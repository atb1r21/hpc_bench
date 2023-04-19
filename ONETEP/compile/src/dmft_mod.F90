! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!         ____  __  __ _____ _____   __  __           _       _               !
!        |  _ \|  \/  |  ___|_   _| |  \/  | ___   __| |_   _| | ___          !
!        | | | | |\/| | |_    | |   | |\/| |/ _ \ / _` | | | | |/ _ \         !
!        | |_| | |  | |  _|   | |   | |  | | (_) | (_| | |_| | |  __/         !
!        |____/|_|  |_|_|     |_|   |_|  |_|\___/ \__,_|\__,_|_|\___|         !
!                                                                             !
! This module contains the main subroutine for performing dynamical mean      !
! field theory calculations, by interfacing with TOSCAM                       !
!                                                                             !
! For further information contact Cedric Weber                                !
!-----------------------------------------------------------------------------!
! Written by Cedric Weber as hubbard_build_mod_dmft.h                         !
! Moved to dmft_mod.F90 by Edward Linscott Mar 2016                           !
! Modified extensively by Edward Linscott 2017                                !
! Modified by Robert Charlton for embedding structures, September 2018.       !
!-----------------------------------------------------------------------------!
! To do:                                                                      !
!  LONG-TERM                                                                  !
!  - revisit parallelism schemes                                              !
!  - re-introduce GPU capabilities                                            !
!  SHORT-TERM                                                                 !
!  - look at the ien variable - I think we can get away with                  !
!    passing it back and forth less frequently                                !
!  - Cedric: go through CWTODO flags                                          !
!=============================================================================!

module dmft
   use constants, only: DP, VERBOSE, STDOUT, PI
   use dense, only: DEM
   use sparse, only: SPAM3
   use sparse_EMBED, only: SPAM3_EMBED
   use rundat, only: pub_output_detail, PUB_1K

implicit none

! ebl: container type to store all information relating to a loop
type DMFT_LOOP_INFO
   integer :: hat
   integer :: theatom
   integer :: channels
   integer :: sp
   integer :: hatb
   integer :: theatomb
   integer :: channelsb
   integer :: spb
end type DMFT_LOOP_INFO

! ebl: container type to store matrices which are projected onto the impurity,
!      namely the green function and the self energy
type DMFT_MATRIX
   type(SPAM3) :: mat
   type(SPAM3) :: mat_v
   type(SPAM3) :: w_mat
   type(SPAM3) :: w_mat_v
   type(SPAM3) :: inv_mat
   type(SPAM3) :: inv_mat_v
   type(SPAM3) :: w_inv_mat
   type(SPAM3) :: w_inv_mat_v
   type(SPAM3) :: backup_spam
   type(SPAM3) :: inv_backup_spam
   complex(kind=DP), allocatable :: backup(:,:)
   complex(kind=DP), allocatable :: inv_backup(:,:)
   logical     :: backup_flag
   !CWTODO: discuss backup and inv_backup matrices - you had both SPAM3 greenf_backup
   ! and array greenfbackup, discuss how to proceed
end type

! ebl: container type to store the various matrices associated with the overlap
! CWTODO Would it be possible to only use the complex parts? The real parts
! seem to be given by rep%overlap etc so there might be duplication going on here
type DMFT_OVERLAP
   type(SPAM3) :: mat
   type(SPAM3) :: tmat
   type(SPAM3) :: inv_mat
   type(SPAM3) :: inv_tmat
   type(SPAM3) :: matc
   type(SPAM3) :: tmatc
   type(SPAM3) :: inv_matc
   type(SPAM3) :: inv_tmatc
end type

type DMFT_SC_MATRICES
   real(kind=DP),allocatable       :: eigenvals(:,:)        ! eigenvalues of ham (was 'vectorAA')
   real(kind=DP),allocatable       :: eigenvals_backup(:,:) ! backup (was 'vectorAAback')
   type(DEM), allocatable          :: eigenvecs(:)          ! eigenvectors of ham
   type(DEM), allocatable          :: eigenvecs_t(:)        ! ...transposed
   type(SPAM3), allocatable        :: ham(:)                ! hamiltonian
   type(SPAM3)                     :: overlap               ! overlap
   type(DEM)                       :: tail                  ! Green function tail contribution
end type

! ebl: container type that allows the size of the occupancy matrix to be
!      dependent on the atom type (or more precisely, l)
type ARRAY_OF_MATRICES
   real(kind=DP), allocatable :: occupancy(:,:,:)
end type

! ebl: container type to store matrices associated with k-point parallelisation
type DMFT_KPOINT_MATRICES
   type(SPAM3) :: ham               ! Hamiltonian (was 'hamK')
   type(SPAM3) :: overlap           ! NGWF overlap in rotated frame (was 'ovK')
   type(SPAM3) :: inv_overlap       ! inverse NGWF overlap in rotated frame (was 'invK')
   type(SPAM3) :: self_energy       ! self energy (was 'selfK')
   type(SPAM3) :: hub_overlap       ! hubbard overlap (was 'O_bar')
   type(SPAM3) :: unrot_overlap     ! unrotated NGWF overlap (was 'Ok')
   type(SPAM3) :: unrot_inv_overlap ! unrotated NGWF inverse overlap (was 'invOk')
end type

   private

   public :: hubbard_dmft_interface

   integer, public :: pub_dmft_mpi_size
   integer, public :: pub_dmft_mpi_rank
   ! ebl: TODO make GPU logicals cleaner
   logical, public :: pub_dmft_use_gpu
   logical, public :: pub_dmft_use_gpu_onlyme
   logical, public :: pub_dmft_use_gpu_partially

   ! ebl: keywords that are, for the moment, only included as devel_code
   real(kind=DP)     :: pub_dmft_chem_shift
   real(kind=DP)     :: pub_dmft_cutoff_tail
   real(kind=DP)     :: pub_dmft_doping
   integer           :: pub_dmft_embed_iter
   real(kind=DP)     :: pub_dmft_embed_mix
   real(kind=DP)     :: pub_dmft_free_green_frequ
   integer           :: pub_dmft_gpu_num
   logical           :: pub_dmft_impose_chem_spin
   logical           :: pub_dmft_impose_same_coeffs
   logical           :: pub_dmft_in_bohr
   logical           :: pub_dmft_integrate_green
   logical           :: pub_dmft_invert_overlap
   logical           :: pub_dmft_kpoints_kernel_gamma
   logical           :: pub_dmft_lin_scaling
   logical           :: pub_dmft_local_scratch
   integer           :: pub_dmft_norm_proj
   real(kind=DP)     :: pub_dmft_optics_x1
   real(kind=DP)     :: pub_dmft_optics_y1
   real(kind=DP)     :: pub_dmft_optics_z1
   logical           :: pub_dmft_plot_all_proj
   logical           :: pub_dmft_plot_real_space_sigma
   real(kind=DP)     :: pub_dmft_scaling_cutoff_h
   integer           :: pub_dmft_scaling_iter
   integer           :: pub_dmft_scaling_maxspace
   real(kind=DP)     :: pub_dmft_scaling_tol
   logical           :: pub_dmft_split
   logical           :: pub_dmft_splitk

contains
 subroutine hubbard_dmft_interface(eigen_en_input, n_occ_input, ham, &
      overlap_input, inv_overlap_input, ngwf_basis, hub_proj_basis, hub, rep, &
      elements, denskern, mdl, nl_projectors, proj_basis, dmft_energy_cor, &
      dmft_self, dmft_kernel, dmft_z)

   !---------------------------------------------------------------------------!
   ! This interface is responsible for half of the DMFT self-consistent cycle. !
   ! It...                                                                     !
   !    1) inverts the impurity Green's function                               !
   !    2) calculates the impurity self-energy                                 !
   !    3) upfolds the self-energy to the full system                          !
   !    4) updates the total Green's function (and potentially the chemical    !
   !       potential and the Hamiltonian)                                      !
   !    5) projects the Green's function, self energy, and Hamiltonian onto    !
   !       the Hubbard subspaces                                               !
   !                                                                           !
   ! Having completed this, it produces the files TOSCAM requires to generate  !
   ! and solve the corresponding Anderson impurity model (the other half of    !
   ! the DMFT self-consistent cycle).                                          !
   !                                                                           !
   ! N.B. for the first loop, the self-energy is assumed to be zero and only   !
   ! step 5 is performed                                                       !
   !                                                                           !
   ! The code itself is structured in three major steps:                       !
   !    1) initialisation of variables                                         !
   !    2) a triple-loop over frequencies, chemical potential, and spin        !
   !    3) post-processing and wiping of variables                             !
   !---------------------------------------------------------------------------!
   ! HISTORY                                                                   !
   !   Written by Cedric Weber                                                 !
   !   Cleaned extensively by Edward Linscott Mar-Aug 2017                     !
   !---------------------------------------------------------------------------!

   use comms, only: comms_barrier, comms_reduce, pub_my_proc_id, pub_on_root, &
        pub_total_num_procs
   use file_handling, only: file_handling_remove_file, &
        file_handling_remove_numbered_files
   use function_basis,      only: FUNC_BASIS
   use hubbard_build,       only: HUBBARD_MODEL
   use geometry,            only: POINT
   use ion,                 only: ELEMENT
   use kernel,              only: DKERN, kernel_create
   use model_type,          only: MODEL
   use ngwf_representation, only: NGWF_REP
   use parallel_strategy,   only: par => pub_par
   use projectors,          only: PROJECTOR_SET
   use rundat,              only: pub_debug, pub_debug_on_root, &
        pub_dmft_complex_freq, pub_dmft_dos_max, pub_dmft_dos_min, &
        pub_dmft_kernel, pub_dmft_nmu_loop, pub_dmft_paramagnetic, &
        pub_dmft_plot_real_space, pub_dmft_points, pub_dmft_read, pub_dmft_sc, &
        pub_dmft_temp, pub_num_kpoints, pub_num_spins, pub_rootname
   use sparse, only: sparse_copy, sparse_create, sparse_scale, &
        sparse_show_matrix
   use sparse_embed,        only: SPAM3_EMBED, sparse_embed_create, &
        sparse_embed_array_scale
   use utils,               only: utils_abort, utils_assert, utils_unit, &
        utils_banner, utils_alloc_check, utils_dealloc_check, utils_int_to_str

   implicit none

   ! Arguments
   ! molecular orbital energies
   real(kind=DP), allocatable, intent(in   ) :: eigen_en_input(:,:)
   ! number of occupied orbitals
   integer,                    intent(in   ) :: n_occ_input(:)
   type(SPAM3),                intent(inout) :: ham(pub_num_spins)
   type(SPAM3_EMBED),          intent(in   ) :: overlap_input
   type(SPAM3_EMBED),          intent(in   ) :: inv_overlap_input
   type(FUNC_BASIS),           intent(in   ) :: ngwf_basis(1)
   type(FUNC_BASIS),           intent(in   ) :: hub_proj_basis
   type(HUBBARD_MODEL),        intent(inout) :: hub
   type(NGWF_REP),             intent(in   ) :: rep
   type(ELEMENT),              intent(in   ) :: elements(:)
   type(DKERN),                intent(inout) :: denskern
   type(MODEL),                intent(in   ) :: mdl
   type(PROJECTOR_SET),        intent(inout) :: nl_projectors(1)
   ! Optional arguments for on-the-fly DMFT calculations
   ! This functionality is not yet implemented
   type(FUNC_BASIS), optional, intent(in   ) :: proj_basis(1)
   ! total DMFT energy
   real(kind = DP), optional,  intent(inout) :: dmft_energy_cor
   ! self energy
   type(SPAM3), optional,      intent(inout) :: dmft_self(pub_num_spins)
   ! density kernel evaluated via trace of Green's function
   type(SPAM3), optional,      intent(inout) :: dmft_kernel(pub_num_spins)
   ! quasiparticle weights
   type(SPAM3), optional,      intent(inout) :: dmft_z(pub_num_spins)

   ! Local variables
   type(DMFT_MATRIX)               :: self_energy
   type(DMFT_MATRIX)               :: greenf
   type(DMFT_OVERLAP)              :: hub_overlap
   type(DMFT_KPOINT_MATRICES)      :: kpoint_data
   type(DMFT_SC_MATRICES)          :: sc_data

   real(kind = DP), allocatable    :: occupancy(:,:,:)
   real(kind = DP), allocatable    :: eigen_en(:,:)
   integer, allocatable            :: n_occ(:)
   integer                         :: ii, num_spins, i
   real(kind = DP), allocatable    :: nabla_square(:,:,:), temp1_real(:,:), &
                                      temp2_real(:,:)
   complex(kind = DP), allocatable :: temp2_complex(:,:)
   type(SPAM3_EMBED)               :: nabla(3)
   type(SPAM3)                     :: ham_embed(pub_num_spins)
   ! rc2013: annoyingly we need this structure for kernel_purify...
   type(SPAM3_EMBED)               :: overlap
   type(SPAM3_EMBED)               :: inv_overlap
   complex(kind = DP), allocatable :: site_greenf_buffer(:,:,:,:)
   type(SPAM3)                     :: overlap_cplx
   integer                         :: is
   integer                         :: ien
   real(kind = DP)                 :: en_start, en_step
   complex(kind = DP)              :: energy
   logical                         :: check
   logical                         :: real_space_exit
   integer                         :: error
   integer                         :: num
   real(kind = DP)                 :: dist
   character(200)                  :: value(100)
   logical                         :: cluster
   integer                         :: j, totn
   integer                         :: hub_atom
   integer                         :: hub_atomb
   real(kind = DP)                 :: fermi_e(2)
   real(kind = DP)                 :: en_step_over_pi
   integer, parameter              :: nmerge = 3
   integer                         :: k, l, k2, l2, i2, jj
   integer, allocatable            :: dimer_table(:,:), merge_table(:,:), &
                                      merge_neigh(:)
   real(kind = DP), allocatable    :: rot_vec_angles(:,:,:), &
                                      rot_vec_angles_split(:,:,:)
   logical, allocatable            :: connection_table(:,:)
   integer                         :: ffirst, llast, iien, energyzero
   integer, allocatable            :: myfrequencies(:)
   integer                         :: myfrequencies_numb
   complex(kind = DP), allocatable :: overlap_hubbard(:,:,:)
   integer, parameter              :: nspecial_frequ = 7
   type(ARRAY_OF_MATRICES)         :: h_atoms(par%nat_hub)
   real(kind = DP)                 :: E_sum_eigenvalues_DFT(2), total_energy, &
                                      lhxc_energy, hubbard_energy, &
                                      dmft_Energy(2), edc_correction
   type(DEM), allocatable          :: matrixSigStat(:), matrixAsig(:)
   logical                         :: force_reopen_file
   real(kind = DP)                 :: dmft_energy_cor_(2)
   type(SPAM3)                     :: self_infinity(pub_num_spins)

   ! Chemical potential, kernel, and total n (self-consistency variables)
   ! Lump into a type?
   integer                         :: kkk1, kkk2, nmu_step
   real(kind = DP)                 :: target_N
   integer                         :: nmu
   real(kind = DP)                 :: chem_shift
   logical                         :: muconverged, muconvergedspin
   real(kind = DP)                 :: vectorMu(1000), vectorNmu(1000)

   ! Kernel
   type(DKERN)                     :: denskern_tmp
   type(DKERN)                     :: denskern_tmpc
   type(SPAM3)                     :: matrixDensKernel
   real(kind = DP)                 :: dmft_kernel_diis_c_out

   ! kpoints
   integer                         :: kstart, kstep
   logical                         :: shifted
   integer                         :: ikpoint, kp_file_shift
   real(kind = DP), allocatable    :: norm_kpoints(:), kpt(:,:)
   integer                         :: lastkpoint, nkpoints

   ! Embedding
   logical                         :: check_embed, check_embed_h
   integer, allocatable            :: maskH1_embed(:,:,:), &
                                      maskH0_embed(:,:,:), &
                                      maskH0_embed_atoms(:,:,:)
   integer                         :: nplanes
   logical, allocatable            :: orb_in_plane(:,:), atom_in_plane(:,:)
   logical                         :: wont_compute_residues
   real(kind = DP)                 :: tempvar, totkpoint
   integer                         :: free_unit
   logical                         :: same_self_for_all

   ! Introduce a split type?
   ! type(DMFT_SPLIT_INFO)             :: split_data
   complex(kind = DP), allocatable :: split_hubt(:,:), split_hub(:,:), &
                                      split_wgv(:,:), split_Sm(:,:), &
                                      split_S(:,:), split_H(:,:)
   integer                         :: splitk_start, splitk_step, splitk_end, &
                                      splitk_last_batch, dmft_splitk_iter, &
                                      splitk_first
   integer                         :: cluster_scheme
   character(2000)                 :: scratchdir

   ! ebl: Run options, think about how to include these
   logical                         :: reopen_sigma_file_each_step ! keep this
   ! local
   logical, parameter              :: one_k_one_file = .false. ! Keep as a
   ! parameter here, rename?, and pass as an argument
   logical, parameter              :: full_prec_dmft_energy = .false.
   logical, parameter              :: use_gpu_eigenvectors = .true. ! Keep as
   ! a parameter, pass as argument
   logical, parameter              :: optimisation_parallel = .false. ! Keep as
   ! a parameter, pass as argument
   logical                         :: inv_projectors
   logical                         :: restrict_window ! pass as arguments
   logical, parameter              :: show_matrices = .false.

   ! ebl: variables for plotting, hopefully these can become local
   ! when the plotting subroutine is restructured
   !CW: this could be a type...: TYPE: DENSITY_FINE(x,y,z), TYPE(DENSITY_FINE) :: sigma(pub_num_spins), ...
   real(kind = DP), allocatable    :: sigma_density_fine(:,:,:,:)
   real(kind = DP), allocatable    :: green_density_fine(:,:,:,:)
   real(kind = DP), allocatable    :: green_density_fine_proj(:,:,:,:)
   real(kind = DP), allocatable    :: green_tot_fine(:,:,:,:)
   real(kind = DP), allocatable    :: Z_density_fine(:,:,:,:)
   ! Convert these to SPAM3_ARRAYS?
   !CW: spam 3 arrays with nk
   type(SPAM3)                     :: sigma_density_spam(2)
   type(SPAM3)                     :: green_density_spam(2,1)
   type(SPAM3)                     :: green_density_spam_proj(2)
   type(SPAM3)                     :: green_tot_spam(2)
   type(SPAM3)                     :: Z_density_spam(2)

   integer                         :: ierr
   character(len=8)                :: spin_label

   !---------------------------------------------------------------------------!
   !                                                                           !
   !                           STEP 1: INITIALISATION                          !
   !                                                                           !
   !---------------------------------------------------------------------------!

   ! ebl: printing header
   if (pub_on_root) write(stdout,'(/a/)') utils_banner('=', ' Dynamical Mean &
        &Field Theory ')

   call utils_assert(pub_num_kpoints == PUB_1K, "DMFT code not compatible with &
        &DFT calculation performed for multiple k-points. Exiting...")

   ! ebl: parse all devel_code variables
   call dmft_load_devel_code_variables()

   ! ebl: Set parallel scheme (faster way to get the energy when running mpi)
   pub_dmft_split = (present(dmft_energy_cor) .and. pub_total_num_procs > 1) &
        .or. (pub_dmft_split .and. pub_total_num_procs > 1)

   if (pub_dmft_split .and. pub_on_root) then
      write(stdout, '(a)') 'Splitting the frequency summation over procs'
   endif

   call dmft_set_kpoints(kpt, norm_kpoints, nkpoints, lastkpoint, totkpoint, &
        mdl)

   ! Catch input arguments and set module variables
   inv_projectors = .false.
   if (pub_debug) then
      inv_projectors = .true.
   end if
   if(inv_projectors)then
      pub_dmft_invert_overlap = .true.
   endif

   call utils_assert(.not.(pub_dmft_kernel == 0 .and. &
        present(dmft_energy_cor)), "Error in hubbard_dmft_interface: cannot &
        &compute the dmft energy but not the dmft kernel")

   nmu = pub_dmft_nmu_loop
   if(nmu == 0 .and. (pub_dmft_kernel /= 0 .or. present(dmft_energy_cor)))then
      nmu = 1
   endif

   if (pub_dmft_lin_scaling) then
      nmu = 0
      pub_dmft_kernel = 0
#ifndef LINEAR_SCALING_INV
      call utils_abort('Error in hubbard_dmft_interface: if &
           &pub_dmft_lin_scaling is to be set to T, please compile with &
           &preprocessing flag LINEAR_SCALING_INV')
#endif
   endif

   ! Allocation of variables
   if(pub_dmft_kernel /= 0)then
      call kernel_create(denskern_tmp, denskern%kern%m(1, 1)%structure // &
           rep%postfix, .false.)
      call sparse_embed_array_scale(denskern_tmp%kern, 0.0_DP)
   endif
   call sparse_embed_create(     overlap,     overlap_input, iscmplx = .false. )
   call sparse_embed_create( inv_overlap, inv_overlap_input, iscmplx = .false. )
   call sparse_create( hub_overlap%mat,   rep%hub_overlap%p, iscmplx = .false. )
   call sparse_create( hub_overlap%tmat, rep%hub_overlap_t%p, iscmplx = .false. )
   do is = 1, pub_num_spins
      call sparse_create( ham_embed(is), ham(is), iscmplx = .true. )
   enddo

   if(nkpoints > 1)then
      call sparse_create(kpoint_data%overlap,     overlap%p , iscmplx = .true.)
      call sparse_create(kpoint_data%ham,         ham(1)  , iscmplx = .true.)
      call sparse_create(kpoint_data%inv_overlap, inv_overlap%p, iscmplx = .true.)
      if(pub_dmft_kernel /= 0)then
         call sparse_create(denskern_tmpc%kern%m(1, PUB_1K)%p, &
              denskern%kern%m(1, PUB_1K)%p, iscmplx = .true.)
      endif
   endif

   ! Where possible, read in variables from file
   inquire(file = trim(pub_rootname) // '.overlap', exist = check)
   if(.not.check .or. .not. pub_dmft_read) then
      call sparse_copy(               overlap%p,     overlap_input%p  )
   endif
   inquire(file = trim(pub_rootname) // '.inv_overlap', exist = check)
   if(.not.check .or. .not. pub_dmft_read)then
      call sparse_copy(           inv_overlap%p, inv_overlap_input%p  )
   endif
   inquire(file = trim(pub_rootname) // '.hub_overlap', exist = check)
   if(.not.check .or. .not. pub_dmft_read)then
      call sparse_copy(   hub_overlap%mat,   rep%hub_overlap%p  )
   endif
   inquire(file = trim(pub_rootname) // '.hub_overlap_t', exist = check)
   if(.not.check .or. .not. pub_dmft_read)then
      call sparse_copy(   hub_overlap%tmat, rep%hub_overlap_t%p  )
   endif

   num_spins = pub_num_spins
   if(pub_dmft_paramagnetic) num_spins = 1

   ! Depreciating gpu_max file in favour of the dmft_gpu_num flag within ONETEP
   ! inquire(file = 'gpu_max', exist = check)
   ! if(check)then
   !    if(pub_on_root) write(stdout, *) 'gpu_max file present'
   !    free_unit = utils_unit()
   !    open(unit = free_unit, file = 'gpu_max')
   !    read(free_unit, *) gpu_max
   !    close(free_unit)
   ! else
   !    gpu_max = 8
   !    write(stdout, *) 'WARNING gpu_max file not present, assuming there are 8 &
   !         &GPUS'
   ! endif

   if(pub_dmft_points <= nspecial_frequ + 1)then
      write(stdout, '(a, i4, a)') 'Please use at least ', nspecial_frequ + 2, &
           ' dmft points (pub_dmft_points)'
      write(stdout, '(a, i4, a)') 'The last ', nspecial_frequ, ' frequencies &
           &are used for testing purposes'
      call utils_abort("Error in hubbard_dmft_interface: dmft_points too small")
   endif

   cluster = .false.
   inquire(file = "mask_dimer", exist = check)
   call comms_barrier
   cluster = check
   if(cluster .and. pub_total_num_procs /= 1)then
      if(pub_on_root) write(stdout, '(a)') 'WARNING, dimer mode was removed &
           &because of mpi use'
      cluster = .false.
      check = .false.
   endif

   ! If performing a cluster DMFT calculation, finds pairs of atoms
   call dmft_define_cluster(rot_vec_angles, rot_vec_angles_split, dimer_table, &
        nmerge, merge_table, merge_neigh, connection_table, hub, &
        hub_proj_basis, mdl%cell, elements, cluster)

   ! Initialises various routine variables
   call dmft_init_routine_variables(ham, overlap%p, h_atoms, hub, hub_overlap, &
        hub_proj_basis, overlap_hubbard, num_spins, en_start, en_step, &
        en_step_over_pi, energyzero, nspecial_frequ, cluster, cluster_scheme, &
        force_reopen_file, inv_overlap%p, reopen_sigma_file_each_step, &
        restrict_window, same_self_for_all, use_gpu_eigenvectors)

   ! Sets up a scratch directory if requested
   call dmft_init_scratch_dir(scratchdir)

   ! Initialisation of sparse and dense matrices
   call dmft_init_matrices(ham, greenf, self_energy, self_infinity, overlap%p, &
        hub_overlap, ngwf_basis(1), kpoint_data, sc_data, overlap_cplx, &
        matrixAsig, matrixSigStat, matrixDensKernel, nkpoints, inv_overlap%p, &
        sigma_density_spam, green_density_spam, green_density_spam_proj, &
        green_tot_spam, Z_density_spam, cluster_scheme, inv_projectors)

   ! Construct the connection tables for dimers
   call dmft_build_my_dimers(dimer_table, cluster, hub, hub_proj_basis, &
        mdl%cell, elements, nmerge, merge_table, merge_neigh)

   ! Initialise quantities required for optics calculations
   call dmft_init_optics(nabla_square, nabla, overlap%p, ngwf_basis, rep, mdl, &
        nl_projectors, proj_basis)

   ! Check if we can load the hamiltonian from a file
   inquire(file = trim(pub_rootname) // '.eigen', exist = check)
   call comms_barrier
   call dmft_read_write_hamiltonian(eigen_en, n_occ, eigen_en_input, &
        n_occ_input, h_atoms, hub, hub_proj_basis, num_spins, check)
   call dmft_print_hubbard_atom_info(eigen_en, n_occ, h_atoms, hub, &
        hub_proj_basis, num_spins, rot_vec_angles)

   ! Fix the fermi energy, and update total energies
   if (present(dmft_energy_cor)) then
      call dmft_fix_fermi_and_total_energies(total_energy, mdl, ngwf_basis(1), &
           nkpoints, E_sum_eigenvalues_DFT, sc_data, fermi_e, eigen_en, n_occ, &
           overlap%p, ham, full_prec_dmft_energy, use_gpu_eigenvectors, &
           dmft_energy_cor)

   else
      call dmft_fix_fermi_and_total_energies(total_energy, mdl, ngwf_basis(1), &
           nkpoints, E_sum_eigenvalues_DFT, sc_data, fermi_e, eigen_en, n_occ, &
           overlap%p, ham, full_prec_dmft_energy, use_gpu_eigenvectors)
   end if
   dmft_energy_cor_ = 0.0_DP

   ! If performing an embedding calculation, set up variables
   call dmft_embedding_prepare_tables(maskH0_embed, maskH1_embed, &
        maskH0_embed_atoms, ngwf_basis(1), ham, overlap%p, nkpoints, check_embed)

   ! Define the Matsubara (or real) frequencies
   call dmft_define_frequencies(myfrequencies, myfrequencies_numb, ffirst, &
        llast, fermi_e, energy, en_start, en_step, check_embed, &
        restrict_window)

   ! In an embedding calculation, define the planes which interface the
   ! calculations with an interface
   call dmft_embedding_setup_planes(orb_in_plane, atom_in_plane, nplanes, &
        ngwf_basis(1), par%nat)

   if(.not.pub_dmft_sc) then

      !------------------------------------------------------------------------!
      !                                                                        !
      !      STEP 2: LOOPS OVER CHEMICAL POTENTIAL, FREQUENCIES, AND SPIN      !
      !                                                                        !
      !------------------------------------------------------------------------!

      nmu_step = 0
      vectorMu = 0.0_DP
      vectorNmu = 0.0_DP
      chem_shift = 0.0_DP

      muconvergedspin = .false.

      if (pub_on_root) then
         write(stdout, "(/a)") repeat(" ",19)//repeat("-",43)
         write(stdout, "(a)") repeat(" ",19)//"| chemical potential | &
              &     occupancy     |"
         write(stdout, "(a)") repeat(" ",19)//"|"//repeat("-",20)//"|"&
              //repeat("-",20)//"|"
      end if

      do while(.not.muconvergedspin)

         if(pub_dmft_complex_freq .and. .not.pub_dmft_impose_chem_spin) then
            ! Update the chemical potential
            call dmft_fix_chemical_potential_loop(n_occ, chem_shift, target_N, &
                 fermi_E, vectorMu, vectorNmu, num_spins, nmu_step, nmu)
         end if

         do is = 1, num_spins
            if (pub_debug_on_root) write(stdout, '(a,i1)') 'DEBUG: Entering &
                 &spin loop for spin ', is

            if(pub_dmft_complex_freq .and. pub_dmft_impose_chem_spin) then
               nmu_step = 0
               vectorMu = 0.0_DP
               vectorNmu = 0.0_DP
               chem_shift = 0.0_DP
            endif

            muconverged = .false.
            do while (.not.muconverged)

               if(pub_dmft_complex_freq .and. pub_dmft_impose_chem_spin) then
                  ! Update chemical potential for this spin channel
                  call dmft_fix_chemical_potential_loop_spin(n_occ, &
                       chem_shift, target_N, fermi_E, is, vectorMu, vectorNmu, &
                       num_spins, nmu_step, nmu)
               end if

               ! Set up k-point parallelism
               if(pub_dmft_splitk)then
                  splitk_start = pub_dmft_mpi_rank-1
                  splitk_step = pub_dmft_mpi_size
                  splitk_end = nkpoints/pub_dmft_mpi_size + 1 !number of batch
                  splitk_end = splitk_end*pub_dmft_mpi_size
                  splitk_last_batch = splitk_end-pub_dmft_mpi_size + 1 ! >=
                  ! splitk_last_batch is in last batch
                  lastkpoint = splitk_last_batch
                  dmft_splitk_iter = 0
                  if(pub_on_root)then
                     write(stdout, *) 'SPLITTING K-POINTS - START/END - my &
                          &rank : ', splitk_start, splitk_end, &
                          pub_dmft_mpi_rank
                     write(stdout, *) 'SPLITTING K-POINTS - STEP : ', &
                          splitk_step
                     write(stdout, *) 'SPLITTING K-POINTS - LAST BATCH START &
                          &AT : ', splitk_last_batch
                     write(stdout, *) 'BREAK NFS COMM : ', pub_dmft_splitk
                  endif
                  splitk_first = nkpoints + 1
                  do ikpoint = splitk_start + 1, splitk_end, splitk_step
                     if(abs(norm_kpoints(ikpoint)) > 1.0E-4_DP) then
                        splitk_first = ikpoint
                        exit
                     endif
                  enddo
                  write(stdout, *) 'SPLITTING K-POINTS - FIRST K POINT : ', &
                       splitk_first
               else
                  splitk_start = 0
                  splitk_step = 1
                  splitk_end = nkpoints
               endif

               do ikpoint = splitk_start + 1, splitk_end, splitk_step

                  ! Reset kernel and energy
                  call sparse_scale(matrixDensKernel, 0.0_DP)
                  dmft_Energy(is) = 0.0_DP

                  ! Prepare the arrays used in the calculation of the energy
                  ! for a new iteration of the chemical potential
                  if(pub_dmft_splitk)then
                     if(abs(norm_kpoints(ikpoint)) < 1.0E-4_DP) then
                        write(stdout, "(a)") 'Skipping k-point due to&
                             & symmetry'
                        if(pub_dmft_complex_freq .and. nmu > 0) then
                           if (present(dmft_energy_cor)) then
                              call dmft_clean_dmft_energy(ham, overlap%p, &
                                   ngwf_basis(1), sc_data, matrixDensKernel, &
                                   matrixAsig, matrixSigStat, denskern_tmp, &
                                   denskern_tmpc, dmft_Energy, total_energy, &
                                   edc_correction, E_sum_eigenvalues_DFT, &
                                   fermi_e, is, num_spins, ikpoint, nkpoints, &
                                   kpt, totkpoint, lastkpoint, nmu, nmu_step, &
                                   chem_shift, vectorNmu, vectorMu, target_N, &
                                   norm_kpoints, dmft_splitk_iter, &
                                   optimisation_parallel, dmft_energy_cor_, &
                                   dmft_energy_cor)
                           else
                              call dmft_clean_dmft_energy(ham, overlap%p, &
                                   ngwf_basis(1), sc_data, matrixDensKernel, &
                                   matrixAsig, matrixSigStat, denskern_tmp, &
                                   denskern_tmpc, dmft_Energy, total_energy, &
                                   edc_correction, E_sum_eigenvalues_DFT, &
                                   fermi_e, is, num_spins, ikpoint, nkpoints, &
                                   kpt, totkpoint, lastkpoint, nmu, nmu_step, &
                                   chem_shift, vectorNmu, vectorMu, target_N, &
                                   norm_kpoints, dmft_splitk_iter, &
                                   optimisation_parallel, dmft_energy_cor_)
                           end if
                        end if
                        cycle
                     endif
                     if(ikpoint > nkpoints)then
                        write(stdout, *) 'DMFT SPLIT K : skipping k-point - &
                             &outbound index'
                        if(pub_dmft_complex_freq .and. nmu > 0) then
                           if (present(dmft_energy_cor)) then
                              call dmft_clean_dmft_energy(ham, overlap%p, &
                                   ngwf_basis(1), sc_data, matrixDensKernel, &
                                   matrixAsig, matrixSigStat, denskern_tmp, &
                                   denskern_tmpc, dmft_Energy, total_energy, &
                                   edc_correction, E_sum_eigenvalues_DFT, &
                                   fermi_e, is, num_spins, ikpoint, nkpoints, &
                                   kpt, totkpoint, lastkpoint, nmu, nmu_step, &
                                   chem_shift, vectorNmu, vectorMu, target_N, &
                                   norm_kpoints, dmft_splitk_iter, &
                                   optimisation_parallel, dmft_energy_cor_, &
                                   dmft_energy_cor)
                           else
                              call dmft_clean_dmft_energy(ham, overlap%p, &
                                   ngwf_basis(1), sc_data, matrixDensKernel, &
                                   matrixAsig, matrixSigStat, denskern_tmp, &
                                   denskern_tmpc, dmft_Energy, total_energy, &
                                   edc_correction, E_sum_eigenvalues_DFT, &
                                   fermi_e, is, num_spins, ikpoint, nkpoints, &
                                   kpt, totkpoint, lastkpoint, nmu, nmu_step, &
                                   chem_shift, vectorNmu, vectorMu, target_N, &
                                   norm_kpoints, dmft_splitk_iter, &
                                   optimisation_parallel, dmft_energy_cor_)
                           end if
                        end if
                        cycle
                     endif
                  endif

                  if(abs(norm_kpoints(ikpoint)) < 1.0E-4_DP) then
                     write(stdout, *) 'SKIPPING K POINT DUE TO SYMMETRY'
                     cycle
                  endif

                  ! write(stdout, "(/a)") ' Performing k-point number '// &
                  !      trim(utils_int_to_str(ikpoint))
                  ! if(ikpoint <= nkpoints)then
                  !    write(*, '(a, 3f14.4)') ' KP : ', kpt(ikpoint, :)
                  ! endif

                  ! Prepare and define the Hamiltonian and overlap
                  ! for multiple k-point calculation
                  call dmft_load_k_ham_and_overlap(ham, greenf, self_energy, &
                       overlap%p, inv_overlap%p, overlap_cplx, h_atoms, hub, &
                       hub_overlap, hub_proj_basis, ngwf_basis(1), kpoint_data, &
                       energy, fermi_e, is, ien, ikpoint, ffirst, llast, &
                       kp_file_shift, hub_atom, hub_atomb, nkpoints, &
                       matrixSigStat, rot_vec_angles, connection_table, &
                       dimer_table, nmerge, merge_table, merge_neigh, &
                       ham_embed, check_embed, check_embed_h, split_hub, &
                       split_hubt, split_wgv, split_S, split_Sm, split_H, &
                       cluster, force_reopen_file, inv_projectors, &
                       one_k_one_file, reopen_sigma_file_each_step, &
                       same_self_for_all, show_matrices, scratchdir)

                  if(pub_dmft_complex_freq .and. nmu > 0) then
                     ! Prepare the arrays used for the calculation of the
                     ! density kernel
                     call dmft_prepare_for_dmft_density_kernel(ham, greenf, &
                          self_energy, self_infinity, overlap%p, h_atoms, hub, &
                          hub_overlap, hub_proj_basis, ngwf_basis(1), &
                          kpoint_data, sc_data, fermi_e, is, ien, ffirst, &
                          llast, hub_atom, hub_atomb, nkpoints, kkk1, kkk2, &
                          en_start, en_step, nspecial_frequ, nmu_step, &
                          matrixSigStat, matrixAsig, overlap_cplx, &
                          rot_vec_angles, connection_table, dimer_table, &
                          cluster, nmerge, merge_table, merge_neigh, &
                          ham_embed, check_embed, check_embed_h, &
                          force_reopen_file, full_prec_dmft_energy, &
                          optimisation_parallel, reopen_sigma_file_each_step, &
                          same_self_for_all, show_matrices, &
                          use_gpu_eigenvectors, scratchdir, dmft_self, dmft_z)

                  end if

                  real_space_exit = .false.

                  ! Compute the NGWF overlap projected onto the Hubbard subspace
                  if (pub_debug) call dmft_compute_projected_overlap(greenf, &
                       hub_overlap, hub, hub_proj_basis, overlap_cplx, &
                       overlap_hubbard, is, nmu_step, ikpoint, ien, energy, &
                       inv_projectors)

                  ! Setting up frequency parallelism
                  if (pub_dmft_split) then
                     kstart = pub_my_proc_id
                     kstep = pub_total_num_procs
                     if(mod(pub_dmft_points, kstep) /= 0) then
                        write(stdout, *) 'Total number of procs : ',  kstep
                        write(stdout, *) 'Total number of frequ : ', &
                             pub_dmft_points
                        call utils_abort('Error in hubbard_dmft_interface: the &
                             &total number of procs and frequencies should be &
                             &commensurate. Exiting...')
                     endif
                  else
                     kstart = 0
                     kstep = 1
                  endif

                  ! Loop over frequencies
                  if (pub_debug_on_root) write(stdout, "(a)") &
                       'DEBUG: Looping over frequencies ... '

                  do iien = kstart + 1, myfrequencies_numb, kstep
                     ien = myfrequencies(iien)

                     if(present(dmft_energy_cor) .and. pub_dmft_split)then
                        if(iien > pub_dmft_points-nspecial_frequ) then
                           write(stdout, *) '...processor done with &
                                &dmft_energy part so cycle now...'
                           exit
                        endif
                     endif
                     if(pub_debug_on_root) write(stdout, '(a, i4, a, i4)') &
                          'DEBUG: Frequency ', iien, ' of ', myfrequencies_numb

                     ! Define the energy for the current loop index
                     call dmft_define_energy(energy, fermi_e, is, ien, &
                          en_start, en_step)

                     if( .not. restrict_window .or. pub_dmft_split .or. &
                          check_embed .or. ( real(energy) - (fermi_e(is) + &
                          pub_dmft_chem_shift) > pub_dmft_dos_min .and. &
                          real(energy) - (fermi_e(is) + pub_dmft_chem_shift) < &
                          pub_dmft_dos_max ) ) then

                        call dmft_embedding_read_potential(ham, ham_embed, is, &
                             ien, fermi_e, ffirst, llast, &
                             reopen_sigma_file_each_step, force_reopen_file, &
                             check_embed, check_embed_h, ngwf_basis(1))

                        ! Build the inverted Green's function
                        call dmft_build_inv_green_matrix(greenf, overlap_cplx, &
                             ham, kpoint_data%ham, ham_embed, is, ien, &
                             nkpoints, energy, check_embed_h, split_S, &
                             split_H)

                        ! Calculate the local self energy and upfold it to the
                        ! NGWF subspace
                        call dmft_wrapper_build_self_energy(self_energy, &
                             h_atoms, hub_overlap, hub, hub_proj_basis, is, &
                             ien, iien, nkpoints, hub_atom, hub_atomb, &
                             cluster, connection_table, nmerge, merge_table, &
                             merge_neigh, dimer_table, kpoint_data, kstart, &
                             kstep, ffirst, llast, rot_vec_angles, &
                             same_self_for_all, show_matrices, &
                             force_reopen_file, reopen_sigma_file_each_step, &
                             scratchdir)

                        ! Add the self energy to the Green's function
                        ! constructor (omega - H - Sigma) and invert to get the
                        ! Green's function
                        call dmft_wrapper_add_self_and_invert_green(greenf, &
                             self_energy, ngwf_basis(1), ham, overlap%p, &
                             overlap_cplx, inv_overlap%p, energy, fermi_e, ien, &
                             is, ikpoint, nkpoints, norm_kpoints, ffirst, &
                             kp_file_shift, splitk_first, kpoint_data, &
                             split_Sm, split_H, matrixSigStat, nplanes, &
                             orb_in_plane, one_k_one_file)

                        if(ien < pub_dmft_points-nspecial_frequ + 1 .and. &
                             pub_dmft_complex_freq .and. nmu > 0)then
                           ! Compute the DMFT energy via the Migdal formula
                           call dmft_compute_dmft_energy(dmft_Energy, greenf, &
                                self_energy, ngwf_basis(1), kpoint_data, sc_data, &
                                energy, fermi_e, nmu, nmu_step, is, ien, &
                                nkpoints, kkk1, kkk2, nspecial_frequ, &
                                ham_embed, check_embed_h, matrixDensKernel, &
                                matrixSigStat, matrixAsig)
                           if(present(dmft_energy_cor)) cycle
                        endif

                        call dmft_dump_data_optics(nabla_square, greenf, &
                             ngwf_basis(1), mdl, energy, fermi_e(is), chem_shift, &
                             ien, is, nkpoints, ikpoint, norm_kpoints, ffirst, &
                             splitk_first, nspecial_frequ, one_k_one_file)

                        ! Project the Green's function onto Hubbard subspaces
                        call dmft_project_green_function(greenf, hub_overlap, &
                             split_hub, split_hubt, split_wgv)

                        if(.not.pub_dmft_split .and. pub_debug .and. &
                             inv_projectors) then
                           call dmft_project_hamiltonian(greenf, hub_overlap, &
                                inv_projectors)
                        end if

                        if(.not.pub_dmft_split) call comms_barrier

                        ! if(.not.pub_dmft_split) then
                        !    if(show_matrices .and. ien == 9)then
                        !       if(pub_on_root) write(stdout, *) '-------showing &
                        !            &greenf in Atom basis--------'
                        !       call sparse_show_matrix(greenf%w_mat_v, &
                        !            show_elems = .true.)
                        !       if(pub_on_root) write(stdout, *) &
                        !            '-------------------------------------------'
                        !    endif
                        !    if(inv_projectors)then
                        !       if(show_matrices .and. ien == 9)then
                        !          if(pub_on_root) write(stdout, *) &
                        !               '-----showing inv greenf in Atom &
                        !               &basis------'
                        !          call sparse_show_matrix(greenf%w_inv_mat_v, &
                        !               show_elems = .true.)
                        !          if(pub_on_root) write(stdout, *) &
                        !               '-------------------------------------------'
                        !       endif
                        !    endif
                        ! end if

                        ! Gathers and distributes the projected Green functions
                        ! for every Hubbard atom
                        call dmft_build_projected_green(greenf, h_atoms, hub, &
                             hub_proj_basis, mdl%cell, elements, energy, &
                             fermi_e, num_spins, is, ien, ikpoint, ffirst, &
                             norm_kpoints, nkpoints, kp_file_shift, &
                             en_step_over_pi, energyzero, rot_vec_angles, &
                             rot_vec_angles_split, splitk_first, split_wgv, &
                             cluster, connection_table, dimer_table, &
                             merge_neigh, merge_table, nmerge, one_k_one_file)

                        ! Prepares the embedding potential for subsequent
                        ! embedding calculations
                        call dmft_embedding_write_sigma(check_embed, &
                             maskH0_embed, maskH0_embed_atoms, self_energy, &
                             energy, fermi_e, is, ien, ffirst)

                        ! Computes various quantities that we want to plot, for
                        ! this particular spin and frequency. These will be
                        ! collated and plotted once the loops have been
                        ! completed
                        call dmft_build_real_space_quantities( &
                             sigma_density_fine, green_density_fine, &
                             green_density_fine_proj, green_tot_fine, &
                             Z_density_fine, sigma_density_spam, &
                             green_density_spam, green_density_spam_proj, &
                             green_tot_spam, Z_density_spam, greenf, &
                             self_energy, hub_overlap, overlap_cplx, mdl, &
                             nkpoints, is, ien, en_start, en_step, &
                             nspecial_frequ, inv_projectors, real_space_exit)

                        if(pub_dmft_plot_real_space .and. real_space_exit) exit
                     endif
                  end do
                  if (pub_debug_on_root) write(stdout, "(a)") 'DEBUG: Finished &
                       &looping over frequencies'

                  if(pub_dmft_complex_freq .and. nmu > 0 .and. pub_dmft_split) &
                       then
                     call comms_reduce('SUM', dmft_Energy(is))
                     if(nkpoints == 1)then
                        call comms_reduce('SUM', matrixDensKernel%dmtx)
                     else
                        call comms_reduce('SUM', matrixDensKernel%zmtx)
                     endif
                  end if

                  if(pub_dmft_complex_freq .and. nmu > 0) then
                     ! Prepare arrays for a new chemical potential iteration
                     if (present(dmft_energy_cor)) then
                        call dmft_clean_dmft_energy(ham, overlap%p, &
                             ngwf_basis(1), sc_data, matrixDensKernel, &
                             matrixAsig, matrixSigStat, denskern_tmp, &
                             denskern_tmpc, dmft_Energy, total_energy, &
                             edc_correction, E_sum_eigenvalues_DFT, fermi_e, &
                             is, num_spins, ikpoint, nkpoints, kpt, totkpoint, &
                             lastkpoint, nmu, nmu_step, chem_shift, vectorNmu, &
                             vectorMu, target_N, norm_kpoints, &
                             dmft_splitk_iter, optimisation_parallel, &
                             dmft_energy_cor_, dmft_energy_cor)
                     else
                        call dmft_clean_dmft_energy(ham, overlap%p, &
                             ngwf_basis(1), sc_data, matrixDensKernel, &
                             matrixAsig, matrixSigStat, denskern_tmp, &
                             denskern_tmpc, dmft_Energy, total_energy, &
                             edc_correction, E_sum_eigenvalues_DFT, fermi_e, &
                             is, num_spins, ikpoint, nkpoints, kpt, totkpoint, &
                             lastkpoint, nmu, nmu_step, chem_shift, vectorNmu, &
                             vectorMu, target_N, norm_kpoints, &
                             dmft_splitk_iter, optimisation_parallel, &
                             dmft_energy_cor_)
                     end if
                  end if

               end do

               ! Check convergence of chemical potential for this spin channel
               muconverged = .true.
               if(pub_dmft_complex_freq .and. pub_dmft_impose_chem_spin)then
                  if(nmu_step < nmu) muconverged = .false.
               endif
            end do
            if (pub_debug_on_root) write(stdout, '(a,i1)') 'DEBUG: Exiting &
                 &spin loop for spin ', is
         end do

         ! Check convergence of chemical potential
         muconvergedspin = .true.
         if(pub_dmft_complex_freq .and. .not.pub_dmft_impose_chem_spin)then
            if(nmu_step < nmu) muconvergedspin = .false.
         endif
      end do
      if (pub_debug_on_root) write(stdout, "(a)") 'DEBUG: Exiting chemical &
           &potential loop'

      ! ebl: close chemical potential/occupancy table
      if (pub_on_root) then
         write(stdout, "(a/)") repeat(" ",19)//repeat("-",43)
      end if

      ! ebl: print energy summary
      if (pub_on_root .and. pub_dmft_complex_freq) then
         write(stdout, "(a)") utils_banner("-", 'Energy summary')
         write(stdout, "(a,f16.9)") ' Double-counting correction  &
              &                         :', edc_correction
         write(stdout, "(a,f16.9)") ' Total E_LDA energy / 2&
              &                               :', total_energy
         do is = 1,pub_num_spins
            if (pub_num_spins == 2) then
               spin_label = '(spin ' // trim(adjustl(utils_int_to_str(is))) // ')'
            else
               spin_label = '        '
            end if
            write(stdout, "(3a,f16.9)") ' Total E_KS energy ', spin_label, &
                 '                           :', E_sum_eigenvalues_DFT(is)
            write(stdout, "(3a,f16.9)") ' Total DMFT energy ', spin_label, &
                 '                           :', dmft_Energy(is)
            write(stdout, "(3a,f16.9)") ' Double-counting-corrected total DMFT &
                 &energy ', spin_label, ' :', &
                 dmft_Energy(is)-edc_correction/2.0_DP
         end do
         write(stdout, "(a,f16.9)") ' Chemical potential &
              &                                  :', vectorMu(nmu_step)
         write(stdout, "(a/)") repeat("-", 80)
      endif


      !------------------------------------------------------------------------!
      !                                                                        !
      !            STEP 3: POST-PROCESSING AND CLEARING OF VARIABLES           !
      !                                                                        !
      !------------------------------------------------------------------------!

      ! Plotting various quantities of interest
      call dmft_plot_real_space_quantities(sigma_density_fine, &
           green_density_fine, green_density_fine_proj, green_tot_fine, &
           Z_density_fine, sigma_density_spam, green_density_spam, &
           green_density_spam_proj, green_tot_spam, Z_density_spam, overlap%p, &
           ngwf_basis(1), rep, mdl, elements, nkpoints, en_start, en_step, &
           real_space_exit)

      wont_compute_residues = pub_dmft_use_gpu_partially .and. .not. &
           (present(dmft_kernel) .or. present(dmft_energy_cor))

      ! Compute the density kernel via the trace of the Green's function
      ! over Matsubara frequencies
      if (present(dmft_energy_cor) .and. present(dmft_kernel)) then
         call dmft_finalize_kernel(denskern, denskern_tmp, ham, overlap, &
              inv_overlap, rep, num_spins, shifted, wont_compute_residues, &
              dmft_energy_cor, dmft_kernel)
      else if (present(dmft_kernel)) then
         call dmft_finalize_kernel(denskern, denskern_tmp, ham, overlap, &
              inv_overlap, rep, num_spins, shifted, wont_compute_residues, &
              dmft_kernel = dmft_kernel)
      else if (present(dmft_energy_cor)) then
         call dmft_finalize_kernel(denskern, denskern_tmp, ham, overlap, &
              inv_overlap, rep, num_spins, shifted, wont_compute_residues, &
              dmft_energy_cor = dmft_energy_cor)
      else
         call dmft_finalize_kernel(denskern, denskern_tmp, ham, overlap, &
              inv_overlap, rep, num_spins, shifted, wont_compute_residues)
      end if

      if(pub_dmft_use_gpu_partially .and. pub_dmft_kernel /= 0 .and. &
           pub_dmft_complex_freq .and. pub_dmft_mpi_size > 1)then
         tempvar = 0.0_DP
         call utils_abort("Reached point in code where nfs sync used to &
              &appear; the code needs to be updated to use these settings.")
         ! if(.not.pub_dmft_splitk)then
         !    call sync_via_nfs(tempvar, 'sum_tempvar_bef_exit', .true.)
         ! else
         !    if(pub_on_root) call sync_via_nfs(tempvar, 'sum_tempvar_bef_exit', &
         !         .true.)
         !    call comms_barrier
         ! endif
      endif

      if(present(dmft_energy_cor)) dmft_energy_cor = sum(dmft_energy_cor_)

      if(pub_on_root .and. pub_dmft_mpi_size > 1 .and. pub_dmft_mpi_rank == 1) &
           then
         call file_handling_remove_file( &
              "onetep_dmft_confirmation_dmft_density_count")
         call file_handling_remove_numbered_files( &
              "onetep_dmft_confirmation_dmft_density_kernel")
      end if

      ! Updates the dimer connection table to match the ONETEP ordering
      ! of atoms
      call dmft_correct_dimer_file(dimer_table, hub, mdl%cell, elements)
   else
      write(stdout, '(a)') ' EXITING ONETEP DMFT INTERFACE '
   end if

   if (pub_on_root) write(stdout,'(a/)') repeat('=',80)

   ! Deallocating variables in prepartion for exiting
   call dmft_prepare_for_exit(ham, greenf, self_energy, self_infinity, &
        hub_overlap, ngwf_basis(1), kpoint_data, sc_data, energy, &
        myfrequencies, fermi_e, nkpoints, ham_embed, overlap_cplx, &
        denskern_tmp, denskern_tmpc, matrixDensKernel, matrixSigStat, &
        matrixAsig, nabla, nabla_square, orb_in_plane, atom_in_plane, &
        split_hub, split_hubt, split_wgv, split_S, split_Sm, split_H, &
        inv_projectors)

  end subroutine hubbard_dmft_interface

 !=============================================================================!
 ! General subroutines for performing DMFT                                     !
 !=============================================================================!
 ! Written by Cedric Weber as hubbard_build_mod_dmft_routines.h                !
 ! March 2016 - split off from hubbard_build mod and revised by Edward         !
 !              Linscott                                                       !
 ! March 2017 - merged with dmft_gen_routines.h                                !
 ! July 2017  - directly added into dmft_mod.F90                               !
 !=============================================================================!

 subroutine dmft_load_devel_code_variables()
   !---------------------------------------------------------------------------!
   ! Loads devel code variables required by any DMFT subroutine                !
   !---------------------------------------------------------------------------!
   ! NOTES                                                                     !
   ! - Once the corresponding functionalities have been properly checked,      !
   !   each keyword will be moved to rundat_mod                                !
   !---------------------------------------------------------------------------!

   use rundat, only: pub_devel_code
   use utils, only: utils_devel_code

   pub_dmft_chem_shift            = utils_devel_code(  0.0_DP, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_cutoff_tail           = utils_devel_code( 10.0_DP, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_doping                = utils_devel_code(  0.0_DP, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_embed_iter            = utils_devel_code(      10, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_embed_mix             = utils_devel_code(  0.7_DP, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_free_green_frequ      = utils_devel_code(200.0_DP, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_gpu_num               = utils_devel_code(       0, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_impose_chem_spin      = utils_devel_code( .false., "DMFT", "dmft_", pub_devel_code)
   pub_dmft_impose_same_coeffs    = utils_devel_code( .false., "DMFT", "dmft_", pub_devel_code)
   pub_dmft_in_bohr               = utils_devel_code( .false., "DMFT", "dmft_", pub_devel_code)
   pub_dmft_integrate_green       = utils_devel_code( .false., "DMFT", "dmft_", pub_devel_code)
   pub_dmft_invert_overlap        = utils_devel_code( .false., "DMFT", "dmft_", pub_devel_code)
   pub_dmft_kpoints_kernel_gamma  = utils_devel_code( .false., "DMFT", "dmft_", pub_devel_code)
   pub_dmft_lin_scaling           = utils_devel_code( .false., "DMFT", "dmft_", pub_devel_code)
   pub_dmft_local_scratch         = utils_devel_code( .false., "DMFT", "dmft_", pub_devel_code)
   pub_dmft_norm_proj             = utils_devel_code(       0, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_optics_x1             = utils_devel_code(  0.0_DP, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_optics_y1             = utils_devel_code(  0.0_DP, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_optics_z1             = utils_devel_code(  0.0_DP, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_plot_all_proj         = utils_devel_code( .false., "DMFT", "dmft_", pub_devel_code)
   pub_dmft_plot_real_space_sigma = utils_devel_code( .false., "DMFT", "dmft_", pub_devel_code)
   pub_dmft_scaling_cutoff_h      = utils_devel_code(1.E-5_DP, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_scaling_iter          = utils_devel_code(    2000, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_scaling_maxspace      = utils_devel_code(      20, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_scaling_tol           = utils_devel_code(1.E-8_DP, "DMFT", "dmft_", pub_devel_code)
   pub_dmft_split                 = utils_devel_code( .false., "DMFT", "dmft_", pub_devel_code)
   pub_dmft_splitk                = utils_devel_code( .false., "DMFT", "dmft_", pub_devel_code)

 end subroutine

 subroutine dmft_finalize_kernel(denskern, denskern_tmp, ham, overlap, &
      inv_overlap, rep, num_spins, shifted, wont_compute_residues, &
      dmft_energy_cor, dmft_kernel)

   !---------------------------------------------------------------------------!
   ! Computes the density kernel from the trace of the Green's function over   !
   ! Matsubara frequencies                                                     !
   !---------------------------------------------------------------------------!
   ! NOTES                                                                     !
   ! - To eventually be depreciated if we can adopt standard kernel_diis       !
   !   routines                                                                !
   ! - TODO: depreciate kernel_info file and ientry flag; instead only         !
   !   allocate arrays up to ientry rather than pub_kernel_diis_size           !
   !---------------------------------------------------------------------------!

   use restart,             only: restart_kernel_read, restart_kernel_write
   use sparse, only: sparse_axpy, sparse_copy, sparse_create, sparse_destroy, &
        sparse_read, sparse_scale, sparse_trace, sparse_write
   use sparse_array, only: SPAM3_ARRAY, sparse_array_create, &
        sparse_array_destroy, sparse_array_write, sparse_array_read
   use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_array_create, SPAM3_EMBED
   use kernel,              only: DKERN, kernel_purify
   use utils,               only: utils_unit, utils_banner, utils_assert
   use comms,               only: pub_on_root
   use rundat, only: pub_debug, pub_debug_on_root, pub_dmft_kernel, &
        pub_dmft_kernel_mix, pub_dmft_paramagnetic, pub_dmft_read, &
        pub_dmft_temp, pub_dmft_write, pub_kernel_diis, pub_kernel_diis_coeff, &
        pub_kernel_diis_size, pub_num_spins, pub_rootname, pub_write_denskern
   use ngwf_representation, only: NGWF_REP

   implicit none

   ! Arguments
   type(DKERN),               intent(inout) :: denskern
   type(DKERN),               intent(inout) :: denskern_tmp
   type(SPAM3),               intent(in   ) :: ham(pub_num_spins)
   type(SPAM3_EMBED),               intent(in   ) :: overlap
   type(SPAM3_EMBED),               intent(in   ) :: inv_overlap
   type(NGWF_REP),            intent(in   ) :: rep
   integer,                   intent(in   ) :: num_spins
   logical,                   intent(inout) :: shifted
   logical,                   intent(in   ) :: wont_compute_residues
   real(kind = DP), optional, intent(in   ) :: dmft_energy_cor
   type(SPAM3), optional,     intent(inout) :: dmft_kernel(pub_num_spins)

   ! Local variables
   real(kind = DP)   :: dmft_kernel_diis_c_out
   type(SPAM3_ARRAY) :: dkn_in
   type(SPAM3_ARRAY) :: dkn_out
   type(SPAM3_ARRAY) :: residues
   type(SPAM3)       :: next_dkn_in(pub_num_spins)
   integer           :: ientry
   integer           :: iterdiis
   integer           :: sizdiis
   logical           :: checkdiis
   type(SPAM3_EMBED_ARRAY) :: pur_denskern_tmp
   integer           :: funit
   integer           :: is
   integer           :: i

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_finalize_kernel'

   dmft_kernel_diis_c_out = 1.0_DP - pub_kernel_diis_coeff

   if (pub_dmft_kernel /= 0) call sparse_embed_array_create(pur_denskern_tmp, &
        denskern%kern, iscmplx = .false.)

   if( pub_dmft_temp > 0.0_DP .and. (.not. wont_compute_residues .or. &
        pub_dmft_mpi_rank == 1) ) then

      if(pub_dmft_kernel == -1 .and. .not.present(dmft_energy_cor) .and. &
           .not.present(dmft_kernel))then
         call kernel_purify(pur_denskern_tmp, denskern_tmp, overlap, &
              inv_overlap, rep%n_occ)
         do is = 1, pub_num_spins
            write(stdout, *) 'Trace PURE-DensKern(is) Overlap : ', &
                 sparse_trace(pur_denskern_tmp%m(is, PUB_1K)%p, overlap%p)
            write(stdout, *) 'Trace PURE-DensKern(is) HAM(is) : ', &
                 sparse_trace(pur_denskern_tmp%m(is, PUB_1K)%p, ham(is))
            call sparse_copy(denskern_tmp%kern%m(is, PUB_1K)%p, &
                 pur_denskern_tmp%m(is, PUB_1K)%p)
         enddo
      endif

      if( (pub_dmft_kernel /= 0 .or. present(dmft_kernel)) .and. &
           .not.present(dmft_energy_cor) )then
         if(pub_dmft_mpi_rank == 1 .and. pub_on_root) write(stdout, "(a)") &
              utils_banner("-", 'Kernel mixing')
         if(.not.present(dmft_kernel)) call restart_kernel_read(denskern%kern)

         !BUG September 2014, to impose paramagnetic calculation
         if(pub_num_spins == 2 .and. pub_dmft_paramagnetic)then
            call sparse_copy(denskern%kern%m(2, PUB_1K)%p, denskern%kern%m(1, &
                 PUB_1K)%p)
         endif
         !END BUG

         if(present(dmft_kernel) .or. pub_dmft_kernel /= 2)then
            if(pub_dmft_mpi_rank == 1 .and. pub_on_root) write(stdout, "(a)") &
                 'Performing linear mixing'
            do is = 1, pub_num_spins
               call sparse_scale(denskern_tmp%kern%m(is, PUB_1K)%p, &
                    pub_dmft_kernel_mix)
               if(pub_num_spins == 1 .and. present(dmft_kernel))then
                  !when DMFT is called in subroutine Hamiltonian_build_dens,
                  ! the kernel has a factor x2
                  call sparse_axpy(denskern_tmp%kern%m(is, PUB_1K)%p, &
                       denskern%kern%m(is, PUB_1K)%p, &
                       0.5_DP*(1.0_DP-pub_dmft_kernel_mix))
                  call sparse_scale(denskern_tmp%kern%m(is, PUB_1K)%p, 2.0_DP)
               else
                  call sparse_axpy(denskern_tmp%kern%m(is, PUB_1K)%p, &
                       denskern%kern%m(is, PUB_1K)%p, 1.0_DP-pub_dmft_kernel_mix)
               endif
            enddo
         endif

         !bug was v1
         !#ifndef debugMIXv1
         !          sizdiis = min(iterdiis, pub_kernel_diis_size)
         !#else
         sizdiis = pub_kernel_diis_size
         !#endif
         if(.not.present(dmft_kernel) .and. pub_dmft_kernel == 2)then

            call utils_assert(pub_dmft_write, "Kernel mixing relies on writing &
                 &kernels to file. Please turn dmft_write to TRUE")

            call utils_assert(pub_dmft_read, "Kernel mixing relies on reading &
                 &kernels from file. Please turn dmft_read to TRUE")

            pub_kernel_diis_coeff  = 1.0_DP - pub_dmft_kernel_mix
            dmft_kernel_diis_c_out =       pub_dmft_kernel_mix
            if(pub_dmft_mpi_rank == 1 .and. pub_on_root) write(stdout, "(a)") &
                 'Performing Pulay mixing'
            ientry = 1
            shifted  = .false.
            iterdiis = 1
            funit = utils_unit()

            ! Create arrays
            call sparse_array_create(dkn_in, n_kpoints = pub_kernel_diis_size, &
                 n_spins = pub_num_spins, structure = 'K')
            call sparse_array_create(dkn_out, n_kpoints = pub_kernel_diis_size, &
                 n_spins = pub_num_spins, structure = 'K')
            call sparse_array_create(residues, n_kpoints = pub_kernel_diis_size, &
                 n_spins = pub_num_spins, structure = 'K')
            do is = 1, pub_num_spins
               next_dkn_in%structure = 'K'
               call sparse_create(next_dkn_in(is))
            end do
            ! call dmft_kernel_diis_sparse_create(dkn_in, dkn_out, next_dkn_in, &
            !      residues, denskern)
            inquire(file = trim(pub_rootname) // '.diis_kernel_info', exist = &
                 checkdiis)
            if(checkdiis)then
               open(file = trim(pub_rootname) // '.diis_kernel_info', unit = &
                    funit, form = 'unformatted')
               read(funit) ientry, shifted, iterdiis
               close(funit)
               if (pub_dmft_mpi_rank == 1) then
                  call sparse_array_read(dkn_in, trim(pub_rootname) // &
                       '.diis_kernel_in')
                  call sparse_array_read(dkn_out, trim(pub_rootname) // &
                       '.diis_kernel_out')
                  call sparse_array_read(residues, trim(pub_rootname) // &
                       '.diis_kernel_res')
               end if
            else
               call dmft_kernel_diis_init(residues, next_dkn_in, dkn_out, &
                    dkn_in, denskern)
            endif

            do is = 1, pub_num_spins
               call sparse_copy(dkn_out%m(is, ientry), denskern_tmp%kern%m(is, &
                    PUB_1K)%p)
            enddo

            call dmft_kernel_diis_residue_inout(residues%m(:, ientry), &
                 dkn_out%m(:, ientry), dkn_in%m(:, ientry))
            call dmft_kernel_diis_mix(next_dkn_in, dkn_in%m, dkn_out%m, residues%m, &
                 overlap%p, iterdiis, ientry, num_spins, wont_compute_residues)
            call dmft_kernel_diis_shift(shifted, dkn_in%m, dkn_out%m, residues%m, &
                 iterdiis)
            call dmft_kernel_diis_find_ientry(ientry, iterdiis, shifted)
            iterdiis = iterdiis + 1
            do is = 1, pub_num_spins
               call sparse_copy(dkn_in%m(is, ientry), next_dkn_in(is))
               call sparse_copy(denskern_tmp%kern%m(is, PUB_1K)%p, &
                    next_dkn_in(is))
            end do

            if(pub_dmft_mpi_rank == 1 .and. pub_on_root) then
               open(file = trim(pub_rootname) // '.updated_diis_kernel', unit = &
                    funit, form = 'unformatted')
               write(funit) ientry, shifted, iterdiis
               close(funit)
            endif

            if(pub_dmft_mpi_rank == 1)then
               if (pub_dmft_mpi_rank == 1) then
                  call sparse_array_write(dkn_in, trim(pub_rootname) // &
                       '.updated_diis_kernel_in')
                  call sparse_array_write(dkn_out, trim(pub_rootname) // &
                       '.updated_diis_kernel_out')
                  call sparse_array_write(residues, trim(pub_rootname) // &
                       '.updated_diis_kernel_res')
               end if
            endif

            call sparse_array_destroy(dkn_in)
            call sparse_array_destroy(dkn_out)
            call sparse_array_destroy(residues)
            do is = 1, pub_num_spins
               call sparse_destroy(next_dkn_in(is))
            end do

            ! call dmft_kernel_diis_sparse_destroy(dkn_in, dkn_out, next_dkn_in, &
            !      residues)
         endif

         if(pub_dmft_mpi_rank == 1 .and. .not.present(dmft_kernel) .and. &
            pub_write_denskern) then
            call restart_kernel_write(denskern_tmp%kern, write_dmft = .true.)
         end if
         ! if(.not.present(dmft_kernel) .and. pub_on_root) write(stdout, "(a/)") &
         !      'WARNING: replacing DFT kernel with DMFT kernel'
         if(.not.present(dmft_kernel))then
            do is = 1, pub_num_spins
               call sparse_copy(denskern%kern%m(is, PUB_1K)%p, &
                    denskern_tmp%kern%m(is, PUB_1K)%p)
            enddo
         else
            do is = 1, pub_num_spins
               call sparse_copy(dmft_kernel(is), denskern%kern%m(is, PUB_1K)%p)
               call sparse_copy(denskern%kern%m(is, PUB_1K)%p, &
                    denskern_tmp%kern%m(is, PUB_1K)%p)
            enddo
         endif
         if(pub_dmft_mpi_rank == 1 .and. pub_on_root) write(stdout, "(a/)") &
              repeat("-", 80)
      endif
   endif

   !    if (pub_dmft_kernel /= 0) call sparse_array_destroy(pur_denskern_tmp)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_finalize_kernel'

 end subroutine

 subroutine dmft_wrapper_add_self_and_invert_green(greenf, self_energy, &
      ngwf_basis, ham, overlap, overlap_cplx, inv_overlap, energy, fermi_e, &
      ien, is, ikpoint, nkpoints, norm_kpoints, ffirst, kp_file_shift, &
      splitk_first, kpoint_data, split_Sm, split_H, matrixSigStat, nplanes, &
      orb_in_plane, one_k_one_file)

   !---------------------------------------------------------------------------!
   ! Add the self energy to the Green's function constructor                   !
   !                           omega - H - Sigma                               !
   ! and then invert to obtain the Green's function                            !
   !---------------------------------------------------------------------------!

   use sparse, only: sparse_axpy, sparse_convert, sparse_index_length, &
        sparse_num_cols
   use linalg,         only: linalg_invert_sym_cmatrix
   use comms,          only: pub_on_root
   use rundat, only: pub_debug, pub_debug_on_root, &
        pub_dmft_points, pub_num_spins
   use function_basis, only: FUNC_BASIS

   implicit none

   ! Arguments
   type(DMFT_MATRIX),               intent(inout) :: greenf
   type(DMFT_MATRIX),               intent(inout) :: self_energy
   type(FUNC_BASIS),                intent(in   ) :: ngwf_basis
   type(SPAM3),                     intent(in   ) :: ham(pub_num_spins)
   type(SPAM3),                     intent(in   ) :: overlap
   type(SPAM3),                     intent(in   ) :: overlap_cplx
   type(SPAM3),                     intent(in   ) :: inv_overlap
   complex(kind = DP),              intent(in   ) :: energy
   real(kind = DP),                 intent(in   ) :: fermi_e(2)
   integer,                         intent(in   ) :: ien
   integer,                         intent(in   ) :: is
   integer,                         intent(in   ) :: ikpoint
   integer,                         intent(in   ) :: nkpoints
   real(kind = DP), allocatable,    intent(inout) :: norm_kpoints(:)
   integer,                         intent(in   ) :: ffirst
   integer,                         intent(in   ) :: kp_file_shift
   integer,                         intent(in   ) :: splitk_first
   type(DMFT_KPOINT_MATRICES),      intent(in   ) :: kpoint_data
   complex(kind = DP), allocatable, intent(inout) :: split_Sm(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_H(:,:)
   type(DEM), allocatable,          intent(inout) :: matrixSigStat(:)
   integer,                         intent(in   ) :: nplanes
   logical, allocatable,            intent(inout) :: orb_in_plane(:,:)
   logical,                         intent(in   ) :: one_k_one_file

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_wrapper_add_self_and_invert_green'

   if (.not.pub_dmft_split) then

      ! Subtract the self energy from G^{-1}
      if (nkpoints == 1) then
         call sparse_axpy(greenf%inv_mat, self_energy%mat, (-1.0_DP, 0.0_DP))
      else
         call sparse_axpy(greenf%inv_mat, kpoint_data%self_energy, &
              (-1.0_DP, 0.0_DP))
      endif
      ! Invert to obtain the Green's function
      call dmft_invert_green_function(greenf, self_energy, ngwf_basis, ham, &
           overlap, overlap_cplx, inv_overlap, energy, fermi_e, ien, is, &
           ikpoint, nkpoints, norm_kpoints, ffirst, kp_file_shift, &
           splitk_first, kpoint_data, matrixSigStat, nplanes, orb_in_plane, &
           one_k_one_file)

   else
      if (pub_dmft_lin_scaling) then
         if (ien /= pub_dmft_points-3 .and. ien /= pub_dmft_points-4 .and. ien &
              /= pub_dmft_points-5 .and. ien /= pub_dmft_points-6)then

            call internal_invert_linear_scaling(greenf, energy, overlap_cplx, &
                 ham(is), self_energy%mat, is, fermi_e, ngwf_basis%num, &
                 self_energy%backup, matrixSigStat, update = .false., &
                 realorcomplexH = .true.)
         else

            if(ien == pub_dmft_points-4)then
               greenf%backup = MATMUL(split_Sm, MATMUL(-split_H, split_Sm))
            elseif(ien == pub_dmft_points-3 .or. ien == pub_dmft_points-6)then
               greenf%backup = split_Sm
            elseif(ien == pub_dmft_points-5)then
               greenf%backup = MATMUL(split_Sm, MATMUL(-self_energy%backup, &
                    split_Sm))
            endif

         endif
      else
         if(ien /= pub_dmft_points-3 .and. ien /= pub_dmft_points-4 .and. ien &
              /= pub_dmft_points-5 .and. ien /= pub_dmft_points-6)then

            greenf%backup = greenf%backup - self_energy%backup

            call linalg_invert_sym_cmatrix(greenf%backup, ngwf_basis%num)
            ! CWTODO: why was this loop here? It doesn't appear to achieve
            ! anything
            ! if(nkpoints == 1)then
            !    call linalg_invert_sym_cmatrix(greenf%backup, ngwf_basis%num)
            ! else
            !    call invert_gen_cmat(ngwfbasis%num, greenf%backup)
            ! endif

         else

            if(ien == pub_dmft_points-4)then
               greenf%backup = MATMUL(split_Sm, MATMUL(-split_H, split_Sm))
            elseif(ien == pub_dmft_points-3 .or. ien == pub_dmft_points-6)then
               greenf%backup = split_Sm
            elseif(ien == pub_dmft_points-5)then
               greenf%backup = MATMUL(split_Sm, MATMUL(-self_energy%backup, &
                    split_Sm))
            endif

         endif
      endif
   endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_wrapper_add_self_and_invert_green'

 end subroutine

 subroutine dmft_wrapper_build_self_energy(self_energy, h_atoms, hub_overlap, &
      hub, hub_proj_basis, is, ien, iien, nkpoints, hub_atom, hub_atomb, &
      cluster, connection_table, nmerge, merge_table, merge_neigh, &
      dimer_table, kpoint_data, kstart, kstep, ffirst, llast, rot_vec_angles, &
      same_self_for_all, show_matrices, force_reopen_file, &
      reopen_sigma_file_each_step, scratchdir)

   !---------------------------------------------------------------------------!
   ! Upfolds the local self energy to the NGWF subspace                        !
   !---------------------------------------------------------------------------!

   use comms,             only: pub_my_proc_id, pub_on_root
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par
   use rundat,            only: pub_debug, pub_debug_on_root
   use sparse,            only: sparse_convert

   implicit none

   ! Arguments
   type(DMFT_MATRIX),            intent(inout) :: self_energy
   type(ARRAY_OF_MATRICES),      intent(in   ) :: h_atoms(par%nat_hub)
   type(DMFT_OVERLAP),           intent(in   ) :: hub_overlap
   type(HUBBARD_MODEL),          intent(in   ) :: hub
   type(FUNC_BASIS),             intent(in   ) :: hub_proj_basis
   integer,                      intent(in   ) :: is
   integer,                      intent(in   ) :: ien
   integer,                      intent(in   ) :: iien
   integer,                      intent(in   ) :: nkpoints
   integer,                      intent(inout) :: hub_atom
   integer,                      intent(inout) :: hub_atomb
   logical,                      intent(in   ) :: cluster
   logical, allocatable,         intent(in   ) :: connection_table(:,:)
   integer,                      intent(in   ) :: nmerge
   integer, allocatable,         intent(in   ) :: merge_neigh(:)
   integer, allocatable,         intent(in   ) :: merge_table(:,:)
   integer, allocatable,         intent(in   ) :: dimer_table(:,:)
   type(DMFT_KPOINT_MATRICES),   intent(inout) :: kpoint_data
   integer,                      intent(in   ) :: kstart
   integer,                      intent(in   ) :: kstep
   integer,                      intent(in   ) :: ffirst
   integer,                      intent(in   ) :: llast
   real(kind = DP), allocatable, intent(in   ) :: rot_vec_angles(:,:,:)
   logical,                      intent(in   ) :: same_self_for_all
   logical,                      intent(in   ) :: show_matrices
   logical,                      intent(in   ) :: force_reopen_file
   logical,                      intent(in   ) :: reopen_sigma_file_each_step
   character(2000),              intent(in   ) :: scratchdir

   ! Local variables
   integer :: kien
   integer :: kkien
   integer :: channels

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_wrapper_build_self_energy'

   if (.not.pub_dmft_split) then
      kien = ien
      call dmft_construct_self_energy(ien, kien, is, self_energy, hub_atom, &
           hub_atomb, connection_table, nmerge, merge_table, merge_neigh, &
           dimer_table, h_atoms, hub, hub_proj_basis, rot_vec_angles, ffirst, &
           llast, cluster, show_matrices, reopen_sigma_file_each_step, &
           force_reopen_file, scratchdir)
      call dmft_upfold_self_energy(self_energy, hub_overlap, nkpoints, &
           kpoint_data, same_self_for_all)

   else

      do kkien = 1, kstep

         kien = iien - kstart + (kkien-1)
         call dmft_construct_self_energy(ien, kien, is, self_energy, hub_atom, &
              hub_atomb, connection_table, nmerge, merge_table, merge_neigh, &
              dimer_table, h_atoms, hub, hub_proj_basis, rot_vec_angles, &
              ffirst, llast, cluster, show_matrices, &
              reopen_sigma_file_each_step, force_reopen_file, scratchdir)
         call dmft_upfold_self_energy(self_energy, hub_overlap, nkpoints, &
              kpoint_data, same_self_for_all)

         if (pub_my_proc_id == kkien-1) then
            if (nkpoints == 1) then
               call sparse_convert(self_energy%backup, self_energy%mat)
            else
               call sparse_convert(self_energy%backup, kpoint_data%self_energy)
            endif
         endif

      enddo

   endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_wrapper_build_self_energy'

 end subroutine

 subroutine dmft_build_inv_green_matrix(greenf, overlap_cplx, ham, ham_kpoint, &
      ham_embed, is, ien, nkpoints, energy, check_embed_h, split_S, split_H)

   !---------------------------------------------------------------------------!
   ! Builds the inverted impurity Green's function                             !
   !                      G^-1(E) = E S - H - (Self - Edc)                     !
   !---------------------------------------------------------------------------!

   use sparse, only: sparse_axpy, sparse_copy, sparse_scale, sparse_show_matrix
   use utils,  only: utils_abort, utils_int_to_str
   use rundat, only: pub_debug, pub_debug_on_root, pub_dmft_points, &
        pub_num_spins

   implicit none

   ! Arguments
   type(DMFT_MATRIX),               intent(inout) :: greenf
   type(SPAM3),                     intent(in   ) :: overlap_cplx
   type(SPAM3),                     intent(in   ) :: ham(pub_num_spins)
   type(SPAM3),                     intent(in   ) :: ham_kpoint
   type(SPAM3),                     intent(in   ) :: ham_embed(pub_num_spins)
   integer,                         intent(in   ) :: is
   integer,                         intent(in   ) :: ien
   integer,                         intent(in   ) :: nkpoints
   complex(kind = DP),              intent(in   ) :: energy
   logical,                         intent(in   ) :: check_embed_h
   complex(kind = DP), allocatable, intent(in   ) :: split_S(:,:)
   complex(kind = DP), allocatable, intent(in   ) :: split_H(:,:)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_build_inv_green_matrix'

   if(.not.pub_dmft_split) call sparse_copy(greenf%inv_mat, overlap_cplx)

   if(ien /= pub_dmft_points-3 .and. ien /= pub_dmft_points-6)then

      if(.not.pub_dmft_split)then
         call sparse_scale(greenf%inv_mat, energy)
         if(ien /= pub_dmft_points-5) then
            if(nkpoints == 1)then
               call sparse_axpy(greenf%inv_mat, ham(is), (-1.0_DP, 0.0_DP))
            else
               call sparse_axpy(greenf%inv_mat, ham_kpoint, (-1.0_DP, 0.0_DP))
            endif
            if(check_embed_h .and. ien /= pub_dmft_points-4)then
               call sparse_axpy(greenf%inv_mat, ham_embed(is), &
                    (-1.0_DP, 0.0_DP))
            endif
         endif
      else
         if(ien /= pub_dmft_points-5 .and. ien /= pub_dmft_points-4) then
            greenf%backup = energy*split_S-split_H
         end if
         if(check_embed_h)then
            call utils_abort('Error in dmft_build_inv_green_matrix: Embedding &
                 &and split are not compatible')
         endif
      endif

   endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_build_inv_green_matrix'

 end subroutine

 subroutine dmft_set_kpoints(kpt, norm_kpoints, nkpoints, lastkpoint, &
      totkpoint, mdl)

   !---------------------------------------------------------------------------!
   ! Selects the k-points in order to make the Gamma point calculations        !
   ! periodic (under the tight binding approximation). The k-points are        !
   ! chosen uniformly in the Brillouin zone.                                   !`
   !---------------------------------------------------------------------------!

   use utils,      only: utils_abort, utils_alloc_check, utils_dealloc_check
   use comms,      only: pub_on_root
   use model_type, only: MODEL
   use rundat, only: pub_debug, pub_debug_on_root, pub_dmft_kernel, &
        pub_dmft_kpoints_sym, pub_dmft_nkpoints

   implicit none

   ! Arguments
   real(kind = DP), allocatable, intent(inout) :: norm_kpoints(:)
   real(kind = DP), allocatable, intent(inout) :: kpt(:,:)
   integer,                      intent(  out) :: nkpoints
   integer,                      intent(  out) :: lastkpoint
   real(kind = DP),              intent(  out) :: totkpoint
   type(MODEL),                  intent(in   ) :: mdl

   ! Local variables
   integer         :: jp, ii1, ii2, ii3, i, j
   real(kind = DP) :: t1_, t2_, t3_, norm
   integer         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_set_kpoints'

   jp = abs(pub_dmft_nkpoints)
   if (jp > 1) then

      nkpoints = 0

      if(pub_dmft_nkpoints < 0)then
         do ii1 = -(jp-1)/2, (jp-1)/2
            do ii2 = -(jp-1)/2, (jp-1)/2
               do ii3 = -(jp-1)/2, (jp-1)/2
                  if(ii3 > 0 .or. (ii3 == 0 .and. ((ii2 > 0) .or. (ii2 == 0 &
                       .and. ii1 >= 0)))) then
                     nkpoints = nkpoints + 1
                  endif
               enddo
            enddo
         enddo
         if(allocated(norm_kpoints)) then
            deallocate(norm_kpoints, stat = ierr)
            call utils_dealloc_check('dmft_set_kpoints', 'norm_kpoints', ierr)
            deallocate(kpt, stat = ierr)
            call utils_dealloc_check('dmft_set_kpoints', 'kpt', ierr)
         end if
         allocate(norm_kpoints(nkpoints), stat = ierr)
         call utils_alloc_check('dmft_set_kpoints', 'norm_kpoints', ierr)
         allocate(kpt(nkpoints, 3), stat = ierr)
         call utils_alloc_check('dmft_set_kpoints', 'kpt', ierr)
         j = 0
         do ii1 = -(jp-1)/2, (jp-1)/2
            do ii2 = -(jp-1)/2, (jp-1)/2
               do ii3 = -(jp-1)/2, (jp-1)/2
                  if( ii3 > 0 .or. (ii3 == 0 .and. ( (ii2 > 0) .or. (ii2 == 0 &
                       .and. ii1 >= 0) ) ) )then
                     j = j + 1
                     t1_ = real(ii1, kind=DP)/real(jp, kind=DP)
                     t2_ = real(ii2, kind=DP)/real(jp, kind=DP)
                     t3_ = real(ii3, kind=DP)/real(jp, kind=DP)
                     kpt(j, 1) = t1_*mdl%cell%b1%x + t2_*mdl%cell%b2%x + &
                          t3_*mdl%cell%b3%x
                     kpt(j, 2) = t1_*mdl%cell%b1%y + t2_*mdl%cell%b2%y + &
                          t3_*mdl%cell%b3%y
                     kpt(j, 3) = t1_*mdl%cell%b1%z + t2_*mdl%cell%b2%z + &
                          t3_*mdl%cell%b3%z
                  endif
               enddo
            enddo
         enddo
         do i = 1, nkpoints
            norm = 2.0_DP
            if(sum(abs(kpt(i, :)**2)) < 1.0E-5_DP) norm = norm/2.0_DP !Gamma Point
            norm_kpoints(i) = norm
         enddo
      else
         nkpoints = jp**3
         if(allocated(norm_kpoints)) then
            deallocate(norm_kpoints, stat = ierr)
            call utils_dealloc_check('dmft_set_kpoints', 'norm_kpoints', ierr)
            deallocate(kpt, stat = ierr)
            call utils_dealloc_check('dmft_set_kpoints', 'kpt', ierr)
         end if
         allocate(norm_kpoints(nkpoints), stat = ierr)
         call utils_alloc_check('dmft_set_kpoints', 'norm_kpoints', ierr)
         allocate(kpt(nkpoints, 3), stat = ierr)
         call utils_alloc_check('dmft_set_kpoints', 'kpt', ierr)
         j = 0
         do ii1 = 1, jp
            do ii2 = 1, jp
               do ii3 = 1, jp
                  j = j + 1
                  t1_ = real(ii1-1, kind=DP)/real(2*jp-1, kind=DP)
                  t2_ = real(ii2-1, kind=DP)/real(2*jp-1, kind=DP)
                  t3_ = real(ii3-1, kind=DP)/real(2*jp-1, kind=DP)
                  kpt(j, 1) = t1_*mdl%cell%b1%x + t2_*mdl%cell%b2%x + &
                       t3_*mdl%cell%b3%x
                  kpt(j, 2) = t1_*mdl%cell%b1%y + t2_*mdl%cell%b2%y + &
                       t3_*mdl%cell%b3%y
                  kpt(j, 3) = t1_*mdl%cell%b1%z + t2_*mdl%cell%b2%z + &
                       t3_*mdl%cell%b3%z
               enddo
            enddo
         enddo
         do i = 1, nkpoints
            norm = 8.0_DP
            do j = 1, 3
               if(abs(kpt(i, j)) < 1.0E-5_DP) norm = norm/2.0_DP !we only keep
               ! mirrors...
            enddo
            norm_kpoints(i) = norm
         enddo
      endif

      lastkpoint = nkpoints

      if(pub_dmft_kpoints_sym .and. (pub_dmft_kernel == 0 .or. &
           pub_dmft_kpoints_kernel_gamma))then
         !WARNING: CUBIC SYM ON KPOINTS NOT COMPATIBLE WITH KERNEL
         do i = 1, nkpoints
            do j = 1, nkpoints
               if(j /= i)then
                  do ii1 = -1, 1, 2
                     do ii2 = -1, 1, 2
                        do ii3 = -1, 1, 2
                           kpt(j, 1) = real(ii1, kind=DP)*kpt(j, 1) ! mirror
                           kpt(j, 2) = real(ii2, kind=DP)*kpt(j, 2) ! mirror
                           kpt(j, 3) = real(ii3, kind=DP)*kpt(j, 3) ! mirror
                           if(internal_norm_vector(kpt(i, :)-kpt(j, :)) < &
                                1.0E-5_DP)then
                              kpt(j, :) = 1.0E15_DP
                              norm_kpoints(i) = norm_kpoints(i) + &
                                   norm_kpoints(j)
                              norm_kpoints(j) = 0.0_DP
                           endif
                           kpt(j, 1) = real(ii1, kind=DP)*kpt(j, 1) ! mirror
                           kpt(j, 2) = real(ii2, kind=DP)*kpt(j, 2) ! mirror
                           kpt(j, 3) = real(ii3, kind=DP)*kpt(j, 3) ! mirror
                        enddo
                     enddo
                  enddo
               endif
            enddo
         enddo
         lastkpoint = 0
         do i = 1, nkpoints
            if(abs(norm_kpoints(i)) > 1.0E-3_DP) lastkpoint = i
         enddo
      endif

      totkpoint = sum(norm_kpoints)
      norm_kpoints = norm_kpoints/totkpoint
      if(pub_on_root) then
         write(stdout, *) 'TOTAL NUMBER OF K POINTS : ', totkpoint
         write(stdout, *) 'TOTAL COMPUTED K POINTS (orthorombic symmetries) : &
              &', nkpoints
         do i = 1, nkpoints
            if(abs(norm_kpoints(i)) > 1.0E-4_DP)then
               write(stdout, *) 'KP - TOT             : ', i, nkpoints
               write(stdout, *) 'WG - TOT : ', norm_kpoints(i), totkpoint
               write(stdout, '(a, 3f10.5)') 'KP        : ', kpt(i, :)
            endif
         enddo
      endif
   else
      nkpoints = 1
      lastkpoint = 1
      if(allocated(norm_kpoints)) then
         deallocate(norm_kpoints, stat=ierr)
         call utils_dealloc_check('dmft_set_kpoints', 'norm_kpoints', ierr)
         deallocate(kpt, stat=ierr)
         call utils_dealloc_check('dmft_set_kpoints', 'kpt', ierr)
      end if
      allocate(norm_kpoints(nkpoints), stat = ierr)
      call utils_alloc_check('dmft_set_kpoints', 'norm_kpoints', ierr)
      allocate(kpt(nkpoints, 3), stat = ierr)
      call utils_alloc_check('dmft_set_kpoints', 'kpt', ierr)
      norm_kpoints = 1.0_DP
      kpt = 0.0_DP
   endif

   if(nkpoints == 1)then
      if(pub_dmft_splitk) then
         call utils_abort("Error in dmft_set_kpoints: pub_dmft_splitk = True &
              &and nkpoints = 1 is incompatible")
      endif
   endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving dmft_set_kpoints'

 end subroutine

 !=============================================================================!
 ! EMBEDDING ROUTINES                                                          !
 ! Embedding merges previous bulk ONETEP calculations with the current         !
 ! calculation (useful for problems such as the deposition of a                !
 ! nano-structure on a surface)                                                !
 !=============================================================================!

 subroutine dmft_embedding_setup_planes(orb_in_plane, atom_in_plane, nplanes, &
      ngwf_basis, nat)

   !---------------------------------------------------------------------------!
   ! Defines the planes which interface the ONETEP calculations with a surface !
   ! or other interface                                                        !
   !---------------------------------------------------------------------------!

   use comms,             only: pub_on_root
   use function_basis,    only: FUNC_BASIS
   use parallel_strategy, only: par => pub_par
   use rundat,            only: pub_debug, pub_debug_on_root
   use utils,             only: utils_unit,  utils_alloc_check, &
        utils_dealloc_check

   implicit none

   ! Arguments
   logical, allocatable, intent(inout) :: orb_in_plane(:,:)
   logical, allocatable, intent(inout) :: atom_in_plane(:,:)
   integer,              intent(  out) :: nplanes
   integer,              intent(in   ) :: nat
   type(FUNC_BASIS),     intent(in   ) :: ngwf_basis

   ! Local variables
   integer :: nn, ufii, i, j, k, ii, jj, kk, iii, jjj, lll
   logical :: planes_check
   integer :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_embedding_setup_planes'

   nplanes = 0

   inquire(file = 'group_of_atoms_label', exist = planes_check)

   if(.not.planes_check) then
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_embedding_setup_planes as group_of_atoms_label not found'
      return
   end if

   ufii = utils_unit()
   open(unit = ufii, file = 'group_of_atoms_label')
   nn = lines_in_file(ufii)-1
   read(ufii, *) nplanes

   if(allocated(orb_in_plane)) then
      deallocate(orb_in_plane, stat = ierr)
      call utils_dealloc_check('dmft_embedding_setup_planes', 'orb_in_plane', ierr)
      allocate(orb_in_plane(nplanes, ngwf_basis%num), stat = ierr)
      call utils_alloc_check('dmft_embedding_setup_planes', 'orb_in_plane', ierr)
   end if

   if(allocated(atom_in_plane)) then
      deallocate(atom_in_plane, stat = ierr)
      call utils_dealloc_check('dmft_embedding_setup_planes', 'atom_in_plane', ierr)
      allocate(atom_in_plane(nplanes, nat), stat = ierr)
      call utils_alloc_check('dmft_embedding_setup_planes', 'atom_in_plane', ierr)
   end if

   orb_in_plane = .false.
   atom_in_plane = .false.

   do i = 1, nn
      read(ufii, *) k, j
      atom_in_plane(j, k) = .true. ! cw: atom_in_plane is in ORIGINAL atom
      ! notations - for popn !
      jj = par%distr_atom(k)
      jjj = 0
      do ii = 1, jj-1
         jjj = jjj + ngwf_basis%num_on_atom(ii)
      enddo
      do ii = 1, ngwf_basis%num_on_atom(jj)
         orb_in_plane(j, jjj + ii) = .true.
      enddo
   enddo
   close(ufii)

   if(pub_dmft_mpi_rank == 1 .and. pub_on_root)then
      open(unit = ufii, file = 'group_of_atoms_label_sorted', form = &
           'unformatted')
      write(ufii) nplanes
      write(ufii) size(atom_in_plane, 1), size(atom_in_plane, 2)
      write(ufii) atom_in_plane
      close(ufii)
   endif


   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_embedding_setup_planes'

 end subroutine

 subroutine dmft_embedding_prepare_tables(maskH0_embed, maskH1_embed, &
      maskH0_embed_atoms, ngwf_basis, ham, overlap, nkpoints, check_embed)

   !---------------------------------------------------------------------------!
   ! Constructing tables required for embedding calculations                   !
   !---------------------------------------------------------------------------!

   use comms,             only: comms_reduce, pub_my_proc_id, pub_on_root
   use function_basis,    only: FUNC_BASIS
   use parallel_strategy, only: par => pub_par
   use rundat, only: pub_debug, pub_debug_on_root, pub_num_spins
   use sparse,            only: sparse_get_element
   use utils,             only: utils_abort, utils_unit, utils_int_to_str, &
                                utils_alloc_check, utils_dealloc_check

   implicit none

   ! Arguments
   integer, allocatable, intent(inout) :: maskH0_embed(:,:,:)
   integer, allocatable, intent(inout) :: maskH1_embed(:,:,:)
   integer, allocatable, intent(inout) :: maskH0_embed_atoms(:,:,:)
   type(FUNC_BASIS),     intent(in   ) :: ngwf_basis
   type(SPAM3),          intent(in   ) :: ham(2)
   type(SPAM3),          intent(in   ) :: overlap
   integer,              intent(in   ) :: nkpoints
   logical,              intent(  out) :: check_embed

   ! Local variables
   integer                      :: i, j, k, l, ufi, ufi1, kk, s1, s2, t1, t2, &
                                   kkk, i1, j1, i0, j0, ii, ii0, jj0, ii1, jj1
   integer                      :: tot_li, isp
   logical                      :: check_embedl, check_embedr
   real(kind = DP), allocatable :: tmpH(:,:), ovH(:,:)
   real(kind = DP)              :: h_element
   integer, allocatable         :: embed_index(:)
   integer, allocatable         :: embed_size(:)
   integer, allocatable         :: embed_pos(:)
   integer, allocatable         :: embed_pos_all(:)
   integer, allocatable         :: maskH1_embed_atoms(:,:,:)
   integer                      :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_embedding_prepare_tables'

   check_embed = .false.
   inquire(file = 'connections_D_L', exist = check_embedl)
   inquire(file = 'connections_D_R', exist = check_embedr)
   check_embed = check_embedl .or. check_embedr

   if(.not.check_embed) then
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_embedding_prepare_tables as connections file not found'
      return
   end if

   if(pub_dmft_split)then
      call utils_abort('Error in dmft_embedding_prepare_tables: split not &
           &compatible with embedding, as embedding requires all atoms to be &
           &on the same proc')
   endif
   if(nkpoints > 1) then
      call utils_abort('Error in dmft_embedding_prepare_tables: embedding and &
           &k-points not yet implemented')
   endif

   if(check_embedl .and. check_embedr)then
      call utils_abort('Error in dmft_embedding_prepare_tables: both right and &
           &left files detected, should only have one or the other')
   endif

   ufi = utils_unit()
   ufi1 = utils_unit()

   if(check_embedl) open(unit = ufi, file = 'connections_D_L')
   if(check_embedr) open(unit = ufi, file = 'connections_D_R')
   if(check_embedl) open(unit = ufi1, file = 'connections_D_L_lead_sites')
   if(check_embedr) open(unit = ufi1, file = 'connections_D_R_lead_sites')

   if(allocated(embed_index)) then
      deallocate(embed_index, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potentials', &
           'embed_index', ierr)
      deallocate(embed_size, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potentials', 'embed_size', &
           ierr)
      deallocate(embed_pos, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potentials', 'embed_pos', &
           ierr)
      deallocate(embed_pos_all, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potentials', &
           'embed_pos_all', ierr)
      deallocate(maskH0_embed, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potentials', &
           'maskH0_embed', ierr)
      deallocate(maskH1_embed, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potentials', &
           'maskH1_embed', ierr)
      deallocate(maskH0_embed_atoms, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potentials', &
           'maskH0_embed_atoms', ierr)
      deallocate(maskH1_embed_atoms, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potentials', &
           'maskH1_embed_atoms', ierr)
   end if

   kk = lines_in_file(ufi1)
   allocate(embed_index(kk), stat = ierr)
   call utils_alloc_check('dmft_embedding_prepare_tables', 'embed_index', ierr)
   allocate(embed_size(kk), stat = ierr)
   call utils_alloc_check('dmft_embedding_prepare_tables', 'embed_size', ierr)
   allocate(embed_pos(kk), stat = ierr)
   call utils_alloc_check('dmft_embedding_prepare_tables', 'embed_pos', ierr)
   allocate(embed_pos_all(par%nat), stat = ierr)
   call utils_alloc_check('dmft_embedding_prepare_tables', 'embed_pos_all', ierr)
   do i = 1, kk
      read(ufi1, *) s2
      embed_index(i) = par%distr_atom(s2)
   enddo

   kkk = 0
   do i = 1, kk
      embed_size(i) = ngwf_basis%num_on_atom( embed_index(i) )
      kkk           = kkk + embed_size(i)
      ! s2 in original positions (as in dat file) notations
      ! reminder: par%orig_atom(iat) gives the original position, it iat is
      ! the NGWFs index
      ! reminder: par%distr_atom is the opposite
   enddo

   do i = 1, kk
      embed_pos(i) = 0
      do j = 1, embed_index(i)-1
         embed_pos(i) = embed_pos(i) + ngwf_basis%num_on_atom( j )
      enddo
   enddo
   do ii = 1, par%nat
      embed_pos_all(ii) = 0
      i = par%distr_atom(ii)
      do j = 1, i-1
         embed_pos_all(ii) = embed_pos_all(ii) + ngwf_basis%num_on_atom( j )
      enddo
   enddo

   allocate(maskH0_embed(kkk, kkk, 2), stat = ierr)
   call utils_alloc_check('dmft_embedding_prepare_tables', 'maskH0_embed', ierr)
   allocate(maskH1_embed(kkk, kkk, 2), stat = ierr)
   call utils_alloc_check('dmft_embedding_prepare_tables', 'maskH1_embed', ierr)
   allocate(maskH0_embed_atoms(kkk, kkk, 2), stat = ierr)
   call utils_alloc_check('dmft_embedding_prepare_tables', 'maskH0_embed_atoms', ierr)
   allocate(maskH1_embed_atoms(kkk, kkk, 2), stat = ierr)
   call utils_alloc_check('dmft_embedding_prepare_tables', 'maskH1_embed_atoms', ierr)

   do i = 1, kk
      do j = 1, kk
         i0 = 0
         if(i > 1) i0 = sum(embed_size(1:i-1))
         j0 = 0
         if(j > 1) j0 = sum(embed_size(1:j-1))
         do i1 = 1, embed_size(i)
            do j1 = 1, embed_size(j)
               maskH0_embed( i0 + i1, j0 + j1, 1 )      = embed_pos(i) + i1
               maskH0_embed( i0 + i1, j0 + j1, 2 )      = embed_pos(j) + j1
               maskH0_embed_atoms(i0 + i1, j0 + j1, 1)  = embed_index(i)
               maskH0_embed_atoms(i0 + i1, j0 + j1, 2)  = embed_index(j)
            enddo
         enddo
      enddo
   enddo

   tot_li = lines_in_file(ufi)
   do ii = 1, tot_li

      read(ufi, *) s1, s2, t1, t2

      if(s1 == 0 .or. s2 == 0 .or. t1 == 0 .or. t2 == 0)then
         write(stdout, '(a)') 'Zero index found:'
         write(stdout, '(a, 4i5)') 's1, s2, t1, t2 : ', s1, s2, t1, t2
         call utils_abort("Error in dmft_embedding_prepare_tables: zero index &
              &detected")
      endif

      ii0 = 0
      if(s1 > 1) ii0 = sum(embed_size(1:s1-1))
      jj0 = 0
      if(s2 > 1) jj0 = sum(embed_size(1:s2-1))

      if(embed_size(s1) /= ngwf_basis%num_on_atom( par%distr_atom(t1)))then
         write(stdout, *) 'BUILD EMBEDDING TABLE ERROR : size mismatch1'
         write(stdout, *) 's1, s2, t1, t2 : ', s1, s2, t1, t2
         write(stdout, *) 'S2 atom index : ', par%orig_atom(embed_index(s2))
         write(stdout, *) 'embed_size(s1) : ', embed_size(s1)
         write(stdout, *) 'ngwf_basis%num_on_atom( par%distr_atom(t1) ) : ', &
              ngwf_basis%num_on_atom( par%distr_atom(t1) )
         call utils_abort("Error in dmft_embedding_prepare_tables: size &
              &mismatch")
      endif

      if(embed_size(s2) /= ngwf_basis%num_on_atom( par%distr_atom(t2)))then
         write(stdout, *) 'BUILD EMBEDDING TABLE ERROR : size mismatch2'
         write(stdout, *) 's1, s2, t1, t2 : ', s1, s2, t1, t2
         write(stdout, *) 'S2 atom index : ', par%orig_atom(embed_index(s2))
         write(stdout, *) 'embed_size(s2) : ', embed_size(s2)
         write(stdout, *) 'NUM_NG(S2 atom index) : ', ngwf_basis%num_on_atom( &
              par%distr_atom(par%orig_atom(embed_index(s2))) )
         write(stdout, *) 'ngwf_basis%num_on_atom( par%distr_atom(t2) ) : ', &
              ngwf_basis%num_on_atom( par%distr_atom(t2))
         call utils_abort("Error in dmft_embedding_prepare_tables: size &
              &mismatch")
      endif

      do i1 = 1, ngwf_basis%num_on_atom( par%distr_atom(t1) )
         do j1 = 1, ngwf_basis%num_on_atom( par%distr_atom(t2)  )
            maskH1_embed(ii0 + i1, jj0 + j1, 1 )       =  embed_pos_all(t1) + i1
            maskH1_embed(ii0 + i1, jj0 + j1, 2 )       =  embed_pos_all(t2) + j1
            maskH1_embed_atoms(ii0 + i1, jj0 + j1, 1 ) = par%distr_atom(t1)
            maskH1_embed_atoms(ii0 + i1, jj0 + j1, 2 ) = par%distr_atom(t2)
         enddo
      enddo
   enddo

   close(ufi)
   close(ufi1)

   allocate(tmpH(kkk, kkk), stat = ierr)
   call utils_alloc_check('dmft_embedding_prepare_tables', 'tmpH', ierr)
   allocate(ovH(kkk, kkk), stat = ierr)
   call utils_alloc_check('dmft_embedding_prepare_tables', 'ovH', ierr)

   do isp = 1, pub_num_spins

      if(pub_dmft_mpi_rank == 1 .and. pub_on_root) open(unit = ufi, file = &
           '_H0_' // trim(adjustl(utils_int_to_str(isp))), form = 'unformatted')

      tmpH = 0.0_DP
      ovH = 0.0_DP

      do i = 1, kkk
         do j = 1, kkk
            i1 = maskH0_embed(i, j, 1)
            j1 = maskH0_embed(i, j, 2)
            ii1 = maskH0_embed_atoms(i, j, 1)
            jj1 = maskH0_embed_atoms(i, j, 2)

            if(par%proc_of_atom(par%orig_atom(ii1)) == pub_my_proc_id .and. &
                 par%proc_of_atom(par%orig_atom(jj1)) == pub_my_proc_id) call &
                 sparse_get_element(h_element, ham(isp), i1, j1)

            tmpH(i, j) = h_element

            if(par%proc_of_atom(par%orig_atom(ii1)) == pub_my_proc_id .and. &
                 par%proc_of_atom(par%orig_atom(jj1)) == pub_my_proc_id) call &
                 sparse_get_element(h_element, overlap, i1, j1)

            ovH(i, j) = h_element
         enddo
      enddo

      call comms_reduce('SUM', tmpH)
      call comms_reduce('SUM', ovH)

      if(pub_dmft_mpi_rank == 1 .and. pub_on_root) then
         write(ufi) kkk
         write(ufi) tmpH
         write(ufi) ovH
      endif

      if(pub_dmft_mpi_rank == 1 .and. pub_on_root) close(ufi)
      if(pub_dmft_mpi_rank == 1 .and. pub_on_root) open(unit = ufi, file = &
           '_H1_' // trim(adjustl(utils_int_to_str(isp))), form = 'unformatted')

      tmpH = 0
      ovH = 0

      do i = 1, kkk
         do j = 1, kkk
            i1 = maskH1_embed(i, j, 1)
            j1 = maskH1_embed(i, j, 2)
            ii1 = maskH1_embed_atoms(i, j, 1)
            jj1 = maskH1_embed_atoms(i, j, 2)

            if(par%proc_of_atom(par%orig_atom(ii1)) == pub_my_proc_id .and. &
                 par%proc_of_atom(par%orig_atom(jj1)) == pub_my_proc_id) call &
                 sparse_get_element(h_element, ham(isp), i1, j1)

            tmpH(i, j) = h_element

            if(par%proc_of_atom(par%orig_atom(ii1)) == pub_my_proc_id .and. &
                 par%proc_of_atom(par%orig_atom(jj1)) == pub_my_proc_id) call &
                 sparse_get_element(h_element, overlap, i1, j1)

            ovH(i, j) = h_element
         enddo
      enddo

      call comms_reduce('SUM', tmpH)
      call comms_reduce('SUM', ovH)

      if(pub_dmft_mpi_rank == 1 .and. pub_on_root) then
         write(ufi) kkk
         write(ufi) tmpH
         write(ufi) ovH
      endif

      if(pub_dmft_mpi_rank == 1 .and. pub_on_root) close(ufi)

   enddo

   deallocate(tmpH, stat = ierr)
   call utils_dealloc_check('dmft_embedding_prepare_tables', 'tmpH', ierr)
   deallocate(ovH, stat = ierr)
   call utils_dealloc_check('dmft_embedding_prepare_tables', 'ovH', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_embedding_prepare_tables'

 end subroutine

 subroutine dmft_embedding_write_sigma(check_embed, maskH0_embed, &
      maskH0_embed_atoms, self_energy, energy, fermi_e, is, ien, ffirst)

   !---------------------------------------------------------------------------!
   ! Prepares the embedding potential for subsequent embedding calculations    !
   !---------------------------------------------------------------------------!

   use sparse,            only: sparse_get_element
   use comms,             only: comms_reduce, pub_my_proc_id, pub_on_root
   use utils,             only: utils_unit, utils_int_to_str, &
                                utils_alloc_check, utils_dealloc_check
   use parallel_strategy, only: par => pub_par
   use rundat, only: pub_debug, pub_debug_on_root, &
        pub_dmft_points, pub_dmft_temp

   implicit none

   ! Arguments
   logical,              intent(in   ) :: check_embed
   integer, allocatable, intent(in   ) :: maskH0_embed(:,:,:)
   integer, allocatable, intent(in   ) :: maskH0_embed_atoms(:,:,:)
   type(DMFT_MATRIX),    intent(in   ) :: self_energy
   complex(kind = DP),   intent(in   ) :: energy
   real(kind = DP),      intent(in   ) :: fermi_e(2)
   integer,              intent(in   ) :: is
   integer,              intent(in   ) :: ien
   integer,              intent(in   ) :: ffirst

   ! Local variables
   complex(kind = DP), allocatable :: tmp_mat(:,:)
   integer                         :: size_mat0, size_mat1, i, j, s, t, ii, jj
   integer                         :: funit
   complex(kind = DP)              :: sig_element
   character(len = 100)            :: filename
   integer                         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_embedding_write_sigma'

   if(.not.check_embed) return

   if(pub_on_root .and. pub_dmft_mpi_rank == 1) write(stdout, *) 'WRITING &
        &EMBEDDING POTENTIALS'

   allocate(tmp_mat(size(maskH0_embed, 1), size(maskH0_embed, 1)), stat = ierr)
   call utils_alloc_check('dmft_embedding_write_sigma', 'tmp_mat', ierr)

   ! ebl: eliminating this check, as maskH0_embed and maskH1_embed are currently
   ! allocated at exactly the same time with exactly the same shape.
   ! maskH1_embed is therefore not needed in this routine at all, although I
   ! suspect it will be at a later date, so I have kept it accessible to this
   ! routine

   size_mat0 = size(maskH0_embed, 1)
   ! size_mat1 = size(maskH1_embed, 1)
   ! if(size_mat0 /= size_mat1)then
   !    write(stdout, *) 'ERROR : embedding, size_mat do not match'
   !    stop
   ! endif

   tmp_mat = 0.0_DP
   do i = 1, size_mat0
      do j = 1, size_mat0
         s  = maskH0_embed(i, j, 1)
         t  = maskH0_embed(i, j, 2)
         ii = maskH0_embed_atoms(i, j, 1)
         jj = maskH0_embed_atoms(i, j, 2)
         if(par%proc_of_atom(par%orig_atom(ii)) == pub_my_proc_id .and. &
              par%proc_of_atom(par%orig_atom(jj)) == pub_my_proc_id) call &
              sparse_get_element(sig_element, self_energy%mat, s, t)
         tmp_mat(i, j) = sig_element
      enddo
   enddo
   call comms_reduce('SUM', tmp_mat)

   if(pub_my_proc_id /= 0) then
      deallocate(tmp_mat, stat = ierr)
      call utils_dealloc_check('dmft_embedding_write_sigma', 'tmp_mat', ierr)
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_embedding_write_sigma'
      return
   endif

   filename = 'embedding_potentials_' // TRIM(ADJUSTL(utils_int_to_str(is)))
   if(pub_dmft_mpi_size > 1)then
      filename = TRIM(ADJUSTL(filename)) // "_rank" // &
           TRIM(ADJUSTL(utils_int_to_str(pub_dmft_mpi_rank)))
   endif

   if(pub_on_root)then
      funit = utils_unit()
      if(ien == ffirst)then
         open(unit = funit, file = trim(adjustl(filename)), form = &
              'unformatted', status = 'replace')
      else
         open(unit = funit, file = trim(adjustl(filename)), form = &
              'unformatted', position = 'append', status = 'old')
      endif
      write(funit)  ien, pub_dmft_points, pub_dmft_temp, size_mat0, size_mat0
      write(funit) fermi_e(is) + pub_dmft_chem_shift, energy, &
           tmp_mat(1:size_mat0, 1:size_mat0)
      close(funit)
   endif

   deallocate(tmp_mat, stat = ierr)
   call utils_dealloc_check('dmft_embedding_write_sigma', 'tmp_mat', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_embedding_write_sigma'

 end subroutine

 subroutine dmft_embedding_read_potential(ham, ham_embed, is, ien, fermi_e, &
      ffirst, llast, reopen_sigma_file_each_step, force_reopen_file, &
      check_embed, check_embed_h, ngwf_basis)

   !---------------------------------------------------------------------------!
   ! Read the embedding potential                                              !
   !---------------------------------------------------------------------------!

   use sparse, only: sparse_convert, sparse_get_element, sparse_put_element, &
        sparse_scale
   use utils,             only: utils_abort, utils_assert, utils_unit, &
                                utils_int_to_str, utils_alloc_check, &
                                utils_dealloc_check
   use comms,             only: pub_my_proc_id, pub_on_root
   use rundat, only: pub_debug, pub_debug_on_root, &
        pub_dmft_points, pub_dmft_temp, pub_num_spins
   use function_basis,    only: FUNC_BASIS
   use parallel_strategy, only: par => pub_par

   implicit none

   ! Arguments
   type(SPAM3),      intent(inout) :: ham(pub_num_spins)
   type(SPAM3),      intent(inout) :: ham_embed(pub_num_spins)
   integer,          intent(in   ) :: is
   integer,          intent(in   ) :: ien
   real(kind = DP),  intent(in   ) :: fermi_e(pub_num_spins)
   integer,          intent(in   ) :: ffirst
   integer,          intent(in   ) :: llast
   logical,          intent(in   ) :: reopen_sigma_file_each_step
   logical,          intent(in   ) :: force_reopen_file
   logical,          intent(in   ) :: check_embed
   logical,          intent(  out) :: check_embed_h
   type(FUNC_BASIS), intent(in   ) :: ngwf_basis

   ! Local variables
   logical                         :: checkl
   logical                         :: checkr
   logical                         :: buildl
   logical                         :: buildr
   integer                         :: i, j, s, t
   complex(kind = DP)              :: sig_element, tmp, ww_
   real(kind = DP)                 :: mmuL, mmuR
   integer                         :: kk, kkk, ll, lll, jj, iii, jjj, ii1, jj1
   integer                         :: l1, l2
   complex(kind = DP)              :: www_
   integer                         :: fopenE, fopenER, fopenEL
   complex(kind = DP), allocatable :: ttt0L(:,:), ttt0R(:,:)
#ifdef savememory
   complex(kind = DP), allocatable :: SSL_(:,:), SSR_(:,:)
#else
   complex(kind = DP), allocatable :: SSL_(:,:,:), SSR_(:,:,:)
#endif
   real(kind = DP), allocatable    :: S0L_(:,:), S1L_(:,:), H0L_(:,:), H1L_(:,:)
   real(kind = DP), allocatable    :: ttL_(:,:), tttL_(:,:), ttR_(:,:), &
                                      tttR_(:,:)
   ! #ifdef GPU_SPEEDUP
   !    complex(8), allocatable         :: ttRc_(:,:), tttRc_(:,:)
   ! #endif
   complex(kind = DP), allocatable :: ttLc_(:,:), tttLc_(:,:)
   complex(kind = DP), allocatable :: GGL_(:,:), GGR_(:,:), GGLback_(:,:), &
                                      GGRback_(:,:)
   integer                         :: size_indexL, size_indexR, size_mat0L, &
                                      size_mat0R
   integer, allocatable            :: maskTL(:,:,:), maskTLatom(:,:,:), &
                                      mask_vectorL(:), maskTR(:,:,:), &
                                      maskTRatom(:,:,:), mask_vectorR(:)
   real(kind = DP), allocatable    :: mmuRar(:), mmuLar(:)
   complex(kind = DP), allocatable :: wwRar(:), wwLar(:)
   real(kind = DP), allocatable    :: S0R_(:,:), S1R_(:,:), H0R_(:,:), H1R_(:,:)
   character(30)                   :: filename
   integer                         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_embedding_read_potential'

   buildr = .false.
   buildl = .false.

   if(pub_dmft_temp > 0.0_DP)then
      filename = 'sigma_embedding'
   else
      filename = 'sigma_embedding_real'
   endif

   inquire(file = trim(adjustl(filename)) // 'L' // &
        trim(adjustl(utils_int_to_str(is))), exist = checkl)
   inquire(file = trim(adjustl(filename)) // 'R' // &
        trim(adjustl(utils_int_to_str(is))), exist = checkr)
   check_embed_h = checkr .or. checkl
   if(.not.check_embed_h) then
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_embedding_read_potential as sigma_embedding file not found'
      return
   end if

   if(pub_on_root .and. pub_dmft_mpi_rank == 1) write(stdout, *) 'READING &
        &EMBEDDING POTENTIALS, L/R PRESENT ? ', checkl, checkr
   if( (checkr .and. .not.checkl) .or. (checkl .and. .not.checkr) )then
      write(stdout, *) 'WARNING : could find file sigma_embedding_(real)L but &
           &NOT sigma_embedding_(real)R'
   endif

   call utils_assert( .not. check_embed, "Error in &
        &dmft_embedding_read_potential: both L/R lead files and &
        &sigma_embedding files found; this should not happen")

   call sparse_scale(ham_embed(is), 0.0_DP)

   if (checkl) then

      if(ien == ffirst .or. reopen_sigma_file_each_step .or. &
           force_reopen_file) then
         if(pub_dmft_mpi_rank == 1 .and. pub_on_root) write(stdout, *) &
              'opening file : ', trim(adjustl(filename)) // 'L' // &
              trim(adjustl(utils_int_to_str(is)))
         fopenEL = utils_unit()
         open(unit = fopenEL, file = trim(adjustl(filename)) // 'L' // &
              trim(adjustl(utils_int_to_str(is))), form = 'unformatted')
         read(fopenEL) size_mat0L, size_indexL
         if(pub_dmft_mpi_rank == 1 .and. pub_on_root) write(stdout, *) 'SIZE L &
              &: ', size_mat0L
         if(.not.allocated(ttt0L))then
            write(stdout, *) 'ALLOCATING EMBEDDING ARRAYS LEFT'
            buildl = .true.
            allocate(H0L_(size_mat0L, size_mat0L), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'H0L_', ierr)
            allocate(H1L_(size_mat0L, size_mat0L), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'H1L_', ierr)
            allocate(S0L_(size_mat0L, size_mat0L), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'S0L_', ierr)
            allocate(S1L_(size_mat0L, size_mat0L), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'S1L_', ierr)
#ifndef savememory
            allocate(SSL_(size_mat0L, size_mat0L, pub_dmft_points), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'SSL_', ierr)
#else
            allocate(SSL_(size_mat0L, size_mat0L), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'SSL_', ierr)
#endif
            allocate(ttL_(size_mat0L, size_mat0L), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'ttL_', ierr)
            allocate(tttL_(size_mat0L, size_mat0L), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'tttL_', ierr)
            ! #ifdef GPU_SPEEDUP
            !             allocate( ttLc_(size_mat0L, size_mat0L))
            !             allocate(tttLc_(size_mat0L, size_mat0L))
            ! #endif
            allocate(GGL_(size_mat0L, size_mat0L), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'GGL_', ierr)
            allocate(GGLback_(size_mat0L, size_mat0L), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'GGLback_', ierr)
            allocate(ttt0L(size_mat0L, size_mat0L), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'ttt0L', ierr)
            allocate(maskTLatom(size_mat0L, size_mat0L, 2), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'maskTLatom', ierr)
            allocate(maskTL(size_mat0L, size_mat0L, 2), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'maskTL', ierr)
            allocate(mask_vectorL(size_indexL), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'mask_vectorL', ierr)
            allocate(mmuLar(pub_dmft_points), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'mmuLar', ierr)
            allocate(wwLar(pub_dmft_points), stat = ierr)
            call utils_alloc_check('dmft_embedding_read_potential', 'wwLar', ierr)
            read(fopenEL) mask_vectorL, H0L_, H1L_, S0L_, S1L_, SSL_, mmuLar, &
                 wwLar
         endif
#ifdef savememory
         read(fopenEL) mask_vectorL, H0L_, H1L_, S0L_, S1L_
#endif
      endif

      ! if(pub_debug_on_root) then
      !    write(stdout, *) 'size      : ', size(mask_vectorL)
      !    write(stdout, *) 'size H0H1 : ', shape(H0L_), shape(H1L_), &
      !         shape(S0L_), shape(S1L_)
      !    write(stdout, *) 'ien       : ', ien
      !    write(stdout, *) 'LOGIC : ', reopen_sigma_file_each_step, &
      !         force_reopen_file
      ! endif

      mmuL = mmuLar(ien)
      ww_ = wwLar(ien)

#ifdef savememory
      if(.not.(reopen_sigma_file_each_step .or. force_reopen_file))then
         read(fopenEL) mmuL, SSL_, ww_
      else
         do kk = 1, ien
            read(fopenEL) mmuL, SSL_, ww_
         enddo
      endif
#endif

      if(ien == llast .or. reopen_sigma_file_each_step .or. &
           force_reopen_file)then
         close(fopenEL)
      endif

      if(pub_debug_on_root) write(stdout, *) 'DEBUG: DYSON EQUATION'
      www_ = (ww_-mmuL) + (fermi_e(is)-mmuL) + fermi_e(is)
      if(pub_debug_on_root) write(stdout, *) 'EMBEDDING SHIFT MU L : ', &
           fermi_e(is)-mmuL
      ttt0L  = 0.0_DP
      ttL_   =           H1L_  - www_*          S1L_
      tttL_  = transpose(H1L_) - www_*transpose(S1L_)
      ! #ifdef GPU_SPEEDUP
      !       ttLc_ = ttL_
      !       tttLc_ = tttL_
      !
      ! #endif
      do kk = 1, pub_dmft_embed_iter
         if(pub_dmft_mpi_rank == 1 .and. pub_on_root .and. (pub_output_detail &
              >= VERBOSE)) write(stdout, *) 'ITER EMBED', kk, real(ttt0L(1, &
              1)), aimag(ttt0L(1, 1))

#ifndef savememory
         GGL_ = www_*S0L_ - H0L_ - ttt0L - SSL_(:, :, ien)
#else
         GGL_ = www_*S0L_ - H0L_ - ttt0L - SSL_
#endif

         if (pub_debug) then
            GGL_ = (GGL_ + transpose(GGL_))/2.0_DP
            ttt0L = GGL_
         end if

         i = size(GGL_, 1)
         call dmft_sym_invert(i, GGL_)

         if (pub_debug) then
            if(pub_on_root) write(stdout, *) 'TESTING INVERSION : ', &
                 maxval(abs(matmul(GGL_, ttt0L))), minval(abs(matmul(GGL_, &
                 ttt0L)))
         end if

         ! #ifdef GPU_SPEEDUP
         ! i = size(ttL_, 1)
         ! if(i /= size(GGL_, 1))then
         !    utils_abort("Error in dmft_embedding_read_potential: dimensions do &
         !         &not match")
         ! endif
         ! if(.not.(kk == 1 .or. abs(pub_dmft_embed_mix-1.0_DP) < 1.0E-4_DP)) GGLback_ &
         !      = ttt0L
         ! #ifdef GPU_SPEEDUP_SINGPREC
         ! call matmulcuda_c_singleprec_bypass(ttLc_, GGL_, ttt0L, i, i, i)
         ! GGL_ = ttt0L
         ! call matmulcuda_c_singleprec_bypass(GGL_, tttLc_, ttt0L, i, i, i)
         ! #else
         ! call matmulcuda_c(ttLc_, GGL_, ttt0L, i, i, i)
         ! GGL_ = ttt0L
         ! call matmulcuda_c(GGL_, tttLc_, ttt0L, i, i, i)
         ! #endif
         ! if(.not.(kk == 1 .or. abs(pub_dmft_embed_mix-1.0_DP) < 1.0E-4_DP)) then
         !    ttt0L = pub_dmft_embed_mix*ttt0L + &
         !         (1.0_DP-pub_dmft_embed_mix)*GGLback_
         ! end if
         ! #else
         if(kk == 1 .or. abs(pub_dmft_embed_mix-1.0_DP) < 1.0E-4_DP)then
            ttt0L = matmul ( matmul( ttL_, GGL_ ), tttL_ )
         else
            ttt0L = pub_dmft_embed_mix*matmul ( matmul( ttL_, GGL_ ), tttL_ ) &
                 + (1.0_DP-pub_dmft_embed_mix)*ttt0L
         endif
         ! #endif

         if(pub_on_root .and. (ien == pub_dmft_points-1 .or. ien == &
              pub_dmft_points))then
            write(stdout, *) 'w = oo, frequ    = ', www_
            write(stdout, *) 'MAX RE SIGMA_L = ', maxval(abs(real(ttt0L)))
            write(stdout, *) 'MAX IM SIGMA_L = ', maxval(abs(aimag(ttt0L)))
         endif

         if (pub_debug_on_root) then
            write(stdout, *) 'DYSON ITERATION L (maxval) : ', kk, &
                 maxval(abs(ttt0L))
            write(stdout, *) 'MAXVAL SELF                : ', maxval(abs(SSL_))
         end if
      enddo

      if (pub_debug) then
         if(pub_debug_on_root) write(stdout, *) 'DEBUG: Building mask L'
         i = 0
         do j = 1, size(mask_vectorL)
            i = i + ngwf_basis%num_on_atom( par%distr_atom(mask_vectorL(j))  )
         enddo
         if(pub_debug_on_root) write(stdout, *) 'TOTAL SIZE OF ARRAYS : ', i
         if( i /= size(maskTL, 1)) then
            write(stdout, *) 'ERROR EMBEDDING'
            write(stdout, *) 'cumulated dimensions of atoms obtained from &
                 &device indices : ', i
            write(stdout, *) 'shape of mask array obtained from file : ', &
                 shape(maskTL)
            write(stdout, *) 'mask_vector : ', mask_vectorL
            call utils_abort('Error in dmft_embedding_read_potential: &
                 &dimensions of atoms does not match mask array')
         end if
      end if

      if(buildl)then
         maskTL = 0
         kk = size_indexL
         do i = 1, kk
            do j = 1, kk
               ll = 0
               do jj = 1, par%distr_atom(mask_vectorL(i))-1
                  ll = ll + ngwf_basis%num_on_atom( jj )
               enddo
               lll = 0
               do jj = 1, par%distr_atom(mask_vectorL(j))-1
                  lll = lll + ngwf_basis%num_on_atom( jj )
               enddo
               iii = 0
               do jj = 1, i-1
                  iii = iii + ngwf_basis%num_on_atom( &
                       par%distr_atom(mask_vectorL(jj)) )
               enddo
               jjj = 0
               do jj = 1, j-1
                  jjj = jjj + ngwf_basis%num_on_atom( &
                       par%distr_atom(mask_vectorL(jj)) )
               enddo
               do l1 = 1, ngwf_basis%num_on_atom( par%distr_atom(mask_vectorL( &
                    i )) )
                  do l2 = 1, ngwf_basis%num_on_atom( &
                       par%distr_atom(mask_vectorL( j )) )
                     maskTL( iii + l1, jjj + l2, 1 ) = ll  + l1
                     maskTL( iii + l1, jjj + l2, 2 ) = lll + l2
                     maskTLatom( iii + l1, jjj + l2, 1 ) = &
                          par%distr_atom(mask_vectorL( i ))
                     maskTLatom( iii + l1, jjj + l2, 2 ) = &
                          par%distr_atom(mask_vectorL( j ))
                  enddo
               enddo
            enddo
         enddo
      endif

      do i = 1, size_mat0L
         do j = 1, size_mat0L
            s = maskTL(i, j, 1)
            t = maskTL(i, j, 2)
            ii1 = maskTLatom(i, j, 1)
            jj1 = maskTLatom(i, j, 2)
            if(par%proc_of_atom(par%orig_atom(ii1)) == pub_my_proc_id .and. &
                 .not. par%proc_of_atom(par%orig_atom(jj1)) == pub_my_proc_id) &
                 then
               call utils_abort("Error in dmft_embedding_read_potential: one &
                    &atom on my proc, but not the other one")
            endif
            if(.not.par%proc_of_atom(par%orig_atom(ii1)) == pub_my_proc_id &
                 .and. par%proc_of_atom(par%orig_atom(jj1)) == pub_my_proc_id) &
                 then
               call utils_abort("Error in dmft_embedding_read_potential: one &
                    &atom one my proc, but not the other one")
            endif
            if(par%proc_of_atom(par%orig_atom(ii1)) == pub_my_proc_id .and. &
                 par%proc_of_atom(par%orig_atom(jj1)) == pub_my_proc_id) then
               sig_element = ttt0L(i, j)
               call sparse_put_element(sig_element, ham_embed(is), s, t)
            endif
         enddo
      enddo
   end if
   if(.not.checkr) return

   if(ien == ffirst .or. reopen_sigma_file_each_step .or. force_reopen_file) &
        then
      fopenER = utils_unit()
      open(unit = fopenER, file = trim(adjustl(filename)) // 'R' // &
           trim(adjustl(utils_int_to_str(is))), form = 'unformatted')
      read(fopenER) size_mat0R, size_indexR
      if(pub_dmft_mpi_rank == 1 .and. pub_on_root) write(stdout, *) 'SIZE R / &
           &size index : ', size_mat0R, size_indexR
      if(.not.allocated(ttt0R))then
         buildr = .true.
         write(stdout, *) 'ALLOCATING EMBEDDING ARRAYS RIGHT'
         allocate(H0R_(size_mat0R, size_mat0R), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'H0R_', ierr)
         allocate(H1R_(size_mat0R, size_mat0R), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'H1R_', ierr)
         allocate(S0R_(size_mat0R, size_mat0R), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'S0R_', ierr)
         allocate(S1R_(size_mat0R, size_mat0R), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'S1R_', ierr)
#ifndef savememory
         allocate(SSR_(size_mat0R, size_mat0R, pub_dmft_points), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'SSR_', ierr)
#else
         allocate(SSR_(size_mat0R, size_mat0R), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'SSR_', ierr)
#endif
         allocate(ttR_(size_mat0R, size_mat0R), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'ttR_', ierr)
         allocate(tttR_(size_mat0R, size_mat0R), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'tttR_', ierr)
         ! #ifdef GPU_SPEEDUP
         ! allocate( ttRc_(size_mat0R, size_mat0R))
         ! allocate(tttRc_(size_mat0R, size_mat0R))
         ! #endif
         allocate(GGR_(size_mat0R, size_mat0R), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'GGR_', ierr)
         allocate(GGRback_(size_mat0R, size_mat0R), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'GGRback_', &
              ierr)
         allocate(ttt0R(size_mat0R, size_mat0R), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'ttt0R', ierr)
         allocate(maskTR(size_mat0R, size_mat0R, 2), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'maskTR', ierr)
         allocate(maskTRatom(size_mat0R, size_mat0R, 2), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'maskTRatom', &
              ierr)
         allocate(mask_vectorR(size_indexR), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', &
              'mask_vectorR', ierr)
         allocate(mmuRar(pub_dmft_points), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'mmuRar', ierr)
         allocate(wwRar(pub_dmft_points), stat = ierr)
         call utils_alloc_check('dmft_embedding_read_potential', 'wwRar', ierr)
         read(fopenER) mask_vectorR, H0R_, H1R_, S0R_, S1R_, SSR_, mmuRar, wwRar
      endif
#ifdef savememory
      read(fopenER) mask_vectorR, H0R_, H1R_, S0R_, S1R_
#endif
   endif

   ! if(pub_debug_on_root) then
   !    write(stdout, *) 'size      : ', size(mask_vectorR)
   !    write(stdout, *) 'size H0H1 : ', shape(H0R_), shape(H1R_), shape(S0R_), &
   !         shape(S1R_)
   !    write(stdout, *) 'ien       : ', ien
   !    write(stdout, *) 'LOGIC : ', reopen_sigma_file_each_step, &
   !         force_reopen_file
   ! endif

   mmuR = mmuRar(ien)
   ww_ = wwRar(ien)

#ifdef savememory
   if(.not.(reopen_sigma_file_each_step .or. force_reopen_file))then
      read(fopenER) mmuR, SSR_, ww_
   else
      do kk = 1, ien
         read(fopenER) mmuR, SSR_, ww_
      enddo
   endif
#endif

   if(ien == llast .or. reopen_sigma_file_each_step .or. force_reopen_file)then
      close(fopenER)
   endif

   if (pub_debug) then
      write(stdout, *) 'ASSYMETRY H0 : ', maxval(abs(H0R_-transpose(H0R_)))
      write(stdout, *) 'ASSYMETRY S0 : ', maxval(abs(S0R_-transpose(S0R_)))
      write(stdout, *) 'ASSYMETRY SS : ', maxval(abs(SSR_(:, :, &
           ien)-transpose(SSR_(:, :, ien))))
      write(stdout, *) 'SYMETRY   H0 : ', maxval(abs(H0R_ + transpose(H0R_)))
      write(stdout, *) 'SYMETRY   S0 : ', maxval(abs(S0R_ + transpose(S0R_)))
      write(stdout, *) 'SYMETRY SS : ', maxval(abs(SSR_(:, :, ien) + &
           transpose(SSR_(:, :, ien))))
   end if

   if(pub_debug_on_root) write(stdout, *) 'DEBUG: DYSON EQUATION'
   www_ = (ww_-mmuR) + (fermi_e(is)-mmuR) + fermi_e(is)
   if(pub_debug_on_root)  write(stdout, *) 'EMBEDDING SHIFT MU R : ', fermi_e(is)-mmuR
   ttt0R  =  0.0_DP
   ttR_   =           H1R_  - www_*          S1R_
   tttR_  = transpose(H1R_) - www_*transpose(S1R_)
   ! #ifdef GPU_SPEEDUP
   ! ttRc_ = ttR_
   ! tttRc_ = tttR_
   ! #endif

   do kk = 1, pub_dmft_embed_iter
      if(pub_dmft_mpi_rank == 1 .and. pub_on_root .and. (pub_output_detail >= &
           VERBOSE)) write(stdout, *) 'ITER EMBED', kk, real(ttt0R(1, 1)), &
           aimag(ttt0R(1, 1))
#ifdef savememory
      GGR_ = www_*S0R_ - H0R_ - ttt0R - SSR_
#else
      GGR_ = www_*S0R_ - H0R_ - ttt0R - SSR_(:, :, ien)
#endif
      if (pub_debug) GGR_ = (GGR_ + transpose(GGR_))/2.0_DP

      i = size(GGR_, 1)
      call dmft_sym_invert(i, GGR_)

      ! #ifdef GPU_SPEEDUP
      ! i = size(ttR_, 1)
      ! call utils_assert(i == size(GGR_, 1), "Error in &
      !      &dmft_embedding_read_potential: gpu speedup dimensions do not match &
      !      &R")
      ! if(.not.(kk == 1 .or. abs(pub_dmft_embed_mix-1.0_DP) < 1.0E-4_DP)) GGRback_ = &
      !      ttt0R
      ! #ifdef GPU_SPEEDUP_SINGPREC
      ! call matmulcuda_c_singleprec_bypass(ttRc_, GGR_, ttt0R, i, i, i)
      ! GGR_ = ttt0R
      ! call matmulcuda_c_singleprec_bypass(GGR_, tttRc_, ttt0R, i, i, i)
      ! #else
      ! call matmulcuda_c(ttRc_, GGR_, ttt0R, i, i, i)
      ! GGR_ = ttt0R
      ! call matmulcuda_c(GGR_, tttRc_, ttt0R, i, i, i)
      ! #endif
      ! if(.not.(kk == 1 .or. abs(pub_dmft_embed_mix-1.0_DP) < 1.0E-4_DP)) ttt0R = &
      !      pub_dmft_embed_mix*ttt0R + (1.0_DP-pub_dmft_embed_mix)*GGRback_
      ! #else
      if(kk == 1 .or. abs(pub_dmft_embed_mix-1.0_DP) < 1.0E-4_DP)then
         ttt0R = matmul ( matmul( ttR_, GGR_ ), tttR_ )
      else
         ttt0R = pub_dmft_embed_mix*matmul ( matmul( ttR_, GGR_ ), tttR_ ) + &
              (1.0_DP-pub_dmft_embed_mix)*ttt0R
      endif
      ! #endif

      if(pub_on_root .and. (ien == pub_dmft_points-1 .or. ien == &
           pub_dmft_points))then
         write(stdout, *) 'w = oo, frequ    = ', www_
         write(stdout, *) 'MAX RE SIGMA_R = ', maxval(abs(real(ttt0R)))
         write(stdout, *) 'MAX IM SIGMA_R = ', maxval(abs(aimag(ttt0R)))
      endif
      if(pub_debug_on_root) then
         write(stdout, *) 'DYSON ITERATION R (max) : ', kk, maxval(abs(ttt0R))
         write(stdout, *) '  G (max)               : ', maxval(abs(GGR_))
      endif
   enddo

   if (pub_debug) then
      if(pub_debug_on_root) write(stdout, *) 'DEBUG: building mask R'
      i = 0
      do j = 1, size(mask_vectorR)
         i = i + ngwf_basis%num_on_atom( par%distr_atom(mask_vectorR(j))  )
      enddo
      if(pub_debug_on_root)  write(stdout, *) 'TOTAL SIZE OF ARRAYS : ', i
      if( i /= size(maskTR, 1)) then
         write(stdout, *) 'ERROR EMBEDDING'
         write(stdout, *) 'cumulated dimensions of atoms obtained from device &
              &indices : ', i
         write(stdout, *) 'shape of mask array obtained from file : ', &
              shape(maskTR)
         write(stdout, *) 'mask_vector : ', mask_vectorR
         call utils_abort('Error in DMFT: dimensions of atoms does not match &
              &mask array')
      endif
   endif

   !erase ham by boundary
   if(buildr)then
      kk  = size_indexR
      kkk = size_indexL
      do i = 1, kk
         do j = 1, kkk
            call utils_assert(mask_vectorR(i) /= mask_vectorL(j), "Error in &
                 &dmft_embedding_read_potential: same site cannot be in both &
                 &right and left leads")
            ll = 0
            do jj = 1, par%distr_atom(mask_vectorR(i))-1
               ll = ll + ngwf_basis%num_on_atom( jj )
            enddo
            lll = 0
            do jj = 1, par%distr_atom(mask_vectorL(j))-1
               lll = lll + ngwf_basis%num_on_atom( jj )
            enddo
            do l1 = 1, ngwf_basis%num_on_atom( par%distr_atom(mask_vectorR( i &
                 )) )
               do l2 = 1, ngwf_basis%num_on_atom( par%distr_atom(mask_vectorL( &
                    j )) )
                  s = ll + l1
                  t = lll + l2
                  call sparse_put_element(0.0_DP, ham(is), s, t)
                  call sparse_put_element(0.0_DP, ham(is), t, s)
               enddo
            enddo
         enddo
      enddo
   endif
   !end erase ham boundary

   if(buildr)then
      if(pub_debug_on_root)  write(stdout, *) 'DEBUG: building mask R'
      maskTR = 0
      kk = size_indexR
      do i = 1, kk
         do j = 1, kk
            ll = 0
            do jj = 1, par%distr_atom(mask_vectorR(i))-1
               ll = ll + ngwf_basis%num_on_atom( jj )
            enddo
            lll = 0
            do jj = 1, par%distr_atom(mask_vectorR(j))-1
               lll = lll + ngwf_basis%num_on_atom( jj )
            enddo
            iii = 0
            do jj = 1, i-1
               iii = iii + ngwf_basis%num_on_atom( &
                    par%distr_atom(mask_vectorR(jj)) )
            enddo
            jjj = 0
            do jj = 1, j-1
               jjj = jjj + ngwf_basis%num_on_atom( &
                    par%distr_atom(mask_vectorR(jj)) )
            enddo
            do l1 = 1, ngwf_basis%num_on_atom( par%distr_atom(mask_vectorR( i &
                 )) )
               do l2 = 1, ngwf_basis%num_on_atom( par%distr_atom(mask_vectorR( &
                    j )) )
                  maskTR( iii + l1, jjj + l2, 1 ) = ll  + l1
                  maskTR( iii + l1, jjj + l2, 2 ) = lll + l2
                  maskTRatom( iii + l1, jjj + l2, 1 ) = &
                       par%distr_atom(mask_vectorR( i ))
                  maskTRatom( iii + l1, jjj + l2, 2 ) = &
                       par%distr_atom(mask_vectorR( j ))
               enddo
            enddo
         enddo
      enddo
   endif

   do i = 1, size_mat0R
      do j = 1, size_mat0R
         s = maskTR(i, j, 1)
         t = maskTR(i, j, 2)
         ii1 = maskTRatom(i, j, 1)
         jj1 = maskTRatom(i, j, 2)
         if(par%proc_of_atom(par%orig_atom(ii1)) == pub_my_proc_id .and. .not. &
              par%proc_of_atom(par%orig_atom(jj1)) == pub_my_proc_id) then
            call utils_abort("Error in dmft_embedding_read_potential: one, but &
                 &not both atoms on proc")
         endif
         if(.not.par%proc_of_atom(par%orig_atom(ii1)) == pub_my_proc_id .and. &
              par%proc_of_atom(par%orig_atom(jj1)) == pub_my_proc_id) then
            call utils_abort("Error in dmft_embedding_read_potential: one, but &
                 &not both atoms on proc")
         endif
         if(par%proc_of_atom(par%orig_atom(ii1)) == pub_my_proc_id .and. &
              par%proc_of_atom(par%orig_atom(jj1)) == pub_my_proc_id) then
            sig_element = ttt0R(i, j)
            tmp = 0.0_DP
            call sparse_get_element(        tmp, ham_embed(is), s, t)
            if (pub_debug) then
               if(abs(tmp) > 1.0E-5_DP)then
                  write(stdout, *) 'SITE 1 : ', par%orig_atom(ii1)
                  write(stdout, *) 'SITE 2 : ', par%orig_atom(jj1)
                  write(stdout, *) 'connected to the L and R leads...'
                  write(stdout, *) 'amplitude : ', tmp
                  write(stdout, *) 'should not happen....'
                  call utils_abort('Error in dmft_embedding_read_potential: &
                       &amplitude exceeds 1e-5, should not happen.')
               endif
            end if
            sig_element = sig_element + tmp
            call sparse_put_element(sig_element, ham_embed(is), s, t)
         endif
      enddo
   enddo

   ! Cleaning
   if(allocated(mmuRar)) then
      deallocate(mmuRar, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'mmuRar', ierr)
      deallocate(wwRar, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'wwRar', ierr)
   end if
   if(allocated(mmuLar)) then
      deallocate(mmuLar, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'mmuLar', ierr)
      deallocate(wwLar, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'wwLar', ierr)
   end if
   if(allocated(ttt0R)) then
      deallocate(ttt0R, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'ttt0R', ierr)
      deallocate(maskTR, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'maskTR', ierr)
      deallocate(maskTRatom, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'maskTRatom', &
           ierr)
      deallocate(mask_vectorR, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', &
           'mask_vectorR', ierr)
   end if
   if(allocated(H0R_)) then
      deallocate(H0R_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'H0R_', ierr)
      deallocate(H1R_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'H1R_', ierr)
      deallocate(S0R_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'S0R_', ierr)
      deallocate(S1R_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'S1R_', ierr)
      deallocate(SSR_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'SSR_', ierr)
   end if
   if(allocated(ttR_)) then
      deallocate(ttR_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'ttR_', ierr)
      deallocate(tttR_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'tttR_', ierr)
      deallocate(GGR_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'GGR_', ierr)
      deallocate(GGRback_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'GGRback_', &
           ierr)
   end if
   if(allocated(ttt0L)) then
      deallocate(ttt0L, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'ttt0L', ierr)
      deallocate(maskTL, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'maskTL', ierr)
      deallocate(maskTLatom, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'maskTLatom', &
           ierr)
      deallocate(mask_vectorL, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', &
           'mask_vectorL', ierr)
   end if
   if(allocated(H0L_)) then
      deallocate(H0L_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'H0L_', ierr)
      deallocate(H1L_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'H1L_', ierr)
      deallocate(S0L_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'S0L_', ierr)
      deallocate(S1L_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'S1L_', ierr)
      deallocate(SSL_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'SSL_', ierr)
   end if
   if(allocated(ttL_)) then
      deallocate(ttL_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'ttL_', ierr)
      deallocate(tttL_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'tttL_', ierr)
      deallocate(GGL_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'GGL_', ierr)
      deallocate(GGLback_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'GGLback_', &
           ierr)
   end if
   ! #ifdef GPU_SPEEDUP
   ! if(allocated(ttLc_))
   !    deallocate(ttLc_, stat = ierr)
   !    call utils_dealloc_check('dmft_embedding_read_potential', 'ttLc_', ierr)
   !    deallocate(tttLc_, stat = ierr)
   !    call utils_dealloc_check('dmft_embedding_read_potential', 'tttLc_', ierr)
   !    deallocate(ttRc_, stat = ierr)
   !    call utils_dealloc_check('dmft_embedding_read_potential', 'ttRc_', ierr)
   !    deallocate(tttRc_, stat = ierr)
   !    call utils_dealloc_check('dmft_embedding_read_potential', 'tttRc_', ierr)
   ! end if
   ! #else
   if(allocated(ttLc_)) then
      deallocate(ttLc_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'ttLc_', ierr)
      deallocate(tttLc_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'tttLc_', ierr)
      deallocate(ttLc_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'ttLc_', ierr)
      deallocate(tttLc_, stat = ierr)
      call utils_dealloc_check('dmft_embedding_read_potential', 'tttLc_', ierr)
   end if
   ! #endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_embedding_read_potential'

 end subroutine

 subroutine dmft_fix_fermi_and_total_energies(total_energy, mdl, ngwf_basis, &
      nkpoints, E_sum_eigenvalues_DFT, sc_data, fermi_e, eigen_en, n_occ, &
      overlap, ham, full_prec_dmft_energy, use_gpu_eigenvectors, &
      dmft_energy_cor)

   !---------------------------------------------------------------------------!
   ! Starting from the traces of the Green's function and self energies, this  !
   ! subroutine fixes (via the Migdal formula) the Fermi energy and updates    !
   ! energies                                                                  !
   !---------------------------------------------------------------------------!

   ! #ifdef GPU_SPEEDUP
   ! use fortran_cuda,   only: eigenvector_gen_cuda_r
   ! #endif
   use linalg,         only: linalg_dsygv_lt
   use dense,          only: dense_convert
   use sparse,         only: sparse_convert, sparse_copy
   use utils,          only: utils_abort, utils_assert, utils_unit, &
        utils_banner, utils_alloc_check, utils_dealloc_check
   use comms,          only: pub_on_root
   use rundat, only: pub_debug, pub_debug_on_root, pub_dmft_read, &
        pub_num_spins, pub_rootname
   use function_basis, only: FUNC_BASIS
   use model_type,     only: MODEL

   implicit none

   ! Arguments
   real(kind = DP),              intent(inout) :: total_energy
   type(MODEL),                  intent(in   ) :: mdl
   type(FUNC_BASIS),             intent(in   ) :: ngwf_basis
   integer,                      intent(in   ) :: nkpoints
   real(kind = DP),              intent(inout) :: E_sum_eigenvalues_DFT(2)
   type(DMFT_SC_MATRICES),       intent(inout) :: sc_data
   real(kind = DP),              intent(inout) :: fermi_e(2)
   real(kind = DP), allocatable, intent(inout) :: eigen_en(:,:)
   integer, allocatable,         intent(inout) :: n_occ(:)
   type(SPAM3),                  intent(inout) :: overlap
   type(SPAM3),                  intent(inout) :: ham(pub_num_spins)
   logical,                      intent(in   ) :: full_prec_dmft_energy
   logical,                      intent(in   ) :: use_gpu_eigenvectors
   real(kind = DP), optional,    intent(in   ) :: dmft_energy_cor

   ! Local variables
   real(kind = DP)              :: lhxc_energy
   real(kind = DP)              :: hubbard_energy
   logical                      :: check_chem, check
   integer                      :: ii, free_unit, is
   real(kind = DP), allocatable :: overlap_square(:,:)
   integer                      :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_fix_fermi_and_total_energies'

   free_unit = utils_unit()

   inquire(file = trim(pub_rootname) // '.tot_energy', exist = check)
   if(check)then
      open(unit = free_unit, file = trim(pub_rootname) // '.tot_energy')
      read(free_unit, *) total_energy, lhxc_energy, hubbard_energy
      if (pub_debug_on_root) then
         write(stdout, "(a)") utils_banner("-", ' Energy summary ')
         write(stdout, "(a,f14.6)") ' Total                    : ', &
              total_energy
         write(stdout, "(a,f14.6)") ' Ewald                    : ', &
              mdl%ewald_energy
         write(stdout, "(a,f14.6)") ' lHxc_fine                : ', &
              lhxc_energy
         write(stdout, "(a,f14.6)") ' Hubbard                  : ', &
              hubbard_energy
      end if
      close(free_unit)
   endif

   if(present(dmft_energy_cor)) total_energy = 0.0_DP

   ii = ngwf_basis%num

   if(present(dmft_energy_cor))then
      call utils_assert(nkpoints <= 1, "Error in &
           &dmft_fix_fermi_and_total_energies: dmft_energy_cor and k-points &
           &not yet implemented")
      do is = 1, pub_num_spins
         if(.not.full_prec_dmft_energy .and. .not.(pub_dmft_use_gpu_onlyme &
              .and. use_gpu_eigenvectors))then
            call dmft_eigenvec_fullmat(ii, ham(is), overlap, &
                 sc_data%eigenvecs(is), sc_data%eigenvals(is, :))
         else
            call dense_convert(sc_data%eigenvecs(is), ham(is))
            call sparse_copy(sc_data%overlap, overlap)
            if(.not.(pub_dmft_use_gpu_onlyme .and. use_gpu_eigenvectors))then
               allocate(overlap_square(ii, ii), stat = ierr)
               call utils_alloc_check('dmft_fix_fermi_and_total_energies', &
                    'overlap_square', ierr)
               call sparse_convert(overlap_square, sc_data%overlap)
               call linalg_dsygv_lt(eigenvectors = sc_data%eigenvecs(is)%dmtx, &
                    eigenvalues = sc_data%eigenvals(is, :), overlap_square = &
                    overlap_square, num = ii)
               deallocate(overlap_square, stat = ierr)
               call utils_dealloc_check('dmft_fix_fermi_and_total_energies', &
                    'overlap_square', ierr)
            else
               ! #ifdef GPU_SPEEDUP
               ! call eigenvector_gen_cuda_r(ii, sc_data%eigenvecs(is, :, :),
               ! sc_data%overlap, sc_data%eigenvals(is, :), tmpr, .false.)
               !           sc_data%eigenvecs(is, :, :) = tmpr
               ! #else
               call utils_abort('Error in dmft_fix_fermi_and_total_energies: &
                    &error gpu 1')
               ! #endif
            endif
         endif
         E_sum_eigenvalues_DFT(is) = sum(sc_data%eigenvals(is, 1:n_occ(is)))
         if(pub_on_root)write(stdout, *) 'SUM EIGENVALUES OF H : ', &
              E_sum_eigenvalues_DFT(is)
         if(pub_on_root)write(stdout, *) 'NOCC                 : ', n_occ(is)
         eigen_en(:, is) = sc_data%eigenvals(is, :)
         fermi_e(is) = 0.5_DP*(eigen_en(n_occ(is), is) + eigen_en(n_occ(is) + &
              1, is))
      enddo
      if(pub_num_spins == 1)then
         !paramagnetic
         if(size(fermi_e)   == 2)    fermi_e(2) = fermi_e(1)
         if(size(eigen_en, 2) == 2) eigen_en(:, 2) = eigen_en(:, 1)
      endif
   else
      E_sum_eigenvalues_DFT(1) = sum(eigen_en(1:n_occ(1), 1))
      if(pub_num_spins == 2) E_sum_eigenvalues_DFT(2) = &
           sum(eigen_en(1:n_occ(2), 2))
      do is = 1, pub_num_spins
         fermi_e(is) = 0.5_DP*(eigen_en(n_occ(is), is) + eigen_en(n_occ(is) + &
              1, is))
      enddo
   endif

   ! if(pub_on_root .and. (pub_output_detail >= VERBOSE))then
   !    write(stdout, *) 'FERMI E IS IN THE Mid-gap, shift it with &
   !         &pub_dmft_chem_pot'
   !    write(stdout, *) 'EIGEN_EN SHAPE      : ', shape(eigen_en)
   !    write(stdout, *) 'MPI RANK            : ', pub_dmft_mpi_rank
   !    write(stdout, *) 'MPI SIZE            : ', pub_dmft_mpi_size
   !    write(stdout, *) 'number of states    : ', n_occ(:)
   !    do is = 1, pub_num_spins
   !       write(stdout, *) 'eigen E last state  : ', eigen_en(n_occ(is), is)
   !       write(stdout, *) 'eigen E last + 1     : ', eigen_en(n_occ(is) + 1, is)
   !    enddo
   !    write(stdout, *) 'sum occ eigenvalues : ', E_sum_eigenvalues_DFT
   ! endif

   inquire(file = 'chem.potential.nmu.iter1', exist = check_chem)
   if(check_chem .and. pub_dmft_read)then
      if (pub_on_root) write(stdout, '(a)') 'Reading chemical potential (spin &
           &up) from file'
      open(unit = free_unit, file = 'chem.potential.nmu.iter1')
      read(free_unit, *) fermi_e(1)
      close(free_unit)
      if(pub_num_spins == 2)then
         inquire(file = 'chem.potential.nmu.iter2', exist = check_chem)
         if(check_chem)then
            if (pub_on_root) write(stdout, '(a)') 'Reading chemical potential &
                 &(spin down) from file'
            open(unit = free_unit, file = 'chem.potential.nmu.iter2')
            read(free_unit, *) fermi_e(2)
            close(free_unit)
         else
            fermi_e(2) = fermi_e(1)
         endif
      endif
   endif

   if (pub_debug_on_root) then
      write(stdout, "(a,f11.6)") ' Fermi energy (spin up)   : ', fermi_e(1)
      write(stdout, "(a,f11.6)") '              (spin down) : ', fermi_e(2)
   end if

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_fix_fermi_and_total_energies'

 end subroutine

 subroutine dmft_fix_chemical_potential_loop(n_occ, chem_shift, target_N, &
      fermi_E, vectorMu, vectorNmu, num_spins, nmu_step, nmu)

   !---------------------------------------------------------------------------!
   ! Calculates a single chemical potential from the (spin-dependent) Fermi    !
   ! levels                                                                    !
   !---------------------------------------------------------------------------!

   use comms,  only: pub_on_root
   use rundat, only: pub_debug_on_root, pub_dmft_mu_diff_max

   implicit none

   ! Arguments
   integer, allocatable, intent(in   ) :: n_occ(:)
   real(kind = DP),      intent(inout) :: chem_shift
   real(kind = DP),      intent(  out) :: target_N
   real(kind = DP),      intent(inout) :: fermi_e(2)
   real(kind = DP),      intent(inout) :: vectorMu(1000)
   real(kind = DP),      intent(inout) :: vectorNmu(1000)
   integer,              intent(in   ) :: num_spins
   integer,              intent(inout) :: nmu_step
   integer,              intent(in   ) :: nmu

   ! Local variables
   real(kind = DP)            :: chem_left
   real(kind = DP)            :: chem_right
   real(kind = DP), parameter :: chem_left_shift = 0.1_DP
   real(kind = DP), parameter :: chem_right_shift = 0.1_DP
   integer                    :: is

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_fix_chemical_potential_loop'

   target_N = sum(n_occ(1:num_spins)) + pub_dmft_doping

   ! if(pub_debug_on_root) &
   !      write(stdout, *) 'target N : ', target_N

   nmu_step = nmu_step + 1
   if(nmu_step > 1)then
      if(abs( vectorNmu(nmu_step-1) - target_N ) < pub_dmft_mu_diff_max)then
         vectorNmu(nmu-1) = vectorNmu(nmu_step-1)
         vectorMu(nmu-1) = vectorMu(nmu_step-1)
         nmu_step = nmu
         chem_shift = 0.0_DP
      endif
   endif

   chem_left  = minval(fermi_e(1:num_spins)) - chem_left_shift
   chem_right = maxval(fermi_e(1:num_spins)) + chem_right_shift

   ! if(pub_debug_on_root) &
   !      write(stdout, *) 'chem pot brackets : ', chem_left, chem_right

   vectorNmu(nmu_step) = 0.0_DP

   do is = 1, num_spins
      fermi_e(is) = fermi_e(is) + chem_shift
      if(fermi_e(is) < chem_left )  fermi_e(is) = chem_left
      if(fermi_e(is) > chem_right)  fermi_e(is) = chem_right
   enddo
   vectorMu(nmu_step) = sum(fermi_e(1:num_spins))/real(num_spins, kind=DP)
   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_fix_chemical_potential_loop'

 end subroutine

 subroutine dmft_fix_chemical_potential_loop_spin(n_occ, chem_shift, target_N, &
      fermi_E, is, vectorMu, vectorNmu, num_spins, nmu_step, nmu)

   !---------------------------------------------------------------------------!
   ! Updates the chemical potential for a single spin channel                  !
   !---------------------------------------------------------------------------!

   use comms,  only: pub_on_root
   use rundat, only: pub_debug_on_root, pub_dmft_mu_diff_max

   implicit none

   ! Arguments
   integer, allocatable, intent(in   ) :: n_occ(:)
   real(kind = DP),      intent(inout) :: chem_shift
   real(kind = DP),      intent(  out) :: target_N
   real(kind = DP),      intent(inout) :: fermi_e(2)
   integer,              intent(in   ) :: is
   real(kind = DP),      intent(inout) :: vectorMu(1000)
   real(kind = DP),      intent(inout) :: vectorNmu(1000)
   integer,              intent(in   ) :: num_spins
   integer,              intent(inout) :: nmu_step
   integer,              intent(in   ) :: nmu

   ! Local variables
   real(kind = DP)            :: chem_left
   real(kind = DP)            :: chem_right
   real(kind = DP), parameter :: chem_left_shift = 0.1_DP
   real(kind = DP), parameter :: chem_right_shift = 0.1_DP

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_fix_chemical_potential_loop_spin'

   target_N = n_occ(is) + pub_dmft_doping/2.0_DP

   if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, *) &
        'target N : ', target_N

   nmu_step = nmu_step + 1
   if(nmu_step > 1)then
      if(abs( vectorNmu(nmu_step-1) - target_N ) < pub_dmft_mu_diff_max)then
         vectorNmu(nmu-1) = vectorNmu(nmu_step-1)
         vectorMu(nmu-1) = vectorMu(nmu_step-1)
         nmu_step = nmu
         chem_shift = 0.0_DP
      endif
   endif

   chem_left  = fermi_e(is) - chem_left_shift
   chem_right = fermi_e(is) + chem_right_shift

   if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, *) &
      'chem pot brackets : ', chem_left, &
        chem_right

   vectorNmu(nmu_step) = 0.0_DP

   fermi_e(is) = fermi_e(is) + chem_shift
   if(fermi_e(is) < chem_left )  fermi_e(is) = chem_left
   if(fermi_e(is) > chem_right)  fermi_e(is) = chem_right

   vectorMu(nmu_step) = fermi_e(is)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_fix_chemical_potential_loop_spin'

 end subroutine

 pure integer function dmft_hub_atom_hub_projs(hat, hub, hub_proj_basis)

   !---------------------------------------------------------------------------!
   ! Returns the number of hubbard projectors on a given hubbard atom          !
   !---------------------------------------------------------------------------!

   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par
   use function_basis,    only: FUNC_BASIS

   implicit none

   ! Arguments
   integer,             intent(in   ) :: hat
   type(HUBBARD_MODEL), intent(in   ) :: hub
   type(FUNC_BASIS),    intent(in   ) :: hub_proj_basis

   ! Local variables
   integer :: theatom

   theatom = par%distr_atom(hub%orig(hat))
   dmft_hub_atom_hub_projs = hub_proj_basis%num_on_atom(theatom)
 end function


 subroutine dmft_define_frequencies(myfrequencies, myfrequencies_numb, ffirst, &
      llast, fermi_e, energy, en_start, en_step, check_embed, restrict_window)

   !---------------------------------------------------------------------------!
   ! Define the Matsubara frequencies                                          !
   !---------------------------------------------------------------------------!

   use utils,  only: utils_abort, utils_unit, utils_int_to_str, &
                     utils_alloc_check
   use comms,  only: pub_on_root
   use rundat, only: pub_debug_on_root, pub_dmft_dos_max, &
        pub_dmft_dos_min, pub_dmft_points

   implicit none

   ! Arguments
   integer, allocatable, intent(  out) :: myfrequencies(:)
   integer,              intent(  out) :: myfrequencies_numb
   integer,              intent(  out) :: ffirst
   integer,              intent(  out) :: llast
   real(kind = DP),      intent(in   ) :: fermi_e(2)
   complex(kind = DP),   intent(  out) :: energy
   real(kind = DP),      intent(in   ) :: en_start
   real(kind = DP),      intent(in   ) :: en_step
   logical,              intent(in   ) :: check_embed
   logical,              intent(in   ) :: restrict_window

   ! Local variables
   real(kind = DP) :: profile_ratio
   integer         :: iien, is, ien, i
   integer         :: iistart, iistep
   integer         :: iostatus
   integer         :: free_unit
   logical         :: check_profile
   integer         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_define_frequencies'

   allocate(myfrequencies(pub_dmft_points), stat = ierr)
   call utils_alloc_check('dmft_define_frequencies', 'myfrequencies', ierr)

   myfrequencies = 0
   myfrequencies_numb = 0
   inquire(file = "mask_frequ" // trim(adjustl(utils_int_to_str( &
        pub_dmft_mpi_rank))), exist = check_profile)

   if(check_profile .and. .not.pub_dmft_lin_scaling .and. .not.pub_dmft_split &
        .and. .not.pub_dmft_splitk)then

      if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, *) &
           'OPENING FREQUENCIES PROFILE'
      free_unit = utils_unit()
      open(unit = free_unit, file = "mask_frequ" // &
           trim(adjustl(utils_int_to_str(pub_dmft_mpi_rank))))
      myfrequencies_numb = 0
      read(free_unit, *) profile_ratio

      if(abs(profile_ratio) < 1.0E-3_DP)then
         do
            read(free_unit, *, iostat = iostatus) iien
            if (iostatus /= 0) exit
            myfrequencies_numb = myfrequencies_numb + 1
            myfrequencies(myfrequencies_numb) = iien
         enddo

      else
         if(pub_dmft_mpi_size /= 2) then
            call utils_abort('Error in dmft_define_frequencies: if not using 2 &
                 &GPUS, the first line of mask_frequ should be "0.00"')
         endif
         if(pub_dmft_mpi_rank == 1)then
            myfrequencies_numb = NINT(profile_ratio &
                 * real(pub_dmft_points, kind=DP))
            myfrequencies(1:myfrequencies_numb) = (/( i, i = 1, &
                 myfrequencies_numb )/)
         else
            myfrequencies_numb = pub_dmft_points - &
                 NINT(profile_ratio*real(pub_dmft_points, kind=DP))
            myfrequencies(1:myfrequencies_numb) = (/( i, i = &
                 pub_dmft_points-myfrequencies_numb + 1, pub_dmft_points )/)
         endif
      endif
      close(free_unit)
   else
      is = 1
      if(.not.pub_dmft_splitk) then
         iistart = (pub_dmft_mpi_rank-1)
         iistep = pub_dmft_mpi_size
      else
         iistart = 0
         iistep = 1
      endif
      do ien = iistart + 1, pub_dmft_points, iistep
         call dmft_define_energy(energy, fermi_e, is, ien, en_start, en_step)
         if( .not. restrict_window .or. check_embed .or. ( real(energy) - &
              (fermi_e(is) + pub_dmft_chem_shift) > pub_dmft_dos_min .and. &
              real(energy) - (fermi_e(is) + pub_dmft_chem_shift) < &
              pub_dmft_dos_max ) ) then
            myfrequencies_numb = myfrequencies_numb + 1
            myfrequencies(myfrequencies_numb) = ien
         endif
      enddo
   endif

   if(pub_debug_on_root) write(stdout, '(a)') 'DEBUG: define first/last &
        &frequencies'

   ffirst = myfrequencies(1)
   llast  = myfrequencies(myfrequencies_numb)

   if(pub_debug_on_root)then
      write(stdout, *) 'FIRST FREQU   :  ', ffirst
      write(stdout, *) 'LAST  FREQU   :  ', llast
      write(stdout, *) 'MY RANK       :  ', pub_dmft_mpi_rank
      write(stdout, *) 'my frequs     :  ', myfrequencies(1:myfrequencies_numb)
      write(stdout, *) 'total number  :  ', myfrequencies_numb
   endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_define_frequencies'

 end subroutine

 subroutine dmft_define_energy(energy, fermi_e, is, ien, en_start, en_step)

   !---------------------------------------------------------------------------!
   ! Defines the energy used in the Green's function for the current Matsubara !
   ! frequency. The final few points are special points used to calculate the  !
   ! moments of the Green's function                                           !
   !---------------------------------------------------------------------------!

   use rundat, only: pub_debug_on_root, &
        pub_dmft_points, pub_dmft_smear, &
        pub_dmft_smear_T, pub_dmft_smear_eta, pub_dmft_smear_shift, &
        pub_dmft_smear_w, pub_dmft_temp

   implicit none

   ! Arguments
   real(kind = DP),    intent(in   ) :: fermi_e(2)
   integer,            intent(in   ) :: is
   integer,            intent(in   ) :: ien
   real(kind = DP),    intent(in   ) :: en_start
   real(kind = DP),    intent(in   ) :: en_step
   complex(kind = DP), intent(  out) :: energy

   ! Local variables
   real(kind = DP) :: ww
   real(kind = DP) :: realfrequ
   real(kind = DP) :: smearing

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_define_energy'

   if (pub_dmft_temp > 0.0_DP) then
      !-----------------------!
      ! Matsubara frequencies !
      !-----------------------!
      energy = cmplx( (fermi_e(is) + pub_dmft_chem_shift), en_start + &
           real(ien-1, kind = DP)*en_step, kind = DP)
      if(ien == pub_dmft_points .or. ien == pub_dmft_points-1) energy = cmplx( &
           fermi_e(is) + pub_dmft_chem_shift, pub_dmft_free_green_frequ, kind &
           = DP)
      if(ien == pub_dmft_points-2) energy = cmplx( fermi_e(is) + &
           pub_dmft_chem_shift, 0.0_DP, kind = DP)
      if(ien == pub_dmft_points-3) energy = (1.0_DP, 0.0_DP)
      if(ien == pub_dmft_points-6) energy = (1.0_DP, 0.0_DP)
      if(ien == pub_dmft_points-4) energy = (0.0_DP, 0.0_DP)
      if(ien == pub_dmft_points-5) energy = (0.0_DP, 0.0_DP)
   else
      !------------------!
      ! Real frequencies !
      !------------------!
      realfrequ = en_start + real(ien-1, kind = DP) * en_step
      ww = realfrequ - pub_dmft_smear_shift
      smearing =  pub_dmft_smear
      smearing = smearing + pub_dmft_smear_eta / ( &
           internal_dexpc(-(ww-pub_dmft_smear_w)/pub_dmft_smear_T) + 1.0_DP )
      smearing = smearing + pub_dmft_smear_eta / ( internal_dexpc( (ww + &
           pub_dmft_smear_w)/pub_dmft_smear_T) + 1.0_DP )

      energy = cmplx( (fermi_e(is) + pub_dmft_chem_shift) + realfrequ, &
           smearing, kind = DP)
      if(ien == pub_dmft_points .or. ien == pub_dmft_points-1) energy = cmplx( &
           pub_dmft_free_green_frequ, smearing, kind = DP)
      if(ien == pub_dmft_points-2) energy = cmplx( fermi_e(is) + &
           pub_dmft_chem_shift, 0.0_DP, kind = DP)
      if(ien == pub_dmft_points-3) energy = (1.0_DP, 0.0_DP)
      if(ien == pub_dmft_points-6) energy = (1.0_DP, 0.0_DP)
      if(ien == pub_dmft_points-4) energy = (0.0_DP, 0.0_DP)
      if(ien == pub_dmft_points-5) energy = (0.0_DP, 0.0_DP)
      if(ien == 1) energy = cmplx( -pub_dmft_free_green_frequ, smearing, kind &
           = DP)
   endif

   if (pub_debug_on_root) then
      if (aimag(energy) >= 0.0_DP) then
         write(stdout, '(a,f10.5," +",f10.5,"i")') &
              "DEBUG: energy = ", real(energy), aimag(energy)
      else
         write(stdout, '(a,f8.5," -",f8.5,"i")') &
              "DEBUG: energy = ", real(energy), -aimag(energy)
      end if
   end if

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_define_energy'

 end subroutine

 subroutine dmft_build_projected_green(greenf, h_atoms, hub, hub_proj_basis, &
      cell, elements, energy, fermi_e, num_spins, is, ien, ikpoint, ffirst, &
      norm_kpoints, nkpoints, kp_file_shift, en_step_over_pi, energyzero, &
      rot_vec_angles, rot_vec_angles_split, splitk_first, split_wgv, cluster, &
      connection_table, dimer_table, merge_neigh, merge_table, nmerge, &
      one_k_one_file)

   !---------------------------------------------------------------------------!
   ! Master routine to drive the calculation of the projected Green's function !
   ! for every Hubbard subspace                                                !
   !---------------------------------------------------------------------------!

   use parallel_strategy, only: par => pub_par
   use simulation_cell,   only: CELL_INFO
   use comms,             only: comms_barrier, pub_my_proc_id, pub_on_root
   use rundat,            only: pub_debug, pub_debug_on_root
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use ion,               only: ELEMENT
   use utils,             only: utils_alloc_check, utils_dealloc_check

   implicit none

   ! Arguments
   type(DMFT_MATRIX),               intent(in   ) :: greenf
   type(ARRAY_OF_MATRICES),         intent(in   ) :: h_atoms(par%nat_hub)
   type(HUBBARD_MODEL),             intent(in   ) :: hub
   type(FUNC_BASIS),                intent(in   ) :: hub_proj_basis
   type(CELL_INFO),                 intent(in   ) :: cell
   type(ELEMENT),                   intent(in   ) :: elements(:)
   complex(kind = DP),              intent(in   ) :: energy
   real(kind = DP),                 intent(in   ) :: fermi_e(2)
   integer,                         intent(in   ) :: num_spins
   integer,                         intent(in   ) :: is
   integer,                         intent(in   ) :: ien
   integer,                         intent(in   ) :: ikpoint
   integer,                         intent(in   ) :: ffirst
   real(kind = DP), allocatable,    intent(in   ) :: norm_kpoints(:)
   integer,                         intent(in   ) :: nkpoints
   integer,                         intent(in   ) :: kp_file_shift
   real(kind = DP),                 intent(in   ) :: en_step_over_pi
   integer,                         intent(in   ) :: energyzero
   real(kind = DP), allocatable,    intent(in   ) :: rot_vec_angles(:,:,:)
   real(kind = DP), allocatable,    intent(in   ) :: rot_vec_angles_split(:,:,:)
   integer,                         intent(in   ) :: splitk_first
   complex(kind = DP), allocatable, intent(in   ) :: split_wgv(:,:)
   logical,                         intent(in   ) :: cluster
   logical, allocatable,            intent(in   ) :: connection_table(:,:)
   integer, allocatable,            intent(in   ) :: dimer_table(:,:)
   integer, allocatable,            intent(in   ) :: merge_neigh(:)
   integer, allocatable,            intent(in   ) :: merge_table(:,:)
   integer,                         intent(in   ) :: nmerge
   logical,                         intent(in   ) :: one_k_one_file

   ! Local variables
   integer                         :: ierr
   integer                         :: i
   complex(kind = DP), allocatable :: site_greenf_buffer(:,:,:,:)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_build_projected_green'

   if(.not.pub_dmft_split) call comms_barrier

   if(pub_dmft_split)then
      i = par%nat_hub
   else
      i = par%num_hub_atoms_on_proc(pub_my_proc_id)
   endif

   allocate(site_greenf_buffer(i, i, 12, 12), stat = ierr)
   call utils_alloc_check('dmft_build_projected_green', 'site_greenf_buffer', &
        ierr)

   site_greenf_buffer(:, :, :, :) = 0.0_DP

   ! Looping over Hubbard atoms on my proc...
   call dmft_extract_green_function_elements(greenf, h_atoms, hub, &
        hub_proj_basis, cell, elements, energy, fermi_e, is, ien, ffirst, &
        ikpoint, norm_kpoints, en_step_over_pi, kp_file_shift, energyzero, &
        rot_vec_angles, rot_vec_angles_split, splitk_first, split_wgv, &
        site_greenf_buffer, cluster, connection_table, one_k_one_file)

   call dmft_write_green_function(h_atoms, hub, hub_proj_basis, num_spins, &
        energy, fermi_e, is, ien, ikpoint, nkpoints, norm_kpoints, ffirst, &
        kp_file_shift, splitk_first, site_greenf_buffer, cluster, dimer_table, &
        merge_neigh, merge_table, nmerge, rot_vec_angles, &
        rot_vec_angles_split, one_k_one_file)

   deallocate(site_greenf_buffer, stat = ierr)
   call utils_dealloc_check('dmft_build_projected_green', &
        'site_greenf_buffer', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_build_projected_green'

 end subroutine

 subroutine dmft_upfold_self_energy(self_energy, hub_overlap, nkpoints, &
      kpoint_data, same_self_for_all, icall)

   !---------------------------------------------------------------------------!
   ! Upfold the self energy to the NGWF basis                                  !
   !---------------------------------------------------------------------------!

   use sparse, only: sparse_convert, sparse_create, sparse_destroy, &
        sparse_num_cols, sparse_num_rows, sparse_product, sparse_show_matrix
   use utils,  only: utils_abort, utils_alloc_check, utils_dealloc_check
   use rundat, only: pub_debug, pub_debug_on_root

   implicit none

   ! Arguments
   type(DMFT_MATRIX),          intent(inout) :: self_energy
   type(DMFT_OVERLAP),         intent(in   ) :: hub_overlap
   integer,                    intent(in   ) :: nkpoints
   type(DMFT_KPOINT_MATRICES), intent(inout) :: kpoint_data
   logical,                    intent(in   ) :: same_self_for_all
   integer, optional,          intent(in   ) :: icall

   ! Local variables
   integer                         :: nn
   complex(kind = DP), allocatable :: mat_square(:,:)
   logical                         :: skip_sigma_imp_calc
   type(SPAM3)                     :: selfK1
   type(SPAM3)                     :: selfK2
   type(SPAM3)                     :: selfK3
   type(SPAM3)                     :: self_energy_vK
   integer                         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_upfold_self_energy'

   skip_sigma_imp_calc = .false.

   nn = sparse_num_rows(hub_overlap%tmatc)

   if(nkpoints > 1)then

      if(same_self_for_all)then
         if(.not.present(icall))then
            call utils_abort('ERROR icall should be present')
         endif
         if(pub_dmft_use_gpu_partially .and. pub_dmft_mpi_size > 1 .and. &
              .not.pub_dmft_splitk)then
            allocate(mat_square(nn, nn), stat = ierr)
            call utils_alloc_check('dmft_upfold_self_energy', 'mat_square', &
                 ierr)
            if( pub_dmft_mpi_rank /= 1)then
               mat_square = 0.0_DP
               skip_sigma_imp_calc = .true.
            endif
         endif
      endif

      call sparse_create(selfK1, kpoint_data%unrot_inv_overlap, iscmplx = &
           .true.)
      call sparse_create(selfK2, kpoint_data%unrot_inv_overlap, iscmplx = &
           .true.)
      call sparse_create(selfK3, kpoint_data%unrot_inv_overlap, iscmplx = &
           .true.)
      call sparse_create(self_energy_vK, hub_overlap%matc, selfK2, iscmplx = &
           .true.)

      ! Sigma_imp - > unrot_inv_overlap x hub_overlap x Sigma_imp x
      ! hub_overlap x unrot_inv_overlap
      if (.not. skip_sigma_imp_calc) then


         ! #ifndef GPU_SPEEDUP_MATMUL
         call sparse_product( selfK1, kpoint_data%hub_overlap, &
              kpoint_data%unrot_inv_overlap ) ! K1 = hub_overlap x
         ! unrot_inv_overlap
         if (pub_debug) then
            write(stdout, *) 'SHOWING kpoint hub_overlap x O_k'
            call sparse_show_matrix(selfK1, show_elems = .true.)
         end if
         call sparse_product( selfK2, self_energy%w_mat_v, selfK1 ) ! K2 =
         ! Sigma_imp x hub_overlap x unrot_inv_overlap
         call sparse_product( selfK1, kpoint_data%unrot_inv_overlap, &
              kpoint_data%hub_overlap ) ! K1 = unrot_inv_overlap x hub_overlap
         call sparse_product( selfK3, selfK1, selfK2 ) ! K3 = K1 x K2
         ! #else
         ! if(pub_dmft_use_gpu_onlyme) then
         !    call dmft_matmul_in_dense_format( selfK1, kpoint_data%hub_overlap &
         !         , kpoint_data%unrot_inv_overlap ) ! K1 = hub_overlap x
         !    ! unrot_inv_overlap
         !    if (pub_debug) then
         !       write(stdout, *) 'SHOWING hub_overlap x O_k'
         !       call sparse_show_matrix(selfK1, show_elems = .true.)
         !    end if
         !    call dmft_matmul_in_dense_format( selfK2, self_energy%w_mat_v, &
         !         selfK1 ) ! K2 = Sigma_imp x hub_overlap x unrot_inv_overlap
         !    call dmft_matmul_in_dense_format( selfK1, &
         !         kpoint_data%unrot_inv_overlap, kpoint_data%hub_overlap ) !
         !    ! K1 = unrot_inv_overlap x hub_overlap
         !    call dmft_matmul_in_dense_format( selfK3, selfK1, selfK2 ) ! K3 =
         !    ! K1 x K2
         ! else
         !    call sparse_product( selfK1, kpoint_data%hub_overlap, &
         !         kpoint_data%unrot_inv_overlap ) ! K1 = hub_overlap x
         !    ! unrot_inv_overlap
         !    if (pub_debug) then
         !       write(stdout, *) 'SHOWING hub_overlap x O_k'
         !       call sparse_show_matrix(selfK1, show_elems = .true.)
         !    end if
         !    call sparse_product( selfK2, self_energy%w_mat_v, selfK1 ) ! K2 =
         !    ! Sigma_imp x hub_overlap x unrot_inv_overlap
         !    call sparse_product( selfK1, kpoint_data%unrot_inv_overlap, &
         !         kpoint_data%hub_overlap ) ! K1 = unrot_inv_overlap x
         !    ! hub_overlap
         !    call sparse_product( selfK3, selfK1, selfK2 ) ! K3 = K1 x K2
         ! end if
         ! #endif
      end if
      if(same_self_for_all)then
         if(pub_dmft_mpi_size > 1 .and. pub_dmft_use_gpu_partially .and. &
              .not.pub_dmft_splitk)then
            if (pub_dmft_mpi_rank == 1) then
               call sparse_convert(mat_square, selfK3)
            endif
            call utils_abort("Reached point in code where nfs sync used to &
                 &appear; the code needs to be updated to use these settings.")
            ! call sync_via_nfs(mat_square, 'upfold_self_imp' // &
            !      trim(adjustl(utils_int_to_str(icall))), nn)
            call sparse_convert(selfK3, mat_square)
            deallocate(mat_square, stat = ierr)
            call utils_dealloc_check('dmft_upfold_self_energy', 'mat_square', &
                 ierr)
         endif
      endif

      ! #ifndef GPU_SPEEDUP_MATMUL
      call sparse_product(self_energy_vK, hub_overlap%matc, selfK3)
      call sparse_product(kpoint_data%self_energy, self_energy_vK, &
           hub_overlap%tmatc)
      ! #else
      ! if(pub_dmft_use_gpu_onlyme)then
      !    call dmft_matmul_in_dense_format(self_energy_vK, hub_overlap%matc, &
      !         selfK3)
      !    call dmft_matmul_in_dense_format(kpoint_data%self_energy, &
      !         self_energy_vK, hub_overlap%tmatc)
      ! else
      !    call sparse_product(self_energy_vK, hub_overlap%matc, selfK3)
      !    call sparse_product(kpoint_data%self_energy, self_energy_vK, &
      !         hub_overlap%tmatc)
      ! endif
      ! #endif

      if (pub_debug) then
         write(stdout, *) 'WARNING : checking : k-point self-energy'
         call dmft_print_spam3_info(kpoint_data%self_energy)
         write(stdout, *) 'WARNING : checking : self_energy'
         call dmft_print_spam3_info(self_energy%mat)
      end if

      call sparse_destroy(selfK1)
      call sparse_destroy(selfK2)
      call sparse_destroy(selfK3)
      call sparse_destroy(self_energy_vK)
   else
      ! GAMMA ONLY VERSION

      ! #ifndef GPU_SPEEDUP_MATMUL
      call sparse_product(self_energy%mat_v, hub_overlap%matc, &
           self_energy%w_mat_v)
      ! #else
      ! if(pub_dmft_use_gpu_onlyme)then
      !    call dmft_matmul_in_dense_format(self_energy%mat_v, hub_overlap%matc, &
      !         self_energy%w_mat_v)
      ! else
      !    call sparse_product(self_energy%mat_v, hub_overlap%matc, &
      !         self_energy%w_mat_v)
      ! endif
      ! #endif

      ! #ifndef GPU_SPEEDUP_MATMUL
      call sparse_product(self_energy%mat, self_energy%mat_v, hub_overlap%tmatc)
      ! #else
      ! if(pub_dmft_use_gpu_onlyme)then
      !    call dmft_matmul_in_dense_format(self_energy%mat, self_energy%mat_v, &
      !         hub_overlap%tmatc)
      ! else
      !    call sparse_product(self_energy%mat, self_energy%mat_v, &
      !         hub_overlap%tmatc)
      ! endif
      ! #endif

   end if

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_upfold_self_energy'

 end subroutine

 subroutine dmft_construct_self_energy(ien, kien, is, self_energy, hub_atom, &
      hub_atomb, connection_table, nmerge, merge_table, merge_neigh, &
      dimer_table, h_atoms, hub, hub_proj_basis, rot_vec_angles, ffirst, &
      llast, cluster, show_matrices, reopen_sigma_file_each_step, &
      force_reopen_file, scratchdir)

   !---------------------------------------------------------------------------!
   ! Constructs the projected self energy                                      !
   !---------------------------------------------------------------------------!

   use sparse,            only: sparse_scale, sparse_show_matrix
   use comms,             only: comms_barrier, pub_my_proc_id, pub_on_root
   use rundat,            only: pub_debug, pub_debug_on_root
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par
   use utils,             only: utils_alloc_check, utils_dealloc_check

   implicit none

   ! Arguments
   integer,                      intent(in   ) :: ien
   integer,                      intent(in   ) :: kien
   integer,                      intent(in   ) :: is
   type(DMFT_MATRIX),            intent(inout) :: self_energy
   integer,                      intent(inout) :: hub_atom
   integer,                      intent(inout) :: hub_atomb
   logical, allocatable,         intent(in   ) :: connection_table(:,:)
   integer,                      intent(in   ) :: nmerge
   integer, allocatable,         intent(in   ) :: merge_table(:,:)
   integer, allocatable,         intent(in   ) :: merge_neigh(:)
   integer, allocatable,         intent(in   ) :: dimer_table(:,:)
   type(ARRAY_OF_MATRICES),      intent(in   ) :: h_atoms(par%nat_hub)
   type(HUBBARD_MODEL),          intent(in   ) :: hub
   type(FUNC_BASIS),             intent(in   ) :: hub_proj_basis
   real(kind = DP), allocatable, intent(in   ) :: rot_vec_angles(:,:,:)
   integer,                      intent(in   ) :: ffirst
   integer,                      intent(in   ) :: llast
   logical,                      intent(in   ) :: cluster
   logical,                      intent(in   ) :: show_matrices
   logical,                      intent(in   ) :: reopen_sigma_file_each_step
   logical,                      intent(in   ) :: force_reopen_file
   character(2000),              intent(in   ) :: scratchdir

   ! Local variables
   complex(kind = DP), allocatable :: site_self_energy_buffer(:,:,:,:)
   integer                         :: i
   integer                         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_construct_self_energy'

   i = par%num_hub_atoms_on_proc(pub_my_proc_id)

   ! if(.not.allocated(site_self_energy_buffer)) then
   !    allocate(site_self_energy_buffer(i, i, 12, 12), stat = ierr)
   !    call utils_alloc_check('dmft_construct_self_energy', &
   !         'site_self_energy_buffer', ierr)
   ! end if
   ! if(size(site_self_energy_buffer, 1) /= i)then
   !    deallocate(site_self_energy_buffer, stat = ierr)
   !    call utils_dealloc_check('dmft_construct_self_energy', &
   !         'site_self_energy_buffer', ierr)
   !    allocate(site_self_energy_buffer(i, i, 12, 12), stat = ierr)
   !    call utils_alloc_check('dmft_construct_self_energy', &
   !         'site_self_energy_buffer', ierr)
   ! endif
   allocate(site_self_energy_buffer(i, i, 12, 12), stat = ierr)
   call utils_alloc_check('dmft_construct_self_energy', &
        'site_self_energy_buffer', ierr)

   site_self_energy_buffer = 0.0_DP
   call dmft_build_self_energy_from_extended_atoms(is, ien, kien, hub_atom, &
        hub_atomb, site_self_energy_buffer, nmerge, merge_table, merge_neigh, &
        dimer_table, hub, hub_proj_basis, ffirst, llast, cluster, &
        reopen_sigma_file_each_step, force_reopen_file, scratchdir)
   call comms_barrier

   ! site_self_energy_buffer is defined from previous section
   call dmft_distribute_self_energy(site_self_energy_buffer, self_energy, &
        connection_table, hub, hub_proj_basis, is, rot_vec_angles, h_atoms, &
        cluster)

   deallocate(site_self_energy_buffer, stat = ierr)
   call utils_dealloc_check('dmft_construct_self_energy', &
        'site_self_energy_buffer', ierr)
   call comms_barrier

   if(pub_debug_on_root) then
      if (show_matrices .and. ien == 4) then
         ! if (pub_on_root) write(stdout, *) '-'*80
         if(pub_debug_on_root) then
            write(stdout, *) 'DEBUG: self energy in atomic basis:'
            call sparse_scale(self_energy%w_mat_v, (1000.0_DP, 0.0_DP))
            call sparse_show_matrix(self_energy%w_mat_v, show_elems = .true.)
            call sparse_scale(self_energy%w_mat_v, &
                 CMPLX(1.0_DP/1000.0_DP, 0.0_DP, kind=DP))
         end if
         ! if (pub_on_root) write(stdout, *) '-'*80
      endif
   endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_construct_self_energy'

 end subroutine

 subroutine internal_renorm_Z(mat, i)

   !---------------------------------------------------------------------------!
   ! Renormalises Z, the quasiparticle weights                                 !
   !---------------------------------------------------------------------------!

   use rundat, only: pub_debug_on_root

   implicit none

   ! Arguments
   real(kind = DP), allocatable, intent(inout) :: mat(:,:,:,:)
   integer,                      intent(in   ) :: i

   ! Local variables
   integer :: k, l, m

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &internal_renorm_Z'

   do k = 1, size(mat, 1)
      do l = 1, size(mat, 2)
         do m = 1, size(mat, 3)
            mat(k, l, m, i) = 1.0_DP/(1.0_DP-mat(k, l, m, i))
            if(mat(k, l, m, i) > 1.0_DP) mat(k, l, m, i) = 0.0_DP
            if(mat(k, l, m, i) < 0.0_DP) mat(k, l, m, i) = 0.0_DP
         enddo
      enddo
   enddo
   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &internal_renorm_Z'

 end subroutine

 subroutine internal_one_minus_Z(mat, i)

   !---------------------------------------------------------------------------!
   ! Calculates one minus the quasiparticle weights                            !
   !---------------------------------------------------------------------------!

   use rundat, only: pub_debug_on_root

   implicit none

   ! Arguments
   real(kind = DP), allocatable, intent(inout) :: mat(:,:,:,:)
   integer,                      intent(in   ) :: i

   ! Local variables
   integer :: k, l, m

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &internal_one_minus_Z'

   do k = 1, size(mat, 1)
      do l = 1, size(mat, 2)
         do m = 1, size(mat, 3)
            mat(k, l, m, i) = 1.0_DP-mat(k, l, m, i)
         enddo
      enddo
   enddo
   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &internal_one_minus_Z'

 end subroutine

 logical function dmft_dimer_connected(loop, dimer_table)

   !---------------------------------------------------------------------------!
   ! Checks if two Hubbard atoms are connected                                 !
   !---------------------------------------------------------------------------!

   use utils, only: utils_assert

   implicit none

   ! Arguments
   type(DMFT_LOOP_INFO), intent(in) :: loop
   integer, allocatable, intent(in) :: dimer_table(:,:)

   dmft_dimer_connected = dimer_table(loop%hat, 2) == loop%hatb
   call utils_assert((dmft_hubbard_on_my_proc(loop%hat) /= 0 .and. &
        dmft_hubbard_on_my_proc(loop%hatb) /= 0), 'Error in &
        &dmft_dimer_connected: connected atoms are not on my proc')

 end function

 integer function dmft_hubbard_on_my_proc(j)

   !---------------------------------------------------------------------------!
   ! Checks if a given Hubbard atom belongs to this proc                       !
   !---------------------------------------------------------------------------!

   use comms,             only: pub_my_proc_id
   use parallel_strategy, only: par => pub_par
   use rundat,            only: pub_debug_on_root

   implicit none

   ! Arguments
   integer, intent(in) :: j

   ! Local variables
   integer :: i, k, hub_atom, hat

   dmft_hubbard_on_my_proc = 0
   do hub_atom = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)
      hat = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
      if(hat == j)then
         dmft_hubbard_on_my_proc = hub_atom
         if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
              &dmft_hubbard_on_my_proc'
         return
      endif
   enddo

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_hubbard_on_my_proc'

 end function

 subroutine dmft_build_real_space_quantities(sigma_density_fine, &
      green_density_fine, green_density_fine_proj, green_tot_fine, &
      Z_density_fine, sigma_density_spam, green_density_spam, &
      green_density_spam_proj, green_tot_spam, Z_density_spam, greenf, &
      self_energy, hub_overlap, overlap_cplx, mdl, nkpoints, is, ien, &
      en_start, en_step, nspecial_frequ, inv_projectors, real_space_exit)

   !---------------------------------------------------------------------------!
   ! Constructs various real space quantities for a given spin and frequency,  !
   ! in preparation for them to be plotted                                     !
   !---------------------------------------------------------------------------!

   use constants,  only: HARTREE_IN_EVS
   use visual,     only: visual_scalarfield
   use sparse, only: sparse_axpy, sparse_copy, sparse_create, sparse_destroy, &
        sparse_product, sparse_scale
   use utils,      only: utils_assert, utils_alloc_check
   use rundat, only: pub_debug_on_root, pub_dmft_complex_freq, &
        pub_dmft_plot_real_space, &
        pub_dmft_points, pub_dmft_temp, pub_num_spins
   use model_type, only: MODEL

   implicit none

   ! Arguments
   real(kind = DP), allocatable, intent(inout) :: sigma_density_fine(:,:,:,:)
   real(kind = DP), allocatable, intent(inout) :: green_density_fine(:,:,:,:)
   real(kind = DP), allocatable, intent(inout) :: &
                                                  green_density_fine_proj(:,:,:,:)
   real(kind = DP), allocatable, intent(inout) :: green_tot_fine(:,:,:,:)
   real(kind = DP), allocatable, intent(inout) :: Z_density_fine(:,:,:,:)
   type(SPAM3),                  intent(inout) :: sigma_density_spam(2)
   type(SPAM3),                  intent(inout) :: green_density_spam(2, 1)
   type(SPAM3),                  intent(inout) :: green_density_spam_proj(2)
   type(SPAM3),                  intent(inout) :: green_tot_spam(2)
   type(SPAM3),                  intent(inout) :: Z_density_spam(2)
   type(DMFT_MATRIX),            intent(inout) :: greenf
   type(DMFT_MATRIX),            intent(inout) :: self_energy
   type(DMFT_OVERLAP),           intent(in   ) :: hub_overlap
   type(SPAM3),                  intent(in   ) :: overlap_cplx
   type(MODEL),                  intent(in   ) :: mdl
   integer,                      intent(in   ) :: nkpoints
   integer,                      intent(in   ) :: is
   integer,                      intent(in   ) :: ien
   real(kind = DP),              intent(in   ) :: en_start
   real(kind = DP),              intent(in   ) :: en_step
   integer,                      intent(in   ) :: nspecial_frequ
   logical,                      intent(in   ) :: inv_projectors
   logical,                      intent(inout) :: real_space_exit

   ! Local variables
   integer     :: i
   integer     :: j
   type(SPAM3) :: tmp_spam_complex
   integer     :: u(1), ii, ierr
   integer     :: funit

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_build_real_space_quantities'

   if(.not. pub_dmft_plot_real_space) then
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_build_real_space_quantities as we are not plotting real space&
           & quantities'
      return
   end if

   if(      nkpoints > 1              ) then
      write(stdout, *) 'ERROR: skipping real_space_quantities plots due to k &
           &points'
      write(stdout, *) ' If you want to plot real space quantities, you will &
           &have to average over k'
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_build_real_space_quantities'
      return
   endif

   call utils_assert(.not. pub_dmft_split, "Error in &
        &dmft_plotting_real_space_quantities: plotting is incompatible with &
        &split")

   !do not use energyzero here...check below...!
   u = minloc ( abs ( (/( en_start + real(i-1, kind = DP)*en_step, i = 1, &
        pub_dmft_points )/)) )
   ii = u(1)

   if(ien == 1 .and. pub_dmft_integrate_green)then
      call sparse_scale(green_tot_spam(is), 0.0_DP)
   endif

   if(((ien == 1 .or. ien == 2 .or. (pub_dmft_integrate_green .and. ien <= &
        pub_dmft_points-nspecial_frequ)) .and. pub_dmft_temp > 0.0_DP) .or. &
        (.not. pub_dmft_complex_freq .and. (ien == ii .or. ien == ii + 1 .or. &
        (pub_dmft_integrate_green .and. ien <= ii .and. ien > 1)) ) ) then

      real_space_exit = ien == 2 .and. .not.pub_dmft_integrate_green .and. &
           pub_dmft_temp > 0.0_DP
      real_space_exit = real_space_exit .or. (ien == &
           pub_dmft_points-nspecial_frequ .and. pub_dmft_temp > 0.0_DP)
      real_space_exit = real_space_exit .or. (ien == ii + 1 .and. .not. &
           pub_dmft_complex_freq)

      if(.not.allocated(sigma_density_fine))then
         allocate(sigma_density_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
              mdl%fine_grid%max_slabs12, pub_num_spins), stat = ierr)
         call utils_alloc_check('dmft_build_real_space_quantities', &
              'sigma_density_fine', ierr)
      endif
      if(.not.allocated(green_density_fine))then
         allocate(green_density_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
              mdl%fine_grid%max_slabs12, pub_num_spins), stat = ierr)
         call utils_alloc_check('dmft_build_real_space_quantities', &
              'green_density_fine', ierr)
      endif
      if(.not.allocated(green_density_fine_proj))then
         allocate(green_density_fine_proj(mdl%fine_grid%ld1, &
              mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins), &
              stat = ierr)
         call utils_alloc_check('dmft_build_real_space_quantities', &
              'green_density_fine_proj', ierr)
      endif
      if(.not.allocated(Z_density_fine))then
         allocate(Z_density_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
              mdl%fine_grid%max_slabs12, pub_num_spins), stat = ierr)
         call utils_alloc_check('dmft_build_real_space_quantities', &
              'Z_density_fine', ierr)
      endif
      if(.not.allocated(green_tot_fine))then
         allocate(green_tot_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
              mdl%fine_grid%max_slabs12, pub_num_spins), stat = ierr)
         call utils_alloc_check('dmft_build_real_space_quantities', &
              'green_tot_fine', ierr)
      endif

      if(pub_dmft_integrate_green)then
         call sparse_create(tmp_spam_complex, greenf%mat, iscmplx = .true.)
         call sparse_copy(tmp_spam_complex, greenf%mat)
         if(pub_dmft_temp > 0.0_DP)then
            if(ien < pub_dmft_points-nspecial_frequ)then
               call sparse_axpy(green_tot_spam(is), tmp_spam_complex, &
                    2.0_DP*abs(pub_dmft_temp))
            endif
            if(ien == pub_dmft_points-nspecial_frequ)then
               call sparse_axpy(green_tot_spam(is), tmp_spam_complex, &
                    2.0_DP*abs(pub_dmft_temp))
               call sparse_scale(tmp_spam_complex, &
                    cmplx(0.0_DP, &
                    0.5_DP*(en_start + real(ien-1, kind = DP)*en_step), &
                    kind = DP))
               call sparse_axpy(green_tot_spam(is), tmp_spam_complex, 1.0_DP)
               call sparse_copy(tmp_spam_complex, green_tot_spam(is))
               write(stdout, *) 'OCCUPATION SPIN : ', is
               call dmft_matmul_in_dense_format(tmp_spam_complex, &
                    tmp_spam_complex, overlap_cplx, show_trace = .true.)
               if(is == 2)then
                  call sparse_copy(tmp_spam_complex, green_tot_spam(1))
                  call sparse_axpy(tmp_spam_complex, green_tot_spam(2), -1.0_DP)
                  call sparse_scale(tmp_spam_complex, (0.5_DP, 0.0_DP))
                  write(stdout, *) 'INTEGRATED MAGNETISM'
                  call dmft_matmul_in_dense_format(tmp_spam_complex, &
                       tmp_spam_complex, overlap_cplx, show_trace = .true.)
               endif
            endif
         endif
         if(.not. pub_dmft_complex_freq)then
            if(ien <= ii .and. ien > 1)then
               call sparse_scale(tmp_spam_complex, &
                    cmplx(0.0_DP, en_step/dacos(-1.0_DP), kind = DP))
               call sparse_axpy(green_tot_spam(is), tmp_spam_complex, 1.0_DP)
            endif
            if(ien == ii)then
               call sparse_copy(tmp_spam_complex, green_tot_spam(is))
               write(stdout, *) 'OCCUPATION SPIN : ', is
               call dmft_matmul_in_dense_format(tmp_spam_complex, &
                    tmp_spam_complex, overlap_cplx, show_trace = .true.)
               if(is == 2)then
                  call sparse_copy(tmp_spam_complex, green_tot_spam(1))
                  call sparse_axpy(tmp_spam_complex, green_tot_spam(2), -1.0_DP)
                  call sparse_scale(tmp_spam_complex, (0.5_DP, 0.0_DP))
                  write(stdout, *) 'INTEGRATED MAGNETISM'
                  call dmft_matmul_in_dense_format(tmp_spam_complex, &
                       tmp_spam_complex, overlap_cplx, show_trace = .true.)
               endif
            endif
         endif
         call sparse_destroy(tmp_spam_complex)
      endif

      if( (ien == 1 .and. pub_dmft_temp > 0.0_DP) .or. (ien == ii .and. .not. &
           pub_dmft_complex_freq) )then
         !ImGreen
         call sparse_create(tmp_spam_complex, greenf%mat, iscmplx = .true.)
         call sparse_copy(tmp_spam_complex, greenf%mat)
         call sparse_scale(tmp_spam_complex, &
              cmplx(0.0_DP, 1.0_DP/dacos(-1.0_DP), kind = DP))
         call sparse_copy(green_density_spam(is, PUB_1K), tmp_spam_complex)
         call sparse_destroy(tmp_spam_complex)
         if(inv_projectors)then
            !ImGreen_proj
            call sparse_create(tmp_spam_complex, greenf%mat_v, &
                 hub_overlap%inv_tmatc, iscmplx = .true.)
            call sparse_product(greenf%mat_v, hub_overlap%matc, greenf%w_mat_v)
            call sparse_product(tmp_spam_complex, greenf%mat_v, &
                 hub_overlap%inv_tmatc)
            call sparse_scale(tmp_spam_complex, &
                 cmplx(0.0_DP, 1.0_DP/dacos(-1.0_DP), kind = DP))
            call sparse_copy(green_density_spam_proj(is), tmp_spam_complex)
            call sparse_destroy(tmp_spam_complex)
            !ImSelf
            call sparse_create(tmp_spam_complex, self_energy%mat_v, &
                 hub_overlap%inv_tmatc, iscmplx = .true.)
            call sparse_product(self_energy%mat_v, hub_overlap%inv_matc, &
                 self_energy%mat_v)
            call sparse_product(tmp_spam_complex, self_energy%mat_v, &
                 hub_overlap%inv_tmatc)
            if(pub_dmft_temp > 0.0_DP)then
               call sparse_scale(tmp_spam_complex, (0.0_DP, -1.0_DP))
            endif
                        ! Calculate the projected Green function for each
                        ! Hubbard atom
            call sparse_copy(sigma_density_spam(is), tmp_spam_complex)
            call sparse_destroy(tmp_spam_complex)
         endif
      endif

      if( (ien == 2 .and. pub_dmft_temp > 0.0_DP) .or. (ien == ii + 1 .and. &
           .not. pub_dmft_complex_freq) )then
         if(inv_projectors)then
            call sparse_create(tmp_spam_complex, self_energy%mat_v, &
                 hub_overlap%inv_tmatc, iscmplx = .true.)
            call sparse_product(self_energy%mat_v, hub_overlap%inv_matc, &
                 self_energy%w_mat_v)
            call sparse_product(tmp_spam_complex, self_energy%mat_v, &
                 hub_overlap%inv_tmatc)
            if(pub_dmft_temp > 0.0_DP)then
               call sparse_scale(tmp_spam_complex, (0.0_DP, -1.0_DP))
            endif
            call sparse_copy(Z_density_spam(is), tmp_spam_complex)
            call sparse_destroy(tmp_spam_complex)
         endif
      endif
   endif
   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_build_real_space_quantities'

 end subroutine

 subroutine dmft_plot_real_space_quantities(sigma_density_fine, &
      green_density_fine, green_density_fine_proj, green_tot_fine, &
      Z_density_fine, sigma_density_spam, green_density_spam, &
      green_density_spam_proj, green_tot_spam, Z_density_spam, overlap, &
      ngwf_basis, rep, mdl, elements, nkpoints, en_start, en_step, &
      real_space_exit)

   !---------------------------------------------------------------------------!
   ! Plots various real space quantities                                       !
   !---------------------------------------------------------------------------!

   use constants,           only: HARTREE_IN_EVS
   use density,             only: density_on_grid
   use function_basis,      only: FUNC_BASIS
   use ion,                 only: ELEMENT
   use model_type,          only: MODEL
   use ngwf_representation, only: NGWF_REP
   use visual,              only: visual_scalarfield
   use sparse,              only: sparse_copy, sparse_destroy
   use utils,               only: utils_assert, utils_unit, utils_int_to_str, &
                                  utils_dealloc_check
   use comms,               only: pub_on_root
   use rundat, only: pub_debug_on_root, &
        pub_dmft_paramagnetic, &
        pub_dmft_plot_real_space, &
        pub_dmft_points, pub_dmft_temp, pub_num_spins

   implicit none

   ! Arguments
   real(kind = DP), allocatable, intent(inout) :: sigma_density_fine(:,:,:,:)
   real(kind = DP), allocatable, intent(inout) :: green_density_fine(:,:,:,:)
   real(kind = DP), allocatable, intent(inout) :: &
                                                  green_density_fine_proj(:,:,:,:)
   real(kind = DP), allocatable, intent(inout) :: green_tot_fine(:,:,:,:)
   real(kind = DP), allocatable, intent(inout) :: Z_density_fine(:,:,:,:)
   type(SPAM3),                  intent(inout) :: sigma_density_spam(2)
   type(SPAM3),                  intent(inout) :: green_density_spam(2, 1)
   type(SPAM3),                  intent(inout) :: green_density_spam_proj(2)
   type(SPAM3),                  intent(inout) :: green_tot_spam(2)
   type(SPAM3),                  intent(inout) :: Z_density_spam(2)
   type(SPAM3),                  intent(in   ) :: overlap
   type(FUNC_BASIS),             intent(in   ) :: ngwf_basis
   type(NGWF_REP),               intent(in   ) :: rep
   type(MODEL),                  intent(in   ) :: mdl
   type(ELEMENT),                intent(in   ) :: elements(:)
   integer,                      intent(in   ) :: nkpoints
   real(kind = DP),              intent(in   ) :: en_start
   real(kind = DP),              intent(in   ) :: en_step
   logical,                      intent(inout) :: real_space_exit

   ! Local variables
   integer     :: i
   integer     :: j
   integer     :: is
   type(SPAM3) :: tmp_spam_complex
   integer     :: u(1), ii, ierr
   integer     :: funit

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_plot_real_space_quantities'

   if(.not. pub_dmft_plot_real_space) then
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_plot_real_space_quantities prematurely since &
           &dmft_plot_real_space : F'
      return
   end if
   if(nkpoints > 1) then
      write(stdout, *) 'ERROR skipping real_space_quantities plots due to k &
           &points'
      write(stdout, *) 'If you want to plot real space quantities, you will &
           &have to average over k'
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_plot_real_space_quantities'
      return
   endif

   call utils_assert(.not. pub_dmft_split, "Error in &
        &dmft_plotting_real_space_quantities: plotting is incompatible with &
        &split")

   ! cw: do not use energyzero here...check below...!
   u = minloc ( abs ( (/( en_start + real(i-1, kind = DP)*en_step, i = 1, &
        pub_dmft_points )/)) )
   ii = u(1)

   if(allocated(sigma_density_fine) .and. allocated(green_density_fine) .and. &
        allocated(Z_density_fine) .and. allocated(green_density_fine_proj) &
        .and. allocated(green_tot_fine) )then
      write(stdout, *) 'building real-space grid for density plots'
      if(pub_dmft_paramagnetic .and. pub_num_spins /= 1)then
         call sparse_copy( green_density_spam_proj(2), &
              green_density_spam_proj(1) )
         call sparse_copy( green_density_spam(2, PUB_1K), &
              green_density_spam(1, PUB_1K) )
         call sparse_copy( sigma_density_spam(2), sigma_density_spam(1) )
         call sparse_copy( Z_density_spam(2), Z_density_spam(1) )
         call sparse_copy( green_tot_spam(2), green_tot_spam(1) )
      endif
      write(stdout, *) '-------------------------------build &
           &density-------------------------------------'
      write(stdout, *) 'go for density plot, green_density'
      if(.not.pub_dmft_plot_real_space_sigma) then
         call density_on_grid(green_density_fine, mdl%fine_grid, mdl%dbl_grid, &
              mdl%cell, mdl%fftbox, green_density_spam, overlap, &
              rep%ngwfs_on_grid(1), ngwf_basis, rep%ngwfs_on_grid(1), ngwf_basis)
      end if
      write(stdout, *) 'go for density plot, green_density_proj'
      if(.not.pub_dmft_plot_real_space_sigma) then
         call density_on_grid(green_density_fine_proj, mdl%fine_grid, &
              mdl%dbl_grid, mdl%cell, mdl%fftbox, green_density_spam_proj, &
              overlap, rep%ngwfs_on_grid(1), ngwf_basis, rep%ngwfs_on_grid(1), &
              ngwf_basis)
      end if
      write(stdout, *) 'go for density plot, green_tot'
      if(.not.pub_dmft_plot_real_space_sigma .and. pub_dmft_integrate_green) &
           then
         call density_on_grid(green_tot_fine, mdl%fine_grid, mdl%dbl_grid, &
              mdl%cell, mdl%fftbox, green_tot_spam, overlap, &
              rep%ngwfs_on_grid(1), ngwf_basis, rep%ngwfs_on_grid(1), ngwf_basis)
      end if
      write(stdout, *) 'go for density plot, sigma_density'
      call density_on_grid(sigma_density_fine, mdl%fine_grid, mdl%dbl_grid, &
           mdl%cell, mdl%fftbox, sigma_density_spam, overlap, &
           rep%ngwfs_on_grid(1), ngwf_basis, rep%ngwfs_on_grid(1), ngwf_basis)
      write(stdout, *) 'go for density plot, Z_density'
      if(.not.pub_dmft_plot_real_space_sigma) then
         call density_on_grid(Z_density_fine, mdl%fine_grid, mdl%dbl_grid, &
              mdl%cell, mdl%fftbox, Z_density_spam, overlap, &
              rep%ngwfs_on_grid(1), ngwf_basis, rep%ngwfs_on_grid(1), ngwf_basis)
      end if
      write(stdout, *) '--------------------------done, now go for &
           &plots--------------------------------'
      if(pub_on_root)then
         funit = utils_unit()
         open(unit = funit, file = 'density_cut_G' // &
              trim(adjustl(utils_int_to_str(1))))
         write(funit, *) maxval(green_density_fine(:, :, 1, 1)), &
              minval(green_density_fine(:, :, 1, 1))
         do i = 1, size(green_density_fine, 1)
            do j = 1, size(green_density_fine, 2)
               write(funit, '(i5, i5, f10.5)') i, j, green_density_fine(i, j, &
                    1, 1)
            enddo
         enddo
         close(funit)
         open(unit = funit, file = 'density_cut_Sigma' // &
              trim(adjustl(utils_int_to_str(1))))
         write(funit, *) maxval(sigma_density_fine(:, :, 1, 1)), &
              minval(sigma_density_fine(:, :, 1, 1))
         do i = 1, size(sigma_density_fine, 1)
            do j = 1, size(sigma_density_fine, 2)
               write(funit, '(i5, i5, f10.5)') i, j, sigma_density_fine(i, j, &
                    1, 1)
            enddo
         enddo
         close(funit)
         open(unit = funit, file = 'density_cut_Z' // &
              trim(adjustl(utils_int_to_str(1))))
         write(funit, *) maxval(Z_density_fine(:, :, 1, 1)), &
              minval(Z_density_fine(:, :, 1, 1))
         do i = 1, size(Z_density_fine, 1)
            do j = 1, size(Z_density_fine, 2)
               write(funit, '(i5, i5, f10.5)') i, j, Z_density_fine(i, j, 1, 1)
            enddo
         enddo
         close(funit)
      endif
      do i = 1, pub_num_spins
         if(.not.pub_dmft_plot_real_space_sigma)then
            write(stdout, *) 'display scalarfield in cube format, go for &
                 &green_proj'
            call visual_scalarfield(green_density_fine_proj(:, :, :, i), &
                 mdl%fine_grid, mdl%cell, 'density_proj_fermi_level_dmft', &
                 'density_proj_fermi_level_dmft' // &
                 trim(adjustl(utils_int_to_str(i))), elements, 1.0_DP)
            write(stdout, *) 'display scalarfield in cube format, go for green'
            call visual_scalarfield(green_density_fine(:, :, :, i), &
                 mdl%fine_grid, mdl%cell, 'density_fermi_level_dmft_for_spin', &
                 'density_fermi_level_dmft_for_spin' // &
                 trim(adjustl(utils_int_to_str(i))), elements, 1.0_DP)
            write(stdout, *) 'display scalarfield in cube format, go for green &
                 &tot'
            if(pub_dmft_integrate_green) then
               call visual_scalarfield(green_tot_fine(:, :, :, i), &
                    mdl%fine_grid, mdl%cell, 'density_dmft_for_spin', &
                    'density_dmft_for_spin' // trim(adjustl(utils_int_to_str( &
                    i))), elements, 1.0_DP)
            end if
         endif
         write(stdout, *) 'display scalarfield in cube format, go for sigma'
         call visual_scalarfield(sigma_density_fine(:, :, :, i), &
              mdl%fine_grid, mdl%cell, 'scattering_fermi_level_dmft', &
              'scattering_fermi_level_dmft' // trim(adjustl(utils_int_to_str( &
              i))), elements, HARTREE_IN_EVS)
         if(.not.pub_dmft_plot_real_space_sigma)then
            call visual_scalarfield(abs(sigma_density_fine(:, :, :, i)), &
                 mdl%fine_grid, mdl%cell, 'abs_scattering_fermi_level_dmft', &
                 'abs_scattering_fermi_level_dmft' // &
                 trim(adjustl(utils_int_to_str(i))), elements, HARTREE_IN_EVS)
            !--------------------------------------!
            Z_density_fine(:, :, :, i) = (Z_density_fine(:, :, :, &
                 i)-sigma_density_fine(:, :, :, i))/en_step
            !--------------------------------------!
            write(stdout, *) 'display scalarfield in cube format, go for slope'
            call visual_scalarfield(Z_density_fine(:, :, :, i), mdl%fine_grid, &
                 mdl%cell, 'Slope_fermi_level_dmft', 'Slope_fermi_level_dmft' &
                 // trim(adjustl(utils_int_to_str(i))), elements, 1.0_DP)
            write(stdout, *) 'display scalarfield in cube format, go for abs &
                 &slope'
            call visual_scalarfield(abs(Z_density_fine(:, :, :, i)), &
                 mdl%fine_grid, mdl%cell, 'Abs_Slope_fermi_level_dmft', &
                 'Abs_Slope_fermi_level_dmft' // trim(adjustl( &
                 utils_int_to_str(i))), elements, 1.0_DP)
            !--------------------------------------!
            call internal_renorm_Z(Z_density_fine, i)
            !--------------------------------------!
            call visual_scalarfield(Z_density_fine(:, :, :, i), mdl%fine_grid, &
                 mdl%cell, 'Z_fermi_level_dmft', 'Z_fermi_level_dmft' // &
                 trim(adjustl(utils_int_to_str(i))), elements, 1.0_DP)
            !--------------------------------------!
            if(pub_dmft_temp > 0.0_DP)then
               Z_density_fine(:, :, :, i) = sigma_density_fine(:, :, :, &
                    i)/en_start
               call visual_scalarfield(Z_density_fine(:, :, :, i), &
                    mdl%fine_grid, mdl%cell, &
                    'Slope_first_only_fermi_level_dmft_for_spin', &
                    'Slope_first_only_fermi_level_dmft_for_spin' // &
                    trim(adjustl(utils_int_to_str(i))), elements, 1.0_DP)
               call internal_renorm_Z(Z_density_fine, i)
               call visual_scalarfield(Z_density_fine(:, :, :, i), &
                    mdl%fine_grid, mdl%cell, &
                    'Z_first_only_fermi_level_dmft_for_spin', &
                    'Z_first_only_fermi_level_dmft_for_spin' // &
                    trim(adjustl(utils_int_to_str(i))), elements, 1.0_DP)
               call internal_one_minus_Z(Z_density_fine, i)
               call visual_scalarfield(Z_density_fine(:, :, :, i), &
                    mdl%fine_grid, mdl%cell, &
                    '1mZ_first_only_fermi_level_dmft_for_spin', &
                    '1mZ_first_only_fermi_level_dmft_for_spin' // &
                    trim(adjustl(utils_int_to_str(i))), elements, 1.0_DP)
            endif
         endif
      enddo

      if(.not.pub_dmft_plot_real_space_sigma .and. .not.pub_dmft_paramagnetic &
           .and. pub_num_spins /= 1)then
         write(stdout, *) 'display scalarfield in cube format, go for magnetic &
              &density at Fermi level'
         call visual_scalarfield((green_density_fine_proj(:, :, :, &
              2)-green_density_fine_proj(:, :, :, 1))/2.0_DP, mdl%fine_grid, &
              mdl%cell, 'mag_density_proj_fermi_level_dmft', &
              'mag_density_proj_fermi_level_dmft', elements, 1.0_DP)
         write(stdout, *) 'display scalarfield in cube format, go for proj &
              &magnetic density at Fermi level'
         call visual_scalarfield((green_density_fine(:, :, :, &
              2)-green_density_fine(:, :, :, 1))/2.0_DP, mdl%fine_grid, &
              mdl%cell, 'mag_density_fermi_level_dmft_for_spin', &
              'mag_density_fermi_level_dmft_for_spin', elements, 1.0_DP)
         write(stdout, *) 'display scalarfield in cube format, go for magnetic &
              &density'
         if(pub_dmft_integrate_green) then
            call visual_scalarfield((green_tot_fine(:, :, :, &
                 2)-green_tot_fine(:, :, :, 1))/2.0_DP, mdl%fine_grid, mdl%cell, &
                 'mag_density_dmft_for_spin', 'mag_density_dmft_for_spin', &
                 elements, 1.0_DP)
         end if
      endif
      deallocate(sigma_density_fine, stat = ierr)
      call utils_dealloc_check('dmft_plot_real_space_quantities', &
           'sigma_density_fine', ierr)
      deallocate(green_density_fine, stat = ierr)
      call utils_dealloc_check('dmft_plot_real_space_quantities', &
           'green_density_fine', ierr)
      deallocate(green_tot_fine, stat = ierr)
      call utils_dealloc_check('dmft_plot_real_space_quantities', &
           'green_tot_fine', ierr)
      deallocate(Z_density_fine, stat = ierr)
      call utils_dealloc_check('dmft_plot_real_space_quantities', &
           'Z_density_fine', ierr)
      deallocate(green_density_fine_proj, stat = ierr)
      call utils_dealloc_check('dmft_plot_real_space_quantities', &
           'green_density_fine_proj', ierr)
   endif
   do is = 1, pub_num_spins
      call sparse_destroy( green_density_spam_proj(i) )
      call sparse_destroy( green_density_spam(i, PUB_1K))
      call sparse_destroy(     green_tot_spam(i)      )
      call sparse_destroy( sigma_density_spam(i)      )
      call sparse_destroy(     Z_density_spam(i)      )
   end do

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_plot_real_space_quantities'

 end subroutine

 subroutine dmft_compute_projected_overlap(greenf, hub_overlap, hub, &
      hub_proj_basis, overlap_cplx, overlap_hubbard, is, nmu_step, ikpoint, &
      ien, energy, inv_projectors)

   !---------------------------------------------------------------------------!
   ! Computes the overlap projected onto the Hubbard subspace                  !
   !---------------------------------------------------------------------------!

   use sparse,            only: sparse_create, sparse_product
   use utils,             only: utils_abort, utils_alloc_check, utils_dealloc_check
   use comms,             only: pub_my_proc_id, pub_on_root, pub_total_num_procs
   use parallel_strategy, only: par => pub_par
   use rundat,            only: pub_debug, pub_debug_on_root
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL

   implicit none

   ! Arguments
   type(DMFT_MATRIX),               intent(inout) :: greenf
   type(DMFT_OVERLAP),              intent(in   ) :: hub_overlap
   type(HUBBARD_MODEL),             intent(in   ) :: hub
   type(FUNC_BASIS),                intent(in   ) :: hub_proj_basis
   type(SPAM3),                     intent(in   ) :: overlap_cplx
   complex(kind = DP), allocatable, intent(inout) :: overlap_hubbard(:,:,:)
   integer,                         intent(in   ) :: is
   integer,                         intent(in   ) :: nmu_step
   integer,                         intent(in   ) :: ikpoint
   integer,                         intent(  out) :: ien
   complex(kind = DP),              intent(  out) :: energy
   logical,                         intent(in   ) :: inv_projectors

   ! Local variables
   integer                         :: i, k1, k2, j, l1, l2, hub_atom, kien
   type(SPAM3)                     :: w_S_v
   complex(kind = DP), allocatable :: matout(:,:)
   complex(kind = DP), allocatable :: fullmat(:,:)
   integer                         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_compute_projected_overlap'

   if(.not.inv_projectors) then
      write(stdout, *) 'Error in dmft_compute_projected_overlap: &
           &inv_projectors must be set to True to compute the projected &
           &overlap matrix'
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_compute_projected_overlap'
      return
   endif

   if(pub_output_detail /= VERBOSE) return

   if(pub_total_num_procs /= 1 .or. is /= 1 .or. nmu_step /= 1 .or. ikpoint /= &
        1) return

   call sparse_create(w_S_v, greenf%w_mat_v, iscmplx = .true.)

   energy = 0.0_DP
   ien = 1
   kien = 1  ! FOR MATMUL ON GPU BELOW - ENFORCE DOUBLE PREC !

   ! #ifndef GPU_SPEEDUP_MATMUL
   if(pub_debug_on_root) write(stdout, '(a)') &
        'DEBUG: compute < \psi^m | S | phi^\alpha > '
   call sparse_product(greenf%w_inv_mat, hub_overlap%inv_tmatc, overlap_cplx)
   ! #else
   ! if(pub_dmft_use_gpu_onlyme)then
   !    if(pub_on_root) write(stdout, *) 'Compute < \psi^m | S | phi^\alpha > '
   !    call dmft_matmul_in_dense_format(greenf%w_inv_mat, &
   !         hub_overlap%inv_tmatc, overlap_cplx)
   ! else
   !    call sparse_product(greenf%w_inv_mat, hub_overlap%inv_tmatc, overlap_cplx)
   ! endif
   ! #endif

   if(pub_debug_on_root) write(stdout, '(a)') &
         'DEBUG: compute < \psi^m | S | psi_m_prime > '
   ! #ifndef GPU_SPEEDUP_MATMUL
   call sparse_product(w_S_v, greenf%w_inv_mat, hub_overlap%inv_matc)
   ! #else
   ! if(pub_dmft_use_gpu_onlyme)then
   !    call dmft_matmul_in_dense_format(w_S_v, greenf%w_inv_mat, &
   !         hub_overlap%inv_matc)
   ! else
   !    call sparse_product(w_S_v, greenf%w_inv_mat, hub_overlap%inv_matc)
   ! endif
   ! #endif

   do hub_atom = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)
      i = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
      k1 = dmft_hub_atom_hub_projs(i, hub, hub_proj_basis)
      k2 = dmft_hub_atom_hub_projs(i, hub, hub_proj_basis)
      allocate(matout(k1, k2), stat = ierr)
      call utils_alloc_check('dmft_compute_projected_overlap', 'matout', ierr)
      if(i == 1)then
         ! ebl: we can't pass a slice of an allocatable array to a subroutine,
         ! hence the following gymnastics with fullmat
         call dmft_show_hubbard_atom(w_S_v, i, matout, hub, hub_proj_basis, &
              showdiag = .true., fullmat = fullmat)
         overlap_hubbard(:, :, is) = fullmat
         deallocate(fullmat, stat = ierr)
         call utils_dealloc_check('dmft_compute_projected_overlap', 'fullmat', &
              ierr)
      else
         call dmft_show_hubbard_atom(w_S_v, i, matout, hub, hub_proj_basis, &
              showdiag = .false.)
      endif
      if (pub_debug .and. maxval(abs(matout)) < 1.0E-8_DP) then
         call utils_abort('Error in DMFT: Hubbard projection matrix is empty')
      end if
      ! if(pub_on_root)then
      !    write(stdout, *) 'Re(S) matrix of ATOM [x] / total [y] : ', i, &
      !         par%nat_hub
      !    do l1 = 1, k1
      !       write(stdout, '(200f10.3)') (real(matout(l1, l2)), l2 = 1, k2)
      !    enddo
      ! endif
      deallocate(matout, stat = ierr)
      call utils_dealloc_check('dmft_compute_projected_overlap', 'matout', ierr)
   enddo

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_compute_projected_overlap'

 end subroutine


 subroutine dmft_project_green_function(greenf, hub_overlap, split_hub, &
      split_hubt, split_wgv, nosplit)

   !---------------------------------------------------------------------------!
   ! Projects the total Green function onto Hubbard subspaces                  !
   !---------------------------------------------------------------------------!

   use sparse, only: sparse_create, sparse_destroy, sparse_num_cols, &
        sparse_num_rows, sparse_product
   use utils,  only: utils_abort
   use comms,  only: pub_on_root
   use rundat, only: pub_debug, pub_debug_on_root

   implicit none

   ! Arguments
   type(DMFT_MATRIX),               intent(inout) :: greenf
   type(DMFT_OVERLAP),              intent(in   ) :: hub_overlap
   complex(kind = DP), allocatable, intent(inout) :: split_hub(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_hubt(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_wgv(:,:)
   logical, optional,               intent(in   ) :: nosplit

   ! Local variables
   type(SPAM3) :: w_mat
   type(SPAM3) :: test

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_project_green_function'

   if(.not.present(nosplit))then
      if(pub_dmft_split)then
         if(pub_on_root) write(stdout, *) 'BUILD PROJECTED GF'
         if (pub_debug) then
            if(size(split_hubt, 1) /= size(split_wgv, 1) .or. size(split_hub, &
                 2) /= size(split_wgv, 2))then
               call utils_abort('ERROR sizes dont match split hubt')
            endif
         end if
         split_wgv = MATMUL( MATMUL( split_hubt, greenf%backup), split_hub )
         if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
              &dmft_project_green_function'
         return
      endif
   endif

   if(pub_debug_on_root)write(stdout, *) 'Compute < \psi^m | \hat{G} | &
        &phi^\alpha > '
   ! #ifndef GPU_SPEEDUP_MATMUL
   call sparse_product(greenf%w_mat, hub_overlap%tmatc, greenf%mat)
   ! #else
   ! if(pub_dmft_use_gpu_onlyme)then
   !    if(pub_debug_on_root)write(stdout, *) 'go for full GPU'
   !    call dmft_matmul_in_dense_format(greenf%w_mat, hub_overlap%tmatc, &
   !         greenf%mat)
   ! else
   !    call sparse_product(greenf%w_mat, hub_overlap%tmatc, greenf%mat)
   !    if (pub_debug) then
   !       write(stdout, *) 'CHECKING SPARSE PRODUCT'
   !       call sparse_create(test, greenf%w_mat, iscmplx = .true.)
   !       call dmft_matmul_in_dense_format(test, hub_overlap%tmatc, greenf%mat)
   !       call internal_sparse_max_diff(test, greenf%w_mat)
   !       call sparse_destroy(test)
   !    end if
   ! endif
   ! #endif

   if(pub_debug_on_root)write(stdout, *) 'Compute < \psi^m | \hat{G} | &
        &psi_m_prime > '

   ! #ifndef GPU_SPEEDUP_MATMUL
   call sparse_product(greenf%w_mat_v, greenf%w_mat, hub_overlap%matc)
   ! #else
   ! if(pub_dmft_use_gpu_onlyme)then
   !    if(pub_debug_on_root)write(stdout, *) 'go for full GPU'
   !    call dmft_matmul_in_dense_format(greenf%w_mat_v, greenf%w_mat, &
   !         hub_overlap%matc)
   ! else
   !    call sparse_product(greenf%w_mat_v, greenf%w_mat, hub_overlap%matc)
   !    if (pub_debug) then
   !       write(stdout, *) 'CHECKING SPARSE PRODUCT'
   !       call sparse_create(test, greenf%w_mat_v, iscmplx = .true.)
   !       call dmft_matmul_in_dense_format(test, greenf%w_mat, hub_overlap%matc)
   !       call internal_sparse_max_diff(test, greenf%w_mat_v)
   !       call sparse_destroy(test)
   !    end if
   ! endif
   ! #endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_project_green_function'

 end subroutine

 subroutine dmft_project_spam(w_mat_v, mat, hub_overlap, split_hub, &
      split_hubt, split_wgv, greenfbackup, nosplit)

   !---------------------------------------------------------------------------!
   ! Projects a generic sparse matrix onto Hubbard subspaces                   !
   !---------------------------------------------------------------------------!
   ! NOTES                                                                     !
   !   This is a modified project_green_function instead of calling it as      !
   !   call project_green_function(Ok, invK, nosplit = .true.)                 !
   !   in dmft_load_k_ham_and_overlap                                          !
   !   This may have some residual behaviour that only ought to be applied in  !
   !   the green function case and not here e.g. the updating of split_wgv     !
   !---------------------------------------------------------------------------!

   use sparse, only: sparse_create, sparse_destroy, sparse_product
   use utils,  only: utils_abort
   use comms,  only: pub_on_root
   use rundat, only: pub_debug, pub_debug_on_root

   implicit none

   ! Arguments
   type(SPAM3),                     intent(  out) :: w_mat_v
   type(SPAM3),                     intent(in   ) :: mat
   type(DMFT_OVERLAP),              intent(in   ) :: hub_overlap
   complex(kind = DP), allocatable, intent(inout) :: split_hub(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_hubt(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_wgv(:,:)
   complex(kind = DP), allocatable, intent(inout) :: greenfbackup(:,:)
   logical, optional,               intent(in   ) :: nosplit

   ! Local variables
   type(SPAM3) :: w_mat
   type(SPAM3) :: test

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_project_spam'

   if(.not.present(nosplit))then
      if(pub_dmft_split)then
         if(pub_on_root) write(stdout, *) 'BUILD PROJECTED GF'
         if (pub_debug) then
            if(size(split_hubt, 1) /= size(split_wgv, 1) .or. size(split_hub, &
                 2) /= size(split_wgv, 2))then
               call utils_abort('ERROR sizes dont match split hubt')
            endif
         end if
         split_wgv = MATMUL( MATMUL( split_hubt, greenfbackup), split_hub )
         if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
              &dmft_project_spam'
         return
      endif
   endif

   if(pub_debug_on_root)write(stdout, *) 'Compute < \psi^m | \hat{M} | &
        &phi^\alpha > '
   ! #ifndef GPU_SPEEDUP_MATMUL
   call sparse_product(w_mat, hub_overlap%tmatc, mat)
   ! #else
   ! if(pub_dmft_use_gpu_onlyme)then
   !    if(pub_debug_on_root)write(stdout, *) 'go for full GPU'
   !    call dmft_matmul_in_dense_format(w_mat, hub_overlap%tmatc, mat)
   ! else
   !    call sparse_product(w_mat, hub_overlap%tmatc, mat)
   !    if (pub_debug) then
   !       write(stdout, *) 'CHECKING SPARSE PRODUCT'
   !       call sparse_create(test, w_mat, iscmplx = .true.)
   !       call dmft_matmul_in_dense_format(test, hub_overlap%tmatc, mat)
   !       call internal_sparse_max_diff(test, w_mat)
   !       call sparse_destroy(test)
   !    end if
   ! endif
   ! #endif

   if(pub_debug_on_root)write(stdout, *) 'Compute < \psi^m | \hat{G} | &
        &psi_m_prime > '

   ! #ifndef GPU_SPEEDUP_MATMUL
   call sparse_product(w_mat_v, w_mat, hub_overlap%matc)
   ! #else
   ! if(pub_dmft_use_gpu_onlyme)then
   !    if(pub_debug_on_root)write(stdout, *) 'go for full GPU'
   !    call dmft_matmul_in_dense_format(w_mat_v, w_mat, hub_overlap%matc)
   ! else
   !    call sparse_product(w_mat_v, w_mat, hub_overlap%matc)
   !    if (pub_debug) then
   !       write(stdout, *) 'CHECKING SPARSE PRODUCT'
   !       call sparse_create(test, w_mat_v, iscmplx = .true.)
   !       call dmft_matmul_in_dense_format(test, w_mat, hub_overlap%matc)
   !       call internal_sparse_max_diff(test, w_mat_v)
   !       call sparse_destroy(test)
   !    end if
   ! endif
   ! #endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_project_spam'

 end subroutine

 subroutine dmft_project_hamiltonian(greenf, hub_overlap, inv_projectors)

   !---------------------------------------------------------------------------!
   ! Projects the Hamiltonian onto Hubbard subspaces                           !
   !---------------------------------------------------------------------------!

   use utils,  only: utils_assert
   use sparse, only: sparse_product
   use rundat, only: pub_debug, pub_debug_on_root

   implicit none

   ! Arguments
   type(DMFT_MATRIX),  intent(inout) :: greenf
   type(DMFT_OVERLAP), intent(inout) :: hub_overlap
   logical,            intent(in   ) :: inv_projectors

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &project_hamiltonian'

   call utils_assert(inv_projectors, "Error in dmft_project_hamiltonian: &
        &inv_projectors must be set to true")

   ! #ifndef GPU_SPEEDUP_MATMUL
   if(pub_debug_on_root) write(stdout, *) 'Compute < \psi^m | E-\hat{H} | &
        &phi^\alpha > '
   call sparse_product(greenf%w_inv_mat, hub_overlap%inv_tmatc, greenf%inv_mat)
   ! #else
   ! if(pub_dmft_use_gpu_onlyme) then
   !    if(pub_debug_on_root)write(stdout, *) 'Compute < \psi^m | E-\hat{H} | &
   !         &phi^\alpha > '
   !    call dmft_matmul_in_dense_format(greenf%w_inv_mat, &
   !         hub_overlap%inv_tmatc, greenf%inv_mat)
   ! else
   !    call sparse_product(greenf%w_inv_mat, hub_overlap%inv_tmatc, &
   !         greenf%inv_mat)
   ! endif
   ! #endif

   if(pub_debug_on_root) write(stdout, *) 'Compute < \psi^m | E-\hat{H} | &
        &psi_m_prime > '
   ! #ifndef GPU_SPEEDUP_MATMUL
   call sparse_product(greenf%w_inv_mat_v, greenf%w_inv_mat, &
        hub_overlap%inv_matc)
   ! #else
   ! if(pub_dmft_use_gpu_onlyme) then
   !    call dmft_matmul_in_dense_format(greenf%w_inv_mat_v, greenf%w_inv_mat, &
   !         hub_overlap%inv_matc)
   ! else
   !    call sparse_product(greenf%w_inv_mat_v, greenf%w_inv_mat, &
   !         hub_overlap%inv_matc)
   ! endif
   ! #endif
   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &project_hamiltonian'

 end subroutine

 subroutine dmft_init_matrices(ham, greenf, self_energy, self_infinity, overlap, &
      hub_overlap, ngwf_basis, kpoint_data, sc_data, overlap_cplx, matrixAsig, &
      matrixSigStat, matrixDensKernel, nkpoints, inv_overlap, &
      sigma_density_spam, green_density_spam, green_density_spam_proj, &
      green_tot_spam, Z_density_spam, cluster_scheme, inv_projectors)

   !---------------------------------------------------------------------------!
   ! Initialisation of most sparse/dense DMFT matrices required throughout the !
   ! entire DMFT loop                                                          !
   !---------------------------------------------------------------------------!
   ! HISTORY                                                                   !
   !   Mar 2017 ebl27 split dmft_init_sparse_and_dimers into                   !
   !                  dmft_init_matrices and dmft_build_my_dimers              !
   !---------------------------------------------------------------------------!

   use comms,          only: pub_on_root
   use dense,          only: dense_create
   use function_basis, only: FUNC_BASIS
   use rundat,         only: pub_debug, pub_debug_on_root, pub_num_spins
   use sparse,         only: sparse_copy, sparse_create, sparse_product
   use utils,          only: utils_alloc_check

   implicit none

   ! Arguments
   type(SPAM3),                intent(inout) :: ham(pub_num_spins)
   type(DMFT_MATRIX),          intent(inout) :: greenf
   type(DMFT_MATRIX),          intent(inout) :: self_energy
   type(SPAM3),                intent(inout) :: self_infinity(pub_num_spins)
   type(SPAM3),                intent(in   ) :: overlap
   type(DMFT_OVERLAP),         intent(inout) :: hub_overlap
   type(FUNC_BASIS),           intent(in   ) :: ngwf_basis
   type(DMFT_KPOINT_MATRICES), intent(inout) :: kpoint_data
   type(DMFT_SC_MATRICES),     intent(inout) :: sc_data
   type(SPAM3),                intent(inout) :: overlap_cplx
   type(DEM), allocatable,     intent(inout) :: matrixAsig(:)
   type(DEM), allocatable,     intent(inout) :: matrixSigStat(:)
   type(SPAM3),                intent(inout) :: matrixDensKernel
   integer,                    intent(in   ) :: nkpoints
   type(SPAM3),                intent(in   ) :: inv_overlap
   type(SPAM3),                intent(inout) :: sigma_density_spam(2)
   type(SPAM3),                intent(inout) :: green_density_spam(2, 1)
   type(SPAM3),                intent(inout) :: green_density_spam_proj(2)
   type(SPAM3),                intent(inout) :: green_tot_spam(2)
   type(SPAM3),                intent(inout) :: Z_density_spam(2)
   integer,                    intent(in   ) :: cluster_scheme
   logical,                    intent(in   ) :: inv_projectors

   ! Local variables
   integer :: i, ii, ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_init_sparse'

   greenf%mat%structure = 'K'
   call sparse_create(greenf%mat, iscmplx = .true.)
   call sparse_create(greenf%inv_mat, ham(1), iscmplx = .true.)
   call sparse_create(greenf%inv_backup_spam, ham(1), iscmplx = .true.)
   call sparse_create(greenf%backup_spam, greenf%mat, iscmplx = .true.)
   call sparse_create(overlap_cplx, overlap, iscmplx = .true.)
   call sparse_create(hub_overlap%matc, hub_overlap%mat, iscmplx = .true.)
   call sparse_create(hub_overlap%tmatc, hub_overlap%tmat, iscmplx = .true.)
   call sparse_create(greenf%w_mat, hub_overlap%tmatc, greenf%mat, iscmplx = &
        .true.) ! WGr

   if(pub_debug_on_root) then
      write(stdout,"(a/)") repeat("-",80)
   end if


   select case(cluster_scheme)
   case(1)
      greenf%w_mat_v%structure = 'G' ! block-diagonal
   case(2)
      greenf%w_mat_v%structure = 'FULL'
   case(3)
      call sparse_create(greenf%w_mat_v, greenf%w_mat, hub_overlap%matc, &
           iscmplx = .true.)
   end select
   if(cluster_scheme /= 3) then
      call sparse_create(greenf%w_mat_v, iscmplx = .true.) ! WGrV- > 'G'
   endif

   if(nkpoints > 1)then
      if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, *) &
           'building arrays for K points'
      kpoint_data%unrot_overlap%structure = 'FULL'
      call sparse_create(       kpoint_data%unrot_overlap, iscmplx = .true.)
      call sparse_create( kpoint_data%unrot_inv_overlap, &
           kpoint_data%unrot_overlap, iscmplx = .true.)
      call sparse_create( kpoint_data%hub_overlap, kpoint_data%unrot_overlap, &
           iscmplx = .true.)
      call sparse_create(kpoint_data%self_energy, greenf%mat, iscmplx = .true.)
   endif

   if(inv_projectors)then
      if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, *) &
           'building arrays for inv. projectors'
      call sparse_create(hub_overlap%inv_mat, inv_overlap, hub_overlap%mat)
      call sparse_create(hub_overlap%inv_matc, hub_overlap%inv_mat, iscmplx = &
           .true.)
      call sparse_create(hub_overlap%inv_tmat, hub_overlap%tmat, inv_overlap)
      call sparse_create(hub_overlap%inv_tmatc, hub_overlap%inv_tmat, iscmplx &
           = .true.)
      call sparse_create(greenf%w_inv_mat, hub_overlap%inv_tmatc, &
           greenf%inv_mat, iscmplx = .true.)
      select case(cluster_scheme)
      case(1)
         greenf%w_inv_mat_v%structure = 'G' ! Blank to block-diagonal 'G'
         ! structure
      case(2)
         greenf%w_inv_mat_v%structure = 'FULL'
      case(3)
         call sparse_create(greenf%w_inv_mat_v, greenf%w_mat, &
              hub_overlap%matc, iscmplx = .true.)
      end select
      if(.not.cluster_scheme == 3) call sparse_create(greenf%w_inv_mat_v, &
           iscmplx = .true.)
   endif

   select case(cluster_scheme)
   case(1)
      self_energy%w_mat_v%structure = 'G' ! Blank to block-diagonal 'G'
      ! structure
   case(2)
      self_energy%w_mat_v%structure = 'FULL'
   case(3)
      call sparse_create(self_energy%w_mat_v, greenf%w_mat, hub_overlap%matc, &
           iscmplx = .true.)
   end select
   if(.not.cluster_scheme == 3) call sparse_create(self_energy%w_mat_v, &
        iscmplx = .true.)

   call sparse_create(greenf%mat_v, hub_overlap%matc, greenf%w_mat_v, iscmplx &
        = .true.)
   call sparse_create(self_energy%mat_v, hub_overlap%matc, &
        self_energy%w_mat_v, iscmplx = .true.)
   call sparse_create(self_energy%mat, self_energy%mat_v, hub_overlap%tmatc, &
        iscmplx = .true.)

   do i = 1, pub_num_spins
      call sparse_create(self_infinity(i), ham(i), iscmplx = (nkpoints > 1))
   enddo

   call sparse_copy(overlap_cplx, overlap)
   call sparse_copy(hub_overlap%matc, hub_overlap%mat)
   call sparse_copy(hub_overlap%tmatc, hub_overlap%tmat)

   if(inv_projectors)then
      ! #ifndef GPU_SPEEDUP_MATMUL_R
      call sparse_product(hub_overlap%inv_mat, inv_overlap, hub_overlap%mat)
      ! #else
      !       if(pub_dmft_use_gpu_onlyme)then
      ! call dmft_matmul_in_dense_format_r(hub_overlap%inv_mat, inv_overlap,
      ! hub_overlap%mat)
      !       else
      ! call sparse_product(hub_overlap%inv_mat, inv_overlap, hub_overlap%mat)
      !       endif
      ! #endif
      call sparse_copy(hub_overlap%inv_matc, hub_overlap%inv_mat)

      ! #ifndef GPU_SPEEDUP_MATMUL_R
      call sparse_product(hub_overlap%inv_tmat, hub_overlap%tmat, inv_overlap)
      ! #else
      !       if(pub_dmft_use_gpu_onlyme)then
      ! if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, *)
      ! 'initialization sparse step [9.1]'
      ! call dmft_matmul_in_dense_format_r(hub_overlap%inv_tmat,
      ! hub_overlap%tmat, inv_overlap)
      !       else
      ! if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, *)
      ! 'initialization sparse step [9.2]'
      ! call sparse_product(hub_overlap%inv_tmat, hub_overlap%tmat,
      ! inv_overlap)
      !       endif
      ! #endif

      call sparse_copy(hub_overlap%inv_tmatc, hub_overlap%inv_tmat)
   endif

   do i = 1, pub_num_spins
      call sparse_create( green_density_spam_proj(i), greenf%mat, iscmplx = &
           .false.)
      call sparse_create( green_density_spam(i, PUB_1K), greenf%mat, iscmplx = &
           .false.)
      call sparse_create( green_tot_spam(i), greenf%mat, iscmplx = .false.)
      call sparse_create( sigma_density_spam(i), self_energy%mat, iscmplx = &
           .false.)
      call sparse_create( Z_density_spam(i), sigma_density_spam(i), iscmplx = &
           .false.)
   enddo

   ii = ngwf_basis%num

   allocate(sc_data%eigenvals(pub_num_spins, ii), stat = ierr)
   call utils_alloc_check('dmft_init_matrices', 'sc_data%eigenvals', ierr)
   if (pub_debug) then
      allocate(sc_data%eigenvals_backup(pub_num_spins, ii), stat = ierr)
      call utils_alloc_check('dmft_init_matrices', 'sc_data%eigenvals_backup', &
           ierr)
   end if

   allocate(sc_data%eigenvecs(pub_num_spins), stat = ierr)
   call utils_alloc_check('dmft_init_matrices', 'sc_data%eigenvecs', ierr)
   allocate(sc_data%eigenvecs_t(pub_num_spins), stat = ierr)
   call utils_alloc_check('dmft_init_matrices', 'sc_data%eigenvecs_t', ierr)
   allocate(sc_data%ham(pub_num_spins), stat = ierr)
   call utils_alloc_check('dmft_init_matrices', 'sc_data%ham', ierr)
   allocate(matrixSigStat(pub_num_spins), stat = ierr)
   call utils_alloc_check('dmft_init_matrices', 'matrixSigStat', ierr)
   allocate(matrixAsig(pub_num_spins), stat = ierr)
   call utils_alloc_check('dmft_init_matrices', 'matrixAsig', ierr)

   do i = 1, pub_num_spins
      call dense_create(sc_data%eigenvecs(i), ngwf_basis%num, ngwf_basis%num, &
           nkpoints > 1)
      call dense_create(sc_data%eigenvecs_t(i), ngwf_basis%num, &
           ngwf_basis%num, nkpoints > 1)
      call sparse_create(sc_data%ham(i), ham(1), iscmplx = nkpoints > 1)
   end do

   call dense_create(sc_data%tail, ngwf_basis%num, ngwf_basis%num, nkpoints > 1)
   call sparse_create(sc_data%overlap, overlap, iscmplx = nkpoints > 1)

   ! if (nkpoints == 1) then
   ! allocate(sc_data%eigenvecs(pub_num_spins, ii, ii), sc_data%tail(ii, ii),
   ! sc_data%overlap(ii, ii), sc_data%ham(pub_num_spins, ii, ii))
   ! else
   ! allocate(sc_data%eigenvecs_c(pub_num_spins, ii, ii), sc_data%tail_c(ii,
   ! ii), sc_data%overlap_c(ii, ii), sc_data%ham_c(pub_num_spins, ii, ii))
   ! endif

   call sparse_create(matrixDensKernel, ham(1), iscmplx = nkpoints > 1)
   ! if (nkpoints == 1) then
   !    allocate(matrixDensKernel(ii, ii))
   ! else
   !    allocate(matrixDensKernelc(ii, ii))
   ! endif
   do i = 1, pub_num_spins
      call dense_create(matrixAsig(i), ngwf_basis%num, ngwf_basis%num, .true.)
      call dense_create(matrixSigStat(i), ngwf_basis%num, ngwf_basis%num, &
           .true.)
   end do

   allocate(greenf%backup(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_init_matrices', 'greenf%backup', ierr)
   allocate(greenf%inv_backup(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_init_matrices', 'greenf%inv_backup', ierr)
   allocate(self_energy%backup(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_init_matrices', 'self_energy%backup', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving dmft_init_sparse'

 end subroutine

 subroutine dmft_init_routine_variables(ham, overlap, h_atoms, hub, &
      hub_overlap, hub_proj_basis, overlap_hubbard, num_spins, en_start, &
      en_step, en_step_over_pi, energyzero, nspecial_frequ, cluster, &
      cluster_scheme, force_reopen_file, inv_overlap, &
      reopen_sigma_file_each_step, restrict_window, same_self_for_all, &
      use_gpu_eigenvectors)

   !---------------------------------------------------------------------------!
   ! Initialise various options for the DMFT calculation                       !
   !---------------------------------------------------------------------------!

   use comms, only: comms_reduce, pub_my_proc_id, pub_on_root, &
        pub_total_num_procs
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL, hubbard_atoms_occupancy
   use parallel_strategy, only: par => pub_par
   use rundat, only: pub_debug_on_root, pub_dmft_complex_freq, &
        pub_dmft_dos_max, pub_dmft_dos_min, pub_dmft_emax, pub_dmft_emin, &
        pub_dmft_nkpoints, pub_dmft_points, pub_dmft_read, pub_dmft_temp, &
        pub_num_spins, pub_rootname
   use sparse,            only: sparse_copy
   use utils,             only: utils_abort, utils_assert, utils_banner, &
        utils_alloc_check, utils_dealloc_check

   implicit none

   ! Arguments
   type(SPAM3),                     intent(inout) :: ham(pub_num_spins)
   type(SPAM3),                     intent(inout) :: overlap
   type(ARRAY_OF_MATRICES),         intent(inout) :: h_atoms(par%nat_hub)
   type(HUBBARD_MODEL),             intent(in   ) :: hub
   type(DMFT_OVERLAP),              intent(inout) :: hub_overlap
   type(FUNC_BASIS),                intent(in   ) :: hub_proj_basis
   complex(kind = DP), allocatable, intent(inout) :: overlap_hubbard(:,:,:)
   type(SPAM3),                     intent(inout) :: inv_overlap
   integer,                         intent(in   ) :: num_spins
   real(kind = DP),                 intent(  out) :: en_start
   real(kind = DP),                 intent(  out) :: en_step
   real(kind = DP),                 intent(  out) :: en_step_over_pi
   integer,                         intent(  out) :: energyzero
   integer,                         intent(in   ) :: nspecial_frequ
   logical,                         intent(in   ) :: cluster
   integer,                         intent(  out) :: cluster_scheme
   logical,                         intent(  out) :: force_reopen_file
   logical,                         intent(  out) :: reopen_sigma_file_each_step
   logical,                         intent(  out) :: restrict_window
   logical,                         intent(  out) :: same_self_for_all
   logical,                         intent(in   ) :: use_gpu_eigenvectors

   ! Local variables
   integer         :: hub_atom
   integer         :: hub_nproj
   logical         :: ham1_exists
   logical         :: ham2_exists
   integer         :: jij
   integer         :: u(1), ii
   integer         :: channels
   integer         :: hat
   integer         :: en_range
   integer         :: length
   integer         :: stat
   real(kind = DP) :: temp_times_pi
   integer         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_init_routine_variables'

   same_self_for_all = .false.
   force_reopen_file = .false.

   hub_nproj = sum( (/( dmft_hub_atom_hub_projs(jij, hub, hub_proj_basis), jij &
        = 1, par%nat_hub )/) )

   if(allocated(overlap_hubbard)) then
      deallocate(overlap_hubbard, stat = ierr)
      call utils_dealloc_check('dmft_init_routine_variables', &
           'overlap_hubbard', ierr)
   end if
   allocate(overlap_hubbard(hub_nproj, hub_nproj, 2), stat = ierr)
   call utils_alloc_check('dmft_init_routine_variables', 'overlap_hubbard', &
        ierr)

   cluster_scheme = 1
   if(cluster) cluster_scheme = 2

   ! ebl: removing this command_agrument_count procedure
   ! CWTODO: check this, I am unsure of the way you parse arguments to ONETEP

   ! For the moment, leave pub_dmft_mpi_size and _rank always to 1
   pub_dmft_mpi_size = 1
   pub_dmft_mpi_rank = 1

   ! num = COMMAND_ARGUMENT_COUNT()
   ! if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, *)
   ! 'NUMBER OF ARGUMENTS IN ONETEP : ', num
   ! do i = 1, num
   !    call GET_COMMAND_ARGUMENT(i, value(i), length, stat)
   ! if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, *)
   ! 'ARGUMENT : ', i, TRIM(ADJUSTL(value(i)))
   ! enddo
   ! if(num >= 2)then
   !    pub_dmft_mpi_size = StrInt2(value(3))
   !    pub_dmft_mpi_rank = StrInt2(value(2))
   ! if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, *)
   ! 'RUNNING ONETEP AS MULTI_PROCS : '
   ! if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, *) 'MY
   ! RANK / SIZE : ', pub_dmft_mpi_rank, pub_dmft_mpi_size
   ! endif
   ! if(num >= 4)then
   !    pub_dmft_points = StrInt2(value(4))
   !    sign_temp       = StrInt2(value(5))
   !    pub_dmft_emin   = StringToReal(value(6))
   !    pub_dmft_emax   = StringToReal(value(7))
   !    if(sign_temp /= 0) pub_dmft_temp = -abs(pub_dmft_temp)

   ! ebl: when performing DOS/optics calculation, work with real frequencies
   ! This is also denoted within the code as a negative temperature
   if (.not. pub_dmft_complex_freq)  pub_dmft_temp = -abs(pub_dmft_temp)

   if (pub_on_root) then
      write(stdout,'(a)') utils_banner('-', ' Calculation parameters ')
      write(stdout, "(a)") " WARNING: not reading the argument given to ONETEP"
      write(stdout, "(a,i9)")      ' Number of frequency points    : ', pub_dmft_points
      write(stdout, "(a,e13.2e2)") ' Temperature (Ha)              : ', pub_dmft_temp
      write(stdout, "(a,f9.2)")    ' Minimum energy (Ha)           : ', pub_dmft_emin
      write(stdout, "(a,f9.2)")    ' Maximum energy (Ha)           : ', pub_dmft_emax
      write(stdout, "(a/)") repeat("-", 80)
   end if

   restrict_window = ((pub_dmft_dos_min > pub_dmft_emin) .or. &
        (pub_dmft_dos_max < pub_dmft_emax)) .and. (.not.pub_dmft_complex_freq)
   reopen_sigma_file_each_step = pub_dmft_mpi_size > 1 .or. restrict_window &
        .or. pub_dmft_split

   pub_dmft_use_gpu = (pub_total_num_procs <= pub_dmft_gpu_num) .and. &
        (pub_dmft_mpi_size <= pub_dmft_gpu_num)
   if(pub_dmft_mpi_size == 1 .and. pub_total_num_procs > 1 .and. &
        pub_dmft_gpu_num == 1 .and. .not.pub_dmft_splitk)then
      write(stdout, *) 'WARNING switching off split due to the presence of a &
           &GPU'
      pub_dmft_use_gpu = .true.
      pub_dmft_split = .false.
   endif
   pub_dmft_use_gpu_onlyme = pub_dmft_use_gpu
   if(pub_dmft_mpi_size /= 1 .and. pub_total_num_procs == 1 .and. &
        pub_dmft_gpu_num >= 1 .and. .not.pub_dmft_splitk)then
      pub_dmft_use_gpu_onlyme = pub_dmft_mpi_rank <= pub_dmft_gpu_num
   endif
   pub_dmft_use_gpu_partially = pub_dmft_mpi_size /= 1 .and. &
        pub_total_num_procs == 1 .and. pub_dmft_gpu_num >= 1 .and. &
        pub_dmft_gpu_num < pub_dmft_mpi_size .and. .not.pub_dmft_splitk

   if(pub_dmft_splitk) then
      if(pub_dmft_use_gpu_onlyme) write(stdout, *) 'WARNING turning off GPU, &
           &because of splitk'
      pub_dmft_use_gpu = .false.
      pub_dmft_use_gpu_partially = .false.
      pub_dmft_use_gpu_onlyme = .false.
   endif

   if (pub_debug_on_root) then
      if(pub_output_detail >= VERBOSE) then
         write(stdout, "(a)")  "============= DMFT CPU/GPU SCHEME ============="
         write(stdout,"(a,i4)")' Number of MPI procs :', pub_total_num_procs
         write(stdout,"(a,i4)")' MPI rank            :', pub_dmft_mpi_rank
         write(stdout,"(a,i4)")' MPI proc ID         :', pub_my_proc_id
         write(stdout,"(a,l4)")' Using GPUs          :', pub_dmft_use_gpu
         write(stdout,"(a,i4)")' Number of GPUs      :', pub_dmft_gpu_num
         write(stdout,"(a,l4)")' This proc is a GPU  :', &
              pub_dmft_use_gpu_onlyme
         write(stdout,"(a,l4)")' Partially using GPUs ', &
              pub_dmft_use_gpu_partially
         write(stdout, "(a)")  "==============================================="
      end if
   end if

   ! excluding since init_gpu_device and choose_gpu are not defined
   ! #ifdef GPU_SPEEDUP
   ! if(pub_dmft_use_gpu_onlyme)then
   !    call init_gpu_device
   !    if(pub_total_num_procs == 1)then
   !       call choose_gpu(pub_dmft_mpi_rank-1)
   !    else
   !       if(pub_my_proc_id < pub_dmft_gpu_num)then
   !          call choose_gpu(pub_my_proc_id)
   !       endif
   !    endif
   !    call cublas_init()
   ! endif
   ! #endif

   do hat = 1, par%nat_hub
      if(allocated(h_atoms(hat)%occupancy)) then
         deallocate(h_atoms(hat)%occupancy, stat = ierr)
         call utils_dealloc_check('dmft_init_routine_variables', 'h_atoms', &
              ierr)
      end if
      channels = dmft_hub_atom_hub_projs(hat, hub, hub_proj_basis)
      allocate(h_atoms(hat)%occupancy(1:channels, 1:channels, 2), stat = ierr)
      call utils_alloc_check('dmft_init_routine_variables', 'h_atoms', ierr)
      h_atoms(hat)%occupancy = 0.0_DP
   enddo

   if (pub_debug_on_root) write(stdout, *) 'DEBUG: building occupancy matrices'

   do hub_atom = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)
      hat      = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
      channels = dmft_hub_atom_hub_projs(hat, hub, hub_proj_basis)
      if(par%proc_of_atom(hub%orig(hat)) /= pub_my_proc_id)then
         write(stdout, *) 'ERROR in getting h_atom, proc id do not match'
         write(stdout, *) ' hat                             : ', hat
         write(stdout, *) ' orig atom                       : ', hub%orig(hat)
         write(stdout, *) ' hub_atom                        : ', hub_atom
         write(stdout, *) ' par%proc_of_atom(hub%orig(hat)) : ', &
              par%proc_of_atom(hub%orig(hat))
         write(stdout, *) ' pub proc id                     : ', pub_my_proc_id
         call utils_abort("Error in dmft_init_routine_variables: proc id for &
              &hubbard atoms does not match")
      endif
      h_atoms(hat)%occupancy = 0.0_DP
      h_atoms(hat)%occupancy(1:channels, 1:channels, 1) = &
           hubbard_atoms_occupancy(hub_proj_basis, hub, hat, 1)
      if(num_spins == 2)then
         h_atoms(hat)%occupancy(1:channels, 1:channels, 2) = &
              hubbard_atoms_occupancy(hub_proj_basis, hub, hat, 2)
      endif
   enddo
   do hat = 1, par%nat_hub
      call comms_reduce('SUM', h_atoms(hat)%occupancy)
   enddo

   ! ebl: checking if .ham files are present

   ! OLD APPROACH - kept to allow comparison to ONETEP 3.1 implementation
   ! call check_sparse_in_file(ham(1), trim(pub_rootname)//'.ham1')
   ! if(num_spins == 2)then
   !    INQUIRE(file = trim(pub_rootname)//'.ham2', EXIST = ham_is_present)
   !    if(ham_is_present)then
   !       call check_sparse_in_file(ham(2), trim(pub_rootname)//'.ham2')
   !    else
   !       call sparse_copy(ham(2), ham(1))
   !    endif
   ! endif

   ! ebl: reading/writing of .ham files
   inquire(file = trim(pub_rootname) // '.ham1', exist = ham1_exists)
   inquire(file = trim(pub_rootname) // '.ham2', exist = ham2_exists)
   ! ebl: If .ham1 exists, read it. If not, write ham(1) to .ham1
   call check_sparse_in_file(ham(1), trim(pub_rootname) // '.ham1')

   ! ebl: in the case of spin-polarised...
   if (num_spins == 2) then
      if ((ham1_exists .and. (.not. ham2_exists)) .and. pub_dmft_read) then
         ! Copy ham(1) over to other spin channel
         write(stdout, *) 'WARNING: ' // trim(pub_rootname) // '.ham2 file is &
              &not present; copying ' // trim(pub_rootname)//'.ham1'
         call sparse_copy(ham(2), ham(1))
      end if
      ! If .ham2 exists, read it. If not, write ham(2) to .ham2
      call check_sparse_in_file(ham(2), trim(pub_rootname) // '.ham2')
   endif

   call check_sparse_in_file(overlap, trim(pub_rootname) // '.overlap')
   call check_sparse_in_file(inv_overlap, trim(pub_rootname) // '.inv_overlap')
   call check_sparse_in_file(hub_overlap%mat, trim(pub_rootname) // &
        '.hub_overlap')
   call check_sparse_in_file(hub_overlap%tmat, trim(pub_rootname) // &
        '.hub_overlap_t')

   if(pub_dmft_invert_overlap)then
      call dmft_invert_sparse_real(overlap, inv_overlap, 1)
   endif

   if(pub_dmft_norm_proj /= 0) then
      call utils_assert(abs(pub_dmft_nkpoints) < 1, "Error in &
           &dmft_init_routine_variables: renormalised Green function and &
           &k-points are not currently compatible")

      if(pub_dmft_norm_proj > 0)then
         call dmft_normalise_projectors(hub_overlap%mat, 2, overlap, &
              use_gpu_eigenvectors, mat_t = hub_overlap%tmat)
      else
         call dmft_normalise_projectors(hub_overlap%mat, 2, overlap, &
              use_gpu_eigenvectors)
         call dmft_normalise_projectors(hub_overlap%tmat, 1, overlap, &
              use_gpu_eigenvectors)
      endif
   endif

   temp_times_pi = ABS(pub_dmft_temp)*PI
   if (pub_dmft_temp > 0.0_DP) then
      en_range = REAL(2*pub_dmft_points, kind = DP)*temp_times_pi
      en_start = temp_times_pi
      en_step  = 2.0_DP*temp_times_pi
   else
      en_range = pub_dmft_emax-pub_dmft_emin
      en_start = pub_dmft_emin
      en_step  = en_range/REAL(pub_dmft_points-1, kind = DP)
   endif

   en_step_over_pi = -1.0_DP * en_step / PI

   ! Define energyzero
   if (pub_dmft_temp > 0.0_DP) then
      energyzero = 1
   else
      u = minloc(abs( (/( en_start + real(ii-1, kind = DP)*en_step, ii = 1, &
           pub_dmft_points-nspecial_frequ )/) ))
      energyzero = u(1)
   endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_init_routine_variables'

 end subroutine

 subroutine dmft_init_scratch_dir(scratchdir)

   !---------------------------------------------------------------------------!
   ! Locates a scratch directory if desired                                    !
   !---------------------------------------------------------------------------!

   use comms,             only: comms_barrier, pub_on_root
   use file_handling,     only: file_handling_copy_file
   use utils,             only: utils_abort, utils_unit, utils_int_to_str
   use parallel_strategy, only: par => pub_par
   use rundat, only: pub_debug_on_root, pub_num_spins

   implicit none

   ! Arguments
   character(2000), intent(  out) :: scratchdir

   ! Local variables
   logical         :: checkscratch
   real(kind = DP) :: lock
   integer         :: hat, sp, is
   integer         :: free_unit
   character(10)   :: filename
   character(100)  :: filename2

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_init_scratch_dir'

   if(pub_dmft_local_scratch)then
      lock = 1.0_DP
      inquire(file = 'scratchdir', exist = checkscratch)
      if(.not.checkscratch)then
         call utils_abort('Error in dmft_init_scratch_dir: file "scratchdir" &
              &is not present')
      endif
      free_unit = utils_unit()
      open(unit = free_unit, file = 'scratchdir')
      read(free_unit, '(a)') scratchdir
      close(free_unit)
      if(pub_dmft_mpi_rank == 1 .and. pub_on_root)then
         write(stdout, *) 'The scratch directory is ', trim(adjustl(scratchdir))
         do hat = 1, par%nat_hub
            do sp = 1, par%num_hub_species
               do is = 1, pub_num_spins
                  filename = 'sigma_output' // trim(adjustl(utils_int_to_str( &
                       hat))) // "_" // trim(adjustl( utils_int_to_str(sp))) &
                       // "_" // trim(adjustl( utils_int_to_str(is)))
                  filename2 = "./" // trim(adjustl(scratchdir)) // "/" // &
                       filename
                  call file_handling_copy_file(filename, filename2, formatted &
                       = .false.)
               end do
            end do
         end do
      endif
      call utils_abort("Reached point in code where nfs sync used to appear; &
           &the code needs to be updated to use these settings.")
      ! if(pub_on_root) call sync_via_nfs(lock, 'copy_sigma_output_scratch_dir', &
      !      longdelay = .true.)
      ! call comms_barrier
   endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_init_scratch_dir'

 end subroutine

 subroutine dmft_distribute_self_energy(site_self_energy_buffer, self_energy, &
      connection_table, hub, hub_proj_basis, is, rot_vec_angles, h_atoms, &
      cluster)

   !---------------------------------------------------------------------------!
   ! Distribute self energy to atomic components                               !
   !---------------------------------------------------------------------------!

   use utils,             only: utils_abort
   use sparse,            only: sparse_put_element
   use comms,             only: pub_my_proc_id
   use rundat,            only: pub_debug_on_root, pub_dmft_rotate_green
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par

   implicit none

   ! Arguments
   complex(kind = DP),      intent(inout) :: site_self_energy_buffer(:,:,:,:)
   type(DMFT_MATRIX),       intent(inout) :: self_energy
   logical, allocatable,    intent(in   ) :: connection_table(:,:)
   type(HUBBARD_MODEL),     intent(in   ) :: hub
   type(FUNC_BASIS),        intent(in   ) :: hub_proj_basis
   integer,                 intent(in   ) :: is
   real(kind = DP),         intent(in   ) :: rot_vec_angles(:,:,:)
   type(ARRAY_OF_MATRICES), intent(in   ) :: h_atoms(par%nat_hub)
   logical,                 intent(in   ) :: cluster

   ! Local variables
   complex(kind = DP)   :: self_energy_element
   type(DMFT_LOOP_INFO) :: loop
   integer              :: hub_atom
   integer              :: hub_atomb
   integer              :: row_counter
   integer              :: row_proj
   integer              :: col_counter
   integer              :: col_proj

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_distribute_self_energy'

   do hub_atom = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)
      do hub_atomb = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

         call dmft_set_loop_variable(loop, hub_atom, hub_atomb, hub, &
              hub_proj_basis)

         if(.not.cluster .and. hub_atom /= hub_atomb) cycle
         if(cluster)then
            if(.not. connection_table(loop%hat, loop%hatb) ) cycle
         endif

         ! Rotate self func
         if (pub_dmft_rotate_green) then

            if (loop%channels /= loop%channelsb) then
               write(stdout, *) 'hub_atom, hub_atomb : ', hub_atom, hub_atomb
               write(stdout, *) 'hat, hatb            : ', loop%hat, loop%hatb
               write(stdout, *) 'table connection : ', &
                    connection_table(loop%hat, loop%hatb)
               write(stdout, *) 'channels, channelsb : ', loop%channels, &
                    loop%channelsb
               call utils_abort('Error in dmft_distribute_self_energy: &
                    &channels, channelsb are different, incompatible with &
                    &local rotation')
            endif

            call dmft_rotate_to_cart_basis( site_self_energy_buffer(hub_atom, &
                 hub_atomb, 1:loop%channels, 1:loop%channelsb), &
                 h_atoms(loop%hat )%occupancy(1:loop%channels, &
                 1:loop%channels, is), &
                 h_atoms(loop%hatb)%occupancy(1:loop%channelsb, &
                 1:loop%channelsb, is), loop%channels, &
                 rot_vec_angles(hub_atom, 1:3, 1:3), rot_vec_angles(hub_atomb, &
                 1:3, 1:3), loop%hat, loop%hatb)
         endif

         row_counter = 0
         do row_proj = hub_proj_basis%first_on_atom(loop%theatom), &
              hub_proj_basis%first_on_atom(loop%theatom) + &
              hub_proj_basis%num_on_atom(loop%theatom) - 1
            row_counter = row_counter + 1
            col_counter = 0
            do col_proj = hub_proj_basis%first_on_atom(loop%theatomb), &
                 hub_proj_basis%first_on_atom(loop%theatomb) + &
                 hub_proj_basis%num_on_atom(loop%theatomb) - 1
               col_counter = col_counter + 1
               self_energy_element = site_self_energy_buffer(hub_atom, &
                    hub_atomb, row_counter, col_counter)
               call sparse_put_element(self_energy_element, &
                    self_energy%w_mat_v, row_proj, col_proj)
            end do
         end do

         if ((row_counter .ne. loop%channels) .or. (col_counter .ne. &
              loop%channelsb) ) then
            write(stdout, *) 'channels = ', loop%channels, loop%channelsb
            write(stdout, *) 'row_counter = ', row_counter
            write(stdout, *) 'col_counter = ', col_counter
            call utils_abort('Error in dmft_distribute_energy: row_counter or &
                 &col_counter is not equal to the number of channels')
         endif
      end do
   end do

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_distribute_self_energy'

 end subroutine

 subroutine dmft_build_self_energy_from_extended_atoms(is, ien, kien, &
      hub_atom, hub_atomb, site_self_energy_buffer, nmerge, merge_table, &
      merge_neigh, dimer_table, hub, hub_proj_basis, ffirst, llast, cluster, &
      reopen_sigma_file_each_step, force_reopen_file, scratchdir)

   !---------------------------------------------------------------------------!
   ! Description:                                                              !
   !---------------------------------------------------------------------------!

   use utils,             only: utils_abort
   use comms,             only: pub_my_proc_id, pub_on_root
   use hubbard_init,      only: h_species
   use rundat,            only: pub_debug, pub_debug_on_root
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par

   implicit none

   ! Arguments
   integer,              intent(in   ) :: is
   integer,              intent(in   ) :: ien
   integer,              intent(in   ) :: kien
   integer,              intent(  out) :: hub_atom
   integer,              intent(  out) :: hub_atomb
   complex(kind = DP),   intent(inout) :: site_self_energy_buffer(:,:,:,:)
   integer,              intent(in   ) :: nmerge
   integer, allocatable, intent(in   ) :: merge_table(:,:)
   integer, allocatable, intent(in   ) :: merge_neigh(:)
   integer, allocatable, intent(in   ) :: dimer_table(:,:)
   type(HUBBARD_MODEL),  intent(in   ) :: hub
   type(FUNC_BASIS),     intent(in   ) :: hub_proj_basis
   integer,              intent(in   ) :: ffirst
   integer,              intent(in   ) :: llast
   logical,              intent(in   ) :: cluster
   logical,              intent(in   ) :: reopen_sigma_file_each_step
   logical,              intent(in   ) :: force_reopen_file
   character(2000),      intent(in   ) :: scratchdir

   ! Local variables
   type(DMFT_LOOP_INFO) :: loop

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_build_self_energy_from_extended_atoms'

   do hub_atom = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)
      do hub_atomb = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

         call dmft_set_loop_variable(loop, hub_atom, hub_atomb, hub, &
              hub_proj_basis)

         if(.not.cluster .and. hub_atom /= hub_atomb) cycle
         if(cluster)then
            if((loop%hat /= loop%hatb .or. dimer_table(loop%hat, 1) == 0) &
                 .and. .not. dmft_dimer_connected(loop, dimer_table)) then
               cycle
            endif
         endif

         call dmft_print_loop_info(loop, hub)

         if (loop%channels .ne. (2*h_species(loop%sp)%hub_ang_mom + 1)) then
            call utils_abort('Error in hubbard_dmft_interface: &
                 &dmft_hub_atom_hub_projs and hub_ang_mom mismatch')
         endif

         if (loop%channels > size(site_self_energy_buffer, 3) .or. &
              loop%channelsb > size(site_self_energy_buffer, 4) .or. hub_atom &
              > size(site_self_energy_buffer, 1) .or. hub_atomb > &
              size(site_self_energy_buffer, 2)) then
            call utils_abort('ERROR onetep, buffer is too small')
         endif

         if(pub_debug_on_root) write(stdout, '(a)') 'DEBUG: reading self energy...'
         call dmft_read_self_energy(hub, hub_proj_basis, is, ien, kien, &
              hub_atom, hub_atomb, loop, site_self_energy_buffer, ffirst, &
              llast, nmerge, merge_table, merge_neigh, cluster, &
              reopen_sigma_file_each_step, force_reopen_file, scratchdir)
         if(pub_debug_on_root) write(stdout, '(a)') '...done'

      end do
   end do

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_build_self_energy_from_extended_atoms'

 end subroutine

 subroutine dmft_show_hubbard_atom(mat, hub_atom, matout, hub, hub_proj_basis, &
      showdiag, fullmat)

   !---------------------------------------------------------------------------!
   ! Prints details associated with a given hubbard atom                       !
   !---------------------------------------------------------------------------!

   use sparse,            only: sparse_convert
   use utils,             only: utils_abort, utils_alloc_check, &
                                utils_dealloc_check
   use comms,             only: pub_my_proc_id
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par
   use rundat,            only: pub_debug_on_root

   implicit none

   ! Arguments
   type(SPAM3),                               intent(in   ) :: mat
   integer,                                   intent(in   ) :: hub_atom
   complex(kind = DP), allocatable,           intent(inout) :: matout(:,:)
   type(HUBBARD_MODEL),                       intent(in   ) :: hub
   type(FUNC_BASIS),                          intent(in   ) :: hub_proj_basis
   logical,                                   intent(in   ) :: showdiag
   complex(kind = DP), allocatable, optional, intent(inout) :: fullmat(:,:)

   ! Local variables
   integer                         :: nn, hat, channels, theatom
   integer                         :: i1, i2, j1, j2, i, j
   complex(kind = DP), allocatable :: mat_square(:,:)
   real(kind = DP)                 :: totdiag
   integer                         :: hub_nproj
   integer                         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_show_hubbard_atom'

   hat = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
   channels = dmft_hub_atom_hub_projs(hat, hub, hub_proj_basis)
   hub_nproj = sum( (/( dmft_hub_atom_hub_projs(j, hub, hub_proj_basis), j = &
        1, par%nat_hub)/) )

   allocate(mat_square(hub_nproj, hub_nproj), stat = ierr)
   call utils_alloc_check('dmft_show_hubbard_atom', 'mat_square', ierr)
   call sparse_convert(mat_square, mat)

   if(showdiag)then
      write(stdout, *) 'total number of hubbard projections : ', hub_nproj
      totdiag = 0.0_DP
      do i = 1, hub_nproj
         totdiag = totdiag + real(mat_square(i, i))
      enddo
      write(stdout, *) 'real part of diagonal elements of hubbard projectors :'
      write(stdout, '(300f7.3)') (real(mat_square(i, i)), i = 1, hub_nproj)
      write(stdout, *) 'SUM IS : ', totdiag
   endif

   if(size(matout, 1) /= channels .or. size(matout, 2) /= channels)then
      write(stdout, *) 'shape matout : ', shape(matout)
      write(stdout, *) 'channels : ', channels
      call utils_abort("Error in dmft_show_hubbard_atom: shape mismatch")
   endif

   theatom = par%distr_atom(hub%orig(hat))
   i1 =   hub_proj_basis%first_on_atom(theatom)
   i2 = i1 + hub_proj_basis%num_on_atom(theatom) - 1
   write(stdout, *) 'first and last projections on atom : ', i1, i2

   matout = mat_square(i1:i2, i1:i2)
   if(present(fullmat)) then
      if (allocated(fullmat)) then
         deallocate(fullmat, stat = ierr)
         call utils_dealloc_check('dmft_show_hubbard_atom', 'fullmat', ierr)
      end if
      allocate(fullmat(hub_nproj, hub_nproj), stat = ierr)
      call utils_alloc_check('dmft_show_hubbard_atom', 'fullmat', ierr)
      fullmat = mat_square
   endif

   deallocate(mat_square, stat = ierr)
   call utils_dealloc_check('dmft_show_hubbard_atom', 'mat_square', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_show_hubbard_atom'

 end subroutine

 subroutine dmft_invert_green_function(greenf, self_energy, ngwf_basis, ham, &
      overlap, overlap_cplx, inv_overlap, energy, fermi_e, ien, is, ikpoint, &
      nkpoints, norm_kpoints, ffirst, kp_file_shift, splitk_first, &
      kpoint_data, matrixSigStat, nplanes, orb_in_plane, one_k_one_file)

   !---------------------------------------------------------------------------!
   ! Inverts omega - H - Sigma to obtain the Green's function of the whole     !
   ! system
   !---------------------------------------------------------------------------!

   use sparse,         only: sparse_copy
   use comms,          only: pub_on_root
   use rundat, only: pub_debug, pub_debug_on_root, &
        pub_dmft_points, pub_dmft_temp, pub_num_spins
   use function_basis, only: FUNC_BASIS

   implicit none

   ! Arguments
   type(DMFT_MATRIX),            intent(inout) :: greenf
   type(DMFT_MATRIX),            intent(inout) :: self_energy
   type(FUNC_BASIS),             intent(in   ) :: ngwf_basis
   type(SPAM3),                  intent(in   ) :: ham(pub_num_spins)
   type(SPAM3),                  intent(in   ) :: overlap
   type(SPAM3),                  intent(in   ) :: overlap_cplx
   type(SPAM3),                  intent(in   ) :: inv_overlap
   complex(kind = DP),           intent(in   ) :: energy
   real(kind = DP),              intent(in   ) :: fermi_e(2)
   integer,                      intent(in   ) :: ien
   integer,                      intent(in   ) :: is
   integer,                      intent(in   ) :: ikpoint
   integer,                      intent(in   ) :: nkpoints
   real(kind = DP), allocatable, intent(inout) :: norm_kpoints(:)
   integer,                      intent(in   ) :: ffirst
   integer,                      intent(in   ) :: kp_file_shift
   integer,                      intent(in   ) :: splitk_first
   type(DMFT_KPOINT_MATRICES),   intent(in   ) :: kpoint_data
   type(DEM), allocatable,       intent(inout) :: matrixSigStat(:)
   integer,                      intent(in   ) :: nplanes
   logical, allocatable,         intent(inout) :: orb_in_plane(:,:)
   logical,                      intent(in   ) :: one_k_one_file

   ! Local
   complex(kind=DP) :: det

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_invert_green_function'

   if(ien == pub_dmft_points-4 .or. ien == pub_dmft_points-5)then
      call sparse_copy(greenf%inv_backup_spam, greenf%inv_mat)
      call sparse_copy(greenf%inv_mat, overlap_cplx)
   endif

   if(pub_dmft_lin_scaling) then
      if(pub_on_root)write(stdout, '(a)') 'DEBUG: Inverting on CPU with linear &
           &scaling scheme'
      if(nkpoints == 1)then
         if(ien /= pub_dmft_points-3 .and. ien /= pub_dmft_points-4 .and. ien &
              /= pub_dmft_points-5) then
            call internal_invert_linear_scaling(greenf, energy, overlap_cplx, &
                 ham(is), self_energy%mat, is, fermi_e, ngwf_basis%num, &
                 self_energy%backup, matrixSigStat, update = .false., &
                 realorcomplexH = .true.)
         else
            call sparse_copy(greenf%mat, inv_overlap)
         endif
      else
         if(ien /= pub_dmft_points-3 .and. ien /= pub_dmft_points-4 .and. ien &
              /= pub_dmft_points-5)then
            call internal_invert_linear_scaling(greenf, energy, overlap_cplx, &
                 kpoint_data%ham, kpoint_data%self_energy, is, fermi_e, &
                 ngwf_basis%num, self_energy%backup, matrixSigStat, update = &
                 .false., realorcomplexH = .false.)
         else
            call sparse_copy(greenf%mat, kpoint_data%inv_overlap)
         endif
      endif
      if(pub_dmft_temp < 0.0_DP) then
         call dmft_show_dos(greenf, overlap, ngwf_basis, energy, is, ikpoint, &
              nkpoints, norm_kpoints, ien, ffirst, fermi_e, kp_file_shift, &
              splitk_first, nplanes, orb_in_plane, one_k_one_file)
      endif
   else
      ! #ifdef GPU_SPEEDUP
      ! if(pub_on_root)write(stdout, *) ' DEBUG: Inverting full matrix on GPU'
      ! if(pub_dmft_use_gpu_onlyme)then
      !    call internal_invert_full_gpu(greenf%mat, greenf%inv_mat, &
      !         overlap_cplx, energy, greenfunc = .true.)
      ! else
      !    call internal_invert_full(greenf%mat, greenf%inv_mat, ngwf_basis%num, &
      !         greenf%backup, greenf%backup_flag)
      !    if(pub_dmft_temp < 0.0_DP) then
      !       call dmft_show_dos(greenf, overlap, ngwf_basis, energy, is, &
      !            ikpoint, nkpoints, norm_kpoints, ien, ffirst, fermi_e, &
      !            kp_file_shift, splitk_first, nplanes, orb_in_plane, &
      !            one_k_one_file)
      !    endif
      ! endif
      ! #else
      ! if(pub_debug_on_root) write(stdout, '(a)') 'DEBUG: inverting full matrix &
      !      on CPU (no Open-MP)'
      call internal_invert_full(greenf%mat, greenf%inv_mat, ngwf_basis%num, &
           greenf%backup, greenf%backup_flag)
      if(pub_dmft_temp < 0.0_DP) then
         call dmft_show_dos(greenf, overlap, ngwf_basis, energy, is, ikpoint, &
              nkpoints, norm_kpoints, ien, ffirst, fermi_e, kp_file_shift, &
              splitk_first, nplanes, orb_in_plane, one_k_one_file)
      endif
      ! #endif
   endif

   if(ien == pub_dmft_points-4 .or. ien == pub_dmft_points-5)then
      call dmft_matmul_in_dense_format(greenf%backup_spam, greenf%mat, &
           greenf%inv_backup_spam, aba = .true.)
      call sparse_copy(greenf%mat, greenf%backup_spam) ! in greenf we have now
      ! (S^-1*Self*S^-1) or (S^-1*H*S^-1)
   endif

   if(nkpoints > 1)then
      if (pub_debug) then
         call dmft_print_spam3_info(greenf%mat)
         write(stdout, *) 'checking overlap_cplx'
         call dmft_print_spam3_info(overlap_cplx)
         write(stdout, *) 'checking k-point overlap'
         call dmft_print_spam3_info(kpoint_data%overlap)
         write(stdout, *) 'checking overlap'
         call dmft_print_spam3_info(overlap)
         write(stdout, *) 'checking self energy'
         call dmft_print_spam3_info(self_energy%mat)
         write(stdout, *) 'checking k-point self-energy'
         call dmft_print_spam3_info(kpoint_data%self_energy)
      end if
   endif


   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_invert_green_function'

 end subroutine

 subroutine dmft_extract_green_function_elements(greenf, h_atoms, hub, &
      hub_proj_basis, cell, elements, energy, fermi_e, is, ien, ffirst, &
      ikpoint, norm_kpoints, en_step_over_pi, kp_file_shift, energyzero, &
      rot_vec_angles, rot_vec_angles_split, splitk_first, split_wgv, &
      site_greenf_buffer, cluster, connection_table, one_k_one_file)

   !---------------------------------------------------------------------------!
   ! Extract Green's function elements from the local projectors to build the  !
   ! Green's function matrix                                                   !
   !---------------------------------------------------------------------------!

   use constants,         only: ANGSTROM, HARTREE_IN_EVS
   use hubbard_init,      only: h_species
   use simulation_cell,   only: CELL_INFO
   use comms, only: comms_barrier, pub_my_proc_id, pub_on_root, &
        pub_total_num_procs
   use sparse,            only: sparse_get_element
   use utils,             only: utils_abort, utils_unit, utils_int_to_str
   use rundat, only: pub_debug_on_root, pub_dmft_complex_freq, pub_dmft_points, &
        pub_dmft_rotate_green, pub_dmft_temp, pub_rootname
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use ion,               only: ELEMENT
   use parallel_strategy, only: par => pub_par

   implicit none

   ! Arguments
   type(DMFT_MATRIX),               intent(in   ) :: greenf
   type(ARRAY_OF_MATRICES),         intent(in   ) :: h_atoms(par%nat_hub)
   type(HUBBARD_MODEL),             intent(in   ) :: hub
   type(FUNC_BASIS),                intent(in   ) :: hub_proj_basis
   type(CELL_INFO),                 intent(in   ) :: cell
   type(ELEMENT),                   intent(in   ) :: elements(:)
   complex(kind = DP),              intent(in   ) :: energy
   real(kind = DP),                 intent(in   ) :: fermi_e(2)
   integer,                         intent(in   ) :: is
   integer,                         intent(in   ) :: ien
   integer,                         intent(in   ) :: ffirst
   integer,                         intent(in   ) :: ikpoint
   real(kind = DP), allocatable,    intent(in   ) :: norm_kpoints(:)
   real(kind = DP),                 intent(in   ) :: en_step_over_pi
   integer,                         intent(in   ) :: kp_file_shift
   integer,                         intent(in   ) :: energyzero
   real(kind = DP), allocatable,    intent(in   ) :: rot_vec_angles(:,:,:)
   real(kind = DP), allocatable,    intent(in   ) :: rot_vec_angles_split(:,:,:)
   integer,                         intent(in   ) :: splitk_first
   complex(kind = DP), allocatable, intent(in   ) :: split_wgv(:,:)
   complex(kind = DP), allocatable, intent(inout) :: site_greenf_buffer(:,:,:,:)
   logical,                         intent(in   ) :: cluster
   logical, allocatable,            intent(in   ) :: connection_table(:,:)
   logical,                         intent(in   ) :: one_k_one_file

   ! Local variables
   integer              :: numberHatom, funit
   logical              :: open_file
   integer              :: kkklll, kkkk, kk
   type(DMFT_LOOP_INFO) :: loop
   complex(kind = DP)   :: greenf_element
   real(kind = DP)      :: mytrace
   integer              :: row_proj
   integer              :: col_proj
   integer              :: row_counter
   integer              :: col_counter
   integer              :: hub_atom
   integer              :: hub_atomb

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_extract_green_function_elements'

   if (pub_dmft_split) then
      numberHatom = par%nat_hub
   else
      numberHatom = par%num_hub_atoms_on_proc(pub_my_proc_id)
   endif

   do hub_atom = 1, numberHatom
      do hub_atomb = 1, numberHatom

         call dmft_set_loop_variable(loop, hub_atom, hub_atomb, hub, &
              hub_proj_basis)

         if(.not.(ien == energyzero .and. pub_on_root .and. loop%channels == &
              loop%channelsb)) then

            if(.not.cluster .and. hub_atom /= hub_atomb) cycle

            if(cluster)then
               if( connection_table(loop%hat, loop%hatb) .or. &
                    (pub_dmft_plot_all_proj .and. hub_atom == hub_atomb) )then
               else
                  cycle
               endif
            endif
         end if

         if(loop%channels .ne. (2 * h_species(loop%sp)%hub_ang_mom + 1) )then
            call utils_abort('Error in dmft_extract_green_function_elements: &
                 &dmft_hub_atom_hub_projs and hub_ang_mom mismatch')
         endif

         if(loop%channels > size(site_greenf_buffer, 3) .or. loop%channelsb > &
              size(site_greenf_buffer, 4) .or. hub_atom > &
              size(site_greenf_buffer, 1) .or. hub_atomb > &
              size(site_greenf_buffer, 2))then
            call utils_abort('Error in dmft_extract_green_function_elements: &
                 &buffer is too small')
         endif

         row_counter = 0

         do row_proj = hub_proj_basis%first_on_atom(loop%theatom), &
              hub_proj_basis%first_on_atom(loop%theatom) + &
              hub_proj_basis%num_on_atom(loop%theatom) - 1

            row_counter = row_counter + 1
            col_counter = 0

            do col_proj = hub_proj_basis%first_on_atom(loop%theatomb), &
                 hub_proj_basis%first_on_atom(loop%theatomb) + &
                 hub_proj_basis%num_on_atom(loop%theatomb) - 1

               col_counter = col_counter + 1

               if(.not.pub_dmft_split)then
                  call sparse_get_element(greenf_element, greenf%w_mat_v, &
                       row_proj, col_proj)
               else
                  greenf_element = split_wgv(row_proj, col_proj)
               endif

               site_greenf_buffer(hub_atom, hub_atomb, row_counter, &
                    col_counter) = greenf_element

            end do
         end do

         !---rotate Green's Func------!
         if(pub_dmft_rotate_green) then
            if(loop%channels /= loop%channelsb) then
               call utils_abort('Error in &
                    &dmft_extract_green_function_elements: channels and &
                    &channelsb are different, and not compatible with local &
                    &rotation')
            endif

            ! Rotate green function back to local basis
            if(pub_dmft_split)then
               call dmft_rotate_to_local_basis( site_greenf_buffer(hub_atom, &
                    hub_atomb, 1:loop%channels, 1:loop%channelsb), &
                    h_atoms(loop%hat )%occupancy(1:loop%channels, &
                    1:loop%channels, is), &
                    h_atoms(loop%hatb)%occupancy(1:loop%channelsb, &
                    1:loop%channelsb, is), loop%channels, &
                    rot_vec_angles_split(hub_atom, 1:3, 1:3), &
                    rot_vec_angles_split(hub_atomb, 1:3, 1:3), loop%hat, &
                    loop%hatb, pub_dmft_rotate_green)
            else
               call dmft_rotate_to_local_basis( site_greenf_buffer(hub_atom, &
                    hub_atomb, 1:loop%channels, 1:loop%channelsb), &
                    h_atoms(loop%hat )%occupancy(1:loop%channels, &
                    1:loop%channels, is), &
                    h_atoms(loop%hatb)%occupancy(1:loop%channelsb, &
                    1:loop%channelsb, is), loop%channels, &
                    rot_vec_angles(hub_atom, 1:3, 1:3), &
                    rot_vec_angles(hub_atomb, 1:3, 1:3), loop%hat, loop%hatb, &
                    pub_dmft_rotate_green)
            endif
         endif
         !----------------------------!

         if ((row_counter .ne. loop%channels) .or. (col_counter .ne. &
              loop%channelsb) ) then
            write(stdout, *) 'channels = ', loop%channels, loop%channelsb
            write(stdout, *) 'row_counter = ', row_counter
            write(stdout, *) 'col_counter = ', col_counter
            call utils_abort('Error in dmft_extract_green_function_elements: &
                 &row_counter or col_counter mismatch with channels')
         endif

         !---dump local orbital DOS---!
         if(.not. pub_dmft_complex_freq) then
            if(hub_atom == hub_atomb)then
               open_file = .false.
               ! BUG CORRECTED JULY 8th 2014
               if ( cluster ) then
                  if(connection_table(loop%hat, loop%hatb)) open_file = .true.
               endif
               ! END BUG
               if ((.not.cluster) .or. open_file) then
                  funit = utils_unit()
                  if(.not.pub_dmft_split)then
                     if(pub_on_root)then
                        if((ien == ffirst .and. one_k_one_file) .or. &
                             (.not.one_k_one_file .and. ien == ffirst .and. &
                             ikpoint == 1))then
                           open(unit = funit, file = trim(pub_rootname) // &
                                '.dmft_ldos_atom_' // TRIM(ADJUSTL( &
                                utils_int_to_str(loop%hat))) // '_spin_' // &
                                TRIM(ADJUSTL(utils_int_to_str(is))), status = &
                                'replace')
                        else
                           open(unit = funit, file = trim(pub_rootname) // &
                                '.dmft_ldos_atom_' // TRIM(ADJUSTL( &
                                utils_int_to_str(loop%hat))) // '_spin_' // &
                                TRIM(ADJUSTL(utils_int_to_str(is))), status = &
                                'old', position = 'append')
                        endif

                        mytrace = 0.0_DP

                        do kk = 1, loop%channels
                           mytrace = mytrace + &
                                AIMAG(site_greenf_buffer(hub_atom, hub_atomb, &
                                kk, kk)*en_step_over_pi*norm_kpoints(ikpoint))
                        end do
                        write(funit, '(3i5, 300f14.6)') ien, loop%channels + &
                             1, pub_dmft_points, (real(energy-(fermi_e(is) + &
                             pub_dmft_chem_shift)))*HARTREE_IN_EVS, &
                             (AIMAG(site_greenf_buffer(hub_atom, hub_atomb, &
                             kk, kk)*en_step_over_pi*norm_kpoints(ikpoint)), &
                             kk = 1, loop%channels), mytrace
                        close(funit)
                     endif
                  else
                     do kkkk = 1, pub_total_num_procs
                        if(kkkk == pub_my_proc_id + 1)then
                           if (ien == ffirst .and. kkkk == 1 .and. &
                                (((.not.pub_dmft_splitk) .and. (ikpoint == 1 &
                                .or. one_k_one_file)) .or. (pub_dmft_splitk &
                                .and. ikpoint == splitk_first)))then
                              open(unit = funit, file = 'LOCAL_DOS_atom_' // &
                                   TRIM(ADJUSTL(utils_int_to_str(loop%hat))) &
                                   // '_spin' // TRIM(ADJUSTL( &
                                   utils_int_to_str(is))) // '_rank' // TRIM( &
                                   ADJUSTL(utils_int_to_str(pub_dmft_mpi_rank &
                                   + kp_file_shift))), STATUS = 'REPLACE')
                           else
                              open(unit = funit, file = 'LOCAL_DOS_atom_' // &
                                   TRIM(ADJUSTL(utils_int_to_str(loop%hat))) &
                                   // '_spin' // TRIM(ADJUSTL( &
                                   utils_int_to_str(is))) // '_rank' // TRIM( &
                                   ADJUSTL(utils_int_to_str(pub_dmft_mpi_rank &
                                   + kp_file_shift))), position = 'append', &
                                   STATUS = 'OLD' )
                           endif

                           mytrace = 0.0_DP

                           do kk = 1, loop%channels
                              mytrace = mytrace + &
                                   AIMAG(site_greenf_buffer(hub_atom, &
                                   hub_atomb, kk, kk) &
                                   *en_step_over_pi*norm_kpoints(ikpoint))
                           end do
                           write(funit, '(3i5, 300f14.6)') ien, loop%channels &
                                + 1, pub_dmft_points, &
                                (real(energy-(fermi_e(is) + &
                                pub_dmft_chem_shift)))*HARTREE_IN_EVS, &
                                (AIMAG(site_greenf_buffer(hub_atom, hub_atomb, &
                                kk, kk)*en_step_over_pi* &
                                norm_kpoints(ikpoint)), kk = 1, &
                                loop%channels), mytrace
                           close(funit)
                        endif
                        call comms_barrier
                     enddo
                  end if
               end if
            end if
         end if
         !----------------------------!
      end do
   end do

   ! if(pub_on_root .and. ien == energyzero)then
   !    funit = utils_unit()
   !    if(pub_dmft_temp > 0.0_DP)then
   !       open(unit = funit, file = &
   !            'projected_rotated_green_function_zero_energy_matsubara', form = &
   !            'unformatted')
   !    else
   !       open(unit = funit, file = &
   !            'projected_rotated_green_function_zero_energy_real', form = &
   !            'unformatted')
   !    endif
   !    write(funit) shape(site_greenf_buffer)
   !    write(funit) site_greenf_buffer
   !    write(funit) par%nat_hub
   !    write(funit) cell%a1%x/ANGSTROM, cell%a1%y/ANGSTROM, cell%a1%z/ANGSTROM
   !    write(funit) cell%a2%x/ANGSTROM, cell%a2%y/ANGSTROM, cell%a2%z/ANGSTROM
   !    write(funit) cell%a3%x/ANGSTROM, cell%a3%y/ANGSTROM, cell%a3%z/ANGSTROM
   !    do kkklll = 1, par%nat_hub
   !       write(funit) internal_coordinates_atom(kkklll, hub, elements)
   !    enddo
   !    close(funit)
   ! endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_extract_green_function_elements'

 end subroutine

 subroutine dmft_write_green_function(h_atoms, hub, hub_proj_basis, num_spins, &
      energy, fermi_e, is, ien, ikpoint, nkpoints, norm_kpoints, ffirst, &
      kp_file_shift, splitk_first, site_greenf_buffer, cluster, dimer_table, &
      merge_neigh, merge_table, nmerge, rot_vec_angles, rot_vec_angles_split, &
      one_k_one_file)

   !---------------------------------------------------------------------------!
   ! Writes the Green's function to file                                       !
   !---------------------------------------------------------------------------!

   use comms,             only: pub_my_proc_id
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par
   use rundat,            only: pub_debug, pub_debug_on_root

   implicit none

   ! Arguments
   type(ARRAY_OF_MATRICES),         intent(in   ) :: h_atoms(par%nat_hub)
   type(HUBBARD_MODEL),             intent(in   ) :: hub
   type(FUNC_BASIS),                intent(in   ) :: hub_proj_basis
   integer,                         intent(in   ) :: num_spins
   complex(kind = DP),              intent(in   ) :: energy
   real(kind = DP),                 intent(in   ) :: fermi_e(2)
   integer,                         intent(in   ) :: is
   integer,                         intent(in   ) :: ien
   integer,                         intent(in   ) :: ikpoint
   integer,                         intent(in   ) :: nkpoints
   real(kind = DP), allocatable,    intent(in   ) :: norm_kpoints(:)
   integer,                         intent(in   ) :: ffirst
   integer,                         intent(in   ) :: kp_file_shift
   integer,                         intent(in   ) :: splitk_first
   complex(kind = DP), allocatable, intent(in   ) :: site_greenf_buffer(:,:,:,:)
   logical,                         intent(in   ) :: cluster
   integer, allocatable,            intent(in   ) :: dimer_table(:,:)
   integer, allocatable,            intent(in   ) :: merge_neigh(:)
   integer, allocatable,            intent(in   ) :: merge_table(:,:)
   integer,                         intent(in   ) :: nmerge
   real(kind = DP), allocatable,    intent(in   ) :: rot_vec_angles(:,:,:)
   real(kind = DP), allocatable,    intent(in   ) :: rot_vec_angles_split(:,:,:)
   logical,                         intent(in   ) :: one_k_one_file

   ! Local variables
   integer              :: hub_atom
   integer              :: hub_atomb
   integer              :: i
   integer              :: numberHatom
   type(DMFT_LOOP_INFO) :: loop

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_write_green_function'

   numberHatom = par%num_hub_atoms_on_proc(pub_my_proc_id)

   if(pub_dmft_split)then
      numberHatom = par%nat_hub
   endif

   do hub_atom  = 1, numberHatom
      do hub_atomb = 1, numberHatom

         call dmft_set_loop_variable(loop, hub_atom, hub_atomb, hub, &
              hub_proj_basis)

         if(.not.cluster .and. hub_atom /= hub_atomb) cycle

         if(     cluster)then
            if( (loop%hat == loop%hatb .and. dimer_table(loop%hat, 1) /= 0) &
                 .or. dmft_dimer_connected(loop, dimer_table) ) then
            else
               cycle
            endif
         endif

         if(pub_debug_on_root)then
            write(stdout, '(a)')         'DEBUG: green function information'
            write(stdout, '(a,i4)')      ' hub_atom            : ', hub_atom
            write(stdout, '(a,i4)')      ' hub_atomb           : ', hub_atomb
            write(stdout, '(a,i4)')      ' hat                 : ', loop%hat
            write(stdout, '(a,i4)')      ' hatb                : ', loop%hatb
            write(stdout, '(a, 5e15.7)') ' green func diagonal : ', &
                 (REAL(site_greenf_buffer(hub_atom, hub_atomb, i, i)), i = 1, 5)
         endif

         call dmft_hubbard_greenf_info(loop, hub_atom, hub_atomb, &
              dmft_totn_atom_merge(hub_atom, hub, hub_proj_basis, merge_neigh, &
              merge_table, nmerge, cluster), dmft_totn_atom_merge(hub_atomb, &
              hub, hub_proj_basis, merge_neigh, merge_table, nmerge, cluster), &
              h_atoms, hub, hub_proj_basis, energy, fermi_e, num_spins, is, &
              ien, ikpoint, norm_kpoints, ffirst, nkpoints, kp_file_shift, &
              splitk_first, site_greenf_buffer, cluster, merge_neigh, &
              merge_table, nmerge, rot_vec_angles, rot_vec_angles_split, &
              one_k_one_file)

      enddo
   enddo

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_write_green_function'

 end subroutine

 subroutine dmft_read_self_energy(hub, hub_proj_basis, is, ien, kien, &
      hub_atom, hub_atomb, loop, site_self_energy_buffer, ffirst, llast, &
      nmerge, merge_table, merge_neigh, cluster, reopen_sigma_file_each_step, &
      force_reopen_file, scratchdir)

   !---------------------------------------------------------------------------!
   ! Reads the self energy from file                                           !
   !---------------------------------------------------------------------------!

   use function_basis, only: FUNC_BASIS
   use hubbard_build,  only: HUBBARD_MODEL
   use rundat,         only: pub_debug, pub_debug_on_root, pub_dmft_points, &
                             pub_dmft_temp
   use utils,          only: utils_int_to_str, utils_alloc_check, &
                             utils_dealloc_check

   implicit none

   ! Arguments
   type(HUBBARD_MODEL),  intent(in   ) :: hub
   type(FUNC_BASIS),     intent(in   ) :: hub_proj_basis
   integer,              intent(in   ) :: is
   integer,              intent(in   ) :: ien
   integer,              intent(in   ) :: kien
   type(DMFT_LOOP_INFO), intent(inout) :: loop
   integer,              intent(in   ) :: hub_atom
   integer,              intent(in   ) :: hub_atomb
   complex(kind = DP),   intent(inout) :: site_self_energy_buffer(:,:,:,:)
   integer,              intent(in   ) :: ffirst
   integer,              intent(in   ) :: llast
   integer,              intent(in   ) :: nmerge
   integer, allocatable, intent(in   ) :: merge_table(:,:)
   integer, allocatable, intent(in   ) :: merge_neigh(:)
   logical,              intent(in   ) :: cluster
   logical,              intent(in   ) :: reopen_sigma_file_each_step
   logical,              intent(in   ) :: force_reopen_file
   character(2000),      intent(in   ) :: scratchdir

   ! Local variables
   integer                         :: funit
   logical                         :: check
   integer                         :: l1, l2
   integer                         :: kk
   complex(kind = DP), allocatable :: self_temp(:,:)
   character(100)                  :: filename
   integer                         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_read_self_energy'

   ! cw: read in self-energy file here
   ! cw: self energy in orthogonalized basis (Self energy is diagonal)
   filename = 'sigma_output' // trim(adjustl(utils_int_to_str(loop%hat))) // &
        "_" // trim(adjustl( utils_int_to_str(loop%sp))) // "_" // &
        trim(adjustl(utils_int_to_str(is)))

   if(pub_dmft_local_scratch) then
      filename = trim(adjustl(scratchdir)) // trim(adjustl(filename))
   endif

   if( loop%hat /= loop%hatb ) filename = trim(adjustl(filename)) // '_dimer'

   inquire(file = trim(adjustl(filename)), exist = check)

   ! cw: No self energy for extrema: to extract in DMFT the non-orthogonal
   ! factor for the tail
   if(check)then

      ! ebl: this is foul, but TOSCAM relies on it so we're forced to follow
      ! suit
      funit = 40010 + loop%theatom + 2000*loop%theatomb + 110*loop%hat + &
           10101*is + 2503*loop%sp

      if(ien == 1 .and. pub_debug_on_root) write(stdout, "(/a)", &
           advance="no") 'Self energy file found'

      if(ien == ffirst .or. reopen_sigma_file_each_step .or. &
           force_reopen_file) then
         if(ien == 1 .and. pub_debug_on_root) write(stdout, *) &
              ': opening ', trim(adjustl(filename))
         open(unit = funit, file = trim(adjustl(filename)), form = &
              'unformatted')
      endif

      l1 = dmft_totn_atom_merge(hub_atom, hub, hub_proj_basis, merge_neigh, &
           merge_table, nmerge, cluster)
      l2 = dmft_totn_atom_merge(hub_atomb, hub, hub_proj_basis, merge_neigh, &
           merge_table, nmerge, cluster)

      allocate(self_temp(l1, l2), stat = ierr)
      call utils_alloc_check('dmft_read_self_energy', 'self_temp', ierr)

      if(.not.(reopen_sigma_file_each_step .or. force_reopen_file))then
         read(funit) self_temp
      else
         do kk = 1, kien
            read(funit) self_temp
         enddo
      endif

      ! cw: distribute self_temp in the blocks of site_self_energy_buffer
      call dmft_distribute(self_temp, hub_atom, hub_atomb, l1, l2, loop, &
           site_self_energy_buffer, hub, hub_proj_basis, nmerge, merge_table, &
           merge_neigh, cluster)

      deallocate(self_temp, stat = ierr)
      call utils_dealloc_check('dmft_read_self_energy', 'self_temp', ierr)

      if(pub_dmft_temp > 0.0_DP .and. ( ien == pub_dmft_points .or. ien == &
           pub_dmft_points-2 .or. ien == pub_dmft_points-3 .or. ien == &
           pub_dmft_points-4)) site_self_energy_buffer = 0.0_DP
      if(pub_dmft_temp < 0.0_DP .and. (ien == 1 .or. ien == pub_dmft_points &
           .or. ien == pub_dmft_points-2 .or. ien == pub_dmft_points-3 .or. &
           ien == pub_dmft_points-4)) site_self_energy_buffer = 0.0_DP
      if(ien == llast .or. reopen_sigma_file_each_step .or. &
           force_reopen_file)then
         if(ien == 1 .and. (pub_output_detail >= VERBOSE)) write(stdout, *) &
              'closing file : ', funit
         close(funit)
      endif
   else
      if(ien == 1 .and. pub_debug_on_root) write(stdout, "(/a)") &
           'No self energy file'
   endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_read_self_energy'

 end subroutine

 function internal_coordinates_atom(jj, hub, elements)

   !---------------------------------------------------------------------------!
   ! Returns the coordinates of hubbard atom jj (in Angstroms)                 !
   !---------------------------------------------------------------------------!

   use constants,     only: ANGSTROM
   use hubbard_build, only: HUBBARD_MODEL
   use ion,           only: ELEMENT

   implicit none

   ! Arguments
   integer,             intent(in   ) :: jj
   type(HUBBARD_MODEL), intent(in   ) :: hub
   type(ELEMENT),       intent(in   ) :: elements(:)

   ! Local variables
   integer         :: j
   real(kind = DP) :: internal_coordinates_atom(3)

   j = hub%orig(jj)
   internal_coordinates_atom(1) = elements(j)%centre%x/ANGSTROM
   internal_coordinates_atom(2) = elements(j)%centre%y/ANGSTROM
   internal_coordinates_atom(3) = elements(j)%centre%z/ANGSTROM

 end function

 integer function internal_atom_number_from_coordinates(coord, hub, elements)

   !---------------------------------------------------------------------------!
   ! Works out the index of an atom based on its coordinates                   !
   !---------------------------------------------------------------------------!

   use hubbard_build,     only: HUBBARD_MODEL
   use ion,               only: ELEMENT
   use parallel_strategy, only: par => pub_par

   implicit none

   ! Arguments
   real(kind = DP),     intent(in   ) :: coord(3)
   type(HUBBARD_MODEL), intent(in   ) :: hub
   type(ELEMENT),       intent(in   ) :: elements(:)

   ! Local variables
   integer         :: i, j
   real(kind = DP) :: v(3)

   j = 0
   do i = 1, par%nat_hub
      v = internal_coordinates_atom(i, hub, elements)
      if( ALL ( abs(coord-v) < 1.0E-3_DP ) )then
         j = i
         exit
      endif
   enddo
   internal_atom_number_from_coordinates = j
 end function


 integer function internal_atom_number_from_coordinates_shifted(coord, hub, &
      elements, cell)

   !---------------------------------------------------------------------------!
   ! Description:                                                              !
   !---------------------------------------------------------------------------!

   use constants,         only: ANGSTROM
   use hubbard_build,     only: HUBBARD_MODEL
   use ion,               only: ELEMENT
   use simulation_cell,   only: CELL_INFO
   use utils,             only: utils_abort
   use parallel_strategy, only: par => pub_par
   use rundat,            only: pub_debug_on_root

   implicit none

   ! Arguments
   real(kind = DP),     intent(in   ) :: coord(3)
   type(HUBBARD_MODEL), intent(in   ) :: hub
   type(ELEMENT),       intent(in   ) :: elements(:)
   type(CELL_INFO),     intent(in   ) :: cell

   ! Local variables
   integer         :: i, j, k, kk, i1, i2, i3
   real(kind = DP) :: v(3)
   real(kind = DP) :: real_lattice(3, 3)

   real_lattice(1, 1) = cell%a1%x/ANGSTROM
   real_lattice(1, 2) = cell%a1%y/ANGSTROM
   real_lattice(1, 3) = cell%a1%z/ANGSTROM
   real_lattice(2, 1) = cell%a2%x/ANGSTROM
   real_lattice(2, 2) = cell%a2%y/ANGSTROM
   real_lattice(2, 3) = cell%a2%z/ANGSTROM
   real_lattice(3, 1) = cell%a3%x/ANGSTROM
   real_lattice(3, 2) = cell%a3%y/ANGSTROM
   real_lattice(3, 3) = cell%a3%z/ANGSTROM

   ! ebl: work out which supercell we're in (note "int" rounds towards zero)
   i1 = INT(internal_scalprod(coord, real_lattice(1, :)) &
        / internal_norm_vector(real_lattice(1, :))**2)
   i2 = INT(internal_scalprod(coord, real_lattice(2, :)) &
        / internal_norm_vector(real_lattice(2, :))**2)
   i3 = INT(internal_scalprod(coord, real_lattice(3, :)) &
        / internal_norm_vector(real_lattice(3, :))**2)

   ! ebl: loop over neigbouring periodic images
   do i = -i1 - 1, -i1 + 1
      do j = -i2 - 1, -i2 + 1
         do k = -i3 - 1, -i3 + 1
            v =  real(i, kind=DP)*real_lattice(1, :)
            v = v + real(j, kind=DP)*real_lattice(2, :)
            v = v + real(k, kind=DP)*real_lattice(3, :)
            kk = internal_atom_number_from_coordinates(coord + v, hub, elements)
            if(kk /= 0) then
               internal_atom_number_from_coordinates_shifted = kk
               if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
                    &internal_atom_number_from_coordinates_shifted'
               return
            endif
         enddo
      enddo
   enddo

   ! ebl: if we get here, we haven't found a matching Hubbard atom
   call utils_abort("Error in &
        &internal_atom_number_from_coordinates_shifted: could not find &
        &a Hubbard atom with coordinates specified by mask_local_rotations")

   return

 end function

 subroutine dmft_correct_dimer_file(dimer_table, hub, cell, elements)

   !---------------------------------------------------------------------------!
   ! Updates the dimer connection table to match ONETEP's internal ordering of !
   ! atoms spread over procs                                                   !
   !---------------------------------------------------------------------------!

   use comms,             only: pub_on_root
   use constants,         only: ANGSTROM
   use file_handling,     only: file_handling_copy_file, file_handling_move_file
   use hubbard_build,     only: HUBBARD_MODEL
   use ion,               only: ELEMENT
   use parallel_strategy, only: par => pub_par
   use rundat,            only: pub_debug_on_root
   use simulation_cell,   only: CELL_INFO
   use utils,             only: utils_unit

   implicit none

   ! Arguments
   integer, allocatable, intent(in   ) :: dimer_table(:,:)
   type(HUBBARD_MODEL),  intent(in   ) :: hub
   type(CELL_INFO),      intent(in   ) :: cell
   type(ELEMENT),        intent(in   ) :: elements(:)

   ! Local variables
   logical         :: check_coord
   integer         :: i, k1, k2, alpha1, alpha2, sss
   real(kind = DP) :: coord1(3), coord2(3)
   integer         :: funit1
   integer         :: funit2
   integer         :: funit3

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_correct_dimer_file'

   funit1 = utils_unit()
   INQUIRE(file = 'mask_dimer_coord', EXIST = check_coord)
   if(check_coord)then
      if(pub_dmft_mpi_rank == 1 .and. pub_on_root) then
         call file_handling_move_file("mask_dimer", &
              "mask_dimer_original_coordinates", formatted = .true.)
         call file_handling_move_file("mask_dimer_coord", &
              "mask_dimer_coord_original_coordinates", formatted = .true.)
         open(unit = funit1, file = 'mask_dimer')
         do i = 1, par%nat_hub
            if(dimer_table(i, 1) > 0) write(funit1, *) dimer_table(i, 1), &
                 dimer_table(i, 2)
         enddo
         close(funit1)
      endif
   endif

   INQUIRE(file = 'mask_uniform_coord', EXIST = check_coord)
   if(check_coord)then
      if(pub_dmft_mpi_rank == 1 .and. pub_on_root) then
         open(unit = funit1, file = 'mask_uniform_coord')
         funit2 = utils_unit()
         open(unit = funit2, file = 'mask_uniform_corrected')
         funit3 = utils_unit()
         open(unit = funit3, file = 'mask_uniform')
         do i = 1, par%nat_hub
            read(funit1, *)   coord1
            read(funit1, *)   coord2
            read(funit3, *)   alpha1, alpha2
            if(pub_dmft_in_bohr)then
               coord1 = coord1/ANGSTROM
               coord2 = coord2/ANGSTROM
            endif
            sss = sign(1, alpha2 )
            write(funit2, *) &
                 internal_atom_number_from_coordinates_shifted(coord1, hub, &
                 elements, cell), &
                 sss*internal_atom_number_from_coordinates_shifted(coord2, &
                 hub, elements, cell)
         enddo
         close(funit1)
         close(funit2)
         close(funit3)
         call file_handling_move_file("mask_uniform", "mask_uniform_original", &
              formatted = .true.)
         call file_handling_move_file("mask_uniform_coord", &
              "mask_uniform_coord_original", formatted = .true.)
         call file_handling_copy_file("mask_uniform_corrected", &
              "mask_uniform", formatted = .true.)
      endif
   endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_correct_dimer_file'

 end subroutine

 subroutine dmft_define_cluster(rot_vec_angles, rot_vec_angles_split, &
      dimer_table, nmerge, merge_table, merge_neigh, connection_table, hub, &
      hub_proj_basis, cell, elements, cluster)

   !---------------------------------------------------------------------------!
   ! Finds pairs of atoms for cluster DMFT                                     !
   !---------------------------------------------------------------------------!

   use comms,             only: comms_barrier, pub_my_proc_id, pub_on_root
   use constants,         only: ANGSTROM
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL, hubbard_cmp_rot
   use ion,               only: ELEMENT
   use parallel_strategy, only: par => pub_par
   use rundat,            only: pub_debug_on_root, pub_debug, pub_dmft_write
   use simulation_cell,   only: CELL_INFO
   use utils,             only: utils_abort, utils_unit, utils_int_to_str, &
                                utils_alloc_check, utils_dealloc_check, &
                                utils_assert

   implicit none

   ! Arguments
   real(kind = DP), allocatable, intent(inout) :: rot_vec_angles(:,:,:)
   real(kind = DP), allocatable, intent(inout) :: rot_vec_angles_split(:,:,:)
   integer, allocatable,         intent(inout) :: dimer_table(:,:)
   integer,                      intent(in   ) :: nmerge
   integer, allocatable,         intent(inout) :: merge_table(:,:)
   integer, allocatable,         intent(inout) :: merge_neigh(:)
   logical, allocatable,         intent(inout) :: connection_table(:,:)
   type(HUBBARD_MODEL),          intent(in   ) :: hub
   type(FUNC_BASIS),             intent(in   ) :: hub_proj_basis
   type(CELL_INFO),              intent(in   ) :: cell
   type(ELEMENT),                intent(in   ) :: elements(:)
   logical,                      intent(in   ) :: cluster

   ! Local variables
   logical         :: allsites, check, check_coord, check_merge_file
   integer         :: i, j, k, l, i2, k2, l2, ii, jj, ijk_, iostatus
   integer         :: ll, lll, chan, lll_, llll_, free_unit
   real(kind = DP) :: v_temp(3), dtemp, v1(3), v2(3)
   real(kind = DP) :: coord(3), mat_temp(3,3)
   integer         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_define_cluster'

   allsites = .true.

   if(allocated(rot_vec_angles)) then
      deallocate(rot_vec_angles, stat = ierr)
      call utils_dealloc_check('dmft_define_cluster', 'rot_vec_angles', ierr)
   end if
   if(allocated(rot_vec_angles_split)) then
      deallocate(rot_vec_angles_split, stat = ierr)
      call utils_dealloc_check('dmft_define_cluster', 'rot_vec_angles_split', &
           ierr)
   end if

   allocate(rot_vec_angles(par%nat_hub, 3, 3), stat = ierr)
   call utils_alloc_check('dmft_define_cluster', 'rot_vec_angles', ierr)
   allocate(rot_vec_angles_split(par%nat_hub, 3, 3), stat = ierr)
   call utils_alloc_check('dmft_define_cluster', 'rot_vec_angles_split', ierr)

   if (pub_debug) then
      write(stdout, '(a, i4, a, i4)') 'There are ', &
           par%num_hub_atoms_on_proc(pub_my_proc_id), ' Hubbard atoms on proc',&
           pub_my_proc_id
   end if

   rot_vec_angles = 0.0d0
   rot_vec_angles_split = 0.0d0

   inquire(file = 'mask_local_rotations', exist = check)

   call comms_barrier

   if(check)then
      free_unit = utils_unit()
      open(unit = free_unit, file = 'mask_local_rotations', form = &
           'unformatted')
      do ijk_ = 1, par%nat_hub

         read(free_unit, IOSTAT = iostatus) lll, coord(:), mat_temp(:,:)
         if (iostatus < 0) call utils_abort('Error in dmft_define_cluster: the &
              &mask_local_rotations file does not have enough entries')
         if (iostatus > 0) call utils_abort('Error in dmft_define_cluster: the &
              &mask_local_rotations file does not contain the expected contents')

         if(pub_dmft_in_bohr)then
            coord = coord/ANGSTROM
         endif
         if(lll == 0)then
            allsites = .false.
            cycle
         endif

         if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, '(a)') &
              'Applying rotation for hubbard atom ' //  trim(utils_int_to_str(ijk_))
         j = internal_atom_number_from_coordinates_shifted(coord, hub, &
              elements, cell)
         if(j /= ijk_)then
            write(stdout, *) 'WARNING, a shift of coordinate was done by &
                 &onetep!!!!'
            write(stdout, *) 'might give problems with rotations when in mpi'
         endif

         call utils_assert(j > 0  .and. j <= par%nat_hub, "Error in &
              &dmft_define_cluster: encountered non-sensical Hubbard index")

         rot_vec_angles_split(j, :, :) = mat_temp

         ! ebl: transforming "j" into hub_atom notation
         j = dmft_hubbard_on_my_proc(j)

         if(j /= 0)then
            rot_vec_angles(j, :, :) = mat_temp

            chan = dmft_hub_atom_hub_projs(par%hub_atoms_on_proc(j, &
                 pub_my_proc_id), hub, hub_proj_basis)
            ! hub_projs take 1..nat argument, so hat notation

            if((pub_output_detail >= VERBOSE) .and. pub_on_root)then
               write(stdout, '(a, i5, a)') 'Atom ', lll, ' has rotation'
               do lll_ = 1, 3
                  write(stdout, '(3f10.3)') (rot_vec_angles(j, lll_, llll_), &
                       llll_ = 1, 3)
               enddo
               write(stdout, *) 'j, pub_my_proc_id : ', j, pub_my_proc_id
               write(stdout, *) 'checking rotations for [x] channels : ', chan
               write(stdout, *) 'TEST, B*B^T (maxval): ', maxval(abs( &
                    matmul(hubbard_cmp_rot(rot_vec_angles(j, :, :), &
                    (chan-1)/2, chan, j), &
                    transpose(hubbard_cmp_rot(rot_vec_angles(j, :, :), &
                    (chan-1)/2, chan, j) )) ))
               write(stdout, *) 'TEST, B*B^T (minval): ', minval(abs( &
                    matmul(hubbard_cmp_rot(rot_vec_angles(j, :, :), &
                    (chan-1)/2, chan, j), &
                    transpose(hubbard_cmp_rot(rot_vec_angles(j, :, :), &
                    (chan-1)/2, chan, j) )) ))
            endif
         endif
      enddo

      do lll_ = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)
         if(maxval(abs(rot_vec_angles(lll_, :, :))) < 1.0E-4_DP .and. allsites)then
            write(stdout, '(a,i3)') ' WARNING: rotation matrix was not defined for &
                 &site ', lll_
            v_temp = internal_coordinates_atom( par%hub_atoms_on_proc(lll_, &
                 pub_my_proc_id), hub, elements)
            write(stdout, *) ' coordinates of this site : ', v_temp
         endif
      enddo
      close(free_unit)

   endif

   if(cluster)then
      if(allocated(dimer_table)) then
         deallocate(dimer_table, stat = ierr)
         call utils_dealloc_check('dmft_define_cluster', 'dimer_table', ierr)
         deallocate(merge_table, stat = ierr)
         call utils_dealloc_check('dmft_define_cluster', 'merge_table', ierr)
         deallocate(merge_neigh, stat = ierr)
         call utils_dealloc_check('dmft_define_cluster', 'merge_neigh', ierr)
      end if
      allocate(dimer_table(par%nat_hub, 2), stat = ierr)
      call utils_alloc_check('dmft_define_cluster', 'dimer_table', ierr)
      allocate(merge_table(par%nat_hub, nmerge), stat = ierr)
      call utils_alloc_check('dmft_define_cluster', 'merge_table', ierr)
      allocate(merge_neigh(par%nat_hub), stat = ierr)
      call utils_alloc_check('dmft_define_cluster', 'merge_neigh', ierr)

      if((pub_output_detail >= VERBOSE) .and. pub_on_root)then
         write(stdout, *) 'NUMBER OF PROJECTIONS ON ATOMS'
         do j = 1, par%nat_hub
            write(stdout, '(a, 2i7)') 'ATOM, NPROJ : ', j, &
                 dmft_hub_atom_hub_projs(j, hub, hub_proj_basis)
         enddo
         write(stdout, *) 'READING MASK_DIMER'
      endif

      free_unit = utils_unit()
      open(unit = free_unit, file = 'mask_dimer')
      dimer_table = 0
      do i = 1, par%nat_hub
         read(free_unit, *, iostat = iostatus) j, jj
         if (iostatus /= 0) exit
         if(j > 0)then
            dimer_table(j, 1) = j
            dimer_table(j, 2) = jj
            if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, &
                 '(a, 200i4)') 'DIMER TABLE : ', (dimer_table(j, jj), jj = 1, &
                 2)
         endif
         if(jj > 0 .and. j > 0 .and. dimer_table(j, 1) > 0 .and. &
              dimer_table(j, 2) > 0)then
            dtemp = &
                 internal_norm_vector(internal_coordinates_atom(dimer_table(j, &
                 1), hub, elements)-internal_coordinates_atom(dimer_table(j, &
                 2), hub, elements))
            if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, &
                 *) 'DIMER BOND LENGTH : ', dtemp
         endif
      enddo
      close(free_unit)

      INQUIRE(file = 'mask_dimer_coord', EXIST = check_coord)
      if(check_coord)then
         dimer_table = 0
         free_unit = utils_unit()
         open(unit = free_unit, file = 'mask_dimer_coord')
         do i = 1, par%nat_hub
            read(free_unit, *) v1
            read(free_unit, *) v2
            if(pub_dmft_in_bohr)then
               v1 = v1/ANGSTROM
               v2 = v2/ANGSTROM
            endif
            j = internal_atom_number_from_coordinates_shifted(v1, hub, &
                 elements, cell)
            dimer_table(j, 1) = j
            dimer_table(j, 2) = &
                 internal_atom_number_from_coordinates_shifted(v2, hub, &
                 elements, cell)
         enddo
         close(free_unit)
      endif

      INQUIRE(FILE = 'mask_merge', EXIST = check_merge_file)
      merge_table = 0
      free_unit = utils_unit()
      if(check_merge_file) open(unit = free_unit, file = 'mask_merge')
      merge_neigh = 0
      do i = 1, par%nat_hub
         if(check_merge_file) then
            read(free_unit, *, iostat = iostatus) (merge_table(i, j), j = 1, &
                 nmerge)
            if (iostatus /= 0) exit
         else
            merge_table(i, :) = 0
            merge_table(i, 1) = i
         endif
         if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, &
              '(a, 200i4)') 'MERGE TABLE : ', (merge_table(i, j), j = 1, &
              nmerge)
         merge_neigh(i) = 0
         do j = 2, nmerge
            if(merge_table(i, j) == 0)then
               merge_neigh(i) = j-1
               exit
            else
               merge_neigh(i) = j
            endif
         enddo
      enddo

      if(check_merge_file) close(free_unit)

      if((pub_output_detail >= VERBOSE) .and. pub_on_root) write(stdout, *) &
           'build connection table'

      if(allocated(connection_table)) then
         deallocate(connection_table, stat = ierr)
         call utils_dealloc_check('dmft_define_cluster', 'connection_table', &
              ierr)
      end if

      allocate(connection_table(par%nat_hub, par%nat_hub), stat = ierr)
      call utils_alloc_check('dmft_define_cluster', 'connection_table', ierr)

      connection_table = .false.

      do i = 1, par%nat_hub
         if(dimer_table(i, 1) /= 0)then
            do i2 = 1, 2
               do k = 1, merge_neigh(i)
                  l = merge_table(i, k)
                  jj = dimer_table(i, i2)
                  if(jj /= 0)then
                     do k2 = 1, merge_neigh(jj)
                        l2 = merge_table(jj, k2)
                        connection_table(l, l2) = .true.
                        connection_table(l2, l) = .true.
                        if(dmft_hubbard_on_my_proc(l) == 0 .or. &
                             dmft_hubbard_on_my_proc(l2) == 0)then
                           call utils_abort('cluster DMFT, connection table &
                                &not compatible with distribution of atoms on &
                                &mpi procs, proc_id : ' // &
                                utils_int_to_str(pub_my_proc_id))
                        endif
                     enddo
                  endif
               enddo
            enddo
         endif
      enddo

      !connection_table = .true. !TESTING-DEBUG
      if((pub_output_detail >= VERBOSE) .and. pub_on_root)then
         write(stdout, *) 'TABLE'
         write(stdout, '(100i0)') (l, l = 1, par%nat_hub)
         do i = 1, par%nat_hub
            write(stdout, *) i, (connection_table(i, l), l = 1, par%nat_hub)
         enddo
      endif

   else
      if((pub_output_detail >= VERBOSE) .and. pub_on_root) then
         write(stdout, "(a/)") ' mask_dimer file not present; performing &
              &single-site DMFT'
      end if
   endif

   ! if((pub_output_detail >= VERBOSE) .and. pub_on_root .and. pub_dmft_write) then
   !    free_unit = utils_unit()
   !    open(unit = free_unit, file = 'positions_hubbard_atoms')
   !    do ii = 1, par%nat_hub
   !       if( hub%orig(ii) /= par%orig_atom(par%distr_atom(hub%orig(ii))) )then
   !          write(stdout, *) ' WARNING: position definition and access to &
   !               &element(j) array is ill defined'
   !          write(stdout, *) 'this is hub%orig(ii) : ', hub%orig(ii)
   !          write(stdout, *) 'this is &
   !               &par%orig_atom(par%distr_atom(hub%orig(ii))) : ', &
   !               par%orig_atom(par%distr_atom(hub%orig(ii)))
   !       endif
   !       v_temp = internal_coordinates_atom(ii, hub, elements)
   !       write(free_unit, '(i5, i5, 10f15.4)') ii, hub%species_number(ii), &
   !            v_temp
   !    enddo
   !    close(free_unit)
   ! endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_define_cluster'

 end subroutine

 subroutine dmft_build_my_dimers(dimer_table, cluster, hub, hub_proj_basis, &
      cell, elements, nmerge, merge_table, merge_neigh)

   !---------------------------------------------------------------------------!
   ! Builds connection tables for dimers                                       !
   !---------------------------------------------------------------------------!

   use comms,             only: pub_my_proc_id, pub_on_root
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use ion,               only: ELEMENT
   use parallel_strategy, only: par => pub_par
   use rundat,            only: pub_debug_on_root, pub_dmft_write
   use simulation_cell,   only: CELL_INFO
   use utils,             only: utils_unit, utils_abort

   implicit none

   ! Arguments
   logical,              intent(in   ) :: cluster
   integer, allocatable, intent(in   ) :: dimer_table(:,:)
   type(HUBBARD_MODEL),  intent(in   ) :: hub
   type(FUNC_BASIS),     intent(in   ) :: hub_proj_basis
   type(CELL_INFO),      intent(in   ) :: cell
   type(ELEMENT),        intent(in   ) :: elements(:)
   integer,              intent(in   ) :: nmerge
   integer, allocatable, intent(in   ) :: merge_table(:,:)
   integer, allocatable, intent(in   ) :: merge_neigh(:)

   ! Local variables
   integer         :: funit_sites, funit_dimers, funit_dists, funit_blocks, &
                      funit_torb
   integer         :: hub_atom, hub_atomb
   integer         :: hat, hatb, nat
   real(kind = DP) :: dist

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_build_my_dimers'

   nat = par%nat_hub

   if (.not. pub_dmft_write) then
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_build_my_dimers immediately since dmft_write is set to FALSE'
      return
   end if

   if(pub_on_root)then
      if(cluster)then
         call utils_abort("Code in dmft_build_my_dimers is under development")

         ! ebl: As implemented currently, multiple MPI tasks may be writing to
         ! the same file at once (if it weren't for the if (pub_on_root), which
         ! means that atoms are being skipped. This needs correcting.

         ! if(pub_debug_on_root) write(stdout, '(a)') &
         !      'DEBUG: writing cluster files'

         ! funit_sites = utils_unit()
         ! open(unit = funit_sites, file = 'compound_sites')
         ! funit_dimers = utils_unit()
         ! open(unit = funit_dimers, file = 'compound_dimers')
         ! funit_dists = utils_unit()
         ! open(unit = funit_dists, file = 'compound_distances')
         ! funit_blocks = utils_unit()
         ! open(unit = funit_blocks, file = 'compound_blocks')

         ! do hub_atom = 1, nat
         !    if (pub_dmft_split) then
         !       hat = hub_atom
         !    else
         !       hat = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
         !    end if

         !    do hub_atomb = 1, nat
         !       if (pub_dmft_split) then
         !          ! CWTODO: should this be hub_atomb?
         !          hatb = hub_atom
         !       else
         !          hatb = par%hub_atoms_on_proc(hub_atomb, pub_my_proc_id)
         !       end if

         !       if( dimer_table(hat, 2) == hatb ) then
         !          write(funit_blocks, *) dmft_hub_atom_hub_projs(hat, hub, &
         !               hub_proj_basis), dmft_hub_atom_hub_projs(hatb, hub, &
         !               hub_proj_basis)
         !       endif
         !    enddo
         ! enddo

         ! do hub_atom = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)
         !    hat   = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
         !    write(funit_sites, '(i3, 10f11.4)') hub_atom, &
         !         internal_coordinates_atom(hat, hub, elements)

         !    do hub_atomb = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)
         !       hatb  = par%hub_atoms_on_proc(hub_atomb, pub_my_proc_id)
         !       dist = internal_distance(hat, hatb, hub, cell, elements)

         !       if(hub_atom /= hub_atomb) write(funit_dists, '(a, i4, i4, &
         !            &f10.4)') 'ATOMS : ', hub_atom, hub_atomb, dist

         !       if( dimer_table(hat, 2) == hatb ) then
         !          write(funit_dimers, '(a, i4, 10f11.4)')' A = : ', hub_atom, &
         !               internal_coordinates_atom(hat, hub, elements)
         !          write(funit_dimers, '(a, i4, 10f11.4)')' B = : ', hub_atomb, &
         !               internal_coordinates_atom(hatb, hub, elements)
         !       endif
         !    enddo
         ! enddo

         ! close(funit_sites)
         ! close(funit_dimers)
         ! close(funit_dists)
         ! close(funit_blocks)
      endif

      if(pub_debug_on_root) write(stdout, '(a)') 'DEBUG: writing &
           &compound_tot_orb file'

      funit_torb = utils_unit()
      open(unit = funit_torb, file = 'compound_tot_orb')
      do hat = 1, nat
         write(funit_torb, *) dmft_hub_atom_hub_projs(hat, hub, hub_proj_basis)
         ! write(funit_torb, *) dmft_totn_atom_merge(hub_atom, hub, &
         !      hub_proj_basis, merge_neigh, merge_table, nmerge, cluster)
      enddo
      close(funit_torb)

   endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_build_my_dimers'

 end subroutine

 subroutine dmft_read_write_hamiltonian(eigen_en, n_occ, eigen_en_input, &
      n_occ_input, h_atoms, hub, hub_proj_basis, num_spins, check)

   !---------------------------------------------------------------------------!
   ! Reads or writes the Hamiltonian to file                                   !
   !---------------------------------------------------------------------------!

   use comms,             only: comms_barrier, pub_my_proc_id, pub_on_root
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par
   use rundat,            only: pub_debug_on_root, pub_dmft_read, &
        pub_dmft_write, pub_rootname
   use utils,             only: utils_alloc_check, utils_unit, &
        utils_int_to_str, utils_alloc_check, utils_dealloc_check

   implicit none

   ! Arguments
   real(kind = DP), allocatable, intent(inout) :: eigen_en(:,:)
   integer, allocatable,         intent(inout) :: n_occ(:)
   real(kind = DP), allocatable, intent(in   ) :: eigen_en_input(:,:)
   integer,                      intent(in   ) :: n_occ_input(:)
   type(ARRAY_OF_MATRICES),      intent(inout) :: h_atoms(par%nat_hub)
   type(HUBBARD_MODEL),          intent(in   ) :: hub
   type(FUNC_BASIS),             intent(in   ) :: hub_proj_basis
   integer,                      intent(in   ) :: num_spins
   logical,                      intent(in   ) :: check

   ! Local variables
   real(kind = DP), allocatable :: temp_alloc(:,:)
   integer                      :: iboundary
   integer                      :: store_unit
   integer                      :: hat
   integer                      :: channels
   integer                      :: n1_eigen
   integer                      :: n2_eigen
   integer                      :: n3_eigen
   integer                      :: ierr
   integer                      :: hub_atom
   integer                      :: ii

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_read_write_hamiltonian'

   store_unit = utils_unit()

   if(.not.check .or. .not. pub_dmft_read)then

      ! Copy over _input files
      if (allocated(eigen_en)) then
         deallocate(eigen_en, stat = ierr)
         call utils_dealloc_check('dmft_read_write_hamiltonian', 'eigen_en', &
              ierr)
      end if

      if(allocated(n_occ)) then
         deallocate(n_occ, stat = ierr)
         call utils_dealloc_check('dmft_read_write_hamiltonian', 'n_occ', ierr)
      end if

      allocate(eigen_en(size(eigen_en_input, 1), size(eigen_en_input, 2)), &
           stat = ierr)
      call utils_alloc_check('dmft_read_write_hamiltonian', 'eigen_en', ierr)
      allocate(n_occ(size(n_occ_input, 1)), stat = ierr)
      call utils_alloc_check('dmft_read_write_hamiltonian', 'n_occ', ierr)
      eigen_en = eigen_en_input
      n_occ = n_occ_input

      if(pub_on_root) then
         if (check) then
            ! File exists (but pub_dmft_read is false)
            write(stdout, '(a)', advance="no") trim(pub_rootname) // &
                 &'.eigen is present; '
         else
            ! File does not exist
            write(stdout, '(a)', advance="no") trim(pub_rootname) // &
                 &'.eigen is not present; '
         end if
      end if

      if (pub_dmft_write .and. .not. check) then
         if (pub_on_root) then
            write(stdout, '(a)') "writing"
            open(unit = store_unit, file = trim(pub_rootname) // &
                 '.eigen', form = 'unformatted')
            write(store_unit) size(eigen_en_input, 1), &
                 size(eigen_en_input, 2), size(n_occ_input, 1)
            write(store_unit) n_occ, eigen_en
            close(store_unit)
         end if
      else
         if (pub_on_root) write(stdout, '(a)') "it will not be (re)generated"
      end if

      call comms_barrier

      if (pub_dmft_write) then
         do hub_atom = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)
            hat = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
            channels = dmft_hub_atom_hub_projs(hat, hub, hub_proj_basis)

            open(unit = store_unit, file = trim(pub_rootname) // '.occupancy' &
                 // trim(adjustl(utils_int_to_str(hat))), form = 'unformatted')
            write(store_unit) channels, num_spins
            write(store_unit) h_atoms(hat)%occupancy(1:channels, 1:channels, 1)
            if(num_spins == 2)then
               write(store_unit) h_atoms(hat)%occupancy(1:channels, &
                    1:channels, 2)
            endif
            close(store_unit)
         enddo
      end if
      call comms_barrier
   else

      if((pub_output_detail >= VERBOSE) .and. pub_on_root) &
           write(stdout, '(a)') trim(pub_rootname) // '.eigen is present; &
           &reading'
      open(unit = store_unit, file = trim(pub_rootname) // '.eigen', form = &
           'unformatted')

      read(store_unit) n1_eigen, n2_eigen, n3_eigen
      if (allocated(eigen_en)) then
         deallocate(eigen_en, stat = ierr)
         call utils_dealloc_check('dmft_read_write_hamiltonian', 'eigen_en', ierr)
         deallocate(n_occ, stat = ierr)
         call utils_dealloc_check('dmft_read_write_hamiltonian', 'n_occ', ierr)
      end if
      allocate(eigen_en(n1_eigen, n2_eigen), stat = ierr)
      call utils_alloc_check('dmft_read_write_hamiltonian', 'eigen_en', ierr)
      allocate(n_occ(n3_eigen), stat = ierr)
      call utils_alloc_check('dmft_read_write_hamiltonian', 'n_occ', ierr)
      read(store_unit) n_occ, eigen_en

      if(size(eigen_en, 2) == 1 .and. num_spins == 2)then
         write(stdout, *) ' WARNING: spin-polarised DMFT calculation requested but the DFT '
         write(stdout, *) '          calculation is not spin polarised'
         write(stdout, *) '          Will proceed with spin-polarised DMFT calculation'
         if (allocated(temp_alloc)) then
            deallocate(temp_alloc, stat = ierr)
            call utils_dealloc_check('dmft_read_write_hamiltonian', 'temp_alloc', ierr)
         end if
         allocate(temp_alloc(size(eigen_en, 1), size(eigen_en, 2)), stat = ierr)
         call utils_alloc_check('dmft_read_write_hamiltonian', 'temp_alloc', ierr)
         temp_alloc = eigen_en
         deallocate(eigen_en, stat = ierr)
         call utils_dealloc_check('dmft_read_write_hamiltonian', 'eigen_en', ierr)
         allocate(eigen_en(size(temp_alloc, 1), 2), stat = ierr)
         call utils_alloc_check('dmft_read_write_hamiltonian', 'eigen_en', ierr)
         eigen_en(:, 1) = temp_alloc(:, 1)
         eigen_en(:, 2) = temp_alloc(:, 1)
         deallocate(temp_alloc, stat = ierr)
         call utils_dealloc_check('dmft_read_write_hamiltonian', 'temp_alloc', ierr)
      endif

      close(store_unit)

      if(pub_dmft_split)then
         iboundary = par%nat_hub
      else
         iboundary = par%num_hub_atoms_on_proc(pub_my_proc_id)
      endif

      do hub_atom = 1, iboundary
         if(pub_dmft_split)then
            hat = hub_atom
         else
            hat = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
         endif
         channels = dmft_hub_atom_hub_projs(hat, hub, hub_proj_basis)
         open(unit = store_unit, file = trim(pub_rootname) // '.occupancy' // &
              trim(adjustl(utils_int_to_str(hat))), form = 'unformatted')
         read(store_unit) channels, ii
         read(store_unit) h_atoms(hat)%occupancy(1:channels, 1:channels, 1)
         if(ii /= num_spins)then
            write(stdout, *) 'not the right number of spin species, should be &
                 &: ', ii
            write(stdout, *) 'check variable dmft_paramagnetic or number of &
                 &spin species in onetep input file'
         endif
         if(num_spins == 2)then
            if(ii == 2)then
               read(store_unit) h_atoms(hat)%occupancy(1:channels, 1:channels, &
                    2)
            else
               h_atoms(hat)%occupancy(1:channels, 1:channels, 2) = &
                    h_atoms(hat)%occupancy(1:channels, 1:channels, 1)
            endif
         endif
         close(store_unit)
      enddo
   endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_read_write_hamiltonian'

 end subroutine

 subroutine dmft_print_hubbard_atom_info(eigen_en, n_occ, h_atoms, hub, &
      hub_proj_basis, num_spins, rot_vec_angles)

   !---------------------------------------------------------------------------!
   ! Prints various information associated with a given hubbard atom           !
   !---------------------------------------------------------------------------!

   use utils,             only: utils_assert, utils_unit
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL, hubbard_diago_occupancy
   use comms,             only: pub_my_proc_id, pub_on_root
   use parallel_strategy, only: par => pub_par
   use rundat, only: pub_debug_on_root, pub_dmft_smear, pub_dmft_smear_T, &
        pub_dmft_smear_eta, pub_dmft_smear_w

   implicit none

   ! Arguments
   real(kind = DP), allocatable, intent(in   ) :: eigen_en(:,:)
   integer, allocatable,         intent(in   ) :: n_occ(:)
   type(ARRAY_OF_MATRICES),      intent(in   ) :: h_atoms(par%nat_hub)
   type(HUBBARD_MODEL),          intent(in   ) :: hub
   type(FUNC_BASIS),             intent(in   ) :: hub_proj_basis
   integer,                      intent(in   ) :: num_spins
   real(kind = DP), allocatable, intent(in   ) :: rot_vec_angles(:,:,:)

   ! Local variables
   real(kind = DP) :: temp_dens(100, 100)
   integer         :: hub_atom, i
   ! integer         :: unit_dens, unit_dens_rot
   integer         :: channels
   integer         :: hat
   integer         :: j

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_print_hubbard_atom_info'

   if (pub_output_detail /= VERBOSE) return

   if (pub_on_root) then
      if (pub_debug_on_root) then
         write(stdout, *) 'SPECTRUM'
         write(stdout, *) 'n_occ              : ', n_occ(:)
         write(stdout, *) 'shape eigen        : ', shape(eigen_en)
         write(stdout, *) 'eigen up           : ', eigen_en(n_occ(1), 1)
         if(num_spins == 2)then
            write(stdout, *) 'eigen dn           : ', eigen_en(n_occ(2), 2)
         endif
         ! write(stdout, *) 'pub_dmft_smear     : ', pub_dmft_smear
         ! write(stdout, *) 'pub_dmft_smear_T   : ', pub_dmft_smear_T
         ! write(stdout, *) 'pub_dmft_smear_w   : ', pub_dmft_smear_w
         ! write(stdout, *) 'pub_dmft_smear_eta : ', pub_dmft_smear_eta
      end if
      ! unit_dens = utils_unit()
      ! unit_dens_rot = utils_unit()
      ! open(unit = unit_dens, file = 'density_spin_up_non_rotated')
      ! open(unit = unit_dens_rot, file = 'density_spin_up_rotated')

      do hub_atom = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)
         hat = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
         channels = dmft_hub_atom_hub_projs(hat, hub, hub_proj_basis)
         call utils_assert(channels <= size(temp_dens), 'Error in &
              &dmft_print_hubbard_atom_info: temp_dens too small')
         if (pub_debug_on_root) then
            write(stdout, '(a,i4,a)') 'Atom ', hat, 'with spin 1 occupancy matrix'
            do i = 1, channels
               write(stdout, '(200f10.5)') ( h_atoms(hat)%occupancy(i, j, 1), &
                    j = 1, channels )
            enddo
         end if
         ! write(unit_dens, '(200f8.4)') ( h_atoms(hat)%occupancy(i, i, 1), i = &
         !      1, channels ), sum( (/( h_atoms(hat)%occupancy(i, i, 1), i = 1, &
         !      channels )/) )
         call hubbard_diago_occupancy(h_atoms(hat)%occupancy(1:channels, &
              1:channels, 1), channels, hat, rot_vec_angle = &
              rot_vec_angles(hub_atom, 1:3, 1:3), tmp_dens = temp_dens)
         ! write(unit_dens_rot, '(200f8.4)') ( temp_dens(i, i), i = 1, channels &
         !      ), sum( (/( temp_dens(i, i), i = 1, channels )/) )
      enddo

      ! close(unit_dens)
      ! close(unit_dens_rot)

      if(num_spins == 2)then

         ! open(unit = unit_dens, file = 'density_spin_dn_non_rotated')
         ! open(unit = unit_dens_rot, file = 'density_spin_dn_rotated')

         do hub_atom = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)
            hat = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
            channels = dmft_hub_atom_hub_projs(hat, hub, hub_proj_basis)
            if (pub_debug_on_root) then
               write(stdout, '(a,i4,a)') 'Atom ', hat, ' with spin 2 occupancy matrix'
               do i = 1, channels
                  write(stdout, '(200f10.5)') (h_atoms(hat)%occupancy(i, j, 2),&
                       j = 1, channels)
               enddo
            end if
            ! write(unit_dens, '(200f8.4)') ( h_atoms(hat)%occupancy(i, i, 2), i &
            !      = 1, channels ), sum( (/( h_atoms(hat)%occupancy(i, i, 2), i &
            !      = 1, channels )/) )
            call hubbard_diago_occupancy(h_atoms(hat)%occupancy(1:channels, &
                 1:channels, 2), channels, hat, rot_vec_angle = &
                 rot_vec_angles(hub_atom, 1:3, 1:3), tmp_dens = temp_dens)
            ! write(unit_dens_rot, '(200f8.4)') ( temp_dens(i, i), i = 1, &
            !      channels ), sum( (/( temp_dens(i, i), i = 1, channels )/) )
         enddo

         ! close(unit_dens)
         ! close(unit_dens_rot)

      endif

   endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_print_hubbard_atom_info'

 end subroutine

 subroutine dmft_clean_dmft_energy(ham, overlap, ngwf_basis, sc_data, &
      matrixDensKernel, matrixAsig, matrixSigStat, denskern_tmp, &
      denskern_tmpc, dmft_Energy, total_energy, edc_correction, &
      E_sum_eigenvalues_DFT, fermi_e, is, num_spins, ikpoint, nkpoints, kpt, &
      totkpoint, lastkpoint, nmu, nmu_step, chem_shift, vectorNmu, vectorMu, &
      target_N, norm_kpoints, dmft_splitk_iter, optimisation_parallel, &
      dmft_energy_cor_, dmft_energy_cor)

   !---------------------------------------------------------------------------!
   ! Prepares arrays used in the calculation of the energy for a new iteration !
   ! of the chemical potential                                                 !
   !---------------------------------------------------------------------------!

   use dense, only: dense_axpy, dense_convert, dense_create, dense_destroy, &
        dense_product, dense_scale, dense_trace, dense_copy, dense_put_element
   use sparse, only: sparse_axpy, sparse_copy, sparse_create, sparse_destroy, &
        sparse_get_element, sparse_index_length, sparse_num_rows, &
        sparse_scale, sparse_trace, sparse_transpose
   ! #ifdef GPU_SPEEDUP
   ! use fortran_cuda,   only: matmulcuda_r
   ! #endif
   use utils,          only: utils_abort, utils_busy_wait, utils_unit, &
        utils_banner, utils_int_to_str, utils_alloc_check, &
        utils_dealloc_check, utils_qc_print
   use comms, only: comms_barrier, comms_bcast, comms_reduce, pub_my_proc_id, &
        pub_on_root, pub_total_num_procs
   use rundat, only: pub_debug, pub_debug_on_root, &
        pub_dmft_kernel, pub_dmft_mu_diff_max, pub_dmft_mu_order, &
        pub_dmft_nkpoints, pub_dmft_paramagnetic, pub_dmft_skip_energy, &
        pub_dmft_temp, pub_num_spins, pub_dmft_ks_shift, pub_dmft_write, &
        pub_print_qc
   use kernel,         only: DKERN
   use function_basis, only: FUNC_BASIS

   implicit none

   ! Arguments
   type(SPAM3),                  intent(in   ) :: ham(pub_num_spins)
   type(SPAM3),                  intent(in   ) :: overlap
   type(FUNC_BASIS),             intent(in   ) :: ngwf_basis
   type(DMFT_SC_MATRICES),       intent(in   ) :: sc_data
   type(SPAM3),                  intent(inout) :: matrixDensKernel
   type(DEM), allocatable,       intent(in   ) :: matrixAsig(:)
   type(DEM), allocatable,       intent(in   ) :: matrixSigStat(:)
   type(DKERN),                  intent(inout) :: denskern_tmp
   type(DKERN),                  intent(inout) :: denskern_tmpc
   real(kind = DP),              intent(inout) :: dmft_Energy(2)
   real(kind = DP),              intent(inout) :: total_energy
   real(kind = DP),              intent(inout) :: edc_correction
   real(kind = DP),              intent(inout) :: E_sum_eigenvalues_DFT(2)
   real(kind = DP),              intent(in   ) :: fermi_e(2)
   integer,                      intent(in   ) :: is
   integer,                      intent(in   ) :: num_spins
   integer,                      intent(in   ) :: ikpoint
   integer,                      intent(in   ) :: nkpoints
   real(kind = DP), allocatable, intent(in   ) :: kpt(:,:)
   real(kind = DP),              intent(in   ) :: totkpoint
   integer,                      intent(in   ) :: lastkpoint
   integer,                      intent(in   ) :: nmu
   integer,                      intent(inout) :: nmu_step
   real(kind = DP),              intent(inout) :: chem_shift
   real(kind = DP),              intent(inout) :: vectorNmu(100)
   real(kind = DP),              intent(inout) :: vectorMu(100)
   real(kind = DP),              intent(in   ) :: target_N
   real(kind = DP), allocatable, intent(in   ) :: norm_kpoints(:)
   integer,                      intent(inout) :: dmft_splitk_iter
   logical,                      intent(in   ) :: optimisation_parallel
   real(kind = DP),              intent(inout) :: dmft_energy_cor_(2)
   real(kind = DP), optional,    intent(in   ) :: dmft_energy_cor

   ! Local variables
   integer                         :: ii, i, j, m, jj, num_denskern, np
   integer                         :: iiistart, iiistep, iistart, iistep, funit
   real(kind = DP)                 :: totN, totE, beta, chem, lambda
   real(kind = DP)                 :: ccback, cc, slope_mu
   real(kind = DP)                 :: ttscreen, mymu
   complex(kind = DP)              :: ccc, ccbackc
   logical                         :: skip_energy_part
   real(kind = DP)                 :: contrib_energy
   real(kind = DP)                 :: contrib_N
   real(kind = DP), parameter      :: slope_mu_mixing = 1.0_DP ! was 0.3
   real(kind = DP), parameter      :: slope_mu_max_step = 0.10_DP ! was 0.01
   real(kind = DP), allocatable    :: tmpr(:,:)
   real(kind = DP), allocatable    :: tmpr2(:,:)
   complex(kind = DP), allocatable :: tmpc(:,:)
   complex(kind = DP), allocatable :: tmpc2(:,:)
   type(SPAM3)                     :: matrixDensKernel_t
   type(SPAM3)                     :: tail_spam
   real(kind = DP)                 :: el_r
   complex(kind = DP)              :: el_c
   logical                         :: check
   integer                         :: ierr
   type(DEM)                       :: prefactor
   type(DEM)                       :: tmp_Asig
   type(DEM)                       :: tmp_dem1
   type(DEM)                       :: tmp_dem2
   type(DEM)                       :: tmp_dem3
   logical, parameter              :: mu_via_bisection = .false.

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_clean_dmft_energy'

   ii = ngwf_basis%num
   ! Allocation of matrices
   allocate(tmpc(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_clean_dmft_energy', 'tmpc', ierr)
   allocate(tmpc2(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_clean_dmft_energy', 'tmpc2', ierr)
   allocate(tmpr(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_clean_dmft_energy', 'tmpr', ierr)
   allocate(tmpr2(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_clean_dmft_energy', 'tmpr2', ierr)

   beta = 1.0_DP/pub_dmft_temp
   chem = fermi_e(is) + pub_dmft_chem_shift

   call sparse_scale(matrixDensKernel, pub_dmft_temp)

   !--------------------------!
   if(pub_dmft_mpi_size /= 1 .and. .not.pub_dmft_splitk)then

      ! if (.not. pub_debug) then
      call utils_abort("Reached point in code where nfs sync used to appear; &
           &the code needs to be updated to use these settings.")
      ! call sync_via_nfs(dmft_Energy(is), 'summing_dmft_energy' // &
      !      trim(adjustl(utils_int_to_str(is))))
      ! call sync_via_nfs(matrixDensKernel, 'sum_matrixDensKernel', &
      !      sparse_index_length(matrixDensKernel))

      ! else
      !    funit = utils_unit()
      ! open(file = 'dens_kern_tmp' //
      ! trim(adjustl(utils_int_to_str(pub_dmft_mpi_rank))), unit = funit, &
      ! status = "replace", position = 'rewind', form = 'unformatted')
      !    if(nkpoints == 1)then
      !       write(funit) matrixDensKernel%dmtx, dmft_Energy(is)
      ! write(funit) matrixDensKernel%dmtx, dmft_Energy(is) !to flush the
      ! file...
      !    else
      !       write(funit) matrixDensKernel%zmtx, dmft_Energy(is)
      ! write(funit) matrixDensKernel%zmtx, dmft_Energy(is) !to flush the
      ! file...
      !    endif
      !    close(funit)
      ! #ifdef DMFT_SAFE_NFS
      !    call utils_busy_wait(1000)
      ! #endif
      !    do
      ! #ifdef DMFT_SAFE_NFS
      !       call utils_busy_wait(500)
      ! #else
      !       call utils_busy_wait(100)
      ! #endif
      ! ! ressys = system(" echo 'I am done' >
      ! onetep_dmft_confirmation_dmft_density_kernel_" //
      ! trim(adjustl(utils_int_to_str(pub_dmft_mpi_rank))))
      !       funit = utils_unit()
      !       if(pub_dmft_mpi_rank == 1) then
      ! open(funit, file = "onetep_dmft_confirmation_dmft_density_kernel_" //
      ! trim(adjustl(utils_int_to_str(pub_dmft_mpi_rank))))
      !          write(funit, *) "I am done"
      !          close(funit)
      !       end if
      ! ! if(pub_dmft_mpi_rank == 1) ressys = system("echo `ls -l
      ! onetep_dmft_confirmation_dmft_density_kernel_* 2 > &1 | grep -v 'No
      ! such' | wc -l` > onetep_dmft_confirmation_dmft_density_count")

      ! open(unit = funit, file =
      ! 'onetep_dmft_confirmation_dmft_density_count')
      !       jj = 0
      !       read(funit, *) jj
      !       close(funit)
      !       write(stdout, *) 'N PROCESS ARE DONE = ', jj
      !       if(jj == pub_dmft_mpi_size) then
      !          write(stdout, *) 'ALL ARE DONE'
      !          exit
      !       endif
      !    enddo
      ! #ifndef DMFT_SAFE_NFS
      !    call utils_busy_wait(1000)
      ! #else
      !    call utils_busy_wait(2000)
      ! #endif!
      !    call dense_scale(matrixDensKernel, 0.0_DP)
      !    dmft_Energy(is) = 0.0_DP

      !    do i = 1, pub_dmft_mpi_size
      !       funit = utils_unit()
      !       open(file = 'dens_kern_tmp' // trim(adjustl( &
      !            utils_int_to_str(i))), unit = funit, action = 'read', &
      !            position = 'rewind', form = 'unformatted')
      !       if(nkpoints == 1)then
      !          read(funit) tmpr, cc
      !          matrixDensKernel%dmtx  = matrixDensKernel%dmtx  + tmpr
      !       else
      !          read(funit) tmpc, cc
      !          matrixDensKernel%zmtx = matrixDensKernel%zmtx + tmpc
      !       endif
      !       close(unit = funit)
      !       write(*, *) cc
      !       dmft_Energy(is)  = dmft_Energy(is)  + cc
      !    enddo
      ! end if
   endif

   !--------------------------!
   if(ikpoint <= nkpoints)then
      if(abs(norm_kpoints(ikpoint)) > 1.0E-4_DP)then
         !first sync frequencies, then add tail for this k point
         call sparse_create(tail_spam, matrixDensKernel)
         call dense_convert(tail_spam, sc_data%tail)
         call sparse_axpy(matrixDensKernel, tail_spam, 1.0_DP)
         call sparse_destroy(tail_spam)
      endif
   endif

   totN = 0.0_DP

   if( optimisation_parallel )Then
      iistart = pub_my_proc_id
      iistep = pub_total_num_procs
   else
      iistart = 0
      iistep = 1
   endif

   totN = sparse_trace(matrixDensKernel, sc_data%overlap)

   ! if (pub_debug) then
   !    if(nkpoints > 1)then
   ! write(stdout, *) 'MINVAL-MAXVAL MATRIXAsig : ', minval(abs(matrixAsig)),
   ! maxval(abs(matrixAsig))
   ! write(stdout, *) 'MINVAL-MAXVAL MATRIXAAC : ',
   ! minval(abs(sc_data%eigenvecs_c)), maxval(abs(sc_data%eigenvecs_c))
   ! write(stdout, *) 'MINVAL-MAXVAL DENSKERNEL : ',
   ! minval(abs(matrixDensKernelc)), maxval(abs(matrixDensKernelc))
   ! write(stdout, *) 'MINVAL-MAXVAL DENSK TAIL : ',
   ! minval(abs(sc_data%tail_c)), maxval(abs(sc_data%tail_c))
   ! write(stdout, *) 'MINVAL-MAXVAL overlapc : ',
   ! minval(abs(sc_data%overlap_c)), maxval(abs(sc_data%overlap_c))
   !    end if
   ! end if

   if(ikpoint <= nkpoints)then
      contrib_N = totN*norm_kpoints(ikpoint)
   else
      contrib_N = 0.0_DP
   endif

   if(pub_dmft_splitk)then
      dmft_splitk_iter = dmft_splitk_iter + 1
      if(pub_on_root) write(stdout, *) 'this is sync number : ', &
           dmft_splitk_iter
      call utils_abort("Reached point in code where nfs sync used to appear; &
           &the code needs to be updated to use these settings.")
      ! if(pub_on_root) call sync_via_nfs(contrib_N, 'summing_up_contribN' // &
      !      trim(adjustl(utils_int_to_str(dmft_splitk_iter))))
      ! call comms_barrier
      ! call comms_bcast(0, contrib_N)
      if(pub_on_root) write(stdout, *) 'TOTAL CHARGE THIS BATCH K POINTS : ', &
           contrib_N
   endif

   vectorNmu(nmu_step) = vectorNmu(nmu_step) + contrib_N

   if(pub_on_root .and. is == num_spins) write(stdout, "(a,f10.6,a,f10.4,a)") &
        repeat(" ",19)//"|     ", vectorMu(nmu_step), "     |     ", &
        vectorNmu(nmu_step), "     |"

   if(nmu_step >= nmu) then

      if (pub_debug_on_root) then
         if (num_spins == 1) then
            write(stdout, "(/a)") utils_banner("-", "Energy summary")
         else
            write(stdout, "(/a)") utils_banner("-", "Energy summary for spin "&
                 //trim(adjustl(utils_int_to_str(is))))
         end if
      end if
      if(pub_debug_on_root) write(stdout, "(a,f16.9)") ' < H_U > from dynamical&
           & effect                : ', dmft_Energy(is)

      ! 1/2 Tr ( Sigma_an * G_an ) !
      if(.not.pub_dmft_use_gpu .and. .not.pub_dmft_use_gpu_partially)then
         ! Create diagonal matrix
         call dense_create(prefactor, ii, ii, nkpoints/=1)
         call dense_scale(prefactor, 0.0_DP)
         do m = 1, ii
            lambda = sc_data%eigenvals(is, m) - chem
            call dense_put_element(1.0_DP / ( -lambda ) * ( 1.0_DP - &
                 internal_dexpc(-beta*lambda) ) / ( 1.0_DP + &
                 internal_dexpc(-beta*lambda) ), prefactor, m, m)
         end do

         ! Creating temporary matrices for product
         call dense_create(tmp_Asig, ii, ii, nkpoints/=1)
         call dense_create(tmp_dem1, ii, ii, nkpoints/=1)
         call dense_create(tmp_dem2, ii, ii, nkpoints/=1)
         call dense_create(tmp_dem3, ii, ii, nkpoints/=1)

         ! Performing trace of product
         if (nkpoints == 1) then
            call dense_product(tmp_dem1, prefactor, sc_data%eigenvecs(is), "T")
         else
            call dense_product(tmp_dem1, prefactor, sc_data%eigenvecs(is), "C")
         end if
         call dense_copy(tmp_Asig, matrixAsig(is), cmplx_to_real=(nkpoints==1))
         call dense_product(tmp_dem2, tmp_dem1, tmp_Asig)
         call dense_product(tmp_dem3, tmp_dem2, sc_data%eigenvecs(is))
         if (nkpoints == 1) then
            cc = dense_trace(tmp_dem3)
         else
            ccc = dense_trace(tmp_dem3)
         end if

         ! Tidying up
         call dense_destroy(prefactor)
         call dense_destroy(tmp_Asig)
         call dense_destroy(tmp_dem1)
         call dense_destroy(tmp_dem2)
         call dense_destroy(tmp_dem3)

         if(pub_dmft_mpi_size /= 1 .and. .not.pub_dmft_splitk) then
            call utils_abort("Reached point in code where nfs sync used to &
                 &appear; the code needs to be updated to use these settings.")
            ! if(nkpoints == 1)then
            !    call sync_via_nfs(cc, 'summing_up_cc' )
            ! else
            !    call sync_via_nfs(ccc, 'summing_up_ccc')
            ! endif
         endif

      else
         call utils_abort("Error in dmft_clean_dmft_energy: GPU not yet &
              &implemented")
         ! write(stdout,*) 'compute 1/2 Tr ( Sigma_an * G_an ) on GPU'
         ! cc = 0.0_DP
         ! ccc = 0.0_DP

         ! if(pub_dmft_mpi_rank == 1 .or. pub_dmft_splitk)then
         !    do m = 1, ii
         !       lambda = sc_data%eigenvals(is, m) - chem
         !       dtempd = 1.0_DP / ( -lambda ) * ( 1.0_DP - &
         !            internal_dexpc(-beta*lambda) ) / ( 1.0_DP + &
         !            internal_dexpc(-beta*lambda) )
         !       if(nkpoints == 1)then
         !          tmpr(m, :) = sc_data%eigenvecs(is)%dmtx(:, m)*dtempd
         !       else
         !          tmpc(m, :) = sc_data%eigenvecs(is)%zmtx(:, m)*dtempd
         !       endif
         !    enddo
         !    write(stdout, *) 'MATMUL CUDA ', ii
         !    ! #ifdef GPU_SPEEDUP
         !    ! if(nkpoints == 1)then
         !    !    call matmulcuda_r( sc_data%eigenvecs(is)%dmtx, tmpr, tmpr2, ii, &
         !    !         ii, ii) ! tmpr2 = sc_data%eigenvecs(is, j, m)*tmpr(m, i)
         !    ! else
         !    !    call matmulcuda_c(conjg(sc_data%eigenvecs(is)%zmtx), tmpc, &
         !    !         tmpc2, ii, ii, ii)
         !    ! endif
         !    ! #else
         !    call utils_abort("Error in dmft_clean_dmft_energy: should not &
         !         &arrive here")
         !    ! #endif

         !    if(nkpoints == 1)then
         !       cc = 0.0_DP
         !       do i = 1, ii
         !          do j = 1, ii
         !             cc = cc + matrixAsig(is)%dmtx(i, j)*tmpr2(j, i)
         !          enddo
         !       enddo
         !       write(stdout, *) 'CC FROM CUDA : ', cc
         !    else
         !       ccc = 0.0_DP
         !       do i = 1, ii
         !          do j = 1, ii
         !             ccc = ccc + matrixAsig(is)%zmtx(i, j)*tmpc2(i, j)
         !          enddo
         !       enddo
         !    endif
         ! endif
         ! if(pub_dmft_mpi_size /= 1 .and. .not.pub_dmft_splitk) then
         !    write(stdout, *) 'I AM DONE RANK : ', pub_dmft_mpi_rank
         !    call utils_abort("Reached point in code where nfs sync used to &
         !         &appear; the code needs to be updated to use these settings.")
         !    ! if(nkpoints == 1)then
         !    !    call sync_via_nfs(cc, 'summing_up_cc')
         !    ! else
         !    !    call sync_via_nfs(ccc, 'summing_up_ccc')
         !    ! endif
         ! endif
      endif

      if(nkpoints == 1)then
         if(pub_debug_on_root) write(stdout, "(a,f16.9)") &
              ' < H_U > Sigma_an * G_an                      :', 0.25_DP*cc
      else
         if(pub_debug_on_root) write(stdout, "(a,f16.9)") &
              ' < H_U > Sigma_an * G_an                      :', 0.25_DP*ccc
      endif
      if(nkpoints == 1)then
         dmft_Energy(is) = dmft_Energy(is) + 0.25_DP * cc
      else
         dmft_Energy(is) = dmft_Energy(is) + 0.25_DP * ccc
      endif

      call dense_create(tmp_dem1, ii, ii, iscmplx=.false.)
      call dense_create(tmp_dem2, ii, ii, iscmplx=.false.)
      call dense_copy(tmp_dem1, sc_data%tail, cmplx_to_real=(nkpoints==1))
      call dense_copy(tmp_dem2, matrixSigStat(is), cmplx_to_real=(nkpoints==1))
      totE = 0.5_DP*dense_trace(tmp_dem1, tmp_dem2)
      call dense_destroy(tmp_dem1)
      call dense_destroy(tmp_dem2)

      if(pub_debug_on_root) write(stdout, "(a,f16.9)") ' < H_U > Sigma(oo) &
           &                           :', totE

      dmft_Energy(is) = dmft_Energy(is) + totE

      if(pub_debug_on_root) write(stdout, "(a,f16.9)") ' < H_U > TOTAL     &
           &                           :', dmft_Energy(is)

      totE = sparse_trace(matrixDensKernel, sc_data%ham(is))
      if(pub_debug_on_root) write(stdout, "(a,f16.9)") ' Trace( G_dmft*H_DFT ) &
           &                       :', totE

      ! IF NO ENERGY FROM DENS KERNEL CORRECTION DUE TO CORRELATIONS
      if(present(dmft_energy_cor))then
         if(.not.pub_dmft_KS_shift)then
            E_sum_eigenvalues_DFT = 0.0_DP
            totE = 0.0_DP
         endif
         total_energy = 0.0_DP
      endif

      dmft_Energy(is) = dmft_Energy(is) + totE - E_sum_eigenvalues_DFT(is) + &
           total_energy/2.0_DP

      edc_correction = 0.0_DP
      inquire(file = 'edc_total', exist = check)
      if(check)then
         funit = utils_unit()
         open(unit = funit, file = 'edc_total')
         read(funit, *) edc_correction
         close(funit)
      endif

      ! ebl: Quality control - print total energy. This is a good test as it
      ! incorporates both the calculation of the Green's function and the
      ! reading/writing of the density kernel
      if (pub_on_root .and. pub_print_qc) then
         call utils_qc_print("E_DMFT("//trim(adjustl(utils_int_to_str(is))) // &
              ")", dmft_Energy(is))
      end if

   end if

   !skipping energy part if nmu_step /= nmu

   ! if(pub_on_root)then
   !    write(stdout, *) 'CHEM POTENTIAL : ', fermi_e(is) + pub_dmft_chem_shift
   !    write(stdout, *) 'TOTAL DENSITY  : ', vectorNmu(nmu_step)
   ! endif

   ! if (pub_debug) then
   !    if(nkpoints == 1)then
   !       num_denskern = sparse_num_rows(denskern%kern%m(1, PUB_1K))
   !       if(num_denskern /= sparse_num_rows(matrixDensKernel))then
   !          write(stdout, *) 'ERROR : dim do not match denskernel in dmft!'
   !          write(stdout, *) 'shape denskern ?      : ', num_denskern
   ! write(stdout, *) 'shape matrix dmft ? : ', shape(matrixDensKernel)
   !          call utils_abort('ERROR : dim do not match denskernel in dmft!')
   !       endif
   !    endif
   ! end if


   if(ikpoint >= lastkpoint .and. (is == num_spins .or. &
        pub_dmft_impose_chem_spin))then
      if(abs( vectorNmu(nmu_step) - target_N ) < pub_dmft_mu_diff_max .and. &
           pub_dmft_skip_energy)then
         if(pub_on_root)then
            write(stdout, *) 'WARNING : TARGET REACHED - SKIPPING OTHER MU &
                 &STEPS'
            write(stdout, *) ' target   : ', target_N
            write(stdout, *) ' obtained : ', vectorNmu(nmu_step)
         endif
         vectorNmu(nmu) = vectorNmu(nmu_step)
         vectorMu(nmu) = vectorMu(nmu_step)
         nmu_step = nmu
         chem_shift = 0.0_DP
      endif
   endif

   if(nmu_step == nmu .and. pub_dmft_kernel /= 0) then
      if(ikpoint <= nkpoints)then
         contrib_energy = &
              (dmft_Energy(is)-edc_correction/2.0_DP)*norm_kpoints(ikpoint)
      else
         contrib_energy = 0.0_DP
      endif
      if(pub_dmft_splitk)then
         dmft_splitk_iter = dmft_splitk_iter + 1
         call utils_abort("Reached point in code where nfs sync used to &
              &appear; the code needs to be updated to use these settings.")
         ! if(pub_on_root)call sync_via_nfs(contrib_energy, &
         !      'summing_contrib_energy' // &
         !      trim(adjustl(utils_int_to_str(dmft_splitk_iter))))
         ! call comms_barrier
      endif
      dmft_energy_cor_(is) = dmft_energy_cor_(is) + contrib_energy

      ! ebl: Last mu step reached, preparing kernel calculations
      if(nkpoints == 1)then
         call sparse_copy(denskern_tmp%kern%m(is, PUB_1K)%p, matrixDensKernel)
      else
         if(pub_dmft_nkpoints > 0 .and. .not.pub_dmft_kpoints_kernel_gamma)then
            write(stdout, *) 'ERROR : if you use symmetries (other than k &
                 &inversion)'
            write(stdout, *) ' the density kernel will be wrong, only fine to &
                 &compute projected quantities and total charge'
            call utils_abort('Error in dmft_clean_dmft_energy')
         else
            if(pub_dmft_kpoints_kernel_gamma)Then
               if(ikpoint <= nkpoints)then
                  if(internal_norm_vector(kpt(ikpoint, :)) > 1.0E-4_DP) call &
                       sparse_scale(matrixDensKernel, 0.0_DP) !Not Gamma point
               end if
               if(pub_dmft_splitk)then
                  dmft_splitk_iter = dmft_splitk_iter + 1
                  call utils_abort("Reached point in code where nfs sync used &
                       &to appear; the code needs to be updated to use these &
                       &settings.")
                  ! if(pub_on_root) call sync_via_nfs(matrixDensKernel, &
                  !      'sum_matrixDensKernel' // &
                  !      trim(adjustl(utils_int_to_str(dmft_splitk_iter))), &
                  !      sparse_index_length(matrixDensKernel))
                  call comms_barrier
                  call comms_bcast(0, matrixDensKernel%zmtx)
               endif
               call sparse_copy(denskern_tmpc%kern%m(is, PUB_1K)%p, &
                    matrixDensKernel)
               call sparse_axpy(denskern_tmp%kern%m(is, PUB_1K)%p, &
                    denskern_tmpc%kern%m(is, PUB_1K)%p, 1.0_DP)
            else
               if(ikpoint <= nkpoints)then
                  if(    abs(norm_kpoints(ikpoint)*totkpoint-1.0_DP) < 1.0E-4_DP)then
                     call sparse_scale(matrixDensKernel, norm_kpoints(ikpoint))
                  elseif(abs(norm_kpoints(ikpoint)*totkpoint-2.0_DP) < 1.0E-4_DP)then
                     call sparse_create(matrixDensKernel_t, matrixDensKernel)
                     call sparse_transpose(matrixDensKernel_t, matrixDensKernel)
                     call sparse_axpy(matrixDensKernel, matrixDensKernel_t, &
                          1.0_DP)
                     call sparse_scale(matrixDensKernel, &
                          norm_kpoints(ikpoint)/2.0_DP)
                     call sparse_destroy(matrixDensKernel_t)
                  elseif( norm_kpoints(ikpoint)*totkpoint > 2.001_DP .or. &
                       norm_kpoints(ikpoint)*totkpoint < 0.99_DP ) then
                     write(stdout, *) ' norm kpoint is : ', &
                          norm_kpoints(ikpoint)*totkpoint
                     call utils_abort('Error in clean_dmft_energy: error in &
                          &k-points')
                  endif
               endif
               if(pub_dmft_splitk)then
                  dmft_splitk_iter = dmft_splitk_iter + 1
                  call utils_abort("Reached point in code where nfs sync used &
                       &to appear; the code needs to be updated to use these &
                       &settings.")
                  ! if(pub_on_root) call sync_via_nfs(matrixDensKernel, &
                  !      'sum_matrixDensKernel' // &
                  !      trim(adjustl(utils_int_to_str(dmft_splitk_iter))), &
                  !      sparse_index_length(matrixDensKernel))
                  call comms_barrier
                  call comms_bcast(0, matrixDensKernel%zmtx)
               endif
               call sparse_copy(denskern_tmpc%kern%m(is, PUB_1K)%p, &
                    matrixDensKernel)
               call sparse_axpy(denskern_tmp%kern%m(is, PUB_1K)%p, &
                    denskern_tmpc%kern%m(is, PUB_1K)%p, 1.0_DP)
            endif
         endif
      endif

      if(ikpoint >= lastkpoint)then
         ttscreen = sparse_trace(denskern_tmp%kern%m(is, PUB_1K)%p, overlap)
         if(nkpoints > 1)then
            if(pub_on_root)write(stdout, *) 'TOTAL CHARGE ESTIMATED WITH &
                 &rho_loc * S : ', ttscreen
            if(pub_on_root)write(stdout, *) 'TOTAL CHARGE with Sum_k rho_k S_k &
                 &: ', vectorNmu(nmu_step)
         endif
         ! if(pub_on_root)write(stdout, *) 'Trace DensKern(is) Overlap : ', &
         !      ttscreen
         ! ttscreen = sparse_trace(denskern_tmp%kern%m(is, PUB_1K), ham(is))
         ! if(pub_on_root)write(stdout, *) 'Trace DensKern(is) HAM(is) : ', &
         !      ttscreen
         if(pub_num_spins == 1)then
            dmft_energy_cor_(1) = dmft_energy_cor_(1)*2.0_DP
            dmft_energy_cor_(2) = 0.0_DP
         endif
         if(pub_dmft_paramagnetic .and. pub_num_spins == 2) then
            if(pub_debug_on_root)write(stdout, *) 'DEBUG: Copying the spin 1 &
                 &density kernel to the spin 2 channel (paramagnetic &
                 &calculation)'
            call sparse_copy(denskern_tmp%kern%m(2, PUB_1K)%p, &
                 denskern_tmp%kern%m(1, PUB_1K)%p)
            dmft_energy_cor_(1) = dmft_energy_cor_(1)*2.0_DP
            dmft_energy_cor_(2) = 0.0_DP
         endif
      endif

   endif

   if(nmu_step < nmu .and. (is == num_spins .or. pub_dmft_impose_chem_spin) &
        .and. ikpoint >= lastkpoint) then

      if (mu_via_bisection) then
         ! ebl: simplified bisection method, currently hard-coded off
         !      (awaiting tweaking of hubbard_dmft_interface -- i.e. removal of nmu
         !       -- to turn it on)
         chem_shift = dmft_mu_updated_by_bisection(vectorMu(1:nmu_step), &
              vectorNmu(1:nmu_step), target_N) - vectorMu(nmu_step)
      else
         ! ebl: TODO
         ! The following two branches (if/if not pub_debug) are NOT equivalent.
         ! This should be
         ! investigated, but for now the debug branch will simply be removed...

         ! if (pub_debug) then
         !    if(nmu_step > 1)then
         !       if(abs(vectorMu(nmu_step)-vectorMu(nmu_step-1)) < 1.0E-8_DP)then
         !          slope_mu =  0.0_DP
         !       else
         ! slope_mu = (vectorNmu(nmu_step)-vectorNmu(nmu_step-1)) /
         ! (vectorMu(nmu_step) - vectorMu(nmu_step-1))
         !       endif
         !       if(pub_on_root) write(stdout, *)  'slope is : ', slope_mu
         !       if(abs(slope_mu) < 50.0_DP)then
         ! chem_shift = sign(1.0_DP, slope_mu) * ( target_N - vectorNmu(nmu_step) )
         ! * 0.01_DP
         !       else
         ! chem_shift = slope_mu_mixing * ( target_N - vectorNmu(nmu_step) ) /
         ! slope_mu
         ! if(chem_shift > slope_mu_max_step) chem_shift = slope_mu_max_step
         ! if(chem_shift < -slope_mu_max_step) chem_shift = -slope_mu_max_step
         !       endif
         !    endif

         !    if(nmu_step == 1)then
         !       if(vectorNmu(nmu_step) > target_N)then
         !          chem_shift = -0.003_DP
         !       else
         !          chem_shift = 0.003_DP
         !       endif
         !    endif
         ! else
         np = min(nmu_step, abs(pub_dmft_mu_order)) ! if min of (nmu_step, 2)
         ! recover method if debug mode
         if (pub_debug_on_root) then
            write(stdout, *) 'np = ', np
            write(stdout, *) 'vectorMu(nmu_step-np + 1:nmu_step) = ', &
                 vectorMu(nmu_step-np + 1:nmu_step)
            write(stdout, *) 'vectorNmu(nmu_step-np + 1:nmu_step) = ', &
                 vectorNmu(nmu_step-np + 1:nmu_step)
            write(stdout, *) 'slope_mu_max_step = ', slope_mu_max_step
            write(stdout, *) 'slope_mu_mixing = ', slope_mu_mixing
            write(stdout, *) 'target_N = ', target_N
            write(stdout, *) 'reshoot = ', pub_dmft_mu_order < 0
            write(stdout, *) 'chemshift = ', chem_shift
         end if
         call dmft_shoot_next_mu(np, vectorMu(nmu_step-np + 1:nmu_step), &
              vectorNmu(nmu_step-np + 1:nmu_step), slope_mu_max_step, &
              slope_mu_mixing, target_N, reshoot = pub_dmft_mu_order < 0, &
              chem_shift = chem_shift)
         if (pub_debug_on_root) then
            write(stdout, *) 'np = ', np
            write(stdout, *) 'vectorMu(nmu_step-np + 1:nmu_step) = ', &
                 vectorMu(nmu_step-np + 1:nmu_step)
            write(stdout, *) 'vectorNmu(nmu_step-np + 1:nmu_step) = ', &
                 vectorNmu(nmu_step-np + 1:nmu_step)
            write(stdout, *) 'slope_mu_max_step = ', slope_mu_max_step
            write(stdout, *) 'slope_mu_mixing = ', slope_mu_mixing
            write(stdout, *) 'target_N = ', target_N
            write(stdout, *) 'reshoot = ', pub_dmft_mu_order < 0
            write(stdout, *) 'chemshift = ', chem_shift
         end if
         ! end if

         if(pub_debug_on_root)then
            write(stdout, *) 'CHEM_SHIFT       : ', chem_shift
            write(stdout, *) 'NMU_STEP         : ', nmu_step
            write(stdout, *) 'TARGET_N         : ', target_N
            write(stdout, *) 'N_(mustep)       : ', vectorNmu(nmu_step)
            write(stdout, *) 'chem(mustep)     : ', vectorMu(nmu_step)
            if(nmu_step > 1)then
               write(stdout, *) 'N_(mustep-1)     : ', vectorNmu(nmu_step-1)
               write(stdout, *) 'chem(mustep-1)   : ', vectorMu(nmu_step-1)
            endif
         endif
      end if

   endif

   if(nmu > 1 .and. ikpoint >= lastkpoint .and. pub_dmft_write)then
      ! if(pub_on_root)then
      !    write(stdout, *) 'CHEM POT / CHARGE OBTAINED UP TO NOW ....'
      !    do i = 1, nmu_step
      !       write(stdout, *) vectorMu(i), vectorNmu(i)
      !    enddo
      !    write(stdout, *) '.................'
      ! endif
      if(pub_on_root .and. pub_dmft_mpi_rank == 1)then
         funit = utils_unit()
         open(unit = funit, file = 'chem.potential.nmu.iter' // &
              trim(adjustl(utils_int_to_str(is))))
         write(funit, *) vectorMu(nmu_step)
         close(funit)
      endif
   endif

   ! Deallocation of matrices
   deallocate(tmpr, stat = ierr)
   call utils_dealloc_check('dmft_clean_dmft_energy', 'tmpr', ierr)
   deallocate(tmpr2, stat = ierr)
   call utils_dealloc_check('dmft_clean_dmft_energy', 'tmpr2', ierr)
   deallocate(tmpc, stat = ierr)
   call utils_dealloc_check('dmft_clean_dmft_energy', 'tmpc', ierr)
   deallocate(tmpc2, stat = ierr)
   call utils_dealloc_check('dmft_clean_dmft_energy', 'tmpc2', ierr)
   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_clean_dmft_energy'

 end subroutine

 subroutine dmft_prepare_for_exit(ham, greenf, self_energy, self_infinity, &
      hub_overlap, ngwf_basis, kpoint_data, sc_data, energy, myfrequencies, &
      fermi_e, nkpoints, ham_embed, overlap_cplx, denskern_tmp, denskern_tmpc, &
      matrixDensKernel, matrixSigStat, matrixAsig, nabla, nabla_square, &
      orb_in_plane, atom_in_plane, split_hub, split_hubt, split_wgv, split_S, &
      split_Sm, split_H, inv_projectors)

   !---------------------------------------------------------------------------!
   ! Deallocates variables in preparation for exiting the program              !
   !---------------------------------------------------------------------------!

   use dense,          only: dense_destroy
   use utils,          only: utils_busy_wait, utils_unit, utils_dealloc_check
   use comms,          only: comms_barrier, pub_my_proc_id, pub_on_root
   use sparse,         only: sparse_destroy
   use sparse_embed,   only: sparse_embed_destroy
   use file_handling, only: file_handling_remove_file, &
        file_handling_remove_numbered_files
   use kernel,         only: DKERN, kernel_destroy
   use rundat, only: pub_debug, pub_debug_on_root, pub_dmft_complex_freq, &
        pub_dmft_kernel, pub_dmft_optics, pub_num_spins
   use function_basis, only: FUNC_BASIS

   implicit none

   ! Arguments
   type(SPAM3),                     intent(in   ) :: ham(pub_num_spins)
   type(DMFT_MATRIX),               intent(inout) :: greenf
   type(DMFT_MATRIX),               intent(inout) :: self_energy
   type(SPAM3),                     intent(inout) :: &
                                                     self_infinity(pub_num_spins)
   type(DMFT_OVERLAP),              intent(inout) :: hub_overlap
   type(FUNC_BASIS),                intent(in   ) :: ngwf_basis
   type(DMFT_KPOINT_MATRICES),      intent(inout) :: kpoint_data
   type(DMFT_SC_MATRICES),          intent(inout) :: sc_data
   integer, allocatable,            intent(inout) :: myfrequencies(:)
   complex(kind = DP),              intent(in   ) :: energy
   real(kind = DP),                 intent(in   ) :: fermi_e(2)
   integer,                         intent(in   ) :: nkpoints
   type(SPAM3),                     intent(inout) :: ham_embed(pub_num_spins)
   type(SPAM3),                     intent(inout) :: overlap_cplx
   type(DKERN),                     intent(inout) :: denskern_tmp
   type(DKERN),                     intent(inout) :: denskern_tmpc
   type(SPAM3),                     intent(inout) :: matrixDensKernel
   type(DEM), allocatable,          intent(inout) :: matrixSigStat(:)
   type(DEM), allocatable,          intent(inout) :: matrixAsig(:)
   type(SPAM3_EMBED),                     intent(inout) :: nabla(3)
   real(kind = DP), allocatable,    intent(inout) :: nabla_square(:,:,:)
   logical, allocatable,            intent(inout) :: orb_in_plane(:,:)
   logical, allocatable,            intent(inout) :: atom_in_plane(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_hub(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_hubt(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_wgv(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_S(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_Sm(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_H(:,:)
   logical,                         intent(in   ) :: inv_projectors

   ! Local variables
   integer :: i
   integer :: j
   integer :: is
   integer :: funit
   integer :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_prepare_for_exit'

   deallocate(myfrequencies, stat = ierr)
   call utils_dealloc_check('hubbard_dmft_interface', 'myfrequencies', ierr)

   if(pub_dmft_split .and. allocated(split_hub)) then
      deallocate(split_hub, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'split_hub', ierr)
      deallocate(split_hubt, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'split_hubt', ierr)
      deallocate(split_wgv, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'split_wgv', ierr)
      deallocate(split_S, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'split_S', ierr)
      deallocate(split_Sm, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'split_Sm', ierr)
      deallocate(split_H, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'split_H', ierr)
   endif

#ifdef LINEAR_SCALING_INV
   if(pub_dmft_lin_scaling)then
      call internal_invert_linear_scaling(greenf, energy, overlap_cplx, &
           ham(is), self_energy%mat, is, fermi_e, ngwf_basis%num, &
           self_energy%backup, matrixSigStat, update = .false., realorcomplexH &
           = .true., finalize = .true.)
   endif
#endif

   if(nkpoints > 1)then
      call sparse_destroy(kpoint_data%unrot_overlap)
      call sparse_destroy(kpoint_data%unrot_inv_overlap)
      call sparse_destroy(kpoint_data%overlap)
      call sparse_destroy(kpoint_data%ham)
      call sparse_destroy(kpoint_data%inv_overlap)
      call sparse_destroy(kpoint_data%self_energy)
      call sparse_destroy(kpoint_data%hub_overlap)
   endif

   if(allocated(orb_in_plane)) then
      deallocate(orb_in_plane, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'orb_in_plane', ierr)
   end if
   if(allocated(atom_in_plane)) then
      deallocate(atom_in_plane, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'atom_in_plane', ierr)
   end if

   call comms_barrier

   if (pub_debug .and. allocated(sc_data%eigenvals_backup)) then
      deallocate(sc_data%eigenvals_backup, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', &
           'sc_data%eigenvals_backup', ierr)
   end if

   if(allocated(sc_data%eigenvals)) then
      do is = 1, pub_num_spins
         call sparse_destroy(sc_data%ham(is))
         call dense_destroy(sc_data%eigenvecs(is))
         call dense_destroy(matrixSigStat(is))
         call dense_destroy(matrixAsig(is))
      end do
      call sparse_destroy(matrixDensKernel)
      call dense_destroy(sc_data%tail)
      call sparse_destroy(sc_data%overlap)
      deallocate(sc_data%eigenvals, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'sc_data%eigenvals', &
           ierr)
      deallocate(matrixSigStat, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'matrixSigStat', ierr)
      deallocate(greenf%backup, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'greenf%backup', ierr)
      deallocate(greenf%inv_backup, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'greenf%inv_backup', &
           ierr)
      deallocate(self_energy%backup, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'self_energy%backup', &
           ierr)
      deallocate(matrixAsig, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'matrixAsig', ierr)
      deallocate(sc_data%eigenvecs, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'sc_data%eigenvecs', &
           ierr)
      deallocate(sc_data%ham, stat = ierr)
      call utils_dealloc_check('dmft_prepare_for_exit', 'sc_data%ham', ierr)
   endif

   do is = 1, pub_num_spins
      call sparse_destroy( ham_embed(is) )
   enddo
   if (pub_dmft_kernel /= 0) then
      call kernel_destroy(denskern_tmp)
   end if
   do is = 1, pub_num_spins
      if(pub_dmft_kernel /= 0)then
         if(nkpoints > 1)then
            if(is == 1)then
               call kernel_destroy( denskern_tmpc )
            endif
         endif
      endif
   enddo
   if(.not. pub_dmft_complex_freq .and. pub_dmft_optics)then
      if(allocated(nabla_square)) then
         deallocate(nabla_square, stat=ierr)
         call utils_dealloc_check('dmft_prepare_for_exit', 'nabla_square', ierr)
      end if
      do is = 1, 3
         call sparse_embed_destroy(nabla(is))
      enddo
   endif

   do is = 1, pub_num_spins
      call sparse_destroy(self_infinity(is))
   enddo

   call sparse_destroy(self_energy%mat)
   ! call sparse_destroy(self_energy%backup_spam)
   call sparse_destroy(self_energy%mat_v)
   call sparse_destroy(self_energy%w_mat_v)

   if( inv_projectors )Then
      call sparse_destroy(greenf%w_inv_mat_v)
      call sparse_destroy(greenf%w_inv_mat)
      call sparse_destroy(hub_overlap%inv_tmatc)
      call sparse_destroy(hub_overlap%inv_tmat)
      call sparse_destroy(hub_overlap%inv_matc)
      call sparse_destroy(hub_overlap%inv_mat)
   endif

   call sparse_destroy(greenf%w_mat_v)
   call sparse_destroy(greenf%w_mat)
   call sparse_destroy(hub_overlap%tmatc)
   call sparse_destroy(hub_overlap%matc)
   call sparse_destroy(overlap_cplx)
   call sparse_destroy(greenf%inv_mat)
   call sparse_destroy(greenf%inv_backup_spam)
   call sparse_destroy(greenf%backup_spam)
   call sparse_destroy(greenf%mat)
   call comms_barrier
   ! funit = utils_unit()
   ! open(unit = funit, file = 'onetep_confirmation_' // trim(adjustl( &
   !    utils_int_to_str(pub_dmft_mpi_rank))), position = 'append')
   ! write(funit, *) 1111
   ! close(funit)

   if(pub_on_root .and. pub_dmft_mpi_size > 1 .and. pub_dmft_mpi_rank == 1)then
#ifdef DMFT_SAFE_NFS
      call utils_busy_wait(4000)
#endif
      if (pub_debug) then
         call file_handling_remove_numbered_files( &
              "onetep_dmft_confirmation_dmft_density_kernel_")
         call file_handling_remove_file( &
              "onetep_dmft_confirmation_dmft_density_count")
      end if
   end if

   ! #ifdef GPU_SPEEDUP
   ! if(pub_dmft_use_gpu_onlyme)then
   !    if(pub_output_detail >= VERBOSE) write(stdout, *) 'SHUTTING DOWN CUBLAS, &
   !         &proc/mpirank : ', pub_my_proc_id, pub_dmft_mpi_rank
   !    call cublas_shutdown()
   ! endif
   ! #endif
   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_prepare_for_exit'

 end subroutine

 subroutine dmft_eigenvec_fullmat(n_eig, aspam, bspam, eigenvecs, eigenvals)

   !---------------------------------------------------------------------------!
   ! Interface to dense_mod to solve for the eigenvectors and values           !
   ! corresponding to the generalised eigenproblem defined by the sparse       !
   ! matrices aspam and bspam                                                  !
   !---------------------------------------------------------------------------!

   use dense, only: DEM, dense_convert, dense_create, dense_destroy, &
        dense_eigensolve, dense_get_element
   use sparse, only: sparse_convert
   use rundat, only: pub_debug_on_root

   implicit none

   ! Arguments
   integer,         intent(in   ) :: n_eig
   type(SPAM3),     intent(in   ) :: aspam
   type(SPAM3),     intent(in   ) :: bspam
   type(DEM),       intent(inout) :: eigenvecs
   real(kind = DP), intent(inout) :: eigenvals(:)

   ! Local variables
   integer   :: i, j
   type(DEM) :: adens, bdens

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_eigenvec_fullmat'


   call dense_create(adens, n_eig, n_eig, iscmplx=eigenvecs%iscmplx)
   call dense_create(bdens, n_eig, n_eig, iscmplx=eigenvecs%iscmplx)
   call dense_convert(adens, aspam)
   call dense_convert(bdens, bspam)
   call dense_eigensolve(n_eig, eigenvals, adens, bdens, 1, eigenvecs)
   call dense_destroy(adens)
   call dense_destroy(bdens)
   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_eigenvec_fullmat'

 end subroutine

 subroutine internal_invert_full_real(mat_inv, mat, nrows, backupflag)

   !---------------------------------------------------------------------------!
   ! Performs the explicit inversion of a sparse (real) matrix)                !
   !---------------------------------------------------------------------------!

   use dense, only: DEM, dense_convert, dense_create, dense_destroy, &
        dense_invert
   use rundat, only: pub_debug_on_root

   implicit none

   ! Arguments
   type(SPAM3),       intent(in   ) :: mat
   type(SPAM3),       intent(inout) :: mat_inv
   integer,           intent(in   ) :: nrows ! Will always be ngwf_basis%num
   logical, optional, intent(inout) :: backupflag

   ! Local variables
   type(DEM) :: mat_dens

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &internal_invert_full_real'

   call dense_create(mat_dens, nrows, nrows, iscmplx = .false.)
   call dense_convert(mat_dens, mat)
   call dense_invert(mat_dens)
   call dense_convert(mat_inv, mat_dens)
   call dense_destroy(mat_dens)

   if(present(backupflag)) backupflag = .false.

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &internal_invert_full_real'

 end subroutine

 subroutine internal_invert_full(mat_inv, mat, nrows, backup, backupflag)

   !---------------------------------------------------------------------------!
   ! Performs the explicit inversion of a sparse matrix                        !
   !---------------------------------------------------------------------------!

   use dense, only: DEM, dense_convert, dense_create, dense_destroy, &
        dense_invert, dense_write
   use sparse, only: sparse_num_cols, sparse_num_rows
   use utils,  only: utils_assert
   use comms,  only: pub_total_num_procs
   use rundat, only: pub_debug, pub_debug_on_root

   implicit none

   ! Arguments
   type(SPAM3),                               intent(in   ) :: mat
   type(SPAM3),                               intent(inout) :: mat_inv
   integer,                                   intent(in   ) :: nrows ! Should
   complex(kind = DP), allocatable, optional, intent(inout) :: backup(:,:)
   logical, optional,                         intent(inout) :: backupflag

   ! Local variables
   type(DEM) :: mat_dens

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &internal_invert_full'

   call utils_assert(present(backup) .eqv. present(backupflag), "Error in &
        &internal_invert_full: one, but not both, of backup and backupflag &
        &provided")

   if(pub_total_num_procs == 1)then
      if (present(backup)) then
         call internal_invert_full_slow(mat_inv, mat, nrows, backup, &
              backupflag) !back to LAPACK
      else
         call internal_invert_full_slow(mat_inv, mat, nrows) !back to LAPACK
      end if
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &internal_invert_full'
      return
   endif
#ifdef INVERT_N_CUBE
   if (present(backup)) then
      ! back to LAPACK
      call internal_invert_full_slow(mat_inv, mat, nrows, backup, backupflag)
   else
      call internal_invert_full_slow(mat_inv, mat, nrows) !back to LAPACK
   end if
   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &internal_invert_full'
   return
#endif

   call utils_assert(sparse_num_rows(mat) == sparse_num_cols(mat), "Error in &
        &internal_invert_full: matrix is not square")

   call dense_create(mat_dens, nrows, nrows, iscmplx = .true.)
   call dense_convert(mat_dens, mat)
   call dense_invert(mat_dens)
   call dense_convert(mat_inv, mat_dens)
   call dense_destroy(mat_dens)

   if(present(backup)) backupflag = .false.

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &internal_invert_full'

 end subroutine

 subroutine internal_invert_linear_scaling(G, energy, ov, H, sig, is, fermi_e, &
      nngwf, sigmabackup, matrixSigStat, update, realorcomplexH, finalize)

   !---------------------------------------------------------------------------!
   ! Performs inversion of a sparse matrix in a linear-scaling fashion         !
   !---------------------------------------------------------------------------!
   ! NOTES                                                                     !
   !   1) this will NOT do the special point correctly, especially Simp Eimp   !
   !      will not be obtained. However, the high frequency limit will be      !
   !      fine, so Simp obtained from the tail will be correct and therefore   !
   !      eimp_from_onetep and simp_from_onetep need to be turned off in the   !
   !      dmft package                                                         !
   !   2) this won't work for nmu > 1 or the density kernel; use nmu = 1 and   !
   !      dmft_kernel = 0                                                      !
   !   3) the chemical potential should hence be the one obtained by the DFT - !
   !      not from the chem.potential file                                     !
   !---------------------------------------------------------------------------!


#ifdef LINEAR_SCALING_INV
   use inversion_with_window, only: invert_window
#endif
   use sparse,                only: sparse_convert
   use linalg,                only: linalg_invert_sym_cmatrix
   use utils,                 only: utils_assert, utils_abort
   use comms,                 only: pub_my_proc_id, pub_total_num_procs
   use rundat, only: pub_debug_on_root, &
        pub_dmft_nval, pub_dmft_sc, &
        pub_dmft_scaling_cutoff, &
        pub_dmft_scaling_meth, pub_dmft_scaling_nmpi, pub_dmft_scaling_tail, &
        pub_dmft_win

   implicit none

   ! Arguments
   type(DMFT_MATRIX),               intent(inout) :: G
   type(SPAM3),                     intent(in   ) :: ov
   type(SPAM3),                     intent(in   ) :: H
   type(SPAM3),                     intent(in   ) :: sig
   integer,                         intent(in   ) :: is
   real(kind = DP),                 intent(in   ) :: fermi_e(2)
   integer,                         intent(in   ) :: nngwf
   complex(kind = DP), allocatable, intent(inout) :: sigmabackup(:,:)
   type(DEM), allocatable,          intent(inout) :: matrixSigStat(:)
   logical,                         intent(in   ) :: update
   logical,                         intent(in   ) :: realorcomplexH
   logical, optional,               intent(in   ) :: finalize

   ! Local variables
   complex(kind = DP)              :: energy, ww
   real(kind = DP)                 :: chem
   integer                         :: lin_inv_nval, lin_inv_ns
   real(kind = DP)                 :: lin_inv_win
   real(kind = DP)                 :: lin_inv_tail, lin_inv_centre, cutoffsig, &
                                      cutoffH
   integer                         :: nmpi, nconv
   logical                         :: skip_invert_window
   integer                         :: lin_nconv_back
   real(kind = DP)                 :: lin_rconv_back
   complex(kind = DP), allocatable :: lin_ovSq(:,:)
   complex(kind = DP), allocatable :: lin_HrSq(:,:)
   complex(kind = DP), allocatable :: lin_HSq(:,:)
   complex(kind = DP), allocatable :: lin_sigSq(:,:)
   complex(kind = DP), allocatable :: lin_GSq(:,:)
   complex(kind = DP), allocatable :: lin_vecs(:,:)
   complex(kind = DP), allocatable :: lin_vec(:)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &internal_invert_linear_scaling'

#ifdef LINEAR_SCALING_INV
   call utils_assert("Error in &
        &internal_invert_linear_scaling: should not have reached here with &
        &dmft_lin_scaling : F")
   lin_inv_ns = 20
   lin_inv_nval = min(pub_dmft_nval, nngwf-1)
   lin_inv_win = pub_dmft_win

   ! Allocate matrices
   allocate(lin_HrSq(nngwf, nngwf), stat = ierr)
   call utils_alloc_check('internal_invert_linear_scaling', 'lin_HrSq', ierr)
   allocate(lin_ovSq(nngwf, nngwf), stat = ierr)
   call utils_alloc_check('internal_invert_linear_scaling', 'lin_ovSq', ierr)
   allocate(lin_HSq(nngwf, nngwf), stat = ierr)
   call utils_alloc_check('internal_invert_linear_scaling', 'lin_HSq', ierr)
   allocate(lin_sigSq(nngwf, nngwf), stat = ierr)
   call utils_alloc_check('internal_invert_linear_scaling', 'lin_sigSq', ierr)
   allocate(lin_GSq(nngwf, nngwf), stat = ierr)
   call utils_alloc_check('internal_invert_linear_scaling', 'lin_GSq', ierr)
   allocate(lin_vecs(nngwf, lin_inv_nval), stat = ierr)
   call utils_alloc_check('internal_invert_linear_scaling', 'lin_vecs', ierr)
   allocate(lin_vec(lin_inv_nval), stat = ierr)
   call utils_alloc_check('internal_invert_linear_scaling', 'lin_vec', ierr)

   if(.not.present(finalize))then
      if(update)then
         call sparse_convert(lin_ovSq, ov)
         if(.not.realorcomplexH)then
            call sparse_convert(lin_HSq, H)  ! Hk
         else
            call sparse_convert(lin_HrSq, H) ! H is real
            lin_HSq = lin_HrSq
         endif
      endif
      chem = fermi_e(is) + pub_dmft_chem_shift
      lin_inv_centre = chem
      lin_inv_tail = pub_dmft_scaling_tail
      ww = energy-chem
      nmpi = pub_dmft_scaling_nmpi
      cutoffH = pub_dmft_scaling_cutoff_h
      cutoffsig = pub_dmft_scaling_cutoff
      if(pub_dmft_split)then
         lin_sigSq = sigmabackup
      else
         call sparse_convert(lin_sigSq, sig)
      endif
      if(update)then
         lin_HSq   = lin_HSq   - real(matrixSigStat(is)%zmtx)
      else
         lin_sigSq = lin_sigSq - real(matrixSigStat(is)%zmtx)
      endif
   endif

   skip_invert_window = .false.
   if(pub_dmft_mpi_rank /= 1 .and. update .and. .not.pub_dmft_splitk)then
      if(.not.present(finalize)) then
         lin_vecs = 0.0_DP
         lin_vec = 0.0_DP
         lin_rconv_back = 0.0_DP
      endif
      skip_invert_window = .true.
   endif

   if (.not. skip_invert_window) then
      if(pub_dmft_mpi_size == 1 .or. pub_total_num_procs > 1 .or. &
           pub_dmft_splitk) nmpi = 1

      call invert_window(pub_my_proc_id, pub_total_num_procs, lin_HSq, ww, &
           lin_GSq, lin_inv_nval, nconv, lin_inv_ns, lin_inv_win, &
           lin_inv_centre, .false., update, lin_sigSq, cutoffsig, &
           pub_dmft_scaling_meth, lin_ovSq, lin_inv_tail, lin_inv_centre, &
           pub_dmft_mpi_rank, cutoffH, nmpi, pub_dmft_scaling_tol, &
           pub_dmft_scaling_maxspace, pub_dmft_scaling_iter, lin_nconv_back, &
           lin_vecs, lin_vec, finalize)
      lin_rconv_back = real(lin_nconv_back, kind=DP)
   end if

   if(pub_dmft_mpi_size > 1 .and. .not.present(finalize) .and. update .and. &
        .not.pub_dmft_splitk)then
      call utils_abort("Error in internal_invert_linear_scaling: depreciated syncing")
      ! call sync_via_nfs(lin_vecs, 'nfs_linear_inversionVECS', nngwf, &
      !      lin_inv_nval, longdelay = .true.)
      ! call sync_via_nfs(lin_vec, 'nfs_linear_inversionVEC', lin_inv_nval)
      ! call sync_via_nfs(lin_rconv_back, 'nfs_linear_inversionSCAL', longdelay &
      !      = .true.)
      ! lin_nconv_back = NINT(lin_rconv_back)
   endif

   if(.not.present(finalize)) then
      if(.not.pub_dmft_split)then
         call sparse_convert(G, lin_GSq)
      endif
      G%backup_flag = .true.
      G%backup = lin_GSq
   endif

   ! Deallocate matrices
   deallocate(lin_HrSq, stat = ierr)
   call utils_dealloc_check('internal_invert_linear_scaling', 'lin_HrSq', ierr)
   deallocate(lin_ovSq, stat = ierr)
   call utils_dealloc_check('internal_invert_linear_scaling', 'lin_ovSq', ierr)
   deallocate(lin_HSq, stat = ierr)
   call utils_dealloc_check('internal_invert_linear_scaling', 'lin_HSq', ierr)
   deallocate(lin_sigSq, stat = ierr)
   call utils_dealloc_check('internal_invert_linear_scaling', 'lin_sigSq', ierr)
   deallocate(lin_GSq, stat = ierr)
   call utils_dealloc_check('internal_invert_linear_scaling', 'lin_GSq', ierr)
   deallocate(lin_vecs, stat = ierr)
   call utils_dealloc_check('internal_invert_linear_scaling', 'lin_vecs', ierr)
   deallocate(lin_vec, stat = ierr)
   call utils_dealloc_check('internal_invert_linear_scaling', 'lin_vec', ierr)

#endif
   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &internal_invert_linear_scaling'

 end subroutine

 subroutine internal_invert_full_slow(mat_inv, mat, nrows, backup, backupflag)

   !---------------------------------------------------------------------------!
   ! Explicitly inverts a sparse matrix                                        !
   !---------------------------------------------------------------------------!

   use sparse, only: sparse_convert
   use linalg, only: linalg_invert_sym_cmatrix
   use utils,  only: utils_alloc_check, utils_assert, utils_dealloc_check
   use comms,  only: pub_on_root
   use rundat, only: pub_debug, pub_debug_on_root

   implicit none

   ! Arguments
   type(SPAM3),                               intent(in)    :: mat
   type(SPAM3),                               intent(inout) :: mat_inv
   integer,                                   intent(in   ) :: nrows
   complex(kind = DP), optional, allocatable, intent(inout) :: backup(:,:)
   logical, optional,                         intent(inout) :: backupflag

   ! Local variables
   complex(kind = DP), allocatable, dimension(:,:) :: mat_square
   integer                                         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &internal_invert_full_slow'

   call utils_assert(present(backup) .eqv. present(backupflag), "Error in &
        &internal_invert_full_slow: one, but not both, of backup and &
        &backupflag provided")

   allocate(mat_square(nrows, nrows), stat = ierr)
   call utils_alloc_check('internal_invert_full_slow', 'mat_square', ierr)

   call sparse_convert(mat_square, mat)

   call linalg_invert_sym_cmatrix(mat_square, nrows)

   ! CWTODO: I'm not entirely sure why this loop was needed, as there is no
   ! k-dependent behaviour?
   ! if(nkpoints == 1)then
   !    call linalg_invert_sym_cmatrix(mat_square, ngwf_basis%num)
   ! else
   !    call invert_gen_cmat(ngwf_basis%num, mat_square)
   ! endif

   call sparse_convert(mat_inv, mat_square)

   if (present(backup)) then
      backupflag = .true.
      backup = mat_square
   endif

   deallocate(mat_square, stat = ierr)
   call utils_dealloc_check('internal_invert_full_slow', 'mat_square', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &internal_invert_full_slow'

 end subroutine

 subroutine dmft_init_optics(nabla_square, nabla, overlap, ngwf_basis, rep, &
      mdl, nl_projectors, proj_basis)

   !---------------------------------------------------------------------------!
   ! Initialises variables associated with optics calculations                 !
   !---------------------------------------------------------------------------!

   use comms,               only: pub_on_root
   use function_basis,      only: FUNC_BASIS
   use integrals,           only: integrals_grad
   use model_type,          only: MODEL
   use ngwf_representation, only: NGWF_REP
   use optics,              only: optics_grad_mat_els
   use paw,                 only: paw_species_init_proj
   use projectors,          only: PROJECTOR_SET, projectors_create_real
   use pseudopotentials,    only: pseudo_species_init_proj
   use rundat, only: pub_any_nl_proj, pub_debug_on_root, pub_dmft_optics, &
        pub_dmft_optics_i1, pub_dmft_optics_i2, pub_dmft_read, pub_dmft_write, &
        pub_paw, pub_realspace_projectors, pub_rootname
   use sparse,              only: sparse_convert, sparse_create
   use utils,               only: utils_abort, utils_assert, utils_unit, &
                                  utils_alloc_check, utils_dealloc_check

   implicit none

   ! Arguments
   real(kind = DP), allocatable, intent(inout) :: nabla_square(:,:,:)
   type(SPAM3_EMBED),                  intent(inout) :: nabla(3)
   type(SPAM3),                  intent(in   ) :: overlap
   type(FUNC_BASIS),             intent(in   ) :: ngwf_basis(1)
   type(NGWF_REP),               intent(in   ) :: rep
   type(MODEL),                  intent(in   ) :: mdl
   type(PROJECTOR_SET),          intent(inout) :: nl_projectors(1)
   type(FUNC_BASIS), optional,   intent(in   ) :: proj_basis(1)

   ! Local variables
   integer         :: i
   logical         :: check_store_nabla
   real(kind = DP) :: norm_optic
   integer         :: free_unit
   integer         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_init_optics'

   if(pub_dmft_optics)then
      write(stdout, *) '--initialize optical conductivity calculations--'
      if (allocated(nabla_square)) then
         deallocate(nabla_square, stat=ierr)
         call utils_dealloc_check('dmft_init_optics', 'nabla_square', ierr)
      end if
      allocate(nabla_square(ngwf_basis(1)%num, ngwf_basis(1)%num, 3), stat = ierr)
      call utils_alloc_check('dmft_init_optics', 'nabla_square', ierr)
      do i = 1, 3
         call sparse_create(nabla(i)%p, overlap, iscmplx = .false.) !BUG corrected
         !call sparse_create(nabla(i), greenf, iscmplx = .false.)
      enddo
      write(stdout, *) 'build < phi_i | d/dx phi_j > '
      inquire(file = trim(pub_rootname) // '.nabla', exist = check_store_nabla)
      !---------------------------------------------------------------!
      if(check_store_nabla .and. pub_dmft_read)then
         free_unit = utils_unit()
         open(unit = free_unit, file = trim(pub_rootname) // '.nabla', form = &
              'unformatted')
         read(free_unit) nabla_square
         close(free_unit)
      else
         ! BUG (?) :
         ! call integrals_grad( nabla, rep%ngwfs_on_grid, ngwf_basis,
         ! rep%ngwfs_on_grid, ngwf_basis )
         call utils_assert(present(proj_basis), "Error in dmft_init_optics: &
              &NLPP argument PROJ_BASIS is not present")

         ! ! ebl
         ! ! Dirty fix to allow access nl_projectors. Will likely need to
         ! place this routine
         ! ! elsewhere (rather than in dmft_mod) so that access is more natural
         ! if (pub_any_nl_proj) then
         ! call pseudo_species_init_proj(nl_projectors, mdl%pseudo_sp,
         ! mdl%elements)
         ! end if
         ! if (pub_paw) then
         ! call paw_species_init_proj(nl_projectors, mdl%paw_sp, mdl%elements)
         ! end if

         ! if (pub_any_nl_proj .and. pub_realspace_projectors) call
         ! projectors_create_real(proj_basis, nl_projectors, mdl%fftbox,
         ! mdl%cell)
         ! if (pub_paw .and. pub_realspace_projectors) call
         ! projectors_create_real(proj_basis, nl_projectors, mdl%fftbox,
         ! mdl%cell)

         call optics_grad_mat_els(nabla, rep, ngwf_basis, proj_basis, &
              nl_projectors, mdl)
         do i = 1, 3
            call sparse_convert( nabla_square(:, :, i), nabla(i)%p )
         enddo
         if(pub_on_root .and. pub_dmft_mpi_rank == 1 .and. pub_dmft_write)then
            free_unit = utils_unit()
            open(unit = free_unit, file = trim(pub_rootname) // '.nabla', &
                 form = 'unformatted')
            write(free_unit) nabla_square
            close(free_unit)
            write(stdout, *) trim(pub_rootname) // '.nabla file written'
         endif
      endif
      if(abs(pub_dmft_optics_x1) > 1.0E-5_DP &
           .or. abs(pub_dmft_optics_y1) > 1.0E-5_DP &
           .or. abs(pub_dmft_optics_z1) > 1.0E-5_DP)then
         norm_optic = sqrt(pub_dmft_optics_x1**2 + pub_dmft_optics_y1**2 + &
              pub_dmft_optics_z1**2)
         pub_dmft_optics_x1 = pub_dmft_optics_x1/norm_optic
         pub_dmft_optics_y1 = pub_dmft_optics_y1/norm_optic
         pub_dmft_optics_z1 = pub_dmft_optics_z1/norm_optic
         nabla_square(:, :, 1) = pub_dmft_optics_x1*nabla_square(:, :, 1) + &
              pub_dmft_optics_y1 *nabla_square(:, :, 2) + &
              pub_dmft_optics_z1*nabla_square(:, :, 3)
         if(abs(pub_dmft_optics_x1**2 + pub_dmft_optics_y1**2 + &
              pub_dmft_optics_z1**2-1.0_DP) > 1.0E-4_DP)then
            call utils_abort("Error in dmft_init_optics: dmft_optics_x1, y1, &
                 &z1 should collectively be a unitary vector")
         endif
         if(pub_dmft_optics_i1 /= 1 .or. pub_dmft_optics_i2 /= 1)then
            write(stdout, *) 'ERROR optical conductivity as defined along a &
                 &vector only works for current correlator along same &
                 &direction'
            write(stdout, *) ' please use pub_dmft_optics_i1 = &
                 &pub_dmft_optics_i2 = 1'
            call utils_abort("Error in DMFT optics calculation")
         endif
      endif
      !---------------------------------------------------------------!
      write(stdout, '(a, 200f8.4)') 'done, testing : nabla_square(1:2, 1:2) is &
           &', nabla_square(1:2, 1:2, 1)
   endif
   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving dmft_init_optics'

 end subroutine

 subroutine dmft_load_k_ham_and_overlap(ham, greenf, self_energy, overlap, &
      inv_overlap, overlap_cplx, h_atoms, hub, hub_overlap, hub_proj_basis, &
      ngwf_basis, kpoint_data, energy, fermi_e, is, ien, ikpoint, ffirst, &
      llast, kp_file_shift, hub_atom, hub_atomb, nkpoints, matrixSigStat, &
      rot_vec_angles, connection_table, dimer_table, nmerge, merge_table, &
      merge_neigh, ham_embed, check_embed, check_embed_h, split_hub, &
      split_hubt, split_wgv, split_S, split_Sm, split_H, cluster, &
      force_reopen_file, inv_projectors, one_k_one_file, &
      reopen_sigma_file_each_step, same_self_for_all, show_matrices, &
      scratchdir)

   !---------------------------------------------------------------------------!
   ! Defines the Hamiltonian and overlap matrices for multiple k-point         !
   ! calculations                                                              !
   !---------------------------------------------------------------------------!

   use comms,             only: pub_on_root
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par
   use rundat, only: pub_debug, pub_debug_on_root, pub_dmft_points, &
        pub_dmft_read, pub_num_spins, pub_rootname
   use sparse, only: sparse_convert, sparse_copy, sparse_num_cols, &
        sparse_num_rows, sparse_read
   use utils,             only: utils_assert, utils_int_to_str, &
                                utils_alloc_check, utils_dealloc_check

   implicit none

   ! Arguments
   type(SPAM3),                     intent(inout) :: ham(pub_num_spins)
   type(DMFT_MATRIX),               intent(inout) :: greenf
   type(DMFT_MATRIX),               intent(inout) :: self_energy
   type(SPAM3),                     intent(in   ) :: overlap
   type(SPAM3),                     intent(in   ) :: inv_overlap
   type(SPAM3),                     intent(inout) :: overlap_cplx
   type(ARRAY_OF_MATRICES),         intent(in   ) :: h_atoms(par%nat_hub)
   type(HUBBARD_MODEL),             intent(in   ) :: hub
   type(DMFT_OVERLAP),              intent(inout) :: hub_overlap
   type(FUNC_BASIS),                intent(in   ) :: hub_proj_basis
   type(FUNC_BASIS),                intent(in   ) :: ngwf_basis
   type(DMFT_KPOINT_MATRICES),      intent(inout) :: kpoint_data
   complex(kind = DP),              intent(in   ) :: energy
   real(kind = DP),                 intent(in   ) :: fermi_e(2)
   integer,                         intent(in   ) :: is
   integer,                         intent(inout) :: ien
   integer,                         intent(in   ) :: ikpoint
   integer,                         intent(in   ) :: ffirst
   integer,                         intent(in   ) :: llast
   integer,                         intent(inout) :: kp_file_shift
   integer,                         intent(inout) :: hub_atom
   integer,                         intent(inout) :: hub_atomb
   integer,                         intent(in   ) :: nkpoints
   type(DEM), allocatable,          intent(inout) :: matrixSigStat(:)
   real(kind = DP), allocatable,    intent(in   ) :: rot_vec_angles(:,:,:)
   logical, allocatable,            intent(inout) :: connection_table(:,:)
   integer, allocatable,            intent(inout) :: dimer_table(:,:)
   integer,                         intent(in   ) :: nmerge
   integer, allocatable,            intent(inout) :: merge_table(:,:)
   integer, allocatable,            intent(in   ) :: merge_neigh(:)
   type(SPAM3),                     intent(inout) :: ham_embed(pub_num_spins)
   logical,                         intent(inout) :: check_embed
   logical,                         intent(inout) :: check_embed_h
   complex(kind = DP), allocatable, intent(inout) :: split_hub(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_hubt(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_wgv(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_S(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_Sm(:,:)
   complex(kind = DP), allocatable, intent(inout) :: split_H(:,:)
   logical,                         intent(in   ) :: cluster
   logical,                         intent(inout) :: force_reopen_file
   logical,                         intent(in   ) :: inv_projectors
   logical,                         intent(in   ) :: one_k_one_file
   logical,                         intent(in   ) :: reopen_sigma_file_each_step
   logical,                         intent(inout) :: same_self_for_all
   logical,                         intent(in   ) :: show_matrices
   character(2000),                 intent(in   ) :: scratchdir

   ! Local variables
   logical                                         :: check_k
   integer                                         :: i, iproj, ingwf, i1, i2, &
                                                      kk, kien
   logical                                         :: verbose
   complex(kind = DP), allocatable, dimension(:,:) :: mat1, mat2, mat3, mat4, &
                                                      mat5, mat6
   real(kind = DP), allocatable                    :: split_Hr(:,:)
   real(kind = DP), allocatable                    :: split_Sr(:,:)
   character(20)                                   :: filek
   integer                                         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_load_k_ham_and_overlap'

   if(pub_dmft_split)then
      if(.not.allocated(split_hub))then
         allocate(split_hub(sparse_num_rows(hub_overlap%matc), &
              sparse_num_cols(hub_overlap%matc)), stat = ierr)
         call utils_alloc_check('dmft_load_k_ham_and_overlap', 'split_hub', ierr)
         allocate(split_hubt(sparse_num_rows(hub_overlap%tmatc), &
              sparse_num_cols(hub_overlap%tmatc)), stat = ierr)
         call utils_alloc_check('dmft_load_k_ham_and_overlap', 'split_hubt', ierr)
         allocate(split_wgv(sparse_num_rows(greenf%w_mat_v), &
              sparse_num_cols(greenf%w_mat_v)), stat = ierr)
         call utils_alloc_check('dmft_load_k_ham_and_overlap', 'split_wgv', ierr)
         allocate(split_S(sparse_num_rows(overlap_cplx), &
              sparse_num_cols(overlap_cplx)), stat = ierr)
         call utils_alloc_check('dmft_load_k_ham_and_overlap', 'split_S', ierr)
         allocate(split_Sm(sparse_num_rows(overlap_cplx), &
              sparse_num_cols(overlap_cplx)), stat = ierr)
         call utils_alloc_check('dmft_load_k_ham_and_overlap', 'split_Sm', ierr)
         allocate(split_H(sparse_num_rows(overlap_cplx), &
              sparse_num_cols(overlap_cplx)), stat = ierr)
         call utils_alloc_check('dmft_load_k_ham_and_overlap', 'split_H', ierr)
      endif
   endif

   if(nkpoints == 1) then
      if(pub_dmft_lin_scaling) then
         call dmft_build_self_infinity(ham, h_atoms, hub, hub_overlap, &
              self_energy, ngwf_basis, hub_proj_basis, fermi_e, hub_atom, &
              hub_atomb, is, nkpoints, cluster, connection_table, dimer_table, &
              nmerge, merge_table, merge_neigh, ffirst, llast, kpoint_data, &
              rot_vec_angles, ham_embed, check_embed, check_embed_h, &
              matrixSigStat, force_reopen_file, same_self_for_all, &
              show_matrices, reopen_sigma_file_each_step, scratchdir)
         call internal_invert_linear_scaling(greenf, energy, overlap_cplx, &
              ham(is), self_energy%mat, is, fermi_e, ngwf_basis%num, &
              self_energy%backup, matrixSigStat, update = .true., &
              realorcomplexH = .true.)
      endif
      kp_file_shift = 0
      if(pub_dmft_split)then
         call sparse_convert(split_hubt, hub_overlap%tmatc)
         call sparse_convert(split_hub, hub_overlap%matc)
         allocate(split_Sr(sparse_num_rows(overlap_cplx), &
              sparse_num_cols(overlap_cplx)), stat = ierr)
         call utils_alloc_check('dmft_load_k_ham_and_overlap', 'split_Sr', ierr)
         call sparse_convert(split_Sr, inv_overlap)
         split_Sm = split_Sr
         call sparse_convert(split_Sr, overlap)
         split_S = split_Sr
         deallocate(split_Sr, stat = ierr)
         call utils_dealloc_check('dmft_load_k_ham_and_overlap', 'split_Sr', ierr)
         allocate(split_Hr(sparse_num_rows(ham(is)), &
              sparse_num_cols(ham(is))), stat = ierr)
         call utils_alloc_check('dmft_load_k_ham_and_overlap', 'split_Hr', ierr)
         call sparse_convert(split_Hr, ham(is))
         split_H = split_Hr
         deallocate(split_Hr, stat = ierr)
         call utils_dealloc_check('dmft_load_k_ham_and_overlap', 'split_Hr', ierr)
      endif
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_load_k_ham_and_overlap'
      return
   endif

   verbose = pub_on_root .and. pub_dmft_mpi_rank == 1

   if(.not.pub_dmft_splitk)then
      if(one_k_one_file)then
         kp_file_shift = pub_dmft_mpi_size*(ikpoint-1)
      else
         kp_file_shift = 0
      endif
   else
      kp_file_shift = 0
   endif

   ! ebl: the following code depends on reading in
   call utils_assert(pub_dmft_read, "Error in dmft_load_k_ham_and_overlap: &
        &dmft_read set to false")

   filek = trim(pub_rootname) // ".hamk_" // trim(adjustl(utils_int_to_str( &
        ikpoint ))) // "_" // trim(adjustl(utils_int_to_str( is )))
   inquire(file = trim(adjustl(filek)), exist = check_k)
   call utils_assert(check_k, "Error in dmft_load_k_ham_and_overlap: " // &
        trim(adjustl(filek)) // " not found")
   if(verbose) write(stdout, *) 'Reading k-dependent Hamiltonian from file '//&
        trim(adjustl(filek))
   if (pub_dmft_mpi_rank == 1) call sparse_read(kpoint_data%ham, &
        trim(adjustl(filek)))
   !-------------!
   filek = trim(pub_rootname) // ".ovk_" // trim(adjustl(utils_int_to_str( &
        ikpoint ))) // "_" // trim(adjustl(utils_int_to_str( is )))
   inquire(file = trim(adjustl(filek)), exist = check_k)
   call utils_assert(check_k, "Error in dmft_load_k_ham_and_overlap: " // &
        trim(adjustl(filek)) // " not found")

   if(verbose) write(stdout, *) 'Reading k-dependent overlap from file '//&
        trim(adjustl(filek))
   if (pub_dmft_mpi_rank == 1) call sparse_read(kpoint_data%overlap, &
        trim(adjustl(filek)))
   call sparse_copy(overlap_cplx, kpoint_data%overlap)
   call dmft_invert_sparse_complex(kpoint_data%overlap, &
        kpoint_data%inv_overlap, 1)

   if(pub_on_root .and. pub_dmft_mpi_rank == 1) write(stdout, *) 'Now update &
        &the projectors'

   filek = trim(pub_rootname) // ".hubk_" // trim(adjustl( &
        utils_int_to_str( ikpoint )))
   inquire(file = trim(adjustl(filek)), exist = check_k)
   call utils_assert(check_k, "Error in dmft_load_k_ham_and_overlap: " // &
        trim(adjustl(filek)) // " not found")
   if(verbose) write(stdout, *) 'Reading k-dependent Hubbard overlap from &
        &file '//trim(adjustl(filek))
   if (pub_dmft_mpi_rank == 1) call sparse_read(hub_overlap%matc, &
        trim(adjustl(filek)))

   filek = trim(pub_rootname) // ".hubTk_" // trim(adjustl( &
        utils_int_to_str( ikpoint )))
   inquire(file = trim(adjustl(filek)), exist = check_k)
   call utils_assert(check_k, "Error in dmft_load_k_ham_and_overlap: " // &
        trim(adjustl(filek)) // " not found")
   if(verbose) write(stdout, *) 'Reading k-dependent Hubbard overlap transpose &
        &from file '//trim(adjustl(filek))
   if (pub_dmft_mpi_rank == 1) call sparse_read(hub_overlap%tmatc, &
        trim(adjustl(filek)))

   if(verbose) write(stdout, *) 'Computing P x {S_k^-1} x P'
   ! take care : hub_overlap and unrot_inv_overlap need to be in the unrotated
   ! basis
   !    Pk Ok Obar^-1 Simp Obar^-1 Ok Pk^\dagger
   !      !                          !
   !       rotation happens here     and here .........
   !   so kpoint_data%unrot_overlap should be in the unrotated basis
   ! - rotation is NOT done in project_green_function (but in
   ! build_green_function)
   ! - Obar was rotated when collected in files green_output, but is rotated
   ! back from sigma_output

   if(pub_debug_on_root)then
      write(stdout, *) 'Checking construction of matrix ovk: '
      call dmft_print_spam3_info(kpoint_data%overlap)
      write(stdout, *) 'Checking construction of matrix invk: '
      call dmft_print_spam3_info(kpoint_data%inv_overlap)
      write(stdout, *) 'Checking construction of matrix tmatrix cplx: '
      call dmft_print_spam3_info(hub_overlap%tmatc)
      write(stdout, *) 'Checking construction of matrix matrix cplx: '
      call dmft_print_spam3_info(hub_overlap%matc)
   endif

   call dmft_project_spam(kpoint_data%unrot_overlap, kpoint_data%inv_overlap, &
        hub_overlap, split_hub, split_hubt, split_wgv, greenf%backup, nosplit &
        = .true.)

   if(pub_debug_on_root)then
      write(stdout, *) 'Checking construction of unrotated overlap matrix: '
      call dmft_print_spam3_info(kpoint_data%unrot_overlap)
   endif

   if(verbose) write(stdout, *) 'Computing inverse of projection'
   call dmft_invert_sparse_complex(kpoint_data%unrot_overlap, &
        kpoint_data%unrot_inv_overlap, 2)

   if(pub_debug_on_root)then
      write(stdout, *) 'Checking construction of matrix invok: '
      call dmft_print_spam3_info(kpoint_data%unrot_inv_overlap)
   endif

   force_reopen_file = .true.
   kien = pub_dmft_points-6
   ien = kien
   call dmft_construct_self_energy(ien, kien, is, self_energy, hub_atom, &
        hub_atomb, connection_table, nmerge, merge_table, merge_neigh, &
        dimer_table, h_atoms, hub, hub_proj_basis, rot_vec_angles, ffirst, &
        llast, cluster, show_matrices, reopen_sigma_file_each_step, &
        force_reopen_file, scratchdir)
   call sparse_copy(kpoint_data%hub_overlap, self_energy%w_mat_v)
   force_reopen_file = .false.
   if(pub_debug_on_root)then
      write(stdout, *) 'Checking construction of matrix kpoint hub_overlap: '
      call dmft_print_spam3_info(kpoint_data%hub_overlap)
   endif

   if(pub_dmft_lin_scaling) then
      call dmft_build_self_infinity(ham, h_atoms, hub, hub_overlap, &
           self_energy, ngwf_basis, hub_proj_basis, fermi_e, hub_atom, &
           hub_atomb, is, nkpoints, cluster, connection_table, dimer_table, &
           nmerge, merge_table, merge_neigh, ffirst, llast, kpoint_data, &
           rot_vec_angles, ham_embed, check_embed, check_embed_h, &
           matrixSigStat, force_reopen_file, same_self_for_all, show_matrices, &
           reopen_sigma_file_each_step, scratchdir)
      call internal_invert_linear_scaling(greenf, energy, overlap_cplx, &
           kpoint_data%ham, kpoint_data%self_energy, is, fermi_e, &
           ngwf_basis%num, self_energy%backup, matrixSigStat, update = .true., &
           realorcomplexH = .false.)
   endif

   call dmft_update_inverse_projectors(hub_overlap, kpoint_data%inv_overlap, &
        inv_projectors)

   if(pub_dmft_split)then
      call sparse_convert(split_hubt, hub_overlap%tmatc)
      call sparse_convert(split_hub, hub_overlap%matc)
      call sparse_convert(split_Sm, kpoint_data%inv_overlap)
      call sparse_convert(split_S, kpoint_data%overlap)
      call sparse_convert(split_H, kpoint_data%ham)
   endif

   if (pub_debug) then
      write(stdout, *) 'DEBUG: running some matrix inversion tests'
      iproj  = sparse_num_rows(hub_overlap%tmatc)
      ingwf  = sparse_num_cols(hub_overlap%tmatc)
      allocate(mat1(iproj, ingwf), stat = ierr)
      call utils_alloc_check('dmft_load_k_ham_and_overlap', 'mat1', ierr)
      allocate(mat2(ingwf, iproj), stat = ierr)
      call utils_alloc_check('dmft_load_k_ham_and_overlap', 'mat2', ierr)
      allocate(mat5(ingwf, ingwf), stat = ierr)
      call utils_alloc_check('dmft_load_k_ham_and_overlap', 'mat5', ierr)
      allocate(mat6(iproj, iproj), stat = ierr)
      call utils_alloc_check('dmft_load_k_ham_and_overlap', 'mat6', ierr)
      call sparse_convert(mat1, hub_overlap%tmatc)
      call sparse_convert(mat2, hub_overlap%matc)
      call sparse_convert(mat5, overlap_cplx)
      allocate(mat3(size(mat1, 1), size(mat2, 2)), stat = ierr)
      call utils_alloc_check('dmft_load_k_ham_and_overlap', 'mat3', ierr)
      allocate(mat4(size(mat2, 1), size(mat1, 2)), stat = ierr)
      call utils_alloc_check('dmft_load_k_ham_and_overlap', 'mat4', ierr)
      mat3 = MATMUL(mat1, mat2)
      mat4 = MATMUL(mat2, mat1)

      !mat1 = Hub^T
      !mat2 = Hub
      !proj space
      write(300 + ikpoint, *) 'Hub^T Hub'
      do i = 1, size(mat3, 1)
         write(300 + ikpoint, *) 'i -', mat3(i, i)
      enddo
      !NGWFs space
      write(200 + ikpoint, *) 'Hub Hub^T'
      do i = 1, size(mat4, 1)
         write(200 + ikpoint, *) 'i -', mat4(i, i)
      enddo

      call sparse_convert(mat1, hub_overlap%tmatc)
      call sparse_convert(mat4, overlap_cplx)
      call sparse_convert(mat5, kpoint_data%inv_overlap)
      call sparse_convert(mat6, kpoint_data%unrot_inv_overlap)

      write(800 + ikpoint, *) 'MAX-MIN S_k*S_k^-1 : ', minval(abs(MATMUL(mat4, &
           mat5))), maxval(abs(MATMUL(mat4, mat5)))

      mat1 = MATMUL(mat1, mat5)
      mat3 = MATMUL(mat1, mat2)
      write(800 + ikpoint, *) 'Hub^T S^-1 Hub = Id_proj'
      do i = 1, size(mat3, 1)
         write(800 + ikpoint, *) 'i -', mat3(i, i)
      enddo

      mat3 = MATMUL(mat6, mat3)
      write(900 + ikpoint, *) 'Hub^T S^-1 Hub * unrot_inv_overlap'
      do i = 1, size(mat3, 1)
         write(900 + ikpoint, *) 'i -', mat3(i, i)
      enddo
      deallocate(mat1, stat = ierr)
      call utils_dealloc_check('dmft_load_k_ham_and_overlap', 'mat1', ierr)
      deallocate(mat2, stat = ierr)
      call utils_dealloc_check('dmft_load_k_ham_and_overlap', 'mat2', ierr)
      deallocate(mat3, stat = ierr)
      call utils_dealloc_check('dmft_load_k_ham_and_overlap', 'mat3', ierr)
      deallocate(mat4, stat = ierr)
      call utils_dealloc_check('dmft_load_k_ham_and_overlap', 'mat4', ierr)
      deallocate(mat5, stat = ierr)
      call utils_dealloc_check('dmft_load_k_ham_and_overlap', 'mat5', ierr)
      deallocate(mat6, stat = ierr)
      call utils_dealloc_check('dmft_load_k_ham_and_overlap', 'mat6', ierr)
   end if

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_load_k_ham_and_overlap'

 end subroutine

 subroutine dmft_update_inverse_projectors(hub_overlap, inv_overlap_kp, &
      inv_projectors)

   !---------------------------------------------------------------------------!
   ! Description:                                                              !
   ! NOTES                                                                     !
   !   For debugging purposes only. As inv_projectors is otherwise set to      !
   !   false, this might be able to be removed                                 !
   !---------------------------------------------------------------------------!

   use sparse, only: sparse_product
   use rundat, only: pub_debug_on_root

   implicit none

   ! Arguments
   type(DMFT_OVERLAP), intent(inout) :: hub_overlap
   type(SPAM3),        intent(inout) :: inv_overlap_kp
   logical,            intent(in   ) :: inv_projectors

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_update_inverse_projectors'

   if(.not.inv_projectors) return

   ! #ifndef GPU_SPEEDUP_MATMUL_R
   call sparse_product(hub_overlap%inv_matc, inv_overlap_kp, hub_overlap%matc)
   ! #else
   ! if(pub_dmft_use_gpu_onlyme)then
   !    call dmft_matmul_in_dense_format_r(hub_overlap%inv_matc, inv_overlap_kp, &
   !         hub_overlap%matc)
   ! else
   !    call sparse_product(hub_overlap%inv_matc, inv_overlap_kp, &
   !         hub_overlap%matc)
   ! endif
   ! #endif


   ! #ifndef GPU_SPEEDUP_MATMUL_R
   call sparse_product(hub_overlap%inv_tmatc, hub_overlap%tmatc, inv_overlap_kp)
   ! #else
   ! if(pub_dmft_use_gpu_onlyme)then
   !    call dmft_matmul_in_dense_format_r(hub_overlap%inv_tmatc, &
   !         hub_overlap%tmatc, inv_overlap_kp)
   ! else
   !    call sparse_product(hub_overlap%inv_tmatc, hub_overlap%tmatc, &
   !         inv_overlap_kp)
   ! endif
   ! #endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_update_inverse_projectors'

 end subroutine

 subroutine dmft_build_self_infinity(ham, h_atoms, hub, hub_overlap, &
      self_energy, ngwf_basis, hub_proj_basis, fermi_e, hub_atom, hub_atomb, &
      is, nkpoints, cluster, connection_table, dimer_table, nmerge, &
      merge_table, merge_neigh, ffirst, llast, kpoint_data, rot_vec_angles, &
      ham_embed, check_embed, check_embed_h, matrixSigStat, force_reopen_file, &
      same_self_for_all, show_matrices, reopen_sigma_file_each_step, &
      scratchdir)

   !---------------------------------------------------------------------------!
   ! Extracts the self energy at the largest frequency considered (stored in   !
   ! the (nfreq - 1)th entry)                                                  !
   !---------------------------------------------------------------------------!

   use dense, only: dense_axpy, dense_convert, dense_create, dense_destroy
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par
   use rundat, only: pub_debug_on_root, pub_dmft_points, pub_num_spins

   implicit none

   ! Arguments
   type(SPAM3),                  intent(inout) :: ham(pub_num_spins)
   type(ARRAY_OF_MATRICES),      intent(in   ) :: h_atoms(par%nat_hub)
   type(HUBBARD_MODEL),          intent(in   ) :: hub
   type(DMFT_OVERLAP),           intent(in   ) :: hub_overlap
   type(DMFT_MATRIX),            intent(inout) :: self_energy
   type(FUNC_BASIS),             intent(in   ) :: ngwf_basis
   type(FUNC_BASIS),             intent(in   ) :: hub_proj_basis
   real(kind = DP),              intent(in   ) :: fermi_e(2)
   integer,                      intent(inout) :: hub_atom
   integer,                      intent(inout) :: hub_atomb
   integer,                      intent(in   ) :: is
   integer,                      intent(in   ) :: nkpoints
   logical,                      intent(in   ) :: cluster
   logical, allocatable,         intent(inout) :: connection_table(:,:)
   integer, allocatable,         intent(inout) :: dimer_table(:,:)
   integer,                      intent(in   ) :: nmerge
   integer, allocatable,         intent(inout) :: merge_table(:,:)
   integer, allocatable,         intent(in   ) :: merge_neigh(:)
   integer,                      intent(in   ) :: ffirst
   integer,                      intent(in   ) :: llast
   type(DMFT_KPOINT_MATRICES),   intent(inout) :: kpoint_data
   real(kind = DP), allocatable, intent(in   ) :: rot_vec_angles(:,:,:)
   type(SPAM3),                  intent(inout) :: ham_embed(pub_num_spins)
   logical,                      intent(in   ) :: check_embed
   logical,                      intent(  out) :: check_embed_h
   type(DEM), allocatable,       intent(inout) :: matrixSigStat(:)
   logical,                      intent(  out) :: force_reopen_file
   logical,                      intent(  out) :: same_self_for_all
   logical,                      intent(in   ) :: show_matrices
   logical,                      intent(in   ) :: reopen_sigma_file_each_step
   character(2000),              intent(in   ) :: scratchdir

   ! Local variables
   integer   :: kien, ien
   type(DEM) :: tmp1

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_build_self_infinity'

   force_reopen_file = .true.
   kien = pub_dmft_points-1
   ien = kien
   call dmft_construct_self_energy(ien, kien, is, self_energy, hub_atom, &
        hub_atomb, connection_table, nmerge, merge_table, merge_neigh, &
        dimer_table, h_atoms, hub, hub_proj_basis, rot_vec_angles, ffirst, &
        llast, cluster, show_matrices, reopen_sigma_file_each_step, &
        force_reopen_file, scratchdir)

   same_self_for_all = .true.
   call dmft_upfold_self_energy(self_energy, hub_overlap, nkpoints, &
        kpoint_data, same_self_for_all, 1)
   same_self_for_all = .false.

   kien = pub_dmft_points-1
   ien = kien
   call dmft_embedding_read_potential(ham, ham_embed, is, ien, fermi_e, &
        ffirst, llast, reopen_sigma_file_each_step, force_reopen_file, &
        check_embed, check_embed_h, ngwf_basis)

   if(nkpoints == 1)then
      call dense_convert(matrixSigStat(is), self_energy%mat )
   else
      call dense_convert(matrixSigStat(is), kpoint_data%self_energy )
   endif
   if(check_embed_h) then
      call dense_create(tmp1, matrixSigStat(is))
      call dense_convert(tmp1, ham_embed(is))
      call dense_axpy(matrixSigStat(is), tmp1, 1.0_DP)
      call dense_destroy(tmp1)
   endif

   force_reopen_file = .false.

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_build_self_infinity'

 end subroutine

 subroutine dmft_prepare_for_dmft_density_kernel(ham, greenf, self_energy, &
      self_infinity, overlap, h_atoms, hub, hub_overlap, hub_proj_basis, &
      ngwf_basis, kpoint_data, sc_data, fermi_e, is, ien, ffirst, llast, &
      hub_atom, hub_atomb, nkpoints, kkk1, kkk2, en_start, en_step, &
      nspecial_frequ, nmu_step, matrixSigStat, matrixAsig, overlap_cplx, &
      rot_vec_angles, connection_table, dimer_table, cluster, nmerge, &
      merge_table, merge_neigh, ham_embed, check_embed, check_embed_h, &
      force_reopen_file, full_prec_dmft_energy, optimisation_parallel, &
      reopen_sigma_file_each_step, same_self_for_all, show_matrices, &
      use_gpu_eigenvectors, scratchdir, dmft_self, dmft_z)

   !---------------------------------------------------------------------------!
   ! Prepare the arrays for the calculation of the density kernel via the      !
   ! trace of the Green's function                                             !
   !---------------------------------------------------------------------------!

   use comms, only: comms_reduce, pub_my_proc_id, pub_on_root, &
        pub_total_num_procs
   use dense, only: dense_axpy, dense_convert, dense_create, dense_destroy, &
        dense_scale, dense_get_element, dense_put_element, dense_product, &
        dense_trace
   use file_handling, only: file_handling_remove_file, &
        file_handling_remove_numbered_files
   ! #ifdef GPU_SPEEDUP
   ! use fortran_cuda,      only: eigenvector_gen_cuda_r
   ! #endif
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use linalg,            only: linalg_dsygv_lt, linalg_zhegv_lt
   use parallel_strategy, only: par => pub_par
   use rundat, only: pub_debug, pub_debug_on_root, &
        pub_dmft_points, pub_dmft_temp, pub_num_spins
   use sparse, only: sparse_axpy, sparse_convert, sparse_copy, sparse_create, &
        sparse_destroy, sparse_product, sparse_scale
   use utils,             only: utils_abort, utils_busy_wait, &
        utils_alloc_check, utils_dealloc_check

   implicit none

   ! Arguments
   type(SPAM3),                  intent(inout) :: ham(pub_num_spins)
   type(DMFT_MATRIX),            intent(in   ) :: greenf
   type(DMFT_MATRIX),            intent(inout) :: self_energy
   type(SPAM3),                  intent(inout) :: self_infinity(pub_num_spins)
   type(SPAM3),                  intent(in   ) :: overlap
   type(ARRAY_OF_MATRICES),      intent(in   ) :: h_atoms(par%nat_hub)
   type(HUBBARD_MODEL),          intent(in   ) :: hub
   type(DMFT_OVERLAP),           intent(in   ) :: hub_overlap
   type(FUNC_BASIS),             intent(in   ) :: hub_proj_basis
   type(FUNC_BASIS),             intent(in   ) :: ngwf_basis
   type(DMFT_KPOINT_MATRICES),   intent(inout) :: kpoint_data
   type(DMFT_SC_MATRICES),       intent(inout) :: sc_data
   real(kind = DP),              intent(in   ) :: fermi_e(2)
   integer,                      intent(in   ) :: is
   integer,                      intent(inout) :: ien
   integer,                      intent(in   ) :: ffirst
   integer,                      intent(in   ) :: llast
   integer,                      intent(inout) :: hub_atom
   integer,                      intent(inout) :: hub_atomb
   integer,                      intent(in   ) :: nkpoints
   integer,                      intent(  out) :: kkk1
   integer,                      intent(  out) :: kkk2
   real(kind = DP),              intent(in   ) :: en_start
   real(kind = DP),              intent(in   ) :: en_step
   integer,                      intent(in   ) :: nspecial_frequ
   integer,                      intent(in   ) :: nmu_step
   type(DEM), allocatable,       intent(inout) :: matrixSigStat(:)
   type(DEM), allocatable,       intent(inout) :: matrixAsig(:)
   type(SPAM3),                  intent(inout) :: overlap_cplx
   real(kind = DP), allocatable, intent(in   ) :: rot_vec_angles(:,:,:)
   logical, allocatable,         intent(in   ) :: connection_table(:,:)
   integer, allocatable,         intent(in   ) :: dimer_table(:,:)
   logical,                      intent(in   ) :: cluster
   integer,                      intent(in   ) :: nmerge
   integer, allocatable,         intent(inout) :: merge_table(:,:)
   integer, allocatable,         intent(in   ) :: merge_neigh(:)
   type(SPAM3),                  intent(inout) :: ham_embed(pub_num_spins)
   logical,                      intent(inout) :: check_embed
   logical,                      intent(inout) :: check_embed_h
   logical,                      intent(inout) :: force_reopen_file
   logical,                      intent(in   ) :: full_prec_dmft_energy
   logical,                      intent(in   ) :: optimisation_parallel
   logical,                      intent(in   ) :: reopen_sigma_file_each_step
   logical,                      intent(inout) :: same_self_for_all
   logical,                      intent(in   ) :: show_matrices
   logical,                      intent(in   ) :: use_gpu_eigenvectors
   character(2000),              intent(in   ) :: scratchdir
   type(SPAM3), optional,        intent(inout) :: dmft_self(pub_num_spins)
   type(SPAM3), optional,        intent(inout) :: dmft_z(pub_num_spins)

   ! Local variables
   integer                         :: ii, i, j, m, iistart, iistep, kien
   real(kind = DP)                 :: chem, beta, occup, mintail, dd, ddd, &
                                      tmpdiff
   type(SPAM3)                     :: scratch
   complex(kind = DP)              :: dddc
   logical                         :: skip_it, skip_eigen_calc
   type(DEM)                       :: tmp1
   type(DEM)                       :: tmp2
   real(kind = DP), allocatable    :: tmpr(:,:)
   complex(kind = DP), allocatable :: tmpc(:,:)
   real(kind = DP), allocatable    :: overlap_square_r(:,:)
   complex(kind = DP), allocatable :: overlap_square_c(:,:)
   integer                         :: ierr
   complex(kind=DP)                :: dense_el
   complex(kind=DP)                :: dense_new_el
   type(DEM)                       :: prefactor
   type(DEM)                       :: dens_tmp
   character                       :: op

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_prepare_for_dmft_density_kernel'

   allocate(tmpr(ngwf_basis%num, ngwf_basis%num), stat = ierr)
   call utils_alloc_check('dmft_prepare_for_dmft_density_kernel', 'tmpr', ierr)
   allocate(tmpc(ngwf_basis%num, ngwf_basis%num), stat = ierr)
   call utils_alloc_check('dmft_prepare_for_dmft_density_kernel', 'tmpc', ierr)

   if(pub_on_root .and. pub_dmft_mpi_size > 1 .and. pub_dmft_mpi_rank == 1) then
#ifdef DMFT_SAFE_NFS
      call utils_busy_wait(4000)
#endif
      if (pub_debug) then
         call file_handling_remove_numbered_files( &
              "onetep_dmft_confirmation_dmft_density_kernel_")
         call file_handling_remove_file( &
              "onetep_dmft_confirmation_dmft_density_count")
      end if
   endif

   beta = 1.0_DP/pub_dmft_temp
   chem = fermi_e(is) + pub_dmft_chem_shift
   ii = ngwf_basis%num

   skip_it = .false.
   if(nkpoints == 1)then
      if(nmu_step > 1) skip_it = .true.
   endif

   if (.not. skip_it) then
      force_reopen_file = .true.

      kien = pub_dmft_points-1
      ien = kien

      call dmft_construct_self_energy(ien, kien, is, self_energy, hub_atom, &
           hub_atomb, connection_table, nmerge, merge_table, merge_neigh, &
           dimer_table, h_atoms, hub, hub_proj_basis, rot_vec_angles, ffirst, &
           llast, cluster, show_matrices, reopen_sigma_file_each_step, &
           force_reopen_file, scratchdir)
      same_self_for_all = .true.

      call dmft_upfold_self_energy(self_energy, hub_overlap, nkpoints, &
           kpoint_data, same_self_for_all, 1)
      same_self_for_all = .false.

      !EMBEDDING POTENTIALS
      kien = pub_dmft_points-1
      ien = kien
      call dmft_embedding_read_potential(ham, ham_embed, is, ien, fermi_e, &
           ffirst, llast, reopen_sigma_file_each_step, force_reopen_file, &
           check_embed, check_embed_h, ngwf_basis)
      !END EMBEDDING

      if(present(dmft_self)) call sparse_copy(dmft_self(is), self_energy%mat)

      if(nkpoints == 1)then
         call dense_convert(matrixSigStat(is), self_energy%mat)
         call sparse_copy(self_infinity(is), self_energy%mat, cmplx_to_real = &
              .true.)
      else
         call dense_convert(matrixSigStat(is), kpoint_data%self_energy )
         call sparse_copy(self_infinity(is), kpoint_data%self_energy)
         call dmft_sparse_realpart(self_infinity(is))
      endif
      if(nkpoints == 1)then
         call sparse_axpy(self_infinity(is), ham(is), 1.0_DP)
      else
         call sparse_axpy(self_infinity(is), kpoint_data%ham, 1.0_DP)
      endif

      !EMBEDDING POTENTIALS
      if(check_embed_h) then
         call dense_create(tmp1, matrixSigStat(is))
         call sparse_axpy(self_infinity(is), ham_embed(is), 1.0_DP)
         call dense_convert(tmp1, ham_embed(is))
         call dense_axpy(matrixSigStat(is), tmp1, 1.0_DP)
         call dense_destroy(tmp1)
      endif
      !END EMBEDDING
      ! if(pub_debug_on_root) write(stdout, *) ' MAX SELF ENERGY STATIC (max, &
      !      &min): ', maxval(real(matrixSigStat(is)%zmtx)), &
      !      minval(real(matrixSigStat(is)%zmtx))

      if(present(dmft_z))then
         kien = 1
         ien = kien
         call sparse_create(scratch, overlap, iscmplx = .false.)
         call dmft_construct_self_energy(ien, kien, is, self_energy, hub_atom, &
              hub_atomb, connection_table, nmerge, merge_table, merge_neigh, &
              dimer_table, h_atoms, hub, hub_proj_basis, rot_vec_angles, &
              ffirst, llast, cluster, show_matrices, &
              reopen_sigma_file_each_step, force_reopen_file, scratchdir)
         call dmft_upfold_self_energy(self_energy, hub_overlap, nkpoints, &
              kpoint_data, same_self_for_all)
         call sparse_scale(scratch, 0.0_DP)
         ! -Z
         call sparse_axpy(scratch, self_energy%mat, (0.0_DP, -1.0_DP))
         call sparse_scale(scratch, -1.0_DP/(pi/beta))
         ! S-Z
         call sparse_axpy(scratch, overlap, 1.0_DP)
         call internal_invert_full_real(dmft_z(is), scratch, ngwf_basis%num)
         ! S (S-Z)^-1
         call sparse_product(scratch, overlap, dmft_z(is))
         call sparse_copy(dmft_z(is), scratch)
         call sparse_destroy(scratch)
      endif

      kien = pub_dmft_points-nspecial_frequ-1
      ien = kien
      call dmft_construct_self_energy(ien, kien, is, self_energy, hub_atom, &
           hub_atomb, connection_table, nmerge, merge_table, merge_neigh, &
           dimer_table, h_atoms, hub, hub_proj_basis, rot_vec_angles, ffirst, &
           llast, cluster, show_matrices, reopen_sigma_file_each_step, &
           force_reopen_file, scratchdir)
      same_self_for_all = .true.
      call dmft_upfold_self_energy(self_energy, hub_overlap, nkpoints, &
           kpoint_data, same_self_for_all, 2)
      same_self_for_all = .false.

      if(nkpoints == 1)then
         call dense_convert(matrixAsig(is), self_energy%mat)
      else
         call dense_convert(matrixAsig(is), kpoint_data%self_energy)
      endif

      !EMBEDDING POTENTIALS
      if(check_embed_h)then
         call dmft_embedding_read_potential(ham, ham_embed, is, ien, fermi_e, &
              ffirst, llast, reopen_sigma_file_each_step, force_reopen_file, &
              check_embed, check_embed_h, ngwf_basis)
         call dense_create(tmp2, matrixAsig(is))
         call dense_convert(tmp2, ham_embed(is))
         call dense_axpy(matrixAsig(is), tmp2, 1.0_DP)
         call dense_destroy(tmp2)
      endif
      !END EMBEDDING

      force_reopen_file = .false.

      ! ebl: building self energy tail
      kkk1 = 1
      kkk2 = 1
      mintail = 0.0_DP
      do i = 1, ii
         do j = 1, ii
            call dense_get_element(dense_el, matrixAsig(is), i, j)
            if(abs(aimag(dense_el)) > mintail .and. i == j) &
                 then
               kkk1 = i
               kkk2 = j
               mintail = abs(aimag(dense_el))
            endif
            dense_new_el = -aimag(dense_el)*(en_start + real(ien-1,kind=DP) * &
                 en_step)
            call dense_put_element(dense_new_el, matrixAsig(is), i, j)
         enddo
      enddo

      ! if(pub_debug_on_root) write(stdout, *) 'MAX SELF ENERGY : ', mintail, &
      !      kkk1, kkk2

      if(nkpoints == 1)then
         call sparse_copy(sc_data%ham(is), ham(is))
         call sparse_copy(sc_data%overlap, overlap)
      else
         if(pub_on_root .and. pub_dmft_mpi_rank == 1) write(stdout, *) &
              'loading k-point Hamiltonian'
         call sparse_copy(sc_data%ham(is), kpoint_data%ham)
         if(pub_on_root .and. pub_dmft_mpi_rank == 1) write(stdout, *) &
              'loading k-point overlap'
         call sparse_copy(sc_data%overlap, kpoint_data%overlap)
      endif

      !H v = lambda S v
      !S = L L^T
      !L^-1 H L^-T  x (L^T v) = lambda (L^T v)
      !   C = L^-1 H L^-T  y = L^T v
      ! C y = lambda y
      ! v   = L^-T y
      ! v1^T . v2 = y1^T L^-1 L^-1^T y2

      if(pub_debug_on_root) write(stdout, '(a)') 'DEBUG: start to diagonalize &
           &Hamiltonian'
      if( .not. full_prec_dmft_energy .and. .not. (pub_dmft_use_gpu_onlyme &
           .and. use_gpu_eigenvectors) .and. .not. (use_gpu_eigenvectors .and. &
           pub_dmft_use_gpu_partially) )then
         if(nkpoints == 1)then
            call dmft_eigenvec_fullmat(ii, self_infinity(is), overlap, &
                 sc_data%eigenvecs(is), sc_data%eigenvals(is, :))
         else
            call dmft_eigenvec_fullmat(ii, self_infinity(is), overlap_cplx, &
                 sc_data%eigenvecs(is), sc_data%eigenvals(is, :))
         endif
      else
         skip_eigen_calc = .false.
         if( pub_dmft_mpi_rank /= 1 .and. pub_dmft_use_gpu_partially .and. &
              .not.pub_dmft_splitk) then
            call dense_scale(sc_data%eigenvecs(is), 0.0_DP)
            sc_data%eigenvals(is, :) = 0.0_DP
            skip_eigen_calc = .true.
         endif

         if (.not. skip_eigen_calc) then
            ! call dense_scale(sc_data%eigenvecs(is), 0.0_DP)
            call dense_convert(sc_data%eigenvecs(is), sc_data%ham(is))
            call dense_axpy(sc_data%eigenvecs(is), matrixSigStat(is), 1.0_DP)

            if(pub_on_root) write(stdout, *) '-USE DSYGV-SIZE = ', ii
            if(.not.(pub_dmft_use_gpu_onlyme .and. use_gpu_eigenvectors))then

               if (nkpoints == 1) then
                  allocate(overlap_square_r(ii, ii), stat = ierr)
                  call utils_alloc_check('dmft_prepare_for_dmft_density_&
                       &kernel', 'overlap_square_r', ierr)
                  call sparse_convert(overlap_square_r, sc_data%overlap)
               else
                  allocate(overlap_square_c(ii, ii), stat = ierr)
                  call utils_alloc_check('dmft_prepare_for_dmft_density_&
                       &kernel', 'overlap_square_c', ierr)
                  call sparse_convert(overlap_square_c, sc_data%overlap)
               end if

               if(nkpoints == 1)then
                  call linalg_dsygv_lt(eigenvectors = &
                       sc_data%eigenvecs(is)%dmtx, eigenvalues = &
                       sc_data%eigenvals(is, :), overlap_square = &
                       overlap_square_r, num = ii)
               else
                  call linalg_zhegv_lt(eigenvectors = &
                       sc_data%eigenvecs(is)%zmtx, eigenvalues = &
                       sc_data%eigenvals(is, :), metric = overlap_square_c, &
                       num = ii)
               endif
            else
               ! #ifdef GPU_SPEEDUP
               ! if(nkpoints == 1)then
               !    call eigenvector_gen_cuda_r( ii, sc_data%eigenvecs(is)%dmtx, &
               !         overlap_square_r, sc_data%eigenvals(is, :), tmpr, &
               !         .false.)
               ! else
               !    call eigenvector_cuda_type_c(1, ii, &
               !         sc_data%eigenvecs(is)%zmtx, overlap_square_c, &
               !         sc_data%eigenvals(is, :), tmpc, .false.)
               !    if (pub_debug) then
               !       write(stdout, *) 'CHECKING EIGENVECTOR - VALUES CUDA'
               !       write(stdout, *) '-- > running on CPU'
               !       call sparse_copy( sc_data%overlap, kpoint_data%overlap )
               !       call sparse_convert(overlap_square_c, sc_data%overlap)
               !       sc_data%eigenvals_backup(is, :) = sc_data%eigenvals(is, :)
               !       call linalg_zhegv_lt(eigenvectors = &
               !            sc_data%eigenvecs(is)%zmtx, eigenvalues = &
               !            sc_data%eigenvals(is, :), overlap_square = &
               !            overlap_square_c, num = ii)
               !       write(stdout, *) 'WARNING : CPU MIN-MAX AB VALUE &
               !            &EIGENVECTORS : ', &
               !            minval(abs(sc_data%eigenvecs(is)%zmtx)), &
               !            maxval(abs(sc_data%eigenvecs(is)%zmtx))
               !       write(stdout, *) 'WARNING : GPU MIN-MAX AB VALUE &
               !            &EIGENVECTORS : ', minval(abs(tmpc)), &
               !            maxval(abs(tmpc))
               !       write(stdout, *) 'WARNING : CPU MIN-MAX RE VALUE &
               !            &EIGENVECTORS : ', &
               !            minval(real(sc_data%eigenvecs(is)%zmtx)), &
               !            maxval(real(sc_data%eigenvecs(is, :, :)))
               !       write(stdout, *) 'WARNING : GPU MIN-MAX RE VALUE &
               !            &EIGENVECTORS : ', minval(real(tmpc)), &
               !            maxval(real(tmpc))
               !       write(stdout, *) 'WARNING : CPU MIN-MAX IM VALUE &
               !            &EIGENVECTORS : ', &
               !            minval(aimag(sc_data%eigenvecs(is)%zmtx)), &
               !            maxval(aimag(sc_data%eigenvecs(is, :, :)))
               !       write(stdout, *) 'WARNING : GPU MIN-MAX IM VALUE &
               !            &EIGENVECTORS : ', minval(aimag(tmpc)), &
               !            maxval(aimag(tmpc))
               !       write(stdout, *) 'WARNING : continue with CPU &
               !            &eigenvectors ... safer...'
               !       tmpc = sc_data%eigenvecs(is)%zmtx
               !       tmpdiff = maxval(abs( sc_data%eigenvals_backup(is, &
               !            :)-sc_data%eigenvals(is, :) ) )
               !       write(stdout, *) 'DIFFERENCE IN EIGENVALUES : ', tmpdiff
               !       if(tmpdiff > 1.0E-5_DP)then
               !          write(stdout, *) 'ERROR GPU EIGENVALUES WRONG, DIFF &
               !               &FROM CPU = ', tmpdiff
               !       endif
               !    endif
               ! endif
               ! #else
               call utils_abort('Error in &
                    &dmft_prepare_for_dmft_density_kernel: error in gpu 2')
               ! #endif
               if(nkpoints == 1)then
                  deallocate(overlap_square_r, stat = ierr)
                  call utils_dealloc_check('dmft_prepare_for_dmft_density_&
                       &kernel', 'overlap_square_r', ierr)
                  sc_data%eigenvecs(is)%dmtx = tmpr
               else
                  deallocate(overlap_square_c, stat = ierr)
                  call utils_dealloc_check('dmft_prepare_for_dmft_density_&
                       &kernel', 'overlap_square_c', ierr)
                  sc_data%eigenvecs(is)%zmtx = tmpc
               endif
            endif
         end if

         if(pub_dmft_use_gpu_partially .and. .not.pub_dmft_splitk)then
            if(pub_dmft_mpi_size > 1) then
               call utils_abort("Reached point in code where nfs sync used to &
                    &appear; the code needs to be updated to use these &
                    &settings.")
               ! write(stdout, *) 'I AM DONE RANK : ', pub_dmft_mpi_rank
               ! call sync_via_nfs(sc_data%eigenvecs(is), 'sumup_AA', ii)
               ! call sync_via_nfs(sc_data%eigenvals(is, :), 'sumup_vectorAA', ii)
            endif
         endif
      endif

      ! if (.true.) then
      !    if(pub_debug_on_root) then
      !       write(stdout, *) 'TEST EIGENVALUES : ', sc_data%eigenvals(is, 1:10)
      !       if(nkpoints == 1)then
      !          write(stdout, *) 'TEST EIGENVECTORS: ', &
      !               sc_data%eigenvecs(is)%dmtx(1, 1:10)
      !          write(stdout, *) 'TRANSPOSE : ', &
      !               sc_data%eigenvecs(is)%dmtx(1:10, 1)
      !          ! write(stdout, *) 'MATRIX H : ', sc_data%ham(is)%dmtx(1, 1:10)
      !          write(stdout, *) 'MIN-MAX EIGEN : ', &
      !               minval(sc_data%eigenvals(is, :)), &
      !               maxval(sc_data%eigenvals(is, :))
      !          write(stdout, *) 'MIN-MAX MATRIX : ', &
      !               minval(sc_data%eigenvecs(is)%dmtx), &
      !               maxval(sc_data%eigenvecs(is)%dmtx)
      !          write(stdout, *) 'MIN-MAX COLUMN 1 : ', &
      !               minval(sc_data%eigenvecs(is)%dmtx(:, 1)), &
      !               maxval(sc_data%eigenvecs(is)%dmtx(:, 1))
      !          write(stdout, *) 'MIN-MAX LINE 1 : ', &
      !               minval(sc_data%eigenvecs(is)%dmtx(1, :)), &
      !               maxval(sc_data%eigenvecs(is)%dmtx(1, :))
      !       else
      !          write(stdout, *) 'TEST EIGENVECTORS: ', &
      !               sc_data%eigenvecs(is)%zmtx(1, 1:10)
      !          write(stdout, *) 'TRANSPOSE : ', &
      !               sc_data%eigenvecs(is)%zmtx(1:10, 1)
      !          ! write(stdout, *) 'MATRIX H : ', sc_data%ham(is)%zmtx(1, 1:10)
      !          write(stdout, *) 'MIN-MAX EIGEN : ', &
      !               minval(sc_data%eigenvals(is, :)), &
      !               maxval(sc_data%eigenvals(is, :))
      !          write(stdout, *) 'MIN-MAX MATRIX : ', &
      !               minval(abs(sc_data%eigenvecs(is)%zmtx(:, :))), &
      !               maxval(abs(sc_data%eigenvecs(is)%zmtx))
      !          write(stdout, *) 'MIN-MAX COLUMN 1 : ', &
      !               minval(abs(sc_data%eigenvecs(is)%zmtx(:, 1))), &
      !               maxval(abs(sc_data%eigenvecs(is)%zmtx(:, 1)))
      !          write(stdout, *) 'MIN-MAX LINE 1 : ', &
      !               minval(abs(sc_data%eigenvecs(is)%zmtx(1, :))), &
      !               maxval(abs(sc_data%eigenvecs(is)%zmtx(1, :)))
      !       endif
      !    endif
      ! end if

      if(nkpoints == 1)then
         call sparse_copy( sc_data%overlap, overlap )
      else
         if (pub_debug) then
            write(stdout, *) 'WARNING : copying k-point overlap to &
                 &overlap_cplx - it should not affect the result, just to make &
                 &sure'
            write(stdout, *) ' that overlap_cplx was not affected up to this &
                 &point'
            call sparse_copy(overlap_cplx, kpoint_data%overlap)
         end if
         call sparse_copy( sc_data%overlap, overlap_cplx )
      endif

      !   A = H_LDA + Sig(oo)
      !  (H + Sig) x A_m = lambda_m S x A_m
      !          G_ab(iw)  =       A_m(a) A_m(b) /     (  iw + mu  - lambda_m  )
      ! Sum_iw (G_ab(iw)) = Sum_m A_m(a) A_m(b) / ( 1 + exp (
      ! beta*(lambda_m-mu) ) )

   end if

   if(pub_debug_on_root) write(stdout, '(a)') "DEBUG: building tail of Green's &
        &function"

   ! GF tail traced analytically
   ! Sum_iw (G_ab(iw)) = Sum_m A_m(a) A_m(b) / (1 + exp(beta*(lambda_m-mu)))

   ! Setting up
   call dense_scale(sc_data%tail, 0.0_DP)
   call dense_create(prefactor, ii, ii, iscmplx=.false.)
   call dense_create(dens_tmp, ii, ii, iscmplx=sc_data%eigenvecs(is)%iscmplx)

   ! Casting the exponential prefactor as a diagonal matrix
   do m = 1, ii
      call dense_put_element(1.0_DP/(1.0_DP + internal_dexpc(beta * &
           (sc_data%eigenvals(is, m)-chem))), prefactor, m, m)
   end do

   ! Computing the three-matrix product in two steps
   if (sc_data%eigenvecs(is)%iscmplx) then
      op = "C"
   else
      op = "T"
   end if
   call dense_product(dens_tmp, sc_data%eigenvecs(is), prefactor)
   call dense_product(sc_data%tail, dens_tmp, sc_data%eigenvecs(is), opB=op)

   if(pub_dmft_mpi_size > 1 .and. .not.pub_dmft_splitk) then
      call utils_abort("Reached point in code where nfs sync used to appear; &
           &the code needs to be updated to use these settings.")
      ! write(stdout, *) 'I AM DONE RANK : ', pub_dmft_mpi_rank
      ! call sync_via_nfs(sc_data%tail, 'sumup_tail', ii)
   endif

   if(pub_debug_on_root) write(stdout, '(a)') 'DEBUG: setting up Fermi occupations&
        & of the tail'

   occup = dense_trace(prefactor)

   ! Tidying up
   call dense_destroy(prefactor)
   call dense_destroy(dens_tmp)

   deallocate(tmpr, stat = ierr)
   call utils_dealloc_check('dmft_prepare_for_dmft_density_kernel', 'tmpr', ierr)
   deallocate(tmpc, stat = ierr)
   call utils_dealloc_check('dmft_prepare_for_dmft_density_kernel', 'tmpc', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_prepare_for_dmft_density_kernel'

 end subroutine

 subroutine dmft_compute_dmft_energy(dmft_Energy, greenf, self_energy, &
      ngwf_basis, kpoint_data, sc_data, energy, fermi_e, nmu, nmu_step, is, &
      ien, nkpoints, kkk1, kkk2, nspecial_frequ, ham_embed, check_embed_h, &
      matrixDensKernel, matrixSigStat, matrixAsig)

   !---------------------------------------------------------------------------!
   ! Computes the DMFT energy via the Migdal formula                           !
   !---------------------------------------------------------------------------!

   ! #ifdef GPU_SPEEDUP
   ! use fortran_cuda,   only: matmulcuda_c, matmulcuda_c_singleprec_bypass
   ! #endif
   use sparse, only: sparse_axpy, sparse_convert, sparse_copy, sparse_create, &
        sparse_destroy, sparse_axpy, sparse_scale
   use dense, only: dense_axpy, dense_convert, dense_create, dense_destroy, &
        dense_scale, dense_product, dense_put_element, dense_copy
   use utils,          only: utils_abort, utils_alloc_check, utils_dealloc_check
   use comms,          only: pub_on_root
   use rundat, only: pub_debug, pub_debug_on_root, pub_num_spins, pub_dmft_temp
   use function_basis, only: FUNC_BASIS

   implicit none

   ! Arguments
   real(kind = DP),            intent(inout) :: dmft_Energy(2)
   type(DMFT_MATRIX),          intent(in   ) :: greenf
   type(DMFT_MATRIX),          intent(in   ) :: self_energy
   type(FUNC_BASIS),           intent(in   ) :: ngwf_basis
   type(DMFT_KPOINT_MATRICES), intent(in   ) :: kpoint_data
   type(DMFT_SC_MATRICES),     intent(in   ) :: sc_data
   complex(kind = DP),         intent(in   ) :: energy
   real(kind = DP),            intent(in   ) :: fermi_e(2)
   integer,                    intent(in   ) :: nmu
   integer,                    intent(in   ) :: nmu_step
   integer,                    intent(in   ) :: is
   integer,                    intent(in   ) :: ien
   integer,                    intent(in   ) :: nkpoints
   integer,                    intent(in   ) :: kkk1
   integer,                    intent(in   ) :: kkk2
   integer,                    intent(in   ) :: nspecial_frequ
   type(SPAM3),                intent(in   ) :: ham_embed(pub_num_spins)
   logical,                    intent(in   ) :: check_embed_h
   type(SPAM3),                intent(inout) :: matrixDensKernel
   type(DEM), allocatable,     intent(in   ) :: matrixSigStat(:)
   type(DEM), allocatable,     intent(in   ) :: matrixAsig(:)

   ! Local variables
   complex(kind = DP)              :: HU
   integer                         :: i, j, k, l
   integer                         :: ii, m, INFO
   real(kind = DP)                 :: chem, beta
   complex(kind = DP)              :: cc
   complex(kind = DP), parameter   :: imi = (0.0_DP, 1.0_DP)
   logical                         :: check
   complex(kind = DP)              :: tempd, tempdd
   complex(kind = DP), allocatable :: matrixGan(:,:)
   complex(kind = DP), allocatable :: tmp1(:,:)
   complex(kind = DP), allocatable :: tmp2(:,:)
   complex(kind = DP), allocatable :: tmp3(:,:)
   complex(kind = DP), allocatable :: tmpc(:,:)
   complex(kind = DP), allocatable :: tmpc2(:,:)
   complex(kind = DP), allocatable :: tmpc3(:,:)
   complex(kind = DP), allocatable :: ham_tmp(:,:)
   complex(kind = DP), allocatable :: overlap_tmp(:,:)
   type(SPAM3)                     :: tmp_spam1
   type(SPAM3)                     :: tmp_spam2
   type(DEM)                       :: prefactor
   type(DEM)                       :: tmp_dem1
   type(DEM)                       :: tmp_dem2
   type(DEM)                       :: eigenvecs_cmplx
   integer                         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_compute_dmft_energy'

   ! TODO: convert all tmp matrices to DEM/SPAM3 type. This will make lots of
   ! the code far simpler

   beta = 1.0_DP/pub_dmft_temp
   chem = fermi_e(is) + pub_dmft_chem_shift
   ii = ngwf_basis%num

   ! Allocation of matrices
   allocate(matrixGan(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_compute_dmft_energy', 'matrixGan', ierr)
   allocate(tmpc3(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_compute_dmft_energy', 'tmpc3', ierr)
   allocate(tmp1(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_compute_dmft_energy', 'tmp1', ierr)
   allocate(tmp2(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_compute_dmft_energy', 'tmp2', ierr)
   allocate(tmp3(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_compute_dmft_energy', 'tmp3', ierr)
   allocate(tmpc(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_compute_dmft_energy', 'tmpc', ierr)
   allocate(tmpc2(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_compute_dmft_energy', 'tmpc2', ierr)
   allocate(ham_tmp(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_compute_dmft_energy', 'ham_tmp', ierr)
   allocate(overlap_tmp(ii, ii), stat = ierr)
   call utils_alloc_check('dmft_compute_dmft_energy', 'overlap_tmp', ierr)

   if(.not.pub_dmft_split)Then
      if(nkpoints == 1)then
         ! call sparse_copy(tmp1, self_energy%mat )
         call sparse_convert(tmp1, self_energy%mat )
      else
         ! call sparse_copy(tmp1, kpoint_data%self_energy)
         call sparse_convert(tmp1, kpoint_data%self_energy)
      endif
   else
      tmp1 = self_energy%backup
      ! call sparse_convert(tmp1, self_energy%backup)
   endif

   if(check_embed_h) then
      call sparse_convert(tmp2, ham_embed(is))
      tmp1 = tmp1 + tmp2
      ! call sparse_axpy(tmp1, ham_embed(is), 1.0_DP)
   endif

   if(.not.greenf%backup_flag .and. .not.pub_dmft_split) then
      ! call sparse_copy(tmp2, greenf%mat)
      call sparse_convert(tmp2, greenf%mat)
   else
      tmp2 = greenf%backup
      ! call sparse_convert(tmp2, greenf%backup)
   endif

   if (pub_debug) then
      if(nkpoints > 1)then
         write(stdout, *) 'COMPUTING G(-iw)'
         ! call dense_convert(tmp3, tmp1)
         ! call dense_scale(tmp3, -1.0_DP)
         ! call dense_axpy(tmp3, sc_data%overlap, energy)
         ! call dense_axpy(tmp3, sc_data%ham(is)%zmtx, -1.0_DP)
         call sparse_convert(ham_tmp, sc_data%ham(i))
         call sparse_convert(overlap_tmp, sc_data%overlap)
         tmpc = energy*overlap_tmp - ham_tmp - tmp1
         call dmft_sym_invert(size(tmpc, 1), tmpc)
         write(stdout, *) ' iw + mu : ', energy
         write(stdout, *) '-iw + mu : ', conjg(energy)
         write(stdout, *) 'MAX DIFFERENCE BETWEEN G(iw) AND G(iw) : '
         write(stdout, *) maxval(abs( tmpc - tmp2 ))
         tmpc = energy*overlap_tmp - ham_tmp
         call dmft_sym_invert(size(tmpc, 1), tmpc)
         tmpc2 = tmpc
         tmpc = conjg(energy)*overlap_tmp - ham_tmp
         call dmft_sym_invert(size(tmpc, 1), tmpc)
         write(stdout, *) 'MAX DIFFERENCE BETWEEN G(iw) AND G(-iw)^ + : '
         write(stdout, *) maxval(abs( transpose(conjg(tmpc)) - tmpc2 ))
         tmpc = conjg(energy)*conjg(ham_tmp) - conjg(ham_tmp)
         call dmft_sym_invert(size(tmpc, 1), tmpc)
         write(stdout, *) 'MAX DIFFERENCE BETWEEN G(iw) AND G(-k, -iw)* : '
         write(stdout, *) maxval(abs( conjg(tmpc) - tmpc2 ))
      endif
   end if

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: computing G_analytic &
        &term'

   if(.not.pub_dmft_use_gpu_onlyme)then
      ! Initialising
      call dense_create(prefactor,  ii, ii,     iscmplx=.true.)
      call dense_create(tmp_dem1,   ii, ii,     iscmplx=.true.)
      call dense_create(tmp_dem2,   ii, ii,     iscmplx=.true.)
      call sparse_create(tmp_spam2, greenf%mat, iscmplx=.true.)
      call dense_create(eigenvecs_cmplx, sc_data%eigenvecs(is), iscmplx=.true.)
      call dense_copy(eigenvecs_cmplx, sc_data%eigenvecs(is))
      do m = 1, ii
         call dense_put_element(1.0_DP/( energy - sc_data%eigenvals(is, m)), &
              prefactor, m, m)
      end do

      ! Performing matrix product eigenvec * prefactor * eigenvec^dag
      call dense_product(tmp_dem1, eigenvecs_cmplx, prefactor)
      if (sc_data%eigenvecs(is)%iscmplx) then
         call dense_product(tmp_dem2, tmp_dem1, eigenvecs_cmplx, opB="C")
      else
         call dense_product(tmp_dem2, tmp_dem1, eigenvecs_cmplx, opB="T")
      end if

      ! Converting result to a standard complex array
      call dense_convert(tmp_spam2, tmp_dem2)
      call sparse_convert(matrixGan, tmp_spam2)

      ! Tidying up
      call dense_destroy(eigenvecs_cmplx)
      call sparse_destroy(tmp_spam2)
      call dense_destroy(prefactor)
      call dense_destroy(tmp_dem1)
      call dense_destroy(tmp_dem2)

      ! Old code:
      ! do m = 1, ii
      !    tempd = 1.0_DP/( energy - sc_data%eigenvals(is, m))
      !    if (nkpoints == 1) then
      !       tmpc(:, m) = sc_data%eigenvecs(is)%dmtx(:, m)*tempd
      !    else
      !       tmpc(:, m) = sc_data%eigenvecs(is)%zmtx(:, m)*tempd
      !    end if
      ! enddo
      ! dense_tmp
      ! if (nkpoints == 1) then
      !
      !    matrixGan = MATMUL(tmpc, sc_data%eigenvecs_t(is)%dmtx)
      ! else
      !    matrixGan = MATMUL(tmpc, sc_data%eigenvecs_t(is)%zmtx)
      ! end if
      if (pub_debug) tmpc2 = tmp2
      tmp2 = tmp2-matrixGan
   else
      ! #ifdef GPU_SPEEDUP
      ! write(stdout, *) 'Ganalytic on GPU'
      ! do m = 1, ii
      !    tempd = 1.0_DP/( energy - sc_data%eigenvals(is, m))
      !    if (nkpoints == 1) then
      !       tmpc(:, m) = sc_data%eigenvecs(is)%dmtx(:, m)*tempd
      !    else
      !       tmpc(:, m) = sc_data%eigenvecs(is)%zmtx(:, m)*tempd
      !    end if
      ! enddo
      ! if( internal_doubleprec(energy, chem, fermi_e, is, ien,
      ! nspecial_frequ) ) then
      ! call matmulcuda_c(tmpc, sc_data%eigenvecs_tc(is, :, :), matrixGan, ii,
      ! ii, ii)
      ! else
      ! call matmulcuda_c_singleprec_bypass(tmpc, sc_data%eigenvecs_tc(is, :,
      ! :), matrixGan, ii, ii, ii)
      ! endif
      ! #else
      call utils_abort("Error in dmft_compute_energy: should not have arrived &
           &here")
      ! #endif
      if (pub_debug) tmpc2 = tmp2
      tmp2 = tmp2-matrixGan
   endif

   ! if (pub_debug) then
   !    tmpc3 = tmp2
   !    tmpc = matrixGan
   !    tmp2 = tmpc2
   !    matrixGan = 0.0_DP
   !    do m = 1, ii
   !       tempd = 1.0_DP/( energy - sc_data%eigenvals(is, m))
   !       do i = 1, ii
   !          if (nkpoints == 1) then
   !             tempdd = tempd * sc_data%eigenvecs(is)%dmtx(i, m)
   !          else
   !             tempdd = tempd * sc_data%eigenvecs(is)%zmtx(i, m)
   !          end if
   !          do j = 1, ii
   !             if (nkpoints == 1) then
   !                cc = sc_data%eigenvecs(is)%dmtx(j, m) * tempdd
   !             else
   !                cc = sc_data%eigenvecs(is)%zmtx(j, m) * tempdd
   !             end if
   !             matrixGan(i, j) = matrixGan(i, j) + cc
   !             tmp2(i, j)      = tmp2(i, j)      - cc
   !          enddo
   !       enddo
   !       ! if(nkpoints == 1)then
   !       !    do i = 1, ii
   !       !       tempdd = tempd * sc_data%eigenvecs(is, i, m)
   !       !       do j = 1, ii
   !       !          cc = sc_data%eigenvecs(is, j, m) * tempdd
   !       !          matrixGan(i, j) = matrixGan(i, j) + cc
   !       !          tmp2(i, j)      = tmp2(i, j)      - cc
   !       !       enddo
   !       !    enddo
   !       ! else
   !       !    do i = 1, ii
   !       !       tempdd = tempd * sc_data%eigenvecs_c(is, i, m)
   !       !       do j = 1, ii
   !       !          cc = conjg(sc_data%eigenvecs_c(is, j, m)) * tempdd
   !       !          matrixGan(i, j) = matrixGan(i, j) + cc
   !       !          tmp2(i, j)      = tmp2(i, j)      - cc
   !       !       enddo
   !       !    enddo
   !       ! endif
   !    enddo

   !    if(maxval(abs( tmp2 - tmpc3 )) > 1.0E-5_DP)then
   !       write(stdout, *) 'ERROR TMP2 in GPU OPTIMIZATION'
   !       write(stdout, *) ' DIFF       : ', maxval(abs( tmp2 - tmpc3 ))
   !       write(stdout, *) ' doubleprec : ', internal_doubleprec(energy, chem, &
   !            fermi_e, is, ien, nspecial_frequ)
   !    end if
   !    if(maxval(abs( matrixGan - tmpc )) > 1.0E-5_DP)then
   !       write(stdout, *) 'ERROR GPU optimization - matrixGan'
   !       write(stdout, *) 'DIFF = ', maxval(abs( matrixGan - tmpc ))
   !    end if
   ! end if

   HU = 0.0_DP

   ! cw: only compute energy once the chemical potential is converged
   if(nmu_step == nmu)then
      ! Sigma * G_num
      do i = 1, ii
         do m = 1, ii
            ! for positive and negative matsubara frequencies at the same time
            ! Tr ( Sigma * G_num )
            HU = HU + 2.0_DP*real( tmp1(i, m)*tmp2(m, i) )
         enddo
      enddo

      ! Calculating Sigma_num
      ! Create workspace
      call dense_create(tmp_dem1, ii, ii, .false.)
      call sparse_create(tmp_spam1, greenf%mat, iscmplx=.false.)
      call sparse_create(tmp_spam2, greenf%mat, iscmplx=.true.)

      ! Take real part
      call dense_copy(tmp_dem1, matrixSigStat(is), cmplx_to_real = .true.)

      ! Convert to sparse matrices
      call dense_convert(tmp_spam1, tmp_dem1)
      call dense_convert(tmp_spam2, matrixAsig(is))

      ! Perform algebra
      call sparse_scale(tmp_spam2, 1.0_DP/(imi*internal_myfrequ(energy)))
      call sparse_axpy(tmp_spam2, tmp_spam1, 1.0_DP)

      ! Convert to array and subtract from Sigma to obtain Sigma_num
      call sparse_convert(tmp3, tmp_spam2)
      tmp1 = tmp1-tmp3

      ! Destroying temporary arrays
      call dense_destroy(tmp_dem1)
      call sparse_destroy(tmp_spam1)
      call sparse_destroy(tmp_spam2)

      ! Sigma_num * G_an
      do i = 1, ii
         do m = 1, ii
            ! for positive and negative matsubara frequencies at the same time
            ! Tr ( Sig_num * G_an )
            HU = HU + 2.0_DP*real( tmp1(i, m)*matrixGan(m, i) )
         enddo
      enddo
   endif

   dmft_Energy(is) = dmft_Energy(is) + HU*0.5_DP*pub_dmft_temp

   ! conjg transpose tmp2 corresponds to negative matsubara frequencies !
   ! if(nkpoints == 1)then
   ! matrixDensKernel%dmtx = matrixDensKernel%dmtx + tmp2 +
   ! conjg(transpose(tmp2))
   ! else
   ! matrixDensKernel%zmtx = matrixDensKernel%zmtx + tmp2 +
   ! conjg(transpose(tmp2))
   ! endif
   call sparse_create(tmp_spam1, matrixDensKernel)
   if (nkpoints == 1) then
      call sparse_convert(tmp_spam1, real(tmp2 + conjg(transpose(tmp2))))
   else
      call sparse_convert(tmp_spam1, tmp2 + conjg(transpose(tmp2)))
   end if
   call sparse_axpy(matrixDensKernel, tmp_spam1, 1.0_DP)
   call sparse_destroy(tmp_spam1)

   ! Deallocation of matrices
   deallocate(matrixGan, stat = ierr)
   call utils_dealloc_check('dmft_compute_dmft_energy', 'matrixGan', ierr)
   deallocate(tmpc3, stat = ierr)
   call utils_dealloc_check('dmft_compute_dmft_energy', 'tmpc3', ierr)
   deallocate(tmp1, stat = ierr)
   call utils_dealloc_check('dmft_compute_dmft_energy', 'tmp1', ierr)
   deallocate(tmp2, stat = ierr)
   call utils_dealloc_check('dmft_compute_dmft_energy', 'tmp2', ierr)
   deallocate(tmp3, stat = ierr)
   call utils_dealloc_check('dmft_compute_dmft_energy', 'tmp3', ierr)
   deallocate(tmpc, stat = ierr)
   call utils_dealloc_check('dmft_compute_dmft_energy', 'tmpc', ierr)
   deallocate(tmpc2, stat = ierr)
   call utils_dealloc_check('dmft_compute_dmft_energy', 'tmpc2', ierr)
   deallocate(ham_tmp, stat = ierr)
   call utils_dealloc_check('dmft_compute_dmft_energy', 'ham_tmp', ierr)
   ! call sparse_destroy(tmp1)
   ! call sparse_destroy(tmp2)
   ! call dense_destroy(tmp2)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_compute_dmft_energy'

 end subroutine

 subroutine dmft_dump_data_optics(nabla_square, greenf, ngwf_basis, mdl, &
      energy, fermi_e, chem_shift, ien, is, nkpoints, ikpoint, norm_kpoints, &
      ffirst, splitk_first, nspecial_frequ, one_k_one_file)

   !---------------------------------------------------------------------------!
   ! Prints optical properties to file                                         !
   !---------------------------------------------------------------------------!

   ! #ifdef GPU_SPEEDUP
   ! use fortran_cuda,   only: matmul_my_cuda, matmulcuda_r_cublas
   ! #endif
   use sparse,         only: sparse_convert
   use model_type,     only: MODEL
   use utils,          only: utils_abort, utils_unit, utils_int_to_str, &
                             utils_alloc_check, utils_dealloc_check
   use comms, only: comms_barrier, comms_bcast, pub_my_proc_id, pub_on_root, &
        pub_total_num_procs
   use rundat, only: pub_debug_on_root, pub_dmft_complex_freq, &
        pub_dmft_optics, pub_dmft_optics_i1, pub_dmft_optics_i2, &
        pub_dmft_optics_window, pub_dmft_points, pub_dmft_temp
   use function_basis, only: FUNC_BASIS

   implicit none

   ! Arguments
   real(kind = DP), allocatable, intent(inout) :: nabla_square(:,:,:)
   type(DMFT_MATRIX),            intent(in   ) :: greenf
   type(FUNC_BASIS),             intent(in   ) :: ngwf_basis
   type(MODEL),                  intent(in   ) :: mdl
   complex(kind = DP),           intent(in   ) :: energy
   real(kind = DP),              intent(in   ) :: fermi_e
   real(kind = DP),              intent(in   ) :: chem_shift
   integer,                      intent(in   ) :: ien
   integer,                      intent(in   ) :: is
   integer,                      intent(in   ) :: nkpoints
   integer,                      intent(in   ) :: ikpoint
   real(kind = DP), allocatable, intent(in   ) :: norm_kpoints(:)
   integer,                      intent(in   ) :: splitk_first
   integer,                      intent(in   ) :: ffirst
   integer,                      intent(in   ) :: nspecial_frequ
   logical,                      intent(in   ) :: one_k_one_file

   ! Local variables
   integer                         :: i, j, k, nn, sizbuffer, kkkk
   integer                         :: funit
   real(kind = DP), allocatable    :: temp1_real(:,:), temp2_real(:,:)
   complex(kind = DP), allocatable :: temp2_complex(:,:)
   integer                         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_dump_data_optics'

   if(pub_dmft_temp > 0.0_DP .or. .not.pub_dmft_optics) then
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_dump_data_optics'
      return
   end if

   allocate(temp1_real(ngwf_basis%num, ngwf_basis%num), stat = ierr)
   call utils_alloc_check('dmft_dump_data_optics', 'temp1_real', ierr)
   allocate(temp2_real(ngwf_basis%num, ngwf_basis%num), stat = ierr)
   call utils_alloc_check('dmft_dump_data_optics', 'temp2_real', ierr)
   allocate(temp2_complex(ngwf_basis%num, ngwf_basis%num), stat = ierr)
   call utils_alloc_check('dmft_dump_data_optics', 'temp2_complex', ierr)

   if(abs(real(energy) - (fermi_e + chem_shift) ) < pub_dmft_optics_window &
        .and. (ien < pub_dmft_points-nspecial_frequ + 1)) then

      if(any(abs(norm_kpoints) < 1.0E-3_DP))then
         call utils_abort("Error in dmft_dump_data_optics: optics won't work &
              &if some k-points are excluded")
      endif
      nn = ngwf_basis%num
      write(stdout, *) ' for frequ / energy : ', ien, energy
      write(stdout, *) 'convert green SPAM to dense'
      if(.not.pub_dmft_split)then
         call sparse_convert(temp2_complex, greenf%mat)
      else
         temp2_complex = greenf%backup
      endif
      temp1_real = (temp2_complex-conjg(transpose(temp2_complex))) / &
           (0.0_DP, 2.0_DP) / PI
      !temp1_real = aimag(temp2_complex)/pi !BUG : faster for symmetric GF
      ! (Gamma = 0 only)
      write(stdout, *) ' G*NABMA_x '
      ! #ifdef GPU_SPEEDUP
      ! if(pub_dmft_use_gpu_onlyme)then
      !    if(pub_on_root .or. pub_dmft_split)then
      !       CALL matmulcuda_r(temp1_real, nabla_square(:, :, &
      !            pub_dmft_optics_i1), temp2_real, nn, nn, nn)
      !    endif
      !    if(pub_total_num_procs > 1 .and. .not.pub_dmft_split) call &
      !         comms_bcast(0, temp2_real)
      ! else
      !    temp2_real = MATMUL(temp1_real, nabla_square(:, :, pub_dmft_optics_i1))
      ! endif
      ! #else
      temp2_real = MATMUL(temp1_real, nabla_square(:, :, pub_dmft_optics_i1))
      ! #endif
      if(.not.pub_dmft_split)then
         if(pub_on_root)then
            funit = utils_unit()
            if((ien == ffirst .and. one_k_one_file) .or. (.not.one_k_one_file &
                 .and. ien == ffirst .and. ikpoint == 1))then
               open(unit = funit, file = 'optics_data_spin_' // &
                    trim(adjustl(utils_int_to_str(is))) // '_k_' // &
                    trim(adjustl(utils_int_to_str(ikpoint))) // '_rank_' // &
                    trim(adjustl(utils_int_to_str(pub_dmft_mpi_rank))), form = &
                    'unformatted')
            else
               open(unit = funit, file = 'optics_data_spin_' // &
                    trim(adjustl(utils_int_to_str(is))) // '_k_' // &
                    trim(adjustl(utils_int_to_str(ikpoint))) // '_rank_' // &
                    trim(adjustl(utils_int_to_str(pub_dmft_mpi_rank))), form = &
                    'unformatted', position = 'append')
            endif
         endif
      endif

      if(.not.pub_dmft_split)then
         if(pub_on_root)then
            if(pub_dmft_optics_i1 == pub_dmft_optics_i2)then
               write(funit) &
                    nkpoints, norm_kpoints(ikpoint), 1, &
                    internal_volume_cell_angstrom(mdl%cell), pub_dmft_temp, &
                    pub_dmft_points, pub_dmft_mpi_size, nn, ien, energy, &
                    fermi_e, chem_shift, temp2_real
            else
               write(funit) &
                    nkpoints, norm_kpoints(ikpoint), 2, &
                    internal_volume_cell_angstrom(mdl%cell), pub_dmft_temp, &
                    pub_dmft_points, pub_dmft_mpi_size, nn, ien, energy, &
                    fermi_e, chem_shift, temp2_real
            endif
         endif
      else
         do kkkk = 1, pub_total_num_procs
            if(kkkk == pub_my_proc_id + 1)then
               funit = utils_unit()
               if(ien == ffirst .and. kkkk == 1 .and. ( (.not.pub_dmft_splitk &
                    .and. (ikpoint == 1 .or. one_k_one_file)) .or. &
                    (pub_dmft_splitk .and. ikpoint == splitk_first) ) )then
                  open(funit, file = 'optics_data_spin_' // &
                       trim(adjustl(utils_int_to_str(is))) // '_k_' // &
                       trim(adjustl(utils_int_to_str(ikpoint))) // '_rank_' // &
                       trim(adjustl(utils_int_to_str(pub_dmft_mpi_rank))), &
                       form = 'unformatted', status = 'replace')
               else
                  open(funit, file = 'optics_data_spin_' // &
                       trim(adjustl(utils_int_to_str(is))) // '_k_' // &
                       trim(adjustl(utils_int_to_str(ikpoint))) // '_rank_' // &
                       trim(adjustl(utils_int_to_str(pub_dmft_mpi_rank))), &
                       form = 'unformatted', position = 'append', status = 'old')
               endif
               if(pub_dmft_optics_i1 == pub_dmft_optics_i2)then
                  write(funit) nkpoints, norm_kpoints(ikpoint), 1, &
                       internal_volume_cell_angstrom(mdl%cell), pub_dmft_temp, &
                       pub_dmft_points, pub_dmft_mpi_size, nn, ien, energy, &
                       fermi_e, chem_shift, temp2_real
               else
                  write(funit) nkpoints, norm_kpoints(ikpoint), 2, &
                       internal_volume_cell_angstrom(mdl%cell), pub_dmft_temp, &
                       pub_dmft_points, pub_dmft_mpi_size, nn, ien, energy, &
                       fermi_e, chem_shift, temp2_real
               endif
               close(funit)

            endif
            call comms_barrier
         enddo
      endif

      if(pub_dmft_optics_i1 /= pub_dmft_optics_i2)then
         ! #ifdef GPU_SPEEDUP
         ! if(pub_dmft_use_gpu_onlyme)then
         !    if(pub_on_root .or. pub_dmft_split)then
         !       CALL matmulcuda_r(temp1_real, nabla_square(:, :, &
         !            pub_dmft_optics_i2), temp2_real, nn, nn, nn)
         !    endif
         !    if(pub_total_num_procs > 1 .and. .not.pub_dmft_split) call &
         !         comms_bcast(0, temp2_real)
         ! else
         !    temp2_real = MATMUL(temp1_real, nabla_square(:, :, &
         !         pub_dmft_optics_i2))
         ! endif
         ! #else
         temp2_real = MATMUL(temp1_real, nabla_square(:, :, pub_dmft_optics_i2))
         ! #endif
      endif

      if(pub_dmft_optics_i1 /= pub_dmft_optics_i2)then
         if(.not.pub_dmft_split)then
            if(pub_on_root)then
               write(funit) &
                    nkpoints, norm_kpoints(ikpoint), 2, &
                    internal_volume_cell_angstrom(mdl%cell), pub_dmft_temp, &
                    pub_dmft_points, pub_dmft_mpi_size, nn, ien, energy, &
                    fermi_e, chem_shift, temp2_real
            endif
         else
            do kkkk = 1, pub_total_num_procs
               if(kkkk == pub_my_proc_id + 1)then
                  funit = utils_unit()
                  if(ien == ffirst .and. kkkk == 1 .and. ( &
                       (.not.pub_dmft_splitk .and. (ikpoint == 1 .or. &
                       one_k_one_file)) .or. (pub_dmft_splitk .and. ikpoint == &
                       splitk_first) ) )then
                     open(unit = funit, file = 'optics_data_spin_' // &
                          trim(adjustl(utils_int_to_str(is))) // '_k_' // &
                          trim(adjustl(utils_int_to_str(ikpoint))) // '_rank_' &
                          // trim(adjustl(utils_int_to_str( &
                          pub_dmft_mpi_rank))), form = 'unformatted', status &
                          = 'replace')
                  else
                     open(unit = funit, file = 'optics_data_spin_' // &
                          trim(adjustl(utils_int_to_str(is))) // '_k_' // &
                          trim(adjustl(utils_int_to_str(ikpoint))) // '_rank_' &
                          // trim(adjustl(utils_int_to_str( &
                          pub_dmft_mpi_rank))), form = 'unformatted', position &
                          = 'append', status = 'old')
                  endif
                  write(funit) nkpoints, norm_kpoints(ikpoint), 2, &
                       internal_volume_cell_angstrom(mdl%cell), pub_dmft_temp, &
                       pub_dmft_points, pub_dmft_mpi_size, nn, ien, energy, &
                       fermi_e, chem_shift, temp2_real
                  close(funit)
               endif
               call comms_barrier
            enddo
         endif
      endif

      if(.not.pub_dmft_split)then
         if(pub_on_root)then
            close(funit)
         endif
      endif

   endif

   deallocate(temp1_real, stat = ierr)
   call utils_dealloc_check('dmft_dump_data_optics', 'temp1_real', ierr)
   deallocate(temp2_real, stat = ierr)
   call utils_dealloc_check('dmft_dump_data_optics', 'temp2_real', ierr)
   deallocate(temp2_complex, stat = ierr)
   call utils_dealloc_check('dmft_dump_data_optics', 'temp2_complex', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_dump_data_optics'
   return

 end subroutine

 subroutine dmft_show_dos(greenf, overlap, ngwf_basis, energy, is, ikpoint, &
      nkpoints, norm_kpoints, ien, ffirst, fermi_e, kp_file_shift, &
      splitk_first, nplanes, orb_in_plane, one_k_one_file)

   !---------------------------------------------------------------------------!
   ! Computes and writes the density of states from the Green's function       !
   !---------------------------------------------------------------------------!

   ! #ifdef GPU_SPEEDUP
   ! use fortran_cuda,   only: matmulcuda_c, matmulcuda_c_cublas
   ! #endif
   use constants,      only: HARTREE_IN_EVS
   use sparse,         only: sparse_convert
   use comms, only: comms_barrier, comms_bcast, pub_my_proc_id, pub_on_root, &
        pub_total_num_procs
   use utils,          only: utils_alloc_check, utils_dealloc_check, &
        utils_unit, utils_int_to_str, utils_qc_print
   use rundat, only: pub_debug, pub_debug_on_root, pub_dmft_points, &
        pub_rootname, pub_print_qc
   use function_basis, only: FUNC_BASIS

   implicit none

   ! Arguments
   type(DMFT_MATRIX),            intent(in   ) :: greenf
   type(SPAM3),                  intent(in   ) :: overlap
   type(FUNC_BASIS),             intent(in   ) :: ngwf_basis
   complex(kind = DP),           intent(in   ) :: energy
   integer,                      intent(in   ) :: is
   integer,                      intent(in   ) :: ikpoint
   integer,                      intent(in   ) :: nkpoints
   real(kind = DP), allocatable, intent(in   ) :: norm_kpoints(:)
   integer,                      intent(in   ) :: ien
   integer,                      intent(in   ) :: ffirst
   real(kind = DP),              intent(in   ) :: fermi_e(2)
   integer,                      intent(in   ) :: kp_file_shift
   integer,                      intent(in   ) :: splitk_first
   integer,                      intent(in   ) :: nplanes
   logical, allocatable,         intent(in   ) :: orb_in_plane(:,:)
   logical,                      intent(in   ) :: one_k_one_file

   ! Local variables
   integer                                         :: i, j, k, nn, kkkk
   complex(kind = DP), allocatable, dimension(:,:) :: mat_square, mat_square2, &
                                                      mat_square3
   real(kind = DP), allocatable                    :: totatom(:)
   real(kind = DP)                                 :: tot
   logical                                         :: computes_full_product
   real(kind = DP), allocatable                    :: totplanes(:)
   integer                                         :: ierr
   integer                                         :: free_unit

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering dmft_show_dos'

   allocate(mat_square2(ngwf_basis%num, ngwf_basis%num), stat = ierr)
   call utils_alloc_check('dmft_show_dos', 'mat_square2', ierr)
   allocate(mat_square3(ngwf_basis%num, ngwf_basis%num), stat = ierr)
   call utils_alloc_check('dmft_show_dos', 'mat_square3', ierr)
   if(.not.pub_dmft_split)then
      call sparse_convert(mat_square2, greenf%mat)
   else
      mat_square2 = greenf%backup
   endif
   call sparse_convert(mat_square3, overlap)
   nn = ngwf_basis%num
   allocate(totatom(nn), stat = ierr)
   call utils_alloc_check('dmft_show_dos', 'totatom', ierr)

   if(pub_on_root)then

      computes_full_product = .false.

      if(computes_full_product)then
         allocate(mat_square(ngwf_basis%num, ngwf_basis%num), stat = ierr)
         call utils_alloc_check('dmft_show_dos', 'mat_square', ierr)

         ! ebl: removing GPU for now
         ! #ifdef GPU_SPEEDUP
         ! if(pub_dmft_use_gpu_onlyme)then
         !    if(pub_on_root .or. pub_dmft_split)then
         !       ! A = Trace(G * S)
         !       i = size(mat_square2, 1)
         !       j = size(mat_square2, 2)
         !       k = size(mat_square3, 2)
         !       write(stdout, *) 'call matmul cuda'
         !       !1)
         ! call matmulcuda_c(mat_square2, mat_square3, mat_square, i, j, k) !
         ! Using Magma, is fastest
         !       !2)
         ! ! call matmulcuda_c_cublas(mat_square2, mat_square3, mat_square, i,
         ! j, k) ! Using Cublas,
         !       !3)
         ! ! call matmul_my_cuda_c(mat_square2, mat_square3, mat_square, i, j,
         ! k)
         !       if(pub_debug_on_root) write(stdout, *) 'done'
         !    endif
         ! if(pub_total_num_procs > 1 .and. .not.pub_dmft_split) call
         ! comms_bcast(0, mat_square)
         ! else
         !    mat_square = MATMUL(mat_square2, mat_square3)
         ! endif
         ! #else
         mat_square = MATMUL(mat_square2, mat_square3)
         ! #endif
         tot = 0.0_DP
         do i = 1, nn
            totatom(i) = -aimag(mat_square(i, i))/PI
            tot        =  tot + totatom(i)
         enddo
         deallocate(mat_square, stat = ierr)
         call utils_dealloc_check('dmft_show_dos', 'mat_square', ierr)
      else
         tot = 0.0_DP
         do i = 1, nn
            totatom(i)  = - aimag(sum(mat_square2(i, :)*mat_square3(:, i))) / PI
            tot         =   tot + totatom(i)
         enddo
      endif

      if(nkpoints > 1)then
         totatom = totatom * norm_kpoints(ikpoint)
         tot    = tot     * norm_kpoints(ikpoint)
      endif

      free_unit = utils_unit()
      if(.not.pub_dmft_split)then
         if(pub_on_root)then
            !in eV
            if((ien == ffirst .and. one_k_one_file) .or. (.not.one_k_one_file &
                 .and. ien == ffirst .and. ikpoint == 1))then
               open(unit = free_unit, file = trim(pub_rootname) // &
                    '.dmft_dos_spin_' // trim(adjustl(utils_int_to_str(is))), &
                    status = 'replace')
            else
               open(unit = free_unit, file = trim(pub_rootname) // &
                    '.dmft_dos_spin_' // trim(adjustl(utils_int_to_str(is))), &
                    status = 'old', position = 'append')
            endif
            write(free_unit, '(i5, i3, i5, 300f14.6)') ien, 1, &
                 pub_dmft_points, (real(energy-(fermi_e(is) + &
                 pub_dmft_chem_shift)))*HARTREE_IN_EVS, tot
            close(free_unit)
         endif
      else
         do kkkk = 1, pub_total_num_procs
            if(kkkk == pub_my_proc_id + 1)then
               if(ien == ffirst .and. kkkk == 1 .and. ( (.not.pub_dmft_splitk &
                    .and. (ikpoint == 1 .or. one_k_one_file)) .or. &
                    (pub_dmft_splitk .and. ikpoint == splitk_first) ) )then
                  open(unit = free_unit, file = 'FULL_DOS_' // &
                       trim(adjustl(utils_int_to_str(is))) // '_rank' // &
                       trim(adjustl(utils_int_to_str(pub_dmft_mpi_rank + &
                       kp_file_shift))), status = 'replace')
               else
                  open(unit = free_unit, file = 'FULL_DOS_' // &
                       trim(adjustl(utils_int_to_str(is))) // '_rank' // &
                       trim(adjustl(utils_int_to_str(pub_dmft_mpi_rank + &
                       kp_file_shift))), position = 'append', status = 'old')
               endif
               write(free_unit, '(i5, i3, i5, 300f14.6)') ien, 1, &
                    pub_dmft_points, (real(energy-(fermi_e(is) + &
                    pub_dmft_chem_shift)))*HARTREE_IN_EVS, tot
               close(free_unit)
            endif
            call comms_barrier
         enddo
      endif

      if(nplanes /= 0)then
         allocate(totplanes(nplanes), stat = ierr)
         call utils_alloc_check('dmft_show_dos', 'totplanes', ierr)
         totplanes = 0.0_DP
         do j = 1, nplanes
            do i = 1, nn
               if(orb_in_plane(j, i)) totplanes(j) = totplanes(j) + totatom(i)
            enddo
         enddo

         do i = 1, nplanes
            if(.not.pub_dmft_split)then
               if(pub_on_root)then
                  if((ien == ffirst .and. one_k_one_file) .or. &
                       (.not.one_k_one_file .and. ien == ffirst .and. ikpoint &
                       == 1))then
                     open(unit = free_unit, file = 'FULL_DOS_PLANE_' // &
                          trim(adjustl(utils_int_to_str(i))) // '_' // &
                          trim(adjustl(utils_int_to_str(is))) // '_rank' // &
                          trim(adjustl(utils_int_to_str(pub_dmft_mpi_rank + &
                          kp_file_shift))), status = 'replace')
                  else
                     open(unit = free_unit, file = 'FULL_DOS_PLANE_' // &
                          trim(adjustl(utils_int_to_str(i))) // '_' // &
                          trim(adjustl(utils_int_to_str(is))) // '_rank' // &
                          trim(adjustl(utils_int_to_str(pub_dmft_mpi_rank + &
                          kp_file_shift))), status = 'old', position = &
                          'append')
                  endif
                  write(free_unit, '(i5, i3, i5, 300f14.6)') ien, 1, &
                       pub_dmft_points, (real(energy-(fermi_e(is) + &
                       pub_dmft_chem_shift)))*HARTREE_IN_EVS, totplanes(i)
                  close(free_unit)
               endif
            else
               do kkkk = 1, pub_total_num_procs
                  if(kkkk == pub_my_proc_id + 1)then
                     if(ien == ffirst .and. kkkk == 1 .and. ( &
                          (.not.pub_dmft_splitk .and. (ikpoint == 1 .or. &
                          one_k_one_file)) .or. (pub_dmft_splitk .and. ikpoint &
                          == splitk_first) ) )then
                        open(unit = free_unit, file = 'FULL_DOS_PLANE_' // &
                             trim(adjustl(utils_int_to_str(i))) // '_' // &
                             trim(adjustl(utils_int_to_str(is))) // '_rank' // &
                             trim(adjustl(utils_int_to_str(pub_dmft_mpi_rank + &
                             kp_file_shift))), status = 'replace')
                     else
                        open(unit = free_unit, file = 'FULL_DOS_PLANE_' // &
                             trim(adjustl(utils_int_to_str(i))) // '_' // &
                             trim(adjustl(utils_int_to_str(is))) // '_rank' // &
                             trim(adjustl(utils_int_to_str(pub_dmft_mpi_rank + &
                             kp_file_shift))), position = 'append', status = &
                             'old')
                     endif
                     write(free_unit, '(i5, i3, i5, 300f14.6)') ien, 1, &
                          pub_dmft_points, (real(energy-(fermi_e(is) + &
                          pub_dmft_chem_shift)))*HARTREE_IN_EVS, totplanes(i)
                     close(free_unit)
                  endif
                  call comms_barrier
               enddo

            endif
         enddo

         deallocate(totplanes, stat = ierr)
         call utils_dealloc_check('dmft_show_dos', 'totplanes', ierr)
      endif

      ! ebl: QC printing for comparing real-frequency DMFT calculations
      !      Only do this for the first 10 energy points
      if ((pub_print_qc) .and. (ien-ffirst < 10)) then
         call utils_qc_print("energy", real(energy-(fermi_e(is) + &
                 pub_dmft_chem_shift))*HARTREE_IN_EVS)
         call utils_qc_print("DOS(spin="//trim(adjustl(utils_int_to_str(is))) // ")", tot)
      end if

   endif

   deallocate(mat_square2, stat = ierr)
   call utils_dealloc_check('dmft_show_dos', 'mat_square2', ierr)
   deallocate(mat_square3, stat = ierr)
   call utils_dealloc_check('dmft_show_dos', 'mat_square3', ierr)
   deallocate(totatom, stat = ierr)
   call utils_dealloc_check('dmft_show_dos', 'totatom', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving dmft_show_dos'

 end subroutine

 subroutine dmft_invert_full_gpu_real(mat_inv, mat, greenfunc)

   !---------------------------------------------------------------------------!
   ! Explicit inversion of a sparse (real) matrix                              !
   !---------------------------------------------------------------------------!
   ! NOTES                                                                     !
   !   ebl: disabling GPU capabilities for now...                              !
   !---------------------------------------------------------------------------!

   use comms,  only: comms_bcast, pub_on_root, pub_total_num_procs
   ! #ifdef GPU_SPEEDUP
   ! use fortran_cuda, only : magma_fortran_double_singleprec_bypass
   ! #endif
   use rundat, only: pub_debug_on_root
   use sparse, only: sparse_convert, sparse_num_rows
   use utils,  only: utils_abort, utils_assert, utils_alloc_check, &
        utils_dealloc_check

   implicit none

   ! Arguments
   type(SPAM3), intent(inout) :: mat_inv
   type(SPAM3), intent(in   ) :: mat
   logical,     intent(in   ) :: greenfunc

   ! Local variables
   integer                                      :: nrows
   real(kind = DP), allocatable, dimension(:,:) :: mat_square
   integer                                      :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_invert_full_gpu_real'

   call utils_assert(greenfunc, "Error in invert_full_gpu_real: this routine &
        &is only intended for green functions")

   nrows = sparse_num_rows(mat)
   allocate(mat_square(nrows, nrows), stat = ierr)
   call utils_alloc_check('dmft_invert_full_gpu_real', 'mat_square', ierr)
   call sparse_convert(mat_square, mat)

   ! #ifdef GPU_SPEEDUP
   !     if(pub_on_root .or. pub_dmft_split)then
   !        ! Single precision
   !        call magma_fortran_double_singleprec_bypass(nn, mat_square)
   !     endif
   ! if(pub_total_num_procs > 1 .and. .not.pub_dmft_split) call comms_bcast(0,
   ! mat_square)
   ! #else
   call utils_abort("Error in invert_full_gpu_real: arrived in GPU function &
        &when GPUs are not being used")
   ! #endif
   call sparse_convert(mat_inv, mat_square)
   deallocate(mat_square, stat = ierr)
   call utils_dealloc_check('dmft_invert_full_gpu_real', 'mat_square', ierr)
   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_invert_full_gpu_real'

 end subroutine

 subroutine internal_invert_full_gpu(mat_inv, mat, overlap, energy, greenfunc)

   !---------------------------------------------------------------------------!
   ! Explicit inversion of a sparse matrix using GPUs                          !
   !---------------------------------------------------------------------------!

   use constants,    only: HARTREE_IN_EVS
   ! #ifdef GPU_SPEEDUP
   ! use fortran_cuda, only: diago_cuda_it_c, magma_fortran_comp_, &
   !      matinv_ge_complex, matinv_magma_complex
   ! #endif
   use sparse,       only: sparse_convert, sparse_num_rows
   use comms, only: comms_barrier, comms_bcast, pub_my_proc_id, pub_on_root, &
        pub_total_num_procs
   use utils, only: utils_abort, utils_alloc_check, utils_assert, &
        utils_dealloc_check, utils_unit
   use rundat, only: pub_debug, pub_debug_on_root, pub_dmft_points, &
        pub_dmft_temp

   implicit none

   ! Arguments
   type(SPAM3),        intent(in   ) :: mat
   type(SPAM3),        intent(inout) :: mat_inv
   type(SPAM3),        intent(in   ) :: overlap
   complex(kind = DP), intent(in   ) :: energy
   logical,            intent(in   ) :: greenfunc

   ! Local variables
   integer                         :: nn
   integer                         :: kkkk
   real(kind = DP)                 :: tot
   logical                         :: computes_full_product
   complex(kind = DP), allocatable :: mat_square(:,:)
   complex(kind = DP), allocatable :: overlap_square(:,:)
   real(kind = DP), allocatable    :: totatom(:)
   real(kind = DP), allocatable    :: totplanes(:)
   integer                         :: funit

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &internal_invert_full_gpu'

   call utils_abort("Error in internal_invert_full_gpu: should not be &
        &attempting to use this function without GPU_SPEEDUP turned on")
   ! #endif

   ! ebl: The following code has been tidied (indented and nominally cleaned)
   ! but still need to work on proper documentation of arguments/local variables
   !
   ! call utils_assert(greenfunc, "Error in invert_full_gpu: this routine is
   ! only intended for green functions")
   !
   !    allocate(mat_square(ngwf_basis%num, ngwf_basis%num), stat = ierr)
   !    call utils_alloc_check('internal_invert_full_gpu', 'mat_square', ierr)
   !    call sparse_convert(mat_square, mat)
   !
   !    nn = sparse_num_rows(mat)
   !    allocate(totatom(nn), stat = ierr)
   !    call utils_alloc_check('internal_invert_full_gpu', 'totatom', ierr)
   !
   ! if(pub_debug_on_root) write(stdout, *) 'SIZE OF MATRIX TO DIAGONALIZE :
   ! ', nn
   !    if(pub_debug_on_root) write(stdout, *) 'start diago NOW'
   !
   ! #ifdef GPU_SPEEDUP
   !    ! Double precision for tail, where Green func is small
   !    if(pub_on_root .or. pub_dmft_split)then
   ! if( internal_doubleprec(energy, chem, fermi_e, is, ien, nspecial_frequ) )
   ! then
   ! #ifdef GPU_HANGS
   !          ! call diago_cuda_it_c(nn, mat_square)   ! WORKS, 5.5-6sec
   !          call matinv_ge_complex(nn, mat_square)   ! WORKS, 2.5-3.5 sec
   ! #else
   ! call magma_fortran_comp_(nn, mat_square) ! WORKS, 2-3 sec (best), but
   ! more memory demanding
   ! #endif
   !       else
   !          ! Single precision
   !          call magma_fortran_comp_singleprec_bypass(nn, mat_square)
   !       endif
   !    endif
   ! if(pub_total_num_procs > 1 .and. .not.pub_dmft_split) call comms_bcast(0,
   ! mat_square)
   ! #endif
   !
   !    if(pub_debug_on_root .and. ien == 1) write(stdout, *) 'done'
   !
   !    ! Density of states
   !    if(pub_dmft_temp < 0.0_DP) then
   !       allocate(overlap_square(ngwf_basis%num, ngwf_basis%num), stat = ierr)
   !       call utils_alloc_check('internal_invert_full_gpu', 'overlap_square', &
   !            ierr)
   !       call sparse_convert(overlap_square, overlap)
   !
   !       if (pub_on_root) then
   !          tot = 0.0_DP
   !          do i = 1, nn
   ! totatom(i) = -aimag(sum(mat_square(i, :)*overlap_square(:, i)))/pi
   !             tot        =  tot + totatom(i)
   !          enddo
   !
   !          if (nkpoints > 1) then
   !            totatom = totatom * norm_kpoints(ikpoint)
   !            tot     = tot     * norm_kpoints(ikpoint)
   !          endif
   !
   !          if (.not.pub_dmft_split) then
   !             if (pub_on_root) then
   !                !in eV
   !                funit = utils_unit()
   ! if ((ien == ffirst .and. one_k_one_file) .or. (.not. one_k_one_file .and.
   ! ien == ffirst .and. ikpoint == 1)) then
   ! open(unit = funit, file = 'FULL_DOS_' // TRIM(ADJUSTL( &
   ! utils_int_to_str(is))) // '_rank' // TRIM(ADJUSTL(utils_int_to_str( &
   ! pub_dmft_mpi_rank + kp_file_shift))), status = 'replace')
   !                else
   ! open(unit = funit, file = 'FULL_DOS_' // TRIM(ADJUSTL( &
   ! utils_int_to_str(is))) // '_rank' // TRIM(ADJUSTL(utils_int_to_str( &
   ! pub_dmft_mpi_rank + kp_file_shift))), status = 'old', position = 'append')
   !                endif
   ! write(funit, '(i5, i3, i5, 300f14.6)') ien, 1, pub_dmft_points,
   ! (real(energy-(fermi_e(is) + pub_dmft_chem_shift))) *HARTREE_IN_EVS, tot
   !                close(funit)
   !             endif
   !          else
   !             funit = utils_unit()
   !             do kkkk = 1, pub_total_num_procs
   !                if (kkkk == pub_my_proc_id + 1) then
   ! if ( ien == ffirst .and. kkkk == 1 .and. ((.not.pub_dmft_splitk .and.
   ! (ikpoint == 1 .or. one_k_one_file)) .or. (pub_dmft_splitk .and. ikpoint
   ! == splitk_first))) then
   ! open(unit = funit, file = 'FULL_DOS_' // TRIM(ADJUSTL( &
   ! utils_int_to_str(is))) // '_rank' // TRIM(ADJUSTL(utils_int_to_str( &
   ! pub_dmft_mpi_rank + kp_file_shift))), status = 'replace')
   !                   else
   ! open(unit = funit, file = 'FULL_DOS_' // TRIM(ADJUSTL( &
   ! utils_int_to_str(is))) // '_rank' // TRIM(ADJUSTL(utils_int_to_str( &
   ! pub_dmft_mpi_rank + kp_file_shift))), position = 'append', status = 'old')
   !                   endif
   ! write(funit, '(i5, i3, i5, 300f14.6)') ien, 1, pub_dmft_points,
   ! (real(energy-(fermi_e(is) + pub_dmft_chem_shift))) *HARTREE_IN_EVS, tot
   !                   close(funit)
   !                endif
   !                call comms_barrier
   !             enddo
   !          endif
   !       endif
   !       deallocate(overlap_square, stat = ierr)
   !       call utils_dealloc_check('internal_invert_full_gpu', 'overlap_square', &
   !            ierr)
   !
   !       if(nplanes /= 0)then
   !          allocate(totplanes(nplanes), stat = ierr)
   !          call utils_dealloc_check('internal_invert_full_gpu', 'totplanes', &
   !               ierr)
   !          totplanes = 0.0_DP
   !
   !          do j = 1, nplanes
   !             do i = 1, nn
   ! if(orb_in_plane(j, i)) totplanes(j) = totplanes(j) + totatom(i)
   !             enddo
   !          enddo
   !
   !          if(.not.pub_dmft_split)then
   !             do i = 1, nplanes
   !                funit = utils_unit()
   !                if(pub_on_root)then
   ! if((ien == ffirst .and. one_k_one_file) .or. (.not.one_k_one_file .and.
   ! ien == ffirst .and. ikpoint == 1)) then
   ! open(unit = funit, file = 'FULL_DOS_PLANE_' // TRIM(ADJUSTL( &
   ! utils_int_to_str(i))) // '_' // TRIM(ADJUSTL(utils_int_to_str(is))) // &
   ! '_rank' // TRIM(ADJUSTL(utils_int_to_str(pub_dmft_mpi_rank + &
   ! kp_file_shift))), status = 'replace')
   !                   else
   ! open(unit = funit, file = 'FULL_DOS_PLANE_' // TRIM(ADJUSTL( &
   ! utils_int_to_str(i))) // '_' // TRIM(ADJUSTL(utils_int_to_str(is))) &
   ! // '_rank' // TRIM(ADJUSTL(utils_int_to_str(pub_dmft_mpi_rank + &
   ! kp_file_shift))), status = 'old', position = 'append')
   !                   endif
   ! write(funit, '(i5, i3, i5, 300f14.6)') ien, 1, pub_dmft_points,
   ! (real(energy-(fermi_e(is) + pub_dmft_chem_shift)))*HARTREE_IN_EVS,
   ! totplanes(i)
   !                close(funit)
   !             endif
   !          enddo
   !          else
   !             do i = 1, nplanes
   !                do kkkk = 1, pub_total_num_procs
   !                   if(kkkk == pub_my_proc_id + 1)then
   !                      funit = utils_unit()
   ! if(ien == ffirst .and. kkkk == 1 .and. ((.not.pub_dmft_splitk .and.
   ! (ikpoint == 1 .or. one_k_one_file)) .or. (pub_dmft_splitk .and. ikpoint
   ! == splitk_first))) then
   ! open(unit = funit, file = 'FULL_DOS_PLANE_' // TRIM(ADJUSTL( &
   ! utils_int_to_str(i)))// '_' // TRIM(ADJUSTL(utils_int_to_str(is))) // &
   ! '_rank' // TRIM(ADJUSTL(utils_int_to_str(pub_dmft_mpi_rank + &
   ! kp_file_shift))), status = 'replace')
   !                      else
   ! open(unit = funit, file = 'FULL_DOS_PLANE_' // TRIM(ADJUSTL( &
   ! utils_int_to_str(i))) // '_' // TRIM(ADJUSTL(utils_int_to_str(is))) // &
   ! '_rank' // TRIM(ADJUSTL(utils_int_to_str(pub_dmft_mpi_rank + &
   ! kp_file_shift))), & position = 'append', status = 'old')
   !                      endif
   ! write(funit, '(i5, i3, i5, 300f14.6)') ien, 1, pub_dmft_points,
   ! (real(energy-(fermi_e(is) + pub_dmft_chem_shift)))*HARTREE_IN_EVS,
   ! totplanes(i)
   !                      close(funit)
   !                   endif
   !                   call comms_barrier
   !                enddo
   !             enddo
   !          endif
   !          deallocate(totplanes, stat=ierr)
   !          call utils_dealloc_check('internal_invert_full_gpu', 'totplanes', &
   !               ierr)
   !       endif
   !    endif
   !    ! DOS calculation complete
   !
   !    greenf%backupflag = .true.
   !    greenf%backup%zmtx = mat_square
   !
   !    call sparse_convert(mat_inv, mat_square)
   !    deallocate(mat_square, stat = ierr)
   !    call utils_dealloc_check('internal_invert_full_gpu', 'mat_square', ierr)
   !    deallocate(totatom, stat = ierr)
   !    call utils_dealloc_check('internal_invert_full_gpu', 'totatom', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &internal_invert_full_gpu'

 end subroutine

 subroutine dmft_matmul_in_dense_format(matc, mata, matb, aba, show_trace)

   !---------------------------------------------------------------------------!
   ! Calculates C = A*B, or C = A*B*A if flag ABA is present                   !
   !---------------------------------------------------------------------------!

   ! #ifdef GPU_SPEEDUP
   ! use fortran_cuda, only: matmulcuda_c, matmulcuda_c_cublas, &
   !      matmulcuda_c_singleprec_bypass
   ! #endif
   use sparse, only: sparse_convert, sparse_num_cols, sparse_num_rows, &
        sparse_show_matrix
   use comms,        only: comms_bcast, pub_on_root, pub_total_num_procs
   use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
        utils_unit
   use ion,          only: ELEMENT
   use rundat,       only: pub_debug, pub_debug_on_root

   implicit none

   ! Arguments
   type(SPAM3),       intent(in   ) :: mata
   type(SPAM3),       intent(in   ) :: matb
   type(SPAM3),       intent(inout) :: matc
   logical, optional, intent(in   ) :: aba
   logical, optional, intent(in   ) :: show_trace

   ! Local variables
   integer                         :: i, n1, n2, n3
   complex(kind = DP), allocatable :: matdensa(:,:)
   complex(kind = DP), allocatable :: matdensb(:,:)
   complex(kind = DP), allocatable :: matdensc(:,:)
   complex(kind = DP)              :: traceit
   integer                         :: funit
   integer                         :: ierr
   character(2), allocatable       :: labels(:)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_matmul_in_dense_format'

   n1 = sparse_num_rows(mata)  ! row A
   n2 = sparse_num_cols(mata)  ! sum
   n3 = sparse_num_cols(matb)  ! col B

   allocate(matdensa(n1, n2), stat = ierr)
   call utils_alloc_check('dmft_matmul_in_dense_format', 'matdensa', ierr)
   allocate(matdensb(n2, n3), stat = ierr)
   call utils_alloc_check('dmft_matmul_in_dense_format', 'matdensb', ierr)
   allocate(matdensc(n1, n3), stat = ierr)
   call utils_alloc_check('dmft_matmul_in_dense_format', 'matdensc', ierr)

   ! Convert matrix to full square
   call sparse_convert(matdensa, mata)
   call sparse_convert(matdensb, matb)

   ! if(pub_debug_on_root .and. (pub_output_detail >= VERBOSE))then
   !    write(stdout, *) '------------------------------------'
   !    write(stdout, *) 'row A     : ', sparse_num_rows(mata)
   !    write(stdout, *) 'col A     : ', sparse_num_cols(mata)
   !    write(stdout, *) 'row B     : ', sparse_num_rows(matb)
   !    write(stdout, *) 'col B     : ', sparse_num_cols(matb)
   !    write(stdout, *) 'row C     : ', sparse_num_rows(matc)
   !    write(stdout, *) 'col C     : ', sparse_num_cols(matc)
   !    write(stdout, *) 'shape A   : ', shape(matdensa)
   !    write(stdout, *) 'shape B   : ', shape(matdensb)
   !    write(stdout, *) 'shape C   : ', shape(matdensc)
   !    write(stdout, *) 'n1n2n3    : ', n1, n2, n3
   !    write(stdout, *) '------------------------------------'
   ! endif

   ! ebl: disabling GPU for now...
   ! #ifdef GPU_SPEEDUP
   ! if(pub_dmft_use_gpu_onlyme)then
   !    if(pub_on_root .or. pub_dmft_split)then
   ! if( internal_doubleprec(energy, chem, fermi_e, is, ien, nspecial_frequ) )
   ! then
   !          if(pub_debug_on_root) write(stdout, *) 'matmucuda_c'
   ! call matmulcuda_c(matdensa, matdensb, matdensc, n1, n2, n3) ! fastest
   ! ! call matmulcuda_c_cublas(matdensa, matdensb, matdensc, n1, n2, n3) !
   ! can be a factor 2 slower
   !       else
   ! if(pub_debug_on_root .and. (pub_output_detail >= VERBOSE)) write(stdout,
   ! *) 'matmucuda_c_singleprec_bypass'
   ! call matmulcuda_c_singleprec_bypass(matdensa, matdensb, matdensc, n1, n2,
   ! n3)
   !       endif
   !       if(present(aba))then
   !          if(n1 /= n2 .or. n1 /= n3)then
   ! call utils_abort("Error in dmft_matmul_in_dense_format: tried to
   ! calculate product A*B*A with rectangular matrices")
   !          endif
   !          call matmulcuda_c(matdensc, matdensa, matdensb, n1, n2, n3)
   !          matdensc = matdensb
   !       endif
   !    endif
   ! if(pub_total_num_procs > 1 .and. .not.pub_dmft_split) call comms_bcast(0,
   ! matdensc)
   ! else
   !    matdensc = matmul(matdensa, matdensb)
   !    if(present(aba))then
   !       matdensc = matmul(matdensc, matdensa)
   !    endif
   ! endif
   ! #else

   matdensc = matmul(matdensa, matdensb)
   if(present(aba))then
      matdensc = matmul(matdensc, matdensa)
   endif

   ! #endif

   if(present(show_trace))then
      if(show_trace)then

         if(allocated(labels)) then
            deallocate(labels, stat = ierr)
            call utils_dealloc_check('dmft_matmul_in_dense_format', 'labels', &
                 ierr)
         end if

         allocate(labels(size(matdensc, 1)), stat = ierr)
         call utils_alloc_check('dmft_matmul_in_dense_format', 'labels', ierr)

         funit = utils_unit()
         call sparse_show_matrix(matc, outunit = funit, show_elems = .true., &
              labels = labels)

         if(sparse_num_cols(matc) /= size(matdensc, 1))then
            write(stdout, *) 'number of labels : ', sparse_num_cols(matc)
            write(stdout, *) 'size of matrix   : ', size(matdensc, 1)
            call utils_abort("Error in dmft_matmul_in_dense_format: matrix &
                 &size mismatch")
         endif

         traceit = 0.0_DP

         if(pub_on_root .and. pub_dmft_mpi_rank == 1)then
            write(stdout, *) 'TOTAL NUMBER OF ELEMENTS : ', size(matdensc, 1)
            do i = 1, size(matdensc, 1)
               traceit = traceit + matdensc(i, i)
               write(stdout, '(a, i5, 4f20.12)') labels(i), i, matdensc(i, i)
            enddo
            write(stdout, '(a, 4f20.12)') 'TRACE IS : ', traceit
         endif

      endif
   endif

   ! Convert matrix back to SPAM3
   call sparse_convert(matc, matdensc)
   deallocate(matdensa, stat = ierr)
   call utils_dealloc_check('dmft_matmul_in_dense_format', 'matdensa', ierr)
   deallocate(matdensb, stat = ierr)
   call utils_dealloc_check('dmft_matmul_in_dense_format', 'matdensb', ierr)
   deallocate(matdensc, stat = ierr)
   call utils_dealloc_check('dmft_matmul_in_dense_format', 'matdensc', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_matmul_in_dense_format'

 end subroutine

 subroutine dmft_matmul_in_dense_format_r(matc, mata, matb, aba, show_trace)

   !---------------------------------------------------------------------------!
   ! Calculates C = A*B, or C = A*B*A if flag ABA is present                   !
   !---------------------------------------------------------------------------!

   ! #ifdef GPU_SPEEDUP
   ! use fortran_cuda, only: matmul_my_cuda, matmulcuda_r_cublas
   ! #endif
   use sparse, only: sparse_convert, sparse_num_cols, sparse_num_rows, &
        sparse_show_matrix
   use comms,        only: comms_bcast, pub_on_root, pub_total_num_procs
   use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
        utils_unit
   use ion,          only: ELEMENT
   use rundat,       only: pub_debug, pub_debug_on_root

   implicit none

   ! Arguments
   type(SPAM3),       intent(in   ) :: mata
   type(SPAM3),       intent(in   ) :: matb
   type(SPAM3),       intent(inout) :: matc
   logical, optional, intent(in   ) :: aba
   logical, optional, intent(in   ) :: show_trace

   ! Local variables
   integer                      :: i, n1, n2, n3
   real(kind = DP), allocatable :: matdensa(:,:)
   real(kind = DP), allocatable :: matdensb(:,:)
   real(kind = DP), allocatable :: matdensc(:,:)
   real(kind = DP)              :: traceit
   integer                      :: ierr
   integer                      :: funit
   character(2), allocatable    :: labels(:)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_matmul_in_dense_format_r'

   n1 = sparse_num_rows(mata)   !row A
   n2 = sparse_num_cols(mata)   !sum
   n3 = sparse_num_cols(matb)   !col B

#ifdef debugMIXv2
   if(sparse_num_rows(matb) /= n2)then
      write(stdout, *) 'rows b, col a : ', sparse_num_rows(matb), n2
      call utils_abort("Error in dmft_matmul_in_dense_format_r: dimensions do &
           &not match")
   endif
   if(sparse_num_rows(matc) /= n1)then
      write(stdout, *) 'rows c, rows a : ', sparse_num_rows(matc), n1
      call utils_abort("Error in dmft_matmul_in_dense_format_r: dimensions do &
           &not match")
   endif
   if(sparse_num_cols(matc) /= n3)then
      write(stdout, *) 'cols c, cols b : ', sparse_num_cols(matc), n3
      call utils_abort("Error in dmft_matmul_in_dense_format_r: dimensions do &
           &not match")
   endif
#endif

   allocate(matdensa(n1, n2), stat = ierr)
   call utils_alloc_check('dmft_matmul_in_dense_format_r', 'matdensa', ierr)
   allocate(matdensb(n2, n3), stat = ierr)
   call utils_alloc_check('dmft_matmul_in_dense_format_r', 'matdensb', ierr)
   allocate(matdensc(n1, n3), stat = ierr)
   call utils_alloc_check('dmft_matmul_in_dense_format_r', 'matdensc', ierr)

   call sparse_convert(matdensa, mata)
   call sparse_convert(matdensb, matb)

   if(pub_debug_on_root .and. (pub_output_detail >= VERBOSE))then
      write(stdout, *) '------------------------------------'
      write(stdout, *) 'row A     : ', sparse_num_rows(mata)
      write(stdout, *) 'col A     : ', sparse_num_cols(mata)
      write(stdout, *) 'row B     : ', sparse_num_rows(matb)
      write(stdout, *) 'col B     : ', sparse_num_cols(matb)
      write(stdout, *) 'row C     : ', sparse_num_rows(matc)
      write(stdout, *) 'col C     : ', sparse_num_cols(matc)
      write(stdout, *) 'shape A   : ', shape(matdensa)
      write(stdout, *) 'shape B   : ', shape(matdensb)
      write(stdout, *) 'shape C   : ', shape(matdensc)
      write(stdout, *) 'n1n2n3    : ', n1, n2, n3
      write(stdout, *) '------------------------------------'
   endif

#ifdef debugMIXv3
   if(pub_on_root)then
      write(stdout, *)
      write(stdout, '(a, 40f15.3)') 'some elements of matrix a : ', &
           matdensa(1, 1:min(20, size(matdensa, 2)))
      write(stdout, *)
      write(stdout, '(a, 40f15.3)') 'some elements of matrix b : ', &
           matdensa(1, 1:min(20, size(matdensa, 2)))
   endif
#endif

   ! ebl: disabling GPU for now...
   ! #ifdef GPU_SPEEDUP
   ! if(pub_dmft_use_gpu_onlyme)then
   !    if(pub_on_root .or. pub_dmft_split)then
   !       call matmulcuda_r(matdensa, matdensb, matdensc, n1, n2, n3)
   !       ! call matmulcuda_r_cublas(matdensa, matdensb, matdensc, n1, n2, n3)
   !       if(pub_debug_on_root)write(stdout, *) 'done/matmulcuda_r'
   !       ! call matmul_my_cuda_r(matdensa, matdensb, matdensc, n1, n2, n3)
   !       if(present(aba))then
   !          if(n1 /= n2 .or. n1 /= n3)then
   ! call utils_abort("Error in dmft_matmul_in_dense_format_r: tried to
   ! calculate product A*B*A with rectangular matrices")
   !          endif
   !          call matmulcuda_r_cublas(matdensc, matdensa, matdensb, n1, n2, n3)
   !          matdensc = matdensb
   !       endif
   !    endif
   ! if(pub_total_num_procs > 1 .and. .not.pub_dmft_split) call comms_bcast(0,
   ! matdensc)
   ! else
   !    matdensc = matmul(matdensa, matdensb)
   !    if(present(aba))then
   !       matdensc = matmul(matdensc, matdensa)
   !    endif
   ! endif
   ! #else
   matdensc = matmul(matdensa, matdensb)
   if(present(aba))then
      matdensc = matmul(matdensc, matdensa)
   endif
   ! #endif

   if(present(show_trace))then
      if(show_trace)then
         if (allocated(labels)) then
            deallocate(labels, stat=ierr)
            call utils_dealloc_check('dmft_matmul_in_dense_format_r', &
                 'labels', ierr)
         end if
         allocate(labels(size(matdensc, 1)), stat = ierr)
         call utils_alloc_check('dmft_matmul_in_dense_format_r', 'labels', ierr)
         if(sparse_num_cols(matc) /= size(matdensc, 1))then
            write(stdout, *) 'number of labels : ', sparse_num_cols(matc)
            write(stdout, *) 'size of matrix   : ', size(matdensc, 1)
            call utils_abort("Error in dmft_matmul_in_dense_format_r: matrix &
                 &size mismatch")
         endif

         funit = utils_unit()
         call sparse_show_matrix(matc, outunit = funit, show_elems = .true., &
              labels = labels)
         traceit = 0.0_DP
         if(pub_on_root .and. pub_dmft_mpi_rank == 1)then
            write(stdout, *) 'TOTAL NUMBER OF ELEMENTS : ', size(matdensc, 1)
            do i = 1, size(matdensc, 1)
               traceit = traceit + matdensc(i, i)
               write(stdout, '(a, i5, 4f20.12)') labels(i), i, matdensc(i, i)
            enddo
            write(stdout, '(a, 4f20.12)') 'TRACE IS : ', traceit
         endif
      endif
   endif

   call sparse_convert(matc, matdensc)
   deallocate(matdensa, stat = ierr)
   call utils_dealloc_check('dmft_matmul_in_dense_format_r', 'matdensa', ierr)
   deallocate(matdensb, stat = ierr)
   call utils_dealloc_check('dmft_matmul_in_dense_format_r', 'matdensb', ierr)
   deallocate(matdensc, stat = ierr)
   call utils_dealloc_check('dmft_matmul_in_dense_format_r', 'matdensc', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_matmul_in_dense_format_r'

 end subroutine


 subroutine dmft_hubbard_greenf_info(loop, hub_atom, hub_atomb, n_, m_, &
      h_atoms, hub, hub_proj_basis, energy, fermi_e, num_spins, is, ien, &
      ikpoint, norm_kpoints, ffirst, nkpoints, kp_file_shift, splitk_first, &
      site_greenf_buffer, cluster, merge_neigh, merge_table, nmerge, &
      rot_vec_angles, rot_vec_angles_split, one_k_one_file)

   !---------------------------------------------------------------------------!
   ! Writes the Green's function to file                                       !
   !---------------------------------------------------------------------------!

   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use hubbard_init,      only: h_species
   use comms, only: comms_reduce, pub_my_proc_id, pub_on_root, &
        pub_total_num_procs
   use parallel_strategy, only: par => pub_par
   use rundat, only: pub_debug, pub_debug_on_root, &
        pub_dmft_points, pub_dmft_temp, pub_dmft_write, pub_rootname
   use utils,             only: utils_abort, utils_assert, utils_unit, &
        utils_int_to_str, utils_alloc_check, utils_dealloc_check

   implicit none

   ! Arguments
   type(DMFT_LOOP_INFO),            intent(in   ) :: loop
   integer,                         intent(in   ) :: hub_atom
   integer,                         intent(in   ) :: hub_atomb
   integer,                         intent(in   ) :: n_
   integer,                         intent(in   ) :: m_
   type(ARRAY_OF_MATRICES),         intent(in   ) :: h_atoms(par%nat_hub)
   type(HUBBARD_MODEL),             intent(in   ) :: hub
   type(FUNC_BASIS),                intent(in   ) :: hub_proj_basis
   complex(kind = DP),              intent(in   ) :: energy
   real(kind = DP),                 intent(in   ) :: fermi_e(2)
   integer,                         intent(in   ) :: num_spins
   integer,                         intent(in   ) :: is
   integer,                         intent(in   ) :: ien
   integer,                         intent(in   ) :: ikpoint
   real(kind = DP), allocatable,    intent(in   ) :: norm_kpoints(:)
   integer,                         intent(in   ) :: ffirst
   integer,                         intent(in   ) :: nkpoints
   integer,                         intent(in   ) :: kp_file_shift
   integer,                         intent(in   ) :: splitk_first
   complex(kind = DP), allocatable, intent(in   ) :: site_greenf_buffer(:,:,:,:)
   logical,                         intent(in   ) :: cluster
   integer, allocatable,            intent(in   ) :: merge_neigh(:)
   integer, allocatable,            intent(in   ) :: merge_table(:,:)
   integer,                         intent(in   ) :: nmerge
   real(kind = DP), allocatable,    intent(in   ) :: rot_vec_angles(:,:,:)
   real(kind = DP), allocatable,    intent(in   ) :: rot_vec_angles_split(:,:,:)
   logical,                         intent(in   ) :: one_k_one_file

   ! Local variables
   integer                         :: channels, channelsb
   complex(kind = DP)              :: ttt(n_, m_)
   complex(kind = DP), allocatable :: ttt_buffer(:,:,:)
   integer, allocatable            :: fre_buffer(:), at_buffer(:)
   complex(kind = DP), allocatable :: ener_buffer(:)
   real(kind = DP), allocatable    :: ener_bufferR(:), ener_bufferI(:)
   integer                         :: species, speciesb
   integer                         :: jj, kkkk, sizbuffer
   integer                         :: kk
   integer                         :: funit, funit2
   logical                         :: check
   integer                         :: l1, l2
   character(100)                  :: filename
   integer                         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_hubbard_greenf_info'

   if (.not. pub_dmft_write) then
      ! Do not generate Green's function files if dmft_write is set to FALSE
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_hubbard_greenf_info'
      return
   end if

   species = loop%sp
   species = loop%spb
   channels = loop%channels
   channelsb = loop%channelsb

   if(pub_dmft_split)then
      allocate(ttt_buffer(pub_total_num_procs, n_, m_), stat = ierr)
      call utils_alloc_check('dmft_hubbard_greenf_info', 'ttt_buffer', ierr)
      allocate(fre_buffer(pub_total_num_procs), stat = ierr)
      call utils_alloc_check('dmft_hubbard_greenf_info', 'fre_buffer', ierr)
      allocate(at_buffer(pub_total_num_procs), stat = ierr)
      call utils_alloc_check('dmft_hubbard_greenf_info', 'at_buffer', ierr)
      allocate(ener_buffer(pub_total_num_procs), stat = ierr)
      call utils_alloc_check('dmft_hubbard_greenf_info', 'ener_buffer', ierr)
      allocate(ener_bufferR(pub_total_num_procs), stat = ierr)
      call utils_alloc_check('dmft_hubbard_greenf_info', 'ener_bufferR', ierr)
      allocate(ener_bufferI(pub_total_num_procs), stat = ierr)
      call utils_alloc_check('dmft_hubbard_greenf_info', 'ener_bufferI', ierr)
      ener_buffer = 0
      ener_bufferR = 0.0_DP
      ener_bufferI = 0.0_DP
      fre_buffer = 0
      ttt_buffer = 0
      at_buffer = 0
   endif

   channels = dmft_totn_atom_merge(hub_atom, hub, hub_proj_basis, merge_neigh, &
        merge_table, nmerge, cluster)
   channelsb = dmft_totn_atom_merge(hub_atomb, hub, hub_proj_basis, &
        merge_neigh, merge_table, nmerge, cluster)

   if(channels /= size(ttt, 1) .or. channelsb /= size(ttt, 2))then
      write(stdout, *) 'ERROR shape of ttt is not matching channels'
      write(stdout, *) 'channels, channelsb : ', channels, channelsb
      write(stdout, *) 'shape ttt          : ', shape(ttt)
      call utils_abort("Mismatching shapes of arrays in dmft_hubbard_greenf_info")
   endif

   ttt = dmft_merged_green_block(hub_atom, hub_atomb, channels, channelsb, &
        site_greenf_buffer, hub, hub_proj_basis, cluster, merge_neigh, &
        merge_table, nmerge)

   filename = 'green_output' // TRIM(ADJUSTL( utils_int_to_str(loop%hat))) // &
        "_" // TRIM(ADJUSTL(utils_int_to_str(species))) // "_" // TRIM( &
        ADJUSTL(utils_int_to_str(is)))
   if(hub_atom /= hub_atomb)then
      ! if(pub_debug_on_root)write(stdout, *) 'dimer with atoms : ', hub_atom, &
      !      hub_atomb
      filename = TRIM(ADJUSTL(filename)) // "_dimer"
   endif

   if(pub_dmft_mpi_size > 1 .or. nkpoints > 1)then
      filename = TRIM(ADJUSTL(filename)) // "_rank" // &
           TRIM(ADJUSTL(utils_int_to_str(pub_dmft_mpi_rank + kp_file_shift)))
   endif

   if(ien == ffirst .and. ( (.not.pub_dmft_splitk .and. (ikpoint == 1 .or. &
        one_k_one_file)) .or. (pub_dmft_splitk .and. ikpoint == splitk_first) &
        ))then
      if(pub_on_root .or. .not.pub_dmft_split) then
         ! if(pub_debug_on_root)then
         !    write(stdout, *) 'RANK 0 OPENING GREEN_OUTPUT FILE : ', &
         !         trim(adjustl(filename))
         !    write(stdout, *) 'CHEM POTENTIAL                   : ', fermi_e(is)
         !    write(stdout, *) 'FREQUENCY                        : ', energy
         !    write(stdout, *) 'FREQUENCY LABEL  / FIRST         : ', ien, ffirst
         !    write(stdout, *) 'KPOINT                           : ', ikpoint
         ! endif
         funit = utils_unit()
         open(unit = funit, file = trim(adjustl(filename)), form = &
              'unformatted', STATUS = 'REPLACE')
      endif
   else
      if (pub_debug) then
         if(pub_on_root .or. .not.pub_dmft_split) then
            inquire(file = trim(adjustl(filename)), exist = check)
            call utils_assert(check, 'Error in dmft_hubbard_greenf_info: ' // &
                 trim(adjustl(filename)) // ' not found')
         end if
      end if
      if(pub_on_root .or. .not.pub_dmft_split) then
         funit = utils_unit()
         open(unit = funit, file = trim(adjustl(filename)), form = &
              'unformatted', position = 'append', status = 'old')
         ! if(pub_debug_on_root)then
         !    write(stdout, *) 'APPENDING INTO FILE : ', trim(adjustl(filename))
         !    write(stdout, *) 'FREQUENCY NUMBER    : ', ien
         !    write(stdout, *) 'KPOINT              : ', ikpoint
         ! endif
      endif
   endif

   if(nkpoints > 1)then
      ttt = ttt * norm_kpoints(ikpoint)
   endif

   if(pub_dmft_split)then
      ttt_buffer(pub_my_proc_id + 1, :, :) = ttt
      fre_buffer(pub_my_proc_id + 1) = ien
      at_buffer(pub_my_proc_id + 1) = loop%hat
      ener_bufferR(pub_my_proc_id + 1) = real(energy)
      ener_bufferI(pub_my_proc_id + 1) = aimag(energy)
      call comms_reduce('SUM', ttt_buffer(:, :, :))
      call comms_reduce('SUM', fre_buffer)
      call comms_reduce('SUM', at_buffer)
      call comms_reduce('SUM', ener_bufferR)
      call comms_reduce('SUM', ener_bufferI)
      ener_buffer = ener_bufferR + (0.0_DP, 1.0_DP)*ener_bufferI
      if(ANY( abs(at_buffer - at_buffer(1)) > 0 ))then
         write(stdout, *) 'atoms on proc : ', at_buffer
         call utils_abort("Error in dmft_hubbard_greenf_info: pub_dmft_split &
              &is True and procs are not doing same atom")
      endif
   endif

   if(pub_dmft_split)then
      sizbuffer = pub_total_num_procs
   else
      sizbuffer = 1
   endif

   !---------------------!

   l1 = dmft_totn_atom_merge(hub_atom, hub, hub_proj_basis, merge_neigh, &
        merge_table, nmerge, cluster)
   l2 = dmft_totn_atom_merge(hub_atom, hub, hub_proj_basis, merge_neigh, &
        merge_table, nmerge, cluster)

   do kkkk = 1, sizbuffer
      if(pub_dmft_split)then
         if(pub_on_root) write(funit) fre_buffer(kkkk), pub_dmft_points, &
              pub_dmft_temp, channels, channelsb, loop%hat, loop%hatb
      else
         write(funit) ien, pub_dmft_points, pub_dmft_temp, channels, &
              channelsb, loop%hat, loop%hatb
      endif

      if(num_spins == 1)then !paramagnetic
         if ((pub_dmft_split .and. pub_on_root) .or. .not. pub_dmft_split) then
            write(funit) dmft_merged_occupancy_block(hub_atom, 1, l1, h_atoms, &
                 hub, hub_proj_basis, merge_neigh, merge_table, nmerge, &
                 rot_vec_angles, rot_vec_angles_split, cluster), &
                 dmft_merged_occupancy_block(hub_atom, 1, l2, h_atoms, hub, &
                 hub_proj_basis, merge_neigh, merge_table, nmerge, &
                 rot_vec_angles, rot_vec_angles_split, cluster)
         endif
      elseif(num_spins == 2)then !spin up and dn
         if ((pub_dmft_split .and. pub_on_root) .or. .not. pub_dmft_split) then
            write(funit) dmft_merged_occupancy_block(hub_atom, 1, l1, h_atoms, &
                 hub, hub_proj_basis, merge_neigh, merge_table, nmerge, &
                 rot_vec_angles, rot_vec_angles_split, cluster), &
                 dmft_merged_occupancy_block(hub_atom, 2, l2, h_atoms, hub, &
                 hub_proj_basis, merge_neigh, merge_table, nmerge, &
                 rot_vec_angles, rot_vec_angles_split, cluster)
         endif
      else
         call utils_abort("Error in dmft_hubbard_greenf_info: spin shoould be &
              &1 or 2")
      endif

      if(ien == 1) then
         funit2 = utils_unit()
         open(unit = funit2, file = trim(pub_rootname) // &
              '.greenf_summary', position = "append")
         write(funit2, '(/a)')'########################################'
         write(funit2, '(a, i6, a, a)') 'DFT + U Greens function of atom ', &
              loop%hat, ' of Hubbard species ', h_species(species)%hub_species
         write(funit2, '(a,2i4)') ' hub_atom, hub_atomb : ', hub_atom, hub_atomb
         write(funit2, '(a,2i4)') ' hat and theatom : ', loop%hat, loop%theatom
         write(funit2, '(a,i4)') ' hat                           : ', loop%hat
         write(funit2, '(a,i4)') ' species                       : ', species
         write(funit2, '(a)')'########################################'
      endif

      if(pub_dmft_split)then
         if(pub_on_root)then
            if (pub_debug) then
               write(stdout, *) 'WRITING TO FILE : ', (fermi_e(is) + &
                    pub_dmft_chem_shift), ener_buffer(kkkk)
               write(stdout, *) 'WRITING TO FILE : ', ttt_buffer(kkkk, &
                    1:channels, 1:channelsb)
            end if
            write(funit) (fermi_e(is) + pub_dmft_chem_shift), &
                 ener_buffer(kkkk), ttt_buffer(kkkk, 1:channels, 1:channelsb)
         endif
         ttt = ttt_buffer(kkkk, 1:channels, 1:channelsb)
      else
         write(funit) (fermi_e(is) + pub_dmft_chem_shift), energy, &
              ttt(1:channels, 1:channelsb)
      endif

      if (ien == 1) then
         if(.not.pub_dmft_split)then
            write(funit2, '(a, i6, a, f15.8, a, f15.8)') 'Greens &
                 &function of spin ', is, ' and energy ', REAL(energy), &
                 ' + i ', AIMAG(energy)
         else
            write(funit2, '(a, i6, a, f15.8, a, f15.8)') 'Greens &
                 &function of spin ', is, ' and energy ', &
                 REAL(ener_buffer(kkkk)), ' + i ', AIMAG(ener_buffer(kkkk))
         endif
         write(funit2, '(a)') '----------------------------------------'
         write(funit2, '(a)') 'REAL PART'
         write(funit2, '(a,2i3)') 'channels, channelsb = ', channels, &
              channelsb
         do jj = 1, channels
            write(funit2, '(300f15.8)')(REAL(ttt(jj, kk)), kk = 1, &
                 channelsb)
         enddo
         write(funit2, '(a)') 'IMAGINARY PART'
         do jj = 1, channels
            write(funit2, '(300f15.8)')(AIMAG(ttt(jj, kk)), kk = 1, &
                 channelsb)
         enddo

         write(funit2, '(a)') '----------occupation matrix-------------'
      end if

      ttt = dmft_merged_occupancy_block(hub_atom, 1, l1, h_atoms, hub, &
           hub_proj_basis, merge_neigh, merge_table, nmerge, rot_vec_angles, &
           rot_vec_angles_split, cluster)
      if (ien == 1) then
         do jj = 1, channels
            write(funit2, '(300f15.8)')(real(ttt(jj, kk)), kk = 1, &
                 channelsb)
         enddo
         write(funit2, '(a)')'########################################'
         close(funit2)
      end if

      ! ! ebl: Quality control - print first row of the projected Green's
      ! ! function
      ! if (ien == 1 .and. is == 1) then
      !    do jj = 1, channels
      !       call utils_qc_print("Real(G_{"//trim(adjustl(utils_int_to_str(jj)))&
      !            //",1})", real(ttt(jj,1)))
      !    end do
      ! end if

   enddo
   !---------------------!

   if(pub_on_root .or. .not.pub_dmft_split) close(funit)
   if(pub_dmft_split) then
      deallocate(ttt_buffer, stat = ierr)
      call utils_dealloc_check('dmft_hubbard_greenf_info', 'ttt_buffer', ierr)
      deallocate(ener_buffer, stat = ierr)
      call utils_dealloc_check('dmft_hubbard_greenf_info', 'ener_buffer', ierr)
      deallocate(fre_buffer, stat = ierr)
      call utils_dealloc_check('dmft_hubbard_greenf_info', 'fre_buffer', ierr)
      deallocate(ener_bufferR, stat = ierr)
      call utils_dealloc_check('dmft_hubbard_greenf_info', 'ener_bufferR', ierr)
      deallocate(ener_bufferI, stat = ierr)
      call utils_dealloc_check('dmft_hubbard_greenf_info', 'ener_bufferI', ierr)
      deallocate(at_buffer, stat = ierr)
      call utils_dealloc_check('dmft_hubbard_greenf_info', 'at_buffer', ierr)
   end if

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_hubbard_greenf_info'

 end subroutine

 subroutine check_sparse_in_file(mat, filename)

   !---------------------------------------------------------------------------!
   ! Checks if a file exists. If it does, read the sparse matrix contained     !
   ! within it; if not, write the sparse matrix to it                          !
   !---------------------------------------------------------------------------!

   use sparse, only: sparse_read, sparse_write
   use comms,  only: comms_barrier, pub_on_root
   use rundat, only: pub_debug_on_root, pub_dmft_read, pub_dmft_write

   implicit none

   ! Arguments
   type(SPAM3),        intent(inout) :: mat
   character(len = *), intent(in   ) :: filename

   ! Local variables
   logical :: check

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &check_sparse_in_file'

   call comms_barrier
   inquire(file = filename, exist = check)
   call comms_barrier
   if(.not.check)then
      if (pub_on_root) write(stdout, '(a)', advance="no") &
           trim(adjustl(filename))//' is not present; '
      if (pub_dmft_write) then
         if (pub_on_root) write(stdout, '(a)') "writing"
         if (pub_dmft_mpi_rank == 1) call sparse_write(mat, filename)
      else
         if (pub_on_root) write(stdout, '(a)') "it will not be generated &
              &since dmft_write is set to FALSE"
      end if
   else
      if (pub_on_root) write(stdout, '(a)', advance="no") &
           trim(adjustl(filename))//' is present; '
      if (pub_dmft_read) then
         if (pub_on_root) write(stdout, '(a)') 'reading'
         if (pub_dmft_mpi_rank == 1) call sparse_read(mat, filename = filename)
      else
         if (pub_on_root) write(stdout, '(a)') "it will not be read &
              &since dmft_read is set to FALSE"
      end if
   endif
   call comms_barrier
   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &check_sparse_in_file'

 end subroutine

 subroutine dmft_sparse_realpart(mat)

   !---------------------------------------------------------------------------!
   ! Isolate the real part of a sparse matrix                                  !
   !---------------------------------------------------------------------------!

   use sparse, only: sparse_convert, sparse_num_cols, sparse_num_rows
   use comms,  only: comms_barrier
   use rundat, only: pub_debug_on_root
   use utils,  only: utils_alloc_check, utils_dealloc_check

   implicit none

   ! Arguments
   type(SPAM3), intent(inout) :: mat

   ! Local variables
   complex(kind = DP), allocatable, dimension(:,:) :: mat_square
   integer                                         :: n1, n2
   integer                                         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_sparse_realpart'

   call comms_barrier

   n1 = sparse_num_rows(mat)
   n2 = sparse_num_cols(mat)

   allocate(mat_square(n1, n2), stat = ierr)
   call utils_alloc_check('dmft_sparse_realpart', 'mat_square', ierr)

   call sparse_convert(mat_square, mat)
   mat_square = real(mat_square)
   call sparse_convert(mat, mat_square)

   deallocate(mat_square, stat = ierr)
   call utils_dealloc_check('dmft_sparse_realpart', 'mat_square', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_sparse_realpart'

 end subroutine

 subroutine dmft_print_spam3_info(mat)

   !---------------------------------------------------------------------------!
   ! Print information associated with a sparse matrix                         !
   !---------------------------------------------------------------------------!

   use sparse, only: sparse_convert, sparse_num_cols, sparse_num_rows
   use comms,  only: comms_barrier
   use rundat, only: pub_debug_on_root
   use utils,  only: utils_alloc_check, utils_dealloc_check

   implicit none

   ! Arguments
   type(SPAM3), intent(in) :: mat

   ! Local variables
   real(kind = DP), allocatable    :: mat_square_real(:,:)
   complex(kind = DP), allocatable :: mat_square_cplx(:,:)
   integer                         :: n1, n2, i
   integer                         :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_print_spam3_info'

   call comms_barrier
   n1 = sparse_num_rows(mat)
   n2 = sparse_num_cols(mat)

   if(n1 == 0 .and. n2 == 0)then
      write(stdout, *) 'WARNING: check sparse matrix, 0-shape matrix'
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_print_spam3_info'
      return
   endif

   if (mat%iscmplx) then
      allocate(mat_square_cplx(n1, n2), stat = ierr)
      call utils_alloc_check('dmft_print_spam3_info', 'mat_square_cplx', ierr)
      call sparse_convert(mat_square_cplx, mat)
   else
      allocate(mat_square_real(n1, n2), stat = ierr)
      call utils_alloc_check('dmft_print_spam3_info', 'mat_square_real', ierr)
      call sparse_convert(mat_square_real, mat)
   end if

   if (mat%iscmplx) then
      if(maxval(abs(mat_square_cplx)) < 1.0E-8_DP)then
         write(stdout, *) 'WARNING: empty matrix'
      endif
      write(stdout, *) 'WARNING: min-max values of abs(matrix) = ', &
           minval(abs(mat_square_cplx)), maxval(abs(mat_square_cplx))
      write(stdout, *) 'WARNING: min-max values of abs(imag(matrix)) = ', &
           minval(abs(aimag(mat_square_cplx))), &
           maxval(abs(aimag(mat_square_cplx)))
      write(stdout, *) 'WARNING: min-max values of abs(real(matrix)) = ', &
           minval(abs(real(mat_square_cplx))), &
           maxval(abs(real(mat_square_cplx)))
      write(stdout, *) 'WARNING: sum(mat_square) = ', sum(mat_square_cplx)
      write(stdout, *) 'WARNING: sum values (xline number) = ', sum( (/( &
           sum(mat_square_cplx(i, :)*real(i, kind=DP)), i = 1, size(mat_square_cplx, 1) &
           )/) )
      write(stdout, *) 'WARNING: sum values (xcol number) = ', sum( (/( &
           sum(mat_square_cplx(:, i)*real(i, kind=DP)), i = 1, size(mat_square_cplx, 1) &
           )/) )
   else
      if(maxval(abs(mat_square_real)) < 1.0E-8_DP)then
         write(stdout, *) 'WARNING: empty matrix'
      endif
      write(stdout, *) 'WARNING: min-max values of abs(matrix) = ', &
           minval(abs(mat_square_real)), maxval(abs(mat_square_real))
      write(stdout, *) 'WARNING: sum(mat_square) = ', sum(mat_square_real)
   end if

   if (mat%iscmplx) then
      deallocate(mat_square_cplx, stat = ierr)
      call utils_dealloc_check('dmft_print_spam3_info', 'mat_square_cplx', ierr)
   else
      deallocate(mat_square_real, stat = ierr)
      call utils_dealloc_check('dmft_print_spam3_info', 'mat_square_real', ierr)
   end if
   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_print_spam3_info'

 end subroutine

 subroutine dmft_normalise_projectors(mat, proj_on_line_or_col, overlap, &
      use_gpu_eigenvectors, mat_t)

   !---------------------------------------------------------------------------!
   ! Reorthogonalise the Hubbard projectors to obtain an identity-protected    !
   ! overlap matrix                                                            !
   !---------------------------------------------------------------------------!

   ! #ifdef GPU_SPEEDUP
   ! use fortran_cuda, only: eigenvector_cuda_r
   ! #endif
   use comms,  only: comms_barrier, comms_bcast, pub_on_root, &
                     pub_total_num_procs
   use linalg, only: linalg_dsyev_lt
   use rundat, only: pub_debug, pub_debug_on_root, pub_dmft_read, &
        pub_dmft_write, pub_rootname
   use sparse, only: sparse_convert, sparse_num_cols, sparse_num_rows
   use utils,  only: utils_abort, utils_assert, utils_unit, &
                     utils_alloc_check, utils_dealloc_check

   implicit none

   ! Arguments
   type(SPAM3),           intent(inout) :: mat
   integer,               intent(in   ) :: proj_on_line_or_col
   type(SPAM3),           intent(in   ) :: overlap
   logical,               intent(in   ) :: use_gpu_eigenvectors
   type(SPAM3), optional, intent(inout) :: mat_t

   ! Local variables
   real(kind = DP), allocatable, dimension(:,:) :: mat_square, Sinv
   real(kind = DP), allocatable, dimension(:,:) :: L_LT
   integer                                      :: nn1, nn2, ierr, n1, n2, i, &
                                                   j, k
   integer                                      :: mat_len
   real(kind = DP)                              :: normhub
   logical                                      :: use_cholesky, u_dagger_u, &
                                                   flip, check
   real(kind = DP), allocatable                 :: vaps(:), temp(:,:)
   integer                                      :: funit

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_normalise_projectors'

   use_cholesky = abs(pub_dmft_norm_proj) == 2 .or. abs(pub_dmft_norm_proj) == &
        3 .or. abs(pub_dmft_norm_proj) == 4
   u_dagger_u = abs(pub_dmft_norm_proj) == 3

   call comms_barrier
   n1 = sparse_num_rows(mat)
   n2 = sparse_num_cols(mat)
   if(n1 == 0 .and. n2 == 0)then
      call utils_abort("Error in dmft_normalize_projectors: projectors are &
           &void")
   endif
   allocate(mat_square(n1, n2), stat = ierr)
   call utils_alloc_check('dmft_normalise_projectors', 'mat_square', ierr)
   call sparse_convert(mat_square, mat)

   write(stdout, *) 'USE CHOLESKY? ', use_cholesky

   if(use_cholesky)then
      write(stdout, *) 'USING ORTHONORMALIZATION OF PROJECTORS'
      if((n1 < n2 .and. u_dagger_u) .or. (n1 > n2 .and. .not.u_dagger_u))then
         nn1 = n1
         nn2 = n2
         flip = .false.
      else
         nn1 = n2
         nn2 = n1
         flip = .true.
      endif

      if (allocated(L_LT)) then
         deallocate(L_LT, stat=ierr)
         call utils_dealloc_check('dmft_normalise_projectors', 'L_LT', ierr)
      end if
      allocate(L_LT(nn1, nn1), stat = ierr)
      call utils_alloc_check('dmft_normalise_projectors', 'L_LT', ierr)
      if(pub_dmft_use_gpu_onlyme)then
         ! #ifdef GPU_SPEEDUP
         ! if(.not.flip)then
         !    i = size(mat_square, 1)
         !    j = size(mat_square, 2)
         !    k = size(L_LT, 2)
         !    if(pub_on_root)then
         !       CALL matmulcuda_r(mat_square, transpose(mat_square), L_LT, i, &
         !            j, k)
         !    endif
         ! else
         !    i = size(mat_square, 2)
         !    j = size(mat_square, 1)
         !    k = size(L_LT, 2)
         !    if(pub_on_root .or. pub_dmft_split)then
         !       CALL matmulcuda_r(transpose(mat_square), mat_square, L_LT, i, &
         !            j, k)
         !    endif
         ! endif
         ! if(pub_total_num_procs > 1 .and. .not.pub_dmft_split) call &
         !      comms_bcast(0, L_LT)
         ! #else
         write(stdout, *) 'error, use GPU but no compiled with libraries'
         ! #endif
      else
         if(.not.flip)then
            L_LT = matmul(mat_square, transpose(mat_square))
         else
            L_LT = matmul(transpose(mat_square), mat_square)
         endif
      endif
      write(stdout, *) 'O symmetric?', maxval(abs((L_LT-transpose(L_LT))/2.0_DP))
      L_LT = (L_LT + transpose(L_LT))/2.0_DP
      if (allocated(vaps)) then
         deallocate(vaps, stat=ierr)
         call utils_dealloc_check('dmft_normalise_projectors', 'vaps', ierr)
         deallocate(temp, stat=ierr)
         call utils_dealloc_check('dmft_normalise_projectors', 'temp', ierr)
      end if
      allocate(vaps(nn1), stat = ierr)
      call utils_alloc_check('dmft_normalise_projectors', 'vaps', ierr)
      allocate(temp(nn1, nn1), stat = ierr)
      call utils_alloc_check('dmft_normalise_projectors', 'temp', ierr)
      temp = L_LT
      do i = 1, size(temp, 1)
         temp(i, i) = 0.0_DP
      enddo
      if(maxval(abs(temp)) > 1.0E-3_DP)then
         write(stdout, *) 'L_LT is not diagonal'
         write(stdout, *) 'max off-diagonal : ', maxval(abs(temp))
         write(stdout, *) 'max     diagonal : ', maxval(abs(diag(L_LT)))
         if(use_gpu_eigenvectors .and. pub_dmft_use_gpu_onlyme)then
            if(pub_on_root .or. pub_dmft_split)then
               ! #ifdef GPU_SPEEDUP
               ! call eigenvector_cuda_r(nn1, L_LT, vaps, temp, .false.)
               ! #else
               call utils_abort("Error in dmft_normalise_projectors: error &
                    &with gpu 3")
               ! #endif
               if (pub_debug) then
                  call linalg_dsyev_lt(L_LT, vaps, nn1)
                  if(maxval(abs(temp-L_LT)) > 1.0E-5_DP)then
                     write(stdout, *) 'ERROR IN EIGENVECTORS LLT!'
                     call utils_abort('Error in DMFT: temp and L_LT are too &
                          &dissimilar')
                  endif
               end if
               L_LT = temp
            endif
            if(pub_total_num_procs > 1 .and. .not.pub_dmft_split) call &
                 comms_bcast(0, L_LT)
            if(pub_total_num_procs > 1 .and. .not.pub_dmft_split) call &
                 comms_bcast(0, vaps)
         else
            call linalg_dsyev_lt(L_LT, vaps, nn1)
         endif
      else
         write(stdout, *) 'L_LT is diagonal'
         vaps = diag(L_LT)
      endif
      vaps = 1.0_DP/sqrt(abs(vaps))
      temp = 0.0_DP
      do i = 1, nn1
         temp(i, i) = vaps(i)
      enddo

      if(pub_dmft_use_gpu_onlyme)then
         ! #ifdef GPU_SPEEDUP
         ! if(pub_on_root .or. pub_dmft_split)then
         !    i = size(temp, 1)
         !    j = size(temp, 2)
         !    k = size(temp, 2)
         !    CALL matmulcuda_r(temp, transpose(L_LT), temp, i, j, k)
         !    i = size(L_LT, 1)
         !    j = size(L_LT, 2)
         !    k = size(L_LT, 2)
         !    CALL matmulcuda_r(L_LT, temp, L_LT, i, j, k)
         !    if(.not.flip)then
         !       i = size(L_LT, 1)
         !       j = size(L_LT, 2)
         !       k = size(mat_square, 2)
         !       CALL matmulcuda_r(L_LT, mat_square, mat_square, i, j, k)
         !    else
         !       i = size(mat_square, 1)
         !       j = size(mat_square, 2)
         !       k = size(mat_square, 2)
         !       CALL matmulcuda_r(mat_square, L_LT, mat_square, i, j, k)
         !    endif
         ! endif
         ! if(pub_total_num_procs > 1 .and. .not.pub_dmft_split) call &
         !      comms_bcast(0, mat_square)
         ! #else
         call utils_abort("Error in dmft_normalise_projectors: attempted to &
              &run a GPU calculation using a binary compiled without GPU &
              &capabilities")
         ! #endif
      else
         L_LT = matmul( matmul(L_LT,  temp), transpose(L_LT) )
         if(.not.flip)then
            mat_square = matmul(L_LT, mat_square)
         else
            mat_square = matmul(mat_square, L_LT)
         endif
      endif

      if(abs(pub_dmft_norm_proj) == 4)then

         if(allocated(Sinv)) then
            deallocate(Sinv, stat = ierr)
            call utils_dealloc_check('dmft_normalise_projectors', 'Sinv', ierr)
         end if

         mat_len = sparse_num_rows(overlap)
         allocate(Sinv(mat_len, mat_len), stat = ierr)
         call utils_alloc_check('dmft_normalise_projectors', 'Sinv', ierr)
         call sparse_convert(Sinv, overlap)

         if(allocated(vaps)) then
            deallocate(vaps, stat = ierr)
            call utils_dealloc_check('dmft_normalize_projectors', 'vaps', ierr)
            deallocate(temp, stat = ierr)
            call utils_dealloc_check('dmft_normalize_projectors', 'temp', ierr)
         end if

         allocate(vaps(mat_len), stat = ierr)
         call utils_alloc_check('dmft_normalise_projectors', 'vaps', ierr)
         allocate(temp(mat_len, mat_len), stat = ierr)
         call utils_alloc_check('dmft_normalise_projectors', 'temp', ierr)

         inquire(file = trim(pub_rootname) // '.sinv', exist = check)
         funit = utils_unit()
         if (check .and. pub_dmft_read) then
            open(unit = funit, file = trim(pub_rootname) // '.sinv', form = &
                 'unformatted')
            read(funit) Sinv, vaps
            close(funit)
         else
            call linalg_dsyev_lt(Sinv, vaps, mat_len)
            if(pub_dmft_mpi_rank == 1 .and. pub_on_root .and. pub_dmft_write) then
               open(unit = funit, file = trim(pub_rootname) // '.sinv', form = &
                    'unformatted')
               write(funit) Sinv, vaps
               close(funit)
            endif
         endif
         vaps = sqrt(abs(vaps))
         temp = 0.0_DP
         do i = 1, mat_len
            temp(i, i) = vaps(i)
         enddo

         if(pub_dmft_use_gpu_onlyme)then
            ! #ifdef GPU_SPEEDUP
            ! if(pub_on_root .or. pub_dmft_split)then
            !    i = size(temp, 1)
            !    j = size(temp, 2)
            !    k = size(temp, 2)
            !    CALL matmulcuda_r(temp, transpose(Sinv), temp, i, j, k)
            !    i = size(Sinv, 1)
            !    j = size(Sinv, 2)
            !    k = size(temp, 2)
            !    CALL matmulcuda_r(Sinv, temp, temp, i, j, k)
            !    if(size(mat_square, 1) == mat_len)then
            !       i = size(temp, 1)
            !       j = size(temp, 2)
            !       k = size(mat_square, 2)
            !       CALL matmulcuda_r(temp, mat_square, mat_square, i, j, k)
            !    elseif(size(mat_square, 2) == mat_len)then
            !       i = size(mat_square, 1)
            !       j = size(mat_square, 2)
            !       k = size(mat_square, 2)
            !       CALL matmulcuda_r(mat_square, temp, mat_square, i, j, k)
            !    else
            !       write(stdout, *) 'shape mat_square : ', shape(mat_square)
            !       write(stdout, *) 'shape Sinv       : ', shape(Sinv)
            !       call utils_abort("Error in dmft_normalise_projectors: &
            !            &mismatching shapes")
            !    endif
            ! endif
            ! if(pub_total_num_procs > 1 .and. .not.pub_dmft_split) call &
            !      comms_bcast(0, mat_square)
            ! #else
            call utils_abort("Error in dmft_normalise_projectors: attempted to &
                 &run a GPU calculation using a binary compiled without GPU &
                 &capabilities")
            ! #endif
         else
            temp = matmul(matmul(Sinv, temp), transpose(Sinv))
            if(size(mat_square, 1) == mat_len)then
               mat_square = matmul(temp, mat_square)
            elseif(size(mat_square, 2) == mat_len)then
               mat_square = matmul(mat_square, temp)
            else
               write(stdout, *) 'shape mat_square : ', shape(mat_square)
               write(stdout, *) 'shape Sinv       : ', shape(Sinv)
               call utils_abort("Error in dmft_normalise_projectors: &
                    &mismatching shapes")
            endif
         endif
      endif
   endif

   if(.not.use_cholesky)then
      if(n1 > n2)then
         call utils_assert(proj_on_line_or_col == 2, "Error in &
              &dmft_normalise_projectors: more projectors than NGWFs")
         write(stdout, *) 'SCANNING HUBOVERLAP MATRIX, COLUMNS'
         write(stdout, *) 'WARNING: please check the normalization of the &
              &projectors...'
         do i = 1, n2
            normhub = sum(abs(mat_square(:, i))**2)
            mat_square(:, i) = mat_square(:, i) / sqrt(normhub)
         enddo
      else
         call utils_assert(proj_on_line_or_col == 1, "Error in &
              &dmft_normalise_projectors: more projectors than NGWFs")
         write(stdout, *) 'SCANNING HUBOVERLAP MATRIX, LINES'
         do i = 1, n1
            normhub = sum(abs(mat_square(i, :))**2)
            mat_square(i, :) = mat_square(i, :) / sqrt(normhub)
         enddo
      endif
   endif

   if(present(mat_t))then
      call sparse_convert(mat_t, transpose(mat_square))
   endif

   call sparse_convert(mat, mat_square)
   deallocate(mat_square, stat = ierr)
   call utils_dealloc_check('dmft_normalise_projectors', 'mat_square', ierr)

   write(stdout, *) 'projectors renormalised, return'

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_normalise_projectors'

 end subroutine

 real(kind = DP) function internal_volume_cell_angstrom(cell)

   !---------------------------------------------------------------------------!
   ! Return the volume of the cell in cubic angstroms                          !
   !---------------------------------------------------------------------------!

   use constants,       only: ANGSTROM
   use geometry,        only: operator(.CROSS.), operator(.DOT.)
   use simulation_cell, only: CELL_INFO

   implicit none

   ! Arguments
   type(CELL_INFO), intent(in) :: cell

   ! Local variables
   real(kind = DP) :: volume_cell

   volume_cell = abs((cell%a1 .CROSS. cell%a2) .DOT. cell%a3)
   internal_volume_cell_angstrom = volume_cell/(ANGSTROM**3)
 end function

 real(kind = DP) function internal_distance(hat, hatb, hub, cell, elements)

   !---------------------------------------------------------------------------!
   ! Description:                                                              !
   !---------------------------------------------------------------------------!

   use constants,       only: ANGSTROM
   use hubbard_build,   only: HUBBARD_MODEL
   use ion,             only: ELEMENT
   use simulation_cell, only: CELL_INFO

   implicit none

   ! Arguments
   integer,             intent(in   ) :: hat, hatb
   type(HUBBARD_MODEL), intent(in   ) :: hub
   type(CELL_INFO),     intent(in   ) :: cell
   type(ELEMENT),       intent(in   ) :: elements(:)

   ! Local variables
   integer         :: i, j, k
   real(kind = DP) :: real_lattice(3, 3), dist(-1:1, -1:1, -1:1), v(3)

   real_lattice(1, 1) = cell%a1%x/ANGSTROM
   real_lattice(1, 2) = cell%a1%y/ANGSTROM
   real_lattice(1, 3) = cell%a1%z/ANGSTROM
   real_lattice(2, 1) = cell%a2%x/ANGSTROM
   real_lattice(2, 2) = cell%a2%y/ANGSTROM
   real_lattice(2, 3) = cell%a2%z/ANGSTROM
   real_lattice(3, 1) = cell%a3%x/ANGSTROM
   real_lattice(3, 2) = cell%a3%y/ANGSTROM
   real_lattice(3, 3) = cell%a3%z/ANGSTROM
   do i = -1, 1
      do j = -1, 1
         do k = -1, 1
            v =  real(i, kind=DP)*real_lattice(1, :)
            v = v + real(j, kind=DP)*real_lattice(2, :)
            v = v + real(k, kind=DP)*real_lattice(3, :)
            dist(i, j, k) = &
                 internal_norm_vector(internal_coordinates_atom(hatb, hub, &
                 elements) + v-internal_coordinates_atom(hat, hub, elements))
         enddo
      enddo
   enddo
   internal_distance = minval(dist(:, :, :))
 end function

 real(kind = DP) function internal_myfrequ(energy)

   !---------------------------------------------------------------------------!
   ! Calculates the frequency omega from the (possibly complex) energy,        !
   ! accounting for sign conventions                                           !
   !---------------------------------------------------------------------------!

   use rundat, only: pub_dmft_temp

   implicit none

   ! Arguments
   complex(kind = DP), intent(in   ) :: energy

   if (pub_dmft_temp .gt. 0.0_DP) then
      internal_myfrequ = abs(aimag(energy))
   else
      internal_myfrequ = abs(real(energy))
   endif
 end function

 pure integer function dmft_totn_atom_merge(hub_atom, hub, hub_proj_basis, &
      merge_neigh, merge_table, nmerge, cluster)

   !---------------------------------------------------------------------------!
   ! Merge atoms to define a large AIM which mixes different orbital character !
   !---------------------------------------------------------------------------!

   use comms,             only: pub_my_proc_id
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par

   implicit none

   ! Arguments
   integer,              intent(in   ) :: hub_atom
   type(HUBBARD_MODEL),  intent(in   ) :: hub
   type(FUNC_BASIS),     intent(in   ) :: hub_proj_basis
   integer, allocatable, intent(in   ) :: merge_neigh(:)
   integer, allocatable, intent(in   ) :: merge_table(:,:)
   integer,              intent(in   ) :: nmerge
   logical,              intent(in   ) :: cluster

   ! Local variables
   integer              :: hat, totn, j
   integer              :: vv(0:nmerge)

   if(cluster)then
      vv = dmft_totn_atom_merge_detail(hub_atom, hub, hub_proj_basis, &
           merge_neigh, merge_table, nmerge)
      dmft_totn_atom_merge = sum(vv(1:))
   else
      hat = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
      dmft_totn_atom_merge = dmft_hub_atom_hub_projs(hat, hub, hub_proj_basis)
   endif

 end function

 pure function dmft_totn_atom_merge_detail(hub_atom, hub, hub_proj_basis, &
      merge_neigh, merge_table, nmerge)

   !---------------------------------------------------------------------------!
   ! Merge atoms to define a large AIM which mixes different orbital character !
   !---------------------------------------------------------------------------!

   use comms,             only: pub_my_proc_id
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par

   implicit none

   ! Arguments
   integer,              intent(in   ) :: hub_atom
   type(HUBBARD_MODEL),  intent(in   ) :: hub
   type(FUNC_BASIS),     intent(in   ) :: hub_proj_basis
   integer, allocatable, intent(in   ) :: merge_neigh(:)
   integer, allocatable, intent(in   ) :: merge_table(:,:)
   integer,              intent(in   ) :: nmerge

   ! Local variables
   integer              :: kk, hat0, j
   integer              :: dmft_totn_atom_merge_detail(0:nmerge)
   integer              :: ierr

   dmft_totn_atom_merge_detail = 0
   dmft_totn_atom_merge_detail(0:0) = 0

   if (pub_dmft_split) then
      hat0 = hub_atom
   else
      hat0 = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
   end if

   do kk = 1, merge_neigh(hat0)
      j = merge_table(hat0, kk)
      if (j > 0) then
         dmft_totn_atom_merge_detail(kk) = dmft_hub_atom_hub_projs(j, hub, &
              hub_proj_basis)
      endif
   enddo

 end function

 function dmft_merged_green_block(hub_atom, hub_atomb, n, m, &
      site_greenf_buffer, hub, hub_proj_basis, cluster, merge_neigh, &
      merge_table, nmerge)

   !---------------------------------------------------------------------------!
   ! Calculates the merged projected Green's function, which spans several     !
   ! atoms                                                                     !
   !---------------------------------------------------------------------------!
   ! NOTES                                                                     !
   !   ebl: consider passing loop and not n and m as an argument               !
   !---------------------------------------------------------------------------!

   use comms,             only: pub_my_proc_id
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par
   use rundat,            only: pub_debug, pub_debug_on_root
   use utils,             only: utils_abort

   implicit none

   ! Arguments
   integer,                         intent(in   ) :: hub_atom
   integer,                         intent(in   ) :: hub_atomb
   integer,                         intent(in   ) :: n
   integer,                         intent(in   ) :: m
   complex(kind = DP), allocatable, intent(in   ) :: site_greenf_buffer(:,:,:,:)
   type(HUBBARD_MODEL),             intent(in   ) :: hub
   type(FUNC_BASIS),                intent(in   ) :: hub_proj_basis
   logical,                         intent(in   ) :: cluster
   integer, allocatable,            intent(in   ) :: merge_neigh(:)
   integer, allocatable,            intent(in   ) :: merge_table(:,:)
   integer,                         intent(in   ) :: nmerge

   ! Local variables
   integer            :: s1, s2, kk, kkb, hat, hatb, totn, totnb, j, jb, &
                         v(0:nmerge), vb(0:nmerge)
   complex(kind = DP) :: dmft_merged_green_block(n, m)

   if(.not.cluster)then
      j = dmft_totn_atom_merge(hub_atom, hub, hub_proj_basis, merge_neigh, &
           merge_table, nmerge, cluster)
      jb = dmft_totn_atom_merge(hub_atomb, hub, hub_proj_basis, merge_neigh, &
           merge_table, nmerge, cluster)
      dmft_merged_green_block(1:j, 1:jb) = site_greenf_buffer(hub_atom, &
           hub_atomb, 1:j, 1:jb)
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_merged_green_block'
      return
   endif

   if(dmft_totn_atom_merge(hub_atom, hub, hub_proj_basis, merge_neigh, &
        merge_table, nmerge, cluster) /= dmft_totn_atom_merge(hub_atomb, hub, &
        hub_proj_basis, merge_neigh, merge_table, nmerge, cluster)) then
      call utils_abort('Error in dmft_merged_green_block: extended merging &
           &scheme not yet implememented for rectangular subblocks')
   endif
   dmft_merged_green_block = 0.0_DP
   v = dmft_totn_atom_merge_detail(hub_atom, hub, hub_proj_basis, merge_neigh, &
        merge_table, nmerge)
   vb = dmft_totn_atom_merge_detail(hub_atomb, hub, hub_proj_basis, &
        merge_neigh, merge_table, nmerge)

   if (pub_dmft_split) then
      hat = hub_atom
      hatb = hub_atomb
   else
      hat = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
      hatb = par%hub_atoms_on_proc(hub_atomb, pub_my_proc_id)
   end if
   do kk = 1, merge_neigh(hat)
      do kkb = 1, merge_neigh(hatb)
         j  = merge_table(hat, kk)
         jb = merge_table(hatb, kkb)
         s1 = sum( v (0:kk-1 ) )
         s2 = sum( vb(0:kkb-1) )
         ! if(pub_debug_on_root)then
         !    write(stdout, *) 'kk (1-nneigh) : ', kk, merge_neigh(hat)
         !    write(stdout, *) 'kkb           : ', kkb, merge_neigh(hatb)
         !    write(stdout, *) 's1-s2         : ', s1, s2
         !    write(stdout, *) 'vkk, vkk2      : ', v(kk), vb(kkb)
         !    write(stdout, *) 'j, jb          : ', j, jb
         ! endif
         dmft_merged_green_block(s1 + 1:s1 + v(kk), s2 + 1:s2 + vb(kkb)) = &
              site_greenf_buffer(dmft_hubbard_on_my_proc(j), &
              dmft_hubbard_on_my_proc(jb), 1:v(kk), 1:vb(kkb))
      enddo
   enddo

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_merged_green_block'
   return

 end function

 subroutine dmft_distribute(big, hub_atom, hub_atomb, n_, m_, loop, &
      site_self_energy_buffer, hub, hub_proj_basis, nmerge, merge_table, &
      merge_neigh, cluster)

   !---------------------------------------------------------------------------!
   ! Distribute the unfolded self energy to local atomic components            !
   !---------------------------------------------------------------------------!

   use comms,             only: pub_my_proc_id
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par
   use utils,             only: utils_assert
   use rundat,            only: pub_debug_on_root

   implicit none

   ! Arguments
   integer,              intent(in   ) :: n_
   integer,              intent(in   ) :: m_
   complex(kind = DP),   intent(in   ) :: big(n_, m_)
   integer,              intent(in   ) :: hub_atom
   integer,              intent(in   ) :: hub_atomb
   type(DMFT_LOOP_INFO), intent(inout) :: loop
   complex(kind = DP),   intent(inout) :: site_self_energy_buffer(:,:,:,:)
   type(HUBBARD_MODEL),  intent(in   ) :: hub
   type(FUNC_BASIS),     intent(in   ) :: hub_proj_basis
   integer,              intent(in   ) :: nmerge
   integer, allocatable, intent(in   ) :: merge_table(:,:)
   integer, allocatable, intent(in   ) :: merge_neigh(:)
   logical,              intent(in   ) :: cluster

   ! Local variables
   integer :: s1, s2, kk, kkb, totn, totnb, j, jb, v(0:nmerge), vb(0:nmerge)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering dmft_distribute'

   j = dmft_totn_atom_merge(hub_atom, hub, hub_proj_basis, merge_neigh, &
        merge_table, nmerge, cluster)
   jb = dmft_totn_atom_merge(hub_atomb, hub, hub_proj_basis, merge_neigh, &
        merge_table, nmerge, cluster)

   if(.not.cluster)then
      site_self_energy_buffer(hub_atom, hub_atomb, 1:j, 1:jb) = big(1:j, 1:jb)
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_distribute'
      return
   endif

   call utils_assert(j == jb, 'ERROR in dmft_distribute: extended merging &
        &scheme not implemented for rectangular subblocks')

   v = dmft_totn_atom_merge_detail(hub_atom, hub, hub_proj_basis, merge_neigh, &
        merge_table, nmerge)
   vb = dmft_totn_atom_merge_detail(hub_atomb, hub, hub_proj_basis, &
        merge_neigh, merge_table, nmerge)
   loop%hat = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
   loop%hatb = par%hub_atoms_on_proc(hub_atomb, pub_my_proc_id)

   do kk = 1, merge_neigh(loop%hat)
      do kkb = 1, merge_neigh(loop%hatb)
         j  = merge_table(loop%hat, kk)
         jb = merge_table(loop%hatb, kkb)
         s1 = sum( v (0:kk-1 ) )
         s2 = sum( vb(0:kkb-1) )
         site_self_energy_buffer(dmft_hubbard_on_my_proc(j), &
              dmft_hubbard_on_my_proc(jb), 1:v(kk), 1:vb(kkb)) = big(s1 + 1:s1 &
              + v(kk), s2 + 1:s2 + vb(kkb))
      enddo
   enddo

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving dmft_distribute'

 end subroutine


 function dmft_merged_occupancy_block(hub_atom, spin, n_, h_atoms, hub, &
      hub_proj_basis, merge_neigh, merge_table, nmerge, rot_vec_angles, &
      rot_vec_angles_split, cluster)

   !---------------------------------------------------------------------------!
   ! Calculates the occupany matrix merged from the components obtained from   !
   ! several atoms                                                             !
   !---------------------------------------------------------------------------!

   use comms,             only: pub_my_proc_id
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par
   use rundat,            only: pub_debug_on_root, pub_dmft_rotate_green
   use utils,             only: utils_alloc_check, utils_dealloc_check
   implicit none

   ! Arguments
   integer,                      intent(in   ) :: hub_atom
   integer,                      intent(in   ) :: spin
   integer,                      intent(in   ) :: n_
   type(ARRAY_OF_MATRICES),      intent(in   ) :: h_atoms(par%nat_hub)
   type(HUBBARD_MODEL),          intent(in   ) :: hub
   type(FUNC_BASIS),             intent(in   ) :: hub_proj_basis
   integer, allocatable,         intent(in   ) :: merge_neigh(:)
   integer, allocatable,         intent(in   ) :: merge_table(:,:)
   integer,                      intent(in   ) :: nmerge
   real(kind = DP), allocatable, intent(in   ) :: rot_vec_angles(:,:,:)
   real(kind = DP), allocatable, intent(in   ) :: rot_vec_angles_split(:,:,:)
   logical,                      intent(in   ) :: cluster

   ! Local variables
   integer                      :: s1, kk, totn, j, jproc
   integer                      :: hat
   integer                      :: v(0:nmerge)
   real(kind = DP), allocatable :: dmft_merged_occupancy_block(:,:)
   integer                      :: ierr

   allocate(dmft_merged_occupancy_block(n_, n_), stat = ierr)
   call utils_alloc_check('dmft_merged_occupancy_block', &
        'dmft_merged_occupancy_block', ierr)

   if (pub_dmft_split) then
      hat = hub_atom
   else
      hat = par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
   end if

   if (.not.cluster) then
      j = dmft_totn_atom_merge(hub_atom, hub, hub_proj_basis, merge_neigh, &
           merge_table, nmerge, cluster)
      dmft_merged_occupancy_block(1:j, 1:j) = &
           diago_occupancy_func(h_atoms(hat)%occupancy(1:j, 1:j, spin), j, &
           pub_dmft_rotate_green, hat, rot_vec_angles(hub_atom, 1:3, 1:3))
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_merged_occupancy_block'
      return
   endif

   dmft_merged_occupancy_block = 0.0_DP
   v = dmft_totn_atom_merge_detail(hub_atom, hub, hub_proj_basis, merge_neigh, &
        merge_table, nmerge)

   do kk = 1, merge_neigh(hat)
      j = merge_table(hat, kk)
      s1 = sum( v(0:kk-1) )

      ! Neaten this when rot_vec_angles and rot_vec_angles_split have been
      ! turned into
      ! one object
      if (pub_dmft_split) then
         dmft_merged_occupancy_block(s1 + 1:s1 + v(kk), s1 + 1:s1 + v(kk)) = &
              diago_occupancy_func(h_atoms(j)%occupancy(1:v(kk), 1:v(kk), &
              spin), v(kk), pub_dmft_rotate_green, hat, &
              rot_vec_angles_split(j, 1:3, 1:3))
      else
         dmft_merged_occupancy_block(s1 + 1:s1 + v(kk), s1 + 1:s1 + v(kk)) = &
              diago_occupancy_func(h_atoms(j)%occupancy(1:v(kk), 1:v(kk), &
              spin), v(kk), pub_dmft_rotate_green, hat, &
              rot_vec_angles(dmft_hubbard_on_my_proc(j), 1:3, 1:3))
      end if
   enddo

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_merged_occupancy_block'
   return
 end function

 subroutine dmft_set_loop_variable(loop, hub_atom, hub_atomb, hub, &
      hub_proj_basis)

   !---------------------------------------------------------------------------!
   ! Updates the "loop" variable to correspond to the hubbard atoms "hub_atom" !
   ! and "hub_atomb"                                                           !
   !---------------------------------------------------------------------------!

   use comms,             only: pub_my_proc_id
   use rundat,            only: pub_debug, pub_debug_on_root
   use function_basis,    only: FUNC_BASIS
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par

   implicit none

   ! Arguments
   type(DMFT_LOOP_INFO), intent(  out) :: loop
   integer,              intent(in   ) :: hub_atom, hub_atomb
   type(HUBBARD_MODEL),  intent(in   ) :: hub
   type(FUNC_BASIS),     intent(in   ) :: hub_proj_basis

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_set_loop_variable'

   if (pub_dmft_split) then
      loop%hat    =  hub_atom
   else
      loop%hat    =  par%hub_atoms_on_proc(hub_atom, pub_my_proc_id)
   end if
   loop%theatom   =  par%distr_atom( hub%orig(loop%hat)  )
   loop%channels = hub_proj_basis%num_on_atom(loop%theatom)
   loop%sp        =  hub%species_number(loop%hat)
   if (pub_dmft_split) then
      loop%hatb   =  hub_atomb
   else
      loop%hatb   =  par%hub_atoms_on_proc(hub_atomb, pub_my_proc_id)
   end if
   loop%theatomb  =  par%distr_atom( hub%orig(loop%hatb) )
   loop%channelsb = hub_proj_basis%num_on_atom(loop%theatomb)
   loop%spb       =  hub%species_number(loop%hatb)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_set_loop_variable'

 end subroutine

 subroutine dmft_print_loop_info(loop, hub)

   !---------------------------------------------------------------------------!
   ! Prints information associated with the primary hubbard atom in "loop"     !
   !---------------------------------------------------------------------------!

   use comms,             only: pub_my_proc_id, pub_on_root
   use hubbard_build,     only: HUBBARD_MODEL
   use parallel_strategy, only: par => pub_par
   use rundat,            only: pub_debug, pub_debug_on_root

   implicit none

   ! Arguments
   type(HUBBARD_MODEL),  intent(in   ) :: hub
   type(DMFT_LOOP_INFO), intent(in   ) :: loop

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_print_loop_info'

   if(pub_output_detail /= VERBOSE) return

   if (pub_debug_on_root) then
      write(stdout, '(a, i2, a, i2)') 'Proc ', pub_my_proc_id, ' is treating atom &
           &', loop%hat
      write(stdout, '(a,i3)') ' Hubbard atom  : ', par%distr_atom( &
           hub%orig(loop%hat) )
      write(stdout, '(a,i3)') ' hub%orig(hat) : ', hub%orig(loop%hat)
      write(stdout, '(a,i3)') ' the atom      : ', loop%theatom
      write(stdout, '(a,i3)') ' species       : ', loop%sp
      write(stdout, '(a,i3)') ' first atom    : ', par%first_atom_on_proc(pub_my_proc_id)
      write(stdout, '(a,i3)') ' num atoms     : ', par%num_atoms_on_proc(pub_my_proc_id)
      write(stdout, '(a,i3)') ' orig atom     : ', par%orig_atom(loop%hat)
   end if

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_print_loop_info'

 end subroutine


 logical function internal_doubleprec(energy, chem, fermi_e, is, ien, &
      nspecial_frequ)

   !---------------------------------------------------------------------------!
   ! Decides whether or not to perform GPU calculations with single or double  !
   ! precision. Double precision will be used at frequencies between 0 and     !
   ! pub_dmft_cutoff_small, or if frequencies above pub_dmft_cutoff_tail are   !
   ! considered.                                                               !
   !---------------------------------------------------------------------------!

   use rundat, only: pub_dmft_complex_freq, &
        pub_dmft_cutoff_small, pub_dmft_points, &
        pub_dmft_temp

   implicit none

   ! Arguments
   complex(kind = DP), intent(in   ) :: energy
   real(kind = DP),    intent(  out) :: chem
   real(kind = DP),    intent(in   ) :: fermi_e(2)
   integer,            intent(in   ) :: is
   integer,            intent(in   ) :: ien
   integer,            intent(in   ) :: nspecial_frequ

   chem = fermi_e(is) + pub_dmft_chem_shift

   internal_doubleprec = pub_dmft_temp > 0.0_DP .and. (internal_myfrequ(energy) < &
        pub_dmft_cutoff_small .or. internal_myfrequ(energy) > &
        pub_dmft_cutoff_tail)

   internal_doubleprec = internal_doubleprec .or. ien >= &
        pub_dmft_points-nspecial_frequ + 1 .or. (ien == 1 .and. .not. &
        pub_dmft_complex_freq)

   internal_doubleprec = internal_doubleprec .or. ((.not. &
        pub_dmft_complex_freq) .and. ( abs(real(energy)-chem) < &
        pub_dmft_cutoff_small .or. abs(real(energy)-chem) > &
        pub_dmft_cutoff_tail ) )
 end function

 subroutine dmft_sym_invert(i, mat)

   !---------------------------------------------------------------------------!
   ! Inverts a symmetric (complex) matrix                                      !
   !---------------------------------------------------------------------------!

   use linalg, only: linalg_sym_invert_serial
   use utils,  only: utils_abort
   use rundat, only: pub_debug_on_root

   implicit none

   ! Arguments
   integer,            intent(in   ) :: i
   complex(kind = DP), intent(inout) :: mat(:,:)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering dmft_sym_invert'

   if(i /= size(mat, 1) .or. i /= size(mat, 2))then
      call utils_abort("Error in dmft_sym_invert: matrix is not square")
   endif

   ! ebl: disabling GPU capabilites for now...
   ! #ifdef GPU_SPEEDUP
   ! if(pub_dmft_use_gpu_onlyme)then
   ! #ifdef GPU_SPEEDUP_SINGPREC
   !    call magma_fortran_comp_singleprec_bypass(i, mat)
   ! #else
   ! #ifdef GPU_HANGS
   !    call matinv_ge_complex(i, mat)
   ! #else
   !    call magma_fortran_comp_(i, mat)
   ! #endif
   ! #endif
   ! else
   !    call linalg_sym_invert_serial(mat, i)
   ! endif
   ! #else
   call linalg_sym_invert_serial(mat, i)
   ! #endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving dmft_sym_invert'

 end subroutine

 subroutine dmft_invert_sparse_complex(mat, mat_inv, iter, backup, backupflag)

   !---------------------------------------------------------------------------!
   ! Directly invert a sparse (complex) matrix                                 !
   !---------------------------------------------------------------------------!

   use sparse, only: sparse_convert, sparse_num_cols, sparse_num_rows
   use linalg, only: linalg_invert_serial
   use utils,  only: utils_assert, utils_abort, utils_alloc_check, &
                     utils_dealloc_check
   use comms,  only: pub_total_num_procs
   use rundat, only: pub_debug_on_root

   implicit none

   ! Arguments
   type(SPAM3),                               intent(in   ) :: mat
   type(SPAM3),                               intent(inout) :: mat_inv
   integer,                                   intent(in   ) :: iter
   complex(kind = DP), optional, allocatable, intent(inout) :: backup(:,:)
   logical, optional,                         intent(inout) :: backupflag

   ! Local variables
   integer                         :: nn, nnn
   complex(kind = DP), allocatable :: mat_square(:,:)
   logical                         :: skip_invert
   integer                         :: ierr
   character(len = 3)              :: errStr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_invert_sparse_complex'

   call utils_assert(present(backup) .eqv. present(backupflag), "Error in &
        &dmft_invert_sparse_complex: one, but not both, of backup and backupflag &
        &provided")

   nn = sparse_num_rows(mat)

   ! #ifndef GPU_SPEEDUP
#ifndef INVERT_N_CUBE
   if (pub_total_num_procs > 1) then
      if (present(backupflag)) then
         call internal_invert_full(mat_inv, mat, nn, backup, backupflag)
      else
         call internal_invert_full(mat_inv, mat, nn)
      end if
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_invert_sparse_complex'
      return
   endif
#endif
   ! #endif

   ! #ifdef GPU_SPEEDUP
   ! #ifndef INVERT_N_CUBE
   ! if(.not.(pub_dmft_use_gpu .or. pub_dmft_use_gpu_partially)) then
   !    if(pub_total_num_procs > 1) then
   !       if (present(backupflag)) then
   !          call internal_invert_full(mat_inv, mat, nn, backup, backupflag)
   !       else
   !          call internal_invert_full(mat_inv, mat, nn)
   !       end if
   !       return
   !    endif
   ! endif
   ! #endif
   ! #endif

   allocate(mat_square(nn, nn), stat = ierr)
   call utils_alloc_check('dmft_invert_sparse_complex', 'mat_square', ierr)

   skip_invert = .false.
   if( pub_dmft_mpi_rank /= 1 .and. pub_dmft_use_gpu_partially .and. .not. &
        pub_dmft_splitk)then
      mat_square = 0.0_DP
      skip_invert = .true.
      ! goto 3131
   endif

   if (.not.skip_invert) then
      call sparse_convert(mat_square, mat)

      ! #ifdef GPU_SPEEDUP
      ! if(pub_dmft_use_gpu_onlyme)then
      !    ! #ifdef GPU_SPEEDUP_SINGPREC
      !    call magma_fortran_comp_singleprec_bypass(nn, mat_square)
      !    ! #else
      !    ! #ifdef GPU_HANGS
      !    call matinv_ge_complex(nn, mat_square)
      !    ! #else
      !    call magma_fortran_comp_(nn, mat_square)
      !    !#endif
      !    !#endif
      ! else
      !    call linalg_invert_serial(mat_square, ierr)
      !    write(errStr, '(i3)') ierr
      ! call utils_assert(ierr == 0, "Error in linalg_invert_serial: exit code
      ! " // adjustl(errStr))
      ! endif
      ! #else
      call linalg_invert_serial(mat_square, ierr)
      write(errStr, '(i3)') ierr
      call utils_assert(ierr == 0, "Error in linalg_invert_serial: exit code " &
           // adjustl(errStr))
      ! #endif
   end if

   if(pub_dmft_mpi_size > 1 .and. pub_dmft_use_gpu_partially .and. .not. &
        pub_dmft_splitk)then
      call utils_abort("Reached point in code where nfs sync used to appear; &
           &the code needs to be updated to use these settings.")
      ! call sync_via_nfs(mat_square, 'invert_matc' // &
      !      adjustl(trim(utils_int_to_str(iter))), nn)
   endif

   if(present(backup))then
      backupflag = .true.
      backup = mat_square
   endif

   call sparse_convert(mat_inv, mat_square)
   deallocate(mat_square, stat = ierr)
   call utils_dealloc_check('dmft_invert_sparse_complex', 'mat_square', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_invert_sparse_complex'

 end subroutine

 subroutine internal_sparse_max_diff(spam1, spam2)

   !---------------------------------------------------------------------------!
   ! Finds the largest difference between the elements of two sparse matrices  !
   !---------------------------------------------------------------------------!

   use sparse, only: sparse_axpy, sparse_copy, sparse_create, sparse_destroy, &
        sparse_max_abs_element
   use rundat, only: pub_debug_on_root

   implicit none

   ! Arguments
   type(SPAM3), intent(in) :: spam1, spam2

   ! Local variables
   type(SPAM3)     :: tmp
   real(kind = DP) :: diff

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &internal_sparse_max_diff'

   call sparse_create(tmp, spam1, spam2)
   call sparse_copy(tmp, spam1)
   call sparse_axpy(tmp, spam2, -1.0_DP)
   diff = sparse_max_abs_element(tmp)
   write(stdout, *) 'MAX DIFFERENCE BETWEEN SPARSE IS : ', diff
   call sparse_destroy(tmp)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &internal_sparse_max_diff'

 end subroutine

 subroutine dmft_invert_sparse_real(mat, mat_inv, iter)

   !---------------------------------------------------------------------------!
   ! Inverts a sparse (real) matrix                                            !
   !---------------------------------------------------------------------------!
   ! NOTES                                                                     !
   !   This is not supposed to be used to compute the Green function           !
   !---------------------------------------------------------------------------!

   use sparse, only: sparse_convert, sparse_num_cols, sparse_num_rows
   use linalg, only: linalg_invert_serial
   use utils,  only: utils_assert, utils_abort, utils_alloc_check, &
                     utils_dealloc_check
   use comms,  only: pub_total_num_procs
   use rundat, only: pub_debug_on_root

   implicit none

   ! Arguments
   integer,     intent(in   ) :: iter
   type(SPAM3), intent(inout) :: mat
   type(SPAM3), intent(inout) :: mat_inv

   ! Local variables
   integer                                      :: nn
   real(kind = DP), allocatable, dimension(:,:) :: mat_square
   logical                                      :: skip_invert
   integer                                      :: ierr
   character(len = 3)                           :: errStr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_invert_sparse_real'

   nn = sparse_num_rows(mat_inv)

   ! #ifndef GPU_SPEEDUP
#ifndef INVERT_N_CUBE
   if(pub_total_num_procs > 1) then
      call internal_invert_full_real(mat, mat_inv, nn)
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_invert_sparse_real'
      return
   endif
#endif
   ! #endif

   ! #ifdef GPU_SPEEDUP
   ! #ifndef INVERT_N_CUBE
   ! if(.not.(pub_dmft_use_gpu .or. pub_dmft_use_gpu_partially)) then
   !    if(pub_total_num_procs > 1) then
   !       call internal_invert_full_real(mat, mat_inv, nn)
   !       return
   !    endif
   ! endif
   ! #endif
   ! #endif

   allocate(mat_square(nn, nn), stat = ierr)
   call utils_alloc_check('dmft_invert_sparse_real', 'mat_square', ierr)

   skip_invert = .false.
   if( pub_dmft_mpi_rank /= 1 .and. pub_dmft_use_gpu_partially .and. .not. &
        pub_dmft_splitk)then
      mat_square = 0.0_DP
      skip_invert = .true.
   endif

   if (.not. skip_invert) then
      call sparse_convert(mat_square, mat)

      ! #ifdef GPU_SPEEDUP
      ! if(pub_dmft_use_gpu_onlyme)then
      !    ! #ifdef GPU_SPEEDUP_SINGPREC
      !    call magma_fortran_double_singleprec_bypass(nn, mat_square)
      !    ! #else
      !    call magma_fortran_double_(nn, mat_square)
      !    ! #endif
      ! else
      !    call linalg_invert_serial(mat_square, ierr)
      !    write(errStr, '(i3)') ierr
      ! call utils_assert(ierr == 0, "Error in linalg_invert_serial: exit code
      ! " // adjustl(errStr))
      ! endif
      ! #else
      call linalg_invert_serial(mat_square, ierr)
      write(errStr, '(i3)') ierr
      call utils_assert(ierr == 0, "Error in linalg_invert_serial: exit code " &
           // adjustl(errStr))
      ! #endif
   end if

   if(pub_dmft_mpi_size > 1 .and. pub_dmft_use_gpu_partially .and. .not. &
        pub_dmft_splitk)then
      call utils_abort("Reached point in code where nfs sync used to appear;&
           & the code needs to be updated to use these settings.")
      ! call sync_via_nfs(mat_square, 'invert_matr' // &
      !      adjustl(trim(utils_int_to_str(iter))), nn)
   endif

   call sparse_convert(mat_inv, mat_square)
   deallocate(mat_square, stat = ierr)
   call utils_dealloc_check('dmft_invert_sparse_real', 'mat_square', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_invert_sparse_real'

 end subroutine

 real(kind = DP) function dmft_mu_updated_by_bisection(muvec, Nvec, &
      target_N)

   !==========================================================================!
   ! # DESCRIPTION                                                            !
   ! Updates the chemical potential using the bisection method                !
   !--------------------------------------------------------------------------!
   ! # ARGUMENTS                                                              !
   ! muvec                        (in)  vector of guesses for the chemical    !
   !                                    potential mu                          !
   ! Nvec                         (in)  vector of resulting values for the    !
   !                                    number of electrons N                 !
   ! target_N                     (in)  the target value for N                !
   ! dmft_mu_updated_by_bisection (out) an updated value for mu               !
   !--------------------------------------------------------------------------!
   ! # AUTHORS & CHANGELOG                                                    !
   ! Author(s):        Edward Linscott                                        !
   ! Date of creation: May 2019                                               !
   !==========================================================================!

   use rundat, only: pub_debug_on_root, pub_dmft_mu_diff_max
   use utils, only: utils_abort, utils_assert

   implicit none

   ! Arguments
   real(kind = DP), intent(in   ) :: muvec(:)   ! vector containing mu
   real(kind = DP), intent(in   ) :: Nvec(:)    ! vector containng N
   real(kind = DP), intent(in   ) :: target_N   ! target N

   ! Local variables
   real(kind = DP)            :: N_upper_bound
   real(kind = DP)            :: N_lower_bound
   real(kind = DP)            :: mu_upper_bound
   real(kind = DP)            :: mu_lower_bound
   integer                    :: i
   logical                    :: found_upper
   logical                    :: found_lower
   real(kind = DP), parameter :: max_mu_step = 0.01_DP

   ! Work out bounds based on muvec and Nvec
   found_upper = .false.
   found_lower = .false.
   do i=1,size(muvec)
      if (abs(Nvec(i) - target_N) < pub_dmft_mu_diff_max) then
         ! If we've find the solution, exit immediately
         if (pub_debug_on_root) then
            write(stdout, '(a)') 'DEBUG: A previous mu gave the target N'
         end if
         dmft_mu_updated_by_bisection = muvec(i)
         return
      else if (Nvec(i) > target_N) then
         ! This mu is a candidate for the upper bound on the true mu
         if (.not. found_upper) then
            found_upper = .true.
         else if (Nvec(i) > N_upper_bound) then
            continue
         end if
         mu_upper_bound = muvec(i)
         N_upper_bound = Nvec(i)
      else
         ! This mu is a candidate for the lower bound on the true mu
         if (.not. found_lower) then
            found_lower = .true.
         else if (Nvec(i) < N_lower_bound) then
            continue
         end if
         mu_lower_bound = muvec(i)
         N_lower_bound = Nvec(i)
      end if
   end do

   if (found_lower .and. found_upper) then
      call utils_assert(N_upper_bound > N_lower_bound, 'Error in &
            &dmft_mu_updated_by_bisection: lower bound for N is greater than upper bound')
      call utils_assert(mu_upper_bound > mu_lower_bound, 'Error in &
            &dmft_mu_updated_by_bisection: lower bound for mu is greater than upper bound')
      dmft_mu_updated_by_bisection = 0.5*(mu_lower_bound + mu_upper_bound)
   else if (found_upper) then
      dmft_mu_updated_by_bisection = mu_upper_bound - max_mu_step
   else if (found_lower) then
      dmft_mu_updated_by_bisection = mu_lower_bound + max_mu_step
   else
      call utils_abort("Error in dmft_mu_updated_by_bisection: should not arrive here")
   end if

   if (pub_debug_on_root) then
      write(stdout, '(a)') 'DEBUG: Updating mu'
      write(stdout, '(a,f8.3)') 'Target N: ', target_N
      if (found_lower) then
         write(stdout, '(a,f11.6,a,f11.6,a)') 'Lower bound: ', mu_lower_bound, &
              ' (N = ', N_lower_bound, ')'
      else
         write(stdout, '(a,f11.6,a,f11.6)') 'Lower bound for mu not yet found'
      end if
      if (found_upper) then
         write(stdout, '(a,f11.6,a,f11.6,a)') 'Upper bound: ', mu_upper_bound, &
              ' (N = ', N_upper_bound, ')'
      else
         write(stdout, '(a,f11.6,a,f11.6)') 'Upper bound for mu not yet found'
      end if
      write(stdout, '(a,f11.6,a,f11.6)') 'Updating mu from ', muvec(size(muvec)), ' to ', &
         dmft_mu_updated_by_bisection
   end if

 end function dmft_mu_updated_by_bisection

 subroutine dmft_shoot_next_mu(np, muvec, Nvec, maxdiff, mixing, target_N, &
      reshoot, chem_shift)

   !---------------------------------------------------------------------------!
   ! Calculates the next chemical potential via Newton's rootfinding method    !
   !   receives : vectorNmu(nmu_step-np + 1:nmu_step)                          !
   !   receives : vectormu(nmu_step-np + 1:nmu_step)                           !
   !   root finding : f = (N(mu) - target_N)                                   !
   !   Newton : x_{n + 1} - x_{n} = - f(x_n) / f'(x_n + 1)                     !
   !---------------------------------------------------------------------------!

   use comms,  only: pub_my_proc_id, pub_on_root
   use linalg, only: linalg_invert_serial
   use rundat, only: pub_debug_on_root
   use utils,  only: utils_alloc_check, utils_dealloc_check

   implicit none

   ! Arguments
   integer,         intent(in   ) :: np
   real(kind = DP), intent(in   ) :: muvec(:)
   real(kind = DP), intent(in   ) :: Nvec(:)
   real(kind = DP), intent(in   ) :: maxdiff
   real(kind = DP), intent(in   ) :: mixing
   real(kind = DP), intent(in   ) :: target_N
   logical,         intent(in   ) :: reshoot
   real(kind = DP), intent(  out) :: chem_shift

   ! Local variables
   integer                      :: nmu_step
   real(kind = DP)              :: x0, fp0
   real(kind = DP), allocatable :: vectorNmu(:), vectorMu(:), vecy(:)
   real(kind = DP), allocatable :: mat(:,:)
   real(kind = DP)              :: f(100), x(100), fp(100)
   real(kind = DP)              :: tt1, tt2, tt
   integer                      :: i, j, id
   real(kind = DP)              :: yn, fpyn, zn, fzn
   integer                      :: info
   integer                      :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_shoot_next_mu'

   nmu_step = np
   allocate(vectorNmu(size(Nvec)), stat = ierr)
   call utils_alloc_check('dmft_shoot_next_mu', 'vectorNmu', ierr)
   allocate(vectorMu(size(muvec)), stat = ierr)
   call utils_alloc_check('dmft_shoot_next_mu', 'vectorMu', ierr)
   allocate(vecy(np), stat = ierr)
   call utils_alloc_check('dmft_shoot_next_mu', 'vecy', ierr)
   allocate(mat(np, np), stat = ierr)
   call utils_alloc_check('dmft_shoot_next_mu', 'mat', ierr)

   vectorNmu = Nvec
   vectorMu = muvec

   if (np == 1) then
      if (vectorNmu(nmu_step) > target_N) then
         chem_shift = -0.003_DP
      else
         chem_shift = 0.003_DP
      end if
      deallocate(vectorNmu, stat = ierr)
      call utils_dealloc_check('dmft_shoot_next_mu', 'vectorNmu', ierr)
      deallocate(vectorMu, stat = ierr)
      call utils_dealloc_check('dmft_shoot_next_mu', 'vectorMu', ierr)
      deallocate(vecy, stat = ierr)
      call utils_dealloc_check('dmft_shoot_next_mu', 'vecy', ierr)
      deallocate(mat, stat = ierr)
      call utils_dealloc_check('dmft_shoot_next_mu', 'mat', ierr)
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_shoot_next_mu'
      return
   end if

   do i = 1, np
      f(i) = vectorNmu(nmu_step-(np-i))-target_N
      x(i) = vectorMu(nmu_step-(np-i))
   enddo

   x0 = x(np)
   fp0 = f(np)

   if (abs(x(2) - x(1)) < 0.00000001_DP) then
      chem_shift = 0.0_DP ! you got there...
      if (pub_on_root) write(stdout, *) 'chem_shift = ', chem_shift
      deallocate(vectorNmu, stat = ierr)
      call utils_dealloc_check('dmft_shoot_next_mu', 'vectorNmu', ierr)
      deallocate(vectorMu, stat = ierr)
      call utils_dealloc_check('dmft_shoot_next_mu', 'vectorMu', ierr)
      deallocate(vecy, stat = ierr)
      call utils_dealloc_check('dmft_shoot_next_mu', 'vecy', ierr)
      deallocate(mat, stat = ierr)
      call utils_dealloc_check('dmft_shoot_next_mu', 'mat', ierr)
      if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
           &dmft_shoot_next_mu'
      return
   endif

   if (np == 2) then
      fp(1) = (f(2)-f(1)) / (x(2)-x(1))
   endif

   if(np > 2)then
      mat = 0.0_DP
      vecy = 0.0_DP
      do i = 1, np
         do j = 1, np
            mat(j, i) = x(j)**(np-i)
         enddo
         vecy(i) = f(i)
      enddo
      call linalg_invert_serial(mat, info)
      vecy = matmul(mat, vecy)
   endif

   if(np > 2)then
      do i = 1, np-1
         fp(i) = internal_derivate(np, vecy, i, x0)
      enddo
      if(fp(1) < -1.0E-4_DP)then
         if (pub_debug_on_root) write(stdout, *) 'WARNING: negative derivative &
              &in N(mu) - only use two points'
         fp = 0.0_DP
         fp(1) = ( f(np)-f(np-1) ) / ( x(np)-x(np-1) )
      endif
   endif

   ! if(pub_on_root)          write(stdout, *)  'Position0: ', x0
   ! if(pub_on_root)          write(stdout, *)  'function0: ', fp0
   ! do i = 1, np
   !    if(pub_on_root)          write(stdout, *)  'Position : ', x(i)
   !    if(pub_on_root)          write(stdout, *)  'Function : ', f(i)
   !    if(i < np .and. pub_on_root) write(stdout, *)  'Der      : ', fp(i)
   ! enddo

   if(np == 2)then
      if(abs(fp(1)) < 50.0_DP)then
         chem_shift = sign(1.0_DP, fp(1)) * ( - fp0 ) * 0.01_DP
      else
         chem_shift =           mixing * ( - fp0 )  / fp(1)
      endif
      ! if (pub_debug_on_root) write(stdout, *) 'np = ', np, '; setting chem_shift = &
      !      &', chem_shift
   endif

   if(np == 3)then
      tt = fp0*fp(2) / (2.0_DP*fp(1))
      chem_shift = - fp0 / ( fp(1) - tt )
      ! if (pub_debug_on_root) write(stdout, *) 'np = ', np, '; setting chem_shift = &
      !      &', chem_shift
   endif

   if(np >= 4)then
      tt2 = (fp0**2) * ( 3.0_DP*(fp(2)**2)-fp(1)*fp(3) ) / ( 6.0_DP* (fp(1)**4) )
      tt1 = 1.0_DP + fp0*fp(2) / ( 2.0_DP*(fp(1)**2) ) + tt2
      chem_shift = - fp0 / fp(1) * tt1
      ! if (pub_debug_on_root) write(stdout, *) 'np = ', np, '; setting chem_shift = &
      !      &', chem_shift
   endif

   if(np > 2 .and. ( (fp0 > 0.0_DP .and. chem_shift > 0.0_DP) .or. &
        (fp0 < 0.0_DP .and.chem_shift < 0.0_DP) ) )then
      if (pub_debug_on_root) then
         write(stdout, *) 'WARNING going in wrong direction, switch back to &
              &first derivative only'
         ! write(stdout, *) ' fp0                 : ', fp0
         ! write(stdout, *) ' obtained chem_shift : ', chem_shift
      end if
      fp = 0.0_DP
      fp(1) = (f(np)-f(np-1)) / (x(np)-x(np-1))
      if(abs(fp(1)) < 50.0_DP)then
         chem_shift = sign(1.0_DP, fp(1)) * ( - fp0 ) * 0.01_DP
      else
         chem_shift =             mixing * ( - fp0 )  / fp(1)
      endif
      ! if (pub_on_root) write(stdout, *) 'corrected chem_shift : ', chem_shift
   endif

   if(reshoot .and. np > 2)then
      yn = x0 - fp0/fp(1)
      fpyn = internal_derivate(np, vecy, 1, yn)
      zn = x0 - 2.0_DP*fp0/( fp(1) + fpyn )
      fzn = internal_derivate(np, vecy, 0, zn)
      chem_shift = zn - (fpyn + fp(1))/(3.0_DP*fpyn - fp(1)) * fzn / fp(1)
      chem_shift = (chem_shift - x0)
      if (pub_on_root) write(stdout, *) 'Reshoot, setting chem_shift = ', &
           chem_shift
   endif

   if(chem_shift > maxdiff) then
      chem_shift = maxdiff
      ! if (pub_debug_on_root) write(stdout, *) 'chemshift > maxdiff, setting &
      !      &chem_shift = maxdiff = ', chem_shift
   end if

   if(chem_shift < -maxdiff) then
      chem_shift = -maxdiff
      ! if (pub_debug_on_root) write(stdout, *) 'chemshift < -maxdiff, setting &
      !      &chem_shift = -maxdiff = ', chem_shift
   end if

   deallocate(vectorNmu, stat = ierr)
   call utils_dealloc_check('dmft_shoot_next_mu', 'vectorNmu', ierr)
   deallocate(vectorMu, stat = ierr)
   call utils_dealloc_check('dmft_shoot_next_mu', 'vectorMu', ierr)
   deallocate(vecy, stat = ierr)
   call utils_dealloc_check('dmft_shoot_next_mu', 'vecy', ierr)
   deallocate(mat, stat = ierr)
   call utils_dealloc_check('dmft_shoot_next_mu', 'mat', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_shoot_next_mu'

 end subroutine

 real(kind = DP) function internal_derivate(n, vec, order, x0)

   !---------------------------------------------------------------------------!
   ! Description:                                                              !
   !---------------------------------------------------------------------------!

   implicit none

   ! Arguments
   integer,         intent(in   ) :: n
   real(kind = DP), intent(in   ) :: vec(n)
   integer,         intent(in   ) :: order
   real(kind = DP), intent(in   ) :: x0

   ! Local variables
   real(kind = DP) :: coef
   integer         :: iorder
   integer         :: derorder
   integer         :: i

   internal_derivate = 0.0_DP
   do i = 1, n
      iorder = n-i
      coef = 1.0_DP
      do derorder = 1, order
         coef = coef*real(iorder-derorder + 1, kind=DP)
      enddo
      if(iorder-order >= 0) internal_derivate = internal_derivate + &
           coef*(x0**(iorder-order))*vec(i)
   enddo

 end function

 integer function lines_in_file(i_unit)

   !---------------------------------------------------------------------------!
   ! Counts the number of lines in a file                                      !
   !---------------------------------------------------------------------------!

   implicit none

   ! Arguments
   integer, intent(in   ) :: i_unit

   ! Local variables
   integer :: iostatus

   rewind(i_unit)
   lines_in_file = 0
   do
      read(i_unit, *, iostat = iostatus)
      if (iostatus < 0) exit
      lines_in_file = lines_in_file + 1
   enddo
   rewind(i_unit)
 end function

 function diag(cc)

   !---------------------------------------------------------------------------!
   ! Returns the diagonal of a real 2D matrix                                  !
   !---------------------------------------------------------------------------!

   implicit none

   ! Arguments
   real(kind = DP), intent(in   ) :: cc(:,:)

   ! Local variables
   real(kind = DP) :: diag(size(cc, 1))
   integer         :: i

   do i = 1, size(cc, 1)
      diag(i) = cc(i, i)
   enddo
 end function

 pure real(kind = DP) function internal_dexpc(rr)

   !---------------------------------------------------------------------------!
   ! A protected exponential                                                   !
   !---------------------------------------------------------------------------!

   implicit none

   ! Arguments
   real(kind = DP), intent(in) :: rr

   ! Local variables
   real(kind = DP), parameter :: MAX_EXP = 700.0_DP
   real(kind = DP), parameter :: MIN_EXP = -700.0_DP

   if(rr < MAX_EXP) then
      if(rr < MIN_EXP)then
         internal_dexpc = 0.0_DP
      else
         internal_dexpc = EXP(rr)
      endif
   else
      internal_dexpc = EXP(MAX_EXP)
   endif

 end function

 real(kind = DP) function internal_norm_vector(xx)

   !---------------------------------------------------------------------------!
   ! Returns the norm of a given vector                                        !
   !---------------------------------------------------------------------------!

   implicit none

   ! Arguments
   real(kind = DP), intent(in   ) :: xx(3)

   internal_norm_vector = sqrt(sum(xx(:)**2))

 end function

 function internal_scalprod(xx, yy)

   !---------------------------------------------------------------------------!
   ! Returns the scalar product of two real vectors                            !
   !---------------------------------------------------------------------------!

   implicit none

   ! Arguments
   real(kind = DP), intent(in), dimension(:) :: xx, yy

   ! Local variables
   real(kind = DP) :: internal_scalprod
   integer         :: j

   internal_scalprod = 0.0_DP
   do j = 1, size(xx)
      internal_scalprod = internal_scalprod + xx(j)*yy(j)
   enddo

 end function

 integer function StrInt2(ch)

   !---------------------------------------------------------------------------!
   ! Converts a character to the corresponding integer                         !
   !---------------------------------------------------------------------------!

   implicit none

   ! Local variables
   integer :: i, j, ilen

   character(len=*), intent(in) :: ch
   j = 0

   ilen = LEN_TRIM(ch)
   do i = 1, ilen ! + 1 BUG corrected
      if (i.eq.ilen) then
         j = j + IACHAR(ch(i:i))-48
      else
         j = j + INT((IACHAR(ch(i:i))-48)*10**(ilen-i))
      end if
   end do
   StrInt2 = j
 end function

 subroutine internal_check_rotation_transpose(chan, rot_vec_angle, hat)

   !---------------------------------------------------------------------------!
   ! Checks a rotation by ensuring the transpose gives the reverse             !
   ! transformation                                                            !
   !---------------------------------------------------------------------------!

   use hubbard_build, only: hubbard_cmp_rot
   use utils,         only: utils_assert, utils_alloc_check, utils_dealloc_check
   use rundat,        only: pub_debug_on_root

   implicit none

   ! Arguments
   integer,         intent(in   ) :: chan
   integer,         intent(in   ) :: hat
   real(kind = DP), intent(in   ) :: rot_vec_angle(3, 3)

   ! Local variables
   integer                      :: i
   real(kind = DP), allocatable :: rotation(:,:)
   real(kind = DP), allocatable :: rotation_t(:,:)
   real(kind = DP), allocatable :: mat2(:,:)
   real(kind = DP)              :: mat(3, 3)
   integer                      :: ierr

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &internal_check_rotation_transpose'

   allocate(rotation(chan, chan), stat = ierr)
   call utils_alloc_check('internal_check_rotation_transpose', 'rotation', ierr)
   allocate(rotation_t(chan, chan), stat = ierr)
   call utils_alloc_check('internal_check_rotation_transpose', 'rotation_t', &
        ierr)
   allocate(mat2(chan, chan), stat = ierr)
   call utils_alloc_check('internal_check_rotation_transpose', 'mat2', ierr)

   rotation = hubbard_cmp_rot( rot_vec_angle, (chan-1)/2, chan, hat)
   rotation_t = hubbard_cmp_rot( transpose(rot_vec_angle), (chan-1)/2, chan, &
        hat)

   mat = matmul( rot_vec_angle, transpose(rot_vec_angle) )

   do i = 1, 3
      call utils_assert(abs(mat(i, i)-1.0_DP) < 1.0E-5_DP, "Error in &
           &internal_check_rotation_transpose: rotation diagonals incorrect")
   enddo

   call utils_assert(maxval(abs(rotation_t-transpose(rotation))) < 1.0E-5_DP, &
        "Error in internal_check_rotation_transpose: transpose incorrect")

   mat2 = matmul( rotation_t, rotation )

   do i = 1, chan
      call utils_assert(abs(mat2(i, i)-1.0_DP) < 1.0E-5_DP, "Error in &
           &internal_check_rotation_transpose: cubic harmonic issues")
   enddo

   deallocate(rotation, stat = ierr)
   call utils_dealloc_check('internal_check_rotation_transpose', 'rotation', &
        ierr)
   deallocate(rotation_t, stat = ierr)
   call utils_dealloc_check('internal_check_rotation_transpose', 'rotation_t', &
        ierr)
   deallocate(mat2, stat = ierr)
   call utils_dealloc_check('internal_check_rotation_transpose', 'mat2', ierr)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &internal_check_rotation_transpose'

 end subroutine

 subroutine dmft_rotate_to_local_basis(green, occupancy1, occupancy2, chan, &
      rot_vec_angle1, rot_vec_angle2, hat, hatb, rotateit)

   !---------------------------------------------------------------------------!
   ! Rotates the green function into the local basis specified by two rotation !
   ! angles, or alternatively in order to diagonalise the occupancy matrices   !
   ! of two Hubbard atoms. Often this will be done for just one Hubbard atom   !
   ! by setting occupancy2 = occupancy1                                        !
   !---------------------------------------------------------------------------!

   use hubbard_build, only: hubbard_cmp_rot
   use linalg,        only: linalg_dsyev_lt
   use utils,         only: utils_unit, utils_int_to_str
   use rundat,        only: pub_debug_on_root, pub_dmft_write

   implicit none

   ! Arguments
   complex(kind = DP), intent(inout) :: green(:,:)
   real(kind = DP),    intent(in   ) :: occupancy1(:,:)
   real(kind = DP),    intent(in   ) :: occupancy2(:,:)
   integer,            intent(in   ) :: chan
   real(kind = DP),    intent(in   ) :: rot_vec_angle1(:,:)
   real(kind = DP),    intent(in   ) :: rot_vec_angle2(:,:)
   integer,            intent(in   ) :: hat
   integer,            intent(in   ) :: hatb
   logical,            intent(in   ) :: rotateit

   ! Local variables
   logical         :: read_in_file
   real(kind = DP) :: rotation1(chan, chan)
   real(kind = DP) :: rotation2(chan, chan)
   real(kind = DP) :: diagdens(chan)
   logical         :: check
   real(kind = DP) :: rotation2_spin(2*chan, 2*chan)
   integer         :: funit

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_rotate_to_local_basis'

   read_in_file = maxval(abs(rot_vec_angle1)) > 1.0E-4_DP .and. &
        maxval(abs(rot_vec_angle2)) > 1.0E-4_DP

   if(.not.read_in_file)then
      rotation1 = occupancy1
      rotation2 = occupancy2
      call linalg_dsyev_lt(rotation1, diagdens, chan)
      call linalg_dsyev_lt(rotation2, diagdens, chan)
   else
      call internal_check_rotation_transpose(chan, rot_vec_angle1, hat )
      call internal_check_rotation_transpose(chan, rot_vec_angle2, hatb)
      rotation1 = hubbard_cmp_rot(rot_vec_angle1, (chan-1)/2, chan, hat)
      rotation2 = hubbard_cmp_rot(rot_vec_angle2, (chan-1)/2, chan, hatb)
   endif

   if(hat == hatb .and. rotateit)then
      inquire(file = 'mask_loc_rot_atom_pm_' // trim(adjustl(utils_int_to_str(&
           hat))), exist = check)

      if(.not.check .and. pub_dmft_write)then
         funit = utils_unit()
         open(unit = funit, file = 'mask_loc_rot_atom_pm_' // &
              trim(adjustl(utils_int_to_str( hat ))), form = 'unformatted')
         write(funit) shape(rotation2)
         write(funit)       rotation2
         close(funit)
      endif

      inquire(file = 'mask_loc_rot_atom_spin_' // trim(adjustl( &
           utils_int_to_str( hat ))), exist = check)
      rotation2_spin                              = 0.0_DP
      rotation2_spin(     1:  chan,     1:  chan) = rotation2
      rotation2_spin(chan + 1:2*chan, chan + 1:2*chan) = rotation2

      if(.not.check .and. pub_dmft_write)then
         open(unit = funit, file = 'mask_loc_rot_atom_spin_' // &
              trim(adjustl(utils_int_to_str( hat ))), form = 'unformatted')
         write(funit) shape(rotation2_spin)
         write(funit)       rotation2_spin
         close(funit)
      endif

   endif

   green = matmul(matmul(transpose(rotation1), green), rotation2)

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_rotate_to_local_basis'

 end subroutine

 subroutine dmft_rotate_to_cart_basis(sigma, occupancy1, occupancy2, chan, &
      rot_vec_angle1, rot_vec_angle2, hat, hatb)

   !---------------------------------------------------------------------------!
   ! Rotates the self energy sigma back to the cartesian basis                 !
   !---------------------------------------------------------------------------!

   use hubbard_build, only: hubbard_cmp_rot
   use linalg,        only: linalg_dsyev_lt
   use rundat,        only: pub_debug_on_root

   implicit none

   ! Arguments
   complex(kind = DP), intent(inout) :: sigma(:,:)
   real(kind = DP),    intent(in   ) :: occupancy1(:,:)
   real(kind = DP),    intent(in   ) :: occupancy2(:,:)
   real(kind = DP),    intent(in   ) :: rot_vec_angle1(:,:)
   real(kind = DP),    intent(in   ) :: rot_vec_angle2(:,:)
   integer,            intent(in   ) :: chan
   integer,            intent(in   ) :: hat
   integer,            intent(in   ) :: hatb

   ! Local variables
   real(kind = DP) :: diagdens(chan)
   real(kind = DP) :: rotation1(chan, chan)
   real(kind = DP) :: rotation2(chan, chan)
   logical         :: read_in_file

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: entering &
        &dmft_rotate_to_cart_basis'

   read_in_file = maxval(abs(rot_vec_angle1)) > 1.0E-4_DP .and. &
        maxval(abs(rot_vec_angle2)) > 1.0E-4_DP

   if (.not.read_in_file) then
      ! ebl: ensuring not to overwrite the occupancy matrices
      rotation1 = occupancy1
      rotation2 = occupancy2
      call linalg_dsyev_lt(rotation1, diagdens, chan)
      call linalg_dsyev_lt(rotation2, diagdens, chan)
      ! ebl: these replace the old calls:
      ! call eigenvector_matrix(chan, occupancy1, diagdens, rotation1)
      ! call eigenvector_matrix(chan, occupancy2, diagdens, rotation2)
   else
      rotation1 = hubbard_cmp_rot(rot_vec_angle1, (chan-1)/2, chan, hat)
      rotation2 = hubbard_cmp_rot(rot_vec_angle2, (chan-1)/2, chan, hatb)
   endif

   sigma = matmul(matmul(rotation1, sigma), transpose(rotation2))

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &dmft_rotate_to_cart_basis'

 end subroutine

 function diago_occupancy_func(tt, chan, rotateit, hat, rot_vec_angle)

   !---------------------------------------------------------------------------!
   ! Return the diagonal of a (potentially rotated) matrix                     !
   !---------------------------------------------------------------------------!

   use hubbard_build, only: hubbard_cmp_rot
   use linalg,        only: linalg_dsyev_lt
   use rundat,        only: pub_debug_on_root

   implicit none

   ! Arguments
   real(kind = DP),           intent(in   ) :: tt(:,:)
   integer,                   intent(in   ) :: chan
   logical,                   intent(in   ) :: rotateit
   integer,                   intent(in   ) :: hat
   real(kind = DP), optional, intent(in   ) :: rot_vec_angle(3, 3)

   ! Local variables
   real(kind = DP) :: diago_occupancy_func(chan, chan)
   real(kind = DP) :: diagdens(chan)
   real(kind = DP) :: rotation(chan, chan)
   logical         :: read_in_file

   read_in_file = .false.

   if(present(rot_vec_angle))then
      read_in_file = maxval(abs(rot_vec_angle)) > 1.0E-4_DP
   endif

   if(.not.read_in_file)then
      rotation = (tt + transpose(tt))/2.0_DP
      call linalg_dsyev_lt(rotation, diagdens, chan)
      ! call eigenvector_matrix(chan, tt_, diagdens, rotation)
   else
      rotation = hubbard_cmp_rot(rot_vec_angle, (chan-1)/2, chan, hat)
   endif

   if(rotateit)then
      diago_occupancy_func = matmul(matmul(transpose(rotation), tt), rotation)
   else
      diago_occupancy_func = tt
   endif

   if (pub_debug_on_root) write(stdout, '(a)') 'DEBUG: leaving &
        &diago_occupancy_func'
   return

 end function

 integer function count_files(filename1, mpi_size)

   !---------------------------------------------------------------------------!
   ! Counts the number of files, assuming naming follows an index              !
   !---------------------------------------------------------------------------!

   use utils, only: utils_int_to_str

   implicit none

   ! Arguments
   character(*), intent(in   ) :: filename1
   integer,      intent(in   ) :: mpi_size

   ! Local variables
   integer :: i
   logical :: check

   count_files = 0
   do i = 1, mpi_size
      inquire(file = trim(adjustl(filename1)) // trim(adjustl( &
           utils_int_to_str(i))), exist = check)
      if(check) count_files = count_files + 1
   enddo
 end function

 !=============================================================================!
 ! KERNEL DIIS ROUTINES                                                        !
 !=============================================================================!
 ! Adapted by Cedric Weber as hubbard_build_mod_dmft_diis.h, as a modified     !
 ! version of the pre-existing routines                                        !
 ! Split off from hubbard_build mod and revised by Edward Linscott March 2016  !
 ! Ultimately, if possible, we would like to be using the standard ONETEP      !
 ! kernel_diis mod in place of these routines                                  !
 !=============================================================================!

 subroutine dmft_kernel_diis_mix(next_dkn_in, dkn_in, dkn_out, residues, &
      overlap, iter, ientry, num_spins, wont_compute_residues)

   !---------------------------------------------------------------------------!
   ! Selects method for kernel DIIS.                                           !
   ! Written by Alvaro Ruiz Serrano in November 2010.                          !
   !---------------------------------------------------------------------------!

   use comms,  only: pub_on_root
   use rundat, only: pub_kernel_diis, pub_kernel_diis_linear_iter, &
        pub_kernel_diis_scheme, pub_kernel_diis_size, pub_num_spins
   use sparse, only: sparse_copy
   use utils,  only: utils_abort

   implicit none

   ! Arguments
   type(SPAM3), intent(inout) :: next_dkn_in(pub_num_spins)
   type(SPAM3), intent(inout) :: dkn_in(:,:)
   type(SPAM3), intent(inout) :: dkn_out(:,:)
   type(SPAM3), intent(inout) :: residues(:,:)
   type(SPAM3), intent(in   ) :: overlap
   integer,     intent(in   ) :: iter
   integer,     intent(in   ) :: ientry
   integer,     intent(in   ) :: num_spins
   logical,     intent(in   ) :: wont_compute_residues

   ! Local variables
   character(len = 40) :: diis_error_message
   integer             :: is

   select case(pub_kernel_diis_scheme)
   case('diag', 'DIAG')
      ! ars: do diagonalisation only
      do is = 1, pub_num_spins
         call sparse_copy(next_dkn_in(is), dkn_out(is, ientry))
      end do

   case('dkn_linear', 'DKN_LINEAR')
      ! ars: do linear mixing only
      call dmft_kernel_diis_linear_mixing(next_dkn_in, dkn_in, dkn_out, ientry)

   case('dkn_pulay', 'DKN_PULAY')
      !ars: do Pulay mixing. Optional linear mixing for the first iterations.
      if(iter.le.pub_kernel_diis_linear_iter) then
         call dmft_kernel_diis_linear_mixing(next_dkn_in, dkn_in, dkn_out, &
              ientry)
      else
         call dmft_kernel_diis_pulay_mixing(next_dkn_in, dkn_out, residues, &
              overlap, ientry, num_spins, wont_compute_residues)
      end if

      case default
      ! if(pub_on_root) write(stdout, *) "Unknown mixing method. ONETEP stops."
      call utils_abort('Unkown mixing scheme for DIIS')

   end select

 end subroutine

 subroutine dmft_kernel_diis_linear_mixing(next_dkn_in, dkn_in, dkn_out, ientry)

   !---------------------------------------------------------------------------!
   ! Performs linear mixing of density kernels.                                !
   ! Written by Alvaro Ruiz Serrano in June 2010.                              !
   !---------------------------------------------------------------------------!

   use comms,  only: pub_on_root
   use rundat, only: pub_kernel_diis, pub_kernel_diis_coeff, &
        pub_kernel_diis_size, pub_num_spins, pub_output_detail

   implicit none

   ! Arguments
   type(SPAM3), intent(inout) :: next_dkn_in(pub_num_spins)
   type(SPAM3), intent(inout) :: dkn_in(:,:)
   type(SPAM3), intent(inout) :: dkn_out(:,:)
   integer,     intent(in   ) :: ientry

   ! Local variables
   real(kind = DP) :: dmft_kernel_diis_c_out

   dmft_kernel_diis_c_out = 1.0_DP - pub_kernel_diis_coeff

   ! if(pub_on_root .and. pub_output_detail.ge.VERBOSE) then
   !    write(stdout, '(/a)') "... PERFORMING LINEAR MIXING ..."
   !    write(stdout, '(a9, f6.4, a10, f6.4/)') " C_in = ", &
   !         pub_kernel_diis_coeff, ", C_out = ", dmft_kernel_diis_c_out
   ! end if

   ! ars: mix kernels linearly
   call dmft_kernel_diis_linear_mix_kernels(next_dkn_in, dkn_out, dkn_in, &
        ientry)

 end subroutine


 subroutine dmft_kernel_diis_linear_mix_kernels(next_dkn_in, dkn_out, dkn_in, &
      ientry)

   !---------------------------------------------------------------------------!
   ! This subroutine performs a linear mixing of the density kernels.          !
   ! Written by Alvaro Ruiz Serrano in June 2010.                              !
   !---------------------------------------------------------------------------!

   use rundat, only: pub_kernel_diis, pub_kernel_diis_coeff, &
        pub_kernel_diis_size, pub_num_spins
   use sparse, only: sparse_axpy, sparse_scale
   use rundat, only: pub_kernel_diis, pub_kernel_diis_coeff, pub_num_spins

   implicit none

   ! Arguments
   type(SPAM3), intent(inout) :: next_dkn_in(pub_num_spins)
   type(SPAM3), intent(inout) :: dkn_out(:,:)
   type(SPAM3), intent(inout) :: dkn_in(:,:)
   integer,     intent(in   ) :: ientry

   ! Local variables
   integer         :: is
   real(kind = DP) :: dmft_kernel_diis_c_out

   dmft_kernel_diis_c_out = 1.0_DP - pub_kernel_diis_coeff

   do is = 1, pub_num_spins
      ! ars: K_in(n + 1) = 0
      call sparse_scale(next_dkn_in(is), 0.0_DP)
      ! ars: K_in(n + 1) = c_out*K_out(n)
      call sparse_axpy(next_dkn_in(is), dkn_out(is, ientry), &
           dmft_kernel_diis_c_out)
      ! ars: K_in(n + 1) = c_out*K_out(n) + c_in*K_in(n)
      call sparse_axpy(next_dkn_in(is), dkn_in(is, ientry), &
           pub_kernel_diis_coeff)
   end do

 end subroutine

 subroutine dmft_kernel_diis_pulay_mixing(next_dkn_in, dkn_out, residues, &
      overlap, ientry, num_spins, wont_compute_residues)

   !---------------------------------------------------------------------------!
   ! Performs Pulay mixing of density kernels.                                 !
   ! Based on Pulay, Chem. Phys. Lett. 73, 2, 1980.                            !
   ! Written by Alvaro Ruiz Serrano in June 2010.                              !
   ! Updated by Alvaro Ruiz Serrano in October 2010 to combine linear + Pulay  !
   ! mixing schemes.                                                           !
   !---------------------------------------------------------------------------!

   use comms,  only: pub_on_root
   use rundat, only: pub_kernel_diis_size, pub_num_spins, pub_output_detail
   use utils,  only: utils_alloc_check, utils_dealloc_check, utils_int_to_str
   use rundat, only: pub_num_spins

   implicit none

   ! Arguments
   type(SPAM3), intent(inout) :: next_dkn_in(pub_num_spins)
   type(SPAM3), intent(inout) :: dkn_out(:,:)
   type(SPAM3), intent(inout) :: residues(:,:)
   type(SPAM3), intent(in)    :: overlap
   integer,     intent(in)    :: ientry
   integer,     intent(in   ) :: num_spins
   logical,     intent(in   ) :: wont_compute_residues

   ! Local variables
   real(kind = DP), allocatable :: Bmat(:,:,:)!ars: linear system to be solved
   real(kind = DP), allocatable :: coeffs(:,:)!ars: Pulay mixing coeffs
   integer                      :: ierr

   if(pub_on_root .and. pub_output_detail.ge.VERBOSE) then
      write(stdout, '(/, a)') "... PERFORMING PULAY MIXING ..."
   end if

   ! ars: allocate Bmat and coeffs
   allocate(Bmat(1:pub_num_spins, 1:ientry + 1, 1:ientry + 1), stat = ierr)
   call utils_alloc_check('dmft_kernel_diis_pulay_mixing', 'Bmat', ierr)
   allocate(coeffs(1:pub_num_spins, 1:ientry + 1), stat = ierr)
   call utils_alloc_check('dmft_kernel_diis_pulay_mixing', 'coeffs', ierr)

   ! ars: calculate B matrix.
   call dmft_kernel_diis_pulay_calc_Bmatrix(Bmat, residues, overlap, ientry, &
        num_spins, wont_compute_residues)

   ! ars: solve linear system Bd = 0 (including Lagrange mult)
   call dmft_kernel_diis_pulay_find_coeffs(coeffs, Bmat, ientry)
   if(pub_output_detail.ge.VERBOSE) call &
        dmft_kernel_diis_pulay_coeffs_summ(coeffs, ientry)

   ! ars: mix kernel
   call dmft_kernel_diis_pulay_mix_kernels(next_dkn_in, dkn_out, coeffs, ientry)

   ! ars: deallocate Bmat and coeffs
   deallocate(coeffs, stat = ierr)
   call utils_dealloc_check('dmft_kernel_diis_pulay_mixing', 'coeffs', ierr)
   deallocate(Bmat, stat = ierr)
   call utils_dealloc_check('dmft_kernel_diis_pulay_mixing', 'Bmat', ierr)

 end subroutine

 subroutine dmft_kernel_diis_pulay_calc_Bmatrix(Bmat, residues, overlap, &
      ientry, num_spins, wont_compute_residues)

   !---------------------------------------------------------------------------!
   ! This subroutine calculates matrix B of the system of linear               !
   ! equations to be solved in order to mix the density kernels                !
   ! according to the Pulay method.                                            !
   ! Written by Alvaro Ruiz Serrano in April 2010.                             !
   !---------------------------------------------------------------------------!

   use comms,  only: comms_barrier, pub_my_proc_id, pub_on_root
   use rundat, only: pub_debug, pub_debug_on_root, pub_kernel_diis_size, &
        pub_num_spins
   use sparse, only: sparse_create, sparse_destroy, sparse_product, &
        sparse_trace, sparse_transpose
   use utils,  only: utils_abort

   implicit none

   ! Arguments
   real(kind = DP), intent(inout) :: Bmat(:,:,:)
   type(SPAM3),     intent(inout) :: residues(:,:)
   type(SPAM3),     intent(in   ) :: overlap
   integer,         intent(in   ) :: ientry
   integer,         intent(in   ) :: num_spins
   logical,         intent(in   ) :: wont_compute_residues

   ! Local variables
   real(kind = DP) :: Bmat_tmp(ientry, ientry)
   integer         :: table(ientry*ientry*num_spins, 3), kkk, kk, iini, istep
   integer         :: is, ii, jj
   type(SPAM3)     :: spam3buffer1, spam3buffer2, spam3buffer3

   call sparse_create(spam3buffer1, residues(1, 1), overlap)
   call sparse_create(spam3buffer2, residues(1, 1))
   call sparse_create(spam3buffer3, spam3buffer1, spam3buffer2)

   kkk = 0
   do jj = 1, ientry
      do ii = 1, ientry
         do is = 1, num_spins
            kkk = kkk + 1
            table(kkk, 1) = jj
            table(kkk, 2) = ii
            table(kkk, 3) = is
         enddo
      enddo
   enddo

   iini = 0
   istep = 1
   Bmat = 0.0_DP

   if(.not. pub_dmft_use_gpu_partially)then
      iini = pub_dmft_mpi_rank-1
      istep = pub_dmft_mpi_size
   else
      if(pub_dmft_mpi_rank /= 1 .and. wont_compute_residues) then
         write(*, *) 'ERROR : WRITING KERNEL WITH 1 GPU : MPI_RANK = ', &
              pub_dmft_mpi_rank
         write(*, *) 'WAS NOT SUPPOSED TO BE HERE'
         call utils_abort('Error: Child MPI thread erroneously attempting to &
              &write')
      endif
      ! if(pub_dmft_mpi_rank /= 1) goto 8888
   endif

   if (.not.(pub_dmft_mpi_rank /= 1 .and. pub_dmft_use_gpu_partially)) then
      do kk = iini + 1, kkk, istep

         if(pub_my_proc_id == 0) write(*, '(a, i4, 2i6)') 'MPI_RANK SCANNING : &
              &', pub_dmft_mpi_rank, kk, kkk

         jj = table(kk, 1)
         ii = table(kk, 2)
         is = table(kk, 3)

#ifndef debugMIXv3
         if(pub_dmft_use_gpu_onlyme) then
            call dmft_matmul_in_dense_format_r(spam3buffer1, residues(is, ii), &
                 overlap)
         else
#endif
            ! ars: spam3buffer1 = R_i x S
            call sparse_product(spam3buffer1, residues(is, ii), overlap)
#ifndef debugMIXv3
         endif
#endif
         ! ars: spam3buffer2 = (R_j)^t

         call sparse_transpose(spam3buffer2, residues(is, jj))

#ifndef debugMIXv3
         if(pub_dmft_use_gpu_onlyme)then
            call dmft_matmul_in_dense_format_r(spam3buffer3, spam3buffer1, &
                 spam3buffer2)
         else
#endif
            ! ars: spam3buffer3 = (R_i) x S x (R_j)^t
            call sparse_product(spam3buffer3, spam3buffer1, spam3buffer2)
#ifndef debugMIXv3
         endif
#endif
         ! ars: Bmat(is, ii, jj) = tr[(R_i) x S x (R_j)^t x S]
         Bmat(is, ii, jj) = sparse_trace(spam3buffer3, overlap)

      end do
   end if

   if(.not.wont_compute_residues)then
      if(pub_dmft_mpi_size /= 1)then
         is = 1
         Bmat_tmp(1:ientry, 1:ientry) = Bmat(is, 1:ientry, 1:ientry)
         write(*, *) 'MPI_RANK CALLING NFS SUM_B : ', pub_dmft_mpi_rank, &
              pub_dmft_mpi_size
         write(*, *) 'IENTRY SIZE                : ', ientry
         call utils_abort("Reached point in code where nfs sync used to &
              &appear; the code needs to be updated to use these settings.")
         ! call sync_via_nfs(Bmat_tmp, 'sumBup', ientry, .true., .true.)
         write(*, *) 'SYNC DONE, RANK : ', pub_dmft_mpi_rank, pub_dmft_mpi_size
         Bmat(is, 1:ientry, 1:ientry) = Bmat_tmp(1:ientry, 1:ientry)
         if(num_spins == 2)then
            is = 2
            Bmat_tmp(1:ientry, 1:ientry) = Bmat(is, 1:ientry, 1:ientry)
            write(*, *) 'MPI_RANK CALLING NFS SUM_B : ', pub_dmft_mpi_rank, &
                 pub_dmft_mpi_size
            write(*, *) 'IENTRY SIZE                : ', ientry
            call utils_abort("Reached point in code where nfs sync used to &
                 &appear; the code needs to be updated to use these settings.")
            ! call sync_via_nfs(Bmat_tmp, 'sumBdown', ientry, .true., .true.)
            write(*, *) 'SYNC DONE, RANK : ', pub_dmft_mpi_rank, &
                 pub_dmft_mpi_size
            Bmat(is, 1:ientry, 1:ientry) = Bmat_tmp(1:ientry, 1:ientry)
         endif
      endif
   endif

   if(num_spins /= pub_num_spins)then
      if(pub_num_spins == 1)then
         write(*, *) 'ERROR,   num_spins = ', num_spins
         write(*, *) 'pub_num_spins  = ', pub_num_spins
         call utils_abort("ERROR in dmft_kernel_diis_pulay_calc_Bmatrix: &
              &num_spins/=pub_num_spins")
      endif
      Bmat(2, :, :) = Bmat(1, :, :)
   endif

   ! ars: destroy SPAM3 buffers
   call sparse_destroy(spam3buffer1)
   call sparse_destroy(spam3buffer2)
   call sparse_destroy(spam3buffer3)

   ! ars: add space for Langrange multiplier
   Bmat(:, 1:ientry, ientry + 1) = -1.0_DP
   Bmat(:, ientry + 1, 1:ientry) = -1.0_DP
   Bmat(:, ientry + 1, ientry + 1) = 0.0_DP

 end subroutine dmft_kernel_diis_pulay_calc_Bmatrix

 subroutine dmft_kernel_diis_pulay_find_coeffs(coeffs, Bmat, ientry)

   !---------------------------------------------------------------------------!
   ! This subroutine finds the coefficients for the Pulay mixing               !
   ! after solving a system of linear equations Bd = 0.                        !
   ! Written by Alvaro Ruiz Serrano in May 2010.                               !
   !---------------------------------------------------------------------------!

   use comms,  only: comms_barrier, pub_my_proc_id, pub_on_root
   use rundat, only: pub_num_spins

   implicit none

   ! Arguments
   real(kind = DP), intent(inout) :: coeffs(:,:)
   real(kind = DP), intent(in   ) :: Bmat(:,:,:)
   integer,         intent(in   ) :: ientry

   ! Local variables
   integer            :: INFO, LDA, LDB, M, N, NRHS
   integer            :: IPIV(1:ientry + 1)
   real(kind = DP)    :: A(ientry + 1, ientry + 1)
   real(kind = DP)    :: B(ientry + 1)
   character(LEN = 1) :: TRANS
   integer            :: is

   ! LAPACK subroutine
   external :: dgetrf
   external :: dgetrs

   ! ars: library wrapper variables
   M = ientry + 1
   N = ientry + 1
   LDA = ientry + 1
   ! ars: set up parameters for LAPACK DGETRS
   TRANS = 'N'
   NRHS = 1
   LDB = ientry + 1

   ! ars: solve linear system
   do is = 1, pub_num_spins

      ! ars: call LAPACK DGETRF
      A = Bmat(is, :, :)
      !#ifdef debug
      ! call write_array( A, ' matrix Bcoef ' // trim(adjustl(&
      !      utils_int_to_str(ientry + 1))) )
      !#endif
      call DGETRF(M, N, A, LDA, IPIV, INFO)
      if((INFO.ne.0) .and. pub_on_root) write(stdout, *) "Error calling &
           &DGETRF. INFO|proc = ", INFO, pub_my_proc_id

      ! ars: call LAPACK DGETRS
      B(1:ientry) = 0.0_DP
      B(ientry + 1) = -1.0_DP
      call DGETRS(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      if((INFO.ne.0) .and. pub_on_root) write(stdout, *) "Error calling &
           &DGETRS. INFO|proc = ", INFO, pub_my_proc_id

      ! ars: set coeffs(is) before exit
      coeffs(is, :) = B(:)

   end do

   call comms_barrier

   if(pub_dmft_impose_same_coeffs .and. pub_num_spins == 2) coeffs(2, :) = &
        coeffs(1, :)

 end subroutine

 subroutine dmft_kernel_diis_pulay_coeffs_summ(coeffs, ientry)

   !---------------------------------------------------------------------------!
   ! This subroutine prints a summary of the coefficients that are             !
   ! used during the density kernel mixing.                                    !
   ! Written by Alvaro Ruiz Serrano in July 2010.                              !
   !---------------------------------------------------------------------------!

   use comms,  only: pub_on_root
   use rundat, only: pub_num_spins

   implicit none

   ! Arguments
   real(kind = DP), intent(in   ) :: coeffs(:,:)
   integer,         intent(in   ) :: ientry

   ! Local variables
   integer :: is, it

   if(pub_on_root) then
      do it = 1, ientry + 1
         do is = 1, pub_num_spins

            if(it.le.ientry) then
               write(stdout, '(a6, i1, a1, i2, a4, f20.12)') "Coeff(", is, ", &
                    &", it, ") = ", coeffs(is, it)
            else
               write(stdout, '(a7, i1, a6, f20.12, /)') "Lambda(", is, ") = ", &
                    coeffs(is, it)
            end if

         end do
      end do
   end if

 end subroutine

 subroutine dmft_kernel_diis_pulay_mix_kernels(next_dkn_in, dkn_out, coeffs, &
      ientry)

   !---------------------------------------------------------------------------!
   ! This subroutine mixes the kernels according to the Pulay method.          !
   ! Written by Alvaro Ruiz Serrano in May 2010.                               !
   !---------------------------------------------------------------------------!

   use rundat, only: pub_kernel_diis_size, pub_num_spins
   use sparse, only: sparse_axpy, sparse_scale

   implicit none

   ! Arguments
   type(SPAM3),     intent(inout) :: next_dkn_in(pub_num_spins)
   type(SPAM3),     intent(in   ) :: dkn_out(:,:)
   real(kind = DP), intent(in   ) :: coeffs(:,:)
   integer,         intent(in   ) :: ientry

   ! Local variables
   integer :: is, counter

   do is = 1, pub_num_spins
      call sparse_scale(next_dkn_in(is), 0.0_DP)
   end do

   ! ars: K_in(n + 1) = sum_i K_out(i) * d_i
   do counter = 1, ientry
      do is = 1, pub_num_spins
         call sparse_axpy(next_dkn_in(is), dkn_out(is, counter), coeffs(is, &
              counter))
      end do
   end do


 end subroutine

 subroutine dmft_kernel_diis_init(residues, next_dkn_in, dkn_out, dkn_in, &
      denskern)

   !---------------------------------------------------------------------------!
   ! This subroutine initialises the first element of dkn_in by copying the    !
   ! latest density kernel onto it (coming from the last NGWF iteration).      !
   ! Written by Alvaro Ruiz Serrano in April 2010.                             !
   ! Adapted by Edward Linscott to use SPAM3_ARRAYS                            !
   !---------------------------------------------------------------------------!

   use comms,        only: pub_on_root
   use rundat,       only: pub_kernel_diis, pub_kernel_diis_size, pub_num_spins
   use sparse,       only: sparse_copy, sparse_scale
   use sparse_array, only: sparse_array_scale, SPAM3_ARRAY
   use kernel,       only: DKERN

   implicit none

   ! Arguments
   type(SPAM3_ARRAY), intent(inout) :: residues
   type(SPAM3),       intent(inout) :: next_dkn_in(pub_num_spins)
   type(SPAM3_ARRAY), intent(inout) :: dkn_out
   type(SPAM3_ARRAY), intent(inout) :: dkn_in
   type(DKERN),       intent(in)    :: denskern

   ! Local variables
   integer :: is

   call sparse_array_scale(residues, 0.0_DP)
   call sparse_array_scale(dkn_in, 0.0_DP)
   call sparse_array_scale(dkn_out, 0.0_DP)
   ! ars: initialise first ientry of dkn_in
   do is = 1, pub_num_spins
      call sparse_scale(next_dkn_in(is), 0.0_DP)
      call sparse_copy(dkn_in%m(is, 1), denskern%kern%m(is, PUB_1K)%p)
   end do

 end subroutine

 subroutine dmft_kernel_diis_residue_inout(residue, kern_out, kern_in)

   !---------------------------------------------------------------------------!
   ! This subroutine calculates the density kernel residue after each          !
   ! kernel DIIS iteration. R_i = K^out_i - K^in_i.                            !
   ! Written by Alvaro Ruiz Serrano in April 2010.                             !
   !---------------------------------------------------------------------------!

   use comms,  only: pub_on_root
   use rundat, only: pub_num_spins
   use sparse, only: sparse_axpy, sparse_copy

   implicit none

   ! Arguments
   type(SPAM3), intent(inout) :: residue(pub_num_spins)
   type(SPAM3), intent(in   ) :: kern_out(pub_num_spins)
   type(SPAM3), intent(in   ) :: kern_in(pub_num_spins)

   ! Local variables
   integer :: is

   do is = 1, pub_num_spins
      call sparse_copy(residue(is), kern_out(is))
      ! calculate [R = K^out - K^in]
      call sparse_axpy(residue(is), kern_in(is), -1.0_DP)
   end do

 end subroutine

 subroutine dmft_kernel_diis_shift(shifted, dkn_in, dkn_out, residues, iter)

   !---------------------------------------------------------------------------!
   ! This subroutine shifts the matrices in arrays in case iter has reached    !
   ! pub_kernel_diis_size value.                                               !
   ! Written by Alvaro Ruiz Serrano in May 2010.                               !
   !---------------------------------------------------------------------------!

   use rundat, only: pub_kernel_diis, pub_kernel_diis_size, pub_num_spins
   use sparse, only: sparse_copy

   implicit none

   ! Arguments
   logical,     intent(  out) :: shifted
   type(SPAM3), intent(inout) :: dkn_in(:,:)
   type(SPAM3), intent(inout) :: dkn_out(:,:)
   type(SPAM3), intent(inout) :: residues(:,:)
   integer,     intent(in   ) :: iter

   ! Local variables
   integer :: scntr, is

   shifted = .false.

   ! ars: move old entries one position back
   if (iter.ge.pub_kernel_diis_size) then
      shifted = .true.
      do scntr = 1, pub_kernel_diis_size-1
         do is = 1, pub_num_spins
            call sparse_copy(dkn_in(is, scntr), dkn_in(is, scntr + 1))
            call sparse_copy(dkn_out(is, scntr), dkn_out(is, scntr + 1))
            call sparse_copy(residues(is, scntr), residues(is, scntr + 1))
         end do
      end do
   end if

 end subroutine

 subroutine dmft_kernel_diis_find_ientry(ientry, iter, shifted)

   !---------------------------------------------------------------------------!
   ! This subroutine finds the next ientry to work with.                       !
   ! Written by Alvaro Ruiz Serrano in May 2010.                               !
   !---------------------------------------------------------------------------!

   use comms,  only: comms_barrier
   use rundat, only: pub_kernel_diis, pub_kernel_diis_size

   implicit none

   ! Arguments
   integer, intent(  out) :: ientry
   integer, intent(in   ) :: iter
   logical, intent(in   ) :: shifted

   ! ars: find next position
   if (shifted) then
      ientry = pub_kernel_diis_size
   else
      ientry = iter + 1
   end if

   call comms_barrier

 end subroutine

end module dmft
