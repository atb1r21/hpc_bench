! -*- mode: F90 ; mode: font-lock ; column-numbero-mode: true -*-

!==============================================================================!
! This module contains the subroutines to perform self-consistent energy       !
! minimisation in the inner loop of ONETEP in the following manner:            !
!                                                                              !
! 1) Diagonalise the Hamiltonian.                                              !
! 2) Build density kernel from the eigenvalues (enforce idempotency).          !
! 3) Perform density kernel DIIS to achieve self-consistency.                  !
!                                                                              !
! Density kernel DIIS is an alternative to LNV in the inner loop of ONETEP.    !
! Hamiltonian diagonalisation only is allowed although is unstable for large   !
! systems.                                                                     !
!------------------------------------------------------------------------------!
! *** Please note that the Hamiltonian diagonalisation and kernel DIIS         !
! functionalities are currently under development.                             !
!------------------------------------------------------------------------------!
! Written by Alvaro Ruiz Serrano in April 2010.                                !
! Modified by Alvaro Ruiz Serrano in April 2012 to add ODA and LiSTi.          !
! Added level shifting by Peter J. Cherry in May 2012.                         !
! Added Hamiltonian mixing by Alvaro Ruiz Serrano in May 2012.                 !
! Modified for embedding by Joseph Prentice, May 2018                          !
!==============================================================================!


module kernel_diis

  use constants, only: DP
  use dense, only: DEM
  use ngwf_representation, only: NGWF_HAM
  use sparse, only: SPAM3
  use sparse_embed, only: SPAM3_EMBED
#ifdef GPU_PGI
#ifdef GPU_SP_TEST
  use fourier_gpu_wrapper_mod
#endif
#endif

  implicit none

  public :: kernel_diis_calculate
  public :: kernel_diis_lagrangian
  public :: kernel_diis_build_pq
  public :: kernel_diis_build_idemp_dkn
  public :: kernel_diis_mu

  ! ars: the following variables are accesible to all routines within the module

  ! ars: logical flags to set task
  logical, private :: do_ls, do_oda, do_list, do_listi, do_listb, &
       do_pulay, allow_linear, mix_dkn, mix_ham, do_diag

  ! ars: logical flags from pub_devel_code
  logical, private :: ls_dkn_eigs

  ! ars: integer flag to set scheme
  integer, private :: scheme

  ! ars: array size of the DIIS history
  integer, private :: diis_size

  ! ars: number of Hamiltonian/dkn eigenstates
  integer, private :: num

  ! ars: call counter
  integer, save, private :: ncalls = 0
  character(len=20), private :: ncalls_char

  ! ars: parameter for level shifting (changes during runtime)
  real(kind=DP) :: beta

  ! ars: integer parameters to define schemes from input file
  integer, private, parameter :: DIAG       = 1
  integer, private, parameter :: DKN_LINEAR = 2
  integer, private, parameter :: DKN_PULAY  = 3
  integer, private, parameter :: DKN_LISTI  = 4
  integer, private, parameter :: DKN_LISTB  = 5
  integer, private, parameter :: HAM_LINEAR = 6
  integer, private, parameter :: HAM_PULAY  = 7
  integer, private, parameter :: HAM_LISTI  = 8
  integer, private, parameter :: HAM_LISTB  = 9

  ! ars: integer parameters for generalised schemes
  integer, private, parameter :: DKN_LIST  = 100
  integer, private, parameter :: HAM_LIST  = 200


  ! ars: and now the arrays necessary for the different schemes - made accesible
  !      for all the routine within the module to avoid code replication.

  ! ars: output density kernel (always used)
  type(SPAM3_EMBED), private, allocatable :: dkn_out(:)

  ! ars: output Hamiltonian (always used)
  type(NGWF_HAM), private :: ham_out

  ! ars: history of dkn matrices (used for Pulay or LIST mixing of dkn)
  type(SPAM3_EMBED), private, allocatable :: dkn_his(:,:)

  ! ars: history of ham matrices (used for Pulay or LIST mixing of ham)
  type(SPAM3_EMBED), private, allocatable :: ham_his(:,:)

  ! ars: history of dkn error matrices (used for Pulay mixing of dkn)
  type(SPAM3_EMBED), private, allocatable :: dkn_err_his(:,:)

  ! ars: history of ham error matrices (used for Pulay mixing of ham)
  type(SPAM3_EMBED), private, allocatable :: ham_err_his(:,:)

  ! ars: history of dkn differentials (used for LIST)
  type(SPAM3_EMBED), private, allocatable :: dkn_diff_his(:,:)

  ! ars: history of ham differentials (used for LIST)
  type(SPAM3_EMBED), private, allocatable :: ham_diff_his(:,:)

  ! ars: Hamiltonian/dkn eigenvectors (used for level shifting)
  type(DEM), private, allocatable :: eigvecs(:)

  ! ars: quantities to account for convergence
  real(kind=DP), private, allocatable :: commut(:), idemp(:), trks(:), residual(:)
  real(kind=DP), private, allocatable :: homo(:), lumo(:), gap(:), delta_gap(:)
  real(kind=DP), private :: deltaE

  ! ars: iteration counter
  integer, private :: iter



contains


  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !
  !                           DRIVER ROUTINE
  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine kernel_diis_calculate(total_energy, mu, pur_denskern, &
       denskern, ham, lhxc_fine, ngwf_basis, hub_proj_basis, hub, rep, &
       mdl, hfxstate)


    use comms, only: pub_on_root, pub_total_num_procs
    use constants, only: stdout, DP, max_spins, VERBOSE
#ifdef GPU_SP_TEST
    use fourier_gpu_wrapper_mod, only: GPU_SP, E_gpu_previous, deltaE_gpu, &
        GPU_SP_switched
#endif
    use function_basis, only: FUNC_BASIS
    use kernel, only: kernel_rescale_spam3
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_num_spins, pub_output_detail, pub_print_qc, &
         pub_num_kpoints, PUB_1K
    use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_copy, &
         sparse_embed_array_copy, sparse_embed_array_num_rows
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_assert, utils_banner

    implicit none

    ! Arguments
    real(kind=DP),    intent(  out) :: total_energy
    real(kind=DP),    intent(  out) :: mu(max_spins)
    type(SPAM3_EMBED_ARRAY), intent(inout) :: pur_denskern
    type(SPAM3_EMBED_ARRAY), intent(inout) :: denskern
    type(NGWF_HAM),   intent(inout) :: ham
    type(MODEL),      intent(in   ) :: mdl
    type(HFX_STATE),  intent(inout), target :: hfxstate
    real(kind=DP),    intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins,mdl%nsub)
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in   ) :: hub_proj_basis(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(NGWF_REP),   intent(in   ) :: rep

    real(kind=dp),dimension(1:pub_num_spins,1:pub_num_kpoints) :: t_occ
    ! ars: local variables
    integer :: is
    ! agrecocmplx
    logical :: loc_cmplx


    !--------------------------------------------------------------------------
    ! ars: start-up

    ! ars: start timer
    call timer_clock('kernel_diis_calculate',1)

    ! ars: print header
    if (pub_on_root) then
       write(stdout,'(/,a80)') utils_banner('-', &
            'Hamiltonian diagonalisation + kernel DIIS')

       ! ars: print name of the linear algebra routine
#ifdef SCALAPACK
       if(pub_output_detail.ge.VERBOSE) &
            write(stdout,'(a21, i4, a7)') "Running ScaLAPACK on ", &
            pub_total_num_procs, " procs."
#else
       if(pub_output_detail.ge.VERBOSE) &
            write(stdout,'(a21, i4, a45)') "Running LAPACK on ", &
            pub_total_num_procs, " procs. Diagonalisation will occur in serial."
#endif
    end if

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(1)%iscmplx

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine kernel_diis_calculate not ready yet for more&
         & than one k-point.')

    ! ars: check input parameters
    call kernel_diis_check_params

    ! ars: set maximum number of eigenstates to be calculated
    num = int(sparse_embed_array_num_rows(denskern))


    ! ars: count how many times we entered this routine
    ncalls = ncalls + 1
    write(ncalls_char,*) ncalls
    ncalls_char=trim(adjustl(ncalls_char))


    ! ars: rescale input kernel to keep the number of electrons constant
    t_occ = real(rep%n_occ,kind=dp)
    call kernel_rescale_spam3(denskern, rep%overlap, t_occ)


    ! ars: allocate workspace
    call kernel_diis_alloc(rep,mdl)


    ! ars: initialise level shifter if required
    if (do_ls) then

       ! ars: allocate workspace
       ! agrecocmplx: call with optional argument in complex case
       call kernel_diis_ls_alloc(iscmplx=loc_cmplx)

       ! ars: build initial eigenvectors
       call kernel_diis_ls_init(beta,ham%ham,rep%overlap,rep%n_occ,num)

    end if


    ! ars: initialise output variables
    mu(:) = 0.0_DP
    total_energy = 0.0_DP

#ifdef GPU_SP_TEST
    GPU_SP = .true.
    GPU_SP_switched = .false.
    E_gpu_previous = 0.0_DP
    deltaE_gpu = 0.0_DP
#endif

    !--------------------------------------------------------------------------
    ! ars: apply appropriate SCF minimisation scheme

    select case (scheme)
    case (DIAG)

       ! ars: diagonalisation only: no mixing involved
       call kernel_diis_diag_scheme(denskern%m(:,PUB_1K), ham, &
            lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, &
            mdl, hfxstate)

    case (DKN_LINEAR)

       ! ars: linear mix of density kernels
       call kernel_diis_dkn_linear_scheme(denskern%m(:,PUB_1K), &
            ham, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate)

    case (DKN_PULAY)

       ! ars: Pulay mixing of density kernels
       call kernel_diis_dkn_pulay_scheme(denskern%m(:,PUB_1K), &
            ham, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate)

    case (DKN_LIST)

       ! ars: LiST mixing of density kernels
       call kernel_diis_dkn_list_scheme(denskern%m(:,PUB_1K), &
            ham, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate)

    case (HAM_LINEAR)

       ! ars: linear mixing of Hamiltonians (%ham)
       call kernel_diis_ham_linear_scheme(denskern%m(:,PUB_1K), &
            ham, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate)

    case (HAM_PULAY)

       ! ars: Pulay mixing of Hamiltonians (%ham)
       call kernel_diis_ham_pulay_scheme(denskern%m(:,PUB_1K), &
            ham, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate)

    case (HAM_LIST)
       ! ars: LiST mixing of Hamiltonians (%ham)
       call kernel_diis_ham_list_scheme(denskern%m(:,PUB_1K), ham, &
            lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate)

    case default
       call utils_abort("Kernel DIIS scheme not recognised.")
    end select

! kaw: Added to reset back to double precision and remove the need for the dEnergy=0.0 check
#ifdef GPU_PGI
#ifdef GPU_SP_TEST
       GPU_SP=.false.
#endif
#endif


    ! ars: copy dkn_out to denskern and get pur_denskern
    do is = 1, pub_num_spins
       call sparse_embed_copy(denskern%m(is,PUB_1K), dkn_out(is))
    end do
    call sparse_embed_array_copy(pur_denskern, denskern)

    ! ars: build final Hamiltonian, calculate final energy
    call kernel_diis_build_ham(total_energy, ham, &
         denskern%m(:,PUB_1K), lhxc_fine, &
         ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
         .true., .false., .true.)

    ! ars: calculate mu
    call kernel_diis_mu(mu, denskern, ham%ham, rep%overlap)

    ! ars: print QC test if required
    if (pub_print_qc) call kernel_diis_qc_output(total_energy, mu, &
         denskern%m(:,PUB_1K) ,ham%ham)

    ! ars: deallocate workspace for level shifting
    call kernel_diis_ls_dealloc

    ! ars: deallocate workspace
    call kernel_diis_dealloc


    ! ars: stop timer
    call timer_clock('kernel_diis_calculate',2)


  end subroutine kernel_diis_calculate



  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !
  !                           END DRIVER ROUTINE
  !_____________________________________________________________________________
  !_____________________________________________________________________________





  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !
  !                      SCF ENERGY MINIMISATION SCHEMES
  !_____________________________________________________________________________
  !_____________________________________________________________________________




  subroutine kernel_diis_diag_scheme(dkn_in, ham, lhxc_fine, &
       ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate)

    !==========================================================================!
    ! Hamiltonian diagonalisation SCF minimisation scheme:                     !
    ! Warning: unstable unless level-shifting is enabled.                      !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2012 from previous routine         !
    ! kernel_diis_calculate.                                                   !
    ! Modified for embedding by Joseph Prentice, May 2018                      !
    !==========================================================================!

    use constants, only: DP
    use function_basis, only: FUNC_BASIS
#ifdef GPU_PGI
#ifdef GPU_SP_TEST
    use fourier_gpu_wrapper_mod, only: GPU_SP, GPU_SP_switched, GPU_SP_switch, &
                                       deltaE_gpu
#endif
#endif
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_kernel_diis_maxit, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: dkn_in(1:pub_num_spins)
    type(NGWF_HAM), intent(inout) :: ham
    type(MODEL),      intent(in   ) :: mdl
    type(HFX_STATE),  intent(inout) ,target :: hfxstate
    real(kind=DP),    intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub)
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in   ) :: hub_proj_basis(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(NGWF_REP),   intent(in   ) :: rep


    ! real
    real(kind=DP) :: total_energy_out
    real(kind=DP), save :: total_energy_old = 0.0_DP
    real(kind=DP), save :: gap_old_1 = 0.0_DP
    real(kind=DP), save :: gap_old_2 = 0.0_DP

    ! logical
    logical :: converged

    ! integer
    integer :: is

    ! character
    character(len=8) :: okmessage


    ! ars: initilise variables
    converged = .false.

    ! ars: build initial dens_dep matrices with dkn_in and calculate total energy
    call kernel_diis_build_ham(total_energy_out, ham, dkn_in, lhxc_fine, &
         ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
         .true., .false., .true.)


    ! ars: start kernel_diis loop (ONETEP inner loop)
    kernel_diis_loop:    do iter=1,pub_kernel_diis_maxit

#ifdef GPU_SP_TEST
    ! kaw: Accuracy switching toggles
    if ((.not.GPU_SP_switched).and.(abs(deltaE_gpu).gt.(GPU_SP_switch))) then
       GPU_SP=.true.                                                    ! Turn on SP code
    else
       GPU_SP=.false.                                                   ! kaw: Switch to DP
       if (.not.(iter.eq.1)) GPU_SP_switched=.true.
    end if
    if (iter.eq.1) GPU_SP=.true.
#endif

       !----------------------------------------------------------
       ! ars: SCF

       ! ars: diagonalise Hamiltonian and build dkn (out)
       if (do_ls) then
          ! pjc: Go into level shifting proceedure if appropriate
          call kernel_diis_level_shifter(dkn_out, dkn_in, ham%ham, &
               rep%overlap, rep%n_occ)
       else
          ! ars: no level shifter -> diag ham and build kernel
          call kernel_diis_build_idemp_dkn(dkn_out, ham%ham, &
               rep%overlap, rep%n_occ, homo, lumo, gap)
       end if

       ! ars: update dens_dep matrices with dkn_out and calculate total energy
       call kernel_diis_build_ham(total_energy_out, ham, dkn_out, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
            .true., .false., .true.)


       !-----------------------------------------------------------
       ! ars: calculate energy, gap, commutator, idemp. and ne.

       ! ars: calculate delta in energy
       deltaE = total_energy_out - total_energy_old
       total_energy_old = total_energy_out

       ! ars: calculate delta in the gap
       delta_gap(1) = gap(1) - gap_old_1
       if (pub_num_spins.eq.2) delta_gap(2) = gap(2) - gap_old_2
       gap_old_1 = gap(1)
       if (pub_num_spins.eq.2) gap_old_2 = gap(2)

       ! ars: calculate residual
       call kernel_diis_residual_inout(residual, dkn_out, dkn_in, rep%overlap)

       ! ars: calculate commutators and number of electrons
       call kernel_diis_constants(commut, idemp, trks, dkn_out, dkn_in, &
            ham%ham, rep%overlap, rep%inv_overlap)


       !-----------------------------------------------------------
       ! ars: check convergence, print line and exit if appropriate

       ! ars: check convergence
       call kernel_diis_conv_check(converged, okmessage)

       ! ars: print info line
       call kernel_diis_info_line(converged, okmessage, total_energy_out)

       ! ars: exit if appropiate
       if(converged.or.(iter.eq.pub_kernel_diis_maxit)) then
          exit kernel_diis_loop
       end if
       !----------------------------------------------------------




       !----------------------------------------------------------
       ! ars: diagonalisation: K_in (i+1) = K_out(i)
       do is=1, pub_num_spins
          call sparse_embed_copy(dkn_in(is), dkn_out(is))
       end do
       !----------------------------------------------------------


    end do kernel_diis_loop


    !----


  end subroutine kernel_diis_diag_scheme



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_dkn_linear_scheme(dkn_in, ham_in, lhxc_fine,&
       ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate)

    !==========================================================================!
    ! Density kernel linear mixing SCF energy minimisation scheme.             !
    ! Warning: potentially unstable unless level-shifting or ODA are enabled.  !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2012 from previous routine         !
    ! kernel_diis_calculate.                                                   !
    ! Modified for embedding by Joseph Prentice, May 2018                      !
    !==========================================================================!

    use constants, only: DP
#ifdef GPU_PGI
#ifdef GPU_SP_TEST
    use fourier_gpu_wrapper_mod, only: GPU_SP, GPU_SP_switched, GPU_SP_switch, &
                                       deltaE_gpu
#endif
#endif
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_kernel_diis_maxit, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3_EMBED),      intent(inout) :: dkn_in(1:pub_num_spins)
    type(NGWF_HAM),   intent(inout) :: ham_in
    type(MODEL),      intent(in   ) :: mdl
    type(HFX_STATE),  intent(inout), target :: hfxstate
    real(kind=DP),    intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub)
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in   ) :: hub_proj_basis(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(NGWF_REP),   intent(in   ) :: rep


    ! real
    real(kind=DP) :: total_energy_in, total_energy_out
    real(kind=DP), save :: total_energy_old = 0.0_DP
    real(kind=DP), save :: gap_old_1 = 0.0_DP
    real(kind=DP), save :: gap_old_2 = 0.0_DP

    ! logical
    logical :: converged

    ! character
    character(len=8) :: okmessage

    ! ars: initilise variables
    converged = .false.


    ! ars: start kernel_diis loop (ONETEP inner loop)
    kernel_diis_loop:    do iter=1,pub_kernel_diis_maxit

#ifdef GPU_SP_TEST
    ! kaw: Accuracy switching toggles
    if ((.not.GPU_SP_switched).and.(abs(deltaE_gpu).gt.(GPU_SP_switch))) then
       GPU_SP=.true.                                                    ! Turn on SP code
    else
       GPU_SP=.false.                                                   ! kaw: Switch to DP
       if (.not.(iter.eq.1)) GPU_SP_switched=.true.
    end if
    if (iter.eq.1) GPU_SP=.true.
#endif

       !----------------------------------------------------------
       ! ars: SCF

       ! ars: update dens_dep matrices with dkn_in and calculate total energy
       call kernel_diis_build_ham(total_energy_in, ham_in, dkn_in, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
            .true., .false., .true.)

       ! ars: diagonalise Hamiltonian and build dkn (out)
       if (do_ls) then
          ! pjc: Go into level shifting proceedure if appropriate
          call kernel_diis_level_shifter(dkn_out, dkn_in, ham_in%ham, &
               rep%overlap, rep%n_occ)
       else
          ! ars: no level shifter -> diag ham and build kernel
          call kernel_diis_build_idemp_dkn(dkn_out, ham_in%ham, &
               rep%overlap, rep%n_occ, homo, lumo, gap)
       end if

       ! ars: update dens_dep matrices with dkn_out and calculate total energy
       call kernel_diis_build_ham(total_energy_out, ham_out, dkn_out, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
            .true., .false., .true.)



       !-----------------------------------------------------------
       ! ars: calculate energy, gap, commutator, idemp. and ne.

       ! ars: calculate delta in energy
       deltaE = total_energy_out - total_energy_old
       total_energy_old = total_energy_out

       ! ars: calculate delta in the gap
       delta_gap(1) = gap(1) - gap_old_1
       if (pub_num_spins.eq.2) delta_gap(2) = gap(2) - gap_old_2
       gap_old_1 = gap(1)
       if (pub_num_spins.eq.2) gap_old_2 = gap(2)

       ! ars: calculate residual
       call kernel_diis_residual_inout(residual, dkn_out, dkn_in, rep%overlap)

       ! ars: calculate commutators and number of electrons
       call kernel_diis_constants(commut, idemp, trks, dkn_out, dkn_in, &
            ham_out%ham, rep%overlap, rep%inv_overlap)


       !----------------------------------------------------------
       ! ars: check convergence, print line and exit if appropriate

       ! ars: check convergence
       call kernel_diis_conv_check(converged, okmessage)

       ! ars: print info line
       call kernel_diis_info_line(converged, okmessage, total_energy_out)

       ! ars: exit if appropiate
       if(converged.or.(iter.eq.pub_kernel_diis_maxit)) then
          ! ars: print banner and exit
          exit kernel_diis_loop
       end if
       !----------------------------------------------------------





       !----------------------------------------------------------
       ! ars: linear mixing: K_in (i+1) = (1-C)*K_in(i) + C*K_out(i)
       if (do_oda) then
          call kernel_diis_linear_mix(dkn_in, dkn_out, dkn_in, dkn_out,&
               ham_in%ham, ham_out%ham, total_energy_in, total_energy_out)
       else
          call kernel_diis_linear_mix(dkn_in, dkn_out)
       end if


    end do kernel_diis_loop


    !----


  end subroutine kernel_diis_dkn_linear_scheme



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_dkn_pulay_scheme(dkn_in, ham_in, lhxc_fine,&
       ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate)

    !==========================================================================!
    ! Density kernel Pulay DIIS mixing SCF energy minimisation scheme.         !
    ! Allows to choose error vector from the following:                        !
    !  1) [K_out - K_in]                                                       !
    ! Warning: unstable for large systems unless level-shifting is enabled.    !
    ! Warning: under development - use with caution.                           !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2012 from previous routine         !
    ! kernel_diis_calculate.                                                   !
    ! Modified for embedding by Joseph Prentice, May 2018                      !
    !==========================================================================!

    use constants, only: DP
#ifdef GPU_PGI
#ifdef GPU_SP_TEST
    use fourier_gpu_wrapper_mod, only: GPU_SP, GPU_SP_switched, GPU_SP_switch, &
                                       deltaE_gpu
#endif
#endif
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_kernel_diis_maxit, pub_kernel_diis_linear_iter, &
         pub_num_spins, pub_inner_loop_iteration
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3_EMBED),      intent(inout) :: dkn_in(1:pub_num_spins)
    type(NGWF_HAM),   intent(inout) :: ham_in
    type(MODEL),      intent(in   ) :: mdl
    type(HFX_STATE),  intent(inout), target :: hfxstate
    real(kind=DP),    intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub)
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in   ) :: hub_proj_basis(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(NGWF_REP),   intent(in   ) :: rep



    ! real
    real(kind=DP) :: total_energy_in, total_energy_out
    real(kind=DP), save :: total_energy_old = 0.0_DP
    real(kind=DP), save :: gap_old_1 = 0.0_DP
    real(kind=DP), save :: gap_old_2 = 0.0_DP

    ! logical
    logical :: converged, shifted


    ! integer
    integer :: is, ientry

    ! character
    character(len=8) :: okmessage

    ! ars: initilise variables
    converged = .false.
    shifted   = .false.
    ientry  = 1




    ! ars: start kernel_diis loop (ONETEP inner loop)
    kernel_diis_loop:    do iter=1,pub_kernel_diis_maxit
       ! jd: Globally accessible copy of iteration number, only used for contro-
       !     lling devel_code dumps of grid quantities from various modules
       pub_inner_loop_iteration = iter


#ifdef GPU_SP_TEST
    ! kaw: Accuracy switching toggles
    if ((.not.GPU_SP_switched).and.(abs(deltaE_gpu).gt.(GPU_SP_switch))) then
       GPU_SP=.true.                                                    ! Turn on SP code
    else
       GPU_SP=.false.                                                   ! kaw: Switch to DP
       if (.not.(iter.eq.1)) GPU_SP_switched=.true.
    end if
    if (iter.eq.1) GPU_SP=.true.
#endif


       !----------------------------------------------------------
       ! ars: SCF

       ! ars: update dens_dep matrices with dkn_in and calculate total energy
       call kernel_diis_build_ham(total_energy_in, ham_in, dkn_in, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
            .true., .false., .true.)

       ! ars: diagonalise Hamiltonian and build dkn (out)
       if (do_ls) then
          ! pjc: Go into level shifting proceedure if appropriate
          call kernel_diis_level_shifter(dkn_out, dkn_in, ham_in%ham, &
               rep%overlap, rep%n_occ)
       else
          ! ars: no level shifter -> diag ham and build kernel
          call kernel_diis_build_idemp_dkn(dkn_out, ham_in%ham, &
               rep%overlap, rep%n_occ, homo, lumo, gap)
       end if

       ! ars: update dens_dep matrices with dkn_out and calculate total energy
       call kernel_diis_build_ham(total_energy_out, ham_out, dkn_out, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
            .true., .false., .true.)



       !----------------------------------------------------------
       ! ars: arrange history matrices for Pulay mixing

       ! ars: add matrices to history
       do is = 1, pub_num_spins
          call sparse_embed_copy(dkn_his(ientry,is), dkn_out(is))
       end do

       ! ars: calculate error vector if Pulay |  R = K_out - K_in
       call kernel_diis_diff_inout(dkn_err_his(ientry,:), dkn_out, dkn_in)



       !-----------------------------------------------------------
       ! ars: calculate energy, gap, commutator, idemp. and ne.

       ! ars: calculate delta in energy
       deltaE = total_energy_out - total_energy_old
       total_energy_old = total_energy_out

       ! ars: calculate delta in the gap
       delta_gap(1) = gap(1) - gap_old_1
       if (pub_num_spins.eq.2) delta_gap(2) = gap(2) - gap_old_2
       gap_old_1 = gap(1)
       if (pub_num_spins.eq.2) gap_old_2 = gap(2)

       ! ars: calculate residual
       call kernel_diis_residual_inout(residual, dkn_out, dkn_in, rep%overlap)

       ! ars: calculate commutators and number of electrons
       call kernel_diis_constants(commut, idemp, trks, dkn_out, dkn_in, &
            ham_out%ham, rep%overlap, rep%inv_overlap)


       !----------------------------------------------------------
       ! ars: check convergence, print line and exit if appropriate

       ! ars: check convergence
       call kernel_diis_conv_check(converged, okmessage)

       ! ars: print info line
       call kernel_diis_info_line(converged, okmessage, total_energy_out)

       ! ars: exit if appropiate
       if(converged.or.(iter.eq.pub_kernel_diis_maxit)) then
          ! ars: print banner and exit
          exit kernel_diis_loop
       end if
       !----------------------------------------------------------



       !----------------------------------------------------------
       ! ars: mix kernels
       if (allow_linear.and.(iter.le.pub_kernel_diis_linear_iter)) then

          ! ars: linear mixing: K_in (i+1) = (1-C)*K_in(i) + C*K_out(i)
          if (do_oda) then
             call kernel_diis_linear_mix(dkn_in, dkn_out, dkn_in, dkn_out,&
                  ham_in%ham, ham_out%ham, total_energy_in, total_energy_out)
          else
             call kernel_diis_linear_mix(dkn_in, dkn_out)
          end if

       else

          ! ars: Pulay mixing: K_in (i+1) = sum_j^i C_j*K_out(j)
          call kernel_diis_pulay_mix(dkn_in, dkn_his, dkn_err_his, &
               rep%overlap, ientry)

       end if


       ! ars: shift arrays if we have reached pub_kernel_diis_size iterations
       if (iter.ge.diis_size) then
          shifted = .true.
          call kernel_diis_shift(dkn_his, diis_size, pub_num_spins)
          call kernel_diis_shift(dkn_err_his, diis_size, pub_num_spins)
       end if

       ! ars: find next position in history to work with
       call kernel_diis_find_ientry(ientry, iter, shifted)
       !----------------------------------------------------------


    end do kernel_diis_loop


    !----


  end subroutine kernel_diis_dkn_pulay_scheme



  !_____________________________________________________________________________
  !_____________________________________________________________________________




  subroutine kernel_diis_dkn_list_scheme(dkn_in, ham_in, lhxc_fine, &
       ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate)

    !==========================================================================!
    ! Density kernel Linear Shooting Techniques (LiSTb and LiSTi) SCF energy   !
    ! minimisation schemes.                                                    !
    ! Warning: under development - use with caution.                           !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2012.                              !
    !==========================================================================!

    use constants, only: DP
#ifdef GPU_PGI
#ifdef GPU_SP_TEST
    use fourier_gpu_wrapper_mod, only: GPU_SP, GPU_SP_switched, GPU_SP_switch, &
                                       deltaE_gpu
#endif
#endif
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_kernel_diis_maxit, pub_kernel_diis_linear_iter, &
         pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy, sparse_embed_trace
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3_EMBED),      intent(inout) :: dkn_in(1:pub_num_spins)
    type(NGWF_HAM),   intent(inout) :: ham_in
    type(MODEL),      intent(in   ) :: mdl
    type(HFX_STATE),  intent(inout), target :: hfxstate
    real(kind=DP),    intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub)
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in   ) :: hub_proj_basis(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(NGWF_REP),   intent(in   ) :: rep


    ! real
    real(kind=DP) :: total_energy_in, total_energy_out
    real(kind=DP), save :: total_energy_old = 0.0_DP
    real(kind=DP), save :: gap_old_1 = 0.0_DP
    real(kind=DP), save :: gap_old_2 = 0.0_DP
    real(kind=DP) :: trace

    ! logical
    logical :: converged, shifted


    ! integer
    integer :: is, ientry

    ! character
    character(len=8) :: okmessage

    ! ars: initilise variables
    converged = .false.
    shifted   = .false.
    ientry  = 1


    ! ars: start kernel_diis loop (ONETEP inner loop)
    kernel_diis_loop:    do iter=1,pub_kernel_diis_maxit

#ifdef GPU_SP_TEST
    ! kaw: Accuracy switching toggles
    if ((.not.GPU_SP_switched).and.(abs(deltaE_gpu).gt.(GPU_SP_switch))) then
       GPU_SP=.true.                                                    ! Turn on SP code
    else
       GPU_SP=.false.                                                   ! kaw: Switch to DP
       if (.not.(iter.eq.1)) GPU_SP_switched=.true.
    end if
    if (iter.eq.1) GPU_SP=.true.
#endif

       !----------------------------------------------------------
       ! ars: SCF

       ! ars: update dens_dep matrices with dkn_in and calculate total energy
       call kernel_diis_build_ham(total_energy_in, ham_in, dkn_in, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
            .true., .false., .true.)

       ! ars: diagonalise Hamiltonian and build dkn (out)
       if (do_ls) then
          ! pjc: Go into level shifting proceedure if appropriate
          call kernel_diis_level_shifter(dkn_out, dkn_in, ham_in%ham, &
               rep%overlap, rep%n_occ)
       else
          ! ars: no level shifter -> diag ham and build kernel
          call kernel_diis_build_idemp_dkn(dkn_out, ham_in%ham, &
               rep%overlap, rep%n_occ, homo, lumo, gap)
       end if

       ! ars: update dens_dep matrices with dkn_out and calculate total energy
       call kernel_diis_build_ham(total_energy_out, ham_out, dkn_out, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
            .true., .false., .true.)



       !----------------------------------------------------------
       ! ars: add matrices to history for LIST mixing

       ! ars: add kernel
       do is = 1, pub_num_spins
          call sparse_embed_copy(dkn_his(ientry,is), dkn_out(is))
       end do

       ! ars: add dkn and ham differentials
       call kernel_diis_diff_inout(dkn_diff_his(ientry,:),dkn_out, dkn_in)
       call kernel_diis_diff_inout(ham_diff_his(ientry,:),ham_out%ham, ham_in%ham)

       ! ars: calculate corrected Hohenberg-Kohn-Sham functional
       do is = 1, pub_num_spins
          call sparse_embed_trace(trace, dkn_diff_his(ientry,is), &
               ham_diff_his(ientry,is))
          total_energy_out = total_energy_out + 0.5_DP * trace
       end do



       !-----------------------------------------------------------
       ! ars: calculate energy, gap, commutator, idemp. and ne.

       ! ars: calculate delta in energy
       deltaE = total_energy_out - total_energy_old
       total_energy_old = total_energy_out

       ! ars: calculate delta in the gap
       delta_gap(1) = gap(1) - gap_old_1
       if (pub_num_spins.eq.2) delta_gap(2) = gap(2) - gap_old_2
       gap_old_1 = gap(1)
       if (pub_num_spins.eq.2) gap_old_2 = gap(2)

       ! ars: calculate residual
       call kernel_diis_residual_inout(residual, dkn_out, dkn_in, rep%overlap)

       ! ars: calculate commutators and number of electrons
       call kernel_diis_constants(commut, idemp, trks, dkn_out, dkn_in, &
            ham_out%ham, rep%overlap, rep%inv_overlap)


       !----------------------------------------------------------
       ! ars: check convergence, print line and exit if appropriate

       ! ars: check convergence
       call kernel_diis_conv_check(converged, okmessage)

       ! ars: print info line
       call kernel_diis_info_line(converged, okmessage, total_energy_out)

       ! ars: exit if appropiate
       if(converged.or.(iter.eq.pub_kernel_diis_maxit)) then
          ! ars: print banner and exit
          exit kernel_diis_loop
       end if
       !----------------------------------------------------------





       !----------------------------------------------------------
       ! ars: mix kernels
       if (allow_linear.and.(iter.le.pub_kernel_diis_linear_iter)) then

          ! ars: linear mixing: K_in (i+1) = (1-C)*K_in(i) + C*K_out(i)
          if (do_oda) then
             call kernel_diis_linear_mix(dkn_in, dkn_out, dkn_in, dkn_out,&
                  ham_in%ham, ham_out%ham, total_energy_in, total_energy_out)
          else
             call kernel_diis_linear_mix(dkn_in, dkn_out)
          end if

       else

          ! ars: LiST mixing: K_in (i+1) = sum_j^i C_j*K_out(j)
          call kernel_diis_list_mix(dkn_in,dkn_his,ham_diff_his,dkn_diff_his,ientry)

       end if


       ! ars: shift arrays if we have reached pub_kernel_diis_size iterations
       if (iter.ge.diis_size) then
          shifted = .true.
          call kernel_diis_shift(dkn_his, diis_size, pub_num_spins)
          call kernel_diis_shift(ham_diff_his, diis_size, pub_num_spins)
          call kernel_diis_shift(dkn_diff_his, diis_size, pub_num_spins)
       end if

       ! ars: find next position in history to work with
       call kernel_diis_find_ientry(ientry, iter, shifted)
       !----------------------------------------------------------


    end do kernel_diis_loop


    !----


  end subroutine kernel_diis_dkn_list_scheme


  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_ham_linear_scheme(dkn_in, ham_in, lhxc_fine, &
       ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate)

    !==========================================================================!
    ! Hamiltonian linear mixing SCF energy minimisation scheme.                !
    ! Warning: potentially unstable unless level shifting or ODA are enabled.  !
    ! Warning: under development - use with caution.                           !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2012.                              !
    ! Modified for embedding by Joseph Prentice, May 2018
    !==========================================================================!

    use constants, only: DP
#ifdef GPU_PGI
#ifdef GPU_SP_TEST
    use fourier_gpu_wrapper_mod, only: GPU_SP, GPU_SP_switched, GPU_SP_switch, &
                                       deltaE_gpu
#endif
#endif
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_kernel_diis_maxit, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy
    use timer, only: timer_clock

    implicit none


    ! Arguments
    type(SPAM3_EMBED),      intent(inout) :: dkn_in(1:pub_num_spins)
    type(NGWF_HAM),   intent(inout) :: ham_in
    type(MODEL),      intent(in   ) :: mdl
    type(HFX_STATE),  intent(inout), target :: hfxstate
    real(kind=DP),    intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub)
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in   ) :: hub_proj_basis(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(NGWF_REP),   intent(in   ) :: rep

    ! real
    real(kind=DP) :: total_energy_in, total_energy_out
    real(kind=DP), save :: total_energy_old = 0.0_DP
    real(kind=DP), save :: gap_old_1 = 0.0_DP
    real(kind=DP), save :: gap_old_2 = 0.0_DP

    ! logical
    logical :: converged

    ! integer
    integer :: is

    ! character
    character(len=8) :: okmessage

    ! ars: initilise variables
    converged = .false.



    ! ars: build initial Ham from rescaled kernel and calculate total energy
    call kernel_diis_build_ham(total_energy_in, ham_in, dkn_in, lhxc_fine, &
         ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
         .true., .false., .true.)


    ! ars: start kernel_diis loop (ONETEP inner loop)
    kernel_diis_loop:    do iter=1,pub_kernel_diis_maxit

#ifdef GPU_SP_TEST
    ! kaw: Accuracy switching toggles
    if ((.not.GPU_SP_switched).and.(abs(deltaE_gpu).gt.(GPU_SP_switch))) then
       GPU_SP=.true.                                                    ! Turn on SP code
    else
       GPU_SP=.false.                                                   ! kaw: Switch to DP
       if (.not.(iter.eq.1)) GPU_SP_switched=.true.
    end if
    if (iter.eq.1) GPU_SP=.true.
#endif

       !----------------------------------------------------------
       ! ars: SCF


       ! ars: calculate total energy with incoming Ham and dkn - do not update ham
       call kernel_diis_build_ham(total_energy_in, ham_in, dkn_in, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
            .true., .false., .false.)

       ! ars: diagonalise Hamiltonian and build dkn (out)
       if (do_ls) then
          ! pjc: Go into level shifting proceedure if appropriate
          call kernel_diis_level_shifter(dkn_out, dkn_in, ham_in%ham, &
               rep%overlap, rep%n_occ)
       else
          ! ars: no level shifter -> diag ham and build kernel
          call kernel_diis_build_idemp_dkn(dkn_out, ham_in%ham, &
               rep%overlap, rep%n_occ, homo, lumo, gap)
       end if

       ! ars: update dens_dep matrices with dkn_out and calculate total energy
       call kernel_diis_build_ham(total_energy_out, ham_out, dkn_out, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
            .true., .false., .true.)



       !-----------------------------------------------------------
       ! ars: calculate energy, gap, commutator, idemp. and ne.

       ! ars: calculate delta in energy
       deltaE = total_energy_out - total_energy_old
       total_energy_old = total_energy_out

       ! ars: calculate delta in the gap
       delta_gap(1) = gap(1) - gap_old_1
       if (pub_num_spins.eq.2) delta_gap(2) = gap(2) - gap_old_2
       gap_old_1 = gap(1)
       if (pub_num_spins.eq.2) gap_old_2 = gap(2)

       ! ars: calculate residual
       call kernel_diis_residual_inout(residual, ham_out%ham, ham_in%ham, &
            rep%inv_overlap)

       ! ars: calculate commutators and number of electrons
       call kernel_diis_constants(commut, idemp, trks, dkn_out, dkn_in, &
            ham_out%ham, rep%overlap, rep%inv_overlap)


       !----------------------------------------------------------
       ! ars: check convergence, print line and exit if appropriate

       ! ars: check convergence
       call kernel_diis_conv_check(converged, okmessage)

       ! ars: print info line
       call kernel_diis_info_line(converged, okmessage, total_energy_out)

       ! ars: exit if appropiate
       if(converged.or.(iter.eq.pub_kernel_diis_maxit)) then
          ! ars: print banner and exit
          exit kernel_diis_loop
       end if
       !----------------------------------------------------------





       !----------------------------------------------------------
       ! ars: linear mixing: H_in (i+1) = (1-C)*H_in(i) + C*H_out(i)
       if (do_oda) then
          call kernel_diis_linear_mix(ham_in%ham, ham_out%ham, dkn_in, dkn_out,&
               ham_in%ham, ham_out%ham, total_energy_in, total_energy_out)
       else
          call kernel_diis_linear_mix(ham_in%ham, ham_out%ham)
       end if

       ! ars: copy dkn_out to dkn_in
       do is = 1, pub_num_spins
          call sparse_embed_copy(dkn_in(is), dkn_out(is))
       end do

       !----------------------------------------------------------


    end do kernel_diis_loop


    !----



  end subroutine kernel_diis_ham_linear_scheme



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_ham_pulay_scheme(dkn_in, ham_in, lhxc_fine, &
       ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate)

    !==========================================================================!
    ! Hamiltonian Pulay DIIS mixing SCF energy minimisation scheme.            !
    ! Allows to choose error vector from the following:                        !
    !  1) H_out - H_in                                                         !
    !  2) H_out x K_out x S - S x K_out x H_out                                !
    ! Warning: potentially unstable for large systems.                         !
    ! Warning: under development - use with caution.                           !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2012.                              !
    ! Modifies by Max Phipps in April 2015 for SCF MI error vector             !
    ! Modified for embedding by Joseph Prentice, May 2018                      !
    !==========================================================================!

    use constants, only: DP
#ifdef GPU_PGI
#ifdef GPU_SP_TEST
    use fourier_gpu_wrapper_mod, only: GPU_SP, GPU_SP_switched, GPU_SP_switch, &
         deltaE_gpu
#endif
#endif
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_kernel_diis_maxit, pub_kernel_diis_linear_iter, &
         pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3_EMBED),      intent(inout) :: dkn_in(1:pub_num_spins)
    type(NGWF_HAM),   intent(inout) :: ham_in
    type(MODEL),      intent(in   ) :: mdl
    type(HFX_STATE),  intent(inout), target :: hfxstate
    real(kind=DP),    intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub)
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in   ) :: hub_proj_basis(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(NGWF_REP),   intent(in   ) :: rep

    ! real
    real(kind=DP) :: total_energy_in, total_energy_out
    real(kind=DP), save :: total_energy_old = 0.0_DP
    real(kind=DP), save :: gap_old_1 = 0.0_DP
    real(kind=DP), save :: gap_old_2 = 0.0_DP

    ! logical
    logical :: converged, shifted

    ! integer
    integer :: is, ientry

    ! character
    character(len=8) :: okmessage


    ! ars: initilise variables
    converged = .false.
    shifted   = .false.
    ientry  = 1

    ! ars: build initial Ham from rescaled kernel and calculate total energy
    call kernel_diis_build_ham(total_energy_in, ham_in, dkn_in, lhxc_fine, &
         ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
         .true., .false., .true.)


    ! ars: start kernel_diis loop (ONETEP inner loop)
    kernel_diis_loop:    do iter=1,pub_kernel_diis_maxit

#ifdef GPU_SP_TEST
    ! kaw: Accuracy switching toggles
    if ((.not.GPU_SP_switched).and.(abs(deltaE_gpu).gt.(GPU_SP_switch))) then
       GPU_SP=.true.                                                    ! Turn on SP code
    else
       GPU_SP=.false.                                                   ! kaw: Switch to DP
       if (.not.(iter.eq.1)) GPU_SP_switched=.true.
    end if
    if (iter.eq.1) GPU_SP=.true.
#endif

       !----------------------------------------------------------
       ! ars: SCF

       ! ars: calculate total energy with incoming Ham and dkn - do not update ham
       call kernel_diis_build_ham(total_energy_in, ham_in, dkn_in, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
            .true., .false., .false.)


       ! ars: diagonalise Hamiltonian and build dkn (out)
       if (do_ls) then
          ! pjc: Go into level shifting proceedure if appropriate
          call kernel_diis_level_shifter(dkn_out, dkn_in, ham_in%ham, &
               rep%overlap, rep%n_occ)
       else
          ! ars: no level shifter -> diag ham and build kernel
          call kernel_diis_build_idemp_dkn(dkn_out, ham_in%ham, &
               rep%overlap, rep%n_occ, homo, lumo, gap)
       end if


       ! ars: update dens_dep matrices with dkn_out and calculate total energy
       call kernel_diis_build_ham(total_energy_out, ham_out, dkn_out, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
            .true., .false., .true.)


       !----------------------------------------------------------
       ! ars: add matrices to history for Pulay mixing
       do is = 1, pub_num_spins
          call sparse_embed_copy(ham_his(ientry,is), ham_out%ham(is))
       end do
       ! ars: calculate error vector
       !  E = H_out - H_in
       call kernel_diis_diff_inout(ham_err_his(ientry,:),ham_out%ham, ham_in%ham)

       !-----------------------------------------------------------
       ! ars: calculate energy, gap, commutator, idemp. and ne.

       ! ars: calculate delta in energy
       deltaE = total_energy_out - total_energy_old
       total_energy_old = total_energy_out

       ! ars: calculate delta in the gap
       delta_gap(1) = gap(1) - gap_old_1
       if (pub_num_spins.eq.2) delta_gap(2) = gap(2) - gap_old_2
       gap_old_1 = gap(1)
       if (pub_num_spins.eq.2) gap_old_2 = gap(2)

       ! ars: calculate residual
       call kernel_diis_residual_inout(residual, ham_out%ham, ham_in%ham, &
            rep%inv_overlap)

       ! ars: calculate commutators and number of electrons
       call kernel_diis_constants(commut, idemp, trks, dkn_out, dkn_in, &
            ham_out%ham, rep%overlap, rep%inv_overlap)


       !----------------------------------------------------------
       ! ars: check convergence, print line and exit if appropriate

       ! ars: check convergence
       call kernel_diis_conv_check(converged, okmessage)

       ! ars: print info line
       call kernel_diis_info_line(converged, okmessage, total_energy_out)

       ! ars: exit if appropiate
       if(converged.or.(iter.eq.pub_kernel_diis_maxit)) then
          ! ars: print banner and exit
          exit kernel_diis_loop
       end if
       !----------------------------------------------------------


       !----------------------------------------------------------
       ! ars: mix Hamiltonians
       if (allow_linear.and.(iter.le.pub_kernel_diis_linear_iter)) then
          ! ars: linear mixing: H_in(i+1) = (1-C)*H_in(i) + C*H_out(i)
          if (do_oda) then
             call kernel_diis_linear_mix(ham_in%ham, ham_out%ham, &
                  dkn_in, dkn_out, ham_in%ham, ham_out%ham, &
                  total_energy_in, total_energy_out)
          else
             call kernel_diis_linear_mix(ham_in%ham, ham_out%ham)
          end if
       else
          ! ars: Pulay mixing: H_in (i+1) = sum_j^i C_j*H_out(j)
          call kernel_diis_pulay_mix(ham_in%ham, ham_his, ham_err_his, &
               rep%overlap, ientry)
       end if

       ! ars: copy dkn_out to dkn_in
       do is = 1, pub_num_spins
          call sparse_embed_copy(dkn_in(is), dkn_out(is))
       end do


       ! ars: shift arrays if we have reached pub_kernel_diis_size iterations
       if (iter.ge.diis_size) then
          shifted = .true.
          call kernel_diis_shift(ham_his, diis_size, pub_num_spins)
          call kernel_diis_shift(ham_err_his, diis_size, pub_num_spins)
       end if


       ! ars: find next position in history to work with
       call kernel_diis_find_ientry(ientry, iter, shifted)
       !----------------------------------------------------------


    end do kernel_diis_loop


    !----


  end subroutine kernel_diis_ham_pulay_scheme

  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine kernel_diis_ham_list_scheme(dkn_in, ham_in, lhxc_fine, &
       ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate)

    !==========================================================================!
    ! Subroutine that performs density kernel mixing and hamiltonian           !
    ! diagonalisation. It minimises the energy with respect to the density     !
    ! kernel. ONETEP inner loop.                                               !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010.                            !
    ! Updated by Alvaro Ruiz Serrano in October 2010 to combine linear+Pulay   !
    ! mixing schemes.                                                          !
    ! Level shifting by Peter Cherry and Alvaro Ruiz Serrano in April 2012.    !
    ! Modified by Alvaro Ruiz Serrano in April 2012 to add ODA and Linear      !
    ! Shooting Techniques (LiSTb and LiSTi).                                   !
    ! Modified for embedding by Joseph Prentice, May 2018                      !
    !==========================================================================!

    use constants, only: DP
#ifdef GPU_PGI
#ifdef GPU_SP_TEST
    use fourier_gpu_wrapper_mod, only: GPU_SP, GPU_SP_switched, GPU_SP_switch, &
                                       deltaE_gpu
#endif
#endif
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_kernel_diis_maxit, pub_kernel_diis_linear_iter, &
         pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy, sparse_embed_trace
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3_EMBED),      intent(inout) :: dkn_in(1:pub_num_spins)
    type(NGWF_HAM),   intent(inout) :: ham_in
    type(MODEL),      intent(in   ) :: mdl
    type(HFX_STATE),  intent(inout), target :: hfxstate
    real(kind=DP),    intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub)
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in   ) :: hub_proj_basis(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(NGWF_REP),   intent(in   ) :: rep

    ! real
    real(kind=DP) :: total_energy_in, total_energy_out
    real(kind=DP), save :: total_energy_old = 0.0_DP
    real(kind=DP), save :: gap_old_1 = 0.0_DP
    real(kind=DP), save :: gap_old_2 = 0.0_DP
    real(kind=DP) :: trace

    ! logical
    logical :: converged, shifted

    ! integer
    integer :: is, ientry

    ! character
    character(len=8) :: okmessage


    ! ars: initilise variables
    converged = .false.
    shifted   = .false.
    ientry  = 1


    ! ars: build initial Ham from rescaled kernel and calculate total energy
    call kernel_diis_build_ham(total_energy_in, ham_in, dkn_in, lhxc_fine, &
         ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
         .true., .false., .true.)


    ! ars: start kernel_diis loop (ONETEP inner loop)
    kernel_diis_loop:    do iter=1,pub_kernel_diis_maxit

#ifdef GPU_SP_TEST
    ! kaw: Accuracy switching toggles
    if ((.not.GPU_SP_switched).and.(abs(deltaE_gpu).gt.(GPU_SP_switch))) then
       GPU_SP=.true.                                                    ! Turn on SP code
    else
       GPU_SP=.false.                                                   ! kaw: Switch to DP
       if (.not.(iter.eq.1)) GPU_SP_switched=.true.
    end if
    if (iter.eq.1) GPU_SP=.true.
#endif

       !----------------------------------------------------------
       ! ars: SCF

       ! ars: calculate total energy with incoming Ham and dkn - do not update ham
       call kernel_diis_build_ham(total_energy_in, ham_in, dkn_in, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
            .true., .false., .false.)


       ! ars: diagonalise Hamiltonian and build dkn (out)
       if (do_ls) then
          ! pjc: Go into level shifting proceedure if appropriate
          call kernel_diis_level_shifter(dkn_out, dkn_in, ham_in%ham, &
               rep%overlap, rep%n_occ)
       else
          ! ars: no level shifter -> diag ham and build kernel
          call kernel_diis_build_idemp_dkn(dkn_out, ham_in%ham, &
               rep%overlap, rep%n_occ, homo, lumo, gap)
       end if


       ! ars: update dens_dep matrices with dkn_out and calculate total energy
       call kernel_diis_build_ham(total_energy_out, ham_out, dkn_out, lhxc_fine, &
            ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, &
            .true., .false., .true.)



       !----------------------------------------------------------
       ! ars: add matrices to history for LIST mixing
       do is = 1, pub_num_spins
          call sparse_embed_copy(ham_his(ientry,is), ham_out%ham(is))
       end do

       ! ars: arrange dkn and ham differentials
       call kernel_diis_diff_inout(dkn_diff_his(ientry,:),dkn_out, dkn_in)
       call kernel_diis_diff_inout(ham_diff_his(ientry,:),ham_out%ham, ham_in%ham)

       ! ars: calculate corrected Hohenberg-Kohn-Sham functional
       do is = 1, pub_num_spins
          call sparse_embed_trace(trace, dkn_diff_his(ientry,is), &
               ham_diff_his(ientry,is))
          total_energy_out = total_energy_out + 0.5_DP * trace
       end do



       !-----------------------------------------------------------
       ! ars: calculate energy, gap, commutator, idemp. and ne.

       ! ars: calculate delta in energy
       deltaE = total_energy_out - total_energy_old
       total_energy_old = total_energy_out

       ! ars: calculate delta in the gap
       delta_gap(1) = gap(1) - gap_old_1
       if (pub_num_spins.eq.2) delta_gap(2) = gap(2) - gap_old_2
       gap_old_1 = gap(1)
       if (pub_num_spins.eq.2) gap_old_2 = gap(2)

       ! ars: calculate residual
       call kernel_diis_residual_inout(residual, ham_out%ham, ham_in%ham, &
            rep%inv_overlap)

       ! ars: calculate commutators and number of electrons
       call kernel_diis_constants(commut, idemp, trks, dkn_out, dkn_in, &
            ham_out%ham, rep%overlap, rep%inv_overlap)


       !----------------------------------------------------------
       ! ars: check convergence, print line and exit if appropriate

       ! ars: check convergence
       call kernel_diis_conv_check(converged, okmessage)

       ! ars: print info line
       call kernel_diis_info_line(converged, okmessage, total_energy_out)

       ! ars: exit if appropiate
       if(converged.or.(iter.eq.pub_kernel_diis_maxit)) then
          ! ars: print banner and exit
          exit kernel_diis_loop
       end if
       !----------------------------------------------------------




       !----------------------------------------------------------
       ! ars: mix Hamiltonians
       if (allow_linear.and.(iter.le.pub_kernel_diis_linear_iter)) then

          ! ars: linear mixing: H_in(i+1) = (1-C)*H_in(i) + C*H_out(i)
          if (do_oda) then
             call kernel_diis_linear_mix(ham_in%ham, ham_out%ham, dkn_in, dkn_out,&
                  ham_in%ham, ham_out%ham, total_energy_in, total_energy_out)
          else
             call kernel_diis_linear_mix(ham_in%ham, ham_out%ham)
          end if

       else

          ! ars: LiST mixing: H_in (i+1) = sum_j^i C_j*H_out(j)
          call kernel_diis_list_mix(ham_in%ham, ham_his, &
               ham_diff_his, dkn_diff_his, ientry)

       end if


       ! ars: copy dkn_out to dkn_in
       do is = 1, pub_num_spins
          call sparse_embed_copy(dkn_in(is), dkn_out(is))
       end do


       ! ars: shift arrays if we have reached pub_kernel_diis_size iterations
       if (iter.ge.diis_size) then
          shifted = .true.
          call kernel_diis_shift(ham_his, diis_size, pub_num_spins)
          call kernel_diis_shift(ham_diff_his, diis_size, pub_num_spins)
          call kernel_diis_shift(dkn_diff_his, diis_size, pub_num_spins)
       end if


       ! ars: find next position in history to work with
       call kernel_diis_find_ientry(ientry, iter, shifted)
       !----------------------------------------------------------


    end do kernel_diis_loop


    !----


  end subroutine kernel_diis_ham_list_scheme





  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !
  !                     END SCF ENERGY MINIMISATION SCHEMES
  !_____________________________________________________________________________
  !_____________________________________________________________________________





  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !
  !                          PRIVATE COMMON ROUTINES
  !_____________________________________________________________________________
  !_____________________________________________________________________________




  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine kernel_diis_check_params

    !========================================================================!
    ! This subroutine checks the correctness of the kernel-DIIS parameters.  !
    !------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in December 2010.                       !
    !========================================================================!

    use comms, only : pub_on_root, comms_barrier, pub_root_proc_id, comms_bcast
    use constants, only: VERBOSE, stdout, DP
    use rundat, only: pub_kernel_diis_size, &
         pub_kernel_diis_maxit, pub_kernel_diis_threshold, &
         pub_kernel_diis_linear_iter, pub_kernel_diis_coeff, &
         pub_kernel_diis_lshift, pub_kernel_diis_ls_iter, &
         pub_kernel_diis_scheme, pub_output_detail, pub_devel_code
    use utils, only: utils_abort, utils_assert

    implicit none


    logical, save :: check_success = .false.
    character(len=200) :: kernel_diis_devel_code
    integer :: start_pos, stop_pos, test_pos


    ! ars: do only the first time that we call this routine
    if (.not.check_success) then

       if (pub_on_root.and.pub_output_detail.ge.VERBOSE) &
            write(stdout,'(/,a35)', advance='no') "Checking kernel-DIIS parameters"


       !-----------------------------------------------------------------------
       ! ars: initialise flags
       do_oda       = .false.
       do_list      = .false.
       do_listi     = .false.
       do_listb     = .false.
       do_pulay     = .false.
       allow_linear = .false.
       mix_dkn      = .false.
       mix_ham      = .false.
       do_diag      = .false.



       !-----------------------------------------------------------------------
       ! ars: check scheme and set flags
       select case (pub_kernel_diis_scheme)

       case ('DIAG')
          scheme = DIAG
          do_diag = .true.
       case ('DKN_LINEAR')
          scheme = DKN_LINEAR
          mix_dkn = .true.
       case ('DKN_PULAY')
          scheme = DKN_PULAY
          mix_dkn = .true.
          do_pulay = .true.
       case ('DKN_LISTI')
          scheme = DKN_LISTI
          mix_dkn = .true.
          do_listi = .true.
       case ('DKN_LISTB')
          scheme = DKN_LISTB
          mix_dkn = .true.
          do_listb = .true.
       case ('HAM_LINEAR')
          scheme = HAM_LINEAR
          mix_ham = .true.
       case ('HAM_PULAY')
          scheme = HAM_PULAY
          mix_ham = .true.
          do_pulay = .true.
       case ('HAM_LISTI')
          scheme = HAM_LISTI
          mix_ham = .true.
          do_listi = .true.
       case ('HAM_LISTB')
          scheme = HAM_LISTB
          mix_ham = .true.
          do_listb = .true.
       case default
          call utils_abort("Kernel DIIS scheme not recognised.")
       end select

       ! ars: set generalised schemes if appropiate
       do_list = (do_listi.or.do_listb)
       if (mix_dkn.and.do_list)  scheme = DKN_LIST
       if (mix_ham.and.do_list)  scheme = HAM_LIST


       !-----------------------------------------------------------------------
       ! ars: check pub_kernel_diis_maxit and abort if negative
       call utils_assert(pub_kernel_diis_maxit > 0, &
            'Invalid kernel-DIIS max iterations')

       !-----------------------------------------------------------------------
       ! ars: check pub_kernel_diis_size and correct if negative
       if (pub_kernel_diis_size.le.0) then

          if(pub_on_root) then
             write(stdout,'(/,a)') &
                  "*** WARNING: Invalid kernel-DIIS max ***************"
             write(stdout,'(a)')    "             kernel_diis_size <= 0"
             write(stdout,'(a,i3)') &
                  "             kernel_diis_size = ", pub_kernel_diis_size
             write(stdout,'(a,/)') &
                  "             Reseting to safe value kernel_diis_size = 5"
          end if

          pub_kernel_diis_size = 5

       end if

       !-----------------------------------------------------------------------
       ! ars: check pub_kernel_diis_size is not greater than pub_kernel_diis_maxit
       if (pub_kernel_diis_size.gt.pub_kernel_diis_maxit) then

          if(pub_on_root) then
             write(stdout,'(/,a)') &
                  "*** WARNING: Invalid kernel-DIIS max ***************"
             write(stdout,'(a)')    "             kernel_diis_size > kernel_diis_maxit"
             write(stdout,'(a,i3)')&
                  "             kernel_diis_size = ", pub_kernel_diis_size
          end if

          if (pub_kernel_diis_maxit.le.5) then
             pub_kernel_diis_size = pub_kernel_diis_maxit
          else
             pub_kernel_diis_size = 5
          end if

          if(pub_on_root) write(stdout,'(a,i1,/)') &
               "             Reseting to safe value kernel_diis_size = ", &
               pub_kernel_diis_size
       end if

       !-----------------------------------------------------------------------
       ! ars: check pub_kernel_diis_linear_iter is not greater than pub_kernel_diis_maxit
       if (pub_kernel_diis_linear_iter.gt.pub_kernel_diis_maxit) then

          if(pub_on_root) then
             write(stdout,'(/,a)')  &
                  "*** WARNING: Invalid kernel-DIIS linear iterations ******"
             write(stdout,'(a,i3)') &
                  "             kernel_diis_liter = ", pub_kernel_diis_linear_iter
          end if

          if (pub_kernel_diis_maxit.le.3) then
             pub_kernel_diis_linear_iter = pub_kernel_diis_maxit
          else
             pub_kernel_diis_linear_iter = 3
          end if


          if(pub_on_root) write(stdout,'(a,i1,/)') &
               "             Reseting to safe value kernel_diis_liter = ",&
               pub_kernel_diis_linear_iter
       end if

       ! ars: set flag if we allow linear iter when Pulay or LiST
       if (pub_kernel_diis_linear_iter.gt.0) allow_linear = .true.




       !-----------------------------------------------------------------------
       ! ars: check C_out for linear mixing and set ODA flag
       if (pub_kernel_diis_coeff.gt.1.0_DP) then

          if(pub_on_root) then
             write(stdout,'(/,a)')     &
                  "*** WARNING: Invalid kernel-DIIS C_out coefficient ******"
             write(stdout,'(a,f22.12)')&
                  "             kernel_diis_coeff = ", pub_kernel_diis_coeff
             write(stdout,'(a,/)')     &
                  "             Negative number for ODA, or else in range [0,1]"
             write(stdout,'(a,/)')     &
                  "             Reseting to safe value kernel_diis_coeff=0.1000"
          end if

          pub_kernel_diis_coeff = 0.1_DP

       else if (pub_kernel_diis_coeff.lt.0.0_DP) then

          ! ars: activate ODA
          do_oda = .true.

       end if


       !-----------------------------------------------------------------------
       ! ars: check DIIS threshold
       if(pub_kernel_diis_threshold.le.0.0_DP) then
          if(pub_on_root) then
             write(stdout,'(/,a)')     &
                  "*** WARNING: Invalid kernel-DIIS theshold **************"
             write(stdout,'(a,f22.12)')&
                  "             kernel_diis_threshold = ", pub_kernel_diis_threshold
             write(stdout,'(a,/)')     &
                  "             Reseting to safe value kernel_diis_threshold=1E-09"
          end if

          pub_kernel_diis_threshold = 0.000000001_DP

       end if


       !-----------------------------------------------------------------------
       ! ars: set history size
       if (do_pulay.or.do_list) then
          diis_size = pub_kernel_diis_size
       else
          diis_size = 1
       end if


       !-----------------------------------------------------------------------
       ! ars: exit nice and happy
       check_success = .true.
       if (pub_on_root.and.pub_output_detail.ge.VERBOSE) &
            write(stdout,'(a,/)') " ... done"


       call comms_barrier


       !-----------------------------------------------------------------------
       !-----------------------------------------------------------------------
       !-----------------------------------------------------------------------
       ! ars: flags due to pub_devel_code
       !-----------------------------------------------------------------------
       ! ars: syntax: devel_code: KERNEL_DIIS:LS_DKN_EIGS=T:KERNEL_DIIS
       !-----------------------------------------------------------------------

       ls_dkn_eigs = .false.
       if (pub_on_root) then
          kernel_diis_devel_code=pub_devel_code
          if (len_trim(kernel_diis_devel_code)>0) then
             start_pos=index(kernel_diis_devel_code,'KERNEL_DIIS:')
             stop_pos=index(kernel_diis_devel_code,':KERNEL_DIIS')
             if (stop_pos<=0) stop_pos=len_trim(kernel_diis_devel_code) !missing end so go to end of string
             if (start_pos>0) then

                ! ars: set finite differences scheme
                test_pos=index(kernel_diis_devel_code,'LS_DKN_EIGS=')
                if (test_pos>start_pos.and.test_pos<stop_pos) then
                   test_pos=test_pos+len('LS_DKN_EIGS=')
                   read(kernel_diis_devel_code(test_pos:test_pos+ &
                        & index(kernel_diis_devel_code(test_pos:stop_pos),':')-2),*) ls_dkn_eigs
                   ! ars: write warning on screen
                   write(stdout,'(a,l1)') "kernel_diis warning: devel_code: ls_dkn_eigs = ", ls_dkn_eigs

                end if
             end if
          end if
       end if
       call comms_bcast(pub_root_proc_id,ls_dkn_eigs)

    end if


    !-----------------------------------------------------------------------
    ! ars: set level shifter
    do_ls = .false.
    if ( (pub_kernel_diis_lshift.gt.0.0_DP).and.&
         (pub_kernel_diis_ls_iter.gt.0) ) do_ls = .true.


  end subroutine kernel_diis_check_params



  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine kernel_diis_alloc(rep,mdl)

    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_HAM, NGWF_REP, &
         ngwf_ham_create
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create
    use utils, only: utils_alloc_check

    implicit none

    type(NGWF_REP), intent(in   ) :: rep
    type(MODEL),          intent(in   ) :: mdl

    integer :: is, ierr, icntr
    ! agrecocmplx
    logical :: loc_cmplx

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(1)%iscmplx

    ! ars: ---- ALLOCATE MATRICES--------------------
    ! ars: allocate dkn_out - always necessary
    allocate(dkn_out(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('kernel_diis_alloc', 'dkn_out', ierr)
    do is = 1, pub_num_spins
       dkn_out(is)%structure='K'
       ! agrecocmplx
       call sparse_embed_create(dkn_out(is), iscmplx=loc_cmplx)
    end do


    ! ars: allocate ham_out (unless we are doing diagonalisation only)
    if (.not.do_diag) call ngwf_ham_create(ham_out,rep)


    ! ars: allocate history of dkn to mix if required for dkn Pulay or LIST
    if (mix_dkn.and.(do_pulay.or.do_list)) then
       allocate(dkn_his(1:diis_size, 1:pub_num_spins), stat=ierr)
       call utils_alloc_check('kernel_diis_alloc', 'dkn_his', ierr)
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             dkn_his(icntr, is)%structure = 'K'
             ! agrecocmplx
             call sparse_embed_create(dkn_his(icntr, is), iscmplx=loc_cmplx)
          end do
       end do
    end if


    ! ars: allocate history of ham to mix if required for ham Pulay or LIST
    if (mix_ham.and.(do_pulay.or.do_list)) then
       allocate(ham_his(1:diis_size, 1:pub_num_spins), stat=ierr)
       call utils_alloc_check('kernel_diis_alloc', 'ham_his', ierr)
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             ham_his(icntr, is)%structure = 'H'
             ! agrecocmplx
             call sparse_embed_create(ham_his(icntr, is), iscmplx=loc_cmplx)
          end do
       end do
    end if

    ! ars: allocate history of dkn error matrices if required for dkn-Pulay
    if (do_pulay.and.mix_dkn) then
       allocate(dkn_err_his(1:diis_size, 1:pub_num_spins), stat=ierr)
       call utils_alloc_check('kernel_diis_alloc', 'dkn_err_his', ierr)
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             dkn_err_his(icntr, is)%structure = 'K'
             ! agrecocmplx
             call sparse_embed_create(dkn_err_his(icntr, is), iscmplx=loc_cmplx)
          end do
       end do
    end if


    ! ars: allocate history of ham error matrices if required for ham-Pulay
    if (do_pulay.and.mix_ham) then
       allocate(ham_err_his(1:diis_size, 1:pub_num_spins), stat=ierr)
       call utils_alloc_check('kernel_diis_alloc', 'ham_err_his', ierr)
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             ham_err_his(icntr, is)%structure = 'H'
             ! agrecocmplx
             call sparse_embed_create(ham_err_his(icntr, is), iscmplx=loc_cmplx)
          end do
       end do
    end if


    ! ars: allocate history of dkn and ham differentials if required for LIST
    if (do_list.and.(mix_dkn.or.mix_ham)) then
       allocate(dkn_diff_his(1:diis_size, 1:pub_num_spins), stat=ierr)
       call utils_alloc_check('kernel_diis_alloc', 'dkn_diff_his', ierr)
       allocate(ham_diff_his(1:diis_size, 1:pub_num_spins), stat=ierr)
       call utils_alloc_check('kernel_diis_alloc', 'ham_diff_his', ierr)
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             dkn_diff_his(icntr, is)%structure = 'K'
             ! agrecocmplx
             call sparse_embed_create(dkn_diff_his(icntr, is), iscmplx=loc_cmplx)
             ham_diff_his(icntr, is)%structure = 'H'
             ! agrecocmplx
             call sparse_embed_create(ham_diff_his(icntr, is), iscmplx=loc_cmplx)
          end do
       end do
    end if


    ! ars: ---- ALLOCATE QUANTITIES --------------
    allocate(commut(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('kernel_diis_alloc', 'commut', ierr)
    allocate(idemp(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('kernel_diis_alloc', 'idemp', ierr)
    allocate(trks(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('kernel_diis_alloc', 'trks', ierr)
    allocate(residual(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('kernel_diis_alloc', 'residual', ierr)
    allocate(homo(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('kernel_diis_alloc', 'homo', ierr)
    allocate(lumo(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('kernel_diis_alloc', 'lumo', ierr)
    allocate(gap(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('kernel_diis_alloc', 'gap', ierr)
    allocate(delta_gap(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('kernel_diis_alloc', 'delta_gap', ierr)


  end subroutine kernel_diis_alloc

  !_____________________________________________________________________________
  !_____________________________________________________________________________




  subroutine kernel_diis_dealloc()

    use ngwf_representation, only: NGWF_HAM, ngwf_ham_destroy
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_destroy
    use utils, only: utils_dealloc_check

    implicit none

    integer :: is, ierr, icntr


    ! ars: ---- DEALLOCATE QUANTITIES --------------
    deallocate(commut, stat=ierr)
    call utils_dealloc_check('kernel_diis_dealloc', 'commut', ierr)
    deallocate(idemp, stat=ierr)
    call utils_dealloc_check('kernel_diis_dealloc', 'idemp', ierr)
    deallocate(trks, stat=ierr)
    call utils_dealloc_check('kernel_diis_dealloc', 'trks', ierr)
    deallocate(residual, stat=ierr)
    call utils_dealloc_check('kernel_diis_dealloc', 'residual', ierr)
    deallocate(homo, stat=ierr)
    call utils_dealloc_check('kernel_diis_dealloc', 'homo', ierr)
    deallocate(lumo, stat=ierr)
    call utils_dealloc_check('kernel_diis_dealloc', 'lumo', ierr)
    deallocate(gap, stat=ierr)
    call utils_dealloc_check('kernel_diis_dealloc', 'gap', ierr)
    deallocate(delta_gap, stat=ierr)
    call utils_dealloc_check('kernel_diis_dealloc', 'delta_gap', ierr)

    ! ars: ---- DEALLOCATE MATRICES --------------
    ! ars: deallocate history of ham differentials if required for ham-LIST
    if (allocated(ham_diff_his)) then
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             call sparse_embed_destroy(ham_diff_his(icntr, is))
          end do
       end do
       deallocate(ham_diff_his, stat=ierr)
       call utils_dealloc_check('kernel_diis_dealloc', 'ham_diff_his', ierr)
    end if


    ! ars: deallocate history of dkn differentials if required for dkn-LIST
    if (allocated(dkn_diff_his)) then
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             call sparse_embed_destroy(dkn_diff_his(icntr, is))
          end do
       end do
       deallocate(dkn_diff_his, stat=ierr)
       call utils_dealloc_check('kernel_diis_dealloc', 'dkn_diff_his', ierr)
    end if

    ! ars: deallocate history of ham error matrices if required for ham-Pulay
    if (allocated(ham_err_his)) then
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             call sparse_embed_destroy(ham_err_his(icntr, is))
          end do
       end do
       deallocate(ham_err_his, stat=ierr)
       call utils_dealloc_check('kernel_diis_dealloc', 'ham_err_his', ierr)
    end if


    ! ars: deallocate history of dkn error matrices if required for dkn-Pulay
    if (allocated(dkn_err_his)) then
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             call sparse_embed_destroy(dkn_err_his(icntr, is))
          end do
       end do
       deallocate(dkn_err_his, stat=ierr)
       call utils_dealloc_check('kernel_diis_dealloc', 'dkn_err_his', ierr)
    end if

    ! ars: deallocate history of ham to mix if required for ham Pulay or LIST
    if (allocated(ham_his)) then
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             call sparse_embed_destroy(ham_his(icntr, is))
          end do
       end do
       deallocate(ham_his, stat=ierr)
       call utils_dealloc_check('kernel_diis_dealloc', 'ham_his', ierr)
    end if


    ! ars: deallocate history of dkn to mix if required for dkn Pulay or LIST
    if (allocated(dkn_his)) then
       do is = 1, pub_num_spins
          do icntr = 1, diis_size
             call sparse_embed_destroy(dkn_his(icntr, is))
          end do
       end do
       deallocate(dkn_his, stat=ierr)
       call utils_dealloc_check('kernel_diis_dealloc', 'dkn_his', ierr)
    end if


    ! ars: deallocate ham_out (unless we are doing diagonalisation only)
    if (.not.do_diag) call ngwf_ham_destroy(ham_out)


    ! ars: deallocate dkn_out - always necessary
    do is = 1, pub_num_spins
       call sparse_embed_destroy(dkn_out(is))
    end do
    deallocate(dkn_out, stat=ierr)
    call utils_dealloc_check('kernel_diis_dealloc', 'dkn_out', ierr)


  end subroutine kernel_diis_dealloc

  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine kernel_diis_ls_alloc(iscmplx)

    use dense, only: dense_create
    use rundat, only: pub_num_spins
    use utils, only: utils_alloc_check

    implicit none

    ! Argument
    ! agrecocmplx
    logical, optional, intent(in) :: iscmplx

    integer :: is, ierr
    ! agrecocmplx
    logical :: loc_cmplx

    ! agrecocmplx
    if (present(iscmplx)) then
       loc_cmplx = iscmplx
    else
       loc_cmplx = .false.
    end if

    if (.not.allocated(eigvecs)) then
       allocate(eigvecs(1:pub_num_spins), stat=ierr)
       call utils_alloc_check('kernel_diis_ls_alloc', 'eigvecs', ierr)
       do is = 1, pub_num_spins
          ! agrecocmplx
          call dense_create(eigvecs(is), num, num, iscmplx=loc_cmplx)
       end do
    end if

  end subroutine kernel_diis_ls_alloc

  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine kernel_diis_ls_dealloc()

    use dense, only: dense_destroy
    use rundat, only: pub_num_spins
    use utils, only: utils_dealloc_check

    implicit none

    integer :: is, ierr

    if (allocated(eigvecs)) then

       do is = 1, pub_num_spins
          call dense_destroy(eigvecs(is))
       end do
       deallocate(eigvecs, stat=ierr)
       call utils_dealloc_check('kernel_diis_ls_dealloc', 'eigvecs', ierr)

    end if

  end subroutine kernel_diis_ls_dealloc


  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_build_ham(total_energy, ham, denskern, lhxc_fine, &
       ngwf_basis, hub_proj_basis, hub, rep, mdl, hfxstate, ham_update, &
       lhxc_fine_fixed, build_ham)

    !====================================================================!
    ! This subroutine calculates the density dependent matrices and      !
    ! builds the Hamiltonian matrix in NGWF representation.              !
    !--------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in June 2010.                       !
    ! Modified by Alvaro Ruiz Serrano in May 2012 to return total_energy.!
    ! Modified for embeddding by Joseph Prentice, May 2018               !
    !====================================================================!

    use constants, only: DP, paw_en_size
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_dens_dep_matrices,&
         hamiltonian_build_matrix
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_num_spins, pub_num_kpoints, PUB_1K, &
         pub_dmft_spoil_kernel, &
         pub_xc_ke_density_required
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_trace, &
         sparse_embed_array_create, sparse_embed_array_copy, &  !!! Needed for workaround
         sparse_embed_array_destroy, sparse_embed_copy !!! Needed for workaround
    use utils, only: utils_assert

    implicit none

    ! ars: arguments
    real(kind=DP),    intent(  out) :: total_energy
    type(NGWF_HAM),   intent(inout) :: ham
    type(SPAM3_EMBED),      intent(inout) :: denskern(pub_num_spins)
    type(MODEL),      intent(in   ) :: mdl
    type(HFX_STATE),  intent(inout), target :: hfxstate
    real(kind=DP),    intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub)
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in   ) :: hub_proj_basis(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(NGWF_REP),   intent(in   ) :: rep
    logical,          intent(in   ) :: ham_update
    logical,          intent(in   ) :: lhxc_fine_fixed
    logical,          intent(in   ) :: build_ham

    ! ars: local variables
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    real(kind=DP) :: lhxc_energy
    real(kind=DP) :: hubbard_energy
    logical, save :: updated=.false.

!!!BEGIN WORKAROUND (PART 1)
    !jme: Temporary workaround (to be removed when hamiltonians
    !     become SPARSE_ARRAYs).
    type(SPAM3_EMBED_ARRAY) :: aux_denskern
    integer :: is

    ! JCW: Abort if KE density required (meta-GGA + DIIS not implemented)
    call utils_assert(.not.pub_xc_ke_density_required, "Error in &
         &kernel_diis_build_ham: pub_xc_ke_density_required is true, but &
         &combination of KE-density-dependent XC functionals with kernel &
         &DIIS has not been implemented/tested")

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine kernel_diis_build_ham not ready yet for more&
         & than one k-point.')

    !jme: create SPARSE_ARRAY
    call sparse_embed_array_create(aux_denskern, denskern(1), &
         n_spins=size(denskern), n_kpoints=PUB_1K)
    !jme: copy contents of denskern into it
    do is = 1, pub_num_spins
       call sparse_embed_copy(aux_denskern%m(is, PUB_1K), denskern(is))
    end do
!!!END WORKAROUND (PART 1)

    ! ars: build hamiltonian:
    !   1: dens_indep_matrices already initialised
    !   2: initiliase density independent matrices
    call hamiltonian_dens_dep_matrices(ham, lhxc_fine, total_energy, &
         lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
         ngwf_basis, hub_proj_basis, hub, aux_denskern, mdl, hfxstate, &
         ham_update, lhxc_fine_fixed)

!!!BEGIN WORKAROUND (PART 2)
    !jme: copy contents of denskern into it
    do is = 1, pub_num_spins
       call sparse_embed_copy(denskern(is), aux_denskern%m(is, PUB_1K))
    end do

!!!END WORKAROUND (PART 2)

    !   3: build hamiltonian matrix
    if (build_ham) call hamiltonian_build_matrix(ham, rep)
    ! ebl note: may need to update arguments of hamiltonian_dens_dep_matrices
    !CW
    if(pub_dmft_spoil_kernel.and..not.updated)then
       updated=.true.
       call hamiltonian_dens_dep_matrices(ham, lhxc_fine, total_energy, &
            lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
            ngwf_basis, hub_proj_basis, hub, aux_denskern, mdl, hfxstate, &
            ham_update, lhxc_fine_fixed,spoil_force=.true.)
       call hamiltonian_build_matrix(ham, rep)
    endif
    !END CW

    call sparse_embed_array_destroy(aux_denskern)


  end subroutine kernel_diis_build_ham


  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_diag_ham(eigs_dens, ham, overlap, n_occ, num, &
       Ehomo, Elumo, Egap)

    !=====================================================================!
    ! This subroutine diagonalises the hamiltonian and builds the output  !
    ! density from the eigenvectors.                                      !
    !---------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010 as an adapation of the !
    ! existing routine kernel_reset_occupancies.                          !
    ! Modified for embedding by Joseph Prentice, May 2018                 !
    !=====================================================================!

    use constants, only: DP, max_spins
    use dense, only: DEM, dense_create, dense_convert, dense_destroy, &
         dense_eigensolve
    use rundat, only: pub_eigensolver_orfac, pub_eigensolver_abstol, &
         pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_trace
    use utils, only : utils_alloc_check, utils_dealloc_check
    use timer, only: timer_clock


    implicit none

    ! ars: arguments
    type(DEM),     intent(inout) :: eigs_dens(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: ham(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: overlap
    integer,       intent(in   ) :: n_occ(max_spins)
    integer,       intent(in   ) :: num
    real(kind=DP), optional, intent(  out) :: Ehomo(pub_num_spins)
    real(kind=DP), optional, intent(  out) :: Elumo(pub_num_spins)
    real(kind=DP), optional, intent(  out) :: Egap(pub_num_spins)

    ! ars: local variables
    integer   :: ierr
    integer   :: is
    real(kind=DP), allocatable :: eigenvalues(:)
    type(DEM) :: overlap_dens, ham_dens
    logical :: do_gap
    ! agrecocmplx
    logical :: loc_cmplx


    ! ars: start timer
    call timer_clock('kernel_diis_diag_ham',1)

    ! ars: check optional arguments
    do_gap=.false.
    if (present(Ehomo).and.present(Elumo).and.present(Egap)) do_gap=.true.

    ! ars: allocate Hamiltonian eigenvalues
    allocate(eigenvalues(num),stat=ierr)
    call utils_alloc_check('kernel_diis_diag_ham','eigenvalues',ierr)

    ! agrecocmplx
    loc_cmplx = overlap%p%iscmplx

    ! ars: DEM buffers
    ! agrecocmplx
    call dense_create(ham_dens,num,num,iscmplx=loc_cmplx) !ars
    call dense_create(overlap_dens,num,num,iscmplx=loc_cmplx)


    ! cks: loop over spins
    do is=1, pub_num_spins

       ! ars: init
       eigenvalues(:) = 0.0_DP

       ! ars: expand denskern, ham and overlap to dense matrices
       call dense_convert(ham_dens,ham(is))
       call dense_convert(overlap_dens,overlap)

       ! ars: solve HC = eSC
       call dense_eigensolve(num,eigenvalues,ham_dens,overlap_dens, &
            1, eigs_dens(is),pub_eigensolver_orfac, pub_eigensolver_abstol)

       ! ars: calculate HOMO, LUMO and energy gap
       if (do_gap) then
          Ehomo(is) = eigenvalues(n_occ(is))
          Elumo(is) = eigenvalues(n_occ(is)+1)
          Egap(is) = Elumo(is) - Ehomo(is)
       end if

    end do

    ! ars: destroy DEM buffers
    call dense_destroy(overlap_dens)
    call dense_destroy(ham_dens) !ars

    ! ars: deallocate
    deallocate(eigenvalues,stat=ierr)
    call utils_dealloc_check('kernel_diis_diag_ham','eigenvalues',ierr)

    ! ars: stop timer
    call timer_clock('kernel_diis_diag_ham',2)


  end subroutine kernel_diis_diag_ham



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_dkn_from_eigs(denskern, eigs_dens, n_occ, num)

    !=====================================================================!
    ! This subroutine diagonalises the hamiltonian and builds the output  !
    ! density from the eigenvectors.                                      !
    !---------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010 as an adapation of the !
    ! existing routine kernel_reset_occupancies.                          !
    ! Modified for embedding by Joseph Prentice, May 2018                 !
    !=====================================================================!

    use constants, only: max_spins
    use dense, only: dense_create, dense_convert, dense_destroy, &
         dense_product
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED
    use timer, only: timer_clock


    implicit none

    ! ars: arguments
    type(SPAM3_EMBED), intent(inout) :: denskern(pub_num_spins)
    type(DEM),   intent(in   ) :: eigs_dens(pub_num_spins)
    integer,     intent(in   ) :: n_occ(max_spins)
    integer,     intent(in   ) :: num

    ! ars: local variables
    integer   :: is
    type(DEM) :: denskern_dens
    ! agrecocmplx
    logical :: loc_cmplx
    character :: opB_loc

    ! ars: start timer
    call timer_clock('kernel_diis_dkn_from_eigs',1)

    ! agrecocmplx
    loc_cmplx = denskern(1)%p%iscmplx

    ! agrecocmplx
    if (loc_cmplx) then
       opB_loc = 'C'
    else
       opB_loc = 'T'
    end if

    ! ars: DEM buffer
    ! agrecocmplx
    call dense_create(denskern_dens,num,num,iscmplx=loc_cmplx)

    ! cks: loop over spins
    do is=1, pub_num_spins

       ! ndmh: construct a density kernel from the first n_occ(is) eigenvectors
       ! agrecocmplx: need to take hermitian of matrix B in complex case?
       ! I think so....
       call dense_product(denskern_dens, eigs_dens(is), eigs_dens(is), &
            opA='N',opB=opB_loc,first_k=1,last_k=n_occ(is))

       ! ndmh: convert back to sparse matrix
       call dense_convert(denskern(is),denskern_dens)

    end do

    ! ars: destroy DEM buffer
    call dense_destroy(denskern_dens)

    ! ars: stop timer
    call timer_clock('kernel_diis_dkn_from_eigs',2)


  end subroutine kernel_diis_dkn_from_eigs



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_mu(mu, denskern, ham, overlap)

    !=====================================================================!
    ! This subroutine calculates mu for kernel-rescaling.                 !
    !---------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in January 2011.                     !
    ! Modified for embedding by Joseph Prentice, May 2018                 !
    !=====================================================================!

    use constants, only: max_spins, DP
    use rundat, only: pub_num_spins, PUB_1K
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_trace

    implicit none

    real(kind=DP), intent(  out) :: mu(max_spins)
    type(SPAM3_EMBED_ARRAY),   intent(in   ) :: denskern
    type(SPAM3_EMBED),   intent(in   ) :: ham(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: overlap

    integer :: is
    real(kind=DP) :: trace

    ! ars: calculate mu
    do is = 1, pub_num_spins
       call sparse_embed_trace(mu(is), denskern%m(is,PUB_1K), ham(is))
       call sparse_embed_trace(trace, denskern%m(is,PUB_1K), overlap)
       mu(is) = mu(is) / trace
    end do

  end subroutine kernel_diis_mu


  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_constants(hks_skh, k_ksk, ks_ne, dkn_out, dkn_in, &
       ham, overlap, inv_overlap)

    !==================================================================!
    ! This subroutine calculates the constants of the calculation      !
    ! related to the density kernel and prints a summary:              !
    !  1) [HKS,SKH]                                                    !
    !  2) [K,KSK]                                                      !
    !  3) Ne = 2tr[KS]                                                 !
    !------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2010.                     !
    ! Modified by D.D. O'Regan for commutator invariance in July 2013. !
    ! Modified for embedding by Joseph Prentice, May 2018              !
    !==================================================================!

    use constants, only: DP, SAFE_DIV_EPS
    use rundat, only: pub_rms_kernel_measure, pub_num_spins, pub_spin_fac
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_rms_element, sparse_embed_trace, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_transpose

    implicit none

    ! ars: arguments
    real(kind=DP),    intent(  out) :: hks_skh(pub_num_spins)
    real(kind=DP),    intent(  out) :: k_ksk(pub_num_spins)
    real(kind=DP),    intent(  out) :: ks_ne(pub_num_spins)
    type(SPAM3_EMBED),      intent(in   ) :: dkn_out(pub_num_spins)
    type(SPAM3_EMBED),      intent(in   ) :: dkn_in(pub_num_spins)
    type(SPAM3_EMBED),      intent(in   ) :: ham(pub_num_spins)
    type(SPAM3_EMBED),      intent(in   ) :: overlap
    type(SPAM3_EMBED),      intent(in   ) :: inv_overlap !ddor

    ! ars: local variables
    integer :: is
    type(SPAM3_EMBED) :: sk, hk, skh, hks, ksk, ksks, mixed_commutator


    ! ars: create matrices
    ! ddor: Changed to HKH structure
    ! ndmh: (a destroyed matrix is fine for structure code use)
    ! jcap: this is not true here for embedding matrices, since
    ! we need the number of regions to allocate the new matrix
    call sparse_embed_create(sk, overlap, dkn_out(1))
    call sparse_embed_create(ksk, dkn_out(1), sk)
    call sparse_embed_create(hk,ham(1),dkn_out(1))
    call sparse_embed_create(skh, hk, ham(1))
    call sparse_embed_destroy(hk)
    call sparse_embed_create(hks, skh)
    if (.not. pub_rms_kernel_measure) then ! ddor
       call sparse_embed_create(mixed_commutator,hks,inv_overlap)
       call sparse_embed_create(ksks,ksk,overlap)
    endif

    do is = 1, pub_num_spins

       ! ars: calculate matrices
       call sparse_embed_product(sk, overlap, dkn_out(is))
       call sparse_embed_product(skh, sk, ham(is))
       call sparse_embed_transpose(hks, skh)

       ! ddor: calculate tr[KS]
       call sparse_embed_trace(ks_ne(is),sk)

       ! ars: calculate [HKS,SKH]
       call sparse_embed_axpy(hks, skh, -1.0_DP)

       ! ddor: Changed for invariant commutator
       if (pub_rms_kernel_measure) then
          hks_skh(is) = sparse_embed_rms_element(hks)
          hks_skh(is) = abs(hks_skh(is))
       else
          ! ddor: Calculate the mixed-index commutator C = (SKH-HKS)
          call sparse_embed_product(mixed_commutator,hks,inv_overlap)
          ! ddor: Trace the squared mixed-index commutator and renormalise
          !     : to give -Tr[C^2] / N^2
          ! ddor: Allow for case of empty density kernel for one spin
          if (abs(ks_ne(is)) < SAFE_DIV_EPS) then
             call sparse_embed_trace(hks_skh(is),mixed_commutator,mixed_commutator)
             hks_skh(is) = -1.0_DP * hks_skh(is) / (ks_ne(is))**2.0_DP
             ! ddor: This is already a per-electron measure, just take the norm
             hks_skh(is) = sqrt(hks_skh(is))
          else
             hks_skh(is) = 0.0_DP
          endif
       endif

       ! ars: calculate [K,KSK]
       call sparse_embed_product(sk, overlap, dkn_in(is))
       call sparse_embed_product(ksk, dkn_in(is), sk)
       call sparse_embed_axpy(ksk, dkn_in(is), -1.0_DP)
       if (pub_rms_kernel_measure) then ! ddor
          k_ksk(is) = sparse_embed_rms_element(ksk)
          k_ksk(is) = abs(k_ksk(is))
       else
          call sparse_embed_product(ksks,ksk,overlap)
          ! ddor: Trace the squared mixed-index idenpotency error
          !     : and renormalise to give -Tr[(rho^2-rho)^2] / N^2
          ! ddor: Allow for case of empty density kernel  for one spin
          if (abs(ks_ne(is)) < SAFE_DIV_EPS) then
             call sparse_embed_trace(k_ksk(is),ksks,ksks)
             k_ksk(is) = k_ksk(is) / (ks_ne(is))**2.0_DP
             ! ddor: This is already a per-electron measure, just take the norm
             k_ksk(is) = sqrt(k_ksk(is))
          else
             k_ksk(is) = 0.0_DP
          endif
       endif

    end do

    ! ars: destroy matrices
    if (.not. pub_rms_kernel_measure) then ! ddor
       call sparse_embed_destroy(ksks)
       call sparse_embed_destroy(mixed_commutator)
    endif
    call sparse_embed_destroy(hks)
    call sparse_embed_destroy(skh)
    call sparse_embed_destroy(ksk)
    call sparse_embed_destroy(sk)

    ! ars: calculate 2tr[KS]
    do is = 1, pub_num_spins
       ks_ne(is) = pub_spin_fac * ks_ne(is)
    end do


  end subroutine kernel_diis_constants


  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_residual_inout(res_inout, mat_out, mat_in, overlap)

    !=====================================================================!
    ! This subroutine checks for convergence of the density kernel.       !
    ! Based on the latest residue matrix only.                            !
    !---------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010.                       !
    ! Modified by Alvaro Ruiz Serrano in May 2012.                        !
    ! Modified for embedding by Joseph Prentice, May 2018                 !
    !=====================================================================!

    use constants, only: DP
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_trace, &
         sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_transpose, sparse_embed_copy, &
         sparse_embed_axpy

    implicit none

    real(kind=DP), intent(  out) :: res_inout(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: mat_out(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: mat_in(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: overlap

    integer :: is
    type(SPAM3_EMBED) :: res_ov, res_t, res_ov_res_t, res



    residual(:) = 0.0_DP

    ! ars: init spam3buffer blocking scheme
    call sparse_embed_create(res, mat_out(1))
    call sparse_embed_create(res_ov, res, overlap)
    call sparse_embed_create(res_t, res)
    call sparse_embed_create(res_ov_res_t, res_ov, res_t)


    do is = 1, pub_num_spins

       ! ars: calculate residue
       call sparse_embed_copy(res, mat_out(is))
       call sparse_embed_axpy(res, mat_in(is), -1.0_DP)

       ! ars: res_ov = R x S
       call sparse_embed_product(res_ov, res, overlap)

       ! ars: res_t = (R)^t
       call sparse_embed_transpose(res_t,res)

       ! ars: res_ov_res_t = (R) x S x (R)^t
       call sparse_embed_product(res_ov_res_t, res_ov, res_t)

       ! ars: Res(is) = + sqrt( tr[(R) x S x (R)^t x S] )
       call sparse_embed_trace(res_inout(is),res_ov_res_t,overlap)
       res_inout(is) = abs(sqrt(res_inout(is)))

    end do

    ! ars: destroy SPAM3 buffers
    call sparse_embed_destroy(res_ov)
    call sparse_embed_destroy(res_t)
    call sparse_embed_destroy(res_ov_res_t)
    call sparse_embed_destroy(res)


  end subroutine kernel_diis_residual_inout



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_conv_check(converged, okmessage)

    use rundat, only: pub_kernel_diis_conv_criteria, pub_kernel_diis_threshold, &
         pub_num_spins

    implicit none

    ! ars: arguments
    logical,          intent(  out) :: converged
    character(len=8), intent(  out) :: okmessage

    ! ars: local variables
    integer :: is
    logical :: conv_criteria(1:4)

    converged = .false.
    conv_criteria(1:4) = .false.
    okmessage(1:8) = '        '


    ! ars: validate convergence criteria:

    ! residual
    if (pub_kernel_diis_conv_criteria(1:1).eq.'1') then
       residual_loop: do is = 1, pub_num_spins
          if (abs(residual(is)).lt.pub_kernel_diis_threshold) then
             okmessage(1:2) = 'ok'
             conv_criteria(1) = .true.
          else
             okmessage(1:2) = 'no'
             conv_criteria(1) = .false.
             exit residual_loop
          end if
       end do residual_loop
    else
       conv_criteria(1) = .true.
    end if

    ! [HKS,SKH]
    if (pub_kernel_diis_conv_criteria(2:2).eq.'1') then
       hks_loop: do is = 1, pub_num_spins
          if (abs(commut(is)).lt.pub_kernel_diis_threshold) then
             okmessage(3:4) = 'ok'
             conv_criteria(2) = .true.
          else
             okmessage(3:4) = 'no'
             conv_criteria(2) = .false.
             exit hks_loop
          end if
       end do hks_loop
    else
       conv_criteria(2) = .true.
    end if

    ! Gap energy variation
    if (pub_kernel_diis_conv_criteria(3:3).eq.'1') then
       gap_loop: do is = 1, pub_num_spins
          if (abs(delta_gap(is)).lt.pub_kernel_diis_threshold) then
             okmessage(5:6) = 'ok'
             conv_criteria(3) = .true.
          else
             conv_criteria(3) = .false.
             okmessage(5:6) = 'no'
             exit gap_loop
          end if
       end do gap_loop
    else
       conv_criteria(3) = .true.
    end if

    ! Energy variations
    if (pub_kernel_diis_conv_criteria(4:4).eq.'1') then
       if(abs(deltaE).lt.pub_kernel_diis_threshold) then
          okmessage(7:8) = 'ok'
          conv_criteria(4) = .true.
       else
          okmessage(7:8) = 'no'
       end if
    else
       conv_criteria(4) = .true.
    end if

    converged = conv_criteria(1).and.conv_criteria(2).and.&
         conv_criteria(3).and.conv_criteria(4)


  end subroutine kernel_diis_conv_check



  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine kernel_diis_info_line(converged, okmessage, total_energy)

    !===================================================================!
    ! This subroutine prints convergence information at each DIIS iter. !
    !-------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in December 2010.                  !
    ! Modified by Alvaro Ruiz Serrano in April 2012 for easy grepping.  !
    !===================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout, VERBOSE
    use rundat, only: pub_output_detail, pub_kernel_diis_maxit, pub_num_spins

    implicit none

    logical,          intent(in   ) :: converged
    character(len=8), intent(in   ) :: okmessage
    real(kind=DP),    intent(in   ) :: total_energy

    integer :: is
    character(len=200) :: fmt

    if (pub_on_root) then

       ! ars: print header at each iteration if VERBOSE, else print at the first
       if ((pub_output_detail.ge.VERBOSE).or.(iter.eq.1)) then    !
          write(stdout,'(/,a80)') repeat('-',80)
          write(stdout,'(a80)')   "  Iter Spin   Resid.  Commut.    H-L Gap   &
               &  DGap                Energy  DEnergy"
       end if


       ! ars: print line
       do is = 1, pub_num_spins

          ! ars: write line with the correct format
          write(fmt,'(a)') &
               '("@ ", i4, 1x, i4, 2(1x,es8.1e2), 1x, es10.3e2, 1x, es8.1e2'

          if (is.eq.1) then

             ! total energy and deltaE
             if (abs(total_energy)<100000.0_DP) then
                write(fmt,'(a,a)') trim(fmt), ', f22.14, 1x, es8.1e2)'
             else
                write(fmt,'(a,a)') trim(fmt), ', f22.12, 1x, es8.1e2)'
             end if

             ! ars: write
             write(stdout,fmt) iter, is, residual(is), commut(is), &
                  gap(is), delta_gap(is), total_energy, deltaE

          else
             ! ars: no energy or deltaE for is.gt.1 - include gap and delta_gap
             write(fmt,'(a,a)') trim(fmt),')'

             ! ars: write
             write(stdout,fmt) iter, is, residual(is), commut(is),  &
                  gap(is), delta_gap(is)

          end if

       end do



       ! ars: write extra information if VERBOSE
       if (pub_output_detail.ge.VERBOSE) then

          ! ars: write OK message
          write(stdout,'("Converged? ",7x,a2,7x,a2,18x,a2,29x,a2)') &
               okmessage(1:2), okmessage(3:4), okmessage(5:6), okmessage(7:8)

          ! ars: secondary header for VERBOSE calculations
          write(stdout,'(/,a80)') "  Iter Spin  Idemp.(in)         Ne&
               &                   HOMO                   LUMO"

          ! ars: write line with the correct format
          write(fmt,'(a)')'("% ", i4, 1x, i4, 2x, es7.1e2'
          if (abs(maxval(trks))<100000.0_DP) then
             write(fmt,'(a,a)') trim(fmt),', 1x, f13.7'
          else
             write(fmt,'(a,a)') trim(fmt),', 1x, f13.5'
          end if
          if (abs(maxval(homo))<100000.0_DP) then
             write(fmt,'(a,a)') trim(fmt),', 1x, f22.14'
          else
             write(fmt,'(a,a)') trim(fmt),', 1x, f22.12'
          end if
          if (abs(maxval(lumo))<100000.0_DP) then
             write(fmt,'(a,a)') trim(fmt),', 1x, f22.14)'
          else
             write(fmt,'(a,a)') trim(fmt),', 1x, f22.12)'
          end if

          ! ars: print lines
          do is = 1, pub_num_spins
             write(stdout,fmt) iter, is, idemp(is), trks(is), homo(is), lumo(is)
          end do

       end if


       ! ars: print banner if this is the last iteration
       if(converged) then
          write(stdout,'(/a14, a30, i4, a12,/)') &
               adjustl("@CONVERGED"//ncalls_char), &
               ": calculation converged after ", iter, " iterations."
       else if ((.not.converged).and.(iter.eq.pub_kernel_diis_maxit)) then
          write(stdout,'(/,a17, a38, i4, a12,/)')&
               adjustl("@NOTCONVERGED"//ncalls_char),&
               ": calculation failed to converge after", iter, " iterations."
       end if

    end if


  end subroutine kernel_diis_info_line



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_diff_inout(residue,kern_out,kern_in)

    !====================================================================!
    ! This subroutine calculates the density kernel residue after each   !
    ! kernel DIIS iteration. R_i = K^out_i - K^in_i.                     !
    !--------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010.                      !
    ! Modified for embedding by Joseph Prentice, May 2018                !
    !====================================================================!

    use constants, only: DP
    use rundat, only: pub_num_spins
    use sparse_embed, only: sparse_embed_copy, sparse_embed_axpy, SPAM3_EMBED


    implicit none

    ! ars: arguments
    type(SPAM3_EMBED), intent(inout) :: residue(:)  ! jd: Dimensions are
    type(SPAM3_EMBED), intent(in   ) :: kern_out(:) !     pub_num_spins.
    type(SPAM3_EMBED), intent(in   ) :: kern_in(:)  !     Leave ':' be.

    ! ars: local variables
    integer :: is

    do is = 1, pub_num_spins
       ! ars: R = K^out
       call sparse_embed_copy(residue(is),kern_out(is))
       ! ars: calculate [R = K^out - K^in]
       call sparse_embed_axpy(residue(is),kern_in(is),-1.0_DP)
    end do


  end subroutine kernel_diis_diff_inout



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_hks_commutator(commutator, ham, dkn, overlap)

    use constants, only: DP
    use rundat, only: pub_num_spins
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_product,&
         sparse_axpy

    implicit none

    ! ars: Arguments
    type(SPAM3), intent(inout) :: commutator(1:pub_num_spins)
    type(SPAM3), intent(in   ) :: ham(1:pub_num_spins)
    type(SPAM3), intent(in   ) :: dkn(1:pub_num_spins)
    type(SPAM3), intent(in   ) :: overlap


    ! ars: Local variables
    integer :: is
    type(SPAM3) :: hk, kh, skh


    ! ars: create buffers
    call sparse_create(kh, dkn(1), ham(1))
    call sparse_create(skh, overlap, kh)
    call sparse_create(hk, ham(1), dkn(1))

    do is = 1, pub_num_spins

       ! ars: calculate matrices
       call sparse_product(kh, dkn(is), ham(is))
       call sparse_product(skh, overlap, kh)
       call sparse_product(hk, ham(is), dkn(is))
       call sparse_product(commutator(is), hk, overlap)
       call sparse_axpy(commutator(is), skh, -1.0_DP)

    end do

    ! ars: destroy buffers
    call sparse_destroy(hk)
    call sparse_destroy(skh)
    call sparse_destroy(kh)


  end subroutine kernel_diis_hks_commutator


  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine kernel_diis_shift(mat_his, his_size, nspins)

    !================================================================!
    ! This subroutine shifts the matrices in arrays in case iter has !
    ! reached pub_kernel_diis_size value.                            !
    !----------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2010.                    !
    ! Modified for embedding by Joseph Prentice, May 2018            !
    !================================================================!

    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy

    implicit none

    ! ars: arguments
    integer,     intent(in   ) :: his_size, nspins
    type(SPAM3_EMBED), intent(inout) :: mat_his(his_size, nspins)

    ! ars: local variables
    integer :: icntr, is


    ! ars: move old entries one position back
    do is = 1, nspins
       do icntr = 1, his_size-1
          call sparse_embed_copy(mat_his(icntr, is), mat_his(icntr+1, is))
       end do
    end do


  end subroutine kernel_diis_shift



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_find_ientry(ientry, iter, shifted)

    !==============================================================!
    ! This subroutine finds the next ientry to work with.          !
    !--------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2010.                  !
    !==============================================================!

    use comms, only: comms_barrier

    implicit none

    ! ars: arguments
    integer, intent(  out) :: ientry
    integer, intent(in   ) :: iter
    logical, intent(in   ) :: shifted


    ! ars: find next position
    if (shifted) then
       ientry = diis_size
    else
       ientry = iter+1
    end if
    call comms_barrier


  end subroutine kernel_diis_find_ientry


  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_qc_output(total_energy, mu, denskern, ham)

    !==============================================================!
    ! This subroutine prints <QC> lines related to kernel-DIIS.    !
    !--------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in December 2010.             !
    ! Modified by Alvaro Ruiz Serrano in May 2012.                 !
    ! Modified for embedding by Joseph Prentice, May 2018          !
    !==============================================================!

    use constants, only: DP, max_spins
    use comms, only: comms_barrier, pub_on_root
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_max_abs_element
    use utils, only: utils_qc_print

    implicit none

    ! ars: arguments
    real(kind=DP), intent(in   ) :: total_energy
    real(kind=DP), intent(in   ) :: mu(max_spins)
    type(SPAM3_EMBED),   intent(in   ) :: denskern(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: ham(pub_num_spins)

    ! ars: local variables
    integer :: is
    real(kind=DP) :: max_denskern(pub_num_spins)
    real(kind=DP) :: max_ham(pub_num_spins)
    character(len=3) :: spin_string


    ! ars: calculate maximum element of each matrix
    do is = 1 , pub_num_spins
       max_denskern(is) = sparse_embed_max_abs_element(denskern(is))
       max_ham(is) = sparse_embed_max_abs_element(ham(is))
    end do


    if(pub_on_root) then

       !       write(stdout, '(a30, i9)')     '<QC> [kernel_diis_iterations]:', &
       !            iter

       call utils_qc_print('total_energy',total_energy)
       call utils_qc_print('deltaE',deltaE)

       do is = 1, pub_num_spins

          write(spin_string,'(a1,i1,a1)') '(', is, ')'
          call utils_qc_print('HKS,SKH'//spin_string,commut(is))
          call utils_qc_print('K,KSK'//spin_string,idemp(is))
          call utils_qc_print('residual'//spin_string,residual(is))
          call utils_qc_print('mu'//spin_string,mu(is))
          call utils_qc_print('max(denskern)'//spin_string,max_denskern(is))
          call utils_qc_print('max(ham)'//spin_string,max_ham(is))

       end do

    end if

    call comms_barrier

  end subroutine kernel_diis_qc_output


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! ars: the following routines are used for linear, Pulay LiST and LS mixing.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




  subroutine kernel_diis_linear_mix(mat_in, mat_out, loc_dkn_in, loc_dkn_out, &
       loc_ham_in, loc_ham_out, energy_in, energy_out)

    !==============================================================!
    ! This subroutine performs linear mixing of matrices.          !
    !--------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010.                !
    ! Added ODA method by Alvaro Ruiz Serrano in April 2012.       !
    ! Modified for embedding by Joseph Prentice, May 2018          !
    !==============================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, DP, VERBOSE
    use function_basis, only: FUNC_BASIS
    use hubbard_build, only: HUBBARD_MODEL
    use rundat, only: pub_kernel_diis_coeff, pub_output_detail, &
         pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_scale, sparse_embed_axpy
    use utils, only: utils_abort

    implicit none

    ! ars: Arguments
    type(SPAM3_EMBED), intent(inout) :: mat_in(pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: mat_out(pub_num_spins)
    type(SPAM3_EMBED), optional, intent(in   ) :: loc_dkn_in(pub_num_spins)
    type(SPAM3_EMBED), optional, intent(in   ) :: loc_dkn_out(pub_num_spins)
    type(SPAM3_EMBED), optional, intent(in   ) :: loc_ham_in(pub_num_spins)
    type(SPAM3_EMBED), optional, intent(in   ) :: loc_ham_out(pub_num_spins)
    real(kind=DP), optional, intent(in   ) :: energy_in
    real(kind=DP), optional, intent(in   ) :: energy_out


    ! ars: Local variables
    real(kind=DP) :: lin_coeff(pub_num_spins)
    integer :: is
    character(len=12) :: mat_name
    logical :: opt_args

    ! ars: check present arguments
    if (present(loc_dkn_in).and.present(loc_dkn_out).and.&
         present(loc_ham_in).and.present(loc_ham_out).and.&
         present(energy_in).and.present(energy_out) ) opt_args = .true.

    ! ars: check inconsistency
    if (do_oda.and.(.not.opt_args)) then
       call utils_abort('Error in kernel_diis_linear_mix(): &
            &Not enough arguments for ODA.')
    end if

    ! ars: find linear mixing coefficient
    if (do_oda) then
       call kernel_diis_oda(lin_coeff, loc_dkn_in, loc_dkn_out, &
            loc_ham_in, loc_ham_out, energy_in, energy_out)
    else
       do is = 1, pub_num_spins
          lin_coeff(is) = pub_kernel_diis_coeff
       end do
    end if


    ! ars: set name of the matrix that we are mixing
    if (mix_dkn) then
       mat_name = 'dens kernels'
    else if (mix_ham) then
       mat_name = 'Hamiltonians'
    end if

    ! ars: banner
    if(pub_on_root.and.pub_output_detail.ge.VERBOSE) then
       if (do_oda) then
          if (pub_num_spins.eq.1) then
             write(stdout,'(/,a25, a12, a14, f8.6, a4)',advance ='no') &
                  "Performing ODA mixing of ", mat_name, " with coeff = ",&
                  lin_coeff(1), " ..."
          else if (pub_num_spins.eq.2) then
             write(stdout,'(/, a25, a12, a14, f8.6, a2, f8.6, a4)',advance ='no') &
                  "Performing ODA mixing of ", mat_name, " with coeff = ",&
                  lin_coeff(1),", ", lin_coeff(2), " ..."
          end if
       else
          if (pub_num_spins.eq.1) then
             write(stdout,'(/,a28, a12, a14, f8.6, a4)',advance ='no') &
                  "Performing linear mixing of ", mat_name, " with coeff = ",&
                  lin_coeff(1), " ..."
          else if (pub_num_spins.eq.2) then
             write(stdout,'(/, a28, a12, a14, f8.6, a2, f8.6, a4)',advance ='no') &
                  "Performing linear mixing of ", mat_name, " with coeff = ",&
                  lin_coeff(1),", ", lin_coeff(2), " ..."
          end if
       end if
    end if


    ! ars: do linear mixing of density kernels
    do is = 1, pub_num_spins
       ! ars: M_in(n+1) = c_in*M_in(n)
       call sparse_embed_scale(mat_in(is), 1.0_DP - lin_coeff(is))
       ! ars: M_in(n+1) = c_out*M_out(n) + c_in*M_in(n)
       call sparse_embed_axpy(mat_in(is), mat_out(is), lin_coeff(is))
    end do

    ! ars: banner
    if(pub_on_root.and.pub_output_detail.ge.VERBOSE) then
       write(stdout,'(a,/)') "done"
    end if



  end subroutine kernel_diis_linear_mix


  !_____________________________________________________________________________
  !_____________________________________________________________________________




  subroutine kernel_diis_oda(lin_coeff, loc_dkn_in, loc_dkn_out, &
       loc_ham_in, loc_ham_out, energy_in, energy_out)

    !=========================================================================!
    ! Calculates Optimal Damping Algorithm (ODA) linear mixing coefficient.   !
    ! See "Cances, J. Chem. Phys. 114(24), 2001".                             !
    ! Warning: under development (use with caution).                          !
    !-------------------------------------------------------------------------!
    ! Created by Alvaro Ruiz Serrano in April 2012.                           !
    ! Modified for embedding by Joseph Prentice, May 2018                     !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_copy, &
         sparse_embed_axpy, sparse_embed_trace


    implicit none

    ! ars: arguments
    real(kind=DP), intent(  out) :: lin_coeff(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: loc_dkn_in(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: loc_dkn_out(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: loc_ham_in(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: loc_ham_out(pub_num_spins)
    real(kind=DP), intent(in   ) :: energy_in
    real(kind=DP), intent(in   ) :: energy_out

    integer :: is
    real(kind=DP) :: afac, bfac, cfac, dfac, l_sqrt, c_pos, c_neg
    type(SPAM3_EMBED) :: dkn_diff

    real(kind=DP), parameter :: threshold = 0.001_DP

    ! ars: create buffers
    call sparse_embed_create(dkn_diff, dkn_out(1))

    do is = 1, pub_num_spins

       ! ars: calculate dkn differential
       call sparse_embed_copy(dkn_diff, loc_dkn_out(is))
       call sparse_embed_axpy(dkn_diff, loc_dkn_in(is), -1.0_DP)

       ! ars: set dfac
       dfac = energy_in

       ! ars: set cfac
       call sparse_embed_trace(cfac,loc_ham_in(is),dkn_diff)
       cfac = 2.0_DP*cfac

       ! ars: set afac
       call sparse_embed_trace(afac,loc_ham_out(is),dkn_diff)
       afac = 2.0_DP*afac - 2.0_DP*energy_out + cfac + 2.0_DP*dfac

       ! ars: set bfac
       bfac = energy_out - afac - cfac - dfac

       ! ars: calculate optimum lin_coeff
       l_sqrt = abs(sqrt(bfac*bfac - 3.0_DP*afac*cfac))

       ! ars: calculate the two possible coefficients
       c_pos = (-bfac + l_sqrt)/(3.0_DP*afac)
       c_neg = (-bfac - l_sqrt)/(3.0_DP*afac)

       ! ars: set lin_coeff
       if ( (c_pos.ge.0.0_DP).and.(c_pos.le.1.0_DP).and.&
            (c_neg.ge.0.0_DP).and.(c_neg.le.1.0_DP) ) then
          ! ars: both are in the range. Take the greatest
          lin_coeff(is) = max(c_pos,c_neg)
       else if ( (c_pos.ge.0.0_DP).and.(c_pos.le.1.0_dP) ) then
          ! ars: c_pos is in range
          lin_coeff(is) = c_pos
       else if ( (c_neg.ge.0.0_DP).and.(c_neg.le.1.0_dP) ) then
          ! ars: c_neg is in range
          lin_coeff(is) = c_neg
       else
          ! ars: none of the coeffs are in range - set to safe value
          if (max(c_pos,c_neg).lt.0.0_DP) then
             lin_coeff(is) = 0.1_DP
             if (pub_on_root) write(stdout,*) &
                  "ODA can't converge - setting coeff = 0.1"
          else if (max(c_pos, c_neg).gt.1.0_DP) then
             lin_coeff(is) = 1.0_DP
          end if
       end if

       ! ars: sanity check against too small ODA coeff
       if (lin_coeff(is).lt.threshold) then
          if (pub_on_root) write(stdout,*) &
               "Coefficient too small - setting to safe value coeff = 0.1"
          lin_coeff(is) = 0.1_DP
       end if

    end do

    ! ars: destroy buffers
    call sparse_embed_destroy(dkn_diff)


  end subroutine kernel_diis_oda

  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine kernel_diis_pulay_mix(mat_in, mat_his, err_his, overlap, ientry)

    !==============================================================!
    ! This subroutine performs Pulay mixing of matrices.           !
    !--------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010.                !
    ! Modified for embedding by Joseph Prentice, May 2018          !
    !==============================================================!


    use constants, only: DP
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! ars: Arguments
    type(SPAM3_EMBED), intent(inout) :: mat_in(pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: mat_his(diis_size,pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: err_his(diis_size,pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: overlap
    integer,     intent(in   ) :: ientry

    ! ars: local variables
    real(kind=DP), allocatable :: diis_mat(:,:,:)
    real(kind=DP), allocatable :: diis_coeffs(:,:)
    integer :: ierr

    ! ars: allocate diis_mat and coeffs
    allocate(diis_mat(1:ientry+1, 1:ientry+1, 1:pub_num_spins),stat=ierr)
    call utils_alloc_check('kernel_diis_pulay_mix', 'diis_mat', ierr)
    allocate(diis_coeffs(1:ientry+1, 1:pub_num_spins), stat=ierr)
    call utils_alloc_check('kernel_diis_pulay_mix', 'diis_coeffs', ierr)


    ! ars: calculate DIIS matrix based on the error vectors.
    call kernel_diis_pulay_matrix(diis_mat, err_his, overlap, ientry)

    ! ars: get DIIS coefficients
    call kernel_diis_coeffs(diis_coeffs, diis_mat, ientry)
    call kernel_diis_coeffs_summ(diis_coeffs, ientry)

    ! ars: mix kernels
    call kernel_diis_mix_kernels(mat_in, mat_his, diis_coeffs, ientry)

    ! ars: deallocate Bmat and coeffs
    deallocate(diis_coeffs, stat=ierr)
    call utils_dealloc_check('kernel_diis_pulay_mix', 'diis_coeffs', ierr)
    deallocate(diis_mat, stat=ierr)
    call utils_dealloc_check('kernel_diis_pulay_mix', 'diis_mat', ierr)



  end subroutine kernel_diis_pulay_mix

  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine kernel_diis_list_mix(mat_in, mat_his, ham_diff, dkn_diff, ientry)


    !====================================================================!
    ! This subroutine performs LiST (LiSTi or LiSTb) mixing of matrices. !
    !--------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2012.                        !
    ! Modified for embedding by Joseph Prentice, May 2018                !
    !====================================================================!

    use constants, only: DP
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! ars: Arguments
    type(SPAM3_EMBED), intent(inout) :: mat_in(pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: mat_his(diis_size,pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: ham_diff(diis_size,pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: dkn_diff(diis_size,pub_num_spins)
    integer,     intent(in   ) :: ientry

    ! ars: local variables
    real(kind=DP), allocatable :: diis_mat(:,:,:)
    real(kind=DP), allocatable :: diis_coeffs(:,:)
    integer :: ierr

    ! ars: allocate diis_mat and coeffs
    allocate(diis_mat(1:ientry+1, 1:ientry+1, 1:pub_num_spins),stat=ierr)
    call utils_alloc_check('kernel_diis_list_mix', 'diis_mat', ierr)
    allocate(diis_coeffs(1:ientry+1, 1:pub_num_spins), stat=ierr)
    call utils_alloc_check('kernel_diis_list_mix', 'diis_coeffs', ierr)


    ! ars: calculate LiST matrix based on the method.
    if (do_listi) then
       call kernel_diis_listi_matrix(diis_mat, dkn_diff, ham_diff, ientry)
    else if (do_listb) then
       call kernel_diis_listb_matrix(diis_mat, dkn_diff, ham_diff, ientry)
    end if

    ! ars: get DIIS coefficients
    call kernel_diis_coeffs(diis_coeffs, diis_mat, ientry)
    call kernel_diis_coeffs_summ(diis_coeffs, ientry)

    ! ars: mix kernels
    call kernel_diis_mix_kernels(mat_in, mat_his, diis_coeffs, ientry)

    ! ars: deallocate Bmat and coeffs
    deallocate(diis_coeffs, stat=ierr)
    call utils_dealloc_check('kernel_diis_list_mix', 'diis_coeffs', ierr)
    deallocate(diis_mat, stat=ierr)
    call utils_dealloc_check('kernel_diis_list_mix', 'diis_mat', ierr)



  end subroutine kernel_diis_list_mix

  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine kernel_diis_pulay_matrix(Bmat,residues, overlap, ientry)

    !==============================================================!
    ! This subroutine calculates matrix B of the system of linear  !
    ! equations to be solved in order to mix the density kernels   !
    ! according to the Pulay method.                               !
    !--------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in April 2010.                 !
    ! Modified for embedding by Joseph Prentice, May 2018          !
    !==============================================================!

    use constants, only: DP
    use comms, only: comms_barrier
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_transpose, &
         sparse_embed_product, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_trace


    implicit none

    ! ars: arguments
    real(kind=DP), intent(  out) :: Bmat(:,:,:)
    type(SPAM3_EMBED),   intent(in   ) :: residues(diis_size, pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: overlap
    integer,       intent(in   ) :: ientry

    ! ars: local variables
    integer     :: is, ii, jj
    type(SPAM3_EMBED) :: res_ov, res_t, res_ov_res_t


    ! ars: init spam3buffer blocking scheme
    call sparse_embed_create(res_ov, residues(1,1), overlap)
    call sparse_embed_create(res_t, residues(1,1))
    call sparse_embed_create(res_ov_res_t, res_ov, res_t)


    do is = 1, pub_num_spins
       do jj = 1, ientry
          do ii = 1, ientry

             ! ars: res_ov = R_i x S
             call sparse_embed_product(res_ov, residues(ii, is), overlap)

             ! ars: res_t = (R_j)^t
             call sparse_embed_transpose(res_t,residues(jj, is))

             ! ars: res_ov_res_t = (R_i) x S x (R_j)^t
             call sparse_embed_product(res_ov_res_t, res_ov, res_t)

             ! ars: Bmat(is,ii,jj) = tr[(R_i) x S x (R_j)^t x S]
             call sparse_embed_trace(Bmat(ii,jj,is), res_ov_res_t,overlap)

          end do
       end do
    end do

    ! ars: destroy SPAM3 buffers
    call sparse_embed_destroy(res_ov)
    call sparse_embed_destroy(res_t)
    call sparse_embed_destroy(res_ov_res_t)

    ! ars: add space for Langrange multiplier
    Bmat(1:ientry,ientry+1,:) = -1.0_DP
    Bmat(ientry+1,1:ientry,:) = -1.0_DP
    Bmat(ientry+1,ientry+1,:) = 0.0_DP



  end subroutine kernel_diis_pulay_matrix

  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine kernel_diis_listi_matrix(Bmat, residues, ham_residues, ientry)

    !==============================================================!
    ! This subroutine calculates matrix B of the system of linear  !
    ! equations to be solved in order to mix the density kernels   !
    ! according to the LiSTi method.                               !
    ! See J. Chem. Phys. (134) 241102, 2011.                       !
    !--------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2012.                !
    ! Modified for embedding by Joseph Prentice, May 2018          !
    !==============================================================!

    use constants, only: DP
    use comms, only: comms_barrier
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_trace


    implicit none

    ! ars: arguments
    real(kind=DP), intent(  out) :: Bmat(:,:,:)
    type(SPAM3_EMBED),   intent(in   ) :: residues(diis_size, pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: ham_residues(diis_size, pub_num_spins)
    integer,       intent(in   ) :: ientry

    ! ars: local variables
    integer     :: is, ii, jj

    ! ars: Bmat(ii,jj,is) = tr[(K_j^out - K_j^in) * (H_i^out - H_i^in)/2]
    do is = 1, pub_num_spins
       do jj = 1, ientry
          do ii = 1, ientry
             call sparse_embed_trace(Bmat(ii,jj,is),residues(jj,is),ham_residues(ii,is))
          end do
       end do
    end do

    ! ars: scale according to LiSTi
    Bmat(:,:,:) = 0.5_DP*Bmat(:,:,:)

    ! ars: add space for Langrange multiplier
    Bmat(1:ientry,ientry+1,:) = -1.0_DP
    Bmat(ientry+1,1:ientry,:) = -1.0_DP
    Bmat(ientry+1,ientry+1,:) = 0.0_DP


  end subroutine kernel_diis_listi_matrix

  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine kernel_diis_listb_matrix(Bmat, residues, ham_residues, ientry)

    !==============================================================!
    ! This subroutine calculates matrix B of the system of linear  !
    ! equations to be solved in order to mix the density kernels   !
    ! according to the LiSTi method.                               !
    ! See J. Chem. Phys. (134) 241102, 2011.                       !
    !--------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2012.                !
    ! Modified for embedding by Joseph Prentice, May 2018          !
    !==============================================================!

    use constants, only: DP
    use comms, only: comms_barrier
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_trace


    implicit none

    ! ars: arguments
    real(kind=DP), intent(  out) :: Bmat(:,:,:)
    type(SPAM3_EMBED),   intent(in   ) :: residues(diis_size, pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: ham_residues(diis_size, pub_num_spins)
    integer,       intent(in   ) :: ientry

    ! ars: local variables
    integer     :: is, ii, jj

    ! ars: Bmat(ii,jj,is) = tr[(K_j^out - K_j^in) * (H_i^out - H_i^in)/2]
    do is = 1, pub_num_spins
       do jj = 1, ientry
          do ii = 1, ientry
             call sparse_embed_trace(Bmat(ii,jj,is),residues(ii,is),ham_residues(jj,is))
          end do
       end do
    end do

    ! ars: scale according to LiSTb
    Bmat(:,:,:) = 0.5_DP*Bmat(:,:,:)

    ! ars: add space for Langrange multiplier
    Bmat(1:ientry,ientry+1,:) = -1.0_DP
    Bmat(ientry+1,1:ientry,:) = -1.0_DP
    Bmat(ientry+1,ientry+1,:) = 0.0_DP


  end subroutine kernel_diis_listb_matrix

  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_coeffs(coeffs,Bmat,ientry)

    !==============================================================!
    ! This subroutine finds the coefficients for the Pulay mixing  !
    ! after solving a system of linear equations Bd=0.             !
    !--------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2010.                  !
    !==============================================================!

    use comms, only: pub_my_proc_id, pub_on_root, comms_barrier
    use constants, only: DP, stdout
    use rundat, only: pub_num_spins


    implicit none

    ! ars: arguments
    real(kind=DP), intent(inout) :: coeffs(:,:)
    real(kind=DP), intent(in   ) :: Bmat(:,:,:)
    integer      , intent(in   ) :: ientry

    ! LAPACK subroutines
    external :: dgetrf, dgetrs

    ! ars: library wrapper variables
    integer :: INFO, LDA, LDB, M, N, NRHS
    integer :: IPIV(1:ientry+1)
    real(kind=DP) :: A(ientry+1,ientry+1)
    real(kind=DP) :: B(ientry+1)
    character(LEN=1) :: TRANS

    ! ars: local variables
    integer :: is



    ! ars: set up parameters for LAPACK DGETRF
    M = ientry+1
    N = ientry+1
    LDA = ientry+1
    ! ars: set up parameters for LAPACK DGETRS
    TRANS ='N'
    NRHS = 1
    LDB = ientry+1

    ! ars: solve linear system
    do is = 1, pub_num_spins

       ! ars: call LAPACK DGETRF
       A = Bmat(:,:,is)
       call DGETRF(M, N, A, LDA, IPIV, INFO)
       if((INFO.ne.0).and.pub_on_root) write(stdout,*) &
            "Error calling DGETRF. INFO|proc = ", INFO, pub_my_proc_id

       ! ars: call LAPACK DGETRS
       B(1:ientry) = 0.0_DP
       B(ientry+1) = -1.0_DP
       call DGETRS(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       if((INFO.ne.0).and.pub_on_root) write(stdout,*) &
            "Error calling DGETRS. INFO|proc = ", INFO, pub_my_proc_id

       ! ars: set coeffs(is) before exit
       coeffs(:, is) = B(:)

    end do

    call comms_barrier

  end subroutine kernel_diis_coeffs



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_coeffs_summ(coeffs, ientry)

    !==================================================================!
    ! This subroutine prints a summary of the coefficients that are    !
    ! used during the density kernel mixing.                           !
    !------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in July 2010.                     !
    !==================================================================!

    use comms, only : pub_on_root
    use constants, only: DP, stdout, VERBOSE
    use rundat, only: pub_output_detail, pub_num_spins

    implicit none

    real(kind=DP), intent(in   ) :: coeffs(:,:)
    integer,       intent(in   ) :: ientry

    integer :: icntr
    character(len=15) :: mat_name
    character(len=5) :: method_name
    character(len=80) :: fmt


    ! ars: prepare format
    write(fmt,'(a)') '(/, a11, a5, a11'

    ! ars: method name
    if (do_pulay) then
       method_name = "Pulay"
    else if (do_listb) then
       method_name = "LISTB"
    else if (do_listi) then
       method_name = "LISTI"
    end if

    ! ars: matrix name
    if (mix_dkn) then
       mat_name = 'density kernels'
       write(fmt,'(a,a)') trim(fmt), ', a15, a19)'
    else if (mix_ham) then
       mat_name = 'Hamiltonians'
       mat_name = trim(adjustl(mat_name))
       write(fmt,'(a,a)') trim(fmt), ', a12, a19)'
    end if

    ! ars: banner
    if (pub_on_root.and.pub_output_detail.ge.VERBOSE) then
       write(stdout,fmt) &
            "Performing ", method_name, " mixing of ", &
            mat_name, " with coefficients:"
       if (pub_num_spins.eq.1) then
          write(stdout,'(10x,a20)') "              Spin 1"
          do icntr = 1, ientry
             write(stdout,'(10x,f20.12)') coeffs(icntr,1)
          end do
          write(stdout,'(a10,f20.12)',advance='no') "Lambda--> ", coeffs(ientry+1,1)
       else if (pub_num_spins.eq.2) then
          write(stdout,'(10x,a20,1x,a20)') "              Spin 1", &
               "              Spin 2"
          do icntr = 1, ientry
             write(stdout,'(10x,f20.12,1x,f20.12)') coeffs(icntr,1), coeffs(icntr,2)
          end do
          write(stdout,'(a10,f20.12,1x,f20.12)',advance='no') &
               "Lambda--> ", coeffs(ientry+1,1), coeffs(ientry+1,2)
       end if
       write(stdout,'(a,/)') "     ... done"
    end if


  end subroutine kernel_diis_coeffs_summ



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_mix_kernels(next_dkn_in, dkn_out, coeffs, ientry)

    !==================================================================!
    ! This subroutine mixes the kernels according to the Pulay method. !
    !------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in May 2010.                      !
    ! Modified for embedding by Joseph Prentice, May 2018              !
    !==================================================================!

    use constants, only: DP
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_axpy, sparse_embed_scale

    implicit none

    type(SPAM3_EMBED),   intent(inout) :: next_dkn_in(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: dkn_out(diis_size, pub_num_spins)
    real(kind=DP), intent(in   ) :: coeffs(:,:)
    integer,       intent(in   ) :: ientry

    integer :: is, counter

    ! ars: K_in(n+1) = 0
    do is = 1, pub_num_spins
       call sparse_embed_scale(next_dkn_in(is),0.0_DP)
    end do

    ! ars: K_in(n+1) = sum_i K_out(i) * d_i
    do is = 1, pub_num_spins
       do counter = 1, ientry
          call sparse_embed_axpy(next_dkn_in(is),dkn_out(counter,is),&
               coeffs(counter,is))
       end do
    end do


  end subroutine kernel_diis_mix_kernels





  !_____________________________________________________________________________
  !_____________________________________________________________________________


  !==========================================================================!
  ! Subroutines for level shifting                                           !
  !==========================================================================!


  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_ls_init(beta, ham, overlap, n_occ, num)

    use comms, only: pub_on_root
    use constants, only: DP, stdout, max_spins, VERBOSE
    use rundat, only: pub_output_detail, pub_kernel_diis_lshift, &
         pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy

    implicit none

    ! ars: Arguments
    real(kind=DP), intent(  out) :: beta
    type(SPAM3_EMBED),   intent(in   ) :: ham(1:pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: overlap
    integer,       intent(in   ) :: n_occ(max_spins)
    integer,       intent(in   ) :: num


    ! ars: local variables
    integer :: is


    if (.not.ls_dkn_eigs) then

       ! ars: diagonalise input ham and store eigenvalues
       call kernel_diis_diag_ham(eigvecs, ham, overlap, n_occ, num, &
            homo, lumo, gap)

       ! ars: print initial homo-lumo gap on screen if required
       if (pub_on_root.and.pub_output_detail.ge.VERBOSE) then
          write(stdout,'(a)') "UNSHIFTED HAM | HOMO     | LUMO     | GAP      "
          do is = 1, pub_num_spins
             write(stdout,'(11x, i2, 1x, 3(1x,es10.3e2))') &
                  is, homo(is), lumo(is), gap(is)
          end do
          write(stdout,'(a)')
       end if

    end if

    ! ars: set beta
    beta = pub_kernel_diis_lshift

  end subroutine kernel_diis_ls_init

  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine kernel_diis_build_ls_mat(beta, ls_matrix, eigs_dens, overlap, &
       n_occ, num)

    !==================================================================!
    ! This 1) Construcuts the levelshifting matrix:                    !
    !         Delta = SMBMS                                            !
    !------------------------------------------------------------------!
    ! Written by Peter John Cherry and Alvaro Ruiz Serrano in May 2012 !
    ! Modified for embedding by Joseph Prentice, May 2018              !
    !==================================================================!


    use comms, only: pub_on_root
    use constants, only: stdout, DP, max_spins, VERBOSE
    use dense, only: dense_product, dense_scale, dense_create, &
         dense_destroy, dense_convert
    use rundat, only: pub_output_detail, pub_num_spins, pub_ngwf_regions_ngroups
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_axpy, sparse_embed_product, sparse_embed_destroy



    implicit none

    !Arguments
    real(kind=DP), intent(inout) :: beta
    type(SPAM3_EMBED),   intent(inout) :: ls_matrix(1:pub_num_spins)
    type(DEM),     intent(in   ) :: eigs_dens(1:pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: overlap
    integer,       intent(in   ) :: n_occ(max_spins)
    integer,       intent(in   ) :: num


    !Local Variables
    integer :: is, nsub
    ! agrecocmplx
    logical :: loc_cmplx
    character :: opB_loc

    ! Level shifting matrix
    type(DEM) :: mm_ls_matrix_dens

    ! SPAM3 level_shift_mat
    type(SPAM3_EMBED) :: mm_ls_matrix, ds


    ! ars: banner
    if(pub_on_root.and.pub_output_detail.ge.VERBOSE) then
       write(stdout,'(a43,f8.6,/)') &
            "Performing level shifting with parameter = ", beta
    end if

    ! agrecocmplx
    loc_cmplx = overlap%p%iscmplx
    ! jcap: set local copy of number of sub regions
    nsub=pub_ngwf_regions_ngroups

    ! agrecocmplx
    if (loc_cmplx) then
       opB_loc = 'C'
    else
       opB_loc = 'T'
    end if

    ! ars: create structures
    ! agrecocmplx
    call dense_create(mm_ls_matrix_dens, num, num,iscmplx=loc_cmplx)
    mm_ls_matrix%structure='K'
    ! agrecocmplx
    call sparse_embed_create(mm_ls_matrix,iscmplx=loc_cmplx)
    call sparse_embed_create(ds, mm_ls_matrix, overlap)

    do is=1, pub_num_spins


       ! cks: ensure the eigenvectors are zero initially
       call dense_scale(mm_ls_matrix_dens, 0.0_DP)

       ! pjc: construct a levelshifting matrix from the last (num-n_occ(is))
       ! pjc: eigenvectors
       ! agrecocmplx: need to take hermitian of matrix B in complex case?
       ! I think so....
       call dense_product(mm_ls_matrix_dens,eigs_dens(is),eigs_dens(is), &
            opA='N',opB=opB_loc,first_k=n_occ(is)+1,last_k=num)

       ! pjc: multiply this matrix by level shift value
       call dense_scale(mm_ls_matrix_dens, beta)

       ! ars: convert back to SPAM3
       call dense_convert(mm_ls_matrix, mm_ls_matrix_dens)

       ! ars: calculate DS
       call sparse_embed_product(ds, mm_ls_matrix, overlap)

       ! calculate sds
       call sparse_embed_product(ls_matrix(is), overlap, ds)
    end do


    ! ars: destroy structures
    call sparse_embed_destroy(ds)
    call sparse_embed_destroy(mm_ls_matrix)
    call dense_destroy(mm_ls_matrix_dens)

  end subroutine kernel_diis_build_ls_mat


  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine kernel_diis_dkn_diag(eigs_dens, dkn, overlap, num)

    !==============================================================!
    ! This subroutine gets the eigenvectors from the kernel to be  !
    ! in the generation of the level shifting matrix               !
    !==============================================================!
    ! Peter John Cherry May 2012                                   !
    ! Modified for embedding by Joseph Prentice, May 2018          !
    !==============================================================!

    use constants, only: DP
    use dense, only: dense_create, dense_destroy, dense_convert,&
         dense_eigensolve, dense_scale
    use sparse_embed, only: SPAM3_EMBED
    use utils, only: utils_alloc_check, utils_dealloc_check
    use rundat, only: pub_num_spins, pub_eigensolver_orfac, &
         pub_eigensolver_abstol

    implicit none

    ! ars: Arguments
    type(DEM),   intent(inout) :: eigs_dens(pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: dkn(pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: overlap
    integer,     intent(in   ) :: num

    ! ars: local variables
    type(DEM) :: overlap_dense, dkn_in_dense
    real(kind=DP), allocatable :: occupancies(:)
    integer :: ierr, is
    ! agrecocmplx
    logical :: loc_cmplx

    ! agrecocmplx
    loc_cmplx = overlap%p%iscmplx

    ! ars: allocate occupancies
    allocate(occupancies(num),stat=ierr)
    call utils_alloc_check('kernel_diis_dkn_diag','occupancies',ierr)


    ! ndmh: create temporary dense matrices
    ! agrecocmplx
    call dense_create(dkn_in_dense,num,num,iscmplx=loc_cmplx)
    call dense_create(overlap_dense,num,num,iscmplx=loc_cmplx)

    ! cks: loop over spins
    do is=1,pub_num_spins

       ! ars: initialise eigenvectors to zero
       call dense_scale(eigs_dens(is), 0.0_DP)

       ! ndmh: expand the density kernel to a dense matrix
       call dense_convert(dkn_in_dense,dkn(is))

       ! ars: multiply by -1 to get eigenvecs in increasing order of energy
       call dense_scale(dkn_in_dense,-1.0_DP)

       ! ndmh: expand the overlap matrix to a dense matrix
       call dense_convert(overlap_dense,overlap)

       ! ndmh: solve the generalised eigenvalue problem
       call dense_eigensolve(num,occupancies,dkn_in_dense,overlap_dense, &
            2,eigs_dens(is),pub_eigensolver_orfac, pub_eigensolver_abstol)

    end do

    ! ndmh: clean up temporary matrices
    call dense_destroy(overlap_dense)
    call dense_destroy(dkn_in_dense)

    ! ars: deallocate occupancies
    deallocate(occupancies,stat=ierr)
    call utils_dealloc_check('kernel_diis_dkn_diag','occupancies',ierr)

  end subroutine kernel_diis_dkn_diag


  !_____________________________________________________________________________
  !_____________________________________________________________________________

  subroutine kernel_diis_level_shifter(dkn_out, dkn_in, ham, overlap, n_occ)


    use comms, only: pub_on_root
    use constants, only: DP, max_spins, stdout
    use rundat, only: pub_kernel_diis_ls_iter, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_copy, sparse_embed_axpy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    type(SPAM3_EMBED),   intent(inout) :: dkn_out(pub_num_spins)
    type(SPAM3_EMBED),   intent(in) :: dkn_in(pub_num_spins)
    type(SPAM3_EMBED),   intent(inout) :: ham(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: overlap
    integer,             intent(in   ) :: n_occ(max_spins)



    ! ars: local variables
    type(SPAM3_EMBED), allocatable :: ls_matrix(:)
    type(SPAM3_EMBED), allocatable :: shifted_ham(:)

    integer :: is, ierr

    ! ars: start timer
    call timer_clock('kernel_diis_level_shifter',1)


    ! ars: matrix to add to ham and shift eigenvalues
    allocate(ls_matrix(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('kernel_diis_level_shifter', 'ls_matrix', ierr)
    do is = 1, pub_num_spins
       call sparse_embed_create(ls_matrix(is), ham(is))
    end do

    ! ars: shifted hamiltonian
    allocate(shifted_ham(1:pub_num_spins), stat=ierr)
    call utils_alloc_check('kernel_diis_level_shifter', 'shifted_ham', ierr)
    do is = 1, pub_num_spins
       call sparse_embed_create(shifted_ham(is), ham(is))
    end do


    ! ars: use dkn_eigs to build ls_matrix if required
    if (ls_dkn_eigs) then
       call kernel_diis_dkn_diag(eigvecs, dkn_in, overlap, num)
       if (pub_on_root) write(stdout,'(a)') &
            "Level shifter: diagonalising density kernel"
    end if

    ! ars: build LS matrix from eigenvectors of previous iteration (or init)
    call kernel_diis_build_ls_mat(beta, ls_matrix, eigvecs, overlap,&
         n_occ, num)

    ! ars: build shifted Hamiltonian
    do is = 1, pub_num_spins
       call sparse_embed_copy(shifted_ham(is), ham(is))
       call sparse_embed_axpy(shifted_ham(is), ls_matrix(is), 1.0_DP)
    end do

    ! ars: diagonalise shifted Hamiltonian and store eigenvectors for next iteration
    call kernel_diis_diag_ham(eigvecs, shifted_ham, overlap, n_occ, num,&
         homo, lumo, gap)

    ! ars: build new density kernel from eigenvectors
    call kernel_diis_dkn_from_eigs(dkn_out, eigvecs, n_occ, num)

    ! ars: destroy shifted ham
    do is = 1, pub_num_spins
       call sparse_embed_destroy(shifted_ham(is))
    end do
    deallocate(shifted_ham, stat=ierr)
    call utils_dealloc_check('kernel_diis_level_shifter', 'shifted_ham', ierr)

    ! ars: destroy ls_matrix
    do is = 1, pub_num_spins
       call sparse_embed_destroy(ls_matrix(is))
    end do
    deallocate(ls_matrix, stat=ierr)
    call utils_dealloc_check('kernel_diis_level_shifter', 'ls_matrix', ierr)


    ! ars: stop LS if reached max iter
    if (iter.ge.pub_kernel_diis_ls_iter) then
       call kernel_diis_ls_dealloc()
       do_ls = .false.
    end if

    ! ars: stop timer
    call timer_clock('kernel_diis_level_shifter',2)


  end subroutine kernel_diis_level_shifter


  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !
  !                         END PRIVATE COMMON ROUTINES
  !_____________________________________________________________________________
  !_____________________________________________________________________________




  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !
  !                            PUBLIC COMMON ROUTINES
  !_____________________________________________________________________________
  !_____________________________________________________________________________




  subroutine kernel_diis_lagrangian(electronic_lagrangian, overlap, denskern, &
       ham)

    !==========================================================================!
    ! Evaluates the electronic Lagrangian according to kernel-DIIS.            !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in January 2011.                          !
    ! Modified for embedding by Joseph Prentice, May 2018                      !
    !==========================================================================!

    use constants, only: DP
    use rundat, only: pub_num_spins, pub_spin_fac
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_trace, &
         sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_copy


    implicit none


    ! Arguments
    real(kind=DP), intent(inout) :: electronic_lagrangian
    type(SPAM3_EMBED),   intent(in   ) :: overlap
    type(SPAM3_EMBED),   intent(in   ) :: denskern(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: ham(pub_num_spins)

    ! Local variable
    integer :: is
    real(kind=DP) :: trhksk
    type(SPAM3_EMBED) :: ksk, ks


    call sparse_embed_create(ks, denskern(1), overlap)
    call sparse_embed_create(ksk, ks, denskern(1))


    ! ars: apply kernel DIIS lagrangian
    !    : at this point, K already contains the rescaling factor
    do is =1, pub_num_spins

       ! ars: calculate tr[H(KSK-K)]
       call sparse_embed_product(ks, denskern(is), overlap)
       call sparse_embed_product(ksk, ks, denskern(is))
       call sparse_embed_axpy(ksk, denskern(is), -1.0_DP)
       call sparse_embed_trace(trhksk, ham(is), ksk)

       ! ars: apply Lagrangian
       electronic_lagrangian = electronic_lagrangian - pub_spin_fac * trhksk

    end do

    call sparse_embed_destroy(ksk)
    call sparse_embed_destroy(ks)


  end subroutine kernel_diis_lagrangian

  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine kernel_diis_build_pq(pmat, qmat, denskern, ham, overlap, n_occ, mu)

    !=====================================================================!
    ! This subroutine returns the P and Q matrices that multiply the H-fb !
    ! and the fb-only parts of the gradient of the kernel DIIS energy as  !
    ! a function of the NGWF expansion coefficients in the psinc basis.   !
    !---------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in January 2011.                     !
    !=====================================================================!

    use constants, only: DP, max_spins
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, sparse_embed_product, &
         sparse_embed_copy, sparse_embed_scale, sparse_embed_axpy, sparse_embed_trace


    implicit none


    ! ars: arguments
    type(SPAM3_EMBED),   intent(inout) :: pmat(pub_num_spins)
    type(SPAM3_EMBED),   intent(inout) :: qmat(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: denskern(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: ham(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: overlap
    integer,       intent(in   ) :: n_occ(max_spins)
    real(kind=DP), intent(in   ) :: mu(max_spins)

    ! ars: local variables
    integer :: is
    type(SPAM3_EMBED) :: kh, ks, ksk
    real(kind=DP) :: nu, rescale, trace


    ! ars: allocate workspace
    call sparse_embed_create(kh,denskern(1), ham(1))
    call sparse_embed_create(ks,denskern(1), overlap)
    call sparse_embed_create(ksk, ks, denskern(1))

    do is = 1, pub_num_spins

       ! ars: build buffers
       call sparse_embed_product(kh, denskern(is), ham(is))
       call sparse_embed_product(ks, denskern(is), overlap)
       call sparse_embed_product(ksk, ks, denskern(is))

       ! ars: build P = [2Ne/tr(KS)] *K - [Ne/tr(KS)]^2 *KSK
       call sparse_embed_trace(rescale, ks)
       rescale = real(n_occ(is), kind=DP) / rescale
       call sparse_embed_copy(pmat(is), denskern(is))
       call sparse_embed_scale(pmat(is), 2*rescale)
       call sparse_embed_axpy(pmat(is), ksk, -rescale*rescale)

       ! ars: build Q = -[Ne/tr(KS)]^2 *KHK - 2(mu-[Ne/tr(KS)]^2 * nu] *K
       call sparse_embed_trace(nu, ham(is), ksk)
       call sparse_embed_trace(trace, ks)
       nu = nu / trace
       call sparse_embed_product(qmat(is), kh, denskern(is))
       call sparse_embed_scale(qmat(is), rescale*rescale)
       call sparse_embed_axpy(qmat(is), denskern(is), 2.0_DP*(mu(is)-rescale**2*nu))
       call sparse_embed_scale(qmat(is), -1.0_DP)

    end do

    ! ars: deallocate workspace
    call sparse_embed_destroy(kh)
    call sparse_embed_destroy(ks)
    call sparse_embed_destroy(ksk)


  end subroutine kernel_diis_build_pq


  !_____________________________________________________________________________
  !_____________________________________________________________________________


  subroutine kernel_diis_build_idemp_dkn(denskern, ham, overlap, n_occ, &
       Ehomo, Elumo, Egap)

    !=====================================================================!
    ! This subroutine diagonalises the hamiltonian and builds the output  !
    ! density from the eigenvectors.                                      !
    !---------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010 as an adapation of the !
    ! existing routine kernel_reset_occupancies.                          !
    ! Modified for embedding by Joseph Prentice, May 2018                 !
    !=====================================================================!

    use constants, only: DP, max_spins
    use dense, only: DEM, dense_create, dense_convert, dense_destroy, &
         dense_eigensolve, dense_product
    use rundat, only: pub_eigensolver_orfac, pub_eigensolver_abstol, &
         pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_trace
    use utils, only : utils_alloc_check, utils_dealloc_check
    use timer, only: timer_clock


    implicit none

    ! ars: arguments
    type(SPAM3_EMBED),   intent(inout) :: denskern(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: ham(pub_num_spins)
    type(SPAM3_EMBED),   intent(in   ) :: overlap
    integer,       intent(in   ) :: n_occ(max_spins)
    real(kind=DP), optional, intent(  out) :: Ehomo(pub_num_spins)
    real(kind=DP), optional, intent(  out) :: Elumo(pub_num_spins)
    real(kind=DP), optional, intent(  out) :: Egap(pub_num_spins)

    ! ars: local variables
    integer   :: ierr
    integer   :: n_eigs
    integer   :: is
    real(kind=DP), allocatable :: eigenvalues(:)
    type(DEM) :: denskern_dens, overlap_dens, ham_dens, eigs_dens
    logical :: do_gap
    ! agrecocmplx
    logical :: loc_cmplx
    character :: opB_loc

    ! ars: start timer
    call timer_clock('kernel_diis_build_idemp_dkn',1)

    ! agrecocmplx
    loc_cmplx = overlap%p%iscmplx

    ! agrecocmplx
    if (loc_cmplx) then
       opB_loc = 'C'
    else
       opB_loc = 'T'
    end if

    ! ars: check optional arguments
    do_gap=.false.
    if (present(Ehomo).and.present(Elumo).and.present(Egap)) do_gap=.true.


    ! ars: allocate Hamiltonian eigenvalues
    allocate(eigenvalues(num),stat=ierr)
    call utils_alloc_check('kernel_diis_build_idemp_dkn','eigenvalues',ierr)


    ! ndmh: create temporary dense matrices
    ! agrecocmplx: complex dense matrices in complex case
    call dense_create(denskern_dens,num,num,iscmplx=loc_cmplx)
    call dense_create(ham_dens,num,num,iscmplx=loc_cmplx) !ars
    call dense_create(overlap_dens,num,num,iscmplx=loc_cmplx)
    call dense_create(eigs_dens,num,num,iscmplx=loc_cmplx)


    ! cks: loop over spins
    do is=1, pub_num_spins

       ! ars: number eigenvalues to calculate
       n_eigs = n_occ(is) + 1
       if (n_occ(is).eq.num) n_eigs = num
       eigenvalues(:) = 0.0_DP


       ! ars: expand denskern, ham and overlap to dense matrices
       call dense_convert(ham_dens,ham(is))
       call dense_convert(overlap_dens,overlap)


       ! ars: solve HC = eSC
       call dense_eigensolve(n_eigs,eigenvalues,ham_dens,overlap_dens, &
            1,eigs_dens,pub_eigensolver_orfac, pub_eigensolver_abstol)


       ! ndmh: construct a density kernel from the first n_occ(is) eigenvectors
       ! agrecocmplx: need to take hermitian of matrix B in complex case?
       ! I think so....
       call dense_product(denskern_dens,eigs_dens,eigs_dens, &
            opA='N',opB=opB_loc, first_k=1,last_k=n_occ(is))


       ! ndmh: convert back to sparse matrix
       call dense_convert(denskern(is),denskern_dens)


       ! ars: calculate HOMO, LUMO and energy gap
       if (do_gap) then
          Ehomo(is) = eigenvalues(n_occ(is))
          Elumo(is) = eigenvalues(n_occ(is)+1)
          Egap(is) = Elumo(is) - Ehomo(is)
       end if

    end do

    ! ndmh: clean up temporary matrices
    call dense_destroy(overlap_dens)
    call dense_destroy(ham_dens) !ars
    call dense_destroy(denskern_dens)
    call dense_destroy(eigs_dens)

    ! ars: deallocate
    deallocate(eigenvalues,stat=ierr)
    call utils_dealloc_check('kernel_diis_build_idemp_dkn','eigenvalues',ierr)

    ! ars: stop timer
    call timer_clock('kernel_diis_build_idemp_dkn',2)


  end subroutine kernel_diis_build_idemp_dkn



  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !
  !                           END PUBLIC COMMON ROUTINES
  !_____________________________________________________________________________
  !_____________________________________________________________________________


end module kernel_diis
