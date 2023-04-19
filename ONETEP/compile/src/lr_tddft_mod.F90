
! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!         Linear-Response Time-Dependent DFT  module             !
!                                                                !
! This module contains routines associated with the use of the   !
! LR-TDDFT method for calculation of excitations of a system,    !
! and the resulting optical spectra.                             !
!----------------------------------------------------------------!
! This module was created by Tim Zuehlsdorff in 2013.            !
!================================================================!

module lr_tddft

  use constants, only: DP

  implicit none

  private

  public :: lr_tddft_calc_excitations

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_calc_excitations(val_denskern, val_rep, &
       val_ngwf_basis, cond_ngwf_basis, mdl, hfxstate, lhxc_fine, &
       proj_basis, hub_proj_basis, total_energy,&
       hub,nl_projectors,joint_ngwf_basis, &
       aux_rep, aux_ngwf_basis, aux_denskern)

    !===================================================================!
    ! Subroutine to calculate the optical transition energies given the !
    ! input of conduction and valence NGWF's, as well as the valence    !
    ! and conduction density kernels. This routine is the main driver   !
    ! for the LR-TDDFT calculation.                                     !
    ! Modified by Andrea Greco on 28/06/2015 to use complex NGWFs.      !
    ! Modified for embedding by Joseph Prentice, July 2018.             !
    ! Modified to support XC with HFX by James Womack in 2019.          !
    !-------------------------------------------------------------------!
    ! Variables:                                                        !
    !                                                                   !
    !===================================================================!

    use datatypes, only: data_functions_alloc
    use comms, only: comms_barrier, pub_on_root
    use dense, only: DEM, dense_create, dense_destroy
    use function_basis, only: FUNC_BASIS
    use constants, only: stdout, NORMAL
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: DKERN, kernel_create, kernel_destroy
    use model_type, only: MODEL
    use ngwfs, only: ngwfs_initialise
    use ngwf_representation, only: NGWF_REP, ngwf_rep_create, NGWF_HAM
    use projectors, only: PROJECTOR_SET
    use qnto, only: qnto_calculate
    use lr_tddft_rpa, only: lr_tddft_rpa_calculate
    use rundat, only: pub_output_detail, pub_lr_tddft_RPA,&
         pub_lr_tddft_num_states, pub_lr_tddft_joint_set, &
         pub_use_aux_ngwfs, pub_debug_on_root, pub_num_spins, &
         pub_num_kpoints, pub_qnto_analysis, pub_emft, pub_emft_follow, &
         PUB_1K
    use services, only: services_flush
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
         sparse_embed_create, sparse_embed_destroy, sparse_embed_product
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type (FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type (FUNC_BASIS), intent(inout) :: cond_ngwf_basis(:)
    type (MODEL), intent(inout) :: mdl
    type (HFX_STATE), intent(inout), target :: hfxstate
    type (NGWF_REP), intent(inout) :: val_rep
    type (SPAM3_EMBED_ARRAY), intent(inout) :: val_denskern
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,mdl%nsub)
    type (FUNC_BASIS), intent(inout) :: hub_proj_basis(:)
    type (FUNC_BASIS), intent(in) :: proj_basis(:)
    real (kind=DP), intent(inout) :: total_energy
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(:)
    type(FUNC_BASIS), intent(inout) :: joint_ngwf_basis(:)
    type(FUNC_BASIS), optional, intent(in) :: aux_ngwf_basis(:)
    type(NGWF_REP), optional, intent(in) :: aux_rep
    type (SPAM3_EMBED_ARRAY), optional, intent(in) :: aux_denskern

    ! Local Variables
    type (NGWF_REP) :: cond_rep
    type (NGWF_HAM) :: cond_ham
    type (NGWF_HAM) :: val_ham
    type (NGWF_REP) :: joint_rep
    type (NGWF_HAM) :: joint_ham
    type (SPAM3_EMBED), allocatable, dimension(:) :: hvkv, kchc
    type(DEM), allocatable, dimension(:) :: cond_evecs, val_evecs, joint_evecs
    type(SPAM3_EMBED), allocatable, dimension(:,:) :: full_evecs
    type(SPAM3_EMBED), allocatable, dimension(:,:,:) :: full_evecs_RPA
    type (DKERN) :: cond_denskern
    type (DKERN) :: joint_denskern
    real (kind=DP), allocatable, dimension(:,:,:,:,:) :: fxc_fine
    real (kind=DP), allocatable, dimension(:,:,:,:) :: dens_fine
    real (kind=DP), allocatable, dimension(:,:,:,:) :: ground_state_dens
    real (kind=DP), allocatable, dimension(:,:,:,:,:,:) :: sub_dens_fine
    real (kind=DP), allocatable, dimension(:,:,:,:,:,:) :: sub_ground_state_dens
    ! The eigenvalues we want to solve for
    real (kind=DP), allocatable, dimension(:) :: lr_tddft_evals
    integer :: icount
    integer :: ierr,is
    real(kind=DP), allocatable, dimension(:,:) :: val_evals,&
         cond_evals, joint_evals
    ! agrecocmplx
    logical :: loc_cmplx
    ! jcap: embedding variables
    integer :: val_ngwf_num, cond_ngwf_num, isub

    if (pub_debug_on_root) write(stdout, '(a)') &
         'DEBUG: Entering LR_TDDFT'

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine lr_tddft_calc_excitations not ready yet for more&
         & than one k-point.')

    ! tjz21: Lock out the possibility of spin-polarised calculations
    ! for the time being as this feature is not fuly tested yet
    !call utils_assert(pub_num_spins==1, &
    !     'Spin-polarised calculations not yet enabled for LR_TDDFT')

    loc_cmplx = val_denskern%m(1,1)%p%iscmplx

    ! tjz07: Allocate storage space for fxc.
    allocate(fxc_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,2,mdl%nsub), stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations','fxc_fine',ierr)


    ! tjz07: Trying to allocate conduction NGWFs
    ! agrecocmplx: allocate using appropriate routine
    ! jcap: loop over regions
    ! jcap: also find total number of NGWFs
    val_ngwf_num=0
    cond_ngwf_num=0
    allocate(cond_rep%ngwfs_on_grid(mdl%nsub),stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations','cond_rep%ngwfs_on_grid',ierr)
    do isub=1,mdl%nsub
       call data_functions_alloc(cond_rep%ngwfs_on_grid(isub), &
            cond_ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
       val_ngwf_num=val_ngwf_num+val_ngwf_basis(isub)%num
       cond_ngwf_num=cond_ngwf_num+cond_ngwf_basis(isub)%num
       call comms_barrier
    end do

    ! tjz07: Creating a Conduction representation for NGWFs
    call ngwf_rep_create(cond_rep, 'c', mdl, is_cmplx=loc_cmplx)

    ! tjz07: Reading converged conduction NGWFs from file
    call ngwfs_initialise(cond_rep,cond_ngwf_basis,mdl,hfxstate, &
         do_not_reexpand = .true.)
    ! JCW: do_not_reexpand == .true. because
    !      (i)  cond_rep%swexes have not been initialised yet (this is done by
    !           hf_exchange_init in lr_tddft_initialise
    !      (ii) ngwfs_initialise would initialise the REP_SWEX_HFX element
    !           of cond_rep%swexes (cond-cond NGWF product expansion),
    !      For evaluating the conduction NGWF Hamiltonian, we need cond-val
    !      NGWF expansions, to be placed in the REP_SWEX_HFX_OTHER element of
    !      cond_rep%swexes by hf_exchange_dkn_indep_stage

    call comms_barrier
    if (pub_on_root .and. (pub_output_detail>=NORMAL)) &
         write(stdout, '(a)') '...Done'
    call services_flush

    ! Creating structures for conduction denskern and joint denskern
    call kernel_create(cond_denskern, 'K'//cond_rep%postfix, &
         is_cmplx=loc_cmplx)

    ! tjz07: Allocate storage space for a density on a fine grid.
    ! tjz07: Allocate storage space for fine denity
    allocate(dens_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations', &
         'dens_fine', ierr)
    allocate(ground_state_dens(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations', &
         'ground_state_dens', ierr)
    allocate(sub_dens_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins,mdl%nsub,mdl%nsub), stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations', &
         'sub_dens_fine', ierr)
    allocate(sub_ground_state_dens(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins,mdl%nsub,mdl%nsub), stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations', &
         'sub_ground_state_dens', ierr)


    ! allocate evals for analysis of TDDFT transitions.
    allocate(val_evecs(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations','val_evecs',ierr)
    allocate(cond_evecs(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations','cond_evecs',ierr)
    allocate(joint_evecs(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations','joint_evecs',ierr)
    do is=1, pub_num_spins
       call dense_create(val_evecs(is), val_ngwf_num, &
            val_ngwf_num, iscmplx=loc_cmplx)
       call dense_create(cond_evecs(is), cond_ngwf_num, &
            cond_ngwf_num, iscmplx=loc_cmplx)
    enddo
    ! allocate evals
    allocate(val_evals(val_ngwf_num,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations','val_evals',&
         ierr)
    allocate(cond_evals(cond_ngwf_num,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations','val_evals',&
         ierr)
    allocate(joint_evals(val_ngwf_num+&
         cond_ngwf_num,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations','val_evals',&
         ierr)
    ! tjz07: Function constructs the ground state and conduction state
    ! density kernels and hamiltonians to be used in the remainder of
    ! the TDDFT code.
    ! function also calculates the TDDFT exchange correlation kernel
    call lr_tddft_initialise(total_energy, val_denskern, val_rep, &
         val_ngwf_basis, proj_basis, hub_proj_basis, cond_ngwf_basis, &
         cond_denskern,cond_rep,mdl,hfxstate,lhxc_fine, cond_ham, val_ham, &
         ground_state_dens, sub_ground_state_dens, &
         fxc_fine,hub, nl_projectors, val_evecs, cond_evecs,&
         joint_ngwf_basis, joint_rep, joint_ham, joint_denskern,&
         joint_evecs,val_evals,cond_evals,joint_evals)

    call comms_barrier

    ! create storage space for lr_tddft_evals and evecs
    allocate(lr_tddft_evals(pub_lr_tddft_num_states), stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations', &
         'lr_tddft_evals',ierr)
    allocate(full_evecs(pub_lr_tddft_num_states,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations','full_evecs',&
         ierr)
    if(pub_lr_tddft_RPA) then
       allocate(full_evecs_RPA(pub_lr_tddft_num_states,2,pub_num_spins),stat=ierr)
       call utils_alloc_check('lr_tddft_calc_excitations','full_evecs_RPA',&
          ierr)
    endif

    ! create correct sparse matrix structure for the TDDFT eigenvectors
    do icount=1, pub_lr_tddft_num_states
       do is=1, pub_num_spins
          if(.not. pub_lr_tddft_RPA) then
             if(pub_lr_tddft_joint_set) then
                full_evecs(icount,is)%structure= 'TDR'//trim(joint_rep%postfix)
             else
                full_evecs(icount,is)%structure= 'TDR'//trim(cond_rep%postfix)
             endif
             call sparse_embed_create(full_evecs(icount,is),iscmplx=loc_cmplx)
             ! allocate full RPA evecs
          else
             if(pub_lr_tddft_joint_set) then
                full_evecs_RPA(icount,1,is)%structure= 'TDR'//trim(joint_rep%postfix)
                full_evecs_RPA(icount,2,is)%structure= 'TDR'//trim(joint_rep%postfix)
             else
                full_evecs_RPA(icount,1,is)%structure= 'TDR'//trim(cond_rep%postfix)
                full_evecs_RPA(icount,2,is)%structure= 'TDR'//trim(cond_rep%postfix)
             endif
             call sparse_embed_create(full_evecs_RPA(icount,1,is),iscmplx=loc_cmplx)
             call sparse_embed_create(full_evecs_RPA(icount,2,is),iscmplx=loc_cmplx)
          endif
       enddo ! end loop over spins
    enddo

    call timer_clock('lr_tddft_calc_excitations', 1)

    ! create kchc and hvkv, the valence and conduction hamiltonian times
    ! their appropriate density kernels
    allocate(kchc(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations','kchc',ierr)
    allocate(hvkv(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations','kchc',ierr)
    do is=1, pub_num_spins
       if (pub_lr_tddft_joint_set) then
          call sparse_embed_create(kchc(is),joint_denskern%kern%m(is,PUB_1K),&
               joint_ham%ham(is))
          call sparse_embed_create(hvkv(is),val_ham%ham(is),&
               val_denskern%m(is,PUB_1K))
          call sparse_embed_product(kchc(is),joint_denskern%kern%m(is,PUB_1K),&
               joint_ham%ham(is))
          call sparse_embed_product(hvkv(is),val_ham%ham(is),&
               val_denskern%m(is,PUB_1K))
       else
          call sparse_embed_create(kchc(is),cond_denskern%kern%m(is,PUB_1K),&
               cond_ham%ham(is))
          call sparse_embed_create(hvkv(is),val_ham%ham(is),&
               val_denskern%m(is,PUB_1K))
          call sparse_embed_product(kchc(is),cond_denskern%kern%m(is,PUB_1K),&
               cond_ham%ham(is))
          call sparse_embed_product(hvkv(is),val_ham%ham(is),&
               val_denskern%m(is,PUB_1K))
       endif
    enddo

    ! call main lr_TDDFT function. Check if we use the joint set
    ! to represent the conduction space or the conduction set.
    if(.not. pub_lr_tddft_RPA) then
       if (pub_lr_tddft_joint_set) then
          call lr_tddft_calculate(full_evecs,lr_tddft_evals, dens_fine, &
               sub_dens_fine, fxc_fine, joint_ngwf_basis, joint_rep, &
               joint_denskern%kern%m(:,PUB_1K), &
               kchc, joint_ham%ham(:), val_ngwf_num+cond_ngwf_num, &
               val_ngwf_basis, val_rep, val_denskern%m(:,PUB_1K), &
               hvkv, val_ham%ham(:), val_ngwf_num, &
               proj_basis, nl_projectors, mdl, hfxstate, &
               val_evecs, joint_evecs, ground_state_dens, &
               sub_ground_state_dens, joint_evals)
       else
          call lr_tddft_calculate(full_evecs,lr_tddft_evals, dens_fine, &
               sub_dens_fine, fxc_fine, cond_ngwf_basis, cond_rep, &
               cond_denskern%kern%m(:,PUB_1K), &
               kchc, cond_ham%ham(:), cond_ngwf_num, &
               val_ngwf_basis, val_rep, val_denskern%m(:,PUB_1K), &
               hvkv, val_ham%ham(:), val_ngwf_num, &
               proj_basis, nl_projectors, mdl, hfxstate, &
               val_evecs, cond_evecs, ground_state_dens, &
               sub_ground_state_dens)
       endif
       if (pub_qnto_analysis) then
          call utils_assert(mdl%nsub.eq.1,'QNTO not ready yet for more&
               & than one subsystem.')
          call qnto_calculate(full_evecs(:,1), lr_tddft_evals, val_rep, &
               val_ngwf_basis(1), val_denskern, cond_rep, cond_ngwf_basis(1), &
               cond_denskern%kern, joint_rep, joint_ngwf_basis(1), &
               joint_denskern%kern, pub_lr_tddft_num_states,&
               mdl, aux_rep, aux_ngwf_basis(1), aux_denskern)
       end if
    else
       if (pub_lr_tddft_joint_set) then
          call lr_tddft_rpa_calculate(full_evecs_RPA,lr_tddft_evals, dens_fine, &
               sub_dens_fine, fxc_fine, joint_ngwf_basis, joint_rep, &
               joint_denskern%kern%m(:,PUB_1K), &
               kchc, joint_ham%ham(:), val_ngwf_num+cond_ngwf_num, &
               val_ngwf_basis, val_rep, val_denskern%m(:,PUB_1K), &
               hvkv, val_ham%ham(:), val_ngwf_num, &
               proj_basis, nl_projectors, mdl, hfxstate, &
               val_evecs, joint_evecs, ground_state_dens, &
               sub_ground_state_dens, joint_evals)
       else
          call lr_tddft_rpa_calculate(full_evecs_RPA,lr_tddft_evals, dens_fine, &
               sub_dens_fine, fxc_fine, cond_ngwf_basis, cond_rep, &
               cond_denskern%kern%m(:,PUB_1K), &
               kchc, cond_ham%ham(:), cond_ngwf_num, &
               val_ngwf_basis, val_rep, val_denskern%m(:,PUB_1K), &
               hvkv, val_ham%ham(:), val_ngwf_num, &
               proj_basis, nl_projectors, mdl, hfxstate, &
               val_evecs,cond_evecs, ground_state_dens, &
               sub_ground_state_dens)
       endif
       if (pub_qnto_analysis) then
          call utils_assert(mdl%nsub.eq.1,'QNTO not ready yet for more&
               & than one subsystem.')
          call qnto_calculate(full_evecs_RPA(:,2,1), lr_tddft_evals, val_rep, &
               val_ngwf_basis(1), val_denskern, cond_rep, cond_ngwf_basis(1), &
               cond_denskern%kern, joint_rep, joint_ngwf_basis(1), &
               joint_denskern%kern, pub_lr_tddft_num_states,&
               mdl, aux_rep, aux_ngwf_basis(1), aux_denskern)
       end if
    endif

    call timer_clock('lr_tddft_calc_excitations', 2)

    ! Clean up: Deallocate data structures
    do icount=1, pub_lr_tddft_num_states
       do is=1, pub_num_spins
          if(.not. pub_lr_tddft_RPA) then
             call sparse_embed_destroy(full_evecs(icount,is))
          else
             call sparse_embed_destroy(full_evecs_RPA(icount,1,is))
             call sparse_embed_destroy(full_evecs_RPA(icount,2,is))
          endif
       enddo
    enddo
    deallocate(dens_fine, stat=ierr)
    call utils_dealloc_check('lr_tddft_calc_excitations', 'dens_fine',&
         ierr)
    deallocate(sub_dens_fine, stat=ierr)
    call utils_dealloc_check('lr_tddft_calc_excitations', 'sub_dens_fine',&
         ierr)
    deallocate(fxc_fine, stat=ierr)
    call utils_dealloc_check('lr_tddft_calc_excitations', 'fxc_fine',&
         ierr)
    deallocate(lr_tddft_evals, stat=ierr)
    call utils_dealloc_check('lr_tddft_calc_excitations', &
         'lr_tddft_evals', ierr)
    deallocate(full_evecs, stat=ierr)
    call utils_dealloc_check('lr_tddft_calc_excitations','full_evecs',&
         ierr)
    if(pub_lr_tddft_RPA) then
       deallocate(full_evecs_RPA, stat=ierr)
       call utils_dealloc_check('lr_tddft_calc_excitations','full_evecs_RPA',&
            ierr)
    endif
    call kernel_destroy(cond_denskern)
    do is=1, pub_num_spins
       call sparse_embed_destroy(kchc(is))
       call sparse_embed_destroy(hvkv(is))
    enddo
    deallocate(kchc, stat=ierr)
    call utils_dealloc_check('lr_tddft_calc_excitations','kchc',ierr)
    deallocate(hvkv, stat=ierr)
    call utils_dealloc_check('lr_tddft_calc_excitations','kchc',ierr)
    do is=1, pub_num_spins
       call dense_destroy(val_evecs(is))
       call dense_destroy(cond_evecs(is))
       if(pub_lr_tddft_joint_set) then
          call dense_destroy(joint_evecs(is))
       endif
    enddo
    deallocate(val_evecs, stat=ierr)
    call utils_dealloc_check('lr_tddft_calc_excitations','val_evecs',ierr)
    deallocate(cond_evecs, stat=ierr)
    call utils_alloc_check('lr_tddft_calc_excitations','cond_evecs',ierr)
    deallocate(joint_evecs, stat=ierr)
    call utils_dealloc_check('lr_tddft_calc_excitations','joint_evecs',ierr)
    deallocate(val_evals, stat=ierr)
    call utils_dealloc_check('lr_tddft_calc_excitations','val_evals',&
         ierr)
    deallocate(cond_evals, stat=ierr)
    call utils_dealloc_check('lr_tddft_calc_excitations','val_evals',&
         ierr)
    deallocate(joint_evals, stat=ierr)
    call utils_dealloc_check('lr_tddft_calc_excitations','val_evals',&
         ierr)
    deallocate(ground_state_dens, stat=ierr)
    call utils_dealloc_check('lr_tddft_calc_excitations', &
         'ground_state_dens', ierr)
    deallocate(sub_ground_state_dens, stat=ierr)
    call utils_dealloc_check('lr_tddft_calc_excitations', &
         'sub_ground_state_dens', ierr)


  end subroutine lr_tddft_calc_excitations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_initialise(total_energy, val_denskern, val_rep, &
       val_ngwf_basis, proj_basis, hub_proj_basis, cond_ngwf_basis, &
       cond_denskern, cond_rep, mdl, hfxstate, lhxc_fine, cond_ham, val_ham, &
       dens_fine, sub_dens_fine, fxc_fine, hub, nl_projectors, val_evecs, &
       cond_evecs, joint_ngwf_basis, joint_rep, joint_ham, joint_denskern,&
       joint_evecs, val_evals, cond_evals, joint_evals)

    !=========================================================================!
    ! Subroutine setting up all appropriate matrices needed for the LR_TDDFT  !
    ! Calculation, such as the conduction and valence Hamiltonian, the        !
    ! projector onto the unoccupied subspace, and the Kohn-Sham matrices of   !
    ! eigenvectors needed for analysing TDDFT transitions.                    !
    ! Modified for embedding by Joseph Prentice, July 2018                    !
    ! Modified to support XC with HFX by James Womack in 2019.                !
    !=========================================================================!

    use augmentation, only: augmentation_density_on_grid
    use datatypes, only: data_functions_alloc
    use comms, only: comms_barrier, pub_on_root
    use constants, only: stdout, paw_en_size, UP, DN, VERBOSE, &
         REP_SWEX_HFX_OTHER
    use dense, only: DEM, dense_create, dense_destroy, &
         dense_convert, dense_eigensolve
    use density, only: density_on_grid
    use electronic_init, only: electronic_init_denskern
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_dens_indep_matrices, &
         hamiltonian_dens_dep_matrices, hamiltonian_build_matrix, &
         hamiltonian_dens_dep_nonsc
    use hf_exchange, only: HFX_STATE, hf_exchange_init, &
         hf_exchange_dkn_indep_stage
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: DKERN, kernel_create
    use model_type, only: MODEL
    use ngwfs, only: ngwfs_merge_sets
    use ngwf_representation, only: NGWF_REP, ngwf_rep_create,&
         NGWF_HAM, ngwf_ham_create
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_output_detail, cond_num_states, pub_spin, &
         pub_aug, pub_debug_on_root, cond_init_shift, pub_lr_tddft_projector, &
         pub_lr_tddft_analysis, pub_num_spins, pub_nhat_in_xc, pub_aug_den_dim,&
         pub_nlcc, pub_spin_fac, pub_paw, pub_num_kpoints, PUB_1K,&
         pub_lr_tddft_xc_finite_diff, &
         pub_xc_ke_density_required, pub_emft, pub_active_region, pub_devel_code, &
         pub_use_hfx, pub_use_activehfx, cond_read_denskern
    use services, only: services_flush
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_copy, sparse_embed_axpy, &
         sparse_embed_product,sparse_embed_transpose_structure, &
         sparse_embed_transpose, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_destroy
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_devel_code
    use xc, only: xc_fxc_potential, xc_gradients, xc_fxc_potential_emft

    implicit none

    ! Arguments
    type (FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type (FUNC_BASIS), intent(inout) :: cond_ngwf_basis(:)
    type (MODEL), intent(inout) :: mdl
    type (HFX_STATE), intent(inout), target :: hfxstate
    type (NGWF_REP), intent(in) :: val_rep
    type (NGWF_REP), intent(inout) :: cond_rep
    type (SPAM3_EMBED_ARRAY), intent(inout) :: val_denskern
    type (DKERN), intent(inout) :: cond_denskern
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,mdl%nsub)
    type (FUNC_BASIS), intent(inout) :: hub_proj_basis(:)
    type (FUNC_BASIS), intent(in) :: proj_basis(:)
    real (kind=DP), intent(inout) :: total_energy
    type (NGWF_HAM), intent(inout) :: cond_ham
    type (NGWF_HAM), intent(inout) :: val_ham
    real (kind=DP), intent(inout) :: sub_dens_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2,mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)
    real (kind=DP), intent(inout) :: dens_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2,mdl%fine_grid%max_slabs12,pub_num_spins)
    real (kind=DP), intent(inout) :: fxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2,mdl%fine_grid%max_slabs12,2,mdl%nsub)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(:)
    type(DEM), intent(inout) :: val_evecs(pub_num_spins)
    type(DEM), intent(inout) :: cond_evecs(pub_num_spins)
    type (NGWF_HAM), intent(inout) :: joint_ham
    type (FUNC_BASIS), intent(inout) :: joint_ngwf_basis(:)
    type (NGWF_REP), intent(inout) :: joint_rep
    type (DKERN), intent(inout) :: joint_denskern
    type(DEM), intent(inout) :: joint_evecs(pub_num_spins)
    real(kind=DP), intent(inout) :: val_evals(:,:)
    real(kind=DP), intent(inout) :: cond_evals(:,:)
    real(kind=DP), intent(inout) :: joint_evals(:,:)

    ! Local Variables
    real (kind=DP) :: paw_sphere_energies(paw_en_size)
    real (kind=DP) :: lhxc_energy
    real (kind=DP) :: hubbard_energy
    real (kind=DP), allocatable, dimension(:,:,:,:,:) :: nhat_den_grad
    integer :: ierr
    integer :: is
    logical :: shift_changed
    integer :: fine_max_slabs12, fine_ld1, fine_ld2 ! jd: shortcut notation

    ! jcap: embedding local variables
    type(SPAM3) :: kern_array(pub_num_spins)
    integer :: isub,jsub,val_ngwf_num,cond_ngwf_num,joint_ngwf_num
    real (kind=DP), allocatable :: sub_nhat_den_grad(:,:,:,:,:,:,:)
    ! jcap: specifically EMFT variables
    real (kind=DP), allocatable :: full_fxc_fine(:,:,:,:)
    real (kind=DP), allocatable :: fxc_fine_emft(:,:,:,:)
    real (kind=DP), allocatable :: dens_buffer(:,:,:,:)
    real (kind=DP), allocatable :: full_dens_grad(:,:,:,:,:)
    real (kind=DP), allocatable :: emft_dens_grad(:,:,:,:,:)
    real (kind=DP), allocatable :: demft_dfull(:,:,:,:)
    complex(kind=DP), allocatable :: recip_work(:,:,:,:)
    real (kind=DP) :: grad1(3), grad2(3), grad1_abs, grad2_abs
    integer :: i1,i2,islab12
    logical :: temp_cond_read_denskern

    ! local variables for diagnostics+generating new denskerns
    type(DEM) :: overlap_dense, ham_dense
    type (SPAM3_EMBED) :: eff_cond_denskern, S_crossS_inv
    type (SPAM3_EMBED) :: S_crossKv, S_invS_crossKv, S_cross_trans
    type (SPAM3_EMBED) :: eff_cond_denskern2

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine lr_tddft_initialise not ready yet for more&
         & than one k-point.')

    call sparse_embed_create(eff_cond_denskern, cond_denskern%kern%m(1,1))
    call sparse_embed_create(S_crossS_inv,cond_rep%cross_overlap,&
         cond_rep%inv_overlap)
    call sparse_embed_transpose_structure(S_cross_trans%structure,&
         cond_rep%cross_overlap)
    call sparse_embed_create(S_cross_trans, iscmplx=cond_rep%cross_overlap%p%iscmplx, &
         mrows=mdl%nsub, ncols=mdl%nsub)
    call sparse_embed_create(S_crossKv, S_cross_trans, val_denskern%m(1,1))
    call sparse_embed_create(S_invS_crossKv, cond_rep%inv_overlap, &
         S_crossKv)
    call sparse_embed_create(eff_cond_denskern2, S_invS_crossKv,&
         S_crossS_inv)

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout, '(a)') &
         'DEBUG: Entering lr_tddft_initialise routine.'

    ! JCW: Abort if KE density required (meta-GGA + LR-TDDFT not
    ! JCW: implemented).
    call utils_assert(.not.pub_xc_ke_density_required, "Error in &
         &lr_tddft_initialise: pub_xc_ke_density_required is true, but &
         &combination of KE-density-dependent XC functionals with LR-TDDFT &
         &has not been implemented/tested.")

    ! tjz07: Try to create the valence hamiltonian
    call ngwf_ham_create(val_ham, val_rep)

    ! JCW: Presumably a call to hamiltonian_dens_indep_matrices is not necessary
    !      here because the density-independent parts of the Hamiltonian (e.g.
    !      val_rep%kinet) are already part of val_rep (passed to
    !      lr_tddft_calc_excitations)

    !ndmh: calculate density dependent energies and matrices for valence NGWFS
    call hamiltonian_dens_dep_matrices(val_ham, lhxc_fine, total_energy, &
         lhxc_energy, hubbard_energy, paw_sphere_energies, val_rep, &
         val_ngwf_basis, hub_proj_basis, hub, val_denskern, &
         mdl, hfxstate, ham_update=.true., lhxc_fixed=.false.)

    call hamiltonian_build_matrix(val_ham, val_rep)

    if (pub_on_root .and. pub_output_detail>=VERBOSE) then
       write(stdout, '(a)') 'Successfully built hamiltonian matrix &
            &for valence ham.'
    endif

    ! ndmh: create cond hamiltonian
    call ngwf_ham_create(cond_ham, cond_rep)

    ! JCW: Initialise HFx for conduction (only in the active region for EMFT)
    if (pub_use_hfx.or.pub_use_activehfx) then
       ! JCW: For EMFT calculations, HFX is only available in the active region
       !      For non-EMFT calculations, pub_active_region == 1 by default
       isub = pub_active_region
       call hf_exchange_init(hfxstate, cond_rep, &
            cond_ham%hfexchange(1)%m(isub,isub), mdl%cell, &
            size(mdl%regions(isub)%elements), init_rep_only = .true.)

       ! JCW: We need to construct the cond NGWF Hamiltonian, which
       !      requires mixed (cond-val) NGWF expansions in SWs. The expansion
       !      coefficients are not computed in ngwfs_initialise so need
       !      to be computed manually, as in conduction_ngwf_optimise. The
       !      SW_EX representing the expansion coeffs is stored in
       !      cond_rep%swexes(REP_SWEX_HFX_OTHER)
       call hf_exchange_dkn_indep_stage(hfxstate, mdl, isub, &
            mdl%regions(isub)%par, cond_rep, cond_ngwf_basis(isub), &
            rep2 = val_rep, ngwf_basis2 = val_ngwf_basis(isub), &
            basis_selector = [ 1,1,2,2,REP_SWEX_HFX_OTHER ])
            ! basis for alpha: 1 (cond)
            ! basis for beta: 1 (cond)
            ! basis for gamma: 2 (val)
            ! basis for delta: 2 (val)
            ! swex selector: HFX_OTHER (of cond_rep)
            ! (corresponds to ket_basis_selector = (/1,2,REP_SWEX_HFX_OTHER/)
            ! which is passed to swx_expand_my_ngwf_pairs inside
            ! hf_exchange_dkn_indep_stage)
            !
            ! Hamiltonian matrix elements are <alpha|H|beta> and the exchange
            ! contribution is <alpha gamma|1/r12|delta beta> P^{delta gamma}
            ! in physicists' notation (electron indices <12|12>), with P a
            ! suitable density matrix to contract with NGWFs gamma and delta.
            !
            ! For the conduction Hamiltonian alpha, beta are conduction
            ! NGWFs and gamma, delta (from the exchange operator) are valence
            ! NGWFs and thus P is the valence density kernel K

       ! JCW: Note that we should not now call ngwfs_initialise or any other
       !      routine that results in a call to ngwf_rep_register_change for
       !      cond_rep until the conduction Hamiltonian is evaluated.
       !      ngwf_rep_register_changes calls swx_invalidate_expansion for
       !      each element of cond_rep%swexes(:), purging the coeffs_ht
       !      hash table for each SW_EX.
    end if

    ! set conduction occupation numbers
    cond_rep%n_occ(UP,PUB_1K) = (cond_num_states - pub_spin)/2
    if (pub_num_spins > 1) then
       cond_rep%n_occ(DN,PUB_1K) = (cond_num_states + pub_spin)/2
    end if

    ! tjz07: Initialise conduction density kernel. Should read from file
    ! tjz07: and calculate overlap inverse.
    ! --------------------------------------------------------------------------
    ! JCW: 26/03/19 We only expect cond_denskern to be read from a file here
    !      When cond_read_denskern == .true. (which is not the default).
    !
    !      The default settings redefine the conduction density kernel as a
    !      projector based on the inverse overlap of the conduction NGWFs, so
    !      the conduction density kernel is not needed later
    ! --------------------------------------------------------------------------
    ! JCW: 26/03/19 @optimize
    !      Avoid pointless initialization of initial guess cond_denskern:
    !      When using the projector with cond_read_densken == .false.,
    !      this call seems to pointlessly initialize an initial guess
    !      cond_denskern, which is not used later. But, electronic_init_denskern
    !      also appears to do a few other things (e.g. call
    !      hamiltonian_dens_indep_matrices to initialize density kernel
    !      independent Hamiltonian matrix components in cond_rep), so it likely
    !      cannot be simply skipped when cond_read_denskern == .false.
    !
    !      The pointless creation of an initial guess cond density
    !      kernel also appears to have the side effect that any
    !      existing .dkn_cond is overwritten by the initial guess
    !      (which is probably not desirable).
    !      --------------------------------------------------------------------------
    !      jcap: ensure that the conduction kernel is read in here if
    !      we are doing an EMFT calculation (over-ride the user
    !      setting of cond_read_denskern), to ensure that the
    !      EMFT-optimised kernel is read in
    if (pub_emft) then
       temp_cond_read_denskern = cond_read_denskern
       cond_read_denskern = .true.
    end if
    cond_ham%cond_shift = cond_init_shift
    call electronic_init_denskern(cond_rep,cond_ham,cond_denskern, &
         lhxc_fine,cond_ngwf_basis, &
         proj_basis, nl_projectors,hub_proj_basis, hub, mdl, hfxstate, &
         val_rep, val_ngwf_basis, val_denskern, val_ham)
    if (pub_emft) then
       cond_read_denskern = temp_cond_read_denskern
    end if
    call comms_barrier

    ! ndmh: copy in dijhat from valence hamiltonian
    if (pub_aug) then
       do is=1,pub_num_spins
          call sparse_embed_copy(cond_ham%dijhat(is),val_ham%dijhat(is))
       end do
    end if

    ! JCW: Presumably a call to hamiltonian_dens_indep_matrices is not necessary
    !      here because the density-independent parts of the Hamiltonian (e.g.
    !      cond_rep%kinet) are already part of cond_rep (evaluated via a call
    !      hamiltonian_dens_indep_matrices in electronic_init_denskern)

    ! ndmh: just calculate the cond hamiltonian for properties
    ! calculation
    shift_changed = .false. ! <-- TODO this looks like a mistake
                            !     If this argument is present, then
                            !     the projected conduction Hamiltonian
                            !     will be computed regardless of the value
                            !     of shift_changed (hamiltonian_dens_dep_nonsc
                            !     only checks for the presence of this argument)
    call hamiltonian_dens_dep_nonsc(cond_ham, cond_rep,cond_ngwf_basis, &
         mdl,hfxstate,lhxc_fine,hub,val_rep,val_ham, &
         val_denskern%m(:,PUB_1K), updated_shift=shift_changed, &
         val_ngwf_basis = val_ngwf_basis, &         ! optional arguments for HFx
         cond_dkn = cond_denskern%kern%m(:,PUB_1K)) ! optional arguments for HFx
    ! JCW: ^--- TODO It seems odd to pass an uninitialized cond density kernel
    !           here (when cond_read_denskern F) but hf_exchange_calculate
    !           requires it. It is used to evaluate the HFX energy and NGWF
    !           gradient, so probably not needed in this context (where we
    !           only require the X matrix)
    ! JCW: Note that cond_denskern%kern is a SPAM3_EMBED_ARRAY object and
    !      the cond_denskern%kern%m component is an array of SPAM3_EMBED
    !      structures indexed m(is,ik), where "is" is the spin index and "ik"
    !      is the kpoint index. PUB_1K=1 in rundat, presumably for single
    !      k-point/gamma-point case.

    ! JCW: hamiltonian_dens_dep_nonsc includes a call to
    !      hamiltonian_build_matrix, so at this point all cond_ham components
    !      have been summed into cond_ham%ham

    if (pub_on_root .and. pub_output_detail>=VERBOSE) then
       write(stdout, '(a)') 'Successfully created conduction &
            &hamiltonian'
    endif

    ! tjz07: Calculate density from converged valence kernel
    ! in case of PAW, need to add a compensation density
    if(pub_paw) then

       ! jcap: need to create regional density arrays
       fine_ld1 = mdl%fine_grid%ld1
       fine_ld2 = mdl%fine_grid%ld2
       fine_max_slabs12 = mdl%fine_grid%max_slabs12
       allocate(sub_nhat_den_grad(fine_ld1,fine_ld2,fine_max_slabs12,pub_num_spins,&
            0:pub_aug_den_dim,mdl%nsub,mdl%nsub),stat=ierr)
       call utils_alloc_check('lr_tddft_initialise','sub_nhat_den_grad',ierr)

       allocate(nhat_den_grad(fine_ld1,fine_ld2,fine_max_slabs12,pub_num_spins,&
            0:pub_aug_den_dim),stat=ierr)
       call utils_alloc_check('lr_tddft_initialise','nhat_den_grad',ierr)

       sub_nhat_den_grad = 0.0_DP
       nhat_den_grad = 0.0_DP
    end if

    ! jcap: need to loop over regions and sum up density across them
    dens_fine = 0.0_DP
    do jsub=1,mdl%nsub
       do isub=1,mdl%nsub

          ! rc2013: HACK! problems with deferred shape arrays using
          ! embedding structures.  In practice this will probably
          ! require denskern to be added properly to density_on_grid
          ! (complex quantities).
          call sparse_embed_extract_from_array(kern_array,val_denskern%m(:,PUB_1K),&
               isub,jsub)

          call density_on_grid(sub_dens_fine(:,:,:,:,isub,jsub), mdl%fine_grid, &
               mdl%dbl_grid, mdl%cell, mdl%fftbox, kern_array, &
               val_rep%ngwf_overlap%m(isub,jsub), &
               val_rep%ngwfs_on_grid(isub), val_ngwf_basis(isub), &
               val_rep%ngwfs_on_grid(jsub), val_ngwf_basis(jsub))

          dens_fine = dens_fine + sub_dens_fine(:,:,:,:,isub,jsub)

          ! in case of PAW, need to add a compensation density
          if (pub_paw) then
             call augmentation_density_on_grid(sub_nhat_den_grad(:,:,:,:,:,isub,jsub), &
                  mdl%fine_grid,mdl%cell,mdl%regions(isub)%pseudo_sp,&
                  mdl%regions(isub)%paw_sp,mdl%aug_box, &
                  kern_array, val_rep%sp_overlap%m(isub,jsub))

             nhat_den_grad = nhat_den_grad + sub_nhat_den_grad(:,:,:,:,:,isub,jsub)
             dens_fine = dens_fine + sub_nhat_den_grad(:,:,:,:,0,isub,jsub)

          endif

          call sparse_embed_destroy_extracted_array(kern_array)

       end do
    end do

    ! add on core density before calculating fxc potential for NLCC
    if (pub_nlcc) then
       do is=1,pub_num_spins
          dens_fine(:,:,:,is) = dens_fine(:,:,:,is) + &
               mdl%core_density_fine * 0.5_DP * pub_spin_fac
       end do
    end if

    ! tjz21: subtract nhat density if it is not required in fxc.
    if (pub_paw .and. .not. pub_nhat_in_xc) then
       dens_fine(:,:,:,:) = dens_fine(:,:,:,:) - nhat_den_grad(:,:,:,:,0)
    end if

    ! calculate the exchange correlation kernel
    ! only do this if fxc term is not evaluated by finite difference
    fxc_fine = 0.0_DP
    if(pub_num_spins==1 .and. .not. pub_lr_tddft_xc_finite_diff) then

       allocate(full_fxc_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
            mdl%fine_grid%max_slabs12,2),stat=ierr)
       call utils_alloc_check('lr_tddft_initialise','full_fxc_fine',ierr)

       call xc_fxc_potential(mdl%fine_grid, full_fxc_fine, dens_fine)

       ! jcap: EMFT correction - taken from
       ! hamiltonian_lhxc_calculate_embed
       if (pub_emft) then

          allocate(full_dens_grad(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
               mdl%fine_grid%max_slabs12,3,pub_num_spins),stat=ierr)
          call utils_alloc_check('lr_tddft_initialise','full_dens_grad',ierr)
          allocate(emft_dens_grad(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
               mdl%fine_grid%max_slabs12,3,pub_num_spins),stat=ierr)
          call utils_alloc_check('lr_tddft_initialise','emft_dens_grad',ierr)
          allocate(demft_dfull(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
               mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
          call utils_alloc_check('lr_tddft_initialise','demft_dfull',ierr)
          allocate(recip_work(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
               mdl%fine_grid%max_slabs12,3),stat=ierr)
          call utils_alloc_check('lr_tddft_initialise','recip_work',ierr)
          allocate(dens_buffer(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
               mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
          call utils_alloc_check('lr_tddft_initialise','dens_buffer',ierr)
          allocate(fxc_fine_emft(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
               mdl%fine_grid%max_slabs12,2),stat=ierr)
          call utils_alloc_check('lr_tddft_initialise','fxc_fine_emft',ierr)

          ! rc2013: calculate the gradient of the full density if required
          if (utils_devel_code(.false., 'HAM','EMFT_POT',pub_devel_code)) then
             call xc_gradients(dens_fine, full_dens_grad, recip_work, &
                  mdl%fine_grid, 0)
          end if

          ! jcap: If necessary, add on the nlcc charges
          dens_buffer = sub_dens_fine(:,:,:,:,pub_active_region,pub_active_region)

          if (pub_nlcc) then
             do is=1,pub_num_spins
                dens_buffer(:,:,:,is) = dens_buffer(:,:,:,is) + &
                     mdl%regions(pub_active_region)%core_density_fine * 0.5_DP * pub_spin_fac
             end do
          end if

          call xc_fxc_potential_emft(mdl%fine_grid, fxc_fine_emft, dens_buffer)

          ! rc2013: In principle we may need to correct the XC potential
          ! by accounting for the derivative of the total system density
          ! wrt the active density
          if (utils_devel_code(.false., 'HAM','EMFT_POT',pub_devel_code)) then
             call xc_gradients(dens_buffer, emft_dens_grad, recip_work, &
                  mdl%fine_grid, 0)
             demft_dfull = 0.0_DP
             do is=1,pub_num_spins
                do islab12=1,mdl%fine_grid%num_my_slabs12
                   do i2=1,mdl%fine_grid%n2
                      do i1=1,mdl%fine_grid%n1
                         grad1(:) = full_dens_grad(i1,i2,islab12,:,is)
                         grad2(:) = emft_dens_grad(i1,i2,islab12,:,is)
                         grad1_abs = sqrt(grad1(1)*grad1(1)+grad1(2)*grad1(2)+grad1(3)*grad1(3))
                         grad2_abs = sqrt(grad2(1)*grad2(1)+grad2(2)*grad2(2)+grad2(3)*grad2(3))
                         ! rc2013: make sure we don't divide by zero
                         if(grad1_abs > 0.0_DP) &
                              demft_dfull(i1,i2,islab12,is) = grad2_abs/grad1_abs
                      end do
                   end do
                end do
             end do
             ! rc2013: correct the EMFT XC potential by multiplying it by the
             ! quotient of the density gradients
             write(*,*) 'Correcting XC potential...'
             fxc_fine_emft = fxc_fine_emft*demft_dfull
          end if

          deallocate(full_dens_grad,stat=ierr)
          call utils_dealloc_check('lr_tddft_initialise','full_dens_grad',ierr)
          deallocate(emft_dens_grad,stat=ierr)
          call utils_dealloc_check('lr_tddft_initialise','emft_dens_grad',ierr)
          deallocate(demft_dfull,stat=ierr)
          call utils_dealloc_check('lr_tddft_initialise','demft_dfull',ierr)
          deallocate(recip_work,stat=ierr)
          call utils_dealloc_check('lr_tddft_initialise','recip_work',ierr)
          deallocate(dens_buffer,stat=ierr)
          call utils_dealloc_check('lr_tddft_initialise','dens_buffer',ierr)

       end if

       do is=1,pub_num_spins
          do isub=1,mdl%nsub
              fxc_fine(:,:,:,is,isub) = fxc_fine(:,:,:,is,isub) &
                   + full_fxc_fine(:,:,:,is)
              ! rc2013: only add the EMFT potential to the active region
              if(isub == pub_active_region .and. pub_emft) &
                   fxc_fine(:,:,:,is,isub) = fxc_fine(:,:,:,is,isub) &
                   + fxc_fine_emft(:,:,:,is)
           end do
        end do

        deallocate(full_fxc_fine,stat=ierr)
        call utils_dealloc_check('lr_tddft_initialise','full_fxc_fine',ierr)
        if (pub_emft) then
           deallocate(fxc_fine_emft,stat=ierr)
           call utils_dealloc_check('lr_tddft_initialise','fxc_fine_emft',ierr)
        end if

    endif

    ! now make sure that the ground state density corresponds to the valence
    ! density only, and no core density is included.
    if (pub_nlcc) then
       do is=1,pub_num_spins
          dens_fine(:,:,:,is) = dens_fine(:,:,:,is) - &
               mdl%core_density_fine * 0.5_DP * pub_spin_fac
       end do
    end if

    ! readd nhat density again if PAW method is used
    if (pub_paw) then
       if(.not. pub_nhat_in_xc) then
          dens_fine(:,:,:,:) = dens_fine(:,:,:,:)+nhat_den_grad(:,:,:,:,0)
       endif

       ! deallocate nhat_den
       deallocate(nhat_den_grad, stat=ierr)
       call utils_dealloc_check('lr_tddft_initialise','nhat_den_grad',ierr)
       deallocate(sub_nhat_den_grad, stat=ierr)
       call utils_dealloc_check('lr_tddft_initialise','sub_nhat_den_grad',ierr)
    endif

    ! Now diagonalise hamiltonians and get their eigenvalues/eigenvectors
    ! start with valence hamiltonian.
    ! This section of the code scales as O(N^3). Only required if we
    ! need to analyse the TDDFT transitions.
    if (pub_lr_tddft_analysis) then

       ! jcap: need to sum up the number of NGWFs over regions
       val_ngwf_num=0
       do isub=1,mdl%nsub
          val_ngwf_num=val_ngwf_num+val_ngwf_basis(isub)%num
       end do

       call dense_create(overlap_dense, val_ngwf_num, &
            val_ngwf_num, iscmplx=val_rep%overlap%p%iscmplx)
       call dense_create(ham_dense, val_ngwf_num, &
            val_ngwf_num, iscmplx=val_ham%ham(1)%p%iscmplx)
       do is=1, pub_num_spins
          call dense_convert(overlap_dense, val_rep%overlap)
          call dense_convert(ham_dense, val_ham%ham(is))
          call dense_eigensolve(val_ngwf_num, val_evals(:,is), ham_dense, &
               overlap_dense, 1,val_evecs(is))
       enddo

       ! clear memory for doing the same with conduction states.
       call dense_destroy(overlap_dense)
       call dense_destroy(ham_dense)

       ! now do the same for the conduction states..

       ! jcap: need to sum up the number of NGWFs over regions
       cond_ngwf_num=0
       do isub=1,mdl%nsub
          cond_ngwf_num=cond_ngwf_num+cond_ngwf_basis(isub)%num
       end do

       call dense_create(overlap_dense, cond_ngwf_num, &
            cond_ngwf_num, iscmplx=cond_rep%overlap%p%iscmplx)
       call dense_create(ham_dense, cond_ngwf_num, &
            cond_ngwf_num, iscmplx=cond_ham%ham(1)%p%iscmplx)
       do is=1, pub_num_spins
          call dense_convert(overlap_dense, cond_rep%overlap)
          call dense_convert(ham_dense, cond_ham%ham(is))
          call dense_eigensolve(cond_ngwf_num, cond_evals(:,is), ham_dense, &
               overlap_dense, 1, cond_evecs(is))
       enddo
    endif

    ! NOW set cond denskern to the projector onto unoccupied subspace
    ! only do so if required
    if (pub_lr_tddft_projector) then
       do is=1, pub_num_spins
          call sparse_embed_product(S_crossS_inv, cond_rep%cross_overlap,&
               cond_rep%inv_overlap)
          call sparse_embed_transpose(S_cross_trans, cond_rep%cross_overlap)
          call sparse_embed_product(S_crossKv, S_cross_trans, &
                val_denskern%m(is,PUB_1K))
          call sparse_embed_product(S_invS_crossKv, cond_rep%inv_overlap,&
               S_crossKv)
          call sparse_embed_product(eff_cond_denskern2, S_invS_crossKv, &
              S_crossS_inv)
          call sparse_embed_copy(eff_cond_denskern, cond_rep%inv_overlap)
          call sparse_embed_axpy(eff_cond_denskern,eff_cond_denskern2, -1.0_DP)

          call sparse_embed_copy(cond_denskern%kern%m(is,PUB_1K), &
                eff_cond_denskern)
       enddo
    endif

    ! NOW do exactly the same for joint set of NGWFs
    ! ndmh: Allocate SPAM3 structures for overlap matrix, kinetic energy,
    ! ndmh: density kernel, inverse overlap
    if (pub_on_root) write(stdout,*)
    call ngwf_rep_create(joint_rep, 'j', mdl, &
         is_cmplx=val_rep%overlap%p%iscmplx)
    call kernel_create(joint_denskern, 'K'//joint_rep%postfix, &
         is_cmplx=val_denskern%m(1,1)%p%iscmplx)

    ! ndmh: Allocate joint rep NGWFs
    ! jcap: loop over regions
    do isub=1,mdl%nsub
       call data_functions_alloc(joint_rep%ngwfs_on_grid(isub), &
            joint_ngwf_basis(isub)%size_on_grid, &
            iscmplx=val_rep%overlap%p%iscmplx)
       call comms_barrier
    end do

    ! tjz07: Initialise create joint ngwf set
    if (pub_on_root) write(stdout,'(a)',advance='no') 'Joint Valence + &
         &Conduction NGWF initialisation ...'

    ! ndmh: create a merged set of NGWFs containing the valence and
    ! conduction
    ! ndmh: NGWFs previously calculated
    ! jcap: loop over regions
    do isub=1,mdl%nsub
       call ngwfs_merge_sets(joint_rep%ngwfs_on_grid(isub),joint_ngwf_basis(isub), &
            mdl%cell,val_rep%ngwfs_on_grid(isub),val_ngwf_basis(isub), &
            cond_rep%ngwfs_on_grid(isub), cond_ngwf_basis(isub), &
            mdl%regions(isub)%par)
       call comms_barrier
    end do
    if (pub_on_root) write(stdout,'(a/)') '... done'
    call services_flush

    ! tjz07: Set the correct occupation numbers for the joint set.
    ! Using the joint set means that ALL conduction states representable
    ! by the joint representation are included in the calculation
    ! jcap: sum up ngwf nums over regions
    joint_ngwf_num=0
    do isub=1,mdl%nsub
       joint_ngwf_num=joint_ngwf_num+joint_ngwf_basis(isub)%num
    end do
    joint_rep%n_occ(:,:) = joint_ngwf_num-val_rep%n_occ(:,:)

    ! ndmh: allocate matrices in the joint hamiltonian
    call ngwf_ham_create(joint_ham,joint_rep)

    ! ndmh: copy in dijhat from valence hamiltonian
    if (pub_aug) then
       do is=1,pub_num_spins
          call sparse_embed_copy(joint_ham%dijhat(is), val_ham%dijhat(is))
       end do
    end if

    ! ndmh: calculate density-kernel-independent parts of the joint
    ! hamiltonian
    call hamiltonian_dens_indep_matrices(joint_rep, joint_ngwf_basis,&
         proj_basis, nl_projectors, hub_proj_basis, hub, mdl, val_rep,&
         val_ngwf_basis)

    ! JCW: Initialise HFx for evaluation of joint NGWF Hamiltonian here
    if (pub_use_hfx.or.pub_use_activehfx) then
       ! JCW: For EMFT calculations, HFX is only available in the active region
       !      For non-EMFT calculations, pub_active_region == 1 by default
       isub = pub_active_region
       ! JCW: Initialize joint_rep for HFX, as with cond_rep above
       call hf_exchange_init(hfxstate, joint_rep, &
            joint_ham%hfexchange(1)%m(isub,isub), mdl%cell, &
            size(mdl%regions(isub)%elements), init_rep_only = .true.)

       ! JCW: We need to construct the joint NGWF Hamiltonian, which
       !      requires mixed (joint-val) NGWF expansions in SWs. The expansion
       !      coefficients are not computed in ngwfs_initialise so need
       !      to be computed manually, as in conduction_ngwf_optimise. The
       !      SW_EX representing the expansion coeffs is stored in
       !      joint_rep%swexes(REP_SWEX_HFX_OTHER)
       call hf_exchange_dkn_indep_stage(hfxstate, mdl, isub, &
            mdl%regions(isub)%par, joint_rep, joint_ngwf_basis(isub), &
            rep2 = val_rep, ngwf_basis2 = val_ngwf_basis(isub), &
            basis_selector = [ 1,1,2,2,REP_SWEX_HFX_OTHER ])
            ! basis for alpha: 1 (joint)
            ! basis for beta: 1 (joint)
            ! basis for gamma: 2 (val)
            ! basis for delta: 2 (val)
            ! swex selector: HFX_OTHER (of joint_rep)
            ! (corresponds to ket_basis_selector = (/1,2,REP_SWEX_HFX_OTHER/)
            ! which is passed to swx_expand_my_ngwf_pairs inside
            ! hf_exchange_dkn_indep_stage)
            !
            ! Hamiltonian matrix elements are <alpha|H|beta> and the exchange
            ! contribution is <alpha gamma|1/r12|delta beta> P^{delta gamma}
            ! in physicists' notation (electron indices <12|12>), with P a
            ! suitable density matrix to contract with NGWFs gamma and delta.
            !
            ! For the joint Hamiltonian alpha, beta are joint NGWFs and gamma,
            ! delta (from the exchange operator) are valence NGWFs and thus
            ! P is the valence density kernel K

       ! JCW: Note that we should not now call ngwfs_initialise or any other
       !      routine that results in a call to ngwf_rep_register_change for
       !      joint_rep until the joint Hamiltonian is evaluated.
       !      ngwf_rep_register_changes calls swx_invalidate_expansion for
       !      each element of joint_rep%swexes(:), purging the coeffs_ht
       !      hash table for each SW_EX.
    end if

    ! ndmh: calculate the rest of the joint hamiltonian
    call hamiltonian_dens_dep_nonsc(joint_ham, joint_rep,&
         joint_ngwf_basis, mdl, hfxstate, lhxc_fine, hub, val_rep, val_ham, &
         val_denskern%m(:,PUB_1K),&
         val_ngwf_basis = val_ngwf_basis, &          ! optional arguments for HFx
         cond_dkn = joint_denskern%kern%m(:,PUB_1K)) ! optional arguments for HFx
    ! JCW: ^--- TODO It seems odd to pass an uninitialized joint density kernel
    !           here but hf_exchange_calculate requires it. It is used to
    !           evaluate the HFX energy and NGWF gradient, so probably not
    !           needed in this context (where we only require the X matrix)

    ! create new effective density kernel for joint ngwf basis. first
    ! start by destroying all old matrices
    if (pub_lr_tddft_projector) then
       call sparse_embed_destroy(eff_cond_denskern)
       call sparse_embed_destroy(S_crossS_inv)
       call sparse_embed_destroy(S_cross_trans)
       call sparse_embed_destroy(S_crossKv)
       call sparse_embed_destroy(S_invS_crossKv)
       call sparse_embed_destroy(eff_cond_denskern2)

       ! reallocating data
       call sparse_embed_create(eff_cond_denskern, joint_denskern%kern%m(1,1))
       call sparse_embed_create(S_crossS_inv,joint_rep%cross_overlap,&
            joint_rep%inv_overlap)
       call sparse_embed_transpose_structure(S_cross_trans%structure,&
            joint_rep%cross_overlap)
       call sparse_embed_create(S_cross_trans, &
            iscmplx=joint_rep%cross_overlap%p%iscmplx)
       call sparse_embed_create(S_crossKv, S_cross_trans, val_denskern%m(1,1))
       call sparse_embed_create(S_invS_crossKv, joint_rep%inv_overlap, &
            S_crossKv)
       call sparse_embed_create(eff_cond_denskern2, S_invS_crossKv,&
            S_crossS_inv)

       ! NOW set joint denskern...
       do is=1, pub_num_spins
          call sparse_embed_product(S_crossS_inv, joint_rep%cross_overlap,&
               joint_rep%inv_overlap)
          call sparse_embed_transpose(S_cross_trans, joint_rep%cross_overlap)
          call sparse_embed_product(S_crossKv, S_cross_trans, &
                val_denskern%m(is,PUB_1K))
          call sparse_embed_product(S_invS_crossKv, joint_rep%inv_overlap,&
               S_crossKv)
          call sparse_embed_product(eff_cond_denskern2, S_invS_crossKv, &
               S_crossS_inv)
          call sparse_embed_copy(eff_cond_denskern, joint_rep%inv_overlap)
          call sparse_embed_axpy(eff_cond_denskern,eff_cond_denskern2, -1.0_DP)

          call sparse_embed_copy(joint_denskern%kern%m(is,PUB_1K), &
                eff_cond_denskern)
       enddo
    endif

    ! now construct joint evecs for conduction space:
    do is=1, pub_num_spins
       call dense_create(joint_evecs(is), joint_ngwf_num, &
            joint_ngwf_num, iscmplx=joint_rep%overlap%p%iscmplx)
    enddo

    ! Again, this diagonalisation scales as O(N^3). It is
    ! only required if we need to analyse the TDDFT eigenstates
    if (pub_lr_tddft_analysis) then
       ! clear memory for calculating joint states.
       call dense_destroy(overlap_dense)
       call dense_destroy(ham_dense)

       ! now do the same for the joint conduction states..
       call dense_create(overlap_dense, joint_ngwf_num, &
            joint_ngwf_num, iscmplx=joint_rep%overlap%p%iscmplx)
       call dense_create(ham_dense, joint_ngwf_num, &
            joint_ngwf_num, iscmplx=joint_ham%ham(1)%p%iscmplx)
       do is=1, pub_num_spins
          call dense_convert(overlap_dense, joint_rep%overlap)
          call dense_convert(ham_dense, joint_ham%ham(is))
          call dense_eigensolve(joint_ngwf_num, joint_evals(:,is), ham_dense,&
               overlap_dense, 1,joint_evecs(is))
       enddo
    endif

    ! JCW: Now the HF exchange contributions to valence, conduction and joint
    !      Hamiltonians are evaluated. We do not need to preserve the SW_EX
    !      components of the valence, conduction or joint NGWF_REP instances.
    !      The SW expansions required for the HF contribution to the SCF
    !      response potential matrix (V_{SCF}^{HF}) are currently evaluated
    !      within linear_response_calc_SCF, via calls to
    !      hf_exchange_dkn_indep_stage, so no further HFX preparation is
    !      required here.
    ! JCW: @optimize
    !      It may be more efficient to precompute the needed SW expansions of
    !      the required NGWF products (val-val and cond/joint-cond/joint) here
    !      and store these in the swexes components of corresponding NGWF_REP
    !      instances.
    !      This would require some modification of the hf_exchange module since
    !      this currently only supports a single expansions_ht hash table
    !      (SWOP-expanded version of potential) tied to a particular NGWF
    !      product expansion. hf_exchange_dkn_indep_stage purges expansions_ht
    !      which is presumably then rebuilt based on the SW expansion of NGWF
    !      products that is subsequently evaluated (by hfx_expand_ngwf_pairs)

    ! destroy storage
    call sparse_embed_destroy(eff_cond_denskern)
    call sparse_embed_destroy(S_crossS_inv)
    call sparse_embed_destroy(S_cross_trans)
    call sparse_embed_destroy(S_crossKv)
    call sparse_embed_destroy(S_invS_crossKv)
    call sparse_embed_destroy(eff_cond_denskern2)
    if(pub_lr_tddft_analysis) then
       call dense_destroy(overlap_dense)
       call dense_destroy(ham_dense)
    endif

  end subroutine lr_tddft_initialise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_gram_schmidt(set_of_vecs,set_of_vecs_cov, num_vecs, &
       cond_overlap, val_overlap)

    !=======================================================================!
    ! Subroutine does a gram_schmidt orthonormalisation such that the trial !
    ! evecs stored in evec_mat are orthonormal with respect to the overlap  !
    ! matrix S_dense and the ORIGINAL ground state evecs. Ie if we denote   !
    ! the ground state evecs as x0 and the evecs in evec mat as x1 and we   !
    ! denote the overlap matrix as S, we require that (x1)^T*S*x0=I, where  !
    ! I is the identity matrix.                                             !
    ! Note that if this is a spin polarised calculation, what we require is !
    ! That the two vectors are orthogonal in EACH spin channel.             !
    ! Modified for embedding by Joseph Prentice, July 2018                  !
    !=======================================================================!

    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_transpose, sparse_embed_transpose_structure, &
         sparse_embed_trace, sparse_embed_copy
    use rundat, only: pub_lr_tddft_num_conv_states,pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: num_vecs  ! defines the number of vectors
    type(SPAM3_EMBED), intent(inout) :: set_of_vecs(num_vecs,pub_num_spins)  ! The total set of
    ! num_vec batches of vectors. gram schmidt will make all
    ! of them orthogonal to each other
    type(SPAM3_EMBED), intent(inout) :: set_of_vecs_cov(num_vecs,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_overlap
    type(SPAM3_EMBED), intent(in) :: val_overlap

    ! Local Variables
    integer :: ierr, is
    integer :: count_i
    integer :: count_j
    type(SPAM3_EMBED), allocatable, dimension(:) :: vec1_trans
    type(SPAM3_EMBED), allocatable, dimension(:) :: vec2_trans
    real(kind=DP) :: temp1, temp2, temp_scaling, trace

    allocate(vec1_trans(pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_gram_schmidt','vec1_trans',ierr)
    allocate(vec2_trans(pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_gram_schmidt','vec2_trans',ierr)

    do is=1, pub_num_spins
       call sparse_embed_transpose_structure(vec1_trans(is)%structure,set_of_vecs(1,is))
       call sparse_embed_create(vec1_trans(is), &
            iscmplx=set_of_vecs(1,is)%p%iscmplx)
       call sparse_embed_create(vec2_trans(is), vec1_trans(is))
    enddo

    ! check if states are already converged. In that case
    ! orthogonalise all unconverged states against them.
    if(pub_lr_tddft_num_conv_states>0) then
       do count_j=1+pub_lr_tddft_num_conv_states, num_vecs
          do is=1, pub_num_spins
             call sparse_embed_transpose(vec1_trans(is),set_of_vecs(count_j,is))
          enddo
          do count_i=1, pub_lr_tddft_num_conv_states
             temp1=0.0_DP
             temp2=0.0_DP
             do is=1, pub_num_spins
                call sparse_embed_trace(trace,vec1_trans(is),set_of_vecs_cov(count_i,is))
                temp1=temp1+trace

                call sparse_embed_transpose(vec2_trans(is),set_of_vecs(count_i,is))
                call sparse_embed_trace(trace,vec2_trans(is),set_of_vecs(count_i,is))
                temp2=temp2+trace
             enddo

             do is=1, pub_num_spins
                call sparse_embed_axpy(set_of_vecs(count_j,is),&
                       set_of_vecs(count_i,is),-1.0_DP*temp1/temp2)
                call sparse_embed_axpy(set_of_vecs_cov(count_j,is),&
                     set_of_vecs_cov(count_i,is),-1.0_DP*temp1/temp2)
             enddo
          enddo
       enddo
    endif

    ! start regular Gram schmidt only for unconverged
    ! vectors
    do count_j=1+pub_lr_tddft_num_conv_states, num_vecs
       do is=1, pub_num_spins
          call sparse_embed_transpose(vec1_trans(is), set_of_vecs(count_j,is))
       enddo

       do count_i=1+pub_lr_tddft_num_conv_states, count_j-1
          ! get vec2
          temp1=0.0_DP
          temp2=0.0_DP
          do is=1, pub_num_spins
             call sparse_embed_transpose(vec2_trans(is), set_of_vecs(count_i,is))
             call sparse_embed_trace(trace, vec1_trans(is), set_of_vecs_cov(count_i,is))
             temp1=temp1+trace
             call sparse_embed_trace(trace, vec2_trans(is), set_of_vecs_cov(count_i,is))
             temp2=temp2+trace
          enddo

          do is=1,pub_num_spins
             temp_scaling=temp1/temp2
             !call sparse_scale(vec2(is),temp_scaling)
             ! Successfully calculated the projector.
             call sparse_embed_axpy(set_of_vecs(count_j,is), &
                   set_of_vecs(count_i,is), -1.0_DP*temp_scaling)
             call sparse_embed_axpy(set_of_vecs_cov(count_j,is),&
                   set_of_vecs_cov(count_i,is),-1.0_DP*temp_scaling)
          enddo
       enddo

       ! normalize
       temp1=0.0_DP
       do is=1, pub_num_spins
          call sparse_embed_transpose(vec1_trans(is), set_of_vecs(count_j,is))
          call sparse_embed_trace(trace, vec1_trans(is), set_of_vecs_cov(count_j,is))
          temp1=temp1+trace
       enddo

       do is=1, pub_num_spins
          call sparse_embed_scale(set_of_vecs(count_j,is), 1.0_DP/sqrt(temp1))
          call sparse_embed_scale(set_of_vecs_cov(count_j,is),1.0_DP/sqrt(temp1))
       enddo
    enddo

    do is=1, pub_num_spins
       call sparse_embed_destroy(vec1_trans(is))
       call sparse_embed_destroy(vec2_trans(is))
    enddo

    deallocate(vec1_trans,stat=ierr)
    call utils_dealloc_check('lr_tddft_gram_schmidt','vec1_trans',ierr)
    deallocate(vec2_trans,stat=ierr)
    call utils_dealloc_check('lr_tddft_gram_schmidt','vec2_trans',ierr)

  end subroutine lr_tddft_gram_schmidt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_calc_omega(omega,evec, evec_cov,y_mat, cond_overlap,&
       val_overlap)

    !=======================================================================!
    ! Subroutine returning the norm omega that is minimised in the process  !
    ! of conjugate gradient calculations.                                   !
    ! Omega is given by the trace of (x*y)/(x*S*x*S)                        !
    ! Modified for embedding by Joseph Prentice, July 2018                  !
    !=======================================================================!

    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_transpose_structure, sparse_embed_trace, sparse_embed_transpose
    use rundat, only: pub_num_spins

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: omega
    type(SPAM3_EMBED), intent(in) :: evec(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: evec_cov(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: y_mat(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_overlap
    type(SPAM3_EMBED), intent(in) :: val_overlap

    ! Local Variables
    real(kind=DP) :: temp_sum, temp_sum2, trace
    type(SPAM3_EMBED) :: evec_trans
    integer :: is

    call sparse_embed_transpose_structure(evec_trans%structure, evec(1))
    call sparse_embed_create(evec_trans, iscmplx=evec(1)%p%iscmplx)

    temp_sum=0.0_DP
    temp_sum2=0.0_DP
    do is=1, pub_num_spins
       call sparse_embed_transpose(evec_trans, evec(is))
       ! Calculate trace of x*y.
       call sparse_embed_trace(trace, evec_trans, y_mat(is))
       temp_sum=temp_sum+trace
       ! calculate the trace of x*x
       call sparse_embed_trace(trace, evec_trans, evec_cov(is))
       temp_sum2=temp_sum2+trace
    enddo

    omega=temp_sum/temp_sum2

    call sparse_embed_destroy(evec_trans)

  end subroutine lr_tddft_calc_omega

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_energy_lambda(lambda, evec, evec_cov,D_vec,D_vec_cov,&
       y_vec,D_on_operator, cond_overlap, val_overlap, num_states,&
       energy)

    !========================================================================!
    ! Routine calculates the energy of the system after a step of step lenght!
    ! lambda into the direction D_vec has been taken. ie x --> x+lambda*D.   !
    ! However, the energy is calculated without calculating the total action !
    ! of x+lambda*D on the 2-particle hamiltonian. The energy is purely      !
    ! calculated by knowing the individual actions of H2p on D_vec and x_vec.!
    ! This routine is needed in the line minimisation, where 2 possible ideal!
    ! step lenghts are calculated and we need to decide which one corresponds!
    ! to the minimum energy.                                                 !
    ! Modified for embedding by Joseph Prentice, July 2018                   !
    !========================================================================!

    use comms, only: comms_barrier
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_transpose, sparse_embed_transpose_structure, &
         sparse_embed_trace, sparse_embed_copy
    use rundat, only: pub_lr_tddft_num_conv_states,pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: num_states
    real(kind=DP), intent(in) :: lambda
    type(SPAM3_EMBED), intent(in) :: evec(num_states,pub_num_spins) ! matrix
    type(SPAM3_EMBED), intent(in) :: evec_cov(num_states,pub_num_spins) !covariant
    type(SPAM3_EMBED), intent(in) :: D_vec(num_states,pub_num_spins) ! matrix
    type(SPAM3_EMBED), intent(in) :: D_vec_cov(num_states,pub_num_spins) !covariant
    type(SPAM3_EMBED), intent(in) :: y_vec(num_states,pub_num_spins) !The result of a batch of x
    ! vectors operating on operator
    type(SPAM3_EMBED), intent(in) :: D_on_operator(num_states,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_overlap
    type(SPAM3_EMBED), intent(in) :: val_overlap
    real(kind=DP), intent(inout) :: energy

    ! Local Variables
    type(SPAM3_EMBED), allocatable, dimension(:) :: k_trans, D_trans
    integer :: counter, ierr, is
    real(kind=DP) :: temp, norm, temp_energy, trace

    allocate(k_trans(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_energy_lambda','k_trans',ierr)
    allocate(D_trans(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_energy_lambda','D_trans',ierr)

    do is=1, pub_num_spins
       call sparse_embed_transpose_structure(k_trans(is)%structure, evec(1,is))
       call sparse_embed_transpose_structure(D_trans(is)%structure, D_vec(1,is))
       call sparse_embed_create(k_trans(is), iscmplx=evec(1,is)%p%iscmplx)
       call sparse_embed_create(D_trans(is), iscmplx=D_vec(1,is)%p%iscmplx)
    enddo

    ! sum over all states
    energy=0.0_DP
    do counter=1+pub_lr_tddft_num_conv_states, num_states
       norm=0.0_DP
       temp_energy=0.0_DP
       do is=1, pub_num_spins
          call sparse_embed_transpose(k_trans(is), evec(counter,is))
          call sparse_embed_transpose(D_trans(is), D_vec(counter,is))
       enddo

       ! start with energy, then do normalisation
       temp=0.0_DP
       do is=1, pub_num_spins
          call sparse_embed_trace(trace, k_trans(is), y_vec(counter,is))
          temp_energy=temp_energy+trace
          call sparse_embed_trace(trace, D_trans(is), y_vec(counter,is))
          temp = temp+trace
       enddo
       temp_energy=temp_energy+temp*lambda
       temp=0.0_DP
       do is=1, pub_num_spins
          call sparse_embed_trace(trace, k_trans(is), D_on_operator(counter,is))
          temp= temp+trace
       enddo
       temp_energy=temp_energy+temp*lambda
       temp=0.0_DP
       do is=1, pub_num_spins
          call sparse_embed_trace(trace, D_trans(is), D_on_operator(counter,is))
          temp=temp+trace
       enddo
       temp_energy=temp_energy+temp*lambda*lambda

       ! do normalisation for this eigenvalue
       temp=0.0_DP
       do is=1, pub_num_spins
          call sparse_embed_trace(trace, k_trans(is), evec_cov(counter,is))
          norm=norm+trace
          call sparse_embed_trace(trace, D_trans(is), evec_cov(counter,is))
          temp=temp+trace
       enddo
       norm=norm+temp*2.0_DP*lambda
       temp=0.0_DP
       do is=1, pub_num_spins
          call sparse_embed_trace(trace, D_trans(is), D_vec_cov(counter,is))
          temp=temp+trace
       enddo
       norm=norm+lambda*lambda*temp

       ! add the energy of this eigenvalue to the total energy
       energy=energy+temp_energy/norm
    enddo

    do is=1, pub_num_spins
       call sparse_embed_destroy(k_trans(is))
       call sparse_embed_destroy(D_trans(is))
    enddo
    deallocate(k_trans,stat=ierr)
    call utils_dealloc_check('lr_tddft_energy_lambda','k_trans',ierr)
    deallocate(D_trans,stat=ierr)
    call utils_dealloc_check('lr_tddft_energy_lambda','D_trans',ierr)

  end subroutine lr_tddft_energy_lambda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_line_min(lambda_opt, evec,evec_cov, D_vec, y_vec, &
       mdl, hfxstate, response_dens, sub_response_dens, fxc_fine, &
       cond_ngwf_basis, cond_rep, cond_denskern, kchc, kcSc, &
       val_ngwf_basis, val_rep, val_denskern, hvkv, Svkv, &
       num_states, energy_conv_states,L_vec,ground_state_dens,&
       sub_ground_state_dens,ground_state_hartree)

    !========================================================================!
    ! Subroutine calculating the optimum lambda by solving a quadratic       !
    ! equation of the form a*lambda^2+b*lambda+c=0. The equation is obtained !
    ! by minimizing the energy with respect to lambda.                       !
    ! Modified for embedding by Joseph Prentice, July 2018                   !
    !========================================================================!

    use comms, only: comms_barrier, pub_on_root
    use function_basis, only: FUNC_BASIS
    use constants, only: stdout, NORMAL, VERBOSE
    use hf_exchange, only: HFX_STATE
    use linear_response, only: linear_response_operator
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_lr_tddft_num_conv_states, pub_num_spins,&
         pub_multigrid_hartree, pub_lr_tddft_preopt,pub_output_detail
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_transpose, sparse_embed_transpose_structure, &
         sparse_embed_trace, sparse_embed_copy
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: num_states
    real(kind=DP), intent(inout) :: energy_conv_states
    real(kind=DP), intent(inout) :: lambda_opt
    type(SPAM3_EMBED), intent(inout) :: evec(num_states,pub_num_spins) ! matrix
    type(SPAM3_EMBED), intent(inout) :: evec_cov(num_states,pub_num_spins) ! covariant vec
    type(SPAM3_EMBED), intent(inout) :: D_vec(num_states,pub_num_spins) ! matrix
    type(SPAM3_EMBED), intent(inout) :: y_vec(num_states,pub_num_spins) !The result of a batch of
    ! x vectors operating on operator
    type(MODEL), intent(in) :: mdl
    type (HFX_STATE), intent(inout), target :: hfxstate
    real(kind=DP), intent(inout) :: response_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_response_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)
    real(kind=DP), intent(in) :: fxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,2,mdl%nsub)
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis(:)
    type(NGWF_REP), intent(inout) :: cond_rep
    type(SPAM3_EMBED), intent(in) :: cond_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: kchc(pub_num_spins)
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type(NGWF_REP), intent(inout) :: val_rep
    type(SPAM3_EMBED), intent(in) :: val_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: hvkv(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: Svkv(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: kcSc(pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: L_vec(num_states)
    real(kind=DP), intent(inout) :: ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)
    real(kind=DP), optional, intent(in) :: ground_state_hartree(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)

    ! Local Variables
    type(SPAM3_EMBED), allocatable, dimension(:,:) :: trial_x1, &
         operator_on_D_vec
    type(SPAM3_EMBED), allocatable, dimension(:,:) :: trial_x2
    type(SPAM3_EMBED), allocatable, dimension(:,:) :: D_vec_cov
    type(SPAM3_EMBED), allocatable, dimension(:) :: kSc,kSv
    type(SPAM3_EMBED) :: temp_transpose, temp_evec, KSvkv
    type(SPAM3_EMBED) :: LSvkv
    real(kind=DP) :: temp1, temp2, a, b, c, trace
    real(kind=DP) :: discriminant, lambda1, lambda2, omega1, omega2
    integer :: icount,jcount, ierr,is

    allocate(D_vec_cov(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_line_min','D_vec_cov',ierr)
    allocate(trial_x1(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_line_min','trial_x1', ierr)
    allocate(trial_x2(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_line_min','trial_x2', ierr)
    allocate(operator_on_D_vec(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_line_min','operator_on_D_vec', &
         ierr)
    allocate(kSc(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_line_min','kSc',ierr)
    allocate(kSv(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_line_min','kSv',ierr)

    call sparse_embed_transpose_structure(temp_transpose%structure, evec(1,1))
    call sparse_embed_create(temp_transpose, iscmplx=evec(1,1)%p%iscmplx)
    do is=1, pub_num_spins
       call sparse_embed_create(kSc(is), temp_transpose, cond_rep%overlap)
       call sparse_embed_create(kSv(is), evec(1,is), val_rep%overlap)
       do icount=1, num_states
          call sparse_embed_create(trial_x1(icount,is), evec(1,is))
          call sparse_embed_create(trial_x2(icount,is), evec(1,is))
          call sparse_embed_create(operator_on_D_vec(icount,is), y_vec(1,is))
          call sparse_embed_create(D_vec_cov(icount,is),y_vec(1,is))
       enddo
    enddo
    ! enforce new trial solution to be only composed of states in
    ! 'physical' subspace
    call sparse_embed_create(temp_evec, evec(1,1))
    call sparse_embed_create(KSvkv, temp_evec, Svkv(1))
    call sparse_embed_create(LSvkv,L_vec(1),Svkv(1))

    ! calculate total operator acting the batch of D vectors
    do icount=1+pub_lr_tddft_num_conv_states, num_states
       ! create correct, truncated, physical gradient...
       do is=1, pub_num_spins
          call sparse_embed_copy(L_vec(icount), D_vec(icount,is))
          call sparse_embed_product(LSvkv,L_vec(icount), Svkv(is))
          call sparse_embed_product(D_vec(icount,is),kcSc(is),LSvkv)
       enddo
       ! now for CG algorithm to work, D_vec(icount) must still be
       ! orthogonal to all previous evecs. This should preserve ortho-
       ! normality of the new eigenvectors to first order.
       do is=1, pub_num_spins
          ! create D_cov_vec
          call sparse_embed_product(ksv(is),D_vec(icount,is), val_rep%overlap)
          call sparse_embed_product(D_vec_cov(icount,is),cond_rep%overlap,ksv(is))
       enddo
       do jcount=1, num_states
          temp1=0.0_DP
          do is=1, pub_num_spins
             call sparse_embed_transpose(temp_transpose, evec(jcount,is))
             !call sparse_product(kSc(is), temp_transpose,cond_rep%overlap)
             call sparse_embed_trace(trace,temp_transpose,D_vec_cov(icount,is))
             temp1=temp1+trace
          enddo

          do is=1, pub_num_spins
             call sparse_embed_axpy(D_vec(icount,is),evec(jcount,is), -1.0_DP*temp1)
             call sparse_embed_axpy(D_vec_cov(icount,is),evec_cov(jcount,is),-1.0_DP*temp1)
          enddo
       enddo

       ! check if ground_state_hartree potential is needed in the routine
       if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
          call linear_response_operator(operator_on_D_vec(icount,:),&
               D_vec(icount,:),  cond_rep, cond_ngwf_basis, cond_denskern,&
               kchc, val_rep, val_ngwf_basis, val_denskern, hvkv,&
               mdl, hfxstate, response_dens, sub_response_dens, fxc_fine, &
               ground_state_dens, sub_ground_state_dens, &
               ground_state_hartree=ground_state_hartree)
       else
          call linear_response_operator(operator_on_D_vec(icount,:),&
               D_vec(icount,:),  cond_rep, cond_ngwf_basis, cond_denskern,&
               kchc, val_rep, val_ngwf_basis, val_denskern, hvkv,&
               mdl, hfxstate, response_dens, sub_response_dens, fxc_fine, &
               ground_state_dens, sub_ground_state_dens)
       endif
    enddo

    ! Now construct c, b and a coefficients.
    ! Start with a:
    a=0.0_DP
    do icount=1+pub_lr_tddft_num_conv_states, num_states
       temp2=0.0_DP
       temp1=0.0_DP
       do is=1, pub_num_spins
          call sparse_embed_transpose(temp_transpose, D_vec(icount,is))
          call sparse_embed_trace(trace,temp_transpose,y_vec(icount,is))
          temp2=temp2+trace

          call sparse_embed_transpose(temp_transpose, evec(icount,is))
          call sparse_embed_trace(trace,temp_transpose,evec_cov(icount,is))
          temp1=temp1+trace
       enddo
       a=a+temp1*temp2

       temp1=0.0_DP
       temp2=0.0_DP
       do is=1, pub_num_spins
          call sparse_embed_transpose(temp_transpose,evec(icount,is))
          call sparse_embed_trace(trace,temp_transpose,y_vec(icount,is))
          temp1=temp1+trace
          call sparse_embed_transpose(temp_transpose, D_vec(icount,is))
          call sparse_embed_trace(trace,temp_transpose,evec_cov(icount,is))
          temp2=temp2+trace
       enddo
       a=a-temp1*temp2
    enddo

    ! now construct b
    b=0.0_DP
    do icount=1+pub_lr_tddft_num_conv_states, num_states
       temp1=0.0_DP
       temp2=0.0_DP
       do is=1, pub_num_spins
          call sparse_embed_transpose(temp_transpose, evec(icount,is))
          call sparse_embed_trace(trace,temp_transpose,evec_cov(icount,is))
          temp1=temp1+trace

          call sparse_embed_transpose(temp_transpose, D_vec(icount,is))
          call sparse_embed_trace(trace,temp_transpose,operator_on_D_vec(icount,is))
          temp2=temp2+trace
       enddo
       b=b+temp1*temp2

       temp1=0.0_DP
       temp2=0.0_DP
       do is=1,pub_num_spins
          call sparse_embed_transpose(temp_transpose, D_vec(icount,is))
          call sparse_embed_trace(trace,temp_transpose,D_vec_cov(icount,is))
          temp1=temp1+trace

          call sparse_embed_transpose(temp_transpose, evec(icount,is))
          call sparse_embed_trace(trace,temp_transpose,y_vec(icount,is))
          temp2=temp2+trace
       enddo
       b=b-temp1*temp2
    enddo

    ! finally construct c:
    c=0.0_DP
    do icount=1+pub_lr_tddft_num_conv_states, num_states
       temp1=0.0_DP
       temp2=0.0_DP
       do is=1, pub_num_spins
          call sparse_embed_transpose(temp_transpose, evec(icount,is))
          call sparse_embed_trace(trace,temp_transpose,D_vec_cov(icount,is))
          temp1=temp1+trace

          call sparse_embed_transpose(temp_transpose, D_vec(icount,is))
          call sparse_embed_trace(trace,temp_transpose,operator_on_D_vec(icount,is))
          temp2=temp2+trace
       enddo
       c=c+temp1*temp2

       temp1=0.0_DP
       temp2=0.0_DP
       do is=1, pub_num_spins
          call sparse_embed_transpose(temp_transpose, D_vec(icount,is))
          call sparse_embed_trace(trace,temp_transpose,D_vec_cov(icount,is))
          temp1=temp1+trace
          call sparse_embed_transpose(temp_transpose, evec(icount,is))
          call sparse_embed_trace(trace,temp_transpose,operator_on_D_vec(icount,is))
          temp2=temp2+trace
       enddo
       c=c-temp1*temp2
    enddo

    if (pub_on_root .and. pub_output_detail>=VERBOSE) then
       write(stdout,'(a,E12.4,a)') '****** Line-min: c coefficient= ',c,' *****'
       write(stdout,'(a,E12.4,a)') '****** Line-min: b coefficient= ',b,' *****'
       write(stdout,'(a,E12.4,a)') '****** Line-min: a coefficient= ',a,' *****'
    end if

    discriminant=(b*b)-4.0_DP*a*c
    if (discriminant<0.0_DP .and. abs(discriminant)>10e-15_DP) then
       if (pub_on_root) write(stdout, '(a)') 'ERROR, discriminant in&
            & lr_tddft_line_min is negative'
       if (pub_on_root) write(stdout, *) discriminant
       lambda_opt=-b/(2.0_DP*c) ! Set it to b/2.0. Corresponds to discriminant =0.
    else if (discriminant<0.0_DP .and. abs(discriminant)<10e-15_DP) then
       lambda_opt=-b/(2.0_DP*c)
    else
       if(b<0.0_DP) then
          lambda1=(-b+sqrt(discriminant))/(2.0_DP*c)
       else
          lambda1=(-b-sqrt(discriminant))/(2.0_DP*c)
       endif

       lambda2=a/(c*lambda1)

       ! construct trial vector. Also ensure that trial_x1 is only
       ! composed of 'physical' evecs

       do icount=1+pub_lr_tddft_num_conv_states, num_states
          do is=1, pub_num_spins
             call sparse_embed_copy(trial_x1(icount,is), evec(icount,is))
             call sparse_embed_axpy(trial_x1(icount,is), D_vec(icount,is), lambda1)
          enddo
       enddo


       ! Calculate the energy corresponding to trial vector batch 1:
       omega1=0.0_DP
       call lr_tddft_energy_lambda(lambda1, evec,evec_cov, D_vec,D_vec_cov,&
            y_vec, operator_on_D_vec, cond_rep%overlap, val_rep%overlap, &
            num_states, omega1)
       omega1=omega1+energy_conv_states

       ! Now construct vector batch 2 and project into physical states
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          do is=1, pub_num_spins
             call sparse_embed_copy(trial_x2(icount,is), evec(icount,is))
             call sparse_embed_axpy(trial_x2(icount,is), D_vec(icount,is), lambda2)
          enddo
       enddo


       ! calculate the energy corresponding to trial vector batch 2:
       omega2=0.0_DP
       call lr_tddft_energy_lambda(lambda2, evec,evec_cov, D_vec,&
            D_vec_cov, y_vec, operator_on_D_vec, cond_rep%overlap,&
            val_rep%overlap, num_states, omega2)
       omega2=omega2+energy_conv_states

       ! check which step size produces smaller energy
       if (omega2<omega1) then
          lambda_opt=lambda2
          ! Set evec
          do icount=1+pub_lr_tddft_num_conv_states, num_states
             do is=1, pub_num_spins
                call sparse_embed_copy(evec(icount,is), trial_x2(icount,is))
                call sparse_embed_axpy(evec_cov(icount,is),D_vec_cov(icount,is),&
                   lambda2)
             enddo
          enddo
          if (pub_on_root .and. pub_output_detail>=NORMAL) &
             write(stdout,'(a,f19.14,a)') &
               '****** Line-min: predicted omega= ',omega2,' *****'

       else
          lambda_opt=lambda1
          do icount=1+pub_lr_tddft_num_conv_states, num_states
             do is=1, pub_num_spins
                call sparse_embed_copy(evec(icount,is), trial_x1(icount,is))
                call sparse_embed_axpy(evec_cov(icount,is),D_vec_cov(icount,is),&
                   lambda1)
             enddo
          enddo
          if (pub_on_root .and. pub_output_detail>=NORMAL) &
             write(stdout,'(a,f19.14,a)') &
               '****** Line-min: predicted omega= ',omega1,' *****'
       endif
    endif

    ! deallocate data structures
    do is=1, pub_num_spins
       call sparse_embed_destroy(kSv(is))
       call sparse_embed_destroy(kSc(is))
       do icount=1, num_states
          call sparse_embed_destroy(trial_x1(icount,is))
          call sparse_embed_destroy(trial_x2(icount,is))
          call sparse_embed_destroy(D_vec_cov(icount,is))
          call sparse_embed_destroy(operator_on_D_vec(icount,is))
       enddo
    enddo
    call sparse_embed_destroy(temp_evec)
    call sparse_embed_destroy(KSvkv)
    deallocate(kSv,stat=ierr)
    call utils_dealloc_check('lr_tddft_line_min','kSv',ierr)
    deallocate(kSc,stat=ierr)
    call utils_dealloc_check('lr_tddft_line_min','kSc',ierr)
    deallocate(trial_x1, stat=ierr)
    call utils_dealloc_check('lr_tddft_line_min', 'trial_x1',ierr)
    deallocate(trial_x2, stat=ierr)
    call utils_dealloc_check('lr_tddft_line_min','trial_x2', ierr)
    deallocate(operator_on_D_vec, stat=ierr)
    call utils_dealloc_check('lr_tddft_line_min','operator_on_D_vec', &
         ierr)
    deallocate(D_vec_cov,stat=ierr)
    call utils_dealloc_check('lr_tddft_line_min','D_vec_cov',ierr)


  end subroutine lr_tddft_line_min

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_check_convergence(x_vec,x_vec_cov,y_vec,subspace_eval,&
         subspace_eval_prev, evec_mat, eval_mat,identity_mat,&
         states_converged,energy_conv_states,num_states)

    !=====================================================================!
    ! Subroutine checks for the convergence of individual eigenstates in  !
    ! the entire converged manifold. If any eigenstates have not changed  !
    ! within some tolerance from the previous values a certain number of  !
    ! iterations ago, these vectors will be removed from the manifold and !
    ! the calculation will be continued with a smaller subset of vectors  !
    ! Modified for embedding by Joseph Prentice, July 2018                !
    !=====================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, NORMAL
    use linalg, only: linalg_dsygv_lt
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_trace, sparse_embed_create, &
         sparse_embed_axpy, sparse_embed_transpose_structure, &
         sparse_embed_transpose, sparse_embed_copy, sparse_embed_scale, &
         sparse_embed_destroy
    use rundat, only: pub_lr_tddft_num_conv_states, &
           pub_lr_tddft_cg_threshold, pub_output_detail,pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: num_states
    type(SPAM3_EMBED), intent(inout) :: x_vec(num_states,pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: x_vec_cov(num_states,pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: y_vec(num_states,pub_num_spins)
    real(kind=DP), intent(inout) :: subspace_eval(num_states-&
        pub_lr_tddft_num_conv_states)
    real(kind=DP), intent(inout) :: subspace_eval_prev(num_states-&
        pub_lr_tddft_num_conv_states)
    real(kind=DP), intent(inout) :: evec_mat(num_states-&
        pub_lr_tddft_num_conv_states,num_states-pub_lr_tddft_num_conv_states)
    real(kind=DP), intent(inout) :: eval_mat(num_states-&
        pub_lr_tddft_num_conv_states,num_states-pub_lr_tddft_num_conv_states)
    real(kind=DP), intent(inout) :: identity_mat(num_states-&
        pub_lr_tddft_num_conv_states,num_states-pub_lr_tddft_num_conv_states)
    integer, intent(inout) :: states_converged
    real(kind=DP), intent(inout) :: energy_conv_states

    ! temp variables:
    integer :: icount, jcount, temp_icount,temp_jcount,subspace_size, ierr
    type(SPAM3_EMBED), allocatable, dimension(:,:) :: temp_x, temp_x_cov
    type(SPAM3_EMBED), allocatable, dimension(:) :: temp_transpose
    real(kind=DP) :: temp_result, trace
    logical :: any_converged
    integer :: is

    subspace_size=num_states-pub_lr_tddft_num_conv_states
    states_converged=0
    any_converged=.false.

    ! allocate space
    allocate(temp_transpose(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_check_convergence','temp_transpose',ierr)
    allocate(temp_x(subspace_size,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_check_convergence','temp_x',ierr)
    allocate(temp_x_cov(subspace_size,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_check_convergence','temp_x_cov',ierr)
    do is=1, pub_num_spins
       call sparse_embed_transpose_structure(temp_transpose(is)%structure,x_vec(1,is))
       call sparse_embed_create(temp_transpose(is), iscmplx=x_vec(1,is)%p%iscmplx)
    enddo
    do icount=1, subspace_size
       do is=1, pub_num_spins
          call sparse_embed_create(temp_x(icount,is), x_vec(1,is))
          call sparse_embed_create(temp_x_cov(icount,is),x_vec_cov(1,is))
       enddo
    enddo

    ! Construct the matrix and diagonalise it explicitly
    do icount=1, subspace_size
       do is=1, pub_num_spins
          call sparse_embed_transpose(temp_transpose(is), x_vec(icount+&
               pub_lr_tddft_num_conv_states,is))
       enddo
       do jcount=icount, subspace_size
          temp_result=0.0_DP
          do is=1, pub_num_spins
             call sparse_embed_trace(trace, temp_transpose(is), &
                  y_vec(jcount+pub_lr_tddft_num_conv_states,is))
             temp_result=temp_result+trace
          enddo
          eval_mat(icount,jcount)=temp_result
          eval_mat(jcount,icount)=temp_result
          if (icount==jcount) then
             identity_mat(icount,icount)=1.0_DP
          else
             identity_mat(icount,jcount)=0.0_DP
             identity_mat(jcount,icount)=0.0_DP
          endif
       enddo
    enddo

    ! call internal eigensolve to diagonalise converged subspace
    evec_mat=eval_mat
    call linalg_dsygv_lt(evec_mat,subspace_eval,identity_mat,subspace_size)

    if(pub_on_root .and. pub_output_detail>=NORMAL) &
        write(stdout,'(a)') '****** Checking for convergence &
        &of individual states *****'
    ! first check if ANY state has converged:
    do icount=1, subspace_size
       temp_result=subspace_eval_prev(icount)-subspace_eval(icount)
       if(abs(temp_result)<pub_lr_tddft_cg_threshold) then
         any_converged=.true.
         energy_conv_states=energy_conv_states+subspace_eval(icount)
       endif
    enddo


    ! only construct new x_vec and reorder etc if any of the states is actually
    ! converged. Otherwise dont
    if(any_converged) then

       ! now need to construct 'real' k from a superposition of Ks, determined
       ! through eval mat.
       do icount=1, subspace_size
          do is=1, pub_num_spins
             call sparse_embed_copy(temp_x(icount,is),x_vec(icount+&
               pub_lr_tddft_num_conv_states,is))
             call sparse_embed_copy(temp_x_cov(icount,is),x_vec_cov(icount+&
               pub_lr_tddft_num_conv_states,is))
          enddo
       enddo


       ! now construct correct response density matrices
       ! seems to be done correctly
       do icount=1, subspace_size
          temp_icount=icount+pub_lr_tddft_num_conv_states
          do jcount=1, subspace_size
             temp_jcount=jcount+pub_lr_tddft_num_conv_states
             do is=1, pub_num_spins
                if (jcount==1) then
                   call sparse_embed_copy(x_vec(temp_icount,is), temp_x(jcount,is))
                   call sparse_embed_copy(x_vec_cov(temp_icount,is),&
                      temp_x_cov(jcount,is))
                   temp_result=evec_mat(jcount,icount)
                   call sparse_embed_scale(x_vec(temp_icount,is), temp_result)
                   call sparse_embed_scale(x_vec_cov(temp_icount,is),temp_result)
                else
                   temp_result=evec_mat(jcount,icount)
                   call sparse_embed_axpy(x_vec(temp_icount,is), temp_x(jcount,is), &
                       temp_result)
                   call sparse_embed_axpy(x_vec_cov(temp_icount,is),&
                        temp_x_cov(jcount,is), temp_result)
                endif
             enddo
          enddo
       enddo

       ! now check for convergence with previous set of evecs:
       do icount=1, subspace_size
          temp_result=subspace_eval_prev(icount)-subspace_eval(icount)
          if(abs(temp_result)<pub_lr_tddft_cg_threshold) then
             ! state is converged. increase the counter
             states_converged=states_converged+1
             ! swap this state if necessary to make sure that the lowest
             ! n x_vecs of the set correspond to the n converged states
             if(icount>states_converged+1) then
               temp_jcount=pub_lr_tddft_num_conv_states+states_converged+1
               temp_icount=pub_lr_tddft_num_conv_states+icount
               do is=1, pub_num_spins
                  call sparse_embed_copy(temp_x(1,is),x_vec(temp_jcount,is))
                  call sparse_embed_copy(x_vec(temp_jcount,is), x_vec(temp_icount,is))
                  call sparse_embed_copy(x_vec(temp_icount,is),temp_x(1,is))
                  call sparse_embed_copy(temp_x_cov(1,is),x_vec_cov(temp_jcount,is))
                  call sparse_embed_copy(x_vec_cov(temp_jcount,is), &
                       x_vec_cov(temp_icount,is))
                  call sparse_embed_copy(x_vec_cov(temp_icount,is),temp_x_cov(1,is))
               enddo
               ! also swap eigenvalues in subspace_eval
               temp_result=subspace_eval(temp_jcount)
               subspace_eval(temp_jcount)=subspace_eval(temp_icount)
               subspace_eval(temp_icount)=temp_result
             endif
          endif
       enddo
    endif

    if(pub_on_root .and. pub_output_detail>=NORMAL) &
      write(stdout,'(a, i5, a)') '****** Number of newly converged &
        &states:', states_converged, '  *****'

    ! deallocate storage space:
    do is=1, pub_num_spins
       call sparse_embed_destroy(temp_transpose(is))
       do icount=1, subspace_size
          call sparse_embed_destroy(temp_x(icount,is))
          call sparse_embed_destroy(temp_x_cov(icount,is))
       enddo
    enddo
    deallocate(temp_x_cov,stat=ierr)
    call utils_dealloc_check('lr_tddft_check_convergence','temp_x_cov',ierr)
    deallocate(temp_transpose,stat=ierr)
    call utils_dealloc_check('lr_tddft_check_convergence','temp_transpose',ierr)
    deallocate(temp_x, stat=ierr)
    call utils_dealloc_check('lr_tddft_check_convergence','temp_x',ierr)

  end subroutine lr_tddft_check_convergence

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_cg(x_vec,lr_tddft_eval,mdl,hfxstate,response_dens,&
       sub_response_dens,fxc_fine,&
       cond_rep, cond_ngwf_basis, cond_denskern,kchc, &
       val_rep, val_ngwf_basis,val_denskern,hvkv,num_states,&
       ground_state_dens,sub_ground_state_dens)


    !=========================================================================!
    ! Subroutine does a cg optimisation of a batch of TDDFT trial vecs of the !
    ! Sternheimer equation. Essentially we solve for the minimum of the func- !
    ! tion (x*A*x)/(x*x), where x represents our batch of trial eigenvectors  !
    ! and A is the effective 2 particle TDDFT Hamiltonian                     !
    ! Modified for embedding by Joseph Prentice, July 2018                    !
    !=========================================================================!

    use comms, only: comms_barrier, pub_on_root
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use constants, only: stdout, NORMAL, BRIEF, VERBOSE
    use model_type, only: MODEL
    use hartree, only: hartree_via_multigrid
    use hf_exchange, only: HFX_STATE
    use is_solvation, only: have_initial_eps
    use linalg, only: linalg_dsygv_lt
    use linear_response, only: linear_response_operator
    use lr_tddft_utils, only: lr_tddft_utils_precond_batch,&
         lr_tddft_utils_penalty_batch
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_output_detail, &
         pub_lr_tddft_cg_threshold, &
         pub_lr_tddft_maxit_cg, pub_lr_tddft_maxit_pen, &
         pub_lr_tddft_penalty_tol, pub_lr_tddft_reset_cg, &
         pub_lr_tddft_write_kernels, pub_lr_tddft_preopt, &
         pub_lr_tddft_preopt_iter, pub_print_qc, pub_lr_tddft_num_conv_states,&
         pub_lr_tddft_check_conv_iter,pub_lr_tddft_precond, &
         pub_lr_tddft_sparse_region, pub_lr_tddft_restart, pub_debug_on_root, &
         pub_num_spins, pub_multigrid_hartree, pub_lr_optical_permittivity,&
         pub_is_bulk_permittivity
    use restart, only: restart_response_kernel_batch_write,&
         restart_response_kernel_batch_read
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_transpose, sparse_embed_trace, &
         sparse_embed_transpose_structure, sparse_embed_copy
    use sparse, only: sparse_show_matrix
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_qc_print, utils_banner

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    integer, intent(in) :: num_states  ! number of excitation energies
    type(SPAM3_EMBED), intent(inout) :: x_vec(num_states,pub_num_spins) ! solution to the
    ! equation. we are optimising several vectors simultaneously
    real(kind=DP), intent(inout) :: lr_tddft_eval(num_states) ! just the
    ! converged eval
    real(kind=DP), intent(inout) :: response_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_response_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)
    real(kind=DP), intent(in) :: fxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,2,mdl%nsub)
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis(:)
    type(NGWF_REP), intent(inout) :: cond_rep
    type(SPAM3_EMBED), intent(in) :: cond_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: kchc(pub_num_spins)
    type(NGWF_REP), intent(inout) :: val_rep
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type(SPAM3_EMBED), intent(in) :: val_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: hvkv(pub_num_spins)
    real(kind=DP), intent(inout) :: ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)

    ! Local Variables
    integer :: maxit
    real(kind=DP) :: tolerance, penalty_value
    real(kind=DP) :: energy_conv_states
    real(kind=DP), allocatable, dimension (:,:,:,:) :: ground_state_hartree
    real(kind=DP), allocatable, dimension (:,:) :: omega_conv
    real(kind=DP) :: gamma1, gamma2, gamma_current, gamma_prev, &
         temp_result, lambda_opt
    real(kind=DP) :: total_omega, total_omega_prev, rms_grad, rms_grad2, trace
    type(SPAM3_EMBED), allocatable, dimension (:,:) :: f_vec, f_prev, a_vec, &
         a_prev, y_vec, g_contravariant, L_vec, x_vec_cov
    type(SPAM3_EMBED), allocatable, dimension (:,:) :: g_vec, D_vec
    type(SPAM3_EMBED) :: fs_inv, ksv, ksc, scksv
    type(SPAM3_EMBED), allocatable, dimension(:) :: Svkv,kcSc,kvSv,Sckc
    type(SPAM3_EMBED), allocatable, dimension(:) :: temp_transpose
    type(SPAM3_EMBED) :: temp_kvSv
    type(SPAM3_EMBED) :: PSv, temp_Svkv, LSvkv
    real(kind=DP), allocatable, dimension(:,:) :: eval_mat, identity,&
         evec_mat
    real(kind=DP),allocatable,dimension(:) :: subspace_eval, subspace_eval_prev
    integer :: counter, reset_counter, subspace_size, subspace_conv_states
    integer :: ierr, is, isub, jsub
    integer :: icount, jcount
    logical :: converged, lr_tddft_precond

    call timer_clock('lr_tddft_cg',1)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entered tddft_cg &
         &routine'


    ! initalise values.
    lr_tddft_precond=pub_lr_tddft_precond
    if (pub_lr_tddft_preopt) then
       maxit=pub_lr_tddft_preopt_iter
       tolerance=pub_lr_tddft_cg_threshold
    else
       maxit=pub_lr_tddft_maxit_cg
       tolerance=pub_lr_tddft_cg_threshold
    endif
    subspace_size=num_states-pub_lr_tddft_num_conv_states
    energy_conv_states=0.0_DP

    ! allocate data structures
    if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
      allocate(ground_state_hartree(mdl%fine_grid%ld1,&
           mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins), &
           stat=ierr)
      call utils_alloc_check('lr_tddft_cg', 'ground_state_hartree', ierr)
    endif
    allocate(omega_conv(maxit+1,3), stat=ierr)
    call utils_alloc_check('lr_tddft_cg', 'omega_conv', ierr)
    omega_conv(:,:) = 0.0_DP ! jd: or valgrind complains
    allocate(x_vec_cov(num_states,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_cg','x_vec_cov',ierr)
    allocate(y_vec(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_cg', 'y_vec', ierr)
    allocate(L_vec(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_cg', 'L_vec', ierr)
    allocate(f_vec(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_cg', 'f_vec', ierr)
    allocate(f_prev(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_cg', 'f_prev', ierr)
    allocate(g_vec(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_cg', 'g_vec', ierr)
    allocate(a_vec(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_cg', 'a_vec', ierr)
    allocate(a_prev(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_cg', 'a_prev', ierr)
    allocate(D_vec(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_cg', 'D_vec', ierr)
    allocate(g_contravariant(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_cg','g_contravariant',ierr)
    allocate(kcSc(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_cg','kcSc',ierr)
    allocate(Sckc(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_cg','Sckc',ierr)
    allocate(kvSv(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_cg','kvSv',ierr)
    allocate(Svkv(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_cg','Svkv',ierr)
    allocate(temp_transpose(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_cg','temp_transpose',ierr)

    ! data structures for reduced eigenvalue problem
    allocate(eval_mat(subspace_size,subspace_size), stat=ierr)
    call utils_alloc_check('lr_tddft_cg','eval_mat',ierr)
    allocate(identity(subspace_size,subspace_size), stat=ierr)
    call utils_alloc_check('lr_tddft_cg','identity',ierr)
    allocate(evec_mat(subspace_size,subspace_size),stat=ierr)
    call utils_alloc_check('lr_tddft_cg','evec_mat', ierr)
    allocate(subspace_eval(subspace_size), stat=ierr)
    call utils_alloc_check('lr_tddft_cg','subspace_eval',ierr)
    allocate(subspace_eval_prev(subspace_size),stat=ierr)
    call utils_alloc_check('lr_tddft_cg','subspace_eval_prev',ierr)

    ! need to declare matrix in the correct size.
    ! PSv is the correct vector to transform from contra to
    ! covariant.
    call sparse_embed_create(PSv, x_vec(1,1), val_rep%overlap)

    do icount=1, num_states
       do is=1, pub_num_spins
          ! Set this to the correct structure as well.
          ! these vectors are covariant vectors
          L_vec(icount,is)%structure='TDRA'//trim(cond_rep%postfix)
          call sparse_embed_create(L_vec(icount,is), iscmplx=val_rep%overlap%p%iscmplx)
          call sparse_embed_create(f_vec(icount,is), cond_rep%overlap, PSv)
          call sparse_embed_create(f_prev(icount,is), f_vec(icount,is))
          call sparse_embed_create(y_vec(icount,is), f_vec(icount,is))
          call sparse_embed_create(x_vec_cov(icount,is),y_vec(icount,is))
          ! these vectors are contravariant vectors
          call sparse_embed_create(g_contravariant(icount,is),x_vec(icount,is))
          call sparse_embed_create(g_vec(icount,is), x_vec(icount,is))
          call sparse_embed_create(a_vec(icount,is), g_vec(icount,is))
          call sparse_embed_create(a_prev(icount,is), a_vec(icount,is))
          call sparse_embed_create(D_vec(icount,is), a_vec(icount,is))
          ! --> vectors have different sparsity patterns
       enddo
    enddo

    call sparse_embed_create(fs_inv, f_vec(1,1), val_rep%inv_overlap)
    ! vectors to enforce orthogonality
    do is=1, pub_num_spins
       call sparse_embed_create(kvSv(is), val_denskern(is), val_rep%overlap)
       call sparse_embed_create(Svkv(is), val_rep%overlap, val_denskern(is))
       call sparse_embed_create(kcSc(is), cond_denskern(is), cond_rep%overlap)
       call sparse_embed_create(Sckc(is), cond_rep%overlap, cond_denskern(is))
       call sparse_embed_transpose_structure(temp_transpose(is)%structure, x_vec(1,1))
       call sparse_embed_create(temp_transpose(is), iscmplx=x_vec(1,1)%p%iscmplx)
    enddo
    call sparse_embed_create(ksv, x_vec(1,1), val_rep%overlap)
    call sparse_embed_create(ksc, temp_transpose(1), cond_rep%overlap)
    call sparse_embed_create(scksv, cond_rep%overlap, ksv)
    call sparse_embed_create(temp_kvSv, x_vec(1,1), kvSv(1))
    call sparse_embed_create(temp_Svkv, x_vec(1,1), Svkv(1))
    call sparse_embed_create(LSvkv, L_vec(1,1),Svkv(1))

    ! construct vectors to enforce orthogonality
    do is=1, pub_num_spins
       call sparse_embed_product(kvSv(is), val_denskern(is), val_rep%overlap)
       call sparse_embed_product(Svkv(is), val_rep%overlap, val_denskern(is))
       call sparse_embed_product(kcSc(is), cond_denskern(is), cond_rep%overlap)
       call sparse_embed_product(Sckc(is), cond_rep%overlap, cond_denskern(is))
    enddo

    ! initialise variables
    converged = .false.
    gamma_prev=1.0_DP
    reset_counter=1
    penalty_value=0.0_DP

    ! OUTPUT HEADER
    if (pub_on_root .and. pub_output_detail >= BRIEF .and. .not. &
       pub_lr_tddft_preopt) then
       write(stdout,'(/a)') repeat('#',80)
       write(stdout,'(a)') utils_banner('#', 'LR-TDDFT CG optimisation')
       write(stdout,'(a)') repeat('#',80)
    end if

    if (pub_on_root .and. pub_output_detail >= BRIEF .and. &
         pub_lr_tddft_preopt) then
       write(stdout,'(/a)') repeat('#',80)
       write(stdout,'(a)') utils_banner('#', 'LR-TDDFT pre-optimisation')
       write(stdout,'(a)') repeat('#',80)
    endif

    ! brief output format
    if (pub_on_root .and. pub_output_detail==BRIEF) write(stdout,'(/a80)')&
       '|Iter     |   Energy (in Ha)  |  Delta E (in Ha) | Step Length | Penalty '

    ! create ground state hartree potential if required
    ! only required if TDDFT is performed within the implicit solvent
    ! model
    if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
      ! change bulk permittivity to have the correct value for optical
      ! frequencies
      pub_is_bulk_permittivity=pub_lr_optical_permittivity
      have_initial_eps=.false.

      if(pub_num_spins==1) then
         ground_state_dens=2.0_DP*ground_state_dens
      endif
      ! jcap: needs to be done over the whole system
      call hartree_via_multigrid(ground_state_hartree,ground_state_dens,&
           mdl%fine_grid,mdl%cell,elements=mdl%elements)
      if(pub_num_spins==1) then
         ground_state_dens=0.5_DP*ground_state_dens
      endif
    endif

    ! First project initial vectors into physical states
    do icount=1+pub_lr_tddft_num_conv_states, num_states
       do is=1, pub_num_spins
          call sparse_embed_copy(L_vec(icount,is),x_vec(icount,is))
          call sparse_embed_product(LSvkv, L_vec(icount,is), Svkv(is))
          call sparse_embed_product(x_vec(icount,is), kcSc(is), LSvkv)
       enddo
    enddo

    ! in order to consistently restart a calculation in case of sparsity, we
    ! need to READ IN x_vec instead of generating it from L_vec.
    if(pub_lr_tddft_sparse_region .and. pub_lr_tddft_restart) then
       if(pub_lr_tddft_num_conv_states>0) then
          call restart_response_kernel_batch_read(x_vec,&
               pub_lr_tddft_num_conv_states)
       else
          call restart_response_kernel_batch_read(x_vec,&
               num_states)
       endif
    endif

    ! compute x_vec_cov for all x_vecs:
    do icount=1, num_states
       do is=1, pub_num_spins
          call sparse_embed_product(PSv,x_vec(icount,is),val_rep%overlap)
          call sparse_embed_product(x_vec_cov(icount,is),cond_rep%overlap,PSv)
       enddo
    enddo

    call lr_tddft_utils_penalty_batch(penalty_value, x_vec, &
            x_vec_cov,cond_denskern, val_denskern, &
            cond_rep%overlap, val_rep%overlap,Svkv,kcSc,&
            pub_lr_tddft_num_conv_states,num_states)

    call lr_tddft_gram_schmidt(x_vec,x_vec_cov,num_states, cond_rep%overlap, &
         val_rep%overlap)

    ! Measure the penalty functional
    call lr_tddft_utils_penalty_batch(penalty_value, x_vec, &
            x_vec_cov,cond_denskern, val_denskern, &
            cond_rep%overlap, val_rep%overlap,Svkv,kcSc,&
            pub_lr_tddft_num_conv_states,num_states)

    counter=1
    ! iteration header for first iteration:
    if (pub_on_root .and. pub_output_detail >= NORMAL) then
       write(stdout,'(/a)') repeat('#',80)
       write(stdout,'(a,i4,a)')'######################### LR_TDDFT &
            &CG iteration ',counter,' ###########################'
       write(stdout,'(a,E12.4,a)') &
            '****** Initial penalty functional value of K^{1}: &
            &Q[K^{1}] = ', penalty_value, ' ******'
    endif

    ! HEADER for penalty func improvement
    if (penalty_value>pub_lr_tddft_penalty_tol) then
       if (pub_on_root.and.pub_output_detail>=VERBOSE) then
          write(stdout,'(a)') '******  Iterative improvement of &
               & K^{1}: ******'
          write(stdout, '(a)') '    |ITER|    Penalty value   &
               &|'
       endif
    endif

    ! Iteratively project into the correct states to improve the
    ! penalty functional value if neccessary.
    icount=1
    do while(icount<pub_lr_tddft_maxit_pen .and. &
         penalty_value>pub_lr_tddft_penalty_tol)
       ! reset penalty val
       penalty_value=0.0_DP
       do jcount=1+pub_lr_tddft_num_conv_states, num_states
          do is=1, pub_num_spins
             call sparse_embed_product(temp_Svkv, x_vec(jcount,is), Svkv(is))
             call sparse_embed_product(x_vec(jcount,is), kcSc(is), temp_Svkv)
          enddo
       enddo

       call lr_tddft_utils_penalty_batch(penalty_value, x_vec, &
            x_vec_cov,cond_denskern, val_denskern, &
            cond_rep%overlap, val_rep%overlap,Svkv,kcSc,&
            pub_lr_tddft_num_conv_states,num_states)

       ! OUTPUT iterative improvement
       if (pub_on_root.and.pub_output_detail >= VERBOSE) write(stdout,&
            '(i6,e12.4)') icount, penalty_value

       ! compute x_vec_cov for all x_vecs:
       do jcount=1, num_states
          do is=1, pub_num_spins
             call sparse_embed_product(PSv,x_vec(jcount,is),val_rep%overlap)
             call sparse_embed_product(x_vec_cov(jcount,is),cond_rep%overlap,PSv)
          enddo
       enddo

       icount=icount+1
    enddo

    ! check for convergence of iterative improvement
    if (penalty_value>pub_lr_tddft_penalty_tol) then
       if (pub_on_root.and.pub_output_detail >=NORMAL) &
            write(stdout,'(a,i4,a)')  'WARNING: maximum number of &
            &iterative improvements for K^{1} (',pub_lr_tddft_maxit_pen,&
            ') exceeded!'
    endif


    ! make sure that that all x_vecs are orthogonal to each other
    call lr_tddft_gram_schmidt(x_vec, x_vec_cov,num_states, cond_rep%overlap,&
         val_rep%overlap)

    if (pub_print_qc) then
       if(pub_on_root) write(stdout,'(a)') 'P_1 after second orthogonalisation:'
       do icount=1, num_states
          do isub=1,mdl%nsub
             do jsub=1,mdl%nsub
                if (pub_on_root .and. mdl%nsub .gt. 1) &
                     write(stdout,'(a,i0,i0)') 'Subregion ', isub, jsub
                call sparse_show_matrix(x_vec(icount,1)%m(isub,jsub),stdout)
             end do
          end do
       enddo
    endif

    ! Calculate y_vec first before calculating Omega
    ! Y=operator acting on x for all n_states x vectors
    ! Y is a covariant vector, g_noortho is its contravariant
    ! counterpart.
    do icount=1+pub_lr_tddft_num_conv_states, num_states
       if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
          call linear_response_operator(y_vec(icount,:), x_vec(icount,:),&
               cond_rep, cond_ngwf_basis, cond_denskern,kchc, val_rep, &
               val_ngwf_basis, val_denskern,hvkv,mdl,hfxstate,response_dens,&
               sub_response_dens, fxc_fine, ground_state_dens, &
               sub_ground_state_dens, g_contravariant(icount,:), &
               ground_state_hartree)
       else
          call linear_response_operator(y_vec(icount,:), x_vec(icount,:),&
               cond_rep, cond_ngwf_basis, cond_denskern, kchc, val_rep, &
               val_ngwf_basis, val_denskern,hvkv,mdl,hfxstate,response_dens,&
               sub_response_dens, fxc_fine, ground_state_dens, &
               sub_ground_state_dens, g_contravariant(icount,:))
       endif
    enddo

    ! Now initialise the set of inidvidual subspace eigenstates by dia-
    ! gonalising the subspace.
    ! Construct the matrix and diagonalise it explicitly
    do icount=1, subspace_size
       do is=1, pub_num_spins
          call sparse_embed_transpose(temp_transpose(is), x_vec(icount+&
                pub_lr_tddft_num_conv_states,is))
       enddo
       do jcount=icount, subspace_size
          temp_result=0.0_DP
          do is=1, pub_num_spins
             call sparse_embed_trace(trace, temp_transpose(is), &
                y_vec(jcount+pub_lr_tddft_num_conv_states,is))
             temp_result=temp_result+trace
          enddo
          eval_mat(icount,jcount)=temp_result
          eval_mat(jcount,icount)=temp_result
          if (icount==jcount) then
             identity(icount,icount)=1.0_DP
          else
             identity(icount,jcount)=0.0_DP
             identity(jcount,icount)=0.0_DP
          endif
       enddo
    enddo


    ! call internal eigensolve to diagonalise converged subspace
    evec_mat=eval_mat
    call linalg_dsygv_lt(evec_mat,subspace_eval_prev,identity,subspace_size)

    ! calculate initial total omega value.
    total_omega=0.0_DP
    do icount=1+pub_lr_tddft_num_conv_states, num_states
       call lr_tddft_calc_omega(omega_conv(2,1),x_vec(icount,:), &
            x_vec_cov(icount,:),y_vec(icount,:), cond_rep%overlap,&
            val_rep%overlap)
       total_omega=total_omega+omega_conv(2,1)
    enddo

    ! store current omega val
    omega_conv(2,1)=total_omega

    if (pub_output_detail >= NORMAL) then
       if (pub_on_root) then
          write(stdout,'(/11X,a)') repeat('.',56)
          write(stdout,'(a)')'           |     *** TDDFT optimisation &
               &not converged ***         |'
          write(stdout,'(a,f19.14,a)') &
               '           |       LR-TDDFT energy  =  ',total_omega,'     &
               &   |'
          write(stdout,'(11X,a)') repeat('=',56)
       endif
    endif

    ! store penalty value in convergence storage
    omega_conv(2,2)=penalty_value

    ! Set the first, initial omega(1)
    total_omega_prev=total_omega+0.1_DP
    ! Initialise the counter


    if (pub_on_root .and. pub_output_detail==BRIEF) then
       write(stdout,'(i16,f18.8,f18.8,f16.8,e15.4)') counter, total_omega, 1.0_DP, &
               1.0_DP,penalty_value
    endif

    counter=2


    ! MAIN LOOP
    do while(abs(total_omega-total_omega_prev)>tolerance.and.counter<maxit)
       ! PRINT iteration header
       if (pub_on_root .and. pub_output_detail >= NORMAL) then
          write(stdout,'(/a)') repeat('#',80)
          write(stdout,'(a,i4,a)')'######################### LR_TDDFT &
               &CG iteration ',counter,' ###########################'
          write(stdout,'(a,E12.4,a)') &
               '****** Initial penalty functional value of K^{1}: &
               &Q[K^{1}] = ', penalty_value, ' ******'
       endif

       ! First calculate f_vec. Make sure it is orthogonal to all current trial
       ! vecs. If we have any converged vecs, f_vec is orthogonal to all of
       ! those

       do icount=1+pub_lr_tddft_num_conv_states, num_states
          do is=1, pub_num_spins
             call sparse_embed_copy(f_vec(icount,is), y_vec(icount,is))
          enddo
          ! orthogonality to all unconverged x_vecs
          do jcount=1, num_states
             temp_result=0.0_DP
             do is=1, pub_num_spins
                call sparse_embed_transpose(temp_transpose(is), x_vec(jcount,is))
                call sparse_embed_trace(trace,temp_transpose(is),y_vec(icount,is))
                temp_result=temp_result+trace
             enddo

             do is=1, pub_num_spins
                call sparse_embed_axpy(f_vec(icount,is),x_vec_cov(jcount,is), &
                     -1.0_DP*temp_result)
             enddo
          enddo
       enddo

       ! g_vec is the tensorially correct version of f_vec. Need to multiply with
       ! inverse overlaps.
       if(lr_tddft_precond .and. .not. pub_lr_tddft_preopt) then
          call lr_tddft_utils_precond_batch(g_contravariant,f_vec,&
               num_states,pub_lr_tddft_num_conv_states,&
               response_dens,sub_response_dens,fxc_fine,&
               cond_rep,cond_ngwf_basis,cond_denskern,kchc,&
               val_rep,val_ngwf_basis,val_denskern,hvkv,&
               mdl,hfxstate,ground_state_dens,sub_ground_state_dens)
       endif

       do icount=1+pub_lr_tddft_num_conv_states, num_states
          do is=1, pub_num_spins
             call sparse_embed_copy(g_vec(icount,is), g_contravariant(icount,is))
          enddo

          ! now need to make sure that g is orthogonal to all current evecs
          do jcount=1, num_states

             temp_result=0.0_DP
             do is=1, pub_num_spins
                call sparse_embed_transpose(temp_transpose(is), g_contravariant(icount,is))
                call sparse_embed_trace(trace,temp_transpose(is),x_vec_cov(jcount,is))
                temp_result=temp_result+trace
             enddo

             do is=1, pub_num_spins
                call sparse_embed_axpy(g_vec(icount,is), x_vec(jcount,is), -1.0_DP*&
                     temp_result)
             enddo
          enddo
       enddo

       ! project g into physical space
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          do is=1, pub_num_spins
             call sparse_embed_product(temp_Svkv, g_vec(icount,is), Svkv(is))
             call sparse_embed_product(g_vec(icount,is), kcSc(is), temp_Svkv)
          enddo
       enddo

       ! calculate the rms gradient from g_vec.
       rms_grad=0.0_DP
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          temp_result=0.0_DP
          do is=1, pub_num_spins
             call sparse_embed_transpose(temp_transpose(is), g_vec(icount,is))
             call sparse_embed_product(ksv, g_vec(icount,is), val_rep%overlap)
             call sparse_embed_product(scksv, cond_rep%overlap, ksv)
             call sparse_embed_trace(trace,temp_transpose(is),scksv)
             temp_result=temp_result+trace
          enddo
          rms_grad=rms_grad+temp_result
       enddo
       temp_result=rms_grad/(1.0_DP*icount)
       rms_grad=sqrt(temp_result)

       ! Print RMS gradient
       if (pub_on_root .and. pub_output_detail >= NORMAL) then
          write(stdout,'(a,E12.4,a)') &
               '****** Calculated RMS gradient of g_vec. &
               &RMS grad = ', rms_grad, ' ******'
       endif

       ! calculate gamma1 and gamma2
       gamma1=0.0_DP
       gamma2=0.0_DP
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          temp_result=0.0_DP
          do is=1, pub_num_spins
             call sparse_embed_transpose(temp_transpose(is), g_vec(icount,is))
             call sparse_embed_trace(trace,temp_transpose(is),f_vec(icount,is))
             temp_result=temp_result+trace
          enddo
          gamma1=gamma1+temp_result

          ! temp transpose is still equal to g_mat_transpose
          ! f_prev is 0 in the first iteration, so leave it out.
          ! also in case of reset counter
          if (counter>2 .and. reset_counter<pub_lr_tddft_reset_cg+1) then
             temp_result=0.0_DP
             do is=1, pub_num_spins
                call sparse_embed_transpose(temp_transpose(is), g_vec(icount,is))
                call sparse_embed_trace(trace,temp_transpose(is),f_prev(icount,is))
                temp_result=temp_result+trace
             enddo
             gamma2=gamma2+temp_result
          else
             ! reset gamma_prev
             gamma_prev=1.0_DP
          endif
       enddo

       if (pub_on_root .and. pub_output_detail >= VERBOSE) then
          write(stdout,'(a,f12.10)') 'Polak-Ribiere coeff. gamma1 &
               &=     ',gamma1
          write(stdout,'(a,f12.10)') 'Polak-Ribiere coeff. gamma2 &
               &=     ',gamma2
       endif


       !calculate the current gamma value
       gamma_current=(gamma1-gamma2)/gamma_prev

       if (pub_on_root .and. pub_output_detail >= NORMAL) then
          write(stdout,'(a,f12.10)') 'Conjugate gradient coefficient&
               &=     ',gamma_current
       endif

       ! Calculate A vec. A= -G+gamma*A_prev. Again, only necessary
       ! if we are not in the first iteration. Otherwise, A_prev=0
       ! A_prev also equal to 0 in case we reset search direction
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          do is=1, pub_num_spins
             call sparse_embed_copy(a_vec(icount,is), g_vec(icount,is))
             call sparse_embed_scale(a_vec(icount,is), -1.0_DP)
          enddo
          if (counter>2 .and. reset_counter<pub_lr_tddft_reset_cg+1) then
             do is=1, pub_num_spins
                call sparse_embed_axpy(a_vec(icount,is), a_prev(icount,is), &
                     gamma_current)
             enddo
          endif
       enddo

       ! now only need to make A vec orthogonal to x again==> D
       ! calculate D_mat:
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          do is=1, pub_num_spins
             call sparse_embed_copy(D_vec(icount,is), a_vec(icount,is))
          enddo
          do jcount=1, num_states
             temp_result=0.0_DP
             do is=1, pub_num_spins
                call sparse_embed_transpose(temp_transpose(is),a_vec(icount,is))
                call sparse_embed_trace(trace,temp_transpose(is),x_vec_cov(jcount,is))
                temp_result=temp_result+trace
             enddo

             do is=1, pub_num_spins
                call sparse_embed_axpy(D_vec(icount,is), x_vec(jcount,is), &
                     -1.0_DP*temp_result)
             enddo
          enddo
       enddo

       ! calculate the rms gradient from D_vec.
       rms_grad2=0.0_DP
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          temp_result=0.0_DP
          do is=1, pub_num_spins
             call sparse_embed_transpose(temp_transpose(is), D_vec(icount,is))
             call sparse_embed_product(ksv, D_vec(icount,is), val_rep%overlap)
             call sparse_embed_product(scksv, cond_rep%overlap, ksv)
             call sparse_embed_trace(trace,temp_transpose(is),scksv)
             temp_result=temp_result+trace
          enddo
          rms_grad2=rms_grad2+temp_result
       enddo
       temp_result=rms_grad2/(1.0_DP*(num_states-pub_lr_tddft_num_conv_states))
       rms_grad2=sqrt(temp_result)

       ! Print RMS gradient
       if (pub_on_root .and. pub_output_detail >= NORMAL) then
          write(stdout,'(a,E12.4,a)') &
               '***** Calculated RMS gradient of search direction. &
               &RMS grad = ', rms_grad2, ' *****'
       endif

       ! Find optimum search direction, ie lambda optimum by minimizing the energy
       if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
          call lr_tddft_line_min(lambda_opt, x_vec,x_vec_cov, D_vec, y_vec, &
               mdl,hfxstate,response_dens,sub_response_dens,fxc_fine,&
               cond_ngwf_basis,cond_rep, cond_denskern,kchc,kcSc, &
               val_ngwf_basis, val_rep, val_denskern, hvkv,Svkv,&
               num_states,energy_conv_states,L_vec,ground_state_dens,&
               sub_ground_state_dens,ground_state_hartree)
       else
          call lr_tddft_line_min(lambda_opt, x_vec, x_vec_cov,D_vec, y_vec, &
               mdl,hfxstate,response_dens,sub_response_dens,fxc_fine,&
               cond_ngwf_basis,cond_rep,cond_denskern,kchc,kcSc, &
               val_ngwf_basis, val_rep, val_denskern, hvkv,Svkv,&
               num_states,energy_conv_states,L_vec,ground_state_dens,&
               sub_ground_state_dens)
       endif

       omega_conv(counter+1,3)=lambda_opt

       ! Print step size
       if (pub_on_root .and. pub_output_detail >= NORMAL) then
          write(stdout,'(a,E12.4,a)') &
               '****** Performed line minimisation. Ideal step size&
               & Lambda = ', lambda_opt, ' ******'
       endif


       ! recalculate pentalty value for new set of x_vec
       penalty_value=0.0_DP
       call lr_tddft_utils_penalty_batch(penalty_value, x_vec, &
            x_vec_cov,cond_denskern, val_denskern, &
            cond_rep%overlap, val_rep%overlap,Svkv,kcSc,&
            pub_lr_tddft_num_conv_states,num_states)

       ! HEADER for penalty func improvement
       if (penalty_value>pub_lr_tddft_penalty_tol) then
          if (pub_on_root.and.pub_output_detail>=VERBOSE) then
             write(stdout,'(a)') '******  Iterative improvement of &
                  & K^{1}: ******'
             write(stdout, '(a)') '     |ITER|    Penalty value   &
                  &|'
          endif
       endif


       ! Iteratively project into the correct states to improve the
       ! penalty functional value if neccessary.
       icount=1
       do while(icount<pub_lr_tddft_maxit_pen .and. &
            penalty_value>pub_lr_tddft_penalty_tol)
          ! reset penalty val
          penalty_value=0.0_DP
          do jcount=1+pub_lr_tddft_num_conv_states, num_states
             do is=1, pub_num_spins
                call sparse_embed_product(temp_Svkv, x_vec(jcount,is), Svkv(is))
                call sparse_embed_product(x_vec(jcount,is), kcSc(is), temp_Svkv)
             enddo
          enddo

          call lr_tddft_utils_penalty_batch(penalty_value, x_vec, &
                 x_vec_cov,cond_denskern, val_denskern, &
                 cond_rep%overlap, val_rep%overlap,Svkv,kcSc,&
                 pub_lr_tddft_num_conv_states,num_states)

          ! OUTPUT iterative improvement
          if (pub_on_root.and.pub_output_detail >= VERBOSE) &
               write(stdout,'(i6,e12.4)') icount, penalty_value
          ! compute x_vec_cov for all x_vecs:
          do jcount=1, num_states
             do is=1, pub_num_spins
               call sparse_embed_product(PSv,x_vec(jcount,is),val_rep%overlap)
               call sparse_embed_product(x_vec_cov(jcount,is),cond_rep%overlap,PSv)
             enddo
          enddo

          icount=icount+1
       enddo

       ! check for convergence of iterative improvement
       if (penalty_value>pub_lr_tddft_penalty_tol) then
          if (pub_on_root.and.pub_output_detail >=NORMAL) &
               write(stdout,'(a,i4,a)')  'WARNING: maximum number of &
               &iterative improvements for K^{1}(',pub_lr_tddft_maxit_pen,&
               ') exceeded!'
       endif


       ! orthogonalise and recalculate penalty value
       call lr_tddft_gram_schmidt(x_vec, x_vec_cov,num_states, cond_rep%overlap,&
            val_rep%overlap)


       penalty_value=0.0_DP
       call lr_tddft_utils_penalty_batch(penalty_value, x_vec, &
            x_vec_cov,cond_denskern, val_denskern, &
            cond_rep%overlap, val_rep%overlap,Svkv,kcSc,&
            pub_lr_tddft_num_conv_states,num_states)

       ! write out the tddft response kernels if required
       if (pub_lr_tddft_write_kernels) then
          call restart_response_kernel_batch_write(x_vec, num_states)
       endif

       ! need to recalculate y_vec
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
             call linear_response_operator(y_vec(icount,:), x_vec(icount,:),&
                  cond_rep, cond_ngwf_basis, cond_denskern, kchc, val_rep, &
                  val_ngwf_basis, val_denskern,hvkv,mdl,hfxstate,response_dens,&
                  sub_response_dens,fxc_fine, ground_state_dens,&
                  sub_ground_state_dens,g_contravariant(icount,:), &
                  ground_state_hartree)
          else
             call linear_response_operator(y_vec(icount,:), x_vec(icount,:),&
                  cond_rep, cond_ngwf_basis, cond_denskern,kchc, val_rep, &
                  val_ngwf_basis, val_denskern,hvkv,mdl,hfxstate,response_dens,&
                  sub_response_dens,fxc_fine, ground_state_dens,&
                  sub_ground_state_dens,g_contravariant(icount,:))
          endif
       enddo


       total_omega_prev=total_omega
       total_omega=0.0_DP
       ! Calculate new omega
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          call lr_tddft_calc_omega(omega_conv(counter+1,1),x_vec(icount,:),&
               x_vec_cov(icount,:),y_vec(icount,:), cond_rep%overlap, &
               val_rep%overlap)
          total_omega=total_omega+omega_conv(counter+1,1)
       enddo
       total_omega=total_omega+energy_conv_states

       ! store final omega:
       omega_conv(counter+1,1)=total_omega

       ! Set all the correct previous values for gamma_prev, a_prev and f_prev
       gamma_prev=gamma1
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          do is=1, pub_num_spins
             call sparse_embed_copy(a_prev(icount,is), a_vec(icount,is))
             call sparse_embed_copy(f_prev(icount,is), f_vec(icount,is))
          enddo
       enddo

       omega_conv(counter+1,2)=penalty_value

       ! BRIEF output
       if (pub_on_root .and. pub_output_detail==BRIEF) then
          write(stdout,'(i16,f18.8,f18.8,f16.8,e15.4)') counter, total_omega, &
             total_omega-total_omega_prev,lambda_opt,penalty_value
       endif


       ! OUTPUT convergence information
       if (pub_output_detail >= NORMAL) then
          if (pub_on_root) then
             write(stdout,'(/a)')'           ............................&
                  &............................'
             if (abs(total_omega-total_omega_prev)>tolerance) then
                write(stdout,'(a)')'           |     *** TDDFT optimisation&
                     & not converged ***         |'
             else
                write(stdout,'(a)')'           |     *** TDDFT optimisation&
                     & converged ***             |'
             endif
             write(stdout,'(a,f19.14,a)') &
                  '           |       LR-TDDFT energy  =  ',total_omega,'     &
                  &   |'
             write(stdout,'(a,f19.14,a)') &
                  '           |       Change in omega  =  ',total_omega-&
                  total_omega_prev,'        |'
             write(stdout,'(a,f19.14,a)') &
                  '           |       RMS gradient     =  ',rms_grad,'        |'
             if (abs(total_omega-total_omega_prev)>tolerance) then
                write(stdout,'(a)')'           | --> Total change in omega &
                     &higher than threshold      |'
             else
                write(stdout,'(a)')'           | --> Total change in omega &
                     &lower than threshold       |'
             endif
             write(stdout,'(a)')'           ============================&
                  &============================'
          endif
       endif

       ! Check for convergence of individual eigenstates in the subset
       ! only check every pub_lr_tddft_num_conv_iter iterations
       if(mod(counter,pub_lr_tddft_check_conv_iter)==0) then
          call lr_tddft_check_convergence(x_vec,x_vec_cov,y_vec,subspace_eval,&
               subspace_eval_prev, evec_mat, eval_mat,identity,&
               subspace_conv_states,energy_conv_states,num_states)

          if(subspace_conv_states==0) then
            do icount=1, subspace_size
               subspace_eval_prev(icount)=subspace_eval(icount)
            enddo
          else
            ! some states have actually converged
            ! set reset counter to the right value so that the cg algorithm
            ! gets reset the next iteration
            reset_counter=pub_lr_tddft_reset_cg

            pub_lr_tddft_num_conv_states=pub_lr_tddft_num_conv_states+&
               subspace_conv_states
            subspace_size=subspace_size-subspace_conv_states
            deallocate(subspace_eval_prev,stat=ierr)
            call utils_dealloc_check('lr_tddft_cg','subspace_eval_prev',ierr)
            deallocate(evec_mat,stat=ierr)
            call utils_dealloc_check('lr_tddft_cg','evec_mat',ierr)
            deallocate(eval_mat,stat=ierr)
            call utils_dealloc_check('lr_tddft_cg','eval_mat',ierr)
            deallocate(identity,stat=ierr)
            call utils_dealloc_check('lr_tddft_cg','identity',ierr)

            ! reallocate matrices at the correct size
            allocate(subspace_eval_prev(subspace_size),stat=ierr)
            call utils_alloc_check('lr_tddft_cg','subspace_eval_prev',ierr)
            allocate(evec_mat(subspace_size,subspace_size),stat=ierr)
            call utils_alloc_check('lr_tddft_cg','evec_mat',ierr)
            allocate(eval_mat(subspace_size,subspace_size),stat=ierr)
            call utils_alloc_check('lr_tddft_cg','eval_mat',ierr)
            allocate(identity(subspace_size,subspace_size),stat=ierr)
            call utils_alloc_check('lr_tddft_cg','identity',ierr)

            ! set subspace_eval_prev and reallocate subspace_eval at correct
            ! size
            do icount=1, subspace_size
               subspace_eval_prev(icount)=subspace_eval(icount+subspace_conv_states)
            enddo

            deallocate(subspace_eval,stat=ierr)
            call utils_dealloc_check('lr_tddft_cg','subspace_eval_prev',ierr)
            allocate(subspace_eval(subspace_size),stat=ierr)
            call utils_alloc_check('lr_tddft_cg','subspace_eval_prev',ierr)

            ! do a gram-schmidt and generate new y_vecs for the set
            ! orthogonalise and recalculate penalty value
            call lr_tddft_gram_schmidt(x_vec, x_vec_cov,num_states, cond_rep%overlap,&
                val_rep%overlap)


            penalty_value=0.0_DP
            call lr_tddft_utils_penalty_batch(penalty_value, x_vec, &
                x_vec_cov,cond_denskern, val_denskern, &
                cond_rep%overlap, val_rep%overlap,Svkv,kcSc,&
                pub_lr_tddft_num_conv_states,num_states)

            ! write out the tddft response kernels if required
            if (pub_lr_tddft_write_kernels) then
               call restart_response_kernel_batch_write(x_vec, num_states)
            endif

            ! need to recalculate y_vec
            do icount=1+pub_lr_tddft_num_conv_states, num_states
               if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
                  call linear_response_operator(y_vec(icount,:), x_vec(icount,:),&
                       cond_rep, cond_ngwf_basis, cond_denskern,kchc, val_rep, &
                       val_ngwf_basis, val_denskern,hvkv,mdl,hfxstate,&
                       response_dens,sub_response_dens,fxc_fine,&
                       ground_state_dens,sub_ground_state_dens,&
                       g_contravariant(icount,:),ground_state_hartree)
               else
                  call linear_response_operator(y_vec(icount,:), x_vec(icount,:),&
                       cond_rep, cond_ngwf_basis, cond_denskern, kchc, val_rep, &
                       val_ngwf_basis, val_denskern,hvkv,mdl,hfxstate,&
                       response_dens,sub_response_dens,fxc_fine,&
                       ground_state_dens,sub_ground_state_dens,&
                       g_contravariant(icount,:))
               endif
            enddo
          endif
       endif

       ! end check


       ! include QC
       if (pub_on_root .and. pub_print_qc .and. .not. pub_lr_tddft_preopt) then
          call utils_qc_print('Omega', total_omega)
          call utils_qc_print('delta_Omega',total_omega-&
               total_omega_prev)
          call utils_qc_print('TDDFT_grad', rms_grad)
       endif

       !in case of reset, set reset counter back to 1
       if (reset_counter>pub_lr_tddft_reset_cg) then
          reset_counter=1
       endif


       ! increase counters
       counter=counter+1
       reset_counter=reset_counter+1
    enddo

    ! End of main loop
    ! check for convergence
    if (abs(total_omega-total_omega_prev)>tolerance) then
       if (pub_on_root) write(stdout,'(a,i4,a)') &
            'WARNING: maximum number of LR-TDDFT CG iterations(',maxit, &
            ') exceeded!'
    endif

    ! Print summary sheet
    if (pub_on_root .and. .not. pub_output_detail==BRIEF) &
       write(stdout,'(/a80)') '|ITER|    Total Energy   &
         &|     Penalty value      |     step       '
    do icount=1, counter-1
       if (pub_on_root .and. .not. pub_output_detail==BRIEF) then
            write(stdout,'(i16,f16.10,e22.10,f21.10)') icount,&
            omega_conv(icount+1,1),omega_conv(icount+1,2),&
            omega_conv(icount+1,3)
       end if
    enddo


    ! now construct final omega
    do icount=1+pub_lr_tddft_num_conv_states, num_states
       call lr_tddft_calc_omega(lr_tddft_eval(icount), x_vec(icount,:), &
            x_vec_cov(icount,:),y_vec(icount,:), cond_rep%overlap, val_rep%overlap)
    enddo

    ! if some states are already converged, create appropriate y_vec for them
    if(pub_lr_tddft_num_conv_states>0) then
       do icount=1, pub_lr_tddft_num_conv_states
          if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
             call linear_response_operator(y_vec(icount,:), x_vec(icount,:),&
                  cond_rep, cond_ngwf_basis, cond_denskern, kchc, val_rep, &
                  val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,response_dens,&
                  sub_response_dens,fxc_fine, ground_state_dens,&
                  sub_ground_state_dens,g_contravariant(icount,:), &
                  ground_state_hartree)
          else
             call linear_response_operator(y_vec(icount,:), x_vec(icount,:),&
                  cond_rep, cond_ngwf_basis, cond_denskern, kchc, val_rep, &
                  val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,response_dens,&
                  sub_response_dens,fxc_fine, ground_state_dens,&
                  sub_ground_state_dens,g_contravariant(icount,:))
          endif
       enddo
    endif

    ! reallocate matrices at correct size
    deallocate(evec_mat,stat=ierr)
    call utils_dealloc_check('lr_tddft_cg','evec_mat',ierr)
    deallocate(eval_mat,stat=ierr)
    call utils_dealloc_check('lr_tddft_cg','eval_mat',ierr)
    deallocate(identity,stat=ierr)
    call utils_dealloc_check('lr_tddft_cg','identity',ierr)

    ! reallocate matrices at the correct size
    allocate(evec_mat(num_states,num_states),stat=ierr)
    call utils_alloc_check('lr_tddft_cg','evec_mat',ierr)
    allocate(eval_mat(num_states,num_states),stat=ierr)
    call utils_alloc_check('lr_tddft_cg','eval_mat',ierr)
    allocate(identity(num_states,num_states),stat=ierr)
    call utils_alloc_check('lr_tddft_cg','identity',ierr)


    ! Construct the matrix and diagonalise it explicitly
    do icount=1, num_states
       do is=1, pub_num_spins
          call sparse_embed_transpose(temp_transpose(is), x_vec(icount,is))
       enddo
       do jcount=icount, num_states
          temp_result=0.0_DP
          do is=1, pub_num_spins
             call sparse_embed_trace(trace,temp_transpose(is),y_vec(jcount,is))
             temp_result=temp_result+trace
          enddo
          eval_mat(icount,jcount)=temp_result
          eval_mat(jcount,icount)=temp_result
          if (icount==jcount) then
             identity(icount,icount)=1.0_DP
          else
             identity(icount,jcount)=0.0_DP
             identity(jcount,icount)=0.0_DP
          endif
       enddo
    enddo

    ! call internal eigensolve to diagonalise converged subspace
    evec_mat=eval_mat
    call linalg_dsygv_lt(evec_mat,lr_tddft_eval,identity,num_states)

    ! now need to construct 'real' k from a superposition of Ks, determined through eval mat.
    ! for this purpose, use D_vec, which we dont need any more, as a temporary storage of K mats.
    ! need to use D_vec instead of y_vec since D_vec has the same
    ! sparsity pattern as x_vec
    do icount=1, num_states
       do is=1, pub_num_spins
          call sparse_embed_copy(D_vec(icount,is), x_vec(icount,is))
       enddo
    enddo


    ! now construct correct response density matrices
    ! seems to be done correctly
    do icount=1, num_states
       do jcount=1, num_states
          if (jcount==1) then
             do is=1,pub_num_spins
                call sparse_embed_copy(x_vec(icount,is), D_vec(jcount,is))
                temp_result=evec_mat(jcount,icount)
                call sparse_embed_scale(x_vec(icount,is), temp_result)
             enddo
          else
             temp_result=evec_mat(jcount,icount)
             do is=1, pub_num_spins
                call sparse_embed_axpy(x_vec(icount,is), D_vec(jcount,is), temp_result)
             enddo
          endif
       enddo
    enddo

    ! write out the final tddft response kernels if required
    if (pub_lr_tddft_write_kernels) then
       call restart_response_kernel_batch_write(x_vec, num_states)
    endif

    ! deallocate storage...
    deallocate(omega_conv, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg', 'omega_conv', ierr)

    call sparse_embed_destroy(fs_inv)
    call sparse_embed_destroy(ksv)
    do is=1, pub_num_spins
       call sparse_embed_destroy(temp_transpose(is))
       call sparse_embed_destroy(Sckc(is))
       call sparse_embed_destroy(kcSc(is))
       call sparse_embed_destroy(Svkv(is))
       call sparse_embed_destroy(kvSv(is))
       do icount=1, num_states
          call sparse_embed_destroy(x_vec_cov(icount,is))
          call sparse_embed_destroy(g_contravariant(icount,is))
          call sparse_embed_destroy(y_vec(icount,is))
          call sparse_embed_destroy(L_vec(icount,is))
          call sparse_embed_destroy(f_vec(icount,is))
          call sparse_embed_destroy(f_prev(icount,is))
          call sparse_embed_destroy(g_vec(icount,is))
          call sparse_embed_destroy(a_vec(icount,is))
          call sparse_embed_destroy(a_prev(icount,is))
       enddo
    enddo
    call sparse_embed_destroy(ksc)
    call sparse_embed_destroy(scksv)

    deallocate(x_vec_cov,stat=ierr)
    call utils_dealloc_check('lr_tddft_cg','x_vec_cov',ierr)
    deallocate(g_contravariant, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg', 'g_contravariant',ierr)
    deallocate(y_vec, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg', 'y_vec', ierr)
    deallocate(L_vec, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg', 'L_vec', ierr)
    deallocate(f_vec, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg', 'f_vec', ierr)
    deallocate(f_prev, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg', 'f_prev', ierr)
    deallocate(g_vec, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg', 'g_vec', ierr)
    deallocate(a_vec, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg', 'a_vec', ierr)
    deallocate(a_prev, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg', 'a_prev', ierr)
    deallocate(D_vec, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg', 'D_vec', ierr)
    deallocate(eval_mat, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg','eval_mat',ierr)
    deallocate(identity, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg','identity',ierr)
    deallocate(evec_mat,stat=ierr)
    call utils_dealloc_check('lr_tddft_cg','evec_mat', ierr)
    deallocate(subspace_eval_prev,stat=ierr)
    call utils_dealloc_check('lr_tddft_cg','subspace_eval_prev',ierr)
    if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
      deallocate(ground_state_hartree, stat=ierr)
      call utils_dealloc_check('lr_tddft_cg', 'ground_state_hartree',ierr)
    endif
    deallocate(kcSc, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg','kcSc',ierr)
    deallocate(Sckc, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg','Sckc',ierr)
    deallocate(kvSv, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg','kvSv',ierr)
    deallocate(Svkv, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg','Svkv',ierr)
    deallocate(temp_transpose, stat=ierr)
    call utils_dealloc_check('lr_tddft_cg','temp_transpose',ierr)


    call timer_clock('lr_tddft_cg',2)


  end subroutine lr_tddft_cg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_calculate(full_evecs,lr_tddft_evals, response_dens, &
       sub_response_dens, fxc_fine, &
       cond_ngwf_basis, cond_rep, cond_denskern, kchc, cond_ham, cond_ngwf_num, &
       val_ngwf_basis, val_rep, val_denskern, hvkv, val_ham, val_ngwf_num, &
       proj_basis, nl_projectors, mdl, hfxstate, val_evecs, cond_evecs, &
       ground_state_dens, sub_ground_state_dens, joint_evals)

    !=====================================================================!
    ! Subroutine successively calls lr_tddft_cg to calculate the first N  !
    ! Excitation energies of the system. Converged eigenvalues are stored !
    ! in an array of matrices and passed down to properties calculations  !
    ! Modified for embedding by Joseph Prentice, July 2018                !
    !=====================================================================!

    use comms, only: comms_barrier, pub_on_root
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
         dense_put_element
    use function_basis, only: FUNC_BASIS
    use constants, only: stdout
    use hf_exchange, only: HFX_STATE
    use lr_tddft_utils, only: lr_tddft_utils_init_tddft_vecs,&
         lr_tddft_utils_calc_mlwfs
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use projectors, only: projector_set
    use qnto, only: qnto_calculate
    use restart, only: restart_response_kernel_batch_read
    use rundat, only: pub_lr_tddft_num_states, pub_lr_tddft_restart, &
         pub_lr_tddft_preopt, pub_qnto_ref_dir, &
         pub_lr_tddft_analysis, pub_print_qc, &
         pub_lr_tddft_num_conv_states, pub_num_spins, &
         pub_lr_tddft_subsystem_coupling,pub_lr_tddft_init_random,&
         pub_lr_tddft_joint_set, pub_lr_tddft_mlwf_analysis
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_transpose, sparse_embed_transpose_structure, &
         sparse_embed_copy
    use sparse, only: sparse_show_matrix
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_qc_print, utils_rand_gen

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: full_evecs(pub_lr_tddft_num_states,pub_num_spins)
    real(kind=DP), intent(inout) :: lr_tddft_evals(pub_lr_tddft_num_states)
    type(MODEL), intent(in) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(:)
    real(kind=DP), intent(inout) :: response_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_response_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)
    real(kind=DP), intent(in) :: fxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,2,mdl%nsub)
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis(:)
    type(NGWF_REP), intent(inout) :: cond_rep
    type(SPAM3_EMBED), intent(in) :: cond_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: kchc(pub_num_spins)
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type(NGWF_REP), intent(inout) :: val_rep
    type(SPAM3_EMBED), intent(in) :: val_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: hvkv(pub_num_spins)
    type(DEM), intent(inout) :: val_evecs(pub_num_spins)
    type(DEM), intent(inout) :: cond_evecs(pub_num_spins)
    integer, intent(in) :: val_ngwf_num, cond_ngwf_num
    real(kind=DP), intent(inout) :: ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)
    type(SPAM3_EMBED) :: cond_ham(pub_num_spins)
    type(SPAM3_EMBED) :: val_ham(pub_num_spins)
    real(kind=DP), optional, intent(in) :: joint_evals(cond_ngwf_num,pub_num_spins)

    ! Local Variables
    integer :: jat_orig, kat_orig, loc_jngwf, loc_kngwf, jngwf, kngwf, &
               jngwf_prev, kngwf_prev
    integer :: icount, ierr, jcount, kcount, num_states, &
               loc_num_conv_states, is, isub, jsub, ksub
    real(kind=DP) :: temp_result
    type(SPAM3_EMBED) :: L_vec    ! auxillary kernel
    type(DEM) :: initial_guess_vec
    real(kind=DP), allocatable, dimension(:) :: rand_array

    num_states=pub_lr_tddft_num_states
    ! allocate data storage for converged_evecs
    L_vec%structure= 'TDRA'//trim(cond_rep%postfix)
    call sparse_embed_create(L_vec, iscmplx=full_evecs(1,1)%p%iscmplx)

    call dense_create(initial_guess_vec, cond_ngwf_num, &
         val_ngwf_num, iscmplx=full_evecs(1,1)%p%iscmplx)


    ! check for restart option
    if (pub_lr_tddft_restart .or. pub_lr_tddft_subsystem_coupling) then
       call restart_response_kernel_batch_read(full_evecs,&
            num_states)
    else
       ! are there already converged states? read them in
       if(pub_lr_tddft_num_conv_states>0) then
          call restart_response_kernel_batch_read(full_evecs,&
               pub_lr_tddft_num_conv_states)
       endif

       ! initialise to KS transitions or random vecs
       if(.not. pub_lr_tddft_init_random .and. pub_lr_tddft_num_conv_states==0) then
         call lr_tddft_utils_init_tddft_vecs(full_evecs,cond_ngwf_basis,&
             val_ngwf_basis,cond_rep,val_rep,cond_ham,val_ham,val_denskern,&
             val_ngwf_num,cond_ngwf_num,mdl,num_states)

      else
          allocate(rand_array(cond_ngwf_num*val_ngwf_num),&
             stat=ierr)
          call utils_alloc_check('lr_tddft_calculate','rand_array',ierr)

          ! generate all other states randomly.
          do icount=1+pub_lr_tddft_num_conv_states, num_states
             ! create an array of random numbers for this kernel
             call utils_rand_gen(rand_array,cond_ngwf_num*&
                val_ngwf_num, icount)

             ! ndmh: new version - uses input file order for random array
             jcount = 0  ! jcap: counts through NGWFs in the order we consider them, for the random number generator
             jngwf_prev = 0  ! jcap: keeps track of number of NGWFs in regions previously dealt with
             ! jcap: loop over regions
             do jsub=1,mdl%nsub
                do jat_orig=1,mdl%regions(jsub)%par%nat
                   do loc_jngwf=1,cond_ngwf_basis(jsub)%num_on_atom(mdl%regions(jsub)%par%distr_atom( &
                        jat_orig))
                      jcount = jcount + 1
                      jngwf = cond_ngwf_basis(jsub)%first_on_atom(mdl%regions(jsub)%par%distr_atom( &
                           jat_orig)) + loc_jngwf - 1
                      kcount = 0  ! jcap: Same as jcount and jngwf_prev above
                      kngwf_prev = 0
                      do ksub=1,mdl%nsub
                         do kat_orig=1,mdl%regions(ksub)%par%nat
                            do loc_kngwf=1,val_ngwf_basis(ksub)%num_on_atom( &
                                 mdl%regions(ksub)%par%distr_atom(kat_orig))
                               kngwf = val_ngwf_basis(ksub)%first_on_atom( &
                                    mdl%regions(ksub)%par%distr_atom(kat_orig)) + loc_kngwf - 1
                               kcount = kcount + 1
                               temp_result = &
                                    rand_array((jcount-1)*val_ngwf_num+kcount)
                               call dense_put_element(temp_result, initial_guess_vec, &
                                    jngwf_prev+jngwf, kngwf_prev+kngwf)
                            end do ! loc_kngwf
                         end do ! kat_orig
                         ! jcap: Add on the total number of functions in this region
                         kngwf_prev = kngwf_prev + val_ngwf_basis(ksub)%num
                      end do ! ksub
                   end do ! loc_jngwf
                end do ! jat_orig
                ! jcap: Add on the total number of functions in this region
                jngwf_prev = jngwf_prev + cond_ngwf_basis(jsub)%num
             end do ! jsub

             ! loop over elements in matrix
             ! old version - did not use input file order
             !do jcount=1, cond_ngwf_num
             !   do kcount=1, val_ngwf_num
             !      temp_result=rand_array((jcount-1)*val_ngwf_num+kcount)
             !      call dense_put_element(temp_result, initial_guess_vec, &
             !          jcount,kcount)
             !   enddo
             !enddo

             call dense_convert(L_vec, initial_guess_vec)
             do is=1, pub_num_spins
                call sparse_embed_copy(full_evecs(icount,is), L_vec)
             enddo
          enddo
          deallocate(rand_array,stat=ierr)
          call utils_dealloc_check('lr_tddft_calculate','rand_array',ierr)
       endif
    endif

    ! print initial response density kernel for QC testing
    if (pub_print_qc) then
       if (pub_on_root) write(stdout,'(a)') 'Initial P_1 matrix:'
       do icount=1, num_states
          ! jcap: loop over regions
          do isub=1,mdl%nsub
             do jsub=1,mdl%nsub
                if (pub_on_root .and. mdl%nsub .gt. 1) &
                     write(stdout,'(a,i0,i0)') 'Subregion ', isub, jsub
                call sparse_show_matrix(full_evecs(icount,1)%m(isub,jsub), stdout)
             end do
          end do
       enddo
       ! extra debugging info in QC mode to help diagnose issues around
       ! sign-flips in kernel
       if (pub_on_root) then
          write(stdout,'(a)') 'Radial NGWFs:'
          do isub=1,mdl%nsub
             do icount=1,mdl%regions(isub)%par%num_species
                do jcount=1,mdl%regions(isub)%radial_ngwfs(icount,1)%nshells
                   if (mdl%nsub .gt. 1) then
                      write(stdout,'(a,i3,a,i3,a,i3,a,3f15.8)') 'Region ',isub,', species ',icount, &
                           ', shell ',jcount,': ', &
                           mdl%regions(isub)%radial_ngwfs(icount,1)%func_real(1:3,jcount)
                   else
                      write(stdout,'(a,i3,a,i3,a,3f15.8)') 'Species ',icount, &
                           ', shell ',jcount,': ', &
                           mdl%regions(isub)%radial_ngwfs(icount,1)%func_real(1:3,jcount)
                   end if
                end do
             end do
          end do
       end if
       do jsub=1,mdl%nsub
          do isub=1,mdl%nsub
             if (pub_on_root .and. mdl%nsub .gt. 1) &
                  write(stdout,'(a,i0,i0)') 'Subregion ', isub, jsub
             if (pub_on_root) write(stdout,'(a)') 'P_v matrix:'
             call sparse_show_matrix(val_denskern(1)%m(isub,jsub), stdout, show_elems=.true.)
             if (pub_on_root) write(stdout,'(a)') 'P_c matrix:'
             call sparse_show_matrix(cond_denskern(1)%m(isub,jsub), stdout, show_elems=.true.)
          end do
       end do
    endif

    ! None of the rest of the routine is relevant for QNTO calculations
    if (pub_qnto_ref_dir=='') then

    if (.not. pub_lr_tddft_subsystem_coupling) then
       if (pub_lr_tddft_preopt) then
          ! first do preoptimisation.
          loc_num_conv_states=pub_lr_tddft_num_conv_states
          call lr_tddft_cg(full_evecs,lr_tddft_evals,mdl,hfxstate,&
               response_dens,sub_response_dens,fxc_fine,&
               cond_rep, cond_ngwf_basis, cond_denskern, kchc, val_rep,&
               val_ngwf_basis, val_denskern, hvkv, num_states,&
               ground_state_dens, sub_ground_state_dens)
          pub_lr_tddft_num_conv_states=loc_num_conv_states
          pub_lr_tddft_preopt=.false.

          ! now do real optimisation
          call lr_tddft_cg(full_evecs,lr_tddft_evals,mdl,hfxstate,&
               response_dens,sub_response_dens,fxc_fine,&
               cond_rep, cond_ngwf_basis, cond_denskern, kchc, val_rep,&
               val_ngwf_basis, val_denskern, hvkv, num_states,&
               ground_state_dens, sub_ground_state_dens)
       else
          ! now call a the TDDFT cg algorithm.
          call lr_tddft_cg(full_evecs,lr_tddft_evals,mdl,hfxstate,&
               response_dens,sub_response_dens,fxc_fine,&
               cond_rep, cond_ngwf_basis, cond_denskern, kchc, val_rep,&
               val_ngwf_basis, val_denskern, hvkv, num_states,&
               ground_state_dens, sub_ground_state_dens)
       endif
    else
       call lr_tddft_couplings(full_evecs,lr_tddft_evals,mdl,hfxstate,&
            response_dens,sub_response_dens, fxc_fine,&
            cond_rep, cond_ngwf_basis, cond_denskern, kchc, val_rep, &
            val_ngwf_basis, val_denskern, hvkv, num_states, &
            ground_state_dens, sub_ground_state_dens)
    endif

    ! Perform some properties calculations. The momentum rep for computing
    ! oscillator strengths can only be used in the joint rep and if lr_tddft_analysis
    ! is set to T, as it requires Kohn-Sham eigenstates
    if(pub_lr_tddft_joint_set .and. pub_lr_tddft_analysis) then
       call lr_tddft_properties(full_evecs, cond_rep, cond_ngwf_basis,&
            val_rep, val_ngwf_basis, proj_basis,nl_projectors,response_dens,&
            num_states, lr_tddft_evals, mdl, cond_evecs,joint_evals)
    else
       call lr_tddft_properties(full_evecs,cond_rep, cond_ngwf_basis,&
            val_rep, val_ngwf_basis,proj_basis,nl_projectors,response_dens,&
            num_states, lr_tddft_evals, mdl)
    endif

    ! This scales very unfavourably with system size (This has a full
    ! O(N^3) scaling). Only do if requested.
    if (pub_lr_tddft_analysis) then
       call lr_tddft_analysis(full_evecs, cond_rep, &
            cond_ngwf_basis, val_rep, val_ngwf_basis, val_evecs,&
            cond_evecs, num_states, val_ngwf_num, cond_ngwf_num)

       if(pub_lr_tddft_mlwf_analysis) then
          ! construct maximally localised wannier functions
          call lr_tddft_utils_calc_mlwfs(full_evecs,cond_ngwf_basis,&
               val_ngwf_basis,cond_rep,val_rep, mdl, proj_basis,&
               lr_tddft_evals,num_states)
       endif

    endif

    ! print QC information
    if (pub_on_root .and. pub_print_qc) then
       call utils_qc_print('omega_1', &
            lr_tddft_evals(1))
    endif

    endif

    call sparse_embed_destroy(L_vec)
    call dense_destroy(initial_guess_vec)

  end subroutine lr_tddft_calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_couplings(x_vec,lr_tddft_evals,mdl,hfxstate,&
       response_dens,sub_response_dens,fxc_fine,&
       cond_rep, cond_ngwf_basis, cond_denskern, kchc,&
       val_rep,val_ngwf_basis,val_denskern,hvkv,num_states,&
       ground_state_dens,sub_ground_state_dens)

     !=======================================================================!
     ! Subroutine taking the converged density kernels from different sub-   !
     ! system TDDFT calculations and calculating the couplings between them  !
     ! in order to construct the effective eigenstates of the entire system  !
     ! Modified for embedding by Joseph Prentice, July 2018                  !
     !=======================================================================!
    use comms, only: pub_on_root
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use constants, only: stdout
    use model_type, only: MODEL
    use hartree, only: hartree_via_multigrid
    use hf_exchange, only: HFX_STATE
    use is_solvation, only: have_initial_eps
    use linalg, only: linalg_dsygv_lt
    use linear_response, only: linear_response_operator
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_lr_tddft_write_kernels, pub_lr_tddft_preopt, &
         pub_num_spins, pub_multigrid_hartree, pub_is_bulk_permittivity,&
         pub_lr_optical_permittivity
    use restart, only: restart_response_kernel_batch_write
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_transpose, sparse_embed_trace, &
         sparse_embed_transpose_structure, sparse_embed_copy
    use sparse, only: sparse_show_matrix
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    integer, intent(in) :: num_states  ! number of excitation energies
    type(SPAM3_EMBED), intent(inout) :: x_vec(num_states,pub_num_spins) ! solution to the
    ! equation. we are optimising several vectors simultaneously
    real(kind=DP), intent(inout) :: lr_tddft_evals(num_states) ! just the
    ! converged eval
    real(kind=DP), intent(inout) :: response_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_response_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)
    real(kind=DP), intent(in) :: fxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,2,mdl%nsub)
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis(:)
    type(NGWF_REP), intent(inout) :: cond_rep
    type(SPAM3_EMBED), intent(in) :: cond_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: kchc(pub_num_spins)
    type(NGWF_REP), intent(inout) :: val_rep
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type(SPAM3_EMBED), intent(in) :: val_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: hvkv(pub_num_spins)
    real(kind=DP), intent(inout) :: ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)


    ! local variables
    integer :: icount, jcount, ierr, is, isub, jsub
    type(SPAM3_EMBED), allocatable, dimension(:,:) :: y_vec, temp_vec
    type(SPAM3_EMBED) :: PSv, temp_trans
    type(SPAM3_EMBED), allocatable, dimension(:) :: ScPSv
    real(kind=DP), allocatable, dimension(:,:,:,:) :: ground_state_hartree
    real(kind=DP) :: temp_val1, temp_val2, trace
    real(kind=DP), allocatable, dimension (:,:) :: coupling_matrix, evecs,&
      overlap_mat

    ! allocate appropriate variables
    allocate(ScPSv(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_couplings','ScPSv',ierr)
    allocate(y_vec(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_couplings','y_vec',ierr)
    allocate(temp_vec(num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_couplings','temp_vec',ierr)
    allocate(coupling_matrix(num_states,num_states), stat=ierr)
    call utils_alloc_check('lr_tddft_couplings','coupling_matrix',ierr)
    allocate(evecs(num_states,num_states), stat=ierr)
    call utils_alloc_check('lr_tddft_couplings','evecs',ierr)
    allocate(overlap_mat(num_states,num_states), stat=ierr)
    call utils_alloc_check('lr_tddft_couplings','overlap_mat',ierr)
    if(pub_multigrid_hartree) then
       allocate(ground_state_hartree(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
         mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
       call utils_alloc_check('lr_tddft_couplings','ground_state_hartree',&
         ierr)
    endif

    ! allocate sparse matrix structures
    call sparse_embed_create(PSv, x_vec(1,1), val_rep%overlap)
    do is=1, pub_num_spins
       call sparse_embed_create(ScPSv(is),cond_rep%overlap,PSv)
    enddo
    call sparse_embed_transpose_structure(temp_trans%structure,x_vec(1,1))
    call sparse_embed_create(temp_trans, iscmplx=x_vec(1,1)%p%iscmplx)
    do icount=1, num_states
       do is=1, pub_num_spins
          call sparse_embed_create(y_vec(icount,is), ScPSv(is))
          call sparse_embed_create(temp_vec(icount,is), x_vec(1,1))
       enddo
    enddo

    ! First determine if ground state hartree potential needs to be
    ! evaluated
    if (pub_multigrid_hartree) then
       pub_is_bulk_permittivity=pub_lr_optical_permittivity
       have_initial_eps=.false.

       if(pub_num_spins==1) then
          ground_state_dens=2.0_DP*ground_state_dens
       endif
       ! jcap: needs to be done over the whole system
       call hartree_via_multigrid(ground_state_hartree,ground_state_dens,&
            mdl%fine_grid,mdl%cell, elements=mdl%elements)
       if(pub_num_spins==1) then
          ground_state_dens=0.5_DP*ground_state_dens
       endif
    endif

    ! calculate y_vec for every single x_vec
    do icount=1, num_states
       if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
          call linear_response_operator(y_vec(icount,:), x_vec(icount,:),&
               cond_rep, cond_ngwf_basis, cond_denskern, kchc, val_rep, &
               val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,response_dens,&
               sub_response_dens,fxc_fine, ground_state_dens, &
               sub_ground_state_dens,ground_state_hartree=ground_state_hartree)
       else
          call linear_response_operator(y_vec(icount,:), x_vec(icount,:),&
               cond_rep, cond_ngwf_basis, cond_denskern, kchc, val_rep, &
               val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,response_dens,&
               sub_response_dens,fxc_fine, ground_state_dens,&
               sub_ground_state_dens)
       endif
    enddo

    ! Now build the coupling and the overlap matrix
    do icount=1, num_states
       do is=1, pub_num_spins
          call sparse_embed_product(PSv, x_vec(icount,is), val_rep%overlap)
          call sparse_embed_product(ScPSv(is), cond_rep%overlap, PSv)
       enddo
       do jcount=icount, num_states
          temp_val1=0.0_DP
          temp_val2=0.0_DP
          do is=1, pub_num_spins
             call sparse_embed_transpose(temp_trans, x_vec(jcount,is))
             call sparse_embed_trace(trace,temp_trans,ScPsv(is))
             temp_val1=temp_val1+trace
             call sparse_embed_trace(trace,temp_trans,y_vec(icount,is))
             temp_val2=temp_val2+trace
          enddo
          overlap_mat(icount,jcount)=temp_val1
          overlap_mat(jcount,icount)=temp_val1
          coupling_matrix(icount,jcount)=temp_val2
          coupling_matrix(jcount,icount)=temp_val2
       enddo
    enddo

    evecs=coupling_matrix

    ! solve the matrix problem in order to gain the couplings between all
    ! excited states of the system
    call linalg_dsygv_lt(evecs,lr_tddft_evals,overlap_mat,num_states)

    ! now need to construct 'real' K from a superposition of Ks, determined through eval mat.
    ! sparsity pattern as x_vec
    do icount=1, num_states
       do is=1, pub_num_spins
          call sparse_embed_copy(temp_vec(icount,is), x_vec(icount,is))
       enddo
    enddo


    ! now construct correct response density matrices
    do icount=1, num_states
       do jcount=1, num_states
          do is=1, pub_num_spins
             if (jcount==1) then
                call sparse_embed_copy(x_vec(icount,is), temp_vec(jcount,is))
                temp_val1=evecs(jcount,icount)
                call sparse_embed_scale(x_vec(icount,is), temp_val1)
             else
                temp_val1=evecs(jcount,icount)
                call sparse_embed_axpy(x_vec(icount,is), temp_vec(jcount,is), temp_val1)
             endif
          enddo
       enddo
    enddo

    ! write them out if required
    if (pub_lr_tddft_write_kernels) then
       call restart_response_kernel_batch_write(x_vec, num_states)
    endif

    ! in principle we are done now. However, before leaving the routine, do
    ! some extra analysis, where the response density matrices of the coupled
    ! system are analysed in terms of the response density matrices of the
    ! uncoupled systems. All the information required is stored within the
    ! eigenvector matrix of the subysystem diagonalisation.
    if(pub_on_root) write(stdout, '(a)') '********  SUBSYSTEM DIAGONALISATION &
             &********'
    if(pub_on_root) write(stdout, '(a)') 'Analysis of global transitions in terms &
       &of original subsystem states:'

    do icount=1, num_states
       if (pub_on_root) write(stdout, '(a,i4,a)') '******Transition: ',&
            icount, '*****************'
       do jcount=1, num_states
          if(abs(evecs(jcount,icount))>0.001) then
             if (pub_on_root) write(stdout,'(a,i4,a,f12.8)') &
                     'Original transition ', jcount, '   =   ',&
                     evecs(jcount,icount)*evecs(jcount,icount)
          endif
       enddo
    enddo

    ! deallocate sparse_matrix structure
    call sparse_embed_destroy(PSv)
    do is=1, pub_num_spins
       call sparse_embed_destroy(ScPSv(is))
    enddo
    call sparse_embed_destroy(temp_trans)
    do icount=1, num_states
       do is=1, pub_num_spins
          call sparse_embed_destroy(y_vec(icount,is))
          call sparse_embed_destroy(temp_vec(icount,is))
       enddo
    enddo

    ! deallocate appropriate variables
    deallocate(ScPSv, stat=ierr)
    call utils_dealloc_check('lr_tddft_couplings','ScPSv',ierr)
    deallocate(y_vec, stat=ierr)
    call utils_dealloc_check('lr_tddft_couplings','y_vec',ierr)
    deallocate(temp_vec, stat=ierr)
    call utils_dealloc_check('lr_tddft_couplings','temp_vec',ierr)
    deallocate(coupling_matrix, stat=ierr)
    call utils_dealloc_check('lr_tddft_couplings','coupling_matrix',ierr)
    deallocate(evecs, stat=ierr)
    call utils_dealloc_check('lr_tddft_couplings','evecs',ierr)
    deallocate(overlap_mat, stat=ierr)
    call utils_dealloc_check('lr_tddft_couplings','overlap_mat',ierr)
    if(pub_multigrid_hartree) then
       deallocate(ground_state_hartree, stat=ierr)
       call utils_dealloc_check('lr_tddft_couplings','ground_state_hartree',&
         ierr)
    endif

  end subroutine lr_tddft_couplings

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_analysis(response_kernel, cond_rep, &
       cond_ngwf_basis, val_rep, val_ngwf_basis, val_evecs,&
       cond_evecs, num_states, val_ngwf_num, cond_ngwf_num)

    !========================================================================!
    ! Subroutine calculates some properties for individual excitations (ie   !
    ! Which transitions make up which response vectors. Also, and very im -  !
    ! portantly: How much percentage of 'forbidden' transitions is contained !
    ! per excitation energy?                                                 !
    ! Modified for embedding by Joseph Prentice, July 2018                   !
    !========================================================================!

    use comms, only: comms_barrier, pub_on_root
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
         dense_put_element,dense_get_element
    use function_basis, only: FUNC_BASIS
    use constants, only: stdout
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_lr_tddft_joint_set,pub_lr_tddft_HOMO_num,&
         pub_lr_tddft_LUMO_num, pub_num_kpoints, PUB_1K, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_transpose, sparse_embed_transpose_structure
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: num_states
    type(SPAM3_EMBED), intent(inout) :: response_kernel(num_states,pub_num_spins)
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis(:)
    type(NGWF_REP), intent(in) :: cond_rep
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type(NGWF_REP), intent(in) :: val_rep
    type(DEM), intent(inout) :: val_evecs(pub_num_spins)
    type(DEM), intent(inout) :: cond_evecs(pub_num_spins)
    integer, intent(in) :: val_ngwf_num, cond_ngwf_num

    ! Local Variables
    type(SPAM3_EMBED):: ksv, scksv, val_sparse, cond_sparse, cond_trans,&
        P1_evec,evec_P1_evec
    type(DEM) :: kernel_dens
    integer :: is, ierr
    integer :: icount, jcount, omega_count
    integer, allocatable, dimension(:) :: val_limit, cond_limit, cond_start
    real(kind=DP) :: temp_sum, temp_check_sum, check_sum_printed

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine lr_tddft_analysis not ready yet for more&
         & than one k-point.')

    allocate(val_limit(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_analysis','val_limit',ierr)
    allocate(cond_limit(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_analysis','cond_limit',ierr)
    allocate(cond_start(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_analysis','cond_start',ierr)


    ! allocate data structures
    call sparse_embed_create(ksv, response_kernel(1,1), val_rep%overlap)
    call sparse_embed_create(scksv,cond_rep%overlap,ksv)
    call sparse_embed_create(cond_sparse,cond_rep%inv_overlap)
    call sparse_embed_create(val_sparse,val_rep%inv_overlap)
    call sparse_embed_create(cond_trans, cond_rep%inv_overlap)
    call sparse_embed_create(P1_evec, scksv, val_sparse)
    call sparse_embed_create(evec_P1_evec, cond_trans, P1_evec)
    call dense_create(kernel_dens,cond_ngwf_num,val_ngwf_num)


    ! First check allowed transitions..
    do is=1, pub_num_spins
       if (pub_lr_tddft_joint_set) then
          cond_start(is)=val_rep%n_occ(is,PUB_1K)+1
          if (val_rep%n_occ(is,PUB_1K)<pub_lr_tddft_HOMO_num) then
             val_limit(is)=val_rep%n_occ(is,PUB_1K)
          else
             val_limit(is)=pub_lr_tddft_HOMO_num
          endif
          if (cond_rep%n_occ(is,PUB_1K)<pub_lr_tddft_LUMO_num) then
             cond_limit(is)=cond_start(is)+cond_rep%n_occ(is,PUB_1K)-1
          else
             cond_limit(is)=cond_start(is)+pub_lr_tddft_LUMO_num-1
          endif
       else
          cond_start(is)=1
          if (val_rep%n_occ(is,PUB_1K)<pub_lr_tddft_HOMO_num) then
             val_limit(is)=val_rep%n_occ(is,PUB_1K)
          else
             val_limit(is)=pub_lr_tddft_HOMO_num
          endif
          if (cond_rep%n_occ(is,PUB_1K)<pub_lr_tddft_LUMO_num) then
             cond_limit(is)=cond_rep%n_occ(is,PUB_1K)
          else
             cond_limit(is)=pub_lr_tddft_LUMO_num
          endif
       endif
    enddo

    do omega_count=1, num_states
       if (pub_on_root) write(stdout, '(a,i4,a)') '******Transition: ',&
            omega_count, '*****************'
       temp_check_sum=0.0_DP
       check_sum_printed=0.0_DP

       do is=1, pub_num_spins
          if(pub_num_spins==2) then
             if(is==1) then
                if (pub_on_root) write(stdout, '(a,i4,a)') '*********  SPIN UP&
                 &    **************'
             else
                if (pub_on_root) write(stdout, '(a,i4,a)') '*********  SPIN DOWN&
                 &  **************'
             endif
          endif
          ! create effective dense response matrix in cond-val evec space
          call sparse_embed_product(ksv, response_kernel(omega_count,is),&
               val_rep%overlap)
          call sparse_embed_product(scksv,cond_rep%overlap,ksv)
          call dense_convert(val_sparse,val_evecs(is))
          call dense_convert(cond_sparse,cond_evecs(is))
          call sparse_embed_transpose(cond_trans,cond_sparse)
          call sparse_embed_product(P1_evec,scksv,val_sparse)
          call sparse_embed_product(evec_P1_evec, cond_trans, P1_evec)
          call dense_convert(kernel_dens,evec_P1_evec)


          ! loop over allowed valence states
          ! we limit ourselves to the most important states
          do icount=val_rep%n_occ(is,PUB_1K)-val_limit(is)+1, &
                    val_rep%n_occ(is,PUB_1K)
             do jcount=cond_start(is), cond_limit(is)
                ! get overlap between vectors
                temp_sum=0.0_DP
                call dense_get_element(temp_sum, kernel_dens,jcount,icount)

                ! square vector element
                temp_sum=temp_sum*temp_sum

                if (temp_sum>0.001_DP) then
                   if (pub_on_root) write(stdout,'(a,i4,a,i4,a,f12.8)') &
                        'KS transition ', icount, '  ----> ',&
                        jcount+val_rep%n_occ(is,PUB_1K)+1-cond_start(is), &
                        '  = ',temp_sum
                   check_sum_printed=check_sum_printed+temp_sum
                endif
                temp_check_sum=temp_check_sum+temp_sum
             enddo
          enddo
        enddo  ! end spin loop
        if (pub_on_root) write(stdout, '(a, i4, a, f15.12)') &
               'Total percentage of transition ',omega_count,&
            '  printed above: ', check_sum_printed


       if (pub_on_root) write(stdout, '(a, i4, a, f15.12)') &
            'Total percentage of transition ',omega_count,&
            '  in allowed space: ', temp_check_sum
    enddo

    ! deallocate all data structures:
    deallocate(val_limit)
    call utils_dealloc_check('lr_tddft_analysis','val_limit',ierr)
    deallocate(cond_start)
    call utils_dealloc_check('lr_tddft_analysis','cond_start',ierr)
    deallocate(cond_limit)
    call utils_dealloc_check('lr_tddft_analysis','cond_limit',ierr)
    call sparse_embed_destroy(ksv)
    call sparse_embed_destroy(scksv)
    call sparse_embed_destroy(cond_sparse)
    call sparse_embed_destroy(val_sparse)
    call sparse_embed_destroy(cond_trans)
    call sparse_embed_destroy(P1_evec)
    call sparse_embed_destroy(evec_P1_evec)
    call dense_destroy(kernel_dens)

  end subroutine lr_tddft_analysis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_properties(response_kernel, cond_rep, &
       cond_ngwf_basis,val_rep, val_ngwf_basis, proj_basis,&
       nl_projectors,response_dens, num_states, energies, mdl,&
       joint_evecs,joint_evals)

    !=========================================================================!
    ! Subroutine plotting the response density of all optimised states and    !
    ! outputs them in a chosen format. Furthermore, the routine calculates    !
    ! and plots projections of the response density into conduction and       !
    ! valence states.                                                         !
    ! Modified for embedding by Joseph Prentice, July 2018                    !
    !=========================================================================!

    use comms, only: comms_barrier, pub_on_root
    use dense, only: DEM
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use constants, only: stdout, ANGSTROM, FINE_STRUCTURE, &
         HARTREE_IN_NS
    use lr_tddft_utils, only: lr_tddft_utils_oscillator
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_lr_tddft_write_densities, &
         pub_lr_tddft_triplet, pub_num_spins, pub_lr_tddft_spectrum_smear
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_transpose, sparse_embed_transpose_structure, sparse_embed_copy
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_destroy
    use utils, only: utils_alloc_check, utils_dealloc_check
    use visual, only: visual_scalarfield

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    real(kind=DP), intent(inout) :: response_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis(:)
    type(NGWF_REP), intent(in) :: cond_rep
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(:)
    type(NGWF_REP), intent(in) :: val_rep
    integer, intent(in) :: num_states ! the number of states that
                                      ! are explicitly optimised
    type(SPAM3_EMBED), intent(inout) :: response_kernel(num_states,pub_num_spins)
    real(kind=DP), intent(inout) :: energies(num_states)
    type(DEM), optional, intent(in) :: joint_evecs(pub_num_spins)
    real(kind=DP), optional,intent(in) :: joint_evals(:,:)

    ! Local Variables
    type(SPAM3_EMBED) :: temp_vec
    type(SPAM3_EMBED) :: Ktrans, KSv, ScK
    type(SPAM3_EMBED), allocatable, dimension(:) :: val_kernel, cond_kernel
    type(SPAM3), allocatable, dimension(:) :: eff_kernel
    integer :: icount, ierr, is
    character(50) :: filename, tempfile
    real(kind=DP), allocatable, dimension(:) :: oscillator
    real(kind=DP), allocatable, dimension(:) :: lifetime

    ! jcap: embedding local variables
    integer :: isub,jsub
    real(kind=DP), allocatable, dimension(:,:,:,:) :: sub_response_dens

    ! allocate data space
    call sparse_embed_create(temp_vec, response_kernel(1,1))
    allocate(eff_kernel(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_properties', 'eff_kernel', ierr)
    allocate(val_kernel(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_properties', 'val_kernel', ierr)
    allocate(cond_kernel(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_properties', 'cond_kernel', ierr)

    do is=1, pub_num_spins
       val_kernel(is)%structure='K'//val_rep%postfix
       call sparse_embed_create(val_kernel(is), &
            iscmplx=response_kernel(1,is)%p%iscmplx)
       cond_kernel(is)%structure='K'//cond_rep%postfix
       call sparse_embed_create(cond_kernel(is), &
            iscmplx=response_kernel(1,is)%p%iscmplx)
    enddo
    call sparse_embed_transpose_structure(Ktrans%structure,&
         response_kernel(1,1))
    call sparse_embed_create(Ktrans, iscmplx=response_kernel(1,1)%p%iscmplx)
    call sparse_embed_create(KSv,response_kernel(1,1),val_rep%overlap)
    call sparse_embed_create(ScK,cond_rep%overlap,response_kernel(1,1))
    allocate(oscillator(num_states), stat=ierr)
    call utils_alloc_check('lr_tddft_properties', 'oscillator', ierr)
    allocate(lifetime(num_states), stat=ierr)
    call utils_alloc_check('lr_tddft_properties', 'lifetime', ierr)


    ! For each excitation energy, calculate its oscillator strength
    if (pub_lr_tddft_triplet) then
       if (pub_on_root) write(stdout,'(/a80)')'|Excitation|    Energy &
            &(in Ha)   |     Oscillator Str.  | Lifetime (in ns)   '

       do icount=1, num_states
          if (pub_on_root) write(stdout,'(i12,f18.8,e21.5,e21.5)') &
               icount, energies(icount), 0.0_DP, 0.0_DP
       enddo
    else

       if(present(joint_evecs) .and. present(joint_evals)) then
          call lr_tddft_utils_oscillator(oscillator, energies, &
               cond_ngwf_basis, cond_rep,val_ngwf_basis, val_rep,&
               proj_basis,nl_projectors,response_kernel,num_states,mdl,&
               joint_evecs,joint_evals)
       else
          call lr_tddft_utils_oscillator(oscillator, energies, &
               cond_ngwf_basis, cond_rep, val_ngwf_basis, val_rep, proj_basis,&
               nl_projectors,response_kernel, num_states, mdl)
       endif

       ! HEADER for excitation energies and oscillator strengths
       if (pub_on_root) write(stdout,'(/a80)')'|Excitation|    Energy &
            &(in Ha)   |     Oscillator Str.  | Lifetime (in ns)   '

       do icount=1, num_states
          lifetime(icount)=((1.0_DP/FINE_STRUCTURE)**3)/&
               (2.0_DP*energies(icount)*energies(icount)*oscillator(icount))
          if (pub_on_root) write(stdout,'(i12,f18.8,e21.5,e21.5)') &
               icount, energies(icount), oscillator(icount), &
               lifetime(icount)*HARTREE_IN_NS
       enddo
    endif

    ! plot spectrum if required
    if (pub_lr_tddft_spectrum_smear > 0.0_DP) then
       call lr_tddft_spectrum(energies, oscillator, lifetime, num_states)
    end if

    ! now print transition densities. Also, print density as projected
    ! into conduction and into valence space.

    if(pub_lr_tddft_write_densities) then

       allocate(sub_response_dens(mdl%fine_grid%ld1,&
            mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('lr_tddft_properties', 'sub_response_dens', ierr)

       do icount=1, num_states
          ! jcap: loop over regions and sum up densities
          do isub=1,mdl%nsub
             do jsub=1,mdl%nsub

                sub_response_dens=0.d0

                do is=1, pub_num_spins
                   call sparse_create(eff_kernel(is), &
                        response_kernel(icount,is)%m(isub,jsub))
                   call sparse_copy(eff_kernel(is), &
                        response_kernel(icount,is)%m(isub,jsub))
                enddo

                call density_on_grid(sub_response_dens, mdl%fine_grid, mdl%dbl_grid, &
                     mdl%cell, mdl%fftbox, eff_kernel, &
                     cond_rep%ngwf_cross_overlap_tr%m(isub,jsub), &
                     cond_rep%ngwfs_on_grid(isub), cond_ngwf_basis(isub), &
                     val_rep%ngwfs_on_grid(jsub), val_ngwf_basis(jsub))

                response_dens=response_dens+sub_response_dens

                do is=1, pub_num_spins
                   call sparse_destroy(eff_kernel(is))
                enddo
             enddo
          enddo


          write(filename,*) icount
          write(tempfile,*) '_TDDFT_response_dens_'
          write(filename,*) trim(adjustl(tempfile))//trim(adjustl(filename))

          if(pub_num_spins==2) then
             response_dens(:,:,:,1)=response_dens(:,:,:,1)+&
               response_dens(:,:,:,2)
          endif

          call visual_scalarfield(response_dens(:,:,:,1),mdl%fine_grid, &
               mdl%cell,'e-h response density (in e/ang^3) for:', filename, &
               mdl%elements, ANGSTROM**3)
       enddo
    endif

    ! transition densities are not very meaningfull. Should also print
    ! the projected transition densities into valence and conduction
    ! states --> Effective electron and hole densities
    if (pub_lr_tddft_write_densities) then
       do icount=1, num_states
          do is=1, pub_num_spins
             call sparse_embed_transpose(Ktrans, response_kernel(icount,is))
             ! start with valence kernel
             call sparse_embed_product(ScK,cond_rep%overlap,&
                  response_kernel(icount,is))
             call sparse_embed_product(val_kernel(is),Ktrans,ScK)
          enddo

          ! jcap: loop over regions and sum up densities
          do isub=1,mdl%nsub
             do jsub=1,mdl%nsub

                sub_response_dens=0.d0

                do is=1, pub_num_spins
                   call sparse_create(eff_kernel(is), val_kernel(is)%m(isub,jsub))
                   call sparse_copy(eff_kernel(is), val_kernel(is)%m(isub,jsub))
                enddo

                ! calculate the density
                call density_on_grid(sub_response_dens, mdl%fine_grid, mdl%dbl_grid, &
                     mdl%cell, mdl%fftbox, eff_kernel,&
                     val_rep%ngwf_overlap%m(isub,jsub), val_rep%ngwfs_on_grid(isub),&
                     val_ngwf_basis(isub), val_rep%ngwfs_on_grid(jsub), &
                     val_ngwf_basis(jsub))

                response_dens=response_dens+sub_response_dens

                do is=1, pub_num_spins
                   call sparse_destroy(eff_kernel(is))
                enddo

             end do
          end do

          ! set filename
          write(filename,*) icount
          write(tempfile,*) '_TDDFT_hole_dens_'
          write(filename,*) trim(adjustl(tempfile))//trim(adjustl(filename))

          if(pub_num_spins==2) then
             response_dens(:,:,:,1)=response_dens(:,:,:,1)+&
               response_dens(:,:,:,2)
          endif

          call visual_scalarfield(response_dens(:,:,:,1),mdl%fine_grid, &
               mdl%cell,'e-h response density (in e/ang^3) for:', filename, &
               mdl%elements, ANGSTROM**3)

          ! now do the same for the electron density
          do is=1, pub_num_spins
             call sparse_embed_transpose(ktrans,response_kernel(icount,is))
             call sparse_embed_product(KSv,response_kernel(icount,is),&
                  val_rep%overlap)
             call sparse_embed_product(cond_kernel(is),KSv,Ktrans)
          enddo

          ! jcap: loop over regions and sum up densities
          do isub=1,mdl%nsub
             do jsub=1,mdl%nsub

                sub_response_dens=0.d0

                do is=1, pub_num_spins
                   call sparse_create(eff_kernel(is), cond_kernel(is)%m(isub,jsub))
                   call sparse_copy(eff_kernel(is), cond_kernel(is)%m(isub,jsub))
                enddo

                ! calculate the density
                call density_on_grid(sub_response_dens, mdl%fine_grid,mdl%dbl_grid, &
                     mdl%cell, mdl%fftbox, eff_kernel, &
                     cond_rep%ngwf_overlap%m(isub,jsub), cond_rep%ngwfs_on_grid(isub),&
                     cond_ngwf_basis(isub), cond_rep%ngwfs_on_grid(jsub), &
                     cond_ngwf_basis(jsub))

                response_dens=response_dens+sub_response_dens

                do is=1, pub_num_spins
                   call sparse_destroy(eff_kernel(is))
                enddo

             end do
          end do

          ! set filename
          write(filename,*) icount
          write(tempfile,*) '_TDDFT_elec_dens_'
          write(filename,*) trim(adjustl(tempfile))//trim(adjustl(filename))

          if(pub_num_spins==2) then
             response_dens(:,:,:,1)=response_dens(:,:,:,1)+&
               response_dens(:,:,:,2)
          endif

          call visual_scalarfield(response_dens(:,:,:,1),mdl%fine_grid, &
               mdl%cell,'e-h response density (in e/ang^3) for:', filename, &
               mdl%elements, ANGSTROM**3)
       enddo

       deallocate(sub_response_dens,stat=ierr)
       call utils_dealloc_check('lr_tddft_properties', 'sub_response_dens', ierr)

    endif

    ! Deallocate temporary data
    call sparse_embed_destroy(temp_vec)
    do is=1, pub_num_spins
       call sparse_embed_destroy(val_kernel(is))
       call sparse_embed_destroy(cond_kernel(is))
    enddo
    deallocate(eff_kernel, stat=ierr)
    call utils_dealloc_check('lr_tddft_properties', 'eff_kernel',ierr)
    deallocate(val_kernel, stat=ierr)
    call utils_dealloc_check('lr_tddft_properties', 'val_kernel',ierr)
    deallocate(cond_kernel, stat=ierr)
    call utils_dealloc_check('lr_tddft_properties','cond_kernel',ierr)
    call sparse_embed_destroy(Ktrans)
    call sparse_embed_destroy(KSv)
    call sparse_embed_destroy(ScK)
    deallocate(oscillator, stat=ierr)
    call utils_dealloc_check('lr_tddft_properties', 'oscillator', ierr)
    deallocate(lifetime, stat=ierr)
    call utils_dealloc_check('lr_tddft_properties', 'lifetime', ierr)

  end subroutine lr_tddft_properties

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_spectrum(energies, oscillator, lifetimes, num_states)

    !========================================================================!
    ! Subroutine plots spectrum of excited states, using the inverse lifetime!
    ! of each state as a realistic broadening factor.                        !
    !========================================================================!

    use comms, only: comms_barrier, pub_on_root
    use constants, only: HARTREE_IN_EVS, PI
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit,&
      utils_open_unit_check, utils_close_unit_check
    use rundat, only: pub_rootname, pub_lr_tddft_spectrum_smear

    implicit none

    ! Arguments
    integer, intent(in) :: num_states
    real(kind=DP), intent(in) :: energies(num_states)
    real(kind=DP), intent(in) :: oscillator(num_states)
    real(kind=DP), intent(in) :: lifetimes(num_states)

    ! Local Variables
    real(kind=DP), allocatable, dimension(:) :: sample_points
    real(kind=DP), allocatable, dimension(:) :: spectrum_points
    integer :: num_points, icount, jcount, ierr
    real(kind=DP) :: start_omega, end_omega, step_size, temp_sum
    real(kind=DP) :: current_step, width, gauss_norm
    character(50) :: filename, suffix
    integer :: output_unit


    ! define output variables
    write(suffix,*) '.tddft_spectrum'
    write(filename,'(a)') trim(adjustl(pub_rootname))//trim(adjustl(suffix))

    ! define variables
    num_points=1000
    allocate(sample_points(num_points), stat=ierr)
    call utils_alloc_check('plot_spectrum', 'sample_points', ierr)
    allocate(spectrum_points(num_points), stat=ierr)
    call utils_alloc_check('plot_stectrum', 'spectrum_points', ierr)
    start_omega=energies(1)-0.02_DP
    end_omega=energies(num_states)+0.02_DP
    step_size=(end_omega-start_omega)/(1.0_DP*num_points)
    current_step=start_omega
    width = pub_lr_tddft_spectrum_smear
    gauss_norm = 1.0_DP / (sqrt(2*PI)*width)

    do icount=1, num_points
       sample_points(icount)=current_step
       temp_sum=0.0_DP
       do jcount=1, num_states
          ! check for infinite lifetimes
          !if (lifetimes(jcount)<1.0e11_DP) then
          !   width=1.0_DP/(lifetimes(jcount))
          !else
          !   width=1.0e-11_DP
          !endif
          ! calculate lorentzian for each excitation energy. sum them up
          !temp_sum=temp_sum+width/((current_step-energies(jcount))**2.0_DP&
          !     +width**2.0_DP)
          ! calculate Gaussian for each excitation energy. sum them up
          temp_sum = temp_sum + oscillator(jcount) * gauss_norm * &
               exp(-0.5_DP*(current_step-energies(jcount))**2/width**2)
       enddo
       spectrum_points(icount)=temp_sum
       current_step=current_step+step_size
    enddo

    ! write into file
    if(pub_on_root) then
       output_unit=utils_unit()
       open(output_unit,file=filename, iostat=ierr)
       call utils_open_unit_check('lr_tddft_spectrum',filename,ierr)
       do icount=1, num_points
          if (pub_on_root) then
             write(output_unit,*) sample_points(icount)*HARTREE_IN_EVS, &
                  spectrum_points(icount)
          endif
       enddo
       close(output_unit, iostat=ierr)
       call utils_close_unit_check('lr_tddft_spectrum',filename,ierr)
    endif

    ! Deallocate data
    deallocate(sample_points, stat=ierr)
    call utils_dealloc_check('plot_spectrum', 'sample_points', ierr)
    deallocate(spectrum_points, stat=ierr)
    call utils_dealloc_check('plot_stectrum', 'spectrum_points', ierr)


  end subroutine lr_tddft_spectrum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module lr_tddft

