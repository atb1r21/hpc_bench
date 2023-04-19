!================================================================!
!                                                                !
!                      NGWF gradients module                     !
!                                                                !
! This module calculates gradients of the NGWFs with respect to  !
! their expansion coefficients in the psinc basis.               !
!----------------------------------------------------------------!
! Originally written by Chris-Kriton Skylaris and                !
! Arash A. Mostofi in 2000 and 2001.                             !
! Subsequent modifications by Chris-Kriton Skylaris,             !
! Arash A. Mostofi, Peter D. Haynes and Nicholas D. M. Hine.     !
! Extended to polarisable embedding by Jacek Dziedzic, June 2017.!
!================================================================!

!-----------------------------------------------------------------------!
! jd: For reference see                                                 !
!     Mostofi, Arash, On Linear-Scaling Methods for Quantum-Mechanical  !
!     First Principles Calculation, PhD thesis, University of           !
!     Cambridge, 2004, pp. 136-137.                                     !
!     and                                                               !
!     Ruiz Serrano, Alvaro, Computational methods for density           !
!     functional theory calculations on insulators and metals based     !
!     on localised orbitals, PhD thesis, University of Soutahmpton,     !
!     2013, pp. 93-100.                                                 !
!-----------------------------------------------------------------------!


module ngwf_gradient

  use constants, only: DP

  implicit none

  private

  real(kind=DP), allocatable, dimension(:,:,:) :: precond_func_real
  real(kind=DP), allocatable, dimension(:,:,:) :: precond_func_recip
  real(kind=DP), allocatable, dimension(:,:,:) :: smooth_func_real
  real(kind=DP), allocatable, dimension(:,:,:) :: smooth_func_recip

  public :: ngwf_gradient_exit
  public :: ngwf_gradient_lnv
  public :: ngwf_gradient_paw_precond_init
  public :: ngwf_gradient_paw_precond_exit

  ! jd: Number of matrices in gradients (parameters shared by ngwf_gradient_lnv,
  !     ngwf_gradient_batch)
  integer, parameter :: nmat=4        ! # entries in coeff_mat


contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_gradient_lnv(contra_grad, cov_grad,      &          ! in/out
       denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, &       ! input
       nl_projectors, lhxc_fine, ham, hub, mu, mdl, hfxstate, ireg, & ! input (ham is in/out)
       muext, val_dkn, val_rep, val_ham, val_ngwf_basis, &            ! opt input
       cond_shift, step, kpt, dfdtau_fine )                           ! opt input

    !==========================================================================!
    ! This subroutine calculates the gradient of the LNV total energy          !
    ! function with respect to the expansion coefficients of the NGWFs.        !
    !==========================================================================!
    ! Arguments:                                                               !
    ! contra_grad (in/out) : Contravariant NGWF gradient in ppd representation !
    !    for the NGWFs of pub_my_proc_id. Initial contents will be irrelevant  !
    !    as they are zeroed on entry.                                          !
    ! cov_grad (in/out)    : Covariant NGWF gradient in ppd representation     !
    !    for the NGWFs of pub_my_proc_id. Initial contents will be irrelevant  !
    !    as they are zeroed on entry.
    ! denskern (input) : Density kernel in DKERN format.                       !
    ! rep (input) : NGWF Representation (functions and matrices).              !
    ! ngwf_basis (input) : Function basis describing the NGWFs                 !
    ! proj_basis (input) : Function basis describing the nonlocal projectors   !
    ! hub_proj_basis (input) : Function basis describing Hubbard projectors    !
    ! nl_projectors (input) : Projector set for nonlocal pseudo projectors     !
    ! lhxc_fine (input) : Total local potential in the fine grid of the        !
    !    whole simulation cell.                                                !
    ! ham (in/out) : Hamiltonian matrix in NGWF_HAM format.                    !
    !                (%hfexchange member is modified).                         !
    ! dijhat (input) : Screened nonlocal energies                              !
    ! hub (input) : HUBBARD_MODEL type defining a Hubbard model                !
    ! mu (input) : Lagrange multiplier (original version) or fudge             !
    !    parameter (Millam-Scuseria version) of LNV function.                  !
    ! mdl (input)     : Container for cell, grid, pseudo_sp, paw_sp.           !
    ! ireg (input)    : Active region index.                                   !
    ! val_dkn (input) : Valence density kernel (optional)                      !
    ! val_rep (input) : Valence NGWF representation (optional)                 !
    ! val_ham (input) : Valence Hamiltonian matrix (optional)                  !
    ! val_ngwf_basis (input) : Function basis describing the valence NGWF      !
    !    basis (optional)                                                      !
    ! cond_shift (input) : Value by which the projected conduction Hamiltonian !
    !    is currently being shifted (optional)                                 !
    ! dfdtau_fine     (inout) : Gradient of XC energy per unit volume wrt      !
    !                           KE density energy (optional)                   !
    ! muext (input) : Chemical potential calculated for the NGWF cycle in the  !
    !                 mermin mode                                              !
    ! step (input) : trial_length of the step used to calculate Mermin gradient!
    ! with non-linear muext contributions                                      !
    !--------------------------------------------------------------------------!
    ! Key internal variables:                                                  !
    !   batch_size: This is set equal to pub_fftbox_batch_size which comes from    !
    !     the rundat module. The value of batch size determines how large is   !
    !     the batch of accumulated fftboxes. Increasing this number decreases  !
    !     the communication per processor but increases the allocated memory   !
    !     per processor.                                                       !
    !==========================================================================!
    ! Originally written by Chris-Kriton Skylaris in January 2001,             !
    ! to use a "pair-box".                                                     !
    ! Rewritten by Arash Mostofi in 2003 so that it uses a "triple-box".       !
    ! Rewritten by Chris-Kriton Skylaris on 23/11/2003 so that it runs on      !
    ! parallel computers.                                                      !
    ! Modified by Chris-Kriton Skylaris on 14/02/2005 to implement mixing      !
    ! of occupancy preconditioning.                                            !
    ! Modified by Peter Haynes to use parallel SPAM 2, July 2006               !
    ! DFT+U added by David O'Regan, April 2009                                 !
    ! Adapted to use SPAM3 matrices by Nicholas Hine, July 2009                !
    ! Modified by Laura Ratcliff for conduction calculations, Oct 2010         !
    ! Modified by Alvaro Ruiz Serrano for kernel DIIS, November 2010.          !
    ! Modified by Jacek Dziedzic for polarisable embedding, June 2017.         !
    ! Misleading intents in header comment fixed by Jacek Dziedzic, June 2017. !
    ! Modified for embedding by Robert Charlton, July 2017.                    !
    ! Modified for embedding and conduction calculations by Joseph Prentice,   !
    ! June 2018                                                                !
    !==========================================================================!

    use comms, only: comms_reduce, pub_on_root
    use constants, only: DP, stdout, max_spins, REP_SWEX_HFX_OTHER, VERBOSE, &
         CRLF
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero, &
         data_fftbox_alloc, data_fftbox_dealloc, data_functions_scale, &
        data_functions_copy, data_functions_alloc, data_functions_axpy, &
        data_functions_dealloc
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_batch_col_start, function_ops_brappd_ketppd, &
        function_ops_sum_ppd_funcs
    use geometry, only: POINT ! agrecokpt
    use hf_exchange, only: HFX_STATE, hf_exchange_calculate
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: DKERN
    use model_type, only: MODEL, REGION
    use polarisable_embedding, only: polarisable_embedding_ngwf_gradient
    use potential, only: potential_input_to_workspace
    use projectors, only: PROJECTOR_SET
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use restart, only: restart_ngwfs_tightbox_output
    use rundat, only: pub_fftbox_batch_size, pub_use_hfx, pub_hubbard, &
         pub_any_nl_proj, pub_aug, pub_cond_calculate, &
         pub_precond_real, pub_precond_recip, pub_debug_on_root, pub_num_spins, &
         pub_num_kpoints, PUB_1K, pub_eda_scfmi, pub_xc_ke_density_required, &
         pub_pol_emb_qmstar, pub_pol_emb_write_vacuum_restart, pub_rootname, &
         pub_output_detail, pub_use_activehfx, pub_active_region, &
         pub_emft, pub_do_fandt, pub_scissor_ngroups, pub_ngwf_gradient_needed,&
         pub_inner_loop_iteration, pub_mermin, pub_mermin_smearing_width, &
         pub_spin_fac
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_transpose, &
         sparse_transpose_structure, sparse_copy, sparse_scale
    use sparse_array, only: SPAM3_ARRAY, sparse_array_create, &
         sparse_array_destroy
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_copy, sparse_embed_transpose, sparse_embed_array_create, &
         sparse_embed_array_destroy, sparse_embed_write, &
         SPAM3_EMBED_ARRAY, sparse_embed_diagnose, sparse_embed_scale, &
         sparse_embed_array_create_sub, sparse_embed_trace, &
         sparse_embed_create_sub, sparse_embed_array_to_sparse_array, &
         sparse_embed_extract_from_array, sparse_embed_axpy, &
         sparse_embed_product, sparse_embed_destroy_extracted_array
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_postfix_to_ngwf_set_name, utils_abort
    use xc, only: xc_embed_swap_functional

    implicit none

    ! Arguments
    type(DKERN), intent(inout)         :: denskern
    type(NGWF_REP), intent(in)         :: rep
    type(FUNC_BASIS), intent(in)       :: ngwf_basis(rep%nsub)
    type(FUNC_BASIS), intent(in)       :: proj_basis(rep%nsub)
    type(FUNC_BASIS), intent(in)       :: hub_proj_basis(rep%nsub)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(rep%nsub)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(MODEL), intent(in)            :: mdl
    type(HFX_STATE), intent(inout), target:: hfxstate
    integer, intent(in)                :: ireg
    real(kind=DP), intent(in)          :: lhxc_fine(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, rep%nsub)
    type(FUNCTIONS), intent(inout)     :: contra_grad(rep%nsub)
    type(FUNCTIONS), intent(inout)     :: cov_grad(rep%nsub)
    type(NGWF_HAM), intent(inout)      :: ham ! jd: %hfexchange is modified
    real(kind=DP), intent(in)          :: mu(max_spins)
    ! lr408: Optional conduction arguments
    type(SPAM3_EMBED_ARRAY), optional, intent(in) :: val_dkn
    type(SPAM3_EMBED), optional, intent(in)       :: val_ham(pub_num_spins)
    type(FUNC_BASIS), optional, intent(in)  :: val_ngwf_basis(rep%nsub)
    type(NGWF_REP), optional, intent(in)    :: val_rep
    real(kind=dp), optional, intent(in)     :: cond_shift
    ! JCW: Optional gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density.
    real(kind=DP), optional, intent(in)     :: dfdtau_fine(:,:,:,:)
    ! JCW: In general, we have
    ! JCW:    dfdtau_fine(mdl%fine_grid%ld1,&
    ! JCW:    mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    ! JCW: or
    ! JCW:    dfdtau_fine(1,1,1,1)
    ! JCW: when it is unneeded.

    ! agrecokpt: optional k-point argument
    type(POINT), optional, intent(in)       :: kpt

    ! Local Variables
    type(SPAM3_EMBED_ARRAY), allocatable :: coeff_mat(:)
    type(SPAM3_EMBED), allocatable :: aux_coeff_mat(:,:)
    type(SPAM3_EMBED), allocatable :: proj_coeff_mat(:,:)
    type(SPAM3_EMBED), allocatable :: hub_proj_coeff_mat(:,:) !ddor
    type(SPAM3_EMBED), allocatable :: cond_coeff_mat(:,:)
    type(SPAM3_EMBED), allocatable :: scissor_coeff_mat(:,:) !ny
    type(SPAM3_EMBED)              :: ps_overlap
    type(FFTBOX_DATA), allocatable :: fftbox_batch(:,:,:)
    real(kind=DP), allocatable :: lhxc_dbl(:,:,:,:,:)
    ! JCW: Gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density on double grid. This has target attribute, because
    ! JCW: a pointer is assigned to this array in ngwf_gradient_batch
    real(kind=DP), allocatable, target :: dfdtau_dbl(:,:,:,:)
    integer :: imat
    integer :: is, ik
    integer :: ierr
    integer :: batch_size
    integer :: batch_count
    integer :: n_batches
    integer :: local_start, local_end, local_len
    integer :: max_current_size ! maximum batch size over all procs
    integer, allocatable :: fa_box_start(:,:)
    integer, allocatable :: fa_start_in_box(:,:)
    character(len=128) :: filename
    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt
    type(POINT) :: loc_kpt
    logical :: loc_is_force_call
    character(len=*), parameter :: myself = 'ngwf_gradient_lnv'
    integer :: isub
    ! rc2013: horrible workaround for spin/embedding issue
    type(SPAM3), allocatable :: kern_array(:), tc_kern_array(:), val_kern_array(:)
    type(SPAM3_ARRAY) :: temp_kern, tc_temp_kern
    !ep: mermin gradient variables
    type(FUNCTIONS) :: sav_grad(rep%nsub)
    type(FUNCTIONS) :: sav_grad_cov(rep%nsub)
    type(FUNCTIONS) :: tmp_direction(rep%nsub)
    type(FUNCTIONS) :: tmp_direction_cov(rep%nsub)
    type(SPAM3_EMBED) :: ngwf_overlapa_dirb
    type(SPAM3_EMBED) :: ngwf_dira_overlapb
    type(SPAM3_EMBED) :: dira_dirb
    type(SPAM3_EMBED_ARRAY) :: mu_mat
    real(kind=DP), optional, intent(in):: step
    real(kind=DP) , optional, intent(inout):: muext(max_spins,pub_num_kpoints)
    real(kind=DP) :: muapp(max_spins,pub_num_kpoints)
    type(FFTBOX_DATA), allocatable :: save_fftbox_batch(:,:,:)
    type(FFTBOX_DATA), allocatable :: save_fftbox_batch_cov(:,:,:)
    integer :: jsub
    !ep: mermin gradient variables

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering '//myself

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         myself//': not ready yet for more than one k-point.')

    ! Start timer
    call timer_clock(myself,1)

    call utils_assert(pub_ngwf_gradient_needed, myself//': Logic error. &
         &Logic in energy_and_force for establishing if NGWF gradient will need&
         & to be calculated is wrong.')

    pub_inner_loop_iteration = -1 ! jd: We're not in the LNV, KDIIS or EDFT loop

    if (.not.pub_xc_ke_density_required) then
       ! JCW: Non-tau-dependent XC functional
       ! JCW: dfdtau_fine is either not present, or has a length of 1
       if ( present(dfdtau_fine ) ) then
          call utils_assert( &
               size(dfdtau_fine) == 1,&
               "Error in ngwf_gradient_lnv: &
               &pub_xc_ke_density_required is false, but dfdtau_fine is present &
               &and has a size > 1.")
       end if
    else
       ! JCW: tau-dependent XC functional
       ! JCW: Check dfdtau_fine is present
       call utils_assert( present(dfdtau_fine),&
            "Error in ngwf_gradient_lnv: &
            &pub_xc_ke_density_required is true, but dfdtau_fine is not &
            &present.")
    end if

    ! ######## INITIALISATIONS #################################
    ! check if this is routine is called from a inside force_mod
    ! while in a cond task.
    loc_is_force_call=.false.
    if(pub_cond_calculate .and. .not. present(val_rep)) then
       loc_is_force_call=.true.
    endif

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(1)%iscmplx

    ! agrecokpt
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    do isub=1,rep%nsub
       call data_set_to_zero(contra_grad(isub))
       call data_set_to_zero(cov_grad(isub))
    end do
    !ep: mermin grad var init
    if (pub_mermin) then
       do isub=1,mdl%nsub
          call data_functions_alloc(sav_grad(isub), &
               ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
          call data_functions_alloc(sav_grad_cov(isub), &
               ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
          call data_set_to_zero(sav_grad(isub))
          call data_set_to_zero(sav_grad_cov(isub))
       end do
       call sparse_embed_create(ngwf_dira_overlapb,rep%overlap)
       call sparse_embed_create(ngwf_overlapa_dirb,rep%overlap)
       call sparse_embed_create(dira_dirb,rep%overlap)
       call sparse_embed_array_create(mu_mat,denskern%kern)
    end if
    !ep: mermin grad var init

    batch_size = pub_fftbox_batch_size

    ! ndmh: projector-ngwf overlap matrix
    if (pub_any_nl_proj.or.pub_aug) then
       ! rc2013: get transpose structure directly from sparse_embed_create
       call sparse_embed_create(ps_overlap, rep%sp_overlap, trans=.true.)
       call sparse_embed_transpose(ps_overlap, rep%sp_overlap)
    else
       ! rc2013: even if we don't need the matrix, allocate the array
       allocate(ps_overlap%m(rep%nsub,rep%nsub), stat=ierr)
       call utils_alloc_check(myself, 'ps_overlap%m', ierr)
    end if

    ! ndmh: allocate storage for coefficient matrices
    allocate(coeff_mat(nmat),stat=ierr)
    call utils_alloc_check(myself,'coeff_mat',ierr)
    allocate(cond_coeff_mat(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check(myself,'cond_coeff_mat',ierr)
    allocate(proj_coeff_mat(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check(myself,'proj_coeff_mat',ierr)
    allocate(hub_proj_coeff_mat(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check(myself,'hub_proj_coeff_mat',ierr)
    allocate(scissor_coeff_mat(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check(myself,'scissor_coeff_mat',ierr)

    call sparse_embed_array_create(coeff_mat(1), denskern%kern, n_kpoints=PUB_1K)
    call sparse_embed_array_create(coeff_mat(2), coeff_mat(1), rep%overlap, &
         n_kpoints=PUB_1K)
    call sparse_embed_array_create(coeff_mat(3), denskern%kern)
    call sparse_embed_array_create(coeff_mat(4), denskern%kern, rep%overlap)

    do is=1,pub_num_spins
       ! ndmh: create elements of proj_coeff_mat array of projector
       ! ndmh: coefficient matrices
       if (pub_aug) then
          call sparse_embed_create(proj_coeff_mat(is,1), ps_overlap, &
               denskern%kern%m(is,PUB_1K))
          call sparse_embed_create(proj_coeff_mat(is,2), ps_overlap, &
               coeff_mat(2)%m(is,PUB_1K))
       end if
       if (pub_any_nl_proj.or.pub_aug) then
          call sparse_embed_create(proj_coeff_mat(is,3), ps_overlap, &
               denskern%kern%m(is,PUB_1K))
          call sparse_embed_create(proj_coeff_mat(is,4), ps_overlap, &
               coeff_mat(2)%m(is,PUB_1K))
       end if
       ! lr408: Conduction coefficient matrices
       if (pub_cond_calculate .and. .not. loc_is_force_call) then
          cond_coeff_mat(is,3)%structure = 'Mc'
          ! agrecocmplx
          call sparse_embed_create(cond_coeff_mat(is,3),iscmplx=loc_cmplx)
          if(trim(adjustl(rep%overlap%structure)) == "S") then
             call sparse_embed_create(cond_coeff_mat(is,4),rep%overlap,cond_coeff_mat(is,3))
          else
             call sparse_embed_create(cond_coeff_mat(is,4),cond_coeff_mat(is,3), &
                  rep%overlap)
          end if
       end if
       ! ddor: Hubbard projector coefficient matrices
       if (pub_hubbard) then
          call sparse_embed_create(hub_proj_coeff_mat(is,3), &
               rep%hub_overlap_t, denskern%kern%m(is,PUB_1K))
          call sparse_embed_create(hub_proj_coeff_mat(is,4), &
               rep%hub_overlap_t,coeff_mat(2)%m(is,PUB_1K))
       endif

       ! ny: Scissor coefficient matrices
       if (pub_scissor_ngroups > 0) then
          call sparse_embed_create(scissor_coeff_mat(is,3), &
               rep%overlap, denskern%kern%m(is,PUB_1K))
          call sparse_embed_create(scissor_coeff_mat(is,4), &
               rep%overlap,coeff_mat(2)%m(is,PUB_1K))
       endif
    end do ! spins

    ! ndmh: Initialise the NGWF coefficients:
    ! ndmh: coeff_mat(:,1) (qmat) = (3LHL -2LSLHL -2LHLSL)
    ! ndmh: coeff_mat(:,2) (tc_qmat) = (3LHL -2LSLHL -2LHLSL)S
    ! ndmh: coeff_mat(:,3) (pur_denskern) = \tilde{K}
    ! ndmh: coeff_mat(:,4) (tc_pur_denskern) = \tilde{K}S

    ! ars: If kernel DIIS the NGWF coefficients are:
    ! ars: coeff_mat(:,1) (qmat) = Q = -[Ne/tr(KS)]^2 *KHK - 2(mu-[Ne/tr(KS)]^2 * nu] *K
    ! ars: coeff_mat(:,2) (tc_qmat) = QS
    ! ars: coeff_mat(:,3) (pur_denskern) = P = [2Ne/tr(KS)] *K - [Ne/tr(KS)]^2 *KSK
    ! ars: coeff_mat(:,4) (tc_pur_denskern) = PS
    ! ars: mu = tr[KH]/tr[KS]
    ! ars: nu = tr[HKSK]/tr[KS]

    ! ars: In EDFT the NGWF coefficients are:
    ! ars: coeff_mat(:,1) = (qmat) = Q = -KHS^.
    ! ars: coeff_mat(:,2) = (tc_qmat) = QS = -KH.
    ! ars: coeff_mat(:,3) = (pur_denskern) = K.
    ! ars: coeff_mat(:,4) = (tc_pur_denskern) = KS.

    ! ndmh: Initialise (optionally) the nonlocal projector coefficients:
    ! ndmh: proj_coeff_mat(:,1) rQ = O<proj|ngwf>*qmat (PAW only)
    ! ndmh: proj_coeff_mat(:,2) tc_rQ  = O<proj|ngwf>*tc_qmat (PAW only)
    ! ndmh: proj_coeff_mat(:,3) rk = D<proj|ngwf>*pur_denskern (NCPP,PAW)
    ! ndmh: proj_coeff_mat(:,4) tc_rk = D<proj|ngwf>*tc_pur_denskern (NCPP,PAW)
    ! ndmh: where O_ij is the block-diagonal projector overlap matrix and
    ! ndmh: D_ij is the matrix of PAW nonlocal energies (or KB coeffs for NCPPs)

    ! ndmh: Initialise (optionally) the Hubbard projector coefficients:
    ! ndmh: hub_proj_coeff_mat(:,3) hub_wk = G<proj|ngwf>*pur_denskern
    ! ndmh: hub_proj_coeff_mat(:,4) tc_hub_wk = G<proj|ngwf>*tc_pur_denskern
    ! ndmh: where G is the Hubbard Hamiltonian (block-diagonal)

    ! ndmh: Initialise (optionally) the valence (jd: sic!) NGWF coefficients
    ! ndmh: in a conduction NGWF optimisation,
    ! ndmh: cond_coeff_mat(:,3) cond_grad, cond_coeff_mat(:,4) tc_cond_grad
    ! ndmh: cond_grad = M<val_ngwf|cond_ngwf>*pur_denskern(cond)
    ! ndmh: tc_cond_grad = M<val_ngwf|cond_ngwf>*tc_pur_denskern(cond)
    ! ndmh: where M = w.Kval.Sval.Kval - Kval.Hval.Kval and w is the
    ! ndmh: conduction state shift parameter in the projection, Kval is
    ! ndmh: the valence density kernel, Sval the valence overlap matrix
    ! ndmh: and Hval the valence Hamiltonian.

    if (.not. pub_eda_scfmi) then
       call ngwf_gradient_coeffs(coeff_mat(1)%m(:,PUB_1K), &
            coeff_mat(2)%m(:,PUB_1K),coeff_mat(3), &
            coeff_mat(4),proj_coeff_mat(:,1),proj_coeff_mat(:,2), &
            proj_coeff_mat(:,3),proj_coeff_mat(:,4), &
            hub_proj_coeff_mat(:,3),hub_proj_coeff_mat(:,4), &
            scissor_coeff_mat(:,3),scissor_coeff_mat(:,4), &
            cond_coeff_mat(:,3),cond_coeff_mat(:,4), &
            denskern, &
            rep%overlap,ham%ham,rep%inv_overlap,ps_overlap,rep%sp_overlap,&
            hub,rep%hub_overlap_t, &
            mu,rep%n_occ,ham%dijhat,mdl,loc_is_force_call,&
            val_dkn, val_rep%overlap, val_rep%sp_overlap, val_ham, &
            rep%cross_overlap, cond_shift, rep, ngwf_basis)
    else
       ! mjsp: SCF-MI: use the supermolecule (full) overlap for
       ! mjsp: calculation of the NGWF gradient
       call ngwf_gradient_coeffs(coeff_mat(1)%m(:,PUB_1K), &
            coeff_mat(2)%m(:,PUB_1K),coeff_mat(3), &
            coeff_mat(4),proj_coeff_mat(:,1),proj_coeff_mat(:,2), &
            proj_coeff_mat(:,3),proj_coeff_mat(:,4), &
            hub_proj_coeff_mat(:,3),hub_proj_coeff_mat(:,4), &
            scissor_coeff_mat(:,3),scissor_coeff_mat(:,4), &
            cond_coeff_mat(:,3),cond_coeff_mat(:,4), &
            denskern, &
            rep%overlap_scfmi_full,ham%ham,rep%inv_overlap_scfmi_full, &
            ps_overlap,rep%sp_overlap,&
            hub,rep%hub_overlap_t, &
            mu,rep%n_occ,ham%dijhat,mdl,loc_is_force_call,&
            val_dkn, val_rep%overlap, val_rep%sp_overlap, val_ham, &
            rep%cross_overlap, cond_shift, rep, ngwf_basis)
    end if

    ! ### END INITIALISATIONS #################################

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Completed initialisation in ngwf_gradient_lnv'

    ! ndmh: allocate preconditioner
    if (pub_precond_recip.or.pub_precond_real) then
       call ngwf_grad_init_precond_recip(mdl%fftbox)
    end if
    if (pub_precond_real) then
       call ngwf_grad_init_precond_real(mdl%fftbox,mdl%cell)
    end if

    ! ndmh: in case the fine grid is denser than the double grid, create a
    ! ndmh: temporary double grid
    allocate(lhxc_dbl(mdl%dbl_grid%ld1,mdl%dbl_grid%ld2,&
         mdl%dbl_grid%max_group_slabs12,pub_num_spins, rep%nsub), stat=ierr)
    call utils_alloc_check(myself,'lhxc_dbl',ierr)
    if (pub_xc_ke_density_required) then
       ! JCW: do the same for the dfdtau_potential
       allocate(dfdtau_dbl(mdl%dbl_grid%ld1,mdl%dbl_grid%ld2,&
            mdl%dbl_grid%max_group_slabs12,pub_num_spins), stat=ierr)
       call utils_alloc_check(myself,'dfdtau_dbl',ierr)
    else
       ! JCW: allocate a single element dummy array
       allocate(dfdtau_dbl(1,1,1,1), stat=ierr)
       call utils_alloc_check(myself,'dfdtau_dbl',ierr)
    end if

    ! ndmh: fb start positions in box and start positions of FFT boxes
    allocate(fa_box_start(3, batch_size), stat=ierr)
    call utils_alloc_check(myself,'fa_box_start',ierr)
    allocate(fa_start_in_box(3, batch_size), stat=ierr)
    call utils_alloc_check(myself,'fa_start_in_box',ierr)

    ! ndmh: allocate storage for fftboxes for this batch
    allocate(fftbox_batch(pub_num_spins, nmat, batch_size), stat=ierr)
    call utils_alloc_check(myself,'fftbox_batch',ierr)
    do batch_count = 1, batch_size
       do imat = 1, nmat
          do is = 1, pub_num_spins
             call data_fftbox_alloc(fftbox_batch(is, imat, batch_count), &
                  mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
                  mdl%fftbox%total_pt3, loc_cmplx)
          end do
       end do
    end do

    !ep: alloc fftboxes for mermin method
    if (pub_mermin) then
       allocate(save_fftbox_batch(pub_num_spins, 1, batch_size), stat=ierr)
       call utils_alloc_check('ngwf_gradient_batch','save_fftbox_batch',ierr)
       allocate(save_fftbox_batch_cov(pub_num_spins, 1, batch_size), stat=ierr)
       call utils_alloc_check('ngwf_gradient_batch','save_fftbox_batch_cov', &
            ierr)
       do batch_count = 1, batch_size
          do is = 1, pub_num_spins
            call data_fftbox_alloc(save_fftbox_batch(is, 1, batch_count), &
                 mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
                 mdl%fftbox%total_pt3, loc_cmplx)
            call data_fftbox_alloc(save_fftbox_batch_cov(is, 1, batch_count), &
                 mdl%fftbox%total_ld1, mdl%fftbox%total_ld2, &
                 mdl%fftbox%total_pt3, loc_cmplx)
          end do
       end do
    end if
   !ep

    ! ndmh: filter (or just copy) the lhxc potential to the double grid
    do is=1,pub_num_spins
       ! rc2013: transfer all LHXC potentials to the double grid
       do isub=1,rep%nsub
          call potential_input_to_workspace(lhxc_dbl(:,:,:,is,isub), &
               lhxc_fine(:,:,:,is,isub),mdl%dbl_grid,mdl%fine_grid)
       end do
       if (pub_xc_ke_density_required) then
          ! JCW: do the same for the dfdtau_potential
          call potential_input_to_workspace(dfdtau_dbl(:,:,:,is), &
               dfdtau_fine(:,:,:,is),mdl%dbl_grid,mdl%fine_grid)
       end if
    end do

    ! rc2013: need to come up with a better way of handling this...
    allocate(kern_array(pub_num_spins),stat=ierr)
    call utils_alloc_check('ngwf_gradient_lnv','kern_array',ierr)
    allocate(tc_kern_array(pub_num_spins),stat=ierr)
    call utils_alloc_check('ngwf_gradient_lnv','tc_kern_array',ierr)
    if (trim(rep%postfix) == 'c') then
       allocate(val_kern_array(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_gradient_lnv','val_kern_array',ierr)
    end if

    ! jme: KPOINTS_DANGER: only one k-point considered here!
    ! Copy the first k-point info of coeff_mat into aux_coeff_mat
    allocate(aux_coeff_mat(pub_num_spins,nmat), stat=ierr)
    call utils_alloc_check(myself, 'aux_coeff_mat', ierr)
    do is=1,pub_num_spins
       ! ndmh: create elements of aux_coeff_mat array of NGWF coefficient matrices
       do imat = 1, nmat
          call sparse_embed_create(aux_coeff_mat(is,imat), &
               coeff_mat(imat)%m(is,PUB_1K))
          call sparse_embed_copy(aux_coeff_mat(is,imat), &
               coeff_mat(imat)%m(is,PUB_1K))
       end do
    end do
    ! jme: KPOINTS_DANGER end

    ! rc2013: loop over regions (if necessary)
    do isub=1,rep%nsub
       ! rc2013: if using freeze-and-thaw optimisation only build the batch for
       ! the active subsystem
       if(pub_do_fandt .and. isub.ne.ireg) cycle

       ! cks: number of row-steps per row-block
       n_batches = ngwf_basis(isub)%max_on_proc / batch_size
       if (mod(ngwf_basis(isub)%max_on_proc,batch_size) > 0) n_batches = n_batches + 1

       ! ndmh: loop over batches of NGWFs
       local_start = 1
       do batch_count=1,n_batches

          if (pub_debug_on_root) write(stdout,'(a,2(i4,a))') 'DEBUG: Batch loop ',&
               batch_count, ' of ',n_batches,' in '//myself

          local_end = min(local_start+batch_size-1,ngwf_basis(isub)%proc_num)
          local_len = local_end - local_start + 1

          ! cks: maximum size of current batch over all procs
          max_current_size = local_len
          call comms_reduce('MAX', max_current_size)

          call function_ops_batch_col_start(fa_box_start,fa_start_in_box, &
               batch_size,local_start,local_end,mdl%fftbox,mdl%cell, &
               ngwf_basis(isub))

          ! cks: zero batch (all spins, all matrices) before accumulation
          call data_set_to_zero(fftbox_batch(:,:,1:local_len))
          ! ep: zero batch (all spins, all matrices) before accumulation
          ! ep: mermin vars
          if (pub_mermin) then
              call data_set_to_zero(save_fftbox_batch(:,:,1:local_len))
              call data_set_to_zero(save_fftbox_batch_cov(:,:,1:local_len))
          end if

          ! ndmh: deposit sums of various functions (the NGWFs, the Hamiltonian
          ! ndmh: acting on the NGWFS, nonlocal projectors, Hubbard projectors,
          ! ndmh: valence NGWFs in a conduction optimisation) to the accumulating
          ! ndmh: FFTboxes, then precondition the gradient, extract it from the
          ! ndmh: FFTBoxes to PPD storage, and then shave them according to the
          ! ndmh: NGWFs radii.
#ifdef GPU_PGI_NGWFGRAD
          call utils_assert(.not.pub_xc_ke_density_required, myself//': Combinatio&
               &n of GPU and tau-dependent XC functionals not implemented/tested')
          ! jcap: GPUs and embedding structures have not been
          ! implemented together yet - this won't be compiled unless
          ! GPUs are asked for, but it serves as a reminder to fix this
          call ngwf_gpu_gradient_batch(contra_grad(isub), cov_grad(isub), & ! output
               fftbox_batch, lhxc_dbl, mdl%dbl_grid, rep, ngwf_basis, proj_basis, &
               hub_proj_basis, nl_projectors, hub,  fa_box_start, fa_start_in_box, &
               batch_size, local_start, local_end, aux_coeff_mat, proj_coeff_mat, &
               hub_proj_coeff_mat, cond_coeff_mat, ps_overlap, mdl%fftbox, mdl%cell, &
               max_current_size, loc_is_force_call,&
               val_ngwf_basis, val_rep)
#else
          call ngwf_gradient_batch(contra_grad(isub), cov_grad(isub), &  ! output
               fftbox_batch, lhxc_dbl, mdl%dbl_grid, rep, ngwf_basis, proj_basis, &
               hub_proj_basis, nl_projectors, hub, fa_box_start, fa_start_in_box, &
               batch_size, local_start, local_end, aux_coeff_mat, proj_coeff_mat, &
               hub_proj_coeff_mat, scissor_coeff_mat, cond_coeff_mat, ps_overlap, &
               mdl%fftbox, mdl%cell, max_current_size, loc_is_force_call,         &
               val_ngwf_basis, val_rep, save_fftbox_batch=save_fftbox_batch,      &
               save_fftbox_batch_cov=save_fftbox_batch_cov,                       &
               sav_grad=sav_grad(isub), sav_grad_cov=sav_grad_cov(isub),          &
               kpt=loc_kpt, dfdtau_dbl = dfdtau_dbl, ireg=isub)
#endif

          local_start = local_start + batch_size
       end do
    end do

    !ep: calculate muext & muext contribution to the grad
    !THERE IS SURELY A BETTER/MORE EFFICIENT WAY TO DO IT!
    !alternative to be explored once the method is finalised !
    if ( pub_mermin ) then
       do isub=1,mdl%nsub
           call data_functions_alloc(tmp_direction(isub), &
                ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
           call data_set_to_zero(tmp_direction(isub))
       enddo
       call sparse_embed_scale(ngwf_dira_overlapb, 0.0_DP)
       call sparse_embed_scale(ngwf_overlapa_dirb, 0.0_DP)
       call sparse_embed_scale(dira_dirb, 0.0_DP)
       !tmp direction without mu
       do isub=1,mdl%nsub
           call data_functions_copy(tmp_direction(isub),cov_grad(isub))
           call data_functions_scale(tmp_direction(isub),-1.0_DP)
        enddo
        do isub=1,mdl%nsub
           do jsub=1,mdl%nsub
              !calculate <psi_a|L_b>
              call function_ops_brappd_ketppd(ngwf_overlapa_dirb%m(isub,jsub), &
                   rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                   tmp_direction(jsub), ngwf_basis(jsub), mdl%cell)
              !calculate <psi_b|L_a>
              call function_ops_brappd_ketppd(ngwf_dira_overlapb%m(isub,jsub), &
                   tmp_direction(isub), ngwf_basis(isub), &
                   rep%ngwfs_on_grid(jsub), ngwf_basis(jsub), mdl%cell)
              if (present(step)) then
                 !calculate <L_a|L_b>
                 call function_ops_brappd_ketppd(dira_dirb%m(isub,jsub), &
                      tmp_direction(isub), ngwf_basis(isub), &
                      tmp_direction(jsub), ngwf_basis(jsub), mdl%cell)
              end if
           end do
       end do
       !calculate <psi_a|L_b> + <psi_b|L_a>
       call sparse_embed_axpy(ngwf_overlapa_dirb, ngwf_dira_overlapb, 1.0_DP)
       if (present(step)) then
          !calculate <psi_a|L_b> + <psi_b|L_a> + step*<L_a|L_b>
          call sparse_embed_axpy(ngwf_overlapa_dirb, dira_dirb, step)
       end if
       do ik = 1, denskern%kern%num_kpoints
          do is=1,pub_num_spins
             !calculate K^ab*(<psi_a|L_b> + <psi_b|L_a>)
             call sparse_embed_product(mu_mat%m(is,ik),denskern%kern%m(is,ik), &
             ngwf_overlapa_dirb)
             !trace K^ab*(<psi_a|L_b> + <psi_b|L_a>)
             call sparse_embed_trace(muapp(is,ik),mu_mat%m(is,ik))
             !calculate -(1/2)*N_e*K^ab*(<psi_a|L_b> + <psi_b|L_a>) i.e. muext
             muext(is,ik)=(-1.0_DP/(2.0_DP*(rep%n_occ(is,ik))))*muapp(is,ik)
          end do
       end do

       ! rc2013: loop over regions (if necessary)
       do isub=1,rep%nsub
          ! rc2013: if using freeze-and-thaw optimisation only build the batch
          ! for
          ! the active subsystem
          if(pub_do_fandt .and. isub.ne.ireg) cycle

         ! cks: number of row-steps per row-block
          n_batches = ngwf_basis(isub)%max_on_proc / batch_size
          if (mod(ngwf_basis(isub)%max_on_proc,batch_size) > 0) n_batches = &
          n_batches + 1

         ! ndmh: loop over batches of NGWFs
          local_start = 1
          do batch_count=1,n_batches
             if (pub_debug_on_root) write(stdout,'(a,2(i4,a))') 'DEBUG: Batch &
                & loop ', batch_count, ' of ',n_batches,' in '//myself
             local_end = min(local_start+batch_size-1,ngwf_basis(isub)%proc_num)
             local_len = local_end - local_start + 1

            ! cks: maximum size of current batch over all procs
             max_current_size = local_len
             call comms_reduce('MAX', max_current_size)

             call function_ops_batch_col_start(fa_box_start,fa_start_in_box, &
                  batch_size,local_start,local_end,mdl%fftbox,mdl%cell, &
                  ngwf_basis(isub))

            ! cks: zero batch (all spins, all matrices) before accumulation
             call data_set_to_zero(fftbox_batch(:,:,1:local_len))

             call ngwf_gradient_batch(contra_grad(isub), cov_grad(isub), &  !output
                  fftbox_batch, lhxc_dbl, mdl%dbl_grid, rep, ngwf_basis, &
                  proj_basis, hub_proj_basis, nl_projectors, hub, fa_box_start, &
                  fa_start_in_box, batch_size, local_start, local_end, &
                  aux_coeff_mat, proj_coeff_mat, hub_proj_coeff_mat, &
                  scissor_coeff_mat, cond_coeff_mat, ps_overlap, &
                  mdl%fftbox, mdl%cell, max_current_size, loc_is_force_call, &
                  val_ngwf_basis, val_rep, save_fftbox_batch=save_fftbox_batch, &
                  save_fftbox_batch_cov=save_fftbox_batch_cov,  &
                  sav_grad=sav_grad(isub), sav_grad_cov=sav_grad_cov(isub), &
                  kpt=loc_kpt, dfdtau_dbl = dfdtau_dbl, ireg=isub, muext=muext)

             local_start = local_start + batch_size
          end do
       end do
    end if
    !ep: calculate muext & muext contribution to the grad
    !THERE IS SURELY A BETTER/MORE EFFICIENT WAY TO DO IT!

    ! pdh: deallocate workspace
    do batch_count = batch_size, 1, -1
       do imat = nmat, 1, -1
          do is = pub_num_spins, 1, -1
             call data_fftbox_dealloc(fftbox_batch(is, imat, batch_count))
          end do
       end do
    end do
    deallocate(fftbox_batch,stat=ierr)
    call utils_dealloc_check(myself,'fftbox_batch',ierr)
    deallocate(fa_start_in_box,stat=ierr)
    call utils_dealloc_check(myself,'fa_start_in_box',ierr)
    deallocate(fa_box_start,stat=ierr)
    call utils_dealloc_check(myself,'fa_box_start',ierr)
    deallocate(dfdtau_dbl, stat=ierr)
    call utils_dealloc_check(myself,'dfdtau_dbl',ierr)
    deallocate(lhxc_dbl, stat=ierr)
    call utils_dealloc_check(myself,'lhxc_dbl',ierr)

    !ep: dellaoc Mermin vars
    if (pub_mermin) then
       call sparse_embed_destroy(ngwf_dira_overlapb)
       call sparse_embed_destroy(ngwf_overlapa_dirb)
       call sparse_embed_destroy(dira_dirb)
       call sparse_embed_array_destroy(mu_mat)
       do isub=1,mdl%nsub
          call data_functions_dealloc(sav_grad_cov(isub))
          call data_functions_dealloc(sav_grad(isub))
       end do
       do batch_count = batch_size, 1, -1
             do is = pub_num_spins, 1, -1
                call data_fftbox_dealloc(save_fftbox_batch(is, 1, batch_count))
             end do
       end do
       deallocate(save_fftbox_batch,stat=ierr)
       call utils_dealloc_check(myself,'save_fftbox_batch',ierr)
       do batch_count = batch_size, 1, -1
             do is = pub_num_spins, 1, -1
                call data_fftbox_dealloc(save_fftbox_batch_cov(is, 1, &
                     batch_count))
             end do
       end do
       deallocate(save_fftbox_batch_cov,stat=ierr)
       call utils_dealloc_check(myself,'save_fftbox_batch_cov',ierr)
    endif
    !ep: dellaoc Mermin vars

    ! qoh: calculate NGWF gradient for Hartree-Fock exchange
    ! agrecocmplx: this is not compatible yet with complex NGWFs
    ! jcap: modified for embedding
    if (pub_use_hfx.or.(pub_use_activehfx.and.pub_emft.and.(ireg==pub_active_region))) then
       if(pub_on_root .and. pub_output_detail >= VERBOSE) then
          write(stdout,'(a,a,a)') CRLF//'HFx: NGWF gradient calculation for ', &
               trim(utils_postfix_to_ngwf_set_name(rep%postfix)), ' rep.'
       end if
       call sparse_embed_extract_from_array(kern_array,coeff_mat(3)%m(:,PUB_1K),&
            ireg,ireg)
       call sparse_embed_extract_from_array(tc_kern_array,coeff_mat(4)%m(:,PUB_1K),&
            ireg,ireg)
       if(trim(rep%postfix) == 'c') then
          call sparse_embed_extract_from_array(val_kern_array,val_dkn%m(:,PUB_1K),&
               ireg,ireg)
       end if
       do is=1,pub_num_spins
          ! rc2013: zero the inactive components of the exchange matrix
          call sparse_embed_scale(ham%hfexchange(is), 0.0_DP)
       end do
       ! jcap: check if we need to swap functionals for EMFT, to get
       ! pub_hfxfraction right
       if (pub_emft.and.pub_use_activehfx.and.(ireg==pub_active_region)) &
            call xc_embed_swap_functional(.true.)

       if(trim(rep%postfix) == '') then
          ! jd: HFx gradient for valence
          call hf_exchange_calculate(hfxstate, ham%hfexchange, &        ! inout
               rep, ireg, kern_array, kern_array, &                     ! input
               ngwf_basis(ireg), mdl%fftbox, mdl%cell, &                ! input
               mdl%regions(ireg)%elements, .true., &                    ! input
               tc_kern_array, &
               cov_grad(ireg), contra_grad(ireg), &                     ! in/out
               precond_func_recip)                                      ! input
       else if(trim(rep%postfix) == 'c') then
          ! jd: HFx gradient for conduction
          !     coeff_mat is conduction
          ! rc2013: EMBED_FIX
          call hf_exchange_calculate(hfxstate, ham%hfexchange, &        ! inout
               rep, ireg, kern_array, val_kern_array, &                 ! input
               ngwf_basis(ireg), mdl%fftbox, mdl%cell, &                ! input
               mdl%regions(ireg)%elements, .true., &                    ! input
               tc_kern_array, &                                         ! input
               cov_grad(ireg), contra_grad(ireg), &                     ! in/out
               precond_func_recip, &                                    ! input
               rep2 = val_rep, ngwf_basis2 = val_ngwf_basis(ireg), &
               basis_selector = (/1,1,2,2,REP_SWEX_HFX_OTHER/), &
               energy_prefactor = 2.0/real(pub_num_spins,kind=DP), &
               grad_prefactor = 0.5_DP * 2.0_DP)
               ! qoh's alpha-beta-gamma-delta conv'n.
               ! '1' is rep and ngwf_basis (so conduction),
               ! '2' is rep2 and ngwf_basis2 (so valence)
       else
          call utils_abort(myself//': NGWF_REP postfix '''//trim(rep%postfix)//&
               ''' is currently unsupported with HFx NGWF gradient.')
       end if
       ! jcap: swap back if we need to
       if (pub_emft.and.pub_use_activehfx.and.(ireg==pub_active_region)) &
            call xc_embed_swap_functional(.false.)

       ! rc2013: destroy the temporary arrays
       call sparse_embed_destroy_extracted_array(kern_array)
       call sparse_embed_destroy_extracted_array(tc_kern_array)
       if(trim(rep%postfix) == 'c') then
          call sparse_embed_destroy_extracted_array(val_kern_array)
          deallocate(val_kern_array,stat=ierr)
          call utils_dealloc_check('ngwf_gradient_lnv','val_kern_array',ierr)
       end if
    end if

    ! jd: Calculate NGWF gradient for polarisable embedding, if QM* rep used
    !     This accounts for terms (1) and (3).
    if(pub_pol_emb_qmstar) then
       ! jcap: convert coeff arrays to sparse_arrays - this is fine as
       ! long as we only have 1 region, since we actually want
       ! m(ireg,1)
       call sparse_embed_array_to_sparse_array(temp_kern,coeff_mat(3),1)
       call sparse_embed_array_to_sparse_array(tc_temp_kern,coeff_mat(4),1)
       call polarisable_embedding_ngwf_gradient(&
            rep, &
            temp_kern, &                                                ! input
            ngwf_basis(ireg), mdl, &                                    ! input
            tc_temp_kern, &                                             ! input
            cov_grad(ireg), contra_grad(ireg), &                        ! in/out
            precond_func_recip)                                         ! input
       ! jcap: destroy the temporary arrays
       call sparse_array_destroy(temp_kern)
       call sparse_array_destroy(tc_temp_kern)

    end if

    ! jd: For polarisable embedding in-vacuum calculation only: save
    !     kernel and NGWFs that will be needed later
    if(pub_pol_emb_write_vacuum_restart) then
       call restart_ngwfs_tightbox_output(rep%ngwfs_on_grid(ireg), ngwf_basis(ireg),&
            mdl%cell, mdl%fftbox, mdl%regions(ireg)%elements, &
            'pol_emb_vac_tightbox_'//ngwf_basis(ireg)%name, mdl%regions(ireg))
       filename = trim(pub_rootname)//'.pol_emb_kdenskern'
       call sparse_embed_write(coeff_mat(3)%m,trim(filename))
    end if

    call ngwf_gradient_exit

    ! ndmh: deallocate coefficient matrices
    do is=pub_num_spins,1,-1
       if (pub_hubbard) then !ddor
          call sparse_embed_destroy(hub_proj_coeff_mat(is,4))
          call sparse_embed_destroy(hub_proj_coeff_mat(is,3))
       endif
       if (pub_any_nl_proj.or.pub_aug) then
          call sparse_embed_destroy(proj_coeff_mat(is,4))
          call sparse_embed_destroy(proj_coeff_mat(is,3))
       end if
       if (pub_aug) then
          call sparse_embed_destroy(proj_coeff_mat(is,2))
          call sparse_embed_destroy(proj_coeff_mat(is,1))
       end if
       if (pub_cond_calculate .and. .not. loc_is_force_call) then
          call sparse_embed_destroy(cond_coeff_mat(is,4))
          call sparse_embed_destroy(cond_coeff_mat(is,3))
       end if
       do imat=4,1,-1
          call sparse_embed_destroy(aux_coeff_mat(is,imat))
       end do
       if (pub_scissor_ngroups > 0) then
          call sparse_embed_destroy(scissor_coeff_mat(is,4))
          call sparse_embed_destroy(scissor_coeff_mat(is,3))
       endif
    end do
    do imat=4,1,-1
       call sparse_embed_array_destroy(coeff_mat(imat))
    end do
    !ep: dealloc mermin var
    if (pub_mermin) then
       do isub=1,mdl%nsub
          call data_functions_dealloc(tmp_direction(isub))
       end do
    endif
    !ep

    deallocate(hub_proj_coeff_mat,stat=ierr)
    call utils_dealloc_check(myself,'hub_proj_coeff_mat',ierr)
    deallocate(proj_coeff_mat,stat=ierr)
    call utils_dealloc_check(myself,'proj_coeff_mat',ierr)
    deallocate(cond_coeff_mat,stat=ierr)
    call utils_dealloc_check(myself,'cond_coeff_mat',ierr)
    deallocate(aux_coeff_mat,stat=ierr)
    call utils_dealloc_check(myself,'aux_coeff_mat',ierr)
    deallocate(coeff_mat,stat=ierr)
    call utils_dealloc_check(myself,'coeff_mat',ierr)
    deallocate(scissor_coeff_mat,stat=ierr)
    call utils_dealloc_check(myself,'scissor_coeff_mat',ierr)
    if (pub_any_nl_proj.or.pub_aug) then
       call sparse_embed_destroy(ps_overlap)
    else
       deallocate(ps_overlap%m,stat=ierr)
       call utils_dealloc_check(myself,'ps_overlap%m',ierr)
    end if
    deallocate(kern_array,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_lnv','kern_array',ierr)
    deallocate(tc_kern_array,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_lnv','tc_kern_array',ierr)

    call timer_clock(myself, 2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving ngwf_gradient_lnv'

  end subroutine ngwf_gradient_lnv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwf_gradient_batch(contra_grad, cov_grad, &  ! output
       fftbox_batch, lhxc_dbl, dbl_grid, rep, &
       ngwf_basis, proj_basis, hub_proj_basis, nl_projectors, hub, &
       fa_box_start, fa_start_in_box, batch_size, local_start, local_end, &
       coeff_mat, proj_coeff_mat, hub_proj_coeff_mat, scissor_coeff_mat, &
       cond_coeff_mat, ps_overlap, fftbox, cell, max_current_size, &
       is_force_call, val_ngwf_basis, val_rep, save_fftbox_batch, &
       save_fftbox_batch_cov, sav_grad, sav_grad_cov, kpt, dfdtau_dbl, &
       ireg, muext)

    !==========================================================================!
    ! This subroutine returns the contravariant and covariant NGWF             !
    ! gradients for the NGWFs of the current batch. It does this               !
    ! by applying the Hamiltonian operator to the functions accumulated        !
    ! in the fftbox of each batch, then applying kinetic energy                !
    ! preconditioning if required and then by extracting the                   !
    ! relevant ppds from the fftboxes and shaving their values                 !
    ! so that they are non-zero only within their spheres.                     !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! contra_grad (output) : Contravariant NGWF gradient in ppd representation !
    !    for the NGWFs of pub_my_proc_id.                                      !
    ! cov_grad (output)    : Covariant NGWF gradient in ppd representation     !
    !    for the NGWFs of pub_my_proc_id.                                      !
    ! fftbox_batch (in/out)                                                    !
    ! lhxc_dbl (input) : Total local potential in the double grid for the      !
    !    whole simulation cell.                                                !
    ! rep (input) : NGWF Representation (functions and matrices).              !
    ! ngwf_basis (input) : Function basis describing the NGWFs                 !
    ! proj_basis (input) : Function basis describing the nonlocal projectors   !
    ! hub_proj_basis (input) : Function basis describing Hubbard projectors    !
    ! nl_projectors (input) : Projector set for nonlocal pseudo projectors     !
    ! hub (input) : Storage for Hubbard Model information                      !
    ! fa_box_start (input) : Start positions of the FFTboxes for each function !
    ! fa_start_in_box (input) : Start positions of functions in FFTboxes       !
    ! batch_size (input) : Size of each batch of NGWFs                         !
    ! local_start (input) : First NGWF in this batch on this proc              !
    ! local_end (input) : Last NGWF in this batch on this proc                 !
    ! coeff_mat (input) : Coefficients of NGWFs to add to gradient             !
    ! proj_coeff_mat (input) : Coefficients of projectors to add to gradient   !
    ! hub_proj_coeff_mat (input) : Coefficients of Hubbard projectors          !
    ! cond_coeff_mat (input) : Coefficients of valence NGWFs for gradient of   !
    !    conduction NGWFs (optional - conduction NGWF optimisations only).     !
    ! val_ngwf_basis (input) : Valence NGWF basis (optional - cond NGWF only). !
    ! val_rep  (input) : NGWF_REP type for valence NGWFs (optional as above.   !
    ! save_fftbox_batch (inout) : fftbox to save contravariant contributions   !
    ! to the gradient for the mermin mode                                      !
    ! save_fftbox_batch (inout) : fftbox to save covariant contributions       !
    ! to the gradient for the mermin mode                                      !
    ! sav_grad (inout) : save the gradient contravariant part to calculate     !
    !                    muext contribution to the final gradient              !
    ! sav_grad_cov (inout) : save the gradient covariant part to calculate     !
    !                    muext contribution to the final gradient              !
    ! muext (input) : Chemical potential calculated for the NGWF cycle in the  !
    !                 mermin mode                                              !
    ! dfdtau_dbl (inout)     : Gradient of XC energy per unit volume wrt       !
    !                          KE density energy on double grid (optional)     !
    ! ireg (input) : Region counter for extracting info from rep.              !
    !--------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 26/12/2003                           !
    ! Modified for various speed enhancements by Nicholas Hine 2007-2009       !
    ! Modified for DFT+U by David O'Regan, 2009.                               !
    ! Modified to not use workspace_mod, Nicholas Hine, November 2009          !
    ! Tidied up, added more comments, reduced memory usage and made flow       !
    ! more logical, Nicholas Hine, October 2010                                !
    ! OpenMP parallelised by Karl Wilkinson and Nicholas Hine, May 2013        !
    ! Modified for embedding by Robert Charlton, 13/09/2017.                   !
    ! Hamiltonian application extracted by Robert Charlton, June 2018.         !
    !==========================================================================!

    use basis, only: basis_extract_function_from_box, basis_clean_function
    use cell_grid, only: GRID_INFO, cell_grid_extract_box
    use constants, only: DP, stdout, EDA_POLFRAGLOC_DEVEL, max_spins
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero, &
         data_fftbox_alloc, data_fftbox_dealloc, data_fftbox_axpy, &
         data_fftbox_copy, data_fftbox_scale
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use fragment_data, only: pub_super2fragid_on_proc
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_sum_fftbox_batch
    use geometry, only: POINT ! agrecokpt
    use hubbard_build, only: HUBBARD_MODEL
    use ngwf_representation, only: NGWF_REP
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_precond_real, pub_precond_recip, &
         pub_any_nl_proj, pub_hubbard, &
         pub_aug, pub_paw, &
         pub_cond_calculate, pub_debug_on_root, &
         pub_smooth_scheme, pub_num_spins, pub_threads_fftbox, pub_spin_fac, &
         pub_frag_counter, pub_eda_mode, pub_num_kpoints, &
         pub_xc_ke_density_required, &  ! JCW
         pub_emft, pub_active_region, pub_scissor_ngroups, pub_mermin
!$  use rundat, only: pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index, &
         sparse_create, sparse_copy, sparse_destroy, sparse_scale
    use sparse_embed, only: SPAM3_EMBED
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: batch_size
    type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(rep%nsub)
    type(FUNC_BASIS), intent(in) :: proj_basis(rep%nsub)
    type(FUNC_BASIS), intent(in) :: hub_proj_basis(rep%nsub)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(rep%nsub)
    type(HUBBARD_MODEL), intent(inout) :: hub
    integer, intent(in) :: fa_box_start(3,batch_size)
    integer, intent(in) :: fa_start_in_box(3,batch_size)
    integer, intent(in) :: local_start, local_end
    integer, intent(in) :: max_current_size
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FUNCTIONS), intent(inout) :: contra_grad
    type(FUNCTIONS), intent(inout) :: cov_grad
    type(FFTBOX_DATA), intent(inout) :: fftbox_batch(pub_num_spins,nmat,batch_size)
    type(GRID_INFO), intent(in) :: dbl_grid
    real(kind=DP), intent(in), target :: lhxc_dbl(dbl_grid%ld1, &
         dbl_grid%ld2, dbl_grid%max_group_slabs12, pub_num_spins, rep%nsub)
    type(SPAM3_EMBED), intent(in) :: coeff_mat(pub_num_spins,nmat)
    type(SPAM3_EMBED), intent(in) :: proj_coeff_mat(pub_num_spins,nmat)
    type(SPAM3_EMBED), intent(in) :: hub_proj_coeff_mat(pub_num_spins,nmat) !ddor
    type(SPAM3_EMBED), intent(in) :: scissor_coeff_mat(pub_num_spins,nmat) !ny
    type(SPAM3_EMBED), intent(in) :: cond_coeff_mat(pub_num_spins,nmat) !lr408
    type(SPAM3_EMBED), intent(in) :: ps_overlap
    logical, intent(in) :: is_force_call ! Check if this is cond task but
      ! function is called from forces_mod
    type(FUNC_BASIS), optional, intent(in) :: val_ngwf_basis(rep%nsub)
    type(NGWF_REP), optional, intent(in) :: val_rep
    ! JCW: Optional gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density, on double grid. A pointer will be assigned to this if
    ! JCW: present
    real(kind=DP), optional, target, intent(in) :: dfdtau_dbl(:,:,:,:)
    ! JCW: In general, we have
    ! JCW:    dfdtau_dbl(dbl_grid%ld1, &
    ! JCW:      dbl_grid%ld2, dbl_grid%max_group_slabs12, pub_num_spins)
    ! JCW: or
    ! JCW:    dfdtau_dbl(1,1,1,1)
    ! JCW: when it is unneeded.

    ! agrecokpt: argument to specify a k-point
    type(POINT), optional, intent(in) :: kpt
    ! rc2013: active region identifier
    integer, intent(in) :: ireg

    !ep: Mermin variables
    type(FFTBOX_DATA), optional, intent(inout) :: save_fftbox_batch(:,:,:)
    type(FFTBOX_DATA), optional, intent(inout) :: save_fftbox_batch_cov(:,:,:)
    type(FUNCTIONS), optional, intent(inout) :: sav_grad
    type(FUNCTIONS), optional, intent(inout) :: sav_grad_cov
    complex(kind=DP), allocatable, dimension(:,:,:) :: sav_coarse_work
    real(kind=DP), optional, intent(in) :: muext(max_spins,pub_num_kpoints)
    !$  integer, external :: omp_get_thread_num, omp_get_num_threads
    !ep

    ! Local Variables
    integer :: local_fa
    integer :: batch_count
    integer :: is, ik
    integer :: ierr
    integer :: idx_len
    integer, parameter :: fb_box = 1
    integer, parameter :: tc_fb_box = 2
    integer, parameter :: ham_box = 3
    integer, parameter :: tc_ham_box = 4
    integer, parameter :: fb_tilde_box = 5     ! jd: for polarisable embedding
    integer, parameter :: tc_fb_tilde_box = 6  ! jd: for polarisable embedding
    real(kind=DP) :: common_fac
    integer, allocatable, dimension(:) :: overlap_idx
    complex(kind=DP), allocatable, dimension(:,:,:) :: coarse_work
    ! JCW: Additional arrays required for tau-dependent contribution
    real(kind=DP), pointer, dimension(:,:,:,:) :: vtau_dbl
    ! rc2013: coefficient matrices for each subsystem to avoid producing
    ! array temporarys
    ! rc2013: Idea: we could use a (generalised) SPAM3_EMBED structure for
    ! these matrices, including everywhere they're used
    type(SPAM3), allocatable :: reg_coeff_mat(:,:), reg_proj_coeff_mat(:,:), &
         reg_cond_coeff_mat(:,:), &
         reg_hub_proj_coeff_mat(:,:), reg_scissor_coeff_mat(:,:)
    ! rc2013: should be allocatable
    type(FFTBOX_DATA), allocatable :: reg_fftbox_batch(:,:,:)
    real(kind=DP), pointer, dimension(:,:,:,:) :: loc_lhxc_dbl

    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt
    type(POINT) :: loc_kpt
    integer :: isub, imat

    integer :: size1,size2,size3


    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering ngwf_gradient_batch'

    call timer_clock('ngwf_gradient_batch',1)

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(ireg)%iscmplx

    ! agrecokpt
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    ! agrecocmplx: should be ok now, need to check!
    !jmecmplx
    !call utils_assert(.not. loc_cmplx, 'Error in&
    !     & ngwf_gradient_batch: not ready yet for complex NGWFs.')

    ! JCW: If dfdtau_dbl is present, and pub_xc_ke_density_required is .true.
    ! JCW: make vtau_dbl point to it, otherwise pointer remains nullified
    nullify(vtau_dbl)
    if (present(dfdtau_dbl)) then
       if (pub_xc_ke_density_required) then
          vtau_dbl => dfdtau_dbl
       else
          call utils_assert( size(dfdtau_dbl) == 1, &
               "Error in ngwf_gradient_batch: &
               &When pub_xc_ke_density_required is false, dfdtau_dbl must not &
               &be present, or have a size == 1.")
       end if
    end if

    ! rc2013: allocate regional matrices
    allocate(reg_coeff_mat(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','reg_coeff_mat',ierr)
    allocate(reg_proj_coeff_mat(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','reg_proj_coeff_mat',ierr)
    allocate(reg_cond_coeff_mat(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','reg_cond_coeff_mat',ierr)
    allocate(reg_hub_proj_coeff_mat(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','reg_hub_proj_coeff_mat',ierr)
    allocate(reg_scissor_coeff_mat(pub_num_spins,nmat),stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','reg_scissor_coeff_mat',ierr)

    ! rc2013: allocate storage for FFTboxes in each region
    allocate(reg_fftbox_batch(pub_num_spins, 4, batch_size), stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','reg_fftbox_batch',ierr)
    do batch_count = 1, batch_size
       do imat = 1,nmat
          do is = 1, pub_num_spins
             call data_fftbox_alloc(reg_fftbox_batch(is, imat, batch_count), &
                  fftbox%total_ld1, fftbox%total_ld2, &
                  fftbox%total_pt3, loc_cmplx)
          end do
       end do
    end do

    ! rc2013: begin iteration over regions
    do isub=1,size(ngwf_basis)

       ! pdh: overlap matrix index
       idx_len = sparse_index_length(rep%ngwf_overlap%m(isub,ireg))
       allocate(overlap_idx(idx_len),stat=ierr)
       call utils_alloc_check('ngwf_gradient_batch','overlap_idx',ierr)

       ! rc2013: extract the relevant coefficient matrix for this region
       do is=1,pub_num_spins
          do imat=1,nmat
             call sparse_create(reg_coeff_mat(is,imat),coeff_mat(is,imat)%m(isub,ireg))
             call sparse_copy(reg_coeff_mat(is,imat),coeff_mat(is,imat)%m(isub,ireg))
             if (pub_paw .or. (pub_any_nl_proj .and. (imat.ge.3))) then
                call sparse_create(reg_proj_coeff_mat(is,imat), &
                     proj_coeff_mat(is,imat)%m(isub,ireg))
                call sparse_copy(reg_proj_coeff_mat(is,imat), &
                     proj_coeff_mat(is,imat)%m(isub,ireg))
             end if
             if (pub_cond_calculate .and. (imat.ge.3) .and. .not. is_force_call) then
                call sparse_create(reg_cond_coeff_mat(is,imat), &
                     cond_coeff_mat(is,imat)%m(isub,ireg))
                call sparse_copy(reg_cond_coeff_mat(is,imat), &
                     cond_coeff_mat(is,imat)%m(isub,ireg))
             end if
             if (pub_hubbard .and. (imat.ge.3)) then
                call sparse_create(reg_hub_proj_coeff_mat(is,imat), &
                     hub_proj_coeff_mat(is,imat)%m(isub,ireg))
                call sparse_copy(reg_hub_proj_coeff_mat(is,imat), &
                     hub_proj_coeff_mat(is,imat)%m(isub,ireg))
             end if
             if (pub_scissor_ngroups > 0 .and. (imat.ge.3)) then
                call sparse_create(reg_scissor_coeff_mat(is,imat), &
                     scissor_coeff_mat(is,imat)%m(isub,ireg))
                call sparse_copy(reg_scissor_coeff_mat(is,imat), &
                     scissor_coeff_mat(is,imat)%m(isub,ireg))
             end if
          end do
       end do

       !ep: Mermin vars to zero
       if (pub_mermin) then
          call data_set_to_zero(reg_fftbox_batch(:,:,:))
          call data_set_to_zero(save_fftbox_batch(:,:,:))
          call data_set_to_zero(save_fftbox_batch_cov(:,:,:))
       else
          call data_set_to_zero(reg_fftbox_batch(:,:,:))
       end if
       !ep

       ! ndmh: calculate \sum_\beta \phi_\beta (r) * C^i_\alpha\beta in
       ! ndmh: FFT boxes for the four C^i matrices.
       ! ndmh: Once this is done:
       ! ndmh: Box 1 (fb_box) contains \sum_b Q^ab phi_b
       ! ndmh: Box 2 (tc_fb_box) contains \sum_b (QS)^a_b phi_b
       ! ndmh: Box 3 (ham_box) contains \sum_b K^ab phi_b
       ! ndmh: Box 4 (tc_ham_box) contains \sum_b (KS)^a_b phi_b
       ! jd:   Box 5 (fb_tilde_box) contains \sum_b Q~^ab phi0_b
       ! jd:   Box 6 (tc_fb_tilde_box) contains \sum_b (Q~S0)^a_b phi0_b
       ! jd:   where phi0 are NGWFs in the absence of embedding, S0 is the overlap
       ! jd:   in the absence of embedding, for definition of Q~ see [*].
       ! jd:   The contents of the boxes are determined by what is passed in
       !       coeff_mat (containing C^i_\alpha\beta).
       !       Boxes 3 and 4 are special in that the final answer is collected
       !       there.
       ! rc2013: get the information for this column
       ! rc2013: i.e. coeff matrices for all NGWF lists
       call sparse_generate_index(overlap_idx,rep%ngwf_overlap%m(isub,ireg))
       common_fac = 2.0_DP * cell%weight * pub_spin_fac
       call function_ops_sum_fftbox_batch(reg_fftbox_batch, fftbox, cell, &
            rep%ngwfs_on_grid(isub), ngwf_basis(isub), fa_box_start, batch_size, &
            local_start, local_end, overlap_idx, idx_len, reg_coeff_mat, &
            fb_box, tc_ham_box, common_fac)
       deallocate(overlap_idx,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_batch','overlap_idx',ierr)

       !ep: copy components in Mermin vars
       if (pub_mermin) then
            do batch_count = 1, local_end - local_start +1 !batch_size
               do is = 1, pub_num_spins
                  call data_fftbox_axpy(save_fftbox_batch(is,1,batch_count), &
                       reg_fftbox_batch(is,3,batch_count))
                  call data_fftbox_axpy(save_fftbox_batch_cov(is,1,batch_count), &
                       reg_fftbox_batch(is,4,batch_count))
               end do
            end do
       end if
       !ep

       ! rc2013: select the lhxc potential required for this region
       if (pub_emft .and. size(ngwf_basis) .gt. 1) then
          if( isub==ireg .and. ireg==pub_active_region) then
             ! rc2013: use the active potential
             loc_lhxc_dbl => lhxc_dbl(:,:,:,:,pub_active_region)
          else if(ireg==pub_active_region) then
             ! use the inactive potential for cross overlap terms
             loc_lhxc_dbl => lhxc_dbl(:,:,:,:,isub)
          else
             ! we're not in the active region; use this region's potential
             loc_lhxc_dbl => lhxc_dbl(:,:,:,:,ireg)
          end if
       else
          ! just use isub
          loc_lhxc_dbl => lhxc_dbl(:,:,:,:,isub)
       end if

       ! rc2013: now apply the Hamiltonian to this part of the batch
       ! rc2013: temporary fix: only apply COND parameters if needed
       ! jcap: if is_force_call is T, valence NGWF quantities and the
       ! cross overlap don't exist, so we need the other version
       if(pub_cond_calculate .and. .not. is_force_call) then
          call ngwf_gradient_build_batch(reg_fftbox_batch, &! output
               loc_lhxc_dbl, dbl_grid, ngwf_basis(ireg), proj_basis(isub), &
               hub_proj_basis(isub), nl_projectors(isub), hub, fa_box_start, fa_start_in_box, &
               batch_size, local_start, local_end, reg_coeff_mat, reg_proj_coeff_mat, &
               reg_cond_coeff_mat, ps_overlap%m(isub,ireg), fftbox, cell, &
               max_current_size, is_force_call,&
               val_ngwf_basis(isub), val_rep%ngwfs_on_grid(isub), &
               rep%ngwf_cross_overlap%m(isub,ireg), kpt=loc_kpt, dfdtau_dbl = dfdtau_dbl, &
               hub_proj_coeff_mat=reg_hub_proj_coeff_mat, &
               hub_overlap_t=rep%hub_overlap_t%p )
       else
          ! rc2013: EMBED_FIX: Hubbard here will only work with 1 subsystem!
          call ngwf_gradient_build_batch(reg_fftbox_batch, &! output
               loc_lhxc_dbl, dbl_grid, ngwf_basis(ireg), proj_basis(isub), &
               hub_proj_basis(isub), nl_projectors(isub), hub, fa_box_start, fa_start_in_box, &
               batch_size, local_start, local_end, reg_coeff_mat, reg_proj_coeff_mat, &
               reg_cond_coeff_mat, ps_overlap%m(isub,ireg), fftbox, cell, &
               max_current_size, is_force_call,&
               kpt=loc_kpt, dfdtau_dbl = dfdtau_dbl, &
               hub_proj_coeff_mat=reg_hub_proj_coeff_mat, &
               hub_overlap_t=rep%hub_overlap_t%p )
               !hub_overlap_t=rep%hub_overlap_t%m(isub,ireg) )
       end if

       ! SSSSSS SCISSOR SHIFT SSSSSSSSSSSSSSSSSSSSSSSSSSSSS
       if (pub_scissor_ngroups > 0) then
          idx_len = sparse_index_length(rep%ngwf_overlap%m(isub,ireg))
          allocate(overlap_idx(idx_len),stat=ierr)
          call utils_alloc_check('ngwf_gradient_batch','overlap_idx',ierr)

          call sparse_generate_index(overlap_idx,rep%ngwf_overlap%m(isub,ireg))
          common_fac = 2.0_DP * cell%weight * pub_spin_fac
          call function_ops_sum_fftbox_batch(reg_fftbox_batch, fftbox, cell, &
               rep%ngwfs_on_grid(isub), ngwf_basis(isub), fa_box_start, batch_size, &
               local_start, local_end, overlap_idx, idx_len, reg_scissor_coeff_mat, &
               ham_box, tc_ham_box, common_fac)

          deallocate(overlap_idx,stat=ierr)
          call utils_dealloc_check('ngwf_gradient_batch','overlap_idx',ierr)
       end if
       ! SSSSSS END SCISSOR SHIFT SSSSSSSSSSSSSSSSSSSSSSSSS

       ! rc2013: add this into the full batch
       do batch_count = 1, local_end - local_start +1 !batch_size
          do imat = 1, nmat
             do is = 1, pub_num_spins
                call data_fftbox_axpy(fftbox_batch(is,imat,batch_count), &
                     reg_fftbox_batch(is,imat,batch_count))
             end do
          end do
       end do

       ! rc2013: deallocate temporary matrices
       do is=1,pub_num_spins
          do imat=1,nmat
             call sparse_destroy(reg_coeff_mat(is,imat))
             if (pub_paw .or. (pub_any_nl_proj .and. (imat.ge.3))) &
                  call sparse_destroy(reg_proj_coeff_mat(is,imat))
             if (pub_hubbard .and. (imat.ge.3)) &
                  call sparse_destroy(reg_hub_proj_coeff_mat(is,imat))
             if (pub_cond_calculate .and. (imat.ge.3) .and. .not. is_force_call) &
                  call sparse_destroy(reg_cond_coeff_mat(is,imat))
             if (pub_scissor_ngroups > 0 .and. (imat.ge.3)) &
                  call sparse_destroy(reg_scissor_coeff_mat(is,imat))
          end do
       end do
    end do
    deallocate(reg_coeff_mat,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','reg_coeff_mat',ierr)
    deallocate(reg_proj_coeff_mat,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','reg_proj_coeff_mat',ierr)
    deallocate(reg_hub_proj_coeff_mat,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','reg_hub_proj_coeff_mat',ierr)
    deallocate(reg_cond_coeff_mat,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','reg_cond_coeff_mat',ierr)
    deallocate(reg_scissor_coeff_mat,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','reg_scissor_coeff_mat',ierr)

    nullify(loc_lhxc_dbl)

    ! pdh: deallocate workspace
    do batch_count = batch_size, 1, -1
       do imat = nmat, 1, -1
          do is = pub_num_spins, 1, -1
             call data_fftbox_dealloc(reg_fftbox_batch(is, imat, batch_count))
          end do
       end do
    end do
    deallocate(reg_fftbox_batch, stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','reg_fftbox_batch',ierr)

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(is,local_fa,batch_count,coarse_work, sav_coarse_work, &
!$OMP       ierr) &
!$OMP SHARED(fa_start_in_box,pub_precond_real,pub_aug, &
!$OMP      pub_eda_mode,pub_super2fragid_on_proc,pub_frag_counter, &
!$OMP      pub_smooth_scheme,pub_precond_recip,local_start,local_end, &
!$OMP      cell,fftbox,fftbox_batch,ngwf_basis,cov_grad,contra_grad, &
!$OMP      pub_num_spins,pub_threads_fftbox,pub_threads_num_fftboxes, &
!$OMP      loc_cmplx, pub_mermin,save_fftbox_batch,sav_grad, &
!$OMP      save_fftbox_batch_cov,sav_grad_cov, stdout, muext, ireg)

    ! ndmh: allocate generic complex FFTbox workspace
    allocate(coarse_work(fftbox%total_ld1,fftbox%total_ld2,&
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','coarse_work',ierr)

    !ep: Mermin vars
    if (pub_mermin) then
       allocate(sav_coarse_work(fftbox%total_ld1,fftbox%total_ld2,&
            fftbox%total_pt3),stat=ierr)
       call utils_alloc_check('ngwf_gradient_batch','sav_coarse_work',ierr)
    end if
    !ep



    ! RRRRRRR RECIPROCAL SPACE PRECONDITIONING RRRRRRRRRRRRRRRR
    if (pub_precond_recip) then
       do is=1,pub_num_spins
!$OMP DO
          do local_fa=local_start,local_end,2
             batch_count = local_fa - local_start + 1

             ! agrecocmplx: in complex case, we can only work with one
             ! fftbox_batch at a time
             ! real case
             if (.not.loc_cmplx) then
                ! Copy the covariant gradient into complex workspace array
                if (local_fa < local_end) then
                   coarse_work = &
                        cmplx(fftbox_batch(is,tc_ham_box,batch_count)%d, &
                        fftbox_batch(is,tc_ham_box,batch_count+1)%d,kind=DP)
                   !ep: mermin switch present in all the reciprocal space
                   !preconditioning steps to save covariant and contravariant
                   !gradient components to calculate muext contributions
                   if (pub_mermin) then
                      sav_coarse_work =&
                           cmplx(save_fftbox_batch_cov(is,1,batch_count)%d, &
                           save_fftbox_batch_cov(is,1,batch_count+1)%d,kind=DP)
                   endif
                else
                   coarse_work = &
                        cmplx(fftbox_batch(is,tc_ham_box,batch_count)%d, &
                        0.0_DP,kind=DP)
                   !ep
                   if (pub_mermin) then
                      sav_coarse_work =&
                           cmplx(save_fftbox_batch_cov(is,1,batch_count)%d, &
                           0.0_DP,kind=DP)
                   endif
                end if

                ! Forward FFT the covariant gradient to reciprocal space
                call fourier_apply_box('Coarse','Forward',coarse_work, &
                                       omp=pub_threads_fftbox)
                !ep
                if (pub_mermin) then
                   call fourier_apply_box('Coarse','Forward',sav_coarse_work, &
                                          omp=pub_threads_fftbox)
                end if

                ! Apply kinetic energy preconditioning to covariant gradient
                call ngwf_gradient_precond_recip(coarse_work,fftbox)

                !ep
                if (pub_mermin) then
                    call ngwf_gradient_precond_recip(sav_coarse_work,fftbox)
                end if


                ! Backward FFT the covariant gradient to real space
                call fourier_apply_box('Coarse','Backward',coarse_work, &
                                       omp=pub_threads_fftbox)

                !ep
                if (pub_mermin) then
                   call fourier_apply_box('Coarse','Backward',sav_coarse_work, &
                                          omp=pub_threads_fftbox)
                end if

                ! Copy the preconditioned covariant gradient out of the complex
                ! workspace array
                fftbox_batch(is,tc_ham_box,batch_count)%d(:,:,:) = &
                     real(coarse_work,kind=DP)
                !ep
                if (pub_mermin) then
                   save_fftbox_batch_cov(is,1,batch_count)%d(:,:,:) = &
                        real(sav_coarse_work,kind=DP)
                end if

                if (local_fa < local_end) then
                     fftbox_batch(is,tc_ham_box,batch_count+1)%d(:,:,:) = &
                     aimag(coarse_work)
                end if
                !ep
                if (pub_mermin) then
                   if (local_fa < local_end) then
                      save_fftbox_batch_cov(is,1,batch_count+1)%d(:,:,:) = &
                      aimag(sav_coarse_work)
                   end if
                end if
             ! complex case
             else
                ! Forward FFT the covariant gradient to reciprocal space
                call fourier_apply_box('Coarse','Forward', &
                        fftbox_batch(is,tc_ham_box,batch_count)%z, &
                        omp=pub_threads_fftbox)
                !ep
                if (pub_mermin) then
                   call fourier_apply_box('Coarse','Forward', &
                        save_fftbox_batch_cov(is,1,batch_count)%z, &
                        omp=pub_threads_fftbox)

                end if

                ! Apply kinetic energy preconditioning to covariant gradient
                call ngwf_gradient_precond_recip( &
                        fftbox_batch(is,tc_ham_box,batch_count)%z,fftbox)

                !ep
                if (pub_mermin) then
                   call ngwf_gradient_precond_recip( &
                        save_fftbox_batch_cov(is,1,batch_count)%z,fftbox)

                end if

                ! Backward FFT the covariant gradient to real space and
                call fourier_apply_box('Coarse','Backward', &
                        fftbox_batch(is,tc_ham_box,batch_count)%z, &
                        omp=pub_threads_fftbox)
                !ep
                if (pub_mermin) then
                   call fourier_apply_box('Coarse','Backward', &
                        save_fftbox_batch_cov(is,1,batch_count)%z, &
                        omp=pub_threads_fftbox)
                end if

                ! agrecocmplx: same for function 2 if we have it
                if (local_fa < local_end) then
                   ! Forward FFT the covariant gradient to reciprocal space
                   call fourier_apply_box('Coarse','Forward', &
                           fftbox_batch(is,tc_ham_box,batch_count+1)%z, &
                           omp=pub_threads_fftbox)
                   !ep
                   if (pub_mermin) then
                      call fourier_apply_box('Coarse','Forward', &
                           save_fftbox_batch_cov(is,1,batch_count+1)%z, &
                           omp=pub_threads_fftbox)
                   end if

                   ! Apply kinetic energy preconditioning to covariant gradient
                   call ngwf_gradient_precond_recip( &
                           fftbox_batch(is,tc_ham_box,batch_count+1)%z,fftbox)
                   !ep
                   if (pub_mermin) then
                      call ngwf_gradient_precond_recip( &
                           save_fftbox_batch_cov(is,1,batch_count+1)%z,fftbox)
                   end if

                   ! Backward FFT the covariant gradient to real space and
                   call fourier_apply_box('Coarse','Backward', &
                           fftbox_batch(is,tc_ham_box,batch_count+1)%z, &
                           omp=pub_threads_fftbox)
                   !ep
                   if (pub_mermin) then
                       call fourier_apply_box('Coarse','Backward', &
                               save_fftbox_batch_cov(is,1,batch_count+1)%z, &
                               omp=pub_threads_fftbox)
                   end if

                end if
             end if
          end do
!$OMP END DO
       end do
    end if
    ! RRRRRRR RECIPROCAL SPACE PRECONDITIONING RRRRRRRRRRRRRRRR

    ! ndmh: APPLY REAL SPACE PRECONDITIONING, APPLY SPIN-AVERAGING,
    ! ndmh: THEN EXTRACT AND SHAVE GRADIENT
!$OMP DO
    do local_fa=local_start,local_end
       batch_count = local_fa - local_start + 1

       if (pub_precond_real) then
          do is=1,pub_num_spins
             ! cks: precondition the covariant gradient in real space
             call ngwf_grad_precond_rspace(  &
                  fftbox_batch(is,tc_ham_box,batch_count), & ! in/out
                  fa_start_in_box(1,batch_count), &
                  fa_start_in_box(2,batch_count), &
                  fa_start_in_box(3,batch_count), &
                  ngwf_basis(ireg)%tight_boxes(local_fa)%tight_pts1, &
                  ngwf_basis(ireg)%tight_boxes(local_fa)%tight_pts2, &
                  ngwf_basis(ireg)%tight_boxes(local_fa)%tight_pts3,fftbox,cell)
             !ep reciprocal space preconditionic for the saved fftbox
             ! used in the mermin method
             if (pub_mermin) then
                 call ngwf_grad_precond_rspace(  &
                      save_fftbox_batch_cov(is,1,batch_count), & ! in/out
                      fa_start_in_box(1,batch_count), &
                      fa_start_in_box(2,batch_count), &
                      fa_start_in_box(3,batch_count), &
                      ngwf_basis(ireg)%tight_boxes(local_fa)%tight_pts1, &
                      ngwf_basis(ireg)%tight_boxes(local_fa)%tight_pts2, &
                      ngwf_basis(ireg)%tight_boxes(local_fa)%tight_pts3,fftbox,cell)

             end if
          end do
       end if

       !ep: calculate muext contributions to the gradient
       !if (pub_mermin) then
       !   call utils_assert(present(muext), &
       !        'Mu needs to be passed if pub_mermin is true')
       !end if
       if (present(muext)) then
          do ik = 1, pub_num_kpoints
             do is=1,pub_num_spins
                call data_fftbox_scale(save_fftbox_batch(is,1,batch_count), &
                     muext(is,ik))
                call data_fftbox_scale(save_fftbox_batch_cov(is,1,batch_count), &
                     muext(is,ik))
                call data_fftbox_axpy(fftbox_batch(is,ham_box,batch_count), &
                     save_fftbox_batch(is,1,batch_count), 1.0_DP)
                call data_fftbox_axpy(fftbox_batch(is,tc_ham_box,batch_count), &
                     save_fftbox_batch_cov(is,1,batch_count), 1.0_DP)
             end do
          end do
       end if

       ! ndmh: Average gradients over spins for spin-polarised calculations.
       if (pub_num_spins == 2) then

          ! ny: use fft_box axpy routine
          call data_fftbox_axpy(fftbox_batch(1,ham_box,batch_count), &
               fftbox_batch(2,ham_box,batch_count))
          call data_fftbox_axpy(fftbox_batch(1,tc_ham_box,batch_count), &
               fftbox_batch(2,tc_ham_box,batch_count))

       end if

       !ep
       if (pub_mermin .and. pub_num_spins == 2) then
          call data_fftbox_axpy(save_fftbox_batch(1,1,batch_count), &
               save_fftbox_batch(2,1,batch_count))
          call data_fftbox_axpy(save_fftbox_batch_cov(1,1,batch_count), &
               save_fftbox_batch_cov(2,1,batch_count))
       end if

       if (pub_smooth_scheme .ne. 'NONE') then
          ! smmd: smooth the covariant gradient in real space
          call ngwf_grad_smooth_rspace(  &
               fftbox_batch(1,tc_ham_box,batch_count), & ! in/out
               fa_start_in_box(1,batch_count), &
               fa_start_in_box(2,batch_count), &
               fa_start_in_box(3,batch_count), &
               ngwf_basis(ireg)%tight_boxes(local_fa), &
               ngwf_basis(ireg)%spheres(local_fa)%centre, &
               ngwf_basis(ireg)%spheres(local_fa)%radius,fftbox,cell)
          !ep
          if (pub_mermin) then
             call ngwf_grad_smooth_rspace(  &
                  save_fftbox_batch_cov(1,1,batch_count), & ! in/out
                  fa_start_in_box(1,batch_count), &
                  fa_start_in_box(2,batch_count), &
                  fa_start_in_box(3,batch_count), &
                  ngwf_basis(ireg)%tight_boxes(local_fa), &
                  ngwf_basis(ireg)%spheres(local_fa)%centre, &
                  ngwf_basis(ireg)%spheres(local_fa)%radius,fftbox,cell)
          end if
       endif

       !ep: shaving mermin grad vars
       if (pub_mermin) then
          call basis_extract_function_from_box(sav_grad, &
               save_fftbox_batch(1,1,batch_count), &
               ngwf_basis(ireg)%spheres(local_fa),&
               ngwf_basis(ireg)%tight_boxes(local_fa), &
                fa_start_in_box(1,batch_count), &
               fa_start_in_box(2,batch_count), fa_start_in_box(3,batch_count), &
               ngwf_basis(ireg)%spheres(local_fa)%offset,cell,fftbox)

          call basis_extract_function_from_box(sav_grad_cov, &
               save_fftbox_batch_cov(1,1,batch_count), &
               ngwf_basis(ireg)%spheres(local_fa),&
               ngwf_basis(ireg)%tight_boxes(local_fa), &
               fa_start_in_box(1,batch_count), &
               fa_start_in_box(2,batch_count), fa_start_in_box(3,batch_count), &
               ngwf_basis(ireg)%spheres(local_fa)%offset,cell,fftbox)

       end if

       ! cks: shaving - stage 1 / extract ppds from the FFT box
       ! extract the contravariant gradient from ham_box to PPDs
       call basis_extract_function_from_box(contra_grad, &
            fftbox_batch(1,ham_box,batch_count), &
            ngwf_basis(ireg)%spheres(local_fa),&
            ngwf_basis(ireg)%tight_boxes(local_fa), fa_start_in_box(1,batch_count), &
            fa_start_in_box(2,batch_count), fa_start_in_box(3,batch_count), &
            ngwf_basis(ireg)%spheres(local_fa)%offset,cell,fftbox)

       ! extract the covariant gradient from ham_box to PPDs
       call basis_extract_function_from_box(cov_grad, &
            fftbox_batch(1,tc_ham_box,batch_count), &
            ngwf_basis(ireg)%spheres(local_fa),&
            ngwf_basis(ireg)%tight_boxes(local_fa), fa_start_in_box(1,batch_count), &
            fa_start_in_box(2,batch_count), fa_start_in_box(3,batch_count), &
            ngwf_basis(ireg)%spheres(local_fa)%offset, cell, fftbox)

       ! cks: shaving - stage 2 / zero points outside NGWF sphere in PPD rep.
       ! smmd: if smoothing scheme applied previously, no needs to clean the
       ! smmd: covariant gradient (nor the contravariant gradient which is
       ! smmd: always used within contra.DOT.cov expressions)
       if (cell%n_pts > 1 .and. pub_smooth_scheme .eq. 'NONE') then
          call basis_clean_function(contra_grad, &
               ngwf_basis(ireg)%spheres(local_fa), cell,fftbox)
          call basis_clean_function(cov_grad, &
               ngwf_basis(ireg)%spheres(local_fa), cell,fftbox)
          !ep
          if (pub_mermin) then
             call basis_clean_function(sav_grad, &
                  ngwf_basis(ireg)%spheres(local_fa), cell,fftbox)
             call basis_clean_function(sav_grad_cov, &
                  ngwf_basis(ireg)%spheres(local_fa), cell,fftbox)
          endif
       endif

       ! mjsp: EDA fragment-wise polarisation: prevent optimisation of specific frozen NGWFs:
       if (pub_eda_mode == EDA_POLFRAGLOC_DEVEL) then

          ! check if this NGWF is on the fragment(s) that the gradient is to be
          ! set to 0 for. (i.e. check if this NGWF is not being polarised).
          ! note: pub_frag_counter=0 is reserved for polarisation of all NGWFs.
          if ((pub_frag_counter .ne. 0) .and. &
            (pub_super2fragid_on_proc(local_fa) /= pub_frag_counter)) then

             ! set the gradient to 0
             ! agrecocmplx
             if (loc_cmplx) then  ! if complex NGWFs
                contra_grad%z(ngwf_basis(ireg)%spheres(local_fa)%offset: &
                     ngwf_basis(ireg)%spheres(local_fa)%offset - 1 &
                     + ngwf_basis(ireg)%spheres(local_fa)%n_ppds_sphere*cell%n_pts) = 0.0_DP
                cov_grad%z(ngwf_basis(ireg)%spheres(local_fa)%offset: &
                     ngwf_basis(ireg)%spheres(local_fa)%offset - 1 &
                     + ngwf_basis(ireg)%spheres(local_fa)%n_ppds_sphere*cell%n_pts) = 0.0_DP
             else
                contra_grad%d(ngwf_basis(ireg)%spheres(local_fa)%offset: &
                     ngwf_basis(ireg)%spheres(local_fa)%offset - 1 &
                     + ngwf_basis(ireg)%spheres(local_fa)%n_ppds_sphere*cell%n_pts) = 0.0_DP
                cov_grad%d(ngwf_basis(ireg)%spheres(local_fa)%offset: &
                     ngwf_basis(ireg)%spheres(local_fa)%offset - 1 &
                     + ngwf_basis(ireg)%spheres(local_fa)%n_ppds_sphere*cell%n_pts) = 0.0_DP
             endif

          end if

       end if

    end do
!$OMP END DO
    ! ndmh: END APPLY REAL SPACE PRECONDITIONING, APPLY SPIN-AVERAGING,
    ! ndmh: THEN EXTRACT AND SHAVE GRADIENT

    ! ndmh: deallocate generic workspace
    deallocate(coarse_work,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','coarse_work',ierr)
!$OMP END PARALLEL

    do batch_count = 1, local_end - local_start +1 !batch_size
       do imat = 1, nmat
          do is = 1, pub_num_spins
             size1=size(fftbox_batch(is,imat,batch_count)%d(:,1,1))
             size2=size(fftbox_batch(is,imat,batch_count)%d(1,:,1))
             size3=size(fftbox_batch(is,imat,batch_count)%d(1,1,:))
          end do
       end do
    end do

    ! JCW: Nullify vtau_dbl pointer, if associated
    if (associated(vtau_dbl)) nullify(vtau_dbl)

    call timer_clock('ngwf_gradient_batch',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving ngwf_gradient_batch'

  end subroutine ngwf_gradient_batch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_gradient_kinetic(out1,out2,in1,in2,zwork_box,fftbox, &
                  kpt) ! agrecokpt

    !=========================================================!
    ! This subroutine applies the kinetic energy operator to  !
    ! an NGWF in an FFTbox.                                   !
    !---------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris as          !
    ! ngwf_gradient_kinloc_pairbox on 15/6/2001.              !
    ! Rewritten by Arash A. Mostofi and modified to use       !
    ! complex-to-complex FFTs on April 2003                   !
    ! Modified by Andrea Greco to work with complex NGWFs and !
    ! k-points, April 2016.                                   !
    !=========================================================!

    use datatypes, only: FFTBOX_DATA ! agrecocmplx
    use constants, only: stdout
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box_pair, fourier_apply_box
    use geometry, only: POINT ! agrecokpt
    use kinetic, only: kinetic_apply_on_box
    ! agrecokpt
    use rundat, only: pub_threads_fftbox, pub_kpoint_method, &
            pub_debug_on_root
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(in) :: fftbox
    ! agrecocmplx: use FFTBOX_DATA types
    !real(kind=DP), intent(out)   :: out1(fftbox%total_ld1, &
    !     fftbox%total_ld2,fftbox%total_pt3)
    !real(kind=DP), intent(out)   :: out2(fftbox%total_ld1, &
    !     fftbox%total_ld2,fftbox%total_pt3)
    !real(kind=DP), intent(inout) :: in1(fftbox%total_ld1, &
    !     fftbox%total_ld2,fftbox%total_pt3)
    !real(kind=DP), intent(inout) :: in2(fftbox%total_ld1, &
    !     fftbox%total_ld2,fftbox%total_pt3)
    type(FFTBOX_DATA), intent(inout) :: out1
    type(FFTBOX_DATA), intent(inout) :: out2
    type(FFTBOX_DATA), intent(inout) :: in1
    type(FFTBOX_DATA), intent(inout) :: in2
    complex(kind=DP), intent(out) :: zwork_box(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    ! agrecokpt
    type(POINT), optional, intent(in) :: kpt

    ! Local variables
    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt: k-point in cartesian coordinates
    real(kind=DP) :: kcart(3)

    ! agrecocmplx
    loc_cmplx = in1%iscmplx
    ! agrecokpt
    if (present(kpt)) then
       kcart(1) = kpt%x
       kcart(2) = kpt%y
       kcart(3) = kpt%z
    else
       kcart(:) = 0.0_DP
    end if

    ! agrecokpt: currently only KP method is being implemented
    call utils_assert(pub_kpoint_method == 'KP', &
         'Subroutine ngwf_gradient_kinetic currently supports&
         & only KP method for BZ sampling')

    ! agrecocmplx: in complex case, can transform only one at a time
    if (loc_cmplx) then

       ! Forward Fourier transform 1 to reciprocal space
       call fourier_apply_box('Coarse', 'Forward', in1%z, &
            gspc=zwork_box, omp=pub_threads_fftbox)

       ! Apply kinetic operator to 1in reciprocal space
       call kinetic_apply_on_box(zwork_box,fftbox)

       ! Backward Fourier transform 1 to real space
       call fourier_apply_box('Coarse', 'Backward', out1%z, &
            gspc=zwork_box, omp=pub_threads_fftbox)

       ! Forward Fourier transform 2 to reciprocal space
       call fourier_apply_box('Coarse', 'Forward', in2%z, &
            gspc=zwork_box, omp=pub_threads_fftbox)

       ! Apply kinetic operator to 1in reciprocal space
       call kinetic_apply_on_box(zwork_box,fftbox)

       ! Backward Fourier transform 1 to real space
       call fourier_apply_box('Coarse', 'Backward', out2%z, &
            gspc=zwork_box, omp=pub_threads_fftbox)

    ! real case: transform 2 at a time
    else
       ! Forward Fourier transform to reciprocal space
       call fourier_apply_box_pair('Coarse','Forward',in1%d,in2%d,zwork_box, &
            omp=pub_threads_fftbox)

       ! Apply kinetic operator in reciprocal space
       call kinetic_apply_on_box(zwork_box,fftbox)

       ! Backward Fourier transform to real space
       call fourier_apply_box_pair('Coarse','Backward',out1%d,out2%d,zwork_box, &
            omp=pub_threads_fftbox)

    end if

    ! agrecokpt: add k-point dependent terms in KP method
    ! only if a non Gamma point is specified
    if ((pub_kpoint_method == 'KP') .and. any(kcart/=0.0_DP)) then
       if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
          &Adding KP terms in ngwf_gradient_kinetic'

       ! agrecokpt: need to update both out1 and out2?
       call ngwf_gradient_update_kinet_kp(out1, in1, &
               zwork_box, fftbox, kcart)

       call ngwf_gradient_update_kinet_kp(out2, in2, &
               zwork_box, fftbox, kcart)

    end if

    return

  end subroutine ngwf_gradient_kinetic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine ngwf_gradient_local(data_box1, data_box2, pot_dbl1, pot_dbl2, &
       dwork_box_dbl,coarse_work,fine_work,fftbox)

    !=========================================================!
    ! This subroutine applies the local potential to a pair   !
    ! of NGWFs in FFTboxes.                                   !
    !---------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris as          !
    ! ngwf_gradient_kinloc_pairbox on 15/6/2001.              !
    ! Rewritten by Arash A. Mostofi and modified to use       !
    ! complex-to-complex FFTs in January 2003.                !
    ! Modified by Chris-Kriton Skylaris on 28/12/2003 to run  !
    ! with the parallel version of ONETEP.                    !
    ! Modified by Andrea Greco to allow compatibility with    !
    ! complex NGWFs, April 2016.                              !
    !=========================================================!

    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_interpolate, fourier_filter
    use rundat, only: pub_threads_per_fftbox
    use utils, only: utils_use_var
    ! agrecocmplx
    use datatypes, only: FFTBOX_DATA

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(in) :: fftbox
    ! agrecocmplx: use FFTBOX_DATA types
    !real(kind=DP), intent(inout) :: data_box1(fftbox%total_ld1, &
    !     fftbox%total_ld2,fftbox%total_pt3)
    !real(kind=DP), intent(inout) :: data_box2(fftbox%total_ld1, &
    !     fftbox%total_ld2,fftbox%total_pt3)
    type(FFTBOX_DATA), intent(inout) :: data_box1
    type(FFTBOX_DATA), intent(inout) :: data_box2
    real(kind=DP), intent(in)    :: pot_dbl1(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl)
    real(kind=DP), intent(in)    :: pot_dbl2(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl)
    real(kind=DP), intent(out)   :: dwork_box_dbl(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl,2)
    complex(kind=DP), intent(out) :: fine_work(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl)
    complex(kind=DP), intent(out) :: coarse_work(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    ! agrecocmplx
    logical :: loc_cmplx

    loc_cmplx = data_box1%iscmplx

    ! agrecocmplx: in complex case, filter/interpolate 1 at a time
    if (loc_cmplx) then
       ! First
       call fourier_interpolate(coarse_work,fine_work,data_box1%z)

!!$OMP PARALLEL WORKSHARE NUM_THREADS(pub_threads_per_fftbox)
       fine_work(:,:,:) = fine_work(:,:,:) * pot_dbl1
!!$OMP END PARALLEL WORKSHARE

       call fourier_filter(coarse_work,fine_work,fine_work,data_box1%z)

       ! Second
       call fourier_interpolate(coarse_work,fine_work,data_box2%z)

!!$OMP PARALLEL WORKSHARE NUM_THREADS(pub_threads_per_fftbox)
       fine_work(:,:,:) = fine_work(:,:,:) * pot_dbl2
!!$OMP END PARALLEL WORKSHARE

       call fourier_filter(coarse_work,fine_work,fine_work,data_box2%z)

    ! real case: fourier/interpolate 2 at a time
    else

       call fourier_interpolate(coarse_work,fine_work,data_box1%d,data_box2%d, &
            dwork_box_dbl(:,:,:,1),dwork_box_dbl(:,:,:,2))

!!$OMP PARALLEL WORKSHARE NUM_THREADS(pub_threads_per_fftbox)
       dwork_box_dbl(:,:,:,1) = dwork_box_dbl(:,:,:,1) * pot_dbl1
       dwork_box_dbl(:,:,:,2) = dwork_box_dbl(:,:,:,2) * pot_dbl2
!!$OMP END PARALLEL WORKSHARE

       call fourier_filter(coarse_work,fine_work,dwork_box_dbl(:,:,:,1), &
            dwork_box_dbl(:,:,:,2),data_box1%d, data_box2%d)

    end if

    ! jd: Silence var unused warning when in non-OMP compiles
    call utils_use_var(pub_threads_per_fftbox)

  end subroutine ngwf_gradient_local

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwf_gradient_div_locpot_grad(data_box1, data_box2, &
       pot_dbl1, pot_dbl2,coarse_work,fine_work,&
       dwork_grad_box,dwork_grad_box_dbl,zwork_grad_box,&
       fftbox)

    !=========================================================!
    ! Evaluates the divergence of a local potential applied   !
    ! to the gradient of a pair of functions in fftboxes, i.e.!
    !      \nabla . ( Vloc \nabla f )                         !
    ! For evaluation of the gradient of the energy wrt NGWFs, !
    ! f is                                                    !
    !   1. \sum_b K^ab phi_b (contravariant part)             !
    !   2. \sum_b (KS)_a^b phi_b (covariant part)             !
    ! for evaluating the contribution to the NGWF gradient    !
    ! from the KE-density-dependent part of the XC potential  !
    ! Vloc is dfdtau_dbl and the result should be multiplied  !
    ! by -0.5_DP when added to the NGWF gradient              !
    !---------------------------------------------------------!
    ! TODO Describe arguments                                 !
    !---------------------------------------------------------!
    ! Based on ngwf_gradient_local                            !
    ! Created by James C. Womack, 02/2016                     !
    !---------------------------------------------------------!
    ! From ngwf_gradient_local:                               !
    ! Originally written by Chris-Kriton Skylaris as          !
    ! ngwf_gradient_kinloc_pairbox on 15/6/2001.              !
    ! Rewritten by Arash A. Mostofi and modified to use       !
    ! complex-to-complex FFTs in January 2003.                !
    ! Modified by Chris-Kriton Skylaris on 28/12/2003 to run  !
    ! with the parallel version of ONETEP.                    !
    !=========================================================!

    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box_pair, fourier_interpolate, fourier_filter
    use kinetic, only: kinetic_grad_on_box, kinetic_grad_on_cart_box
    use rundat, only: pub_threads_fftbox, pub_threads_per_fftbox
    use utils, only: utils_use_var

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(inout) :: data_box1(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    real(kind=DP), intent(inout) :: data_box2(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    real(kind=DP), intent(in)    :: pot_dbl1(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl)
    real(kind=DP), intent(in)    :: pot_dbl2(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl)
    complex(kind=DP), intent(inout) :: fine_work(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl)
    complex(kind=DP), intent(inout) :: coarse_work(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    real(kind=DP), intent(inout) :: dwork_grad_box(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3,3,2)
    real(kind=DP), intent(inout) :: dwork_grad_box_dbl(fftbox%total_ld1_dbl,&
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl,3,2)
    complex(kind=DP), intent(inout) :: zwork_grad_box(fftbox%total_ld1,&
         fftbox%total_ld2, fftbox%total_pt3,3,2)

    ! Local variables
    integer :: dim


    ! 1. Apply grad operator to f
    ! Forward Fourier transform to reciprocal space
    call fourier_apply_box_pair('Coarse','Forward',data_box1,data_box2,&
         zwork_grad_box(:,:,:,1,2), &
         omp=pub_threads_fftbox)
    ! use one unused complex FFTbox in zwork_grad_box to store transformed data,
    ! avoiding need for additional allocatable array

    ! Apply kinetic operator in reciprocal space
    call kinetic_grad_on_box(zwork_grad_box(:,:,:,1,2),&
         zwork_grad_box(:,:,:,:,1),fftbox)

    ! Backward Fourier transform to real space
    do dim = 1, 3
       call fourier_apply_box_pair('Coarse','Backward',dwork_grad_box(:,:,:,dim,1),&
            dwork_grad_box(:,:,:,dim,2),zwork_grad_box(:,:,:,dim,1), &
            omp=pub_threads_fftbox)
    end do
    ! ...have \nabla f (x, y, z components) in dwork_grad_box(:,:,:,1:3,1:2)

    ! 2. Apply potential to \nabla f
    do dim = 1, 3
       call fourier_interpolate(coarse_work,fine_work,&
            dwork_grad_box(:,:,:,dim,1),dwork_grad_box(:,:,:,dim,2), &
            dwork_grad_box_dbl(:,:,:,dim,1),dwork_grad_box_dbl(:,:,:,dim,2))

       dwork_grad_box_dbl(:,:,:,dim,1) = dwork_grad_box_dbl(:,:,:,dim,1) * pot_dbl1
       dwork_grad_box_dbl(:,:,:,dim,2) = dwork_grad_box_dbl(:,:,:,dim,2) * pot_dbl2

       call fourier_filter(coarse_work,fine_work,&
            dwork_grad_box_dbl(:,:,:,dim,1),dwork_grad_box_dbl(:,:,:,dim,2),&
            dwork_grad_box(:,:,:,dim,1), dwork_grad_box(:,:,:,dim,2))
    end do
    ! ...have Vloc \nabla f (x, y, z components) in dwork_grad_box(:,:,:,1:3,1:2)

    ! 3. Apply div operator to (Vloc \nabla f)
    ! JCW: Copy real space FFTboxes containing (Vloc \nabla f)_{x,y,z} to complex
    ! JCW: array and Fourier transform to reciprocal space
    do dim = 1, 3
       call fourier_apply_box_pair('Coarse','Forward',dwork_grad_box(:,:,:,dim,1),&
            dwork_grad_box(:,:,:,dim,2),zwork_grad_box(:,:,:,dim,1),&
            omp=pub_threads_fftbox)
    end do

    ! JCW: Apply {x,y,z} components of grad operator to corresponding Cartesian
    ! JCW: components of (Vloc \nabla f)_{x,y,z} to give \nabla_{i} (Vloc \nabla f)_{i},
    ! JCW: where i = x, y, z.
    call kinetic_grad_on_cart_box(zwork_grad_box(:,:,:,:,1),&
         zwork_grad_box(:,:,:,:,2),fftbox)

    do dim=1,3
       ! JCW: Reverse transform FFTboxes containing \nabla_{i} (Vloc \nabla f)_{i} f
       ! JCW: to real space following application of gradient operator
       call fourier_apply_box_pair('Coarse','Backward',dwork_grad_box(:,:,:,dim,1),&
            dwork_grad_box(:,:,:,dim,2),zwork_grad_box(:,:,:,dim,2), &
            omp=pub_threads_fftbox)
    end do

    ! ...have \nabla . (Vloc \nabla f)

    ! 4. Accumulate result in output arrays
    data_box1 = 0.0_DP; data_box2 = 0.0_DP ! zero before accumulating
    do dim=1,3
       ! JCW: Sum components of dwork_grad_box(:,:,:,:,{1,2}) into data_box{1,2}
       ! JCW: to give final Vtau = -0.5 * \nabla . (Vloc \nabla f )
       data_box1(:,:,:) = data_box1(:,:,:) + dwork_grad_box(:,:,:,dim,1)
       data_box2(:,:,:) = data_box2(:,:,:) + dwork_grad_box(:,:,:,dim,2)
    end do

    ! jd: Silence var unused warning when in non-OMP compiles
    call utils_use_var(pub_threads_per_fftbox)

  end subroutine ngwf_gradient_div_locpot_grad

  subroutine test_ngwf_gradient_div_locpot_grad(&
             data_box1, data_box2,&
             test_data_box1, test_data_box2,&
             pot_dbl1, pot_dbl2,coarse_work,fine_work,&
             dwork_grad_box,dwork_grad_box_dbl,zwork_grad_box,&
             fftbox)
    !=========================================================!
    ! Testing routine for ngwf_gradient_div_locpot_grad.      !
    !                                                         !
    ! Checks that the output of ngwf_gradient_div_locpot_grad !
    ! is consistent with expectations.                        !
    !                                                         !
    ! Test 1: Comparison against ngwf_gradient_kinetic, by    !
    !         setting pot_dbl1 = pot_dbl2 = 1.0_DP            !
    !                                                         !
    ! Note, currently pot_dbl1 and pot_dbl2 are unused but    !
    ! may be used in future tests. test_data_box1, and        !
    ! test_data_box2 are also unused at present.              !
    !                                                         !
    ! Test results are output to stdout.                      !
    !---------------------------------------------------------!
    ! TODO Describe arguments                                 !
    !---------------------------------------------------------!
    ! James C. Womack, 2016                                   !
    !=========================================================!

    use constants, only: stdout, DP, file_maxsize
    use comms, only: pub_on_root, pub_total_num_procs
    use datatypes, only: FFTBOX_DATA, &
         data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box_pair, fourier_interpolate, fourier_filter
    use kinetic, only: kinetic_grad_on_box, kinetic_grad_on_cart_box, &
         kinetic_apply_on_box
    use rundat, only: pub_rootname
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_use_var, &
         utils_alloc_check, utils_dealloc_check, &
         utils_unit, utils_open_unit_check, utils_close_unit_check, &
         utils_use_var
    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(in) :: fftbox
    !real(kind=DP), intent(inout) :: data_box1(fftbox%total_ld1, &
    !     fftbox%total_ld2,fftbox%total_pt3)
    !real(kind=DP), intent(inout) :: data_box2(fftbox%total_ld1, &
    !     fftbox%total_ld2,fftbox%total_pt3)
    !real(kind=DP), intent(in) :: test_data_box1(fftbox%total_ld1, &
    !     fftbox%total_ld2,fftbox%total_pt3)
    !real(kind=DP), intent(in) :: test_data_box2(fftbox%total_ld1, &
    !     fftbox%total_ld2,fftbox%total_pt3)
    type(FFTBOX_DATA) :: data_box1
    type(FFTBOX_DATA) :: data_box2
    type(FFTBOX_DATA) :: test_data_box1
    type(FFTBOX_DATA) :: test_data_box2
    real(kind=DP), intent(in)    :: pot_dbl1(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl)
    complex(kind=DP), intent(inout) :: fine_work(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl)
    complex(kind=DP), intent(inout) :: coarse_work(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    real(kind=DP), intent(in)    :: pot_dbl2(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl)
    real(kind=DP), intent(inout) :: dwork_grad_box(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3,3,2)
    real(kind=DP), intent(inout) :: dwork_grad_box_dbl(fftbox%total_ld1_dbl,&
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl,3,2)
    complex(kind=DP), intent(inout) :: zwork_grad_box(fftbox%total_ld1,&
         fftbox%total_ld2, fftbox%total_pt3,3,2)

    ! Local variables
    !real(kind=DP),allocatable :: tmp_data_box1(:,:,:,:)
    !real(kind=DP),allocatable :: tmp_data_box2(:,:,:,:)
    ! JCW: Use FFTBOX_DATA type for compatibility with
    ! JCW: ngwf_gradient_kinetic, but always assume real NGWFs
    type(FFTBOX_DATA) :: tmp_data_box1(2)
    type(FFTBOX_DATA) :: tmp_data_box2(2)
    real(kind=DP),allocatable :: tmp_pot_dbl1(:,:,:)
    real(kind=DP),allocatable :: tmp_pot_dbl2(:,:,:)
    integer :: ierr !, dim <-- unused
    integer :: i1, i2, i3
    real(kind=DP) :: abs_diff1, abs_diff2
    real(kind=DP) :: max_diff1
    real(kind=DP) :: max_diff2
    real(kind=DP) :: varA1, varB1
    real(kind=DP) :: varA2, varB2
    real(kind=DP) :: rmsd1, rmsd2
    real(kind=DP) :: proc_test_data(8)
    ! proc_test_data(1):  max_diff1 (between set A1 and B1)
    ! proc_test_data(2):  max_diff2 (between set A2 and B2)
    ! proc_test_data(3):  varA1 at max_diff1
    ! proc_test_data(4):  varB1 at max_diff1
    ! proc_test_data(5):  varA2 at max_diff2
    ! proc_test_data(6):  varB2 at max_diff2
    ! proc_test_data(7):  rmsd1
    ! proc_test_data(8):  rmsd2
    integer       :: npoints
    integer,save :: itest = 0
    logical,save :: file_created = .false.
    character(len=file_maxsize) :: output_file
    integer      :: iunit
    integer      :: nthreads
!$  integer, external :: omp_get_thread_num, omp_get_num_threads

    ! Assert a single MPI process and OpenMP thread, since multiples of
    ! these confuse the output.
    call utils_assert(pub_total_num_procs.eq.1,&
         "ERROR in ngwf_gradient_div_locpot_grad: &
         &This routine does not support more than one MPI process. &
         &Please run ONETEP with only a single MPI process.")
    nthreads = 1
!$  nthreads = omp_get_num_threads()
    call utils_assert(nthreads.eq.1,"ERROR in ngwf_gradient_div_locpot_grad: &
         &This routine does not support OpenMP threading. Please set &
         &OMP_NUM_THREADS=1 in your shell.")

    if (pub_on_root) then
       ! Should only ever be called in debug mode
       write(stdout,'(a)') "DEBUG: Entering test_ngwf_gradient_div_locpot_grad"
    end if

    ! JCW: Use unused potential input arguments to avoid compiler warnings
    call utils_use_var(pot_dbl1(1,1,1))
    call utils_use_var(pot_dbl2(1,1,1))

    ! Open unit for writing out test data (rather than using output file)
    if (pub_on_root) then
       iunit = utils_unit()
       output_file = trim(pub_rootname)//"_ngwf_gradient_test.out"
       ! Should only ever be called in debug mode
       write(stdout,'(a,a,a,i5)') "DEBUG: Writing test output to file ", &
            trim(output_file),", unit: ", iunit
       if (.not.file_created) then
          ! Start by overwriting any pre-existing file of the same name
          open(unit=iunit,file=output_file,action="write",status="replace",&
               iostat=ierr)
          file_created = .true.
          ! Write header to file
          write(iunit,*) "JCW test code: ngwf_gradient_div_locpot_grad"
          ! Write intro to test
          write(iunit,*)
          write(iunit,'(a,i5)') "Test sequence ",itest
          itest = itest + 1
       else
          ! File already exists, so simply append to this
          open(unit=iunit,file=output_file,action="write",status="old",&
               position="append",iostat=ierr)
          ! Write intro for new test
          write(iunit,*)
          write(iunit,'(a,i5)') "Test sequence ",itest
          itest = itest + 1
       end if
       call utils_open_unit_check('test_ngwf_gradient_div_locpot_grad',&
            output_file,ierr)
    else
      iunit = -100 ! Set to nonsensical value -- should not be used on non-root
                   ! proc
    end if

    ! Start timer
    call timer_clock('test_ngwf_gradient_div_locpot_grad',1)

    ! Allocate temporary internal data arrays
    !allocate(tmp_data_box1(fftbox%total_ld1,fftbox%total_ld2,&
    !     fftbox%total_pt3,2),stat=ierr)
    !call utils_alloc_check('test_ngwf_gradient_div_locpot_grad',&
    !     'tmp_data_box1',ierr)
    !allocate(tmp_data_box2(fftbox%total_ld1,fftbox%total_ld2,&
    !     fftbox%total_pt3,2),stat=ierr)
    !call utils_alloc_check('test_ngwf_gradient_div_locpot_grad',&
    !     'tmp_data_box2',ierr)
    do i1 = 1, 2
       ! This routine does not support complex NGWFs, so
       ! we can assume iscmplx=.false.
       call data_fftbox_alloc(tmp_data_box1(i1),&
            fftbox%total_ld1,fftbox%total_ld2,&
            fftbox%total_pt3,iscmplx=.false.)
       call data_fftbox_alloc(tmp_data_box2(i1),&
            fftbox%total_ld1,fftbox%total_ld2,&
            fftbox%total_pt3,iscmplx=.false.)
    end do
    allocate(tmp_pot_dbl1(fftbox%total_ld1_dbl,fftbox%total_ld2_dbl,&
         fftbox%total_pt3_dbl),stat=ierr)
    call utils_alloc_check('test_ngwf_gradient_div_locpot_grad',&
         'tmp_pot_dbl1',ierr)
    allocate(tmp_pot_dbl2(fftbox%total_ld1_dbl,fftbox%total_ld2_dbl,&
         fftbox%total_pt3_dbl),stat=ierr)
    call utils_alloc_check('test_ngwf_gradient_div_locpot_grad',&
         'tmp_pot_dbl2',ierr)

    ! Copy input data arrays to temporary internal arrays
    !tmp_data_box1(:,:,:,1) = data_box1
    !tmp_data_box1(:,:,:,2) = data_box1
    !tmp_data_box2(:,:,:,1) = data_box2
    !tmp_data_box2(:,:,:,2) = data_box2
    tmp_data_box1(1)%d(:,:,:) = data_box1%d
    tmp_data_box1(2)%d(:,:,:) = data_box1%d
    tmp_data_box2(1)%d(:,:,:) = data_box2%d
    tmp_data_box2(2)%d(:,:,:) = data_box2%d

    ! Compare with kinetic energy operator, by setting
    !    Vloc = 1.0_DP.
    ! [ ngwf_gradient_div_locpot_grad ]
    tmp_pot_dbl1 = 1.0_DP
    tmp_pot_dbl2 = 1.0_DP
    !call ngwf_gradient_div_locpot_grad( &
    !     tmp_data_box1(:,:,:,1), &                   ! input-output
    !     tmp_data_box2(:,:,:,1), &                   ! input-output
    !     tmp_pot_dbl1, tmp_pot_dbl2, &               ! input-output
    !     coarse_work,fine_work, &                    ! workspace
    !     dwork_grad_box,dwork_grad_box_dbl,&         ! workspace
    !     zwork_grad_box, &                           ! workspace
    !     fftbox)
    call ngwf_gradient_div_locpot_grad( &
         tmp_data_box1(1)%d, &                   ! input-output
         tmp_data_box2(1)%d, &                   ! input-output
         tmp_pot_dbl1, tmp_pot_dbl2, &               ! input-output
         coarse_work,fine_work, &                    ! workspace
         dwork_grad_box,dwork_grad_box_dbl,&         ! workspace
         zwork_grad_box, &                           ! workspace
         fftbox)

    ! Multiply by -0.5_DP to obtain kinetic energy operator
    !tmp_data_box1(:,:,:,1) = -0.5_DP * tmp_data_box1(:,:,:,1)
    !tmp_data_box2(:,:,:,1) = -0.5_DP * tmp_data_box2(:,:,:,1)
    tmp_data_box1(1)%d(:,:,:) = -0.5_DP * tmp_data_box1(1)%d
    tmp_data_box2(1)%d(:,:,:) = -0.5_DP * tmp_data_box2(1)%d

    ! [ kinetic energy operator ]
    !call ngwf_gradient_kinetic(&
    !     tmp_data_box1(:,:,:,2), tmp_data_box2(:,:,:,2), &
    !     data_box1, data_box2, &
    !     coarse_work,fftbox)
    call ngwf_gradient_kinetic(&
         tmp_data_box1(2), tmp_data_box2(2), &
         data_box1, data_box2, &
         coarse_work,fftbox)


    proc_test_data = 0.0_DP
    ! [ Comparison ]

    ! max_difference
    max_diff1 = 0.0_DP
    max_diff2 = 0.0_DP
    do i3 = 1, fftbox%total_pt3
       do i2 = 1, fftbox%total_ld2
          do i1 = 1, fftbox%total_ld1
             !abs_diff1 = abs( tmp_data_box1(i1,i2,i3,1) - tmp_data_box1(i1,i2,i3,2) )
             !abs_diff2 = abs( tmp_data_box2(i1,i2,i3,1) - tmp_data_box2(i1,i2,i3,2) )
             abs_diff1 = abs( tmp_data_box1(1)%d(i1,i2,i3) - tmp_data_box1(2)%d(i1,i2,i3) )
             abs_diff2 = abs( tmp_data_box2(1)%d(i1,i2,i3) - tmp_data_box2(2)%d(i1,i2,i3) )
             if (abs_diff1.gt.max_diff1) then
               max_diff1 = abs_diff1
               !varA1 = tmp_data_box1(i1,i2,i3,1)
               !varB1 = tmp_data_box1(i1,i2,i3,2)
               varA1 = tmp_data_box1(1)%d(i1,i2,i3)
               varB1 = tmp_data_box1(2)%d(i1,i2,i3)
             end if
             if (abs_diff2.gt.max_diff2) then
               max_diff2 = abs_diff2
               !varA2 = tmp_data_box2(i1,i2,i3,1)
               !varB2 = tmp_data_box2(i1,i2,i3,2)
               varA2 = tmp_data_box2(1)%d(i1,i2,i3)
               varB2 = tmp_data_box2(2)%d(i1,i2,i3)
             end if
          end do
       end do
    end do

    ! RMSD
    rmsd1 = 0.0_DP
    rmsd2 = 0.0_DP
    do i3 = 1, fftbox%total_pt3
       do i2 = 1, fftbox%total_ld2
          do i1 = 1, fftbox%total_ld1
             ! Sum all mean squared differences
             !rmsd1 = rmsd1 + ( tmp_data_box1(i1,i2,i3,1) &
             !                - tmp_data_box1(i1,i2,i3,2) )**2
             !rmsd2 = rmsd2 + ( tmp_data_box2(i1,i2,i3,1) &
             !                - tmp_data_box2(i1,i2,i3,2) )**2
             rmsd1 = rmsd1 + ( tmp_data_box1(1)%d(i1,i2,i3) &
                             - tmp_data_box1(2)%d(i1,i2,i3) )**2
             rmsd2 = rmsd2 + ( tmp_data_box2(1)%d(i1,i2,i3) &
                             - tmp_data_box2(2)%d(i1,i2,i3) )**2
          end do
       end do
    end do
    ! Calculate RMSD
    npoints = fftbox%total_ld1*fftbox%total_ld2*fftbox%total_pt3
    rmsd1 = sqrt(rmsd1/real(npoints,kind=DP))
    rmsd2 = sqrt(rmsd2/real(npoints,kind=DP))

    ! Initialize per-proc data array:
    ! proc_test_data(1):  max_diff1 (between set A1 and B1)
    ! proc_test_data(2):  max_diff2 (between set A2 and B2)
    ! proc_test_data(3):  varA1 at max_diff1
    ! proc_test_data(4):  varB1 at max_diff1
    ! proc_test_data(5):  varA2 at max_diff2
    ! proc_test_data(6):  varB2 at max_diff2
    ! proc_test_data(7):  rmsd1
    ! proc_test_data(8):  rmsd2
    proc_test_data(1) = max_diff1
    proc_test_data(2) = max_diff2
    proc_test_data(3) = varA1
    proc_test_data(4) = varB1
    proc_test_data(5) = varA2
    proc_test_data(6) = varB2
    proc_test_data(7) = rmsd1
    proc_test_data(8) = rmsd2

    call output_test_data(iunit,1,&
         "Vloc = 1.0, comparison against ngwf_gradient_kinetic",proc_test_data)

    ! Deallocate temporary data arrays
    !deallocate(tmp_data_box1,stat=ierr)
    !call utils_dealloc_check('test_ngwf_gradient_div_locpot_grad',&
    !     'tmp_data_box1',ierr)
    !deallocate(tmp_data_box2,stat=ierr)
    !call utils_dealloc_check('test_ngwf_gradient_div_locpot_grad',&
    !     'tmp_data_box2',ierr)
    do i1 = 1, 2
       call data_fftbox_dealloc(tmp_data_box1(i1))
       call data_fftbox_dealloc(tmp_data_box2(i1))
    end do
    deallocate(tmp_pot_dbl1,stat=ierr)
    call utils_dealloc_check('test_ngwf_gradient_div_locpot_grad',&
         'tmp_pot_dbl1',ierr)
    deallocate(tmp_pot_dbl2,stat=ierr)
    call utils_dealloc_check('test_ngwf_gradient_div_locpot_grad',&
         'tmp_pot_dbl2',ierr)

    ! Close
    if (pub_on_root) then
       ! Should only ever be called in debug mode
       write(stdout,'(a,a,a,i5)') "DEBUG: Closing file ", &
            trim(output_file),", unit: ", iunit
       close(iunit,iostat=ierr)
       call utils_close_unit_check('test_ngwf_gradient_div_locpot_grad',&
            output_file,ierr)
    end if

    ! Stop timer
    call timer_clock('test_ngwf_gradient_div_locpot_grad',2)

    if (pub_on_root) then
       ! Should only ever be called in debug mode
       write(stdout,'(a)') "DEBUG: Leaving test_ngwf_gradient_div_locpot_grad"
    end if

    contains
      subroutine output_test_data(iunit,test_number,test_desc,proc_test_data)
        !=========================================================!
        ! Output per-proc results of test contained in array      !
        ! proc_test_data to unit iunit in human-readable format.  !
        ! See below for details of contents of proc_test_data.    !
        !=========================================================!
        ! James C. Womack, 2016                                   !
        !=========================================================!

        use comms, only: pub_my_proc_id, comms_barrier
        implicit none

        ! Arguments
        integer, intent(in)      :: iunit
        integer, intent(in)      :: test_number
        character(len=*),intent(in) :: test_desc
        real(kind=DP),intent(inout) :: proc_test_data(8)
        ! proc_test_data(1):  max_diff1 (between set A1 and B1)
        ! proc_test_data(2):  max_diff2 (between set A2 and B2)
        ! proc_test_data(3):  varA1 at max_diff1
        ! proc_test_data(4):  varB1 at max_diff1
        ! proc_test_data(5):  varA2 at max_diff2
        ! proc_test_data(6):  varB2 at max_diff2
        ! proc_test_data(7):  rmsd1
        ! proc_test_data(8):  rmsd2

        ! Local variables
        real(kind=DP) :: max_diff1
        real(kind=DP) :: max_diff2
        real(kind=DP) :: varA1, varB1
        real(kind=DP) :: varA2, varB2
        real(kind=DP) :: rmsd1, rmsd2
        integer :: iproc

        iproc = pub_my_proc_id
        ! Unpack data into human-readable variables
        ! proc_test_data(1):  max_diff1 (between set A1 and B1)
        ! proc_test_data(2):  max_diff2 (between set A2 and B2)
        ! proc_test_data(3):  varA1 at max_diff1
        ! proc_test_data(4):  varB1 at max_diff1
        ! proc_test_data(5):  varA2 at max_diff2
        ! proc_test_data(6):  varB2 at max_diff2
        ! proc_test_data(7):  rmsd1
        ! proc_test_data(8):  rmsd2
        max_diff1 = proc_test_data(1)
        max_diff2 = proc_test_data(2)
        varA1     = proc_test_data(3)
        varB1     = proc_test_data(4)
        varA2     = proc_test_data(5)
        varB2     = proc_test_data(6)
        rmsd1     = proc_test_data(7)
        rmsd2     = proc_test_data(8)
        ! Output test data
        write(iunit,'(a,i3)') "Test: ",test_number
        write(iunit,'(a)') test_desc
        write(iunit,'(a,i4,a)') "data_box1 (proc ",iproc,"):"
        write(iunit,'(a,es24.16)') "Maximum per-element difference = ", &
             max_diff1
        write(iunit,'(a,es24.16,a)') "(div_locpot_grad value         = ", &
             varA1,")"
        write(iunit,'(a,es24.16,a)') "(comparison value              = ", &
             varB1,")"
        write(iunit,'(a,es24.16)') "Root mean-square deviation     = ", &
             rmsd1
        write(iunit,'(a,i4,a)') "data_box2 (proc ",iproc,"):"
        write(iunit,'(a,es24.16)') "Maximum per-element difference = ", &
             max_diff2
        write(iunit,'(a,es24.16,a)') "(div_locpot_grad value         = ", &
             varA2,")"
        write(iunit,'(a,es24.16,a)') "(comparison value              = ", &
             varB2,")"
        write(iunit,'(a,es24.16)') "Root mean-square deviation     = ", &
             rmsd2

     end subroutine output_test_data

  end subroutine test_ngwf_gradient_div_locpot_grad


  subroutine ngwf_grad_smooth_rspace( &
       fftbox_grad, fa_start1, fa_start2, fa_start3, &
       fa_tight_box, sphere_centre, sphere_radius, fftbox, cell)

    !======================================================================!
    ! Smooth a function inside an FFTbox by convolution in real space      !
    ! with a localized low-pass filter                                     !
    !----------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on January 2011                        !
    !     largely based on ngwf_grad_precond_rspace originally written by  !
    !     Chris-Kriton Skylaris on 29/4/2003.                              !
    !======================================================================!

    use basis, only: FUNCTION_TIGHT_BOX, basis_copy_tightbox_to_fftbox, &
         basis_location_func_wrt_cell
    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use datatypes, only: FFTBOX_DATA, data_set_to_zero, &
         data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use geometry, only: POINT, operator(*), operator(.DOT.), operator(+)
    use rundat, only: pub_smooth_scheme, pub_r_smooth
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(CELL_INFO), intent(in) :: cell
    integer, intent(in) :: fa_start1, fa_start2, fa_start3
    type(FUNCTION_TIGHT_BOX), intent(in) :: fa_tight_box
    type(POINT), intent(in)   :: sphere_centre
    real(kind=DP), intent(in) :: sphere_radius
    type(FFTBOX_DATA), intent(inout) :: fftbox_grad

    ! Local Variables
    integer :: fa_tight_pts1, fa_tight_pts2, fa_tight_pts3
    integer :: row1,row2,row3,step1,step2,step3
    integer :: npt1,npt2,npt3
    integer :: row1_start, row2_start, row3_start
    integer :: period_pt1, period_pt2, period_pt3
    integer :: row1_end, row2_end, row3_end
    integer :: tc1, tc2, tc3
    type(FFTBOX_DATA) :: smooth_tightbox_grad
    integer :: ierr
    ! agrecocmplx
    logical :: loc_cmplx

    call timer_clock("ngwf_grad_smooth_rspace",1)

    ! agrecocmplx
    loc_cmplx = fftbox_grad%iscmplx

    ! agrecocmplx: should be ok now, need to check!
    !jmecmplx
    !call utils_assert(.not. loc_cmplx, 'Error in&
    !     & ngwf_grad_smooth_rspace: not yet ready for complex NGWFs.')

    !==========================================================================!
    !===== Preliminaries ======================================================!

    fa_tight_pts1 = fa_tight_box%tight_pts1
    fa_tight_pts2 = fa_tight_box%tight_pts2
    fa_tight_pts3 = fa_tight_box%tight_pts3

    ! smmd: determine number of grid points in smoothing radius
    npt1 = nint( pub_r_smooth/cell%d1 )
    npt2 = nint( pub_r_smooth/cell%d2 )
    npt3 = nint( pub_r_smooth/cell%d3 )

    ! smmd : debug
    if (pub_on_root) then
       write(stdout,*) "Smoothing npts : ", npt1, npt2, npt3
    endif

    ! smmd: allocate tightbox for smooth gradient
    call data_fftbox_alloc(smooth_tightbox_grad, fa_tight_pts1, &
         fa_tight_pts2, fa_tight_pts3, iscmplx=loc_cmplx)
    call data_set_to_zero(smooth_tightbox_grad)

    ! smmd: Initialise smoothing function in real space
    if ( .not. allocated(smooth_func_real) ) then
       call internal_init_smooth()
    endif

    ! smmd: Shave gradient in FFTbox
    call internal_shave_in_fftbox(sphere_radius-pub_r_smooth)

    ! smmd : debug
    if (pub_on_root) then
       write(stdout,*) "Smoothing shave : done ! "
    endif

    !==========================================================================!
    !===== Real space convolution ( gradient & smoothing function ) ===========!

    ! smmd: case without periodic boundary conditions on the fftbox
    if (.not.(fftbox%coin3 .or. &
        fftbox%coin2 .or. fftbox%coin1)) then

       ! cks: loop over all points selected for convolution
       do step3 = -npt3, npt3
          do step2 = -npt2, npt2
             do step1 = -npt1, npt1

                ! cks: convolve only in a spherical region around central point
                ! cks: and only for the points of the (tight) fft_box of fa
                ! cks: which is much smaller than the fft box
                if ( smooth_func_real(step1,step2,step3) .ne. 0.0_dp ) then

                   row3_start =max(1 -step3, fa_start3)
                   row2_start =max(1 -step2, fa_start2)
                   row1_start =max(1 -step1, fa_start1)

                   row3_end = min(fftbox%total_pt3 - step3, &
                        fa_start3 + fa_tight_pts3 -1)
                   row2_end = min(fftbox%total_pt2 - step2, &
                        fa_start2 + fa_tight_pts2 -1)
                   row1_end = min(fftbox%total_pt1 - step1, &
                        fa_start1 + fa_tight_pts1 -1)

                   tc3 = row3_start - fa_start3                      ! smmd
                   do row3= row3_start, row3_end
                      tc3 =tc3+1

                      tc2 = row2_start - fa_start2                   ! smmd
                      do row2= row2_start, row2_end
                         tc2 =tc2 +1

                         tc1 = row1_start - fa_start1                ! smmd
                         do row1= row1_start, row1_end
                            tc1 =tc1 +1

                            ! agrecocmplx: smooth_func_real should stay real
                            ! in any case, but need to check this
                            if (loc_cmplx) then
                               smooth_tightbox_grad%z(tc1, tc2, tc3) = &
                                    smooth_tightbox_grad%z(tc1, tc2, tc3) &
                                    + fftbox_grad%z(row1+step1,row2+step2,row3+step3) &
                                    * smooth_func_real(step1,step2,step3)

                            else
                               smooth_tightbox_grad%d(tc1, tc2, tc3) = &
                                    smooth_tightbox_grad%d(tc1, tc2, tc3) &
                                    + fftbox_grad%d(row1+step1,row2+step2,row3+step3) &
                                    * smooth_func_real(step1,step2,step3)

                            end if

                         enddo
                      enddo
                   enddo

                endif

             enddo
          enddo
       enddo

    ! smmd: case with periodic boundary conditions on the fftbox
    else

       ! cks: loop over all points selected for convolution
       do step3 = -npt3, npt3
          do step2 = -npt2, npt2
             do step1 = -npt1, npt1

                ! cks: convolve only in a spherical region around central point
                ! cks: and only for the points of the tight fft_box of fa
                ! cks: which is much smaller than the fft box
                if ( smooth_func_real(step1,step2,step3) .ne. 0.0_dp ) then

                   if (fftbox%coin3) then
                      row3_start = fa_start3
                      row3_end = fa_start3 + fa_tight_pts3 -1
                   else
                      row3_start =max(1 -step3, fa_start3)
                      row3_end = min(fftbox%total_pt3 - step3, &
                           fa_start3 + fa_tight_pts3 -1)
                   endif

                   if (fftbox%coin2) then
                      row2_start = fa_start2
                      row2_end = fa_start2 + fa_tight_pts2 -1
                   else
                      row2_start =max(1 -step2, fa_start2)
                      row2_end = min(fftbox%total_pt2 - step2, &
                           fa_start2 + fa_tight_pts2 -1)
                   endif

                   if (fftbox%coin1) then
                      row1_start = fa_start1
                      row1_end = fa_start1 + fa_tight_pts1 -1
                   else
                      row1_start =max(1 -step1, fa_start1)
                      row1_end = min(fftbox%total_pt1 - step1, &
                           fa_start1 + fa_tight_pts1 -1)
                   endif

                   tc3 = 0
                   do row3 = row3_start, row3_end
                      tc3 = tc3 + 1
                      if (fftbox%coin3) then
                         period_pt3 = modulo(row3+step3-1,fftbox%total_pt3)+1
                      else
                         period_pt3 = row3 + step3
                      endif

                      tc2 = 0
                      do row2 = row2_start, row2_end
                         tc2 = tc2 + 1
                         if (fftbox%coin2) then
                            period_pt2 = modulo(row2+step2-1,fftbox%total_pt2)+1
                         else
                            period_pt2 = row2 + step2
                         endif

                         tc1 = 0
                         do row1 = row1_start, row1_end
                            tc1 = tc1 + 1
                            if (fftbox%coin1) then
                               period_pt1 = modulo(row1+step1-1,fftbox%total_pt1)+1
                            else
                               period_pt1 = row1 + step1
                            endif

                            ! agrecocmplx: smooth_func_real should stay real
                            ! in any case, but need to check this
                            if (loc_cmplx) then
                               smooth_tightbox_grad%z(tc1, tc2, tc3) = &
                                    smooth_tightbox_grad%z(tc1, tc2, tc3) &
                                    + fftbox_grad%z(period_pt1,period_pt2,period_pt3) &
                                    * smooth_func_real(step1,step2,step3)

                            else
                               smooth_tightbox_grad%d(tc1, tc2, tc3) = &
                                    smooth_tightbox_grad%d(tc1, tc2, tc3) &
                                    + fftbox_grad%d(period_pt1,period_pt2,period_pt3) &
                                    * smooth_func_real(step1,step2,step3)

                            end if
                         enddo
                      enddo
                   enddo

                endif

             enddo
          enddo
       enddo

    endif

    ! smmd : debug
    if (pub_on_root) then
       write(stdout,*) "Smoothing convolution : done ! "
    endif

    ! cks: copy the tightbox with the smooth gradient to
    ! cks: the appropriate position in the fftbox that contained the
    ! cks: original gradient.
    call basis_copy_tightbox_to_fftbox(fftbox_grad, &                ! output
         fa_start1, fa_start2, fa_start3, smooth_tightbox_grad, &    ! input
          fa_tight_pts1, fa_tight_pts2, fa_tight_pts3)   ! input


    ! cks: free up memory
    call data_fftbox_dealloc(smooth_tightbox_grad)

    ! smmd : debug
    if (pub_on_root) then
       write(stdout,*) "Smoothing complete : done ! "
    endif

    call timer_clock("ngwf_grad_smooth_rspace", 2)


    contains

    !=========================================================================!
    !=========================================================================!

    subroutine internal_shave_in_fftbox(rcut)

      use basis, only: basis_func_centre_wrt_fftbox
      use geometry, only: local_displacement, geometry_distance

      implicit none

      ! Argument
      real(kind=DP), intent(in) :: rcut

      ! Local Variables
      type(POINT) :: fa_centre
      type(POINT) :: current_point
      type(POINT) :: periodic_point
      real(kind=DP) :: distance
      real(kind=DP) :: pt1, pt2, pt3
      integer :: fftbox_start1, fftbox_start2, fftbox_start3
      integer :: period1, period2, period3
      integer :: ip1, ip2, ip3
      logical :: in_sphere

      ! agrecocmplx: should be ok now
      !call utils_assert(.not. loc_cmplx, 'Error in&
      !     & internal_shave_in_fftbox: not yet ready for complex NGWFs.')

      ! Centre of the current function
      call basis_location_func_wrt_cell(fftbox_start1, fftbox_start2, &
           fftbox_start3, fa_tight_box, cell)
      fa_centre = basis_func_centre_wrt_fftbox(sphere_centre, &
           fa_start1,fa_start2,fa_start3,fftbox_start1,fftbox_start2, &
           fftbox_start3,cell)

      ! Check for periodic boundary conditions
      period1 = 0
      period2 = 0
      period3 = 0
      if (fftbox%coin1) period1 = 1
      if (fftbox%coin2) period2 = 1
      if (fftbox%coin3) period3 = 1

      ! Loop over the points in tight box
      do tc3 = 1, fa_tight_pts3
         pt3 = real(fa_start3 + tc3 - 1, DP)*cell%d3

         do tc2 = 1, fa_tight_pts2
            pt2 = real(fa_start2 + tc2 - 1, DP)*cell%d2

            do tc1 = 1, fa_tight_pts1
               pt1 = real(fa_start1 + tc1 - 1, DP)*cell%d1

               in_sphere = .false.
               current_point = local_displacement(cell%a1_unit,&
                       cell%a2_unit,cell%a3_unit,pt1,pt2,pt3)

               ! Account for the periodic images of current_point in neighboring cells
               do ip3 = -period3, period3
                  do ip2 = -period2, period2
                     do ip1 = -period1, period1

                        periodic_point = current_point &
                              + real(ip1,DP)*cell%a1 &
                              + real(ip2,DP)*cell%a2 &
                              + real(ip3,DP)*cell%a3
                        distance = geometry_distance(periodic_point,fa_centre)
                        if (distance .lt. rcut) in_sphere = .true.

                     enddo
                  enddo
               enddo

               if (.not. in_sphere) then
                  ! agrecocmplx
                  if (loc_cmplx) then
                     fftbox_grad%z(fa_start1+tc1-1,fa_start2+tc2-1,fa_start3+tc3-1) = (0.0_DP,0.0_DP)
                  else
                     fftbox_grad%d(fa_start1+tc1-1,fa_start2+tc2-1,fa_start3+tc3-1) = 0.0_DP
                  end if
               end if
            enddo
         enddo
      enddo

    end subroutine internal_shave_in_fftbox

    !=========================================================================!
    !=========================================================================!

    subroutine internal_init_smooth()

      use constants, only: two_pi

      implicit none

      ! Local Variables
      type(POINT) :: distance_vec
      complex(kind=DP), allocatable :: smooth_complex(:,:,:)
      real(kind=DP) :: rr
      integer :: s1, s2, s3

      allocate(smooth_func_real(-npt1:npt1,-npt2:npt2,-npt3:npt3),stat=ierr)
      call utils_alloc_check('ngwf_grad_smooth_rspace','smooth_func_real',&
           ierr)

      if (pub_smooth_scheme == 'BG' .or. &
          pub_smooth_scheme == 'MAURI' .or. &
          pub_smooth_scheme == 'TETER' ) then

         allocate(smooth_complex(fftbox%total_ld1, fftbox%total_ld2, &
              fftbox%total_pt3),stat=ierr)
         call utils_alloc_check('ngwf_grad_smooth_rspace','smooth_complex',&
              ierr)

         if (.not. allocated(smooth_func_recip)) &
              call ngwf_grad_init_smooth_recip(fftbox)

         smooth_complex = cmplx(smooth_func_recip,0.0_DP,kind=DP)

         ! aam: use new fourier routines
         call fourier_apply_box('Coarse','Backward',smooth_complex)

      endif

      do step3 = -npt3, npt3
         s3 = abs(step3)
         do step2 = -npt2, npt2
            s2 = abs(step2)
            do step1 = -npt1, npt1
               s1 = abs(step1)

               distance_vec = real(step1,kind=DP)*cell%d1*cell%a1_unit&
                    + real(step2,kind=DP)*cell%d2*cell%a2_unit &
                    + real(step3,kind=DP)*cell%d3*cell%a3_unit

               ! spherical convolution cut-off radius
               rr = sqrt( distance_vec.DOT.distance_vec  )

               if ( rr.le.pub_r_smooth ) then
                  if (pub_smooth_scheme .eq. 'GAUSSIAN') then
                     smooth_func_real(step1,step2,step3) &
                           = (1/(sqrt(two_pi)*(pub_r_smooth/6)))*exp(-rr**2/(2*(pub_r_smooth/6)**2))
                          != exp(-0.025)*bessel_in(rr,0.025)
                  else
                     smooth_func_real(step1,step2,step3) &
                          = real(smooth_complex(s1+1,s2+1,s3+1),kind=DP)
                  endif
               else
                  smooth_func_real(step1,step2,step3) = 0.0_DP
               endif

            enddo
         enddo
      enddo

      if (pub_smooth_scheme .ne. 'GAUSSIAN') then

         ! free unused memory
         deallocate(smooth_complex,stat=ierr)
         call utils_dealloc_check('ngwf_grad_smooth_rspace','smooth_complex',&
              ierr)
      endif

    end subroutine internal_init_smooth

    !=========================================================================!
    !=========================================================================!

  end subroutine ngwf_grad_smooth_rspace


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_grad_init_precond_real(fftbox,cell)

    !======================================================================!
    ! Initialises the kinetic energy preconditioner in real space.         !
    !----------------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris on 29/4/2003.            !
    ! Modified by Arash A. Mostofi on 26/6/2003, 10/7/2003 and 8/8/2003.   !
    ! Modified for speed by Chris-Kriton Skylaris on 1/3/2004.             !
    ! Modified for accuracy by Simon M.-M. Dubois on 15/04/2011            !
    ! Split to own subroutine by Nicholas Hine on 05/08/2011.              !
    !======================================================================!

    use constants, only: DP
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use geometry, only: POINT, operator(*), operator(.DOT.), operator(+)
    use rundat, only: pub_r_precond
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(CELL_INFO), intent(in) :: cell

    ! Local Variables
    integer :: step1,step2,step3
    integer :: s1,s2,s3,npt1,npt2,npt3
    integer :: ierr
    real(kind=DP) :: rr

    complex(kind=DP), dimension(:,:,:), allocatable :: prec_complex
    type(POINT) :: distance_vec

    ! determine number of grid points in preconditioning radius
    npt1 = nint( pub_r_precond/cell%d1 )
    npt2 = nint( pub_r_precond/cell%d2 )
    npt3 = nint( pub_r_precond/cell%d3 )

    ! ===============================================================
    ! ========= INIT REAL SPACE PRECONDITIONER ======================
    if ( .not. allocated(precond_func_real) ) then

       allocate(precond_func_real(-npt1:npt1,-npt2:npt2,-npt3:npt3),stat=ierr)
       call utils_alloc_check('ngwf_grad_init_precond_real','precond_func_real',&
            ierr)
       allocate(prec_complex(fftbox%total_ld1, fftbox%total_ld2, &
            fftbox%total_pt3),stat=ierr)
       call utils_alloc_check('ngwf_grad_init_precond_real','prec_complex',&
            ierr)

       prec_complex = cmplx(precond_func_recip,0.0_DP,kind=DP)

       ! aam: use new fourier routines
       call fourier_apply_box('Coarse','Backward',prec_complex)

       do step3 = -npt3, npt3
          s3 = abs(step3)
          do step2 = -npt2, npt2
             s2 = abs(step2)
             do step1 = -npt1, npt1
                s1 = abs(step1)

                distance_vec = real(step1,kind=DP)*cell%d1*cell%a1_unit&
                     + real(step2,kind=DP)*cell%d2*cell%a2_unit &
                     + real(step3,kind=DP)*cell%d3*cell%a3_unit

                ! spherical convolution cut-off radius
                rr = sqrt( distance_vec.DOT.distance_vec  )

                if ( rr.le.pub_r_precond ) then
                   precond_func_real(step1,step2,step3) &
                        = real(prec_complex(s1+1,s2+1,s3+1),kind=DP)
                else
                   precond_func_real(step1,step2,step3) = 0.0_DP
                endif

             enddo
          enddo
       enddo


       ! cks: free unused memory
       deallocate(prec_complex,stat=ierr)
       call utils_dealloc_check('ngwf_grad_init_precond_real','prec_complex',&
            ierr)

       ! smmd: no need to scale the preconditioning function
       ! smmd: the following code lines are therefore commented
       !precond_func_real(-npt1:npt1,-npt2:npt2,-npt3:npt3) = &
       !     1.0_DP / precond_func_real(0,0,0) &
       !     * precond_func_real(-npt1:npt1,-npt2:npt2,-npt3:npt3)

    endif
    ! ===== END INIT REAL SPACE PRECONDITIONER ======================
    ! ===============================================================

  end subroutine ngwf_grad_init_precond_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_grad_precond_rspace( &
       fftbox_grad, fa_start1, fa_start2, fa_start3, &
       fa_tight_pts1, fa_tight_pts2, fa_tight_pts3,fftbox,cell)

    !======================================================================!
    ! Performs kinetic energy preconditioning in real space on a function  !
    ! inside an FFTbox.                                                    !
    !----------------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris on 29/4/2003.            !
    ! Modified by Arash A. Mostofi on 26/6/2003, 10/7/2003 and 8/8/2003.   !
    ! Modified for speed by Chris-Kriton Skylaris on 1/3/2004.             !
    ! Modified for accuracy by Simon M.-M. Dubois on 15/04/2011            !
    ! Data parallel OpenMP added by Karl Wilkinson in July 2013            !
    !======================================================================!

    use basis, only: basis_copy_tightbox_to_fftbox
    use constants, only: DP
    use datatypes, only: FFTBOX_DATA, data_set_to_zero, &
         data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use rundat, only: pub_r_precond, pub_threads_per_fftbox
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_use_var, &
         utils_assert

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    integer, intent(in) :: fa_start1, fa_start2, fa_start3
    integer, intent(in) :: fa_tight_pts1
    integer, intent(in) :: fa_tight_pts2
    integer, intent(in) :: fa_tight_pts3
    type(FFTBOX_DATA), intent(inout) :: fftbox_grad

    ! Local Variables
    integer :: row1,row2,row3
    integer :: period_pt1,period_pt2,period_pt3
    integer :: step1,step2,step3
    integer :: npt1,npt2,npt3
    integer :: row1_start
    integer :: row2_start
    integer :: row3_start
    integer :: row1_end
    integer :: row2_end
    integer :: row3_end
    integer :: tc1
    integer :: tc2
    integer :: tc3
    type(FFTBOX_DATA) :: prec_tightbox_grad
    ! agrecocmplx
    logical :: loc_cmplx

    ! -------------------------------------------------------------------------

    call timer_clock("ngwf_grad_precond_rspace",1)

    ! agrecocmplx
    loc_cmplx = fftbox_grad%iscmplx

    ! agrecocmplx: should be ok now, need to make sure precond_func_real
    ! stays real even in the complex case
    !jmecmplx
    !call utils_assert(.not. loc_cmplx, 'Error in&
    !     & ngwf_grad_precond_rspace: not ready yet for complex NGWFs.')

    ! jd: Silence var unused warning when in non-OMP compiles
    call utils_use_var(pub_threads_per_fftbox)

    call data_fftbox_alloc(prec_tightbox_grad, fa_tight_pts1, fa_tight_pts2, &
         fa_tight_pts3, iscmplx=loc_cmplx)

    ! determine number of grid points in preconditioning radius
    npt1 = nint( pub_r_precond/cell%d1 )
    npt2 = nint( pub_r_precond/cell%d2 )
    npt3 = nint( pub_r_precond/cell%d3 )

    call data_set_to_zero(prec_tightbox_grad)

!!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) DEFAULT(none) &
!!$OMP PRIVATE(step3,step2,step1,row3_start,row2_start,row1_start, &
!!$OMP      row3_end,row2_end,row1_end,tc3,tc2,tc1,row3,row2,row1, &
!!$OMP      period_pt3,period_pt2,period_pt1) &
!!$OMP SHARED(npt3,npt2,npt1,fa_start3,fa_start2,fa_start1, &
!!$OMP      prec_tightbox_grad,fftbox,precond_func_real,fa_tight_pts3, &
!!$OMP      fa_tight_pts2,fa_tight_pts1,fftbox_grad,pub_threads_per_fftbox,
!!$OMP      loc_cmplx)

    ! smmd: case without periodic boundary conditions on the fftbox
    if (.not.(fftbox%coin3 .or. &
        fftbox%coin2 .or. fftbox%coin1)) then

       ! cks: loop over all points selected for convolution
!!$OMP DO
       do step3 = -npt3, npt3
          do step2 = -npt2, npt2
             do step1 = -npt1, npt1

                ! cks: convolve only in a spherical region around central point
                ! cks: and only for the points of the (tight) fft_box of fa
                ! cks: which is much smaller than the fft box
                ! cks: ( = the whole precond_grad array).
                if ( precond_func_real(step1,step2,step3) .ne. 0.0_dp ) then

                   row3_start =max(1 -step3, fa_start3)
                   row2_start =max(1 -step2, fa_start2)
                   row1_start =max(1 -step1, fa_start1)

                   row3_end = min(fftbox%total_pt3 - step3, &
                        fa_start3 + fa_tight_pts3 -1)
                   row2_end = min(fftbox%total_pt2 - step2, &
                        fa_start2 + fa_tight_pts2 -1)
                   row1_end = min(fftbox%total_pt1 - step1, &
                        fa_start1 + fa_tight_pts1 -1)

                   tc3 = row3_start - fa_start3                      ! smmd
                   do row3= row3_start, row3_end
                      tc3 =tc3+1

                      tc2 = row2_start - fa_start2                   ! smmd
                      do row2= row2_start, row2_end
                         tc2 =tc2 +1

                         tc1 = row1_start - fa_start1                ! smmd
                         do row1= row1_start, row1_end
                            tc1 =tc1 +1

                            ! agrecocmplx: precond_func_real should stay real in
                            ! any case, but need to check this
                            if (loc_cmplx) then
                               prec_tightbox_grad%z(tc1, tc2, tc3) = &
                                    prec_tightbox_grad%z(tc1, tc2, tc3) &
                                    + fftbox_grad%z(row1+step1,row2+step2,row3+step3) &
                                    * precond_func_real(step1,step2,step3)

                            else
                               prec_tightbox_grad%d(tc1, tc2, tc3) = &
                                    prec_tightbox_grad%d(tc1, tc2, tc3) &
                                    + fftbox_grad%d(row1+step1,row2+step2,row3+step3) &
                                    * precond_func_real(step1,step2,step3)

                            end if


                         enddo
                      enddo
                   enddo

                endif

             enddo
          enddo
       enddo
!!$OMP END DO

    ! smmd: apply the correct periodic boundary conditions on the fftbox
    else

       ! cks: loop over all points selected for convolution
!!$OMP DO
       do step3 = -npt3, npt3
          do step2 = -npt2, npt2
             do step1 = -npt1, npt1

                ! cks: convolve only in a spherical region around central point
                ! cks: and only for the points of the tight fft_box of fa
                ! cks: which is much smaller than the fft box
                if ( precond_func_real(step1,step2,step3) .ne. 0.0_dp ) then

                   if (fftbox%coin3) then
                      row3_start = fa_start3
                      row3_end = fa_start3 + fa_tight_pts3 -1
                   else
                      row3_start =max(1 -step3, fa_start3)
                      row3_end = min(fftbox%total_pt3 - step3, &
                           fa_start3 + fa_tight_pts3 -1)
                   endif

                   if (fftbox%coin2) then
                      row2_start = fa_start2
                      row2_end = fa_start2 + fa_tight_pts2 -1
                   else
                      row2_start =max(1 -step2, fa_start2)
                      row2_end = min(fftbox%total_pt2 - step2, &
                           fa_start2 + fa_tight_pts2 -1)
                   endif

                   if (fftbox%coin1) then
                      row1_start = fa_start1
                      row1_end = fa_start1 + fa_tight_pts1 -1
                   else
                      row1_start =max(1 -step1, fa_start1)
                      row1_end = min(fftbox%total_pt1 - step1, &
                           fa_start1 + fa_tight_pts1 -1)
                   endif

                   tc3 = 0
                   do row3 = row3_start, row3_end
                      tc3 = tc3 + 1
                      if (fftbox%coin3) then
                         period_pt3 = modulo(row3+step3-1,fftbox%total_pt3)+1
                      else
                         period_pt3 = row3 + step3
                      endif

                      tc2 = 0
                      do row2 = row2_start, row2_end
                         tc2 = tc2 + 1
                         if (fftbox%coin2) then
                            period_pt2 = modulo(row2+step2-1,fftbox%total_pt2)+1
                         else
                            period_pt2 = row2 + step2
                         endif

                         tc1 = 0
                         do row1 = row1_start, row1_end
                            tc1 = tc1 + 1
                            if (fftbox%coin1) then
                               period_pt1 = modulo(row1+step1-1,fftbox%total_pt1)+1
                            else
                               period_pt1 = row1 + step1
                            endif

                            ! agrecocmplx: precond_func_real should stay real in
                            ! any case, but need to check
                            if (loc_cmplx) then
                               prec_tightbox_grad%z(tc1, tc2, tc3) = &
                                    prec_tightbox_grad%z(tc1, tc2, tc3) &
                                    + fftbox_grad%z(period_pt1,period_pt2,period_pt3) &
                                    * precond_func_real(step1,step2,step3)

                            else
                               prec_tightbox_grad%d(tc1, tc2, tc3) = &
                                    prec_tightbox_grad%d(tc1, tc2, tc3) &
                                    + fftbox_grad%d(period_pt1,period_pt2,period_pt3) &
                                    * precond_func_real(step1,step2,step3)

                            end if


                         enddo
                      enddo
                   enddo

                endif

             enddo
          enddo
       enddo
!!$OMP END DO
    endif
!!$OMP END PARALLEL


    ! cks: copy the tightbox with the preconditioned gradient to
    ! cks: the appropriate position in the fftbox that contained the
    ! cks: original gradient.
    call basis_copy_tightbox_to_fftbox(fftbox_grad, &                    !output
         fa_start1, fa_start2, fa_start3, &                               !input
         prec_tightbox_grad, fa_tight_pts1, fa_tight_pts2, fa_tight_pts3) !input


    ! cks: free up memory
    call data_fftbox_dealloc(prec_tightbox_grad)

    call timer_clock("ngwf_grad_precond_rspace", 2)


  end subroutine ngwf_grad_precond_rspace


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_gradient_precond_recip(data_box,fftbox)

    !======================================================================!
    ! Performs kinetic energy preconditioning in reciprocal space on a     !
    ! function inside an FFTbox.                                           !
    !----------------------------------------------------------------------!
    ! Modified on 08/08/2003 by aam.                                       !
    ! Modified on 28/12/2003 by cks.                                       !
    ! Modified on 15/7/2004 by pdh.                                        !
    !======================================================================!

    use fft_box, only: FFTBOX_INFO
    use rundat, only: pub_threads_per_fftbox
    use utils, only: utils_use_var

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(in) :: fftbox
    complex(kind=DP), intent(inout) :: data_box(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)

!!$OMP PARALLEL WORKSHARE NUM_THREADS(pub_threads_per_fftbox)
    data_box = data_box * precond_func_recip
!!$OMP END PARALLEL WORKSHARE

    ! jd: Silence var unused warning when in non-OMP compiles
    call utils_use_var(pub_threads_per_fftbox)

  end subroutine ngwf_gradient_precond_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_grad_init_precond_recip(fftbox)

    !======================================================================!
    ! Initialises the preconditioning function in reciprocal space.        !
    !----------------------------------------------------------------------!
    ! Originally written in 2003.                                          !
    !======================================================================!

    use constants, only: DP
    use fft_box, only: FFTBOX_INFO
    use rundat, only: pub_precond_scheme, pub_smooth_scheme, pub_k_zero
    use utils, only: utils_alloc_check, utils_assert, utils_abort

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local variables
    integer :: ierr
    integer :: i1,i2,i3,i1start
    real(kind=DP) :: scale
    real(kind=DP) :: x,temp

    ! cks: test to avoid divisions by zero and other problems
    call utils_assert(pub_k_zero > 0.0_DP, &
         'Error in ngwf_gradient_precond_recip: invalid k_zero =', pub_k_zero)

    allocate(precond_func_recip(fftbox%total_ld1,fftbox%total_ld2, &
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('ngwf_grad_init_precond_recip', &
         'precond_func_recip',ierr)

    scale = 1.0_DP / (0.5_DP * pub_k_zero * pub_k_zero)
    precond_func_recip = 0.0_DP

    select case (pub_precond_scheme)
    case ('BG')    ! B-G method
       do i3=1,fftbox%total_pt3
          do i2=1,fftbox%total_pt2
             do i1=1,fftbox%total_pt1
                x = scale * fftbox%recip_grid(5,i1,i2,i3)
                if (pub_smooth_scheme .eq. 'NONE' .or. pub_smooth_scheme .eq. 'GAUSSIAN') then
                   precond_func_recip(i1,i2,i3) = 1.0_DP / (1.0_DP + x)
                else
                   precond_func_recip(i1,i2,i3) = sqrt(1.0_DP / (1.0_DP + x))
                endif
             end do
          end do
       end do
    case ('MAURI') ! Mauri method
       precond_func_recip(1,1,1) = 1.0_DP   ! G=0
       i1start = 2                          ! miss G=0
       do i3=1,fftbox%total_pt3
          do i2=1,fftbox%total_pt2
             do i1=i1start,fftbox%total_pt1
                x = scale * fftbox%recip_grid(5,i1,i2,i3)
                if (pub_smooth_scheme .eq. 'NONE' .or. pub_smooth_scheme .eq. 'GAUSSIAN') then
                   precond_func_recip(i1,i2,i3) = min(1.0_DP,1.0_DP/x)
                else
                   precond_func_recip(i1,i2,i3) = sqrt(min(1.0_DP,1.0_DP/x))
                endif
             end do
             i1start = 1
          end do
       end do
    case ('TETER') ! Teter-Allen-Payne method
       do i3=1,fftbox%total_pt3
          do i2=1,fftbox%total_pt2
             do i1=1,fftbox%total_pt1
                x = scale * fftbox%recip_grid(5,i1,i2,i3)
                temp = 27.0_DP + x * (18.0_DP + x * (12.0_DP + 8.0_DP * x))
                if (pub_smooth_scheme .eq. 'NONE' .or. pub_smooth_scheme .eq. 'GAUSSIAN') then
                   precond_func_recip(i1,i2,i3) = temp / (temp + 16.0_DP * x**4)
                else
                   precond_func_recip(i1,i2,i3) = sqrt(temp / (temp + 16.0_DP * x**4))
                endif
             end do
          end do
       end do
    case ('NONE')
       precond_func_recip = 1.0_DP
    case default
       call utils_abort('Error in ngwf_grad_init_precond_recip: preconditioning&
            & scheme "'//trim(pub_precond_scheme)//'" not recognised')
    end select

  end subroutine ngwf_grad_init_precond_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_grad_init_smooth_recip(fftbox)

    !======================================================================!
    ! Initialises the preconditioning function in reciprocal space.        !
    !----------------------------------------------------------------------!
    ! Originally written in 2003.                                          !
    !======================================================================!

    use constants, only: DP
    use fft_box, only: FFTBOX_INFO
    use rundat, only: pub_precond_scheme, pub_smooth_scheme, pub_k_smooth
    use utils, only: utils_alloc_check, utils_abort

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local variables
    integer :: ierr
    integer :: i1,i2,i3,i1start
    real(kind=DP) :: scale
    real(kind=DP) :: x,temp

    if (pub_smooth_scheme .ne. 'SHAVE') then
       allocate(smooth_func_recip(fftbox%total_ld1,fftbox%total_ld2, &
            fftbox%total_pt3),stat=ierr)
       call utils_alloc_check('ngwf_grad_init_smooth_recip', &
            'smooth_func_recip',ierr)

       scale = 1.0_DP / (0.5_DP * pub_k_smooth * pub_k_smooth)
       smooth_func_recip = 0.0_DP

       select case (pub_smooth_scheme)
       case ('BG')    ! B-G method
          do i3=1,fftbox%total_pt3
             do i2=1,fftbox%total_pt2
                do i1=1,fftbox%total_pt1
                   x = scale * fftbox%recip_grid(5,i1,i2,i3)
                   smooth_func_recip(i1,i2,i3) = sqrt(1.0_DP / (1.0_DP + x))
                end do
             end do
          end do
       case ('MAURI') ! Mauri method
          smooth_func_recip(1,1,1) = 1.0_DP   ! G=0
          i1start = 2                          ! miss G=0
          do i3=1,fftbox%total_pt3
             do i2=1,fftbox%total_pt2
                do i1=i1start,fftbox%total_pt1
                   x = scale * fftbox%recip_grid(5,i1,i2,i3)
                   smooth_func_recip(i1,i2,i3) = sqrt(min(1.0_DP,1.0_DP/x))
                end do
                i1start = 1
             end do
          end do
       case ('TETER') ! Teter-Allen-Payne method
          do i3=1,fftbox%total_pt3
             do i2=1,fftbox%total_pt2
                do i1=1,fftbox%total_pt1
                   x = scale * fftbox%recip_grid(5,i1,i2,i3)
                   temp = 27.0_DP + x * (18.0_DP + x * (12.0_DP + 8.0_DP * x))
                   smooth_func_recip(i1,i2,i3) = sqrt(temp / (temp + 16.0_DP * x**4))
                end do
             end do
          end do
       case ('NONE')
          smooth_func_recip = 1.0_DP
       case default
          call utils_abort('Error in ngwf_grad_init_smooth_recip: precondition&
               &ing scheme "'//trim(pub_precond_scheme)//'" not recognised')
       end select
    endif

  end subroutine ngwf_grad_init_smooth_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_gradient_exit

    !===========================================!
    ! Free up allocated memory from the module. !
    !-------------------------------------------!
    ! Written by Arash A. Mostofi on 10/7/2003  !
    !===========================================!

    use utils, only: utils_dealloc_check

    implicit none

    integer :: ierr

    if (allocated(precond_func_real)) then
       deallocate(precond_func_real,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_exit','precond_func_real',ierr)
    end if

    if (allocated(precond_func_recip)) then
       deallocate(precond_func_recip,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_exit','precond_func_recip',ierr)
    end if

  end subroutine ngwf_gradient_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_gradient_coeffs(qmat,tc_qmat,pur_denskern,tc_pur_denskern, &
       rq,tc_rq,rk,tc_rk,hub_wk,tc_hub_wk,scissor_mat,tc_scissor_mat, &
       cond_grad,tc_cond_grad, denskern,overlap,ham,inv_overlap,ps_overlap, &
       sp_overlap,hub,hub_overlap_t,mu,n_occ,dijhat,mdl, &
       is_force_call, val_dkn, val_overlap, val_sp_overlap, val_ham,&
       cross_overlap, cond_shift, rep, ngwf_basis)

    !========================================================================!
    ! This subroutine returns the various matrices which multiply the        !
    ! NGWFs in row-sums required for the NGWF gradient                       !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    !  qmat (inout) : Q matrix in SPAM3 storage.                             !
    !  tc_qmat (inout) : Tensor-corrected Q matrix, QS, in SPAM3 storage.    !
    !  pur_denskern (inout) : Purified density kernel K in SPAM3 storage     !
    !  tc_pur_denskern (inout) : Tensor-corrected purified density kernel    !
    !    KS in SPAM3 storage.                                                !
    !  rk (inout) : Coefficients of nonlocal projectors in gradient in SPAM3 !
    !    storage.                                                            !
    !  rsk (inout) : Tensor-corrected coefficients of nonlocal projectors in !
    !    gradient in SPAM3 storage.                                          !
    !  hub_wk (inout) : Coefficients of DFT+U projectors in gradient in      !
    !  SPAM3 storage.                                                        !
    !  tc_hub_wk (inout) : Coefficients of DFT+U projectors in gradient in   !
    !  SPAM3 storage.                                                        !
    !  ham (input) : Hamiltonian matrix in SPAM3 format.                     !
    !  denskern (input) : Density kernel in DKERN format.                    !
    !  overlap (input) : Overlap matrix in SPAM3 format.                     !
    !  inv_overlap (input) : Inverse overlap matrix in SPAM3 format.         !
    !  mu (input) : Lagrange multiplier (original version) or fudge          !
    !    parameter (Millam-Scuseria version) of LNV function.                !
    !  n_occ (input) : number of occupied orbitals for each spin channel.    !
    !  dijhat (input) : Screened part of nonlocal projector energies         !
    !  val_dkn (input) : Valence density kernel (optional)                   !
    !  val_overlap (input) : Valence overlap matrix (optional)               !
    !  val_sp_overlap (input) : Valence NGWF-proj overlap matrix (optional)  !
    !  val_ham (input) : Valence Hamiltonian matrix (optional)               !
    !  cross_overlap (input) : Valence-conduction overlap matrix (optional)  !
    !  cond_shift (input) : Value by which the projected conduction          !
    !    Hamiltonian is currently being shifted (optional)                   !
    !------------------------------------------------------------------------!
    ! Written by Nicholas D.M. Hine on 26/07/2009 reusing some code from     !
    !   previous versions of ngwf_gradient_lnv written by Chris-Kriton       !
    !   Skylaris in 2000-2007.                                               !
    ! DFT+U added by David D. O'Regan in September 2009.                     !
    ! Modified for conduction calculation by Laura Ratcliff in October 2010. !
    ! Modified by Alvaro Ruiz Serrano for kernel DIIS in November 2010.      !
    ! Modified for embedding by Robert Charlton, September 2017.             !
    !========================================================================!

    use augmentation, only: aug_projector_denskern
    use constants, only: DP, max_spins, stdout
    use ensemble_dft, only: edft_ngwf_gradient
    use fragment_scfmi, only: scfmi_ngwf_gradient
    use function_basis, only: FUNC_BASIS
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: DKERN, kernel_purify, kernel_rescale_spam3
    use kernel_diis, only: kernel_diis_build_pq
    use mermin, only: mermin_ngwf_gradient
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use paw, only: paw_nonlocal_energies, paw_projector_overlap
    use pseudopotentials, only: PSEUDO_SPECIES, pseudo_get_dij, &
         pseudo_aug_Q_matrix
    use rundat, only: pub_exact_lnv, pub_occ_mix, pub_hubbard, pub_any_nl_proj, &
         pub_paw, pub_usp, pub_aug, pub_cond_calculate, pub_kernel_diis, &
         pub_edft, pub_num_spins, pub_spin_fac, pub_num_kpoints, PUB_1K, &
         pub_eda_scfmi, pub_debug_on_root, pub_scissor_ngroups, pub_mermin, &
         pub_mermin_smearing_width
    use scissor, only: scissor_shift_denskern
    use sparse, only: SPAM3, sparse_destroy, &
         sparse_create, sparse_product, sparse_scale, sparse_copy, &
         sparse_axpy, sparse_transpose_structure, sparse_transpose
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_destroy, &
         sparse_embed_create, sparse_embed_product, sparse_embed_scale, &
         sparse_embed_copy, sparse_embed_axpy, sparse_embed_transpose_structure, &
         sparse_embed_transpose, sparse_embed_array_product, &
         sparse_embed_array_scale, sparse_embed_array_copy, SPAM3_EMBED_ARRAY, &
         sparse_embed_array_extract_sub, sparse_embed_array_create, &
         sparse_embed_extract_sub, sparse_embed_array_destroy, &
         sparse_embed_array_create_sub, sparse_embed_create_sub, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: qmat(pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: tc_qmat(pub_num_spins)
    type(SPAM3_EMBED_ARRAY), intent(inout) :: pur_denskern
    type(SPAM3_EMBED_ARRAY), intent(inout) :: tc_pur_denskern
    type(SPAM3_EMBED), intent(inout) :: rq(pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: tc_rq(pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: rk(pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: tc_rk(pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: hub_wk(pub_num_spins)  !ddor
    type(SPAM3_EMBED), intent(inout) :: tc_hub_wk(pub_num_spins) !ddor
    type(SPAM3_EMBED), intent(inout) :: scissor_mat(pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: tc_scissor_mat(pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: cond_grad(pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: tc_cond_grad(pub_num_spins)
    type(DKERN), intent(inout) :: denskern
    type(SPAM3_EMBED), intent(in) :: overlap
    type(SPAM3_EMBED), intent(in) :: ham(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: inv_overlap
    type(SPAM3_EMBED), intent(in) :: ps_overlap
    type(SPAM3_EMBED), intent(in) :: sp_overlap
    type(SPAM3_EMBED), intent(in) :: hub_overlap_t
    real(kind=DP), intent(in) :: mu(max_spins)
    integer, intent(in) :: n_occ(:,:)
    type(SPAM3_EMBED), intent(in) :: dijhat(:)
    type(HUBBARD_MODEL), intent(in) :: hub
    ! rc2013: need pseudo_sp info from here; maybe there's a better way?
    type(MODEL), intent(in) :: mdl
    !type(PSEUDO_SPECIES), intent(in) :: pseudo_sp(:)
    logical, intent(in) :: is_force_call
    ! lr408: Optional conduction input arguments
    type(SPAM3_EMBED_ARRAY), optional, intent(in) :: val_dkn
    type(SPAM3_EMBED), optional, intent(in) :: val_ham(pub_num_spins)
    type(SPAM3_EMBED), optional, intent(in) :: val_overlap
    type(SPAM3_EMBED), optional, intent(in) :: val_sp_overlap
    type(SPAM3_EMBED), optional, intent(in) :: cross_overlap
    real(kind=dp), optional, intent(in) :: cond_shift
    type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(rep%nsub)

    ! Local Variables
    type(SPAM3_EMBED) :: tc_buffer
    type(SPAM3_EMBED) :: oij
    type(SPAM3_EMBED), allocatable :: rhoij(:)
    type(SPAM3_EMBED), allocatable :: dij(:)
    type(SPAM3_EMBED), allocatable :: dij_ps_overlap(:)
    type(SPAM3_EMBED), allocatable :: rk_cond(:), tc_rk_cond(:)
    type(SPAM3_EMBED) :: oijr_scissor
    type(SPAM3_EMBED) :: hub_ps_overlap
    type(SPAM3_EMBED) :: val_ps_overlap, oij_val_ps_overlap
    type(SPAM3_EMBED) :: ktm, kt, khv, ksv
    type(SPAM3) :: shifted_denskern, sk
    real(kind=DP),dimension(pub_num_spins, pub_num_kpoints) :: norm_fac, t_occ
    integer :: is
    integer :: ierr
    ! agrecocmplx
    logical :: loc_cmplx
    type(SPAM3) :: hub_projector_ham_cmplx
    character(len=*), parameter :: myself = 'ngwf_gradient_coeffs'
    integer :: mrows, ncols, isub, jsub
    ! jcap: work around for spins
    type(SPAM3), allocatable :: dij_tmp(:),rhoij_tmp(:),kern_array(:)

    ! -------------------------------------------------------------------------

    call timer_clock(myself, 1)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering '//myself//'.'

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, myself//&
         ': not ready yet for more than one k-point.')

    ! agrecocmplx
    loc_cmplx = overlap%iscmplx

    ! rc2013: get subsystem info
    mrows = overlap%mrows
    ncols = overlap%ncols

    ! ndmh: Choose means of determining Q matrix according to the procedure
    ! ndmh: used to optimise the density kernel
    if ((.not.pub_kernel_diis).and.(.not.pub_edft).and.(.not.pub_mermin).and. &
       (.not.(pub_eda_scfmi))) then

       ! cks: Construct purified density kernel K from density kernel L
       call kernel_purify(pur_denskern, denskern, overlap, inv_overlap, n_occ)

       ! Form NGWF gradient matrix
       if (pub_exact_lnv) then

          t_occ = real(n_occ,kind=dp)
          call kernel_rescale_spam3(pur_denskern, overlap, t_occ, &
               norm_fac=norm_fac, silent=.true.)

          call ngwf_grad_matrix_lnv_robust(qmat, &
               pur_denskern%m(:,PUB_1K), &
               denskern%kern%m(:,PUB_1K), ham, overlap, mu, &
               norm_fac(:,PUB_1K))
       else
          call ngwf_grad_matrix_ms_robust(qmat, &
               denskern%kern%m(:,PUB_1K), ham, overlap, mu)
       end if

    else if (pub_kernel_diis) then

       ! Form NGWF gradient matrix for kernel DIIS
       call kernel_diis_build_pq(pur_denskern%m(:,PUB_1K), qmat, &
            denskern%kern%m(:,PUB_1K), ham, overlap, &
            n_occ(:,PUB_1K), mu)

    else if (pub_edft) then

       ! ars: Form NGWF gradient matrix for EDFT
       call edft_ngwf_gradient(pur_denskern%m(:,PUB_1K), qmat, &
            denskern%kern%m(:,PUB_1K), ham, inv_overlap)

    else if (pub_mermin) then

       !ep: Form mermin gradient matrix for MERMIN method
       call mermin_ngwf_gradient(pur_denskern%m(:,PUB_1K), qmat, &
            denskern, inv_overlap, overlap )
       !ep

    else if (pub_eda_scfmi) then

       ! mjsp: Form NGWF gradient matrix for SCF-MI
       call scfmi_ngwf_gradient(pur_denskern%m(:,PUB_1K), qmat, &
            denskern%kern%m(:,PUB_1K), ham, inv_overlap)

    end if

    ! lr408: ++++ CONDUCTION GRADIENT COEFFS (SPAM3 MATRICES) ++++++++++++++++++
    if (pub_cond_calculate .and. .not. is_force_call) then

       call sparse_embed_create(kt,val_dkn%m(1,1),cross_overlap)
       call sparse_embed_create(ktm, kt, pur_denskern%m(1,1))
       call sparse_embed_create(khv,val_dkn%m(1,1),val_ham(1))
       call sparse_embed_create(ksv,val_dkn%m(1,1),val_overlap)

       do is=1,pub_num_spins

          call sparse_embed_product(kt,val_dkn%m(is,PUB_1K),cross_overlap)
          call sparse_embed_product(ktm, kt, pur_denskern%m(is,PUB_1K))
          call sparse_embed_product(khv,val_dkn%m(is,PUB_1K),val_ham(is))
          call sparse_embed_product(ksv,val_dkn%m(is,PUB_1K),val_overlap)

          call sparse_embed_scale(khv,-1.0_dp)
          call sparse_embed_axpy(khv,ksv,cond_shift)

          ! lr408: -KHVKTM + wKSVKTM
          call sparse_embed_product(cond_grad(is),khv,ktm)

       end do

       call sparse_embed_destroy(khv)
       call sparse_embed_destroy(ksv)
       call sparse_embed_destroy(ktm)
       call sparse_embed_destroy(kt)

    end if
    ! lr408: ++++ END CONDUCTION GRADIENT COEFFS (SPAM3 MATRICES) ++++++++++++++


    ! Apply tensor correction: the result is not symmetric!
    ! cks: Determine matrices depending on if occ-precond happens
    ! cks: or not.
    call sparse_embed_array_product(tc_pur_denskern, pur_denskern, overlap)
    do is=1,pub_num_spins
       call sparse_embed_product(tc_qmat(is),qmat(is),overlap)
       if (pub_cond_calculate .and. .not. is_force_call) then ! lr408
          call sparse_embed_product(tc_cond_grad(is),cond_grad(is),overlap)
       end if
    end do

    ! cks: ++++ OCC-MIX (SPAM3 MATRICES) ++++++++++++++++++++++++++++++++++++
    ! cks: Mix fraction of occupation preconditioned covariant NGWF gradient
    ! cks: into covariant NGWF gradient.
    if (pub_occ_mix /= 0.0_DP) then

       call sparse_embed_create(tc_buffer,tc_qmat(1))


       ! cks: add pub_occ_mix*I to (1.0_DP-pub_occ_mix)*K.S
       call sparse_embed_array_scale(tc_pur_denskern, 1.0_DP-pub_occ_mix, pub_occ_mix)

       do is=1,pub_num_spins

          ! cks: Set tc_buffer to pub_occ_mix*(-1.0_DP)*S^-1.H
          call sparse_embed_product(tc_buffer,inv_overlap,ham(is))
          call sparse_embed_scale(tc_buffer,-pub_occ_mix)

          ! cks: set tc_qmat to pub_occ_mix*(-1.0_DP)*S^-1.H + (1.0_DP-pub_occ_mix)*Q_tc
          call sparse_embed_axpy(tc_buffer,tc_qmat(is),1.0_DP-pub_occ_mix)

          call sparse_embed_copy(tc_qmat(is),tc_buffer)

       end do

       call sparse_embed_destroy(tc_buffer)

       ! Occ-mixing for conduction NGWFs
       if (pub_cond_calculate .and. .not. is_force_call) then
          call sparse_embed_create(tc_buffer,tc_cond_grad(1))

          call sparse_embed_create(kt,val_dkn%m(1,1),cross_overlap)
          call sparse_embed_create(khv,val_dkn%m(1,1),val_ham(1))
          call sparse_embed_create(ksv,val_dkn%m(1,1),val_overlap)

          do is=1,pub_num_spins

             call sparse_embed_product(kt,val_dkn%m(is,PUB_1K),cross_overlap)
             call sparse_embed_product(khv,val_dkn%m(is,PUB_1K),val_ham(is))
             call sparse_embed_product(ksv,val_dkn%m(is,PUB_1K),val_overlap)
             call sparse_embed_scale(khv,-pub_occ_mix)
             call sparse_embed_axpy(khv,ksv,pub_occ_mix*cond_shift)

             ! lr408: Set tc_buffer to pub_occ_mix*(-KHvKT + wKSvKT)
             call sparse_embed_product(tc_buffer,khv,kt)

             ! lr408: set tc_qmat to
             ! lr408: pub_occ_mix*(-KHvKT + wKSvKT) + (1.0_DP-pub_occ_mix)*Q_cond_tc
             call sparse_embed_axpy(tc_buffer,tc_cond_grad(is),1.0_DP-pub_occ_mix)
             call sparse_embed_copy(tc_cond_grad(is),tc_buffer)
          end do

          call sparse_embed_destroy(khv)
          call sparse_embed_destroy(ksv)
          call sparse_embed_destroy(kt)

          call sparse_embed_destroy(tc_buffer)
       end if

    end if
    ! cks: ++++ END OCC-MIX (SPAM3 MATRICES) ++++++++++++++++++++++++++++++++

    ! ny: ++++++ SCISSOR GRADIENT COEFFS (SPAM3 MATRICES) ++++++++++++++++
    if (pub_scissor_ngroups > 0) then
       allocate(kern_array(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_gradient_coeffs','kern_array',ierr)

       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub
             call sparse_embed_extract_from_array(kern_array, &
                  pur_denskern%m(:,PUB_1K),isub,jsub)
             call sparse_create(shifted_denskern,kern_array(1))
             call sparse_create(sk,overlap%m(isub,jsub),kern_array(1))

             do is=1,pub_num_spins
                call scissor_shift_denskern(shifted_denskern, &
                     kern_array(is),ngwf_basis(isub),mdl%regions(isub), &
                     inv_overlap%m(isub,jsub))
                call sparse_product(sk,overlap%m(isub,jsub),kern_array(is))
                call sparse_product(scissor_mat(is)%m(isub,jsub), &
                     shifted_denskern,sk)
                call sparse_product(tc_scissor_mat(is)%m(isub,jsub), &
                     scissor_mat(is)%m(isub,jsub),overlap%m(isub,jsub))
             end do

             call sparse_destroy(shifted_denskern)
             call sparse_destroy(sk)
          enddo
       enddo

       deallocate(kern_array,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_coeffs','kern_array',ierr)
    end if
    ! ny: ++++++ END SCISSOR GRADIENT COEFFS (SPAM3 MATRICES) ++++++++++++

    ! ndmh: ++++ NONLOCAL GRADIENT COEFFS (SPAM3 MATRICES) ++++++++++++++++++
    if (pub_any_nl_proj.or.pub_aug) then

       ! Create the matrix of nonlocal energies D_ij
       if (pub_aug) then

          allocate(dij_ps_overlap(pub_num_spins),stat=ierr)
          call utils_alloc_check(myself,'dij_ps_overlap',ierr)
          allocate(dij(pub_num_spins),stat=ierr)
          call utils_alloc_check(myself,'dij',ierr)
          allocate(rhoij(pub_num_spins),stat=ierr)
          call utils_alloc_check(myself,'rhoij',ierr)
          if (pub_cond_calculate .and. .not. is_force_call) then
             allocate(rk_cond(pub_num_spins),stat=ierr)
             call utils_alloc_check(myself,'rk_cond',ierr)
             allocate(tc_rk_cond(pub_num_spins),stat=ierr)
             call utils_alloc_check(myself,'tc_rk_cond',ierr)
          end if
       else
          allocate(dij_ps_overlap(1),stat=ierr)
          call utils_alloc_check(myself,'dij_ps_overlap',ierr)
          allocate(dij(1),stat=ierr)
          call utils_alloc_check(myself,'dij',ierr)
       end if

       ! Norm-Conserving PSP version
       if (pub_any_nl_proj.and.(.not.pub_aug)) then

          ! Create matrix structures for D_ij and D_ij.<Proj_j|NGWF_b>
          call sparse_embed_create(dij_ps_overlap(1),ps_overlap)
          dij(1)%structure = 'E'
          ! agrecocmplx: initialise to complex in complex case
          call sparse_embed_create(dij(1),iscmplx=loc_cmplx)

          ! Get Kleinman-Bylander Denominators D_i (nonzero on diagonal only)
          do isub=1,mrows
             call pseudo_get_dij(dij(1)%m(isub,isub),mdl%regions(isub)%pseudo_sp)
          end do

          ! Calculate the matrix D_i.<Proj_i|NGWF_b>
          call sparse_embed_product(dij_ps_overlap(1),dij(1),ps_overlap)

          ! Calculate the matrices 4*(D_ij.<Proj_j|NGWF_b>).K^ab
          ! and 4*(D_ij.<Proj_j|NGWF_b>).K^bc.S_ca
          do is=1,pub_num_spins
             call sparse_embed_product(rk(is), dij_ps_overlap(1), &
                  pur_denskern%m(is,PUB_1K))
             call sparse_embed_scale(rk(is),2.0_DP*pub_spin_fac)
             call sparse_embed_product(tc_rk(is),dij_ps_overlap(1), &
                  tc_pur_denskern%m(is,PUB_1K))
             call sparse_embed_scale(tc_rk(is),2.0_DP*pub_spin_fac)
          end do

          ! Clean up temporary matrices
          call sparse_embed_destroy(dij(1))
          call sparse_embed_destroy(dij_ps_overlap(1))

       end if

       ! PAW or USP version
       if (pub_aug) then

          ! Create matrix structures for rho_ij, D_ij and D_ij.<Proj_j|NGWF_b>
          do is=1,pub_num_spins
             call sparse_embed_create(dij_ps_overlap(is),ps_overlap)
             rhoij(is)%structure = 'E'
             ! agrecocmplx: rhoij real in any case, need to stay real here
             ! because routine paw_nonlocal_energies does not work
             ! with complex matrices; dij can be complex, we need it to be
             ! complex here because of the call to sparse_product afterwards
             call sparse_embed_create(rhoij(is))
             call sparse_embed_create(dij(is),rhoij(is),iscmplx=loc_cmplx)
          end do
          ! agrecocmplx
          call sparse_embed_create(oij,rhoij(1),iscmplx=loc_cmplx)

          ! Calculate structures for valence NGWF contributions in COND mode
          if (pub_cond_calculate .and. .not. is_force_call) then
             call sparse_embed_transpose_structure(val_ps_overlap%structure, &
                  val_sp_overlap)
             ! agrecocmplx
             call sparse_embed_create(val_ps_overlap,iscmplx=loc_cmplx)
             call sparse_embed_transpose(val_ps_overlap,val_sp_overlap)
             call sparse_embed_create(oij_val_ps_overlap,oij,val_ps_overlap)

             do is=1,pub_num_spins
                call sparse_embed_create(rk_cond(is),oij_val_ps_overlap,cond_grad(is))
                call sparse_embed_create(tc_rk_cond(is),oij_val_ps_overlap, &
                     tc_cond_grad(is))
             end do
          end if

          ! Create projector density kernel
          allocate(rhoij_tmp(pub_num_spins),stat=ierr)
          call utils_alloc_check('ngwf_gradient_coeffs',&
               'rhoij_tmp',ierr)
          allocate(kern_array(pub_num_spins),stat=ierr)
          call utils_alloc_check('ngwf_gradient_coeffs',&
               'kern_array',ierr)
          call sparse_embed_extract_from_array(rhoij_tmp,rhoij)
          if (pub_cond_calculate .and. .not. is_force_call) then
             call sparse_embed_extract_from_array(kern_array,val_dkn%m(:,PUB_1K))
             call aug_projector_denskern(rhoij_tmp,kern_array,val_sp_overlap%p)
          else
             call sparse_embed_extract_from_array(kern_array,pur_denskern%m(:,PUB_1K))
             call aug_projector_denskern(rhoij_tmp,kern_array,sp_overlap%p)
          end if
          call sparse_embed_destroy_extracted_array(rhoij_tmp,rhoij,.true.)
          call sparse_embed_destroy_extracted_array(kern_array)
          deallocate(kern_array,stat=ierr)
          call utils_dealloc_check('ngwf_gradient_coeffs',&
               'kern_array',ierr)

          do is=1,pub_num_spins
             call sparse_embed_scale(rhoij(is),pub_spin_fac)
          end do

          ! Calculate nonlocal energies (block diagonal matrix)
          ! and add pre-calculated screened part
          ! jcap: EMBED_FIX!
          allocate(dij_tmp(pub_num_spins),stat=ierr)
          call utils_alloc_check('ngwf_gradient_coeffs',&
               'dij_tmp',ierr)
          do isub=1,mdl%nsub
             call sparse_embed_extract_from_array(dij_tmp,dij,isub,isub)
             call sparse_embed_extract_from_array(rhoij_tmp,rhoij,isub,isub)
             if (pub_paw) call paw_nonlocal_energies(dij_tmp,&
                  rhoij_tmp,mdl%regions(isub)%paw_sp, &
                  mdl%regions(isub)%par,show_matrices=.false.)
             call sparse_embed_destroy_extracted_array(dij_tmp,dij,.true.,&
                  isub,isub)
             call sparse_embed_destroy_extracted_array(rhoij_tmp,rhoij,.true.,&
                  isub,isub)
          end do
          deallocate(dij_tmp,stat=ierr)
          call utils_dealloc_check('ngwf_gradient_coeffs',&
               'dij_tmp',ierr)
          deallocate(rhoij_tmp,stat=ierr)
          call utils_dealloc_check('ngwf_gradient_coeffs',&
               'rhoij_tmp',ierr)
          do is=1,pub_num_spins
             call sparse_embed_axpy(dij(is),dijhat(is),1.0_DP)
          end do

          ! Calculate projector overlap matrix O_ij
          do isub=1,mdl%nsub
             if (pub_paw) then
                call paw_projector_overlap(oij%m(isub,isub),&
                     mdl%regions(isub)%paw_sp)
             else if (pub_usp) then
                call pseudo_aug_Q_matrix(oij%m(isub,isub),&
                     mdl%regions(isub)%pseudo_sp)
             end if
          end do

          ! Calculate the matrix O_ij.<Proj_i|NGWF_b>
          call sparse_embed_product(dij_ps_overlap(1),oij,ps_overlap)
          call sparse_embed_scale(dij_ps_overlap(1),2.0_DP*pub_spin_fac)
          do is=1,pub_num_spins
             ! Calculate the matrices 4*(O_ij.<Proj_j|NGWF_b>).Q^ab
             ! and 4*(O_ij.<Proj_j|NGWF_b>).Q^bc.S_ca
             call sparse_embed_product(rq(is),dij_ps_overlap(1),qmat(is))
             call sparse_embed_product(tc_rq(is),dij_ps_overlap(1),tc_qmat(is))
          end do

          if (pub_cond_calculate .and. .not. is_force_call) then
             ! Calculate the matrix O_ij.<Proj_i|val_NGWF_b>
             call sparse_embed_product(oij_val_ps_overlap,oij,val_ps_overlap)
             call sparse_embed_scale(oij_val_ps_overlap,2.0_DP*pub_spin_fac)
             do is=1,pub_num_spins
                ! Calculate matrices 4*(O_ij.<Proj_j|val_NGWF_b>).K(-HV+wSV)KTM
                ! and 4*(O_ij.<Proj_j|NGWF_b>).K(-HV+wSV)KTM
                call sparse_embed_product(rk_cond(is),oij_val_ps_overlap, &
                     cond_grad(is))
                call sparse_embed_product(tc_rk_cond(is),oij_val_ps_overlap, &
                     tc_cond_grad(is))
             end do
          end if

          do is=1,pub_num_spins
             ! Calculate the matrix D_ij.<Proj_j|NGWF_b>
             call sparse_embed_product(dij_ps_overlap(is),dij(is),ps_overlap)
             call sparse_embed_scale(dij_ps_overlap(is),2.0_DP*pub_spin_fac)

             ! Calculate the matrices 4*(D_ij.<Proj_j|NGWF_b>).K^ab
             ! and 4*(D_ij.<Proj_j|NGWF_b>).K^bc.S_ca
             call sparse_embed_product(rk(is), dij_ps_overlap(is), &
                  pur_denskern%m(is,PUB_1K))
             call sparse_embed_product(tc_rk(is),dij_ps_overlap(is), &
                  tc_pur_denskern%m(is,PUB_1K))
          end do

          ! All sets of terms will be evaluated together in one fftbox sum, so
          ! add 4*(O_ij.<Proj_j_NGWF_b>).Q^ab to rk(:) and
          ! add 4*(O_ij.<Proj_j|NGWF_b>).Q^bc.S_ca to tc_rk(:)
          do is=1,pub_num_spins
             ! ny: Add species scissor specific term
             if (pub_scissor_ngroups > 0) then
                call sparse_embed_create(oijr_scissor,rq(is))
                call sparse_embed_product(dij_ps_overlap(is),oij,ps_overlap)
                call sparse_embed_scale(dij_ps_overlap(is),2.0_DP*pub_spin_fac)
                call sparse_embed_product(oijr_scissor,dij_ps_overlap(is),scissor_mat(is))
                call sparse_embed_axpy(rq(is),oijr_scissor,1.0_DP)
                call sparse_embed_product(oijr_scissor,dij_ps_overlap(is),tc_scissor_mat(is))
                call sparse_embed_axpy(tc_rq(is),oijr_scissor,1.0_DP)
                call sparse_embed_destroy(oijr_scissor)
             end if

             call sparse_embed_axpy(rk(is),rq(is),1.0_DP)
             call sparse_embed_axpy(tc_rk(is),tc_rq(is),1.0_DP)
             if (pub_cond_calculate .and. .not. is_force_call) then
                call sparse_embed_axpy(rk(is),rk_cond(is),1.0_DP)
                call sparse_embed_axpy(tc_rk(is),tc_rk_cond(is),1.0_DP)
             end if
          end do

          ! Clean up temporary matrices
          if (pub_cond_calculate .and. .not. is_force_call) then
             call sparse_embed_destroy(oij_val_ps_overlap)
             call sparse_embed_destroy(val_ps_overlap)
          end if
          call sparse_embed_destroy(oij)
          do is=pub_num_spins,1,-1
             if (pub_cond_calculate .and. .not. is_force_call) then
                call sparse_embed_destroy(tc_rk_cond(is))
                call sparse_embed_destroy(rk_cond(is))
             end if
             call sparse_embed_destroy(dij(is))
             call sparse_embed_destroy(rhoij(is))
             call sparse_embed_destroy(dij_ps_overlap(is))
          end do
          if (pub_cond_calculate .and. .not. is_force_call) then
             deallocate(rk_cond,stat=ierr)
             call utils_dealloc_check(myself,'rk_cond',ierr)
             deallocate(tc_rk_cond,stat=ierr)
             call utils_dealloc_check(myself,'tc_rk_cond',ierr)
          end if
          deallocate(rhoij,stat=ierr)
          call utils_dealloc_check(myself,'rhoij',ierr)

       end if

       ! Deallocate SPAM3 arrays
       deallocate(dij,stat=ierr)
       call utils_dealloc_check(myself,'dij',ierr)
       deallocate(dij_ps_overlap,stat=ierr)
       call utils_dealloc_check(myself,'dij_ps_overlap',ierr)

    end if
    ! ndmh: ++++ END NONLOCAL GRADIENT COEFFS (SPAM3 MATRICES) ++++++++++++++

    ! ddor: ++++ DFT+U GRADIENT COEFFS (SPAM3 MATRICES) ++++++++++++++++++
    if (pub_hubbard) then

       ! Create matrix structure for [H^j_k O^ki + O^jk H_k^i] <Proj_i|NGWF_b>
       call sparse_embed_create(hub_ps_overlap,hub_overlap_t)

       ! Calculate the matrices 4*(G <Proj_i|NGWF_b>).K^ab
       ! and 4*(G <Proj_i|NGWF_b>).S_ac.K^cb

       ! agrecocmplx: need complex copy of hub_projector_ham
       ! to make the product in complex case
       if (loc_cmplx) then
          call sparse_create(hub_projector_ham_cmplx,hub%projector_ham(1),&
                  iscmplx=loc_cmplx)
       end if

       do is=1,pub_num_spins

          ! agrecocmplx
          if (loc_cmplx) then
             ! convert real to complex
             call sparse_copy(hub_projector_ham_cmplx,hub%projector_ham(is))
             call sparse_product(hub_ps_overlap%p,hub_projector_ham_cmplx, &
               hub_overlap_t%p)
          else
             ! Calculate the matrix G <Proj_i|NGWF_b>
             call sparse_product(hub_ps_overlap%p,hub%projector_ham(is), &
                  hub_overlap_t%p)
          end if

          ! ddor: Multiply on the right with Kernel and scale
          call sparse_embed_product(hub_wk(is), hub_ps_overlap, &
               pur_denskern%m(is,PUB_1K))
          call sparse_embed_scale(hub_wk(is),2.0_DP*pub_spin_fac)
          call sparse_embed_product(tc_hub_wk(is),hub_ps_overlap, &
               tc_pur_denskern%m(is,PUB_1K))
          call sparse_embed_scale(tc_hub_wk(is),2.0_DP*pub_spin_fac)
       end do

       ! agrecocmplx: destroy temporary matrix
       if (loc_cmplx) then
          call sparse_destroy(hub_projector_ham_cmplx)
       end if

       ! Clean up temporary matrices
       call sparse_embed_destroy(hub_ps_overlap)

    endif
    ! ddor: ++++ END DFT+U GRADIENT COEFFS (SPAM3 MATRICES) ++++++++++++++

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving '//myself//'.'

    call timer_clock(myself, 2)

  end subroutine ngwf_gradient_coeffs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwf_grad_matrix_ms_robust(gradient, &
       denskern, ham, overlap, mu)

    !=====================================================================!
    ! This subroutine returns the matrix that multiplies the fb-only      !
    ! part of the gradient of the LNV-MS energy as a function             !
    ! of the NGWF expansion coefficients in the psinc basis.              !
    !---------------------------------------------------------------------!
    ! Written by Arash A. Mostofi in March 2003, based on the subroutine  !
    ! ngwf_gradient_plain_matrix_for_fb_robust which was written by       !
    ! Chris-Kriton Skylaris on 21/6/2001 for the ONES code.               !
    ! Rewritten by Peter Haynes for block sparse matrices in spring 2004. !
    ! Modified for conduction calculations by Laura Ratcliff Oct 2010.    !
    ! Modified for embedding by Robert Charlton, 16/08/2018.              !
    !=====================================================================!

    use constants, only: DP, max_spins
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_copy, sparse_embed_transpose, sparse_embed_product, &
         sparse_embed_axpy, sparse_embed_scale

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: gradient(pub_num_spins)
    type(SPAM3_EMBED), intent(in)  :: denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in)  :: ham(pub_num_spins)
    type(SPAM3_EMBED), intent(in)  :: overlap
    real(kind=DP), intent(in)  :: mu(max_spins)

    ! Local Variables
    integer :: is
    type(SPAM3_EMBED) :: sl,ls,hl,lhl,ltmp

    ! Allocate workspace
    call sparse_embed_create(sl,overlap,denskern(1))
    call sparse_embed_create(ls,denskern(1),overlap)
    call sparse_embed_create(hl,ham(1),denskern(1))
    call sparse_embed_create(lhl,denskern(1),hl)
    call sparse_embed_create(ltmp,denskern(1))

    do is=1,pub_num_spins

       ! Calculate 'L' := 3*L.H.L - 2*(L.H.L.S.L + L.S.L.H.L) - mu*L
       call sparse_embed_product(sl,overlap,denskern(is))
       call sparse_embed_transpose(ls,sl)
       call sparse_embed_product(hl,ham(is),denskern(is))
       call sparse_embed_product(lhl,denskern(is),hl)
       call sparse_embed_copy(gradient(is),denskern(is))
       call sparse_embed_scale(gradient(is),-mu(is))
       call sparse_embed_product(ltmp,lhl,sl)
       call sparse_embed_axpy(gradient(is),ltmp,-2.0_DP)
       call sparse_embed_product(ltmp,ls,lhl)
       call sparse_embed_axpy(gradient(is),ltmp,-2.0_DP)
       call sparse_embed_axpy(gradient(is),lhl,3.0_DP)

    end do

    ! Deallocate workspace
    call sparse_embed_destroy(ltmp)
    call sparse_embed_destroy(lhl)
    call sparse_embed_destroy(hl)
    call sparse_embed_destroy(ls)
    call sparse_embed_destroy(sl)

  end subroutine ngwf_grad_matrix_ms_robust



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine ngwf_grad_matrix_lnv_robust(gradient, &
       purkern, denskern, ham, overlap, mu, norm_fac)

    !=====================================================================!
    ! This subroutine returns the matrix that multiplies the fb-only      !
    ! part of the gradient of the original-LNV energy as a function       !
    ! of the NGWF expansion coefficients in the psinc basis.              !
    ! @docme: args.                                                       !
    !---------------------------------------------------------------------!
    ! Written by Arash A. Mostofi in March 2003, based on the subroutine  !
    ! ngwf_gradient_plain_matrix_for_fb_robust which was written by       !
    ! Chris-Kriton Skylaris on 21/6/2001 for the ONES code.               !
    ! Rewritten by Peter Haynes for block sparse matrices in spring 2004. !
    ! Adjusted to allow independent rescaling of different spin channels  !
    ! of LNV gradient to match other parts of the code, by Nicholas Hine  !
    ! in April 2010.                                                      !
    ! Modified for conduction calculations by Laura Ratcliff Oct 2010.    !
    ! Modified for embedding by Robert Charlton, 20/07/17.                !
    !=====================================================================!

    use constants, only: DP, max_spins
    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_copy, sparse_embed_transpose, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: gradient(pub_num_spins)
    type(SPAM3_EMBED), intent(in)    :: purkern(pub_num_spins)
    type(SPAM3_EMBED), intent(in)    :: denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in)    :: ham(pub_num_spins)
    type(SPAM3_EMBED), intent(in)    :: overlap
    real(kind=DP), intent(in)        :: mu(max_spins)
    real(kind=DP), intent(in)        :: norm_fac(max_spins)

    ! Local Variables
    integer :: is
    type(SPAM3_EMBED) :: sl,hl,lhl,ktmp,htmp

    ! Allocate workspace
    call sparse_embed_create(sl,overlap,denskern(1))
    call sparse_embed_create(hl,ham(1),denskern(1))
    call sparse_embed_create(lhl,denskern(1),hl)
    call sparse_embed_create(ktmp,purkern(1))
    call sparse_embed_create(htmp,ham(1))

    ! pdh: new LNV requires gradient to be scaled by Nocc / Tr(KS)
    ! ndmh: alteration to allow separate rescaling of each spin channel
    do is=1,pub_num_spins

       ! Calculate 'L' := 3*L.H'.L - 2*(L.H'.L.S.L + L.S.L.H'.L) - mu*K
       ! ndmh: Removed two sparse_products and replaced them with transposes
       call sparse_embed_copy(htmp,ham(is))
       call sparse_embed_axpy(htmp,overlap,-mu(is))
       call sparse_embed_product(sl,overlap,denskern(is))
       call sparse_embed_product(hl,htmp,denskern(is))
       call sparse_embed_product(lhl,denskern(is),hl)
       call sparse_embed_product(gradient(is),lhl,sl)
       call sparse_embed_scale(gradient(is),-2.0_DP)
       call sparse_embed_transpose(ktmp,gradient(is))
       call sparse_embed_axpy(gradient(is),ktmp,1.0_DP)
       call sparse_embed_axpy(gradient(is),lhl,3.0_DP)

       ! pdh: re-scale gradient
       call sparse_embed_scale(gradient(is),norm_fac(is))

       ! ndmh: add mu*K term (already rescaled)
       call sparse_embed_axpy(gradient(is),purkern(is),-mu(is))

    end do

    ! Deallocate workspace
    call sparse_embed_destroy(htmp)
    call sparse_embed_destroy(ktmp)
    call sparse_embed_destroy(lhl)
    call sparse_embed_destroy(hl)
    call sparse_embed_destroy(sl)

  end subroutine ngwf_grad_matrix_lnv_robust


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef GPU_PGI_NGWFGRAD
  subroutine ngwf_gpu_gradient_batch(contra_grad, cov_grad, &  ! output
       fftbox_batch, lhxc_dbl, dbl_grid, rep, &
       ngwf_basis, proj_basis, hub_proj_basis, nl_projectors, hub, &
       fa_box_start, fa_start_in_box, batch_size, local_start, local_end, &
       coeff_mat, proj_coeff_mat, hub_proj_coeff_mat, cond_coeff_mat, &
       ps_overlap, max_current_size,is_force_call, &
       val_ngwf_basis, val_rep)

    !==========================================================================!
    ! This subroutine returns the contravariant and covariant NGWF             !
    ! gradients for the NGWFs of the current batch. It does this               !
    ! by applying the Hamiltonian operator to the functions accumulated        !
    ! in the fftbox of each batch, then applying kinetic energy                !
    ! preconditioning if required and then by extracting the                   !
    ! relevant ppds from the fftboxes and shaving their values                 !
    ! so that they are non-zero only within their spheres.                     !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! contra_grad (output) : Contravariant NGWF gradient in ppd representation !
    !    for the NGWFs of pub_my_proc_id.                                      !
    ! cov_grad (output)    : Covariant NGWF gradient in ppd representation     !
    !    for the NGWFs of pub_my_proc_id.                                      !
    ! fftbox_batch (input) :                                                   !
    ! lhxc_dbl (input) : Total local potential in the double grid for the      !
    !    whole simulation cell.                                                !
    ! rep (input) : NGWF Representation (functions and matrices).              !
    ! ngwf_basis (input) : Function basis describing the NGWFs                 !
    ! proj_basis (input) : Function basis describing the nonlocal projectors   !
    ! hub_proj_basis (input) : Function basis describing Hubbard projectors    !
    ! batch_size (input) : Size of each batch of NGWFs                         !
    ! local_start (input) : First NGWF in this batch on this proc              !
    ! local_end (input) : Last NGWF in this batch on this proc                 !
    ! coeff_mat (input) : Coefficients of NGWFs to add to gradient             !
    ! proj_coeff_mat (input) : Coefficients of projectors to add to gradient   !
    ! hub_proj_coeff_mat (input) : Coefficients of Hubbard projectors          !
    ! cond_coeff_mat (input) : Coefficients of valence NGWFs for gradient of   !
    !    conduction NGWFs (optional - conduction NGWF optimisations only).     !
    ! val_ngwf_basis (input) : Valence NGWF basis (optional - cond NGWF only). !
    ! val_rep  (input) : NGWF_REP type for valence NGWFs (optional as above.   !
    !--------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 26/12/2003                           !
    ! Modified for various speed enhancements by Nicholas Hine 2007-2009       !
    ! Modified for DFT+U by David O'Regan, 2009.                               !
    ! Modified to not use workspace_mod, Nicholas Hine, November 2009          !
    ! Tidied up, added more comments, reduced memory usage and made flow       !
    ! more logical, Nicholas Hine, October 2010                                !
    ! Branched from the original version to take advantage of GPUs by Karl     !
    ! Wilkinson, October 2012                                                  !
    !==========================================================================!

    use basis, only: basis_location_func_wrt_cell, &
         basis_extract_function_from_box, basis_clean_function
    use cell_grid, only: GRID_INFO, cell_grid_extract_box
    use comms, only : pub_on_root
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS
    use fourier, only: fourier_apply_box, fourier_gpu_interpolate_grad, &
         fourier_gpu_filter_grad
    use projectors, only: PROJECTOR_SET, projectors_gradient_batch
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_sum_fftbox_batch
    use hubbard_build, only: HUBBARD_MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_precond_real, pub_precond_recip, pub_task, &
         pub_any_nl_proj, pub_hubbard, pub_hubbard_restart, &
         pub_hubbard_atomsolve, pub_hub_on_the_fly, pub_paw, &
         pub_cond_calculate, pub_dbl_grid_scale, pub_smooth_scheme, &
         pub_debug_on_root, pub_num_spins, pub_spin_fac
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use fourier_gpu_wrapper_mod
    use accel_lib !! External dependency
    use cudafor   !! External dependency

    implicit none

    ! Arguments
    integer, intent(in) :: batch_size
    type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(PROJECTOR_SET), intent(inout) :: nl_projectors
    type(HUBBARD_MODEL), intent(inout) :: hub
    integer, intent(in) :: fa_box_start(3,batch_size)
    integer, intent(in) :: fa_start_in_box(3,batch_size)
    integer, intent(in) :: local_start, local_end
    integer, intent(in) :: max_current_size
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FUNCTIONS), intent(inout) ::contra_grad
    type(FUNCTIONS), intent(inout) :: cov_grad
    real(kind=DP), intent(inout) :: fftbox_batch(fftbox%total_ld1, &
         fftbox%total_ld2, fftbox%total_pt3, pub_num_spins, &
         4, batch_size)
    type(GRID_INFO), intent(in) :: dbl_grid
    real(kind=DP), intent(in) :: lhxc_dbl(dbl_grid%ld1, &
         dbl_grid%ld2, dbl_grid%max_group_slabs12, pub_num_spins)
    type(SPAM3), intent(in) :: coeff_mat(pub_num_spins,4)
    type(SPAM3), intent(in) :: proj_coeff_mat(pub_num_spins,4)
    type(SPAM3), intent(in) :: hub_proj_coeff_mat(pub_num_spins,4) !ddor
    type(SPAM3), intent(in) :: cond_coeff_mat(pub_num_spins,4) !lr408
    type(SPAM3), intent(in) :: ps_overlap
    logical, intent(in) :: is_force_call ! check if it is a cond task
     ! but called from forces mod
    type(FUNC_BASIS), optional, intent(in) :: val_ngwf_basis
    type(NGWF_REP), optional, intent(in) :: val_rep

    ! Local Variables
    logical :: i_need_potential
    integer :: prev_start1, prev_start2, prev_start3
    integer :: fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl
    integer :: local_fa
    integer :: batch_count
    integer :: is
    integer :: ierr
    integer :: idx_len
    integer, parameter :: fb_box = 1
    integer, parameter :: tc_fb_box = 2
    integer, parameter :: ham_box = 3
    integer, parameter :: tc_ham_box = 4
    real(kind=DP) :: common_fac
    integer, allocatable, dimension(:) :: overlap_idx
    real(kind=DP), allocatable, dimension(:,:,:) :: lhxc_fftbox_dbl
    real(kind=DP), allocatable, dimension(:,:,:) :: kin_buffer
    real(kind=DP), allocatable, dimension(:,:,:) :: tc_kin_buffer
    real(kind=DP), allocatable, dimension(:,:,:) :: buffer_dbl
    real(kind=DP), allocatable, dimension(:,:,:,:) :: dwork_box_dbl
    complex(kind=DP), allocatable, dimension(:,:,:) :: zwork_box

    real(kind=DP), allocatable, dimension(:,:,:,:) :: potential_box_dbl_pair
    real(kind=DP), allocatable, dimension(:,:,:), device :: kin_buffer_gpu
    real(kind=DP), allocatable, dimension(:,:,:), device :: tc_kin_buffer_gpu
    real(kind=DP) :: scalefac
    integer :: n1,n2,n3,i1,i2,i3,imat
    real(kind=DP) :: n1_dist, n2_dist, n3_dist, dist, dist_cutoff
    real(kind=DP) :: cell_d3, cell_d2, cell_d1
    integer :: n3max, n2max, n1max

    real(kind=DP) test_r(3)
    complex(kind=DP) test_c(3)
    integer :: npt, k1, k2, k3, istat

    ! Kaw: These arrays are a crude fix to a bug relating to the population of a
    !      complex array from two subarrays of the same array.
    real(kind=DP), allocatable, dimension(:,:,:,:,:), device :: fftbox_batch_gpu_ham
    real(kind=DP), allocatable, dimension(:,:,:,:,:), device :: fftbox_batch_gpu_tc_ham



    ! kaw: Allocate GPU arrays
    n1 = fftbox%total_ld1
    n2 = fftbox%total_ld2
    n3 = fftbox%total_pt3
    allocate(fftbox_batch_gpu(fftbox%total_ld1, fftbox%total_ld2, &
                              fftbox%total_pt3, pub_num_spins, 4, batch_size),stat=ierr)
    call utils_alloc_check('ngwf_gpu_gradient_batch','fftbox_batch_gpu',ierr)
    ! kaw: Used due to compiler error, see above. Should be done outside of loops
    !      over batches.
    allocate(fftbox_batch_gpu_ham(fftbox%total_ld1, fftbox%total_ld2, &
                                  fftbox%total_pt3, pub_num_spins, batch_size),stat=ierr)
    call utils_alloc_check('ngwf_gpu_gradient_batch','fftbox_batch_gpu_ham',ierr)
    allocate(fftbox_batch_gpu_tc_ham(fftbox%total_ld1, fftbox%total_ld2, &
                                  fftbox%total_pt3, pub_num_spins, batch_size),stat=ierr)
    call utils_alloc_check('ngwf_gpu_gradient_batch','fftbox_batch_gpu_tc_ham',ierr)
    scalefac = 1.0_DP / (fftbox%total_ld1*fftbox%total_ld2*fftbox%total_pt3)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering ngwf_gpu_gradient_batch'

    ! pdh: overlap matrix index
    idx_len = sparse_index_length(rep%ngwf_overlap)
    allocate(overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('ngwf_gpu_gradient_batch','overlap_idx',ierr)

    ! ndmh: calculate \sum_\beta \phi_\beta (r) * C^i_\alpha\beta in
    ! ndmh: FFT boxes for the four C^i matrices.
    ! ndmh: Once this is done:
    ! ndmh: Box 1 (fb_box) contains \sum_b Q^ab phi_b
    ! ndmh: Box 2 (tc_fb_box) contains \sum_b (QS)^a_b phi_b
    ! ndmh: Box 3 (ham_box) contains \sum_b K^ab phi_b
    ! ndmh: Box 4 (tc_ham_box) contains \sum_b (KS)^a_b phi_b
    call sparse_generate_index(overlap_idx,rep%ngwf_overlap)
    common_fac = 2.0_DP * cell%weight * pub_spin_fac
    call function_ops_sum_fftbox_batch(fftbox_batch, &
         rep%ngwfs_on_grid, ngwf_basis, ngwf_basis, fa_box_start, batch_size, &
         local_start, local_end, overlap_idx, idx_len, coeff_mat, &
         fb_box, tc_ham_box, common_fac)


    ! kaw: Copy fftbox batch onto GPU
    fftbox_batch_gpu = fftbox_batch
    fftbox_batch_gpu_ham = fftbox_batch(:,:,:,:,3,:)
    fftbox_batch_gpu_tc_ham = fftbox_batch(:,:,:,:,4,:)

    ! kaw: Copy kinetic data to GPU - this would ideally be done outside the
    !      loop over batches
    fftbox_recip_grid_gpu = fftbox%recip_grid(5,:,:,:)

    deallocate(overlap_idx,stat=ierr)
    call utils_dealloc_check('ngwf_gpu_gradient_batch','overlap_idx',ierr)

    ! ndmh: allocate generic complex FFTbox workspace
    allocate(zwork_box(fftbox%total_ld1,fftbox%total_ld2,&
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('ngwf_gpu_gradient_batch','zwork_box',ierr)


    ! KKKKKKKKKKK KINETIC KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

    allocate(kin_buffer_gpu(fftbox%total_ld1,fftbox%total_ld2,&
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('ngwf_gpu_gradient_batch','kin_buffer_gpu',ierr)
    allocate(tc_kin_buffer_gpu(fftbox%total_ld1,fftbox%total_ld2, &
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('ngwf_gpu_gradient_batch','tc_kin_buffer_gpu',ierr)

!$acc data region
    do is=1,pub_num_spins
       batch_count = 0
       do local_fa=local_start,local_start+max_current_size+1

          batch_count = batch_count + 1

          if (local_fa <= local_end) then

!$acc region
             !populate complex array
             ! kaw: This is the specific code relating to the compiler error mentioned above.
             !coarse_work_gpu(1:n1,1:n2,1:n3) = &
             !   cmplx(fftbox_batch_gpu(1:n1,1:n2,1:n3,is,ham_box,batch_count), &
             !         fftbox_batch_gpu(1:n1,1:n2,1:n3,is,tc_ham_box,batch_count),kind=DP)
             ! kaw: And this is the temporary fix:
             coarse_work_gpu(1:n1,1:n2,1:n3) = &
                cmplx(fftbox_batch_gpu_ham(1:n1,1:n2,1:n3,is,batch_count), &
                      fftbox_batch_gpu_tc_ham(1:n1,1:n2,1:n3,is,batch_count),kind=DP)

!$acc end region

             call cufftExec(cufftplan_coarse,planType,coarse_work_gpu, &
                            coarse_work_gpu,CUFFT_FORWARD)

             ! apply the kinetic energy operator
!$acc region
             do i3=1,fftbox%total_pt3
                do i2=1,fftbox%total_pt2
                   do i1=1,fftbox%total_pt1
                      coarse_work_gpu(i1,i2,i3) = coarse_work_gpu(i1,i2,i3) * &
                           fftbox_recip_grid_gpu(i1,i2,i3) * scalefac
                   end do
                end do
             end do
!$acc end region

             call cufftExec(cufftplan_coarse,planType,coarse_work_gpu, &
                            coarse_work_gpu,CUFFT_INVERSE)

!$acc region
             do i3=1,fftbox%total_pt3
                do i2=1,fftbox%total_pt2
                   do i1=1,fftbox%total_pt1
                      kin_buffer_gpu(i1,i2,i3) = real(coarse_work_gpu(i1,i2,i3),DP)
                      tc_kin_buffer_gpu(i1,i2,i3) = aimag(coarse_work_gpu(i1,i2,i3))
                   end do
                end do
             end do
!$acc end region

             ! ndmh: add T \sum_b fb to the fb_box
!$acc region
             fftbox_batch_gpu(:,:,:,is,fb_box,batch_count) = &
                  fftbox_batch_gpu(:,:,:,is,fb_box,batch_count) + &
                  kin_buffer_gpu(:,:,:)

             fftbox_batch_gpu(:,:,:,is,tc_fb_box,batch_count) = &
                  fftbox_batch_gpu(:,:,:,is,tc_fb_box,batch_count) + &
                  tc_kin_buffer_gpu(:,:,:)
!$acc end region

          end if

       end do
    end do
!$acc end data region

    deallocate(kin_buffer_gpu,stat=ierr)
    call utils_dealloc_check('ngwf_gpu_gradient_batch','kin_buffer_gpu',ierr)
    deallocate(tc_kin_buffer_gpu,stat=ierr)
    call utils_dealloc_check('ngwf_gpu_gradient_batch','tc_kin_buffer_gpu',ierr)

    ! KKKKKKKK END KINETIC KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

    ! LLLLLL LHXC POTENTIAL LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    allocate(potential_box_dbl_pair(fftbox%total_ld1_dbl, &
             fftbox%total_ld2_dbl, fftbox%total_pt3_dbl, 2),stat=ierr)
    call utils_alloc_check('ngwf_gpu_gradient_batch','potential_box_dbl_pair',ierr)
    allocate(lhxc_fftbox_dbl(fftbox%total_ld1_dbl,&
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl),stat=ierr)
    call utils_alloc_check('ngwf_gpu_gradient_batch','lhxc_fftbox_dbl',ierr)
    allocate(buffer_dbl(fftbox%total_ld1_dbl,fftbox%total_ld2_dbl,&
         dbl_grid%max_group_slabs12),stat=ierr)
    call utils_alloc_check('ngwf_gpu_gradient_batch','buffer_dbl',ierr)

    ! cks: initialisation
    buffer_dbl = 0.0_DP

    ! pdh: loop over spins
    do is=1,pub_num_spins

       prev_start1= -1111111
       prev_start2= -2222222
       prev_start3= -3333333
       batch_count = 0
       do local_fa=local_start,local_start+max_current_size+1
          batch_count = batch_count + 1

          ! ndmh: find FFTbox start position on double grid
          if (local_fa <= local_end) then
             if (pub_dbl_grid_scale > 1.0_DP) then
                fftbox_start1_dbl = 2*fa_box_start(1,batch_count) - 1
                fftbox_start2_dbl = 2*fa_box_start(2,batch_count) - 1
                fftbox_start3_dbl = 2*fa_box_start(3,batch_count) - 1
             else ! ndmh: dbl_grid is scale 1.0
                fftbox_start1_dbl = fa_box_start(1,batch_count)
                fftbox_start2_dbl = fa_box_start(2,batch_count)
                fftbox_start3_dbl = fa_box_start(3,batch_count)
             end if
          else
             fftbox_start1_dbl = -1234
             fftbox_start2_dbl = -1234
             fftbox_start3_dbl = -1234
          end if

          ! Extract potential from cell to FFT box if position of FFT box
          ! is different from last time.
          i_need_potential = .false.
          if ((fftbox_start1_dbl /= prev_start1 .or. &
               fftbox_start2_dbl /= prev_start2 .or. &
               fftbox_start3_dbl /= prev_start3) .and. local_fa <= local_end) &
               i_need_potential = .true.

          call cell_grid_extract_box(lhxc_fftbox_dbl,&
               buffer_dbl, lhxc_dbl(:,:,:,is), dbl_grid, &
               fftbox%total_pt1_dbl, fftbox%total_pt2_dbl, &
               fftbox%total_pt3_dbl, fftbox%total_ld1_dbl, &
               fftbox%total_ld2_dbl, fftbox_start1_dbl, &
               fftbox_start2_dbl, fftbox_start3_dbl, i_need_potential, .true.)

          prev_start1 = fftbox_start1_dbl
          prev_start2 = fftbox_start2_dbl
          prev_start3 = fftbox_start3_dbl

          if (local_fa <= local_end) then

             ! kaw: There is a potential saving here as the same value is used
             !      for both...
             potential_box_dbl_pair(:,:,:,1) = lhxc_fftbox_dbl
             potential_box_dbl_pair(:,:,:,2) = lhxc_fftbox_dbl

             ! Calculate V_loc on sum of phi_beta
             if (pub_dbl_grid_scale>1.0_DP) then
                ! ndmh: interpolate function sums to fine grid, multiply
                ! ndmh: by fine grid potential, then filter back to coarse grid

                !kaw: These routines need to be modified in order to allow them
                !     to use the "input" array that is already on the GPU
                call fourier_gpu_interpolate_grad(potential_box_dbl_pair, &
                      .true., is, batch_count)
                call fourier_gpu_filter_grad(is, batch_count, .true.)
             else

                ! ndmh: dbl_grid_scale:1, so no need to interpolate
!$acc data region
!$acc region
                fftbox_batch_gpu(:,:,:,is,ham_box,batch_count) = &
                     fftbox_batch_gpu(:,:,:,is,ham_box,batch_count) * &
                     potential_box_dbl_gpu(:,:,:,1)
                fftbox_batch_gpu(:,:,:,is,tc_ham_box,batch_count) = &
                     fftbox_batch_gpu(:,:,:,is,tc_ham_box,batch_count) * &
                     potential_box_dbl_gpu(:,:,:,2)
!$acc end region
!$acc end data region
             end if

!$acc data region
!$acc region
             ! cks: grad = Vloc*fb + fb
             fftbox_batch_gpu(:,:,:,is,ham_box,batch_count) = &
                  fftbox_batch_gpu(:,:,:,is,ham_box,batch_count) + &
                  fftbox_batch_gpu(:,:,:,is,fb_box,batch_count)
             fftbox_batch_gpu(:,:,:,is,tc_ham_box,batch_count) = &
                  fftbox_batch_gpu(:,:,:,is,tc_ham_box,batch_count) + &
                  fftbox_batch_gpu(:,:,:,is,tc_fb_box,batch_count)
!$acc end region
!$acc end data region

          end if

       end do

    end do


    deallocate(buffer_dbl,stat=ierr)
    call utils_dealloc_check('ngwf_gpu_gradient_batch','buffer_dbl',ierr)
    deallocate(lhxc_fftbox_dbl,stat=ierr)
    call utils_dealloc_check('ngwf_gpu_gradient_batch','lhxc_fftbox_dbl',ierr)

    deallocate(potential_box_dbl_pair,stat=ierr)
    call utils_dealloc_check('ngwf_gpu_gradient_batch','potential_box_dbl_pair',&
         ierr)
    ! LLLLLL END LHXC POTENTIAL LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    ! kaw: Copy fftbox batch back from the GPU, this needs to be shifted as other
    !      areas of the code are accelerated.
    fftbox_batch = fftbox_batch_gpu

    ! NNNNNN NON-LOCAL POTENTIAL NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    if (pub_any_nl_proj.or.pub_aug) then
       call projectors_gradient_batch(fftbox_batch, ham_box, tc_ham_box, &
            local_start, local_end, batch_size, fa_start_in_box, &
            ngwf_basis, proj_basis, proj_coeff_mat, ps_overlap, nl_projectors)
    end if
    ! NNNNNN END NON-LOCAL POTENTIAL NNNNNNNNNNNNNNNNNNNNNNNNNNNN

    ! HHHHHH HUBBARD DFT+U POTENTIAL HHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if (pub_hubbard) then
       if (((pub_task == 'HUBBARDSCF') .and. (hub%consistency_iteration > 1)) .or.&
            & pub_hubbard_restart .or. pub_hubbard_atomsolve .or. &
            & pub_hub_on_the_fly) then

          idx_len = sparse_index_length(rep%hub_overlap_t)
          allocate(overlap_idx(idx_len),stat=ierr)
          call utils_alloc_check('ngwf_gpu_gradient_batch','overlap_idx',ierr)

          call sparse_generate_index(overlap_idx,rep%hub_overlap_t)
          common_fac = cell%weight
          call function_ops_sum_fftbox_batch(fftbox_batch, &
               hub%consistency_projs, hub_proj_basis, ngwf_basis, &
               fa_box_start, batch_size, local_start, &
               local_end, overlap_idx, idx_len, hub_proj_coeff_mat, &
               ham_box, tc_ham_box, common_fac)

          deallocate(overlap_idx,stat=ierr)
          call utils_dealloc_check('ngwf_gpu_gradient_batch','overlap_idx',ierr)
       else
          call projectors_gradient_batch(fftbox_batch, ham_box, tc_ham_box, &
               local_start, local_end, batch_size, fa_start_in_box, &
               ngwf_basis, hub_proj_basis, hub_proj_coeff_mat, &
               rep%hub_overlap_t, hub%projectors)
       endif
    endif
    ! HHHHHH END HUBBARD DFT+U POTENTIAL HHHHHHHHHHHHHHHHHHHHHHHHH


    ! CCCCCCC CONDUCTION TERM CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    if (pub_cond_calculate .and. .not. is_force_call) then
       ! lr408: Indexing for cross overlap matrix
       idx_len = sparse_index_length(rep%ngwf_cross_overlap)
       allocate(overlap_idx(idx_len),stat=ierr)
       call utils_alloc_check('ngwf_gpu_gradient_batch','overlap_idx',ierr)

       call sparse_generate_index(overlap_idx,rep%ngwf_cross_overlap)
       common_fac = 2.0_DP * cell%weight * pub_spin_fac
       call function_ops_sum_fftbox_batch(fftbox_batch, &
            val_rep%ngwfs_on_grid, val_ngwf_basis, ngwf_basis, fa_box_start, &
            batch_size, local_start, local_end, overlap_idx, idx_len, &
            cond_coeff_mat, ham_box, tc_ham_box, common_fac)

       deallocate(overlap_idx,stat=ierr)
       call utils_dealloc_check('ngwf_gpu_gradient_batch','overlap_idx',ierr)
    end if
    ! CCCCCCC END CONDUCTION TERM CCCCCCCCCCCCCCCCCCCCCCCCCCCCC


    ! RRRRRRR RECIPROCAL SPACE PRECONDITIONING RRRRRRRRRRRRRRRR
    if (pub_precond_recip) then
       do is=1,pub_num_spins
          batch_count = 1
          do local_fa=local_start,local_end,2

             ! Copy the covariant gradient into complex workspace array
             if (local_fa < local_end) then
                zwork_box = &
                     cmplx(fftbox_batch(:,:,:,is,tc_ham_box,batch_count), &
                     fftbox_batch(:,:,:,is,tc_ham_box,batch_count+1),kind=DP)
             else
                zwork_box = &
                     cmplx(fftbox_batch(:,:,:,is,tc_ham_box,batch_count), &
                     0.0_DP,kind=DP)
             end if

             ! Forward FFT the covariant gradient to reciprocal space
             call fourier_apply_box('Coarse','Forward',zwork_box)

             ! Apply kinetic energy preconditioning to covariant gradient
             call ngwf_gradient_precond_recip(zwork_box)

             ! Backward FFT the covariant gradient to real space
             call fourier_apply_box('Coarse','Backward',zwork_box)

             ! Copy the preconditioned covariant gradient out of the complex
             ! workspace array
             fftbox_batch(:,:,:,is,tc_ham_box,batch_count) = &
                  real(zwork_box,kind=DP)
             if (local_fa < local_end) &
                  fftbox_batch(:,:,:,is,tc_ham_box,batch_count+1) = &
                  aimag(zwork_box)

             batch_count = batch_count + 2
          end do
       end do

    end if
    ! RRRRRRR RECIPROCAL SPACE PRECONDITIONING RRRRRRRRRRRRRRRR

    ! ndmh: APPLY REAL SPACE PRECONDITIONING, APPLY SPIN-AVERAGING,
    ! ndmh: THEN EXTRACT AND SHAVE GRADIENT
    batch_count = 0
    do local_fa=local_start,local_end
       batch_count = batch_count + 1

       if (pub_precond_real) then
          do is=1,pub_num_spins
             ! cks: precondition the covariant gradient in real space
             call ngwf_grad_precond_rspace(  &
                  fftbox_batch(:,:,:,is,tc_ham_box,batch_count), & ! in/out,
                  fa_start_in_box(1,batch_count), &
                  fa_start_in_box(2,batch_count), &
                  fa_start_in_box(3,batch_count), &
                  ngwf_basis%tight_boxes(local_fa)%tight_pts1, &
                  ngwf_basis%tight_boxes(local_fa)%tight_pts2, &
                  ngwf_basis%tight_boxes(local_fa)%tight_pts3)
          end do
       end if

       ! ndmh: Average gradients over spins for spin-polarised calculations.
       if (pub_num_spins == 2) then
          fftbox_batch(:,:,:,1,ham_box,batch_count) = &
               fftbox_batch(:,:,:,1,ham_box,batch_count) + &
               fftbox_batch(:,:,:,2,ham_box,batch_count)
          fftbox_batch(:,:,:,1,tc_ham_box,batch_count) = &
               fftbox_batch(:,:,:,1,tc_ham_box,batch_count) + &
               fftbox_batch(:,:,:,2,tc_ham_box,batch_count)
       end if

       if (pub_smooth_scheme .ne. 'NONE') then
          ! smmd: smooth the covariant gradient in real space
          call ngwf_grad_smooth_rspace(  &
               fftbox_batch(:,:,:,1,tc_ham_box,batch_count), & ! in/out,
               fa_start_in_box(1,batch_count), &
               fa_start_in_box(2,batch_count), &
               fa_start_in_box(3,batch_count), &
               ngwf_basis%tight_boxes(local_fa), &
               ngwf_basis%spheres(local_fa)%centre, &
               ngwf_basis%spheres(local_fa)%radius)
       endif



       ! cks: shaving - stage 1 / extract ppds from the FFT box
       ! extract the contravariant gradient from ham_box to PPDs
       call basis_extract_function_from_box(contra_grad, &
            fftbox_batch(:,:,:,1,ham_box,batch_count), &
            ngwf_basis%spheres(local_fa),&
            ngwf_basis%tight_boxes(local_fa), fa_start_in_box(1,batch_count), &
            fa_start_in_box(2,batch_count), fa_start_in_box(3,batch_count), &
            ngwf_basis%spheres(local_fa)%offset,cell,fftbox)

       ! extract the covariant gradient from ham_box to PPDs
       call basis_extract_function_from_box(cov_grad, &
            fftbox_batch(:,:,:,1,tc_ham_box,batch_count), &
            ngwf_basis%spheres(local_fa),&
            ngwf_basis%tight_boxes(local_fa), fa_start_in_box(1,batch_count), &
            fa_start_in_box(2,batch_count), fa_start_in_box(3,batch_count), &
            ngwf_basis%spheres(local_fa)%offset,cell,fftbox)

       ! cks: shaving - stage 2 / zero points outside NGWF sphere in PPD rep.
       ! smmd: if smoothing scheme applied previously, no needs to clean the
       ! smmd: covariant gradient (nor the contravariant gradient which is
       ! smmd: always used within contra.DOT.cov expressions)
       if (cell%n_pts > 1 .and. pub_smooth_scheme .eq. 'NONE') then
          call basis_clean_function(contra_grad, &
               ngwf_basis%spheres(local_fa), ngwf_basis%n_ppds,cell,fftbox)
          call basis_clean_function(cov_grad, &
               ngwf_basis%spheres(local_fa), ngwf_basis%n_ppds,cell,fftbox)
       endif

    end do
    ! ndmh: END APPLY REAL SPACE PRECONDITIONING, APPLY SPIN-AVERAGING,
    ! ndmh: THEN EXTRACT AND SHAVE GRADIENT

    ! ndmh: deallocate generic workspace
    deallocate(zwork_box,stat=ierr)
    call utils_dealloc_check('ngwf_gpu_gradient_batch','zwork_box',ierr)


    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving ngwf_gpu_gradient_batch'

  deallocate(fftbox_batch_gpu,stat=ierr)
  call utils_dealloc_check('ngwf_gpu_gradient_batch','fftbox_batch_gpu',ierr)

  end subroutine ngwf_gpu_gradient_batch
#else
  subroutine dummy_gpu
  use fourier_gpu_wrapper_mod
  !kaw: This is a dummy use statement to make the check_dependencies script happy.
  end subroutine dummy_gpu
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_gradient_paw_precond_init(proj_set,paw_sp,cell,fftbox)

    !==================================================================!
    ! This subroutine returns the matrix used in the PAW-correction to !
    ! the kinetic energy preconditioning:                              !
    ! proj_set%O_prec = -O*(I+C*O)^-1                                  !
    ! with O - partial wave overlap difference                         !
    ! and C - <p^i|(1+T)^-1|p^j>                                       !
    !------------------------------------------------------------------!
    ! This subroutine was written by Gabriel Constantinescu in May'14. !
    !==================================================================!

    use fft_box, only: FFTBOX_INFO
    use paw, only: PAW_SPECIES, paw_species_calc_proj_prec_mat
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_precond_recip, pub_precond_real, pub_aug, pub_paw
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(PROJECTOR_SET), intent(inout) :: proj_set
    type(PAW_SPECIES), intent(in) :: paw_sp(:)

    if (.not.pub_aug) return

    ! ndmh: allocate preconditioner
    if (pub_precond_recip.or.pub_precond_real) then
       call ngwf_grad_init_precond_recip(fftbox)
    end if

    if (pub_paw) then
       call paw_species_calc_proj_prec_mat(proj_set,cell,fftbox, &
            precond_func_recip,paw_sp)
    end if

    call ngwf_gradient_exit

  end subroutine ngwf_gradient_paw_precond_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_gradient_paw_precond_exit(proj_set)

    !==================================================================!
    ! This subroutine deallocates storage for the O_prec matrix in     !
    ! the projector set.                                               !
    !------------------------------------------------------------------!
    ! This subroutine was written by Gabriel Constantinescu in May'14. !
    !==================================================================!

    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_aug
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(PROJECTOR_SET), intent(inout) :: proj_set

    ! Local Variables
    integer :: ierr

    if (.not.pub_aug) return

    ! deallocate O_prec matrix
    if (allocated(proj_set%O_prec)) then
       deallocate(proj_set%O_prec,stat=ierr)
       call utils_dealloc_check('paw_species_init_proj', &
            'proj_set%O_prec',ierr)
    end if

  end subroutine ngwf_gradient_paw_precond_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwf_gradient_update_kinet_kp(kin_buf, box_data, &
      zwork_box, fftbox, kcart) ! agrecokpt

    !=========================================================================!
    ! This subroutine adds the extra terms of the kinetic energy arising in   !
    ! the KP method for non-Gamma k-points.                                   !
    ! Currently works for a single k-point.                                   !
    !                                                                         !
    ! Arguments:                                                              !
    !     kin_buf : on-entry, original kinetic operator applied to |kets>     !
    !     box_data : NGWFs in FFTBOX_DATA format                              !
    !     kcart : current k-point in cartesian components                     !
    !                                                                         !
    ! Written by Andrea Greco on 19/04/2016.                                  !
    !=========================================================================!

    use datatypes, only: FFTBOX_DATA
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use kinetic, only: kinetic_grad_on_box
    use rundat, only: pub_threads_fftbox
    use utils, only: utils_assert, utils_alloc_check, &
            utils_dealloc_check

    ! Argument
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FFTBOX_DATA), intent(inout) :: kin_buf
    type(FFTBOX_DATA), intent(inout) :: box_data
    complex(kind=DP), intent(out) :: zwork_box(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    ! agrecokpt
    real(kind=DP), intent(in) :: kcart(3)

    ! Local variables
    integer :: ierr
    integer :: dim_cart
    complex(kind=DP), dimension(:,:,:,:), allocatable :: zwork_grad_box

    ! agrecokpt: only complex NGWFs are supported in this case
    call utils_assert(kin_buf%iscmplx, &
         'Subroutine ngwf_gradient_update_kinet_kp supports&
         & only complex NGWFs')

    allocate(zwork_grad_box(fftbox%total_ld1,fftbox%total_ld2, &
         fftbox%total_pt3,3),stat=ierr)
    call utils_alloc_check('ngwf_gradient_update_kinet_kp','zwork_grad_box',ierr)

    ! FFT to reciprocal space
    call fourier_apply_box('Coarse', 'Forward', box_data%z, &
            gspc=zwork_box,omp=pub_threads_fftbox)

    ! apply grad operator in reciprocal space
    call kinetic_grad_on_box(zwork_box,zwork_grad_box,fftbox)

    do dim_cart=1,3
       ! FFT back to real space
       call fourier_apply_box('Coarse','Backward', &
            zwork_grad_box(:,:,:,dim_cart), omp=pub_threads_fftbox)
    end do

    ! KE(Gamma)|ket> - i k * grad|ket>
    do dim_cart=1,3
       ! add KP term
       kin_buf%z(:,:,:) = kin_buf%z(:,:,:) + &
           cmplx(0.0_DP,-kcart(dim_cart),kind=DP) * &
           zwork_grad_box(:,:,:,dim_cart)
    end do

    ! KE(Gamma)|ket> - i k * grad|ket> + 0.5*|k|^2 |ket>
    kin_buf%z(:,:,:) = kin_buf%z(:,:,:) + &
        0.5_DP*(kcart(1)*kcart(1) + kcart(2)*kcart(2) + &
        kcart(3)*kcart(3))*box_data%z(:,:,:)

    ! deallocate workspace
    deallocate(zwork_grad_box,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_update_kinet_kp','zwork_grad_box',ierr)

  end subroutine ngwf_gradient_update_kinet_kp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwf_gradient_build_batch(fftbox_batch, &! output
       lhxc_dbl, dbl_grid, &!rep, &
       ngwf_basis, proj_basis, hub_proj_basis, nl_projectors, hub, &
       fa_box_start, fa_start_in_box, batch_size, local_start, local_end, &
       coeff_mat, proj_coeff_mat, cond_coeff_mat, &
       ps_overlap, fftbox, cell, max_current_size, is_force_call, &
       val_ngwf_basis, val_ngwfs_on_grid, ngwf_cross_overlap, kpt, dfdtau_dbl, &
       hub_proj_coeff_mat, hub_overlap_t)

    !==========================================================================!
    ! This subroutine applies the Hamiltonian operator to the functions        !
    ! accumulated in the fftbox of each batch, accounting for the different    !
    ! Hamiltonians that may exist in different regions.                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! fftbox_batch (output) :                                                  !
    ! lhxc_dbl (input) : Total local potential in the double grid for the      !
    !    whole simulation cell.                                                !
    ! ngwf_basis (input) : Function basis describing the NGWFs                 !
    ! proj_basis (input) : Function basis describing the nonlocal projectors   !
    ! hub_proj_basis (input) : Function basis describing Hubbard projectors    !
    ! nl_projectors (input) : Projector set for nonlocal pseudo projectors     !
    ! hub (input) : Storage for Hubbard Model information                      !
    ! fa_box_start (input) : Start positions of the FFTboxes for each function !
    ! fa_start_in_box (input) : Start positions of functions in FFTboxes       !
    ! batch_size (input) : Size of each batch of NGWFs                         !
    ! local_start (input) : First NGWF in this batch on this proc              !
    ! local_end (input) : Last NGWF in this batch on this proc                 !
    ! coeff_mat (input) : Coefficients of NGWFs to add to gradient             !
    ! proj_coeff_mat (input) : Coefficients of projectors to add to gradient   !
    ! hub_proj_coeff_mat (input) : Coefficients of Hubbard projectors          !
    ! cond_coeff_mat (input) : Coefficients of valence NGWFs for gradient of   !
    !    conduction NGWFs (optional - conduction NGWF optimisations only).     !
    ! val_ngwf_basis (input) : Valence NGWF basis (optional - cond NGWF only). !
    ! val_rep  (input) : NGWF_REP type for valence NGWFs (optional as above.   !
    ! dfdtau_dbl (inout)     : Gradient of XC energy per unit volume wrt       !
    !                          KE density energy on double grid (optional)     !
    !--------------------------------------------------------------------------!
    ! Extracted from ngwf_gradient_batch by Robert Charlton, June 2018.        !
    !==========================================================================!

    use basis, only: basis_extract_function_from_box, basis_clean_function
    use cell_grid, only: GRID_INFO, cell_grid_extract_box
    use constants, only: DP, stdout, EDA_POLFRAGLOC_DEVEL
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_copy, &
         data_fftbox_alloc, data_fftbox_dealloc, data_fftbox_axpy
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_sum_fftbox_batch
    use geometry, only: POINT ! agrecokpt
    use hubbard_build, only: HUBBARD_MODEL
    use ngwf_representation, only: NGWF_REP
    use projectors, only: PROJECTOR_SET, projectors_gradient_batch, &
         projectors_grad_precond_batch
    use rundat, only: pub_precond_recip, pub_task, &
         pub_any_nl_proj, pub_hubbard, pub_hubbard_restart, &
         pub_hubbard_atomsolve, pub_hub_on_the_fly, pub_aug, pub_paw, &
         pub_cond_calculate, pub_debug, pub_debug_on_root, &
         pub_dbl_grid_scale, pub_num_spins, &
         pub_threads_fftbox, pub_spin_fac, &
         pub_devel_code, &
         pub_xc_ke_density_required ! JCW
!$  use rundat, only: pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index, &
        sparse_create, sparse_copy, sparse_destroy
    use sparse_embed, only: SPAM3_EMBED
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: batch_size
    !type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(PROJECTOR_SET), intent(inout) :: nl_projectors
    type(HUBBARD_MODEL), intent(inout) :: hub
    integer, intent(in) :: fa_box_start(3,batch_size)!,size(ngwf_basis))
    integer, intent(in) :: fa_start_in_box(3,batch_size)!,size(ngwf_basis))
    integer, intent(in) :: local_start, local_end
    !integer, intent(in) :: local_start_reg(size(ngwf_basis)), &
    !  local_end_reg(size(ngwf_basis))
    integer, intent(in) :: max_current_size
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FFTBOX_DATA), intent(inout) :: fftbox_batch(pub_num_spins,4,batch_size)
    type(GRID_INFO), intent(in) :: dbl_grid
    real(kind=DP), intent(in) :: lhxc_dbl(dbl_grid%ld1, &
         dbl_grid%ld2, dbl_grid%max_group_slabs12, pub_num_spins)
    type(SPAM3), intent(in) :: coeff_mat(pub_num_spins,4)
    type(SPAM3), intent(in) :: proj_coeff_mat(pub_num_spins,4)
    type(SPAM3), intent(in) :: cond_coeff_mat(pub_num_spins,4) !lr408
    ! rc2013: make this optional so we don't need to mess around with hub
    type(SPAM3), intent(in), optional :: hub_proj_coeff_mat(pub_num_spins,4) !ddor
    type(SPAM3), intent(in), optional :: hub_overlap_t
    type(SPAM3), intent(in) :: ps_overlap
    type(SPAM3), intent(in), optional :: ngwf_cross_overlap
    logical, intent(in) :: is_force_call ! Check if this is cond task but
      ! function is called from forces_mod
    type(FUNC_BASIS), optional, intent(in) :: val_ngwf_basis
    type(FUNCTIONS), optional, intent(in)  :: val_ngwfs_on_grid
    ! JCW: Optional gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density, on double grid. A pointer will be assigned to this if
    ! JCW: present.
    real(kind=DP), optional, target, intent(in) :: dfdtau_dbl(dbl_grid%ld1, &
         dbl_grid%ld2, dbl_grid%max_group_slabs12, pub_num_spins)
    ! agrecokpt: argument to specify a k-point
    type(POINT), optional, intent(in) :: kpt

    ! Local Variables
    logical :: i_need_potential
    integer :: prev_start1, prev_start2, prev_start3
    integer :: fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl
    integer :: local_fa!, local_start, local_end
    integer :: batch_count
    integer :: is
    integer :: ierr
    integer :: idx_len
    integer, parameter :: fb_box = 1
    integer, parameter :: tc_fb_box = 2
    integer, parameter :: ham_box = 3
    integer, parameter :: tc_ham_box = 4
    real(kind=DP) :: common_fac
    integer, allocatable, dimension(:) :: overlap_idx
    real(kind=DP), allocatable, dimension(:,:,:) :: lhxc_fftbox_dbl
    real(kind=DP), allocatable, dimension(:,:,:) :: vtau_fftbox_dbl
    ! agrecocmplx: use FFTBOX_DATA to switch between real and complex case
    !real(kind=DP), allocatable, dimension(:,:,:) :: kin_buffer
    !real(kind=DP), allocatable, dimension(:,:,:) :: tc_kin_buffer
    type(FFTBOX_DATA) :: kin_buffer
    type(FFTBOX_DATA) :: tc_kin_buffer
    !real(kind=DP), allocatable, dimension(:,:,:) :: vtau_buffer
    !real(kind=DP), allocatable, dimension(:,:,:) :: tc_vtau_buffer
    type(FFTBOX_DATA) :: vtau_buffer
    type(FFTBOX_DATA) :: tc_vtau_buffer
    real(kind=DP), allocatable, dimension(:,:,:) :: buffer_dbl
    real(kind=DP), allocatable, dimension(:,:,:,:) :: dwork_box_dbl
    complex(kind=DP), allocatable, dimension(:,:,:) :: coarse_work,fine_work
    ! JCW: Additional arrays required for tau-dependent contribution
    real(kind=DP), pointer, dimension(:,:,:,:) :: vtau_dbl
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: dwork_grad_box
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: dwork_grad_box_dbl
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: zwork_grad_box

    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt
    type(POINT) :: loc_kpt

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering ngwf_gradient_build_batch'

    call timer_clock('ngwf_gradient_build_batch',1)

    ! agrecocmplx
    loc_cmplx = coeff_mat(1,1)%iscmplx

    ! agrecokpt
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    ! agrecocmplx: should be ok now, need to check!
    !jmecmplx
    !call utils_assert(.not. loc_cmplx, 'Error in&
    !     & ngwf_gradient_batch: not ready yet for complex NGWFs.')

    ! JCW: If dfdtau_dbl is present, make vtau_dbl point to it
    if (present(dfdtau_dbl)) then
       vtau_dbl => dfdtau_dbl
    else
       nullify(vtau_dbl)
    end if

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(local_fa, batch_count, i_need_potential, is, &
!$OMP      fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl, &
!$OMP      coarse_work, dwork_box_dbl, &
!$OMP      lhxc_fftbox_dbl, vtau_fftbox_dbl, buffer_dbl, ierr, &
!$OMP      prev_start1, prev_start2, prev_start3,fine_work, &
!$OMP      dwork_grad_box, dwork_grad_box_dbl, zwork_grad_box) &
!$OMP FIRSTPRIVATE(kin_buffer, tc_kin_buffer, vtau_buffer, &
!$OMP      tc_vtau_buffer) &
!$OMP SHARED(fa_start_in_box, cell, local_start, fa_box_start, &
!$OMP      local_end, max_current_size, fftbox, dbl_grid, &
!$OMP      ngwf_basis, fftbox_batch, pub_dbl_grid_scale, lhxc_dbl, &
!$OMP      pub_num_spins, pub_threads_num_fftboxes,loc_cmplx,loc_kpt, &
!$OMP      pub_xc_ke_density_required,pub_devel_code,pub_debug_on_root, &
!$OMP      vtau_dbl)

    if (pub_xc_ke_density_required) then
       ! JCW: Check that vtau_dbl pointer is associated when a tau-dependent
       ! JCW: XC functional is used. This should be set to point to optional
       ! JCW: dfdtau_dbl argument, if present.
       call utils_assert(associated(vtau_dbl),"Error in ngwf_gradient_lnv: &
            &vtau_dbl pointer should be associated if pub_xc_ke_density_required &
            &is .true.")
    end if

    ! ndmh: allocate generic complex FFTbox workspace
    allocate(coarse_work(fftbox%total_ld1,fftbox%total_ld2,&
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('ngwf_gradient_build_batch','coarse_work',ierr)


    ! KKKKKKKKKKK KINETIC KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

    ! agrecocmplx
    !allocate(kin_buffer(fftbox%total_ld1,fftbox%total_ld2,&
    !     fftbox%total_pt3),stat=ierr)
    !call utils_alloc_check('ngwf_gradient_build_batch','kin_buffer',ierr)
    !allocate(tc_kin_buffer(fftbox%total_ld1,fftbox%total_ld2, &
    !     fftbox%total_pt3),stat=ierr)
    !call utils_alloc_check('ngwf_gradient_build_batch','tc_kin_buffer',ierr)
    call data_fftbox_alloc(kin_buffer,fftbox%total_ld1,fftbox%total_ld2,&
            fftbox%total_pt3,iscmplx=loc_cmplx)
    call data_fftbox_alloc(tc_kin_buffer,fftbox%total_ld1,fftbox%total_ld2,&
            fftbox%total_pt3,iscmplx=loc_cmplx)

    do is=1,pub_num_spins
!$OMP DO
       do local_fa=local_start,local_start+max_current_size-1
          batch_count = local_fa - local_start + 1

          if (local_fa <= local_end) then
             ! Calculate T on sum of phi_beta
             ! agrecocmplx: this is now compatible with new types
             ! agrecokpt: added extra term for KP method
             call ngwf_gradient_kinetic(kin_buffer, tc_kin_buffer, &
                  fftbox_batch(is,ham_box,batch_count),  &
                  fftbox_batch(is,tc_ham_box,batch_count), &
                  coarse_work,fftbox,kpt=loc_kpt)

             ! ndmh: add T \sum_b fb to the fb_box
             ! agrecocmplx: this is now compatible with complex types

             ! ny: use fft_box axpy routine
             call data_fftbox_axpy(fftbox_batch(is,fb_box,batch_count), &
                  kin_buffer)
             call data_fftbox_axpy(fftbox_batch(is,tc_fb_box,batch_count), &
                  tc_kin_buffer)
          end if

       end do
!$OMP END DO
    end do

    ! agrecocmplx
    call data_fftbox_dealloc(tc_kin_buffer)
    call data_fftbox_dealloc(kin_buffer)

    ! KKKKKKKK END KINETIC KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK



    ! LLLLLL LHXC POTENTIAL LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    ! JCW: Includes tau-dependent part of XC potential, if required

    allocate(lhxc_fftbox_dbl(fftbox%total_ld1_dbl,&
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl),stat=ierr)
    call utils_alloc_check('ngwf_gradient_build_batch','lhxc_fftbox_dbl',ierr)
    if (pub_xc_ke_density_required) then
       ! JCW: Allocate space for gradient component containing tau-dependent
       ! JCW: part of XC potential
       call utils_assert(.not.loc_cmplx,"Error in ngwf_gradient_batch: &
            &complex NGWFs not yet supported for tau-dependent functionals")
       allocate(vtau_fftbox_dbl(fftbox%total_ld1_dbl,&
            fftbox%total_ld2_dbl,fftbox%total_pt3_dbl),stat=ierr)
       call utils_alloc_check('ngwf_gradient_batch','vtau_fftbox_dbl',ierr)
       ! JCW: Allocate temporary space for vtau contribution to gradient
       call data_fftbox_alloc(vtau_buffer,fftbox%total_ld1,fftbox%total_ld2,&
               fftbox%total_pt3,iscmplx=loc_cmplx)
       call data_fftbox_alloc(tc_vtau_buffer,fftbox%total_ld1,fftbox%total_ld2,&
               fftbox%total_pt3,iscmplx=loc_cmplx)
       ! JCW: Allocate workspace for vtau contribution
       allocate(dwork_grad_box(fftbox%total_ld1, fftbox%total_ld2,&
            fftbox%total_pt3,3,2),stat=ierr)
       call utils_alloc_check('ngwf_gradient_batch','dwork_grad_box',ierr)
       allocate(dwork_grad_box_dbl(fftbox%total_ld1_dbl, fftbox%total_ld2_dbl,&
            fftbox%total_pt3_dbl,3,2),stat=ierr)
       call utils_alloc_check('ngwf_gradient_batch','dwork_grad_box',ierr)
       allocate(zwork_grad_box(fftbox%total_ld1, fftbox%total_ld2,&
            fftbox%total_pt3,3,2),stat=ierr)
       call utils_alloc_check('ngwf_gradient_batch','zwork_grad_box',ierr)
    end if
    allocate(buffer_dbl(fftbox%total_ld1_dbl,fftbox%total_ld2_dbl,&
         dbl_grid%max_group_slabs12),stat=ierr)
    call utils_alloc_check('ngwf_gradient_build_batch','buffer_dbl',ierr)
    allocate(dwork_box_dbl(fftbox%total_ld1_dbl,fftbox%total_ld2_dbl,&
         fftbox%total_pt3_dbl,2),stat=ierr)
    call utils_alloc_check('ngwf_gradient_build_batch','dwork_box_dbl',ierr)
    allocate(fine_work(fftbox%total_ld1_dbl,fftbox%total_ld2_dbl,&
         fftbox%total_pt3_dbl),stat=ierr)
    call utils_alloc_check('ngwf_gradient_build_batch','fine_work',ierr)

    ! cks: initialisation
    buffer_dbl = 0.0_DP

    ! pdh: loop over spins
    do is=1,pub_num_spins

       prev_start1= -1111111
       prev_start2= -2222222
       prev_start3= -3333333

!$OMP DO
       do local_fa=local_start,local_start+max_current_size-1
          batch_count = local_fa - local_start + 1

          ! ndmh: find FFTbox start position on double grid
          if (local_fa <= local_end) then
             if (pub_dbl_grid_scale > 1.0_DP) then
                fftbox_start1_dbl = 2*fa_box_start(1,batch_count) - 1
                fftbox_start2_dbl = 2*fa_box_start(2,batch_count) - 1
                fftbox_start3_dbl = 2*fa_box_start(3,batch_count) - 1
             else ! ndmh: dbl_grid is scale 1.0
                fftbox_start1_dbl = fa_box_start(1,batch_count)
                fftbox_start2_dbl = fa_box_start(2,batch_count)
                fftbox_start3_dbl = fa_box_start(3,batch_count)
             end if
          else
             fftbox_start1_dbl = -1234
             fftbox_start2_dbl = -1234
             fftbox_start3_dbl = -1234
          end if

          ! Extract potential from cell to FFT box if position of FFT box
          ! is different from last time.
          i_need_potential = .false.
          if ((fftbox_start1_dbl /= prev_start1 .or. &
               fftbox_start2_dbl /= prev_start2 .or. &
               fftbox_start3_dbl /= prev_start3) .and. local_fa <= local_end) &
               i_need_potential = .true.

!$OMP CRITICAL
          call cell_grid_extract_box(lhxc_fftbox_dbl,&
               buffer_dbl, lhxc_dbl(:,:,:,is), dbl_grid, &
               fftbox%total_pt1_dbl, fftbox%total_pt2_dbl, &
               fftbox%total_pt3_dbl, fftbox%total_ld1_dbl, &
               fftbox%total_ld2_dbl, fftbox_start1_dbl, &
               fftbox_start2_dbl, fftbox_start3_dbl, i_need_potential, .true.)
          if (pub_xc_ke_density_required) then
             ! JCW: Extract FFTbox of vtau (dfdtau) potential
             call cell_grid_extract_box(vtau_fftbox_dbl,&
                  buffer_dbl, vtau_dbl(:,:,:,is), dbl_grid, &
                  fftbox%total_pt1_dbl, fftbox%total_pt2_dbl, &
                  fftbox%total_pt3_dbl, fftbox%total_ld1_dbl, &
                  fftbox%total_ld2_dbl, fftbox_start1_dbl, &
                  fftbox_start2_dbl, fftbox_start3_dbl, i_need_potential, .true.)
          end if
!$OMP END CRITICAL

          prev_start1 = fftbox_start1_dbl
          prev_start2 = fftbox_start2_dbl
          prev_start3 = fftbox_start3_dbl

          if (local_fa <= local_end) then

             ! Calculate V_loc on sum of phi_beta
             if (pub_dbl_grid_scale>1.0_DP) then
                if (pub_xc_ke_density_required) then
                   ! JCW: Copy \sum_b K^ab phi_b and \sum_b (KS)^a_b phi_b from
                   ! JCW: fftbox_batch to vtau_buffer and tc_vtau_buffer prior
                   ! JCW: to the modification of fftbox_batch by
                   ! JCW: ngwf_gradient_local.
                   ! TODO Implement support for complex NGWFs
                   call utils_assert(.not.loc_cmplx,"Error in ngwf_gradient_batch: &
                        &ngwf_gradient_div_locpot_grad does not support complex &
                        &NGWFs.")
                   call data_fftbox_copy(vtau_buffer,fftbox_batch(is,ham_box,batch_count))
                   call data_fftbox_copy(tc_vtau_buffer,fftbox_batch(is,tc_ham_box,batch_count))
                   !vtau_buffer%d(:,:,:) = fftbox_batch(is,ham_box,batch_count)%d(:,:,:)
                   !tc_vtau_buffer%d(:,:,:) = fftbox_batch(is,tc_ham_box,batch_count)%d(:,:,:)

                   ! JCW: Calculate contribution for tau-dependent part of XC potential
                   call ngwf_gradient_div_locpot_grad( &
                      vtau_buffer%d, &                              ! input-output
                      tc_vtau_buffer%d,&                            ! input-output
                      vtau_fftbox_dbl, vtau_fftbox_dbl, &         ! input-output
                      coarse_work,fine_work,&                     ! workspace
                      dwork_grad_box,dwork_grad_box_dbl,&         ! workspace
                      zwork_grad_box,&                            ! workspace
                      fftbox)
                   if (pub_debug_on_root) then
                      ! Test only works for single MPI process, so use
                      ! pub_debug_on_root
                      ! To activate test, need
                      !   devel_code: NGWFGRAD:TAUTEST=T:NGWFGRAD
                      ! in input file.
                      ! Directly check pub_devel_code, rather than using
                      ! utils_devel_code, since this does a comms_bcast
                      call utils_assert(.not.loc_cmplx,"Error in ngwf_gradient_batch: &
                           &Devel code NGWFGRAD:TAUTEST=T:NGWFGRAD does not &
                           &support complex NGWFs")
                      ! TODO Implement support for complex NGWFs
                      if (index(&
                         pub_devel_code,"NGWFGRAD:TAUTEST=T:NGWFGRAD")>0) then
                         call test_ngwf_gradient_div_locpot_grad(&
                              fftbox_batch(is,ham_box,batch_count),&
                              fftbox_batch(is,tc_ham_box,batch_count),&
                              vtau_buffer,tc_vtau_buffer, &
                              vtau_fftbox_dbl, vtau_fftbox_dbl, &
                              coarse_work,fine_work,&
                              dwork_grad_box,dwork_grad_box_dbl,&
                              zwork_grad_box,&
                              fftbox)
                      end if
                   end if
                end if
                ! ndmh: interpolate function sums to fine grid, multiply
                ! ndmh: by fine grid potential, then filter back to coarse grid
                ! agrecocmplx: now compatible with complex NGWFs
                call ngwf_gradient_local( &
                     fftbox_batch(is,ham_box,batch_count), &   ! input-output
                     fftbox_batch(is,tc_ham_box,batch_count),& ! input-output
                     lhxc_fftbox_dbl, lhxc_fftbox_dbl, &             ! input-output
                     dwork_box_dbl,coarse_work,fine_work,fftbox)     ! workspace
             else
                ! ndmh: dbl_grid_scale:1, so no need to interpolate
                call utils_assert(.not.pub_xc_ke_density_required,"Error in &
                     &ngwf_gradient_build_batch: dbl_grid_scale==1.0 not implemented/&
                     &tested with tau-dependent XC functionals")
                ! agrecocmplx
                if (loc_cmplx) then
                   fftbox_batch(is,ham_box,batch_count)%z(:,:,:) = &
                        fftbox_batch(is,ham_box,batch_count)%z(:,:,:) * &
                        lhxc_fftbox_dbl(:,:,:)
                   fftbox_batch(is,tc_ham_box,batch_count)%z(:,:,:) = &
                        fftbox_batch(is,tc_ham_box,batch_count)%z(:,:,:) * &
                        lhxc_fftbox_dbl(:,:,:)
                else
                   fftbox_batch(is,ham_box,batch_count)%d(:,:,:) = &
                        fftbox_batch(is,ham_box,batch_count)%d(:,:,:) * &
                        lhxc_fftbox_dbl(:,:,:)
                   fftbox_batch(is,tc_ham_box,batch_count)%d(:,:,:) = &
                        fftbox_batch(is,tc_ham_box,batch_count)%d(:,:,:) * &
                        lhxc_fftbox_dbl(:,:,:)
                end if
             end if

             ! cks: grad = Vloc*fb + fb
             ! ny: use fft_box axpy routine
             call data_fftbox_axpy(fftbox_batch(is,ham_box,batch_count), &
                  fftbox_batch(is,fb_box,batch_count))
             call data_fftbox_axpy(fftbox_batch(is,tc_ham_box,batch_count), &
                  fftbox_batch(is,tc_fb_box,batch_count))
             if (pub_xc_ke_density_required) then
                ! TODO Implement support for complex NGWFs
                call utils_assert(.not.loc_cmplx,"Error in ngwf_gradient_build_batch: &
                     &complex NGWFs not yet supported for tau-dependent functionals")

                ! JCW: Add contribution for tau-dependent part of XC potential
                ! ny: use fft_box axpy routine
                call data_fftbox_axpy(fftbox_batch(is,ham_box,batch_count), &
                     vtau_buffer,-0.5_DP)
                call data_fftbox_axpy(fftbox_batch(is,tc_ham_box,batch_count),&
                     tc_vtau_buffer,-0.5_DP)

             end if
          end if

       end do
!$OMP END DO
!$OMP BARRIER
    end do

    deallocate(fine_work,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_build_batch','fine_work',ierr)
    deallocate(dwork_box_dbl,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_build_batch','dwork_box_dbl',ierr)
    deallocate(buffer_dbl,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_build_batch','buffer_dbl',ierr)
    if (pub_xc_ke_density_required) then
       ! JCW: Deallocate temporary space for vtau contribution to gradient
       call data_fftbox_dealloc(tc_vtau_buffer)
       call data_fftbox_dealloc(vtau_buffer)
       ! JCW: Deallocate space for gradient component containing tau-dependent
       ! JCW: part of XC potential
       deallocate(vtau_fftbox_dbl,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_batch','vtau_fftbox_dbl',ierr)
       ! JCW: Deallocate workspace for vtau contribution
       deallocate(dwork_grad_box,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_batch','dwork_grad_box',ierr)
       deallocate(dwork_grad_box_dbl,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_batch','dwork_grad_box',ierr)
       deallocate(zwork_grad_box,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_batch','zwork_grad_box',ierr)
    end if
    deallocate(lhxc_fftbox_dbl,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_build_batch','lhxc_fftbox_dbl',ierr)
    deallocate(coarse_work,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_build_batch','coarse_work',ierr)
!$OMP END PARALLEL

    ! LLLLLL END LHXC POTENTIAL LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    ! NNNNNN NON-LOCAL POTENTIAL NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    if (pub_any_nl_proj.or.pub_aug) then
       ! agrecokpt: call with specified k-point
       call projectors_gradient_batch(fftbox_batch, ham_box, tc_ham_box, &
            local_start, local_end, batch_size, fa_box_start, &
            ngwf_basis, proj_basis, cell, fftbox, proj_coeff_mat, &
            ps_overlap, nl_projectors,kpt=loc_kpt)
    end if
    ! NNNNNN END NON-LOCAL POTENTIAL NNNNNNNNNNNNNNNNNNNNNNNNNNNN


    ! HHHHHH HUBBARD DFT+U POTENTIAL HHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if (pub_hubbard) then
       if (((pub_task == 'HUBBARDSCF') .and. (hub%consistency_iteration > 1)) .or.&
            & pub_hubbard_restart .or. pub_hubbard_atomsolve .or. &
            & pub_hub_on_the_fly) then

          idx_len = sparse_index_length(hub_overlap_t)
          allocate(overlap_idx(idx_len),stat=ierr)
          call utils_alloc_check('ngwf_gradient_batch','overlap_idx',ierr)

          call sparse_generate_index(overlap_idx,hub_overlap_t)
          common_fac = cell%weight
          call function_ops_sum_fftbox_batch(fftbox_batch, fftbox, cell, &
               hub%consistency_projs, hub_proj_basis, &
               fa_box_start, batch_size, local_start, &
               local_end, overlap_idx, idx_len, hub_proj_coeff_mat, &
               ham_box, tc_ham_box, common_fac)

          deallocate(overlap_idx,stat=ierr)
          call utils_dealloc_check('ngwf_gradient_batch','overlap_idx',ierr)
       else
          ! agrecokpt: k-point dependence needs to be included here as well?
          ! probably yes...
          call projectors_gradient_batch(fftbox_batch, ham_box, tc_ham_box, &
               local_start, local_end, batch_size, fa_box_start, &
               ngwf_basis, hub_proj_basis, cell, fftbox, &
               hub_proj_coeff_mat, hub_overlap_t, hub%projectors, &
               kpt=loc_kpt)
       endif
    endif
    ! HHHHHH END HUBBARD DFT+U POTENTIAL HHHHHHHHHHHHHHHHHHHHHHHHH


    ! CCCCCCC CONDUCTION TERM CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    if (pub_cond_calculate .and. .not. is_force_call) then

       ! lr408: Indexing for cross overlap matrix
       idx_len = sparse_index_length(ngwf_cross_overlap)
       allocate(overlap_idx(idx_len),stat=ierr)
       call utils_alloc_check('ngwf_gradient_batch','overlap_idx',ierr)

       call sparse_generate_index(overlap_idx,ngwf_cross_overlap)
       common_fac = 2.0_DP * cell%weight * pub_spin_fac
       call function_ops_sum_fftbox_batch(fftbox_batch, fftbox, cell, &
            val_ngwfs_on_grid, val_ngwf_basis, fa_box_start, &
            batch_size, local_start, local_end, overlap_idx, idx_len, &
            cond_coeff_mat, ham_box, tc_ham_box, common_fac)

       deallocate(overlap_idx,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_batch','overlap_idx',ierr)

    end if
    ! CCCCCCC END CONDUCTION TERM CCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    ! PPPPPPP PAW PRECONDITIONER PROJ-OVERLAP TERM PPPPPPPPPPPP
    if (pub_paw.and.pub_precond_recip) then
       ! Calculate overlap of this batch of NGWF covariant gradients
       ! with preconditioned projectors
       ! Multiply <g_a|q_i> matrix by O^prec_ij to get coefficients
       ! of preconditioned projectors to add back into gradient
       ! agrecokpt: does this need to be called with kpt argument?
       ! probably yes since we are using non-local projectors which
       ! are k-point dependent in the KP method
       call projectors_grad_precond_batch(fftbox_batch,tc_ham_box, &
            local_start,local_end,batch_size,fa_box_start, &
            cell,fftbox,ngwf_basis,ps_overlap,nl_projectors, &
            precond_func_recip, kpt=loc_kpt)
    end if
    ! PPPPPPP END PAW PRECONDITIONER PROJ-OVERLAP TERM PPPPPPPP

    ! JCW: Nullify vtau_dbl pointer, if associated
    if (associated(vtau_dbl)) nullify(vtau_dbl)

    call timer_clock('ngwf_gradient_build_batch',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving ngwf_gradient_build_batch'

  end subroutine ngwf_gradient_build_batch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module ngwf_gradient
