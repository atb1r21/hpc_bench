! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                        Hamiltonian Matrix Module                            !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   Module for evaluating the Hamiltonian matrix and related other matrices,  !
!   eg overlap and inverse overlap, in various formats, in the NGWF           !
!   Representation. Uses the NGWF_HAM and NGWF_REP types defined in           !
!   ngwf_representation_mod.                                                  !
!                                                                             !
!   The following public routines are available:                              !
!                                                                             !
!   hamiltonian_dens_indep_matrices: assembles the DENSity INDEPendent terms  !
!       of the Hamiltonian matrix, namely the kinetic matrix, nonlocal matrix !
!       (for NCPPs), overlap, overlaps with any other representations, and    !
!       inverse overlap.                                                      !
!                                                                             !
!   hamiltonian_dens_dep_matrices: assembles the DENSity DEPendent terms of   !
!       the Hamiltonian matrix. Takes in density kernel and calculates the    !
!       electronic density, and thence the effective potential (local+hartree !
!       +xc) matrix, and also calculates HFX exchange and DFT+U contributions !
!       (should only be called for the valence NGWFs)                         !
!                                                                             !
!   hamiltonian_dens_dep_nonsc: assembles the DENSity DEPendent terms in a    !
!       NON Self-Consistent way. Takes in an already-calculated LHXC          !
!       effective potential, and calculates the LHXC matrix and any other     !
!       kernel-dependent terms (assumes fixed valence density, hence ideal    !
!       for Conduction NGWF calculations etc)                                 !
!                                                                             !
!   hamiltonian_lhxc_calculate: standalone routine for calculating lhxc       !
!       potential                                                             !
!                                                                             !
!   hamiltonian_build_matrix: takes the matrices calculated above, and adds   !
!       together the relevant terms to produce the total Hamiltonian matrix   !
!                                                                             !
!   hamiltonian_proj_cond_matrix: takes the conduction NGWF Hamiltonian and   !
!       shifts-and-scales the valence states so that they are raised above    !
!       all the conduction states in energy.                                  !
!                                                                             !
!   hamiltonian_energy_components: evaluates all the individual contributions !
!       to the total energy of the system.                                    !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   The subroutines in this file were written by:                             !
!                                                                             !
!   Chris-Kriton Skylaris, Nicholas D.M. Hine, Jacek Dziedzic,                !
!   David O'Regan, Laura Ratcliff, and Andrea Greco.                          !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module hamiltonian

  use constants, only: DP
  use rundat, only: pub_debug_on_root

  private

  public :: hamiltonian_dens_indep_matrices
  public :: hamiltonian_dens_dep_matrices
  public :: hamiltonian_dens_dep_nonsc
  public :: hamiltonian_lhxc_calculate
  public :: hamiltonian_build_matrix
  public :: hamiltonian_proj_cond_matrix
  public :: hamiltonian_energy_components
  ! agrecokpt: to be used in KP method
  public :: hamiltonian_update_kinet_kp
  ! agrecokpt: to be used in TB method
  public :: hamiltonian_apply_phases_tb_method
  public :: hamiltonian_soc_matrices
  public :: hamiltonian_proj_embed_matrix

contains

  subroutine hamiltonian_dens_indep_matrices(rep, &                       ! out
       ngwf_basis,projector_basis,nl_projectors,hub_proj_basis,hub,mdl, & ! in
       val_rep, val_ngwf_basis, kpt, allow_pseuinvS)                      ! in, opt

    !=============================================================!
    ! This subroutine initialises the matrices which make up the  !
    ! Hamiltonian and are independent of the density kernel.      !
    !-------------------------------------------------------------!
    ! Arguments:                                                  !
    ! rep             (inout) : NGWF Representation (functions    !
    !                           and matrices).                    !
    ! ngwf_basis      (input) : Function basis type for the NGWFs !
    ! projector_basis (input) : Function basis type for nonlocal  !
    !                           pseudopotential projectors        !
    ! nl_projectors   (input) : Type describing the nonlocal      !
    !                           projectors                        !
    ! hub             (input) : Type describing the Hubbard Model !
    ! mdl             (input) : Container type containing general !
    !                           information about the system      !
    ! hub_proj_basis  (input) : Function basis type for Hubbard   !
    !                           projectors                        !
    ! val_rep         (input) : Valence NGWF Representation       !
    !                           !optional)                        !
    ! val_ngwf_basis  (input) : Function basis type for the       !
    !                           valence NGWFs (optional)          !
    !-------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 3/7/2001.               !
    ! DFT+U added by David O'Regan, April 2008                    !
    ! Adapted for SPAM3, function basis and new pseudopotential   !
    ! routines by Nicholas Hine, May-July 2009.                   !
    ! Moved to hamiltonian_mod, April 2010.                       !
    ! Calculation of valence-conduction overlap matrix added,     !
    ! with necessary optional arguments by Laura Ratcliff,        !
    ! Oct 2010.                                                   !
    ! Modified by Andrea Greco on 21/06/2015 to use complex NGWFs.!
    ! Modified for embedding by Robert Charlton, 26/05/2017.      !
    !=============================================================!

    use augmentation, only: augmentation_overlap
    use comms, only: pub_on_root
    use constants, only: DP, stdout, VERBOSE
    use dense, only: DEM,dense_create,dense_destroy,dense_convert, &
         dense_invert
    use dftb, only: dftb_analytical_gaussian_overlap
    use function_basis, only: FUNC_BASIS
    use fragment_matrix_ops, only: fmo_construct_overlap_mats
    use fragment_scfmi, only: overlap_f_dens
    use hubbard_build, only: HUBBARD_MODEL, hubbard_projection_mtx, &
         hubbard_projector_update
    use function_ops, only: function_ops_brappd_ketppd
    ! agrecokpt: needed for kpt argument in projectors_func_ovlp_box
    use geometry, only: POINT
    use integrals, only: integrals_kinetic
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use projectors, only: PROJECTOR_SET, projectors_func_ovlp_box
    use pseudopotentials, only: pseudopotentials_nonlocal_mat, PSEUDO_SPECIES, &
         pseudo_get_dij
    use rundat, only: pub_any_nl_proj, pub_hubbard, &
         pub_task, pub_hubbard_restart, pub_hubbard_atomsolve, &
         pub_hub_on_the_fly, pub_paw, pub_maxit_hotelling, &
         pub_max_resid_hotelling, pub_output_detail, pub_usp, pub_aug, &
         pub_realspace_projectors, pub_eda_scfmi, pub_build_bo, &
         pub_kpoint_method, pub_check_hermitian, pub_imag_thr, pub_debug, &
         pub_dftb_overlap_analytical
    use sparse, only: sparse_check_hermitian, sparse_rms_element
    use sparse_embed, only: sparse_embed_create, sparse_embed_destroy, SPAM3_EMBED, &
         sparse_embed_trace, sparse_embed_hotelling_invert, &
         sparse_embed_hotelling_init, sparse_embed_check_hermitian
    use utils, only: utils_abort, utils_assert

    implicit none

    ! Arguments
    type(NGWF_REP), intent(inout) :: rep
    type(FUNC_BASIS), intent(inout) :: ngwf_basis(rep%nsub)
    type(FUNC_BASIS), intent(in) :: projector_basis(rep%nsub)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(rep%nsub)
    type(FUNC_BASIS), intent(inout) :: hub_proj_basis(rep%nsub)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(MODEL), intent(in) :: mdl
    type(NGWF_REP), optional, intent(in) :: val_rep
    type(FUNC_BASIS), optional, intent(in) :: val_ngwf_basis(rep%nsub)
    real(kind=DP), optional, intent(in) :: kpt(3)
    logical, optional, intent(in) :: allow_pseuinvS ! ja531-> allow inverse S to be a pseudoinverse

    ! Local Variables
    type(DEM) :: inv_overlap_dens
    real(kind=DP) :: final_max_resid
    logical :: loc_allow_pseuinvS
    ! agrecokpt
    type(POINT) :: loc_kpt
    ! agrecocmplx
    logical :: loc_cmplx
    integer :: isub,jsub
    integer :: mrows, ncols
    ! rc2013: nonlocal projectors matrix
    type(SPAM3_EMBED) :: dij
    real(kind=DP) :: rms_hermitian_check

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering hamiltonian_dens_indep_matrices'

    ! agrecokpt: currently only KP method is being implemented
    ! now both methods are implemented, even though some routines
    ! have not been checked yet... but warnings/error messages are
    ! in place!
    !call utils_assert(pub_kpoint_method == 'KP', &
    !     'Subroutine hamiltonian_dens_indep_matrices currently supports&
    !     & only KP method for BZ sampling')

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(1)%iscmplx

    ! agrecokpt: default is Gamma point
    if (present(kpt)) then
       loc_kpt%x = kpt(1)
       loc_kpt%y = kpt(2)
       loc_kpt%z = kpt(3)
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    mrows = rep%overlap%mrows
    ncols = rep%overlap%ncols

    ! rc2013: if we need to block-orthogonalise the overlap do this first
    if(pub_build_bo) call internal_build_bo_matrices

    ! rc2013: begin regions loop
    do isub=1,ncols
        do jsub=1,mrows
            ! rc2013: build each matrix block separately
            call internal_build_matrices(jsub,isub)
        end do
    end do

    ! rc2013: calculate the nonlocal pseudopotentials matrix
    call internal_build_nonlocal

    ! rc2013: initialise inverse overlap
    call internal_build_inverse(rep%inv_overlap, rep%overlap, &
         rep%inv_overlap_init)

    ! jcap: calculate SCF-MI overlap matrices (if required)
    call internal_build_scfmi

    ! rc2013: diagnostics
    call internal_matrix_check

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hamiltonian_dens_indep_matrices'

    return

contains

    subroutine internal_build_matrices(jsub,isub)
      ! rc2013: calculate all the required matrix components

      implicit none

      ! Arguments
      integer, intent(in) :: isub, jsub ! rc2013: reguin indices
      real(kind=DP) :: tmp_val

      ! agrecocmplx
      loc_cmplx = rep%ngwfs_on_grid(isub)%iscmplx


      ! %%%%%%%%%%%%%%%%%%%%%%%%% OVERLAP MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Starting overlap matrix'

      if (pub_dftb_overlap_analytical) then
         call dftb_analytical_gaussian_overlap(rep%overlap%m(jsub,isub), &
              mdl, ngwf_basis(jsub))
      else
         call function_ops_brappd_ketppd(rep%overlap%m(jsub,isub), &
              rep%ngwfs_on_grid(jsub), ngwf_basis(jsub), rep%ngwfs_on_grid(isub), &
              ngwf_basis(isub), mdl%cell)
      end if

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Finished overlap matrix'

      ! %%%%%%%%%%%%%%%%%%%%%%% END OVERLAP MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ! ========================= KINETIC MATRIX ===============================
      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Starting kinetic matrix'

      call integrals_kinetic(rep%kinet%m(jsub,isub), rep%ngwfs_on_grid(jsub), &
          ngwf_basis(jsub), rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
          mdl%cell, mdl%fftbox)

      ! agrecokpt: add extra terms only if a non Gamma point is specified
      if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP) then
         ! rc2013: abort if trying to do this with embedding
         call utils_assert(mdl%nsub == 1, &
              'Error in hamiltonian_dens_indep_matrices: non-Gamma &
              &k-points not allowed with more than 1 subsystem.')
         select case(pub_kpoint_method)
            case('KP')
               if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
                  &Adding KP terms in kinetic matrix'

               ! agrecokpt: should we take into account the augmented overlap
               ! matrix in this case, when using PAW?
               call hamiltonian_update_kinet_kp(rep%kinet%p, &
                    rep%overlap%p, rep%ngwfs_on_grid(1), ngwf_basis(1), &
                    mdl, loc_kpt)

            case('TB')
               if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
                  &Adding TB terms in kinetic matrix'

               call hamiltonian_apply_phases_tb_method(rep%kinet%p, &
                    ngwf_basis(1), mdl, loc_kpt)

            case default
              call utils_abort('Illegal k-point method specified')

         end select
      end if

      ! agrecocmplx: check kinetic matrix is hermitian if requested
      if (pub_check_hermitian) then
         rms_hermitian_check = sparse_embed_check_hermitian(rep%kinet)

         if (pub_on_root) write(stdout,'(a,f8.5)') 'RMS kinetic hermitian check: ', &
            rms_hermitian_check

         call utils_assert(rms_hermitian_check < pub_imag_thr, &
              'WARNING: hermitian check for kinetic matrix above threshold')

      end if

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Finished kinetic matrix'
      !! ======================= END KINETIC MATRIX =============================

      ! ================= NGWF-PROJECTOR OVERLAP MATRIX =====================
      ! calculate the ngwf-projector overlap matrix
      if (pub_any_nl_proj) then
         if (pub_debug_on_root) write(stdout,'(a)') &
              'DEBUG: Starting NGWF-Proj overlap matrix'

         if (.not.pub_realspace_projectors) then
            ! agrecokpt: use specified kpoint, but only in KP method; in TB
            ! method the whole hamiltonian is multiplied by k-point dependent
            ! terms, so no need to have k-point dependence at this stage
            select case (pub_kpoint_method)
               ! KP method
               case('KP')
                  call projectors_func_ovlp_box(rep%sp_overlap%m(jsub,isub), &
                       rep%ngwfs_on_grid(jsub),ngwf_basis(jsub),projector_basis(isub), &
                       nl_projectors(isub), mdl%fftbox,mdl%cell, kshift=loc_kpt)
               ! TB method
               case('TB')
                  call projectors_func_ovlp_box(rep%sp_overlap%m(jsub,isub), &
                       rep%ngwfs_on_grid(jsub),ngwf_basis(jsub),projector_basis(isub), &
                       nl_projectors(isub), mdl%fftbox,mdl%cell)
               case default
                  call utils_abort('Illegal k-point method specified')
            end select
         else
            ! agrecokpt: k-point dependence already included in the
            ! definition of projs_on_grid
            call function_ops_brappd_ketppd(rep%sp_overlap%m(jsub,isub), &
                 rep%ngwfs_on_grid(jsub), ngwf_basis(jsub), &
                 nl_projectors(isub)%projs_on_grid, projector_basis(isub), mdl%cell)
         end if

         ! rc2013: nonlocal matrix extracted to internal_build_nonlocal

      else if (pub_paw) then
         if (pub_debug_on_root) write(stdout,'(a)') &
              'DEBUG: Starting sp overlap matrix'
         ! ndmh: Calculate the overlap matrix of the NGWFs and PAW projectors
         ! agrecokpt: kpoint dependence added, but only in KP method; in TB
         ! method, k-point dependent terms are added to whole hamiltonian,
         ! no need to include them in the SP overlap now
         select case (pub_kpoint_method)
            ! KP method
            case('KP')
               call projectors_func_ovlp_box(rep%sp_overlap%m(jsub,isub), &
                    rep%ngwfs_on_grid(jsub),ngwf_basis(jsub),projector_basis(isub), &
                    nl_projectors(isub), mdl%fftbox,mdl%cell, kshift=loc_kpt)
            ! TB method
            case('TB')
               call projectors_func_ovlp_box(rep%sp_overlap%m(jsub,isub), &
                    rep%ngwfs_on_grid(jsub),ngwf_basis(jsub),projector_basis(isub), &
                    nl_projectors(isub), mdl%fftbox,mdl%cell)
            case default
               call utils_abort('Illegal k-point method specified')
         end select

         ! ndmh: NB - the nonlocal matrix cannot be calculated at this stage in
         ! ndmh: PAW, since it depends on the density kernel!
         if (pub_debug_on_root) write(stdout,'(a)') &
              'DEBUG: Finished sp overlap matrix'
      end if

    ! ================= END NGWF-PROJECTOR OVERLAP MATRIX =====================

    ! %%%%%%%%%%%%%%%%%%%%%%%%% CROSS OVERLAP MATRIX %%%%%%%%%%%%%%%%%%%%%%%%
    ! lr408: calculate cross overlap matrix if necessary
    if (present(val_rep)) then

       if (pub_debug_on_root) write(stdout,'(a)') &
            'DEBUG: Starting cross overlap matrix'

       ! agrecokpt: not sure here what to do,
       ! my guess is that we would store a different kernel
       ! and NGWF set for each k-point, then compute the cross
       ! overlap by explicitly looking at the previosuly converged
       ! kernel/NGWF set for that k-point; not sure if the k-point
       ! dependence requires inclusion of extra factors in the
       ! matrix elements
       if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP) then
          if (pub_debug_on_root) write(stdout,'(a)') &
             'WARNING: computation of cross overlap matrix in &
             &hamiltonian_dens_indep_matrices not yet checked for &
             &non-Gamma k-point'
       end if

       call function_ops_brappd_ketppd(rep%cross_overlap%m(jsub,isub), &
            val_rep%ngwfs_on_grid(jsub), val_ngwf_basis(jsub), &
            rep%ngwfs_on_grid(isub), ngwf_basis(isub), mdl%cell)

       if (pub_aug) then
          ! ndmh: Calculate the augmentation of the cross overlap matrix due to
          ! ndmh: the augmentation region part of the overlap operator
          call augmentation_overlap(rep%cross_overlap%m(jsub,isub), &
               mdl%regions(jsub)%pseudo_sp, mdl%regions(isub)%paw_sp, &
               val_rep%sp_overlap%m(jsub,isub), rep%sp_overlap%m(jsub,isub))
       end if

       ! agrecocmplx: check cross overlap matrix is hermitian if requested
       ! rc2013: this is not set up with embedding yet...should be outside loop
       if (pub_check_hermitian) then
          rms_hermitian_check = &
               sparse_check_hermitian(rep%cross_overlap%m(jsub,isub))

          if (pub_on_root) write(stdout,'(a,f8.5)') 'RMS cross overlap&
             & hermitian check: ', rms_hermitian_check

          call utils_assert(rms_hermitian_check < pub_imag_thr, &
               'WARNING: hermitian check for cross overlap above threshold')

       end if

       if (pub_debug_on_root) write(stdout,'(a)') &
            'DEBUG: Finished cross overlap matrix'

    end if

    ! %%%%%%%%%%%%%%%%%%%%%%% END CROSS OVERLAP MATRIX %%%%%%%%%%%%%%%%%%%%%%



    ! ====================== DFT+U CORRECTION MATRIX =========================

    ! ddor: DFT+U calculation
    if (pub_hubbard) then

       ! ddor: DFT+U calculation with self-consistency over Hubbard projectors
       !     : (either HUBBARDSCF or on-the fly), or just a fixed-projector
       !     : DFT+U calculation using projectors in NGWF format.
       consist: if ( ((pub_task == 'HUBBARDSCF') .and.  &
            & (hub%consistency_iteration > 1)) .or. &
            & pub_hubbard_restart .or. pub_hubbard_atomsolve .or. &
            & pub_hub_on_the_fly) then

          if (pub_debug_on_root) write(stdout,'(a,i4)') &
               'DEBUG: Hubbard correction matrix on HUBBARDSCF iteration', &
               hub%consistency_iteration

          ! ddor: If carrying out on-the-fly HUBBARDSCF, renew the projectors
          !     : here each time the NGWFs change. This affects NGWF CG etc.
          if (pub_hub_on_the_fly) then
             ! agrecokpt: kpoint dependence needs to be checked in this routine
             ! in KP method, no need for k-point dependence in TB method
             select case(pub_kpoint_method)
                ! KP method
                case('KP')
                   call hubbard_projector_update( &
                        rep%ngwfs_on_grid(1), ngwf_basis(1), &
                        nl_projectors(1), projector_basis(1), &
                        hub_proj_basis(1), hub, &
                        rep%inv_overlap%p, rep%hub_overlap%p, &
                        rep%hub_overlap_t%p, &
                        rep%sp_overlap%p, rep%hub_proj_paw_overlap%p, &
                        rep%hub_ngwf_paw_overlap%p, kpt=loc_kpt)
                case('TB')
                   call hubbard_projector_update( &
                        rep%ngwfs_on_grid(1), ngwf_basis(1), &
                        nl_projectors(1), projector_basis(1), &
                        hub_proj_basis(1), hub, &
                        rep%inv_overlap%p, rep%hub_overlap%p, &
                        rep%hub_overlap_t%p, &
                        rep%sp_overlap%p, rep%hub_proj_paw_overlap%p, &
                        rep%hub_ngwf_paw_overlap%p)
                case default
                   call utils_abort('Illegal k-point method specified')
             end select

          endif

          ! ddor: Using the NGWFs from the previous HUBBARDSCF
          !     : iteration calculate the new hub_overlap.
          call function_ops_brappd_ketppd(rep%hub_overlap%p, &    ! input-output
               rep%ngwfs_on_grid(1), ngwf_basis(1), hub%consistency_projs, & !input
               hub_proj_basis(1), mdl%cell)                               !input

          ! ddor: Augment the <NGWF|proj> matrix for PAW+U or USP calculations
          if (pub_aug) then

             if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP) then
                if (pub_debug_on_root) write(stdout,'(a)') &
                   'WARNING: computation of PAW-U overlap matrix in &
                   &hamiltonian_dens_indep_matrices not yet checked for &
                   &non-Gamma k-point'
             end if

             ! ddor: Calculate the overlap matrix of the
             ! ddor: Hubbard and PAW projectors
             ! agrecokpt: not sure if kpt dependence is needed here
             call projectors_func_ovlp_box(rep%hub_proj_paw_overlap%p, &
                  hub%consistency_projs,hub_proj_basis(1),&
                  projector_basis(1),nl_projectors(1),mdl%fftbox,mdl%cell)!,kpt=loc_kpt

             ! ddor: Calculate the augmentation of the
             ! ddor: overlap matrix due to the
             ! ddor: augmentation region part of the overlap operator
             call augmentation_overlap(rep%hub_overlap%p,&
                  mdl%regions(1)%pseudo_sp,mdl%regions(1)%paw_sp, &
                  rep%sp_overlap%p,rep%hub_proj_paw_overlap%p)

          end if

          ! calculate the Hubbard projection operator for each site,
          ! using hub_overlap and the tensorial correction matrix
          call hubbard_projection_mtx(.true.,.false., rep%hub_overlap%p, &
               rep%hub_overlap_t%p, hub%o_matrix)

          ! ddor: DFT+U calculation with projectors in non-NGWF format
       elseif ( ((.not. pub_hubbard_restart) .and. &  ! This is a must
            & (.not. pub_hubbard_atomsolve) .and. &
            & (.not. pub_hub_on_the_fly)) .and. &
            & ((pub_task .ne. 'HUBBARDSCF') .or. &        ! Then either standard
            & ((pub_task == 'HUBBARDSCF') .and. &         ! or HUBBARDSCF first iter
            & (hub%consistency_iteration .eq. 1) .and. &
            & hub%dftu_on_first_hubbardscf )) ) then

          if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Starting &
               &hubbard_sp_ovlp_box'

          if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP) then
             if (pub_debug_on_root) write(stdout,'(a)') &
                'WARNING: computation of Hubbard-SP overlap matrix in &
                &hamiltonian_dens_indep_matrices not yet checked for &
                &non-Gamma k-point'
          end if

          ! ddor: Calculate the Hubbard on-site <ngwf|atomic> overlap matrix
          ! agrecokpt: not sure if kpt dependence is needed here
          call projectors_func_ovlp_box(rep%hub_overlap%p, &             ! out
               rep%ngwfs_on_grid(1),ngwf_basis(1),hub_proj_basis(1),&
               hub%projectors,mdl%fftbox,mdl%cell)!,kpt=loc_kpt

          if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Finished &
               &hubbard_sp_ovlp_box'

          ! ddor: calculate the Hubbard projection operator for each site, with
          ! ddor: no tensorial correction (.false.) since we use atomic
          ! ddor: projectors
          call hubbard_projection_mtx(.false.,.false., rep%hub_overlap%p, &
               rep%hub_overlap_t%p, hub%o_matrix)

       endif consist

       ! agrecocmplx: check hubbard overlap matrix is hermitian if requested
       if (pub_check_hermitian) then
          rms_hermitian_check = sparse_check_hermitian(rep%hub_overlap%p)

          if (pub_on_root) write(stdout,'(a,f8.5)') 'RMS hubbard overlap&
             & hermitian check: ', rms_hermitian_check

          call utils_assert(rms_hermitian_check < pub_imag_thr, &
               'WARNING: hermitian check for hubbard overlap matrix above threshold')

       end if

       if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Finished Hubbard &
            &correction matrix'

    endif ! DFT+U

    end subroutine internal_build_matrices

    ! ================= NON-LOCAL PSEUDOPOTENTIAL MATRIX =====================

    ! rc2013: a separate subroutine to set up matrices that can't be done on
    ! a "block-by-block" basis
    subroutine internal_build_nonlocal

      implicit none

      if(pub_any_nl_proj) then
         if (.not.pub_usp) then
            if (pub_debug_on_root) write(stdout,'(a)') &
                'DEBUG: Starting nonlocal matrix'

            ! rc2013: create dij to hold nonlocal terms
            dij%structure = 'E'
            call sparse_embed_create(dij)
            ! rc2013: build diagonal terms of dij -- this allows us
            ! to avoid passing all pseudo_sp's for nonlocpot
            do isub=1,mdl%nsub
               call pseudo_get_dij(dij%m(isub,isub), mdl%regions(isub)%pseudo_sp)
            enddo

            call pseudopotentials_nonlocal_mat(rep%nonlocpot, &
                 rep%sp_overlap, dij)
            ! rc2013: clean-up
            call sparse_embed_destroy(dij)

            if (pub_debug_on_root) write(stdout,'(a)') &
                 'DEBUG: Finished nonlocal matrix'

            if ((pub_kpoint_method == 'TB') .and. &
                 (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP)) then

               if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
                    &Adding TB terms in nonlocpot matrix'

               call hamiltonian_apply_phases_tb_method(rep%kinet%p, &
                    ngwf_basis(1), mdl, loc_kpt)

            end if

            ! agrecocmplx: check nonlocal matrix is hermitian if requested
            if (pub_check_hermitian) then
               rms_hermitian_check = sparse_embed_check_hermitian(rep%nonlocpot)

               if (pub_on_root) write(stdout,'(a,f8.5)') 'RMS nonlocal hermitian&
                  & check: ', rms_hermitian_check

               call utils_assert(rms_hermitian_check < pub_imag_thr, &
                    'WARNING: hermitian check for nonlocal matrix above threshold')

            end if
        end if
      end if

      if (pub_aug) then
         ! ndmh: Calculate the augmentation of the overlap matrix due to the
         ! ndmh: augmentation region part of the overlap operator
         ! agrecokpt: in KP method, the k-point dependence is already included
         ! in the sp_overlap due to the k-point dependent projectors;
         ! do we need to use the augmented overlap when we add the extra term
         ! in the kinetic energy?
         ! in TB method, augment the Gamma point overlap matrix, and only then
         ! add the extra k-point dependent terms
         ! rc2013: this has not been tested with embedding structures yet
         do jsub=1,mdl%nsub
            do isub=1,mdl%nsub
               call augmentation_overlap(rep%overlap%m(isub,jsub), &
                    mdl%regions(isub)%pseudo_sp,mdl%regions(jsub)%paw_sp, &
                    rep%sp_overlap%m(isub,jsub))
            end do
         end do

      end if

      ! agrecokpt: in TB method, need to include appropriate extra terms in the
      ! overlap matrix, only if a non Gamma point is specified; the whole overlap
      ! matrix is multiplied by the appropriate extra-terms, regardless if it is
      ! the standard or the augmented one
      if ((pub_kpoint_method == 'TB') .and. &
         (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP)) then
         if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
            &Adding TB terms in overlap matrix'

         ! rc2013: abort if trying to do this with embedding
         call utils_assert(mdl%nsub == 1, &
              'Error in hamiltonian_dens_indep_matrices: non-Gamma &
              &k-points not allowed with more than 1 subsystem.')
         call hamiltonian_apply_phases_tb_method(rep%overlap%p, &
              ngwf_basis(1), mdl, loc_kpt)

      end if

      ! agrecocmplx: check final overlap matrix is hermitian if requested
      if (pub_check_hermitian) then
         rms_hermitian_check = sparse_embed_check_hermitian(rep%overlap)

         if (pub_on_root) write(stdout,'(a,f8.5)') 'RMS overlap hermitian check: ', &
            rms_hermitian_check

         call utils_assert(rms_hermitian_check < pub_imag_thr, &
              'WARNING: hermitian check for overlap matrix above threshold')

      end if

    end subroutine internal_build_nonlocal

    ! =============== END NON-LOCAL PSEUDOPOTENTIAL MATRIX ===================

    subroutine internal_build_inverse(inv_overlap, overlap, inv_overlap_init, &
         ireg)

      use sparse, only: sparse_create, sparse_product, sparse_destroy, &
           sparse_scale, sparse_rms_element, sparse_max_abs_element, &
           sparse_num_element, SPAM3

      implicit none

      ! Arguments
      type(SPAM3_EMBED), intent(inout) :: inv_overlap
      type(SPAM3_EMBED), intent(in)    :: overlap
      logical, intent(inout)           :: inv_overlap_init
      integer, intent(in), optional    :: ireg

      ! Local variables
      integer :: isub, jsub
      real(kind=DP) :: max_resid  ! maximum value of residual
      real(kind=DP) :: frob_norm
      type(SPAM3)   :: s_sinv

      ! ======================== INVERSE OVERLAP MATRIX ========================

      if (.not.inv_overlap_init) then
          ! cks: initialise inverse overlap guess if Hotelling recursion
          ! cks: will be used to approximate it.
          ! rc2013: call sparse_embed_mod version of hotelling's algorithm initialisation
          if (pub_maxit_hotelling > 0) &
              call sparse_embed_hotelling_init(inv_overlap, overlap)
          inv_overlap_init = .true.
      end if


      ! cks : approximate inverse overlap by Hotelling's recursion
      if (pub_maxit_hotelling > 0) then

        if (pub_output_detail>=VERBOSE .and. pub_on_root) then
            write(stdout,'(a)')'============ Calculation of NGWF S^-1 using Hotelling algorithm ================ '
            write(stdout,'(a)')'         Iteration  |  Frobenius norm and Maximum &
                &abs value    (of I-S*S_n^-1)   '
        end if

        call sparse_embed_hotelling_invert(inv_overlap, & !inout
             overlap,show_output=(pub_output_detail>=VERBOSE), &
             max_resid_converged=pub_max_resid_hotelling, &
             num_iter=pub_maxit_hotelling,final_max_resid=final_max_resid, &
             ireg=ireg)


        ! ndmh: test if Hotelling failed to invert the matrix
        if (final_max_resid > 0.95_DP) then

           if (pub_output_detail>=VERBOSE.and.pub_on_root) then
              write(stdout,'(15x,a)')'Resetting Hotelling algorithm'
           end if

           ! ndmh: re-initialise Hotelling algorithm from scratch
           call sparse_embed_hotelling_init(inv_overlap, overlap)

           call sparse_embed_hotelling_invert(inv_overlap, & !inout
               overlap,show_output=(pub_output_detail>=VERBOSE), &
               max_resid_converged=pub_max_resid_hotelling, &
               num_iter=pub_maxit_hotelling,final_max_resid=final_max_resid)


           loc_allow_pseuinvS=.false.
           if(present(allow_pseuinvS)) then
              loc_allow_pseuinvS=allow_pseuinvS
           end if
           if(.not.loc_allow_pseuinvS) then
              ! ndmh: check if it is still broken
              if (final_max_resid > 0.95_DP) then
                 call utils_abort('Error in hamiltonian_dens_indep_matrices: &
                      &Inversion of overlap matrix failed')
              end if
           end if

        end if

       if (pub_output_detail>=VERBOSE.and.pub_on_root) then
          write(stdout,'(a)')'===================================&
               &============================================='
       end if

    else if (pub_maxit_hotelling == 0) then
       call utils_assert(mdl%nsub .eq. 1, &
            'Error in hamiltonian_dens_indep_matrices: &
            &maxit_hotelling = 0 but dense inversion is not &
            &implemtented/tested with more than 1 subsystem.')

       ! ndmh: allocate storage for dense inverse
       ! agrecocmplx: inv_overlap_dens is complex if
       ! inv_overlap is complex
       call dense_create(inv_overlap_dens,ngwf_basis(1)%num,ngwf_basis(1)%num, &
                         iscmplx=loc_cmplx)

       ! ndmh: copy current overlap matrix into dense inverse overlap
       ! agrecokpt: in TB method, k-point dependence already included
       ! in overlap matrix
       call dense_convert(inv_overlap_dens, overlap)

       ! ndmh: invert dense matrix
       call dense_invert(inv_overlap_dens)

       ! ndmh: copy dense inverse overlap back to sparse matrix
       call dense_convert(inv_overlap,inv_overlap_dens)

       ! ndmh: deallocate storage for dense inverse
       call dense_destroy(inv_overlap_dens)
    else
       call utils_abort('Error in hamiltonian_dens_indep_matrices: &
            &negative pub_maxit_hotelling ',pub_maxit_hotelling)
    end if

    ! agrecocmplx: check inverse overlap matrix is hermitian if requested
    if (pub_check_hermitian) then
       rms_hermitian_check = sparse_embed_check_hermitian(inv_overlap)

       if (pub_on_root) write(stdout,'(a,f8.5)') 'RMS inverse overlap&
          & hermitian check: ', rms_hermitian_check

       call utils_assert(rms_hermitian_check < pub_imag_thr, &
            'WARNING: hermitian check for inverse overlap matrix above threshold')

    end if

    ! ==================== END INVERSE OVERLAP MATRIX ========================

    end subroutine internal_build_inverse

    ! jcap: calculate the SCF-MI overlap matrices if required
    subroutine internal_build_scfmi

    implicit none

    ! ======================== SCF-MI OVERLAP MATRICES =======================

    ! mjsp: if SCF-MI:
    if (pub_eda_scfmi) then

       if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP) then
          if (pub_debug_on_root) write(stdout,'(a)') &
             'WARNING: computation of fragment and supermolecule overlap &
             &matrix in hamiltonian_dens_indep_matrices not yet checked for &
             &non-Gamma k-point'
       end if

       ! mjsp: construct fragment and supermolecule overlap and inverse
       ! mjsp: overlap matrices from rep%overlap and rep%inv_overlap
       ! mjsp: (this involves rep%overlap and rep%inv_overlap
       ! mjsp: being converted from the full supermolecule to block-diagonal
       ! mjsp: representation)
       ! agrecokpt: not checked for k-points yet, I assume dependence in
       ! TB method propagates from overlap and inv_overlap?
       call fmo_construct_overlap_mats(rep, ngwf_basis)

       ! mjsp: backup fragment overlap in dense form for later use
       ! mjsp: in SCF-MI module
       call dense_convert(overlap_f_dens,rep%overlap_scfmi_full)

    end if

    ! ==================== END SCF-MI OVERLAP MATRICES =======================

    end subroutine internal_build_scfmi

    ! %%%%%%%%%%%%%%%%%%%%%% BLOCK-ORTHOGONALISED MATRICES %%%%%%%%%%%%%%%%%%%%%

    ! rc2013: construct block-orthogonalised matrices and NGWFs
    subroutine internal_build_bo_matrices

      use datatypes, only: FUNCTIONS, data_set_to_zero, data_functions_axpy,  &
           data_functions_alloc, data_functions_dealloc
      use function_ops, only: function_ops_sum_ppd_funcs
      use rundat, only: pub_emft, pub_active_region
      use sparse, only: sparse_copy, SPAM3, sparse_create, sparse_destroy
      use sparse_embed, only: sparse_embed_scale, sparse_embed_product, &
           sparse_embed_extract_sub

      implicit none

      ! Local variables
      type(SPAM3_EMBED) :: reg_overlap(rep%nsub), reg_inv_overlap(rep%nsub)
      logical :: inv_overlap_init
      type(FUNCTIONS) :: proj_ngwfs_on_grid(1)
      type(SPAM3) :: proj_matrix(1)
      integer :: ireg

      ! agrecocmplx
      loc_cmplx = rep%ngwfs_on_grid(1)%iscmplx

      ! rc2013: construct non-orthogonalised NGWF overlap
      call sparse_embed_scale(rep%inv_overlap, 0.0_DP)

      if (pub_output_detail>=VERBOSE .and. pub_on_root) &
           write(stdout,'(a)') 'Block orthogonalising NGWFs...'
      do isub=1,mrows

         inv_overlap_init = .false.
         if (pub_output_detail>=VERBOSE .and. pub_on_root) &
              write(stdout,'(a,i1)') '... constructing inverse overlap &
              &for subsystem ', isub
         do jsub=1,ncols
            call function_ops_brappd_ketppd(rep%overlap%m(jsub,isub), &
                 rep%ngwfs_on_grid(jsub), ngwf_basis(jsub), &
                 rep%ngwfs_on_grid(isub), ngwf_basis(isub), mdl%cell)
         end do

         ! rc2013: now construct the diagonal matrices
         call sparse_embed_create(reg_overlap(isub),rep%overlap, &
              arow=isub,bcol=isub)
         call sparse_embed_create(reg_inv_overlap(isub),rep%inv_overlap, &
              arow=isub,bcol=isub)
         call sparse_embed_extract_sub(reg_overlap(isub),rep%overlap,isub,isub)
         ! rc2013: for testing, always reset inv_overlap -- probably unnecessary
         call internal_build_inverse(reg_inv_overlap(isub), reg_overlap(isub), &
              inv_overlap_init, isub)
         ! Copy the diagonal inverse overlap's back to full matrix
         call sparse_copy(rep%inv_overlap%m(isub,isub), reg_inv_overlap(isub)%m(1,1))
      end do

      ! rc2013: now let's build the projector matrices
      call sparse_embed_product(rep%bo_projector, rep%inv_overlap, rep%overlap)
      call sparse_embed_product(rep%bo_projector_t, rep%overlap, rep%inv_overlap)

      ireg = pub_active_region
      do isub=1,ncols
         if(isub .ne. ireg) then
            call data_functions_alloc(proj_ngwfs_on_grid(1), &
                 ngwf_basis(isub)%size_on_grid, iscmplx=rep%ngwfs_on_grid(isub)%iscmplx)
            call sparse_create(proj_matrix(1), rep%bo_projector%m(ireg,isub))

            ! Initialise
            call data_set_to_zero(proj_ngwfs_on_grid(1))
            call sparse_copy(proj_matrix(1), rep%bo_projector%m(ireg,isub))

            ! rc2013: calculate |\chi_A>P^A_B
            call function_ops_sum_ppd_funcs(proj_ngwfs_on_grid(1:1), &   ! inout
                 ngwf_basis(isub), &!proj_matrix(1:1), 1, 1, &
                 rep%bo_projector%m(ireg,isub), 1, 1, &           ! input
                 rep%overlap%m(ireg,isub), &
                 rep%ngwfs_on_grid(ireg), ngwf_basis(ireg))

            ! rc2013: finally add the projector term to the NGWFs
            call data_functions_axpy(rep%ngwfs_on_grid(isub), &
                 proj_ngwfs_on_grid(1), -1.0_DP)

            ! rc2013: clean-up
            call data_functions_dealloc(proj_ngwfs_on_grid(1))
            call sparse_destroy(proj_matrix(1))
         end if
      end do
      if (pub_output_detail>=VERBOSE .and. pub_on_root) &
           write(stdout,'(a)') '... block orthogonalisation done.'

      ! %%%%%%%%%%%%%%%%%%% END BLOCK-ORTHOGONALISED MATRICES %%%%%%%%%%%%%%%%%%

    end subroutine internal_build_bo_matrices


    ! rc2013: print matrix traces as check
    subroutine internal_matrix_check

      use rundat, only: pub_debug_on_root, pub_debug
      use sparse_embed, only: sparse_embed_trace
      use sparse, only: sparse_trace

      implicit none

      ! Local variable
      real(kind=DP) :: debug_trace

      if(pub_debug) then
          call sparse_embed_trace(debug_trace, rep%kinet, rep%kinet)
          if(pub_on_root) write(stdout,'(a,f24.14)') &
                'DEBUG: tr[Kin*Kin]:', debug_trace
          call sparse_embed_trace(debug_trace, rep%overlap, rep%overlap)
          if(pub_on_root) write(stdout,'(a,f24.14)') &
                'DEBUG: tr[S*S]:', debug_trace
          call sparse_embed_trace(debug_trace, rep%overlap, rep%inv_overlap)
          if(pub_on_root) write(stdout,'(a,f24.14)') &
                'DEBUG: tr[S*S^-1]:', debug_trace
          call sparse_embed_trace(debug_trace, rep%inv_overlap, rep%inv_overlap)
          if(pub_on_root) write(stdout,'(a,f24.14)') &
                'DEBUG: tr[S^-1*S^-1]:', debug_trace
      end if

    end subroutine internal_matrix_check

  end subroutine hamiltonian_dens_indep_matrices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hamiltonian_dens_dep_matrices(ham, lhxc_fine, total_energy, &
       lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
       ngwf_basis, hub_proj_basis, hub, denskern, mdl, hfxstate, &
       ham_update, lhxc_fixed, spoil_force, dfdtau_fine, kpt)

    !===================================================================!
    ! This subroutine calculates density-dependent contributions to the !
    ! energy and components of the Hamiltonian matrix.                  !
    ! Updating of the Hamiltonian matrix components is optional, and is !
    ! controlled by the logical flag ham_update                         !
    !-------------------------------------------------------------------!
    ! Written by Nicholas Hine on 22/04/2010, assembled from bits of    !
    ! code formerly in electronic_mod.                                  !
    !===================================================================!

    use augmentation, only: aug_projector_denskern, aug_nonlocal_mat
    use comms, only: pub_on_root
    use constants, only: DP, paw_en_size, paw_en_etxc, &
         paw_en_ehart, paw_en_exc, paw_en_dijhat, paw_en_dijxc, &
         paw_en_exc_core, stdout
#ifdef GPU_SP_TEST
    use fourier_gpu_wrapper_mod, only: deltaE_gpu, E_gpu_previous
#endif
    use function_basis, only: FUNC_BASIS
    ! agrecokpt
    use geometry, only: POINT
    use hf_exchange, only: HFX_STATE, hf_exchange_calculate
    use hubbard_build, only: HUBBARD_MODEL, hubbard_ham_matrix, &
         hubbard_energy_total, hubbard_projector_ham
    use integrals, only: integrals_locpot, integrals_div_locpot_grad_functions, &
         integrals_confinement_pot, &
         test_integrals_div_locpot_grad_functions
    use is_smeared_ions, only: smeared_ion_E_self, smeared_ion_E_smeared
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_HAM, NGWF_REP
    use polarisable_embedding, only: polarisable_embedding_calculate
    use rundat, only: pub_use_hfx, pub_hubbard, pub_paw, pub_any_nl_proj, &
         pub_is_smeared_ion_rep, pub_aug, pub_usp, pub_num_spins, &
         pub_spin_fac,pub_confined_ngwfs, pub_pol_emb_pot, &
         pub_num_kpoints, PUB_1K, pub_dmft_points, &
         pub_dmft_fully_sc, pub_dmft_fully_sc_h, pub_imag_thr, &
         pub_debug, pub_devel_code, & ! JCW
         pub_xc_ke_density_required, pub_inner_loop_iteration, &
         pub_kpoint_method, pub_check_hermitian, pub_rootname, &
         pub_use_activehfx, pub_active_region, pub_emft, pub_emft_follow, & ! jcap
         pub_pol_emb_vacuum_qmstar, pub_scissor_ngroups
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_trace, &
        sparse_product, sparse_copy, sparse_take_real_part
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
         sparse_embed_create, sparse_embed_destroy, sparse_embed_scale, &
         sparse_embed_trace, sparse_embed_product, sparse_embed_copy, &
         sparse_embed_array2mat, sparse_embed_array_scale, &
         sparse_embed_check_hermitian, sparse_embed_array_to_sparse_array, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use sparse_array, only: SPAM3_ARRAY, sparse_array_destroy
    use file_handling, only: file_handling_remove_file, file_handling_remove_numbered_files
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_devel_code, utils_int_to_str
    use vdwcorrection, only: pub_dispersion_energy
    use visual, only: visual_scalarfield_ngwf_basis
    use xc, only: xc_embed_swap_functional

    implicit none

    ! Arguments
    type(NGWF_REP), intent(in) :: rep
    type(NGWF_HAM), intent(inout) :: ham
    type(FUNC_BASIS), intent(in) :: ngwf_basis(rep%nsub)
    type(FUNC_BASIS), intent(in) :: hub_proj_basis(rep%nsub)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(SPAM3_EMBED_ARRAY), intent(inout)   :: denskern !could be intent(in) if rescalings avoided
    real(kind=DP), intent(out)  :: total_energy
    ! jd: It's important that lhxc_energy is inout. When lhxc_fixed is .true., this subroutine
    !     *DOES NOT* calculate lhxc_energy (see [*] below). In this case we want to reuse the
    !     previous value for lhxc_energy. This only happens in the call from ensemble_dft_mod:671.
    !     Previously this argument was intent(out), which led to it being used uninitialised.
    real(kind=DP), intent(inout)  :: lhxc_energy
    real(kind=DP), intent(out)  :: hubbard_energy
    real(kind=DP), intent(out)  :: paw_sphere_energies(paw_en_size)
    type(MODEL), intent(in)     :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    ! rc2013: each region can have its own lhxc potential
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, rep%nsub)
    ! JCW: Optional gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density. Required to evaluate lhxc matrix elements for meta-GGAs.
    real(kind=DP), optional, intent(inout) :: dfdtau_fine(:,:,:,:)
    ! JCW: In general, we have
    ! JCW:    dfdtau_fine(mdl%fine_grid%ld1,&
    ! JCW:    mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    ! JCW: or
    ! JCW:    dfdtau_fine(1,1,1,1)
    ! JCW: when it is unneeded.

    logical, intent(in) :: ham_update
    logical, intent(in) :: lhxc_fixed
    logical,optional :: spoil_force
    ! agrecokpt: kpoint needed for TB method; currently the extra
    ! terms are included to each term of the Hamiltonian (testing/debug),
    ! but if it is sufficient to add them to the whole Hamiltonian, do that
    ! instead
    real(kind=DP), optional, intent(in) :: kpt(3)

    ! Local Variables
    type(SPAM3_EMBED), allocatable   :: projector_denskern(:)
    real(kind=DP), allocatable :: density_fine(:,:,:,:,:,:)
    real(kind=DP) :: polemb_energy_mm_valence
    real(kind=DP) :: polemb_energy_mm_vdw
    real(kind=DP) :: polemb_energy_mm_es
    real(kind=DP) :: polemb_energy_qmstar_elec_mm_perm
    real(kind=DP) :: polemb_energy_qmstar_elec_mm_ind
    real(kind=DP) :: polemb_energy_qmfull_elec_mm_perm
    real(kind=DP) :: polemb_energy_qmfull_elec_mm_ind
    real(kind=DP) :: polemb_energy_qm_core_mm_perm
    real(kind=DP) :: polemb_energy_qm_core_mm_ind
    real(kind=DP) :: polemb_energy_qmmm_vdw
    real(kind=DP) :: polemb_energy_other
    integer, save :: lhxc_output_iter = 1
    logical, save :: lhxc_output_done = .false.
    integer :: is, imat
    integer :: ierr
    logical :: dmft_ham_update
    real(kind=DP), save        :: dmft_correction
    !real(kind=DP)       :: eigen_en(ngwf_basis%num,pub_num_spins)
    type(SPAM3_EMBED) :: dmft_kernel(pub_num_spins),dmft_self(pub_num_spins), &
         dmft_z(pub_num_spins),backupH(pub_num_spins)
    integer          :: godmft
    logical          :: loop_ham_update
    ! agrecokpt
    type(POINT) :: loc_kpt
    real(kind=DP) :: rms_hermitian_check
    integer :: mrows, ncols, isub, jsub, lhxc_index
    real(kind=DP)    :: ewald_energy
    real(kind=DP) :: trace  ! Matrix trace
    ! rc2013: horrible workaround for spin/embedding issue
    type(SPAM3), allocatable :: kern_array(:), aug_array(:), ham_array(:)
    type(SPAM3_ARRAY) :: temp_kern

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hamiltonian_dens_dep_matrices'

    ! jme: KPOINTS_DANGER
    ! Some parts of this function enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine hamiltonian_dens_dep_matrices not ready yet &
         &for more than one k-point.')

    ! agrecokpt: default is Gamma point
    if (present(kpt)) then
       loc_kpt%x = kpt(1)
       loc_kpt%y = kpt(2)
       loc_kpt%z = kpt(3)
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    if (.not.pub_xc_ke_density_required) then
       ! JCW: Non-tau-dependent XC functional
       ! JCW: Check dfdtau_fine is either not present, or has a length of 1
       if ( present(dfdtau_fine) ) then
          call utils_assert( &
               size(dfdtau_fine) == 1,&
               "Error in hamiltonian_dens_dep_matrices: &
               &pub_xc_ke_density_required is false, but dfdtau_fine is present &
               &or has a size > 1.")
       end if
    else
       ! JCW: tau-dependent XC functional
       ! JCW: Check dfdtau_fine is present
       call utils_assert( present(dfdtau_fine),&
            "Error in hamiltonian_dens_dep_matrices: &
            &pub_xc_ke_density_required is true, but dfdtau_fine is not &
            &present.")
    end if

    ! rc2013: set size of embedding matrices
    mrows = denskern%mrows
    ncols = denskern%ncols

    allocate(kern_array(pub_num_spins),stat=ierr)
    call utils_alloc_check('HDDM','kern_array',ierr)

    ! CW
    dmft_ham_update = (pub_dmft_points>0.and.pub_dmft_fully_sc)
    if (dmft_ham_update) then
       do is=1,pub_num_spins
          if (.not.ham_update) call sparse_embed_create(backupH(is), &
                ham%ham(is),iscmplx=.false.)
          ! agrecocmplx: use safe conversion from complex to real
          ! when ham is complex
          if (.not.ham_update) then
             if (ham%ham(is)%iscmplx) then
                ! rc2013: need embedding version
                do jsub=1,ncols
                   do isub=1,mrows
                      call sparse_take_real_part(backupH(is)%m(isub,jsub), &
                           ham%ham(is)%m(isub,jsub),pub_imag_thr)
                   end do
                end do
             ! use standard copy for real to real case
             else
                call sparse_embed_copy(backupH(is), ham%ham(is))
             end if
          end if
       end do
       do is=1,pub_num_spins
          call sparse_embed_create(dmft_self(is)   ,ham%ham(is),iscmplx=.false.)
          call sparse_embed_create(dmft_z(is)      ,ham%ham(is),iscmplx=.false.)
          if (pub_dmft_fully_sc_h) call sparse_embed_create(dmft_kernel(is), &
               denskern%m(is,PUB_1K),iscmplx=.false.)
       end do
    end if
    godmft=1

    ! ebl: In a DMFT calculation with self-consistency, loop until converged
    ! ebl: In normal calculations, just loop once.
    loop_ham_update = .true.
    do while (loop_ham_update)
       loop_ham_update = .false.
       dmft_correction = 0.0d0

       ! ndmh: include spin-degeneracy in density kernel
       if (pub_num_spins /= 2) then
          call sparse_embed_array_scale(denskern, pub_spin_fac)
       end if

       ! ndmh: calculate the lhxc potential
       if(.not.lhxc_fixed) then ! [*]
          ! jd: In the absence of polarisable embedding
          if(.not. pub_pol_emb_pot) then
             ! agrecokpt: check that in the TB method, it is correct to use
             ! the k-dependent overlap matrix and not the original one at Gamma!
             call hamiltonian_lhxc_calculate(lhxc_fine,lhxc_energy,ham%dijhat, &
                  mdl,rep%ngwfs_on_grid, &
                  ngwf_basis,denskern,rep%ngwf_overlap,rep%sp_overlap,&
                  dfdtau_fine = dfdtau_fine )
          else
             ! agrecokpt
             if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP) then
                if (pub_debug_on_root) write(stdout,'(a)') &
                   'WARNING: computation of polemb matrix in &
                   &hamiltonian_dens_dep_matrices not yet checked for &
                   &non-Gamma k-point'
             end if

             ! jd: For polarisable embedding potential briefly store density_fine
             !     calculated by lhxc_calculate. This avoids re-calculating it and
             !     keeping module-wide state.
             allocate(density_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
                  mdl%fine_grid%max_slabs12,pub_num_spins,mdl%nsub,mdl%nsub), stat=ierr)
             call utils_alloc_check('hamiltonian_dens_dep_matrices','density_fine',&
                  ierr)
             call hamiltonian_lhxc_calculate(lhxc_fine,lhxc_energy,ham%dijhat, &
                  mdl,rep%ngwfs_on_grid, &
                  ngwf_basis,denskern,rep%ngwf_overlap,rep%sp_overlap, &
                  return_density_fine = density_fine)
             ! jd: Then calculate and include the effect of polarisable embedding
             !     NB: this also *modifies* lhxc_fine and lhxc_energy to account
             !     for non-QM* embedding multipoles and/or MM repulsive
             !     potentials, if any.
             ! jcap: convert denskern first, extracting region 1
             call sparse_embed_array_to_sparse_array(temp_kern,denskern,1)
             ! jcap: also create temporary arrays for hamiltonian
             allocate(ham_array(6),stat=ierr)
             call utils_alloc_check('HDDM','ham_array',ierr)
             call sparse_embed_extract_from_array(ham_array(1:2),&
                  ham%polemb(1:2))
             if (pub_pol_emb_vacuum_qmstar) then
                call sparse_embed_extract_from_array(ham_array(3:4),&
                     ham%vacuum_polemb(1:2))
                call sparse_embed_extract_from_array(ham_array(5:6),&
                     ham%vacuum_aux_polemb(1:2))
             end if
             call polarisable_embedding_calculate(ham_array(1:2), &      ! inout
                  ham_array(3:4),ham_array(5:6),&                        ! inout
                  polemb_energy_mm_valence, &                            ! out
                  polemb_energy_mm_vdw, &                                ! out
                  polemb_energy_mm_es, &                                 ! out
                  polemb_energy_qmstar_elec_mm_perm, &                   ! out
                  polemb_energy_qmstar_elec_mm_ind, &                    ! out
                  polemb_energy_qmfull_elec_mm_perm, &                   ! out
                  polemb_energy_qmfull_elec_mm_ind, &                    ! out
                  polemb_energy_qm_core_mm_perm, &                       ! out
                  polemb_energy_qm_core_mm_ind, &                        ! out
                  polemb_energy_qmmm_vdw, &                              ! out
                  polemb_energy_other, &                                 ! out
                  lhxc_fine(:,:,:,:,1), lhxc_energy, &                   ! inout
                  density_fine(:,:,:,:,1,1), mdl, rep, ngwf_basis(1), &  ! in
                  temp_kern)                                             ! in

             ! jcap: destroy temporary kernel
             call sparse_array_destroy(temp_kern)
             ! jcap: copy results back over and destroy
             call sparse_embed_destroy_extracted_array(ham_array(1:2),&
                  ham%polemb(1:2),.true.)
             if (pub_pol_emb_vacuum_qmstar) then
                call sparse_embed_destroy_extracted_array(ham_array(3:4),&
                     ham%vacuum_polemb(1:2),.true.)
                call sparse_embed_destroy_extracted_array(ham_array(5:6),&
                     ham%vacuum_aux_polemb(1:2),.true.)
             end if

             deallocate(ham_array,stat=ierr)
             call utils_dealloc_check('HDDM','ham_array',ierr)

             ! jd: Clean up the temporary density buffer
             deallocate(density_fine,stat=ierr)
             call utils_dealloc_check('hamiltonian_dens_dep_matrices', &
                  'density_fine', ierr)
          end if
       end if

       ! ndmh: only calculate the lhxc matrix if we need to update the Hamiltonian
       if (ham_update.or.dmft_ham_update)then
          ! agrecokpt: need to sum over k-points, i.e. need ham%lhxc(is,ik) computed
          ! from ngwfs_on_grid(ik), then weighted sum over k-points; lhxc_fine is
          ! k-point NON-dependent in KP method and if using summed over density
          ! cks: calculate the lhxc matrix
          do jsub=1,ncols
             do isub=1,mrows
                do is=1,pub_num_spins
                   ! rc2013: the LHXC potential required will depend on the region we're in
                   lhxc_index = isub
                   if((isub .ne. jsub) .and. (isub == pub_active_region)) &
                        lhxc_index = jsub
                   ! rc2013: TEST: use "average" LHXC in cross regions
                   ! rc2013: add the EMFT potential to the active region only
                   if (utils_devel_code(.false., 'HAM','LHXC_AVG',pub_devel_code) &
                        .and. (isub.ne.jsub)) then
                      call integrals_locpot(ham%lhxc(is)%m(isub,jsub), &
                           rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                           rep%ngwfs_on_grid(jsub), ngwf_basis(jsub), mdl%fine_grid, &
                           mdl%dbl_grid, mdl%cell, mdl%fftbox, &
                           (lhxc_fine(:,:,:,is,isub)+lhxc_fine(:,:,:,is,jsub))/2.0_DP)
                   else
                      call integrals_locpot(ham%lhxc(is)%m(isub,jsub), &
                           rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                           rep%ngwfs_on_grid(jsub), ngwf_basis(jsub), mdl%fine_grid, &
                           mdl%dbl_grid, mdl%cell, mdl%fftbox, &
                           lhxc_fine(:,:,:,is,lhxc_index))
                   end if
                enddo
             enddo
          enddo
          do is=1,pub_num_spins
             ! agrecokpt: add TB method extra terms if requested
             if ((pub_kpoint_method == 'TB') .and. &
                (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP)) then
                if (pub_debug_on_root) write(stdout,'(a,i3)') 'DEBUG: &
                   &Adding TB terms in lhxc matrix spin ', is

                call hamiltonian_apply_phases_tb_method(ham%lhxc(is)%p, &
                     ngwf_basis(1), mdl, loc_kpt)

             end if

             ! agrecocmplx: check final lhxc matrix is hermitian if requested
             if (pub_check_hermitian) then
                rms_hermitian_check = sparse_embed_check_hermitian(ham%lhxc(is))

                if (pub_on_root) write(stdout,'(a,i3,a,f8.5)') 'RMS lhxc spin ',&
                   is,' hermitian check: ', rms_hermitian_check

                call utils_assert(rms_hermitian_check < pub_imag_thr, &
                     'WARNING: hermitian check for lhxc matrix above threshold')

             end if

             ! jd: Devel code hook to see what LHXC looks like after it has been
             !     projected to the NGWF basis. Comparison between original
             !     realspace LHXC and a realspace version reconstructed from
             !     ham%lhxc can reveal tricky issues with representing external
             !     potentials that are difficult to capture in the NGWF basis
             if(utils_devel_code(.false.,'DEBUG','DUMP_LHXC', pub_devel_code, &
                  no_warn = .true.)) then
                if(pub_inner_loop_iteration == 1 .and. .not. lhxc_output_done) then
                   ! rc2013: EMBED_FIX!
                   call visual_scalarfield_ngwf_basis(ham%lhxc(is)%p, rep, mdl, 1, &
                       ngwf_basis(1), &
                       'LHXC pot as represented in NGWF basis (internal units)', &
                        '_reconstructed_lhxc_pot_spin_'//&
                        trim(utils_int_to_str(is))//'_'//&
                        trim(utils_int_to_str(lhxc_output_iter)))
                   if(is == pub_num_spins) then
                      lhxc_output_iter = lhxc_output_iter + 1
                      lhxc_output_done = .true.
                   end if
                else
                   lhxc_output_done = .false.
                end if
             end if

          end do

          ! JCW: Calculate the dfdtau matrix, if KE-density-dependent functional
          ! JCW: is used
          if (pub_xc_ke_density_required) then
             ! agrecokpt
             if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP) then
                if (pub_debug_on_root) write(stdout,'(a)') &
                   'WARNING: computation of dfdtau matrix in &
                   &hamiltonian_dens_dep_matrices not yet checked for &
                   &non-Gamma k-point'
             end if

             ! JCW: For testing, set dfdtau_fine = -0.5_DP. The trace of the integrals
             ! JCW: with the density kernel should be equal to the kinetic energy
             ! JCW:  dfdtau_fine = -0.5_DP
             do jsub=1,ncols
                do isub=1,mrows
                   do is=1,pub_num_spins
                      ! rc2013: we probably need a separate dfdtau_fine for each
                      ! rc2013: region for embedding
                      call integrals_div_locpot_grad_functions( &
                           ham%dfdtau(is)%m(isub,jsub), &
                           rep%ngwfs_on_grid(isub),ngwf_basis(isub), &
                           rep%ngwfs_on_grid(jsub),ngwf_basis(jsub), &
                           mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
                           dfdtau_fine(:,:,:,is))
                   end do
                   if (pub_debug) then
                      ! JCW: Internal testing routine for div_locpot_grad integrals:
                      ! JCW: Requires dfdtau integral SPAM3 object (ham%dfdtau) to
                      ! JCW: set sparsity of internal SPAM3 objects. Any data in
                      ! JCW: the SPAM3 object passed in is ignored, however.
                      ! JCW: as well as arguments necessary to evaluate new sets of
                      ! JCW: dfdtau integrals.
                      ! To activate test, need
                      !   devel_code: HAM:TAUTEST=T:HAM
                      ! in input file.
                      if (utils_devel_code(&
                         .false.,'HAM','TAUTEST',pub_devel_code)) then
                         call test_integrals_div_locpot_grad_functions( &
                              ham%dfdtau(1)%m(isub,jsub), &
                              rep%ngwfs_on_grid(isub),ngwf_basis(isub), &
                              rep%ngwfs_on_grid(jsub),ngwf_basis(jsub), &
                              mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
                              dfdtau_fine(:,:,:,:),denskern)
                      end if
                   end if
                end do
             end do
          end if
       end if

       ! gcc32:
       if(pub_confined_ngwfs) then
          if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP) then
             if (pub_debug_on_root) write(stdout,'(a)') &
                'WARNING: computation of confinement matrix in &
                &hamiltonian_dens_dep_matrices not yet checked for &
                &non-Gamma k-point'
          end if

          do is=1,pub_num_spins
             do jsub=1,mdl%nsub
                do isub=1,mdl%nsub
                   ! rc2013: EMBED_FIX! this is designed for only 1 subsystem
                   call integrals_confinement_pot( &
                        ham%confinement(is)%m(isub,jsub), &
                        rep%ngwfs_on_grid(isub),ngwf_basis(isub), &
                        rep%ngwfs_on_grid(jsub),ngwf_basis(jsub), &
                        mdl%dbl_grid, mdl%cell, mdl%fftbox)
                end do
             end do
          end do
       end if

       ! qoh: Calculate and include HF exchange if necessary
       ! jcap: this part only works for doing a non-embedding calculation
       ! rc2013: this format is obsolete, will need updating
       if (pub_use_hfx) then
          ! rc2013: EMBED_FIX!
          do isub=1,mdl%nsub
             call sparse_embed_extract_from_array(kern_array,&
                  denskern%m(:,PUB_1K),isub,isub)
             call hf_exchange_calculate(hfxstate, ham%hfexchange, rep, &
                  isub, kern_array, kern_array, ngwf_basis(isub), &
                  mdl%fftbox, mdl%cell, mdl%regions(isub)%elements, .false.)
             call sparse_embed_destroy_extracted_array(kern_array)
          end do
          ! jcap: if we only have exact exchange in the active
          ! subregion (and we are doing an EMFT calculation), pass
          ! appropriate arguments
       else if (pub_use_activehfx.and.pub_emft) then
          ! rc2013: EMBED_FIX!
          call sparse_embed_extract_from_array(kern_array,&
               denskern%m(:,PUB_1K),pub_active_region,pub_active_region)
          do is=1,pub_num_spins
             ! rc2013: zero the inactive components of the exchange matrix
             call sparse_embed_scale(ham%hfexchange(is), 0.0_DP)
          end do
          ! jcap: swap functionals for EMFT, to get pub_hfxfraction
          ! right
          call xc_embed_swap_functional(.true.)
          call hf_exchange_calculate(hfxstate, ham%hfexchange, rep, &
               pub_active_region, kern_array, kern_array, &
               ngwf_basis(pub_active_region), mdl%fftbox, mdl%cell, &
               mdl%regions(pub_active_region)%elements, .false.)
          call xc_embed_swap_functional(.false.)
          call sparse_embed_destroy_extracted_array(kern_array)
       else if (pub_use_activehfx.and.pub_emft_follow) then
          ! rc2013: if we'll be doing an EMFT calculation after a low-level calculation,
          ! explicitly zero the hfexchange matrix
          do is=1,pub_num_spins
             call sparse_embed_scale(ham%hfexchange(is), 0.0_DP)
          end do
       end if

       ! ndmh: calculate Hubbard energy and Hamiltonian matrix if necessary
       if (pub_hubbard) then
          if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP) then
             if (pub_debug_on_root) write(stdout,'(a)') &
                  'WARNING: computation of hubbard matrix in &
                  &hamiltonian_dens_dep_matrices not yet checked for &
                  &non-Gamma k-point'
          end if

          if (pub_num_spins /= 2) then
             call sparse_embed_array_scale(denskern, 0.5_DP*pub_num_spins)
          end if
          ! jcap: create temporary arrays
          allocate(ham_array(pub_num_spins),stat=ierr)
          call utils_alloc_check('HDDM','ham_array',ierr)
          call sparse_embed_extract_from_array(kern_array,denskern%m(:,PUB_1K))
          call sparse_embed_extract_from_array(ham_array,ham%hubbard_ham)
          call hubbard_energy_total(hub,hubbard_energy, &
               kern_array,hub_proj_basis(1),rep%hub_overlap%p,&
               rep%hub_overlap_t%p)
          if (ham_update.or.dmft_ham_update) then
             call hubbard_projector_ham(hub, &
                  kern_array,rep%hub_overlap%p,rep%hub_overlap_t%p,&
                  hub_proj_basis(1))
             call hubbard_ham_matrix(hub,ham_array, &
                  rep%hub_overlap%p,rep%hub_overlap_t%p)
          end if
          call sparse_embed_destroy_extracted_array(kern_array,&
               denskern%m(:,PUB_1K),.true.)
          call sparse_embed_destroy_extracted_array(ham_array,&
               ham%hubbard_ham,.true.)
          deallocate(ham_array,stat=ierr)
          call utils_dealloc_check('HDDM','ham_array',ierr)
          if (pub_num_spins /= 2) then
             call sparse_embed_array_scale(denskern, pub_spin_fac)
          end if
       else
          hubbard_energy = 0.0_DP
       end if

       ! ndmh: calculate nonlocal potential matrix
       if (pub_aug) then
          ! agrecokpt
          if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP) then
             if (pub_debug_on_root) write(stdout,'(a)') &
                'WARNING: computation of nonlocal matrix in &
                &hamiltonian_dens_dep_matrices not yet checked for &
                &non-Gamma k-point'
          end if

          ! Create temporary matrix for projector denskern
          allocate(projector_denskern(pub_num_spins),stat=ierr)
          call utils_alloc_check('hamiltonian_dens_dep_matrices', &
               'projector_denskern',ierr)
          do is=1,pub_num_spins
             projector_denskern(is)%structure = 'E'
             call sparse_embed_create(projector_denskern(is))
          end do

          ! Create the PAW projector density kernel
          ! jcap: need temporary arrays to pass in
          allocate(ham_array(pub_num_spins),stat=ierr)
          call utils_alloc_check('hamiltonian_dens_dep_matrices',&
               'ham_array',ierr)
          allocate(aug_array(pub_num_spins),stat=ierr)
          call utils_alloc_check('hamiltonian_dens_dep_matrices',&
               'aug_array',ierr)
          call sparse_embed_extract_from_array(aug_array,projector_denskern)
          call sparse_embed_extract_from_array(kern_array,denskern%m(:,PUB_1K))

          call aug_projector_denskern(aug_array, kern_array, &
               rep%sp_overlap%p)

          call sparse_embed_destroy_extracted_array(aug_array,projector_denskern,&
               .true.)
          call sparse_embed_destroy_extracted_array(kern_array)

          ! Calculate the PAW nonlocal matrix
          call sparse_embed_extract_from_array(ham_array,ham%nonlocpot)
          call sparse_embed_extract_from_array(aug_array,ham%dijhat)
          call sparse_embed_extract_from_array(kern_array,projector_denskern)

          call aug_nonlocal_mat(ham_array,aug_array,kern_array,&
               rep%sp_overlap%p,mdl%regions(1)%pseudo_sp, &
               mdl%regions(1)%paw_sp,paw_sphere_energies,.false.)

          call sparse_embed_destroy_extracted_array(ham_array,ham%nonlocpot,&
               .true.)
          call sparse_embed_destroy_extracted_array(aug_array)
          call sparse_embed_destroy_extracted_array(kern_array)
          deallocate(ham_array,stat=ierr)
          call utils_dealloc_check('hamiltonian_dens_dep_matrices',&
               'ham_array',ierr)
          deallocate(aug_array,stat=ierr)
          call utils_dealloc_check('hamiltonian_dens_dep_matrices',&
               'aug_array',ierr)

          ! agrecokpt: add TB method extra terms if requested
          ! check this is correct in PAW case
          do is=1,pub_num_spins
             if ((pub_kpoint_method == 'TB') .and. (loc_kpt%x/=0.0_DP.or.&
                  loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP)) then
                if (pub_debug_on_root) write(stdout,'(a,i3)') 'DEBUG: &
                   &Adding TB terms in nonlocpot matrix spin ', is

                call hamiltonian_apply_phases_tb_method(ham%nonlocpot(is)%p, &
                     ngwf_basis(1), mdl, loc_kpt)

             end if

             ! agrecocmplx: check nonlocpot matrix is hermitian if requested
             if (pub_check_hermitian) then
                rms_hermitian_check = sparse_embed_check_hermitian(ham%nonlocpot(is))

                if (pub_on_root) write(stdout,'(a,i3,a,f8.5)') 'RMS nonlocpot spin ',&
                   is,' hermitian check: ', rms_hermitian_check

                call utils_assert(rms_hermitian_check < pub_imag_thr, &
                     'WARNING: hermitian check for nonlocpot matrix above threshold')

             end if
          end do

          ! Clean up temporary matrices
          do is=pub_num_spins,1,-1
             call sparse_embed_destroy(projector_denskern(is))
          end do
          deallocate(projector_denskern,stat=ierr)
          call utils_dealloc_check('hamiltonian_dens_dep_matrices', &
               'projector_denskern',ierr)
       end if

       ! Add up local potential, Hartree, exchange-correlation, Ewald, Hubbard
       ! and dispersion energies.
       ! NB: If polarisable embedding in effect, lhxc_energy already contains
       !     those polemb contributions use the local potential model.
       total_energy = lhxc_energy + mdl%ewald_energy + hubbard_energy + &
            pub_dispersion_energy

       ! jd: If polarisable embedding, add relevant energy terms.
       !     Note that:
       !     a) Non-QM* (locpot) energy terms involving MM multipoles interac-
       !        ting with electrons have already been included in lhxc_energy
       !        and lhxc_fine and are not added again, except that there is an
       !        indirect contribution in mm_es from MM+ pol.
       !     b) QM* <-> MM energy terms give rise to qmstar terms.
       !     c) Regardless of scheme, there are QM core <-> MM terms.
       if(pub_pol_emb_pot) then
          total_energy = total_energy + &
               polemb_energy_mm_valence + polemb_energy_mm_vdw + &
               polemb_energy_mm_es + polemb_energy_qmstar_elec_mm_perm + &
               polemb_energy_qmstar_elec_mm_ind + &
               polemb_energy_qm_core_mm_perm + polemb_energy_qm_core_mm_ind + &
               polemb_energy_qmmm_vdw + polemb_energy_other
       end if

       ! Add kinetic, nonlocal, hfexchange and polarisable embedding energies
       do is=1,pub_num_spins
          call sparse_embed_trace(trace,denskern%m(is,PUB_1K),rep%kinet)
          total_energy = total_energy + trace
          ! rc2013: print energy components for debugging
          if(pub_debug) then
              if(pub_on_root) write(stdout,'(a17,f15.6)') 'kinetic_energy = ', &
                  trace
          endif

          ! ndmh: add nonlocal pseudopotential energy if necessary
          if (pub_any_nl_proj.and.(.not.pub_usp)) then
              call sparse_embed_trace(trace,denskern%m(is,PUB_1K),rep%nonlocpot)
              total_energy = total_energy + trace
              if(pub_debug) then
                  if(pub_on_root) write(stdout,'(a15,f15.6)') 'Nonlocal energy = ', &
                      trace
              endif
          end if
          if (pub_aug) then
             call sparse_embed_trace(trace,denskern%m(is,PUB_1K),ham%nonlocpot(is))
             total_energy = total_energy + trace
             if(pub_debug) then
                if(pub_on_root) write(stdout,'(a15,f15.6)') 'Augmentation energy = ', &
                     trace
             endif
          end if

          ! qoh: include HF exchange if necessary
          ! jcap: for cases where we've calculated it
          if (pub_use_hfx.or.(pub_use_activehfx.and.pub_emft)) then
              call sparse_embed_trace(trace,denskern%m(is,PUB_1K),ham%hfexchange(is))
              total_energy = total_energy - 0.5_DP*trace
           end if

       end do
       if(pub_debug_on_root) write(stdout,'(a15,f15.6)') &
            'DEBUG: total_energy = ', total_energy

       ! Add up PAW energy terms
       if (pub_paw) then
          ! ndmh: subtract unnecessary bits of nonlocal energy, add AE sphere term
          total_energy = total_energy &
                 - paw_sphere_energies(paw_en_dijhat) &
                 - paw_sphere_energies(paw_en_ehart) &
                 + paw_sphere_energies(paw_en_exc) &
                 - paw_sphere_energies(paw_en_etxc) &
                 - paw_sphere_energies(paw_en_dijxc) &
                 - paw_sphere_energies(paw_en_exc_core)
       end if

       ! jd: Include the smeared ion energy corrections, if necessary
       if (pub_is_smeared_ion_rep) then
          total_energy = total_energy + smeared_ion_E_self + smeared_ion_E_smeared
       end if

       ! ebl
       if (dmft_ham_update.and.godmft==1) then

          call hamiltonian_build_matrix(ham, rep)

          if(pub_on_root) then
             call file_handling_remove_file(trim(pub_rootname)//".eigen")
             call file_handling_remove_file(trim(pub_rootname)//".ham1")
             call file_handling_remove_file(trim(pub_rootname)//".ham2")
             call file_handling_remove_file(trim(pub_rootname)//".hub_overlap")
             call file_handling_remove_file(trim(pub_rootname)//".hub_overlap_t")
             call file_handling_remove_file(trim(pub_rootname)//".inv_overlap")
             call file_handling_remove_file(trim(pub_rootname)//".overlap")
             call file_handling_remove_file(trim(pub_rootname)//".tot_energy")
             call file_handling_remove_file(trim(pub_rootname)//".nabla")
             call file_handling_remove_numbered_files("chem.potential.nmu.iter")
             call file_handling_remove_numbered_files(trim(pub_rootname)//".occupancy")
          endif
          ! ebl
          ! The idea here was to replace all together H with [H+Sigma(oo)]*Z in
          ! the onetep energy minimisation. You would have to do DMFT
          ! calculations for every single evaluation of the hamiltonian, which
          ! would be quite useful as it would affect the NGWF optimisation as well.
          ! However, it is quite expensive and is not a priority to implement
          ! if(.not.pub_dmft_fully_sc_h)then
          !      call hubbard_dmft_interface(eigen_en,rep%n_occ,                           &
          !           ham%ham,rep%overlap,rep%inv_overlap,ngwf_basis,    &
          !           hub_proj_basis,hub,rep,elements,denskern,mdl,      &
          !           dmft_energy_cor=dmft_correction,dmft_self=dmft_self,dmft_z=dmft_z)
          ! else
          !      call hubbard_dmft_interface(eigen_en,rep%n_occ,                           &
          !           ham%ham,rep%overlap,rep%inv_overlap,ngwf_basis,    &
          !           hub_proj_basis,hub,rep,elements,denskern,mdl,      &
          !           dmft_energy_cor=dmft_correction,dmft_kernel=dmft_kernel)
          ! endif
          ! if(pub_my_proc_id==0) then
          !    ressys=system(" rm store*                   ")
          !    ressys=system(" rm chem.potential.nmu.iter* ")
          ! endif
          ! if(pub_my_proc_id==0) write(*,*) 'DMFT CORRECTION : ', dmft_correction
          !
          ! H -> Z * (H+Sig(oo))
          !
          ! do is=1,pub_num_spins
          !    if(ham_update.and..not.pub_dmft_fully_sc_h) then
          !        call sparse_axpy(ham%ham(is),dmft_self(is),1.d0)
          !        call sparse_copy(dmft_self(is),ham%ham(is))
          !        call sparse_product(ham%ham(is),dmft_self(is),dmft_z(is))
          !    endif
          call sparse_embed_destroy(dmft_z(is))
          call sparse_embed_destroy(dmft_self(is))
          ! enddo
          if(pub_dmft_fully_sc_h.and.godmft==1)then
             godmft=2
             loop_ham_update = .true.
          endif
       endif
    end do

    ! ny: species-dependent scissor operator
    if (pub_scissor_ngroups > 0) then
       call hamiltonian_species_scissor(ham%scissor,denskern%m(:,PUB_1K), &
            rep,ngwf_basis,mdl)
       do is=1,pub_num_spins
          call sparse_embed_trace(trace,denskern%m(is,PUB_1K),ham%scissor(is))
          total_energy = total_energy + trace
       enddo
    end if

    if (pub_dmft_fully_sc.and.pub_dmft_fully_sc_h)then
      do is=1,pub_num_spins
         call sparse_embed_copy(denskern%m(is,PUB_1K),dmft_kernel(is))
         call sparse_embed_destroy(dmft_kernel(is))
      enddo
    endif

    if (.not.ham_update.and.dmft_ham_update) then
      do is=1,pub_num_spins
         call sparse_embed_copy(ham%ham(is),backupH(is))
         call sparse_embed_destroy(backupH(is))
      enddo
    endif

    total_energy = total_energy + dmft_correction

    ! gcc32
    if (pub_confined_ngwfs) then
      ! add confinement "energy"
       do is=1,pub_num_spins
          call sparse_embed_trace(trace,denskern%m(is,PUB_1K),ham%confinement(is))
          total_energy = total_energy + trace
      end do
    end if

    ! ndmh: remove spin-degeneracy from density kernel
    if (pub_num_spins /= 2) then
       call sparse_embed_array_scale(denskern, 0.5_DP*pub_num_spins)
    end if

#ifdef GPU_SP_TEST
    ! kaw: Load values for precision switching
    deltaE_gpu = abs(E_gpu_previous-total_energy)
    E_gpu_previous = total_energy
#endif

    deallocate(kern_array,stat=ierr)
    call utils_dealloc_check('HDDM','kern_array',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hamiltonian_dens_dep_matrices'

  end subroutine hamiltonian_dens_dep_matrices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hamiltonian_dens_dep_nonsc(ham, rep, ngwf_basis, mdl, hfxstate, &
       lhxc_fine, hub,val_rep,val_ham,val_dkn,updated_shift,dfdtau_fine, &
       val_ngwf_basis, cond_dkn)

    !====================================================================!
    ! Updates the density-matrix dependent parts of the Hamiltonian      !
    ! when the NGWFs have changed, without recalculating the density     !
    ! (ie non-self-consistently).                                        !
    ! @docme                                                             !
    ! val_ngwf_basis(in, opt): Only needed for HFx with mixed NGWF bases.!
    ! cond_dkn(in, opt):       Only needed for HFx with mixed NGWF bases.!
    !--------------------------------------------------------------------!
    ! Written by Laura Ratcliff June 2010 as hamiltonian_cond_ham_update !
    ! Moved to hamiltonian_mod and rearranged by Nicholas Hine in April  !
    ! 2011.                                                              !
    !====================================================================!

    use augmentation, only: aug_projector_denskern, aug_nonlocal_mat
    use constants, only: paw_en_size, stdout, REP_SWEX_HFX_OTHER
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE, hf_exchange_calculate
    use hubbard_build, only: HUBBARD_MODEL, hubbard_ham_matrix
    use integrals, only: integrals_locpot, integrals_div_locpot_grad_functions
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_aug, pub_hubbard, pub_use_hfx, pub_num_spins, &
         pub_spin_fac, pub_xc_ke_density_required, pub_emft, &
         pub_use_activehfx, pub_active_region
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_destroy
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_scale, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_assert
    use xc, only: xc_embed_swap_functional

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(inout)   :: ham
    type(NGWF_REP), intent(in)      :: rep
    type(FUNC_BASIS), intent(in)    :: ngwf_basis(rep%nsub)
    type(MODEL), intent(in)         :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    real(kind=DP), intent(inout)    :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins,rep%nsub)
    type(HUBBARD_MODEL), intent(in) :: hub
    type(NGWF_HAM), intent(in)      :: val_ham
    type(NGWF_REP), intent(in)      :: val_rep
    type(SPAM3_EMBED), intent(in)   :: val_dkn(pub_num_spins)
    logical, optional, intent(out)  :: updated_shift
    ! JCW: Optional gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density. Required to evaluate lhxc matrix elements for meta-GGAs.
    real(kind=DP), optional, intent(in) :: dfdtau_fine(:,:,:,:)
    ! JCW: In general, we have
    ! JCW:    dfdtau_fine(mdl%fine_grid%ld1,&
    ! JCW:    mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    ! JCW: or
    ! JCW:    dfdtau_fine(1,1,1,1)
    ! JCW: when it is unneeded.

    ! jd: Optionals for conduction with HFx
    type(FUNC_BASIS), optional, intent(in) :: val_ngwf_basis(val_rep%nsub)
    type(SPAM3_EMBED), optional, intent(in) :: cond_dkn(pub_num_spins)

    ! Local Variables
    type(SPAM3_EMBED), allocatable :: projector_denskern(:)
    type(SPAM3),dimension(pub_num_spins) :: val_kern_array,cond_kern_array,ham_array,aug_array
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    integer :: ierr
    integer :: is,isub,jsub
    integer :: lhxc_index

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hamiltonian_dens_dep_nonsc'

    if(pub_use_hfx) then
       call utils_assert(present(val_ngwf_basis), &
            'hamiltonian_dens_dep_nonsc(): The optional argument val_ngwf_basis&
            & needs to be passed when HFx is enabled.')
       call utils_assert(present(cond_dkn), &
            'hamiltonian_dens_dep_nonsc(): The optional argument cond_dkn&
            & needs to be passed when HFx is enabled.')
    end if

    ! lr408: calculate the lhxc matrix
    ! jcap: loop over regions
    do jsub=1,mdl%nsub
       do isub=1,mdl%nsub
          do is=1,pub_num_spins
             ! rc2013: the LHXC potential required will depend on the region we're in
             lhxc_index = isub
             if((isub .ne. jsub) .and. (isub == pub_active_region)) &
                  lhxc_index = jsub
             call integrals_locpot(ham%lhxc(is)%m(isub,jsub), &
                  rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                  rep%ngwfs_on_grid(jsub), ngwf_basis(jsub), mdl%fine_grid, &
                  mdl%dbl_grid, mdl%cell, mdl%fftbox, &
                  lhxc_fine(:,:,:,is,lhxc_index))
          enddo
       enddo
    enddo

    ! JCW: Abort if tau-dependent XC functional requested, since this has not
    ! JCW: been implemented in hamiltonian_dens_dep_nonsc (yet)
    if (pub_xc_ke_density_required) then
      call utils_abort("Error in hamiltonian_dens_dep_nonsc: &
           &A tau-dependent functional has been requested, but support for &
           &this has not been implemented in hamiltonian_dens_dep_nonsc.")
    end if
    ! JCW: ... also abort if optional argument dfdtau_fine is present (and with
    ! JCW: size > 1, since we accept a size == 1 dummy array)
    if ( present(dfdtau_fine) ) then
       call utils_assert( size(dfdtau_fine) == 1, &
            "Error in hamiltonian_dens_dep_nonsc: &
            &Optional argument dfdtau_fine has been passed in and has size /= 1 &
            &but support for this has not been implemented in &
            &hamiltonian_dens_dep_nonsc.")
    end if
    ! ndmh: Calculate and include HF exchange if necessary
    if (pub_use_hfx) then
        ! rc2013: EMBED_FIX!
       if(rep%postfix=='j' .or. rep%postfix=='c') then
          do isub=1,mdl%nsub
             call sparse_embed_extract_from_array(val_kern_array,val_dkn,isub,isub)
             call sparse_embed_extract_from_array(cond_kern_array,cond_dkn,isub,isub)
             ! jd: For 'cond' calculation we have cond_dkn as denskern_ab
             !     and val_dkn as denskern_cd. For 'joint' calculation cond_dkn
             !     is actually joint_dkn (passed from conduction_properties), and
             !     val_dkn is the valence kernel. Currently conduction_properties
             !     does not initialise joint_dkn (it is all zeroes).
             call hf_exchange_calculate(hfxstate, ham%hfexchange, rep, &
                  isub, cond_kern_array, val_kern_array, ngwf_basis(isub), &
                  mdl%fftbox, mdl%cell, mdl%regions(isub)%elements, .false., &
                  rep2 = val_rep, &
                  ngwf_basis2 = val_ngwf_basis(isub), &
                  basis_selector = (/1,1,2,2,REP_SWEX_HFX_OTHER/), &
                  !                  ^ qoh's alpha-beta-gamma-delta conv'n
                  energy_prefactor = 2.0/real(pub_num_spins,kind=DP))
             call sparse_embed_destroy_extracted_array(val_kern_array)
             call sparse_embed_destroy_extracted_array(cond_kern_array)
          end do
       else
          call utils_abort('Unsupported postfix in HFx call: "'//trim(rep%postfix)//'"')
       end if
       ! jcap: if we only have exact exchange in the active
       ! subregion (and we are doing an EMFT calculation), pass
       ! appropriate arguments
    else if (pub_use_activehfx.and.pub_emft) then
       if(rep%postfix=='j' .or. rep%postfix=='c') then
          ! rc2013: EMBED_FIX!
          call sparse_embed_extract_from_array(val_kern_array,val_dkn,&
               pub_active_region,pub_active_region)
          call sparse_embed_extract_from_array(cond_kern_array,cond_dkn,&
               pub_active_region,pub_active_region)
          do is=1,pub_num_spins
             ! rc2013: zero the inactive components of the exchange matrix
             call sparse_embed_scale(ham%hfexchange(is), 0.0_DP)
          end do
          ! jcap: swap functionals for EMFT, to get pub_hfxfraction
          ! right
          call xc_embed_swap_functional(.true.)
          call hf_exchange_calculate(hfxstate, ham%hfexchange, rep, &
               pub_active_region, cond_kern_array, val_kern_array, &
               ngwf_basis(pub_active_region), mdl%fftbox, mdl%cell, &
               mdl%regions(pub_active_region)%elements, .false., &
               rep2 = val_rep, &
               ngwf_basis2 = val_ngwf_basis(pub_active_region), &
               basis_selector = (/1,1,2,2,REP_SWEX_HFX_OTHER/), &
               !                  ^ qoh's alpha-beta-gamma-delta conv'n
               energy_prefactor = 2.0/real(pub_num_spins,kind=DP))
          call xc_embed_swap_functional(.false.)
          call sparse_embed_destroy_extracted_array(val_kern_array)
          call sparse_embed_destroy_extracted_array(cond_kern_array)
       else
          call utils_abort('Unsupported postfix in HFx call: "'//trim(rep%postfix)//'"')
       end if
    end if

    ! ndmh: Calculate and include Hubbard Hamiltonian if necessary
    if (pub_hubbard) then
       call sparse_embed_extract_from_array(ham_array,ham%hubbard_ham)
       call hubbard_ham_matrix(hub,ham_array, &
            rep%hub_overlap%p,rep%hub_overlap_t%p)
       call sparse_embed_destroy_extracted_array(ham_array,ham%hubbard_ham,.true.)
    end if

    ! ndmh: in PAW, recalculate the nonlocal matrix
    if (pub_aug) then

       ! Create temporary matrix for projector denskern
       allocate(projector_denskern(pub_num_spins),stat=ierr)
       call utils_alloc_check('hamiltonian_dens_dep_nonsc', &
            'projector_denskern',ierr)
       do is=1,pub_num_spins
          projector_denskern(is)%structure = 'E'
          call sparse_embed_create(projector_denskern(is))
       end do

       ! Create the PAW projector density kernel
       ! agrecocmplx: projector_denskern is real in any case
       ! jcap: need temporary arrays to pass in
       call sparse_embed_extract_from_array(aug_array,projector_denskern)
       call sparse_embed_extract_from_array(val_kern_array,val_dkn)

       call aug_projector_denskern(aug_array,val_kern_array, &
            val_rep%sp_overlap%p)

       call sparse_embed_destroy_extracted_array(aug_array,projector_denskern,&
               .true.)
       call sparse_embed_destroy_extracted_array(val_kern_array)
       do is=1,pub_num_spins
          call sparse_embed_scale(projector_denskern(is),pub_spin_fac)
       end do

       ! Calculate the PAW nonlocal matrix
       call sparse_embed_extract_from_array(ham_array,ham%nonlocpot)
       call sparse_embed_extract_from_array(aug_array,ham%dijhat)
       call sparse_embed_extract_from_array(val_kern_array,projector_denskern)

       call aug_nonlocal_mat(ham_array,aug_array,val_kern_array,&
            rep%sp_overlap%p,mdl%regions(1)%pseudo_sp, &
            mdl%regions(1)%paw_sp,paw_sphere_energies,.false.)

       call sparse_embed_destroy_extracted_array(ham_array,ham%nonlocpot,.true.)
       call sparse_embed_destroy_extracted_array(aug_array)
       call sparse_embed_destroy_extracted_array(val_kern_array)

       ! Clean up temporary matrices
       do is=pub_num_spins,1,-1
          call sparse_embed_destroy(projector_denskern(is))
       end do
       deallocate(projector_denskern,stat=ierr)
       call utils_dealloc_check('hamiltonian_dens_dep_nonsc', &
            'projector_denskern',ierr)
    end if

    ! lr408: Build unprojected conduction hamiltonian
    call hamiltonian_build_matrix(ham, rep)

    ! lr408: Finally build projected matrix
    if (present(updated_shift)) then
       call hamiltonian_proj_cond_matrix(rep, ham, &
            val_rep, val_ham, val_dkn, shift_changed=updated_shift)
    end if

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hamiltonian_dens_dep_nonsc'

  end subroutine hamiltonian_dens_dep_nonsc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hamiltonian_lhxc_calculate(lhxc_fine, lhxc_energy, dijhat, &
       mdl, ngwfs_on_grid, ngwf_basis, &
       denskern, overlap, sp_overlap, add_xc_pot, dfdtau_fine, &
       return_density_fine)

    !=====================================================================!
    ! This subroutine calculates and returns the sum of the               !
    ! local pseudopotential, hartree potential and exchange-correlation   !
    ! potentials and their energy.                                        !
    !---------------------------------------------------------------------!
    ! @docme arguments                                                    !
    !                                                                     !
    ! If implicit solvation is used with a self-consistently adapting     !
    ! cavity, the returned potential also includes a gradient correction  !
    ! that reflects the dependence of dielectric permittivity on density. !
    !                                                                     !
    ! jd: If the optional argument 'return_density_fine' is passed, it is !
    ! used as a buffer for density_fine. It is then assumed to have been  !
    ! allocated prior to the call and the responsibility for freeing it   !
    ! lies with the caller. In normal operation this argument is omitted  !
    ! and the buffer is allocated and deallocated here.                   !
    !---------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 18/7/2001.                      !
    ! Modified by Chris-Kriton Skylaris on 14/02/2004 so that it          !
    ! is less memory-hungry.                                              !
    ! Modified by Peter Haynes on 1/7/2004 to use fourier parallelisation.!
    ! Modified by Nicholas Hine on 6/2/2009 for NLCC core charges         !
    ! Modified by David D. O'Regan on 6/5/2009 for TDDFT                  !
    ! Modified by Nicholas Hine in July 2009 for function basis type and  !
    ! for SPAM3                                                           !
    ! Modified by Nicholas Hine in July 2010 for PAW                      !
    ! Modified by Nicholas Hine in February 2011 for augmentation_mod     !
    ! Modified by Nicholas Hine in March 2011 to make logic clearer and   !
    ! easier to read.                                                     !
    ! Modified by Jacek Dziedzic in May-Dec 2015 to allow devel_codes to  !
    ! control dumps of density and potential. Also to optionally return   !
    ! the calculated density.                                             !
    ! Modified for subsystem calculations and embedded mean-field theory  !
    ! by Robert Charlton and Joseph Prentice, 2017-18.                    !
    !=====================================================================!

    use augmentation, only: augmentation_density_on_grid, &
         augmentation_screen_dij
    use datatypes, only: FUNCTIONS
    use comms, only: comms_bcast, comms_barrier, pub_on_root, pub_my_proc_id
    use constants, only: DP, stdout, VERBOSE, max_spins
    use cutoff_coulomb, only: cutoff_coulomb_hartree
    use density, only: density_on_grid
    use enthalpy, only: enthalpy_terms
    use function_basis, only: FUNC_BASIS
    use hartree, only: hartree_on_grid, hartree_via_multigrid
    use integrals, only: integrals_product_on_grid
    use is_poisson, only: SOLVATION_ENERGY_TERMS, zero_solvation_terms
    use is_solvation_boltzmann, only: implicit_solvent_boltzmann_energy_contrib
    ! JCW: Add support for kinetic-energy-density-dependent XC functionals
    use ke_density, only: ke_density_on_grid
    use model_type, only: MODEL
    use rundat, only: pub_coulomb_cutoff, pub_nlcc, pub_paw, &
         pub_multigrid_hartree, pub_aug_den_dim, pub_aug, &
         pub_nhat_in_xc, pub_external_pressure, pub_num_spins, &
         pub_turn_off_hartree, pub_spin_fac, pub_num_kpoints, PUB_1K, &
         pub_xc_ke_density_required, pub_hhf_nstates, pub_hhf_factor, &
         pub_debug, pub_emft, pub_active_region, pub_devel_code, pub_inner_loop_iteration
    use sparse, only: sparse_copy, sparse_create, sparse_destroy, SPAM3
    use sparse_array, only: SPAM3_ARRAY, sparse_array_create, sparse_array_copy, &
        sparse_array_destroy
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_abort, utils_devel_code
    use xc, only: xc_energy_potential, xc_init, xc_exit, xc_hfxinit, &
         xc_gradients, xc_emft_calculate

    implicit none

    ! Arguments
    type(MODEL), intent(in)             :: mdl
    type(FUNC_BASIS), intent(in)        :: ngwf_basis(mdl%nsub)
    type(SPAM3_EMBED_ARRAY), intent(in) :: denskern
    type(SPAM3_EMBED), intent(in)       :: overlap
    type(SPAM3_EMBED), intent(in)       :: sp_overlap
    type(FUNCTIONS), intent(in)   :: ngwfs_on_grid(mdl%nsub)
    real(kind=DP), intent(out)    :: lhxc_energy
    real(kind=DP), intent(out)    :: lhxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub)
    type(SPAM3_EMBED), intent(inout)    :: dijhat(:)
    logical, intent(in), optional :: add_xc_pot
    ! JCW: Optional gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density. Required to evaluate lhxc matrix elements for meta-GGAs.
    real(kind=DP), optional, intent(out) :: dfdtau_fine(:,:,:,:)
    ! JCW: In general, we have
    ! JCW:    dfdtau_fine(mdl%fine_grid%ld1,&
    ! JCW:      mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    ! JCW: or
    ! JCW:    dfdtau_fine(1,1,1,1)
    ! JCW: when it is unneeded.
    real(kind=DP), intent(inout), optional, target :: return_density_fine(:,:,:,:,:,:)

    ! Local Variables
    integer :: ierr   ! error flag
    integer :: is
    logical :: loc_add_xc_pot
    integer :: fine_max_slabs12, fine_ld1, fine_ld2 ! jd: shortcut notation
    ! JCW: evaluate ke_density for use in meta-GGAs, if needed
    real(kind=DP), pointer, dimension(:,:,:,:)       :: ke_density_fine
    real(kind=DP), pointer, dimension(:,:,:,:)       :: density_fine
    real(kind=DP), pointer, dimension(:,:,:,:)       :: pot_fine
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: nhat_den_grad
    real(kind=DP) :: xc_energy
    real(kind=DP) :: enthalpy_energy
    real(kind=DP) :: locpspot_energy(mdl%nsub, mdl%nsub)   ! ndmh: energy of density in local pspot
    real(kind=DP) :: hartree_energy(mdl%nsub, mdl%nsub)    ! jd: Multigrid Hartree energy
    type(SOLVATION_ENERGY_TERMS) :: solvation_terms
    integer       :: isub, jsub, mrows, ncols, ik, kk, ll
    ! rc2013: temporary matrix to hold kernel info; avoids array temporaries
    type(SPAM3), allocatable :: kern_array(:)
    ! rc2013: placeholder for calculating XC energies
    ! jcap: placeholders for calculating XC energies within EMFT
    real(kind=DP) :: xc_energy_emft
    ! rc2013: for calculating gradient
    real(kind=DP), allocatable, target, dimension(:,:,:,:,:,:) :: sub_ke_density_fine
    real(kind=DP), pointer, dimension(:,:,:,:,:,:) :: sub_density_fine
    real(kind=DP), allocatable, target, dimension(:,:,:,:,:,:) :: sub_pot_fine
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: full_density_grad
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: emft_density_grad
    real(kind=DP), allocatable, dimension(:,:,:,:) :: demft_dfull
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: recip_work
    integer :: i1,i2,islab12                       ! Grid counters
    real(kind=DP) :: grad1(3), grad2(3), grad1_abs, grad2_abs
    ! rc2013: pass this to XC_mod to calculate the subsystem potential for EMFT
    real(kind=DP), allocatable, dimension(:,:,:,:) :: lhxc_fine_emft

    ! -------------------------------------------------------------------------

    call timer_clock('hamiltonian_lhxc_calculate',1)

    ! jme: KPOINTS_DANGER
    ! Some parts of this function enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine hamiltonian_lhxc_calculate not ready yet &
         &for more than one k-point.')

    ! rc2013: set no. of rows and columns in system
    mrows = overlap%mrows
    ncols = overlap%ncols

    if (.not.pub_xc_ke_density_required) then
       ! JCW: Non-tau-dependent XC functional
       ! JCW: Check dfdtau_fine is either not present, or has a length of 1
       if ( present(dfdtau_fine) ) then
          call utils_assert( &
               size(dfdtau_fine) == 1,&
               "Error in hamiltonian_lhxc_calculate: &
               &pub_xc_ke_density_required is false, but dfdtau_fine is present &
               &or has a size > 1.")
       end if
    else
       ! JCW: tau-dependent XC functional
       ! JCW: Check dfdtau_fine is present
       call utils_assert( present(dfdtau_fine),&
            "Error in hamiltonian_lhxc_calculate: &
            &pub_xc_ke_density_required is true, but dfdtau_fine is not &
            &present.")
    end if

    ! Optional argument to suppress xc potential being added to lhxc_fine
    if (present(add_xc_pot)) then
       loc_add_xc_pot = add_xc_pot
    else
       loc_add_xc_pot = .true.
    end if

    ! Allocate workspace
    fine_ld1 = mdl%fine_grid%ld1
    fine_ld2 = mdl%fine_grid%ld2
    fine_max_slabs12 = mdl%fine_grid%max_slabs12
    if(.not. present(return_density_fine)) then
       allocate(sub_density_fine(fine_ld1,fine_ld2,fine_max_slabs12, &
            pub_num_spins, mrows, ncols),stat=ierr)
       call utils_alloc_check('hamiltonian_lhxc_calculate','sub_density_fine',ierr)
    else
       sub_density_fine => return_density_fine
    end if
    allocate(sub_pot_fine(fine_ld1,fine_ld2,fine_max_slabs12, &
         pub_num_spins, mrows, ncols), stat=ierr)
    call utils_alloc_check('hamiltonian_lhxc_calculate','sub_pot_fine',ierr)
    if ((mrows.eq.1).and.(ncols.eq.1)) then
       density_fine => sub_density_fine(:,:,:,:,1,1)
       pot_fine => sub_pot_fine(:,:,:,:,1,1)
    else
       allocate(density_fine(fine_ld1,fine_ld2,fine_max_slabs12,pub_num_spins),&
            stat=ierr)
       call utils_alloc_check('hamiltonian_lhxc_calculate','density_fine',ierr)
       allocate(pot_fine(fine_ld1,fine_ld2,fine_max_slabs12,pub_num_spins), &
            stat=ierr)
       call utils_alloc_check('hamiltonian_lhxc_calculate','pot_fine',ierr)
    end if
    if (pub_aug) then
       allocate(nhat_den_grad(fine_ld1,fine_ld2,fine_max_slabs12,pub_num_spins,&
            0:pub_aug_den_dim),stat=ierr)
       call utils_alloc_check('hamiltonian_lhxc_calculate','nhat_den_grad',ierr)
    end if
    ! JCW: Allocate space for kinetic energy density
    ! (pub_xc_ke_density_required set in xc_init)
    if (pub_xc_ke_density_required) then
       allocate(sub_ke_density_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub, mdl%nsub), stat=ierr)
       call utils_alloc_check('hamiltonian_lhxc_calculate','sub_ke_density_fine',&
            ierr)
       if ((mrows.eq.1).and.(ncols.eq.1)) then
          ke_density_fine => sub_ke_density_fine(:,:,:,:,1,1)
       else
          allocate(ke_density_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
               mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
          call utils_alloc_check('hamiltonian_lhxc_calculate','ke_density_fine',&
               ierr)
       end if
    else
       ! JCW: Allocate a single element dummy array
       allocate(sub_ke_density_fine(1,1,1,1,1,1), stat=ierr)
       call utils_alloc_check('hamiltonian_lhxc_calculate','sub_ke_density_fine',&
            ierr)
       ke_density_fine => sub_ke_density_fine(:,:,:,:,1,1)
    end if

    ! rc2013: allocate space for density gradient and EMFT potential
    if(pub_emft) then
       allocate(lhxc_fine_emft(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
       call utils_alloc_check('hamiltonian_lhxc_calculate','lhxc_fine_emft',ierr)
       if (utils_devel_code(.false., 'HAM','EMFT_POT',pub_devel_code)) then
          allocate(full_density_grad(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
               mdl%fine_grid%max_slabs12, 3, pub_num_spins), stat=ierr)
          call utils_alloc_check('hamiltonian_lhxc_calculate','full_density_grad',ierr)
          allocate(emft_density_grad(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
               mdl%fine_grid%max_slabs12, 3, pub_num_spins), stat=ierr)
          call utils_alloc_check('hamiltonian_lhxc_calculate','emft_density_grad',ierr)
          allocate(demft_dfull(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
               mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
          !mdl%fine_grid%max_slabs12, 3, pub_num_spins), stat=ierr)
          call utils_alloc_check('hamiltonian_lhxc_calculate','demft_dfull',ierr)
          allocate(recip_work(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
               mdl%fine_grid%max_slabs12, 3), stat=ierr)
          call utils_alloc_check('hamiltonian_lhxc_calculate','recip_work',ierr)
       end if
    end if

    lhxc_fine = 0.0_DP
    pot_fine = 0.0_DP

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%% DENSITY ON GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! agrecokpt: need to replace denskern%m(:,PUB_1K) with full array
    ! denskern: sum over k-point component will be done automatically
    ! in density_on_grid routine

    ! rc2013: allocate the array to hold the DK components
    allocate(kern_array(pub_num_spins),stat=ierr)
    call utils_alloc_check('hamiltonian_lhxc_calculate','kern_array',ierr)

    sub_density_fine = 0.0_DP
    do jsub=1,ncols
       do isub=1,mrows
          ! rc2013: HACK! problems with deferred shape arrays using embedding structures.
          ! In practice this will probably require denskern to be added properly to
          ! density_on_grid (complex quantities).
          call sparse_embed_extract_from_array(kern_array,denskern%m(:,PUB_1K),&
               isub,jsub)
          call density_on_grid(sub_density_fine(:,:,:,:,isub,jsub), &   ! output
               mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
               kern_array, overlap%m(isub,jsub), &
               ngwfs_on_grid(isub), ngwf_basis(isub), &
               ngwfs_on_grid(jsub), ngwf_basis(jsub))

          ! ndmh: for PAW, get augmentation density nhat
          ! agrecokpt: need to replace denskern%m(:,PUB_1K) with full array
          ! with all k-points: sum over k-points will be done automatically
          ! in augmentation_density_on_grid routine
          ! rc2013: EMBED_FIX! not suitable for embedding structures!
          if (pub_aug) then
             nhat_den_grad = 0.0_DP
             call augmentation_density_on_grid(nhat_den_grad, &
                  mdl%fine_grid,mdl%cell,mdl%regions(isub)%pseudo_sp, &
                  mdl%regions(isub)%paw_sp,mdl%aug_box, &
                  kern_array, sp_overlap%m(isub,jsub))
             sub_density_fine(:,:,:,:,isub,jsub) = sub_density_fine(:,:,:,:,isub,jsub) &
                  + nhat_den_grad(:,:,:,:,0)
          end if

          !! JCW: if kinetic energy density required for XC, evaluate
          if (pub_xc_ke_density_required) then
             if (pub_debug_on_root) write(stdout,'(a)') &
                  "DEBUG: Evaluating kinetic energy density..."
             ! Evaluate kinetic energy density on grid
             ! rc2013: for now this is not ready to handle off-diagonal matrices
             if(isub == jsub) &
                  call ke_density_on_grid(sub_ke_density_fine(:,:,:,:,isub,isub), & ! output
                  mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
                  kern_array, overlap%m(isub,isub), &
                  ngwfs_on_grid(isub), ngwf_basis(isub))
             ! Also set dfdtau_fine to 0.0_DP
             dfdtau_fine = 0.0_DP
          end if

          ! rc2013: destroy kern_array
          call sparse_embed_destroy_extracted_array(kern_array)
       end do
    end do

    ! %%%%%%%%%%%%%%%%%%%%%%%%%% END DENSITY ON GRID %%%%%%%%%%%%%%%%%%%%%%%%%%

    ! jd: Write out current density to a file if a suitable devel code specified
    call internal_visualise_density

    ! =============== LOCAL PSEUDOPOTENTIAL ENERGY AND POTENTIAL ===============

    locpspot_energy = 0.0_DP
    do is=1,pub_num_spins
       ! rc2013: use full potential
       do jsub=1,ncols
          do isub=1,mrows
             locpspot_energy(isub,jsub) = locpspot_energy(isub,jsub) + &
                  integrals_product_on_grid( &
                  mdl%fine_grid,mdl%localpseudo_fine(:,:,:), &
                  sub_density_fine(:,:,:,is,isub,jsub))
          end do
          lhxc_fine(:,:,:,is,jsub) = mdl%localpseudo_fine(:,:,:)
       end do
    end do

    ! jd: Write out pot. component to a file if a suitable devel code specified
    ! rc2013: only use 1st (active) subsystem
    call internal_visualise_potential_and_efield(lhxc_fine(:,:,:,1,1),'locps',1)

    ! ============= END LOCAL PSEUDOPOTENTIAL ENERGY AND POTENTIAL =============

    ! ====================== HARTREE ENERGY AND POTENTIAL ======================

    ! jd: Calculate either:
    !     (1) the usual Hartree potential
    !     (2) Hartree potential with cutoff Coulomb
    !     (3) Hartree potential with multigrid, which could be
    !         a) electronic Hartree potential, with open BCs,
    !         b) molecular Hartree potential, with smeared-ions, in vacuo,
    !         c) molecular Hartree potential, in implicit solvent.
    !     All of the above are mutually exclusive.
    !     Iff (3c) and the dielectric is not fixed, but instead depends on
    !     the electronic density, an extra term caused by this dependence
    !     (V_eps(r)+V_apolar(r)) is included in the potential to get the correct
    !     derivative. The Hartree energy then needs to be
    !     calculated differently, and hartree_via_multigrid() takes care of
    !     that, as well as calculating the apolar energy to solvation.
    hartree_energy = 0.0_DP
    call zero_solvation_terms(solvation_terms, zero_apolar = .true., &
         zero_vacuum = .false.)

    if (.not. pub_turn_off_hartree) then
       if (pub_coulomb_cutoff) then
          do jsub=1,ncols
             do isub=1,mrows
                call cutoff_coulomb_hartree(sub_pot_fine(:,:,:,:,isub,jsub), &
                     sub_density_fine(:,:,:,:,isub,jsub), mdl)
             end do
          end do
       ! jcap: The multigrid Hartree calculation needs to be done
       ! for the whole system at once, to avoid complaints about
       ! non-integer charges
       else if (pub_multigrid_hartree) then
          ! jcap: calculate full density if we have multiple regions
          if (.not.((mrows.eq.1).and.(ncols.eq.1))) then
             density_fine = 0.0_DP
             do jsub=1,ncols
                do isub=1,mrows
                   density_fine = density_fine + &
                        sub_density_fine(:,:,:,:,isub,jsub)
                end do
             end do
          end if
          call hartree_via_multigrid(pot_fine, density_fine, &             ! (3)
               mdl%fine_grid, mdl%cell, hartree_energy(1,1), &
               solvation_terms, elements = mdl%elements)
          ! jcap: reset density_fine to 0 if required
          if (.not.((mrows.eq.1).and.(ncols.eq.1))) density_fine = 0.0_DP
       else
          do jsub=1,ncols
             do isub=1,mrows
                call hartree_on_grid(sub_pot_fine(:,:,:,:,isub,jsub), &
                     sub_density_fine(:,:,:,:,isub,jsub),  mdl%fine_grid, &   ! (1)
                     mdl%cell)
             end do
          end do
       end if
    else
       sub_pot_fine=0.0_DP
    endif

    do jsub=1,ncols
       do isub=1,mrows
          ! cks: scale Hartree potential for hyper Hartree-Fock
          if (pub_hhf_nstates > 0) then
             do is=1,pub_num_spins
                sub_pot_fine(:,:,:,is,isub,jsub) = &
                     pub_hhf_factor * sub_pot_fine(:,:,:,is,isub,jsub)
             enddo
             call utils_assert(.not. pub_multigrid_hartree, &
                  'Hyper Hartree-Fock (hhf_nstates) is currently incompatible with &
                  &multigrid Hartree calculation -- definition of hartree_energy &
                  &in hamiltonian_mod would have to be updated')
          endif
       end do
    end do

    ! jd: If working with the multigrid Hartree, hartree_energy is already
    !     calculated, no need to calculate hartree_energy directly.
    do jsub=1,ncols
       do isub=1,mrows
          if (.not.pub_multigrid_hartree) then
             do is=1,pub_num_spins
                ! rc2013: need to sum up all contributions from the Hartree potential
                do ll=1,ncols
                   do kk=1,mrows
                      hartree_energy(isub,jsub) = hartree_energy(isub,jsub) + &
                           0.5_DP*integrals_product_on_grid(mdl%fine_grid, &
                           sub_pot_fine(:,:,:,is,kk,ll), sub_density_fine(:,:,:,is,isub,jsub))
                   end do
                end do
             end do
             ! rc2013: add up the regional potentials
             ! jcap: This now only needs to be done if we aren't using
             ! multigrid
             if (.not.((mrows.eq.1).and.(ncols.eq.1))) then
                do is=1,pub_num_spins
                   pot_fine(:,:,:,is) = pot_fine(:,:,:,is) + &
                        sub_pot_fine(:,:,:,is,isub,jsub)
                end do
             end if
          end if
       end do
    end do

    ! rc2013: now add the full Hartree potential
    do is=1,pub_num_spins
       do isub=1,mrows
          lhxc_fine(:,:,:,is,isub) = lhxc_fine(:,:,:,is,isub) + pot_fine(:,:,:,is)
       end do
    end do

    ! jd: Write out pot. component to file if a suitable devel code specified
    call internal_visualise_potential_and_efield(pot_fine(:,:,:,1),'hartree',5)

    ! ==================== END HARTREE ENERGY AND POTENTIAL ====================

    ! ndmh: screen nonlocal projector energies now if the augmentation density
    ! ndmh: is not going to be included in the exchange energy calculation
    if (pub_paw.and.(.not.pub_nhat_in_xc)) then
       do isub=1,mdl%nsub
          call sparse_embed_extract_from_array(kern_array,dijhat,isub,isub)
          call augmentation_screen_dij(kern_array, &
               lhxc_fine(:,:,:,:,isub),mdl%aug_box,mdl%cell,mdl%fine_grid, &
               mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp)
          call sparse_embed_destroy_extracted_array(kern_array,dijhat,.true.,&
               isub,isub)
       end do
    end if

    ! ========================= XC ENERGY AND POTENTIAL ========================

    ! ndmh: add on core density before calculating xc potential for NLCC
    sub_pot_fine = 0.0_DP
    if((mrows .gt. 1) .or. (ncols .gt. 1)) then
       density_fine = 0.0_DP
       ke_density_fine = 0.0_DP
    end if
    if (pub_nlcc) then
       do is=1,pub_num_spins
          density_fine(:,:,:,is) = density_fine(:,:,:,is) + &
               mdl%core_density_fine * 0.5_DP * pub_spin_fac
       end do
    end if

    ! rc2013: get the total density for XC calculation if required
    if (.not.((mrows.eq.1).and.(ncols.eq.1))) then
       do jsub=1,ncols
          do isub=1,mrows
             do is=1,pub_num_spins
                density_fine(:,:,:,is) = density_fine(:,:,:,is) + &
                     sub_density_fine(:,:,:,is,isub,jsub)
                if (pub_xc_ke_density_required) &
                     ke_density_fine(:,:,:,is) = ke_density_fine(:,:,:,is) + &
                     sub_ke_density_fine(:,:,:,is,isub,jsub)
             end do
          end do
       end do
    endif

    ! ndmh: calculate XC energy and potential
    if (pub_aug) then
       ! JCW: Exit with error if KE density required, since KE density +
       ! augmentation is untested (temporary)
       if(pub_xc_ke_density_required) then
          call utils_abort('Error in hamiltonian_lhxc_calculate: &
               &KE density evaluation and augmentation not tested.')
       end if
       ! rc2013: also exit if EMFT is requested
       if(pub_emft) then
          call utils_abort('Error in hamiltonian_lhxc_calculate: &
               &augmentation not tested with EMFT.')
       end if
       ! ndmh: if we do not want the nhat density in the XC calculation,
       ! ndmh: remove it now and calculate xc potential without it
       if (pub_paw) then
          if (.not.pub_nhat_in_xc) then
             density_fine = density_fine - nhat_den_grad(:,:,:,:,0)
          end if
       end if
       call xc_energy_potential(density_fine, xc_energy, pot_fine, &
            mdl%fine_grid, mdl%cell, pub_aug_den_dim, nhat_den_grad)
    else
       call xc_energy_potential(density_fine, xc_energy, pot_fine, &
            mdl%fine_grid, mdl%cell, 0, &
            ke_density_fine = ke_density_fine, &
            dfdtau_fine = dfdtau_fine )

       ! jcap: If we are using embedded mean field theory, we need to
       ! call this twice more, in order to get the xc energy and
       ! potential for the active region
       if (pub_emft) then

          ! rc2013: calculate the gradient of the full density if required
          if (utils_devel_code(.false., 'HAM','EMFT_POT',pub_devel_code)) then
             call xc_gradients(density_fine, full_density_grad, recip_work, &
                  mdl%fine_grid, 0)
          end if

          ! jcap: If necessary, add on the nlcc charges
          density_fine = sub_density_fine(:,:,:,:,pub_active_region,pub_active_region)
          if (pub_nlcc) then
             do is=1,pub_num_spins
                density_fine(:,:,:,is) = density_fine(:,:,:,is) + &
                     mdl%regions(pub_active_region)%core_density_fine * 0.5_DP * pub_spin_fac
             end do
          end if

          ! rc2013: Now get the EMFT correction
          call xc_emft_calculate(density_fine, xc_energy_emft, &
               lhxc_fine_emft, mdl%fine_grid, mdl%cell)

          ! rc2013: add the EMFT correction to the XC energy
          xc_energy = xc_energy + xc_energy_emft

          ! rc2013: In principle we may need to correct the XC potential
          ! by accounting for the derivative of the total system density
          ! wrt the active density
          if (utils_devel_code(.false., 'HAM','EMFT_POT',pub_devel_code)) then
             call xc_gradients(density_fine, emft_density_grad, recip_work, &
                  mdl%fine_grid, 0)
             demft_dfull = 0.0_DP
             do is=1,pub_num_spins
                do islab12=1,mdl%fine_grid%num_my_slabs12
                   do i2=1,mdl%fine_grid%n2
                      do i1=1,mdl%fine_grid%n1
                         grad1(:) = full_density_grad(i1,i2,islab12,:,is)
                         grad2(:) = emft_density_grad(i1,i2,islab12,:,is)
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
             lhxc_fine_emft = lhxc_fine_emft*demft_dfull
          end if
       end if
    end if

    if (loc_add_xc_pot) then
       do is=1,pub_num_spins
          do isub=1,mdl%nsub
             lhxc_fine(:,:,:,is,isub) = lhxc_fine(:,:,:,is,isub) &
                  + pot_fine(:,:,:,is)
             ! rc2013: only add the EMFT potential to the active region
             if(isub == pub_active_region .and. pub_emft) &
                  lhxc_fine(:,:,:,is,isub) = lhxc_fine(:,:,:,is,isub) &
                  + lhxc_fine_emft(:,:,:,is)
          end do
       end do
       ! jd: Write out pot. component to file if a suitable devel code specified
       !call internal_visualise_potential(pot_fine(:,:,:,1),'xc',6)
    end if

    ! ====================== END XC ENERGY AND POTENTIAL =======================

    ! ndmh: screen nonlocal projector energies now if the augmentation density
    ! ndmh: is going to be included in the exchange energy calculation
    if (pub_paw.and.pub_nhat_in_xc) then
       do isub=1,mdl%nsub
          call sparse_embed_extract_from_array(kern_array,dijhat,isub,isub)
          call augmentation_screen_dij(kern_array, &
               lhxc_fine(:,:,:,:,isub),mdl%aug_box,mdl%cell,mdl%fine_grid, &
               mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp)
          call sparse_embed_destroy_extracted_array(kern_array,dijhat,.true.,&
               isub,isub)
       end do
    end if

    !==========================ELECTRONIC ENTHALPY==============================

    enthalpy_energy = 0.0_DP
    if (pub_external_pressure>0.0_DP) then
       call enthalpy_terms(pot_fine,enthalpy_energy,density_fine,mdl%fine_grid)
       do isub=1,mdl%nsub
          lhxc_fine(:,:,:,:,isub) = lhxc_fine(:,:,:,:,isub) + pot_fine(:,:,:,:)
       end do

       ! jd: Write out pot. component to file if a suitable devel code specified
       call internal_visualise_potential_and_efield(pot_fine(:,:,:,1),'enthalpy',7)
    end if

    ! ADD UP ENERGIES
    lhxc_energy = xc_energy + sum(locpspot_energy) + sum(hartree_energy) + &
         solvation_terms%E_apolar_cavitation + solvation_terms%E_apolar_disrep+&
         implicit_solvent_boltzmann_energy_contrib(solvation_terms) + &
         enthalpy_energy

    if(pub_debug_on_root) then
       write(stdout,'(a,f15.6)') 'Total xc_energy       = ', xc_energy
       write(stdout,'(a,f15.6)') 'Total locpspot_energy = ', sum(locpspot_energy)
       write(stdout,'(a,f15.6)') 'Total hartree_energy  = ', sum(hartree_energy)
       write(stdout,'(a,f15.6)') 'lhxc_energy           = ', lhxc_energy
       write(stdout,'(a,f15.6)') 'IS_apolar_energy      = ', &
            solvation_terms%E_apolar_cavitation
       write(stdout,'(a,f15.6)') 'IS_disrep_energy      = ', &
            solvation_terms%E_apolar_disrep
       write(stdout,'(a,f15.6)') 'IS_boltzmann_terms    = ', &
            implicit_solvent_boltzmann_energy_contrib(solvation_terms)
       write(stdout,'(a,f15.6)') 'enthalpy_energy       = ', enthalpy_energy
    endif

    ! Deallocate workspace
    ! rc2013: nullify and deallocate fine_grid pointers if necessary
    if ((mrows.eq.1).and.(ncols.eq.1)) then
       nullify(ke_density_fine)
       nullify(density_fine)
       nullify(pot_fine)
    else
       if (pub_xc_ke_density_required) then
          deallocate(ke_density_fine,stat=ierr)
          call utils_dealloc_check('hamiltonian_lhxc_calculate','ke_density_fine',ierr)
       end if
       deallocate(pot_fine,stat=ierr)
       call utils_dealloc_check('hamiltonian_lhxc_calculate','pot_fine',ierr)
       deallocate(density_fine,stat=ierr)
       call utils_dealloc_check('hamiltonian_lhxc_calculate','density_fine',ierr)
    end if
    ! JCW: Deallocate ke_density_fine. Always allocated since we allocate a
    ! JCW: single element dummy array if pub_xc_ke_density_required is false.
    deallocate(sub_ke_density_fine,stat=ierr)
    call utils_dealloc_check('hamiltonian_lhxc_calculate','sub_ke_density_fine',&
         ierr)
    if (pub_aug) then
       deallocate(nhat_den_grad,stat=ierr)
       call utils_dealloc_check('hamiltonian_lhxc_calculate','nhat_den_grad', &
            ierr)
    end if
    deallocate(sub_pot_fine,stat=ierr)
    call utils_dealloc_check('hamiltonian_lhxc_calculate','sub_pot_fine',ierr)
    if(pub_emft) then
       deallocate(lhxc_fine_emft,stat=ierr)
       call utils_dealloc_check('hamiltonian_lhxc_calculate','lhxc_fine_emft',ierr)
       if (utils_devel_code(.false., 'HAM','EMFT_POT',pub_devel_code)) then
          deallocate(full_density_grad,stat=ierr)
          call utils_dealloc_check('hamiltonian_lhxc_calculate','full_density_grad',ierr)
          deallocate(emft_density_grad,stat=ierr)
          call utils_dealloc_check('hamiltonian_lhxc_calculate','emft_density_grad',ierr)
          deallocate(demft_dfull,stat=ierr)
          call utils_dealloc_check('hamiltonian_lhxc_calculate','demft_dfull',ierr)
          deallocate(recip_work,stat=ierr)
          call utils_dealloc_check('hamiltonian_lhxc_calculate','recip_work',ierr)
       end if
    end if
    if(.not. present(return_density_fine)) then
       deallocate(sub_density_fine,stat=ierr)
       call utils_dealloc_check('hamiltonian_lhxc_calculate','sub_density_fine',ierr)
    end if
    deallocate(kern_array,stat=ierr)
    call utils_dealloc_check('hamiltonian_lhxc_calculate','kern_array',ierr)
    call timer_clock('hamiltonian_lhxc_calculate',2)

    return

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_visualise_density
      !=====================================================================!
      ! If 'DEVEL_CODE DEBUG: DUMP_DENSITY=T :DEBUG' is specified, writes   !
      ! the current density on the fine grid to a file (dx, cube, grd).     !
      ! This is useful for debugging runaway densities (eg. overpolarisation!
      ! due to charge sucking).                                             !
      !---------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in May 2015.                              !
      !---------------------------------------------------------------------!

      use constants, only: ANGSTROM, UP, DN
      use rundat, only: pub_devel_code, pub_inner_loop_iteration
      use utils, only: utils_devel_code
      use visual, only: visual_scalarfield

      implicit none

      ! jd: Local variables
      character(len=80) :: dens_output_filename_suffix
      integer, save :: dens_output_iter = 1
      logical, save :: dens_output_done = .false.

      ! -----------------------------------------------------------------------

      if(utils_devel_code(.false.,'DEBUG','DUMP_DENSITY', pub_devel_code, &
           no_warn = .true.)) then
         if(pub_inner_loop_iteration == 1 .and. .not. dens_output_done) then
            write(dens_output_filename_suffix,'(a,i0)') '_density_', &
                 dens_output_iter
            if(pub_num_spins == 1) then
               call visual_scalarfield( &
                    density_fine(:,:,:,1), mdl%fine_grid, mdl%cell, &
                    'Electronic density (in e/ang^3)', &
                    dens_output_filename_suffix, mdl%elements, ANGSTROM**3)
            else
               call visual_scalarfield( &
                    density_fine(:,:,:,UP) + density_fine(:,:,:,DN), &
                    mdl%fine_grid, mdl%cell, &
                    'Electronic density (in e/ang^3)', &
                    dens_output_filename_suffix, mdl%elements, ANGSTROM**3)
            end if
            dens_output_iter = dens_output_iter + 1
            dens_output_done = .true.
         else
            dens_output_done = .false.
         end if
      end if

    end subroutine internal_visualise_density

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_visualise_potential_and_efield(&
         pot_array, pot_string, pot_id)
      !=====================================================================!
      ! If 'DEVEL_CODE DEBUG: DUMP_POTENTIAL=T :DEBUG' is specified, writes !
      ! a given potential on the fine grid to a file (dx, cube, grd).       !
      ! This is useful for debugging runaway potentials (eg. overpolari-    !
      ! sation due to charge sucking).                                      !
      !                                                                     !
      ! Similarly, if 'DEVEL_CODE DEBUG: DUMP_EFIELD_{PBC}=T :DEBUG' is     !
      ! specified, calculates and dumps the electric field due to the given !
      ! potential on the fine grid to a file (dx, cube, grd). For PBC the   !
      ! electric fields are calculated in reciprocal space, with the usual  !
      ! caveats (ringing). For OBC the electric fields are calculated with  !
      ! finite differences, with an order 'pub_finite_difference_order'.    !
      ! Both calculations can be performed at the same time.                !
      !                                                                     !
      ! To save time and disk space we only output in the 1st inner loop    !
      ! iteration and in properties. We do not output in subsequent inner   !
      ! loop iterations or in the NGWF gradient stage. Output files are     !
      ! numbered sequentially. The index is incremented after each call.    !
      !---------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in October 2016.                          !
      !---------------------------------------------------------------------!

      use constants, only: HARTREE_PER_BOHR_TO_MV_PER_CM
      use rundat, only: pub_devel_code, pub_finite_difference_order, &
           pub_inner_loop_iteration
      use utils, only: utils_devel_code
      use visual, only: visual_scalarfield, visual_efield_on_grid_obc, &
           visual_efield_on_grid_pbc

      implicit none

      ! jd: Arguments
      real(kind=DP), intent(in)    :: pot_array(:,:,:)
      character(len=*), intent(in) :: pot_string
      integer, intent(in)          :: pot_id

      ! jd: Local variables
      real(kind=DP), allocatable :: efield_fine(:,:,:,:)
      character(len=80) :: output_filename_suffix
      logical :: dump_efield_obc, dump_efield_pbc
      integer, save :: output_iter(10) = (/1,1,1,1,1,1,1,1,1,1/)
      logical, save :: pot_output_done(10) = (/.false.,.false.,.false.,.false.,&
           .false.,.false.,.false.,.false.,.false.,.false./)

      ! -----------------------------------------------------------------------

      ! NB: Output
      if((pub_inner_loop_iteration == 0 .or. pub_inner_loop_iteration == 1) &
           .and. .not. pot_output_done(pot_id)) then
         if(utils_devel_code(.false.,'DEBUG','DUMP_POTENTIAL',pub_devel_code, &
              no_warn = .true.)) then
            write(output_filename_suffix,'(a,a,a,i0)') '_potential_', &
                 pot_string, '_', output_iter(pot_id)
            call visual_scalarfield(pot_array, mdl%fine_grid, mdl%cell, &
                 'Potential (internal units)', &
                 output_filename_suffix, mdl%elements)
         end if

         dump_efield_obc = utils_devel_code(.false.,'DEBUG','DUMP_EFIELD_OBC', &
              pub_devel_code)
         dump_efield_pbc = utils_devel_code(.false.,'DEBUG','DUMP_EFIELD_PBC', &
              pub_devel_code)

         if(dump_efield_obc .or. dump_efield_pbc) then
            allocate(efield_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
                 mdl%fine_grid%max_slabs12, 3), stat=ierr)
            call utils_alloc_check('internal_visualise_potential &
                 &(hamiltonian_lhxc_calculate)','efield_fine', ierr)

            if(dump_efield_obc) then
               write(output_filename_suffix,'(a,a,a,i0)') '_efield_obc_', &
                    pot_string, '_', output_iter(pot_id)
               call visual_efield_on_grid_obc(efield_fine, pot_array, &
                    mdl%fine_grid, mdl%cell, pub_finite_difference_order)
               call visual_scalarfield(efield_fine(:,:,:,1), &
                    mdl%fine_grid, mdl%cell, &
                    'Electric field '//pot_string//' [MV/cm] for:', &
                    trim(output_filename_suffix)//'_x', mdl%elements, &
                    -HARTREE_PER_BOHR_TO_MV_PER_CM)
               call visual_scalarfield(efield_fine(:,:,:,2), &
                    mdl%fine_grid, mdl%cell, &
                    'Electric field '//pot_string//' [MV/cm] for:', &
                    trim(output_filename_suffix)//'_y', mdl%elements, &
                    -HARTREE_PER_BOHR_TO_MV_PER_CM)
               call visual_scalarfield(efield_fine(:,:,:,3), &
                    mdl%fine_grid, mdl%cell, &
                    'Electric field '//pot_string//' [MV/cm] for:', &
                    trim(output_filename_suffix)//'_z', mdl%elements, &
                    -HARTREE_PER_BOHR_TO_MV_PER_CM)
            end if
            if(dump_efield_pbc) then
               write(output_filename_suffix,'(a,a,a,i0)') '_efield_pbc_', &
                    pot_string, '_', output_iter(pot_id)
               call visual_efield_on_grid_pbc(efield_fine, pot_array, &
                    mdl%fine_grid, mdl%cell)
               call visual_scalarfield(efield_fine(:,:,:,1), &
                    mdl%fine_grid, mdl%cell, &
                    'Electric field '//pot_string//' [MV/cm] for:', &
                    trim(output_filename_suffix)//'_x', mdl%elements, &
                    -HARTREE_PER_BOHR_TO_MV_PER_CM)
               call visual_scalarfield(efield_fine(:,:,:,2), &
                    mdl%fine_grid, mdl%cell, &
                    'Electric field '//pot_string//' [MV/cm] for:', &
                    trim(output_filename_suffix)//'_y', mdl%elements, &
                    -HARTREE_PER_BOHR_TO_MV_PER_CM)
               call visual_scalarfield(efield_fine(:,:,:,3), &
                    mdl%fine_grid, mdl%cell, &
                    'Electric field '//pot_string//' [MV/cm] for:', &
                    trim(output_filename_suffix)//'_z', mdl%elements, &
                    -HARTREE_PER_BOHR_TO_MV_PER_CM)
            end if

            deallocate(efield_fine, stat=ierr)
            call utils_dealloc_check('internal_visualise_potential &
                 &(hamiltonian_lhxc_calculate)','efield_fine', ierr)
         end if

         output_iter(pot_id) = output_iter(pot_id) + 1
         pot_output_done(pot_id) = .true.
      else
         pot_output_done(pot_id) = .false.
      end if

    end subroutine internal_visualise_potential_and_efield

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine hamiltonian_lhxc_calculate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hamiltonian_species_scissor(scissor_ham, denskern, rep, ngwf_basis, &
       mdl)

    !========================================================================!
    ! This subroutine gets the scissor shifted Hamiltonian:                  !
    ! H_scissor = S(vK + c(S^-1-K))S                                         !
    !------------------------------------------------------------------------!
    ! scissor_ham        (inout) : Hamiltonian                               !
    ! denskern              (in) : Density kernel                            !
    ! ngwf_basis            (in) : NGWF basis                                !
    ! mdl                   (in) : Model                                     !
    ! rep                   (in) : NGWF representation                       !
    !------------------------------------------------------------------------!
    ! Written by Nelson Yeung in January 2019.                               !
    !========================================================================!

    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_num_spins, pub_spin_fac, PUB_1K
    use scissor, only: scissor_shift_denskern
    use sparse, only: SPAM3, sparse_create, sparse_product, sparse_destroy
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
         sparse_embed_extract_from_array, sparse_embed_scale
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout)       :: scissor_ham(pub_num_spins)
    type(SPAM3_EMBED), intent(inout)       :: denskern(pub_num_spins)
    type(NGWF_REP), intent(in)             :: rep
    type(FUNC_BASIS), intent(in)           :: ngwf_basis(rep%nsub)
    type(MODEL), intent(in)                :: mdl

    ! Local variables
    integer                        :: ierr
    integer                        :: is
    integer                        :: isub, jsub
    type(SPAM3)                    :: shifted_denskern, sk

    ! Add spin degeneracy to kernel
    if (pub_num_spins /= 2) then
       call sparse_embed_scale(denskern(1), 1.0_DP/pub_spin_fac)
    end if

    do isub=1,mdl%nsub
       do jsub=1,mdl%nsub
          call sparse_create(shifted_denskern,denskern(1)%m(isub,jsub))
          call sparse_create(sk,rep%overlap%m(isub,jsub),denskern(1)%m(isub,jsub))

          do is=1,pub_num_spins
             ! H = S(vK + c(S^-1-K))S
             call scissor_shift_denskern(shifted_denskern,denskern(is)%m(isub,jsub), &
                  ngwf_basis(isub),mdl%regions(isub), &
                  rep%inv_overlap%m(isub,jsub))
             call sparse_product(sk,rep%overlap%m(isub,jsub),shifted_denskern)
             call sparse_product(scissor_ham(is)%m(isub,jsub),sk, &
                  rep%overlap%m(isub,jsub))
          enddo

          call sparse_destroy(shifted_denskern)
          call sparse_destroy(sk)
       enddo
    enddo

    ! Remove spin degeneracy
    if (pub_num_spins /= 2) then
       call sparse_embed_scale(denskern(1), pub_spin_fac)
    end if

  end subroutine hamiltonian_species_scissor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hamiltonian_build_matrix(ham,rep)

    !========================================================================!
    ! This subroutine builds the full Hamiltonian matrix in SPAM3 format out !
    ! of its various components, which have already been calculated as SPAM3 !
    ! matrices.                                                              !
    !------------------------------------------------------------------------!
    ! ham        (inout) : Hamiltonian Wrapper (contains lhxc,hubbard_ham,   !
    !                      hfexchange and PAW nonlocpot)                     !
    ! rep           (in) : NGWF Representation (contains kinet and normal    !
    !                      nonlocpot)                                        !
    !------------------------------------------------------------------------!
    ! Written by Nicholas Hine in April 2010.                                !
    ! Modified by James C. Womack to add contribution from matrix elements   !
    ! over tau-dependent part of XC potential when a tau-dependent XC        !
    ! functional is used, 2015                                               !
    ! Adapted for embedding structures by Robert Charlton, August 2017.      !
    !========================================================================!

    use constants, only: DP, stdout
    use comms, only: pub_on_root ! JCW: For debugging
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_any_nl_proj, pub_use_hfx, pub_hubbard, &
         pub_aug, pub_usp, pub_num_spins, pub_confined_ngwfs, &
         pub_pol_emb_qmstar, pub_eda_mode, pub_eda_scfmi, pub_debug, &
         pub_xc_ke_density_required, pub_use_activehfx, pub_scissor_ngroups, &
         pub_emft
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_axpy, sparse_embed_copy, &
        sparse_embed_trace ! JCW: For debugging
    use constants, only: EDA_POLFRAGLOC_DEVEL_OFF, EDA_POLSIMUL_OFF, &
         EDA_CTFRAGLOC_DEVEL_OFF
    use fragment_scfmi, only: scfmi_hamiltonian_projection
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(inout) :: ham
    type(NGWF_REP), intent(in) :: rep

    ! Local variable
    integer :: is    ! spin loop counter
    real(kind=DP) :: trace ! JCW: For debugging

    ! -------------------------------------------------------------------------

    do is=1,pub_num_spins

       ! ndmh: copy lhxc matrix to Hamiltonian
       call sparse_embed_copy(ham%ham(is),ham%lhxc(is))

       if (pub_debug) then
          ! JCW: Output trace of matrix for debugging XC functionals
          call sparse_embed_trace(trace, ham%lhxc(is))
          if (pub_on_root) then
               write(stdout,'(a,es20.10)') "DEBUG: tr[lhxc(is)] = ", trace
          end if
       end if

       if (pub_xc_ke_density_required) then
          call utils_assert(allocated(ham%dfdtau),"Error in &
               &hamiltonian_build_matrix: tau-dependent XC functional &
               &requested but ham%dfdtau not allocated.")
          ! JCW: add dfdtau matrix (correction to lhxc for tau-dependent
          ! JCW: functionals, when XC potential is evaluated using
          ! JCW: functional derivatives wrt orbitals, i.e. FDO method).
          if (pub_debug_on_root) then
             write(stdout,'(a)') "DEBUG: Adding dfdtau matrix to Hamiltonian"
          end if
          call sparse_embed_axpy(ham%ham(is),ham%dfdtau(is),1.0_DP)
          if (pub_debug) then
             ! JCW: Output trace of matrix for debugging XC functionals
             call sparse_embed_trace(trace, ham%dfdtau(is))
             if (pub_on_root) then
                write(stdout,'(a,es20.10)') "DEBUG: tr[dfdtau(is)] = ", trace
             end if
          end if
       end if

       ! nmdh: add kinetic energy matrix
       call sparse_embed_axpy(ham%ham(is),rep%kinet,1.0_DP)

       ! ndmh: add the nonlocal potential matrix if necessary
       if (pub_any_nl_proj.and.(.not.pub_usp)) &
            call sparse_embed_axpy(ham%ham(is),rep%nonlocpot,1.0_DP)
       if (pub_aug) &
            call sparse_embed_axpy(ham%ham(is),ham%nonlocpot(is),1.0_DP)

       ! ny: add the scissor shifted matrix if necessary
       if (pub_scissor_ngroups > 0) &
            call sparse_embed_axpy(ham%ham(is),ham%scissor(is),1.0_DP)

       ! qoh: add the hfexchange matrix if necessary
       if (pub_use_hfx.or.(pub_use_activehfx.and.pub_emft)) &
            call sparse_embed_axpy(ham%ham(is),ham%hfexchange(is),-1.0_DP)

       ! jd: add the polemb matrix if necessary
       if (pub_pol_emb_qmstar) then
          call sparse_embed_axpy(ham%ham(is),ham%polemb(1),1.0_DP) ! jd: induced
          call sparse_embed_axpy(ham%ham(is),ham%polemb(2),1.0_DP) ! jd: permanent
          ! jd: NB the vacuum_polemb and vacuum_aux_polemb matrices do NOT get
          !     added, since they do not contribute to tr[KH].
       end if

       ! ddor: add the hubbard hamiltonian if necessary
       if (pub_hubbard) call sparse_embed_axpy(ham%ham(is),ham%hubbard_ham(is),1.0_DP)

       ! gcc32: add NGWF confinement hamiltonian if necessary
       if (pub_confined_ngwfs) call sparse_embed_axpy(ham%ham(is),ham%confinement(is), &
            1.0_DP)

    end do

    ! mjsp: Electronic polarisation via SCF-MI:
    if (pub_eda_scfmi) then

       ! mjsp:     Simultanous polarisation    (pub_eda_mode = EDA_POLSIMUL)
       ! mjsp:        Fragment polarisation    (pub_eda_mode = EDA_POLFRAGLOC_DEVEL)
       ! mjsp: (Inter)Fragment charge transfer (pub_eda_mode = EDA_CTFRAGLOC_DEVEL)

       ! mjsp: During NGWF gradient calculation, the full and normalised
       ! density representation is used and so projection is unwarranted.
       if ((pub_eda_mode /= EDA_POLSIMUL_OFF) .and. &
           (pub_eda_mode /= EDA_POLFRAGLOC_DEVEL_OFF) .and. &
           (pub_eda_mode /= EDA_CTFRAGLOC_DEVEL_OFF)) then

          ! mjsp: SCF-MI projection equations
          call scfmi_hamiltonian_projection(ham, rep)

       end if

    end if

  end subroutine hamiltonian_build_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hamiltonian_proj_cond_matrix(cond_rep, cond_ham, &
       val_rep, val_ham, val_dkn, cond_dkn, shift_changed)

    !=========================================================================!
    ! This subroutine calculates the projected conduction Hamiltonian matrix  !
    ! for a given shift w such that:                                          !
    ! H_proj = H_c - wS_c - T+K_vH_vK_vT + wT+K_vS_vK_vT                      !
    ! See L.E. Ratcliff, N.D.M. Hine and P.D. Haynes, Phys. Rev. B 84, 165131 !
    ! (2011) for details, specifically Sections II.A to II.D                  !
    !-------------------------------------------------------------------------!
    ! cond_rep       (in) : conduction NGWF representation                    !
    ! val_rep        (in) : valence NGWF representation                       !
    ! cond_ham    (inout) : conduction Hamiltonian matrices                   !
    ! val_ham        (in) : conduction Hamiltonian matrices                   !
    ! val_kernel     (in) : valence density kernel                            !
    ! shift_changed (out) : optional logical output for changing shifts       !
    !-------------------------------------------------------------------------!
    ! Written by Laura Ratcliff in June 2010.                                 !
    ! Modified by Laura Ratcliff in October 2010 to include NGWF_REP and      !
    ! NGWF_HAM types.                                                         !
    ! Modified by Nicholas Hine in January 2012 to add comments and ensure    !
    ! proper working with spin polarisation.                                  !
    ! Modified for embedding by Joseph Prentice, June 2018                    !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout, VERBOSE, max_spins
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_output_detail, cond_fixed_shift, cond_calc_max_eigen, &
         cond_shift_buffer, pub_debug, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_axpy, sparse_embed_copy, &
         sparse_embed_create, sparse_embed_destroy, sparse_embed_extremal_eigenvalue, &
         sparse_embed_product, sparse_embed_transpose, sparse_embed_trace, &
         sparse_embed_transpose_structure

    implicit none

    ! Arguments
    type(NGWF_REP), intent(in)     :: cond_rep
    type(NGWF_REP), intent(in)     :: val_rep
    type(NGWF_HAM), intent(inout)  :: cond_ham
    type(NGWF_HAM), intent(in)     :: val_ham
    type(SPAM3_EMBED), intent(in)        :: val_dkn(pub_num_spins)
    type(SPAM3_EMBED), optional, intent(in)        :: cond_dkn(pub_num_spins)
    logical, optional, intent(out) :: shift_changed

    ! Local variables
    integer :: is    ! spin loop counter
    real(kind=DP),parameter :: delta_e_thresh =1.0E-12_DP ! energy/atom threshold
    real(kind=DP) :: max_en(max_spins)      ! maximum hamiltonian eigenvalue

    type(SPAM3_EMBED) :: sinv_ham         ! inverse overlap times hamiltonian
    type(SPAM3_EMBED) :: trans_cross_olap ! the transpose of the cross overlap
    type(SPAM3_EMBED) :: uk, kt, ukh, uks, ukhkt, ukskt ! matrix products

    real(kind=DP) :: debug_trace    ! matrix trace for debugging purposes
    ! agrecocmplx
    logical :: loc_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
         &Entering hamiltonian_proj_cond_matrix'

    max_en=0.0_DP

    ! agrecocmplx
    loc_cmplx = val_dkn(1)%p%iscmplx

    if (present(shift_changed)) shift_changed = .false.

    ! Create transpose of cross_overlap matrix
    call sparse_embed_transpose_structure(trans_cross_olap%structure, &
         cond_rep%cross_overlap)
    call sparse_embed_create(trans_cross_olap,iscmplx=loc_cmplx)
    call sparse_embed_transpose(trans_cross_olap,cond_rep%cross_overlap)

    do is=1,pub_num_spins

       ! lr408: copy unprojected ham to cond_non_proj_ham
       call sparse_embed_copy(cond_ham%cond_non_proj_ham(is),cond_ham%ham(is))

       ! lr408: If required, calculate maximum eigenvalue of the conduction
       ! lr408: Hamiltonian
       if (cond_calc_max_eigen) then

          call sparse_embed_create(sinv_ham, cond_rep%inv_overlap, &
               cond_ham%cond_non_proj_ham(is))

          ! lr408: S^-1.H
          call sparse_embed_product(sinv_ham, cond_rep%inv_overlap, &
               cond_ham%cond_non_proj_ham(is))

          ! lr408: maximum orbital energy
          call sparse_embed_extremal_eigenvalue(sinv_ham, cond_rep%overlap, &
               max_en(is), delta_e_thresh)

          if ((pub_output_detail >= VERBOSE).and.(pub_on_root)) then
             write(stdout,'(/a,f20.14,f20.14)') &
                  'Conduction Hamiltonian max eigenvalue, shift:', max_en(is), &
                  cond_ham%cond_shift
          end if

          call sparse_embed_destroy(sinv_ham)

       end if

    end do

    ! lr408: If the shift is not being kept constant and the maximum
    ! lr408: eigenvalue has gone above the current shift, increase
    ! lr408: the shift accordingly
    if (.not. cond_fixed_shift) then
       if (maxval(max_en) > cond_ham%cond_shift) then
          cond_ham%cond_shift = maxval(max_en) + cond_shift_buffer
          if (pub_output_detail >= VERBOSE) then
             if (pub_on_root) write(stdout,'(/a,f24.14)') &
                  'Conduction shift updated:', cond_ham%cond_shift
          end if
          if (present(shift_changed)) shift_changed = .true.
       end if
    end if

    ! Allocate storage for matrices required for projection
    call sparse_embed_create(uk,trans_cross_olap,val_dkn(1))
    call sparse_embed_create(kt,val_dkn(1),cond_rep%cross_overlap)
    call sparse_embed_create(ukh,uk,val_ham%ham(1))
    call sparse_embed_create(uks,uk,val_rep%overlap)
    call sparse_embed_create(ukhkt,ukh,kt)
    call sparse_embed_create(ukskt,uks,kt)

    ! Now project out valence Hamiltonian and shift valence states up in energy
    ! (see Eq.11 of PRB 84 165131 (2011))
    do is=1,pub_num_spins

       ! Calculate ukhkt and ukhks
       call sparse_embed_product(uk,trans_cross_olap,val_dkn(is))
       call sparse_embed_product(kt,val_dkn(is),cond_rep%cross_overlap)
       call sparse_embed_product(ukh,uk,val_ham%ham(is))
       call sparse_embed_product(uks,uk,val_rep%overlap)
       call sparse_embed_product(ukhkt,ukh,kt)
       call sparse_embed_product(ukskt,uks,kt)

       ! lr408: add -TKHKT to conduction hamiltonian
       call sparse_embed_axpy(cond_ham%ham(is),ukhkt,-1.0_dp)

       ! lr408: add wTKSKT to conduction hamiltonian
       call sparse_embed_axpy(cond_ham%ham(is),ukskt,cond_ham%cond_shift)

       ! lr408: Calculate Tr[T+KSKTM] as a measure of orthogonality with valence
       ! lr408: states
       if (present(cond_dkn) .and. pub_debug) then
          call sparse_embed_trace(debug_trace,ukskt,cond_dkn(is))
          if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
               &Tr[T+KS_vKTM] = ',debug_trace
          call sparse_embed_trace(debug_trace,ukskt,cond_rep%inv_overlap)
          if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
               &Tr[T+KS_vKTS^-1] = ',debug_trace
       end if

    end do

    ! Clean up temporary matrices
    call sparse_embed_destroy(ukskt)
    call sparse_embed_destroy(ukhkt)
    call sparse_embed_destroy(uks)
    call sparse_embed_destroy(ukh)
    call sparse_embed_destroy(kt)
    call sparse_embed_destroy(uk)

    call sparse_embed_destroy(trans_cross_olap)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
         &Leaving hamiltonian_proj_cond_matrix'

  end subroutine hamiltonian_proj_cond_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! agrecokpt: do we need to include k-point dependence for TB method?
  subroutine hamiltonian_energy_components(pur_denskern, &    ! in/out
       rep, mdl, ngwf_basis, hub_proj_basis, &                ! input
       hub, &                                                 ! in/out
       hfexchange, mu, edft)                                  ! input

    !=====================================================================!
    ! This subroutine calculates and prints the components that make      !
    ! up the total energy.                                                !
    !---------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 17/4/2001 as                    !
    ! hamiltonian_density_ded_debug and hamiltonian_energy_terms.         !
    ! Modified by Peter Haynes 1/7/2004 for Fourier parallelisation.      !
    ! Modified by Chris-Kriton Skylaris on 7/10/2004 to return the        !
    ! Hartree energy as calculated directly from the charge density.      !
    ! Renamed to hamiltonian_energy_components and rewritten by           !
    ! Chris-Kriton Skylaris on 8/10/2004 in order to print more details   !
    ! using less memory.                                                  !
    ! Converted to SPAM2 by Peter Haynes, July 2006                       !
    ! Added HF exchange (optional) by Quintin Hill 23/10/2008.            !
    ! Modified by Nicholas Hine on 6/2/2009 for NLCC core charges         !
    ! Modified for DFT+U by D. D. O'Regan in April 2009                   !
    ! Modified by Nicholas Hine in July 2009 for function basis type      !
    ! Modified by Jacek Dziedzic in May 2015 to include local polarisable !
    ! embedding potential.                                                !
    ! Modified for embedding by Robert Charlton, 27/07/17.                !
    !=====================================================================!

    use augmentation, only: augmentation_density_on_grid, &
         aug_projector_denskern, augmentation_screen_dij, aug_nonlocal_mat
    use comms, only: comms_reduce, pub_on_root
    use constants, only: DP, UP, DN, stdout, max_spins, VERBOSE, paw_en_size, &
         paw_en_dij0, paw_en_ehart, paw_en_exc, paw_en_etxc, paw_en_exc_core, &
         k_B
    use cutoff_coulomb, only: cutoff_coulomb_hartree
    use density, only: density_on_grid
    use ensemble_dft_type, only: EDFT_MODEL
    use function_basis, only: FUNC_BASIS
    use hartree, only: hartree_on_grid, hartree_via_multigrid
    use hubbard_build, only: HUBBARD_MODEL, hubbard_energy_total
    use integrals, only: integrals_trace_on_grid, integrals_product_on_grid, &
         integrals_locpot
    use is_poisson, only: SOLVATION_ENERGY_TERMS, zero_solvation_terms
    use is_smeared_ions, only: smeared_ion_E_self, &
         smeared_ion_E_smeared
    use is_solvation_boltzmann, only: implicit_solvent_boltzmann_energy_contrib
    use ke_density, only: ke_density_on_grid, ke_density_kinetic_energy
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use polarisable_embedding, only: polarisable_embedding_calculate
    use rundat, only: pub_debug, pub_edft, &
         pub_output_detail, pub_dispersion, pub_use_hfx, &
         pub_coulomb_cutoff, pub_hubbard, pub_spin, pub_any_nl_proj, &
         pub_nlcc, pub_paw, &
         pub_pspot_bc_is_periodic, pub_ion_ion_bc_is_periodic, &
         pub_is_smeared_ion_rep, pub_aug, pub_usp, pub_is_pbe, &
         pub_is_include_apolar, pub_multigrid_hartree, &
         pub_multigrid_bc_is_periodic, pub_aug_den_dim, &
         pub_nhat_in_xc, pub_cdft, pub_num_spins, pub_turn_off_hartree, &
         pub_hubbard_compute_u_or_j, pub_spin_fac, pub_pol_emb_pot, pub_1K, &
         pub_eda_scfmi, pub_charge, pub_xc_ke_density_required, &
         pub_hhf_nstates, pub_hhf_factor, pub_scissor_ngroups, &
         pub_emft, pub_active_region, pub_use_activehfx, pub_inner_loop_iteration, &
         !ep
         pub_mermin, pub_mermin_smearing_width, pub_mermin_cheb, PUB_1K, &
         pub_num_kpoints

    use sparse, only: SPAM3, sparse_trace, sparse_create, &
         sparse_destroy, sparse_scale, sparse_copy
    use sparse_array, only: SPAM3_ARRAY, sparse_array_create, &
         sparse_array_destroy
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_trace, sparse_embed_scale, sparse_embed_copy, SPAM3_EMBED_ARRAY, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    !ep
    use smearing_operator, only: calculate_fermi_entropy_mermin
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, &
         utils_assert
    use vdwcorrection, only: pub_dispersion_energy
    use xc, only: xc_energy_potential,xc_init,xc_exit,xc_hfxinit,xc_emft_calculate

    implicit none

    ! Arguments
    type(NGWF_REP), intent(in) :: rep
    type(SPAM3_EMBED), intent(inout) :: pur_denskern(:) ! jd: inout because of scaling
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: hub_proj_basis(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(MODEL), intent(in) :: mdl
    type(SPAM3_EMBED), intent(in) :: hfexchange(:) ! Hartree-Fock exchange matrix
    real(kind=DP), optional, intent(in) :: mu(:)  !mu for MERMIN
    type(EDFT_MODEL), intent(in), optional :: edft ! jd: for printing entropic
                                                   !     contrib. & free energy

    ! Local Variables
    integer :: is
    integer :: ierr
    integer :: i1,i2,islab12                       ! Grid counters
    integer :: fine_max_slabs12, fine_ld1, fine_ld2 ! jd: shortcut notation
    real(kind=DP), allocatable, dimension(:,:,:,:) :: density_fine
    ! JCW: evaluate ke_density for use in meta-GGAs, if needed
    real(kind=DP), allocatable, dimension(:,:,:,:) :: ke_density_fine
    ! JCW: evaluate gradient of XC energy per unit volume wrt tau, if needed
    real(kind=DP), allocatable, dimension(:,:,:,:) :: dfdtau_fine
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: nhat_den_grad
    real(kind=DP), allocatable, dimension(:,:,:,:) :: buffer_fine
    ! rc2013: storage for densities when doing embedding
    real(kind=DP), allocatable, dimension(:,:,:,:,:,:) :: sub_density_fine
    real(kind=DP), allocatable, dimension(:,:,:,:,:,:) :: sub_buffer_fine
    real(kind=DP), allocatable, dimension(:,:,:,:,:,:) :: sub_ke_density_fine
    real(kind=DP), allocatable, dimension(:,:,:,:) :: lhxc_fine
    real(kind=DP) :: integrated_density
    real(kind=DP) :: integrated_ps_val_density
    real(kind=DP) :: integrated_spin
    real(kind=DP) :: integrated_field(max_spins)
    real(kind=DP) :: integrated_nhat(max_spins)
    real(kind=DP) :: tr_integrated_density
    real(kind=DP) :: tr_integrated_spin
    real(kind=DP) :: integrated_mod_spin
    real(kind=DP) :: integrated_depleted_spin
    real(kind=DP) :: tr_field(max_spins)
    real(kind=DP) :: hartree_energy(mdl%nsub,mdl%nsub)
    real(kind=DP) :: nhat_hartree_energy
    real(kind=DP) :: tr_hartree_energy
    real(kind=DP) :: tr_nhat_hartree_energy
    real(kind=DP) :: total_hartree_energy  ! rc2013: test of subsystem Hartree
    real(kind=DP) :: locpot_energy
    real(kind=DP) :: nhat_locpot_energy
    real(kind=DP) :: tr_locpot_energy
    real(kind=DP) :: tr_nhat_locpot_energy
    real(kind=DP) :: locps_gzero_energy
    real(kind=DP) :: kinetic_energy
    real(kind=DP) :: kedens_kinetic_energy ! JCW: to test KE density evaluation
    real(kind=DP) :: nonlocpot_energy
    type(SOLVATION_ENERGY_TERMS) :: solvation_terms
    real(kind=DP) :: scissor_energy
    real(kind=DP) :: xc_energy
    real(kind=DP) :: xc_denpot_energy
    real(kind=DP) :: xc_nhat_denpot_energy
    real(kind=DP) :: hfx_energy ! Hartree-Fock exchange energy
    real(kind=DP) :: hubbard_energy ! DFT+U correction energy
    real(kind=DP) :: hubbard_e_alpha ! DFT+alpha correction energy
    real(kind=DP) :: hubbard_e_spin  ! DFT+spin  correction energy
    real(kind=DP) :: spin_squared_expectation
    real(kind=DP) :: factor
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    real(kind=DP) :: total_energy
    real(kind=DP) :: total_sphere_energy
    real(kind=DP) :: trace
    type(SPAM3_EMBED) :: localmat
    type(SPAM3_EMBED), allocatable :: projector_denskern(:)
    type(SPAM3_EMBED), allocatable :: local_nonlocpot(:)
    type(SPAM3_EMBED), allocatable :: local_dijhat(:)
    type(SPAM3_EMBED), allocatable :: local_scissor_ham(:)
    real(kind=DP) :: smeared_ion_self_energy    ! jd: Extra terms resulting from
    real(kind=DP) :: smeared_ion_nonself_energy ! the smeared ion representation
    ! agrecocmplx
    logical :: loc_cmplx
    ! jcap: further embedding variables
    real(kind=DP) :: xc_energy_emft,xc_energy_emft_corr
    real(kind=DP), allocatable, dimension(:,:,:,:) :: buffer_fine_emft

    ! jd: ~~~ Energy terms arising from polarisable embedding ~~~
    character(len=27), allocatable :: polemb_energy_term_names(:)
    real(kind=DP), allocatable     :: polemb_energy_term_values(:)
    integer                        :: polemb_n_energy_terms
    real(kind=DP) :: polemb_energy_mm_valence
    real(kind=DP) :: polemb_energy_mm_vdw
    real(kind=DP) :: polemb_energy_mm_es
    real(kind=DP) :: polemb_energy_qmstar_elec_mm_perm
    real(kind=DP) :: polemb_energy_qmstar_elec_mm_ind
    real(kind=DP) :: polemb_energy_qmfull_elec_mm_perm
    real(kind=DP) :: polemb_energy_qmfull_elec_mm_ind
    real(kind=DP) :: polemb_energy_qm_core_mm_perm
    real(kind=DP) :: polemb_energy_qm_core_mm_ind
    real(kind=DP) :: polemb_energy_qmmm_vdw
    real(kind=DP) :: polemb_energy_other
    real(kind=DP) :: dummy_delta_lhxc_energy
    real(kind=DP) :: polemb_energy_sum1
    real(kind=DP) :: polemb_energy_sum2
    type(SPAM3)   :: scratch(6) ! scratch space for P matrices
    integer       :: i_term
    integer       :: imat
    !ep : mermin var
    real(kind=DP), &
         dimension(pub_num_spins,pub_num_kpoints) :: trace_mer
    real(kind=DP), dimension(pub_num_spins) :: smermin
    real(kind=DP) :: entr_cont
    real(kind=DP) :: muext_cont


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    type(SPAM3_ARRAY)              :: pur_denskern_array
    type(SPAM3), allocatable       :: kern_array(:),aug_array(:),ham_array(:)
    integer       :: mrows, ncols, isub, jsub, ll, kk

    ! -------------------------------------------------------------------------

    ! Start Timer
    call timer_clock('hamiltonian_energy_components',1)

    ! --------------------------------------------------------------------------
    ! rc2013: set embedding parameters
    mrows = pur_denskern(1)%mrows
    ncols = pur_denskern(1)%ncols

    allocate(kern_array(pub_num_spins),stat=ierr)
    call utils_alloc_check('hamiltonian_energy_components', 'kern_array',ierr)

    ! ndmh: add spin degeneracy to kernel
    if (pub_num_spins /= 2) then
       call sparse_embed_scale(pur_denskern(1), pub_spin_fac)
    end if

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(1)%iscmplx

    ! Allocate density workspace

    fine_ld1 = mdl%fine_grid%ld1
    fine_ld2 = mdl%fine_grid%ld2
    fine_max_slabs12 = mdl%fine_grid%max_slabs12
    allocate(density_fine(fine_ld1,fine_ld2,fine_max_slabs12,pub_num_spins),&
         stat=ierr)
    call utils_alloc_check('hamiltonian_energy_components','density_fine',ierr)
    allocate(sub_density_fine(fine_ld1,fine_ld2,fine_max_slabs12,pub_num_spins, &
         mrows, ncols), stat=ierr)
    call utils_alloc_check('hamiltonian_energy_components','sub_density_fine',ierr)
    if (pub_aug) then
       allocate(nhat_den_grad(fine_ld1,fine_ld2,fine_max_slabs12,pub_num_spins,&
            0:pub_aug_den_dim),stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components', &
            'nhat_den_grad',ierr)
    end if
    ! JCW: Allocate space for kinetic energy density and gradient of XC energy
    ! JCW: per unit volume wrt tau
    ! (pub_xc_ke_density_required set in xc_init)
    if (pub_xc_ke_density_required) then
       allocate(ke_density_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components',&
            'ke_density_fine',ierr)
       allocate(sub_ke_density_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12, pub_num_spins, mrows, ncols), stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components',&
            'sub_ke_density_fine',ierr)
       allocate(dfdtau_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components',&
            'dfdtau_fine',ierr)
    else
       ! JCW: Allocate a single element dummy array
       allocate(ke_density_fine(1,1,1,1), stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components',&
            'ke_density_fine',ierr)
       allocate(sub_ke_density_fine(1,1,1,1,mrows,ncols), stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components',&
            'sub_ke_density_fine',ierr)
       allocate(dfdtau_fine(1,1,1,1), stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components',&
            'dfdtau_fine',ierr)
    end if

    mrows = rep%overlap%mrows
    ncols = rep%overlap%ncols

    ! ndmh: allocate SPAM3 matrices for projector density kernel and local
    ! ndmh: workspace with dimensions of PAW projectors, and create projector
    ! ndmh: density kernel
    if (pub_aug) then
       allocate(projector_denskern(pub_num_spins),stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components', &
            'projector_denskern',ierr)
       allocate(local_nonlocpot(pub_num_spins),stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components', &
            'local_nonlocpot',ierr)
       allocate(local_dijhat(pub_num_spins),stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components', &
            'local_dijhat',ierr)
       do is=1,pub_num_spins
          projector_denskern(is)%structure='E'
          ! agrecocmplx: projector_denskern is real even when using complex NGWFs
          call sparse_embed_create(projector_denskern(is))
          call sparse_embed_create(local_dijhat(is),projector_denskern(is))
       end do
       ! jcap: need temporary arrays to pass in
       allocate(aug_array(pub_num_spins),stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components',&
            'aug_array',ierr)
       call sparse_embed_extract_from_array(aug_array,projector_denskern)
       call sparse_embed_extract_from_array(kern_array,pur_denskern)

       call aug_projector_denskern(aug_array,kern_array,rep%sp_overlap%p)

       call sparse_embed_destroy_extracted_array(aug_array,projector_denskern,&
            .true.)
       call sparse_embed_destroy_extracted_array(kern_array)
       deallocate(aug_array,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components',&
            'aug_array',ierr)
    end if

    ! ################## CHARGE DENSITY #############################

    ! cks: initialise density_fine
    sub_density_fine = 0.0_DP
    density_fine = 0.0_DP
    integrated_density = 0.0_DP
    integrated_field = 0.0_DP
    do jsub=1,ncols
       do isub=1,mrows
          ! rc2013: extract kernel component to avoid array temporarys
          call sparse_embed_extract_from_array(kern_array,pur_denskern,&
               isub,jsub)
          call density_on_grid(sub_density_fine(:,:,:,:,isub,jsub), &    ! output
               mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
               kern_array, rep%ngwf_overlap%m(isub,jsub), &
               rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
               rep%ngwfs_on_grid(jsub), ngwf_basis(jsub))

          if (pub_aug) then
             nhat_den_grad = 0.0_DP
             ! rc2013: this has not been set up with embedding structures
             call augmentation_density_on_grid(nhat_den_grad,mdl%fine_grid, &
                  mdl%cell,mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp, &
                  mdl%aug_box,kern_array,rep%sp_overlap%m(isub,jsub))
          end if

          ! JCW: if kinetic energy density required for XC, evaluate
          if (pub_xc_ke_density_required) then
             ! Evaluate kinetic energy density on grid
             ! rc2013: EMBED_FIX!
             ! rc2013: for now this is not ready to handle off-diagonal matrices
             if(isub == jsub) &
                  call ke_density_on_grid(sub_ke_density_fine(:,:,:,:,isub,isub), & ! output
                  mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
                  kern_array, rep%ngwf_overlap%m(isub,isub), &
                  rep%ngwfs_on_grid(isub), ngwf_basis(isub))
             ! Also set dfdtau_fine to 0.0_DP
             dfdtau_fine = 0.0_DP
          end if
          call sparse_embed_destroy_extracted_array(kern_array)

          ! cks: find number of electrons from valence charge density
          do is=1,pub_num_spins
             integrated_field(is) = integrated_field(is) + &
                  integrals_trace_on_grid(sub_density_fine(:,:,:,is,isub,jsub),mdl%fine_grid)
             if(pub_debug_on_root) write(stdout,'(a40,f15.9)') 'integrated_field        = ', &
                  integrated_field(is)
          end do
          ! ndmh: find number of electrons from compensation density
          if (pub_aug) then
             do is=1,pub_num_spins
                integrated_nhat(is) = &
                     integrals_trace_on_grid(nhat_den_grad(:,:,:,is,0), &
                     mdl%fine_grid)
             end do
          else
             integrated_nhat(:) = 0.0_DP
          end if
          if (pub_num_spins == 1) then
             integrated_ps_val_density = integrated_field(1)
             integrated_density = integrated_field(1) + integrated_nhat(1)
             integrated_spin = 0.0_DP
          else
             integrated_density = integrated_field(UP) + integrated_field(DN)
             integrated_ps_val_density = integrated_density
             integrated_density = integrated_density + &
                  integrated_nhat(UP) + integrated_nhat(DN)
             integrated_spin = integrated_field(UP) - integrated_field(DN)
             integrated_spin = integrated_spin + &
                  integrated_nhat(UP) - integrated_nhat(DN)
             ! pdh: calculate integrated |spin density|
             integrated_mod_spin = 0.0_DP
             integrated_depleted_spin = 0.0_DP
             do islab12=1,mdl%fine_grid%num_my_slabs12
                do i2=1,mdl%fine_grid%n2
                   do i1=1,mdl%fine_grid%n1
                      integrated_mod_spin = integrated_mod_spin + &
                           abs(sub_density_fine(i1,i2,islab12,UP,isub,jsub) - &
                           sub_density_fine(i1,i2,islab12,DN,isub,jsub))
                      !ddor: Estimate <S^2>
                      if (pub_output_detail >= VERBOSE) then
                         integrated_depleted_spin = &
                              integrated_depleted_spin + &
                              DIM(sub_density_fine(i1,i2,islab12,DN,isub,jsub),&
                              &sub_density_fine(i1,i2,islab12,UP,isub,jsub))
                      endif
                   end do
                end do
             end do
             call comms_reduce('SUM',integrated_mod_spin)
             integrated_mod_spin = integrated_mod_spin * mdl%fine_grid%weight
          end if
          ! rc2013: get total density for XC
          do is=1,pub_num_spins
              density_fine(:,:,:,is) = &
                  density_fine(:,:,:,is) +  sub_density_fine(:,:,:,is,isub,jsub)
              if (pub_xc_ke_density_required) ke_density_fine(:,:,:,is) = &
                   ke_density_fine(:,:,:,is) + sub_ke_density_fine(:,:,:,is,isub,jsub)
          end do
       end do
    end do

    if (pub_output_detail >= VERBOSE) then
       do is=1,pub_num_spins
          if (.not. (pub_eda_scfmi)) then
             call sparse_embed_trace(tr_field(is),pur_denskern(is),rep%overlap)
          else
             ! mjsp: if SCF-MI, then use supermolecule overlap matrix
             call sparse_embed_trace(tr_field(is),pur_denskern(is),rep%overlap_scfmi_full)
          end if
       end do
       if (pub_num_spins == 1) then
          tr_integrated_density = tr_field(1)
          tr_integrated_spin = 0.0_DP
       else
          tr_integrated_density = tr_field(UP) + tr_field(DN)
          tr_integrated_spin = tr_field(UP) - tr_field(DN)
          ! ddor: An estimate of the expectation value of the S^2 operator, used
          !       to correct calculations of splitting between spin-states where
          !       spin contamination has taken place.
          !       Based on Local Density Approximation and derived in
          !       Wang, Becke, Smith. J. Chem. Phys 102 (8) 1995.
          !       The integral of the magnetisation density in the down-spin
          !       excessive areas is subtracted from S(S+1)
          call comms_reduce('SUM',integrated_depleted_spin)
          integrated_depleted_spin = integrated_depleted_spin * &
               mdl%fine_grid%weight
          spin_squared_expectation = (0.25_DP * REAL(pub_spin,kind=DP) * &
               &(REAL(pub_spin,kind=DP) + 2.0_DP) ) + integrated_depleted_spin
       end if
    end if

    ! ############### END CHARGE DENSITY ############################


    ! *********** LOCAL PSEUDOPOTENTIAL ENERGIES **********************
    locpot_energy = 0.0_DP
    nhat_locpot_energy = 0.0_DP
    do is=1,pub_num_spins
       locpot_energy = locpot_energy + &
            integrals_product_on_grid(mdl%fine_grid, mdl%localpseudo_fine, &
            density_fine(:,:,:,is))
       if (pub_aug) then
          nhat_locpot_energy = nhat_locpot_energy + &
               integrals_product_on_grid(mdl%fine_grid,mdl%localpseudo_fine, &
               nhat_den_grad(:,:,:,is,0))
       end if
    end do
    if (pub_aug) locpot_energy = locpot_energy + nhat_locpot_energy

    if (pub_output_detail >= VERBOSE) then
       call sparse_embed_create(localmat, rep%ngwf_overlap)
       do jsub=1,ncols
          do isub=1,mrows
             call integrals_locpot(localmat%m(isub,jsub), &             ! output
                  rep%ngwfs_on_grid(isub), ngwf_basis(isub), &          ! input
                  rep%ngwfs_on_grid(jsub), ngwf_basis(jsub), &          ! input
                  mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
                  mdl%localpseudo_fine)
          end do
       end do
       tr_locpot_energy = 0.0_DP
       do is=1,pub_num_spins
          call sparse_embed_trace(trace, pur_denskern(is), localmat)
          tr_locpot_energy = tr_locpot_energy + trace
       end do
       call sparse_embed_destroy(localmat)
    end if ! verbose

    locps_gzero_energy = mdl%locps_gzero_term * pub_charge

    ! ********END LOCAL PSEUDOPOTENTIAL ENERGIES **********************

    ! cks: allocate simulation cell fine workspace array
    allocate(buffer_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
    call utils_alloc_check('hamiltonian_energy_components','buffer_fine',ierr)
    if (pub_aug .or. pub_pol_emb_pot) then ! jd: polarisable_embedding_calculate
                                           ! relies on lhxc_fine being allocated
       allocate(lhxc_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components','lhxc_fine',ierr)
       do is=1,pub_num_spins
          lhxc_fine(:,:,:,is) = mdl%localpseudo_fine
       end do
    end if

    ! ndmh: Add on augmentation density contribution:
    ! ndmh: \int v_ion \hat{n} dr
    if (pub_aug .and. pub_output_detail >= VERBOSE) then
       ! rc2013: nasty hack to avoid temporary arrays
       do isub=1,mdl%nsub
          call sparse_embed_extract_from_array(kern_array,local_dijhat,&
               isub,isub)
          do is=1,pub_num_spins
             call sparse_scale(kern_array(is),0.0_DP)
          end do
          call augmentation_screen_dij(kern_array, lhxc_fine, &
               mdl%aug_box,mdl%cell,mdl%fine_grid, &
               mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp)
          call sparse_embed_destroy_extracted_array(kern_array,local_dijhat,&
               .true.,isub,isub)
       end do
       tr_nhat_locpot_energy = 0.0_DP
       do is=1,pub_num_spins
          call sparse_embed_trace(trace,projector_denskern(is),local_dijhat(1))
          tr_nhat_locpot_energy = tr_nhat_locpot_energy + trace
       end do
       tr_locpot_energy = tr_locpot_energy + tr_nhat_locpot_energy
    end if

    ! rc2013: subsystem fine grid workspace
    allocate(sub_buffer_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins, mrows, ncols), stat=ierr)
    call utils_alloc_check('hamiltonian_energy_components','sub_buffer_fine',ierr)
    ! jcap: allocate fine grid workspace for EMFT
    if (pub_emft) then
       allocate(buffer_fine_emft(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components','buffer_fine_emft',ierr)
    end if

    ! =================== LOCAL POLARISABLE EMBEDDING ENERGY ===================
    ! jd: Calculate and include the effect of polarisable embedding if necessary
    !     NB: this also *modifies* lhxc_fine to account for QM interactions
    !     with embedding multipoles handled via locpot.
    if(pub_pol_emb_pot) then

       ! jd: Temporary hack for getting a SPAM3_ARRAY representation of
       !     denskern. I assume hamiltonian_energy_components() will be
       !     migrated to SPAM3_ARRAY one day.
       ! jcap: Modified to convert from an array of SPAM3_EMBEDs
       call sparse_array_create(pur_denskern_array, pur_denskern(1)%p, &
            n_spins=pub_num_spins, n_kpoints=PUB_1K)
       do is = 1, pub_num_spins
          call sparse_copy(pur_denskern_array%m(is,PUB_1K),pur_denskern(is)%p)
       end do

       ! jd: polarisable_embedding_calculate() fills polemb and vacuum_polemb
       !     matrices and uses them as scratch space. Here we are not
       !     interested in the matrices, so a dummy is sufficient for scratch.
       scratch(:)%structure='S'
       do imat = 1, 6
          call sparse_create(scratch(imat))
       end do

       dummy_delta_lhxc_energy = 0D0
       call polarisable_embedding_calculate(scratch(1:2), scratch(3:4),&! inout
            scratch(5:6), &                                             ! inout
            polemb_energy_mm_valence, &                                 ! out
            polemb_energy_mm_vdw, &                                     ! out
            polemb_energy_mm_es, &                                      ! out
            polemb_energy_qmstar_elec_mm_perm, &                        ! out
            polemb_energy_qmstar_elec_mm_ind, &                         ! out
            polemb_energy_qmfull_elec_mm_perm, &                        ! out
            polemb_energy_qmfull_elec_mm_ind, &                         ! out
            polemb_energy_qm_core_mm_perm, &                            ! out
            polemb_energy_qm_core_mm_ind, &                             ! out
            polemb_energy_qmmm_vdw, &                                   ! out
            polemb_energy_other, &                                      ! out
            lhxc_fine, dummy_delta_lhxc_energy, &                       ! inout
            sub_density_fine(:,:,:,:,1,1), mdl, rep, ngwf_basis, &      ! in
            pur_denskern_array, &                                       ! in
            polemb_energy_term_names, polemb_energy_term_values, &      ! out
            polemb_n_energy_terms)                                      ! out

       do imat = 1, 6
          call sparse_destroy(scratch(imat))
       end do

       call sparse_array_destroy(pur_denskern_array)

       ! jd: Sanity check on the partitioning of energies
       polemb_energy_sum1 = polemb_energy_mm_valence + polemb_energy_mm_vdw + &
            polemb_energy_mm_es + polemb_energy_qmstar_elec_mm_perm + &
            polemb_energy_qmstar_elec_mm_ind + &
            polemb_energy_qmfull_elec_mm_perm + &
            polemb_energy_qmfull_elec_mm_ind + &
            polemb_energy_qm_core_mm_perm + polemb_energy_qm_core_mm_ind + &
            polemb_energy_qmmm_vdw + polemb_energy_other
       polemb_energy_sum2 = &
            sum(polemb_energy_term_values(1:polemb_n_energy_terms))
       call utils_assert(abs(polemb_energy_sum1-polemb_energy_sum2) < 1D-9, &
            'hamiltonian_energy_components: Energy components do not add up', &
            polemb_energy_sum1, polemb_energy_sum2)
    end if
    ! ================== END LOCAL POLARISABLE EMBEDDING ENERGY ================

    ! %%%%%%%%%%%%% HARTREE ENERGIES %%%%%%%%%%%%%%%%%%%%%%%
    ! cks: initialise hartree_fine
    buffer_fine = 0.0_DP
    sub_buffer_fine = 0.0_DP

    ! ndmh: Add on compensation density before Hartree calculation
    if (pub_aug) then
       density_fine = density_fine + nhat_den_grad(:,:,:,:,0)
       ! rc2013: add this to the subsystem contributions too
       do jsub=1,mdl%nsub
          do isub=1,mdl%nsub
             sub_density_fine(:,:,:,:,isub,jsub) = sub_density_fine(:,:,:,:,isub,jsub) &
                  + nhat_den_grad(:,:,:,:,0)
          end do
       end do
    end if

    ! jd: Calculate either:
    !     (1) the usual Hartree potential
    !     (2) Hartree potential with cutoff Coulomb
    !     (3) Hartree potential with multigrid, which could be
    !         a) electronic Hartree potential, with open BCs,
    !         b) molecular Hartree potential, with smeared-ions, in vacuo,
    !         c) molecular Hartree potential, in implicit solvent.
    !     All of the above are mutually exclusive.
    !     Iff (3c) and the dielectric is not fixed, but instead depends on
    !     the electronic density, an extra term caused by this dependence
    !     (V_eps(r)+V_apolar(r)) is included in the potential to get the correct
    !     derivative. The Hartree energy then needs to be
    !     calculated differently, and hartree_via_multigrid() takes care of
    !     that, as well as calculating the apolar terms to solvation energy.
    hartree_energy = 0.0_DP
    call zero_solvation_terms(solvation_terms, zero_apolar = .true., &
         zero_vacuum = .false.)

    if (.not. pub_turn_off_hartree) then
       ! rc2013: calculate subsystem hartree potentials
       if (pub_coulomb_cutoff) then
          do jsub=1,ncols
             do isub=1,mrows
                call cutoff_coulomb_hartree(sub_buffer_fine(:,:,:,:,isub,jsub), &
                     sub_density_fine(:,:,:,:,isub,jsub), mdl)
                buffer_fine = buffer_fine + sub_buffer_fine(:,:,:,:,isub,jsub)
             end do
          end do
       else if (pub_multigrid_hartree) then
          ! rc2013: EMBED_FIX: this won't work for calculating
          ! subsystem forces
          ! jcap: The multigrid Hartree calculation needs to be done
          ! for the whole system at once, to avoid complaints about
          ! non-integer charges
          call hartree_via_multigrid(buffer_fine, density_fine, &
               mdl%fine_grid, mdl%cell, hartree_energy(1,1), &
               solvation_terms, elements = mdl%elements)
       else
          do jsub=1,ncols
             do isub=1,mrows
                call hartree_on_grid(sub_buffer_fine(:,:,:,:,isub,jsub), &
                     sub_density_fine(:,:,:,:,isub,jsub), mdl%fine_grid, mdl%cell)
                buffer_fine = buffer_fine + sub_buffer_fine(:,:,:,:,isub,jsub)
             end do
          end do
       end if
    else
       buffer_fine=0.0_DP
       sub_buffer_fine=0.0_DP
    endif

    ! cks: scale Hartree potential for hyper Hartree-Fock
    if (pub_hhf_nstates > 0) then
       do is=1,pub_num_spins
          buffer_fine(:,:,:,:) = pub_hhf_factor * buffer_fine(:,:,:,:)
       enddo
       call utils_assert(.not. pub_multigrid_hartree, &
            'Hyper Hartree-Fock (hhf_nstates) is currently incompatible with &
            &multigrid Hartree calculation -- definition of hartree_energy &
            &in hamiltonian_mod would have to be updated')
    endif

    ! jd: If working with the multigrid Hartree, hartree_energy is already
    !     calculated, no need to calculate hartree_energy directly.
    if (.not.pub_multigrid_hartree) then
       do jsub=1,ncols
          do isub=1,mrows
             do is=1,pub_num_spins
                do ll=1,ncols
                   do kk=1,mrows
                      hartree_energy(isub,jsub) = hartree_energy(isub,jsub) + &
                           0.5_DP*integrals_product_on_grid(mdl%fine_grid, &
                           sub_buffer_fine(:,:,:,is,kk,ll), &
                           sub_density_fine(:,:,:,is,isub,jsub))
                   end do
                end do
             end do
          end do
       end do
    end if

    nhat_hartree_energy = 0.0_DP
    if (pub_aug) then

       ! ndmh: add the Hartree potential to the LHXC potential
       lhxc_fine = lhxc_fine + buffer_fine
       do jsub=1,ncols

          ! ndmh: subract off the compensation density to get just the band
          ! ndmh: energy part of the Hartree energy
          do isub=1,mrows
             sub_density_fine(:,:,:,:,isub,jsub) = sub_density_fine(:,:,:,:,isub,jsub) &
                  - nhat_den_grad(:,:,:,:,0)
          end do
       end do
       density_fine = density_fine - nhat_den_grad(:,:,:,:,0)

       ! ndmh: calculate the energy of the compensation density in the Hartree
       ! ndmh: potential of the total soft density
       nhat_hartree_energy = 0.0_DP
       do is=1,pub_num_spins
          nhat_hartree_energy = nhat_hartree_energy + 0.5_DP* &
               integrals_product_on_grid(mdl%fine_grid, buffer_fine(:,:,:,is), &
               density_fine(:,:,:,is))
       end do

    end if

    if (pub_output_detail >= VERBOSE) then
       call sparse_embed_create(localmat, rep%ngwf_overlap)
       tr_hartree_energy = 0.0_DP
       do is=1,pub_num_spins
          do jsub=1,ncols
             do isub=1,mrows
                call integrals_locpot(localmat%m(isub,jsub), &            ! output
                     rep%ngwfs_on_grid(isub), ngwf_basis(isub), &         ! input
                     rep%ngwfs_on_grid(jsub), ngwf_basis(jsub), &         ! input
                     mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, & ! input
                     buffer_fine(:,:,:,is))                              ! input
             end do
          end do
          call sparse_embed_trace(trace, pur_denskern(is), localmat)
          tr_hartree_energy = tr_hartree_energy + 0.5_DP * trace
       end do
       call sparse_embed_destroy(localmat)
       ! ndmh: Calculate energy of augmentation density in Hartree potential
       ! ndmh: (1/2) \int v_H[\tilde{n}+\hat{n}] \hat{n} dr
       if (pub_aug) then
          do isub=1,mrows
             ! rc2013: horrible workaround for temporary array issue!
             call sparse_embed_extract_from_array(kern_array,local_dijhat,&
                  isub,isub)
             do is=1,pub_num_spins
                call sparse_scale(kern_array(is),0.0_DP)
             end do
             call augmentation_screen_dij(kern_array,sub_buffer_fine(:,:,:,:,isub,isub), &
                  mdl%aug_box,mdl%cell,mdl%fine_grid, &
                  mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp)
             call sparse_embed_destroy_extracted_array(kern_array,local_dijhat,&
                  .true.,isub,isub)
          end do
          tr_nhat_hartree_energy = 0.0_DP
          do is=1,pub_num_spins
             call sparse_embed_trace(trace, projector_denskern(is), local_dijhat(is))
             tr_nhat_hartree_energy = tr_nhat_hartree_energy + 0.5_DP * trace
          end do
          tr_hartree_energy = tr_hartree_energy + tr_nhat_hartree_energy
       end if
    end if
    ! %%%%%%%%%% END HARTREE ENERGIES %%%%%%%%%%%%%%%%%%%%%%%

    ! ============== XC MATRIX & ENERGY =================================

    ! ndmh: in PAW, we need to add on the compensation density before finding
    ! ndmh: the XC potential
    if (pub_aug.and.(pub_nhat_in_xc.or.pub_usp)) then
       !sub_density_fine = sub_density_fine + nhat_den_grad(:,:,:,:,0)
       density_fine = density_fine + nhat_den_grad(:,:,:,:,0)
    end if

    ! ndmh: add on core density before calculating xc potential for NLCC
    if (pub_nlcc) then
       factor = 1.0_DP/pub_num_spins
       do is=1,pub_num_spins
          ! jcap: calculate for the full system
          density_fine(:,:,:,is) = density_fine(:,:,:,is) + &
               factor*mdl%core_density_fine
       end do
    end if

    ! cks: initialise
    buffer_fine = 0.0_DP
    sub_buffer_fine = 0.0_DP

    ! ndmh: calculate XC energy and potential
    ! ddor: include optional TDDFT functional switch if necessary
    if (pub_aug) then
       ! JCW: Exit with error if KE density required, since KE density +
       ! augmentation is untested (temporary)
       if(pub_xc_ke_density_required) then
          call utils_abort('Error in hamiltonian_energy_components: &
               &combination of KE density evaluation and augmentation not tested.')
       end if
       call xc_energy_potential(density_fine, xc_energy, buffer_fine, &
            mdl%fine_grid, mdl%cell, pub_aug_den_dim, nhat_den_grad)
    else
       call xc_energy_potential(density_fine, xc_energy, buffer_fine, &
            mdl%fine_grid, mdl%cell, 0, &
            ke_density_fine = ke_density_fine, &
            dfdtau_fine = dfdtau_fine )

       if (pub_emft) then
          ! jcap: If necessary, add on the nlcc charges
          density_fine = sub_density_fine(:,:,:,:,pub_active_region,pub_active_region)
          if (pub_nlcc) then
             do is=1,pub_num_spins
                density_fine(:,:,:,is) = density_fine(:,:,:,is) + &
                     mdl%regions(pub_active_region)%core_density_fine * 0.5_DP * pub_spin_fac
             end do
          end if

          ! rc2013: Now get the EMFT correction
          call xc_emft_calculate(density_fine, xc_energy_emft_corr, &
               buffer_fine_emft, mdl%fine_grid, mdl%cell)

          ! rc2013: add the EMFT correction to the XC energy
          xc_energy = xc_energy + xc_energy_emft_corr
       end if
    end if

    ! ndmh: now subtract off the compensation density and calculate the
    ! ndmh: integral of the pseudo-density \tilde{n} in the XC potential
    if (pub_aug) then

       ! jd: @Possible issue with combining PAW with IS. If both are used,
       !     lhxc_fine is not just a potential, but a potential with a
       !     gradient correction added, i.e. a proper dE/drho, which in these
       !     cases is *not* V_lhxc.
       !     I'm not sure which of these augmentation_screen_dij() expects.

       ! ndmh: zero the dijhat matrices
       do is=1,pub_num_spins
          call sparse_embed_scale(local_dijhat(is),0.0_DP)
       end do

       ! ndmh: calculate screening of dij now if nhat not included in XC
       if (.not.pub_nhat_in_xc) then
          do isub=1,mrows
             ! rc2013: horrible workaround for temporary array issue!
             call sparse_embed_extract_from_array(kern_array,local_dijhat,&
                  isub,isub)
             call augmentation_screen_dij(kern_array,lhxc_fine, &
                  mdl%aug_box,mdl%cell,mdl%fine_grid, &
                  mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp)
             call sparse_embed_destroy_extracted_array(kern_array,local_dijhat,&
                  .true.,isub,isub)
          end do
       end if

       ! Add the XC potential to the LHXC potential
       lhxc_fine = lhxc_fine + buffer_fine

       ! ndmh: calculate screening of dij now if nhat is included in XC
       if (pub_nhat_in_xc) then
          do isub=1,mrows
             ! rc2013: horrible workaround for temporary array issue!
             call sparse_embed_extract_from_array(kern_array,local_dijhat,&
                  isub,isub)
             call augmentation_screen_dij(kern_array,lhxc_fine, &
                  mdl%aug_box,mdl%cell,mdl%fine_grid, &
                  mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp)
             call sparse_embed_destroy_extracted_array(kern_array,local_dijhat,&
                  .true.,isub,isub)
          end do
       end if

       ! ndmh: subract off the compensation density to get just the band
       ! ndmh: energy part of the XC energy
       density_fine = density_fine - nhat_den_grad(:,:,:,:,0)

       ! ndmh: calculate the energy of the PS valence density and of the
       ! ndmh: augmentation density in the XC potential of the total density
       xc_denpot_energy = 0.0_DP
       xc_nhat_denpot_energy = 0.0_DP
       do is=1,pub_num_spins
          xc_denpot_energy = xc_denpot_energy + &
               integrals_product_on_grid(mdl%fine_grid,buffer_fine(:,:,:,is), &
               density_fine)
          xc_nhat_denpot_energy = xc_nhat_denpot_energy  + &
               integrals_product_on_grid(mdl%fine_grid,buffer_fine(:,:,:,is), &
               nhat_den_grad(:,:,:,is,0))
       end do

    end if

    ! ========== END XC MATRIX & ENERGY =================================

    ! ========== PAW NONLOCAL MATRIX ====================================
    paw_sphere_energies(:) = 0.0_DP
    if (pub_aug) then
       ! ndmh: Create workspace to hold nonlocpot
       do is=1,pub_num_spins
          local_nonlocpot(is)%structure = 'H'//rep%postfix
          ! agrecocmplx: needs to be complex for complex NGWFs
          call sparse_embed_create(local_nonlocpot(is),iscmplx=loc_cmplx)
       end do

       ! jcap: need temporary arrays to pass in
       allocate(ham_array(pub_num_spins),stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components',&
            'ham_array',ierr)
       allocate(aug_array(pub_num_spins),stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components',&
            'aug_array',ierr)
       call sparse_embed_extract_from_array(ham_array,local_nonlocpot)
       call sparse_embed_extract_from_array(aug_array,local_dijhat)
       call sparse_embed_extract_from_array(kern_array,projector_denskern)
       if (pub_paw) then
          ! agrecokpt: need to include k-point dependent terms in TB method here?
          call aug_nonlocal_mat(ham_array,aug_array, &
               kern_array,rep%sp_overlap%p,mdl%regions(1)%pseudo_sp, &
               mdl%regions(1)%paw_sp,paw_sphere_energies,show_matrices=.true.)
       else if (pub_usp) then
          call aug_nonlocal_mat(ham_array,aug_array, &
               kern_array,rep%sp_overlap%p,mdl%regions(1)%pseudo_sp, &
               mdl%regions(1)%paw_sp,paw_sphere_energies,show_matrices=.true.)
       end if
       call sparse_embed_destroy_extracted_array(ham_array,local_nonlocpot,.true.)
       call sparse_embed_destroy_extracted_array(aug_array)
       call sparse_embed_destroy_extracted_array(kern_array)
       deallocate(ham_array,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components',&
            'ham_array',ierr)
       deallocate(aug_array,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components',&
            'aug_array',ierr)
    end if
    ! ========== END PAW NONLOCAL MATRIX ================================

    ! cks: deallocate fine grid workspace
    deallocate(sub_buffer_fine,stat=ierr)
    call utils_dealloc_check('hamiltonian_energy_components','sub_buffer_fine',ierr)
    deallocate(buffer_fine,stat=ierr)
    call utils_dealloc_check('hamiltonian_energy_components','buffer_fine',ierr)

    ! jcap: deallocate fine grid space for EMFT
    if(pub_emft) then
       deallocate(buffer_fine_emft,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components','buffer_fine_emft',ierr)
    end if

    if(pub_aug .or. pub_pol_emb_pot) then
       deallocate(lhxc_fine,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components','lhxc_fine',ierr)
    end if

    ! ndmh: deallocate augmentation workspace
    if (pub_aug) then
       do is=pub_num_spins,1,-1
          call sparse_embed_destroy(local_dijhat(is))
          call sparse_embed_destroy(local_nonlocpot(is))
          call sparse_embed_destroy(projector_denskern(is))
       end do
       deallocate(local_dijhat,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components', &
            'local_dijhat',ierr)
       deallocate(local_nonlocpot,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components', &
            'local_nonlocpot',ierr)
       deallocate(projector_denskern,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components', &
            'projector_denskern',ierr)
    end if

    ! JCW: If ke_density_fine evaluated, and in debug mode, compute kinetic
    ! JCW: energy using ke_density_fine, for comparison to energy calculated
    ! JCW: using normal method
    if (pub_xc_ke_density_required.and.pub_debug) then
      kedens_kinetic_energy = &
           ke_density_kinetic_energy(ke_density_fine,mdl%fine_grid)
    end if

    ! cks: Free up density workspace
    ! JCW: Deallocate ke_density_fine and dfdtau_fine. Always allocated since
    ! JCW: we allocate a single element dummy array if
    ! JCW: pub_xc_ke_density_required is false.
    deallocate(dfdtau_fine, stat=ierr)
    call utils_dealloc_check('hamiltonian_energy_components', &
         'dfdtau_fine',ierr)
    deallocate(sub_ke_density_fine, stat=ierr)
    call utils_dealloc_check('hamiltonian_energy_components', &
         'sub_ke_density_fine',ierr)
    ! jd: Seems like this was missed
    deallocate(ke_density_fine, stat=ierr)
    call utils_dealloc_check('hamiltonian_energy_components', &
         'ke_density_fine',ierr)

    if (pub_aug) then
       deallocate(nhat_den_grad,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components', &
            'nhat_den_grad',ierr)
    end if
    deallocate(sub_density_fine,stat=ierr)
    call utils_dealloc_check('hamiltonian_energy_components','sub_density_fine', &
         ierr)
    deallocate(density_fine,stat=ierr)
    call utils_dealloc_check('hamiltonian_energy_components','density_fine', &
         ierr)

    if (pub_scissor_ngroups > 0) then

       allocate(local_scissor_ham(pub_num_spins),stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components', &
            'local_scissor_ham',ierr)

       ! ndmh: Create workspace to hold scissor operator hamiltonian
       do is=1,pub_num_spins
          call sparse_embed_create(local_scissor_ham(is), rep%overlap, &
               iscmplx=loc_cmplx)
       end do

       call hamiltonian_species_scissor(local_scissor_ham, &
            pur_denskern(:),rep,ngwf_basis,mdl)

    end if

    kinetic_energy = 0.0_DP
    nonlocpot_energy = 0.0_DP
    scissor_energy = 0.0_DP
    do is=1,pub_num_spins

       ! cks: calculate kinetic energy
       call sparse_embed_trace(trace, pur_denskern(is), rep%kinet)
       kinetic_energy = kinetic_energy + trace

       ! cks: calculate non-local potential energy
       ! ndmh: only if there are nonlocal projectors
       if (pub_any_nl_proj.and.(.not.(pub_usp))) then
          call sparse_embed_trace(trace, pur_denskern(is), rep%nonlocpot)
          nonlocpot_energy = nonlocpot_energy + trace
       end if
       if (pub_any_nl_proj.and.pub_usp) then
          call sparse_embed_trace(trace, pur_denskern(is), local_nonlocpot(is))
          nonlocpot_energy = nonlocpot_energy + trace
       end if
       if (pub_scissor_ngroups > 0) then
          call sparse_embed_trace(trace, pur_denskern(is), local_scissor_ham(is))
          scissor_energy = scissor_energy + trace
       end if

    end do

    if (pub_scissor_ngroups > 0) then

       ! ndmh: Destroy workspace to hold scissor operator hamiltonian
       do is=pub_num_spins,1,-1
          call sparse_embed_destroy(local_scissor_ham(is))
       end do

       deallocate(local_scissor_ham,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components', &
            'local_scissor_ham',ierr)

    end if

    !qoh: Calculate and include HF exchange if necessary
    hfx_energy = 0.0_DP
    if (pub_use_hfx.or.pub_use_activehfx) then
       do is=1,pub_num_spins
          call sparse_embed_trace(trace,pur_denskern(is),hfexchange(is))
          hfx_energy = hfx_energy - 0.5_DP * trace
       end do
    end if

    ! jd: Include the smeared ion energy corrections, if necessary
    smeared_ion_self_energy = 0.0_DP
    smeared_ion_nonself_energy = 0.0_DP
    if (pub_is_smeared_ion_rep) then
       smeared_ion_self_energy = smeared_ion_E_self
       smeared_ion_nonself_energy = smeared_ion_E_smeared
    end if

    ! =============== DFT+U contribution =======================
    ! ddor: Calculate DFT+U contribution to the total energy
    !       For this we need to remove spin degeneracy factor
    hubbard_energy = 0.0_DP
    if (pub_hubbard_compute_u_or_j) then
       hubbard_e_alpha = 0.0_DP
       hubbard_e_spin = 0.0_DP
    endif
    if (pub_hubbard) then
       ! ddor: Temporarily remove spin degeneracy factor
       if (pub_num_spins /= 2) then
          call sparse_embed_scale(pur_denskern(1), 0.5_DP*pub_num_spins)
       end if
       ! rc2013: extract kernel component to avoid array temporarys
       call sparse_embed_extract_from_array(kern_array,pur_denskern)
       if (.not.pub_hubbard_compute_u_or_j) then
          call hubbard_energy_total(hub, hubbard_energy, kern_array, &
               hub_proj_basis(1), rep%hub_overlap%p, rep%hub_overlap_t%p)
       else
          call hubbard_energy_total(hub, hubbard_energy, kern_array, &
               hub_proj_basis(1), rep%hub_overlap%p, rep%hub_overlap_t%p, &
               hubbard_e_alpha, hubbard_e_spin) !gibo:optional for DFT+alpha/spin
       endif
       call sparse_embed_destroy_extracted_array(kern_array,pur_denskern,.true.)
       if (pub_num_spins /= 2) then
          call sparse_embed_scale(pur_denskern(1), pub_spin_fac)
       end if
    endif
    ! ================= end DFT+U contribution ==================

    !ep : centr_contalculate mermin contributions where needed
    if (pub_mermin) then
       call utils_assert(present(mu), &
            'Mu needs to be passed if pub_mermin is true')
    end if
    entr_cont=0.0_DP
    muext_cont=0.0_DP
    if (pub_num_spins /= 2) then
       call sparse_embed_scale(pur_denskern(1), 0.5_DP*pub_num_spins)
    end if
    if (pub_mermin) then
       do is  = 1, pub_num_spins
          call sparse_embed_trace(trace_mer(is,PUB_1K), pur_denskern(is), &
               rep%overlap)
          call calculate_fermi_entropy_mermin(pur_denskern(is),rep%overlap,1,smermin(is), pub_mermin_cheb)
          entr_cont= entr_cont + smermin(is)*pub_spin_fac* &
               (pub_mermin_smearing_width/k_B)
          muext_cont= muext_cont+ mu(is) * ((trace_mer(is,PUB_1K)*pub_spin_fac) - &
               (rep%n_occ(is,PUB_1K)*pub_spin_fac))
       enddo
    end if
    if (pub_num_spins /= 2) then
       call sparse_embed_scale(pur_denskern(1), pub_spin_fac)
    end if

    ! rc2013: deallocate extracted kernel
    deallocate(kern_array,stat=ierr)
    call utils_dealloc_check('hamiltonian_energy_components', 'kern_array',ierr)

    ! ndmh: add up total energy
    ! rc2013: add up all the subsystem contributions
    total_hartree_energy = 0.0_DP
    do isub=1,mrows
       do jsub=1,ncols
          total_hartree_energy = total_hartree_energy + hartree_energy(isub,jsub)
       end do
    end do
    total_energy = kinetic_energy + locpot_energy + nonlocpot_energy &
            + total_hartree_energy + xc_energy + hfx_energy + hubbard_energy &
            + mdl%ewald_energy + pub_dispersion_energy &
            + smeared_ion_self_energy + smeared_ion_nonself_energy &
            + solvation_terms%E_apolar_cavitation &
            + solvation_terms%E_apolar_disrep &
            + implicit_solvent_boltzmann_energy_contrib(solvation_terms) &
            + scissor_energy + entr_cont + muext_cont
    if (pub_paw) then
       total_sphere_energy = paw_sphere_energies(paw_en_dij0) &
            + paw_sphere_energies(paw_en_ehart) &
            + paw_sphere_energies(paw_en_exc) &
            - paw_sphere_energies(paw_en_etxc) &
            - paw_sphere_energies(paw_en_exc_core)
       total_energy = total_energy + total_sphere_energy
    end if

    if(pub_pol_emb_pot) then
       ! jd: Add the polarisable embedding potential energy terms, if any.
       !     The ones tagged as ignored have been zeroed already.
       do i_term = 1, polemb_n_energy_terms
          total_energy = total_energy + polemb_energy_term_values(i_term)
       end do
    end if

    if (pub_on_root) then
       write(stdout,'(/a)') '===================================================&
            &============================='
       write(stdout,'(11x,a)') '---------------- ENERGY COMPONENTS (Eh) &
            &----------------'
       write(stdout,'(11x,a,f24.14,a)') &
            '| Kinetic                    :',            kinetic_energy,' |'
       if (pub_xc_ke_density_required.and.pub_debug) then
          ! JCW: If pub_debug and KE-density-dependent functional used, then
          ! JCW: output kinetic energy evaluated by integrating ke_density_fine
          ! JCW: for comparison to kinetic energy evaluated in standard way
          write(stdout,'(11x,a,f24.14,a)') &
               '| Kinetic (using KE density) :',  kedens_kinetic_energy,' |'
       end if
       if (.not.pub_paw) then

          if (pub_is_smeared_ion_rep) then
             if (.not.any(pub_pspot_bc_is_periodic(1:3))) then
                ! Fully open BCs
                write(stdout,'(a,f24.14,a)') &
                     '           | Pseudo. (local,OBC,corr''d) :', &
                     locpot_energy,' |'
             else if (all(pub_pspot_bc_is_periodic(1:3))) then
                ! Fully periodic BCs
                write(stdout,'(a,f24.14,a)') &
                     '           | Pseudo. (local,PBC,corr''d) :', &
                     locpot_energy,' |'
                ! TODO Determine whether the locps_gzero_energy is meaningful
                !      when the local pseudopotential is corrected for smeared
                !      ionic cores.
                !write(stdout,'(a,f24.14,a)') &
                !     '           | Pseudo (non-coul chg cor)  :', &
                !locps_gzero_energy
             else
                ! Mixed open/periodic BCs
                write(stdout,'(a,f24.14,a)') &
                     '           | Pseudo. (loc.,mixed BC,corr):', &
                     locpot_energy,' |'
             end if
          else
             if (.not.any(pub_pspot_bc_is_periodic(1:3))) then
                ! Fully open BCs
                write(stdout,'(a,f24.14,a)') &
                     '           | Pseudopot. (local, OBC):    ', &
                     locpot_energy,' |'
             else if (all(pub_pspot_bc_is_periodic(1:3))) then
                ! Fully periodic BCs
                write(stdout,'(a,f24.14,a)') &
                     '           | Pseudopotential (local)    :', &
                     locpot_energy,' |'
                write(stdout,'(a,f24.14,a)') &
                     '           | Pseudo (non-coul chg cor)  :', &
                locps_gzero_energy
             else
                ! Mixed open/periodic BCs
                write(stdout,'(a,f24.14,a)') &
                     '           | Pseudo. (local,mixed BC):   ', &
                     locpot_energy,' |'
             end if
          end if

          write(stdout,'(11x,a,f24.14,a)') &
               '| Pseudopotential (non-local):',       nonlocpot_energy,' |'
          if (pub_is_smeared_ion_rep) then
             write(stdout,'(a,f24.14,a)') &
                  '           | Hartree (molecular)        :', &
                  total_hartree_energy,' |'
          else
             if (pub_multigrid_hartree) then
                if (.not.any(pub_multigrid_bc_is_periodic)) then  ! Fully open BCs
                   write(stdout,'(a,f24.14,a)') &
                        '           | Hartree (OBC)              :', &
                        hartree_energy(1,1),' |'
                else if (all(pub_multigrid_bc_is_periodic)) then  ! Fully periodic BCs
                   write(stdout,'(a,f24.14,a)') &
                        '           | Hartree (PBC)              :', &
                        hartree_energy(1,1),' |'
                else                                ! Mixed BCs
                   write(stdout,'(a,f24.14,a)') &
                        '           | Hartree (Mixed BCs)        :', &
                        hartree_energy(1,1),' |'
                end if
             else
                if(pub_output_detail >= VERBOSE .and. mdl%nsub .gt. 1) then
                   do isub=1,mrows
                      do jsub=1,ncols
                         write(stdout,'(a,i1,a2,i1,a,f24.14,a)') &
                              '           |           Hartree (',isub,', ', &
                              jsub, ')   :', hartree_energy(isub,jsub),' |'
                      end do
                   end do
                endif
                write(stdout,'(a,f24.14,a)') &
                     '           | Hartree                    :', &
                     total_hartree_energy,' |'
             end if
          end if
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exchange-correlation       :',              xc_energy,' |'
          ! jcap: Print out EMFT correction
          if (pub_emft) then
             write(stdout,'(11x,a,f24.14,a)') &
                  '| Exchange-correlation (EMFT):',              xc_energy_emft_corr,' |'
          end if
       else
          write(stdout,'(11x,a,f24.14,a)') &
               '| Core Hartree               :',          locpot_energy,' |'

          write(stdout,'(11x,a,f24.14,a)') &
               '| CoreHart (non-coul chg cor):', &
          locps_gzero_energy
          write(stdout,'(11x,a,f24.14,a)') &
               '| Hartree                    :',         total_hartree_energy,' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exchange-correlation       :',              xc_energy,' |'
       end if
       if (pub_pol_emb_pot) then
          ! jd: Print the polarisable embedding energy terms
          !     Ignore the ones marked with a hash
          do i_term = 1, polemb_n_energy_terms
             if(index(polemb_energy_term_names(i_term),'#') == 0) then
                write(stdout,'(11x,a,a,a,f24.14,a)') '| ', &
                     polemb_energy_term_names(i_term), ':', &
                     polemb_energy_term_values(i_term),' |'
             end if
          end do
       end if
       if (pub_use_hfx.or.pub_use_activehfx) write(stdout,'(11x,a,f24.14,a)') &
            '| Hartree-Fock exchange      :',       hfx_energy,' |'
       if ((pub_hubbard).and.(.not.pub_cdft)) then
            write(stdout,'(11x,a,f24.14,a)') &
            '| Hubbard DFT+U correction   :',   hubbard_energy,' |'
       elseif (pub_cdft) then
            write(stdout,'(11x,a,f24.14,a)') &
            '| cDFT [+(DFT+U)] correction :',   hubbard_energy,' |'
       endif
       !gibo: for U-(J-)-calculate runs, print also DFT+alpha/spin energy
       if (pub_hubbard_compute_u_or_j) then
            write(stdout,'(11x,a,f24.14,a)') &
            '| \Sum alpha*Tr[n]           :',  hubbard_e_alpha,' |'
            write(stdout,'(11x,a,f24.14,a)') &
            '| \Sum s*(Tr[n_UP]-Tr[n_DN]) :',   hubbard_e_spin,' |'
            write(stdout,'(11x,a,f24.14,a)') &
            '| Total for computing U or J :',&
              total_energy - hubbard_e_spin - hubbard_e_alpha,' |'
       endif
       if (.not.any(pub_ion_ion_bc_is_periodic(1:3)) ) then
          ! JCW: Fully open BCs
          write(stdout,'(11x,a,f24.14,a)') &
               '| Ion-ion (open BC)          :',  mdl%ewald_energy,' |'
       else if (all(pub_ion_ion_bc_is_periodic(1:3)) ) then
          ! JCW: Fully periodic BCs
          write(stdout,'(11x,a,f24.14,a)') &
               '| Ewald                      :',  mdl%ewald_energy,' |'
       else
          ! JCW: Mixed open/periodic BCs
          write(stdout,'(11x,a,f24.14,a)') &
               '| Ion-ion (mixed BC)         :',  mdl%ewald_energy,' |'
       end if
       if (pub_paw) then
          write(stdout,'(11x,a,f24.14,a)') &
               '| Dij0                       :',  paw_sphere_energies(2),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Dijhartree                 :',  paw_sphere_energies(3),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exc [AE core]              :',  paw_sphere_energies(4),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exc DC [AE core]           :',  paw_sphere_energies(5),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exc [PS core]              :',  paw_sphere_energies(6),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exc DC [PS core]           :',  paw_sphere_energies(7),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exc [Core]                 :',  paw_sphere_energies(9),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Total [Sphere]             :',     total_sphere_energy,' |'
       end if
       !ep: print entropy contribution if needed
       if (pub_mermin) then
          write(stdout,'(11x,a,f24.14,a)') &
               '| Entropy Contribution (-T*S):',     entr_cont,' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Mu Contribution            :',     muext_cont,' |'
       end if
       !ep: print entropy contribution if needed
       if (pub_dispersion /= '0') write(stdout,'(11x,a,f24.14,a)') &
            '| Dispersion Correction      :', pub_dispersion_energy,&
            ' |'
       if (pub_is_smeared_ion_rep) then
          write(stdout,'(a,f24.14,a)') &
               '           | Smeared ion non-self corr. :', &
               smeared_ion_nonself_energy,' |'
          write(stdout,'(a,f24.14,a)') &
               '           | Smeared ion self corr.     :', &
               smeared_ion_self_energy,' |'
       end if
       if (pub_is_include_apolar) then
          write(stdout,'(a,f24.14,a)') &
               '           | Solvent cavitation energy  :', &
               solvation_terms%E_apolar_cavitation,' |'
          write(stdout,'(a,f24.14,a)') &
               '           | Solute-solvent dis-rep en. :', &
               solvation_terms%E_apolar_disrep,' |'
       end if

       if (pub_is_pbe /= "NONE") then
          write(stdout,'(a,f24.14,a)') &
               '           | Elect. mobile ion energy   :', &
               solvation_terms%E_elec_mob,' |'
          write(stdout,'(a,f24.14,a)') &
               '           | Osmotic pressure energy    :', &
               solvation_terms%E_osmo,' |'
          write(stdout,'(a,f24.14,a)') &
               '           | Acc. corr. (steric) energy :', &
               solvation_terms%E_acc,' |'
          write(stdout,'(a,f24.14,a)') &
               '           | Ionic atmo. rearrang. en.  :', &
               solvation_terms%E_atmo,' |'
          write(stdout,'(a,f24.14,a)') &
               '           | Chemical pot. contribution :', &
               solvation_terms%E_chempot,' |'
       end if

       write(stdout,'(11x,a,f24.14,a)') '| Total                      :',&
            total_energy,' |'

       ! jd: If edft container passed and we are actually doing EDFT,
       !     also print the entropic contribution and total free energy
       if(present(edft) .and. pub_edft) then
          write(stdout,'(11x,a)') '|---------------------------------------&
            &---------------|'
          write(stdout,'(11x,a,f24.14,a)') '| Entropic contribution      :',&
               edft%free_energy - total_energy, ' |'
          write(stdout,'(11x,a,f24.14,a)') '| Total free energy          :',&
               edft%free_energy,' |'
       end if

       write(stdout,'(11x,a)') '----------------------------------------&
            &----------------'

       if (pub_output_detail >= VERBOSE) then
          write(stdout,'(11x,a)') '------ LOCAL ENERGY COMPONENTS FROM &
               &MATRIX TRACES ------'
          if (.not.pub_paw) then
             write(stdout,'(11x,a,f24.14,a)') &
                  '| Pseudopotential (local)    :', tr_locpot_energy,' |'
             write(stdout,'(11x,a,f24.14,a)') &
                  '| Hartree                    :',tr_hartree_energy,' |'
          else
             write(stdout,'(11x,a,f24.14,a)') &
                  '| Core Hartree               :', tr_locpot_energy,' |'
             write(stdout,'(11x,a,f24.14,a)') &
                  '| Hartree (Valence)          :',tr_hartree_energy,' |'
          end if
          write(stdout,'(11x,a)') '-------------------------------------&
               &-------------------'
       end if

       if (pub_aug) then
          write(stdout,'(11x,a,f24.14)')  &
               'Integrated compensation chrg.:', sum(integrated_nhat( &
               1:pub_num_spins))
       end if

       write(stdout,'(11x,a,f24.14)')  &
            'Integrated density           :', integrated_density
       if (pub_num_spins > 1) then
          write(stdout,'(11x,a,f24.14)') &
               'Integrated spin density      :', integrated_spin
          write(stdout,'(11x,a,f24.14)') &
               'Integrated |spin density|    :', integrated_mod_spin
       end if

       if (pub_output_detail >= VERBOSE) then
          if (pub_num_spins > 1) write(stdout,'(11x,a,f24.14)') &
               'Local density approx. <S^2>  :', &
               spin_squared_expectation
          write(stdout,'(11x,a,f24.14)')  &
               'Integrated density tr(KS)    :', tr_integrated_density
          if (pub_num_spins > 1) write(stdout,'(11x,a,f24.14)') &
               'Integrated spin tr(KS)       :',tr_integrated_spin
       end if

       write(stdout,'(a/)') '===================================================&
            &============================='

    end if

    if (pub_pol_emb_pot) then
       ! jd: Deallocate the energy terms allocated in
       !     polarisable_embedding_calculate() and calls therein
       deallocate(polemb_energy_term_values, stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components', &
            'polemb_energy_term_values', ierr)
       deallocate(polemb_energy_term_names, stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components', &
            'polemb_energy_term_names', ierr)
    end if

    ! ndmh: remove spin-degeneracy from density kernel
    if (pub_num_spins /= 2) then
       call sparse_embed_scale(pur_denskern(1), 0.5_DP*pub_num_spins)
    end if

    ! Stop Timer
    call timer_clock('hamiltonian_energy_components',2)

  end subroutine hamiltonian_energy_components

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hamiltonian_update_kinet_kp(kinet, overlap, ngwfs_on_grid, &
             ngwf_basis, mdl, kpt)

    !=========================================================================!
    ! This subroutine adds the extra terms of the kinetic energy arising in   !
    ! the KP method for non-Gamma k-points.                                   !
    ! Currently works for a single k-point.                                   !
    !                                                                         !
    ! Arguments:                                                              !
    !     kinet : on-entry, original kinetic matrix (non k-dependent)         !
    !     overlap : overlap matrix <bra|ket>                                  !
    !     ngwfs_on_grid : NGWFs of basis ngwfs_basis, on grid                 !
    !     kpt : KPOINT type for the current k-point                           !
    !                                                                         !
    ! Written by Andrea Greco on 08/02/2016.                                  !
    !=========================================================================!

    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_axpy
    use integrals, only: integrals_grad
    use datatypes, only: FUNCTIONS
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use geometry, only: POINT

    implicit none

    ! Argument
    type(SPAM3), intent(inout)    :: kinet           ! kinetic matrix
    type(SPAM3), intent(in)       :: overlap         ! overlap matrix
    type(FUNCTIONS), intent(in)   :: ngwfs_on_grid   ! NGWFs
    type(FUNC_BASIS), intent(in)  :: ngwf_basis      ! NGWF basis type
    type(POINT), intent(in)       :: kpt             ! k-point
    type(MODEL), intent(in)       :: mdl

    ! Local variables
    type(SPAM3)                   :: grad(3)         ! Grad matrix to compute KP term
    real(kind=DP)                 :: kcart(3)        ! k-point in Cartesians
    integer                       :: dim_cart        ! dimension counter

    ! Cartesian coordinates of k-point
    kcart(1) = kpt%x
    kcart(2) = kpt%y
    kcart(3) = kpt%z

    ! Calculate grad matrix elements
    do dim_cart=1,3
       call sparse_create(grad(dim_cart), overlap)
    end do

    call integrals_grad(grad, ngwfs_on_grid, ngwf_basis, &
         ngwfs_on_grid, ngwf_basis, mdl%cell, mdl%fftbox)

    ! KE(Gamma) - i k*<bra| grad | ket>
    do dim_cart=1,3
       call sparse_axpy(kinet ,grad(dim_cart), &
            cmplx(0.0_DP,-kcart(dim_cart),kind=DP))
    end do

    ! KE(Gamma) - i k*<bra| grad | ket> + 0.5 * |k|^2 <bra|ket>
    call sparse_axpy(kinet, overlap, &
         0.5_DP*(kcart(1)*kcart(1) + kcart(2)*kcart(2) + &
         kcart(3)*kcart(3)))

    ! deallocate temp grad matrix
    do dim_cart=3,1,-1
       call sparse_destroy(grad(dim_cart))
    end do

  end subroutine hamiltonian_update_kinet_kp

  subroutine hamiltonian_apply_phases_tb_method(in_matrix, ngwf_basis, mdl, &
    kpt, out_matrix)
    !=========================================================================!
    ! This subroutine includes the k-point dependent phases in a given matrix,!
    ! effectively overlap/hamiltonian matrices, in the TB method.             !
    ! Currently works for a single k-point.                                   !
    !                                                                         !
    ! Arguments:                                                              !
    !     in_matrix: on entry, the original matrix at Gamma point             !
    !     kpt : KPOINT type for the current k-point                           !
    !     out_matrix: if present, the matrix with k-dependent terms           !
    !                                                                         !
    ! Written by Andrea Greco on 05/09/2016, adapted from similar routine in  !
    ! bandstructure module.                                                   !
    !=========================================================================!
    ! rc2013: PAR WARNING: this subroutine only works with 1 subsystem

    use comms, only: pub_my_proc_id
    use constants, only: TWO_PI
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(.DOT.), operator(-), &
        operator(*)
    use model_type, only: MODEL
    use sparse, only: SPAM3, sparse_index_length, sparse_copy, &
        sparse_put_block, sparse_get_block, sparse_generate_index
    use utils, only: utils_alloc_check, utils_dealloc_check, &
        utils_assert

    implicit none

    ! Argument
    type(SPAM3), intent(inout)           :: in_matrix      ! input matrix
    type(FUNC_BASIS), intent(in)         :: ngwf_basis     ! NGWF basis type
    type(MODEL), intent(in)              :: mdl            ! model
    type(POINT), intent(in)              :: kpt            ! k-point
    type(SPAM3), intent(inout), optional :: out_matrix     ! output matrix

    ! Local variables
    integer                       :: ierr
    integer                       :: idxlen             ! index length
    integer, allocatable          :: idx(:)             ! sparse index
    integer                       :: max_ngwfs_on_atom  ! max NGWFs on any atom
    integer                       :: loc_iat            ! atom counter on local proc
    integer                       :: iat, jat           ! atom counters
    integer                       :: orig_iat, orig_jat ! atom counters input file order
    integer                       :: jdx                ! sparse matrix counter
    integer                       :: dim                ! dimension counter
    real(kind=DP)                 :: kfrac(3)           ! k-point in frac coordinates
    real(kind=DP)                 :: dispfrac(3)        ! Displacement vector (fractional)
    type(POINT)                   :: atdisp             ! atom displacement
    complex(kind=DP)              :: i2pi               ! 2i * PI
    complex(kind=DP)              :: zphase             ! Complex translation phase
    complex(kind=DP), allocatable :: block(:,:)         ! Sparse matrix block
    real(kind=DP), parameter      :: recip_twopi = 1.0_DP / TWO_PI

    ! only for complex matrices
    call utils_assert(in_matrix%iscmplx, 'Error in &
         &hamiltonian_apply_phases_tb_method: in_matrix must be complex.')

    i2pi = cmplx(0.0_DP,-TWO_PI,kind=DP)

    ! generate index of sparse matrix
    idxlen = sparse_index_length(in_matrix)
    allocate(idx(idxlen),stat=ierr)
    call utils_alloc_check('hamiltonian_apply_phases_tb_method', &
         'idx',ierr)
    call sparse_generate_index(idx, in_matrix)

    ! block size
    max_ngwfs_on_atom = maxval(ngwf_basis%num_on_atom)
    allocate(block(max_ngwfs_on_atom,max_ngwfs_on_atom),stat=ierr)
    call utils_alloc_check('hamiltonian_apply_phases_tb_method', &
         'block',ierr)
    ! initialise block
    block = (0.0_DP,0.0_DP)
    ! k-point in fractional coordinates
    kfrac(1) = kpt .DOT. mdl%cell%a1
    kfrac(2) = kpt .DOT. mdl%cell%a2
    kfrac(3) = kpt .DOT. mdl%cell%a3
    kfrac = kfrac * recip_twopi

    ! copy input to output matrix if
    ! output matrix is specified
    if (present(out_matrix)) then
       call utils_assert(out_matrix%iscmplx, 'Error in &
            &hamiltonian_apply_phases_tb_method: out_matrix must be complex.')

       call sparse_copy(out_matrix, in_matrix)
    end if

    ! loop over first atom
    do loc_iat=1,mdl%par%num_atoms_on_proc(pub_my_proc_id)
       iat = loc_iat + mdl%par%first_atom_on_proc(pub_my_proc_id) - 1
       orig_iat = mdl%par%orig_atom(iat)
       ! loop over second atom
       do jdx=idx(loc_iat),idx(loc_iat+1)-1
          jat = idx(jdx)
          orig_jat = mdl%par%orig_atom(jat)

          atdisp = mdl%elements(orig_iat)%centre - &
             mdl%elements(orig_jat)%centre

          dispfrac(1) = atdisp .DOT. mdl%cell%b1
          dispfrac(2) = atdisp .DOT. mdl%cell%b2
          dispfrac(3) = atdisp .DOT. mdl%cell%b3

          dispfrac = dispfrac * recip_twopi
          ! initialize extra phase
          zphase = (1.0_DP,0.0_DP)
          ! loop over directions
          do dim=1,3
             if (dispfrac(dim) > 0.5_DP) &
                zphase = zphase * exp(i2pi * kfrac(dim))
             if (dispfrac(dim) < -0.5_DP) &
                zphase = zphase * exp(-i2pi * kfrac(dim))
          end do

          if (zphase /= (1.0_DP,0.0_DP)) then
             call sparse_get_block(block,in_matrix,jat,iat)
             block = block * zphase
             ! modified matrix is a different matrix than input
             if (present(out_matrix)) then
                call sparse_put_block(block,out_matrix,jat,iat)
             ! input matrix is modified directly
             else
                call sparse_put_block(block,in_matrix,jat,iat)
             end if
          end if
       ! end loop over second atom
       end do
    ! end loop over first atom
    end do

    deallocate(block,stat=ierr)
    call utils_dealloc_check('hamiltonian_apply_phases_tb_method', &
         'block',ierr)

    deallocate(idx,stat=ierr)
    call utils_dealloc_check('hamiltonian_apply_phases_tb_method', &
         'idx',ierr)

  end subroutine hamiltonian_apply_phases_tb_method


  subroutine hamiltonian_soc_matrices(hamso,ham,rhoij,rep,mdl)

    !=========================================================================!
    ! This subroutine calculates the spin-orbit terms in the Hamiltonian for  !
    ! post-processing calculation (bandstructure projection, DOS) purposes    !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use paw, only: paw_dij_so
    use rundat, only: pub_num_spins
    use sparse, only: SPAM3, sparse_copy, sparse_create, sparse_destroy
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_transpose, sparse_embed_transpose_structure, &
         sparse_embed_product, sparse_embed_axpy
    use utils, only: utils_alloc_check, utils_dealloc_check, &
        utils_assert

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: hamso(4)          ! Hamiltonian matrices
    type(SPAM3_EMBED), intent(in) :: rhoij(pub_num_spins) ! Projector denskern
    type(NGWF_HAM), intent(in) :: ham
    type(NGWF_REP), intent(in) :: rep               ! NGWF Representation
    type(MODEL), intent(in)    :: mdl

    ! Local variables
    integer                       :: ierr
    integer                       :: isp, isub
    type(SPAM3_EMBED), allocatable :: dijso(:)  ! Spin-orbit terms
    type(SPAM3_EMBED) :: sp_overlap
    type(SPAM3_EMBED) :: ps_overlap, spdijso
    type(SPAM3), allocatable :: reg_dijso(:)  ! Spin-orbit terms
    type(SPAM3), allocatable :: reg_rhoij(:)  ! Spin-orbit terms

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering hamiltonian_soc_matrices'

    ! Prepare perturbative SOC
    allocate(dijso(4), stat=ierr)
    call utils_alloc_check('hamiltonian_soc_matrices','dijso',ierr)
    allocate(reg_dijso(4), stat=ierr)
    call utils_alloc_check('hamiltonian_soc_matrices','reg_dijso',ierr)
    allocate(reg_rhoij(pub_num_spins), stat=ierr)
    call utils_alloc_check('hamiltonian_soc_matrices','reg_dijso',ierr)
    do isp=1,4
       dijso(isp)%structure = 'E'
       call sparse_embed_create(dijso(isp), iscmplx=.true.)
    end do
    call sparse_embed_create(sp_overlap, rep%sp_overlap, iscmplx = .true.)
    call sparse_embed_copy(sp_overlap, rep%sp_overlap) ! copy real to complex
    call sparse_embed_transpose_structure(ps_overlap%structure, sp_overlap)
    call sparse_embed_create(ps_overlap, iscmplx = .true.)
    call sparse_embed_transpose(ps_overlap, sp_overlap)
    call sparse_embed_create(spdijso, sp_overlap, dijso(1), iscmplx = .true.)

    ! Create Hamiltonian matrix for this k-point and spin
    ! rc2013: loop over all regions
    ! This is a hideous structure but the alternative is going directly into PAW
    do isub=1,rep%nsub
       do isp=1,4
          call sparse_create(reg_dijso(isp), dijso(isp)%m(isub,isub))
          call sparse_copy(reg_dijso(isp), dijso(isp)%m(isub,isub))
       end do
       do isp=1,pub_num_spins
          call sparse_create(reg_rhoij(isp), rhoij(isp)%m(isub,isub))
          call sparse_copy(reg_rhoij(isp), rhoij(isp)%m(isub,isub))
       end do
       call paw_dij_so(reg_dijso, 4, reg_rhoij, mdl%regions(isub)%paw_sp)
       do isp=1,4
          call sparse_copy(dijso(isp)%m(isub,isub), reg_dijso(isp))
          call sparse_destroy(reg_dijso(isp))
       end do
       do isp=1,pub_num_spins
          call sparse_destroy(reg_rhoij(isp))
       end do
    end do
    do isp=1,4
       call sparse_embed_product(spdijso, sp_overlap, dijso(isp))
       call sparse_embed_product(hamso(isp), spdijso, ps_overlap)
    end do
    do isp=1,pub_num_spins
       call sparse_embed_axpy(hamso(isp), ham%ham(isp), 1.0_DP)
    end do
    if (pub_num_spins == 1) then
       call sparse_embed_axpy(hamso(2), ham%ham(1), 1.0_DP)
    end if

    ! Cleanup perturbative SOC
    call sparse_embed_destroy(spdijso)
    call sparse_embed_destroy(ps_overlap)
    call sparse_embed_destroy(sp_overlap)
    do isp=4,1,-1
       call sparse_embed_destroy(dijso(isp))
    end do
    deallocate(dijso, stat=ierr)
    call utils_dealloc_check('hamiltonian_soc_matrices','dijso',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving hamiltonian_soc_matrices'

  end subroutine hamiltonian_soc_matrices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hamiltonian_proj_embed_matrix(ham, unproj_ham, overlap, denskern)

    !=========================================================================!
    ! This subroutine calculates the projected embedding Hamiltonian matrix:  !
    ! H_{AA}^proj = H_{AA} - S_{AB}K_{BB}H_{BA} - H_{AB}K_{BB}S_{BA}          !
    !                      + S_{AB}K_{BB}H_{BB}K_{BB}S_{BA}                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! unproj_ham (inout) :: Projected Hamiltonian matrix to be calculated.    !
    !      ham (inout)   :: Unprojected  Hamiltonian matrix.                  !
    !  overlap (in)      :: Overlap matrix for all systems.                   !
    ! denskern (in)      :: Density kernel for all systems.                   !
    !-------------------------------------------------------------------------!
    ! Written by Robert Charlton, 14/12/2017.                                 !
    ! Based on hamiltonain_proj_cond_matrix by Laura Ratcliff.                !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout, VERBOSE, max_spins
    use rundat, only: pub_output_detail, cond_fixed_shift, cond_calc_max_eigen, &
         cond_shift_buffer, pub_debug, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_axpy, sparse_embed_copy, &
        sparse_embed_create, sparse_embed_destroy, sparse_embed_product, &
        sparse_embed_trace, sparse_embed_diagnose
    use sparse, only: sparse_trace, sparse_copy, sparse_axpy

    implicit none

    ! Arguments -- for now don't pass anything back out; this is just for testing
    type(SPAM3_EMBED), intent(inout)  :: unproj_ham
    type(SPAM3_EMBED), intent(inout)  :: ham
    type(SPAM3_EMBED), intent(in)     :: overlap
    type(SPAM3_EMBED), intent(in)     :: denskern

    ! Local variables
    integer :: is    ! spin loop counter
    type(SPAM3_EMBED) :: sab_kbb, kbb_sba, hab_pba, pab_hba, pab_hbb, &
                         pab_hbb_pba  ! matrix products
    real(kind=DP)     :: debug_trace  ! matrix trace for debugging purposes
    integer           :: isub, jsub   ! region iterators

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
         &Entering hamiltonian_proj_embed_matrix'

    ! agrecocmplx
    !loc_cmplx = dkn(1)%iscmplx

    ! rc2013: copy original Hamiltonian into unprojected form
    ! This should really be done separately as otherwise we make things worse in the SCF cycle
    call sparse_embed_copy(unproj_ham, ham)

    ! rc2013: loop over all regions
    do isub=1,ham%mrows
        do jsub=1,ham%ncols
            ! rc2013: only build projectors from other regions
            if(jsub .ne. isub) then
                ! Allocate storage for matrices required for projection
                ! rc2013: specify the row/col we want to extract
                call sparse_embed_create(sab_kbb, overlap, denskern, arow=isub, bcol=jsub)
                call sparse_embed_create(kbb_sba, denskern, overlap, arow=jsub, bcol=isub)
                ! rc2013: need to specify the column structure index
                call sparse_embed_create(hab_pba, unproj_ham, kbb_sba, arow=isub, acol=jsub, cstruc=isub)
                call sparse_embed_create(pab_hba, sab_kbb, unproj_ham, brow=jsub, bcol=isub, rstruc=isub)
                call sparse_embed_create(pab_hbb, sab_kbb, unproj_ham, brow=jsub, bcol=jsub, rstruc=isub)
                ! rc2013: the structures should already be set here
                call sparse_embed_create(pab_hbb_pba, pab_hbb, kbb_sba, rstruc=isub, cstruc=isub)

                ! rc2013: now construct projected Hamiltonian
                call sparse_embed_product(sab_kbb, overlap%m(isub,jsub), denskern%m(jsub,jsub))
                call sparse_embed_product(kbb_sba, denskern%m(jsub,jsub), overlap%m(jsub,isub))
                call sparse_embed_product(hab_pba, unproj_ham%m(isub,jsub), kbb_sba)
                call sparse_embed_product(pab_hba, sab_kbb, unproj_ham%m(jsub,isub))
                call sparse_embed_product(pab_hbb, sab_kbb, unproj_ham%m(jsub,jsub))
                call sparse_embed_product(pab_hbb_pba, pab_hbb, kbb_sba)

                ! rc2013: now add up all the matrices into projected Hamiltonian
                call sparse_axpy(ham%m(isub,isub), pab_hba%p, -1.0_DP)
                call sparse_axpy(ham%m(isub,isub), hab_pba%p, -1.0_DP)
                call sparse_axpy(ham%m(isub,isub), pab_hbb_pba%p, 1.0_DP)

                ! lr408: Calculate Tr[T+KSKTM] as a measure of orthogonality with valence
                ! lr408: states
                if (pub_debug) then
                    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: print traces'
                    debug_trace = sparse_trace(unproj_ham%m(isub,isub))
                    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
                        &H_{AA}               = ',debug_trace
                    debug_trace = sparse_trace(ham%m(isub,isub))
                    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
                        &H_{AA}^proj          = ',debug_trace
                    call sparse_embed_trace(debug_trace,hab_pba)
                    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
                        &H_{AB}*P_{BA}        = ',debug_trace
                    call sparse_embed_trace(debug_trace,pab_hba)
                    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
                        &P_{AB}*H_{BA}        = ',debug_trace
                    call sparse_embed_trace(debug_trace,pab_hbb_pba)
                    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
                        &P_{AB}*H_{BB}*P_{BA} = ',debug_trace
                    call sparse_embed_trace(debug_trace,denskern)
                    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
                        &K                    = ',debug_trace

                    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: more traces...'
                    debug_trace = sparse_trace(denskern%m(isub,isub), unproj_ham%m(isub,isub))
                    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
                        &K_AA*H_{AA}               = ',debug_trace
                    debug_trace = sparse_trace(denskern%m(isub,isub), ham%m(isub,isub))
                    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
                        &K_AA*H_{AA}^proj          = ',debug_trace
                    debug_trace = sparse_trace(hab_pba%p)
                    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
                        &H_{AB}*P_{BA}             = ',debug_trace
                    debug_trace = sparse_trace(pab_hba%p)
                    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
                        &P_{AB}*H_{BA}             = ',debug_trace
                    debug_trace = sparse_trace(denskern%m(isub,isub), pab_hbb_pba%p)
                    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
                        &K_AA*P_{AB}*H_{BB}*P_{BA} = ',debug_trace
                    debug_trace = sparse_trace(denskern%m(isub,isub))
                    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
                        &K_{AA}                    = ',debug_trace
                    debug_trace = sparse_trace(denskern%m(jsub,jsub))
                    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
                        &K_{BB}                    = ',debug_trace
                end if

                ! Clean up temporary matrices
                call sparse_embed_destroy(sab_kbb)
                call sparse_embed_destroy(kbb_sba)
                call sparse_embed_destroy(hab_pba)
                call sparse_embed_destroy(pab_hba)
                call sparse_embed_destroy(pab_hbb)
                call sparse_embed_destroy(pab_hbb_pba)
            end if
        end do
    end do

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
         &Leaving hamiltonian_proj_embed_matrix'

  end subroutine hamiltonian_proj_embed_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end module hamiltonian
