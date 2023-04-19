! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!            Electronic energy optimisation module               !
!                                                                !
! This module contains routines used in the optimisation of the  !
! electronic energy with respect to the NGWFs and the density    !
! kernel.                                                        !
!----------------------------------------------------------------!
! Originally written by Chris-Kriton Skylaris in 2000.           !
! Subsequent modifications by Chris-Kriton Skylaris,             !
! Arash A. Mostofi, Peter D. Haynes and Nicholas D.M. Hine.      !
! Modified by Nicholas Hine to use the NGWF_REP container for    !
! NGWF-related data and matrices, October 2010.                  !
!================================================================!


module electronic

  use constants, only: DP
  use rundat, only: pub_debug_on_root

  implicit none

  private

  public :: electronic_energy
  public :: electronic_lagrangian

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  real(kind=DP) function electronic_energy(denskern, pur_denskern, ham, &
       lhxc_fine, mu, edft, rep, ngwf_basis, hub_proj_basis, hub, &
       mdl, hfxstate, lnv_threshold, current_maxit_lnv, kernel_update, &
       conv_status, dfdtau_fine, kpt, force_no_IO, mermin_threshold, &
       current_maxit_mermin)

    !=============================================================!
    ! This function, given a set of NGWFs, returns the total      !
    ! energy after first optimising the density kernel in the     !
    ! LNV scheme (if required).                                   !
    !-------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in June 2003.              !
    ! Spin polarised by Peter Haynes, July 2006                   !
    ! Modified for NLCC by Nicholas Hine, January 2009            !
    ! DFT+U added by David O'Regan, April 2008                    !
    ! Adapted for SPAM3 and function basis by Nicholas Hine,      !
    ! May-July 2009.                                              !
    ! Adapted for NGWF_REP by Nicholas Hine in October 2010.      !
    ! Modified to calculated the conduction energy by Laura       !
    ! Ratcliff in October 2010.                                   !
    ! Kernel DIIS added by Alvaro Ruiz Serrano, November 2010.    !
    ! EDA SCF-MI added by Max Phipps, November 2015.              !
    !=============================================================!

    use constants, only: DP, max_spins, paw_en_size, stdout, k_B
    use ensemble_dft, only: edft_calculate, edft_fermi_level, &
         edft_reorthogonalise_mo, edft_matrix_from_eigs, edft_diag_ngwf_ham
    use ensemble_dft_type, only: EDFT_MODEL
    use ensemble_dftb, only: edftb_calculate
    use fragment_data, only: eda_dfdtau_fine
    use function_basis, only: FUNC_BASIS
    use fragment_scfmi, only: scfmi_construct_nonorth_kernel, denskern_R
    use hamiltonian, only: hamiltonian_dens_dep_matrices, &
         hamiltonian_build_matrix, hamiltonian_proj_embed_matrix
    use hf_exchange, only: HFX_STATE
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: DKERN, kernel_purify, kernel_rescale, kernel_occ_check
    use kernel_diis, only: kernel_diis_calculate, &
         kernel_diis_build_idemp_dkn
    use lnv, only: lnv_denskernel_optimise_cg, lnv_calculate_mu
    use mermin, only: mermin_denskernel_optimise_cg, mermin_calculate_mu
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_exact_lnv, pub_cond_calculate, pub_kernel_diis, &
         pub_eda_scfmi, pub_dftb, &
         pub_edft, pub_edft_maxit, pub_maxit_kernel_occ_check, &
         pub_edft_smearing_width, pub_num_spins, pub_spin_fac, pub_foe, &
         pub_num_kpoints, PUB_1K, pub_xc_ke_density_required, &
         pub_edft_spin_fix, pub_mermin_smearing_width, pub_mermin_cheb, &
         pub_mermin
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_trace, sparse_embed_copy, &
         sparse_embed_scale, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_transpose, sparse_embed_product, sparse_embed_array_copy, &
         SPAM3_EMBED_ARRAY, sparse_embed_array_trace, &
         sparse_embed_array_create, sparse_embed_array_product
    use smearing_operator, only: calculate_fermi_entropy_mermin
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments (inputs)
    type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: hub_proj_basis(:)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(MODEL), intent(in) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    real(kind=DP), intent(in) :: lnv_threshold
    integer, intent(in) :: current_maxit_lnv
    logical, intent(in) :: kernel_update
    ! agrecokpt: needed to add extra terms in TB method
    real(kind=DP), optional, intent(in) :: kpt(3)

    ! Arguments (outputs)
    type(DKERN), intent(inout) :: denskern
    type(DKERN), intent(inout) :: pur_denskern
    type(NGWF_HAM), intent(inout) :: ham
    real(kind=DP), intent(inout) :: lhxc_fine(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, rep%nsub)
    real(kind=DP), intent(inout) :: mu(pub_num_spins, pub_num_kpoints)
    type(EDFT_MODEL), intent(inout) :: edft
    integer, optional, intent(out) :: conv_status
    ! JCW: Optional gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density. Required to evaluate lhxc matrix elements for meta-GGAs.
    real(kind=DP), optional, intent(inout) :: dfdtau_fine(:,:,:,:)
    ! JCW: In general, we have
    ! JCW:    dfdtau_fine(mdl%fine_grid%ld1,&
    ! JCW:      mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    ! JCW: or
    ! JCW:    dfdtau_fine(1,1,1,1)
    ! JCW: when it is unneeded.
    logical, optional, intent(in) :: force_no_IO
    real(kind=DP), optional, intent(in) :: mermin_threshold
    integer, optional, intent(in) :: current_maxit_mermin

    ! Local variables
    real(kind=DP) :: hubbard_energy
    real(kind=DP) :: lhxc_energy
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    real(kind=DP) :: trace
    integer :: is
    integer :: cvstat
    integer :: occ_check_iter
    logical :: kernel_reset
    logical :: ham_update
    type(SPAM3_EMBED),allocatable :: tmp_spam(:)
    integer :: ierr
    ! agrecokpt
    real(kind=DP) :: loc_kpt(3)
    ! kkbd: to hold non-integer occs for dkern rescale
    real(kind=DP),dimension(pub_num_spins, PUB_1K) :: t_occ
    logical :: loc_force_no_IO
    !ep
    real(kind=DP) :: smermin
    real(kind=DP) :: tmp


    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Function electronic_energy not ready yet for more&
         & than one k-point.')

    ! agrecokpt: default is Gamma point
    if (present(kpt)) then
       loc_kpt(:) = kpt(:)
    else
       loc_kpt(:) = 0.0_DP
    end if

    loc_force_no_IO=.false.
    if(present(force_no_IO)) then
       loc_force_no_IO=force_no_IO
    end if

    if (.not.pub_xc_ke_density_required) then
       ! JCW: Non-tau-dependent XC functional
       ! JCW: dfdtau_fine is either not present, or has a length of 1
       if ( present(dfdtau_fine ) ) then
          call utils_assert( &
               size(dfdtau_fine) == 1,&
               "Error in electronic_energy: &
               &pub_xc_ke_density_required is false, but dfdtau_fine is present &
               &and has a size > 1.")
       end if
    else
       ! JCW: tau-dependent meta-GGA XC functional
       ! JCW: Check dfdtau_fine has been passed in.
       call utils_assert(present(dfdtau_fine),"Error in &
            &electronic_energy: pub_xc_ke_density_required is true &
            &but no dfdtau_fine array has been passed in.")
    end if

    ! Old version: smmd: initialise convergence status to -1
    ! dhpt: default should be that convergence is not reached
    ! dhpt: need to check if this works in all situations
    !if (present(conv_status)) conv_status = -1
    if (present(conv_status)) conv_status = 1

    ! ndmh: energy evaluation for COND calculations
    if (pub_cond_calculate) then
       ! JCW: Abort if tau-dependent XC functional requested
       call utils_assert(.not.pub_xc_ke_density_required, "Error in &
            &electronic energy: pub_xc_ke_density_required and &
            &pub_cond_calculate are true, but this combination is not &
            &implemented.")
       if ((current_maxit_lnv > 0).and.(.not.pub_kernel_diis).and.(.not.pub_edft)) then

          if (kernel_update) then

             ! ndmh: loop for kernel occupancy checking if required
             occ_check_iter = 1
             do

                ! agrecokpt: in TB method, need to include k-dependent terms
                ! in internal_energy when updating density-dependent matrices
                call lnv_denskernel_optimise_cg( &
                     denskern, pur_denskern, ham, lhxc_fine, mu, &
                     electronic_energy, rep, ngwf_basis, hub_proj_basis, hub, &
                     mdl, hfxstate, lnv_threshold, current_maxit_lnv, &
                     cvstat, kpt=loc_kpt, force_no_IO=loc_force_no_IO)

                ! ndmh: quit if we have already done more cycles than required
                if (occ_check_iter > pub_maxit_kernel_occ_check) &
                   exit

                ! ndmh: kernel not reset unless kernel_occ_check says so for at
                ! ndmh: least one spin
                kernel_reset = .false.
                if (pub_maxit_kernel_occ_check>0) then
                   ! ndmh: compare occupancy of kernel to that based on
                   ! ndmh: canonical purification
                   call kernel_occ_check(denskern%kern, kernel_reset, &
                        ham%ham, rep%overlap, rep%inv_overlap, rep%n_occ)
                end if

                ! ndmh: no need to continue if kernel was not reset. If it was, go
                ! ndmh: round again and re-do LNV.
                if (.not.kernel_reset) exit
                occ_check_iter = occ_check_iter + 1

             end do

          end if

          ! smmd: If lnv_threshold is reached then conv_status = 0,
          !       else conv_status .gt. 0
          if (present(conv_status)) conv_status = cvstat

          ! cks: purify denskern and scale before Lagrangian evaluation
          call kernel_purify(pur_denskern%kern, denskern, rep%overlap, &
               rep%inv_overlap, rep%n_occ)

          ! Set normalisation factor for revised LNV
          if (pub_exact_lnv) call kernel_rescale(pur_denskern, rep%overlap, &
               rep%n_occ, silent=.true.)

       else if (pub_kernel_diis.or.(current_maxit_lnv < 0)) then

          ! ars: build idempotent kernel - scale appropriately
          call kernel_diis_build_idemp_dkn(pur_denskern%kern%m(:,PUB_1K), &
               ham%ham, rep%overlap, rep%n_occ)

       else if (pub_edft) then

          ! KPOINTS_DANGER: we are only considering 1 k-point here!!
          ! ars: if edft -> build kernel with smeared occupancies
          do is = 1, pub_num_spins

             ! ars: diagonalise ham and store e_i and M^\alpha_i
             call edft_diag_ngwf_ham(ham%ham(is), rep%overlap, edft%num, &
                  edft%mo(is), edft%h_evals(:,is))

             ! ars: check orthogonality of MOs and re-orthogonalise if necessary
             call edft_reorthogonalise_mo(edft%mo(is), denskern%kern%m(is,PUB_1K), &
                  ham%ham(is), rep%overlap, edft%num, edft%orth_residue(is))

          end do

          ! ars: build smeared occupancies (diagonal)
          call edft_fermi_level(edft%occ, edft%fermi, &
               edft%integrated_ne, edft%h_evals, edft%s_occ(:), &
               edft%num, edft%nbands, edft%nelec, &
               pub_edft_smearing_width, pub_edft_spin_fix)

          do is = 1, pub_num_spins
             ! ars: build density kernel
             call edft_matrix_from_eigs(denskern%kern%m(is,PUB_1K), edft%mo(is), &
                  edft%occ(:,is), edft%num)
          end do

          call sparse_embed_array_copy(pur_denskern%kern, denskern%kern)

       end if

       ! lr408: Take trace of kernel with the projected Hamiltonian
       ! KPOINTS_DANGER: only 1 k-point considered!!
       electronic_energy = 0.0_DP
       do is = 1, pub_num_spins
          call sparse_embed_trace(trace, ham%ham(is), pur_denskern%kern%m(is,PUB_1K))
          electronic_energy = electronic_energy + trace * pub_spin_fac
       end do

       ! rab207: print total energy obtained with final kernel
       ! only print if kernel has been updated
       if (kernel_update) call internal_print_final_state

       ! ndmh: skip the rest of this routine (non-COND code)
       return

    end if

    ! cks: optimise density kernel L using an LNV variant
    if ((current_maxit_lnv > 0).and.(kernel_update)&
         .and.(.not.pub_kernel_diis).and.(.not.pub_edft) &
         .and. (.not.pub_mermin)) then

       ! ndmh: loop for kernel occupancy checking if required
       occ_check_iter = 1
       do
          ! agrecokpt: in TB method, need to include k-dependent terms
          ! in internal_energy when updating density-dependent matrices
          call lnv_denskernel_optimise_cg( &
               denskern, pur_denskern, ham, lhxc_fine, mu, &
               electronic_energy, rep, ngwf_basis, hub_proj_basis, hub, &
               mdl, hfxstate, lnv_threshold, current_maxit_lnv, &
               cvstat, dfdtau_fine = dfdtau_fine, kpt=loc_kpt, &
               force_no_IO=loc_force_no_IO)

          ! ndmh: quit if we have already done more cycles than required
          if (occ_check_iter > pub_maxit_kernel_occ_check) exit

          ! ndmh: kernel not reset unless kernel_occ_check says so for at
          ! ndmh: least one spin
          kernel_reset = .false.
          if (pub_maxit_kernel_occ_check>0) then
             ! ndmh: ccompare occupancy of kernel to that based on
             ! ndmh: canonical purification
             call kernel_occ_check(denskern%kern, kernel_reset, &
                  ham%ham, rep%overlap, rep%inv_overlap, rep%n_occ)
          end if

          ! ndmh: no need to continue if kernel was not reset. If it was, go
          ! ndmh: round again and re-do LNV.
          if (.not.kernel_reset) exit
          occ_check_iter = occ_check_iter + 1

       end do

       ! smmd: If lnv_threshold is reached then conv_status = 0,
       !       else conv_status .gt. 0
       if (present(conv_status)) conv_status = cvstat

       ! cks: purify denskern and scale before Lagrangian evaluation
       call kernel_purify(pur_denskern%kern, denskern, rep%overlap, &
            rep%inv_overlap, rep%n_occ)

       ! Set normalisation factor for revised LNV
       if (pub_exact_lnv) call kernel_rescale(pur_denskern,rep%overlap, &
            rep%n_occ,silent=.true.)

    !ep
    else if (pub_mermin.and.kernel_update) then

       call mermin_denskernel_optimise_cg( &
            denskern, pur_denskern, ham, lhxc_fine, mu, &
            electronic_energy, rep, ngwf_basis, hub_proj_basis, hub, &
            mdl, hfxstate, mermin_threshold, current_maxit_mermin, &
            cvstat, dfdtau_fine = dfdtau_fine, kpt=loc_kpt, &
            force_no_IO=loc_force_no_IO)
       call sparse_embed_array_copy(pur_denskern%kern, denskern%kern)

    else if (pub_kernel_diis.and.kernel_update) then

       ! JCW: Abort if KE density required (meta-GGA + DIIS untested).
       call utils_assert(.not.pub_xc_ke_density_required, "Error in &
            &electronic energy: pub_xc_ke_density_required and &
            &pub_kernel_diis are true, but this combination is not &
            &implemented.")

       ! agrecokpt
       if (any(loc_kpt/=0.0_DP)) then
          if (pub_debug_on_root) write(stdout,'(a)') &
               'WARNING: kernel optimisation with DIIS in &
               &electronic_energy not yet checked for &
               &non-Gamma k-point'
       end if

       ! ars: optimise the density kernel using kernel DIIS
       electronic_energy = 0.0_DP
       call kernel_diis_calculate(electronic_energy, mu, pur_denskern%kern, &
            denskern%kern, ham, lhxc_fine, ngwf_basis, hub_proj_basis, hub, &
            rep, mdl, hfxstate)

    else if (pub_edft.and.kernel_update) then

       ! JCW: Abort if KE density required (meta-GGA + eDFT untested).
       call utils_assert(.not.pub_xc_ke_density_required, "Error in &
            &electronic energy: pub_xc_ke_density_required and &
            &pub_edft are true, but this combination is not &
            &implemented.")

       ! agrecokpt
       if (any(loc_kpt/=0.0_DP)) then
          if (pub_debug_on_root) write(stdout,'(a)') &
             'WARNING: kernel optimisation with edft in &
             &electronic_energy not yet checked for &
             &non-Gamma k-point'
       end if

       ! ars: optimise the density kernel using ensemble-DFT
       electronic_energy = 0.0_DP
       if (pub_dftb) then
          call edftb_calculate(edft, denskern%kern, ham, rep, &
               ngwf_basis, mdl, pub_edft_maxit,&
               force_no_IO=loc_force_no_IO)
       else
          call edft_calculate(edft, denskern%kern, ham, rep, lhxc_fine, &
               ngwf_basis, hub_proj_basis, hub, mdl, hfxstate, pub_edft_maxit,&
               force_no_IO=loc_force_no_IO)
       end if

       call sparse_embed_array_copy(pur_denskern%kern, denskern%kern)

       ! ars: we are minimising the free energy functional
       electronic_energy = edft%free_energy

    else ! ndmh: no density kernel optimisation

       ! ars: set flag for ham update
       ham_update = .false.
       if ((pub_kernel_diis).or. &
            ((current_maxit_lnv==0).and.kernel_update) .or. &
             (pub_eda_scfmi .eqv. .true.)) ham_update = .true.

       if ((.not.pub_kernel_diis).and.(.not.pub_edft) &
          .and.(.not.pub_mermin)) then

          ! cks: purify denskern
          call kernel_purify(pur_denskern%kern, denskern, rep%overlap, &
               rep%inv_overlap, rep%n_occ, fixed_denskern=.true.)

          ! Set normalisation factor for revised LNV
          if (pub_exact_lnv) call kernel_rescale(pur_denskern, rep%overlap, &
               rep%n_occ, silent=.true.)

       else if (pub_kernel_diis) then

          ! ars: rescale density kernel and pur_denskern
          call kernel_rescale(denskern, rep%overlap, rep%n_occ)
          call sparse_embed_array_copy(pur_denskern%kern, denskern%kern)

       else if (pub_edft) then !(ars: assumes ngwf_cg_rotate=T)
          if(pub_foe) then
             allocate(tmp_spam(2),stat=ierr)
             call utils_alloc_check('electronic_energy','tmp_spam',ierr)
             call sparse_embed_create(tmp_spam(1),ham%ham(1), edft%rot)
             call sparse_embed_create(tmp_spam(2),edft%rot)
          end if


          ! ars: rescale rotated kernel, copy to pur_denskern and rebuild Hamiltonian
          ! WARNING will only work at gamma
          t_occ(:,PUB_1K) = edft%integrated_ne(:)
          call kernel_rescale(denskern, rep%overlap, t_occ)
          call sparse_embed_array_copy(pur_denskern%kern, denskern%kern)
          do is = 1, pub_num_spins
             if(.not.pub_foe) then
                ! ars: rebuild the Hamiltonian
                call edft_matrix_from_eigs(ham%ham(is), edft%mo(is), &
                     edft%h_evals(:,is), edft%num, rep%overlap)
             else
!                edft%rot^T * ham * edft%rot
!                call sparse_product(tmp_spam(1), ham%ham(is), edft%rot)
!                call sparse_transpose(tmp_spam(2), edft%rot)
!                call sparse_product(ham%ham(is), tmp_spam(2), tmp_spam(1))
                call sparse_embed_product(tmp_spam(1), edft%old_ham(is), edft%rot)
                call sparse_embed_transpose(tmp_spam(2), edft%rot)
                call sparse_embed_product(ham%ham(is), tmp_spam(2), tmp_spam(1))
             end if
          end do

          if(pub_foe) then
             call sparse_embed_destroy(tmp_spam(1))
             call sparse_embed_destroy(tmp_spam(2))
             deallocate(tmp_spam,stat=ierr)
             call utils_dealloc_check('electronic_energy','tmp_spam',ierr)
          end if

       !ep
       else if (pub_mermin) then
          call sparse_embed_array_copy(pur_denskern%kern, denskern%kern)

       end if

       ! mjsp: If fragment SCF-MI
       if (pub_eda_scfmi) then
          ! rc2013: nor is embedding
          call utils_assert(mdl%nsub .eq. 1, "Error in &
               &electronic energy: more than 1 region has been specified in &
               &species_ngwf_regions with pub_eda_scfmi, &
               &but this combination is not tested.")

          ! mjsp: update denskern_R: proper representation of the density kernel in the

          ! full overlap
          call scfmi_construct_nonorth_kernel(pur_denskern,rep)

          ! calculate the density using the denskern_R kernel
          ! agrecokpt: include k-point dependence here?
          call hamiltonian_dens_dep_matrices(ham, lhxc_fine, &
               electronic_energy, lhxc_energy, hubbard_energy, &
               paw_sphere_energies, rep, ngwf_basis, &
               hub_proj_basis, hub, denskern_R%kern, mdl, hfxstate, &
               ham_update, lhxc_fixed=.false., kpt=loc_kpt, &
               dfdtau_fine = eda_dfdtau_fine)

       else
          ! ndmh: calculate density dependent energies and matrices
          ! KPOINTS_DANGER: only 1st k-point passed!!
          ! agrecokpt: call with k-point dependence to include terms
          ! in TB method
          call hamiltonian_dens_dep_matrices(ham, lhxc_fine, &
               electronic_energy, lhxc_energy, hubbard_energy, &
               paw_sphere_energies, rep, ngwf_basis, &
               hub_proj_basis, hub, pur_denskern%kern, mdl, hfxstate, &
               ham_update, lhxc_fixed=.false.,dfdtau_fine=dfdtau_fine,&
               kpt=loc_kpt)
       end if

       ! ars: update Hamiltonian in case of kernel-DIIS or EDFT
       ! agrecokpt: each component already k-point dependent at this stage
       if (ham_update) call hamiltonian_build_matrix(ham,rep)

       ! ndmh: to ensure first NGWF iteration uses correct mu, we
       ! ndmh: need to recalculate it if restarting
       if (((.not.pub_kernel_diis).and.(.not.pub_edft)).and. &
           ((.not.pub_kernel_diis).and.(.not.pub_mermin)).and. &
            ((current_maxit_lnv==0).and.kernel_update)) then

          ! Calculate chemical potential mu for current LNV scheme
          call lnv_calculate_mu(mu, denskern, ham, rep, ngwf_basis)

       end if

       ! ars: use free energy in EDFT calculations
       if (pub_edft) then
          edft%energy = electronic_energy
          electronic_energy = edft%energy + edft%entropy
       end if

    endif

    ! rab207: print total energy obtained with final kernel
    ! only print if kernel has been updated
    if (kernel_update) call internal_print_final_state

    contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_print_final_state

      use constants, only: NORMAL, stdout
      use comms, only: pub_on_root
      use kernel, only: kernel_rms_err, kernel_commutator, &
           kernel_occupancy_bounds, kernel_middle_occupancy,&
           kernel_bandgap
      use rundat, only: pub_output_detail, pub_kernel_track_mid_occ, &
           pub_num_kpoints

      implicit none

      ! internal
      real(kind=DP), &
           dimension(pub_num_spins, pub_num_kpoints) :: max_occ, mid_occ, min_occ
      real(kind=DP) :: rms_occ_err
      real(kind=DP) :: commutator
      real(kind=DP) :: bandgap(pub_num_spins)
      integer :: ik, isp
      character(len=22) :: msg(pub_num_spins, pub_num_kpoints)
      logical :: is_lnv

      ! jme: KPOINTS_DANGER
      ! Some parts of this subroutine enforce a single k-point.
      call utils_assert(pub_num_kpoints == PUB_1K, &
           'Subroutine internal_print_final_state not checked yet for more&
           & than one k-point.')

      is_lnv = .not.(pub_edft.or.pub_kernel_diis)

      msg(:,:) = '(OK)'

      ! compute commutator
      commutator = kernel_commutator(denskern, ham%ham, rep%overlap, &
           rep%inv_overlap)

      ! compute bandgap and occupancy information only when LNV
      if (is_lnv) then

         ! compute estimated bandgap
         bandgap = kernel_bandgap(denskern, ham%ham, rep%overlap,&
           rep%inv_overlap)

         rms_occ_err = kernel_rms_err(denskern, rep%overlap)
         call kernel_occupancy_bounds(max_occ, min_occ, denskern, rep%overlap)
         if (pub_kernel_track_mid_occ) then
            call kernel_middle_occupancy(mid_occ, denskern, rep%overlap)
            do ik = 1, pub_num_kpoints
               do isp = 1, pub_num_spins
                  if ((mid_occ(isp,ik)>0.15_DP).and.(mid_occ(isp,ik)<0.85_DP)) &
                    msg(isp,ik) = '(possible degeneracy)'
               end do
            end do
         end if
      end if


      ! rab207: report summary of calculation
      if (pub_output_detail >= NORMAL .and. pub_on_root) then
         write(stdout,'(//a)') &
              '>>> Density kernel optimised for the current NGWF basis:'
         write(stdout,'(5x,a)') repeat('~',59)
         if (pub_edft) then
            write(stdout,'(5x,a,es22.14,a)') &
                 'Total free energy           = ', electronic_energy, '  Eh'
         else
            write(stdout,'(5x,a,es22.14,a)') &
                 'Total energy                = ', electronic_energy, '  Eh'
         endif
         if (is_lnv) then
            if (.not.pub_cond_calculate) then
               if (pub_num_spins == 1) then
                  write(stdout,'(5x,a,es12.4,a)') &
                          'Estimated bandgap           = ', bandgap(1), '  Eh'
               else
                  write(stdout,'(5x,a,es12.4,a)') &
                          'Estimated bandgap (UP spin) = ', bandgap(1), '  Eh'
                  write(stdout,'(5x,a,es12.4,a)') &
                          'Estimated bandgap (DN spin) = ', bandgap(2), '  Eh'
               end if
            else
               if (pub_num_spins == 1) then
                  write(stdout,'(5x,a,es12.4,a)') &
                          'Estimated gap above Kc      = ', bandgap(1), '  Eh'
               else
                  write(stdout,'(5x,a,es12.4,a)') &
                          'Estimated gap above Kc (UP) = ', bandgap(1), '  Eh'
                  write(stdout,'(5x,a,es12.4,a)') &
                          'Estimated gap above Kc (DN) = ', bandgap(2), '  Eh'
               end if
            end if
            write(stdout,'(5x,a,es12.4)') &
                 'RMS occupancy error         = ', rms_occ_err
         endif
         write(stdout,'(5x,a,es12.4)') &
              '[H,K] commutator            = ', commutator

         if (is_lnv) then
            ! KPOINTS_DANGER: check that the information printed out makes sense!!
            do ik = 1, pub_num_kpoints
               if (pub_num_kpoints > 1) then
                  write(stdout, '(5x, a, I3)') 'k-point number', ik
               end if
               if (pub_num_spins == 1) then
                  write(stdout,'(5x,a,f7.3,a,f7.3,a)') &
                       'Occupancy bounds            = [',min_occ(1,ik),':',max_occ(1,ik),']'
                  if (pub_kernel_track_mid_occ) then
                     write(stdout,'(5x,a,f6.3,2x,a)') &
                       'Occupancy closest to 0.5    = ', mid_occ(1,ik), trim(msg(1,ik))
                  end if
               else
                  write(stdout,'(5x,a,f7.3,a,f7.3,a)') &
                       'UP spin occupancy bounds    = [',min_occ(1,ik),':',max_occ(1,ik),']'
                  if (pub_kernel_track_mid_occ) then
                     write(stdout,'(5x,a,f6.3,2x,a)') &
                       'UP spin occ. closest to 0.5 = ', mid_occ(1,ik), trim(msg(1,ik))
                  end if

                  write(stdout,'(5x,a,f7.3,a,f7.3,a)') &
                       'DN spin occupancy bounds    = [',min_occ(2,ik),':',max_occ(2,ik),']'
                  if (pub_kernel_track_mid_occ) then
                     write(stdout,'(5x,a,f6.3,2x,a)') &
                       'DN spin occ. closest to 0.5 = ', mid_occ(2,ik), trim(msg(2,ik))
                  endif
               end if
            end do
         end if
         write(stdout,'(5x,a/)') repeat('~',59)
      end if

    end subroutine internal_print_final_state

  end function electronic_energy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function electronic_lagrangian(Einput, overlap, denskern, &
       ham, mu, n_occ, edft, muext)

    !========================================================!
    ! This function returns the value of the LNV Lagrangian. !
    !--------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in June 2003.         !
    ! Spin polarised by Peter Haynes, July 2006              !
    ! Kernel DIIS by Alvaro Ruiz Serrano, November 2010.     !
    !========================================================!

    use constants, only: DP, max_spins
    use ensemble_dft, only: edft_lagrangian
    use ensemble_dft_type, only: EDFT_MODEL
    use mermin, only: mermin_lagrangian
    use kernel_diis, only: kernel_diis_lagrangian
    use rundat, only: pub_exact_lnv, pub_kernel_diis, pub_edft, pub_num_spins, &
         pub_spin_fac, pub_mermin
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_trace
    use utils, only: utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: Einput
    type(SPAM3_EMBED),   intent(in) :: overlap
    type(SPAM3_EMBED),   intent(in) :: denskern(:)
    type(SPAM3_EMBED),   intent(in) :: ham(:)
    real(kind=DP), intent(in) :: mu(max_spins)
    integer,       intent(in) :: n_occ(max_spins)
    type(EDFT_MODEL), intent(in) :: edft ! ars
    real(kind=DP), optional, intent(in) :: muext(max_spins)

    ! Local variables
    integer :: is
    real(kind=DP) :: trace


    ! ars: W = electronic_lagrangian
    electronic_lagrangian = Einput

    if (pub_kernel_diis) then

       ! ars: W = E(K') - f*tr[H(K'SK'-K')], where K' = [Ne/tr(KS)]*K
       call kernel_diis_lagrangian(electronic_lagrangian, overlap, denskern, &
            ham)

    else if (pub_edft) then

       ! ars: A = E -TS - tr[HKS(XSX-X)] - mu*[\sum_i f_i - Ne]
       ! ars: note that edft%energy has been updated in electronic_energy
       electronic_lagrangian = edft_lagrangian(edft)

    !ep
    elseif (pub_mermin) then

       ! ep: A= E -T*s[{(K^**S_**)}]+ mu*(tr(K^ab S_ab)-N_e)
       electronic_lagrangian = mermin_lagrangian(electronic_lagrangian, &
                 denskern(:), overlap, muext, n_occ(:))


    elseif (pub_exact_lnv) then

       ! ars: W = E(K'), where K' = [Ne/tr(KS)]*K
       electronic_lagrangian = electronic_lagrangian

    elseif (.not.pub_exact_lnv) then

       ! ars: W = E(K) - f*mu*(tr(LS)-Nocc)
       do is=1,pub_num_spins
          call sparse_embed_trace(trace,denskern(is),overlap)
          electronic_lagrangian = electronic_lagrangian - &
               pub_spin_fac * mu(is) * (trace - real(n_occ(is),kind=DP))
       end do

    else

       ! ars: incorrect logic for kernel minimisation - ONETEP stop
       call utils_abort('Error in electronic_lagrangian: incorrect kernel &
            &optimisation method.')
    end if

  end function electronic_lagrangian


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module electronic
