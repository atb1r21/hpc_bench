! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=================================================================!
!                                                                 !
!                          Optics Module                          !
!                                                                 !
! This module calculates optical absorption spectra.              !
!-----------------------------------------------------------------!
! Relocated from properties and slightly re-organised by Laura    !
! Ratcliff in November 2011.                                      !
! Maintained and partly re-written to use dense matrices where    !
! possible and to calculate local decomposition of spectra, by    !
! Nicholas Hine in July 2012.                                     !
!=================================================================!

module optics

  use constants, only: DP
  use rundat, only: pub_debug_on_root

  implicit none

  private

  public :: optics_calculate_spectra
  public :: optics_pos_mat_els
  public :: optics_grad_mat_els

contains

  subroutine optics_calculate_spectra(cur_spin,eigs_dens,eigen_en, &
       rep,ngwf_basis,nl_projectors,proj_basis,mdl,num_opt_states)

    !==================================================================!
    ! This subroutine calculates optical absorption spectra using      !
    ! either position or momentum matrix elements (including the       !
    ! commutator between the non-local potential and the position      !
    ! operator).  Both the matrix elements and the imaginary component !
    ! of the dielectric function are written to file.                  !
    !------------------------------------------------------------------!
    ! Moved from properties_spectra, which was written by Laura        !
    ! Ratcliff in December 2010 and updated to include the output of   !
    ! the imaginary component of the dielectric function by Laura      !
    ! Ratcliff in October 2011.                                        !
    ! Significantly re-worked by Nicholas Hine in July 2012, including !
    ! addition of local spectra and removal of explicit trans_weights  !
    ! and trans_energies arrays (which could become very large!)       !
    !==================================================================!

    use constants, only: stdout
    use dense, only: DEM, dense_product, dense_create, &
         dense_convert, dense_destroy, dense_copy
    use function_basis, only: FUNC_BASIS
    use geometry, only: operator(.CROSS.), operator(.DOT.)
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_calc_mom_mat_els, pub_any_nl_proj, &
         pub_spectra_print_mat_els, pub_opt_smear, pub_paw, pub_spec_scissor_op, &
         pub_num_kpoints, PUB_1K
    use sparse_embed, only: SPAM3_EMBED,sparse_embed_create,sparse_embed_destroy
    use utils, only: utils_abort, utils_assert

    implicit none

    ! Arguments
    type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(rep%nsub)
    type(FUNC_BASIS), intent(in) :: proj_basis(rep%nsub)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(rep%nsub)
    type(MODEL), intent(in) :: mdl
    type(DEM), intent(in) :: eigs_dens
    real(kind=DP), intent(inout) :: eigen_en(:)
    integer, intent(in) :: cur_spin
    integer, intent(in), optional :: num_opt_states ! lr408: total num states to include

    ! Local variables
    type(DEM) :: r_eigs_dens(3)
    type(DEM) :: s_eigs_dens(1)
    type(DEM) :: pos_dens
    type(DEM) :: eigs_dens_c
    type(SPAM3_EMBED) :: pos_elements(3)
    integer :: homo, ieig, num_ngwfs
    integer :: xyz
    integer :: loc_num_opt_states
    character(len=6) :: file_type
    ! agrecocmplx
    logical :: loc_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &optics_calculate_spectra'

    ! jme: KPOINTS_DANGER
    ! Two instructions of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine optics_calculate_spectra not ready yet for more&
         & than one k-point.')

    ! Determine prefix from NGWF_REP postfix
    if (rep%postfix=='') then
       file_type='val'
    else if (rep%postfix=='j') then
       file_type='joint'
    else
       return
    end if

    num_ngwfs = sum(ngwf_basis(:)%num)

    ! agrecocmplx: need to take care of complex case as well
    loc_cmplx = pub_calc_mom_mat_els.or.eigs_dens%iscmplx

    ! Calculate highest occupied orbital for this spin (else use first)
    homo = rep%n_occ(cur_spin,PUB_1K)
    if (homo==0) homo = 1

    ! ndmh: create duplicate of eigs_dens which may be complex if
    ! ndmh: calculating momentum matrix elements
    ! agrecocmplx
    call dense_create(eigs_dens_c,num_ngwfs,num_ngwfs, &
         iscmplx=loc_cmplx)
    call dense_copy(eigs_dens_c,eigs_dens)

    ! lr408: Create sparse matrices for matrix elems in NGWF representation
    do xyz=1,3
       if (pub_calc_mom_mat_els) then
          if (pub_any_nl_proj) then
             call sparse_embed_create(pos_elements(xyz),rep%nonlocpot, &
                  iscmplx=.true.)
          else if (pub_paw) then
             call utils_abort('Error in optics_calculate_spectra: PAW grad &
                  &terms not yet implemented')
             call sparse_embed_create(pos_elements(xyz),rep%nonlocpot, &
                  iscmplx=.true.)
          else ! non-local commutator not needed, but still complex matrix elems
             call sparse_embed_create(pos_elements(xyz),rep%overlap,iscmplx=.true.)
          end if
       else ! position matrix elements only (all real)
          ! agrecocmplx: need to be complex when using complex NGWFs
          ! because of the routine used to actually compute them
          call sparse_embed_create(pos_elements(xyz),rep%ngwf_overlap, &
               iscmplx=loc_cmplx)
       end if
    end do

    ! Calculate position or momentum matrix elems in NGWF representation
    if (.not. pub_calc_mom_mat_els) then
       call optics_pos_mat_els(pos_elements,rep,ngwf_basis,rep,ngwf_basis,&
            rep%overlap,rep%ngwf_overlap,proj_basis, mdl)
    else
       call optics_grad_mat_els(pos_elements,rep,ngwf_basis,proj_basis, &
            nl_projectors,mdl)
    end if

    ! Calculate S|psi> in NGWF representation
    call dense_create(s_eigs_dens(1),num_ngwfs,num_ngwfs, &
         iscmplx=loc_cmplx)
    call dense_create(pos_dens,num_ngwfs,num_ngwfs, &
         iscmplx=loc_cmplx)
    call dense_convert(pos_dens,rep%overlap)
    call dense_product(s_eigs_dens(1),pos_dens,eigs_dens_c, &
         opA='T')

    ! Apply Scissor Operator
    if (pub_spec_scissor_op /= 0.0_DP) then
       do ieig=homo+1,num_ngwfs
          eigen_en(ieig) = eigen_en(ieig) + pub_spec_scissor_op
       end do
    end if

    ! Calculate joint DOS
    if (present(num_opt_states)) then
       loc_num_opt_states = num_opt_states
    else
       loc_num_opt_states = num_ngwfs - rep%n_occ(cur_spin,PUB_1K)
    end if
    if (pub_opt_smear > 0.0_DP) then
       ! rc2013: EMBED_FIX -- this requires 1 subsystem
       call optics_ldos_gp(eigen_en,homo,loc_num_opt_states,eigs_dens, &
            s_eigs_dens,cur_spin,mdl%par,mdl%elements,ngwf_basis(1),file_type)
    end if
    call dense_destroy(pos_dens)
    call dense_destroy(s_eigs_dens(1))

    ! Calculate r|psi> in NGWF representation
    do xyz=1,3
       call dense_create(r_eigs_dens(xyz),num_ngwfs,num_ngwfs, &
            iscmplx=loc_cmplx)
    end do
    call dense_create(pos_dens,num_ngwfs,num_ngwfs, &
         iscmplx=loc_cmplx)
    do xyz=1,3
       call dense_convert(pos_dens,pos_elements(xyz))
       call dense_product(r_eigs_dens(xyz),pos_dens,eigs_dens_c, &
            opA='T')
    end do

    ! Calculate imaginary part of dielectric function
    if (pub_opt_smear > 0.0_DP) then
       ! rc2013: EMBED_FIX -- this requires 1 subsystem
       call optics_ldos_gp(eigen_en,homo,loc_num_opt_states,eigs_dens, &
            r_eigs_dens(1:3),cur_spin,mdl%par,mdl%elements,ngwf_basis(1),file_type)
    end if
    call dense_destroy(pos_dens)

    ! Write out matrix elements if required
    if (pub_spectra_print_mat_els) call optics_print_mat_els(r_eigs_dens, &
         eigs_dens_c,eigen_en,mdl%cell,  cur_spin,num_ngwfs, &
         loc_num_opt_states,file_type)

    ! Cancel Scissor Operator
    if (pub_spec_scissor_op /= 0.0_DP) then
       do ieig=homo+1,num_ngwfs
          eigen_en(ieig) = eigen_en(ieig) - pub_spec_scissor_op
       end do
    end if

    ! Clean up sparse and dense matrices
    do xyz=3,1,-1
       call dense_destroy(r_eigs_dens(xyz))
    end do
    do xyz=3,1,-1
       call sparse_embed_destroy(pos_elements(xyz))
    end do
    call dense_destroy(eigs_dens_c)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving optics_calculate_spectra'

  end subroutine optics_calculate_spectra


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine optics_pos_mat_els(r_elements, bra_rep, bra_ngwf_basis, &
       ket_rep, ket_ngwf_basis, bra_ket_overlap, ngwf_bra_ket_overlap, &
       proj_basis, mdl, &
       first_order, bra_siGp_overlap, ket_siGp_overlap, direction, weight, axis)

    !==================================================================!
    ! This subroutine calculates position matrix elements for a given  !
    ! NGWF representation. This is only expected to be useful for      !
    ! systems where no NGWFs overlap with others from neigbouring      !
    ! periodic images of the simulation cell.                          !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff in December 2010, based on             !
    ! properties_polarisation.                                         !
    ! Modified for PAW by Nicholas Hine in December 2011.              !
    ! Modified by Tim Zuehlsdorff in November 2014 to allow for the    !
    ! use in LR-TDDFT.                                                 !
    ! Modified for embedding by Joseph Prentice, July 2018             !
    ! Updated to include changes for lr_phonons by Joseph Prentice,    !
    ! October 2018                                                     !
    !==================================================================!

    use augmentation, only: augmentation_pos, augmentation_overlap
    use basis, only: basis_location_func_wrt_cell, basis_ket_start_wrt_fftbox, &
         basis_point_wrt_box, basis_start_of_box_wrt_cell
    use comms, only: pub_my_proc_id
    use constants, only: stdout
    use datatypes, only: COEF
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(*), operator(+)
    use integrals, only: integrals_pos
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_aug
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_copy, sparse_embed_axpy, &
         sparse_embed_scale
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_destroy, &
         sparse_element_exists, sparse_get_element, sparse_put_element
    use utils, only: utils_alloc_check, utils_dealloc_check
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout)   :: r_elements(3)
    type(FUNC_BASIS), intent(in) :: bra_ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: ket_ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(NGWF_REP), intent(in) :: bra_rep
    type(NGWF_REP), intent(in) :: ket_rep
    type(SPAM3_EMBED), intent(in) :: bra_ket_overlap
    type(SPAM3_EMBED), intent(in) :: ngwf_bra_ket_overlap
    type(MODEL), intent(in) :: mdl
    logical, intent(in), optional     :: first_order
    type(SPAM3_EMBED), intent(in), optional :: bra_siGp_overlap
    type(SPAM3_EMBED), intent(in), optional :: ket_siGp_overlap
    integer, intent(in), optional :: direction(mdl%nat)
    real(kind=DP), intent(in), optional :: weight(mdl%nat)
    integer, intent(in), optional :: axis

    ! Local Variables
    type(SPAM3)         :: r_elements_bare(3)
    type(SPAM3_EMBED)   :: r_elements_aug(3)
    type(SPAM3_EMBED)   :: overlap
    type(POINT)   :: a1,a2,a3
    type(POINT)   :: r
    real(kind=DP) :: R_fft(3)
    ! agrecocmplx: use COEF type to switch between real and complex case
    type(COEF) :: r_el,o_el
    real(kind=DP) :: cbg1,cbg2,cbg3
#ifdef FFTBOX_OLD_POS
    integer       :: cs1,cs2,cs3
#endif
    integer       :: bs1,bs2,bs3
    integer       :: axmin,axmax
    integer       :: xyz,jngwf,loc_ingwf,ingwf
    ! agrecocmplx
    logical :: loc_cmplx, loc_first_order
    integer :: isub,jsub,ierr,iat
    integer, allocatable :: tmp_direction(:)
    real(kind=DP), allocatable :: tmp_weight(:)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering optics_pos_mat_els'

    ! Start timer
    call timer_clock('optics_pos_mat_els',1)

    if (present(first_order)) then
       loc_first_order = first_order
    else
       loc_first_order = .false.
    end if

    ! Allocate only one matrix if axis is specified
    if (present(axis)) then
       axmin = axis
       axmax = axis
    else
       axmin = 1
       axmax = 3
    end if

    ! agrecocmplx
    loc_cmplx = bra_ket_overlap%p%iscmplx
    r_el%iscmplx = loc_cmplx
    o_el%iscmplx = loc_cmplx

    ! jcap: loop over regions
    do isub=1,mdl%nsub
       do jsub=1,mdl%nsub

          ! Create temporary unaugmented matrices to hold result
          do xyz=1,3
             call sparse_create(r_elements_bare(xyz),ngwf_bra_ket_overlap%m(isub,jsub))
          end do

          ! calculate matrix elements < phi_a | r_op - R_fft | phi_b >
          ! make sure that it works for different NGWF species used as bras and kets
          call integrals_pos(r_elements_bare,&
               bra_rep%ngwfs_on_grid(isub),bra_ngwf_basis(isub),&
               ket_rep%ngwfs_on_grid(jsub),ket_ngwf_basis(jsub), &
               mdl%dbl_grid,mdl%cell,mdl%fftbox,1,axis)

          do xyz=1,3
             call sparse_copy(r_elements(xyz)%m(isub,jsub),r_elements_bare(xyz))
             call sparse_destroy(r_elements_bare(xyz))
          end do

       end do
    end do

    ! Add R_fft * S_ab to matrix elements
    ! NOT IMPLEMENTED YET: accounting for atoms moving across cell boundary
    ! ndmh: NB this involves an O(N^2) check - could use the
    ! ndmh: index of S to calculate which elements need shifting
    ! ndmh: but this is a one-time calculation and the prefactor is tiny

    ! Calculate vectors between grid points
    a1 = (1.0_DP / mdl%cell%total_pt1) * mdl%cell%a1
    a2 = (1.0_DP / mdl%cell%total_pt2) * mdl%cell%a2
    a3 = (1.0_DP / mdl%cell%total_pt3) * mdl%cell%a3

    ! Create bare (un-augmented) NGWF overlap
    call sparse_embed_create(overlap,bra_ket_overlap)
    if (pub_aug) then
       ! Create -1*<fa|pi>oij<pj|fb> and add it to <fa|S|fb>
       ! to give <fa|fb>
       ! jcap: loop over regions here - not sure this is right, but it
       ! will do for now
       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub
             call augmentation_overlap(overlap%m(isub,jsub),&
                  mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp, &
                  bra_rep%sp_overlap%m(isub,jsub),ket_rep%sp_overlap%m(isub,jsub))
          end do
       end do
       call sparse_embed_scale(overlap,-1.0_DP)
       call sparse_embed_axpy(overlap,bra_ket_overlap,1.0_DP)
    else
       call sparse_embed_copy(overlap,bra_ket_overlap)
    end if

    ! jcap: loop over regions
    do jsub=1,mdl%nsub
       ! Loop over all ngwfs on this proc
       do loc_ingwf=1,ket_ngwf_basis(jsub)%num_on_proc(pub_my_proc_id)
          ingwf = loc_ingwf + ket_ngwf_basis(jsub)%first_on_proc(pub_my_proc_id) - 1

#ifdef FFTBOX_OLD_POS
          ! Determine where tightbox begins wrt fftbox
          call basis_ket_start_wrt_fftbox(bs1,bs2,bs3,mdl%fftbox%total_pt1,&
               mdl%fftbox%total_pt2,mdl%fftbox%total_pt3,mdl%fftbox)

          ! Find start of tightbox of ingwf
          call basis_location_func_wrt_cell(cs1,cs2,cs3, &
               ket_ngwf_basis(jsub)%tight_boxes(loc_ingwf),mdl%cell)

          ! Find vector to origin of FFTbox
          r = real(cs1-bs1,kind=DP) * a1 &
               + real(cs2-bs2,kind=DP) * a2 &
               + real(cs3-bs3,kind=DP) * a3

#else
          ! Centre of function wrt fftbox in terms of grid spacings
          call basis_point_wrt_box(cbg1, cbg2, cbg3, &  ! output
               ket_ngwf_basis(jsub)%spheres(loc_ingwf)%centre, &
               mdl%fftbox%total_pt1, mdl%fftbox%total_pt2, &
               mdl%fftbox%total_pt3, mdl%cell, mdl%fftbox)

          ! Start of fftbox wrt cell in terms of grid-point number
          call basis_start_of_box_wrt_cell(bs1,bs2,bs3, &
               ket_ngwf_basis(jsub)%spheres(loc_ingwf)%centre,&
               cbg1,cbg2,cbg3,mdl%cell)

          ! Find vector to origin of FFTbox
          r = real(bs1-1,kind=DP) * a1 &
               + real(bs2-1,kind=DP) * a2 &
               + real(bs3-1,kind=DP) * a3
#endif

          R_fft(1) = r%X
          R_fft(2) = r%Y
          R_fft(3) = r%Z

          ! jcap: second loop over regions
          do isub=1,mdl%nsub
             ! Loop over all row ngwfs
             do jngwf=1,bra_ngwf_basis(isub)%num

                ! Test if these NGWFs overlap
                if (.not.sparse_element_exists(bra_ket_overlap%m(isub,jsub),&
                     jngwf,ingwf)) cycle

                ! Extract overlap element
                ! agrecocmplx
                if (loc_cmplx) then
                   call sparse_get_element(o_el%z,overlap%m(isub,jsub),&
                        jngwf,ingwf)
                else
                   call sparse_get_element(o_el%d,overlap%m(isub,jsub),&
                        jngwf,ingwf)
                end if

                ! Extract element from r_elements and shift by R_fft*o_el
                do xyz=axmin,axmax
                   ! agrecocmplx
                   if (loc_cmplx) then
                      call sparse_get_element(r_el%z,&
                           r_elements(xyz)%m(isub,jsub),jngwf,ingwf)
                      r_el%z = R_fft(xyz) * o_el%z + r_el%z
                      call sparse_put_element(r_el%z,&
                           r_elements(xyz)%m(isub,jsub),jngwf,ingwf)
                   else
                      call sparse_get_element(r_el%d,&
                           r_elements(xyz)%m(isub,jsub),jngwf,ingwf)
                      r_el%d = R_fft(xyz) * o_el%d + r_el%d
                      call sparse_put_element(r_el%d,&
                           r_elements(xyz)%m(isub,jsub),jngwf,ingwf)
                   end if
                end do

             end do
          end do
       end do
    end do

    ! ndmh: in PAW/USP formalism, add contribution from sphere part
    if (pub_aug) then

       ! ndmh: create position elements matrix with augmented overlap sparsity
       ! ndmh: and copy existing matrix into it, then recreate old matrix with
       ! ndmh: augmented sparsity, copy back to that and destroy temporary
       ! ndmh: matrix
       do xyz=axmin,axmax
          call sparse_embed_create(r_elements_aug(xyz),bra_ket_overlap)
          call sparse_embed_copy(r_elements_aug(xyz),r_elements(xyz))
          call sparse_embed_destroy(r_elements(xyz))
          call sparse_embed_create(r_elements(xyz),r_elements_aug(xyz))
          call sparse_embed_copy(r_elements(xyz),r_elements_aug(xyz))
          call sparse_embed_destroy(r_elements_aug(xyz))
       end do

       ! ndmh: augment the position matrix
       ! jcap: loop over regions - again not sure this is correct, but
       ! it will do for the moment
       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub

             do xyz=1,3
                call sparse_create(r_elements_bare(xyz),r_elements(xyz)%m(isub,jsub))
                call sparse_copy(r_elements_bare(xyz),r_elements(xyz)%m(isub,jsub))
             end do

             ! jcap: need temporary arrays to only pass in the
             ! appropriate directions and weights
             ! jcap: only if we have been given a direction/weight
             if (present(direction)) then
                allocate(tmp_direction(mdl%regions(isub)%par%nat),stat=ierr)
                call utils_alloc_check('optics_pos_mat_els','tmp_direction',ierr)
                do iat=1,mdl%regions(isub)%par%nat
                   tmp_direction(iat) = direction( &
                        mdl%regions(isub)%elements(iat)%global_atom_number)
                end do
             end if
             if (present(weight)) then
                allocate(tmp_weight(mdl%regions(isub)%par%nat),stat=ierr)
                call utils_alloc_check('optics_pos_mat_els','tmp_weight',ierr)
                do iat=1,mdl%regions(isub)%par%nat
                   tmp_weight(iat) = weight( &
                        mdl%regions(isub)%elements(iat)%global_atom_number)
                end do
             end if


             if (loc_first_order) then
                if (present(direction).and.present(weight)) then
                   call augmentation_pos(r_elements_bare,proj_basis(isub), &
                        mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp, &
                        bra_rep%sp_overlap%m(isub,jsub),ket_rep%sp_overlap%m(isub,jsub),&
                        first_order = .true., &
                        direction = tmp_direction, weight = tmp_weight)
                else
                   call augmentation_pos(r_elements_bare,proj_basis(isub), &
                        mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp, &
                        bra_rep%sp_overlap%m(isub,jsub),ket_rep%sp_overlap%m(isub,jsub),&
                        first_order = .true.)
                end if
                call augmentation_pos(r_elements_bare,proj_basis(isub), &
                     mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp,&
                     bra_rep%sp_overlap%m(isub,jsub), ket_siGp_overlap%m(isub,jsub))
                call augmentation_pos(r_elements_bare,proj_basis(isub), &
                     mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp,&
                     bra_siGp_overlap%m(isub,jsub), ket_rep%sp_overlap%m(isub,jsub))
             else
                call augmentation_pos(r_elements_bare,proj_basis(isub),&
                     mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp, &
                     bra_rep%sp_overlap%m(isub,jsub),ket_rep%sp_overlap%m(isub,jsub))
             end if

             do xyz=1,3
                call sparse_copy(r_elements(xyz)%m(isub,jsub),r_elements_bare(xyz))
                call sparse_destroy(r_elements_bare(xyz))
             end do

             if (present(direction)) then
                deallocate(tmp_direction,stat=ierr)
                call utils_dealloc_check('calc_FO_ham','tmp_direction',ierr)
             end if
             if (present(weight)) then
                deallocate(tmp_weight,stat=ierr)
                call utils_dealloc_check('calc_FO_ham','tmp_weight',ierr)
             end if
          end do
       end do

    end if

    call sparse_embed_destroy(overlap)

    ! Stop timer
    call timer_clock('optics_pos_mat_els',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving optics_pos_mat_els'

  end subroutine optics_pos_mat_els


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine optics_grad_mat_els(pos_elements,rep,ngwf_basis,&
         proj_basis,nl_projectors,mdl)

    !==================================================================!
    ! This subroutine calculates momentum matrix elements (including   !
    ! the commutator between the non-local potential and the position  !
    ! operator) between NGWFs for optical absorption spectra           !
    ! calculations.                                                    !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff in November 2011, based on             !
    ! properties_spectra, which was written by Laura Ratcliff in       !
    ! December 2010.                                                   !
    ! Modified for PAW by Nicholas Hine in December 2011.              !
    ! Modified for embedding by Joseph Prentice, July 2018             !
    !==================================================================!

    use augmentation, only: augmentation_grad!,aug_nonlocal_commutator_mat
    use constants, only: stdout
    use function_basis, only: FUNC_BASIS
    use geometry, only: operator(.CROSS.), operator(.DOT.)
    use integrals, only: integrals_grad
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use projectors, only: PROJECTOR_SET
    use pseudopotentials, only: pseudo_nonlocal_commutator_mat, &
         pseudo_nonlocal_com_mat_fd, pseudo_get_dij
    use rundat, only: pub_any_nl_proj, &
         pub_spec_nonloc_fin_diff, pub_calc_nonloc_comm, &
         pub_spec_cont_deriv, pub_aug
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_scale, sparse_embed_axpy, &
         sparse_embed_copy
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_scale, &
         sparse_copy, sparse_axpy
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout)   :: pos_elements(3)
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(NGWF_REP), intent(in) :: rep
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(:)
    type(MODEL), intent(in) :: mdl

    ! Local variables
    type(SPAM3)   :: pos_elements_rc(3), temp_spam(3)
    type(SPAM3_EMBED)   :: nonloc_com_elements(3)
    type(SPAM3_EMBED) :: dij
    real(kind=DP) :: delta
    integer :: xyz, isub, jsub
    ! agrecocmplx
    logical :: loc_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering optics_grad_mat_els'

    ! Start timer
    call timer_clock('optics_grad_mat_els',1)

    delta = pub_spec_nonloc_fin_diff

    ! agrecocmplx
    loc_cmplx = rep%overlap%p%iscmplx

    ! lr408: Use relationship q.<f|r|i> = q.{<f|p|i>/(im*w_fi) + <f|[V_NL,r]|i>/(hbar*w_fi)}
    ! lr408: to calculate momentum matrix elements rather than position

    ! lr408: Calculate grad mat els then commutator
    ! lr408: (NB both need dividing by w_fi and we are in atomic units, so can
    ! lr408: ignore constants)
    ! agrecocmplx: case when rep%overlap is real
    if (.not.loc_cmplx) then
       ! jcap: loop over regions
       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub

             do xyz=1,3
                call sparse_create(pos_elements_rc(xyz),rep%overlap%m(isub,jsub),&
                     iscmplx=.false.)
             end do

             call integrals_grad(pos_elements_rc, rep%ngwfs_on_grid(isub), &
                  ngwf_basis(isub), rep%ngwfs_on_grid(jsub), ngwf_basis(jsub), &
                  mdl%cell, mdl%fftbox)

             if (pub_calc_nonloc_comm) then
                do xyz=1,3
                   call sparse_scale(pos_elements_rc(xyz),-1.0_DP)
                end do
             end if
             do xyz=1,3
                call sparse_copy(pos_elements(xyz)%m(isub,jsub),pos_elements_rc(xyz))
             end do

             do xyz=3,1,-1
                call sparse_destroy(pos_elements_rc(xyz))
             end do

          end do
       end do
    ! agrecocmplx: compute directly in complex case
    else
       ! jcap: loop over regions
       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub

             do xyz=1,3
                call sparse_create(temp_spam(xyz),pos_elements(xyz)%m(isub,jsub))
                call sparse_copy(temp_spam(xyz),pos_elements(xyz)%m(isub,jsub))
             end do

             call integrals_grad(temp_spam, &
                  rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                  rep%ngwfs_on_grid(jsub), ngwf_basis(jsub), &
                  mdl%cell, mdl%fftbox)

             do xyz=1,3
                call sparse_copy(temp_spam(xyz),pos_elements(xyz)%m(isub,jsub))
                call sparse_destroy(temp_spam(xyz))
             end do

          end do
       end do
       if (pub_calc_nonloc_comm) then
          do xyz=1,3
             call sparse_embed_scale(pos_elements(xyz),-1.0_DP)
          end do
       end if
    end if

    ! ndmh: for PAW/USP, augment the grad operator matrix
    if (pub_aug) then
       ! rc2013: augmentation not set up with embedding
       if(mdl%nsub .gt. 1) call utils_abort('Error in optics_grad_mat_els: &
            &augmentation requested with more than 1 subsystem, but this &
            &combination is not implemented/tested.')
       !call augmentation_grad(pos_elements,mdl%pseudo_sp, &
       !     mdl%paw_sp,rep%sp_overlap%p)
       ! jcap: loop over regions - I'm not sure this is the correct
       ! thing to do, but it'll do for now
       !do isub=1,mdl%nsub
       !   do jsub=1,mdl%nsub

       ! rc2013: EMBED_FIX!
       do xyz=1,3
          call sparse_create(temp_spam(xyz),pos_elements(xyz)%m(isub,jsub))
          call sparse_copy(temp_spam(xyz),pos_elements(xyz)%m(isub,jsub))
       end do

       call augmentation_grad(temp_spam,mdl%regions(isub)%pseudo_sp, &
            mdl%regions(isub)%paw_sp,rep%sp_overlap%m(isub,jsub))

       do xyz=1,3
          call sparse_copy(temp_spam(xyz),pos_elements(xyz)%m(isub,jsub))
          call sparse_destroy(temp_spam(xyz))
       end do

       !   end do
       !end do
    end if

    if ((pub_any_nl_proj.or.pub_aug).and.(pub_calc_nonloc_comm)) then

       do xyz=1,3
          call sparse_embed_create(nonloc_com_elements(xyz),rep%nonlocpot, &
               iscmplx=.true.)
       end do

       ! jcap: need to precalculate dij for each set of species
       dij%structure = 'E'
       call sparse_embed_create(dij, iscmplx=.true.)
       do isub=1,mdl%nsub
          call pseudo_get_dij(dij%m(isub,isub),mdl%regions(isub)%pseudo_sp)
       end do

       if (.not. pub_spec_cont_deriv) then ! i.e. finite diff

          if (pub_any_nl_proj) then
             call pseudo_nonlocal_com_mat_fd(nonloc_com_elements, &
                  proj_basis,nl_projectors,ngwf_basis,rep%ngwfs_on_grid, &
                  rep%sp_overlap,mdl%cell,mdl%fftbox,dij, &
                  rep%nonlocpot,delta)
          else if (pub_aug) then
             call utils_abort('Error in optics_grad_mat_els: &
                  &aug_nonlocal_com_mat_fd does not exist')
          end if

       else ! continuous derivative

          if (pub_any_nl_proj) then
             call pseudo_nonlocal_commutator_mat(nonloc_com_elements, &
                  proj_basis,nl_projectors,ngwf_basis,rep%ngwfs_on_grid, &
                  rep%sp_overlap,mdl%cell,mdl%fftbox,dij,delta)
          else if (pub_aug) then
             !call aug_nonlocal_commutator_mat(nonloc_com_elements, &
             !     proj_basis, elements, dijhat, rho_ij, ngwf_basis, &
             !     rep%ngwfs_on_grid, rep%sp_overlap, delta)
          end if

       end if ! finite/cont diff

       call sparse_embed_destroy(dij)

       do xyz=1,3
          call sparse_embed_axpy(pos_elements(xyz), &
               nonloc_com_elements(xyz), (1.0_DP,0.0_DP))
       end do

       do xyz=3,1,-1
          call sparse_embed_destroy(nonloc_com_elements(xyz))
       end do

    end if ! non-local commutator present or not

    ! Stop timer
    call timer_clock('optics_grad_mat_els',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving optics_grad_mat_els'

  end subroutine optics_grad_mat_els


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine optics_ldos_gp(eigen_en, homo, norb_max, eigenvecs_dens, &
       op_eigenvecs_dens, is, par, elements, ngwf_basis, ham_type)

    !===================================================================!
    ! This subroutine prepares and outputs in simple txt format a       !
    ! Local Joint Density of States (LDOS) plot or a plot of the        !
    ! imaginary dielectric function, which are been generated           !
    ! with the Gaussian smearing method, using gamma point only         !
    ! energies. The groups over which to sum the LDOS are defined by    !
    ! pub_ldos_groups, which is set through the species_ldos_groups     !
    ! block in the input file. LDOS is calculated as a histogram of the !
    ! sum of all Gaussian-broadened energies multiplied by a coefficient!
    ! expressing the contribution of the NGWFs on each atom in the LDOS !
    ! group to that orbital. Energies in the output are in eV and so is !
    ! the half-width of the smearing Gaussians which is set by the      !
    ! input parameter ldos_smear.                                       !
    !-------------------------------------------------------------------!
    ! Written by Nicholas Hine on 15th February 2010, based on parts of !
    ! eigenstates_dos_gp which was written by Chris-Kriton Skylaris on  !
    ! 10/05/2006.                                                       !
    ! Modified to include Hubbard subspaces by David O'Regan, 10/2/2011.!
    ! Moved to eigenstates_mod by Nicholas Hine in November 2011.       !
    ! Adapted for spectra and JDOS by Nicholas Hine in July 2012.       !
    !===================================================================!

    use bibliography, only: bibliography_cite
    use comms, only: comms_barrier, pub_on_root
    use constants, only: DP, UP, DN, stdout, HARTREE_IN_EVS, PI
    use datatypes, only: FUNCTIONS, COEF, data_functions_alloc, &
         data_functions_dealloc, data_functions_dot
    use dense, only: DEM, dense_get_col
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_rootname, pub_opt_smear, pub_ldos_ngroups, &
         pub_ldos_group_nsp, pub_ldos_groups, pub_cond_calculate, &
         pub_calc_mom_mat_els, pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_open_unit_check, utils_close_unit_check, utils_abort

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis   ! Function basis for NGWFs
    real(kind=DP), intent(in) :: eigen_en(ngwf_basis%num)
    type(DEM), intent(in) :: eigenvecs_dens
    type(DEM), intent(in) :: op_eigenvecs_dens(:)
    integer, intent(in) :: norb_max
    integer, intent(in) :: homo
    type(PARAL_INFO), intent(in) :: par
    type(ELEMENT), intent(in) :: elements(par%nat)
    integer, intent(in) :: is                    ! Current spin
    character(len=*), intent(in) :: ham_type

    ! Local Variables
    real(kind=DP) :: alpha        ! Gaussian exponent
    real(kind=DP) :: delta_e      ! energy increment
    real(kind=DP) :: e_point      ! energy point
    real(kind=DP) :: cfac         ! common factor
    real(kind=DP) :: gnorm        ! Gaussian normalisation factor
    real(kind=DP) :: orb_coeff(3) ! Coefficient of orbital for this NGWF
    real(kind=DP) :: orb_coeff_tot(3) ! Coefficient of orbital for all NGWFs
    real(kind=DP), allocatable, dimension(:) :: energies ! buffer for energies in eV
    real(kind=DP), allocatable, dimension(:,:,:) :: histo_val ! histogram values
    real(kind=DP) :: en_offset    ! left and right offset for energy scale
    ! agrecocmplx: use FUNCTIONS type to switch between real and complex cases
    type(FUNCTIONS) :: eig_col_i
    type(FUNCTIONS), allocatable :: op_eig_col_f(:)
    integer, parameter :: histonum=2001 ! number of histogram points
    integer :: row  ! row counting index
    integer :: output_unit ! fortran output unit number
    integer :: ierr ! memory allocation error flag
    integer :: igroup
    integer :: iat
    integer :: isp
    integer :: ingwf
    integer :: icart,ncart
    integer :: io_status ! file access error flag
    integer :: iorb,forb ! orbital counting indices
    logical :: iat_in_group
    character(len=256) :: output_file  ! output file name buffer
    character(len=4) :: iat_orig_id
    ! agrecocmplx
    logical :: loc_cmplx
    type(COEF) :: temp_coef

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering optics_ldos_gp'

    call bibliography_cite('LDOS')

    ncart = size(op_eigenvecs_dens)

    ! agrecocmplx
    loc_cmplx = eigenvecs_dens%iscmplx

    ! ndmh: Allocate temporary arrays
    allocate(energies(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('optics_ldos_gp','energies',ierr)
    allocate(histo_val(histonum,1:ncart,0:pub_ldos_ngroups+1),stat=ierr)
    call utils_alloc_check('optics_ldos_gp','histo_val',ierr)
    ! agrecocmplx: allocate using appropriate routines
    call data_functions_alloc(eig_col_i,ngwf_basis%num,iscmplx=loc_cmplx)
    allocate(op_eig_col_f(ncart),stat=ierr)
    call utils_alloc_check('optics_ldos_gp','op_eig_col_f',ierr)
    do icart=1,ncart
       call data_functions_alloc(op_eig_col_f(icart),ngwf_basis%num, &
            iscmplx=loc_cmplx)
    end do

    ! cks: store eigen energies in buffer and convert to eV
    energies = eigen_en
    energies = energies*HARTREE_IN_EVS
    ! qoh: Initialise to avoid compiler warning
    output_unit = -1

    if (pub_on_root) then

       if (ncart==1) then
          if (pub_num_spins == 1) then
             write(stdout,'(a)')&
                  '================== Joint density of States (JDOS) calculation &
                  &=================='
          elseif (is == UP) then
             write(stdout,'(a)')&
                  '============ Joint density of States (JDOS) calculation &
                  &for UP spin ============='
          elseif (is == DN) then
             write(stdout,'(a)')&
                  '=========== Joint density of States (JDOS) calculation &
                  &for DOWN spin ============'
          endif
       else
          if (pub_num_spins == 1) then
             write(stdout,'(a)')&
                  '================ Imaginary component of the dielectric function &
                  &================'
          elseif (is == UP) then
             write(stdout,'(a)')&
                  '========== Imaginary component of the dielectric function &
                  &for UP spin ===='
          elseif (is == DN) then
             write(stdout,'(a)')&
                  '========== Imaginary component of the dielectric function &
                  &for DOWN spin ===='
          endif
       end if

       ! ndmh: find name for output files
       ! ndmh: Include possibility of various COND task ham_types
       write(output_file,*) trim(pub_rootname)
       if (pub_cond_calculate) write(output_file,*) &
            trim(output_file)//'_'//trim(ham_type)
       if ((pub_num_spins==2).and.(is == UP)) &
            write(output_file,*) trim(output_file)//'_'//'UP'
       if ((pub_num_spins==2).and.(is == DN)) &
            write(output_file,*) trim(output_file)//'_'//'DN'
       if (ncart==1) then
          write(output_file,*) trim(output_file)//'_JDOS.txt'
       else
          write(output_file,*) trim(output_file)//'_OPT_SPEC'//'.txt'
       end if
       output_file = adjustl(output_file)

       ! ndmh: print output notification
       write(stdout,'(3a)',advance ='no') 'Writing "', trim(output_file),'" ...'

       ! ndmh: get a free unit number
       output_unit = utils_unit()

       ! ndmh: open the output file
       open(unit=output_unit, form="formatted" ,file=trim(output_file), &
            action="write",iostat=io_status)
       call utils_open_unit_check('optics_ldos_gp','output_file', &
            io_status)

       ! ndmh: write key to column entries
       if (ncart==3) then
          write(output_unit,'(a,i4,a)',err=100) '# Column      ',1, &
               ': Transition Energy'
       else
          write(output_unit,'(a,i4,a)',err=100) '# Column ',1, &
               ': Transition Energy'
       end if
       do igroup=1,pub_ldos_ngroups
          if (ncart==3) then
             write(output_unit,'(a,i4,a,i4,a)',advance='no',err=100) &
                  '# Columns',(igroup-1)*3+2,'-',igroup*3+1,': Imag diel fn ('
          else
             write(output_unit,'(a,i4,a)',advance='no',err=100) &
                  '# Column ',igroup+1,': JDOS ('
          end if
          do isp=1,pub_ldos_group_nsp(igroup)
             write(output_unit,'(a)',advance='no',err=100) &
                  pub_ldos_groups(isp,igroup)
             if (isp<pub_ldos_group_nsp(igroup)) &
                  write(output_unit,'(a)',advance='no',err=100) ', '
          end do
          write(output_unit,'(a)')')'
       end do
       if (ncart==3) then
          write(output_unit,'(a,i4,a,i4,a)',advance='no',err=100) &
               '# Columns',pub_ldos_ngroups*3+2,'-',pub_ldos_ngroups*3+4, &
               ': Imag diel fn ('
       else
          write(output_unit,'(a,i4,a)',advance='no',err=100) &
               '# Column ',pub_ldos_ngroups+2,': JDOS ('
       end if
       write(output_unit,'(a)') 'total)'
       write(output_unit,'(a)') '#'

       ! ndmh: write column header line
       write(output_unit,'(a)',advance='no',err =100)'#  Energy (eV)  |'
       do igroup=1,pub_ldos_ngroups
          if (ncart==3) then
             write(output_unit,'(a,i3,a)',advance='no',err=100) &
                  'Imag diel fn grp',igroup,':1                    :2                    :3|'
          else
             write(output_unit,'(a,i4,a)',advance='no',err=100) &
                  '    JDOS group ',igroup,'  |'
          end if
       end do

       if (ncart==3) then
          write(output_unit,'(a)',err=100) ' Total Imag diel fn:1 &
               &                   :2                    :3|'
       else
          write(output_unit,'(a)',err=100) '          Total JDOS |'
       end if

    end if

    ! ndmh: initialisations

    ! ndmh: smearing Gaussian exponent (pub_opt_smear is the half-width in eVs)
    alpha = log(2.0_DP) / ((pub_opt_smear*HARTREE_IN_EVS)**2)
    gnorm = sqrt(alpha/PI)
    en_offset = pub_opt_smear * 30.0_dp
    delta_e = (2.0_DP*en_offset + energies(norb_max) - energies(1)) / &
         real(histonum-1, kind=DP)
    histo_val(:,:,:) = 0.0_DP
    ! ndmh: common factor (au's) is 2 pi / epsilon_0 + conversion to per eV
    cfac = (8.0_DP * PI**2) * HARTREE_IN_EVS

    ! ndmh: accumulate contributions to JLDOS for each relevant orbital pair
    do iorb=1,homo

       ! agrecocmplx: distinguish between real and complex case
       if (loc_cmplx) then
          call dense_get_col(eig_col_i%z,eigenvecs_dens,iorb)
       else
          call dense_get_col(eig_col_i%d,eigenvecs_dens,iorb)
       end if

       do forb=homo+1,norb_max

          ! agrecocmplx: distinguish between real and complex case
          if (loc_cmplx) then
             do icart=1,ncart
                call dense_get_col(op_eig_col_f(icart)%z,op_eigenvecs_dens(icart),forb)
             end do
          else
             do icart=1,ncart
                call dense_get_col(op_eig_col_f(icart)%d,op_eigenvecs_dens(icart),forb)
             end do
          end if

          ! ndmh: Loop over LDOS histogram points to set up exp(-alpha*(e-(ei-ej))^2)
          ! ndmh: array for each point, where ei is the current transition energy
          e_point = -en_offset
          do row=1,histonum
             histo_val(row,1,0) = &
                  gnorm*exp(-alpha*(e_point - (energies(forb)-energies(iorb)))**2)
             e_point = e_point + delta_e
          end do

          ! ndmh: Calculate contribution to total optical spectrum or total
          ! ndmh: JDOS. If we are doing the local optical spectrum then
          ! ndmh: we need to calculate the |<psi_c|r|psi_v>| element here
          ! ndmh: as half of the squared matrix element, since we do not
          ! ndmh: want to insert the projector twice
          if (ncart==3) then
             do icart=1,ncart
                ! agrecocmplx: deal with both real and complex cases
                temp_coef = data_functions_dot(eig_col_i,op_eig_col_f(icart))
                if (temp_coef%iscmplx) then
                   orb_coeff_tot(icart) = abs(temp_coef%z)
                else
                   orb_coeff_tot(icart) = temp_coef%d
                end if
                !orb_coeff_tot(icart) = sum(eig_col_i(1:ngwf_basis%num) * &
                !        op_eig_col_f(1:ngwf_basis%num,icart))

                ! ndmh: momentum matrix elements
                if (pub_calc_mom_mat_els) then

                   ! ndmh: calculate contribution to total spectrum
                   do row=1,histonum
                      histo_val(row,icart,pub_ldos_ngroups+1) = &
                           histo_val(row,icart,pub_ldos_ngroups+1) &
                           + histo_val(row,1,0) * orb_coeff_tot(icart)**2 &
                           / (eigen_en(forb)-eigen_en(iorb))**2 * cfac
                   end do

                   ! ndmh: include factors in orb_coeff_tot for local decomp.
                   orb_coeff_tot(icart) =  orb_coeff_tot(icart) &
                        / (eigen_en(forb)-eigen_en(iorb))**2 * cfac

                ! ndmh: position matrix elements
                else

                   orb_coeff_tot(icart) =  orb_coeff_tot(icart) * cfac
                   ! ndmh: calculate contribution to total spectrum
                   do row=1,histonum
                      histo_val(row,icart,pub_ldos_ngroups+1) = &
                           histo_val(row,icart,pub_ldos_ngroups+1) &
                           + histo_val(row,1,0) * orb_coeff_tot(icart)**2
                   end do

                end if
             end do
          else
             ! ndmh: add this i,f pair to total JDOS
             do row=1,histonum
                histo_val(row,1:ncart,pub_ldos_ngroups+1) = &
                        histo_val(row,1:ncart,pub_ldos_ngroups+1) + &
                        histo_val(row,1,0)
             end do
          end if

          ! ndmh: cycle if there are no groups other than total
          if (pub_ldos_ngroups==0) cycle

          ! ndmh: loop over atoms to add up contributions to the local DOS
          do iat=1,par%nat
             iat_orig_id = elements(par%orig_atom(iat))%species_id

             ! ndmh: go to next atom if this atom does not occur in any LDOS groups
             iat_in_group = .false.
             do igroup=1,pub_ldos_ngroups
                if (any(pub_ldos_groups(1:pub_ldos_group_nsp(igroup),igroup) &
                     ==iat_orig_id)) iat_in_group=.true.
             end do

             ! ndmh: loop over NGWFs on this atom
             do ingwf=ngwf_basis%first_on_atom(iat), &
                  ngwf_basis%first_on_atom(iat)+ngwf_basis%num_on_atom(iat)-1

                ! ndmh: evaluate orbital coefficient for NGWF a
                ! ndmh: given by M(dagger)_an * sum_b(S_ab M_nb) (LJDOS)
                ! ndmh: or by M(dagger)_an * sum_b(R_ab M_nb) (local spectra)
                ! agrecocmplx: deal with both real and complex cases
                if (loc_cmplx) then
                   do icart=1,ncart
                      orb_coeff(icart) = abs(conjg(eig_col_i%z(ingwf)) * &
                          op_eig_col_f(icart)%z(ingwf))
                   end do
                else
                   do icart=1,ncart
                      orb_coeff(icart) = eig_col_i%d(ingwf) * &
                          op_eig_col_f(icart)%d(ingwf)
                   end do
                end if

                ! ndmh: add contribution of this eigenvalue to any of the
                ! ndmh: groups in which this atom is listed
                if (iat_in_group) then
                   do igroup=1,pub_ldos_ngroups

                      if (any(pub_ldos_groups(1:pub_ldos_group_nsp(igroup),igroup) &
                           ==iat_orig_id)) then
                         if (ncart==3) then
                            do icart=1,ncart
                               histo_val(:,icart,igroup) = &
                                    histo_val(:,icart,igroup) + &
                                    orb_coeff_tot(icart) * orb_coeff(icart) * &
                                    histo_val(:,1,0)
                            end do
                         else
                            histo_val(:,1,igroup) = histo_val(:,1,igroup) + &
                                 orb_coeff(1) * histo_val(:,1,0)
                         end if
                      end if
                   end do  ! igroup
                endif ! iat_in_group

             end do  ! ingwf

          end do  ! iat
       end do ! forb
    end do  ! iorb

    if (pub_on_root) then

       ! ndmh: loop over histogram points writing output file
       e_point = -en_offset
       do row=1,histonum

          ! ndmh: write to file current histo-point for each group
          write(output_unit,'(f16.8)',advance='no',err=100) e_point
          do igroup=1,pub_ldos_ngroups
             do icart=1,ncart
                write(output_unit,'(f22.10)',advance='no',err=100) &
                     histo_val(row,icart,igroup)
              end do
          end do
          ! ndmh: write to file current histo-point for total
          do icart=1,ncart
             write(output_unit,'(f22.10)',advance='no',err=100) &
                  histo_val(row,icart,pub_ldos_ngroups+1)
          end do
          ! ndmh: end the line
          write(output_unit,*)

          e_point = e_point + delta_e

       end do

       ! ndmh: close output file
       close(unit=output_unit,iostat=io_status)
       call utils_close_unit_check('optics_ldos_gp','output_unit', &
            io_status)

       ! ndmh: finish writing message to stdout
       write(stdout,*)' done'
       write(stdout,'(a)')'================================&
            &================================================'

    endif

    ! ndmh: deallocate temporary arrays
    ! agrecocmplx: deallocate using appropriate routines
    do icart=1,ncart
       call data_functions_dealloc(op_eig_col_f(icart))
    end do
    deallocate(op_eig_col_f,stat=ierr)
    call utils_dealloc_check('optics_ldos_gp','op_eig_col_f',ierr)
    call data_functions_dealloc(eig_col_i)
    deallocate(histo_val,stat=ierr)
    call utils_dealloc_check('optics_ldos_gp','histo_val',ierr)
    deallocate(energies,stat=ierr)
    call utils_dealloc_check('optics_ldos_gp','energies',ierr)

    call comms_barrier

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving optics_ldos_gp'

    return

100 call utils_abort('Problem writing to file in optics_ldos_gp.')

  end subroutine optics_ldos_gp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine optics_print_mat_els(r_eigs_dens,eigs_dens, &
       eigen_en,cell,cur_spin,num_basis,num_opt_states,file_type)

    !==================================================================!
    ! This subroutine writes the optical matrix elements to a file.    !
    !------------------------------------------------------------------!
    ! Routine created by Nicholas Hine in July 2012 based on code      !
    ! by Laura Ratcliff, Nicholas Hine and Edgar Engel, originally     !
    ! written in December 2010 onward                                  !
    !==================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, HARTREE_IN_EVS, UP, DN
    use dense, only: DEM, dense_get_element, dense_product, dense_create, &
         dense_destroy
    use geometry, only: operator(.CROSS.), operator(.DOT.)
    use simulation_cell, only: CELL_INFO
    use rundat, only: pub_rootname, pub_calc_mom_mat_els, pub_num_spins
    use utils, only:  utils_unit, utils_open_unit_check, &
         utils_close_unit_check, utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(in) :: r_eigs_dens(3)
    type(DEM), intent(in) :: eigs_dens
    type(CELL_INFO), intent(in) :: cell
    real(kind=DP), intent(in) :: eigen_en(:) ! hamiltonian eigenvalues
    integer, intent(in) :: cur_spin
    integer, intent(in) :: num_basis
    integer, intent(in) :: num_opt_states
    character(len=6), intent(in) :: file_type

    ! Local variables
    type(DEM) :: opt_mat_elements ! Optical matrix elements of eigenvectors
    real(kind=DP) :: cell_volume
    real(kind=DP) :: opt_mat_el
    complex(kind=DP) :: opt_mat_el_cmplx
    integer :: xyz   ! cartesian direction
    integer :: i, f  ! initial and final state indices
    integer :: output_unit
    integer :: io_status
    character(len=256) :: output_file  ! output file name buffer

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering optics_print_mat_els'

    ! lr408: calculate cell volume
    cell_volume = abs((cell%a1 .CROSS. cell%a2) .DOT. cell%a3)
    call utils_assert(cell_volume /= 0.0_DP,'Cell volume is zero')

    if (pub_on_root) then

       write(stdout,'(a)') ''

       write(stdout,'(a)')&
            '================ Writing optical matrix elements &
            &================'
       if (pub_num_spins == 1) then
          write(output_file,*) trim(pub_rootname)//'_'// &
               trim(file_type)//'_OPT_MAT_ELS.txt'
       else if (cur_spin == UP) then
          write(output_file,*)trim(pub_rootname)//'_'// &
                trim(file_type)//'_OPT_MAT_ELS_UP.txt'
       else if (cur_spin == DN) then
          write(output_file,*)trim(pub_rootname)//'_'// &
               trim(file_type)//'_OPT_MAT_ELS_DN.txt'
       end if

       output_file = adjustl(output_file)

       ! cks: print output warning
       write(stdout,'(3a)') 'Writing "', trim(output_file),'" ...  done'

       write(stdout,'(a)')'================================&
            &================================================'
       write(stdout,'(a)') ''

       ! cks: get a unit number that is free
       output_unit = utils_unit()

       open(unit=output_unit, form="formatted" ,file=trim(output_file), &
            action="write",iostat=io_status)
       call utils_open_unit_check('optics_print_mat_els','output_file', &
            io_status)

       ! cks: write first line
       write(output_unit,'(a)')'# initial | Energy (eV) | final  | Energy (eV) &
           &|            Matrix el.  |        Trans. E (eV)   | Complex Matrix el.'

       write(output_unit,'(i5,2x,i5)') num_opt_states,num_opt_states,1
       write(output_unit,'(f24.12)') cell_volume

    end if

    ! Loop over cartesian directions
    do xyz=1,3

       if (pub_on_root) write(output_unit,'(a,i1)') &
            '# Cartesian component ',xyz

       ! Create storage for optical matrix elements in eigenstate basis
       call dense_create(opt_mat_elements,num_basis,num_basis, &
            iscmplx=pub_calc_mom_mat_els)

       ! Calculate optical matrix elements in eigenstate basis
       call dense_product(opt_mat_elements,eigs_dens,r_eigs_dens(xyz), &
            opA='T')

       ! Double loop over states
       do i=1,num_opt_states ! initial state
          do f=1,num_opt_states ! final state

             ! Write matrix element to file
             if (.not.pub_calc_mom_mat_els) then
                call dense_get_element(opt_mat_el,opt_mat_elements,i,f)
                if (pub_on_root) then
                   write(output_unit,'(2x,2(i5,4x,f12.6,2x),2(f24.12,2x))') &
                        i,eigen_en(i)*HARTREE_IN_EVS,&
                        f,eigen_en(f)*HARTREE_IN_EVS,&
                        opt_mat_el,&
                        (eigen_en(f)-eigen_en(i))*HARTREE_IN_EVS
                end if
             else ! momentum matrix elements
                call dense_get_element(opt_mat_el_cmplx,opt_mat_elements,i,f)
                if (pub_on_root) then
                   if (i /= f) then
                      ! eae32: set matrix elements to zero for like states
                      write(output_unit,'(2x,2(i5,4x,f12.6,2x),4(f24.12,2x))') &
                           i,eigen_en(i)*HARTREE_IN_EVS,&
                           f,eigen_en(f)*HARTREE_IN_EVS,&
                           abs(opt_mat_el_cmplx/(eigen_en(f)-eigen_en(i))),&
                           (eigen_en(f)-eigen_en(i))*HARTREE_IN_EVS,&
                           opt_mat_el_cmplx/(eigen_en(f)-eigen_en(i))
                   else
                      write(output_unit,'(2x,2(i5,4x,f12.6,2x),4(f24.12,2x))') &
                           i,eigen_en(i)*HARTREE_IN_EVS,&
                           f,eigen_en(f)*HARTREE_IN_EVS,&
                           0.0_DP,(eigen_en(f)-eigen_en(i))*HARTREE_IN_EVS,&
                           0.0_DP,0.0_DP
                   end if
                end if

             end if ! position / momentum matrix elements

          end do
       end do

       ! Clean up storage for optical matrix elements
       call dense_destroy(opt_mat_elements)

    end do ! loop over xyz

    ! Close output file for matrix elements
    if (pub_on_root) then
       close(unit=output_unit,iostat=io_status)
       call utils_close_unit_check('optics_print_mat_els','output_unit', &
            io_status)
    end if

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving optics_print_mat_els'

  end subroutine optics_print_mat_els


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




end module optics
