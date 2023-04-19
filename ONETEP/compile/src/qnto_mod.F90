! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!   Quantified Natural Transition Orbital Analysis Module        !
!                                                                !
! This module contains routines associated with the use of the   !
! Quantified Natural Transition Orbital method.                  !
!----------------------------------------------------------------!
! This module was created by Jian-Hao Li in 2013, with help from !
! Nicholas Hine.                                                 !
! Modified for embedding and to remove pub_par by Joseph         !
! Prentice, September 2018                                       !
!----------------------------------------------------------------!
! Default setting of input parameters:                           !
! pub_qnto_analysis : T                                          !
! pub_qnto_num_transition : 2 (limited to <= 5)                  !
! pub_qnto_write_orbitals : F                                    !
! pub_qnto_svd_method : 2                                        !
! pub_qnto_nbo_proj : F                                          !
! pub_qnto_ref_dir : ''                                          !
! pub_qnto_num_ref_states : 1                                    !
! pub_qnto_num_core_atoms : 0                                    !
!================================================================!

module qnto

  use constants, only: DP

  implicit none

  private

  public :: qnto_calculate

  type, public :: NBO_info

    integer :: nbo_idx
    character(len=3) :: type
    integer :: count
    character(len=2) :: a1e
    character(len=3) :: a1id
    character(len=2) :: a2e
    character(len=3) :: a2id
    character(len=2) :: a3e
    character(len=3) :: a3id
    real(kind=DP) :: occ
    real(kind=DP) :: energy

  end type NBO_info

contains

  subroutine qnto_calculate(response_kernel, exc_energies, val_rep, &
       val_ngwf_basis, val_denskern, cond_rep, cond_ngwf_basis, cond_denskern, &
       joint_rep, joint_ngwf_basis, joint_denskern, num_states, mdl, aux_rep, &
       aux_ngwf_basis, aux_denskern)

    !======================================================================!
    ! Subroutine that performs the QNTO analysis for each electronic       !
    ! excitation using the optimised response kernel obtained from running !
    ! the lr_tddft_cg subroutine.                                          !
    !----------------------------------------------------------------------!
    ! Originally written by Jian-Hao Li in August 2013.                    !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,     !
    ! September 2018                                                       !
    !======================================================================!

    use comms, only: comms_barrier, comms_bcast, pub_on_root, &
         pub_root_proc_id, pub_my_proc_id
    use constants, only: stdout, HARTREE_IN_EVS
    use datatypes, only: FUNCTIONS, data_functions_alloc, data_functions_dealloc
    use dense, only: DEM, dense_create, dense_destroy, dense_copy, &
         dense_convert, dense_eigensolve, dense_axpy, dense_scale, &
         dense_get_col, dense_put_col, &
         dense_product, dense_point2point_element, dense_get_element, &
         dense_put_element, dense_transpose, dense_mat_sqrt, dense_invert, &
         dense_svdsolve
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketppd
    use ion, only: ELEMENT
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, ngwf_rep_create
    use npa, only: npa_read_nbodata
    use restart, only: restart_response_kernel_batch_read, &
         restart_ngwfs_tightbox_input
    use rundat, only: pub_rootname, pub_use_aux_ngwfs, pub_qnto_svd_method, &
         pub_qnto_num_transition, pub_qnto_ref_dir, pub_lr_tddft_joint_set, &
         pub_qnto_num_ref_states, pub_qnto_nbo_proj, pub_num_spins, PUB_1K
    use sparse, only: sparse_get_element, sparse_put_element, sparse_get_col, &
         sparse_put_col, sparse_proc_of_elem
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product, sparse_embed_trace, &
         sparse_embed_axpy, sparse_embed_copy, sparse_embed_transpose, &
         sparse_embed_transpose_structure, sparse_embed_num_rows, sparse_embed_scale
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit,&
         utils_heapsort, utils_open_unit_check, utils_close_unit_check, &
         utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: num_states
    type(SPAM3_EMBED), intent(inout) :: response_kernel(num_states)
    real(kind=DP), intent(in) :: exc_energies(num_states)
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis
    type(SPAM3_EMBED_ARRAY), intent(inout) :: val_denskern
    type(NGWF_REP), intent(in) :: cond_rep
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis
    type(SPAM3_EMBED_ARRAY), intent(in) :: cond_denskern
    type(NGWF_REP), intent(in) :: joint_rep
    type(FUNC_BASIS), intent(in) :: joint_ngwf_basis
    type(SPAM3_EMBED_ARRAY), intent(in) :: joint_denskern
    type(MODEL), intent(inout) :: mdl
    type(NGWF_REP), optional, intent(in) :: aux_rep
    type(FUNC_BASIS), optional, intent(in) :: aux_ngwf_basis
    type(SPAM3_EMBED_ARRAY), intent(in) :: aux_denskern

    ! Local variables for calculating SVD to obtain NTOs
    type(SPAM3_EMBED) :: response_kernel_trans
    type(SPAM3_EMBED) :: ScK, SvKt, ScKSvKt, SvKtScK, ScKSvKtSc, SvKtScKSv
    type(SPAM3_EMBED) :: Tv, Tc
    type(DEM) :: dbpt_dens, dbmt_dens
    type(DEM) :: dbpSv_dens, dbptSc_dens, dbmSv_dens, dbmtSc_dens
    type(DEM) :: dbpSvdbptSc_dens, dbmSvdbmtSc_dens
    type(DEM) :: SvKtScKSv_dens, Sv_dens, ntoh_sys_dens
    type(DEM) :: ScKSvKtSc_dens, Sc_dens, ntoe_sys_dens
    type(DEM) :: Tv_dens, Tc_dens
    type(DEM) :: response_kernel_dens
    type(DEM) :: dbp_dens, dbm_dens
    type(DEM) :: db_prev_dens, tmp_dens
    real(kind=DP), allocatable :: lambda_t(:), lambda_tv(:), lambda_tc(:)
    real(kind=DP), allocatable :: trans_coeff(:)
    real(kind=DP) :: Dbp_norm=0.0, Dbm_norm=0.0
    real(kind=DP) :: tot_perc
    integer :: I
    integer :: num_diag, n
    integer :: ierr
    character(len=256) :: nbo_output_file, qnto_output_file
    character(len=256) :: ntoh_output_file, ntoe_output_file, ntoh_in_filename
    character(len=256) :: ntoe_in_filename, label_io_filename
    character(len=2) :: name_n

    ! Local variables for projecting NTOs to NBOs & for comparing NTOs from
    ! different molecules or the same molecule with different geometries
    type(DEM), allocatable :: nbo_sys_dens(:), stdnbo_sys_dens(:)
    logical :: pre_ao
    integer :: is, nspin ! # spins of orbitals in NBO's output files
    integer :: num_ao, num_nbo, num_stdnbo
    integer :: nbo_idx
    character(len=64) :: nbo_trfile
    type(DEM) :: Tvntoh_sys_dens, ntoh2nbo_dens
    type(DEM) :: Tcntoe_sys_dens, ntoe2nbo_dens
    type(NBO_info), allocatable :: nbo_list(:), stdnbo_list(:)
    character(len=4) :: apx
    type(SPAM3_EMBED) :: nbo, nbo_std
    logical :: std
    integer :: qnto_output_unit
    integer :: it, irow, icol
    type(DEM) :: nbo_reftgt_dens
    type(DEM) :: ntoh_systgt_dens, ntoh_reftgt_dens
    type(DEM) :: ntoe_systgt_dens, ntoe_reftgt_dens
    type(DEM) :: ntoh_sctgt_dens, ntoh_rctgt_dens
    type(DEM) :: ntoe_sctgt_dens, ntoe_rctgt_dens
    type(DEM) :: nto_sch2sch, nto_sce2sce, nto_sch2sce, nto_sce2sch
    type(DEM) :: nto_sh2sch, nto_se2sce, nto_sh2sce, nto_se2sch
    type(DEM) :: nto_rch2rch, nto_rce2rce, nto_rch2rce, nto_rce2rch
    type(DEM) :: nto_sch2rch, nto_sce2rce, nto_sch2rce, nto_sce2rch
    type(DEM) :: Svrs_dens, Scrs_dens, Stcrvr_dens, Stvrcr_dens
    integer :: systgt_size, reftgt_size, nto_io_idx, state_idx_sys
    integer :: atom_idx, ngwf_idx
    integer, allocatable :: val_ngwf_1st_idx_inp(:), val_ngwf_1st_idx_nod(:)
    integer, allocatable :: cond_ngwf_1st_idx_inp(:), cond_ngwf_1st_idx_nod(:)
    integer, allocatable :: joint_ngwf_1st_idx_inp(:), joint_ngwf_1st_idx_nod(:)
    integer, allocatable :: aux_ngwf_1st_idx_inp(:), aux_ngwf_1st_idx_nod(:)
    integer :: label_io_unit
    character(len=80) :: string, dir_sub
    integer :: num_dir, sub_it, idx_start, idx_end
    integer :: pos, pos_sub
    character(len=80), allocatable :: ref_dir(:)
    integer :: num_nto_trans
    integer :: num_min_nto_trans
    type(NGWF_REP) :: ref_nbo_rep, val_s2r_rep, val_s2s_rep, val_r2r_rep
    type(NGWF_REP) :: cj_s2r_rep, cj_s2s_rep, cj_r2r_rep, cond_s2r_rep
    real(kind=DP), allocatable :: global_match_weight_sys(:)
    character(len=2), allocatable :: global_label_sys(:), global_label_ref(:)
    character(len=2), allocatable :: old_label_sys(:)
    character(len=5), allocatable :: tmp_string(:)
    character(len=1) :: mark
    logical :: f_exist, diff_flag



    if(pub_qnto_ref_dir /= '') then
    !-------------------------------------------------------------------------!
    ! Construct (multiple) ref_dir from the pub_qnto_ref_dir                  !
    !-------------------------------------------------------------------------!
       string = pub_qnto_ref_dir
       num_dir = 0
       pos = 1
       do while(pos /= 0)
          pos = index(string, ',')
          num_dir = num_dir + 1
          if(pos /= 0) then
             dir_sub = string(1:pos-1)
          else
             dir_sub = string(pos+1:)
          end if
          ! Test if there are subdirectories
          pos_sub = index(dir_sub, '|')
          if(pos_sub /= 0) then
             read(dir_sub(pos_sub-2:pos_sub-1),'(i2)') idx_start
             read(dir_sub(pos_sub+1:pos_sub+2),'(i2)') idx_end
             num_dir = num_dir + abs(idx_end - idx_start)
          end if
          string = string(pos+1:)
       end do

       string = pub_qnto_ref_dir
       allocate(ref_dir(num_dir), stat=ierr)
       call utils_alloc_check('qnto/qnto_calculate', 'ref_dir', ierr)
       it = 1
       do while(it <= num_dir)
          pos = index(string, ',')
          if(pos /= 0) then
             dir_sub = string(1:pos-1)
          else
             dir_sub = string(pos+1:)
          end if
          ! Test if there are subdirectories
          pos_sub = index(dir_sub, '|')
          if(pos_sub /= 0) then
             read(dir_sub(pos_sub-2:pos_sub-1),'(i2)') idx_start
             read(dir_sub(pos_sub+1:pos_sub+2),'(i2)') idx_end
             if(idx_end-idx_start>=0) then
                do sub_it=0, idx_end-idx_start
                   write(ref_dir(it+sub_it),'(a,i2.2)') dir_sub(1:pos_sub-3), &
                        idx_start+sub_it
                   if(pub_on_root) write(stdout,'(/a,a)') 'dir:', ref_dir(it+sub_it)
                end do
             else
                do sub_it=0, idx_end-idx_start, -1
                   write(ref_dir(it+abs(sub_it)),'(a,i2.2)') &
                        dir_sub(1:pos_sub-3), idx_start+sub_it
                   if(pub_on_root) write(stdout,'(/a,a)') 'dir:', ref_dir(it+abs(sub_it))
                end do
             end if
             it = it + abs(idx_end - idx_start) + 1
          else
             ref_dir(it) = dir_sub
             if(pub_on_root) write(stdout,'(/a,a)') 'dir:', ref_dir(it)
             it = it + 1
          end if

          string = string(pos+1:)
       end do
    end if !pub_qnto_ref_dir


    !-------------------------------------------------------------------------!
    ! Construct val_ngwf_1st_idx_inp and val_ngwf_1st_idx_nod for saving the  !
    ! 1st idx of VAL-NGWFs of atoms according to the atom order in input file !
    ! and procs, respectively.                                                !
    ! Similarly, cond_ngwf_1st_idx_inp and cond_ngwf_1st_idx_nod for          !
    ! COND-NGWFs; aux_ngwf_1st_idx_inp and aux_ngwf_1st_idx_nod for AUX-NGWFs.!
    !-------------------------------------------------------------------------!
    allocate(val_ngwf_1st_idx_inp(mdl%nat), stat=ierr)
    call utils_alloc_check('qnto/qnto_calculate', 'val_ngwf_1st_idx_inp', ierr)
    allocate(val_ngwf_1st_idx_nod(mdl%nat), stat=ierr)
    call utils_alloc_check('qnto/qnto_calculate', 'val_ngwf_1st_idx_nod', ierr)
    val_ngwf_1st_idx_inp(1) = 1
    val_ngwf_1st_idx_nod(1) = 1
    do atom_idx=2, mdl%nat
        val_ngwf_1st_idx_inp(atom_idx) = val_ngwf_1st_idx_inp(atom_idx-1) &
                     + mdl%elements(atom_idx-1)%nfunctions
        val_ngwf_1st_idx_nod(atom_idx) = val_ngwf_1st_idx_nod(atom_idx-1) &
                     + mdl%elements(mdl%par%orig_atom(atom_idx-1))%nfunctions
    end do

    if(pub_lr_tddft_joint_set) then
       allocate(joint_ngwf_1st_idx_inp(mdl%nat), stat=ierr)
       call utils_alloc_check('qnto/qnto_calculate', 'joint_ngwf_1st_idx_inp', ierr)
       allocate(joint_ngwf_1st_idx_nod(mdl%nat), stat=ierr)
       call utils_alloc_check('qnto/qnto_calculate', 'joint_ngwf_1st_idx_nod', ierr)
       joint_ngwf_1st_idx_inp(1) = 1
       joint_ngwf_1st_idx_nod(1) = 1
       do atom_idx=2, mdl%nat
          joint_ngwf_1st_idx_inp(atom_idx) = joint_ngwf_1st_idx_inp(atom_idx-1) &
                     + mdl%elements(atom_idx-1)%nfunctions_cond &
                     + mdl%elements(atom_idx-1)%nfunctions
          joint_ngwf_1st_idx_nod(atom_idx) = joint_ngwf_1st_idx_nod(atom_idx-1) &
                     + mdl%elements(mdl%par%orig_atom(atom_idx-1) )%nfunctions_cond&
                     + mdl%elements(mdl%par%orig_atom(atom_idx-1) )%nfunctions
       end do
    else
       allocate(cond_ngwf_1st_idx_inp(mdl%nat), stat=ierr)
       call utils_alloc_check('qnto/qnto_calculate', 'cond_ngwf_1st_idx_inp', ierr)
       allocate(cond_ngwf_1st_idx_nod(mdl%nat), stat=ierr)
       call utils_alloc_check('qnto/qnto_calculate', 'cond_ngwf_1st_idx_nod', ierr)
       cond_ngwf_1st_idx_inp(1) = 1
       cond_ngwf_1st_idx_nod(1) = 1
       do atom_idx=2, mdl%nat
          cond_ngwf_1st_idx_inp(atom_idx) = cond_ngwf_1st_idx_inp(atom_idx-1) &
                     + mdl%elements(atom_idx-1)%nfunctions_cond
          cond_ngwf_1st_idx_nod(atom_idx) = cond_ngwf_1st_idx_nod(atom_idx-1) &
                     + mdl%elements(mdl%par%orig_atom(atom_idx-1) )%nfunctions_cond
       end do
    end if !pub_lr_tddft_joint_set
    if(pub_use_aux_ngwfs) then
       allocate(aux_ngwf_1st_idx_inp(mdl%nat), stat=ierr)
       call utils_alloc_check('qnto/qnto_calculate', 'aux_ngwf_1st_idx_inp', ierr)
       allocate(aux_ngwf_1st_idx_nod(mdl%nat), stat=ierr)
       call utils_alloc_check('qnto/qnto_calculate', 'aux_ngwf_1st_idx_nod', ierr)
       aux_ngwf_1st_idx_inp(1) = 1
       aux_ngwf_1st_idx_nod(1) = 1
       do atom_idx=2, mdl%nat
           aux_ngwf_1st_idx_inp(atom_idx) = aux_ngwf_1st_idx_inp(atom_idx-1) &
                     + mdl%elements(atom_idx-1)%nfunctions_aux
           aux_ngwf_1st_idx_nod(atom_idx) = aux_ngwf_1st_idx_nod(atom_idx-1) &
                     + mdl%elements(mdl%par%orig_atom(atom_idx-1) )%nfunctions_aux
       end do
    end if !pub_use_aux_ngwfs


    !-------------------------------------------------------------------------!
    ! Read in NBOs from previous NBO generation cal.  If AUX-NGWFs are used,  !
    ! also read in STD-NBOs.  If pub_qnto_ref_dir == '', output NGWF-coeff of !
    ! (STD-)NBOs to *_NBO.dat(_std).  If pub_qnto_ref_dir /= '', Sys-NBOs are !
    ! projected to Ref-NBOs and saved in *_NBO_s2r_*.  '*' denotes any string.!
    !-------------------------------------------------------------------------!
    if(pub_qnto_nbo_proj) call internal_readin_nbo


    !-------------------------------------------------------------------------!
    ! Set up the number of NTO transitions.                                   !
    !-------------------------------------------------------------------------!
    num_nto_trans = pub_qnto_num_transition
    if(val_rep%n_occ(1,PUB_1K)>5) then
       num_diag = 5
    else
       num_diag = val_rep%n_occ(1,PUB_1K)
    end if
    if(pub_qnto_ref_dir == '') then
       if(num_nto_trans /= num_diag) then
          num_nto_trans = num_diag
          if(pub_on_root) &
             write(stdout,'(/a,i1)') 'Number of NTO transitions is reset to ', num_diag
       end if
    else !pub_qnto_ref_dir /= ''
       if(num_nto_trans > num_diag) then
          num_nto_trans = num_diag
          if(pub_on_root) &
             write(stdout,'(/a,i1)') 'Number of NTO transitions is reset to ', num_diag
       end if
    end if


    if(pub_qnto_ref_dir /= '') then
    !-------------------------------------------------------------------------!
    ! Construct matrices for saving NGWF-coeff of NTOn-H(E) of the Sys as well!
    ! as matrices saving those of the Ref NTOn-H(E) of the Sys will later be  !
    ! projected to those of the Ref                                           !
    !-------------------------------------------------------------------------!
       systgt_size = num_nto_trans * num_states
       reftgt_size = num_nto_trans * pub_qnto_num_ref_states
       call dense_create(ntoh_systgt_dens, val_ngwf_basis%num, systgt_size)
       call dense_create(ntoh_sctgt_dens, val_ngwf_basis%num, systgt_size)
       if(pub_lr_tddft_joint_set) then
          call dense_create(ntoe_systgt_dens, joint_ngwf_basis%num, systgt_size)
          call dense_create(ntoe_sctgt_dens, joint_ngwf_basis%num, systgt_size)
       else
          call dense_create(ntoe_systgt_dens, cond_ngwf_basis%num, systgt_size)
          call dense_create(ntoe_sctgt_dens, cond_ngwf_basis%num, systgt_size)
       end if
    end if !pub_qnto_ref_dir


    !-------------------------------------------------------------------------!
    ! Perform QNTO analysis for each excitation I                             !
    !-------------------------------------------------------------------------!
    if(pub_on_root) then
       qnto_output_file = trim(pub_rootname)//'_QNTO_info.dat'
       qnto_output_file = adjustl(qnto_output_file)
       qnto_output_unit = utils_unit()

       open(unit=qnto_output_unit,form="formatted", &
            file=trim(qnto_output_file),action="write",iostat=ierr)
       call utils_open_unit_check('qnto/qnto_calculate',trim(qnto_output_file),ierr)
    end if
    if(val_rep%n_occ(1,PUB_1K)>5) then
       num_diag = 5
    else
       num_diag = val_rep%n_occ(1,PUB_1K)
    end if
    do I=1, num_states
       if(pub_qnto_ref_dir /= '') then
          call mo_coeff_i(I, num_nto_trans, val_ngwf_basis, mdl, &
               val_ngwf_1st_idx_inp, val_ngwf_1st_idx_nod, 'sys-val', '', &
               ntoh_systgt_dens, ntoh_sctgt_dens)
          if(pub_lr_tddft_joint_set) then
             call mo_coeff_i(I, num_nto_trans, joint_ngwf_basis, mdl, &
                  joint_ngwf_1st_idx_inp, joint_ngwf_1st_idx_nod, 'sys-joint', '', &
                  ntoe_systgt_dens, ntoe_sctgt_dens)
          else
             call mo_coeff_i(I, num_nto_trans, cond_ngwf_basis, mdl, &
                  cond_ngwf_1st_idx_inp, cond_ngwf_1st_idx_nod, 'sys-cond', '', &
                  ntoe_systgt_dens, ntoe_sctgt_dens)
          end if
       else
          !-------------------------------------------------------------------------!
          ! Solve for the valence and conduction NTO bases                          !
          !-------------------------------------------------------------------------!
          if(pub_qnto_svd_method == 1) then
             ! trans_coeff, ntoh_sys_dens, and ntoe_sys_dens are constructed
             if(pub_lr_tddft_joint_set) then
                call solve_nto_svd1(num_nto_trans, num_min_nto_trans, trans_coeff, &
                     tot_perc, ntoh_sys_dens, ntoe_sys_dens, response_kernel(I), &
                     val_rep, val_ngwf_basis, joint_rep, joint_ngwf_basis)
             else
                call solve_nto_svd1(num_nto_trans, num_min_nto_trans, trans_coeff, &
                     tot_perc, ntoh_sys_dens, ntoe_sys_dens, response_kernel(I), &
                     val_rep, val_ngwf_basis, cond_rep, cond_ngwf_basis)
             end if
          else if(pub_qnto_svd_method == 2) then
             if(pub_lr_tddft_joint_set) then
                call solve_nto_svd2(num_nto_trans, num_min_nto_trans, trans_coeff, &
                     tot_perc, ntoh_sys_dens, ntoe_sys_dens, response_kernel(I), &
                     val_rep, val_ngwf_basis, joint_rep, joint_ngwf_basis)
             else
                call solve_nto_svd2(num_nto_trans, num_min_nto_trans, trans_coeff, &
                     tot_perc, ntoh_sys_dens, ntoe_sys_dens, response_kernel(I), &
                     val_rep, val_ngwf_basis, cond_rep, cond_ngwf_basis)
             end if
          else
             call utils_abort('ERROR in qnto_calculate: qnto_svd_method has to &
                  &be set to 1 or 2')
          end if !pub_qnto_svd_method

          !-------------------------------------------------------------------------!
          ! Project NTOn-H and NTOn-E to NBOs                                       !
          !-------------------------------------------------------------------------!
          if(pub_qnto_nbo_proj) call internal_projto_nbo

          !-------------------------------------------------------------------------!
          ! Output NTO transition vector components to stdout and QNTO_info file    !
          !-------------------------------------------------------------------------!
          if(pub_on_root) then
             write(stdout, '(a,i4,a)') '===================  Excitation: ', &
                  I, '  ==================='
             write(qnto_output_unit, '(a,i4,a)') '===================  Excitation: ',&
                  I, '  ==================='
             write(qnto_output_unit, '(a,f12.8)') 'Eigenvalue (in eV) : ', &
                  exc_energies(I)*HARTREE_IN_EVS
             do n=1, num_nto_trans
                write(stdout, '(1x,a,i1,a,f7.4)') 'NTO', n, &
                     ' Trans. Coe. = ', trans_coeff(n)
                write(qnto_output_unit, '(1x,a,i1,a,f7.4)') 'NTO', n, &
                     ' Trans. Coe. = ', trans_coeff(n)
             end do
             write(stdout, '(1x,a,f5.2,a)') 'Total Trans. Dens. &
                  &accounted for: ', tot_perc*100.0_DP, '%'
             if(num_min_nto_trans == num_diag+1) then
                num_min_nto_trans = num_diag
                write(stdout,'(a,i1)') ' No. of NTO pairs to &
                     &reach 95% Trans. Dens.: >', num_min_nto_trans
                write(qnto_output_unit,'(a,i1)') ' No. of NTO pairs to &
                     &reach 95% Trans. Dens.: >', num_min_nto_trans
             else
                write(stdout,'(a,i3)') ' No. of NTO pairs to &
                     &reach 95% Trans. Dens.:', num_min_nto_trans
                write(qnto_output_unit,'(a,i3)') ' No. of NTO pairs to &
                     &reach 95% Trans. Dens.:', num_min_nto_trans
             end if
             if( num_nto_trans /= 1 ) then
                write(stdout,'(a,i3,a,i3,a)') ' Outputting',&
                                              num_nto_trans, ' NTO Trans. ...'
             else
                write(stdout,'(a,i3,a,i3,a)') ' Outputting',&
                                              num_nto_trans, ' NTO Trans. ...'
             end if
             write(qnto_output_unit, '(1x,a,f5.2,a)') 'Total Trans. Dens. &
                  &accounted for: ', tot_perc*100.0_DP, '%'
          end if !pub_on_root
          deallocate(trans_coeff, stat=ierr)
          call utils_dealloc_check('qnto/qnto_calculate', 'trans_coeff', ierr)

          if(pub_qnto_nbo_proj) then
          !-------------------------------------------------------------------------!
          ! Output NTO-to-NBO projection results to stdout and QNTO_info file       !
          !-------------------------------------------------------------------------!
             do n=1, num_nto_trans
                call output_nbo(ntoh2nbo_dens, nbo_list, num_nbo, n, qnto_output_unit, 'H')
                call output_nbo(ntoe2nbo_dens, nbo_list, num_nbo, n, qnto_output_unit, 'E')
             end do !n
             call dense_destroy(ntoh2nbo_dens)
             call dense_destroy(ntoe2nbo_dens)
          end if !pub_qnto_nbo_proj

          !-------------------------------------------------------------------------!
          ! Output NGWF-coeff of NTOn-H(E) to *_NTOn-H(E)_#.dat   #: integer        !
          ! and save NTOn-H(E) in *_NTO#-H(E)_#.cube                                !
          !-------------------------------------------------------------------------!
          call mo_coeff_o(I, num_nto_trans, val_ngwf_basis, ntoh_sys_dens,&
               val_rep, mdl, val_ngwf_1st_idx_inp, val_ngwf_1st_idx_nod, 'val', val_denskern)
          if(pub_lr_tddft_joint_set) then
             call mo_coeff_o(I, num_nto_trans, joint_ngwf_basis, &
                  ntoe_sys_dens, joint_rep, mdl, joint_ngwf_1st_idx_inp, &
                  joint_ngwf_1st_idx_nod, 'joint', joint_denskern)
          else
             call mo_coeff_o(I, num_nto_trans, cond_ngwf_basis, &
                  ntoe_sys_dens, cond_rep, mdl, cond_ngwf_1st_idx_inp, &
                  cond_ngwf_1st_idx_nod, 'cond', cond_denskern)
          end if
          if(pub_on_root) then
             write(stdout,*)
             write(qnto_output_unit,*)
          end if

          call dense_destroy(ntoh_sys_dens)
          call dense_destroy(ntoe_sys_dens)
      end if !pub_qnto_ref_dir
    end do !I
    if(pub_on_root) then
       close(unit=qnto_output_unit,iostat=ierr)
       call utils_close_unit_check('qnto/qnto_calculate', 'qnto_output_unit',ierr)
    end if


    if(pub_qnto_ref_dir /= '') then
    !-------------------------------------------------------------------------!
    ! Read in NGWF-coeff of NTOn-H and NTOn-E of Ref                          !
    !-------------------------------------------------------------------------!
       ! Construct needed matrices and arrays
       call internal_ini_clo_ref_proj(1)

       allocate(global_label_sys(num_states),stat=ierr)
       call utils_alloc_check('qnto/qnto_calculate','global_label_sys', ierr)
       allocate(old_label_sys(num_states),stat=ierr)
       call utils_alloc_check('qnto/qnto_calculate','old_label_sys', ierr)
       allocate(global_label_ref(num_states),stat=ierr)
       call utils_alloc_check('qnto/qnto_calculate','global_label_ref', ierr)
       allocate(global_match_weight_sys(num_states),stat=ierr)
       call utils_alloc_check('qnto/qnto_calculate','global_match_weight_sys', ierr)
       allocate(tmp_string(num_states), stat=ierr)
       call utils_alloc_check('qnto/qnto_calculate','tmp_string', ierr)
       ! Set up global_match_weight_sys(:), global_label_sys(:), global_label_ref(:)
       if(pub_on_root) then
          write(label_io_filename,'(a)') trim(pub_rootname)//'_state_labellings.dat.tmp'
          label_io_filename = adjustl(label_io_filename)
          label_io_unit = utils_unit()
          inquire(FILE=trim(label_io_filename), EXIST=f_exist)
          if(.not. f_exist) then
             do state_idx_sys=1, num_states
                global_match_weight_sys(state_idx_sys) = 0.0_DP
                global_label_sys(state_idx_sys) = '  '
                global_label_ref(state_idx_sys) = '  '
             end do
          else
             open(unit=label_io_unit,form="formatted", file=trim(label_io_filename),&
                  action="read", iostat=ierr)
             call utils_open_unit_check('qnto/qnto_calculate',&
                  trim(label_io_filename),ierr)
             n = 0
             do
                read(label_io_unit,'(a)',advance='no') mark
                if(mark == '-') then
                   rewind(label_io_unit)
                   do i=1, n-3
                      read(label_io_unit,*)
                   end do
                   read(label_io_unit,'(3x)',advance='no')
                   do state_idx_sys=1, num_states
                      read(label_io_unit,'(f6.5)',advance='no') &
                          global_match_weight_sys(state_idx_sys)
                   end do
                   read(label_io_unit,*)
                   read(label_io_unit,'(4x)',advance='no')
                   do state_idx_sys=1, num_states
                      read(label_io_unit,'(a,4x)',advance='no') &
                          global_label_ref(state_idx_sys)
                   end do
                   read(label_io_unit,*)
                   read(label_io_unit,*)
                   read(label_io_unit,'(4x)',advance='no')
                   do state_idx_sys=1, num_states
                      read(label_io_unit,'(a,4x)',advance='no') &
                          global_label_sys(state_idx_sys)
                   end do
                   exit
                else
                   read(label_io_unit,*)
                   n = n + 1
                end if
             end do
          end if
          close(unit=label_io_unit,iostat=ierr)
          call utils_close_unit_check('qnto/qnto_calculate',&
               trim(label_io_filename),ierr)
       end if !pub_on_root

       do it=1, num_dir
          do I=1, pub_qnto_num_ref_states
             call mo_coeff_i(I, num_nto_trans, val_ngwf_basis, mdl, &
                  val_ngwf_1st_idx_inp, val_ngwf_1st_idx_nod, 'val', ref_dir(it), &
                  ntoh_reftgt_dens, ntoh_rctgt_dens)

             if(pub_lr_tddft_joint_set) then
                call mo_coeff_i(I, num_nto_trans, joint_ngwf_basis, &
                     mdl, joint_ngwf_1st_idx_inp, joint_ngwf_1st_idx_nod, &
                     'joint', ref_dir(it), ntoe_reftgt_dens, ntoe_rctgt_dens)
             else
                call mo_coeff_i(I, num_nto_trans, cond_ngwf_basis, &
                     mdl, cond_ngwf_1st_idx_inp, cond_ngwf_1st_idx_nod, &
                     'cond', ref_dir(it), ntoe_reftgt_dens, ntoe_rctgt_dens)
             end if
          end do !I

          !-------------------------------------------------------------------------!
          ! Perform Sys-to-Ref projection for NTOn-H and NTOn-E                     !
          !-------------------------------------------------------------------------!
          if(pub_lr_tddft_joint_set) then
             call nto_sys_ref_proj( &
                  nto_sch2sch, nto_sce2sce, nto_sch2sce, nto_sce2sch, &
                  nto_sh2sch, nto_se2sce, nto_sh2sce, nto_se2sch, &
                  nto_rch2rch, nto_rce2rce, nto_rch2rce, nto_rce2rch, &
                  nto_sch2rch, nto_sce2rce, nto_sch2rce, nto_sce2rch, &
                  ntoh_rctgt_dens, ntoh_sctgt_dens, ntoe_rctgt_dens, &
                  ntoe_sctgt_dens, ntoh_reftgt_dens, ntoh_systgt_dens, &
                  ntoe_reftgt_dens, ntoe_systgt_dens, val_ngwf_basis, &
                  joint_ngwf_basis, val_rep, joint_rep, systgt_size, &
                  reftgt_size, mdl, ref_dir(it), &
                  val_s2r_rep, cj_s2r_rep, val_r2r_rep, cj_r2r_rep, &
                  val_s2s_rep, cj_s2s_rep, cond_s2r_rep, cond_ngwf_basis)
          else
             call nto_sys_ref_proj( &
                  nto_sch2sch, nto_sce2sce, nto_sch2sce, nto_sce2sch, &
                  nto_sh2sch, nto_se2sce, nto_sh2sce, nto_se2sch, &
                  nto_rch2rch, nto_rce2rce, nto_rch2rce, nto_rce2rch, &
                  nto_sch2rch, nto_sce2rce, nto_sch2rce, nto_sce2rch, &
                  ntoh_rctgt_dens, ntoh_sctgt_dens, ntoe_rctgt_dens, &
                  ntoe_sctgt_dens, ntoh_reftgt_dens, ntoh_systgt_dens, &
                  ntoe_reftgt_dens, ntoe_systgt_dens, val_ngwf_basis, &
                  cond_ngwf_basis, val_rep, cond_rep, systgt_size, &
                  reftgt_size, mdl, ref_dir(it), &
                  val_s2r_rep, cj_s2r_rep, val_r2r_rep, cj_r2r_rep, &
                  val_s2s_rep, cj_s2s_rep, cond_s2r_rep)
          end if

          !-------------------------------------------------------------------------!
          ! Output Sys-to-Ref projection results of NTOn-H and NTOn-E to files:     !
          ! *_sc2sc.dat, *_s2sc(_HE).dat, *_rc2rc_*.dat, *_r2rc_*.dat, and          !
          ! *_sc2rc_*(HE).dat                                                       !
          !-------------------------------------------------------------------------!
          call nto_output_proj_coeff( &
               nto_sch2sch, nto_sce2sch, nto_sch2sce, nto_sce2sce, &
               nto_sh2sch, nto_se2sch, nto_sh2sce, nto_se2sce, &
               nto_rch2rch, nto_rce2rch, nto_rch2rce, nto_rce2rce, &
               nto_sch2rch, nto_sce2rch, nto_sch2rce, nto_sce2rce, &
               num_states, systgt_size, reftgt_size, ref_dir(it), &
               global_label_sys, global_label_ref, global_match_weight_sys, &
               val_rep%n_occ(1,PUB_1K), num_nto_trans)
       end do !it

       !-------------------------------------------------------------------------!
       ! Read in labels of states if *_state_labellings.dat exists               !
       ! Compare these labels with newly generated lables and write both of them !
       ! in *_state_labellings.dat.tmp for checking self-consistency             !
       !-------------------------------------------------------------------------!
       if(pub_on_root) then
          write(label_io_filename,'(a)') trim(pub_rootname)//'_state_labellings.dat'
          label_io_filename = adjustl(label_io_filename)
          label_io_unit = utils_unit()
          inquire(FILE=trim(label_io_filename), EXIST=f_exist)
          if(.not. f_exist) then
             ! Set up old_label_sys(:)
             do state_idx_sys=1, num_states
                old_label_sys(state_idx_sys) = '  '
             end do
          else
             ! Read in old_label_sys(:)
             open(unit=label_io_unit,form="formatted", file=trim(label_io_filename),&
                  action="read", iostat=ierr)
             call utils_open_unit_check('qnto/nto_output_proj_coeff',&
                  trim(label_io_filename),ierr)
             do
                read(label_io_unit,'(a)',advance='no') mark
                if(mark == '-') then
                   read(label_io_unit,'(3x)',advance='no')
                   do state_idx_sys=1, num_states
                      read(label_io_unit,'(a,4x)',advance='no') &
                          old_label_sys(state_idx_sys)
                   end do
                   exit
                else
                   read(label_io_unit,*)
                end if
             end do
             close(unit=label_io_unit,iostat=ierr)
             call utils_close_unit_check('qnto/nto_output_proj_coeff',&
                  trim(label_io_filename),ierr)
          end if !f_exist

          ! Write old_label_sys(:) and global_label_sys(:)
          write(label_io_filename,'(a)') trim(pub_rootname)//'_state_labellings.dat.tmp'
          label_io_filename = adjustl(label_io_filename)
          label_io_unit = utils_unit()
          open(unit=label_io_unit,form="formatted", file=trim(label_io_filename),&
               action="readwrite", iostat=ierr)
          call utils_open_unit_check('qnto/nto_output_proj_coeff',&
               trim(label_io_filename),ierr)

          n = 0
          do
             read(label_io_unit,'(a)',advance='no') mark
             if(mark == '-') then
                read(label_io_unit,'(3x)',advance='no')
                do state_idx_sys=1, num_states
                   read(label_io_unit,'(a,1x)',advance='no') tmp_string(state_idx_sys)
                end do
                rewind(label_io_unit)
                do i=1, n
                   read(label_io_unit,*)
                end do
                write(label_io_unit,'(a,2x,a)',advance='no') '%', '|'
                do state_idx_sys=1, num_states
                   write(label_io_unit,'(a,a)',advance='no') tmp_string(state_idx_sys), '|'
                end do
                write(label_io_unit,*)
                write(label_io_unit,'(a)',advance='no') '%'
                exit
             else
                read(label_io_unit,*)
                n = n + 1
             end if
          end do

          write(label_io_unit,*)
          write(label_io_unit,'(a)',advance='no') '%  '
          do state_idx_sys=1, num_states
             write(label_io_unit,'(f6.5)',advance='no') global_match_weight_sys(state_idx_sys)
          end do
          write(label_io_unit,'(a)') '|Global Match Weight'
          write(label_io_unit,'(a)',advance='no') '%  |'
          do state_idx_sys=1, num_states
             write(label_io_unit,'(a,a)',advance='no') global_label_ref(state_idx_sys), '   |'
          end do
          write(label_io_unit,'(a)') 'Global Ref label'
          write(label_io_unit,'(a)',advance='no') '%  |'
          do state_idx_sys=1, num_states
             write(label_io_unit,'(a,a)',advance='no') old_label_sys(state_idx_sys), '   |'
          end do
          write(label_io_unit,'(a)') 'Old Sys label'
          write(label_io_unit,'(a)',advance='no') '-> |'
          diff_flag = .FALSE.
          do state_idx_sys=1, num_states
             if( global_match_weight_sys(state_idx_sys) > sqrt(1.0_DP/3.0_DP) &
                 .OR. global_label_sys(state_idx_sys) == '--' ) then
                write(label_io_unit,'(a,a)',advance='no') &
                     global_label_sys(state_idx_sys), '   |'
                if(old_label_sys(state_idx_sys) /= global_label_sys(state_idx_sys)) &
                  diff_flag = .TRUE.
             else
                write(label_io_unit,'(a)',advance='no') '     |'
                if(old_label_sys(state_idx_sys) /= '  ') diff_flag = .TRUE.
             end if
          end do
          if (diff_flag .eqv. .TRUE.) then
             write(label_io_unit,'(a)') 'New Sys label(Diff)'
          else
             write(label_io_unit,'(a)') 'New Sys label'
          end if

          close(unit=label_io_unit,iostat=ierr)
          call utils_close_unit_check('qnto/nto_output_proj_coeff',&
               trim(label_io_filename),ierr)
       end if !pub_on_root
       deallocate(global_label_sys,stat=ierr)
       call utils_dealloc_check('qnto/qnto_calculate','global_label_sys', ierr)
       deallocate(old_label_sys,stat=ierr)
       call utils_dealloc_check('qnto/qnto_calculate','old_label_sys', ierr)
       deallocate(global_label_ref,stat=ierr)
       call utils_dealloc_check('qnto/qnto_calculate','global_label_ref', ierr)
       deallocate(global_match_weight_sys,stat=ierr)
       call utils_dealloc_check('qnto/qnto_calculate','global_match_weight_sys', ierr)
       deallocate(tmp_string,stat=ierr)
       call utils_dealloc_check('qnto/qnto_calculate','tmp_string', ierr)

       ! Destruct matrices and arrays
       call internal_ini_clo_ref_proj(0)
    end if !pub_qnto_ref_dir


    !--------------------------------------------------------------------!
    ! Deallocate remaining matrices and arrays used                      !
    !--------------------------------------------------------------------!
    deallocate(val_ngwf_1st_idx_inp, stat=ierr)
    call utils_dealloc_check('qnto/qnto_calculate', 'val_ngwf_1st_idx_inp', ierr)
    deallocate(val_ngwf_1st_idx_nod, stat=ierr)
    call utils_dealloc_check('qnto/qnto_calculate', 'val_ngwf_1st_idx_nod', ierr)
    if(pub_lr_tddft_joint_set) then
       deallocate(joint_ngwf_1st_idx_inp, stat=ierr)
       call utils_dealloc_check('qnto/qnto_calculate', 'joint_ngwf_1st_idx_inp', ierr)
       deallocate(joint_ngwf_1st_idx_nod, stat=ierr)
       call utils_dealloc_check('qnto/qnto_calculate', 'joint_ngwf_1st_idx_inp', ierr)
    else
       deallocate(cond_ngwf_1st_idx_inp, stat=ierr)
       call utils_dealloc_check('qnto/qnto_calculate', 'cond_ngwf_1st_idx_inp', ierr)
       deallocate(cond_ngwf_1st_idx_nod, stat=ierr)
       call utils_dealloc_check('qnto/qnto_calculate', 'cond_ngwf_1st_idx_inp', ierr)
    end if
    if(pub_use_aux_ngwfs) then
       deallocate(aux_ngwf_1st_idx_inp, stat=ierr)
       call utils_dealloc_check('qnto/qnto_calculate', 'aux_ngwf_1st_idx_inp', ierr)
       deallocate(aux_ngwf_1st_idx_nod, stat=ierr)
       call utils_dealloc_check('qnto/qnto_calculate', 'aux_ngwf_1st_idx_nod', ierr)
    end if !pub_use_aux_ngwfs

    if(pub_qnto_nbo_proj) then
       deallocate(nbo_sys_dens, stat=ierr)
       call utils_dealloc_check('qnto/qnto_calculate', 'nbo_sys_dens', ierr)
       if(pub_use_aux_ngwfs) then
          deallocate(stdnbo_sys_dens, stat=ierr)
          call utils_dealloc_check('qnto/qnto_calculate', 'stdnbo_sys_dens', ierr)
          deallocate(stdnbo_list, stat=ierr)
          call utils_dealloc_check('qnto/qnto_calculate', 'stdnbo_list', ierr)
       end if ! pub_use_aux_ngwfs
       deallocate(nbo_list, stat=ierr)
       call utils_dealloc_check('qnto/qnto_calculate', 'nbo_list', ierr)
    end if !pub_qnto_nbo_proj


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

  subroutine internal_readin_nbo

    implicit none

    !--------------------------------------------------------------------!
    ! Read in a natural bond orbitals (NBOs) set (S1) for interpreting   !
    ! the NTOn-H and NTOn-E in QNTO analysis.  S1 may be generated from  !
    ! the auxiliary NGWFs set in a PROPERTIES calculation.               !
    ! "xx_inittr_nao_nbo.dat" and "xx_nao.37" are needed,                !
    ! the former storing the expression of NAOs in NGWFs, whereas        !
    ! the latter storing the expression of NBOs in NAOs.                 !
    ! Reference: properties_mod.F90/internal_plot_nbo                    !
    !--------------------------------------------------------------------!
    pre_ao = .false.
    nbo_trfile = '.37'
    nbo_trfile = adjustl(nbo_trfile)
    write(nbo_trfile,'(a)') trim(pub_rootname)//'_nao'//trim(nbo_trfile)
    nbo_trfile = adjustl(nbo_trfile)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Inquires if input file from nbo_trfile is spin-polarized
    if(pub_use_aux_ngwfs) then
       call npa_read_nbodata(nbo_sys_dens, aux_ngwf_basis, trim(nbo_trfile),&
            pre_ao, num_ao, spin_inquire=.true., nspin=nspin, &
            par=mdl%par)
    else
       call npa_read_nbodata(nbo_sys_dens, val_ngwf_basis, trim(nbo_trfile),&
            pre_ao, num_ao, spin_inquire=.true., nspin=nspin, &
            par=mdl%par)
    end if !pub_use_aux_ngwfs
    allocate(nbo_sys_dens(nspin), stat=ierr) ! nspin = 1 or 2
    call utils_alloc_check('qnto/internal_readin_nbo', 'nbo_sys_dens', ierr)
    do is=1, nspin
       if(pub_use_aux_ngwfs) then
          call dense_create(nbo_sys_dens(is), aux_ngwf_basis%num, &
               aux_ngwf_basis%num)
       else
          call dense_create(nbo_sys_dens(is), val_ngwf_basis%num, &
               val_ngwf_basis%num)
       endif !pub_use_aux_ngwfs
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Retrieve expressions of NBOs in NGWFs, stored in nbo_sys_dens
    if(pub_use_aux_ngwfs) then
       call npa_read_nbodata(nbo_sys_dens, aux_ngwf_basis, trim(nbo_trfile),&
            pre_ao, num_ao, spin_inquire=.false., nspin=nspin, &
            par=mdl%par)
    else
       call npa_read_nbodata(nbo_sys_dens, val_ngwf_basis, trim(nbo_trfile),&
            pre_ao, num_ao, spin_inquire=.false., nspin=nspin, &
            par=mdl%par)
    end if !pub_use_aux_ngwfs
    ! Set up num_nbo, num_stdnbo
    if(pub_use_aux_ngwfs) then
       num_nbo = aux_ngwf_basis%num
    else
       num_nbo = val_ngwf_basis%num
    endif !pub_use_aux_ngwfs
    num_stdnbo = val_ngwf_basis%num


    if(pub_use_aux_ngwfs) then ! S2 /= S1
       !--------------------------------------------------------------------!
       ! Read in a standard NBOs set (S2) for calibrating the sign and order!
       ! of the orbitals in S1.  S2 may be generated from the optimised     !
       ! NGWFs in a valence calculation.                                    !
       ! "xx_inittr_nao_nbo.dat_std" and "xx_nao.37_std" are needed.        !
       !--------------------------------------------------------------------!
       write(nbo_trfile,'(a)') trim(nbo_trfile)//'_std'
       nbo_trfile = adjustl(nbo_trfile)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Inquires if input file from nbo_trfile is spin-polarized
       call npa_read_nbodata(stdnbo_sys_dens, val_ngwf_basis, trim(nbo_trfile),&
            pre_ao, num_ao, spin_inquire=.true., nspin=nspin, &
            par=mdl%par)
       allocate(stdnbo_sys_dens(nspin), stat=ierr) ! nspin = 1 or 2
       call utils_alloc_check('qnto/internal_readin_nbo', 'stdnbo_sys_dens',ierr)
       do is=1, nspin
          call dense_create(stdnbo_sys_dens(is), val_ngwf_basis%num, &
               val_ngwf_basis%num)
       end do
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Retrieve expressions of NBOs in NGWFs, stored in stdnbo_sys_dens
       call npa_read_nbodata(stdnbo_sys_dens, val_ngwf_basis, trim(nbo_trfile),&
            pre_ao, num_ao, spin_inquire=.false., nspin=nspin, &
            par=mdl%par, std=.true.)
     end if ! pub_use_aux_ngwfs


    !--------------------------------------------------------------------!
    ! Read in NBO info from "outNBO" and "outNBO_std", stored in         !
    ! nbo_list and stdnbo_list, respectively.  Calibrate the S1 NBOs set !
    ! (re-ordering and re-signing orbitals) according to set S2.         !
    !--------------------------------------------------------------------!
    apx = ''
    call read_nbo_info(nbo_list, apx, num_nbo)
    if(pub_use_aux_ngwfs) then ! S2 /= S1
       apx = '_std'
       call read_nbo_info(stdnbo_list, apx, num_stdnbo)
          call calibrate_nbo(nbo_sys_dens(1), stdnbo_sys_dens(1), &
               nbo_list, stdnbo_list, aux_rep%cross_overlap, &
               num_stdnbo, num_nbo)
    end if ! pub_use_aux_ngwfs


    if(pub_use_aux_ngwfs) then
       !--------------------------------------------------------------------!
       ! Output NGWF-coeff of NBOs and STD-NBOs.  If pub_qnto_ref_dir /= '',!
       ! read in corresponding NGWF-coeff of NBOs of the Ref and perform    !
       ! Sys-to-Ref projection                                              !
       !--------------------------------------------------------------------!
       if(pub_qnto_ref_dir == '') then
          call mo_coeff_o(1, num_nbo, aux_ngwf_basis, nbo_sys_dens(1), &
               aux_rep, mdl, aux_ngwf_1st_idx_inp, aux_ngwf_1st_idx_nod, &
               'nbo', aux_denskern)
          call mo_coeff_o(1, num_stdnbo, val_ngwf_basis, stdnbo_sys_dens(1), &
               val_rep, mdl, val_ngwf_1st_idx_inp, val_ngwf_1st_idx_nod, &
               'std-nbo', val_denskern)
       end if ! pub_qnto_ref_dir
    else
       if(pub_qnto_ref_dir == '') then
          call mo_coeff_o(1, num_nbo, val_ngwf_basis, nbo_sys_dens(1), val_rep, &
               mdl, val_ngwf_1st_idx_inp, val_ngwf_1st_idx_nod, 'nbo', val_denskern)
       end if ! pub_qnto_ref_dir
    end if !pub_use_aux_ngwfs
  end subroutine internal_readin_nbo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine internal_projto_nbo

    implicit none

    !--------------------------------------------------------------------!
    ! Project the NTOn-H and NTOn-E to the NBOs set read in              !
    !--------------------------------------------------------------------!
    if(pub_use_aux_ngwfs) then
       ! Construct expressions of NTOn-H in NBOs (ntoh2nbo_dens)
       !==================== Tv = <AUX-NGWFs|VAL-NGWFs> ====================!
       call sparse_embed_transpose_structure(Tv%structure, aux_rep%cross_overlap) !
       call sparse_embed_create(Tv)                                               !
       call sparse_embed_transpose(Tv, aux_rep%cross_overlap)                     !
       call dense_create(Tv_dens, aux_ngwf_basis%num, val_ngwf_basis%num)         !
       call dense_convert(Tv_dens, Tv)                                            !
       !====================================================================!
       call dense_create(Tvntoh_sys_dens, aux_ngwf_basis%num, num_nto_trans)
       call dense_create(ntoh2nbo_dens, aux_ngwf_basis%num, num_nto_trans)
       call dense_product(Tvntoh_sys_dens, Tv_dens, ntoh_sys_dens)
       call dense_product(ntoh2nbo_dens, nbo_sys_dens(1), Tvntoh_sys_dens, &
            opA='T')
       call dense_destroy(Tvntoh_sys_dens)
       call dense_destroy(Tv_dens)
       ! Construct expressions of NTOn-E in NBOs (ntoe2nbo_dens)
       !==================== Tc = <AUX-NGWFs|COND-NGWFs> ===================!
       call sparse_embed_transpose_structure( &                             !
            response_kernel_trans%structure, response_kernel(I))            !
       call sparse_embed_create(response_kernel_trans)                      !
       call sparse_embed_create(Tc, Tv, response_kernel_trans)              !
       call sparse_embed_destroy(response_kernel_trans)                     !
       call sparse_embed_destroy(Tv)                                        !
       if(pub_lr_tddft_joint_set) then                                      !
          call function_ops_brappd_ketppd(Tc%p, &                           !
               aux_rep%ngwfs_on_grid(1), aux_ngwf_basis, &                  !
               joint_rep%ngwfs_on_grid(1), joint_ngwf_basis, mdl%cell)      !
       else                                                                 !
          call function_ops_brappd_ketppd(Tc%p, &                           !
               aux_rep%ngwfs_on_grid(1), aux_ngwf_basis, &                  !
               cond_rep%ngwfs_on_grid(1), cond_ngwf_basis, mdl%cell)        !
       end if                                                               !
       if(pub_lr_tddft_joint_set) then                                      !
          call dense_create(Tc_dens,aux_ngwf_basis%num,joint_ngwf_basis%num)!
       else                                                                 !
          call dense_create(Tc_dens,aux_ngwf_basis%num,cond_ngwf_basis%num) !
       end if                                                               !
       call dense_convert(Tc_dens, Tc)                                      !
       call sparse_embed_destroy(Tc)                                        !
       !====================================================================!
       if(pub_lr_tddft_joint_set) then
          call dense_create(Tcntoe_sys_dens, aux_ngwf_basis%num, num_nto_trans)
          call dense_create(ntoe2nbo_dens, aux_ngwf_basis%num, num_nto_trans)
       else
          call dense_create(Tcntoe_sys_dens, aux_ngwf_basis%num, num_nto_trans)
          call dense_create(ntoe2nbo_dens, aux_ngwf_basis%num, num_nto_trans)
       end if
       call dense_product(Tcntoe_sys_dens, Tc_dens, ntoe_sys_dens)
       call dense_product(ntoe2nbo_dens, nbo_sys_dens(1), &
                          Tcntoe_sys_dens, opA='T')
       call dense_destroy(Tcntoe_sys_dens)
       call dense_destroy(Tc_dens)
    else
       ! Construct expressions of NTOn-H in NBOs (ntoh2nbo_dens)
       !==================== Tv = <VAL-NGWFs|VAL-NGWFs> ====================!
       call dense_create(Tv_dens, val_ngwf_basis%num, val_ngwf_basis%num)   !
       call dense_convert(Tv_dens, val_rep%overlap)                         !
       !====================================================================!
       call dense_create(Tvntoh_sys_dens, val_ngwf_basis%num, num_nto_trans)
       call dense_create(ntoh2nbo_dens, val_ngwf_basis%num, num_nto_trans)
       call dense_product(Tvntoh_sys_dens, Tv_dens, ntoh_sys_dens)
       call dense_product(ntoh2nbo_dens, nbo_sys_dens(1), Tvntoh_sys_dens, &
            opA='T')
       call dense_destroy(Tvntoh_sys_dens)
       call dense_destroy(Tv_dens)
       ! Construct expressions of NTOn-E in NBOs (ntoe2nbo_dens)
       !==================== Tc = <VAL-NGWFs|COND-NGWFs> ===================!
       if(pub_lr_tddft_joint_set) then                                      !
          call dense_create(Tc_dens,val_ngwf_basis%num,joint_ngwf_basis%num)!
          call dense_convert(Tc_dens, joint_rep%cross_overlap)              !
       else                                                                 !
          call dense_create(Tc_dens,val_ngwf_basis%num,cond_ngwf_basis%num) !
          call dense_convert(Tc_dens, cond_rep%cross_overlap)               !
       end if                                                               !
       !====================================================================!
       if(pub_lr_tddft_joint_set) then
          call dense_create(Tcntoe_sys_dens, val_ngwf_basis%num, num_nto_trans)
          call dense_create(ntoe2nbo_dens, val_ngwf_basis%num, num_nto_trans)
       else
          call dense_create(Tcntoe_sys_dens, val_ngwf_basis%num, num_nto_trans)
          call dense_create(ntoe2nbo_dens, val_ngwf_basis%num, num_nto_trans)
       end if
       call dense_product(Tcntoe_sys_dens, Tc_dens, ntoe_sys_dens)
       call dense_product(ntoe2nbo_dens, nbo_sys_dens(1), &
                          Tcntoe_sys_dens, opA='T')
       call dense_destroy(Tcntoe_sys_dens)
       call dense_destroy(Tc_dens)
    endif ! aux
  end subroutine internal_projto_nbo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine internal_ini_clo_ref_proj(flag)

    implicit none

    integer, intent(in) :: flag

    if( flag == 1 ) then ! initialise
       !-------------------------------------------------------------------------!
       ! Construct matrices for saving NGWF-coeff of NTOn-H(E) of the Sys and    !
       ! matrices saving those of the Ref                                        !
       ! NTOn-H(E) of the Sys will be projected to those of the Ref              !
       !-------------------------------------------------------------------------!
       call dense_create(ntoh_reftgt_dens, reftgt_size, val_ngwf_basis%num)
       call dense_create(ntoh_rctgt_dens, reftgt_size, val_ngwf_basis%num)
       if(pub_lr_tddft_joint_set) then
          call dense_create(ntoe_reftgt_dens, reftgt_size, joint_ngwf_basis%num)
          call dense_create(ntoe_rctgt_dens, reftgt_size, joint_ngwf_basis%num)
       else
          call dense_create(ntoe_reftgt_dens, reftgt_size, cond_ngwf_basis%num)
          call dense_create(ntoe_rctgt_dens, reftgt_size, cond_ngwf_basis%num)
       end if

       ! For <Sc-NTO|Sc-NTO>
       call dense_create(nto_sch2sch, systgt_size, systgt_size)
       call dense_create(nto_sce2sce, systgt_size, systgt_size)
       call dense_create(nto_sch2sce, systgt_size, systgt_size)
       call dense_create(nto_sce2sch, systgt_size, systgt_size)

       ! For <Sc-NTO|S-NTO>
       call dense_create(nto_sh2sch, systgt_size, systgt_size)
       call dense_create(nto_se2sce, systgt_size, systgt_size)
       call dense_create(nto_sh2sce, systgt_size, systgt_size)
       call dense_create(nto_se2sch, systgt_size, systgt_size)

       ! For <Rc-NTO|Rc-NTO>
       call dense_create(nto_rch2rch, reftgt_size, reftgt_size)
       call dense_create(nto_rce2rce, reftgt_size, reftgt_size)
       call dense_create(nto_rch2rce, reftgt_size, reftgt_size)
       call dense_create(nto_rce2rch, reftgt_size, reftgt_size)

       ! For <Rc-NTO|Sc-NTO>
       call dense_create(nto_sch2rch, reftgt_size, systgt_size)
       call dense_create(nto_sce2rce, reftgt_size, systgt_size)
       call dense_create(nto_sch2rce, reftgt_size, systgt_size)
       call dense_create(nto_sce2rch, reftgt_size, systgt_size)

       ! Allocate SPAM3 structures for ref overlap matrix, kinetic energy,
       ! density kernel, inverse overlap
       call ngwf_rep_create(val_s2r_rep,'',mdl)
       call ngwf_rep_create(val_s2s_rep,'',mdl)
       call ngwf_rep_create(val_r2r_rep,'',mdl)
       if(.not.pub_lr_tddft_joint_set) then
          call ngwf_rep_create(cj_s2r_rep,'c',mdl)
          call ngwf_rep_create(cj_s2s_rep,'c',mdl)
          call ngwf_rep_create(cj_r2r_rep,'c',mdl)
       else
          call ngwf_rep_create(cj_s2r_rep,'j',mdl)
          call ngwf_rep_create(cond_s2r_rep,'c',mdl)
          call ngwf_rep_create(cj_s2s_rep,'j',mdl)
          call ngwf_rep_create(cj_r2r_rep,'j',mdl)
       end if
       ! Allocate storage for NGWFs
       call data_functions_alloc(val_s2r_rep%ngwfs_on_grid(1), val_ngwf_basis%&
            n_ppds*mdl%cell%n_pts, iscmplx=val_rep%ngwfs_on_grid(1)%iscmplx)
       if(pub_lr_tddft_joint_set) then
          call data_functions_alloc(cj_s2r_rep%ngwfs_on_grid(1), joint_ngwf_basis%&
               n_ppds*mdl%cell%n_pts, iscmplx=val_rep%ngwfs_on_grid(1)%iscmplx)
          call data_functions_alloc(cond_s2r_rep%ngwfs_on_grid(1), cond_ngwf_basis%&
               n_ppds*mdl%cell%n_pts, iscmplx=val_rep%ngwfs_on_grid(1)%iscmplx)
       else
          call data_functions_alloc(cj_s2r_rep%ngwfs_on_grid(1), cond_ngwf_basis%&
               n_ppds*mdl%cell%n_pts, iscmplx=val_rep%ngwfs_on_grid(1)%iscmplx)
       end if
       !--------------------------------------------------------------------!
       ! Build up needed (cross-)overlap matrices for Sys-to-Sys projection !
       !--------------------------------------------------------------------!
       if(pub_lr_tddft_joint_set) then
          call function_ops_brappd_ketppd(val_s2s_rep%overlap%p, &
               val_rep%ngwfs_on_grid(1), val_ngwf_basis, &
               val_rep%ngwfs_on_grid(1), val_ngwf_basis, mdl%cell)
          call function_ops_brappd_ketppd(cj_s2s_rep%overlap%p, &
               joint_rep%ngwfs_on_grid(1), joint_ngwf_basis, &
               joint_rep%ngwfs_on_grid(1), joint_ngwf_basis, mdl%cell)
          call function_ops_brappd_ketppd(cj_s2s_rep%cross_overlap_tr%p, &
               joint_rep%ngwfs_on_grid(1), joint_ngwf_basis, &
               val_rep%ngwfs_on_grid(1), val_ngwf_basis, mdl%cell)
          call function_ops_brappd_ketppd(cj_s2s_rep%cross_overlap%p, &
               val_rep%ngwfs_on_grid(1), val_ngwf_basis, &
               joint_rep%ngwfs_on_grid(1), joint_ngwf_basis, mdl%cell)
       else
          call function_ops_brappd_ketppd(val_s2s_rep%overlap%p, &
               val_rep%ngwfs_on_grid(1), val_ngwf_basis, &
               val_rep%ngwfs_on_grid(1), val_ngwf_basis, mdl%cell)
          call function_ops_brappd_ketppd(cj_s2s_rep%overlap%p, &
               cond_rep%ngwfs_on_grid(1), cond_ngwf_basis, &
               cond_rep%ngwfs_on_grid(1), cond_ngwf_basis, mdl%cell)
          call function_ops_brappd_ketppd(cj_s2s_rep%cross_overlap_tr%p, &
               cond_rep%ngwfs_on_grid(1), cond_ngwf_basis, &
               val_rep%ngwfs_on_grid(1), val_ngwf_basis, mdl%cell)
          call function_ops_brappd_ketppd(cj_s2s_rep%cross_overlap%p, &
               val_rep%ngwfs_on_grid(1), val_ngwf_basis, &
               cond_rep%ngwfs_on_grid(1), cond_ngwf_basis, mdl%cell)
       end if !pub_lr_tddft_joint_set
    else ! close
       call data_functions_dealloc(val_s2r_rep%ngwfs_on_grid(1))
       call data_functions_dealloc(cj_s2r_rep%ngwfs_on_grid(1))
       if(pub_lr_tddft_joint_set) then
          call data_functions_dealloc(cond_s2r_rep%ngwfs_on_grid(1))
       end if

       call dense_destroy(ntoh_reftgt_dens)
       call dense_destroy(ntoh_systgt_dens)
       call dense_destroy(ntoe_reftgt_dens)
       call dense_destroy(ntoe_systgt_dens)

       call dense_destroy(ntoh_rctgt_dens)
       call dense_destroy(ntoh_sctgt_dens)
       call dense_destroy(ntoe_rctgt_dens)
       call dense_destroy(ntoe_sctgt_dens)

       call dense_destroy(nto_sch2sch)
       call dense_destroy(nto_sch2sce)
       call dense_destroy(nto_sce2sch)
       call dense_destroy(nto_sce2sce)

       call dense_destroy(nto_sh2sch)
       call dense_destroy(nto_se2sce)
       call dense_destroy(nto_sh2sce)
       call dense_destroy(nto_se2sch)

       call dense_destroy(nto_rch2rch)
       call dense_destroy(nto_rch2rce)
       call dense_destroy(nto_rce2rch)
       call dense_destroy(nto_rce2rce)

       call dense_destroy(nto_sch2rch)
       call dense_destroy(nto_sch2rce)
       call dense_destroy(nto_sce2rch)
       call dense_destroy(nto_sce2rce)
    end if
  end subroutine internal_ini_clo_ref_proj
  end subroutine qnto_calculate


!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------- Supporting subroutines -----------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Solve for NTOs set by running SVD decomposition                      !
!----------------------------------------------------------------------!
  subroutine solve_nto_svd1(num_nto_trans, num_min_nto_trans, trans_coeff, &
             tot_perc, ntoh_sys_dens, ntoe_sys_dens, & !output
             response_kernel, val_rep, val_ngwf_basis, cj_rep, cj_ngwf_basis)

    use comms, only: pub_on_root
    use constants, only: stdout
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
        dense_product, dense_transpose, dense_mat_sqrt, dense_invert, &
        dense_get_col, dense_put_col, dense_svdsolve
    use function_basis, only: FUNC_BASIS
    use ngwf_representation, only: NGWF_REP
    use rundat, only: PUB_1K, pub_debug_on_root
    use sparse_embed, only: SPAM3_EMBED
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: num_nto_trans
    integer, intent(inout) :: num_min_nto_trans
    real(kind=DP), allocatable, intent(inout) :: trans_coeff(:)
    real(kind=DP), intent(inout) :: tot_perc
    type(DEM), intent(inout) :: ntoh_sys_dens
    type(DEM), intent(inout) :: ntoe_sys_dens
    type(SPAM3_EMBED), intent(in) :: response_kernel
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis
    type(NGWF_REP), intent(in) :: cj_rep
    type(FUNC_BASIS), intent(in) :: cj_ngwf_basis

    ! Local variables for calculating SVD to obtain NTOs
    type(DEM) :: K_dens, BK_dens, BKAt_dens
    type(DEM) :: Sv_dens, Sc_dens
    type(DEM) :: leftvecs, rightvecs, rightvecs_tr
    type(DEM) :: tmp_h_dens, tmp_e_dens
    real(kind=DP), allocatable :: singulars(:)
    real(kind=DP), allocatable :: ntoh_sys_vec(:), ntoe_sys_vec(:)
    integer :: n
    integer :: ierr

    if (pub_debug_on_root) write(stdout,'(a)') 'In solve_nto_svd1.'

    ! Start timer
    call timer_clock('qnto_SVD1', 1)

    call dense_create(K_dens, cj_ngwf_basis%num, val_ngwf_basis%num)
    call dense_convert(K_dens, response_kernel)

    call dense_create(Sv_dens, val_ngwf_basis%num, val_ngwf_basis%num)
    call dense_create(Sc_dens, cj_ngwf_basis%num, cj_ngwf_basis%num)
    call dense_convert(Sv_dens, val_rep%overlap)
    call dense_convert(Sc_dens, cj_rep%overlap)

    ! Solve Sv = At * A and Sc = Bt * B
    call dense_mat_sqrt(Sv_dens)
    call dense_mat_sqrt(Sc_dens)

    ! Construct B * K * At to be SVD-decomposed as BU * Sig * (AV)t
    ! (BU)t * (BU) = 1 and (AV)t * (AV) = 1
    call dense_create(BK_dens, cj_ngwf_basis%num, val_ngwf_basis%num)
    call dense_product(BK_dens, Sc_dens, K_dens)
    call dense_destroy(K_dens)
    call dense_create(BKAt_dens, cj_ngwf_basis%num, val_ngwf_basis%num)
    call dense_product(BKAt_dens, BK_dens, Sv_dens)
    call dense_destroy(BK_dens)

    call dense_create(leftvecs, cj_ngwf_basis%num, val_ngwf_basis%num)
    call dense_create(rightvecs, val_ngwf_basis%num, val_ngwf_basis%num)
    allocate(singulars(val_ngwf_basis%num), stat=ierr)
    call utils_alloc_check('qnto/solve_nto_svd1', 'singulars', ierr)

    call dense_svdsolve(singulars, BKAt_dens, leftvecs, rightvecs)
    call dense_destroy(BKAt_dens)

    ! Calculate A^-1 and B^-1
    call dense_invert(Sv_dens)
    call dense_invert(Sc_dens)

    call dense_create(tmp_h_dens, val_ngwf_basis%num, val_ngwf_basis%num)
    call dense_create(tmp_e_dens, cj_ngwf_basis%num, val_ngwf_basis%num)

    ! U = B^-1 * BU; V = A^-1 * AV
    call dense_product(tmp_e_dens, Sc_dens, leftvecs)
    call dense_create(rightvecs_tr, val_ngwf_basis%num, val_ngwf_basis%num)
    call dense_transpose(rightvecs_tr, rightvecs)
    call dense_product(tmp_h_dens, Sv_dens, rightvecs_tr)
    call dense_destroy(Sv_dens)
    call dense_destroy(Sc_dens)

    tot_perc = 0.0_DP
    do n=1, num_nto_trans
       tot_perc = tot_perc + singulars(n)**2
       if(tot_perc>=0.95_DP) then
          num_min_nto_trans = n
          exit
       end if
    end do

    ! calculate tot_perc for output
    tot_perc = 0.0_DP
    do n=1, num_nto_trans
       tot_perc = tot_perc + singulars(n)**2
    end do

    allocate(trans_coeff(num_nto_trans), stat=ierr)
    call utils_dealloc_check('qnto/solve_nto_svd1', 'trans_coeff', ierr)
    do n=1, num_nto_trans
       trans_coeff(n) = singulars(n)
       if (pub_debug_on_root) write(stdout,*) 'Singular value:', trans_coeff(n)
    end do

    deallocate(singulars, stat=ierr)
    call utils_dealloc_check('qnto/solve_nto_svd1', 'singulars', ierr)
    call dense_destroy(rightvecs)
    call dense_destroy(rightvecs_tr)
    call dense_destroy(leftvecs)

    call dense_create(ntoh_sys_dens, val_ngwf_basis%num, num_nto_trans)
    allocate(ntoh_sys_vec(val_ngwf_basis%num), stat=ierr)
    call utils_dealloc_check('qnto/solve_nto_svd1', 'ntoh_sys_vec', ierr)
    do n=1, num_nto_trans
       call dense_get_col(ntoh_sys_vec, tmp_h_dens, n)
       call dense_put_col(ntoh_sys_vec, ntoh_sys_dens, n)
    end do
    call dense_destroy(tmp_h_dens)
    deallocate(ntoh_sys_vec, stat=ierr)
    call utils_dealloc_check('qnto/solve_nto_svd1', 'ntoh_sys_vec', ierr)

    call dense_create(ntoe_sys_dens, cj_ngwf_basis%num, num_nto_trans)
    allocate(ntoe_sys_vec(cj_ngwf_basis%num), stat=ierr)
    call utils_dealloc_check('qnto/solve_nto_svd1', 'ntoe_sys_vec', ierr)
    do n=1, num_nto_trans
       call dense_get_col(ntoe_sys_vec, tmp_e_dens, n)
       call dense_put_col(ntoe_sys_vec, ntoe_sys_dens, n)
    end do
    call dense_destroy(tmp_e_dens)
    deallocate(ntoe_sys_vec, stat=ierr)
    call utils_dealloc_check('qnto/solve_nto_svd1', 'ntoe_sys_vec', ierr)

    call timer_clock('qnto_SVD1', 2)
  end subroutine solve_nto_svd1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------!
! Solve for NTOs set by diagonalising SvKtScKSv_dens and ScKSvKtSc_dens!
! and determine the relative phases among NTOs of interest             !
!----------------------------------------------------------------------!
  subroutine solve_nto_svd2(num_nto_trans, num_min_nto_trans, trans_coeff, &
             tot_perc, ntoh_sys_dens, ntoe_sys_dens, & !output
             response_kernel, val_rep, val_ngwf_basis, cj_rep, cj_ngwf_basis)

    use comms, only: pub_on_root
    use constants, only: stdout
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
        dense_get_col, dense_put_col, dense_eigensolve
    use function_basis, only: FUNC_BASIS
    use ngwf_representation, only: NGWF_REP
    use rundat, only: PUB_1K, pub_debug_on_root
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_transpose_structure, &
         sparse_embed_create, sparse_embed_product, sparse_embed_destroy, &
         sparse_embed_scale, sparse_embed_transpose
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: num_nto_trans
    integer, intent(inout) :: num_min_nto_trans
    real(kind=DP), allocatable, intent(inout) :: trans_coeff(:)
    real(kind=DP), intent(inout) :: tot_perc
    type(DEM), intent(inout) :: ntoh_sys_dens
    type(DEM), intent(inout) :: ntoe_sys_dens
    type(SPAM3_EMBED), intent(in) :: response_kernel
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis
    type(NGWF_REP), intent(in) :: cj_rep
    type(FUNC_BASIS), intent(in) :: cj_ngwf_basis

    ! Local variables for calculating SVD to obtain NTOs
    type(SPAM3_EMBED) :: response_kernel_trans
    type(SPAM3_EMBED) :: SvKt, ScK, SvKtScK, ScKSvKt, SvKtScKSv, ScKSvKtSc
    type(DEM) :: SvKtScKSv_dens, ScKSvKtSc_dens
    type(DEM) :: Sv_dens, Sc_dens
    type(DEM) :: tmp_h_dens, tmp_e_dens
    real(kind=DP), allocatable :: lambda_tv(:), lambda_tc(:)
    real(kind=DP), allocatable :: ntoh_sys_vec(:), ntoe_sys_vec(:)
    integer :: n
    integer :: ierr


    if (pub_debug_on_root) write(stdout,'(a)') 'In solve_nto_svd2.'

    ! Start timer
    call timer_clock('qnto_SVD2', 1)
    !--------------------------------------------------------------------!
    ! Construct sparse matrices and allocate lambda array used to solve  !
    ! for the Val and Cond NTO bases                                     !
    !--------------------------------------------------------------------!
    call sparse_embed_transpose_structure( &
         response_kernel_trans%structure, response_kernel)
    call sparse_embed_create(response_kernel_trans)
    call sparse_embed_transpose(response_kernel_trans, response_kernel)

    call sparse_embed_create(SvKt, val_rep%overlap, response_kernel_trans)
    call sparse_embed_product(SvKt, val_rep%overlap, response_kernel_trans)
    call sparse_embed_destroy(response_kernel_trans)
    call sparse_embed_create(ScK, cj_rep%overlap, response_kernel)
    call sparse_embed_product(ScK, cj_rep%overlap, response_kernel)

    call sparse_embed_create(SvKtScK, SvKt, ScK)
    call sparse_embed_product(SvKtScK, SvKt, ScK)
    call sparse_embed_create(ScKSvKt, ScK, SvKt)
    call sparse_embed_product(ScKSvKt, ScK, SvKt)
    call sparse_embed_destroy(SvKt)
    call sparse_embed_destroy(ScK)

    call sparse_embed_create(SvKtScKSv, SvKtScK, val_rep%overlap)
    call sparse_embed_product(SvKtScKSv, SvKtScK, val_rep%overlap)
    call sparse_embed_destroy(SvKtScK)

    call sparse_embed_create(ScKSvKtSc, ScKSvKt, cj_rep%overlap)
    call sparse_embed_product(ScKSvKtSc, ScKSvKt, cj_rep%overlap)
    call sparse_embed_destroy(ScKSvKt)


    !--------------------------------------------------------------------!
    ! Solve for the Val NTO basis                                        !
    !--------------------------------------------------------------------!
    allocate(lambda_tv(val_ngwf_basis%num), stat=ierr)
    call utils_alloc_check('qnto/solve_nto_svd2', 'lambda_tv', ierr)
    call dense_create(SvKtScKSv_dens, val_ngwf_basis%num, val_ngwf_basis%num)
    call dense_create(Sv_dens, val_ngwf_basis%num, val_ngwf_basis%num)
    call dense_create(tmp_h_dens, val_ngwf_basis%num, val_ngwf_basis%num)
    call sparse_embed_scale(SvKtScKSv, -1.0_DP)
    call dense_convert(SvKtScKSv_dens, SvKtScKSv)
    call sparse_embed_destroy(SvKtScKSv)
    call dense_convert(Sv_dens, val_rep%overlap)

    call dense_eigensolve(num_nto_trans, lambda_tv, SvKtScKSv_dens, &
         Sv_dens, 1, tmp_h_dens)
    lambda_tv(:) = lambda_tv(:) * (-1.0_DP)
    call dense_destroy(SvKtScKSv_dens)
    call dense_destroy(Sv_dens)

    tot_perc = 0.0_DP
    do n=1, num_nto_trans !val_rep%n_occ(1,PUB_1K)
       tot_perc = tot_perc + lambda_tv(n)
       if(tot_perc>=0.95_DP) then
          num_min_nto_trans = n
          exit
       end if
    end do

    ! calculate tot_perc for output
    tot_perc = 0.0_DP
    do n=1, num_nto_trans
       tot_perc = tot_perc + lambda_tv(n)
    end do

    allocate(trans_coeff(num_nto_trans), stat=ierr)
    call utils_dealloc_check('qnto/solve_nto_svd2', 'trans_coeff', ierr)
    do n=1, num_nto_trans
       trans_coeff(n) = sqrt(lambda_tv(n))
    end do

    deallocate(lambda_tv)
    call utils_dealloc_check('qnto/solve_nto_svd2', 'lambda_tv', ierr)


    !--------------------------------------------------------------------!
    ! Solve for the Cond NTO basis                                       !
    !--------------------------------------------------------------------!
    allocate(lambda_tc(cj_ngwf_basis%num), stat=ierr)
    call utils_alloc_check('qnto/solve_nto_svd2', 'lambda_tc', ierr)
    call dense_create(ScKSvKtSc_dens, cj_ngwf_basis%num, cj_ngwf_basis%num)
    call dense_create(Sc_dens, cj_ngwf_basis%num, cj_ngwf_basis%num)
    call dense_create(tmp_e_dens, cj_ngwf_basis%num, cj_ngwf_basis%num)

    call sparse_embed_scale(ScKSvKtSc, -1.0_DP)
    call dense_convert(ScKSvKtSc_dens, ScKSvKtSc)
    call sparse_embed_destroy(ScKSvKtSc)
    call dense_convert(Sc_dens, cj_rep%overlap)

    call dense_eigensolve(num_nto_trans, lambda_tc, ScKSvKtSc_dens, &
         Sc_dens, 1, tmp_e_dens)
    lambda_tc(:) = lambda_tc(:) * (-1.0_DP)
    call dense_destroy(ScKSvKtSc_dens)
    call dense_destroy(Sc_dens)

    if (pub_debug_on_root) then
       write(stdout,*) '************lambda(n)^2************'
       do n=1, num_nto_trans
          write(stdout,*) lambda_tc(n)
       enddo
    end if

    deallocate(lambda_tc)
    call utils_dealloc_check('qnto/solve_nto_svd2', 'lambda_tc', ierr)


!    call internal_solve_nto_phase1 !norm check dense
    call internal_solve_nto_phase2 !norm check sparse
!    call internal_solve_nto_phase3 !element check


    ! Construct ntoh_sys_dens
    call dense_create(ntoh_sys_dens, val_ngwf_basis%num, num_nto_trans)
    allocate(ntoh_sys_vec(val_ngwf_basis%num), stat=ierr)
    call utils_alloc_check('qnto/solve_nto_svd2', 'ntoh_sys_vec', ierr)
    do n=1, num_nto_trans
       call dense_get_col(ntoh_sys_vec, tmp_h_dens, n)
       call dense_put_col(ntoh_sys_vec, ntoh_sys_dens, n)
    end do
    call dense_destroy(tmp_h_dens)
    deallocate(ntoh_sys_vec)
    call utils_dealloc_check('qnto/solve_nto_svd2', 'ntoh_sys_vec', ierr)

    ! Construct ntoe_sys_dens
    call dense_create(ntoe_sys_dens, cj_ngwf_basis%num, num_nto_trans)
    allocate(ntoe_sys_vec(cj_ngwf_basis%num), stat=ierr)
    call utils_alloc_check('qnto/solve_nto_svd2', 'ntoe_sys_vec', ierr)
    do n=1, num_nto_trans
       call dense_get_col(ntoe_sys_vec, tmp_e_dens, n)
       call dense_put_col(ntoe_sys_vec, ntoe_sys_dens, n)
    end do
    call dense_destroy(tmp_e_dens)
    deallocate(ntoe_sys_vec)
    call utils_dealloc_check('qnto/solve_nto_svd2', 'ntoe_sys_vec', ierr)
    call timer_clock('qnto_SVD2', 2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

  subroutine internal_solve_nto_phase1

    use dense, only: dense_copy, dense_axpy, dense_product, &
        dense_rank1_update, dense_transpose, dense_get_element

    implicit none

    ! Local variables
    type(DEM) :: db_prev_dens, response_kernel_dens, dbp_dens, dbm_dens
    type(DEM) :: tmp_dens, dbpt_dens, dbmt_dens, dbpSv_dens, dbmSv_dens
    type(DEM) :: dbptSc_dens, dbmtSc_dens, dbpSvdbptSc_dens, dbmSvdbmtSc_dens
    real(kind=DP) :: Dbp_norm, Dbm_norm
    integer :: it, num_e
    real(kind=DP) :: ele

    call timer_clock('dbp_method1',1)
    !--------------------------------------------------------------------!
    ! Determine the correct transition vector of NTOn that will depend   !
    ! on the results of previous NTO pairs, i.e. NTO(n-1), ..., 1        !
    !--------------------------------------------------------------------!
    call dense_create(db_prev_dens, cj_ngwf_basis%num, val_ngwf_basis%num)
    do n=1, num_nto_trans
       ! Construct dbp (d^bar^plus) and dbm (d^bar^minus) matrices
       ! and their transposes, dbpt and dbpm, respectively
       call dense_create(response_kernel_dens, cj_ngwf_basis%num, val_ngwf_basis%num)
       call dense_convert(response_kernel_dens, response_kernel)
       call dense_create(dbp_dens, cj_ngwf_basis%num, val_ngwf_basis%num)
       call dense_create(dbm_dens, cj_ngwf_basis%num, val_ngwf_basis%num)

       call dense_copy(dbp_dens, response_kernel_dens) !step1
       call dense_destroy(response_kernel_dens)
       call dense_axpy(dbp_dens, db_prev_dens, -1.0_DP) !step2-p
       call dense_copy(dbm_dens, dbp_dens) !step2-m
       call dense_create(tmp_dens, cj_ngwf_basis%num, val_ngwf_basis%num)
       call dense_rank1_update(tmp_dens, tmp_e_dens, tmp_h_dens, n, n, &
            trans_coeff(n),skip_size_check=.true.) !step3
       call dense_axpy(dbp_dens, tmp_dens, 1.0_DP) !step4-p
       call dense_axpy(dbm_dens, tmp_dens, -1.0_DP) !step4-m
       call dense_create(dbpt_dens, val_ngwf_basis%num, cj_ngwf_basis%num)
       call dense_create(dbmt_dens, val_ngwf_basis%num, cj_ngwf_basis%num)
       call dense_transpose(dbpt_dens, dbp_dens) !step5-p
       call dense_transpose(dbmt_dens, dbm_dens) !step5-m

       ! Calculate Dbp_norm and Dbm_norm that is the norm of
       ! D^bar^plus and D^bar^minus, respectively
       call dense_create(Sv_dens, val_ngwf_basis%num, val_ngwf_basis%num)
       call dense_convert(Sv_dens, val_rep%overlap)
       call dense_create(dbpSv_dens, cj_ngwf_basis%num, val_ngwf_basis%num)
       call dense_create(dbmSv_dens, cj_ngwf_basis%num, val_ngwf_basis%num)
       call dense_product(dbpSv_dens, dbp_dens, Sv_dens) !step6-p
       call dense_product(dbmSv_dens, dbm_dens, Sv_dens) !step6-m
       call dense_destroy(Sv_dens)
       call dense_destroy(dbp_dens)
       call dense_destroy(dbm_dens)
       call dense_create(Sc_dens, cj_ngwf_basis%num, cj_ngwf_basis%num)
       call dense_convert(Sc_dens, cj_rep%overlap)
       call dense_create(dbptSc_dens, val_ngwf_basis%num, cj_ngwf_basis%num)
       call dense_create(dbmtSc_dens, val_ngwf_basis%num, cj_ngwf_basis%num)
       call dense_product(dbptSc_dens, dbpt_dens, Sc_dens) !step7-p
       call dense_product(dbmtSc_dens, dbmt_dens, Sc_dens) !step7-m
       call dense_destroy(Sc_dens)
       call dense_destroy(dbpt_dens)
       call dense_destroy(dbmt_dens)
       Dbp_norm = 0.0_DP
       Dbm_norm = 0.0_DP
       ! dense_trace(dbpSv_dens, dbptSc_dens)
       call dense_create(dbpSvdbptSc_dens, cj_ngwf_basis%num, cj_ngwf_basis%num)
       call dense_create(dbmSvdbmtSc_dens, cj_ngwf_basis%num, cj_ngwf_basis%num)
       call dense_product(dbpSvdbptSc_dens, dbpSv_dens, dbptSc_dens)
       call dense_product(dbmSvdbmtSc_dens, dbmSv_dens, dbmtSc_dens)
       call dense_destroy(dbptSc_dens)
       call dense_destroy(dbmtSc_dens)
       call dense_destroy(dbpSv_dens)
       call dense_destroy(dbmSv_dens)
       num_e = cj_ngwf_basis%num
       do it=1, num_e
          call dense_get_element(ele, dbpSvdbptSc_dens, it, it) !step8-p
          Dbp_norm = Dbp_norm + ele
          call dense_get_element(ele, dbmSvdbmtSc_dens, it, it) !step8-m
          Dbm_norm = Dbm_norm + ele
       end do
       call dense_destroy(dbpSvdbptSc_dens)
       call dense_destroy(dbmSvdbmtSc_dens)

       if (pub_debug_on_root) then
          write(stdout,'(a,i2.2,a,f18.15)') 'NTO', n, ' --> Dbp_norm:', Dbp_norm
          write(stdout,'(a,i2.2,a,f18.15)') 'NTO', n, ' --> Dbm_norm:', Dbm_norm
          write(stdout,*) '(Dbp_norm-Dbm_norm)/4.0 = ', &
               ((Dbp_norm-Dbm_norm)/4.0_DP)
       end if

       ! Obtain the correct singular vector pair
       allocate(ntoh_sys_vec(val_ngwf_basis%num), stat=ierr)
       call utils_alloc_check('qnto/internal_solve_nto_phase1', 'ntoh_sys_vec', ierr)
       if (Dbp_norm > Dbm_norm) then
          call dense_axpy(db_prev_dens, tmp_dens, 1.0_DP)
       else
          call dense_get_col(ntoh_sys_vec, tmp_h_dens, n)
          ntoh_sys_vec(:) = -1.0_DP * ntoh_sys_vec(:)
          call dense_put_col(ntoh_sys_vec, tmp_h_dens, n)
          call dense_axpy(db_prev_dens, tmp_dens, -1.0_DP)
       endif
       call dense_destroy(tmp_dens)
       deallocate(ntoh_sys_vec)
       call utils_dealloc_check('qnto/internal_solve_nto_phase1', 'ntoh_sys_vec', ierr)
    enddo !n
    call dense_destroy(db_prev_dens)

    call timer_clock('dbp_method1',2)
  end subroutine internal_solve_nto_phase1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine internal_solve_nto_phase2

    use comms, only: pub_on_root
    use dense, only: dense_rank1_update
    use sparse_embed, only: sparse_embed_copy, sparse_embed_axpy, &
        sparse_embed_trace

    implicit none

    ! Local variables
    type(SPAM3_EMBED) :: db_prev, dbp, dbm, tmp, dbpt, dbmt
    type(SPAM3_EMBED) :: dbpSv, dbmSv, dbptSc, dbmtSc
    type(DEM) :: tmp_dens
    real(kind=DP) :: Dbp_norm, Dbm_norm

    call timer_clock('dbp_method2',1)
    !--------------------------------------------------------------------!
    ! Determine the correct transition vector of NTOn that will depend   !
    ! on the results of previous NTO pairs, i.e. NTO(n-1), ..., 1        !
    !--------------------------------------------------------------------!
    call sparse_embed_create(db_prev, response_kernel)
    do n=1, num_nto_trans
       ! Construct dbp (d^bar^plus) and dbm (d^bar^minus) matrices
       ! and their transposes, dbpt and dbpm, respectively
       call sparse_embed_create(dbp, response_kernel)
       call sparse_embed_copy(dbp, response_kernel) !step1
       call sparse_embed_create(dbm, response_kernel)
       call sparse_embed_copy(dbm, response_kernel)
       call sparse_embed_axpy(dbp, db_prev, -1.0_DP) !step2-p
       call sparse_embed_copy(dbm, dbp) !step2-m

       call dense_create(tmp_dens, cj_ngwf_basis%num, val_ngwf_basis%num)
       call dense_rank1_update(tmp_dens, tmp_e_dens, tmp_h_dens, n, n, &
            trans_coeff(n),skip_size_check=.true.) !step3
       call sparse_embed_create(tmp, response_kernel)
       call dense_convert(tmp, tmp_dens)
       call dense_destroy(tmp_dens)
       call sparse_embed_axpy(dbp, tmp, 1.0_DP) !step4-p
       call sparse_embed_axpy(dbm, tmp, -1.0_DP) !step4-m
       call sparse_embed_transpose_structure(dbpt%structure, dbp)
       call sparse_embed_create(dbpt)
       call sparse_embed_transpose(dbpt, dbp) !step5-p
       call sparse_embed_transpose_structure(dbmt%structure, dbm)
       call sparse_embed_create(dbmt)
       call sparse_embed_transpose(dbmt, dbm) !step5-m

       ! Calculate Dbp_norm and Dbm_norm that is the norm of
       ! D^bar^plus and D^bar^minus, respectively
       call sparse_embed_create(dbpSv, dbp, val_rep%overlap)
       call sparse_embed_create(dbmSv, dbm, val_rep%overlap)
       call sparse_embed_product(dbpSv, dbp, val_rep%overlap) !step6-p
       call sparse_embed_product(dbmSv, dbm, val_rep%overlap) !step6-m
       call sparse_embed_destroy(dbp)
       call sparse_embed_destroy(dbm)
       call sparse_embed_create(dbptSc, dbpt, cj_rep%overlap)
       call sparse_embed_create(dbmtSc, dbmt, cj_rep%overlap)
       call sparse_embed_product(dbptSc, dbpt, cj_rep%overlap) !step7-p
       call sparse_embed_product(dbmtSc, dbmt, cj_rep%overlap) !step7-m
       call sparse_embed_destroy(dbpt)
       call sparse_embed_destroy(dbmt)
       Dbp_norm = 0.0_DP
       Dbm_norm = 0.0_DP
       call sparse_embed_trace(Dbp_norm, dbpSv, dbptSc) !step8-p
       call sparse_embed_trace(Dbm_norm, dbmSv, dbmtSc) !step8-m
       call sparse_embed_destroy(dbpSv)
       call sparse_embed_destroy(dbmSv)
       call sparse_embed_destroy(dbptSc)
       call sparse_embed_destroy(dbmtSc)

       if (pub_debug_on_root) then
          write(stdout,'(a,i2.2,a,f18.15)') 'NTO', n, ' --> Dbp_norm:', Dbp_norm
          write(stdout,'(a,i2.2,a,f18.15)') 'NTO', n, ' --> Dbm_norm:', Dbm_norm
          write(stdout,*) '(Dbp_norm-Dbm_norm)/4.0 = ', &
               ((Dbp_norm-Dbm_norm)/4.0_DP)
       end if


       ! Obtain the correct singular vector pair
       allocate(ntoh_sys_vec(val_ngwf_basis%num), stat=ierr)
       call utils_alloc_check('qnto/internal_solve_nto_phase2', 'ntoh_sys_vec', ierr)
       if (Dbp_norm > Dbm_norm) then
          call sparse_embed_axpy(db_prev, tmp, 1.0_DP)
       else
          call dense_get_col(ntoh_sys_vec, tmp_h_dens, n)
          ntoh_sys_vec(:) = -1.0_DP * ntoh_sys_vec(:)
          call dense_put_col(ntoh_sys_vec, tmp_h_dens, n)
          call sparse_embed_axpy(db_prev, tmp, -1.0_DP)
       endif
       call sparse_embed_destroy(tmp)
       deallocate(ntoh_sys_vec)
       call utils_dealloc_check('qnto/internal_solve_nto_phase2', 'ntoh_sys_vec', ierr)
    enddo !n
    call sparse_embed_destroy(db_prev)

    call timer_clock('dbp_method2',2)
  end subroutine internal_solve_nto_phase2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine internal_solve_nto_phase3

    use comms, only: pub_my_proc_id, comms_bcast
    use sparse, only: sparse_proc_of_elem, sparse_get_element, sparse_put_element
    use sparse_embed, only: sparse_embed_num_rows, sparse_embed_copy, &
         sparse_embed_axpy

    implicit none

    ! Local variables
    type(DEM) :: tmp2_h_dens, tmp2_e_dens
    type(SPAM3_EMBED) :: ntoh_sys, ntoe_sys, tmp, db
    real(kind=DP) :: el_tmp, el_h, el_e, el_db
    real(kind=DP) :: dev_p, dev_n
    integer :: num_h, num_e, ie, ih

    call timer_clock('dbp_method3',1)

    call dense_create(tmp2_h_dens, val_ngwf_basis%num, cj_ngwf_basis%num)
    allocate(ntoh_sys_vec(val_ngwf_basis%num), stat=ierr)
    call utils_alloc_check('qnto/internal_solve_nto_phase3', 'ntoh_sys_vec', ierr)
    do n=1, num_nto_trans
       call dense_get_col(ntoh_sys_vec, tmp_h_dens, n)
       call dense_put_col(ntoh_sys_vec, tmp2_h_dens, n)
    end do
    deallocate(ntoh_sys_vec)
    call utils_dealloc_check('qnto/internal_solve_nto_phase3', 'ntoh_sys_vec', ierr)

    call dense_create(tmp2_e_dens, cj_ngwf_basis%num, val_ngwf_basis%num)
    allocate(ntoe_sys_vec(cj_ngwf_basis%num), stat=ierr)
    call utils_alloc_check('qnto/internal_solve_nto_phase3', 'ntoe_sys_vec', ierr)

    do n=1, num_nto_trans
       call dense_get_col(ntoe_sys_vec, tmp_e_dens, n)
       call dense_put_col(ntoe_sys_vec, tmp2_e_dens, n)
    end do
    deallocate(ntoe_sys_vec)
    call utils_dealloc_check('qnto/internal_solve_nto_phase3', 'ntoe_sys_vec', ierr)

    call sparse_embed_transpose_structure(ntoh_sys%structure, response_kernel)
    call sparse_embed_create(ntoh_sys)
    call dense_convert(ntoh_sys, tmp2_h_dens)
    call sparse_embed_create(ntoe_sys, response_kernel)
    call dense_convert(ntoe_sys, tmp2_e_dens)

    num_h = int(sparse_embed_num_rows(ntoh_sys))
    num_e = int(sparse_embed_num_rows(ntoe_sys))

    call sparse_embed_create(tmp, response_kernel)
    call sparse_embed_create(db, response_kernel)
    call sparse_embed_copy(db, response_kernel)
    do n=1, num_nto_trans
       dev_p = 0.0_DP
       dev_n = 0.0_DP
       do ih=1, num_h
          if(sparse_proc_of_elem(n,ntoh_sys%p,'C') == pub_my_proc_id) then
             call sparse_get_element(el_h, ntoh_sys%p, ih, n)
          end if
          call comms_bcast(pub_my_proc_id, el_h)
          do ie=1, num_e
             if(sparse_proc_of_elem(n,ntoe_sys%p,'C') == pub_my_proc_id) then
                call sparse_get_element(el_e, ntoe_sys%p, ie, n)
             end if
             call comms_bcast(pub_my_proc_id, el_e)
             el_tmp = el_h * el_e * trans_coeff(n)
             ! ----------------------------------------------------------------
             !sparse_rank1_update(tmp,ntoe_sys,ntoh_sys,n,n,trans_coeff(n))
             if(sparse_proc_of_elem(ih,tmp%p,'C') == pub_my_proc_id) then
                call sparse_put_element(el_tmp, tmp%p, ie, ih)
             end if
             !-----------------------------------------------------------------

             if(sparse_proc_of_elem(ih,db%p,'C') == pub_my_proc_id) then
                call sparse_get_element(el_db, db%p, ie, ih)
             end if
             call comms_bcast(pub_my_proc_id, el_db)
             dev_p = dev_p + (el_db-el_tmp)**2
             dev_n = dev_n + (el_db+el_tmp)**2
          end do
       end do

       allocate(ntoh_sys_vec(val_ngwf_basis%num), stat=ierr)
       call utils_alloc_check('qnto/internal_solve_nto_phase3', 'ntoh_sys_vec', ierr)
       if(dev_p < dev_n) then
          call sparse_embed_axpy(db, tmp, -1.0_DP)
          call dense_get_col(ntoh_sys_vec, tmp_h_dens, n)
          ntoh_sys_vec(:) = -1.0_DP * ntoh_sys_vec(:)
          call dense_put_col(ntoh_sys_vec, tmp_h_dens, n)
          if(pub_on_root) then
             write(stdout,'(a,i3.3,a,i5,a,f18.15,a)') 'NTO', n, &
                  ' --> RMS of Dev_p for ', num_h*num_e, ' elements:', &
                  sqrt(dev_p/real(num_h*num_e)),'*'
             write(stdout,'(a,i3.3,a,i5,a,f18.15)') 'NTO', n, &
                  ' --> RMS of Dev_n for ', num_h*num_e, ' elements:', &
                  sqrt(dev_n/real(num_h*num_e))
             write(stdout,*)
          end if
       else
          call sparse_embed_axpy(db, tmp, 1.0_DP)
          if(pub_on_root) then
             write(stdout,'(a,i3.3,a,i5,a,f18.15)') 'NTO', n, &
                  ' --> RMS of Dev_p for ', num_h*num_e, ' elements:', &
                  sqrt(dev_p/real(num_h*num_e))
             write(stdout,'(a,i3.3,a,i5,a,f18.15,a)') 'NTO', n, &
                  ' --> RMS of Dev_n for ', num_h*num_e, ' elements:', &
                  sqrt(dev_n/real(num_h*num_e)),'*'
             write(stdout,*)
          end if
       end if
       deallocate(ntoh_sys_vec)
       call utils_dealloc_check('qnto/internal_solve_nto_phase3', 'ntoh_sys_vec', ierr)

       ! Prepare for the calculation of NTO(n+1)
       call sparse_embed_scale(tmp, 0.0_DP)
    end do !n

    call dense_destroy(tmp2_h_dens)
    call dense_destroy(tmp2_e_dens)
    call sparse_embed_destroy(tmp)
    call sparse_embed_destroy(ntoh_sys)
    call sparse_embed_destroy(ntoe_sys)
    call sparse_embed_destroy(db)

    call timer_clock('dbp_method3',2)
  end subroutine internal_solve_nto_phase3
  end subroutine solve_nto_svd2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_nbo_info(nbo_list, apx, num_nbo)

    use comms, only: pub_on_root
    use constants, only: stdout
    use rundat, only: pub_rootname, pub_debug_on_root
    use utils, only: utils_open_unit_check, utils_close_unit_check, &
        utils_alloc_check, utils_unit

    implicit none

    ! Argument
    type(NBO_info), allocatable, intent(inout) :: nbo_list(:)
    character(len=4) :: apx
    integer, intent(in) :: num_nbo

    ! Local variable
    integer :: num_tmp, nbo_idx
    integer :: i, j, idx, ierr, nbo_input_unit
    character(len=256) :: filename
    character(len=32) :: str_tmp1
    character(len=3) :: str_tmp2
    character(len=3) :: type
    character(len=1) :: end_line

    100 format(i1,1x,a2,a3,1x,a2,a3,1x,a2,a3,1x,f8.5,2x,f10.5)
!    write(filename,'(a)') trim(pub_rootname)//'_inittr_nao_nbo.dat'//apx
!    nbo_input_unit = utils_unit()
!    open(unit=nbo_input_unit,file=trim(filename),action='read',iostat=ierr)
!    call utils_open_unit_check('qnto/read_nbo_info',trim(filename),ierr)
!    read(nbo_input_unit,*)
!    read(nbo_input_unit,*)
!    read(nbo_input_unit,'(i6)') num_nbo
!    close(nbo_input_unit,iostat=ierr)
!    call utils_close_unit_check('qnto/read_nbo_info',trim(filename),ierr)

    if (pub_debug_on_root) write(stdout,'(a,i0)') 'num_nbo', num_nbo

    if(pub_on_root) then
       write(filename,'(a)') 'outNBO'//apx
       nbo_input_unit = utils_unit()
       open(unit=nbo_input_unit,file=trim(filename),action='read',iostat=ierr)
       call utils_open_unit_check('qnto/read_nbo_info',trim(filename),ierr)
       read(nbo_input_unit,'(a22)') str_tmp1
       do while(trim(str_tmp1) /= ' NATURAL BOND ORBITALS')
          read(nbo_input_unit,'(a22)') str_tmp1
       enddo
       do i=1, 5
          read(nbo_input_unit,*)
       enddo
    end if

    allocate(nbo_list(num_nbo), stat=ierr)
    call utils_alloc_check('qnto/read_nbo_info', 'nbo_list', ierr)
    if (pub_on_root) then
       i = 1
       do while (i <= num_nbo)
          read(nbo_input_unit, '(i4,2x,a3,1x,a1)',advance='no') nbo_idx, type, end_line
          if(type /= ' ') then
             nbo_list(i)%nbo_idx = nbo_idx
             nbo_list(i)%type = type
             read(nbo_input_unit,100) nbo_list(i)%count, &
                 nbo_list(i)%a1e, nbo_list(i)%a1id, nbo_list(i)%a2e, &
                 nbo_list(i)%a2id, nbo_list(i)%a3e, nbo_list(i)%a3id, &
                 nbo_list(i)%occ, nbo_list(i)%energy
             i = i + 1
          else if(end_line /= '-') then
             read(nbo_input_unit,*)
          else
             do j=1,9
                read(nbo_input_unit,*)
             end do
          endif

          if (pub_debug_on_root) &
               write(stdout,'(a,i0,a,i0)') 'i:', i, 'nbo_idx:', nbo_idx

       end do
       close(nbo_input_unit,iostat=ierr)
       call utils_close_unit_check('qnto/read_nbo_info',trim(filename),ierr)
    end if !pub_on_root

  end subroutine read_nbo_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------!
! calibrate_nbo                                                        !
!----------------------------------------------------------------------!
  subroutine calibrate_nbo(nbo_dens, nbo_std_dens, nbo_list, stdnbo_list,&
             Tv, num_stdnbo, num_nbo)

    use sparse_embed, only: SPAM3_EMBED
    use dense, only: DEM, dense_convert, dense_product, dense_create, &
        dense_destroy, dense_get_element
    use constants, only: stdout
    use comms, only: pub_on_root
    use rundat, only: pub_debug_on_root, pub_debug

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(in) :: Tv
    type(NBO_info), intent(inout) :: nbo_list(:), stdnbo_list(:)
    type(DEM), intent(inout) :: nbo_dens
    type(DEM), intent(in) :: nbo_std_dens
    integer, intent(in) :: num_stdnbo, num_nbo

    ! Local variables
    type(DEM) :: tmp_dens, norm_dens, Tv_dens
    integer :: std, tgt
    integer :: ierr
    real(kind=DP) :: norm
    integer :: i

    call dense_create(norm_dens, num_stdnbo, num_nbo)
    call dense_create(Tv_dens, num_stdnbo, num_nbo)
    call dense_convert(Tv_dens, Tv)
    call dense_create(tmp_dens, num_stdnbo, num_nbo)
    call dense_product(tmp_dens, Tv_dens, nbo_dens)
    call dense_product(norm_dens, nbo_std_dens, tmp_dens, opA='T')

    if (pub_debug) then
       300 format(i4,a1,a3,1x,i2,2x,a2,a3,1x,a2,a3,1x,a2,a3,1x,f8.5,4x, &
       f8.5)
       if(pub_on_root) write(stdout,*) 'DEBUG: Before NBOs calibration'
          do i=1, num_nbo
             if(pub_on_root) write(stdout,300) i, '.', nbo_list(i)%type, &
             nbo_list(i)%count, nbo_list(i)%a1e, nbo_list(i)%a1id, &
             nbo_list(i)%a2e, nbo_list(i)%a2id, nbo_list(i)%a3e, &
             nbo_list(i)%a3id, nbo_list(i)%occ, nbo_list(i)%energy
          end do
       do std=1, num_stdnbo
          do tgt=1, num_nbo
             if(nbo_list(tgt)%type/='RY*') then
                call dense_get_element(norm, norm_dens, std, tgt)
                if(pub_on_root) write(stdout,'(f5.1)', advance='no') norm
             endif
          end do
          if(pub_on_root) write(stdout,*)
       end do
    end if

    do std=1, num_stdnbo
       do tgt=std, num_nbo
          if( (nbo_list(tgt)%type==stdnbo_list(std)%type) .and. &
              (nbo_list(tgt)%a1e==stdnbo_list(std)%a1e) .and. &
             (nbo_list(tgt)%a1id==stdnbo_list(std)%a1id) .and. &
              (nbo_list(tgt)%a2e==stdnbo_list(std)%a2e) .and. &
             (nbo_list(tgt)%a2id==stdnbo_list(std)%a2id) ) then
             call dense_get_element(norm, norm_dens, std, tgt)
          else
             cycle
          end if
          if(norm >= 0.5) then
             if(tgt /= std) then
                call internal_swap_scale_nbo(tgt, std, 1.0_DP)
                exit
             end if
          else if(norm < -0.5) then
             call internal_swap_scale_nbo(tgt, std, -1.0_DP)
             exit
          end if !norm
       end do !tgt
    end do !std

    if (pub_debug_on_root) write(stdout,*) 'DEBUG: After NBOs calibration'

    if (pub_debug) then
       do i=1, num_nbo
          if (pub_on_root) write(stdout,300) i, '.', nbo_list(i)%type, &
          nbo_list(i)%count, nbo_list(i)%a1e, nbo_list(i)%a1id, &
          nbo_list(i)%a2e, nbo_list(i)%a2id, nbo_list(i)%a3e, &
          nbo_list(i)%a3id, nbo_list(i)%occ, nbo_list(i)%energy
       end do
       do std=1, num_stdnbo
          do tgt=1, num_nbo
             if (nbo_list(tgt)%type/='RY*') then
                call dense_get_element(norm, norm_dens, std, tgt)
                if(pub_on_root) write(stdout,'(f5.1)', advance='no') norm
             end if
          end do
          if (pub_on_root) write(stdout,*)
       end do
    end if

    call dense_destroy(norm_dens)
    call dense_destroy(Tv_dens)

  contains

    !----------------------------------------------------------------------!
    ! swap_scale_nbo                                                       !
    !----------------------------------------------------------------------!
    subroutine internal_swap_scale_nbo(tgt, std, scale)

      use dense, only: dense_get_col, dense_put_col

      implicit none

      integer, intent(in) :: tgt, std
      real(kind=DP), intent(in) :: scale

      ! Local variables
      real(kind=DP) :: vec1(num_nbo), vec2(num_nbo)
      real(kind=DP) :: vec3(num_stdnbo), vec4(num_stdnbo)
      type(NBO_info) :: nbo_tmp_info

      if(std /= tgt) then
         call dense_get_col(vec1, nbo_dens, tgt)
         call dense_get_col(vec2, nbo_dens, std)
         vec1(:) = vec1(:) * scale
         call dense_put_col(vec1, nbo_dens, std)
         call dense_put_col(vec2, nbo_dens, tgt)

         call dense_get_col(vec3, norm_dens, tgt)
         call dense_get_col(vec4, norm_dens, std)
         vec3(:) = vec3(:) * scale
         call dense_put_col(vec3, norm_dens, std)
         call dense_put_col(vec4, norm_dens, tgt)

         nbo_tmp_info = nbo_list(tgt)
         nbo_list(tgt) = nbo_list(std)
         nbo_list(std) = nbo_tmp_info
         nbo_list(std)%count = stdnbo_list(std)%count
      else
         call dense_get_col(vec1, nbo_dens, tgt)
         vec1(:) = vec1(:) * scale
         call dense_put_col(vec1, nbo_dens, tgt)

         call dense_get_col(vec3, norm_dens, tgt)
         vec3(:) = vec3(:) * scale
         call dense_put_col(vec3, norm_dens, tgt)
      endif

    end subroutine internal_swap_scale_nbo
  end subroutine calibrate_nbo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_nbo(nto2nbo_dens, nbo_list, num_nbo, n, qnto_output_unit, orb)
    use comms, only: pub_root_proc_id, pub_on_root
    use constants, only: stdout
    use dense, only: DEM, dense_point2point_element
    use rundat, only: pub_use_aux_ngwfs
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_heapsort

    implicit none

    ! Argument
    type(DEM), intent(in) :: nto2nbo_dens
    type(NBO_info), intent(in) :: nbo_list(:)
    integer, intent(in) :: num_nbo
    integer, intent(in) :: n
    integer, intent(in) :: qnto_output_unit
    character(LEN=1), intent(in) :: orb

    ! Local variable
    integer :: nbo_idx
    integer :: ry_skip
    real(kind=DP) :: coeff_nbo, sum_nbo
    real(kind=DP), allocatable :: abs_coeff(:), coeff_sorted(:), coeff(:)
    integer, allocatable :: idx(:)
    type(NBO_info), allocatable :: nbo_list_sorted(:)
    integer :: ierr

    200 format(i4,a1,a3,2x,a2,a3,1x,a2,a3,1x,a2,a3,1x,f8.5,4x, f8.5,4x,f8.5)
    allocate(abs_coeff(num_nbo), stat=ierr)
    call utils_alloc_check('qnto/output_nbo', 'abs_coeff', ierr)
    allocate(idx(num_nbo), stat=ierr)
    call utils_alloc_check('qnto/output_nbo', 'idx', ierr)
    allocate(nbo_list_sorted(num_nbo), stat=ierr)
    call utils_alloc_check('qnto/output_nbo', 'nbo_list_sorted', ierr)
    allocate(coeff_sorted(num_nbo), stat=ierr)
    call utils_alloc_check('qnto/output_nbo', 'coeff_sorted', ierr)
    allocate(coeff(num_nbo), stat=ierr)
    call utils_alloc_check('qnto/output_nbo', 'coeff', ierr)

    ! Sort NBOs according to projection coefficients of NTOn-H or NTOn-E
    do nbo_idx=1, num_nbo
       call dense_point2point_element(coeff_nbo, nto2nbo_dens, &
            nbo_idx, n, pub_root_proc_id)
       if(pub_on_root) then
          abs_coeff(nbo_idx) = abs(coeff_nbo)
          coeff(nbo_idx) = coeff_nbo
       end if
    end do
    if(pub_on_root) then
       call utils_heapsort(num_nbo, abs_coeff, idx)
       do nbo_idx=1, num_nbo
          nbo_list_sorted(nbo_idx) = nbo_list(idx(num_nbo+1-nbo_idx))
          coeff_sorted(nbo_idx) = coeff(idx(num_nbo+1-nbo_idx))
       end do

       ! Output coeff to stdout and file
       sum_nbo = 0.0_DP
       ry_skip = 0
       write(stdout,'(a)') '###########################################&
            &###################'
       write(stdout,'(a,i1,a)') 'NTO',n,'-'//orb//' projected to NBOs:'
       write(stdout,'(10x,a,15x,a,4x,a,4x,a)') 'NBO', &
            'Occupancy', 'Energy', 'proj-coeff'
       write(qnto_output_unit,'(a)') '###########################################&
            &###################'
       write(qnto_output_unit,'(a,i1,a)') 'NTO',n,'-'//orb//' projected to NBOs:'
       write(qnto_output_unit,'(10x,a,15x,a,4x,a,4x,a)') 'NBO', &
            'Occupancy', 'Energy', 'proj-coeff'
       do nbo_idx=1, 10
          do while (nbo_list_sorted(nbo_idx+ry_skip)%type == 'RY*')
             ry_skip = ry_skip + 1
          end do
          write(stdout,200) &
               nbo_list_sorted(nbo_idx+ry_skip)%nbo_idx, '.', &
               nbo_list_sorted(nbo_idx+ry_skip)%type, &
               nbo_list_sorted(nbo_idx+ry_skip)%a1e, &
               nbo_list_sorted(nbo_idx+ry_skip)%a1id, &
               nbo_list_sorted(nbo_idx+ry_skip)%a2e, &
               nbo_list_sorted(nbo_idx+ry_skip)%a2id, &
               nbo_list_sorted(nbo_idx+ry_skip)%a3e, &
               nbo_list_sorted(nbo_idx+ry_skip)%a3id, &
               nbo_list_sorted(nbo_idx+ry_skip)%occ, &
               nbo_list_sorted(nbo_idx+ry_skip)%energy, &
               coeff_sorted(nbo_idx+ry_skip)
          write(qnto_output_unit,200) &
               nbo_list_sorted(nbo_idx+ry_skip)%nbo_idx, '.', &
               nbo_list_sorted(nbo_idx+ry_skip)%type, &
               nbo_list_sorted(nbo_idx+ry_skip)%a1e, &
               nbo_list_sorted(nbo_idx+ry_skip)%a1id, &
               nbo_list_sorted(nbo_idx+ry_skip)%a2e, &
               nbo_list_sorted(nbo_idx+ry_skip)%a2id, &
               nbo_list_sorted(nbo_idx+ry_skip)%a3e, &
               nbo_list_sorted(nbo_idx+ry_skip)%a3id, &
               nbo_list_sorted(nbo_idx+ry_skip)%occ, &
               nbo_list_sorted(nbo_idx+ry_skip)%energy, &
               coeff_sorted(nbo_idx+ry_skip)
          sum_nbo = sum_nbo + coeff_sorted(nbo_idx+ry_skip)* &
                              coeff_sorted(nbo_idx+ry_skip)
       end do
       if(pub_use_aux_ngwfs .AND. orb=='E') then
          write(stdout,'(a,f8.4)') &
               'Total % accounted for by first 10 dominant NBOs (-RY*):', &
               sum_nbo
          write(qnto_output_unit,'(a,f8.4)') &
               'Total % accounted for by first 10 dominant NBOs (-RY*):', &
               sum_nbo
          write(stdout,'(1x,a)') repeat('-', 60)
          write(qnto_output_unit,'(1x,a)') repeat('-', 60)
       else
          write(stdout,'(a,f8.4)') &
               'Total % accounted for by first 10 dominant NBOs (-RY*):', &
               sum_nbo
          write(qnto_output_unit,'(a,f8.4)') &
               'Total % accounted for by first 10 dominant NBOs (-RY*):', &
               sum_nbo
       end if

       if(pub_use_aux_ngwfs .AND. orb=='E') then
          sum_nbo = 0.0_DP
          ry_skip = 0
          do nbo_idx=1, 10
             do while (nbo_list_sorted(nbo_idx+ry_skip)%type /= 'RY*')
                ry_skip = ry_skip + 1
             end do
             write(stdout,200) &
                  nbo_list_sorted(nbo_idx+ry_skip)%nbo_idx, '.', &
                  nbo_list_sorted(nbo_idx+ry_skip)%type, &
                  nbo_list_sorted(nbo_idx+ry_skip)%a1e, &
                  nbo_list_sorted(nbo_idx+ry_skip)%a1id, &
                  nbo_list_sorted(nbo_idx+ry_skip)%a2e, &
                  nbo_list_sorted(nbo_idx+ry_skip)%a2id, &
                  nbo_list_sorted(nbo_idx+ry_skip)%a3e, &
                  nbo_list_sorted(nbo_idx+ry_skip)%a3id, &
                  nbo_list_sorted(nbo_idx+ry_skip)%occ, &
                  nbo_list_sorted(nbo_idx+ry_skip)%energy, &
                  coeff_sorted(nbo_idx+ry_skip)
             write(qnto_output_unit,200) &
                  nbo_list_sorted(nbo_idx+ry_skip)%nbo_idx, '.', &
                  nbo_list_sorted(nbo_idx+ry_skip)%type, &
                  nbo_list_sorted(nbo_idx+ry_skip)%a1e, &
                  nbo_list_sorted(nbo_idx+ry_skip)%a1id, &
                  nbo_list_sorted(nbo_idx+ry_skip)%a2e, &
                  nbo_list_sorted(nbo_idx+ry_skip)%a2id, &
                  nbo_list_sorted(nbo_idx+ry_skip)%a3e, &
                  nbo_list_sorted(nbo_idx+ry_skip)%a3id, &
                  nbo_list_sorted(nbo_idx+ry_skip)%occ, &
                  nbo_list_sorted(nbo_idx+ry_skip)%energy, &
                  coeff_sorted(nbo_idx+ry_skip)
             sum_nbo = sum_nbo + coeff_sorted(nbo_idx+ry_skip)* &
                                 coeff_sorted(nbo_idx+ry_skip)
          end do
          write(stdout,'(a,f8.4)') &
               'Total % accounted for by first 10 dominant RY* NBOs:', &
               sum_nbo
          write(qnto_output_unit,'(a,f8.4)') &
               'Total % accounted for by first 10 dominant RY* NBOs:', &
               sum_nbo
       end if !pub_use_aux_ngwfs .AND. orb

       sum_nbo = 0.0_DP
       do nbo_idx=1, num_nbo
          sum_nbo = sum_nbo + coeff_sorted(nbo_idx)*coeff_sorted(nbo_idx)
       end do
       write(stdout,'(a,i6,1x,a,f8.4)') &
            'Total % accounted for by all', num_nbo, 'NBOs:', sum_nbo
       write(stdout,*)
       write(qnto_output_unit,'(a,i6,1x,a,f8.4)') &
            'Total % accounted for by all', num_nbo, 'NBOs:', sum_nbo
       write(qnto_output_unit,*)
       if(orb=='E') then
          write(stdout,'(a)') '###########################################&
               &###################'
          write(qnto_output_unit,'(a)') '###########################################&
               &###################'
       end if !orb
    end if !pub_on_root
    deallocate(coeff, stat=ierr)
    call utils_dealloc_check('qnto/output_nbo', 'coeff', ierr)
    deallocate(abs_coeff, stat=ierr)
    call utils_dealloc_check('qnto/output_nbo', 'abs_coeff', ierr)
    deallocate(idx, stat=ierr)
    call utils_dealloc_check('qnto/output_nbo', 'idx', ierr)
    deallocate(nbo_list_sorted, stat=ierr)
    call utils_dealloc_check('qnto/output_nbo', 'nbo_list_sorted', ierr)
    deallocate(coeff_sorted, stat=ierr)
    call utils_dealloc_check('qnto/output_nbo', 'coeff_sorted', ierr)
  end subroutine output_nbo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mo_coeff_o(I, row_num, ngwf_basis, mo_sys_dens, rep, mdl, &
             ngwf_1st_idx_inp, ngwf_1st_idx_nod, str, denskern)

    use augmentation, only: augmentation_overlap
    use constants, only: max_spins
    use comms, only: pub_on_root
    use datatypes, only: FUNCTIONS, data_functions_alloc, data_functions_dealloc
    use dense, only: DEM, dense_get_col
    use eigenstates, only: eigenstates_plot_mo
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_rootname, pub_qnto_write_orbitals, pub_aug, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_scale
    use utils, only: utils_open_unit_check, utils_close_unit_check, &
        utils_alloc_check, utils_dealloc_check, utils_unit

    implicit none

    ! Argument
    integer, intent(in) :: I, row_num
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(DEM), intent(in) :: mo_sys_dens
    type(NGWF_REP), intent(in) :: rep
    type(MODEL), intent(in) :: mdl
    integer, intent(in) :: ngwf_1st_idx_inp(:), ngwf_1st_idx_nod(:)
    character(len=*), intent(in) :: str
    type(SPAM3_EMBED_ARRAY), intent(in) :: denskern

    ! Local variable
    character(len=3) :: name_I
    character(len=256) :: mo_output_file
    integer :: mo_io_unit
    integer :: col_num
    integer :: mo_idx
    integer :: ierr
    type(FUNCTIONS) :: mo_sys_vec
    type(FUNCTIONS) :: mo_ref_vec
    type(SPAM3_EMBED) :: mo_kern(max_spins), aug_overlap


    col_num = ngwf_basis%num
    ! Allocate mo_sys_vec
    call data_functions_alloc(mo_sys_vec,col_num,iscmplx=rep%ngwfs_on_grid(1)%iscmplx)
    write(name_I,'(i3)') I
    !--------------------------------------------------------------------!
    ! Open files for outputting NGWF-coeff of NTOn-H(E) or (STD-)NBOs    !
    ! Allocate arrays used for saving these coeff later                  !
    !--------------------------------------------------------------------!
    if(pub_on_root) then
       if(str=='val' .OR. str=='cond' .OR. str=='joint') then
          write(mo_output_file,'(a)') trim(pub_rootname)//&
               '_NTOs_coeff_'//trim(adjustl(name_I))//'.dat'
       else if(str=='std-nbo') then
          write(mo_output_file,'(a)') trim(pub_rootname)//'_NBO.dat_std'
       else if(str=='nbo') then
          write(mo_output_file,'(a)') trim(pub_rootname)//'_NBO.dat'
       end if !str

       mo_output_file = adjustl(mo_output_file)
       mo_io_unit = utils_unit()

       ! Write header
       if(str=='val') then
          open(unit=mo_io_unit,form="formatted", &
               file=trim(mo_output_file),action="write",iostat=ierr)
          call utils_open_unit_check('qnto/mo_coeff_o',trim(mo_output_file),ierr)
          write(mo_io_unit,'(a,/a)') trim(pub_rootname), &
               'NTOn-H_'//trim(adjustl(name_I))// &
               ' in NGWF basis: (NTO to NGWF transformation)'
       else if(str=='cond' .OR. str=='joint') then
          open(unit=mo_io_unit,form="formatted",status="old",position="append", &
               file=trim(mo_output_file),action="write",iostat=ierr)
          call utils_open_unit_check('qnto/mo_coeff_o',trim(mo_output_file),ierr)
          write(mo_io_unit,'(/a,/a)') trim(pub_rootname), &
               'NTOn-E_'//trim(adjustl(name_I))// &
               ' in NGWF basis: (NTO to NGWF transformation)'
       else if(str=='std-nbo') then
          open(unit=mo_io_unit,form="formatted", &
               file=trim(mo_output_file),action="write",iostat=ierr)
          call utils_open_unit_check('qnto/mo_coeff_o',trim(mo_output_file),ierr)
          write(mo_io_unit,'(a,/a)') trim(pub_rootname), &
               'STD-NBO in NGWF basis: (STD-NBO to NGWF transformation)'
       else if(str=='nbo') then
          open(unit=mo_io_unit,form="formatted", &
               file=trim(mo_output_file),action="write",iostat=ierr)
          call utils_open_unit_check('qnto/mo_coeff_o',trim(mo_output_file),ierr)
          write(mo_io_unit,'(a,/a)') trim(pub_rootname), &
               'NBO in NGWF basis: (NBO to NGWF transformation)'
       end if !str
       ! Write matrix size
       write(mo_io_unit,'(1x,I5,1x,I5)',advance='no') row_num, col_num
    end if !pub_on_root


    ! 'kernel' and aug overlap part for dealing with augmentation
    if (pub_aug) then
       call sparse_embed_create(mo_kern(1), denskern%m(1,1))
       call sparse_embed_create(aug_overlap,rep%overlap)
       call augmentation_overlap(aug_overlap%p,mdl%pseudo_sp,mdl%paw_sp, &
            rep%sp_overlap%p)
    end if
    do mo_idx=1, row_num
       call dense_get_col(mo_sys_vec%d, mo_sys_dens, mo_idx)
       !--------------------------------------------------------------------!
       ! Plot NTOn-H(E) and output their NGWF-coeff                         !
       !--------------------------------------------------------------------!
       if( pub_qnto_write_orbitals .AND. mo_idx == 1 .AND. &
          (str=='val'.OR.str=='cond'.OR.str=='joint') ) then
          call internal_plot_mo(mo_idx)
       end if !str
       call internal_output_mo_coeff
    end do ! mo_idx
    if (pub_aug) then
       call sparse_embed_destroy(aug_overlap)
       call sparse_embed_destroy(mo_kern(1))
    end if


    !--------------------------------------------------------------------!
    ! Close files to which NGWF-coeff of NTOn-H(E)_I or (STD-)NBOs have  !
    ! been outputted; Destroy arrays and matrices saving these coeff     !
    ! If (read_in)==.true., close files where NGWF-coeff of Ref were     !
    ! read from, and deallocate arrays that were used                    !
    !--------------------------------------------------------------------!

    call data_functions_dealloc(mo_sys_vec)
    if(pub_on_root) then
       close(unit=mo_io_unit,iostat=ierr)
       call utils_close_unit_check('qnto/mo_coeff_o','MO NGWF-coeff output',ierr)
    end if

  contains

    subroutine internal_plot_mo(n)

      use constants, only: stdout

      implicit none

      ! Argument
      integer, intent(in) :: n

      ! Local variable
      character(len=2) :: name_n
      character(len=12) :: name_file
      character(len=256) :: outfile_mo  ! file names
      character(len=256) :: title_line   ! title line in output file

      !--------------------------------------------------------------------!
      ! Plot NTOn-H(E) or (STD)-NBOs                                       !
      !--------------------------------------------------------------------!
      write(name_n,'(i2)') n
      if(str == 'val') then
         write(name_file,*) 'NTO'//trim(adjustl(name_n))//'-H'//'_'&
              //trim(adjustl(name_I))
         ! NTOn-H of excitation I, outputted as "xxx_NTO<n>-H_<I>"
         write(outfile_mo,*) '_'//trim(adjustl(name_file))
         write(title_line,*) trim(adjustl(name_file))
      else if(str == 'cond' .OR. str == 'joint') then
         write(name_file,*) 'NTO'//trim(adjustl(name_n))//'-E'//'_'&
              //trim(adjustl(name_I))
         ! NTOn-E of excitation I, outputted as "xxx_NTO<n>-E_<I>"
         write(outfile_mo,*) '_'//trim(adjustl(name_file))
         write(title_line,*) trim(adjustl(name_file))
      else if(str == 'std-nbo') then
         write(name_file,*) 'NBO_STD'//trim(adjustl(name_n))
         write(outfile_mo,*) '_NBO_STD'//trim(adjustl(name_n))
         write(title_line,*) 'NBO_STD'//trim(adjustl(name_n))
      else if(str == 'nbo') then
         write(name_file,*) 'NBO'//trim(adjustl(name_n))
         write(outfile_mo,*) '_NBO'//trim(adjustl(name_n))
         write(title_line,*) 'NBO'//trim(adjustl(name_n))
      end if !str
      if(pub_on_root) write(stdout,'(a)',advance ='no') adjustl('| '//name_file//'  ')
      call eigenstates_plot_mo(mo_sys_vec, rep%ngwfs_on_grid(1), &
           ngwf_basis, mdl, title_line, outfile_mo, aug_overlap%p, mo_kern(1)%p)
    end subroutine internal_plot_mo

    subroutine internal_output_mo_coeff

      implicit none

      ! Local variable
      integer :: atom_idx, ngwf_idx, it, mo_idx
      integer :: ierr
      real(kind=DP), allocatable :: mo_tmp_vec(:)


      !--------------------------------------------------------------------!
      ! Write NGWF-coeff of NTOn-H(E) or (STD-)NBOs of the Sys to file     !
      !--------------------------------------------------------------------!
      allocate(mo_tmp_vec(ngwf_basis%num), stat=ierr)
      call utils_alloc_check('qnto/output_nto', 'mo_tmp_vec', ierr)
      ! Reorder NTOn-H(E)_I of system according to the atom order in the input
      ! file atom_idx refers to atom order in procs, orig_atom saves the
      ! original order of atom_idx_th atom in the input file
      if(str == 'nbo') then
         do atom_idx=1, mdl%nat
            do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions_aux
               mo_tmp_vec(ngwf_1st_idx_inp(mdl%par%orig_atom(atom_idx))+ngwf_idx-1) &
               = mo_sys_vec%d(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1)
            end do
         end do
      else if(str == 'val' .OR. str == 'std-nbo') then
         do atom_idx=1, mdl%par%nat
            do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions
               mo_tmp_vec(ngwf_1st_idx_inp(mdl%par%orig_atom(atom_idx))+ngwf_idx-1) &
               = mo_sys_vec%d(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1)
            end do
         end do
      else if(str == 'joint') then
        do atom_idx=1, mdl%par%nat
           do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions_cond &
                         +mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions
              mo_tmp_vec(ngwf_1st_idx_inp(mdl%par%orig_atom(atom_idx))+ngwf_idx-1 ) &
              = mo_sys_vec%d(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1 )
           end do
        end do
      else if(str == 'cond') then
         do atom_idx=1, mdl%par%nat
            do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions_cond
               mo_tmp_vec(ngwf_1st_idx_inp(mdl%par%orig_atom(atom_idx))+ngwf_idx-1) &
               = mo_sys_vec%d(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1)
            end do
         end do
      end if !str

      ! Write to file
      it = 0
      if(pub_on_root) then
         do mo_idx=1, ngwf_basis%num
            if(MOD(it,5) == 0) write(mo_io_unit,'(a)') ''
            write(mo_io_unit,'(1x,E21.14)',advance='no') mo_tmp_vec(mo_idx)
            it = it + 1
         end do
         write(mo_io_unit,*)
      end if
      deallocate(mo_tmp_vec, stat=ierr)
      call utils_dealloc_check('qnto/internal_output_mo_coeff', 'mo_tmp_vec', ierr)
    end subroutine internal_output_mo_coeff
  end subroutine mo_coeff_o

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mo_coeff_i(I, row_num, ngwf_basis, mdl, ngwf_1st_idx_inp, &
        ngwf_1st_idx_nod, str, ref_dir, mo_tgt_dens, mo_ctgt_dens)

    use comms, only: pub_on_root
    use constants, only: stdout
    use dense, only: DEM, dense_get_col, dense_put_element
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use model_type, only: MODEL
    use rundat, only: pub_rootname
    use utils, only: utils_open_unit_check, utils_close_unit_check, &
        utils_alloc_check, utils_dealloc_check, utils_unit
    implicit none

    ! Argument
    integer, intent(in) :: I, row_num
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(MODEL), intent(in) :: mdl
    integer, intent(in) :: ngwf_1st_idx_inp(:), ngwf_1st_idx_nod(:)
    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: ref_dir
    type(DEM), intent(inout) :: mo_tgt_dens
    type(DEM), optional, intent(inout) :: mo_ctgt_dens

    ! Local variable
    character(len=3) :: name_I
    character(len=256) :: mo_input_file
    integer :: mo_io_unit
    integer :: col_num
    integer :: exc_idx, mo_idx, icol
    integer :: ierr
    real(kind=DP), allocatable :: mo_vec(:)
    integer :: read_ngwf_num, mo_num


    col_num = ngwf_basis%num
    write(name_I,'(i3)') I

    !--------------------------------------------------------------------!
    ! Open Ref file for reading in NGWF-coeff of NTOn-H(E) or (STD-)NBOs !
    ! Allocate arrays used for building mo_reftgt_dens later             !
    !--------------------------------------------------------------------!
    if(str=='val' .OR. str=='cond' .OR. str=='joint') then
       write(mo_input_file,'(a)') '../'//trim(ref_dir)//'/'//trim(pub_rootname)&
            //'_NTOs_coeff_'//trim(adjustl(name_I))//'.dat'
    else if(str=='std-nbo') then
       write(mo_input_file,'(a)') '../'//trim(ref_dir)//'/'//trim(pub_rootname)&
            //'_NBO.dat_std'
    else if(str=='nbo') then
       write(mo_input_file,'(a)') '../'//trim(ref_dir)//'/'//trim(pub_rootname)&
            //'_NBO.dat'
    else if(str=='sys-val') then
       write(mo_input_file,'(a)') trim(pub_rootname)&
            //'_NTOs_coeff_'//trim(adjustl(name_I))//'.dat'
    else if(str=='sys-cond' .OR. str=='sys-joint') then
       write(mo_input_file,'(a)') trim(pub_rootname)&
            //'_NTOs_coeff_'//trim(adjustl(name_I))//'.dat'
    end if !str

    mo_input_file = adjustl(mo_input_file)
    mo_io_unit = utils_unit()
    open(unit=mo_io_unit,form="formatted", &
         file=trim(mo_input_file),action="read",iostat=ierr)
    call utils_open_unit_check('qnto/mo_coeff_i',trim(mo_input_file),ierr)

    read(mo_io_unit,*)
    read(mo_io_unit,*)
    read(mo_io_unit,'(i6,i6)',advance='no') mo_num, read_ngwf_num

    if(str=='sys-cond' .OR. str=='sys-joint' .OR. &
       str=='cond' .OR. str=='joint') then
       read(mo_io_unit,*)
       do mo_idx=1, mo_num
          do icol=0,read_ngwf_num/5
             read(mo_io_unit,*)
          end do
          read(mo_io_unit,*)
       end do
       read(mo_io_unit,*)
       read(mo_io_unit,*)
       read(mo_io_unit,'(6x,i6)',advance='no') read_ngwf_num
    end if

    if(str=='sys-val' .OR. str=='sys-cond' .OR. str=='sys-joint') then
       ! Allocate mo_sys_vec
       allocate(mo_vec(col_num), stat=ierr)
       call utils_alloc_check('qnto/mo_coeff_i', 'mo_sys_vec', ierr)
    else
       ! Allocate mo_ref_vec
       if(str /= 'nbo') then
          allocate(mo_vec(col_num), stat=ierr)
          call utils_alloc_check('qnto/mo_coeff_i', 'mo_ref_vec', ierr)
       else
          allocate(mo_vec(read_ngwf_num), stat=ierr)
          call utils_alloc_check('qnto/mo_coeff_i', 'mo_ref_vec', ierr)
       end if
    end if

    exc_idx = (I-1)*row_num

    if(str=='sys-val' .OR. str=='sys-cond' .OR. str=='sys-joint') then
       !--------------------------------------------------------------------!
       ! Build mo_systgt_dens and mo_sctgt_dens for NTOn-H(E)               !
       !--------------------------------------------------------------------!
       call internal_systgt_sctgt_build
    else
       !--------------------------------------------------------------------!
       ! Build mo_reftgt_dens for NTOn-H(E) or (STD-)NBOs of the Ref        !
       ! and build mo_rctgt_dens for NTOn-H(E)                              !
       !--------------------------------------------------------------------!
! For NBO, <AUX1-NGWFs|AUX2-NGWFs> has to be developed
! Not done yet
       call internal_reftgt_rctgt_build
    end if

    !--------------------------------------------------------------------!
    ! Close files to which NGWF-coeff of NTOn-H(E)_I or (STD-)NBO have   !
    ! been outputted; Destroy arrays and matrices saving these coeff     !
    ! If (read_in)==.true., close files where NGWF-coeff of the Ref      !
    ! were read from, and deallocate arrays that were used               !
    !--------------------------------------------------------------------!
    close(unit=mo_io_unit,iostat=ierr)
    call utils_close_unit_check('qnto/mo_coeff_i',trim(mo_input_file),ierr)
    deallocate(mo_vec, stat=ierr)
    call utils_dealloc_check('qnto/mo_coeff_i', 'mo_vec', ierr)

  contains

    subroutine internal_systgt_sctgt_build

       use constants, only: stdout
       use dense, only: DEM, dense_get_col, dense_put_col
       use function_basis, only: FUNC_BASIS
       use ion, only: ELEMENT
       use rundat, only: pub_qnto_num_core_atoms
       use utils, only: utils_alloc_check, utils_dealloc_check

       implicit none

       ! Local variable
       integer :: mo_io_idx, atom_idx, ngwf_idx, it, icol
       integer :: ierr
       real(kind=DP), allocatable :: mo_tmp_vec(:)


       do mo_idx=1, row_num
          ! Read NGWF coeff of NTOn-H(E) of Sys. geometry from file
          allocate(mo_tmp_vec(ngwf_basis%num), stat=ierr)
          call utils_alloc_check('qnto/internal_systgt_sctgt_build', 'mo_tmp_vec', ierr)
          it = 0
          do mo_io_idx=1, ngwf_basis%num
             if(MOD(it,5) == 0) read(mo_io_unit,*)
             read(mo_io_unit,'(1x,E21.14)',advance='no') mo_tmp_vec(mo_io_idx)
             it = it + 1
          end do
          read(mo_io_unit,*)

          ! Reorder NGWF-coeff of NTOn-H(E) or NBOs of the Sys according to the
          ! atom order saved in calculation procs
          if(str == 'sys-val') then
             do atom_idx=1, mdl%par%nat
                do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions
                   mo_vec(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1) &
                   = mo_tmp_vec(ngwf_1st_idx_inp(mdl%par%orig_atom(atom_idx))+ngwf_idx-1)
                end do
             end do
          else if(str == 'sys-joint') then
             do atom_idx=1, mdl%par%nat
                do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions_cond &
                             + mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions
                   mo_vec(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1) &
                   = mo_tmp_vec(ngwf_1st_idx_inp(mdl%par%orig_atom(atom_idx))+ngwf_idx-1)
                end do
             end do
          else if(str == 'sys-cond') then
             do atom_idx=1, mdl%nat
                do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions_cond
                   mo_vec(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1 ) &
                   = mo_tmp_vec(ngwf_1st_idx_inp(mdl%par%orig_atom(atom_idx))+ngwf_idx-1)
                end do
             end do
          end if !str
          deallocate(mo_tmp_vec, stat=ierr)
          call utils_dealloc_check('qnto/internal_systgt_sctgt_build', 'mo_tmp_vec', ierr)

          call dense_put_col(mo_vec, mo_tgt_dens, exc_idx+mo_idx)

          ! Set the NGWF-coeff of non-projected atoms to zero
          if(str == 'sys-val') then
             do atom_idx=1, mdl%par%nat
                if(mdl%par%orig_atom(atom_idx) > pub_qnto_num_core_atoms) then
                   do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions
                      mo_vec(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1) = 0.0_DP
                   end do
                end if
             end do
          else if(str == 'sys-joint') then
             do atom_idx=1, mdl%par%nat
                if(mdl%par%orig_atom(atom_idx) > pub_qnto_num_core_atoms) then
                   do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions_cond&
                                 +mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions
                      mo_vec(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1) = 0.0_DP
                   end do
                end if
             end do
          else if(str == 'sys-cond') then
             do atom_idx=1, mdl%par%nat
                if(mdl%par%orig_atom(atom_idx) > pub_qnto_num_core_atoms) then
                   do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions_cond
                      mo_vec(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1) = 0.0_DP
                   end do
                end if
             end do
          end if !str
          !--------------------------------------------------------------------!
          ! Build mo_systgt_dens for NTOn-H(E)                                 !
          !--------------------------------------------------------------------!
          call dense_put_col(mo_vec, mo_ctgt_dens, exc_idx+mo_idx)
       end do !mo_idx
    end subroutine internal_systgt_sctgt_build

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine internal_reftgt_rctgt_build

       use constants, only: stdout
       use rundat, only: pub_qnto_num_core_atoms
       use dense, only: dense_put_element

       implicit none

       ! Local variable
       integer :: mo_io_idx, atom_idx, ngwf_idx, it, icol, mo_num
       integer :: ierr
       real(kind=DP), allocatable :: mo_tmp_vec(:)


       if(str /= 'nbo' .AND. str /= 'std-nbo') then
          if( ngwf_basis%num <= read_ngwf_num ) then
             mo_num = ngwf_basis%num
          else
             mo_num = read_ngwf_num
          end if
       else
          mo_num = read_ngwf_num
       end if
       do mo_idx=1, row_num
          ! Read NGWF-coeff of NTOn-H(E) or NBOs of the Ref from file
          allocate(mo_tmp_vec(mo_num), stat=ierr)
          call utils_alloc_check('qnto/internal_rctgt_build', 'mo_tmp_vec', ierr)
          it = 0
          do mo_io_idx=1, mo_num
             if(MOD(it,5) == 0) read(mo_io_unit,*)
             read(mo_io_unit,'(1x,E21.14)',advance='no') mo_tmp_vec(mo_io_idx)
             it = it + 1
          end do
          read(mo_io_unit,*)

          if(str=='val' .OR. str=='cond' .OR. str=='joint') then
             do icol=1,read_ngwf_num/5-ngwf_basis%num/5
                read(mo_io_unit,*)
             end do
          end if


          ! Reorder NGWF-coeff of NTOn-H(E) or NBOs of the Ref according to the
          ! atom order saved in calculation procs
          if(str == 'nbo') then
             do atom_idx=1, mdl%par%nat
                do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions_aux
                   mo_vec(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1) &
                   = mo_tmp_vec(ngwf_1st_idx_inp(mdl%par%orig_atom(atom_idx))+ngwf_idx-1)
                end do
             end do
          else if(str == 'val' .OR. str == 'std-nbo') then
             do atom_idx=1, mdl%par%nat
                do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions
                   mo_vec(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1 ) &
                   = mo_tmp_vec(ngwf_1st_idx_inp(mdl%par%orig_atom(atom_idx))+ngwf_idx-1)
                end do
             end do
          else if(str == 'joint') then
             do atom_idx=1, mdl%par%nat
                do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions_cond &
                             + mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions
                   mo_vec(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1) &
                   = mo_tmp_vec(ngwf_1st_idx_inp(mdl%par%orig_atom(atom_idx))+ngwf_idx-1)
                end do
             end do
          else if(str == 'cond') then
             do atom_idx=1, mdl%par%nat
                do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions_cond
                   mo_vec(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1) &
                   = mo_tmp_vec(ngwf_1st_idx_inp(mdl%par%orig_atom(atom_idx))+ngwf_idx-1)
                end do
             end do
          end if !str
          deallocate(mo_tmp_vec, stat=ierr)
          call utils_dealloc_check('qnto/internal_rctgt_build', 'mo_tmp_vec', ierr)

!          ! Save NGWF-coeff of NTOn-H(E) or NBOs of the Ref in mo_reftgt_dens
!          do icol=1, ngwf_basis%num ! dense_put_row
!            call dense_put_element(mo_vec(icol), mo_tgt_dens, &
!                 exc_idx+mo_idx, icol)
!          end do

          if(str /= 'nbo' .AND. str /= 'std-nbo') then
             ! Set the NGWF-coeff of non-projected atoms to zero
             if(str == 'val') then
                do atom_idx=1, mdl%nat
                   if(mdl%par%orig_atom(atom_idx) > pub_qnto_num_core_atoms) then
                      do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions
                         mo_vec(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1) = 0.0_DP
                      end do
                   end if
                end do
             else if(str == 'joint') then
                do atom_idx=1, mdl%par%nat
                   if(mdl%par%orig_atom(atom_idx) > pub_qnto_num_core_atoms) then
                      do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions_cond &
                                    +mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions
                         mo_vec(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1) = 0.0_DP
                      end do
                   end if
                end do
             else if(str == 'cond') then
                do atom_idx=1, mdl%par%nat
                   if(mdl%par%orig_atom(atom_idx) > pub_qnto_num_core_atoms) then
                      do ngwf_idx=1, mdl%elements(mdl%par%orig_atom(atom_idx))%nfunctions_cond
                         mo_vec(ngwf_1st_idx_nod(atom_idx)+ngwf_idx-1) = 0.0_DP
                      end do
                   end if
                end do
             end if !str

             ! Save NGWF-coeff of NTOn-H(E) or NBOs of the Ref in mo_rctgt_dens
             do icol=1, ngwf_basis%num ! dense_put_row
               call dense_put_element(mo_vec(icol), mo_ctgt_dens, &
                    exc_idx+mo_idx, icol)
             end do
          end if
       end do
    end subroutine internal_reftgt_rctgt_build
  end subroutine mo_coeff_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! jcap: this routine doesn't appear to be used anywhere... commented out

!!$  subroutine nbo_sys_ref_proj(nbo_reftgt_dens, nbo_systgt_dens, ngwf_basis, &
!!$             sys_rep, ref_nbo_rep, mdl, ref_dir)
!!$
!!$    use comms, only: pub_on_root
!!$    use dense, only: DEM, dense_create, dense_destroy, dense_product, &
!!$        dense_convert, dense_get_element
!!$    use function_basis, only: FUNC_BASIS
!!$    use function_ops, only: function_ops_brappd_ketppd
!!$    use ion, only: ELEMENT
!!$    use model_type, only : MODEL
!!$    use ngwf_representation, only: NGWF_REP, ngwf_rep_create
!!$    use restart, only: restart_ngwfs_tightbox_input
!!$    use rundat, only: pub_rootname, pub_use_aux_ngwfs
!!$    use utils, only: utils_alloc_check, utils_unit, utils_open_unit_check, &
!!$        utils_close_unit_check
!!$
!!$    implicit none
!!$
!!$    ! Argument
!!$    type(DEM), intent(in) :: nbo_reftgt_dens, nbo_systgt_dens
!!$    type(FUNC_BASIS), intent(in) :: ngwf_basis
!!$    type(NGWF_REP), intent(in) :: sys_rep
!!$    type(NGWF_REP), intent(inout) :: ref_nbo_rep
!!$    type(MODEL), intent(inout) :: mdl
!!$    character(len=*), intent(in) :: ref_dir
!!$
!!$    ! Local variable
!!$    type(DEM) :: nbo_proj_dens
!!$    integer :: ierr
!!$    type(DEM) :: Svrs_dens, Srs_tmp
!!$    character(len=256) :: sr_out_filename
!!$    integer :: sr_out_unit
!!$    integer :: nbo_idx
!!$    integer :: irow, icol
!!$    real(kind=DP) :: mtx_el
!!$
!!$
!!$    call dense_create(nbo_proj_dens, ngwf_basis%num, ngwf_basis%num)
!!$    call dense_create(Svrs_dens, ngwf_basis%num, ngwf_basis%num)
!!$    !--------------------------------------------------------------------!
!!$    ! Initialise REF-NGWFs                                               !
!!$    !--------------------------------------------------------------------!
!!$    if(.not.pub_use_aux_ngwfs) then
!!$       ! Initialise the REF-NGWF values by loading from a file
!!$       call restart_ngwfs_tightbox_input(ref_nbo_rep%ngwfs_on_grid, & ! output
!!$            ngwf_basis, mdl%cell, mdl%fftbox, mdl%elements, &
!!$            'tightbox_ngwfs_', mdl%regions(1), '../'//ref_dir)
!!$       call function_ops_brappd_ketppd(ref_nbo_rep%overlap, &
!!$            ref_nbo_rep%ngwfs_on_grid, ngwf_basis, sys_rep%ngwfs_on_grid, &
!!$            ngwf_basis, mdl%cell)
!!$       ! Svrs = <VAL_REF-NGWFs|VAL_SYS-NGWFs>
!!$       call dense_convert(Svrs_dens, ref_nbo_rep%overlap)
!!$    else
!!$       ! Svrs = <AUX_REF-NGWFs|AUX_SYS-NGWFs>
!!$       call dense_convert(Svrs_dens, sys_rep%overlap)
!!$    end if
!!$
!!$    call dense_create(Srs_tmp, ngwf_basis%num, ngwf_basis%num)
!!$    call dense_product(Srs_tmp, nbo_reftgt_dens, Svrs_dens)
!!$    call dense_product(nbo_proj_dens, Srs_tmp, nbo_systgt_dens)
!!$    call dense_destroy(Srs_tmp)
!!$    call dense_destroy(Svrs_dens)
!!$
!!$
!!$    !-------------------------------------------------------------------------!
!!$    ! Perform the Sys-to-Ref projection for NBOs and output to file           !
!!$    !-------------------------------------------------------------------------!
!!$    if(pub_on_root) then
!!$       ! Open files for outputting NBO NGWF-coeff
!!$       sr_out_filename = trim(pub_rootname)//'_NBO_s2r_'//&
!!$                         trim(ref_dir)//'.dat'
!!$       sr_out_filename = adjustl(sr_out_filename)
!!$       sr_out_unit = utils_unit()
!!$
!!$       open(unit=sr_out_unit,form="formatted", &
!!$            file=trim(sr_out_filename),action="write",iostat=ierr)
!!$       call utils_open_unit_check('qnto/nbo_sys_ref_proj', &
!!$            trim(sr_out_filename),ierr)
!!$       ! Write header
!!$       write(sr_out_unit,'(a,/a)') trim(pub_rootname), &
!!$            'NBO_sys(col idx) to NBO_ref(row idx):'
!!$       write(sr_out_unit,'(5x)',advance='no')
!!$       do nbo_idx=1,ngwf_basis%num
!!$          write(sr_out_unit,'(3x,i5.5)',advance='no') nbo_idx
!!$       end do
!!$       write(sr_out_unit,*)
!!$    end if !pub_on_root
!!$    ! Output the NBO projection results to file
!!$    do irow=1,ngwf_basis%num
!!$       if(pub_on_root) write(sr_out_unit,'(i5.5)',advance='no') irow
!!$       do icol=1,ngwf_basis%num
!!$          call dense_get_element(mtx_el,nbo_proj_dens,irow,icol)
!!$          if(pub_on_root) write(sr_out_unit,'(2x,f6.2)', advance='no') mtx_el
!!$       end do
!!$       if(pub_on_root) write(sr_out_unit,*)
!!$    end do
!!$    if(pub_on_root) then
!!$       close(unit=sr_out_unit,iostat=ierr)
!!$       call utils_close_unit_check('qnto/nbo_sys_ref_proj', &
!!$            trim(sr_out_filename),ierr)
!!$    end if
!!$
!!$    call dense_destroy(nbo_proj_dens)
!!$  end subroutine nbo_sys_ref_proj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nto_sys_ref_proj( &
             nto_sch2sch, nto_sce2sce, nto_sch2sce, nto_sce2sch, &
             nto_sh2sch, nto_se2sce, nto_sh2sce, nto_se2sch, &
             nto_rch2rch, nto_rce2rce, nto_rch2rce, nto_rce2rch, &
             nto_sch2rch, nto_sce2rce, nto_sch2rce, nto_sce2rch, &
             ntoh_rctgt_dens, ntoh_sctgt_dens, ntoe_rctgt_dens, &
             ntoe_sctgt_dens, ntoh_reftgt_dens, ntoh_systgt_dens, &
             ntoe_reftgt_dens, ntoe_systgt_dens, val_ngwf_basis, &
             cj_ngwf_basis, val_rep, cj_rep, systgt_size, reftgt_size, &
             mdl, ref_dir, &
             val_s2r_rep, cj_s2r_rep, val_r2r_rep, cj_r2r_rep, &
             val_s2s_rep, cj_s2s_rep, cond_s2r_rep, cond_ngwf_basis)

    use augmentation, only: augmentation_overlap
    use comms, only: pub_on_root
    use constants, only: stdout
    use dense, only: DEM, dense_create, dense_product, dense_convert, &
        dense_transpose, dense_destroy!, dense_copy
    use function_basis, only: FUNC_BASIS, function_basis_allocate, &
        function_basis_distribute, function_basis_init_spheres, &
        function_basis_init_tight_boxes, function_basis_deallocate
    use function_ops, only: function_ops_brappd_ketppd
    use ion, only: ELEMENT
    use model_type, only: MODEL
    use ngwfs, only: ngwfs_merge_sets
    use ngwf_representation, only: NGWF_REP
    use restart, only: restart_ngwfs_tightbox_input
    use rundat, only: pub_lr_tddft_joint_set, pub_rootname, pub_aug
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_axpy
    use sparse, only: sparse_init_blocking_scheme, BLKS_COND
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
                     utils_open_unit_check, utils_close_unit_check

    implicit none

    ! Argument
    type(DEM), intent(inout) :: nto_sch2sch, nto_sce2sce, nto_sch2sce, nto_sce2sch
    type(DEM), intent(inout) :: nto_sh2sch, nto_se2sce, nto_sh2sce, nto_se2sch
    type(DEM), intent(inout) :: nto_rch2rch, nto_rce2rce, nto_rch2rce, nto_rce2rch
    type(DEM), intent(inout) :: nto_sch2rch, nto_sce2rce, nto_sch2rce, nto_sce2rch
    type(DEM), intent(in) :: ntoh_rctgt_dens, ntoh_sctgt_dens
    type(DEM), intent(in) :: ntoe_rctgt_dens, ntoe_sctgt_dens
    type(DEM), intent(in) :: ntoh_reftgt_dens, ntoh_systgt_dens
    type(DEM), intent(in) :: ntoe_reftgt_dens, ntoe_systgt_dens
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis
    type(FUNC_BASIS), intent(in) :: cj_ngwf_basis
    type(NGWF_REP), intent(in) :: val_rep, cj_rep
    integer, intent(in) :: systgt_size, reftgt_size
    type(MODEL), intent(inout) :: mdl
    character(len=*), intent(in) :: ref_dir
    type(FUNC_BASIS), optional, intent(in) :: cond_ngwf_basis
    type(NGWF_REP), intent(inout) :: val_s2r_rep, cj_s2r_rep
    type(NGWF_REP), intent(inout) :: val_r2r_rep, cj_r2r_rep
    type(NGWF_REP), intent(inout) :: val_s2s_rep, cj_s2s_rep
    type(NGWF_REP), intent(inout) :: cond_s2r_rep

    ! Local variable
    integer :: ierr
    type(DEM) :: Svrs_dens, Scrs_dens, Stcrvs_dens, Stvrcs_dens, Srs_tmp
    type(DEM) :: Svrr_dens, Scrr_dens, Stcrvr_dens, Stvrcr_dens, Srr_tmp
    type(DEM) :: Svss_dens, Scss_dens, Stcsvs_dens, Stvscs_dens, Sss_tmp
    type(DEM) :: ntoh_sctgt_tr_dens, ntoe_sctgt_tr_dens
    character(len=256) :: ngwfs_i_name, ngwfs_o_name
    integer :: ngwfs_i_unit, ngwfs_o_unit, i, read_num, nx, ny, nz
    real(kind=DP), allocatable :: box(:,:,:)
    real(kind=DP) :: pos(3)
    type(SPAM3_EMBED) :: aug_overlap



    !--------------------------------------------------------------------!
    ! Initialise NGWFs for Sys-to-Ref projection                         !
    !--------------------------------------------------------------------!
    ! Initialise the REF-NGWF values by loading from a file
    call restart_ngwfs_tightbox_input(val_s2r_rep%ngwfs_on_grid(1), &
         val_ngwf_basis, mdl%cell, mdl%fftbox, mdl%elements, 'tightbox_ngwfs', &
         mdl%regions(1), '../'//ref_dir)

    if(.not.pub_lr_tddft_joint_set) then
       call restart_ngwfs_tightbox_input(cj_s2r_rep%ngwfs_on_grid(1), &
            cj_ngwf_basis, mdl%cell, mdl%fftbox, mdl%elements, &
            'tightbox_ngwfs_cond', mdl%regions(1), '../'//ref_dir)
    else ! The passed-in cj_ngwf_basis is the joint_ngwf_basis
       call restart_ngwfs_tightbox_input(cond_s2r_rep%ngwfs_on_grid(1), &
            cond_ngwf_basis, mdl%cell, mdl%fftbox, mdl%elements, &
            'tightbox_ngwfs_cond', mdl%regions(1), '../'//ref_dir)
       call ngwfs_merge_sets(cj_s2r_rep%ngwfs_on_grid(1), cj_ngwf_basis, mdl%cell, &
            val_s2r_rep%ngwfs_on_grid(1), val_ngwf_basis, &
            cond_s2r_rep%ngwfs_on_grid(1), cond_ngwf_basis, &
            mdl%par)
    end if ! pub_lr_tddft_joint_set

    !--------------------------------------------------------------------!
    ! Build up needed (cross-)overlap matrices and do the projection     !
    !--------------------------------------------------------------------!
    call function_ops_brappd_ketppd(val_s2r_rep%overlap%p, &
         val_s2r_rep%ngwfs_on_grid(1), val_ngwf_basis, val_rep%ngwfs_on_grid(1), &
         val_ngwf_basis, mdl%cell)
    call function_ops_brappd_ketppd(cj_s2r_rep%overlap%p, &
         cj_s2r_rep%ngwfs_on_grid(1), cj_ngwf_basis, cj_rep%ngwfs_on_grid(1), &
         cj_ngwf_basis, mdl%cell)
    call function_ops_brappd_ketppd(cj_s2r_rep%cross_overlap_tr%p, &
         cj_s2r_rep%ngwfs_on_grid(1), cj_ngwf_basis, val_rep%ngwfs_on_grid(1), &
         val_ngwf_basis, mdl%cell)
    call function_ops_brappd_ketppd(cj_s2r_rep%cross_overlap%p, &
         val_s2r_rep%ngwfs_on_grid(1), val_ngwf_basis, cj_rep%ngwfs_on_grid(1), &
         cj_ngwf_basis, mdl%cell)

    ! nto_sch2rch
    ! Svrs = <VAL_REF-NGWFs|VAL_SYS-NGWFs>
    call dense_create(Svrs_dens, val_ngwf_basis%num, val_ngwf_basis%num)
    if(pub_aug) then
       call sparse_embed_create(aug_overlap, val_s2r_rep%overlap)
       call augmentation_overlap(aug_overlap%p, mdl%pseudo_sp, mdl%paw_sp, &
            val_rep%sp_overlap%p)
       call sparse_embed_axpy(aug_overlap, val_s2r_rep%overlap, 1.0_DP)
       call dense_convert(Svrs_dens, aug_overlap)
       call sparse_embed_destroy(aug_overlap)
    else
       call dense_convert(Svrs_dens, val_s2r_rep%overlap)
    end if
    call dense_create(Srs_tmp, reftgt_size, val_ngwf_basis%num)
    call dense_product(Srs_tmp, ntoh_rctgt_dens, Svrs_dens)
    call dense_product(nto_sch2rch, Srs_tmp, ntoh_sctgt_dens)
    call dense_destroy(Srs_tmp)
    call dense_destroy(Svrs_dens)

    ! nto_sce2rce
    ! Scrs = <COND_REF-NGWFs|COND_SYS-NGWFs>
    call dense_create(Scrs_dens, cj_ngwf_basis%num, cj_ngwf_basis%num)
    if(pub_aug) then
       call sparse_embed_create(aug_overlap, cj_s2r_rep%overlap)
       call augmentation_overlap(aug_overlap%p, mdl%pseudo_sp, mdl%paw_sp, &
            cj_rep%sp_overlap%p)
       call sparse_embed_axpy(aug_overlap, cj_s2r_rep%overlap, 1.0_DP)
       call dense_convert(Scrs_dens, aug_overlap)
       call sparse_embed_destroy(aug_overlap)
    else
       call dense_convert(Scrs_dens, cj_s2r_rep%overlap)
    end if
    call dense_create(Srs_tmp, reftgt_size, cj_ngwf_basis%num)
    call dense_product(Srs_tmp, ntoe_rctgt_dens, Scrs_dens)
    call dense_product(nto_sce2rce, Srs_tmp, ntoe_sctgt_dens)
    call dense_destroy(Srs_tmp)
    call dense_destroy(Scrs_dens)

    ! nto_sch2rce
    ! Stcrvs = <COND_REF-NGWFs|VAL_SYS-NGWFs>
    call dense_create(Stcrvs_dens, cj_ngwf_basis%num, val_ngwf_basis%num)
    if(pub_aug) then
       call sparse_embed_create(aug_overlap, cj_s2r_rep%cross_overlap_tr)
       call augmentation_overlap(aug_overlap%p, mdl%pseudo_sp, mdl%paw_sp, &
            cj_rep%sp_overlap%p, val_rep%sp_overlap%p)
       call sparse_embed_axpy(aug_overlap, cj_s2r_rep%cross_overlap_tr, 1.0_DP)
       call dense_convert(Stcrvs_dens, aug_overlap)
       call sparse_embed_destroy(aug_overlap)
    else
       call dense_convert(Stcrvs_dens, cj_s2r_rep%cross_overlap_tr)
    end if
    call dense_create(Srs_tmp, reftgt_size, val_ngwf_basis%num)
    call dense_product(Srs_tmp, ntoe_rctgt_dens, Stcrvs_dens)
    call dense_product(nto_sch2rce, Srs_tmp, ntoh_sctgt_dens)
    call dense_destroy(Srs_tmp)
    call dense_destroy(Stcrvs_dens)

    ! nto_sce2rch
    ! Stvrcs = <VAL_REF-NGWFs|COND_SYS-NGWFs>
    call dense_create(Stvrcs_dens, val_ngwf_basis%num, cj_ngwf_basis%num)
    if(pub_aug) then
       call sparse_embed_create(aug_overlap, cj_s2r_rep%cross_overlap)
       call augmentation_overlap(aug_overlap%p, mdl%pseudo_sp, mdl%paw_sp, &
            val_rep%sp_overlap%p, cj_rep%sp_overlap%p)
       call sparse_embed_axpy(aug_overlap, cj_s2r_rep%cross_overlap, 1.0_DP)
       call dense_convert(Stvrcs_dens, aug_overlap)
       call sparse_embed_destroy(aug_overlap)
    else
       call dense_convert(Stvrcs_dens, cj_s2r_rep%cross_overlap)
    end if
    call dense_create(Srs_tmp, reftgt_size, cj_ngwf_basis%num)
    call dense_product(Srs_tmp, ntoh_rctgt_dens, Stvrcs_dens)
    call dense_product(nto_sce2rch, Srs_tmp, ntoe_sctgt_dens)
    call dense_destroy(Srs_tmp)
    call dense_destroy(Stvrcs_dens)

    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------!
    ! Initialise NGWFs for Sys-to-Sys projection                         !
    !--------------------------------------------------------------------!
    call dense_create(ntoh_sctgt_tr_dens, systgt_size, val_ngwf_basis%num)
    call dense_create(ntoe_sctgt_tr_dens, systgt_size, cj_ngwf_basis%num)
    call dense_transpose(ntoh_sctgt_tr_dens, ntoh_sctgt_dens)
    call dense_transpose(ntoe_sctgt_tr_dens, ntoe_sctgt_dens)
    ! nto_sch2sch, nto_sh2sch
    ! Svss = <VAL_SYS-NGWFs|VAL_SYS-NGWFs>
    call dense_create(Svss_dens, val_ngwf_basis%num, val_ngwf_basis%num)
    if(pub_aug) then
       call sparse_embed_create(aug_overlap, val_s2s_rep%overlap)
       call augmentation_overlap(aug_overlap%p, mdl%pseudo_sp, mdl%paw_sp, &
            val_rep%sp_overlap%p)
       call sparse_embed_axpy(aug_overlap, val_s2s_rep%overlap, 1.0_DP)
       call dense_convert(Svss_dens, aug_overlap)
       call sparse_embed_destroy(aug_overlap)
    else
       call dense_convert(Svss_dens, val_s2s_rep%overlap)
    end if
    call dense_create(Sss_tmp, systgt_size, val_ngwf_basis%num)
    call dense_product(Sss_tmp, ntoh_sctgt_tr_dens, Svss_dens)
    call dense_product(nto_sch2sch, Sss_tmp, ntoh_sctgt_dens)
    call dense_product(nto_sh2sch, Sss_tmp, ntoh_systgt_dens)
    call dense_destroy(Sss_tmp)
    call dense_destroy(Svss_dens)

    ! nto_sce2sce, nto_se2sce
    ! Scss = <COND_SYS-NGWFs|COND_SYS-NGWFs>
    call dense_create(Scss_dens, cj_ngwf_basis%num, cj_ngwf_basis%num)
    if(pub_aug) then
       call sparse_embed_create(aug_overlap, cj_s2s_rep%overlap)
       call augmentation_overlap(aug_overlap%p, mdl%pseudo_sp, mdl%paw_sp, &
            cj_rep%sp_overlap%p)
       call sparse_embed_axpy(aug_overlap, cj_s2s_rep%overlap, 1.0_DP)
       call dense_convert(Scss_dens, aug_overlap)
       call sparse_embed_destroy(aug_overlap)
    else
       call dense_convert(Scss_dens, cj_s2s_rep%overlap)
    end if
    call dense_create(Sss_tmp, systgt_size, cj_ngwf_basis%num)
    call dense_product(Sss_tmp, ntoe_sctgt_tr_dens, Scss_dens)
    call dense_product(nto_sce2sce, Sss_tmp, ntoe_sctgt_dens)
    call dense_product(nto_se2sce, Sss_tmp, ntoe_systgt_dens)
    call dense_destroy(Sss_tmp)
    call dense_destroy(Scss_dens)

    ! nto_sch2sce, nto_sh2sce
    ! Stcsvs = <COND_SYS-NGWFs|VAL_SYS-NGWFs>
    call dense_create(Stcsvs_dens, cj_ngwf_basis%num, val_ngwf_basis%num)
    if(pub_aug) then
       call sparse_embed_create(aug_overlap, cj_s2s_rep%cross_overlap_tr)
       call augmentation_overlap(aug_overlap%p, mdl%pseudo_sp, mdl%paw_sp, &
            cj_rep%sp_overlap%p, val_rep%sp_overlap%p)
       call sparse_embed_axpy(aug_overlap, cj_s2s_rep%cross_overlap_tr, 1.0_DP)
       call dense_convert(Stcsvs_dens, aug_overlap)
       call sparse_embed_destroy(aug_overlap)
    else
       call dense_convert(Stcsvs_dens, cj_s2s_rep%cross_overlap_tr)
    end if
    call dense_create(Sss_tmp, systgt_size, val_ngwf_basis%num)
    call dense_product(Sss_tmp, ntoe_sctgt_tr_dens, Stcsvs_dens)
    call dense_product(nto_sch2sce, Sss_tmp, ntoh_sctgt_dens)
    call dense_product(nto_sh2sce, Sss_tmp, ntoh_systgt_dens)
    call dense_destroy(Sss_tmp)
    call dense_destroy(Stcsvs_dens)

    ! nto_sce2sch, nto_se2sch
    ! Stvscs = <VAL_SYS-NGWFs|COND_SYS-NGWFs>
    call dense_create(Stvscs_dens, val_ngwf_basis%num, cj_ngwf_basis%num)
    if(pub_aug) then
       call sparse_embed_create(aug_overlap, cj_s2s_rep%cross_overlap)
       call augmentation_overlap(aug_overlap%p, mdl%pseudo_sp, mdl%paw_sp, &
            val_rep%sp_overlap%p, cj_rep%sp_overlap%p)
       call sparse_embed_axpy(aug_overlap, cj_s2s_rep%cross_overlap, 1.0_DP)
       call dense_convert(Stvscs_dens, aug_overlap)
       call sparse_embed_destroy(aug_overlap)
    else
       call dense_convert(Stvscs_dens, cj_s2s_rep%cross_overlap)
    end if
    call dense_create(Sss_tmp, systgt_size, cj_ngwf_basis%num)
    call dense_product(Sss_tmp, ntoh_sctgt_tr_dens, Stvscs_dens)
    call dense_product(nto_sce2sch, Sss_tmp, ntoe_sctgt_dens)
    call dense_product(nto_se2sch, Sss_tmp, ntoe_systgt_dens)
    call dense_destroy(Sss_tmp)
    call dense_destroy(Stvscs_dens)
    call dense_destroy(ntoh_sctgt_tr_dens)
    call dense_destroy(ntoe_sctgt_tr_dens)

    !--------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------!
    ! Ref-to-Ref: Build up needed (cross-)overlap matrices and do the projection    !
    !-------------------------------------------------------------------------------!
    call function_ops_brappd_ketppd(val_r2r_rep%overlap%p, &
         val_s2r_rep%ngwfs_on_grid(1), val_ngwf_basis, &
         val_s2r_rep%ngwfs_on_grid(1), val_ngwf_basis, mdl%cell)
    call function_ops_brappd_ketppd(cj_r2r_rep%overlap%p, &
         cj_s2r_rep%ngwfs_on_grid(1), cj_ngwf_basis, &
         cj_s2r_rep%ngwfs_on_grid(1), cj_ngwf_basis, mdl%cell)
    call function_ops_brappd_ketppd(cj_r2r_rep%cross_overlap_tr%p, &
         cj_s2r_rep%ngwfs_on_grid(1), cj_ngwf_basis, &
         val_s2r_rep%ngwfs_on_grid(1), val_ngwf_basis, mdl%cell)
    call function_ops_brappd_ketppd(cj_r2r_rep%cross_overlap%p, &
         val_s2r_rep%ngwfs_on_grid(1), val_ngwf_basis, &
         cj_s2r_rep%ngwfs_on_grid(1), cj_ngwf_basis, mdl%cell)

    ! nto_rch2rch
    ! Svrr = <VAL_REF-NGWFs|VAL_REF-NGWFs>
    call dense_create(Svrr_dens, val_ngwf_basis%num, val_ngwf_basis%num)
    if(pub_aug) then
       call sparse_embed_create(aug_overlap, val_r2r_rep%overlap)
       call augmentation_overlap(aug_overlap%p, mdl%pseudo_sp, mdl%paw_sp, &
            val_rep%sp_overlap%p)
       call sparse_embed_axpy(aug_overlap, val_r2r_rep%overlap, 1.0_DP)
       call dense_convert(Svrr_dens, aug_overlap)
       call sparse_embed_destroy(aug_overlap)
    else
       call dense_convert(Svrr_dens, val_r2r_rep%overlap)
    end if
    call dense_create(Srr_tmp, reftgt_size, val_ngwf_basis%num)
    call dense_product(Srr_tmp, ntoh_rctgt_dens, Svrr_dens)
    call dense_product(nto_rch2rch, Srr_tmp, ntoh_rctgt_dens, opB='T')
    call dense_destroy(Srr_tmp)
    call dense_destroy(Svrr_dens)

    ! nto_rce2rce
    ! Scrr = <COND_REF-NGWFs|COND_REF-NGWFs>
    call dense_create(Scrr_dens, cj_ngwf_basis%num, cj_ngwf_basis%num)
    if(pub_aug) then
       call sparse_embed_create(aug_overlap, cj_r2r_rep%overlap)
       call augmentation_overlap(aug_overlap%p, mdl%pseudo_sp, mdl%paw_sp, &
            cj_rep%sp_overlap%p)
       call sparse_embed_axpy(aug_overlap, cj_r2r_rep%overlap, 1.0_DP)
       call dense_convert(Scrr_dens, aug_overlap)
       call sparse_embed_destroy(aug_overlap)
    else
       call dense_convert(Scrr_dens, cj_r2r_rep%overlap)
    end if
    call dense_create(Srr_tmp, reftgt_size, cj_ngwf_basis%num)
    call dense_product(Srr_tmp, ntoe_rctgt_dens, Scrr_dens)
    call dense_product(nto_rce2rce, Srr_tmp, ntoe_rctgt_dens, opB='T')
    call dense_destroy(Srr_tmp)
    call dense_destroy(Scrr_dens)

    ! nto_rch2rce
    ! Stcrvr = <COND_REF-NGWFs|VAL_REF-NGWFs>
    call dense_create(Stcrvr_dens, cj_ngwf_basis%num, val_ngwf_basis%num)
    if(pub_aug) then
       call sparse_embed_create(aug_overlap, cj_r2r_rep%cross_overlap_tr)
       call augmentation_overlap(aug_overlap%p, mdl%pseudo_sp, mdl%paw_sp, &
            cj_rep%sp_overlap%p, val_rep%sp_overlap%p)
       call sparse_embed_axpy(aug_overlap, cj_r2r_rep%cross_overlap_tr, 1.0_DP)
       call dense_convert(Stcrvr_dens, aug_overlap)
       call sparse_embed_destroy(aug_overlap)
    else
       call dense_convert(Stcrvr_dens, cj_r2r_rep%cross_overlap_tr)
    end if
    call dense_create(Srr_tmp, reftgt_size, val_ngwf_basis%num)
    call dense_product(Srr_tmp, ntoe_rctgt_dens, Stcrvr_dens)
    call dense_product(nto_rch2rce, Srr_tmp, ntoh_rctgt_dens, opB='T')
    call dense_destroy(Srr_tmp)
    call dense_destroy(Stcrvr_dens)

    ! nto_rce2rch
    ! Stvrcr = <VAL_REF-NGWFs|COND_REF-NGWFs>
    call dense_create(Stvrcr_dens, val_ngwf_basis%num, cj_ngwf_basis%num)
    if(pub_aug) then
       call sparse_embed_create(aug_overlap, cj_r2r_rep%cross_overlap)
       call augmentation_overlap(aug_overlap%p, mdl%pseudo_sp, mdl%paw_sp, &
            val_rep%sp_overlap%p, cj_rep%sp_overlap%p)
       call sparse_embed_axpy(aug_overlap, cj_r2r_rep%cross_overlap, 1.0_DP)
       call dense_convert(Stvrcr_dens, aug_overlap)
       call sparse_embed_destroy(aug_overlap)
    else
       call dense_convert(Stvrcr_dens, cj_r2r_rep%cross_overlap)
    end if
    call dense_create(Srr_tmp, reftgt_size, cj_ngwf_basis%num)
    call dense_product(Srr_tmp, ntoh_rctgt_dens, Stvrcr_dens)
    call dense_product(nto_rce2rch, Srr_tmp, ntoe_rctgt_dens, opB='T')
    call dense_destroy(Srr_tmp)
    call dense_destroy(Stvrcr_dens)

  end subroutine nto_sys_ref_proj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nto_output_proj_coeff( &
             nto_sch2sch, nto_sce2sch, nto_sch2sce, nto_sce2sce, &
             nto_sh2sch, nto_se2sch, nto_sh2sce, nto_se2sce, &
             nto_rch2rch, nto_rce2rch, nto_rch2rce, nto_rce2rce, &
             nto_sch2rch, nto_sce2rch, nto_sch2rce, nto_sce2rce, &
             num_states, systgt_size, reftgt_size, ref_dir, &
             global_label_sys, global_label_ref, global_match_weight_sys, &
             n_occ, num_nto_trans)

    use comms, only: pub_on_root, pub_my_proc_id, comms_bcast
    use constants, only: stdout
    use dense, only: DEM, dense_get_element
    use rundat, only: pub_rootname, pub_qnto_num_ref_states, &
        pub_qnto_num_core_atoms, pub_debug_on_root
    use utils, only: utils_open_unit_check, utils_close_unit_check, &
        utils_unit, utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Argument
    type(DEM), intent(in) :: nto_sch2sch, nto_sce2sch, nto_sch2sce, nto_sce2sce
    type(DEM), intent(in) :: nto_sh2sch, nto_se2sch, nto_sh2sce, nto_se2sce
    type(DEM), intent(in) :: nto_rch2rch, nto_rce2rch, nto_rch2rce, nto_rce2rce
    type(DEM), intent(in) :: nto_sch2rch, nto_sce2rch, nto_sch2rce, nto_sce2rce
    integer, intent(in) :: num_states, systgt_size, reftgt_size
    character(len=80), intent(in) :: ref_dir
    character(len=2), intent(inout) :: global_label_sys(:), global_label_ref(:)
    real(kind=DP), intent(inout) :: global_match_weight_sys(:)
    integer, intent(in) :: n_occ
    integer, intent(in) :: num_nto_trans

    ! Local variable
    character(len=256) :: t_coe_input_name
    character(len=256) :: sr_out_filename, ss_out_filename, rr_io_filename
    character(len=256) :: ss2_out_filename
    character(len=256) :: sr2_out_filename
    character(len=256) :: label_io_filename
    character(len=256) :: last_ref_filename
    character(len=15) :: title
    character(len=4) :: name_num_core
    character(len=2) :: name_n
    integer :: t_coe_input_unit
    integer :: sr_out_unit, ss_out_unit, rr_io_unit, ss2_out_unit, sr2_out_unit
    integer :: label_io_unit
    integer :: state_idx, state_idx_sys, state_idx2_sys, state_idx_ref
    integer :: nto_idx_sys, nto_idx2_sys, nto_idx_ref
    integer :: irow, icol, i, n
    integer :: num_match
    real(kind=DP), allocatable :: trans_coeff_sys(:), trans_coeff_ref(:)
    integer, allocatable :: I_match(:)
    integer :: ierr, reclen
    real(kind=DP) :: mtx_el, el_h, el_e
    real(kind=DP), allocatable :: nor_h_sc(:), nor_e_sc(:)
    real(kind=DP) :: h_s2sc
    real(kind=DP), allocatable :: e_s2sc(:)
    real(kind=DP), allocatable :: nor_h_s2sc(:), nor_e_s2sc(:)
    real(kind=DP), allocatable :: nor_h_rc(:), nor_e_rc(:)
    integer :: nor_h_r2rc_int, nor_e_r2rc_int
    real(kind=DP), allocatable :: nor_h_r2rc(:), nor_e_r2rc(:)
    real(kind=DP), allocatable :: h_sc2rc(:,:), e_sc2rc(:,:)
    real(kind=DP), allocatable :: match_weight(:), match_weight_sys(:)
    character(len=2), allocatable :: label_sys(:), label_ref(:)
    character(len=1), allocatable :: hash_flag_ref(:), hash_flag_sys(:)
    character(len=2), allocatable :: tmp_string(:)
    logical, allocatable :: label_flag_sys(:), label_flag_ref(:)
    real(kind=DP), allocatable :: exc_energies_sys(:), exc_energies_ref(:)
    real(kind=DP), allocatable :: match_energy_diff(:)
    character(len=1) :: mark
    character(len=4) :: str
    logical :: f_exist, ex, op


    !-------------------------------------------------------------------------!
    ! Read excitation energies and transition vector components of NTOn       !
    !-------------------------------------------------------------------------!
    ! systgt_size = num_nto_trans * num_states
    allocate(trans_coeff_sys(systgt_size),stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','trans_coeff_sys', ierr)
    allocate(exc_energies_sys(num_states),stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','exc_energies_sys', ierr)
    if(pub_on_root) then
       t_coe_input_name = trim(pub_rootname)//'_QNTO_info.dat'
       t_coe_input_unit = utils_unit()

       open(unit=t_coe_input_unit,file=trim(t_coe_input_name),action='read',iostat=ierr)
       call utils_open_unit_check('qnto/nto_output_proj_coeff',trim(t_coe_input_name),ierr)

       do state_idx_sys=1, num_states
          do
            read(t_coe_input_unit,'(a)') mark
            if( mark == '=' ) exit
          end do
          read(t_coe_input_unit,'(21x,f12.8)') exc_energies_sys(state_idx_sys)
          do nto_idx_sys=1, num_nto_trans
             n = (state_idx_sys-1)*num_nto_trans + nto_idx_sys
             read(t_coe_input_unit,'(20x,f12.8)') trans_coeff_sys(n)
          end do
       end do

       close(unit=t_coe_input_unit,iostat=ierr)
       call utils_close_unit_check('qnto/nto_output_proj_coeff',trim(t_coe_input_name),ierr)
    end if !pub_on_root

    ! reftgt_size = num_nto_trans * pub_qnto_num_ref_states
    allocate(trans_coeff_ref(reftgt_size),stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','trans_coeff_ref', ierr)
    allocate(exc_energies_ref(pub_qnto_num_ref_states),stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','exc_energies_ref', ierr)
    if(pub_on_root) then
       t_coe_input_name = '../'//trim(ref_dir)//'/'//trim(pub_rootname)//'_QNTO_info.dat'
       t_coe_input_unit = utils_unit()

       open(unit=t_coe_input_unit,file=trim(t_coe_input_name),action='read',iostat=ierr)
       call utils_open_unit_check('qnto/nto_output_proj_coeff',trim(t_coe_input_name),ierr)

       do state_idx_ref=1, pub_qnto_num_ref_states
          do
            read(t_coe_input_unit,'(a)') mark
            if( mark == '=' ) exit
          end do
          read(t_coe_input_unit,'(21x,f12.8)') exc_energies_ref(state_idx_ref)
          do nto_idx_ref=1, num_nto_trans
             n = (state_idx_ref-1)*num_nto_trans + nto_idx_ref
             read(t_coe_input_unit,'(20x,f12.8)') trans_coeff_ref(n)
          end do
       end do

       close(unit=t_coe_input_unit,iostat=ierr)
       call utils_close_unit_check('qnto/nto_output_proj_coeff',trim(t_coe_input_name),ierr)
    end if !pub_on_root


    ! Construct normalisation factors for hole and electron projection results
    allocate(nor_h_sc(systgt_size),stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','nor_h_sc', ierr)
    allocate(nor_e_sc(systgt_size),stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','nor_e_sc', ierr)
    !------------------------------------------------------------------------------------!
    ! Output cNTOn-H(E)_Sys to cNTOn-H(E)_Sys projection results                         !
    !------------------------------------------------------------------------------------!
    if(pub_on_root) then
       write(name_num_core,'(i4)') pub_qnto_num_core_atoms
       write(name_n,'(i2)') num_nto_trans
       write(ss_out_filename,'(a)') trim(pub_rootname)//'_QNTO_'//&
            trim(adjustl(name_num_core))//'a'//trim(adjustl(name_n))//'t_s2s.dat'
       ss_out_filename = adjustl(ss_out_filename)
       ss_out_unit = utils_unit()
       open(unit=ss_out_unit,form="formatted", &
            file=trim(ss_out_filename),action="write",iostat=ierr)
       call utils_open_unit_check('qnto/nto_output_proj_coeff',trim(ss_out_filename),ierr)
       ! Write header
       write(ss_out_unit,'(a,//a,/a,/a,i4,a,//a,/a,/a,/a/)') &
            'Sc2Sc:', & !'Projection of cNTOn-H(E)a_S to cNTOn-H(E)b_S:'
            &'   "<cNTOb_S|cNTOa_S>"', & ! / (<cNTOa_S|cNTOa_S>*<cNTOb_S|cNTOb_S>)^0.5"', &
            &'   using <NGWF_S|NGWF_S> on S-geom', &
            '   (qnto_num_core_atoms = ', pub_qnto_num_core_atoms, ')', &
            'Map format:', &
            '           cNTOa_S', &
            '   cNTOb_S| H2H; E2H', &
            '          | H2E; E2E'

       write(ss_out_unit,'(12x)',advance='no')
       do state_idx_sys=1, num_states
          do nto_idx_sys=1, num_nto_trans
             write(ss_out_unit,'(1x,i1,a,i3.3,1x,i1,a,i3.3)',advance='no') &
                  nto_idx_sys,'H_',state_idx_sys, nto_idx_sys,'E_',state_idx_sys
          end do
          write(ss_out_unit,'(a)',advance='no') '|'
       end do
       write(ss_out_unit,*)

       ! Output transition results of NTOs_sys
       write(ss_out_unit,'(11x,a)',advance='no') '|'
       do state_idx_sys=1, num_states
          do nto_idx_sys=1, num_nto_trans
             n = (state_idx_sys-1)*num_nto_trans + nto_idx_sys
             write(ss_out_unit,'(6x,i3.2,5x)',advance='no') &
                  nint(trans_coeff_sys(n)*100.0_DP)
          end do
          write(ss_out_unit,'(a)',advance='no') '|'
       end do
       write(ss_out_unit,*)
    end if !pub_on_root

    ! Plot the map
    do irow=1, systgt_size
       if(pub_on_root) write(ss_out_unit,'(1x,i1,a,i3.3,1x,i3.2,a)',advance='no') &
       irow-(irow-1)/num_nto_trans*num_nto_trans,'H_', &
       (irow-1)/num_nto_trans+1, nint(trans_coeff_sys(irow)*100.0_DP),'|'
       do icol=1, systgt_size
          call dense_get_element(mtx_el,nto_sch2sch,irow,icol)
          if(pub_on_root) then
             ! Construct nor_h_sc
             if(icol == irow) nor_h_sc(irow) = sqrt(mtx_el)
             if(abs(mtx_el) >= sqrt(1.0_DP/2.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '*', nint(mtx_el*100.0_DP)
             else if(abs(mtx_el) >= sqrt(1.0_DP/3.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '#', nint(mtx_el*100.0_DP)
             else if(mtx_el>-0.005_DP .AND. mtx_el<0.0_DP) then
                write(ss_out_unit,'(3x,a)',advance='no') ' -00'
             else
                write(ss_out_unit,'(3x,i4.2)',advance='no') nint(mtx_el*100.0_DP)
             end if
          end if

          call dense_get_element(mtx_el,nto_sce2sch,irow,icol)
          if(pub_on_root) then
             if(abs(mtx_el) >= sqrt(1.0_DP/2.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '*', nint(mtx_el*100.0_DP)
             else if(abs(mtx_el) >= sqrt(1.0_DP/3.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '#', nint(mtx_el*100.0_DP)
             else if(mtx_el>-0.005_DP .AND. mtx_el<0.0_DP) then
                write(ss_out_unit,'(3x,a)',advance='no') ' -00'
             else
                write(ss_out_unit,'(3x,i4.2)',advance='no') nint(mtx_el*100.0_DP)
             end if
             if(MOD(icol,num_nto_trans)==0) &
             write(ss_out_unit,'(a)',advance='no') '|'
          end if
       end do

       if(pub_on_root) write(ss_out_unit,*)
       if(pub_on_root) write(ss_out_unit,'(1x,i1,a,i3.3,4x,a)',advance='no') &
       irow-(irow-1)/num_nto_trans*num_nto_trans,'E_', &
       (irow-1)/num_nto_trans+1,'|'
       do icol=1, systgt_size
          call dense_get_element(mtx_el,nto_sch2sce,irow,icol)
          if(pub_on_root) then
             if(abs(mtx_el) >= sqrt(1.0_DP/2.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '*', nint(mtx_el*100.0_DP)
             else if(abs(mtx_el) >= sqrt(1.0_DP/3.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '#', nint(mtx_el*100.0_DP)
             else if(mtx_el>-0.005_DP .AND. mtx_el<0.0_DP) then
                write(ss_out_unit,'(3x,a)',advance='no') ' -00'
             else
                write(ss_out_unit,'(3x,i4.2)',advance='no') nint(mtx_el*100.0_DP)
             end if
          end if

          call dense_get_element(mtx_el,nto_sce2sce,irow,icol)
          if(pub_on_root) then
             ! Construct nor_e_sc
             if(icol == irow) nor_e_sc(irow) = sqrt(mtx_el)
             if(abs(mtx_el) >= sqrt(1.0_DP/2.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '*', nint(mtx_el*100.0_DP)
             else if(abs(mtx_el) >= sqrt(1.0_DP/3.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '#', nint(mtx_el*100.0_DP)
             else if(mtx_el>-0.005_DP .AND. mtx_el<0.0_DP) then
                write(ss_out_unit,'(3x,a)',advance='no') ' -00'
             else
                write(ss_out_unit,'(3x,i4.2)',advance='no') nint(mtx_el*100.0_DP)
             end if
             if(MOD(icol,num_nto_trans)==0) &
               write(ss_out_unit,'(a)',advance='no') '|'
          end if
       end do
       if(pub_on_root) then
          write(ss_out_unit,*)
          if(MOD(irow,num_nto_trans)==0) write(ss_out_unit,*)
       end if
    end do
    if(pub_on_root) then
       do irow=1, 1
         do i=1, 13+15*num_states
           write(ss_out_unit,'(a)',advance='no') '='
         end do
         write(ss_out_unit,*)
       end do
       close(unit=ss_out_unit,iostat=ierr)
       call utils_close_unit_check('qnto/nto_output_proj_coeff',trim(ss_out_filename),ierr)
    end if


    !------------------------------------------------------------------------------------!
    ! Output NTOn-H(E)_Sys to cNTOn-H(E)_Sys projection results                          !
    !------------------------------------------------------------------------------------!
    if(pub_on_root) then
       write(ss_out_filename,'(a)') trim(pub_rootname)//'_QNTO_'//&
            trim(adjustl(name_num_core))//'a'//trim(adjustl(name_n))//'t_s2s.dat'
       ss_out_filename = adjustl(ss_out_filename)
       ss_out_unit = utils_unit()
       open(unit=ss_out_unit,form="formatted",status="old",position="append", &
            file=trim(ss_out_filename),action="write",iostat=ierr)
       call utils_open_unit_check('qnto/nto_output_proj_coeff',trim(ss_out_filename),ierr)

       ! Write header
       write(ss_out_unit,'(a,//a,/a,/a,i4,a,//a,/a,/a,/a/)') &
            'S2Sc:', & !'Projection of NTOn-H(E)_S to cNTOn-H(E)_S:'
            &'   "<cNTO_S|NTO_S> / (<cNTO_S|cNTO_S>)^0.5"', &
            '   using <NGWF_S|NGWF_S> on S-geom', &
            '   (qnto_num_core_atoms = ', pub_qnto_num_core_atoms, ')', &
            'Map format:', &
            '             NTO_S', &
            '   cNTO_S| H2H; E2H', &
            '         | H2E; E2E'

       write(ss_out_unit,'(12x)',advance='no')
       do state_idx_sys=1, num_states
          do nto_idx_sys=1, num_nto_trans
             write(ss_out_unit,'(1x,i1,a,i3.3,1x,i1,a,i3.3)',advance='no') &
                  nto_idx_sys,'H_',state_idx_sys, nto_idx_sys,'E_',state_idx_sys
          end do
          write(ss_out_unit,'(a)',advance='no') '|'
       end do
       write(ss_out_unit,*)

       ! Output transition coeff of NTOs_sys
       write(ss_out_unit,'(11x,a)',advance='no') '|'
       do state_idx_sys=1, num_states
          do nto_idx_sys=1, num_nto_trans
             n = (state_idx_sys-1)*num_nto_trans + nto_idx_sys
             write(ss_out_unit,'(6x,i3.2,5x)',advance='no') &
                  nint(trans_coeff_sys(n)*100.0_DP)
          end do
          write(ss_out_unit,'(a)',advance='no') '|'
       end do
       write(ss_out_unit,*)
    end if !pub_on_root

    ! Plot the map
    allocate(e_s2sc(systgt_size), stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','e_s2sc', ierr)
    allocate(nor_h_s2sc(systgt_size), stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','nor_h_s2sc', ierr)
    allocate(nor_e_s2sc(systgt_size), stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','nor_e_s2sc', ierr)
    allocate(label_flag_sys(num_states), stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','label_flag_sys', ierr)
    allocate(I_match(systgt_size),stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','I_match', ierr)
    if(pub_on_root) then
       do state_idx_sys=1, num_states
          label_flag_sys(state_idx_sys) = .FALSE.
       end do
    end if
    do irow=1, systgt_size
       if(pub_on_root) then
          write(ss_out_unit,'(1x,i1,a,i3.3,1x,i3.2,a)',advance='no') &
               irow-(irow-1)/num_nto_trans*num_nto_trans,'H_', &
               (irow-1)/num_nto_trans+1,&
               nint(trans_coeff_sys(irow)*100.0_DP),'|'
          num_match = 0
          I_match(:) = 0
       end if

       do icol=1, systgt_size
          call dense_get_element(mtx_el,nto_sh2sch,irow,icol)
          if(pub_on_root) then
             h_s2sc = mtx_el / nor_h_sc(irow)
             if(icol==irow) nor_h_s2sc(icol) = h_s2sc
          end if
          call dense_get_element(mtx_el,nto_se2sce,irow,icol)
          if(pub_on_root) then
             e_s2sc(icol) = mtx_el / nor_e_sc(irow)
             if(icol==irow) nor_e_s2sc(icol) = e_s2sc(icol)

             if(abs(h_s2sc) >= sqrt(1.0_DP/2.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '*', nint(h_s2sc*100.0_DP)
             else if(abs(h_s2sc) >= sqrt(1.0_DP/3.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '#', nint(h_s2sc*100.0_DP)
             else if(h_s2sc>-0.005_DP .AND. h_s2sc<0.0_DP) then
                write(ss_out_unit,'(3x,a)',advance='no') ' -00'
             else
                write(ss_out_unit,'(3x,i4.2)',advance='no') nint(h_s2sc*100.0_DP)
             end if
          end if !pub_on_root

          call dense_get_element(mtx_el,nto_se2sch,irow,icol)
          if(pub_on_root) then
             mtx_el = mtx_el / nor_h_sc(irow)
             if(abs(mtx_el) >= sqrt(1.0_DP/2.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '*', nint(mtx_el*100.0_DP)
             else if(abs(mtx_el) >= sqrt(1.0_DP/3.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '#', nint(mtx_el*100.0_DP)
             else if(mtx_el>-0.005_DP .AND. mtx_el<0.0_DP) then
                write(ss_out_unit,'(3x,a)',advance='no') ' -00'
             else
                write(ss_out_unit,'(3x,i4.2)',advance='no') nint(mtx_el*100.0_DP)
             end if
             if(MOD(icol,num_nto_trans)==0) &
               write(ss_out_unit,'(a)',advance='no') '|'
          end if !pub_on_root
       end do

       if(pub_on_root) then
          nto_idx_sys = MOD(irow+num_nto_trans-1, num_nto_trans) + 1
          state_idx_sys = (irow+num_nto_trans-1) / num_nto_trans
          write(ss_out_unit,*)
          write(ss_out_unit,'(1x,i1,a,i3.3,4x,a)',advance='no') &
               irow-(irow-1)/num_nto_trans*num_nto_trans,'E_', &
               (irow-1)/num_nto_trans+1,'|'
       end if !pub_on_root

       do icol=1, systgt_size
          call dense_get_element(mtx_el,nto_sh2sce,irow,icol)
          if(pub_on_root) then
             mtx_el = mtx_el / nor_e_sc(irow)
             if(abs(mtx_el) >= sqrt(1.0_DP/2.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '*', nint(mtx_el*100.0_DP)
             else if(abs(mtx_el) >= sqrt(1.0_DP/3.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '#', nint(mtx_el*100.0_DP)
             else if(mtx_el>-0.005_DP .AND. mtx_el<0.0_DP) then
                write(ss_out_unit,'(3x,a)',advance='no') ' -00'
             else
                write(ss_out_unit,'(3x,i4.2)',advance='no') nint(mtx_el*100.0_DP)
             end if
             if(abs(e_s2sc(icol)) >= sqrt(1.0_DP/2.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '*', &
                     nint(e_s2sc(icol)*100.0_DP)
             else if(abs(e_s2sc(icol)) >= sqrt(1.0_DP/3.0_DP)) then
                write(ss_out_unit,'(2x,a,i4.2)',advance='no') '#', &
                     nint(e_s2sc(icol)*100.0_DP)
             else if(e_s2sc(icol)>-0.005_DP .AND. e_s2sc(icol)<0.0_DP) then
                write(ss_out_unit,'(3x,a)',advance='no') ' -00'
             else
                write(ss_out_unit,'(3x,i4.2)',advance='no') nint(e_s2sc(icol)*100.0_DP)
             end if
             if(MOD(icol,num_nto_trans)==0) &
               write(ss_out_unit,'(a)',advance='no') '|'
          end if !pub_on_root
       end do
       if(pub_on_root) then
          write(ss_out_unit,*)
          if(MOD(irow,num_nto_trans)==0) write(ss_out_unit,*)
       end if
    end do
    if(pub_on_root) then
       do irow=1, 1
         do i=1, 13+15*num_states
           write(ss_out_unit,'(a)',advance='no') '='
         end do
         write(ss_out_unit,*)
       end do
       close(unit=ss_out_unit,iostat=ierr)
       call utils_close_unit_check('qnto/nto_output_proj_coeff',trim(ss_out_filename),ierr)
    end if


    if(pub_on_root) then
       write(ss2_out_filename,'(a)') trim(pub_rootname)//'_QNTO_'//&
            trim(adjustl(name_num_core))//'a'//trim(adjustl(name_n))//'t_s2s.dat'
       ss2_out_filename = adjustl(ss2_out_filename)
       ss2_out_unit = utils_unit()
       open(unit=ss2_out_unit,form="formatted",status="old",position="append", &
            file=trim(ss2_out_filename),action="write",iostat=ierr)
       call utils_open_unit_check('qnto/nto_output_proj_coeff',trim(ss2_out_filename),ierr)

       ! Write header
       write(ss2_out_unit,'(a,//a,/a,/a,i4,a,//a,/a,/a/)') &
            'S2Sc_H2H_E2E:', & !'(Projection of NTOn-H(E)_S to cNTOn-H(E)_S):'
            &'   "<cNTO_S|NTO_S> / (<cNTO_S|cNTO_S>)^0.5"', &
            '   using <NGWF_S|NGWF_S> on S-geom', &
            '   (qnto_num_core_atoms = ', pub_qnto_num_core_atoms, ')', &
            'Map format:', &
            '             NTO_S', &
            '   cNTO_S| H2H; E2E'

       write(ss2_out_unit,'(13x)',advance='no')
       do state_idx_sys=1, num_states
          do nto_idx_sys=1, num_nto_trans
             write(ss2_out_unit,'(3x,i1,a,i3.3)',advance='no') &
                  nto_idx_sys,'HE_',state_idx_sys
          end do
          write(ss2_out_unit,'(a)',advance='no') '|'
       end do
       write(ss2_out_unit,'(3x,a)') 'Matchings'

       ! Output transition coeff of NTOs_sys
       write(ss2_out_unit,'(12x,a)',advance='no') '|'
       do state_idx_sys=1, num_states
          do nto_idx_sys=1, num_nto_trans
             n = (state_idx_sys-1)*num_nto_trans + nto_idx_sys
             write(ss2_out_unit,'(5x,i3.2,2x)',advance='no') &
                  nint(trans_coeff_sys(n)*100.0_DP)
          end do
          write(ss2_out_unit,'(a)',advance='no') '|'
       end do
       write(ss2_out_unit,*)
    end if
    do irow=1, systgt_size
       if(pub_on_root) then
          write(ss2_out_unit,'(1x,i1,a,i3.3,1x,i3.2,a)',advance='no') &
               irow-(irow-1)/num_nto_trans*num_nto_trans,'HE_', &
               (irow-1)/num_nto_trans+1, &
               nint(trans_coeff_sys(irow)*100.0_DP),'|'
          num_match = 0
          I_match(:) = 0
       end if

       do icol=1, systgt_size
          call dense_get_element(mtx_el,nto_sh2sch,irow,icol)
          if(pub_on_root) then
             h_s2sc = mtx_el / nor_h_sc(irow)
             if(icol==irow) nor_h_s2sc(icol) = h_s2sc
          end if
          call dense_get_element(mtx_el,nto_se2sce,irow,icol)
          if(pub_on_root) then
             e_s2sc(icol) = mtx_el / nor_e_sc(irow)
             if(icol==irow) nor_e_s2sc(icol) = e_s2sc(icol)

             if(abs(h_s2sc)>=sqrt(1.0_DP/2.0_DP) .AND. &
                abs(e_s2sc(icol))>=sqrt(1.0_DP/2.0_DP) ) then
                write(ss2_out_unit,'(a1)',advance='no') '*'
                num_match = num_match + 1
                I_match(num_match) = icol
             else if(abs(h_s2sc)>=sqrt(1.0_DP/3.0_DP) .AND. &
                abs(e_s2sc(icol))>=sqrt(1.0_DP/3.0_DP) ) then
                write(ss2_out_unit,'(a1)',advance='no') '#'
                num_match = num_match + 1
                I_match(num_match) = icol-systgt_size
             else
                write(ss2_out_unit,'(a1)',advance='no') ' '
             end if

             if(h_s2sc>-0.005_DP .AND. h_s2sc<0.0_DP) then
                write(ss2_out_unit,'(a)',advance='no') ' -00'
             else
                write(ss2_out_unit,'(i4.2)',advance='no') nint(h_s2sc*100.0_DP)
             end if
             write(ss2_out_unit,'(a1)',advance='no') ';'

             if(e_s2sc(icol)>-0.005_DP .AND. e_s2sc(icol)<0.0_DP) then
                write(ss2_out_unit,'(a)',advance='no') ' -00'
             else
                write(ss2_out_unit,'(i4.2)',advance='no') nint(e_s2sc(icol)*100.0_DP)
             end if

             if(MOD(icol,num_nto_trans)==0) &
               write(ss2_out_unit,'(a)',advance='no') '|'
          end if !pub_on_root
       end do

       if(pub_on_root) then
          nto_idx_sys = MOD(irow+num_nto_trans-1, num_nto_trans) + 1
          state_idx_sys = (irow+num_nto_trans-1) / num_nto_trans
          write(ss2_out_unit,'(2x,i3,1x,a,i1,a,i3.3,a)',advance='no') &
               num_match, 'Sc', nto_idx_sys, '_', state_idx_sys, ' - S'
          do i=1, num_match
             if(I_match(i) > 0) then
                nto_idx2_sys = MOD(I_match(i)+num_nto_trans-1, &
                               num_nto_trans) + 1
                state_idx2_sys = (I_match(i)+num_nto_trans-1) / &
                                 num_nto_trans
                if(I_match(i) == irow) then
                   write(ss2_out_unit,'(a,i1,a,i3.3)',advance='no') &
                        '**', nto_idx2_sys, '_', state_idx2_sys
                   if(nto_idx2_sys == 1) label_flag_sys(state_idx_sys) = .TRUE.
                else
                   write(ss2_out_unit,'(a,i1,a,i3.3)',advance='no') &
                        ' *', nto_idx2_sys, '_', state_idx2_sys
                end if
             else
                nto_idx2_sys = MOD(I_match(i)+systgt_size+num_nto_trans-1, &
                               num_nto_trans) + 1
                state_idx2_sys = (I_match(i)+systgt_size+num_nto_trans-1) / &
                                 num_nto_trans
                if(I_match(i)+systgt_size == irow) then
                   write(ss2_out_unit,'(a,i1,a,i3.3)',advance='no') &
                        '##', nto_idx2_sys, '_', state_idx2_sys
                else
                   write(ss2_out_unit,'(a,i1,a,i3.3)',advance='no') &
                        ' #', nto_idx2_sys, '_', state_idx2_sys
                end if
             end if
          end do
          write(ss2_out_unit,*)
       end if !pub_on_root
    end do
    deallocate(e_s2sc,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','e_s2sc', ierr)
    if(pub_on_root) then
       close(unit=ss2_out_unit,iostat=ierr)
       call utils_close_unit_check('qnto/nto_output_proj_coeff',trim(ss2_out_filename),ierr)
    end if


    ! Construct normalisation factors for hole and electron projection results
    allocate(nor_h_rc(reftgt_size),stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','nor_h_rc', ierr)
    allocate(nor_e_rc(reftgt_size),stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','nor_e_rc', ierr)
    !------------------------------------------------------------------------------------!
    ! Output cNTOn-H(E)_Ref to cNTOn-H(E)_Ref projection results                         !
    !------------------------------------------------------------------------------------!
    if(pub_on_root) then
       write(rr_io_filename,'(a)') trim(pub_rootname)//'_QNTO_'//&
            trim(adjustl(name_num_core))//'a'//trim(adjustl(name_n))//'t_s2r_'//trim(ref_dir)//'.dat'
       rr_io_filename = adjustl(rr_io_filename)
       rr_io_unit = utils_unit()
       open(unit=rr_io_unit,form="formatted", &
            file=trim(rr_io_filename),action="write",iostat=ierr)
       call utils_open_unit_check('qnto/nto_output_proj_coeff',trim(rr_io_filename),ierr)
       ! Write header
       write(rr_io_unit,'(a,//a,/a,/a,i4,a,//a,/a,/a,/a/)') &
            'Rc2Rc:', & !'Projection of cNTOn-H(E)a_R to cNTOn-H(E)b_R:'
            &'   "<cNTOb_R|cNTOa_R>"', &
            '   using <NGWF_R|NGWF_R> on S-geom', &
            '   (qnto_num_core_atoms = ', pub_qnto_num_core_atoms, ')', &
            'Map format:', &
            '            cNTOa_R', &
            '   cNTOb_R| H2H; E2H', &
            '          | H2E; E2E'

       write(rr_io_unit,'(12x)',advance='no')
       do state_idx_ref=1, pub_qnto_num_ref_states
          do nto_idx_ref=1, num_nto_trans
             write(rr_io_unit,'(1x,i1,a,i3.3,1x,i1,a,i3.3)',advance='no') &
                  nto_idx_ref,'H_',state_idx_ref, nto_idx_ref,'E_',state_idx_ref
          end do
          write(rr_io_unit,'(a)',advance='no') '|'
       end do
       write(rr_io_unit,*)

       ! Output transition coeff of NTOs_ref
       write(rr_io_unit,'(11x,a)',advance='no') '|'
       do state_idx_ref=1, pub_qnto_num_ref_states
          do nto_idx_ref=1, num_nto_trans
             n = (state_idx_ref-1)*num_nto_trans + nto_idx_ref
             write(rr_io_unit,'(6x,i3.2,5x)',advance='no') &
                  nint(trans_coeff_ref(n)*100.0_DP)
          end do
          write(rr_io_unit,'(a)',advance='no') '|'
       end do
       write(rr_io_unit,*)
    end if !pub_on_root

    ! Plot the map
    do irow=1, reftgt_size
       if(pub_on_root) write(rr_io_unit,'(1x,i1,a,i3.3,1x,i3.2,a)',advance='no') &
       irow-(irow-1)/num_nto_trans*num_nto_trans,'H_', &
       (irow-1)/num_nto_trans+1,nint(trans_coeff_ref(irow)*100.0_DP),'|'
       do icol=1, reftgt_size
          call dense_get_element(mtx_el,nto_rch2rch,irow,icol)
          if(pub_on_root) then
             ! Construct nor_h_rc
             if(icol == irow) nor_h_rc(irow) = sqrt(mtx_el)
             if(abs(mtx_el) >= sqrt(1.0_DP/2.0_DP)) then
                write(rr_io_unit,'(2x,a,i4.2)',advance='no') '*',nint(mtx_el*100.0_DP)
             else if(abs(mtx_el) >= sqrt(1.0_DP/3.0_DP)) then
                write(rr_io_unit,'(2x,a,i4.2)',advance='no') '#',nint(mtx_el*100.0_DP)
             else if(mtx_el>-0.005_DP .AND. mtx_el<0.0_DP) then
                write(rr_io_unit,'(3x,a)',advance='no') ' -00'
             else
                write(rr_io_unit,'(3x,i4.2)',advance='no') nint(mtx_el*100.0_DP)
             end if
          end if

          call dense_get_element(mtx_el,nto_rce2rch,irow,icol)
          if(pub_on_root) then
             if(abs(mtx_el) >= sqrt(1.0_DP/2.0_DP)) then
                write(rr_io_unit,'(2x,a,i4.2)',advance='no') '*',nint(mtx_el*100.0_DP)
             else if(abs(mtx_el) >= sqrt(1.0_DP/3.0_DP)) then
                write(rr_io_unit,'(2x,a,i4.2)',advance='no') '#',nint(mtx_el*100.0_DP)
             else if(mtx_el>-0.005_DP .AND. mtx_el<0.0_DP) then
                write(rr_io_unit,'(3x,a)',advance='no') ' -00'
             else
                write(rr_io_unit,'(3x,i4.2)',advance='no') nint(mtx_el*100.0_DP)
             end if
             if(MOD(icol,num_nto_trans)==0) &
               write(rr_io_unit,'(a)',advance='no') '|'
          end if
       end do

       if(pub_on_root) write(rr_io_unit,*)
       if(pub_on_root) write(rr_io_unit,'(1x,i1,a,i3.3,4x,a)',advance='no') &
       irow-(irow-1)/num_nto_trans*num_nto_trans,'E_', &
       (irow-1)/num_nto_trans+1,'|'
       do icol=1, reftgt_size
          call dense_get_element(mtx_el,nto_rch2rce,irow,icol)
          if(pub_on_root) then
             if(abs(mtx_el) >= sqrt(1.0_DP/2.0_DP)) then
                write(rr_io_unit,'(2x,a,i4.2)',advance='no') '*',nint(mtx_el*100.0_DP)
             else if(abs(mtx_el) >= sqrt(1.0_DP/3.0_DP)) then
                write(rr_io_unit,'(2x,a,i4.2)',advance='no') '#',nint(mtx_el*100.0_DP)
             else if(mtx_el>-0.005_DP .AND. mtx_el<0.0_DP) then
                write(rr_io_unit,'(3x,a)',advance='no') ' -00'
             else
                write(rr_io_unit,'(3x,i4.2)',advance='no') nint(mtx_el*100.0_DP)
             end if
          end if

          call dense_get_element(mtx_el,nto_rce2rce,irow,icol)
          if(pub_on_root) then
             ! Construct nor_e_rc
             if(icol == irow) nor_e_rc(irow) = sqrt(mtx_el)
             if(abs(mtx_el) >= sqrt(1.0_DP/2.0_DP)) then
                write(rr_io_unit,'(2x,a,i4.2)',advance='no') '*',nint(mtx_el*100.0_DP)
             else if(abs(mtx_el) >= sqrt(1.0_DP/3.0_DP)) then
                write(rr_io_unit,'(2x,a,i4.2)',advance='no') '#',nint(mtx_el*100.0_DP)
             else if(mtx_el>-0.005_DP .AND. mtx_el<0.0_DP) then
                write(rr_io_unit,'(3x,a)',advance='no') ' -00'
             else
                write(rr_io_unit,'(3x,i4.2)',advance='no') nint(mtx_el*100.0_DP)
             end if
             if(MOD(icol,num_nto_trans)==0) &
               write(rr_io_unit,'(a)',advance='no') '|'
          end if
       end do
       if(pub_on_root) then
          write(rr_io_unit,*)
          if(MOD(irow,num_nto_trans)==0) write(rr_io_unit,*)
       end if
    end do
    if(pub_on_root) then
       do irow=1, 1
         do i=1, 13+15*num_states
           write(rr_io_unit,'(a)',advance='no') '='
         end do
         write(rr_io_unit,*)
       end do
       close(unit=rr_io_unit,iostat=ierr)
       call utils_close_unit_check('qnto/nto_output_proj_coeff',trim(rr_io_filename),ierr)
    end if


    !------------------------------------------------------------------------------------!
    ! Initialise label_sys(:), label_ref(:), label_flag_ref(:), nor_h_r2rc(:),           !
    ! and nor_e_r2rc(:)                                                                  !
    !------------------------------------------------------------------------------------!
    allocate(label_flag_ref(pub_qnto_num_ref_states), stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','label_flag_ref', ierr)
    allocate(nor_h_r2rc(reftgt_size), stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','nor_h_r2rc', ierr)
    allocate(nor_e_r2rc(reftgt_size), stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','nor_e_r2rc', ierr)
    allocate(label_sys(num_states), stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','label_sys', ierr)
    allocate(tmp_string(num_states), stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','tmp_string', ierr)
    allocate(hash_flag_sys(num_states), stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','hash_flag_sys', ierr)
    allocate(label_ref(pub_qnto_num_ref_states), stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','label_ref', ierr)
    allocate(hash_flag_ref(pub_qnto_num_ref_states), stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','hash_flag_ref', ierr)

    if(pub_on_root) then
       ! Set up label_flag_ref
       write(rr_io_filename,'(a)') '../'//trim(ref_dir)//'/'//&
            trim(pub_rootname)//'_QNTO_'//trim(adjustl(name_num_core))//&
            'a'//trim(adjustl(name_n))//'t_s2s.dat'
       rr_io_filename = adjustl(rr_io_filename)
       rr_io_unit = utils_unit()
       inquire(FILE=trim(rr_io_filename), EXIST=f_exist)
       if(.not. f_exist) then
          call utils_abort('ERROR in qnto_calculate: '//trim(rr_io_filename)//' does not exist.')

       end if
       open(unit=rr_io_unit,form="formatted", &
            file=trim(rr_io_filename),action="read",iostat=ierr)
       call utils_open_unit_check('qnto/nto_output_proj_coeff',trim(rr_io_filename),ierr)
       do
          read(rr_io_unit,'(a)') title
          if( trim(adjustl(title)) == 'S2Sc_H2H_E2E:') exit
       end do
       do i=1,11
          read(rr_io_unit,*)
       end do
       do state_idx_ref=1, pub_qnto_num_ref_states
          do nto_idx_ref=1, num_nto_trans
          ! The referenced file has to be generated with the same num_nto_trans
             read(rr_io_unit,'(13x)',advance='no')
             do i=1, (state_idx_ref-1)*(10*num_nto_trans+1)+10*(nto_idx_ref-1)
                read(rr_io_unit,'(1x)',advance='no')
             end do
             if(nto_idx_ref==1) then
                read(rr_io_unit,'(a1)',advance='no') mark
                if(mark == '*') then
                   label_flag_ref(state_idx_ref) = .TRUE.
                else
                   label_flag_ref(state_idx_ref) = .FALSE.
                end if
                read(rr_io_unit,'(a4)',advance='no') str
                read(str,'(i4)') nor_h_r2rc_int
                nor_h_r2rc((state_idx_ref-1)*num_nto_trans+nto_idx_ref) = &
                nor_h_r2rc_int / 100.0_DP
                read(rr_io_unit,'(1x,a4)') str
                read(str,'(i4)') nor_e_r2rc_int
                nor_e_r2rc((state_idx_ref-1)*num_nto_trans+nto_idx_ref) = &
                nor_e_r2rc_int / 100.0_DP
             else
                read(rr_io_unit,'(1x,i4)',advance='no') nor_h_r2rc_int
                nor_h_r2rc((state_idx_ref-1)*num_nto_trans+nto_idx_ref) = &
                nor_h_r2rc_int / 100.0_DP
                read(rr_io_unit,'(1x,i4)') nor_e_r2rc_int
                nor_e_r2rc((state_idx_ref-1)*num_nto_trans+nto_idx_ref) = &
                nor_e_r2rc_int / 100.0_DP
             end if
             if (pub_debug_on_root) then
                write(stdout,*)'nor_he_r2rc:',(state_idx_ref-1)* &
                     num_nto_trans+nto_idx_ref, &
                     nor_h_r2rc((state_idx_ref-1)*num_nto_trans+nto_idx_ref), &
                     nor_e_r2rc((state_idx_ref-1)*num_nto_trans+nto_idx_ref)
             end if
          end do
       end do
       close(unit=rr_io_unit,iostat=ierr)
       call utils_close_unit_check('qnto/nto_output_proj_coeff',trim(rr_io_filename),ierr)

       ! Read in label_ref
       write(label_io_filename,'(a)') '../'//trim(ref_dir)//'/'//trim(pub_rootname)//&
            '_state_labellings.dat'
       label_io_filename = adjustl(label_io_filename)
       label_io_unit = utils_unit()
       inquire(FILE=trim(label_io_filename), EXIST=f_exist)
       if(.not. f_exist) then
          write(stdout,*) 'Warning:', trim(label_io_filename), ' does not exist.'
          do state_idx_ref=1, pub_qnto_num_ref_states
             label_ref(state_idx_ref) = '  '
             hash_flag_ref(state_idx_ref) = ' '
          end do
       else
          open(unit=label_io_unit,form="formatted", file=trim(label_io_filename),&
               action="read", iostat=ierr)
          call utils_open_unit_check('qnto/nto_output_proj_coeff',&
               trim(label_io_filename),ierr)
          do
             read(label_io_unit,'(a1)',advance='no') mark
             if(mark == '-') then
                read(label_io_unit,'(3x)',advance='no')
                do state_idx_ref=1, pub_qnto_num_ref_states
                   read(label_io_unit,'(a,a,3x)',advance='no') label_ref(state_idx_ref), &
                       hash_flag_ref(state_idx_ref)
                end do
                exit
             else
                read(label_io_unit,*)
             end if
          end do
          close(unit=label_io_unit,iostat=ierr)
          call utils_close_unit_check('qnto/nto_output_proj_coeff',&
               trim(label_io_filename),ierr)
       end if !f_exist

       ! Read in label_sys
       write(label_io_filename,'(a)') trim(pub_rootname)//'_state_labellings.dat.tmp'
       label_io_filename = adjustl(label_io_filename)
       label_io_unit = utils_unit()

       inquire(file=label_io_filename, exist=f_exist)
       if(.not. f_exist) then
          do state_idx_sys=1, num_states
             label_sys(state_idx_sys) = '  '
          end do
       else
          open(unit=label_io_unit,form="formatted", file=trim(label_io_filename),&
               action="readwrite", iostat=ierr)
          call utils_open_unit_check('qnto/nto_output_proj_coeff',&
               trim(label_io_filename),ierr)
          n = 0
          do
             read(label_io_unit,'(a)',advance='no') mark
             if(mark == '-') then
                read(label_io_unit,'(3x)',advance='no')
                do state_idx_sys=1, num_states
                   read(label_io_unit,'(a,a,a,1x)',advance='no') label_sys(state_idx_sys), &
                       hash_flag_sys(state_idx_sys), tmp_string(state_idx_sys)
                end do
                rewind(label_io_unit)
                do i=1, n
                   read(label_io_unit,*)
                end do
                write(label_io_unit,'(a,2x,a)',advance='no') '%', '|'
                do state_idx_sys=1, num_states
                   write(label_io_unit,'(a,a,a,a)',advance='no') label_sys(state_idx_sys), &
                        hash_flag_sys(state_idx_sys), tmp_string(state_idx_sys), '|'
                   label_sys(state_idx_sys) = '  '
                end do
                write(label_io_unit,*)
                write(label_io_unit,'(a)',advance='no') '%'
                exit
             else
                read(label_io_unit,*)
                n = n + 1
             end if
          end do
       end if !f_exist
       close(unit=label_io_unit,iostat=ierr)
       call utils_close_unit_check('qnto/nto_output_proj_coeff',trim(label_io_filename),ierr)
    end if !pub_on_root
    deallocate(tmp_string,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','tmp_string', ierr)


    !------------------------------------------------------------------------------------!
    ! Output cNTOn-H(E)_Sys to cNTOn-H(E)_Ref projection results                         !
    !------------------------------------------------------------------------------------!
    if(pub_on_root) then
       write(sr_out_filename,'(a)') trim(pub_rootname)//'_QNTO_'//&
            trim(adjustl(name_num_core))//'a'//trim(adjustl(name_n))//'t_s2r_'//trim(ref_dir)//'.dat'
       sr_out_filename = adjustl(sr_out_filename)
       sr_out_unit = utils_unit()
       open(unit=sr_out_unit,form="formatted",status="old",position="append", &
            file=trim(sr_out_filename),action="write",iostat=ierr)
       call utils_open_unit_check('qnto/nto_output_proj_coeff',trim(sr_out_filename),ierr)

       ! Write header
       write(sr_out_unit,'(a,//a,/a,/a,i4,a,//a,/a,/a,/a/)') &
            'Sc2Rc:', & !Projection of cNTOn-H(E)_S to cNTOn-H(E)_R:'
            &'   "<cNTO_R|cNTO_S> / (<cNTO_R|cNTO_R>)^0.5 * <cNTO_S|cNTO_S>^0.5)"', &
            '   using <NGWF_R|NGWF_R>, <NGWF_S|NGWF_S>, and <NGWF_R|NGWF_S> on S-geom', &
            '   (qnto_num_core_atoms = ', pub_qnto_num_core_atoms, ')', &
            'Map format:', &
            '          cNTO_S', &
            '   cNTO_R| H2H; E2H', &
            '         | H2E; E2E'

       write(sr_out_unit,'(12x)',advance='no')
       do state_idx_sys=1, num_states
          do nto_idx_sys=1, num_nto_trans
             write(sr_out_unit,'(1x,i1,a,i3.3,1x,i1,a,i3.3)',advance='no') &
                  nto_idx_sys,'H_',state_idx_sys, nto_idx_sys,'E_',state_idx_sys
          end do
          write(sr_out_unit,'(a)',advance='no') '|'
       end do
       write(sr_out_unit,*)

       ! Output transition coeff of NTOs_sys
       write(sr_out_unit,'(11x,a)',advance='no') '|'
       do state_idx_sys=1, num_states
          do nto_idx_sys=1, num_nto_trans
             n = (state_idx_sys-1)*num_nto_trans + nto_idx_sys
             write(sr_out_unit,'(6x,i3.2,5x)',advance='no') &
                  nint(trans_coeff_sys(n)*100.0_DP)
          end do
          write(sr_out_unit,'(a)',advance='no') '|'
       end do
       write(sr_out_unit,*)
    end if !pub_on_root


    allocate(h_sc2rc(reftgt_size,systgt_size),stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','h_sc2rc', ierr)
    allocate(e_sc2rc(reftgt_size,systgt_size),stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','e_sc2rc', ierr)
    allocate(match_weight(num_states),stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','match_weight', ierr)
    allocate(match_weight_sys(num_states),stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','match_weight_sys', ierr)
    allocate(match_energy_diff(num_states),stat=ierr)
    call utils_alloc_check('qnto/nto_output_proj_coeff','match_energy_diff', ierr)
    match_weight(:) = 0.0_DP
    match_weight_sys(:) = 0.0_DP

    ! Plot the map
    do irow=1, reftgt_size
       if(pub_on_root) then
          write(sr_out_unit,'(1x,i1,a,i3.3,1x,i3.2,a)',advance='no') &
               irow-(irow-1)/num_nto_trans*num_nto_trans,'H_', &
               (irow-1)/num_nto_trans+1,&
               nint(trans_coeff_ref(irow)*100.0_DP),'|'
          nto_idx_ref = MOD(irow+num_nto_trans-1, num_nto_trans) + 1
          state_idx_ref = (irow+num_nto_trans-1) / num_nto_trans
       end if !pub_on_root

       ! 1st icol run
       do icol=1, systgt_size
          call dense_get_element(mtx_el,nto_sch2rch,irow,icol)
          if(pub_on_root) h_sc2rc(irow,icol) = mtx_el / (nor_h_rc(irow)*nor_h_sc(icol))
          call dense_get_element(mtx_el,nto_sce2rce,irow,icol)
          if(pub_on_root) then
             e_sc2rc(irow,icol) = mtx_el / (nor_e_rc(irow)*nor_e_sc(icol))

             if(abs(h_sc2rc(irow,icol)) >= sqrt(1.0_DP/2.0_DP)) then
                write(sr_out_unit,'(2x,a,i4.2)',advance='no') '*', nint(h_sc2rc(irow,icol)*100.0_DP)
             else if(abs(h_sc2rc(irow,icol)) >= sqrt(1.0_DP/3.0_DP)) then
                write(sr_out_unit,'(2x,a,i4.2)',advance='no') '#', nint(h_sc2rc(irow,icol)*100.0_DP)
             else if(h_sc2rc(irow,icol)>-0.005_DP .AND. h_sc2rc(irow,icol)<0.0_DP) then
                write(sr_out_unit,'(3x,a)',advance='no') ' -00'
             else
                write(sr_out_unit,'(3x,i4.2)',advance='no') nint(h_sc2rc(irow,icol)*100.0_DP)
             end if
          end if !pub_on_root

          call dense_get_element(mtx_el,nto_sce2rch,irow,icol)
          if(pub_on_root) then
             mtx_el = mtx_el / (nor_h_rc(irow)*nor_e_sc(icol))
             if(abs(mtx_el) >= sqrt(1.0_DP/2.0_DP)) then
                write(sr_out_unit,'(2x,a,i4.2)',advance='no') '*', nint(mtx_el*100.0_DP)
             else if(abs(mtx_el) >= sqrt(1.0_DP/3.0_DP)) then
                write(sr_out_unit,'(2x,a,i4.2)',advance='no') '#', nint(mtx_el*100.0_DP)
             else if(mtx_el>-0.005_DP .AND. mtx_el<0.0_DP) then
                write(sr_out_unit,'(3x,a)',advance='no') ' -00'
             else
                write(sr_out_unit,'(3x,i4.2)',advance='no') nint(mtx_el*100.0_DP)
             end if
             if(MOD(icol,num_nto_trans)==0) &
             write(sr_out_unit,'(a)',advance='no') '|'
          end if !pub_on_root
       end do
       ! 1st icol run finished

       if(pub_on_root) then
          write(sr_out_unit,*)
          write(sr_out_unit,'(1x,i1,a,i3.3,4x,a)',advance='no') &
               irow-(irow-1)/num_nto_trans*num_nto_trans,'E_', &
               (irow-1)/num_nto_trans+1,'|'
       end if !pub_on_root

       ! 2nd icol run
       do icol=1, systgt_size
          call dense_get_element(mtx_el,nto_sch2rce,irow,icol)
          if(pub_on_root) then
             mtx_el = mtx_el / (nor_e_rc(irow)*nor_h_sc(icol))
             if(abs(mtx_el) >= sqrt(1.0_DP/2.0_DP)) then
                write(sr_out_unit,'(2x,a,i4.2)',advance='no') '*', nint(mtx_el*100.0_DP)
             else if(abs(mtx_el) >= sqrt(1.0_DP/3.0_DP)) then
                write(sr_out_unit,'(2x,a,i4.2)',advance='no') '#', nint(mtx_el*100.0_DP)
             else if(mtx_el>-0.005_DP .AND. mtx_el<0.0_DP) then
                write(sr_out_unit,'(3x,a)',advance='no') ' -00'
             else
                write(sr_out_unit,'(3x,i4.2)',advance='no') nint(mtx_el*100.0_DP)
             end if

             if(abs(e_sc2rc(irow,icol)) >= sqrt(1.0_DP/2.0_DP)) then
                write(sr_out_unit,'(2x,a,i4.2)',advance='no') '*', &
                     nint(e_sc2rc(irow,icol)*100.0_DP)
             else if(abs(e_sc2rc(irow,icol)) >= sqrt(1.0_DP/3.0_DP)) then
                write(sr_out_unit,'(2x,a,i4.2)',advance='no') '#', &
                     nint(e_sc2rc(irow,icol)*100.0_DP)
             else if(e_sc2rc(irow,icol)>-0.005_DP .AND. e_sc2rc(irow,icol)<0.0_DP) then
                write(sr_out_unit,'(3x,a)',advance='no') ' -00'
             else
                write(sr_out_unit,'(3x,i4.2)',advance='no') nint(e_sc2rc(irow,icol)*100.0_DP)
             end if
             if(MOD(icol,num_nto_trans)==0) &
               write(sr_out_unit,'(a)',advance='no') '|'
          end if !pub_on_root
       end do
       ! 2nd icol run finished
       if(pub_on_root) then
          write(sr_out_unit,*)
          if(MOD(irow,num_nto_trans)==0) write(sr_out_unit,*)
       end if
    end do !irow
    if(pub_on_root) then
       do irow=1, 1
         do i=1, 13+15*num_states
           write(sr_out_unit,'(a)',advance='no') '='
         end do
         write(sr_out_unit,*)
       end do
       close(unit=sr_out_unit,iostat=ierr)
       call utils_close_unit_check('qnto/nto_output_proj_coeff',trim(sr_out_filename),ierr)
    end if



    if(pub_on_root) then
       write(sr2_out_filename,'(a)') trim(pub_rootname)//'_QNTO_'//&
            trim(adjustl(name_num_core))//'a'//trim(adjustl(name_n))//'t_s2r_'//trim(ref_dir)//'.dat'
       sr2_out_filename = adjustl(sr2_out_filename)
       sr2_out_unit = utils_unit()
       open(unit=sr2_out_unit,form="formatted",status="old",position="append", &
            file=trim(sr2_out_filename),action="write",iostat=ierr)

       ! Write header
       write(sr2_out_unit,'(a,//a,/a,/a,i4,a,//a,/a,/a/)') &
            'Sc2Rc_H2H_E2E:', & !Projection of cNTOn-H(E)_S to cNTOn-H(E)_R:'
            &'   "<cNTO_R|cNTO_S> / (<cNTO_R|cNTO_R>)^0.5 * <cNTO_S|cNTO_S>^0.5)"', &
            '   using <NGWF_R|NGWF_R>, <NGWF_S|NGWF_S>, and <NGWF_R|NGWF_S> on S-geom', &
            '   (qnto_num_core_atoms = ', pub_qnto_num_core_atoms, ')', &
            'Map format:', &
            '            cNTO_S', &
            '   cNTO_R| H2H; E2E'

       write(sr2_out_unit,'(13x)',advance='no')
       do state_idx_sys=1, num_states
          do nto_idx_sys=1, num_nto_trans
             write(sr2_out_unit,'(3x,i1,a,i3.3)',advance='no') &
                  nto_idx_sys,'HE_',state_idx_sys
          end do
          write(sr2_out_unit,'(a)',advance='no') '|'
       end do
       write(sr2_out_unit,'(3x,a)') 'Matchings'

       ! Output transition coeff of NTOs_sys
       write(sr2_out_unit,'(12x,a)',advance='no') '|'
       do state_idx_sys=1, num_states
          do nto_idx_sys=1, num_nto_trans
             n = (state_idx_sys-1)*num_nto_trans + nto_idx_sys
             write(sr2_out_unit,'(5x,i3.2,2x)',advance='no') &
                  nint(trans_coeff_sys(n)*100.0_DP)
          end do
          write(sr2_out_unit,'(a)',advance='no') '|'
       end do
       write(sr2_out_unit,*)
    end if !pub_on_root


    ! Plot the map
    do irow=1, reftgt_size
       if(pub_on_root) then
          write(sr2_out_unit,'(1x,i1,a,i3.3,1x,i3.2,a)',advance='no') &
               irow-(irow-1)/num_nto_trans*num_nto_trans,'HE_', &
               (irow-1)/num_nto_trans+1,&
               nint(trans_coeff_ref(irow)*100.0_DP),'|'
          num_match = 0
          I_match(:) = 0
          nto_idx_ref = MOD(irow+num_nto_trans-1, num_nto_trans) + 1
          state_idx_ref = (irow+num_nto_trans-1) / num_nto_trans
       end if !pub_on_root

       ! 1st icol run
       do icol=1, systgt_size
          if(pub_on_root) then
             if(abs(h_sc2rc(irow,icol))>=sqrt(1.0_DP/2.0_DP) .AND. &
                abs(e_sc2rc(irow,icol))>=sqrt(1.0_DP/2.0_DP) ) then
                write(sr2_out_unit,'(a1)',advance='no') '*'
                num_match = num_match + 1
                I_match(num_match) = icol
             else if(abs(h_sc2rc(irow,icol))>=sqrt(1.0_DP/3.0_DP) .AND. &
                     abs(e_sc2rc(irow,icol))>=sqrt(1.0_DP/3.0_DP) ) then
                write(sr2_out_unit,'(a1)',advance='no') '#'
                num_match = num_match + 1
                I_match(num_match) = icol
             else
                write(sr2_out_unit,'(a1)',advance='no') ' '
             end if

             ! Finding the best match between sys and ref
             nto_idx_sys = MOD(icol+num_nto_trans-1, num_nto_trans) + 1
             state_idx_sys = (icol+num_nto_trans-1) / num_nto_trans
             match_weight(state_idx_sys) = match_weight(state_idx_sys) + &
               trans_coeff_sys(icol)*trans_coeff_ref(irow)*h_sc2rc(irow,icol)*e_sc2rc(irow,icol)* &
               nor_h_s2sc(icol)*nor_h_r2rc(irow)*nor_e_s2sc(icol)*nor_e_r2rc(irow)
             if (pub_debug_on_root) then
                if(state_idx_ref==state_idx_sys) then
                     write(stdout,*) 'state_idx_sys:', state_idx_sys, trans_coeff_sys(icol),&
                     trans_coeff_ref(irow), h_sc2rc(irow,icol), e_sc2rc(irow,icol), nor_h_s2sc(icol), nor_h_r2rc(irow),&
                     nor_e_s2sc(icol), nor_e_r2rc(irow), match_weight(state_idx_sys)
                end if
             end if
             if( nto_idx_sys == num_nto_trans .AND. &
                 nto_idx_ref == num_nto_trans ) then
                if( hash_flag_ref(state_idx_ref) /= '!' &
                  .AND. (label_flag_ref(state_idx_ref) .eqv. .TRUE.) &
                  .AND. (label_flag_sys(state_idx_sys) .eqv. .TRUE.) &
                  .AND. abs(match_weight(state_idx_sys)) > &
                        match_weight_sys(state_idx_sys) ) then
                   label_sys( state_idx_sys ) = label_ref( state_idx_ref )
                   match_energy_diff( state_idx_sys ) = &
                                    abs( exc_energies_sys(state_idx_sys)-&
                                    exc_energies_ref(state_idx_ref) )/( n_occ*0.2_DP )
                   match_weight_sys(state_idx_sys) = abs(match_weight(state_idx_sys))
                end if
                match_weight(state_idx_sys) = 0.0_DP
             end if

             if(h_sc2rc(irow,icol)>-0.005_DP .AND. h_sc2rc(irow,icol)<0.0_DP) then
                write(sr2_out_unit,'(a)',advance='no') ' -00'
             else
                write(sr2_out_unit,'(i4.2)',advance='no') nint(h_sc2rc(irow,icol)*100.0_DP)
             end if
             write(sr2_out_unit,'(a1)',advance='no') ';'

             if(e_sc2rc(irow,icol)>-0.005_DP .AND. e_sc2rc(irow,icol)<0.0_DP) then
                write(sr2_out_unit,'(a)',advance='no') ' -00'
             else
                write(sr2_out_unit,'(i4.2)',advance='no') nint(e_sc2rc(irow,icol)*100.0_DP)
             end if

             if(MOD(icol,num_nto_trans)==0) &
             write(sr2_out_unit,'(a)',advance='no') '|'
          end if !pub_on_root
       end do
       ! 1st icol run finished

       if(pub_on_root) then
          write(sr2_out_unit,'(2x,i3,1x,a,i1,a,i3.3,a)',advance='no') &
               num_match, 'Rc', nto_idx_ref, '_', state_idx_ref, ' - Sc'
          do i=1, num_match
             if(I_match(i) > 0) then
                nto_idx_sys = MOD(I_match(i)+num_nto_trans-1, &
                              num_nto_trans) + 1
                state_idx_sys = (I_match(i)+num_nto_trans-1) / &
                            num_nto_trans
                write(sr2_out_unit,'(a,i1,a,i3.3)',advance='no') ' *', &
                     nto_idx_sys, '_', state_idx_sys
             else
                nto_idx_sys = MOD(I_match(i)+systgt_size+num_nto_trans-1, &
                              num_nto_trans) + 1
                state_idx_sys = (I_match(i)+systgt_size+num_nto_trans-1) / &
                                num_nto_trans
                write(sr2_out_unit,'(a,i1,a,i3.3)',advance='no') ' #', nto_idx_sys, &
                     '_', state_idx_sys
             end if
          end do
          write(sr2_out_unit,*)
       end if !pub_on_root
    end do !irow
    if(pub_on_root) then
       close(unit=sr2_out_unit,iostat=ierr)
       call utils_close_unit_check('qnto/nto_output_proj_coeff',trim(sr2_out_filename),ierr)
    end if


    deallocate(I_match,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','I_match', ierr)
    deallocate(exc_energies_sys,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','exc_energies_sys', ierr)
    deallocate(exc_energies_ref,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','exc_energies_ref', ierr)
    deallocate(match_weight,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','match_weight', ierr)
    deallocate(nor_h_rc,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','nor_h_rc', ierr)
    deallocate(nor_e_rc,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','nor_e_rc', ierr)
    deallocate(h_sc2rc,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','h_sc2rc', ierr)
    deallocate(e_sc2rc,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','e_sc2rc', ierr)
    deallocate(nor_h_s2sc,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','nor_h_s2sc', ierr)
    deallocate(nor_e_s2sc,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','nor_e_s2sc', ierr)
    deallocate(nor_h_r2rc,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','nor_h_r2rc', ierr)
    deallocate(nor_e_r2rc,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','nor_e_r2rc', ierr)
    deallocate(nor_h_sc,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','nor_h_sc', ierr)
    deallocate(nor_e_sc,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','nor_e_sc', ierr)
    deallocate(trans_coeff_sys,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','trans_coeff_sys', ierr)
    deallocate(trans_coeff_ref,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','trans_coeff_ref', ierr)

    !------------------------------------------------------------------------------------!
    ! Output label_sys(:)                                                                !
    !------------------------------------------------------------------------------------!
    if(pub_on_root) then
       write(label_io_filename,'(a)') trim(pub_rootname)//'_state_labellings.dat.tmp'
       label_io_filename = adjustl(label_io_filename)
       label_io_unit = utils_unit()
       open(unit=label_io_unit,form="formatted", file=trim(label_io_filename),&
            action="write", position="append", iostat=ierr)
       call utils_open_unit_check('qnto/nto_output_proj_coeff',&
            trim(label_io_filename),ierr)

       if(num_states > pub_qnto_num_ref_states) then
          n = num_states
       else
          n = pub_qnto_num_ref_states
       end if
       write(label_io_unit,'(a,2x,a)',advance='no') '%', '|'
       do state_idx=1, n
          write(label_io_unit,'(i5.5,a)',advance='no') state_idx, '|'
       end do
       write(label_io_unit,*)
       write(label_io_unit,'(a,2x,a)',advance='no') '%', '|'
       if(pub_on_root)then!debug_on_root) then
          write(stdout,*) 'match_weight_sys(1):', match_weight_sys(1)
          write(stdout,*) 'label_ref(1):', label_ref(1)
          write(stdout,*) 'label_sys(1):', label_sys(1)
          write(stdout,*) 'global_match_weight_sys(1):', global_match_weight_sys(1)
          write(stdout,*) 'global_label_ref(1):', global_label_ref(1)
          write(stdout,*) 'global_label_sys(1):', global_label_sys(1)
          write(stdout,*)
       end if

       do state_idx_ref=1, pub_qnto_num_ref_states
          if( hash_flag_ref(state_idx_ref) == '1' .OR. &
              hash_flag_ref(state_idx_ref) == ' ' ) then
             if(label_flag_ref(state_idx_ref) .eqv. .FALSE.) then
                write(label_io_unit,'(a)',advance='no') '--   |'
             else
                write(label_io_unit,'(a,3x,a)',advance='no') label_ref(state_idx_ref), '|'
             end if
          else
             write(label_io_unit,'(a,a,2x,a)',advance='no') label_ref(state_idx_ref), &
                  hash_flag_ref(state_idx_ref), '|'
          end if
       end do
       write(label_io_unit,'(a,1x,a,a,a,a,a,a)') 'Ref:', trim(ref_dir), &
            '[', trim(adjustl(name_num_core)), 'at;', trim(adjustl(name_n)), 'tr]'

       write(label_io_unit,'(a)',advance='no') '-> |'
       do state_idx_sys=1, num_states
          if( label_flag_sys(state_idx_sys) .eqv. .FALSE. ) then
             write(label_io_unit,'(a)',advance='no') '--   |'
             global_label_sys(state_idx_sys) = '--'
          else if( match_weight_sys(state_idx_sys) > sqrt(1.0_DP/2.0_DP) ) then
             write(label_io_unit,'(a,i3.2,a)',advance='no') label_sys(state_idx_sys), &
                  nint(match_weight_sys(state_idx_sys)*100.0_DP), '|'
          else if( match_weight_sys(state_idx_sys) > sqrt(1.0_DP/3.0_DP) ) then
             write(label_io_unit,'(a,a,i2.2,a)',advance='no') label_sys(state_idx_sys), &
                  '#', nint(match_weight_sys(state_idx_sys)*100.0_DP), '|'
          else
             write(label_io_unit,'(a,a,i2.2,a)',advance='no') label_sys(state_idx_sys), &
                  '!', nint(match_weight_sys(state_idx_sys)*100.0_DP), '|'
          end if

          ! Set up global_label_sys(:)
          if(state_idx_sys <= pub_qnto_num_ref_states) then
             if( match_weight_sys(state_idx_sys) > sqrt(1.0_DP/3.0_DP) ) then
                if( label_sys(state_idx_sys) /= label_ref(state_idx_sys) ) then
                   ! Global sys and ref having the same label is given the lowest priority
                   if( global_label_sys(state_idx_sys) == global_label_ref(state_idx_sys) ) then
                      global_match_weight_sys(state_idx_sys) = match_weight_sys(state_idx_sys)
                      global_label_sys(state_idx_sys) = label_sys(state_idx_sys)
                      global_label_ref(state_idx_sys) = label_ref(state_idx_sys)
                   else !global_label_sys(state_idx_sys) /= global_label_ref(state_idx_sys)
                     if(global_match_weight_sys(state_idx_sys) < &
                                match_weight_sys(state_idx_sys)) then
                      global_match_weight_sys(state_idx_sys) = match_weight_sys(state_idx_sys)
                      global_label_sys(state_idx_sys) = label_sys(state_idx_sys)
                      global_label_ref(state_idx_sys) = label_ref(state_idx_sys)
                     end if
                   end if
                else !label_sys(state_idx_sys) == label_ref(state_idx_sys)
                   if( global_label_sys(state_idx_sys) == global_label_ref(state_idx_sys) ) then
                      if(global_match_weight_sys(state_idx_sys) < &
                                match_weight_sys(state_idx_sys)) then
                      global_match_weight_sys(state_idx_sys) = match_weight_sys(state_idx_sys)
                      global_label_sys(state_idx_sys) = label_sys(state_idx_sys)
                      global_label_ref(state_idx_sys) = label_ref(state_idx_sys)
                      end if
                   end if
                end if
             end if !match_weight_sys(state_idx_sys)>sqrt(1.0_DP/3.0_DP)
          else !state_idx_sys > pub_qnto_num_ref_states
             if( match_weight_sys(state_idx_sys) > sqrt(1.0_DP/3.0_DP) ) then
                if(global_match_weight_sys(state_idx_sys) < match_weight_sys(state_idx_sys)) then
                   global_match_weight_sys(state_idx_sys) = match_weight_sys(state_idx_sys)
                   global_label_sys(state_idx_sys) = label_sys(state_idx_sys)
                end if
             end if
          end if
       end do
       if(pub_on_root)then!debug_on_root) then
          write(stdout,*) 'new global_match_weight_sys(1):', global_match_weight_sys(1)
          write(stdout,*) 'new global_label_ref(1):', global_label_ref(1)
          write(stdout,*) 'new global_label_sys(1):', global_label_sys(1)
       end if

       write(label_io_unit,*)
       write(label_io_unit,'(a)',advance='no') '%'
       close(unit=label_io_unit,iostat=ierr)
       call utils_close_unit_check('qnto/nto_output_proj_coeff',&
            trim(label_io_filename),ierr)
    end if !pub_on_root

    deallocate(label_flag_sys,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','label_flag_sys', ierr)
    deallocate(label_flag_ref,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','label_flag_ref', ierr)
    deallocate(label_sys,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','label_sys', ierr)
    deallocate(label_ref,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','label_ref', ierr)
    deallocate(hash_flag_sys,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','hash_flag_sys', ierr)
    deallocate(hash_flag_ref,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','hash_flag_ref', ierr)
    deallocate(match_weight_sys,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','match_weight_sys', ierr)
    deallocate(match_energy_diff,stat=ierr)
    call utils_dealloc_check('qnto/nto_output_proj_coeff','match_energy_diff', ierr)

  end subroutine nto_output_proj_coeff
end module qnto
