! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-   !

!==================================================================!
!               Module: sparse_initialise_mod.f90                  !
!                                                                  !
! This module was created by Nicholas Hine in April 2011 based on  !
! routines from sparse_mod.F90, created by Peter Haynes and later  !
! extended by Nicholas Hine from 2007.                             !
! Further subroutines have been written by Nicholas Hine,          !
! Jacek Dziedzic, Quintin Hill, David O'Regan, Laura Ratcliff      !
! and Tim Zuehlsdorff between April 2011 and April 2014.           !
! Modified for embedding structures by Robert Charlton, 2017.      !
!                                                                  !
!==================================================================!
module sparse_initialise

  use constants, only: DP, LONG, stdout
  use ion, only: ELEMENT

  implicit none

  public :: sparse_init_rep

  ! rc2013: wrapper for subsystem info
  type SFC_INFO
     type(ELEMENT), allocatable   ::   elems_sfc(:)
  end type SFC_INFO

  private

contains

  !===========================================================================!
  ! This subroutine performs all the initialisation routines involving the    !
  ! sparse matrices for a given set of NGWFs (eg valence, conduction, joint). !
  !---------------------------------------------------------------------------!
  ! Arguments:                                                                !
  !   elements (input) : A list of atoms in the cell                          !
  !   suffix (input)   : The suffix for this NGWF representation.             !
  !---------------------------------------------------------------------------!
  ! SPAM3 version written by Nicholas Hine, May 2009.                         !
  !---------------------------------------------------------------------------!
  ! SPAM2 version originally written by Peter Haynes, March 2004.             !
  ! Revised for distributed data, June 2006.                                  !
  ! Indexing for dense matrices added by Nicholas Hine, Dec 2007              !
  ! Modified Hamiltonian sparsity for HF exchange, Quintin Hill Feb 2008      !
  ! Modified to include conduction matrix types by Laura Ratcliff, Oct 2010   !
  ! Massive re-write by Nicholas Hine, April 2011.                            !
  ! Modified to handle cross-overlap matrices by Robert Charlton, Sept 2018.  !                                                  !
  !===========================================================================!

  subroutine sparse_init_rep(mdl, suffix)

    use comms, only: pub_on_root, pub_my_proc_id, &
         pub_total_num_procs, comms_barrier
    use constants, only: VERBOSE
    use ion, only : ELEMENT
    use model_type, only: MODEL
    use parallel_strategy, only: PARAL_INFO, OVERLAP_LIST, &
         parallel_strategy_list_cross_overlaps, parallel_strategy_exit_overlaps
    use rundat, only: pub_use_swx, pub_kernel_cutoff, &
         pub_any_nl_proj, pub_paw, pub_aug, pub_hubbard, pub_output_detail, &
         cond_kernel_cutoff, pub_hfx_nlpp_for_exchange, &
         pub_eels_calculate, pub_lr_tddft_calculate, pub_devel_code, &
         pub_lr_tddft_kernel_cutoff, pub_use_hfx, pub_lr_tddft_sparse_region, &
         pub_dos_smear, pub_pdos_max_l, pub_lr_phonons_calculate, &
         pub_lr_phonons_kernel_cutoff, pub_contracoham_radmult, &
         pub_H2denskern_sparsity, &
         pub_use_activehfx, pub_active_region, pub_use_activeswx
    use simulation_cell, only: CELL_INFO
    use sparse, only: BLKS_NGWF, BLKS_COND, BLKS_JOINT, BLKS_PROJ, &
         BLKS_HUB_PROJ, BLKS_SW, BLKS_CORE, BLKS_AUX, sparse_fill_fac_denom, &
         sparse_num_elems_on_atom, sparse_count_ss, sparse_index_ss, &
         sparse_count_union, sparse_index_union, sparse_num_element, &
         sparse_show_segment_filling, pub_sparse_allow_new_matrices, &
         BLKS_PDOS, BLKS_JOINTPH
    use sparse_embed, only: SPAM3_EMBED,sparse_embed_create,sparse_embed_destroy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_devel_code, utils_banner !, utils_unit
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Argument
    type(MODEL), intent(in), target :: mdl
    character(*), intent(in)        :: suffix

    ! Local variables
    integer :: ierr                    ! Error flag
    integer :: iat                     ! Atom loop counter
    integer :: proc                    ! Proc loop counter
    integer :: iatproc                 ! Atom on proc loop counter
    integer :: iblk                    ! Block counter
    integer :: nzb,my_nze,my_nzb       ! Various sizes of current matrix
    integer(kind=LONG) :: nze          ! Number of elems of current matrix
    integer(kind=LONG) :: nze_q,nze_r  ! Number of elems of q and r matrices
    integer(kind=LONG) :: nze_b,nze_c  ! Number of elems of b and c matrices
    integer :: bslib(mdl%nsub,mdl%nsub) ! Library entry for basic NGWF overlap
    integer :: slib(mdl%nsub,mdl%nsub) ! Library entry for overlap
    integer :: occhlib(mdl%nsub,mdl%nsub) ! Library entry for basic NGWF overlap (for Contra-Covariant Ham) <-ja531
    integer :: scchlib(mdl%nsub,mdl%nsub) ! Library entry for overlap (for Contra-Covariant Ham)
    integer, dimension(mdl%nsub,mdl%nsub) :: klib       ! Library entry for K
    integer, dimension(mdl%nsub,mdl%nsub) :: qrcchlib
    integer, dimension(mdl%nsub,mdl%nsub) :: qrlib,vwlib ! Library entries for H intermediates
    integer, dimension(mdl%nsub,mdl%nsub) :: qcrlib,qrclib ! Library entries for Cond-NGWF crossterms
    integer, dimension(mdl%nsub,mdl%nsub) :: qprlib,prqlib           ! Library entries for Proj-SW intermediates
    integer, dimension(mdl%nsub,mdl%nsub) :: swngwflib,ngwfswlib     ! Library entries for NGWF-SW intermediates
    integer :: hlib(mdl%nsub,mdl%nsub)   ! Library entry for H
    integer :: cchlib(mdl%nsub,mdl%nsub) ! Library entry for Contra-Covariant H
    integer, dimension(mdl%nsub,mdl%nsub) :: kslib,kskslib           ! Library entries for kernel products
    integer :: lr_tddft_proj_kernel_lib(mdl%nsub,mdl%nsub)! Library entry for projected TDDFT kernel
    integer :: K1SvKv_lib(mdl%nsub,mdl%nsub)  ! Library entry for partially projected kernel
    integer, dimension(mdl%nsub,mdl%nsub) :: ksksklib
    integer, dimension(mdl%nsub,mdl%nsub) :: tlib, ulib
    integer, dimension(mdl%nsub,mdl%nsub) :: btlib, bulib
    integer, dimension(mdl%nsub,mdl%nsub) :: cond_proj_ham_lib
    integer, dimension(mdl%nsub,mdl%nsub) :: proj_embed_ham_lib
    integer, dimension(mdl%nsub,mdl%nsub) :: cond_grad_lib1, cond_grad_lib2
    integer :: embed_grad_lib1, embed_grad_lib2
    integer :: blks
    character(len=10) :: pscode
    real(kind=DP) :: fill_fac          ! Denominator for filling fractions
    integer,allocatable :: seg_nzb(:)  ! Number of nonzero blocks on each seg
    integer,allocatable :: seg_nze(:)  ! Number of nonzero elements on each seg
    type(ELEMENT), allocatable :: elems_sfc1(:), elems_sfc2(:)
    type(OVERLAP_LIST) :: overlaps(mdl%nsub,mdl%nsub)
    logical :: any_proj
    logical :: show_sparsity
    type(SFC_INFO), allocatable :: system_info(:)
    type(PARAL_INFO), pointer :: row_par, col_par
    integer :: first_atom_on_proc, last_atom_on_proc
    integer :: isub, jsub
    character(len=1) :: isub_str, jsub_str
    character(len=2) :: reg_struc, reg_struc_t

    ! Start Timer (and initially synchronise all procs)
    call comms_barrier
    call timer_clock('sparse_init',1)
#ifdef ITC_TRACE
    call VTBEGIN(vt_sparse_init,vt_err)
#endif

    pub_sparse_allow_new_matrices = .true.

    any_proj = (pub_any_nl_proj.or.pub_paw)

    show_sparsity = (pub_output_detail>=VERBOSE)

    ! Allocate workspace for this routine
    allocate(system_info(mdl%nsub),stat=ierr)
    call utils_alloc_check('sparse_init','system_info',ierr)
    do isub=1,mdl%nsub
       ! rc2013: now allocate space for elements in each subsystem
       allocate(system_info(isub)%elems_sfc(mdl%regions(isub)%par%nat),stat=ierr)
       call utils_alloc_check('sparse_init','system_info%elems_sfc',ierr)
    end do
    allocate(seg_nzb(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('sparse_init','seg_nzb',ierr)
    allocate(seg_nze(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('sparse_init','seg_nze',ierr)

    ! Transfer atoms to space-filling curve order
    do isub=1,mdl%nsub
       ! rc2013: reset iblk for each subsystem
       iblk = 1
       do proc=0,pub_total_num_procs-1
          do iatproc=1,mdl%regions(isub)%par%num_atoms_on_proc(proc)
             iat = mdl%regions(isub)%par%atoms_on_proc(iatproc,proc)
             system_info(isub)%elems_sfc(iblk) = mdl%regions(isub)%elements(iat)
             iblk = iblk + 1
          end do
       end do
    end do

    ! Sparse matrix information header
    if (pub_on_root.and.show_sparsity) then
       if (suffix=='c') then
          write(stdout,'(/a)') utils_banner('=', &
               'Conduction Sparse matrix information')
       else if (suffix=='j') then
          write(stdout,'(/a)') utils_banner('=', &
               'Joint Valence + Conduction Sparse matrix information')
       else if (suffix=='z') then
          write(stdout,'(/a)') utils_banner('=', &
               'JointPH Sparse matrix information')
       else if (suffix=='a') then
          write(stdout,'(/a)') utils_banner('=', &
               'Auxiliary NGWF Sparse matrix information')
       else
          write(stdout,'(/a)') utils_banner('=', 'Sparse matrix information')
       end if
    end if

    !---------------------------------------------------------------------------
    ! NGWF overlaps and kernel
    !---------------------------------------------------------------------------

    if (suffix=='') then
       blks = BLKS_NGWF
       pscode = 'Region'
    else if (suffix=='c') then
       blks = BLKS_COND
       pscode = 'A'
    else if (suffix=='j') then
       blks = BLKS_JOINT
       pscode = 'A'
    else if (suffix=='z') then
       blks = BLKS_JOINTPH
       pscode = 'Region'
    else if (suffix=='a') then
       blks = BLKS_AUX
       pscode = 'Region'
    else if(suffix=='p') then
       blks = BLKS_PDOS
       pscode = 'Region'
    end if

    call internal_create_representation

    if(.not.pub_H2denskern_sparsity) then
       ! Create other matrices: ks,ksk,ksks
       if (suffix/='a') call internal_create_ksksk
    end if

    !---------------------------------------------------------------------------
    ! Hamiltonian matrices
    !---------------------------------------------------------------------------

    call internal_create_hamiltonian



    if(pub_H2denskern_sparsity) then
       ! Create other matrices: ks,ksk,ksks
       if (suffix/='a') call internal_create_ksksk
    end if



    call internal_create_ham_strucs
    if (any_proj) call internal_create_nl_strucs

    ! Create Hubbard structures if required
    if (pub_hubbard.and.(suffix=='')) then
       call internal_create_hubbard_strucs
    end if

    if ((suffix=='c').or.(suffix=='j').or.(suffix=='a')) then
       call internal_create_cond_structures
       ! tjz07: ADD option for TDDFT calculations
       if (pub_lr_tddft_calculate) then
          call internal_create_tddft_structures
       endif
    end if


    ! gcc32: LR_PHONONS
    if (suffix=='z') then
       call internal_create_resp_structures
    end if
    if (pub_lr_phonons_calculate) then
       call internal_create_lr_phonons_structures
    endif

    !---------------------------------------------------------------------------
    ! Block diagonal matrices
    !---------------------------------------------------------------------------

    ! rc2013: initialise projector matrix etc.
    call internal_create_block_diag

    !---------------------------------------------------------------------------
    ! Display filling factors of other structures
    !---------------------------------------------------------------------------

    ! rc2013: EMBED_FIX
    ! rc2013: should rewrite this for full matrix rather than sub-blocks
    ! PROBLEM: we use the library entries for ks and ksks, which is a mess...
    do isub=1,mdl%nsub
       do jsub=1,mdl%nsub
          col_par => mdl%regions(isub)%par
          row_par => mdl%regions(jsub)%par

          ! Calculate filling factor denominator for fraction
          fill_fac = sparse_fill_fac_denom(blks,blks,row_par,col_par)

          ! Find the number of nonzero elements in the patterns
          ! of commonly-used matrices and print results as part of filling
          ! factor summary
          if (pub_on_root.and.show_sparsity.and.(suffix/='a')) then
             if (mdl%nsub.eq.1) then
                write(stdout,'(29x,a,f7.2,2a)') 'KS matrix filling:', &
                     fill_fac * sparse_num_element(kslib(isub,jsub)),'% '
                write(stdout,'(27x,a,f7.2,2a)') 'KSKS matrix filling:', &
                     fill_fac * sparse_num_element(kskslib(isub,jsub)),'% '
             else
                write(stdout,'(23x,a,i0,a,i0,a,f7.2,2a)') &
                     'KS matrix filling (',isub,',',jsub,'):', &
                     fill_fac * sparse_num_element(kslib(isub,jsub)),'% '
                write(stdout,'(21x,a,i0,a,i0,a,f7.2,2a)') &
                     'KSKS matrix filling (',isub,',',jsub,'):', &
                     fill_fac * sparse_num_element(kskslib(isub,jsub)),'% '
             end if
          end if
       enddo
    enddo

    ! jd: Create V and/or O matrix for SW expansion
    if (pub_use_swx.or.pub_use_activeswx) call internal_create_metric_matrix
    ! jd: Create X (exchange) matrix for HFx
    if (pub_use_hfx.or.pub_use_activehfx) call internal_create_exchange_matrix

    ! Sparse matrix information footer
    if (pub_on_root.and.show_sparsity) &
         write(stdout,'(a/)') repeat('=',80)

    if (utils_devel_code(.false.,'SPARSE','SHOW_SEGS',pub_devel_code)) then
       call sparse_show_segment_filling
    end if

    pub_sparse_allow_new_matrices = .false.

    call comms_barrier

    ! Deallocate workspace for this routine
    do isub=1,mdl%nsub
       do jsub=1,mdl%nsub
          call parallel_strategy_exit_overlaps(overlaps(jsub,isub))
       end do
       deallocate(system_info(isub)%elems_sfc,stat=ierr)
       call utils_dealloc_check('sparse_init','system_info%elems_sfc',ierr)
    end do
    deallocate(seg_nze,stat=ierr)
    call utils_dealloc_check('sparse_init','seg_nze',ierr)
    deallocate(seg_nzb,stat=ierr)
    call utils_dealloc_check('sparse_init','seg_nzb',ierr)
    deallocate(system_info,stat=ierr)
    call utils_dealloc_check('sparse_init','system_info',ierr)

    ! Stop Timer
    call timer_clock('sparse_init',2)
#ifdef ITC_TRACE
    call VTEND(vt_sparse_init,vt_err)
#endif

  contains


    !==========================================================================!
    ! This subroutine creates the structures of an NGWF representation: the    !
    ! overlap, density kernel, and the sp- and ps- overlaps with the nonlocal  !
    ! and Hubbard projectors.                                                  !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in April 2011 based on existing code.           !
    !==========================================================================!

    subroutine internal_create_representation

      implicit none

      ! Arguments

      ! Local variables
      real(kind=DP) :: range
      character(len=30) :: sstruc

      do isub=1,mdl%nsub
         call internal_allocate_sfc(elems_sfc1, col_par, isub)
         do jsub=1,mdl%nsub
            call internal_allocate_sfc(elems_sfc2, row_par, jsub)
            call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

            ! Calculate filling factor denominator for fraction
            fill_fac = sparse_fill_fac_denom(blks,blks,row_par,col_par)

            !---------------------------------------------------------------------------
            ! Overlap matrix
            !---------------------------------------------------------------------------

            ! Obtain a list of overlaps between spherical regions
            call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub), &
                 elems_sfc1, elems_sfc2, mdl%cell,pscode,pscode)

            ! Count number of corresponding nonzero NGWF matrix elements
            call sparse_count_ss(overlaps(jsub,isub),blks,blks,nze,nzb, &
                 my_nze,my_nzb,seg_nze,seg_nzb,row_par,col_par)

            ! Create sparse overlap matrix
            if (.not.pub_aug) then
               sstruc = 'S'//suffix
            else
               sstruc = 'bS'//suffix ! records library index of "bare" NGWF overlap
            end if
            call sparse_index_ss(overlaps(jsub,isub),blks,blks, &
                 trim(sstruc)//reg_struc, &
                 nze,nzb,my_nze, my_nzb,seg_nze,seg_nzb,slib(jsub,isub), &
                 trim(sstruc)//reg_struc_t, row_par, col_par)
            if (pub_aug) bslib = slib

            if (pub_on_root.and.show_sparsity) then
               if (mdl%nsub.eq.1) then
                  write(stdout,'(24x,a,f7.2,a)') 'Overlap matrix filling:', &
                       fill_fac * nze,'%'
               else
                  write(stdout,'(18x,a,i0,a,i0,a,f7.2,a)') &
                       'Overlap matrix filling (',jsub,',',isub,'):', &
                       fill_fac * nze,'%'
               end if
            end if

            ! Synchronise
            call comms_barrier


            if(abs(pub_contracoham_radmult-1.0_dp)>epsilon(1.0_dp)) then
               !---------------------------------------------------------------
               ! Overlap matrix for Contra-Covariant Ham
               !---------------------------------------------------------------

               ! Obtain a list of overlaps between spherical regions
               call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub), &
                    elems_sfc1, elems_sfc2, mdl%cell, pscode, pscode, &
                    multiplier=pub_contracoham_radmult)

               ! Count number of corresponding nonzero NGWF matrix elements
               call sparse_count_ss(overlaps(jsub,isub),blks,blks,nze,nzb, &
                    my_nze,my_nzb,seg_nze,seg_nzb,row_par,col_par)

               ! Create sparse overlap matrix
               if (.not.pub_aug) then
                  call sparse_index_ss(overlaps(jsub,isub),blks,blks, &
                       'SccH'//suffix//reg_struc, &
                       nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb,rlib=scchlib(jsub,isub), &
                       transpose_name='SccH'//suffix//reg_struc, &
                       row_par=row_par, col_par=col_par)
               else
                  call sparse_index_ss(overlaps(jsub,isub),blks,blks, &
                       'OccH'//suffix//reg_struc, &
                       nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb,rlib=occhlib(jsub,isub), &
                       transpose_name='OccH'//suffix//reg_struc, &
                       row_par=row_par, col_par=col_par)
               end if

            end if

            ! Synchronise
            call comms_barrier



            !------------------------------------------------------------------
            ! Density kernel
            !------------------------------------------------------------------

            !ja531-> Non-H^2 Sparsity
            if(.not.pub_H2denskern_sparsity) then

               ! Obtain a list of overlaps of fixed cutoff
               if ((suffix=='c').or.(suffix=='j')) then
                  range = cond_kernel_cutoff
               else
                  range = pub_kernel_cutoff
               end if
               call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub), &
                    elems_sfc1, elems_sfc2, mdl%cell, 'Fixed','Fixed',0.5_DP * range)

               ! Count number of corresponding nonzero NGWF matrix elements
               call sparse_count_ss(overlaps(jsub,isub),blks,blks,nze,nzb,my_nze,my_nzb, &
                    seg_nze,seg_nzb,row_par,col_par)

               ! Create density kernel matrix
               call sparse_index_ss(overlaps(jsub,isub),blks,blks, &
                    'K'//suffix//reg_struc, &
                    nze,nzb,my_nze, my_nzb, seg_nze,seg_nzb,klib(jsub,isub), &
                    'K'//suffix//reg_struc_t, row_par, col_par)

               if (pub_on_root.and.show_sparsity) then
                  if (mdl%nsub.eq.1) then
                     write(stdout,'(24x,a,f7.2,a)') 'Density kernel filling:', &
                          fill_fac * nze,'%'
                  else
                     write(stdout,'(18x,a,i0,a,i0,a,f7.2,a)') &
                          'Density kernel filling (',jsub,',',isub,'):', &
                          fill_fac * nze,'%'
                  end if
               end if


               ! Create W and V matrices if required
               if ((pub_hubbard).and.(suffix=='')) then

                  ! Count number of corresponding nonzero HUB_PROJ  matrix elements
                  call sparse_count_ss(overlaps(jsub,isub),BLKS_HUB_PROJ,BLKS_HUB_PROJ, &
                       nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb, &
                       row_par=row_par,col_par=col_par)

                  ! Create density kernel matrix pattern in Hubbard Projector block scheme
                  call sparse_index_ss(overlaps(jsub,isub),BLKS_HUB_PROJ,BLKS_HUB_PROJ, &
                       'I'//suffix//reg_struc_t,nze,nzb,my_nze,my_nzb, &
                       seg_nze,seg_nzb, row_par=row_par, col_par=col_par)

               end if

            end if

            ! Synchronise
            call comms_barrier
            ! rc2013: deallocate temporary arrays
            deallocate(elems_sfc2,stat=ierr)
            call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
         end do
         deallocate(elems_sfc1,stat=ierr)
         call utils_dealloc_check('sparse_init_rep', 'elems_sfc', ierr)
      end do

      !---------------------------------------------------------------------------
      ! Projector-NGWF matrices
      !---------------------------------------------------------------------------

      ! Set up these matrices only if there are any projectors
      if (any_proj) then
         do isub=1,mdl%nsub
            call internal_allocate_sfc(elems_sfc1, col_par, isub)
            do jsub=1,mdl%nsub
               call internal_allocate_sfc(elems_sfc2, row_par, jsub)
               call internal_write_embed_struc(reg_struc,reg_struc_t,jsub,isub)

               ! Create Projector-NGWF overlap matrix
               call internal_create_structure(overlaps(jsub,isub),elems_sfc1, &
                    elems_sfc2,pscode,'Core',blks,BLKS_PROJ,nze_r,nzb, &
                    my_nze,my_nzb,seg_nze,seg_nzb,col_par,row_par, &
                    'R'//suffix//reg_struc,'Q'//suffix//reg_struc_t)

               ! Synchronise
               call comms_barrier

               ! Create NGWF-Projector overlap matrix
               call internal_create_structure(overlaps(jsub,isub),elems_sfc1, &
                    elems_sfc2,'Core',pscode,BLKS_PROJ,blks,nze_q,nzb, &
                    my_nze,my_nzb,seg_nze,seg_nzb,col_par,row_par, &
                    'Q'//suffix//reg_struc,'R'//suffix//reg_struc_t)

               ! Calculate filling factor denominator for fraction
               fill_fac = sparse_fill_fac_denom(blks,BLKS_PROJ,row_par,col_par)

               if (pub_on_root.and.show_sparsity) then
                  if (mdl%nsub.eq.1) then
                     write(stdout,'(14x,a,f7.2,a)') 'NGWF-Proj overlap matrix &
                          &filling:', fill_fac * nze_r,'%'
                  else
                     write(stdout,'(8x,a,i0,a,i0,a,f7.2,a)') &
                          'NGWF-Proj overlap matrix filling (', &
                          jsub,',',isub,'):', fill_fac * nze_r,'%'
                  end if
               end if

               ! rc2013: deallocate elements lists
               deallocate(elems_sfc2,stat=ierr)
               call utils_dealloc_check('sparse_init_rep', 'elems_sfc', ierr)
            end do ! jsub
            deallocate(elems_sfc1,stat=ierr)
            call utils_dealloc_check('sparse_init_rep', 'elems_sfc', ierr)
         end do ! isub
         ! Synchronise
         call comms_barrier

         ! Create matrix product of Q and R matrices as nonlocal psp matrix
         call internal_create_product(qrlib,'Q'//suffix,'R'//suffix)

         ! Augmented overlap 'S' is union of direct sphere-sphere S and QR
         if (pub_aug) then
            do isub=1,mdl%nsub
               col_par => mdl%regions(isub)%par
               do jsub=1,mdl%nsub
                  row_par => mdl%regions(jsub)%par
                  call internal_write_embed_struc(reg_struc,reg_struc_t,jsub,isub)

                  ! Count elements in union of bS and QR
                  call sparse_count_union(bslib(jsub,isub),qrlib(jsub,isub), &
                       nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)

                  ! S structure becomes this union
                  call sparse_index_union(bslib(jsub,isub),qrlib(jsub,isub), &
                       'S'//suffix//reg_struc,nze,nzb,my_nze, &
                       my_nzb,seg_nze,seg_nzb,rlib=slib(jsub,isub))

                  ! Calculate filling factor denominator for fraction
                  fill_fac = sparse_fill_fac_denom(blks,blks,row_par,col_par)

                  if (pub_on_root.and.show_sparsity) then
                     if (mdl%nsub.eq.1) then
                        write(stdout,'(20x,a,f7.2,a)') 'Aug Overlap matrix filling:', &
                             fill_fac * nze,'%'
                     else
                        write(stdout,'(14x,a,i0,a,i0,a,f7.2,a)') &
                             'Aug Overlap matrix filling (',jsub,',',isub,'):', &
                             fill_fac * nze,'%'
                     end if
                  end if

               end do
            end do
         end if

      else
         call utils_assert(.not. pub_aug, 'Error in internal_create_&
              &representation(): Cannot have no projectors in augmented &
              &formalisms.')
         qrlib = 0
      end if

      ! Set up a version of slib for contra-co variant Ham
      if(abs(pub_contracoham_radmult-1.0_dp)>epsilon(1.0_dp)) then

         ! Set up these matrices only if there are any projectors
         if (any_proj) then
            do isub=1,mdl%nsub
               call internal_allocate_sfc(elems_sfc1, col_par, isub)
               do jsub=1,mdl%nsub
                  call internal_allocate_sfc(elems_sfc2, row_par, jsub)
                  call internal_write_embed_struc(reg_struc, reg_struc_t, &
                       jsub, isub)

                  ! Create Projector-NGWF overlap matrix
                  call internal_create_structure(overlaps(jsub,isub), elems_sfc1, &
                       elems_sfc2, pscode, 'Core', blks, BLKS_PROJ, nze_r, nzb, &
                       my_nze, my_nzb, seg_nze, seg_nzb, col_par, row_par, &
                       'ccR'//suffix//reg_struc, &
                       'ccQ'//suffix//reg_struc_t,multiplier=pub_contracoham_radmult)

                  ! Synchronise
                  call comms_barrier

                  ! Create NGWF-Projector overlap matrix
                  call internal_create_structure(overlaps(jsub,isub), elems_sfc1, &
                       elems_sfc2, 'Core', pscode, BLKS_PROJ, blks, nze_q, nzb, &
                       my_nze, my_nzb, seg_nze, seg_nzb, col_par, row_par, &
                       'ccQ'//suffix//reg_struc, &
                       'ccR'//suffix//reg_struc_t,multiplier=pub_contracoham_radmult)

!                  ! Calculate filling factor denominator for fraction
!                  fill_fac = sparse_fill_fac_denom(blks,BLKS_PROJ)
!
!                  if (pub_on_root.and.show_sparsity) &
!                       write(stdout,'(14x,a,f7.2,a)') 'NGWF-Proj overlap matrix &
!                       &filling:', fill_fac * nze_r,'%'

                  ! Synchronise
                  call comms_barrier
                  ! rc2013: deallocate temporary arrays
                  deallocate(elems_sfc2,stat=ierr)
                  call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
               end do
               deallocate(elems_sfc1,stat=ierr)
               call utils_dealloc_check('sparse_init_rep', 'elems_sfc', ierr)
            end do

            ! rc2013: don't see why this'd be the case with cross overlaps?
            call utils_assert(nze_q == nze_r, 'Error in internal_create_&
                 &representation(): NGWF-Projector and Projector-NGWF overlap &
                 &matrices must have same number of nonzero elements.')

            ! Create matrix product of Q and R matrices as nonlocal psp matrix
            call internal_create_product(qrcchlib,'ccQ'//suffix,'ccR'//suffix)

            ! Augmented overlap 'S' is union of direct sphere-sphere S and QR
            if (pub_aug) then
               do isub=1,mdl%nsub
                  col_par => mdl%regions(isub)%par
                  do jsub=1,mdl%nsub
                     row_par => mdl%regions(jsub)%par
                     call internal_write_embed_struc(reg_struc,reg_struc_t,jsub,isub)

                     ! Count elements in union of O and QR
                     call sparse_count_union(occhlib(jsub,isub), &
                          qrcchlib(jsub,isub),nze,nzb,my_nze,my_nzb, &
                          seg_nze,seg_nzb)

                     ! S structure becomes this union
                     call sparse_index_union(occhlib(jsub,isub), &
                          qrcchlib(jsub,isub),'SccH'//suffix//reg_struc, &
                          nze,nzb,my_nze, &
                          my_nzb,seg_nze,seg_nzb,rlib=scchlib(jsub,isub))

                     ! Calculate filling factor denominator for fraction
!                     fill_fac = sparse_fill_fac_denom(blks,blks)
!
!                     if (pub_on_root.and.show_sparsity) &
!                          write(stdout,'(20x,a,f7.2,a)') 'Aug Overlap matrix filling:', &
!                          fill_fac * nze,'%'

                     ! rc2013: deallocate temporary arrays
                     deallocate(elems_sfc2,stat=ierr)
                     call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
                  end do
                  deallocate(elems_sfc1,stat=ierr)
                  call utils_dealloc_check('sparse_init_rep', 'elems_sfc', ierr)
               end do
            end if

         else
            call utils_assert(.not. pub_aug, 'Error in internal_create_&
                 &representation(): Cannot have no projectors in augmented &
                 &formalisms.')
            qrlib = 0
         end if

      end if


      ! rc2013: loop over the remaining structures
      do isub=1,mdl%nsub
         call internal_allocate_sfc(elems_sfc1, col_par, isub)
         do jsub=1,mdl%nsub
            call internal_allocate_sfc(elems_sfc2, row_par, jsub)
            call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

            ! Set up these matrices only if there are any projectors ! lr408
            if (pub_eels_calculate) then

               !----------------------------------------------------------------
               ! Core-NGWF matrices
               !----------------------------------------------------------------

               ! Obtain a list of overlaps between spherical regions & ion cores
               call internal_create_structure(overlaps(jsub,isub),elems_sfc1, &
                    elems_sfc2,pscode,'Laura',blks,BLKS_CORE,nze_c,nzb, &
                    my_nze,my_nzb,seg_nze,seg_nzb,col_par,row_par, &
                    'C'//suffix//reg_struc,'B'//suffix//reg_struc_t)

               ! Synchronise
               call comms_barrier

               ! Obtain a list of overlaps between ion cores and spherical regions
               call internal_create_structure(overlaps(jsub,isub),elems_sfc1, &
                    elems_sfc2,'Laura',pscode,BLKS_CORE,blks,nze_b,nzb, &
                    my_nze,my_nzb,seg_nze,seg_nzb,col_par,row_par, &
                    'B'//suffix//reg_struc,'C'//suffix//reg_struc_t)

               ! Calculate filling factor denominator for fraction
               fill_fac = sparse_fill_fac_denom(blks,BLKS_CORE,row_par,col_par)

               if (pub_on_root.and.show_sparsity) then
                  if (mdl%nsub.eq.1) then
                     write(stdout,'(14x,a,f7.2,a)') 'NGWF-Core overlap matrix filling:', &
                          fill_fac * nze_c,'%'
                  else
                     write(stdout,'(8x,a,i0,a,i0,a,f7.2,a)') &
                          'NGWF-Core overlap matrix filling (',jsub,',',isub,'):', &
                          fill_fac * nze_c,'%'
                  end if
               end if

               ! Synchronise
               call comms_barrier

               call utils_assert(nze_b == nze_c, 'Error in internal_create_&
                    &representation(): NGWF-Core and Core-NGWF overlap matrices &
                    &must have same number of nonzero elements.')

            end if

            !---------------------------------------------------------------------------
            ! Hubbard Projector-NGWF overlaps
            !---------------------------------------------------------------------------

            ! Create W and V matrices if required
            if (pub_hubbard) then

               ! Obtain a list of overlaps between spherical regions and ion cores
               call internal_create_structure(overlaps(jsub,isub),elems_sfc1, &
                    elems_sfc2,pscode,'Region',blks,BLKS_HUB_PROJ,nze,nzb, &
                    my_nze,my_nzb,seg_nze,seg_nzb,row_par,col_par, &
                    'W'//suffix//reg_struc,'V'//suffix//reg_struc_t)

               ! Obtain a list of overlaps between spherical regions and ion cores
               call internal_create_structure(overlaps(jsub,isub),elems_sfc1, &
                    elems_sfc2,'Region',pscode,BLKS_HUB_PROJ,blks,nze,nzb, &
                    my_nze,my_nzb,seg_nze,seg_nzb,col_par,row_par, &
                    'V'//suffix//reg_struc,'W'//suffix//reg_struc_t)

               ! Calculate filling factor denominator for fraction
               fill_fac = sparse_fill_fac_denom(blks,BLKS_HUB_PROJ,row_par,col_par)

               if (pub_on_root.and.show_sparsity) then
                  if (mdl%nsub.eq.1) then
                     write(stdout,'(6x,a,f7.2,a)') 'NGWF-Hubbard Proj overlap &
                          &matrix filling:', fill_fac * nze,'%'
                  else
                     write(stdout,'(a,i0,a,i0,a,f7.2,a)') &
                          'NGWF-Hubbard Proj overlap matrix filling (',&
                          jsub,',',isub,'):', fill_fac * nze,'%'
                  end if
               end if

            end if
            deallocate(elems_sfc2,stat=ierr)
            call utils_dealloc_check('sparse_init_rep', 'elems_sfc', ierr)
         end do
         deallocate(elems_sfc1,stat=ierr)
         call utils_dealloc_check('sparse_init_rep', 'elems_sfc', ierr)
      end do

      if (pub_hubbard) then
         ! Create matrix product of V and W matrices as Hubbard Hamiltonian
         call internal_create_product(vwlib,'V'//suffix,'W'//suffix)

         ! Calculate filling factor denominator for fraction
         do isub=1,mdl%nsub
            do jsub=1,mdl%nsub

               ! rc2013: assign pointers for row and column parallel strategies
               col_par => mdl%regions(isub)%par
               row_par => mdl%regions(jsub)%par

               fill_fac = sparse_fill_fac_denom(blks,blks,row_par,col_par)
               if (pub_on_root.and.show_sparsity) then
                  if (mdl%nsub.eq.1) then
                     write(stdout,'(29x,a,f7.2,a)') 'VW matrix filling:', &
                          fill_fac * sparse_num_element(vwlib(jsub,isub)),'%'
                  else
                     write(stdout,'(23x,a,i0,a,i0,a,f7.2,a)') &
                          'VW matrix filling (',jsub,',',isub,'):', &
                          fill_fac * sparse_num_element(vwlib(jsub,isub)),'%'
                  end if
               end if
            end do
         end do
      else
         vwlib = 0
      end if

      do isub=1,mdl%nsub
         call internal_allocate_sfc(elems_sfc1, col_par, isub)
         do jsub=1,mdl%nsub
            call internal_allocate_sfc(elems_sfc2, row_par, jsub)
            call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

            !---------------------------------------------------------------------------
            ! Hubbard Projector - Projector overlaps
            !---------------------------------------------------------------------------

            ! ddor: Create Y and X matrices if required, where these are akin
            ! ddor: to W and V, respectively, but with non-local projectors
            ! ddor: replacing NGWFs.
            if (pub_hubbard .and. any_proj) then

               ! Create NGWF-Hubbard Projector overlap matrix
               call internal_create_structure(overlaps(jsub,isub),elems_sfc1, &
                    elems_sfc2,'Core','Region',BLKS_PROJ,BLKS_HUB_PROJ,nze,nzb, &
                    my_nze,my_nzb,seg_nze,seg_nzb,col_par,row_par, &
                    'Y'//suffix//reg_struc,'X'//suffix//reg_struc_t)

               ! Create Hubbard Projector-NGWF overlap matrix
               call internal_create_structure(overlaps(jsub,isub),elems_sfc1, &
                    elems_sfc2,'Region','Core',BLKS_HUB_PROJ,BLKS_PROJ,nze,nzb, &
                    my_nze,my_nzb,seg_nze,seg_nzb,col_par,row_par, &
                    'X'//suffix//reg_struc,'Y'//suffix//reg_struc_t)

            end if

            ! ja531 --> Add pdos SW projection sparsity patterns.

            !---------------------------------------------------------------------------
            ! pDOS : Spherical wave - Spherical wave overlaps
            !---------------------------------------------------------------------------

            if(pub_dos_smear > 0.0_DP .and. pub_pdos_max_l >=0) then

               call internal_create_structure(overlaps(jsub,isub),elems_sfc1, &
                    elems_sfc2,'Region','Region',BLKS_PDOS,BLKS_PDOS, &
                    nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb,col_par,row_par, &
                    'L$'//suffix//reg_struc, &
                    'L$'//suffix//reg_struc_t)

               fill_fac = sparse_fill_fac_denom(BLKS_PDOS,BLKS_PDOS,row_par,col_par)

!               if (pub_on_root.and.show_sparsity) &
!                    write(stdout,'(13x,a,f7.2,a)') 'PDOS SW-SW Overlap matrix filling:', &
!                    fill_fac * nze,'%'

               ! Synchronise
               call comms_barrier

               !---------------------------------------------------------------------------
               ! pDOS : Spherical wave - NGWF overlaps
               !---------------------------------------------------------------------------

               ! Create Projector-NGWF overlap matrix
               call internal_create_structure(overlaps(jsub,isub),elems_sfc1, &
                    elems_sfc2,'Region',pscode,BLKS_PDOS,blks, &
                    nze_r,nzb,my_nze,my_nzb,seg_nze,seg_nzb,col_par,row_par, &
                    'P%'//suffix//reg_struc, &
                    'P$'//suffix//reg_struc_t, rlib=ngwfswlib(jsub,isub))

               ! Synchronise
               call comms_barrier

               ! Create NGWF-Projector overlap matrix
               call internal_create_structure(overlaps(jsub,isub),elems_sfc1, &
                    elems_sfc2,pscode,'Region',BLKS,BLKS_PDOS, &
                    nze_q,nzb,my_nze,my_nzb,seg_nze,seg_nzb,col_par,row_par, &
                    'P$'//suffix//reg_struc, &
                    'P%'//suffix//reg_struc_t,rlib=swngwflib(jsub,isub))

               ! Calculate filling factor denominator for fraction
               fill_fac = sparse_fill_fac_denom(BLKS_PDOS,BLKS,row_par,col_par)

!               if (pub_on_root.and.show_sparsity) &
!                    write(stdout,'(11x,a,f7.2,a)') 'PDOS SW-NGWF overlap matrix &
!                    &filling:', fill_fac * nze_r,'%'


               ! Set up these matrices only if there are any projectors
               if (any_proj) then

                  ! Create Projector-SW overlap matrix
                  call internal_create_structure(overlaps(jsub,isub),elems_sfc1, &
                       elems_sfc2,'Region','Core',BLKS_PDOS,BLKS_PROJ, &
                       nze_r,nzb,my_nze,my_nzb,seg_nze,seg_nzb,col_par,row_par, &
                       'Pr%'//suffix//reg_struc, &
                       'Pr$'//suffix//reg_struc_t)

                  ! Synchronise
                  call comms_barrier

                  ! Create NGWF-Projector overlap matrix
                  call internal_create_structure(overlaps(jsub,isub),elems_sfc1, &
                       elems_sfc2,'Core','Region',BLKS_PROJ,BLKS_PDOS, &
                       nze_q,nzb,my_nze,my_nzb,seg_nze,seg_nzb,col_par,row_par, &
                       'Pr$'//suffix//reg_struc, &
                       'Pr%'//suffix//reg_struc_t)

                  ! Calculate filling factor denominator for fraction
                  fill_fac = sparse_fill_fac_denom(BLKS_PDOS,BLKS_PROJ,row_par,col_par)

!                  if (pub_on_root.and.show_sparsity) &
!                       write(stdout,'(14x,a,f7.2,a)') 'NGWF-Proj overlap matrix &
!                       &filling:', fill_fac * nze_r,'%'

                  ! Synchronise
                  call comms_barrier

                  ! rc2013: don't see why this'd be the case with cross overlaps?
                  if(isub == jsub) call utils_assert(nze_q == nze_r, &
                       'Error in internal_create_&
                       &representation(): SW-Projector and Projector-SW overlap &
                       &matrices must have same number of nonzero elements.')

               end if
            end if
            deallocate(elems_sfc2,stat=ierr)
            call utils_dealloc_check('sparse_init_rep', 'elems_sfc', ierr)
         end do
         deallocate(elems_sfc1,stat=ierr)
         call utils_dealloc_check('sparse_init_rep', 'elems_sfc', ierr)
      end do


      if(pub_dos_smear > 0.0_DP .and. pub_pdos_max_l >=0) then

         ! Set up these matrices only if there are any projectors
         if (any_proj) then

            ! Create matrix product of Q and R matrices as nonlocal psp matrix
!            call internal_create_product(prqlib,'Pr$','R'//suffix)
!            call internal_create_product(qprlib,'Q'//suffix,'Pr%')

!write(*,*) pub_my_proc_id,"prqlib"
            call internal_create_product(prqlib,'Pr$'//suffix,'R'//suffix) !! 117, 156
!write(*,*) pub_my_proc_id,"qprlib"
            call internal_create_product(qprlib,'Q'//suffix,'Pr%'//suffix)             !! 156, 117

            ! Augmented overlap 'S' is union of direct sphere-sphere S and QR
            if (pub_aug) then

               do isub=1,mdl%nsub
                  col_par => mdl%regions(isub)%par
                  do jsub=1,mdl%nsub
                     row_par => mdl%regions(jsub)%par
                     call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

                     ! Count elements in union of O and QR
                     call sparse_count_union(prqlib(jsub,isub), &
                          swngwflib(jsub,isub),nze,nzb,my_nze,my_nzb, &
                          seg_nze,seg_nzb)
                     ! S structure becomes this union
                     call sparse_index_union(prqlib(jsub,isub), &
                          swngwflib(jsub,isub),'P$'//suffix//reg_struc, &
                          nze,nzb,my_nze, & !117 156
                          my_nzb,seg_nze,seg_nzb)!,rlib=slib)

                     ! Calculate filling factor denominator for fraction
                     fill_fac = sparse_fill_fac_denom(blks,blks,row_par,col_par)

!                     if (pub_on_root.and.show_sparsity) &
!                          write(stdout,'(20x,a,f7.2,a)') 'Aug Overlap matrix filling:', &
!                          fill_fac * nze,'%'

                     ! Count elements in union of O and QR
                     call sparse_count_union(ngwfswlib(jsub,isub), &
                          qprlib(jsub,isub),nze,nzb,my_nze,my_nzb, &
                          seg_nze,seg_nzb)
                     ! S structure becomes this union
                     call sparse_index_union(ngwfswlib(jsub,isub), &
                          qprlib(jsub,isub),'P%'//suffix//reg_struc, &
                          nze,nzb,my_nze, &
                          my_nzb,seg_nze,seg_nzb)!,rlib=slib)

                     ! Calculate filling factor denominator for fraction
                     fill_fac = sparse_fill_fac_denom(blks,blks,row_par,col_par)


!                     if (pub_on_root.and.show_sparsity) &
!                          write(stdout,'(20x,a,f7.2,a)') 'Aug Overlap matrix filling:', &
!                          fill_fac * nze,'%'
                  end do
               end do


            end if
         end if


      end if

    end subroutine internal_create_representation


    !==========================================================================!
    ! This subroutine creates the structures of the Hamiltonian, as the union  !
    ! of various structures that might contribute nonzero elements.            !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in April 2011 based on existing code.           !
    !==========================================================================!

    subroutine internal_create_hamiltonian

      implicit none

      integer :: klib_tmp(mdl%nsub,mdl%nsub)

      ! Set up Contra-Covariant Ham pattern
      do isub=1,mdl%nsub
         col_par => mdl%regions(isub)%par
         do jsub=1,mdl%nsub
            row_par => mdl%regions(jsub)%par
            call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

            if(abs(pub_contracoham_radmult-1.0_dp)>epsilon(1.0_dp)) then
               ! Start by copying the overlap matrix into Hamiltonian
               call sparse_count_union(scchlib(jsub,isub),scchlib(jsub,isub), &
                    nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
               call sparse_index_union(scchlib(jsub,isub),scchlib(jsub,isub), &
                    'ccH'//suffix//reg_struc,nze,nzb, &
                    my_nze,my_nzb,seg_nze,seg_nzb,rlib=cchlib(jsub,isub))

               if (any_proj) then

                  ! Count elements in union of H and QR
                  call sparse_count_union(cchlib(jsub,isub),qrlib(jsub,isub), &
                       nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
                  ! Hamiltonian structure becomes this union
                  call sparse_index_union(cchlib(jsub,isub),qrlib(jsub,isub), &
                       'ccH'//suffix//reg_struc,nze,nzb, &
                       my_nze,my_nzb,seg_nze,seg_nzb)
               end if

               if (pub_hubbard) then

                  ! Count elements in union of H and VW
                  call sparse_count_union(cchlib(jsub,isub),vwlib(jsub,isub), &
                       nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
                  ! Hamiltonian structure becomes this union
                  call sparse_index_union(cchlib(jsub,isub),vwlib(jsub,isub), &
                       'ccH'//suffix//reg_struc,nze,nzb, &
                       my_nze,my_nzb,seg_nze,seg_nzb)

               end if

               if (pub_use_hfx .and. .not. pub_hfx_nlpp_for_exchange) then

                  ! Count elements in union of H and KSKSK
                  call sparse_count_union(cchlib(jsub,isub),ksksklib(jsub,isub), &
                       nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
                  ! Hamiltonian structure becomes this union
                  call sparse_index_union(cchlib(jsub,isub),ksksklib(jsub,isub), &
                       'ccH'//suffix//reg_struc,nze,nzb, &
                       my_nze,my_nzb,seg_nze,seg_nzb)
               else if(pub_use_activehfx .and. .not. pub_hfx_nlpp_for_exchange) then
                  ! jcap: if we only need HF in the active subregion, only
                  ! do this for the active subregion (exchange within
                  ! subregion only)
                  if ((isub==jsub) .and. (isub==pub_active_region)) then

                     ! Count elements in union of H and KSKSK
                     call sparse_count_union(cchlib(jsub,isub),ksksklib(jsub,isub), &
                          nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
                     ! Hamiltonian structure becomes this union
                     call sparse_index_union(cchlib(jsub,isub),ksksklib(jsub,isub), &
                          'ccH'//suffix//reg_struc, &
                          nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
                  end if
               end if

               ! Calculate filling factor and display
               fill_fac = sparse_fill_fac_denom(blks,blks,row_par,col_par)
               if (pub_on_root.and.show_sparsity) then
                  if (mdl%nsub.eq.1) then
                     write(stdout,'(20x,a,f7.2,a)') 'Contra-Covariant &
                          &Hamiltonian matrix filling:', fill_fac * nze,'%'
                  else
                     write(stdout,'(14x,a,i0,a,i0,a,f7.2,a)') &
                          'Contra-Covariant Hamiltonian matrix filling (',&
                          jsub,',',isub,'):', fill_fac * nze,'%'
                  end if
               end if

            end if

            ! And now standard Ham:

            ! Start by copying the overlap matrix into Hamiltonian
            call sparse_count_union(slib(jsub,isub),slib(jsub,isub), &
                 nze,nzb,my_nze,my_nzb, seg_nze,seg_nzb)
            call sparse_index_union(slib(jsub,isub),slib(jsub,isub), &
                 'H'//suffix//reg_struc, &
                 nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb,rlib=hlib(jsub,isub))

            if (any_proj) then

               ! Count elements in union of H and QR
               call sparse_count_union(hlib(jsub,isub),qrlib(jsub,isub), &
                    nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
               ! Hamiltonian structure becomes this union
               call sparse_index_union(hlib(jsub,isub),qrlib(jsub,isub), &
                    'H'//suffix//reg_struc, &
                    nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
            end if

            if (pub_hubbard) then

               ! Count elements in union of H and VW
               call sparse_count_union(hlib(jsub,isub),vwlib(jsub,isub), &
                    nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
               ! Hamiltonian structure becomes this union
               call sparse_index_union(hlib(jsub,isub),vwlib(jsub,isub), &
                    'H'//suffix//reg_struc, &
                    nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)

            end if

            if (pub_use_hfx .and. .not. pub_hfx_nlpp_for_exchange) then

               ! Count elements in union of H and KSKSK
               call sparse_count_union(hlib(jsub,isub),ksksklib(jsub,isub), &
                    nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
               ! Hamiltonian structure becomes this union
               call sparse_index_union(hlib(jsub,isub),ksksklib(jsub,isub), &
                    'H'//suffix//reg_struc, &
                    nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
            else if(pub_use_activehfx .and. .not. pub_hfx_nlpp_for_exchange) then
               ! jcap: if we only need HF in the active subregion, only
               ! do this for the active subregion (exchange within
               ! subregion only)
               if ((isub==jsub) .and. (isub==pub_active_region)) then

                  ! Count elements in union of H and KSKSK
                  call sparse_count_union(hlib(jsub,isub),ksksklib(jsub,isub), &
                       nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
                  ! Hamiltonian structure becomes this union
                  call sparse_index_union(hlib(jsub,isub),ksksklib(jsub,isub), &
                       'H'//suffix//reg_struc, &
                       nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
               end if
            end if

            ! Calculate filling factor and display
            fill_fac = sparse_fill_fac_denom(blks,blks,row_par,col_par)
            if (pub_on_root.and.show_sparsity) then
               if (mdl%nsub.eq.1) then
                  write(stdout,'(20x,a,f7.2,a)') 'Hamiltonian matrix filling:', &
                       fill_fac * nze,'%'
               else
                  write(stdout,'(14x,a,i0,a,i0,a,f7.2,a)') &
                       'Hamiltonian matrix filling (',jsub,',',isub,'):', &
                       fill_fac * nze,'%'
               end if
            end if

         end do
      end do


      ! ja531-> If H2denskern sparsity denskern -> H^2
      if(pub_H2denskern_sparsity) then

         if(abs(pub_contracoham_radmult-1.0_dp)>epsilon(1.0_dp)) then
            call internal_create_product(klib,'ccH'//suffix,'ccH'//suffix)
         else
            call internal_create_product(klib,'H'//suffix,'H'//suffix)
         end if


         do isub=1,mdl%nsub
            col_par => mdl%regions(isub)%par
            do jsub=1,mdl%nsub
               row_par => mdl%regions(jsub)%par

               call sparse_count_union(klib(jsub,isub),klib(jsub,isub), &
                    nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)

               call sparse_index_union(klib(jsub,isub),klib(jsub,isub), &
                    'K'//suffix//reg_struc,nze,nzb,my_nze, &
                    my_nzb,seg_nze,seg_nzb,rlib=klib_tmp(jsub,isub))

               klib(jsub,isub) = klib_tmp(jsub,isub)

               fill_fac = sparse_fill_fac_denom(blks,blks,row_par,col_par)
               if (pub_on_root.and.show_sparsity) then
                  if (mdl%nsub.eq.1) then
                     write(stdout,'(24x,a,f7.2,a)') 'Density kernel filling:', &
                          fill_fac * nze,'%'
                  else
                     write(stdout,'(18x,a,i0,a,i0,a,f7.2,a)') &
                          'Density kernel filling (',jsub,',',isub,'):', &
                          fill_fac * nze,'%'
                  end if
               end if
            end do
         end do

      end if


    end subroutine internal_create_hamiltonian


    !==========================================================================!
    ! This subroutine creates a structure by counting overlaps for a given     !
    ! pair of region types and a pair of block lists.                          !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in April 2011 based on existing code.           !
    ! Modified by Robert Charlton, 25/04/2018, to handle multiple regions.     !
    !==========================================================================!

    subroutine internal_create_structure(overlaps,elems_sfc1,elems_sfc2, &
         row_ovlp_type,col_ovlp_type,col_blks,row_blks,nze,nzb,my_nze,my_nzb, &
         seg_nze,seg_nzb,col_par,row_par,name,transpose_name, &
         radius,rlib,tddft_region,multiplier)

      implicit none

      ! Arguments
      type(OVERLAP_LIST), intent(inout) :: overlaps
      integer, intent(in) :: col_blks
      integer, intent(in) :: row_blks
      character, intent(in) :: col_ovlp_type
      character, intent(in) :: row_ovlp_type
      character(*), intent(in) :: name
      integer(kind=LONG), intent(out) :: nze
      integer, intent(out) :: nzb
      integer, intent(out) :: my_nze, my_nzb
      integer, intent(out) :: seg_nze(0:pub_total_num_procs)
      integer, intent(out) :: seg_nzb(0:pub_total_num_procs)
      character(*), intent(in) :: transpose_name
      type(PARAL_INFO), pointer, intent(in) :: col_par
      type(PARAL_INFO), pointer, intent(in) :: row_par
      type(ELEMENT), intent(in) :: elems_sfc1(col_par%nat)
      type(ELEMENT), intent(in) :: elems_sfc2(row_par%nat)
      real(kind=DP), intent(in), optional :: radius
      integer, intent(out), optional :: rlib
      logical, optional :: tddft_region
      real(kind=dp), intent(in), optional :: multiplier

      ! Obtain a list of overlaps between spherical regions
      call parallel_strategy_list_cross_overlaps(overlaps, elems_sfc1, &
           elems_sfc2, mdl%cell, row_ovlp_type, col_ovlp_type, &
           radius, tddft_region, multiplier)

      ! Count number of corresponding nonzero NGWF matrix elements
      call sparse_count_ss(overlaps,col_blks,row_blks,nze,nzb,my_nze,my_nzb, &
           seg_nze, seg_nzb, row_par, col_par)

      ! Create sparse overlap matrix
      call sparse_index_ss(overlaps,col_blks,row_blks,name,nze,nzb,my_nze, &
           my_nzb, seg_nze, seg_nzb, transpose_name=transpose_name, &
           row_par=row_par, col_par=col_par, rlib=rlib)

    end subroutine internal_create_structure


    subroutine internal_create_cross_overlap(blks1,blks2,suffix1,suffix2, &
         pscode1,pscode2)

      ! Arguments
      integer, intent(in) :: blks1,blks2
      character(len=*), intent(in) :: suffix1,suffix2
      character(len=*), intent(in) :: pscode1,pscode2

      ! Local variables
      character(len=30) :: ustruc, tstruc

      !----------------------------------------------------------
      ! Cross overlap matrices (eg between Valence and COND NGWFs
      !----------------------------------------------------------

      do isub=1,mdl%nsub
         call internal_allocate_sfc(elems_sfc1, col_par, isub)
         do jsub=1,mdl%nsub
            call internal_allocate_sfc(elems_sfc2, row_par, jsub)
            call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

            if (.not.pub_aug) then
               ustruc = 'U'//suffix2//suffix1
               tstruc = 'T'//suffix1//suffix2
            else
               ustruc = 'bU'//suffix2//suffix1
               tstruc = 'bT'//suffix1//suffix2
            end if

            ! Obtain a list of overlaps between NGWF set 1 and set 2 regions
            call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub), &
                 elems_sfc1, elems_sfc2, mdl%cell, pscode1, pscode2)

            ! Count number of corresponding nonzero matrix elements
            call sparse_count_ss(overlaps(jsub,isub),blks1,blks2,nze,nzb,my_nze,my_nzb, &
                 seg_nze,seg_nzb,row_par,col_par)

            ! Create bare NGWF overlap matrix
            call sparse_index_ss(overlaps(jsub,isub),blks1,blks2, &
                 trim(ustruc)//reg_struc,nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb, &
                 transpose_name=trim(tstruc)//reg_struc_t,rlib=ulib(jsub,isub), &
                 row_par=row_par,col_par=col_par)
            if (pub_aug) bulib(jsub,isub) = ulib(jsub,isub)

            ! Obtain a list of overlaps between NGWF set 2 and set 1 regions
            call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub), &
                 elems_sfc1, elems_sfc2, mdl%cell, pscode2, pscode1)

            ! Count number of corresponding nonzero matrix elements
            call sparse_count_ss(overlaps(jsub,isub),blks2,blks1,nze,nzb,my_nze,my_nzb, &
                 seg_nze,seg_nzb,row_par,col_par)

            ! Create bare NGWF overlap matrix
            call sparse_index_ss(overlaps(jsub,isub),blks2,blks1, &
                 trim(tstruc)//reg_struc,nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb, &
                 transpose_name=trim(ustruc)//reg_struc_t,rlib=tlib(jsub,isub), &
                 row_par=row_par,col_par=col_par)
            if (pub_aug) btlib = tlib

            ! Calculate filling factor denominator for fraction
            fill_fac = sparse_fill_fac_denom(blks1,blks2,row_par,col_par)

            if (pub_on_root.and.show_sparsity) then
               if (mdl%nsub.eq.1) then
                  write(stdout,'(14x,a,f7.2,a)') 'NGWF-COND overlap matrix filling:', &
                       fill_fac * nze,'%'
               else
                  write(stdout,'(8x,a,i0,a,i0,a,f7.2,a)') &
                       'NGWF-COND overlap matrix filling (',jsub,',',isub,'):', &
                       fill_fac * nze,'%'
               end if
            end if

            ! jcap: check if we need to do this
            if (any_proj) then
               ! Create matrix product of Q and R matrices as nonlocal psp matrix
               call internal_create_product(qcrlib,'Q'//suffix2,'R'//suffix1)

               ! Create matrix product of Q and R matrices as nonlocal psp matrix
               call internal_create_product(qrclib,'Q'//suffix1,'R'//suffix2)
            end if

            ! Augmented cross overlap 'T' is union of direct sphere-sphere T and QR
            if (pub_aug) then

               ustruc = 'U'//suffix2//suffix1
               tstruc = 'T'//suffix1//suffix2

               ! Count elements in union of U and RQ
               call sparse_count_union(bulib(jsub,isub),qcrlib(jsub,isub), &
                    nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)

               ! U structure becomes this union
               call sparse_index_union(bulib(jsub,isub),qcrlib(jsub,isub), &
                    trim(ustruc)//reg_struc,nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb, &
                    rlib=ulib(jsub,isub), transpose_name=trim(tstruc)//reg_struc_t)

               ! Calculate filling factor denominator for fraction
               fill_fac = sparse_fill_fac_denom(blks1,blks2,row_par,col_par)

               if (pub_on_root.and.show_sparsity) then
                  if (mdl%nsub.eq.1) then
                     write(stdout,'(9x,a,f7.2,a)') 'Aug Trans Cross Overlap &
                          &matrix filling:', fill_fac * nze,'%'
                  else
                     write(stdout,'(3x,a,i0,a,i0,a,f7.2,a)') &
                          'Aug Trans Cross Overlap matrix filling (',&
                          jsub,',',isub,'):', fill_fac * nze,'%'
                  end if
               end if

               ! Count elements in union of T and QR
               call sparse_count_union(btlib(jsub,isub),qrclib(jsub,isub), &
                    nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)

               ! T structure becomes this union
               call sparse_index_union(btlib(jsub,isub),qrclib(jsub,isub), &
                    trim(tstruc)//reg_struc,nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb, &
                    rlib=tlib(jsub,isub),transpose_name=trim(ustruc)//reg_struc_t)

               ! Calculate filling factor denominator for fraction
               fill_fac = sparse_fill_fac_denom(blks2,blks1,row_par,col_par)

               if (pub_on_root.and.show_sparsity) then
                  if (mdl%nsub.eq.1) then
                     write(stdout,'(14x,a,f7.2,a)') 'Aug Cross Overlap &
                          &matrix filling:', fill_fac * nze,'%'
                  else
                     write(stdout,'(8x,a,i0,a,i0,a,f7.2,a)') &
                          'Aug Cross Overlap matrix filling (',&
                          jsub,',',isub,'):', fill_fac * nze,'%'
                  end if
               end if
            end if

            deallocate(elems_sfc2,stat=ierr)
            call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
         end do
         deallocate(elems_sfc1,stat=ierr)
         call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
      end do

    end subroutine internal_create_cross_overlap

    !==========================================================================!
    ! This subroutine creates the structures of various matrices associated    !
    ! with the calculation of conduction NGWFs                                 !
    !--------------------------------------------------------------------------!
    ! Created by Nicholas Hine in April 2011, based on code by Laura Ratcliff, !
    ! from June 2010.                                                          !
    !==========================================================================!

    subroutine internal_create_cond_structures

      ! Local variables
      character(len=30) :: ustruc, tstruc

#if 1
      call internal_create_cross_overlap(BLKS_NGWF,blks,'',suffix, &
           'Region',pscode)

      ! Structures for overlap between joint and cond NGWFs
      if (suffix=='j') call internal_create_cross_overlap(BLKS_COND,blks, &
           'c',suffix,'A',pscode)
#else

      ! Create NGWF-COND overlap matrix
      if (.not.pub_aug) then
         ustruc = 'U'//suffix
         tstruc = 'T'//suffix
      else
         ustruc = 'bU'//suffix
         tstruc = 'bT'//suffix
      end if

      do isub=1,mdl%nsub
         call internal_allocate_sfc(elems_sfc1, col_par, isub)
         do jsub=1,mdl%nsub
            call internal_allocate_sfc(elems_sfc2, row_par, jsub)
            call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

            !------------------------------------------
            ! Cross overlap matrices with Valence NGWFs
            !------------------------------------------

            ! Create NGWF-COND overlap matrix
            call internal_create_structure(overlaps(jsub,isub),elems_sfc1,elems_sfc2, &
                 'Region', pscode, BLKS_NGWF, blks, nze, nzb, my_nze, my_nzb, &
                 seg_nze, seg_nzb, col_par, row_par, trim(ustruc)//reg_struc, &
                 trim(tstruc)//reg_struc_t, rlib=ulib(jsub,isub))
            if (pub_aug) bulib(jsub,isub) = ulib(jsub,isub)

            ! Create COND-NGWF blks matrix
            call internal_create_structure(overlaps(jsub,isub),elems_sfc1,elems_sfc2, &
                 pscode, 'Region', blks, BLKS_NGWF, nze, nzb, my_nze, my_nzb, &
                 seg_nze, seg_nzb, col_par, row_par, trim(tstruc)//reg_struc, &
                 trim(ustruc)//reg_struc_t, rlib=tlib(jsub,isub))

            ! Calculate filling factor denominator for fraction
            fill_fac = sparse_fill_fac_denom(BLKS_NGWF,blks,row_par,col_par)

            if (pub_on_root.and.show_sparsity) then
               if (mdl%nsub.eq.1) then
                  write(stdout,'(14x,a,f7.2,a)') 'NGWF-COND overlap &
                       &matrix filling:', fill_fac * nze,'%'
               else
                  write(stdout,'(8x,a,i0,a,i0,a,f7.2,a)') &
                       'NGWF-COND overlap matrix filling (',&
                       jsub,',',isub,'):', fill_fac * nze,'%'
               end if
            end if

            deallocate(elems_sfc2,stat=ierr)
            call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
         end do
         deallocate(elems_sfc1,stat=ierr)
         call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
      end do

      ! jcap: check if we need to do this
      if (any_proj) then
         ! Create matrix product of Q and R matrices as nonlocal psp matrix
         call internal_create_product(qcrlib,'Q'//suffix,'R')

         ! Create matrix product of Q and R matrices as nonlocal psp matrix
         call internal_create_product(qrclib,'Q','R'//suffix)
      end if

      ! Augmented cross overlap 'T' is union of direct sphere-sphere T and QR
      if (pub_aug) then

         ustruc = 'U'//suffix
         tstruc = 'T'//suffix

         do isub=1,mdl%nsub
            call internal_allocate_sfc(elems_sfc1, col_par, isub)
            do jsub=1,mdl%nsub
               call internal_allocate_sfc(elems_sfc2, row_par, jsub)
               call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

               ! Count elements in union of U and RQ
               call sparse_count_union(bulib(jsub,isub),qcrlib(jsub,isub), &
                    nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)

               ! U structure becomes this union
               call sparse_index_union(bulib(jsub,isub),qcrlib(jsub,isub), &
                    trim(ustruc)//reg_struc,nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb, &
                    rlib=ulib(jsub,isub),transpose_name=trim(tstruc)//reg_struc_t)

               ! Calculate filling factor denominator for fraction
               fill_fac = sparse_fill_fac_denom(BLKS_NGWF,blks,row_par,col_par)

               if (pub_on_root.and.show_sparsity) then
                  if (mdl%nsub.eq.1) then
                     write(stdout,'(8x,a,f7.2,a)') 'Aug Trans Cross Overlap &
                          &matrix filling:', fill_fac * nze,'%'
                  else
                     write(stdout,'(2x,a,i0,a,i0,a,f7.2,a)') &
                          'Aug Trans Cross Overlap matrix filling (',&
                          jsub,',',isub,'):', fill_fac * nze,'%'
                  end if
               end if

               ! Count elements in union of T and QR
               call sparse_count_union(btlib(jsub,isub),qrclib(jsub,isub), &
                    nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)

               ! T structure becomes this union
               call sparse_index_union(btlib(jsub,isub),qrclib(jsub,isub), &
                    trim(tstruc)//reg_struc,nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb, &
                    rlib=tlib(jsub,isub),transpose_name=trim(ustruc)//reg_struc_t)

               ! Calculate filling factor denominator for fraction
               fill_fac = sparse_fill_fac_denom(blks,BLKS_NGWF,row_par,col_par)

               if (pub_on_root.and.show_sparsity) then
                  if (mdl%nsub.eq.1) then
                     write(stdout,'(14x,a,f7.2,a)') 'Aug Cross Overlap &
                          &matrix filling:', fill_fac * nze,'%'
                  else
                     write(stdout,'(8x,a,i0,a,i0,a,f7.2,a)') &
                          'Aug Cross Overlap matrix filling (',&
                          jsub,',',isub,'):', fill_fac * nze,'%'
                  end if
               end if
               deallocate(elems_sfc2,stat=ierr)
               call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
            end do
            deallocate(elems_sfc1,stat=ierr)
            call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
         end do
      end if

      ! Structures for overlap between joint and cond NGWFs
      if (suffix=='j') then
         do isub=1,mdl%nsub
            call internal_allocate_sfc(elems_sfc1, col_par, isub)
            do jsub=1,mdl%nsub
               call internal_allocate_sfc(elems_sfc2, row_par, jsub)
               call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

               ! Create JOINT-COND overlap matrix
               call internal_create_structure(overlaps(jsub,isub),elems_sfc1,elems_sfc2, &
                    pscode, 'A', blks, BLKS_COND, nze, nzb, my_nze, my_nzb, &
                    seg_nze, seg_nzb, col_par, row_par, &
                    'T'//'c'//trim(suffix)//reg_struc, &
                    'U'//trim(suffix)//'c'//reg_struc_t)

               ! Create COND-JOINT overlap matrix
               call internal_create_structure(overlaps(jsub,isub),elems_sfc1,elems_sfc2, &
                    'A', pscode, BLKS_COND, blks, nze, nzb, my_nze, my_nzb, &
                    seg_nze, seg_nzb, col_par, row_par, &
                    'U'//trim(suffix)//'c'//reg_struc, &
                    'T'//'c'//trim(suffix)//reg_struc_t)
               deallocate(elems_sfc2,stat=ierr)
               call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
            end do
            deallocate(elems_sfc1,stat=ierr)
            call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
         end do
      end if
#endif

      ! Rest is not relevant to joint matrices - only during COND optimisation
      if (suffix/='c') return

      !-------------------------------
      ! Conduction-projector overlaps
      !-------------------------------

      ! Create COND-PROJ and PROJ-COND overlaps
      if (any_proj) then

         ! Create matrix product of Q and R matrices as nonlocal psp matrix
         !call internal_create_product(tlib,'Qc','R')

         ! Create matrix product of Q and R matrices as nonlocal psp matrix
         !call internal_create_product(tlib,'Q','Rc')
      end if

      call internal_cond_ham_products

      !---------------------------------------------------------------------------
      ! Projected conduction Hamiltonian
      !---------------------------------------------------------------------------

      do isub=1,mdl%nsub
       col_par => mdl%regions(isub)%par
         do jsub=1,mdl%nsub
            row_par => mdl%regions(jsub)%par
            call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

            call sparse_count_union(cond_proj_ham_lib(jsub,isub),hlib(jsub,isub), &
                 nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)

            call sparse_index_union(cond_proj_ham_lib(jsub,isub),hlib(jsub,isub), &
                 'L'//suffix//reg_struc,nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)

            fill_fac = sparse_fill_fac_denom(BLKS_COND,BLKS_COND,row_par,col_par)

            ! Find the number of nonzero elements in the patterns
            ! of commonly-used matrices and print results as part of filling
            ! factor summary
            if (pub_on_root.and.show_sparsity) then
               if (mdl%nsub.eq.1) then
                  write(stdout,'(8x,a,f7.2,a)') 'Proj. Cond. Hamiltonian &
                       &matrix filling:', fill_fac * nze,'%'
               else
                  write(stdout,'(2x,a,i0,a,i0,a,f7.2,a)') &
                       'Proj. Cond. Hamiltonian matrix filling (',&
                       jsub,',',isub,'):', fill_fac * nze,'%'
               end if
            end if

            ! Gradient term !k...k_cond sparsity
            call sparse_count_union(cond_grad_lib1(jsub,isub),cond_grad_lib2(jsub,isub), &
                 nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)

            call sparse_index_union(cond_grad_lib1(jsub,isub),cond_grad_lib2(jsub,isub), &
                 'M'//suffix//reg_struc,nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)

         end do ! jsub
      end do ! isub

      call internal_cond_extra_products

    end subroutine internal_create_cond_structures

    !==========================================================================!
    ! This subroutine creates the structures of various matrices associated    !
    ! with the calculation of lr_phonons response NGWFs                        !
    !--------------------------------------------------------------------------!
    ! Adapted for embedding structures by Robert Charlton, 14/08/2018.         !
    !==========================================================================!

    subroutine internal_create_resp_structures

      do isub=1,mdl%nsub
         call internal_allocate_sfc(elems_sfc1, col_par, isub)
         do jsub=1,mdl%nsub
            call internal_allocate_sfc(elems_sfc2, row_par, jsub)
            call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

            !------------------------------------------
            ! Cross overlap matrices with Valence NGWFs
            !------------------------------------------

            ! Create NGWF-RESP overlap matrix
            call internal_create_structure(overlaps(jsub,isub),elems_sfc1,elems_sfc2, &
                 'Region', pscode, BLKS_NGWF, blks, nze, nzb, my_nze, my_nzb, &
                 seg_nze, seg_nzb, col_par, row_par, &
                 'U'//suffix//reg_struc, 'T'//suffix//reg_struc_t)

            ! Create RESP-NGWF blks matrix
            call internal_create_structure(overlaps(jsub,isub),elems_sfc1,elems_sfc2, &
                 pscode, 'Region', blks, BLKS_NGWF, nze, nzb, my_nze, my_nzb, &
                 seg_nze, seg_nzb, col_par, row_par, &
                 'T'//suffix//reg_struc, 'U'//suffix//reg_struc_t)

            ! Calculate filling factor denominator for fraction
            fill_fac = sparse_fill_fac_denom(BLKS_NGWF,blks,row_par,col_par)

            if (pub_on_root.and.show_sparsity) then
               if (mdl%nsub.eq.1) then
                  write(stdout,'(14x,a,f7.2,a)') 'NGWF-RESP overlap &
                       &matrix filling:', fill_fac * nze,'%'
               else
                  write(stdout,'(8x,a,i0,a,i0,a,f7.2,a)') &
                       'NGWF-RESP overlap matrix filling (',&
                       jsub,',',isub,'):', fill_fac * nze,'%'
               end if
            end if

   !        ! Structures for overlap between joint val-resp and resp NGWFs
   !        if (suffix=='z') then

   !           ! Obtain a list of overlaps between RESP NGWF region and current region
   !           call parallel_strategy_list_overlaps(overlaps,par%nat,elems_sfc, &
   !                cell,'A',pscode)

   !           ! Count number of corresponding nonzero matrix elements
   !           call sparse_count_ss(overlaps,BLKS_NGWF,blks,nze,nzb,my_nze,my_nzb, &
   !                seg_nze,seg_nzb, col_par=par)

   !           ! Create NGWF-RESP overlap matrix
   !           call sparse_index_ss(overlaps,BLKS_NGWF,blks,'U'//trim(suffix)//'c', &
   !                nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb, &
   !                transpose_name='T'//'c'//trim(suffix), row_par=par, col_par=par)

   !           ! Obtain a list of overlaps between current region and conduction
   !           ! NGWF region
   !           call parallel_strategy_list_overlaps(overlaps,par%nat,elems_sfc, &
   !                cell,pscode,'A')

   !           ! Count number of corresponding nonzero matrix elements
   !           call sparse_count_ss(overlaps,blks,BLKS_NGWF,nze,nzb,my_nze,my_nzb, &
   !                seg_nze,seg_nzb, col_par=par)

   !           ! Create RESP-NGWF blks matrix
   !           call sparse_index_ss(overlaps,blks,BLKS_NGWF,'T'//'c'//trim(suffix), &
   !                nze, nzb,my_nze,my_nzb,seg_nze,seg_nzb, &
   !                transpose_name='U'//trim(suffix)//'c', row_par=par, col_par=par)

   !        end if

            deallocate(elems_sfc2,stat=ierr)
            call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
         end do
         deallocate(elems_sfc1,stat=ierr)
         call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
      end do

    end subroutine internal_create_resp_structures


    !==========================================================================!
    ! This subroutine creates filling structures for the TDDFT calculation,    !
    ! most notably the effective response density kernel.                      !
    !--------------------------------------------------------------------------!
    ! Written by Tim Zuehlsdorff, November 2012                                !
    ! Copied and modified for embedding by Joseph Prentice, August 2018        !
    !==========================================================================!

    subroutine internal_create_tddft_structures

      ! jcap: loop over regions
      do isub=1,mdl%nsub
         call internal_allocate_sfc(elems_sfc1, col_par, isub)
         do jsub=1,mdl%nsub
            call internal_allocate_sfc(elems_sfc2, row_par, jsub)
            call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

            !--------------------------------------
            ! Response density kernel
            !--------------------------------------
           ! determine the overlaps of fixed cutoff. Might need to change kernel
           ! cutoff to a special response kernel cutoff. Include the option of a
           ! cutoff region rather than a fixed radius
            if(pub_lr_tddft_sparse_region) then
               call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub),&
                    elems_sfc1,elems_sfc2,mdl%cell,'Fixed','Fixed',&
                    0.5_DP * pub_lr_tddft_kernel_cutoff,.true.)
            else
               call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub),&
                    elems_sfc1,elems_sfc2,mdl%cell,'Fixed','Fixed',&
                    0.5_DP * pub_lr_tddft_kernel_cutoff)
            endif

            ! counts the number of non-zero matrix elements between cond and val
            ! ngwfs. Also works if the joint matrix is used to represent the
            ! conduction space
            call sparse_count_ss(overlaps(jsub,isub),BLKS_NGWF,blks,nze,nzb,&
                 my_nze,my_nzb,seg_nze,seg_nzb,row_par=row_par,col_par=col_par)

            ! create effective response density kernel
            call sparse_index_ss(overlaps(jsub,isub),BLKS_NGWF,blks,&
                 'TDRA'//suffix//reg_struc,nze,nzb,&
                 my_nze,my_nzb,seg_nze,seg_nzb,&
                 transpose_name='TDRAt'//suffix//reg_struc_t,&
                 row_par=row_par,col_par=col_par)
            ! allocate the transpose as well
            call sparse_count_ss(overlaps(jsub,isub),blks,BLKS_NGWF,nze,nzb,&
                 my_nze,my_nzb,seg_nze,seg_nzb,row_par=row_par,col_par=col_par)

            call sparse_index_ss(overlaps(jsub,isub),blks,BLKS_NGWF,&
                 'TDRAt'//suffix//reg_struc,nze,nzb, &
                 my_nze,my_nzb,seg_nze,seg_nzb,&
                 transpose_name='TDRA'//suffix//reg_struc_t,&
                 row_par=row_par,col_par=col_par)

            ! Calculate filling factor denominator for fraction
            fill_fac = sparse_fill_fac_denom(blks,BLKS_NGWF,row_par,col_par)

            if (pub_on_root.and.show_sparsity) then
               if (mdl%nsub.eq.1) then
                  write(stdout,'(14x,a,f7.2,a)') '   TDDFT response &
                       &matrix filling:', fill_fac * nze,'%'
               else
                  write(stdout,'(8x,a,i0,a,i0,a,f7.2,a)') &
                       '   TDDFT response matrix filling (',&
                       jsub,',',isub,'):', fill_fac * nze,'%'
               end if
            end if

            ! create sparse data structures for an effective TDDFT overlap matrix
            if(pub_lr_tddft_sparse_region) then
               call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub),&
                    elems_sfc1,elems_sfc2,mdl%cell,pscode,'Region',&
                    tddft_region=.true.)
            else
               call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub),&
                    elems_sfc1,elems_sfc2,mdl%cell,pscode,'Region')
            endif

            ! Count number of corresponding nonzero matrix elements
            call sparse_count_ss(overlaps(jsub,isub),blks,BLKS_NGWF,nze,nzb,&
                 my_nze,my_nzb,seg_nze,seg_nzb,row_par=row_par,col_par=col_par)

            ! Create TDDFT kernel overlap matrix
            call sparse_index_ss(overlaps(jsub,isub),blks,BLKS_NGWF,&
                 'TDOAt'//suffix//reg_struc,nze,nzb,my_nze,&
                 my_nzb,seg_nze,seg_nzb,&
                 transpose_name='TDOA'//suffix//reg_struc_t,&
                 row_par=row_par,col_par=col_par)

            ! create sparse data structures for the transpose of TDDFT overlap matrix
            if(pub_lr_tddft_sparse_region) then
               call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub),&
                    elems_sfc1,elems_sfc2,mdl%cell,'Region',pscode,&
                    tddft_region=.true.)
            else
               call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub),&
                    elems_sfc1,elems_sfc2,mdl%cell,'Region',pscode)
            endif

            ! Count number of corresponding nonzero matrix elements
            call sparse_count_ss(overlaps(jsub,isub),BLKS_NGWF,blks,nze,nzb,&
                 my_nze,my_nzb,seg_nze,seg_nzb,row_par=row_par,col_par=col_par)

            ! Create NGWF-COND overlap matrix
            call sparse_index_ss(overlaps(jsub,isub),BLKS_NGWF,blks,&
                 'TDOA'//suffix//reg_struc,nze,nzb,my_nze,&
                 my_nzb,seg_nze,seg_nzb,&
                 transpose_name='TDOAt'//suffix//reg_struc_t,&
                 row_par=row_par,col_par=col_par)

            ! Calculate filling factor denominator for fraction
            fill_fac = sparse_fill_fac_denom(BLKS_NGWF,blks,row_par,col_par)

            ! now if this is a sparse region tddft calculation, create the structures for
            ! fully dense matrix as well: In Reality this would have to be done
            ! properly. Need the correct sparsity AND transpose structure if
            ! PcScLSvPv is NOT fully dense.
            ! here, we estimate the sparsity of PcScLSvPv by taking a kernel cutoff equal
            ! to pub_lr_tddft_kernel_cutoff+kernel_cutoff+cond_kernel_cutoff)
            call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub),&
                 elems_sfc1,elems_sfc2,mdl%cell,'Fixed','Fixed',&
                 0.5_DP*(pub_lr_tddft_kernel_cutoff+pub_kernel_cutoff+cond_kernel_cutoff))

            ! counts the number of non-zero matrix elements between cond and val
            ! ngwfs. Also works if the joint matrix is used to represent the
            ! conduction space
            call sparse_count_ss(overlaps(jsub,isub),BLKS_NGWF,blks,nze,nzb,&
                 my_nze,my_nzb,seg_nze,seg_nzb,row_par=row_par,col_par=col_par)

            ! create effective response density kernel
            call sparse_index_ss(overlaps(jsub,isub),BLKS_NGWF,blks,&
                 'TDR'//suffix//reg_struc,nze,nzb,&
                 my_nze,my_nzb,seg_nze,seg_nzb,&
                 transpose_name='TDRt'//suffix//reg_struc_t,&
                 row_par=row_par,col_par=col_par)
            ! allocate the transpose as well
            call sparse_count_ss(overlaps(jsub,isub),blks,BLKS_NGWF,nze,nzb,&
                 my_nze,my_nzb,seg_nze,seg_nzb,row_par=row_par,col_par=col_par)

            call sparse_index_ss(overlaps(jsub,isub),blks,BLKS_NGWF,&
                 'TDRt'//suffix//reg_struc,nze,nzb,&
                 my_nze,my_nzb,seg_nze,seg_nzb,&
                 transpose_name='TDR'//suffix//reg_struc_t,&
                 row_par=row_par,col_par=col_par)

            ! create sparse data structures for an effective TDDFT
            ! overlap matrix
            call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub),&
                 elems_sfc1,elems_sfc2,mdl%cell,pscode,'Region')

            ! Count number of corresponding nonzero matrix elements
            call sparse_count_ss(overlaps(jsub,isub),blks,BLKS_NGWF,nze,nzb,&
                 my_nze,my_nzb,seg_nze,seg_nzb,row_par=row_par,col_par=col_par)

            ! Create TDDFT kernel overlap matrix
            call sparse_index_ss(overlaps(jsub,isub),blks,BLKS_NGWF,&
                 'TDOt'//suffix//reg_struc,nze,nzb,my_nze,&
                 my_nzb,seg_nze,seg_nzb,&
                 transpose_name='TDO'//suffix//reg_struc_t,&
                 row_par=row_par,col_par=col_par)

            ! create sparse data structures for the transpose of TDDFT overlap matrix
            call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub),&
                 elems_sfc1,elems_sfc2,mdl%cell,'Region',pscode)

            ! Count number of corresponding nonzero matrix elements
            call sparse_count_ss(overlaps(jsub,isub),BLKS_NGWF,blks,nze,nzb,&
                 my_nze,my_nzb,seg_nze,seg_nzb,row_par=row_par,col_par=col_par)

            ! Create NGWF-COND overlap matrix
            call sparse_index_ss(overlaps(jsub,isub),BLKS_NGWF,blks,&
                 'TDO'//suffix//reg_struc,nze,nzb,my_nze,&
                 my_nzb,seg_nze,seg_nzb,&
                 transpose_name='TDOt'//suffix//reg_struc_t,&
                 row_par=row_par,col_par=col_par)

            deallocate(elems_sfc2,stat=ierr)
            call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
         end do
         deallocate(elems_sfc1,stat=ierr)
         call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
      end do

      call internal_lr_tddft_products

      do isub=1,mdl%nsub
         do jsub=1,mdl%nsub
            col_par => mdl%regions(isub)%par
            row_par => mdl%regions(jsub)%par

            fill_fac = sparse_fill_fac_denom(BLKS_NGWF,blks,row_par,col_par)
            if (pub_on_root.and.show_sparsity) then
               if (mdl%nsub.eq.1) then
                  write(stdout,'(14x,a,f7.2,a)') '           K1SvKv &
                       &matrix filling:', fill_fac*&
                       sparse_num_element(K1SvKv_lib(isub,jsub)),'%'
               else
                  write(stdout,'(8x,a,i0,a,i0,a,f7.2,a)') &
                       '           K1SvKv matrix filling (',&
                       isub,',',jsub,'):', fill_fac*&
                       sparse_num_element(K1SvKv_lib(isub,jsub)),'%'
               end if
            end if

            if (pub_on_root.and.show_sparsity) then
               if (mdl%nsub.eq.1) then
                  write(stdout,'(14x,a,f7.2,a)') '       KcScK1SvKv &
                       &matrix filling:', fill_fac*&
                       sparse_num_element(lr_tddft_proj_kernel_lib(isub,jsub)),'%'
               else
                  write(stdout,'(8x,a,i0,a,i0,a,f7.2,a)') &
                       '       KcScK1SvKv matrix filling (',&
                       isub,',',jsub,'):', fill_fac*&
                       sparse_num_element(lr_tddft_proj_kernel_lib(isub,jsub)),'%'
               end if
            end if
         end do
      end do


    end subroutine internal_create_tddft_structures


    !==========================================================================!
    ! This subroutine creates filling structures for a LR PHONONS calculation  !
    !--------------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu, December 2016                         !
    ! Adapted for embedding structures by Robert Charlton, 14/08/2018.         !
    !==========================================================================!

    subroutine internal_create_lr_phonons_structures

      do isub=1,mdl%nsub
         call internal_allocate_sfc(elems_sfc1, col_par, isub)
         do jsub=1,mdl%nsub
            call internal_allocate_sfc(elems_sfc2, row_par, jsub)
            call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

            !--------------------------------------
            ! Response density kernel
            !--------------------------------------
            ! create effective response density kernel
            call internal_create_structure(overlaps(jsub,isub),elems_sfc1,elems_sfc2, &
                 'Fixed', 'Fixed', BLKS_NGWF, blks, nze, nzb, my_nze, my_nzb, &
                 seg_nze, seg_nzb, col_par, row_par, &
                 'TDRA'//suffix//reg_struc, 'TDRAt'//suffix//reg_struc_t, &
                 0.5_DP*pub_lr_phonons_kernel_cutoff)

            ! allocate the transpose as well -- no need to recount overlaps
            call sparse_count_ss(overlaps(jsub,isub),blks,BLKS_NGWF, nze,nzb,my_nze,my_nzb, &
                  seg_nze, seg_nzb, row_par, col_par)

            call sparse_index_ss(overlaps(jsub,isub),blks,BLKS_NGWF, &
                 'TDRAt'//suffix//reg_struc,nze,nzb, &
                 my_nze,my_nzb,seg_nze,seg_nzb, &
                 transpose_name='TDRA'//suffix//reg_struc_t, &
                 row_par=row_par, col_par=col_par)

            ! Calculate filling factor denominator for fraction
            fill_fac = sparse_fill_fac_denom(blks,BLKS_NGWF,row_par,col_par)

            if (pub_on_root.and.show_sparsity) then
               if (mdl%nsub.eq.1) then
                  write(stdout,'(14x,a,f7.2,a)') '   LR PHONONS response &
                       &matrix filling:', fill_fac * nze,'%'
               else
                  write(stdout,'(8x,a,i0,a,i0,a,f7.2,a)') &
                       '   LR PHONONS response matrix filling (',&
                       jsub,',',isub,'):', fill_fac * nze,'%'
               end if
            end if

            ! Create TDDFT kernel overlap matrix
            call internal_create_structure(overlaps(jsub,isub),elems_sfc1,elems_sfc2, &
                 pscode, 'Region', blks, BLKS_NGWF, nze, nzb, my_nze, my_nzb, &
                 seg_nze, seg_nzb, col_par, row_par, &
                 'TDOAt'//suffix//reg_struc, 'TDOA'//suffix//reg_struc_t)

            ! Create NGWF-RESPNSE overlap matrix
            call internal_create_structure(overlaps(jsub,isub),elems_sfc1,elems_sfc2, &
                 'Region', pscode, BLKS_NGWF, blks, nze, nzb, my_nze, my_nzb, &
                 seg_nze, seg_nzb, col_par, row_par, &
                 'TDOA'//suffix//reg_struc, 'TDOAt'//suffix//reg_struc_t)

             ! Calculate filling factor denominator for fraction
             fill_fac = sparse_fill_fac_denom(BLKS_NGWF,blks,row_par,col_par)

             ! Need the correct sparsity AND transpose structure if PcScLSvPv is NOT
             ! fully dense. Here, we estimate the sparsity of PcScLSvPv by taking a
             ! kernel cutoff equal to pub_lr_phonons_kernel_cutoff+kernel_cutoff)
            call internal_create_structure(overlaps(jsub,isub),elems_sfc1,elems_sfc2, &
                 'Fixed', 'Fixed', BLKS_NGWF, blks, nze, nzb, my_nze, my_nzb, &
                 seg_nze, seg_nzb, col_par, row_par, &
                 'TDR'//suffix//reg_struc, 'TDRt'//suffix//reg_struc_t, &
                 0.5_DP*(pub_lr_phonons_kernel_cutoff+pub_kernel_cutoff))

            ! allocate the transpose as well -- no need to recount overlaps
            call sparse_count_ss(overlaps(jsub,isub),blks,BLKS_NGWF, nze,nzb,my_nze,my_nzb, &
                  seg_nze, seg_nzb, row_par, col_par)

            call sparse_index_ss(overlaps(jsub,isub),blks,BLKS_NGWF,'TDRt'//suffix//reg_struc,nze,nzb, &
                 my_nze,my_nzb,seg_nze,seg_nzb,transpose_name='TDR'//suffix//reg_struc_t, &
                 row_par=row_par, col_par=col_par)

            ! Create LR PHONONS kernel overlap matrix
            call internal_create_structure(overlaps(jsub,isub),elems_sfc1,elems_sfc2, &
                 pscode, 'Region', blks, BLKS_NGWF, nze, nzb, my_nze, my_nzb, &
                 seg_nze, seg_nzb, col_par, row_par, &
                 'TDOt'//suffix//reg_struc, 'TDO'//suffix//reg_struc_t)

            ! Create NGWF-RESP overlap matrix
            call internal_create_structure(overlaps(jsub,isub),elems_sfc1,elems_sfc2, &
                 'Region', pscode, BLKS_NGWF, blks, nze, nzb, my_nze, my_nzb, &
                 seg_nze, seg_nzb, col_par, row_par, &
                 'TDO'//suffix//reg_struc, 'TDOt'//suffix//reg_struc_t)

            deallocate(elems_sfc2,stat=ierr)
            call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
         end do
         deallocate(elems_sfc1,stat=ierr)
         call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
      end do

      ! rc2013: setup the LR-TDDFT matrices once we're out of the regions loop
      call internal_lr_tddft_products

      do isub=1,mdl%nsub
         do jsub=1,mdl%nsub
            col_par => mdl%regions(isub)%par
            row_par => mdl%regions(jsub)%par

            fill_fac = sparse_fill_fac_denom(BLKS_NGWF,blks,row_par,col_par)
            if (pub_on_root.and.show_sparsity) then
               if (mdl%nsub.eq.1) then
                  write(stdout,'(14x,a,f7.2,a)') '           K1SvKv &
                       &matrix filling:', fill_fac*&
                       sparse_num_element(K1SvKv_lib(jsub,isub)),'%'
               else
                  write(stdout,'(8x,a,i0,a,i0,a,f7.2,a)') &
                       '           K1SvKv matrix filling (',&
                       jsub,',',isub,'):', fill_fac*&
                       sparse_num_element(K1SvKv_lib(jsub,isub)),'%'
               end if
            end if

            if (pub_on_root.and.show_sparsity) then
               if (mdl%nsub.eq.1) then
                  write(stdout,'(14x,a,f7.2,a)') '       KcScK1SvKv &
                       &matrix filling:', fill_fac*&
                       sparse_num_element(lr_tddft_proj_kernel_lib(jsub,isub)),'%'
               else
                  write(stdout,'(8x,a,i0,a,i0,a,f7.2,a)') &
                       '       KcScK1SvKv matrix filling (',&
                       jsub,',',isub,'):', fill_fac*&
                       sparse_num_element(lr_tddft_proj_kernel_lib(jsub,isub)),'%'
               end if
            end if
         end do
      end do


     end subroutine internal_create_lr_phonons_structures

     !=========================================================================!
     ! This subroutine creates the filling structure of the projected TDDFT    !
     ! density kernel with the correct sparsity pattern                        !
     !-------------------------------------------------------------------------!
     ! Written by Tim Zuehlsdorff, February 2014                               !
     ! Copied and modified for embedding by Joseph Prentice, August 2018       !
     !=========================================================================!

     subroutine internal_lr_tddft_products

       implicit none

       ! local variables
       type(SPAM3_EMBED) :: k1,kc, sc, kv, sv, kcsc,svkv
       type(SPAM3_EMBED) :: k1svkv, kcsck1svkv

       k1%structure = 'TDRA'//suffix ! TDDFT auxillary response density kernel
       kc%structure = 'K'//suffix   ! Conduction density kernel
       sc%structure = 'S'//suffix   ! conduction NGWF overlap
       kv%structure = 'K'           ! Valence density kernel
       sv%structure = 'S'           ! Valence NGWF overlap

       ! jcap: Create block matrices - we need the SPAM3_EMBED info
       ! for matrix products, so we can't destroy them until we're
       ! finished
       call sparse_embed_create(k1, mrows=mdl%nsub, ncols=mdl%nsub)
       call sparse_embed_create(kc, mrows=mdl%nsub, ncols=mdl%nsub)
       call sparse_embed_create(sc, mrows=mdl%nsub, ncols=mdl%nsub)
       call sparse_embed_create(kv, mrows=mdl%nsub, ncols=mdl%nsub)
       call sparse_embed_create(sv, mrows=mdl%nsub, ncols=mdl%nsub)

       call sparse_embed_create(kcsc,kc,sc)
       call sparse_embed_create(svkv,sv,kv)
       call sparse_embed_create(k1svkv,k1,svkv,rlib=K1SvKv_lib)
       call sparse_embed_create(kcsck1svkv,kcsc,k1svkv,rlib=lr_tddft_proj_kernel_lib)

       ! jcap: now destroy everything
       call sparse_embed_destroy(k1)
       call sparse_embed_destroy(kc)
       call sparse_embed_destroy(sc)
       call sparse_embed_destroy(kv)
       call sparse_embed_destroy(sv)
       call sparse_embed_destroy(kcsc)
       call sparse_embed_destroy(svkv)
       call sparse_embed_destroy(k1svkv)

     end subroutine internal_lr_tddft_products


    !==========================================================================!
    ! This subroutine creates the filling structures of the projected          !
    ! conduction Hamiltonian and associated gradient structure                 !
    !--------------------------------------------------------------------------!
    ! Written by Laura Ratcliff, June 2010.                                    !
    !==========================================================================!

    subroutine internal_cond_ham_products

      implicit none

      ! Local variables
      type(SPAM3_EMBED) :: k, h, u, t, s, m
      type(SPAM3_EMBED) :: kt, kh, ks, uk
      type(SPAM3_EMBED) :: ukh, uks
      type(SPAM3_EMBED) :: kskt, khkt
      type(SPAM3_EMBED) :: ksktm, ukhkt, ukskt, khktm

      k%structure = 'K'         ! valence kernel
      h%structure = 'H'         ! valence hamiltonian
      u%structure = 'U'//suffix ! cross overlap <\chi_\alpha|\phi_\beta>
      t%structure = 'T'//suffix ! cross overlap <\phi_\alpha|\chi_\beta>
      s%structure = 'S'         ! valence overlap
      m%structure = 'K'//suffix ! cond kernel

      ! rc2013: we need the SPAM3_EMBED info for matrix products, so
      ! we can't destroy them until we're finished

      ! lr408: Products for projected conduction Hamiltonian
      call sparse_embed_create(m, mrows=mdl%nsub, ncols=mdl%nsub)
      call sparse_embed_create(s, mrows=mdl%nsub, ncols=mdl%nsub)
      call sparse_embed_create(k, mrows=mdl%nsub, ncols=mdl%nsub)
      call sparse_embed_create(h, mrows=mdl%nsub, ncols=mdl%nsub)
      call sparse_embed_create(u, mrows=mdl%nsub, ncols=mdl%nsub)
      call sparse_embed_create(t, mrows=mdl%nsub, ncols=mdl%nsub)
      call sparse_embed_create(kt,k,t)
      call sparse_embed_create(uk,u,k)
      call sparse_embed_create(ukh,uk,h)
      call sparse_embed_create(uks,uk,s)
      call sparse_embed_create(kh,k,h)
      call sparse_embed_create(khkt,kh,kt)
      call sparse_embed_create(ukhkt,u,khkt,rlib=cond_proj_ham_lib)
      call sparse_embed_create(ukskt,uks,kt)

      ! lr408:  products for projected conduction Hamiltonian gradient
      call sparse_embed_create(khktm,khkt,m,rlib=cond_grad_lib1)
      call sparse_embed_create(ks,k,s)
      call sparse_embed_create(kskt,ks,kt)
      call sparse_embed_create(ksktm,kskt,m,rlib=cond_grad_lib2)

      ! rc2013: now destroy everything
      call sparse_embed_destroy(m)
      call sparse_embed_destroy(s)
      call sparse_embed_destroy(k)
      call sparse_embed_destroy(h)
      call sparse_embed_destroy(u)
      call sparse_embed_destroy(t)
      call sparse_embed_destroy(kt)
      call sparse_embed_destroy(uk)
      call sparse_embed_destroy(ukh)
      call sparse_embed_destroy(uks)
      call sparse_embed_destroy(kh)
      call sparse_embed_destroy(khkt)
      call sparse_embed_destroy(ukhkt)
      call sparse_embed_destroy(ukskt)
      call sparse_embed_destroy(khktm)
      call sparse_embed_destroy(ks)
      call sparse_embed_destroy(kskt)
      call sparse_embed_destroy(ksktm)

    end subroutine internal_cond_ham_products


    !==========================================================================!
    ! This subroutine creates the structures associated with products of the   !
    ! projected conduction Hamiltonian and other matrices.                     !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in November 2011.                               !
    !==========================================================================!

    subroutine internal_cond_extra_products

      implicit none

      ! Local variables
      type(SPAM3_EMBED) :: k, l, t, s, kv, m
      type(SPAM3_EMBED) :: kl, ks, lk, tk, ms
      type(SPAM3_EMBED) :: skl, lks, lkl, klk, kvtk

      k%structure = 'K'//suffix ! cond kernel
      l%structure = 'L'//suffix ! projected cond hamiltonian
      t%structure = 'T'//suffix ! cross overlap <\phi_\alpha|\chi_\beta>
      s%structure = 'S'//suffix ! cond overlap
      m%structure = 'M'//suffix ! cond gradient term
      kv%structure = 'K'        ! valence kernel

      ! rc2013: we need the SPAM3_EMBED info for matrix products, so
      ! we can't destroy them until we're finished

      ! ndmh: Products involving projected conduction Hamiltonian
      call sparse_embed_create(k, mrows=mdl%nsub, ncols=mdl%nsub)
      call sparse_embed_create(l, mrows=mdl%nsub, ncols=mdl%nsub)
      call sparse_embed_create(s, mrows=mdl%nsub, ncols=mdl%nsub)
      call sparse_embed_create(kl,k,l)
      call sparse_embed_create(lk,l,k)
      call sparse_embed_create(ks,k,s)
      call sparse_embed_create(skl,s,kl)
      call sparse_embed_create(lks,l,ks)
      call sparse_embed_create(lkl,l,kl)
      call sparse_embed_create(klk,k,lk)

      ! ndmh: Products involving gradient term m
      call sparse_embed_create(m, mrows=mdl%nsub, ncols=mdl%nsub)
      call sparse_embed_create(ms,m,s)

      ! ndmh: other cond matrices
      call sparse_embed_create(t, mrows=mdl%nsub, ncols=mdl%nsub)
      call sparse_embed_create(kv, mrows=mdl%nsub, ncols=mdl%nsub)
      call sparse_embed_create(tk,t,k)
      call sparse_embed_create(kvtk,kv,tk)

      ! rc2013: now destroy everything
      call sparse_embed_destroy(k)
      call sparse_embed_destroy(l)
      call sparse_embed_destroy(s)
      call sparse_embed_destroy(kl)
      call sparse_embed_destroy(lk)
      call sparse_embed_destroy(ks)
      call sparse_embed_destroy(skl)
      call sparse_embed_destroy(lks)
      call sparse_embed_destroy(lkl)
      call sparse_embed_destroy(klk)
      call sparse_embed_destroy(m)
      call sparse_embed_destroy(ms)
      call sparse_embed_destroy(t)
      call sparse_embed_destroy(kv)
      call sparse_embed_destroy(tk)
      call sparse_embed_destroy(kvtk)

    end subroutine internal_cond_extra_products

    !==========================================================================!
    ! This subroutine creates the structures of the ks, ... ,ksksk matrices    !
    ! so that their filling factors can be displayed at init time              !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine, July 2008.                                     !
    ! Modified by Quintin Hill, 02/02/2009 to include LSLSL for HF exchange.   !
    ! Updated for SPAM3 by Nicholas Hine, May 2009.                            !
    ! Modified for extra memory contiguousness by Nicholas Hine in Jan 2010.   !
    ! Modified to allow reuse for conduction matrices by Laura Ratcliff        !
    ! in Oct 2010.                                                             !
    ! Modified for embedding by Robert Charlton, 14/07/17.                     !
    !==========================================================================!

    subroutine internal_create_ksksk

      implicit none

      ! Local variables
      type(SPAM3_EMBED) :: k,s,ks,ksk,sk,sks,ksks,ksksk,kk

      k%structure = 'K'//suffix
      s%structure = 'S'//suffix

      ! Create the matrices and immediately destroy them to free up the data
      ! part... the structures will be kept in each case (this allows contiguous
      ! blocks of memory to remain after the creation of the structures)
      call sparse_embed_create(k,mrows=mdl%nsub,ncols=mdl%nsub)
      call sparse_embed_create(s,mrows=mdl%nsub,ncols=mdl%nsub)
      call sparse_embed_create(ks,k,s,rlib=kslib)
      call sparse_embed_create(ksk,ks,k)
      call sparse_embed_create(ksks,ksk,s,rlib=kskslib)
      call sparse_embed_create(sk,s,k)
      call sparse_embed_create(sks,sk,s)
      call sparse_embed_create(kk,k,k)

      ! With Hartree-Fock exchange, sparsity pattern of Hamiltonian is KSKSK
      if (pub_use_hfx.or.pub_use_activehfx) then
         call sparse_embed_create(ksksk,ksks,k,rlib=ksksklib)
         call sparse_embed_destroy(ksksk)
       end if
      ! rc2013: don't destroy matrices until we're finished
      ! so that blocking info is correct
      call sparse_embed_destroy(k)
      call sparse_embed_destroy(s)
      call sparse_embed_destroy(ks)
      call sparse_embed_destroy(ksk)
      call sparse_embed_destroy(ksks)
      call sparse_embed_destroy(sk)
      call sparse_embed_destroy(sks)
      call sparse_embed_destroy(kk)

    end subroutine internal_create_ksksk


    !==========================================================================!
    ! This subroutine creates the structures of all required matrices relating !
    ! to the Hamiltonian, to prevent creation of new structures at run time.   !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in Jan 2010.                                    !
    ! Modified to allow reuse for conduction matrices by Laura Ratcliff        !
    ! in Oct 2010.                                                             !
    ! Modified to allow embedding projector matrices by Robert Charlton,       !
    ! 15/12/17.                                                                !
    !==========================================================================!

    subroutine internal_create_ham_strucs

      implicit none

      ! Local variables
      type(SPAM3_EMBED) :: k,s,h,kh,hk,hkh,khk,skh,hks,cch,cchcch,cchk,f

      k%structure = 'K'//suffix
      h%structure = 'H'//suffix
      s%structure = 'S'//suffix

      ! Create the matrices and immediately destroy them to free up the data
      ! part... the structures will be kept in each case (this allows contiguous
      ! blocks of memory to remain after the creation of the structures)
      call sparse_embed_create(k,mrows=mdl%nsub,ncols=mdl%nsub)
      call sparse_embed_create(s,mrows=mdl%nsub,ncols=mdl%nsub)
      call sparse_embed_create(h,mrows=mdl%nsub,ncols=mdl%nsub)
      call sparse_embed_create(kh,k,h)
      call sparse_embed_create(hk,h,k)
      call sparse_embed_create(hkh,h,kh)
      call sparse_embed_create(khk,kh,k)
      call sparse_embed_create(skh,s,kh)
      call sparse_embed_create(hks,hk,s)

      if(abs(pub_contracoham_radmult-1.0_dp)>epsilon(1.0_dp)) then
         cch%structure = 'ccH'//suffix
         call sparse_embed_create(cch,mrows=mdl%nsub,ncols=mdl%nsub)
         call sparse_embed_create(cchcch,cch,cch)
         call sparse_embed_create(cchk,cch,k)
         call sparse_embed_destroy(cch)
         call sparse_embed_destroy(cchcch)
         call sparse_embed_destroy(cchk)
      end if

      ! rc2013: don't destroy matrices until we're finished
      ! so that blocking info is correct
      call sparse_embed_destroy(k)
      call sparse_embed_destroy(s)
      call sparse_embed_destroy(h)
      call sparse_embed_destroy(kh)
      call sparse_embed_destroy(hk)
      call sparse_embed_destroy(hkh)
      call sparse_embed_destroy(khk)
      call sparse_embed_destroy(skh)
      call sparse_embed_destroy(hks)

    end subroutine internal_create_ham_strucs


    !==========================================================================!
    ! This subroutine creates the structures of all required matrices relating !
    ! to nonlocal psps, to prevent creation of new structures at run time.     !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in Jan 2010.                                    !
    ! Modified to allow reuse for conduction matrices by Laura Ratcliff        !
    ! in Oct 2010.                                                             !
    !==========================================================================!

    subroutine internal_create_nl_strucs

      implicit none

      ! Local variables
      type(SPAM3_EMBED) :: r,q,k,s,rk,rks,kq

      r%structure = 'R'//suffix
      q%structure = 'Q'//suffix
      k%structure = 'K'//suffix
      s%structure = 'S'//suffix

      ! Create the matrices and immediately destroy them to free up the data
      ! part... the structures will be kept in each case (this allows contiguous
      ! blocks of memory to remain after the creation of the structures)
      call sparse_embed_create(r,mrows=mdl%nsub,ncols=mdl%nsub)
      call sparse_embed_create(q,mrows=mdl%nsub,ncols=mdl%nsub)
      call sparse_embed_create(k,mrows=mdl%nsub,ncols=mdl%nsub)
      call sparse_embed_create(s,mrows=mdl%nsub,ncols=mdl%nsub)
      call sparse_embed_create(rk,r,k)
      call sparse_embed_create(rks,rk,s)
      call sparse_embed_create(kq,k,q)

      ! rc2013: don't destroy matrices until we're finished
      ! so that blocking info is correct
      call sparse_embed_destroy(r)
      call sparse_embed_destroy(q)
      call sparse_embed_destroy(k)
      call sparse_embed_destroy(s)
      call sparse_embed_destroy(rk)
      call sparse_embed_destroy(rks)
      call sparse_embed_destroy(kq)

    end subroutine internal_create_nl_strucs

    !==========================================================================!
    ! This subroutine creates the structures of all required matrices relating !
    ! to Hubbard projectors, to prevent creation of new structures at run time.!
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in Jan 2010.                                    !
    !==========================================================================!

    subroutine internal_create_hubbard_strucs

      implicit none

      ! Local variables
      type(SPAM3_EMBED) :: v,w,k,s,wk,wks,kv

      v%structure = 'V'//suffix
      w%structure = 'W'//suffix
      k%structure = 'K'//suffix
      s%structure = 'S'//suffix

      ! Create the matrices and immediately destroy them to free up the data
      ! part... the structures will be kept in each case (this allows contiguous
      ! blocks of memory to remain after the creation of the structures)
      ! rc2013: specify the number of rows/cols in the matrix arrays
      call sparse_embed_create(v,mrows=mdl%nsub,ncols=mdl%nsub)
      call sparse_embed_create(w,mrows=mdl%nsub,ncols=mdl%nsub)
      call sparse_embed_create(k,mrows=mdl%nsub,ncols=mdl%nsub)
      call sparse_embed_create(s,mrows=mdl%nsub,ncols=mdl%nsub)
      call sparse_embed_create(wk,w,k)
      call sparse_embed_create(wks,wk,s)
      call sparse_embed_create(kv,k,v)
      ! rc2013: don't destroy matrices until we're finished
      ! so that blocking info is correct
      call sparse_embed_destroy(v)
      call sparse_embed_destroy(w)
      call sparse_embed_destroy(k)
      call sparse_embed_destroy(s)
      call sparse_embed_destroy(wk)
      call sparse_embed_destroy(wks)
      call sparse_embed_destroy(kv)

    end subroutine internal_create_hubbard_strucs

    !==========================================================================!
    ! This subroutine creates the structure of a matrix product from two given !
    ! structure codes. Matrix product indexing is normally done automatically  !
    ! upon creation of a new matrix with an unknown code, but for the          !
    ! creation of the Hamiltonian matrix we need the patterns of some extra    !
    ! matrices such as QR and VW, so this routine can be used to force the     !
    ! necessary call to sparse_embed_create.                                   !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine, May 2009 as internal_create_qr.                !
    ! Generalised to any matrix product by Nicholas Hine in September 2009.    !
    ! Modified for extra memory contiguousness by Nicholas Hine in Jan 2010.   !
    ! Rewritten with embedding structures by Robert Charlton, 25/09/2017.      !
    !==========================================================================!

    subroutine internal_create_product(newlib,astruc,bstruc)

      implicit none

      ! Arguments
      integer, intent(out) :: newlib(mdl%nsub, mdl%nsub)
      character(len=*), intent(in) :: astruc,bstruc

      ! Local variables
      type(SPAM3_EMBED) :: amat,bmat,cmat

      amat%structure = astruc
      bmat%structure = bstruc

      ! Create temporary matrices and immediately destroy them to free up the
      ! memory allocate to the data before creating new structures
      ! rc2013: for embedding matrices we need all the system info
      call sparse_embed_create(amat, mrows=mdl%nsub, ncols=mdl%nsub)
      call sparse_embed_create(bmat, mrows=mdl%nsub, ncols=mdl%nsub)

      ! Create the product and store its library entry
      call sparse_embed_create(cmat,amat,bmat,rlib=newlib)

      ! Destroy temporary matrices
      call sparse_embed_destroy(amat)
      call sparse_embed_destroy(bmat)
      call sparse_embed_destroy(cmat)

      ! Synchronise
      call comms_barrier

    end subroutine internal_create_product

    !==========================================================================!
    ! This subroutine creates the structure of the vmatrix/omatrix used in the !
    ! sw_expansion module.  The vmatrix is the SW-SW two electron integral     !
    ! metric matrix under the electrostatic metric. The omatrix is the SW-SW   !
    ! two electron integral metric matrix under the usual (overlap) metric.    !
    ! This matrix has the same atom block sparsity pattern as the purified     !
    ! density kernel (but with different sized blocks).                        !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill, October 2009.                                   !
    ! Modified by Nicholas Hine and Jacek Dziedzic in April 2012 for SW metric !
    ! Modified by Jacek Dziedzic in April 2014 to work with SW expansion even  !
    ! in the absence of HFx.                                                   !
    ! Modified for embedding by Joseph Prentice, June 2018                     !
    !==========================================================================!

    subroutine internal_create_metric_matrix

      use rundat, only: pub_use_swx, pub_use_activeswx
      use utils, only: utils_assert, utils_dealloc_check

      implicit none

      call utils_assert(pub_use_swx.or.pub_use_activeswx, &
           'Inconsistency in sparse_initialise::internal_create_metric_matrix.')

      ! jcap: loop over regions
      do isub=1,mdl%nsub
         call internal_allocate_sfc(elems_sfc1, col_par, isub)
         do jsub=1,mdl%nsub
            call internal_allocate_sfc(elems_sfc2, row_par, jsub)
            call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

            ! rc2013: with embedding HFx can only be used in 1 region, so
            ! only initialise that matrix structure
            if(isub==pub_active_region .and. jsub==pub_active_region) then
               ! qoh: Create SW-SW overlap matrix (O matrix)
               call internal_create_structure(overlaps(jsub,isub), elems_sfc1, &
                    elems_sfc2, pscode, pscode, BLKS_SW, BLKS_SW, nze, nzb, &
                    my_nze, my_nzb, seg_nze, seg_nzb, col_par, row_par, &
                    'Ss'//suffix//reg_struc, &
                    'Ss'//suffix//reg_struc_t)

               ! qoh: Create SW-SW overlap matrix (V matrix)
               call sparse_index_ss(overlaps(jsub,isub),BLKS_SW,BLKS_SW, &
                    'Vs'//suffix//reg_struc,nze,nzb,my_nze,my_nzb, &
                    seg_nze,seg_nzb,transpose_name='Vs'//suffix//reg_struc_t, &
                    row_par=row_par,col_par=col_par)
            end if

            deallocate(elems_sfc2,stat=ierr)
            call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
         end do
         deallocate(elems_sfc1,stat=ierr)
         call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
      end do

    end subroutine internal_create_metric_matrix

    !==========================================================================!
    ! This subroutine creates the structure of the exchange (X) matrix used in !
    ! Hartree-Fock exchange.                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill, October 2009.                                   !
    ! Modified by Nicholas Hine and Jacek Dziedzic in April 2012 for SW metric !
    ! Separated from internal_create_metric_matrix by Jacek Dziedzic in April  !
    ! 2014.                                                                    !
    ! Modified for embedding by Joseph Prentice, June 2018                     !
    !==========================================================================!

    subroutine internal_create_exchange_matrix

      use rundat, only: pub_use_hfx, pub_use_activehfx, pub_hfx_cutoff
      use utils, only: utils_assert, utils_dealloc_check

      implicit none

      call utils_assert(pub_use_hfx.or.pub_use_activehfx, &
           'Inconsistency in sparse_initialise::internal_create_exchange_matrix.')

      ! jcap: loop over regions
      do isub=1,mdl%nsub
         call internal_allocate_sfc(elems_sfc1, col_par, isub)
         do jsub=1,mdl%nsub
            call internal_allocate_sfc(elems_sfc2, row_par, jsub)
            call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

            ! Create X matrix
            call internal_create_structure(overlaps(jsub,isub), elems_sfc1, &
                 elems_sfc2, 'Fixed', 'Fixed', blks, blks, nze, nzb, &
                 my_nze, my_nzb, seg_nze, seg_nzb, col_par, row_par, &
                 'HFx'//suffix//reg_struc, &
                 'HFx'//suffix//reg_struc_t, radius=0.5_DP * pub_hfx_cutoff)

            deallocate(elems_sfc2,stat=ierr)
            call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
         end do
         deallocate(elems_sfc1,stat=ierr)
         call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
      end do

    end subroutine internal_create_exchange_matrix

    !==========================================================================!
    ! This subroutine creates the structures for the block-diagonal matrices:  !
    ! projector diagonal, Hubbard projectors, NGWF diagonal matrix.            !
    !--------------------------------------------------------------------------!
    ! Extracted to separate subroutine by Robert Charlton, October 2017.       !
    !==========================================================================!

    subroutine internal_create_block_diag

      implicit none

      ! Arguments

      ! Local variables
      real(kind=DP) :: range

      do isub=1,mdl%nsub
         call internal_allocate_sfc(elems_sfc1, col_par, isub)
         ! rc2013: if a diagonal block of the overlaps_list is empty
         ! e.g. no projectors reallocate it now
         if(size(overlaps(isub,isub)%overlap_list,dim=1) == 0) then
            deallocate(overlaps(isub,isub)%overlap_list,stat=ierr)
            call utils_dealloc_check('sparse_init_rep', &
                 'overlaps%overlap_list',ierr)
            allocate(overlaps(isub,isub)%overlap_list(1,col_par%nat),stat=ierr)
            call utils_alloc_check('sparse_init_rep', &
                 'overlaps%overlap_list',ierr)
         end if
         do jsub=1,mdl%nsub
            call internal_allocate_sfc(elems_sfc2, row_par, jsub)
            call internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

            first_atom_on_proc = col_par%first_atom_on_proc(pub_my_proc_id)
            last_atom_on_proc = first_atom_on_proc + &
                 col_par%num_atoms_on_proc(pub_my_proc_id) - 1
            if (suffix=='') then

               ! Projector diagonal matrix
               if (any_proj) then
                  do iblk=first_atom_on_proc,last_atom_on_proc
                     ! rc2013: all off-diagonal blocks will be zero, so omit
                     if (sparse_num_elems_on_atom(iblk,BLKS_PROJ,col_par)>0 &
                          .and. isub==jsub) then
                        overlaps(jsub,isub)%num_overlaps(iblk) = 1
                        overlaps(jsub,isub)%overlap_list(1,iblk) = iblk
                     else
                        overlaps(jsub,isub)%num_overlaps(iblk) = 0
                     end if
                  end do
                  call sparse_count_ss(overlaps(jsub,isub),BLKS_PROJ,BLKS_PROJ,nze,nzb, &
                       my_nze,my_nzb,seg_nze,seg_nzb, row_par, col_par)
                  call sparse_index_ss(overlaps(jsub,isub),BLKS_PROJ,BLKS_PROJ, &
                       'E'//suffix//reg_struc,nze,nzb,my_nze,my_nzb, &
                       seg_nze,seg_nzb,transpose_name='E'//suffix//reg_struc_t, &
                       row_par=row_par, col_par=col_par)
               end if

               ! Hubbard Projector diagonal matrix
               if (pub_hubbard) then
                  do iblk=first_atom_on_proc, last_atom_on_proc
                     ! rc2013: all off-diagonal blocks will be zero, so omit
                     if (sparse_num_elems_on_atom(iblk,BLKS_HUB_PROJ,col_par)>0 .and. isub==jsub) then
                        overlaps(jsub,isub)%num_overlaps(iblk) = 1
                        overlaps(jsub,isub)%overlap_list(1,iblk) = iblk
                     else
                        overlaps(jsub,isub)%num_overlaps(iblk) = 0
                     end if
                  end do
                  call sparse_count_ss(overlaps(jsub,isub),BLKS_HUB_PROJ,BLKS_HUB_PROJ,&
                       nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb,row_par,col_par)
                  call sparse_index_ss(overlaps(jsub,isub),BLKS_HUB_PROJ,BLKS_HUB_PROJ,&
                       'G'//reg_struc,nze,nzb,my_nze,my_nzb, &
                       seg_nze,seg_nzb,transpose_name='G'//reg_struc_t, &
                       row_par=row_par, col_par=col_par)

                  ! CW: FULL PROJ
                  ! ebl: This code is associated with a (presently) commented line
                  !      in parallel_strategy mod
                  !      It requires further thought before uncommenting
                  ! my_first_blk = par%first_atom_on_proc(pub_my_proc_id)
                  ! my_last_blk = par%first_atom_on_proc(pub_my_proc_id + 1) - 1
                  ! do iblk=my_first_blk,my_last_blk
                  !    if (sparse_num_elems_on_atom(iblk,BLKS_HUB_PROJ,par)>0) then
                  !       j=1
                  !       overlaps%overlap_list(1,iblk) = iblk

                  !       do i=my_first_blk,my_last_blk
                  !         if (i/=iblk.and.sparse_num_elems_on_atom(i,BLKS_HUB_PROJ,par)>0.and. &
                  !                       & sparse_num_elems_on_atom(i,BLKS_HUB_PROJ,par)== &
                  !                       & sparse_num_elems_on_atom(iblk,BLKS_HUB_PROJ,par)) then
                  !           j=j+1
                  !           overlaps%overlap_list(j,iblk) = i
                  !         endif
                  !       enddo
                  !       overlaps%num_overlaps(iblk) = j
                  !       write(*,*) 'hubbard atom [x] has [y] buddies : ', iblk,j
                  !       write(*,*) 'his pals are                     : ', &
                  !            overlaps%overlap_list(1:j,iblk)
                  !       write(*,*) 'NUMS elements on atoms           : ', &
                  !            sparse_num_elems_on_atom(iblk,BLKS_HUB_PROJ,par)
                  !    else
                  !       overlaps%num_overlaps(iblk) = 0
                  !    end if
                  ! end do
                  ! call sparse_count_ss(overlaps, BLKS_HUB_PROJ,BLKS_HUB_PROJ, &
                  !      nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb,row_par,col_par)
                  ! call sparse_index_ss(overlaps, BLKS_HUB_PROJ,BLKS_HUB_PROJ, &
                  !      'FULL',nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb, row_par=par, col_par=par)
                  ! if (pub_on_root) then
                  !    funit = utils_unit()
                  !    open(unit=funit,file='positions_cor_atoms_in_onetep')
                  !       do iblk=my_first_blk,my_last_blk
                  !          if (sparse_num_elems_on_atom(iblk,BLKS_HUB_PROJ,par)>0) then
                  !             write(funit,'(a,3f10.4)') &
                  !                  trim(adjustl(par%elements_on_proc(iblk-my_first_blk+1)%symbol)), &
                  !                  par%elements_on_proc(iblk-my_first_blk+1)%centre%x*0.529177249, &
                  !                  par%elements_on_proc(iblk-my_first_blk+1)%centre%y*0.529177249, &
                  !                  par%elements_on_proc(iblk-my_first_blk+1)%centre%z*0.529177249
                  !          endif
                  !       enddo
                  !    close(funit)
                  ! endif
                  ! END CW

               end if


               if (pub_eels_calculate) then ! lr408
                  do iblk=first_atom_on_proc,last_atom_on_proc
                     if (sparse_num_elems_on_atom(iblk,BLKS_PROJ,col_par)>0 .and. &
                          sparse_num_elems_on_atom(iblk,BLKS_CORE,col_par)>0 &
                          .and. isub==jsub) then
                        overlaps(jsub,isub)%num_overlaps(iblk) = 1
                        overlaps(jsub,isub)%overlap_list(1,iblk) = iblk
                     else
                        overlaps(jsub,isub)%num_overlaps(iblk) = 0
                     end if
                  end do
                  call sparse_count_ss(overlaps(jsub,isub),BLKS_CORE,BLKS_PROJ,&
                       nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb, row_par, col_par)
                  call sparse_index_ss(overlaps(jsub,isub),BLKS_CORE,BLKS_PROJ,&
                       'J',nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb, &
                       transpose_name='P'//reg_struc_t, &
                       row_par=row_par,col_par=col_par)
                  call sparse_count_ss(overlaps(jsub,isub),BLKS_PROJ,BLKS_CORE,&
                       nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb,row_par,col_par)
                  call sparse_index_ss(overlaps(jsub,isub),BLKS_PROJ,BLKS_CORE,&
                       'P',nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb, &
                       transpose_name='J'//reg_struc_t, &
                       row_par=row_par, col_par=col_par)

               end if

            end if

            ! NGWF diagonal matrix for this rep
            ! rc2013: build the structures but only diagonal blocks will
            ! have non-zero elements
            fill_fac = sparse_fill_fac_denom(blks,blks,row_par,col_par)
            do iblk=first_atom_on_proc,last_atom_on_proc
               if (isub==jsub) then
                  overlaps(jsub,isub)%num_overlaps(iblk) = 1
                  overlaps(jsub,isub)%overlap_list(1,iblk) = iblk
               endif
            end do
            call sparse_count_ss(overlaps(jsub,isub),blks,blks,nze,nzb,my_nze,my_nzb, &
                 seg_nze,seg_nzb,row_par,col_par)
            call sparse_index_ss(overlaps(jsub,isub),blks,blks, &
                 'D'//reg_struc,nze,nzb,my_nze,my_nzb, &
                 seg_nze,seg_nzb,transpose_name='D'//reg_struc_t, &
                 row_par=row_par,col_par=col_par)
            if (pub_on_root.and.show_sparsity) then
               if (mdl%nsub.eq.1) then
                  write(stdout,'(18x,a,f7.2,a)') 'NGWF diagonal &
                       &matrix filling:', fill_fac * nze,'%'
               else
                  write(stdout,'(12x,a,i0,a,i0,a,f7.2,a)') &
                       'NGWF diagonal matrix filling (',&
                       jsub,',',isub,'):', fill_fac * nze,'%'
               end if
            end if
            deallocate(elems_sfc2,stat=ierr)
            call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
         end do
         deallocate(elems_sfc1,stat=ierr)
         call utils_dealloc_check('sparse_init_rep','elems_sfc',ierr)
      end do

    end subroutine internal_create_block_diag

    !==========================================================================!
    ! This subroutine creates the structures of various matrices associated    !
    ! with the calculation of conduction NGWFs                                 !
    !--------------------------------------------------------------------------!
    ! Created by Robert Charlton, 11/01/2018.                                  !
    !==========================================================================!

    subroutine internal_create_embed_structures

      type(ELEMENT), allocatable :: elems_sfc1(:)
      type(ELEMENT), allocatable :: elems_sfc2(:)

      !---------------------------------------------------------------------------
      ! Projected embedding Hamiltonian
      !---------------------------------------------------------------------------

      do isub=1,mdl%nsub
         do jsub=1,mdl%nsub
            write(isub_str,'(i1)') isub
            write(jsub_str,'(i1)') jsub
            call sparse_count_union(proj_embed_ham_lib(jsub,isub),hlib(jsub,isub),nze,nzb,my_nze,my_nzb, &
                 seg_nze,seg_nzb)

            call sparse_index_union(proj_embed_ham_lib(jsub,isub),hlib(jsub,isub),'F'//suffix,nze,nzb, &
                 my_nze,my_nzb,seg_nze,seg_nzb)

            fill_fac = sparse_fill_fac_denom(BLKS_NGWF,BLKS_NGWF,row_par,col_par)

            ! Find the number of nonzero elements in the patterns
            ! of commonly-used matrices and print results as part of filling
            ! factor summary
            if (pub_on_root.and.show_sparsity) then
               if (mdl%nsub.eq.1) then
                  write(stdout,'(8x,a,f7.2,a)') 'Proj. embedding Hamiltonian &
                       &matrix filling:', fill_fac * nze,'%'
               else
                  write(stdout,'(2x,a,i0,a,i0,a,f7.2,a)') &
                       'Proj. embedding Hamiltonian matrix filling (',&
                       jsub,',',isub,'):', fill_fac * nze,'%'
               end if
            end if

            ! Gradient term !k...k_cond sparsity
            call sparse_count_union(embed_grad_lib1,embed_grad_lib2,nze,nzb,my_nze, &
                 my_nzb,seg_nze,seg_nzb)

            call sparse_index_union(embed_grad_lib1,embed_grad_lib2,'M'//suffix,nze, &
                 nzb,my_nze,my_nzb,seg_nze,seg_nzb)

         end do
      end do

    end subroutine internal_create_embed_structures

    !==========================================================================!
    ! This subroutine allocates the elements list in space-filling curve order,!
    ! assigns the related parallel strategy and writes the region number to a  !
    ! string for later use in sparse_embed_mod.                                !
    !--------------------------------------------------------------------------!
    ! Created by Robert Charlton, 25/04/2018.                                  !
    !==========================================================================!

    subroutine internal_allocate_sfc(elems_sfc, this_par, isub)

      implicit none

      ! Arguments
      type(ELEMENT), allocatable, intent(inout) :: elems_sfc(:)
      type(PARAL_INFO), pointer, intent(inout)  :: this_par
      integer, intent(in)                       :: isub

      ! rc2013: allocate the array
      allocate(elems_sfc(mdl%regions(isub)%par%nat),stat=ierr)
      call utils_alloc_check('sparse_init_rep', 'elems_sfc', ierr)

      this_par => mdl%regions(isub)%par
      do iat=1,this_par%nat
         elems_sfc(iat) = system_info(isub)%elems_sfc(iat)
      end do

    end subroutine internal_allocate_sfc

    !==========================================================================!
    ! This subroutine sets the embedding sub-block structure for a given SPAM3 !
    ! and its transpose. If there's only 1 region then the structure is blank. !
    !--------------------------------------------------------------------------!
    ! Created by Robert Charlton, 27/09/2018.                                  !
    !==========================================================================!

    subroutine internal_write_embed_struc(reg_struc, reg_struc_t, jsub, isub)

      implicit none

      ! Arguments
      character(len=*), intent(out) :: reg_struc
      character(len=*), intent(out) :: reg_struc_t
      integer, intent(in)           :: jsub
      integer, intent(in)           :: isub

      if(mdl%nsub .gt. 1) then
         write(isub_str,'(i1)') isub
         write(jsub_str,'(i1)') jsub
         reg_struc   = trim(jsub_str//isub_str)
         reg_struc_t = trim(isub_str//jsub_str)
      else
         reg_struc   = ''
         reg_struc_t = ''
      end if

    end subroutine internal_write_embed_struc

  end subroutine sparse_init_rep

end module sparse_initialise
