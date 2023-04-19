! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!        LR-TDDFT utils module                                   !
!                                                                !
! The module contains a number of utility routines that are      !
! used in several other modules related to TDDFT, namely         !
! lr_tddft_mod, lr_tddft_RPA_mod and lr_tddft_forces_mod.        !
!----------------------------------------------------------------!
! This module was created by Tim Zuehlsdorff in 2015.            !
!================================================================!

module lr_tddft_utils

  use constants, only: DP

  implicit none

  private

  public :: lr_tddft_utils_oscillator
  public :: lr_tddft_utils_precond
  public :: lr_tddft_utils_precond_batch
  public :: lr_tddft_utils_precond_batch_RPA
  public :: lr_tddft_utils_eh_denskern_RPA
  public :: lr_tddft_utils_init_tddft_vecs
  public :: lr_tddft_utils_penalty_batch
  public :: lr_tddft_utils_penalty_batch_RPA
  public :: lr_tddft_utils_calc_mlwfs

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_utils_oscillator(oscillator, omega, &
       cond_ngwf_basis, cond_rep, val_ngwf_basis,val_rep,proj_basis,&
       nl_projectors,kernel,num_states, mdl,joint_evecs,joint_evals)

    !========================================================================!
    ! Routine calculates oscillator strength of a given excitation from its  !
    ! response_density and its energy. This routine can be called for both   !
    ! RPA and TDA, but for RPA the (X+Y) kernel is used, which reduces to    !
    ! just X in the Tamm-Dancoff approximation. Furthermore, the routine is  !
    ! capable of computing the oscillator strength both in momentum and in   !
    ! position representation. However, currently the momentum               !
    ! representation only works if this is not a PAW calculation and the     !
    ! joint_rep is used to represent the conduction space.                   !
    !------------------------------------------------------------------------!
    ! Written by Tim Zuehlsdorff in March 2016 using earlier routines from   !
    ! lr_tddft_mod.                                                          !
    ! Modified for embedding by Joseph Prentice, July 2018                   !
    !========================================================================!

    use comms, only: pub_on_root
    use function_basis, only: FUNC_BASIS
    use constants, only: stdout, HARTREE_IN_DEBYE
    use dense, only: DEM,dense_get_element,dense_create,dense_destroy,&
        dense_product,dense_copy,dense_convert,&
        dense_put_element
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_num_spins, pub_lr_tddft_mom_mat_els,&
         pub_any_nl_proj, pub_lr_tddft_joint_set, pub_paw
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_copy, sparse_embed_trace
    use optics, only: optics_pos_mat_els, optics_grad_mat_els
    use utils, only: utils_alloc_check, utils_dealloc_check,&
         utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: num_states
    real(kind=DP), intent(inout) :: oscillator(num_states)
    real(kind=DP), intent(in) :: omega(num_states)
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis(:)
    type(NGWF_REP), intent(in) :: cond_rep
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(:)
    type(SPAM3_EMBED), intent(inout) :: kernel(num_states,pub_num_spins)
    type(DEM),optional, intent(in) :: joint_evecs(pub_num_spins)
    real(kind=DP), optional, intent(in) :: joint_evals(:,:)!cond_ngwf_num,pub_num_spins)
    type(MODEL), intent(in) :: mdl

    ! Local Variables
    type(SPAM3_EMBED), allocatable, dimension(:) :: rmat,eff_joint_kernel
    integer :: ierr, icount, jcount, is,cond_count,val_count
    real(kind=DP) :: dipole_mom, temp_val
    real(kind=DP) :: dipole_xyz(3)
    type(SPAM3_EMBED) :: joint_proj, SP, eff_kernel,eff_rmat
    type(DEM), allocatable, dimension(:) :: evecs_dense_c
    type(DEM), allocatable, dimension(:,:) :: rmat_dense,rmat_dense_real
    type(DEM) :: overlap_dense,temp_dense
    real(kind=DP) :: temp, trace
    COMPLEX(kind=DP) :: temp_z
    ! jcap: embedding local variables
    integer :: cond_ngwf_num

    ! calculate oscillator strength in position representation
    if(.not. (pub_lr_tddft_mom_mat_els .and. present(joint_evecs) &
          .and. present(joint_evals))) then
       ! allocate data
       allocate(rmat(3), stat=ierr)
       call utils_alloc_check('lr_tddft_utils_oscillator','rmat', ierr)
       do icount=1, 3
          call sparse_embed_create(rmat(icount), cond_rep%cross_overlap)
       enddo

       ! calculate rmat
       call optics_pos_mat_els(rmat, val_rep, val_ngwf_basis, &
            cond_rep, cond_ngwf_basis, cond_rep%cross_overlap, &
            cond_rep%ngwf_cross_overlap,proj_basis,mdl)

       ! header for xyz transition dipole moment
       if (pub_on_root) write(stdout,'(/a80)')'|Excitation|    d_x (in D) &
               &  |     d_y (in D)  |    d_z (in D)   '

       ! construct oscillator strength from the dipole moment
       ! of each transition density
       do jcount=1, num_states

          dipole_mom=0.0_DP
          do icount=1, 3
             temp_val=0.0_DP
             do is=1, pub_num_spins
                call sparse_embed_trace(trace, kernel(jcount,is), rmat(icount))
                temp_val=temp_val+trace
             enddo
             dipole_mom=dipole_mom+temp_val*temp_val
             ! convert to Debye units
             dipole_xyz(icount)=temp_val*HARTREE_IN_DEBYE
          enddo
          ! print transition dipole moment
          if (pub_on_root) write(stdout,'(i24,e18.8,e18.8,e18.8)') &
                jcount,dipole_xyz(1),dipole_xyz(2),dipole_xyz(3)

          ! factor of 2 because of spin.
          if(pub_num_spins==1) then
             oscillator(jcount)=4.0_DP/3.0_DP*dipole_mom*omega(jcount)
          else
             oscillator(jcount)=2.0_DP/3.0_DP*dipole_mom*omega(jcount)
          endif
       enddo
    else ! momentum representation
       ! this only works if joint rep is used for cond rep
       ! build in a couple of utils assert checks here

       call utils_assert(pub_lr_tddft_joint_set, &
            'Computation of oscillator strength in momentum rep is&
             &only possible if pub_lr_tddft_joint_set=T')
       call utils_assert(.not. pub_paw, &
            'Computation of oscillator strenth in momentum rep not &
            &implemented for PAW')

       ! jcap: set number of conduction NGWFs in total
       if (present(joint_evals)) cond_ngwf_num=size(joint_evals(:,1))

       allocate(rmat(3), stat=ierr)
       call utils_alloc_check('lr_tddft_utils_oscillator','rmat', ierr)
       allocate(eff_joint_kernel(pub_num_spins),stat=ierr)
       call utils_alloc_check('lr_tddft_utils_oscillator',&
           'eff_joint_kernel',ierr)
       allocate(evecs_dense_c(pub_num_spins),stat=ierr)
       call utils_alloc_check('lr_tddft_utils_oscillator','evecs_dense_c',ierr)
       allocate(rmat_dense(pub_num_spins,3), stat=ierr)
       call utils_alloc_check('lr_tddft_utils_oscillator','rmat_dense',ierr)
       allocate(rmat_dense_real(pub_num_spins,3), stat=ierr)
       call utils_alloc_check('lr_tddft_utils_oscillator','rmat_dense_real',ierr)

       do icount=1, 3
          if(pub_any_nl_proj) then
             call sparse_embed_create(rmat(icount), cond_rep%nonlocpot,&
                  iscmplx=.true.)
          else
             call sparse_embed_create(rmat(icount),cond_rep%overlap,iscmplx=.true.)
          endif
       enddo
       ! ceate data structures to generate effective joint response kernel
       call sparse_embed_create(joint_proj,cond_rep%cross_overlap,cond_rep%inv_overlap)
       do is=1, pub_num_spins
          call sparse_embed_create(eff_joint_kernel(is),kernel(1,is),joint_proj)
       enddo

       call sparse_embed_product(joint_proj,cond_rep%cross_overlap,cond_rep%inv_overlap)
       call sparse_embed_create(SP,cond_rep%overlap,eff_joint_kernel(1))
       call sparse_embed_create(eff_kernel,eff_joint_kernel(1))
       call sparse_embed_create(eff_rmat,eff_kernel)

       ! create mom matrix elements
       call optics_grad_mat_els(rmat,cond_rep,cond_ngwf_basis,proj_basis,&
            nl_projectors,mdl)

       do is=1, pub_num_spins
          call dense_create(evecs_dense_c(is),cond_ngwf_num,cond_ngwf_num,&
               iscmplx=.true.)
          call dense_copy(evecs_dense_c(is),joint_evecs(is))
       enddo

       ! convert rmat into KS space.
       do is=1, pub_num_spins
          do icount=1,3
             call dense_create(rmat_dense(is,icount),cond_ngwf_num,cond_ngwf_num,&
                   iscmplx=.true.)
             call dense_create(rmat_dense_real(is,icount),cond_ngwf_num,cond_ngwf_num,&
                   iscmplx=.false.)
          enddo
       enddo
       call dense_create(temp_dense,cond_ngwf_num,cond_ngwf_num,&
                iscmplx=.true.)
       call dense_create(overlap_dense,cond_ngwf_num,cond_ngwf_num,&
                iscmplx=.true.)

       ! construct effective rmat in KS space
       do is=1, pub_num_spins
          do icount=1,3
             call dense_convert(temp_dense,rmat(icount))
             call dense_product(overlap_dense,temp_dense,evecs_dense_c(is),opA='T')
             call dense_product(rmat_dense(is,icount),evecs_dense_c(is),overlap_dense,&
                  opA='T')
          enddo
       enddo
       call dense_destroy(temp_dense)
       call dense_destroy(overlap_dense)
       call dense_create(temp_dense,cond_ngwf_num,cond_ngwf_num,iscmplx=.false.)
       call dense_create(overlap_dense,cond_ngwf_num,cond_ngwf_num,iscmplx=.false.)

       ! now create effective R_mat matrix elements, just as in optics mod
       do icount=1,3
          do is=1, pub_num_spins
             ! loop over individual matrix elements:
             do val_count=1,cond_ngwf_num
                do cond_count=1,cond_ngwf_num
                  ! if(val_count<val_rep%n_occ(is,1)+1 .and. cond_count>val_rep%n_occ(is,1)) then
                   if(val_count/=cond_count) then
                      call dense_get_element(temp_z,rmat_dense(is,icount),&
                         val_count,cond_count)
                      temp=real(temp_z)/(joint_evals(cond_count,is)-&
                          joint_evals(val_count,is))
                   else
                      temp=0.0_DP
                   endif
                   call dense_put_element(temp,rmat_dense_real(is,icount),&
                       val_count,cond_count) ! transpose
                enddo
             enddo
          enddo
       enddo

       ! convert rmat_dense_real back into NGWF space. This can be done by multiplying
       ! from left and right by the Kohn-Sham eigenvectors in joint rep.
       do icount=1,3
          do is=1, pub_num_spins
             call dense_product(temp_dense,rmat_dense_real(is,icount),joint_evecs(is),opB='T')
             call dense_product(rmat_dense_real(is,icount),joint_evecs(is),temp_dense)
          enddo
       enddo

       ! loop over all states of which we want to compute the oscillator strength
       do icount=1, num_states
          ! loop over spins
          do is=1, pub_num_spins
             call sparse_embed_product(eff_joint_kernel(is),kernel(icount,is),&
                  joint_proj)
          enddo

          dipole_mom=0.0_DP
          do jcount=1,3
             temp_val=0.0_DP
             do is=1, pub_num_spins
                ! project SPS into KS space.
                call sparse_embed_copy(eff_kernel,eff_joint_kernel(is))
                call dense_convert(eff_rmat,rmat_dense_real(is,jcount))
                ! eff rmat is a contravariant quantity. Convert into covariant
                call sparse_embed_product(SP,cond_rep%overlap,eff_rmat)
                call sparse_embed_product(eff_rmat,SP,cond_rep%overlap)

                call sparse_embed_trace(trace,eff_kernel,eff_rmat)
                temp_val=temp_val+trace
             enddo
             dipole_mom=dipole_mom+temp_val*temp_val
          enddo

          ! successfully calculated matrix element for this transition
          if(pub_num_spins==1) then
             oscillator(icount)=4.0_DP/3.0_DP*dipole_mom*omega(icount)
          else
             oscillator(icount)=2.0_DP/3.0_DP*dipole_mom*omega(icount)
          endif
       enddo


       ! deallocate data
       do is=1, pub_num_spins
          call sparse_embed_destroy(eff_joint_kernel(is))
          call dense_destroy(evecs_dense_c(is))
          do icount=1,3
             call dense_destroy(rmat_dense(is,icount))
             call dense_destroy(rmat_dense_real(is,icount))
          enddo
       enddo
       call sparse_embed_destroy(eff_kernel)
       call sparse_embed_destroy(SP)
       call sparse_embed_destroy(eff_rmat)
       call dense_destroy(temp_dense)
       call dense_destroy(overlap_dense)
       call sparse_embed_destroy(joint_proj)
       deallocate(eff_joint_kernel,stat=ierr)
       call utils_dealloc_check('lr_tddft_utils_oscillator',&
           'eff_joint_kernel',ierr)
       deallocate(evecs_dense_c,stat=ierr)
       call utils_dealloc_check('lr_tddft_utils_oscillator',&
           'evecs_dense_c',ierr)
       deallocate(rmat_dense,stat=ierr)
       call utils_dealloc_check('lr_tddft_utils_oscillator',&
           'rmat_dense',ierr)
       deallocate(rmat_dense_real,stat=ierr)
       call utils_dealloc_check('lr_tddft_utils_oscillator',&
           'rmat_dense_real',ierr)
    endif

    ! deallocate data
    do icount=1, 3
       call sparse_embed_destroy(rmat(icount))
    enddo
    deallocate(rmat, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_oscillator','rmat', ierr)

  end subroutine lr_tddft_utils_oscillator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_utils_calc_mlwfs(vecs,cond_ngwf_basis,val_ngwf_basis,&
     cond_rep,val_rep, mdl, proj_basis, excitations,num_states,vecs_q)

     !====================================================================!
     ! Subroutine computes a unitary transformation matrix that rotates   !
     ! a set of converged (RPA) excitons into a maximally localised rep.  !
     ! Furthermore, it expresses the excitation energies in this site     !
     ! basis. The construction of the localised representation is         !
     ! equivalent to the construction of maximally localised wannier      !
     ! functions in real space for isolated systems. Thus the algorithm   !
     ! presented here is only appropriate for excited states in isolated  !
     ! systems.                                                           !
     ! Modified for embedding by Joseph Prentice, July 2018               !
     !====================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, ANGSTROM, HARTREE_IN_DEBYE
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use linalg, only: linalg_mat_mul_serial
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_product, sparse_embed_create,&
         sparse_embed_destroy, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_copy, sparse_embed_transpose_structure, sparse_embed_trace
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_destroy, &
         sparse_scale, sparse_trace
    use rundat, only: pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check
    use visual, only: visual_scalarfield
    use optics, only: optics_pos_mat_els

    ! Arguments
    integer, intent(in) :: num_states
    type(SPAM3_EMBED), intent(inout) :: vecs(num_states,pub_num_spins)
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type(NGWF_REP), intent(in) :: cond_rep
    type(NGWF_REP), intent(in) :: val_rep
    real(kind=DP), intent(in) :: excitations(num_states)
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(MODEL), intent(in) :: mdl
    type(SPAM3_EMBED), optional, intent(inout) :: vecs_q(num_states,pub_num_spins)
      ! vecs q holds the q vector. this is only present in a RPA calculation

    ! local variables
    real(kind=DP), allocatable, dimension(:,:) :: dipole_xyz
    character(50) :: filename, tempfile
    integer :: icount, jcount, counter, MAX_IT,ierr,rpa_counter,&
       rpa_num, is, isub, jsub
    real(kind=DP) :: step_length, tolerance, TOL, trace
    real(kind=DP) :: f0,f1,f2,a,b,c,lambda, f_new
    real(kind=DP) :: trial_step1,trial_step2, max_step_length
    type(SPAM3_EMBED), allocatable, dimension(:) :: temp_dkn
    real(kind=DP), allocatable, dimension(:,:,:) :: xyz_0_mat
    real(kind=DP), allocatable, dimension(:,:) :: grad_mat
    real(kind=DP), allocatable, dimension(:,:) :: delta_W
    real(kind=DP), allocatable, dimension(:,:) :: temp_u
    real(kind=DP), allocatable, dimension(:,:) :: u_mat
    real(kind=DP), allocatable, dimension(:,:) :: temp_mat
    real(kind=DP), allocatable, dimension(:,:,:,:,:,:) :: trans_dens_array
    real(kind=DP), allocatable, dimension(:,:,:,:) :: sub_trans_dens_array
    type(SPAM3_EMBED), allocatable, dimension(:) :: rmat
    type(SPAM3), allocatable, dimension(:) :: temp_spam

    ! check if this is an RPA or TDA calculation
    if(present(vecs_q)) then
       rpa_num=2
    else
       rpa_num=1
    endif

    ! first allocate trans_dens_array and xyz_0_mat and compute it
    allocate(trans_dens_array(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
         mdl%fine_grid%max_slabs12,pub_num_spins,num_states,rpa_num), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_calc_mlwfs',&
         'trans_dens_array',ierr)
    allocate(sub_trans_dens_array(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
         mdl%fine_grid%max_slabs12,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_calc_mlwfs',&
         'sub_trans_dens_array',ierr)
    allocate(xyz_0_mat(num_states,num_states,3), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_calc_mlwfs','xyz_0_mat',ierr)
    allocate(temp_dkn(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_calc_mlwfs','temp_dkn',ierr)
    !do is=1, pub_num_spins
    !   call sparse_create(temp_dkn(is),vecs(1,is))
    !enddo
    allocate(temp_spam(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_calc_mlwfs','temp_spam',ierr)

    ! calculate all transition densities for all spins
    trans_dens_array=0.d0
    do icount=1, num_states
       if(rpa_num==1) then ! this is a TDA calculation and vecs contains X

          ! jcap: loop over regions and sum up densities
          do isub=1,mdl%nsub
             do jsub=1,mdl%nsub

                sub_trans_dens_array=0.d0

                do is=1, pub_num_spins
                   call sparse_create(temp_spam(is),vecs(icount,is)%m(isub,jsub))
                   call sparse_copy(temp_spam(is),vecs(icount,is)%m(isub,jsub))
                enddo

                call density_on_grid(sub_trans_dens_array,mdl%fine_grid,&
                     mdl%dbl_grid,mdl%cell,mdl%fftbox,temp_spam,&
                     cond_rep%ngwf_cross_overlap_tr%m(isub,jsub),&
                     cond_rep%ngwfs_on_grid(isub),cond_ngwf_basis(isub),&
                     val_rep%ngwfs_on_grid(jsub),val_ngwf_basis(jsub))

                trans_dens_array(:,:,:,:,icount,rpa_num)=&
                     trans_dens_array(:,:,:,:,icount,rpa_num)+sub_trans_dens_array

                do is=1,pub_num_spins
                   call sparse_destroy(temp_spam(is))
                end do

             end do
          end do
       else ! this is an RPA calculation

          ! jcap: loop over regions and sum up densities
          do isub=1,mdl%nsub
             do jsub=1,mdl%nsub

                sub_trans_dens_array=0.d0

                do is=1, pub_num_spins
                   call sparse_create(temp_spam(is),vecs(icount,is)%m(isub,jsub))
                   call sparse_copy(temp_spam(is),vecs(icount,is)%m(isub,jsub))
                enddo

                call density_on_grid(sub_trans_dens_array,mdl%fine_grid,&
                     mdl%dbl_grid,mdl%cell,mdl%fftbox,temp_spam,&
                     cond_rep%ngwf_cross_overlap_tr%m(isub,jsub),&
                     cond_rep%ngwfs_on_grid(isub),cond_ngwf_basis(isub),&
                     val_rep%ngwfs_on_grid(jsub),val_ngwf_basis(jsub))

                trans_dens_array(:,:,:,:,icount,1)=&
                     trans_dens_array(:,:,:,:,icount,1)+sub_trans_dens_array

                do is=1,pub_num_spins
                   call sparse_destroy(temp_spam(is))
                end do

                sub_trans_dens_array=0.d0

                ! now construct Y
                do is=1, pub_num_spins
                   call sparse_create(temp_spam(is),vecs_q(icount,is)%m(isub,jsub))
                   call sparse_copy(temp_spam(is),vecs_q(icount,is)%m(isub,jsub))
                enddo

                call density_on_grid(sub_trans_dens_array,mdl%fine_grid,&
                     mdl%dbl_grid,mdl%cell,mdl%fftbox,temp_spam,&
                     cond_rep%ngwf_cross_overlap_tr%m(isub,jsub),&
                     cond_rep%ngwfs_on_grid(isub),cond_ngwf_basis(isub),&
                     val_rep%ngwfs_on_grid(jsub),val_ngwf_basis(jsub))

                trans_dens_array(:,:,:,:,icount,2)=&
                     trans_dens_array(:,:,:,:,icount,2)+sub_trans_dens_array

                do is=1,pub_num_spins
                   call sparse_destroy(temp_spam(is))
                end do
             end do
          end do

       endif
    enddo

    ! now calculate xyz_0_mat
    do icount=1, num_states
       do jcount=1, num_states
          if(rpa_num==1) then
             call internal_calc_pos_mat_int(xyz_0_mat(icount,jcount,:),&
                  trans_dens_array(:,:,:,:,icount,1),&
                  trans_dens_array(:,:,:,:,jcount,1),mdl%fine_grid)
          else
             call internal_calc_pos_mat_int(xyz_0_mat(icount,jcount,:),&
                  trans_dens_array(:,:,:,:,icount,1),&
                  trans_dens_array(:,:,:,:,jcount,2),mdl%fine_grid)
          endif
       enddo
    enddo

    ! get rid of all unnecessary dara strucutres
    deallocate(trans_dens_array, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_calc_mlwfs','trans_dens_array',ierr)
    deallocate(sub_trans_dens_array, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_calc_mlwfs','sub_trans_dens_array',ierr)

    ! now can start actual optimisation procedure of u_mat. First, allocate data
    ! structures
    allocate(u_mat(num_states,num_states), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_calc_mlwfs','u_mat',ierr)
    allocate(temp_u(num_states,num_states), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_calc_mlwfs','temp_u',ierr)
    allocate(grad_mat(num_states,num_states), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_calc_mlwfs','grad_mat',ierr)
    allocate(delta_W(num_states,num_states), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_calc_mlwfs','delta_W',ierr)
    allocate(temp_mat(num_states,num_states), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_calc_mlwfs','u_mat',ierr)

    ! initialise optimisation parameters
    step_length=0.2_DP
    max_step_length=60.0_DP
    MAX_IT=100000
    counter=1
    TOL=10E-12_DP
    tolerance=1.0_DP

    ! initialise u_mat to identity matrix
    do icount=1, num_states
       do jcount=1,num_states
          if(icount==jcount) then
             u_mat(icount,jcount)=1.0_DP
             temp_u(icount,jcount)=1.0_DP
          else
             u_mat(icount,jcount)=0.0_DP
             temp_u(icount,jcount)=0.0_DP
          endif
       enddo
    enddo

    ! compute intital localisation functional value
    call internal_calc_functional(f0,u_mat,xyz_0_mat,num_states)

    ! compute initial gradient
    call internal_calc_gradient(grad_mat,u_mat,xyz_0_mat,num_states)
    if(pub_on_root) write(stdout,'(a)') 'Computing maximally localised wannier funcs'
    if(pub_on_root) write(stdout,'(a)') 'Iter  |  Func_val   |   change    '

    ! start optimisation loop
    do while (counter<MAX_IT .and. abs(tolerance)>TOL .and. tolerance>0.0_DP)
       ! take a trial step. compute delta W for trial step and set temp_mat=u
       trial_step1=step_length
       do icount=1, num_states
          do jcount=1, num_states
             delta_W(icount,jcount)=trial_step1*grad_mat(icount,jcount)
             temp_mat(icount,jcount)=u_mat(icount,jcount)
          enddo
       enddo
       ! evaluate u_mat for that step and func_val for that step
       call internal_update_u_mat(temp_mat,delta_W,num_states)
       call internal_calc_functional(f1,temp_mat,xyz_0_mat,num_states)

       ! if functional increases in value, make sure to step into correct
       ! direction
       if(f1>f0) then
          trial_step1=-1.0_DP*step_length
          do icount=1, num_states
             do jcount=1, num_states
                delta_W(icount,jcount)=trial_step1*grad_mat(icount,jcount)
                temp_mat(icount,jcount)=u_mat(icount,jcount)
             enddo
          enddo
          ! evaluate u_mat for that step and func_val for that step
          call internal_update_u_mat(temp_mat,delta_W,num_states)
          call internal_calc_functional(f1,temp_mat,xyz_0_mat,num_states)
       endif

       if(trial_step1>0.0_DP) then
         trial_step2=2.0_DP*step_length
       else
         trial_step2=-2.0_DP*step_length
       endif

       ! take a second trial step. compute delta@ and set temp_mat=u
       do icount=1, num_states
          do jcount=1, num_states
             delta_W(icount,jcount)=trial_step2*grad_mat(icount,jcount)
             temp_mat(icount,jcount)=u_mat(icount,jcount)
          enddo
       enddo
       ! evaluate u_mat for that step and func_val for that step
       call internal_update_u_mat(temp_mat,delta_W,num_states)
       call internal_calc_functional(f2,temp_mat,xyz_0_mat,num_states)

       ! have evaluated the functional at 3 points. Fit a parabola and
       ! find minimum to calculate ideal step length
       c=f0
       a=(f1-f2*trial_step1/trial_step2-f0*(1.0_DP-trial_step1/trial_step2))/&
         (trial_step1*trial_step1-trial_step1*trial_step2)
       b=(f1-f0)/trial_step1-a*trial_step1
       ! now compute ideal step length
       ! first check if parabola has the correct sign. Are we looking for
       ! a minimum?
       if(a<0.0_DP) then
          ! check if f1 or f2 are smaller than f0
          if(f1<f0 .and. f2<f0) then
             if(f2<f1) then
                lambda=trial_step2
             else
                lambda=trial_step1
             endif
          else if(f1<f0) then
             lambda=trial_step1
          else if(f2<f0) then
             lambda=trial_step2
          else
             if(pub_on_root) write(stdout,'(a)') 'Error in line min. Func val increases'
             lambda=0.5_DP*trial_step1
          endif
       else
          lambda=-b/(2.0_DP*a)
          if(abs(lambda)>max_step_length) then
             if(lambda<0.0_DP) then
                lambda=-1.0_DP*max_step_length
             else
                lambda=max_step_length
             endif
          endif
       endif

       ! compute correct delta_W
       do icount=1, num_states
          do jcount=1, num_states
             delta_W(icount,jcount)=lambda*grad_mat(icount,jcount)
          enddo
       enddo
       ! update u_mat
       call internal_update_u_mat(u_mat,delta_W,num_states)
       call internal_calc_functional(f_new,u_mat,xyz_0_mat,num_states)

       ! compute new gradient for updated u_mat
       call internal_calc_gradient(grad_mat,u_mat,xyz_0_mat,num_states)

       ! if this new step actually yields an increase in functional value
       ! check f1 and f2
       if(f0<f_new) then
          if(pub_on_root) write(stdout,'(a)') 'Func val increased!!!!'
          ! do a safe step
          if(f2<f0) then
             lambda=trial_step2
          else if(f1<f0) then
             lambda=trial_step1
          else
             if(pub_on_root) write(stdout,'(a)') 'Failed to fix step!'
          endif

          do icount=1, num_states
             do jcount=1, num_states
                u_mat(icount,jcount)=temp_u(icount,jcount)
                delta_W(icount,jcount)=lambda*grad_mat(icount,jcount)
             enddo
          enddo
          call internal_update_u_mat(u_mat,delta_W,num_states)
          call internal_calc_functional(f_new,u_mat,xyz_0_mat,num_states)

          ! compute new gradient for updated u_mat
          call internal_calc_gradient(grad_mat,u_mat,xyz_0_mat,num_states)

       endif
       ! set temp_u
       do icount=1, num_states
          do jcount=1, num_states
             temp_u(icount,jcount)=u_mat(icount,jcount)
          enddo
       enddo

       ! set tolerance
       ! and new f0 and update counter
       tolerance=f0-f_new
       f0=f_new
       counter=counter+1
    enddo

    ! print final iteration number, residue etc.
    if(pub_on_root) write(stdout,'(i6,f12.8,f12.8)') counter,f0,tolerance

    ! now ideal u_mat should have been found. As a first step, compute
    ! the TDDFT hamiltonian in the representation of the maximally localised
    ! states. This can be done by computing U^T*diag(omega)*U, where
    ! diag(omega) is the diagonal matrix of eigenvalues.
    do icount=1,num_states
       do jcount=1, num_states
          if(icount==jcount) then
             delta_W(icount,jcount)=excitations(icount)
          else
             delta_W(icount,jcount)=0.0_DP
          endif
       enddo
    enddo

    call linalg_mat_mul_serial(grad_mat,delta_W,&
                u_mat)
    call linalg_mat_mul_serial(temp_mat,u_mat,grad_mat,opA='T')

    if(pub_on_root) write(stdout,'(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    if(pub_on_root) write(stdout,'(a)') 'Maximally Localised Wannier function analysis:'
    if(pub_on_root) write(stdout,'(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    do icount=1,num_states
       do jcount=1,num_states
          if(pub_on_root) write(stdout,'(f12.8)', advance='no') temp_mat(icount,jcount)
       enddo
       if(pub_on_root) write(stdout,*) '' ! line break
    enddo
    if(pub_on_root) write(stdout,'(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

    ! now construct effective density kernels for each localised state and print
    ! the transition density associated with it
    ! reallocate trans_mat_array for this
    ! Also compute transition dipole moment of localised states
    ! If this is an RPA calculation, only do this for Q, which is the
    ! effective X+Y vector
    allocate(trans_dens_array(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
         mdl%fine_grid%max_slabs12,pub_num_spins,1,1), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_calc_mlwfs',&
         'trans_dens_array',ierr)
    allocate(sub_trans_dens_array(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
         mdl%fine_grid%max_slabs12,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_calc_mlwfs',&
         'sub_trans_dens_array',ierr)
    allocate(dipole_xyz(3,num_states), stat=ierr)
    call utils_alloc_check('lr_tddft_uitls_calc_mlwfs','dipole_xyz',ierr)
    allocate(rmat(3), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_calc_mlwfs','rmat',ierr)
    do icount=1, 3
       call sparse_embed_create(rmat(icount), cond_rep%cross_overlap)
    enddo

    ! calculate rmat. This only works in position rep
    call optics_pos_mat_els(rmat, val_rep, val_ngwf_basis, &
         cond_rep, cond_ngwf_basis, cond_rep%cross_overlap, &
         cond_rep%ngwf_cross_overlap,proj_basis,mdl)

    trans_dens_array=0.d0
    dipole_xyz=0.d0
    do icount=1,num_states
       do jcount=1,num_states
          do is=1, pub_num_spins
             if(jcount==1) then
                if(rpa_num==1) then
                   ! jcap: only create if jcount is 1
                   call sparse_embed_create(temp_dkn(is),vecs(jcount,is))
                   call sparse_embed_copy(temp_dkn(is),vecs(jcount,is))
                   call sparse_embed_scale(temp_dkn(is),u_mat(jcount,icount))
                else
                   ! jcap: only create if jcount is 1
                   call sparse_embed_create(temp_dkn(is),vecs_q(jcount,is))
                   call sparse_embed_copy(temp_dkn(is),vecs_q(jcount,is))
                   call sparse_embed_scale(temp_dkn(is),u_mat(jcount,icount))
                endif
             else
                if(rpa_num==1) then
                   call sparse_embed_axpy(temp_dkn(is),vecs(jcount,is),&
                        u_mat(jcount,icount))
                else
                   call sparse_embed_axpy(temp_dkn(is),vecs_q(jcount,is),&
                        u_mat(jcount,icount))
                endif
             endif
          enddo
       enddo

       ! jcap: need to loop over regions here to sum up densities
       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub

             sub_trans_dens_array=0.d0

             do is=1,pub_num_spins
                call sparse_create(temp_spam(is),temp_dkn(is)%m(isub,jsub))
                call sparse_copy(temp_spam(is),temp_dkn(is)%m(isub,jsub))
             end do

             ! now temp_dkn(1) contains denskern of localised exciton icount.
             ! construct its density and print it.
             call density_on_grid(sub_trans_dens_array,mdl%fine_grid,&
                  mdl%dbl_grid,mdl%cell,mdl%fftbox,temp_spam,&
                  cond_rep%ngwf_cross_overlap_tr%m(isub,jsub),&
                  cond_rep%ngwfs_on_grid(isub),cond_ngwf_basis(isub),&
                  val_rep%ngwfs_on_grid(jsub),val_ngwf_basis(jsub))

             trans_dens_array(:,:,:,:,1,1)=trans_dens_array(:,:,:,:,1,1)+&
                  sub_trans_dens_array

             ! jcap: destroy temporary density kernel here
             do is=1,pub_num_spins
                call sparse_destroy(temp_spam(is))
             end do

          end do
       end do

       ! jcap: need to construct the new density kernel

       ! jcap: calculate dipole here
       do jcount=1, 3
          do is=1, pub_num_spins
             call sparse_embed_trace(trace,rmat(jcount),temp_dkn(is))
             dipole_xyz(jcount,icount)=dipole_xyz(jcount,icount)+trace
          enddo
          dipole_xyz(jcount,icount)=dipole_xyz(jcount,icount)*HARTREE_IN_DEBYE
       enddo

       if(pub_num_spins==2) then
          trans_dens_array(:,:,:,1,1,1)=trans_dens_array(:,:,:,1,1,1)+&
               trans_dens_array(:,:,:,2,1,1)  ! add up spin channels
       endif

       ! write density
       write(filename,*) icount
       write(tempfile,*) '_TDDFT_MLWF_dens_'
       write(filename,*) trim(adjustl(tempfile))//trim(adjustl(filename))

       call visual_scalarfield(trans_dens_array(:,:,:,1,1,1),mdl%fine_grid,&
          mdl%cell, 'MLWF density (in e/ang^3) for:',filename,mdl%elements,&
          ANGSTROM**3)

    enddo

    ! print dipole mom
    ! header for xyz transition dipole moment
    if (pub_on_root) write(stdout,'(/a80)')'|Site |    d_x (in D) &
         &  |     d_y (in D)  |    d_z (in D)   '
    do icount=1, num_states
       ! print transition dipole moment
       if (pub_on_root) write(stdout,'(i24,e18.8,e18.8,e18.8)') &
             icount,dipole_xyz(1,icount),dipole_xyz(2,icount),&
             dipole_xyz(3,icount)
    enddo

    do icount=1,3
       call sparse_embed_destroy(rmat(icount))
    end do
    deallocate(rmat,stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_calc_mlwfs','rmat',ierr)
    deallocate(dipole_xyz,stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_calc_mlwfs',&
      'dipole_xyz',ierr)
    deallocate(trans_dens_array,stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_calc_mlwfs','trans_dens_array',ierr)
    deallocate(sub_trans_dens_array,stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_calc_mlwfs','sub_trans_dens_array',ierr)
    if(pub_on_root) write(stdout,'(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

    ! at the end of the routine, deallocate remaining data structures
    do is=1, pub_num_spins
       call sparse_embed_destroy(temp_dkn(is))
    enddo
    deallocate(temp_dkn,stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_calc_mlwfs','temp_dkn',ierr)
    deallocate(u_mat, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_calc_mlwfs','u_mat',ierr)
    deallocate(grad_mat, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_calc_mlwfs','grad_mat',ierr)
    deallocate(delta_W, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_calc_mlwfs','delta_W',ierr)
    deallocate(temp_mat, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_calc_mlwfs','u_mat',ierr)

  contains

  !-----------------------------------------------------------------------!
  subroutine internal_calc_functional(func_val,u_mat,xyz_0_mat,mat_dim)

     !====================================================================!
     ! Subroutine computes the value of the localisation functional,      !
     ! which is given by f(X,Y,Z)=Tr[(X')^2+(Y')^2+(Z')^2], where X' is   !
     ! given by X'=X-diag(X).                                             !
     !====================================================================!

     use linalg, only: linalg_mat_mul_serial
     use utils, only: utils_alloc_check, utils_dealloc_check

     ! arguments
     integer, intent(in) :: mat_dim
     real(kind=DP), intent(inout) :: func_val
     real(kind=DP), intent(in) :: u_mat(mat_dim,mat_dim)
     real(kind=DP), intent(in) :: xyz_0_mat(mat_dim,mat_dim,3)

     ! local variables
     integer :: icount, jcount, xyz_count, ierr, is, rpa_counter
     real(kind=DP), allocatable, dimension(:,:,:) :: xyz_mat
     real(kind=DP), allocatable, dimension(:,:,:) :: xyz_dash_mat
     real(kind=DP), allocatable, dimension(:,:) :: temp_mat

     ! allocate variables.
     allocate(xyz_mat(mat_dim,mat_dim,3), stat=ierr)
     call utils_alloc_check('internal_calc_functional','xyz_mat',ierr)
     allocate(xyz_dash_mat(mat_dim,mat_dim,3), stat=ierr)
     call utils_alloc_check('internal_calc_functional','xyz_dash_mat',ierr)
     allocate(temp_mat(mat_dim,mat_dim), stat=ierr)
     call utils_alloc_check('internal_calc_functional','temp_mat',ierr)

     func_val=0.0_DP

     ! create xyz_mat
     do xyz_count=1,3
        call linalg_mat_mul_serial(temp_mat,xyz_0_mat(:,:,xyz_count),&
             u_mat)
        call linalg_mat_mul_serial(xyz_mat(:,:,xyz_count),u_mat,temp_mat,&
             opA='T')
     enddo
     ! now create xyz_dash
     do xyz_count=1,3
        do icount=1,mat_dim
           do jcount=1, mat_dim
              if(icount==jcount) then
                 xyz_dash_mat(icount,jcount,xyz_count)=0.0_DP
              else
                 xyz_dash_mat(icount,jcount,xyz_count)=xyz_mat(icount,jcount,xyz_count)
              endif
           enddo
        enddo
     enddo

     ! now evaluate functional
     do xyz_count=1,3
        ! compute square of xyz_dahs
        call linalg_mat_mul_serial(temp_mat,xyz_dash_mat(:,:,xyz_count),&
                xyz_dash_mat(:,:,xyz_count))

        ! compute trace of result
        do icount=1,num_states
           func_val=func_val+temp_mat(icount,icount)
        enddo
     enddo

     ! deallocate variables
     deallocate(xyz_mat, stat=ierr)
     call utils_dealloc_check('internal_calc_functional','xyz_mat',ierr)
     deallocate(xyz_dash_mat, stat=ierr)
     call utils_dealloc_check('internal_calc_functional','xyz_dash_mat',ierr)
     deallocate(temp_mat, stat=ierr)
     call utils_dealloc_check('internal_calc_functional','temp_mat',ierr)

  end subroutine internal_calc_functional

  !-----------------------------------------------------------------------!
  subroutine internal_calc_gradient(grad_mat,u_mat,xyz_0_mat,mat_dim)

     !====================================================================!
     ! Subroutine computes the gradient of the localisation functional    !
     ! with respect to a change in the unitary matrix u. The gradient is  !
     ! computed as g=2([X',X_D]+[Y',Y_D]+[Z',Z_D]) where X_D=Diag(X) and  !
     ! X'=X-X_D                                                           !
     !====================================================================!

     use linalg, only: linalg_mat_mul_serial
     use utils, only: utils_alloc_check, utils_dealloc_check

     ! arguments
     integer, intent(in) :: mat_dim
     real(kind=DP), intent(inout) :: grad_mat(mat_dim,mat_dim)
     real(kind=DP), intent(in) :: u_mat(mat_dim,mat_dim)
     real(kind=DP), intent(in) :: xyz_0_mat(mat_dim,mat_dim,3)

     ! local variables
     integer :: icount, jcount, xyz_count,ierr
     real(kind=DP), allocatable, dimension(:,:,:) :: xyz_mat
     real(kind=DP), allocatable, dimension(:,:,:) :: xyz_D_mat
     real(kind=DP), allocatable, dimension(:,:,:) :: xyz_dash_mat
     real(kind=DP), allocatable, dimension(:,:) :: temp_mat

     ! allocate local variables
     allocate(xyz_mat(mat_dim,mat_dim,3), stat=ierr)
     call utils_alloc_check('internal_calc_gradient','xyz_mat',ierr)
     allocate(xyz_D_mat(mat_dim,mat_dim,3), stat=ierr)
     call utils_alloc_check('internal_calc_gradient','xyz_D_mat',ierr)
     allocate(xyz_dash_mat(mat_dim,mat_dim,3), stat=ierr)
     call utils_alloc_check('internal_calc_gradient','xyz_dash_mat',ierr)
     allocate(temp_mat(mat_dim,mat_dim), stat=ierr)
     call utils_alloc_check('internal_calc_gradient','temp_mat',ierr)

     ! first create xyz from xyz_0 matrix and overlap matrices
     do xyz_count=1,3
        call linalg_mat_mul_serial(temp_mat,xyz_0_mat(:,:,xyz_count),&
                u_mat)
        call linalg_mat_mul_serial(xyz_mat(:,:,xyz_count),u_mat,temp_mat,&
              opA='T')
     enddo

     ! now create xyz_D and xyz_dash
     do xyz_count=1,3
        do icount=1,mat_dim
           do jcount=1, mat_dim
              if(icount==jcount) then
                 xyz_D_mat(icount,jcount,xyz_count)=xyz_mat(icount,jcount,xyz_count)
                 xyz_dash_mat(icount,jcount,xyz_count)=0.0_DP
              else
                 xyz_D_mat(icount,jcount,xyz_count)=0.0_DP
                 xyz_dash_mat(icount,jcount,xyz_count)=xyz_mat(icount,jcount,xyz_count)
              endif
           enddo
        enddo
     enddo

     ! set gradient matrix to zero
     do icount=1,mat_dim
        do jcount=1, mat_dim
           grad_mat(icount,jcount)=0.0_DP
        enddo
     enddo

     ! now have all ingredients to evaluate gradient
     do xyz_count=1,3
        call linalg_mat_mul_serial(temp_mat,xyz_dash_mat(:,:,xyz_count),&
            xyz_D_mat(:,:,xyz_count))

        do icount=1,mat_dim
           do jcount=1, mat_dim
              grad_mat(icount,jcount)=grad_mat(icount,jcount)+2.0_DP*temp_mat(icount,jcount)
           enddo
        enddo

        call linalg_mat_mul_serial(temp_mat,xyz_D_mat(:,:,xyz_count),&
            xyz_dash_mat(:,:,xyz_count))

        do icount=1,mat_dim
           do jcount=1, mat_dim
              grad_mat(icount,jcount)=grad_mat(icount,jcount)-2.0_DP*temp_mat(icount,jcount)
           enddo
        enddo
     enddo

     ! deallocate mat
     deallocate(xyz_mat, stat=ierr)
     call utils_dealloc_check('internal_calc_gradient','xyz_mat',ierr)
     deallocate(xyz_D_mat, stat=ierr)
     call utils_dealloc_check('internal_calc_gradient','xyz_D_mat',ierr)
     deallocate(xyz_dash_mat, stat=ierr)
     call utils_dealloc_check('internal_calc_gradient','xyz_dash_mat',ierr)
     deallocate(temp_mat, stat=ierr)
     call utils_dealloc_check('internal_calc_gradient','temp_mat',ierr)

  end subroutine internal_calc_gradient

  !-----------------------------------------------------------------------!
  subroutine internal_update_u_mat(u_mat,delta_W,mat_dim)

     !====================================================================!
     ! Internal subroutine updates the unitary matrix u_mat according to  !
     ! U-->U*exp(deltaW). Note that u_mat and delta_W are both real       !
     ! matrices and delta_W is anti-hermitian. Thus to construct          !
     ! the exponential of delta_W, we solve for the eigenvalues and       !
     ! eigenvectors of the hermitian matrix H=i*deltaW. Using the eigen-  !
     ! vectors and eigenvalues it is then possible to compute exp(deltaW).!
     !====================================================================!

     use linalg, only: linalg_mat_mul_serial
     use utils, only: utils_alloc_check, utils_dealloc_check

     ! arguments
     integer, intent(in) :: mat_dim
     real(kind=DP), intent(inout) :: u_mat(mat_dim,mat_dim)
     real(kind=DP), intent(in) :: delta_W(mat_dim,mat_dim)

     ! local variables:
     complex(kind=DP), allocatable,dimension(:,:) :: H_mat
     complex(kind=DP), allocatable,dimension(:,:) :: diag_mat
     real(kind=DP), allocatable,dimension(:) :: evals
     complex(kind=DP), allocatable,dimension(:,:) :: expW_cmplx
     real(kind=DP), allocatable, dimension(:,:) :: expW_real
     real(kind=DP), allocatable, dimension(:,:) :: buffer_mat_real
     real(kind=DP), allocatable, dimension(:,:) :: buffer_mat_real2
     complex(kind=DP), allocatable, dimension(:,:) :: buffer_mat_cmplx
     integer :: icount, jcount, ierr, counter
     real(kind=DP) :: factorial, current_factor
     real(kind=DP) :: temp_real
     real(kind=DP) :: temp_cmplx

     ! create data structures
     allocate(H_mat(mat_dim,mat_dim), stat=ierr)
     call utils_alloc_check('internal_update_u_mat','H_mat',ierr)
     allocate(diag_mat(mat_dim,mat_dim), stat=ierr)
     call utils_alloc_check('internal_update_u_mat','diag_mat',ierr)
     allocate(expW_cmplx(mat_dim,mat_dim), stat=ierr)
     call utils_alloc_check('internal_update_u_mat','expW_cmplx',ierr)
     allocate(expW_real(mat_dim,mat_dim), stat=ierr)
     call utils_alloc_check('internal_update_u_mat','expW_realt',ierr)
     allocate(buffer_mat_real(mat_dim,mat_dim), stat=ierr)
     call utils_alloc_check('internal_update_u_mat','buffer_mat_real',ierr)
     allocate(buffer_mat_real2(mat_dim,mat_dim), stat=ierr)
     call utils_alloc_check('internal_update_u_mat','buffer_mat_real2',ierr)
     allocate(buffer_mat_cmplx(mat_dim,mat_dim), stat=ierr)
     call utils_alloc_check('internal_update_u_mat','buffer_mat_cmplx',ierr)
     allocate(evals(mat_dim), stat=ierr)
     call utils_alloc_check('internal_update_u_mat','evals',ierr)

     ! set exp_W to identity mat plus delta W to identity
     do icount=1,mat_dim
        do jcount=1,mat_dim
           if(icount==jcount) then
             expW_real(icount,jcount)=1.0_DP+delta_W(icount,jcount)
           else
             expW_real(icount,jcount)=delta_W(icount,jcount)
           endif
           buffer_mat_real(icount,jcount)=delta_W(icount,jcount)
        enddo
     enddo

     ! alternative construction: Power series
     factorial=1.0_DP
     current_factor=2.0_DP
     counter=2

     do while (counter<25)
        factorial=current_factor*factorial
        call linalg_mat_mul_serial(buffer_mat_real2,buffer_mat_real,delta_W)

        do icount=1,num_states
           do jcount=1,num_states
              expW_real(icount,jcount)=expW_real(icount,jcount)+&
                  buffer_mat_real2(icount,jcount)/factorial
           enddo
        enddo

        current_factor=current_factor+1.0_DP

        do icount=1,num_states
           do jcount=1,num_states
              buffer_mat_real(icount,jcount)=buffer_mat_real2(icount,jcount)
           enddo
        enddo

        counter=counter+1
     enddo

     ! compute update matrix
     call linalg_mat_mul_serial(buffer_mat_real,u_mat,expW_real)

     ! copy to u_mat
     do icount=1,mat_dim
        do jcount=1,mat_dim
           u_mat(icount,jcount)=buffer_mat_real(icount,jcount)
        enddo
     enddo

     ! destroy data structures
     deallocate(H_mat, stat=ierr)
     call utils_dealloc_check('internal_update_u_mat','H_mat',ierr)
     deallocate(diag_mat, stat=ierr)
     call utils_dealloc_check('internal_update_u_mat','diag_mat',ierr)
     deallocate(expW_cmplx, stat=ierr)
     call utils_dealloc_check('internal_update_u_mat','expW_cmplx',ierr)
     deallocate(expW_real, stat=ierr)
     call utils_dealloc_check('internal_update_u_mat','expW_realt',ierr)
     deallocate(buffer_mat_real, stat=ierr)
     call utils_dealloc_check('internal_update_u_mat','buffer_mat_real',ierr)
     deallocate(buffer_mat_cmplx, stat=ierr)
     call utils_dealloc_check('internal_update_u_mat','buffer_mat_cmplx',ierr)
     deallocate(evals, stat=ierr)
     call utils_dealloc_check('internal_update_u_mat','evals',ierr)

  end subroutine internal_update_u_mat

  !-----------------------------------------------------------------------!
  subroutine internal_calc_pos_mat_int(resulting_vec,den1,den2,grid)
    !=====================================================================!
    ! Subroutine computes the integral of two densities and the position  !
    ! operator directly on the simulation grid. Note that this is only    !
    ! valid in isolated systems.                                          !
    !=====================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_real_pt
    use comms, only: comms_reduce
    use rundat, only: pub_num_spins

    ! jd: Arguments
    real(kind=DP), intent(inout)        :: resulting_vec(3)
    type(GRID_INFO), intent(in)         :: grid
    real(kind=DP), intent(in)           :: den1(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(in)           :: den2(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_num_spins)

    ! jd: Internal variables
    integer :: i1, i2, i3,is
    real(kind=DP) :: integral
    real(kind=DP) :: real_point(3)
    integer :: p1, p2, p3
    integer :: axis

    !------------------------------------------------------------------------

    ! jd: Default upper bounds of integration
    p1 = grid%n1
    p2 = grid%n2
    p3 = grid%num_my_slabs12

    do axis=1,3
       resulting_vec(axis) = 0.0_DP
    enddo

    do is=1, pub_num_spins
       do i3=1, p3
          do i2=1, p2
             do i1=1, p1
                call cell_grid_real_pt(real_point,i1,i2,i3,grid)
                resulting_vec(1)=resulting_vec(1)+real_point(1)&
                       *den1(i1,i2,i3,is)*den2(i1,i2,i3,is)
                resulting_vec(2)=resulting_vec(2)+real_point(2)&
                       *den1(i1,i2,i3,is)*den2(i1,i2,i3,is)
                resulting_vec(3)=resulting_vec(3)+real_point(3)&
                       *den1(i1,i2,i3,is)*den2(i1,i2,i3,is)
             enddo
          enddo
       end do
    end do

    call comms_reduce('SUM',resulting_vec(1:3))

    do axis=1,3
       resulting_vec(axis)=resulting_vec(axis)*grid%weight
    enddo

  end subroutine internal_calc_pos_mat_int

  !-----------------------------------------------------------------------!

  end subroutine lr_tddft_utils_calc_mlwfs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_utils_penalty_batch(value,vecs,vecs_cov,cond_denskern, &
      val_denskern, cond_overlap,val_overlap,Svkv,kcSc,num_conv_states,num_states)
    !========================================================================!
    ! Subroutine computes total penalty value for a batch of vectors in TDA. !
    ! Modified for embedding by Joseph Prentice, July 2018                   !
    !========================================================================!

    use rundat, only: pub_num_spins,pub_lr_tddft_penalty_func
    use sparse_embed, only: SPAM3_EMBED

    ! Arguments
    real(kind=DP), intent(out) :: value
    integer, intent(in) :: num_states
    type(SPAM3_EMBED), intent(in) :: vecs(num_states,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: vecs_cov(num_states,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: val_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: Svkv(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: kcSc(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_overlap
    type(SPAM3_EMBED), intent(in) :: val_overlap
    integer, intent(in) :: num_conv_states

    ! local variables
    real(kind=DP) :: temp_val, temp_val2
    integer :: icount

    temp_val=0.0_DP
    temp_val2=0.0_DP

    if(pub_lr_tddft_penalty_func) then

       do icount=1+num_conv_states, num_states
          call lr_tddft_utils_penalty_func(temp_val2,vecs(icount,:),&
             vecs_cov(icount,:),cond_denskern,val_denskern,&
             cond_overlap,val_overlap,Svkv,kcSc)

          temp_val=temp_val+temp_val2
       enddo

       value=temp_val

    else
       value=0.0_DP
    endif

  end subroutine lr_tddft_utils_penalty_batch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_utils_penalty_batch_RPA(value,vecs,vecs_cov,cond_denskern, &
      val_denskern, cond_overlap,val_overlap,Svkv,kcSc,num_conv_states,num_states)
    !========================================================================!
    ! Subroutine computes total penalty value for a batch of vectors in RPA. !
    ! Modified for embedding by Joseph Prentice, July 2018                   !
    !========================================================================!

    use rundat, only: pub_num_spins, pub_lr_tddft_penalty_func
    use sparse_embed, only: SPAM3_EMBED

    ! Arguments
    real(kind=DP), intent(out) :: value
    integer, intent(in) :: num_states
    type(SPAM3_EMBED), intent(in) :: vecs(num_states,2,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: vecs_cov(num_states,2,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: val_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: Svkv(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: kcSc(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_overlap
    type(SPAM3_EMBED), intent(in) :: val_overlap
    integer, intent(in) :: num_conv_states

    ! local variables
    real(kind=DP) :: temp_val, temp_val2, temp_val3
    integer :: icount

    temp_val=0.0_DP
    temp_val2=0.0_DP
    temp_val3=0.0_DP

    if(pub_lr_tddft_penalty_func) then

       do icount=1+num_conv_states, num_states
          call lr_tddft_utils_penalty_func(temp_val2,vecs(icount,1,:),&
               vecs_cov(icount,1,:),cond_denskern,val_denskern,&
               cond_overlap,val_overlap,Svkv,kcSc)

          call lr_tddft_utils_penalty_func(temp_val3,vecs(icount,2,:),&
               vecs_cov(icount,2,:),cond_denskern,val_denskern,&
               cond_overlap,val_overlap,Svkv,kcSc)

          temp_val=temp_val+(temp_val2+temp_val3)/2.0_DP
       enddo

       value=temp_val

    else

       value=0.0_DP

    endif

  end subroutine lr_tddft_utils_penalty_batch_RPA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_utils_penalty_func(value1, vec, vec_cov,cond_denskern, &
       val_denskern, cond_overlap, val_overlap,Svkv,kcSc)
    !========================================================================!
    ! Subroutine is calculating the value of the penalty norm Q[K]. For a    !
    ! fully dense response kernel, this is 0 by construction. However,       !
    ! as soon as a cutoff is introduced, Q[K] is a measure of how much       !
    ! the kernel deviates from being a valid response kernel.                !
    ! Modified for embedding by Joseph Prentice, July 2018                   !
    !========================================================================!

    use rundat, only: pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_transpose, &
         sparse_embed_transpose_structure, sparse_embed_trace

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: value1
    type(SPAM3_EMBED), intent(in) :: vec(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: vec_cov(pub_num_spins) ! covariant version of vec
    type(SPAM3_EMBED), intent(in) :: cond_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: val_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_overlap
    type(SPAM3_EMBED), intent(in) :: val_overlap
    type(SPAM3_EMBED), intent(in) :: Svkv(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: kcSc(pub_num_spins)

    ! Local Variables
    type(SPAM3_EMBED) :: kSvkv
    type(SPAM3_EMBED) :: k_primed, vec_trans
    type(SPAM3_EMBED) :: kSv, SckSv, kSckSv, temp
    real(kind=DP) :: trace
    integer :: is

    ! allocate data
    call sparse_embed_create(kSvkv, vec(1), Svkv(1))
    ! change sparsity pattern of kprimed to that of k
    call sparse_embed_create(k_primed, vec(1))
    call sparse_embed_transpose_structure(vec_trans%structure, vec(1))
    call sparse_embed_create(vec_trans, iscmplx=vec(1)%p%iscmplx)
    call sparse_embed_create(kSv, vec(1), val_overlap)
    call sparse_embed_create(SckSv, cond_overlap, kSv)
    call sparse_embed_create(kSckSv, vec_trans, SckSv)
    call sparse_embed_create(temp, kSckSv)

    value1=0.0_DP

    ! sum over spin channels:
    do is=1, pub_num_spins
       ! calculate k_primed
       call sparse_embed_product(kSvkv, vec(is), Svkv(is))
       call sparse_embed_product(k_primed, kcSc(is), kSvkv)

       ! evaluate the functional
       call sparse_embed_transpose(vec_trans, vec(is))
       call sparse_embed_product(kSckSv, vec_trans, vec_cov(is))

       call sparse_embed_product(kSv, k_primed, val_overlap)
       call sparse_embed_product(SckSv, cond_overlap, kSv)
       call sparse_embed_transpose(vec_trans, k_primed)
       call sparse_embed_product(temp, vec_trans, SckSv)

       call sparse_embed_axpy(kSckSv, temp, -1.0_DP)

       ! now square the penalty functional and find the trace
       call sparse_embed_trace(trace, kSckSv, kSckSv)
       value1 =value1+trace
    enddo

    ! deallocate data
    call sparse_embed_destroy(kSvkv)
    call sparse_embed_destroy(k_primed)
    call sparse_embed_destroy(vec_trans)
    call sparse_embed_destroy(kSv)
    call sparse_embed_destroy(SckSv)
    call sparse_embed_destroy(kSckSv)
    call sparse_embed_destroy(temp)

  end subroutine lr_tddft_utils_penalty_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_utils_init_tddft_vecs(vecs,cond_ngwf_basis,&
    val_ngwf_basis,cond_rep,val_rep,cond_ham,val_ham,&
    val_denskern,val_ngwf_num,cond_ngwf_num,mdl,num_states)

    !======================================================================!
    ! Subroutine computes the lowest N conduction and highest N valence    !
    ! states using iterative processes that do not require any direct      !
    ! diagonalisation. It then initialises the N TDDFT kernels to kernels  !
    ! constructed from the NxN set of possible vectors. Either the         !
    ! initialisation is done such that those transitions with the smallest !
    ! Kohn-Sham eigenvalue differences are picked, or the initialisation   !
    ! is performed such that the highest overlap between electron and      !
    ! hole is achieved.                                                    !
    ! Modified for embedding by Joseph Prentice, July 2018                 !
    !======================================================================!

    use datatypes, only: FUNCTIONS, data_functions_alloc, data_functions_dealloc
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_product_on_grid
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_product, &
         sparse_embed_create, sparse_embed_destroy, sparse_embed_axpy, &
         sparse_embed_copy, sparse_embed_transpose_structure, &
         sparse_embed_transpose, sparse_embed_extremal_eigenvalue, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_destroy
    use rundat, only: pub_num_spins, pub_lr_tddft_init_max_overlap
    use utils, only: utils_alloc_check, utils_dealloc_check

    ! Arguments
    integer, intent(in) :: num_states
    type(SPAM3_EMBED), intent(inout) :: vecs(num_states,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_ham(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: val_ham(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: val_denskern(pub_num_spins)
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis(:)
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type(NGWF_REP), intent(in) :: cond_rep
    type(NGWF_REP), intent(in) :: val_rep
    type(MODEL), intent(in) :: mdl
    integer, intent(in) :: val_ngwf_num, cond_ngwf_num

    ! local variables
    integer :: is
    integer :: loc_num_states
    integer :: icount,jcount,counter,ierr
    integer :: index_pair(2)
    real(kind=DP) :: shift,temp_val
    type(FUNCTIONS), dimension(:,:), allocatable :: val_evecs, cond_evecs
    real(kind=DP), allocatable, dimension (:,:) :: cond_eval
    real(kind=DP), allocatable, dimension (:,:) :: val_eval
    real(kind=DP), allocatable, dimension(:,:,:) :: energy_list ! list of
      ! effective Kohn-Sham excitation energies
    type(SPAM3_EMBED) :: temp_proj_ham, HPS,PS, SPS, SP,proj, SinvH,temp_dkn
    type(SPAM3_EMBED) :: SvalcondSinv,trans,PSvalcondSinv,SinvScondval
    logical :: maximise_overlap ! will be a user defined variable in future
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: cond_dens_batch
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: val_dens_batch
    type(SPAM3_EMBED), allocatable, dimension(:,:) :: temp_cond_dkn,temp_val_dkn
    logical :: loc_cmplx

    ! jcap: embedding local variables
    type(SPAM3) :: val_kern_array(pub_num_spins)
    type(SPAM3) :: cond_kern_array(pub_num_spins)
    integer :: isub,jsub
    real (kind=DP), allocatable :: sub_val_dens_batch(:,:,:,:)
    real (kind=DP), allocatable :: sub_cond_dens_batch(:,:,:,:)


    ! jmecmplx
    loc_cmplx = vecs(1,1)%p%iscmplx

    maximise_overlap=pub_lr_tddft_init_max_overlap
    loc_num_states=num_states
    if(loc_num_states>minval(val_rep%n_occ(:,1))) then
      loc_num_states=minval(val_rep%n_occ(:,1))
    endif

    ! allocate and create data structures
    allocate(val_evecs(loc_num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_init_tddft_vecs','val_evecs',ierr)
    allocate(cond_evecs(loc_num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_init_tddft_vecs','cond_evecs',ierr)
    do is = 1, pub_num_spins
       do icount = 1, loc_num_states
          call data_functions_alloc(val_evecs(icount,is), &
               val_ngwf_num, iscmplx=loc_cmplx)
          call data_functions_alloc(cond_evecs(icount,is), &
               cond_ngwf_num, iscmplx=loc_cmplx)
       end do
    end do
    allocate(energy_list(loc_num_states,loc_num_states,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_utils_init_tddft_vecs','energy_list',ierr)
    call sparse_embed_create(temp_proj_ham,val_ham(1))
    call sparse_embed_create(proj,val_denskern(1))
    call sparse_embed_create(SP,val_rep%overlap,proj)
    call sparse_embed_create(SinvH,val_rep%inv_overlap,val_ham(1))
    call sparse_embed_create(PS,proj,val_rep%overlap)
    call sparse_embed_create(HPS,val_ham(1),PS)
    call sparse_embed_create(SPS,val_ham(1))
    call sparse_embed_create(temp_dkn,val_denskern(1))
    allocate(val_eval(loc_num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_init_tddft_vecs','val_eval',ierr)
    allocate(cond_eval(loc_num_states,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_init_tddft_vecs','cond_eval',ierr)

    ! Start with computing the num_states highest occupied states of the valence
    ! hamiltonian.
    ! first compute the smallest eval of the valence hamiltonian to determine
    ! shift value.
    ! do so for each spin channel
    do is=1, pub_num_spins
       call sparse_embed_product(SinvH,val_rep%inv_overlap,val_ham(is))
       call sparse_embed_extremal_eigenvalue(SinvH,val_rep%overlap,shift,&
            tol=0.00001_DP,min_val=.true.)
       shift=shift-0.1_DP

       ! create a projected hamiltonian, with all valence states projected to the
       ! shift value
       ! create an effective projector onto all conduction states
       call sparse_embed_copy(temp_dkn,val_rep%inv_overlap)
       call sparse_embed_axpy(temp_dkn,val_denskern(is),-1.0_DP)
       ! use this as projector to project out all conduction states from ham.
       call sparse_embed_product(PS,temp_dkn,val_rep%overlap)
       call sparse_embed_product(SP,val_rep%overlap,temp_dkn)
       call sparse_embed_product(SPS, SP, val_rep%overlap)
       call sparse_embed_product(HPS,val_ham(is),PS)
       call sparse_embed_product(proj,SP,HPS)
       call sparse_embed_axpy(proj,SPS,-1.0_DP*shift)
       call sparse_embed_copy(temp_proj_ham,val_ham(is))
       call sparse_embed_axpy(temp_proj_ham,proj,-1.0_DP)
       call sparse_embed_product(SinvH,val_rep%inv_overlap,temp_proj_ham)

       ! now compute the num_states highest eigenvalues of proj hamiltonian
       ! these correspond to the highest occupied Kohn-Sham states
       do icount=1, loc_num_states
          ! compute highest eval of projected ham
          call sparse_embed_extremal_eigenvalue(SinvH,val_rep%overlap, &
               val_eval(icount,is), tol=0.00001_DP,evec=val_evecs(icount,is))

          ! project out the newly converged state from Hamiltonian
          call internal_sparse_outer_product(proj, val_evecs(icount,is), &
               val_evecs(icount,is))
          temp_val=-1.0_DP*(val_eval(icount,is)-shift)
          call sparse_embed_product(PS,proj,val_rep%overlap)
          call sparse_embed_axpy(SinvH,PS,temp_val)
       enddo
    enddo ! end of spin loop

    ! we have successfuly computed num_states valence states. Now do the
    ! same for conduction hamiltonian and conduction states.
    call sparse_embed_destroy(temp_proj_ham)
    call sparse_embed_destroy(proj)
    call sparse_embed_destroy(SP)
    call sparse_embed_destroy(SinvH)
    call sparse_embed_destroy(PS)
    call sparse_embed_destroy(HPS)
    call sparse_embed_destroy(SPS)
    call sparse_embed_destroy(temp_dkn)
    ! reallocate data structures for cond ham
    call sparse_embed_create(temp_proj_ham,cond_ham(1))
    call sparse_embed_create(proj,cond_rep%inv_overlap)
    call sparse_embed_create(SP,cond_rep%overlap,proj)
    call sparse_embed_create(SinvH,cond_rep%inv_overlap,cond_ham(1))
    call sparse_embed_create(PS,proj,cond_rep%overlap)
    call sparse_embed_create(HPS,cond_ham(1),PS)
    call sparse_embed_create(SPS,cond_ham(1))
    call sparse_embed_create(temp_dkn,cond_rep%inv_overlap)
    ! Data structures for constructing temp_dkn=SinvSchiphiPSphichiSinv
    call sparse_embed_create(SvalcondSinv,cond_rep%cross_overlap,&
         cond_rep%inv_overlap)
    call sparse_embed_transpose_structure(trans%structure,&
         cond_rep%cross_overlap)
    call sparse_embed_create(trans)
    call sparse_embed_create(SinvScondval,cond_rep%inv_overlap,trans)
    call sparse_embed_create(PSvalcondSinv,val_denskern(1),SvalcondSinv)

    ! loop over spin indices
    do is=1, pub_num_spins
       ! compute maximum eval of cond ham for projection
       call sparse_embed_product(SinvH,cond_rep%inv_overlap,cond_ham(is))
       call sparse_embed_extremal_eigenvalue(SinvH,cond_rep%overlap,shift,&
            tol=0.00001_DP)
       shift=shift+0.1_DP

       ! create an effective projected conduction hamiltonian, with all
       ! valence states projected above the shift
       call sparse_embed_product(SvalcondSinv,cond_rep%cross_overlap,cond_rep%inv_overlap)
       call sparse_embed_transpose(trans,cond_rep%cross_overlap)
       call sparse_embed_product(SinvScondval,cond_rep%inv_overlap,trans)
       call sparse_embed_product(PSvalcondSinv,val_denskern(is),SvalcondSinv)
       call sparse_embed_product(temp_dkn,SinvScondval,PSvalcondSinv)
       call sparse_embed_product(PS,temp_dkn,cond_rep%overlap)
       call sparse_embed_product(SP,cond_rep%overlap,temp_dkn)
       call sparse_embed_product(SPS, SP, cond_rep%overlap)
       call sparse_embed_product(HPS,cond_ham(is),PS)
       call sparse_embed_product(proj,SP,HPS)
       call sparse_embed_axpy(proj,SPS,-1.0_DP*shift)
       call sparse_embed_copy(temp_proj_ham,cond_ham(is))
       call sparse_embed_axpy(temp_proj_ham,proj,-1.0_DP)
       call sparse_embed_product(SinvH,cond_rep%inv_overlap,temp_proj_ham)

       ! now compute the num_states lowest eigenvalues of proj hamiltonian
       ! these correspond to the lowest unoccupied Kohn-Sham states
       do icount=1, loc_num_states
          ! compute highest eval of projected ham
          call sparse_embed_extremal_eigenvalue(SinvH, cond_rep%overlap, &
               cond_eval(icount,is), tol=0.00001_DP, min_val=.true., &
               evec=cond_evecs(icount,is))

          ! project out the newly converged state from Hamiltonian
          call internal_sparse_outer_product(proj, cond_evecs(icount,is),&
             cond_evecs(icount,is))
          temp_val=-1.0_DP*(cond_eval(icount,is)-shift)
          call sparse_embed_product(PS,proj,cond_rep%overlap)
          call sparse_embed_axpy(SinvH,PS,temp_val)
       enddo
    enddo ! end loop over spin indices

    ! successfully computed the lowest num_states conduction and highest
    ! num states valence Kohn-Sham states in linear scaling effort.
    ! now fill initial states: If maximise overlaps=true we aim to initialise
    ! states to those that have the maximum overlap between effective electron
    ! and hole densities. If not, we aim to have the minimum Kohn sham gaps
    ! if this is a spin polarised calculation, both spin channels are filled
    ! by only considering overlaps or min KS gaps in their spin channel. This
    ! might not be the ideal way of doing it but will undoubtedly work.
    if(.not. maximise_overlap) then
       ! compute Kohn-Sham energies
       do icount=1, loc_num_states
          do jcount=1, loc_num_states
             do is=1, pub_num_spins
                energy_list(icount,jcount,is)=cond_eval(icount,is)-val_eval(jcount,is)
             enddo
          enddo
       enddo
       ! find the smallest num_states eigenvalues and initialise the batch
       ! of vectors to Kohn-Sham transitions corresponding to them.
       ! do so for both spin channels
       do is=1, pub_num_spins
          do counter=1, num_states
             temp_val=1000.0_DP
             do icount=1,loc_num_states
                do jcount=1,loc_num_states
                   if (energy_list(icount,jcount,is)<temp_val) then
                      temp_val=energy_list(icount,jcount,is)
                      index_pair(1)=icount
                      index_pair(2)=jcount
                   endif
                enddo
             enddo
             ! found the index pair corresponding to the smallest KS transition
             ! now set TDDFT vector
             call internal_sparse_outer_product(vecs(counter,is), &
                  cond_evecs(index_pair(1),is), val_evecs(index_pair(2),is))
             ! make sure each state is only found once
             energy_list(index_pair(1),index_pair(2),is)=1000_DP
          enddo
       enddo ! end spin
    else
       ! instead of minimizing the KS gap of the inital TDDFT vectors, maximise
       ! the overlap of effective electron and hole wavefunctions.
       allocate(temp_cond_dkn(pub_num_spins,1), stat=ierr)
       call utils_alloc_check('lr_tddft_utils_init_tddft_vecs','temp_cond_dkn',ierr)
       allocate(temp_val_dkn(pub_num_spins,1), stat=ierr)
       call utils_alloc_check('lr_tddft_utils_init_tddft_vecs','temp_val_dkn',ierr)
       do is=1, pub_num_spins
          call sparse_embed_create(temp_val_dkn(is,1),val_rep%inv_overlap)
          call sparse_embed_create(temp_cond_dkn(is,1),cond_rep%inv_overlap)
       enddo
       allocate(val_dens_batch(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
          mdl%fine_grid%max_slabs12,pub_num_spins,loc_num_states),stat=ierr)
       call utils_alloc_check('lr_tddft_utils_init_tddft_vecs',&
            'val_dens_batch',ierr)
       allocate(cond_dens_batch(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
            mdl%fine_grid%max_slabs12,pub_num_spins,loc_num_states),stat=ierr)
       call utils_alloc_check('lr_tddft_utils_init_tddft_vecs',&
            'cond_dens_batch',ierr)

       allocate(sub_val_dens_batch(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
            mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('lr_tddft_utils_init_tddft_vecs',&
            'sub_val_dens_batch',ierr)
       allocate(sub_cond_dens_batch(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
            mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('lr_tddft_utils_init_tddft_vecs',&
            'sub_cond_dens_batch',ierr)

       ! construct densities associated with all electrons and holes
       ! built from the pure KS states converged
       do icount=1, loc_num_states
          val_dens_batch(:,:,:,:,icount)=0.d0
          cond_dens_batch(:,:,:,:,icount)=0.d0
          do is=1, pub_num_spins
             call internal_sparse_outer_product(temp_val_dkn(is,1), &
                  val_evecs(icount,is), val_evecs(icount,is))
             call internal_sparse_outer_product(temp_cond_dkn(is,1), &
                  cond_evecs(icount,is), cond_evecs(icount,is))
          enddo

          ! jcap: Loop over regions and sum up densities
          do isub=1,mdl%nsub
             do jsub=1,mdl%nsub
                sub_val_dens_batch=0.d0
                sub_cond_dens_batch=0.d0

                ! rc2013: EMBED_FIX! problems with deferred shape arrays using
                ! embedding structures.  In practice this will probably
                ! require denskern to be added properly to density_on_grid
                ! (complex quantities).
                call sparse_embed_extract_from_array(val_kern_array,&
                     temp_val_dkn(:,1),isub,jsub)
                call sparse_embed_extract_from_array(cond_kern_array,&
                     temp_cond_dkn(:,1),isub,jsub)

                call density_on_grid(sub_cond_dens_batch,&
                     mdl%fine_grid,mdl%dbl_grid,mdl%cell,mdl%fftbox,&
                     cond_kern_array,cond_rep%ngwf_overlap%m(isub,jsub),&
                     cond_rep%ngwfs_on_grid(isub), cond_ngwf_basis(isub),&
                     cond_rep%ngwfs_on_grid(jsub), cond_ngwf_basis(jsub))

                call density_on_grid(sub_val_dens_batch,&
                     mdl%fine_grid,mdl%dbl_grid,mdl%cell,mdl%fftbox,&
                     val_kern_array,val_rep%ngwf_overlap%m(isub,jsub),&
                     val_rep%ngwfs_on_grid(isub), val_ngwf_basis(isub),&
                     val_rep%ngwfs_on_grid(jsub), val_ngwf_basis(jsub))

                cond_dens_batch(:,:,:,:,icount)=cond_dens_batch(:,:,:,:,icount)+&
                     sub_cond_dens_batch
                val_dens_batch(:,:,:,:,icount)=val_dens_batch(:,:,:,:,icount)+&
                     sub_val_dens_batch

                call sparse_embed_destroy_extracted_array(val_kern_array)
                call sparse_embed_destroy_extracted_array(cond_kern_array)

             end do
          end do

       enddo

       deallocate(sub_val_dens_batch,stat=ierr)
       call utils_dealloc_check('lr_tddft_utils_init_tddft_vecs',&
            'sub_val_dens_batch',ierr)
       deallocate(sub_cond_dens_batch,stat=ierr)
       call utils_dealloc_check('lr_tddft_utils_init_tddft_vecs',&
            'sub_cond_dens_batch',ierr)


       ! fill effective matrix of overlap integrals
       do icount=1, loc_num_states
          do jcount=1, loc_num_states
             do is=1,pub_num_spins
                energy_list(icount,jcount,is)=integrals_product_on_grid(mdl%fine_grid,&
                cond_dens_batch(:,:,:,is,icount),val_dens_batch(:,:,:,is,jcount))
             enddo
          enddo
       enddo

       ! now fill states according to their highest overlaps
       ! loop over spins
       do is=1, pub_num_spins
          if(num_states>loc_num_states) then
             do counter=1,num_states
                temp_val=0.0_DP
                do icount=1,loc_num_states
                   do jcount=1,loc_num_states
                      if (energy_list(icount,jcount,is)>temp_val) then
                         temp_val=energy_list(icount,jcount,is)
                         index_pair(1)=icount
                         index_pair(2)=jcount
                      endif
                   enddo
                enddo
                ! found the index pair corresponding to the smallest KS transition
                ! now set TDDFT vector
                call internal_sparse_outer_product(vecs(counter,is), &
                     cond_evecs(index_pair(1),is), val_evecs(index_pair(2),is))
                ! make sure each state is only found once
                energy_list(index_pair(1),index_pair(2),is)=0.0_DP
            enddo
          else
             do icount=1,num_states
                temp_val=0.0_DP
                index_pair(1)=icount
                do jcount=1,loc_num_states
                   if (energy_list(icount,jcount,is)> temp_val) then
                      temp_val=energy_list(icount,jcount,is)
                      index_pair(2)=jcount
                   endif
                enddo
                call internal_sparse_outer_product(vecs(icount,is), &
                     cond_evecs(index_pair(1),is), val_evecs(index_pair(2),is))
             enddo
          endif
       enddo ! end of spin loop

       ! deallocate data structure
       do is=1, pub_num_spins
          call sparse_embed_destroy(temp_val_dkn(is,1))
          call sparse_embed_destroy(temp_cond_dkn(is,1))
       enddo
       deallocate(temp_cond_dkn,stat=ierr)
       call utils_dealloc_check('lr_tddft_utils_init_tddft_vecs',&
          'temp_cond_dkn',ierr)
       deallocate(temp_val_dkn,stat=ierr)
       call utils_dealloc_check('lr_tddft_utils_init_tddft_vecs',&
          'temp_val_dkn',ierr)
       deallocate(val_dens_batch,stat=ierr)
       call utils_dealloc_check('lr_tddft_utils_init_tddft_vecs',&
          'val_dens_batch',ierr)
       deallocate(cond_dens_batch,stat=ierr)
       call utils_dealloc_check('lr_tddft_utils_init_tddft_vecs',&
          'cond_dens_batch',ierr)
    endif

    ! deallocate data structures
    call sparse_embed_destroy(temp_proj_ham)
    call sparse_embed_destroy(proj)
    call sparse_embed_destroy(SP)
    call sparse_embed_destroy(SinvH)
    call sparse_embed_destroy(PS)
    call sparse_embed_destroy(HPS)
    call sparse_embed_destroy(SPS)
    call sparse_embed_destroy(temp_dkn)
    call sparse_embed_destroy(SvalcondSinv)
    call sparse_embed_destroy(trans)
    call sparse_embed_destroy(SinvScondval)
    call sparse_embed_destroy(PSvalcondSinv)
    do is = pub_num_spins, 1, -1
       do icount = loc_num_states, 1, -1
          call data_functions_dealloc(cond_evecs(icount,is))
          call data_functions_dealloc(val_evecs(icount,is))
       end do
    end do
    deallocate(cond_evecs,stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_init_tddft_vecs','cond_evecs',ierr)
    deallocate(val_evecs,stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_init_tddft_vecs','val_evecs',ierr)
    deallocate(energy_list,stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_init_tddft_vecs','energy_list',ierr)
    deallocate(val_eval, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_init_tddft_vecs','val_eval',ierr)
    deallocate(cond_eval, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_init_tddft_vecs','cond_eval',ierr)

  contains

     subroutine internal_sparse_outer_product(mat, uvec, vvec)

        !=====================================================================!
        ! Wrapper for sparse_outer_product, mat_ji = u_j * conjg(v_i)         !
        ! Written by JM Escartin, August 2016.                                !
        !=====================================================================!

        use datatypes, only: FUNCTIONS
        use sparse_embed, only: SPAM3_EMBED, sparse_embed_outer_product
        implicit none

        ! Arguments
        type(SPAM3_EMBED), intent(inout) :: mat
        type(FUNCTIONS), intent(in) :: uvec ! row vector u
        type(FUNCTIONS), intent(in) :: vvec ! column vector v

        if (loc_cmplx) then ! loc_cmplx from the containing subroutine
           call sparse_embed_outer_product(mat, uvec%z, vvec%z)
        else
           call sparse_embed_outer_product(mat, uvec%d, vvec%d)
        end if

     end subroutine internal_sparse_outer_product

  end subroutine lr_tddft_utils_init_tddft_vecs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_utils_eh_denskern_RPA(eh_kernel,vec,cond_overlap,&
     val_overlap,cond_inv_overlap,val_inv_overlap,val_cond_overlap)

    !=====================================================================!
    ! This subroutine constructs the effective electron-hole difference   !
    ! density kernel in the RPA. The routine assumes that the conduction  !
    ! space manifold is represented by the joint rep.                     !
    ! Modified for embedding by Joseph Prentice, July 2018                !
    !=====================================================================!

    use sparse_embed, only: SPAM3_EMBED, sparse_embed_product, &
         sparse_embed_create, sparse_embed_destroy, sparse_embed_axpy, &
         sparse_embed_scale, sparse_embed_transpose_structure, sparse_embed_transpose
    use rundat, only: pub_num_spins

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: eh_kernel(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: vec(2,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_overlap
    type(SPAM3_EMBED), intent(in) :: val_overlap
    type(SPAM3_EMBED), intent(in) :: cond_inv_overlap
    type(SPAM3_EMBED), intent(in) :: val_inv_overlap
    type(SPAM3_EMBED), intent(in) :: val_cond_overlap

    ! local variables
    type(SPAM3_EMBED) :: e_kernel, h_kernel, h_kernel_temp, PSval,ScondP,&
                 SvcSinv,PvSvcSinv, trans, SinvScv
    integer :: is

    ! allocate data structures
    call sparse_embed_create(e_kernel,cond_inv_overlap)
    call sparse_embed_create(h_kernel,val_inv_overlap)
    call sparse_embed_create(h_kernel_temp,val_inv_overlap)
    call sparse_embed_create(PSval,vec(1,1),val_overlap)
    call sparse_embed_create(ScondP, cond_overlap,vec(1,1))
    call sparse_embed_create(SvcSinv,val_cond_overlap,cond_inv_overlap)
    call sparse_embed_create(PvSvcSinv,h_kernel,SvcSinv)
    call sparse_embed_transpose_structure(trans%structure,vec(1,1))
    call sparse_embed_create(trans)

    do is=1,pub_num_spins
       ! construct electron part first
       call sparse_embed_product(PSval,vec(1,is), val_overlap)
       call sparse_embed_transpose(trans, vec(1,is))
       call sparse_embed_product(eh_kernel(is), PSval,trans)
       call sparse_embed_transpose(trans,vec(2,is))
       call sparse_embed_product(PSval,vec(2,is),val_overlap)
       call sparse_embed_product(e_kernel,PSval,trans)
       call sparse_embed_axpy(eh_kernel(is),e_kernel,1.0_DP)
       call sparse_embed_scale(eh_kernel(is),0.5_DP)

       ! now construct hole part
       call sparse_embed_product(ScondP,cond_overlap,vec(1,is))
       call sparse_embed_transpose(trans, vec(1,is))
       call sparse_embed_product(h_kernel, trans,ScondP)
       call sparse_embed_product(ScondP,cond_overlap,vec(2,is))
       call sparse_embed_transpose(trans,vec(2,is))
       call sparse_embed_product(h_kernel_temp, trans,ScondP)
       call sparse_embed_axpy(h_kernel,h_kernel_temp,1.0_DP)
       call sparse_embed_scale(h_kernel,0.5_DP)
    enddo

    call sparse_embed_destroy(trans)
    call sparse_embed_transpose_structure(trans%structure,val_cond_overlap)
    call sparse_embed_create(trans)
    call sparse_embed_create(SinvScv,cond_inv_overlap,trans)

    do is=1, pub_num_spins
       ! project hole part into conduction space and subtract
       ! from electron part to get electron-hole difference
       call sparse_embed_product(SvcSinv,val_cond_overlap,cond_inv_overlap)
       call sparse_embed_product(PvSvcSinv,h_kernel,SvcSinv)
       call sparse_embed_transpose(trans, val_cond_overlap)
       call sparse_embed_product(SinvScv,cond_inv_overlap,trans)
       call sparse_embed_product(e_kernel,SinvScv,PvSvcSinv)
       call sparse_embed_axpy(eh_kernel(is),e_kernel,-1.0_DP)
    enddo

    ! deallocate data structures
    call sparse_embed_destroy(SinvScv)
    call sparse_embed_destroy(e_kernel)
    call sparse_embed_destroy(h_kernel)
    call sparse_embed_destroy(h_kernel_temp)
    call sparse_embed_destroy(PSval)
    call sparse_embed_destroy(ScondP)
    call sparse_embed_destroy(SvcSinv)
    call sparse_embed_destroy(PvSvcSinv)
    call sparse_embed_destroy(trans)

  end subroutine lr_tddft_utils_eh_denskern_RPA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_utils_precond_batch_RPA(g_batch,f_batch,num_states,&
       num_conv_states,response_dens,sub_response_dens,fxc_fine,&
       cond_rep,cond_ngwf_basis,cond_denskern,kchc,val_rep,&
       val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,ground_state_dens,&
       sub_ground_state_dens)

    !====================================================================!
    ! Subroutine applying the preconditioner to a batch of f_vecs, to ob-!
    ! tain a batch of preconditioned g_vecs. Only vectors that are not   !
    ! converged yet get the preconditioning applied to them. Routine     !
    ! also prints out the final residue of vectors obtained.             !
    ! Modified for embedding by Joseph Prentice, July 2018               !
    !====================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, NORMAL
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_num_spins, pub_output_detail
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    integer, intent(in) :: num_states
    integer, intent(in) :: num_conv_states
    type(SPAM3_EMBED), intent(inout) :: f_batch(num_states,2,pub_num_spins) ! gradient of the system
      ! with respect to a solution vector x_vec.
    type(SPAM3_EMBED), intent(inout) :: g_batch(num_states,2,pub_num_spins) ! preconditioned version of f
    type(MODEL), intent(in) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
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
    integer :: icount,is,ierr
    real(kind=DP) :: final_tol
    type(SPAM3_EMBED) :: PSv
    type(SPAM3_EMBED), allocatable, dimension(:,:) :: eff_fvec

    allocate(eff_fvec(2,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_utils_precond_batch_RPA',&
      'eff_fvec',ierr)

    call sparse_embed_create(PSv, f_batch(1,1,1), val_rep%overlap)
    do is=1,pub_num_spins
       call sparse_embed_create(eff_fvec(1,is), cond_rep%overlap, PSv)
       call sparse_embed_create(eff_fvec(2,is), cond_rep%overlap, PSv)
    enddo

    ! first output some printing information
    ! output preconditioner final residue
    if (pub_on_root .and. pub_output_detail>=NORMAL) then
       write(stdout,'(a)') '******  Iteratively applying preconditioner&
        & on G^{1}: ******'
       write(stdout, '(a)') '    |STATE|    Final residue  &
         &|'
    endif

    do icount=1+num_conv_states, num_states
       do is=1,pub_num_spins
          call sparse_embed_product(PSv,f_batch(icount,1,is), val_rep%overlap)
          call sparse_embed_product(eff_fvec(1,is),cond_rep%overlap, PSv)
          call sparse_embed_product(PSv,f_batch(icount,2,is), val_rep%overlap)
          call sparse_embed_product(eff_fvec(2,is),cond_rep%overlap, PSv)
       enddo

       call lr_tddft_utils_precond_RPA(g_batch(icount,:,:),eff_fvec,&
       response_dens,sub_response_dens,fxc_fine,&
       cond_rep,cond_ngwf_basis,cond_denskern,kchc,val_rep,&
       val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,ground_state_dens,&
       sub_ground_state_dens,final_tol)

       if(pub_on_root .and. pub_output_detail>=NORMAL) &
          write(stdout,'(i7,f16.10)') icount, final_tol
    enddo

    ! deallocate temporary data structures
    call sparse_embed_destroy(PSv)
    do is=1, pub_num_spins
       call sparse_embed_destroy(eff_fvec(1,is))
       call sparse_embed_destroy(eff_fvec(2,is))
    enddo
    deallocate(eff_fvec, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_precond_batch_RPA',&
      'eff_fvec',ierr)

  end subroutine lr_tddft_utils_precond_batch_RPA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_utils_precond_batch(g_batch,f_batch,num_states,&
       num_conv_states,response_dens,sub_response_dens,fxc_fine,&
       cond_rep,cond_ngwf_basis,cond_denskern,kchc,val_rep,&
       val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,ground_state_dens,&
       sub_ground_state_dens)

    !====================================================================!
    ! Subroutine applying the preconditioner to a batch of f_vecs, to ob-!
    ! tain a batch of preconditioned g_vecs. Only vectors that are not   !
    ! converged yet get the preconditioning applied to them. Routine     !
    ! also prints out the final residue of vectors obtained.             !
    ! Modified for embedding by Joseph Prentice, July 2018               !
    !====================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, NORMAL
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_num_spins, pub_output_detail
    use sparse_embed, only: SPAM3_EMBED

    implicit none

    integer, intent(in) :: num_states
    integer, intent(in) :: num_conv_states
    type(SPAM3_EMBED), intent(inout) :: f_batch(:,:) ! gradient of the full system
      ! with respect to a solution vector x_vec. ! jd: Dims: num_states, pub_num_spins. Leave ':,:' be.
    type(SPAM3_EMBED), intent(inout) :: g_batch(:,:) ! preconditioned version of f
                                                 ! jd: Dims: num_states, pub_num_spins. Leave ':,:' be.
    type(MODEL), intent(in) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
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
    integer :: icount
    real(kind=DP) :: final_tol

    ! first output some printing information
    ! output preconditioner final residue
    if (pub_on_root .and. pub_output_detail>=NORMAL) then
       write(stdout,'(a)') '******  Iteratively applying preconditioner&
        & on G^{1}: ******'
       write(stdout, '(a)') '    |STATE|    Final residue  &
         &|'
    endif

    do icount=1+num_conv_states, num_states
       call lr_tddft_utils_precond(g_batch(icount,:),f_batch(icount,:),&
            response_dens,sub_response_dens,fxc_fine,&
            cond_rep,cond_ngwf_basis,cond_denskern,kchc,val_rep,&
            val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,ground_state_dens,&
            sub_ground_state_dens,final_tol)

       if(pub_on_root .and. pub_output_detail>=NORMAL) &
        write(stdout,'(i7,f16.10)') icount,final_tol
    enddo

  end subroutine lr_tddft_utils_precond_batch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_utils_precond_RPA(g_vec,f_vec,response_dens,&
       sub_response_dens,fxc_fine,&
       cond_rep,cond_ngwf_basis,cond_denskern,kchc,val_rep,&
       val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,ground_state_dens,&
       sub_ground_state_dens,final_tol)

    !====================================================================!
    ! Subroutine applies the diagonal preconditioner of Kohn sham eigen- !
    ! value differences onto a single vector f_vec. Since the precon-    !
    ! ditioner is only diagonal in KS-eigenstate space, its action has   !
    ! to be computed by solving a linear system using a conjugate        !
    ! gradient algorithm.                                                !
    ! Modified for embedding by Joseph Prentice, July 2018               !
    !====================================================================!

    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use linear_response, only: linear_response_RPA_operator
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_lr_tddft_preopt, pub_lr_tddft_precond_iter, &
         pub_lr_tddft_precond_tol, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_transpose, sparse_embed_trace, &
         sparse_embed_transpose_structure, sparse_embed_copy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    type(SPAM3_EMBED), intent(inout) :: f_vec(2,pub_num_spins) ! gradient of the full system
      ! with respect to a solution vector x_vec.
    type(SPAM3_EMBED), intent(inout) :: g_vec(2,pub_num_spins) ! preconditioned version of f
    type(MODEL), intent(in) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
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
    real(kind=DP), intent(inout) :: final_tol

    ! temporary variables
    integer :: MAX_IT, icount
    real(kind=DP) :: tolerance, TOL, alpha, beta, temp1, temp2, trace
    type(SPAM3_EMBED) :: rSinv, temp_trans
    type(SPAM3_EMBED) :: rSc,rSv
    integer :: is, ierr
    type(SPAM3_EMBED), allocatable, dimension(:,:) :: grad, grad_contra, r_contra,&
       r_vec, p_vec, r_prev

    call timer_clock('lr_tddft_utils_precond_RPA',1)

    MAX_IT=pub_lr_tddft_precond_iter
    TOL=pub_lr_tddft_precond_tol
    tolerance=1.0_DP
    icount=1
    pub_lr_tddft_preopt=.true.

    allocate(grad(2,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_precond_RPA','grad',ierr)
    allocate(grad_contra(2,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_precond_RPA','grad_contra',ierr)
    allocate(r_contra(2,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_precond_RPA','r_contra',ierr)
    allocate(r_vec(2,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_precond_RPA','r_vec',ierr)
    allocate(p_vec(2,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_precond_RPA','p_vec',ierr)
    allocate(r_prev(2,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_precond_RPA','r_prev',ierr)

    do is=1, pub_num_spins
       call sparse_embed_create(r_contra(1,is),g_vec(1,is))
       call sparse_embed_create(r_contra(2,is),g_vec(2,is))
       call sparse_embed_create(r_vec(1,is),f_vec(1,is))
       call sparse_embed_create(r_vec(2,is),f_vec(2,is))
       call sparse_embed_create(p_vec(1,is),g_vec(1,is))
       call sparse_embed_create(p_vec(2,is),g_vec(2,is))
       call sparse_embed_create(grad(1,is),f_vec(1,is))
       call sparse_embed_create(grad(2,is),f_vec(2,is))
       call sparse_embed_create(grad_contra(1,is),r_contra(1,is))
       call sparse_embed_create(grad_contra(2,is),r_contra(2,is))
       call sparse_embed_create(r_prev(1,is),r_vec(1,is))
       call sparse_embed_create(r_prev(2,is),r_vec(2,is))
    enddo
    call sparse_embed_create(rSinv, f_vec(1,1), val_rep%inv_overlap)
    call sparse_embed_transpose_structure(temp_trans%structure, g_vec(1,1))
    call sparse_embed_create(temp_trans, iscmplx=g_vec(1,1)%p%iscmplx)
    call sparse_embed_create(rSc,temp_trans,cond_rep%overlap)
    call sparse_embed_create(rSv, r_vec(1,1),val_rep%overlap)

    ! create initial g_vec as the contravariant version of f
    do is=1,pub_num_spins
       call sparse_embed_product(rSinv, f_vec(1,is), val_rep%inv_overlap)
       call sparse_embed_product(g_vec(1,is),cond_rep%inv_overlap,rSinv)
       call sparse_embed_product(rSinv, f_vec(2,is), val_rep%inv_overlap)
       call sparse_embed_product(g_vec(2,is),cond_rep%inv_overlap,rSinv)
    enddo
    ! calc operator on g_vec
    call linear_response_RPA_operator(grad, g_vec, cond_rep, &
         cond_ngwf_basis,cond_denskern, kchc, val_rep, &
         val_ngwf_basis, val_denskern, hvkv, mdl, hfxstate, response_dens,&
         sub_response_dens,fxc_fine,ground_state_dens,&
         sub_ground_state_dens,grad_contra)

    ! initialise r_0 and p_0:
    do is=1, pub_num_spins
       call sparse_embed_copy(r_vec(1,is), f_vec(1,is))
       call sparse_embed_axpy(r_vec(1,is), grad(1,is),-1.0_DP)
       call sparse_embed_copy(r_vec(2,is), f_vec(2,is))
       call sparse_embed_axpy(r_vec(2,is), grad(2,is),-1.0_DP)
       call sparse_embed_product(rSinv,r_vec(1,is),val_rep%inv_overlap)
       call sparse_embed_product(r_contra(1,is),cond_rep%inv_overlap,rSinv)
       call sparse_embed_copy(p_vec(1,is),r_contra(1,is))
       call sparse_embed_product(rSinv,r_vec(2,is),val_rep%inv_overlap)
       call sparse_embed_product(r_contra(2,is),cond_rep%inv_overlap,rSinv)
       call sparse_embed_copy(p_vec(2,is),r_contra(2,is))
    enddo

    ! start cg loop:
    do while(icount<MAX_IT .and. abs(tolerance)>TOL)
       ! calculate Pvec acting on operator:
       call linear_response_RPA_operator(grad, p_vec, cond_rep, &
            cond_ngwf_basis,cond_denskern, kchc, val_rep, &
            val_ngwf_basis, val_denskern, hvkv, mdl, hfxstate, response_dens, &
            sub_response_dens,fxc_fine, ground_state_dens,&
            sub_ground_state_dens,grad_contra)
       temp2=0.0_DP
       temp1=0.0_DP
       do is=1, pub_num_spins
          call sparse_embed_transpose(temp_trans,p_vec(1,is))
          call sparse_embed_trace(trace,temp_trans,grad(1,is))
          temp2=temp2+trace
          call sparse_embed_transpose(temp_trans,p_vec(2,is))
          call sparse_embed_trace(trace,temp_trans,grad(2,is))
          temp2=temp2+trace
          call sparse_embed_transpose(temp_trans,r_contra(1,is))
          call sparse_embed_trace(trace,temp_trans,r_vec(1,is))
          temp1=temp1+trace
          call sparse_embed_transpose(temp_trans,r_contra(2,is))
          call sparse_embed_trace(trace,temp_trans,r_vec(2,is))
          temp1=temp1+trace
       enddo
       alpha=temp1/temp2

       ! update g_vec and r_vec. Save old r_vec in r_prev.
       ! also construct r_contra
       do is=1,pub_num_spins
          call sparse_embed_copy(r_prev(1,is),r_vec(1,is))
          call sparse_embed_axpy(g_vec(1,is),p_vec(1,is),alpha)
          call sparse_embed_axpy(r_vec(1,is),grad(1,is),-1.0_DP*alpha)
          call sparse_embed_copy(r_prev(2,is),r_vec(2,is))
          call sparse_embed_axpy(g_vec(2,is),p_vec(2,is),alpha)
          call sparse_embed_axpy(r_vec(2,is),grad(2,is),-1.0_DP*alpha)
          call sparse_embed_product(rSinv,r_vec(1,is),val_rep%inv_overlap)
          call sparse_embed_product(r_contra(1,is),cond_rep%inv_overlap,rSinv)
          call sparse_embed_product(rSinv,r_vec(2,is),val_rep%inv_overlap)
          call sparse_embed_product(r_contra(2,is),cond_rep%inv_overlap,rSinv)
       enddo
       ! calculate stopping condition
       temp2=0.0_DP
       do is=1,pub_num_spins
          call sparse_embed_transpose(temp_trans, r_contra(1,is))
          call sparse_embed_trace(trace,temp_trans,r_vec(1,is))
          temp2=temp2+trace
          call sparse_embed_transpose(temp_trans, r_contra(2,is))
          call sparse_embed_trace(trace,temp_trans,r_vec(2,is))
          temp2=temp2+trace
       enddo
       tolerance=sqrt(abs(temp2))
       beta=temp2/temp1

       ! update p_vec and increment iteration number
       do is=1,pub_num_spins
          call sparse_embed_scale(p_vec(1,is),beta)
          call sparse_embed_axpy(p_vec(1,is),r_contra(1,is),1.0_DP)
          call sparse_embed_scale(p_vec(2,is),beta)
          call sparse_embed_axpy(p_vec(2,is),r_contra(2,is),1.0_DP)
       enddo
       icount=icount+1
    enddo

    ! once we have left the main loop, g_vec is the contravariant
    ! preconditioned version of f_vec

    final_tol=tolerance
    pub_lr_tddft_preopt=.false.

    ! deallocate data
    do is=1,pub_num_spins
       call sparse_embed_destroy(r_contra(1,is))
       call sparse_embed_destroy(r_vec(1,is))
       call sparse_embed_destroy(p_vec(1,is))
       call sparse_embed_destroy(grad(1,is))
       call sparse_embed_destroy(grad_contra(1,is))
       call sparse_embed_destroy(r_prev(1,is))
       call sparse_embed_destroy(r_contra(2,is))
       call sparse_embed_destroy(r_vec(2,is))
       call sparse_embed_destroy(p_vec(2,is))
       call sparse_embed_destroy(grad(2,is))
       call sparse_embed_destroy(grad_contra(2,is))
       call sparse_embed_destroy(r_prev(2,is))
    enddo
    call sparse_embed_destroy(rSinv)
    call sparse_embed_destroy(temp_trans)
    call sparse_embed_destroy(rSc)
    call sparse_embed_destroy(rSv)

    deallocate(grad, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_precond_RPA','grad',ierr)
    deallocate(grad_contra, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_precond_RPA','grad_contra',ierr)
    deallocate(r_contra, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_precond_RPA','r_contra',ierr)
    deallocate(r_vec, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_precond_RPA','r_vec',ierr)
    deallocate(p_vec, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_precond_RPA','p_vec',ierr)
    deallocate(r_prev, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_precond_RPA','r_prev',ierr)

    call timer_clock('lr_tddft_utils_precond_RPA',2)

  end subroutine lr_tddft_utils_precond_RPA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_utils_precond(g_vec,f_vec,response_dens,&
       sub_response_dens,fxc_fine,&
       cond_rep,cond_ngwf_basis,cond_denskern,kchc,val_rep,&
       val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,ground_state_dens,&
       sub_ground_state_dens,final_tol)

    !====================================================================!
    ! Subroutine applies the diagonal preconditioner of Kohn sham eigen- !
    ! value differences onto a single vector f_vec. Since the precon-    !
    ! ditioner is only diagonal in KS-eigenstate space, its action has   !
    ! to be computed by solving a linear system using a conjugate        !
    ! gradient algorithm.                                                !
    !====================================================================!

    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use linear_response, only: linear_response_operator
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_lr_tddft_preopt, pub_lr_tddft_precond_iter, &
         pub_lr_tddft_precond_tol, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_transpose, sparse_embed_trace,&
         sparse_embed_transpose_structure, sparse_embed_copy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check


    implicit none

    type(SPAM3_EMBED), intent(inout) :: f_vec(:) ! gradient of the full system
      ! with respect to a solution vector x_vec. jd: Dim: pub_num_spins. Leave ':' be.
    type(SPAM3_EMBED), intent(inout) :: g_vec(:) ! preconditioned version of f
                                               ! jd: Dim: pub_num_spins. Leave ':' be.
    type(MODEL), intent(in) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
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
    real(kind=DP), intent(inout) :: final_tol

    ! temporary variables
    integer :: MAX_IT, icount
    real(kind=DP) :: tolerance, TOL, alpha, beta, temp1, temp2, trace
    type(SPAM3_EMBED) :: rSinv, temp_trans
    type(SPAM3_EMBED) :: rSc,rSv
    integer :: is, ierr
    type(SPAM3_EMBED), allocatable, dimension(:) :: grad, grad_contra, r_contra,&
       r_vec, p_vec, r_prev

    call timer_clock('lr_tddft_utils_precond',1)

    MAX_IT=pub_lr_tddft_precond_iter
    TOL=pub_lr_tddft_precond_tol
    tolerance=1.0_DP
    icount=1
    pub_lr_tddft_preopt=.true.

    allocate(grad(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_precond','grad',ierr)
    allocate(grad_contra(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_precond','grad_contra',ierr)
    allocate(r_contra(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_precond','r_contra',ierr)
    allocate(r_vec(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_precond','r_vec',ierr)
    allocate(p_vec(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_precond','p_vec',ierr)
    allocate(r_prev(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_utils_precond','r_prev',ierr)

    do is=1, pub_num_spins
       call sparse_embed_create(r_contra(is),g_vec(is))
       call sparse_embed_create(r_vec(is),f_vec(is))
       call sparse_embed_create(p_vec(is),g_vec(is))
       call sparse_embed_create(grad(is),f_vec(is))
       call sparse_embed_create(grad_contra(is),r_contra(is))
       call sparse_embed_create(r_prev(is),r_vec(is))
    enddo
    call sparse_embed_create(rSinv, f_vec(1), val_rep%inv_overlap)
    call sparse_embed_transpose_structure(temp_trans%structure, g_vec(1))
    call sparse_embed_create(temp_trans, iscmplx=g_vec(1)%p%iscmplx)
    call sparse_embed_create(rSc,temp_trans,cond_rep%overlap)
    call sparse_embed_create(rSv, r_vec(1),val_rep%overlap)

    ! create initial g_vec as the contravariant version of f
    do is=1,pub_num_spins
       call sparse_embed_product(rSinv, f_vec(is), val_rep%inv_overlap)
       call sparse_embed_product(g_vec(is),cond_rep%inv_overlap,rSinv)
    enddo
    ! calc operator on g_vec
    call linear_response_operator(grad, g_vec, cond_rep, &
         cond_ngwf_basis,cond_denskern, kchc, val_rep, &
         val_ngwf_basis, val_denskern, hvkv, mdl, hfxstate, response_dens, &
         sub_response_dens, fxc_fine, ground_state_dens, &
         sub_ground_state_dens, grad_contra)

    ! initialise r_0 and p_0:
    do is=1, pub_num_spins
       call sparse_embed_copy(r_vec(is), f_vec(is))
       call sparse_embed_axpy(r_vec(is), grad(is),-1.0_DP)
       call sparse_embed_product(rSinv,r_vec(is),val_rep%inv_overlap)
       call sparse_embed_product(r_contra(is),cond_rep%inv_overlap,rSinv)
       call sparse_embed_copy(p_vec(is),r_contra(is))
    enddo

    ! start cg loop:
    do while(icount<MAX_IT .and. abs(tolerance)>TOL)
       ! calculate Pvec acting on operator:
       call linear_response_operator(grad, p_vec, cond_rep, &
            cond_ngwf_basis,cond_denskern, kchc, val_rep, &
            val_ngwf_basis, val_denskern, hvkv, mdl, hfxstate, response_dens, &
            sub_response_dens, fxc_fine, ground_state_dens, &
            sub_ground_state_dens, grad_contra)
       temp2=0.0_DP
       temp1=0.0_DP
       do is=1, pub_num_spins
          call sparse_embed_transpose(temp_trans,p_vec(is))
          call sparse_embed_trace(trace,temp_trans,grad(is))
          temp2=temp2+trace
          call sparse_embed_transpose(temp_trans,r_contra(is))
          call sparse_embed_trace(trace,temp_trans,r_vec(is))
          temp1=temp1+trace
       enddo
       alpha=temp1/temp2

       ! update g_vec and r_vec. Save old r_vec in r_prev.
       ! also construct r_contra
       do is=1,pub_num_spins
          call sparse_embed_copy(r_prev(is),r_vec(is))
          call sparse_embed_axpy(g_vec(is),p_vec(is),alpha)
          call sparse_embed_axpy(r_vec(is),grad(is),-1.0_DP*alpha)
          call sparse_embed_product(rSinv,r_vec(is),val_rep%inv_overlap)
          call sparse_embed_product(r_contra(is),cond_rep%inv_overlap,rSinv)
       enddo
       ! calculate stopping condition
       temp2=0.0_DP
       do is=1,pub_num_spins
          call sparse_embed_transpose(temp_trans, r_contra(is))
          call sparse_embed_trace(trace,temp_trans,r_vec(is))
          temp2=temp2+trace
       enddo
       tolerance=sqrt(abs(temp2))
       beta=temp2/temp1

       ! update p_vec and increment iteration number
       do is=1,pub_num_spins
          call sparse_embed_scale(p_vec(is),beta)
          call sparse_embed_axpy(p_vec(is),r_contra(is),1.0_DP)
       enddo
       icount=icount+1
    enddo

    ! once we have left the main loop, g_vec is the contravariant
    ! preconditioned version of f_vec

    final_tol=tolerance
    pub_lr_tddft_preopt=.false.

    ! deallocate data
    do is=1,pub_num_spins
       call sparse_embed_destroy(r_contra(is))
       call sparse_embed_destroy(r_vec(is))
       call sparse_embed_destroy(p_vec(is))
       call sparse_embed_destroy(grad(is))
       call sparse_embed_destroy(grad_contra(is))
       call sparse_embed_destroy(r_prev(is))
    enddo
    call sparse_embed_destroy(rSinv)
    call sparse_embed_destroy(temp_trans)
    call sparse_embed_destroy(rSc)
    call sparse_embed_destroy(rSv)

    deallocate(grad, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_precond','grad',ierr)
    deallocate(grad_contra, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_precond','grad_contra',ierr)
    deallocate(r_contra, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_precond','r_contra',ierr)
    deallocate(r_vec, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_precond','r_vec',ierr)
    deallocate(p_vec, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_precond','p_vec',ierr)
    deallocate(r_prev, stat=ierr)
    call utils_dealloc_check('lr_tddft_utils_precond','r_prev',ierr)

    call timer_clock('lr_tddft_utils_precond',2)

  end subroutine lr_tddft_utils_precond

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module lr_tddft_utils
