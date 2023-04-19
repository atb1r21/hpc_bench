! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!         Linear Response TDDFT RPA (full TDDFT) module          !
!                                                                !
! This module contains routines to calculate the linear res-     !
! ponse TDDFT eigenvalues of a given system. Here the full TDDFT !
! eigenvalue problem, rather than the simplified hermitian       !
! CIS eigenvalue problem in the Tamm-Dancoff approximation is    !
! solved. This requires the representation of both particle-hole !
! and hole-particle subspaces and the solution to a non-hermitian!
! eigenvalue problem.                                            !
!----------------------------------------------------------------!
! This module was created by Tim Zuehlsdorff in 2015.            !
!================================================================!

module lr_tddft_rpa

  use constants, only: DP

  implicit none

  private

  public :: lr_tddft_rpa_calculate

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine lr_tddft_rpa_calculate(vector_batch,lr_tddft_evals,&
       response_dens,sub_response_dens,fxc_fine,&
       cond_ngwf_basis,cond_rep,cond_denskern,kchc,cond_ham,cond_ngwf_num,&
       val_ngwf_basis,val_rep,val_denskern,hvkv,val_ham,val_ngwf_num,&
       proj_basis,nl_projectors,mdl,hfxstate,val_evecs,cond_evecs,&
       ground_state_dens,sub_ground_state_dens,joint_evals)
    !=================================================================!
    ! Main driver routine for calculating the RPA TDDFT energies of an!
    ! excitation. Initialises a TDDFT density matrix, minimises the   !
    ! Rayleigh-Ritz value and then computes properties of the         !
    ! excitation such as the oscillator strength. Returns the         !
    ! converged vector batch.                                         !
    ! Modified for embedding by Joseph Prentice, July 2018            !
    !=================================================================!

    use comms, only: pub_on_root
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
         dense_put_element
    use function_basis, only: FUNC_BASIS
    use constants, only: stdout, FINE_STRUCTURE,&
        HARTREE_IN_NS, ANGSTROM
    use hf_exchange, only: HFX_STATE
    use lr_tddft_utils, only: lr_tddft_utils_init_tddft_vecs,&
         lr_tddft_utils_calc_mlwfs,lr_tddft_utils_oscillator
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use projectors, only: PROJECTOR_SET
    use restart, only: restart_RPA_kernel_read
    use rundat, only: pub_lr_tddft_num_states, pub_lr_tddft_restart, &
         pub_lr_tddft_preopt, pub_lr_tddft_triplet,&
         pub_lr_tddft_analysis, pub_lr_tddft_write_densities,&
         pub_lr_tddft_num_conv_states, pub_num_spins, &
         pub_lr_tddft_subsystem_coupling, pub_lr_tddft_init_random,&
         pub_lr_tddft_mlwf_analysis, pub_lr_tddft_joint_set
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_copy
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_rand_gen
    use visual, only: visual_scalarfield

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: vector_batch(pub_lr_tddft_num_states,2,&
         pub_num_spins)
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
    type(SPAM3_EMBED), intent(in) :: cond_ham(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: val_ham(pub_num_spins)
    real(kind=DP), optional, intent(in) :: joint_evals(cond_ngwf_num,pub_num_spins)

    ! local variables:
    character(50) :: filename, tempfile
    integer :: icount, jcount, kcount, num_vecs, ierr, is
    integer :: num_conv_states
    type(SPAM3_EMBED) :: L_vec ! auxillary kernel
    type(DEM) :: initial_guess_vec
    real(kind=DP) :: temp_result
    real(kind=DP), allocatable, dimension(:) :: rand_array
    real(kind=DP), allocatable, dimension(:) :: oscillator
    real(kind=DP), allocatable, dimension(:) :: lifetime

    num_vecs=pub_lr_tddft_num_states
    num_conv_states=0
    ! allocate data storage for converged_evecs
    L_vec%structure= 'TDRA'//trim(cond_rep%postfix)
    call sparse_embed_create(L_vec, iscmplx=vector_batch(1,1,1)%p%iscmplx)

    call dense_create(initial_guess_vec, cond_ngwf_num, val_ngwf_num)

    ! check for restart option:
    if(pub_lr_tddft_restart .or. pub_lr_tddft_subsystem_coupling) then
      do icount=1, num_vecs
         call restart_RPA_kernel_read(vector_batch(icount,:,:),icount)
      enddo
    else
      if(pub_lr_tddft_num_conv_states>0) then
         do icount=1,pub_lr_tddft_num_conv_states
            call restart_RPA_kernel_read(vector_batch(icount,:,:),icount)
         enddo

      endif

      ! generate all other states randomly (if pub_lr_tddft_init_random: T)
      if(.not. pub_lr_tddft_init_random .and. pub_lr_tddft_num_conv_states==0) then
         ! initialise to KS states rather than random vecs
         call lr_tddft_utils_init_tddft_vecs(vector_batch(:,1,:),cond_ngwf_basis,&
              val_ngwf_basis,cond_rep,val_rep,cond_ham,val_ham,val_denskern,&
              val_ngwf_num,cond_ngwf_num,mdl,num_vecs)

         do icount=1,num_vecs
            do is=1, pub_num_spins
               call sparse_embed_copy(vector_batch(icount,2,is),&
                    vector_batch(icount,1,is))
            enddo
         enddo

      else
         allocate(rand_array(cond_ngwf_num*val_ngwf_num),&
               stat=ierr)
         call utils_alloc_check('lr_tddft_RPA_calculate','rand_array',ierr)

         do icount=1+pub_lr_tddft_num_conv_states, num_vecs
            ! create an array of random numbers for this kernel
            call utils_rand_gen(rand_array,cond_ngwf_num*&
               val_ngwf_num, icount)

            ! loop over elements in matrix
            do jcount=1, cond_ngwf_num
               do kcount=1, val_ngwf_num
                  temp_result=rand_array((jcount-1)*val_ngwf_num+kcount)
                  call dense_put_element(temp_result, initial_guess_vec, &
                       jcount,kcount)
               enddo
            enddo
            call dense_convert(L_vec, initial_guess_vec)
            do is=1, pub_num_spins
               call sparse_embed_copy(vector_batch(icount,1,is), L_vec)
               call sparse_embed_copy(vector_batch(icount,2,is), L_vec)
            enddo
            ! initialise p and q to the same vector (equivalent to setting Y to 0).
            ! thus first iteration step is a Tamm-Dancoff step)
         enddo
         deallocate(rand_array,stat=ierr)
         call utils_dealloc_check('lr_tddft_RPA_calculate','rand_array',ierr)
      endif
    endif

    ! if required, handle preoptimisation:
    if(pub_lr_tddft_preopt) then
       call lr_tddft_RPA_cg(vector_batch,lr_tddft_evals,mdl,hfxstate,&
            response_dens,sub_response_dens,fxc_fine,&
            cond_rep,cond_ngwf_basis,cond_denskern,kchc,&
            val_rep,val_ngwf_basis,val_denskern,hvkv,num_vecs,&
            ground_state_dens,sub_ground_state_dens)
       pub_lr_tddft_preopt=.false.

    endif

    ! now perform a cg optimisation
    call lr_tddft_RPA_cg(vector_batch,lr_tddft_evals,mdl,hfxstate,&
         response_dens,sub_response_dens,fxc_fine,&
         cond_rep,cond_ngwf_basis,cond_denskern,kchc,&
         val_rep,val_ngwf_basis,val_denskern,hvkv,num_vecs,&
         ground_state_dens,sub_ground_state_dens)


    allocate(oscillator(num_vecs),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_calculate','oscillator',ierr)
    allocate(lifetime(num_vecs),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_calculate','lifetime',ierr)


    ! now calculate oscillator strength
    ! compting in momentum rep only works if this is a joint rep
    ! calculation and lr_tddft_analysis is set to T, as we require
    ! KS eigenstates
    if(pub_lr_tddft_joint_set .and. pub_lr_tddft_analysis) then
       call lr_tddft_utils_oscillator(oscillator,lr_tddft_evals,&
         cond_ngwf_basis,cond_rep,val_ngwf_basis,val_rep,&
         proj_basis,nl_projectors,vector_batch(:,2,:),num_vecs,mdl,&
         cond_evecs,joint_evals)
    else
       call lr_tddft_utils_oscillator(oscillator,lr_tddft_evals,&
         cond_ngwf_basis,cond_rep,val_ngwf_basis,val_rep,&
         proj_basis,nl_projectors,vector_batch(:,2,:),num_vecs,mdl)
    endif

    ! For each excitation energy, calculate its oscillator strength
    if (pub_lr_tddft_triplet) then
       if (pub_on_root) write(stdout,'(/a80)')'|Excitation|    Energy &
            &(in Ha)   |     Oscillator Str.  | Lifetime (in ns)   '

       do icount=1, num_vecs
          lifetime(icount)=((1.0_DP/FINE_STRUCTURE)**3)/&
               (2.0_DP*lr_tddft_evals(icount)*lr_tddft_evals(icount)&
               *oscillator(icount))
          if (pub_on_root) write(stdout,'(i12,f18.8,e21.5,e21.5)') &
               icount, lr_tddft_evals(icount), 0.0_DP, 0.0_DP
       enddo
    else
       ! HEADER for excitation energies and oscillator strengths
       if (pub_on_root) write(stdout,'(/a80)')'|Excitation|    Energy &
            &(in Ha)   |     Oscillator Str.  | Lifetime (in ns)   '

       do icount=1, num_vecs
          lifetime(icount)=((1.0_DP/FINE_STRUCTURE)**3)/&
               (2.0_DP*lr_tddft_evals(icount)*lr_tddft_evals(icount)&
                *oscillator(icount))
          if (pub_on_root) write(stdout,'(i12,f18.8,e21.5,e21.5)') &
               icount, lr_tddft_evals(icount), oscillator(icount), &
               lifetime(icount)*HARTREE_IN_NS
       enddo
    endif

    if(pub_lr_tddft_analysis) then
       call lr_tddft_RPA_analysis(vector_batch, cond_rep, &
          cond_ngwf_basis, val_rep, val_ngwf_basis, val_evecs,&
          cond_evecs, num_vecs, val_ngwf_num, cond_ngwf_num)

       if(pub_lr_tddft_mlwf_analysis) then
          ! construct maximally localised wannier functions
          call lr_tddft_utils_calc_mlwfs(vector_batch(:,1,:),cond_ngwf_basis,&
               val_ngwf_basis,cond_rep,val_rep, mdl, proj_basis,&
               lr_tddft_evals,num_vecs,vector_batch(:,2,:))
       endif
    endif

    ! if required, write the e-h difference density
    if(pub_lr_tddft_write_densities) then

      do icount=1,num_vecs
         call lr_tddft_RPA_eh_dens(response_dens,vector_batch(icount,:,:), &
            cond_rep,cond_ngwf_basis,val_rep,val_ngwf_basis,mdl)
         write(filename,*) icount
         write(tempfile,*) '_TDDFT_eh_diff_dens_'
         write(filename,*) trim(adjustl(tempfile))//trim(adjustl(filename))

         if (pub_num_spins==2) then
            response_dens(:,:,:,1)=response_dens(:,:,:,1)+response_dens(:,:,:,2)
         endif

         call visual_scalarfield(response_dens(:,:,:,1),mdl%fine_grid,&
               mdl%cell, 'eh difference density (in e/ang^3) for:',filename, &
               mdl%elements, ANGSTROM**3)

      enddo
    endif

    ! deallocate data
    call dense_destroy(initial_guess_vec)
    call sparse_embed_destroy(L_vec)
    deallocate(oscillator,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_calculate','oscillator',ierr)
    deallocate(lifetime,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_calculate','lifetime',ierr)

  end subroutine lr_tddft_RPA_calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_RPA_eh_dens(dens,vector, cond_rep,&
      cond_ngwf_basis,val_rep,val_ngwf_basis,mdl)

    !===================================================================!
    ! Subroutine constructs the appropriate e-h difference density      !
    ! for an excitation described by 'vector'.                          !
    ! Modified for embedding by Joseph Prentice, July 2018              !
    !===================================================================!

    use constants, only: stdout
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use lr_tddft_utils, only: lr_tddft_utils_eh_denskern_RPA
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_num_spins, pub_lr_tddft_joint_set
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_transpose, sparse_embed_transpose_structure, sparse_embed_copy, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_destroy
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    real(kind=DP), intent(inout) :: dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: vector(2,pub_num_spins)
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis(:)
    type(NGWF_REP), intent(in) :: cond_rep
    type(NGWF_REP), intent(in) :: val_rep
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)

    ! allocatable variables
    integer :: ierr, is, isub, jsub
    type(SPAM3_EMBED) :: temp_trans, PSphi, SchiP,elec2,hole2
    type(SPAM3_EMBED), allocatable, dimension(:) :: eh_kernel,elec,hole
    type(SPAM3), allocatable, dimension(:) :: kern_array
    real(kind=DP), allocatable, dimension(:,:,:,:) :: dens_hole, sub_dens

    ! allocate variables
    allocate(eh_kernel(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_eh_dens','eh_kernel',ierr)
    allocate(hole(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_eh_dens','eh_kernel',ierr)
    allocate(elec(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_eh_dens','eh_kernel',ierr)
    allocate(dens_hole(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_eh_dens','dens_hole',ierr)
    allocate(sub_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_eh_dens','sub_dens',ierr)
    allocate(kern_array(pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_eh_dens','kern_array',ierr)
    call sparse_embed_transpose_structure(temp_trans%structure,vector(1,1))
    call sparse_embed_create(temp_trans, iscmplx=vector(1,1)%p%iscmplx)
    call sparse_embed_create(PSphi,vector(1,1),val_rep%overlap)
    call sparse_embed_create(SchiP,cond_rep%overlap,vector(1,1))
    do is=1, pub_num_spins
       call sparse_embed_create(elec(is),PSphi,temp_trans)
       call sparse_embed_create(hole(is),temp_trans,SchiP)
       call sparse_embed_create(eh_kernel(is),elec(is))
    enddo
    call sparse_embed_create(hole2,hole(1))
    call sparse_embed_create(elec2,elec(1))

    ! if the joint set is used for representing the conduction space, this
    ! is simple:
    dens_hole=0.d0
    if(pub_lr_tddft_joint_set) then
       call lr_tddft_utils_eh_denskern_RPA(eh_kernel,vector,cond_rep%overlap,&
            val_rep%overlap,cond_rep%inv_overlap,val_rep%inv_overlap,cond_rep%cross_overlap)
       ! jcap: loop over regions
       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub

             sub_dens=0.d0

             call sparse_embed_extract_from_array(kern_array,eh_kernel,isub,jsub)

             call density_on_grid(sub_dens,mdl%fine_grid,mdl%dbl_grid, &
                  mdl%cell, mdl%fftbox, kern_array, &
                  cond_rep%ngwf_overlap%m(isub,jsub), &
                  cond_rep%ngwfs_on_grid(isub), cond_ngwf_basis(isub), &
                  cond_rep%ngwfs_on_grid(jsub), cond_ngwf_basis(jsub))

             dens=dens+sub_dens

             call sparse_embed_destroy_extracted_array(kern_array)

          end do
       end do
    else
       ! now construct elec and hole denskerns assuming cond rep is used
       ! for conduction space rather than joint rep
       do is=1, pub_num_spins
          call sparse_embed_transpose(temp_trans,vector(1,is))
          call sparse_embed_product(PSphi,vector(1,is),val_rep%overlap)
          call sparse_embed_product(SchiP,cond_rep%overlap,vector(1,is))
          call sparse_embed_product(elec(is),PSphi,temp_trans)
          call sparse_embed_product(hole(is),temp_trans,SchiP)
          call sparse_embed_transpose(temp_trans,vector(2,is))
          call sparse_embed_product(PSphi,vector(2,is),val_rep%overlap)
          call sparse_embed_product(SchiP,cond_rep%overlap,vector(2,is))
          call sparse_embed_product(elec2,PSphi,temp_trans)
          call sparse_embed_product(hole2,temp_trans,SchiP)
          call sparse_embed_axpy(elec(is),elec2,1.0_DP)
          call sparse_embed_axpy(hole(is),hole2,1.0_DP)
          call sparse_embed_scale(elec(is),0.5_DP)
          call sparse_embed_scale(hole(is),0.5_DP)
          call sparse_embed_copy(eh_kernel(is),elec(is))
       enddo

       ! jcap: loop over regions
       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub

             sub_dens=0.d0

             call sparse_embed_extract_from_array(kern_array,eh_kernel,isub,jsub)

             ! create e-h kernel. Now construct elec density
             call density_on_grid(sub_dens,mdl%fine_grid,mdl%dbl_grid, &
                  mdl%cell, mdl%fftbox, kern_array, &
                  cond_rep%ngwf_overlap%m(isub,jsub), &
                  cond_rep%ngwfs_on_grid(isub), cond_ngwf_basis(isub), &
                  cond_rep%ngwfs_on_grid(jsub), cond_ngwf_basis(jsub))

             dens=dens+sub_dens

             call sparse_embed_destroy_extracted_array(kern_array)

          end do
       end do

       do is=1, pub_num_spins
          call sparse_embed_destroy(eh_kernel(is))
          call sparse_embed_create(eh_kernel(is),hole(is))
          call sparse_embed_copy(eh_kernel(is),hole(is))
       enddo

       ! jcap: loop over regions
       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub

             sub_dens=0.d0

             call sparse_embed_extract_from_array(kern_array,eh_kernel,isub,jsub)

             ! now construct hole density
             call density_on_grid(sub_dens,mdl%fine_grid,mdl%dbl_grid, &
                  mdl%cell, mdl%fftbox, kern_array, &
                  val_rep%ngwf_overlap%m(isub,jsub), &
                  val_rep%ngwfs_on_grid(isub), val_ngwf_basis(isub), &
                  val_rep%ngwfs_on_grid(jsub), val_ngwf_basis(jsub))

             dens_hole=dens_hole+sub_dens

             call sparse_embed_destroy_extracted_array(kern_array)

          end do
       end do

       ! difference density is electron density minus hole density
       dens=dens-dens_hole
    endif

    ! done. Deallocate all data
    do is=1, pub_num_spins
       call sparse_embed_destroy(eh_kernel(is))
       call sparse_embed_destroy(elec(is))
       call sparse_embed_destroy(hole(is))
    enddo
    deallocate(eh_kernel, stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_eh_dens','eh_kernel',ierr)
    deallocate(hole, stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_eh_dens','hole',ierr)
    deallocate(elec, stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_eh_dens','elec',ierr)
    deallocate(dens_hole,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_eh_dens','dens_hole',ierr)
    call sparse_embed_destroy(elec2)
    call sparse_embed_destroy(hole2)
    call sparse_embed_destroy(temp_trans)
    call sparse_embed_destroy(PSphi)
    call sparse_embed_destroy(SchiP)

  end subroutine lr_tddft_RPA_eh_dens

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_RPA_cg(vector_batch,lr_tddft_evals,mdl,hfxstate,&
       response_dens,sub_response_dens,fxc_fine,&
       cond_rep,cond_ngwf_basis,cond_denskern,kchc,&
       val_rep,val_ngwf_basis,val_denskern,hvkv,num_states,&
       ground_state_dens,sub_ground_state_dens)
    !================================================================!
    ! Routine minimises the Rayleigh-Ritz value of the RPA problem   !
    ! and returns both the individual LR-TDDFT eigenvalues and the   !
    ! corresponing eigenvectors.                                     !
    ! Modified for embedding by Joseph Prentice, July 2018           !
    !================================================================!
    use comms, only: pub_on_root
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
         dense_put_element
    use function_basis, only: FUNC_BASIS
    use constants, only: stdout, NORMAL, BRIEF, VERBOSE
    use model_type, only: MODEL
    use hartree, only: hartree_via_multigrid
    use hf_exchange, only: HFX_STATE
    use is_solvation, only: have_initial_eps
    use linalg, only: linalg_dsygv_lt
    use linear_response, only: linear_response_RPA_operator
    use lr_tddft_utils, only: lr_tddft_utils_precond_batch_RPA,&
         lr_tddft_utils_penalty_batch_RPA
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_output_detail, &
         pub_lr_tddft_cg_threshold, &
         pub_lr_tddft_maxit_cg, pub_lr_tddft_reset_cg, &
         pub_lr_tddft_write_kernels, pub_lr_tddft_preopt, &
         pub_lr_tddft_preopt_iter, pub_lr_tddft_num_conv_states,&
         pub_lr_tddft_precond, &
         pub_lr_tddft_sparse_region, pub_lr_tddft_restart, &
         pub_num_spins, pub_multigrid_hartree, pub_lr_optical_permittivity,&
         pub_is_bulk_permittivity,pub_lr_tddft_maxit_pen,&
         pub_lr_tddft_penalty_tol, pub_print_qc
    use restart, only: restart_RPA_kernel_write, restart_RPA_kernel_read
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_transpose, sparse_embed_transpose_structure, sparse_embed_copy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_banner, &
         utils_qc_print

    ! Arguments
    type(MODEL), intent(in) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    integer, intent(in) :: num_states  ! number of converged excitation energies
    type(SPAM3_EMBED), intent(inout) :: vector_batch(num_states,2,pub_num_spins)
    real(kind=DP), intent(inout) :: lr_tddft_evals(num_states) ! converged eval
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
    real(kind=DP), allocatable,dimension(:,:) :: identity_mat,subspace_mat
    real(kind=DP) :: omega, omega_prev, lambda, temp1,&
      temp_total,penalty_value
    real(kind=DP) :: temp_val,temp_val2,temp_val3
    real(kind=DP), allocatable, dimension (:,:,:,:,:,:) :: ground_state_hartree
    type(SPAM3_EMBED) :: eff_fvec
    type(SPAM3_EMBED), allocatable, dimension (:,:,:) :: y_vec, y_contra_vec, &
        grad_vec, grad_prev_vec,g_vec, f_prev_vec, f_vec,y_prev_vec,&
        vector_batch_cov ! all these are contravariant,
        ! apart from y batch and vector_batch_cov
    type(SPAM3_EMBED), allocatable, dimension (:) :: L_vec ! auxillary kernel
    type(SPAM3_EMBED) :: LSvkv, PSv, temp_Svkv
    type(SPAM3_EMBED), allocatable, dimension(:) :: kcSc, Svkv
    type(SPAM3_EMBED), allocatable, dimension(:) :: temp_g
    integer :: maxit, icount,jcount,counter, ierr, is, isub, jsub
    real(kind=DP) :: tolerance,gamma_current, gamma_prev, gamma1,gamma2
    logical :: converged, lr_tddft_precond
    integer :: reset_counter


    lr_tddft_precond=pub_lr_tddft_precond
    if(pub_lr_tddft_preopt) then
       maxit=pub_lr_tddft_preopt_iter
       tolerance=pub_lr_tddft_cg_threshold
    else
       maxit=pub_lr_tddft_maxit_cg
       tolerance=pub_lr_tddft_cg_threshold
    endif

    ! allocate storage space
    if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
      allocate(ground_state_hartree(mdl%fine_grid%ld1,&
           mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,&
           mdl%nsub,mdl%nsub), stat=ierr)
      call utils_alloc_check('lr_tddft_RPA_cg', 'ground_state_hartree', ierr)
    endif
    allocate(Svkv(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','Svkv',ierr)
    allocate(kcSc(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','kcSc',ierr)
    allocate(L_vec(2),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','L_vec',ierr)
    allocate(y_vec(num_states,2,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','y_vec',ierr)
    allocate(y_contra_vec(num_states,2,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','y_contra_vec',ierr)
    allocate(grad_vec(num_states,2,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','grad_vec',ierr)
    allocate(grad_prev_vec(num_states,2,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','grad_prev_vec',ierr)
    allocate(g_vec(num_states,2,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','g_vec',ierr)
    allocate(f_prev_vec(num_states,2,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','f_prev_vec',ierr)
    allocate(y_prev_vec(num_states,2,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','y_prev_vec',ierr)
    allocate(f_vec(num_states,2,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','f_vec',ierr)
    allocate(vector_batch_cov(num_states,2,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','vector_batch_cov',ierr)
    allocate(temp_g(2),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','temp_g',ierr)

    ! allocate matrices:
    call sparse_embed_create(PSv,vector_batch(1,1,1),val_rep%overlap)
    L_vec(1)%structure='TDRA'//trim(cond_rep%postfix)
    L_vec(2)%structure='TDRA'//trim(cond_rep%postfix)
    call sparse_embed_create(L_vec(1), iscmplx=vector_batch(1,1,1)%p%iscmplx)
    call sparse_embed_create(L_vec(2), iscmplx=vector_batch(1,1,1)%p%iscmplx)
    do icount=1, num_states
       do is=1, pub_num_spins
          call sparse_embed_create(y_vec(icount,1,is),cond_rep%overlap,PSv)
          call sparse_embed_create(y_vec(icount,2,is),cond_rep%overlap,PSv)
          call sparse_embed_create(vector_batch_cov(icount,1,is),y_vec(1,1,is))
          call sparse_embed_create(vector_batch_cov(icount,2,is),y_vec(1,1,is))
          call sparse_embed_create(y_contra_vec(icount,1,is),vector_batch(1,1,is))
          call sparse_embed_create(y_contra_vec(icount,2,is),vector_batch(1,1,is))
          call sparse_embed_create(y_prev_vec(icount,1,is),vector_batch(1,1,is))
          call sparse_embed_create(y_prev_vec(icount,2,is),vector_batch(1,1,is))
          call sparse_embed_create(grad_vec(icount,1,is),vector_batch(1,1,is))
          call sparse_embed_create(grad_vec(icount,2,is),vector_batch(1,1,is))
          call sparse_embed_create(grad_prev_vec(icount,1,is),vector_batch(1,1,is))
          call sparse_embed_create(grad_prev_vec(icount,2,is),vector_batch(1,1,is))
          call sparse_embed_create(g_vec(icount,1,is),vector_batch(1,1,is))
          call sparse_embed_create(g_vec(icount,2,is),vector_batch(1,1,is))
          call sparse_embed_create(f_prev_vec(icount,1,is),vector_batch(1,1,is))
          call sparse_embed_create(f_prev_vec(icount,2,is),vector_batch(1,1,is))
          call sparse_embed_create(f_vec(icount,1,is),vector_batch(1,1,is))
          call sparse_embed_create(f_vec(icount,2,is),vector_batch(1,1,is))
       enddo
    enddo
    call sparse_embed_create(eff_fvec,y_vec(1,1,1))
    do is=1, pub_num_spins
       call sparse_embed_create(Svkv(is),val_rep%overlap,val_denskern(1))
       call sparse_embed_create(kcSc(is),cond_denskern(1),cond_rep%overlap)
    enddo
    call sparse_embed_create(LSvkv,L_vec(1),Svkv(1))
    call sparse_embed_create(temp_Svkv,vector_batch(1,1,1),Svkv(1))
    do is=1, pub_num_spins
       call sparse_embed_product(Svkv(is),val_rep%overlap,val_denskern(is))
       call sparse_embed_product(kcSc(is),cond_denskern(is),cond_rep%overlap)
    enddo
    call sparse_embed_create(temp_g(1),g_vec(1,1,1))
    call sparse_embed_create(temp_g(2),g_vec(1,1,1))

    call timer_clock('lr_tddft_RPA_cg',1)

    converged=.false.
    reset_counter=1

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
    end if


    if (pub_on_root .and. pub_output_detail==BRIEF) write(stdout,'(/a80)')&
       '|Iter    |  Energy (in Ha) |  Delta E (in Ha)  |  Step Length | Penalty '

    ! create ground state hartree potential if required
    ! only required if TDDFT is performed within the implicit solvent
    ! model
    if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
      ! change bulk permittivity to have the correct value for optical
      ! frequencies
      pub_is_bulk_permittivity=pub_lr_optical_permittivity
      have_initial_eps=.false.

      if(pub_num_spins==1) then
         sub_ground_state_dens=2.0_DP*sub_ground_state_dens
      endif
      ! jcap: loop over regions
      do isub=1,mdl%nsub
         do jsub=1,mdl%nsub
            call hartree_via_multigrid(ground_state_hartree(:,:,:,:,isub,jsub),&
                 sub_ground_state_dens(:,:,:,:,isub,jsub),&
                 mdl%fine_grid,mdl%cell,elements=mdl%elements)
         end do
      end do
      if(pub_num_spins==1) then
         sub_ground_state_dens=0.5_DP*sub_ground_state_dens
      endif
    endif

    ! project initial vector into physical space
    if(.not. (pub_lr_tddft_sparse_region .and. pub_lr_tddft_restart)) then
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          do is=1, pub_num_spins
             call sparse_embed_copy(L_vec(1),vector_batch(icount,1,is))
             call sparse_embed_copy(L_vec(2),vector_batch(icount,2,is))
             call sparse_embed_product(LSvkv, L_vec(1), Svkv(is))
             call sparse_embed_product(vector_batch(icount,1,is), kcSc(is), LSvkv)
             call sparse_embed_product(LSvkv, L_vec(2), Svkv(is))
             call sparse_embed_product(vector_batch(icount,2,is), kcSc(is), LSvkv)
          enddo
       enddo
    endif

    ! create covariant vector batch for all vectors (also already converged ones)
    do icount=1, num_states
       do is=1, pub_num_spins
          call sparse_embed_product(PSv,vector_batch(icount,1,is),val_rep%overlap)
          call sparse_embed_product(vector_batch_cov(icount,1,is),cond_rep%overlap,PSv)
          call sparse_embed_product(PSv,vector_batch(icount,2,is),val_rep%overlap)
          call sparse_embed_product(vector_batch_cov(icount,2,is),cond_rep%overlap,PSv)
       enddo
    enddo

    ! orthonormalise vector batch
    call lr_tddft_RPA_gram_schmidt(vector_batch,vector_batch_cov, num_states,&
         cond_rep%overlap, val_rep%overlap)

    ! Measure the penalty functional
    call lr_tddft_utils_penalty_batch_RPA(penalty_value, vector_batch, &
            vector_batch_cov,cond_denskern, val_denskern, &
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
       if (pub_on_root.and.pub_output_detail==VERBOSE) then
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
             call sparse_embed_product(temp_Svkv, vector_batch(jcount,1,is), Svkv(is))
             call sparse_embed_product(vector_batch(jcount,1,(is)), kcSc(is), temp_Svkv)
             call sparse_embed_product(temp_Svkv, vector_batch(jcount,2,is), Svkv(is))
             call sparse_embed_product(vector_batch(jcount,2,is), kcSc(is), temp_Svkv)
             call sparse_embed_product(PSv,vector_batch(jcount,1,is),val_rep%overlap)
             call sparse_embed_product(vector_batch_cov(jcount,1,is),cond_rep%overlap,PSv)
             call sparse_embed_product(PSv,vector_batch(jcount,2,is),val_rep%overlap)
             call sparse_embed_product(vector_batch_cov(jcount,2,is),cond_rep%overlap,PSv)
          enddo
       enddo

       call lr_tddft_utils_penalty_batch_RPA(penalty_value,vector_batch, &
            vector_batch_cov,cond_denskern, val_denskern, &
            cond_rep%overlap, val_rep%overlap,Svkv,kcSc,&
            pub_lr_tddft_num_conv_states,num_states)

       ! OUTPUT iterative improvement
       if (pub_on_root.and.pub_output_detail == VERBOSE) write(stdout,&
            '(i6,e12.4)') icount, penalty_value
       icount=icount+1
    enddo

    ! check for convergence of iterative improvement
    if (penalty_value>pub_lr_tddft_penalty_tol) then
       if (pub_on_root.and.pub_output_detail >=NORMAL) &
            write(stdout,'(a,i4,a)')  'WARNING: maximum number of &
            &iterative improvements for K^{1} (',pub_lr_tddft_maxit_pen,&
            ') exceeded!'
    endif

    ! calculate intial y_vec and y_contra_vec
    do icount=1+pub_lr_tddft_num_conv_states, num_states
       if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
          call linear_response_RPA_operator(y_vec(icount,:,:), vector_batch(icount,:,:),&
               cond_rep, cond_ngwf_basis, cond_denskern, kchc, val_rep, &
               val_ngwf_basis, val_denskern,hvkv,mdl,hfxstate,response_dens,&
               sub_response_dens,fxc_fine, ground_state_dens,&
               sub_ground_state_dens,y_contra_vec(icount,:,:), &
               ground_state_hartree=ground_state_hartree)
       else
          call linear_response_RPA_operator(y_vec(icount,:,:), vector_batch(icount,:,:),&
               cond_rep, cond_ngwf_basis, cond_denskern, kchc, val_rep, &
               val_ngwf_basis, val_denskern,hvkv,mdl,hfxstate,response_dens,&
               sub_response_dens,fxc_fine, ground_state_dens,&
               sub_ground_state_dens,y_contra_vec(icount,:,:))
       endif
    enddo

    ! calculate omega
    call lr_tddft_RPA_calc_omega(omega,vector_batch,vector_batch_cov,y_vec,&
        cond_rep%overlap,val_rep%overlap,num_states)

    ! set omega_prev:
    omega_prev=omega+0.1_DP ! add 0.1_DP to ensure tolerance is not
    ! met in first iteration

    if (pub_output_detail >= NORMAL) then
       if (pub_on_root) then
          write(stdout,'(/a)')'           ............................&
               &............................'
          write(stdout,'(a)')'           |     *** TDDFT optimisation &
               &not converged ***         |'
          write(stdout,'(a,f19.14,a)') &
               '           |       LR-TDDFT energy  =  ',omega,'     &
               &   |'
          write(stdout,'(a)') '           ============================&
               &============================'
       endif
    endif

    if (pub_on_root .and. pub_output_detail==BRIEF) then
       write(stdout,'(i16,f18.8,f18.8,f16.8,e15.4)') counter, omega,&
           1.0_DP, 1.0_DP, penalty_value
    endif

    counter=2

    ! start main loop:
    do while(abs(omega-omega_prev)>tolerance.and.counter<maxit)
       ! PRINT iteration header
       if (pub_on_root .and. pub_output_detail >= NORMAL) then
          write(stdout,'(/a)') repeat('#',80)
          write(stdout,'(a,i4,a)')'######################### LR_TDDFT &
               &CG iteration ',counter,' ###########################'
          write(stdout,'(a,E12.4,a)') &
               '****** Initial penalty functional value of K^{1}: &
               &Q[K^{1}] = ', penalty_value, ' ******'
       endif

       ! first calculate f_vec. fvec is the full gradient:
       !f_p=y_p-Tr[]
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          do is=1, pub_num_spins
             call sparse_embed_copy(f_vec(icount,1,is),y_contra_vec(icount,1,is))
             call sparse_embed_copy(f_vec(icount,2,is),y_contra_vec(icount,2,is))
          enddo
          call lr_tddft_RPA_inner_product(temp_val,vector_batch(icount,1,:),&
               y_vec(icount,1,:),cond_rep%overlap,val_rep%overlap,.false.)
          call lr_tddft_RPA_inner_product(temp_val2,vector_batch(icount,2,:),&
               y_vec(icount,2,:),cond_rep%overlap,val_rep%overlap,.false.)
          call lr_tddft_RPA_inner_product(temp_val3,vector_batch(icount,1,:),&
               vector_batch_cov(icount,2,:),cond_rep%overlap,val_rep%overlap,.false.)

          temp_total=(temp_val+temp_val2)/(2.0_DP*temp_val3)

          do is=1, pub_num_spins
             call sparse_embed_axpy(f_vec(icount,1,is),vector_batch(icount,2,is),&
                 -1.0_DP*temp_total)
             call sparse_embed_axpy(f_vec(icount,2,is),vector_batch(icount,1,is),&
                 -1.0_DP*temp_total)
          enddo
       enddo


       ! now orthogonalise f_vec to all current x vecs. Note that this is
       ! an RPA orthonormalisation. Thus it is an orthonormalisation with
       ! respect to the effective RPA normalisation condition
       call lr_tddft_RPA_orthogonalise(f_vec,vector_batch,vector_batch_cov,&
             num_states,cond_rep%overlap,val_rep%overlap)


       ! now precondition g_batch
       if(lr_tddft_precond .and. .not. pub_lr_tddft_preopt) then
          call lr_tddft_utils_precond_batch_RPA(g_vec,f_vec,&
               num_states, pub_lr_tddft_num_conv_states,&
               response_dens,sub_response_dens,fxc_fine,&
               cond_rep,cond_ngwf_basis,cond_denskern,kchc,&
               val_rep, val_ngwf_basis,val_denskern,hvkv,&
               mdl,hfxstate,ground_state_dens,sub_ground_state_dens)
       else
          do icount=1+pub_lr_tddft_num_conv_states, num_states
             do is=1, pub_num_spins
                call sparse_embed_copy(g_vec(icount,1,is),f_vec(icount,1,is))
                call sparse_embed_copy(g_vec(icount,2,is),f_vec(icount,2,is))
             enddo
          enddo
       endif

       ! reorthogonalise wrt to all current vectors in batch
       ! only necessary if the gradient is preconditioned
       if(lr_tddft_precond .and. .not. pub_lr_tddft_preopt) then
          call lr_tddft_RPA_orthogonalise(g_vec,vector_batch,vector_batch_cov,&
               num_states,cond_rep%overlap,val_rep%overlap)
       endif

       ! calculate gamma1 and gamma2 of conjugate gradients algorithm
       gamma1=0.0_DP
       gamma2=0.0_DP
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          call lr_tddft_RPA_inner_product(temp_val,g_vec(icount,1,:),&
               f_vec(icount,1,:),cond_rep%overlap,val_rep%overlap,.true.)
          do is=1, pub_num_spins
             gamma1=gamma1+temp_val
          enddo
          call lr_tddft_RPA_inner_product(temp_val,g_vec(icount,2,:),&
               f_vec(icount,2,:),cond_rep%overlap,val_rep%overlap,.true.)
          do is=1, pub_num_spins
             gamma1=gamma1+temp_val
          enddo

          ! now calculate gamma2 unless it is first iteration or conjugate
          ! gradient reset
          if(counter>2 .and. reset_counter<pub_lr_tddft_reset_cg+1) then
             call lr_tddft_RPA_inner_product(temp_val,g_vec(icount,1,:),&
               f_prev_vec(icount,1,:),cond_rep%overlap,val_rep%overlap,.true.)
               do is=1, pub_num_spins
                  gamma2=gamma2+temp_val
               enddo
             call lr_tddft_RPA_inner_product(temp_val,g_vec(icount,2,:),&
               f_prev_vec(icount,2,:),cond_rep%overlap,val_rep%overlap,.true.)
               do is=1, pub_num_spins
                  gamma2=gamma2+temp_val
               enddo
          else
             gamma_prev=1.0_DP
          endif
       enddo

       if (pub_on_root .and. pub_output_detail == VERBOSE) then
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

       do icount=1+pub_lr_tddft_num_conv_states, num_states
          do is=1, pub_num_spins
             call sparse_embed_copy(grad_vec(icount,1,is), g_vec(icount,1,is))
             call sparse_embed_copy(grad_vec(icount,2,is), g_vec(icount,2,is))
             call sparse_embed_scale(grad_vec(icount,1,is),-1.0_DP)
             call sparse_embed_scale(grad_vec(icount,2,is),-1.0_DP)

             ! completely comment out conjugate gradient for time being
             if (counter>2 .and. reset_counter< pub_lr_tddft_reset_cg+1) then
                call sparse_embed_axpy(grad_vec(icount,1,is),grad_prev_vec(icount,1,is), &
                     gamma_current)
                call sparse_embed_axpy(grad_vec(icount,2,is),grad_prev_vec(icount,2,is),&
                     gamma_current)
             endif
          enddo
       enddo


       ! now project the conjugate gradient vector into physical space
       ! This is needed for the subsystem tddft approach
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          do is=1, pub_num_spins
             call sparse_embed_copy(L_vec(1),grad_vec(icount,1,is))
             call sparse_embed_copy(L_vec(2),grad_vec(icount,2,is))
             call sparse_embed_product(LSvkv,L_vec(1),Svkv(is))
             call sparse_embed_product(grad_vec(icount,1,is),kcSc(is), LSvkv)
             call sparse_embed_product(LSvkv,L_vec(2),Svkv(is))
             call sparse_embed_product(grad_vec(icount,2,is),kcSc(is),LSvkv)
          enddo
       enddo

       ! reorthogonalise grad batch
       call lr_tddft_RPA_orthogonalise(grad_vec,vector_batch,vector_batch_cov,&
              num_states,cond_rep%overlap,val_rep%overlap)

       ! now do a line minimisation
       if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
          call lr_tddft_RPA_line_min(lambda, vector_batch,vector_batch_cov,&
            y_vec,grad_vec,cond_rep,cond_ngwf_basis,cond_denskern,kchc,val_rep,&
            val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,response_dens,&
            sub_response_dens,fxc_fine,L_vec,num_states,ground_state_dens,&
            sub_ground_state_dens,ground_state_hartree)
       else
          call lr_tddft_RPA_line_min(lambda, vector_batch,vector_batch_cov,&
            y_vec,grad_vec,cond_rep,cond_ngwf_basis,cond_denskern,kchc,val_rep,&
            val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,response_dens,&
            sub_response_dens,fxc_fine,L_vec,num_states,ground_state_dens,&
            sub_ground_state_dens)
       endif

       call lr_tddft_utils_penalty_batch_RPA(penalty_value,vector_batch, &
               vector_batch_cov,cond_denskern, val_denskern, &
               cond_rep%overlap, val_rep%overlap,Svkv,kcSc,&
               pub_lr_tddft_num_conv_states,num_states)

       ! HEADER for penalty func improvement
       if (penalty_value>pub_lr_tddft_penalty_tol) then
          if (pub_on_root.and.pub_output_detail==VERBOSE) then
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
                call sparse_embed_product(temp_Svkv, vector_batch(jcount,1,is), &
                     Svkv(is))
                call sparse_embed_product(vector_batch(jcount,1,is), kcSc(is), &
                     temp_Svkv)
                call sparse_embed_product(temp_Svkv, vector_batch(jcount,2,is), &
                     Svkv(is))
                call sparse_embed_product(vector_batch(jcount,2,is), kcSc(is), &
                     temp_Svkv)
                call sparse_embed_product(PSv,vector_batch(jcount,1,is), &
                     val_rep%overlap)
                call sparse_embed_product(vector_batch_cov(jcount,1,is), &
                     cond_rep%overlap,PSv)
                call sparse_embed_product(PSv,vector_batch(jcount,2,is), &
                     val_rep%overlap)
                call sparse_embed_product(vector_batch_cov(jcount,2,is), &
                     cond_rep%overlap,PSv)
             enddo
          enddo

          call lr_tddft_utils_penalty_batch_RPA(penalty_value,vector_batch, &
               vector_batch_cov,cond_denskern, val_denskern, &
               cond_rep%overlap, val_rep%overlap,Svkv,kcSc,&
               pub_lr_tddft_num_conv_states,num_states)

          ! OUTPUT iterative improvement
          if (pub_on_root.and.pub_output_detail == VERBOSE) write(stdout,&
               '(i6,e12.4)') icount, penalty_value
          icount=icount+1
       enddo

       ! check for convergence of iterative improvement
       if (penalty_value>pub_lr_tddft_penalty_tol) then
          if (pub_on_root.and.pub_output_detail >=NORMAL) &
               write(stdout,'(a,i4,a)')  'WARNING: maximum number of &
               &iterative improvements for K^{1} (',pub_lr_tddft_maxit_pen,&
               ') exceeded!'
       endif

       ! ortonormalise the new set of vectors and calcualte the new y batch
       call lr_tddft_RPA_gram_schmidt(vector_batch,vector_batch_cov,num_states,&
           cond_rep%overlap, val_rep%overlap)

       ! calculate intial y_batch and y_contra_batch
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
             call linear_response_RPA_operator(y_vec(icount,:,:), vector_batch(icount,:,:),&
                  cond_rep, cond_ngwf_basis, cond_denskern,kchc, val_rep, &
                  val_ngwf_basis, val_denskern,hvkv,mdl,hfxstate,response_dens,&
                  sub_response_dens,fxc_fine, ground_state_dens,&
                  sub_ground_state_dens,y_contra_vec(icount,:,:), &
                  ground_state_hartree=ground_state_hartree)
          else
             call linear_response_RPA_operator(y_vec(icount,:,:), vector_batch(icount,:,:),&
                  cond_rep, cond_ngwf_basis, cond_denskern,kchc, val_rep, &
                  val_ngwf_basis, val_denskern,hvkv,mdl,hfxstate,response_dens,&
                  sub_response_dens,fxc_fine, ground_state_dens,&
                  sub_ground_state_dens,y_contra_vec(icount,:,:))
          endif
       enddo

       ! set new omega_prev and calculate new omega
       ! also set g_prev and grad_prev
       omega_prev=omega
       gamma_prev=gamma1
       do icount=1+pub_lr_tddft_num_conv_states, num_states
          do is=1, pub_num_spins
             call sparse_embed_copy(grad_prev_vec(icount,1,is),grad_vec(icount,1,is))
             call sparse_embed_copy(grad_prev_vec(icount,2,is),grad_vec(icount,2,is))
             call sparse_embed_copy(f_prev_vec(icount,1,is),f_vec(icount,1,is))
             call sparse_embed_copy(f_prev_vec(icount,2,is),f_vec(icount,2,is))
          enddo
       enddo

       call lr_tddft_RPA_calc_omega(omega,vector_batch,vector_batch_cov,y_vec,&
           cond_rep%overlap,val_rep%overlap, num_states)

       if(pub_lr_tddft_write_kernels) then
         do icount=1+pub_lr_tddft_num_conv_states, num_states
            call restart_RPA_kernel_write(vector_batch(icount,:,:),icount)
         enddo
       endif

       if (pub_on_root .and. pub_output_detail==BRIEF) then
          write(stdout,'(i16,f18.8,f18.8,f16.8,e15.4)') counter, omega, &
           omega-omega_prev, lambda, penalty_value
       endif

       ! OUTPUT convergence information
       if (pub_output_detail >= NORMAL) then
          if (pub_on_root) then
             write(stdout,'(/a)')'           ............................&
                  &............................'
             if (abs(omega-omega_prev)>tolerance) then
                write(stdout,'(a)')'           |     *** TDDFT optimisation&
                     & not converged ***         |'
             else
                write(stdout,'(a)')'           |     *** TDDFT optimisation&
                     & converged ***             |'
             endif
             write(stdout,'(a,f19.14,a)') &
                  '           |       LR-TDDFT energy  =  ',omega,'     &
                  &   |'
             write(stdout,'(a,f19.14,a)') &
                  '           |       Change in omega  =  ',omega-&
                  omega_prev,'        |'
             if (abs(omega-omega_prev)>tolerance) then
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

       if(reset_counter>pub_lr_tddft_reset_cg) then
         reset_counter=1
       endif
       reset_counter=reset_counter+1

       counter=counter+1

       ! include QC
       if (pub_on_root .and. pub_print_qc .and. .not. pub_lr_tddft_preopt) then
          call utils_qc_print('Omega',omega)
          call utils_qc_print('delta_Omega',omega-omega_prev)
       endif

       ! done. Close the outer loop
    enddo

    ! calculation is converged. Now we require a diagonalisation of the subspace
    allocate(identity_mat(num_states,num_states),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','identity_mat',ierr)
    allocate(subspace_mat(num_states,num_states),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_cg','subspace_mat',ierr)

    ! if some vectors are already converged, need to recalculate y for those vecs

    if(pub_lr_tddft_num_conv_states>0) then
       do icount=1, pub_lr_tddft_num_conv_states
          if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
             call linear_response_RPA_operator(y_vec(icount,:,:), vector_batch(icount,:,:),&
                  cond_rep, cond_ngwf_basis, cond_denskern,kchc, val_rep, &
                  val_ngwf_basis, val_denskern,hvkv,mdl,hfxstate,response_dens,&
                  sub_response_dens,fxc_fine, ground_state_dens,&
                  sub_ground_state_dens,y_contra_vec(icount,:,:), &
                  ground_state_hartree=ground_state_hartree)
          else
             call linear_response_RPA_operator(y_vec(icount,:,:), vector_batch(icount,:,:),&
                  cond_rep, cond_ngwf_basis, cond_denskern, kchc, val_rep, &
                  val_ngwf_basis, val_denskern,hvkv,mdl,hfxstate,response_dens,&
                  sub_response_dens,fxc_fine, ground_state_dens,&
                  sub_ground_state_dens,y_contra_vec(icount,:,:))
          endif
       enddo
    endif

    ! now fill identity mat and subspace mat:
    do icount=1, num_states
       do jcount=icount, num_states
          if(icount==jcount) then
            identity_mat(icount,jcount)=1.0_DP
          else
            identity_mat(icount,jcount)=0.0_DP
            identity_mat(jcount,icount)=0.0_DP
          endif

          call lr_tddft_RPA_inner_product(temp_val,vector_batch(icount,1,:),&
               y_vec(jcount,1,:),cond_rep%overlap,val_rep%overlap,.false.)
          call lr_tddft_RPA_inner_product(temp_val2,vector_batch(icount,2,:),&
               y_vec(jcount,2,:),cond_rep%overlap,val_rep%overlap,.false.)
          temp_total=(temp_val+temp_val2)/(2.0_DP)

          subspace_mat(icount,jcount)=temp_total
          subspace_mat(jcount,icount)=temp_total

       enddo
    enddo

    ! now diagonalise subspace:
    call linalg_dsygv_lt(subspace_mat,lr_tddft_evals,identity_mat,num_states)

    ! construct 'correct' eigenvectors
    do icount=1, num_states
       do is=1, pub_num_spins
          call sparse_embed_copy(grad_vec(icount,1,is), vector_batch(icount,1,is))
          call sparse_embed_copy(grad_vec(icount,2,is), vector_batch(icount,2,is))
       enddo
    enddo


    ! now construct correct response density matrices
    ! seems to be done correctly
    do icount=1, num_states
       do jcount=1, num_states
          do is=1, pub_num_spins
             if (jcount==1) then
                call sparse_embed_copy(vector_batch(icount,1,is), &
                     grad_vec(jcount,1,is))
                call sparse_embed_copy(vector_batch(icount,2,is), &
                     grad_vec(jcount,2,is))
                temp1=subspace_mat(jcount,icount)
                call sparse_embed_scale(vector_batch(icount,1,is), temp1)
                call sparse_embed_scale(vector_batch(icount,2,is), temp1)
             else
                temp1=subspace_mat(jcount,icount)
                call sparse_embed_axpy(vector_batch(icount,1,is), &
                     grad_vec(jcount,1,is), temp1)
                call sparse_embed_axpy(vector_batch(icount,2,is), &
                     grad_vec(jcount,2,is), temp1)
             endif
          enddo
       enddo
    enddo

    ! write the new eigenvectors to file if necessary
    if(pub_lr_tddft_write_kernels) then
      do icount=1, num_states
         call restart_RPA_kernel_write(vector_batch(icount,:,:),icount)
      enddo
    endif

    call timer_clock('lr_tddft_RPA_cg',2)

    ! done here. clean up and leave the routine
    call sparse_embed_destroy(temp_Svkv)
    call sparse_embed_destroy(PSv)
    call sparse_embed_destroy(L_vec(1))
    call sparse_embed_destroy(L_vec(2))
    do icount=1, num_states
       do is=1, pub_num_spins
          call sparse_embed_destroy(vector_batch_cov(icount,1,is))
          call sparse_embed_destroy(vector_batch_cov(icount,2,is))
          call sparse_embed_destroy(y_vec(icount,1,is))
          call sparse_embed_destroy(y_vec(icount,2,is))
          call sparse_embed_destroy(grad_vec(icount,1,is))
          call sparse_embed_destroy(grad_vec(icount,2,is))
          call sparse_embed_destroy(g_vec(icount,1,is))
          call sparse_embed_destroy(g_vec(icount,2,is))
          call sparse_embed_destroy(f_vec(icount,1,is))
          call sparse_embed_destroy(f_vec(icount,2,is))
          call sparse_embed_destroy(y_contra_vec(icount,1,is))
          call sparse_embed_destroy(y_contra_vec(icount,2,is))
          call sparse_embed_destroy(y_prev_vec(icount,1,is))
          call sparse_embed_destroy(y_prev_vec(icount,2,is))
          call sparse_embed_destroy(grad_prev_vec(icount,1,is))
          call sparse_embed_destroy(grad_prev_vec(icount,2,is))
          call sparse_embed_destroy(f_prev_vec(icount,1,is))
          call sparse_embed_destroy(f_prev_vec(icount,2,is))
       enddo
    enddo
    do is=1, pub_num_spins
       call sparse_embed_destroy(kcSc(is))
       call sparse_embed_destroy(Svkv(is))
    enddo
    call sparse_embed_destroy(LSvkv)
    call sparse_embed_destroy(temp_g(1))
    call sparse_embed_destroy(temp_g(2))

    deallocate(Svkv, stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','Svkv',ierr)
    deallocate(kcSc, stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','kcSc',ierr)
    deallocate(subspace_mat, stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','subspace_mat',ierr)
    deallocate(identity_mat, stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','identity_mat',ierr)
    if(pub_multigrid_hartree .and. .not. pub_lr_tddft_preopt) then
      deallocate(ground_state_hartree, stat=ierr)
      call utils_dealloc_check('lr_tddft_RPA_cg', 'ground_state_hartree', ierr)
    endif
    deallocate(L_vec,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','L_vec',ierr)
    deallocate(y_vec,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','y_vec',ierr)
    deallocate(grad_vec,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','grad_vec',ierr)
    deallocate(g_vec,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','g_vec',ierr)
    deallocate(f_vec, stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','f_vec',ierr)
    deallocate(y_prev_vec,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','y_prev_vec',ierr)
    deallocate(y_contra_vec,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','y_contra_vec',ierr)
    deallocate(grad_prev_vec,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','grad_prev_vec',ierr)
    deallocate(f_prev_vec,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','f_prev_vec',ierr)
    deallocate(temp_g,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','temp_g',ierr)
    deallocate(vector_batch_cov,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_cg','vector_batch_cov',ierr)

  end subroutine lr_tddft_RPA_cg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine lr_tddft_RPA_calc_omega_lambda(omega,lambda,&
     y_vec,operator_on_grad,x_trial,x_trial_cov,cond_overlap,&
     val_overlap, num_vecs)

     !===============================================================!
     ! Subroutine computes a predicted omega for a given line step   !
     ! and returns that number                                       !
     ! Modified for embedding by Joseph Prentice, July 2018          !
     !===============================================================!
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_copy, sparse_embed_transpose, sparse_embed_transpose_structure
    use rundat, only: pub_lr_tddft_num_conv_states, pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer :: num_vecs
    real(kind=DP), intent(inout) :: omega
    real(kind=DP), intent(in) :: lambda
    type(SPAM3_EMBED), intent(in) :: y_vec(num_vecs,2,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: operator_on_grad(num_vecs,2,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: x_trial(num_vecs,2,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: x_trial_cov(num_vecs,2,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_overlap
    type(SPAM3_EMBED), intent(in) :: val_overlap

    ! variables
    real(kind=DP) :: temp_val, temp_val2, &
       temp_normalise
    integer :: icount, ierr, is
    type(SPAM3_EMBED), allocatable, dimension(:,:) :: temp_cov

    omega=0.0_DP
    ! allocate data
    allocate(temp_cov(2,pub_num_spins),stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_calc_omega_lambda','temp_cov',ierr)
    do is=1, pub_num_spins
       call sparse_embed_create(temp_cov(1,is),y_vec(1,1,is))
       call sparse_embed_create(temp_cov(2,is),y_vec(1,2,is))
    enddo

    do icount=1+pub_lr_tddft_num_conv_states,num_vecs
       do is=1, pub_num_spins
          call sparse_embed_copy(temp_cov(1,is),y_vec(icount,1,is))
          call sparse_embed_copy(temp_cov(2,is),y_vec(icount,2,is))
          call sparse_embed_axpy(temp_cov(1,is),operator_on_grad(icount,1,is),lambda)
          call sparse_embed_axpy(temp_cov(2,is),operator_on_grad(icount,2,is),lambda)
       enddo
       call lr_tddft_RPA_inner_product(temp_val,x_trial(icount,1,:),&
            temp_cov(1,:),cond_overlap,val_overlap,.false.)
       call lr_tddft_RPA_inner_product(temp_val2,x_trial(icount,2,:),&
            temp_cov(2,:),cond_overlap,val_overlap,.false.)
       call lr_tddft_RPA_normalisation(temp_normalise,x_trial(icount,:,:),&
             x_trial_cov(icount,:,:), cond_overlap,val_overlap,.false.)

       omega=omega+(temp_val+temp_val2)/(abs(temp_normalise))

    enddo

    ! deallocate data
    do is=1, pub_num_spins
       call sparse_embed_destroy(temp_cov(1,is))
       call sparse_embed_destroy(temp_cov(2,is))
    enddo
    deallocate(temp_cov,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_calc_omega_lambda',&
       'temp_cov',ierr)

  end subroutine lr_tddft_RPA_calc_omega_lambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_RPA_calc_omega(omega,vector,vector_cov,y_vec,&
    cond_overlap,val_overlap, num_vecs)

    !===============================================================!
    ! Subroutine computes the Rayleigh quotient that is minimised   !
    ! in the conjugate gradient approach.                           !
    ! Modified for embedding by Joseph Prentice, July 2018          !
    !===============================================================!
    use sparse_embed, only: SPAM3_EMBED
    use rundat, only: pub_lr_tddft_num_conv_states,pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer :: num_vecs
    real(kind=DP), intent(inout) :: omega
    type(SPAM3_EMBED), intent(in) :: vector(num_vecs,2,pub_num_spins)
      ! contains (p,q).
    type(SPAM3_EMBED), intent(in) :: vector_cov(num_vecs,2,pub_num_spins)
      ! contains covariant version of (p,q)
    type(SPAM3_EMBED), intent(in) :: y_vec(num_vecs,2,pub_num_spins)
      ! contains [ (A-B)p, (A+B)q
    type(SPAM3_EMBED), intent(in) :: cond_overlap
    type(SPAM3_EMBED), intent(in) :: val_overlap

    ! variables
    real(kind=DP) :: temp_normalise, &
        temp_val1,temp_val2
    integer :: icount

    omega=0.0_DP
    do icount=1+pub_lr_tddft_num_conv_states, num_vecs
       call lr_tddft_RPA_inner_product(temp_val1,vector(icount,1,:),&
              y_vec(icount,1,:),cond_overlap,val_overlap,.false.)
       call lr_tddft_RPA_inner_product(temp_val2,vector(icount,2,:),&
              y_vec(icount,2,:),cond_overlap,val_overlap,.false.)
       call lr_tddft_RPA_normalisation(temp_normalise,vector(icount,:,:),&
                vector_cov(icount,:,:),cond_overlap,val_overlap,.false.)

       omega=omega+(temp_val1+temp_val2)/temp_normalise

    enddo

  end subroutine lr_tddft_RPA_calc_omega

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine lr_tddft_RPA_line_min(lambda, vector,vector_cov,y_vec,&
      grad_vec,cond_rep,cond_ngwf_basis,cond_denskern,kchc,val_rep,&
      val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,response_dens,&
      sub_response_dens,fxc,L_vec,num_vecs,ground_state_dens,&
      sub_ground_state_dens,ground_state_hartree)

    !==============================================================!
    ! Subroutine computes the ideal conjugate gradient line step   !
    ! length for the RPA eigenvalue problem.                       !
    ! The subroutine also updates the vector batch accordingly,    !
    ! using the ideal line step length and the gradient.           !
    ! Modified for embedding by Joseph Prentice, July 2018         !
    !==============================================================!

    use comms, only: pub_on_root
    use function_basis, only: FUNC_BASIS
    use constants, only: stdout, VERBOSE, NORMAL
    use hf_exchange, only: HFX_STATE
    use linear_response, only: linear_response_RPA_operator
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_lr_tddft_num_conv_states, pub_num_spins,&
         pub_output_detail
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_copy, sparse_embed_transpose, sparse_embed_transpose_structure
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer :: num_vecs
    real(kind=DP), intent(inout) :: lambda
    type(SPAM3_EMBED), intent(inout) :: vector(num_vecs,2,pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: vector_cov(num_vecs,2,pub_num_spins)
      ! covariant version of vector
    type(SPAM3_EMBED), intent(in) :: y_vec(num_vecs,2,pub_num_spins) ! y_batch
      ! is the result of vector_batch acting on the TDDFT operator
    type(SPAM3_EMBED), intent(inout) :: grad_vec(num_vecs,2,pub_num_spins) ! grad_batch
      ! is the batch of conjugate gradient search directions
    type(SPAM3_EMBED), intent(inout) :: L_vec(2) ! Auxillary
      ! response kernel. These terms are used to construct the correct
      ! kernel in case a kernel truncation is used.
    type(NGWF_REP), intent(inout) :: cond_rep
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis(:)
    type(SPAM3_EMBED), intent(in) :: kchc(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_denskern(pub_num_spins)
    type(NGWF_REP), intent(inout) :: val_rep
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type(SPAM3_EMBED), intent(in) :: hvkv(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: val_denskern(pub_num_spins)
    type(MODEL), intent(in) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    real(kind=DP), intent(inout) :: response_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_response_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)
    real(kind=DP), intent(in) :: fxc(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,2,mdl%nsub)
    real(kind=DP), intent(inout) :: ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)
    real(kind=DP), optional,intent(in) :: ground_state_hartree(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, &
         mdl%nsub, mdl%nsub)


    ! local variables
    real(kind=DP) :: temp_val, temp_val2, &
          temp_val3,temp_val4
    real(kind=DP) :: a, b, c, discriminant,&
           lambda1,lambda2
    integer :: ierr, icount, is
    ! values that can be reused in line min
    real(kind=DP) :: xx,xy,xD,DopD,DD,Dy
    type(SPAM3_EMBED) :: PSv
    type(SPAM3_EMBED), allocatable, dimension(:,:,:) :: operator_on_grad
    type(SPAM3_EMBED), allocatable, dimension(:,:,:) :: operator_on_grad_contra
    type(SPAM3_EMBED), allocatable, dimension(:,:,:) :: grad_vec_cov
    type(SPAM3_EMBED), allocatable, dimension(:,:,:) :: trial_x,trial_x_cov
    real(kind=DP) :: omega1,omega2

    call timer_clock('lr_tddft_RPA_line_min',1)

    ! allocate storage space
    allocate(trial_x(num_vecs,2,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_line_min','trial_x',ierr)
    allocate(trial_x_cov(num_vecs,2,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_line_min','trial_x_cov',ierr)
    allocate(operator_on_grad(num_vecs,2,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_line_min','operator_on_grad',ierr)
    allocate(operator_on_grad_contra(num_vecs,2,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_line_min','operator_on_grad_contra',ierr)
    allocate(grad_vec_cov(num_vecs,2,pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_line_min','grad_vec_cov',ierr)

    ! create sparse matrices
    do icount=1,num_vecs
       do is=1,pub_num_spins ! loop over spin variables
          call sparse_embed_create(grad_vec_cov(icount,1,is),y_vec(1,1,1))
          call sparse_embed_create(grad_vec_cov(icount,2,is),y_vec(1,2,1))
          call sparse_embed_create(operator_on_grad_contra(icount,1,is),&
              vector(1,1,1))
          call sparse_embed_create(operator_on_grad_contra(icount,2,is),&
              vector(1,2,1))
          call sparse_embed_create(operator_on_grad(icount,1,is),&
              y_vec(1,1,1))
          call sparse_embed_create(operator_on_grad(icount,2,is),&
              y_vec(1,2,1))
          call sparse_embed_create(trial_x(icount,1,is),vector(1,1,1))
          call sparse_embed_create(trial_x(icount,2,is),vector(1,2,1))
          call sparse_embed_create(trial_x_cov(icount,1,is),vector(1,1,1))
          call sparse_embed_create(trial_x_cov(icount,2,is),vector(1,2,1))
       enddo
    enddo
    call sparse_embed_create(PSv,vector(1,1,1),val_rep%overlap)

    ! now calculate the correct action.
    do icount=1+pub_lr_tddft_num_conv_states, num_vecs
       if(present(ground_state_hartree)) then
          call linear_response_RPA_operator(operator_on_grad(icount,:,:),&
             grad_vec(icount,:,:),cond_rep,cond_ngwf_basis,cond_denskern,&
             kchc,val_rep,val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,&
             response_dens,sub_response_dens,fxc,ground_state_dens,&
             sub_ground_state_dens,operator_on_grad_contra(icount,:,:),&
             ground_state_hartree=ground_state_hartree)
       else
          call linear_response_RPA_operator(operator_on_grad(icount,:,:),&
             grad_vec(icount,:,:),cond_rep,cond_ngwf_basis,cond_denskern,&
             kchc,val_rep,val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,&
             response_dens,sub_response_dens,fxc,ground_state_dens,&
             sub_ground_state_dens,operator_on_grad_contra(icount,:,:))
       endif

       do is=1, pub_num_spins
          call sparse_embed_product(PSv,grad_vec(icount,1,is),val_rep%overlap)
          call sparse_embed_product(grad_vec_cov(icount,1,is),cond_rep%overlap,PSv)
          call sparse_embed_product(PSv,grad_vec(icount,2,is),val_rep%overlap)
          call sparse_embed_product(grad_vec_cov(icount,2,is),cond_rep%overlap,PSv)
       enddo
    enddo

    a=0.0_DP
    b=0.0_DP
    c=0.0_DP
    ! now calculate appropriate a, b and c coefficients. Start with a
    do icount=1+pub_lr_tddft_num_conv_states, num_vecs

       call lr_tddft_RPA_inner_product(temp_val, grad_vec(icount,1,:),&
             y_vec(icount,1,:),cond_rep%overlap,val_rep%overlap,.false.)
       call lr_tddft_RPA_inner_product(temp_val2,grad_vec(icount,2,:),&
             y_vec(icount,2,:),cond_rep%overlap,val_rep%overlap,.false.)
       call lr_tddft_RPA_inner_product(temp_val3,vector(icount,1,:),&
             vector_cov(icount,2,:),cond_rep%overlap,val_rep%overlap,.false.)

       xx=temp_val3
       Dy=temp_val+temp_val2

       a=a+2.0_DP*xx*Dy


       call lr_tddft_RPA_inner_product(temp_val,grad_vec(icount,1,:),&
              vector_cov(icount,2,:),cond_rep%overlap,val_rep%overlap,.false.)
       call lr_tddft_RPA_inner_product(temp_val2,grad_vec(icount,2,:),&
              vector_cov(icount,1,:),cond_rep%overlap,val_rep%overlap,.false.)

       call lr_tddft_RPA_inner_product(temp_val3,vector(icount,1,:),&
             y_vec(icount,1,:),cond_rep%overlap,val_rep%overlap,.false.)
       call lr_tddft_RPA_inner_product(temp_val4,vector(icount,2,:),&
             y_vec(icount,2,:),cond_rep%overlap,val_rep%overlap,.false.)

       xy=temp_val3+temp_val4
       xD=temp_val+temp_val2

       a=a-1.0_DP*xD*xy

       ! now b component
       call lr_tddft_RPA_inner_product(temp_val, grad_vec(icount,1,:),&
             operator_on_grad(icount,1,:),cond_rep%overlap,val_rep%overlap,.false.)
       call lr_tddft_RPA_inner_product(temp_val2,grad_vec(icount,2,:),&
             operator_on_grad(icount,2,:),cond_rep%overlap,val_rep%overlap,.false.)

       DopD=temp_val+temp_val2

       b=b+2.0_DP*(DopD)*(xx)

       call lr_tddft_RPA_inner_product(temp_val3,grad_vec(icount,1,:),&
             grad_vec_cov(icount,2,:),cond_rep%overlap,val_rep%overlap,.false.)

       DD=temp_val3

       b=b-2.0_DP*(xy)*(DD)
       ! finally c component

       c=c+1.0_DP*(DopD)*(xD)

       c=c-2.0_DP*(DD)*(Dy)
    enddo

    discriminant=b*b-4.0_DP*a*c

    if(abs(discriminant)>1.0e-20_DP .and. discriminant>0.0_DP) then
      if(b<0.0_DP) then
         lambda1=(-b+sqrt(discriminant))/(2.0_DP*c)
      else
         lambda1=(-b-sqrt(discriminant))/(2.0_DP*c)
      endif

      lambda2=a/(c*lambda1)

      ! check which of the 2 solutions is the ideal one:
      do icount=1+pub_lr_tddft_num_conv_states, num_vecs
         do is=1, pub_num_spins ! loop over spin variables
            call sparse_embed_copy(trial_x(icount,1,is),vector(icount,1,is))
            call sparse_embed_axpy(trial_x(icount,1,is),grad_vec(icount,1,is),lambda1)
            call sparse_embed_copy(trial_x(icount,2,is),vector(icount,2,is))
            call sparse_embed_axpy(trial_x(icount,2,is),grad_vec(icount,2,is),lambda1)
            call sparse_embed_copy(trial_x_cov(icount,1,is),vector_cov(icount,1,is))
            call sparse_embed_axpy(trial_x_cov(icount,1,is),grad_vec_cov(icount,1,is),&
                 lambda1)
            call sparse_embed_copy(trial_x_cov(icount,2,is),vector_cov(icount,2,is))
            call sparse_embed_axpy(trial_x_cov(icount,2,is),grad_vec_cov(icount,2,is),&
                 lambda1)
         enddo
      enddo

      call lr_tddft_RPA_calc_omega_lambda(omega1,lambda1,&
          y_vec,operator_on_grad,trial_x,trial_x_cov,cond_rep%overlap,&
          val_rep%overlap, num_vecs)

      do icount=1+pub_lr_tddft_num_conv_states, num_vecs
         do is=1, pub_num_spins ! loop over spin variables
            call sparse_embed_copy(trial_x(icount,1,is),vector(icount,1,is))
            call sparse_embed_axpy(trial_x(icount,1,is),grad_vec(icount,1,is),lambda2)
            call sparse_embed_copy(trial_x(icount,2,is),vector(icount,2,is))
            call sparse_embed_axpy(trial_x(icount,2,is),grad_vec(icount,2,is),lambda2)
            call sparse_embed_copy(trial_x_cov(icount,1,is),vector_cov(icount,1,is))
            call sparse_embed_axpy(trial_x_cov(icount,1,is),grad_vec_cov(icount,1,is),&
                 lambda2)
            call sparse_embed_copy(trial_x_cov(icount,2,is),vector_cov(icount,2,is))
            call sparse_embed_axpy(trial_x_cov(icount,2,is),grad_vec_cov(icount,2,is),&
                 lambda2)
         enddo
      enddo

      call lr_tddft_RPA_calc_omega_lambda(omega2,lambda2,&
          y_vec,operator_on_grad,trial_x,trial_x_cov,cond_rep%overlap,&
          val_rep%overlap, num_vecs)

      if(pub_on_root .and. pub_output_detail==VERBOSE) write(stdout,'(a,E12.4,a)') &
         '******** Line-min: omega1=',omega1,'*********'
      if(pub_on_root .and. pub_output_detail==VERBOSE) write(stdout,'(a,E12.4,a)') &
        '******** Line-min: omega2=',omega2,'*********'

      if(omega1<omega2) then
        lambda=lambda1
      else
        omega1=omega2
        lambda=lambda2
      endif

    else
       if (pub_on_root) write(stdout, '(a)') 'WARNING, discriminant in&
            & lr_tddft_line_min is negative or close to zero'
       if (pub_on_root) write(stdout, *) discriminant
       lambda=-b/(2.0_DP*c)
    endif

    ! successfully calculated the ideal line step.
    ! update vector batch

    if(pub_on_root .and. pub_output_detail>=VERBOSE) then
       write(stdout,'(a,E12.4,a)') &
       '******** Line-min: c coefficient=',c,'*********'
       write(stdout,'(a,E12.4,a)') &
       '******** Line-min: b coefficient=',b,'*********'
       write(stdout,'(a,E12.4,a)') &
       '******** Line-min: a coefficient=',a,'*********'
    end if
    if(pub_on_root .and. pub_output_detail>=NORMAL) write(stdout,'(a,E12.4,a)') &
       '******** Line-min: lambda coefficient=',lambda,'*********'

    if(pub_on_root .and. pub_output_detail>=NORMAL) write(stdout,'(a,f14.10,a)') &
       '******** Line-min: Predicted omega=',omega1,'*********'

    do icount=1+pub_lr_tddft_num_conv_states, num_vecs
       do is=1, pub_num_spins
          call sparse_embed_axpy(vector(icount,1,is),grad_vec(icount,1,is),lambda)
          call sparse_embed_axpy(vector(icount,2,is),grad_vec(icount,2,is),lambda)
          call sparse_embed_axpy(vector_cov(icount,1,is),grad_vec_cov(icount,1,is),&
               lambda)
          call sparse_embed_axpy(vector_cov(icount,2,is),grad_vec_cov(icount,2,is),&
               lambda)
       enddo
    enddo

    ! deallocate data
    ! destroy
    do icount=1, num_vecs
       do is=1, pub_num_spins
          call sparse_embed_destroy(grad_vec_cov(icount,1,is))
          call sparse_embed_destroy(grad_vec_cov(icount,2,is))
          call sparse_embed_destroy(operator_on_grad_contra(icount,1,is))
          call sparse_embed_destroy(operator_on_grad_contra(icount,2,is))
          call sparse_embed_destroy(operator_on_grad(icount,1,is))
          call sparse_embed_destroy(operator_on_grad(icount,2,is))
          call sparse_embed_destroy(trial_x(icount,1,is))
          call sparse_embed_destroy(trial_x(icount,2,is))
          call sparse_embed_destroy(trial_x_cov(icount,1,is))
          call sparse_embed_destroy(trial_x_cov(icount,2,is))
       enddo
    enddo
    call sparse_embed_destroy(PSv)

    deallocate(grad_vec_cov,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_line_min','grad_vec_cov',ierr)
    deallocate(operator_on_grad,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_line_min','operator_on_grad',ierr)
    deallocate(operator_on_grad_contra,stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_line_min','operator_on_grad_contra',ierr)
    deallocate(trial_x_cov, stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_line_min','trial_x_cov',ierr)
    deallocate(trial_x, stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_line_min','trial_x',ierr)

    call timer_clock('lr_tddft_RPA_line_min',2)

  end subroutine lr_tddft_RPA_line_min

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_RPA_orthogonalise(vector_batch,orthogonal_vectors,&
          orthogonal_vectors_cov,num_vecs,cond_overlap,val_overlap)

     !=============================================================!
     ! Subroutine takes a contravariant vector batch and makes it  !
     ! orthogonal to a contravariant vector batch called           !
     ! orthogonal_vectors.                                         !
     ! In case the calculation is not spin-degenerate, we make     !
     ! sure that the orthogonalisation holds for each spin channel !
     ! separately.                                                 !
     ! Modified for embedding by Joseph Prentice, July 2018        !
     !=============================================================!

    use sparse_embed, only: SPAM3_EMBED, sparse_embed_axpy
    use rundat, only: pub_lr_tddft_num_conv_states,pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: num_vecs
    type(SPAM3_EMBED), intent(inout) :: vector_batch(num_vecs,2,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: orthogonal_vectors(num_vecs,2,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: orthogonal_vectors_cov(num_vecs,2,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_overlap
    type(SPAM3_EMBED), intent(in) :: val_overlap

    ! Local Variables
    integer :: jcount, icount, is, ierr
    real(kind=DP) :: temp1, temp2, temp_scaling
    type(SPAM3_EMBED) :: PSv

    ! loop over all unconverged states of vector batch
    do icount=1+pub_lr_tddft_num_conv_states, num_vecs
       ! loop over all states of orthogonal vectors

       do jcount=1, num_vecs
          call lr_tddft_RPA_normalisation(temp1,vector_batch(icount,:,:),&
             orthogonal_vectors_cov(jcount,:,:),cond_overlap,val_overlap,.false.)

          call lr_tddft_RPA_normalisation(temp2,orthogonal_vectors(jcount,:,:),&
             orthogonal_vectors_cov(jcount,:,:),cond_overlap,val_overlap,.false.)

             temp_scaling=temp1/temp2

          do is=1,pub_num_spins
             call sparse_embed_axpy(vector_batch(icount,1,is),orthogonal_vectors(jcount,1,is),&
                   -1.0_DP*temp_scaling)
             call sparse_embed_axpy(vector_batch(icount,2,is),orthogonal_vectors(jcount,2,is),&
                   -1.0_DP*temp_scaling)
          enddo

       enddo
    enddo

  end subroutine lr_tddft_RPA_orthogonalise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_RPA_gram_schmidt(vector_batch,vector_batch_cov,&
       num_vecs,cond_overlap,val_overlap)

    !============================================================!
    ! Subroutine performs an orthonormalisation procedure of the !
    ! p and q vectors by enforcing that Tr[pq]=1                 !
    ! Note that if p*q is negative, this routine also reverts the!
    ! sign of the q vector (this is an arbitrary choice)         !
    ! In case that the system is spin-polarised, the Gram-       !
    ! Schmidt orthonormalisation guarantees that all vectors are !
    ! orthogonal to each other in each spin channel.             !
    ! Modified for embedding by Joseph Prentice, July 2018       !
    !============================================================!

    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_copy, sparse_embed_transpose, sparse_embed_transpose_structure
    use rundat, only: pub_lr_tddft_num_conv_states,pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: num_vecs
    type(SPAM3_EMBED), intent(inout) :: vector_batch(num_vecs,2,pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: vector_batch_cov(num_vecs,2,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_overlap
    type(SPAM3_EMBED), intent(in) :: val_overlap

    ! Local Variables
    integer :: count_i
    integer :: count_j, ierr, is
    real(kind=DP) :: temp1, temp2, temp_scaling

    ! check if RPA states are already converged. In that case
    ! orthogonalise all unconverged states against them.
    if(pub_lr_tddft_num_conv_states>0) then
       do count_j=1+pub_lr_tddft_num_conv_states, num_vecs

          do count_i=1, pub_lr_tddft_num_conv_states

             call lr_tddft_RPA_normalisation(temp1,&
                  vector_batch(count_i,:,:), vector_batch_cov(count_j,:,:),&
                  cond_overlap,val_overlap,.false.)

             call lr_tddft_RPA_normalisation(temp2,&
                  vector_batch(count_i,:,:),vector_batch_cov(count_i,:,:),&
                  cond_overlap,val_overlap,.false.)

             temp_scaling=temp1/temp2

             do is=1, pub_num_spins
                call sparse_embed_axpy(vector_batch(count_j,1,is),&
                     vector_batch(count_i,1,is),-1.0_DP*temp_scaling)
                call sparse_embed_axpy(vector_batch(count_j,2,is),&
                     vector_batch(count_i,1,is),-1.0_DP*temp_scaling)
                call sparse_embed_axpy(vector_batch_cov(count_j,1,is),&
                     vector_batch_cov(count_i,1,is),-1.0_DP*temp_scaling)
                call sparse_embed_axpy(vector_batch_cov(count_j,2,is),&
                     vector_batch_cov(count_i,1,is),-1.0_DP*temp_scaling)
             enddo
          enddo
       enddo
    endif

    ! start regular Gram schmidt only for unconverged
    ! vectors
    do count_i=1+pub_lr_tddft_num_conv_states, num_vecs
       ! normalize vec. Make sure that the normalisation reduces the
       ! result of the inner product to 2!
       call lr_tddft_RPA_normalisation(temp1,vector_batch(count_i,:,:),&
          vector_batch_cov(count_i,:,:),cond_overlap, val_overlap,.false.)

       ! check for negative p*q norm in each spin channel
       if(temp1<0.0_DP) then
          do is=1, pub_num_spins
              call sparse_embed_scale(vector_batch(count_i,2,is),-1.0_DP)
              call sparse_embed_scale(vector_batch_cov(count_i,2,is), -1.0_DP)
          enddo
          call lr_tddft_RPA_normalisation(temp1,vector_batch(count_i,:,:),&
             vector_batch_cov(count_i,:,:),cond_overlap,val_overlap,.false.)
       endif

       ! normalise
       do is=1, pub_num_spins
          call sparse_embed_scale(vector_batch(count_i,1,is),sqrt(2.0_DP)/sqrt(abs(temp1)))
          call sparse_embed_scale(vector_batch(count_i,2,is),sqrt(2.0_DP)/sqrt(abs(temp1)))
          call sparse_embed_scale(vector_batch_cov(count_i,1,is),sqrt(2.0_DP)/sqrt(abs(temp1)))
          call sparse_embed_scale(vector_batch_cov(count_i,2,is),sqrt(2.0_DP)/sqrt(abs(temp1)))
       enddo

       do count_j=count_i+1, num_vecs

          call lr_tddft_RPA_normalisation(temp1,vector_batch(count_j,:,:),&
              vector_batch_cov(count_i,:,:),cond_overlap,val_overlap,.false.)
          call lr_tddft_RPA_normalisation(temp2,vector_batch(count_i,:,:),&
              vector_batch_cov(count_i,:,:),cond_overlap,val_overlap,.false.)

          temp_scaling=temp1/temp2

          do is=1, pub_num_spins
             call sparse_embed_axpy(vector_batch(count_j,1,is),vector_batch(count_i,1,is),&
                  -1.0_DP*temp_scaling)
             call sparse_embed_axpy(vector_batch(count_j,2,is),vector_batch(count_i,2,is),&
                 -1.0_DP*temp_scaling)
             call sparse_embed_axpy(vector_batch_cov(count_j,1,is),vector_batch_cov(count_i,1,is),&
                  -1.0_DP*temp_scaling)
             call sparse_embed_axpy(vector_batch_cov(count_j,2,is),vector_batch_cov(count_i,2,is),&
                 -1.0_DP*temp_scaling)
          enddo
       enddo

       ! normalize vec again. Make sure that the normalisation reduces the
       ! result of the inner product to 2! Only necessary if calculation is
       ! spin polarised
       if(pub_num_spins==2) then
          call lr_tddft_RPA_normalisation(temp1,vector_batch(count_i,:,:),&
             vector_batch_cov(count_i,:,:),cond_overlap, val_overlap,.false.)

          ! check for negative p*q norm in each spin channel
          if(temp1<0.0_DP) then
             do is=1, pub_num_spins
                call sparse_embed_scale(vector_batch(count_i,2,is),-1.0_DP)
                call sparse_embed_scale(vector_batch_cov(count_i,2,is),-1.0_DP)
             enddo
             call lr_tddft_RPA_normalisation(temp1,vector_batch(count_i,:,:),&
                vector_batch_cov(count_i,:,:),cond_overlap,val_overlap,.false.)
          endif

          ! normalise
          do is=1, pub_num_spins
             call sparse_embed_scale(vector_batch(count_i,1,is),&
                    sqrt(2.0_DP)/sqrt(abs(temp1)))
             call sparse_embed_scale(vector_batch(count_i,2,is),&
                    sqrt(2.0_DP)/sqrt(abs(temp1)))
             call sparse_embed_scale(vector_batch_cov(count_i,1,is),&
                    sqrt(2.0_DP)/sqrt(abs(temp1)))
             call sparse_embed_scale(vector_batch_cov(count_i,2,is),&
                    sqrt(2.0_DP)/sqrt(abs(temp1)))
          enddo
       endif
    enddo


  end subroutine lr_tddft_RPA_gram_schmidt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_RPA_normalisation(value, bra_vec_contra,ket_vec,&
   cond_overlap, val_overlap,is_contra)

    !==============================================================!
    ! Subroutine computes the effective normalisation between two  !
    ! 2-component RPA vectors. In symbolic notation, the           !
    ! Normalisation factor is given by:                            !
    !         (x1 y1)(0 1)(x2)                                     !
    !                (1 0)(y2)                                     !
    ! Modified for embedding by Joseph Prentice, July 2018         !
    !==============================================================!

    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_copy, sparse_embed_transpose, sparse_embed_transpose_structure
    use rundat, only: pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: value
    type(SPAM3_EMBED), intent(in) :: bra_vec_contra(2,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: ket_vec(2,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_overlap
    type(SPAM3_EMBED), intent(in) :: val_overlap
    logical, optional, intent(in) :: is_contra

    ! local variables
    logical :: loc_is_contra
    real(kind=DP) :: temp1, temp2

    loc_is_contra=.false.
    if(present(is_contra)) then
      loc_is_contra=is_contra
    endif

    call lr_tddft_RPA_inner_product(temp1,bra_vec_contra(1,:),&
          ket_vec(2,:),cond_overlap,val_overlap,loc_is_contra)

    call lr_tddft_RPA_inner_product(temp2,bra_vec_contra(2,:),&
          ket_vec(1,:),cond_overlap,val_overlap,loc_is_contra)

    value=temp1+temp2

  end subroutine lr_tddft_RPA_normalisation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_RPA_inner_product(value1, bra_vec_contra,&
       ket_vec,cond_overlap,val_overlap,is_contra)

    !==============================================================!
    ! This routine calculates the inner product between two vectors!
    ! the bra_vec is always a contravariant quantity, while the ket!
    ! vec is either contravariant or covariant. If it is contra-   !
    ! variant, it is transformed by multiplying the overlap        !
    ! matrices from the left and right.                            !
    ! Modified for embedding by Joseph Prentice, July 2018         !
    !==============================================================!

    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_copy, sparse_embed_transpose, &
         sparse_embed_transpose_structure, sparse_embed_trace
    use rundat, only: pub_num_spins

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: value1
    type(SPAM3_EMBED), intent(in) :: bra_vec_contra(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: ket_vec(pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: cond_overlap
    type(SPAM3_EMBED), intent(in) :: val_overlap
    logical, optional, intent(in) :: is_contra

    ! local variables
    type(SPAM3_EMBED) :: temp_cov, xSv
    type(SPAM3_EMBED) :: temp_trans
    real(kind=DP) :: trace
    logical :: loc_is_contra
    integer :: is

    loc_is_contra=.false.
    if(present(is_contra)) then
      loc_is_contra=is_contra
    endif

    ! allocate sparse matrices
    call sparse_embed_create(xSv,bra_vec_contra(1),val_overlap)
    call sparse_embed_create(temp_cov,cond_overlap,xSv)
    call sparse_embed_transpose_structure(temp_trans%structure,&
         bra_vec_contra(1))
    call sparse_embed_create(temp_trans, iscmplx=bra_vec_contra(1)%p%iscmplx)

    value1=0.0_DP
    ! sum over spins
    do is=1, pub_num_spins
       if(loc_is_contra) then
         ! create covariant quantity out of contravariant ket_vec
         call sparse_embed_product(xSv,ket_vec(is),val_overlap)
         call sparse_embed_product(temp_cov,cond_overlap,xSv)
       else
         call sparse_embed_copy(temp_cov,ket_vec(is))
       endif

       ! calculate the inner product.
       call sparse_embed_transpose(temp_trans,bra_vec_contra(is))
       call sparse_embed_trace(trace,temp_trans,temp_cov)
       value1=value1+trace
    enddo

    ! deallocate storage space
    call sparse_embed_destroy(xSv)
    call sparse_embed_destroy(temp_cov)
    call sparse_embed_destroy(temp_trans)

  end subroutine lr_tddft_RPA_inner_product

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_tddft_RPA_analysis(response_kernel, cond_rep, &
       cond_ngwf_basis, val_rep, val_ngwf_basis, val_evecs,&
       cond_evecs, num_states, val_ngwf_num, cond_ngwf_num)

    !========================================================================!
    ! Subroutine calculates some properties for individual excitations (ie   !
    ! Which transitions make up which response vectors. Also, and very im -  !
    ! portantly: How much percentage of 'forbidden' transitions is contained !
    ! per excitation energy?                                                 !
    ! Modified for embedding by Joseph Prentice, July 2018                   !
    !========================================================================!

    use comms, only: pub_on_root
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
         dense_put_element,dense_get_element
    use function_basis, only: FUNC_BASIS
    use constants, only: stdout
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_lr_tddft_joint_set,pub_lr_tddft_HOMO_num,&
         pub_lr_tddft_LUMO_num, pub_num_kpoints, PUB_1K, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_transpose, sparse_embed_transpose_structure, sparse_embed_copy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: num_states
    type(SPAM3_EMBED), intent(inout) :: response_kernel(num_states,2,pub_num_spins)
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis(:)
    type(NGWF_REP), intent(in) :: cond_rep
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type(NGWF_REP), intent(in) :: val_rep
    type(DEM), intent(inout) :: val_evecs(pub_num_spins)
    type(DEM), intent(inout) :: cond_evecs(pub_num_spins)
    integer, intent(in) :: val_ngwf_num, cond_ngwf_num

    ! Local Variables
    type(SPAM3_EMBED):: ksv, scksv, val_sparse, cond_sparse, cond_trans,&
        P1_evec,evec_P1_evec, x_vec, y_vec
    type(DEM) :: kernel_dens_x, kernel_dens_y
    integer :: icount, jcount, omega_count
    real(kind=DP) :: temp_sum_x, temp_sum_y,temp_check_sum_x,temp_check_sum_y,&
        check_sum_printed_x, check_sum_printed_y
    integer, allocatable, dimension(:) :: cond_limit, val_limit, cond_start
    integer :: is, ierr

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine lr_tddft_RPA_analysis not ready yet for more&
         & than one k-point.')

    call timer_clock('lr_tddft_RPA_analysis',1)

    ! allocate data structures
    allocate(cond_limit(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_analysis','cond_limit',ierr)
    allocate(val_limit(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_analysis','val_limit',ierr)
    allocate(cond_start(pub_num_spins), stat=ierr)
    call utils_alloc_check('lr_tddft_RPA_analysis','cond_start',ierr)
    call sparse_embed_create(ksv, response_kernel(1,1,1), val_rep%overlap)
    call sparse_embed_create(scksv,cond_rep%overlap,ksv)
    call sparse_embed_create(cond_sparse,cond_rep%inv_overlap)
    call sparse_embed_create(val_sparse,val_rep%inv_overlap)
    call sparse_embed_create(cond_trans, cond_rep%inv_overlap)
    call sparse_embed_create(P1_evec, scksv,val_sparse)
    call sparse_embed_create(evec_P1_evec, cond_trans,P1_evec)
    call dense_create(kernel_dens_x,cond_ngwf_num,val_ngwf_num)
    call dense_create(kernel_dens_y,cond_ngwf_num,val_ngwf_num)
    call sparse_embed_create(x_vec,response_kernel(1,1,1))
    call sparse_embed_create(y_vec,response_kernel(1,1,1))

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
    enddo ! end spin loop

    do omega_count=1, num_states
       ! start spin loop
       if (pub_on_root) write(stdout, '(a,i4,a)') '******Transition: ',&
            omega_count, '*****************'

       temp_check_sum_x=0.0_DP
       temp_check_sum_y=0.0_DP
       check_sum_printed_x=0.0_DP
       check_sum_printed_y=0.0_DP

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
          ! first create effective x, then effective y matrix:
          call dense_convert(val_sparse,val_evecs(is))
          call dense_convert(cond_sparse,cond_evecs(is))
          call sparse_embed_transpose(cond_trans,cond_sparse)
          call sparse_embed_copy(x_vec,response_kernel(omega_count,1,is))
          call sparse_embed_axpy(x_vec,response_kernel(omega_count,2,is),1.0_DP)
          call sparse_embed_scale(x_vec,0.5_DP)
          call sparse_embed_copy(y_vec,response_kernel(omega_count,1,is))
          call sparse_embed_axpy(y_vec,response_kernel(omega_count,2,is),-1.0_DP)
          call sparse_embed_scale(y_vec,0.5_DP)
          call sparse_embed_product(ksv,x_vec,val_rep%overlap)
          call sparse_embed_product(scksv,cond_rep%overlap,ksv)
          call sparse_embed_product(P1_evec,scksv,val_sparse)
          call sparse_embed_product(evec_P1_evec, cond_trans, P1_evec)
          call dense_convert(kernel_dens_x,evec_P1_evec)
          call sparse_embed_product(ksv,y_vec,val_rep%overlap)
          call sparse_embed_product(scksv,cond_rep%overlap,ksv)
          call dense_convert(val_sparse,val_evecs(is))
          call sparse_embed_product(P1_evec,scksv,val_sparse)
          call sparse_embed_product(evec_P1_evec, cond_trans, P1_evec)
          call dense_convert(kernel_dens_y,evec_P1_evec)

          ! loop over allowed valence states
          ! we limit ourselves to the most important states
          do icount = val_rep%n_occ(is,PUB_1K)-val_limit(is)+1, &
                      val_rep%n_occ(is,PUB_1K)
             do jcount=cond_start(is), cond_limit(is)
                ! get overlap between vectors
                temp_sum_x=0.0_DP
                temp_sum_y=0.0_DP
                call dense_get_element(temp_sum_x, kernel_dens_x, jcount, icount)
                call dense_get_element(temp_sum_y, kernel_dens_y, jcount, icount)

                ! square vector element
                temp_sum_x=temp_sum_x*temp_sum_x
                temp_sum_y=temp_sum_y*temp_sum_y

                if (temp_sum_x>0.001_DP) then
                   if (pub_on_root) write(stdout,'(a,i4,a,i4,a,f12.8,a)') &
                        'KS transition ', icount, '  ----> ',&
                        jcount+val_rep%n_occ(is,PUB_1K)+1-cond_start(is), &
                        '  = ',temp_sum_x,'   X'
                   check_sum_printed_x=check_sum_printed_x+temp_sum_x
                endif
                if (temp_sum_y>0.001_DP) then
                   if (pub_on_root) write(stdout,'(a,i4,a,i4,a,f12.8,a)') &
                        'KS transition ', icount, '  ----> ',&
                        jcount+val_rep%n_occ(is,PUB_1K)+1-cond_start(is), &
                        '  = ',temp_sum_y,'   Y'
                   check_sum_printed_y=check_sum_printed_y+temp_sum_y
                endif

                temp_check_sum_x=temp_check_sum_x+temp_sum_x
                temp_check_sum_y=temp_check_sum_y+temp_sum_y
             enddo
          enddo
       enddo ! end spin loop

       if (pub_on_root) write(stdout, '(a, i4, a, f15.12)') &
            'Total contribution of transition X',omega_count,&
            '  printed above: ', check_sum_printed_x


       if (pub_on_root) write(stdout, '(a, i4, a, f15.12)') &
            'Total contribution of transition X',omega_count,&
            '  in allowed space: ', temp_check_sum_x

       if (pub_on_root) write(stdout, '(a, i4, a, f15.12)') &
            'Total contribution of transition Y',omega_count,&
            '  printed above: ', check_sum_printed_y


       if (pub_on_root) write(stdout, '(a, i4, a, f15.12)') &
            'Total contribution of transition Y',omega_count,&
            '  in allowed space: ', temp_check_sum_y

    enddo

    ! deallocate all data structures:
    call sparse_embed_destroy(ksv)
    call sparse_embed_destroy(scksv)
    call sparse_embed_destroy(cond_sparse)
    call sparse_embed_destroy(val_sparse)
    call sparse_embed_destroy(cond_trans)
    call sparse_embed_destroy(P1_evec)
    call sparse_embed_destroy(evec_P1_evec)
    call dense_destroy(kernel_dens_x)
    call dense_destroy(kernel_dens_y)
    call sparse_embed_destroy(x_vec)
    call sparse_embed_destroy(y_vec)
    deallocate(cond_limit, stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_analysis','cond_limit',ierr)
    deallocate(val_limit, stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_analysis','val_limit',ierr)
    deallocate(cond_start, stat=ierr)
    call utils_dealloc_check('lr_tddft_RPA_analysis','cond_start',ierr)


    call timer_clock('lr_tddft_RPA_analysis',2)

  end subroutine lr_tddft_RPA_analysis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module lr_tddft_RPA
