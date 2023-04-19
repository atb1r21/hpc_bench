! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!         General Linear-Response module                         !
!                                                                !
! This module contains routines associated with calculating the  !
! first order response of a Kohn-Sham system to a perturbation   !
! of the ground state density. The routines are mainly needed in !
! the LR-TDDFT module, but are expected to be of more general    !
! use, for example in calculating phonons from response theory   !
!----------------------------------------------------------------!
! This module was created by Tim Zuehlsdorff in 2015.            !
! Modified to add embedding routines by Joseph Prentice, July    !
! 2018                                                           !
! Modified to add initial support for HF-contribution to         !
! SCF response potential matrix by James Womack, 2019            !
!================================================================!

module linear_response

  use constants, only: DP

  implicit none

  private

  public :: linear_response_calc_SCF_2nd
  public :: linear_response_calc_SCF
  public :: linear_response_operator
  public :: linear_response_RPA_operator
  public :: linear_response_finite_diff

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine linear_response_kxc(fxc_buffer,rho_perturbed,&
       ground_state_dens,grid,cell,nhat_den_grad,nhat_den_grad_gs)

    !==============================================================!
    ! Subroutine approximates the second-order selfconsistent      !
    ! response potential by a finite differences approximation.    !
    ! For semi-local xc functionals, V^{2}_{SCF}(r) can be written !
    ! as k_xc(r)*rho^{1}(r)*rho^{1}(r), where k_xc is the 3rd      !
    ! derivative of the exchange correlation energy. Since k_xc    !
    ! is difficult to construct on the standard fine grid for GGAs,!
    ! it is computed here via a finite difference approximation.   !
    ! This routine currently only works for Singlets!              !
    !==============================================================!

    use cell_grid, only: GRID_INFO
    use rundat, only: pub_num_spins,&
        pub_aug_den_dim
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check
    use xc, only: xc_fxc_finite_diff

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    real(kind=DP), intent(inout) :: fxc_buffer(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: rho_perturbed(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: ground_state_dens(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(in) :: nhat_den_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)
    real(kind=DP), optional, intent(inout) :: nhat_den_grad_gs(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)

    ! local variables
    real(kind=DP) :: epsilon_val
    integer :: ierr
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: temp_fxc
    real(kind=DP), allocatable, dimension(:,:,:,:) :: eff_dens

    ! finite difference parameter
    epsilon_val=5.0e-4_DP
    allocate(temp_fxc(grid%ld1,grid%ld2,&
          grid%max_slabs12,pub_num_spins,2), stat=ierr)
    call utils_alloc_check('linear_response_kxc','temp_fxc',ierr)
    allocate(eff_dens(grid%ld1,grid%ld2,&
          grid%max_slabs12,pub_num_spins), stat=ierr)
    call utils_alloc_check('linear_response_kxc','eff_dens',ierr)

    ! use a 5 point central difference approximation to express the 2nd
    ! derivative:
    eff_dens=2.0_DP*ground_state_dens+2.0_DP*epsilon_val*rho_perturbed

    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
      nhat_den_grad_gs=2.0_DP*nhat_den_grad_gs+2.0_DP*epsilon_val*nhat_den_grad

      call xc_fxc_finite_diff(eff_dens,temp_fxc, grid, cell, ground_state_dens,&
           pub_aug_den_dim,nhat_den_grad_gs)
    else
       call xc_fxc_finite_diff(eff_dens,temp_fxc,&
            grid,cell,ground_state_dens,0)
    endif

    fxc_buffer(:,:,:,:)=-1.0_DP*(temp_fxc(:,:,:,:,1))

    ! 2nd point:
    eff_dens=eff_dens-1.0_DP*epsilon_val*rho_perturbed

    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
      nhat_den_grad_gs=nhat_den_grad_gs-1.0_DP*epsilon_val*nhat_den_grad

      call xc_fxc_finite_diff(eff_dens,temp_fxc, grid, cell,ground_state_dens,&
            pub_aug_den_dim,nhat_den_grad_gs)
    else
       call xc_fxc_finite_diff(eff_dens,temp_fxc,&
            grid,cell,ground_state_dens,0)
    endif

    fxc_buffer(:,:,:,:)=fxc_buffer(:,:,:,:)+ &
             16.0_DP*(temp_fxc(:,:,:,:,1))

    ! 3rd point:
    eff_dens=eff_dens-1.0_DP*epsilon_val*rho_perturbed

    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
      nhat_den_grad_gs=nhat_den_grad_gs-1.0_DP*epsilon_val*nhat_den_grad

      call xc_fxc_finite_diff(eff_dens,temp_fxc, grid, cell,ground_state_dens,&
            pub_aug_den_dim,nhat_den_grad_gs)
    else
       call xc_fxc_finite_diff(eff_dens,temp_fxc,&
            grid,cell,ground_state_dens,0)
    endif

    fxc_buffer(:,:,:,:)=fxc_buffer(:,:,:,:)-&
             30.0_DP*(temp_fxc(:,:,:,:,1))

    ! 4th point:
    eff_dens=eff_dens-1.0_DP*epsilon_val*rho_perturbed

    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
      nhat_den_grad_gs=nhat_den_grad_gs-1.0_DP*epsilon_val*nhat_den_grad

      call xc_fxc_finite_diff(eff_dens,temp_fxc, grid, cell,ground_state_dens,&
           pub_aug_den_dim, nhat_den_grad_gs)
    else
       call xc_fxc_finite_diff(eff_dens,temp_fxc,&
            grid,cell,ground_state_dens,0)
    endif

      fxc_buffer(:,:,:,:)=fxc_buffer(:,:,:,:)+ &
             16.0_DP*(temp_fxc(:,:,:,:,1))

    ! 5th point:
    eff_dens=eff_dens-1.0_DP*epsilon_val*rho_perturbed

    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
      nhat_den_grad_gs=nhat_den_grad_gs-1.0_DP*epsilon_val*nhat_den_grad

      call xc_fxc_finite_diff(eff_dens,temp_fxc, grid, cell,ground_state_dens,&
           pub_aug_den_dim, nhat_den_grad_gs)

      ! restore nhat_den_grad_gs:
      nhat_den_grad_gs=nhat_den_grad_gs+2.0_DP*epsilon_val*nhat_den_grad
      nhat_den_grad_gs=nhat_den_grad_gs/2.0_DP
    else
       call xc_fxc_finite_diff(eff_dens,temp_fxc,&
            grid,cell,ground_state_dens,0)
    endif

    fxc_buffer(:,:,:,:)=fxc_buffer(:,:,:,:)- &
             1.0_DP*(temp_fxc(:,:,:,:,1))

    ! scale the resulting value according to finite diff approx:

    fxc_buffer=fxc_buffer/(12.0_DP*epsilon_val*epsilon_val)

    ! deallocate data
    deallocate(temp_fxc,stat=ierr)
    call utils_dealloc_check('linear_response_kxc','temp_fxc',ierr)
    deallocate(eff_dens,stat=ierr)
    call utils_dealloc_check('linear_response_kxc','eff_dens',ierr)

  end subroutine linear_response_kxc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linear_response_calc_SCF_2nd(pot_SCF,rho_perturbed,&
      response_denskern,basis_bra,bra_rep,basis_ket,ket_rep,&
      ground_state_dens, val_denskern, mdl, Vscf_fine)
    !=======================================================================!
    ! The subroutine calculates the second order self consistent field  re- !
    ! sponse of the perturbation described by rho_perturbed in bra/ket      !
    ! matrix representation.                                                !
    ! NOTE THAT THIS ROUTINE IS NOT YET CONSISTENT WITH THE PAW FORMALISM   !
    !=======================================================================!

    use cell_grid, only: GRID_INFO
    use integrals, only: integrals_locpot
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_num_spins, pub_paw
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_scale
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    real(kind=DP), intent(inout) :: rho_perturbed(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2,mdl%fine_grid%max_slabs12,pub_num_spins)
    type(SPAM3), intent(in) :: response_denskern(pub_num_spins)
    type(SPAM3), intent(inout) :: pot_SCF
    type(FUNC_BASIS), intent(in) :: basis_bra
    type(NGWF_REP), intent(in) :: bra_rep
    type(FUNC_BASIS), intent(in) :: basis_ket
    type(NGWF_REP), intent(in) :: ket_rep
    real(kind=DP), intent(inout) :: ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2,mdl%fine_grid%max_slabs12,pub_num_spins)
    type(SPAM3), intent(in) :: val_denskern(pub_num_spins)
    real(kind=DP), optional, intent(inout) ::Vscf_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2,mdl%fine_grid%max_slabs12,pub_num_spins)

    ! local variables
    real(kind=DP), allocatable, dimension (:,:,:,:) :: SCF_buffer
    integer :: ierr

    call timer_clock('linear_response_calc_SCF_2nd',1)

    call utils_assert(pub_paw .eqv. .false., &
         'Subroutine linear_response_calc_SCF_2nd not compatible with PAW.')

    ! allocate temporary storage space
    allocate(SCF_buffer(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
    call utils_alloc_check('linear_response_calc_SCF_2nd', &
         'SCF_buffer', ierr)

    ! calculate SCF buffer
    call linear_response_kxc(SCF_buffer,rho_perturbed,&
       ground_state_dens,mdl%fine_grid,mdl%cell)

    if(present(Vscf_fine)) then
      ! scale for spin degeneracy
      Vscf_fine=2.0_DP*SCF_buffer
    endif

    ! ignore PAW corrections for the time being and directly
    ! construct the matrix elements of pot_SCF
    call integrals_locpot(pot_SCF, bra_rep%ngwfs_on_grid(1), &
         basis_bra, ket_rep%ngwfs_on_grid(1), basis_ket, &
         mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
         SCF_buffer)

    ! scale according to spin degeneracy
    ! does this scaling make sense?
    call sparse_scale(pot_SCF, 2.0_DP)

    ! deallocate storage space
    deallocate(SCF_buffer, stat=ierr)
    call utils_dealloc_check('linear_response_calc_SCF_2nd',&
         'SCF_buffer',ierr)

    call timer_clock('linear_response_calc_SCF_2nd',2)

  end subroutine linear_response_calc_SCF_2nd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linear_response_finite_diff_sp(fxc_buffer,rho_perturbed,&
       ground_state_dens,grid,cell,nsub,nhat_den_grad,nhat_den_grad_gs,&
       sub_rho_perturbed,sub_ground_state_dens,sub_nhat_den_grad,sub_nhat_den_grad_gs)

    !==============================================================!
    ! Subroutine calculates rho_perturbed*fxc in a finite diff     !
    ! approximation in the case the calculation is spin-polarised. !
    ! Only needed for triplet states.                              !
    ! This is to ensure that no higher order de-                   !
    ! rivatives of the exchange-correlation functional are needed  !
    ! which can be troublesome to compute on the standard grid for !
    ! gradient-corrected functionals.                              !
    ! Modified for embedding by Joseph Prentice, July 2018         !
    !==============================================================!
    ! Arguments:                                                   !
    ! fxc_buffer(inout) : rho_perturbed*fxc in the finite diff     !
    !                      approximation                           !
    ! rho_perturbed(in) : response density associated with P^{1}   !
    ! ground_state_dens : ground state density of the system       !
    !==============================================================!

    use cell_grid, only: GRID_INFO
    use rundat, only: pub_num_spins,pub_aug_den_dim,pub_emft,pub_active_region
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check
    use xc, only: xc_fxc_finite_diff, xc_fxc_finite_diff_emft

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    integer, intent(in)         :: nsub
    real(kind=DP), intent(inout) :: fxc_buffer(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins,nsub)
    real(kind=DP), intent(inout) :: rho_perturbed(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: ground_state_dens(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(in) :: nhat_den_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)
    real(kind=DP), optional, intent(inout) :: nhat_den_grad_gs(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)
    real(kind=DP), optional, intent(inout) :: sub_rho_perturbed(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(inout) :: sub_ground_state_dens(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(in) :: sub_nhat_den_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)
    real(kind=DP), optional, intent(inout) :: sub_nhat_den_grad_gs(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)

    ! local variables
    real(kind=DP) :: epsilon_val
    integer :: ierr, is, isub
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: temp_fxc
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: temp_fxc_emft
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: fxc_up_down
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: fxc_up_down_emft
    real(kind=DP), allocatable, dimension(:,:,:,:) :: eff_dens
    real(kind=DP), allocatable, dimension(:,:,:,:) :: sub_eff_dens
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: eff_nhat_den
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: sub_eff_nhat_den

    ! finite difference parameter
    epsilon_val=1.0e-4_DP

    allocate(temp_fxc(grid%ld1,grid%ld2,&
          grid%max_slabs12,pub_num_spins,2), stat=ierr)
    call utils_alloc_check('linear_response_finite_diff_sp','temp_fxc',ierr)
    allocate(fxc_up_down(grid%ld1,grid%ld2,&
          grid%max_slabs12,pub_num_spins,pub_num_spins), stat=ierr)
    call utils_alloc_check('linear_response_finite_diff_sp','fxc_up_down',ierr)
    allocate(eff_dens(grid%ld1,grid%ld2,&
          grid%max_slabs12,pub_num_spins), stat=ierr)
    call utils_alloc_check('linear_response_finite_diff_sp','eff_dens',ierr)
    if (pub_emft) then
       allocate(temp_fxc_emft(grid%ld1,grid%ld2,&
            grid%max_slabs12,pub_num_spins,2), stat=ierr)
       call utils_alloc_check('linear_response_finite_diff_sp',&
            'temp_fxc_emft',ierr)
       allocate(fxc_up_down_emft(grid%ld1,grid%ld2,&
            grid%max_slabs12,pub_num_spins,pub_num_spins), stat=ierr)
       call utils_alloc_check('linear_response_finite_diff_sp',&
            'fxc_up_down_emft',ierr)
       allocate(sub_eff_dens(grid%ld1,grid%ld2,&
            grid%max_slabs12,pub_num_spins), stat=ierr)
       call utils_alloc_check('linear_response_finite_diff_sp',&
            'eff_dens_emft',ierr)
    end if
    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
       allocate(eff_nhat_den(grid%ld1,grid%ld2,grid%max_slabs12,&
            pub_num_spins,pub_aug_den_dim), stat=ierr)
       call utils_alloc_check('linear_response_finite_diff_sp','eff_nhat_den',&
            ierr)
       if (pub_emft) then
          allocate(sub_eff_nhat_den(grid%ld1,grid%ld2,grid%max_slabs12,&
               pub_num_spins,pub_aug_den_dim), stat=ierr)
          call utils_alloc_check('linear_response_finite_diff_sp',&
               'sub_eff_nhat_den',ierr)
       end if
    endif

    ! use a central difference scheme: rho0+eps*rho1
    ! in both spin channels
    eff_dens=ground_state_dens
    if (pub_emft) sub_eff_dens=sub_ground_state_dens
    do is=1, pub_num_spins
       eff_dens(:,:,:,is)=eff_dens(:,:,:,is)+epsilon_val*&
            rho_perturbed(:,:,:,is)
       if (pub_emft) sub_eff_dens(:,:,:,is)=sub_eff_dens(:,:,:,is)+epsilon_val*&
            sub_rho_perturbed(:,:,:,is)
    enddo

    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
       eff_nhat_den=nhat_den_grad_gs
       if (pub_emft) sub_eff_nhat_den=sub_nhat_den_grad_gs
       do is=1, pub_num_spins
          eff_nhat_den(:,:,:,is,:)=eff_nhat_den(:,:,:,is,:)+&
               epsilon_val*nhat_den_grad(:,:,:,is,:)
          if (pub_emft) sub_eff_nhat_den(:,:,:,is,:)=sub_eff_nhat_den(:,:,:,is,:)+&
               epsilon_val*sub_nhat_den_grad(:,:,:,is,:)
       enddo
       call xc_fxc_finite_diff(eff_dens,temp_fxc,grid,cell,ground_state_dens,&
            pub_aug_den_dim,eff_nhat_den,nhat_den_grad_gs)

       if (pub_emft) call xc_fxc_finite_diff_emft(sub_eff_dens,temp_fxc_emft,&
            grid,cell,sub_ground_state_dens,pub_aug_den_dim,sub_eff_nhat_den,&
            sub_nhat_den_grad_gs)

    else

       call xc_fxc_finite_diff(eff_dens,temp_fxc,grid,cell,ground_state_dens,0)

       if (pub_emft) call xc_fxc_finite_diff_emft(sub_eff_dens,temp_fxc_emft,&
            grid,cell,sub_ground_state_dens,0)

    endif

    fxc_up_down(:,:,:,:,1)=temp_fxc(:,:,:,:,1)
    if (pub_emft) fxc_up_down_emft(:,:,:,:,1)=temp_fxc_emft(:,:,:,:,1)

    ! second part of finite difference scheme: rho0-eps*rho1
    do is=1, pub_num_spins
       eff_dens(:,:,:,is)=eff_dens(:,:,:,is)-2.0_DP*epsilon_val*&
            rho_perturbed(:,:,:,is)
       if (pub_emft) sub_eff_dens(:,:,:,is)=sub_eff_dens(:,:,:,is)-2.0_DP*epsilon_val*&
            sub_rho_perturbed(:,:,:,is)
    enddo

    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
       do is=1, pub_num_spins
          eff_nhat_den(:,:,:,is,:)=eff_nhat_den(:,:,:,is,:)-&
               2.0_DP*epsilon_val*nhat_den_grad(:,:,:,is,:)
          if (pub_emft) sub_eff_nhat_den(:,:,:,is,:)=sub_eff_nhat_den(:,:,:,is,:)-&
               2.0_DP*epsilon_val*sub_nhat_den_grad(:,:,:,is,:)
       enddo
       call xc_fxc_finite_diff(eff_dens,temp_fxc,grid,cell,ground_state_dens,&
            pub_aug_den_dim,eff_nhat_den,nhat_den_grad_gs)

       if (pub_emft) call xc_fxc_finite_diff_emft(sub_eff_dens,temp_fxc_emft,&
            grid,cell,sub_ground_state_dens,pub_aug_den_dim,sub_eff_nhat_den,&
            sub_nhat_den_grad_gs)

    else
       call xc_fxc_finite_diff(eff_dens,temp_fxc,grid,cell,ground_state_dens,0)

       if (pub_emft) call xc_fxc_finite_diff_emft(sub_eff_dens,temp_fxc_emft,&
            grid,cell,sub_ground_state_dens,0)

    endif

    fxc_up_down(:,:,:,:,1)=fxc_up_down(:,:,:,:,1)-temp_fxc(:,:,:,:,1)
    if (pub_emft) fxc_up_down_emft(:,:,:,:,1)=fxc_up_down_emft(:,:,:,:,1)-&
         temp_fxc_emft(:,:,:,:,1)

    ! now apply finite difference scaling
    fxc_up_down=0.5_DP*fxc_up_down/epsilon_val
    fxc_up_down_emft=0.5_DP*fxc_up_down_emft/epsilon_val

    ! now build correct fxc in both spin channels
    do is=1, pub_num_spins
       ! jcap: loop over regions
       do isub=1,nsub
          fxc_buffer(:,:,:,is,isub)=fxc_up_down(:,:,:,is,1)
          if (isub == pub_active_region .and. pub_emft) &
               fxc_buffer(:,:,:,is,isub)=fxc_buffer(:,:,:,is,isub)+&
               fxc_up_down_emft(:,:,:,is,1)
       end do
    enddo

    deallocate(temp_fxc,stat=ierr)
    call utils_dealloc_check('linear_response_finite_diff_sp','temp_fxc',ierr)
    deallocate(eff_dens,stat=ierr)
    call utils_dealloc_check('linear_response_finite_diff_sp','eff_dens',ierr)
    deallocate(fxc_up_down,stat=ierr)
    call utils_dealloc_check('linear_response_finite_diff_sp',&
         'fxc_up_down',ierr)
    if (pub_emft) then
       deallocate(temp_fxc_emft,stat=ierr)
       call utils_dealloc_check('linear_response_finite_diff_sp',&
            'temp_fxc_emft',ierr)
       deallocate(sub_eff_dens,stat=ierr)
       call utils_dealloc_check('linear_response_finite_diff_sp',&
            'sub_eff_dens',ierr)
       deallocate(fxc_up_down_emft,stat=ierr)
       call utils_dealloc_check('linear_response_finite_diff_sp',&
            'fxc_up_down_emft',ierr)
    end if
    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
       deallocate(eff_nhat_den, stat=ierr)
       call utils_dealloc_check('linear_response_finite_diff_sp',&
            'eff_nhat_den',ierr)
       if (pub_emft) then
          deallocate(sub_eff_nhat_den, stat=ierr)
          call utils_dealloc_check('linear_response_finite_diff_sp',&
               'sub_eff_nhat_den',ierr)
       end if
    endif

  end subroutine linear_response_finite_diff_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linear_response_finite_diff(fxc_buffer,rho_perturbed,&
       ground_state_dens,grid,cell,nsub,nhat_den_grad,nhat_den_grad_gs,&
       sub_rho_perturbed,sub_ground_state_dens,sub_nhat_den_grad,sub_nhat_den_grad_gs)

    !==============================================================!
    ! Subroutine calculates rho_perturbed*fxc in a finite diff     !
    ! approximation. This is to ensure that no higher order de-    !
    ! rivatives of the exchange-correlation functional are needed  !
    ! which can be troublesome to compute on the standard grid for !
    ! gradient-corrected functionals.                              !
    ! Modified for embedding by Joseph Prentice, July 2018         !
    !==============================================================!
    ! Arguments:                                                   !
    ! fxc_buffer(inout) : rho_perturbed*fxc in the finite diff     !
    !                      approximation                           !
    ! rho_perturbed(in) : response density associated wiht P^{1}   !
    ! ground_state_dens : ground state density of the system       !
    !==============================================================!

    use cell_grid, only: GRID_INFO
    use rundat, only: pub_lr_tddft_triplet, pub_num_spins,&
        pub_aug_den_dim, pub_emft, pub_active_region
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check
    use xc, only: xc_fxc_finite_diff, xc_fxc_finite_diff_emft

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    integer, intent(in)         :: nsub
    real(kind=DP), intent(inout) :: fxc_buffer(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins,nsub)
    real(kind=DP), intent(inout) :: rho_perturbed(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: ground_state_dens(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(in) :: nhat_den_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)
    real(kind=DP), optional, intent(inout) :: nhat_den_grad_gs(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)
    real(kind=DP), optional, intent(inout) :: sub_rho_perturbed(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(inout) :: sub_ground_state_dens(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(in) :: sub_nhat_den_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)
    real(kind=DP), optional, intent(inout) :: sub_nhat_den_grad_gs(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)


    ! local variables
    real(kind=DP) :: epsilon_val
    integer :: ierr, isub
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: temp_fxc
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: temp_fxc_emft
    real(kind=DP), allocatable, dimension(:,:,:,:) :: eff_dens
    real(kind=DP), allocatable, dimension(:,:,:,:) :: sub_eff_dens
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: eff_nhat_den
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: sub_eff_nhat_den

    ! finite difference parameter
    epsilon_val=1.0e-4_DP
    allocate(temp_fxc(grid%ld1,grid%ld2,&
          grid%max_slabs12,pub_num_spins,2), stat=ierr)
    call utils_alloc_check('linear_response_finite_diff','temp_fxc',ierr)
    ! @optimize
    ! temp_fxc always has final dimension 2, even though temp_fxc(:,:,:,:,2)
    ! seems to only be used when pub_lr_tddft_triplet is true. It is otherwise
    ! a waste of memory.
    allocate(eff_dens(grid%ld1,grid%ld2,&
          grid%max_slabs12,pub_num_spins), stat=ierr)
    call utils_alloc_check('linear_response_finite_diff','eff_dens',ierr)
    if (pub_emft) then
       allocate(temp_fxc_emft(grid%ld1,grid%ld2,&
            grid%max_slabs12,pub_num_spins,2), stat=ierr)
       call utils_alloc_check('linear_response_finite_diff',&
            'temp_fxc_emft',ierr)
       allocate(sub_eff_dens(grid%ld1,grid%ld2,&
            grid%max_slabs12,pub_num_spins), stat=ierr)
       call utils_alloc_check('linear_response_finite_diff',&
            'eff_dens_emft',ierr)
    end if
    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
       allocate(eff_nhat_den(grid%ld1,grid%ld2,grid%max_slabs12,&
         pub_num_spins,pub_aug_den_dim), stat=ierr)
       call utils_alloc_check('linear_response_finite_diff','eff_nhat_den',&
         ierr)
       if (pub_emft) then
          allocate(sub_eff_nhat_den(grid%ld1,grid%ld2,grid%max_slabs12,&
               pub_num_spins,pub_aug_den_dim), stat=ierr)
          call utils_alloc_check('linear_response_finite_diff',&
               'sub_eff_nhat_den',ierr)
       end if
    endif

    ! use a central difference scheme: rho0+eps*rho1
    ! double ground state dens due to spin degeneracy
    eff_dens=2.0_DP*ground_state_dens+epsilon_val*rho_perturbed
    if (pub_emft) sub_eff_dens=2.0_DP*sub_ground_state_dens+&
         epsilon_val*sub_rho_perturbed

    ! if nhat_den_grad is supplied, the finite difference approximation
    ! has to be applied to it as well. It is stored temporarily in
    ! nhat_den_grad_gs
    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
       eff_nhat_den=2.0_DP*nhat_den_grad_gs+epsilon_val*nhat_den_grad
       if (pub_emft) sub_eff_nhat_den=2.0_DP*sub_nhat_den_grad_gs+&
            epsilon_val*sub_nhat_den_grad

      call xc_fxc_finite_diff(eff_dens,temp_fxc,grid,cell,ground_state_dens,&
           pub_aug_den_dim,eff_nhat_den,nhat_den_grad_gs)

      if (pub_emft) call xc_fxc_finite_diff_emft(sub_eff_dens,temp_fxc_emft,&
           grid,cell,sub_ground_state_dens,pub_aug_den_dim,sub_eff_nhat_den,&
           sub_nhat_den_grad_gs)

    else
       call xc_fxc_finite_diff(eff_dens,temp_fxc,grid,cell,ground_state_dens,0)

       if (pub_emft) call xc_fxc_finite_diff_emft(sub_eff_dens,temp_fxc_emft,&
            grid,cell,sub_ground_state_dens,0)

    endif

    ! jcap: loop over regions
    do isub=1,nsub
       if(pub_lr_tddft_triplet) then
          fxc_buffer(:,:,:,:,isub)=(temp_fxc(:,:,:,:,1)-temp_fxc(:,:,:,:,2))
          if (pub_emft.and.(isub==pub_active_region)) then
             fxc_buffer(:,:,:,:,isub)=fxc_buffer(:,:,:,:,isub)+&
                  (temp_fxc_emft(:,:,:,:,1)-temp_fxc_emft(:,:,:,:,2))
          end if
       else
          fxc_buffer(:,:,:,:,isub)=temp_fxc(:,:,:,:,1)
          if (pub_emft.and.(isub==pub_active_region)) then
             fxc_buffer(:,:,:,:,isub)=fxc_buffer(:,:,:,:,isub)+&
                  temp_fxc_emft(:,:,:,:,1)
          end if
       end if
    end do

    ! second part of finite difference scheme: rho0-eps*rho1
    eff_dens=2.0_DP*ground_state_dens-epsilon_val*rho_perturbed
    if (pub_emft) sub_eff_dens=2.0_DP*sub_ground_state_dens-&
         epsilon_val*sub_rho_perturbed

    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
       eff_nhat_den=2.0_DP*nhat_den_grad_gs-epsilon_val*nhat_den_grad
       if (pub_emft) sub_eff_nhat_den=2.0_DP*sub_nhat_den_grad_gs-&
            epsilon_val*sub_nhat_den_grad

       call xc_fxc_finite_diff(eff_dens,temp_fxc,grid,cell,ground_state_dens,&
            pub_aug_den_dim,eff_nhat_den,nhat_den_grad_gs)

       if (pub_emft) call xc_fxc_finite_diff_emft(sub_eff_dens,temp_fxc_emft,&
            grid,cell,sub_ground_state_dens,pub_aug_den_dim,sub_eff_nhat_den,&
            sub_nhat_den_grad_gs)

    else
       call xc_fxc_finite_diff(eff_dens,temp_fxc,grid,cell,ground_state_dens,0)

       if (pub_emft) call xc_fxc_finite_diff_emft(sub_eff_dens,temp_fxc_emft,&
            grid,cell,sub_ground_state_dens,0)

    endif

    do isub=1,nsub
       if(pub_lr_tddft_triplet) then
          fxc_buffer(:,:,:,:,isub)=fxc_buffer(:,:,:,:,isub)-&
               (temp_fxc(:,:,:,:,1)-temp_fxc(:,:,:,:,2))
          if (pub_emft.and.(isub==pub_active_region)) then
             fxc_buffer(:,:,:,:,isub)=fxc_buffer(:,:,:,:,isub)-&
                  (temp_fxc_emft(:,:,:,:,1)-temp_fxc_emft(:,:,:,:,2))
          end if
       else
          fxc_buffer(:,:,:,:,isub)=fxc_buffer(:,:,:,:,isub)-&
               temp_fxc(:,:,:,:,1)
          if (pub_emft.and.(isub==pub_active_region)) then
             fxc_buffer(:,:,:,:,isub)=fxc_buffer(:,:,:,:,isub)-&
                  temp_fxc_emft(:,:,:,:,1)
          end if
       end if
    end do

    ! scale by epsilon to compute the finite difference approximation
    fxc_buffer=0.5_DP*fxc_buffer/epsilon_val

    ! sum over spins
    fxc_buffer=2.0_DP*fxc_buffer

    deallocate(temp_fxc,stat=ierr)
    call utils_dealloc_check('linear_response_finite_diff','temp_fxc',ierr)
    deallocate(eff_dens,stat=ierr)
    call utils_dealloc_check('linear_response_finite_diff','eff_dens',ierr)
    if (pub_emft) then
       deallocate(temp_fxc_emft,stat=ierr)
       call utils_dealloc_check('linear_response_finite_diff',&
            'temp_fxc_emft',ierr)
       deallocate(sub_eff_dens,stat=ierr)
       call utils_dealloc_check('linear_response_finite_diff',&
            'sub_eff_dens',ierr)
    end if
    if(present(nhat_den_grad) .and. present(nhat_den_grad_gs)) then
       deallocate(eff_nhat_den, stat=ierr)
       call utils_dealloc_check('linear_response_finite_diff','eff_nhat_den',&
            ierr)
       if (pub_emft) then
          deallocate(sub_eff_nhat_den, stat=ierr)
          call utils_dealloc_check('linear_response_finite_diff',&
               'sub_eff_nhat_den',ierr)
       end if
    endif

  end subroutine linear_response_finite_diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linear_response_calc_SCF(pot_SCF,fxc_fine, rho_perturbed, &
       sub_rho_perturbed, response_denskern, basis_bra, bra_rep, &
       basis_ket, ket_rep, ground_state_dens, sub_ground_state_dens, &
       val_denskern, mdl, hfxstate, ground_state_hartree)
    !========================================================================!
    ! The subroutine calculates the selfconsistent response of a Kohn-Sham   !
    ! System to a perturbation in the density described by rho_perturbed, a  !
    ! density associated with a generalised response_denskern. The routine   !
    ! calculates V_SCF(r)=int rho^1(r')/|r-r'|+fxc(r)*rho^1(r) and           !
    ! returns pot_SCF in matrix form in <bra_NGWF|V_SCF|ket_NGWF> represen-  !
    ! tation and is written in a flexible way to ensure that the evaluation  !
    ! is compatible with the use of an implicit solvent or PAW. Note that    !
    ! for the moment, this routine only works for spin-degenerate systems.   !
    ! Modified for embedding by Joseph Prentice, July 2018                   !
    ! Modified to support XC with HFX by James Womack in 2019.               !
    !========================================================================!
    ! Arguments:                                                             !
    ! pot_SCF (inout) : V_SCF in bra/ket NGWF representation                 !
    ! fxc_fine (in)   : exchange-correlation kernel on fine grid. Only used  !
    !               for the LDA functional. Otherwise the finite difference  !
    !               approximation is used.                                   !
    ! rho_perturbed   : response density associated with P^{1}               !
    ! response_denskern : response_denskern P^{1}                            !
    ! hfxstate (inout): The HFX_STATE container that stores all HFx data     !
    !                   that needs to persist across calls.                  !
    ! ground_state_hartree : Only needed if the implicit solvent model is    !
    !               used for the hartree potential. Contains the hartree     !
    !               potential of the ground state density                    !
    !========================================================================!

    use augmentation, only: augmentation_density_on_grid,&
        augmentation_screen_dij, aug_projector_denskern
    use cell_grid, only: GRID_INFO
    use comms, only: comms_barrier
    use constants, only: REP_SWEX_HFX
    use function_basis, only: FUNC_BASIS
    use hartree, only: hartree_on_grid, hartree_via_multigrid
    use hf_exchange, only: HFX_STATE, hf_exchange_calculate, &
         hf_exchange_dkn_indep_stage
    use integrals, only: integrals_locpot
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use paw, only: paw_dij_hartree, paw_dij_xc_tddft
    use rundat, only: pub_lr_tddft_triplet,pub_lr_tddft_xc_finite_diff, &
         pub_num_spins, pub_turn_off_hartree, pub_multigrid_hartree,&
         pub_nlcc, pub_paw, pub_nhat_in_xc, pub_aug_den_dim, pub_aug,&
         pub_spin_fac, pub_emft, pub_active_region, &
         pub_use_hfx, pub_use_activehfx, pub_devel_code
    use simulation_cell, only: CELL_INFO
    use sparse, only: sparse_create, sparse_destroy, sparse_transpose, &
         sparse_transpose_structure
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_transpose, sparse_embed_axpy, sparse_embed_copy, &
         sparse_embed_product, sparse_embed_transpose_structure, sparse_embed_scale, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use sparse, only: SPAM3
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_devel_code
    use xc, only : xc_embed_swap_functional
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    real(kind=DP), intent(in) :: fxc_fine(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2,mdl%fine_grid%max_slabs12,2,mdl%nsub)
    real(kind=DP), intent(inout) :: rho_perturbed(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2,mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_rho_perturbed(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2,mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)
    type(SPAM3_EMBED), intent(in) :: response_denskern(pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: pot_SCF(pub_num_spins)
    ! JCW: bra_rep and ket_rep must be inout, since (when computing the HF
    !      contribution to pot_SCF, we use hf_exchange_dkn_indep_stage
    !      to compute and stores SW expansions of NGWF products in
    !      each instance of NGWF_REP (bra_rep and ket_rep)
    type(FUNC_BASIS), intent(in) :: basis_bra(:)
    type(NGWF_REP), intent(inout) :: bra_rep
    type(FUNC_BASIS), intent(in) :: basis_ket(:)
    type(NGWF_REP), intent(inout) :: ket_rep
    real(kind=DP), intent(inout) :: ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2,mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2,mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)
    type(SPAM3_EMBED), intent(in) :: val_denskern(pub_num_spins)
    type(HFX_STATE), intent(inout), target :: hfxstate
    real(kind=DP), optional, intent(in) :: ground_state_hartree(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2,mdl%fine_grid%max_slabs12,pub_num_spins)

    ! Parameters
    ! JCW: Factor to scale matrix returned by hf_exchange_calculate (already
    !      includes c_HF)
    real(kind=DP) :: Vhf_scale_factor = -2.0_DP

    ! Local Variables
    real(kind=DP), allocatable, dimension (:,:,:,:,:) :: SCF_buffer
    real(kind=DP), allocatable, dimension (:,:,:,:) :: hartree_buffer
    real(kind=DP), allocatable, dimension (:,:,:,:,:,:) :: sub_hartree_buffer
    real(kind=DP), allocatable, dimension (:,:,:,:,:) :: fxc_buffer
    real(kind=DP), allocatable, dimension (:,:,:,:,:) :: nhat_den_grad_gs
    real(kind=DP), allocatable, dimension (:,:,:,:,:) :: nhat_den_grad
    real(kind=DP), allocatable, dimension (:,:,:,:,:,:,:) :: sub_nhat_den_grad_gs
    real(kind=DP), allocatable, dimension (:,:,:,:,:,:,:) :: sub_nhat_den_grad
    type(SPAM3_EMBED), allocatable, dimension(:) :: temp_pot_SCF1, temp_pot_SCF2
    type(SPAM3_EMBED) :: sp_overlap_Dij, ps_overlap, SCF_nonloc
    type(SPAM3_EMBED), allocatable, dimension(:) :: Dij, Dij_hartree,Dij_hat, rho_ij
    type(SPAM3_EMBED), allocatable, dimension(:) :: Dij_xc, rho_0_ij
    type(SPAM3), allocatable, dimension(:) :: temp_array, temp_array2, &
         temp_array3, temp_array4
    real(kind=DP) :: epsilon_val

    integer :: fine_max_slabs12, fine_ld1, fine_ld2 ! jd: shortcut notation

    integer :: ierr, is, isub, jsub, ireg, fxc_index

    logical :: calc_Vhf ! set to .false. to disable evaluation of HF
                        ! contribution to pot_SCF using devel_code
                        ! TDDFT:CALC_VHF=[T/F]:TDDFT

    ! -------------------------------------------------------------------------

    call timer_clock('linear_response_calc_SCF', 1)

    epsilon_val=1.0e-4_DP ! finite difference parameter

    if (pub_emft) ireg=pub_active_region   ! local version

    allocate(hartree_buffer(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins), stat=ierr)
    call utils_alloc_check('linear_response_calc_SCF', 'hartree_buffer', ierr)
    allocate(sub_hartree_buffer(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub, mdl%nsub), stat=ierr)
    call utils_alloc_check('linear_response_calc_SCF', 'sub_hartree_buffer', ierr)
    allocate(SCF_buffer(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub), stat=ierr)
    call utils_alloc_check('linear_response_calc_SCF', 'SCF_buffer', ierr)
    allocate(fxc_buffer(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, pub_num_spins, mdl%nsub), stat=ierr)
    call utils_alloc_check('linear_response_calc_SCF', 'fxc_buffer', ierr)
    allocate(temp_array(pub_num_spins),stat=ierr)
    call utils_alloc_check('linear_response_calc_SCF','temp_array',ierr)

    ! if this is a PAW calculation, we need to construct nhat_grad for the
    ! response density. If this is a multigrid Hartree calculation, we also
    ! need to do this for the ground state density.

    ! Note that this routine, due to the way the lr_tddft calculation works, assumes that
    ! on entering the routine, ground_state_dens IS augmented by nhat, while rho_perturbed
    ! IS NOT
    if(pub_aug) then
       ! allocate and create data structures for nhat_den_grad
       fine_ld1 = mdl%fine_grid%ld1
       fine_ld2 = mdl%fine_grid%ld2
       fine_max_slabs12 = mdl%fine_grid%max_slabs12
       allocate(nhat_den_grad(fine_ld1,fine_ld2,fine_max_slabs12,pub_num_spins,&
            0:pub_aug_den_dim),stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF','nhat_den_grad',ierr)
       allocate(sub_nhat_den_grad(fine_ld1,fine_ld2,fine_max_slabs12,pub_num_spins,&
            0:pub_aug_den_dim,mdl%nsub,mdl%nsub),stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF','sub_nhat_den_grad',ierr)

       ! if nhat is NOT to be included in the xc term, we need to constrct nhat_den_grad for the
       ! ground state density kernel again.
       allocate(nhat_den_grad_gs(fine_ld1,fine_ld2,fine_max_slabs12,&
            pub_num_spins,0:pub_aug_den_dim),stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF','nhat_den_grad_gs',ierr)
       allocate(sub_nhat_den_grad_gs(fine_ld1,fine_ld2,fine_max_slabs12,&
            pub_num_spins,0:pub_aug_den_dim,mdl%nsub,mdl%nsub),stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF','sub_nhat_den_grad_gs',ierr)

       nhat_den_grad = 0.0_DP
       nhat_den_grad_gs = 0.0_DP
       sub_nhat_den_grad = 0.0_DP
       sub_nhat_den_grad_gs = 0.0_DP

       ! jcap: loop over regions and sum up densities
       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub

             call sparse_embed_extract_from_array(temp_array,response_denskern,&
                  isub,jsub)
             call augmentation_density_on_grid(sub_nhat_den_grad(:,:,:,:,:,isub,jsub), &
                  mdl%fine_grid, mdl%cell, mdl%regions(isub)%pseudo_sp, &
                  mdl%regions(isub)%paw_sp, mdl%aug_box, &
                  temp_array, bra_rep%sp_overlap%m(isub,jsub), &
                  ket_rep%sp_overlap%m(isub,jsub))
             call sparse_embed_destroy_extracted_array(temp_array)

             call sparse_embed_extract_from_array(temp_array,val_denskern,&
                  isub,jsub)
             call augmentation_density_on_grid(sub_nhat_den_grad_gs(:,:,:,:,:,isub,jsub), &
                  mdl%fine_grid, mdl%cell, mdl%regions(isub)%pseudo_sp, &
                  mdl%regions(isub)%paw_sp, mdl%aug_box, &
                  temp_array, ket_rep%sp_overlap%m(isub,jsub))
             call sparse_embed_destroy_extracted_array(temp_array)

             nhat_den_grad=nhat_den_grad+sub_nhat_den_grad(:,:,:,:,:,isub,jsub)
             nhat_den_grad_gs=nhat_den_grad_gs+sub_nhat_den_grad_gs(:,:,:,:,:,isub,jsub)

          end do
       end do

       ! add nhat to response density
       rho_perturbed = rho_perturbed + nhat_den_grad(:,:,:,:,0)
       sub_rho_perturbed = sub_rho_perturbed + sub_nhat_den_grad(:,:,:,:,0,:,:)

    endif

    ! Now first construct the hartree potential before computing the xc_term
    if (.not. pub_turn_off_hartree .and. .not. pub_lr_tddft_triplet) then
       sub_hartree_buffer=0.0_DP
       if(pub_multigrid_hartree) then
          ! in the multigrid method, the hartree potential is evaluated as
          ! the difference between the ground state hartree potential and
          ! the hartree potential of rho_0+rho_perturbed
          ! check for spin degeneracy:
          ! jcap: need to do this calculation over the whole system
          if(pub_num_spins==1) then
             ground_state_dens=2.0_DP*ground_state_dens
             rho_perturbed=ground_state_dens+2.0_DP*rho_perturbed
          else
             rho_perturbed=ground_state_dens+rho_perturbed
          endif

          call hartree_via_multigrid(hartree_buffer, &
               rho_perturbed, mdl%fine_grid, mdl%cell, &
               elements=mdl%elements)

         hartree_buffer=hartree_buffer-ground_state_hartree
         rho_perturbed=rho_perturbed-ground_state_dens

         ! again check for spin degeneracy:
         if(pub_num_spins==1) then
            rho_perturbed=rho_perturbed*0.5_DP

            ground_state_dens=0.5_DP*ground_state_dens
            hartree_buffer=0.5_DP*hartree_buffer
         endif
      else
         ! jcap: loop over regions
         do isub=1,mdl%nsub
            do jsub=1,mdl%nsub
               call hartree_on_grid(sub_hartree_buffer(:,:,:,:,isub,jsub), &
                    sub_rho_perturbed(:,:,:,:,isub,jsub), mdl%fine_grid, mdl%cell)
            end do
         end do

         ! jcap: sum up over regions if we aren't using multigrid hartree
         hartree_buffer=0.0_DP
         do isub=1,mdl%nsub
            do jsub=1,mdl%nsub
               hartree_buffer=hartree_buffer+sub_hartree_buffer(:,:,:,:,isub,jsub)
            end do
         end do

      endif

    else
       hartree_buffer = 0.0_DP
       sub_hartree_buffer = 0.0_DP
    endif

    ! now compute xc part
    if(pub_lr_tddft_xc_finite_diff) then
       ! add core density to ground state density if necessary
       if(pub_nlcc) then
         do is=1,pub_num_spins
            ground_state_dens(:,:,:,is)= ground_state_dens(:,:,:,is) + &
                 mdl%core_density_fine * 0.5_DP * pub_spin_fac
            if (pub_emft) then
               sub_ground_state_dens(:,:,:,is,ireg,ireg) = &
                    sub_ground_state_dens(:,:,:,is,ireg,ireg) + &
                    mdl%regions(ireg)%core_density_fine * 0.5_DP *pub_spin_fac
            end if
         enddo
       endif

       ! see if nhat is included in the calculation of fxc or not
       if(pub_aug .and. .not. pub_nhat_in_xc) then
          rho_perturbed=rho_perturbed-nhat_den_grad(:,:,:,:,0)
          ground_state_dens=ground_state_dens-nhat_den_grad_gs(:,:,:,:,0)
          if (pub_emft) then
             sub_rho_perturbed(:,:,:,:,ireg,ireg)=sub_rho_perturbed(:,:,:,:,ireg,ireg)-&
                  sub_nhat_den_grad(:,:,:,:,0,ireg,ireg)
             sub_ground_state_dens(:,:,:,:,ireg,ireg)=&
                  sub_ground_state_dens(:,:,:,:,ireg,ireg)-&
                  sub_nhat_den_grad_gs(:,:,:,:,0,ireg,ireg)
          end if
       endif

       if(pub_num_spins==2) then
          if (pub_emft) then
             if (pub_aug) then
                call linear_response_finite_diff_sp(fxc_buffer,rho_perturbed,&
                     ground_state_dens,mdl%fine_grid,mdl%cell,mdl%nsub,&
                     nhat_den_grad=nhat_den_grad,nhat_den_grad_gs=nhat_den_grad_gs,&
                     sub_rho_perturbed=sub_rho_perturbed(:,:,:,:,ireg,ireg),&
                     sub_ground_state_dens=sub_ground_state_dens(:,:,:,:,ireg,ireg),&
                     sub_nhat_den_grad=sub_nhat_den_grad(:,:,:,:,:,ireg,ireg),&
                     sub_nhat_den_grad_gs=sub_nhat_den_grad_gs(:,:,:,:,:,ireg,ireg))
             else
                call linear_response_finite_diff_sp(fxc_buffer,rho_perturbed,&
                     ground_state_dens,mdl%fine_grid,mdl%cell,mdl%nsub,&
                     sub_rho_perturbed=sub_rho_perturbed(:,:,:,:,ireg,ireg),&
                     sub_ground_state_dens=sub_ground_state_dens(:,:,:,:,ireg,ireg))
             endif
          else
             if (pub_aug) then
                call linear_response_finite_diff_sp(fxc_buffer,rho_perturbed,&
                     ground_state_dens,mdl%fine_grid,mdl%cell,mdl%nsub,&
                     nhat_den_grad=nhat_den_grad,nhat_den_grad_gs=nhat_den_grad_gs)
             else
                call linear_response_finite_diff_sp(fxc_buffer,rho_perturbed,&
                     ground_state_dens,mdl%fine_grid,mdl%cell,mdl%nsub)
             endif
          end if
       else
          if (pub_emft) then
             if (pub_aug) then
                call linear_response_finite_diff(fxc_buffer,rho_perturbed,&
                     ground_state_dens,mdl%fine_grid,mdl%cell,mdl%nsub,&
                     nhat_den_grad=nhat_den_grad,nhat_den_grad_gs=nhat_den_grad_gs,&
                     sub_rho_perturbed=sub_rho_perturbed(:,:,:,:,ireg,ireg),&
                     sub_ground_state_dens=sub_ground_state_dens(:,:,:,:,ireg,ireg),&
                     sub_nhat_den_grad=sub_nhat_den_grad(:,:,:,:,:,ireg,ireg),&
                     sub_nhat_den_grad_gs=sub_nhat_den_grad_gs(:,:,:,:,:,ireg,ireg))
             else
                call linear_response_finite_diff(fxc_buffer,rho_perturbed,&
                     ground_state_dens,mdl%fine_grid,mdl%cell,mdl%nsub,&
                     sub_rho_perturbed=sub_rho_perturbed(:,:,:,:,ireg,ireg),&
                     sub_ground_state_dens=sub_ground_state_dens(:,:,:,:,ireg,ireg))
             endif
          else
             if (pub_aug) then
                call linear_response_finite_diff(fxc_buffer,rho_perturbed,&
                     ground_state_dens,mdl%fine_grid,mdl%cell,mdl%nsub,&
                     nhat_den_grad=nhat_den_grad,nhat_den_grad_gs=nhat_den_grad_gs)
             else
                call linear_response_finite_diff(fxc_buffer,rho_perturbed,&
                     ground_state_dens,mdl%fine_grid,mdl%cell,mdl%nsub)
             endif
          end if
       endif

       ! subtract core density from ground state density if necessary
       if(pub_nlcc) then
         do is=1,pub_num_spins
            ground_state_dens(:,:,:,is)= ground_state_dens(:,:,:,is) - &
                 mdl%core_density_fine * 0.5_DP * pub_spin_fac
            if (pub_emft) then
               sub_ground_state_dens(:,:,:,is,ireg,ireg)= &
                    sub_ground_state_dens(:,:,:,is,ireg,ireg) - &
                    mdl%regions(ireg)%core_density_fine * 0.5_DP * pub_spin_fac
            end if
         enddo
       endif

       ! add nhat density back to ground state density to be consistent
       if(pub_aug .and. .not. pub_nhat_in_xc) then
          ground_state_dens=ground_state_dens+nhat_den_grad_gs(:,:,:,:,0)
          if (pub_emft) then
             sub_ground_state_dens(:,:,:,:,ireg,ireg)=&
                  sub_ground_state_dens(:,:,:,:,ireg,ireg)+&
                  sub_nhat_den_grad_gs(:,:,:,:,0,ireg,ireg)
          end if
       endif
    else
       ! subtract nhat density from response density if necessary
       if(pub_aug .and. .not. pub_nhat_in_xc) then
          rho_perturbed=rho_perturbed-nhat_den_grad(:,:,:,:,0)
          if (pub_emft) then
             sub_rho_perturbed(:,:,:,:,ireg,ireg)=&
                  sub_rho_perturbed(:,:,:,:,ireg,ireg)-&
                  sub_nhat_den_grad(:,:,:,:,0,ireg,ireg)
          end if
       endif

       ! this subroutine uses precalculated fxc_fine instead of finite diff
       ! approximation. This only works for the spin degenerate case.
       ! jcap: loop over regions
       do isub=1,mdl%nsub
          fxc_buffer(:,:,:,1,isub)=(fxc_fine(:,:,:,1,isub))*&
               rho_perturbed(:,:,:,1)
       end do
    endif

    ! make sure that response_density leaving this routine is nhat free
    if(pub_aug .and. pub_nhat_in_xc) then
       rho_perturbed=rho_perturbed-nhat_den_grad(:,:,:,:,0)
       if (pub_emft) then
          sub_rho_perturbed(:,:,:,:,ireg,ireg)=&
               sub_rho_perturbed(:,:,:,:,ireg,ireg)-&
               sub_nhat_den_grad(:,:,:,:,0,ireg,ireg)
       end if
    endif

    ! combine fxc part and hartree part to get the full selfconsistent field
    if(.not. pub_lr_tddft_triplet) then
       if(pub_num_spins==1) then
          ! Factor of 2 originates from the fact that P^{1}=
          ! 1/sqrt(2)*(P^{1}_UP+P^{1}_DOWN), with up and down
          ! channel being identical in the spin degenerate system.
          ! net result is a factor of 2, once VSCF is multiplied by P^{1}_up
          ! and P^{1}_down
          ! jcap: need to loop over regions
          do isub=1,mdl%nsub
             SCF_buffer(:,:,:,1,isub)=2.0_DP*hartree_buffer(:,:,:,1)+&
                  fxc_buffer(:,:,:,1,isub)
          end do
       else
          ! in spin polarised calculation, hartree buffer is given by
          ! Vhartree up + Vhartree down. This is already done in the
          ! hartree_on_grid routine.
          ! jcap: need to loop over regions
          do isub=1,mdl%nsub
             SCF_buffer(:,:,:,:,isub)=hartree_buffer+fxc_buffer(:,:,:,:,isub)
          end do
       endif
    else
       ! jcap: need to loop over regions
       do isub=1,mdl%nsub
          SCF_buffer(:,:,:,:,isub)=fxc_buffer(:,:,:,:,isub)
       end do
    endif

    call comms_barrier

    allocate(temp_pot_SCF1(pub_num_spins),stat=ierr)
    call utils_alloc_check('linear_response_calc_SCF','temp_pot_SCF1',ierr)
    do is=1,pub_num_spins
       ! Create un-augmented "bare" SCF potential with direct sphere
       ! sphere overlap sparsity
       call sparse_embed_create(temp_pot_SCF1(is),bra_rep%ngwf_cross_overlap_tr)
       ! jcap: loop over regions
       ! jcap: here, follow calculation of ham%lhxc in
       ! hamiltonian_dens_dep_matrices
       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub
             ! jcap: the SCF_buffer region required will depend on the
             ! region we're in
             fxc_index=isub
             if ((isub.ne.jsub).and.(isub==ireg)) fxc_index=jsub
             call integrals_locpot(temp_pot_SCF1(is)%m(isub,jsub), &
                  bra_rep%ngwfs_on_grid(isub), basis_bra(isub), &
                  ket_rep%ngwfs_on_grid(jsub), basis_ket(jsub), &
                  mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
                  SCF_buffer(:,:,:,is,fxc_index))
          end do
       end do
       ! Copy to version where augmented SCF potential will be constructed
       ! if required
       call sparse_embed_copy(pot_SCF(is),temp_pot_SCF1(is))
       call sparse_embed_destroy(temp_pot_SCF1(is))
    enddo
    deallocate(temp_pot_SCF1,stat=ierr)
    call utils_dealloc_check('linear_response_calc_SCF','temp_pot_SCF1',ierr)

    ! JCW: pot_SCF now contains contribution to the SCF response potential
    !      matrix from semi-local parts of the XC functional.

    ! JCW: If needed, we now compute and add the contribution from Hartree-Fock
    !      exchange, as described in Eq. 26 of T.J. Zuehlsdorff et al.,
    !      J. Chem. Phys. 139, 064104 (2013).
    !
    !      Use hf_exchange module to evaluate the HF contribution to V_SCF, i.e.
    !        [V_{SCF}^{HF}]_ab = -2 c_HF (ad|bc) P^{dc}
    !      where the response density matrix P is contracted with electron
    !      repulsion integrals (ERI) (ad|bc) (with electron indices
    !      (11|22)). Here, a, b, c and d are shorthand for alpha, beta, gamma
    !      and delta, labelling conduction NGWFs alpha, delta and valence
    !      NGWFs beta, gamma.
    !
    !      When evaluating the exchange (X) matrix for ground state
    !      calculations, the symmetry of the ERIs and density kernels involved
    !      allows a balanced density fitting approximation of the matrix
    !      elements by averaging X with its transpose. This is equivalent to
    !      averaging over matrix elements in which ERIs are approximated by
    !      fitting the bra and ket NGWF products.
    !
    !      For the V_{SCF}^{HF} matrix in LR-TDDFT, averaging with the transpose
    !      cannot be used to obtain a balanced density-fitting approximation.
    !      This is because the ERIs are no longer symmetric wrt interchange
    !      of the c and d indices. Instead, we evaluate two full V_{SCF}^{HF}
    !      matrices with (i) val-val and (ii) cond/joint-cond/joint NGWF
    !      products fitted and average over these.

    ! JCW: Check for devel_code TDDFT:CALC_VHF=[T/F]:TDDFT (default .true.)
    calc_Vhf = utils_devel_code(.true.,'TDDFT','CALC_VHF',pub_devel_code)

    if (calc_Vhf.and.(pub_use_hfx.or.pub_use_activehfx)) then

       ! jcap: If we are using EMFT with HFX in the active region, we
       ! need to switch the functional here to get the right HFX
       ! fraction in hf_exchange_calculate
       if (pub_emft.and.pub_use_activehfx) call xc_embed_swap_functional(.true.)

       ! JCW: Need an additional temporary array for storing transpose of
       !      response density kernel
       allocate(temp_array2(pub_num_spins),stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF','temp_array2',ierr)

       ! JCW: For EMFT calculations, HFX is only available in the active region
       !      For non-EMFT calculations, pub_active_region == 1 by default
       isub = pub_active_region

       ! JCW: Create temporary SPAM3_EMBED structures with sparsity of pot_SCF
       !      to store Hartree-Fock exchange contribution to to SCF response
       !      potential matrix. pot_SCF has cond/joint NGWFs as row index and
       !      val NGWFs as column index (it has the structure of
       !      cond_rep%cross_overlap_tr).
       allocate(temp_pot_SCF1(pub_num_spins),stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF','temp_pot_SCF1',ierr)
       allocate(temp_pot_SCF2(pub_num_spins),stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF','temp_pot_SCF2',ierr)
       do is=1,pub_num_spins
          ! temp_pot_SCF1 has the structure of the transpose of pot_SCF
          !    --> row index: val NGWFs, col index: cond/joint NGWFs
          ! Holds output of hf_exchange_calculate for this structure
          !    --> (b)eta:    val NGWFs, (a)lpha:   cond/joint NGWFs
          call sparse_embed_create(temp_pot_SCF1(is),pot_SCF(is),trans=.true.)
          ! temp_pot_SCF2 has the structure of pot_SCF
          !   --> row index: cond/joint NGWFs, col index: val NGWFs
          ! Holds output of hf_exchange_calculate for this structure
          !    --> (b)eta:   cond/joint NGWFs, (a)lpha:   val NGWFs
          call sparse_embed_create(temp_pot_SCF2(is),pot_SCF(is))
       end do

       ! JCW: For calculation of HFX contribution to SCF response potential
       !      matrix, use response density kernel as denskern_cd, extracted
       !      to temp_array from SPAM3_EMBED instance response_denskern
       call sparse_embed_extract_from_array(temp_array,response_denskern,isub,isub)

       ! JCW: Obtain transpose of response density kernel
       !      P^{1}     : temp_array
       !      --> row index: c/j NGWFs, col index: val NGWFs
       !      [P^{1}]^T : temp_array2
       !      --> row index: val NGWFs, col index: c/j NGWFs
       do is=1,pub_num_spins
          call sparse_transpose_structure(temp_array2(is)%structure,temp_array(is))
          call sparse_create(temp_array2(is))
          call sparse_transpose(temp_array2(is),temp_array(is))
       end do

       ! JCW: EVALUATION OF V_{SCF}^{HF}
       !
       ! NGWF_REP and FUNC_BASIS instances:
       ! * bra_rep and basis_bra should be for the conduction NGWFs (or,
       !   more generally, the NGWF set representing the unoccupied states)
       ! * ket_rep and basis_ket should be for the valence NGWFs
       !
       ! We want to evaluate
       !    V_ab = -2 c_HF 0.5 [(ad|{bc}) + ({ad}|bc)] P^{dc}
       ! where {xy} indicates SW-fitting of NGWF product xy and the NGWF indices
       ! are drawn from the following NGWF sets
       !
       !  Index      Set
       !  -----   ----------
       !    a     cond/joint
       !    b     valence
       !    c     valence
       !    d     cond/joint
       !
       ! Since hf_exchange_calculate can only fit the ket of the ERI, |bc), we
       ! actually evaluate
       !    V = 0.5 [ V^1_ab + V^2_ab ]
       ! with
       !    V^1_ab = -2 c_HF (ad|{bc}) P^{dc}
       ! and
       !    V^2_ab = -2 c_HF ({ad}|bc) P^{dc}
       !           = -2 c_HF transpose(W)_ab
       !
       !  Index      Set
       !  -----   ----------
       !    a     cond/joint
       !    b     valence
       !    c     valence
       !    d     cond/joint
       !
       ! V^2_ab is evaluated via
       !    W_ab   = (ad|{bc}) P^{dc}
       ! with alternative NGWF index-set assignments to ensure only ket-fitting
       ! is needed.
       !
       !  Index      Set
       !  -----   ----------
       !    a     valence
       !    b     cond/joint
       !    c     cond/joint
       !    d     valence
       !
       ! Note the use of a transpose response density kernel P with c cond/joint
       ! and d valence to evaluate W.

       ! JCW: MATRIX INDEXING FOR hf_exchange MODULE
       ! hf_exchange_calculate expects the following matrices as arguments:
       !
       ! mat (SPAM3_EMBED, inout)
       ! denskern_ab (SPAM3, in)
       ! denskern_cd (SPAM3, in)
       !
       ! In ground-state HFX evaluation, attention does not need to be paid to
       ! the indexing of rows and columns, since all matrices are symmetric
       ! and have identical NGWF sets along rows and columns.
       !
       ! For evaluation of V_{SCF}^{HF} in LR-TDDFT, care must be taken to
       ! ensure that the matrices are conformable for multiplication, i.e.
       ! their row and column indices must be compatible.
       !
       ! Using the alpha-beta-gamma-delta convention of the basis_selector
       ! argument (see header documentation for hf_exchange_calculate),
       ! hf_exchange_calculate expects
       !
       !     Matrix        ROW index     COL index
       !  -----------      ---------     ---------
       !  mat              b             a
       !  denskern_ab      a             b
       !  denskern_cd      c             d

       ! JCW: EVALUATION OF V^1_ab
       !    V^1_ab = -2 c_HF (ad|{bc}) P^{dc}
       !
       !  Index      Set
       !  -----   ----------
       !    a     cond/joint
       !    b     valence
       !    c     valence
       !    d     cond/joint
       !
       !   Dummy arg       ROW index     COL index      Actual arg
       !  -----------      ---------     ---------    -------------
       !  mat              b (v)         a (c/j)      temp_pot_SCF1
       !  denskern_ab      a (c/j)       b (v)        temp_array
       !  denskern_cd      c (v)         d (c/j)      temp_array2
       !
       !  Transpose mat to produce V^1_ab which is conformable for addition to
       !  pot_SCF (row: c/j, col: v).
       !
       ! The result should already be scaled by the fraction of HF exchange
       ! in the chosen XC functional, c_HF, within hf_exchange_calculate.
       ! -----------------------------------------------------------------------
       ! JCW: @optimize
       !      The NGWF sets are not changing, but rather we are switching
       !      between two fixed sets (val and joint). Since we only have one set
       !      of hash tables in hfxstate, we need to refresh these. This could
       !      possibly be optimized by
       !      1. Allowing NGWF expansions to persist in NGWF_REP instances,
       !         rather than re-expanding in hf_exchange_dkn_indep_stage
       !      2. Having hash tables for more than one NGWF product in hfxstate
       !      For now we will always purge the hash tables and re-expand the
       !      NGWF products in SWs
       !
       ! JCW: The SW_EX representing the val-val expansion coeffs is stored in
       !      ket_rep%swexes(REP_SWEX_HFX)
       call hf_exchange_dkn_indep_stage(hfxstate, mdl, isub, &
            mdl%regions(isub)%par, ket_rep, basis_ket(isub), &
            rep2 = bra_rep, ngwf_basis2 = basis_bra(isub), &
            basis_selector = [ 2,1,1,2,REP_SWEX_HFX ])
            ! basis for alpha: 2 (joint)
            ! basis for beta:  1 (val)
            ! basis for gamma: 1 (val)
            ! basis for delta: 2 (joint)
            ! swex selector: HFX (of val_rep)
            ! (corresponds to ket_basis_selector = (/1,1,REP_SWEX_HFX/)
            ! which is passed to swx_expand_my_ngwf_pairs inside
            ! hf_exchange_dkn_indep_stage)
       ! -----------------------------------------------------------------------

       ! Call 1, where rep (from which SW expansion of NGWFs is
       ! obtained) is valence NGWF rep (ket_rep)
       call hf_exchange_calculate(hfxstate, temp_pot_SCF1, ket_rep, &
            isub, temp_array, temp_array2, basis_ket(isub), &
            mdl%fftbox, mdl%cell, mdl%regions(isub)%elements, .false., &
            rep2 = bra_rep, &
            ngwf_basis2 = basis_bra(isub), &
            basis_selector = (/2,1,1,2,REP_SWEX_HFX/), &
            !                  ^ qoh's alpha-beta-gamma-delta conv'n
            !                    For X_ab = (ad|bg) P^{dg}
            !                    (a)lpha: conduction (bra_rep)
            !                    (b)eta:  valence (ket_rep)
            !                    (g)amma: valence (ket_rep)
            !                    (d)elta: conduction (bra_rep)
            !                    REP_SWEX_HFX: val-val from ket_rep
            symmetrise_xmatrix = .false. )
            ! ^ don't symmetrise by averaging the X matrix with its
            !   transpose, as this is not valid when the alpha-beta /
            !   gamma-delta NGWF pairs are comprised of mixtures of
            !   functions from different NGWF bases

       ! Get transpose of temp_pot_SCF1, which can then be added to pot_SCF
       do is=1,pub_num_spins
          call sparse_embed_transpose(temp_pot_SCF2(is),temp_pot_SCF1(is))
       end do

       ! Add val-val-fitted part of HF-contribution to pot_SCF
       do is=1, pub_num_spins
          call sparse_embed_axpy(pot_SCF(is), temp_pot_SCF2(is), &
               0.5_DP*Vhf_scale_factor)
               ! ^--- scaled by 0.5 since we are averaging over val-val-fitted
               !      and cond/joint-cond/joint-fitted matrices
       end do

       ! Clearing temp_pot_SCF2 should not be necessary, since the result
       ! of hf_exchange_calculate should fully overwrite the contents of the
       ! matrix when it is sparse_copy'd into temp_pot_SCF2.
       ! For debugging: Destroying and recreating temp_pot_SCF2 should
       ! absolutely ensure that all elements are zeroed before the second call
       ! to hf_exchange_calculate
       !do is=1, pub_num_spins
       !   call sparse_embed_destroy(temp_pot_SCF2(is))
       !   call sparse_embed_create(temp_pot_SCF2(is),pot_SCF(is))
       !end do

       ! JCW: EVALUATION OF V^2_ab
       !    V^2_ab = -2 c_HF ({ad}|bc) P^{dc}
       !           = -2 c_HF transpose(W)_ab
       !
       !    W_ab   = (ad|{bc}) P^{dc}
       !
       !  Index      Set
       !  -----   ----------
       !    a     valence
       !    b     cond/joint
       !    c     cond/joint
       !    d     valence
       !
       !   Dummy arg       ROW index     COL index      Actual arg
       !  -----------      ---------     ---------    -------------
       !  mat              b (c/j)       a (v)        temp_pot_SCF2
       !  denskern_ab      a (v)         b (c/j)      temp_array2
       !  denskern_cd      c (c/j)       d (v)        temp_array
       !
       !  W_ab is already transposed on output of hf_exchange_calculate, so is
       !  conformable for addition to pot_SCF (row: c/j, col: v).
       !
       ! The result should already be scaled by the fraction of HF exchange
       ! in the chosen XC functional, c_HF, within hf_exchange_calculate.
       ! -----------------------------------------------------------------------
       ! JCW: @optimize
       !      The NGWF sets are not changing, but rather we are switching
       !      between two fixed sets (val and joint). Since we only have one set
       !      of hash tables in hfxstate, we need to refresh these. This could
       !      possibly be optimized by
       !      1. Allowing NGWF expansions to persist in NGWF_REP instances,
       !         rather than re-expanding in hf_exchange_dkn_indep_stage
       !      2. Having hash tables for more than one NGWF product in hfxstate
       !      For now we will always purge the hash tables and re-expand the
       !      NGWF products in SWs
       !
       ! JCW: The SW_EX representing the joint-joint expansion coeffs is stored
       !      in bra_rep%swexes(REP_SWEX_HFX)
       call hf_exchange_dkn_indep_stage(hfxstate, mdl, isub, &
            mdl%regions(isub)%par, bra_rep, basis_bra(isub), &
            rep2 = ket_rep, ngwf_basis2 = basis_ket(isub), &
            basis_selector = [ 2,1,1,2,REP_SWEX_HFX ])
            ! basis for alpha: 2 (val)
            ! basis for beta:  1 (cond/joint)
            ! basis for gamma: 1 (cond/joint)
            ! basis for delta: 2 (val)
            ! swex selector: HFX (of cond/joint rep)
            ! (corresponds to ket_basis_selector = (/1,1,REP_SWEX_HFX/)
            ! which is passed to swx_expand_my_ngwf_pairs inside
            ! hf_exchange_dkn_indep_stage)
       ! -----------------------------------------------------------------------

       ! Call 2, where rep (from which SW expansion of NGWFs is
       ! obtained) is cond/joint NGWF rep (bra_rep)
       call hf_exchange_calculate(hfxstate, temp_pot_SCF2, bra_rep, &
            isub, temp_array2, temp_array, basis_bra(isub), &
            mdl%fftbox, mdl%cell, mdl%regions(isub)%elements, .false., &
            rep2 = ket_rep, &
            ngwf_basis2 = basis_ket(isub), &
            basis_selector = (/2,1,1,2,REP_SWEX_HFX/), &
            !                  ^ qoh's alpha-beta-gamma-delta conv'n
            !                    For X_ab = (ad|bg) P^{dg}
            !                    (a)lpha: valence (ket_rep)
            !                    (b)eta:  conduction (bra_rep)
            !                    (g)amma: conduction (bra_rep)
            !                    (d)elta: valence (ket_rep)
            !                    REP_SWEX_HFX: cond/joint-cond/joint from
            !                    bra_rep
            symmetrise_xmatrix = .false. )
            ! ^ don't symmetrise by averaging the X matrix with its
            !   transpose, as this is not valid when the alpha-beta /
            !   gamma-delta NGWF pairs are comprised of mixtures of
            !   functions from different NGWF bases

       ! Add cond/joint-cond/joint-fitted part of HF-contribution to pot_SCF
       do is=1, pub_num_spins
          call sparse_embed_axpy(pot_SCF(is), temp_pot_SCF2(is), &
               0.5_DP*Vhf_scale_factor)
               ! ^--- scaled by 0.5 since we are averaging over val-val-fitted
               !      and cond/joint-cond/joint-fitted matrices
       end do

       ! Destroy temporary SPAM3 instances for holding response density kernel
       do is=1,pub_num_spins
          call sparse_destroy(temp_array(is))
          call sparse_destroy(temp_array2(is))
       end do

       deallocate(temp_array2,stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF','temp_array2',ierr)

       ! Destroy temporary SPAM3_EMBED instances for holding V^1_ab and V^2_ab
       do is=1,pub_num_spins
          call sparse_embed_destroy(temp_pot_SCF1(is))
          call sparse_embed_destroy(temp_pot_SCF2(is))
       end do

       deallocate(temp_pot_SCF1,stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF','temp_pot_SCF1',ierr)

       deallocate(temp_pot_SCF2,stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF','temp_pot_SCF2',ierr)

       ! jcap: If we are using EMFT with HFX in the active region, we
       ! need to switch the functional back now
       if (pub_emft.and.pub_use_activehfx) call xc_embed_swap_functional(.false.)

    end if ! pub_use_hfx.or.pub_use_activehfx

    ! pot_SCF now contains soft real space grid contribution. For PAW, we need
    ! to augment it accordingly
    ! jcap: this needs to be seriously looked at for embedding
    if(pub_paw) then
       allocate(Dij(pub_num_spins), stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF','Dij',ierr)
       allocate(Dij_hartree(pub_num_spins), stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF','Dij_hartree',ierr)
       allocate(Dij_hat(pub_num_spins),stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF','Dij_hat',ierr)
       allocate(rho_ij(pub_num_spins), stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF', 'rho_ij',ierr)
       allocate(Dij_xc(pub_num_spins),stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF','Dij_xc',ierr)
       allocate(rho_0_ij(pub_num_spins), stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF', 'rho_0_ij',ierr)
       allocate(temp_array2(pub_num_spins),stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF','temp_array2',ierr)
       allocate(temp_array3(pub_num_spins),stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF','temp_array3',ierr)
       allocate(temp_array4(pub_num_spins),stat=ierr)
       call utils_alloc_check('linear_response_calc_SCF','temp_array4',ierr)

       do is=1, pub_num_spins
          Dij_hat(is)%structure='E'
          rho_ij(is)%structure='E'
          rho_0_ij(is)%structure='E'
          Dij_xc(is)%structure='E'
          Dij(is)%structure='E'
          Dij_hartree(is)%structure='E'
          call sparse_embed_create(Dij_hat(is))
          call sparse_embed_create(rho_ij(is))
          call sparse_embed_create(Dij(is))
          call sparse_embed_create(Dij_hartree(is))
          call sparse_embed_create(Dij_xc(is))
          call sparse_embed_create(rho_0_ij(is))
       enddo

       call sparse_embed_create(sp_overlap_Dij,bra_rep%sp_overlap,Dij(1))
       call sparse_embed_transpose_structure(ps_overlap%structure,&
            ket_rep%sp_overlap)
       call sparse_embed_create(ps_overlap, &
            iscmplx=ket_rep%sp_overlap%p%iscmplx)
       call sparse_embed_create(SCF_nonloc,pot_SCF(1))

       ! if nhat is not included in Vxc and this is a triplet calculation, then
       ! there is NO nhat contribution, since the hartree term vanishes for
       ! triplets

       if(.not. pub_nhat_in_xc) then
         if (.not. pub_lr_tddft_triplet) then
            if(pub_num_spins==1) then
               hartree_buffer=2.0_DP*hartree_buffer
            endif

            ! jcap: loop over regions - not sure this is the right
            ! thing to do, but it'll do for now
            do isub=1,mdl%nsub

               call sparse_embed_extract_from_array(temp_array,Dij_hat,isub,isub)
               call augmentation_screen_dij(temp_array,hartree_buffer,&
                    mdl%aug_box,mdl%cell,mdl%fine_grid,&
                    mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp)
               call sparse_embed_destroy_extracted_array(temp_array,Dij_hat,&
                    .true.,isub,isub)

            end do

            if(pub_num_spins==1) then
               hartree_buffer=1.0_DP/2.0_DP*hartree_buffer
            endif
         endif
      else

         ! jcap: loop over regions - not sure this is the right
         ! thing to do, but it'll do for now
         do isub=1,mdl%nsub

            call sparse_embed_extract_from_array(temp_array,Dij_hat,isub,isub)
            call augmentation_screen_dij(temp_array,SCF_buffer,&
                 mdl%aug_box,mdl%cell,mdl%fine_grid,&
                 mdl%regions(isub)%pseudo_sp,mdl%regions(isub)%paw_sp)
            call sparse_embed_destroy_extracted_array(temp_array,Dij_hat,&
                    .true.,isub,isub)

         end do
       endif
       ! calculate effective denskerns inside augmentation spheres
       ! jcap: need temporary arrays to pass in
       call sparse_embed_extract_from_array(temp_array,rho_ij)
       call sparse_embed_extract_from_array(temp_array2,response_denskern)
       call aug_projector_denskern(temp_array,temp_array2,&
            bra_rep%sp_overlap%p,ket_rep%sp_overlap%p)
       call sparse_embed_destroy_extracted_array(temp_array,rho_ij,&
            .true.,1,1)
       call sparse_embed_destroy_extracted_array(temp_array2)

       call sparse_embed_extract_from_array(temp_array,rho_0_ij)
       call sparse_embed_extract_from_array(temp_array2,val_denskern)
       call aug_projector_denskern(temp_array,temp_array2,ket_rep%sp_overlap%p)
       call sparse_embed_destroy_extracted_array(temp_array,rho_0_ij,&
            .true.)
       call sparse_embed_destroy_extracted_array(temp_array2)

       ! calculate dij_xc and dij_hartree (dij_hartree only necessary
       ! if this is not a triplet calculation)
       if(pub_lr_tddft_triplet) then
          ! jcap: loop over regions - not sure this is the right
          ! thing to do, but it'll do for now
          do isub=1,mdl%nsub

             call sparse_embed_extract_from_array(temp_array,Dij_xc,isub,isub)
             call sparse_embed_extract_from_array(temp_array2,rho_ij,isub,isub)
             call sparse_embed_extract_from_array(temp_array3,rho_0_ij,isub,isub)
             call paw_dij_xc_tddft(temp_array,temp_array2,temp_array3,&
                  mdl%regions(isub)%paw_sp,epsilon_val,triplet=.true.)
             call sparse_embed_destroy_extracted_array(temp_array,Dij_xc,&
                  .true.,isub,isub)
             call sparse_embed_destroy_extracted_array(temp_array2)
             call sparse_embed_destroy_extracted_array(temp_array3)

          end do
       else
          ! jcap: loop over regions - not sure this is the right
          ! thing to do, but it'll do for now
          do isub=1,mdl%nsub

             call sparse_embed_extract_from_array(temp_array,Dij_xc,isub,isub)
             call sparse_embed_extract_from_array(temp_array2,rho_ij,isub,isub)
             call sparse_embed_extract_from_array(temp_array3,rho_0_ij,isub,isub)
             call sparse_embed_extract_from_array(temp_array4,Dij_hartree,isub,isub)
             call paw_dij_hartree(temp_array4,temp_array2,&
                  mdl%regions(isub)%paw_sp)
             call paw_dij_xc_tddft(temp_array,temp_array2,temp_array3,&
                  mdl%regions(isub)%paw_sp,epsilon_val,triplet=.false.)
             call sparse_embed_destroy_extracted_array(temp_array,Dij_xc,&
                  .true.,isub,isub)
             call sparse_embed_destroy_extracted_array(temp_array2)
             call sparse_embed_destroy_extracted_array(temp_array3)
             call sparse_embed_destroy_extracted_array(temp_array4,Dij_hartree,&
                  .true.,isub,isub)

          end do
          ! deal with spin degeneracy if necessary
          if(pub_num_spins==1) then
             call sparse_embed_scale(Dij_hartree(1),2.0_DP)
          endif
       endif

       ! now form the correct Dij term. For triplets that do not contain
       ! nhat_in_xc, this only consists of the xc term. For other triplets,
       ! it also includes Dij_hat. For singlets, it includes Dij_hartree,
       ! Dij_hat and Dij_xc.
       if(pub_lr_tddft_triplet) then
          do is=1, pub_num_spins
             if(pub_nhat_in_xc) then
                call sparse_embed_copy(Dij(is),Dij_hat(is))
                call sparse_embed_axpy(Dij(is),Dij_xc(is),1.0_DP)
             else
                call sparse_embed_copy(Dij(is),Dij_xc(is))
             endif
          enddo
       else
          do is=1, pub_num_spins
             ! sum over spins already performed inside Dij_hartree
             ! routine (if spin=2).
             call sparse_embed_copy(Dij(is),Dij_hartree(is))
             call sparse_embed_axpy(Dij(is),Dij_hat(is),1.0_DP)
             call sparse_embed_axpy(Dij(is),Dij_xc(is),1.0_DP)
          enddo
       endif

       ! NOW CALCULATE NONLOCAL MATRIX
       call sparse_embed_transpose(ps_overlap,ket_rep%sp_overlap)
       do is=1, pub_num_spins
          call sparse_embed_product(sp_overlap_Dij,bra_rep%sp_overlap,Dij(is))
          call sparse_embed_product(SCF_nonloc,sp_overlap_Dij,ps_overlap)

          ! now add nonlocal part to previous matrix part
          call sparse_embed_axpy(pot_SCF(is), SCF_nonloc, 1.0_DP)
       enddo

       ! deallocate appropriate data structures
       do is=1, pub_num_spins
          call sparse_embed_destroy(Dij_hat(is))
          call sparse_embed_destroy(rho_ij(is))
          call sparse_embed_destroy(Dij(is))
          call sparse_embed_destroy(Dij_xc(is))
          call sparse_embed_destroy(rho_0_ij(is))
       enddo
       call sparse_embed_destroy(SCF_nonloc)
       call sparse_embed_destroy(sp_overlap_Dij)
       call sparse_embed_destroy(ps_overlap)
       deallocate(Dij, stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF','Dij',ierr)
       deallocate(Dij_hat,stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF','Dij_hat',ierr)
       deallocate(rho_ij, stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF', 'rho_ij',ierr)
       deallocate(Dij_xc,stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF','Dij_xc',ierr)
       deallocate(rho_0_ij, stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF', 'rho_0_ij',ierr)

       deallocate(nhat_den_grad,stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF',&
            'nhat_den_grad',ierr)
       deallocate(sub_nhat_den_grad,stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF',&
            'sub_nhat_den_grad',ierr)
       deallocate(nhat_den_grad_gs,stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF',&
            'nhat_den_grad_gs',ierr)
       deallocate(sub_nhat_den_grad_gs,stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF',&
            'sub_nhat_den_grad_gs',ierr)
       deallocate(temp_array2,stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF',&
            'temp_array2',ierr)
       deallocate(temp_array3,stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF',&
            'temp_array3',ierr)
       deallocate(temp_array4,stat=ierr)
       call utils_dealloc_check('linear_response_calc_SCF',&
            'temp_array4',ierr)

    endif


    deallocate(hartree_buffer, stat=ierr)
    call utils_dealloc_check('linear_response_calc_SCF','hartree_buffer',ierr)
    deallocate(sub_hartree_buffer, stat=ierr)
    call utils_dealloc_check('linear_response_calc_SCF','sub_hartree_buffer',ierr)
    deallocate(SCF_buffer, stat=ierr)
    call utils_dealloc_check('linear_response_calc_SCF','SCF_buffer',ierr)
    deallocate(fxc_buffer, stat=ierr)
    call utils_dealloc_check('linear_response_calc_SCF','fxc_buffer',ierr)
    deallocate(temp_array,stat=ierr)
    call utils_dealloc_check('linear_response_calc_SCF','temp_array',ierr)

    call timer_clock('linear_response_calc_SCF', 2)

  end subroutine linear_response_calc_SCF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linear_response_RPA_operator(resulting_vec,vec,cond_rep,&
      cond_ngwf_basis,cond_denskern,kchc,val_rep,&
      val_ngwf_basis,val_denskern,hvkv,mdl,hfxstate,response_dens,&
      sub_response_dens,fxc,ground_state_dens,sub_ground_state_dens,&
      contravariant_vec,ground_state_hartree)

    !======================================================================!
    ! This subroutine calculates the full TDDFT operator acting on a two   !
    ! component tranisiton vector. Following Tsiper, the RPA energy can    !
    ! be written as (p(a-B)p+q(A+B)q)/2. The ingoing vector vec is a two-  !
    ! component vector containing (p,q). The outoging vector resulting_vec !
    ! contains [(A-B)p, (A+B)q]                                            !
    ! Modified for embedding by Joseph Prentice, July 2018                 !
    !======================================================================!

    use comms, only: comms_barrier
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_lr_tddft_preopt, pub_num_spins,&
         pub_multigrid_hartree
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product,sparse_embed_axpy, sparse_embed_copy, &
         sparse_embed_transpose, sparse_embed_transpose_structure, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use sparse, only: SPAM3
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: resulting_vec(2,pub_num_spins)
    type(SPAM3_EMBED), intent(in) :: vec(2,pub_num_spins) ! holds the p and q vectors,
      ! where p=X-Y, q=X+Y, X is the excitation vector and Y is the
      ! deexcitation vecortr
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
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,mdl%nsub)
    real(kind=DP), intent(inout) :: ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)
    type(SPAM3_EMBED), optional, intent(inout) :: contravariant_vec(2,pub_num_spins)
    real(kind=DP), optional,intent(in) :: ground_state_hartree(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, &
         mdl%nsub,mdl%nsub)

    ! local variables
    type(SPAM3_EMBED), allocatable, dimension(:) :: Vscf_q ! created in ph space
    type(SPAM3_EMBED) :: k1sv
    type(SPAM3_EMBED) :: SCFkv
    type(SPAM3_EMBED), allocatable, dimension(:) :: response_denskern
    type(SPAM3), allocatable, dimension(:) :: temp_dkn
    type(SPAM3_EMBED), allocatable, dimension(:) :: temp_contravariant_vec,&
         temp_contravariant_vec2
    integer :: ierr, is, isub, jsub
    type(SPAM3_EMBED) :: buffer_vector

    call timer_clock('linear_response_RPA_operator',1)

    ! allocate sparse matrices
    allocate(Vscf_q(pub_num_spins), stat=ierr)
    call utils_alloc_check('linear_response_RPA_operator','Vscf_q',ierr)
    allocate(temp_contravariant_vec(pub_num_spins), stat=ierr)
    call utils_alloc_check('linear_response_RPA_operator',&
         'temp_contravariant_vec',ierr)
    allocate(temp_contravariant_vec2(pub_num_spins), stat=ierr)
    call utils_alloc_check('linear_response_RPA_operator',&
         'temp_contravariant_vec2',ierr)
    allocate(response_denskern(pub_num_spins),stat=ierr)
    call utils_alloc_check('linear_response_RPA_operator',&
         'response_denskern',ierr)
    allocate(temp_dkn(pub_num_spins),stat=ierr)
    call utils_alloc_check('linear_response_RPA_operator',&
         'temp_dkn',ierr)
    do is=1, pub_num_spins
       call sparse_embed_create(response_denskern(is),vec(1,is))
       call sparse_embed_create(temp_contravariant_vec(is),vec(1,is))
       call sparse_embed_create(temp_contravariant_vec2(is),vec(1,is))
    enddo
    do is=1, pub_num_spins
       call sparse_embed_create(Vscf_q(is),cond_rep%cross_overlap_tr)
    enddo
    call sparse_embed_create(k1sv, vec(1,1), val_rep%overlap)
    call sparse_embed_create(SCFkv, Vscf_q(1), val_denskern(1))
    call sparse_embed_create(buffer_vector, vec(1,1))

    ! comput KS eigenvalue difference first
    ! start with p component
    do is=1,pub_num_spins
       call sparse_embed_product(temp_contravariant_vec(is), kchc(is), vec(1,is))
       call sparse_embed_product(buffer_vector, vec(1,is), hvkv(is))
       call sparse_embed_axpy(temp_contravariant_vec(is), buffer_vector, -1.0_DP)
    enddo

    ! now do q component
    do is=1,pub_num_spins
       call sparse_embed_product(temp_contravariant_vec2(is), kchc(is), vec(2,is))
       call sparse_embed_product(buffer_vector, vec(2,is), hvkv(is))
       call sparse_embed_axpy(temp_contravariant_vec2(is), buffer_vector, -1.0_DP)
    enddo

    ! construct Vscf_q. If no hybrid functionals are used, the operator
    ! (A-B) just corresponds to the diagonal part of the TDDFT operator, so
    ! that no Vscf needs to be calculated for vec(1).
    if (.not. pub_lr_tddft_preopt) then
       ! start with xcomponent (hole-particle)
       do is=1,pub_num_spins
          call sparse_embed_copy(response_denskern(is),vec(2,is))
       end do

       ! jcap: loop over regions and sum up densities
       response_dens = 0.0_DP
       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub

             call sparse_embed_extract_from_array(temp_dkn,response_denskern,&
                  isub,jsub)
             call density_on_grid(sub_response_dens(:,:,:,:,isub,jsub),&
                  mdl%fine_grid,mdl%dbl_grid, &
                  mdl%cell, mdl%fftbox, temp_dkn,&
                  cond_rep%ngwf_cross_overlap_tr%m(isub,jsub), &
                  cond_rep%ngwfs_on_grid(isub), cond_ngwf_basis(isub), &
                  val_rep%ngwfs_on_grid(jsub), val_ngwf_basis(jsub))
             call sparse_embed_destroy_extracted_array(temp_dkn)

             response_dens=response_dens+sub_response_dens(:,:,:,:,isub,jsub)

          end do
       end do

       ! Calculate the selfconsistent field
       if(pub_multigrid_hartree) then
          call linear_response_calc_SCF(Vscf_q,fxc,response_dens,&
               sub_response_dens, response_denskern, &
               cond_ngwf_basis, cond_rep, val_ngwf_basis,&
               val_rep, ground_state_dens, sub_ground_state_dens, &
               val_denskern,mdl,hfxstate,ground_state_hartree)
       else
          call linear_response_calc_SCF(Vscf_q,fxc,response_dens,&
               sub_response_dens, response_denskern, &
               cond_ngwf_basis, cond_rep, val_ngwf_basis,&
               val_rep, ground_state_dens, sub_ground_state_dens, &
               val_denskern,mdl,hfxstate)
       endif

       do is=1,pub_num_spins
          call sparse_embed_product(SCFkv, Vscf_q(is), val_denskern(is))
          call sparse_embed_product(buffer_vector, cond_denskern(is), SCFkv)
          ! add x_component of selfconsistent field to ycomponent of diagonal term
          ! note that we need to calculate (A+B)q=(D+4*Vscf_q), where D is the
          ! diagonal KS part
          call sparse_embed_axpy(temp_contravariant_vec2(is), buffer_vector, 2.0_DP)
       enddo

    endif

    if(present(contravariant_vec)) then
      do is=1, pub_num_spins
         call sparse_embed_copy(contravariant_vec(1,is),temp_contravariant_vec(is))
         call sparse_embed_copy(contravariant_vec(2,is), temp_contravariant_vec2(is))
      enddo
    endif

    ! now temp_contravariant_vec and temp_contravariant_vec2 hold the correct
    ! quantities for the p and q vectors. Convert
    ! into covariant quantities and return.

    do is=1,pub_num_spins
       call sparse_embed_product(K1sv, temp_contravariant_vec(is), val_rep%overlap)
       call sparse_embed_product(resulting_vec(1,is),cond_rep%overlap, K1sv)
       call sparse_embed_product(K1sv, temp_contravariant_vec2(is), val_rep%overlap)
       call sparse_embed_product(resulting_vec(2,is),cond_rep%overlap, K1sv)
    enddo

    ! Done. Only thing left to do is deallocating data space
    do is=1,pub_num_spins
       call sparse_embed_destroy(Vscf_q(is))
       call sparse_embed_destroy(temp_contravariant_vec(is))
       call sparse_embed_destroy(temp_contravariant_vec2(is))
    enddo

    deallocate(response_denskern,stat=ierr)
    call utils_dealloc_check('linear_response_RPA_operator',&
         'response_denskern',ierr)
    deallocate(temp_dkn,stat=ierr)
    call utils_dealloc_check('linear_response_RPA_operator',&
         'temp_dkn',ierr)
    deallocate(Vscf_q,stat=ierr)
    call utils_dealloc_check('linear_response_RPA_operator','Vscf_q',ierr)
    deallocate(temp_contravariant_vec,stat=ierr)
    call utils_dealloc_check('linear_response_RPA_operator',&
     'temp_contravariant_vec',ierr)
    deallocate(temp_contravariant_vec2,stat=ierr)
    call utils_dealloc_check('linear_response_RPA_operator',&
     'temp_contravariant_vec2',ierr)
    call sparse_embed_destroy(buffer_vector)
    call sparse_embed_destroy(k1sv)
    call sparse_embed_destroy(SCFkv)

    call timer_clock('linear_response_RPA_operator',2)

  end subroutine linear_response_RPA_operator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linear_response_operator(resulting_vec, vec, cond_rep, &
       cond_ngwf_basis,cond_denskern, kchc, val_rep, &
       val_ngwf_basis, val_denskern, hvkv, mdl, hfxstate, response_dens, &
       sub_response_dens, fxc, ground_state_dens, sub_ground_state_dens, &
       contravariant_vec, ground_state_hartree)

    !=======================================================================!
    ! This subroutine calculates the result of the total linear-response    !
    ! operator acting on a trial response_density kernel 'vec' and returns  !
    ! the result of the action as 'resulting_vec'. Note that the total      !
    ! action of the linear response operator results is a covariant vector. !
    ! Optionally, a contravariant version of the result is also returned in !
    ! contravariant_vec.                                                    !
    ! Modified for embedding by Joseph Prentice, July 2018                  !
    !=======================================================================!

    use comms, only: comms_barrier
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_lr_tddft_preopt, pub_num_spins,&
         pub_multigrid_hartree
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_axpy, sparse_embed_copy, &
         sparse_embed_transpose, sparse_embed_transpose_structure, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use sparse, only: SPAM3
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: resulting_vec(:) ! jd: pub_num_spins. Leave ':' be.
    type(SPAM3_EMBED), intent(in) :: vec(pub_num_spins)
    type(NGWF_REP), intent(inout) :: cond_rep
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis(:)
    type(SPAM3_EMBED), intent(in) :: kchc(pub_num_spins) ! conduction density kernel * cond_ham
    type(SPAM3_EMBED), intent(in) :: cond_denskern(pub_num_spins)
    type(NGWF_REP), intent(inout) :: val_rep
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis(:)
    type(SPAM3_EMBED), intent(in) :: hvkv(pub_num_spins) ! val_ham*val_denskern
    type(SPAM3_EMBED), intent(in) :: val_denskern(pub_num_spins)
    type(MODEL), intent(in) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    real(kind=DP), intent(inout) :: response_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_response_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)
    real(kind=DP), intent(in) :: fxc(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,mdl%nsub)
    real(kind=DP), intent(inout) :: ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(inout) :: sub_ground_state_dens(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins,&
         mdl%nsub,mdl%nsub)
    type(SPAM3_EMBED), optional, intent(inout) :: contravariant_vec(:) ! jd: pub_num_spins. Leave ':' be.
    real(kind=DP), optional,intent(in) :: ground_state_hartree(mdl%fine_grid%ld1,&
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins, &
         mdl%nsub, mdl%nsub)

    ! Local Variables
    type(SPAM3_EMBED) :: k1sv, temp_vec
    type(SPAM3_EMBED), allocatable, dimension(:) :: temp_contravariant_vec
    type(SPAM3_EMBED), allocatable, dimension(:) :: SCF_field
    type(SPAM3_EMBED) :: SCFkv
    type(SPAM3), allocatable, dimension(:) :: response_denskern
    integer :: ierr, is, isub, jsub

    call timer_clock('linear_response_operator', 1)

    ! Data allocation
    allocate(temp_contravariant_vec(pub_num_spins), stat=ierr)
    call utils_alloc_check('linear_response_operator',&
      'temp_contravariant_vec',ierr)
    allocate(SCF_field(pub_num_spins), stat=ierr)
    call utils_alloc_check('linear_response_operator','SCF_field',ierr)
    do is=1,pub_num_spins
       call sparse_embed_create(SCF_field(is),cond_rep%cross_overlap_tr)
       call sparse_embed_create(temp_contravariant_vec(is),vec(is))
    enddo
    call sparse_embed_create(temp_vec, vec(1))
    call sparse_embed_create(k1sv, vec(1), val_rep%overlap)

    allocate(response_denskern(pub_num_spins), stat=ierr)
    call utils_alloc_check('linear_response_operator', &
         'response_denskern', ierr)
    call sparse_embed_create(SCFkv, SCF_field(1), val_denskern(1))

    ! calculate first part of the action (ie. the conduction
    ! valence Kohn-Sham eigenvalue difference
    do is=1, pub_num_spins
       call sparse_embed_product(temp_contravariant_vec(is), kchc(is), vec(is))
       call sparse_embed_product(temp_vec, vec(is), hvkv(is))
       call sparse_embed_axpy(temp_contravariant_vec(is), temp_vec, -1.0_DP)
    enddo

    ! Calculate coupling term. Only if we do not do a preoptimisation
    if (.not. pub_lr_tddft_preopt) then

       ! jcap: loop over regions and sum up densities
       response_dens = 0.0_DP
       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub

             call sparse_embed_extract_from_array(response_denskern,vec,&
                  isub,jsub)
             call density_on_grid(sub_response_dens(:,:,:,:,isub,jsub),&
                  mdl%fine_grid,mdl%dbl_grid, &
                  mdl%cell, mdl%fftbox, response_denskern,&
                  cond_rep%ngwf_cross_overlap_tr%m(isub,jsub), &
                  cond_rep%ngwfs_on_grid(isub), cond_ngwf_basis(isub), &
                  val_rep%ngwfs_on_grid(jsub), val_ngwf_basis(jsub))
             call sparse_embed_destroy_extracted_array(response_denskern)

             response_dens=response_dens+sub_response_dens(:,:,:,:,isub,jsub)

          end do
       end do

       ! Calculate the selfconsistent field
       ! jcap: replace response_denskern with vec, since that's
       ! where we copied it from
       if(pub_multigrid_hartree) then
          call linear_response_calc_SCF(SCF_field,fxc,response_dens,&
               sub_response_dens, vec, & !response_denskern, &
               cond_ngwf_basis, cond_rep, val_ngwf_basis,&
               val_rep, ground_state_dens, sub_ground_state_dens, &
               val_denskern,mdl,hfxstate,ground_state_hartree)
       else
          call linear_response_calc_SCF(SCF_field,fxc,response_dens,&
               sub_response_dens, vec, & !response_denskern, &
               cond_ngwf_basis, cond_rep, val_ngwf_basis,&
               val_rep, ground_state_dens, sub_ground_state_dens, &
               val_denskern,mdl,hfxstate)
       endif

       do is=1,pub_num_spins
          call sparse_embed_product(SCFkv, SCF_field(is), val_denskern(is))
          call sparse_embed_product(temp_vec, cond_denskern(is), SCFkv)
          ! add selfconsistent field to eigenvalue contributions
          call sparse_embed_axpy(temp_contravariant_vec(is), temp_vec, 1.0_DP)
       enddo
    endif

    ! return the tensorially correct g_vec as well
    if (present(contravariant_vec)) then
       do is=1,pub_num_spins
          call sparse_embed_copy(contravariant_vec(is), temp_contravariant_vec(is))
       enddo
    endif

    do is=1,pub_num_spins
       call sparse_embed_product(K1sv, temp_contravariant_vec(is), val_rep%overlap)
       call sparse_embed_product(resulting_vec(is),cond_rep%overlap, K1sv)
    enddo

    ! deallocate data
    call sparse_embed_destroy(temp_vec)
    call sparse_embed_destroy(k1sv)
    call sparse_embed_destroy(SCFkv)
    do is=1, pub_num_spins
       call sparse_embed_destroy(SCF_field(is))
       call sparse_embed_destroy(temp_contravariant_vec(is))
    enddo
    deallocate(response_denskern, stat=ierr)
    call utils_dealloc_check('linear_response_operator',&
         'response_denskern', ierr)
    deallocate(SCF_field, stat=ierr)
    call utils_dealloc_check('linear_response_operator',&
         'SCF_field', ierr)
    deallocate(temp_contravariant_vec, stat=ierr)
    call utils_dealloc_check('linear_response_operator',&
         'temp_contravariant_vec', ierr)

    call timer_clock('linear_response_operator', 2)

  end subroutine linear_response_operator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module linear_response
