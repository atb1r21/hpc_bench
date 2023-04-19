! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                           Hartree module                       !
!----------------------------------------------------------------!
! Originally written by Chris-Kriton Skylaris in June 2000.      !
! Rewritten by Peter Haynes, 28/6/04                             !
! Updated for variable grids by Nicholas Hine, June 2010.        !
! Further update by Jacek Dziedzic, 10/12/2010                   !
!================================================================!

module hartree

  use constants, only: DP
  implicit none

  private

  public :: hartree_on_grid
  public :: hartree_via_multigrid

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hartree_on_grid(hartree_potential,density,grid,cell)

    use cell_grid, only: GRID_INFO
    use constants, only: DP, UP, DN, PI
    use fourier, only: fourier_apply_cell_forward, fourier_apply_cell_backward
    use pbc_corrections, only: pbc_corr_hartree
    use rundat, only: pub_mt_cutoff, pub_num_spins
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    ! Grid definition
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    ! Density in real space:
    real(kind=DP), intent(inout) :: density(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    ! Potential in real space:
    real(kind=DP), intent(out) :: hartree_potential(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)

    ! Local variables
    integer :: ierr                                ! Error flag
    real(kind=DP), parameter :: fourpi = 4.0_DP*PI ! 4 pi
    complex(kind=DP), allocatable :: zwork(:,:,:)  ! Workspace

    ! Start timer
    call timer_clock('hartree_on_grid',1)

    ! Allocate workspace
    allocate(zwork(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('hartree_on_grid','zwork',ierr)

    if (pub_num_spins == 2) density(:,:,:,1) = &
         density(:,:,:,1) + density(:,:,:,2)

    ! Fourier transform density to reciprocal space
    call fourier_apply_cell_forward(density(:,:,:,1),zwork,grid)

    ! jd: Apply the Martyna-Tuckerman correction, if requested
    if(pub_mt_cutoff == 0.0_DP) then
       ! Multiply by 4pi/g^2 to obtain the usual Hartree potential
       zwork = zwork * grid%coulomb_recip(:,:,:) * fourpi
    else
       ! jd: Use the PBC-corrected formula
       call pbc_corr_hartree(zwork,grid,cell)
    end if

    ! Fourier transform potential to real space
    call fourier_apply_cell_backward(hartree_potential(:,:,:,1),zwork,grid)

    if (pub_num_spins == 2) then
       hartree_potential(:,:,:,DN) = hartree_potential(:,:,:,UP)
       density(:,:,:,1) = density(:,:,:,1) - density(:,:,:,2)
    end if

    ! Deallocate workspace
    deallocate(zwork,stat=ierr)
    call utils_dealloc_check('hartree_on_grid','zwork',ierr)

    ! Stop timer
    call timer_clock('hartree_on_grid',2)

  end subroutine hartree_on_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hartree_via_multigrid(phi, rho_elec, grid, cell, &
      E_Hartree, solvation_terms, uniform_eps, elements, no_dump)
    !=========================================================================!
    ! This is a front-end subroutine for the calculation of the Hartree energy!
    ! and Hartree potential using a multigrid solver. There are three,        !
    ! mutually-exclusive cases when this is in order:                         !
    !   a) A calculation without smeared ions, in vacuum.                     !
    !   b) A calculation, with smeared ions, in vacuum.                       !
    !   c) A calculation, with smeared ions, in solvent.                      !
    ! In a) the Hartree energy and potential are only due to electrons.       !
    ! In b) and c) the Hartree energy and potential are due to the whole      !
    !       molecule, i.e. electrons and smeared ions.                        !
    ! With b) and c), because of its molecular nature, E_Hartree cannot be    !
    ! calculated in the usual fashion, by integrating the Hartree potential   !
    ! with the electronic density, hence this subroutine calculates it on its !
    ! own. Futhermore, in the subcase of c) with a self-consistently changing !
    ! dielectric cavity, extra terms appear in the gradient, which are most   !
    ! conveniently treated by adding them to the potential. This subroutine   !
    ! does this.                                                              !
    ! For convenience, in c) this subroutine also calculates the cavitation   !
    ! energy and solute-solvent dispersion-repulsion energy, because it's     !
    ! easiest to do this at the same time. In cases a) and b) these terms     !
    ! will be zero.                                                           !
    ! Because the multigrid only works with the fine grid, there is no need   !
    ! to pass a grid argument to this subroutine.                             !
    !                                                                         !
    ! Open (Dirichlet) BCs and PBCs are supported.                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   phi,              intent=out, the calculated Hartree potential        !
    !   rho_elec,         intent=inout,  the input _electronic_ density       !
    !                     NB: inout, as it's temporarily summed over spins    !
    !   E_Hartree, (opt)  intent=out, calculated Hartree energy, if requested !
    !   solvation_terms, (opt, inout), if passed, its components will be      !
    !                                  populated with a breakdown of energy   !
    !                                  terms relevant in solvation.           !
    !   uniform_eps, (opt) intent=in,  the dielectric permittivity (scalar)   !
    !                                  to use for approximate BCs.            !
    !   elements, (opt) intent=in, passed on to implicit_solvent_hartree for  !
    !                              calculating polarization correction energy !
    !   no_dump, (opt)     intent=in,  logical parameter to avoid 3D dumps    !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 10/12/2010.                                  !
    ! Updated by Jacek Dziedzic, 20/11/2019 for Boltzmann solvation.          !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: DP, DN, UP
    use is_poisson, only: SOLVATION_ENERGY_TERMS, zero_solvation_terms
    use is_smeared_ions, only: smeared_ion_hartree
    use is_solvation, only: implicit_solvent_hartree, &
         implicit_solvent_energy_terms
    use multigrid_methods, only: multigrid_calculate_hartree, &
         multigrid_initialise
    use rundat, only: pub_is_smeared_ion_rep, pub_is_implicit_solvent, &
         pub_multigrid_hartree, pub_num_spins
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort
    use ion, only: ELEMENT

    implicit none

    ! jd: Arguments
    type(GRID_INFO), intent(in)                           :: grid
    type(CELL_INFO), intent(in)                           :: cell
    real(kind=DP), intent(out)                            :: phi(grid%ld1, &
         grid%ld2, grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(inout)                          :: rho_elec(grid%ld1,&
         grid%ld2, grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(out), optional                  :: E_Hartree
    type(SOLVATION_ENERGY_TERMS), intent(inout), optional :: solvation_terms
    real(kind=DP), intent(in), optional                   :: uniform_eps
    type(ELEMENT), intent(in), optional                   :: elements(:)
    logical, intent(in), optional                         :: no_dump

    ! jd: Local variables
    real(kind=DP) :: E_Hartree_loc
    type(SOLVATION_ENERGY_TERMS) :: loc_solvation_terms

    ! ------------------------------------------------------------------------

    if(.not. pub_is_implicit_solvent) then
       call zero_solvation_terms(loc_solvation_terms, zero_apolar = .true., &
            zero_vacuum = .true.)
    end if

    if(present(E_Hartree)) E_Hartree = 0D0

    ! jd: Use total electronic density, if we have two spins
    if (pub_num_spins == 2) rho_elec(:,:,:,UP) = &
         rho_elec(:,:,:,UP) + rho_elec(:,:,:,DN)

    ! jd: Need to initialise multigrid first, because the hartree routines
    !     dimension arrays with variables that are set up in the initialiser
    call multigrid_initialise(grid)

    ! jd: Case c)
    if(pub_is_implicit_solvent) then
       call implicit_solvent_hartree(phi, rho_elec(:,:,:,UP), &
            grid, cell, E_Hartree, loc_solvation_terms, elements, no_dump=no_dump)

    ! jd: Case b)
    else if(pub_is_smeared_ion_rep) then
       call smeared_ion_hartree(phi, rho_elec(:,:,:,UP), grid, cell, &
            E_Hartree_loc, no_dump=no_dump)
       if(present(E_Hartree)) E_Hartree = E_Hartree_loc

    ! jd: Case a)
    else if(pub_multigrid_hartree) then
       E_Hartree_loc = multigrid_calculate_hartree(phi, rho_elec(:,:,:,UP), &
            grid, cell, uniform_eps=uniform_eps, no_dump=no_dump)
       if(present(E_Hartree)) E_Hartree = E_Hartree_loc
    else
       call utils_abort('Internal error in hartree_via_multigrid.')
    end if

    if(.not. pub_is_implicit_solvent) then
       ! jd: Store the most recent in-vacuum Hartree energy in ugly module-wide
       !     state, so that it is accessible at the end of
       !     energy_and_force_calculate() for pretty-printing solvation
       !     energy components in auto solvation.
       implicit_solvent_energy_terms%E_elec_fixed_vac = E_Hartree_loc
    end if

    if(present(solvation_terms)) then
       solvation_terms = loc_solvation_terms
    end if

    ! jd: Undo the spin summation, update second spin component of potential
    if (pub_num_spins == 2) then
       phi(:,:,:,DN) = phi(:,:,:,UP)
       rho_elec(:,:,:,UP) = rho_elec(:,:,:,UP) - rho_elec(:,:,:,DN)
    end if

  end subroutine hartree_via_multigrid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module hartree

