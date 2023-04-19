! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!               Electron Localisation descriptors                !
!                                                                !
!                      Rebecca J. Clements,                      !
!      James C. Womack and Chris-Kriton Skylaris, July 2019      !
!                                                                !
!----------------------------------------------------------------!
! Electron localisation descriptors introduced by Becke and      !
! Edgecombe in 1990 and 2000, and applied to DFT by Savin et al. !
! in 1992.                                                       !
!================================================================!
module eld

  use constants, only: DP, SAFE_DIV_EPS, stdout, stderr
  use rundat, only: pub_debug_on_root

  implicit none

  public :: eld_on_grid
  private :: elf_on_grid
  private :: lol_on_grid

contains

  subroutine eld_on_grid(eld_fine, fine_grid, density_fine, ke_density_fine, &
       eld_str)

    !==========================================================================!
    ! This subroutine calculates the chosen electron localisation descriptor,  !
    ! either the Electron Localisation Function, ELF, or the Localised Orbital !
    ! Locator, LOL, specified by eld_str.                                      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! eld_fine (output) : The electron localisation descriptor on the fine     !
    !  grid of the simulation cell                                             !
    ! fine_grid (input) : GRID_INFO describing fine grid on which the electron !
    !  localisation descriptor is eventually constructed                       !
    ! density_fine (input) : The total charge density of the system on the     !
    !  fine grid of the simulation cell                                        !
    ! ke_density_fine (input) : The total kinetic energy density on the fine   !
    !  grid of the simulation cell                                             !
    ! eld_str (input) : The keywords "ELF" and "LOL" to be used to choose      !
    !  which electron localisation descriptor is to be calculated              !
    !--------------------------------------------------------------------------!
    ! Created by Rebecca J. Clements, 2019.                                    !
    !==========================================================================!

    use bibliography, only: bibliography_cite
    use cell_grid, only: GRID_INFO
    use constants, only: DP
    use rundat, only: pub_num_spins

    implicit none

    ! rjc: Arguments
    type(GRID_INFO), intent(in) :: fine_grid
    real(kind=DP), intent(out)  :: eld_fine(fine_grid%ld1, &
         fine_grid%ld2, fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in)  :: ke_density_fine(fine_grid%ld1, &
         fine_grid%ld2, fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in)  :: density_fine(fine_grid%ld1, &
         fine_grid%ld2, fine_grid%max_slabs12, pub_num_spins)
    character(len=*),intent(in) :: eld_str

    ! rjc: Select subroutine corresponding to eld_str value
    select case (eld_str)

    case ('ELF')
       call elf_on_grid(eld_fine, fine_grid, density_fine, ke_density_fine)

       ! rjc: cite reference
       call bibliography_cite('ELF')
       call bibliography_cite('ELF_DFT')

    case ('LOL')
       call lol_on_grid(eld_fine, fine_grid, density_fine, ke_density_fine)

       ! rjc: cite reference
       call bibliography_cite('LOL')

    end select

  end subroutine eld_on_grid

  subroutine elf_on_grid(eld_fine, fine_grid, density_fine, ke_density_fine)

    !==========================================================================!
    ! This subroutine calculates the Electron Localisation Function, ELF,      !
    ! introduced by Becke and Edgecombe                                        !
    ! (Becke, A. D. and Edgecombe, K. E., J. Chem. Phys., 1990, 92(9),5397-    !
    ! 5403) and later applied to DFT by Savin et al                            !
    ! (Savin, A., Jepsen, O., Flad, J., Andersen, O. K., Preuss, H. and        !
    ! von Schnering, H. G., Angew. Chem. Int. Ed. Engl., 1992, 31(2), 187-188) !
    ! This measures the probability of finding a second like-spin electron     !
    ! near the reference point. It takes the form:                             !
    !                                                                          !
    ! ELF = (1 + \chi_{\sigma}^{2})^{-1}, leading to values between 0 and 1.   !
    !                                                                          !
    ! where \chi = D_{\sigma} / D_{\sigma}^{0}                                 !
    !                                                                          !
    ! and D_{\sigma}^{0} = (3/5) * (6 * \pi^{2})^{2/3} * n_{\sigma}^{5/3},     !
    ! the uniform electron gas with spin-density equal to the local value of n !
    !                                                                          !
    ! and D_{\sigma} =                                                         !
    !    \tau_{\sigma} - (1/4) * |\grad(n_{\sigma})^{2}| / n_{\sigma}.         !
    !                                                                          !
    ! \tau_{\sigma}, the kinetic energy density, is defined without the        !
    ! standard factor of 1/2 by Becke                                          !
    ! (standard: Ernzerhof, M. & Scuseria, J. Chem. Phys. 111, 911-915 (1999)),!
    ! so here we double the kinetic energy density calculated in ONETEP.       !
    ! Becke's version is of the form:                                          !
    ! \tau =                                                                   !
    ! \sum_{\alpha\beta} K^{\alpha\beta}                                       !
    ! (\nabla\phi_{\alpha}).(\nabla\phi_{\beta})                               !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! eld_fine (output) : The electron localisation descriptor on the fine     !
    !  grid of the simulation cell                                             !
    ! fine_grid (input) : GRID_INFO describing fine grid on which ELF          !
    !  is eventually constructed                                               !
    ! density_fine (input) : The total charge density of the system on the     !
    !  fine grid of the simulation cell                                        !
    ! ke_density_fine (input) : The total kinetic energy density on the fine   !
    !  grid of the simulation cell                                             !
    !--------------------------------------------------------------------------!
    ! Created by Rebecca J. Clements, 2019.                                    !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: DP
    use rundat, only: pub_num_spins
    use services, only: services_gradient_on_grid
    use utils, only: utils_alloc_check, utils_dealloc_check

    ! rjc: Arguments
    type(GRID_INFO), intent(in) :: fine_grid
    real(kind=DP), intent(in)  :: ke_density_fine(fine_grid%ld1, &
         fine_grid%ld2, fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in)  :: density_fine(fine_grid%ld1, &
         fine_grid%ld2, fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(out) :: eld_fine(fine_grid%ld1, &
         fine_grid%ld2, fine_grid%max_slabs12, pub_num_spins)

    ! rjc: Local Variables
    integer :: ierr
    integer :: i1, i2, islab12, is ! Loop counters over grid points
    real(kind=DP), allocatable :: func_grad(:,:,:,:,:)
    real(kind=DP) :: w_ke_density
    real(kind=DP) :: D
    real(kind=DP) :: ratio
    real(kind=DP), parameter :: threshold = 1.0e-8_DP

    ! rjc: Allocate per-proc fine grid real space work array.
    allocate(func_grad(fine_grid%ld1,fine_grid%ld2,&
         fine_grid%max_slabs12,3,pub_num_spins),stat=ierr)
    call utils_alloc_check('elf_on_grid','func_grad',ierr)

    ! rjc: Loop over spins for services_gradient_on_grid
    do is=1,pub_num_spins
       call services_gradient_on_grid(func_grad(:,:,:,:,is),&
            density_fine(:,:,:,is),fine_grid)
    end do

    ! rjc: For square of density gradient, dot product func_grad with itself
    ! rjc: for each square and sum Cartesian components of func_grad
    func_grad(:,:,:,1,:) = &
         func_grad(:,:,:,1,:)*func_grad(:,:,:,1,:) + &
         func_grad(:,:,:,2,:)*func_grad(:,:,:,2,:) + &
         func_grad(:,:,:,3,:)*func_grad(:,:,:,3,:)
    ! ...func_grad(:,:,:,1,:) now contains |grad(rho)|^{2}

    ! rjc: Take care of num/max padding
    eld_fine = 0.0_DP

    ! rjc: Loop over spins for ELF formula
    do is=1,pub_num_spins
       do islab12=1,fine_grid%num_my_slabs12
          do i2=1,fine_grid%n2
             do i1=1,fine_grid%n1

                ! rjc: Apply threshold to charge density to avoid numerical
                ! rjc: errors
                if (abs(density_fine(i1,i2,islab12,is)).gt.threshold.and.&
                     abs(ke_density_lsda(density_fine(i1,i2,islab12,is))).gt.&
                     threshold) then

                   ! rjc: Calculate the von Weizsacker kinetic energy density
                   ! rjc: using previous jcw code
                   ! jcw: The von Weizsacker kinetic energy functional was
                   ! jcw: developed in the 1930s by Weizsacker, and is described
                   ! jcw: in many more recent papers, e.g. Acharya, P. K.,
                   ! jcw: Bartolotti, L. J., Sears, S. B. & Parr, R. G.,
                   ! jcw: PNAS 77, 6978-6982 (1980).
                   w_ke_density = &
                        0.25_DP * func_grad(i1,i2,islab12,1,is) &
                        / density_fine(i1,i2,islab12,is)

                   ! rjc: Calculate the Pauli kinetic energy density. Electron
                   ! rjc: localisation relates to the smallness of this
                   ! rjc: expression. Double kinetic energy density according
                   ! rjc: to Becke's formula.
                   D = 2.0_DP*ke_density_fine(i1,i2,islab12,is) - &
                        w_ke_density

                   ! rjc: Compare against a uniform electron gas with
                   ! rjc: spin-density equal to the local value of n_{\sigma}(r)
                   ratio = D/ke_density_lsda(density_fine(i1,i2,islab12,is))

                   ! rjc: Convert values to range 0 to 1
                   eld_fine(i1,i2,islab12,is) = 1.0_DP/(1.0_DP+ratio**2)

                else

                   eld_fine(i1,i2,islab12,is) = 0.0_DP

                end if
             end do
          end do
       end do
    end do

    ! rjc: Deallocate memory
    deallocate(func_grad,stat=ierr)
    call utils_dealloc_check('elf_on_grid','func_grad',ierr)

  end subroutine elf_on_grid


  subroutine lol_on_grid(eld_fine,fine_grid, density_fine, ke_density_fine)

    !==========================================================================!
    ! This subroutine calculates the Localised Orbital Locator, LOL,           !
    ! introduced by Becke and Schmider                                         !
    ! (Schmider, H. L. and Becke, A. D., J. Mol. Struct. (Theochem), 2000,     !
    ! 527,51-61)                                                               !
    ! which takes the form:                                                    !
    !                                                                          !
    ! LOL = t_{\sigma} / (1 + t_{\sigma}, leading to values between 0 and 1.   !
    !                                                                          !
    ! where t_{\sigma} = \tau_{\sigma}^{LSDA} / \tau_{\sigma}^{exact}          !
    !                                                                          !
    ! and \tau_{\sigma}^{LSDA} =                                               !
    ! (3/5) * (6 * \pi^{2})^{2/3} * n_\sigma^{5/3}, the uniform electron       !
    ! as for the ELF subroutine, with spin-density equal to the local value of !
    ! n,                                                                       !
    !                                                                          !
    ! and t_{\sigma}^{exact}, the exact kinetic energy density, is defined     !
    ! without the standard factor of 1/2 by Becke                              !
    ! (standard: Ernzerhof, M. & Scuseria, J. Chem. Phys. 111, 911-915 (1999)),!
    ! so here we double the kinetic energy density calculated in ONETEP.       !
    ! Becke's version is of the form:                                          !
    ! \tau =                                                                   !
    ! \sum_{\alpha\beta} K^{\alpha\beta}                                       !
    ! (\nabla\phi_{\alpha}).(\nabla\phi_{\beta})                               !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! eld_fine (output) : The electron localisation descriptor on the fine     !
    !  grid of the simulation cell                                             !
    ! fine_grid (input) : GRID_INFO describing fine grid on which LOL          !
    !  is eventually constructed                                               !
    ! density_fine (input) : The total charge density of the system on the     !
    !  fine grid of the simulation cell                                        !
    ! ke_density_fine (input) : The total kinetic energy density on the fine   !
    !  grid of the simulation cell                                             !
    !--------------------------------------------------------------------------!
    ! Created by Rebecca J. Clements, 2019.                                    !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use rundat, only: pub_num_spins

    implicit none

    ! rjc: Arguments
    type(GRID_INFO), intent(in) :: fine_grid
    real(kind=DP), intent(in)  :: ke_density_fine(fine_grid%ld1, &
         fine_grid%ld2, fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in)  :: density_fine(fine_grid%ld1, &
         fine_grid%ld2, fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(out)  :: eld_fine(fine_grid%ld1, &
         fine_grid%ld2, fine_grid%max_slabs12, pub_num_spins)

    ! rjc: Local Variables
    real(kind=DP) :: ratio
    integer :: i1, i2, islab12, is ! Loop counters over grid points
    real(kind=DP), parameter :: threshold = 1.0e-8_DP

    ! rjc: Take care of num/max padding
    eld_fine = 0.0_DP

    ! rjc: Get KE density LSDA function and use to calculate LOL density
    do is=1,pub_num_spins
       do islab12=1,fine_grid%num_my_slabs12
          do i2=1,fine_grid%n2
             do i1=1,fine_grid%n1

                ! rjc: Apply threshold to charge density to avoid numerical
                ! rjc: errors
                if (abs(density_fine(i1,i2,islab12,is)).gt.threshold.and.&
                     abs(ke_density_fine(i1,i2,islab12,is)).gt.threshold) then

                   !rjc: Becke's LOL formula. Double kinetic energy density
                   ratio = ke_density_lsda(density_fine(i1,i2,islab12,is)) &
                        / (2.0_DP*ke_density_fine(i1,i2,islab12,is))

                   ! rjc: Convert values to range 0 to 1
                   eld_fine(i1,i2,islab12,is) = ratio / (ratio + 1.0_DP)

                else

                   eld_fine(i1,i2,islab12,is) = 0.0_DP

                end if
             end do
          end do
       end do
    end do

  end subroutine lol_on_grid


  elemental function ke_density_lsda(den) result(ked_lsda)

    !==========================================================================!
    ! This function calculates the local spin density approximation            !
    ! \tau_{\sigma}^{LSDA} = (3/5) * (6 * \pi^{2})^{2/3} * n_\sigma^{5/3}.     !
    ! This serves as the uniform electron gas, with spin-density equal to the  !
    ! local value of n.                                                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! den (input) : The total charge density of the system on the fine grid of !
    !  the simulation cell                                                     !
    ! ked_lsda (output) : The local spin density approximation                 !
    !--------------------------------------------------------------------------!
    ! Created by Rebecca J. Clements, 2019.                                    !
    ! Reference:                                                               !
    ! Becke, A. D. and Edgecombe, K. E., J. Chem. Phys., 1990, 92(9),5397-5403 !
    !==========================================================================!

    use constants, only: DP, PI

    implicit none

    ! rjc: Arguments
    real(kind=DP), intent(in) :: den !density

    ! rjc: Parameters
    real(kind=DP), parameter :: prefactor = (3.0_DP/5.0_DP)*(6.0_DP*PI**2.0_DP)&
         **(2.0_DP/3.0_DP)

    ! rjc: Result
    real(kind=DP) :: ked_lsda

    ! rjc: Function formula
    ! jd: Fix to filter out negative densities
    if(den > SAFE_DIV_EPS) then
       ked_lsda = prefactor*den**(5.0_DP/3.0_DP)
    else
       ked_lsda = 0.0_DP
    end if

  end function ke_density_lsda

end module eld
