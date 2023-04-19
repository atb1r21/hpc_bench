! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!==============================================================================!
!                                                                              !
!  ONETEP exchange-correlation module: xc_mod.F90                              !
!                                                                              !
!  The subroutines in this file were written by Peter Haynes, February 2004    !
!                                                                              !
!  Theory of Condensed Matter, Cavendish Laboratory, University of Cambridge   !
!  Madingley Road, Cambridge CB3 0HE, UK                                       !
!                                                                              !
!  Addition of hybrid functionals, BLYP and XLYP, and major restructuring of   !
!  module by Quintin Hill in February/March 2009 under the supervision of      !
!  Chris-Kriton Skylaris.                                                      !
!  Note for hybrid functionals the Hartree-Fock exchange contribution to the   !
!  energy and potential is calculated and applied outside this subroutine.     !
!                                                                              !
!  Modified for TDDFT by David D. O'Regan in May 2009.                         !
!  Added support for libxc functionals, Nicholas Hine, July 2010.              !
!  Added support for radial xc calculations, Nicholas Hine, September 2010.    !
!------------------------------------------------------------------------------!
!                                                                              !
! The following exchange-correlation functionals (specified by the input file  !
! keyword XC_FUNCTIONAL) are available:                                        !
!                                                                              !
!  Local density approximation functionals:                                    !
!   CAPZ   - L(S)DA based on Ceperley-Alder Monte Carlo data parametrised by   !
!            Perdew & Zunger                                                   !
!   VWN    - L(S)DA based on parametrisation by Vosko, Wilk & Nusair           !
!   PW92   - L(S)DA based on Ceperley-Alder Monte Carlo data parametrised by   !
!            Perdew & Wang                                                     !
!                                                                              !
!  Generalised gradient approximation functionals:                             !
!   BLYP   - GG(S)A based on Becke 88 + LYP (Lee, Yang, Parr)                  !
!   PBE    - GG(S)A based on 1996 functional of Perdew, Burke Ernzerhof        !
!   REVPBE - GG(S)A based on 1998 revised PBE by Zhang & Yang                  !
!   RPBE   - GG(S)A based on 1999 RPBE by Hammer, Hansen & Norskov             !
!   PBESOL - GG(S)A based on 2008 functional of Perdew et al.                  !
!   PW91   - GG(S)A based on 1991 functional of Perdew & Wang                  !
!   XLYP   - GG(S)A based on 2004 functional of Xu and Goddard                 !
!   WC     - GG(S)A based on 2006 functional of Wu and Cohen                   !
!                                                                              !
!  Hybrid functionals:                                                         !
!   B1LYP  - hybrid based on 1997 functional of Adamo and Barone               !
!   B1PW91 - hybrid based on 1997 functional of Adamo and Barone               !
!   B3LYP  - hybrid based on 1993 functional of Becke using LYP instead of PW91!
!   B3PW91 - hybrid based on 1993 functional of Becke                          !
!   PBE0   - hybrid based on 1999 functional of Adamo, Cossi and Barone        !
!   X3LYP  - hybrid based on 2004 functional of Xu and Goddard                 !
!                                                                              !
!  Libxc functionals:                                                          !
!   LIBXC - specify actual functional id with the libxc_x/cfunc_id parameters  !
!                                                                              !
!  Non-functionals:                                                            !
!   HF     - only Hartree-Fock exchange (for solving Hartree-Fock equations)   !
!   NONE   - no exchange-correlation (for solving Hartree equations)           !
!                                                                              !
!   In addition, the following defaults apply:                                 !
!                                                                              !
!   LDA -> CAPZ                                                                !
!   GGA -> RPBE                                                                !
!                                                                              !
!==============================================================================!

module xc

  use constants, only: DP
#ifdef LIBXC
  use xc_f03_lib_m !! External dependency
#endif

  implicit none

  private

  ! Public subroutines
  public :: xc_energy_potential
  public :: xc_hfxinit
  public :: xc_init
  public :: xc_exit
  public :: xc_radial
  public :: xc_fxc_potential
  public :: xc_fxc_finite_diff
  !public :: xc_test
  public :: xc_gradients ! used by ke_density module
  public :: xc_emft_calculate
  public :: xc_embed_swap_functional
  public :: xc_fxc_potential_emft
  public :: xc_fxc_finite_diff_emft

  ! Parsed version of public "pub_xc_functional" variable
  character(len=80) :: functional

  ! JCW: Alternative functionals which can be used for the initial guess /
  ! JCW: radial densities.
  ! JCW: Functional to be used by xc_radial (can be different to pub_xc_functional),
  ! JCW: set in xc_init.
  character(len=80) :: radial_functional
  ! JCW: Functional to be used by xc_energy_potential when optional variable
  ! JCW: "initial" is true, set in xc_init.
  character(len=80) :: initial_functional

  ! Fraction of Hartree-Fock exchange
  real(kind=DP), save, public :: pub_hfxfraction

  ! Whether gradient corrections are required
  logical, save, public :: pub_xc_gradient_corrected
  ! JCW: Whether kinetic energy density is required
  !logical, save, public :: pub_xc_ke_density_required   <-- moved to rundat
  ! JCW: Whether Laplacian of density is required
  !logical, save, public :: pub_xc_lapl_density_required <-- moved to rundat

  ! Private internal variables used by the module
#ifdef LIBXC
  logical :: use_libxc
  logical :: libxc_gga  ! true if GGA or hybrid GGA, false otherwise
  logical :: libxc_mgga ! true if meta-GGA or hybrid meta-GGA, false otherwise
  logical :: libxc_nlc_vv10  ! true if Libxc functional requires non-local
                             ! correlation using the VV10 functional implemented
                             ! in nlxc module
  real(kind=DP) :: libxc_vv10_bval ! Set in xc_init, if required
  real(kind=DP) :: libxc_vv10_Cval ! Set in xc_init, if required
  type(xc_f03_func_t) :: libxc_func(2)
  type(xc_f03_func_info_t) :: libxc_info(2)
#endif


contains

  subroutine xc_fxc_potential(grid,fxc_pot_fine, density_fine,parallel)
    !=========================================================================!
    !Subroutine getting fxc in fine grid representation from the ground state !
    !density of the converged system. Returns fxc_pot_fine. Currently only    !
    !capable of doing ALDA in spin degenerate systems.                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !    density_fine (in) : Input density on the fine grid                   !
    !    grid (in)         : GRID_INFO grid definition                        !
    !    fxc_pot_fine (out) : Exchange-correlation kernel on grid             !
    !-------------------------------------------------------------------------!
    ! Written by Tim Zuehlsdorff, March 2013                                  !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_barrier, pub_on_root
    use constants, only: stdout
    use rundat, only: pub_tddft_xc_functional, pub_num_spins
    use xc_funcs, only: xc_fxc_capz_point

    implicit none

    type (GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: fxc_pot_fine(grid%ld1,&
         grid%ld2,grid%max_slabs12,2)
    real(kind=DP), intent(inout) :: density_fine(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    logical, optional, intent(in) :: parallel

    integer :: islab12, i1, i2 ! loop counters
    real(kind=DP) :: xc_energy
    real(kind=DP) :: pot1, pot2, den
    logical :: loc_parallel

    loc_parallel = .true.
    if(present(parallel)) then
       loc_parallel = parallel
    end if

    if (functional /= 'LDA') then
       if(pub_on_root) write(stdout, '(a)') 'Overriding xc_func &
       &to LDA, since nothing else is implemented for tddft yet'
    endif

    ! tjz07: Start fxc calculation on a double grid. Reset
    ! tjz07: xc_energy.
    xc_energy=0.0_DP

    ! now, calculate fxc by calling the appropriate function for
    ! getting fxc at each point
    fxc_pot_fine=0.0_DP
    a3: do islab12=1, grid%num_my_slabs12
       a2: do i2=1, grid%n2
            a1: do i1=1, grid%n1
              ! get f_xc at that point
              den=density_fine(i1,i2,islab12,1)
              call xc_fxc_capz_point(den, pot1, pot2)

              fxc_pot_fine(i1,i2,islab12,1)=pot1
              fxc_pot_fine(i1,i2,islab12,2)=pot2

          enddo a1
       enddo a2
    enddo a3

    if(loc_parallel) then
       call comms_barrier
    end if

  end subroutine xc_fxc_potential


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_fxc_finite_diff(density_fine,fxc_pot_fine,&
       grid,cell,ground_state_dens,nhat_dim,eff_nhat_den,nhat_den_grad_gs,&
       parallel)

    !===========================================================================!
    ! Given the electronic ground state density and the perturbed density on    !
    ! fine grid, this routine calculates fxc*rho_perturbed using only Vxc       !
    ! and a finite difference technique. This prevents us from having to code   !
    ! up the actual second derivatives of the exchange correlation energy,      !
    ! which for GGA codes is cumbersome and numerically unstable.               !
    !---------------------------------------------------------------------------!
    ! Modified to add initial support for local part of hybrid functionals      !
    ! by James Womack in 2019.                                                  !
    !---------------------------------------------------------------------------!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root
    use constants, only: stdout, DP
    use nlxc, only: nlxc_vdw_energy
    use rundat, only: pub_libxc_x_func_id, pub_libxc_c_func_id,pub_num_spins,&
       pub_lr_tddft_triplet
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_use_var
    use xc_funcs

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid    ! Grid definition
    type(CELL_INFO), intent(in) :: cell    ! Cell definition
    real(kind=DP), intent(inout) :: density_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(out) :: fxc_pot_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,2)
    ! @optimize ------------------------^
    ! fxc_pot_fine always has final dimension 2, even though
    ! fxc_pot_fine(:,:,:,:,2) seems to be only used when pub_lr_tddft_triplet
    ! is true. It is otherwise a waste of memory.
    real(kind=DP), intent(in) :: ground_state_dens(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    integer, intent(in) :: nhat_dim
    real(kind=DP), optional, intent(in) :: eff_nhat_den(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:nhat_dim)
    real(kind=DP), optional, intent(in) :: nhat_den_grad_gs(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:nhat_dim)

    logical, optional, intent(in) :: parallel

    ! la: new local variables used for VDWDF
    integer :: ierr
    real(kind=DP) :: xc_energy,exchange,correlation,&
            nl_energy,xc_potint
    real(kind=DP), allocatable :: nlxc_pot_fine(:,:,:,:)
    real(kind=DP), allocatable :: x_pot_fine(:,:,:,:)
    real(kind=DP), allocatable :: c_pot_fine(:,:,:,:)
    real(kind=DP), allocatable :: xc_pot_fine(:,:,:,:)
    real(kind=DP), allocatable :: triplet_dens(:,:,:,:)
    real(kind=DP), allocatable :: nhat_den_grad_triplet(:,:,:,:,:)
#ifdef WIN32
  real(kind=DP), allocatable, target :: density_grad(:,:,:,:,:)
  real(kind=DP), allocatable, target :: density_aux(:,:,:,:)
  complex(kind=DP), allocatable, target :: recip_work(:,:,:,:)
  real(kind=DP), allocatable, target :: conditioned_density_fine(:,:,:,:)
#else
  real(kind=DP), allocatable :: density_grad(:,:,:,:,:)
  real(kind=DP), allocatable :: density_aux(:,:,:,:)
  complex(kind=DP), allocatable :: recip_work(:,:,:,:)
  real(kind=DP), allocatable :: conditioned_density_fine(:,:,:,:)
#endif
  logical :: loc_parallel

    ! -------------------------------------------------------------------------

    call timer_clock('xc_fxc_finite_diff',1)

    loc_parallel = .true.
    if(present(parallel)) then
       loc_parallel = parallel
    end if



    call utils_use_var(pub_libxc_x_func_id)
    call utils_use_var(pub_libxc_c_func_id)

    if(pub_lr_tddft_triplet .and. pub_num_spins==1 .and.&
        present(eff_nhat_den) .and. present(nhat_den_grad_gs)) then
       allocate(nhat_den_grad_triplet(grid%ld1,grid%ld2,grid%max_slabs12,&
         pub_num_spins,nhat_dim), stat=ierr)
       call utils_alloc_check('xc_fxc_finite_diff','nhat_den_grad_triplet',ierr)
    endif
    if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
       allocate(xc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12,&
       2), stat=ierr)
    else
       allocate(xc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12,&
       pub_num_spins), stat=ierr)
    endif
    call utils_alloc_check('xc_fxc_finite_diff','xc_pot_fine',ierr)
    allocate(x_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12, &
         pub_num_spins),stat=ierr)
    call utils_alloc_check('xc_fxc_finite_diff','x_pot_fine',ierr)
    allocate(c_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12, &
         pub_num_spins),stat=ierr)
    call utils_alloc_check('xc_fxc_finite_diff','c_pot_fine',ierr)
    ! treat triplets separately in spin degenerate case to deal
    ! correctly with non-vanishing terms due to derivatives of the
    ! fully spin polarised case.
    if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
       allocate(triplet_dens(grid%ld1,grid%ld2,grid%max_slabs12, &
        2), stat=ierr)
       call utils_alloc_check('xc_fxc_finite_diff','triplet_dens',ierr)
    endif
    if(pub_xc_gradient_corrected) then
       ! Allocate workspace for calculating gradients of the density
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          allocate(density_grad(grid%ld1,grid%ld2,grid%max_slabs12,3, &
               2),stat=ierr)
          call utils_alloc_check('xc_fxc_finite_diff', &
               'density_grad',ierr)
       else
          allocate(density_grad(grid%ld1,grid%ld2,grid%max_slabs12,3, &
               pub_num_spins),stat=ierr)
          call utils_alloc_check('xc_fxc_finite_diff', &
               'density_grad',ierr)
       endif
       ! Allocate workspace for the auxiliary density array
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          allocate(density_aux(grid%ld1,grid%ld2,grid%max_slabs12, &
               2),stat=ierr)
          call utils_alloc_check('xc_fxc_finite_diff', &
               'density_aux',ierr)
       else
          allocate(density_aux(grid%ld1,grid%ld2,grid%max_slabs12, &
               pub_num_spins),stat=ierr)
          call utils_alloc_check('xc_fxc_finite_diff', &
               'density_aux',ierr)
       endif
       ! Allocate reciprocal space workspace
       allocate(recip_work(grid%ld3,grid%ld2,grid%max_slabs23,3),stat=ierr)
       call utils_alloc_check('xc_fxc_finite_diff', &
            'recip_work',ierr)
    endif

    ! Check if this is a PAW GGA triplet calculation. Currently Triplet
    ! calculations are not enabled in conjunction with PAW and GGA functionals
    if (pub_lr_tddft_triplet .and. (pub_num_spins==1) .and. &
         pub_xc_gradient_corrected .and. present(eff_nhat_den)) &
       call utils_abort('Error in xc_mod: &
            &TDDFT triplet calculations currently not implemented with PAW &
            &and GGA functionals.')

    ! construct correct finite difference density in case we have a triplet
    ! state and pub num spins=1
    ! triplet_dens1=rho_0_up+epsilon_rho1_up
    ! triplet_dens2=rho_0_up
    if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
       triplet_dens(:,:,:,1)=0.5_DP*density_fine(:,:,:,1)
       triplet_dens(:,:,:,2)=ground_state_dens(:,:,:,1)
    endif

    ! condition density to filter out low values
    allocate(conditioned_density_fine(grid%ld1,grid%ld2,grid%max_slabs12,&
         2),stat=ierr)
    call utils_alloc_check('xc_fxc_finite_diff','conditioned_density_fine',&
         ierr)
    if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
       call xc_condition_density(conditioned_density_fine(:,:,:,1),&
          triplet_dens(:,:,:,1),grid)
       call xc_condition_density(conditioned_density_fine(:,:,:,2),&
          triplet_dens(:,:,:,2),grid)
    else
       call xc_condition_density(conditioned_density_fine,density_fine,grid)
    endif

    ! Calculate gradients of the density if needed
    ! ndmh: supplying compensation density and its gradient if necessary
    if (pub_xc_gradient_corrected) then
       if (present(eff_nhat_den) .and. present(nhat_den_grad_gs)) then
          if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
             nhat_den_grad_triplet=0.5_DP*eff_nhat_den
             call xc_gradients(conditioned_density_fine(:,:,:,1),&
              density_grad(:,:,:,:,1),&
              recip_work,grid,nhat_dim,nhat_den_grad_triplet(:,:,:,:,:))
             nhat_den_grad_triplet=nhat_den_grad_gs
             call xc_gradients(conditioned_density_fine(:,:,:,2),&
              density_grad(:,:,:,:,2),&
              recip_work,grid,nhat_dim,nhat_den_grad_triplet(:,:,:,:,:))
          else
             call xc_gradients(conditioned_density_fine,density_grad,recip_work,grid, &
                  nhat_dim,eff_nhat_den)
          endif
       else
          if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
             call xc_gradients(conditioned_density_fine(:,:,:,1),&
              density_grad(:,:,:,:,1),recip_work,grid,0)
             call xc_gradients(conditioned_density_fine(:,:,:,2),&
              density_grad(:,:,:,:,2),&
              recip_work,grid,0)
          else
             call xc_gradients(conditioned_density_fine,density_grad,recip_work,grid,0)
          endif
       end if
    end if

    ! Zero potential
    fxc_pot_fine = 0.0_DP
    x_pot_fine = 0.0_DP
    c_pot_fine = 0.0_DP

    ! Call appropriate routine according to functional type
    select case (functional)
    case ('CAPZ')
       !----------------------------------------------------------------------!
       ! Local density approximation - from Ceperley-Alder Monte Carlo data,  !
       !  and Gell-Mann-Brueckner expansion, parameterised by Perdew & Zunger !
       !----------------------------------------------------------------------!
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          call xc_lda_triplet(conditioned_density_fine,xc_pot_fine,grid,xc_capz_c_point_sp)
          fxc_pot_fine(:,:,:,1,1)=xc_pot_fine(:,:,:,1)
          fxc_pot_fine(:,:,:,1,2)=xc_pot_fine(:,:,:,2)
       else
          call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
               xc_capz_c_point,xc_capz_c_point_sp,xc_potint,exchange,correlation,&
               parallel=loc_parallel)
          fxc_pot_fine(:,:,:,:,1)=xc_pot_fine
       endif

    case ('VWN')
       !----------------------------------------------------------------------!
       ! Local density approximation - Vosko, Wilk & Nusair                   !
       !----------------------------------------------------------------------!
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          call xc_lda_triplet(conditioned_density_fine,xc_pot_fine,grid,&
            xc_vwn_c_point_sp)
          fxc_pot_fine(:,:,:,1,1)=xc_pot_fine(:,:,:,1)
          fxc_pot_fine(:,:,:,1,2)=xc_pot_fine(:,:,:,2)
       else
          call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid,&
               xc_vwn_c_point,xc_vwn_c_point_sp,xc_potint,exchange,correlation,&
               parallel=loc_parallel)
          fxc_pot_fine(:,:,:,:,1)=xc_pot_fine(:,:,:,:)
       endif

    case ('PW92')
       !----------------------------------------------------------------------!
       ! Local density approximation - Perdew & Wang 1992 functional          !
       !----------------------------------------------------------------------!
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          call xc_lda_triplet(conditioned_density_fine,xc_pot_fine,grid,&
             xc_pw92_c_point_sp)
          fxc_pot_fine(:,:,:,1,1)=xc_pot_fine(:,:,:,1)
          fxc_pot_fine(:,:,:,1,2)=xc_pot_fine(:,:,:,2)
       else
          call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid,&
               xc_pw92_c_point, xc_pw92_c_point_sp,xc_potint,exchange,correlation,&
               parallel=loc_parallel)
          fxc_pot_fine(:,:,:,:,1)=xc_pot_fine(:,:,:,:)
       endif

    case ('BLYP')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: Becke88 + LYP            !
       !----------------------------------------------------------------------!
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          call xc_gc_triplet(conditioned_density_fine,density_grad,density_aux,&
               recip_work,grid,xc_pot_fine,xc_b88_x_point_sp,xc_lyp_c_point_sp)
          fxc_pot_fine(:,:,:,1,1)=xc_pot_fine(:,:,:,1)
          fxc_pot_fine(:,:,:,1,2)=xc_pot_fine(:,:,:,2)
       else
          call xc_gc(conditioned_density_fine,density_grad,density_aux,&
               recip_work,grid, xc_energy, xc_pot_fine, xc_b88_x_point,&
               xc_b88_x_point_sp,xc_lyp_c_point,xc_lyp_c_point_sp,xc_potint,&
               exchange,correlation,parallel=loc_parallel)

          fxc_pot_fine(:,:,:,:,1)=xc_pot_fine(:,:,:,:)
       endif

    case ('PBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1996 version of          !
       ! Perdew, Burke & Ernzerhof                                            !
       !----------------------------------------------------------------------!
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          call xc_gc_triplet(conditioned_density_fine,density_grad,density_aux,&
               recip_work,grid,xc_pot_fine,xc_pbe_x_point_sp,xc_pbe_c_point_sp)
          fxc_pot_fine(:,:,:,1,1)=xc_pot_fine(:,:,:,1)
          fxc_pot_fine(:,:,:,1,2)=xc_pot_fine(:,:,:,2)
       else
          call xc_gc(conditioned_density_fine,density_grad,density_aux,&
               recip_work,grid, xc_energy,xc_pot_fine,xc_pbe_x_point,&
               xc_pbe_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp,xc_potint,&
               exchange,correlation,parallel=loc_parallel)
          fxc_pot_fine(:,:,:,:,1)=xc_pot_fine(:,:,:,:)
       endif
    case ('PBEX')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1996 version of          !
       ! Perdew, Burke & Ernzerhof                                            !
       !----------------------------------------------------------------------!
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          call xc_gc_triplet(conditioned_density_fine,density_grad,density_aux,&
               recip_work,grid,xc_pot_fine,xc_pbe_x_point_sp,xc_none_c_point_sp)
          fxc_pot_fine(:,:,:,1,1)=xc_pot_fine(:,:,:,1)
          fxc_pot_fine(:,:,:,1,2)=xc_pot_fine(:,:,:,2)
       else
          call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,grid, &
               xc_energy,xc_pot_fine,xc_pbe_x_point,&
               xc_pbe_x_point_sp,xc_none_c_point,xc_none_c_point_sp,xc_potint,&
               exchange,correlation,parallel=loc_parallel)
          fxc_pot_fine(:,:,:,:,1)=xc_pot_fine(:,:,:,:)
       endif

    case ('REVPBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1998 revised version of  !
       ! PBE due to Zhang & Yang                                              !
       !----------------------------------------------------------------------!
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          call xc_gc_triplet(conditioned_density_fine,density_grad,density_aux,&
               recip_work,grid,xc_pot_fine,xc_revpbe_x_point_sp,xc_pbe_c_point_sp)
          fxc_pot_fine(:,:,:,1,1)=xc_pot_fine(:,:,:,1)
          fxc_pot_fine(:,:,:,1,2)=xc_pot_fine(:,:,:,2)
       else
          call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,grid, &
               xc_energy,xc_pot_fine,xc_revpbe_x_point,&
               xc_revpbe_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp,xc_potint,&
               exchange,correlation,parallel=loc_parallel)
          fxc_pot_fine(:,:,:,:,1)=xc_pot_fine(:,:,:,:)
       endif

    case ('RPBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1999 revised version of  !
       ! PBE due to Hammer, Hansen & Norskov                                  !
       !----------------------------------------------------------------------!
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          call xc_gc_triplet(conditioned_density_fine,density_grad,density_aux,&
               recip_work,grid,xc_pot_fine,xc_rpbe_x_point_sp,xc_pbe_c_point_sp)
          fxc_pot_fine(:,:,:,1,1)=xc_pot_fine(:,:,:,1)
          fxc_pot_fine(:,:,:,1,2)=xc_pot_fine(:,:,:,2)
       else
          call xc_gc(conditioned_density_fine,density_grad,density_aux,&
               recip_work,grid, xc_energy,xc_pot_fine,xc_rpbe_x_point,&
               xc_rpbe_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp,xc_potint,&
               exchange,correlation,parallel=loc_parallel)
          fxc_pot_fine(:,:,:,:,1)=xc_pot_fine(:,:,:,:)
       endif

    case ('PBESOL')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 2008 revised version of  !
       ! PBE for solids due to Perdew et al.                                  !
       !----------------------------------------------------------------------!
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          call xc_gc_triplet(conditioned_density_fine,density_grad,density_aux,&
               recip_work,grid,xc_pot_fine,xc_pbesol_x_point_sp,xc_pbesol_c_point_sp)
          fxc_pot_fine(:,:,:,1,1)=xc_pot_fine(:,:,:,1)
          fxc_pot_fine(:,:,:,1,2)=xc_pot_fine(:,:,:,2)
       else
          call xc_gc(conditioned_density_fine,density_grad,density_aux,&
               recip_work,grid, xc_energy,xc_pot_fine,xc_pbesol_x_point,&
               xc_pbesol_x_point_sp,xc_pbesol_c_point,xc_pbesol_c_point_sp,xc_potint,&
               exchange,correlation,parallel=loc_parallel)
          fxc_pot_fine(:,:,:,:,1)=xc_pot_fine(:,:,:,:)
       endif

    case ('AM05')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 2008 revised version of  !
       ! PBE for solids due to Perdew et al.                                  !
       !----------------------------------------------------------------------!

       if (pub_lr_tddft_triplet) &
            call utils_abort('Error in xc_mod: &
            &TDDFT triplet states not currently supported in AM05')

       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy,xc_pot_fine,xc_am05_x_point,&
            xc_am05_x_point_sp,xc_am05_c_point,xc_am05_c_point_sp,xc_potint,&
            exchange,correlation,parallel=loc_parallel)
       fxc_pot_fine(:,:,:,:,1)=xc_pot_fine(:,:,:,:)

    case ('PW91')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1991 version of          !
       ! Perdew and Wang                                                      !
       !----------------------------------------------------------------------!
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          call xc_gc_triplet(triplet_dens,density_grad,density_aux,recip_work,&
               grid,xc_pot_fine,xc_pw91_x_point_sp,xc_pw91_c_point_sp)
          fxc_pot_fine(:,:,:,1,1)=xc_pot_fine(:,:,:,1)
          fxc_pot_fine(:,:,:,1,2)=xc_pot_fine(:,:,:,2)
       else
          call xc_gc(conditioned_density_fine,density_grad,density_aux,&
               recip_work,grid, xc_energy,xc_pot_fine,xc_pw91_x_point, &
               xc_pw91_x_point_sp,xc_pw91_c_point,xc_pw91_c_point_sp,xc_potint,&
            exchange,correlation,parallel=loc_parallel)
          fxc_pot_fine(:,:,:,:,1)=xc_pot_fine(:,:,:,:)
       endif
    case ('XLYP')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: X + LYP                  !
       !----------------------------------------------------------------------!
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          call xc_gc_triplet(conditioned_density_fine,density_grad,&
               density_aux,recip_work,grid,xc_pot_fine,xc_x_x_point_sp,&
               xc_lyp_c_point_sp)
          fxc_pot_fine(:,:,:,1,1)=xc_pot_fine(:,:,:,1)
          fxc_pot_fine(:,:,:,1,2)=xc_pot_fine(:,:,:,2)
       else
          call xc_gc(conditioned_density_fine,density_grad,density_aux,&
               recip_work,grid, xc_energy, xc_pot_fine, xc_x_x_point, &
               xc_x_x_point_sp, xc_lyp_c_point, xc_lyp_c_point_sp,xc_potint,&
               exchange,correlation,parallel=loc_parallel)
          fxc_pot_fine(:,:,:,:,1)=xc_pot_fine(:,:,:,:)
       endif

    case ('WC')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 2006 functional of Wu &  !
       ! Cohen                                                                !
       !----------------------------------------------------------------------!
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          call xc_gc_triplet(conditioned_density_fine,density_grad,density_aux,&
               recip_work,grid,xc_pot_fine,xc_wc_x_point_sp,xc_pbe_c_point_sp)
          fxc_pot_fine(:,:,:,1,1)=xc_pot_fine(:,:,:,1)
          fxc_pot_fine(:,:,:,1,2)=xc_pot_fine(:,:,:,2)
       else
          call xc_gc(conditioned_density_fine,density_grad,density_aux,&
               recip_work,grid,xc_energy, xc_pot_fine, xc_wc_x_point, &
               xc_wc_x_point_sp, xc_pbe_c_point, xc_pbe_c_point_sp,xc_potint,&
               exchange,correlation,parallel=loc_parallel)
          fxc_pot_fine(:,:,:,:,1)=xc_pot_fine(:,:,:,:)
       endif
    ! --------------------------------------------------------------------------
    ! JCW: Hybrid functionals
    !      The contribution from the semi-local component of these
    !      functionals is evaluated here. The non-local Hartree-Fock exchange
    !      contribution to fxc is evaluated elsewhere.
    !
    ! JCW: @optimize
    !      While the existing system of selecting XC functionals via a single
    !      string persists, this select case statement must be manually updated
    !      when new hybrid functionals are added to ONETEP. Similar select case
    !      statements must be updated in xc_energy_potential and xc_radial.
    !      Ideally this will be replaced by a less-cumbersome system in future,
    !      e.g. utilizing derived types to represent and store data associated
    !      with XC functionals.
    !
    ! JCW: As of 27/03/19 the list of
    !      supported hybrid functionals is (extracted from xc_hfxinit
    !      routine): B1LYP, B1PW91, PBE0, B3LYP, B3PW91, X3LYP and HF (i.e.
    !      100% Hartree-Fock exchange).
    ! --------------------------------------------------------------------------
    case ('B1LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional B1LYP:                                             !
       ! 1/4 HF_X + 3/4 B88_X - LDA_X) + LYP_C                                !
       !----------------------------------------------------------------------!
       call utils_abort('Error in xc_fxc_finite_diff: &
            &Functional '//trim(functional)//' not implemented for LR-TDDFT yet.')
       !call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
       !     grid, xc_energy, xc_pot_fine, xc_b1_x_point, &
       !     xc_b1_x_point_sp, xc_lyp_c_point, xc_lyp_c_point_sp)
    case ('B1PW91')
       !----------------------------------------------------------------------!
       ! Hybrid functional B1PW91:                                             !
       ! 1/4 HF_X + 3/4 B88_X + PW91_C                                        !
       !----------------------------------------------------------------------!
       call utils_abort('Error in xc_fxc_finite_diff: &
            &Functional '//trim(functional)//' not implemented for LR-TDDFT yet.')
       !call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
       !     grid, xc_energy, xc_pot_fine, xc_b1_x_point, &
       !     xc_b1_x_point_sp, xc_pw91_c_point, xc_pw91_c_point_sp)
    case ('B3LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional B3LYP:                                             !
       ! LDA + a(HF_X - LDA_X) + b (B88_X - LDA_X) + c(LYP_C - VWN_C)         !
       !----------------------------------------------------------------------!
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          call utils_abort('Error in xc_fxc_finite_diff: &
               &Functional '//trim(functional)//' not implemented for triplet &
               &excitations in LR-TDDFT yet.')
          !call xc_gc_triplet(conditioned_density_fine,density_grad,density_aux,&
          !     recip_work,grid,xc_pot_fine,xc_pbe_x_point_sp,xc_pbe_c_point_sp)
          !fxc_pot_fine(:,:,:,1,1)=xc_pot_fine(:,:,:,1)
          !fxc_pot_fine(:,:,:,1,2)=xc_pot_fine(:,:,:,2)
       else
          call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
               grid, xc_energy, xc_pot_fine, xc_b3_x_point, &
               xc_b3_x_point_sp, xc_b3lyp_c_point, xc_b3lyp_c_point_sp,&
               parallel=loc_parallel)
          ! JCW: Optional xc_potint, exchange and correlation arguments for
          !      xc_gc are unneeded (results not used in xc_fxc_finite_diff)
          fxc_pot_fine(:,:,:,:,1)=xc_pot_fine(:,:,:,:)
       endif
    case ('B3PW91')
       !----------------------------------------------------------------------!
       ! Hybrid functional B3PW91:                                            !
       ! LDA + a(HF_X - LDA_X) + b (B88_X - LDA_X) + c(PW91_C - VWN_C)        !
       !----------------------------------------------------------------------!
       call utils_abort('Error in xc_fxc_finite_diff: &
            &Functional '//trim(functional)//' not implemented for LR-TDDFT yet.')
       !call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
       !     grid, xc_energy, xc_pot_fine, xc_b3_x_point, &
       !     xc_b3_x_point_sp, xc_b3pw91_c_point, xc_b3pw91_c_point_sp)
    case ('PBE0')
       !----------------------------------------------------------------------!
       ! Hybrid functional PBE0: PBE + 1/4 (HF_X - PBE_X)                     !
       !----------------------------------------------------------------------!
       if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
          call utils_abort('Error in xc_fxc_finite_diff: &
               &Functional '//trim(functional)//' not implemented for triplet &
               &excitations in LR-TDDFT yet.')
          !call xc_gc_triplet(conditioned_density_fine,density_grad,density_aux,&
          !     recip_work,grid,xc_pot_fine,xc_pbe_x_point_sp,xc_pbe_c_point_sp)
          !fxc_pot_fine(:,:,:,1,1)=xc_pot_fine(:,:,:,1)
          !fxc_pot_fine(:,:,:,1,2)=xc_pot_fine(:,:,:,2)
       else
          call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, xc_energy, xc_pot_fine, xc_pbe0_x_point,&
            xc_pbe0_x_point_sp, xc_pbe_c_point, xc_pbe_c_point_sp)
          ! JCW: Optional xc_potint, exchange and correlation arguments for
          !      xc_gc are unneeded (results not used in xc_fxc_finite_diff)
          fxc_pot_fine(:,:,:,:,1)=xc_pot_fine(:,:,:,:)
       end if
    case ('X3LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional X3LYP:                                             !
       ! LDA + a(HF_X - LDA_X) + b ((d B88_X + e PW91_X - LDA_X)              !
       !     + c(LYP_C - VWN_C)                                               !
       !----------------------------------------------------------------------!
       call utils_abort('Error in xc_fxc_finite_diff: &
            &Functional '//trim(functional)//' not implemented for LR-TDDFT yet.')
       !call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
       !     grid, xc_energy, xc_pot_fine, xc_x3_x_point, &
       !     xc_x3_x_point_sp, xc_x3lyp_c_point, xc_x3lyp_c_point_sp)
    case ('HF')
       !----------------------------------------------------------------------!
       ! 100% Hartree-Fock exchange (Hartree-Fock approximation).             !
       !----------------------------------------------------------------------!
       call utils_abort('Error in xc_fxc_finite_diff: &
            &Functional '//trim(functional)//' not implemented for LR-TDDFT yet.')
       xc_energy = 0.0_DP
       xc_pot_fine = 0.0_DP

    ! --------------------------------------------------------------------------
    ! JCW: End of hybrid functionals
    ! --------------------------------------------------------------------------

    case('VDWDF')
       !----------------------------------------------------------------------!
       ! vdW-DF functional                                                    !
       !----------------------------------------------------------------------!

       if (pub_lr_tddft_triplet) &
            call utils_abort('Error in xc_mod: &
            &TDDFT triplet states not currently supported in VDWDF')

       allocate(nlxc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

       ! Calculate the LDA correlation potential and GGA exchange potential
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
            xc_pw92_c_point, xc_pw92_c_point_sp,correlation=correlation,&
            c_pot_fine=c_pot_fine,parallel=loc_parallel)
       call xc_gc(conditioned_density_fine,density_grad,density_aux,&
            recip_work,grid, exchange,x_pot_fine,xc_revpbe_x_point,&
            xc_revpbe_x_point_sp,xc_none_c_point,xc_none_c_point_sp,&
            parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(density_fine,density_grad,recip_work,grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,density_fine,density_grad, &
            grid,cell,'VDWDF1')
       fxc_pot_fine(:,:,:,1,1) = x_pot_fine(:,:,:,1) + c_pot_fine(:,:,:,1) +&
            nlxc_pot_fine(:,:,:,1)

       call timer_clock('nlxc_vdw_energy',2)

       deallocate(nlxc_pot_fine,stat=ierr)
       call utils_dealloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

    case('VDWDF2')
       !----------------------------------------------------------------------!
       ! vdW-DF2 functional                                                   !
       !----------------------------------------------------------------------!
       if (pub_lr_tddft_triplet) &
            call utils_abort('Error in xc_mod: &
            &TDDFT triplet states not currently supported in VDWDF2')

       allocate(nlxc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

       ! JA: Changed xc_revpbe_x_point_sp -> xc_rpw86_x_point_sp.
       !     This matches the spec & the non spin pol version.
       ! Calculate the LDA correlation potential and GGA exchange potential
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,grid, &
            exchange,x_pot_fine,xc_rpw86_x_point,&
            xc_rpw86_x_point_sp,xc_none_c_point,xc_none_c_point_sp,&
            parallel=loc_parallel)
       call xc_lda(density_fine,xc_energy, xc_pot_fine, grid, xc_pw92_c_point, &
            xc_pw92_c_point_sp,correlation=correlation,c_pot_fine=c_pot_fine,&
            parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work,grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,density_fine,density_grad, &
            grid,cell,'VDWDF2')

       fxc_pot_fine(:,:,:,1,1) = x_pot_fine(:,:,:,1) + c_pot_fine(:,:,:,1)+&
            nlxc_pot_fine(:,:,:,1)
       call timer_clock('nlxc_vdw_energy',2)

       deallocate(nlxc_pot_fine,stat=ierr)
       call utils_dealloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

   case('VDWC09')
       !----------------------------------------------------------------------!
       ! vdW-DF-C09x functional                                               !
       !----------------------------------------------------------------------!
       if (pub_lr_tddft_triplet) &
            call utils_abort('Error in xc_mod: &
            &TDDFT triplet states not currently supported in VDWC09')

       allocate(nlxc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

       ! JA: FIXME: This has the wrong spin polarized exchange functional
       ! Calculate the LDA correlation potential and GGA exchange potential
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,grid, &
            exchange,x_pot_fine,xc_c09_x_point,&
            xc_c09_x_point_sp,xc_none_c_point,xc_none_c_point_sp,parallel=loc_parallel)
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, xc_pw92_c_point, &
            xc_pw92_c_point_sp,correlation=correlation,c_pot_fine=c_pot_fine,&
            parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work,grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,density_fine,density_grad, &
            grid,cell,'VDWDF1')
       fxc_pot_fine(:,:,:,1,1) = x_pot_fine(:,:,:,1) + c_pot_fine(:,:,:,1) +&
            nlxc_pot_fine(:,:,:,1)

       call timer_clock('nlxc_vdw_energy',2)

       deallocate(nlxc_pot_fine,stat=ierr)
       call utils_dealloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

   case('C09-2')
       !----------------------------------------------------------------------!
       ! vdW-DF-C09x functional with vdW-DF2 non-local correlation            !
       !----------------------------------------------------------------------!
       if (pub_lr_tddft_triplet) &
            call utils_abort('Error in xc_mod: &
            &TDDFT triplet states not currently supported in C09-2')

       allocate(nlxc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

       ! Calculate the LDA correlation potential and GGA exchange potential
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,grid, &
            exchange,x_pot_fine,xc_c09_x_point,&
            xc_c09_x_point_sp,xc_none_c_point,xc_none_c_point_sp,parallel=loc_parallel)
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
            xc_pw92_c_point,xc_pw92_c_point_sp,correlation=correlation,&
            c_pot_fine=c_pot_fine,parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work,grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,density_fine,density_grad, &
            grid,cell,'VDWDF2')
       fxc_pot_fine(:,:,:,1,1) = x_pot_fine(:,:,:,1) + c_pot_fine(:,:,:,1) + &
            nlxc_pot_fine(:,:,:,1)

       call timer_clock('nlxc_vdw_energy',2)
       deallocate(nlxc_pot_fine,stat=ierr)
       call utils_dealloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

    case('OPTPBE')
       !----------------------------------------------------------------------!
       ! vdW-DF-optPBE functional                                             !
       !----------------------------------------------------------------------!
       if (pub_lr_tddft_triplet) &
            call utils_abort('Error in xc_mod: &
            &TDDFT triplet states not currently supported in OPTPBE')

       allocate(nlxc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

       ! Calculate the LDA correlation potential and GGA exchange potential
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,grid, &
            exchange,x_pot_fine,xc_vdwoptpbe_x_point,&
            xc_vdwoptpbe_x_point_sp,xc_none_c_point,xc_none_c_point_sp,parallel=loc_parallel)
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, xc_pw92_c_point, &
            xc_pw92_c_point_sp,correlation=correlation,c_pot_fine=c_pot_fine,&
            parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work,grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,density_fine,density_grad, &
            grid,cell,'VDWDF1')
       fxc_pot_fine(:,:,:,1,1) = x_pot_fine(:,:,:,1) + c_pot_fine(:,:,:,1)+ &
            nlxc_pot_fine(:,:,:,1)

       call timer_clock('nlxc_vdw_energy',2)
       deallocate(nlxc_pot_fine,stat=ierr)
       call utils_dealloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

    case('OPTB88')
       !----------------------------------------------------------------------!
       ! vdW-DF-optB88 functional                                             !
       !----------------------------------------------------------------------!
       if (pub_lr_tddft_triplet) &
            call utils_abort('Error in xc_mod: &
            &TDDFT triplet states not currently supported in OPTB88')

       if (pub_num_spins>1) &
            call utils_abort('Error in xc_mod: &
            &OPTB88 does not yet support spin-polarization.') ! JA: needs x_point

       allocate(nlxc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

       !Calculate the LDA correlation potential and GGA exchange potential
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,grid, &
            exchange,x_pot_fine,xc_vdwoptb88_x_point,&
            xc_vdwoptb88_x_point_sp,xc_none_c_point,xc_none_c_point_sp,parallel=loc_parallel)
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
            xc_pw92_c_point,xc_pw92_c_point_sp,correlation=correlation,&
            c_pot_fine=c_pot_fine,parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work,grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,density_fine,density_grad, &
            grid,cell,'VDWDF1')
       fxc_pot_fine(:,:,:,1,1) = x_pot_fine(:,:,:,1) + c_pot_fine(:,:,:,1) +&
            nlxc_pot_fine(:,:,:,1)

       call timer_clock('nlxc_vdw_energy',2)
       deallocate(nlxc_pot_fine,stat=ierr)
       call utils_dealloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

   case('PBEK1')
       !----------------------------------------------------------------------!
       ! vdW-DF-optPBEK=1 functional                                          !
       !----------------------------------------------------------------------!
       if (pub_lr_tddft_triplet) &
            call utils_abort('Error in xc_mod: &
            &TDDFT triplet states not currently supported in PBEK1')

       allocate(nlxc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,grid, &
            exchange,x_pot_fine,xc_vdwpbek1_x_point,&
            xc_vdwpbek1_x_point_sp,xc_none_c_point,xc_none_c_point_sp,parallel=loc_parallel)
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
            xc_pw92_c_point, xc_pw92_c_point_sp,correlation=correlation,&
            c_pot_fine=c_pot_fine,parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work,grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,density_fine,density_grad, &
            grid,cell,'VDWDF1')
       fxc_pot_fine(:,:,:,1,1) = x_pot_fine(:,:,:,1) + c_pot_fine(:,:,:,1) + &
            nlxc_pot_fine(:,:,:,1)

       call timer_clock('nlxc_vdw_energy',2)
       deallocate(nlxc_pot_fine,stat=ierr)
       call utils_dealloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

   case('VV10')
       !----------------------------------------------------------------------!
       ! rVV10 functional, as implemented by Sabatini et al.                  !
       !----------------------------------------------------------------------!
       if (pub_lr_tddft_triplet) &
            call utils_abort('Error in xc_mod: &
            &TDDFT triplet states not currently supported in VV10')

       allocate(nlxc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy,xc_pot_fine,xc_rpw86_x_point,&
            xc_rpw86_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp,xc_potint,&
            exchange,correlation,parallel=loc_parallel)
       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work,grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,density_fine, &
            density_grad, grid,cell,'VV10')

       fxc_pot_fine(:,:,:,1,1) = xc_pot_fine(:,:,:,1) +&
            nlxc_pot_fine(:,:,:,1)

       call timer_clock('nlxc_vdw_energy',2)
       deallocate(nlxc_pot_fine,stat=ierr)
       call utils_dealloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)
   case('VV10BC')
       !----------------------------------------------------------------------!
       ! rVV10 functional, as implemented by Sabatini et al., but with b and  !
       ! C parameters from Vydrov and Van Voorhis original paper, i.e. b=5.9, !
       ! C=0.0093.                                                            !
       !----------------------------------------------------------------------!
       if (pub_lr_tddft_triplet) &
            call utils_abort('Error in xc_mod: &
            &TDDFT triplet states not currently supported in VV10')

       allocate(nlxc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy,xc_pot_fine,xc_rpw86_x_point,&
            xc_rpw86_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp,xc_potint,&
            exchange,correlation,parallel=loc_parallel)
       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work,grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,density_fine, &
            density_grad, grid,cell,'VV10BC')

       fxc_pot_fine(:,:,:,1,1) = xc_pot_fine(:,:,:,1) +&
            nlxc_pot_fine(:,:,:,1)

       call timer_clock('nlxc_vdw_energy',2)
       deallocate(nlxc_pot_fine,stat=ierr)
       call utils_dealloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

   case('AVV10S')
       !----------------------------------------------------------------------!
       ! AM05-VV10sol functional                                              !
       !----------------------------------------------------------------------!
       if (pub_lr_tddft_triplet) &
            call utils_abort('Error in xc_mod: &
            &TDDFT triplet states not currently supported in AVV10S')

       allocate(nlxc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy,xc_pot_fine,xc_am05_x_point,&
            xc_am05_x_point_sp,xc_am05_c_point,xc_am05_c_point_sp, &
            xc_potint, exchange,correlation,parallel=loc_parallel)
       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work,grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,density_fine, &
            density_grad, grid,cell,'VV10-S')

       fxc_pot_fine(:,:,:,1,1) = xc_pot_fine(:,:,:,1)  +&
            nlxc_pot_fine(:,:,:,1)

       call timer_clock('nlxc_vdw_energy',2)
       deallocate(nlxc_pot_fine,stat=ierr)
       call utils_dealloc_check('xc_fxc_energy_finite_diff','nlxc_pot_fine',ierr)

    case default
       !----------------------------------------------------------------------!
       ! Unrecognised functional type                                         !
       !----------------------------------------------------------------------!
       if (pub_on_root) write(stdout,'(2a)') 'WARNING: unknown exchange-&
            &correlation functional in xc_fxc_finite_diff: ', trim(functional)
       xc_energy = 0.0_DP
       xc_pot_fine = 0.0_DP
       fxc_pot_fine =0.0_DP
    end select

    ! deallocate data
    deallocate(conditioned_density_fine,stat=ierr)
    call utils_dealloc_check('xc_fxc_finite_diff','conditioned_density_fine',&
         ierr)
    deallocate(xc_pot_fine, stat=ierr)
    call utils_dealloc_check('xc_fxc_finite_diff','xc_pot_fine',ierr)
    deallocate(x_pot_fine,stat=ierr)
    call utils_dealloc_check('xc_fxc_finite_diff','x_pot_fine',ierr)
    deallocate(c_pot_fine,stat=ierr)
    call utils_dealloc_check('xc_fxc_finite_diff','c_pot_fine',ierr)
    if(pub_lr_tddft_triplet .and. pub_num_spins==1) then
       deallocate(triplet_dens,stat=ierr)
       call utils_dealloc_check('xc_fxc_finite_diff','triplet_dens',ierr)
    endif
    if(pub_xc_gradient_corrected) then
       deallocate(density_grad,stat=ierr)
       call utils_dealloc_check('xc_fxc_finite_diff', &
            'density_grad',ierr)
       deallocate(density_aux,stat=ierr)
       call utils_dealloc_check('xc_fxc_finite_diff', &
            'density_aux',ierr)
       deallocate(recip_work,stat=ierr)
       call utils_dealloc_check('xc_fxc_finite_diff', &
            'recip_work',ierr)
    endif
    if(pub_lr_tddft_triplet .and. pub_num_spins==1 .and.&
        present(eff_nhat_den) .and. present(nhat_den_grad_gs)) then
       deallocate(nhat_den_grad_triplet, stat=ierr)
       call utils_dealloc_check('xc_fxc_finite_diff',&
          'nhat_den_grad_triplet',ierr)
    endif

    call timer_clock('xc_fxc_finite_diff',2)


  end subroutine xc_fxc_finite_diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine xc_energy_potential(density_fine,xc_energy,xc_pot_fine,grid,cell, &
       nhat_dim,nhat_den_grad,tddftxc,ke_density_fine,dfdtau_fine,initial,&
       parallel)

    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the exchange-correlation energy and the exchange-correlation potential,   !
    ! according to the functional defined by the input file.                    !
    !---------------------------------------------------------------------------!
    ! Arguments:                                                                !
    !    density_fine (in) : Input density on the fine grid                     !
    !    grid (in)         : GRID_INFO grid definition                          !
    !    xc_energy (out)   : Exchange-correlation energy                        !
    !    xc_pot_fine (out) : Exchange-correlation potential on grid             !
    !    nhat_dim (in)     : Number of components of compensation density (if   !
    !                        present - ignored if not present)                  !
    !    nhat_den_grad (in): Compensation density and Compensation density      !
    !                        gradient (if required).                            !
    !    tddftxc (in)      : Switch to override exchange-correlation functional !
    !                        with TDDFT exchange-correlation functional         !
    !    ke_density_fine (in)                                                   !
    !                      : Optional input KE density on the fine grid         !
    !    dfdtau_fine (out) : Optional output gradient of XC energy per unit     !
    !                        volume wrt kinetic energy density                  !
    !    initial (in)      : Optional, if .true. then use initial_functional    !
    !                        rather than functional                             !
    !---------------------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004                         !
    ! Addition of hybrid functionals and modification to use new subroutines by !
    ! Quintin Hill 11/03/2009.                                                  !
    ! Kinetic-energy-density-dependence added by James C. Womack 2015           !
    !===========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root
    use constants, only: stdout, DP, VERBOSE
    use nlxc, only: nlxc_vdw_energy
    use rundat, only: pub_output_detail, pub_libxc_x_func_id, &
         pub_libxc_c_func_id, pub_num_spins, pub_debug_on_root, &
         pub_xc_ke_density_required, pub_xc_lapl_density_required
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_assert, utils_use_var
    use xc_funcs

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid    ! Grid definition
    type(CELL_INFO), intent(in) :: cell    ! Cell definition
    real(kind=DP), intent(in)  :: density_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(out) :: xc_energy
    real(kind=DP), intent(out) :: xc_pot_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    integer, intent(in) :: nhat_dim
    real(kind=DP), optional, intent(in) :: nhat_den_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:nhat_dim)
    logical, optional, intent(in) :: tddftxc
    ! JCW: optional KE density on grid
    real(kind=DP), optional, intent(in) :: ke_density_fine(:,:,:,:)
    ! JCW: optional gradient of XC energy per unit volume wrt KE density
    real(kind=DP), optional, intent(out) :: dfdtau_fine(:,:,:,:)
    ! JCW: Expect the following dimensions (or otherwise a dummy array with
    ! JCW: size == 1).
    ! JCW:    ke_density_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins)
    ! JCW:    dfdtau_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins)
    logical, optional, intent(in) :: initial ! JCW: optional argument, if true
    ! JCW: then we use initial_functional for this call only
    logical, optional, intent(in) :: parallel

    ! Local variables
    integer :: ierr
    real(kind=DP), allocatable :: conditioned_density_fine(:,:,:,:)
    ! JCW: Conditioning may also be necessary for ke_density
    real(kind=DP), allocatable :: conditioned_ke_density_fine(:,:,:,:)

    ! la: new local variables used for VDWDF
    real(kind=DP) :: exchange
    real(kind=DP) :: correlation
    real(kind=DP) :: nl_energy
    real(kind=DP), allocatable :: nlxc_pot_fine(:,:,:,:)
    real(kind=DP), allocatable :: x_pot_fine(:,:,:,:)
    real(kind=DP), allocatable :: c_pot_fine(:,:,:,:)
#ifdef WIN32
    real(kind=DP), allocatable, target :: density_grad(:,:,:,:,:)
    real(kind=DP), allocatable, target :: density_aux(:,:,:,:)
    complex(kind=DP), allocatable, target :: recip_work(:,:,:,:)
#else
    real(kind=DP), allocatable :: density_grad(:,:,:,:,:)
    real(kind=DP), allocatable :: density_aux(:,:,:,:)
    complex(kind=DP), allocatable :: recip_work(:,:,:,:)
#endif
    logical :: loc_tddftxc ! ddor: A local copy of tddftxc
    character(len=80) :: loc_functional ! JCW: local copy of functional (or
    ! initial_functional, if initial==.true.)
    logical :: ke_density_required ! JCW: local copy of pub_xc_ke_density_required
    logical :: loc_parallel


    ! -------------------------------------------------------------------------

    loc_parallel = .true.
    if(present(parallel)) then
       loc_parallel = parallel
    end if


    ! Set local copy of functional and pub_xc_ke_density_required
    ! Defaults
    loc_functional = functional
    ke_density_required = pub_xc_ke_density_required
    if (present(initial)) then
       ! Override local copy of functional if initial==.true.
       if (initial) then
          ! JCW: Use an initial functional, different to the functional used
          ! JCW: elsewhere (only available for meta-GGAs).
          ! JCW: NB. loc_functional must be a GGA in this case
          ! JCW: (enforced in xc_init).
          call utils_assert(pub_xc_ke_density_required,"Error in &
               &xc_energy_potential: initial==.true. is only allowed for &
               &meta-GGA XC functionals.")
          loc_functional = initial_functional
          ! JCW: If initial==.true., we have an alternative GGA functional,
          ! JCW: so do not need to deal with ke_density or dfdtau_fine
          ke_density_required = .false.
          if (pub_on_root) write(stdout,*) "WARNING: xc_energy_potential using &
             &alternative initial XC functional "//trim(loc_functional)//"."
       end if
    end if

    ! JCW: Use locally set ke_density_required variable to verify arguments
    ! JCW: have been correctly passed to this routine (whether "dfdtau_fine" and
    ! JCW: "ke_density_fine" are required depends on the value of "initial").
    if (.not.ke_density_required) then
       ! JCW: Non-tau-dependent XC functional
       ! JCW: dfdtau_fine is either not present, or has a length of 1
       if ( present(dfdtau_fine) ) then
          call utils_assert( &
               size(dfdtau_fine) == 1,&
               "Error in xc_energy_potential: &
               &pub_xc_ke_density_required is false, but dfdtau_fine is present &
               &and has a size > 1.")
       end if
       if ( present(ke_density_fine) ) then
          call utils_assert( &
               size(ke_density_fine) == 1,&
               "Error in xc_energy_potential: &
               &pub_xc_ke_density_required is false, but ke_density_fine is present &
               &and has a size > 1.")
       end if
    else
       ! JCW: tau-dependent XC functional
       ! JCW: Check dfdtau_fine is present
       call utils_assert( present(dfdtau_fine),&
            "Error in xc_energy_potential: &
            &pub_xc_ke_density_required is true, but dfdtau_fine is not &
            &present.")
       ! JCW: Check ke_density_fine is present
       call utils_assert( present(ke_density_fine),&
            "Error in xc_energy_potential: &
            &pub_xc_ke_density_required is true, but ke_density_fine is not &
            &present.")
    end if

    ! JCW: XC functionals dependent on the Laplacian of the density are not
    ! JCW: supported. Abort.
    call utils_assert(.not.pub_xc_lapl_density_required,&
         "Error in xc_energy_potential: Dependence of XC &
         &functionals on the Laplacian of the density is not &
         &supported.")

    if (pub_debug_on_root) write(stdout,'(a)') "DEBUG: &
       &xc_energy_potential, loc_functional = "//trim(loc_functional)

    ! Start timer
    call timer_clock('xc_energy_potential',1)

    call utils_use_var(pub_libxc_x_func_id)
    call utils_use_var(pub_libxc_c_func_id)

    ! ddor: Set default for loc_tddftxc
    loc_tddftxc = .false.
    if (present(tddftxc)) loc_tddftxc = tddftxc

    ! Check input density
    if (pub_output_detail >= VERBOSE) call internal_check_density(loc_parallel)

    ! jd: Filter out extremely low density values (and negative densities),
    !     so that BLYP and friends do not become ill-conditioned.
    allocate(conditioned_density_fine(grid%ld1,grid%ld2,grid%max_slabs12,&
         pub_num_spins),stat=ierr)
    call utils_alloc_check('xc_energy_potential','conditioned_density_fine',&
         ierr)
    call xc_condition_density(conditioned_density_fine,density_fine,grid)

    ! JCW: Also filtering ke_density_fine, to avoid problems with meta-GGAs
    ! JCW: when very small / negative ke_density_fine values arise
    if (ke_density_required) then
      call utils_assert(present(ke_density_fine),&
           "Error in xc_energy_potential: ke_density_required is true, &
           &but ke_density_fine has not been passed as an argument.")
      allocate(conditioned_ke_density_fine(grid%ld1,grid%ld2, &
           grid%max_slabs12,pub_num_spins),stat=ierr)
      call utils_alloc_check('xc_energy_potential',&
           'conditioned_ke_density_fine',ierr)
      ! JCW: Note that xc_condition_ke_density only removes negative
      ! JCW: ke_density values, setting these to zeroes.
      call xc_condition_ke_density(conditioned_ke_density_fine,ke_density_fine,grid)
    else
      ! JCW: Allocate a single element dummy array
      allocate(conditioned_ke_density_fine(1,1,1,1),stat=ierr)
      call utils_alloc_check('xc_energy_potential',&
           'conditioned_ke_density_fine',ierr)
    end if

    ! Allocate workspace if required
#ifdef LIBXC
    if (pub_xc_gradient_corrected.or.use_libxc) call internal_allocate_workspace
#else
    if (pub_xc_gradient_corrected) call internal_allocate_workspace
#endif

    ! Calculate gradients of the density if needed
    ! ndmh: supplying compensation density and its gradient if necessary
    if (pub_xc_gradient_corrected) then
       if (present(nhat_den_grad)) then
          call xc_gradients(conditioned_density_fine,density_grad,recip_work, &
               grid,nhat_dim,nhat_den_grad)
       else
          call xc_gradients(conditioned_density_fine,density_grad,recip_work,&
               grid,0)
       end if
    end if

    ! Zero potential
    xc_pot_fine = 0.0_DP

    ! Call appropriate routine according to functional type
    select case (loc_functional)
    case ('CAPZ')
       !----------------------------------------------------------------------!
       ! Local density approximation - from Ceperley-Alder Monte Carlo data,  !
       !  and Gell-Mann-Brueckner expansion, parameterised by Perdew & Zunger !
       !----------------------------------------------------------------------!
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
            xc_capz_c_point, xc_capz_c_point_sp,parallel=loc_parallel)
    case ('VWN')
       !----------------------------------------------------------------------!
       ! Local density approximation - Vosko, Wilk & Nusair                   !
       !----------------------------------------------------------------------!
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
            xc_vwn_c_point, xc_vwn_c_point_sp,parallel=loc_parallel)
    case ('PW92')
       !----------------------------------------------------------------------!
       ! Local density approximation - Perdew & Wang 1992 functional          !
       !----------------------------------------------------------------------!
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
            xc_pw92_c_point, xc_pw92_c_point_sp,parallel=loc_parallel)
    case ('BLYP')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: Becke88 + LYP            !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, xc_energy, xc_pot_fine, xc_b88_x_point,xc_b88_x_point_sp, &
            xc_lyp_c_point, xc_lyp_c_point_sp,parallel=loc_parallel)
    case ('PBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1996 version of          !
       ! Perdew, Burke & Ernzerhof                                            !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, xc_energy,xc_pot_fine,xc_pbe_x_point,&
            xc_pbe_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp,&
            parallel=loc_parallel)
    case ('PBEX')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1996 version of          !
       ! Perdew, Burke & Ernzerhof                                            !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, xc_energy,xc_pot_fine,xc_pbe_x_point,&
            xc_pbe_x_point_sp,xc_none_c_point,xc_none_c_point_sp,&
            parallel=loc_parallel)
    case ('REVPBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1998 revised version of  !
       ! PBE due to Zhang & Yang                                              !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid,xc_energy,xc_pot_fine,xc_revpbe_x_point,&
            xc_revpbe_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp,&
            parallel=loc_parallel)
    case ('RPBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1999 revised version of  !
       ! PBE due to Hammer, Hansen & Norskov                                  !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid,xc_energy,xc_pot_fine,xc_rpbe_x_point,&
            xc_rpbe_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp,&
            parallel=loc_parallel)
    case ('PBESOL')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 2008 revised version of  !
       ! PBE for solids due to Perdew et al.                                  !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid,xc_energy,xc_pot_fine,xc_pbesol_x_point,&
            xc_pbesol_x_point_sp,xc_pbesol_c_point,xc_pbesol_c_point_sp,&
            parallel=loc_parallel)
    case ('AM05')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 2005 AM05                !
       !----------------------------------------------------------------------!

       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid,xc_energy,xc_pot_fine,xc_am05_x_point,&
            xc_am05_x_point_sp,xc_am05_c_point,xc_am05_c_point_sp,&
            parallel=loc_parallel)
    case ('PW91')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1991 version of          !
       ! Perdew and Wang                                                      !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid,xc_energy,xc_pot_fine,xc_pw91_x_point, &
            xc_pw91_x_point_sp,xc_pw91_c_point,xc_pw91_c_point_sp,&
            parallel=loc_parallel)
    case ('XLYP')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: X + LYP                  !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, xc_energy, xc_pot_fine, xc_x_x_point, &
            xc_x_x_point_sp, xc_lyp_c_point, xc_lyp_c_point_sp,&
            parallel=loc_parallel)
    case ('WC')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 2006 functional of Wu &  !
       ! Cohen                                                                !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, xc_energy, xc_pot_fine, xc_wc_x_point, &
            xc_wc_x_point_sp, xc_pbe_c_point, xc_pbe_c_point_sp,&
            parallel=loc_parallel)
    case ('B1LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional B1LYP:                                             !
       ! 1/4 HF_X + 3/4 B88_X - LDA_X) + LYP_C                                !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, xc_energy, xc_pot_fine, xc_b1_x_point, &
            xc_b1_x_point_sp, xc_lyp_c_point, xc_lyp_c_point_sp,&
            parallel=loc_parallel)
    case ('B1PW91')
       !----------------------------------------------------------------------!
       ! Hybrid functional B1PW91:                                             !
       ! 1/4 HF_X + 3/4 B88_X + PW91_C                                        !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, xc_energy, xc_pot_fine, xc_b1_x_point, &
            xc_b1_x_point_sp, xc_pw91_c_point, xc_pw91_c_point_sp,&
            parallel=loc_parallel)
    case ('B3LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional B3LYP:                                             !
       ! LDA + a(HF_X - LDA_X) + b (B88_X - LDA_X) + c(LYP_C - VWN_C)         !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, xc_energy, xc_pot_fine, xc_b3_x_point, &
            xc_b3_x_point_sp, xc_b3lyp_c_point, xc_b3lyp_c_point_sp,&
            parallel=loc_parallel)
    case ('B3PW91')
       !----------------------------------------------------------------------!
       ! Hybrid functional B3PW91:                                            !
       ! LDA + a(HF_X - LDA_X) + b (B88_X - LDA_X) + c(PW91_C - VWN_C)        !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, xc_energy, xc_pot_fine, xc_b3_x_point, &
            xc_b3_x_point_sp, xc_b3pw91_c_point, xc_b3pw91_c_point_sp,&
            parallel=loc_parallel)
    case ('PBE0')
       !----------------------------------------------------------------------!
       ! Hybrid functional PBE0: PBE + 1/4 (HF_X - PBE_X)                     !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, xc_energy, xc_pot_fine, xc_pbe0_x_point,&
            xc_pbe0_x_point_sp, xc_pbe_c_point, xc_pbe_c_point_sp,&
            parallel=loc_parallel)
    case ('X3LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional X3LYP:                                             !
       ! LDA + a(HF_X - LDA_X) + b ((d B88_X + e PW91_X - LDA_X)              !
       !     + c(LYP_C - VWN_C)                                               !
       !----------------------------------------------------------------------!
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, xc_energy, xc_pot_fine, xc_x3_x_point, &
            xc_x3_x_point_sp, xc_x3lyp_c_point, xc_x3lyp_c_point_sp,&
            parallel=loc_parallel)
    case ('NONE','HF')
       !----------------------------------------------------------------------!
       ! No exchange-correlation functional (Hartree approximation) or only   !
       ! Hartree-Fock exchange (Hartree-Fock approximation).                  !
       !----------------------------------------------------------------------!
       xc_energy = 0.0_DP
       xc_pot_fine = 0.0_DP
    case ('LIBXC')
       !----------------------------------------------------------------------!
       ! Functional from LIBXC library - type specified by pub_libxc_func_id  !
       !----------------------------------------------------------------------!
#ifdef LIBXC
       if (pub_xc_gradient_corrected) then
          ! JCW: Function is a GGA, or meta-GGA
          call xc_libxc(pub_libxc_x_func_id,pub_libxc_c_func_id, &
               libxc_gga, libxc_mgga, &
               conditioned_density_fine,recip_work,&
               xc_energy, xc_pot_fine,grid, &
               density_grad = density_grad, &
               density_aux = density_aux, &
               ke_density_fine = conditioned_ke_density_fine, &
               dfdtau_fine = dfdtau_fine,parallel=loc_parallel)
       else
          ! JCW: LDA, no gradient or higher-order terms
          call xc_libxc(pub_libxc_x_func_id,pub_libxc_c_func_id, &
               libxc_gga, libxc_mgga, &
               conditioned_density_fine,recip_work,&
               xc_energy, xc_pot_fine,grid,parallel=loc_parallel)
       end if

       if (libxc_nlc_vv10) then
          ! JCW: NLC (VV10 required)

          ! JCW: Allocate workspace for NLC potential
          call internal_allocate_vdwdf

          ! Calculate nonlocal energy and potential
          call timer_clock('nlxc_vdw_energy',1)
          ! JCW: Note: In the future, we may want to stop reusing density_grad
          ! JCW: in xc_libxc, so it is not necessary to regenerate the density
          ! JCW: gradient for nlxc_vdw_energy.
          call xc_gradients(conditioned_density_fine,density_grad,recip_work, &
               grid,0)
          call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,conditioned_density_fine, &
               density_grad, grid,cell,'VV10-LIBXC',&
               vv10_bval=libxc_vv10_bval,vv10_Cval=libxc_vv10_Cval)
          xc_energy = xc_energy + nl_energy
          xc_pot_fine = xc_pot_fine + nlxc_pot_fine

          call timer_clock('nlxc_vdw_energy',2)

          ! JCW: Deallocate workspace for NLC potential
          call internal_deallocate_vdwdf
       end if
#else
       call utils_abort('Error in xc_energy_potential: Functional "LIBXC" &
            &specified, but LIBXC not present in build')
#endif
    ! la: The vdW-DF functional
    case('VDWDF')
       !----------------------------------------------------------------------!
       ! vdW-DF functional                                                    !
       !----------------------------------------------------------------------!

       call internal_allocate_vdwdf

       ! Calculate the LDA correlation potential and GGA exchange potential
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
            xc_pw92_c_point, xc_pw92_c_point_sp,correlation=correlation, &
            c_pot_fine=c_pot_fine,parallel=loc_parallel)
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, exchange,x_pot_fine,xc_revpbe_x_point,&
            xc_revpbe_x_point_sp,xc_none_c_point,xc_none_c_point_sp,&
            parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work,&
            grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,conditioned_density_fine,&
            density_grad, grid,cell,'VDWDF1')
       xc_energy = exchange  + correlation + nl_energy
       xc_pot_fine = x_pot_fine + c_pot_fine + nlxc_pot_fine

       call timer_clock('nlxc_vdw_energy',2)

       call internal_deallocate_vdwdf

    case('VDWDF2')
       !----------------------------------------------------------------------!
       ! vdW-DF2 functional                                                   !
       !----------------------------------------------------------------------!

       call internal_allocate_vdwdf

       ! JA: Changed xc_revpbe_x_point_sp -> xc_rpw86_x_point_sp.
       !     This matches the spec & the non spin pol version.
       ! Calculate the LDA correlation potential and GGA exchange potential
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, exchange,x_pot_fine,xc_rpw86_x_point,&
            xc_rpw86_x_point_sp,xc_none_c_point,xc_none_c_point_sp,&
            parallel=loc_parallel)
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
            xc_pw92_c_point, xc_pw92_c_point_sp,correlation=correlation, &
            c_pot_fine=c_pot_fine,parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work,&
            grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,conditioned_density_fine, &
            density_grad, grid,cell,'VDWDF2')

       xc_energy = exchange + correlation + nl_energy
       xc_pot_fine = x_pot_fine + c_pot_fine + nlxc_pot_fine

       call timer_clock('nlxc_vdw_energy',2)

       call internal_deallocate_vdwdf

   case('VDWC09')
       !----------------------------------------------------------------------!
       ! vdW-DF-C09x functional                                               !
       !----------------------------------------------------------------------!

       call internal_allocate_vdwdf


       ! Calculate the LDA correlation potential and GGA exchange potential
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, exchange,x_pot_fine,xc_c09_x_point,&
            xc_c09_x_point_sp,xc_none_c_point,xc_none_c_point_sp,&
            parallel=loc_parallel)
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
            xc_pw92_c_point, xc_pw92_c_point_sp,correlation=correlation, &
            c_pot_fine=c_pot_fine,parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work, &
            grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,conditioned_density_fine, &
            density_grad, grid,cell,'VDWDF1')
       xc_energy = exchange + correlation + nl_energy
       xc_pot_fine = x_pot_fine + c_pot_fine + nlxc_pot_fine

       call timer_clock('nlxc_vdw_energy',2)

       call internal_deallocate_vdwdf

   case('C09-2')
       !----------------------------------------------------------------------!
       ! vdW-DF-C09x functional with vdW-DF2 non-local correlation            !
       !----------------------------------------------------------------------!

       call internal_allocate_vdwdf


       ! Calculate the LDA correlation potential and GGA exchange potential
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, exchange,x_pot_fine,xc_c09_x_point,&
            xc_c09_x_point_sp,xc_none_c_point,xc_none_c_point_sp,&
            parallel=loc_parallel)
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
            xc_pw92_c_point, xc_pw92_c_point_sp,correlation=correlation, &
            c_pot_fine=c_pot_fine,parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work, &
            grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,conditioned_density_fine, &
            density_grad, grid,cell,'VDWDF2')
       xc_energy = exchange + correlation + nl_energy
       xc_pot_fine = x_pot_fine + c_pot_fine + nlxc_pot_fine

       call timer_clock('nlxc_vdw_energy',2)

       call internal_deallocate_vdwdf

    case('OPTPBE')
       !----------------------------------------------------------------------!
       ! vdW-DF-optPBE functional                                             !
       !----------------------------------------------------------------------!

       call internal_allocate_vdwdf

       ! Calculate the LDA correlation potential and GGA exchange potential
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, exchange,x_pot_fine,xc_vdwoptpbe_x_point,&
            xc_vdwoptpbe_x_point_sp,xc_none_c_point,xc_none_c_point_sp,&
            parallel=loc_parallel)
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
            xc_pw92_c_point, xc_pw92_c_point_sp,correlation=correlation, &
            c_pot_fine=c_pot_fine,parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work, &
            grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,conditioned_density_fine, &
            density_grad, grid,cell,'VDWDF1')
       xc_energy = exchange + correlation + nl_energy
       xc_pot_fine = x_pot_fine + c_pot_fine + nlxc_pot_fine

       call timer_clock('nlxc_vdw_energy',2)

       call internal_deallocate_vdwdf

    case('OPTB88')
       !----------------------------------------------------------------------!
       ! vdW-DF-optB88 functional                                             !
       !----------------------------------------------------------------------!

       call internal_allocate_vdwdf

       !Calculate the LDA correlation potential and GGA exchange potential
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, exchange,x_pot_fine,xc_vdwoptb88_x_point,&
            xc_vdwoptb88_x_point_sp,xc_none_c_point,xc_none_c_point_sp,&
            parallel=loc_parallel)
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
            xc_pw92_c_point, xc_pw92_c_point_sp,correlation=correlation, &
            c_pot_fine=c_pot_fine,parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work, &
            grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,conditioned_density_fine, &
            density_grad, grid,cell,'VDWDF1')
       xc_energy = exchange + correlation + nl_energy
       xc_pot_fine = x_pot_fine + c_pot_fine + nlxc_pot_fine

       call timer_clock('nlxc_vdw_energy',2)

       call internal_deallocate_vdwdf

   case('PBEK1')
       !----------------------------------------------------------------------!
       ! vdW-DF-optPBEK=1 functional                                          !
       !----------------------------------------------------------------------!

       call internal_allocate_vdwdf

       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid, exchange,x_pot_fine,xc_vdwpbek1_x_point,&
            xc_vdwpbek1_x_point_sp,xc_none_c_point,xc_none_c_point_sp,&
            parallel=loc_parallel)
       call xc_lda(conditioned_density_fine,xc_energy, xc_pot_fine, grid, &
            xc_pw92_c_point, xc_pw92_c_point_sp,correlation=correlation, &
            c_pot_fine=c_pot_fine,parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work, &
            grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,conditioned_density_fine, &
            density_grad, grid,cell,'VDWDF1')
       xc_energy = exchange + correlation + nl_energy
       xc_pot_fine = x_pot_fine + c_pot_fine + nlxc_pot_fine

       call timer_clock('nlxc_vdw_energy',2)

       call internal_deallocate_vdwdf

    case('VV10')
       !----------------------------------------------------------------------!
       ! rVV10 functional, as implemented by Sabatini et al.                  !
       !----------------------------------------------------------------------!

       call internal_allocate_vdwdf

       ! JA: Changed xc_revpbe_x_point_sp -> xc_rpw86_x_point_sp.
       !     This matches the spec & the non spin pol version.
       !Calculate the GGA exchange/correlation potential
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid,xc_energy,xc_pot_fine,xc_rpw86_x_point,&
            xc_rpw86_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp,&
            parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work, &
            grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,conditioned_density_fine, &
            density_grad, grid,cell,'VV10')
       xc_energy = xc_energy + nl_energy
       xc_pot_fine = xc_pot_fine + nlxc_pot_fine

       call timer_clock('nlxc_vdw_energy',2)

       call internal_deallocate_vdwdf

    case('VV10BC')
       !----------------------------------------------------------------------!
       ! rVV10 functional, as implemented by Sabatini et al., but with b and  !
       ! C parameters from Vydrov and Van Voorhis original paper, i.e. b=5.9, !
       ! C=0.0093.                                                            !
       !----------------------------------------------------------------------!

       call internal_allocate_vdwdf
       ! JA: Changed xc_revpbe_x_point_sp -> xc_rpw86_x_point_sp.
       !     This matches the spec & the non spin pol version.

       !Calculate the GGA exchange/correlation potential
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid,xc_energy,xc_pot_fine,xc_rpw86_x_point,&
            xc_rpw86_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp,&
            parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work, &
            grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,conditioned_density_fine, &
            density_grad, grid,cell,'VV10BC')
       xc_energy = xc_energy + nl_energy
       xc_pot_fine = xc_pot_fine + nlxc_pot_fine

       call timer_clock('nlxc_vdw_energy',2)

       call internal_deallocate_vdwdf

    case('AVV10S')
       !----------------------------------------------------------------------!
       ! AM05-VV10sol functional                                              !
       !----------------------------------------------------------------------!
       call internal_allocate_vdwdf

       !Calculate the GGA exchange/correlation potential
       call xc_gc(conditioned_density_fine,density_grad,density_aux,recip_work,&
            grid,xc_energy,xc_pot_fine,xc_am05_x_point,&
            xc_am05_x_point_sp,xc_am05_c_point,xc_am05_c_point_sp,&
            parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)
       call xc_gradients(conditioned_density_fine,density_grad,recip_work, &
            grid,0)
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,conditioned_density_fine, &
            density_grad, grid,cell,'VV10-S')
       xc_energy = xc_energy + nl_energy
       xc_pot_fine = xc_pot_fine + nlxc_pot_fine

       call timer_clock('nlxc_vdw_energy',2)

       call internal_deallocate_vdwdf

    !!!!!!!!!!!!!!!! KE-DENSITY-DEPENDENT META-GGA FUNCTIONALS !!!!!!!!!!!!!!!!!
    case('LTA')
       !----------------------------------------------------------------------!
       ! Local tau functional (JCW)                                           !
       ! Ernzerhof, M. & Scuseria, G. E. JCP 111, 911-915 (1999).             !
       ! A simple kinetic-energy-density-dependent functional.                !
       ! WARNING: EXPERIMENTAL CODE                                           !
       !----------------------------------------------------------------------!
       ! Kinetic-energy-dependent functional: requires ke_density_fine to be
       ! passed in.
       if (pub_on_root) then
          if (.not.(present(ke_density_fine).and.present(dfdtau_fine))) then
             call utils_abort("xc_mod: LTA functional requires &
                  &ke_density_fine and dfdtau_fine to be passed to &
                  &xc_energy_potential.")
          end if
          if (.not.ke_density_required) then
             call utils_abort("xc_mod: LTA functional requires &
                  &ke_density_required to be .true.")
          end if
       end if
       call xc_mgga(conditioned_density_fine,density_grad,&
            conditioned_ke_density_fine,&
            recip_work,grid,xc_energy,xc_pot_fine,dfdtau_fine,&
            xc_lta_x_point,xc_lta_x_point_sp,&
            xc_none_mgga_c_point,xc_none_mgga_c_point_sp,&
            parallel=loc_parallel)
    case('PKZB')
       !----------------------------------------------------------------------!
       ! PKZB meta-GGA, exchange + correlation (JCW)                          !
       !    JP Perdew, S Kurth, A Zupan & P Blaha,                            !
       !    Phys. Rev. Lett. 82, 2544 (1999)                                  !
       !----------------------------------------------------------------------!
       if (.not.(present(ke_density_fine).and.present(dfdtau_fine))) then
          call utils_abort("xc_mod: PKZB functional requires &
               &ke_density_fine and dfdtau_fine to be passed to &
               &xc_energy_potential.")
       end if
       if (.not.ke_density_required) then
          call utils_abort("xc_mod: PKZB functional requires &
               &ke_density_required to be .true.")
       end if
       call xc_mgga(conditioned_density_fine,density_grad,&
            conditioned_ke_density_fine,&
            recip_work,grid,xc_energy,xc_pot_fine,dfdtau_fine,&
            xc_pkzb_x_point,xc_pkzb_x_point_sp,&
            xc_pkzb_c_point,xc_pkzb_c_point_sp,parallel=loc_parallel)
            !xc_none_mgga_c_point,xc_none_mgga_c_point_sp)
    case('PKZB-C')
       !----------------------------------------------------------------------!
       ! PKZB correlation only (JCW)                                          !
       !    JP Perdew, S Kurth, A Zupan & P Blaha,                            !
       !    Phys. Rev. Lett. 82, 2544 (1999)                                  !
       !----------------------------------------------------------------------!
       if (.not.(present(ke_density_fine).and.present(dfdtau_fine))) then
          call utils_abort("xc_mod: PKZB functional requires &
               &ke_density_fine and dfdtau_fine to be passed to &
               &xc_energy_potential.")
       end if
       if (.not.ke_density_required) then
          call utils_abort("xc_mod: PKZB functional requires &
               &ke_density_required to be .true.")
       end if
       call xc_mgga(conditioned_density_fine,density_grad,&
            conditioned_ke_density_fine,&
            recip_work,grid,xc_energy,xc_pot_fine,dfdtau_fine,&
            xc_none_mgga_x_point,xc_none_mgga_x_point_sp,&
            xc_pkzb_c_point,xc_pkzb_c_point_sp,parallel=loc_parallel)
    case('PKZB-X')
       !----------------------------------------------------------------------!
       ! PKZB exchange only (JCW)                                             !
       !    JP Perdew, S Kurth, A Zupan & P Blaha,                            !
       !    Phys. Rev. Lett. 82, 2544 (1999)                                  !
       !----------------------------------------------------------------------!
       if (.not.(present(ke_density_fine).and.present(dfdtau_fine))) then
          call utils_abort("xc_mod: PKZB functional requires &
               &ke_density_fine and dfdtau_fine to be passed to &
               &xc_energy_potential.")
       end if
       if (.not.ke_density_required) then
          call utils_abort("xc_mod: PKZB functional requires &
               &ke_density_required to be .true.")
       end if
       call xc_mgga(conditioned_density_fine,density_grad,&
            conditioned_ke_density_fine,&
            recip_work,grid,xc_energy,xc_pot_fine,dfdtau_fine,&
            xc_pkzb_x_point,xc_pkzb_x_point_sp,&
            xc_none_mgga_c_point,xc_none_mgga_c_point_sp,parallel=loc_parallel)
    !!!!!!!!!!!!!!!! KE-DENSITY-DEPENDENT META-GGA FUNCTIONALS !!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!! WITH NON-LOCAL COMPONENT !!!!!!!!!!!!!!!!!!!!!!!!!!!
    case('B97M-V')
       !----------------------------------------------------------------------!
       ! B97M-V functional, exchange + correlation (JCW)                      !
       ! This uses rVV10 non-local correlation, as implemented by Sabatini    !
       ! et al., but with parameters b = 6.0_DP, C = 0.01_DP, as described in !
       !       N. Mardirossian & M. Head-Gordon,                              !
       !            J. Chem. Phys. 142, 074111 (2015)                         !
       ! See nlxc module for further details regarding NLC component.         !
       !----------------------------------------------------------------------!
       ! Allocate dynamic memory for non-local potential
       call internal_allocate_vdwdf

       if (.not.(present(ke_density_fine).and.present(dfdtau_fine))) then
          call utils_abort("xc_mod: B97M-V functional requires &
               &ke_density_fine and dfdtau_fine to be passed to &
               &xc_energy_potential.")
       end if
       if (.not.ke_density_required) then
          call utils_abort("xc_mod: B97M-V functional requires &
               &ke_density_required to be .true.")
       end if
       call xc_mgga(conditioned_density_fine,density_grad,&
            conditioned_ke_density_fine,&
            recip_work,grid,xc_energy,xc_pot_fine,dfdtau_fine,&
            xc_b97mv_x_point,xc_b97mv_x_point_sp,&
            xc_b97mv_c_point,xc_b97mv_c_point_sp,parallel=loc_parallel)

       ! Calculate nonlocal energy and potential
       call timer_clock('nlxc_vdw_energy',1)

       ! xc_gradients not required for xc_mgga, since density_grad is
       ! not changed by xc_mgga (unlike xc_gc).
       call nlxc_vdw_energy(nl_energy,nlxc_pot_fine,conditioned_density_fine, &
            density_grad, grid,cell,'VV10-B97M-V')
       xc_energy = xc_energy + nl_energy
       xc_pot_fine = xc_pot_fine + nlxc_pot_fine

       call timer_clock('nlxc_vdw_energy',2)

       ! Deallocate dynamic memory for non-local potential
       call internal_deallocate_vdwdf

    case default
       !----------------------------------------------------------------------!
       ! Unrecognised functional type                                         !
       !----------------------------------------------------------------------!
       if (pub_on_root) write(stdout,'(2a)') 'WARNING: unknown exchange-&
            &correlation functional: ', trim(loc_functional)
       xc_energy = 0.0_DP
       xc_pot_fine = 0.0_DP
    end select

    ! Free workspace if required
#ifdef LIBXC
    if (pub_xc_gradient_corrected.or.use_libxc) &
         call internal_deallocate_workspace
#else
    if (pub_xc_gradient_corrected) call internal_deallocate_workspace
#endif
    ! JCW: Deallocate conditioned version of ke_density_fine. Always allocated
    ! JCW: since we allocate a single element dummy array if ke_density_required
    ! JCW: is .false.
    deallocate(conditioned_ke_density_fine,stat=ierr)
    call utils_dealloc_check('xc_energy_potential','conditioned_ke_density_fine',&
         ierr)
    deallocate(conditioned_density_fine,stat=ierr)
    call utils_dealloc_check('xc_energy_potential','conditioned_density_fine', &
         ierr)

    ! Stop timer
    call timer_clock('xc_energy_potential',2)

  contains

    subroutine internal_allocate_workspace

      !=========================================================================!
      ! This subroutine allocates workspace for the gradients of the density if !
      ! required.                                                               !
      !-------------------------------------------------------------------------!
      ! Originally written by Peter Haynes in February 2004.                    !
      ! Modified by Quintin Hill in November 2008 to check allocs.              !
      ! Made subroutine-internal by Nicholas Hine, October 2010.                !
      !=========================================================================!

      use rundat, only: pub_num_spins

      implicit none

      ! Local Variables
      integer :: ierr    ! Error flag

      ! Allocate workspace for calculating gradients of the density
      allocate(density_grad(grid%ld1,grid%ld2,grid%max_slabs12,3, &
           pub_num_spins),stat=ierr)
      call utils_alloc_check('internal_allocate_workspace (xc)', &
           'density_grad',ierr)

      ! Allocate workspace for the auxiliary density array
      allocate(density_aux(grid%ld1,grid%ld2,grid%max_slabs12, &
           pub_num_spins),stat=ierr)
      call utils_alloc_check('internal_allocate_workspace (xc)', &
           'density_aux',ierr)

      ! Allocate reciprocal space workspace
      allocate(recip_work(grid%ld3,grid%ld2,grid%max_slabs23,3),stat=ierr)
      call utils_alloc_check('internal_allocate_workspace (xc)', &
           'recip_work',ierr)

    end subroutine internal_allocate_workspace

    subroutine internal_deallocate_workspace

      !========================================================================!
      ! This subroutine frees workspace for the gradients of the density.      !
      !------------------------------------------------------------------------!
      ! Originally written by Peter Haynes in February 2004.                   !
      ! Modified by Quintin Hill in November 2008 to use dealloc check.        !
      ! Made subroutine-internal by Nicholas Hine, October 2010.               !
      !========================================================================!

      implicit none

      ! Local Variables
      integer :: ierr    ! Error flag

      ! Deallocate reciprocal space workspace
      if (allocated(recip_work)) then
         deallocate(recip_work,stat=ierr)
         call utils_dealloc_check('internal_deallocate_workspace (xc)', &
              'recip_work', ierr)
      end if

      ! Deallocate workspace for calculating gradients of the density
      if (allocated(density_aux)) then
         deallocate(density_aux,stat=ierr)
         call utils_dealloc_check('internal_deallocate_workspace (xc)', &
              'density_aux', ierr)
      end if

      ! Deallocate workspace for calculating gradients of the density
      if (allocated(density_grad)) then
         deallocate(density_grad,stat=ierr)
         call utils_dealloc_check('internal_deallocate_workspace (xc)', &
              'density_grad', ierr)
      end if

    end subroutine internal_deallocate_workspace

    subroutine internal_allocate_vdwdf

      !========================================================================!
      ! This subroutine allocates variables used for the vdW-DF functionals    !
      !------------------------------------------------------------------------!
      !                                                                        !
      !========================================================================!

       use rundat, only: pub_num_spins

       implicit none

       integer :: ierr

       ! Allocate storage for nonlocal potential calculation
       allocate(nlxc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('xc_energy_potential','nlxc_pot_fine',ierr)
       allocate(x_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12, &
            pub_num_spins),stat=ierr)
       call utils_alloc_check('xc_energy_potential','x_pot_fine',ierr)
       allocate(c_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12, &
            pub_num_spins),stat=ierr)
       call utils_alloc_check('xc_energy_potential','c_pot_fine',ierr)

    end subroutine internal_allocate_vdwdf

    subroutine internal_deallocate_vdwdf

      !========================================================================!
      ! This subroutine deallocates variables used for the vdW-DF functionals  !
      !------------------------------------------------------------------------!
      !                                                                        !
      !========================================================================!

       implicit none

       integer :: ierr

       ! Deallocate storage for nonlocal potential calculation
       deallocate(c_pot_fine,stat=ierr)
       call utils_dealloc_check('xc_energy_potential','c_pot_fine',ierr)
       deallocate(x_pot_fine,stat=ierr)
       call utils_dealloc_check('xc_energy_potential','x_pot_fine',ierr)
       deallocate(nlxc_pot_fine,stat=ierr)
       call utils_dealloc_check('xc_energy_potential','nlxc_pot_fine',ierr)

    end subroutine internal_deallocate_vdwdf

    subroutine internal_check_density(parallel)

      !------------------------------------------------------------------------!
      ! This subroutine prints a warning if the charge density of either spin  !
      ! component is negative at any point.                                    !
      !------------------------------------------------------------------------!

      use comms, only: comms_reduce, pub_on_root
      use rundat, only: pub_num_spins

      implicit none

      logical, optional, intent(in) :: parallel

      !------------------------
      ! Local variables:
      !------------------------

      integer :: i1,i2,islab12 ! Loop counters over real-space fine grid
      integer :: is            ! Loop counter over spin
      real(kind=DP) :: denmin  ! Minimum values of density
      logical :: loc_parallel

      loc_parallel = .true.
      if(present(parallel)) then
         loc_parallel = parallel
      end if


      spins: do is=1, pub_num_spins

         denmin = density_fine(1,1,1,is)

         a3: do islab12=1,grid%num_my_slabs12
            a2: do i2=1,grid%n2
               a1: do i1=1,grid%n1

                  denmin = min(denmin,density_fine(i1,i2,islab12,is))

               end do a1
            end do a2
         end do a3

         if(loc_parallel) then
            call comms_reduce('MIN',denmin)
         end if

         if (pub_num_spins == 1) then
            if (denmin <= -0.0001_DP .and. pub_on_root) &
                 write(stdout,'(a,f7.4)') 'WARNING: minimum value of charge &
                 &density: ',denmin
         else
            if (denmin <= -0.0001_DP .and. pub_on_root) &
                 write(stdout,'(a,i1,a,f7.4)') 'WARNING: minimum value of charge&
                 & density for spin ',is,': ',denmin
         end if

      end do spins

    end subroutine internal_check_density

  end subroutine xc_energy_potential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_init(use_tddftxc,use_activexc)

    !=========================================================================!
    ! This subroutine initialises the internal variables in the module.       !
    !-------------------------------------------------------------------------!
    ! Originally written by Peter Haynes in February 2004.                    !
    ! Modified by Quintin Hill to remove unnecessary internal uppercase       !
    ! subroutine in November 2008.                                            !
    ! Modified by Quintin Hill to cope with additional functionals and remove !
    ! mention of nspins.                                                      !
    ! Modified by Quintin Hill on 06/04/2009 to remove unnecessary            !
    ! initialisations, workspace allocation is already done elsewhere.        !
    ! Modified for TDDFT by David D. O'Regan in May 2009.                     !
    ! Modified to be called at start only by Nicholas Hine in October 2010.   !
    ! Modified by James C. Womack in 2015/2016 to support KE-density-dependent!
    ! XC functionals (meta-GGAs).                                             !
    ! Modified by James C. Womack in 2015/2016 to support Libxc F2003         !
    ! interface introduced in version 3.0.0 of the library.
    ! Modified by James C. Womack in 2015/2016 to support calling KE-density- !
    ! dependent meta-GGAs from Libxc.                                         !
    !=========================================================================!

#ifdef LIBXC
    use iso_c_binding, only: c_int, c_double
#endif
    use bibliography, only: bibliography_cite
    ! JCW: comms pub_on_root required to print message only on root proc
    use comms, only: pub_on_root
    ! JCW: constants stdout used for temporary meta-GGA warning message
    use constants, only: stdout
    use utils, only: utils_abort, utils_assert, utils_use_var
    use rundat, only: pub_xc_functional, pub_tddft_xc_functional, &
         pub_active_xc_functional, pub_debug_on_root, &
         pub_libxc_x_func_id, pub_libxc_c_func_id, pub_aug, &
         pub_xc_min_tau, pub_xc_initial_functional, &
         pub_xc_ke_density_required, pub_xc_lapl_density_required, &
         pub_cdft, pub_cond_calculate, &              ! JCW: MGGA compatibility checks
         pub_dbl_grid_scale, pub_do_properties, &     ! JCW: MGGA compatibility checks
         pub_eda, pub_edft, &                         ! JCW: MGGA compatibility checks
         pub_hubbard, &                               ! JCW: MGGA compatibility checks
         pub_kernel_diis, pub_lr_phonons_calculate, & ! JCW: MGGA compatibility checks
         pub_lr_tddft_calculate, pub_maxit_pen, &         ! JCW: MGGA compatibility checks
         pub_pol_emb_pot, &                           ! JCW: MGGA compatibility checks
         pub_spin_polarised, pub_use_aux_ngwfs, &     ! JCW: MGGA compatibility checks
#ifdef LIBXC
         pub_paw,pub_num_spins
    use xc_f03_lib_m !! External dependency
#else
         pub_paw
#endif

    implicit none

    ! Arguments
    ! ddor: Switch between DFT xc functional and TDDFT xc functional
    logical, intent(in), optional :: use_tddftxc

    ! jcap: Switch between full system functional and active subregion
    ! functional in EMFT
    logical, intent(in), optional :: use_activexc

    ! Local Variables
    ! la: vdW-DF functional added to the gcfunctionals, as it requires the
    ! la: density gradient for its calculation
    ! JCW: meta-GGAs should be included in this list, since they also require
    ! JCW: the gradient of the density
    character(len=6), parameter :: gcfunctionals(31) &
         = (/'BLYP  ','PBE   ','PW91  ','REVPBE','RPBE  ','PBESOL','XLYP  ',&
         'WC    ','B1LYP ','B1PW91', 'B3LYP ','B3PW91','X3LYP ','PBE0  ',&
         'PBEX  ','VDWDF ','VDWDF2','OPTPBE','OPTB88','PBEK1 ','VDWC09', &
         'C09-2 ','VV10  ','VV10BC', 'AVV10S','AM05  ',& ! meta-GGA functionals follow...
         &'LTA   ','PKZB  ','PKZB-C','PKZB-X','B97M-V' /)
    ! JCW: Functionals requiring evaluation of kinetic energy density
    character(len=6), parameter :: ke_density_functionals(5) &
          = (/'LTA   ','PKZB  ','PKZB-C','PKZB-X','B97M-V' /)
    ! JCW: Functionals requiring evaluation of Laplacian of density
    character(len=6), parameter :: lapl_density_functionals(1) &
          = (/'######'/) ! placeholder
#ifdef LIBXC
    ! JCW: Recognized Libxc functional IDs for Laplacian of density dependent
    ! JCW: meta-GGAs
    integer(kind=c_int), parameter :: lapl_density_functionals_libxc_ids(1) &
          = (/ -123456789 /) ! placeholder
    ! JCW: Recognized Libxc functional IDs for KE-dependent meta-GGAs
    integer(kind=c_int), parameter :: ke_density_functionals_libxc_ids(14) &
          = (/ XC_MGGA_X_LTA, &
               XC_MGGA_X_PKZB,  XC_MGGA_C_PKZB, &
               XC_MGGA_X_TPSS,  XC_MGGA_C_TPSS, &
               XC_MGGA_X_REVTPSS, XC_MGGA_C_REVTPSS, &
               XC_MGGA_X_M06_L, XC_MGGA_C_M06_L, &
               XC_MGGA_X_M11_L, XC_MGGA_C_M11_L, &
               XC_MGGA_C_VSXC,  XC_MGGA_C_BC95, &
               XC_MGGA_XC_B97M_V /)
    ! JCW: List of Libxc functional IDs which use non-local correlation
    ! JCW: These functionals will use the nlxc module. Currently only
    ! JCW: the VV10 NLC functional is supported, with mixing coefficients
    ! JCW: passed to the nlxc_vdw_energy subroutine.
    integer(kind=c_int), parameter :: nlc_vv10_functionals_libxc_ids(2) &
         = (/ XC_GGA_XC_VV10, XC_MGGA_XC_B97M_V /)
    integer :: libxc_family
    integer :: libxc_polarized
    ! JCW: Temporary local variables with kind=c_double to satisfy the
    ! JCW: interface of xc_f03_nlc_coef
    real(kind=c_double) :: bval_tmp
    real(kind=c_double) :: Cval_tmp
#endif

    ! -------------------------------------------------------------------------

    call utils_use_var(pub_libxc_x_func_id)
    call utils_use_var(pub_libxc_c_func_id)

    ! Make a local copy of the xc_functional parameter from the input file

    ! jcap: switch between xc_functional and active_xc_functional
    ! jcap: activexc switch overrides tddft switch
    if (present(use_activexc)) then
       if (use_activexc) then
          functional = pub_active_xc_functional
       else
          ! ddor: loc_tddftxc switches between functionals defined in
          !       xc_functional and tddft_xc_functional
          if (present(use_tddftxc)) then
             if (use_tddftxc) then
                functional = pub_tddft_xc_functional
             else
                functional = pub_xc_functional
             end if
          else
             functional = pub_xc_functional
          endif
       end if
    else
       ! ddor: loc_tddftxc switches between functionals defined in
       !       xc_functional and tddft_xc_functional
       if (present(use_tddftxc)) then
          if (use_tddftxc) then
             functional = pub_tddft_xc_functional
          else
             functional = pub_xc_functional
          end if
       else
          functional = pub_xc_functional
       endif
    end if

    ! Apply defaults to functional type
    if (functional == 'LDA') then
       if (.not.pub_paw) then
          functional = 'CAPZ'
       else
          functional = 'PW92'
       end if
    end if
    if (functional == 'GGA') functional = 'RPBE'

    ! JCW: Print warning about incomplete meta-GGA implementation
    if ( any(functional == ke_density_functionals).and.pub_on_root) then
       ! Attempt to perform calculation with meta-GGA
       write(stdout,'(1x,9a,a)') "WARNING: ","meta-GGA functionality under development."
       write(stdout,'(1x,9x,a)')             "Using this functional is experimental."
       write(stdout,'(1x,9x,a,es16.8)')      "Minimum tau threshold, xc_min_tau = ", &
            pub_xc_min_tau
    end if

    ! Cite appropriate paper in bibliography
    select case (functional)
    case ('CAPZ')
       call bibliography_cite('PZ-LDA')
    case ('PBE')
       call bibliography_cite('PBE')
    end select

    ! ndmh: initialisations for LIBXC 'functional'
    if (functional == 'LIBXC') then

#ifdef LIBXC
       use_libxc = .true.

       ! ndmh: determined polarization
       if (pub_num_spins==1) then
          libxc_polarized = int(XC_UNPOLARIZED,kind=1)
       else
          libxc_polarized = int(XC_POLARIZED,kind=1)
       end if

       ! JCW: Determine whether NLC (VV10) is required
       libxc_nlc_vv10 = any(nlc_vv10_functionals_libxc_ids == &
            int(pub_libxc_c_func_id,kind=c_int) )

       if (pub_debug_on_root) write(stdout,'(a,l1)') &
            "DEBUG: xc_init, libxc_nlc_vv10 = ", libxc_nlc_vv10

       pub_xc_gradient_corrected    = .false. ! gradient corrected functional
       pub_xc_ke_density_required   = .false. ! KE density dependent functional
       pub_xc_lapl_density_required = .false. ! Laplacian of density dependent
                                              ! functional

       if (pub_libxc_x_func_id > 0) then
          call xc_f03_func_init(libxc_func(1), &
               pub_libxc_x_func_id, libxc_polarized)
          libxc_info(1) = xc_f03_func_get_info(libxc_func(1))
          libxc_family = xc_f03_func_info_get_family(libxc_info(1))
          libxc_gga = (int(libxc_family,kind=c_int)==XC_FAMILY_GGA) .or. &
               (int(libxc_family,kind=c_int)==XC_FAMILY_HYB_GGA)
          ! JCW: Additional logical variable for identifying meta-GGAs
          ! JCW: If an X functional is present (id /= 0) then C functional
          ! JCW: and X functional must both be meta-GGAs
          libxc_mgga = (int(libxc_family,kind=c_int)==XC_FAMILY_MGGA) .or. &
               (int(libxc_family,kind=c_int)==XC_FAMILY_HYB_MGGA)
       end if
       if (pub_libxc_c_func_id > 0) then
          call xc_f03_func_init(libxc_func(2), &
               pub_libxc_c_func_id, libxc_polarized)
          libxc_info(2) = xc_f03_func_get_info(libxc_func(2))
          libxc_family = xc_f03_func_info_get_family(libxc_info(2))
          if (pub_libxc_x_func_id == 0) then
             ! JCW: If a null X functional is used (pub_libxc_x_func_id == 0),
             ! JCW: then libxc_{m}gga is not yet set. Set libxc_{m}gga based on C
             ! JCW: functional only in this case.
             libxc_gga = (int(libxc_family,kind=c_int)==XC_FAMILY_GGA) .or. &
                  (int(libxc_family,kind=c_int)==XC_FAMILY_HYB_GGA)
             libxc_mgga = (int(libxc_family,kind=c_int)==XC_FAMILY_MGGA) .or. &
                  (int(libxc_family,kind=c_int)==XC_FAMILY_HYB_MGGA)
          end if
       end if

       if ((pub_libxc_x_func_id/=0).and.(pub_libxc_c_func_id/=0)) then
          ! JCW: Restrict combinations of exchange and correlation functionals
          ! JCW: such that only {meta-}GGA-type functionals may be combined
          ! JCW: (unless one or other of the functions is equal to 0 (null),
          ! JCW: which may be combined with a {meta-}GGA functional).
          if (libxc_gga.neqv.((int(libxc_family,kind=c_int)==XC_FAMILY_GGA) .or. &
               (int(libxc_family,kind=c_int)==XC_FAMILY_HYB_GGA))) then
             call utils_abort('Error in xc_init: Mixing of LIBXC functionals &
                  &with and without gradient corrections is not supported')
          end if
          if (libxc_mgga.neqv.(int(libxc_family,kind=c_int)==XC_FAMILY_MGGA)) then
             call utils_abort('Error in xc_init: Mixing of LIBXC functionals &
                 &that are meta-GGAs with other types not supported')
          end if
       else if ((pub_libxc_x_func_id==0).and.(pub_libxc_c_func_id==0 )) then
          call utils_abort('Error in xc_init: LIBXC X and C functionals &
                           &cannot both be null (func_id == 0)')
       end if

       ! JCW: If Libxc correlation functional includes VV10 non-local
       ! JCW: correlationthen set libxc_vv10_bval and libxc_vv10_Cval
       if (libxc_nlc_vv10) then
          call xc_f03_nlc_coef(libxc_func(2),bval_tmp,Cval_tmp)
          libxc_vv10_bval = real(bval_tmp, kind=DP)
          libxc_vv10_Cval = real(Cval_tmp, kind=DP)
       else
          libxc_vv10_bval = 0.0_DP
          libxc_vv10_Cval = 0.0_DP
       end if

       ! JCW: meta-GGAs may depend on KE density and/or Laplacian of density
       ! JCW: so we must determine which of these quantities needs evaluating.
       ! JCW: We assume that the gradient of the density is always required, in
       ! JCW: line with the conventional definition of meta-GGA functionals,
       ! JCW: see e.g. Perdew, J. P. et al. JCP 123, 062201 (2005).

       ! JCW: Check if functional requires KE density evaluation.
       ! JCW: We would like to be able to obtain this information from Libxc
       ! JCW: based on the functional id, but this is not currently supported
       ! JCW: by Libxc.
       ! JCW: Temporary alternative --- use hard-coded list
       pub_xc_ke_density_required = &
            (any(ke_density_functionals_libxc_ids == int(pub_libxc_x_func_id,kind=c_int)).or.&
            any(ke_density_functionals_libxc_ids == int(pub_libxc_c_func_id,kind=c_int)))

       ! JCW: Check if functional required Laplacian of density evaluation.
       ! JCW: As above, Libxc does not currently provide this information.
       ! JCW: Temporary alternative --- use hard-coded list
       pub_xc_lapl_density_required = &
            (any(lapl_density_functionals_libxc_ids == int(pub_libxc_x_func_id,kind=c_int)).or.&
            any(lapl_density_functionals_libxc_ids == int(pub_libxc_c_func_id,kind=c_int)))

       ! ndmh: check if this functional is gradient corrected
       pub_xc_gradient_corrected = libxc_gga.or.libxc_mgga

       ! JCW: if functional is a meta-GGA, then libxc_mgga should be true and
       ! JCW: libxc_gga should be false (the categories are exclusive).
       call utils_assert(.not.(libxc_gga.and.libxc_mgga),&
            "Error in xc_init: Functional must be either GGA or meta-GGA &
            &not both, i.e. libxc_gga and libxc_mgga cannot both be true")
       ! JCW: Abort if hybrid meta-GGA requested, since these have not been
       ! JCW: tested.
       call utils_assert(.not.(int(libxc_family,kind=c_int)==XC_FAMILY_HYB_MGGA),&
            "Error in xc_init: hybrid meta-GGAs have not yet been &
            &implemented.")
       ! JCW: If libxc_mgga is .true., then pub_libxc_x_func_id and
       ! JCW: pub_libxc_c_func_id must either be zero, or present in
       ! JCW: ke_density_functionals_libxc_ids, or
       ! JCW: lapl_density_functionals_libxc_ids
       if (libxc_mgga) then
          if (pub_libxc_x_func_id > 0) then
             call utils_assert(&
                  any(ke_density_functionals_libxc_ids == &
                  int(pub_libxc_x_func_id,kind=c_int)).or.&
                  any(lapl_density_functionals_libxc_ids == &
                  int(pub_libxc_x_func_id,kind=c_int)),&
                  "Error in xc_init: If a meta-GGA is requested from Libxc, &
                  &then the exchange functional must be dependent on the &
                  &kinetic energy density, the Laplacian of the density, or &
                  &both.")
          end if
          if (pub_libxc_c_func_id > 0) then
             call utils_assert(&
                  any(ke_density_functionals_libxc_ids == &
                  int(pub_libxc_c_func_id,kind=c_int)).or.&
                  any(lapl_density_functionals_libxc_ids == &
                  int(pub_libxc_c_func_id,kind=c_int)),&
                  "Error in xc_init: If a meta-GGA is requested from Libxc, &
                  &then the correlation functional must be dependent on the &
                  &kinetic energy density, the Laplacian of the density, or &
                  &both.")
          end if
       end if

       if (pub_debug_on_root) then
          write(stdout,'(a)') "DEBUG: xc_init -- Libxc functional information:"
          write(stdout,'(a,i5)') "DEBUG: pub_libxc_x_func_id = ", &
               pub_libxc_x_func_id
          write(stdout,'(a,i5)') "DEBUG: pub_libxc_c_func_id = ", &
               pub_libxc_c_func_id
          write(stdout,'(a,l3)') "DEBUG: pub_xc_gradient_corrected = ", &
               pub_xc_gradient_corrected
          write(stdout,'(a,l3)') "DEBUG: pub_xc_ke_density_required = ", &
               pub_xc_ke_density_required
          write(stdout,'(a,l3)') "DEBUG: pub_xc_lapl_density_required = ", &
               pub_xc_lapl_density_required
          write(stdout,'(a,l3)') "DEBUG: libxc_nlc_vv10 = ", libxc_nlc_vv10
          if (libxc_nlc_vv10) then
             write(stdout,'(a,es20.10)') "DEBUG: libxc_vv10_bval = ", &
                  libxc_vv10_bval
             write(stdout,'(a,es20.10)') "DEBUG: libxc_vv10_Cval = ", &
                  libxc_vv10_Cval
          end if
       end if

#else
       call utils_abort('Error in xc_init: Functional "LIBXC" &
            &specified, but LIBXC not present in build')
#endif

    else

#ifdef LIBXC
       use_libxc = .false.
#endif

       ! Determine whether functional requires gradients to be calculated
       pub_xc_gradient_corrected = (any(gcfunctionals == functional))

       ! JCW: meta-GGAs may depend on KE density and/or Laplacian of density
       ! JCW: so we must determine which of these quantities needs evaluating
       ! JCW: We assume that the gradient of the density is always required, in
       ! JCW: line with the conventional definition of meta-GGA functionals,
       ! JCW: see e.g. Perdew, J. P. et al. JCP 123, 062201 (2005).

       ! JCW: Determine whether functional requires kinetic energy density
       ! JCW: to be evaluated.
       pub_xc_ke_density_required = (any(ke_density_functionals == functional))

       ! JCW: Determine whether functional requires Laplacian of density to
       ! JCW: be evaluated
       pub_xc_lapl_density_required = (any(lapl_density_functionals == &
           functional))

       ! JCW: Laplacian of density not (yet) supported -- abort.
       call utils_assert(.not.pub_xc_lapl_density_required, &
            "Error in xc_init: Functionals dependent on the Laplacian of the &
            &density are not supported at present.")
    end if

    ! JCW: Set radial_functional and initial_functional
    if (pub_xc_ke_density_required) then
       ! If meta-GGA:
       if (pub_xc_initial_functional /= "NONE") then
          ! If xc_initial_functional keyword is set, use alternative
          ! functional for radial densities (in xc_radial)
          radial_functional = trim(pub_xc_initial_functional)
          ! JCW: Currently only GGAs are support as alternative radial functionals
          ! TODO Allow any type of non-meta-GGA to be used for xc_radial
          call utils_assert(any(gcfunctionals==radial_functional),&
               "Error in xc_init: &
               &XC_INITIAL_FUNCTIONAL must specify a non-LIBXC GGA &
               &functional to accompany a meta-GGA.")
          call utils_assert(.not.(any(ke_density_functionals==radial_functional)&
               .or.any(lapl_density_functionals==radial_functional)),&
               "Error in xc_init: &
               &XC_INITIAL_FUNCTIONAL cannot be a meta-GGA.")
          ! Print warning
          if (pub_on_root) then
             write(stdout,'(1x,9a,a)') "WARNING: ",&
                  "Using alternative XC functional "//trim(radial_functional)//" for "
             write(stdout,'(1x,9x,a)') "radial densities (xc_radial)."
          end if
       else
          ! If no xc_initial_functional has been specified, neglect exchange and
          ! correlation in xc_radial when a meta-GGA is used
          radial_functional = 'NONE'
          if (pub_on_root) then
             write(stdout,'(1x,9a,a)') "WARNING: ",&
                  "Neglecting exchange and correlation for radial densities."
          end if
       end if
    else ! pub_xc_ke_density_required == F
       ! If not meta-GGA:
       ! Use same functional as used elsewhere
       radial_functional = functional
    end if
    ! Initial functional is always the same as radial functional
    initial_functional = trim(radial_functional)

    if (pub_xc_gradient_corrected.and.pub_aug) then
       !call utils_abort('Error in xc_init: Augmentation charges with &
       !     &gradient-corrected functionals are not yet supported')
    end if
    if (pub_debug_on_root) then
       write(stdout,'(a)') "DEBUG: xc_init -- functional information:"
       write(stdout,'(a,i5)') "DEBUG: pub_libxc_x_func_id = ", &
            pub_libxc_x_func_id
       write(stdout,'(a,i5)') "DEBUG: pub_libxc_c_func_id = ", &
            pub_libxc_c_func_id
       write(stdout,'(a,l3)') "DEBUG: pub_xc_gradient_corrected = ", &
            pub_xc_gradient_corrected
       write(stdout,'(a,l3)') "DEBUG: pub_xc_ke_density_required = ", &
            pub_xc_ke_density_required
       write(stdout,'(a,l3)') "DEBUG: pub_xc_lapl_density_required = ", &
            pub_xc_lapl_density_required
       write(stdout,'(a,a)') "DEBUG: functional = ", &
            trim(functional)
       write(stdout,'(a,a)') "DEBUG: radial_functional = ", &
            trim(radial_functional)
       write(stdout,'(a,a)') "DEBUG: initial_functional = ", &
            trim(initial_functional)
    end if

    ! JCW: Meta-GGA (in)-compatibility checks
    ! JCW: These cannot be performed inside rundat_check_inputs, since
    ! JCW: we need to avoid accessing information from xc module in rundat
    ! JCW: to avoid circular dependencies between modules.
    if (pub_xc_ke_density_required) then
       ! JCW: tau-dependent XC functional requested, perform compatibility
       ! JCW: checks and abort if compatibility unimplemented or untested.

       ! JCW: Double grid scale <= 1.0_DP not supported
       call utils_assert(pub_dbl_grid_scale>1.0_DP,"Error in xc_init: &
            &pub_xc_ke_density_required is true, but combination of &
            &KE-density-dependent XC functionals with pub_dbl_grid_scale &
            &<= 1.0_DP is not supported.")

       ! JCW: Auxiliary NGWFs, not implemented/untested
       call utils_assert(.not.pub_use_aux_ngwfs,"Error in xc_init: &
            &pub_xc_ke_density_required is true, but combination of &
            &KE-density-dependent XC functionals with PAW &
            &has not been implemented.")

       ! JCW: Projector augmented waves (PAW), not implemented
       call utils_assert(.not.pub_paw,"Error in xc_init: &
            &pub_xc_ke_density_required is true, but combination of &
            &KE-density-dependent XC functionals with PAW &
            &has not been implemented.")

       ! JCW: Conduction NGWFs, not implemented/untested
       call utils_assert(.not.pub_cond_calculate,&
            "Error in xc_init: &
            &pub_xc_ke_density_required is true, but combination of &
            &KE-density-dependent XC functionals with valence NGWFs &
            &has not been implemented/tested.")

       ! JCW: Constrained DFT, not implemented/untested
       call utils_assert(.not.pub_cdft,"Error in xc_init: &
            &pub_xc_ke_density_required is true, but combination of &
            &KE-density-dependent XC functionals with constrained DFT &
            &has not been implemented/tested.")

       ! JCW: Ensemble DFT, not implemented/untested
       call utils_assert(.not.pub_edft,"Error in xc_init: &
            &pub_xc_ke_density_required is true, but combination of &
            &KE-density-dependent XC functionals with ensemble DFT &
            &has not been implemented/tested.")

       ! JCW: Hubbard SCF / DFT+U, not implemented/tested
       call utils_assert(.not.pub_hubbard,"Error in xc_init &
            &pub_xc_ke_density_required is true, but combination of &
            &KE-density-dependent XC functionals with Hubbard SCF / DFT+U &
            &has not been implemented/tested.")

       ! JCW: Kernel DIIS, not implemented/untested
       call utils_assert(.not.pub_kernel_diis,"Error in xc_init: &
            &pub_xc_ke_density_required is true, but combination of &
            &KE-density-dependent XC functionals with kernel DIIS &
            &has not been implemented/tested.")

       ! JCW: Linear-response phonons, not implemented/untested
       call utils_assert(.not.pub_lr_phonons_calculate,&
            "Error in xc_init: &
            &pub_xc_ke_density_required is true, but combination of &
            &KE-density-dependent XC functionals with linear-response phonons &
            &has not been implemented/tested.")

       ! JCW: LR-TDDFT, not implemented/untested
       call utils_assert(.not.pub_lr_tddft_calculate,&
            "Error in xc_init: &
            &pub_xc_ke_density_required is true, but combination of &
            &KE-density-dependent XC functionals with LR-TDDFT &
            &has not been implemented/tested.")

       ! JCW: Penalty functional (instead of LNV), not implemented/untested
       call utils_assert(pub_maxit_pen == 0,&
            "Error in xc_init: &
            &pub_xc_ke_density_required is true, but combination of &
            &KE-density-dependent XC functionals with penalty functional &
            &has not been implemented/tested.")

       ! JCW: Properties not yet supported
       call utils_assert(.not.pub_do_properties,"Error in xc_init: &
            &pub_xc_ke_density_required is true, but properties calculation &
            &has not yet been implemented or tested with KE-density-dependent &
            &XC functionals.")

!       ! JCW: Spin-polarisation not yet supported
!       call utils_assert(.not.pub_spin_polarised,"Error in xc_init: &
!            &pub_xc_ke_density_required is true, but spin-polarization &
!            &has not yet been implemented or tested with KE-density-dependent &
!            &XC functionals.")

    end if

  end subroutine xc_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_exit

    !=========================================================================!
    ! This subroutine deallocates the internal variables in the module.       !
    !-------------------------------------------------------------------------!
    ! Originally written by Nicholas Hine in October 2010.                    !
    ! Libxc interface updated for version 3 / F2003  by James C. Womack,      !
    ! 2015/2016                                                               !
    !=========================================================================!

    use rundat, only: pub_libxc_x_func_id, pub_libxc_c_func_id
    use utils, only: utils_abort, utils_assert, utils_use_var
#ifdef LIBXC
    use xc_f03_lib_m !! External dependency
#endif

    implicit none

#ifdef LIBXC
    ! ndmh: de-initialisations for LIBXC 'functional'
    if (use_libxc) then
       if (pub_libxc_x_func_id > 0) then
          call xc_f03_func_end(libxc_func(1))
       end if
       if (pub_libxc_c_func_id > 0) then
          call xc_f03_func_end(libxc_func(2))
       end if
    end if
#else
    call utils_use_var(pub_libxc_x_func_id)
    call utils_use_var(pub_libxc_c_func_id)
#endif

  end subroutine xc_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine xc_hfxinit(active_switch)

    !===========================================================!
    ! This subroutine sets pub_use_hfx and pub_hfxfraction.     !
    !-----------------------------------------------------------!
    ! Written by Quintin Hill in October 2008.                  !
    ! Quintin Hill added hybrid functionals on 11/03/2009.      !
    ! Edited to deal with active region functionals by Joseph   !
    ! Prentice, 23/03/18                                        !
    !===========================================================!

    ! jd: Any changes to the list of hybrid functionals below need to be refle-
    !     cted in xc::xc_fxc_finite_diff() and in rundat::internal_print_bc().

    use rundat, only: pub_xc_functional, rundat_set_pub_use_hfx, &
         pub_active_xc_functional, rundat_set_pub_use_activehfx

    implicit none

    logical, intent(in), optional :: active_switch

    logical :: loc_switch

    if (present(active_switch)) then
       loc_switch=active_switch
    else
       loc_switch=.false.
    end if

    select case (pub_xc_functional)

    case ('HF')
       call rundat_set_pub_use_hfx(.true.)
       if (.not.loc_switch) pub_hfxfraction = 1.0_DP
    case ('B1LYP','B1PW91','PBE0')
       call rundat_set_pub_use_hfx(.true.)
       if (.not.loc_switch) pub_hfxfraction = 0.25_DP
    case ('B3LYP','B3PW91')
       call rundat_set_pub_use_hfx(.true.)
       if (.not.loc_switch) pub_hfxfraction = 0.2_DP
    case ('X3LYP')
       call rundat_set_pub_use_hfx(.true.)
       if (.not.loc_switch) pub_hfxfraction = 0.218_DP
    case default
       call rundat_set_pub_use_hfx(.false.)
       if (.not.loc_switch) pub_hfxfraction = 0.0_DP
    end select

    select case (pub_active_xc_functional)

    case ('HF')
       call rundat_set_pub_use_activehfx(.true.)
       if (loc_switch) pub_hfxfraction = 1.0_DP
    case ('B1LYP','B1PW91','PBE0')
       call rundat_set_pub_use_activehfx(.true.)
       if (loc_switch) pub_hfxfraction = 0.25_DP
    case ('B3LYP','B3PW91')
       call rundat_set_pub_use_activehfx(.true.)
       if (loc_switch) pub_hfxfraction = 0.2_DP
    case ('X3LYP')
       call rundat_set_pub_use_activehfx(.true.)
       if (loc_switch) pub_hfxfraction = 0.218_DP
    case default
       call rundat_set_pub_use_activehfx(.false.)
       if (loc_switch) pub_hfxfraction = 0.0_DP
    end select

  end subroutine xc_hfxinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if 0
  subroutine xc_test(grid,cell,parallel)

    !=========================================================================!
    ! This subroutine tests routines associated with the module.              !
    !-------------------------------------------------------------------------!
    ! Originally written by Peter Haynes in February 2004.                    !
    ! Modified by Quintin Hill in November 2008 to use utils_(de)alloc_check. !
    ! Modified by Quintin Hill in March 2009 to use new subroutines and add   !
    ! test of BLYP.                                                           !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce, pub_on_root, pub_my_proc_id
    use constants, only: DP, stdout
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_banner
    use simulation_cell, only: CELL_INFO
    use xc_funcs

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    logical, optional, intent(in) :: parallel

    ! Local variables
    integer :: ierr                                     ! Error flag
    integer :: i1,i2,i3                                 ! Grid loop counters
    integer :: islab12                                  ! Slab loop counter
    integer :: m1,m2,m3                                 ! Grid midpoints
    integer :: npts                                     ! Number of points
    logical :: grad_alloc                               ! Flag for workspace
    real(kind=DP), parameter :: threshold = 0.005_DP    ! Threshold
    real(kind=DP) :: xc_energy                          ! XC energy
    real(kind=DP) :: xc_potint                          ! XC potential integral
    real(kind=DP) :: denint                             ! Density integral
    real(kind=DP) :: normfac                            ! Normalising factor
    real(kind=DP) :: r1,r2,r3                           ! Distances
    real(kind=DP) :: r3sq,r23sq,rsq                     ! Squared distances
    real(kind=DP) :: gamma                              ! Parameter for test
    real(kind=DP) :: grad(3)                            ! Gradient
    real(kind=DP) :: mod_grad                           ! Modulus of gradient
    real(kind=DP) :: ana_grad                           ! Analytic gradient
    real(kind=DP) :: num_grad                           ! Numerical gradient
    real(kind=DP) :: err,errsq,rmserr                   ! Errors in gradient
    real(kind=DP), allocatable :: density_fine(:,:,:,:) ! Density array
    real(kind=DP), allocatable :: xc_pot_fine(:,:,:,:)  ! Potential array
    logical :: loc_parallel

    loc_parallel = .true.
    if(present(parallel)) then
       loc_parallel = parallel
    end if

    if (pub_on_root) then
       write(stdout,'(/a/)') utils_banner('=', &
            'Exchange-correlation module test')
       write(stdout,'(a/)') '  N.B. Orthorhombic cells only!'
    end if

    ! Initialise workspace for gradients if necessary
    grad_alloc = allocated(density_grad)
    if (.not. grad_alloc) call xc_workspace(grid)

    ! Allocate workspace
    allocate(density_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('xc_test', 'density_fine', ierr)
    allocate(xc_pot_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('xc_test', 'xc_pot_fine', ierr)

    ! Calculate grid midpoints
    m1 = grid%n1 / 2
    m2 = grid%n2 / 2
    m3 = grid%n3 / 2

    ! Calculate optimum exponent
    gamma = 1.0_DP / sqrt(pi * min(grid%n1, &
         grid%n2, grid%n3))

    ! Create a Gaussian density to test the routines
    normfac = par%nat * 2 * (0.5_DP*gamma/pi)**1.5_DP
    denint = 0.0_DP
    do islab12=1,grid%num_my_slabs12
       i3=grid%first_slab12(pub_my_proc_id) + islab12 - 1
       r3 = (i3-m3) * 0.5_DP * cell%d3
       r3sq = r3*r3
       do i2=1,grid%n2
          r2 = (i2-m2) * 0.5_DP * cell%d2
          r23sq = r2*r2 + r3sq
          do i1=1,grid%n1
             r1 = (i1-m1) * 0.5_DP * cell%d1
             rsq = r1*r1 + r23sq
             density_fine(i1,i2,islab12,1) = normfac * &
                  exp(-0.5_DP * gamma * rsq)
             denint = denint + density_fine(i1,i2,i3,1)
          end do
       end do
    end do
    call comms_reduce('SUM',denint)
    denint = denint * grid%weight

    ! Calculate gradients of the density
    call xc_gradients(density_fine,grid)

    ! Compare gradients against analytic result
    npts = 0
    errsq = 0.0_DP
    do islab12=1,grid%num_my_slabs12
       i3=grid%first_slab12(pub_my_proc_id) + islab12 - 1
       r3 = (i3-m3) * 0.5_DP * cell%d3
       r3sq = r3*r3
       do i2=1,grid%n2
          r2 = (i2-m2) * 0.5_DP * cell%d2
          r23sq = r2*r2 + r3sq
          do i1=1,grid%n1
             r1 = (i1-m1) * 0.5_DP * cell%d1
             rsq = r1*r1 + r23sq
             if (gamma * rsq < threshold) then
                grad = density_grad(i1,i2,islab12,:,1)
                mod_grad = sqrt(grad(1)*grad(1) + grad(2)*grad(2) + &
                     grad(3)*grad(3))
                num_grad = mod_grad / density_fine(i1,i2,islab12,1)
                ana_grad = gamma * sqrt(rsq)
                err = num_grad - ana_grad
                errsq = errsq + err*err
                npts = npts + 1
             end if
          end do
       end do
    end do
    call comms_reduce('SUM',errsq)
    call comms_reduce('SUM',npts)
    rmserr = sqrt(errsq / npts)
    if (pub_on_root) then
       write(stdout,'(a)') '  Gradient calculation test:'
       write(stdout,'(a,f12.7)') '    Gaussian density variance  : ',&
            1.0_DP / gamma
       write(stdout,'(a,f12.7)') '    Gaussian density integral  : ',denint
       write(stdout,'(a,e12.4)') '    Error in density integral  : ', &
            abs(denint-2.0_DP*par%nat)
       write(stdout,'(a,e12.4)') '    RMS error in gradient      : ',rmserr
    end if

    ! Call appropriate routines for each functional type available:

    ! Ceperley-Alder LDA:
    call xc_lda(density_fine, xc_energy, xc_pot_fine, grid, xc_capz_c_point, &
         xc_capz_c_point_sp, xc_potint,parallel=loc_parallel)
    if (pub_on_root) then
       write(stdout,'(/a)') '  Ceperley-Alder LDA:'
       write(stdout,'(a,f16.8)') '    Exchange-correlation energy: ',xc_energy
       write(stdout,'(a,f16.8)') '    Potential-density integral : ',xc_potint
    end if

    ! Vosko-Wilk-Nusair LDA:
    call xc_lda(density_fine, xc_energy, xc_pot_fine, grid, xc_vwn_c_point, &
         xc_vwn_c_point_sp, xc_potint,parallel=loc_parallel)
    if (pub_on_root) then
       write(stdout,'(/a)') '  Vosko-Wilk-Nusair LDA:'
       write(stdout,'(a,f16.8)') '    Exchange-correlation energy: ',xc_energy
       write(stdout,'(a,f16.8)') '    Potential-density integral : ',xc_potint
    end if

    ! Perdew-Burke-Ernzerhof GGA:
    call xc_gc(density_fine, xc_energy, xc_pot_fine, grid, xc_pbe_x_point,&
         xc_pbe_x_point_sp, xc_pbe_c_point, xc_pbe_c_point_sp, xc_potint,&
         parallel=loc_parallel)
    if (pub_on_root) then
       write(stdout,'(/a)') '  Perdew-Burke-Ernzerhof GGA:'
       write(stdout,'(a,f16.8)') '    Exchange-correlation energy: ',xc_energy
       write(stdout,'(a,f16.8)') '    Potential-density integral : ',xc_potint
    end if

    ! Perdew-Wang 91 GGA:
    call xc_gc(density_fine, xc_energy, xc_pot_fine, grid, xc_pw91_x_point,&
         xc_pw91_x_point_sp, xc_pw91_c_point, xc_pw91_c_point_sp,xc_potint,&
         parallel=loc_parallel)
    if (pub_on_root) then
       write(stdout,'(/a)') '  Perdew-Wang ''91 GGA:'
       write(stdout,'(a,f16.8)') '    Exchange-correlation energy: ',xc_energy
       write(stdout,'(a,f16.8)') '    Potential-density integral : ',xc_potint
    end if

    ! BLYP GGA:
    call xc_gc(density_fine, xc_energy, xc_pot_fine, grid, xc_b88_x_point, &
         xc_b88_x_point_sp, xc_lyp_c_point, xc_lyp_c_point_sp,xc_potint,&
         parallel=loc_parallel)
    if (pub_on_root) then
       write(stdout,'(/a)') '  BLYP GGA:'
       write(stdout,'(a,f16.8)') '    Exchange-correlation energy: ',xc_energy
       write(stdout,'(a,f16.8)') '    Potential-density integral : ',xc_potint
    end if

    ! Deallocate workspace
    deallocate(xc_pot_fine,stat=ierr)
    call utils_dealloc_check('xc_test', 'xc_pot_fine', ierr)
    deallocate(density_fine,stat=ierr)
    call utils_dealloc_check('xc_test', 'density_fine', ierr)

    ! Free workspace if required
    if (.not. grad_alloc) call xc_free

    if (pub_on_root) write(stdout,'(/a/)') repeat('=',80)

  end subroutine xc_test
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_radial(npts,npts_max,ns,den_rad,grad_rad,vxc_rad, &
       dvxc_dmgd_rad,exc_rad,rad,rmax,rab,prt)

    !==========================================================================!
    ! This subroutine calculates the exchange-correlation energy and potential !
    ! of a radial density distribution, for atoms and PAW calculations.        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !    npts (in)     : Number of actual radial points used in arrays         !
    !    npts_max (in) : Dimension of radial grids used to store arrays        !
    !    rmax (in)     : Maximum of value on radial grid                       !
    !    ns (in)       : Number of spin channels                               !
    !    den_rad (in)  : Density on radial grid (no r^2 factor included)       !
    !    grad_rad (in) : Gradient dn/dr on radial grid (no r^2 factor)         !
    !    vxc_rad (in)  : Potential on radial grid                              !
    !    dvxc_dmgd_rad (inout) : Workspace for d(v_xc)/d|grad(n)|              !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in June 2010.                                   !
    ! Updated to allow XC fuctional used in xc_radial to differ from the XC    !
    ! functional used elsewhere, by James C. Womack, 2015/2016                 !
    ! Updated to use F2003 interface introduced in Libxc 3.0.0 by              !
    ! James C. Womack, 2015/2016                                               !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: max_spins, stdout
    use services, only: services_radial_derivative
    use utils, only: utils_abort
    use xc_funcs
    use rundat, only: pub_num_spins, pub_devel_code
    use utils, only: utils_devel_code

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    integer, intent(in) :: npts_max
    integer, intent(in) :: ns
    real(kind=DP), intent(in) :: den_rad(npts_max,ns)
    real(kind=DP), intent(inout) :: grad_rad(npts_max,ns)
    real(kind=DP), intent(out) :: vxc_rad(npts_max,ns)
    real(kind=DP), intent(out), optional :: dvxc_dmgd_rad(npts_max,ns)
    real(kind=DP), intent(out), optional :: exc_rad(npts_max,ns)
    real(kind=DP), intent(in), optional :: rad(npts_max)
    real(kind=DP), intent(in), optional :: rmax
    real(kind=DP), intent(in), optional :: rab(npts_max)
    logical, intent(in), optional :: prt

    ! Local Variables
    integer :: ipt, is
    real(kind=DP) :: x_energy, x_pot(max_spins)
    real(kind=DP) :: c_energy, c_pot(max_spins)
    logical :: calc_x, calc_c
    logical :: loc_prt

    if (index(pub_devel_code,'NO_XC_INIT')>0) then
       ! JCW: If devel_code: NO_XC_INIT then, skip XC potential and energy
       ! JCW: evaluation in initial guess
       if (pub_on_root) write(stdout,*) "WARNING: NO_XC_INIT, &
          &xc_radial setting vxc_rad = exc_rad = dvxc_dmgd_rad = 0.0_DP"
       vxc_rad = 0.0_DP
       if (present(exc_rad)) exc_rad = 0.0_DP
       if (present(dvxc_dmgd_rad)) dvxc_dmgd_rad = 0.0_DP
       return
    end if
    loc_prt = .false.
    if (present(prt)) loc_prt = prt

    ! Calculate gradient of density wrt r
    if (pub_xc_gradient_corrected) then
       ! Different approaches depending on whether this is a regular radial grid
       ! or a logarithmic one
       if (present(rmax).and.present(rab)) then
          call utils_abort('Error in xc_radial: Both rmax and rab supplied!')
       end if
       if (.not.(present(rmax).or.present(rab))) then
          print *,present(rmax),present(rab)
          call utils_abort('Error in xc_radial: No rmax or rab supplied!')
       end if
       if (.not.present(rad)) then
          call utils_abort('Error in xc_radial: No radial grid supplied!')
       end if

       do is=1,ns

          ! rmax supplied, hence regular grid
          if (present(rmax)) then
             call services_radial_derivative(grad_rad(:,is), &
                  den_rad(:,is),npts,rmax)
          end if
          ! rab supplied, hence logarithmic grid
          if (present(rab)) then
             call services_radial_derivative(grad_rad(:,is), &
                  den_rad(:,is),npts,real(npts,kind=DP))
             grad_rad(1:npts,is) = grad_rad(1:npts,is) / rab(1:npts)
          end if

       end do

    else
       grad_rad = 0.0_DP
    end if

    calc_x = .true.
    calc_c = .true.
    vxc_rad(:,:) = 0.0_DP
    if (present(exc_rad)) exc_rad(:,:) = 0.0_DP
    if (present(dvxc_dmgd_rad)) dvxc_dmgd_rad(:,:) = 0.0_DP

    select case (radial_functional)
    case ('CAPZ')
       !----------------------------------------------------------------------!
       ! Local density approximation - from Ceperley-Alder Monte Carlo data,  !
       !  and Gell-Mann-Brueckner expansion, parameterised by Perdew & Zunger !
       !----------------------------------------------------------------------!
       call internal_xc_radial_lda(xc_capz_c_point,xc_capz_c_point_sp)
    case ('VWN')
       !----------------------------------------------------------------------!
       ! Local density approximation - Vosko, Wilk & Nusair                   !
       !----------------------------------------------------------------------!
       call internal_xc_radial_lda(xc_vwn_c_point,xc_vwn_c_point_sp)
    case ('PW92')
       !----------------------------------------------------------------------!
       ! Local density approximation - Perdew-Wang 1992 functional            !
       !----------------------------------------------------------------------!
       call internal_xc_radial_lda(xc_pw92_c_point,xc_pw92_c_point_sp)
    case ('BLYP')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: Becke88 + LYP            !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_b88_x_point,xc_b88_x_point_sp, &
            xc_lyp_c_point, xc_lyp_c_point_sp)
    case ('PBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1996 version of          !
       ! Perdew, Burke & Ernzerhof                                            !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_pbe_x_point,xc_pbe_x_point_sp, &
            xc_pbe_c_point,xc_pbe_c_point_sp)
    ! la: Added VDWDF so that it is recognised whenever xc_radial is called
    case ('REVPBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1998 revised version of  !
       ! PBE due to Zhang & Yang                                              !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_revpbe_x_point,xc_revpbe_x_point_sp, &
            xc_pbe_c_point,xc_pbe_c_point_sp)
    case ('RPBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1999 revised version of  !
       ! PBE due to Hammer, Hansen & Norskov                                  !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_rpbe_x_point,&
            xc_rpbe_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp)
    case ('PBESOL')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 2008 revised version of  !
       ! PBE for solids due to Perdew et al.                                  !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_pbesol_x_point,&
            xc_pbesol_x_point_sp,xc_pbesol_c_point,xc_pbesol_c_point_sp)
    case ('AM05')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 2005 AM05                !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_am05_x_point,&
            xc_am05_x_point_sp,xc_am05_c_point,xc_am05_c_point_sp)
    case ('PW91')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1991 version of          !
       ! Perdew and Wang                                                      !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_pw91_x_point, &
            xc_pw91_x_point_sp,xc_pw91_c_point,xc_pw91_c_point_sp)
    case ('XLYP')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: X + LYP                  !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_x_x_point, &
            xc_x_x_point_sp, xc_lyp_c_point, xc_lyp_c_point_sp)
    case ('WC')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 2006 functional of Wu &  !
       ! Cohen                                                                !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_wc_x_point, &
            xc_wc_x_point_sp, xc_pbe_c_point, xc_pbe_c_point_sp)
    case ('B1LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional B1LYP:                                             !
       ! 1/4 HF_X + 3/4 B88_X - LDA_X) + LYP_C                                !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_b1_x_point, &
            xc_b1_x_point_sp, xc_lyp_c_point, xc_lyp_c_point_sp)
    case ('B1PW91')
       !----------------------------------------------------------------------!
       ! Hybrid functional B1PW91:                                             !
       ! 1/4 HF_X + 3/4 B88_X + PW91_C                                        !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_b1_x_point, &
            xc_b1_x_point_sp, xc_pw91_c_point, xc_pw91_c_point_sp)
    case ('B3LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional B3LYP:                                             !
       ! LDA + a(HF_X - LDA_X) + b (B88_X - LDA_X) + c(LYP_C - VWN_C)         !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_b3_x_point, &
            xc_b3_x_point_sp, xc_b3lyp_c_point, xc_b3lyp_c_point_sp)
    case ('B3PW91')
       !----------------------------------------------------------------------!
       ! Hybrid functional B3PW91:                                            !
       ! LDA + a(HF_X - LDA_X) + b (B88_X - LDA_X) + c(PW91_C - VWN_C)        !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_b3_x_point, &
            xc_b3_x_point_sp, xc_b3pw91_c_point, xc_b3pw91_c_point_sp)
    case ('PBE0')
       !----------------------------------------------------------------------!
       ! Hybrid functional PBE0: PBE + 1/4 (HF_X - PBE_X)                     !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_pbe0_x_point,&
            xc_pbe0_x_point_sp, xc_pbe_c_point, xc_pbe_c_point_sp)
    case ('X3LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional X3LYP:                                             !
       ! LDA + a(HF_X - LDA_X) + b ((d B88_X + e PW91_X - LDA_X)              !
       !     + c(LYP_C - VWN_C)                                               !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_x3_x_point, &
            xc_x3_x_point_sp, xc_x3lyp_c_point, xc_x3lyp_c_point_sp)
    case ('LIBXC')
       !----------------------------------------------------------------------!
       ! Functional from LIBXC library - type specified by pub_libxc_func_id  !
       !----------------------------------------------------------------------!
#ifdef LIBXC
       call internal_xc_radial_libxc
#else
       call utils_abort('Error in xc_radial: Functional "LIBXC" &
            &specified, but LIBXC not present in build')
#endif
    case ('NONE','HF')
       !----------------------------------------------------------------------!
       ! No exchange-correlation functional (Hartree approximation) or only   !
       ! Hartree-Fock exchange (Hartree-Fock approximation).                  !
       !----------------------------------------------------------------------!
       vxc_rad = 0.0_DP
       if (present(exc_rad)) exc_rad = 0.0_DP
    ! la: The vdW-DF functional
    case('VDWDF')
       !----------------------------------------------------------------------!
       ! vdW-DF functional                                                    !
       !----------------------------------------------------------------------!
       ! Calculate the LDA correlation potential and GGA exchange potential
       calc_c = .true.; calc_x = .false.
       call internal_xc_radial_lda(xc_pw92_c_point,xc_pw92_c_point_sp)
       calc_c = .false.; calc_x = .true.
       call internal_xc_radial_gc(xc_revpbe_x_point,xc_revpbe_x_point_sp, &
            xc_pbe_c_point,xc_pbe_c_point_sp)
    case('VDWDF2')
       !----------------------------------------------------------------------!
       ! vdW-DF2 functional                                                   !
       !----------------------------------------------------------------------!
       ! Calculate the LDA correlation potential and GGA exchange potential
       calc_c = .true.; calc_x = .false.
       call internal_xc_radial_lda(xc_pw92_c_point,xc_pw92_c_point_sp)
       calc_c = .false.; calc_x = .true.
       call internal_xc_radial_gc(xc_rpw86_x_point,xc_rpw86_x_point_sp, &
            xc_pbe_c_point,xc_pbe_c_point_sp)
   case('VDWC09')
       !----------------------------------------------------------------------!
       ! vdW-DF-C09x functional                                               !
       !----------------------------------------------------------------------!
       ! Calculate the LDA correlation potential and GGA exchange potential
       calc_c = .true.; calc_x = .false.
       call internal_xc_radial_lda(xc_pw92_c_point,xc_pw92_c_point_sp)
       calc_c = .false.; calc_x = .true.
       call internal_xc_radial_gc(xc_c09_x_point,xc_c09_x_point_sp, &
            xc_pbe_c_point,xc_pbe_c_point_sp)
   case('C09-2')
       !----------------------------------------------------------------------!
       ! vdW-DF-C09x functional with vdW-DF2 non-local correlation            !
       !----------------------------------------------------------------------!
       ! Calculate the LDA correlation potential and GGA exchange potential
       calc_c = .true.; calc_x = .false.
       call internal_xc_radial_lda(xc_pw92_c_point,xc_pw92_c_point_sp)
       calc_c = .false.; calc_x = .true.
       call internal_xc_radial_gc(xc_c09_x_point,xc_c09_x_point_sp, &
            xc_pbe_c_point,xc_pbe_c_point_sp)
    case('OPTPBE')
       !----------------------------------------------------------------------!
       ! vdW-DF-optPBE functional                                             !
       !----------------------------------------------------------------------!
       ! Calculate the LDA correlation potential and GGA exchange potential
       calc_c = .true.; calc_x = .false.
       call internal_xc_radial_lda(xc_pw92_c_point,xc_pw92_c_point_sp)
       calc_c = .false.; calc_x = .true.
       call internal_xc_radial_gc(xc_vdwoptpbe_x_point,xc_vdwoptpbe_x_point_sp, &
            xc_pbe_c_point,xc_pbe_c_point_sp)
    case('OPTB88')
       !----------------------------------------------------------------------!
       ! vdW-DF-optB88 functional                                             !
       !----------------------------------------------------------------------!
       ! Calculate the LDA correlation potential and GGA exchange potential
       calc_c = .true.; calc_x = .false.
       call internal_xc_radial_lda(xc_pw92_c_point,xc_pw92_c_point_sp)
       calc_c = .false.; calc_x = .true.
       call internal_xc_radial_gc(xc_vdwoptb88_x_point,xc_vdwoptb88_x_point_sp, &
            xc_pbe_c_point,xc_pbe_c_point_sp)
   case('PBEK1')
       !----------------------------------------------------------------------!
       ! vdW-DF-optPBEK=1 functional                                          !
       !----------------------------------------------------------------------!
       calc_c = .true.; calc_x = .false.
       call internal_xc_radial_lda(xc_pw92_c_point,xc_pw92_c_point_sp)
       calc_c = .false.; calc_x = .true.
       call internal_xc_radial_gc(xc_vdwpbek1_x_point,xc_vdwpbek1_x_point_sp, &
            xc_pbe_c_point,xc_pbe_c_point_sp)
    case('VV10')
       !----------------------------------------------------------------------!
       ! rVV10 functional, as implemented by Sabatini et al.                  !
       !----------------------------------------------------------------------!
       ! Calculate the GGA exchange/correlation potentials
       call internal_xc_radial_gc(xc_rpw86_x_point,xc_rpw86_x_point_sp, &
            xc_pbe_c_point,xc_pbe_c_point_sp)
    case('VV10BC')
       !----------------------------------------------------------------------!
       ! rVV10 functional, as implemented by Sabatini et al., but with b and  !
       ! C parameters from Vydrov and Van Voorhis original paper, i.e. b=5.9, !
       ! C=0.0093.                                                            !
       !----------------------------------------------------------------------!
       ! Calculate the GGA exchange/correlation potentials
       call internal_xc_radial_gc(xc_rpw86_x_point,xc_rpw86_x_point_sp, &
            xc_pbe_c_point,xc_pbe_c_point_sp)
    case('AVV10S')
       !----------------------------------------------------------------------!
       ! AM05-VV10sol functional                                              !
       !----------------------------------------------------------------------!
       ! Calculate the GGA exchange/correlation potentials
       call internal_xc_radial_gc(xc_am05_x_point,xc_am05_x_point_sp, &
            xc_am05_c_point,xc_am05_c_point_sp)
    case default
       !----------------------------------------------------------------------!
       ! Unrecognised functional type                                         !
       !----------------------------------------------------------------------!
       if (pub_on_root) write(stdout,'(2a)') 'WARNING: unknown exchange-&
            &correlation functional: ', trim(radial_functional)
    end select

    ! Now add term from dv/dmgd according to radial formula:
    ! v_xc[n(r),g(r)] = df_xc / dn - (d/dr(df/dg) + (2/r)(df/dg))
    if (pub_xc_gradient_corrected) then

       do is=1,ns

          do ipt=1,npts
             if (abs(grad_rad(ipt,is))>0.0_DP) then
                dvxc_dmgd_rad(ipt,is) = -dvxc_dmgd_rad(ipt,is) * &
                     grad_rad(ipt,is) / abs(grad_rad(ipt,is))
             else
                dvxc_dmgd_rad(ipt,is) = 0.0_DP
             end if
          end do
          if (loc_prt) then
          do ipt=1,npts
             write(81,'(i3,5f28.12)') &
                 is,rad(ipt),den_rad(ipt,is),grad_rad(ipt,is),vxc_rad(ipt,is),&
                 dvxc_dmgd_rad(ipt,is)
          end do
          end if

          ! rmax supplied, hence regular grid
          if (present(rmax)) then
             call services_radial_derivative(grad_rad(:,is), &
                  dvxc_dmgd_rad(:,is),npts,rmax)
          end if
          ! rab supplied, hence logarithmic grid
          if (present(rab)) then
             call services_radial_derivative(grad_rad(:,is), &
                  dvxc_dmgd_rad(:,is),npts,real(npts,kind=DP))
             grad_rad(1:npts,is) = grad_rad(1:npts,is) / rab(1:npts)
          end if
          if (loc_prt) then
          do ipt=1,npts
             write(91,'(i3,5f28.12)') &
                 is,rad(ipt),den_rad(ipt,is),grad_rad(ipt,is),vxc_rad(ipt,is),&
                 2.0_DP*dvxc_dmgd_rad(ipt,is)/rad(ipt)
          end do
          end if
          grad_rad(2:npts,is) = grad_rad(2:npts,is) + &
               2.0_DP*dvxc_dmgd_rad(2:npts,is)/rad(2:npts)
          vxc_rad(:npts,is) = vxc_rad(:npts,is) + grad_rad(:npts,is)

          if (loc_prt) write(81,*)
          if (loc_prt) write(91,*)

       end do

    end if

contains

    subroutine internal_xc_radial_lda(c_functional,c_functional_sp)

      interface
         subroutine c_functional(den,c_energy,c_pot)
           use constants, only: DP
           real(kind=DP), intent(in)  :: den
           real(kind=DP), intent(out) :: c_energy
           real(kind=DP), intent(out) :: c_pot
         end subroutine c_functional
      end interface

      interface
         subroutine c_functional_sp(den,den1,den2,c_energy,c_pot1,c_pot2)
           use constants, only: DP
           real(kind=DP), intent(in)  :: den
           real(kind=DP), intent(in)  :: den1
           real(kind=DP), intent(in)  :: den2
           real(kind=DP), intent(out) :: c_energy
           real(kind=DP), intent(out) :: c_pot1
           real(kind=DP), intent(out) :: c_pot2
         end subroutine c_functional_sp
      end interface

      if (ns==1) then
         do ipt=1,npts
            x_pot = 0.0_DP
            if (calc_x) then
               call xc_lda_x_point(den_rad(ipt,1),x_energy,x_pot(1))
            else
               x_energy = 0.0_DP; x_pot(1) = 0.0_DP
            end if
            if (calc_c) then
               call c_functional(den_rad(ipt,1),c_energy,c_pot(1))
            else
               c_energy = 0.0_DP; c_pot(1) = 0.0_DP
            end if

            vxc_rad(ipt,1) = vxc_rad(ipt,1) + x_pot(1) + c_pot(1)
            if (present(exc_rad)) &
                 exc_rad(ipt,1) = exc_rad(ipt,1) + x_energy + c_energy
         end do
      else
         do ipt=1,npts
            if (calc_x) then
               call xc_lsda_x_point(den_rad(ipt,1),den_rad(ipt,2),x_energy, &
                    x_pot(1),x_pot(2))
            else
               x_energy = 0.0_DP; x_pot(:) = 0.0_DP
            end if
            if (calc_c) then
               call c_functional_sp(sum(den_rad(ipt,:)),den_rad(ipt,1), &
                    den_rad(ipt,2),c_energy,c_pot(1),c_pot(2))
            else
               c_energy = 0.0_DP; c_pot(:) = 0.0_DP
            end if
            vxc_rad(ipt,:) = vxc_rad(ipt,:) + x_pot(:) + c_pot(:)
            if (present(exc_rad)) &
                 exc_rad(ipt,1) = exc_rad(ipt,1) + x_energy + c_energy
         end do
      end if

    end subroutine internal_xc_radial_lda

    subroutine internal_xc_radial_gc(x_functional,x_functional_sp, &
         c_functional,c_functional_sp)

      interface
         subroutine x_functional(den,mgd,x_energy,x_pot,x_dfdmgd)
           use constants, only: DP
           real(kind=DP), intent(in)  :: den
           real(kind=DP), intent(in)  :: mgd
           real(kind=DP), intent(out) :: x_energy
           real(kind=DP), intent(out) :: x_pot
           real(kind=DP), intent(out) :: x_dfdmgd
         end subroutine x_functional
      end interface

      interface
         subroutine x_functional_sp(den1,den2,mgd1,mgd2,x_energy,x_pot1,x_pot2,&
              x_dfdmgd1,x_dfdmgd2)
           use constants, only: DP
           real(kind=DP), intent(in)  :: den1
           real(kind=DP), intent(in)  :: den2
           real(kind=DP), intent(in)  :: mgd1
           real(kind=DP), intent(in)  :: mgd2
           real(kind=DP), intent(out) :: x_energy
           real(kind=DP), intent(out) :: x_pot1
           real(kind=DP), intent(out) :: x_pot2
           real(kind=DP), intent(out) :: x_dfdmgd1
           real(kind=DP), intent(out) :: x_dfdmgd2
         end subroutine x_functional_sp
      end interface

      interface
         subroutine c_functional(den,mgd,c_energy,c_pot,c_dfdmgd)
           use constants, only: DP
           real(kind=DP), intent(in)  :: den
           real(kind=DP), intent(in)  :: mgd
           real(kind=DP), intent(out) :: c_energy
           real(kind=DP), intent(out) :: c_pot
           real(kind=DP), intent(out) :: c_dfdmgd
         end subroutine c_functional
      end interface

      interface
         subroutine c_functional_sp(den,den1,den2,mgd,mgd1,mgd2,c_energy,&
              c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)
           use constants, only: DP
           real(kind=DP), intent(in)  :: den
           real(kind=DP), intent(in)  :: den1
           real(kind=DP), intent(in)  :: den2
           real(kind=DP), intent(in)  :: mgd
           real(kind=DP), intent(in)  :: mgd1
           real(kind=DP), intent(in)  :: mgd2
           real(kind=DP), intent(out) :: c_energy
           real(kind=DP), intent(out) :: c_pot1
           real(kind=DP), intent(out) :: c_pot2
           real(kind=DP), intent(out) :: c_dfdmgd
           real(kind=DP), intent(out) :: c_dfdmgd1
           real(kind=DP), intent(out) :: c_dfdmgd2
           real(kind=DP), intent(out) :: c_dfdmgd12
         end subroutine c_functional_sp
      end interface

      ! Local Variables
      real(kind=DP) :: mgd1, mgd2, mgd
      real(kind=DP) :: x_dfdmgd(2)
      real(kind=DP) :: c_dfdmgd(0:3)

      if (ns==1) then
         do ipt=1,npts
            mgd = abs(grad_rad(ipt,1))
            if (calc_x) then
               call x_functional(den_rad(ipt,1),mgd, x_energy, &
                    x_pot(1), x_dfdmgd(1))
            else
               x_energy = 0.0_DP; x_pot(1) = 0.0_DP; x_dfdmgd(1) = 0.0_DP
            end if
            if (calc_c) then
               call c_functional(den_rad(ipt,1),mgd, c_energy, &
                    c_pot(1), c_dfdmgd(1))
            else
               c_energy = 0.0_DP; c_pot(1) = 0.0_DP; c_dfdmgd(1) = 0.0_DP
            end if

            vxc_rad(ipt,1) = vxc_rad(ipt,1) + x_pot(1) + c_pot(1)

            if (present(dvxc_dmgd_rad)) then
               dvxc_dmgd_rad(ipt,1) = dvxc_dmgd_rad(ipt,1) &
                   + x_dfdmgd(1) + c_dfdmgd(1)
            end if
            if (present(exc_rad)) &
                 exc_rad(ipt,1) = exc_rad(ipt,1) + x_energy + c_energy
         end do
      else
         do ipt=1,npts
            mgd1 = abs(grad_rad(ipt,1))
            mgd2 = abs(grad_rad(ipt,2))
            mgd = abs(grad_rad(ipt,1)+grad_rad(ipt,2))
            if (calc_x) then
               call x_functional_sp(den_rad(ipt,1),den_rad(ipt,2),mgd1, &
                    mgd2,x_energy,x_pot(1),x_pot(2),x_dfdmgd(1), &
                    x_dfdmgd(2))
            else
               x_energy = 0.0_DP; x_pot(:) = 0.0_DP; x_dfdmgd(:) = 0.0_DP
            end if
            if (calc_c) then
               call c_functional_sp(sum(den_rad(ipt,:)),den_rad(ipt,1), &
                    den_rad(ipt,2),mgd,mgd1,mgd2,c_energy,c_pot(1),c_pot(2), &
                    c_dfdmgd(0),c_dfdmgd(1),c_dfdmgd(2),c_dfdmgd(3))
            else
               c_energy = 0.0_DP; c_pot(:) = 0.0_DP; c_dfdmgd(:) = 0.0_DP
            end if
            vxc_rad(ipt,:) = vxc_rad(ipt,:) + x_pot(:) + c_pot(:)
            if (present(dvxc_dmgd_rad)) then
               dvxc_dmgd_rad(ipt,1) = dvxc_dmgd_rad(ipt,1) &
                    + x_dfdmgd(1) + c_dfdmgd(0) + c_dfdmgd(1)
               dvxc_dmgd_rad(ipt,2) = dvxc_dmgd_rad(ipt,2) &
                    + x_dfdmgd(2) + c_dfdmgd(0) + c_dfdmgd(2)
               ! TODO: deal with c_dfdmgd(3)
               if (c_dfdmgd(3)>0.0_DP) call utils_abort('Error in xc_radial: &
                    &nonzero c_dfdmgd(3): Functional not yet supported by &
                    &xc_radial')
            end if
            if (present(exc_rad)) &
                 exc_rad(ipt,1) = exc_rad(ipt,1) + x_energy + c_energy
         end do
      end if

    end subroutine internal_xc_radial_gc

#ifdef LIBXC
    subroutine internal_xc_radial_libxc
      ! Updated for F2003 interface introduced in Libxc 3.0.0,
      ! James C. Womack, 2016
      use iso_c_binding, only: c_double
      use rundat, only: pub_libxc_x_func_id, pub_libxc_c_func_id, pub_num_spins
      use xc_f03_lib_m !! External dependency
      use utils, only: utils_assert

      ! Local Variables
      ! JCW: Arrays with length of 1 required to satisfy libxc F2003 interface
      real(kind=c_double) :: den(max_spins) ! Density
      real(kind=c_double) :: sigma(3) ! Modulus of density gradient
      real(kind=c_double) :: lapl_rho(max_spins) ! Laplacian of density
      real(kind=c_double) :: tau(max_spins)      ! Kinetic energy density
      real(kind=c_double) :: ex(1)               ! Accumulated xc energy
      real(kind=c_double) :: ec(1)               ! Accumulated xc energy
      real(kind=c_double) :: lx_pot(max_spins)   ! xc potential at point
      real(kind=c_double) :: lc_pot(max_spins)   ! xc potential at point
      real(kind=c_double) :: x_dfdsigma(3)       ! Derivative df_{xc}/dsigma
      real(kind=c_double) :: x_dfdlapl_rho(1)    ! Derivative df_{xc}/dlapl_rho
      real(kind=c_double) :: x_dfdtau(1)         ! Derivative df_{xc}/dtau
      real(kind=c_double) :: c_dfdsigma(3)       ! Derivative df_{xc}/dsigma
      real(kind=c_double) :: c_dfdlapl_rho(1)    ! Derivative df_{xc}/dlapl_rho
      real(kind=c_double) :: c_dfdtau(1)         ! Derivative df_{xc}/dtau
      real(kind=c_double), parameter :: two = real(2,kind=c_double)
      !
      integer :: libxc_polarized_rad
      type(xc_f03_func_t) :: libxc_func_rad(2)
      type(xc_f03_func_info_t) :: libxc_info_rad(2)

      if (ns.eq.1) then
         ! JCW: Create local versions of libxc_polarized, libxc_func, libxc_info
         ! JCW: for radial
         libxc_polarized_rad = XC_UNPOLARIZED
         if (pub_libxc_x_func_id > 0) then
            call xc_f03_func_init(libxc_func_rad(1), &
                 pub_libxc_x_func_id, libxc_polarized_rad)
            libxc_info_rad(1) = xc_f03_func_get_info(libxc_func_rad(1))
         end if
         if (pub_libxc_c_func_id > 0) then
            call xc_f03_func_init(libxc_func_rad(2), &
                 pub_libxc_c_func_id, libxc_polarized_rad)
            libxc_info_rad(2) = xc_f03_func_get_info(libxc_func_rad(2))
         end if
      else
         ! JCW: Assert that ns==1. The current interface is not set up to deal
         ! JCW: with the spin-polarized case correctly, and the atomic solver
         ! JCW: only ever calls xc_radial with ns==1 (at present).
         ! JCW: If ns>1, then c_dfdsigma(2) may be non-zero, and xc_radial does
         ! JCW: not currently support evaluating the gradient-corrected XC
         ! JCW: potential when this is the case.
         call utils_abort("Error in xc_radial: &
              &nspins > 1 not supported by radial Libxc interface.")
      end if

      ! JCW: Initialize local variables which are modified by Libxc to 0.0
      ! JCW: for safety. Fixes a bug where a check for c_dfdsigma(2)>0.0_DP
      ! JCW: was used to abort program, but setting this
      ! JCW: variable was left to Libxc (which does not necessarily zero the
      ! JCW: variable as expected).
      den(:)           = 0.0_c_double
      sigma(:)         = 0.0_c_double
      lapl_rho(:)      = 0.0_c_double
      tau(:)           = 0.0_c_double
      ex(1)            = 0.0_c_double
      ec(1)            = 0.0_c_double
      lx_pot(:)        = 0.0_c_double
      lc_pot(:)        = 0.0_c_double
      x_dfdsigma(:)    = 0.0_c_double
      x_dfdlapl_rho(1) = 0.0_c_double
      x_dfdtau(1)      = 0.0_c_double
      c_dfdsigma(:)    = 0.0_c_double
      c_dfdlapl_rho    = 0.0_c_double
      c_dfdtau(1)      = 0.0_c_double


      do ipt=1,npts

         den(1:ns) = real( den_rad(ipt,1:ns), kind=c_double)

         ! ndmh: Calculate energy and potential at this point
         if (.not.(libxc_gga.or.libxc_mgga)) then  ! Non-Gradient Corrected

            ! ndmh: Exchange part
            if (pub_libxc_x_func_id > 0) then
               call xc_f03_lda_exc_vxc(libxc_func_rad(1), 1, den(1), &
                    ex(1), lx_pot(1))
            else
               ex(1) = 0.0_c_double; lx_pot(:) = 0.0_c_double
            end if

            ! ndmh: Correlation part
            if (pub_libxc_c_func_id > 0) then
               call xc_f03_lda_exc_vxc(libxc_func_rad(2), 1, den(1), &
                    ec(1), lc_pot(1))
            else
               ec(1) = 0.0_c_double; lc_pot(:) = 0.0_c_double
            end if

         else if (libxc_gga.and.(.not.libxc_mgga)) then
            ! Gradient Corrected (GGA, not meta-GGA)

            ! Calculate reduced gradients
            if (ns==2) then
               sigma(1) = real( abs(grad_rad(ipt,1))**2, kind=c_double )
               sigma(2) = real( grad_rad(ipt,1)*grad_rad(ipt,2), kind=c_double )
               sigma(3) = real( abs(grad_rad(ipt,2))**2, kind=c_double )
            else
               sigma(1) = real( abs(grad_rad(ipt,1))**2, kind=c_double )
            end if

            ! ndmh: Exchange part
            if (pub_libxc_x_func_id > 0) then
               call xc_f03_gga_exc_vxc(libxc_func_rad(1), 1, den(1), sigma(1), &
                    ex(1), lx_pot(1), x_dfdsigma(1))
            else
               ex(1) = 0.0_c_double; lx_pot(:) = 0.0_c_double; x_dfdsigma(:) = 0.0_c_double
            end if

            ! ndmh: Correlation part
            if (pub_libxc_c_func_id > 0) then
               call xc_f03_gga_exc_vxc(libxc_func_rad(2), 1, den(1), sigma(1), &
                    ec(1), lc_pot(1), c_dfdsigma(1))
            else
               ec(1) = 0.0_c_double; lc_pot(:) = 0.0_c_double; c_dfdsigma(:) = 0.0_c_double
            end if

            if (present(dvxc_dmgd_rad)) then
               dvxc_dmgd_rad(ipt,1) = real( x_dfdsigma(1) + c_dfdsigma(1), kind=DP) &
                    * ( 2.0_DP * abs(grad_rad(ipt,1)) )
               if (ns==2) then
                 call utils_abort("Error in xc_radial: &
                      &nspins > 1 not supported by Libxc interface in &
                      &xc_radial")
               end if
               ! TODO: deal with c_dfdsigma(2)
               if (c_dfdsigma(2)>0.0_c_double) call utils_abort('Error in xc_radial: &
                    &nonzero c_dfdsigma(2): Functional not yet supported by &
                    &xc_radial')
            end if

         else if (libxc_mgga.and.(.not.libxc_gga)) then
            ! Gradient Corrected (meta-GGA, not GGA)
            call utils_abort(&
                 "Error in internal_xc_radial_libxc: meta-GGAs not implemented yet.")

            ! JCW: temporary debugging
            ! JCW: Set lapl_rho, tau to 0.0
            !lapl_rho = 0.0_c_double
            !tau = 0.0_c_double

            !! Calculate reduced gradients
            !if (ns==2) then
            !   sigma(1) = real( abs(grad_rad(ipt,1))**2, kind=c_double )
            !   sigma(2) = real( grad_rad(ipt,1)*grad_rad(ipt,2), kind=c_double )
            !   sigma(3) = real( abs(grad_rad(ipt,2))**2, kind=c_double )
            !else
            !   sigma(1) = real( abs(grad_rad(ipt,1))**2, kind=c_double )
            !end if

            !! JCW: Exchange functional
            !if (pub_libxc_x_func_id > 0) then
            !   call xc_f03_mgga_exc_vxc(libxc_func_rad(1), 1, den(1), &
            !        sigma(1), lapl_rho(1), tau(1), ex(1), lx_pot(1), &
            !        x_dfdsigma(1), x_dfdlapl_rho(1), x_dfdtau(1) )
            !else
            !   ex(1) = 0.0_c_double; lx_pot(:) = 0.0_c_double; x_dfdsigma(:) = 0.0_c_double
            !   x_dfdlapl_rho(1) = 0.0_c_double; x_dfdtau(1) = 0.0_c_double
            !end if

            !! JCW: Correlation functional
            !if (pub_libxc_c_func_id > 0) then
            !   call xc_f03_mgga_exc_vxc(libxc_func_rad(2), 1, den(1), &
            !        sigma(1), lapl_rho(1), tau(1), ec(1), lc_pot(1), &
            !        c_dfdsigma(1), c_dfdlapl_rho(1), c_dfdtau(1) )
            !else
            !   ec(1) = 0.0_c_double; lc_pot(:) = 0.0_c_double; c_dfdsigma(:) = 0.0_c_double
            !   x_dfdlapl_rho(1) = 0.0_c_double; x_dfdtau(1) = 0.0_c_double
            !end if

            !if (present(dvxc_dmgd_rad)) then
            !   ! Calculate d(v_xc)/d|grad(n)| using d(v_xc)/d(sigma) for
            !   ! singlet (S=0) case:
            !   ! d(v_xc)/d(sigma) = d(v_xc)/d|grad(n)| . d|grad(n)|/d(sigma)
            !   ! and
            !   ! d|grad(n)|/d(sigma) = 0.5 |grad(n)|
            !   ! thus
            !   ! d(v_xc)/d|grad(n)| = d(v_xc)/d(sigma) * 2.0 * |grad(n)|
            !   dvxc_dmgd_rad(ipt,1) = real(x_dfdsigma(1) + c_dfdsigma(1), kind=DP) &
            !        * ( 2.0_DP * abs(grad_rad(ipt,1)) )
            !   if (ns==2) then
            !     call utils_abort("Error in xc_radial: &
            !          &nspins > 1 not supported by Libxc interface in &
            !          &xc_radial")
            !   end if
            !   ! TODO: deal with c_dfdsigma(2)
            !   if (c_dfdsigma(2)>0.0_c_double) call utils_abort('Error in xc_radial: &
            !        &nonzero c_dfdsigma(2): Functional not yet supported by &
            !        &xc_radial')
            !end if
         else
            ! JCW: Neither or both libxc_gga and libxc_mgga true.
            ! JCW: Undefined state, so abort.
            call utils_abort(&
                 "Error in internal_xc_radial_libxc: unknown functional category.")
         end if

         vxc_rad(ipt,1) = vxc_rad(ipt,1) + real(lx_pot(1) + lc_pot(1),kind=DP)
         ! JCW: ex + ec multiplied by density because Libxc routines
         ! JCW: return XC energy per particle, rather than per unit
         ! JCW: volume. Here we assume ns == 1.
         if (present(exc_rad)) then
            exc_rad(ipt,1) = exc_rad(ipt,1) + real((ex(1) + ec(1))*den(1), kind=DP)
         end if

      end do

      ! JCW: de-initialize local libxc_func
      if (pub_libxc_x_func_id > 0) then
         call xc_f03_func_end(libxc_func_rad(1))
      end if
      if (pub_libxc_c_func_id > 0) then
         call xc_f03_func_end(libxc_func_rad(2))
      end if

    end subroutine internal_xc_radial_libxc
#endif

  end subroutine xc_radial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lda_triplet(density_fine,xc_pot_fine,grid,c_functional_sp)

    !==========================================================================!
    ! This subroutine calculates the exchange correlation potential in a spin  !
    ! polarised formulation, where rho(UP,DOWN)=(rho_0+epsilon*rho1,rho0).     !
    ! This is required when trying to calculate fxc in the finite difference   !
    ! approximation for triplet states in the spin-degenerate case. This       !
    ! routine is only needed for TDDFT triplet calculations in spin-degenerate !
    ! systems.                                                                 !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP
!$  use rundat, only: pub_threads_max
    use xc_funcs, only: xc_lda_x_point, xc_lsda_x_point

    implicit none

    ! Arguments
    ! Grid definition
    type(GRID_INFO), intent(in) :: grid

    ! Input density on the fine grid
    real(kind=DP), intent(in) :: density_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,2)

    ! Exchange-correlation potential
    real(kind=DP), intent(out) :: xc_pot_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,2)

    interface
       subroutine c_functional_sp(den,den1,den2,c_energy,c_pot1,c_pot2)
         use constants, only: DP
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(in)  :: den1
         real(kind=DP), intent(in)  :: den2
         real(kind=DP), intent(out) :: c_energy
         real(kind=DP), intent(out) :: c_pot1
         real(kind=DP), intent(out) :: c_pot2
       end subroutine c_functional_sp
    end interface

    ! Local variables:
    integer :: i1,i2,islab12   ! Fine grid cartesian indices
    integer :: ipt             ! Fine grid loop counter
    real(kind=DP) :: ex        ! Exchange energy in Ha
    real(kind=DP) :: ec        ! Correlation energy in Ha
    real(kind=DP) :: potint    ! Potential-density integral
    real(kind=DP) :: den       ! Charge density
    real(kind=DP) :: den1      ! Charge densities for spin 1
    real(kind=DP) :: den2      ! Charge densities for spin 2
    real(kind=DP) :: x_pot1    ! Exchange potential spin 1
    real(kind=DP) :: x_pot2    ! Exchange potential spin2
    real(kind=DP) :: x_energy  ! Exchange energy at point
    real(kind=DP) :: c_energy  ! Correlation energy at point
    real(kind=DP) :: c_pot1    ! Correlation potential spin 1
    real(kind=DP) :: c_pot2    ! Correlation potential spin 2

    ex = 0.0_DP ; ec = 0.0_DP
    potint = 0.0_DP

    xc_pot_fine = 0.0_DP ! jd: To clear margins between n1 and ld1 etc.

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE (ipt,i1,i2,islab12,den,den1,den2, &
!$OMP      x_energy,c_energy,x_pot1,x_pot2,c_pot1,c_pot2) &
!$OMP SHARED (grid,density_fine,xc_pot_fine,pub_threads_max)
    do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
       i1 = modulo(ipt-1,grid%n1) + 1
       i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
       islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1

       den1 = density_fine(i1,i2,islab12,1)
       den2 = density_fine(i1,i2,islab12,2)
       den = den1 + den2

       call xc_lsda_x_point(den1,den2,x_energy,x_pot1,x_pot2)
       call c_functional_sp(den,den1,den2,c_energy,c_pot1,c_pot2)

       xc_pot_fine(i1,i2,islab12,1) = x_pot1 + c_pot1
       xc_pot_fine(i1,i2,islab12,2) = x_pot2 + c_pot2

    end do
!$OMP END PARALLEL DO

  end subroutine xc_lda_triplet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lda(density_fine,xc_energy,xc_pot_fine,grid,c_functional,&
       c_functional_sp,xc_potint,exchange,correlation,x_pot_fine,c_pot_fine,&
       parallel)

    !==========================================================================!
    ! This subroutine calculates the exchange-correlation energy and potential !
    ! according to the local density approximation (LDA), given the electronic !
    ! density on the fine grid. The correlation functional to be used is given !
    ! as an argument.                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in March 2009 based on existing xc_capz          !
    ! subroutine written by Peter Haynes in February 2004.                     !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP, EDA_ISOLATED, EDA_FROZENIDEM, EDA_FROZENNONIDEM
    use rundat, only: pub_num_spins, pub_eda, pub_eda_mode, pub_frag_counter, &
         pub_eda_isol_x, pub_eda_isol_c, &
         pub_eda_frz_x, pub_eda_rep_x, &
         pub_eda_frz_c, pub_eda_rep_c, pub_eda_read_super
!$  use rundat, only: pub_threads_max
    use xc_funcs, only: xc_lda_x_point, xc_lsda_x_point

    implicit none

    ! Arguments
    ! Grid definition
    type(GRID_INFO), intent(in) :: grid

    ! Input density on the fine grid
    real(kind=DP), intent(in) :: density_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)

    ! Exchange-correlation energy
    real(kind=DP), intent(out) :: xc_energy

    ! Exchange-correlation potential
    real(kind=DP), intent(out) :: xc_pot_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)

    ! Exchange-correlation potential-density integral
    real(kind=DP), optional, intent(out) :: exchange
    real(kind=DP), optional, intent(out) :: correlation
    real(kind=DP), optional, intent(out) :: x_pot_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(out) :: c_pot_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(out) :: xc_potint
    logical,       optional, intent(in)  :: parallel

    interface
       subroutine c_functional(den,c_energy,c_pot)
         use constants, only: DP
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(out) :: c_energy
         real(kind=DP), intent(out) :: c_pot
       end subroutine c_functional
    end interface

    interface
       subroutine c_functional_sp(den,den1,den2,c_energy,c_pot1,c_pot2)
         use constants, only: DP
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(in)  :: den1
         real(kind=DP), intent(in)  :: den2
         real(kind=DP), intent(out) :: c_energy
         real(kind=DP), intent(out) :: c_pot1
         real(kind=DP), intent(out) :: c_pot2
       end subroutine c_functional_sp
    end interface

    ! Local variables:
    integer :: i1,i2,islab12   ! Fine grid cartesian indices
    integer :: ipt             ! Fine grid loop counter
    real(kind=DP) :: ex        ! Exchange energy in Ha
    real(kind=DP) :: ec        ! Correlation energy in Ha
    real(kind=DP) :: potint    ! Potential-density integral
    real(kind=DP) :: den       ! Charge density
    real(kind=DP) :: den1      ! Charge densities for spin 1
    real(kind=DP) :: den2      ! Charge densities for spin 2
    real(kind=DP) :: x_pot     ! Exchange potential
    real(kind=DP) :: x_pot1    ! Exchange potential spin 1
    real(kind=DP) :: x_pot2    ! Exchange potential spin2
    real(kind=DP) :: x_energy  ! Exchange energy at point
    real(kind=DP) :: c_energy  ! Correlation energy at point
    real(kind=DP) :: c_pot     ! Correlation potential
    real(kind=DP) :: c_pot1    ! Correlation potential spin 1
    real(kind=DP) :: c_pot2    ! Correlation potential spin 2
    logical       :: loc_parallel

    loc_parallel = .true.
    if(present(parallel)) then
       loc_parallel = parallel
    end if

    ex = 0.0_DP ; ec = 0.0_DP
    potint = 0.0_DP

    xc_pot_fine = 0.0_DP ! jd: To clear margins between n1 and ld1 etc.

    if (pub_num_spins == 1) then

       ! qoh: Spin-unpolarised case:
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE (ipt,i1,i2,islab12,den, &
!$OMP      x_energy,c_energy,x_pot,c_pot) &
!$OMP SHARED (grid,density_fine,xc_pot_fine,x_pot_fine, &
!$OMP      c_pot_fine,pub_threads_max) &
!$OMP REDUCTION(+:ex,ec,potint)
      do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
          i1 = modulo(ipt-1,grid%n1) + 1
          i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
          islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1

          den = density_fine(i1,i2,islab12,1)

          call xc_lda_x_point(den,x_energy,x_pot)
          call c_functional(den,c_energy,c_pot)

          ex = ex + x_energy
          ec = ec + c_energy

          xc_pot_fine(i1,i2,islab12,1) = x_pot + c_pot
          if(present(x_pot_fine)) then
              x_pot_fine(i1,i2,islab12,1) = x_pot
          endif
          if(present(c_pot_fine)) then
              c_pot_fine(i1,i2,islab12,1) = c_pot
          endif
          potint = potint + den * xc_pot_fine(i1,i2,islab12,1)

       end do
!$OMP END PARALLEL DO

    else

       ! qoh: Spin-polarised case:
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE (ipt,i1,i2,islab12,den,den1,den2, &
!$OMP      x_energy,c_energy,x_pot1,x_pot2,c_pot1,c_pot2) &
!$OMP SHARED (grid,density_fine,xc_pot_fine,x_pot_fine, &
!$OMP      c_pot_fine,pub_threads_max) &
!$OMP REDUCTION(+:ex,ec,potint)
       do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
          i1 = modulo(ipt-1,grid%n1) + 1
          i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
          islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1

          den1 = density_fine(i1,i2,islab12,1)
          den2 = density_fine(i1,i2,islab12,2)
          den = den1 + den2

          call xc_lsda_x_point(den1,den2,x_energy,x_pot1,x_pot2)
          call c_functional_sp(den,den1,den2,c_energy,c_pot1,c_pot2)


           if(present(x_pot_fine)) then
               x_pot_fine(i1,i2,islab12,1) = x_pot1
               x_pot_fine(i1,i2,islab12,2) = x_pot2
           endif
           if(present(c_pot_fine)) then
               c_pot_fine(i1,i2,islab12,1) = c_pot1
               c_pot_fine(i1,i2,islab12,2) = c_pot2
           endif

          ex = ex + x_energy
          ec = ec + c_energy

          xc_pot_fine(i1,i2,islab12,1) = x_pot1 + c_pot1
          xc_pot_fine(i1,i2,islab12,2) = x_pot2 + c_pot2
          potint = potint + den1 * xc_pot_fine(i1,i2,islab12,1)
          potint = potint + den2 * xc_pot_fine(i1,i2,islab12,2)

       end do
!$OMP END PARALLEL DO

    end if

    ! Global reductions
    if(loc_parallel) then
       call comms_reduce('SUM',ex)
       call comms_reduce('SUM',ec)
       call comms_reduce('SUM',potint)
    end if

    if(present(exchange)) then
    exchange  = ex * grid%weight
    endif
    if(present(correlation)) then
    correlation  = ec * grid%weight
    endif

    ! mjsp: if this is an isolated fragment calculation then
    ! update the fragment's stored X,C energies
    if (pub_eda) then
       if (pub_eda_mode == EDA_ISOLATED) then
          pub_eda_isol_x(pub_frag_counter) = (ex)*grid%weight
          pub_eda_isol_c(pub_frag_counter) = (ec)*grid%weight
       ! mjsp: if this is the frozen energy calculation
       ! then store the frozen state's X,C contributions
       ! (if pub_eda_read_super then this is contained in the .eda file)
       else if (pub_eda_mode == EDA_FROZENNONIDEM .and. .not.(pub_eda_read_super)) then
          pub_eda_frz_x = (ex)*grid%weight
          pub_eda_frz_c = (ec)*grid%weight
       ! mjsp: if this is the repulsion energy calculation
       ! then store the Pauli repulsion state's X,C contributions
       else if (pub_eda_mode == EDA_FROZENIDEM) then
          pub_eda_rep_x = (ex)*grid%weight
          pub_eda_rep_c = (ec)*grid%weight
       end if
    end if

    xc_energy = (ex + ec) * grid%weight

    if (present(xc_potint)) xc_potint = potint * grid%weight

  end subroutine xc_lda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine xc_gc_triplet(density_fine,density_grad,density_aux,&
       recip_work,grid,xc_pot_fine,x_functional_sp,c_functional_sp)

    !==========================================================================!
    ! This subroutine calculates the gradient corrected exchange-correlation   !
    ! potential in a spin                                                      !
    ! polarised formulation, where rho(UP,DOWN)=(rho_0+epsilon*rho1,rho0).     !
    ! This is required when trying to calculate fxc in the finite difference   !
    ! approximation for triplet states in the spin-degenerate case. This       !
    ! routine is only needed for TDDFT triplet calculations in spin-degenerate !
    ! systems.                                                                 !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP
!$  use rundat, only: pub_threads_max

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(in) :: density_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,2)
    real(kind=DP), intent(out) :: xc_pot_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,2)
    real(kind=DP),intent(inout) :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,2)
    real(kind=DP),intent(inout) :: density_aux(grid%ld1,grid%ld2, &
         grid%max_slabs12,2)
    complex(kind=DP),intent(inout) :: recip_work(grid%ld3,grid%ld2, &
         grid%max_slabs23,3)

    interface
       subroutine x_functional_sp(den1,den2,mgd1,mgd2,x_energy,x_pot1,x_pot2,&
            x_dfdmgd1,x_dfdmgd2)
         use constants, only: DP
         real(kind=DP), intent(in)  :: den1
         real(kind=DP), intent(in)  :: den2
         real(kind=DP), intent(in)  :: mgd1
         real(kind=DP), intent(in)  :: mgd2
         real(kind=DP), intent(out) :: x_energy
         real(kind=DP), intent(out) :: x_pot1
         real(kind=DP), intent(out) :: x_pot2
         real(kind=DP), intent(out) :: x_dfdmgd1
         real(kind=DP), intent(out) :: x_dfdmgd2
       end subroutine x_functional_sp
    end interface

    interface
       subroutine c_functional_sp(den,den1,den2,mgd,mgd1,mgd2,c_energy,&
            c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)
         use constants, only: DP
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(in)  :: den1
         real(kind=DP), intent(in)  :: den2
         real(kind=DP), intent(in)  :: mgd
         real(kind=DP), intent(in)  :: mgd1
         real(kind=DP), intent(in)  :: mgd2
         real(kind=DP), intent(out) :: c_energy
         real(kind=DP), intent(out) :: c_pot1
         real(kind=DP), intent(out) :: c_pot2
         real(kind=DP), intent(out) :: c_dfdmgd
         real(kind=DP), intent(out) :: c_dfdmgd1
         real(kind=DP), intent(out) :: c_dfdmgd2
         real(kind=DP), intent(out) :: c_dfdmgd12
       end subroutine c_functional_sp
    end interface

    real(kind=DP) :: den      ! Density
    real(kind=DP) :: den1     ! Density spin 1
    real(kind=DP) :: den2     ! Density spin 2
    real(kind=DP) :: mgd      ! Modulus of density gradient
    real(kind=DP) :: mgd1     ! Modulus of density gradient for spin 1
    real(kind=DP) :: mgd2     ! Modulus of density gradient for spin 2
    real(kind=DP) :: grad(3)  ! Density gradient
    real(kind=DP) :: grad1(3) ! Density gradient for spin 1
    real(kind=DP) :: grad2(3) ! Density gradient for spin 2
    !qoh: Exchange variables
    real(kind=DP) :: x_energy  ! Exchange energy at point
    real(kind=DP) :: x_pot1    ! Exchange potential at point for spin 1
    real(kind=DP) :: x_pot2    ! Exchange potential at point for spin 2
    real(kind=DP) :: x_dfdmgd! Derivative df_{x}/d|grad n|
    real(kind=DP) :: x_dfdmgd1! Derivative df_{x}/d|grad n| spin 1
    real(kind=DP) :: x_dfdmgd2! Derivative df_{x}/d|grad n| spin 2
    !qoh: Correlation variables
    real(kind=DP) :: c_energy  ! Correlation energy at point
    real(kind=DP) :: c_pot1    ! Correlation potential at point for spin 1
    real(kind=DP) :: c_pot2    ! Correlation potential at point for spin 2
    real(kind=DP) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=DP) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n| spin 1
    real(kind=DP) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n| spin 2
    real(kind=DP) :: c_dfdmgd12 ! df_{c}/d(grad n_1.grad n_2) component
    real(kind=DP) :: gd1_dot_gd2 ! (grad n_1).(grad n_2)
    integer :: i1,i2,islab12   ! Fine grid cartesian indices
    integer :: ipt             ! Fine grid loop counter

    xc_pot_fine = 0.0_DP ! jd: To clear margins between n1 and ld1 etc.

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE (ipt,i1,i2,islab12,den,den1,den2,grad,grad1,grad2, &
!$OMP      mgd,mgd1,mgd2,gd1_dot_gd2,x_energy,c_energy,x_pot1,x_pot2, &
!$OMP      c_pot1,c_pot2,x_dfdmgd,x_dfdmgd1,x_dfdmgd2,c_dfdmgd,c_dfdmgd1, &
!$OMP      c_dfdmgd2,c_dfdmgd12) &
!$OMP SHARED (grid,density_fine,density_grad,xc_pot_fine,pub_threads_max)
    do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
       i1 = modulo(ipt-1,grid%n1) + 1
       i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
       islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1

       den1 = density_fine(i1,i2,islab12,1)
       den2 = density_fine(i1,i2,islab12,2)
       den = den1 + den2

       grad1 = density_grad(i1,i2,islab12,:,1)
       grad2 = density_grad(i1,i2,islab12,:,2)
       grad = grad1 + grad2

       mgd1 = sqrt(grad1(1)*grad1(1) + grad1(2)*grad1(2) + &
            grad1(3)*grad1(3))
       mgd2 = sqrt(grad2(1)*grad2(1) + grad2(2)*grad2(2) + &
            grad2(3)*grad2(3))
       mgd = sqrt(grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3))
            gd1_dot_gd2 = grad1(1)*grad2(1) + grad1(2)*grad2(2) + &
            grad1(3)*grad2(3)
       ! qoh: Calculate energy and potential at this point
       call x_functional_sp(den1,den2,mgd1,mgd2,&
            x_energy, x_pot1, x_pot2, x_dfdmgd1, x_dfdmgd2)
       call c_functional_sp(den,den1,den2,mgd,mgd1,mgd2,&
            c_energy, c_pot1, c_pot2, &
            c_dfdmgd, c_dfdmgd1, c_dfdmgd2, c_dfdmgd12)

       ! qoh: Store this point's energy and potential
       xc_pot_fine(i1,i2,islab12,1) = x_pot1 + c_pot1
       xc_pot_fine(i1,i2,islab12,2) = x_pot2 + c_pot2

       ! qoh: Make sure we don't divide by zero
       if (mgd1 > 0.0_DP) then
          density_grad(i1,i2,islab12,:,1) = &
               -(x_dfdmgd1+c_dfdmgd1) * grad1(:) / mgd1
       else
          density_grad(i1,i2,islab12,:,1) = 0.0_DP
       end if
       if (mgd2 > 0.0_DP) then
          density_grad(i1,i2,islab12,:,2) = &
               -(x_dfdmgd2+c_dfdmgd2) * grad2(:) / mgd2

       else
          density_grad(i1,i2,islab12,:,2) = 0.0_DP
       end if

       if (abs(gd1_dot_gd2)>0.0_DP) then
          density_grad(i1,i2,islab12,:,1) = &
               density_grad(i1,i2,islab12,:,1) &
               -c_dfdmgd12 * grad2(:) / gd1_dot_gd2
          density_grad(i1,i2,islab12,:,2) = &
               density_grad(i1,i2,islab12,:,2) &
               -c_dfdmgd12 * grad1(:) / gd1_dot_gd2
       end if
       if (mgd > 0.0_DP) then
          density_grad(i1,i2,islab12,:,1) = &
               density_grad(i1,i2,islab12,:,1) &
               -c_dfdmgd * grad(:) / mgd
          density_grad(i1,i2,islab12,:,2) = &
               density_grad(i1,i2,islab12,:,2) &
               -c_dfdmgd * grad(:) / mgd
       end if
    end do
!$OMP END PARALLEL DO

    ! Calculate divergence of gradient dependent part and add to potential
    call xc_divergence(density_grad(:,:,:,:,1),density_aux(:,:,:,1),&
          recip_work,grid)
    call xc_divergence(density_grad(:,:,:,:,2),density_aux(:,:,:,2),&
          recip_work,grid)
    xc_pot_fine = xc_pot_fine + density_aux

  end subroutine xc_gc_triplet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
       xc_energy,xc_pot_fine,x_functional,x_functional_sp, &
       c_functional,c_functional_sp,xc_potint,exchange,correlation,parallel)

    !==========================================================================!
    ! This subroutine calculates the exchange-correlation energy and potential !
    ! for a gradient corrected functional (including hybrids), given the       !
    ! electronic density on the fine grid.  The exchange and correlation       !
    ! functionals to be used are given as arguments.                           !
    !--------------------------------------------------------------------------!
    ! Arguments                                                                !
    !   grid (in) : Grid definition                                            !
    !   density_fine (in) : Input density on the fine grid                     !
    !   xc_energy (out) : Exchange-correlation energy                          !
    !   xc_pot_fine (out) : Exchange-correlation potential                     !
    !   xc_potint (out) : Exchange-correlation potential-density integral      !
    !   x_functional,x_functional_sp (in) : X subroutines given as arguements  !
    !   c_functional,c_functional_sp (in) : C subroutines given as arguements  !
    !   density_grad (inout) : Density gradient on the fine grid               !
    !   density_aux (inout) : Auxiliary density on the fine grid               !
    !   recip_work (inout) : Fine Grid Recip-space workspace                   !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in February / March 2009 based on existing       !
    ! xc_pbe subroutine written by Peter Haynes in February 2004.              !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce, pub_on_root
    use constants, only: DP, EDA_ISOLATED, EDA_FROZENIDEM, EDA_FROZENNONIDEM
    use rundat, only: pub_num_spins, pub_eda, pub_eda_mode, pub_frag_counter, &
         pub_eda_isol_x, pub_eda_isol_c, &
         pub_eda_frz_x, pub_eda_rep_x, &
         pub_eda_frz_c, pub_eda_rep_c, pub_eda_read_super
!$  use rundat, only: pub_threads_max

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(in) :: density_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(out) :: xc_energy
    real(kind=DP), intent(out) :: xc_pot_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(out) :: xc_potint
    real(kind=DP), optional, intent(out) :: exchange
    real(kind=DP), optional, intent(out) :: correlation
    real(kind=DP),intent(inout) :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_num_spins)
    real(kind=DP),intent(inout) :: density_aux(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    complex(kind=DP),intent(inout) :: recip_work(grid%ld3,grid%ld2, &
         grid%max_slabs23,3)
    logical,       optional, intent(in) :: parallel

    interface
       subroutine x_functional(den,mgd,x_energy,x_pot,x_dfdmgd)
         use constants, only: DP
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(in)  :: mgd
         real(kind=DP), intent(out) :: x_energy
         real(kind=DP), intent(out) :: x_pot
         real(kind=DP), intent(out) :: x_dfdmgd
       end subroutine x_functional
    end interface

    interface
       subroutine x_functional_sp(den1,den2,mgd1,mgd2,x_energy,x_pot1,x_pot2,&
            x_dfdmgd1,x_dfdmgd2)
         use constants, only: DP
         real(kind=DP), intent(in)  :: den1
         real(kind=DP), intent(in)  :: den2
         real(kind=DP), intent(in)  :: mgd1
         real(kind=DP), intent(in)  :: mgd2
         real(kind=DP), intent(out) :: x_energy
         real(kind=DP), intent(out) :: x_pot1
         real(kind=DP), intent(out) :: x_pot2
         real(kind=DP), intent(out) :: x_dfdmgd1
         real(kind=DP), intent(out) :: x_dfdmgd2
       end subroutine x_functional_sp
    end interface

    interface
       subroutine c_functional(den,mgd,c_energy,c_pot,c_dfdmgd)
         use constants, only: DP
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(in)  :: mgd
         real(kind=DP), intent(out) :: c_energy
         real(kind=DP), intent(out) :: c_pot
         real(kind=DP), intent(out) :: c_dfdmgd
       end subroutine c_functional
    end interface

    interface
       subroutine c_functional_sp(den,den1,den2,mgd,mgd1,mgd2,c_energy,&
            c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)
         use constants, only: DP
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(in)  :: den1
         real(kind=DP), intent(in)  :: den2
         real(kind=DP), intent(in)  :: mgd
         real(kind=DP), intent(in)  :: mgd1
         real(kind=DP), intent(in)  :: mgd2
         real(kind=DP), intent(out) :: c_energy
         real(kind=DP), intent(out) :: c_pot1
         real(kind=DP), intent(out) :: c_pot2
         real(kind=DP), intent(out) :: c_dfdmgd
         real(kind=DP), intent(out) :: c_dfdmgd1
         real(kind=DP), intent(out) :: c_dfdmgd2
         real(kind=DP), intent(out) :: c_dfdmgd12
       end subroutine c_functional_sp
    end interface

    real(kind=DP) :: den      ! Density
    real(kind=DP) :: den1     ! Density spin 1
    real(kind=DP) :: den2     ! Density spin 2
    real(kind=DP) :: mgd      ! Modulus of density gradient
    real(kind=DP) :: mgd1     ! Modulus of density gradient for spin 1
    real(kind=DP) :: mgd2     ! Modulus of density gradient for spin 2
    real(kind=DP) :: grad(3)  ! Density gradient
    real(kind=DP) :: grad1(3) ! Density gradient for spin 1
    real(kind=DP) :: grad2(3) ! Density gradient for spin 2
    real(kind=DP) :: potint   ! Potential integral
    !qoh: Exchange variables
    real(kind=DP) :: ex        ! Accumulated exchange energy
    real(kind=DP) :: x_energy  ! Exchange energy at point
    real(kind=DP) :: x_pot     ! Exchange potential at point
    real(kind=DP) :: x_pot1    ! Exchange potential at point for spin 1
    real(kind=DP) :: x_pot2    ! Exchange potential at point for spin 2
    real(kind=DP) :: x_dfdmgd! Derivative df_{x}/d|grad n|
    real(kind=DP) :: x_dfdmgd1! Derivative df_{x}/d|grad n| spin 1
    real(kind=DP) :: x_dfdmgd2! Derivative df_{x}/d|grad n| spin 2
    !qoh: Correlation variables
    real(kind=DP) :: ec        ! Accumulated correlation energy
    real(kind=DP) :: c_energy  ! Correlation energy at point
    real(kind=DP) :: c_pot     ! Correlation potential at point
    real(kind=DP) :: c_pot1    ! Correlation potential at point for spin 1
    real(kind=DP) :: c_pot2    ! Correlation potential at point for spin 2
    real(kind=DP) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=DP) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n| spin 1
    real(kind=DP) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n| spin 2
    real(kind=DP) :: c_dfdmgd12 ! df_{c}/d(grad n_1.grad n_2) component
    real(kind=DP) :: gd1_dot_gd2 ! (grad n_1).(grad n_2)
    integer :: i1,i2,islab12   ! Fine grid cartesian indices
    integer :: ipt             ! Fine grid loop counter
    integer :: is ! Spin counter
    logical :: loc_parallel


    real(kind=dp) :: dummy

    loc_parallel = .true.
    if(present(parallel)) then
       loc_parallel = parallel
    end if

    ex = 0.0_DP
    ec = 0.0_DP

    xc_pot_fine = 0.0_DP ! jd: To clear margins between n1 and ld1 etc.

    spin: if ( pub_num_spins == 1 ) then
       !qoh: Spin unpolarised case
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE (ipt,i1,i2,islab12,den,grad,mgd,x_energy,c_energy, &
!$OMP      x_pot,c_pot,x_dfdmgd,c_dfdmgd) &
!$OMP SHARED (grid,density_fine,density_grad,xc_pot_fine, &
!$OMP      pub_threads_max) &
!$OMP REDUCTION(+:ex,ec)
       do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
          i1 = modulo(ipt-1,grid%n1) + 1
          i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
          islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1

          den = density_fine(i1,i2,islab12,1)
          grad = density_grad(i1,i2,islab12,:,1)
          mgd = sqrt(grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3))

          ! qoh: Calculate energy and potential at this point
          call x_functional(den,mgd,x_energy, x_pot, x_dfdmgd)
          call c_functional(den,mgd,c_energy, c_pot, c_dfdmgd)

          ! qoh: Store this point's energy and potential
          ex = x_energy + ex
          ec = c_energy + ec
          xc_pot_fine(i1,i2,islab12,1) = x_pot + c_pot
          ! qoh: Make sure we don't divide by zero
          if (mgd > 0.0_DP) then
             density_grad(i1,i2,islab12,:,1) = &
                  -(c_dfdmgd + x_dfdmgd) * grad(:) / mgd
          else
             density_grad(i1,i2,islab12,:,1) = 0.0_DP
          end if

       end do
!$OMP END PARALLEL DO

    else
       !qoh: Spin polarised case
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE (ipt,i1,i2,islab12,den,den1,den2,grad,grad1,grad2, &
!$OMP      mgd,mgd1,mgd2,gd1_dot_gd2,x_energy,c_energy,x_pot1,x_pot2, &
!$OMP      c_pot1,c_pot2,x_dfdmgd,x_dfdmgd1,x_dfdmgd2,c_dfdmgd,c_dfdmgd1, &
!$OMP      c_dfdmgd2,c_dfdmgd12) &
!$OMP SHARED (grid,density_fine,density_grad,xc_pot_fine, &
!$OMP      pub_threads_max) &
!$OMP REDUCTION(+:ex,ec)
       do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
          i1 = modulo(ipt-1,grid%n1) + 1
          i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
          islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1

          den1 = density_fine(i1,i2,islab12,1)
          den2 = density_fine(i1,i2,islab12,2)
          den = den1 + den2

          grad1 = density_grad(i1,i2,islab12,:,1)
          grad2 = density_grad(i1,i2,islab12,:,2)
          grad = grad1 + grad2

          mgd1 = sqrt(grad1(1)*grad1(1) + grad1(2)*grad1(2) + &
               grad1(3)*grad1(3))
          mgd2 = sqrt(grad2(1)*grad2(1) + grad2(2)*grad2(2) + &
               grad2(3)*grad2(3))
          mgd = sqrt(grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3))
          gd1_dot_gd2 = grad1(1)*grad2(1) + grad1(2)*grad2(2) + &
               grad1(3)*grad2(3)
          ! qoh: Calculate energy and potential at this point
          call x_functional_sp(den1,den2,mgd1,mgd2,&
               x_energy, x_pot1, x_pot2, x_dfdmgd1, x_dfdmgd2)
          call c_functional_sp(den,den1,den2,mgd,mgd1,mgd2,&
               c_energy, c_pot1, c_pot2, &
               c_dfdmgd, c_dfdmgd1, c_dfdmgd2, c_dfdmgd12)

          ! qoh: Store this point's energy and potential
          ex = x_energy + ex
          ec = c_energy + ec
          xc_pot_fine(i1,i2,islab12,1) = x_pot1 + c_pot1
          xc_pot_fine(i1,i2,islab12,2) = x_pot2 + c_pot2

          ! qoh: Make sure we don't divide by zero
          if (mgd1 > 0.0_DP) then
             density_grad(i1,i2,islab12,:,1) = &
                  -(x_dfdmgd1+c_dfdmgd1) * grad1(:) / mgd1
          else
             density_grad(i1,i2,islab12,:,1) = 0.0_DP
          end if
          if (mgd2 > 0.0_DP) then
             density_grad(i1,i2,islab12,:,2) = &
                  -(x_dfdmgd2+c_dfdmgd2) * grad2(:) / mgd2

          else
             density_grad(i1,i2,islab12,:,2) = 0.0_DP

          end if
          if (abs(gd1_dot_gd2)>0.0_DP) then
             density_grad(i1,i2,islab12,:,1) = &
                  density_grad(i1,i2,islab12,:,1) &
                  -c_dfdmgd12 * grad2(:) / gd1_dot_gd2
             density_grad(i1,i2,islab12,:,2) = &
                  density_grad(i1,i2,islab12,:,2) &
                  -c_dfdmgd12 * grad1(:) / gd1_dot_gd2
          end if
          if (mgd > 0.0_DP) then
             density_grad(i1,i2,islab12,:,1) = &
                  density_grad(i1,i2,islab12,:,1) &
                  -c_dfdmgd * grad(:) / mgd
             density_grad(i1,i2,islab12,:,2) = &
                  density_grad(i1,i2,islab12,:,2) &
                  -c_dfdmgd * grad(:) / mgd
          end if
       end do
!$OMP END PARALLEL DO
    end if spin


    ! Calculate divergence of gradient dependent part and add to potential
    call xc_divergence(density_grad,density_aux,recip_work,grid)
    xc_pot_fine = xc_pot_fine + density_aux

    ! Global reductions
    if(loc_parallel) then
       call comms_reduce('SUM',ex)
       call comms_reduce('SUM',ec)
    end if

    if(present(exchange)) then
    exchange  = ex * grid%weight
    endif
    if(present(correlation)) then
    correlation  = ec * grid%weight
    endif

    ! mjsp: if this is an isolated fragment calculation then
    ! update the fragment's stored X,C energies
    if (pub_eda) then
       if (pub_eda_mode == EDA_ISOLATED) then
          pub_eda_isol_x(pub_frag_counter) = (ex)*grid%weight
          pub_eda_isol_c(pub_frag_counter) = (ec)*grid%weight
       ! mjsp: if this is the frozen energy calculation
       ! then store the frozen state's X,C contributions
       ! (if pub_eda_read_super then this is contained in the .eda file)
       else if (pub_eda_mode == EDA_FROZENNONIDEM .and. .not.(pub_eda_read_super)) then
          pub_eda_frz_x = (ex)*grid%weight
          pub_eda_frz_c = (ec)*grid%weight
       ! mjsp: if this is the repulsion energy calculation
       ! then store the Pauli repulsion state's X,C contributions
       else if (pub_eda_mode == EDA_FROZENIDEM) then
          pub_eda_rep_x = (ex)*grid%weight
          pub_eda_rep_c = (ec)*grid%weight
       end if
    end if

    xc_energy = (ex+ec) * grid%weight

    ! qoh: Calculate potential if necessary

    if (present(xc_potint)) then
       potint = 0.0_DP
       do is=1,pub_num_spins
          do islab12=1,grid%num_my_slabs12
             do i2=1,grid%n2
                do i1=1,grid%n1
                   potint = potint + xc_pot_fine(i1,i2,islab12,is) * &
                        density_fine(i1,i2,islab12,is)
                end do
             end do
          end do
       end do
       if(loc_parallel) then
          call comms_reduce('SUM',potint)
       end if
       xc_potint = potint * grid%weight
    end if

  end subroutine xc_gc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_mgga(density_fine,density_grad,ke_density_fine, &
       recip_work,grid, &
       xc_energy,xc_pot_fine,dfdtau_fine,x_functional,x_functional_sp, &
       c_functional,c_functional_sp,xc_potint,exchange,correlation,x_pot_fine,&
       c_pot_fine, parallel)

    !==========================================================================!
    ! This subroutine calculates the exchange-correlation energy and potential !
    ! for a meta-gradient corrected functional, given the electronic density,  !
    ! gradient of the electronic density and kinetic energy density (tau) on   !
    ! the fine grid.  The exchange and correlation functionals to be used are  !
    ! given as arguments.                                                      !
    !--------------------------------------------------------------------------!
    ! Arguments                                                                !
    !   grid (in) : Grid definition                                            !
    !   density_fine (in) : Input density on the fine grid                     !
    !   xc_energy (out) : Exchange-correlation energy                          !
    !   xc_pot_fine (out) : Exchange-correlation potential (GGA-like part)     !
    !   dfdtau_fine (out) : Gradient of XC energy per unit volume with respect !
    !                       to KE-density (required to apply XC potential under!
    !                       FDO method)                                        !
    !   xc_potint (out) : Exchange-correlation potential-density integral      !
    !                     (does not include the KE-density-dependent part of   !
    !                     the potential)                                       !
    !   x_functional,x_functional_sp (in) : X subroutines given as arguments   !
    !   c_functional,c_functional_sp (in) : C subroutines given as arguments   !
    !   density_grad (inout) : Density gradient on the fine grid               !
    !   ke_density_fine (in) : Kinetic energy density on the fine grid         !
    !   recip_work (inout) : Fine Grid Recip-space workspace                   !
    !--------------------------------------------------------------------------!
    ! Modified version of xc_gc, created by James C. Womack, 2016.             !
    !--------------------------------------------------------------------------!
    ! Added spin polarization support, Jolyon Aarons, 2021.                    !
    !==========================================================================!
    ! Comments from original xc_gc routine:                                    !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in February / March 2009 based on existing       !
    ! xc_pbe subroutine written by Peter Haynes in February 2004.              !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce, pub_on_root
    use constants, only: DP, EDA_ISOLATED, EDA_FROZENIDEM, EDA_FROZENNONIDEM
    use rundat, only: pub_num_spins, pub_eda, pub_eda_mode, pub_frag_counter, &
         pub_eda_isol_x, pub_eda_isol_c, &
         pub_eda_frz_x, pub_eda_rep_x, &
         pub_eda_frz_c, pub_eda_rep_c, pub_eda_read_super
!$  use rundat, only: pub_threads_max
    use utils, only: utils_assert, utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(in) :: density_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP),intent(in) :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_num_spins)
    real(kind=DP),intent(in) :: ke_density_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    complex(kind=DP),intent(inout) :: recip_work(grid%ld3,grid%ld2, &
         grid%max_slabs23,3)
    real(kind=DP), intent(out) :: xc_energy
    real(kind=DP), intent(out) :: xc_pot_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(out) :: dfdtau_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(out) :: x_pot_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(out) :: c_pot_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(out) :: xc_potint
    real(kind=DP), optional, intent(out) :: exchange
    real(kind=DP), optional, intent(out) :: correlation
    logical,       optional, intent(in)  :: parallel

    interface
       subroutine x_functional(den,grad,mgd,tau,x_energy,x_pot,x_dfdmgd,&
            x_dfdgrad,x_dfdtau,use_dfdgrad)
         use constants, only: DP
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(in)  :: grad(3)
         real(kind=DP), intent(in)  :: mgd
         real(kind=DP), intent(in)  :: tau
         real(kind=DP), intent(out) :: x_energy
         real(kind=DP), intent(out) :: x_pot
         real(kind=DP), intent(out) :: x_dfdmgd
         real(kind=DP), intent(out) :: x_dfdgrad(3)
         real(kind=DP), intent(out) :: x_dfdtau
         logical,intent(out)        :: use_dfdgrad
       end subroutine x_functional
    end interface

    interface
       subroutine x_functional_sp(den1,den2,grad1,grad2,mgd1,mgd2,tau1,tau2,&
            x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2,&
            x_dfdgrad1,x_dfdgrad2,x_dfdtau1,x_dfdtau2,use_dfdgrad)
         use constants, only: DP
         real(kind=DP), intent(in)  :: den1
         real(kind=DP), intent(in)  :: den2
         real(kind=DP), intent(in)  :: grad1(3)
         real(kind=DP), intent(in)  :: grad2(3)
         real(kind=DP), intent(in)  :: mgd1
         real(kind=DP), intent(in)  :: mgd2
         real(kind=DP), intent(in)  :: tau1
         real(kind=DP), intent(in)  :: tau2
         real(kind=DP), intent(out) :: x_energy
         real(kind=DP), intent(out) :: x_pot1
         real(kind=DP), intent(out) :: x_pot2
         real(kind=DP), intent(out) :: x_dfdmgd1
         real(kind=DP), intent(out) :: x_dfdmgd2
         real(kind=DP), intent(out) :: x_dfdgrad1(3)
         real(kind=DP), intent(out) :: x_dfdgrad2(3)
         real(kind=DP), intent(out) :: x_dfdtau1
         real(kind=DP), intent(out) :: x_dfdtau2
         logical,intent(out)        :: use_dfdgrad
       end subroutine x_functional_sp
    end interface

    interface
       subroutine c_functional(den,grad,mgd,tau,c_energy,c_pot,c_dfdmgd,&
            c_dfdgrad,c_dfdtau,use_dfdgrad)
         use constants, only: DP
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(in)  :: grad(3)
         real(kind=DP), intent(in)  :: mgd
         real(kind=DP), intent(in)  :: tau
         real(kind=DP), intent(out) :: c_energy
         real(kind=DP), intent(out) :: c_pot
         real(kind=DP), intent(out) :: c_dfdmgd
         real(kind=DP), intent(out) :: c_dfdgrad(3)
         real(kind=DP), intent(out) :: c_dfdtau
         logical,intent(out)        :: use_dfdgrad
       end subroutine c_functional
    end interface

    interface
       subroutine c_functional_sp(den,den1,den2,grad1,grad2,mgd,mgd1,mgd2,tau1,tau2,&
            c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12,&
            c_dfdgrad1,c_dfdgrad2,c_dfdtau1,c_dfdtau2,use_dfdgrad)
         use constants, only: DP
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(in)  :: den1
         real(kind=DP), intent(in)  :: den2
         real(kind=DP), intent(in)  :: grad1(3)
         real(kind=DP), intent(in)  :: grad2(3)
         real(kind=DP), intent(in)  :: mgd
         real(kind=DP), intent(in)  :: mgd1
         real(kind=DP), intent(in)  :: mgd2
         real(kind=DP), intent(in)  :: tau1
         real(kind=DP), intent(in)  :: tau2
         real(kind=DP), intent(out) :: c_energy
         real(kind=DP), intent(out) :: c_pot1
         real(kind=DP), intent(out) :: c_pot2
         real(kind=DP), intent(out) :: c_dfdmgd
         real(kind=DP), intent(out) :: c_dfdmgd1
         real(kind=DP), intent(out) :: c_dfdmgd2
         real(kind=DP), intent(out) :: c_dfdmgd12
         real(kind=DP), intent(out) :: c_dfdgrad1(3)
         real(kind=DP), intent(out) :: c_dfdgrad2(3)
         real(kind=DP), intent(out) :: c_dfdtau1
         real(kind=DP), intent(out) :: c_dfdtau2
         logical,intent(out)        :: use_dfdgrad
       end subroutine c_functional_sp
    end interface

    ! Allocatable work arrays (store dfdgrad and divergence of dfdgrad)
    real(kind=DP),allocatable :: dwork(:,:,:,:)
    real(kind=DP),allocatable :: dwork_cart(:,:,:,:,:)

    real(kind=DP) :: den      ! Density
    real(kind=DP) :: den1     ! Density spin 1
    real(kind=DP) :: den2     ! Density spin 2
    real(kind=DP) :: mgd      ! Modulus of density gradient
    real(kind=DP) :: mgd1     ! Modulus of density gradient for spin 1
    real(kind=DP) :: mgd2     ! Modulus of density gradient for spin 2
    real(kind=DP) :: tau      ! Kinetic energy density
    real(kind=DP) :: tau1     ! Kinetic energy density spin 1
    real(kind=DP) :: tau2     ! Kinetic energy density spin 2
    real(kind=DP) :: grad(3)  ! Density gradient
    real(kind=DP) :: grad1(3) ! Density gradient for spin 1
    real(kind=DP) :: grad2(3) ! Density gradient for spin 2
    real(kind=DP) :: potint   ! Potential integral
    !qoh: Exchange variables
    real(kind=DP) :: ex        ! Accumulated exchange energy
    real(kind=DP) :: x_energy  ! Exchange energy at point
    real(kind=DP) :: x_pot     ! Exchange potential at point
    real(kind=DP) :: x_pot1    ! Exchange potential at point for spin 1
    real(kind=DP) :: x_pot2    ! Exchange potential at point for spin 2
    real(kind=DP) :: x_dfdmgd! Derivative df_{x}/d|grad n|
    real(kind=DP) :: x_dfdmgd1! Derivative df_{x}/d|grad n| spin 1
    real(kind=DP) :: x_dfdmgd2! Derivative df_{x}/d|grad n| spin 2
    real(kind=DP) :: x_dfdgrad1(3)! Derivative df_{x}/d(grad n_1)
    real(kind=DP) :: x_dfdgrad2(3)! Derivative df_{x}/d(grad n_2)
    real(kind=DP) :: x_dfdtau ! Derivative df_{x}/dtau
    real(kind=DP) :: x_dfdtau1! Derivative df_{x}/dtau spin 1
    real(kind=DP) :: x_dfdtau2! Derivative df_{x}/dtau spin 2
    logical       :: x_use_dfdgrad
    !qoh: Correlation variables
    real(kind=DP) :: ec        ! Accumulated correlation energy
    real(kind=DP) :: c_energy  ! Correlation energy at point
    real(kind=DP) :: c_pot     ! Correlation potential at point
    real(kind=DP) :: c_pot1    ! Correlation potential at point for spin 1
    real(kind=DP) :: c_pot2    ! Correlation potential at point for spin 2
    real(kind=DP) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=DP) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n| spin 1
    real(kind=DP) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n| spin 2
    real(kind=DP) :: c_dfdmgd12 ! df_{c}/d(grad n_1.grad n_2) component
    real(kind=DP) :: c_dfdgrad1(3)! Derivative df_{c}/d(grad n_1)
    real(kind=DP) :: c_dfdgrad2(3)! Derivative df_{c}/d(grad n_2)
    real(kind=DP) :: c_dfdtau  ! Derivative df_{c}/dtau
    real(kind=DP) :: c_dfdtau1 ! Derivative df_{c}/dtau spin 1
    real(kind=DP) :: c_dfdtau2 ! Derivative df_{c}/dtau spin 2
    logical       :: c_use_dfdgrad
    real(kind=DP) :: gd1_dot_gd2 ! (grad n_1).(grad n_2)
    integer :: i1,i2,islab12   ! Fine grid cartesian indices
    integer :: ipt             ! Fine grid loop counter
    integer :: is ! Spin counter
    integer :: ierr
    logical :: loc_parallel

    real(kind=dp) :: dummy

    loc_parallel = .true.
    if(present(parallel)) then
       loc_parallel = parallel
    end if



    ! Allocate simulation cell workspace arrays
    allocate(dwork(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),&
         stat=ierr)
    call utils_alloc_check('xc_mgga','dwork',ierr)
    allocate(dwork_cart(grid%ld1,grid%ld2,grid%max_slabs12,3,pub_num_spins),&
         stat=ierr)
    call utils_alloc_check('xc_mgga','dwork_cart',ierr)

    ex = 0.0_DP
    ec = 0.0_DP

    xc_pot_fine = 0.0_DP ! jd: To clear margins between n1 and ld1 etc.
    dfdtau_fine = 0.0_DP


    spin: if ( pub_num_spins == 1 ) then
       !qoh: Spin unpolarised case
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE (ipt,i1,i2,islab12,den,grad,mgd,tau,x_energy,c_energy, &
!$OMP      x_pot,c_pot,x_dfdmgd,x_dfdtau,x_dfdgrad1,x_use_dfdgrad,&
!$OMP      c_dfdmgd,c_dfdgrad1,c_dfdtau,c_use_dfdgrad) &
!$OMP SHARED (grid,density_fine,density_grad,dwork_cart,ke_density_fine,xc_pot_fine, &
!$OMP      x_pot_fine,c_pot_fine,dfdtau_fine,pub_threads_max) &
!$OMP REDUCTION(+:ex,ec)
       do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
          i1 = modulo(ipt-1,grid%n1) + 1
          i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
          islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1

          den  = density_fine(i1,i2,islab12,1)
          grad = density_grad(i1,i2,islab12,:,1)
          mgd  = sqrt(grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3))
          tau  = ke_density_fine(i1,i2,islab12,1)

          ! qoh: Calculate energy and potential at this point
          call x_functional(den, grad, mgd, tau, x_energy, x_pot, x_dfdmgd, &
               x_dfdgrad1, x_dfdtau, x_use_dfdgrad)
          call c_functional(den, grad, mgd, tau, c_energy, c_pot, c_dfdmgd, &
               c_dfdgrad1, c_dfdtau, c_use_dfdgrad)

          ! qoh: Store this point's energy and potential
          ex = x_energy + ex
          ec = c_energy + ec
          xc_pot_fine(i1,i2,islab12,1) = x_pot + c_pot
          ! JCW: Store dfdtau for this point
          dfdtau_fine(i1,i2,islab12,1) = x_dfdtau + c_dfdtau
          if(present(x_pot_fine)) then
             x_pot_fine(i1,i2,islab12,1) = x_pot
          endif
          if(present(c_pot_fine)) then
             c_pot_fine(i1,i2,islab12,1) = c_pot
          endif

          dwork_cart(i1,i2,islab12,:,1) = 0.0_DP
          if (x_use_dfdgrad) then
             ! Use dfdgrad directly
             dwork_cart(i1,i2,islab12,:,1) = -x_dfdgrad1(:)
          else if (mgd > 0.0_DP) then
             ! Form dfdgrad using dfdmgd
             dwork_cart(i1,i2,islab12,:,1) = &
                  -(x_dfdmgd) * grad(:) / mgd
          end if
          if (c_use_dfdgrad) then
             ! Use dfdgrad directly
             dwork_cart(i1,i2,islab12,:,1) = dwork_cart(i1,i2,islab12,:,1) &
                  -c_dfdgrad1(:)
          else if (mgd > 0.0_DP) then
             ! Form dfdgrad using dfdmgd
             dwork_cart(i1,i2,islab12,:,1) = dwork_cart(i1,i2,islab12,:,1) &
                  -(c_dfdmgd) * grad(:) / mgd
          end if

       end do
!$OMP END PARALLEL DO

    else
       !JA-> Spin polarised case

       c_dfdgrad2 = 0.0_DP
       x_dfdgrad2 = 0.0_DP
       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP
       x_dfdmgd1 = 0.0_DP
       x_dfdmgd2 = 0.0_DP
       c_dfdtau1 = 0.0_DP
       c_dfdtau2 = 0.0_DP
       x_dfdtau1 = 0.0_DP
       x_dfdtau2 = 0.0_DP
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP
       x_pot1 = 0.0_DP
       x_pot2 = 0.0_DP
       den1 = 0.0_DP
       den2 = 0.0_DP
       gd1_dot_gd2 = 0.0_DP
       grad1 = 0.0_DP
       grad2 = 0.0_DP
       mgd1 = 0.0_DP
       mgd2 = 0.0_DP
       tau1 = 0.0_DP
       tau2 = 0.0_DP
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE (ipt,i1,i2,islab12,den,den1,den2,grad,grad1,grad2, &
!$OMP      mgd,mgd1,mgd2,gd1_dot_gd2,x_energy,c_energy,x_pot1,x_pot2, &
!$OMP      c_pot1,c_pot2,x_dfdmgd,x_dfdmgd1,x_dfdmgd2,c_dfdmgd,c_dfdmgd1, &
!$OMP      x_dfdgrad1,x_dfdgrad2,tau1,tau2,c_dfdtau1,c_dfdtau2,c_dfdgrad1,c_dfdgrad2, &
!$OMP      c_dfdmgd2,c_dfdmgd12,x_dfdtau1,x_dfdtau2,x_use_dfdgrad,c_use_dfdgrad) &
!$OMP SHARED (grid,density_fine,density_grad,dwork_cart,xc_pot_fine,x_pot_fine, &
!$OMP      c_pot_fine,pub_threads_max,ke_density_fine,dfdtau_fine) &
!$OMP REDUCTION(+:ex,ec)
       do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
          i1 = modulo(ipt-1,grid%n1) + 1
          i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
          islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1

          den1 = density_fine(i1,i2,islab12,1)
          den2 = density_fine(i1,i2,islab12,2)
          den = den1 + den2

          grad1 = density_grad(i1,i2,islab12,:,1)
          grad2 = density_grad(i1,i2,islab12,:,2)
          grad = grad1 + grad2

          mgd1 = sqrt(grad1(1)*grad1(1) + grad1(2)*grad1(2) + &
               grad1(3)*grad1(3))
          mgd2 = sqrt(grad2(1)*grad2(1) + grad2(2)*grad2(2) + &
               grad2(3)*grad2(3))
          mgd = sqrt(grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3))
          gd1_dot_gd2 = grad1(1)*grad2(1) + grad1(2)*grad2(2) + &
               grad1(3)*grad2(3)

          tau1  = ke_density_fine(i1,i2,islab12,1)
          tau2  = ke_density_fine(i1,i2,islab12,2)

          ! JA-> Calculate energy and potential at this point
          call x_functional_sp(den1,den2,grad1,grad2,mgd1,mgd2,tau1,tau2,&
            x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2,&
            x_dfdgrad1,x_dfdgrad2,x_dfdtau1,x_dfdtau2,x_use_dfdgrad)

          call c_functional_sp(den,den1,den2,grad1,grad2,mgd,mgd1,mgd2,tau1,tau2,&
            c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12,&
            c_dfdgrad1,c_dfdgrad2,c_dfdtau1,c_dfdtau2,c_use_dfdgrad)


          ! JA-> Store this point's energy and potential
          ex = x_energy + ex
          ec = c_energy + ec
          xc_pot_fine(i1,i2,islab12,1) = x_pot1 + c_pot1
          xc_pot_fine(i1,i2,islab12,2) = x_pot2 + c_pot2

          ! JA-> Store dfdtau for this point
          dfdtau_fine(i1,i2,islab12,1) = x_dfdtau1 + c_dfdtau1
          dfdtau_fine(i1,i2,islab12,2) = x_dfdtau2 + c_dfdtau2
          if(present(x_pot_fine)) then
             x_pot_fine(i1,i2,islab12,1) = x_pot1
             x_pot_fine(i1,i2,islab12,2) = x_pot2
          endif
          if(present(c_pot_fine)) then
             c_pot_fine(i1,i2,islab12,1) = c_pot1
             c_pot_fine(i1,i2,islab12,2) = c_pot2
          endif

          dwork_cart(i1,i2,islab12,:,1) = 0.0_DP
          dwork_cart(i1,i2,islab12,:,2) = 0.0_DP
          if (x_use_dfdgrad) then
             ! JA -> Use dfdgrad directly
             dwork_cart(i1,i2,islab12,:,1) = -x_dfdgrad1(:)
             dwork_cart(i1,i2,islab12,:,2) = -x_dfdgrad2(:)
          else if (mgd > 0.0_DP) then
             ! JA->Form dfdgrad using dfdmgd
             dwork_cart(i1,i2,islab12,:,1) = &
                  -(x_dfdmgd1) * grad1(:) / mgd1
             dwork_cart(i1,i2,islab12,:,2) = &
                  -(x_dfdmgd2) * grad2(:) / mgd2
          end if
          if (c_use_dfdgrad) then
             ! JA-> Use dfdgrad directly
             dwork_cart(i1,i2,islab12,:,1) = dwork_cart(i1,i2,islab12,:,1) &
                  -c_dfdgrad1(:)
             dwork_cart(i1,i2,islab12,:,2) = dwork_cart(i1,i2,islab12,:,2) &
                  -c_dfdgrad2(:)
          else if (mgd > 0.0_DP) then
             ! JA-> Form dfdgrad using dfdmgd
             dwork_cart(i1,i2,islab12,:,1) = dwork_cart(i1,i2,islab12,:,1) &
                  -(c_dfdmgd1) * grad2(:) / mgd1
             dwork_cart(i1,i2,islab12,:,2) = dwork_cart(i1,i2,islab12,:,2) &
                  -(c_dfdmgd2) * grad2(:) / mgd2
          end if

       end do
!$OMP END PARALLEL DO
    end if spin


    ! Calculate divergence of gradient dependent part and add to potential
    call xc_divergence(dwork_cart,dwork,recip_work,grid)
    xc_pot_fine = xc_pot_fine + dwork


    ! Global reductions
    if(loc_parallel) then
       call comms_reduce('SUM',ex)
       call comms_reduce('SUM',ec)
    end if

    if(present(exchange)) then
    exchange  = ex * grid%weight
    endif
    if(present(correlation)) then
    correlation  = ec * grid%weight
    endif

    ! mjsp: if this is an isolated fragment calculation then
    ! update the fragment's stored X,C energies
    if (pub_eda) then
       if (pub_eda_mode == EDA_ISOLATED) then
          pub_eda_isol_x(pub_frag_counter) = (ex)*grid%weight
          pub_eda_isol_c(pub_frag_counter) = (ec)*grid%weight
       ! mjsp: if this is the frozen energy calculation
       ! then store the frozen state's X,C contributions
       ! (if pub_eda_read_super then this is contained in the .eda file)
       else if (pub_eda_mode == EDA_FROZENNONIDEM .and. .not.(pub_eda_read_super)) then
          pub_eda_frz_x = (ex)*grid%weight
          pub_eda_frz_c = (ec)*grid%weight
       ! mjsp: if this is the repulsion energy calculation
       ! then store the Pauli repulsion state's X,C contributions
       else if (pub_eda_mode == EDA_FROZENIDEM) then
          pub_eda_rep_x = (ex)*grid%weight
          pub_eda_rep_c = (ec)*grid%weight
       end if
    end if

    xc_energy = (ex+ec) * grid%weight

    ! qoh: Calculate potential if necessary


    if (present(xc_potint)) then
       potint = 0.0_DP
       do is=1,pub_num_spins
          do islab12=1,grid%num_my_slabs12
             do i2=1,grid%n2
                do i1=1,grid%n1
                   potint = potint + xc_pot_fine(i1,i2,islab12,is) * &
                        density_fine(i1,i2,islab12,is)
                end do
             end do
          end do
       end do
       if(loc_parallel) then
          call comms_reduce('SUM',potint)
       end if
       xc_potint = potint * grid%weight
    end if

    ! Deallocate simulation cell workspace arrays
    deallocate(dwork)
    call utils_dealloc_check('xc_mgga','dwork',ierr)
    deallocate(dwork_cart)
    call utils_dealloc_check('xc_mgga','dwork_cart',ierr)

  end subroutine xc_mgga

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef LIBXC
  subroutine xc_libxc(libxc_x_func_id,libxc_c_func_id,libxc_gga,libxc_mgga,&
       density_fine,recip_work,xc_energy,xc_pot_fine,grid,density_grad,&
       density_aux,ke_density_fine,dfdtau_fine,xc_potint,parallel)

    !==========================================================================!
    ! This subroutine calculates the exchange-correlation energy and potential !
    ! for a gradient corrected functional (including hybrids), given the       !
    ! electronic density on the fine grid.  The exchange and correlation       !
    ! functionals to be used are given as arguments.                           !
    !--------------------------------------------------------------------------!
    ! Arguments                                                                !
    !   grid (in) : Grid definition                                            !
    !   libxc_gga (in)  : .true. if functional is GGA, .false. otherwise       !
    !   libxc_mgga (in) : .true. if functional is meta-GGA, .false. otherwise  !
    !   density_fine (in) : Input density on the fine grid                     !
    !   xc_energy (out) : Exchange-correlation energy                          !
    !   xc_pot_fine (out) : Exchange-correlation potential                     !
    !   xc_potint (out) : Exchange-correlation potential-density integral      !
    !   density_grad (inout) : Density gradient on the fine grid               !
    !   density_aux (inout) : Auxiliary density on the fine grid               !
    !   ke_density_fine (in) : Kinetic energy density on the fine grid         !
    !   dfdtau_fine (out) : Gradient of XC energy per unit volume with respect !
    !                       to KE-density (required to apply XC potential under!
    !                       FDO method)                                        !
    !   recip_work (inout) : Fine Grid Recip-space workspace                   !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in July 2010, based partly on xc_gc by Quintin  !
    ! Hill.                                                                    !
    ! 21/10/15: As noted by James Womack, the arguments for x_func_id and      !
    ! c_func_id were originally reversed. Now fixed: hopefully this was        !
    ! harmless as both are treated identically within the code.                !
    ! Updated to Libxc F2003 interface introduced in version 3.0.0 by          !
    ! James C. Womack, 2015/2016.                                              !
    ! Updated to provide support for meta-GGA exchange-correlation functionals !
    ! by James C. Womack, 2015/2016.                                           !
    !==========================================================================!

    use iso_c_binding, only: c_double
    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce, pub_on_root
    use constants, only: DP, max_spins, stdout
#ifdef LIBXC
    use rundat, only: pub_libxc_x_func_id, pub_libxc_c_func_id
    use xc_f03_lib_m !! External dependency
#endif
    use rundat, only: pub_num_spins, pub_debug_on_root
    use utils, only: utils_abort, utils_assert
    implicit none

    ! Arguments
    integer, intent(in) :: libxc_c_func_id
    integer, intent(in) :: libxc_x_func_id
    logical, intent(in) :: libxc_gga  ! is GGA-type functional
    logical, intent(in) :: libxc_mgga ! is meta-GGA-type functional
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(in) :: density_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    complex(kind=DP),intent(inout) :: recip_work(grid%ld3,grid%ld2, &
         grid%max_slabs23,3)
    real(kind=DP), intent(out) :: xc_energy
    real(kind=DP), intent(out) :: xc_pot_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(inout) :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_num_spins)
    real(kind=DP), optional, intent(inout) :: density_aux(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(in) :: ke_density_fine(:,:,:,:)
    real(kind=DP), optional, intent(out) :: dfdtau_fine(:,:,:,:)
    ! JCW: Expect the following dimensions (or otherwise a dummy array with
    ! JCW: size == 1.
    ! JCW:    ke_density_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins)
    ! JCW:    dfdtau_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), optional, intent(out) :: xc_potint
    logical,       optional, intent(in)  :: parallel

    ! Local Variables
    real(kind=DP) :: potint   ! Potential integral
    ! JCW: Arrays with length of 1 required to satisfy libxc F2003 interface
    real(kind=c_double) :: den(max_spins) ! Density
    real(kind=c_double) :: grad(3,max_spins)  ! Density gradient
    real(kind=c_double) :: sigma(3) ! Modulus of density gradient
    real(kind=c_double) :: lapl_rho(max_spins) ! Laplacian of density
    real(kind=c_double) :: tau(max_spins)      ! Kinetic energy density
    real(kind=c_double) :: ex(1)               ! Accumulated xc energy
    real(kind=c_double) :: ec(1)               ! Accumulated xc energy
    real(kind=c_double) :: x_pot(max_spins)    ! xc potential at point
    real(kind=c_double) :: c_pot(max_spins)    ! xc potential at point
    real(kind=c_double) :: x_dfdsigma(3)       ! Derivative df_{xc}/dsigma
    real(kind=c_double) :: x_dfdlapl_rho(max_spins)
                                                  ! Derivative df_{xc}/dlapl_rho
    real(kind=c_double) :: x_dfdtau(max_spins) ! Derivative df_{xc}/dtau
    real(kind=c_double) :: c_dfdsigma(3)       ! Derivative df_{xc}/dsigma
    real(kind=c_double) :: c_dfdlapl_rho(max_spins)
                                                  ! Derivative df_{xc}/dlapl_rho
    real(kind=c_double) :: c_dfdtau(max_spins) ! Derivative df_{xc}/dtau
    real(kind=c_double), parameter :: two = real(2,kind=c_double)
    integer :: islab12 ! Loop counter
    integer :: i1 ! Loop counter
    integer :: i2 ! Loop counter
    integer :: is ! Spin counter
    integer :: ns ! Spin counter
    logical :: loc_parallel

    ! JCW: Check that libxc_gga and libxc_mgga are not both true (these should
    ! JCW: be mutually exclusive).
    call utils_assert(.not.(libxc_gga.and.libxc_mgga), &
         "Error in xc_libxc: Functional cannot be both GGA-type and meta-GGA-&
         &type (libxc_gga and libxc_mgga cannot both be true).")

    loc_parallel = .true.
    if(present(parallel)) then
       loc_parallel = parallel
    end if



    ! JCW: Ensure that optional arguments are consistent with value
    ! JCW: of libxc_mgga.
    if (libxc_mgga) then
       ! Meta-GGA (tau-dependent)
       ! JCW: TODO add support for lapl_rho dependent functionals
       call utils_assert(present(ke_density_fine),"Error in xc_libxc: &
            &For a tau-dependent meta-GGA, ke_density_fine must be present.")
       call utils_assert(present(dfdtau_fine),"Error in xc_libxc: &
            &For a tau-dependent meta-GGA, dfdtau_fine must be present.")
    else
       ! LDA or GGA
       if ( present(ke_density_fine) ) then
          call utils_assert( &
               size(ke_density_fine) == 1,&
               "Error in xc_libxc: &
               &For a LDA or GGA, ke_density_fine must not be present &
               &or have size == 1.")
       end if
       if ( present(dfdtau_fine ) ) then
          call utils_assert( &
               size(dfdtau_fine) == 1,&
               "Error in xc_libxc: &
               &For a LDA or GGA, dfdtau_fine must not be present &
               &or have size == 1.")
       end if
    end if

    ns = pub_num_spins

    ! JCW: Initialize local variables which are modified by Libxc to 0.0
    ! JCW: for safety.
    den(:)           = 0.0_c_double
    grad(:,:)        = 0.0_c_double
    sigma(:)         = 0.0_c_double
    lapl_rho(:)      = 0.0_c_double
    tau(:)           = 0.0_c_double
    ex(1)            = 0.0_c_double
    ec(1)            = 0.0_c_double
    x_pot(:)         = 0.0_c_double
    c_pot(:)         = 0.0_c_double
    x_dfdsigma(:)    = 0.0_c_double
    x_dfdlapl_rho(:) = 0.0_c_double
    x_dfdtau(:)      = 0.0_c_double
    c_dfdsigma(:)    = 0.0_c_double
    c_dfdlapl_rho(:) = 0.0_c_double
    c_dfdtau(:)      = 0.0_c_double

    xc_energy = 0.0_DP
    xc_pot_fine = 0.0_DP ! jd: To clear margins between n1 and ld1 etc.

    a3: do islab12=1,grid%num_my_slabs12
       a2: do i2=1,grid%n2
          a1: do i1=1,grid%n1

             den(1:ns) = real(density_fine(i1,i2,islab12,1:ns),kind=c_double)
             do is=1,ns
                if (den(is)<0.0_c_double) den(is) = 0.0_c_double
             end do

             ! ndmh: Calculate energy and potential at this point
             if (.not.(libxc_gga.or.libxc_mgga)) then  ! Non-Gradient Corrected

                ! ndmh: Exchange part
                if (libxc_x_func_id > 0) then
                   call xc_f03_lda_exc_vxc(libxc_func(1), 1, den(1), &
                        ex(1), x_pot(1))
                else
                   ex(1) = 0.0_c_double; x_pot(:) = 0.0_c_double
                end if

                ! ndmh: Correlation part
                if (libxc_c_func_id > 0) then
                   call xc_f03_lda_exc_vxc(libxc_func(2), 1, den(1), &
                        ec(1), c_pot(1))
                else
                   ec(1) = 0.0_c_double; c_pot(:) = 0.0_c_double
                end if

             else if (libxc_gga.and.(.not.libxc_mgga)) then
                ! Gradient Corrected (GGA, not meta-GGA)
                grad(1:3,1:ns) = real(density_grad(i1,i2,islab12,1:3,1:ns),&
                    kind=c_double)

                ! Calculate reduced gradients
                if (ns==2) then
                   sigma(1) = grad(1,1)*grad(1,1) + &
                        grad(2,1)*grad(2,1)+grad(3,1)*grad(3,1)
                   sigma(2) = grad(1,1)*grad(1,2) + &
                        grad(2,1)*grad(2,2)+grad(3,1)*grad(3,2)
                   sigma(3) = grad(1,2)*grad(1,2) + &
                        grad(2,2)*grad(2,2)+grad(3,2)*grad(3,2)
                else
                   sigma(1) = grad(1,1)*grad(1,1) + &
                        grad(2,1)*grad(2,1)+grad(3,1)*grad(3,1)
                end if

                ! ndmh: Exchange part
                if (libxc_x_func_id > 0) then
                   call xc_f03_gga_exc_vxc(libxc_func(1), 1, den(1), &
                        sigma(1), ex(1), x_pot(1), x_dfdsigma(1))
                else
                   ex(1) = 0.0_c_double; x_pot(:) = 0.0_c_double; x_dfdsigma(:) = 0.0_c_double
                end if

                ! ndmh: Correlation part
                if (libxc_c_func_id > 0) then
                   call xc_f03_gga_exc_vxc(libxc_func(2), 1, den(1), &
                        sigma(1), ec(1), c_pot(1), c_dfdsigma(1))
                else
                   ec(1) = 0.0_c_double; c_pot(:) = 0.0_c_double; c_dfdsigma(:) = 0.0_c_double
                end if

                if (ns==2) then
                   !x_dfdsigma(1) = x_dfdsigma(1)*2.0_DP*sqrt(sigma(1)) &
                   !     + x_dfdsigma(2)*sqrt(sigma(3))
                   !x_dfdsigma(3) = x_dfdsigma(3)*2.0_DP*sqrt(sigma(3)) &
                   !     + x_dfdsigma(2)*sqrt(sigma(1))
                   !c_dfdsigma(1) = c_dfdsigma(1)*2.0_DP*sqrt(sigma(1)) &
                   !     + c_dfdsigma(2)*sqrt(sigma(3))
                   !c_dfdsigma(3) = c_dfdsigma(3)*2.0_DP*sqrt(sigma(3)) &
                   !     + c_dfdsigma(2)*sqrt(sigma(1))

                   if (sigma(1) > 0.0_c_double) then
                      density_grad(i1,i2,islab12,:,1) = real( &
                           -two*(x_dfdsigma(1)+c_dfdsigma(1)) * grad(:,1) &
                           -(x_dfdsigma(2)+c_dfdsigma(2)) * grad(:,2), kind=DP)
                   else
                      density_grad(i1,i2,islab12,:,1) = 0.0_DP
                   end if
                   if (sigma(3) > 0.0_c_double) then
                      density_grad(i1,i2,islab12,:,2) = real( &
                           -two*(x_dfdsigma(3)+c_dfdsigma(3)) * grad(:,2) &
                           -(x_dfdsigma(2)+c_dfdsigma(2)) * grad(:,1), kind=DP)
                   else
                      density_grad(i1,i2,islab12,:,2) = 0.0_DP
                   end if

                else
                   density_grad(i1,i2,islab12,:,1) = real( &
                        -two*(x_dfdsigma(1)+c_dfdsigma(1)) * grad(:,1), kind=DP)
                end if

             else if (libxc_mgga.and.(.not.libxc_gga)) then
                ! Gradient Corrected (meta-GGA, not GGA)

                ! JCW: Set lapl_rho 0.0_DP (functionals dependent on Laplacian
                ! JCW: of the density are not yet supported).
                lapl_rho(1:ns) = 0.0_c_double

                ! JCW: Get ke_density at point i1,i2,islab12
                tau(1:ns) = real(ke_density_fine(i1,i2,islab12,1:ns),&
                    kind=c_double)

                ! JCW: Get density gradient at point i1,i2,islab12
                grad(1:3,1:ns) = real(density_grad(i1,i2,islab12,1:3,1:ns),&
                    kind=c_double)

                ! Calculate reduced gradients
                if (ns==2) then
                   sigma(1) = grad(1,1)*grad(1,1) + &
                        grad(2,1)*grad(2,1)+grad(3,1)*grad(3,1)
                   sigma(2) = grad(1,1)*grad(1,2) + &
                        grad(2,1)*grad(2,2)+grad(3,1)*grad(3,2)
                   sigma(3) = grad(1,2)*grad(1,2) + &
                        grad(2,2)*grad(2,2)+grad(3,2)*grad(3,2)
                else
                   sigma(1) = grad(1,1)*grad(1,1) + &
                        grad(2,1)*grad(2,1)+grad(3,1)*grad(3,1)
                end if

                ! JCW: Exchange functional
                if (libxc_x_func_id > 0) then
                   call xc_f03_mgga_exc_vxc(libxc_func(1), 1, den(1), &
                        sigma(1), lapl_rho(1), tau(1), ex(1), x_pot(1), &
                        x_dfdsigma(1), x_dfdlapl_rho(1), x_dfdtau(1) )
                else
                   ex(1) = 0.0_c_double
                   x_pot(:) = 0.0_c_double
                   x_dfdsigma(:) = 0.0_c_double
                   x_dfdlapl_rho(:) = 0.0_c_double
                   x_dfdtau(:) = 0.0_c_double
                end if

                ! JCW: Correlation functional
                if (libxc_c_func_id > 0) then
                   call xc_f03_mgga_exc_vxc(libxc_func(2), 1, den(1), &
                        sigma(1), lapl_rho(1), tau(1), ec(1), c_pot(1), &
                        c_dfdsigma(1), c_dfdlapl_rho(1), c_dfdtau(1) )
                else
                   ec(1) = 0.0_c_double
                   c_pot(:) = 0.0_c_double
                   c_dfdsigma(:) = 0.0_c_double
                   c_dfdlapl_rho(:) = 0.0_c_double
                   c_dfdtau(:) = 0.0_c_double
                end if

                if (ns==2) then
                   !x_dfdsigma(1) = x_dfdsigma(1)*2.0_DP*sqrt(sigma(1)) &
                   !     + x_dfdsigma(2)*sqrt(sigma(3))
                   !x_dfdsigma(3) = x_dfdsigma(3)*2.0_DP*sqrt(sigma(3)) &
                   !     + x_dfdsigma(2)*sqrt(sigma(1))
                   !c_dfdsigma(1) = c_dfdsigma(1)*2.0_DP*sqrt(sigma(1)) &
                   !     + c_dfdsigma(2)*sqrt(sigma(3))
                   !c_dfdsigma(3) = c_dfdsigma(3)*2.0_DP*sqrt(sigma(3)) &
                   !     + c_dfdsigma(2)*sqrt(sigma(1))

                   if (sigma(1) > 0.0_c_double) then
                      density_grad(i1,i2,islab12,:,1) = real(&
                           -two*(x_dfdsigma(1)+c_dfdsigma(1)) * grad(:,1) &
                           -(x_dfdsigma(2)+c_dfdsigma(2)) * grad(:,2),kind=DP)
                   else
                      density_grad(i1,i2,islab12,:,1) = 0.0_DP
                   end if
                   if (sigma(3) > 0.0_c_double) then
                      density_grad(i1,i2,islab12,:,2) = &
                           real(-two*(x_dfdsigma(3)+c_dfdsigma(3)) * grad(:,2) &
                           -(x_dfdsigma(2)+c_dfdsigma(2)) * grad(:,1),kind=DP)
                   else
                      density_grad(i1,i2,islab12,:,2) = 0.0_DP
                   end if

                else
                   density_grad(i1,i2,islab12,:,1) = &
                        real(-two*(x_dfdsigma(1)+c_dfdsigma(1)) * grad(:,1),&
                        kind=DP)
                end if

                ! JCW: Accumulate x_dfdtau and c_dfdtau into dfdtau_fine
                dfdtau_fine(i1,i2,islab12,1:ns) = real(&
                     x_dfdtau(1:ns) + c_dfdtau(1:ns), kind=DP)
             else
                ! JCW: Neither or both libxc_gga and libxc_mgga true.
                ! JCW: Undefined state, so abort.
                call utils_abort(&
                     "Error in xc_libxc: unknown functional category.")

             end if

             ! ndmh: Store this point's energy and potential
             xc_energy = xc_energy + (real(ex(1),kind=DP) + real(ec(1),kind=DP))*&
                  real(sum(den(1:ns)),kind=DP)
             !c_energy = c_energy + ec*sum(den(1:ns))
             !x_energy = x_energy + ex*sum(den(1:ns))
             xc_pot_fine(i1,i2,islab12,1:ns) = real(x_pot(1:ns),kind=DP) + &
                  real(c_pot(1:ns),kind=DP)
          end do a1
       end do a2
    end do a3

    if (libxc_gga.or.libxc_mgga) then
       ! Calculate divergence of gradient dependent part and add to potential
       call xc_divergence(density_grad,density_aux,recip_work,grid)
       xc_pot_fine = xc_pot_fine + density_aux
    end if

    ! Global reductions
    if(loc_parallel) then
       call comms_reduce('SUM',xc_energy)
    end if
    xc_energy = xc_energy * grid%weight

    ! qoh: Calculate potential if necessary
    if (present(xc_potint)) then
       potint = 0.0_DP
       do is=1,pub_num_spins
          do islab12=1,grid%num_my_slabs12
             do i2=1,grid%n2
                do i1=1,grid%n1
                   potint = potint + xc_pot_fine(i1,i2,islab12,is) * &
                        density_fine(i1,i2,islab12,is)
                end do
             end do
          end do
       end do
       if(loc_parallel) then
          call comms_reduce('SUM',potint)
       end if
       xc_potint = potint * grid%weight
    end if

  end subroutine xc_libxc
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_gradients(density_fine,density_grad,recip_work,grid, &
       nhat_dim,nhat_den_grad)

    !========================================================================!
    ! This subroutine calculates the gradients of the density in preparation !
    ! for a gradient-corrected functional.                                   !
    !------------------------------------------------------------------------!
    ! Arguments                                                              !
    !   grid (in) : Grid definition                                          !
    !   density_fine (in) : Input density on the fine grid                   !
    !   density_grad (inout) : Density gradient on the fine grid             !
    !   recip_work (inout) : Fine grid recip-space workspace                 !
    !   nhat_den_grad (in) : Compensation density (and gradient) on fine grid!
    !------------------------------------------------------------------------!
    ! Originally written by Peter Haynes in February 2004.                   !
    ! Modified to avoid module-level arrays by Nicholas Hine in October 2010.!
    !========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_proc_id
    use constants, only: cmplx_i
    use fourier, only: fourier_apply_cell_forward, fourier_apply_cell_backward
    use rundat, only: pub_nhat_in_xc, pub_num_spins

    implicit none

    ! Arguments
    type(GRID_INFO),intent(in) :: grid
#ifdef WIN32
    real(kind=DP), intent(inout),target :: density_fine(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    complex(kind=DP),intent(inout),target :: recip_work(grid%ld3,grid%ld2, &
         grid%max_slabs23,3)
    real(kind=DP),intent(inout),target :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_num_spins)
#else
    real(kind=DP), intent(inout) :: density_fine(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    complex(kind=DP),intent(inout) :: recip_work(grid%ld3,grid%ld2, &
         grid%max_slabs23,3)
    real(kind=DP),intent(inout) :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_num_spins)
#endif
    integer, intent(in) :: nhat_dim
    real(kind=DP), optional, intent(in) :: nhat_den_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:nhat_dim)

    ! Local Variables
    integer :: is                      ! Loop variable over spin component
    integer :: i3,i2,islab23           ! Loop variables over fine grids
    integer :: k                       ! Loop variables for Cartesian components
    complex(kind=DP) :: dengtmp        ! Temporary copies of reciprocal density
    real(kind=DP) :: gvec(3)
    ! Introduce pointers to avoid copies to stack on Win32
#ifdef WIN32
    real(kind=DP), pointer :: dens_real(:,:,:)
    complex(kind=DP), pointer :: dens_cmpl(:,:,:)
#endif

    spin: do is=1, pub_num_spins

       ! ndmh: subtract off the compensation density (if present) before
       ! ndmh: taking gradient: analytic gradient will be added afterwards
       if (present(nhat_den_grad).and.pub_nhat_in_xc) then
          density_fine(:,:,:,is) = density_fine(:,:,:,is) - &
               nhat_den_grad(:,:,:,is,0)
       end if

       ! Fourier transform to reciprocal space
#ifdef WIN32
       dens_real => density_fine(:,:,:,is)
       dens_cmpl => recip_work(:,:,:,1)
       call fourier_apply_cell_forward(dens_real,dens_cmpl,grid)
#else
       call fourier_apply_cell_forward(density_fine(:,:,:,is), &
            recip_work(:,:,:,1),grid)
#endif

       ! Apply grad in reciprocal space (multiply by iG)
       b1: do islab23=1,grid%num_slabs23
          b2: do i2=1,grid%n2
             b3: do i3=1,grid%n3

                ! Apply grad
                dengtmp = recip_work(i3,i2,islab23,1)

                call cell_grid_recip_pt(gvec,islab23 + &
                     grid%first_slab23(pub_my_proc_id) - 1,i2,i3,grid)

                cart_comps1: do k=1,3
                   recip_work(i3,i2,islab23,k) = dengtmp * cmplx_i * gvec(k)
                end do cart_comps1

             end do b3
          end do b2
       end do b1

       ! Transform gradient back to real space
       cart_comps2: do k=1,3
#ifdef WIN32
          dens_real => density_grad(:,:,:,k,is)
          dens_cmpl => recip_work(:,:,:,k)
          call fourier_apply_cell_backward(dens_real,dens_cmpl,grid)
#else
          call fourier_apply_cell_backward(density_grad(:,:,:,k,is), &
               recip_work(:,:,:,k),grid)
#endif
       end do cart_comps2

       ! ndmh: add back on the compensation density itself and add on
       ! ndmh: the compensation density gradient (if present)
       if (present(nhat_den_grad).and.pub_nhat_in_xc) then
          density_fine(:,:,:,is) = density_fine(:,:,:,is) + &
               nhat_den_grad(:,:,:,is,0)
          do k=1,3
             density_grad(:,:,:,k,is) = density_grad(:,:,:,k,is) + &
                  nhat_den_grad(:,:,:,is,k)
          end do
       end if

    end do spin

  end subroutine xc_gradients

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_divergence(density_grad,density_aux,recip_work,grid)

    !=========================================================================!
    ! This subroutine calculates the divergence of a vector field (stored in  !
    ! the module array density_grad) which is returned in the array           !
    ! density_aux.                                                            !
    !-------------------------------------------------------------------------!
    ! Originally written by Peter Haynes in February 2004.                    !
    ! Modified to avoid module-level arrays by Nicholas Hine in October 2010. !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_proc_id
    use constants, only: DP, cmplx_i
    use fourier, only: fourier_apply_cell_forward, fourier_apply_cell_backward
    use rundat, only: pub_num_spins

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
#ifdef WIN32
    real(kind=DP),intent(inout),target :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_num_spins)
    complex(kind=DP),intent(inout),target :: recip_work(grid%ld3,grid%ld2, &
         grid%max_slabs23,3)
    real(kind=DP),intent(inout),target :: density_aux(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
#else
    real(kind=DP),intent(inout) :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_num_spins)
    complex(kind=DP),intent(inout) :: recip_work(grid%ld3,grid%ld2, &
         grid%max_slabs23,3)
    real(kind=DP),intent(inout) :: density_aux(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
#endif

    ! Local Variables
    integer :: is                      ! Loop variable over spin component
    integer :: i2,i3,islab23           ! Loop variables over fine grids
    integer :: k                       ! Loop variables for Cartesian components
    complex(kind=DP) :: div
    real(kind=DP) :: gvec(3)
#ifdef WIN32
    real(kind=DP), pointer :: dens_real(:,:,:)
    complex(kind=DP), pointer :: dens_cmpl(:,:,:)
#endif

    spin: do is=1, pub_num_spins

       ! Loop over Cartesian components of the vector field
       cart_comps1: do k=1,3

          ! Fourier transform to reciprocal space
#ifdef WIN32
          dens_real => density_grad(:,:,:,k,is)
          dens_cmpl => recip_work(:,:,:,k)
          call fourier_apply_cell_forward(dens_real,dens_cmpl,grid)
#else
          call fourier_apply_cell_forward(density_grad(:,:,:,k,is), &
               recip_work(:,:,:,k),grid)
#endif
       end do cart_comps1

       ! Apply divergence in reciprocal space (multiply by iG)

       b1: do islab23=1,grid%num_slabs23
          b2: do i2=1,grid%n2
             b3: do i3=1,grid%n3

                ! Calculate divergence
                div = 0.0_DP

                call cell_grid_recip_pt(gvec,islab23 + &
                     grid%first_slab23(pub_my_proc_id) - 1,i2,i3,grid)

                ! Loop over Cartesian components of the vector field
                cart_comps2: do k=1,3

                   div = div + recip_work(i3,i2,islab23,k) * cmplx_i * gvec(k)

                end do cart_comps2

                recip_work(i3,i2,islab23,1) = div

             end do b3
          end do b2
       end do b1

       ! Transform divergence back to real space
#ifdef WIN32
       dens_real => density_aux(:,:,:,is)
       dens_cmpl => recip_work(:,:,:,1)
       call fourier_apply_cell_backward(dens_real,dens_cmpl,grid)
#else
       call fourier_apply_cell_backward(density_aux(:,:,:,is),&
            recip_work(:,:,:,1),grid)
#endif
    end do spin

  end subroutine xc_divergence

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_condition_density(conditioned_density_out,density_in,grid)

    !=========================================================================!
    ! This subroutine conditions density_fine by zeroing all densities below  !
    ! d0 (1E-9), including negative densities, leaving all densities above d1 !
    ! (1E-7) intact, and applying a gradual transition in between.            !
    ! This helps avoid ill-conditioning in certain XC functionals (e.g. BLYP),!
    ! where extreme noise in densities close to dentol (1E-15) leads to large !
    ! values of gradient terms.                                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   conditioned_density_out (out): Density after conditioning.            !
    !   density_in (in): Density to be conditioned.                           !
    !   grid (in): Describes the slab-parallel distribution of the density.   !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 2014/07/01.                                !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: DP
    use rundat, only: pub_num_spins

    implicit none

    ! jd: Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(out)  :: conditioned_density_out(&
         grid%ld1, grid%ld2, grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in)   :: density_in(&
         grid%ld1, grid%ld2, grid%max_slabs12, pub_num_spins)

    ! jd: Local variables
    integer       :: i1, i2, islab12, is
    real(kind=DP) :: den_in, den_out

    ! jd: Parameters of the filter. The filter is a linear function with a value
    !     of 0 at d0 and d1 at d1. Curiously, a cubic polynomial with zero
    !     derivatives at d0 and d1 did not work as well.
    real(kind=DP) :: d0
    real(kind=DP) :: d1
    d0 = 1D-9  / real(pub_num_spins,kind=DP)
    d1 = 1D-7  / real(pub_num_spins,kind=DP)

    ! -------------------------------------------------------------------------

    conditioned_density_out = 0D0 ! jd: Takes care of ld and num/max padding

    do is=1, pub_num_spins
       do islab12=1, grid%num_my_slabs12
          do i2=1, grid%n2
             do i1=1, grid%n1
                den_in = density_in(i1,i2,islab12,is)
                if(den_in > d1) then
                   den_out = den_in
                else if(den_in <= d0) then
                   den_out = 0D0
                else
                   den_out = d1/(d1-d0)*den_in + d0*d1/(d0-d1)
                end if
                conditioned_density_out(i1,i2,islab12,is) = den_out
             end do
          end do
       end do
    end do

  end subroutine xc_condition_density

  subroutine xc_condition_ke_density(conditioned_ke_density_out,ke_density_in,grid)
    !=========================================================================!
    ! This subroutine conditions ke_density_fine by zeroing all negative KE   !
    ! densities below (<0.0_DP)                                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   conditioned_ke_density_out (out): Density after conditioning.         !
    !   ke_density_in (in): Density to be conditioned.                        !
    !   grid (in): Describes the slab-parallel distribution of the density.   !
    !-------------------------------------------------------------------------!
    ! James C. Womack, 2015/2016                                              !
    !=========================================================================!
    use cell_grid, only: GRID_INFO
    use constants, only: DP
    use rundat, only: pub_num_spins
    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(out)  :: conditioned_ke_density_out(&
         grid%ld1, grid%ld2, grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(in)   :: ke_density_in(&
         grid%ld1, grid%ld2, grid%max_slabs12, pub_num_spins)

    ! Local variables
    integer       :: i1, i2, islab12, is
    real(kind=DP) :: den_in, den_out
    conditioned_ke_density_out = 0.0_DP

    do is=1, pub_num_spins
       do islab12=1, grid%num_my_slabs12
          do i2=1, grid%n2
             do i1=1, grid%n1
                den_in = ke_density_in(i1,i2,islab12,is)
                if (den_in<0.0_DP) then
                   den_out = 0.0_DP
                else
                   den_out = den_in
                end if
                conditioned_ke_density_out(i1,i2,islab12,is) = den_out
             end do
          end do
       end do
    end do

  end subroutine xc_condition_ke_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_embed_swap_functional(switch,use_tddft)
    !=========================================================================!
    ! This subroutine swaps functional from the main functional to the active !
    ! region functional (or vice versa) if switch is true (false)             !
    !-------------------------------------------------------------------------!
    ! Written by Joseph Prentice, June 2018                                   !
    ! Added in switch for swapping TDDFT functional, July 2018                !
    !=========================================================================!

    implicit none

    ! Arguments
    logical, intent(in) :: switch
    logical, optional, intent(in) :: use_tddft

    logical :: loc_tddft

    loc_tddft=.false.
    if (present(use_tddft)) loc_tddft=use_tddft

    call xc_exit
    if (loc_tddft) then
       call xc_init(use_tddftxc=loc_tddft,use_activexc=switch)
    else
       call xc_init(use_activexc=switch)
    end if
    call xc_hfxinit(active_switch=switch)

  end subroutine xc_embed_swap_functional

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_emft_calculate(density_fine, xc_energy_emft, &
       xc_emft_potential, grid, cell)

    !=========================================================================!
    ! This subroutine calculates the correction to the XC potential and energy!
    ! with Embedded Mean Field Theory (EMFT), based on the work of Fornace    !
    ! et al, 2015 (doi:10.1021/ct5011032).                                    !
    !-------------------------------------------------------------------------!
    ! Written by Robert Charlton, June 2018.                                  !
    ! Based on original code in hamiltonian_mod by Joseph Prentice.           !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use rundat, only: pub_num_spins
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    real(kind=DP), intent(in)   :: density_fine(&
         grid%ld1, grid%ld2, grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(out)  :: xc_emft_potential(&
         grid%ld1, grid%ld2, grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(out)  :: xc_energy_emft

    ! Local variables
    real(kind=DP) :: xc_energy_low
    real(kind=DP), allocatable :: buffer_fine_pot(:,:,:,:)
    integer :: ierr

    ! rc2013: allocate space for the low-level potential
    allocate(buffer_fine_pot(grid%ld1, grid%ld2, grid%max_slabs12, &
         pub_num_spins), stat=ierr)
    call utils_alloc_check('xc_emft_calculate','buffer_fine_pot',ierr)

    ! jcap: xc part of the hamiltonian is the only part that
    ! changes between different levels of DFT theory. First,
    ! calculate xc for low level theory for active subregion

    call xc_energy_potential(density_fine, xc_energy_low, &
         buffer_fine_pot, grid, cell, 0)

    ! jcap: need to swap to high level xc functional, so call
    ! xc_exit and xc_init again, with optional argument to use
    ! active subregion functional (bit of a hack)
    call xc_embed_swap_functional(.true.)

    ! jcap: calculate xc for high level theory for active subregion
    call xc_energy_potential(density_fine, xc_energy_emft, &
         xc_emft_potential, grid, cell, 0)

    ! jcap: need to swap the xc functional back
    call xc_embed_swap_functional(.false.)

    ! rc2013: calculate the EMFT energy correction
    xc_energy_emft = xc_energy_emft - xc_energy_low

    ! rc2013: subtract the low-level potential from the higher level
    xc_emft_potential = xc_emft_potential - buffer_fine_pot

    ! rc2013: clean-up
    deallocate(buffer_fine_pot, stat=ierr)
    call utils_dealloc_check('xc_emft_calculate','buffer_fine_pot',ierr)

  end subroutine xc_emft_calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_fxc_potential_emft(grid,fxc_pot_fine,density_fine)

    !=========================================================================!
    !Subroutine getting the EMFT correction to fxc in fine grid representation!
    !from the ground state density of the converged system. Returns           !
    !fxc_pot_fine. Currently only capable of doing ALDA in spin degenerate    !
    !systems, as this is all xc_fxc_potential can do                          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !    density_fine (in) : Input active region density on the fine grid     !
    !    grid (in)         : GRID_INFO grid definition                        !
    !    fxc_pot_fine (out) : EMFT correction to exchange-correlation kernel  !
    !                         on grid                                         !
    !-------------------------------------------------------------------------!
    ! Written by Joseph Prentice, July 2018                                   !
    ! Based on xc_fxc_potential, written by Tim Zuhlsdorff, and               !
    ! xc_emft_calculate, written by Rob Charlton                              !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use rundat, only: pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    type (GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: fxc_pot_fine(grid%ld1,&
         grid%ld2,grid%max_slabs12,2)
    real(kind=DP), intent(inout) :: density_fine(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)

    real(kind=DP), allocatable :: buffer_fine_pot(:,:,:,:)
    integer :: ierr

    ! rc2013: allocate space for the low-level potential
    allocate(buffer_fine_pot(grid%ld1, grid%ld2, grid%max_slabs12, &
         pub_num_spins), stat=ierr)
    call utils_alloc_check('xc_fxc_potential_emft','buffer_fine_pot',ierr)

    ! jcap: xc part of the hamiltonian is the only part that
    ! changes between different levels of DFT theory. First,
    ! calculate xc for low level theory for active subregion

    call xc_fxc_potential(grid, buffer_fine_pot, density_fine)

    ! jcap: need to swap to high level xc functional, so call
    ! xc_exit and xc_init again, with optional argument to use
    ! active subregion functional (bit of a hack)
    call xc_embed_swap_functional(.true.)

    ! jcap: calculate xc for high level theory for active subregion
    call xc_fxc_potential(grid, fxc_pot_fine, density_fine)

    ! jcap: need to swap the xc functional back
    call xc_embed_swap_functional(.false.)

    ! rc2013: subtract the low-level potential from the higher level
    fxc_pot_fine = fxc_pot_fine - buffer_fine_pot

    ! rc2013: clean-up
    deallocate(buffer_fine_pot, stat=ierr)
    call utils_dealloc_check('xc_fxc_potential_emft','buffer_fine_pot',ierr)

  end subroutine xc_fxc_potential_emft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_fxc_finite_diff_emft(density_fine,fxc_pot_fine,&
         grid,cell,ground_state_dens,nhat_dim,eff_nhat_den,nhat_den_grad_gs)

    !===========================================================================!
    ! Given the electronic ground state density and the perturbed density on    !
    ! finge grid, this routine calculates the EMFT correction to fxc *          !
    ! rho_perturbed using only Vxc and a finite difference technique. This      !
    ! prevents us from having to code up the actual second derivatives of the   !
    ! exchange correlation energy, which for GGA codes is cumbersome and        !
    ! numerically unstable.                                                     !
    ! Written by Joseph Prentice, July 2018                                     !
    ! Based on xc_fxc_finite_diff, written by Tim Zuhlsdorff, and               !
    ! xc_emft_calculate, written by Rob Charlton                                !
    !---------------------------------------------------------------------------!

    use cell_grid, only: GRID_INFO
    use simulation_cell, only: CELL_INFO
    use rundat, only: pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid    ! Grid definition
    type(CELL_INFO), intent(in) :: cell    ! Cell definition
    real(kind=DP), intent(inout) :: density_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(out) :: fxc_pot_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,2)
    real(kind=DP), intent(in) :: ground_state_dens(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    integer, intent(in) :: nhat_dim
    real(kind=DP), optional, intent(in) :: eff_nhat_den(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:nhat_dim)
    real(kind=DP), optional, intent(in) :: nhat_den_grad_gs(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins,0:nhat_dim)

    real(kind=DP), allocatable :: buffer_fine_pot(:,:,:,:,:)
    integer :: ierr

    ! rc2013: allocate space for the low-level potential
    allocate(buffer_fine_pot(grid%ld1, grid%ld2, grid%max_slabs12, &
         pub_num_spins, 2), stat=ierr)
    call utils_alloc_check('xc_fxc_finite_diff_emft','buffer_fine_pot',ierr)

    ! jcap: xc part of the hamiltonian is the only part that
    ! changes between different levels of DFT theory. First,
    ! calculate xc for low level theory for active subregion

    if (present(eff_nhat_den).and.present(nhat_den_grad_gs)) then
       call xc_fxc_finite_diff(density_fine,buffer_fine_pot,grid,cell,&
            ground_state_dens,nhat_dim,eff_nhat_den,nhat_den_grad_gs)
    else
       call xc_fxc_finite_diff(density_fine,buffer_fine_pot,grid,cell,&
            ground_state_dens,nhat_dim)
    end if

    ! jcap: need to swap to high level xc functional, so call
    ! xc_exit and xc_init again, with optional argument to use
    ! active subregion functional (bit of a hack)
    call xc_embed_swap_functional(.true.)

    ! jcap: calculate xc for high level theory for active subregion
    if (present(eff_nhat_den).and.present(nhat_den_grad_gs)) then
       call xc_fxc_finite_diff(density_fine,fxc_pot_fine,grid,cell,&
            ground_state_dens,nhat_dim,eff_nhat_den,nhat_den_grad_gs)
    else
       call xc_fxc_finite_diff(density_fine,fxc_pot_fine,grid,cell,&
            ground_state_dens,nhat_dim)
    end if

    ! jcap: need to swap the xc functional back
    call xc_embed_swap_functional(.false.)

    ! rc2013: subtract the low-level potential from the higher level
    fxc_pot_fine = fxc_pot_fine - buffer_fine_pot

    ! rc2013: clean-up
    deallocate(buffer_fine_pot, stat=ierr)
    call utils_dealloc_check('xc_fxc_finite_diff_emft','buffer_fine_pot',ierr)

  end subroutine xc_fxc_finite_diff_emft


end module xc
