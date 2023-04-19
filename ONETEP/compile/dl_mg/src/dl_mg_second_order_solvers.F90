module dl_mg_second_order_solvers
  implicit none
  private

!  logical, save :: use_cg
  
  public dlmg_2nd_order_solver 
contains

  !> main driver for the second order solver, parameters are documented in the wrappers: 
  !> dl_mg::dl_mg_solver_poisson(), dl_mg::dl_mg_solver_pbe()
  subroutine dlmg_2nd_order_solver(eq_type, alpha, rho, &
       pot, use_pot_in, remove_source_zero_mode, res, der_pot, ierror)
    ! main solver subroutine
    ! finds the solution of the discretised version of Poisson Equation
    !
    ! \nabla (c(x,y,z) \grad p(x,y,z)) = alpha * f(x,y,z) + fnl(p(x,y,z))
    !
    ! for electrostatic problems alpha = -4 * pi, or -1 according to the units system
    !
    ! fnl(p) = lambda  steric_weight  [\sum c_i q_i exp(-beta ( q_i  p ) )]
    !
    ! JCW, 30/01/17:
    ! If optional argument steric_weight is not present present it is assummed
    ! to be 1.0 for all space.
    use dl_mg_mpi_header
    use dl_mg_types
    use dl_mg_common_data
    use dl_mg_convergence_params, only : use_cg => use_cg_
    use dl_mg_info, only  : mg_info_write
    use dl_mg_grids, only : mg
    use dl_mg_timer
    use dl_mg_errors
    use dl_mg_params, only : DL_MG_SUCCESS
    use dl_mg_multigrid_method, only : mg_solver
    use dl_mg_conjugate_gradient_method, only : conj_grad_solver
    use dl_mg_newton_method, only : newton_solver
    use dl_mg_utils, only : set_mass_term
    ! debug
    use dl_mg_kernels
    use dl_mg_alloc
    implicit none

    integer, intent(in)                 :: eq_type
    real(wp), intent(in)                :: alpha, rho(isx:, isy:, isz:)
    real(wp), intent(inout)             :: pot(isx:, isy:, isz:)
    real(wp), optional, intent(out)     :: res(isx:, isy:, isz:)
    real(wp), optional, intent(in)      :: der_pot(isx:, isy:, isz:) !< value
    !! of potential at with the derivative of Boltzmann term is computed

    logical, optional, intent(in)       :: use_pot_in !< if absent or false
    !! the inital solution is set 0 inside domain
    integer, optional, intent(in)   :: remove_source_zero_mode
    integer, optional, intent(out)      :: ierror

    integer ierr, ierrt
    !debug
    integer t


    call mg_timer(START, TMG, TIGNORE, TIGNORE)

    if (present(ierror)) then
       ierror = DL_MG_SUCCESS
    endif

    !$OMP PARALLEL DEFAULT(NONE) SHARED(mg, eq_type, &
    !$OMP& alpha, rho, pot, res, der_pot, &
    !$OMP& use_cg, use_pot_in, remove_source_zero_mode, ierr,t) &
    !$OMP& FIRSTPRIVATE(ierrt)

    if ( eq_type == EQ_PBE_NEWTON .or. &
         eq_type == EQ_PBE_NEWTON_NOSTERIC ) then
       call newton_solver(mg, eq_type, alpha, rho, pot, res, use_pot_in, ierrt)
    else
       if (eq_type == EQ_LINEAR_PBE_POT) then
            call set_mass_term(mg, der_pot, ierr=ierr)
         end if
       if (use_cg) then
          ! needs to set res to 0 on boundaries as mg_solver does
          call conj_grad_solver(mg, eq_type, alpha, rho, pot, res,&
               use_pot_in, ierr=ierrt)
       else
          call mg_solver(mg, eq_type, alpha, rho, pot, res, &
            use_pot_in, ierr=ierrt)
       endif
    endif

    !$OMP MASTER
    ierr = ierrt
    !$OMP END MASTER

    !$OMP END PARALLEL

    call mg_info_write(eq_type)
    
    if ( .not. check_assertion(ierr == DL_MG_SUCCESS)) then
       call handle_error(ierr, ierror)
    endif

    call mg_timer(STOP, TMG, TIGNORE, TIGNORE)

  end subroutine dlmg_2nd_order_solver
  
end module dl_mg_second_order_solvers
