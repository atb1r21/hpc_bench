!> \brief Store of convergence and iteration control parameters:
!!
!!  - defect correction
!!  - conjugate gradient
!!  - multigrid V cycle
!!  - multigirid coarse level ieration
!!  - Newton iteration
!!
!! The governing convergence criteria is as follow:
!!
!! A. if defect correction order is > 2 the following applies:
!!
!!    defect norm < MAX(tol_rel_res * (initial defect norm), tol_res_abs)
!!
!!    AND
!!
!!    ||pot^(n) -pot^(n-1)|| < MAX(tol_res_pot * ||pot^(n-1)||, tol_abs_pot)
!!
!!    These two conditions try to ensure that the iteration process does
!!    not stops too early because of a norm reduction stall or small defect.
!!    The tol_abs_{res,pot} values ensure stopping in the case of zero solution
!!    or account for roundoff errors.
!!
!! B. defect correction == 2
!!
!!    In this case the standard the multigrid convergence algorithm is used
!!
!!    residual norm < MAX(tol_rel_res * || f||, tol_res_abs)
!!
!!   where tol_res_rel, and tol_res_abs are transferrend to
!!   the selected equation (Poisson or or PB flavor) and their tolerance
!!   parameters are ignored if not present in the argument list.
!!   Also, tol_rel_pot and tol_abs_pot are ignored in this cases.
!!
!!  When the PBE equation is solved with periodic boundary conditions the chemical potential convergence
!!  is tested with the following condition:
!!
!!  abs(mu^(i) - mu^(i-1)) < MAX(tol_mu_rel * abs(mu^(i), tol_mu_abs)
!!
!! Convergece is reached if the above condition and the Newton converge test are satisfied

module dl_mg_convergence_params
  use dl_mg_params, only : wp, EQ_PBE_FAS, EQ_PBE_FAS_NOSTERIC, &
       EQ_PBE_NEWTON, EQ_PBE_NEWTON_NOSTERIC
  implicit none

  type conv_params_t
     real(wp) rel
     real(wp) abs
     integer maxiter
  end type conv_params_t

  ! defaults
  real(wp) :: tol_res_rel_            = 1.0e-2_wp ! defect reduction need not be very tight
  real(wp) :: tol_res_abs_            = 1.0e-2_wp
  real(wp) :: tol_pot_rel_            = 1.0e-6_wp
  real(wp) :: tol_pot_abs_            = 1.0e-6_wp
  real(wp) :: tol_newton_res_rel_     = 1.0e-8_wp
  real(wp) :: tol_newton_res_abs_     = 1.0e-8_wp
  real(wp) :: tol_cg_res_rel_         = 1.0e-8_wp
  real(wp) :: tol_cg_res_abs_         = 1.0e-8_wp
  real(wp) :: tol_vcyc_res_rel_       = 1.0e-8_wp
  real(wp) :: tol_vcyc_res_abs_       = 1.0e-8_wp
  real(wp) :: tol_level_1_res_rel_ !< initialisation is done in the subroutine from below
  real(wp) :: tol_level_1_res_abs_    = 1.0e-12_wp
  real(wp) :: tol_mu_rel_             = 1.0e-3_wp
  real(wp) :: tol_mu_abs_             = 1.0e-3_wp

  integer :: max_iters_defco_        =  30 !< maximum number of iterations for the higher order correction
  integer :: max_iters_newton_       =  30 !< .... Newton iteration
  integer :: max_iters_cg_            = 100
  integer :: max_iters_vcyc_         =  30 !< maximum number of V cycles (recommended value not more that 50)
  !! if conjugate gradient is on this is reset to 1, unless set
  !! explicitly via input argument or envinroment variable 
  integer :: max_iters_level_1_      = 100 !< max number of iteration
                                        !! at the coarsest level 

  !> Number of smoother iteration to be applied at each level.
  !! First value is for the descending branch (coarsening),
  !! the second is for the ascending branch (prolongation).
  integer :: v_iterations_(2)        = (/2, 1/)

  real(wp) :: omega_sor_             = 1.15 !< SOR relaxation parameter,
  !! the default is set tp 1.15 as recommended in Trottenberg's book

  !> if consecutive residuals ratio is larger than the value below
  !! it is assumed that
  !! round-off limit was reached. i.e. the computation cannot
  !! progress becase of roundoff errors, or perhaps the problem
  !! is too hard for this algorithm
  real(wp) :: mg_max_conv_rate_ = 0.85_wp

  !! switch for CG, is this right place for it?
  logical, save :: use_cg_ = .false.
  
  !> convergence tags used in subroutine  mg_convergence_test:
  !!  MG coarse level, MG V-cycle, Cojugate Gradient, Newton
  integer, parameter :: tconv_level_1=1, tconv_vcycle=2, &
        tconv_cg=3, tconv_newton=4, nconv_tags=4

  !> used in mg_convergence_test
  type(conv_params_t), save :: conv_params(nconv_tags)


contains

  subroutine dl_mg_set_conv_params(eq_type, fd_order, &
       tol_res_rel, tol_res_abs,                 &
       tol_pot_rel, tol_pot_abs,                 &
       tol_newton_res_rel, tol_newton_res_abs,   &
       tol_cg_res_rel, tol_cg_res_abs,           &
       tol_vcyc_res_rel, tol_vcyc_res_abs,       &
       tol_level_1_res_rel, tol_level_1_res_abs, &
       tol_mu_rel,          tol_mu_abs,          &
       mg_max_conv_rate, &
       max_iters_defco, &
       max_iters_newton, &
       max_iters_cg,&
       max_iters_vcyc, &
       max_iters_level_1, &
       v_iterations, &
       omega_sor, &
       use_cg)
    use dl_mg_params, only : DL_MG_ERR_UNSPECIFIED
    use dl_mg_errors, only : handle_error, check_assertion
    use dl_mg_common_data, only : mg_levels
    implicit none

    integer,  intent(in)            :: eq_type
    integer,  intent(in)            :: fd_order
    real(wp), intent(in),  optional :: tol_res_rel
    real(wp), intent(in),  optional :: tol_res_abs
    real(wp), intent(in),  optional :: tol_pot_rel
    real(wp), intent(in),  optional :: tol_pot_abs
    real(wp), intent(in),  optional :: tol_newton_res_rel
    real(wp), intent(in),  optional :: tol_newton_res_abs
    real(wp), intent(in),  optional :: tol_cg_res_rel
    real(wp), intent(in),  optional :: tol_cg_res_abs
    real(wp), intent(in),  optional :: tol_vcyc_res_rel
    real(wp), intent(in),  optional :: tol_vcyc_res_abs
    real(wp), intent(in),  optional :: tol_level_1_res_rel
    real(wp), intent(in),  optional :: tol_level_1_res_abs
    real(wp), intent(in), optional :: tol_mu_rel
    real(wp), intent(in), optional :: tol_mu_abs
    real(wp), intent(in),  optional :: mg_max_conv_rate
    integer,  intent(in),  optional :: max_iters_defco
    integer,  intent(in),  optional :: max_iters_newton
    integer,  intent(in),  optional :: max_iters_cg
    integer,  intent(in),  optional :: max_iters_vcyc
    integer,  intent(in),  optional :: max_iters_level_1
    integer,  intent(in),  optional :: v_iterations(2)
    real(wp), intent(in),  optional :: omega_sor
    logical, intent(in), optional   :: use_cg 

    integer istat, len, pos
    character(len=32) val

    call get_environment_variable("DL_MG_USE_CG",value=val, length=len, status=istat)
    if ( len > 0) then
       select case (val(1:1))
          case ('T','t','F','f')
             read(val(1:1),*) use_cg_
          case default
             call handle_error(DL_MG_ERR_UNSPECIFIED, &
                  msg="The value of DL_MG_USE_CG must be one of the following: T,t,F,f")
          end select
    else
       if (present(use_cg)) use_cg_ = use_cg
    end if

    if (present(tol_res_rel)) then
       if (fd_order > 2) then
          tol_res_rel_ = tol_res_rel
       else
          ! the top tolerance is atributed to the top iteration
          ! but this it is superseded if the specific variable is present,
          ! see below
          if (eq_type == EQ_PBE_NEWTON .or. &
             eq_type == EQ_PBE_NEWTON_NOSTERIC) then
             tol_newton_res_rel_ =  tol_res_rel
          else
             if (use_cg) then
                tol_cg_res_rel_ = tol_res_rel
             else
                tol_vcyc_res_rel_ = tol_res_rel
             endif
          endif
       endif
    endif

    if (present(tol_res_abs)) then
       if (fd_order > 2) then
          tol_res_abs_ = tol_res_abs
       else
          ! the top tolerance is atributed to the top iteration
          if (eq_type == EQ_PBE_NEWTON .or. &
               eq_type == EQ_PBE_NEWTON_NOSTERIC) then
             ! do this only if the rel conterpart is not
             ! present, otherwise keep the default
             ! this is to avoid setting abs tol to the generic valure
             ! when only tol rel is provided
             ! for newton or v cycle
             if (.not. present(tol_newton_res_rel)) then
                tol_newton_res_abs_ =  tol_res_abs
             end if
          else
             ! the same as above
             if (use_cg) then
                if(.not. present(tol_cg_res_rel)) then
                   tol_cg_res_abs_ = tol_res_abs
                end if
             else
                if(.not. present(tol_vcyc_res_rel)) then
                   tol_vcyc_res_abs_ = tol_res_abs
                end if
             end if
          endif
       endif
    endif

    !> set more convergence parameters if provided
    if (present(tol_pot_rel)) then
       tol_pot_rel_ = tol_pot_rel
    endif

    if (present(tol_pot_abs)) then
       tol_pot_abs_ = tol_pot_abs
    endif

    if (present(tol_newton_res_rel)) then
       tol_newton_res_rel_ = tol_newton_res_rel
    endif

    if (present(tol_newton_res_abs)) then
       tol_newton_res_abs_ = tol_newton_res_abs
    endif

    if (present(tol_cg_res_rel)) then
       tol_cg_res_rel_ = tol_cg_res_rel
    endif

    if (present(tol_cg_res_abs)) then
       tol_cg_res_abs_ = tol_cg_res_abs
    endif
    
    if (present(tol_vcyc_res_rel)) then
       tol_vcyc_res_rel_ = tol_vcyc_res_rel
    endif

    if (present(tol_vcyc_res_abs)) then
       tol_vcyc_res_abs_ = tol_vcyc_res_abs
    endif

    if (present(tol_level_1_res_rel)) then
       tol_level_1_res_rel_ = tol_level_1_res_rel
    else
       ! Experimentatally I found that a higher precision is needed at the coarse level for FAS.
       ! Set the FAS value here only if it was not set before by the corresponding environment variable
       if ( eq_type == EQ_PBE_FAS .or. eq_type == EQ_PBE_FAS_NOSTERIC) then
          tol_level_1_res_rel_ = 0.1_wp * tol_vcyc_res_rel_
       else
          if (use_cg_) then 
            tol_level_1_res_rel_ = 0.1_wp ! from other experiments
          else
            tol_level_1_res_rel_ = 10.0_wp * tol_vcyc_res_rel_ ! seems to be enough in the studied cases
          endif
       end if
    endif

    if (present(tol_mu_rel)) then
       tol_mu_rel_ = tol_mu_rel
    endif

    if (present(tol_mu_abs)) then
       tol_mu_abs_ = tol_mu_abs
    endif


    if (present(tol_level_1_res_abs)) then
       tol_level_1_res_abs_ = tol_level_1_res_abs
    endif

    if (present(max_iters_defco))   max_iters_defco_   = max_iters_defco
    if (present(max_iters_newton))  max_iters_newton_  = max_iters_newton
    if (present(max_iters_cg))      max_iters_cg_    = max_iters_cg

    if ( use_cg_) then
       max_iters_vcyc_  = 1
    endif
    call get_environment_variable("DL_MG_NITER_VCYC",value=val, length=len, status=istat)
    if ( len > 0) then
        read(val,*) max_iters_vcyc_
    else
       if (present(max_iters_vcyc)) max_iters_vcyc_ = max_iters_vcyc
    endif

    ! the setting via environment variable are  mainly for debugging
    call get_environment_variable("DL_MG_NITER_LEVEL1",value=val, length=len, status=istat)
    if ( len > 0) then
       read(val,*) max_iters_level_1_
    else
       if (present(max_iters_level_1)) max_iters_level_1_ = max_iters_level_1
    endif
    
    ! now we can populate the in conv_params array
    conv_params(tconv_level_1)%abs     = tol_level_1_res_abs_
    conv_params(tconv_level_1)%rel     = tol_level_1_res_rel_
    conv_params(tconv_level_1)%maxiter = max_iters_level_1_
    conv_params(tconv_vcycle)%abs      = tol_vcyc_res_abs_
    conv_params(tconv_vcycle)%rel      = tol_vcyc_res_rel_
    conv_params(tconv_vcycle)%maxiter  = max_iters_vcyc_
    conv_params(tconv_cg)%abs          = tol_cg_res_abs_
    conv_params(tconv_cg)%rel          = tol_cg_res_rel_
    conv_params(tconv_cg)%maxiter      = max_iters_cg_
    conv_params(tconv_newton)%abs      = tol_newton_res_abs_
    conv_params(tconv_newton)%rel      = tol_newton_res_rel_
    conv_params(tconv_newton)%maxiter  = max_iters_newton_

    if (present(mg_max_conv_rate)) then
       mg_max_conv_rate_ = mg_max_conv_rate
    endif

    call get_environment_variable("DL_MG_V_ITER",value=val, length=len, status=istat)
    if ( len > 0) then
       pos = index(val,",")
       read(val(1:pos-1),*) v_iterations_(1)
       read(val(pos+1:),*) v_iterations_(2)
    else
       if (present(v_iterations)) v_iterations_ = v_iterations
    end if

    call get_environment_variable("DL_MG_OMEGA_SOR",value=val, length=len, status=istat)
    !if (istat /= 0 ) then  to be continued with warning
    if ( len > 0) then
       read(val,*) omega_sor_
    else
       if (present(omega_sor)) omega_sor_ = omega_sor
    end if

    if ( mg_levels == 1 .and. use_cg) then
       ! CG doesn't converge with red-black preconditioner with
       ! omega_sor != 1
       omega_sor_ = 1.0_wp
    endif
    !if (check_assertion( use_cg_ .and. omega_sor_ /= 1.0_wp)) then
    !   call handle_error(DL_MG_ERR_UNSPECIFIED, msg="dl_mg_convergence_params:&
    !        & Conjugate gradient works only with omega_sor=1.0")
    !end if
    
    call write_convergence_params(eq_type, fd_order, &
       tol_res_rel_, tol_res_abs_,                 &
       tol_pot_rel_, tol_pot_abs_,                 &
       tol_newton_res_rel_, tol_newton_res_abs_,   &
       tol_cg_res_rel_, tol_cg_res_abs_,           &
       tol_vcyc_res_rel_, tol_vcyc_res_abs_,       &
       tol_level_1_res_rel_, tol_level_1_res_abs_, &
       tol_mu_rel_,           tol_mu_abs_,         &
       mg_max_conv_rate_, &
       max_iters_defco_, &
       max_iters_newton_, &
       max_iters_cg_, &
       max_iters_vcyc_, &
       max_iters_level_1_, &
       v_iterations_, &
       omega_sor_, &
       use_cg_)

  end subroutine dl_mg_set_conv_params

   !>\todo print only the params used in present calculation.
  !! i.e. don't print info on Newton for Poisson case
  subroutine write_convergence_params(eq_type, fd_order, &
       tol_res_rel, tol_res_abs, &
       tol_pot_rel, tol_pot_abs, &
       tol_newton_res_rel, tol_newton_res_abs,&
       tol_cg_res_rel, tol_cg_res_abs,&
       tol_vcyc_res_rel, tol_vcyc_res_abs,&
       tol_level_1_res_rel, tol_level_1_res_abs, &
       tol_mu_rel,          tol_mu_abs,          &
       mg_max_conv_rate, &
       max_iters_defco, &
       max_iters_newton, &
       max_iters_cg, &
       max_iters_vcyc, &
       max_iters_level_1, &
       v_iterations, &
       omega_sor, &
       use_cg)

    use dl_mg_mpi_header
    use dl_mg_params, only : EQ_PBE_NEWTON, EQ_PBE_NEWTON_NOSTERIC, &
         EQ_POISSON, DL_MG_BC_PERIODIC
    use dl_mg_common_data, only : report_unit, report_file, bc, mg_levels
    use dl_mg_mpi, only : get_myid
    use dl_mg_errors
    implicit none

    integer,  intent(in) :: eq_type
    integer,  intent(in) :: fd_order
    real(wp), intent(in) :: tol_res_rel
    real(wp), intent(in) :: tol_res_abs
    real(wp), intent(in) :: tol_pot_rel
    real(wp), intent(in) :: tol_pot_abs
    real(wp), intent(in) :: tol_newton_res_rel
    real(wp), intent(in) :: tol_newton_res_abs
    real(wp), intent(in) :: tol_cg_res_rel
    real(wp), intent(in) :: tol_cg_res_abs
    real(wp), intent(in) :: tol_vcyc_res_rel
    real(wp), intent(in) :: tol_vcyc_res_abs
    real(wp), intent(in) :: tol_level_1_res_rel
    real(wp), intent(in) :: tol_level_1_res_abs
    real(wp), intent(in) :: tol_mu_rel
    real(wp), intent(in) :: tol_mu_abs
    real(wp), intent(in) :: mg_max_conv_rate
    integer,  intent(in) :: max_iters_defco
    integer,  intent(in) :: max_iters_newton
    integer,  intent(in) :: max_iters_cg
    integer,  intent(in) :: max_iters_vcyc
    integer,  intent(in) :: max_iters_level_1
    integer,  intent(in) :: v_iterations(2)
    real(wp), intent(in) :: omega_sor
    logical, intent(in)  :: use_cg


   integer myid, i
   logical lop


   myid = get_myid()

   if (myid == 0) then
      inquire (unit=report_unit,opened=lop)
      if (.not. lop) then
         open (unit=report_unit,file=report_file,status="old",position="append")
      endif

      ! the PBE details are reported in write_nonlinear_params
      if (eq_type == EQ_POISSON) then
         write(report_unit, '(/,a,/)') "Solving Poisson Equation"
      end if

      if (use_cg) then
         write(report_unit, '(a)') "Conjugate gradient enabled"
         if (mg_levels == 1) then
            write(report_unit,'(a/)') "  Omega_sor = 1 for the red-black preconditioner"
         else
            write(report_unit,*)
         endif
      end if

      write(report_unit,'(a,1x,I3)') 'Finite difference approximation order:', fd_order

      write(report_unit,'(/,a)') 'convergence tols     rel     abs'
      write(report_unit,'("-------------------------------------------------------")')
      if ( fd_order > 2) then
         write(report_unit, 100) 'defco residual ',  tol_res_rel,  tol_res_abs
         write(report_unit, 100) 'defco potential',  tol_pot_rel,  tol_pot_abs
      end if
      if (eq_type == EQ_PBE_NEWTON .or. &
           eq_type == EQ_PBE_NEWTON_NOSTERIC) then
         write(report_unit, 100) 'Newton         ',  tol_newton_res_rel, tol_newton_res_abs
         if (all(bc == DL_MG_BC_PERIODIC)) then
            write(report_unit, 100) '  mu (PBC)     ', tol_mu_rel, tol_mu_abs
         end if
      end if
      if (use_cg) then
         write(report_unit, 100) 'CG             ',  tol_cg_res_rel, tol_cg_res_abs
      end if
      write(report_unit, 100) 'V cycle        ',  tol_vcyc_res_rel, tol_vcyc_res_abs
      write(report_unit, 100) 'Level 1        ',  tol_level_1_res_rel,tol_level_1_res_abs
100   format(a15, 2E9.2)


      write(report_unit,'(/,"Max number of iterations:")')
      write(report_unit,'("-------------------------------------------------------")')
      if ( fd_order > 2) then
         write(report_unit, 200) "defect correction: ", max_iters_defco
      end if
      if (eq_type == EQ_PBE_NEWTON .or. &
           eq_type == EQ_PBE_NEWTON_NOSTERIC) then
         write(report_unit, 200) "Newton:            ", max_iters_newton
      end if
      if (use_cg) then
         write(report_unit, 200) "CG:                ", max_iters_cg
      end if
      write(report_unit, 200) "Vcycle:            ", max_iters_vcyc
      write(report_unit, 200) "MG level 1:        ", max_iters_level_1
200   format (a,1x,I6)
      write(report_unit, '(/,a,1x,2(I4,1x))') "V cycle smoother iterations :", &
           v_iterations
      write(report_unit, '(a,1x,g10.3)')    "multigrid max convergence rate:", &
           mg_max_conv_rate
      write(report_unit, '(a,1x,E10.3,/)')    "Omega SOR                   :", &
           omega_sor

       close(report_unit)
       open (unit=report_unit,file=report_file,status="old", position="append")
    end if

  end subroutine write_convergence_params

  
end module dl_mg_convergence_params
