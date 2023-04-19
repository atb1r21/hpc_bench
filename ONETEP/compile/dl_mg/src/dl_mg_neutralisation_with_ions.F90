!> this is essentially a submodule of dl_mg_nonlinear_model
!! which handles the computation of neutralisation ratios
module dl_mg_neutralisation_with_ions
  use dl_mg_params, only : wp
  use dl_mg_nonlinear_model, only : nion, q => qion, &
       c => cion, cin => cion_in, mu, &
       ntrl_m => neutralisation_method, &
       steric_weight_sum, rho_sum

  implicit none

  real(wp), private, parameter :: chrg_tol=1.e-14 ! tolerance below which
                                         ! two charges are considered equal

  logical, private, save :: use_linear_approx = .false.
  logical, private, save :: use_qint = .false. ! use integer q in powers, it should
                                      ! be faster
  logical, private, save :: debug_print = .false.

  !> the algoritm does not need this array to be shared but
  !! it is useful for returning the ratios at the end of calculation
  real(wp), allocatable, save :: counterion_ratios(:)
  !> integer ion charges set up when use_qint is true
  integer, allocatable, private, save :: qi(:)

  public  :: init_counterion_neutralisation, compute_shifted_ion_concentrations, &
       free_neutralisation
  private :: func, der_func

contains

  subroutine init_counterion_neutralisation( ntrl_ratios, ierr)
    use dl_mg_errors
    use dl_mg_params
    implicit none

    real(wp), optional, intent(in) :: ntrl_ratios(:)
    integer, optional, intent(inout) :: ierr

    integer i
    character(len=1) :: val, caux*12

    ! check if we have any bussines to do here
    if (ntrl_m == dl_mg_neutralise_with_jellium_uniform &
         .or. &
         ntrl_m == dl_mg_neutralise_with_jellium_vacc &
         .or. &
         ntrl_m == dl_mg_neutralise_none) then
       c = cin
       return
    end if

    ! collect possible environment variables settings for
    ! debug writing
    call get_environment_variable("DL_MG_NEUTRALISE_DBG_PRINT",value=val)

    if (len_trim(val) > 0) then
       select case(val)
       case ('0','F','f')
          debug_print = .false.
       case default
          debug_print = .true.
       end select
    end if

    ! use integer charges if possible
    if (all(abs(q -nint(q)) < chrg_tol)) then
       use_qint = .true.
       allocate(qi(nion))
       qi = nint(q)
    end if

    allocate(counterion_ratios(nion))

    select case(ntrl_m)
    case(dl_mg_neutralise_with_ions_fixed)
       counterion_ratios &
            = ntrl_ratios/sum(ntrl_ratios)

       c = cin - counterion_ratios &
            * rho_sum/(q * steric_weight_sum)

       do i = 1, nion
          if (check_assertion(c(i) <= 0.0_wp,&
               DL_MG_ERR_CION_NEG)) then
             write(caux, '(i10)') i
             call handle_error(DL_MG_ERR_CION_NEG, ierr, msg="shifted c" &
                 //trim(adjustl(caux))//" < 0")
             exit
          end if
       end do

    case(dl_mg_neutralise_with_ions_auto)
       call compute_shifted_ion_concentrations
    case(dl_mg_neutralise_with_ions_auto_linear)
       use_linear_approx = .true.
       call compute_shifted_ion_concentrations
    end select

  end subroutine init_counterion_neutralisation


  subroutine compute_shifted_ion_concentrations(ierr)
    use dl_mg_params
    use dl_mg_errors
    implicit none

    integer, optional, intent(inout) :: ierr

    integer i
    real(wp) cs
    character(len=100) :: caux

    ! do we have any business here ?
    if (.not. (ntrl_m == dl_mg_neutralise_with_ions_auto &
         .or. &
         ntrl_m == dl_mg_neutralise_with_ions_auto_linear) ) then
       return
    end if

    cs = rho_sum/steric_weight_sum

    if (use_linear_approx) then
       call x_by_lin_approx(cs)
    else

       if (nion == 2 .and. abs(q(1) + q(2)) < chrg_tol) then
          call x_by_quadratic(cs)
       else
          call x_by_newton(cs)
       end if
    end if


    c = cin - counterion_ratios * cs / q
    do i = 1, nion
       if (c(i) < 0.0_wp ) then
          write(caux,'(a,i6,1x,e13.5)') ' x fraction', i, counterion_ratios(i)
          call handle_error(DL_MG_ERR_CION_NEG, ierr, &
               msg='c(i) < 0 in compute_counter_ion_fractions'//trim(caux))
          exit
       endif
    end do

  contains

    subroutine x_by_quadratic(cs)
      implicit none

      real(wp), intent(in) :: cs

      real(wp) expmu, a1, a0, s

      expmu = exp(-mu(1)-mu(2)) !exp((cqb(2,1) * mu(2) - cqb(2,2) * mu(1))/cqb(2,2))

    ! the equation to solve is
    !
    !  (1 - x*cs/(q1*b1) ) * ( 1 - (1-x)*cs/(q2*b2) ) = expmu
    !  or
    !  (x - (q1*b1)/cs) * ( x + q2 * b2 / cs -1) = - expmu * c1 * b1 * q2 * b2 /cs**2

    a1 = q(2) * cin(2)/cs - 1.0_wp - q(1) * cin(1)/ cs
    a0 = expmu * q(1) * cin(1) * q(2) * cin(2) / cs**2 &
         - q(1) * cin(1)/ cs * (q(2) * cin(2)/cs -1.0_wp)

    if ( a1**2 - 4.0_wp * a0 < 0) then
       call handle_error(DL_MG_ERR_UNSPECIFIED, &
            msg="delta < 0 in dl_mg_counter_error")
    else
       ! pick up the solution that has a finite limit when cs -> 0
       if ( a1 < 0.0_wp ) then
          s = -1.0_wp
       else
          s = 1.0_wp
       end if
       counterion_ratios(1) = 2.0_wp * a0 &
            / ( - a1 - s * sqrt (a1**2 - 4.0_wp * a0))
       ! numerically stable version of 0.5_wp * ( - a1 + s * sqrt (a1**2 - 4.0_wp * a0))
       counterion_ratios(2) = 1.0_wp - counterion_ratios(1)
    end if

    if (debug_print) then
       write(*,*)'neutralise x from quadratic', &
            0.5_wp * ( - a1 + sqrt (a1**2 - 4.0_wp * a0)), &
            0.5_wp * ( - a1 - sqrt (a1**2 - 4.0_wp * a0))
    end if
    end subroutine x_by_quadratic


    subroutine x_by_newton(cs)
      implicit none

      real(wp), intent(in) :: cs

      real(wp), parameter :: tol_res = 1.e-3
      real(wp), parameter :: tol_err_rel = 1.e-3
      real(wp), parameter :: tol_err_abs = 1.e-3

      integer i, imax
      real(wp) x, error, res, jac

      imax = 30
      i=1
      ! start formn the linear solution
      x = exp ((sum(mu * q * cin) - cs) &
           / sum(q**2 * cin))

      do
         res = func(x, cs)
         jac = der_func(x)
         !write(*,*) 'nwt cntios', n, res, jac, x, error
         if (abs(jac) < 1.e-8) then
            call handle_error(DL_MG_ERR_UNSPECIFIED, &
                 msg="jacobian less that 1.e-5 in compute_counter_ion_fractions")
         end if
         error = - res / jac
         x = x + error

         if ( abs(res) < tol_res .and. abs(error) < max(tol_err_rel * x, tol_err_abs)) then
            exit
         end if

         i = i+1
         if (i > imax) then
            call handle_error(DL_MG_ERR_UNSPECIFIED, &
                 msg='Newton method did not converge in compute_counter_ion_fractions')
         end if
      end do

      if (use_qint) then
         counterion_ratios = q * cin * (1.0_wp - exp(-mu) * x**qi) /cs
      else
         counterion_ratios = q * cin * (1.0_wp - exp(-mu) * x**q) /cs
      end if

    end subroutine x_by_newton


    subroutine x_by_lin_approx(cs)
      implicit none
      real(wp), intent (in) :: cs

       real(wp) b(nion)

       counterion_ratios = q**2 * cin / sum(q**2 * cin) &
            * (1.0_wp - sum(mu * q * cin)/cs)         &
            + mu * q * cin / cs

    end subroutine x_by_lin_approx

  end subroutine compute_shifted_ion_concentrations


  Pure function func(x, cs)
    use dl_mg_params, only : wp
    implicit none

    real(wp), intent(in) :: x, cs
    real(wp) func

    if (use_qint) then
       func = sum( qi * c * exp(-mu) * x**qi) + cs
    else
       func = sum( q  * c * exp(-mu) * x**q ) + cs
    end if

  end function func


  pure function der_func(x)
    use dl_mg_params, only : wp
    implicit none

    real(wp), intent(in) :: x
    real(wp) der_func

     if (use_qint) then
        der_func = sum( qi**2 * c * exp(-mu)  * x**(qi-1) )
     else
        der_func = sum( q**2  * c * exp(-mu)  * x**(q-1)  )
     endif
  end function der_func


  !> free module arrays
  subroutine free_neutralisation
    implicit none

    if (allocated(counterion_ratios)) deallocate(counterion_ratios)
    if (allocated(qi)) deallocate(qi)

  end subroutine free_neutralisation

end module dl_mg_neutralisation_with_ions
