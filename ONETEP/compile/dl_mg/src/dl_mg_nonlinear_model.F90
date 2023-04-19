!> \brief Computes the Boltzmann term and its derivative
!!
!! uses capped exponential to avoid overflow
!!
!! Lucian Anton
!!
!! \todo Warn the user if the provided expcap is overruled by max_excap; should max_expcap be located in parameters module?

module dl_mg_nonlinear_model
  use dl_mg_params, only : wp
  implicit none

  integer, save :: nion !< number of ion types
  real (wp), save :: lambda, temp, kappa2, expcap, beta
  ! kappa2 is the  the used for linerized version of PBE
  ! kappa2 = -lambda  beta  \sum c_i q_i^2

  real (wp), allocatable, save :: qion(:)  !< ions charges;  array of dimension (nion)
  real (wp), allocatable, save :: cion_in(:) !< the input concentrations
  real(wp), allocatable, save :: cion(:) !< used concentrations, they might differ from cion_in 
                                         !< if neutralisation by concentration shift is used

  real (wp), allocatable, save :: mu(:) !< chemical potential * beta
  !! computed only if full PBE and full PBC
  real (wp), allocatable, save :: mu_prev(:) !< previous chem pot, needed in convergence test; this is ugly !

  real(wp), parameter :: max_expcap = 50.0_wp

  logical :: has_steric_weight = .false. !< true if steric_weight array is provided

  !> steric weight sum over the whole domain ( needed when using  PBC)
  real(wp), save :: steric_weight_sum
  real(wp), save :: rho_sum = huge(0.0_wp)   !< sum charge density
  !! over internal grid points.
  !! Needed in some neutralisations

  !> neutralisation method
  integer, save :: neutralisation_method


  contains

    function fnl(x, steric_w)
      implicit none
      real(wp), intent(in) :: x, steric_w
      real(wp) fnl

      integer  i
      real(wp) y, z


      if (lambda == 0.0_wp ) then
         fnl = 0.0_wp
      else
         fnl = 0.0_wp
         z = beta * x
         do i = 1, nion
            y = - qion(i) * z + mu(i)
            fnl = fnl + cion(i) * qion(i) * steric_w * capexp(y)
         enddo
         fnl = lambda * fnl
      endif
    end function fnl


    function derfnl(x, steric_w)
      implicit none
      real(wp), intent(in) :: x, steric_w
      real(wp) derfnl

      integer   i
      real(wp)  y,z

      if ( lambda == 0.0_wp) then
         derfnl = 0.0_wp
      else
         derfnl = 0.0_wp
         z = beta * x
         do i = 1, nion
            y = - qion(i) * z + mu(i)
            if ( y < expcap ) then
               derfnl = derfnl - cion(i) * qion(i)**2  * beta * &
                    steric_w * exp(y)
            endif
         enddo
         derfnl = lambda * derfnl
      end if

    end function derfnl


    elemental function capexp(x)
      implicit none
      real(wp), intent(in) ::x
      real(wp) capexp

      capexp = exp(min(x,expcap))
    end function capexp

    !> initialise non-linear module data
    subroutine init_nonlinear(lineq, lam, t, n, q, c, &
         has_stw, expcap_in, ntrl_method, ntrl_ion_ratios)
      use dl_mg_params, only : wp, EQ_LINEAR_PBE, EQ_LINEAR_PBE_NOSTERIC, &
           dl_mg_neutralise_with_jellium_vacc, &
           dl_mg_neutralise_none, &
           hartree
      use dl_mg_common_data, only : fullpbc
      use dl_mg_errors
      implicit none
      integer, intent(in)  :: n
      logical, intent(in) :: lineq
      real(wp), intent(in) :: lam, t, q(n), c(n)
      logical, intent(in) :: has_stw
      real(wp), optional, intent(in) :: expcap_in
      integer, optional, intent(in)  :: ntrl_method
      real(wp), optional, intent(in) ::ntrl_ion_ratios(n)


      temp = t
      nion = n
      lambda =lam
      beta = Hartree / t
      allocate(qion(nion), cion_in(nion), cion(nion))
      qion = q
      cion_in = c
      allocate(mu(nion), mu_prev(nion))
      mu      = 0.0_wp
      mu_prev = 0.0_wp
      has_steric_weight = has_stw

      if(present(expcap_in)) then
         expcap = min(expcap_in, max_expcap)
      else
         expcap = max_expcap
      endif

      ! if c_bulk == 0 only  linearised PBE has a well defined
      ! potential shift
      if ( .not. lineq .and. all(c == 0.0_wp)) then
         call handle_error(DL_MG_ERR_CION_ZERO)
      end if

      if (present(ntrl_method)) then
         neutralisation_method = ntrl_method
      else
       ! set the "natural value" neutralisation
       ! on BC
       if (fullpbc) then
          neutralisation_method = dl_mg_neutralise_with_jellium_vacc
       else
          neutralisation_method = dl_mg_neutralise_none
       end if
    end if

  end subroutine init_nonlinear

    !> set the inverse square of Debye length with the used
    ! concentration
    subroutine set_debye_length(eq_type)
      use dl_mg_params, only : EQ_LINEAR_PBE, EQ_LINEAR_PBE_NOSTERIC
      use dl_mg_errors
      implicit none

      integer, intent(in) :: eq_type

      select case(eq_type)
      case(EQ_LINEAR_PBE, EQ_LINEAR_PBE_NOSTERIC)
         kappa2 = -lambda * beta * sum ( cion(1:nion) &
              * qion(1:nion)**2 )
         if ( .not. check_assertion(kappa2 >= 0.0_wp, DL_MG_ERR_DEBYE)) then
            call handle_error(DL_MG_ERR_DEBYE)
         end if
      case default
         ! the code should never use this value
         kappa2 = -HUGE(kappa2)
      end select

    end subroutine set_debye_length


    subroutine free_nonlinear
      implicit none

      if (allocated(qion))     deallocate(qion)
      if (allocated(cion_in))  deallocate(cion_in)
      if (allocated(cion))     deallocate(cion)
      if (allocated(mu))       deallocate(mu)
      if (allocated(mu_prev))  deallocate(mu_prev)

    end subroutine free_nonlinear

end module dl_mg_nonlinear_model
