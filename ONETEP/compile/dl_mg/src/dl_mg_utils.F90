!> auxiliary subroutines
module dl_mg_utils
  implicit none

contains

  !> Initialises the non-linear parameters, called by pbe_params%lam.
  !! It can be called multiple times, e.g. if beta needs to be changed in continuation calculation.
  !!
  !! This subroutine needes to be called also in the case of linearised PBE.
  !!
  !! Note: from physical point of view it is sensible to assume that
  !! the solvent is electrically neutral, i.e. \f$ \sum_i c_i q_i = 0 \f$,
  !! however the neutrality condition is not checked.
  !! \todo Test for charge neutrality?
  subroutine dl_mg_init_nonlin(eq_type, lambda, temp, n, c, q, &
       rho, steric_weight, der_pot, expcap, &
       linearised, use_fas, &
       neutralisation_method,     &
       neutralisation_ion_ratios, &
        ierror)
    use dl_mg_params
    use dl_mg_common_data, only : fullpbc
    use dl_mg_nonlinear_model, only : init_nonlinear, set_debye_length
    use dl_mg_neutralisation_with_ions, only : init_counterion_neutralisation
    use dl_mg_info, only :  write_nonlinear_params, &
         charge_neutralisation_info
    use dl_mg_errors
    implicit none

    integer, intent(out) :: eq_type
    integer, intent(in)  :: n
    real(wp), intent(in) :: lambda, temp, c(n), q(n)
    real(wp), intent(in) :: rho(:,:,:)
    real(wp), intent(in), optional :: steric_weight(:,:,:), der_pot(:,:,:)
    real(wp), intent(in), optional :: expcap
    logical, intent(in),  optional :: linearised, use_fas
    integer,  intent(in), optional  :: neutralisation_method
    real(wp), intent(in), optional  :: neutralisation_ion_ratios(n)
    integer, optional, intent(out)    :: ierror
    !> for arguments description see \ref dl_mg_solver_pbe

    logical lineq, ufas, has_steric_weight

    ierror = DL_MG_SUCCESS

    ! sets has_steric_weight in to be passed to dl_mg_nonlinear_model below
    call compute_steric_weight_sum(has_steric_weight, steric_weight)
    ! sets rho_sum in dl_mg_nonlinear_model
    call compute_rho_sum(rho)

    ! compute eq_type
    lineq = .false.
    if (present(linearised)) lineq = linearised

    ufas = .false.
    if ( present(use_fas)) ufas = use_fas
    if (lineq) then
       if (has_steric_weight) then
          eq_type = EQ_LINEAR_PBE
       else
          eq_type = EQ_LINEAR_PBE_NOSTERIC
       endif

       ! linearise around der_pot
       ! applies for steric and no-steric cases
       if (present(der_pot)) then
          eq_type = EQ_LINEAR_PBE_POT
       endif

    else
       if (has_steric_weight) then
          if (ufas) then
             eq_type = EQ_PBE_FAS
          else
             eq_type = EQ_PBE_NEWTON
          endif
       else
          if (ufas) then
             eq_type = EQ_PBE_FAS_NOSTERIC
          else
             eq_type = EQ_PBE_NEWTON_NOSTERIC
          endif
       endif
    end if

    if (present(neutralisation_method)) then
       ! check for consistency with boundary conditions
       if (fullpbc .and. &
            neutralisation_method == &
            dl_mg_neutralise_none) then
          call handle_error(DL_MG_ERR_UNSPECIFIED, &
               msg='full PBC is not consistent with un-neutralised source')
       else if (.not. fullpbc .and. &
            neutralisation_method /= &
            dl_mg_neutralise_none) then
          call handle_error(DL_MG_ERR_UNSPECIFIED, &
               msg='neutralisation other than none must be used only with full PBCs')
          if ( ierror /= DL_MG_SUCCESS ) return
       end if
    end if

    select case (neutralisation_method)
    case (dl_mg_neutralise_with_ions_fixed)
       if (.not. present(neutralisation_ion_ratios)) then
          call handle_error(DL_MG_ERR_UNSPECIFIED, &
               msg='neutralisation with fixed ions needs the neutralisation_ion_ratios array')
       else
          if (any(neutralisation_ion_ratios < 0.0_wp)) then
             call handle_error(DL_MG_ERR_UNSPECIFIED, &
                  msg='negative value in the neutralisation_ion_ratios')
          end if
       end if
    end select

    call init_nonlinear(lineq, lambda, temp, n, q, c, &
         has_steric_weight, expcap, &
         neutralisation_method, neutralisation_ion_ratios)

    ! set initial concentration shift when needed
    call init_counterion_neutralisation(neutralisation_ion_ratios, ierror)
    ! needed if dl_mg return control to the callee  app
    ! instead of aborting
    if ( ierror /= DL_MG_SUCCESS ) return

    ! this needs to be done after setting neutralisation because
    ! the working ion concentration might be changed by neutralisation
    ! initialisation
    call set_debye_length(eq_type)

    call write_nonlinear_params(eq_type, n, c, q, &
            temp, lambda)

    if (fullpbc) then
       call charge_neutralisation_info
    end if

  end subroutine dl_mg_init_nonlin


  subroutine compute_rho_sum(rho)
    use dl_mg_params, only : wp
    use dl_mg_common_data, only : isx, isy, isz, fullpbc, mg_levels
    use dl_mg_grids, only : mg
    use dl_mg_types, only : mg_get_idx
     use dl_mg_nonlinear_model, only : rho_sum
    use dl_mg_kernels, only : blocked_dotproduct
    implicit none

    real(wp), intent(in) :: rho(isx:,isy:,isz:)

    integer sx, ex, sy, ey, sz, ez

    if (fullpbc) then
       call mg_get_idx(mg(mg_levels), sx, ex, sy, ey, sz, ez)
       rho_sum = blocked_dotproduct(rho(sx:, sy:, sz:), exyz=[ex,ey,ez])
    end if

  end subroutine compute_rho_sum

  !> computes the accesible volume ~ \int_V steric_weight(r) dv / (V/(nx*ny*nz))
  !! Attention: this sum is not multipliplied with dv = dx * dy * dz
  !! because it's used as a normalisation
  !! factor for quatites which are also not multiplied with dv
  subroutine compute_steric_weight_sum(has_stw, steric_weight)
    use dl_mg_params, only : wp
    use dl_mg_common_data, only : isx, isy, isz, mg_levels, &
         nx, ny, nz
    use dl_mg_grids, only : mg
    use dl_mg_types, only : mg_get_idx
    use dl_mg_kernels, only : blocked_dotproduct
    use dl_mg_nonlinear_model, only : steric_weight_sum
    implicit none

    logical, intent(out) :: has_stw
    real(wp), optional, intent(in) :: steric_weight(isx:, isy:, isz:)

    integer sx, sy, sz, ex, ey, ez

    has_stw = .false.
    if (present(steric_weight)) then
       ! sometimes it is more convenient for the calling app
       ! to pass a zero size steric_weight instead of no argument
       if (size(steric_weight) > 0) then
          has_stw = .true.
       end if
    end if

    if (has_stw)  then
       call mg_get_idx(mg(mg_levels), sx, ex, sy, ey, sz, ez)
       steric_weight_sum = blocked_dotproduct( &
            steric_weight(sx:,sy:,sz:), exyz=[ex,ey,ez])
    else
       steric_weight_sum = real(nx,wp) * ny * nz
    end if

  end subroutine compute_steric_weight_sum

  !> called by pbe solver
  subroutine test_charge_neutrality(temp, rho, pot, steric_weight)
    use dl_mg_mpi_header
    use dl_mg_mpi, only : get_myid
    use dl_mg_params
    use dl_mg_common_data, only : mg_comm, dx, dy, dz, bc, &
         isx,iex,isy,iey,isz,iez
    use dl_mg_nonlinear_model, only : mu, c=>cion, q=>qion, cin => cion_in, &
         capexp, lambda, steric_weight_sum, rho_sum, &
         beta, has_steric_weight, &
         n => nion, neutralisation_method
    use dl_mg_neutralisation_with_ions, only : xion => counterion_ratios
    use dl_mg_info, only : charge_neutrality_report
    use dl_mg_mpi, only : get_mpi_grid
    implicit none

    real(wp), intent(in) :: temp
    real(wp), dimension(:,:,:), intent(in) :: rho, pot, steric_weight
    optional steric_weight

    integer i, myid, ierr, proc_coords(3), proc_dims(3), &
         s1,e1,s2,e2,s3,e3
    real(wp) lsum(3), gsum_fixed(3), gsum_ion(3), r, z
    real(wp) obc_deviation(n)
    logical periodic(3)


    !if (any(bc == DL_MG_BC_DIRICHLET)) return
    if (all(cin == 0.0_wp)) return

    call get_mpi_grid(mg_comm, coords=proc_coords, periods=periodic, dims=proc_dims)

    ! use only internal grid  points for sums
    s1 = 1; e1 = iex-isx+1
    if (.not. periodic(1)) then
       if (proc_coords(1) == 0)               s1 = s1+1
       if (proc_coords(1) == proc_dims(1) -1) e1 = e1-1
    endif

    s2 = 1; e2 = iey-isy+1
    if (.not. periodic(2)) then
       if (proc_coords(2) == 0)               s2 = s2+1
       if (proc_coords(2) == proc_dims(2) -1) e2 = e2-1
    endif

    s3 = 1; e3 = iez-isz+1
    if (.not. periodic(3)) then
       if (proc_coords(3) == 0)               s3 = s3+1
       if (proc_coords(3) == proc_dims(3) -1) e3 = e3-1
    endif

    myid = get_myid()

    ! fixed charge integrals

    lsum(1:2) = [ sum(rho(s1:e1, s2:e2, s3:e3)), sum(abs(rho(s1:e1, s2:e2, s3:e3)))]

#ifdef MPI
    call mpi_allreduce(lsum, gsum_fixed, 2, MPI_DOUBLE_PRECISION, MPI_SUM, mg_comm, ierr)
#else
    gsum_fixed(1:2) = lsum(1:2)
#endif

    gsum_fixed(1:2) = gsum_fixed(1:2) * dx * dy * dz

    ! the ions now

    lsum = 0.0_wp
    do i = 1, n
       if (has_steric_weight) then
          r = sum (steric_weight(s1:e1, s2:e2, s3:e3) &
               * capexp(- beta * q(i) * pot(s1:e1,s2:e2,s3:e3) + mu(i)))
       else
          r =  sum (capexp(- beta * q(i) * pot(s1:e1,s2:e2,s3:e3) + mu(i)))
       end if
       lsum(1) = lsum(1) + c(i) * q(i) * r
       lsum(2) = lsum(2) + c(i) * r
    end do

#ifdef MPI
    call mpi_allreduce(lsum, gsum_ion, 2, MPI_DOUBLE_PRECISION, MPI_SUM, mg_comm, ierr)
#else
    gsum_ion(1:2) = lsum(1:2)
#endif

    ! ratio between the total number of ions given by the solution and the number of ions put in
    gsum_ion(3) = gsum_ion(2) / (steric_weight_sum*sum(c(:))) ! concentrations ratio
    ! total charge and total ion number
    gsum_ion(1:2) = gsum_ion(1:2) * dx * dy * dz

    ! deviation from OBC condition when using counter ion neutralisation
    select case(neutralisation_method)
    case (dl_mg_neutralise_with_ions_fixed, &
         dl_mg_neutralise_with_ions_auto, &
         dl_mg_neutralise_with_ions_auto_linear)
       ! compute the ratios a_i/a_i+1 where
       ! a_i = (c_shifted_i e^mu_i / c_infty_i)**(1/q_i)
       ! and compare them to 1
       ! the last element in the aaray is the sum of xions
       ! which should be close to 1
       do i = 1,n-1
          obc_deviation(i) = ((1.0_wp - xion(i) * rho_sum &
               /(steric_weight_sum * q(i)) * cin(i)) &
               * exp(mu(i)))**(1.0_wp/q(i)) &
               / ((1.0_wp - xion(i+1) * rho_sum &
               /(steric_weight_sum * q(i+1) * cin(i+1))) &
               * exp(q(i)* mu(i+1)))**(1.0_wp/q(i+1)) &
               -1.0_wp
       end do
       obc_deviation(n) = sum(xion)
    end select

    call charge_neutrality_report(gsum_fixed, gsum_ion,&
         obc_deviation, myid)

  end subroutine test_charge_neutrality


  !> removes zero modes from source term. Needed when using full PBC
  subroutine remove_zero_mode_source(mg, eq_type, remove_mode, rho, steric_weight, der_pot)
    use dl_mg_params, only : wp, EQ_POISSON, EQ_LINEAR_PBE, &
         EQ_LINEAR_PBE_NOSTERIC, &
         EQ_PBE_NEWTON, EQ_PBE_NEWTON_NOSTERIC, EQ_LINEAR_PBE_POT
    use dl_mg_mpi_header
    use dl_mg_types, only : mg_t, mg_get_idx
    use dl_mg_common_data, only : bc, nx, ny, nz, blk, &
         fullpbc, isx, iex, isy, iey, isz, iez
    use dl_mg_nonlinear_model, only : lambda, has_steric_weight, &
         steric_weight_sum,  rho_sum, neutralisation_method
    use dl_mg_kernels, only : blocked_update2, blocked_dotproduct, z_vtype
    use dl_mg_errors
    implicit none

    type(mg_t), target, intent(inout) :: mg
    integer, intent(in) :: eq_type
    integer, intent(in) :: remove_mode !< 1 (first call) or
                                       !! 2 for call in inner loops (defco or newton)
    real(wp), target, intent(inout), optional :: rho(isx:, isy:, isz:)
    real(wp), target, intent(in), optional :: steric_weight(isx:, isy:, isz:)
    real(wp), intent(in), optional         :: der_pot(isx:, isy:, isz:)

    integer sx, ex, sy, ey, sz, ez, ierr, n
    real(wp) lsum(2), gsum(2)
    real(wp), allocatable, save, target :: derp_stw(:,:,:)
    real(wp), pointer :: rho_pt(:,:,:), stw_pt(:,:,:), derp_stw_pt(:,:,:)
    real(wp) qlam

    ! this is needed only with full PBC
    if (.not. fullpbc) return

    if (present(rho)) then
       rho_pt => rho
       qlam = 1.0_wp
    else
       rho_pt => mg%f
       qlam = lambda
    endif

    if (has_steric_weight) then
       if (present(steric_weight)) then
          stw_pt => steric_weight
       else
          stw_pt => mg%d
       endif
    end if

    call mg_get_idx(mg, sx, ex, sy, ey, sz, ez)

    select case (remove_mode)
    case(1)
       select case(eq_type)
       case (EQ_POISSON)
          call add_jellium

       case(EQ_LINEAR_PBE, EQ_LINEAR_PBE_NOSTERIC)
          select case(neutralisation_method)
          case (dl_mg_neutralise_with_ions_fixed, &
               dl_mg_neutralise_with_ions_auto, &
               dl_mg_neutralise_with_ions_auto_linear)
             ! the linearised pbe with PBC has the zero order term from the
             ! exponential expansion
             ! which is like a jellium but only in the accesible volume
             !!!! probalbe the expression form below are wrong !!!!
             if ( has_steric_weight) then
                !rho_pt(sx:ex, sy:ey, sz:ez) = rho_pt(sx:ex, sy:ey, sz:ez) &
                !     - qlam * rho_sum / steric_weight_sum &
                !     * stw_pt(sx:ex, sy:ey, sz:ez)
                call blocked_update2(mg, &
                     rho_pt(sx:, sy:, sz:), z_vtype, &
                     - qlam * rho_sum / steric_weight_sum, &
                     stw_pt(sx:, sy:, sz:), z_vtype)
             else
                !rho_pt(sx:ex, sy:ey, sz:ez) = rho_pt(sx:ex, sy:ey, sz:ez) &
                !     - qlam * rho_sum
                call blocked_update2(mg, &
                     rho_pt(sx:, sy:, sz:), z_vtype, &
                     - qlam * rho_sum)
             end if

          case default
             call add_jellium
          end select

          !write(*,*) 'remove mode linPBE', remove_mode, sum(rho_pt(sx:ex, sy:ey, sz:ez))/steric_weight_sum

       case(EQ_PBE_NEWTON, EQ_PBE_NEWTON_NOSTERIC)
          select case (neutralisation_method)
          case(dl_mg_neutralise_with_jellium_uniform)
             call add_jellium(smart=.false.)

          case (dl_mg_neutralise_with_jellium_vacc)
             call add_jellium(smart=.true.)
          end select

       case (EQ_LINEAR_PBE_POT)
          ! Not sure about this
          ! we don't solve linear PBE at a non-zero pot
          !>\todo {analyse furter eq_linear_pbe_pot in as primary input}
          !
          ! this is needed in defect correction and newton
          ! both cases have to solve the following linearised equation
          ! for v
          ! A v -N'(u) v = f + N(u) - Au
          ! with full PBC the following condition needs to be satisfied
          ! -\int_V N'(u) v = \int_V ( f+N(u) -Au)
          ! if we pick v such that lhs = 0
          ! then one has to add to rhs
          ! - N'(u) * (\int_V (f + N(u) -Au)) / \int_V N'(u)
          !$OMP MASTER
          if (present(der_pot)) then
             ! need to evaluate the coefficent N'(u)
             ! this is for Newton loop
             allocate(derp_stw(isx:iex,isy:iey,isz:iez))
             call derpot_and_steric_weight
             derp_stw_pt => derp_stw
          else
             ! use the one stored
             ! this is for defco
             derp_stw_pt => mg%w
          end if
          lsum(1) = sum(rho_pt(sx:ex, sy:ey, sz:ez))
          lsum(2) = sum(derp_stw_pt(sx:ex, sy:ey, sz:ez))
#ifdef MPI
          call mpi_allreduce(lsum, gsum, 2, MPI_DOUBLE_PRECISION, MPI_SUM, mg%comm, ierr)
#else
          gsum = lsum
#endif
          rho_pt(sx:ex, sy:ey, sz:ez) = rho_pt(sx:ex, sy:ey, sz:ez) &
               - derp_stw_pt(sx:ex,sy:ey,sz:ez) * gsum(1) / gsum(2)

          !$omp end master
          !$omp barrier
       end select

    case(2)
       ! inner loops in Newton and defect correction
       select case(eq_type)
       case (EQ_POISSON)
          call add_jellium

       case(EQ_LINEAR_PBE, EQ_LINEAR_PBE_NOSTERIC)
          call add_jellium

       case(EQ_PBE_NEWTON, EQ_PBE_NEWTON_NOSTERIC)
          ! mode 2 shouldn't hit this, set a test
          select case (neutralisation_method)
          case (dl_mg_neutralise_with_jellium_uniform, &
               dl_mg_neutralise_with_jellium_vacc)
             call add_jellium
          end select

       case (EQ_LINEAR_PBE_POT)
          ! this is needed in defect correction and Newton
          ! both cases have to solve the following linearised equation
          ! for v
          ! A v -N'(u) v = f + N(u) - Au
          ! with full PBC the following condition needs to be satisfied
          ! -\int_V N'(u) v = \int_V ( f+N(u) -Au)
          ! if we pick v such that lhs = 0
          ! then one has to add to rhs
          ! - N'(u) * (\int_V (f + N(u) -Au)) / \int_V N'(u)

          call add_jellium

          ! the variant from bellow caused
          ! diffs with 2.1.4 in ONETEP regression test
          ! I don't think that we need to be so smart
          ! to be revised
       end select

    end select

  contains

    subroutine add_jellium(smart)
!$    use omp_lib      
      implicit none

      logical, optional, intent(in) :: smart

      real(wp)  gsum
      logical smar_

      if (present(smart)) then
         smar_ = smart
      else
         smar_ =.false.
      endif

      gsum = blocked_dotproduct(rho_pt(sx:, sy:, sz:), scale=1.0_wp)
      if (smar_ .and. has_steric_weight) then
         call blocked_update2(mg, &
              rho_pt(sx:, sy:, sz:), z_vtype, &
              - gsum / steric_weight_sum, stw_pt(sx:, sy:, sz:))
      else
         call blocked_update2(mg, &
              rho_pt(sx:, sy:, sz:), z_vtype, &
              - gsum / (real(nx,wp)*ny*nz))
      endif

    end subroutine add_jellium


    ! used in master region
    subroutine derpot_and_steric_weight
      use dl_mg_nonlinear_model, only : derfnl
      implicit none

      integer i,j,k

      if (has_steric_weight)then
         do k= sz, ez
            do j = sy, ey
               do i = sx, ex
                  derp_stw(i,j,k) = derfnl(der_pot(i,j,k), mg%d(i,j,k))
               enddo
            enddo
         enddo
      else
         do k= sz, ez
            do j = sy, ey
               do i = sx, ex
                  derp_stw(i,j,k) = derfnl(der_pot(i,j,k), 1.0_wp)
               enddo
            enddo
         enddo
      end if

    end subroutine derpot_and_steric_weight

  end subroutine remove_zero_mode_source


  !> removes zero modes from potential. Needed when using full PBC
  !! at the moment this is called inside vcycle for mg%p
  !! it migh need to be extend for arbitrary arrays with optional
  !! args
  subroutine remove_zero_mode_potential(mg, eq_type)
    use dl_mg_params, only :wp,  DL_MG_BC_DIRICHLET, DL_MG_BC_PERIODIC, &
         EQ_PBE_NEWTON, EQ_PBE_NEWTON_NOSTERIC, EQ_LINEAR_PBE_POT, EQ_POISSON, &
         EQ_LINEAR_PBE, EQ_LINEAR_PBE_NOSTERIC
    use dl_mg_mpi_header
    use dl_mg_types, only : mg_t, mg_get_idx
    use dl_mg_nonlinear_model, only : derfnl, has_steric_weight, &
         steric_weight_sum
    use dl_mg_common_data, only : bc, blk, t => mg_levels, nx, ny, nz

    use dl_mg_errors
    implicit none

    type(mg_t), target, intent(inout) :: mg
    integer, intent(in) :: eq_type
    !real(wp), target, optional, intent(inout) :: v(:,:,:)

    integer sx, ex, sy, ey, sz, ez, ierr
    integer n, i, j, k, d, iblk
    integer s(3), e(3), si(3), ei(3)
    real(wp) :: lsum(2), gsum(2), q, pot_shift
    save lsum, gsum ! for OpenMP reduction
    real(wp), pointer :: pp(:,:,:)
#ifdef HAVE_CONTIGUOUS
    contiguous pp
#endif
    
    if ( any(bc /= DL_MG_BC_PERIODIC)) return

    call mg_get_idx(mg, sx, ex, sy, ey, sz, ez)
    si = [sx,sy,sz]
    ei = [ex,ey,ez]

    ! not need of external array so far
    !if (present(v)) then
       ! needs halo ? not used so far
       !pp(mg%psx:, mg%psy:, mg%psz:) => v
    !   continue
    !else
       pp => mg%p
    !end if
    
    !$OMP SINGLE
    lsum = 0.0_wp
    !$OMP END SINGLE

    if (eq_type == EQ_LINEAR_PBE_POT) then
       ! used in defect correction
       n = 2
       !$OMP DO SCHEDULE(STATIC,1) REDUCTION(+:lsum)
       do iblk = 1, size(blk(t)%start, dim=2)
          s(:) = blk(t)%start(:,iblk)
          e(:) = blk(t)%end(:,iblk)
          do k = s(3), e(3)
             do j = s(2), e(2)
                do i = s(1), e(1)
                   q = mg%w(i,j,k)
                   lsum(1) = lsum(1) + pp(i, j, k) * q
                   lsum(2) = lsum(2) + q
                enddo
             end do
          end do
       enddo
       !$OMP ENDDO
    else
       if (has_steric_weight) then
          if (eq_type == EQ_PBE_NEWTON ) then
             ! inside the Newton solver, NB: PBE is linearised around mg%w
             n = 2
             !$OMP DO SCHEDULE(STATIC,1) REDUCTION(+:lsum)
             do iblk = 1, size(blk(t)%start, dim=2)
                s(:) = blk(t)%start(:,iblk)
                e(:) = blk(t)%end(:,iblk)
                do k = s(3), e(3)
                   do j = s(2), e(2)
                      do i = s(1), e(1)
                         q = derfnl(mg%w(i,j,k), mg%d(i,j,k))
                         lsum(1) = lsum(1) + pp(i, j, k) * q
                         lsum(2) = lsum(2) + q
                      enddo
                   end do
                end do
             end do
             !$OMP ENDDO
          else
             ! this is for linearised PBE, steric_weight_sum is computed only
             ! once in dl_mg_init_nonlinear
             n = 1
             !$OMP DO SCHEDULE(STATIC,1) REDUCTION(+:lsum)
             do iblk = 1, size(blk(t)%start, dim=2)
                s(:) = blk(t)%start(:,iblk)
                e(:) = blk(t)%end(:,iblk)
                do k = s(3), e(3)
                   do j = s(2), e(2)
                      do i = s(1), e(1)
                         lsum(1) = lsum(1) + mg%d(i,j,k) * pp(i, j, k)
                      enddo
                   end do
                end do
             end do
             !$OMP ENDDO
          end if
       else
          if (eq_type == EQ_PBE_NEWTON_NOSTERIC) then
             n = 2
             !$OMP DO SCHEDULE(STATIC,1) REDUCTION(+:lsum)
             do iblk = 1, size(blk(t)%start, dim=2)
                s(:) = blk(t)%start(:,iblk)
                e(:) = blk(t)%end(:,iblk)
                do k = s(3), e(3)
                   do j = s(2), e(2)
                      do i = s(1), e(1)
                         q = derfnl(mg%w(i,j,k), 1.0_wp)
                         lsum(1) = lsum(1) + pp(i, j, k) * q
                         lsum(2) = lsum(2) + q
                      enddo
                   end do
                end do
             end do
             !$OMP ENDDO
          else
             ! linearised PBE without steric weight or simply Poisson Eq
             n = 1
             !$OMP DO SCHEDULE(STATIC,1) REDUCTION(+:lsum)
             do iblk = 1, size(blk(t)%start, dim=2)
                s(:) = blk(t)%start(:,iblk)
                e(:) = blk(t)%end(:,iblk)
                do k = s(3), e(3)
                   do j = s(2), e(2)
                      do i = s(1), e(1)
                         lsum(1) = lsum(1) + pp(i, j, k)
                      end do
                   end do
                end do
             enddo
             !$OMP ENDDO
          end if
       end if
    end if

    !$OMP MASTER
#ifdef MPI
    call mpi_allreduce(lsum, gsum, n, MPI_DOUBLE_PRECISION, MPI_SUM, mg%comm, ierr)
#else
    gsum = lsum
#endif
    !$OMP END MASTER
    !$OMP BARRIER

    if ( n == 2) then
       pot_shift = gsum(1) / gsum(2)
    else
       if ( eq_type == EQ_POISSON) then
          pot_shift = gsum(1) / (real(nx,wp) * ny *nz)
       else
          ! if no steric the steric_weight_sum = nx * ny * nz
          pot_shift = gsum(1) / steric_weight_sum
       end if
    end if

    !mg%p(sx:ex, sy:ey, sz:ez) = mg%p(sx:ex, sy:ey, sz:ez) - pot_shift

    !$OMP DO SCHEDULE(STATIC,1)
    do iblk = 1, size(blk(t)%start, dim=2)
       s(:) = blk(t)%start(:,iblk)
       e(:) = blk(t)%end(:, iblk)
       ! we need to cover the halos as well at least in multigrid
       do d=1,3
          if ( s(d) == si(d) ) s(d) = s(d) -1
          if ( e(d) == ei(d) ) e(d) = e(d) +1
       end do
       do k = s(3), e(3)
          do j = s(2), e(2)
             do i = s(1), e(1)
                pp(i, j, k) = pp(i, j, k) - pot_shift
             enddo
          enddo
       enddo
    enddo
    !$OMP ENDDO


    !if ( eq_type == EQ_LINEAR_PBE_POT .and. use_counterions ) then
!!$    if (.false.)then
!!$    ! need to shift the potential with the avarage of lhs
!!$       !$OMP DO SCHEDULE(STATIC,1) REDUCTION(+:lsum)
!!$       do iblk = 1, size(blk(t)%start, dim=2)
!!$          s(:) = blk(t)%start(:,iblk)
!!$          d(:) = blk(t)%dims(:) - 1
!!$          do k = s(3), min(s(3)+d(3), ez)
!!$             do j = s(2), min(s(2)+d(2), ey)
!!$                do i = s(1), min(s(1)+d(1), ex)
!!$                   lsum(1) = lsum(1) + mg%f(i, j, k)
!!$                enddo
!!$             end do
!!$          end do
!!$       enddo
!!$       !$OMP ENDDO
!!$
!!$       !$OMP MASTER
!!$#ifdef MPI
!!$       call mpi_allreduce(lsum, gsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mg%comm, ierr)
!!$#else
!!$       gsum = lsum
!!$#endif
!!$       !$OMP END MASTER
!!$       !$OMP BARRIER
!!$
!!$       pot_shift = gsum(1) / gsum(2)
!!$
!!$       write(*,*) 'pot shift', pot_shift
!!$
!!$       !$OMP DO SCHEDULE(STATIC,1)
!!$       do iblk = 1, size(blk(t)%start, dim=2)
!!$          s(:) = blk(t)%start(:,iblk)
!!$          d(:) = blk(t)%dims(:) - 1
!!$          do k = s(3), min(s(3)+d(3), ez)
!!$             do j = s(2), min(s(2)+d(2), ey)
!!$                do i = s(1), min(s(1)+d(1), ex)
!!$                   pp(i, j, k) = pp(i, j, k) + pot_shift
!!$                enddo
!!$             enddo
!!$          enddo
!!$       enddo
!!$       !$OMP ENDDO
!!$
!!$    end if

  end subroutine remove_zero_mode_potential


  real(wp) function compute_norm(mg, component, scale)
    use dl_mg_types, only : mg_t
    use dl_mg_mpi
    use dl_mg_common_data, only : mg_levels, dx, dy, dz, blk
    use dl_mg_types, only : mg_get_idx
    use dl_mg_kernels, only : blocked_dotproduct, z_vtype
    use dl_mg_errors, only : handle_error, DL_MG_ERR_UNSPECIFIED

    type(mg_t), intent(in), target :: mg
    character(len=1), intent(in) :: component
    logical, optional, intent(in) :: scale ! flag for scaled norm

    integer level, sx, ex, sy, ey, sz, ez, ierr
    real(wp) xscale
    real(wp), pointer :: a_ptr(:,:,:)
    real(wp) :: sumtot

    !> \todo{a_ptr can have contiguous attribute}

    select case ( component(1:1) )
    case ("f")
       a_ptr => mg%f
    case ("r")
       a_ptr => mg%r
    case("p")
       a_ptr => mg%p
    case("w")
       a_ptr => mg%w
    case default
       call handle_error(DL_MG_ERR_UNSPECIFIED, msg="wrong filed component in compute_norm")
    end select

    call mg_get_idx (mg, sx, ex, sy, ey, sz, ez)
    level = mg%level

    ! I think that it is correct to use identical arguments below because they
    ! have intent in
    sumtot = blocked_dotproduct(mg, a_ptr(sx:,sy:,sz:), z_vtype, &
         a_ptr(sx:,sy:,sz:))

!!$    !$OMP DO SCHEDULE(STATIC,1) REDUCTION(+:sumloc)
!!$    do iblk = 1, size(blk(level)%start, dim=2)
!!$       s(:) = blk(level)%start(:,iblk)
!!$       d(:) = blk(level)%dims(:) - 1
!!$       !if ( level == 1) then
!!$       !  write(0,*)'norm ', s, d, sumloc, component, mg%f(2:2,2:2,2:2), a_ptr
!!$       ! endif
!!$       do k = s(3), min(s(3) + d(3), ez)
!!$          do j = s(2), min(s(2) + d(2), ey)
!!$             do i = s(1), min(s(1) + d(1), ex)
!!$                !write(0,*) 'norm ', i,j,k, sumloc, a_ptr(i,j,k)
!!$                sumloc = sumloc + a_ptr(i, j, k)**2
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo
!!$    !$OMP ENDDO
!!$
!!$
!!$    !$OMP MASTER
!!$#ifdef MPI
!!$    call mpi_allreduce(sumloc, sumtot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mg%comm, ierr)
!!$#else
!!$    sumtot = sumloc
!!$#endif
!!$    !$OMP END MASTER
!!$    !$OMP BARRIER

    xscale = 1.0_wp
    if ( present(scale) ) then
       if (scale) then
          xscale = ( 2**(mg_levels - level) )**3
       end if
    endif

    compute_norm = sqrt(xscale * dx * dy * dz * sumtot)

  end function compute_norm


  !> needed in dl_mg_solver_pbe
  !! if the arguments betamu or gammas are present set them to the values
  !! used in code
  !! added the ion concentrartions used internally (it might differ from the input
  !! value for counter-ion neutralisation, when using PBC)
  subroutine return_betamu_gamma_ion_concentrations(n, betamu, ion_concentration, &
       neutralisation_ratios, steric_weight_average)
    use dl_mg_params, only : wp
    use dl_mg_common_data, only     : nx, ny, nz
    use dl_mg_nonlinear_model, only : mu, cion, steric_weight_sum
    use dl_mg_neutralisation_with_ions, only : counterion_ratios
    implicit none
    integer, intent(in)             :: n
    real(wp), optional, intent(out) :: betamu(n), ion_concentration(n)
    real(wp), optional, intent(out) :: steric_weight_average
    real(wp), optional, intent(out) :: neutralisation_ratios(:)

    if (present(betamu)) then
       betamu = mu
    end if

    if (present(ion_concentration)) then
       ion_concentration = cion
    end if

    if (present(steric_weight_average)) then
       steric_weight_average = steric_weight_sum / (real(nx,wp) * ny * nz)
    end if

    if (present(neutralisation_ratios) &
         .and. allocated(counterion_ratios)) then
       neutralisation_ratios = counterion_ratios
    end if

  end subroutine return_betamu_gamma_ion_concentrations

  !> checks if the first error code bubbled up to the top level or if it
  !! was replaced
  subroutine check_error_code(ierror)
    use dl_mg_params, only : DL_MG_SUCCESS
    use dl_mg_common_data, only : report_unit
    use dl_mg_mpi, only : get_myid
    use dl_mg_errors, only : first_error, abort_on_failure
    implicit none
    integer, optional , intent(inout) :: ierror

    if ( present(ierror)) then
       if ( .not. abort_on_failure           &
            .and. ierror /= DL_MG_SUCCESS    &
            .and. ierror /= first_error      &
            .and. get_myid() == 0) then
          write (report_unit, *) "WARNING: The were other errors found on the way up."
          write (report_unit,*) "first error code : ", first_error
          write (report_unit,*) "returned ierror:   ", ierror
          ierror = first_error
       end if
     end if
   end subroutine check_error_code


   !> convergence tests, for details see dl_mg_convergence_params
   subroutine mg_convergence_test(mg, level, eq_type, test_tag, &
       iter, converged, f_norm, res_norm, sol)
    ! test convergences with criterion
    ! || r_k || < max (tol * || f ||, tolmin)
    !  the CFD criterion
    ! || r_k || < tol * || r_0||
    ! is deemed as dubious
    ! see paper by Gereald Recktenwald Google docs
    !
    ! WARNING : all iterations that call this subroutine must start from iter = 1
    !$    use omp_lib
    use dl_mg_info
    use dl_mg_params, only : wp, DL_MG_SUCCESS, DL_MG_ERR_NITER, &
         DL_MG_ERR_RES_RATIO
    use dl_mg_types, only : mg_t
    use dl_mg_common_data, only : mg_levels, blk, fullpbc, isx, isy,isz
    use dl_mg_kernels, only : blocked_dotProduct
    use dl_mg_alloc, only : p_vtype, f_vtype, r_vtype, z_vtype
    use dl_mg_timer
    use dl_mg_convergence_params, only :&
         tconv_level_1, tconv_vcycle, &
         tconv_newton, tconv_cg, nconv_tags, conv_params, &
         mg_max_conv_rate => mg_max_conv_rate_ , &
         tol_mu_rel_, tol_mu_abs_
    use dl_mg_nonlinear_model, only :  mu, mu_prev, neutralisation_method

    implicit none


    type (mg_t), target, intent(inout) :: mg
    integer, intent(in)        :: level, eq_type, test_tag, iter
    integer, intent(out)       :: converged
    real(wp), optional, intent(in) :: f_norm, res_norm
    real(wp), optional, target, intent(in) :: sol(isx:,isy:,isz:) ! current solution

    real(wp), save :: fnorm(nconv_tags)
    real(wp), save :: rnorm_m1(nconv_tags) ! previous residual norm
    !$omp THREADPRIVATE(fnorm, rnorm_m1)
    real(wp) rnorm, pnorm, tol_used, res_ratio, dmu
    integer info_eq_tag, maxiter
    logical in_newton_loop, in_vcycle_loop, in_cg_loop, collect_data, mu_converged
    real(wp), pointer :: sol_pt(:,:,:)

    call mg_timer(start, tconv, level, tcompute)

    converged = -1 ! working, no success (0) or error (>0)

    in_newton_loop = test_tag == tconv_newton
    in_vcycle_loop = test_tag == tconv_vcycle
    in_cg_loop     = test_tag == tconv_cg
    collect_data   = in_newton_loop .or. in_vcycle_loop .or. in_cg_loop

    if ( in_newton_loop) then
       info_eq_tag = info_newton_tag
    else if( in_vcycle_loop) then
       info_eq_tag = info_mg_tag
    else if ( test_tag == tconv_level_1) then
    else if (test_tag == tconv_cg) then
       info_eq_tag = info_cg_tag
    endif

    ! for reporting and debug
    if (present(sol)) then
       sol_pt =>  sol(mg%sx:,mg%sy:,mg%sz:)
    else
       sol_pt => mg%p(mg%sx:,mg%sy:,mg%sz:)
    endif

    
    if ( iter == 1 ) then
       if (present(f_norm)) then
          fnorm(test_tag) = f_norm
       else
          fnorm(test_tag) = sqrt(blocked_dotproduct(mg, mg%f, f_vtype, mg%f, f_vtype))
       endif
    endif
    tol_used =  max(conv_params(test_tag)%abs, conv_params(test_tag)%rel * fnorm(test_tag) )
    maxiter = conv_params(test_tag)%maxiter

    if (present(res_norm)) then
       rnorm = res_norm
    else
       ! write(0,*) 'convergence test 1', iter, omp_get_thread_num(), rnorm_m1
       !rnorm = compute_norm(mg, "r")
       rnorm = sqrt( blocked_dotproduct(mg, mg%r, r_vtype, mg%r, r_vtype))
    endif
    !write(0,*) 'convergence test', iter, fnorm_coarse, tol_used, rnorm, accuracy_limit
    !write(0,*) 'convergence test 2', iter, omp_get_thread_num(), rnorm_m1

    if ( iter == 1 ) rnorm_m1(test_tag) = 10.0_wp*rnorm

    !write(0,*) 'convergence test 3', iter, omp_get_thread_num(), rnorm_m1

    if ( rnorm_m1(test_tag) > 0.0_wp) then
       res_ratio = rnorm/rnorm_m1(test_tag)
    else
       res_ratio = -1.0_wp
    endif

    ! sort out mu convergence if is needed
    if (fullpbc .and. in_newton_loop) then
       if ( all( abs(mu_prev-mu) < max(tol_mu_rel_ * mu, tol_mu_abs_))) then
          mu_converged = .true.
       else
          mu_converged = .false.
       end if
       !$OMP MASTER
       dmu = maxval(abs(mu-mu_prev)) ! for reporting
       !$OMP END MASTER
    else
       mu_converged = .true. ! i.e. we don't care
       dmu = -1.0_wp
    end if

    !write(0,*) 'convergence test 4', iter, test_tag, rnorm, rnorm_m1(test_tag), tol_used, fnorm(test_tag), tol
    if ( debug_write ) then
       pnorm = sqrt(blocked_dotproduct(mg, sol_pt, z_vtype, sol_pt, z_vtype))
       !$OMP MASTER
       call mg_info_iteration_details(test_tag, level, iter, rnorm, res_ratio, pnorm, tol_used)
       !$OMP END MASTER
    endif

    if (collect_data) then
       pnorm = sqrt(blocked_dotproduct(mg, sol_pt, z_vtype, sol_pt, z_vtype))
       !$OMP MASTER
       call accumulate_residual_ratio(info_eq_tag, iter, rnorm, pnorm, res_ratio, dmu)
       !$OMP END MASTER
    endif


    !if ( in_newton_loop) then
    !   write(*,*) 'convergence test', iter, omp_get_thread_num(), rnorm, tol_used, mu_converged, mu, mu_prev!, rnorm/rnorm_m1
    !end if

    if ( rnorm < tol_used .and. mu_converged) then
       converged = DL_MG_SUCCESS

       if( collect_data) then
          !pnorm = compute_norm(mg, field_component)
          !$OMP MASTER
          call mg_info_set(DL_MG_SUCCESS, iter, rnorm, res_ratio,&
               pnorm, tol_used, maxiter, info_eq_tag)
          !$OMP END MASTER
       end if
    else

       ! test for some error conditions
    
       ! for multigrid algorithm:
       ! if residual ratio is larger than mg_max_conv_rate is assumed that
       ! the desired tolerance target cannot be achieved due to
       ! round-off errors or due to algorithm unsuitability
       ! for the problem
       ! the comparison id done after 3 steps to allow the residual ratio to settle
       if ( test_tag == tconv_vcycle .and. &
            iter > 3 .and. res_ratio > mg_max_conv_rate) then
          converged = DL_MG_ERR_RES_RATIO
          if ( collect_data) then
             !pnorm = compute_norm(mg, field_component)
             !$OMP MASTER
             call mg_info_set(DL_MG_ERR_RES_RATIO, iter, rnorm, res_ratio,&
                  pnorm, tol_used, maxiter, info_eq_tag)
             !$OMP END MASTER
          endif
       endif

       if (iter == maxiter ) then
          converged = DL_MG_ERR_NITER
          if ( collect_data ) then
             !pnorm = compute_norm(mg, field_component)
             !$OMP MASTER
             call mg_info_set(DL_MG_ERR_NITER, iter, rnorm, res_ratio, &
                  pnorm, tol_used, maxiter, info_eq_tag)
             !$OMP END MASTER
          endif
       endif
    end if

    rnorm_m1(test_tag) = rnorm

    if (fullpbc .and. test_tag == tconv_newton) then
       !$OMP SINGLE
       mu_prev = mu
       !$OMP END SINGLE
    end if


    call test_thread_consistency

    call mg_timer(stop, tconv, level, tcompute)

  contains

    !> check if the converged flag is the same across threads
    subroutine test_thread_consistency
      use dl_mg_common_data, only : nthreads
      use dl_mg_errors, only : handle_error
      use dl_mg_params, only : DL_MG_ERR_UNSPECIFIED
      implicit none
      integer, save :: cmax, cmin
      integer tid

      if ( nthreads == 1) return
      tid = 0
      !$ tid = omp_get_thread_num()
      !$OMP MASTER
      cmax = converged
      cmin = converged
      !$OMP END MASTER
      !$OMP BARRIER
      !$OMP CRITICAL
      if (tid > 0) then
         cmax = MAX(cmax, converged)
         cmin = MIN(cmin, converged)
      endif
      !$OMP END CRITICAL
      !$OMP MASTER
      if (cmax /= cmin ) then
         call handle_error(DL_MG_ERR_UNSPECIFIED, msg="mg_convergence_test:&
              & convergence flags differ on different&
              & threads at top level")
      endif
      !$OMP END MASTER
      !$OMP BARRIER

    end subroutine test_thread_consistency

  end subroutine mg_convergence_test    

  
  subroutine set_boundary_values(mg, v, v_type, val, aval, av_type)
    use dl_mg_mpi, only : get_mpi_grid
    use dl_mg_params, only : wp, DL_MG_BC_DIRICHLET, DL_MG_ERR_UNSPECIFIED
    use dl_mg_alloc, only : get_lbound
    use dl_mg_types, only : mg_t, mg_get_idx
    use dl_mg_common_data, only : bc
    use dl_mg_errors, only : handle_error
    implicit none

    type(mg_t), intent(inout) :: mg
    real(wp), intent(inout) :: v(:,:,:)
    integer, intent(in) :: v_type
    real(wp), optional, intent(in) :: val
    real(wp), optional, intent(in) :: aval(:,:,:)
    integer,  optional, intent(in) :: av_type

    ! set field to 0 on boundary layers.
    ! it may be generalised later

    integer dims(3), coords(3), sx, ex, sy, ey, sz, ez
    integer  ba(3), bv(3)
    
    ! needed only for Dirichlet BC

    if (present(val) .and. present(aval)) then
       call handle_error(DL_MG_ERR_UNSPECIFIED, &
            msg=" set_boundary_values: scalar and array values provided")
    else if (.not. (present(val) .or. present(aval))) then
       call handle_error(DL_MG_ERR_UNSPECIFIED, &
            msg=" set_boundary_values: needs one of val, aval  arguments")
    end if
    if (present(aval) .and. .not. present(av_type)) then
       call handle_error(DL_MG_ERR_UNSPECIFIED, &
            msg=" set_boundary_values: aval argument provided without av_vtype")
    end if
    
    !$omp master
    call get_mpi_grid(mg%comm, dims=dims, coords=coords)
    call mg_get_idx(mg, sx, ex, sy, ey, sz, ez)

    call get_lbound(mg, v_type, bv)
    if (present(aval)) then
       call get_lbound(mg,av_type, ba)
    end if

    call set_bd(v,aval)
    !$omp end master
    !$omp barrier

  contains

    subroutine set_bd(v, av)
      implicit none
      real(wp) v(bv(1):,bv(2):,bv(3):)
      real(wp), optional :: av(ba(1):,ba(2):,ba(3):)

      
      if (present(av)) then
         if (bc(1) == DL_MG_BC_DIRICHLET) then
            if (coords(1) == 0) then
               v(sx-1, sy:ey, sz:ez ) = av(sx-1, sy:ey, sz:ez )
            endif
            if (coords(1) == dims(1)-1) then
               v(ex+1, sy:ey, sz:ez ) = av(ex+1, sy:ey, sz:ez )
            endif
         endif
         if (bc(2) == DL_MG_BC_DIRICHLET) then
            if (coords(2) == 0) then
               v(sx:ex, sy-1, sz:ez ) = av(sx:ex, sy-1, sz:ez )
            endif
            if (coords(2) == dims(2)-1) then
               v(sx:ex, ey+1, sz:ez ) = av(sx:ex, ey+1, sz:ez )
            endif
         endif
         if (bc(3) == DL_MG_BC_DIRICHLET) then
            if (coords(3) == 0) then
               v(sx:ex, sy:ey, sz-1 ) = av(sx:ex, sy:ey, sz-1 )
            endif
            if (coords(3) == dims(3)-1) then
               v(sx:ex, sy:ey, ez+1 ) = av(sx:ex, sy:ey, ez+1 )
            endif
         endif
      else
         if (coords(1) == 0 .and. bc(1) == DL_MG_BC_DIRICHLET) then
            v(sx-1, sy:ey, sz:ez ) = val
         endif
         
         if (coords(1) == dims(1) -1 .and. bc(1) == DL_MG_BC_DIRICHLET) then
            v(ex+1, sy:ey, sz:ez ) = val
         endif
         
         if (coords(2) == 0 .and. bc(2) == DL_MG_BC_DIRICHLET) then
            v(sx:ex, sy-1, sz:ez ) = val
         endif
         
         if (coords(2) == dims(2) -1 .and. bc(2) == DL_MG_BC_DIRICHLET) then
            v(sx:ex, ey+1, sz:ez ) = val
         endif
         
         if (coords(3) == 0 .and. bc(3) == DL_MG_BC_DIRICHLET) then
            v(sx:ex, sy:ey, sz-1 ) = val
         endif
         
         if (coords(3) == dims(3) -1 .and. bc(3) == DL_MG_BC_DIRICHLET) then
            v(sx:ex, sy:ey, ez+1 ) = val
         endif
      end if
      
    end subroutine set_bd

  end subroutine set_boundary_values


  !> computes mass term at non-zero potential
  !! used in linearized PBEs (Newton and deferred correction)
  subroutine set_mass_term(mg, der_pot, ierr)
    use dl_mg_common_data, only :  mg_levels, blk, isx, isy, isz
    use dl_mg_nonlinear_model, only : has_steric_weight, derfnl
    use dl_mg_types, only : mg_t
    use dl_mg_errors
    implicit none

    type(mg_t), intent(inout)       :: mg(:)
    real(wp), intent(in)  :: der_pot(isx:, isy:, isz:)
    integer, optional, intent(out)  :: ierr

    integer i, j, k, iblk, &
         s1, e1, s2, e2, s3, e3, t
    

    if (present(ierr)) then
       ierr = DL_MG_SUCCESS
    endif

    t = mg_levels
    
    !call mg_get_idx(mg(t), sx, ex, sy, ey, sz, ez)
    !$omp do schedule(static,1)
    do iblk=1, blk(t)%nblks(2)
       call set_blk_idx(iblk)
       if (has_steric_weight)then
          do k= s3, e3
             do j = s2, e2
                do i = s1, e1
                   mg(t)%w(i,j,k) = derfnl(der_pot(i,j,k), &
                        mg(t)%d(i,j,k))
                enddo
             enddo
          enddo
       else
          do k= s3, e3
             do j = s2, e2
                do i = s1, e1
                   mg(t)%w(i,j,k) = derfnl(der_pot(i,j,k), &
                        1.0_wp)
                enddo
             enddo
          enddo
       end if
    end do
    !$omp enddo
       
    do t = mg_levels-1, 1, -1
       if ( .not. mg(t)%active) exit
       !$omp do schedule(static,1)
       do iblk=1, blk(t)%nblks(2)
          call set_blk_idx(iblk)
          do k = s3, e3
             do j = s2, e2
                do i = s1, e1
                   mg(t)%w(i,j,k) = mg(t+1)%w(2*i-1, 2*j-1, 2*k-1)
                enddo
             enddo
          enddo
       enddo
       !$omp enddo
    end do

    contains

      subroutine set_blk_idx(iblk)
        implicit none

        integer, intent(in) :: iblk
        
        s1 = blk(t)%start(1, iblk)
        s2 = blk(t)%start(2, iblk)
        s3 = blk(t)%start(3, iblk)
        
        e1 = blk(t)%end(1,iblk)
        e2 = blk(t)%end(2,iblk)
        e3 = blk(t)%end(3,iblk)
        
      end subroutine set_blk_idx
      
  end subroutine set_mass_term
  
end module dl_mg_utils
