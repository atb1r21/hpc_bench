!> \brief  Higher order derivatives defect correction driver and its subroutines
!!
!> James C. Womack, Lucian Anton 2016-17
module dl_mg_defco
  use dl_mg_params, only: wp
  implicit none

  private

  public :: dl_mg_defco_defect_corr_solver

contains

  !> This subroutine computes the solution for the Poisson or Poisson-Boltzmann Equation
  !! with higher order approximation for the finite difference derivatives.
  !! It uses multigrid to obtain the solution for the second order approximation
  !! and defect correction method to improve it to the specified order.
  !> This subroutine uses the defect correction method with multigrid
  !! to obtain a solution to the Poisson or Poisson-Boltzmann equation which
  !! uses higher order approximation for the finite difference derivatives.
  !! Typically this increases solution accuracy but don't forget that
  !! higher order means does not imply higher precision in all cases.
  subroutine dl_mg_defco_defect_corr_solver(eq_type, eps_full, eps_half, alpha,  & ! out
       rho, fd_order, pot, der_pot, steric_weight, &
       use_pot_in, use_damping, res, ierror)

    use dl_mg_common_data, only: mg_comm, mg_levels, fullpbc, bc, &
         dx, dy, dz, &
         grid_weight, &
         isx, iex, isy, iey, isz, iez, &
         mx, my, mz, &
         nxb, nyb, nzb, &
         report_unit, &
         DL_MG_SUCCESS, DL_MG_BC_DIRICHLET
    use dl_mg_params, only : EQ_PBE_NEWTON, EQ_PBE_NEWTON_NOSTERIC, &
         EQ_PBE_FAS, EQ_PBE_FAS_NOSTERIC, &
         EQ_LINEAR_PBE_POT
    use dl_mg_grids, only: mg, set_fd, fd
    use dl_mg_defco_utils, only: dl_mg_defco_utils_alloc_check, &
         dl_mg_defco_utils_assert, dl_mg_defco_utils_dealloc_check, &
         dl_mg_defco_utils_integrate_product_on_grid
    use dl_mg_convergence_params, only : &
         tol_res_rel_, &
         tol_res_abs_, &
         tol_pot_rel_, &
         tol_pot_abs_, &
         tol_mu_rel_,  &
         tol_mu_abs_,  &
         defco_maxiters => max_iters_defco_
    use dl_mg_second_order_solvers, only : dlmg_2nd_order_solver
    use dl_mg_nonlinear_model, only : mu, mu_prev, &
         neutralisation_method
    use dl_mg_alloc, only : f_vtype, z_vtype
    use dl_mg_utils, only : remove_zero_mode_source, set_mass_term
    use dl_mg_kernels
    use dl_mg_mpi_header
    use dl_mg_mpi, only : get_mpi_grid, get_myid
    use dl_mg_errors
    use dl_mg_timer
    use dl_mg_info, only : info_get_res_norm, defco_status
    ! debug
    !use dl_mg_grids, only : mg

    implicit none

    integer,  intent(in)           :: eq_type !< selector for the algorithm to be used
    real(wp), intent(in)           :: eps_full(:,:,:) !< relative permitivity on the grid points
                                                      !! ignored if the finite difference order is 2
    real(wp), intent(in)           :: eps_half(:,:,:,:) !< relative permitivity, the values correspond
                                                        !! to the points located halfway between the grid points
                                                        !! in each direction

    real(wp), intent(in)           :: alpha !< multiplicative constant for \f$\rho\f$
                                            !! defined by the units used
    real(wp), target, intent(in)   :: rho(:,:,:) !< r.h.s. term (charge density)
    integer,  intent(in)           :: fd_order !< finite difference order used in defect correction
    real(wp), intent(inout)        :: pot(:,:,:) !< the potential we seek
    real(wp), intent(in), optional :: der_pot(:,:,:) !< potential function at which the
                                                     !! derivative of the Boltzmann term is
                                                     !! computed (if not present defaults to 0)
    real(wp), intent(in), optional :: steric_weight(:,:,:) !< steric weight in Boltzmann term, 
                                                           !< i.e. \f$ e^{-\beta V_{steric}(r)} \f$
    logical,  intent(in), optional :: use_pot_in !< start iterating from this potential,
                                                 !! if not present starts from 0 inside the domain
    logical,  intent(in), optional :: use_damping !< if present and true use error damping in defect correction loop
    real(wp), intent(out),optional :: res(:,:,:) !< residual
    integer, intent(out), optional :: ierror !< error flag

    target :: res

    ! Internal variables
    integer       :: ierr  ! jd: Error flag
    integer       :: ierr1 ! JCW: to satisfy MPI_ABORT
    integer       :: i1, i2, i3, islab12, dims(3)
    integer       :: owner_node
    real(kind=wp) :: verbose_actual_y, verbose_actual_z
    real(kind=wp) :: res_norm, error_norm, pot_norm
    real(kind=wp) :: defect_norm, defect_norm_stpt, defect_norm_bt, defect_norm_ord2
    integer       :: damping_iterations
    integer       :: damping_iterations_total
    real(kind=wp) :: s_damping ! damping factor
    integer       :: eq_type_defco
    character(len=*), parameter :: myself = 'dl_mg_defco_defect_corr_solver'
    character(len=*), parameter :: explanation = 'This could mean that the&
         & multigrid is too coarse to represent the quickly changing&
         & quantities (such as the permittivity).'

    ! Defect correction method variables
    real(kind=wp), allocatable :: defect(:,:,:)      ! current defect - total
    real(kind=wp), allocatable :: defect_stpt(:,:,:) ! - src and Poisson term
    real(kind=wp), allocatable :: defect_bt(:,:,:)   ! - Boltzmann term
    real(kind=wp), allocatable :: error(:,:,:)  ! hhh: curr. error approx.
    real(kind=wp), allocatable :: eps_full_1(:,:,:) ! jd: Eps of unity
    real(kind=wp), allocatable, target :: res_loc_data(:,:,:) ! workspace
    real(kind=wp), allocatable, target :: rho_loc_data(:,:,:)
    real(kind=wp), pointer     :: res_loc(:,:,:)
    real(kind=wp), pointer     :: rho_loc(:,:,:)
    logical                    :: converged
    logical                    :: test_conv_mu
    logical                    :: spokesman
    integer                    :: myid
    logical                    :: i_am_root
    logical                    :: do_defco           ! If .true., do defect correction
    logical                    :: is_damping_loc     ! Local version of optional is_damping
    logical                    :: use_mu             ! compute or display mu info if in PBE with PBC
    !------------------------------------------------------------------------

    ! JCW: Everything is awesome, unless someone tells us it isn't
    ierror = DL_MG_SUCCESS


    ! JCW: If mg_comm == MPI_COMM_NULL, then there is nothing to do on this rank
#ifdef MPI
    if ( mg_comm == MPI_COMM_NULL) then
       return
    endif
#endif

    call mg_timer(start, tsolver, tignore, tignore)

    myid = get_myid()

    ! JCW: Determine if current MPI process is root (id == 0)
    i_am_root = (myid == 0)

    ! JCW: Determine whether defect correction is needed
    if (fd_order == 2) then
       do_defco = .false.
    else if (fd_order > 2) then
       do_defco = .true.
       call set_fd(fd_order, ierr)
    else
      ! Unsupported fd_order /=2 and not > 2
      ! ==> Abort or return with error code
      ierr = DL_MG_ERR_FD_ORDER
      call handle_error(ierr, ierror, " failure in "//myself)
      return
    end if

    ! JCW: In the case of the defect correction (fd_order > 2) then der_pot
    ! JCW: should not have been passed, since der_pot is generated directly
    ! JCW: within the defect correction loop
    if ( do_defco .and. present(der_pot) ) then
       ierr = DL_MG_ERR_DEFCO_DERPOT
       call handle_error(ierr, ierror, " failure in "//myself)
       return
    end if

    ! JCW: Determine whether error damping will be used for the defect correction
    if (present(use_damping) ) then
       is_damping_loc = use_damping
    else
       is_damping_loc = .false.
    end if

    ! If not doing defect correction, stay quiet.
    ! As in dl_mg_info, only the root (myid == 0) MPI rank should output to
    ! the report unit.
    spokesman = i_am_root.and.do_defco

    if(spokesman) then
       ! Only output this information if defect correction is being attempted
       write(report_unit,'(a,i0,a)') 'DEFCO Getting initial solution using 2nd order solver'
       flush(report_unit)
    end if

    ! jd: Allocate workspace for residual
    if (present(res)) then
       res_loc => res
    else
       allocate(res_loc_data(isx:iex, isy:iey, isz:iez),stat=ierr)
       call dl_mg_defco_utils_alloc_check(ierr,myself,'res_loc')
       res_loc_data = 0.0_wp ! jd: Takes care of max-num, ld-pt padding
       res_loc => res_loc_data
    endif

    ! allocate mg data structure for the mulitgrid solver and
    ! populate the data arrays for the fields that do not change
    ! during the computation: eps and steric weight
    call set_mg_grids_eps_and_steric(eq_type, eps_half, steric_weight)

    use_mu = .false.
    if (fullpbc) then
       ! need to shift the source term when full pbc
       allocate(rho_loc_data(isx:iex, isy:iey, isz:iez))
       !rho_loc_data = rho
       !$omp parallel
       ! N.B. z_vtype is good for pbc as there is no boundary subgrid
       call blocked_scaled_copy2(mg(mg_levels), 1.0_wp, &
            rho, z_vtype, &
            rho_loc_data, z_vtype)
       !$omp end parallel
       rho_loc => rho_loc_data
       call remove_zero_mode_source(mg(mg_levels), eq_type, 1, &
            rho_loc, steric_weight)
       if (eq_type == EQ_PBE_NEWTON .or. &
            eq_type == EQ_PBE_NEWTON_NOSTERIC) then
          use_mu = .true.
       end if
    else
       rho_loc(isx:,isy:,isz:) => rho
       i1=size(mg)
       !write(0,*)'defco ', lbound(rho_loc), ubound(rho_loc), &
       !     blocked_dotproduct(mg(i1),&
       !     rho_loc(mg(i1)%sx:,mg(i1)%sy:,mg(i1)%sz:), &
       !     z_vtype, rho_loc(mg(i1)%sx:,mg(i1)%sy:,mg(i1)%sz:), z_vtype)
    endif
    
    ! ************************************************************************
    ! *** CALL THE SOLVER (FOR PBE EQUATION PROPER)
    ! ************************************************************************
    call dlmg_2nd_order_solver(eq_type, alpha, rho_loc, &
         pot, use_pot_in, res=res_loc, der_pot=der_pot, &
         ierror=ierr)
     if (ierr /= DL_MG_SUCCESS) then
        call handle_error(ierr, ierror, " failure in first call in defco")
        return ! relevant only if abort_on_error is false
     endif
    ! ************************************************************************

    ! ************************************************************************
    ! *** DEFECT CORRECTION PROCEDURE
    ! ************************************************************************
    ! jd: If higher order was requested, employ defect correction
    if (.not.do_defco) then
       ! 2nd order solution requested, do not perform defect correction
       ! ==> Return ierror from dl_mg_2nd_order_solver
       if (present(ierror)) ierror = ierr
    else
       ! >2nd order solution requested, perform defect correction

       damping_iterations_total = 0

       ! Allocate memory for local arrays needed for defect correction
       allocate(defect(isx:iex,isy:iey,isz:iez),stat=ierr)
       call dl_mg_defco_utils_alloc_check(ierr,myself,'defect')
       allocate(defect_stpt(isx:iex,isy:iey,isz:iez),stat=ierr)
       call dl_mg_defco_utils_alloc_check(ierr,myself,'defect_stpt')
       allocate(defect_bt(isx:iex,isy:iey,isz:iez),stat=ierr)
       call dl_mg_defco_utils_alloc_check(ierr,myself,'defect_bt')
       allocate(error(isx:iex,isy:iey,isz:iez),stat=ierr)
       call dl_mg_defco_utils_alloc_check(ierr,myself,'error')

       ! hhh: Error initialized to zero -this is fine for possibly arbitrary BCs
       !      set earlier by optional bound argument in multigrid_poisson_solver
       error(:,:,:) = 0.0_wp

       ! hhh: We now have all the inputs ready to start the defect correction
       !      iteration
       converged = .false.
       defco_status%iter = 0

       if(spokesman) then
          write(report_unit,'(a,i0)') 'DEFCO Starting defect correction with &
               &a discretization order of: ', fd_order
          if (use_mu) then
             ! in this case we print mu as well
                write(report_unit,'(a)') "DEFCO iter |ST+PT defect|   |BT defect|  &
                     &|tot defect|   | pot |        |mu|           &
                     &|error|       |err.eqn res.|"
             else
                write(report_unit,'(a)') "DEFCO iter |ST+PT defect|   |BT defect|  &
                     &|tot defect|   | pot |        |error|       |err.eqn res.|"
             end if
             flush(report_unit)
       end if


       ! **********************************************************************
       ! *** DEFECT CORRECTION LOOP
       ! **********************************************************************
       do while (.not. converged)

          defco_status%iter = defco_status%iter + 1

          ! hhh: Compute the defect using high order finite difference
          !      Level of accuracy is selected through the input param
          !      discretization_order
          !      defect = source + \nabla \cdot [\epsilon \grad potmg]
          call compute_current_defect(fd_order, eq_type, alpha, &
               rho_loc, eps_full, pot, &
               defect, defect_norm, defect_stpt, defect_norm_stpt, &
               defect_bt, defect_norm_bt, steric_weight=steric_weight )

          if (defco_status%iter == 1) then
             defect_norm_ord2 = defect_norm
             defco_status%target_res_norm = &
                  max(tol_res_abs_, tol_res_rel_ * defect_norm_ord2)
          end if

!         call internal_multigrid_verbose_input_defco

         ! hhh: Now that we have the high order defect we very approximately
         !      solve the error equation using second order multigrid method
         !      with a low convergence threshold
         !      - \nabla \cdot [epsilon \grad error] = defect
         !defect = defect/(FOUR_PI) ! jd: (*) solver expects 'rho', not 'source'
         !error = pot

         ! *****************************************************************
         ! *** CALL THE SOLVER (DEFECT EQUATION)
         ! *****************************************************************
         ! JCW: Error should be zeroed --> zero BCs for solving defect equation
         error = 0.0_wp

         ! needs to solve linearised PBE at pot
         select case(eq_type)
         case(EQ_PBE_NEWTON, EQ_PBE_NEWTON_NOSTERIC, EQ_PBE_FAS, EQ_PBE_FAS_NOSTERIC)
            eq_type_defco = EQ_LINEAR_PBE_POT
         case default
            eq_type_defco = eq_type
         end select

         call remove_zero_mode_source(mg(mg_levels), eq_type_defco, 2, &
            defect)
         call dlmg_2nd_order_solver(eq_type_defco, -1.0_wp, defect, &
              error, use_pot_in, remove_source_zero_mode=2, &
              res=res_loc, der_pot=pot, &
              ierror=ierr)
         if (ierr /= DL_MG_SUCCESS) then
            call handle_error(ierr, ierror)
            exit ! relevant only when abort_on_eror = false
         endif
         ! *****************************************************************

         ! jd: Compute error and error equation residual norms
         error_norm = sqrt(dl_mg_defco_utils_integrate_product_on_grid( &
              grid_weight,error,error,mx,my,mz,comm=fd%comm))
         defco_status%error_norm = error_norm

         res_norm = sqrt(dl_mg_defco_utils_integrate_product_on_grid( &
              grid_weight,res_loc,res_loc,mx,my,mz,comm=fd%comm))
         defco_status%res_norm = res_norm
         ! save some time, res_norm is computed in dlmg_2nd_order_solver
         !res_norm = info_get_res_norm("mg") * sqrt(grid_weight)
         pot_norm = sqrt(dl_mg_defco_utils_integrate_product_on_grid( &
              grid_weight,pot,pot,mx,my,mz,comm=fd%comm))
         defco_status%sol_norm = pot_norm

         if(spokesman) then
            if (use_mu) then
               write(report_unit,'(a,i4,7(a,e13.7))') &
                    "DEFCO ", defco_status%iter,"  ", &
                    defect_norm_stpt," ",defect_norm_bt," ",defect_norm," ", &
                    pot_norm, " ", sum(abs(mu))," -> ", error_norm," ",res_norm
            else
               write(report_unit,'(a,i4,6(a,e13.7))') &
                    "DEFCO ", defco_status%iter,"  ", &
                    defect_norm_stpt," ",defect_norm_bt," ",defect_norm," ", &
                    pot_norm," -> ", error_norm," ",res_norm
            end if
            flush(report_unit)
         end if

         ! la: damp the error if needed
         damping_iterations = 0
         if (is_damping_loc) then
            call internal_error_damping(error, s_damping, pot, defect_norm, ierr, steric_weight)
            if (ierr /= DL_MG_SUCCESS) then
               call handle_error(ierr,ierror)
               exit
            end if
            damping_iterations_total = &
                 damping_iterations_total + damping_iterations
         end if

         ! Correct the current approximation to get the new approximation
         pot(1:mx,1:my,1:mz) = pot(1:mx,1:my,1:mz) + error(isx:iex,isy:iey,isz:iez)

         ! for full PBC and defco correction to newton
         if (use_mu) then
            mu_prev = mu
            call compute_mu_defco(pot, steric_weight)
         end if

         if ( damping_iterations > 1) then
            error_norm = sqrt( s_damping ) * error_norm
         endif

         ! jd: Sanity checks for runaway solutions
         if (error_norm > 1D6) then
            call handle_error(DL_MG_ERR_DEFCO_UNPHYSICAL,ierror,&
                 msg='The error norm is absurdly high (>1.0e6). '//trim(explanation))
            exit
         end if
         if (defect_norm > 1D6) then
            call handle_error(DL_MG_ERR_DEFCO_UNPHYSICAL,ierror,&
                 msg='The defect norm is absurdly high (>1.0e6). '//trim(explanation))
            exit
         end if
         if (res_norm > 1D6) then
            call handle_error(DL_MG_ERR_DEFCO_UNPHYSICAL,ierror,&
                 msg='The residual norm is absurdly high (>1.0e6). '//trim(explanation))
            exit
         end if

         if (defco_status%iter > defco_maxiters) then
            ! Maximum number of defect-correction iterations exceeded, handle error
            ! and exit loop (if handle_error returns, rather than aborting execution)
            call handle_error(DL_MG_ERR_DEFCO_ITER,ierror,msg=explanation)
            exit
         end if

         if (use_mu) then
            test_conv_mu = all( abs(mu_prev-mu) < max(tol_mu_rel_ * abs(mu), tol_mu_abs_))
         else
            test_conv_mu = .true.
         end if
         if (error_norm  < max(tol_pot_abs_, tol_pot_rel_ * pot_norm) .and. &
              defect_norm < defco_status%target_res_norm .and. &
              test_conv_mu ) &
              converged = .true.

      end do
       ! **********************************************************************
       ! *** END OF DEFECT CORRECTION LOOP
       ! **********************************************************************

       if (use_mu) then
          call shift_pot_and_mu_defco(eq_type_defco, pot, steric_weight)
       endif

       defco_status%ierror = ierror

       if(spokesman) then
          write(report_unit,*)
          if (converged) then
             write(report_unit,'(a)') 'DEFCO Defect corrected solution obtained'
          else
             write(report_unit,'(a)') 'DEFCO Defect correction loop failed: last iteration data:'
          end if
          write(report_unit,'(a,e12.4)') 'DEFCO  - final defect norm:            ', &
               defect_norm
          write(report_unit,'(a,e12.4)') 'DEFCO  - final error norm:             ', &
               error_norm
          write(report_unit,'(a,i0)')    'DEFCO  - defect-correction iterations: ', &
               defco_status%iter
          if(is_damping_loc) then
             write(report_unit,'(a,i0)') 'DEFCO  - error-damping iterations:     ', &
                  damping_iterations_total
          end if
          write(report_unit,*)
          flush(report_unit)
       end if

       ! jd: Clean up
       deallocate(error,stat=ierr)
       call dl_mg_defco_utils_dealloc_check(ierr,myself,'error')
       deallocate(defect_bt,stat=ierr)
       call dl_mg_defco_utils_dealloc_check(ierr,myself,'defect_bt')
       deallocate(defect_stpt,stat=ierr)
       call dl_mg_defco_utils_dealloc_check(ierr,myself,'defect_stpt')
       deallocate(defect,stat=ierr)
       call dl_mg_defco_utils_dealloc_check(ierr,myself,'defect')
    end if ! if defect correction used

    ! free mg memory
    call free_mg
    
    call mg_timer(stop, tsolver, tignore, tignore)
    call report_timings

  contains

    subroutine internal_error_damping(error, s, pot, defect_norm, ierror, steric_weight )
      ! -----------------------------------------------------------------------
      ! Error damping facility for full PBE, by Lucian Anton
      ! Modified for ONETEP conventions by Jacek Dziedzic
      ! -----------------------------------------------------------------------

      use dl_mg_common_data, only: mx, my, mz
      use dl_mg_defco_utils, only: dl_mg_defco_utils_alloc_check, dl_mg_defco_utils_assert, dl_mg_defco_utils_dealloc_check
      use dl_mg_errors
      implicit none

      real(kind=wp), intent(inout) :: error(:,:,:)
      real(kind=wp), intent(out)   :: s ! damping factor , needed out for norm rescaling
      real(kind=wp), intent(in)    :: pot(:,:,:)
      real(kind=wp)                :: defect_norm
      integer, intent(out)         :: ierror
      real(wp), intent(in), optional :: steric_weight(:,:,:)

      ! jd: Internal variables
      real(kind=wp), parameter     :: s_min = 1.d-5   ! minimum damping
      real(kind=wp), parameter     :: q_damp =0.75_wp ! damping factor
      integer                      :: l1, l2, l3
      real(kind=wp)                :: defect_norm_aux
      real(kind=wp)                :: defect_norm_aux_stpt
      real(kind=wp)                :: defect_norm_aux_bt
      logical                      :: found
      real(kind=wp), allocatable   :: pot_aux(:,:,:), defect_aux(:,:,:,:)
      integer                      :: ierr
      character(len=*), parameter  :: myself = 'internal_error_damping'

      ! ----------------------------------------------------------------------

      ierror = DL_MG_SUCCESS

      s = 1.0_wp
      damping_iterations = 0
      found = .false.

      ! JCW: Get extents of input array
      l1 = size(pot,1); l2 = size(pot,2); l3 = size(pot,3)
      call dl_mg_defco_utils_assert(&
           l1>=mx.and.l2>=my.and.l3>=mz,&
           "Error in "//trim(myself)//": &
           &extent of pot array is inconsistent with size of data on this MPI rank.")

      allocate(pot_aux(mx,my,mz), stat=ierr)
      call dl_mg_defco_utils_alloc_check(ierr,myself,'pot_aux')
      allocate(defect_aux(mx,my,mz,3), stat=ierr)
      call dl_mg_defco_utils_alloc_check(ierr,myself,'defect_aux')

      do while ( s > s_min )

         pot_aux(1:mx,1:my,1:mz) = pot(1:mx,1:my,1:mz) + s * error(1:mx,1:my,1:mz)

         call compute_current_defect(fd_order, eq_type, alpha, rho_loc, &
              eps_full, pot_aux, &
              defect_aux(:,:,:,1), defect_norm_aux, &
              defect_aux(:,:,:,2), defect_norm_aux_stpt, &
              defect_aux(:,:,:,3), defect_norm_aux_bt, &
              steric_weight=steric_weight )

         damping_iterations = damping_iterations + 1

         if(spokesman) then
            write(report_unit,'(a,i4,a,e13.7,a,e13.7,a,e13.7,a,i2,f11.8,a)') &
                 "DEFCO ", defco_status%iter,"  ", &
                 defect_norm_aux_stpt," ",defect_norm_aux_bt," ",&
                 defect_norm_aux,"        DAMPING (",damping_iterations,s,")"
         end if

         if (defect_norm_aux < defect_norm) then
            error(1:mx,1:my,1:mz) = s * error(1:mx,1:my,1:mz)
            found = .true.
            exit
         else
            s = q_damp * s
         endif

      end do

      if (.not.found) then
         ! Error code to return, if error damping was not successful
         ierror = DL_MG_ERR_DEFCO_DAMPING
      end if

      deallocate(defect_aux, stat=ierr)
      call dl_mg_defco_utils_dealloc_check(ierr,myself,'defect_aux')
      deallocate(pot_aux, stat=ierr)
      call dl_mg_defco_utils_dealloc_check(ierr,myself,'pot_aux')

    end subroutine internal_error_damping


#if 0

    subroutine get_eps_half_from_eps_full
      ! -----------------------------------------------------------------------
      ! computes eps_half from eps_full by linear interpolation from the nearest
      ! grid points
      ! used a grid type to take advantage of halo exchange subroutines
      ! form dl_mg_mpi
      ! A cubic interpolation should be tried  instead of linear
      ! -----------------------------------------------------------------------
      use dl_mg_types, only : mg_t
      use dl_mg_common_data, only : mg_levels, mg_comm
      use dl_mg_mpi
      implicit none

      type(mg_t) g
      integer i, j, k, l1, l2, l3

      l1 = size(rho,1); l2 = size(rho,2); l3 = size(rho,3)
      g%comm = mg_comm
      g%level = mg_levels
      g%active = .true.
      g%fsx = 0
      g%sx  = 1
      g%fex = l1 +1
      g%ex  = l1
      g%fsy = 0
      g%sy  = 1
      g%fey = l2 +1
      g%ey  = l2
      g%fsz = 0
      g%sz  = 1
      g%fez = l3 +1
      g%ez  = l3
      allocate(g%f(0:l1+1,0:l2+1,0:l3+1))
      g%f(1:l1,1:l2,1:l3) = eps_full(:,:,:)
      call find_mpi_neighbors(g)
      do i =1, 3
         ! each rank needs only data from right and it has to send to left ranks
         call exchange_halo_begin(i, exch_sl+exch_rr, 1, g,"f")
         call exchange_halo_end(i, exch_sl+exch_rr, 1, g,"f")
      enddo

      do k = 1, l3
         do j = 1, l2
            do i = 1, l1
               eps_half_loc(i,j,k,1) = 0.5_wp*(g%f(i,j,k) &
                    + g%f(i+1,j,k))
               eps_half_loc(i,j,k,2) = 0.5_wp*(g%f(i,j,k) &
                    + g%f(i,j+1,k))
               eps_half_loc(i,j,k,3) = 0.5_wp*(g%f(i,j,k) &
                    + g%f(i,j,k+1))
            end do
         end do
      end do

    end subroutine get_eps_half_from_eps_full
#endif
  end subroutine dl_mg_defco_defect_corr_solver


  subroutine compute_current_defect(fd_order, eq_type, alpha, rho, eps, v_i,  & ! in
    defect, defect_norm, defect_stpt, defect_norm_stpt, &         ! out
    defect_bt, defect_norm_bt, &                                  ! out
    steric_weight )                                               ! optional, in
    !=========================================================================!
    ! @updateme
    ! This subroutine will compute the current defect for the Poisson equation!
    ! in a dielectric medium.  The defect is defined as the amount that the   !
    ! current approximation fails to satisfy the equation or in other words:  !
    !                                                                         !
    !    defect = source + \nabla \cdot [eps \grad v_i],                      !
    ! where source = - alpha \rho                                             !
    !                                                                         !
    ! The discretization of the operators is selected by the                  !
    ! discretization_order parameter.  This subroutine is intended for use    !
    ! within the defect correction method which will produce a higher         !
    ! accuracy solution than is possible with normal multigrid iteration.     !
    !                                                                         !
    ! JCW, 30/01/17:                                                          !
    ! LA, 4/17:                                                               !
    ! If optional argument steric_weight is not present, then the defect      !
    ! is evaluated assuming steric weight is 1.0  for all space               !
    ! (i.e. zero steric potential).                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !    order       (input): finite difference order to use                  !
    !    rho         (input): fixed charge density                            !
    !    eps         (input): the diel. functional on the full grid points    !
    !    v_i         (input): the current approximation of the potential      !
    !    grid        (input): the defect is computed on the mg grid, but the  !
    !                         host grid is needed here and there.             !
    !    defect      (output): defined as above                               !
    !    defect_norm (output): discrete L2 norm of the defect function        !
    !    steric_weight (input): optional, steric weight to be used in PB      !
    !                         equation                                        !
    ! Arrays are assumed to be in the distributed slab representation.        !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal 30/04/2010                                     !
    ! Mended for onetep by HH Helal 7/6/2010                                  !
    ! Adapted for parallel operation by Jacek Dziedzic, 10-11/06/2010.        !
    !=========================================================================!

    use dl_mg_params, only : EQ_POISSON, EQ_LINEAR_PBE, EQ_LINEAR_PBE_NOSTERIC, EQ_LINEAR_PBE_POT,&
         EQ_PBE_FAS, EQ_PBE_FAS_NOSTERIC, EQ_PBE_NEWTON, EQ_PBE_NEWTON_NOSTERIC, &
         DL_MG_BC_PERIODIC
    use dl_mg_common_data, only: grid_weight, isx, iex, isy, iey, isz, iez, mx, my, mz, &
         nx, ny, nz, bc, fullpbc
    use dl_mg_grids, only: fd
    use dl_mg_defco_utils, only: dl_mg_defco_utils_alloc_check, dl_mg_defco_utils_dealloc_check, &
         dl_mg_defco_utils_assert, dl_mg_defco_utils_integrate_product_on_grid
    use dl_mg_defco_fd, only: dl_mg_defco_fd_gradient, dl_mg_defco_fd_laplacian
    use dl_mg_nonlinear_model, only : neutralisation_method

    implicit none

    ! jd: Arguments
    integer,       intent(in)   :: fd_order, eq_type
    real(kind=wp), intent(in)   :: alpha
    real(kind=wp), intent(in)   :: rho(1:,1:,1:)
    real(kind=wp), intent(in)   :: eps(1:,1:,1:)
    real(kind=wp), intent(in)   :: v_i(1:,1:,1:)
    real(kind=wp), intent(out)  :: defect(1:,1:,1:)
    real(kind=wp), intent(out)  :: defect_norm
    real(kind=wp), intent(out)  :: defect_stpt(1:,1:,1:)
    real(kind=wp), intent(out)  :: defect_norm_stpt
    real(kind=wp), intent(out)  :: defect_bt(1:,1:,1:)
    real(kind=wp), intent(out)  :: defect_norm_bt
    real(kind=wp), optional, intent(in) :: steric_weight(1:,1:,1:)

    ! jd: Internal variables
    real(kind=wp), allocatable     :: lap_v(:,:,:)         ! nabla^2 v
    real(kind=wp), allocatable     :: grad_eps(:,:,:,:)    ! nabla eps
    real(kind=wp), allocatable     :: grad_v(:,:,:,:)      ! nabla v
    integer                        :: ierr                 ! jd: Error flag
    logical                        :: loc_use_linear
    character(len=*), parameter    :: myself = 'compute_current_defect'

    real(kind=wp) :: defect_av ! JCW: Average of defect over entire grid, for
                               ! JCW: subtracting k=0 mode in PBCs

    !------------------------------------------------------------------------

    ! JCW: Sanity check the dimensions:
    ! JCW: Dimensions 1, 2, 3 should be greater than or equal to dimensions
    ! JCW: of arrays held on this MPI rank (as set by dl_mg_init)
    ! JCW: [ The dimensions of the defect_ arrays are mx, my, mz ]
    call dl_mg_defco_utils_assert(&
         size(rho,1)>=mx.and.size(rho,2)>=my.and.size(rho,3)>=mz,&
         "Error in "//trim(myself)//": &
         &extent of rho array is inconsistent with size of data on this MPI rank.")
    call dl_mg_defco_utils_assert(&
         size(eps,1)>=mx.and.size(eps,2)>=my.and.size(eps,3)>=mz,&
         "Error in "//trim(myself)//": &
         &extent of eps array is inconsistent with size of data on this MPI rank.")
    call dl_mg_defco_utils_assert(&
         size(v_i,1)>=mx.and.size(v_i,2)>=my.and.size(v_i,3)>=mz,&
         "Error in "//trim(myself)//": &
         &extent of v_i array is inconsistent with size of data on this MPI rank.")

    ! Allocate memory for local variables
    allocate(lap_v(mx,my,mz),stat=ierr)
    call dl_mg_defco_utils_alloc_check(ierr,myself,'lap_v')
    allocate(grad_eps(mx,my,mz,3),stat=ierr)
    call dl_mg_defco_utils_alloc_check(ierr,myself,'grad_eps')
    allocate(grad_v(mx,my,mz,3),stat=ierr)
    call dl_mg_defco_utils_alloc_check(ierr,myself,'grad_v')

    ! JCW: Use DL_MG-specific high-order FD routines
    !
    ! IMPORANT NOTE: Care should be taken to ensure that the input array
    !                for dl_mg_defco_fd_{gradient,laplacian} is contiguous.
    !                If non-contiguous, this can cause problems with non-blocking
    !                MPI calls downstream.
    call dl_mg_defco_fd_gradient(grad_v(1:mx,1:my,1:mz,1:3),v_i,fd_order)
    call dl_mg_defco_fd_laplacian(lap_v(1:mx,1:my,1:mz),v_i,fd_order)
    call dl_mg_defco_fd_gradient(grad_eps(1:mx,1:my,1:mz,1:3),eps,fd_order)

    ! Now have everything we need to compute the current defect and defect_norm
    defect_stpt(1:mx,1:my,1:mz) = -alpha * rho(1:mx,1:my,1:mz) + &
         eps(1:mx,1:my,1:mz) * lap_v(1:mx,1:my,1:mz) + &
         grad_eps(1:mx,1:my,1:mz,1) * grad_v(1:mx,1:my,1:mz,1) + &
         grad_eps(1:mx,1:my,1:mz,2) * grad_v(1:mx,1:my,1:mz,2) + &
         grad_eps(1:mx,1:my,1:mz,3) * grad_v(1:mx,1:my,1:mz,3)

!!$    if ( fullpbc ) then !.and. (eq_type == EQ_POISSON .or. &
!!$         !(eq_type /= EQ_POISSON .and. .not. use_counterions))) then
!!$       defect_av = dl_mg_defco_utils_integrate_product_on_grid(1.0_wp/(real(nx,wp)*ny*nz), &
!!$            defect_stpt, m1=mx, m2=my, m3=mz,comm=fd%comm)
!!$       defect_stpt(1:mx,1:my,1:mz) = defect_stpt(1:mx,1:my,1:mz) - defect_av
!!$    endif

    defect_norm_stpt = sqrt(dl_mg_defco_utils_integrate_product_on_grid( &
         grid_weight,defect_stpt,defect_stpt,mx,my,mz,comm=fd%comm) )

    if (eq_type /= EQ_POISSON) then

       ! JCW: Set loc_use_linear
       if ( eq_type == EQ_PBE_FAS .or.           &
            eq_type == EQ_PBE_FAS_NOSTERIC .or.  &
            eq_type == EQ_PBE_NEWTON .or.        &
            eq_type == EQ_PBE_NEWTON_NOSTERIC) then
          ! JCW: Full non-linear equation
          loc_use_linear = .false.
       else
          ! JCW: Linearised equation
          loc_use_linear = .true.
       end if

       call compute_boltzmann_defect_term(defect_bt, v_i, &
            mx, my, mz,loc_use_linear,steric_weight)

!!$       if (fullpbc ) then !.and. .not. use_counterions) then
!!$          !jellium
!!$          defect_av = dl_mg_defco_utils_integrate_product_on_grid(1.0_wp/(real(nx,wp)*ny*nz), &
!!$               defect_bt, m1=mx, m2=my, m3=mz,comm=fd%comm)
!!$          defect_bt(1:mx,1:my,1:mz) = defect_bt(1:mx,1:my,1:mz) - defect_av
!!$       endif

       defect_norm_bt = sqrt(dl_mg_defco_utils_integrate_product_on_grid( &
            grid_weight,defect_bt,defect_bt,mx,my,mz,comm=fd%comm))
       defect = defect_stpt - defect_bt
       defect_norm = sqrt(dl_mg_defco_utils_integrate_product_on_grid( &
            grid_weight,defect,defect,mx,my,mz,comm=fd%comm))
    else
       defect_bt(1:mx,1:my,1:mz) = 0.0_wp
       defect_norm_bt = 0.0_wp
       defect = defect_stpt
       defect_norm = defect_norm_stpt
    end if

    ! jd: Clean up
    deallocate(grad_v,stat=ierr)
    call dl_mg_defco_utils_dealloc_check(ierr,myself,'grad_v')
    deallocate(grad_eps,stat=ierr)
    call dl_mg_defco_utils_dealloc_check(ierr,myself,'grad_eps')
    deallocate(lap_v,stat=ierr)
    call dl_mg_defco_utils_dealloc_check(ierr,myself,'lap_v')

  end subroutine compute_current_defect


  !> Compute the Boltzmann term to add to the defect for the Poisson equation
  !!
  !! For the full non-linear PB equation the additional term in the defect is
  !! \f[ -\lambda \sum_{i} c_{i} q_{i} \beta W \exp( -\beta q_{i} \phi^{(i)})\f]
  !! and for the linearized PB equation, the additional term is
  !!  \f[ \lambda {\lambda_{D}^{2}} W \phi^{(i)} \f]
  !! where \f$ c_{i}\f$ and \f$ q_{i} \f$ are the concentration and charge for ion type \f$i\f$,
  !! \f$\beta = 1/kT\f$, \f$ W\f$ is the steric weight \f$\exp(-\beta V)]\f$,
  !! and \f$\lambda_{D}\f$ is
  !! the Debye length.
  !!
  !! If optional argument steric_weight is present, then the defect is evaluated
  !! using the steric potential, if not, then the steric potential is assumed
  !! to be 0.0 everywhere.
  !!
  !! Adapted from implicit_solvent_boltzmann_defect from ONETEP's
  !! is_solvation_boltzmann module for inclusion in DL_MG. See
  !! original documentation in source.
  !!
  !! J. C. Womack, 2016
  subroutine compute_boltzmann_defect_term(bt, pot, n1, n2, n3, use_linear, steric_weight)
    !=========================================================================!
    ! Computes the Boltzmann term to defect correction.                       !
    !=========================================================================!
    use dl_mg_defco_utils, only: dl_mg_defco_utils_assert
    use dl_mg_nonlinear_model, only: temp, nion, q => qion, &
         c => cion, lambda, &
         capexp, kappa2, mu, has_steric_weight, beta, &
         steric_weight_sum, rho_sum
    use dl_mg_grids, only : mg
    use dl_mg_common_data, only : isx, isy, isz, &
         iex, iey, iez, t => mg_levels
    ! Variables imported from dl_mg_nonlinear_model:
    ! temp:        temperature, in Kelvin
    ! nion:           number of ionic species
    ! c(1:n):   array containing concentrations of ionic species in atomic
    !              units, i.e. 1/a_{0}^{-3}
    ! q(1:n):   array containing charges of ionic species in atomic units
    ! lambda_o_beta: multiplicative constant for non-linear Boltzmann term
    !              (-4.0*\pi for atomic units)
    ! lambda:      multiplicative constant for linear Boltzmann term, defined as
    !              z / kT, where z is a multiplicative constant based on unit
    !              system (-4.0*\pi for atomic units).
    ! expcap:      cap of the exponential argument in Boltzmann term
    ! max_expcap:  pre-defined maximum for exponential cap
    ! kappa2:      kappa2 = -lamda * sum ( cq(1, :) * cq(2, :)**2 )
    !              i.e. (1/(Debye length))^{2}, where rel permittivity is 1.0

    implicit none

    ! jd: Arguments
    real(kind=wp), intent(out) :: bt(isx:, isy:, isz:)
    real(kind=wp), intent(in)  :: pot(isx:, isy:, isz:)
    integer, intent(in)        :: n1, n2, n3
    logical, intent(in)        :: use_linear
    real(kind=wp), optional, intent(in) :: steric_weight(isx:, isy:, isz:)

    ! jd: Internal variables
    integer                    :: i, j, k, ion
    real(kind=wp)              :: x
    character(len=*),parameter :: myself = "compute_boltzmann_defect_term"

    ! -----------------------------------------------------------------------



    !$OMP PARALLEL DEFAULT(none) SHARED(mg, bt, has_steric_weight, t, pot, &
    !$OMP steric_weight, beta,  c, q,  nion, &
    !$OMP rho_sum, steric_weight_sum, &
    !$OMP isx, isy, isz, iex, iey, iez, use_linear, mu, lambda, kappa2) &
    !$OMP PRIVATE(i, j, k, x, ion)

    if (.not.use_linear) then
       ! Full non-linear PB equation
       !$OMP WORKSHARE
       bt(:, :, :) = 0.0_wp
       !$OMP END WORKSHARE
       if (has_steric_weight) then
          do ion = 1, nion
             !$OMP DO COLLAPSE(2)
             do k = isz, iez
                do j = isy, iey
                   do i = isx, iex
                      x =  -beta * q(ion) * pot(i,j,k) +mu(ion)
                      bt(i,j,k) = bt(i,j,k) + &
                           lambda * c(ion) * q(ion) * &
                           steric_weight(i, j, k) * capexp(x)
                   enddo
                enddo
             enddo
             !$OMP ENDDO NOWAIT
          end do
       else
          do ion = 1, nion
             !$OMP DO COLLAPSE(2)
             do k = isz, iez
                do j = isy, iey
                   do i = isx, iex
                      x = -beta * q(ion) * pot(i,j,k) +mu(ion)
                      bt(i,j,k) = bt(i,j,k) + &
                           lambda * c(ion) * q(ion) *  capexp(x)
                   enddo
                enddo
             enddo
             !$OMP ENDDO NOWAIT
          end do
       endif
    else
       ! Linearized PB equation
       if (has_steric_weight) then
          !bt(1:n1, 1:n2, 1:n3) = kappa2 * &
          !     steric_weight(1:n1, 1:n2, 1:n3) * pot(1:n1, 1:n2, 1:n3)
          !$OMP DO COLLAPSE(2)
          do k = isz, iez
             do j = isy, iey
                do i = isx, iex
                   bt(i,j,k) = kappa2 * steric_weight(i, j, k) * pot(i,j,k)
                enddo
             enddo
          enddo
          !$OMP ENDDO NOWAIT
       else
          !bt(1:n1, 1:n2, 1:n3) = kappa2 * pot(1:n1, 1:n2, 1:n3)
          !$OMP DO COLLAPSE(2)
          do k = isz, iez
             do j = isy, iey
                do i = isx, iex
                   bt(i,j,k) = kappa2 * pot(i,j,k)
                enddo
             enddo
          enddo
          !$OMP ENDDO NOWAIT
       endif
    end if
    !$OMP END PARALLEL

  end subroutine compute_boltzmann_defect_term


  !> computes the correction to chemical pontential in defco loop
  !! relevant when using PBE with PBC
  subroutine compute_mu_defco(pot, steric_weight)
    use dl_mg_params, only : hartree
    use dl_mg_nonlinear_model, only : capexp, c => cion, &
         q => qion, mu, &
         has_steric_weight, steric_weight_sum, nion, temp, &
         beta
    use dl_mg_common_data, only : isx, isy, isz, blk, t => mg_levels
    use dl_mg_grids, only : mg
    use dl_mg_types, only : mg_get_idx
    use dl_mg_mpi
    implicit none


    real(wp), intent(inout) :: pot(isx:, isy:, isz:)
    real(wp), intent(in), optional :: steric_weight(isx:, isy:, isz:)

    integer i,j,k, ion_type, iblk, ierr
    integer, pointer :: s(:), e(:)
    real(wp) x, stw
    real(wp), allocatable :: lsum(:), w(:)


    allocate(lsum(nion), w(nion))

    lsum = 0.0_wp

    !$OMP PARALLEL DEFAULT(none) SHARED(has_steric_weight, blk, t, pot, &
    !$OMP steric_weight, beta, lsum, c, q, nion) &
    !$OMP PRIVATE(i,j,k,iblk, x, stw, s, e, ion_type)
    if (has_steric_weight) then
       !$OMP DO SCHEDULE(STATIC,1) REDUCTION(+:lsum)
       do iblk = 1, size(blk(t)%start, dim=2)
          s => blk(t)%start(:,iblk)
          e => blk(t)%end(:,iblk)
          do k = s(3), e(3)
             do j = s(2), e(2)
                do i = s(1), e(1)
                   stw = steric_weight(i,j,k)
                   do ion_type = 1, nion
                      x = -beta * q(ion_type) * pot(i,j,k)
                      lsum(ion_type) = lsum(ion_type) + stw * capexp(x)
                   end do
                end do
             end do
          end do
       end do
       !$OMP ENDDO NOWAIT
    else
       !$OMP DO SCHEDULE(STATIC,1) REDUCTION(+:lsum)
       do iblk = 1, size(blk(t)%start, dim=2)
          s => blk(t)%start(:,iblk)
          e => blk(t)%end(:, iblk)
          do k = s(3), e(3)
             do j = s(2), e(2)
                do i = s(1), e(1)
                   do ion_type = 1, nion
                      x = -beta * q(ion_type) *pot(i,j,k)
                      lsum(ion_type) = lsum(ion_type) + capexp(x)
                   enddo
                end do
             end do
          end do
       end do
       !$OMP ENDDO NOWAIT
    end if
    !$OMP END PARALLEL

    ! reduce accros MPI ranks
#ifdef MPI
    call MPI_ALLREDUCE(lsum, w, nion, MPI_DOUBLE_PRECISION, MPI_SUM, mg(t)%comm, ierr)
#else
    w(:) = lsum(:)
#endif

    mu(:) = - log (w(:) / steric_weight_sum)

    deallocate(lsum, w)

  end subroutine compute_mu_defco


  subroutine shift_pot_and_mu_defco(eq_type_defco, pot, steric_weight)
    use dl_mg_nonlinear_model, only : has_steric_weight, steric_weight_sum, &
         beta, q => qion, mu, neutralisation_method
    use dl_mg_common_data, only : isx, isy, isz, blk, t => mg_levels, &
         fullpbc
    use dl_mg_params, only : EQ_LINEAR_PBE_POT, dl_mg_neutralise_with_ions_auto
    use dl_mg_grids, only : mg
    use dl_mg_types, only : mg_get_idx
    use dl_mg_mpi
    use dl_mg_kernels
    use dl_mg_neutralisation_with_ions, only : compute_shifted_ion_concentrations
    implicit none

    integer, intent(in) :: eq_type_defco
    real(wp), intent(inout) :: pot(isx:, isy:, isz:)
    real(wp), intent(in), optional :: steric_weight(isx:, isy:, isz:)

    integer sx, sy, sz, ex, ey, ez
    real(wp) pave, s1, x, cs


    call mg_get_idx (mg(t), sx, ex, sy, ey, sz, ez)

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(s1, pave)
    if (has_steric_weight) then
       s1 = blocked_dotproduct(mg(t), pot(sx:,sy:,sz:), z_vtype, &
            steric_weight(sx:,sy:,sz:), z_vtype)
    else
       s1 = blocked_dotproduct(mg(t), pot(sx:,sy:,sz:), z_vtype)
    end if

    pave = s1/steric_weight_sum
    call blocked_update2(mg(t), pot(sx:,sy:,sz:), z_vtype, -pave)

    ! shift mu_i -> mu_i - beta * q_i s1
    !$OMP MASTER
    mu(:) = mu(:) - beta * q(:) * pave

    ! get the new concentrations after the mu was computed
    if (neutralisation_method == dl_mg_neutralise_with_ions_auto &
       .or. &
        neutralisation_method == dl_mg_neutralise_with_ions_auto_linear) then
       call compute_shifted_ion_concentrations
    end if
    !$OMP END MASTER
    !$OMP END PARALLEL

  end subroutine shift_pot_and_mu_defco


  subroutine set_mg_grids_eps_and_steric(eq_type, eps_half, steric_weight)
    use dl_mg_params, only : DL_MG_BC_DIRICHLET, &
         EQ_PBE_FAS, EQ_PBE_FAS_NOSTERIC, &
         EQ_PBE_NEWTON, EQ_PBE_NEWTON_NOSTERIC, &
         EQ_LINEAR_PBE_POT
    use dl_mg_common_data, only : t => mg_levels, bc, isx, isy, isz
    use dl_mg_grids, only : mg
    use dl_mg_stencil_coefficients, only : coarse_stencils
    use dl_mg_mpi
    use dl_mg_nonlinear_model, only : has_steric_weight
    use dl_mg_alloc, only : mg_allocate
    use dl_mg_types, only : mg_get_idx
    implicit none

    integer, intent(in) :: eq_type
    real(wp), intent(in)           :: eps_half(isx:,isy:,isz:,:)
    real(wp), intent(in), optional :: steric_weight(isx:,isy:,isz:) 
    
    integer i,sx,sy,sz,ex,ey,ez
    integer s(3),e(3),coords(3),dims(3)
    logical use_nonlinear
    
    use_nonlinear = eq_type == EQ_PBE_FAS .or. eq_type == EQ_PBE_FAS_NOSTERIC .or. &
         eq_type == EQ_PBE_NEWTON .or. eq_type == EQ_PBE_NEWTON_NOSTERIC .or. &
         eq_type == EQ_LINEAR_PBE_POT

    
    call mg_allocate(mg, use_nonlinear, has_steric_weight)

    call mg_get_idx(mg(t), sx, ex, sy, ey, sz, ez)
    s(:) = [ sx, sy, sz ]
    e(:) = [ ex, ey, ez ]
    
#ifdef MPI

    call get_mpi_grid(mg(t)%comm, dims=dims, coords=coords)

    do i=1,3
    if (bc(i) ==  DL_MG_BC_DIRICHLET) then
       if (coords(i) == 0) then
          s(i) = s(i)-1;
       endif
       if (coords(i) == dims(i)-1) then
          e(i) = e(i)+1
       end if
    endif
 end do

    mg(t)%c(s(1):e(1), s(2):e(2), s(3):e(3), 1:3) = &
         eps_half(s(1):e(1), s(2):e(2), s(3):e(3), 1:3)

    do i =1, 3
         !write(0,*) 'c', i
         call exchange_halo_begin(i, exch_sr+exch_rl, 1, mg(t),"c")
         call exchange_halo_end(i, exch_sr+exch_rl, 1, mg(t),"c")
      enddo

#else
    ! serial version, just copy f,cof, and bc to top mg
      ! OpenMP to be added

      if (bc(1) ==  DL_MG_BC_DIRICHLET) then
         s(1) = s(1)-1; e(1) = e(1)+1
      endif
      if (bc(2) ==  DL_MG_BC_DIRICHLET) then
         s(2) = s(2)-1; e(2) = e(2)+1
      endif
      if (bc(3) ==  DL_MG_BC_DIRICHLET) then
         s(3) = s(3)-1; e(3) = e(3)+1
      endif

      ! for PBC we copy only the interior points. No need for ds<> de<>, hence use ds<>
      mg(t)%c(s(1):e(1), s(2):e(2), s(3):e(3), :) = &
           eps_half(s(1):e(1), s(2):e(2), s(3):e(3), :)

      ! no need for pbc halo as it is called in steric_coef, but check again. ???
      call pbc_halo(mg(t),"c")

#endif

      call coarse_stencils(mg)
      
      call transfer_steric_weight
    
  contains

    subroutine transfer_steric_weight
      use dl_mg_common_data, only : mg_levels
      use dl_mg_nonlinear_model, only :  has_steric_weight
      use dl_mg_mpi
      use dl_mg_types, only : mg_get_idx
      implicit none

      integer t, i, j, k
      
      if (has_steric_weight) then
         t = mg_levels
         mg(t)%d(sx:ex, sy:ey, sz:ez) = steric_weight(sx:ex, sy:ey, sz:ez)

         do t = mg_levels-1, 1, -1
            if ( .not. mg(t)%active) exit
            call mg_get_idx(mg(t), sx, ex, sy, ey, sz, ez)
            do k = sz, ez
               do j = sy, ey
                  do i = sx, ex
                     mg(t)%d(i,j,k) = mg(t+1)%d(2*i-1, 2*j-1, 2*k-1)
                  enddo
               enddo
            enddo
         end do
      endif

    end subroutine transfer_steric_weight
    
  end subroutine set_mg_grids_eps_and_steric

  
  subroutine free_mg
    use dl_mg_grids, only : mg
    use dl_mg_alloc, only : mg_deallocate
    implicit none
    
    call mg_deallocate(mg)

  end subroutine free_mg
  
end module dl_mg_defco
