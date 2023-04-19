module dl_mg_conjugate_gradient_method
  implicit none

contains

  ! implements the following algorithm taken from wikipedia
  ! r[0] = f - A p[0]
  ! z[0] = M^-1 r[0] # preconditioning
  ! y[0] = z[0]
  ! k = 0
  ! DO
  !   alpha = r[k] . z[k] / y[k] . A y[k]
  !   p[k+1] = p[k] + alpha y[k]
  !   r[k+1] = r[k] - alpha A y[k]
  !   test || r[+1] || and exit if small enough
  !   z[k+1] = M^-1 r[k+1]
  !   beta = r[k+1] . z[k+1] / r[k] . z[k] # one can try
  !                                        # z[k+1] -> z[k+1]-z[k]
  !   y[k+1] = z[k+1] + beta y[k]
  !   k = k + 1
  ! END DO
  !
  ! Notes:
  !   1. DL_MG uses a semi-negative de fined operator but
  !      the algorithm works the same in this case.
  !      formally one should solve -Av = -b
  !   2. The redblack preconditioner (mg_relax) must be used with
  !      omega_sor = 1
  !
  subroutine conj_grad_solver(mg, eq_type, qrho, rho, pot, res,&
        use_pot_in, remove_source_zero_mode, ierr)
!$  use omp_lib
    use dl_mg_params, only : wp
    use dl_mg_types
    use dl_mg_common_data, only : t => mg_levels, blk, &
         isx, isy, isz, fullpbc !, nx, ny, nz
    use dl_mg_kernels
    use dl_mg_convergence_params, only : tconv_cg, &
         max_iters => max_iters_cg_, use_cg => use_cg_, &
         max_iters_vcyc_
    use dl_mg_utils, only : mg_convergence_test, &
         remove_zero_mode_potential, &
         remove_zero_mode_source, set_boundary_values
    use dl_mg_alloc, only : galloc
    use dl_mg_mpi, only : exchange_full_halos, get_myid
    use dl_mg_errors
    implicit none

    type(mg_t), target, intent(inout) :: mg(:)
    integer, intent(in) :: eq_type
    real(wp), intent(in)            :: qrho, rho(isx:, isy:, isz:)
    real(wp), target, intent(inout) :: pot(isx:, isy:, isz:)
    real(wp), optional, intent(out) :: res(isx:, isy:, isz:)
    logical, optional, intent(in)   :: use_pot_in !< if absent or false
    !! the inital solution is set 0 inside domain
    integer, optional, intent(in)   :: remove_source_zero_mode
    integer, intent(out) :: ierr
    
    type(mg_t), pointer :: gp
    integer converged, iter, iblk, myid, myth
    integer sx, sy, sz, ex, ey, ez
    real(wp) fnorm, rnorm, alpha, beta, rz_prod, rz_prod_prev, xy
    real(wp), allocatable, save :: r(:,:,:), y(:,:,:), ay(:,:,:),&
         fa(:,:,:), pa(:,:,:)
    integer r_t, y_t, ay_t, fa_t, pa_t, z_t !array types for kernels
    logical use_precond
    real(wp), pointer :: z(:,:,:) ! preconditioned residual
    target y, pa
#ifdef HAVE_CONTIGUOUS
    contiguous z
#endif
    ! debug
    !real(wp) dy, sumy

    
    if (present(use_pot_in)) then
       if (use_pot_in) then
          call handle_error(DL_MG_ERR_UNSPECIFIED,&
               msg="use input potential is not implemented in CG")
       end if
    end if
    
    ierr = DL_MG_SUCCESS

    use_precond = max_iters_vcyc_ > 0
    
    gp => mg(t)
    z => gp%p

    call mg_get_idx(gp, sx, ex, sy, ey, sz, ez)
    ay_t = p_vtype    
    y_t  = p_vtype
    fa_t = f_vtype
    r_t  = f_vtype
    pa_t = p_vtype
    z_t  = p_vtype
    !$omp master
    call galloc(gp, ay, ay_t)
    call galloc(gp,  y, y_t)
    call galloc(gp, fa, fa_t)
    call galloc(gp,  r, r_t)
    call galloc(gp, pa, pa_t)
    !$omp end master
    !$omp barrier
    myid = get_myid()
    myth = 0
    !$  myth = omp_get_thread_num()

    ! set boundary to 0 for y which are needed in Ay computation
    call set_boundary_values(gp, y, y_t, val=0.0_wp)

    
    ! set the potential to zero inside the domain
    call blocked_update2(gp, pot(sx:,sy:,sz:), z_vtype, 0.0_wp, sa=0.0_wp)
    call blocked_update2(gp, pa, pa_t, 0.0_wp, sa=0.0_wp)
    ! copy boundary values to pa, need for residual calculation
    call set_boundary_values(gp, pa, pa_t, aval=pot, av_type=y_vtype)
    ! step 0
    call blocked_scaled_copy2(gp, qrho, rho, y_vtype,&
         fa, fa_t)
    ! fnorm is needed for the convergence test
    fnorm=sqrt(blocked_dotproduct(gp, &
         fa, fa_t, fa, fa_t))
    ! compute residual r0 = b - A p0
    ! Note: for full PBC p0=0 and there are no boundary values
    !       hence we can save some computation
    if ( fullpbc) then
       call blocked_scaled_copy2(gp, 1.0_wp, fa, fa_t,&
         r, r_t)
    else
       call blocked_residual(gp, eq_type, pa, pa_t, 1.0_wp,&
            fa, fa_t, r, r_t)
    endif
    call set_boundary_values(gp, r, r_t, val=0.0_wp)

    ! precondition: z0 = M^-1 r0
    ! y0 = z0 = mgl%p
    if (use_precond) then
       call precon
       rz_prod = blocked_dotproduct(gp, r, r_t, z, z_t)
       call blocked_scaled_copy2(gp, 1.0_wp, z, z_t, y, y_t)
    else
       rz_prod = blocked_dotproduct(gp, r, r_t, r, r_t)
       call blocked_scaled_copy2(gp, 1.0_wp, r, r_t, y, y_t)
    end if

    do iter = 1, max_iters
       ! compute -Ay, r is a dummy
       call blocked_residual(gp, eq_type, y, y_t, 0.0_wp, r, r_t, &
            ay, ay_t)
       
       ! `-` sign is form -A y
       alpha = - rz_prod / blocked_dotproduct(gp, y, y_t, ay, ay_t)
       
       ! compute new pot
       ! pot[k+1] = pot[k] + alpha * y[k]
       call blocked_update2(gp, pot(sx:,sy:,sz:), z_vtype, &
            alpha, y, y_t)
       
       ! r[k+1] = r[k] + alpha * ( - A y[k])
       call blocked_update2(gp, r, r_t, alpha, ay, ay_t)
       if (.not. use_precond) then
          ! K=0 is removed automatically in preconditioner
          ! from the math I think that r_{K=0} must be set to zero
          ! during iterarion without preconditioning
          !
          ! but numerical test that same solution norm as MG
          ! is reached whithout any mode removal
          ! hence the continue
          !
          ! Anyway this is just for testing, without preconditioning
          ! CG is very slow
          !
          !call remove_zero_mode_potential(gp, eq_type, r(sx:, sy:, sz:))
          continue
       endif
       
       ! test
       rnorm=sqrt(blocked_dotproduct(gp, r, r_t, r, r_t))

       call mg_convergence_test(gp, t, eq_type, tconv_cg, iter, &
            converged, fnorm, rnorm, pot)
       if (converged >= DL_MG_SUCCESS) exit

       rz_prod_prev=rz_prod
       ! compute cg beta = r[k+1] . z[k+1] / r[k] . z[k]
       ! precondition: get z in mgl%p
       if (use_precond) then
          call precon
          rz_prod = blocked_dotproduct(gp, r, r_t, z, z_t)
       else
          rz_prod = blocked_dotproduct(gp, r, r_t, r, r_t)
       end if
          
       beta = rz_prod / rz_prod_prev
       ! compute new y[k+1] = beta * y[k] + z[k+1]
       if (use_precond) then
          call blocked_update2(gp, y, y_t, 1.0_wp, z, z_t, beta)
       else
          call blocked_update2(gp, y, y_t, 1.0_wp, r, r_t, beta)
       endif

    end do

    ! set ierr
    ierr = converged

    if (present(res)) then
       call blocked_scaled_copy2(gp, 1.0_wp, r, r_t, &
            res, y_vtype, with_boundary=.true.)
    endif

    !$omp barrier
    !$OMP MASTER
    deallocate(r,ay,y,fa,pa)
    !$OMP END MASTER
    
  contains
    
    subroutine precon
      use dl_mg_multigrid_method, only : mg_solver, mg_relax
      implicit none

      real(wp), allocatable, save :: pa(:,:,:)
      real(wp), pointer :: vin(:,:,:)
      real(wp) c7
      integer sv(3),ev(3),i,j,k, iter, ierr


      call mg_solver(mg, eq_type, 1.0_wp, r(isx:,isy:,isz:), &
           z(isx:,isy:,isz:), ierr=ierr)

      return
      
! Jacobi iterator
!!$      vin => mgl%p
!!$
!!$      !$omp master
!!$      allocate( pa(mgl%psx:mgl%pex, &
!!$           mgl%psy:mgl%pey, &
!!$           mgl%psz:mgl%pez))
!!$      !$omp end master
!!$      !$omp barrier
!!$
!!$      vin = 0.0_wp
!!$      
!!$      do iter = 1, precon_niter
!!$         call exchange_full_halos(mgl,"p")
!!$         !$OMP DO SCHEDULE(STATIC,1)
!!$         do iblk = 1, blk(t)%nblks(2)
!!$            do i =1, 3
!!$               sv(i) = blk(t)%start(i, iblk)
!!$               ev(i) = blk(t)%end(i, iblk)
!!$            enddo
!!$
!!$            do k = sv(3), ev(3)
!!$               do j = sv(2), ev(2)
!!$                  do i = sv(1), ev(1)
!!$                     c7 = mgl%c(i,j,k,1) + mgl%c(i-1,j,   k,1) + &
!!$                          mgl%c(i,j,k,2) + mgl%c(i  ,j-1 ,k,2) + &
!!$                          mgl%c(i,j,k,3) + mgl%c(i  ,j   ,k-1,3)
!!$
!!$                     pa(i,j,k) =  (-mgl%f(i,j,k) &
!!$                          +mgl%c(i-1,j,k,  1) * vin(i-1,j,  k  ) &
!!$                          +mgl%c(i  ,j,k,  1) * vin(i+1,j,  k  ) &
!!$                          +mgl%c(i,  j-1,k,2) * vin(i  ,j-1,k  ) &
!!$                          +mgl%c(i,  j,  k,2) * vin(i  ,j+1,k  ) &
!!$                          +mgl%c(i,  j,k-1,3) * vin(i  ,j,  k-1) &
!!$                          +mgl%c(i,  j,k,  3) * vin(i  ,j,  k+1))/c7
!!$                  enddo
!!$               enddo
!!$            enddo
!!$         enddo
!!$         !$OMP ENDDO NOWAIT
!!$
!!$         call blocked_scaled_copy2(mgl, 1.0_wp, pa, p_vtype, vin, p_vtype)
!!$      end do
!!$      !$omp barrier
!!$      !$omp master
!!$      deallocate(pa)
!!$      !$omp end master
!!$      !$omp barrier

    end subroutine precon

    
    subroutine test_dirichlet_bc_consistency(mg, v, vtype, ierr)
      ! test if the Dirichlet boundary condition
      ! have been changed during computations
      use dl_mg_common_data, only : mg_levels, bc
      use dl_mg_mpi, only : get_mpi_grid
      use dl_mg_errors
      use dl_mg_types, only : mg_t, mg_get_idx
      implicit none

      type(mg_t), intent(in) :: mg
      real(wp), target, intent(in) :: v(:,:,:)
      integer, intent(in) :: vtype

      integer, optional, intent(out) :: ierr

      integer dims(3), t, count, i, a(3,6), bv(3)
      integer sx, ex, sy, ey, sz, ez
      character(len=255) errmsg
      real(wp), pointer :: pv(:,:,:)

      if(present(ierr)) ierr = DL_MG_SUCCESS

      ! the other equation types scale the boundary
      ! a smarter test is need for the other cases
      !if ( eq_type /= EQ_POISSON ) return

      !$omp master

      count = 0
      t = mg_levels
      call get_mpi_grid(mg%comm, dims=dims)
      call mg_get_idx(mg, sx, ex, sy, ey, sz, ez)
      call get_lbound(mg, vtype, bv)
      pv(bv(1):,bv(2):,bv(3):) => v

      if (mg%coords(1) == 0 .and. bc(1) == DL_MG_BC_DIRICHLET) then
         if (any(pv(sx-1, sy:ey, sz:ez) /= 0.0_wp)) then
            count = count+1
            a(:,count) =  mg%coords(:)
         endif
      endif

      if (mg%coords(1) == dims(1) - 1 .and. bc(1) == DL_MG_BC_DIRICHLET ) then
         if (any(pv(ex+1, sy:ey, sz:ez) /= 0.0_wp )) then
            count = count+1
            a(:,count) =  mg%coords(:)
         endif
      endif

      if (mg%coords(2) == 0 .and. bc(2) == DL_MG_BC_DIRICHLET) then
         if(any(pv(sx:ex, sy-1, sz:ez) /= 0.0_wp )) then
            count = count+1
            a(:,count) =  mg%coords(:)
         endif
      endif

      if (mg%coords(2) == dims(2) - 1 .and. bc(2) == DL_MG_BC_DIRICHLET) then
         if(any(pv(sx:ex, ey+1, sz:ez) /= 0.0_wp )) then
            count = count+1
            a(:,count) =  mg%coords(:)
         endif
      endif

      if (mg%coords(3) == 0 .and. bc(3) == DL_MG_BC_DIRICHLET) then
         if (any(pv(sx:ex, sy:ey, sz-1) /= 0.0 )) then
            count = count+1
            a(:,count) =  mg%coords(:)
         endif
      endif

      if (mg%coords(3) == dims(3) - 1 .and. bc(3) == DL_MG_BC_DIRICHLET) then
         if(any(pv(sx:ex, sy:ey, ez+1) /= 0.0 )) then
            count = count+1
            a(:,count) =  mg%coords(:)
         endif
      endif

      if ( count /= 0 ) then
        write(*,*)  "boundary modified by MPI ranks with coordinates : ", &
              ( a(1:3,i), i=1,count)
      endif
      !if (.not. check_assertion(count == 0, DL_MG_ERR_MOD_BC)) then
      !   write(errmsg,'(a,18I5)') "boundary modified by MPI ranks with coordinates : ", &
      !        ( a(1:3,i), i=1,count)
      !   call handle_error(DL_MG_ERR_MOD_BC, ierr, errmsg)
      !endif

      !$omp end master

    end subroutine test_dirichlet_bc_consistency

    
  end subroutine conj_grad_solver

  
!!$  subroutine set_params
!!$    implicit none
!!$
!!$    character(len=20) :: val
!!$    integer len, idx
!!$
!!$    use_precond=.true.
!!$    call get_environment_variable("DL_MG_USE_CG_PRECON",val, len)
!!$    if (len > 0) then
!!$       idx = index(val,',')
!!$       select case (val(1:1))
!!$       case('f','F','n','N')
!!$          use_precond = .false.
!!$       end select
!!$       if (idx > 0) then
!!$          read(val(idx+1:),*)precon_niter
!!$       end if
!!$    end if
!!$    
!!$  end subroutine set_params
!!$    
!!$
!!$  subroutine report_params
!!$    use dl_mg_mpi
!!$    use dl_mg_common_data, only : report_unit
!!$    implicit none
!!$
!!$    integer myid
!!$
!!$    myid = get_myid()
!!$
!!$    if (myid == 0) then
!!$       if (use_precond) then
!!$          write(report_unit,'(a,1x,I5,/)') 'CG with preconditioner, iterations:',precon_niter
!!$       else
!!$          write(report_unit,'(a,/)') 'CG without preconditioner'
!!$       end if
!!$    end if
!!$    
!!$  end subroutine report_params
  
end module dl_mg_conjugate_gradient_method
