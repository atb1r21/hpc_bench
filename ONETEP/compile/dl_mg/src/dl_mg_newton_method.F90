module dl_mg_newton_method
  implicit none
  
contains

  subroutine newton_solver(mg, eq_type, alpha, rho, pot, res, &
       use_pot_in, ierr)
    use dl_mg_mpi, only : get_myid
    use dl_mg_common_data, only : mg_levels, nthreads, &
         DL_MG_SUCCESS, blk, mg_comm, isx, isy, isz
    use dl_mg_kernels
    use dl_mg_errors
    use dl_mg_mpi_header
    use dl_mg_convergence_params, only : tconv_newton, niter_newton => max_iters_newton_, use_cg => use_cg_
    use dl_mg_utils, only : compute_norm, remove_zero_mode_source, &
         mg_convergence_test, set_boundary_values, set_mass_term
    use dl_mg_nonlinear_model, only : neutralisation_method, mu
    use dl_mg_mpi, only : exchange_full_halos
    use dl_mg_alloc, only : galloc
    use dl_mg_types

    implicit none

    type(mg_t), target, intent(inout) :: mg(:)
    integer, intent(in)             :: eq_type
    real(wp), intent(in)            :: alpha, rho(isx:, isy:, isz:)
    real(wp), target, intent(inout) :: pot(isx:, isy:, isz:)
    real(wp), optional, intent(out) :: res(isx:, isy:, isz:)
    logical, optional, intent(in)   :: use_pot_in !< if absent or false
    !! the inital solution is set 0 inside domain
    integer, intent(out)      :: ierr

    real(wp), allocatable, save :: fma(:,:,:),&
         linsol(:,:,:) , &
         rn(:,:,:), rl(:,:,:), upv(:,:,:), &
         dv(:,:,:)
    type(mg_t), pointer :: gp
    integer i, t, iblk, ierr_lin, myid, converged
    integer sidx(3), eidx(3), p_t, r_t, f_t, b_t
    real(wp) res_u, res_upv, normb, res_utest, fnorm
    target linsol
    
      !> switch for verbose reporting from Newton iterations
  logical, save :: newton_print=.false.

    !write(0,*) ' in newton solver'

    ierr = DL_MG_SUCCESS

    call get_environment_variable("DL_MG_NEWTON_PRINT", length=i)
    if (i > 0) newton_print = .true.

    t = mg_levels
    gp => mg(t)

    if ( newton_print) then
       myid = get_myid()
    endif

    ! attetion at dealloc, only pf0 is tested
    !$OMP MASTER
    allocate(fma(gp%psx:gp%pex,gp%psy:gp%pey,gp%psz:gp%pez), & ! f - A u
         rn(gp%rsx:gp%rex,gp%rsy:gp%rey,gp%rsz:gp%rez), & ! non-linear residual
         rl(gp%rsx:gp%rex,gp%rsy:gp%rey,gp%rsz:gp%rez), & ! linear residual
         upv(gp%rsx:gp%rex,gp%rsy:gp%rey,gp%rsz:gp%rez))  ! u + lam * v
    call galloc(gp, linsol, p_vtype)
    !$OMP END MASTER
    !$OMP BARRIER
    p_t = p_vtype
    r_t = r_vtype
    f_t = f_vtype
    b_t = y_vtype ! i.e. with boundary, input array
    
    ! start and end indices for ranks on the domain boundaries
    !call set_boundary_indices

    ! set the potential to zero inside the domain
    call blocked_update2(gp, pot, b_t, 0.0_wp, sa=0.0_wp)
    ! we know the error on the boundary :)
    call set_boundary_values(gp, linsol, p_t, val=0.0_wp) 

    ! need for the convergence test
    fnorm=abs(alpha) * sqrt(blocked_dotproduct(gp, rho, b_t, rho , b_t))

    do i=1, niter_newton

       ! called here because diagmatvec uses mass terms
       call set_mass_term(mg, pot, ierr)
       
       ! compute f - Au + N(u)  and store it in rt
       ! f - Au is stored in fma, needed in damping procedure
       call blocked_residual(gp, EQ_POISSON, pot, b_t, alpha, rho, y_vtype, fma, p_t)
       !call blocked_scaled_copy2(gp, 1.0_wp, sol, s_t, rn, r_t)
       call blocked_diagmatvec(gp, eq_type, rn, r_t, vi=pot, vi_type=b_t) ! u -> N(u)
       call blocked_update2(gp, rn, r_t, 1.0_wp, fma, p_t ) ! f -Au +N(u)

       res_u = blocked_dotproduct(gp, rn, r_t, rn, r_t) ! NB this is norm**2

       if ( newton_print) then
          !$OMP MASTER
          if ( myid == 0 ) then
             write(*,*) 'newton iter', i, fnorm, sqrt(res_u) !, (res_u-sqrt(res_utest))/res_u
          endif
          !$OMP END MASTER
       endif

       ! for full pbc one needs to satisfy
       ! - \int_V N' v = \int_V  (f -Au + N(u))
       call remove_zero_mode_source(gp, EQ_LINEAR_PBE_POT, 2, &
            rho=rn(isx:, isy:, isz:))

       ! solve J v = r
       call solve_linear(i, ierr_lin)

       if ( newton_print) then
          normb= blocked_dotproduct(gp, gp%r, r_vtype, rn , r_vtype)
          if ( myid == 0) then
             !$OMP MASTER
             write(*,*) '|| F(u) ||^2 - (lin_res, F(u))', i, res_u - normb
             !$OMP END MASTER
          end if
       end if

       if (ierr_lin /= DL_MG_SUCCESS) then
          !$OMP MASTER
          !converged = DL_MG_ERR_SMOOTHER_NEWTON_SOR
          if ( newton_print) then
             if ( myid == 0) then
                write(*,*)"WARNING: Multigrid in newton_solver fail to converge, ..."
             endif
          endif
          !$OMP END MASTER
          !test if (F(u), r) < (F(u), F(u))
          normb = blocked_dotproduct(gp, rl, r_vtype, rn , r_vtype)
          !$OMP MASTER
          if ( myid == 0 .and. normb >= res_u) then
             write(*,*)"WARNING: (F(u), r) >= (F(u), F(u))"
             write(*,*)"         the solution of the linear system is not guarantee to be a descending direction."
          endif
          !$OMP END MASTER
          !exit
       endif

       call compute_damping_parameter(ierr)
       if (ierr /= DL_MG_SUCCESS) exit

       ! does work only when fullpbc
       call compute_mu(i)


       ! no need to have u in gp%p becase the residual and source are passed
       ! in an separate arrays
       call mg_convergence_test(gp, t, eq_type, tconv_newton, &
            i, converged, fnorm, res_upv, pot)
       if (converged >=0 ) ierr = converged ! converged is set to -1 while working
       !> \todo what should happen here if converged has the value DL_MG_ERR_RES_RATIO?
       !! exit or keep going
       if (converged == DL_MG_SUCCESS) then
          exit
       endif
    enddo

    ! write the solution in p array even if convergence was not reached
!!$    !$OMP MASTER
!!$    gp%w(sidx(1):eidx(1),sidx(2):eidx(2),sidx(3):eidx(3)) = &
!!$         gp%w(sidx(1):eidx(1),sidx(2):eidx(2),sidx(3):eidx(3)) - &
!!$         sum(gp%w(sidx(1):eidx(1),sidx(2):eidx(2),sidx(3):eidx(3))&
!!$         *gp%d(sidx(1):eidx(1),sidx(2):eidx(2),sidx(3):eidx(3)))&
!!$         /sum(gp%d(sidx(1):eidx(1),sidx(2):eidx(2),sidx(3):eidx(3)))
!!$    !$OMP END MASTER
!!$    !$OMP BARRIER
!!$    call compute_mu(i)
    call shift_potential_and_mu(gp, ierr)
    ! no return when ierr /= DL_MG_SUCCESS because
    ! we copy the current approximation in gp%p
    !call blocked_scaled_copy2(gp, 1.0_wp, sol, p_t, pot, y_vtype)
    if ( newton_print) then
       normb=blocked_dotproduct(gp, pot, y_vtype, pot, y_vtype)
       if ( myid == 0) then
          !$OMP MASTER
          write(*,*) 'newton iter, p norm at exit', i, sqrt(normb)!, &
               !sum(gp%p(sidx(1):eidx(1),sidx(2):eidx(2),sidx(3):eidx(3))&
               !*gp%d(sidx(1):eidx(1),sidx(2):eidx(2),sidx(3):eidx(3)))&
               !,sum(gp%d(sidx(1):eidx(1),sidx(2):eidx(2),sidx(3):eidx(3))), mu
          !$OMP END MASTER
       end if
    end if

    !$OMP BARRIER
    !$OMP MASTER
    if (allocated(fma)) deallocate(fma, rn, rl, upv, linsol)
    !$OMP END MASTER

  contains

    subroutine solve_linear(iter_newton, ierr_lin)
      !use dl_mg_convergence_params, only : tconv_vcycle, maxcycles=> max_iters_vcyc_
      use dl_mg_multigrid_method, only : mg_solver
      use dl_mg_conjugate_gradient_method, only : conj_grad_solver
      implicit none
      integer, intent(in) :: iter_newton
      integer, intent(inout) :: ierr_lin

      integer iter, conv_mg


      if ( newton_print) then
         !$OMP MASTER
         if (myid == 0) then
            write(*,*) 'linear solve'
         endif
         !$OMP END MASTER
      endif
      ierr_lin = DL_MG_SUCCESS

      if (use_cg) then
         !call conj_grad_solver(mg,eq_type,ierr_lin)
         call conj_grad_solver(mg, EQ_LINEAR_PBE_POT, 1.0_wp,&
              rn(isx:,isy:,isz:), linsol(isx:,isy:,isz:), &
              rl(isx:,isy:,isz:), &
              use_pot_in=.false., ierr=ierr_lin)
         !continue
      else
         call mg_solver(mg, EQ_LINEAR_PBE_POT, 1.0_wp, &
              rn(isx:,isy:,isz:), linsol(isx:,isy:,isz:), &
              rl(isx:,isy:,isz:), &
              use_pot_in=.false., ierr=ierr_lin)
!!$         do iter = 1, maxcycles
!!$            call mg_vcycle(mg, EQ_LINEAR_PBE_POT)
!!$            !$OMP BARRIER
!!$            call mg_convergence_test(mg(mg_levels), mg_levels, &
!!$                 EQ_LINEAR_PBE_POT, tconv_vcycle, &
!!$                 iter, conv_mg)
!!$            if ( newton_print) then
!!$               normb =  compute_norm(gp, "r")
!!$               !$OMP MASTER
!!$               if ( myid == 0) then
!!$               write(*,*) 'linear conv test', iter, normb
!!$            end if
!!$            !$OMP END MASTER
      end if
      
    end subroutine solve_linear


    subroutine compute_damping_parameter(ierr)
      implicit none

      integer, intent(inout) :: ierr
      real(wp), parameter :: zv = 0.5_wp, zd = 0.5_wp ! line backtraking parameters
      integer count
      real(wp) lam, derfnc_v
      logical found

      !$OMP BARRIER
      !$OMP MASTER
      call galloc(gp, dv, p_t)
      !$OMP END MASTER
      !$OMP BARRIER

      ! Compute `-1 * lap v` and store it in dv
      !gp%z => linsol
      !call exchange_full_halos(gp,"z")
      !call blocked_lap(gp, 1.0_wp, linsol, p_t, dv, p_t)
      call blocked_residual(gp, EQ_POISSON, linsol, p_t,&
           0.0_wp, gp%r, r_t, dv, p_t)
      
      count = 1
      lam = 1.0_wp
      found = .false.

      do while ( lam > 1.e-5_wp .and. count <= 1000 )
         ! compute || f - Au - s Av + N(u + s v) ||
         ! use gp%r, upv for intermediate storage
         !call blocked_scaled_copy2(gp, 1.0_wp, fma, p_t, gp%r, r_t)
         !call blocked_update2(gp, gp%r, r_vtype, lam, dv, p_t)
         !call blocked_scaled_copy2(gp, lam, linsol, p_t, upv, r_t)
         !call blocked_update2(gp, upv, r_t, 1.0_wp, sol, p_t)
         !call blocked_diagmatvec(gp, eq_fas, upv, r_vtype)
         !call blocked_update2(gp, gp%r, r_t, 1.0_wp, upv, r_t)
         call blocked_scaled_copy2(gp, 1.0_wp, pot, b_t, gp%r, r_t)
         call blocked_update2(gp, gp%r, r_vtype, lam, linsol, p_t)
         call blocked_diagmatvec(gp, eq_type, gp%r, r_vtype)
         call blocked_update2(gp, gp%r, r_vtype, lam, dv, p_t) ! dv = -lap linsol
         call blocked_update2(gp, gp%r, r_vtype, 1.0_wp, fma, p_t) 
         
         res_upv = blocked_dotproduct(gp, gp%r, r_vtype, gp%r , r_vtype)
         if (newton_print) then
            !$OMP MASTER
            if (myid == 0) then
               write(*,*) 'damping: residuals u, u+v', i, res_u, res_upv, lam
            endif
            !$OMP END MASTER
         endif

         if ( count == 1) then
            if ( res_u > res_upv) then
               found = .true.
               exit
            else
               ! compute
               ! d (F(u) + s * v), F(u + s * v))/ds | s=0, i.e.
               ! -2 ( -Av + N'(u) v, f - Au + N(u) )
               ! for line back tracking
               !call blocked_scaled_copy2(gp, 1.0_wp, linsol, p_t, gp%r, r_t)
               call blocked_diagmatvec(gp, EQ_LINEAR_PBE_POT, gp%r, r_t, vi=linsol, vi_type=p_t)
               call blocked_update2(gp, gp%r, r_t, 1.0_wp, dv, p_t)
               derfnc_v = 2.0_wp * blocked_dotproduct(gp, gp%r, r_t, rn, r_t)
               if (newton_print) then
                  !$OMP MASTER
                  if (myid == 0) then
                     write(*,*) 'damping: functional derivative in v dir', derfnc_v
                  endif
                  !$OMP END MASTER
               endif
               if (derfnc_v >= 0.0_wp) then
                  ierr = DL_MG_ERR_NEWTON_DAMP_DERFUN
                  exit
               endif
               count = count + 1
            end if
         else
            ! use line backtracking to damp v
            if (res_u - res_upv  < - zd * derfnc_v * lam ) then
               lam = zv * lam
               count = count + 1
            else
               found = .true.
               exit
            endif
         endif ! count == 1
      enddo ! search damping parameter loop

      !$OMP BARRIER
      !$OMP MASTER
      deallocate(dv)
      !$OMP END MASTER
      !$OMP BARRIER

      ! store the updated u  in gp%w
      if ( found ) then
         call blocked_update2(gp, pot, b_t, lam, linsol, p_t)
      else
         if (ierr == DL_MG_SUCCESS) then
            ! that is if no other error occured
            if ( .not. check_assertion(.false., DL_MG_ERR_NEWTON_DAMP)) then
               ierr = DL_MG_ERR_NEWTON_DAMP
            endif
         end if
         !call error_abort(" compute_damping_parameter : failed to find a damping parameter for Newton method!")
      endif

    end subroutine compute_damping_parameter


    subroutine set_boundary_indices
      use dl_mg_common_data, only : npx, npy, npz, bc
      use dl_mg_mpi, only : get_mpi_grid
      implicit none

      integer i, coord(3), pe(3)

      sidx = (/ gp%sx,  gp%sy,  gp%sz  /)
      eidx = (/ gp%ex,  gp%ey,  gp%ez  /)
      pe   = (/ npx -1, npy -1, npz -1 /)
      call get_mpi_grid(gp%comm, coords=coord)

      do i=1,3
         if (coord(i) == 0  .and. bc(i) == DL_MG_BC_DIRICHLET    ) sidx(i) = sidx(i) - 1
         if (coord(i) == pe(i) .and. bc(i) == DL_MG_BC_DIRICHLET ) eidx(i) = eidx(i) + 1
      enddo

    end subroutine set_boundary_indices


    !> computes the chemical potential from the conservation of number ions equation
    !! for a given potential
    subroutine compute_mu(iter)
      use dl_mg_nonlinear_model, only : capexp, mu, q => qion, c => cion, &
           has_steric_weight, steric_weight_sum, nion_types => nion, beta
      use dl_mg_common_data, only : nx, ny, nz, dx, dy, dz, fullpbc
      implicit none

      integer, intent(in) :: iter

      integer i,j,k, ion_type, iblk, ierr
      integer, pointer :: s(:), e(:)
      real(wp) x, stw, dv
      real(wp), allocatable, save :: lsum(:), w(:)

      if (.not. fullpbc) return

      !$OMP SINGLE
      allocate(lsum(nion_types), w(nion_types))
      lsum = 0.0_wp
      !$OMP END SINGLE

      if (has_steric_weight) then
         !$OMP DO SCHEDULE(STATIC,1) REDUCTION(+:lsum)
         do iblk = 1, size(blk(t)%start, dim=2)
            s => blk(t)%start(:,iblk)
            e => blk(t)%end(:,iblk)
            do k = s(3), e(3)
               do j = s(2), e(2)
                  do i = s(1), e(1)
                     stw = gp%d(i,j,k)
                     do ion_type = 1, nion_types
                        x = -beta * q(ion_type) * pot(i,j,k)
                        lsum(ion_type) = lsum(ion_type) + stw * capexp(x)
                     end do
                  end do
               end do
            end do
         end do
         !$OMP ENDDO
      else
         !$OMP DO SCHEDULE(STATIC,1) REDUCTION(+:lsum)
         do iblk = 1, size(blk(t)%start, dim=2)
            s => blk(t)%start(:,iblk)
            e => blk(t)%end(:,iblk)
            do k = s(3), e(3)
               do j = s(2), e(2)
                  do i = s(1), e(1)
                     do ion_type = 1, nion_types
                        x = -beta * q(ion_type) * pot(i,j,k)
                        lsum(ion_type) = lsum(ion_type) + capexp(x)
                     enddo
                  end do
               end do
            end do
         end do
         !$OMP ENDDO
      end if

      ! reduce accros MPI ranks
      !$OMP MASTER
#ifdef MPI
      call MPI_ALLREDUCE(lsum, w, nion_types, MPI_DOUBLE_PRECISION, MPI_SUM, gp%comm,ierr)
#else
      w(:) = lsum(:)
#endif

      mu(:) = - log (w(:) / steric_weight_sum)

      if (newton_print .and. myid == 0) then
         dv = dx*dy*dz
         write(*,*) 'iter mu     ', iter, mu
         write(*,*) 'total charge', sum(c(:)*q(:)*w(:)*exp(mu(:)))* dv
         write(*,*) ' total ion c', sum(c(:)*w(:)*exp(mu(:))) * dv
      end if

      deallocate(lsum, w)
      !$OMP END MASTER
      !$OMP BARRIER

    end subroutine compute_mu

    !> ensures that at the end of calculation
    !! \int_V steric_weight(r) pot(r) =0
    !! computes the corresponding chemical potential
    subroutine shift_potential_and_mu(mg, ierr)
!$    use omp_lib
      use dl_mg_params, only : dl_mg_neutralise_with_ions_auto
      use dl_mg_types, only : mg_t
      use dl_mg_kernels, only : blocked_update2
      use dl_mg_nonlinear_model, only : steric_weight_sum, beta, &
           q => qion, mu, has_steric_weight
      use dl_mg_common_data, only : fullpbc
      use dl_mg_neutralisation_with_ions, only : compute_shifted_ion_concentrations
      implicit none

      type(mg_t), intent(inout) :: mg
      integer, intent(inout) :: ierr

      integer i
      real(wp) s1, pot_ave, x, cs

      integer, save :: ierr_shared

!     if (fd_order > 2) return ! pot shift is done after defect correction
!     not possible because fd_order is not passed this deep

      if (.not. fullpbc) return

      if (has_steric_weight) then
         s1 = blocked_dotproduct(mg, pot, b_t, gp%d, f_vtype)
      else
         s1 = blocked_dotproduct(mg, pot, b_t)
      end if


      pot_ave = s1/steric_weight_sum
      call blocked_update2(mg, pot, b_t, -pot_ave)

      ! shift mu_i -> mu_i - beta * q_i s1
      !$OMP MASTER
      mu(:) = mu(:) - beta * q(:) * pot_ave


      ! get a new concentration after the mu was computed
      if (neutralisation_method == dl_mg_neutralise_with_ions_auto &
           .or. &
          neutralisation_method == dl_mg_neutralise_with_ions_auto_linear) then
         call compute_shifted_ion_concentrations(ierr)
         ierr_shared = ierr
      end if
      !$OMP END MASTER
      !$OMP BARRIER
!$      if (omp_get_thread_num() > 0) then
!$        ierr = ierr_shared
!$      endif
    end subroutine shift_potential_and_mu

  end subroutine newton_solver

end module dl_mg_newton_method
