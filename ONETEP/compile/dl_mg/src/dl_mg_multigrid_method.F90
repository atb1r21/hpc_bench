!> \brief Top routines for the multigrid algorith: restrict, relax, prolungate and some ancillaries
!!
!! Lucian Anton
!!
module dl_mg_multigrid_method
  !> contains the subroutines used in V cycle steps ( restrict, relax, prolongation ) and auxiliaries
  use dl_mg_params, only : wp
  use dl_mg_types
  implicit none

  !> buffers use to exchage domain halos at each mg level
  type halo_buffers_t
     real(wp), allocatable :: sbuff_ns(:,:), rbuff_ns(:,:), sbuff_we(:,:), &
          rbuff_we(:,:), sbuff_bt(:,:), rbuff_bt(:,:)
     integer np_recv(6, 0:1), np_send(6, 0:1)
  end type halo_buffers_t

  type(halo_buffers_t), allocatable, target, private :: halo_buff(:)

  logical, save :: test_dbc_corruption
  
contains

  ! this function is first in the module to help a compiler (bizarre!)
  
  real(wp) function residual_norm(mg, eq_type, task_in)
    use dl_mg_types
    use dl_mg_common_data, only : blk
    use dl_mg_kernels
    use dl_mg_utils, only : compute_norm
    implicit none
    type(mg_t), intent(inout) :: mg
    integer,  intent(in) :: eq_type
    integer, optional, intent(in) :: task_in

    integer iblk, task, level, eqtag

    level = mg%level

    if (present(task_in)) then
       task = task_in
    else
       task = RESIDUAL
    endif

    !$OMP DO SCHEDULE(STATIC,1)
    do iblk = 1, size(blk(level)%start, dim=2)
       call residual_kernel(task, blk(level), iblk, mg, eq_type)
    enddo
    !$OMP ENDDO

    residual_norm = compute_norm(mg, "r", scale=.true.)

  end function residual_norm


  subroutine mg_solver(mg, eq_type, alpha, rho, pot, res, &
       use_pot_in, remove_source_zero_mode, ierr)
    use dl_mg_params, only : DL_MG_SUCCESS
    use dl_mg_common_data, only : mg_levels, isx, iex, isy, iey, isz, iez, &
         fullpbc
    use dl_mg_convergence_params, only : tconv_vcycle, &
         maxcycles => max_iters_vcyc_
    use dl_mg_utils, only : mg_convergence_test
    use dl_mg_mpi, only : exchange_full_halos
    use dl_mg_errors
    implicit none

    type(mg_t), intent(inout)       :: mg(:)
    integer, intent(in)             :: eq_type
    real(wp), intent(in)            :: alpha, rho(isx:, isy:, isz:)
    real(wp), intent(inout)         :: pot(isx:, isy:, isz:)
    real(wp), optional, intent(out) :: res(isx:, isy:, isz:)
    logical, optional, intent(in)   :: use_pot_in !< if absent or false
    !! the inital solution is set 0 inside domain
    integer, optional, intent(in)   :: remove_source_zero_mode !< defect correction removes k=0 mode, avoids repetition
    integer, intent(out)            :: ierr

    integer iter, converged

    !$omp master
    call allocate_halo_buffers(mg)
    !$omp end master
    !$omp barrier
    call transfer_data_to_mg_solver
 

    do iter = 1, maxcycles
       call mg_vcycle(mg, eq_type)
       !$OMP BARRIER
       call mg_residual(mg(mg_levels), eq_type, halo_exch=.false.)
       call mg_convergence_test(mg(mg_levels), mg_levels, eq_type, tconv_vcycle, &
            iter, converged)
       !write(0,*) 'did conv test', iter, eq_type, converged

       if (converged == DL_MG_SUCCESS .or. converged == DL_MG_ERR_RES_RATIO) exit
    enddo

    call transfer_data_to_output
    !$omp master
    call free_halo_buffers(mg)
    !$omp end master
    !$omp barrier
    
    ierr = converged

  contains

    subroutine transfer_data_to_mg_solver
      use dl_mg_mpi
      use dl_mg_grids, only : transfer_map_recv,  transfer_map_send
      use dl_mg_kernels
      use dl_mg_utils, only : remove_zero_mode_source
      use dl_mg_common_data, only : bc
      use dl_mg_alloc, only : p_vtype, f_vtype, z_vtype
      implicit none

      integer i, j, k, t, iblk
      integer  dims(3), coords(3), left
      integer dsx, dex, dsy, dey, dsz, dez
      integer sx, ex, sy, ey, sz, ez
      integer s1,e1,s2,e2,s3,e3
      logical use_pot

      call mg_get_idx(mg(mg_levels), sx, ex, sy, ey, sz, ez)

      ! use pot input values ?
      use_pot = .false.
      if ( present(use_pot_in) ) then
         use_pot = use_pot_in
      endif
 
      t = mg_levels
      
#ifdef MPI


      if (.not. use_pot) then
         !mg(t)%p(sx:ex, sy:ey, sz:ez) = 0.0_wp
         call blocked_update2(mg(t), mg(t)%p, p_vtype, 0.0_wp, sa=0.0_wp)
      else
         !mg(t)%p(sx:ex, sy:ey, sz:ez) = pot(sx:ex, sy:ey, sz:ez)
         call blocked_scaled_copy2(mg(t), 1.0_wp, pot(sx:, sy:, sz:), &
           z_vtype, mg(t)%p, p_vtype)
      endif
      
      !mg(t)%f(sx:ex, sy:ey, sz:ez) = &
      !     alpha * rho(sx:ex, sy:ey, sz:ez)
      call blocked_scaled_copy2(mg(t), alpha, rho(sx:, sy:, sz:), &
           z_vtype, mg(t)%f, f_vtype)

      call get_mpi_grid(mg(t)%comm, dims=dims)

      coords = mg(t)%coords
      ! shifts for ranks on boundary
      dsx = 0; dex = 0
      dsy = 0; dey = 0
      dsz = 0; dez = 0

      ! transfer boundary conditions for Dirichlet BC
      
      if (coords(1) == 0 .and. bc(1) == DL_MG_BC_DIRICHLET ) then
         !write(*,*) "N"
         !$omp single
         dsx = -1
         mg(t)%p(sx+dsx, sy:ey, sz:ez) = &
              pot (isx,  sy:ey, sz:ez)
         !$omp end single nowait
      endif

      if (coords(1) == dims(1) -1 .and. bc(1) == DL_MG_BC_DIRICHLET ) then
         !write(*,*) "S"
         !$omp single
         dex = 1
         mg(t)%p(ex+dex, sy:ey, sz:ez) = &
              pot (iex,  sy:ey, sz:ez)
         !$omp end single nowait
      end if
      
      if (coords(2) == 0 .and. bc(2) == DL_MG_BC_DIRICHLET ) then
         !write(*,*) "W"
         !$omp single
         dsy = -1
         mg(t)%p(sx:ex, sy+dsy, sz:ez) = &
              pot (sx:ex,  isy, sz:ez)
         !$omp end single nowait
      endif

      if (coords(2) == dims(2) -1 .and. bc(2) == DL_MG_BC_DIRICHLET) then
         !write(*,*) "E"
         !$omp single
         dey = 1
         mg(t)%p(sx:ex, ey+dey, sz:ez) = &
              pot (sx:ex,  iey, sz:ez)
         !$omp end single nowait
      end if

      if (coords(3) == 0 .and. bc(3) == DL_MG_BC_DIRICHLET) then
         !write(*,*) "B"
         !$omp single
         dsz = -1
         mg(t)%p(sx:ex, sy:ey, sz+dsz) = &
              pot(sx:ex, sy:ey, isz)
         !$omp end single nowait
      endif

      if ( coords(3) == dims(3) -1 .and. bc(3) == DL_MG_BC_DIRICHLET) then
         !write(*,*) "T"
         !$omp single
         dez = 1
         mg(t)%p(sx:ex, sy:ey, ez+dez) = &
              pot(sx:ex, sy:ey, iez)
         !$omp end single nowait
      endif
      !$omp barrier

      ! the above function must have a barrier
      !!$omp barrier
      !$omp master
      ! not clear to me why this is needed, to investigate it
      do i =1, 3
         !write(0,*) 'f', i, t
         call exchange_halo_begin(i, exch_full, 1, mg(t),"f")
         call exchange_halo_end(i, exch_full, 1, mg(t),"f")
      enddo
      !$omp end master
      !$omp barrier

#else
      ! serial version, just copy f,cof, and bc to top mg

      !mg(mg_levels)%f(sx:ex, sy:ey, sz:ez) = alpha * rho(sx:ex, sy:ey, sz:ez)
      call blocked_scaled_copy2(mg(t), alpha, rho(sx:, sy:, sz:), &
           z_vtype, mg(t)%f, f_vtype)

      
      if (use_pot) then
         !mg(t)%p(sx:ex, sy:ey, sz:ez) = &
         !     pot(sx:ex, sy:ey, sz:ez)
         call blocked_scaled_copy2(mg(t), 1.0_wp, pot(sx:, sy:, sz:), &
           z_vtype, mg(t)%p, p_vtype)
      else
         !mg(t)%p(sx:ex, sy:ey, sz:ez) = 0.0_wp
         call blocked_update2(mg(t), mg(t)%p, p_vtype, 0.0_wp, sa=0.0_wp)
      endif

      s1 = sx; e1 = ex
      s2 = sy; e2 = ey
      s3 = sz; e3 = ez

      if (bc(1) ==  DL_MG_BC_DIRICHLET) then
         s1 = s1-1; e1 = e1+1
      endif
      if (bc(2) ==  DL_MG_BC_DIRICHLET) then
         s2 = s2-1; e2 = e2+1
      endif
      if (bc(3) ==  DL_MG_BC_DIRICHLET) then
         s3 = s3-1; e3 = e3+1
      endif

      if (bc(1) == DL_MG_BC_DIRICHLET) then
         !$omp single
         mg(t)%p(sx-1, s2:e2, s3:e3) = &
              pot(sx-1, s2:e2, s3:e3)
         !$omp end single nowait
         !$omp single
         mg(t)%p(ex+1, s2:e2, s3:e3) = &
              pot(ex+1, s2:e2, s3:e3)
         !$omp end single nowait
      endif
      if (bc(2) == DL_MG_BC_DIRICHLET) then
         !$omp single
         mg(t)%p(s1:e1, sy-1, s3:e3) = &
              pot(s1:e1, sy-1, s3:e3)
         !$omp end single nowait
         !$omp single
         mg(t)%p(s1:e1, ey+1, s3:e3) = &
              pot(s1:e1, ey+1, s3:e3)
         !$omp end single nowait
      endif
      if (bc(3) == DL_MG_BC_DIRICHLET) then
         !$omp single
         mg(t)%p(s1:e1, s2:e2, sz-1) = &
              pot(s1:e1, s2:e2, sz-1)
         !$omp end single nowait
         !$omp single
         mg(t)%p(s1:e1, s2:e2, ez+1) = &
              pot(s1:e1, s2:e2, ez+1)
         !$omp end single nowait
      endif
      !$omp barrier
#endif
      if (present(remove_source_zero_mode)) then
         if (remove_source_zero_mode == 2) then
            call remove_zero_mode_source(mg(mg_levels), eq_type, remove_source_zero_mode )
         else
            ! It should not get here
            call handle_error(DL_MG_ERR_UNSPECIFIED,&
                 msg="remove source k=0 mode key /= 2")
         end if
      end if
      
    end subroutine transfer_data_to_mg_solver

    
    subroutine transfer_data_to_output
      use dl_mg_mpi_header
      use dl_mg_grids, only : mg
      use dl_mg_common_data, only : t=> mg_levels
      use dl_mg_alloc, only : p_vtype, r_vtype, z_vtype
      use dl_mg_kernels
      implicit none

      integer sx, ex, sy, ey, sz, ez

      call mg_get_idx(mg(mg_levels), sx, ex, sy, ey, sz, ez)

      !pot(sx:ex, sy:ey, sz:ez) = mg(mg_levels)%p(sx:ex, sy:ey, sz:ez)
      call blocked_scaled_copy2(mg(t), 1.0_wp, mg(t)%p, p_vtype, &
           pot(sx:, sy:, sz:), z_vtype)
      if (present(res)) then
         !res(sx:ex, sy:ey, sz:ez) = mg(mg_levels)%r(sx:ex, sy:ey, sz:ez)
         call blocked_scaled_copy2(mg(t), 1.0_wp, mg(t)%r, r_vtype, &
           res(sx:, sy:, sz:), z_vtype)
         ! set residual 0 on boundaries
         call set_residual_boundary_values
      endif

    end subroutine transfer_data_to_output


    subroutine set_residual_boundary_values
      use dl_mg_mpi, only : get_mpi_grid
      use dl_mg_grids, only : mg
      use dl_mg_params, only : DL_MG_BC_DIRICHLET
      use dl_mg_common_data, only : bc
      implicit none
      ! set residual to 0 on boundary layers.
      ! it may be generalised later

      integer dims(3), coords(3), sx, ex, sy, ey, sz, ez

      ! needed only for Dirichlet BC

      call get_mpi_grid(mg(mg_levels)%comm, dims=dims, coords=coords)
      call mg_get_idx(mg(mg_levels), sx, ex, sy, ey, sz, ez)

      if (coords(1) == 0 .and. bc(1) == DL_MG_BC_DIRICHLET) then
         !$omp single
         res(sx-1, :, : ) = 0.0_wp
         !$omp end single nowait
      endif
      if (coords(1) == dims(1) -1 .and. bc(1) == DL_MG_BC_DIRICHLET) then
         !$omp single
         res(ex+1, :, : ) = 0.0_wp
         !$omp end single nowait
      endif

      if (coords(2) == 0 .and. bc(2) == DL_MG_BC_DIRICHLET) then
         !$omp single
         res(:, sy-1, : ) = 0.0_wp
         !$omp end single nowait
      endif
      if (coords(2) == dims(2) -1 .and. bc(2) == DL_MG_BC_DIRICHLET) then
         !$omp single
         res(:, ey+1, : ) = 0.0_wp
         !$omp end single nowait
      endif

      if (coords(3) == 0 .and. bc(3) == DL_MG_BC_DIRICHLET) then
         !$omp single
         res(:, :, sz-1 ) = 0.0_wp
         !$omp end single nowait
      endif
      if (coords(3) == dims(3) -1 .and. bc(3) == DL_MG_BC_DIRICHLET) then
         !$omp single
         res(:, :, ez+1 ) = 0.0_wp
         !$omp end single nowait
      endif

      !$omp barrier

    end subroutine set_residual_boundary_values


    subroutine test_dirichlet_bc_consistency(ierr)
      ! test if the Dirichlet boundary condition
      ! have been changed during computations
      use dl_mg_common_data, only : mg_levels, bc
      use dl_mg_mpi, only : get_mpi_grid
      use dl_mg_errors
      implicit none

      integer, optional, intent(out) :: ierr

      integer dims(3), t, count, i, a(3,6)
      integer sx, ex, sy, ey, sz, ez
      character(len=255) errmsg

      if(present(ierr)) ierr = DL_MG_SUCCESS

      ! the other equation types scale the boundary
      ! a smarter test is need for the other cases
      if ( eq_type /= EQ_POISSON ) return

      count = 0
      t = mg_levels
      call get_mpi_grid(mg(t)%comm, dims=dims)
      call mg_get_idx(mg(t), sx, ex, sy, ey, sz, ez)

      if (mg(t)%coords(1) == 0 .and. bc(1) == DL_MG_BC_DIRICHLET) then
         if (any(pot(sx-1, sy:ey, sz:ez) /= mg(t)%p(sx-1, sy:ey, sz:ez))) then
            count = count+1
            a(:,count) =  mg(t)%coords(:)
         endif
      endif

      if (mg(t)%coords(1) == dims(1) - 1 .and. bc(1) == DL_MG_BC_DIRICHLET ) then
         if (any(pot(ex+1, sy:ey, sz:ez) /= mg(t)%p(ex+1, sy:ey, sz:ez) )) then
            count = count+1
            a(:,count) =  mg(t)%coords(:)
         endif
      endif

      if (mg(t)%coords(2) == 0 .and. bc(2) == DL_MG_BC_DIRICHLET) then
         if(any(pot(sx:ex, sy-1, sz:ez) /= mg(t)%p(sx:ex, sy-1, sz:ez) )) then
            count = count+1
            a(:,count) =  mg(t)%coords(:)
         endif
      endif

      if (mg(t)%coords(2) == dims(2) - 1 .and. bc(2) == DL_MG_BC_DIRICHLET) then
         if(any(pot(sx:ex, ey+1, sz:ez) /= mg(t)%p(sx:ex, ey+1, sz:ez) )) then
            count = count+1
            a(:,count) =  mg(t)%coords(:)
         endif
      endif

      if (mg(t)%coords(3) == 0 .and. bc(3) == DL_MG_BC_DIRICHLET) then
         if (any(pot(sx:ex, sy:ey, sz-1) /= mg(t)%p(sx:ex, sy:ey, sz-1) )) then
            count = count+1
            a(:,count) =  mg(t)%coords(:)
         endif
      endif

      if (mg(t)%coords(3) == dims(3) - 1 .and. bc(3) == DL_MG_BC_DIRICHLET) then
         if(any(pot(sx:ex, sy:ey, ez+1) /= mg(t)%p(sx:ex, sy:ey, ez+1) )) then
            count = count+1
            a(:,count) =  mg(t)%coords(:)
         endif
      endif

      if (.not. check_assertion(count == 0, DL_MG_ERR_MOD_BC)) then
         write(errmsg,'(a,18I5)') "boundary modified by MPI ranks with coordinates : ", &
              ( a(1:3,i), i=1,count)
         call handle_error(DL_MG_ERR_MOD_BC, ierr, errmsg)
         return
      endif

    end subroutine test_dirichlet_bc_consistency

  end subroutine mg_solver


  !> run the V-cycle
  !!  \todo needs to handle better the convergence flag
  !!   in this version the coarse level does not check for convergence
  subroutine mg_vcycle(mg,  eq_type, niter_level1, test_conv_level1)
    !$ use omp_lib
    use dl_mg_params, only : DL_MG_SUCCESS
    use dl_mg_common_data, only : nthreads, mg_levels
    use dl_mg_convergence_params, only: &
         niter_level1_def => max_iters_level_1_, tconv_level_1, v_iteration => v_iterations_
    use dl_mg_utils, only : remove_zero_mode_potential, mg_convergence_test
    implicit none
    type(mg_t), intent(inout) :: mg(:)
    integer, intent(in) :: eq_type
    integer, optional, intent(in) :: niter_level1
    logical, optional, intent(in) :: test_conv_level1

    integer i, level, niter_newton, niter_smoother, converged, niter_level1_
    logical test_conv_level1_
    
    if (present(niter_level1)) then
       niter_level1_ = niter_level1
    else
       niter_level1_ = niter_level1_def
    end if

    test_conv_level1_ = .true.
    if (present(test_conv_level1)) then
       test_conv_level1_ = test_conv_level1
    end if
    
    do level = mg_levels, 2, -1
       call mg_relax(mg, level, v_iteration(1), eq_type)
       !$OMP BARRIER

       call mg_restrict(mg, level, eq_type)
       !$OMP BARRIER
    enddo
    ! solve at the coarsest level
    ! for level 1 iterations the book recommends v(1)+v(2)
    ! however for large grids we need quite a bit more
    if ( mg(1)%active ) then
       !write(202,*) 'norms', maxval(abs(mg(1)%f))
       do i = 1, niter_level1_
          call mg_relax(mg, 1, 3, eq_type)
          !$OMP BARRIER
          ! if monogrid (i.e. mg_levels == 1) convergence test is done in v-cycle loop
          if ( test_conv_level1_ .and. mg_levels > 1) then
             ! halo exchange is done always at the end of relax
             call mg_residual(mg(1), eq_type, halo_exch=.false.)
             call mg_convergence_test(mg(1), 1, eq_type, tconv_level_1, &
                  i,  converged)
             !$OMP BARRIER
             !write(0,*) 'lv 1 ', i, converged
             if (converged == DL_MG_SUCCESS) exit
          endif
       enddo
    endif

    do level = 1, mg_levels - 1, 1
       call mg_prolongate(mg,level, eq_type)
       !$OMP BARRIER
       call mg_relax(mg, level + 1, v_iteration(2), eq_type)
       !$OMP BARRIER
    enddo
    call remove_zero_mode_potential(mg(mg_levels), eq_type)

  end subroutine mg_vcycle


  subroutine mg_relax(mg, level, iters, eq_type)
    !$    use omp_lib
    use dl_mg_mpi, only : get_mpi_grid, minmax_grid
    use dl_mg_errors
    use dl_mg_common_data, only : nthreads, blk, mg_levels
    use dl_mg_kernels, only : redblack
    use dl_mg_convergence_params, only : omega => omega_sor_
    use dl_mg_mpi, only : exchange_full_halos
    use dl_mg_timer
    implicit none
    type(mg_t), target, intent(inout) :: mg(:)
    integer,            intent(in)    :: level
    integer,            intent(in)    :: iters
    integer,            intent(in)    :: eq_type !< used to switch the smoother
    !< in red black subroutine

    integer i, j, k, ic, it, rb, iblk, &
         t, sx, ex, sy, ey, sz, ez, coords(3), dims(3), comm

#ifdef MPI
    integer, save ::  request(12)
    integer, pointer :: np_send(:, :), np_recv(:, :)

    real(wp), pointer :: sbuff_ns(:,:), rbuff_ns(:,:), &
         sbuff_we(:,:), rbuff_we(:,:), &
         sbuff_bt(:,:), rbuff_bt(:,:)
#endif
    real(wp), pointer :: p_ptr(:,:,:)

    if ( .not. mg(level)%active ) return

    call mg_timer(start, trelax, t, tcompute)

    !set boundary conditions
    !call set_boundary_condition(mg,level)

    call mg_get_idx(mg(level), sx, ex, sy, ey, sz, ez)
    call get_mpi_grid(mg(level)%comm, dims=dims, coords=coords) ! Attention, it calls mpi under
    comm = mg(level)%comm
    t = level

#ifdef MPI

    sbuff_ns => halo_buff(t)%sbuff_ns
    rbuff_ns => halo_buff(t)%rbuff_ns
    sbuff_we => halo_buff(t)%sbuff_we
    rbuff_we => halo_buff(t)%rbuff_we
    sbuff_bt => halo_buff(t)%sbuff_bt
    rbuff_bt => halo_buff(t)%rbuff_bt
    np_recv => halo_buff(t)%np_recv
    np_send => halo_buff(t)%np_send

#endif
    p_ptr => mg(t)%p

    ! useful for quick debug
#ifdef DLMG_RELAX_PRINT
    call mg_relax_debug_print("top")
#endif
    do it = 1, iters

       do ic = 0, 1

          ! exchange halos
          call mg_timer(start, trelax, t, tcomm)

          if (.false.) then
             call exchange_full_halos(mg(t),"p")
          else
             call post_recv( halo_buff(t)%np_recv(:, 1-ic))
             call  buffs_halos_transfers(1, ic)
             !$OMP BARRIER
             call post_send(halo_buff(t)%np_send(:, 1-ic))
             !$OMP BARRIER
             call buffs_halos_transfers(-1, ic)
             !$OMP BARRIER
          end if

          call mg_timer(stop, trelax, t, tcomm)

          !$OMP DO SCHEDULE(STATIC,1)
          do iblk = 1, size(blk(t)%start, dim=2)
             call redblack(mg(t), blk(t), iblk, omega, &
                  mg(t)%f, mg(t)%c, mg(t)%w, p_ptr, ic, ex, ey, ez, eq_type)
          enddo
          !$OMP ENDDO

          !  set boundary conditions
          !call set_boundary_condition(mg,level)
!!$            if ( t == 2) then
!!$               if ( ipass == 1) then
!!$                  call dump_data(mg(t), "p ", "p1_")
!!$               else
!!$                  call dump_data(mg(t), "p ", "p2_")
!!$               endif
!!$            endif
       enddo
    enddo

    ! This ensures very close residuals for different number of MPI ranks
    ! Needed for correct parallel prolongation ( check!!!)
    ! 1 D topolgy needs only color 1 exchage, which should be faster
    if ( .false. .and. ((dims(1) == 1 .and. dims(2) == 1) .or. &
         (dims(1) == 1 .and. dims(3) == 1) .or. &
         (dims(2) == 1 .and. dims(3) == 1)) ) then
       ! only color 1 needs updating, pass 0 because buffs_halo_transfers swaps it

       call mg_timer(start, trelax, t, tcomm)

       call post_recv( halo_buff(t)%np_recv(:, 1 ))
       call  buffs_halos_transfers(1, 0)
       !$OMP BARRIER
       call post_send(halo_buff(t)%np_send(:, 1))
       !$OMP BARRIER
       call buffs_halos_transfers(-1, 0)
       !$OMP BARRIER

       call mg_timer(stop, trelax, t, tcomm)
    else
       call exchange_full_halos(mg(t),"p")
    endif

    ! useful for quick debug
#ifdef DLMG_RELAX_PRINT
    call mg_relax_debug_print("bottom")
#endif

    call mg_timer(stop, trelax, t, tcompute)

  contains

    subroutine post_recv(np_recv)
      use dl_mg_mpi_header
      implicit none
      integer, intent(in) :: np_recv(6)

      integer npoints, ierr

#ifdef MPI
      !$OMP MASTER
      request(:) = MPI_REQUEST_NULL

      !write(0,*) "post_recv ", t, np_recv

      ! N-S
      ! receive ghost points for left face
      ! Post receives the next color on the left face
      call MPI_Irecv(rbuff_ns(:,1), np_recv(1), MPI_DOUBLE_PRECISION, mg(t)%neighbor(north), &
           711, comm, request(1), ierr)

      ! receive ghost points for left face
      !  receiving from right
      call MPI_Irecv(rbuff_ns(:,2), np_recv(2), MPI_DOUBLE_PRECISION, mg(t)%neighbor(south), &
           712, comm, request(2), ierr)
      ! W-E recv

      ! Post receives for the next color on the left face
      !           write(0,*) 'recv on left ' , myid,icol, npoints
      call MPI_Irecv(rbuff_we(:,1),np_recv(3), MPI_DOUBLE_PRECISION, mg(t)%neighbor(west),&
           721, comm, request(3), ierr)

      ! receive ghost points for left face
      !  receiving from right
      !            write(0,*) 'recv on right ' , myid,icol, npoints
      call MPI_Irecv(rbuff_we(:,2), np_recv(4), MPI_DOUBLE_PRECISION, mg(t)%neighbor(east),&
           722, comm, request(4), ierr)

      ! B-T recv

      ! Post receives for the next color on the left face
      call MPI_Irecv(rbuff_bt(:,1), np_recv(5), MPI_DOUBLE_PRECISION, mg(t)%neighbor(bottom) ,&
           731, comm, request(5), ierr)

      ! receive ghost points for left face
      !  receiving from right
      call MPI_Irecv(rbuff_bt(:,2), np_recv(6), MPI_DOUBLE_PRECISION, mg(t)%neighbor(top),&
           732, comm, request(6), ierr)

      !write(0,*) "np_recv ", coords, np_recv

      !$OMP END MASTER
#endif

    end subroutine post_recv


    subroutine post_send(np_send)
      use dl_mg_mpi_header
      implicit none

      integer, intent(in) :: np_send(6)
#ifdef MPI
      integer status_sedrecv(MPI_STATUS_SIZE, 12), ierr, i, bufflen
      character(len=2047) err_buffer
      ! send boundary data

      ! N-S send

      !$OMP MASTER
      !write(0,*)'post_send ', t, np_send
      ! send to the right (i.e. s,e,t)
      call MPI_ISend(sbuff_ns(:,2), np_send(2), MPI_DOUBLE_PRECISION, mg(t)%neighbor(south), &
           711, comm, request(7), ierr)

      ! send to the left
      call MPI_ISend(sbuff_ns(:,1), np_send(1), MPI_DOUBLE_PRECISION, mg(t)%neighbor(north), &
           712, comm, request(8), ierr)

      ! W-E send

      ! send to the right (s,e,t)
      !           write(0,*) 'send to right ', myid, icol, npoints
      call MPI_ISend(sbuff_we(:,2), np_send(4), MPI_DOUBLE_PRECISION, mg(t)%neighbor(east), &
           721, comm, request(9), ierr)

      ! send to the left
      !           write(0,*) 'send to left ', myid, icol, npoints
      call MPI_ISend(sbuff_we(:,1), np_send(3), MPI_DOUBLE_PRECISION, mg(t)%neighbor(west), &
           722, comm, request(10), ierr)

      ! B-T send

      ! send to the right (s,e,t)
      call MPI_ISend(sbuff_bt(:,2), np_send(6), MPI_DOUBLE_PRECISION,  mg(t)%neighbor(top), &
           731, comm, request(11), ierr)

      ! send to the left
      call MPI_ISend(sbuff_bt(:,1), np_send(5), MPI_DOUBLE_PRECISION, mg(t)%neighbor(bottom), &
           732, comm, request(12), ierr)

      !write(0,*) 'np_send', coords, np_send, request, MPI_REQUEST_NULL

      call MPI_Waitall(12, request, status_sedrecv, ierr)
      if ( ierr /= MPI_SUCCESS) then
         write(0,*) 'error smoother waitall', coords, t, ierr
         do i =1,12
            call MPI_Error_string(status_sedrecv(MPI_ERROR,i),err_buffer,bufflen,ierr);
            write(0,*) i, err_buffer(1:min(bufflen,len(err_buffer)))
         enddo
      endif

      !$OMP END MASTER
#endif
    end subroutine post_send


    subroutine buffs_halos_transfers(dir, colour)
      use dl_mg_common_data, only : bc
      implicit none
      integer, intent(in) :: dir !  1 grid -> buff
      ! -1 buff -> halo
      integer, intent(in) :: colour

#ifdef MPI
      integer cshift, i, j, k, oc

      ! fill the transfer buffers, for black the update is done in update face
      ! send and receive is for the opposite colour
      oc = 1 - colour
      ! N-S

      !$OMP SINGLE
      if ( coords(1) > 0 .or. bc(1) == DL_MG_BC_PERIODIC) then
         if ( dir > 0 ) then
            i = 0
            do k = sz, ez
               cshift = mod(sx + sy + k + oc, 2)
               do j = sy + cshift, ey, 2
                  i = i+1
                  sbuff_ns(i,1) = p_ptr(sx, j, k)
               enddo
            enddo
         else
            i=0
            do k = sz, ez
               cshift = mod(sx - 1 + sy + k + oc, 2)
               do j = sy + cshift, ey, 2
                  i = i + 1
                  p_ptr(sx-1,j,k) = rbuff_ns(i,1)
               enddo
            enddo
         endif
      endif
      !$OMP END SINGLE NOWAIT

      !$OMP SINGLE
      if( coords(1) < dims(1)-1 .or. bc(1) == DL_MG_BC_PERIODIC) then
         if ( dir > 0 ) then
            i = 0
            do k = sz, ez
               cshift = mod(ex + sy + k + oc, 2)
               do j = sy + cshift, ey, 2
                  i = i+1
                  sbuff_ns(i, 2) = p_ptr(ex, j, k)
               enddo
            enddo
         else
            i=0
            do k = sz, ez
               cshift = mod(ex + 1 + sy + k + oc, 2)
               do j = sy + cshift, ey, 2
                  i = i + 1
                  p_ptr(ex+1, j, k) = rbuff_ns(i, 2)
               enddo
            enddo
         endif
      endif
      !$OMP END SINGLE NOWAIT

      ! W-E

      !$OMP SINGLE
      if( coords(2) > 0 .or. bc(2) == DL_MG_BC_PERIODIC) then
         if ( dir > 0 ) then
            j = 0
            do k = sz, ez
               cshift = mod(sx + sy + k + oc, 2)
               do i = sx + cshift, ex, 2
                  j = j + 1
                  sbuff_we(j,1) = p_ptr(i, sy, k)
               enddo
            enddo
         else
            j = 0
            do k = sz, ez
               cshift = mod(sx + sy - 1  + k + oc, 2)
               do i = sx + cshift, ex, 2
                  j = j + 1
                  p_ptr(i,sy-1,k) = rbuff_we(j,1)
               enddo
            enddo
         endif
      endif
      !$OMP  END SINGLE NOWAIT

      !$OMP SINGLE
      if( coords(2) < dims(2)-1 .or. bc(2) == DL_MG_BC_PERIODIC) then
         if ( dir > 0 ) then
            j = 0
            do k = sz, ez
               cshift = mod(sx + ey + k + oc, 2)
               do i = sx + cshift, ex, 2
                  j = j + 1
                  sbuff_we(j, 2) = p_ptr(i, ey, k)
               enddo
            enddo
         else
            j = 0
            do k = sz, ez
               cshift = mod(sx + ey + 1  + k + oc, 2)
               do i = sx + cshift, ex, 2
                  j = j + 1
                  p_ptr(i, ey + 1, k) = rbuff_we(j,2)
               enddo
            enddo
         endif
      endif
      !$OMP END SINGLE NOWAIT

      ! B-T send

      !$OMP SINGLE
      if ( coords(3) > 0 .or. bc(3) == DL_MG_BC_PERIODIC) then
         if ( dir > 0 ) then
            k = 0
            do j = sy, ey
               cshift = mod(sx + j + sz + oc, 2)
               do i = sx + cshift, ex, 2
                  k = k + 1
                  sbuff_bt(k,1) = p_ptr(i, j, sz)
               enddo
            enddo
         else
            k = 0
            do j = sy, ey
               cshift = mod(sx + j  + sz - 1 + oc, 2)
               do i = sx + cshift, ex, 2
                  k = k + 1
                  p_ptr(i,j,sz-1) = rbuff_bt(k,1)
               enddo
            enddo
         endif
      endif
      !$OMP END SINGLE NOWAIT

      !$OMP SINGLE
      if( coords(3) < dims(3)-1 .or. bc(3) == DL_MG_BC_PERIODIC) then
         if ( dir > 0 ) then
            k = 0
            do j = sy, ey
               cshift = mod(sx + j + ez + oc, 2)
               do i = sx + cshift, ex, 2
                  k = k + 1
                  sbuff_bt(k,2) = p_ptr(i, j, ez)
               enddo
            enddo
         else
            k = 0
            do j = sy, ey
               cshift = mod(sx + j  + ez + 1 + oc, 2)
               do i = sx + cshift, ex, 2
                  k = k + 1
                  p_ptr(i, j, ez+1) = rbuff_bt(k,2)
               enddo
            enddo
         endif
      endif
      !$OMP  END SINGLE NOWAIT
#else
      ! for no MPI code the data must be shifted around in PBC directions
      integer cshift, i, j, k, oc

      ! it can be done in one go
      if ( dir < 0) return

      oc = 1 - colour
      ! N-S

      !$OMP SINGLE
      if (bc(1) == DL_MG_BC_PERIODIC) then
         do k = sz, ez
            cshift = mod(sx - 1 + sy + k + oc, 2)
            do j = sy + cshift, ey, 2
               p_ptr(sx-1,j,k) = p_ptr(ex,j,k)
            enddo
            !cshift = mod(ex + 1 + sy + k + oc, 2)
            cshift = 1 - cshift
            do j = sy + cshift, ey, 2
               p_ptr(ex+1,j,k) =  p_ptr(sx,j,k)
            enddo
         enddo
      endif
      !$OMP END SINGLE NOWAIT

      !$OMP SINGLE
      if (bc(2) == DL_MG_BC_PERIODIC) then
         do k = sz, ez
            cshift = mod(sx + sy - 1 + k + oc, 2)
            do i = sx + cshift, ex, 2
               p_ptr(i,sy-1,k) = p_ptr(i,ey,k)
            enddo
            cshift = 1- cshift
            do i = sx + cshift, ex, 2
               p_ptr(i,ey+1,k) =  p_ptr(i,sy,k)
            enddo
         enddo
      endif
      !$OMP END SINGLE NOWAIT

      !$OMP SINGLE
      if (bc(3) == DL_MG_BC_PERIODIC) then
         do j = sy, ey
            cshift = mod(sx + j + sz - 1 + oc, 2)
            do i = sx + cshift, ex, 2
               p_ptr(i,j,sz-1) = p_ptr(i,j,ez)
            enddo
            cshift = 1 - cshift
            do i = sx + cshift, ex, 2
               p_ptr(i,j,ez+1) =  p_ptr(i,j,sz)
            enddo
         enddo
      endif
      !$OMP END SINGLE NOWAIT

#endif
      ! there is a omp barrier after the call of this subroutine
    end subroutine buffs_halos_transfers


#ifdef DLMG_RELAX_PRINT

    subroutine mg_relax_debug_print(position)
      use dl_mg_mpi, only :  minmax_grid
      use dl_mg_kernels
      use dl_mg_utils, only : compute_norm
#ifdef DUMP_ARRAYS
      use dl_mg_mpi, only : dump_data
#endif
      implicit none
      character(len=*),   intent(in)  ::  position

      real(wp) xn, yn, zn, cmin(3), cmax(3) !  for debugging
      character(len=1) p_component
      integer comp_task


      comp_task = residual
      p_component = "p"

      select case(position)
      case("top")

         xn = compute_norm(mg(level),"f")
         yn = compute_norm(mg(level),p_component)
         zn = residual_norm(mg(level),eq_type,comp_task)
         call minmax_grid(mg(t),"c",cmin,cmax)
         if ( mg(t)%coords(1) == 0 .and. &
              mg(t)%coords(2) == 0 .and. &
              mg(t)%coords(3) == 0 ) then
            write(0,*) 'relax norm r,p,f', level, zn, yn, xn
            write(0,*) 'relax min  c    ', cmin
            write(0,*) 'relax max  c    ', cmax
         endif
         if ( t == 2) then
!!$         call dump_data(mg(t), "c1", "c1_")
!!$         call dump_data(mg(t), "c3", "c2_")
!!$         call dump_data(mg(t), "c3", "c3_")
#ifdef DUMP_ARRAYS
            call dump_data(mg(t), "p ", "ps_")
#endif
!!$
         end if
!!$      if ( t == 1) then
!!$            write(0,*) 'p(1)222 ',mg(t)%p(2,2,2), mg(t)%w(2,2,2)
!!$         endif

!!$
!!$      write(0,*) 'bc values', level, coords(3), maxval(abs(mg(t)%p(:,:,sz-1))),&
!!$           maxval(abs(mg(t)%p(:,:,sz))),  maxval(abs(mg(t)%p(:,:,ez+1))), maxval(abs(mg(t)%p(:,:,ez)))
!!$            maxval(abs(mg(t)%p(:,sy-1,:))),  maxval(abs(mg(t)%p(:,ey+1,:))),&
!!$            maxval(abs(mg(t)%p(:,:,sz-1))),  maxval(abs(mg(t)%p(:,:,ez+1)))

!!$      do it=1,3
!!$         write(0,*) level, coords(3), it, minval( mg(t)%c(sx:ex,sy:ey,sz:ez+1,it) ),  maxval( mg(t)%c(sx:ex,sy:ey,sz:ez+1,it) )
!!$      enddo
      case("bottom")
         xn = compute_norm(mg(level),p_component)
         yn = residual_norm(mg(level), eq_type, comp_task)
         if ( mg(t)%coords(1) == 0 .and. &
              mg(t)%coords(2) == 0 .and. &
              mg(t)%coords(3) == 0 ) then
            write(0,*) 'relax norm r,p ' , level, yn, xn
         endif
         if (t == 2) then
#ifdef DUMP_ARRAYS
            call dump_data(mg(t), "p ", "pe_")
#endif
            !call error_abort("stop")
         endif
!!$         if ( t == 2) then
!!$            write(0,*) 'p(2)333 ',mg(t)%p(3,3,3)
!!$         endif

      end select
    end subroutine mg_relax_debug_print
#endif

  end subroutine mg_relax


  subroutine mg_restrict(mg, level, eq_type)
    use dl_mg_common_data, only : restriction_weight, blk
    use dl_mg_errors
    use dl_mg_kernels
    use dl_mg_timer
    implicit none
    type(mg_t), intent(inout) :: mg(:)
    integer,    intent(in)    :: level
    integer, intent(in) :: eq_type

    integer  t, iblk, task
    integer sidx(3), eidx(3)

    ! compute residual at current level
    t = level
    if (.not. mg(t)%active ) return

    call mg_timer(start, trestrict, t, tcompute)

    task = RESIDUAL

    !$OMP DO SCHEDULE(STATIC,1)
    do iblk = 1, size(blk(t)%start, dim=2)
       call residual_kernel(task, blk(t), iblk, mg(t), eq_type)
    enddo
    !$OMP ENDDO

    
    if ( mg(t)%aggregate == 1) then
       !$OMP MASTER
       call aggregate_data
       !$OMP END MASTER
       !$OMP BARRIER
    else
       call exchange_halo(t, "r")
    end if
    

    if ( .not. mg(t-1)%active) then
       call mg_timer(stop, trestrict, t, tcompute)
       return
    endif


    if ( eq_type == EQ_PBE_FAS .or. eq_type == EQ_PBE_FAS_NOSTERIC .or. &
         eq_type == EQ_PBE_NEWTON .or. eq_type == EQ_PBE_NEWTON_NOSTERIC ) then
       ! careless initialisation, good for testing, probable not needed
       !$OMP DO SCHEDULE(STATIC,1)
       do iblk = 1, size(blk(t-1)%start, dim=2)
          call blocked_set(mg(t-1), blk(t-1), iblk, 0.0_wp, &
               mg(t-1)%psx, mg(t-1)%psy, mg(t-1)%psz, mg(t-1)%w, &
               lbound( mg(t-1)%w), ubound(mg(t-1)%w))
       enddo
       !$OMP ENDDO NOWAIT
    endif
    !$OMP DO SCHEDULE(STATIC,1)
    do iblk = 1, size(blk(t-1)%start, dim=2)
       call restrict_kernel("r", dl_mg_half_weight_restriction, mg, t, blk, iblk, eq_type)
    enddo
    !$OMP ENDDO

    if (eq_type == EQ_PBE_NEWTON .or. eq_type == EQ_PBE_NEWTON_NOSTERIC ) then
       ! transfer current p_h -> w_H
       !$OMP DO SCHEDULE(STATIC,1)
       do iblk = 1, size(blk(t-1)%start, dim=2)
          call restrict_kernel("w", dl_mg_injection_restriction, mg, t, blk, iblk, eq_type)
       enddo
    endif

    if ( eq_type == EQ_PBE_FAS .or. eq_type == EQ_PBE_FAS_NOSTERIC) then
       ! transfer current p_h -> w_H
       !$OMP DO SCHEDULE(STATIC,1)
       do iblk = 1, size(blk(t-1)%start, dim=2)
          call restrict_kernel("w", dl_mg_injection_restriction, mg, t, blk, iblk, eq_type)
       enddo
       !$OMP ENDDO
       ! probable this halo exchage can be avoided by using p from halos
       ! but needs careful treatment of the boundaries
       !$OMP MASTER
       call exchange_halo(t-1, "w")
       !$OMP END MASTER
       !$OMP BARRIER
       ! build f_H = r_H + N_H(w_H)
       sidx = lbound(mg(t-1)%p)
       eidx = ubound(mg(t-1)%p)
       !$OMP DO SCHEDULE(STATIC,1)
       do iblk = 1, size(blk(t-1)%start, dim=2)
          call residual_kernel(F_NONLINEAR, blk(t-1), iblk, mg(t-1), eq_type)
          call blocked_scaled_copy(mg(t-1), blk(t-1), iblk, &
               mg(t-1)%psx, mg(t-1)%psy, mg(t-1)%psz, 1.0_wp, mg(t-1)%w,&
               mg(t-1)%psx, mg(t-1)%psy, mg(t-1)%psz, mg(t-1)%p, &
               sidx, eidx )
       enddo
       ! the commented out operation from below  is done in the above loop
       !$OMP ENDDO
!!$OMP WORKSHARE
       !mg(t-1)%p(:,:,:) = mg(t-1)%w(:,:,:)
!!$OMP END WORKSHARE
    else
       sidx = lbound(mg(t-1)%p)
       eidx = ubound(mg(t-1)%p)
       !$OMP DO SCHEDULE(STATIC,1)
       do iblk = 1, size(blk(t-1)%start, dim=2)
          call blocked_set(mg(t-1), blk(t-1), iblk, 0.0_wp, &
               mg(t-1)%psx, mg(t-1)%psy, mg(t-1)%psz, mg(t-1)%p, &
               sidx, eidx)
       enddo
       !$OMP ENDDO
!!$         !$OMP WORKSHARE
!!$         mg(t-1)%p(:,:,:) = 0.0_wp
!!$         !$OMP END WORKSHARE
    endif
!!$      if ( t == 2) then
!!$         write(0,*) 'restrict p', t-1, mg(t-1)%p
!!$         write(0,*) '          ', t  , mg(t)%p(3,3,3)
!!$         write(0,*) 'restrict r', t-1, mg(t-1)%r
!!$         write(0,*) "restrict f", t -1, mg(t-1)%f
!!$      endif
    ! initial solution is set to 0 in linear case
!!$      if( .not. use_nonlinear) then
!!$         !$OMP WORKSHARE
!!$         mg(t-1)%p(:,:,:) = 0.0_wp
!!$         !$OMP END WORKSHARE
!!$      endif

    call mg_timer(stop, trestrict, t, tcompute)

  contains

    subroutine aggregate_data
      use dl_mg_mpi
      implicit none
#ifdef MPI
      integer comm, ierr
      integer i, source, dest, n, root
      integer asx, aex, asy, aey, asz, aez
      integer sx, ex, sy, ey, sz, ez
      integer, allocatable :: recvcnts(:), recvdisp(:)
      real(wp), allocatable :: buff(:)

      call mg_timer(start, trestrict, t, tcomm)

      if ( mg(t)%active) then
         call mg_get_idx(mg(t), sx, ex, sy, ey, sz, ez)
         comm = mg(t)%comm
         call mpi_cart_rank(comm,(/0, 0, 0/), root, ierr)
         if ( mg(t)%agg_mast) then
            call mpi_comm_size(comm,n,ierr)
            allocate(recvcnts(0:n-1),recvdisp(0:n-1))
            recvcnts(0) = 0! no need of (ex - sx + 1) * ( ey - sy + 1) * (ez - sz + 1)
            recvdisp(0) = 0
            do i = 1, n - 1
               asx = mg(t)%agg_map(1,i)
               aex = mg(t)%agg_map(2,i)
               asy = mg(t)%agg_map(3,i)
               aey = mg(t)%agg_map(4,i)
               asz = mg(t)%agg_map(5,i)
               aez = mg(t)%agg_map(6,i)
               recvcnts(i) = (aex - asx + 1) * ( aey - asy + 1) * (aez - asz + 1)
               recvdisp(i) = recvdisp(i-1) + recvcnts(i-1)
            enddo
            allocate(buff(sum(recvcnts)))
            call mpi_gatherv(MPI_IN_PLACE,0,MPI_DOUBLE_PRECISION,&
                 buff,recvcnts,recvdisp,MPI_DOUBLE_PRECISION,&
                 root,comm,ierr)
            if ( ierr /= MPI_SUCCESS) then
               write(0,*) 'error in restriction aggregation on root'
            endif
            do i = 1, n - 1
               asx = mg(t)%agg_map(1,i)
               aex = mg(t)%agg_map(2,i)
               asy = mg(t)%agg_map(3,i)
               aey = mg(t)%agg_map(4,i)
               asz = mg(t)%agg_map(5,i)
               aez = mg(t)%agg_map(6,i)
               mg(t)%r(asx:aex, asy:aey, asz:aez) = reshape (buff(recvdisp(i) +1: recvdisp(i)+recvcnts(i)), &
                    (/ aex - asx + 1,  aey - asy + 1, aez - asz + 1 /))
            enddo
!!$              do i = 1, size(mg(t)%agg_map, 2)
!!$                 st = mg(t)%agg_map(1,i)
!!$                 et = mg(t)%agg_map(2,i)
!!$                 n = (ex - sx + 1) * ( ey - sy + 1) * (et - st + 1)
!!$                 source = mg(t)%agg_map(3,i)
!!$                 call mpi_recv(mg(t)%r(sx:ex,sy:ey,st:et), n, &
!!$                      MPI_DOUBLE_PRECISION, source, 33, comm, MPI_STATUS_IGNORE, ierr)
!!$              enddo
         else
            n = (ex - sx + 1) * (ey - sy + 1) * (ez - sz + 1)
            call mpi_gatherv(mg(t)%r(sx:ex,sy:ey,sz:ez),n,MPI_DOUBLE_PRECISION,&
                 mg(t)%p(sx,sy,sz),(/0/),(/0/),MPI_DOUBLE_PRECISION,&
                 root,comm,ierr)
            if ( ierr /= MPI_SUCCESS) then
               write(0,*) 'error in restriction aggregation on sender'
            endif
!!$              dest=mg(t)%agg_map(3,1)
!!$              call mpi_send(mg(t)%r(sx:ex,sy:ey,sz:ez),n,    &
!!$                   MPI_DOUBLE_PRECISION,dest,33,comm,        &
!!$                   ierr)
         endif

!!$           comm = mg(t)%comm
!!$           call get_mpi_grid(comm,rank=myid,coords=coords)
!!$           n=size(mg(t)%r(sx:ex,sy:ey,sz:ez))
!!$
!!$           if (mod(coords(3),2) == 0) then
!!$              call mpi_cart_rank(comm,coords+(/0,0,1/),source,ierr)
!!$              call mpi_recv(mg(t)%r(sx:ex,sy:ey,ez+1:2*ez-sz+1),n, &
!!$                   MPI_DOUBLE_PRECISION,source,33,comm,MPI_STATUS_IGNORE, &
!!$                   ierr)
!!$           else
!!$              call mpi_cart_rank(comm,coords+(/0,0,-1/),dest,ierr)
!!$              call mpi_send(mg(t)%r(sx:ex,sy:ey,sz:ez),n,    &
!!$                   MPI_DOUBLE_PRECISION,dest,33,comm,        &
!!$                   ierr)
!!$           endif
!!$           !         write(0,*)'mgdresrt, aggregate myid ', myid,source,dest
      endif
      call mg_timer(stop, trestrict, t, tcomm)
#endif

    end subroutine aggregate_data


    subroutine exchange_halo(t,component)
      use dl_mg_mpi
      implicit none

      integer, intent(in) :: t
      character(len=1) :: component

#ifdef MPI

      integer i

      call mg_timer(start, trestrict, t, tcomm)

      !$omp master
      do i = 1, 3
         call exchange_halo_begin(i, exch_full, 1, mg(t), component)
         call exchange_halo_end(i, exch_full, 1, mg(t),  component)
      enddo
      !$omp end master
      call mg_timer(stop, trestrict, t, tcomm)
      !$omp barrier
#else
      call pbc_halo(mg(t),component)
#endif
    end subroutine exchange_halo

  end subroutine mg_restrict


  subroutine mg_prolongate(mg, level, eq_type)
    use dl_mg_mpi
    use dl_mg_common_data, only : blk
    use dl_mg_kernels, only : blocked_update
    use dl_mg_timer

    implicit none
    type(mg_t), target, intent(inout) :: mg(:)
    integer,    intent(in)    :: level, eq_type

    integer t, sx, ex, sy, ey, sz, ez, sxc, exc, syc, eyc, szc, ezc, &
         i, j, k, ic, jc, kc, s1, s2, s3, dsx, dex, dsy, dey, dsz, dez
    real(wp), allocatable :: pa(:, :, :)
    real(wp), pointer :: pc_ptr(:,:,:), pf_ptr(:,:,:)
    logical accumulate

    if ( .not. mg(level+1)%active ) return

    pc_ptr => null(); pf_ptr => null()

    t = level
    call mg_timer(start, tprolong, t+1, tcompute)

    call mg_get_idx(mg(t+1), sx, ex, sy, ey, sz, ez)
    call mg_get_idx(mg(t), sxc, exc, syc, eyc, szc, ezc)

    ! one can skip the boundary values for fine grid.  They are fixed in relax subroutine
    select  case ( mg(t+1)%aggregate )
    case(0)

       if ( mg(t)%active) then

!!$          write(0,*) 'prolong', t, sxc, exc, syc, eyc, szc, ezc
!!$          !write(0,*) ' prolong ', use_nonlinear
!!$          if (use_nonlinear) then
!!$             pc_ptr => mg(t)%p
!!$             pf_ptr => mg(t+1)%p
!!$             s1 = mg(t+1)%psx
!!$             s2 = mg(t+1)%psy
!!$             s3 = mg(t+1)%psz
!!$             accumulate = .true.
!!$             !write(0,*) 'before workshare', t, sxc, exc, syc, eyc, szc, ezc
!!$             !$OMP WORKSHARE
!!$
!!$             !$OMP END WORKSHARE
!!$             !write(0,*) 'prolong correction norm', compute_norm(mg(t),"p")
!!$          else

          pc_ptr => mg(t)%p
          pf_ptr  => mg(t+1)%r
          s1 = mg(t+1)%rsx
          s2 = mg(t+1)%rsy
          s3 = mg(t+1)%rsz
          accumulate = .false.
          !write(0,*) 'prolong correction norm', compute_norm(mg(t),"p")

          ! prolongate only the error in FAS
          if (eq_type == EQ_PBE_FAS .or. eq_type == EQ_PBE_FAS_NOSTERIC) then
             call compute_halo_extend(t, dsx, dex, dsy, dey, dsz, dez)
             !$OMP DO SCHEDULE(STATIC,1)
             do i = 1, size(blk(t)%start, dim=2)
                call blocked_update(mg(t), blk(t), i, &
                     mg(t)%psx, mg(t)%psy, mg(t)%psz, -mg(t)%w,&
                     mg(t)%psx, mg(t)%psy, mg(t)%psz, mg(t)%p, &
                     (/ sxc+dsx, syc+dsy, szc+dsz /), &
                     (/ exc+dex, eyc+dey, ezc+dez /))
             enddo
             !$OMP ENDDO
!!$                 !$OMP WORKSHARE
!!$                 mg(t)%p(sxc-dsx:exc+dex, syc-dsy:eyc+dey, szc-dsz:ezc+dez) =&
!!$                      mg(t)%p(sxc-dsx:exc+dex, syc-dsy:eyc+dey, szc-dsz:ezc+dez) &
!!$                      - mg(t)%w(sxc-dsx:exc+dex, syc-dsy:eyc+dey, szc-dsz:ezc+dez)
!!$                 !$OMP END WORKSHARE
          Endif

          call prolongate_loop(s1, s2, s3, mg(t)%p, mg(t+1)%aggregate, accumulate, mg(t+1)%r)

          !if (t == 1) write(44,'(5e15.5)') mg(t+1)%r(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)

!!$           if (use_nonlinear) then
!!$              !$OMP WORKSHARE
!!$              mg(t+1)%p(sx:ex, sy:ey, sz:ez) = p_ptr(sx:ex, sy:ey, sz:ez) + mg(t+1)%p(sx:ex, sy:ey, sz:ez)
!!$              !$OMP END WORKSHARE
!!$           end if


       endif

       !call dump_data(mg(t+1),"r ","r2_")
       !call dump_data(mg(t),"p ","r1_")

       call exchange_halo(t+1)



       !$OMP DO SCHEDULE(STATIC,1)
       do i = 1, size(blk(t+1)%start, dim=2)
          call blocked_update(mg(t+1), blk(t+1), i,&
               mg(t+1)%rsx, mg(t+1)%rsy, mg(t+1)%rsz, mg(t+1)%r,&
               mg(t+1)%psx, mg(t+1)%psy, mg(t+1)%psz, mg(t+1)%p)
       enddo
       !$OMP ENDDO
!!$OMP WORKSHARE
       !mg(t+1)%p(sx:ex, sy:ey, sz:ez) = mg(t+1)%p(sx:ex, sy:ey, sz:ez) + &
       !      mg(t+1)%r(sx:ex, sy:ey, sz:ez)
!!$OMP END WORKSHARE
       !if (t == 1) write(600+mg(t+1)%coords(2)*10+mg(t+1)%coords(3),'(1e15.5)') mg(t+1)%r(sx:ex,sy:ey,sz:ez)


    case(1)

!!! LEFT BEHIND

       !j=maxval(mg(t+1)%agg_map(2,:))
       !allocate (pa(sx-1 : ex+1, sy-1 : ey+1, sz-1 : j + 1) )
       if (mg(t)%active) then

          !mg(t+1)%r = 0.0_wp
          call prolongate_loop(mg(t+1)%rsx, mg(t+1)%rsy, mg(t+1)%rsz, mg(t)%p, mg(t+1)%aggregate,&
               .false., mg(t+1)%r)
          !$OMP WORKSHARE
          mg(t+1)%p(sx:ex,sy:ey,sz:ez) = mg(t+1)%p(sx:ex,sy:ey,sz:ez) + &
               mg(t+1)%r(sx:ex,sy:ey,sz:ez)
          !$OMP END WORKSHARE

       endif

       call distribute_data

    end select
    call mg_timer(stop, tprolong, t+1, tcompute)
  contains

    subroutine prolongate_loop(psx, psy, psz, pc, agg, accumulate, p)
      ! compute correction prolongation ( linear and non-linear )
      use dl_mg_common_data, only : blk
      use dl_mg_kernels

      implicit none
      integer, intent(in) :: psx, psy, psz
      real(wp), intent(in) :: pc(mg(t)%psx:, mg(t)%psy:, mg(t)%psz:)
      integer, intent(in) :: agg
      logical, intent(in) :: accumulate
      real(wp), intent(inout):: p(psx:, psy:, psz:)

      integer i, j, k, ic, jc, kc, dsx, dex, dsy, dey, dsz, dez
      integer iblk, s1, s2, s3, dx, dy, dz, ds, lba(3), lbp(3)
      integer sidx(3), eidx(3)
      real(wp), allocatable, save  :: a(:,:,:)

      ! we need to start form sdc-1 to catch the fine point that might be just before
      ! local boundary. This could be optimised later
      ds = 1

      !$OMP MASTER
      select case (agg)
      case(0)
         allocate(a(sx-2:ex+1, sy-2:ey+1, sz-2:ez+1))
         !a = 0.0_wp Take care, this omission poluates halos with NaN but it shouldn't be
         ! a problem as halos are updated by exchanges in mg_relax
      case(1)
         allocate(a(mg(t+1)%rsx:mg(t+1)%rex, mg(t+1)%rsy:mg(t+1)%rey, mg(t+1)%rsz:mg(t+1)%rez))
      end select
      !$OMP END MASTER
      !$OMP BARRIER

      select case (agg)
      case(0)
         sidx = lbound(a)
         eidx = ubound(a)
         !$OMP DO SCHEDULE(STATIC,1)
         do iblk = 1, size(blk(t+1)%start, dim=2)
            call blocked_set(mg(t+1), blk(t+1), iblk, 0.0_wp, &
                 sx-2, sy-2, sz-2, a, &
                 sidx, eidx)
         enddo
         !$OMP ENDDO
      end select

      !$OMP DO SCHEDULE(STATIC,1)
      do iblk = 1, size(blk(t)%start, dim = 2)
         if ( blk(t)%start(1,iblk)  == sxc) then
            s1 = sxc - ds; dx = ds
         else
            s1 = blk(t)%start(1,iblk); dx = 0
         endif
         if ( blk(t)%start(2,iblk)  == syc) then
            s2 = syc - ds; dy = ds
         else
            s2 = blk(t)%start(2,iblk); dy = 0
         endif
         if ( blk(t)%start(3,iblk)  == szc) then
            s3 = szc - ds; dz = ds
         else
            s3 = blk(t)%start(3,iblk); dz = 0
         endif

         do kc = s3, min(blk(t)%end(3,iblk) + dz, ezc)
            do jc = s2, min(blk(t)%end(2,iblk) + dy, eyc)
               do ic = s1, min(blk(t)%end(1,iblk) + dx, exc)
                  !do kc = szc-1, ezc
                  !do jc = syc-1, eyc
                  !do ic = sxc-1, exc
                  k = 2 * kc - 1;  j = 2 * jc - 1; i = 2 * ic - 1

                  a(i, j, k) =  pc(ic, jc, kc)
                  a(i+1, j, k) =  &
                       0.5_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc))
                  a(i, j+1, k) = &
                       0.5_wp * (pc(ic, jc, kc) + pc(ic, jc+1, kc))
                  a(i, j, k+1) = &
                       0.5_wp * (pc(ic, jc, kc) + pc(ic, jc, kc+1))
                  a(i+1, j+1, k) =  &
                       0.25_wp * (pc(ic, jc, kc) + pc(ic+1, jc, kc) &
                       + pc(ic, jc+1, kc) + pc(ic+1, jc+1, kc))

                  a(i+1, j, k+1) = &
                       0.25_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc) &
                       +pc(ic, jc, kc+1) + pc(ic+1, jc, kc+1))
                  a(i, j+1, k+1) = &
                       0.25_wp*(pc(ic, jc, kc) + pc(ic, jc+1, kc) &
                       +pc(ic, jc, kc+1) + pc(ic, jc+1, kc+1))
                  a(i+1, j+1, k+1) =  &
                       0.125_wp*( pc(ic,   jc,   kc)   + pc(ic+1, jc,   kc) &
                       + pc(ic,   jc+1, kc)   + pc(ic,   jc,   kc+1) &
                       + pc(ic+1, jc+1, kc)   + pc(ic+1, jc,   kc+1) &
                       + pc(ic,   jc+1, kc+1) + pc(ic+1, jc+1, kc+1))

               enddo
            enddo
         enddo
      enddo
      !$OMP ENDDO

      select case (agg)
      case (0)

         call compute_halo_extend(t+1, dsx, dex, dsy, dey, dsz, dez)

         if ( accumulate ) then
            !$OMP WORKSHARE
            p(sx+dsx:ex+dex, sy+dsy:ey+dey, sz+dsz:ez+dez) = &
                 p(sx+dsx:ex+dex, sy+dsy:ey+dey, sz+dsz:ez+dez) + &
                 a(sx+dsx:ex+dex, sy+dsy:ey+dey, sz+dsz:ez+dez)
            !$OMP END WORKSHARE
         else
            lba = lbound(a)
            lbp = lbound(p)
            sidx = (/sx+dsx, sy+dsy, sz+dsz/)
            eidx = (/ex+dex, ey+dey, ez+dez/)
            !$OMP DO SCHEDULE(STATIC,1)
            do i = 1, size(blk(t+1)%start, dim=2)
               call blocked_scaled_copy(mg(t+1), blk(t+1), i,&
                    lba(1), lba(2), lba(3), 1.0_wp, a,&
                    lbp(1), lbp(2), lbp(3), p, &
                    sidx, eidx)
            enddo
            !$OMP ENDDO
!!$            !$OMP WORKSHARE
!!$            p(sx-dsx:ex+dex, sy-dsy:ey+dey, sz-dsz:ez+dez) =  &
!!$                 a(sx-dsx:ex+dex, sy-dsy:ey+dey, sz-dsz:ez+dez)
!!$            !$OMP END WORKSHARE
         endif
      case (1)
         !$OMP WORKSHARE
         p = a
         !$OMP END WORKSHARE
      end select

      !$OMP BARRIER
      !$OMP MASTER
      deallocate(a)
      !$OMP END MASTER
      !$OMP BARRIER

      ! in some cases (inactive processes) the left coarse ghost prolungates to the rank-2
      ! therefore we treat them separately

      ! North face
!!$      ic = sxc - 1
!!$      i  = 2 * ic - 1
!!$      if ( i == sx - 1 ) then
!!$         ! north face ( inside  points )
!!$         do kc = szc-1, ezc
!!$            k = 2 * kc - 1
!!$            do jc = syc, eyc
!!$               j = 2 * jc - 1
!!$                  p(i, j, k) = p(i, j, k) + pc(ic, jc, kc)
!!$               p(i, j+1, k) = p(i, j+1, k) &
!!$                    + 0.5_wp * (pc(ic, jc, kc) + pc(ic, jc+1, kc))
!!$               p(i, j, k+1) = p(i, j, k+1) &
!!$                    + 0.5_wp * (pc(ic, jc, kc) + pc(ic, jc, kc+1))
!!$               p(i, j+1, k+1) = p(i, j+1, k+1) &
!!$                    + 0.25_wp*(pc(ic, jc, kc) + pc(ic, jc+1, kc) &
!!$                    +pc(ic, jc, kc+1) + pc(ic, jc+1, kc+1))
!!$            enddo
!!$         enddo
!!$
!!$         ! the west - north edge
!!$         jc = syc -1
!!$         j = 2 * jc -1
!!$         if ( j == sy - 1 ) then
!!$            ! edge points
!!$            do kc = szc-1, ezc
!!$               k = 2 * kc - 1
!!$               p(i, j, k) = p(i, j, k) + pc(ic, jc, kc)
!!$               p(i, j+1, k) = p(i, j+1, k) &
!!$                    + 0.5_wp * (pc(ic, jc, kc) + pc(ic, jc+1, kc))
!!$               p(i, j, k+1) = p(i, j, k+1) &
!!$                    + 0.5_wp * (pc(ic, jc, kc) + pc(ic, jc, kc+1))
!!$               p(i, j+1, k+1) = p(i, j+1, k+1) &
!!$                    + 0.25_wp*(pc(ic, jc, kc) + pc(ic, jc+1, kc) &
!!$                    +pc(ic, jc, kc+1) + pc(ic, jc+1, kc+1))

!!$               ! we put here the west face bits ( because we are on an edge )
!!$               p(i+1, j, k) = p(i+1, j, k) &
!!$                    + 0.5_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc))
!!$               p(i+1, j, k+1) = p(i+1, j, k+1) &
!!$                    + 0.25_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc) &
!!$                    +pc(ic, jc, kc+1) + pc(ic+1, jc, kc+1))
!!$
!!$               ! and the inside points which are not covered in last loop of this face
!!$               p(i+1, j+1, k) = p(i+1, j+1, k) &
!!$                    + 0.25_wp * (pc(ic, jc, kc) + pc(ic+1, jc, kc) &
!!$                    + pc(ic, jc+1, kc) + pc(ic+1, jc+1, kc))
!!$
!!$               p(i+1, j+1, k+1) = p(i+1, j+1, k+1) &
!!$                    + 0.125_wp*( pc(ic,   jc,   kc)   + pc(ic+1, jc,   kc) &
!!$                    + pc(ic,   jc+1, kc)   + pc(ic,   jc,   kc+1) &
!!$                    + pc(ic+1, jc+1, kc)   + pc(ic+1, jc,   kc+1) &
!!$                    + pc(ic,   jc+1, kc+1) + pc(ic+1, jc+1, kc+1))

!!$            enddo
!!$         endif
!!$
!!$         ! the north - bottom edge
!!$         kc = szc - 1
!!$         k = 2 * kc - 1
!!$         if ( k == sz - 1) then
!!$            do jc = syc, eyc
!!$               j = 2 * jc - 1
!!$               p(i, j, k) = p(i, j, k) + pc(ic, jc, kc)
!!$               p(i, j+1, k) = p(i, j+1, k) &
!!$                    + 0.5_wp * (pc(ic, jc, kc) + pc(ic, jc+1, kc))
!!$               p(i, j, k+1) = p(i, j, k+1) &
!!$                    + 0.5_wp * (pc(ic, jc, kc) + pc(ic, jc, kc+1))
!!$               p(i, j+1, k+1) = p(i, j+1, k+1) &
!!$                    + 0.25_wp*(pc(ic, jc, kc) + pc(ic, jc+1, kc) &
!!$                    +pc(ic, jc, kc+1) + pc(ic, jc+1, kc+1))
!!$
!!$               ! bottom face bits
!!$                p(i+1, j, k) = p(i+1, j, k) &
!!$                    + 0.5_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc))
!!$               p(i+1, j+1, k) = p(i+1, j+1, k) &
!!$                    + 0.25_wp * (pc(ic, jc, kc) + pc(ic+1, jc, kc) &
!!$                    + pc(ic, jc+1, kc) + pc(ic+1, jc+1, kc))
!!$
!!$               ! and the inside points which  are not covered in last loop of this face
!!$               p(i+1, j, k+1) = p(i+1, j, k+1) &
!!$                    + 0.25_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc) &
!!$                    +pc(ic, jc, kc+1) + pc(ic+1, jc, kc+1))
!!$               p(i+1, j+1, k+1) = p(i+1, j+1, k+1) &
!!$                    + 0.125_wp*( pc(ic,   jc,   kc)   + pc(ic+1, jc,   kc) &
!!$                    + pc(ic,   jc+1, kc)   + pc(ic,   jc,   kc+1) &
!!$                    + pc(ic+1, jc+1, kc)   + pc(ic+1, jc,   kc+1) &
!!$                    + pc(ic,   jc+1, kc+1) + pc(ic+1, jc+1, kc+1))
!!$
!!$            enddo
!!$         endif
!!$
!!$         ! the NWB corner
!!$
!!$         jc = syc - 1
!!$         j = 2 * jc - 1
!!$         kc = szc - 1
!!$         k = 2 * kc - 1
!!$         if ( j == sy - 1 .and. k == sz - 1 ) then
!!$            ! NWB point
!!$            p(i, j, k) = p(i, j, k) + pc(ic, jc, kc)
!!$            p(i, j+1, k) = p(i, j+1, k) &
!!$                 + 0.5_wp * (pc(ic, jc, kc) + pc(ic, jc+1, kc))
!!$            p(i, j, k+1) = p(i, j, k+1) &
!!$                 + 0.5_wp * (pc(ic, jc, kc) + pc(ic, jc, kc+1))
!!$            p(i, j+1, k+1) = p(i, j+1, k+1) &
!!$                 + 0.25_wp*(pc(ic, jc, kc) + pc(ic, jc+1, kc) &
!!$                 +pc(ic, jc, kc+1) + pc(ic, jc+1, kc+1))
!!$
!!$            ! bits from west face
!!$            p(i+1, j, k) = p(i+1, j, k) &
!!$                 + 0.5_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc))
!!$            p(i+1, j, k+1) = p(i+1, j, k+1) &
!!$                 + 0.25_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc) &
!!$                 +pc(ic, jc, kc+1) + pc(ic+1, jc, kc+1))
!!$
!!$            ! last bit from bottom face
!!$            p(i+1, j+1, k) = p(i+1, j+1, k) &
!!$                 + 0.25_wp * (pc(ic, jc, kc) + pc(ic+1, jc, kc) &
!!$                 + pc(ic, jc+1, kc) + pc(ic+1, jc+1, kc))
!!$
!!$            ! and the inside point wh is ccovered in last loop of this face
!!$            p(i+1, j+1, k+1) = p(i+1, j+1, k+1) &
!!$                 + 0.125_wp*( pc(ic,   jc,   kc)   + pc(ic+1, jc,   kc) &
!!$                 + pc(ic,   jc+1, kc)   + pc(ic,   jc,   kc+1) &
!!$                 + pc(ic+1, jc+1, kc)   + pc(ic+1, jc,   kc+1) &
!!$                 + pc(ic,   jc+1, kc+1) + pc(ic+1, jc+1, kc+1))
!!$
!!$         endif
!!$      endif
!!$
!!$      ! the rest of the points must be done always
!!$
!!$      do kc = szc-1, ezc
!!$         k = 2 * kc - 1
!!$         do jc = syc, eyc
!!$            j = 2 * jc - 1
!!$            p(i+1, j, k) = p(i+1, j, k) &
!!$                 + 0.5_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc))
!!$            p(i+1, j+1, k) = p(i+1, j+1, k) &
!!$                 + 0.25_wp * (pc(ic, jc, kc) + pc(ic+1, jc, kc) &
!!$                 + pc(ic, jc+1, kc) + pc(ic+1, jc+1, kc))
!!$            p(i+1, j, k+1) = p(i+1, j, k+1) &
!!$                 + 0.25_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc) &
!!$                 +pc(ic, jc, kc+1) + pc(ic+1, jc, kc+1))
!!$            p(i+1, j+1, k+1) = p(i+1, j+1, k+1) &
!!$                 + 0.125_wp*( pc(ic,   jc,   kc)   + pc(ic+1, jc,   kc) &
!!$                 + pc(ic,   jc+1, kc)   + pc(ic,   jc,   kc+1) &
!!$                 + pc(ic+1, jc+1, kc)   + pc(ic+1, jc,   kc+1) &
!!$                 + pc(ic,   jc+1, kc+1) + pc(ic+1, jc+1, kc+1))
!!$         enddo
!!$      enddo
!!$
!!$      ! West face
!!$      jc = syc - 1
!!$      j  = 2 * jc - 1
!!$      if ( j == sy - 1 ) then
!!$         do kc = szc, ezc
!!$            k = 2 * kc - 1
!!$            do ic = sxc, exc
!!$               i = 2 * ic - 1
!!$               p(i, j, k) = p(i, j, k) + pc(ic, jc, kc)
!!$               p(i+1, j, k) = p(i+1, j, k) &
!!$                    + 0.5_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc))
!!$               p(i, j, k+1) = p(i, j, k+1) &
!!$                    + 0.5_wp * (pc(ic, jc, kc) + pc(ic, jc, kc+1))
!!$               p(i+1, j, k+1) = p(i+1, j, k+1) &
!!$                    + 0.25_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc) &
!!$                    + pc(ic, jc, kc+1) + pc(ic+1, jc, kc+1))
!!$            enddo
!!$         enddo
!!$
!!$         ! the west - bottom edge
!!$         kc = szc - 1
!!$         k = 2 * kc - 1
!!$         if ( k == sz - 1 ) then
!!$            do ic = sxc, exc
!!$               i = 2 * ic - 1
!!$               p(i, j, k) = p(i, j, k) + pc(ic, jc, kc)
!!$               p(i+1, j, k) = p(i+1, j, k) &
!!$                    + 0.5_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc))
!!$               p(i, j, k+1) = p(i, j, k+1) &
!!$                    + 0.5_wp * (pc(ic, jc, kc) + pc(ic, jc, kc+1))
!!$               p(i+1, j, k+1) = p(i+1, j, k+1) &
!!$                    + 0.25_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc) &
!!$                    +pc(ic, jc, kc+1) + pc(ic+1, jc, kc+1))
!!$               ! bottom face bits as well
!!$               p(i, j+1, k) = p(i, j+1, k) &
!!$                    + 0.5_wp * (pc(ic, jc, kc) + pc(ic, jc+1, kc))
!!$               p(i+1, j+1, k) = p(i+1, j+1, k) &
!!$                    + 0.25_wp * (pc(ic, jc, kc) + pc(ic+1, jc, kc) &
!!$                    + pc(ic, jc+1, kc) + pc(ic+1, jc+1, kc))
!!$               ! and the inside points
!!$               p(i, j+1, k+1) = p(i, j+1, k+1) &
!!$                    + 0.25_wp*(pc(ic, jc, kc) + pc(ic, jc+1, kc) &
!!$                    +pc(ic, jc, kc+1) + pc(ic, jc+1, kc+1))
!!$               p(i+1, j+1, k+1) = p(i+1, j+1, k+1) &
!!$                    + 0.125_wp*( pc(ic,   jc,   kc)   + pc(ic+1, jc, kc) &
!!$                    + pc(ic,   jc+1, kc)   + pc(ic,   jc,   kc+1) &
!!$                    + pc(ic+1, jc+1, kc)   + pc(ic+1, jc,   kc+1) &
!!$                    + pc(ic,   jc+1, kc+1) + pc(ic+1, jc+1, kc+1))
!!$            enddo
!!$         endif
!!$
!!$         ! Note: W-N edge and NWB corner are already done if North face section
!!$
!!$      endif
!!$
!!$      ! points right to 2* jc -1 face
!!$      do kc = szc, ezc
!!$         k = 2 * kc - 1
!!$         do ic = sxc, exc
!!$            i = 2 * ic - 1
!!$            p(i, j+1, k) = p(i, j+1, k) &
!!$                 + 0.5_wp * (pc(ic, jc, kc) + pc(ic, jc+1, kc))
!!$            p(i+1, j+1, k) = p(i+1, j+1, k) &
!!$                 + 0.25_wp * (pc(ic, jc, kc) + pc(ic+1, jc, kc) &
!!$                 + pc(ic, jc+1, kc) + pc(ic+1, jc+1, kc))
!!$            p(i, j+1, k+1) = p(i, j+1, k+1) &
!!$                 + 0.25_wp*(pc(ic, jc, kc) + pc(ic, jc+1, kc) &
!!$                 +pc(ic, jc, kc+1) + pc(ic, jc+1, kc+1))
!!$            p(i+1, j+1, k+1) = p(i+1, j+1, k+1) &
!!$                 + 0.125_wp*( pc(ic,   jc,   kc)   + pc(ic+1, jc,   kc) &
!!$                 + pc(ic,   jc+1, kc)   + pc(ic,   jc,   kc+1) &
!!$                 + pc(ic+1, jc+1, kc)   + pc(ic+1, jc,   kc+1) &
!!$                 + pc(ic,   jc+1, kc+1) + pc(ic+1, jc+1, kc+1))
!!$         enddo
!!$      enddo
!!$
!!$       ! Bottom face
!!$
!!$      kc = szc - 1
!!$      k = 2 * kc - 1
!!$      if ( k == sz-1) then
!!$         do jc = syc, eyc
!!$            j = 2 * jc - 1
!!$            do ic = sxc, exc
!!$               i = 2 * ic - 1
!!$
!!$               p(i, j, k) = p(i, j, k) + pc(ic, jc, kc)
!!$               p(i+1, j, k) = p(i+1, j, k) &
!!$                    + 0.5_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc))
!!$               p(i, j+1, k) = p(i, j+1, k) &
!!$                    + 0.5_wp * (pc(ic, jc, kc) + pc(ic, jc+1, kc))
!!$               p(i+1, j+1, k) = p(i+1, j+1, k) &
!!$                    + 0.25_wp * (pc(ic, jc, kc) + pc(ic+1, jc, kc) &
!!$                    + pc(ic, jc+1, kc) + pc(ic+1, jc+1, kc))
!!$            enddo
!!$         enddo
!!$      end if
!!$

      ! right to 2*kc-1 all point must be done
!!$      do jc = syc-1, eyc
!!$         j = 2 * jc -1
!!$         do ic = sxc-1, exc
!!$            i = 2 * ic - 1
!!$            p(i, j, k+1) = p(i, j, k+1) &
!!$                    + 0.5_wp * (pc(ic, jc, kc) + pc(ic, jc, kc+1))
!!$            p(i+1, j, k+1) = p(i+1, j, k+1) &
!!$                 + 0.25_wp*(pc(ic, jc, kc) + pc(ic+1, jc, kc) &
!!$                 +pc(ic, jc, kc+1) + pc(ic+1, jc, kc+1))
!!$            p(i, j+1, k+1) = p(i, j+1, k+1) &
!!$                 + 0.25_wp*(pc(ic, jc, kc) + pc(ic, jc+1, kc) &
!!$                 +pc(ic, jc, kc+1) + pc(ic, jc+1, kc+1))
!!$            p(i+1, j+1, k+1) = p(i+1, j+1, k+1) &
!!$                 + 0.125_wp*( pc(ic,   jc,   kc)   + pc(ic+1, jc,   kc) &
!!$                 + pc(ic,   jc+1, kc)   + pc(ic,   jc,   kc+1) &
!!$                 + pc(ic+1, jc+1, kc)   + pc(ic+1, jc,   kc+1) &
!!$                 + pc(ic,   jc+1, kc+1) + pc(ic+1, jc+1, kc+1))
!!$         enddo
!!$      enddo

    end subroutine prolongate_loop


    subroutine exchange_halo(t)
      use dl_mg_mpi
      implicit none

      integer, intent(in) :: t
#ifdef MPI

      ! not active in this version
      integer sx, ex, sy, ey, sz, ez, dsx, dex, dsy, dey, dsz, dez,&
           s1, e1, s2, e2, n, stage, source, dest, ierr


      ! nothing to do
      if ( .not. allocated(mg(t)%prolong_transfers)) return

      call mg_timer(start, tprolong, t, tcomm)

      !$OMP MASTER
      call mg_get_idx(mg(t), sx, ex, sy, ey, sz, ez)

      call compute_halo_extend(t, dsx, dex, dsy, dey, dsz, dez)

      !write(0,*) 'prolong exchange:', t, mg(t)%coords,  mg(t)%prolong_transfers
      !call mpi_barrier(mpi_comm_world,ierr)

      do stage = 1,3

!!$          if ( t == 2) then
!!$             write(0,*) 'stage ', stage, mg(t)%coords, &
!!$                  mg(t)%prolong_transfers(:,stage)
!!$          endif

         if (mg(t)%prolong_transfers(1,stage) == 1) then
            s1 = sy + dsy
            e1 = ey + dey
            s2 = sz + dsz
            e2 = ez + dez
            n = (e1 - s1 + 1)*(e2 - s2 +1)
            dest =mg(t)%neighbor(1)
            call mpi_send(mg(t)%r(sx-1, s1:e1, s2:e2 ), n,&
                 mpi_double_precision, dest, 511,mg(t)%comm ,ierr)
         endif

         if (mg(t)%prolong_transfers(2,stage) == 1) then
            s1 = sy + dsy
            e1 = ey + dey
            s2 = sz + dsz
            e2 = ez + dez
            n = (e1 - s1 + 1)*(e2 - s2 +1)
            dest =mg(t)%neighbor(2)
            call mpi_send(mg(t)%r(ex+1, s1:e1, s2:e2), n, &
                 mpi_double_precision, dest, 512,mg(t)%comm ,ierr)
         endif

         if (mg(t)%prolong_transfers(3,stage) == 1) then
            s1 = sx + dsx
            e1 = ex + dex
            s2 = sz + dsz
            e2 = ez + dez
            n = (e1 - s1 + 1)*(e2 - s2 +1)
            dest =mg(t)%neighbor(3)
            call mpi_send(mg(t)%r(s1:e1, sy-1, s2:e2), n, &
                 mpi_double_precision, dest, 521,mg(t)%comm ,ierr)
         endif

         if (mg(t)%prolong_transfers(4,stage) == 1) then
            s1 = sx + dsx
            e1 = ex + dex
            s2 = sz + dsz
            e2 = ez + dez
            n = (e1 - s1 + 1)*(e2 - s2 + 1)
            dest =mg(t)%neighbor(4)
            call mpi_send(mg(t)%r(s1:e1, ey+1, s2:e2),n, &
                 mpi_double_precision, dest, 522,mg(t)%comm ,ierr)
         endif

         if (mg(t)%prolong_transfers(5,stage) == 1) then
            s1 = sx + dsx
            e1 = ex + dex
            s2 = sy + dsy
            e2 = ey + dey
            n = (e1 - s1 + 1)*(e2 - s2 + 1)
            dest =mg(t)%neighbor(5)
            call mpi_send(mg(t)%r(s1:e1, s2:e2, sz-1),n, &
                 mpi_double_precision, dest, 531,mg(t)%comm ,ierr)
         endif

         if (mg(t)%prolong_transfers(6,stage) == 1) then
            s1 = sx + dsx
            e1 = ex + dex
            s2 = sy + dsy
            e2 = ey + dey
            n = (e1 - s1 + 1)*(e2 - s2 + 1)
            dest =mg(t)%neighbor(6)
            call mpi_send(mg(t)%r(s1:e1, s2:e2, ez+1),n, &
                 mpi_double_precision, dest, 532,mg(t)%comm, ierr)
         endif

         ! receive section
         ! swap the order in each dimension to avoid blocking.
         ! This should be moved to a non-blocking version

         if (mg(t)%prolong_transfers(2,stage) == -1) then
            s1 = sy + dsy
            e1 = ey + dey
            s2 = sz + dsz
            e2 = ez + dez
            n = (e1 - s1 + 1)*(e2 - s2 + 1)
            source = mg(t)%neighbor(2)
            call mpi_recv(mg(t)%r(ex, s1:e1, s2:e2), n,&
                 mpi_double_precision, source,&
                 511,mg(t)%comm, MPI_STATUS_IGNORE, ierr)
         endif

         if (mg(t)%prolong_transfers(1, stage) == -1) then
            s1 = sy + dsy
            e1 = ey + dey
            s2 = sz + dsz
            e2 = ez + dez
            n = (e1 - s1 + 1)*(e2 - s2 + 1)
            source = mg(t)%neighbor(1)
            call mpi_recv(mg(t)%r(sx, s1:e1, s2:e2), n, &
                 mpi_double_precision, source,&
                 512,mg(t)%comm, MPI_STATUS_IGNORE, ierr)
         endif

         if (mg(t)%prolong_transfers(4, stage) == -1) then
            s1 = sx + dsx
            e1 = ex + dex
            s2 = sz + dsz
            e2 = ez + dez
            n = (e1 - s1 + 1)*(e2 - s2 +1)
            source = mg(t)%neighbor(4)
            call mpi_recv(mg(t)%r(s1:e1, ey, s2:e2), n, &
                 mpi_double_precision, source,&
                 521, mg(t)%comm, MPI_STATUS_IGNORE, ierr)
         endif

         if (mg(t)%prolong_transfers(3, stage) == -1) then
            s1 = sx + dsx
            e1 = ex + dex
            s2 = sz + dsz
            e2 = ez + dez
            n = (e1 - s1 + 1)*(e2 - s2 +1)
            source = mg(t)%neighbor(3)
            call mpi_recv(mg(t)%r(s1:e1, sy, s2:e2), n, &
                 mpi_double_precision, source,&
                 522,mg(t)%comm, MPI_STATUS_IGNORE, ierr)
         endif

         if (mg(t)%prolong_transfers(6,stage) == -1) then
            s1 = sx + dsx
            e1 = ex + dex
            s2 = sy + dsy
            e2 = ey + dey
            n = (e1 - s1 + 1)*(e2 - s2 + 1)
            source = mg(t)%neighbor(6)
            call mpi_recv(mg(t)%r(s1:e1, s2:e2, ez), n,&
                 mpi_double_precision, source,&
                 531,mg(t)%comm, MPI_STATUS_IGNORE, ierr)
         endif

         if (mg(t)%prolong_transfers(5,stage) == -1) then
            s1 = sx + dsx
            e1 = ex + dex
            s2 = sy + dsy
            e2 = ey + dey
            n = (e1 - s1 + 1)*(e2 - s2 + 1)
            source = mg(t)%neighbor(5)
            call mpi_recv(mg(t)%r(s1:e1, s2:e2, sz), n, &
                 mpi_double_precision, source,&
                 532,mg(t)%comm , MPI_STATUS_IGNORE, ierr)
         endif

      enddo
      !$OMP END MASTER
      call mg_timer(stop, tprolong, t, tcomm)
      !$OMP BARRIER
#endif
    end subroutine exchange_halo


    subroutine distribute_data
      implicit none
#ifdef MPI
      integer i, root, n, comm, ierr
      integer asx, aex, asy, aey, asz, aez
      integer, allocatable :: sendcnts(:), senddisp(:)
      real(wp), allocatable :: buff(:)

      call mg_timer(start, tprolong, t+1, tcomm)

      !$OMP MASTER
      comm = mg(t+1)%comm
      if ( mg(t+1)%agg_mast ) then
         call mpi_cart_rank(comm,(/0, 0, 0/), root, ierr)
         call mpi_comm_size(comm, n, ierr)
         allocate(sendcnts(0:n-1),senddisp(0:n-1))
         sendcnts(0) = 0 !(ex - sx + 1) * ( ey - sy + 1) * (ez - sz + 1)
         senddisp(0) = 0
         do i=1,n-1
            asx = mg(t+1)%agg_map(1,i)
            aex = mg(t+1)%agg_map(2,i)
            asy = mg(t+1)%agg_map(3,i)
            aey = mg(t+1)%agg_map(4,i)
            asz = mg(t+1)%agg_map(5,i)
            aez = mg(t+1)%agg_map(6,i)
            sendcnts(i) = (aex - asx + 1) * ( aey - asy + 1) * (aez - asz + 1)
            senddisp(i) = senddisp(i-1) + sendcnts(i-1)
         enddo
         allocate(buff(0:sum(sendcnts)-1))
         do i =1, n-1
            asx = mg(t+1)%agg_map(1,i)
            aex = mg(t+1)%agg_map(2,i)
            asy = mg(t+1)%agg_map(3,i)
            aey = mg(t+1)%agg_map(4,i)
            asz = mg(t+1)%agg_map(5,i)
            aez = mg(t+1)%agg_map(6,i)
            buff(senddisp(i):senddisp(i)+sendcnts(i)-1) = &
                 reshape( mg(t+1)%r(asx:aex, asy:aey, asz:aez), &
                 (/ size( mg(t+1)%r(asx:aex, asy:aey, asz:aez)) /) )
         enddo
         call mpi_scatterv(buff, sendcnts, senddisp, MPI_DOUBLE_PRECISION,&
              MPI_IN_PLACE, 0, MPI_DOUBLE_PRECISION, root, comm, ierr)

!!$         do i = 1, n !size(mg(t+1)%agg_map, 2)
!!$            st = mg(t+1)%agg_map(1,i)
!!$            et = mg(t+1)%agg_map(2,i)
!!$            n = (ex - sx + 1) * ( ey - sy + 1) * (et - st + 1)
!!$            dest = mg(t+1)%agg_map(3,i)
!!$            call mpi_send(pa(sx:ex,sy:ey,st:et), n, &
!!$                 MPI_DOUBLE_PRECISION, dest, 43, comm, ierr)
!!$         enddo

      else
         n = (ex - sx + 1) * (ey - sy + 1) * (ez - sz + 1)
         root = mg(t+1)%agg_map(1,1)
         call mpi_scatterv(i, (/0/), (/0/), MPI_DOUBLE_PRECISION,&
              mg(t+1)%r(sx:ex,sy:ey,sz:ez), n, MPI_DOUBLE_PRECISION, root, comm, ierr)
!!$         !allocate(pa(sx:ex, sy:ey, sz:ez))
!!$         source = mg(t+1)%agg_map(3,1) ! can be generalised later
!!$         n = (ex - sx + 1) * (ey - sy + 1) * (ez - sz + 1)
!!$         call mpi_recv(pa(sx:ex,sy:ey,sz:ez),n, &
!!$              MPI_DOUBLE_PRECISION,source,43,comm,MPI_STATUS_IGNORE,ierr)
         mg(t+1)%p(sx:ex,sy:ey,sz:ez) =  &
              mg(t+1)%p(sx:ex,sy:ey,sz:ez)  &
              + mg(t+1)%r(sx:ex,sy:ey,sz:ez)
      end if

!!$        call get_mpi_grid(comm, rank=myid)
!!$
!!$        n=size( mg(t+1)%p(sx:ex,sy:ey,sz:ez) )
!!$
!!$        if (mod(coords(3),2) == 0) then
!!$           call mpi_cart_rank(comm, coords+(/0,0,1/), dest, ierr)
!!$           call mpi_send(pa(sx:ex,sy:ey,ez+1:2*ez-sz+1),n, &
!!$                MPI_DOUBLE_PRECISION, dest, 43, comm, ierr)
!!$        else
!!$           allocate(pa(sx:ex, sy:ey, sz:ez))
!!$           call mpi_cart_rank(comm, coords+(/0,0,-1/), source,ierr)
!!$           call mpi_recv(pa(sx:ex,sy:ey,sz:ez),n, &
!!$                MPI_DOUBLE_PRECISION,source,43,comm,MPI_STATUS_IGNORE,ierr)
!!$           mg(t+1)%p(sx:ex,sy:ey,sz:ez) =  &
!!$                mg(t+1)%p(sx:ex,sy:ey,sz:ez)  &
!!$                + pa(sx:ex,sy:ey,sz:ez)
!!$        endif
!!$

      !$OMP END MASTER
      call mg_timer(stop, tprolong, t, tcomm)
      !$OMP BARRIER
#endif
    end subroutine distribute_data

    subroutine compute_halo_extend(t, dsx, dex, dsy, dey, dsz, dez)
      use dl_mg_common_data, only : bc
      implicit none

      integer, intent(in) :: t
      integer, intent(out) :: dsx, dex, dsy, dey, dsz, dez

      integer dims(3)

      dsx = -1; dex = 1
      dsy = -1; dey = 1
      dsz = -1; dez = 1
      ! for DBC we need the halo points but no boundaries
      call get_mpi_grid(mg(t)%comm,dims=dims)
      if ( mg(t)%coords(1) == 0 .and. &
           bc(1) == DL_MG_BC_DIRICHLET) dsx = 0
      if ( mg(t)%coords(1) == dims(1) -1 .and. &
           bc(1) == DL_MG_BC_DIRICHLET ) dex = 0
      if ( mg(t)%coords(2) == 0 .and. &
           bc(2) == DL_MG_BC_DIRICHLET) dsy = 0
      if ( mg(t)%coords(2) == dims(2) -1 .and. &
           bc(2) == DL_MG_BC_DIRICHLET) dey = 0
      if ( mg(t)%coords(3) == 0 .and. &
           bc(3) == DL_MG_BC_DIRICHLET) dsz = 0
      if ( mg(t)%coords(3) == dims(3) -1 .and. &
           bc(3) == DL_MG_BC_DIRICHLET) dez = 0

    end subroutine compute_halo_extend
  end subroutine mg_prolongate

    subroutine allocate_halo_buffers(mg)
      use dl_mg_common_data, only : mg_levels
      implicit none
      type(mg_t), intent(in) :: mg(:)

      integer t, sx, ex, sy, ey, sz, ez

      allocate(halo_buff(mg_levels))

      do t = mg_levels, 1 , -1
         if ( .not. mg(t)%active) exit
         call mg_get_idx(mg(t), sx, ex, sy, ey, sz, ez)

         allocate(halo_buff(t)%sbuff_ns(((ey - sy + 1)*(ez - sz + 1))/2 + 1, 2), &
              halo_buff(t)%rbuff_ns(((ey - sy + 1)*(ez - sz + 1))/2 + 1, 2), &
              halo_buff(t)%sbuff_we(((ex - sx + 1)*(ez - sz + 1))/2 + 1, 2), &
              halo_buff(t)%rbuff_we(((ex - sx + 1)*(ez - sz + 1))/2 + 1, 2), &
              halo_buff(t)%sbuff_bt(((ey - sy + 1)*(ex - sx + 1))/2 + 1, 2), &
              halo_buff(t)%rbuff_bt(((ey - sy + 1)*(ex - sx + 1))/2 + 1, 2))
         call compute_n_send_recv
      enddo

      contains

      subroutine compute_n_send_recv
          ! computes the number of red and black points on local grid faces
          implicit none

#ifdef MPI
          integer n, ic

          ! North - South
          n = (ey - sy + 1) * (ez - sz +1)
          halo_buff(t)%np_send(1:2,0) = n/2;  halo_buff(t)%np_send(1:2,1) = n/2
          halo_buff(t)%np_recv(1:2,0) = n/2;  halo_buff(t)%np_recv(1:2,1) = n/2
          if ( mod(n , 2) == 1 ) then
             ic = mod(sx+sy+sz,2)
             if ( ic == 0 ) then
                halo_buff(t)%np_send(1,0) = halo_buff(t)%np_send(1,0) + 1
                halo_buff(t)%np_recv(1,1) = halo_buff(t)%np_recv(1,1) + 1
             else
                halo_buff(t)%np_send(1,1) = halo_buff(t)%np_send(1,1) + 1
                halo_buff(t)%np_recv(1,0) = halo_buff(t)%np_recv(1,0) + 1
             endif

             ic = mod(ex+sy+sz,2)
             if ( ic == 0 ) then
                halo_buff(t)%np_send(2,0) =  halo_buff(t)%np_send(2,0) + 1
                halo_buff(t)%np_recv(2,1) =  halo_buff(t)%np_recv(2,1) + 1
             else
                halo_buff(t)%np_send(2,1) =  halo_buff(t)%np_send(2,1) + 1
                halo_buff(t)%np_recv(2,0) =  halo_buff(t)%np_recv(2,0) + 1
             endif

          endif

                   ! West -East
          n = (ex - sx + 1) * (ez - sz +1)
          halo_buff(t)%np_send(3:4,0) = n/2;  halo_buff(t)%np_send(3:4,1) = n/2
          halo_buff(t)%np_recv(3:4,0) = n/2;  halo_buff(t)%np_recv(3:4,1) = n/2
          if ( mod(n, 2) == 1 ) then
             ic = mod(sx+sy+sz,2)
             if ( ic == 0 ) then
                halo_buff(t)%np_send(3,0) = halo_buff(t)%np_send(3,0) + 1
                halo_buff(t)%np_recv(3,1) = halo_buff(t)%np_recv(3,1) + 1
             else
                halo_buff(t)%np_send(3,1) =  halo_buff(t)%np_send(3,1) + 1
                halo_buff(t)%np_recv(3,0) =  halo_buff(t)%np_recv(3,0) + 1
             endif

             ic = mod(sx+ey+sz,2)
             if ( ic == 0 ) then
                halo_buff(t)%np_send(4,0) =  halo_buff(t)%np_send(4,0) + 1
                halo_buff(t)%np_recv(4,1) =  halo_buff(t)%np_recv(4,1) + 1
             else
                halo_buff(t)%np_send(4,1) =  halo_buff(t)%np_send(4,1) + 1
                halo_buff(t)%np_recv(4,0) =  halo_buff(t)%np_recv(4,0) + 1
             endif

          endif

          ! Bottom - Top
          n = (ex - sx + 1) * (ey - sy + 1)
          halo_buff(t)%np_send(5:6,0) = n/2;  halo_buff(t)%np_send(5:6,1) = n/2
          halo_buff(t)%np_recv(5:6,0) = n/2;  halo_buff(t)%np_recv(5:6,1) = n/2
          if ( mod(n , 2) == 1 ) then
             ic = mod(sx+sy+sz,2)
             if ( ic == 0 ) then
                halo_buff(t)%np_send(5,0) = halo_buff(t)%np_send(5,0) + 1
                halo_buff(t)%np_recv(5,1) = halo_buff(t)%np_recv(5,1) + 1
             else
                halo_buff(t)%np_send(5,1) = halo_buff(t)%np_send(5,1) + 1
                halo_buff(t)%np_recv(5,0) = halo_buff(t)%np_recv(5,0) + 1
             endif

             ic = mod(sx+sy+ez,2)
             if ( ic == 0 ) then
                halo_buff(t)%np_send(6,0) = halo_buff(t)%np_send(6,0) + 1
                halo_buff(t)%np_recv(6,1) = halo_buff(t)%np_recv(6,1) + 1
             else
                halo_buff(t)%np_send(6,1) = halo_buff(t)%np_send(6,1) + 1
                halo_buff(t)%np_recv(6,0) = halo_buff(t)%np_recv(6,0) + 1
             endif
          endif
#endif
        end subroutine compute_n_send_recv

   end subroutine allocate_halo_buffers


    subroutine free_halo_buffers(mg)
      use dl_mg_common_data, only : mg_levels
      implicit none

      type(mg_t), intent(in) :: mg(:)
      integer t

      do t = mg_levels, 1 , -1
         if ( .not. mg(t)%active) exit
         deallocate(halo_buff(t)%sbuff_ns, &
              halo_buff(t)%rbuff_ns, &
              halo_buff(t)%sbuff_we, &
              halo_buff(t)%rbuff_we, &
              halo_buff(t)%sbuff_bt, &
              halo_buff(t)%rbuff_bt)
      enddo

      deallocate(halo_buff)

    end subroutine free_halo_buffers

    !> specialisation of the blocked residual for multigrid
    !! no need to exchange halos as this op is done at the end of mg_relax
    subroutine mg_residual(mg, eq_type, halo_exch)
      use dl_mg_types, only : mg_t
      use dl_mg_kernels, only : blocked_residual
      use dl_mg_alloc, only : p_vtype, f_vtype, r_vtype
      implicit none
      type(mg_t), intent (inout) :: mg
      integer, intent(in) :: eq_type
      logical, optional, intent(in) :: halo_exch

       call blocked_residual(mg, eq_type, &
            mg%p, p_vtype, 1.0_wp, mg%f, f_vtype, mg%r, r_vtype, &
            halo_exch=halo_exch)
      
    end subroutine mg_residual
    
  end module dl_mg_multigrid_method


 !! unused bit of code dating most probable from FAS epoch for PBE
!!$    subroutine newton_sor(mg, eq_type, niter_newton, niter_smoother, converged)
!!$      use dl_mg_common_data, only : blk, mg_levels, DL_MG_SUCCESS, nthreads
!!$      use dl_mg_kernels
!!$      implicit none
!!$      type(mg_t), target, intent(inout) :: mg(:)
!!$      real(wp), intent(in) :: tol
!!$      integer, intent(in) :: niter_newton, niter_smoother, eq_type
!!$      integer, intent(out) :: converged
!!$
!!$      integer i, j, t, iblk, eq_tag, conv_smth
!!$      real(wp), allocatable, save :: f0(:,:,:), w0(:,:,:), p0(:,:,:)
!!$      type(mg_t), pointer :: gp => null()
!!$      real(wp)  res_u, res_upv,  res_lin, lam
!!$
!!$      !write(*,*) 'in newton_sor', niter_newton, niter_smoother
!!$      t=1
!!$      gp => mg(t)
!!$
!!$      !$OMP MASTER
!!$      allocate(f0(gp%fsx:gp%fex,gp%fsy:gp%fey,gp%fsz:gp%fez), &
!!$               w0(gp%psx:gp%pex,gp%psy:gp%pey,gp%psz:gp%pez), &
!!$               p0(gp%psx:gp%pex,gp%psy:gp%pey,gp%psz:gp%pez))
!!$      !$OMP END MASTER
!!$      !$OMP BARRIER
!!$      !f0 = mg%f; p0 = mg%p
!!$      !$OMP DO SCHEDULE(STATIC,1)
!!$      do iblk = 1, size(blk(t)%start, dim=2)
!!$         call blocked_scaled_copy(gp, blk(t), iblk, gp%fsx, gp%fsy, gp%fsz,&
!!$              1.0_wp, gp%f, gp%fsx, gp%fsy, gp%fsz, f0)
!!$         if ( mg_levels == 1) then
!!$            call blocked_scaled_copy(gp, blk(t), iblk,&
!!$                 gp%psx, gp%psy, gp%psz, 1.0_wp, gp%p,&
!!$                 gp%psx, gp%psy, gp%psz, gp%w, &
!!$                 (/gp%psx, gp%psy, gp%psz/), (/gp%pex, gp%pey, gp%pez/))
!!$         else
!!$            call blocked_scaled_copy(gp, blk(t), iblk,&
!!$                 gp%psx, gp%psy, gp%psz, 1.0_wp, gp%w,&
!!$                 gp%psx, gp%psy, gp%psz, w0, &
!!$                 (/gp%psx, gp%psy, gp%psz/), (/gp%pex, gp%pey, gp%pez/))
!!$         endif
!!$      enddo
!!$      !$OMP ENDDO
!!$
!!$      do i=1, niter_newton
!!$         ! compute residual
!!$         !$OMP DO SCHEDULE(STATIC,1)
!!$         do iblk = 1, size(blk(t)%start, dim=2)
!!$            call residual_kernel(RESIDUAL, blk(t), iblk, gp, eq_type)
!!$         enddo
!!$         !$OMP ENDDO
!!$
!!$         res_u = compute_norm(gp, "r")
!!$         !write(0,*) 'newton iter', i, compute_norm(gp, "f"), compute_norm(gp, "r")
!!$
!!$         ! move residual in f for mg_relax
!!$         !$OMP DO SCHEDULE(STATIC,1)
!!$         do iblk = 1, size(blk(t)%start, dim=2)
!!$            call blocked_scaled_copy(gp, blk(t), iblk, gp%rsx, gp%rsy, gp%rsz, 1.0_wp, &
!!$              gp%r, gp%fsx, gp%fsy, gp%fsz, gp%f)
!!$
!!$            !zero the whole p array
!!$            call blocked_set(gp, blk(t), iblk, 0.0_wp, gp%psx, gp%psy, gp%psz, gp%p, &
!!$                 (/gp%psx, gp%psy, gp%psz/), (/gp%pex, gp%pey, gp%pez/))
!!$         enddo
!!$         !$OMP ENDDO
!!$
!!$         !write(0,*) 'norms', compute_norm(gp, "f")
!!$
!!$         !write(*,*) 'w', maxval(abs(gp%w))
!!$
!!$         ! solve J e = r
!!$         select case (eq_type)
!!$         case (EQ_PBE_FAS)
!!$            eq_tag = EQ_PBE_NEWTON
!!$         case(EQ_PBE_FAS_NOSTERIC)
!!$            eq_tag = EQ_PBE_NEWTON_NOSTERIC
!!$         end select
!!$
!!$         do j=1, niter_smoother
!!$
!!$            call mg_relax(mg, t, 30, eq_tag) ! redblack uses newton iteration rule
!!$            call mg_convergence_test(gp, t, eq_tag, tconv_newton_sor, j, &
!!$                 niter_smoother, conv_smth)
!!$
!!$            res_lin = compute_norm(gp, "r")
!!$            !if(j == 1) then
!!$!            if ( gp%coords(1) == 0 .and. &
!!$!                 gp%coords(2) == 0 .and. &
!!$!                 gp%coords(3) == 0 ) then
!!$!               write(*,*) 'residuals u, lin', i, j, res_u, res_lin
!!$!            endif
!!$            !endif
!!$
!!$            !write(0,*) 'newton sor smoother', j, conv_smth, smth_tag
!!$            !$OMP BARRIER
!!$
!!$            select case(conv_smth)
!!$            case(DL_MG_SUCCESS)
!!$               exit
!!$            end select
!!$         enddo
!!$
!!$
!!$!         if (conv_smth /= DL_MG_SUCCESS) then
!!$!            !$OMP MASTER
!!$!            if (myid == 0) then
!!$!               write(0,*)"WARNING: smoother in newton_sor fail to converge, ..."
!!$!            endif
!!$!            !$OMP END MASTER
!!$!            !exit
!!$!         endif
!!$
!!$         !$OMP DO SCHEDULE(STATIC,1)
!!$         do iblk = 1, size(blk(t)%start, dim=2)
!!$            call blocked_scaled_copy(gp, blk(t), iblk, gp%fsx, gp%fsy, gp%fsz, 1.0_wp, &
!!$                 f0, gp%fsx, gp%fsy, gp%fsz, gp%f)
!!$            call blocked_scaled_copy(gp, blk(t), iblk, gp%psx, gp%psy, gp%psz, 1.0_wp, &
!!$                 gp%p, gp%psx, gp%psy, gp%psz, p0)
!!$         end do
!!$         !$OMP ENDDO
!!$
!!$         lam = 1.0_wp
!!$         do while ( lam > 1.e-2_wp )
!!$
!!$            !$OMP DO SCHEDULE(STATIC,1)
!!$            do iblk = 1, size(blk(t)%start, dim=2)
!!$               call blocked_scaled_copy(gp, blk(t), iblk, gp%psx, gp%psy, gp%psz, lam, &
!!$                    p0, gp%psx, gp%psy, gp%psz, gp%p)
!!$               call blocked_update(gp, blk(t), iblk, gp%psx, gp%psy, gp%psz, gp%w, &
!!$                    gp%psx, gp%psy, gp%psz, gp%p)
!!$            end do
!!$            !$OMP ENDDO
!!$
!!$
!!$            !$OMP MASTER
!!$            call exchange_full_halos
!!$            !$OMP END MASTER
!!$            !$OMP BARRIER
!!$
!!$            !$OMP DO SCHEDULE(STATIC,1)
!!$            do iblk = 1, size(blk(t)%start, dim=2)
!!$               call residual_kernel(RESIDUAL, blk(t), iblk, gp, eq_type)
!!$            enddo
!!$            !$OMP ENDDO
!!$
!!$            res_upv = compute_norm(gp, "r")
!!$
!!$            !if(j == 1) then
!!$!            res_lin = compute_norm(gp, "r")
!!$!            if ( gp%coords(1) == 0 .and. &
!!$!                 gp%coords(2) == 0 .and. &
!!$!                    gp%coords(3) == 0 ) then
!!$!               write(*,*) 'residuals u, u+v', i, res_u, res_upv, lam
!!$!            endif
!!$            !endif
!!$
!!$            if (res_upv > res_u) then
!!$               lam = 0.75 *lam
!!$               !$OMP DO SCHEDULE(STATIC,1)
!!$               do iblk = 1, size(blk(t)%start, dim=2)
!!$                  call blocked_scaled_copy(gp, blk(t), iblk, gp%psx, gp%psy, gp%psz, 1.0_wp, &
!!$                       p0, gp%psx, gp%psy, gp%psz, gp%p)
!!$               end do
!!$               !$OMP ENDDO
!!$            else
!!$               exit
!!$            endif
!!$         enddo ! search damping parameter loop
!!$
!!$            !$OMP DO SCHEDULE(STATIC,1)
!!$            do iblk = 1, size(blk(t)%start, dim=2)
!!$               call blocked_scaled_copy(gp, blk(t), iblk, gp%psx, gp%psy, gp%psz, 1.0_wp, &
!!$                    gp%p, gp%psx, gp%psy, gp%psz, gp%w)
!!$            end do
!!$            !$OMP ENDDO
!!$
!!$
!!$         if (mg_levels > 1) then
!!$            call mg_convergence_test(gp, t, eq_type, tconv_newton_sor, i, niter_newton, converged)
!!$
!!$            !$OMP BARRIER
!!$            if (converged == DL_MG_SUCCESS) exit
!!$         endif
!!$      enddo
!!$
!!$      if ( mg_levels > 1) then
!!$         !$OMP DO SCHEDULE(STATIC,1)
!!$         do iblk = 1, size(blk(t)%start, dim=2)
!!$            call blocked_scaled_copy(gp, blk(t), iblk,&
!!$                 gp%psx, gp%psy, gp%psz, 1.0_wp, w0,&
!!$                 gp%psx, gp%psy, gp%psz, gp%w, &
!!$                 (/gp%psx, gp%psy, gp%psz/), (/gp%pex, gp%pey, gp%pez/))
!!$         enddo
!!$         !$OMP ENDDO
!!$      endif
!!$
!!$      !$OMP MASTER
!!$      if (allocated(f0)) deallocate(f0, p0)
!!$      if (allocated(w0)) deallocate(w0)
!!$      !$OMP END MASTER
!!$
!!$      contains
!!$
!!$        subroutine exchange_full_halos
!!$          ! this is needed if MPI grid has dimension > 1
!!$          use dl_mg_mpi
!!$          use dl_mg_common_data, only : bc
!!$          implicit none
!!$
!!$#ifdef MPI
!!$
!!$          integer i
!!$
!!$          !call mg_timer(start,trelax,t,tcomm)
!!$
!!$          do i = 1, 3
!!$             call exchange_halo_begin(i,exch_full,1,mg(t),"p")
!!$             call exchange_halo_end(i,exch_full,1,mg(t),"p")
!!$          enddo
!!$
!!$          !call mg_timer(stop,trelax,t,tcomm)
!!$#else
!!$           call pbc_halo(mg(t), "p")
!!$#endif
!!$
!!$
!!$      end subroutine exchange_full_halos
!!$
!!$    end subroutine newton_sor
