!> \brief linear algebra operation over the grids
!!
!! Lucian Anton
!!
module dl_mg_kernels
  use dl_mg_params, only : wp
  use dl_mg_alloc, only : r_vtype, f_vtype, p_vtype, &
       z_vtype, y_vtype, get_lbound 
  implicit none

  ! flags for residual computation
  integer, parameter :: residual=901,  f_nonlinear = 902

  interface blocked_dotproduct
     module procedure blocked_dotproduct_with_mg, blocked_dotproduct_no_mg
  end interface

contains

  subroutine redblack(mg, blk, iblk, omega, f, c, w, p, ic, ex, ey, ez, eq_type)
    use dl_mg_types
    use dl_mg_nonlinear_model, only : kappa2, fnl, derfnl
    use dl_mg_errors
    implicit none

    type(mg_t), intent(in) :: mg
    type(block_list_t), intent(in) :: blk
    real(wp), intent(in) :: omega, f(mg%fsx:, mg%fsy:, mg%fsz:), &
         c(mg%csx:, mg%csy:, mg%csz:, :)
    real(wp), intent(inout) :: p(mg%psx:, mg%psy:, mg%psz:)
    real(wp), intent(in) :: w(mg%psx:, mg%psy:, mg%psz:)
    integer, intent(in) :: iblk, ic, ex, ey, ez
    integer, intent(in) :: eq_type

    integer rb, i, j, k, s1, e1, s2, e2, s3, e3
    real(wp) c7, pa, v

    s1 = blk%start(1, iblk)
    s2 = blk%start(2, iblk)
    s3 = blk%start(3, iblk)

    e1 = min(s1 + blk%dims(1) - 1, ex)
    e2 = min(s2 + blk%dims(2) - 1, ey)
    e3 = min(s3 + blk%dims(3) - 1, ez)

    select case (eq_type)
    case(EQ_POISSON)
       do k = s3, e3
          do j = s2, e2
             rb = mod(s1 + j + k + ic, 2)
             do i = s1 + rb, e1 , 2
                pa = f(i,j,k) &
                   - c(i-1,j,  k,  1) * p(i-1,j  ,k) &
                   - c(i  ,j,  k,  1) * p(i+1,j  ,k) &
                   - c(i,  j-1,k,  2) * p(i  ,j-1,k) &
                   - c(i,  j,  k,  2) * p(i  ,j+1,k) &
                   - c(i,  j,  k-1,3) * p(i  ,j,  k-1) &
                   - c(i,  j,  k,  3) * p(i  ,j,  k+1)

                c7 = c(i-1,j,  k,   1) + c(i  ,j,  k,  1) + &
                     c(i,  j-1,k,   2) + c(i  ,j,  k,  2) + &
                     c(i,  j  ,k-1, 3) + c(i  ,j,  k,  3)
                p(i, j, k) = - omega * pa / c7 + (1.0_wp - omega) * p(i, j, k)
             enddo
          enddo
       enddo
    case (EQ_LINEAR_PBE)
       do k = s3, e3
          do j = s2, e2
             rb = mod(s1 + j + k + ic, 2)
             do i = s1 + rb, e1 , 2
                pa = f(i,j,k) &
                     - c(i-1,j,  k,  1) * p(i-1,j  ,k) &
                     - c(i  ,j,  k,  1) * p(i+1,j  ,k) &
                     - c(i,  j-1,k,  2) * p(i  ,j-1,k) &
                     - c(i,  j,  k,  2) * p(i  ,j+1,k) &
                     - c(i,  j,  k-1,3) * p(i  ,j,  k-1) &
                     - c(i,  j,  k,  3) * p(i  ,j,  k+1)


                c7 = c(i-1,j,  k,   1) + c(i  ,j,  k,  1) + &
                     c(i,  j-1,k,   2) + c(i  ,j,  k,  2) + &
                     c(i,  j  ,k-1, 3) + c(i  ,j,  k,  3) + &
                     kappa2 * mg%d(i,j,k)
                p(i, j, k) = - omega * pa / c7 + (1.0_wp - omega) * p(i,j,k) !- p(i, j, k) / c7
             enddo
          enddo
       enddo
    case(EQ_LINEAR_PBE_NOSTERIC)
       do k = s3, e3
          do j = s2, e2
             rb = mod(s1 + j + k + ic, 2)
             do i = s1 + rb, e1 , 2
                pa = f(i,j,k) &
                     - c(i-1,j,  k,  1) * p(i-1,j  ,k) &
                     - c(i  ,j,  k,  1) * p(i+1,j  ,k) &
                     - c(i,  j-1,k,  2) * p(i  ,j-1,k) &
                     - c(i,  j,  k,  2) * p(i  ,j+1,k) &
                     - c(i,  j,  k-1,3) * p(i  ,j,  k-1) &
                     - c(i,  j,  k,  3) * p(i  ,j,  k+1)


                c7 = c(i-1,j,  k,   1) + c(i  ,j,  k,  1) + &
                     c(i,  j-1,k,   2) + c(i  ,j,  k,  2) + &
                     c(i,  j  ,k-1, 3) + c(i  ,j,  k,  3)   &
                     + kappa2
                p(i, j, k) = - omega * pa / c7 + (1.0_wp - omega) * p(i,j,k) !- p(i, j, k) / c7
             enddo
          enddo
       enddo
    case (EQ_LINEAR_PBE_POT)
       do k = s3, e3
          do j = s2, e2
             rb = mod(s1 + j + k + ic, 2)
             do i = s1 + rb, e1 , 2
                pa = f(i,j,k) &
                     - c(i-1,j,  k,  1) * p(i-1,j  ,k) &
                     - c(i  ,j,  k,  1) * p(i+1,j  ,k) &
                     - c(i,  j-1,k,  2) * p(i  ,j-1,k) &
                     - c(i,  j,  k,  2) * p(i  ,j+1,k) &
                     - c(i,  j,  k-1,3) * p(i  ,j,  k-1) &
                     - c(i,  j,  k,  3) * p(i  ,j,  k+1)


                c7 = c(i-1,j,  k,   1) + c(i  ,j,  k,  1) + &
                     c(i,  j-1,k,   2) + c(i  ,j,  k,  2) + &
                     c(i,  j  ,k-1, 3) + c(i  ,j,  k,  3)   &
                     + mg%w(i,j,k)
                p(i, j, k) = - omega * pa / c7 + (1.0_wp - omega) * p(i,j,k) !- p(i, j, k) / c7
             enddo
          enddo
       enddo
    case(EQ_PBE_FAS)
       do k = s3, e3
          do j = s2, e2
             rb = mod(s1 + j + k + ic, 2)
             do i = s1 + rb, e1 , 2
                v = mg%d(i,j,k)
                pa =  -derfnl(p(i,j,k), v)*p(i,j,k) &
                     + fnl(p(i,j,k), v) &
                     + f(i,j,k) &
                     - c(i-1,j,  k,  1) * p(i-1,j  ,k) &
                     - c(i  ,j,  k,  1) * p(i+1,j  ,k) &
                     - c(i,  j-1,k,  2) * p(i  ,j-1,k) &
                     - c(i,  j,  k,  2) * p(i  ,j+1,k) &
                     - c(i,  j,  k-1,3) * p(i  ,j,  k-1) &
                     - c(i,  j,  k,  3) * p(i  ,j,  k+1)

                c7 = -(c(i-1,j,  k,   1) + c(i  ,j,  k,  1) + &
                     c(i,  j-1,k,   2) + c(i  ,j,  k,  2) + &
                     c(i,  j  ,k-1, 3) + c(i  ,j,  k,  3) + &
                     derfnl(p(i,j,k), v))
                p(i, j, k) = omega * pa / c7 + (1.0_wp - omega) * p(i,j,k)
             enddo
          enddo
       enddo
    case (EQ_PBE_FAS_NOSTERIC)
       do k = s3, e3
          do j = s2, e2
             rb = mod(s1 + j + k + ic, 2)
             do i = s1 + rb, e1 , 2
                !write(0,*) 'fnl ', i,j,k, fnl(p(i,j,k)), derfnl(p(i,j,k))
                pa =  -derfnl(p(i,j,k), 1.0_wp)*p(i,j,k) &
                     + fnl(p(i,j,k), 1.0_wp) &
                     + f(i,j,k) &
                     - c(i-1,j,  k,  1) * p(i-1,j  ,k) &
                     - c(i  ,j,  k,  1) * p(i+1,j  ,k) &
                     - c(i,  j-1,k,  2) * p(i  ,j-1,k) &
                     - c(i,  j,  k,  2) * p(i  ,j+1,k) &
                     - c(i,  j,  k-1,3) * p(i  ,j,  k-1) &
                     - c(i,  j,  k,  3) * p(i  ,j,  k+1)

                c7 = -(c(i-1,j,  k,   1) + c(i  ,j,  k,  1) + &
                     c(i,  j-1,k,   2) + c(i  ,j,  k,  2) + &
                     c(i,  j  ,k-1, 3) + c(i  ,j,  k,  3) + &
                     derfnl(p(i,j,k), 1.0_wp))
                p(i, j, k) = omega * pa / c7 + (1.0_wp - omega) * p(i,j,k)
             enddo
          enddo
       enddo
    case(EQ_PBE_NEWTON)
       do k = s3, e3
          do j = s2, e2
             rb = mod(s1 + j + k + ic, 2)
             do i = s1 + rb, e1 , 2
                v = mg%d(i,j,k)
                pa = + f(i,j,k) &
                     - c(i-1,j,  k,  1) * p(i-1,j  ,k) &
                     - c(i  ,j,  k,  1) * p(i+1,j  ,k) &
                     - c(i,  j-1,k,  2) * p(i  ,j-1,k) &
                     - c(i,  j,  k,  2) * p(i  ,j+1,k) &
                     - c(i,  j,  k-1,3) * p(i  ,j,  k-1) &
                     - c(i,  j,  k,  3) * p(i  ,j,  k+1)

                c7 = -(c(i-1,j,  k,   1) + c(i  ,j,  k,  1) + &
                     c(i,  j-1,k,   2) + c(i  ,j,  k,  2) + &
                     c(i,  j  ,k-1, 3) + c(i  ,j,  k,  3) + &
                     derfnl(w(i,j,k), v))
                p(i, j, k) = omega * pa / c7 + (1.0_wp - omega) * p(i,j,k)
             enddo
          enddo
       enddo
    case(EQ_PBE_NEWTON_NOSTERIC)
       do k = s3, e3
          do j = s2, e2
             rb = mod(s1 + j + k + ic, 2)
             do i = s1 + rb, e1 , 2
                pa = + f(i,j,k) &
                     - c(i-1,j,  k,  1) * p(i-1,j  ,k) &
                     - c(i  ,j,  k,  1) * p(i+1,j  ,k) &
                     - c(i,  j-1,k,  2) * p(i  ,j-1,k) &
                     - c(i,  j,  k,  2) * p(i  ,j+1,k) &
                     - c(i,  j,  k-1,3) * p(i  ,j,  k-1) &
                     - c(i,  j,  k,  3) * p(i  ,j,  k+1)

                c7 = -(c(i-1,j,  k,   1) + c(i  ,j,  k,  1) + &
                     c(i,  j-1,k,   2) + c(i  ,j,  k,  2) + &
                     c(i,  j  ,k-1, 3) + c(i  ,j,  k,  3) + &
                     derfnl(w(i,j,k), 1.0_wp))
                p(i, j, k) = omega * pa / c7 + (1.0_wp - omega) * p(i,j,k)
             enddo
          enddo
       enddo
    case default
       call handle_error(DL_MG_ERR_UNSPECIFIED, msg = "unknown equation type")
    end select

  end subroutine redblack

  !> computes residual for MG, CG or Newton
  !! halo_exch=.false. is meant for calls inside multigrid modules
  !! where we can avoid some copies and halo exchanges (which are already done
  !! by mg_relax)
  subroutine blocked_residual(mg, eq_type, p, p_type, qf, f, f_type, r, r_type, halo_exch)
    use dl_mg_types, only : mg_t
    use dl_mg_mpi, only : exchange_full_halos
    use dl_mg_alloc, only : get_lbound
    use dl_mg_params, only : EQ_POISSON
    implicit none
    type(mg_t), intent(inout) :: mg
    integer, intent(in) :: eq_type, p_type, f_type, r_type
    real(wp), intent(in) :: qf
    real(wp), intent(in) :: p(:,:,:), f(:,:,:)
    real(wp), intent(out) :: r(:,:,:)
    logical,optional, intent(in) :: halo_exch

    real(wp), allocatable, save :: ra(:,:,:), pa(:,:,:)
    target ra
    logical halo_exch_

    if (present(halo_exch)) then
       halo_exch_ = halo_exch
    else
       halo_exch_ = .true.
    end if

    !$OMP MASTER
    allocate(ra(mg%psx:mg%pex, mg%psy:mg%pey, mg%psz:mg%pez), &
         pa(mg%psx:mg%pex, mg%psy:mg%pey, mg%psz:mg%pez))
    mg%z => ra ! for halos exchange
    !$OMP END MASTER
    !$OMP BARRIER

    if (qf == 0.0_wp) then
       call blocked_update2(mg, r, r_type, 0.0_wp, sa=0.0_wp)
    ELSE
       call blocked_scaled_copy2(mg, qf, f, f_type, r, r_type)
    endif

    if (halo_exch_) then
       call blocked_scaled_copy2(mg, 1.0_wp, p, p_type, ra, p_vtype, &
            with_boundary=.true.)
       call exchange_full_halos(mg,"z")
       call blocked_lap(mg, 1.0_wp, ra, p_vtype, pa, p_vtype)
    else
       ! p is mg%p so no copy or halo exchange needed
       ! to be used only inside multigrid
       call blocked_lap(mg, 1.0_wp, p, p_vtype, pa, p_vtype)
    end if
       
    if (eq_type /= EQ_POISSON) then
       if (halo_exch_)then
          call blocked_diagmatvec(mg, eq_type, ra, p_vtype)
       else
          call blocked_diagmatvec(mg, eq_type, ra, p_vtype, p, p_vtype)
       end if
       call blocked_update2(mg, r, r_type,  1.0_wp, ra, p_vtype)
    endif
    call blocked_update2(mg, r, r_type, -1.0_wp, pa, p_vtype)

    !$OMP BARRIER
    !$OMP MASTER
    deallocate(ra, pa)
    !$OMP END MASTER

  end subroutine blocked_residual



  subroutine residual_kernel(compute_task, blk, iblk, mg, eq_type)

    ! computes:
    !  1)  residual for the linear equation
    !      r = f - div [c grad p ]
    !  2)  .........        non-linear equation
    !       r = f - div [c grad p ] + fnl[p]
    !  2') linearised version of 2)
    !  2'') Jacobian for Newton method
    !      r = f - div [c grad p ] + derfnl[p_old] p
    !  3) right hand side at coarse level for non-linear case (FAS)
    !       f = r + ( - alfa div [ c grad p ] + fnl[p] )
    use dl_mg_types
    use dl_mg_nonlinear_model, only : kappa2, fnl, derfnl
    use dl_mg_errors
    implicit none
    integer, intent(in) :: compute_task
    type(mg_t), target, intent(inout) :: mg
    type(block_list_t), intent(in) :: blk
    integer, intent(in) :: iblk
    integer, intent(in) :: eq_type

    integer s1, e1, s2, e2, s3, e3, i, j, k
    real(wp) c7, v
    real(wp), pointer :: p_ptr(:,:,:) !< in parallel region implicit save
                                      !! when initialising will make the
                                      !! the pointer shared

    p_ptr => null()


    p_ptr => mg%p

    s1 = blk%start(1, iblk)
    s2 = blk%start(2, iblk)
    s3 = blk%start(3, iblk)

    e1 = blk%end(1,iblk)
    e2 = blk%end(2,iblk)
    e3 = blk%end(3,iblk)

    select case(compute_task)
    case (RESIDUAL)
       select case(eq_type)
       case(EQ_POISSON)
          do k = s3, e3
             do j = s2, e2
                do i = s1, e1
                   c7 = mg%c(i,j,k,1) + mg%c(i-1,j,   k,1) + &
                        mg%c(i,j,k,2) + mg%c(i  ,j-1 ,k,2) + &
                        mg%c(i,j,k,3) + mg%c(i  ,j   ,k-1,3)

                   mg%r(i,j,k) = mg%f(i,j,k) &
                        + c7               * mg%p(i  ,j,  k  ) &
                        -mg%c(i-1,j,k,  1) * mg%p(i-1,j,  k  ) &
                        -mg%c(i  ,j,k,  1) * mg%p(i+1,j,  k  ) &
                        -mg%c(i,  j-1,k,2) * mg%p(i  ,j-1,k  ) &
                        -mg%c(i,  j,  k,2) * mg%p(i  ,j+1,k  ) &
                        -mg%c(i,  j,k-1,3) * mg%p(i  ,j,  k-1) &
                        -mg%c(i,  j,k,  3) * mg%p(i  ,j,  k+1)
                   !write(0,*) 'p ', i,j,k, mg(t)%p(i,j,k), mg(t)%r(i,j,k)
                enddo
             enddo
          enddo
       case (EQ_LINEAR_PBE)
          do k = s3, e3
             do j = s2, e2
                do i = s1, e1
                   c7 = mg%c(i,j,k,1) + mg%c(i-1,j,   k,1) + &
                        mg%c(i,j,k,2) + mg%c(i  ,j-1 ,k,2) + &
                        mg%c(i,j,k,3) + mg%c(i  ,j   ,k-1,3) + &
                        kappa2 * mg%d(i, j, k)

                   mg%r(i,j,k) = mg%f(i,j,k) &
                        + c7               * mg%p(i  ,j,  k  ) &
                        -mg%c(i-1,j,k,  1) * mg%p(i-1,j,  k  ) &
                        -mg%c(i  ,j,k,  1) * mg%p(i+1,j,  k  ) &
                        -mg%c(i,  j-1,k,2) * mg%p(i  ,j-1,k  ) &
                        -mg%c(i,  j,  k,2) * mg%p(i  ,j+1,k  ) &
                        -mg%c(i,  j,k-1,3) * mg%p(i  ,j,  k-1) &
                        -mg%c(i,  j,k,  3) * mg%p(i  ,j,  k+1)
                   !write(0,*) 'p ', i,j,k, mg(t)%p(i,j,k), mg(t)%r(i,j,k)
                enddo
             enddo
          enddo
       case (EQ_LINEAR_PBE_NOSTERIC)
          do k = s3, e3
             do j = s2, e2
                do i = s1, e1
                   c7 = mg%c(i,j,k,1) + mg%c(i-1,j,   k,1) + &
                        mg%c(i,j,k,2) + mg%c(i  ,j-1 ,k,2) + &
                        mg%c(i,j,k,3) + mg%c(i  ,j   ,k-1,3) + &
                        kappa2

                   mg%r(i,j,k) = mg%f(i,j,k) &
                        + c7               * mg%p(i  ,j,  k  ) &
                        -mg%c(i-1,j,k,  1) * mg%p(i-1,j,  k  ) &
                        -mg%c(i  ,j,k,  1) * mg%p(i+1,j,  k  ) &
                        -mg%c(i,  j-1,k,2) * mg%p(i  ,j-1,k  ) &
                        -mg%c(i,  j,  k,2) * mg%p(i  ,j+1,k  ) &
                        -mg%c(i,  j,k-1,3) * mg%p(i  ,j,  k-1) &
                        -mg%c(i,  j,k,  3) * mg%p(i  ,j,  k+1)
                   !write(0,*) 'p ', i,j,k, mg(t)%p(i,j,k), mg(t)%r(i,j,k)
                enddo
             enddo
          enddo
       case (EQ_LINEAR_PBE_POT)
          do k = s3, e3
             do j = s2, e2
                do i = s1, e1
                   c7 = mg%c(i,j,k,1) + mg%c(i-1,j,   k,1) + &
                        mg%c(i,j,k,2) + mg%c(i  ,j-1 ,k,2) + &
                        mg%c(i,j,k,3) + mg%c(i  ,j   ,k-1,3) + &
                        mg%w(i,j,k)

                   mg%r(i,j,k) = mg%f(i,j,k) &
                        + c7               * mg%p(i  ,j,  k  ) &
                        -mg%c(i-1,j,k,  1) * mg%p(i-1,j,  k  ) &
                        -mg%c(i  ,j,k,  1) * mg%p(i+1,j,  k  ) &
                        -mg%c(i,  j-1,k,2) * mg%p(i  ,j-1,k  ) &
                        -mg%c(i,  j,  k,2) * mg%p(i  ,j+1,k  ) &
                        -mg%c(i,  j,k-1,3) * mg%p(i  ,j,  k-1) &
                        -mg%c(i,  j,k,  3) * mg%p(i  ,j,  k+1)
                   !write(0,*) 'p ', i,j,k, mg(t)%p(i,j,k), mg(t)%r(i,j,k)
                enddo
             enddo
          enddo
       case(EQ_PBE_FAS, EQ_PBE_NEWTON)
          do k = s3, e3
             do j = s2, e2
                do i = s1, e1
                   c7 = mg%c(i,j,k,1) + mg%c(i-1,j,   k,1) + &
                        mg%c(i,j,k,2) + mg%c(i  ,j-1 ,k,2) + &
                        mg%c(i,j,k,3) + mg%c(i  ,j   ,k-1,3)

                   v = mg%d(i,j,k)
                   mg%r(i,j,k) = mg%f(i,j,k) + fnl(p_ptr(i,j,k), v)  &
                        + c7               * p_ptr(i  ,j,  k  ) &
                        -mg%c(i-1,j,k,  1) * p_ptr(i-1,j,  k  ) &
                        -mg%c(i  ,j,k,  1) * p_ptr(i+1,j,  k  ) &
                        -mg%c(i,  j-1,k,2) * p_ptr(i  ,j-1,k  ) &
                        -mg%c(i,  j,  k,2) * p_ptr(i  ,j+1,k  ) &
                        -mg%c(i,  j,k-1,3) * p_ptr(i  ,j,  k-1) &
                        -mg%c(i,  j,k,  3) * p_ptr(i  ,j,  k+1)
                   !write(0,*) 'p ', i,j,k, mg%p(i,j,k), mg%r(i,j,k)
                enddo
             enddo
          enddo
       case (EQ_PBE_FAS_NOSTERIC, EQ_PBE_NEWTON_NOSTERIC)
          do k = s3, e3
             do j = s2, e2
                do i = s1, e1
                   c7 = mg%c(i,j,k,1) + mg%c(i-1,j,   k,1) + &
                        mg%c(i,j,k,2) + mg%c(i  ,j-1 ,k,2) + &
                        mg%c(i,j,k,3) + mg%c(i  ,j   ,k-1,3)

                   mg%r(i,j,k) = mg%f(i,j,k) + fnl(p_ptr(i,j,k), 1.0_wp)  &
                        + c7               * p_ptr(i  ,j,  k  ) &
                        -mg%c(i-1,j,k,  1) * p_ptr(i-1,j,  k  ) &
                        -mg%c(i  ,j,k,  1) * p_ptr(i+1,j,  k  ) &
                        -mg%c(i,  j-1,k,2) * p_ptr(i  ,j-1,k  ) &
                        -mg%c(i,  j,  k,2) * p_ptr(i  ,j+1,k  ) &
                        -mg%c(i,  j,k-1,3) * p_ptr(i  ,j,  k-1) &
                        -mg%c(i,  j,k,  3) * p_ptr(i  ,j,  k+1)
                   !write(0,*) 'p ', i,j,k, mg%p(i,j,k), mg%r(i,j,k)
                enddo
             enddo
          enddo
       case default
             call handle_error(DL_MG_ERR_UNSPECIFIED, &
                  msg ="wrong task parameter in dl_mg_kernels:residual_kernel: &
                  &residual task")
       end select
    case(F_NONLINEAR)
       select case (eq_type)
       case (EQ_PBE_FAS)
          do k = s3, e3
             do j = s2, e2
                do i = s1, e1
                   c7 = mg%c(i,j,k,1) + mg%c(i-1,j,   k,1) + &
                        mg%c(i,j,k,2) + mg%c(i  ,j-1 ,k,2) + &
                        mg%c(i,j,k,3) + mg%c(i  ,j   ,k-1,3)

                   v = mg%d(i,j,k)
                   mg%f(i,j,k) = mg%f(i,j,k) - fnl(mg%w(i,j,k), v)   &
                        - c7               * mg%w(i  ,j,  k  ) &
                        +mg%c(i-1,j,k,  1) * mg%w(i-1,j,  k  ) &
                        +mg%c(i  ,j,k,  1) * mg%w(i+1,j,  k  ) &
                        +mg%c(i,  j-1,k,2) * mg%w(i  ,j-1,k  ) &
                        +mg%c(i,  j,  k,2) * mg%w(i  ,j+1,k  ) &
                        +mg%c(i,  j,k-1,3) * mg%w(i  ,j,  k-1) &
                        +mg%c(i,  j,k,  3) * mg%w(i  ,j,  k+1)
                   !write(0,*) 'p ', i,j,k, mg%p(i,j,k), mg%r(i,j,k)
                enddo
             enddo
          enddo
       case (EQ_PBE_FAS_NOSTERIC)
          do k = s3, e3
             do j = s2, e2
                do i = s1, e1
                   c7 = mg%c(i,j,k,1) + mg%c(i-1,j,   k,1) + &
                        mg%c(i,j,k,2) + mg%c(i  ,j-1 ,k,2) + &
                        mg%c(i,j,k,3) + mg%c(i  ,j   ,k-1,3)

                   mg%f(i,j,k) = mg%f(i,j,k) - fnl(mg%w(i,j,k), 1.0_wp)   &
                        - c7               * mg%w(i  ,j,  k  ) &
                        +mg%c(i-1,j,k,  1) * mg%w(i-1,j,  k  ) &
                        +mg%c(i  ,j,k,  1) * mg%w(i+1,j,  k  ) &
                        +mg%c(i,  j-1,k,2) * mg%w(i  ,j-1,k  ) &
                        +mg%c(i,  j,  k,2) * mg%w(i  ,j+1,k  ) &
                        +mg%c(i,  j,k-1,3) * mg%w(i  ,j,  k-1) &
                        +mg%c(i,  j,k,  3) * mg%w(i  ,j,  k+1)
                   !write(0,*) 'p ', i,j,k, mg%p(i,j,k), mg%r(i,j,k)
                enddo
             enddo
          enddo
       case default
          call handle_error(DL_MG_ERR_UNSPECIFIED, &
               msg = "wrong eq_type for F_NONLINEAR task")
       end select
    case default
          call handle_error(DL_MG_ERR_UNSPECIFIED, &
          msg = "wrong task parameter in dl_mg_kernels:residual_kernel")
    end select

  end subroutine residual_kernel


  subroutine restrict_kernel(component, weight, mg, level, blk, iblk, eq_type)
    use dl_mg_types
    use dl_mg_errors
    implicit none
    type(mg_t), target, intent(inout) :: mg(:)
    type(block_list_t), intent(in) :: blk(:)
    integer, intent(in) :: level, iblk, weight
    character(len=1), intent(in) :: component
    integer,          intent(in) :: eq_type

    integer s1, e1, s2, e2, s3, e3
    integer k1, k2, k3, m1, m2, m3

    s1 = blk(level-1)%start(1, iblk)
    s2 = blk(level-1)%start(2, iblk)
    s3 = blk(level-1)%start(3, iblk)

    e1 = blk(level-1)%end(1, iblk)
    e2 = blk(level-1)%end(2, iblk)
    e3 = blk(level-1)%end(3, iblk)

    select case (component)
    case("r")
       k1=mg(level)%rsx
       k2=mg(level)%rsy
       k3=mg(level)%rsz
       m1=mg(level-1)%fsx
       m2=mg(level-1)%fsy
       m3=mg(level-1)%fsz
       call do_loops(mg(level)%r, k1, k2, k3, &
          mg(level-1)%f, m1, m2, m3 )
    case("w")
       k1=mg(level)%psx
       k2=mg(level)%psy
       k3=mg(level)%psz
       m1=mg(level-1)%psx
       m2=mg(level-1)%psy
       m3=mg(level-1)%psz
       if ( eq_type == EQ_PBE_NEWTON .or. eq_type == EQ_PBE_NEWTON_NOSTERIC) then
          call do_loops(mg(level)%w,k1,k2,k3,&
               mg(level-1)%w,m1,m2,m3)
       else
          call do_loops(mg(level)%p,k1,k2,k3,&
               mg(level-1)%w,m1,m2,m3)
       endif
    case default
          call handle_error(DL_MG_ERR_UNSPECIFIED, &
               msg = "wrong compute parameter in dl_mg_kernels:restrict_kernel")
    end select

  contains

    subroutine do_loops (ain, i1, i2, i3, aout, j1, j2, j3)
      implicit none

      integer, intent(in) :: i1, i2, i3, j1, j2, j3
      real(wp), intent(in) :: ain(i1:,i2:,i3:)
      real(wp), intent(out) :: aout(j1:,j2:,j3:)

      real(wp), parameter :: inv12=1.0_wp/12.0_wp

      integer i,j,k,ic,jc,kc

      select case (weight)
      case (dl_mg_injection_restriction)
         do kc = s3, e3
            k = 2 * kc - 1
            do jc = s2, e2
               j = 2 * jc - 1
               do ic = s1, e1
                  i = 2 * ic - 1
                  aout(ic, jc, kc) = ain(i, j, k)
               enddo
            enddo
         enddo
      case(dl_mg_half_weight_restriction)
         do kc = s3, e3
            k = 2 * kc - 1
            do jc = s2, e2
               j = 2 * jc - 1
               do ic = s1, e1
                  i = 2 * ic - 1
                  aout(ic, jc, kc) = inv12 &
                       * (6.0_wp * ain(i,j,k) &
                       + ain(i-1,j,k) + ain(i+1,j,k) &
                       + ain(i,j-1,k) + ain(i,j+1,k) &
                       + ain(i,j,k-1) + ain(i,j,k+1) )
               enddo
            enddo
         enddo
      case default
         call handle_error(DL_MG_ERR_UNSPECIFIED, &
              msg = "unknown resctriction weight")
      end select

    end subroutine do_loops

  end subroutine restrict_kernel


  subroutine blocked_scaled_copy(mg, blk, iblk, i1, i2, i3, scale, ain, j1, j2, j3, aout, sa, ea)
    use dl_mg_types
    implicit none
    type(mg_t), intent(inout)        :: mg
    type(block_list_t), intent(in)   :: blk
    integer, intent(in)              :: iblk, i1, i2, i3, j1, j2, j3
    real(wp), intent(in)             :: scale, ain(i1:, i2:, i3:)
    real(wp), intent(out)            :: aout(j1:, j2:, j3:)
    integer, optional, intent(in)    :: sa(3), ea(3)

    !locals
    integer s1, s2, s3, e1, e2, e3

    s1 = blk%start(1, iblk)
    s2 = blk%start(2, iblk)
    s3 = blk%start(3, iblk)

    e1 = min(s1 + blk%dims(1) - 1, mg%ex)
    e2 = min(s2 + blk%dims(2) - 1, mg%ey)
    e3 = min(s3 + blk%dims(3) - 1, mg%ez)

    if (present(sa) ) then
       if( s1 == mg%sx) then
          s1 = sa(1)
       endif
       if ( s2 == mg%sy) then
          s2 = sa(2)
       endif
       if ( s3 == mg%sz) then
          s3 = sa(3)
       endif
    endif

    if (present(ea)) then
       if (e1 == mg%ex ) then
          e1 = ea(1)
       endif
       if (e2 == mg%ey) then
          e2 = ea(2)
       endif
       if (e3 == mg%ez) then
          e3 = ea(3)
       endif
    endif

!!$    if (present(ea) .and.present(sa)) then
!!$       write(0,*) mg%level, mg%coords(3), sa, ea
!!$       write(0,*) mg%level, mg%coords(3), s1,e1,s2,e2,s3,e3,mg%sx,mg%ex,mg%sy,mg%ey,mg%sz,mg%ez
!!$    endif

    aout(s1:e1, s2:e2, s3:e3) = scale * ain(s1:e1, s2:e2, s3:e3)

  end subroutine blocked_scaled_copy


  subroutine blocked_scaled_copy2(mg, s, vin, vin_type, vout, vout_type, si, ei, nobarrier, with_boundary)
    use dl_mg_types, only : mg_t
    use dl_mg_errors, only : handle_error
    use dl_mg_params, only : DL_MG_ERR_UNSPECIFIED, DL_MG_BC_DIRICHLET
    use dl_mg_common_data, only : bc
    use dl_mg_mpi, only : get_mpi_grid
    implicit none
    type(mg_t), intent(inout) :: mg
    integer, intent(in) :: vin_type, vout_type
    real(wp), intent(in) :: vin(:,:,:), s
    real(wp), intent(out):: vout(:,:,:)
    integer, optional, intent(in)    :: si(3), ei(3)
    logical, optional, intent(in) :: nobarrier
    logical, optional, intent(in) :: with_boundary !< if true copies
    !! the boundary values in case if DBC

    integer vi(3), vo(3), vi1, vi2, vi3, vo1, vo2, vo3, i
    integer as(3), ae(3), ds(3), de(3), coords(3), dims(3)
    logical nobarr, with_bd

    if (present(si) .or. present(ei)) then
       call handle_error(DL_MG_ERR_UNSPECIFIED, &
       msg="blocked_scaled_copy2: arguments si, ei are dezactivated")
    end if
    
    nobarr = .false.
    if (present(nobarrier)) then
       nobarr = nobarrier
    end if
    
    call get_lbound(mg, vin_type, vi)
    call get_lbound(mg, vout_type, vo)

    vi1 = vi(1); vi2=vi(2); vi3=vi(3)
    vo1 = vo(1); vo2=vo(2); vo3=vo(3)

    with_bd = .false.
    if (present(with_boundary)) then
       with_bd = with_boundary
       if (with_bd) then
          call get_mpi_grid(mg%comm, coords=coords, dims=dims)
          as = [ mg%sx, mg%sy, mg%sz ]
          ae = [ mg%ex, mg%ey, mg%ez ]
          ds = 0
          de = 0
          do i = 1,3
             if ( bc(i) == DL_MG_BC_DIRICHLET) then
                if (coords(i) == 0 ) then
                   ds(i) = -1
                end if
                if (coords(i) == dims(i) -1) then
                   de(i) =  1
                end if
             end if
          end do
       end if
    end if
    
    call do_scaled_copy(vi1,vi2,vi3,vin,vo1,vo2,vo3, vout)

    if (nobarr) return
    
    !$omp barrier

  contains

    subroutine do_scaled_copy(vi1,vi2,vi3,vin, vo1,vo2,vo3,vout)
      use dl_mg_common_data, only : blk
      implicit none
      integer, intent(in) :: vi1, vi2, vi3, vo1, vo2, vo3
      real(wp), intent(in) :: vin(vi1:, vi2:, vi3:)
      real(wp), intent(out) :: vout(vo1:, vo2:, vo3:)

      integer r(3), e(3), t, iblk

      t = mg%level

      !$OMP DO SCHEDULE(STATIC,1)
      do iblk = 1,  blk(t)%nblks(2)
         r = blk(t)%start(:, iblk)
         e = blk(t)%end(:,iblk)

         if (with_bd) then
            do i = 1, 3
               if (r(i) == as(i)) then
                  r(i) = r(i) + ds(i)
               endif
               if (e(i) == ae(i)) then
                  e(i) = e(i) + de(i)
               end if
            end do
         end if

!!$         ! capture the boundary in case DBC 
!!$         if (present(si) ) then
!!$            if( s1 == mg%sx) then
!!$               s1 = si(1)
!!$            endif
!!$            if ( s2 == mg%sy) then
!!$               s2 = si(2)
!!$            endif
!!$            if ( s3 == mg%sz) then
!!$               s3 = si(3)
!!$            endif
!!$         endif
!!$
!!$         if (present(ei)) then
!!$            if (e1 == mg%ex ) then
!!$               e1 = ei(1)
!!$            endif
!!$            if (e2 == mg%ey) then
!!$               e2 = ei(2)
!!$            endif
!!$            if (e3 == mg%ez) then
!!$               e3 = ei(3)
!!$            endif
!!$         endif

            vout(r(1):e(1), r(2):e(2), r(3):e(3)) = &
                 s * vin(r(1):e(1), r(2):e(2), r(3):e(3))
      enddo
      !$OMP ENDDO NOWAIT

     end subroutine do_scaled_copy
   end subroutine blocked_scaled_copy2


  subroutine blocked_update(mg, blk, iblk, i1, i2, i3, ain, j1, j2, j3, aout, sa, ea)
    use dl_mg_types
    implicit none
    type(mg_t), intent(inout)      :: mg
    type(block_list_t), intent(in) :: blk
    integer, intent(in)            :: iblk, i1, i2, i3, j1, j2, j3
    integer, optional, intent(in)  :: sa(3), ea(3)
    real(wp), intent(in)           :: ain(i1:, i2:, i3:)
    real(wp), intent(inout)        :: aout(j1:, j2:, j3:)

    !locals
    integer s1, s2, s3, e1, e2, e3

    s1 = blk%start(1, iblk)
    s2 = blk%start(2, iblk)
    s3 = blk%start(3, iblk)
    if (present(sa) ) then
       if( s1 == mg%sx) then
          s1 = sa(1)
       endif
       if ( s2 == mg%sy) then
          s2 = sa(2)
       endif
       if ( s3 == mg%sz) then
          s3 = sa(3)
       endif
    endif

    e1 = min(blk%start(1, iblk) + blk%dims(1) - 1, mg%ex)
    e2 = min(blk%start(2, iblk) + blk%dims(2) - 1, mg%ey)
    e3 = min(blk%start(3, iblk) + blk%dims(3) - 1, mg%ez)

    if (present(ea)) then
       if (e1 == mg%ex ) then
          e1 = ea(1)
       endif
       if (e2 == mg%ey) then
          e2 = ea(2)
       endif
       if (e3 == mg%ez) then
          e3 = ea(3)
       endif
    endif

!!$    if (present(ea) .and.present(sa)) then
!!$       write(0,*) mg%level, mg%coords(3), sa, ea
!!$       write(0,*) mg%level, mg%coords(3), s1,e1,s2,e2,s3,e3,mg%sx,mg%ex,mg%sy,mg%ey,mg%sz,mg%ez
!!$    endif

    aout(s1:e1, s2:e2, s3:e3) = aout(s1:e1, s2:e2, s3:e3) + ain(s1:e1, s2:e2, s3:e3)

  end subroutine blocked_update

  !> coputes a = sa * a + s * b 
  subroutine blocked_update2(mg, a, a_type, s, b, b_type, sa, nobarrier)
    !$ use omp_lib
    use dl_mg_types, only : mg_t
    use dl_mg_common_data, only : blk
    implicit none

    type(mg_t), intent(inout)      :: mg
    integer, intent(in)            :: a_type
    real(wp), intent(inout)        :: a(:,:,:)
    real(wp), intent(in)           :: s
    real(wp), intent(in), optional :: b(:,:,:)
    integer,  intent(in), optional :: b_type
    real(wp), intent(in), optional :: sa
    logical, optional, intent(in)  :: nobarrier 
    
    integer ba(3), bb(3)
    logical in_parallel, nobarr
    real(wp) qa

    call get_lbound(mg, a_type, ba)
    if (present(b)) then
       if (present(b_type)) then
          call get_lbound(mg, b_type, bb)
       else
          bb = ba
       end if
    end if

    if (present(sa)) then
       qa=sa
    else
       qa=1.0_wp
    endif

    nobarr = .false.
    if (present(nobarrier)) then
       nobarr = nobarrier
    end if    
    
    in_parallel = .false.
    !$ in_parallel = omp_in_parallel()

    if ( in_parallel) then
       call loops(mg, a, s, b, qa)
       if (nobarr) return
       !$omp barrier
    else
       !$omp parallel default(shared)
       call loops(mg, a, s, b, qa)
       !$omp end parallel
    end if

  contains

    subroutine loops(mg, a, q, b, qa)
      implicit none
      type(mg_t), intent(inout)      :: mg
      real(wp), intent(inout)        :: a(ba(1):,ba(2):,ba(3):)
      real(wp), intent(in)           :: q, qa
      real(wp), intent(in), optional :: b(bb(1):,bb(2):,bb(3):)

      integer iblk, t
      integer, pointer :: s(:), e(:)

      t = mg%level

      if (present(b)) then
         !$OMP DO SCHEDULE(STATIC,1)
         do iblk = 1, size(blk(t)%start, dim=2)
            s => blk(t)%start(:, iblk)
            e => blk(t)%end(:, iblk)
            a(s(1):e(1), s(2):e(2), s(3):e(3)) &
                 = qa * a(s(1):e(1), s(2):e(2), s(3):e(3)) &
                 + q * b(s(1):e(1), s(2):e(2), s(3):e(3))

         enddo
         !$OMP ENDDO NOWAIT
      else
         if (qa == 0.0_wp) then
            !$OMP DO SCHEDULE(STATIC,1)
            do iblk = 1, size(blk(t)%start, dim=2)
               s => blk(t)%start(:, iblk)
               e => blk(t)%end(:, iblk)
               a(s(1):e(1), s(2):e(2), s(3):e(3)) = q
            enddo
            !$OMP ENDDO NOWAIT
         else
            !$OMP DO SCHEDULE(STATIC,1)
            do iblk = 1, size(blk(t)%start, dim=2)
               s => blk(t)%start(:, iblk)
               e => blk(t)%end(:, iblk)
               a(s(1):e(1), s(2):e(2), s(3):e(3)) &
                    = qa * a(s(1):e(1), s(2):e(2), s(3):e(3)) + q
               
            enddo
            !$OMP ENDDO NOWAIT
         end if
      endif

    end subroutine loops

  end subroutine blocked_update2


  subroutine blocked_scale(mg, s, a, a_type)
    use dl_mg_types, only : mg_t
    use dl_mg_common_data, only : blk
    implicit none

    type(mg_t), intent(inout) :: mg
    integer, intent(in) :: a_type
    real(wp), intent(inout) :: a(:,:,:)
    real(wp), intent(in) :: s

    integer, dimension(3) ::  exyz, sa, ea, ba
    integer i, iblk, t

    call get_lbound(mg, a_type, ba)
    exyz = (/ mg%ex, mg%ey, mg%ez /)

    t = mg%level

    !$OMP DO SCHEDULE(STATIC,1)
    do iblk = 1, size(blk(t)%start, dim=2)
       do i =1, 3
          sa(i) = blk(t)%start(i, iblk) - ba(i) +1
          ea(i) = min(blk(t)%start(i, iblk) + blk(t)%dims(i) - 1, exyz(i)) -ba(i)+1
       enddo
       a(sa(1):ea(1), sa(2):ea(2), sa(3):ea(3)) = &
            s * a(sa(1):ea(1), sa(2):ea(2), sa(3):ea(3))

    enddo
    !$OMP ENDDO NOWAIT

  end subroutine blocked_scale


  subroutine blocked_set(mg, blk, iblk, val, i1, i2, i3, aout, sa, ea)
    use dl_mg_types
    implicit none
    type(mg_t), intent(inout)        :: mg
    type(block_list_t), intent(in)   :: blk
    integer, intent(in)              :: iblk, i1, i2, i3
    real(wp), intent(in)             :: val
    real(wp), intent(out)            :: aout(i1:, i2:, i3:)
    integer, optional, intent(in)    :: sa(3), ea(3)

    !locals
    integer s1, s2, s3, e1, e2, e3

    s1 = blk%start(1, iblk)
    s2 = blk%start(2, iblk)
    s3 = blk%start(3, iblk)

    e1 = blk%end(1, iblk)
    e2 = blk%end(2, iblk)
    e3 = blk%end(3, iblk)

    if (present(sa) ) then
       if( s1 == mg%sx) then
          s1 = sa(1)
       endif
       if ( s2 == mg%sy) then
          s2 = sa(2)
       endif
       if ( s3 == mg%sz) then
          s3 = sa(3)
       endif
    endif

    if (present(ea)) then
       if (e1 == mg%ex ) then
          e1 = ea(1)
       endif
       if (e2 == mg%ey) then
          e2 = ea(2)
       endif
       if (e3 == mg%ez) then
          e3 = ea(3)
       endif
    endif

    aout(s1:e1, s2:e2, s3:e3) = val

  end subroutine blocked_set


  function blocked_dotproduct_with_mg(mg, a, atype, b, btype, scale)
    !$ use omp_lib
    use dl_mg_mpi_header
    use dl_mg_types, only : mg_t
    use dl_mg_common_data, only : blk
    implicit none
    type(mg_t), intent(in) :: mg
    integer, intent(in) :: atype
    integer, intent(in), optional :: btype
    real(wp), intent(in) :: a(:,:,:)
    real(wp), intent(in), optional :: b(:,:,:)
    real(wp), intent(in), optional :: scale

    real(wp) blocked_dotproduct_with_mg

    real(wp), save :: xloc, xworld
    integer, dimension(3) ::  exyz_, ba, bb
    real(wp) scl
    integer i, iblk, t, ierr, is, ie, comm
    logical in_parallel

    call get_lbound(mg, atype, ba)
    if (present(b)) then
       if (present(btype)) then
          call get_lbound(mg, btype, bb)
       else
          bb = ba
       end if
    endif

    exyz_ = (/ mg%ex, mg%ey, mg%ez /)

    t = mg%level
    comm = mg%comm

#include "implementations/blocked_dotproduct_1.inc"

    blocked_dotproduct_with_mg = xworld

  contains

#include "implementations/blocked_dotproduct_2.inc"

  end function blocked_dotproduct_with_mg


  !> simplification of blocked_dotproduct to be used with top grid
  !! assumes that the array start are provided via the array indexes
  !!  such as s(sx:,sy:,sz:)
  !! if exyz is not provided the upper bound of a array will be used
  function blocked_dotproduct_no_mg(a, b, exyz, scale)
    !$ use omp_lib
    use dl_mg_common_data, only : mg_comm, mg_levels
    use dl_mg_mpi_header
    use dl_mg_common_data, only : blk
    implicit none
    real(wp), intent(in) :: a(:,:,:)
    real(wp), intent(in), optional :: b(:,:,:)
    integer, intent(in), optional :: exyz(:) !< end indices for both arrays
    real(wp), intent(in), optional :: scale
    real(wp) blocked_dotproduct_no_mg

    real(wp), save :: xloc, xworld
    integer, dimension(3) ::  ba(3), bb(3), exyz_(3)
    integer i, iblk, ierr, is, ie, t, comm
    real(wp) scl
    logical in_parallel

    !> \todo test for size of exyz

    if (present(exyz))then
       exyz_ = exyz(1:3)
    else
       exyz_ = ubound(a)
    endif

    if (present(scale)) then
       scl = scale
    else
       scl = 1.0_wp
    end if

    t = mg_levels
    comm = mg_comm
    ba(:) = blk(t)%start(:,1)
    exyz_ = ba + exyz_ -1
    if (present(b)) then
       bb = ba
    end if

#include "implementations/blocked_dotproduct_1.inc"

    blocked_dotproduct_no_mg = xworld

  contains

#include "implementations/blocked_dotproduct_2.inc"

  end function blocked_dotproduct_no_mg


  subroutine blocked_lap(mg, scale, vin, vin_type, vout, vout_type, nobarrier)
    use dl_mg_types, only : mg_t
    use dl_mg_common_data, only : blk
    implicit none
    type(mg_t), intent(inout) :: mg
    integer, intent(in) :: vin_type, vout_type
    real(wp), intent(in) :: scale, vin(:,:,:)
    real(wp), intent(out) :: vout(:,:,:)
    logical, optional, intent(in) :: nobarrier

    integer bi(3), bo(3), vi1, vi2, vi3, vo1, vo2, vo3
    logical nobarr

    nobarr=.false.
    if (present(nobarrier)) then
      nobarr=nobarrier
    endif

    call get_lbound(mg, vin_type, bi)
    call get_lbound(mg, vout_type, bo)

    vi1 = bi(1); vi2 = bi(2); vi3 = bi(3)
    vo1 = bo(1); vo2 = bo(2); vo3 = bo(3)

    ! should exchange the halos here ?

    call compute_matvec(vi1, vi2, vi3, vin, vo1,vo2,vo3, vout)
    
    if (nobarr) return
    !$omp barrier

  contains

    subroutine compute_matvec(vi1,vi2,vi3, vin, vo1, vo2, vo3, vout)
      implicit none
      integer, intent(in)   :: vi1, vi2, vi3, vo1, vo2, vo3
      real(wp), intent(in ) :: vin(vi1:, vi2:, vi3:)
      real(wp), intent(out) :: vout(vo1:, vo2:, vo3:)

      integer sv(3), ev(3), exyz(3), i, j, k, iblk, t
      real(wp) c7, vnew

      exyz = (/ mg%ex, mg%ey, mg%ez /)
      t = mg%level

      !$OMP DO SCHEDULE(STATIC,1)
      do iblk = 1, size(blk(t)%start, dim=2)
         do i =1, 3
            sv(i) = blk(t)%start(i, iblk)
            ev(i) = blk(t)%end(i, iblk)
         enddo

         do k = sv(3), ev(3)
            do j = sv(2), ev(2)
               do i = sv(1), ev(1)
                  c7 = mg%c(i,j,k,1) + mg%c(i-1,j,   k,1) + &
                       mg%c(i,j,k,2) + mg%c(i  ,j-1 ,k,2) + &
                       mg%c(i,j,k,3) + mg%c(i  ,j   ,k-1,3)

                  vnew =  -c7 * vin(i  ,j,  k  ) &
                       +mg%c(i-1,j,k,  1) * vin(i-1,j,  k  ) &
                       +mg%c(i  ,j,k,  1) * vin(i+1,j,  k  ) &
                       +mg%c(i,  j-1,k,2) * vin(i  ,j-1,k  ) &
                       +mg%c(i,  j,  k,2) * vin(i  ,j+1,k  ) &
                       +mg%c(i,  j,k-1,3) * vin(i  ,j,  k-1) &
                       +mg%c(i,  j,k,  3) * vin(i  ,j,  k+1)
                  vout(i,j,k) = scale * vnew
                  enddo
               enddo
            enddo
         enddo
         !$OMP ENDDO NOWAIT

    end subroutine compute_matvec

  end subroutine blocked_lap

  
  subroutine blocked_diagmatvec(mg, eq_type, v, v_type, vi, vi_type)
    use dl_mg_types, only : mg_t
    use dl_mg_params, only : EQ_POISSON, EQ_LINEAR_PBE_NOSTERIC, &
         EQ_LINEAR_PBE, EQ_LINEAR_PBE_POT, EQ_PBE_FAS, EQ_PBE_FAS_NOSTERIC, &
         EQ_PBE_NEWTON_NOSTERIC, EQ_PBE_NEWTON
    use dl_mg_errors
    implicit none

    type(mg_t), intent(inout) :: mg
    integer, intent(in) :: eq_type, v_type
    real(wp), intent(inout) :: v(:,:,:)
    integer, optional, intent(in) :: vi_type
    real(wp), optional, intent(in) :: vi(:,:,:) 

    integer bv(3), bi(3)
    logical use_vi

    use_vi=.false.
    if (present(vi)) then
       use_vi = .true.
    end if

!!$    !$OMP MASTER
!!$    if ( check_assertion((eq_type == EQ_PBE_NEWTON .or. &
!!$         eq_type == EQ_PBE_NEWTON_NOSTERIC) .and. &
!!$         .not. (present(u) .and. present(u_type)))) then
!!$       call handle_error(DL_MG_ERR_UNSPECIFIED, &
!!$            msg = "kernel diagmatvec : missing optional argument needed for &
!!$            & Newton case")
!!$    endif
!!$    !$OMP END MASTER

    call get_lbound(mg, v_type, bv)
    !if ( present(u)) then
    !   call get_lbound(mg, u_type, bu)
    !end if
    if ( use_vi) then
       call get_lbound(mg, vi_type, bi)
    end if

    call compute_diagmatvec(v,vi)

  contains

    subroutine compute_diagmatvec(v,vi)
      use dl_mg_common_data, only : blk
      use dl_mg_nonlinear_model
      implicit none
      real(wp), intent(inout) ::         v(bv(1):,bv(2):,bv(3):)
      real(wp), optional, intent(in) :: vi(bi(1):,bi(2):,bi(3):)

      integer i, j,k, t, iblk, exyz(3), sv(3), ev(3)
      real(wp) s

      exyz = (/ mg%ex, mg%ey, mg%ez /)
      t = mg%level

      !$OMP DO SCHEDULE(STATIC,1)
      do iblk = 1, blk(t)%nblks(2)
         do i =1, 3
            sv(i) = blk(t)%start(i, iblk)
            ev(i) = blk(t)%end(i, iblk)
         enddo

         select case(eq_type)
         case(EQ_POISSON)
            do k = sv(3), ev(3)
               do j = sv(2), ev(2)
                  do i = sv(1), ev(1)
                     v(i,j,k) =  0.0_wp
                  enddo
               enddo
            enddo
         case (EQ_LINEAR_PBE)
            if (use_vi)then
               do k = sv(3), ev(3)
                  do j = sv(2), ev(2)
                     do i = sv(1), ev(1)
                        v(i,j,k) =  kappa2 * mg%d(i, j, k) * vi(i  ,j,  k  )
                     enddo
                  enddo
               enddo
            else
               do k = sv(3), ev(3)
                  do j = sv(2), ev(2)
                     do i = sv(1), ev(1)
                        v(i,j,k) =  kappa2 * mg%d(i, j, k) * v(i  ,j,  k  )
                     enddo
                  enddo
               enddo
            end if
         case (EQ_LINEAR_PBE_NOSTERIC)
            if (use_vi)then
               do k = sv(3), ev(3)
                  do j = sv(2), ev(2)
                     do i = sv(1), ev(1)
                        v(i,j,k) = kappa2 * vi(i  ,j,  k  )
                     enddo
                  enddo
               enddo
            else
               do k = sv(3), ev(3)
                  do j = sv(2), ev(2)
                     do i = sv(1), ev(1)
                        v(i,j,k) = kappa2 * v(i  ,j,  k  )
                     enddo
                  enddo
               enddo
            endif
         case (EQ_LINEAR_PBE_POT)
            if(use_vi)then
               do k = sv(3), ev(3)
                  do j = sv(2), ev(2)
                     do i = sv(1), ev(1)
                        v(i,j,k) = mg%w(i,j,k) * vi(i  ,j,  k  )
                     enddo
                  enddo
               enddo
            else
               do k = sv(3), ev(3)
                  do j = sv(2), ev(2)
                     do i = sv(1), ev(1)
                        v(i,j,k) = mg%w(i,j,k) * v(i  ,j,  k  )
                     enddo
                  enddo
               enddo
            endif
         case(EQ_PBE_FAS, EQ_PBE_NEWTON)
            if(use_vi)then
               do k = sv(3), ev(3)
                  do j = sv(2), ev(2)
                     do i = sv(1), ev(1)
                        s = mg%d(i,j,k)
                        v(i,j,k) =  fnl(vi(i,j,k), s)
                     enddo
                  enddo
               enddo
            else
               do k = sv(3), ev(3)
                  do j = sv(2), ev(2)
                     do i = sv(1), ev(1)
                        s = mg%d(i,j,k)
                        v(i,j,k) =  fnl(v(i,j,k), s)
                     enddo
                  enddo
               enddo
            end if
         case (EQ_PBE_FAS_NOSTERIC,EQ_PBE_NEWTON_NOSTERIC)
            if(use_vi)then
               do k = sv(3), ev(3)
                  do j = sv(2), ev(2)
                     do i = sv(1), ev(1)
                        s = 1.0_wp
                        v(i,j,k) = fnl(vi(i,j,k), s)
                     enddo
                  enddo
               enddo
            else
               do k = sv(3), ev(3)
                  do j = sv(2), ev(2)
                     do i = sv(1), ev(1)
                        s = 1.0_wp
                        v(i,j,k) = fnl(v(i,j,k), s)
                     enddo
                  enddo
               enddo
            endif
         case default
            call handle_error(DL_MG_ERR_UNSPECIFIED, &
                 msg = "wrong task parameter in dl_mg_kernels:residual_kernel:&
                 & residual task")
         end select
      enddo
      !$OMP ENDDO

    end subroutine compute_diagmatvec

  end subroutine blocked_diagmatvec



end module dl_mg_kernels


! could be useful later. Full weight restriction didn't work with
! the current setup
!!$#if 0
!!$    case(dl_mg_full_weight_restriction)
!!$       ! This is not used in current version
!!$       t = level
!!$       do kc = s3, e3
!!$          k = 2 * kc - 1
!!$          do jc = s2, e2
!!$             j = 2 * jc - 1
!!$             do ic = s1, e1
!!$                i = 2 * ic - 1
!!$                !                mg(t-1)%f(ic,jc,kc) = mg(t)%r(i,j,k)
!!$
!!$                mg(t-1)%f(ic,jc,kc) = 0.015625_wp * ( 8.0_wp*mg(t)%r(i,j,k)  &
!!$                     + mg(t)%r(i-1,j-1,k-1) &
!!$                     + mg(t)%r(i+1,j-1,k-1) &
!!$                     + mg(t)%r(i-1,j+1,k-1) &
!!$                     + mg(t)%r(i-1,j-1,k+1) &
!!$                     + mg(t)%r(i+1,j+1,k-1) &
!!$                     + mg(t)%r(i+1,j-1,k+1) &
!!$                     + mg(t)%r(i-1,j+1,k+1) &
!!$                     + mg(t)%r(i+1,j+1,k+1) & ! end corners
!!$                     + 2.0_wp*( mg(t)%r(i,j-1,k-1) & ! lower edges
!!$                     + mg(t)%r(i,j+1,k-1) &
!!$                     + mg(t)%r(i-1,j,k-1) &
!!$                     + mg(t)%r(i+1,j,k-1) )&
!!$                     + 2.0_wp*( mg(t)%r(i,j-1,k+1) & ! upper edges
!!$                     + mg(t)%r(i,j+1,k+1) &
!!$                     + mg(t)%r(i-1,j,k+1) &
!!$                     + mg(t)%r(i+1,j,k+1) )&
!!$                     + 2.0_wp*( mg(t)%r(i-1,j-1,k) & ! middle corners
!!$                     + mg(t)%r(i-1,j+1,k) &
!!$                     + mg(t)%r(i+1,j-1,k) &
!!$                     + mg(t)%r(i+1,j+1,k) ) &
!!$                     + 4.0_wp*(mg(t)%r(i-1,j,k) &  ! face centers
!!$                     + mg(t)%r(i+1,j,k) &
!!$                     + mg(t)%r(i,j-1,k) &
!!$                     + mg(t)%r(i,j+1,k) &
!!$                     + mg(t)%r(i,j,k-1) &
!!$                     + mg(t)%r(i,j,k+1) ) )

!!$                !  write(0,*) ic,i,jc,j,kc,k, mg(t)%r(i,j,k), mg(t)%r(i-1,j,k)
!!$             enddo
!!$          enddo
!!$       enddo
!!$#endif
