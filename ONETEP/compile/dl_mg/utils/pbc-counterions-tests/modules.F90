module params
!  use dl_mg_params, only : wp, pi, elcharge, kB => kboltz
  implicit none

  integer, parameter:: wp = kind(0.d0)
  real(wp), parameter :: pi=4.0_wp*atan(1.0_wp), &
       twopi = 2.0_wp * pi, s2pi = sqrt(twopi), &
       fourpi = 4.0_wp * pi
  !> hartree energy devided by kB in Gauss units, i.e. \f$ e^2/(r_B * K_B ) \f$
  real(wp), parameter :: hartree = 3.1577464761980001721e5_wp

  real(wp), parameter :: max_exp_arg = 70.0_wp

  type pbe_params_t
     integer n ! number of ion species
     real(wp), allocatable :: c(:), q(:) ! concentrations and charges
     real(wp) lam, temp! ...
     real(wp) :: r_steric_weight
     logical fixed_nions, use_steric, linearised, use_fas, use_damping
  end type pbe_params_t

  type(pbe_params_t), target :: pbe_params

  real(wp) xl, yl, zl, dx, dy, dz
  real(wp) :: q_counterions
  real(wp) :: delta_q_cntions

  logical :: periodic(3), use_counterions, use_pbe, use_smart_jellium
  logical :: use_damping ! for defco

  integer nproc, mg_comm
  integer ngx, ngy, ngz ! global grid sizes
  integer sx,ex,sy,ey,sz,ez ! local start end indices
  integer isx, iex, isy,iey, isz,iez ! internal domain ranges (for Dirichlet BC)

  integer fd_order, fd_maxiters, niters

  real(wp) tol_res_rel, tol_res_abs, tol_pot_rel, tol_pot_abs, &
       tol_newton_rel, tol_mg_rel

  logical dump_fields

  integer, parameter :: data_unit=25
end module params


module model_functions
  use params
  implicit none

  type gauss_params_t
     real(wp) sig, r0, scale
  end type gauss_params_t

  type permittivity_params_t
     real(wp) eps0, d, delta
  end type permittivity_params_t


  type(gauss_params_t), target :: gauss_params
  type(permittivity_params_t), target :: eps_params

contains

  function gauss(x,y,z)
    implicit none
    real(wp), intent(in) :: x,y,z
    real(wp) gauss

    real(wp) r2, s
    type(gauss_params_t), pointer :: p

    p => gauss_params

    r2=(x-p%r0)**2+(y-p%r0)**2+(z-p%r0)**2

    s = 1.25_wp * p%sig
    !write(*,*) 'in gauss', p,x,y,z, s2pi

    gauss = p%scale * exp(-r2/(2.0_wp * p%sig**2))/(p%sig*s2pi)**3 !&
         !- p%scale * exp(-r2/(2.0_wp * s**2))/(s*s2pi)**3

  end function gauss


  function eps(x,y,z)
    implicit none
    real(wp), intent(in) :: x,y,z

    real(wp) :: eps

    real(wp) r, r0
    type(permittivity_params_t), pointer :: p !pointer to parameters

    p => eps_params
    r0 = gauss_params%r0
    r = sqrt((x-r0)**2+(y-r0)**2+(z-r0)**2)

    eps = 1.0_wp + (p%eps0-1) * 0.5_wp * (1 + erf((r-p%d)/p%delta))

  end function eps


  subroutine compute_steric_weight(i,j,k,x,y,z,stw)
    use params, only : wp, sx,ex, sy,ey, sz,ez, &
         pbe => pbe_params
    implicit none

    integer, intent(in)  ::  i,j,k
    real(wp), intent(in) :: x,y,z
    real(wp), intent(out):: stw(sx:,sy:,sz:)

    real(wp) r,r0,r_stw

    if (.not. pbe%use_steric .or. .not. use_pbe) return

    r0 = gauss_params%r0
    r_stw = pbe%r_steric_weight
    r = sqrt((x-r0)**2 + (y-r0)**2 + (z-r0)**2)

    stw(i,j,k) = 0.5_wp * (1.0_wp + erf((r-r_stw)/eps_params%delta)) ! it makes sent to have the same steepness


  end subroutine compute_steric_weight


  subroutine grid_sum(a, asum)
#ifdef MPI
    use mpi
#endif
    use params
    implicit none

    real(wp), intent(in) :: a(sx:, sy:, sz:)
    real(wp), intent(out) :: asum

    integer i,j,k, ierr
    real(wp) lsum

    lsum = 0.0_wp
    do k = isz, iez ! use internal points only
       do j = isy, iey
          do i = isx, iex
             lsum =lsum + a(i,j,k)
          enddo
       enddo
    enddo

#ifdef MPI
    call mpi_allreduce(lsum, asum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mg_comm, ierr)
#else
    asum = lsum
#endif

  end subroutine grid_sum


  function capexp(x)
    use params, only : max_exp_arg
    implicit none
    real(wp), intent(in) ::x
    real(wp) capexp

    capexp = exp(min(x, max_exp_arg))

  end function capexp


  subroutine abort(msg)
#ifdef MPI
    use mpi
#endif
    implicit none

    character(len=*), intent(in) :: msg

    integer ierr

    write(0,*) msg

#ifdef MPI
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
    STOP
#endif

  end subroutine abort

end module model_functions


module model_data
  use params, only : mg_comm, wp, sx,ex,sy,ey,sz,ez, use_pbe, &
        pbe_params, periodic
  implicit none

  real(wp), allocatable :: rho(:,:,:), stw(:,:,:), eps_half(:,:,:,:),&
       eps_full(:,:,:), pot(:,:,:), betamu(:)

  real(wp) gamma ! accessible volume fraction

  target stw

contains

  subroutine init_arrays
#ifdef MPI
    use mpi
#endif
    implicit none

    integer myid, ierr, pdims(3), my_coords(3)
    logical dummy(3)

#ifdef MPI
    call mpi_comm_rank(mg_comm,myid,ierr)
    call mpi_cart_get(mg_comm, 3, pdims, dummy, my_coords, ierr)
#else
    myid  = 0
    pdims = (/ 1, 1, 1/)
    my_coords = (/ 0, 0, 0 /)
#endif

    call set_grid_partition

    allocate(pot(sx:ex, sy:ey, sz:ez),&
         rho(sx:ex, sy:ey, sz:ez),&
         eps_half(sx:ex, sy:ey, sz:ez,3),&
         eps_full(sx:ex, sy:ey, sz:ez))

    if (use_pbe) then
       if (pbe_params%use_steric)then
          allocate(stw(sx:ex, sy:ey, sz:ez))
       endif
       allocate(betamu(pbe_params%n))
    endif

  end subroutine init_arrays


  subroutine set_grid_partition
#ifdef MPI
    use mpi
#endif
    use params, only : mg_comm, ngx, ngy, ngz,&
    sx, ex, sy, ey, sz, ez, isx, iex, isy, iey, isz, iez
    implicit none

    integer pgrid_dims(3),  my_coords(3), myid, ierr
    integer nl(3), ngxyz(3), r(3), local_shift(3)
    integer nlx, nly, nlz, i
    logical pgrid_periods(3)

#ifdef MPI
    call mpi_cart_get(mg_comm, 3, pgrid_dims, pgrid_periods, my_coords, ierr)
    call mpi_comm_rank(mg_comm,myid,ierr)
#else
    pgrid_dims = (/ 1, 1, 1/)
    pgrid_periods = (/ .false.,  .false.,  .false. /)
    my_coords = (/ 0, 0, 0 /)
    myid = 0
#endif

    ngxyz = (/ ngx, ngy, ngz /)

    nl = ngxyz/pgrid_dims

    ! test for silly values ?

    ! distribute the remainder to the first r ranks
    do i=1,3

       r(i) = mod(ngxyz(i),pgrid_dims(i))

       if ( my_coords(i) < r(i) ) then
          nl(i) = nl(i)+1
          local_shift(i) = my_coords(i)*nl(i)
       else
          local_shift(i) = r(i) * (nl(i)+1) + (my_coords(i) - r(i)) * nl(i)
       endif

    end do

    !write(0,*) 'local starts', myid, local_shift
    nlx = nl(1); nly = nl(2); nlz = nl(3)

#ifdef MPI

    sx = local_shift(1) + 1
    ex = local_shift(1) + nlx
    sy = local_shift(2) + 1
    ey = local_shift(2) + nly
    sz = local_shift(3) + 1
    ez = local_shift(3) + nlz

#else
    sx = 1; ex = ngx
    sy = 1; ey = ngy
    sz = 1; ez = ngz
#endif

    isx = sx; iex = ex
    if (.not. periodic(1)) then
       if (my_coords(1) == 0) then
          isx = sx+1
       endif
       if (my_coords(1) == pgrid_dims(1) -1) then
          iex = ex-1
       end if
    endif

    isy = sy; iey = ey
    if (.not. periodic(2)) then
       if (my_coords(2) == 0) then
          isy = sy+1
       endif
       if (my_coords(2) == pgrid_dims(2) -1) then
          iey = ey-1
       end if
    endif

    isz = sz; iez = ez
    if (.not. periodic(3)) then
       if (my_coords(3) == 0) then
          isz = sz+1
       endif
       if (my_coords(3) == pgrid_dims(3) -1) then
          iez = ez-1
       end if
    endif

  end subroutine set_grid_partition

end module model_data
