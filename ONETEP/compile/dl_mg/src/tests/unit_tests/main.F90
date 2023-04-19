module tests
  implicit none

contains

  !> contains tests for kernels that
  !! don't need communication between
  !! MPI ranks except for reduce operation
  subroutine test_local_kernels
    implicit none

    call init

    call test_dotproduct

    call test_update2

    call finalise

  end subroutine test_local_kernels

  subroutine init
    use dl_mg, only : dl_mg_init
    use grid_data
    implicit none

    integer bc(3)

    where(periodic)
       bc = DL_MG_BC_PERIODIC
    elsewhere
       bc = DL_MG_BC_DIRICHLET
    end where

        call dl_mg_init(ngx, ngy, ngz, dx, dy, dz, bc, gstart, gend, mg_comm, &
             121,'log_test')

  end subroutine init


  subroutine test_dotproduct
    !$ use omp_lib
    use dl_mg_params
    use dl_mg_grids, only : mg
    use dl_mg_mpi, only : get_myid, get_mpi_grid
    use grid_data
    use dl_mg_kernels
    use dl_mg_types, only : mg_allocate
    implicit none

    real(wp), parameter :: key=1.0_wp
    integer nlvls, t, sx,ex,sy,ey,sz,ez, pdims(3)
    real(wp) dotp(2), q
    real(wp), allocatable :: a(:,:,:), b(:,:,:)

    call mg_allocate(mg, .false., .false.)

    nlvls = size(mg)
    t = nlvls
    sx = mg(t)%sx
    ex = mg(t)%ex
    sy = mg(t)%sy
    ey = mg(t)%ey
    sz = mg(t)%sz
    ez = mg(t)%ez

    mg(t)%p = 0.0_wp
    mg(t)%p(sx:ex,sy:ey,sz:ez) = key

    !write(*,*) lbound(mg(t)%p), ubound(mg(t)%p), mg(t)%sx, mg(t)%ex

    call get_mpi_grid(mg(t)%comm, dims=pdims)

    write(*,*) 'grid global dims', ngx, ngy,ngz
    write(*,*) 'using', omp_get_max_threads(), 'threads'
    write(*,*) 'using', pdims, 'mpi ranks'

     mg(t)%p = 0.0_wp
     mg(t)%p(sx:ex,sy:ey,sz:ez) = key

    call dotproduct_wrapper(key, dotp, &
         a     = mg(t)%p, &
         atype = p_vtype, &
         b     = mg(t)%p, &
         btype = p_vtype, &
         mg    = mg(t))
    if ( get_myid() == 0) then
       write(*,*) 'out omp p * p dot product' , dotp(1) - (key**2 * nix*niy*niz)
       write(*,*) 'in  omp p * p dot product' , dotp(2) - (key**2 * nix*niy*niz)
    end if

    allocate(a(sx-1:ex,sy-3:ey,sz:ez))
    a = 0.0_wp
    a(sx:ex,sy:ey,sz:ez) = key


  end subroutine test_dotproduct


  subroutine dotproduct_wrapper(key, rslt, a, atype, b, btype, scale, exyz, mg)
    use dl_mg_params,  only : wp
    use dl_mg_types,   only : mg_t
    use dl_mg_kernels, only : blocked_dotproduct
    implicit none

    real(wp), intent(in)  :: key
    real(wp), intent(out) :: rslt(2)
    type(mg_t), intent(inout), optional :: mg
    real(wp), intent(in) :: a(:,:,:)
    integer, intent(in), optional :: atype
    real(wp), intent(in), optional :: b(:,:,:)
    integer, intent(in), optional :: btype
    real(wp), intent(in), optional :: scale
    integer, intent(in), optional  :: exyz(:)

    ! call with inside and outside openmp
    if (present(mg)) then
       !$omp parallel default(shared)
       rslt(1) = blocked_dotproduct(mg, a, atype, b=b, btype=btype, scale=scale)
       !$omp end parallel
       rslt(2) = blocked_dotproduct(mg, a, atype, b=b, btype=btype, scale=scale)
    else
       !$omp parallel default(shared)
       rslt(1) = blocked_dotproduct(a, b=b, exyz=exyz, scale=scale)
       !$omp end parallel
       rslt(2) = blocked_dotproduct(a, b=b, exyz=exyz, scale=scale)
    end if
  end subroutine dotproduct_wrapper


  subroutine test_update2
    implicit none

  end subroutine test_update2


  subroutine finalise
    use dl_mg
    implicit none

    call dl_mg_free

  end subroutine finalise

end module tests

!> driver program
program kernels_driver
  use tests
  implicit none

  call initialize

  call test_local_kernels

  call finalize

contains

  subroutine initialize
    use dl_mg_mpi_header
    use dl_mg, only : dl_mg_init
    use grid_data
    implicit none

    integer provided, bc(3), ierr

#ifdef MPI
    call MPI_init_thread(MPI_THREAD_FUNNELED,provided,ierr)
    call MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN, ierr);
    if ( provided /= MPI_THREAD_FUNNELED) then
       write(0,*) "warning, required MPI thread safety level not provided"
    end if
#endif

    call read_input

    call set_grid_partition

  end subroutine initialize


  subroutine read_input
    use dl_mg_mpi_header
    use dl_mg_mpi
    use dl_mg_errors
    use grid_data, only : mg_comm, periodic, &
         ngx, ngy, ngz, periodic
    implicit none
    integer i, ic, j, jj, narg, myid, nproc, ierr, pxyz(3), ngxyz(3)
    character(len=64),allocatable :: file(:)
    character(len=64) :: argval, argval2
    real(wp) dxyz(3)
    logical pperiods(3)

    myid = get_myid(mpi_comm_world)

    narg= command_argument_count()
    ic=1
    i=1
    do
       call get_command_argument(ic, argval)
       select case(argval)
       case('-g','-p','-b')
          jj = 1
          do j=ic+1,ic+3
             call get_command_argument(j, argval2)
             select case(argval)
             case('-g')
                read(argval2,*) ngxyz(jj)
             case('-p')
                read(argval2,*) pxyz(jj)
             case('-b')
                read(argval2,*) periodic(jj)
             end select
             jj = jj+1
          end do
          ic=ic+4
       case default
          call handle_error(DL_MG_ERR_UNSPECIFIED, msg='wrong option')
       end select
       if (ic > narg) exit
    end do


#ifdef MPI
    call mpi_comm_size(mpi_comm_world, nproc, ierr)
    if ( nproc /= pxyz(1) * pxyz(2) * pxyz(3) ) then
       write(0,*) "nproc /= npx * npy * npz !, Aborting ..."
       call mpi_abort(mpi_comm_world,1,ierr)
    end if
#endif


#ifdef MPI
  call mpi_cart_create(mpi_comm_world,3,pxyz,periodic,.false.,mg_comm,ierr)
#endif

  ngx = ngxyz(1); ngy = ngxyz(2); ngz = ngxyz(3)

  end subroutine read_input


  subroutine dotproduct_tests
    !$ use omp_lib
    use dl_mg_params
    use dl_mg_grids, only : mg
    use dl_mg_mpi, only : get_myid, get_mpi_grid
    use grid_data
    use dl_mg_kernels
    use dl_mg_types, only : mg_allocate
    implicit none

    real(wp), parameter :: key=7.0_wp
    integer nlvls, t, sx,ex,sy,ey,sz,ez, pdims(3)
    real(wp) dotp, q
    real(wp), allocatable :: a(:,:,:), b(:,:,:)

    call mg_allocate(mg, .false., .false.)

    nlvls = size(mg)
    t = nlvls
    sx = mg(t)%sx
    ex = mg(t)%ex
    sy = mg(t)%sy
    ey = mg(t)%ey
    sz = mg(t)%sz
    ez = mg(t)%ez

    mg(t)%p = 0.0_wp
    mg(t)%p(sx:ex,sy:ey,sz:ez) = key

    !write(*,*) lbound(mg(t)%p), ubound(mg(t)%p), mg(t)%sx, mg(t)%ex

    call get_mpi_grid(mg(t)%comm, dims=pdims)

    write(*,*) 'grid global dims', ngx, ngy,ngz
    write(*,*) 'using', omp_get_max_threads(), 'threads'
    write(*,*) 'using', pdims, 'mpi ranks'

    !$omp parallel default(shared)
    dotp = blocked_dotproduct(mg(t), mg(t)%p, p_vtype, mg(t)%p, p_vtype)
    !$omp end parallel
    if ( get_myid() == 0) then
       write(*,*) ' p * p dot product' , dotp - (key**2 * nix*niy*niz)
    end if

    dotp = blocked_dotproduct(mg(t), mg(t)%p, p_vtype, mg(t)%p, p_vtype)
    if ( get_myid() == 0) then
       write(*,*) ' p * p dot product' , dotp - (key**2 * nix*niy*niz)
    end if

    q=0.5_wp
    dotp = blocked_dotproduct(mg(t), mg(t)%p, p_vtype, mg(t)%p, p_vtype, scale=q)
    if ( get_myid() == 0) then
       write(*,*) ' p * p dot product' , dotp - (q* key**2 * nix*niy*niz)
    end if

    dotp = blocked_dotproduct(mg(t), mg(t)%p, p_vtype)
    if ( get_myid() == 0) then
       write(*,*) ' p * 1 dot product' , dotp/(key * nix*niy*niz) - 1.0_wp
    end if

    mg(t)%f = 0.0_wp
    mg(t)%f(sx:ex,sy:ey,sz:ez) = key
    dotp = blocked_dotproduct(mg(t), mg(t)%f, f_vtype)
    if ( get_myid() == 0) then
       write(*,*) ' f * 1 dot product' , dotp/(key * nix*niy*niz) - 1.0_wp
    end if

    allocate(a(sx-1:ex,sy-3:ey,sz:ez))
    a = 0.0_wp
    a(sx:ex,sy:ey,sz:ez) = key
    dotp = blocked_dotproduct(a(sx:,sy:,sz:),a(sx:,sy:,sz:),[ex,ey,ez])
    if ( get_myid() == 0) then
       write(*,*) ' a * a dot product' , dotp/(key**2 * nix*niy*niz) - 1.0_wp
    end if

    dotp = blocked_dotproduct(a(sx:,sy:,sz:))
    if ( get_myid() == 0) then
       write(*,*) ' a * 1 dot product' , dotp/(key * nix*niy*niz) - 1.0_wp
    end if

    q=0.9_wp
    !$omp parallel firstprivate(q)
    dotp = blocked_dotproduct(a(sx:,sy:,sz:), exyz=[ex,ey,ez], scale=q)
    !$omp end parallel
    if ( get_myid() == 0) then
       write(*,*) ' q *  a * 1 dot product' , dotp/(q * key * nix*niy*niz) - 1.0_wp
    end if

    mg(t)%p(sx:ex,sy:ey,sz:ez) = 1.0_wp
    mg(t)%f(sx:ex,sy:ey,sz:ez) = 2.0_wp
    call blocked_update2(mg(t), mg(t)%p, p_vtype, 1.0_wp, mg(t)%f, f_vtype)
    dotp = blocked_dotproduct(mg(t), mg(t)%p, p_vtype, mg(t)%f, f_vtype)
    if ( get_myid() == 0) then
       write(*,*) ' a = a + s * b' , dotp/(6.0_wp * nix*niy*niz) - 1.0_wp
    end if

    mg(t)%p(sx:ex,sy:ey,sz:ez) = 1.0_wp
    mg(t)%f(sx:ex,sy:ey,sz:ez) = 2.0_wp
    !$omp parallel
    call blocked_update2(mg(t), mg(t)%p, p_vtype, 1.0_wp, mg(t)%f, f_vtype)
    dotp = blocked_dotproduct(mg(t), mg(t)%p, p_vtype, mg(t)%f, f_vtype)
    !$omp end parallel
    if ( get_myid() == 0) then
       write(*,*) ' omp a = a + s * b' , dotp/(6.0_wp * nix*niy*niz) - 1.0_wp
    end if

  end subroutine dotproduct_tests


  subroutine finalize
    implicit none

    integer ierr

    call mpi_finalize(ierr)

  end subroutine finalize


end program kernels_driver
