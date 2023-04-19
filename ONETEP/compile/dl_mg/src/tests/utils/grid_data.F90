module grid_data
  use dl_mg_params
  implicit none
  integer mg_comm
  integer ngx, ngy, ngz, gstart(3), gend(3)
  integer nix, niy, niz ! global inner points
  real(wp) dx,dy,dz
  logical periodic(3)
  target ngx,ngy,ngz,nix,niy,niz

contains

  subroutine set_grid_partition
    use dl_mg_mpi_header


    integer pgrid_dims(3),  my_coords(3), myid, ierr
    integer nl(3), ngxyz(3), r(3), local_shift(3)
    integer nlx, nly, nlz, i, sx, ex, sy, ey, sz, ez
    logical pgrid_periods(3)
    integer, pointer :: nip, ngp
    real(wp) dxyz(3)

#ifdef MPI
    call mpi_cart_get(mg_comm, 3, pgrid_dims, pgrid_periods, my_coords, ierr)
    call mpi_comm_rank(mg_comm,myid,ierr)
#else
    pgrid_dims = (/ 1, 1, 1/)
    pgrid_periods = periodic
    my_coords = (/ 0, 0, 0 /)
    myid = 0
#endif

    do i=1,3
       select case(i)
       case (1)
          nip => nix
          ngp => ngx
       case (2)
          nip => niy
          ngp => ngy
       case (3)
          nip => niz
          ngp => ngz
       end select
       if (periodic(i)) then
          nip = ngp
       else
          nip = ngp -2
       end if
    end do
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


    gstart(:) = [ sx, sy, sz ]
    gend(:)   = [ ex, ey, ez ]

    ! might need improvement
    ! box length might be needed
    do i=1, 3
       if (periodic(i)) then
          dxyz(i) = 1.0_wp/ngxyz(i)
       else
          dxyz(i) = 1.0_wp/(ngxyz(i) -1)
       end if
    end do

    dx = dxyz(1)
    dy = dxyz(2)
    dz = dxyz(3)

  end subroutine set_grid_partition

end module grid_data
