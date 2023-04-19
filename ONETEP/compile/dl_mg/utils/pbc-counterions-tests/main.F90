program pbc_counterions_test
  use params
  implicit none

  call initialize

  call setup_rho_steric

  ! explore the ion mixture space
  call pbe_driver

  call finish_run

contains

  subroutine initialize
#ifdef MPI
    use mpi
#endif
    use model_data, only : init_arrays
    implicit none

    integer provided, ierr

#ifdef MPI
    call MPI_init_thread(MPI_THREAD_FUNNELED,provided,ierr)
    call MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN, ierr);
    if ( provided /= MPI_THREAD_FUNNELED) then
       write(0,*) "warning, required MPI thread safety level not provided"
    end if
#endif

    call read_input

    call init_arrays

  end subroutine initialize


  subroutine setup_rho_steric
    use params
    use model_functions, only : gauss, eps, compute_steric_weight
    use model_data
    implicit none

    integer i,j,k
    real(wp) x,y,z

    pot = 0.0_wp ! good for homogenous Dirichlet or PBC

    ! epsilon, rho and steric weight
    do k = sz, ez
       do j = sy, ey
          do i = sx, ex
             x=(i-1.0_wp)*dx
             y=(j-1.0_wp)*dy
             z=(k-1.0_wp)*dz

             eps_half(i,j,k,1) = eps(x+0.5_wp*dx,y,z)
             eps_half(i,j,k,2) = eps(x,y+0.5_wp*dy,z)
             eps_half(i,j,k,3) = eps(x,y,z+0.5_wp*dz)
             eps_full(i,j,k)   = eps(x,y,z)

             rho(i,j,k) = gauss(x,y,z)

             call compute_steric_weight(i,j,k,x,y,z, stw)
          enddo
       enddo
    enddo

  end subroutine setup_rho_steric


  subroutine pbe_driver
    use dl_mg
    use params
    implicit none

    real(wp) q_cnt
    integer boundary_condition(3), ierr

    where(periodic)
       boundary_condition = DL_MG_BC_PERIODIC
    elsewhere
       boundary_condition = DL_MG_BC_DIRICHLET
    end where

    call dl_mg_init(ngx, ngy, ngz, &
         dx, dy, dz, &
         boundary_condition, &
         [sx, sy, sz], [ex, ey, ez], &
         mg_comm, &
         98, 'log',&
         full_aggregation_size=0, full_aggregation_level=0, &
         errors_return = .false., ierror=ierr)



    q_cnt = q_counterions
    do while (q_cnt <= 1.0_wp .and. q_cnt >= 0.0_wp &
         .or. .not. use_counterions)
       ! if using jellium do one computation irrespective of q_cnt value
       call solve_pbe(q_cnt)

       call compute_energy(q_cnt)

       if (delta_q_cntions == 0 .or. .not. use_counterions) exit
       q_cnt = q_cnt + delta_q_cntions

       ! are the new concentrations positive ?
       if ( .not. positive_concentrations(q_cnt))then
          exit
       end if
    end do

  end subroutine pbe_driver


  logical function positive_concentrations(y)
    use model_data, only : rho, stw
    use model_functions, only : grid_sum
    implicit none

    real(wp), intent(in) :: y

    real(wp) rho_sum, stw_sum, c1, c2

    positive_concentrations = .true.

    call grid_sum(rho, rho_sum)
    if (pbe_params%use_steric) then
       call grid_sum(stw, stw_sum)
    else
       stw_sum = real(ngx,wp)*ngy*ngz
    end if
    if (use_counterions) then
       c1 = pbe_params%c(1) - y * rho_sum/(pbe_params%q(1) * stw_sum)
       c2 = pbe_params%c(2) - (1.0_wp -y) * rho_sum /(pbe_params%q(2) * stw_sum)
       if ( c1 <= 0.0_wp .or. c2 <= 0.0_wp) then
          positive_concentrations =.false.
       end if
    end if

  end function positive_concentrations


  subroutine finish_run
#ifdef MPI
    use mpi
#endif
    use params, only : data_unit
    implicit none
    integer myid, ierr

#ifdef MPI
    call mpi_comm_rank(mpi_comm_world,myid,ierr)
#else
    myid  = 0
#endif

    if (myid == 0) then
       close(data_unit)
    end if

    call mpi_finalize(ierr)

  end subroutine finish_run


  subroutine read_input
#ifdef MPI
    use mpi
#endif
    use params, only : nproc, mg_comm, periodic, &
         ngx, ngy, ngz, use_pbe, delta_q_cntions, &
         periodic, xl,yl, zl, dump_fields
    use model_functions, only : pbe => pbe_params, &
         gauss => gauss_params, eps => eps_params, abort
    implicit none
    integer i, myid, ierr, pdims(3)
    character(len=64),allocatable :: file(:)
    character(len=64) :: line
    real(wp) dxyz(3)
    logical pperiods(3)

#ifdef MPI
    call mpi_comm_rank(mpi_comm_world,myid,ierr)
#else
    myid  = 0
    nproc = 1
#endif

    if (myid == 0) then

       !find how many lines of data are in input
       open(44,file="input",status="old")

       i=0
       do
          read(44,*,end=101)line
          line = adjustl(line)
          ! skip lines starting with # or empty
          if ( .not. (line(1:1) == '#' .or. len_trim(line) == 0) ) i = i + 1
       end do
101    continue

       allocate(file(i))

       ! read the file
       rewind(44)

       !print*, 'input file'

       i=0
       do
          read(44,'(a)',end=102)line
          line = adjustl(line)
          print*,  trim(line)
          if ( .not. (line(1:1) == '#' .or. len_trim(line) == 0) ) then
             i = i + 1
             file(i) = line
          endif

       end do
102    close(44)

#ifdef MPI
       call mpi_bcast(size(file),1,MPI_INTEGER,0,mpi_comm_world,ierr)
       call mpi_bcast(file,size(file)*len(file(1)),MPI_CHARACTER,0,mpi_comm_world,ierr)
#endif
    else
#ifdef MPI
       call mpi_bcast(i,1,MPI_INTEGER,0,mpi_comm_world,ierr)
       allocate (file(i))
       call mpi_bcast(file,size(file)*len(file(1)),MPI_CHARACTER,0,mpi_comm_world,ierr)
#endif
    end if

    !everybody reads the data

    do i=1, size(file)
       select case(i)
       case(1)
          read(file(i),*) ngx,ngy,ngz, xl, yl, zl
          !write(*,*) ' using grid sizes ', myid,' :',ngx,ngy,ngz
       case(2)
          read(file(i),*) pdims
          !write(*,*) ' using processor grid ', myid,' :',pdims
       case(3)
          read(file(i),*) pperiods
          !write(*,*) ' using processor grid periodicity ', myid,' :',pperiods
       case(4)
          ! PBE data : temperature, number of ion species, fixed_nions
          read(file(i), *) use_pbe, pbe%temp, pbe%n
       case(5)
          ! PBE concentrations
          if (use_pbe) then
             allocate(pbe%c(pbe%n), pbe%q(pbe%n))
             read(file(i), *) pbe%c(:)
             ! convert from molar to rB**-3
             pbe%c(:) = pbe%c(:)* 8.9236832078390854265458e-5_wp
          end if
       case (6)
          ! PBE ion charges
          if (use_pbe) then
             read(file(i),*) pbe%q(:)
          endif
       case (7)
          if (use_pbe) then
             read(file(i),*) pbe%linearised,  pbe%fixed_nions, &
                  pbe%use_steric, pbe%r_steric_weight
          endif
       case (8)
          if (use_pbe) then
             read(file(i), *) use_counterions, q_counterions, delta_q_cntions, use_smart_jellium
          endif
       case(9)
          ! finite diference order
          read(file(i),*) fd_order, use_damping
       case (10)
          read(file(i),*) fd_maxiters, niters
       case (11)
          ! convergence tolerance
          read(file(i),*) tol_res_rel, tol_res_abs, tol_pot_rel, tol_pot_abs, &
               tol_newton_rel, tol_mg_rel

       case (12)
          ! rho dist params
          read(file(i),*) gauss%sig, gauss%r0, gauss%scale
       case (13)
          ! eps params
          ! eps(r) =  1 + (eps0-1) * 0.5_wp * (1 + erf((r-d)/delta))
          read(file(i),*) eps%eps0, eps%d, eps%delta
       case (14)
          read(file(i),*) dump_fields
       end select

    end do

#ifdef MPI
    call mpi_comm_size(mpi_comm_world, nproc, ierr)
    if ( nproc /= pdims(1) * pdims(2) * pdims(3) ) then
       write(0,*) "nproc /= npx * npy * npz !, Aborting ..."
       call mpi_abort(mpi_comm_world,1,ierr)
    end if
#endif

    if (.not. all(pperiods) .and. pbe%fixed_nions) then
       call abort("fixed number of ions works with full periodic")
    end if

#ifdef MPI
  call mpi_cart_create(mpi_comm_world,3,pdims,pperiods,.false.,mg_comm,ierr)
#endif

  periodic = pperiods

  where(periodic)
     dxyz = [xl, yl, zl ]/ [ngx, ngy, ngz]
  elsewhere
     dxyz = [xl, yl, zl ]/ [ngx-1, ngy-1, ngz-1]
  end where

  dx = dxyz(1); dy = dxyz(2); dz = dxyz(3)


  end subroutine read_input


  subroutine solve_pbe(x)
    use params, pbe => pbe_params
    use model_data
    use dl_mg
    implicit none

    real(wp), intent(in) :: x

    integer np, ierror

    real(wp), allocatable, target :: dummy(:,:,:)
    real(wp), pointer :: stw_pt(:,:,:), c(:), q(:)

    np = pbe_params%n

    if (use_pbe) then
       if (.not. pbe%use_steric) then
          allocate(dummy(1:0,0:0,0:0)) ! zero size array
          stw_pt => dummy
       else
          stw_pt => stw
       end if
       c => pbe%c
       q => pbe%q


       call dl_mg_solver(eps_full, eps_half, -fourpi, &
            rho, &
            -fourpi, pbe%temp, np, c, q, pot, fd_order, &
            tol_res_rel = tol_res_rel, tol_res_abs = tol_res_abs, &
            tol_pot_rel = tol_pot_rel, tol_pot_abs = tol_pot_abs, &
            tol_newton_res_rel=tol_newton_rel, &
            tol_newton_res_abs=1.e-8_wp, &
            tol_vcycle_res_rel=tol_mg_rel, &
            tol_vcycle_res_abs=1.e-8_wp, &
            linearised=pbe_params%linearised, &
            steric_weight=stw_pt, use_pot_in = .false., &
            use_damping=pbe%use_damping, use_fas=.false., &
            conserve_ions_number = pbe%fixed_nions, &
            use_counterions = use_counterions, &
            q_counterions = x, &
            use_smart_jellium = use_smart_jellium, &
            betamu_elec = betamu, steric_weight_average = gamma, &
            ierror=ierror)

    else

       call dl_mg_solver(eps_full, eps_half, -fourpi,  rho,  pot, fd_order, &
            tol_res_rel = 1.e-8_wp, tol_res_abs= 1.e-8_wp, &
             use_damping=use_damping,&
            ierror=ierror)

    end if

  end subroutine solve_pbe


  subroutine compute_energy(y)
#ifdef MPI
    use mpi
#endif
    use params
    use model_data
    use model_functions, only : grid_sum, capexp
    implicit none

    real(wp), intent(in) :: y

    integer i,j,k, p, myid, ierr
    real(wp) en(3), gen(3), rmob, csum, beta, &
         c(pbe_params%n), c0(pbe_params%n), q(pbe_params%n), &
          sw, zeta, dv, rho_sum, stw_sum
    real(wp), allocatable :: rho_shift(:,:,:)

    logical, save :: first_call = .true.

    allocate(rho_shift(isx:iex,isy:iey,isz:iez))
    rho_shift = 0.0_wp


    c0 = pbe_params%c
    q = pbe_params%q

    beta = hartree / pbe_params%temp

    dv = dx*dy*dz

    call grid_sum(rho, rho_sum)
    if (pbe_params%use_steric) then
       call grid_sum(stw, stw_sum)
    else
       stw_sum = real(ngx,wp)*ngy*ngz
    end if

    !write(*,*) 'rho stw sum', rho_sum, stw_sum, beta

    if (use_counterions) then
       c(1) = c0(1) - y * rho_sum/(q(1) * stw_sum)
       c(2) = c0(2)- (1.0_wp -y) * rho_sum /(q(2) * stw_sum)
    else
       c = c0
       if (use_smart_jellium .and. use_pbe) then
          rho_shift(isx:iex,isy:iey,isz:iez) = stw(isx:iex,isy:iey,isz:iez) * rho_sum / stw_sum
       else
          rho_shift = rho_sum / (real(ngx,wp)*ngy*ngz)
       endif
    end if

    en = 0.0_wp
    do k = isz, iez
       do j = isy, iey
          do i = isx, iex
             if (pbe_params%use_steric) then
                sw = stw(i,j,k)
             else
                sw =1.0_wp
             end if
             rmob = 0.0
             csum = 0.0
             do p=1,pbe_params%n
                zeta = capexp( - beta * q(p) * pot(i,j,k) + betamu(p))
                rmob = rmob + c(p) * q(p) * sw * zeta
                csum = csum + c(p)  * sw * zeta
             end do
             en(1) = en(1) + (rho(i,j,k) - rho_shift(i,j,k)) * pot(i,j,k)
             en(2) = en(2) + rmob  * pot(i,j,k)
             en(3) = en(3) + csum
          enddo
       enddo
    enddo

#ifdef MPI
    call mpi_reduce(en, gen, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mg_comm, ierr)
#else
    gen = en
#endif

#ifdef MPI
    call mpi_comm_rank(mg_comm,myid,ierr)
#else
    myid  = 0
#endif
    if (myid == 0) then

       if (first_call) then
          first_call = .false.
          open(data_unit,file="GrandPot.csv", status="replace", position="append")
          ! write a header
          write(data_unit,'(11a19)') 'x', 'cb1 ', 'cb2 ', 'c1 ', 'c2 ', 'mu1 ', '  mu2     ', &
               '0.5\int rho*pot', '0.5\int c*pot', 'kBT \int sum c', 'kBT sum cb'
       endif

       write(data_unit,'(1x,E17.10,10(1x,",",E17.10))') y, pbe_params%c, c, betamu, 0.5_wp * gen(1)*dv, &
            0.5_wp * gen(2)*dv, 1.0_wp/beta * gen(3) * dv, 1.0_wp/beta * sum(c) * stw_sum * dv
       call flush(data_unit)
    end if
    call write_fields(beta, c)

  end subroutine compute_energy


  subroutine write_fields(beta, cb)
    use params, only : dump_fields, ngz, &
         sx, ex, sy, ey,  sz, ez, &
         use_pbe, &
         pbe=>pbe_params
    use model_data, only : stw, betamu, rho, pot
    use model_functions, only : capexp
    use mg_utils
    implicit none

    real(wp), intent(in) :: beta, cb(:)

    integer i,j,k,p
    character(len=30) fname
    real(wp), allocatable :: c(:,:,:)
    real(wp) sw

    if (dump_fields) then

       call mg_utils_init(mg_comm, ngx, ngy, ngz, [sx,sy,sz])

       call write_slice("xy",ngz/2, sx, sy, sz, pot, "pot")
       call write_slice("xy",ngz/2, sx, sy, sz, rho, "rho")
       if ( use_pbe) then
          allocate(c(sx:ex,sy:ey,sz:ez))
          do p=1,pbe%n
             do k=sz,ez
                do j=sy,ey
                   do i=sx,ex
                      if (pbe%use_steric) then
                         sw = stw(i,j,k)
                      else
                         sw =1.0_wp
                      end if
                      c(i,j,k) = sw * cb(p) &
                           * capexp(-beta * pbe%q(p) * pot(i,j,k) - betamu(p))
                   enddo
                enddo
             enddo
             write(fname,'(a1,i0)')"c",p
             call write_slice("xy",ngz/2, sx, sy, sz, c, fname)
          end do
       end if
    end if
  end subroutine write_fields

end program pbc_counterions_test
