! driver for testing DL_MG
! initiailises MPI, reads input and calls the tests
program main
!  use dl_mg_mpi
  use dl_mg_params, only : wp, pi
  use mg_tests, only : pbe => pbe_params, periodic, poisson_test
  use mg_utils
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif

  integer provided, mg_comm, pdims(3),ngx,ngy,ngz,nlx,nly,nlz, fd_maxiters, &
       niter, local_shift(3),ierr
  integer slice_index, slice_section, fd_order
  logical pperiods(3), use_pbelinear, use_steric, use_fas, &
       use_damping, errors_return, use_cg
  character(len=64) test_name, model_name, neutr_method
  real(wp) tol_res_rel, tol_res_abs, tol_pot_rel, tol_pot_abs, &
       tol_newton_rel, tol_cg_rel, tol_mg_rel, temp
  real(wp), allocatable :: x_ratios(:)

#ifdef MPI
  call MPI_init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  call MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN, ierr);
  if ( provided /= MPI_THREAD_FUNNELED) then
     write(0,*) "warning, required MPI thread safety level not provided"
  end if
#endif

  call read_input

#ifdef MPI
  !write(*,*) 'WARNING: in this verison the topology is 1D, nzp=nproc'
  !write(*,*) '         the pdims line from input file is neglected  '
  !call mpi_comm_size(mpi_comm_world,pdims(3),ierr)
  !pdims(1)=1; pdims(2)=1;
  call mpi_cart_create(mpi_comm_world,3,pdims,pperiods,.false.,mg_comm,ierr)
#else
  periodic = pperiods
  write(*,*) 'periodic', periodic
#endif

  call grid_partition(nlx,nly,nlz,local_shift)

  select case(test_name)
  case("poisson")
     call poisson_test(mg_comm,ngx,ngy,ngz,local_shift,nlx,nly,nlz, use_cg, &
          fd_maxiters, niter,model_name,&
          tol_res_rel, tol_res_abs, tol_pot_rel, tol_pot_abs, &
          tol_newton_rel, tol_cg_rel, tol_mg_rel, &
          fd_order, use_damping, errors_return, &
          neutr_method, x_ratios)
!  case("laplace")
!     call laplace_test(mg_comm,ngx,ngy,ngz,local_shift,nlx,nly,nlz,niter,tol)
  case default
     write(0,*)'unknown test name in main, quitting ....'
#ifdef MPI
     call mpi_abort(mpi_comm_world,1,ierr)
#else
     STOP
#endif
  end select

#ifdef MPI
  call MPI_Finalize(ierr)
#endif

  contains

    subroutine read_input
      implicit none

      integer i, myid, nproc, ierr
      character(len=64),allocatable :: file(:)
      character(len=64) :: line
      real(wp) aux
      logical ltmp

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
101      continue

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
102      close(44)

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
            read(file(1),*) ngx,ngy,ngz
            !write(*,*) ' using grid sizes ', myid,' :',ngx,ngy,ngz
         case(2)
            read(file(2),*) pdims
            !write(*,*) ' using processor grid ', myid,' :',pdims
         case(3)
            read(file(3),*) pperiods
            !write(*,*) ' using processor grid periodicity ', myid,' :',pperiods
         case(4)
            read(file(4),*) use_cg
         case (5)
            read(file(5),*) fd_maxiters, niter
         case (6)
            ! write slice if index > 0
            read(file(6),*) slice_index, slice_section
         case (7)
            read(file(7),*) test_name, model_name
            !write(*,*) ' test ', trim(test_name), 'model name ', trim(model_name)
         case(8)
            ! convergence tolerance
            read(file(8),*) tol_res_rel, tol_res_abs, tol_pot_rel, tol_pot_abs, &
                 tol_newton_rel, tol_cg_rel, tol_mg_rel
         case(9)
            ! finite diference order
            read(file(9),*) fd_order, use_damping
         case (10)
            ! PBE data : temperature, number of ion species, fixed_nions
            read(file(10), *) pbe%temp, pbe%n
         case (11)
            ! PBE concentrations
            allocate(pbe%c(pbe%n), pbe%q(pbe%n))
            read(file(11), *) pbe%c(:)
            ! convert from molar to rB**-3
            pbe%c(:) = pbe%c(:)* 8.9236832078390854265458e-5_wp
         case (12)
            ! PBE ion charges
            read(file(12),*) pbe%q(:)
         case (13)
            ! PBE lambda prefactor
            read(file(13),*) ltmp, aux
            if (ltmp) then
               pbe%lam = aux
            else
               pbe%lam = -4.0_wp * pi
            end if
         case (14)
            read(file(14),*) pbe%linearised, pbe%use_steric, errors_return, pbe%use_fas
         case (15)
            read(file(15), *) neutr_method
         case (16)
            allocate(x_ratios(pbe%n))
            read(file(16), *) x_ratios
         end select
      end do

#ifdef MPI
      call mpi_comm_size(mpi_comm_world, nproc, ierr)
      if ( nproc /= pdims(1) * pdims(2) * pdims(3) ) then
         write(0,*) "nproc /= npx * npy * npz !, Aborting ..."
         call mpi_abort(mpi_comm_world,1,ierr)
      end if
#endif

    end subroutine read_input


    subroutine grid_partition(nlx,nly,nlz,local_shift)
      implicit none

      integer, intent(out) :: nlx,nly,nlz, local_shift(3)

      integer pgrid_dims(3),my_coords(3), ngxyz(3),myid,nproc,ierr
      integer r(3), nl(3), i, nlmax(3),nlmin(3)
      logical pgrid_periods(3)
      logical :: use_input_partition
      integer, allocatable :: shifts_in(:),local_sizes(:)

#ifdef MPI
      call mpi_cart_get(mg_comm, 3, pgrid_dims, pgrid_periods, my_coords, ierr)
      call mpi_comm_rank(mg_comm,myid,ierr)
#else
      pgrid_dims = (/ 1, 1, 1/)
      pgrid_periods = (/ .false.,  .false.,  .false. /)
      my_coords = (/ 0, 0, 0 /)
      myid = 0
#endif

      call get_environment_variable("USE_INPUT_PARTITION",length=i)
      if ( i> 0) then
        use_input_partition = .true.
     else
        use_input_partition = .false.
     endif


      if (use_input_partition) then
         if (myid == 0) then
            write(*, *) 'using input rank partition'
            ! only for z direction
            allocate(shifts_in(0:pgrid_dims(3)-1),local_sizes(0:pgrid_dims(3)-1))
            open(397,file="input_mpi_partition_128",status="old")
            read(397,*) shifts_in
            close(397)
            ! we need the shifts
            local_sizes(0) = shifts_in(0)
            do i=1,pgrid_dims(3)-1
               local_sizes(i)= shifts_in(i)-shifts_in(i-1)
            enddo
            local_sizes(pgrid_dims(3)-1)=local_sizes(pgrid_dims(3)-1)+1 !boundary

            do i=pgrid_dims(3)-1,1,-1
               shifts_in(i)=shifts_in(i-1)

            enddo
            shifts_in(0) = 0

            do i=0,pgrid_dims(3)-1
               write(0,*) 'shift local size', shifts_in(i),local_sizes(i)
            enddo
         else
            allocate(shifts_in(1),local_sizes(1))
         endif

#ifdef MPI
         call mpi_scatter(shifts_in, 1, mpi_integer, &
              local_shift(3),1, mpi_integer, &
              0,mg_comm,ierr)
         call mpi_scatter(local_sizes, 1, mpi_integer, &
              nl(3),1, mpi_integer, &
              0,mg_comm,ierr)
#endif

         nl(1) = ngx/pgrid_dims(1)
         nl(2) = ngy/pgrid_dims(2)
         local_shift(1:2) = 0
      else
         nl(1) = ngx/pgrid_dims(1)
         nl(2) = ngy/pgrid_dims(2)
         nl(3) = ngz/pgrid_dims(3)

      ! test for silly values ?

         ngxyz = (/ ngx, ngy, ngz /)

         do i=1,3

            r(i) = mod(ngxyz(i),pgrid_dims(i))

            if ( my_coords(i) < r(i) ) then
               nl(i) = nl(i)+1
               local_shift(i) = my_coords(i)*nl(i)
            else
               local_shift(i) = r(i) * (nl(i)+1) + (my_coords(i) - r(i)) * nl(i)
            endif

         end do
   endif

      !write(0,*) 'local starts', myid, local_shift
      nlx = nl(1); nly = nl(2); nlz = nl(3)

#ifdef MPI
      ! some info on grid partition
      call mpi_reduce(nl,nlmax,3,mpi_integer,mpi_max,0,mg_comm,ierr)
      call mpi_reduce(nl,nlmin,3,mpi_integer,mpi_min,0,mg_comm,ierr)

      if (myid == 0) then
         write(*,*)
         write(*,*) 'min local domain sizes ', nlmin
         write(*,*) 'max local domain sizes ', nlmax
         write(*,*)
      endif
#endif

    end subroutine grid_partition


end program main
