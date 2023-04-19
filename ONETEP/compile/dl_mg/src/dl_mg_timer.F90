!> \brief internal timer
!!
!!  measures the computational and communication total times and allso
!! at each multigrid level node 0 the writes average time, min and max
!! ( accros active ranks for each level)
!!
!! If the environment variable DL_MG_TIMER_MG_LEVELS is defined times for multigrid components (per thread) are collected
!!
!! Lucian Anton
!!
module dl_mg_timer
  implicit none
  public

  integer, parameter :: INIT=-1, START=1, STOP=0, TCOMPUTE=1, TCOMM=2, &
       TSOLVER = -2, TDEFCO = -3, TDEFCO_MPI=-4, TMG=-5, &
       TRELAX = 1, TRESTRICT = 2, TPROLONG = 3, &  ! multgrid components
       TCONV = 4, TIGNORE = -111, max_components = 4
  logical, save :: collect_mg_levels =.false.
#ifdef USE_TIMER
  real(kind(0.d0)), allocatable, private :: tdata(:,:,:,:), tstart(:,:, :)
  real(kind(0.d0)), private :: solver_time(TMG:TSOLVER), solver_start(TMG:TSOLVER) ! accumulates the total time spent in dl_mg_solver
#endif

contains

   subroutine mg_timer(switch, component, level, sector )
!$   use omp_lib
     use dl_mg_mpi_header
     use dl_mg_common_data, only : mg_levels, nthreads
     implicit none
     integer, intent (in) :: switch, & !< start/stop
          component, & !< collects time for solver component: tsolver, tdefco, tdefco_mpi,
                       !! tmg, trelax, trestrict, tprolong, ...
          level, &     !< multigrid level
          sector       !< coomunicationor computation, relevant only for multigrid components
#ifdef USE_TIMER

     integer thid, i
     real (kind(0.d0)) t2
#ifndef MPI
     real rt ! to use with cpu_time for pure serial code
#endif

     select case(switch)
     case(START)
        if (component < 0) then
#ifdef MPI
           solver_start(component) = mpi_wtime()
#else
!$         if (.true.) then
!$         solver_start(component) = omp_get_wtime()
!$         else
           call cpu_time(rt)
           solver_start(component) = rt
!$         endif
#endif

        else
           if(collect_mg_levels) then
              thid = 0
!$            thid = omp_get_thread_num()
!$            if (.true. ) then
!$              tstart(sector, component, thid) = omp_get_wtime()
!$            else
#ifdef MPI
                tstart(sector, component, thid) = mpi_wtime()
#else
                call cpu_time(rt)
                tstart(sector, component, thid) = rt
#endif
!$            endif
             endif
        endif
     case(STOP)
        if (component < 0) then
#ifdef MPI
           solver_time(component) = solver_time(component) + &
                mpi_wtime() - solver_start(component)
#else
!$         if (.true.) then
!$         solver_time(component) = solver_time(component) + &
!$             omp_get_wtime() - solver_start(component)
!$         else
           call cpu_time(rt)
           solver_time(component) = solver_time(component) + &
                rt - solver_start(component)
!$         endif
#endif
        else
           if (collect_mg_levels) then
              thid = 0
!$      thid = omp_get_thread_num()
!$      if (.true. ) then
!$         t2 = omp_get_wtime()
!$      else
#ifdef MPI
              t2 = mpi_wtime()
#else
              call cpu_time(rt)
              t2 = rt
#endif
!$      endif
              if (tstart(sector, component, thid) < 0.d0 ) then
                 write(0,*) "dl_mg_timer: internal storage t1 < 0,"
                 write(0,*) "possible inconsistent order call of the timer."
                 write(0,*) "Component", component, " sector ", sector
              else
                 tdata(sector, component, thid, level) = tdata(sector, component, thid, level) &
                      + (t2-tstart(sector, component, thid))
                 tstart(sector, component, thid)=-666.d0
              endif
           endif
        endif
    CASE(INIT)
       ! tdata layout
       ! dim1 =>sector ( computation or communication )
       ! dim2 => components: relax, restrict, prolongate, convergence
       ! dim3 => OpenMP thread
       ! dim4 => mg levels
       call get_environment_variable("DL_MG_TIMER_MG_LEVELS",length=i)
       if ( i > 0)  collect_mg_levels =.true.
       if (collect_mg_levels) then
          allocate(tdata( 2, max_components, 0:nthreads-1, mg_levels), tstart(2, max_components, 0:nthreads-1))
          tdata(:,:,:,:) = 0.d0
       end if
       solver_time(:) = 0.d0
    end select
#endif
   end subroutine mg_timer


   subroutine report_timings
     use dl_mg_common_data, only : unit => report_unit, mg_levels, mg_comm, nthreads
     use dl_mg_grids, only : mg
     use dl_mg_mpi_header
     implicit none
#ifdef USE_TIMER
     integer n, i, j, k, t, mg_id, level_id, ndata, kmax, ierr
#ifdef MPI
     integer leveltag(mg_levels), rq(mg_levels)
     integer waitall_status(MPI_STATUS_SIZE, mg_levels)
#endif
     integer nproc_level(mg_levels)
     real(kind(0.d0)) tstat(2, max_components, 0:nthreads-1, mg_levels,3)
     real(kind(0.d0)), dimension(TMG:TSOLVER) :: st_sum, st_min, st_max
     real(kind(0.d0)), allocatable :: bbuff(:,:), buff(:)
     character(len=100) str


#ifdef MPI

     call mpi_comm_rank(mg_comm, mg_id, ierr)
     call mpi_comm_size(mg(mg_levels)%comm, nproc_level(mg_levels), ierr)
     n = size(solver_time)
     call mpi_reduce(solver_time, st_sum, n, MPI_DOUBLE_PRECISION, &
          MPI_SUM, 0, mg_comm,ierr)
     call mpi_reduce(solver_time, st_min, n, MPI_DOUBLE_PRECISION, &
          MPI_MIN, 0, mg_comm,ierr)
     call mpi_reduce(solver_time, st_max, n, MPI_DOUBLE_PRECISION, &
          MPI_MAX, 0, mg_comm,ierr)

     if (collect_mg_levels) then
        ndata = size(tstat(:,:,:,1:1,:)) ! amount of stat data per mg level
        allocate(buff(ndata+1))  ! needed for normalisation
        kmax = 3 ! print average, min max
        ! stat for total solver time
        leveltag(:) = 0 ! mark the ranks that collect data at each level

        ! accumulate stats at each level
        do t = 1, mg_levels
           if ( .not. mg(t)%active ) cycle
           call mpi_comm_rank(mg(t)%comm, level_id, ierr)
           if ( level_id == 0 ) then
              leveltag(t) = 1
              ! mg level size it obtained above the if block
              if (t < mg_levels) &
                   call mpi_comm_size(mg(t)%comm, nproc_level(t), ierr)
           endif

           n = size(tdata(:,:,:,t))
           call mpi_reduce(tdata(:,:,:,t), tstat(:,:,:,t,1), n, MPI_DOUBLE_PRECISION, &
                MPI_SUM, 0, mg(t)%comm,ierr)
           call mpi_reduce(tdata(:,:,:,t), tstat(:,:,:,t,2), n, MPI_DOUBLE_PRECISION, &
                MPI_MIN, 0, mg(t)%comm,ierr)
           call mpi_reduce(tdata(:,:,:,t), tstat(:,:,:,t,3), n, MPI_DOUBLE_PRECISION, &
                MPI_MAX, 0, mg(t)%comm,ierr)
        enddo

        do t = 1, mg_levels
           ! this rank has stat data but is not rank 0 in mg_comm, hence it must send its data
           if ( leveltag(t) > 0 .and. mg_id > 0) then
              buff(1:ndata) = reshape(tstat(:,:,:,t,:), (/ ndata /))
              buff(ndata+1) = nproc_level(t)
              call mpi_send(buff, ndata+1 , mpi_double_precision, 0, 287+t, mg_comm, ierr)
           end if
        enddo
     endif
#else
     mg_id = 0
     nproc_level(:) = 1
     kmax = 1 ! no need of min max in time stats
#endif

     if ( mg_id == 0) then
#ifdef MPI
        if (collect_mg_levels) then
           allocate(bbuff(ndata+1,mg_levels))
           rq(:) = MPI_REQUEST_NULL
           do t = 1, mg_levels
              ! mg_comm rank 0 does not have the data for level t; someone must send them
              if ( leveltag(t) == 0) then
                 call mpi_irecv(bbuff(:,t), ndata+1, mpi_double_precision, MPI_ANY_SOURCE, 287+t, mg_comm, &
                      rq(t), ierr)
              endif
           enddo

           call mpi_waitall(mg_levels, rq, waitall_status, ierr)

           do t=1, mg_levels
              if ( leveltag(t) == 0) then
                 tstat(:,:,:,t,:) = reshape(bbuff(1:ndata,t), (/2,max_components,nthreads,3/))
                 nproc_level(t) = bbuff(ndata+1,t)
              endif
           enddo
           deallocate(bbuff)
        endif
#else
        st_sum = solver_time
        st_max = -1.0d0
        st_min = -1.0d0

        if (collect_mg_levels) then
           tstat(:,:,:,:,1) = tdata
           tstat(:,:,:,:,2:3) = -1.0d0
        endif
#endif


        ! report data:

        ! A. top level times
#ifdef MPI
        write(unit, '(19x,3(a10))') ' average ', '   min   ', '    max   '
        write(unit, '(a,3(E10.4,1x))') 'solver time:       ', &
             st_sum(TSOLVER)/nproc_level(mg_levels), st_min(TSOLVER), st_max(TSOLVER)
        write(unit, '(a,3(E10.4,1x))') 'derivatives total: ', &
             st_sum(TDEFCO)/nproc_level(mg_levels), st_min(TDEFCO), st_max(TDEFCO)
        write(unit, '(a,3(E10.4,1x))') 'derivatives  MPI:  ', &
             st_sum(TDEFCO_MPI)/nproc_level(mg_levels), st_min(TDEFCO_MPI), st_max(TDEFCO_MPI)
        write(unit, '(a,3(E10.4,1x))') 'multigrid:         ', &
             st_sum(TMG)/nproc_level(mg_levels), st_min(TMG), st_max(TMG)
        write(unit,*)
#else
        write(unit, '(19x,(a10))') '   time   '
        write(unit, '(a,3(E10.4,1x))') 'solver time:       ', st_sum(TSOLVER)
        write(unit, '(a,3(E10.4,1x))') 'derivatives total: ', st_sum(TDEFCO)
        write(unit, '(a,3(E10.4,1x))') 'derivatives  MPI:  ', st_sum(TDEFCO_MPI)
        write(unit, '(a,3(E10.4,1x))') 'multigrid:         ', st_sum(TMG)
        write(unit,*)
#endif
        ! B. total time (across mg levels) for the main components
        if (collect_mg_levels) then
           ! average the time
        do t = 1, mg_levels
           tstat(:,:,:,t,1) = tstat(:,:,:,t,1)/nproc_level(t)
        enddo
           write(unit, '(a)') 'total time spent over relax, restrict, prolungate components: '
           write(unit, '(a)') 'odd columns are totals, even are the communication: '
           do k = 1, kmax ! average, min`, max
              write(unit, '(3x,1000(E10.4,1x))') ((sum(tstat(j,1:3,i,:,k)), j=1,2),i=0,nthreads-1)
           end do
           write(unit,'(1x/1x)')

           do i = 1, max_components ! over mg components: relax, restric, prolongate and convergence
              select case(i)
              case (TRELAX)
                 write(unit,'(a)') "component: relax"
              case(TRESTRICT)
                 write(unit,'(a)') "component: restrict"
              case(TPROLONG)
                 write(unit,'(a)') "component: prolongation"
              case(TCONV)
                 write(unit,'(a)') "component: convergence test"
              case default
                 write(unit,'(a, I3)') "component:", i
              end select
              ! total over level first
              write(unit,'(a)')    '=========================================='
              do k = 1, kmax
                 write(unit, '(3x,1000(E10.4,1x))') ((sum(tstat(j,i,t,:,k)), j=1,2), t=0,nthreads-1)
              enddo
              write(unit,'(a)')    '-----------------------------------------'
              do t=mg_levels, 1, -1
                 write(unit,'(a,I4)') "level", t
                 do k = 1, kmax ! average, min, max
                    write(unit, '(3x,1000(E10.4,1x))') ((tstat(j,i,n,t,k), j=1,2), n=0,nthreads-1)
                 enddo
                 write(unit,*)
              enddo
              write(unit,*)
           enddo
        endif
     end if

#endif
   end subroutine report_timings


   subroutine timer_free
     implicit none
#ifdef USE_TIMER
     if (allocated(tdata)) deallocate(tdata)
     if (allocated(tstart)) deallocate(tstart)
#endif
   end subroutine timer_free

end module dl_mg_timer
