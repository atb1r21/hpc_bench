program test_onetep_ls
  use mpi
  use dl_mg
  implicit none

  real(kind(0.d0)), parameter :: pi = 4.0d0 * atan(1.0d0)
  integer pq1f, pq2f, pq3f, num_slabs12_pq, nproc, myid, mg_comm, ierr
  integer idx_s(3), idx_e(3), nfiles
  integer, allocatable :: zsizes(:)
  real(kind(0.d0)) d1f, d2f, d3f
  real(kind(0.d0)), allocatable :: pot(:,:,:), rho(:,:,:), eps_full(:,:,:), res(:,:,:), eps_half(:,:,:,:)
  character(len=128) model_name
  character(len=10) timetag
  real(kind(0.d0)) t1, t2, tsolver_min, tsolver_max, tsolver_sum, pot_norm
  integer agg_level

  call mpi_init(ierr)

  !model_name="onetep_t27"
  !model_name = "onetep_t28"
  !model_name="data_onetep_ls/vacuum"
  !model_name="data_onetep_ls/solvent"

  call mpi_comm_size(mpi_comm_world, nproc, ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)

  call get_command_argument(1,model_name)
  if ( myid == 0) then
     write(*,*)
     write(*,*) 'using model ', model_name
     write(*,*)
  endif

  call get_command_argument(2,timetag)
  if ( len_trim(timetag) == 0) then
     agg_level = 0
     if ( myid == 0) then
        write(*,*) "No full aggregation level provided!"
        write(*,*) "Using default value 0."
     endif
  else
     read(timetag,*) agg_level
     if ( myid == 0) then
        write(*,*) "Aggregation level :", agg_level
        write(*,*)
     endif
  endif

  select case(model_name)
  case("onetep_t27") ! QC test 27
     pq1f = 161
     pq2f = 161
     pq3f = 161
     d1f  = 0.285714285714286
     d2f  = 0.285714285714286
     d3f  = 0.285714285714286
     nfiles = 2
  case("onetep_t28")
     pq1f = 193
     pq2f = 193
     pq3f = 193
     d1f  = 0.221938775510204
     d2f  = 0.221938775510204
     d3f  = 0.221938775510204
     nfiles = 3
  case("data_onetep_ls/vacuum")
    pq1f = 577
    pq2f = 577
    pq3f = 673
    d1f  = 0.25d0
    d2f  = 0.25d0
    d3f  = 0.25d0
    nfiles = 2
  case("data_onetep_ls/solvent")
     pq1f = 641
     pq2f = 641
     pq3f = 737
     d1f  = 0.225155279503106
     d2f  = 0.225155279503106
     d3f  = 0.224734042553191
     nfiles = 3
   case("lysozyme")
      pq1f = 449
      pq2f = 545
      pq3f = 609
      d1f  = 0.225155279503106
      d2f  = 0.225155279503106
      d3f  = 0.224734042553191
      nfiles = 2
   end select

  call mpi_cart_create(mpi_comm_world, 3, (/ 1, 1, nproc /), &
       (/ .false., .false., .false. /), .false., mg_comm,ierr)

  num_slabs12_pq = (pq3f - 2)/ nproc
  if ( myid < mod(pq3f-2,nproc)) num_slabs12_pq = num_slabs12_pq + 1

  if ( myid == 0 )        num_slabs12_pq =  num_slabs12_pq + 1
  if ( myid == nproc - 1) num_slabs12_pq =  num_slabs12_pq + 1

  allocate (zsizes(0:nproc-1))
  call mpi_allgather(num_slabs12_pq,1,mpi_integer,zsizes,1,mpi_integer,mg_comm,ierr)

  idx_s = (/ 1,       1, sum(zsizes(0:myid - 1)) + 1 /)
  idx_e   = (/ pq1f, pq2f, idx_s(3) + num_slabs12_pq - 1 /)

  call date_and_time(time=timetag)
  call dl_mg_init(pq1f, pq2f, pq3f, &
       d1f, d2f, d3f, &
       (/ DL_MG_BC_DIRICHLET, DL_MG_BC_DIRICHLET, DL_MG_BC_DIRICHLET /), &
       idx_s, idx_e, &
       150, 98, "dl_mg_report_"//timetag//".txt", mg_comm,&
       full_aggregation_level=agg_level)
       !40, 98, "dl_mg_report_"//timetag//".txt", mg_comm,&
       !full_aggregation_level=agg_level)

  allocate(pot(pq1f, pq2f, num_slabs12_pq), rho(pq1f, pq2f, num_slabs12_pq), &
       res(pq1f, pq2f, num_slabs12_pq))

  allocate(eps_half(pq1f, pq2f, num_slabs12_pq,3))
  if ( nfiles == 2) then
     eps_half(:,:,:,:) = 1.d0
  endif

  call read_onetep_data

  t1 = mpi_wtime()
  !call dl_mg_solver_poisson(pot, res, -4.0d0*pi, rho, eps_half=eps_half, &
  !        tol=1.d-5)
  ! JCW: Update to function interface
  !call dl_mg_solver_poisson(eps_half, -4.0d0*pi, rho, pot, tol=1.d-5, res=res)
  t2 = mpi_wtime()

  call mpi_reduce(sum(pot(1:pq1f,1:pq2f,1:num_slabs12_pq)**2),pot_norm,1,mpi_double_precision,MPI_SUM,0,mg_comm,ierr)

  call mpi_reduce(t2-t1,tsolver_max, 1, mpi_double_precision, MPI_MAX, 0, mg_comm,ierr)
  call mpi_reduce(t2-t1,tsolver_min, 1, mpi_double_precision, MPI_MIN, 0, mg_comm, ierr)
  call mpi_reduce(t2-t1,tsolver_sum, 1, mpi_double_precision, MPI_SUM, 0, mg_comm, ierr)


  if (myid == 0)  then
     write(*,*) 'done, norm pot ', sqrt(pot_norm)
     write(*,*) 'solver times:    max         min        ave'
     write(*,'(I4,11x,3(E10.4,1x))') nproc, tsolver_max, tsolver_min, tsolver_sum/nproc
  end if

   call mpi_finalize(ierr)

  contains

    subroutine read_onetep_data
      !use test_onetep_utils
      implicit none
      integer i, j, fh, nt, gsizes(3), lstart(3), lsizes(3), ierr
      integer ftype, ftype_rho, ftype_eps
      integer status(MPI_STATUS_SIZE), nr

      integer(kind=MPI_OFFSET_KIND) disp
      character(len=128) fname
      real(kind(0.d0)) xmin, xmax

      lsizes = (/ pq1f, pq2f, num_slabs12_pq /)

      nt = pq1f * pq2f * num_slabs12_pq

      gsizes =(/ pq1f, pq2f, pq3f /)

      !write(0,*) 'zsizes ', zsizes, myid

      lstart(1) = 0
      lstart(2) = 0
      lstart(3) = sum(zsizes(0:myid-1))

      call mpi_type_create_subarray(3,gsizes,lsizes,lstart,&
           MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,ftype_rho,ierr)
      call mpi_type_commit(ftype_rho,ierr)
      if ( nfiles == 3) then
         call mpi_type_create_subarray(4, (/ pq1f, pq2f, pq3f, 3/), &
              (/ pq1f, pq2f, num_slabs12_pq, 3 /), &
              (/ 0   , 0   , lstart(3)     , 0 /), &
              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,ftype_eps,ierr)
         call mpi_type_commit(ftype_eps,ierr)
      endif
!!$        if (myid == 0) then
!!$          print*, 'array sizes '
!!$          print*, 'gzises', gsizes
!!$          print*, ' x, y ', size(pot,1), size(pot,2)
!!$          print*, 'z ', zsize
!!$        endif

          do i = 1, nfiles ! 3 if eps_half is not constant

             select case(i)
             case (1)
                fname=trim(model_name)//'/pot.dat'
                ftype = ftype_rho
             case(2)
                fname=trim(model_name)//'/rho.dat'
                ftype =ftype_rho
             case(3)
                fname=trim(model_name)//'/eps_half.dat'
                ftype =ftype_eps
             end select

             call mpi_file_open(mg_comm,fname,&
                  MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
             if (ierr /= 0) then
                write(0,*) 'MPI IO error in opening file on rank ', myid
             end if
             disp = 0
             call mpi_file_set_view(fh,disp,MPI_DOUBLE_PRECISION,ftype,&
                  "native",MPI_INFO_NULL,IERR)

             select case (i)
             case (1)
                call mpi_file_read_all(fh,pot,  &
                     nt, MPI_DOUBLE_PRECISION, status,ierr)
                if ( ierr /= 0 ) then
                   write(0,*) 'error in mpi read pot'
                endif
             call mpi_get_count(status, mpi_double_precision, nr, ierr)
             write(0,*) 'read pot', nr
             !call minmax_grid(pot,mg_comm,xmin,xmax)
             !write(0,*) 'pot min max ', xmin, xmax
             case (2)
                call mpi_file_read_all(fh,rho,  &
                     nt,MPI_DOUBLE_PRECISION, status,ierr)
             call mpi_get_count(status, mpi_double_precision, nr, ierr)
             write(0,*) 'read rho', nr
             !call minmax_grid(rho,mg_comm,xmin,xmax)
             !write(0,*) 'rho min max ', xmin, xmax
          case(3)
             call mpi_file_read_all(fh, eps_half, &
                  3*nt, MPI_DOUBLE_PRECISION, status,ierr)
             call mpi_get_count(status, mpi_double_precision, nr, ierr)
             write(0,*) 'read eps_half', nr
             !do j=1, 3
             !   call minmax_grid(eps_half(:,:,:,j),mg_comm,xmin,xmax)
             !   write(0,*) 'eps min max ', j, xmin, xmax
             !enddo
          end select

          call mpi_file_close(fh,ierr)
          call mpi_barrier(mg_comm,ierr)

       end do

       call mpi_type_free(ftype,ierr)

     end subroutine read_onetep_data

end program test_onetep_ls
