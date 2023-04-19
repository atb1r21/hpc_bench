  if (present(scale)) then
     scl = scale
  else
     scl = 1.0_wp
  end if


  in_parallel = .false.
  !$ in_parallel = omp_in_parallel()


  !$OMP MASTER
  xloc = 0.0_wp
  !$OMP END MASTER
  !$OMP BARRIER

  if (in_parallel)then
     call loops_ab(a,b,xloc)
  else
     !$omp parallel
     call loops_ab(a,b,xloc)
     !$omp end parallel
  end if

  !$OMP MASTER
#ifdef MPI
  call mpi_allreduce(xloc, xworld, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
#else
  xworld = xloc
#endif
  !$OMP END MASTER
  !$OMP BARRIER