! tests point by point the differences
! between 2 data files dumped with MPI IO
program check_q
 implicit none

 integer, parameter :: wp=kind(0.d0)

 ! data read from input file
 integer nx1,nz1,ny1,nx2,ny2,nz2
 character(len=100) :: fn1,fn2
 logical :: verbose=.false., print_all = .false.

 ! locals
 integer  i,j,k,n,ic,imax,jmax,kmax,nmax
 real(wp) x1,x2,maxdiff, avdiff

 open(22,file="check_q.in",status="old")
 read(22,'(a)') fn1
 read(22,'(a)') fn2
 read(22,*,end=101) verbose, print_all !if not present use default
 101 close(22)

 print*, trim(fn1)
 print*, trim(fn2)

! allocate(q1(nx,nz,ny,nvar),q2(nx,nz,ny,nvar))

 open(23,file=fn1,status="old",form="unformatted",access="stream")
 read(23) nx1, ny1, nz1
 open(24,file=fn2,status="old",form="unformatted",access="stream")
 read(24) nx2, ny2, nz2

  print*,  trim(fn1), " nx,ny,nz",nx1,ny1,nz1
  print*,  trim(fn2), " nx,ny,nz",nx2,ny2,nz2

  if ( nx1 /= nx2 .or. &
       ny1 /= ny2 .or. &
       nz1 /= nz2 ) then
     write(0,*) 'problem with the input file grid sizes'
     stop
  endif

 ic=0
 maxdiff=0.0_wp
 avdiff=0.0_wp
 do n=1,1
  do k=1,nz1
   do j=1,ny1
    do i=1,nx1
     read(23)x1
     read(24)x2
     if ( verbose) then
       print*, i,j,k,x1,x2
     endif
     if(x1 /= x2) then
      if (abs(x1-x2) > maxdiff) then
       maxdiff=abs(x1-x2)
       imax=i; jmax=j; kmax=k; nmax=n
      endif
      if (print_all .and. .not. verbose) then
         print*, i,j,k,x1,x2
      end if
      avdiff =avdiff + abs(x1-x2)
      ic=ic+1
   endif
    enddo
   enddo
  enddo
 enddo

 close(23)
 close(24)

  if(ic > 0) then
     print*, 'array size ',nx1*ny1*nz1, ' of which' , ic, ' points differ '
     print*, 'average max diff ', avdiff/ic
     print*,'largest difference ', maxdiff, 'at ',imax,jmax,kmax
  endif

end program check_q
