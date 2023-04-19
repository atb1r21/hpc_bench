program grid_diagrams
  implicit none

  integer i, ii, j, jj, ip, jp, t, glevels, nx, ny, npx, npy, rx,ry
  integer, allocatable :: gpart_x(:), gpart_y(:)
  character(len=1) c, active, dead
  character(len=2) cc
  character(len=256) line

  active="X"
  dead="o"

  call read_input

  call grid_levels

  print*, 'glevels ', glevels

  open(34,file='gdia.out', status='unknown')

  do t = glevels, 1, -1

     c=active

     ip =-1

     ! ranks
     rx = 1

     do i = 1, nx

        ii = (i-1) / 2**(glevels-t)
        if( ii == ip ) then
           c=dead
        else
           c=active
           ip= ii
        endif

        write(line(1:1),'(a1)') c

        jp = 0
        ry = 1
        do j = 2, ny

           if ( line(1:1) == active) then
              jj = (j-1)/ 2**(glevels-t)
              if(jj == jp ) then
                 c=dead
              else
                 c=active
                 jp=jj
              endif
           endif

           if ( j-1 == gpart_y(ry)) then
              cc="| "
              ry = ry+1
           else
              cc="  "
           endif
           write(line(3*(j-1)-1:3*(j-1)+2), '(a2,a1)') cc,c
        enddo
        write(34,'(a)') line(1:3*(ny-1)+1)
        ! blank line or rank border
        if ( i == gpart_x(rx) ) then
           c="-"
           rx = rx+1
        else
           c=" "
        endif
        ry = 1
        do j = 1, 3*(ny-1)+1
           if ( j -1 == 3*gpart_y(ry) -2 ) then
              line(j:j)="|"
              ry=ry+1
           else
              line(j:j) = c
           endif
        enddo
        write(34,'(a)') line(1:3*(ny-1)+1)

     enddo

     write(34, *)
     write(34, *)

  enddo

  close(34)

contains

  subroutine grid_levels
    implicit none

    integer gridc(2), k(2), i

! find the number of levels for the global grid
        gridc=(/nx-1,ny-1/)
        do i=1,2
           k(i)=1
          do while (mod(gridc(i),2) == 0 .and. gridc(i) > 2)
            gridc(i)=gridc(i)/2
            k(i)=k(i)+1
          enddo
        enddo

        glevels = minval(k)

  end subroutine grid_levels


  subroutine read_input
    implicit none
    open(33,file="gdia.in",status="old")

    read(33,*)
    read(33,*) nx,ny
    read(33,*)
    read(33,*) npx, npy
    read(33,*)
    allocate( gpart_x(npx), gpart_y(npy))

    if ( npx > 1) then
       read(33,*) gpart_x(1:npx-1)
    else
       read(33,*)
    endif
    if ( npy > 1) then
       read(33,*) gpart_y(1:npy-1)
    else
       read(33,*)
    endif
    close(33)

!!$    if (sum(gpart_x(:)) /= npx) then
!!$       write(0,*) " x partition inconsistent"
!!$       STOP
!!$    endif
!!$
!!$    if (sum(gpart_x(:)) /= npy) then
!!$        write(0,*) " x partition inconsistent"
!!$       STOP
!!$    end if

  end subroutine read_input

end program grid_diagrams
