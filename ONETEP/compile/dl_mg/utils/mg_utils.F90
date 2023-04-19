! module containing useful functions for code
! testing

! 1) MPI IO dump of the grid data

! 2) 3d eigenfunction for red-black matrix

! 3) more to follow

! the eigen functions dx (i.e. 1/(nx-1) ) might need some revision
! as this subroutine was used with various grid setups.
! Lucian Anton, July 2013


module mg_utils
   !USE mg_parameters, only : wp
   implicit none
   private

   integer, parameter :: wp=kind(1.d0)
   real(wp),parameter :: pi = 4.d0*atan(1.d0)

   integer mg_comm,ngx,ngy,ngz
   ! global indexes corresponding to local 1,1,1
   integer igs,jgs,kgs

   ! eigen value index (set in test subroutines)
   integer kx,ky,kz

   ! write slice parameters set in the input file
   !integer slice_index
   !character(len=2) slice_section

   logical :: no_write=.false. ! useful when one does not want
                               ! to dump a lot of data

   public :: kx,ky,kz,mg_utils_init,write_vect, write_slice,eigenfun_rb3d, &
        eigenval_rb3d

contains

  subroutine mg_utils_init(mg_comm_,ngx_,ngy_,ngz_,local_start)
    implicit none

    integer, intent(in) :: mg_comm_,ngx_,ngy_,ngz_,local_start(3)
    integer ierr

    mg_comm = mg_comm_
    ngx = ngx_; ngy = ngy_; ngz = ngz_
    igs = local_start(1)
    jgs = local_start(2)
    kgs = local_start(3)

    ! better use local calls
    !call mpi_comm_rank(mg_comm,myid,ierr)
    !call mpi_cart_coords(mg_comm,myid,3,my_coords,ierr)
    !call mpi_cart_get(smoother_comm, 3, pgrid_dims, pgrid_periods, my_coords, ierr)

  end subroutine mg_utils_init



  subroutine write_vect(vect,gsizes,lbvec,ijkmin,ijkmax,fname_prefix,commin)
    !    use dl_mg_mpi_header

!$  use omp_lib
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif

    integer, dimension(3), intent(in) :: gsizes,lbvec, ijkmin,ijkmax
    real(wp),intent(in) :: vect(lbvec(1):,lbvec(2):,lbvec(3):)
    character(len=*),intent(in) :: fname_prefix
    integer, optional, intent(in) :: commin
    !
    !     MPI IO bits
    !
#ifdef MPI
    integer lsizes(3),lstart(3),&
            nt,nlx,nly,nlz,ftype,fh,&
            status(MPI_STATUS_SIZE),ierr
    integer myid, pdims(3),my_coords(3),comm, comm_x, comm_y, comm_z
    integer,allocatable ::  lsizes_ax(:), lsizes_ay(:), lsizes_az(:)
    logical pgrid_periods(3)

    integer(kind=MPI_OFFSET_KIND) disp !,ofs,fp_head
#endif
    character(len=128) fname

    if(no_write) return

    !
#ifdef MPI

    if(present(commin))then
       comm=commin
    else
       comm=mg_comm
    endif

    call mpi_comm_rank(comm,myid,ierr)
    call mpi_cart_get(comm, 3, pdims, pgrid_periods, my_coords, ierr)

    write(fname,'(a,i0,a,i0,a,i0)') 'npx',pdims(1),'_npy',pdims(2),'_npz',pdims(3)
!$  write(fname(len_trim(fname)+1:len(fname)),'(a,i0)') '_nth', omp_get_max_threads()
    fname=trim(fname_prefix)//adjustl(trim(fname))

    disp=0
    nlx=ijkmax(1) - ijkmin(1) + 1
    nly=ijkmax(2) - ijkmin(2) + 1
    nlz=ijkmax(3) - ijkmin(3) + 1
    nt=nlx*nly*nlz
    lsizes(:) = (/ nlx, nly, nlz /)

    ! get 3 1d subcommunicators
    call mpi_cart_sub(comm,(/ .true., .false., .false. /), comm_x, ierr)
    call mpi_cart_sub(comm,(/ .false., .true., .false. /), comm_y, ierr)
    call mpi_cart_sub(comm,(/ .false., .false., .true. /), comm_z, ierr)

    allocate (lsizes_ax(0:pdims(1)-1), lsizes_ay(0:pdims(2)-1), lsizes_az(0:pdims(3)-1))

    call mpi_allgather(lsizes(1),1,mpi_integer,lsizes_ax,1,mpi_integer,comm_x,ierr)
    call mpi_allgather(lsizes(2),1,mpi_integer,lsizes_ay,1,mpi_integer,comm_y,ierr)
    call mpi_allgather(lsizes(3),1,mpi_integer,lsizes_az,1,mpi_integer,comm_z,ierr)

    lstart(1) = sum(lsizes_ax(0:my_coords(1)-1))
    lstart(2) = sum(lsizes_ay(0:my_coords(2)-1))
    lstart(3) = sum(lsizes_az(0:my_coords(3)-1))

    write(0,*) 'write_vector ', myid, ijkmin(3), ijkmax(3), gsizes, lsizes, lstart

    call mpi_type_create_subarray(3,gsizes,lsizes,lstart,&
         MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,ftype,ierr)
    call mpi_type_commit(ftype,ierr)

    if (myid == 0) then
       call mpi_file_delete(fname,MPI_INFO_NULL,ierr)
    end if
    call mpi_barrier(comm,ierr)

! write the global sizes at the beginig of the file
    if (myid == 0) then
       call mpi_file_open(MPI_COMM_SELF,fname,&
            MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
       if (ierr /= 0) then
          write(0,*) 'MPI IO error in openning file on rank ', myid
       end if
       call mpi_file_write(fh,gsizes,3,mpi_integer, MPI_STATUS_IGNORE,ierr)
       call mpi_file_close(fh,ierr)
    end if
    call mpi_barrier(comm,ierr)

    call mpi_file_open(comm,fname,&
         MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    if (ierr /= 0) then
       write(0,*) 'MPI IO error in openning file on rank ', myid
    end if

    disp = 3*4 ! dodgy !!! get size of int
    call mpi_file_set_view(fh,disp,MPI_DOUBLE_PRECISION,ftype,&
         "native",MPI_INFO_NULL,IERR)

    call mpi_file_write_all(fh,vect(ijkmin(1):ijkmax(1),ijkmin(2):ijkmax(2),ijkmin(3):ijkmax(3)),&
         nt,MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)

    call mpi_file_close(fh,ierr)
    call mpi_type_free(ftype,ierr)
#else
    write(fname,'(a)') 'serial'
!$  write(fname(len_trim(fname)+1:len(fname)),'(a,i0)') '_nth', omp_get_max_threads()

    fname=trim(fname_prefix)//adjustl(trim(fname))
    open(55,file=fname,status="replace",form="unformatted",access="stream")
    write(55) vect(ijkmin(1):ijkmax(1),ijkmin(2):ijkmax(2),ijkmin(3):ijkmax(3))
    close(55)
#endif
  end subroutine write_vect


  ! write a slice of data in (xy,yz,xz)
  subroutine write_slice(slice_section,slice_index,sx,sy,sz,u,prefix_fn)
!    use dl_mg_mpi_header
!$    use omp_lib
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    real(wp), intent(in) :: u(:,:,:)
    integer, intent(in) ::  slice_index, sx, sy, sz
    character(len=*),intent(in) :: slice_section
    character(len=*),optional,intent(in) :: prefix_fn

    integer plane_comm, plane_dims(2),plane_coords(2),nl1,nl2,nl3,&
         ng1,ng2,ng3, perp_dir, nlxyz(3), us1,ue1,us2,ue2,us3,ue3, &
         n1,n2,plane_rank,source, plane_00, my_coords(3),pdims(3)
    integer i,j, i1,i2,i3, nlx,nly,nlz, nth, slice, ierr
    integer anl1,anl2,ags1,ags2,gs1,gs2,ijkgs(3),plane_nproc
    integer, allocatable :: plane_loc_info(:,:)
    logical subdims(3), plane_periods(2), pgrid_periods(3)
    real(wp), allocatable :: buff_plane(:,:), buff(:,:)
    real(wp) q
    character fn*128, fn_aux*16, plane*2

    if(no_write .or. slice_index < 0) return

    slice = slice_index
    plane = slice_section

#ifdef MPI
    call mpi_cart_get(mg_comm, 3, pdims, pgrid_periods, my_coords, ierr)
#else
    pdims = (/ 1, 1, 1 /)
    pgrid_periods = (/ .false., .false., .false. /)
    my_coords = (/ 0, 0, 0 /)
#endif

    nlx = size(u,1); nly = size(u,2); nlz = size(u,3)
    nlxyz = (/ nlx, nly, nlz /)
    ijkgs = (/ sx, sy, sz /)

    select case(plane)
    case( "xy","yx")
       subdims=(/.true., .true., .false. /)
       perp_dir = 3
       nl1 = nlx; nl2 = nly
       us1 = 1; ue1 = nlx
       us2 = 1; ue2 = nly
       us3 = slice - sz + 1
       ue3 = us3
       gs1 = sx; gs2 = sy
    case("yz","zy")
       subdims=(/.false., .true., .true. /)
       perp_dir = 1
       nl1 = nly; nl2 = nlz
       us1 = slice - sx + 1
       ue1 = us1
       us2 = 1; ue2 = nly
       us3 = 1; ue3 = nlz
       gs1 = sy; gs2 = sz
    case("xz","zx")
       subdims=(/.true., .false., .true. /)
       perp_dir = 2
       nl1 = nlx; nl2 = nlz
       us1 = 1; ue1 = nlx
       us2 = slice - sy
       ue2 = us2
       us3 = 1; ue3 = nlz
       gs1 = sx; gs2 = sz
    end select
#ifdef MPI
    call MPI_Cart_sub(mg_comm,subdims,plane_comm,ierr)
#endif
    if ( .not. (slice >= ijkgs(perp_dir) .and. &
         slice < ijkgs(perp_dir)+nlxyz(perp_dir))) return

    ! only the ranks containing the slice work
#ifdef MPI
    call mpi_cart_get(plane_comm,2,plane_dims,plane_periods,&
         plane_coords, ierr)
    call mpi_comm_rank(plane_comm,plane_rank,ierr)
    call mpi_comm_size(plane_comm,plane_nproc,ierr)

    ! set the ranks that collects the data
    plane_00 = -1
    if (plane_coords(1) == 0 .and. plane_coords(2) == 0 ) then
       plane_00 = plane_rank
    endif

    call MPI_Allreduce(MPI_IN_PLACE,plane_00,1,MPI_INTEGER,MPI_MAX, plane_comm,ierr)
#else
    plane_dims = (/1, 1/)
    plane_periods = (/ .false., .false./)
    plane_coords = (/ 0, 0 /)
    plane_rank = 0; plane_00 = 0; plane_nproc=1
#endif
    if (plane_rank == plane_00) then

       allocate (plane_loc_info(4,0:plane_nproc-1))
#ifdef MPI
       call mpi_gather( (/nl1,nl2,gs1,gs2/),4,mpi_integer,plane_loc_info,4,mpi_integer,plane_00,plane_comm,ierr)
#else
     plane_loc_info(:,0) = (/nl1,nl2,gs1,gs2/)
#endif

     !write(0,*) 'write slice info ', plane_loc_info

     n1 = 0
     do i = 0, plane_dims(1)-1
#ifdef MPI
        call mpi_cart_rank(plane_comm,(/i,0/),source,ierr)
#else
        source = 0
#endif
        n1 = n1 + plane_loc_info(1, source)
     enddo
     n2 = 0
     do i = 0, plane_dims(2)-1
#ifdef MPI
        call mpi_cart_rank(plane_comm,(/0,i/),source,ierr)
#else
        source = 0
#endif
        n2 = n2 + plane_loc_info(2, source)
     enddo

     allocate (buff_plane(n1,n2))
       ! get cart dims
       !write(0,*) 'plane_dims ', n1, n2, plane_dims, plane_00

     !!! WARNING buff_plane assumes that min(sx)=1, and min(sy)=1
       buff_plane(1:nl1, 1:nl2) = &
            reshape(u(us1:ue1,us2:ue2,us3:ue3), (/ nl1, nl2/))

       do j = 0, plane_dims(2)-1
          do i = 0, plane_dims(1)-1
             if ( i == 0 .and. j == 0 ) cycle
#ifdef MPI
             call mpi_cart_rank(plane_comm,(/i,j/),source,ierr)
#endif
             !write(0,*) 'source  ', source
             anl1 = plane_loc_info(1,source)
             anl2 = plane_loc_info(2,source)
             allocate(buff(anl1,anl2))
             ags1 = plane_loc_info(3,source)
             ags2 = plane_loc_info(4,source)
#ifdef MPI
             call mpi_recv(buff,anl1*anl2,MPI_DOUBLE_PRECISION,source,1,&
                  plane_comm,mpi_status_ignore,ierr)
#endif

             !write(0,*) 'ags..', ags1, ags2, anl1, anl2
             buff_plane(ags1:ags1+anl1-1,ags2:ags2+anl2-1 ) = &
                  buff(:,:)
             deallocate(buff)
          enddo
       enddo

       ! build file name for slice data
       ! slice coordinates
       write(fn_aux,*) 'slice_'; fn = adjustl(fn_aux)
       write(fn_aux,*) plane;    fn = fn(1:len_trim(fn))//adjustl(fn_aux)
       write(fn_aux,*) slice;    fn = fn(1:len_trim(fn))//adjustl(fn_aux)
!!$         write(fn_aux,*) us1; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
!!$         write(fn_aux,*) 'x'; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
!!$         write(fn_aux,*) ue1; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
!!$         write(fn_aux,*) '_'; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
!!$         write(fn_aux,*) us2; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
!!$         write(fn_aux,*) 'x'; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
!!$         write(fn_aux,*) ue2; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
!!$         write(fn_aux,*) '_'; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
!!$         write(fn_aux,*) us3; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
!!$         write(fn_aux,*) 'x'; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
!!$         write(fn_aux,*) ue3; fn = fn(1:len_trim(fn))//adjustl(fn_aux)

       !processor grid info
       write(fn_aux,*) '_npx'; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
       write(fn_aux,*) pdims(1); fn = fn(1:len_trim(fn))//adjustl(fn_aux)
       write(fn_aux,*) '_npy'; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
       write(fn_aux,*) pdims(2); fn = fn(1:len_trim(fn))//adjustl(fn_aux)
       write(fn_aux,*) '_npz'; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
       write(fn_aux,*) pdims(3); fn = fn(1:len_trim(fn))//adjustl(fn_aux)

       ! iteration
       !write(fn_aux,*) '_iter'; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
       !write(fn_aux,*) iter; fn = fn(1:len_trim(fn))//adjustl(fn_aux)

       if (present(prefix_fn)) then
          fn=trim(prefix_fn)//fn
       endif

       ! add number of threads if using OpenMP
       !$       nth = omp_get_max_threads()
       !$       write(fn_aux,*) '_nth'; fn = fn(1:len_trim(fn))//adjustl(fn_aux)
       !$       write(fn_aux,*) nth; fn = fn(1:len_trim(fn))//adjustl(fn_aux)

       ! analytical solution prefactor

       !q = (1.0_wp/norm_eigenfun) * eigenval_rb3d()**iter

       ! use layout used by gnuplot splot command
       open (56,file=fn,status='replace')
       do j = 1, n2
          do i = 1, n1
             select case(plane)
             case("xy")
                i1 = i; i2 = j; i3 = slice
             case("yz")
                i1 = slice; i2 = i; i3 = j
             case("xz")
                i1 = i; i2 = slice; i3 = j
             end select

             write(56,'(1000E24.15)') buff_plane(i,j) !, q*eigenfun_rb3d(i1,i2,i3)
          enddo
          write(56,*) ! empty record for gnuplot
       end do
       close(56)

    else
#ifdef MPI
       call mpi_gather( (/nl1,nl2,gs1,gs2/),4,mpi_integer,i,4,mpi_integer,plane_00,plane_comm,ierr)

       call mpi_send(u(us1:ue1,us2:ue2,us3:ue3),&
            size(u(us1:ue1,us2:ue2,us3:ue3)),&
            MPI_DOUBLE_PRECISION,plane_00,1,plane_comm,ierr)
#endif
    endif
#ifdef MPI
    call mpi_comm_free(plane_comm,ierr)
#endif

  end subroutine write_slice


!
!----------------------------------------------------------------------
!
! egeinfunction for red-black operator and zero boundary condition
  function eigenfun_rb3d(i,j,k) result(u)
    !USE mg_parameters, only : wp, pi
    implicit none

    integer, intent(in) :: i,j,k!,kx,ky,kz,ngx,ngy,ngz

    real(wp)  u ! retun value

    ! locals
    real(wp) rx,ry,rz

    rx = cos(kx*pi/real(ngx+1,wp))
    ry = cos(ky*pi/real(ngy+1,wp))
    rz = cos(kz*pi/real(ngz+1,wp))

    u = (rx * (1.0_wp + rx) * sin(real(kx*i,wp)*pi/real(ngx+1,wp)) + &
         rx * (1.0_wp - rx) * sin(real((ngx+1-kx)*i,wp)*pi/real(ngx+1,wp))) * &
        (ry * (1.0_wp + ry) * sin(real(ky*j,wp)*pi/real(ngy+1,wp)) + &
         ry * (1.0_wp - ry) * sin(real((ngy+1-ky)*j,wp)*pi/real(ngy+1,wp))) * &
        (rz * (1.0_wp + rz) * sin(real(kz*k,wp)*pi/real(ngz+1,wp)) + &
         rz * (1.0_wp - rz) * sin(real((ngz+1-kz)*k,wp)*pi/real(ngz+1,wp)))


  end function eigenfun_rb3d


! egeinvalues  for red-black operator and zero boundary condition
  function eigenval_rb3d() result(z)
    !USE mg_parameters, only : wp, pi
    implicit none

    !integer, intent(in) :: kx,ky,kz,ngx,ngy,ngz

    real(wp)  z ! return value

    z = (1.d0/3.d0)*(cos((kx*pi)/real(ngx+1,wp))**2 + &
         cos((ky*pi)/real(ngy+1,wp))**2 + cos((kz*pi)/real(ngz+1,wp))**2)

  end function eigenval_rb3d

  end module mg_utils
