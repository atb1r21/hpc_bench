!
!> \brief  utilities for  MPI communications, provides
!! specialised versions for serial case
!!
!
!! Lucian Anton
!! February 2013

module dl_mg_mpi
  use dl_mg_params
  use dl_mg_mpi_header
  implicit none

  ! mode constants for halo exchage
  integer, parameter :: exch_rl=1000, exch_rr=100, exch_sl=10, exch_sr=1, &
  exch_full = exch_rl + exch_rr + exch_sl + exch_sr

  ! arrays used in halo exchage
  real(wp), allocatable, private :: buffrecv_xl(:), buffrecv_xr(:),&
                                    buffsend_xl(:), buffsend_xr(:),&
                                    buffrecv_yl(:), buffrecv_yr(:),&
                                    buffsend_yl(:), buffsend_yr(:),&
                                    buffrecv_zl(:), buffrecv_zr(:),&
                                    buffsend_zl(:), buffsend_zr(:)

  ! request array used in halo exchange
  integer, private :: exrq(12), nrq

contains

  subroutine set_mpi_comm(commin, ierror)
    use dl_mg_common_data, only : mg_comm, bc
    use dl_mg_errors
    implicit none

    integer, intent(in) :: commin
    integer, optional, intent(inout) :: ierror

    integer topo, nproc, dims(3), coords(3), i, ierr
    logical periods(3)
    character(len=255) errmsg

    if (present(ierror)) ierror = DL_MG_SUCCESS
#ifdef MPI
    call mpi_topo_test(commin,topo,ierr)
    if (.not. check_assertion(topo == MPI_CART, DL_MG_ERR_MPI_TOPO)) then
       call handle_error(DL_MG_ERR_MPI_TOPO, ierror)
        return ! if handle_error does not abort the execution control must be
               !return to the caller

        !assume 1 d decomposition in z direction, following onetep
        !call mpi_comm_size(commin,nproc,ierr)
        !call mpi_cart_create(commin, 3, (/1,1,nproc/), (/.false., .false., .false./), .false.,mg_comm,ierr)
     else
        call mpi_comm_dup(commin,mg_comm,ierr)
        call mpi_topo_test(mg_comm,topo,ierr)
        if ( .not. check_assertion(topo == MPI_CART, DL_MG_ERR_MPI_COMM_DUP)) then
           call handle_error(DL_MG_ERR_MPI_COMM_DUP, ierror)
           return
        endif
     endif

     ! test if BC are consistent with MPI topology periodicity
      call mpi_cart_get(mg_comm,3,dims,periods,coords,ierr)
      do i = 1, 3
         if ( .not. check_assertion((periods(i) .and. (bc(i) == DL_MG_BC_PERIODIC)) .or. &
              (.not. periods(i) .and. (bc(i) == DL_MG_BC_DIRICHLET)), &
              DL_MG_ERR_MPI_TOPO_BC) ) then
               write(errmsg, '(a, a, I3, a, L1, a, I3)') 'boundary condition inconsistent', &
                    "with MPI topology in direction",&
                    i, " MPI periods ", periods(i),&
                    " BC", bc(i)
               call handle_error( DL_MG_ERR_MPI_TOPO_BC, ierror, errmsg)
               return
            endif
      enddo

#endif
  end subroutine set_mpi_comm


  subroutine set_nprocs
     use dl_mg_common_data, only : mg_comm,npx,npy,npz
     implicit none

#ifdef MPI
     integer dims(3), ierr

     call get_mpi_grid(mg_comm, dims=dims)
     npx =dims(1); npy = dims(2); npz = dims(3)
#else
      npx=1 ; npy = 1; npz = 1
#endif
    end subroutine set_nprocs


    integer function get_myid(comm) result(myid)
      use dl_mg_common_data, only : mg_comm
      implicit none

      integer, optional, intent(in) :: comm

      integer comm_, ierr

      if (present(comm)) then
         comm_ =comm
      else
         comm_ = mg_comm ! use the app passed communicator
      end if

#ifdef MPI
    call mpi_comm_rank(comm_, myid, ierr)
#else
    myid = 0
#endif

    end function get_myid

   subroutine set_grid_comm(commin,active,comm)
     use dl_mg_errors
     implicit none
     integer, intent(in)     :: commin
     logical, intent(in)     :: active
     integer, intent(out)    :: comm

#ifdef MPI

     integer color,min_color,id,aux4(4),i,j,k,p,np,px,py,pz, &
          comm0,dims(3),coords_in(3),ierr
     integer, allocatable :: col_array(:,:), grid_in(:,:,:)
     logical periods(3)


     comm = MPI_COMM_NULL

     ! this rank was eliminated at a higher level
     if (commin == MPI_COMM_NULL) then
        return
     endif

     if (active) then
        color = 1
     else
        color = 0
     endif

     call mpi_allreduce(color, min_color, 1, mpi_integer, MPI_MIN, commin, ierr)
     if (min_color == 1) then
        ! no spliting needed
        ! is this really necessary?
        ! is split with one coour equivalent with dup?
        call mpi_comm_dup(commin,comm,ierr)
     else
        call mpi_comm_rank(commin,id,ierr)
        if (color == 0) color = MPI_UNDEFINED
        call mpi_comm_split(commin,color,id,comm0,ierr)
        if (color == MPI_UNDEFINED) color = 0

        ! set a cartesian topology at this level
        call mpi_comm_size(commin,np,ierr)
        allocate(col_array(4,0:np-1))
        call get_mpi_grid(commin,periods=periods,dims=dims,coords=coords_in)

        ! I wonder if this version is scalable
        ! alternative : use subcomms ?
        allocate(grid_in(0:dims(1)-1,0:dims(2)-1,0:dims(3)-1))
        aux4(1:3) = coords_in
        aux4(4)   = color
        call mpi_allgather(aux4,4,mpi_integer,col_array,4,mpi_integer,commin,ierr)

        do p = 0, np - 1
           i = col_array(1,p)
           j = col_array(2,p)
           k = col_array(3,p)
           grid_in(i, j, k) = col_array(4,p)

        enddo

        ! find a site with color == 1
        px = 0; py = 0; pz = 0
        ext_loop : do k = 0, dims(3) -1
           do j = 0, dims(2) -1
              do i = 0, dims(1) -1
                 if (grid_in(i,j,k) == 1) then
                    px = sum(grid_in(:, j, k))
                    py = sum(grid_in(i, :, k))
                    pz = sum(grid_in(i, j, :))
                    exit ext_loop
                 endif
              enddo
           enddo
        enddo ext_loop

        if ( .not. check_assertion(px > 0 .and. py > 0 .and. pz > 0)) then
           call handle_error(DL_MG_ERR_UNSPECIFIED, &
                msg = "set_grid_comm: wrong number of ranks in new topology")
        endif

        if (color == 1 ) then
           call mpi_cart_create(comm0, 3, (/ px, py, pz /), periods, .false., comm, ierr)
        endif
     endif
#else
     comm = -1
#endif

   end subroutine set_grid_comm

   !> return info on MPI topology
   subroutine get_mpi_grid(comm,rank,dims,periods,coords)
     implicit none
     integer, intent(in) :: comm
     integer, intent(out), optional :: rank, coords(3), dims(3)
     logical, intent(out), optional :: periods(3)

     integer ierr, dims_(3), coords_(3)
     logical periods_(3)

     if (present(rank)) then
#ifdef MPI
        call mpi_comm_rank(comm,rank,ierr)
#else
        rank=0
#endif
     endif

     if (present(coords) .or. &
          present (dims)  .or. &
          present(periods)) then
#ifdef MPI
        call mpi_cart_get(comm,3,dims_,periods_,coords_,ierr)
#else
        dims_ = (/ 1, 1, 1/)
        periods_ = (/ .false., .false., .false. /)
        coords_ = (/ 0, 0, 0 /)
#endif
        if (present(coords))  coords  = coords_
        if (present(dims))    dims    = dims_
        if (present(periods)) periods = periods_
     endif
   end subroutine get_mpi_grid


   subroutine find_mpi_neighbors(mg)
     use dl_mg_types
     implicit none
     type(mg_t), intent(inout) :: mg

#ifdef MPI
     integer ngb(3),dims(3),coords(3),ierr
     logical periods(3)

     !
     ! neighbors which share a plane with 'myid'
     !

     if (mg%comm == MPI_COMM_NULL ) return

      call MPI_CART_GET(mg%comm,3,dims,periods,coords,ierr)

      ! north
      ngb(1)=coords(1)-1
      ngb(2)=coords(2)
      ngb(3)=coords(3)
      if (coords(1) == 0 .and. .not. periods(1)) then
         mg%neighbor(north)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(north),ierr)
      endif

      ! south
      ngb(1)=coords(1)+1
      ngb(2)=coords(2)
      ngb(3)=coords(3)
      if (coords(1)==dims(1)-1 .and. .not. periods(1)) then
         mg%neighbor(south)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(south),ierr)
      endif

      ! west
      ngb(1)=coords(1)
      ngb(2)=coords(2)-1
      ngb(3)=coords(3)
      if ( coords(2) == 0 .and. .not. periods(2)) then
         mg%neighbor(west)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(west),ierr)
      endif

      ! east
      ngb(1)=coords(1)
      ngb(2)=coords(2)+1
      ngb(3)=coords(3)
      if (coords(2) == dims(2)-1 .and. .not. periods(2)) then
         mg%neighbor(east) = MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(east),ierr)
      endif

      ! bottom
      ngb(1)=coords(1)
      ngb(2)=coords(2)
      ngb(3)=coords(3)-1
      if (coords(3) == 0 .and. .not. periods(3)) then
         mg%neighbor(bottom) = MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(bottom),ierr)
      endif


      ! top
      ngb(1)=coords(1)
      ngb(2)=coords(2)
      ngb(3)=coords(3)+1
      if(coords(3) == dims(3)-1 .and. .not. periods(3))then
         mg%neighbor(top)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(top),ierr)
      endif

!      write(0,*) 'did planes'
!
! neighbors which share a line with 'myid'
!
      ! nw
      ngb(1)=coords(1)-1
      ngb(2)=coords(2)-1
      ngb(3)=coords(3)
       if ( (coords(1) == 0 .and..not. periods(1)) &
            .or. (coords(2) == 0 .and..not. periods(2))) then
         mg%neighbor(nw)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(nw),ierr)
      endif

      ! ne
      ngb(1)=coords(1)-1
      ngb(2)=coords(2)+1
      ngb(3)=coords(3)
      if ( (coords(1) == 0 .and..not.periods(1))&
           .or. (coords(2) == dims(2)-1 .and..not.periods(2))) then
         mg%neighbor(ne)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(ne),ierr)
      endif

      ! nb
      ngb(1)=coords(1)-1
      ngb(2)=coords(2)
      ngb(3)=coords(3)-1
      if ( (coords(1) == 0 .and..not. periods(1)) &
           .or. (coords(3) == 0 .and..not. periods(3)) )then
         mg%neighbor(nb)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(nb),ierr)
      endif

      ! nt
      ngb(1)=coords(1)-1
      ngb(2)=coords(2)
      ngb(3)=coords(3)+1
      if ( (coords(1) == 0 .and..not. periods(1)) &
           .or. (coords(3) == dims(3)-1 .and..not. periods(3)))then
         mg%neighbor(nt)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(nt),ierr)
      endif

      ! sw
      ngb(1)=coords(1)+1
      ngb(2)=coords(2)-1
      ngb(3)=coords(3)
      if (coords(1) == dims(1)-1 .and..not. periods(1) &
           .or. (coords(2) == 0 .and..not. periods(2))) then
         mg%neighbor(sw)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(sw),ierr)
      endif

      ! se
      ngb(1)=coords(1)+1
      ngb(2)=coords(2)+1
      ngb(3)=coords(3)
      if ( (coords(1)==dims(1)-1.and..not.periods(1))&
           .or. (coords(2)==dims(2)-1.and..not.periods(2)))then
         mg%neighbor(se)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(se),ierr)
      endif

      ! sb
      ngb(1)=coords(1)+1
      ngb(2)=coords(2)
      ngb(3)=coords(3)-1
       if ( (coords(1) == dims(1)-1 .and..not. periods(1))&
           .or. (coords(3) == 0 .and..not. periods(3)))then
         mg%neighbor(sb)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(sb),ierr)
      endif

      ! st
      ngb(1)=coords(1)+1
      ngb(2)=coords(2)
      ngb(3)=coords(3)+1
      if ( (coords(1)==dims(1)-1.and..not.periods(1))&
           .or. (coords(3)==dims(3)-1.and..not.periods(3)))then
         mg%neighbor(st)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(st),ierr)
      endif

      ! wb
      ngb(1)=coords(1)
      ngb(2)=coords(2)-1
      ngb(3)=coords(3)-1
      if ( (coords(2) == 0 .and..not. periods(2))&
           .or. (coords(3) == 0 .and..not. periods(3))) then
         mg%neighbor(wb)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(wb),ierr)
      endif

      ! wt
      ngb(1)=coords(1)
      ngb(2)=coords(2)-1
      ngb(3)=coords(3)+1
      if ( (coords(2) == 0 .and..not. periods(2))&
           .or. (coords(3) == dims(3)-1 .and..not. periods(3))) then
         mg%neighbor(wt)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(wt),ierr)
      endif

      ! eb
      ngb(1)=coords(1)
      ngb(2)=coords(2)+1
      ngb(3)=coords(3)-1
      if ( (coords(2) == dims(2)-1 .and..not. periods(2)) &
           .or. (coords(3) == 0 .and..not. periods(3))) then
         mg%neighbor(eb)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(eb),ierr)
      endif

      ! et
      ngb(1)=coords(1)
      ngb(2)=coords(2)+1
      ngb(3)=coords(3)+1
      if ( (coords(2) == dims(2)-1 .and..not. periods(2))&
           .or. (coords(3) == dims(3)-1 .and..not. periods(3))) then
         mg%neighbor(et)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(et),ierr)
      endif

!      write(0,*) 'did  lines'

!
! neighbors which share a corner with 'myid'
!
      ! nwb
      ngb(1)=coords(1)-1
      ngb(2)=coords(2)-1
      ngb(3)=coords(3)-1
      if (      (coords(1)==0        .and..not.periods(1))&
           .or. (coords(2)==0         .and..not.periods(2))&
           .or. (coords(3)==0        .and..not.periods(3))) then
         mg%neighbor(nwb)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(nwb),ierr)
      endif

      ! neb
      ngb(1)=coords(1)-1
      ngb(2)=coords(2)+1
      ngb(3)=coords(3)-1
      if (      (coords(1)==0        .and..not.periods(1))&
           .or. (coords(2)==dims(2)-1.and..not.periods(2))&
           .or. (coords(3)==0        .and..not.periods(3))) then
         mg%neighbor(neb)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(neb),ierr)
      endif

      ! nwt
      ngb(1)=coords(1)-1
      ngb(2)=coords(2)-1
      ngb(3)=coords(3)+1
      if (      (coords(1)==0        .and..not.periods(1))&
           .or. (coords(2)==0        .and..not.periods(2))&
           .or. (coords(3)==dims(3)-1.and..not.periods(3))) then
         mg%neighbor(nwt)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(nwt),ierr)
      endif

      !net
      ngb(1)=coords(1)-1
      ngb(2)=coords(2)+1
      ngb(3)=coords(3)+1
      if (      (coords(1)==0        .and..not.periods(1))&
           .or. (coords(2)==dims(2)-1.and..not.periods(2))&
           .or. (coords(3)==dims(3)-1.and..not.periods(3))) then
         mg%neighbor(net)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(net),ierr)
      endif

      ! swb
      ngb(1)=coords(1)+1
      ngb(2)=coords(2)-1
      ngb(3)=coords(3)-1
      if (      (coords(1)==dims(1)-1.and..not.periods(1))&
           .or. (coords(2)==0        .and..not.periods(2))&
           .or. (coords(3)==0        .and..not.periods(3))) then
         mg%neighbor(swb)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(swb),ierr)
      endif

      ! seb
      ngb(1)=coords(1)+1
      ngb(2)=coords(2)+1
      ngb(3)=coords(3)-1
      if (      (coords(1)==dims(1)-1.and..not.periods(1))&
           .or. (coords(2)==dims(2)-1.and..not.periods(2))&
           .or. (coords(3)==0        .and..not.periods(3))) then
         mg%neighbor(seb)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(seb),ierr)
      endif

      ! swt
      ngb(1)=coords(1)+1
      ngb(2)=coords(2)-1
      ngb(3)=coords(3)+1
      if (      (coords(1)==dims(1)-1.and..not.periods(1))&
           .or. (coords(2)==0        .and..not.periods(2))&
           .or. (coords(3)==dims(3)-1.and..not.periods(3))) then
         mg%neighbor(swt)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(swt),ierr)
      endif

      ! set
      ngb(1)=coords(1)+1
      ngb(2)=coords(2)+1
      ngb(3)=coords(3)+1
      if (      (coords(1)==dims(1)-1.and..not.periods(1))&
           .or. (coords(2)==dims(2)-1.and..not.periods(2))&
           .or. (coords(3)==dims(3)-1.and..not.periods(3))) then
         mg%neighbor(set)=MPI_PROC_NULL
      else
         call MPI_CART_RANK(mg%comm,ngb,mg%neighbor(set),ierr)
      endif

!      if(coords(3) == dims(3)-1)then
!         write(0,*) 'set topology', coords(3), bd
!      endif
#endif
! ifdef MPI

   end subroutine find_mpi_neighbors


   subroutine exchange_halo_begin(dir, mode, thickness, mg, component)
     use dl_mg_types
     use dl_mg_errors
     implicit none
     integer, intent(in)          :: dir, mode, thickness
     type(mg_t), target, intent(in)       :: mg
     character(len=1), intent(in) :: component

#ifdef MPI
     integer asx, sx, aex, ex, asy,sy, aey, ey, asz, sz, aez, ez
     integer left, right ! MPI neighbours
     real(wp), pointer :: pp(:,:,:) ! for p and z components
!#ifdef HAVE_CONTIGUOUS
!     contiguous pp
!#endif

     ! when starting new exchange there must not be another one in process
     ! need to think more about it
!!$     do i =1, 12
!!$        if ( exrq /= MPI_REQUEST_NULL) then
!!$           call error_abort("exchange_halo_begin: found another halo exchage in progress")
!!$        endif
!!$     enddo

     ! current index in request array
     nrq = 0
     call mg_get_idx(mg,sx,ex,sy,ey,sz,ez)
     select case(dir)
        case (1)
           left  = mg%neighbor(north)
           right = mg%neighbor(south)
        case(2)
           left  = mg%neighbor(west)
           right = mg%neighbor(east)
        case(3)
           left  = mg%neighbor(bottom)
           right = mg%neighbor(top)
        case default
           if (check_assertion(.true.))then
              call handle_error(DL_MG_ERR_UNSPECIFIED, &
                   msg="wrong value of dir in exchange_halo_begin")
           endif
        end select

     select case ( component)
     case ("p","z")
        asx = mg%psx
        aex = mg%pex
        asy = mg%psy
        aey = mg%pey
        asz = mg%psz
        aez = mg%pez
        select case(component)
        case("p")
           pp => mg%p
        case("z")
           pp => mg%z
        end select
        call exchange_halo_begin_base(dir, left, right, mg%comm, mode, thickness, &
             asx, asy, asz, aex, aey, aez, 1,&
             sx, ex, sy, ey, sz, ez, pp)
     case("w") ! nonlinear only in restrict
        asx = mg%psx
        aex = mg%pex
        asy = mg%psy
        aey = mg%pey
        asz = mg%psz
        aez = mg%pez
        call exchange_halo_begin_base(dir, left, right, mg%comm, mode, thickness, &
             asx, asy, asz, aex, aey, aez, 1,&
             sx, ex, sy, ey, sz, ez, mg%w)
     case("f")
        asx = mg%fsx
        aex = mg%fex
        asy = mg%fsy
        aey = mg%fey
        asz = mg%fsz
        aez = mg%fez
        call exchange_halo_begin_base(dir, left, right, mg%comm, mode, thickness,&
             asx, asy, asz, aex, aey, aez, 1, &
             sx, ex, sy, ey, sz, ez, mg%f)
     case("c")
        asx = mg%csx
        aex = mg%cex
        asy = mg%csy
        aey = mg%cey
        asz = mg%csz
        aez = mg%cez
        call exchange_halo_begin_base(dir, left, right, mg%comm, mode, thickness,&
             asx, asy, asz, aex, aey, aez, 3, &
             sx, ex, sy, ey, sz, ez, mg%c)
     case("r")
        asx = mg%rsx
        aex = mg%rex
        asy = mg%rsy
        aey = mg%rey
        asz = mg%rsz
        aez = mg%rez
        call exchange_halo_begin_base(dir, left, right, mg%comm, mode, thickness,&
             asx, asy, asz, aex, aey, aez, 1, &
             sx, ex, sy, ey, sz, ez, mg%r)
     case default
        if (check_assertion(.true.)) then
           call handle_error(DL_MG_ERR_UNSPECIFIED, &
                msg = "wrong mg component value in exchange_halo_begin")
        endif
     end select
#endif
   end subroutine exchange_halo_begin


   subroutine exchange_halo_begin_base(dir, left, right, comm, mode, thickness, &
        asx, asy, asz, aex, aey, aez, ncomp,&
        sx, ex, sy, ey, sz, ez, a)
     use dl_mg_types
     use dl_mg_errors
     implicit none
     integer, intent(in) ::dir, left, right, comm, mode, thickness
     integer, intent(in) :: asx, asy, asz, aex, aey, aez, ncomp, &
          sx, ex, sy, ey, sz, ez
     real(wp), intent(in) :: a(asx:aex, asy:aey, asz:aez, ncomp)

#ifdef MPI
     integer n, ierr, dsx,dex,dsy,dey,dsz,dez
     character(len=4) cmode

     !write(0,*) 'halo begin', dir, left, right, MPI_PROC_NULL

     ! mode is any combination of 0,1 in 4 positions
     write(cmode,'(I4)') mode

     select case(dir)
     case(1)
        dsx = thickness
        dex = thickness
        dsy = min(thickness,sy-asy)
        dey = min(thickness,aey-ey)
        dsz = min(thickness,sz-asz)
        dez = min(thickness,aez-ez)

        if (cmode(1:1) == "1" ) then
           nrq = nrq + 1
           if (left /= MPI_PROC_NULL) then
              ! receive from left
              n =  dsx * (ey - sy + 1 + dsy + dey) * &
                         ( ez -sz + 1 + dsz + dez) * ncomp
              allocate(buffrecv_xl(n))
              call mpi_irecv(buffrecv_xl, n, MPI_DOUBLE_PRECISION, left, 701,&
                   comm, exrq(nrq),ierr)
           else
              exrq(nrq) = MPI_REQUEST_NULL
           endif
        endif

        if (cmode(2:2) == "1") then
           nrq = nrq + 1
           if (right /= MPI_PROC_NULL) then
              ! receive from right
              n =  dex * (ey - sy + 1 + dsy + dey) * &
                         (ez - sz + 1 + dsz + dez) * ncomp
              allocate(buffrecv_xr(n))
              call mpi_irecv(buffrecv_xr, n, MPI_DOUBLE_PRECISION, right, 702,&
                   comm, exrq(nrq), ierr)
           else
              exrq(nrq) = MPI_REQUEST_NULL
           endif
        endif

        if (cmode(3:3) == "1") then
           nrq = nrq + 1
           if (left /= MPI_PROC_NULL) then
              ! send left
              n =  dsx * (ey - sy + 1 + dsy + dey) * &
                         (ez - sz + 1 + dsz + dez) * ncomp
              allocate(buffsend_xl(n))
              buffsend_xl(:) = reshape(a(sx:sx+dsx-1, sy-dsy:ey+dey, sz-dsz:ez+dez, 1:ncomp), (/ n /))
              call mpi_isend(buffsend_xl, n, MPI_DOUBLE_PRECISION, left, 702,&
                   comm, exrq(nrq), ierr)
              else
                 exrq(nrq) = MPI_REQUEST_NULL
              endif
        endif

        if (cmode(4:4) == "1") then
           nrq = nrq + 1
           if (right /= MPI_PROC_NULL) then
              ! send right
              n =  dex * (ey - sy + 1 + dsy + dey) * &
                         (ez - sz + 1 + dsz + dez) * ncomp
              allocate(buffsend_xr(n))
              buffsend_xr(:) = reshape(a(ex-dex+1:ex, sy-dsy:ey+dey, sz-dsz:ez+dez, 1:ncomp), (/ n /))
              call mpi_isend(buffsend_xr, n, MPI_DOUBLE_PRECISION, right, 701,&
                   comm, exrq(nrq), ierr)
              else
                 exrq(nrq) = MPI_REQUEST_NULL
              endif
        endif

     case(2)
        dsy = thickness
        dey = thickness
        dsx = min(thickness,sx-asx)
        dex = min(thickness,aex-ex)
        dsz = min(thickness,sz-asz)
        dez = min(thickness,aez-ez)

        if ( cmode(1:1) == "1"  ) then
           nrq = nrq + 1
           if (left /= MPI_PROC_NULL ) then
              ! receive from left
              n = dsy * (ex - sx + 1 + dsx + dex) * &
                        (ez - sz + 1 + dsz + dez) * ncomp
              allocate(buffrecv_yl(n))
              call mpi_irecv(buffrecv_yl, n, MPI_DOUBLE_PRECISION, left, 705,&
                   comm, exrq(nrq), ierr)
           else
              exrq(nrq) = MPI_REQUEST_NULL
          endif
        endif

        if ( cmode(2:2) == "1" ) then
           nrq = nrq + 1
           if ( right /= MPI_PROC_NULL ) then
              ! receive from right
              n = dey * (ex - sx + 1 + dsx + dex) * &
                        (ez - sz + 1 + dsz + dez) * ncomp
              allocate(buffrecv_yr(n))
              call mpi_irecv(buffrecv_yr, n, MPI_DOUBLE_PRECISION,  right, 706, &
                   comm, exrq(nrq),ierr)
           else
              exrq(nrq) = MPI_REQUEST_NULL
           endif
        endif

        if (cmode(3:3) == "1" ) then
           nrq = nrq + 1
           if (left /= MPI_PROC_NULL) then
              ! send to left
              n =  dsy * (ex - sx + 1 + dsx + dex) * &
                         (ez - sz + 1 + dsz + dez) * ncomp
              allocate(buffsend_yl(n))
              buffsend_yl(:) = reshape( a(sx-dsx:ex+dex, sy:sy+dsy-1, sz-dsz:ez+dez, 1:ncomp), (/ n /))
              call mpi_isend(buffsend_yl, n, MPI_DOUBLE_PRECISION, left, 706,&
                   comm, exrq(nrq), ierr)
           else
              exrq(nrq) = MPI_REQUEST_NULL
           endif
        endif


        if (cmode(4:4) == "1" ) then
           nrq = nrq + 1
           if (right /= MPI_PROC_NULL) then
              ! send to right
              n =  dey * (ex - sx + 1 + dsx + dex) * &
                         (ez - sz + 1 + dsz + dez) * ncomp
              allocate(buffsend_yr(n))
              buffsend_yr(:) = reshape(a(sx-dsx:ex+dex, ey-dey+1:ey, sz-dsz:ez+dez, 1:ncomp), (/ n /))
              call mpi_isend(buffsend_yr, n, MPI_DOUBLE_PRECISION, right, 705, &
                   comm, exrq(nrq), ierr)
           else
              exrq(nrq) = MPI_REQUEST_NULL
           endif
        endif

     case(3)
        dsz = thickness
        dez = thickness
        dsx = min(thickness,sx-asx)
        dex = min(thickness,aex-ex)
        dsy = min(thickness,sy-asy)
        dey = min(thickness,aey-ey)

        if ( cmode(1:1) == "1" ) then
           nrq = nrq + 1
           if (left /= MPI_PROC_NULL ) then
              ! receive from left
              n = dsz * (ex - sx + 1 + dsx + dex) * &
                        (ey - sy + 1 + dsy + dey) * ncomp
              allocate(buffrecv_zl(n))
              call mpi_irecv(buffrecv_zl, n, MPI_DOUBLE_PRECISION, left, 709,&
                   comm, exrq(nrq),ierr)
           else
               exrq(nrq) = MPI_REQUEST_NULL
           endif

        endif

        if ( cmode(2:2) == "1" ) then
           nrq = nrq + 1
           if ( right /= MPI_PROC_NULL ) then
              ! receive from right
              n = dez * (ex - sx + 1 + dsx + dex) * &
                        (ey - sy + 1 + dsy + dey) * ncomp
              allocate(buffrecv_zr(n))
              call mpi_irecv(buffrecv_zr, n, MPI_DOUBLE_PRECISION, right, 710,&
                   comm, exrq(nrq),ierr)
              else
                 exrq(nrq) = MPI_REQUEST_NULL
              endif
        endif

        if (cmode(3:3) == "1" ) then
           nrq = nrq + 1
           if ( left /= MPI_PROC_NULL) then
              ! send to left
              n =  dsz * (ex - sx + 1 + dsx + dex) * &
                         (ey - sy + 1 + dsy + dey) * ncomp
              allocate(buffsend_zl(n))
              buffsend_zl(:) = reshape(a(sx-dsx:ex+dex, sy-dsy:ey+dey, sz:sz+dsz-1, 1:ncomp), (/ n /))
              call mpi_isend(buffsend_zl, n, MPI_DOUBLE_PRECISION, left, 710,&
                   comm, exrq(nrq), ierr)
           else
              exrq(nrq) = MPI_REQUEST_NULL
           endif
        endif

        if (cmode(4:4) == "1")then
           nrq = nrq + 1
           if (right /= MPI_PROC_NULL) then
              ! send to right
              n =  dez * (ex - sx + 1 + dsx + dex) * &
                         (ey - sy + 1 + dsy + dey) * ncomp
              allocate(buffsend_zr(n))
              buffsend_zr(:) = reshape(a(sx-dsx:ex+dex, sy-dsy:ey+dey, ez-dez+1:ez, 1:ncomp), (/ n /))
              call mpi_isend(buffsend_zr, n, MPI_DOUBLE_PRECISION, right, 709,&
                   comm, exrq(nrq), ierr)
           else
              exrq(nrq) = MPI_REQUEST_NULL
           endif
        endif
        case default
           if (check_assertion(.true.))then
              call handle_error(DL_MG_ERR_UNSPECIFIED, &
                   msg="wrong dimension value in exchage_halo_begin_base")
           endif
        end select
#endif
   end subroutine exchange_halo_begin_base


   subroutine exchange_halo_end(dir,mode,thickness,mg,component)
     use dl_mg_types
     use dl_mg_errors
     implicit none

     integer, intent(in) :: dir, mode, thickness
     type(mg_t), target, intent(inout) :: mg
     character(len=1), intent(in) :: component

#ifdef MPI
     integer i, ierr
     integer asx, sx, aex, ex, asy,sy, aey, ey, asz, sz, aez, ez
     integer left, right
     integer ex_status(MPI_STATUS_SIZE,12)
     real(wp), pointer :: pp(:,:,:)
!#ifdef HAVE_CONTIGUOUS
!     contiguous pp
!#endif


     ! copy recv to arrays, deallocate auxiliary arrays

     call mg_get_idx(mg,sx,ex,sy,ey,sz,ez)

     !write(0,*) 'exch ', component, dir, nrq, exrq, MPI_REQUEST_NULL, mg%coords(3)
     !call barrier

     select case (dir)
     case(1)
        left = mg%neighbor(north)
        right = mg%neighbor(south)
     case(2)
        left = mg%neighbor(west)
        right = mg%neighbor(east)
     case(3)
        left = mg%neighbor(bottom)
        right = mg%neighbor(top)
     case default
        if (check_assertion(.true.)) then
           call handle_error(DL_MG_ERR_UNSPECIFIED,&
                msg = "wrong value for dir (direction) in  exchange_halo_end")
        endif
     end select

     select case (component)
     case ("p","z")
        asx = mg%psx
        aex = mg%pex
        asy = mg%psy
        aey = mg%pey
        asz = mg%psz
        aez = mg%pez
        select case(component)
        case("p")
           pp => mg%p
        case("z")
           pp => mg%z
        end select
        call mpi_waitall(nrq,exrq,ex_status,ierr)
        call exchange_halo_end_base(dir, left, right, mode, thickness,&
             asx, asy, asz, aex, aey, aez, 1, &
             sx, ex, sy, ey, sz, ez, pp)
     case ("w") ! nonlinear only
        asx = mg%psx
        aex = mg%pex
        asy = mg%psy
        aey = mg%pey
        asz = mg%psz
        aez = mg%pez
        call mpi_waitall(nrq,exrq,ex_status,ierr)
        call exchange_halo_end_base(dir, left, right, mode, thickness,&
             asx, asy, asz, aex, aey, aez, 1, &
             sx, ex, sy, ey, sz, ez, mg%w)
     case("f")
        asx = mg%fsx
        aex = mg%fex
        asy = mg%fsy
        aey = mg%fey
        asz = mg%fsz
        aez = mg%fez
        call mpi_waitall(nrq,exrq,ex_status,ierr)
        call exchange_halo_end_base(dir, left, right, mode, thickness,&
             asx, asy, asz, aex, aey, aez, 1,  &
             sx, ex, sy, ey, sz, ez, mg%f)
        case("c")
        asx = mg%csx
        aex = mg%cex
        asy = mg%csy
        aey = mg%cey
        asz = mg%csz
        aez = mg%cez
        call mpi_waitall(nrq,exrq,ex_status,ierr)
        call exchange_halo_end_base(dir, left, right, mode, thickness,&
             asx, asy, asz, aex, aey, aez, 3, &
             sx, ex, sy, ey, sz, ez, mg%c)

        case("r")
        asx = mg%rsx
        aex = mg%rex
        asy = mg%rsy
        aey = mg%rey
        asz = mg%rsz
        aez = mg%rez
        call mpi_waitall(nrq,exrq,ex_status,ierr)
        call exchange_halo_end_base(dir,  left, right, mode, thickness, &
             asx, asy, asz, aex, aey, aez, 1, &
             sx, ex, sy, ey, sz, ez, mg%r)
     case default
        if (check_assertion(.true.))then
           call handle_error(DL_MG_ERR_UNSPECIFIED, &
                msg="wrong mg component value in exchange_halo_end")
           endif
     end select
#endif
   end subroutine exchange_halo_end


   subroutine exchange_halo_end_base(dir, left, right, mode, thickness,&
        asx, asy, asz, aex, aey, aez, ncomp,&
        sx, ex, sy, ey, sz, ez, a)
     use dl_mg_types
     use dl_mg_errors
     implicit none
     integer, intent(in) ::dir, left, right, mode, thickness
     integer, intent(in) :: asx, asy, asz, aex, aey, aez, ncomp,&
          sx, ex, sy, ey, sz, ez
     real(wp), intent(inout) :: a(asx:aex, asy:aey, asz:aez,ncomp)

#ifdef MPI

     integer ierr, dsx, dex, dsy, dey, dsz, dez
     character(len=4) cmode

     write(cmode,'(I4)') mode

     select case (dir)
     case(1)
        dsx = thickness
        dex = thickness
        dsy = min(thickness, sy-asy)
        dey = min(thickness, aey-ey)
        dsz = min(thickness, sz-asz)
        dez = min(thickness, aez-ez)
        if (cmode(1:1) == "1" .and. left /= MPI_PROC_NULL) then
           ! transfer to left halo
           a(sx-dsx:sx-1,sy-dsy:ey+dey, sz-dsz:ez+dez, 1:ncomp) = &
                reshape(buffrecv_xl,(/ dsx, ey + dey - (sy - dsy) + 1,  ez + dez - (sz - dsz) + 1, ncomp /))
           deallocate(buffrecv_xl)
        endif

        if (cmode(2:2) == "1" .and. right /= MPI_PROC_NULL) then
         a(ex+1:ex+dex,sy-dsy:ey+dey, sz-dsz:ez+dez, 1:ncomp) = &
                reshape(buffrecv_xr,(/ dex, ey + dey - (sy - dsy) + 1,  ez + dez - (sz - dsz) + 1, ncomp /))
           deallocate(buffrecv_xr)
        endif

        if (cmode(3:3) == "1" .and. left /= MPI_PROC_NULL) then
            deallocate(buffsend_xl)
         endif
         if  (cmode(4:4) == "1" .and. right /= MPI_PROC_NULL ) then
            deallocate(buffsend_xr)
         endif
      case(2)
        dsx = min(thickness, sx-asx)
        dex = min(thickness, aex-ex)
        dsy = thickness
        dey = thickness
        dsz = min(thickness, sz-asz)
        dez = min(thickness, aez-ez)
        if (cmode(1:1) == "1" .and. left /= MPI_PROC_NULL) then
           ! transfer to left halo
           a(sx-dsx:ex+dex, sy-dsy:sy-1, sz-dsz:ez+dez, 1:ncomp) = &
                reshape(buffrecv_yl, (/ ex + dex - (sx - dsx) + 1, dsy,  ez + dez - (sz - dsz) + 1, ncomp /))
           deallocate(buffrecv_yl)
        endif

        if (cmode(2:2) == "1" .and. right /= MPI_PROC_NULL) then
           a(sx-dsx:ex+dex, ey+1:ey+dey, sz-dsz:ez+dez, 1:ncomp) = &
                reshape(buffrecv_yr,(/ ex + dex - (sx - dsx) + 1, dey,   ez + dez - (sz - dsz) + 1, ncomp /))
           deallocate(buffrecv_yr)
        endif

        if (cmode(3:3) == "1" .and. left /= MPI_PROC_NULL ) then
            deallocate(buffsend_yl)
         endif
         if  (cmode(4:4) == "1" .and. right /= MPI_PROC_NULL ) then
            deallocate(buffsend_yr)
         endif
       case(3)
        dsx = min(thickness, sx-asx)
        dex = min(thickness, aex-ex)
        dsy = min(thickness, sy-asy)
        dey = min(thickness, aey-ey)
        dsz = thickness
        dez = thickness
        if (cmode(1:1) == "1" .and. left /= MPI_PROC_NULL) then
           ! transfer to left halo
           a(sx-dsx:ex+dex, sy-dsy:ey+dey, sz-dsz:sz-1, 1:ncomp) = &
                reshape(buffrecv_zl, (/ ex + dex - (sx - dsx) + 1, ey +dey - (sy - dsy) + 1, dsz, ncomp /))
           deallocate(buffrecv_zl)
        endif

        if (cmode(2:2) == "1" .and. right /= MPI_PROC_NULL) then
           a(sx-dsx:ex+dex, sy-dsy:ey+dey, ez+1:ez+dez, 1:ncomp) = &
                reshape(buffrecv_zr, (/ ex + dex - (sx - dsx) + 1, ey + dey - (sy - dsy) + 1, dez, ncomp /))
           deallocate(buffrecv_zr)
        endif

        if (cmode(3:3) == "1" .and. left /= MPI_PROC_NULL ) then
            deallocate(buffsend_zl)
         endif
         if  (cmode(4:4) == "1" .and. right /= MPI_PROC_NULL) then
            deallocate(buffsend_zr)
         endif
      case default
         if(check_assertion(.true.))then
            call handle_error(DL_MG_ERR_UNSPECIFIED, &
                 msg = "wrong dimension value in exchage_halo_end_base")
         endif
     end select
#endif
   end subroutine exchange_halo_end_base


   !> useful for lapacian computation and other bits
  subroutine exchange_full_halos(mg, component)
    use dl_mg_types
    !use dl_mg_timer
    implicit none

    type(mg_t), intent(inout) :: mg
    character(len=1), intent(in) :: component

#ifdef MPI

      integer i

      !call mg_timer(start,trelax,t,tcomm)
      !$omp master
      do i = 1, 3
         call exchange_halo_begin(i,exch_full,1,mg,component(1:1))
         call exchange_halo_end(i,exch_full,1,mg,component(1:1))
      enddo
      !$omp end master
      !$omp barrier
      !call mg_timer(stop,trelax,t,tcomm)
#else
      call pbc_halo(mg, component(1:1))

#endif

    end subroutine exchange_full_halos

   
! if no MPI the halo data needs to be shifted for PBC
#ifndef MPI
   subroutine pbc_halo(mg, component)
     use dl_mg_types, only : mg_t, mg_get_idx
     use dl_mg_common_data, only : bc
     use dl_mg_errors
     implicit none
     type(mg_t), target, intent(inout) :: mg
     character(len=1) :: component

     integer i, ii, j, k, n, sx, ex, sy, ey, sz, ez
     real(wp), pointer :: a(:,:,:), b(:,:,:,:)
     logical found_pbc

     a => null(); b => null()

     found_pbc=.false.

     do i=1,3
        if (bc(i) == DL_MG_BC_PERIODIC) then
           found_pbc =.true.
           exit
        end if
     end do

     if ( .not. found_pbc) return

     n = 1

    call mg_get_idx(mg,sx,ex,sy,ey,sz,ez)
     select case(component)
     case("p")
        a => mg%p
     case("r")
        a => mg%r
     case("f")
        a => mg%f
     case("w")
        a => mg%w
     case("z")
        a => mg%z
     case("c")
        b => mg%c
        n = 3
     case default
        call handle_error(DL_MG_ERR_UNSPECIFIED, &
             msg="wrong array component in pbc_halo")
     end select

! Note: the halo needs to be exchange in full, i.e. ex-1 : ez+1, etc
!       for correct prolungation
!       below some data are transfered needlessly but is simple.
!       To be revised.

     do ii = 1, n

        if ( n == 3 ) a => remap (b(:,:,:,ii))

        !$OMP SINGLE
        if (bc(1) == DL_MG_BC_PERIODIC) then
           do k = sz-1, ez+1
              do j = sy-1, ey+1
                 a(sx-1,j,k) = a(ex,j,k)
                 a(ex+1,j,k) = a(sx,j,k)
              enddo
           enddo
        endif
        !$OMP END SINGLE NOWAIT

        !$OMP SINGLE
        if (bc(2) == DL_MG_BC_PERIODIC) then
           do k = sz-1, ez+1
              do i = sx-1, ex+1
                 a(i,sy-1,k) = a(i,ey,k)
                 a(i,ey+1,k) = a(i,sy,k)
              enddo
           enddo
        endif
        !$OMP END SINGLE NOWAIT

        !$OMP SINGLE
        if (bc(3) == DL_MG_BC_PERIODIC) then
           do j = sy-1, ey+1
              do i = sx-1, ex+1
                 a(i,j,sz-1) = a(i,j,ez)
                 a(i,j,ez+1) = a(i,j,sz)
              enddo
           enddo
        endif
        !$OMP END SINGLE NOWAIT
     end do
     !$omp barrier
     
     contains

       function remap(p)
         implicit none
         real(wp), target, intent(in) :: p(mg%csx:, mg%csy:, mg%csz:)
         real(wp), pointer :: remap(:,:,:)

         remap => p

       end function remap


   end subroutine pbc_halo
#endif


! useful for debugging
   subroutine minmax_grid(mg,comp,xmin,xmax)
     use dl_mg_types
     implicit none
     type(mg_t),intent(in) :: mg
     character(len=1), intent(in) :: comp
     real(wp), intent(out) :: xmin(3), xmax(3)

     integer i, ierr
     real(wp) ymin(3), ymax(3)

     select case (comp)
     case ("c")
        do i = 1, 3
           ymin(i) = minval(mg%c(mg%sx:mg%ex,mg%sy:mg%ey,mg%sz:mg%ez,i))
           ymax(i) = maxval(mg%c(mg%sx:mg%ex,mg%sy:mg%ey,mg%sz:mg%ez,i))
         enddo
#ifdef MPI
         call mpi_allreduce(ymin,xmin,3,mpi_double_precision,&
              MPI_MIN,mg%comm,ierr)
         call mpi_allreduce(ymax,xmax,3,mpi_double_precision,&
              MPI_MAX,mg%comm,ierr)
#else
      xmin = ymin; xmax = ymax
#endif
      case default
         write(0,*) ' minmax_grid not implemented for component', comp
      end select

   end subroutine minmax_grid


   subroutine barrier(comm)
     implicit none
     integer, optional, intent(in) :: comm

#ifdef MPI
     integer c, ierr

     if (present(comm)) then
        c=comm
     else
        c=MPI_COMM_WORLD
     endif

     call mpi_barrier(c,ierr)
#endif
   end subroutine barrier

#ifdef DUMP_ARRAYS
   subroutine dump_data(mg, component, fname)
     use dl_mg_types
     use mg_utils
     implicit none

     type(mg_t), target, intent(in) :: mg
     character(len=2), intent(in) :: component
     character(len=*), intent(in) :: fname


     integer comm, comm_x, comm_y, comm_z, gsx, gex, gsy, gey, gsz, gez, ierr
     integer gsizes(3), lbvec(3), ijkmin(3), ijkmax(3), dims(3)
     integer, allocatable :: bx(:,:), by(:,:), bz(:,:)
     real(wp), pointer ::vect(:,:,:) => null()

     comm = mg%comm
     ! get 3 1d subcommunicators
     call mpi_cart_sub(comm,(/ .true., .false., .false. /), comm_x, ierr)
     call mpi_cart_sub(comm,(/ .false., .true., .false. /), comm_y, ierr)
     call mpi_cart_sub(comm,(/ .false., .false., .true. /), comm_z, ierr)

     call get_mpi_grid(comm, dims=dims)
     allocate (bx(2,0:dims(1)-1), by(2,0:dims(2)-1), bz(2,0:dims(3)-1))
     call mpi_allgather((/ mg%sx, mg%ex /), 2,mpi_integer,bx,2,mpi_integer,comm_x,ierr)
     call mpi_allgather((/ mg%sy, mg%ey /), 2,mpi_integer,by,2,mpi_integer,comm_y,ierr)
     call mpi_allgather((/ mg%sz, mg%ez /), 2,mpi_integer,bz,2,mpi_integer,comm_z,ierr)

     gsx = minval(bx(1,:)); gex = maxval(bx(2,:))
     gsy = minval(by(1,:)); gey = maxval(by(2,:))
     gsz = minval(bz(1,:)); gez = maxval(bz(2,:))

     gsizes =(/ gex - gsx +1, gey - gsy +1, gez - gsz +1 /)
     ijkmin = (/ mg%sx, mg%sy, mg%sz /)
     ijkmax = (/ mg%ex, mg%ey, mg%ez /)

     select case (component)
     case ("p")
        vect => mg%p
        lbvec =(/ mg%psx, mg%psy, mg%psz /)
     case("r")
        vect => mg%r
        lbvec =(/ mg%rsx, mg%rsy, mg%rsz /)
     case( "c1")
        vect => mg%c(:,:,:,1)
        lbvec =(/ mg%csx, mg%csy, mg%csz /)
     case("c2")
        vect => mg%c(:,:,:,2)
        lbvec =(/ mg%csx, mg%csy, mg%csz /)
     case ("c3")
        vect => mg%c(:,:,:,3)
        lbvec =(/ mg%csx, mg%csy, mg%csz /)
     case default
        write(0,*) "dump vect not implemented for component :", component
        return
     end select

     call write_vect(vect, gsizes, lbvec, ijkmin, ijkmax, fname, comm)

   end subroutine dump_data
#endif
end module dl_mg_mpi
