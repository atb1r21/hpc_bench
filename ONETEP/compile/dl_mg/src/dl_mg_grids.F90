!> \brief contains subroutines that
!! i) distritute grid blocks to MPI ranks at all levels,
!! ii) compute the communication patterns
!!
!! Lucian Anton


module dl_mg_grids
  use dl_mg_types, only : mg_t, fd_t
  implicit none
  public

  ! the main data structure of the solver
  type(mg_t), allocatable, target :: mg(:)

  ! supplementary data structure for high order finite differences
  type(fd_t) :: fd

  ! flag to decide wheter to redistribute grid data for better load balance or
  ! lighter communication
  ! N.B : redistrubution is not implemented in this version
  logical, save :: redistribute_grid = .false.


  ! variables to assist data transfes between application and dl_mg grid
  integer, allocatable :: transfer_map_recv(:,:), transfer_map_send(:,:),&
       transfer_map_recv_coll(:,:), transfer_map_send_coll(:,:)

  ! auxiliary arrays of pointers
  type ptr_array_t
     integer, pointer :: p => null()
  end type ptr_array_t

  integer dbc_shift(3) ! useful to to distiguish between PBC and DBC

contains

  subroutine set_mg_levels(ierror)
    use dl_mg_common_data, only : mg_levels, full_agg_level, full_agg_size,&
         nx, ny, nz, nxc, nyc, nzc, npx, npy, npz, bc
    use dl_mg_params, only :  DL_MG_BC_PERIODIC, DL_MG_BC_DIRICHLET
    use dl_mg_errors
    implicit none

    integer, optional, intent(out) :: ierror

    integer i, j, nlev, np
    integer zpart(0:npz-1), mpart(0:npz-1)
    integer gridc(3), k(3), nxyzc(3)
    character(len=10) val


    do i=1,3
       if (bc(i) == DL_MG_BC_PERIODIC ) then
          dbc_shift(i) = 0
       elseif  (bc(i) == DL_MG_BC_DIRICHLET ) then
          ! Dirichlet BC keep also the BC layer on their grid
          ! hence the internal points start with a shift of 1
          ! at left boundary and -1 at right boundary
          dbc_shift(i) = 1
       else
          ! check_assertion call is used here for testing
          if (.not. check_assertion( bc(i) == DL_MG_BC_PERIODIC .or. &
               bc(i) == DL_MG_BC_DIRICHLET, DL_MG_ERR_UNKNOWN_BC)) then
             call handle_error(DL_MG_ERR_UNKNOWN_BC, ierror, &
                  "dl_mg_grids:set_mg_levels - wrong value of boundary tag")
             return
          endif
       endif
    enddo

    ! find the number of levels for the global grid
    gridc=(/nx, ny, nz/) + dbc_shift
    do i=1,3
       k(i)=1
       do while (mod(gridc(i),2) == 0 .and. gridc(i) > 2)
          gridc(i)=gridc(i)/2
          k(i)=k(i)+1
       enddo
    enddo

    mg_levels = minval(k)
    call get_environment_variable("DL_MG_MAX_LVL",val,j)
    if (j > 0) then
      read(val,*) i
      mg_levels = min(mg_levels, i)
    endif
    
    ! PBC needs an even number of grid points for red-black smoother
    ! hence we need to lower the number of level if any coarse grid
    ! size is odd in periodic directions
    gridc=(/nx, ny, nz/) + dbc_shift
    levels: do i = mg_levels, 1, -1
       do j = 1, 3
          if ( bc(j) == DL_MG_BC_PERIODIC .and. mod(gridc(j),2) == 1) then
             mg_levels = mg_levels -1
             exit levels
          endif
       enddo
       gridc(:) = gridc(:)/2
    enddo levels

    if (check_assertion(mg_levels == 0)) then
       call handle_error(DL_MG_ERR_UNSPECIFIED, &
            msg = "dl_mg_grids:set_mg_level - mg_level = 0 &
            & probable the number of grid points is odd in a periodic&
            & direction")
    endif

    gridc=(/nx, ny, nz/) + dbc_shift
    do i = mg_levels, 2, -1
       gridc(:) = gridc(:)/2
    enddo
    do i = 1,3
       if(bc(i) ==  DL_MG_BC_PERIODIC) then
          nxyzc(i) = gridc(i)
       elseif (bc(i) ==  DL_MG_BC_DIRICHLET) then
          ! internal number of points for Dirichlet BC
          nxyzc(i) = gridc(i)-1
       else
          call handle_error(DL_MG_ERR_UNSPECIFIED, &
               msg = "dl_mg_grids:set_mg_levels - wrong value of boundary tag")
       endif
    end do

    nxc = nxyzc(1); nyc = nxyzc(2); nzc = nxyzc(3)

#ifdef MPI

    if ( npx == 1 .and. npy == 1 .and. npz == 1 ) then
       full_agg_level = 0
       return
    endif

    ! however full aggregation level can be set higher if the
    ! global grid is small enough and  full_agg_size is set in dl_mg_init
    ! NB TO BE REVISED for MIXED BC
    if ( full_agg_size > 0 ) then
       gridc = (/nx, ny, nz/)
       do i = mg_levels,  2, -1
          if ( gridc(1) * gridc(2) * gridc(3)  <= full_agg_size ) then
             full_agg_level = max(i,2)
             exit
          endif
          gridc(:) = gridc(:)/2 + dbc_shift
       enddo
    endif
#endif

  end subroutine set_mg_levels


  subroutine set_mg_grids(ierror)
    use dl_mg_types
    use dl_mg_mpi
    use dl_mg_common_data, only : mg_levels, mg_comm, npz, full_agg_level, &
         isx, iex, isy, iey, isz, iez, nx, ny, nz, bc
    use dl_mg_errors
    use dl_mg_params, only : DL_MG_BC_DIRICHLET, DL_MG_BC_PERIODIC
    implicit none

    integer, optional, intent(inout) :: ierror

    integer coords(3), dims(3), iaux(3), myid
    integer  sx, ex, sy, ey, sz, ez, nlx, nly, nlz, left, right
    integer  d, j, t, deltahalo(6), sxyz(3), exyz(3), ierr
    type(ptr_array_t) :: sf(3), ef(3), sc(3), ec(3)

    if (present(ierror)) ierror = DL_MG_SUCCESS

    call set_mg_levels

    allocate(mg(mg_levels))

    ! we go twice over multigrid levels
    ! 1. in the first sweep over the levels compute the grids indices
    ! 2. sweep again the levels and build the communication pattern
    !   for prolongation at transition between inactive to active MPI ranks
    ! 3. aggregation ?

    ! start from the top level

    t = mg_levels

    ! all ranks are active at top level ( this might change)
    mg(t)%active = .true.
    mg(t)%level  = t

    call set_grid_comm(mg_comm, mg(t)%active, mg(t)%comm)

    call find_mpi_neighbors(mg(t))

    call get_mpi_grid(mg(t)%comm, coords=coords, dims=dims, rank=myid)
    mg(t)%coords = coords

    ! shorthands
    sxyz = (/ isx, isy, isz /)
    exyz = (/ iex, iey, iez /)

    deltahalo = (/ 1, 1, 1, 1, 1, 1 /) ! we stay with one layer halo for the time being
    ! for Dirichlet BC set s<> e<> from on internal points
    do d = 1, 3
       if ( coords(d) == 0 .and. bc(d) == DL_MG_BC_DIRICHLET) then
          sxyz(d) = sxyz(d) + 1 ! cover only inside grid points for DBC
          deltahalo(2*d-1) = 1
       endif
       if ( coords(d) == dims(d) -1 .and. bc(d) == DL_MG_BC_DIRICHLET) then
          exyz(d) =exyz(d) -1
          deltahalo(2*d) = 1
       endif
    enddo

    ! check for empty inner grids
    if ( .not. check_assertion(exyz(1) - sxyz(1) >= 0 .and. &
         exyz(2) - sxyz(2) >= 0 .and. &
         exyz(3) - sxyz(3) >= 0, DL_MG_ERR_GRID_EMPTY)) then
       call handle_error(DL_MG_ERR_GRID_EMPTY, ierror)
       return
    endif

    mg(t)%sx = sxyz(1)
    mg(t)%ex = exyz(1)
    mg(t)%sy = sxyz(2)
    mg(t)%ey = exyz(2)
    mg(t)%sz = sxyz(3)
    mg(t)%ez = exyz(3)

    !write(0,*) 'set grids', t, mg(t)%sx,  mg(t)%ex, mg(t)%sy, mg(t)%ey, mg(t)%sz,  mg(t)%ez

    call set_array_limits

    if ( t > 1 )  then
       mg(t-1)%active = .true.
       do d = 1,3
          if ( sxyz(d) == exyz(d) .and. mod(sxyz(d), 2) == 0) then
             mg(t-1)%active = .false.
          endif
       enddo
    endif

    if ( t == full_agg_level ) then
       mg(t)%aggregate = 1
       call build_aggregation_map
    else
       mg(t)%aggregate = 0
    end if

    !compute the index ranges for the coarse levels

    do t = mg_levels-1, 1, -1
       mg(t)%level = t
       ! get a communicator and topology for the active ranks at this level
       call set_grid_comm(mg(t+1)%comm, mg(t)%active, mg(t)%comm)

       ! propagate downwards the inactive site
       ! in this version once a rank is tagged as inactive it stays
       ! inactive at lower levels
       if ( .not. mg(t)%active ) then
          mg(t)%sx = 0
          mg(t)%ex = -1 ! -1 in order to have ex-sx+1=0
          mg(t)%sy = 0
          mg(t)%ey = -1
          mg(t)%sz = 0
          mg(t)%ez = -1
          if ( t > 1 ) then
             mg(t-1)%active=.false.
          endif
          cycle
       endif

       call find_mpi_neighbors(mg(t))

       call get_mpi_grid(mg(t)%comm,coords=coords,dims=dims,rank=myid)

       mg(t)%coords = coords

       sf(1)%p => mg(t+1)%sx
       sf(2)%p => mg(t+1)%sy
       sf(3)%p => mg(t+1)%sz
       ef(1)%p => mg(t+1)%ex
       ef(2)%p => mg(t+1)%ey
       ef(3)%p => mg(t+1)%ez

       sc(1)%p => mg(t)%sx
       sc(2)%p => mg(t)%sy
       sc(3)%p => mg(t)%sz
       ec(1)%p => mg(t)%ex
       ec(2)%p => mg(t)%ey
       ec(3)%p => mg(t)%ez

       if ( t < full_agg_level ) then
          ! there is only one rank holding the whole grid at this level
          iaux = (/ nx, ny, nz /) + dbc_shift
          do d = 1, 3
             sc(d)%p = 1 + dbc_shift(d)
             ec(d)%p = iaux(d)/2**(mg_levels - t)
          enddo
       else
          do d = 1, 3
             sc(d)%p = sf(d)%p / 2 + 1
             ec(d)%p = (ef(d)%p + mod(ef(d)%p, 2))/ 2
          enddo
       endif

       !write(0,*) 'set grids', t, mg(t)%active, mg(t)%sx,  mg(t)%ex, mg(t)%sy, mg(t)%ey, mg(t)%sz,  mg(t)%ez

       deltahalo = (/ 1, 1, 1, 1, 1, 1 /)
       do d = 1, 3
          if ( coords(d) == 0 ) then
             deltahalo(2*d-1) = 1
          endif
          if ( coords(d) == dims(d) -1 ) then
             deltahalo(2*d) = 1
          endif
       enddo

       call set_array_limits
       ! set activity at the next level
       if ( t > 1 )  then
          mg(t-1)%active = .true.
          sxyz(:) = (/ mg(t)%sx, mg(t)%sy, mg(t)%sz /)
          exyz(:) = (/ mg(t)%ex, mg(t)%ey, mg(t)%ez /)
          do d = 1,3
             if ( sxyz(d) == exyz(d) .and. mod(sxyz(d), 2) == 0) then
                mg(t-1)%active = .false.
                exit
             endif
          enddo
       endif

       if ( t == full_agg_level ) then
          mg(t)%aggregate = 1
          call build_aggregation_map
       else
          mg(t)%aggregate = 0
       end if

    end do

    call build_prolong_map(ierror)
    if (present(ierror)) then
      if (ierror /= DL_MG_SUCCESS) return
    endif

    ! communication needed for prologation on revived ranks is
    ! left apart for the time being

  contains

    subroutine set_array_limits
      implicit none

      ! set the limits for work arrays. N.B. includes halos or extra space for aggregation
      mg(t)%psx = mg(t)%sx - deltahalo(1)
      mg(t)%pex = mg(t)%ex + deltahalo(2)
      mg(t)%psy = mg(t)%sy - deltahalo(3)
      mg(t)%pey = mg(t)%ey + deltahalo(4)
      mg(t)%psz = mg(t)%sz - deltahalo(5)
      mg(t)%pez = mg(t)%ez + deltahalo(6)

      mg(t)%csx = mg(t)%sx - 2 !deltahalo(1)
      mg(t)%cex = mg(t)%ex + 2 !deltahalo(2)
      mg(t)%csy = mg(t)%sy - 2 !deltahalo(3)
      mg(t)%cey = mg(t)%ey + 2 !deltahalo(4)
      mg(t)%csz = mg(t)%sz - 2 !deltahalo(5)
      mg(t)%cez = mg(t)%ez + 2 !deltahalo(6)

      mg(t)%fsx = mg(t)%sx - 1
      mg(t)%fex = mg(t)%ex + 1
      mg(t)%fsy = mg(t)%sy - 1
      mg(t)%fey = mg(t)%ey + 1
      mg(t)%fsz = mg(t)%sz - 1
      mg(t)%fez = mg(t)%ez + 1

      mg(t)%rsx = mg(t)%sx - 1
      mg(t)%rex = mg(t)%ex + 1
      mg(t)%rsy = mg(t)%sy - 1
      mg(t)%rey = mg(t)%ey + 1
      mg(t)%rsz = mg(t)%sz - 1
      mg(t)%rez = mg(t)%ez + 1

    end subroutine set_array_limits


    subroutine build_aggregation_map
      use dl_mg_errors
      implicit none

      integer i, j, k, mastrank, myid, nproc, ierr
      logical on000

      on000 = mg(t)%coords(1) == 0 .and. &
           mg(t)%coords(2) == 0 .and. &
           mg(t)%coords(3) == 0

      if ( on000) then

         ! data aggregation relays on the assumption that
         ! rank 0 has (0,0,0) position in grid
         ! if this is not so stop the whole thing
         call get_mpi_grid(mg(t)%comm, rank=i)
         if ( .not. check_assertion(i == 0) ) then
            call handle_error(DL_MG_ERR_UNSPECIFIED, &
            msg = "aggregation map error: 000 grid rank is not rank 0 &
                 & in communicator; aggregation operation won't work")
         endif

         mg(t)%rsx = 1
         mg(t)%rex = (nx+1) / 2**(mg_levels - t) + 1
         mg(t)%rsy = 1
         mg(t)%rey = (ny+1) / 2**(mg_levels - t) + 1
         mg(t)%rsz = 1
         mg(t)%rez = (nz+1) / 2**(mg_levels - t) + 1

         mg(t)%agg_mast = .true.
         mg(t-1)%active = .true.
      else

         mg(t)%agg_mast = .false.
         if ( t > 1) then
            mg(t-1)%active = .false.
         endif

      endif
#ifdef MPI
      if ( on000) then
         call mpi_comm_size(mg(t)%comm, nproc, ierr)
         allocate(mg(t)%agg_map(7, 0:nproc-1))
      else
         allocate(mg(t)%agg_map(1,1))
      endif

      call mpi_cart_rank(mg(t)%comm, (/ 0, 0, 0 /), mastrank, ierr)
      if ( .not. on000) mg(t)%agg_map(1,1) = mastrank

      call mpi_comm_rank(mg(t)%comm, myid, ierr)
      call mpi_gather( (/ mg(t)%sx, mg(t)%ex, mg(t)%sy, mg(t)%ey, mg(t)%sz, mg(t)%ez, myid/), &
           7, mpi_integer, mg(t)%agg_map, 7, mpi_integer, mastrank, mg(t)%comm, ierr)
#endif
    end subroutine build_aggregation_map


    subroutine build_prolong_map(ierror)
      ! collects data about the active status of the 6 neighboring ranks
      ! set the type of activated ranks
      implicit none

      integer, optional, intent(inout) :: ierror
#ifdef MPI
      integer t, stage, d, i, n, dims(3), tag(6), myid, rq(12), thickness(3), ierr
      integer waitall_status(MPI_STATUS_SIZE,12)
      logical i_am_active, neigh_active(6), any_active, any_inactive
      logical, allocatable :: test_activity(:)

      if (present(ierror)) ierror = DL_MG_SUCCESS

      do t = max(full_agg_level+1, 2), mg_levels
         if (.not. mg(t)%active ) cycle
         allocate(mg(t)%prolong_transfers(6,3))
         mg(t)%prolong_transfers(:,:)=0
         i_am_active = mg(t-1)%active
         neigh_active = .false.
         call get_mpi_grid(mg(t)%comm, dims=dims, rank=myid)
         do stage = 1,3
            ! correcteness check
!!$                  thickness(1) = mg(t)%ex - mg(t)%sx
!!$                  thickness(2) = mg(t)%ey - mg(t)%sy
!!$                  thickness(3) = mg(t)%ez - mg(t)%sz
!!$                  n=0
!!$                  do i = 1, 3
!!$                     if ( thickness(i) == 0 ) n=n+1
!!$                  enddo
!!$                  select case (stage)
!!$                  case (1)
!!$                     if ( n == 0 ) then
!!$                        call error_abort("build_prolong_map: inactive site thicker &
!!$                             &than one layer in all dimensions at stage 1")
!!$                     endif
!!$                  case(2)
!!$                     if ( n < 1) then
!!$                        call error_abort("build_prolong_map: inactive site does not &
!!$                             &have thickness 1 in at least 2 dimensions at stage 2")
!!$                     endif
!!$                  case(3 )
!!$                     if  ( n < 3 ) then
!!$                        call error_abort("build_prolong_map: inactive site is not a &
!!$                             &unit cube at stage 3")
!!$                     endif
!!$                  end select

            do d = 1,3
               ! collect neighbors activity
               rq(:) = MPI_REQUEST_NULL
               do i = 1, 6
                  call mpi_irecv(neigh_active(i), 1, MPI_LOGICAL, mg(t)%neighbor(i), 950+i, &
                       mg(t)%comm, rq(i), ierr)
               enddo
               tag(:) = (/ 952, 951, 954, 953, 956, 955 /)
               do i = 6, 1, -1
                  call mpi_isend(i_am_active, 1, MPI_LOGICAL,mg(t)%neighbor(i),tag(i), &
                       mg(t)%comm, rq(7+(6-i)), ierr)
               enddo
               !write(0,*) 'build prolong map: before wait all ', t, stage, MPI_PROC_NULL, mg(t)%neighbor(1:6),  rq
               call mpi_waitall(12, rq, waitall_status, ierr)

               !write(0,*) 'build prolong map: after wait all ', t, stage, d, mg(t)%coords, neigh_active, i_am_active


               if ( i_am_active) then
                  ! send to right if not at right boundary (DBC)
                  if ( .not. neigh_active(2*d) .and. &
                       (mg(t)%neighbor(2*d) /= MPI_PROC_NULL)) then
                     mg(t)%prolong_transfers(2*d,stage) = 1
                  endif

                  ! for DBC one needs to send to left if a rank is
                  ! activated at left side
                  if ( (bc(d) == DL_MG_BC_DIRICHLET) .and. &
                       (.not. neigh_active(2*d-1)) .and. &
                       mg(t)%coords(d) == 1 ) then
                     mg(t)%prolong_transfers(2*d-1, stage) = 1
                  endif

               else

                  ! get data from left
                  if ( neigh_active(2*d-1) ) then
                     mg(t)%prolong_transfers(2*d-1, stage) = -1
                     !i_am_active = .true.
                  endif
                  !this should take care of DBC
                  if ( neigh_active(2*d) .and. &
                       (mg(t)%neighbor(2*d-1) == MPI_PROC_NULL))then
                     mg(t)%prolong_transfers(2*d,stage) = -1
                     !i_am_active = .true.
                  endif

               endif
               !write(t*1000+100*mg(t)%coords(3)+10*mg(t)%coords(2)+mg(t)%coords(1),*) &
               !     'prolong transfer ', t, stage, d,  i_am_active, mg(t)%prolong_transfers
            enddo
            if ( .not. i_am_active .and. minval(mg(t)%prolong_transfers(:,stage)) == -1 ) then
               i_am_active = .true.
            endif
         enddo

         if ( myid == 0) then
            allocate(test_activity(dims(1)*dims(2)*dims(3)))
         else
            allocate(test_activity(1))
         end if
         call mpi_gather(i_am_active, 1, MPI_LOGICAL, test_activity, 1, MPI_LOGICAL,&
              0, mg(t)%comm,ierr)
         if ( myid == 0) then
            do i=1, size(test_activity)
               if ( .not. check_assertion(test_activity(i), DL_MG_ERR_PROLONG_MAP)) then
                  call handle_error(DL_MG_ERR_PROLONG_MAP, ierror)
                  return
               endif
            enddo
         endif
         deallocate(test_activity)

         if ( maxval(mg(t)%prolong_transfers) == 0 .and. &
              minval(mg(t)%prolong_transfers) == 0 ) then
            deallocate(mg(t)%prolong_transfers)
         endif

         !if (allocated(mg(t)%prolong_transfers)) then
         !   write(0,*) 'build_prolong_map ', t, mg(t)%coords, mg(t)%prolong_transfers
         !endif

      enddo
#else
      if (present(ierror)) ierror = DL_MG_SUCCESS
#endif
    end subroutine build_prolong_map

  end subroutine set_mg_grids

  !> set the grid block size for openmp threads
  subroutine set_loop_blocks(min_thread_chunk_in, block_size_in)
    use dl_mg_common_data, only : mg_levels, blk, nthreads
    use dl_mg_types, only :  mg_get_idx
    use dl_mg_info, only : write_loop_block_info
    use dl_mg_errors
    implicit none

    integer, optional, intent(in) :: min_thread_chunk_in, block_size_in

    ! blocking parameters for loops over local grids
    integer  min_thread_chunk, block_size

    integer sx, ex, sy, ey, sz, ez, nx, ny, nz, nb, &
         sweep, i, j, k, t, nt, thchnk
    integer btemp, iy, iz, last_active_lvl, lval
    character(len=10) val

    allocate (blk(mg_levels))

    ! simple version that cuts the domain in blocks along y
    ! then each block is divided equally to the threads
    ! if the local domain is too small some threads will do nothing

    !
    ! need to study more the literature (Rivera, ...)
    ! to write a general parallel blocking

    ! minimum number of xz layers below which it does not make sense
    ! to use more than one thread. Default is set to 3.
     min_thread_chunk = 3
     if (present(min_thread_chunk_in)) then
        ! for usability the optional argument could be set <= 0
        ! to signal to use the default value
        if ( min_thread_chunk_in > 0) then
           min_thread_chunk = min_thread_chunk_in
        end if
    endif

    block_size = 0
    call get_environment_variable("DL_MG_BLK_SIZE",val,lval)
    if (lval > 0) then
       read(val,*) block_size
    else
       if (present(block_size_in) ) then
          ! for usability the optional argument could be set <= 0
          ! to signal to use the default value
          block_size = block_size_in
       end if
    end if

    do t = mg_levels, 1, -1
       if ( .not. mg(t)%active ) then
          cycle
       endif
       last_active_lvl = t

       call mg_get_idx(mg(t), sx, ex, sy, ey, sz, ez)
       nx = ex - sx + 1; ny = ey - sy + 1; nz = ez - sz + 1

       if (block_size <= 0) then
          block_size = ny
       else
          block_size = min(block_size, ny)
       end if
       

       ! compute a thread chunk
       if (  block_size <= min_thread_chunk * nthreads ) then
          ! use fewer threads if the domain is very small
          thchnk =  min_thread_chunk
       else
          ! if block_size  is not a multiple of nthreads go the next multiple
          ! and split that number to the number of threads
          ! in this case the last thread is lucky, less work to do
          thchnk = (block_size + min(mod(block_size, nthreads),1)*nthreads) / nthreads
       endif

       call compute_block_y(thchnk)

       ! block partition only in y plane so far
       blk(t)%dims(1) = nx
       blk(t)%nblks(1) = 1
       blk(t)%dims(3) = nz
       blk(t)%nblks(3) = 1



       ! lazy counting
       do sweep = 1, 2
          nb = 0
          do k = sz, ez, blk(t)%dims(3)
             do j = sy, ey, blk(t)%dims(2)
                do i = sx, ex, blk(t)%dims(1)
                   if ( sweep == 1) then
                      nb = nb + 1
                   else
                      nb = nb + 1
                      blk(t)%start(1, nb) = i
                      blk(t)%end  (1, nb)   = min(ex, i + blk(t)%dims(1) -1)
                      blk(t)%start(2, nb) = j
                      blk(t)%end  (2, nb)   = min(ey, j + blk(t)%dims(2) -1)
                      blk(t)%start(3, nb) = k
                      blk(t)%end  (3, nb)   = min(ez, k + blk(t)%dims(3) -1)
                   end if
                enddo
             enddo
          enddo
          if (sweep == 1) then
             allocate (blk(t)%start(3,nb), blk(t)%end(3,nb))
             blk(t)%nblks(2)=nb
          endif
       enddo
       !write(0,*) 'block ', t, sx,ex,sy,ey,sz,ez, blk(t)%dims, blk(t)%start
    end do

    call write_loop_block_info(blk, mg_levels, last_active_lvl, min_thread_chunk, block_size)

  contains

    !stub subroutine, to be expanded in more general versions
    subroutine compute_block_y(blks)
      implicit none

      integer, intent(in) :: blks

      blk(t)%dims(2) = max(1, blks)

      if ( check_assertion(blk(t)%dims(2) == 0)) then
         call handle_error(DL_MG_ERR_UNSPECIFIED, &
              msg = "set openmp blocks: block size is 0 in y direction")
      endif

    end subroutine compute_block_y

  end subroutine set_loop_blocks



!!$      subroutine compute_transfer_maps
!!$        use dl_mg_types
!!$        use dl_mg_mpi_header
!!$        use dl_mg_common_data, only : mg_levels, npz, isz, iez
!!$        implicit none
!!$
!!$#ifdef MPI
!!$        integer i, t, jrecv,jsend,so,eo,sm,em,sml,eml, sz, ez, myid, ierr
!!$        integer ts(3,npz),tr(3,npz), ot_zranges(2, 0:npz-1), mg_zranges(2,0:npz-1)
!!$
!!$          t=mg_levels
!!$          call mpi_allgather( (/isz, iez/), 2, mpi_integer, ot_zranges, 2, MPI_INTEGER, mg(t)%comm, ierr)
!!$          call mpi_allgather( (/mg(t)%sz, mg(t)%ez/), 2, mpi_integer, mg_zranges, 2, MPI_INTEGER, mg(t)%comm, ierr)
!!$          call mpi_comm_rank(mg(t)%comm, myid, ierr)
!!$
!!$          jrecv=0
!!$          jsend=0
!!$          so=ot_zranges(1,myid)
!!$          eo=ot_zranges(2,myid)
!!$
!!$          sz = mg(t)%sz; ez = mg(t)%ez
!!$
!!$          ! grab the boundary on z faces
!!$
!!$          do i = 0, npz-1
!!$             if (ot_zranges(2,i) >= sz-1 .and. ot_zranges(1,i) <= ez+1 ) then
!!$                jrecv=jrecv+1
!!$                tr(1,jrecv) = i
!!$                tr(2,jrecv) = max(ot_zranges(1,i),sz-1)
!!$                tr(3,jrecv) = min(ot_zranges(2,i),ez+1)
!!$             endif
!!$
!!$             if (mg_zranges(2,i) + 1 >= so .and. mg_zranges(1,i) - 1 <= eo ) then
!!$                jsend = jsend + 1
!!$                ts(1,jsend) = i
!!$                ts(2,jsend) = max(mg_zranges(1,i) - 1 ,so)
!!$                ts(3,jsend) = min(mg_zranges(2,i) + 1 ,eo)
!!$             endif
!!$
!!$          end do
!!$
!!$           allocate(transfer_map_recv(3,jrecv),transfer_map_send(3,jsend))
!!$
!!$           transfer_map_recv(:,1:jrecv)=tr(:,1:jrecv)
!!$           transfer_map_send(:,1:jsend)=ts(:,1:jsend)
!!$
!!$           !write(0,*) 'transfer_map_recv ', myid, transfer_map_recv
!!$           !write(0,*) 'transfer_map_send ', myid, transfer_map_send
!!$
!!$          ! now for the collection step, no halo of *_mg array needed
!!$
!!$          jrecv=0
!!$          jsend=0
!!$          so=ot_zranges(1,myid)
!!$          eo=ot_zranges(2,myid)
!!$
!!$          ! grab the boundary on z faces
!!$
!!$          sml=mg_zranges(1,myid)
!!$          eml=mg_zranges(2,myid)
!!$          if (myid == 0) sml=sml-1
!!$          if (myid == npz -1) eml=eml+1
!!$
!!$          do i = 0, npz - 1
!!$
!!$             sm=mg_zranges(1,i)
!!$             em=mg_zranges(2,i)
!!$
!!$             if (i == 0) sm = sm-1
!!$             if (i == npz - 1) em = em+1
!!$
!!$             if (em >= so .and. sm  <= eo ) then
!!$                jrecv=jrecv+1
!!$                tr(1,jrecv) = i
!!$                tr(2,jrecv) = max(sm ,so)
!!$                tr(3,jrecv) = min(em ,eo)
!!$             endif
!!$
!!$             if (ot_zranges(2,i) >= sml .and. ot_zranges(1,i) <= eml ) then
!!$                jsend=jsend+1
!!$                ts(1,jsend) = i
!!$                ts(2,jsend) = max(sml,ot_zranges(1,i))
!!$                ts(3,jsend) = min(eml,ot_zranges(2,i))
!!$             endif
!!$
!!$          end do
!!$
!!$           allocate(transfer_map_recv_coll(3,jrecv),transfer_map_send_coll(3,jsend))
!!$
!!$           transfer_map_recv_coll(:,1:jrecv)=tr(:,1:jrecv)
!!$           transfer_map_send_coll(:,1:jsend)=ts(:,1:jsend)
!!$
!!$           !write(0,*) 'transfer_map_recv_coll ', myid, transfer_map_recv_coll
!!$           !write(0,*) 'transfer_map_send_coll ', myid, transfer_map_send_coll
!!$
!!$#endif
!!$      end subroutine compute_transfer_maps

  !> Setup global instance of the higer order finite diference fd_t
  subroutine set_fd(fd_order, ierror)
    use dl_mg_common_data, only: mg_levels, isx, iex, isy, iey, isz, iez, &
         mx, my, mz, nxb, nyb, nzb, bc
    use dl_mg_errors, only: handle_error, DL_MG_ERR_ALLOCATION, &
         DL_MG_ERR_DEALLOCATION, DL_MG_ERR_MPI_FAIL, DL_MG_ERR_UNSPECIFIED
    use dl_mg_mpi_header
    use dl_mg_mpi
    use dl_mg_params, only : DL_MG_SUCCESS
    use dl_mg_defco_utils, only: dl_mg_defco_utils_assert, dl_mg_defco_utils_abort
    implicit none

    ! Arguments
    integer, intent(in) :: fd_order              !< finite diference order for derivatives
                                                 !! Cartesian directions (x=1,y=2,z=3)
    integer, optional, intent(inout) :: ierror   !< DL_MG error code to return

    ! Parameters
    character(len=*), parameter :: myself = "set_fd"
    integer, parameter          :: root_id = 0
    integer, parameter          :: left=0, right=1

    ! Local variables
    integer :: line_comm
    integer :: dims(3)
    integer :: coords(3)
    integer :: gstart, gend
    Integer :: nlayers
    integer :: isxyz(3), iexyz(3)
    integer :: d, dir, i, k, chnk, count, rsrc, rdest, ierr
    integer :: i1s, i1e, i2s, i2e, i3s, i3e
    integer, allocatable :: all_se(:,:), map_recv(:), &
                            buff_halos(:,:), gmap(:,:,:), halo(:,:)
    logical :: periodic(3), lg_aux(3)

    character(len=256) :: char_buffer

    if (present(ierror)) ierror = DL_MG_SUCCESS

    if (fd%initialised .and. fd%order == fd_order) return

    if (.not. (fd_order > 0 .and. fd_order <=12 .and. mod(fd_order,2)==0) ) then
       char_buffer = ""
       write(char_buffer,'(i32)') fd_order
       call dl_mg_defco_utils_assert(.false.,&
            'Unsupported FD order in '//myself//", order="//trim(char_buffer))
    end if

    ! JCW: Removed sanity checks for moving code into DL_MG (these scanned
    ! JCW: Removed sanity checks for moving code into DL_MG (these scanned
    ! JCW: arrays for occurences of NaN)
    ! TODO Consider replacing sanity checks with DL_MG native routines
    !      Santity checks performed for:
    !         fcn
    call dl_mg_defco_utils_assert(nxb >= fd_order + 1, &
         'The X dimension of the fine grid is smaller than the FD stencil')
    call dl_mg_defco_utils_assert(nyb >= fd_order + 1, &
         'The Y dimension of the fine grid is smaller than the FD stencil')
    call dl_mg_defco_utils_assert(nzb >= fd_order + 1, &
         'The Z dimension of the fine grid is smaller than the FD stencil')

    fd%order = fd_order

    ! Set grid geometry and indexing information in per-MPI rank fd_t instance
    fd%sx = isx; fd%sy = isy; fd%sz = isz
    fd%ex = iex; fd%ey = iey; fd%ez = iez
    fd%mx = mx;  fd%my = my;  fd%mz = mz

    isxyz(1) = isx
    isxyz(2) = isy
    isxyz(3) = isz
    iexyz(1) = iex
    iexyz(2) = iey
    iexyz(3) = iez

    fd%comm = mg(mg_levels)%comm
    ! Construct a per-MPI rank map of neighbouring MPI ranks along x, y and z
    ! directions which may contain halo elements for high order 1D finite
    ! differences

    ! For high-order 1D finite differences, each MPI rank needs to know the
    ! number of elements available in neighbouring ranks along the x, y and z
    ! directions up to the maximum halo size (or to the local grid boundary).


    ! main steps of the algorithm
    ! 1. loop through dimensions (axis)
    ! 2. find how many layers this ranks needs form each neighbor (left and right)
    ! 3. find the global indices on those neighbors, store them in fd together with
    !    mpi_rank (needed for recv operation)
    ! 3. find the send neighbors,  the amount of data and adresses (rank) to be send
    !      -this can done easily by aggregating the recv data on the rows
    !       of a global matrix and read the send information on columns


#ifdef MPI

    call get_mpi_grid(fd%comm, dims = dims, periods = periodic, coords = coords)

    DO d = 1, 3
       ! select the current axis and build a comm
       lg_aux(:) = .false.
       lg_aux(d) = .true.
       call mpi_cart_sub(fd%comm, lg_aux, line_comm, ierr)

       if (allocated(all_se)) deallocate(all_se)
       allocate(all_se(2, 0:dims(d)-1))
       call mpi_allgather((/ isxyz(d), iexyz(d)/), 2, mpi_integer, &
            all_se, 2, mpi_integer, line_comm, ierr)

       fd%gs(d) = minval(all_se(1,:))

       if (allocated(map_recv)) deallocate(map_recv)
       allocate(map_recv(0:dims(d)-1))
       if (allocated(halo)) deallocate(halo)
       allocate(halo(3,dims(d)))
       halo = -1
       if (allocated(gmap)) deallocate(gmap)
       if (allocated(buff_halos)) deallocate(buff_halos)
       allocate(gmap(2, 0:dims(d)-1, 0:dims(d)-1),  buff_halos(2, 0:dims(d) -1))
       !write(0,*) d, isxyz, iexyz, dims, all_se

       gstart = fd%gs(d)
       gend   = maxval(all_se(2,:))

       ! need to do left and right separately to avoid overlaps for PBC
       ! this would cost 3 more alltoall but it shouldn't be significant

       do dir = left, right
          ! compute the number of planes needed for each halo (left and right).
          ! Keep in mind that the
          ! source point  for derivarive computation shifts near boundary
          ! remember that the stencil has always fd_order + 1 points
          nlayers = fd_order/2
          if (.not. periodic(d)) then
             select case(dir)
             case(left)
                if ( isxyz(d) - gstart < fd_order/2 ) then
                   nlayers = isxyz(d) - gstart
                else if(gend - isxyz(d) < fd_order/2) then
                   ! close to the right grid boundary
                   ! left halo might need to be extended because there isn't
                   ! enough space for the right halo
                   nlayers = fd_order - (gend - isxyz(d))
                endif
             case(right)
                if (gend - iexyz(d) < fd_order/2) then
                   nlayers = gend -iexyz(d)
                else if (iexyz(d) - gstart < fd_order/2) then
                   nlayers = fd_order -(iexyz(d) - gstart)
                endif
             end select
          endif ! periodic bc
          ! find the coordinates of the ranks from which I need halo layers
          map_recv=0
          select case(dir)
          case(left)
             count = nlayers
             ! find how many neighbours are covered by halo to the left
             i = coords(d) -1
             do while (count > 0)
                if ( i < 0 ) i = dims(d) -1 ! for PBC
                chnk = min(all_se(2,i) - all_se(1,i) + 1, count)
                map_recv(i) = chnk
                count = count - chnk
                if ( count < 0) then
                   call handle_error(DL_MG_ERR_UNSPECIFIED, ierror, &
                        "dl_mg_grids:set_fd  count < 0 in halo allocation to the left")
                endif
                i = i -1
             end do
          case(right)
             count = nlayers
             i = coords(d) +1
             do while (count > 0)
                if ( i > dims(d) -1) i = 0 ! PBC
                chnk = min(all_se(2,i) - all_se(1,i) +1, count)
                map_recv(i) = chnk
                count = count - chnk
                if ( count < 0) then
                   call handle_error(DL_MG_ERR_UNSPECIFIED, ierror, &
                        "dl_mg_grids:set_fd  count < 0 in halo allocation to the right")
                endif
                i = i+1
             end do
          end select

          ! what we need for the halo exchange?
          ! for receive: source and the global indices of the received buffer
          ! for send   : destination and the global indices of the buffer to be send
          ! typically more than one sources or destination are expected.
          ! Find the number of sends/recv and the index range for each of them

          select case(dir)
          case(left)
             i = coords(d) -1
             k = 0

             do while (nlayers > 0)
                k=k+1
                if ( i < 0) i = dims(d) -1
                ! the end
                halo(2,k) = all_se(2,i)
                ! the start
                halo(1,k) = halo(2,k) - map_recv(i) +1
                ! the source
                call mpi_cart_shift(fd%comm, d-1, k, rsrc, rdest,  ierr)
                ! test rsrc is MPI_RANK_NULL
                halo(3, k) = rsrc
                nlayers = nlayers - map_recv(i)
                i = i-1
             enddo
             fd%num_neighbours_recv_l(d) = k
             if (allocated(fd%halo_map_recv_l(d)%m)) &
                  deallocate(fd%halo_map_recv_l(d)%m)
             allocate(fd%halo_map_recv_l(d)%m(3,k))
             fd%halo_map_recv_l(d)%m = halo(:,1:k)

          case(right)
             i = coords(d) +1
             k = 0
             halo = -1
             do while (nlayers > 0)
                k=k+1
                if ( i > dims(d) -1) i = 0
                ! the start
                halo(1, k) = all_se(1, i)
                ! the end
                halo(2,k) = halo(1,k) + map_recv(i) -1
                ! the source
                call mpi_cart_shift(fd%comm, d-1, -k, rsrc, rdest,  ierr)
                ! test rsrc is MPI_RANK_NULL
                halo(3, k) = rsrc
                nlayers = nlayers - map_recv(i)
                i = i+1
             enddo
             fd%num_neighbours_recv_r(d) = k
             if (allocated(fd%halo_map_recv_r(d)%m)) &
                  deallocate(fd%halo_map_recv_r(d)%m)
             allocate(fd%halo_map_recv_r(d)%m(3,k))
             fd%halo_map_recv_r(d)%m = halo(:,1:k)

          end select

          ! each rank builds a global map; will be used to derive the coordinates and the slices to send
          buff_halos(1,:) = 0
          buff_halos(2,:) = -1  !this so such that e - s + 1 = 0
          ! fill to the left
          select case (dir)
          case (left)
             i = coords(d) -1
             do k = 1, fd%num_neighbours_recv_l(d)
                if ( i < 0 ) i = dims(d) -1
                buff_halos(:, i) = fd%halo_map_recv_l(d)%m(1:2,k)
                i = i -1
             enddo
          case(right)
             ! fill to the right
             i = coords(d) +1
             do k = 1, fd%num_neighbours_recv_r(d)
                if ( i > dims(d) -1 ) i = 0
                buff_halos(:, i) = fd%halo_map_recv_r(d)%m(1:2,k)
                i = i +1
             enddo
          end select
          call mpi_allgather(buff_halos, 2 * dims(d), mpi_integer, &
               gmap, 2 * dims(d), mpi_integer, line_comm, ierr)

          select case(dir)
          case(left)
             ! send to right (lower diagonal part of gmap matrix)
             i = coords(d) +1
             k = 0
             halo = -1
             do
                if (i > dims(d) -1) then
                   if (periodic(d)) then
                      i = 0
                   else
                      exit
                   endif
                endif
                if (gmap(2, coords(d), i) - gmap(1, coords(d), i) + 1 == 0 &
                     .or. k == dims(d)) exit
                k=k+1
                ! the start
                halo(1,k) = gmap(1, coords(d), i)
                ! the end
                halo(2,k) = gmap(2, coords(d), i)
                ! the source
                call mpi_cart_shift(fd%comm, d-1, k, rsrc, rdest,  ierr)
                ! test rsrc is MPI_RANK_NULL
                halo(3, k) = rdest
                if ( k > dims(d)) then
                   call handle_error(DL_MG_ERR_UNSPECIFIED, ierror, &
                        "dl_mg_grids:set_fd  too may iterations in send to right map built")
                endif
                i = i+1
             enddo
             fd%num_neighbours_send_r(d) = k
             if (allocated(fd%halo_map_send_r(d)%m)) &
                  deallocate(fd%halo_map_recv_r(d)%m)
             allocate(fd%halo_map_send_r(d)%m(3,k))
             fd%halo_map_send_r(d)%m = halo(:,1:k)
          case(right)
             ! send to left (uper part of global matrix)
             i = coords(d) -1
             k = 0
             halo = -1
             do
                if (i < 0) then
                   if (periodic(d)) then
                      ! this will spin for ever if dim(d) == 1, hence the condition on k below
                      i = dims(d) -1
                   else
                      exit
                   endif
                endif
                if (gmap(2, coords(d), i) - gmap(1, coords(d), i) + 1 == 0 &
                     .or. k == dims(d)) exit
                k = k+1
                ! the start
                halo(1,k) = gmap(1,coords(d),i)
                ! the end
                halo(2,k) = gmap(2,coords(d),i)
                ! the source
                call mpi_cart_shift(fd%comm, d-1, -k, rsrc, rdest,  ierr)
                ! test rsrc is MPI_RANK_NULL
                halo(3, k) = rdest
                if ( k > dims(d)) then
                   call handle_error(DL_MG_ERR_UNSPECIFIED, ierror, &
                        "dl_mg_grids:set_fd too many iteration  in send to left map built")
                endif
                i = i-1
             enddo
             fd%num_neighbours_send_l(d) = k
             if (allocated(fd%halo_map_send_l(d)%m)) &
                  deallocate(fd%halo_map_send_l(d)%m)
             allocate(fd%halo_map_send_l(d)%m(3, k))
             fd%halo_map_send_l(d)%m = halo(:, 1:k)
          end select
       end do
       call mpi_comm_free(line_comm, ierr)

    enddo ! d = 1,3

!!$    call mpi_comm_rank(fd%comm, i, ierr)
!!$    write(200+i, '(4I4)') fd_order, coords
!!$    write(200+i, '(6I4)') isx, iex, isy, iey, isz, iez
!!$
!!$    do d = 1,3
!!$       write(200+i,*) '--------------------------'
!!$       write(200+i, '(200I4)') fd%num_neighbours_recv_l(d), fd%halo_map_recv_l(d)%m
!!$       write(200+i, '(200I4)') fd%num_neighbours_recv_r(d), fd%halo_map_recv_r(d)%m
!!$       write(200+i, '(200I4)') fd%num_neighbours_send_l(d), fd%halo_map_send_l(d)%m
!!$       write(200+i, '(200I4)') fd%num_neighbours_send_r(d), fd%halo_map_send_r(d)%m
!!$    enddo
!!$
!!$    call mpi_finalize(ierr)
!!$    stop

    ! find the halos ranges from fd maps
    do d = 1,3
       if (fd%num_neighbours_recv_l(d) > 0) then
          select case (d)
          case(1)
             if (bc(d) == DL_MG_BC_PERIODIC) then
                i1s = fd%sx - fd_order/2
                i1e = fd%sx - 1
             else
                i1s = fd%halo_map_recv_l(d)%m(1,fd%num_neighbours_recv_l(d))
                i1e = fd%halo_map_recv_l(d)%m(2,1)
             endif
             i2s = fd%sy
             i2e = fd%ey
             i3s = fd%sz
             i3e = fd%ez
          case(2)
             i1s = fd%sx
             i1e = fd%ex
             if (bc(d) == DL_MG_BC_PERIODIC) then
                i2s = fd%sy - fd_order/2
                i2e = fd%sy - 1
             else
                i2s = fd%halo_map_recv_l(d)%m(1,fd%num_neighbours_recv_l(d))
                i2e = fd%halo_map_recv_l(d)%m(2,1)
             endif
             i3s = fd%sz
             i3e = fd%ez
          case(3)
             i1s = fd%sx
             i1e = fd%ex
             i2s = fd%sy
             i2e = fd%ey
             if (bc(d) == DL_MG_BC_PERIODIC) then
                i3s = fd%sz - fd_order/2
                i3e = fd%sz -1
             else
                i3s = fd%halo_map_recv_l(d)%m(1,fd%num_neighbours_recv_l(d))
                i3e = fd%halo_map_recv_l(d) %m(2,1)
             end if
          end select

          fd%halo_indices(:,1,d) = (/ i1s, i1e, i2s, i2e, i3s, i3e /)
       else
          fd%halo_indices(:,1,d) = (/ 1, 0, 1, 0, 1, 0 /)
       endif

       if (fd%num_neighbours_recv_r(d) > 0) then
          select case(d)
          case(1)
             if (bc(d) == DL_MG_BC_PERIODIC) then
                i1s = fd%ex + 1
                i1e = fd%ex + fd_order/2
             else
                i1s = fd%halo_map_recv_r(d)%m(1,1)
                i1e = fd%halo_map_recv_r(d)%m(2,fd%num_neighbours_recv_r(d))
             endif
             i2s = fd%sy
             i2e = fd%ey
             i3s = fd%sz
             i3e = fd%ez
          case(2)
             i1s = fd%sx
             i1e = fd%ex
             if (bc(d) == DL_MG_BC_PERIODIC) then
                i2s = fd%ey + 1
                i2e = fd%ey + fd_order/2
             else
                i2s = fd%halo_map_recv_r(d)%m(1,1)
                i2e = fd%halo_map_recv_r(d)%m(2, fd%num_neighbours_recv_r(d))
             endif
             i3s = fd%sz
             i3e = fd%ez
          case(3)
             i1s = fd%sx
             i1e = fd%ex
             i2s = fd%sy
             i2e = fd%ey
             if (bc(d) == DL_MG_BC_PERIODIC) then
                i3s = fd%ez + 1
                i3e = fd%ez + fd_order/2
             else
                i3s = fd%halo_map_recv_r(d)%m(1, 1)
                i3e = fd%halo_map_recv_r(d)%m(2, fd%num_neighbours_recv_r(d))
             endif
          end select
          fd%halo_indices(:,2,d) = (/ i1s, i1e, i2s, i2e, i3s, i3e /)
       else
          fd%halo_indices(:,2,d) = (/ 1, 0, 1, 0, 1, 0 /)
       endif
    end do

#else
    ! No need to allocate fd%neighbour_map_* or halo_map_*, since code is
    ! compiled in serial and there is no parallel data decomposition.
    fd%num_neighbours_recv_l(:) = 0
    fd%num_neighbours_recv_r(:) = 0
    fd%num_neighbours_send_l(:) = 0
    fd%num_neighbours_send_r(:) = 0

#endif

    fd%initialised = .true.

  end subroutine set_fd

  !> reset fd data for a new calculation
  subroutine free_fd
    implicit none

    integer d

    fd%initialised = .false.
    do d=1,3
       if (allocated(fd%halo_map_recv_l(d)%m)) deallocate(fd%halo_map_recv_l(d)%m)
       if (allocated(fd%halo_map_recv_r(d)%m)) deallocate(fd%halo_map_recv_r(d)%m)
       if (allocated(fd%halo_map_send_l(d)%m)) deallocate(fd%halo_map_send_l(d)%m)
       if (allocated(fd%halo_map_send_r(d)%m)) deallocate(fd%halo_map_send_r(d)%m)
    enddo

  end subroutine free_fd

!!$  !> Writes the per-MPI-rank fd_t instance to out_unit in human-readable
!!$  !! format
!!$  subroutine debug_write_fd(out_unit)
!!$    use dl_mg_errors, only: handle_error, DL_MG_ERR_MPI_FAIL
!!$    use dl_mg_mpi_header
!!$    implicit none
!!$
!!$    ! Arguments
!!$    integer,intent(in) :: out_unit !< Per-MPI-rank output unit
!!$
!!$    ! Parameters
!!$    character(len=*), parameter :: myself = "debug_write_fd"
!!$
!!$    ! Local variables
!!$    integer :: my_id
!!$    integer :: num_ranks
!!$    integer :: ii, jj
!!$    integer :: ierr
!!$
!!$#ifdef MPI
!!$    call MPI_COMM_RANK(fd%comm, my_id, ierr)
!!$    if (ierr /= MPI_SUCCESS) then
!!$       call handle_error(DL_MG_ERR_MPI_FAIL, &
!!$            msg="Error in "//myself//": MPI_COMM_RANK failed.")
!!$    end if
!!$
!!$    call MPI_COMM_SIZE(fd%comm, num_ranks, ierr)
!!$    if (ierr /= MPI_SUCCESS) then
!!$       call handle_error(DL_MG_ERR_MPI_FAIL, &
!!$            msg="Error in "//myself//": MPI_COMM_SIZE failed.")
!!$    end if
!!$#else
!!$    my_id = 0
!!$    num_ranks = 1
!!$    ! In serial case, fd%num_neighbours_{l,r}(:) = 0 (Dirichlet BCs)
!!$#endif
!!$
!!$    write(out_unit,'(a,i5,a)') "fd output for ",num_ranks," MPI ranks"
!!$
!!$    ! Output contents of fd instance on each MPI rank in sequence
!!$    ! Output contents of fd instance on each MPI rank in sequence
!!$    write(out_unit,'(a5,i5,10x,3a5)') "Rank",my_id,"x","y","z"
!!$    write(out_unit,'(a20,3i5)') "Max halo thickness:",fd%max_halo_thickness(:)
!!$    write(out_unit,'(a20,3i5)') "Num. neighbours L:",fd%num_neighbours_l(:)
!!$    write(out_unit,'(a20,3i5)') "Num. neighbours R:",fd%num_neighbours_r(:)
!!$    write(out_unit,'(a20,3i5)') "Local grid start:", fd%sx, fd%sy, fd%sz
!!$    write(out_unit,'(a20,3i5)') "Local grid end:",   fd%ex, fd%ey, fd%ez
!!$    write(out_unit,'(a20,3i5)') "Local grid size:",  fd%mx, fd%my, fd%mz
!!$    do ii = 1, 3
!!$      write(out_unit,'(a,i3)') "L neighbours, dir.:",ii
!!$      if (fd%num_neighbours_l(ii) > 0) then
!!$        do jj = 1, size(fd%halo_map_l(ii)%m)
!!$          write(out_unit,'(2x,a,i3,a,i6,a,i6)') "Neighbour ",jj,":",&
!!$               fd%halo_map_l(ii)%m(jj)," elements, rank",fd%neighbour_map_l(ii)%m(jj)
!!$        end do
!!$      else
!!$        write(out_unit,'(2x,a)') "No neighbours"
!!$      end if
!!$      write(out_unit,'(a,i3)') "R neighbours, dir.:",ii
!!$      if (fd%num_neighbours_r(ii) > 0) then
!!$        do jj = 1, size(fd%halo_map_r(ii)%m)
!!$          write(out_unit,'(2x,a,i3,a,i6,a,i6)') "Neighbour ",jj,":",&
!!$               fd%halo_map_r(ii)%m(jj)," elements, rank",fd%neighbour_map_r(ii)%m(jj)
!!$        end do
!!$      else
!!$        write(out_unit,'(2x,a)') "No neighbours"
!!$      end if
!!$    end do
!!$
!!$    flush(out_unit)
!!$  end subroutine debug_write_fd
!!$
!!$  !> Un-setup global instance of fd_t (dealllocate arrays)
!!$  subroutine unset_fd(ierror)
!!$    use dl_mg_errors, only: handle_error, DL_MG_ERR_DEALLOCATION
!!$    implicit none
!!$
!!$    ! Arguments
!!$    integer, optional, intent(inout) :: ierror   !< DL_MG error code to return
!!$
!!$    ! Parameters
!!$    character(len=*), parameter :: myself = "unset_fd"
!!$
!!$    ! Local variables
!!$    integer :: icart
!!$    integer :: ierr
!!$
!!$    ! Deallocate map_array_t instances in per-MPI-rank fd instance
!!$    do icart = 1, 3
!!$      if ( fd%num_neighbours_l(icart) > 0 ) then
!!$        ! Only allocated if num_neighbours > 0
!!$        deallocate( fd%neighbour_map_l(icart)%m, stat=ierr )
!!$        if (ierr /= 0) then
!!$           call handle_error(DL_MG_ERR_DEALLOCATION, ierror=ierror, &
!!$           msg = myself//" fd%neighbour_map_l(icart)%m")
!!$        end if
!!$        deallocate( fd%halo_map_l(icart)%m, stat=ierr )
!!$        if (ierr /= 0) then
!!$           call handle_error(DL_MG_ERR_DEALLOCATION, ierror=ierror, &
!!$           msg = myself//" fd%halo_map_l(icart)%m")
!!$        end if
!!$      end if
!!$      if ( fd%num_neighbours_r(icart) > 0 ) then
!!$        deallocate( fd%neighbour_map_r(icart)%m, stat=ierr  )
!!$        if (ierr /= 0) then
!!$           call handle_error(DL_MG_ERR_DEALLOCATION, ierror=ierror, &
!!$           msg = myself//" fd%neighbour_map_r(icart)%m")
!!$        end if
!!$        deallocate( fd%halo_map_r(icart)%m, stat=ierr )
!!$        if (ierr /= 0) then
!!$           call handle_error(DL_MG_ERR_DEALLOCATION, ierror=ierror, &
!!$           msg = myself//" fd%halo_map_r(icart)%m")
!!$        end if
!!$      end if
!!$    end do
!!$
!!$  end subroutine unset_fd

end module dl_mg_grids
