!> \brief Transfers dielectric function (matrix coefficients) to the coarse grids
!!
!! Lucian Anton

module dl_mg_stencil_coefficients
  implicit none

contains


!!$  subroutine fine_stencil_from_fc_face(mg,fc)
!!$    use dl_mg_params
!!$    use dl_mg_types
!!$    use dl_mg_common_data, only : dx, dy, dz
!!$    implicit none
!!$    type(mg_t), intent(inout) :: mg
!!$    real(wp), intent(in)      :: fc(mg%sx-1:, mg%sy-1:, mg%sz-1:, :)
!!$
!!$    integer i,j,k
!!$
!!$    do k = mg%sz, mg%ez
!!$       do j = mg%sy, mg%ey
!!$          do i = mg%sx, mg%ex
!!$             mg%c(i,j,k,1) = fc(i-1,j,k,1) / (dx * dx)
!!$             mg%c(i,j,k,2) = fc(i  ,j,k,1) / (dx * dx)
!!$             mg%c(i,j,k,3) = fc(i,j-1,k,2) / (dy * dy)
!!$             mg%c(i,j,k,4) = fc(i  ,j,k,2) / (dy * dy)
!!$             mg%c(i,j,k,5) = fc(i,j,k-1,3) / (dz * dz)
!!$             mg%c(i,j,k,6) = fc(i,j,k  ,3) / (dz * dz)
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$  end subroutine fine_stencil_from_fc_face


  subroutine coarse_stencils(mg)
    use dl_mg_params
    use dl_mg_types
    use dl_mg_common_data, only : dx, dy, dz, bc, blk
    use dl_mg_errors
    implicit none

    type(mg_t), intent(inout) :: mg(:)

    integer sxc, exc, syc, eyc, szc, ezc, i, j, k, is, js, ks, t, if, jf, kf, &
         si, sj, sk, ds(3)
    real(wp), allocatable :: cofa(:,:,:,:)

    ! set index shifts between PBC and DBC
    do i=1,3
           if (bc(i) == DL_MG_BC_PERIODIC ) then
              ds(i) = 0
           elseif  (bc(i) == DL_MG_BC_DIRICHLET ) then
              ! Dirichlet BC keep also the BC layer on their grid
              ! hence the internal points start with a shift of 1
              ! at left boundary and -1 at right boundary
              ds(i) = 1
           else
              if (check_assertion(.true.)) then
                 call handle_error(DL_MG_ERR_UNSPECIFIED, &
                     msg = "dl_mg_stencil_coefficients:coarse_stencils &
                      &- wrong value of boundary tag")
              endif
           endif
        enddo

    ! rescale top level coefficients
    t = size(mg)

    call compute_blocked_loop(mg(t), blk(t))

!!$    !$OMP DO
!!$    do k = mg(t)%csz, mg(t)%cez
!!$       do j = mg(t)%csy, mg(t)%cey
!!$          do i = mg(t)%csx, mg(t)%cex
!!$             mg(t)%c(i,j,k,1) = mg(t)%c(i,j,k,1) / (dx * dx)
!!$             mg(t)%c(i,j,k,2) = mg(t)%c(i,j,k,2) / (dy * dy)
!!$             mg(t)%c(i,j,k,3) = mg(t)%c(i,j,k,3) / (dz * dz)
!!$          enddo
!!$       enddo
!!$    enddo
!!$    !$OMP ENDDO

    !$OMP MASTER
    do t = size(mg) - 1, 1, -1

       call mg_get_idx(mg(t), sxc, exc, syc, eyc, szc, ezc)

       !write(0,*) "stencil ", t, mg(t+1)%aggregate, sxc, exc

#ifdef MPI
       if (mg(t+1)%aggregate == 1) then

          !write(0,*) 'stencil  before aggregate', t, sxc, exc

          call aggregate_data(t)

          !write(0,*) 'stencil after aggregate', t, sxc, exc, mg(t)%active
          if (.not. mg(t)%active) exit

         do k=szc-ds(3),ezc
            kf=(2*k-1)
            do j=syc-ds(2),eyc
               jf=(2*j-1)
               do i=sxc-ds(1),exc
                  if=2*i-1

                   mg(t)%c(i,j,k,1) = 0.125_wp*(cofa(if,jf,kf,1) + cofa(if+1,jf,kf,1))
                   mg(t)%c(i,j,k,2) = 0.125_wp*(cofa(if,jf,kf,2) + cofa(if,jf+1,kf,2))
                   mg(t)%c(i,j,k,3) = 0.125_wp*(cofa(if,jf,kf,3) + cofa(if,jf,kf+1,3))

                enddo
             enddo
          enddo

          ! for PBC we need the left-1 layer populated with right layer
          ! becase there is no help from MPI
          if (bc(1) == DL_MG_BC_PERIODIC ) then
             mg(t)%c(sxc-1, syc:eyc, szc:ezc, :) = mg(t)%c(exc, syc:eyc, szc:ezc, :)
          endif
          if (bc(2) == DL_MG_BC_PERIODIC ) then
             mg(t)%c(sxc:exc, syc-1, szc:ezc, :) = mg(t)%c(sxc:exc, eyc, szc:ezc, :)
          endif
          if (bc(3) == DL_MG_BC_PERIODIC ) then
             mg(t)%c(sxc:exc, syc:eyc, szc-1, :) = mg(t)%c(sxc:exc, syc:eyc, ezc, :)
          endif

      else
#endif

         if (.not. mg(t+1)%active) exit
         call exchange_halo(t+1)
         if (.not. mg(t)%active) exit

         ! when not at boundary the value at s<d>c-1 is provided by the above exchange
         ! in this way one layer of halo memory is saved for c
         if ( mg(t)%coords(1) == 0 ) then
            si = sxc - ds(1)
         else
            si = sxc
         endif

         if ( mg(t)%coords(2) == 0 ) then
            sj = syc - ds(2)
         else
            sj = syc
         endif

         if ( mg(t)%coords(3) == 0 ) then
            sk = szc - ds(3)
         else
            sk = szc
         endif

         do k = sk, ezc
            kf = 2*k-1
            do j = sj, eyc
               jf = 2*j-1
               do i = si, exc
                  if = 2*i-1

                   mg(t)%c(i,j,k,1)= 0.125_wp*(mg(t+1)%c(if,jf,kf,1) + mg(t+1)%c(if+1,jf,kf,1))
                   mg(t)%c(i,j,k,2)= 0.125_wp*(mg(t+1)%c(if,jf,kf,2) + mg(t+1)%c(if,jf+1,kf,2))
                   mg(t)%c(i,j,k,3)= 0.125_wp*(mg(t+1)%c(if,jf,kf,3) + mg(t+1)%c(if,jf,kf+1,3))

                enddo
             enddo
          enddo
          ! we need to exchange halos in the cases from below
          ! becase next iteration does not take this branch
          if ( t == 1 .or. mg(t)%aggregate == 1) then
             call exchange_halo(t)
          endif
#ifdef MPI
       endif
#endif

    enddo
    !$OMP END MASTER

  contains

    subroutine compute_blocked_loop(mg, blk)
      implicit none

      type(mg_t), intent(inout) :: mg
      type(block_list_t), intent(in) :: blk

      integer s1, e1, s2, e2, s3, e3, i, j, k, iblk

      !$OMP DO SCHEDULE(STATIC,1)
      do iblk = 1, size(blk%start, dim=2)

         s1 = blk%start(1, iblk)
         e1 = min(s1 + blk%dims(1) - 1, mg%cex)
         ! blocks are defined for p grid, cover if necessary the c grid
         if (s1 == mg%sx) s1 = min( s1, mg%csx)

         s2 = blk%start(2, iblk)
         e2 =  min(s2 + blk%dims(2) - 1, mg%cey)
         if (s2 == mg%sy) s2 = min (s2, mg%csy)

         s3 = blk%start(3, iblk)
         e3 = min(s3 + blk%dims(3) - 1, mg%cez)
         if (s3 == mg%sz) s3 = min(s3, mg%csz )

         do k = s3, e3
            do j = s2, e2
               do i = s1, e1
                  mg%c(i,j,k,1) = mg%c(i,j,k,1) / (dx * dx)
                  mg%c(i,j,k,2) = mg%c(i,j,k,2) / (dy * dy)
                  mg%c(i,j,k,3) = mg%c(i,j,k,3) / (dz * dz)
               enddo
            enddo
         enddo
      enddo
      !$OMP ENDDO
    end subroutine compute_blocked_loop


    subroutine aggregate_data(t)
      use dl_mg_mpi
      use dl_mg_errors
      implicit none

      integer, intent(in) :: t

#ifdef MPI

      integer comm, coords(3), ierr
      integer i,j, slice, root, id, n
      integer asx, aex, asy, aey, asz, aez
      integer sx, ex, sy, ey, sz, ez
      integer, allocatable  :: recvcnts(:), recvdisp(:)
      real(wp), allocatable :: buff(:)

      if( .not.  mg(t+1)%active ) return

      comm = mg(t+1)%comm
      call mpi_cart_rank(comm,(/0, 0, 0/), root, ierr)
      call mg_get_idx(mg(t+1), sx, ex, sy, ey, sz, ez)
      call get_mpi_grid(comm, rank=id)

      ! mg(k+1)%agg_mast is true on (/ 0 0 0/)
      ! this is not very elegant at the moment
      if (mg(t+1)%agg_mast) then
!         write(0,*) ' agg mast ', t, mg(t+1)%coords, id
         if (allocated( cofa)) deallocate(cofa)
         allocate(cofa(mg(t+1)%rsx:mg(t+1)%rex, &
              mg(t+1)%rsy:mg(t+1)%rey, &
              mg(t+1)%rsz:mg(t+1)%rez, 3))

         call mpi_comm_size(comm,n,ierr)
         allocate(recvcnts(0:n-1),recvdisp(0:n-1))
!          write(0,*) ' agg mast ', n, size(mg(t+1)%agg_map)
         ! s<> -1 is for the left walls
         ! assumes that we aggregate on 0 rank of communicator
         ! save some memory
         recvcnts(0) = 0! 3 * (ex - (sx -1)  + 1) * ( ey - (sy - 1) + 1 ) &
                        ! * (ez - (sz -1) + 1)
         recvdisp(0) = 0
         do i = 1, n - 1
            call mpi_cart_coords(comm, i, 3, coords, ierr)
            asx = mg(t+1)%agg_map(1,i)
            if ( coords(1) == 0) asx = asx - 1
            aex = mg(t+1)%agg_map(2,i)
            asy = mg(t+1)%agg_map(3,i)
            if ( coords(2) == 0 ) asy = asy - 1
            aey = mg(t+1)%agg_map(4,i)
            asz = mg(t+1)%agg_map(5,i)
            if ( coords(3) == 0 ) asz = asz - 1
            aez = mg(t+1)%agg_map(6,i)
            recvcnts(i) = 3 * (aex - asx + 1) * ( aey - asy + 1) * (aez - asz + 1)
            recvdisp(i) = recvdisp(i-1) + recvcnts(i-1)
         enddo
         allocate(buff(sum(recvcnts)))
         call mpi_gatherv(MPI_IN_PLACE, 0, MPI_DOUBLE_PRECISION,&
              buff,recvcnts,recvdisp,MPI_DOUBLE_PRECISION,&
              root,comm,ierr)
         if ( .not. check_assertion(ierr == MPI_SUCCESS)) then
            call handle_error(DL_MG_ERR_UNSPECIFIED, &
                msg = 'error in stencil aggregation on root')
         endif
         ! this assumes that we use 0,0,0 rank for aggregation !
         cofa(sx-1:ex, sy-1:ey, sz-1:ez, 1:3) = mg(t+1)%c(sx-1:ex, sy-1:ey, sz-1:ez, 1:3)
         do i = 1, n - 1
            call mpi_cart_coords(comm, i, 3, coords,ierr)
 !           write(0,*) "agg map ", i, mg(t+1)%agg_map(1:6,i)
            asx = mg(t+1)%agg_map(1,i)
            if ( coords(1) == 0) asx = asx - 1
            aex = mg(t+1)%agg_map(2,i)
            asy = mg(t+1)%agg_map(3,i)
            if ( coords(2) == 0 ) asy = asy - 1
            aey = mg(t+1)%agg_map(4,i)
            asz = mg(t+1)%agg_map(5,i)
            if ( coords(3) == 0 ) asz = asz - 1
            aez = mg(t+1)%agg_map(6,i)
            slice = (aex - asx + 1) * (aey - asy + 1) * (aez - asz + 1)
            !do j = 1, 3
            cofa(asx:aex, asy:aey, asz:aez, 1:3) = &
                 reshape (buff(recvdisp(i) + 1 : recvdisp(i) + recvcnts(i)), &
                 (/ aex - asx + 1,  aey - asy + 1, aez - asz + 1, 3 /))
!            write(70+i, *) buff(recvdisp(i) + 1 : recvdisp(i) + recvcnts(i))
            !enddo
         enddo
      else
         root = mg(t+1)%agg_map(1,1)
         coords = mg(t+1)%coords
         asx = sx
         if (coords(1) == 0) asx = asx - 1
         asy = sy
         if (coords(2) == 0) asy = asy - 1
         asz = sz
         if (coords(3) == 0) asz = asz - 1
         n = 3 * (ex - asx + 1) * (ey - asy + 1) * (ez - asz + 1)
!         write(80+coords(3),*) mg(t+1)%c(asx:ex,asy:ey,asz:ez,3)
         call mpi_gatherv(mg(t+1)%c(asx:ex,asy:ey,asz:ez,1:3), n, MPI_DOUBLE_PRECISION,&
              i,(/0/),(/0/),MPI_DOUBLE_PRECISION,&
              root,comm,ierr)
         if ( .not. check_assertion(ierr == MPI_SUCCESS)) then
            call handle_error(DL_MG_ERR_UNSPECIFIED, &
                msg = 'error in stencil aggregation on sender')
         endif
      endif

#endif
    end subroutine aggregate_data


    subroutine exchange_halo(t)
      use dl_mg_mpi
      implicit none

      integer, intent(in) :: t
#ifdef MPI

      integer i

      do i =1, 3
         call exchange_halo_begin(i, exch_full, 2, mg(t),"c")
         call exchange_halo_end(i, exch_full, 2, mg(t),"c")
      enddo
#else
      call pbc_halo(mg(t),"c")
#endif
     end subroutine exchange_halo

  end subroutine coarse_stencils

end module dl_mg_stencil_coefficients
