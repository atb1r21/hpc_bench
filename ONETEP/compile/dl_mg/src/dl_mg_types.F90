!
!> \brief derived data types that encapsulate information on grid setup, to be used at each multigrid level
!!
!! For MG calculation we need the following bits of information
!!
!! start/end in local domain in terms of global index space -> sx,ex,sy,ye,sz,ez
!! local grid start/end points  must allow for halo or extra space needed when
!! restricting the number of MPI ranks sxc,exc,syc,eyc,szc,ezc  ( cg(3,2) )
!! data arrays p(sxc:exc,syc:eyc,szc:ezc), cof(:,:,:,8), f(...), r(...)
!! boundary arrays ( save memory when BC=0 with a flag)
!! communicator ( MPI_COMM_NULL if core is not active)
!! should we have one more flag to signal active, inactive rank and a given MG level ?
!! MPI coords
!! myid
!! global grid boundaries to take care of :
!!  N, S, W, E, B, T, NW, NE, NB, NT, SW, SE, SB, ST, WB, WT, EB,
!!  ET, NWB, NEB, NWT, NET, SWB, SEB, SWT, SET
!! flag to mark that restriction at this level needs aggregation, i. e.,
!! the next coarse level has a reduced number of ranks, hence some ranks will receive data
!! from the discarded ranks.
!! timer storage usefull for performace tunning
!! computation, communication, threads(?)
!
!
!
!! Lucian Anton February 2013
!

module dl_mg_types
  use dl_mg_params
  implicit none

  !> boundaries codes
  integer, parameter :: north =  1,&
                        south =  2,&
                        west  =  3,&
                        east  =  4,&
                        bottom=  5,&
                        top   =  6,&
                        nw    =  7,&
                        ne    =  8,&
                        nb    =  9,&
                        nt    = 10,&
                        sw    = 11,&
                        se    = 12,&
                        sb    = 13,&
                        st    = 14,&
                        wb    = 15,&
                        wt    = 16,&
                        eb    = 17,&
                        et    = 18,&
                        nwb   = 19,&
                        neb   = 20,&
                        nwt   = 21,&
                        net   = 22,&
                        swb   = 23,&
                        seb   = 24,&
                        swt   = 25,&
                        set   = 26

  type boundary_values_t
     real(wp), allocatable, dimension(:,:) :: n,s,w,e,b,t
     real(wp), allocatable, dimension(:)   :: nw, ne, nb, nt, sw, se, &
          sb, st, wb, wt, eb, et
     real(wp) c(nwb:set) ! corners
  end type boundary_values_t

  
  type mg_t
     integer sx, ex, sy, ey, sz, ez       !< global indices range of the local grid on this rank
     !> active boundary list
     integer active_bd(26)
     type(boundary_values_t) bd

     integer comm        !< MPI communicator
     integer coords(3)   !< MPI coords

     integer neighbor(26) !< nearest (touching) neighbors in MPI topology

     integer level        !< multigrid level

     logical active      !< wheter this raks is active at this level

!> map for data communicating at prolongations stage
!! points : newly activated rank needs data from the halo of an active neighbor
!!          in three this is done in three stages
!! dimensions of prolong_trnasfers are (6,3)
!! for each stage a rank can send data (1), receive data (-1) or do nothing (0)
!! the first dimension is for the neighborghs n,s,w,e,b,t; second dimension id for the
!! stage
!! if this array is not allocated it means that no prolongation transfers are necessary on this rank

    integer, allocatable :: prolong_transfers(:,:)

     integer aggregate   !< aggregate before reduce MPI ranks
                         !! convention   1 z direction
                         !!             10 y direction
                         !!            100 x direction
                         !!            110 xy direction
                         !!            ...
     logical agg_mast !< if true receives data in aggregation
     integer, allocatable :: agg_map(:,:) !< holds the range of data to be aggregated
                                          !! and the rank which has the data before
                                          !! aggregation in the following order:
                                          !! {start x, endx, start y, end y, start z, end z, rank} ...


     integer psx, pex, psy, pey, psz, pez !< indices for the extent of compute grid and local array bounds
                                          !! these indices could cover more that logical one
                                          !! e.g. in case of aggregation

     real(wp), allocatable :: p(:,:,:)   !< solution field
     real(wp), allocatable :: w(:,:,:)   !< correction field, used only in nonlinear case
     integer csx, cex, csy, cey, csz, cez
     real(wp), allocatable :: c(:,:,:,:) !< permittivity coefficents
     integer fsx, fex, fsy, fey, fsz, fez
     real(wp), allocatable :: f(:,:,:)   !< source term
     real(wp), allocatable :: d(:,:,:)   !< steric potential, used for PBE
     integer rsx, rex, rsy, rey, rsz, rez
     real(wp), allocatable :: r(:,:,:)   !< residual
     real(wp), pointer :: z(:,:,:) ! used to hook local arrays and
                                   ! pass them to the communication routines (e.g. full hallo exchange)
  end type mg_t


  !> Used in fd_t to allow variable length arrays which map distance between neighbouring
  !! MPI ranks to rank id / halo elements available
  !! dim (2,n_of_neighbourd)
  !! for fixed neighbour first element contains the rank of the neighbor,
  !! the second the size of the buffers to be transffered
  !! i.e. a%m(1,1) -> is the rank of the nearest neighbour
  !!      a%m(2,1) -> the buffer size of nearest neighbour
  type map_array1_t
     integer, allocatable :: m(:)
  end type map_array1_t
  type map_array2_t
     integer, allocatable :: m(:,:)
  end type map_array2_t

  !> Type containing data for high-order finite differences, used in defect
  !! correction
  type fd_t
     logical :: initialised
     integer :: sx, ex, sy, ey, sz, ez       !< global indices range of the local grid on this rank
     integer :: gs(3)                        !< global start indices, needed for PBC
     integer :: mx, my, mz                   !< per-MPI rank fine grid size
     integer :: comm                         !< MPI communicator for high-order FD
     integer :: order                        !< approximation order for finite difference formula
                                             !! direction (x=1,y=2,z=3).
     integer :: num_neighbours_recv_l(3)     !< Number of neighbours in the "left" direction, or
     !! direction of decreasing global indexes in x, y, z for
     integer :: num_neighbours_send_l(3)
                                             !! each Cartesian direction.
     integer :: num_neighbours_recv_r(3)          !< Number of neighbours in the "right" direction, or
     !! direction of increasing global indexes in x, y, z for
     integer :: num_neighbours_send_r(3)
                                             !! each Cartesian direction.
     type(map_array2_t) :: halo_map_recv_l(3), &
          halo_map_recv_r(3), &
          halo_map_send_l(3), &
          halo_map_send_r(3) !< Contains the  global index ranges (1,2) and
                             !! rank id (3) of the halo
                             !! elements to requestor send  from/to each rank, along
                             !! each Cartesian direction (x=1, y=2, z=3).
                             !! Array index for the m(:) component of the type
                             !! refers to neighbouring MPI rank
                             !! distance (nearest is 1, increasing value is
                             !! increasing distance) from this rank along
                             !! Cartesian direction.
     integer :: halo_indices(6,2,3) !< halo's indices in the following order
                                    !! start-end (xs,xe,ys,ye,zs,ze),
                                    !! left-right halo,
                                    !!dimension
                                    !! x ->1, y->2, z->3
  end type fd_t


  type block_list_t
     integer dims(3)
     integer nblks(3)
     integer, allocatable :: start(:,:)
     integer, allocatable :: end(:,:)
  end type block_list_t

  contains


    subroutine mg_get_idx(mg, sx, ex, sy, ey, sz, ez)
      implicit none
      type(mg_t), intent(in) :: mg
      integer, intent(out)   :: sx, ex, sy, ey, sz, ez

      sx = mg%sx; ex = mg%ex
      sy = mg%sy; ey = mg%ey
      sz = mg%sz; ez = mg%ez

    end subroutine mg_get_idx

    
end module dl_mg_types
