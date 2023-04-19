  !>\brief collection of subroutines that allocate/deallocate
  !! memory for grid arrays
module dl_mg_alloc
  implicit none
  !> flags for the array types
  !! the vector type are practically identical becase all 3d array
  !! use for multgrid have a halo of size 1
  !! most probable this won't change in the forseable future
  integer, parameter :: p_vtype = 201, & !potetial array
       r_vtype = 202, & ! residual array
       f_vtype = 203, & ! source array
       c_vtype = 204, & ! permittivity type
       z_vtype = 999, & ! assume no halo, useful when an array,
                     ! is passed with explicit lower bonds
                     !  a(sx:, sy:, sz:)
                     ! could be useful in particular cases
                     ! Note: mg%w is p_vtype, dl_mg_types
       y_vtype = 666 ! grabs the boundary layer for Dirichlet BC
  interface mg_allocate
     module procedure mg_allocate_s, mg_allocate_a
  end interface mg_allocate
  
  interface mg_deallocate
     module procedure mg_deallocate_s, mg_deallocate_a
  end interface mg_deallocate
  
  private mg_allocate_s, mg_allocate_a, mg_deallocate_s, mg_deallocate_a
  
contains

  ! allocate mg arrays
      subroutine mg_allocate_s(mg, nonlinear, have_steric)
      use dl_mg_types, only : mg_t
      implicit none
      type(mg_t), intent(inout) :: mg
      logical, intent(in) :: nonlinear, have_steric

      if (.not. mg%active) return

      allocate(mg%p(mg%psx:mg%pex, mg%psy:mg%pey, mg%psz:mg%pez))
      allocate(mg%c(mg%csx:mg%cex, mg%csy:mg%cey, mg%csz:mg%cez,3))
      allocate(mg%f(mg%fsx:mg%fex, mg%fsy:mg%fey, mg%fsz:mg%fez))
      allocate(mg%r(mg%rsx:mg%rex, mg%rsy:mg%rey, mg%rsz:mg%rez))
      !> Attention: this subroutine allocates mg%w at top level were is
      !! not needed.  The aray version avoid this
      if (nonlinear) then
         allocate(mg%w(mg%psx:mg%pex, mg%psy:mg%pey, mg%psz:mg%pez))
      else
         ! w and d components will be passed to some kernels subroutines
         ! and have to be allocated acording to the standard
         allocate(mg%w(1:0, 1:0, 1:0))
      endif
      if (have_steric) then
         allocate(mg%d(mg%fsx:mg%fex, mg%fsy:mg%fey, mg%fsz:mg%fez))
      else
         allocate(mg%d(1:0, 1:0, 1:0))
      endif
    end subroutine mg_allocate_s


    subroutine mg_allocate_a(mg, nonlinear, have_steric)
      use dl_mg_types, only : mg_t
      implicit none
      type(mg_t), intent(inout) :: mg(:)
      logical, intent(in) :: nonlinear, have_steric

      integer i


      ! w array not need at top level in non-linear case
      do i = size(mg), 1, -1

         if (.not. mg(i)%active) cycle

         allocate(mg(i)%p(mg(i)%psx:mg(i)%pex, mg(i)%psy:mg(i)%pey, mg(i)%psz:mg(i)%pez))
         allocate(mg(i)%c(mg(i)%csx:mg(i)%cex, mg(i)%csy:mg(i)%cey, mg(i)%csz:mg(i)%cez,3))
         allocate(mg(i)%f(mg(i)%fsx:mg(i)%fex, mg(i)%fsy:mg(i)%fey, mg(i)%fsz:mg(i)%fez))
         allocate(mg(i)%r(mg(i)%rsx:mg(i)%rex, mg(i)%rsy:mg(i)%rey, mg(i)%rsz:mg(i)%rez))
         if (nonlinear) then
            allocate(mg(i)%w(mg(i)%psx:mg(i)%pex, mg(i)%psy:mg(i)%pey, mg(i)%psz:mg(i)%pez))
         else
            ! w and d components will be passed to some kernels subroutines
            ! and have to be allocated acording to the standard
            allocate(mg(i)%w(1:0, 1:0, 1:0))
         endif
         if ( have_steric) then
            allocate(mg(i)%d(mg(i)%fsx:mg(i)%fex, mg(i)%fsy:mg(i)%fey,&
                 mg(i)%fsz:mg(i)%fez))
         else
            allocate(mg(i)%d(1:0, 1:0, 1:0))
         endif
      enddo

    end subroutine mg_allocate_a


    subroutine mg_deallocate_s(mg)
      use dl_mg_types, only : mg_t
      implicit none
      type(mg_t) mg

      if(.not. mg%active) return

      !deallocate(mg%p, mg%c, mg%f, mg%r, mg%w, mg%d)
         if (allocated(mg%p)) deallocate(mg%p)
         if (allocated(mg%f)) deallocate(mg%f)
         if (allocated(mg%c)) deallocate(mg%c)
         if (allocated(mg%r)) deallocate(mg%r)
         if (allocated(mg%w)) deallocate(mg%w)
         if (allocated(mg%d)) deallocate(mg%d)
    end subroutine mg_deallocate_s


    subroutine mg_deallocate_a(mg)
      use dl_mg_types, only : mg_t
      implicit none
      type(mg_t) mg(:)

      integer i

      do i = size(mg), 1, -1
         if (.not. mg(i)%active) cycle
         !deallocate(mg(i)%p, mg(i)%c, mg(i)%f, mg(i)%r, &
         !     mg(i)%w, mg(i)%d)
         if (allocated(mg(i)%p)) deallocate(mg(i)%p)
         if (allocated(mg(i)%f)) deallocate(mg(i)%f)
         if (allocated(mg(i)%c)) deallocate(mg(i)%c)
         if (allocated(mg(i)%r)) deallocate(mg(i)%r)
         if (allocated(mg(i)%w)) deallocate(mg(i)%w)
         if (allocated(mg(i)%d)) deallocate(mg(i)%d)
      enddo

    end subroutine mg_deallocate_a


    !> allocates 3d real array with halos according to the type
    !! an global indices
    !! the type is irrelevant becase all type have a halo of size 1
    !! most probable this won't change in the forseable future
    subroutine galloc(mg, a, type, ierr)
      use dl_mg_params, only : wp, DL_MG_ERR_UNSPECIFIED
      use dl_mg_types, only : mg_t
      use dl_mg_common_data, only : isx, isy, isz, &
           iex, iey, iez
      use dl_mg_errors
      implicit none

      type(mg_t), intent(inout) :: mg
      real(wp), allocatable, intent(out) :: a(:,:,:)
      integer, optional, intent(in) :: type ! p,f,r,c,z
      integer,optional, intent(out) :: ierr
      
      integer :: t
      
      
      if (present(type)) then
         t=type
      else
         t=p_vtype
      end if
      
      select case (t)
      case(p_vtype)
         allocate(a(mg%psx:mg%pex, mg%psy:mg%pey, mg%psz:mg%pez))
      case(r_vtype)
         allocate(a(mg%rsx:mg%rex, mg%rsy:mg%rey, mg%rsz:mg%rez))
      case(f_vtype)
         allocate(a(mg%fsx:mg%fex, mg%fsy:mg%fey, mg%fsz:mg%fez))
      case(c_vtype)
         allocate(a(mg%csx:mg%cex, mg%csy:mg%cey, mg%csz:mg%cez))
      case(z_vtype)
         allocate(a(isx:iex,isy:iey,isz:iez))
      case default
         if (check_assertion(.true.) ) then
            call handle_error(DL_MG_ERR_UNSPECIFIED, ierr, &
                 msg="galloc: wrong array type")
         endif
      end select
      
    end subroutine galloc


    ! for symmetry, one can use dealloc as well
    subroutine gdealloc(a, ierr)
      use dl_mg_params, only : wp, DL_MG_SUCCESS
      implicit none
      real(wp), allocatable, intent(inout) :: a(:,:,:)
      integer, optional, intent(out) :: ierr

      if(present(ierr)) then
         ierr = DL_MG_SUCCESS
      end if
      
      if (allocated(a)) then
         deallocate(a)
      else
         continue
         ! error here ?
      end if
      
    end subroutine gdealloc
    
  
  subroutine get_lbound(mg, vtype, si)
    use dl_mg_types, only : mg_t
    use dl_mg_common_data, only : isx, isy, isz
    use dl_mg_mpi, only : get_mpi_grid
    use dl_mg_errors
    implicit none
    type(mg_t), intent(in) :: mg
    integer, intent(in)    :: vtype
    integer, intent(out)   :: si(3)

    integer coords(3)
    
    select case(vtype)
    case(p_vtype)
       si(1) = mg%psx
       si(2) = mg%psy
       si(3) = mg%psz
    case(r_vtype)
       si(1) = mg%rsx
       si(2) = mg%rsy
       si(3) = mg%rsz
    case(f_vtype)
       si(1) = mg%fsx
       si(2) = mg%fsy
       si(3) = mg%fsz
    case(c_vtype)
       si(1) = mg%csx
       si(2) = mg%csy
       si(3) = mg%csz
    case(z_vtype)
       si(1) = mg%sx
       si(2) = mg%sy
       si(3) = mg%sz
    case(y_vtype)
       ! grab the boundary layer for DBC
       ! useful to hadndle the input arrays
       call get_mpi_grid(mg%comm, coords=coords)
       si(1) = mg%sx
       si(2) = mg%sy
       si(3) = mg%sz
       if (coords(1) == 0) then
          si(1) = isx
       endif
       if (coords(2) == 0) then
          si(2) = isy
       endif
       if (coords(3) == 0) then
          si(3) = isz
       endif
    case default
       call handle_error(DL_MG_ERR_UNSPECIFIED, &
            msg = "get_lbound: unkown array type in get_lbound")
    end select
  end subroutine get_lbound

end module dl_mg_alloc
