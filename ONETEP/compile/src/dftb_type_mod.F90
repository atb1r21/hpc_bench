
! ab: quantities computed and reused in DFTB

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

module dftb_type

  use dense, only: DEM
  use neighbour_list, only: NL_NEIGHBOUR_LIST

  implicit none

  private

  public :: DFTB_STATE

  type DFTB_STATE
     type(NL_NEIGHBOUR_LIST)  :: nl
     type(DEM) :: dqdr(3)
     type(DEM) :: dcdr(3)
     type(DEM) :: dxdr(3)
  end type DFTB_STATE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module dftb_type
