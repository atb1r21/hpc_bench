! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by Chris-Kriton Skylaris
!
!   February 2000
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module geometry

  use constants, only: DP

  implicit none

  private

  type POINT
     real(kind=DP) :: X,Y,Z
  end type POINT


  interface operator(+)
     module procedure add_points
  end interface

  interface operator(-)
     module procedure subtract_points
  end interface

  interface operator(*)
     module procedure scale_point_real
     module procedure scale_point_int
  end interface

  interface operator(.CROSS.)
     module procedure cross_product
  end interface

  interface operator(.DOT.)
     module procedure dot_prod
  end interface


  public :: POINT

  public :: geometry_distance
  public :: local_displacement
  public :: magnitude
  public :: unit_vector
  public :: geometry_magnitude

  public :: operator(+), operator(-), operator(*)
  public :: operator(.CROSS.)
  public :: operator(.DOT.)
! jme: workaround for cray compiler bug 843178
#ifdef _CRAYFTN
  public :: add_points
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure elemental type (POINT) function unit_vector(A)

    implicit none

    type(POINT), intent(in) :: A

    unit_vector=(1.0_DP/(magnitude(A)))*A


  end function unit_vector


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  pure elemental function geometry_magnitude(A)

    use constants, only: DP
    implicit none

    real(kind=DP) :: geometry_magnitude

    type(POINT), intent(in) :: A

    geometry_magnitude=SQRT(A%X**2+A%Y**2+A%Z**2)

  end function geometry_magnitude

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  pure elemental function magnitude(A)

    use constants, only: DP
    implicit none

    real(kind=DP) :: magnitude

    type(POINT), intent(in) :: A

    magnitude=SQRT(A%X**2+A%Y**2+A%Z**2)

  end function magnitude


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  pure elemental function geometry_DISTANCE(point1,point2)

    use constants, only: DP
    implicit none

    real(kind=DP) :: geometry_distance

    type(POINT), intent(in) :: point1,point2


    geometry_DISTANCE=sqrt((point1%x-point2%x)**2 &
         + (point1%y-point2%y)**2 + (point1%z-point2%z)**2 )

  end function geometry_DISTANCE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  pure elemental type(POINT) function LOCAL_DISPLACEMENT(&
       A1_UNIT,A2_UNIT,A3_UNIT,LOCAL_1,LOCAL_2,LOCAL_3)

    implicit none

    type(POINT),   intent(in) :: A1_UNIT,A2_UNIT,A3_UNIT

    real(kind=DP), intent(in) :: LOCAL_1,LOCAL_2,LOCAL_3


    LOCAL_DISPLACEMENT=LOCAL_1*A1_UNIT+LOCAL_2*A2_UNIT+LOCAL_3*A3_UNIT



  end function LOCAL_DISPLACEMENT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  pure elemental type(POINT) function normal_unit_vector(A,B)

    implicit none

    type(POINT), intent(in) :: A,B
    type(POINT)             :: normal

    ! cks: normal_unit_vector = AxB/|AxB|

    normal%X=A%Y*B%Z-A%Z*B%Y
    normal%Y=A%Z*B%X-A%X*B%Z
    normal%Z=A%X*B%Y-A%Y*B%X


    normal_unit_vector=unit_vector(normal)


  end function normal_unit_vector


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! cks: for definition of new operators

  pure elemental type(POINT) function add_points(A,B)

    implicit none

    type(POINT), intent(in) :: A,B

    add_points%X=A%X+B%X
    add_points%Y=A%Y+B%Y
    add_points%Z=A%Z+B%Z

  end function add_points


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  pure elemental type(POINT) function subtract_points(A,B)

    implicit none

    type(POINT), intent(in) :: A,B

    subtract_points%X=A%X-B%X
    subtract_points%Y=A%Y-B%Y
    subtract_points%Z=A%Z-B%Z

  end function subtract_points


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  pure elemental type(POINT) function scale_point_real(scale,A)

    use constants, only: DP
    implicit none

    real(kind=DP), intent(in) :: scale

    type(POINT), intent(in) :: A

    scale_point_real%X=scale*A%X
    scale_point_real%Y=scale*A%Y
    scale_point_real%Z=scale*A%Z

  end function scale_point_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  pure elemental type(POINT) function scale_point_int(scale,A)

    use constants, only: DP
    implicit none

    integer, intent(in) :: scale

    type(POINT), intent(in) :: A

    scale_point_int%X=real(scale,kind=DP)*A%X
    scale_point_int%Y=real(scale,kind=DP)*A%Y
    scale_point_int%Z=real(scale,kind=DP)*A%Z

  end function scale_point_int


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  pure elemental type(POINT) function cross_product(A,B)

    implicit none

    type(POINT), intent(in) :: A,B

    cross_product%X=A%Y*B%Z-A%Z*B%Y
    cross_product%Y=A%Z*B%X-A%X*B%Z
    cross_product%Z=A%X*B%Y-A%Y*B%X

  end function cross_product


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure elemental real(kind=DP) function dot_prod(A,B)

    implicit none

    type(POINT), intent(in) :: A,B

    dot_prod=A%X*B%X+A%Y*B%Y+A%Z*B%Z

  end function dot_prod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module geometry












