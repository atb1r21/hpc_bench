! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   The subroutines in this file were written by:                             !
!                                                                             !
!   Jose M Escartin, Andrea Greco, and Nicholas D.M. Hine.                    !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   This module was spun off in August 2016 from previously-existing          !
!   types and procedures within basis_mod.F90, by JM Escartin.                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module datatypes

  use constants, only: DP
  implicit none

  private

!------------------------------------------------------------------------------

  ! Data types.

  type COEF
     real(kind=DP)    :: d
     complex(kind=DP) :: z
     logical          :: iscmplx
  end type COEF

  ! ndmh_pointerfun
  type FUNCTIONS
#ifndef POINTERFUN
     real(kind=DP), allocatable    :: d(:)
     complex(kind=DP), allocatable :: z(:)
#else
     real(kind=DP), pointer    :: d(:) => NULL()
     complex(kind=DP), pointer :: z(:) => NULL()
#endif
     logical          :: iscmplx
  end type FUNCTIONS

  type FFTBOX_DATA
#ifndef POINTERFUN
     real(kind=DP), allocatable    :: d(:,:,:)
     complex(kind=DP), allocatable :: z(:,:,:)
#else
     real(kind=DP), pointer    :: d(:,:,:) => NULL()
     complex(kind=DP), pointer :: z(:,:,:) => NULL()
#endif
     logical          :: iscmplx
  end type FFTBOX_DATA

!------------------------------------------------------------------------------

  ! List of public types and procedures.

  public :: COEF, FUNCTIONS, FFTBOX_DATA

  public :: data_check_alloc
  public :: data_check_unalloc
  public :: data_set_to_zero
  public :: data_size

  public :: data_coef_abs
  !public :: data_coef_reduce
  public :: data_coef_scale
  public :: data_coef_assign

  public :: data_functions_alloc
  public :: data_functions_dealloc
  public :: data_functions_dot
  public :: data_functions_copy
  public :: data_functions_scale
  public :: data_functions_axpy
  public :: data_functions_mix
  public :: data_functions_abs

  public :: data_fftbox_alloc
  public :: data_fftbox_dealloc
  public :: data_fftbox_scale
  public :: data_fftbox_axpy
  public :: data_fftbox_copy

  ! ndmh_pointerfun
  !public :: operator(+), operator(-), operator(*), operator(/)

  ! Interfaces
  interface data_check_alloc
     module procedure data_functions_check_alloc
     module procedure data_fftbox_check_alloc
  end interface data_check_alloc

  interface data_check_unalloc
     module procedure data_functions_check_unalloc
     module procedure data_fftbox_check_unalloc
  end interface data_check_unalloc

  interface data_functions_axpy
     module procedure data_functions_axpy_real
     module procedure data_functions_axpy_complex
     module procedure data_functions_axpy_coef
  end interface data_functions_axpy

  interface data_set_to_zero
     module procedure data_set_coef_to_zero
     module procedure data_set_functions_to_zero
     module procedure data_set_fftboxdata_to_zero
  end interface data_set_to_zero

  interface data_size
     module procedure data_size_functions
     module procedure data_size_fftboxdata
  end interface data_size

  interface data_functions_mix
     module procedure data_functions_mix_real
  end interface data_functions_mix

  interface data_coef_scale
     module procedure data_coef_scale_real
     module procedure data_coef_scale_complex
  end interface data_coef_scale

  interface data_functions_scale
     module procedure data_functions_scale_real
     module procedure data_functions_scale_complex
     module procedure data_functions_scale_coef
  end interface data_functions_scale

  interface data_fftbox_scale
     module procedure data_fftbox_scale_real
     module procedure data_fftbox_scale_complex
     module procedure data_fftbox_scale_coef
  end interface data_fftbox_scale

  interface data_coef_assign
     module procedure  data_coef_assign_real
     module procedure  data_coef_assign_complex
  end interface data_coef_assign

  interface data_fftbox_axpy
     module procedure data_fftbox_axpy_real
  end interface data_fftbox_axpy

  interface data_functions_dot
     module procedure data_functions_dot
     module procedure data_functions_dot_array
  end interface data_functions_dot

! ndmh_pointerfun
!  interface operator(+)
!     module procedure data_coef_sum
!     module procedure data_functions_sum
!  end interface

!  interface operator(-)
!     module procedure data_coef_diff
!     module procedure data_functions_diff
!  end interface

!  interface operator(*)
!     module procedure data_product_real_coef
!     module procedure data_product_complex_coef
!     module procedure data_product_coef
!     module procedure data_product_real_functions
!     module procedure data_product_complex_functions
!     module procedure data_product_coef_functions
!     module procedure data_product_functions
!  end interface

!  interface operator(/)
!     module procedure data_division_coef
!  end interface

  ! Auxiliary string
#ifdef POINTERFUN
  character(len=2), parameter :: ptr_str = '=>'
#else
  character(len=0), parameter :: ptr_str = ''
#endif

contains


  subroutine data_functions_alloc(funcs,fsize,iscmplx)

    use utils, only: utils_alloc_check, utils_assert

    implicit none

    ! Arguments
    type(FUNCTIONS), intent(inout) :: funcs
    integer, intent(in) :: fsize
    logical, intent(in) :: iscmplx

    ! Locals
    integer :: ierr

    funcs%iscmplx = iscmplx

    call utils_assert(data_functions_check_unalloc(funcs), &
       'ERROR in data_functions_alloc: array component already allocated.')

    if (funcs%iscmplx) then
       allocate(funcs%z(fsize),stat=ierr)
       call utils_alloc_check('data_functions_alloc',ptr_str//'funcs%z',ierr)
    else
       allocate(funcs%d(fsize),stat=ierr)
       call utils_alloc_check('data_functions_alloc',ptr_str//'funcs%d',ierr)
    end if

  end subroutine data_functions_alloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine data_functions_dealloc(funcs)

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNCTIONS), intent(inout) :: funcs

    ! Locals
    integer :: ierr

    if (funcs%iscmplx) then
       deallocate(funcs%z,stat=ierr)
       call utils_dealloc_check('data_functions_dealloc',ptr_str//'funcs%z',ierr)
#ifdef POINTERFUN
       nullify(funcs%z)
#endif
    else
       deallocate(funcs%d,stat=ierr)
       call utils_dealloc_check('data_functions_dealloc',ptr_str//'funcs%d',ierr)
#ifdef POINTERFUN
       nullify(funcs%d)
#endif
    end if

  end subroutine data_functions_dealloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  logical function data_functions_check_alloc(funcs) result(check_alloc)
    ! jme: returns true if the appropriate component of func is allocated

    implicit none

    ! Arguments
    type(FUNCTIONS), intent(in) :: funcs

#ifndef POINTERFUN
    if (funcs%iscmplx) then
       check_alloc = allocated(funcs%z)
    else
       check_alloc = allocated(funcs%d)
    end if
#else
    if (funcs%iscmplx) then
       check_alloc = associated(funcs%z)
    else
       check_alloc = associated(funcs%d)
    end if
#endif

  end function data_functions_check_alloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  logical function data_functions_check_unalloc(funcs) result(check_unalloc)
    ! jme: returns true if no component of func is allocated

    implicit none

    ! Arguments
    type(FUNCTIONS), intent(in) :: funcs

#ifndef POINTERFUN
    check_unalloc = .not. (allocated(funcs%z) .or. allocated(funcs%d))
#else
    check_unalloc = .not. (associated(funcs%z) .or. associated(funcs%d))
#endif

  end function data_functions_check_unalloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine data_fftbox_alloc(fftbox,n1,n2,n3,iscmplx)

    use utils, only: utils_alloc_check, utils_assert

    implicit none

    ! Arguments
    type(FFTBOX_DATA), intent(inout) :: fftbox
    integer, intent(in) :: n1,n2,n3
    logical, intent(in) :: iscmplx

    ! Locals
    integer :: ierr

    fftbox%iscmplx = iscmplx

    call utils_assert(data_fftbox_check_unalloc(fftbox), &
       'ERROR in data_fftbox_alloc: array component already allocated.')

    if (fftbox%iscmplx) then
       allocate(fftbox%z(n1,n2,n3),stat=ierr)
       call utils_alloc_check('data_fftbox_alloc',ptr_str//'fftbox%z',ierr)
    else
       allocate(fftbox%d(n1,n2,n3),stat=ierr)
       call utils_alloc_check('data_fftbox_alloc',ptr_str//'fftbox%d',ierr)
    end if

  end subroutine data_fftbox_alloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine data_fftbox_dealloc(fftbox)

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(FFTBOX_DATA), intent(inout) :: fftbox

    ! Locals
    integer :: ierr

    if (fftbox%iscmplx) then
       deallocate(fftbox%z,stat=ierr)
       call utils_dealloc_check('data_fftbox_dealloc',ptr_str//'fftbox%z',ierr)
#ifdef POINTERFUN
       nullify(fftbox%z)
#endif
    else
       deallocate(fftbox%d,stat=ierr)
       call utils_dealloc_check('data_fftbox_dealloc',ptr_str//'fftbox%d',ierr)
#ifdef POINTERFUN
       nullify(fftbox%d)
#endif
    end if

  end subroutine data_fftbox_dealloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  logical function data_fftbox_check_alloc(fftbd) result(check_alloc)
    ! jme: returns true if the appropriate component of fftbd is allocated

    implicit none

    ! Arguments
    type(FFTBOX_DATA), intent(in) :: fftbd

#ifndef POINTERFUN
    if (fftbd%iscmplx) then
       check_alloc = allocated(fftbd%z)
    else
       check_alloc = allocated(fftbd%d)
    end if
#else
    if (fftbd%iscmplx) then
       check_alloc = associated(fftbd%z)
    else
       check_alloc = associated(fftbd%d)
    end if
#endif

  end function data_fftbox_check_alloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  logical function data_fftbox_check_unalloc(fftbd) result(check_unalloc)
    ! jme: returns true if no component of fftbd is allocated

    implicit none

    ! Arguments
    type(FFTBOX_DATA), intent(in) :: fftbd

#ifndef POINTERFUN
    check_unalloc = .not. (allocated(fftbd%z) .or. allocated(fftbd%d))
#else
    check_unalloc = .not. (associated(fftbd%z) .or. associated(fftbd%d))
#endif

  end function data_fftbox_check_unalloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  type(COEF) function data_functions_dot(bras,kets,startbra,startket,length)

    use linalg, only: linalg_ddot, linalg_zdotc
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(FUNCTIONS), intent(in) :: bras
    type(FUNCTIONS), intent(in) :: kets
    integer, intent(in), optional :: startbra
    integer, intent(in), optional :: startket
    integer, intent(in), optional :: length

    ! Local variables
    integer :: start1, start2, end1, end2, totlen1, totlen2

    ! agreco: dot product between two real or two complex functions
    call utils_assert(bras%iscmplx.eqv.kets%iscmplx,'Error in &
       &data_functions_dot: incompatible argument types')

    totlen1 = data_size_functions(bras)
    totlen2 = data_size_functions(kets)

    if (present(startbra)) then
       start1 = startbra
       if (present(length)) then
          end1 = start1 + length - 1
       else
          end1 = totlen1
       end if
    else
       start1 = 1
       end1 = totlen1
    end if

    if (present(startket)) then
       start2 = startket
       if (present(length)) then
          end2 = start2 + length - 1
       else
          end2 = totlen2
       end if
    else
       start2 = 1
       end2 = totlen2
    end if

    ! use zdotc to take conjugate of first vector if complex functions
    if (bras%iscmplx) then
       data_functions_dot%iscmplx = .true.
       data_functions_dot%z = linalg_zdotc(end1-start1+1, &
          bras%z(start1:end1),1,kets%z(start2:end2),1)
    else
       data_functions_dot%iscmplx = .false.
       data_functions_dot%d = linalg_ddot(end1-start1+1, &
          bras%d(start1:end1),1,kets%d(start2:end2),1)
    end if

  end function data_functions_dot



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine data_functions_copy(yfun,xfun,starty,startx,length)

    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(FUNCTIONS), intent(inout) :: yfun
    type(FUNCTIONS), intent(in) :: xfun
    integer, intent(in), optional :: starty
    integer, intent(in), optional :: startx
    integer, intent(in), optional :: length

    ! Local variables
    integer :: start1, start2, end1, end2, totlen1, totlen2

    call utils_assert(yfun%iscmplx.eqv.xfun%iscmplx,'Error in &
       &data_functions_copy: incompatible argument types')

    totlen1 = data_size_functions(yfun)
    totlen2 = data_size_functions(xfun)

    if (present(starty)) then
       start1 = starty
       if (present(length)) then
          end1 = start1 + length - 1
       else
          end1 = totlen1
       end if
    else
       start1 = 1
       end1 = totlen1
    end if

    if (present(startx)) then
       start2 = startx
       if (present(length)) then
          end2 = start2 + length - 1
       else
          end2 = totlen2
       end if
    else
       start2 = 1
       end2 = totlen2
    end if


    if (yfun%iscmplx) then
       yfun%z(start1:end1) = xfun%z(start2:end2)
    else
       yfun%d(start1:end1) = xfun%d(start2:end2)
    end if

  end subroutine data_functions_copy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine data_functions_axpy_real(yfun,xfun,alpha,starty,startx,length)

    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(FUNCTIONS), intent(inout) :: yfun
    type(FUNCTIONS), intent(in) :: xfun
    real(kind=DP) :: alpha
    integer, intent(in), optional :: starty
    integer, intent(in), optional :: startx
    integer, intent(in), optional :: length

    ! Local variables
    integer :: start1, start2, end1, end2, totlen1, totlen2

    call utils_assert(yfun%iscmplx.eqv.xfun%iscmplx,'Error in &
       &data_functions_axpy_real: incompatible argument types')

    totlen1 = data_size_functions(yfun)
    totlen2 = data_size_functions(xfun)

    if (present(starty)) then
       start1 = starty
       if (present(length)) then
          end1 = start1 + length - 1
       else
          end1 = totlen1
       end if
    else
       start1 = 1
       end1 = totlen1
    end if

    if (present(startx)) then
       start2 = startx
       if (present(length)) then
          end2 = start2 + length - 1
       else
          end2 = totlen2
       end if
    else
       start2 = 1
       end2 = totlen2
    end if


    if (yfun%iscmplx) then
       yfun%z(start1:end1) = yfun%z(start1:end1) + alpha * xfun%z(start2:end2)
    else
       yfun%d(start1:end1) = yfun%d(start1:end1) + alpha * xfun%d(start2:end2)
    end if

  end subroutine data_functions_axpy_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine data_functions_axpy_complex(yfun,xfun,alpha,starty,startx,length)

    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(FUNCTIONS), intent(inout) :: yfun
    type(FUNCTIONS), intent(in) :: xfun
    complex(kind=DP) :: alpha
    integer, intent(in), optional :: starty
    integer, intent(in), optional :: startx
    integer, intent(in), optional :: length

    ! Local variables
    integer :: start1, start2, end1, end2, totlen1, totlen2

    !jmecmplx: is this a problem?  I would just check that yfun is indeed complex
    call utils_assert(yfun%iscmplx.eqv.xfun%iscmplx,'Error in &
       &data_functions_axpy_complex: incompatible argument types')

    totlen1 = data_size_functions(yfun)
    totlen2 = data_size_functions(xfun)

    if (present(starty)) then
       start1 = starty
       if (present(length)) then
          end1 = start1 + length - 1
       else
          end1 = totlen1
       end if
    else
       start1 = 1
       end1 = totlen1
    end if

    if (present(startx)) then
       start2 = startx
       if (present(length)) then
          end2 = start2 + length - 1
       else
          end2 = totlen2
       end if
    else
       start2 = 1
       end2 = totlen2
    end if

    if (yfun%iscmplx) then
       yfun%z(start1:end1) = yfun%z(start1:end1) + alpha * xfun%z(start2:end2)
    else
       ! jd: Is this branch needed in this subroutine at all?
       call utils_assert(.false., 'data_functions_axpy_complex: &
            &Attempt to axpy two real functions with a complex alpha.')
!       yfun%d(start1:end1) = yfun%d(start1:end1) + real(alpha) * xfun%d(start2:end2)
    end if

  end subroutine data_functions_axpy_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine data_functions_axpy_coef(yfun,xfun,alpha,starty,startx,length)
    implicit none

    ! Arguments
    type(FUNCTIONS), intent(inout) :: yfun
    type(FUNCTIONS), intent(in)    :: xfun
    type(COEF), intent(in)         :: alpha
    integer, intent(in), optional  :: starty
    integer, intent(in), optional  :: startx
    integer, intent(in), optional  :: length

    if (alpha%iscmplx) then
       call data_functions_axpy_complex(yfun, xfun, alpha%z, starty, startx, length)
    else
       call data_functions_axpy_real(yfun, xfun, alpha%d, starty, startx, length)
    end if

  end subroutine data_functions_axpy_coef



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine data_functions_mix_real(yfun,xfun,alpha,starty,startx,length)

    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(FUNCTIONS), intent(inout) :: yfun
    type(FUNCTIONS), intent(in) :: xfun
    real(kind=DP) :: alpha
    integer, intent(in), optional :: starty
    integer, intent(in), optional :: startx
    integer, intent(in), optional :: length

    ! Local variables
    integer :: start1, start2, end1, end2, totlen1, totlen2

    call utils_assert(yfun%iscmplx.eqv.xfun%iscmplx,'Error in &
       &data_functions_mix_real: incompatible argument types.')

    totlen1 = data_size_functions(yfun)
    totlen2 = data_size_functions(xfun)

    if (present(starty)) then
       start1 = starty
       if (present(length)) then
          end1 = start1 + length - 1
       else
          end1 = totlen1
       end if
    else
       start1 = 1
       end1 = totlen1
    end if

    if (present(startx)) then
       start2 = startx
       if (present(length)) then
          end2 = start2 + length - 1
       else
          end2 = totlen2
       end if
    else
       start2 = 1
       end2 = totlen2
    end if


    if (yfun%iscmplx) then
       yfun%z(start1:end1) = alpha * yfun%z(start1:end1) + &
                             (1.0_DP - alpha) * xfun%z(start2:end2)
    else
       yfun%d(start1:end1) = alpha * yfun%d(start1:end1) + &
                             (1.0_DP - alpha) * xfun%d(start2:end2)
    end if

  end subroutine data_functions_mix_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function data_coef_abs(mycoef)
    implicit none
    type(COEF), intent(in) :: mycoef
    if (mycoef%iscmplx) then
       data_coef_abs = abs(mycoef%z)
    else
       data_coef_abs = abs(mycoef%d)
    end if
  end function data_coef_abs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(COEF) pure elemental function data_coef_sum(coef1, coef2) result(csum)
    implicit none

    type(COEF), intent(in) :: coef1, coef2

    csum%iscmplx = coef1%iscmplx .or. coef2%iscmplx

    if (csum%iscmplx) then
       if (coef1%iscmplx) then
          csum%z = coef1%z
       else
          csum%z = cmplx(coef1%d, 0.0_DP, DP)
       end if
       if (coef2%iscmplx) then
          csum%z = csum%z + coef2%z
       else
          csum%z = csum%z + coef2%d
       end if
    else
       csum%d = coef1%d + coef2%d
    end if

  end function data_coef_sum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(COEF) pure elemental function data_coef_diff(coef1, coef2) result(cdiff)
    implicit none

    type(COEF), intent(in) :: coef1, coef2

    cdiff%iscmplx = coef1%iscmplx .or. coef2%iscmplx

    if (cdiff%iscmplx) then
       if (coef1%iscmplx) then
          cdiff%z = coef1%z
       else
          cdiff%z = cmplx(coef1%d, 0.0_DP, DP)
       end if
       if (coef2%iscmplx) then
          cdiff%z = cdiff%z - coef2%z
       else
          cdiff%z = cdiff%z - coef2%d
       end if
    else
       cdiff%d = coef1%d - coef2%d
    end if

  end function data_coef_diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(COEF) pure elemental function data_product_real_coef(x, mycoef) result(cprod)
    implicit none

    real(kind=DP), intent(in)   :: x
    type(COEF), intent(in)      :: mycoef

    cprod%iscmplx = mycoef%iscmplx

    if (mycoef%iscmplx) then
       cprod%z = x * mycoef%z
    else
       cprod%d = x * mycoef%d
    end if

  end function data_product_real_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(COEF) pure elemental function data_product_complex_coef(x, mycoef) result(cprod)
    implicit none

    complex(kind=DP), intent(in) :: x
    type(COEF), intent(in)       :: mycoef

    cprod%iscmplx = .True.

    if (mycoef%iscmplx) then
       cprod%z = x * mycoef%z
    else
       cprod%z = x * mycoef%d
    end if

  end function data_product_complex_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(COEF) pure elemental function data_product_coef(coef1, coef2) result(cprod)
    implicit none

    type(COEF), intent(in) :: coef1, coef2

    cprod%iscmplx = coef1%iscmplx .or. coef2%iscmplx

    if (cprod%iscmplx) then
       if (coef1%iscmplx) then
          cprod%z = coef1%z
       else
          cprod%z = cmplx(coef1%d, 0.0_DP, DP)
       end if
       if (coef2%iscmplx) then
          cprod%z = cprod%z * coef2%z
       else
          cprod%z = cprod%z * coef2%d
       end if
    else
       cprod%d = coef1%d * coef2%d
    end if

  end function data_product_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(COEF) pure elemental function data_division_coef(coef1, coef2) result(cdiv)
    implicit none

    type(COEF), intent(in) :: coef1, coef2

    cdiv%iscmplx = coef1%iscmplx .or. coef2%iscmplx

    if (cdiv%iscmplx) then
       if (coef1%iscmplx) then
          cdiv%z = coef1%z
       else
          cdiv%z = cmplx(coef1%d, 0.0_DP, DP)
       end if
       if (coef2%iscmplx) then
          cdiv%z = cdiv%z / coef2%z
       else
          cdiv%z = cdiv%z / coef2%d
       end if
    else
       cdiv%d = coef1%d / coef2%d
    end if

  end function data_division_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  type(COEF) function data_coef_sum_wopt(coef1, coef2, iscmplx) result(csum)
!    use utils, only: utils_assert
!    implicit none
!
!    type(COEF), intent(in) :: coef1, coef2
!    logical, optional, intent(in) :: iscmplx
!
!    csum%iscmplx = coef1%iscmplx .or. coef2%iscmplx
!    if (present(iscmplx)) then
!       call utils_assert(iscmplx .or. (.not. csum%iscmplx), 'ERROR in&
!            & data_coef_sum: you are casting a complex COEF into a real one.)')
!       csum%iscmplx = csum%iscmplx .or. iscmplx
!    end if
!
!    if (csum%iscmplx) then
!       if (coef1%iscmplx) then
!          csum%z = coef1%z
!       else
!          csum%z = cmplx(coef1%d, 0.0_DP, DP)
!       end if
!       if (coef2%iscmplx) then
!          csum%z = csum%z + coef2%z
!       else
!          csum%z = csum%z + coef2%d
!       end if
!    else
!       csum%d = coef1%d + coef2%d
!    end if
!
!  end function data_coef_sum_wopt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  type(COEF) function data_coef_product_wopt(coef1, coef2, iscmplx) result(cprod)
!    use utils, only: utils_assert
!    implicit none
!
!    type(COEF), intent(in) :: coef1, coef2
!    logical, optional, intent(in) :: iscmplx
!
!    cprod%iscmplx = coef1%iscmplx .or. coef2%iscmplx
!    if (present(iscmplx)) then
!       call utils_assert(iscmplx .or. (.not. cprod%iscmplx), 'ERROR in&
!            & data_coef_sum: you are casting a complex COEF into a real one.)')
!       cprod%iscmplx = cprod%iscmplx .or. iscmplx
!    end if
!
!    if (cprod%iscmplx) then
!       if (coef1%iscmplx) then
!          cprod%z = coef1%z
!       else
!          cprod%z = cmplx(coef1%d, 0.0_DP, DP)
!       end if
!       if (coef2%iscmplx) then
!          cprod%z = cprod%z * coef2%z
!       else
!          cprod%z = cprod%z * coef2%d
!       end if
!    else
!       cprod%d = coef1%d * coef2%d
!    end if
!
!  end function data_coef_product_wopt
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(FUNCTIONS) function data_functions_sum(func1, func2) result(fsum)
    use utils, only: utils_assert
    implicit none

    type(FUNCTIONS), intent(in) :: func1, func2

    call utils_assert(data_size(func1) == data_size(func2), 'Error in&
         & data_functions_sum: sizes of functions do not match.')

    call data_functions_alloc(fsum, data_size(func1), &
         func1%iscmplx .or. func2%iscmplx)

    if (fsum%iscmplx) then
       if (func1%iscmplx) then
          fsum%z(:) = func1%z
       else
          fsum%z(:) = cmplx(func1%d, 0.0_DP, DP)
       end if
       if (func2%iscmplx) then
          fsum%z(:) = fsum%z + func2%z
       else
          fsum%z(:) = fsum%z + func2%d
       end if
    else
       fsum%d(:) = func1%d + func2%d
    end if

  end function data_functions_sum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(FUNCTIONS) function data_functions_diff(func1, func2) result(fdiff)
    use utils, only: utils_assert
    implicit none

    type(FUNCTIONS), intent(in) :: func1, func2

    call utils_assert(data_size(func1) == data_size(func2), 'Error in&
         & data_functions_diff: sizes of functions do not match.')

    call data_functions_alloc(fdiff, data_size(func1), &
         func1%iscmplx .or. func2%iscmplx)

    if (fdiff%iscmplx) then
       if (func1%iscmplx) then
          fdiff%z(:) = func1%z
       else
          fdiff%z(:) = cmplx(func1%d, 0.0_DP, DP)
       end if
       if (func2%iscmplx) then
          fdiff%z(:) = fdiff%z - func2%z
       else
          fdiff%z(:) = fdiff%z - func2%d
       end if
    else
       fdiff%d(:) = func1%d - func2%d
    end if

  end function data_functions_diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(FUNCTIONS) function data_product_real_functions(x, func) result(fprod)
    implicit none

    real(kind=DP), intent(in)   :: x
    type(FUNCTIONS), intent(in) :: func

    call data_functions_alloc(fprod, data_size(func), func%iscmplx)

    if (func%iscmplx) then
       fprod%z(:) = x * func%z
    else
       fprod%d(:) = x * func%d
    end if

  end function data_product_real_functions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(FUNCTIONS) function data_product_complex_functions(x, func) result(fprod)
    implicit none

    complex(kind=DP), intent(in) :: x
    type(FUNCTIONS), intent(in)  :: func

    call data_functions_alloc(fprod, data_size(func), .True.)

    if (func%iscmplx) then
       fprod%z(:) = x * func%z
    else
       fprod%z(:) = x * func%d
    end if

  end function data_product_complex_functions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(FUNCTIONS) function data_product_coef_functions(mycoef, func) result(fprod)
    implicit none

    type(COEF), intent(in)      :: mycoef
    type(FUNCTIONS), intent(in) :: func

    call data_functions_alloc(fprod, data_size(func), &
         mycoef%iscmplx .or. func%iscmplx)

    if (fprod%iscmplx) then
       if (func%iscmplx) then
          fprod%z(:) = func%z
       else
          fprod%z(:) = cmplx(func%d, 0.0_DP, DP)
       end if
       if (mycoef%iscmplx) then
          fprod%z(:) = mycoef%z * fprod%z
       else
          fprod%z(:) = mycoef%d * fprod%z
       end if
    else
       fprod%d(:) = mycoef%d * func%d
    end if

  end function data_product_coef_functions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(FUNCTIONS) function data_product_functions(func1, func2) result(fprod)
    use utils, only: utils_assert
    implicit none

    type(FUNCTIONS), intent(in) :: func1, func2

    call utils_assert(data_size(func1) == data_size(func2), 'Error in&
         & data_product_functions: sizes of functions do not match.')

    call data_functions_alloc(fprod, data_size_functions(func1), &
         func1%iscmplx .or. func2%iscmplx)

    if (fprod%iscmplx) then
       if (func1%iscmplx) then
          fprod%z(:) = func1%z
       else
          fprod%z(:) = cmplx(func1%d, 0.0_DP, DP)
       end if
       if (func2%iscmplx) then
          fprod%z(:) = fprod%z * func2%z
       else
          fprod%z(:) = fprod%z * func2%d
       end if
    else
       fprod%d(:) = func1%d * func2%d
    end if

  end function data_product_functions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure subroutine data_coef_assign_real(c, x)
    implicit none
    type(COEF), intent(inout) :: c
    real(kind=DP), intent(in) :: x
    if (c%iscmplx) then
       c%z = cmplx(x, 0.0_DP, DP)
    else
       c%d = x
    end if
  end subroutine data_coef_assign_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine data_coef_assign_complex(c, x)
    use utils, only: utils_abort
    implicit none
    type(COEF), intent(inout) :: c
    complex(kind=DP), intent(in) :: x
    if (c%iscmplx) then
       c%z = x
    else
       call utils_abort('Error in data_coef_assign_complex: you cannot assign&
            & a complex value to a real coefficient.')
    end if
  end subroutine data_coef_assign_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function data_functions_abs(myfunc)
    use utils, only: utils_alloc_check
    implicit none
    ! output
    real(kind=DP), allocatable :: data_functions_abs(:)
    ! input
    type(FUNCTIONS), intent(in) :: myfunc
    ! local variable
    integer :: ierr

    allocate(data_functions_abs(data_size_functions(myfunc)), stat=ierr)
    call utils_alloc_check('data_functions_abs', 'data_functions_abs', ierr)

    if (myfunc%iscmplx) then
       data_functions_abs(:) = abs(myfunc%z)
    else
       data_functions_abs(:) = abs(myfunc%d)
    end if
  end function data_functions_abs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  subroutine data_coef_reduce(op, mycoef)
!    use comms, only: comms_reduce
!    implicit none
!
!    character(len=*), intent(in) :: op
!    type(COEF), intent(inout) :: mycoef
!    if (mycoef%iscmplx) then
!       call comms_reduce(op, mycoef%z)
!    else
!       call comms_reduce(op, mycoef%d)
!    end if
!  end subroutine data_coef_reduce

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure elemental subroutine data_coef_scale_real(mycoef, alpha)
    use constants, only: DP
    implicit none
    type(COEF), intent(inout) :: mycoef
    real(kind=DP), intent(in) :: alpha
    if (mycoef%iscmplx) then
       mycoef%z = mycoef%z * alpha
    else
       mycoef%d = mycoef%d * alpha
    end if
  end subroutine data_coef_scale_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  jme: impure only supported by the Intel Fortran Compiler from v16
!  impure elemental subroutine data_coef_scale_complex(mycoef, alpha)
  subroutine data_coef_scale_complex(mycoef, alpha)
    use constants, only: DP
    use utils, only: utils_abort
    implicit none
    type(COEF), intent(inout) :: mycoef
    complex(kind=DP), intent(in) :: alpha
    if (mycoef%iscmplx) then
       mycoef%z = mycoef%z * alpha
    else
       call utils_abort('Error in data_coef_scale_complex: product of &
            &real coefficient and complex factor not allowed.')
    end if
  end subroutine data_coef_scale_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure elemental subroutine data_functions_scale_real(myfunc, alpha)
    use constants, only: DP
    implicit none
    type(FUNCTIONS), intent(inout) :: myfunc
    real(kind=DP), intent(in) :: alpha
    if (myfunc%iscmplx) then
       myfunc%z(:) = myfunc%z * alpha
    else
       myfunc%d(:) = myfunc%d * alpha
    end if
  end subroutine data_functions_scale_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  jme: impure only supported by the Intel Fortran Compiler from v16
!  impure elemental subroutine data_functions_scale_complex(myfunc, alpha)
  subroutine data_functions_scale_complex(myfunc, alpha)
    use constants, only: DP
    use utils, only: utils_abort
    implicit none
    type(FUNCTIONS), intent(inout) :: myfunc
    complex(kind=DP), intent(in) :: alpha
    if (myfunc%iscmplx) then
       myfunc%z(:) = myfunc%z * alpha
    else
       call utils_abort('Error in data_coef_scale_complex: product of &
            &real function and complex factor not allowed.')
    end if
  end subroutine data_functions_scale_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  jme: impure only supported by the Intel Fortran Compiler from v16
!  impure elemental subroutine data_functions_scale_coef(myfunc, mycoef, same_type)
  subroutine data_functions_scale_coef(myfunc, mycoef, same_type)
    use utils, only: utils_assert
    implicit none
    type(FUNCTIONS), intent(inout) :: myfunc
    type(COEF), intent(in) :: mycoef
    logical, optional, intent(in) :: same_type

    ! jme: if required, check that either both are real or both are complex
    if (present(same_type)) then
       if (same_type) then
          call utils_assert(myfunc%iscmplx .eqv. mycoef%iscmplx, &
               'Error in data_functions_scale_coef: same type required, but &
               &complex character of function and coef arguments are:', &
               myfunc%iscmplx, mycoef%iscmplx)
       end if
    end if

    if (mycoef%iscmplx) then
       call data_functions_scale_complex(myfunc, mycoef%z)
    else
       call data_functions_scale_real(myfunc, mycoef%d)
    end if
  end subroutine data_functions_scale_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure elemental subroutine data_fftbox_scale_real(fftbd, alpha)
    use constants, only: DP
    implicit none
    type(FFTBOX_DATA), intent(inout) :: fftbd
    real(kind=DP), intent(in) :: alpha
    if (fftbd%iscmplx) then
       fftbd%z(:,:,:) = fftbd%z * alpha
    else
       fftbd%d(:,:,:) = fftbd%d * alpha
    end if
  end subroutine data_fftbox_scale_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  jme: impure only supported by the Intel Fortran Compiler from v16
!  impure elemental subroutine data_fftbox_scale_complex(fftbd, alpha)
  subroutine data_fftbox_scale_complex(fftbd, alpha)
    use constants, only: DP
    use utils, only: utils_abort
    implicit none
    type(FFTBOX_DATA), intent(inout) :: fftbd
    complex(kind=DP), intent(in) :: alpha
    if (fftbd%iscmplx) then
       fftbd%z(:,:,:) = fftbd%z * alpha
    else
       call utils_abort('Error in data_coef_scale_complex: product of &
            &real fftbox_data and complex factor not allowed.')
    end if
  end subroutine data_fftbox_scale_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  jme: impure only supported by the Intel Fortran Compiler from v16
!  impure elemental subroutine data_fftbox_scale_coef(fftbd, mycoef, same_type)
  subroutine data_fftbox_scale_coef(fftbd, mycoef, same_type)
    use utils, only: utils_assert
    implicit none
    type(FFTBOX_DATA), intent(inout) :: fftbd
    type(COEF), intent(in) :: mycoef
    logical, optional, intent(in) :: same_type

    ! jme: if required, check that either both are real or both are complex
    if (present(same_type)) then
       if (same_type) then
          call utils_assert(fftbd%iscmplx .eqv. mycoef%iscmplx, &
               'Error in data_fftbox_scale_coef: same type required, but &
               &complex character of fftbox_data and coef arguments are:', &
               fftbd%iscmplx, mycoef%iscmplx)
       end if
    end if

    if (mycoef%iscmplx) then
       call data_fftbox_scale_complex(fftbd, mycoef%z)
    else
       call data_fftbox_scale_real(fftbd, mycoef%d)
    end if
  end subroutine data_fftbox_scale_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure elemental subroutine data_set_coef_to_zero(mycoef)
    use constants, only: DP, cmplx_0
    implicit none
    type(COEF), intent(inout) :: mycoef
    if (mycoef%iscmplx) then
       mycoef%z = cmplx_0
    else
       mycoef%d = 0.0_DP
    end if
  end subroutine data_set_coef_to_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure elemental subroutine data_set_functions_to_zero(myfunc)
    use constants, only: DP, cmplx_0
    implicit none
    type(FUNCTIONS), intent(inout) :: myfunc
    if (myfunc%iscmplx) then
       myfunc%z(:) = cmplx_0
    else
       myfunc%d(:) = 0.0_DP
    end if
  end subroutine data_set_functions_to_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure elemental subroutine data_set_fftboxdata_to_zero(mydata)
    use constants, only: DP, cmplx_0
    implicit none
    type(FFTBOX_DATA), intent(inout) :: mydata
    if (mydata%iscmplx) then
       mydata%z(:,:,:) = cmplx_0
    else
       mydata%d(:,:,:) = 0.0_DP
    end if
  end subroutine data_set_fftboxdata_to_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer pure elemental function data_size_functions(myfunc) result(mysize)
    implicit none
    type(FUNCTIONS), intent(in) :: myfunc
    if (myfunc%iscmplx) then
       mysize = size(myfunc%z, 1)
    else
       mysize = size(myfunc%d, 1)
    end if
  end function data_size_functions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function data_size_fftboxdata(fftbd) result(mysize)
    implicit none
    integer :: mysize(3)
    type(FFTBOX_DATA), intent(in) :: fftbd
    if (fftbd%iscmplx) then
       mysize = size(fftbd%z)
    else
       mysize = size(fftbd%d)
    end if
  end function data_size_fftboxdata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine data_fftbox_copy(destfftbox, srcfftbox)
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(FFTBOX_DATA), intent(inout) :: destfftbox
    type(FFTBOX_DATA), intent(in) :: srcfftbox

    ! Trap real added to complex
    call utils_assert(srcfftbox%iscmplx.eqv.destfftbox%iscmplx, &
         "ERROR in data_fftbox_copy: FFTBox real/complex flags &
         &do not match.")

    ! Perform copy operation
    if (destfftbox%iscmplx) then
       destfftbox%z(:,:,:) = srcfftbox%z
    else
       destfftbox%d(:,:,:) = srcfftbox%d
    end if

  end subroutine data_fftbox_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine data_fftbox_axpy_real(yfftbox, xfftbox, alpha)
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(FFTBOX_DATA), intent(inout) :: yfftbox
    type(FFTBOX_DATA), intent(in) :: xfftbox
    real(kind=DP), optional :: alpha

    ! Local variables
    real(kind=DP) :: a

    a = 1.0_DP
    if (present(alpha)) a = alpha
    if (a==0.0_DP) return

    ! Trap real added to complex
    call utils_assert(xfftbox%iscmplx.eqv.yfftbox%iscmplx, &
         "ERROR in data_fftbox_axpy_real: FFTBox real/complex flags &
         &do not match.")

    ! Perform axpy operation
    if (yfftbox%iscmplx) then
       yfftbox%z(:,:,:) = a * xfftbox%z + yfftbox%z
    else
       yfftbox%d(:,:,:) = a * xfftbox%d + yfftbox%d
    end if

  end subroutine data_fftbox_axpy_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !============================================================================!
  ! This function calculates the dot product for an 2 equally sized arrays  of !
  ! FUNCTIONS.                                                                 !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, April 2018.                                    !
  ! Slightly modified by Joseph Prentice, April 2019
  !============================================================================!

  type(COEF) function data_functions_dot_array(bras,kets,startbra,startket,length)

    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNCTIONS), intent(in)   :: bras(:)
    type(FUNCTIONS), intent(in)   :: kets(:)
    integer, intent(in), optional :: startbra
    integer, intent(in), optional :: startket
    integer, intent(in), optional :: length

    ! Local variables
    integer :: isub, ierr
    type(COEF), allocatable :: funcs(:)

    ! rc2013: make sure the arrays are the same size
    call utils_assert(size(bras) == size(kets),'Error in &
         &data_functions_dot_array: bras and kets are different sizes.')
    allocate(funcs(size(bras)), stat=ierr)
    call utils_alloc_check('data_functions_dot_array', 'funcs', ierr)

    ! rc2013: I'm not sure if iscmplx has been correctly allocated yet...
    ! jcap: set iscmplx correctly
    if (any(bras(:)%iscmplx) .or. any(kets(:)%iscmplx)) then
       data_functions_dot_array%iscmplx=.true.
    else
       data_functions_dot_array%iscmplx=.false.
    end if
    call data_set_to_zero(data_functions_dot_array)
    do isub = 1,size(bras)
       funcs(isub) = data_functions_dot(bras(isub), kets(isub))
       if (data_functions_dot_array%iscmplx) then
          if (funcs(isub)%iscmplx) then
             data_functions_dot_array%z = &
                  data_functions_dot_array%z + funcs(isub)%z
          else
             data_functions_dot_array%z = &
                  data_functions_dot_array%z + cmplx(funcs(isub)%d,0.d0)
          end if
       else
          ! jcap: there should be no way that funcs can be complex but
          ! data_functions_dot_array is real
          data_functions_dot_array%d = &
               data_functions_dot_array%d + funcs(isub)%d
       end if
    end do
    deallocate(funcs, stat=ierr)
    call utils_dealloc_check('data_functions_dot_array', 'funcs', ierr)

  end function data_functions_dot_array


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module datatypes
