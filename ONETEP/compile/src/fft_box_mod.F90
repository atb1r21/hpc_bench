! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   This file was created in June 2019 by Nicholas Hine
!
!   The subroutines in this file were written by Chris-Kriton Skylaris,
!   Arash A. Mostofi, and Nicholas Hine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module fft_box

  use geometry, only: point
  use constants, only: DP

  implicit none

  private

  type FFTBOX_INFO
     ! lattice vectors
     type(POINT)      :: a1,a2,a3
     ! lattice vectors normalised to unit length
     type(POINT)      :: a1_unit,a2_unit,a3_unit
     ! reciprocal lattice vectors
     type(POINT)      :: b1,b2,b3
     ! regular grid spacing in atomic units
     real(kind=DP)    :: d1,d2,d3
     ! weight for integrals on the standard grid
     real(kind=DP)    :: weight
     ! total number of points in each lattice direction on the standard grid
     integer          :: total_pt1,total_pt2,total_pt3
     ! the total number of points in each lattice direction on the double grid
     integer          :: total_pt1_dbl,total_pt2_dbl,total_pt3_dbl
     ! pdh: additions to improve FFT efficiency
     integer          :: total_ld1,total_ld2
     integer          :: total_ld1_dbl,total_ld2_dbl
     ! reciprocal lattice vector grid
     real(kind=DP), allocatable, dimension(:,:,:,:) :: recip_grid
     ! cks: true if FFT box length coincides with sim cell along a1
     logical          :: coin1
     ! cks: true if FFT box length coincides with sim cell along a2
     logical          :: coin2
     ! cks: true if FFT box length coincides with sim cell along a3
     logical          :: coin3
     ! Index of the entry for this grid in the array of serial_fft3d_info
     ! structures in fourier_mod (needed to perform box FFTs)
     integer          :: fft_index
  end type FFTBOX_INFO

  public :: FFTBOX_INFO
  public :: fftbox_init
  public :: fftbox_exit
  public :: fftbox_fourier_size
  public :: fftbox_find_size
  public :: fftbox_find_size_basic

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fftbox_init(fftbox, n1, n2, n3, cell)

    !============================================================!
    ! This subroutine initialises the components of an           !
    ! FFTBOX_INFO  type variable which groups together           !
    ! information about the FFTbox.                              !
    !------------------------------------------------------------!
    ! WARNING: This subroutine should be called only AFTER       !
    !          the cell variable has been initialised!           !
    !------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 13/5/2003.             !
    ! Reciprocal grid array added by Peter Haynes on 15/7/2004   !
    ! Reciprocal grid fine array added by Quintin Hill on        !
    ! 24/04/2008.                                                !
    !============================================================!

    use geometry, only: operator(*)
    use rundat, only: pub_dbl_grid_scale
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: n1, n2, n3
    type(FFTBOX_INFO), intent(inout) :: fftbox
    type(CELL_INFO), intent(in) :: cell

    ! Local variables
    integer :: ierr
    integer :: i1,i2,i3
    integer :: k1,k2,k3
    integer :: n1half,n2half,n3half
    real(kind=DP) :: b1(3),b2(3),b3(3)
    real(kind=DP) :: g3(3),g23(3),g(3)
    real(kind=DP) :: gsq
    real(kind=DP) :: coulomb_cutoff

    if (mod(n1,2) == 0 .or. mod(n2,2) == 0 .or. mod(n3,2) == 0) then
       call utils_abort('Error in fftbox_init: at least one &
            &dimension of the FFT box is even. n1, n2, n3 follow',n1,n2,n3)
    end if

    fftbox%total_pt1 = n1
    fftbox%total_pt2 = n2
    fftbox%total_pt3 = n3

    ! cks: set flags if FFT box coincides with sim cell along any lattice vector
    fftbox%coin1 =.false.
    if (fftbox%total_pt1 == cell%total_pt1) fftbox%coin1 =.true.
    fftbox%coin2 =.false.
    if (fftbox%total_pt2 == cell%total_pt2) fftbox%coin2 =.true.
    fftbox%coin3 =.false.
    if (fftbox%total_pt3 == cell%total_pt3) fftbox%coin3 =.true.

    if (pub_dbl_grid_scale > 1.0_DP) then
       fftbox%total_pt1_dbl = 2*n1
       fftbox%total_pt2_dbl = 2*n2
       fftbox%total_pt3_dbl = 2*n3
    else
       fftbox%total_pt1_dbl = n1
       fftbox%total_pt2_dbl = n2
       fftbox%total_pt3_dbl = n3
    end if

    fftbox%a1 =( real(n1, kind=DP)/real(cell%total_pt1, kind=DP) ) &
         * cell%a1
    fftbox%a2 =( real(n2, kind=DP)/real(cell%total_pt2, kind=DP) ) &
         * cell%a2
    fftbox%a3 =( real(n3, kind=DP)/real(cell%total_pt3, kind=DP) ) &
         * cell%a3

    fftbox%a1_unit =cell%a1_unit
    fftbox%a2_unit =cell%a2_unit
    fftbox%a3_unit =cell%a3_unit

    fftbox%b1 =( real(cell%total_pt1, kind=DP)/real(n1, kind=DP) ) &
         * cell%b1
    fftbox%b2 =( real(cell%total_pt2, kind=DP)/real(n2, kind=DP) ) &
         * cell%b2
    fftbox%b3 =( real(cell%total_pt3, kind=DP)/real(n3, kind=DP) ) &
         * cell%b3

    fftbox%d1 =cell%d1
    fftbox%d2 =cell%d2
    fftbox%d3 =cell%d3

    fftbox%weight =cell%weight

    ! pdh: improve efficiency of FFTs - assumed complex-to-complex here
    ! aam: only ever do C2C FFTs on the FFT box, therefore don't require
    !      the extra two rows in the first dimension of either the coarse
    !      or fine FFT boxes.
    ! pdh: And FFTw will not work if ld1 /= n1 and ld2 /= n2!
    !    fftbox%total_ld1=fftbox%total_pt1+2
    fftbox%total_ld1=fftbox%total_pt1
    fftbox%total_ld2=fftbox%total_pt2
    !    fftbox%total_ld1_dbl=fftbox%total_pt1_dbl+2
    fftbox%total_ld1_dbl=fftbox%total_pt1_dbl
    fftbox%total_ld2_dbl=fftbox%total_pt2_dbl

    ! ndmh: initialise index to internal storage in fourier_mod
    fftbox%fft_index = 0

    ! pdh: add reciprocal grid
    allocate(fftbox%recip_grid(6,fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('fftbox_init','fftbox%recip_grid', ierr)

    ! local copies of FFT box reciprocal lattice vectors
    b1(1) = fftbox%b1%x ; b1(2) = fftbox%b1%y ; b1(3) = fftbox%b1%z
    b2(1) = fftbox%b2%x ; b2(2) = fftbox%b2%y ; b2(3) = fftbox%b2%z
    b3(1) = fftbox%b3%x ; b3(2) = fftbox%b3%y ; b3(3) = fftbox%b3%z

    coulomb_cutoff = 0.48_DP*min(&
         (real(fftbox%total_pt1,kind=DP)*fftbox%d1),&
         (real(fftbox%total_pt2,kind=DP)*fftbox%d2),&
         (real(fftbox%total_pt3,kind=DP)*fftbox%d3))

    ! loop over FFT box reciprocal grid
    fftbox%recip_grid = 0.0_DP
    n1half = fftbox%total_pt1/2+1
    n2half = fftbox%total_pt2/2+1
    n3half = fftbox%total_pt3/2+1
    do i3=1,fftbox%total_pt3
       if (i3 > n3half) then
          k3 = i3 - fftbox%total_pt3 - 1
       else
          k3 = i3 - 1
       end if
       g3 = k3 * b3
       do i2=1,fftbox%total_pt2
          if (i2 > n2half) then
             k2 = i2 - fftbox%total_pt2 - 1
          else
             k2 = i2 - 1
          end if
          g23 = g3 + k2 * b2
          do i1=1,fftbox%total_pt1
             if (i1 > n1half) then
                k1 = i1 - fftbox%total_pt1 - 1
             else
                k1 = i1 - 1
             end if
             g = g23 + k1 * b1
             fftbox%recip_grid(1:3,i1,i2,i3) = g
             gsq = g(1)*g(1) + g(2)*g(2) + g(3)*g(3)
             fftbox%recip_grid(4,i1,i2,i3) = sqrt(gsq)
             fftbox%recip_grid(5,i1,i2,i3) = 0.5_DP * gsq
             ! qoh: Treat case of G=0 specially
             if (i1 == 1 .and. i2 == 1 .and. i3 == 1) then
                fftbox%recip_grid(6,1,1,1)=0.5_DP*coulomb_cutoff**2
             else
                fftbox%recip_grid(6,i1,i2,i3) &
                     = (1.0_DP - cos(sqrt(gsq)*coulomb_cutoff)) / gsq
             end if
          end do
       end do
    end do

  end subroutine fftbox_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fftbox_exit(fftbox)

    !============================================================!
    ! This subroutine deallocates pointers in fftbox arrays.     !
    !------------------------------------------------------------!
    ! Written by Victor Milman on 03/11/2006.                    !
    ! Modified by Quintin Hill on 24/04/2008.                    !
    !============================================================!

    use utils, only : utils_dealloc_check

    implicit none

    type(FFTBOX_INFO), intent(inout) :: fftbox

    integer :: ierr

    if (allocated(fftbox%recip_grid)) then
       deallocate(fftbox%recip_grid,stat=ierr)
       call utils_dealloc_check('fftbox_exit','fftbox%recip_grid',ierr)
    end if

  end subroutine fftbox_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fftbox_fourier_size(n,flag)

    implicit none

    ! Arguments
    integer, intent(inout) :: n           ! The FFT grid size
    integer, intent(in), optional :: flag ! Flag specifying whether result
    !                                       must be even or odd

    ! The following (FFT library dependent) parameters determine which
    ! factors are allowed:
    !   num_factors : the number of factors available
    !   factor(:)   : the factors themselves (in ascending order)
    !   limit(:)    : the limits on the power of the factors (-ve -> unlimited)

#ifdef MKL_FFTW3
#define FFTW3
#endif
#ifdef FFTW3_NO_OMP
#ifndef FFTW3
#define FFTW3
#endif
#endif
#ifdef FFTW
    integer, parameter :: num_factors = 6
    integer, parameter :: factor(num_factors) = (/ 2,3,5,7,11,13 /)
    integer, parameter :: limit(num_factors) = (/ -1,-1,-1,-1,1,1 /)
#endif
#ifdef FFTW3
    integer, parameter :: num_factors = 6
    integer, parameter :: factor(num_factors) = (/ 2,3,5,7,11,13 /)
    integer, parameter :: limit(num_factors) = (/ -1,-1,-1,-1,1,1 /)
#endif
#ifdef ACML
    integer, parameter :: num_factors = 6
    integer, parameter :: factor(num_factors) = (/ 2,3,5,7,11,13 /)
    integer, parameter :: limit(num_factors) = (/ -1,-1,-1,-1,1,1 /)
#endif
#ifdef VENDOR
#ifdef ALPHA
    integer, parameter :: num_factors = 5
    integer, parameter :: factor(num_factors) = (/ 2,3,5,7,11 /)
    integer, parameter :: limit(num_factors) = (/ -1,-1,-1,1,1 /)
#endif
#ifdef SUN
    integer, parameter :: num_factors = 7
    integer, parameter :: factor(num_factors) = (/ 2,3,4,5,7,11,13 /)
    integer, parameter :: limit(num_factors) = (/ -1,-1,-1,-1,-1,-1,-1 /)
#endif
#endif

#if !defined(FFTW) && !defined(FFTW3) && !defined(ACML) && !defined(SUN)
#error "ONETEP compilation ERROR: One of FFTW, FFTW3, MKL_FFTW3, ACML or SUN must be defined."
#endif

    ! Local variables
    integer :: ifac
    integer :: fac
    integer :: power
    integer :: ntest
    integer :: stride

    ! Set starting point
    if (present(flag)) then
       if (mod(flag,2) == 0) then ! force n to be even
          n = n+mod(n,2)
       else                       ! force n to be odd
          n = n-mod(n,2)+1
       end if
    end if

    ! Set stride
    if (present(flag)) then
       stride = 2
    else
       stride = 1
    end if

    ! Check whether n is an allowed product of factors
    do
       ntest = n
       do ifac=num_factors,1,-1 ! take out largest factors first
          fac = factor(ifac)
          power = limit(ifac)
          do
             if (mod(ntest,fac) /= 0 .or. power == 0) exit
             ntest = ntest / fac
             power = power - 1
          end do
       end do
       if (ntest == 1) exit
       n = n + stride
    end do

  end subroutine fftbox_fourier_size


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine fftbox_find_size_basic(n1, n2, n3,  &              ! output
       max_rad, multiplier, extend, halo, a1, a2, a3, &         ! input
       b1, b2, b3, d1, d2, d3, total_pt1, total_pt2, total_pt3) ! input

    !==================================================================!
    ! This subroutine returns the number of grid points in each        !
    ! lattice vector direction for the FFTbox.                         !
    !------------------------------------------------------------------!
    ! WARNING: At present, the number of points returned by this       !
    !          subroutine is correct only when orthorhombic simulation !
    !          cells are used. For other types of simulation cells,    !
    !          where either of the alpha or beta or gamma angles are   !
    !          not 90 degrees, this subroutine will underestimate      !
    !          the size of the FFTbox, and the only way to set its     !
    !          size right is by the user specifying it in the          !
    !          input file.                                             !
    !------------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris in 2001 as           !
    ! basis_pair_box_specification to work with a "pair-box".          !
    ! Modified by Arash Mostofi in April 2003 to work with             !
    ! a "triple-box".                                                  !
    ! Modified by Chris-Kriton Skylaris on 18/11/2003 so that it       !
    ! runs in parallel.                                                !
    ! Modified by Peter D. Haynes in June 2004 so that it accepts      !
    ! larger, user-defined FFTbox dimensions.                          !
    ! Modified by Chris-Kriton Skylaris on 06/11/2004 in order to      !
    ! work as part of ppd_strategy.                                    !
    ! Modified by Chris-Kriton Skylaris on 11/02/2008 in order to      !
    ! work with non-orthogonal lattice vectors.                        !
    !==================================================================!

    use constants, only: DP
    use geometry, only: POINT, geometry_magnitude, operator(.DOT.)

    implicit none

    ! Arguments
    integer, intent(out)       :: n1, n2, n3
    type(POINT), intent(in)    :: a1, a2, a3  ! lattice vectors
    type(POINT), intent(in)    :: b1, b2, b3  ! reciprocal vectors
    real (kind=dp), intent(in) :: d1, d2, d3  ! grid spacing
    integer, intent(in)        :: total_pt1, total_pt2, & ! total grid points
                                  total_pt3
    real(kind=DP), intent(in)  :: max_rad     ! maximum sphere radius
    integer, intent(in)        :: multiplier
    real(kind=DP), intent(in)  :: halo
    logical, intent(in)        :: extend(3)

    ! Local variables
    integer :: same_n ! common value for pairs of n1, n2, n3
    integer :: odd_grid_pt1 ! highest odd number <= grid%total_pt1
    integer :: odd_grid_pt2 ! highest odd number <= grid%total_pt2
    integer :: odd_grid_pt3 ! highest odd number <= grid%total_pt3
    real(kind=DP)  :: proj1 ! cosine of angle between a1 and a2 x a3
    real(kind=DP)  :: proj2 ! cosine of angle between a2 and a3 x a1
    real(kind=DP)  :: proj3 ! cosine of angle between a3 and a2 x a1

    call internal_check_rad_vs_grid(max_rad, d1)
    call internal_check_rad_vs_grid(max_rad, d2)
    call internal_check_rad_vs_grid(max_rad, d3)

    proj1 = (a1 .DOT. b1) / ( geometry_magnitude(a1) * geometry_magnitude(b1) )
    proj2 = (a2 .DOT. b2) / ( geometry_magnitude(a2) * geometry_magnitude(b2) )
    proj3 = (a3 .DOT. b3) / ( geometry_magnitude(a3) * geometry_magnitude(b3) )

    call internal_np_init(n1, odd_grid_pt1, proj1, d1, total_pt1)
    call internal_np_init(n2, odd_grid_pt2, proj2, d2, total_pt2)
    call internal_np_init(n3, odd_grid_pt3, proj3, d3, total_pt3)

    ! cks: if grid spacing is the same along two directions make sure FFT box
    ! cks: points are the same
    ! ndmh: added protection against dimensions becoming even
    ! jd: added missing abs that ndmh had spotted
    if ( abs(d1 - d2) < epsilon(1.0_DP) ) then
       same_n = max(n1, n2)
       if (same_n <= odd_grid_pt1) n1 = same_n
       if (same_n <= odd_grid_pt2) n2 = same_n
    endif
    if ( abs(d1 - d3) < epsilon(1.0_DP) ) then
       same_n = max(n1, n3)
       if (same_n <= odd_grid_pt1) n1 = same_n
       if (same_n <= odd_grid_pt3) n3 = same_n
    endif
    if ( abs(d3 - d2) < epsilon(1.0_DP) ) then
       same_n = max(n3, n2)
       if (same_n <= odd_grid_pt3) n3 = same_n
       if (same_n <= odd_grid_pt2) n2 = same_n
    endif

    ! agreco: if extended NGWFs along a lattice vector,
    ! set FFT box equal to the grid along that vector
    ! jme: previously forced to be odd
    if (extend(1)) n1 = total_pt1
    if (extend(2)) n2 = total_pt2
    if (extend(3)) n3 = total_pt3

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_check_rad_vs_grid(max_rad_xx, d_xx)

      use utils, only: utils_assert

      implicit none
      real(kind=DP), intent(in) :: max_rad_xx  ! maximum ngwf radius
      real(kind=DP), intent(in) :: d_xx        ! psinc grid spacing


      call utils_assert(max_rad_xx >= d_xx, 'Maximum NGWF radius is&
           & smaller than psinc grid spacing.')

    end subroutine internal_check_rad_vs_grid

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_np_init(np,odd_grid_pt,proj,d_xx,total_pt)

       implicit none

       ! Arguments
       integer, intent(out) :: np          ! number of points along direction
       integer, intent(out) :: odd_grid_pt ! max. odd num. of points
       real (kind=dp), intent(in) :: proj  ! cosine of angle (a_i, b_i)
       real (kind=dp), intent(in) :: d_xx  ! psinc grid spacing
       integer, intent(in) :: total_pt ! total number of points

       ! Auxiliary variables
       integer ::  np_temp

#ifdef OLD_FFTBOX_SIZE
       np = multiplier*( int(2.0_DP*real(max_rad, kind=DP)/(abs(proj)*d_xx)) +1)
#else
       np = multiplier*( int(2.0_DP*real(max_rad, kind=DP)/(abs(proj)*d_xx)) +2)
#endif
       ! jd: Rationale for '2': email exchange between JD, NDMH, CKS, JMA, JCW
       !     on 2018.09.28. Fixes bug reported by Alexander Perlov on 2018.09.18.
       !     This changes the default FFT box size, so caveat emptor.

       ! cks: increase points to account for halo region if there is one
       if (halo > 0.0_DP) then
          np = np + int(multiplier*2.0_DP*halo/d_xx)
       end if

       ! ensures that all dimensions are odd
       if ( mod(np,2) == 0 ) np = np + 1

       ! ndmh: find odd number equal to or below grid dimensions
       odd_grid_pt = total_pt - 1 + mod(total_pt,2)

       ! cks: Increase dimensions until Fourier transform efficiency
       ! cks: provided it does not exceed the grid
       ! ndmh: added protection against grid size becoming even
       ! jme: ensure here that the FFT box is smaller than the grid
       np_temp = np
       call fftbox_fourier_size(np, 1)
       if (np > odd_grid_pt) np = min(np_temp, odd_grid_pt)

    end subroutine internal_np_init

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  end subroutine fftbox_find_size_basic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine fftbox_find_size(n1, n2, n3, &          ! output
       max_rad, multiplier, extend, halo, pref_size, &
       a1,a2,a3,b1,b2,b3,d1,d2,d3,total_pt1,total_pt2,total_pt3) ! input

    use constants, only: stdout, VERBOSE
    use comms, only: pub_on_root
    use rundat, only: pub_output_detail

    implicit none

    ! Arguments
    integer, intent(out)        :: n1, n2, n3
    real(kind=DP), intent(in)   :: max_rad
    real(kind=DP), intent(in)   :: halo
    integer, intent(in)         :: multiplier
    integer, intent(in)         :: pref_size(3)
    logical, intent(in)         :: extend(3)
    type(POINT), intent(in)    :: a1, a2, a3  ! lattice vectors
    type(POINT), intent(in)    :: b1, b2, b3  ! reciprocal vectors
    real (kind=dp), intent(in) :: d1, d2, d3  ! grid spacing
    integer, intent(in)        :: total_pt1, total_pt2, & ! total grid points
                                  total_pt3

    call fftbox_find_size_basic(n1, n2, n3, &
         max_rad, multiplier, extend, halo, &
         a1, a2, a3, b1, b2, b3, d1, d2, d3, total_pt1, total_pt2, total_pt3)

    ! pdh: increase FFT box dimensions to specified preference
    if (pub_on_root .and. pub_output_detail >= VERBOSE .and. &
         (pref_size(1) > n1 .or. pref_size(2) > n2 .or. &
         pref_size(3) > n3)) &
         write(stdout,'(/a)') 'WARNING in fftbox_find_size:'
    if (pref_size(1) > n1) then
       if (pub_on_root .and. pub_output_detail >= VERBOSE) &
            write(stdout,'(2(a,i3))') '  FFT-box dimension 1 increased from ',n1, &
                 ' to preferred size ',pref_size(1)
       n1 = pref_size(1)
    end if
    if (pref_size(2) > n2) then
       if (pub_on_root.and. pub_output_detail >= VERBOSE) &
            write(stdout,'(2(a,i3))') '  FFT-box dimension 2 increased from ',n2, &
                 ' to preferred size ',pref_size(2)
       n2 = pref_size(2)
    end if
    if (pref_size(3) > n3) then
       if (pub_on_root.and. pub_output_detail >= VERBOSE) &
            write(stdout,'(2(a,i3))') '  FFT-box dimension 3 increased from ',n3, &
                 ' to preferred size ',pref_size(3)
       n3 = pref_size(3)
    end if

  end subroutine fftbox_find_size

end module fft_box
