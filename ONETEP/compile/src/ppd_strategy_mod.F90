! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                        PPD strategy module                     !
!----------------------------------------------------------------!
! Written by Chris-Kriton Skylaris, 6/11/2004                    !
!================================================================!

module ppd_strategy

  implicit none

  private

  public :: ppd_strategy_determine
  public :: ppd_strategy_check_and_print

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ppd_strategy_determine(d1, d2, d3, n_ppd_pt1, n_ppd_pt2, n_ppd_pt3, & !output
       a1, a2, a3, elements, grid_scale)                                           !input
    !====================================================================!
    ! This subroutine determines the PSINC grid spacing and the number   !
    ! of points per PPD.                                                 !
    !--------------------------------------------------------------------!
    ! NOTE. The grid spacing is determined according to the following    !
    ! algorithm: First it is calculated according to the given kinetic   !
    ! energy cutoff based on the condition of equal volumes between the  !
    ! ONETEP cube of plane wave vectors and the sphere of wave vectors   !
    ! of conventional plane wave codes. Then the grid spacing in each    !
    ! lattice direction is made finer until it corresponds to an odd     !
    ! number of points with as many prime factors (for efficient FFTs)   !
    ! as possible. If a nonzero value is given in the input for any      !
    ! lattice vector direction, this value is used instead and overrides !
    ! what is selected by the previous procedure. The number of points   !
    ! per ppd is determined so that an integer number of ppds fits in    !
    ! each lattice vector direction and again there, the value obtained  !
    ! can be overriden by a non-zero value specified in the input.       !
    !--------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 06/11/2004.                    !
    ! Modified by Chris-Kriton Skylaris on 07/08/2005 to add even psinc  !
    ! grid capability.                                                   !
    !====================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, VERBOSE, stdout, PI
    use fft_box, only: fftbox_find_size_basic, fftbox_fourier_size
    use geometry, only: POINT, operator(.cross.), operator(.dot.),operator(*)
    use ion, only : ELEMENT
    use rundat, only: pub_cutoff_energy, pub_output_detail, pub_ppd_npoints, &
        pub_psinc_spacing, pub_cond_calculate_any_task, pub_extend_ngwf, &
        pub_ngwf_halo
    use utils, only: utils_abort, utils_assert

    implicit none


    real(kind=DP), intent(out) :: d1 ! Standard grid spacing in direction 1
    real(kind=DP), intent(out) :: d2 ! Standard grid spacing in direction 2
    real(kind=DP), intent(out) :: d3 ! Standard grid spacing in direction 3
    integer, intent(out) :: n_ppd_pt1 ! Number of points per ppd in direction 1
    integer, intent(out) :: n_ppd_pt2 ! Number of points per ppd in direction 2
    integer, intent(out) :: n_ppd_pt3 ! Number of points per ppd in direction 3
    type(POINT), intent(in) :: a1  ! Lattice vector 1
    type(POINT), intent(in) :: a2  ! Lattice vector 2
    type(POINT), intent(in) :: a3  ! Lattice vector 3
    type(ELEMENT), intent(in) :: elements(:)
    real(kind=dp), optional, intent(in) :: grid_scale

    ! Local variables
    real(kind=DP) :: ddd    ! Grid-spacing corresponding to KE cutoff
    real(kind=DP) :: max_rad, min_fftbox_width
    integer :: n_lat_pt1     ! Number of grid points along a1
    integer :: n_lat_pt2     ! Number of grid points along a2
    integer :: n_lat_pt3     ! Number of grid points along a3
    logical :: ext1          ! true if extended NGWFs along a1
    logical :: ext2          ! true if extended NGWFs along a2
    logical :: ext3          ! true if extended NGWFs along a3
    type(POINT) :: b1, b2, b3 ! Reciprocal vectors
    logical :: force_odd_size(3)   ! Force odd sizes in cell.
    integer :: attempt
    integer :: i, nfft(3), nfft_bak(3)

    ! cks: make output info more clear
    if ( pub_on_root .and. pub_output_detail >= VERBOSE) write(stdout,*)'   '

    ext1 = pub_extend_ngwf(1)
    ext2 = pub_extend_ngwf(2)
    ext3 = pub_extend_ngwf(3)

    ! Reciprocal vectors
    b1 = (2.0_DP*PI/( a1.DOT.(a2.CROSS.a3) ) ) * (a2.CROSS.a3)
    b2 = (2.0_DP*PI/( a1.DOT.(a2.CROSS.a3) ) ) * (a3.CROSS.a1)
    b3 = (2.0_DP*PI/( a1.DOT.(a2.CROSS.a3) ) ) * (a1.CROSS.a2)

    ! cks: grid spacing in 1D corresponding to KE energy cutoff
    if (pub_cutoff_energy > 0.0_DP) then
       ddd = ( (6.0 *PI*PI)**(1.0_DP/3.0_DP) ) / sqrt(2.0_DP*pub_cutoff_energy)
    else
       call utils_abort(&
            'Non-positive kinetic energy cutoff in ppd_strategy_determine.')
    endif

    ! cks: minimum required width of FFT-box along any direction
    max_rad = maxval(elements(:)%radius)
    if (pub_cond_calculate_any_task) then
       max_rad = max(max_rad, maxval(elements(:)%radius_cond))
    end if
    min_fftbox_width = 6.0_DP * max_rad

    ! Try to get cells that will lead to acceptable FFT sizes without forcing
    ! odd numbers of points along directions of the cell where ngwfs are not
    ! extended.  Force odd numbers if unsuccessful.
    force_odd_size = pub_extend_ngwf
    do attempt = 1, 2

       ! Determine grid spacing and number of points in each direction of the cell.
       call internal_find_dpsinc_in_latvec(d1, n_lat_pt1, &
            a1, b1, pub_psinc_spacing(1), "a1", force_odd_size(1))

       call internal_find_dpsinc_in_latvec(d2, n_lat_pt2, &
            a2, b2, pub_psinc_spacing(2), "a2", force_odd_size(2))

       call internal_find_dpsinc_in_latvec(d3, n_lat_pt3, &
            a3, b3, pub_psinc_spacing(3), "a3", force_odd_size(3))

       ! Simulate computation of FFT box size.
       call fftbox_find_size_basic(nfft(1), nfft(2), nfft(3), &
            max_rad, 3, pub_extend_ngwf, pub_ngwf_halo, a1, a2, a3, b1, b2, b3, &
            d1, d2, d3, n_lat_pt1, n_lat_pt2, n_lat_pt3)

       ! At the end of the first attempt, check whether the FFTs will be
       ! computed efficiently and try a second attempt otherwise.
       if (attempt == 1) then
          nfft_bak = nfft
          do i = 1, 3
             call fftbox_fourier_size(nfft(i), 1)
          end do
          if ( ALL ( nfft == nfft_bak ) ) then
             exit
          else
             do i = 1, 3
                if (nfft(i) /= nfft_bak(i)) force_odd_size(i) = .true.
             end do
          end if
        end if

    end do

    ! Determine number of points per ppd.
    ! jd: 2017.01: Changed 'i2' to 'i0', as occasionally a single huge PPD
    !              commensurate with the cell is useful for debugging, and
    !              with i2 'output conversion error' would ensue.

    ! Along a1 (n_ppd_pt1).
    ! override number of points per ppd if extended NGWFs are used
    if ((pub_ppd_npoints(1) > 0) .and. (.not.ext1)) then
       n_ppd_pt1 =pub_ppd_npoints(1)
       if (pub_on_root .and. pub_output_detail >= VERBOSE) &
            write(stdout,'(a,i0,a)') '  Number of PPD points along a1 set to ',&
                 n_ppd_pt1,' from input file.'
    else
       call internal_find_npts_per_ppd_1d(n_ppd_pt1, &
           n_lat_pt1, ext1)
       if (pub_on_root.and.pub_output_detail >= VERBOSE.and.ext1) &
            write(stdout,'(a)') '  Using extended NGWFs along a1.'
    endif
    call utils_assert(n_ppd_pt1 /= 0, &
         'Failed to initialise n_ppd_pt1 in ppd_strategy_determine.')

    ! Along a2 (n_ppd_pt2).
    ! override number of points per ppd if extended NGWFs are used
    if ((pub_ppd_npoints(2) > 0) .and. (.not.ext2)) then
       n_ppd_pt2 =pub_ppd_npoints(2)
       if (pub_on_root .and. pub_output_detail >= VERBOSE) &
            write(stdout,'(a,i0,a)') '  Number of PPD points along a2 set to ',&
                 n_ppd_pt2,' from input file.'
    else
       call internal_find_npts_per_ppd_1d(n_ppd_pt2, &
            n_lat_pt2, ext2)
       if (pub_on_root.and.pub_output_detail >= VERBOSE.and.ext2) &
            write(stdout,'(a)') '  Using extended NGWFs along a2.'
    endif
    call utils_assert(n_ppd_pt2 /= 0, &
         'Failed to initialise n_ppd_pt2 in ppd_strategy_determine.')

    ! Along a3 (n_ppd_pt3).
    ! override number of points per ppd if extended NGWFs are used
    if ((pub_ppd_npoints(3) > 0) .and. (.not.ext3)) then
       n_ppd_pt3 =pub_ppd_npoints(3)
       if (pub_on_root .and. pub_output_detail >= VERBOSE) then
          if (n_ppd_pt3 == 1) then
             write(stdout,'(a,i0,a,/,a)') '  Number of PPD points along a3 set to ',&
                  n_ppd_pt3, '.  This is ', '  the default value, it may have been &
                  &provided also in the input file.'
          else
             write(stdout,'(a,i0,a)') '  Number of PPD points along a3 set to ',&
                  n_ppd_pt3,' from input file.'
          end if
       end if
    else
       call internal_find_npts_per_ppd_1d(n_ppd_pt3, &
            n_lat_pt3, ext3)
       if (pub_on_root.and.pub_output_detail >= VERBOSE.and.ext3) &
            write(stdout,'(a)') '  Using extended NGWFs along a3.'
    end if
    call utils_assert(n_ppd_pt3 /= 0, &
         'Failed to initialise n_ppd_pt3 in ppd_strategy_determine.')

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_find_dpsinc_in_latvec(dxx, n_lat_pt, &
         a_lat, b_recip, pub_psinc_spacing, lat_text, force_odd)
      ! cks: Fixed bugs in the manual specification of psinc-grid, 23/04/2010

      use constants, only: TWO_PI
      use fft_box, only: fftbox_fourier_size
      use geometry, only: geometry_magnitude
      use rundat, only: pub_even_psinc_grid, pub_odd_psinc_grid

      implicit none

      ! Arguments
      real(kind=DP),    intent(out)   :: dxx
      integer,          intent(out)   :: n_lat_pt
      real(kind=DP),    intent(in)    :: pub_psinc_spacing
      type(POINT),      intent(in)    :: a_lat
      type(POINT),      intent(in)    :: b_recip
      character(len=*), intent(in)    :: lat_text
      logical,          intent(inout) :: force_odd

      ! Local variables
      real(kind=DP) :: len_a
      real(kind=DP) :: plane_dist


      len_a = geometry_magnitude(a_lat)
      n_lat_pt = nint(len_a/ddd)


      ! cks: distance between consecutive planes defined by other two lattice vectors
      plane_dist = TWO_PI / geometry_magnitude(b_recip)
      ! cks: decide whether FFT box should coincide with simulation cell
      ! cks: and force odd number of points if true
      if ( plane_dist <= min_fftbox_width ) force_odd =.true.

      call utils_assert(.not. force_odd .or. .not. pub_even_psinc_grid, &
           'FFT box needs to coincide with simulation cell along '//&
           trim(lat_text)//' but even_psinc_grid is specified. &
           &These are incompatible.')

      if ( (pub_psinc_spacing > 0.0_DP) ) then
         ! cks: apply pub_psinc_spacing settings from input file

         dxx =pub_psinc_spacing
         if (pub_on_root .and. pub_output_detail >= VERBOSE) &
              write(stdout,'(3a,f16.12,a)') &
              '  PSINC grid spacing along ',lat_text,' set to ',dxx,'a0 from input file.'

         n_lat_pt =nint(len_a/dxx)

         if ( (mod(n_lat_pt, 2) == 0) .and. force_odd) then

            call utils_abort('FFT box needs to coincide with simulation cell &
                 &along '//trim(lat_text)//'. Adjust value of pub_psinc_spacing to &
                 &produce an odd number of psinc points.')
         endif

      else

         ! cks: determine optimum number of psinc-grid points
         if (pub_odd_psinc_grid .or. force_odd) then
            call fftbox_fourier_size(n_lat_pt, 1)
         else if (pub_even_psinc_grid) then
            call fftbox_fourier_size(n_lat_pt, 0)
         else
            call fftbox_fourier_size(n_lat_pt)
         endif

      endif

      ! set psinc grid spacing
      dxx = len_a / real(n_lat_pt, kind=DP)

    end subroutine internal_find_dpsinc_in_latvec

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_find_npts_per_ppd_1d(n_pt_in_ppd, &  ! output
         n_pt_in_lat, ext)                                   ! input

      implicit none
      !====================================================================!
      ! This subroutine subdivides the number of psinc grid points along   !
      ! a lattice vector into an odd-integer number of ppds. The odd number!
      ! of grid points in the ppd is returned.                             !
      !--------------------------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 10/11/2004.                    !
      ! Modified by Chris-Kriton Skylaris on 07/08/2005 to add even psinc  !
      ! grid capability.                                                   !
      ! agreco: added flags to use extended NGWFs in one or more directions!
      ! along those directions only one PPD is used (16/10/2014)           !
      !====================================================================!
      integer, intent(out) :: n_pt_in_ppd
      integer, intent(in)  :: n_pt_in_lat
      logical, intent(in)  :: ext


      ! cks: <<local variables>>
      integer :: row
      integer, parameter :: num_fac =11
      ! cks: PPD sizes in order of preference
      integer, parameter :: primes(num_fac) = (/ 5, 6, 7, 8, 9, 10, 11, 12, 4, 3, 13 /)

      ! agreco: only 1 PPD along this direction if
      ! extended NGWFs are required
      if (ext) then
          n_pt_in_ppd = n_pt_in_lat
      else
          n_pt_in_ppd =0
          ppd_npt_finder: do row= 1, num_fac

             if ( mod( n_pt_in_lat, primes(row)) == 0 ) then
                n_pt_in_ppd =primes(row)
                exit ppd_npt_finder
             endif

          enddo ppd_npt_finder
      endif

    end subroutine internal_find_npts_per_ppd_1d

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  end subroutine ppd_strategy_determine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ppd_strategy_check_and_print(fftbox,cell)
    !=============================================================!
    ! This subroutine prints information on the PSINC grids used  !
    ! by ONETEP. It also checks the sizes of the FFTbox and PPDs  !
    ! and stops if they are not sensible.                         !
    !-------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris in 2001 for the !
    ! ONES code.                                                  !
    ! Modified by Chris-Kriton Skylaris on 6/11/2004 to become    !
    ! part of the ppd_strategy module.                            !
    !=============================================================!

    use comms, only : pub_on_root
    use constants, only: DP, stdout, PI, NORMAL, HARTREE_IN_eVs
    use fft_box, only: FFTBOX_INFO
    use rundat, only: pub_output_detail
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_banner

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! cks: internal declarations
    real(kind=DP) :: ke1    ! Kinetic energy cutoff along direction a1
    real(kind=DP) :: ke2    ! Kinetic energy cutoff along direction a2
    real(kind=DP) :: ke3    ! Kinetic energy cutoff along direction a3
    real(kind=DP) :: ke_fac ! Grid space to kinetic energy conversion factor


    ! cks: print PSINC-grids related info
    if (pub_on_root) then

       write(stdout,'(a)') utils_banner('=','PSINC grid sizes')
       write(stdout,'(3(a,i4))') '                      Simulation cell: ',&
            cell%total_pt1 ,' x',cell%total_pt2 ,' x',cell%total_pt3

       write(stdout,'(3(a,i4))') '                              FFT-box: ',&
            fftbox%total_pt1 ,' x',fftbox%total_pt2 ,' x',fftbox%total_pt3

       write(stdout,'(3(a,i4))') '                                  PPD: ',&
            cell%n_pt1 ,' x',cell%n_pt2 ,' x',cell%n_pt3

       if (pub_output_detail >= NORMAL) then
          ! cks: kinetic energy in Eh corresponding to each grid space
          ke_fac =( (6.0_DP*(PI**2))**(2.0_DP/3.0_DP) ) / 2.0_DP
          ke1 = ke_fac/(cell%d1 **2)
          ke2 = ke_fac/(cell%d2 **2)
          ke3 = ke_fac/(cell%d3 **2)
          write(stdout,'(a,f16.12,a, f10.5,a,f10.5,a)') &
               'Grid space d1=',cell%d1,'a0 (KE cutoff=',ke1, &
               'Eh =',ke1*HARTREE_IN_eVs,'eV)'
          write(stdout,'(a,f16.12,a, f10.5,a,f10.5,a)') &
               'Grid space d2=',cell%d2,'a0 (KE cutoff=',ke2, &
               'Eh =',ke2*HARTREE_IN_eVs,'eV)'
          write(stdout,'(a,f16.12,a, f10.5,a,f10.5,a)') &
               'Grid space d3=',cell%d3,'a0 (KE cutoff=',ke3, &
               'Eh =',ke3*HARTREE_IN_eVs,'eV)'
       endif

       write(stdout,'(a)') '================================================================================'

    end if


    ! cks: make sure the size of the fftbox does not exceed the simulation cell
    call internal_check_fftbox_size


    ! cks: make sure that an integer number of grid points and ppds fit
    ! cks: along each lattice vector.
    call internal_check_ppd_shape( &
         cell%a1, cell%d1, cell%n_pt1, 1)
    call internal_check_ppd_shape( &
         cell%a2, cell%d2, cell%n_pt2, 2)
    call internal_check_ppd_shape( &
         cell%a3, cell%d3, cell%n_pt3, 3)



  contains


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_check_fftbox_size
      !================================================================!
      ! This subroutine checks that the size of the fftbox does not    !
      ! exceed the size of the simulation cell.                        !
      !----------------------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 6/11/2004.                 !
      !================================================================!

      use utils, only: utils_abort

      implicit none

      if ( fftbox%total_pt1 > cell%total_pt1 .or. &
           fftbox%total_pt2 > cell%total_pt2 .or. &
           fftbox%total_pt3 > cell%total_pt3) then

         call utils_abort('Error in ppd_strategy_check_and_print: &
                 &one or more of the dimensions of the FFT-box is larger &
                 &than the simulation cell dimensions.')
      end if


      ! cks: make sure that the FFT box has an odd number of points
      ! cks: along each lattice vector direction
     if ( (mod(fftbox%total_pt1, 2) == 0) .or. &
          (mod(fftbox%total_pt2, 2) == 0) .or. &
          (mod(fftbox%total_pt3, 2) == 0) ) then

         call utils_abort('Error in ppd_strategy_check_and_print: &
                 &the number of psinc grid points along one or more of the &
                 &lattice vectors of the FFT-box is even.')
      end if

    end subroutine internal_check_fftbox_size

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_check_ppd_shape(lat_vec, dx, ppd_npt, &  ! input
         direction_counter)                                      ! input
      !================================================================!
      ! This subroutine checks in a lattice vector direction that the  !
      ! similation cell grid contains an integer number of grid points !
      ! and an integer number of ppds.                                 !
      !----------------------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 21/6/2001.                 !
      !================================================================!

      use geometry, only: point, geometry_magnitude
      use constants, only: DP
      use utils, only: utils_assert

      implicit none

      type(POINT), intent(in) :: lat_vec
      real(kind=DP), intent(in) :: dx
      integer, intent(in) :: ppd_npt, direction_counter

      ! cks: internal declarations
      real(kind=DP) :: total_npt_real
      integer :: total_npt

      call utils_assert(ppd_npt >= 0, 'Error in internal_check_ppd_shape: &
           &the number of points per ppd in at least one direction is negative.&
           & The direction is ', direction_counter)

      total_npt_real = geometry_magnitude(lat_vec) / dx
      total_npt = nint(total_npt_real)

      ! cks: make sure that an essentially integral number of grid points fits
      ! cks: in the current edge of the simulation cell
      call utils_assert(abs(total_npt_real-real(total_npt,kind=DP))<=1.0e-5_DP,&
           'Error in internal_check_ppd_shape: the number of points in &
           &at least one direction of the simulation cell is not an integer &
           & -- adjust grid spacing. The direction and number of points were ',&
           direction_counter, total_npt)

      ! cks: make sure that an integer number of ppds fits in the current edge of
      !      the simulation cell.
      call utils_assert(mod(total_npt, ppd_npt) == 0, &
           'Error in internal_check_ppd_shape: the number of ppds in &
           &at least one direction of the simulation cell is not an integer &
           & -- adjust points in ppd. The direction and number of ppds were ',&
           real(direction_counter,kind=DP), &
           real(total_npt,kind=DP) / real(ppd_npt,kind=DP))

    end subroutine internal_check_ppd_shape

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  end subroutine ppd_strategy_check_and_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module ppd_strategy




