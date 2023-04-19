! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   The subroutines in this file were written by                              !
!                                                                             !
!   Chris-Kriton Skylaris, Arash A. Mostofi, Peter D. Haynes and              !
!   Nicholas D.M. Hine                                                        !
!                                                                             !
!   TCM Group, Cavendish laboratory                                           !
!   Madingley Road                                                            !
!   Cambridge CB3 0HE                                                         !
!   UK                                                                        !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module basis

  use constants, only: DP
  use geometry, only: POINT

  implicit none

  private

  ! This structure holds all the information for the localised fft grid
  ! (contained in the "tight box") that will be used for each NGWF.
  ! It also holds all the information necessary for extracting a NGWF
  ! function from the ppd-in-standard-grid representation and placing it in
  ! the localised fft grid and vice-versa.
  type FUNCTION_TIGHT_BOX

     ! the number of zero padding points in each lattice vector direction
     integer :: pad1,pad2,pad3

     ! the very first and very last ppd number in each lattice vector
     ! direction that contains NGWF values.
     integer :: start_ppds1,start_ppds2,start_ppds3
     integer :: finish_ppds1,finish_ppds2,finish_ppds3

     ! the very first and the very last NGWF value point on the
     ! very first and very last ppd with function values in each lattice
     ! vector direction.
     integer :: start_pts1,start_pts2,start_pts3
     integer :: finish_pts1,finish_pts2,finish_pts3

     ! the number of tight points in the localised tight box grid in each
     ! lattice vector direction.
     integer :: tight_pts1, tight_pts2, tight_pts3

  end type FUNCTION_TIGHT_BOX


  type SPHERE

     ! the cartesian coordinates of the atom to which sphere belongs in a.u.
     type(POINT)       :: centre

     ! the radius of the sphere in atomic units
     ! cks: For NGWFs with halos, note that this is equal to the
     ! cks: NGWF radius (without the halo) and therefore it is different
     ! cks: from element%radius which includes the halo
     real(kind=DP) :: radius

     ! the number of ppds belonging to the sphere
     integer :: n_ppds_sphere

     ! ppd_list(1,:) contains a list of the indices (in the global ppd
     ! counting scheme) of the ppds that belong to the sphere
     ! ppd_list(2,:) contains a list of integers from -13 to 13 for each
     ! ppd belonging to the sphere. This number shows if the contribution
     ! to the sphere comes from the simulation cell or from a periodic image
     ! of the sphere in one of the neighbouring simulation cells
     integer, allocatable, dimension(:,:) :: ppd_list

     ! This is an offset showing the starting index (for the function inside
     ! the sphere) in the large array where all the functions are stored.
     integer :: offset

     ! agreco: flag to specify if the NGWF is extended along a given direction
     ! currently this is set to be equal to pub_extend_ngwf, in the future we
     ! may want to use this to allow different functions to have different
     ! localisation properties within the same simulation
     logical :: extended(3)

     ! =========================================================================
     ! jd: Any additions or changes to this datatype should be reflected
     !     (at least) in remote_mod (search for [*]) and in
     !     ngwfs_mod::internal_send_ngwf() and internal_recv_ngwf().
     ! =========================================================================

  end type SPHERE


  public :: SPHERE, FUNCTION_TIGHT_BOX

  ! ndmh: Location routines for finding points/vectors etc in cells & boxes
  public :: basis_function_origin_wrt_tb
  public :: basis_func_centre_wrt_fftbox
  public :: basis_find_function_wrt_box
  public :: basis_ket_start_wrt_fftbox
  public :: basis_location_func_wrt_cell
  public :: basis_point_wrt_box
  public :: basis_start_of_box_wrt_cell
  public :: basis_box_origin_to_atom

  ! ndmh: Sphere initialisation routines
  public :: basis_initialise_sphere
  public :: basis_sphere_deallocate
  public :: basis_ppd_location
  public :: basis_copy_sphere

  ! ndmh: Function-in-box manipulation routines
  public :: basis_copy_function_to_box
  public :: basis_add_function_to_box
  public :: basis_extract_function_from_box
  public :: basis_dot_function_with_box
  public :: basis_multiply_function_by_box
  public :: basis_clean_function
  public :: basis_put_tightbox_in_fftbox
  public :: basis_copy_tightbox_to_fftbox
  public :: basis_phase_on_fftbox_recip

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_function_origin_wrt_tb(tb_orig1, tb_orig2, tb_orig3,  &
       orig_wrt_cell, tight_box, cell)

    !=======================================================================!
    ! This subroutine returns the coordinates of the atom centre of an NGWF !
    ! from the origin of its tightbox expressed in terms of                 !
    ! (in general non-integer) numbers of standard grid points              !
    ! along each lattice vector direction.                                  !
    !-----------------------------------------------------------------------!
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 12/4/2004                         !
    !=======================================================================!

    use constants, only: DP, two_pi
    use geometry, only : POINT, operator(.DOT.)
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: tb_orig1
    real(kind=DP), intent(out) :: tb_orig2
    real(kind=DP), intent(out) :: tb_orig3
    type(POINT), intent(in) :: orig_wrt_cell
    type(FUNCTION_TIGHT_BOX), intent(in) :: tight_box
    type(CELL_INFO), intent(in) :: cell


    ! cks: << local variables>>
    integer :: start1
    integer :: start2
    integer :: start3
    real(kind=DP) :: ff
    real(kind=DP) :: inv_two_pi

    inv_two_pi = 1.0_DP / TWO_PI

    ! cks: origin of tightbox wrt cell in terms of integer numbers of grid
    !      points

    start1 = (tight_box%start_ppds1-1)*cell%n_pt1 + tight_box%start_pts1 - 1
    start2 = (tight_box%start_ppds2-1)*cell%n_pt2 + tight_box%start_pts2 - 1
    start3 = (tight_box%start_ppds3-1)*cell%n_pt3 + tight_box%start_pts3 - 1

    ! cks: origin of function wrt tightbox in terms of fractional numbers of
    !      grid points
    ff = (orig_wrt_cell .DOT. cell%b1)
    tb_orig1 = ff*inv_two_pi*real(cell%total_pt1, kind=DP) - &
         real(start1, kind=DP)

    ff = (orig_wrt_cell .DOT. cell%b2)
    tb_orig2 = ff*inv_two_pi*real(cell%total_pt2, kind=DP) - &
         real(start2, kind=DP)

    ff = (orig_wrt_cell .DOT. cell%b3)
    tb_orig3 = ff*inv_two_pi*real(cell%total_pt3, kind=DP) - &
         real(start3, kind=DP)

  end subroutine basis_function_origin_wrt_tb


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function basis_func_centre_wrt_fftbox(ngwf_centre, &
       ngwf_start1, ngwf_start2, ngwf_start3, &
       ngwf_cell_start1, ngwf_cell_start2, ngwf_cell_start3, cell)

    !=============================================================!
    ! Returns a vector from the origin of an FFTbox to the centre !
    ! of a function.                                              !
    !-------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 2/5/2003.               !
    ! Modified by Quintin Hill on 10/09/2008 to use cell      !
    !=============================================================!

    use geometry, only: point, operator(+), operator(*)
    use simulation_cell, only: CELL_INFO

    implicit none

    type(POINT) :: basis_func_centre_wrt_fftbox

    ! Arguments
    type(POINT), intent(in) :: ngwf_centre
    integer, intent(in) :: ngwf_start1, ngwf_start2, ngwf_start3
    integer, intent(in) :: ngwf_cell_start1, ngwf_cell_start2, ngwf_cell_start3
    type(CELL_INFO), intent(in) :: cell

    ! cks: internal
    type(POINT) :: buffvec

    buffvec = ngwf_centre + &
         real((ngwf_start1 -ngwf_cell_start1), kind=DP)*&
         &cell%d1*cell%a1_unit + &
         real((ngwf_start2 -ngwf_cell_start2), kind=DP)*&
         &cell%d2*cell%a2_unit + &
         real((ngwf_start3 -ngwf_cell_start3), kind=DP)*&
         &cell%d3*cell%a3_unit

    basis_func_centre_wrt_fftbox = buffvec

  end function basis_func_centre_wrt_fftbox


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_find_function_wrt_box( &
       row_start1, row_start2, row_start3, &           ! output
       box_start1, box_start2, box_start3, tight_box, cell, fftbox)  ! input

    !============================================================!
    ! This subroutine returns the starting gridpoint, in each    !
    ! lattice vector direction, of the current function with     !
    ! respect to the start of the FFTbox.                        !
    !------------------------------------------------------------!
    ! Written by Arash A. Mostofi, September 2002.               !
    ! Modified by Chris-Kriton Skylaris on 11/08/2007 to work in !
    ! the case where the FFT box coincides with the simulation   !
    ! cell.                                                      !
    ! Moved to basis_mod by Nicholas Hine, 29/06/2009.           !
    !============================================================!

    use fft_box, only: FFTBOX_INFO
    use simulation_cell, only: CELL_INFO

    implicit none

    type(FUNCTION_TIGHT_BOX), intent(in) :: tight_box
    integer, intent(out) :: row_start1,row_start2,row_start3
    integer, intent(in)  :: box_start1,box_start2,box_start3
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox


    if (fftbox%coin1) then
       row_start1 = 1
    else
       row_start1 = (tight_box%start_ppds1 - 1)*cell%n_pt1 &
            + tight_box%start_pts1 - box_start1 + 1
       if (row_start1 .lt. 1) then
          row_start1 = row_start1 + cell%total_pt1
       elseif (row_start1 .gt. cell%total_pt1) then
          row_start1 = row_start1 - cell%total_pt1
       endif
    endif

    if (fftbox%coin2) then
       row_start2 = 1
    else
       row_start2 = (tight_box%start_ppds2 - 1)*cell%n_pt2 &
            + tight_box%start_pts2 - box_start2 + 1
       if (row_start2 .lt. 1) then
          row_start2 = row_start2 + cell%total_pt2
       elseif (row_start2 .gt. cell%total_pt2) then
          row_start2 = row_start2 - cell%total_pt2
       endif
    endif

    if (fftbox%coin3) then
       row_start3 = 1
    else
       row_start3 = (tight_box%start_ppds3 - 1)*cell%n_pt3 &
            + tight_box%start_pts3 - box_start3 + 1
       if (row_start3 .lt. 1) then
          row_start3 = row_start3 + cell%total_pt3
       elseif (row_start3 .gt. cell%total_pt3) then
          row_start3 = row_start3 - cell%total_pt3
       endif
    endif

  end subroutine basis_find_function_wrt_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_ket_start_wrt_fftbox(row_start1,row_start2,row_start3, &
       n1, n2, n3, fftbox)

    !===============================================================!
    ! Returns the starting position of the tightbox of the          !
    ! 'ket' function with respect to the beginning of the FFT box.  !
    !---------------------------------------------------------------!
    ! This subroutine also works when the length of the FFT box     !
    ! coincides with that of the simulation cell along one of more  !
    ! lattice vector directions.                                    !
    !---------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 6/4/2007, based on the    !
    ! basis_location_fa_wrt_box subroutine by Arash Mostofi.        !
    !===============================================================!

    use fft_box, only: FFTBOX_INFO
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    integer, intent(out) :: row_start1
    integer, intent(out) :: row_start2
    integer, intent(out) :: row_start3
    integer, intent(in)  :: n1,n2,n3
    type(FFTBOX_INFO), intent(in) :: fftbox


    if (fftbox%coin1) then
       row_start1 =1
    else
       row_start1 =int(n1/3) + 1
    endif


    if (fftbox%coin2) then
       row_start2 =1
    else
       row_start2 = int(n2/3) + 1
    endif


    if (fftbox%coin3) then
       row_start3 =1
    else
       row_start3 = int(n3/3) + 1
    endif


  end subroutine basis_ket_start_wrt_fftbox



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_location_func_wrt_cell( &
       start1,start2,start3,tight_box,cell)

    !==========================================================!
    ! Returns the position of a function wrt simulation cell.  !
    !----------------------------------------------------------!
    ! Written by Arash A. Mostofi, September 2002.             !
    !==========================================================!

    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    type(FUNCTION_TIGHT_BOX), intent(in) :: tight_box
    integer, intent(out)               :: start1
    integer, intent(out)               :: start2
    integer, intent(out)               :: start3
    type(CELL_INFO), intent(in) :: cell

    start1 = (tight_box%start_ppds1 - 1)*cell%n_pt1 + tight_box%start_pts1
    start2 = (tight_box%start_ppds2 - 1)*cell%n_pt2 + tight_box%start_pts2
    start3 = (tight_box%start_ppds3 - 1)*cell%n_pt3 + tight_box%start_pts3

  end subroutine basis_location_func_wrt_cell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_point_wrt_box(dg1, dg2, dg3, &
       box_centre_point, box_pts1, box_pts2, box_pts3, cell, fftbox)

    !===================================================================!
    ! This subroutine returns the distance along each lattice vector    !
    ! direction, in gridspacings, from the origin of the FFT box to     !
    ! a given point. The point is positioned as close to the middle of  !
    ! the box as possible                                               !
    !-------------------------------------------------------------------!
    ! Written by Arash A. Mostofi, September 2002, originally named     !
    ! pseudo_centre_of_proj_wrt_box.                                    !
    ! Updated by Arash A. Mostofi on  10/07/2003 for non-orthogonal     !
    ! lattice geometries.                                               !
    ! Modified by Chris-Kriton Skylaris on 11/08/2007 to work in the    !
    ! case where the FFTbox coincides with the simulation cell.         !
    ! Moved to basis_mod by Nicholas Hine, 29/06/2009                   !
    !===================================================================!

    use constants, only: PI
    use fft_box, only: FFTBOX_INFO
    use geometry, only: POINT, OPERATOR(.dot.)
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    type(POINT), intent(in)     :: box_centre_point
    integer, intent(in)         :: box_pts1,box_pts2,box_pts3
    real(kind=DP), intent(out)  :: dg1,dg2,dg3 !{D}istance in {G}rid-spacings
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! internal declarations
    real(kind=DP), parameter :: inv_two_pi = 0.5_DP/PI
    real(kind=DP) :: dd,xx
    real(kind=DP) :: tp1,tp2,tp3

    tp1 = real(cell%total_pt1,kind=DP)
    tp2 = real(cell%total_pt2,kind=DP)
    tp3 = real(cell%total_pt3,kind=DP)

    if (fftbox%coin1) then
       dg1=tp1*inv_two_pi*(box_centre_point.DOT.cell%b1)
    else
       dd=real(box_pts1/2-1,kind=DP)
       xx=MOD(tp1*inv_two_pi*(box_centre_point.DOT.cell%b1),1.0_DP)
       dg1 = dd + xx
    endif

    if (fftbox%coin2) then
       dg2=tp2*inv_two_pi*(box_centre_point.DOT.cell%b2)
    else
       dd=real(box_pts2/2-1,kind=DP)
       xx=MOD(tp2*inv_two_pi*(box_centre_point.DOT.cell%b2),1.0_DP)
       dg2 = dd + xx
    endif

    if (fftbox%coin3) then
       dg3=tp3*inv_two_pi*(box_centre_point.DOT.cell%b3)
    else
       dd=real(box_pts3/2-1,kind=DP)
       xx=MOD(tp3*inv_two_pi*(box_centre_point.DOT.cell%b3),1.0_DP)
       dg3 = dd + xx
    endif

  end subroutine basis_point_wrt_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_start_of_box_wrt_cell( &
       box_start1, box_start2, box_start3, &         ! output
       box_centre_point, dg1, dg2, dg3, cell)       ! input

    !=================================================================!
    ! This subroutine returns the starting gridpoint in each lattice  !
    ! vector direction of the FFTbox in the simulation cell.          !
    !-----------------------------------------------------------------!
    ! Written by Arash A. Mostofi, September 2002.                    !
    ! Updated by Arash A. Mostofi on  10/07/2003 for non-orthogonal   !
    ! lattice geometries.                                             !
    ! Moved to basis_mod by Nicholas Hine, 29/06/2009                 !
    !=================================================================!

    use constants, only: PI
    use geometry, only: POINT, OPERATOR(.dot.)
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    type(POINT), intent(in) :: box_centre_point
    integer, intent(out) :: box_start1,box_start2,box_start3
    ! !{D}istance in {G}rid-spacings of box centre point from origin
    ! of box in each *lattice vector direction*. Not necessarily integer values
    real(kind=DP), intent(in) :: dg1,dg2,dg3
    type(CELL_INFO), intent(in) :: cell

    ! internal declarations
    real(kind=DP), parameter :: inv_two_pi = 0.5_DP/PI
    real(kind=DP) :: aa
    real(kind=DP) :: tp1,tp2,tp3

    tp1=real(cell%total_pt1,kind=DP)
    tp2=real(cell%total_pt2,kind=DP)
    tp3=real(cell%total_pt3,kind=DP)

    aa=tp1*inv_two_pi*(box_centre_point.DOT.cell%b1) - dg1
    box_start1 = NINT(aa) + 1

    aa=tp2*inv_two_pi*(box_centre_point.DOT.cell%b2) - dg2
    box_start2 = NINT(aa) + 1

    aa=tp3*inv_two_pi*(box_centre_point.DOT.cell%b3) - dg3
    box_start3 = NINT(aa) + 1

  end subroutine basis_start_of_box_wrt_cell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_box_origin_to_atom(box_to_atom,box_origin,atom_origin, &
       cell_start1,cell_start2,cell_start3,da1,da2,da3,cell)

    use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
         OPERATOR(*)
    use simulation_cell, only: CELL_INFO

    ! Arguments
    type(POINT),intent(out) :: box_origin
    type(POINT),intent(out) :: box_to_atom
    type(POINT),intent(in) :: atom_origin
    integer,intent(in) :: cell_start1,cell_start2,cell_start3
    type(POINT),intent(in) :: da1, da2, da3
    type(CELL_INFO), intent(in) :: cell

    ! Find vector to origin of box
    box_origin = (cell_start1-1)*da1 + (cell_start2-1)*da2 + &
         (cell_start3-1)*da3

    ! ndmh: vector from origin of box to centre of atom
    box_to_atom = atom_origin - box_origin

    ! ndmh: check components of this vector in terms of lattice vectors
    ! ndmh: are all positive. If not, the box has been looped back into the
    ! ndmh: cell, so shift the origin back accordingly
    if ((box_to_atom.DOT.cell%b1) < 0.0_DP) then
       box_origin = box_origin - cell%a1
       box_to_atom = box_to_atom + cell%a1
    end if
    if ((box_to_atom.DOT.cell%b2) < 0.0_DP) then
       box_origin = box_origin - cell%a2
       box_to_atom = box_to_atom + cell%a2
    end if
    if ((box_to_atom.DOT.cell%b3) < 0.0_DP) then
       box_origin = box_origin - cell%a3
       box_to_atom = box_to_atom + cell%a3
    end if

  end subroutine basis_box_origin_to_atom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_phase_on_fftbox_recip(f_out, &
       nn1,nn2,nn3,ld1,ld2,len1,len2,len3,f_in)

    !====================================================================!
    ! Application of any desired phase factor to a quantity (f_in)       !
    ! represented on the reciprocal space coarse grid pair box.          !
    ! The phase factor is characterised by the three quantities          !
    ! len1, len2 and len3, which give the phase vector in terms of       !
    ! numbers of gridpoints (not necessarily integer) along each lattice !
    ! direction. n1, n2 and n3 can be even or odd.                       !
    !--------------------------------------------------------------------!
    ! Written by Arash A. Mostofi, August 2002.                          !
    ! Modified January 2003 - got rid of make_even.                      !
    ! Modified April 2003 - complex-to-complex FFTs throughout.          !
    ! Modified April 2008 and August 2008 by Nick Hine for speed         !
    !====================================================================!

    use constants, only: DP, two_pi, cmplx_i

    implicit none

    ! Arguments
    integer, intent(in)   :: nn1,nn2,nn3,ld1,ld2
    real(kind=DP), intent(in)  :: len1,len2,len3
    complex(kind=DP), intent(inout)  :: f_out(ld1,ld2,nn3)
    complex(kind=DP), intent(in),optional  :: f_in(ld1,ld2,nn3)

    ! Local Variables
    complex(kind=DP) :: phase1,phase2,phase3
    complex(kind=DP) :: phase_neg1,phase_neg2,phase_neg3
    complex(kind=DP) :: z1(nn1)
    complex(kind=DP) :: z2(nn2)
    complex(kind=DP) :: z3(nn3)
    complex(kind=DP) :: z2z3
    integer  :: i1,i2,i3
    real(kind=DP) :: r_n1, r_n2, r_n3

    r_n1 = real(nn1,kind=DP)
    r_n2 = real(nn2,kind=DP)
    r_n3 = real(nn3,kind=DP)

    ! the elementary phase factors for the three directions in space
    phase1 = exp(two_pi * cmplx_i * len1 / r_n1)
    phase2 = exp(two_pi * cmplx_i * len2 / r_n2)
    phase3 = exp(two_pi * cmplx_i * len3 / r_n3)

    ! for the negative reciprocal lattice points:
    ! e^(-i*2*pi*len2) and e^(-i*2*pi*len3)
    phase_neg1 = exp(-two_pi * cmplx_i * len1)
    phase_neg2 = exp(-two_pi * cmplx_i * len2)
    phase_neg3 = exp(-two_pi * cmplx_i * len3)

    ! initialise phase factor arrays
    z1(1) = cmplx(1.0_DP,0.0_DP,kind=DP)
    do i1=2,nn1/2+1
       z1(i1) = z1(i1-1)*phase1
    end do
    z1(nn1/2+2) = z1(nn1/2+1)*phase_neg1*phase1
    do i1=nn1/2+3,nn1
       z1(i1) = z1(i1-1)*phase1
    end do

    z2(1) = cmplx(1.0_DP,0.0_DP,kind=DP)
    do i2=2,nn2/2+1
       z2(i2) = z2(i2-1)*phase2
    end do
    z2(nn2/2+2) = z2(nn2/2+1)*phase_neg2*phase2
    do i2=nn2/2+3,nn2
       z2(i2) = z2(i2-1)*phase2
    end do

    z3(1) = cmplx(1.0_DP,0.0_DP,kind=DP)
    do i3=2,nn3/2+1
       z3(i3) = z3(i3-1)*phase3
    end do
    z3(nn3/2+2) = z3(nn3/2+1)*phase_neg3*phase3
    do i3=nn3/2+3,nn3
       z3(i3) = z3(i3-1)*phase3
    end do

    ! Now copy over the f_in array while multiplying by phase
    ! (or multiply f_out if f_in not supplied)
    if (present(f_in)) then
       do i3=1,nn3
          do i2=1,nn2
             z2z3 = z2(i2) * z3(i3)
             f_out(1:nn1,i2,i3) = z1(1:nn1)*z2z3*f_in(1:nn1,i2,i3)
          end do
       end do
    else
       do i3=1,nn3
          do i2=1,nn2
             z2z3 = z2(i2) * z3(i3)
             f_out(1:nn1,i2,i3) = z1(1:nn1)*z2z3*f_out(1:nn1,i2,i3)
          end do
       end do
    end if

  end subroutine basis_phase_on_fftbox_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_copy_sphere(current_sphere,old_sphere,offset, extended)

    !=======================================================================!
    ! This subroutine initialises sphere_a by copying its values from       !
    ! sphere_b                                                              !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! sphere_a       (output)                                               !
    ! sphere_b       (input)                                                !
    ! offset         (input)
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine, 16/09/2008                                  !
    !=======================================================================!

    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(SPHERE), intent(out) :: current_sphere ! sphere to be initialised
    type(SPHERE), intent(in) :: old_sphere      ! sphere to copy from
    integer, intent(in) :: offset               ! current offset of NGWF data
    ! agreco: override extended property of the sphere to copy from
    ! this allows to specify a different localisation for another
    ! function on the same atom
    logical, intent(in), optional :: extended(3)

    ! Local Variables
    integer :: ierr             ! allocation error flag
    integer :: ippd             ! loop counter

    current_sphere%centre%x = old_sphere%centre%x
    current_sphere%centre%y = old_sphere%centre%y
    current_sphere%centre%z = old_sphere%centre%z
    current_sphere%radius = old_sphere%radius
    current_sphere%n_ppds_sphere = old_sphere%n_ppds_sphere
    current_sphere%offset = offset

    ! agreco: override extended property if required
    ! otherwise copy extended value from origin sphere
    if (present(extended)) then
       current_sphere%extended = extended
    else
       current_sphere%extended = old_sphere%extended
    end if

    allocate(current_sphere%ppd_list(2,current_sphere%n_ppds_sphere),stat=ierr)
    call utils_alloc_check('basis_copy_sphere','current_sphere%ppd_list',ierr)

    do ippd=1,current_sphere%n_ppds_sphere
       current_sphere%ppd_list(:,ippd) = old_sphere%ppd_list(:,ippd)
    end do

  end subroutine basis_copy_sphere


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_initialise_sphere(current_sphere,centre,radius,&
       current_offset,ppd_list,ppd_loc,cell,fftbox, extended)

    !=======================================================================!
    ! This subroutine initialises a single sphere for one function on       !
    ! one atom.                                                             !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! centre         (input)                                                !
    ! radius         (input)                                                !
    ! cell           (input)                                                !
    ! current_offset (input)                                                !
    ! current_sphere (output)                                               !
    !-----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.                             !
    ! Updated to work when the FFT box coincides with the simulation cell   !
    ! by Chris-Kriton Skylaris on 131/08/2007.                              !
    ! Modified by Nicholas Hine on 24/06/2009 to combine ppd_loc            !
    ! with ppd_list for MPI efficiency.                                     !
    ! Modified to use function basis_ppd_location rather than ppd_location  !
    ! array by Nicholas Hine, October 2009.                                 !
    !=======================================================================!

    use constants, only: DP, TWO_PI
    use fft_box, only: FFTBOX_INFO
    use geometry, only: point, operator(+), operator(*), geometry_distance, &
         operator(.DOT.), magnitude
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_abort

    implicit none

    ! Arguments
    type(SPHERE), intent(out) :: current_sphere ! sphere to be initialised
    type(POINT), intent(in) :: centre           ! centre of sphere
    real(kind=DP), intent(in) :: radius         ! radius of sphere
    integer, intent(in) :: current_offset       ! psinc-counting offset for spere (local proc)
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    integer, intent(inout) :: ppd_list(cell%n_ppds) ! temporary ppd index list
    integer, intent(inout) :: ppd_loc(cell%n_ppds)  ! temporary ppd loc list
    ! agreco: optional argument to specify extension of the NGWF for this sphere
    logical, intent(in), optional :: extended(3)

    ! cks: internal variables
    integer :: ierr             ! allocation error flag
    integer :: a3_step          ! ppd counter along lattice vector a3
    integer :: a2_step          ! ppd counter along lattice vector a2
    integer :: a1_step          ! ppd counter along lattice vector a1
    integer :: a1_neighbour     ! neighbour simcell counter along a1
    integer :: a2_neighbour     ! neighbour simcell counter along a2
    integer :: a3_neighbour     ! neighbour simcell counter along a3
    integer :: ppd              ! counter of all ppds in sim cell
    integer :: one_or_zero      ! ppd belongs to sphere flag
    type(POINT) :: periodic_centre  ! centre od sphere in neighbour of simc ell
    logical :: ppd_near_sphere  ! result of quick test to select ppds near sphere

    real(kind=DP) :: a1_centre  ! ppd start along a1 (fractional coords)
    real(kind=DP) :: a2_centre  ! ppd start along a2 (fractional coords)
    real(kind=DP) :: a3_centre  ! ppd start along a3 (fractional coords)
    real(kind=DP) :: a1_radius  ! sphere radius expressed in ppds along a1
    real(kind=DP) :: a2_radius  ! sphere radius expressed in ppds along a2
    real(kind=DP) :: a3_radius  ! sphere radius expressed in ppds along a3
    integer :: a1_start         ! ppd start along a1 (periodic)
    integer :: a2_start         ! ppd start along a2 (periodic)
    integer :: a3_start         ! ppd start along a3 (periodic)
    integer :: a1_end           ! ppd finish along a1 (periodic)
    integer :: a2_end           ! ppd finish along a2 (periodic)
    integer :: a3_end           ! ppd finish along a3 (periodic)
    integer :: a1_ppds, a2_ppds, a3_ppds ! shorthand variables
    integer :: ppd_count        ! counter for ppd_list and ppd_loc
    integer :: ppd_check        ! counter for avoiding double inclusion
    integer :: ppd_lowest
    logical :: ppd_already_found

    ! cks: initialise radius of sphere
    current_sphere%radius = radius

    ! cks: initialise the centre of the sphere
    current_sphere%centre = centre

    ! cks: initialise the offset
    current_sphere%offset = current_offset

    ! agreco: initialise the extended property
    ! if not specified, initialise to fully localised spheres
    if (present(extended)) then
       current_sphere%extended = extended
    else
       current_sphere%extended = (/ .false., .false., .false. /)
    end if

    a1_ppds = cell%n_ppds_a1
    a2_ppds = cell%n_ppds_a2
    a3_ppds = cell%n_ppds_a3

    ! find centre in ppd coordinates along a1
    a1_centre = (centre .dot. cell%b1) * real(a1_ppds,kind=DP) / TWO_PI
    ! find radius projected onto a1 axis in terms of a1 ppd lengths
    a1_radius = radius * magnitude(cell%b1) * real(a1_ppds,kind=DP) / TWO_PI
    ! find range of ppds to check along a1
    a1_start = floor(a1_centre - a1_radius + 1_DP)
    a1_end =   floor(a1_centre + a1_radius + 1_DP)
    ! ensure we are not checking last slab twice
    if ((a1_radius>=1_DP).and.(modulo(a1_end + a1_ppds - 1, a1_ppds) == &
         modulo(a1_start + a1_ppds - 1, a1_ppds))) a1_end = a1_end - 1

    ! find centre in ppd coordinates along a2
    a2_centre = (centre .dot. cell%b2) * real(a2_ppds,kind=DP) / TWO_PI
    a2_radius = radius * magnitude(cell%b2) * real(a2_ppds,kind=DP) / TWO_PI
    a2_start = floor(a2_centre - a2_radius + 1_DP)
    a2_end =   floor(a2_centre + a2_radius + 1_DP)
    if ((a2_radius>=1_DP).and.(modulo(a2_end + a2_ppds - 1, a2_ppds) == &
         modulo(a2_start + a2_ppds - 1, a2_ppds))) a2_end = a2_end - 1

    ! find centre in ppd coordinates along a3
    a3_centre = (centre .dot. cell%b3) * real(a3_ppds,kind=DP) / TWO_PI
    a3_radius = radius * magnitude(cell%b3) * real(a3_ppds,kind=DP) / TWO_PI
    a3_start = floor(a3_centre - a3_radius + 1_DP)
    a3_end =   floor(a3_centre + a3_radius + 1_DP)
    if ((a3_radius>=1_DP).and.(modulo(a3_end + a3_ppds - 1, a3_ppds) == &
         modulo(a3_start + a3_ppds - 1, a3_ppds))) a3_end = a3_end - 1

    ppd_count = 0
    ! cks: loop over all ppds in simulation cell first
    do a3_step=a3_start,a3_end
       do a2_step=a2_start,a2_end
          do a1_step=a1_start,a1_end

             ppd = (modulo(a3_step + a3_ppds - 1, a3_ppds)*a2_ppds &
                  + modulo(a2_step + a2_ppds - 1, a2_ppds))*a1_ppds &
                  + modulo(a1_step + a1_ppds - 1, a1_ppds) + 1

             do a1_neighbour = -1,1
                do a2_neighbour = -1,1
                   do a3_neighbour = -1,1

                      ! this is the centre of the NGWF_sphere in one of the
                      ! neighbouring simulation cells
                      periodic_centre = current_sphere%centre + &
                           real(a1_neighbour,kind=DP)*cell%a1 + &
                           real(a2_neighbour,kind=DP)*cell%a2 + &
                           real(a3_neighbour,kind=DP)*cell%a3

                      ! consider if a ppd belongs to the current NGWF
                      ! function (or its periodic image in one of the
                      ! surrounding simulation cells) only if the distance of
                      ! one of its vertices from the NGWF centre is
                      ! smaller than the sum of the NGWF radius and the
                      ! length of the three ppd edges. This test selects, for
                      ! all the spheres, a number of ppds that scales linearly
                      ! with system size.
                      ppd_near_sphere = (geometry_distance(periodic_centre, &
                           basis_ppd_location(ppd,cell)) <= current_sphere%radius + &
                           real(cell%n_pt1,kind=DP)*cell%d1 + &
                           real(cell%n_pt2,kind=DP)*cell%d2 + &
                           real(cell%n_pt3,kind=DP)*cell%d3)

                      if (ppd_near_sphere) then

                         one_or_zero = basis_ppd_belongs_to_sphere( &
                              periodic_centre,basis_ppd_location(ppd,cell), &
                              current_sphere%radius,cell)

                         if (one_or_zero == 1) then

                            ppd_already_found = .false.
                            do ppd_check = 1,ppd_count
                               if (ppd_list(ppd_check)==ppd) then
                                  ppd_already_found = .true.
                                  exit
                               end if
                            end do
                            if (ppd_already_found .and. &
                                 ((.not.fftbox%coin1) .and. &
                                 (.not.fftbox%coin2) .and. &
                                 (.not.fftbox%coin3))) then
                               call utils_abort('Error in basis_initialise_&
                                    &sphere: an NGWF sphere and its periodic &
                                    &image have values on ppd ', ppd)
                            end if

                            if (ppd_already_found) then
                               ppd_loc(ppd_check) = (9*a1_neighbour+ &
                                    3*a2_neighbour+a3_neighbour) * &
                                    one_or_zero - 14*(1-one_or_zero)
                            end if
                            if (.not.ppd_already_found) then
                               ppd_count = ppd_count + 1
                               ppd_list(ppd_count) = ppd
                               ppd_loc(ppd_count) = (9*a1_neighbour+ &
                                    3*a2_neighbour+a3_neighbour) * &
                                    one_or_zero - 14*(1-one_or_zero)
                            end if

                         end if

                      end if

                   end do
                end do
             end do

          end do
       end do
    end do

    ! cks: now determine the other quantities relevant to the sphere
    current_sphere%n_ppds_sphere = ppd_count

    allocate(current_sphere%ppd_list(2,current_sphere%n_ppds_sphere),stat=ierr)
    call utils_alloc_check('basis_initialise_sphere', &
         'current_sphere%ppd_list',ierr)

    ! ndmh: sort the ppds in ascending order
    do ppd_check=1,ppd_count
       ppd_lowest = 1
       do ppd=1,ppd_count
          if (ppd_list(ppd)<ppd_list(ppd_lowest)) ppd_lowest=ppd
       end do
       current_sphere%ppd_list(1,ppd_check) = ppd_list(ppd_lowest)
       current_sphere%ppd_list(2,ppd_check) = ppd_loc(ppd_lowest)
       ppd_list(ppd_lowest) = cell%n_ppds + 1
    end do

  end subroutine basis_initialise_sphere


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  integer function basis_ppd_belongs_to_sphere(centre,ppd_start, &
       NGWF_radius,cell)

    !===================================================================!
    ! This function returns the value 1 if the current ppd belongs to   !
    ! the current sphere and zero otherwise.                            !
    !-------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.                         !
    ! Updated by Chris-Kriton Skylaris on 13/08/2007.                   !
    !===================================================================!

    use constants, only: DP, TWO_PI
    use geometry, only: point, operator(.DOT.), geometry_magnitude, &
         operator(-), local_displacement, operator(+), geometry_distance
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: NGWF_radius
    TYPE(POINT), intent(in) :: centre,ppd_start
    type(CELL_INFO), intent(in) :: cell

    ! cks: internal variable declarations
    real(kind=DP) :: first_plane_a1_a2, second_plane_a1_a2
    real(kind=DP) :: first_plane_a2_a3, second_plane_a2_a3
    real(kind=DP) :: first_plane_a3_a1, second_plane_a3_a1
    integer :: more_than_two, loc_1, loc_2, loc_3
    real(kind=DP) :: local_1, local_2, local_3, length
    TYPE(POINT) :: current_displacement, current_point

    ! ======================= SPHERE CENTRE INSIDE PPD? ======================

    basis_ppd_belongs_to_sphere = 0

    ! cks: find the distance of the sphere from the a1_a2 planes defined by
    !      boundary points of the current ppd
    first_plane_a1_a2 = ((ppd_start - centre) .DOT. cell%b3) / &
         geometry_MAGNITUDE(cell%b3)

    second_plane_a1_a2 = (((ppd_start - centre) .DOT. cell%b3) + &
         TWO_PI*(real(cell%n_pt3-1,kind=DP)*cell%d3) / &
         geometry_MAGNITUDE(cell%a3)) / geometry_MAGNITUDE(cell%b3)

    ! cks: find the distance of the sphere from the a2_a3 planes defined
    !      by boundary points of the current ppd
    first_plane_a2_a3 = ((ppd_start - centre) .DOT. cell%b1) / &
         geometry_MAGNITUDE(cell%b1)

    second_plane_a2_a3 = (((ppd_start - centre) .DOT. cell%b1) + &
         TWO_PI*(real(cell%n_pt1-1,kind=DP)*cell%d1) / &
         geometry_MAGNITUDE(cell%a1)) / geometry_MAGNITUDE(cell%b1)

    ! cks: find the distance of the sphere from the a3_a1 planes defined
    !      by boundary points of the current ppd
    first_plane_a3_a1 = ((ppd_start - centre) .DOT. cell%b2) / &
         geometry_MAGNITUDE(cell%b2)

    second_plane_a3_a1 = (((ppd_start - centre) .DOT. cell%b2) + &
         TWO_PI*(real(cell%n_pt2-1,kind=DP)*cell%d2) / &
         geometry_MAGNITUDE(cell%a2)) / geometry_MAGNITUDE(cell%b2)

    ! cks: check if ppd contains the centre of sphere and therefore belongs
    !      to the sphere
    more_than_two = 0
    if (first_plane_a1_a2*second_plane_a1_a2 <= 0.0_DP) &
         more_than_two = more_than_two + 1
    if (first_plane_a2_a3*second_plane_a2_a3 <= 0.0_DP) &
         more_than_two = more_than_two + 1
    if (first_plane_a3_a1*second_plane_a3_a1 <= 0.0_DP) &
         more_than_two = more_than_two + 1

    if (more_than_two == 3) then
       basis_ppd_belongs_to_sphere = 1
       return
    end if


    ! =================== END SPHERE CENTRE INSIDE PPD? =======================


    !xxxxxxxxxxxx CHECK IF PPD BORDER POINTS BELONG TO SPHERE xxxxxxxxxxxxxxxxx

    ! cks: reaching this point means that the centre of the sphere is
    !      outside the ppd

    ! cks: Loop over all the boundary grid points of the current ppd and
    !      examine if any of them belongs to the current sphere.
    do loc_3=0, cell%n_pt3-1
       local_3=real(loc_3,kind=DP)*cell%d3
       do loc_2=0,cell%n_pt2-1
          local_2=real(loc_2,kind=DP)*cell%d2

          current_displacement = local_displacement(cell%a1_unit, &
               cell%a2_unit,cell%a3_unit,0.0_DP,local_2,local_3)

          current_point = ppd_start + current_displacement
          length = geometry_distance(current_point,centre)

          ! cks: test if the points at the plane at the beginning of the a1
          !      direction belong
          if (length <= NGWF_radius) then
             basis_ppd_belongs_to_sphere = 1
             return
          end if

          current_displacement = local_displacement(cell%a1_unit, &
               cell%a2_unit,cell%a3_unit, &
               real(cell%n_pt1-1,kind=DP)*cell%d1, &
               local_2,local_3)

          current_point = ppd_start + current_displacement
          length = geometry_distance(current_point,centre)

          ! cks: test if the points at the plane at the end of the a1 direction
          !      belong
          if (length <= NGWF_radius) then
             basis_ppd_belongs_to_sphere = 1
             return
          end if

          if (loc_3 == 0 .or. loc_3 == cell%n_pt3-1 .or. loc_2 == 0 .or. &
               loc_2 == cell%n_pt2-1) then

             do loc_1=1, cell%n_pt1-2
                local_1=real(loc_1,kind=DP)*cell%d1

                current_displacement = local_displacement(cell%a1_unit,&
                     cell%a2_unit, cell%a3_unit,local_1,local_2,local_3)

                current_point = ppd_start + current_displacement
                length = geometry_distance(current_point,centre)

                ! cks: test if the points at the planes at the end and at the
                !      beginning of the a2 and a3 directions belong
                if (length <= NGWF_radius) then
                   basis_ppd_belongs_to_sphere = 1
                   return
                end if

             end do

          end if
       end do
    end do

    !xxxxxxxxxxxx END CHECK IF PPD BORDER POINTS BELONG TO SPHERE xxxxxxxxxxxxx

  end function basis_ppd_belongs_to_sphere


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_sphere_deallocate(current_sphere)

    !===================================================!
    ! Deallocate pointers in the sphere object.         !
    !---------------------------------------------------!
    ! Written by Victor Milman on 03/11/2006.           !
    ! Modified by Nicholas Hine on 24/06/2009 to        !
    ! combine ppd_loc with ppd_list for MPI efficiency. !
    !===================================================!

    use utils, only : utils_dealloc_check
    implicit none
    type(SPHERE), intent(inout) :: current_sphere

    integer :: ierr

    ! deallocate "sphere" pointers
    if (allocated(current_sphere%ppd_list)) then
       deallocate(current_sphere%ppd_list,stat=ierr)
       call utils_dealloc_check('basis_sphere_deallocate', &
            'current_sphere%ppd_list',ierr)
    end if

  end subroutine basis_sphere_deallocate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  type(POINT) function basis_ppd_location(ppd,cell)

    !====================================================================!
    ! This function returns the cartesian coordinates of the origin of a !
    ! given ppd.                                                         !
    !--------------------------------------------------------------------!
    ! Written by Nicholas Hine on 19/10/2009 to replace the previous     !
    ! array of ppd locations as that was O(N) memory and not O(1/Nproc). !
    !====================================================================!

    use geometry, only: operator(*), operator(+)
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    integer, intent(in) :: ppd
    type(CELL_INFO), intent(in) :: cell

    ! Local variables
    integer :: d1,d2,d3
    real(kind=DP) :: rd1, rd2, rd3

    ! ndmh: find integer coordinates of this ppd along each lattice direction
    d1 = modulo(ppd-1,cell%n_ppds_a1)
    d2 = modulo((ppd-d1-1)/cell%n_ppds_a1, cell%n_ppds_a2)
    d3 = (ppd-d1-d2*cell%n_ppds_a1-1)/(cell%n_ppds_a2*cell%n_ppds_a1)

    ! ndmh: find distance along each lattice direction of this ppd
    rd1 = real(d1,kind=DP)*real(cell%n_pt1,kind=DP)*cell%d1
    rd2 = real(d2,kind=DP)*real(cell%n_pt2,kind=DP)*cell%d2
    rd3 = real(d3,kind=DP)*real(cell%n_pt3,kind=DP)*cell%d3

    ! ndmh: find cartesian coordinates of this ppd in terms of unit vectors
    ! ndmh: along each lattice vector (avoiding geometry_mod function calls)
    basis_ppd_location%x = rd1*cell%a1_unit%x + &
                           rd2*cell%a2_unit%x + &
                           rd3*cell%a3_unit%x
    basis_ppd_location%y = rd1*cell%a1_unit%y + &
                           rd2*cell%a2_unit%y + &
                           rd3*cell%a3_unit%y
    basis_ppd_location%z = rd1*cell%a1_unit%z + &
                           rd2*cell%a2_unit%z + &
                           rd3*cell%a3_unit%z

  end function basis_ppd_location


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_copy_function_to_box(fa_box, &                       ! in-out
       offset1, offset2, offset3, fa_tightbox, fa_on_grid, fa_sphere, & ! input
       cell, fftbox, opt_offset_override)                               ! input

    !=================================================================!
    ! This subroutine copies an NGWF stored in its ppd representation !
    ! to a three-dimensional array with the size of the function's    !
    ! tightbox.                                                       !
    !-----------------------------------------------------------------!
    ! @docme                                                          !
    ! opt_offset_override (in, opt): Optional override of the offset  !
    !                                used to set fa_point.            !
    !-----------------------------------------------------------------!
    ! Copied from basis_copy_function_to_fftbox by Nicholas Hine on   !
    ! 23/07/2009 to replace old double-loop version of same. See      !
    ! comments on above.                                              !
    ! Amalgamated all the copy_function_to_*box routines into one     !
    ! multi-purpose routine, Nicholas Hine, 28/02/2011.               !
    ! Modified to add optional offset override by Jacek Dziedzic.     !
    ! Modified by Andrea Greco on 07/04/2015 to make it compatible    !
    ! with new FUNCTIONS and FFTBOX_DATA types                        !
    !=================================================================!

    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(FUNCTION_TIGHT_BOX), intent(in) :: fa_tightbox ! function tightbox
    ! integer, intent(in) :: box_n1,box_n2,box_n3
    integer, intent(in) :: offset1,offset2,offset3
    ! real(kind=DP), intent(out) :: fa_box(box_n1,box_n2,box_n3) ! box data
    type(FFTBOX_DATA), intent(inout) :: fa_box ! box data
    ! real(kind=DP), intent(in) :: fa_on_grid(:) ! function data to copy from
    type(FUNCTIONS), intent(in) :: fa_on_grid ! function data to copy from
    type(SPHERE), intent(in) :: fa_sphere ! sphere of function's PPDs
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    integer, intent(in), optional :: opt_offset_override

    ! Local Variables
    integer :: fa_ppd      ! PPD counter within sphere
    integer :: fa_point    ! point counter within fa_on_grid
    integer :: a1_pos      ! absolute 1-coordinate of current PPD
    integer :: a2_pos      ! absolute 2-coordinate of current PPD
    integer :: a3_pos      ! absolute 3-coordinate of current PPD
    integer :: start1      ! start of points within box along direction 1 of current PPD
    integer :: start2      ! start of points within box along direction 2 of current PPD
    integer :: start3      ! start of points within box along direction 3 of current PPD
    integer :: finish1     ! end of points within box along direction 1 of current PPD
    integer :: finish2     ! end of points within box along direction 2 of current PPD
    integer :: finish3     ! end of points within box along direction 3 of current PPD
    integer :: box_offset1 ! offset to express 1-start in box
    integer :: box_offset2 ! offset to express 2-start in box
    integer :: box_offset3 ! offset to express 3-start in box
    integer :: box_pt1     ! current point in box along direction-1
    integer :: box_pt2     ! current point in box along direction-2
    integer :: box_pt3     ! current point in box along direction-3
    integer :: box_start1  ! start point of box points along direction-1
    integer :: box_end1    ! end point of box points along direction-1
    integer :: point2      ! point counter wrt current PPD along direction-2
    integer :: point3      ! point counter wrt current PPD along direction-3
    integer :: fa_offset

    call timer_clock('basis_copy_function_to_box',1)

    ! agreco: check we copy real to real or complex to complex
    call utils_assert(fa_on_grid%iscmplx.eqv.fa_box%iscmplx,'Error in &
       &basis_copy_function_to_box: incompatible argument types')

    call data_set_to_zero(fa_box)

    if(present(opt_offset_override)) then
       fa_offset = opt_offset_override
    else
       fa_offset = fa_sphere%offset
    end if

    if (cell%n_pts == 1) then
       ! cks: case of only one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               cell%n_ppds_a1, cell%n_ppds_a2, cell%n_ppds_a3, fftbox)

          box_pt1 = a1_pos - fa_tightbox%start_ppds1 + offset1
          box_pt2 = a2_pos - fa_tightbox%start_ppds2 + offset2
          box_pt3 = a3_pos - fa_tightbox%start_ppds3 + offset3

          fa_point = fa_offset + fa_ppd -1

          ! ndmh: copy ngwf value to box
          if (fa_box%iscmplx) then
             fa_box%z(box_pt1, box_pt2, box_pt3) = fa_on_grid%z(fa_point)
          else
             fa_box%d(box_pt1, box_pt2, box_pt3) = fa_on_grid%d(fa_point)
          endif

       enddo

    else

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               cell%n_ppds_a1, cell%n_ppds_a2, cell%n_ppds_a3, fftbox)

          call basis_lims_1d_in_ppd_in_tight(start1, finish1, box_offset1, &
               a1_pos, fa_tightbox%start_ppds1, fa_tightbox%finish_ppds1, &
               fa_tightbox%start_pts1, fa_tightbox%finish_pts1, cell%n_pt1)

          call basis_lims_1d_in_ppd_in_tight(start2, finish2, box_offset2, &
               a2_pos, fa_tightbox%start_ppds2, fa_tightbox%finish_ppds2, &
               fa_tightbox%start_pts2, fa_tightbox%finish_pts2, cell%n_pt2)

          call basis_lims_1d_in_ppd_in_tight(start3, finish3, box_offset3, &
               a3_pos, fa_tightbox%start_ppds3, fa_tightbox%finish_ppds3, &
               fa_tightbox%start_pts3, fa_tightbox%finish_pts3, cell%n_pt3)

          ! ndmh: find the first point in the current ppd that needs to be copied
          fa_point = fa_offset + (fa_ppd-1)*cell%n_pts &
               + (start3-1)*cell%n_pt2*cell%n_pt1 &
               + (start2-1)*cell%n_pt1 + (start1-1)

          ! ndmh: find range of points in '1' direction
          box_start1 = box_offset1 + offset1
          box_end1 = box_start1 + (finish1 - start1)

          do point3=start3,finish3
             ! ndmh: coordinates of points with respect to the box,
             !      i.e. the first point inside the box in a certain
             !      dimension starts from 1, etc.
             box_pt3 = point3 - start3 + box_offset3 + offset3
             do point2=start2,finish2

                box_pt2 = point2 - start2 + box_offset2 + offset2

                ! ndmh: copy '1' direction line of function values to box
                if (fa_box%iscmplx) then
                   fa_box%z(box_start1:box_end1,box_pt2,box_pt3) = &
                      fa_on_grid%z(fa_point:fa_point+(finish1-start1))
                else
                   fa_box%d(box_start1:box_end1,box_pt2,box_pt3) = &
                      fa_on_grid%d(fa_point:fa_point+(finish1-start1))
                endif

                ! ndmh: skip to next line in '2' direction
                fa_point = fa_point + cell%n_pt1

             end do

             ! ndmh: skip to next line in '3' direction
             fa_point = fa_point + (cell%n_pt2-finish2+start2-1)* &
                  cell%n_pt1
          end do

       end do

    end if

    call timer_clock('basis_copy_function_to_box',2)

  end subroutine basis_copy_function_to_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_add_function_to_box(fa_box, &                  ! in-out
       fa_start1, fa_start2, fa_start3, fa_tightbox, &            ! input
       fa_on_grid, fa_sphere, factor, cell, fftbox)             ! input

    !==================================================================!
    ! This subroutine adds an NGWF multiplied by some factor to a      !
    ! specified position with a three-dimensional array with the size  !
    ! of the fftbox, reading the ngwf directly from its ppd            !
    ! representation.                                                  !
    !------------------------------------------------------------------!
    ! Adapted from code originally written by Chris-Kriton Skylaris    !
    ! on 11/6/2001 for the ONES code, for extracting ppds from a box,  !
    ! called "basis_extract_ppds_from_pairbox"                         !
    ! Renamed by Arash A. Mostofi in August 2003 and removed the       !
    ! "make_pair_even" increment from the FFTbox dimensions.           !
    ! Modified for speed by Chris-Kriton Skylaris on 08/09/2005.       !
    ! Modified to work faster when the ppds contain only one point     !
    ! by Chris-Kriton Skylaris on 31/12/2006.                          !
    ! Modified to work in reverse to deposit ppds to the fftbox, and   !
    ! loop-unrolled for speed by Nicholas Hine on 09/05/2008.          !
    ! Modified by Nicholas Hine on 24/06/2009 to combine ppd_loc       !
    ! with ppd_list for MPI efficiency.                                !
    ! Renamed by Nicholas Hine on 28/02/2011 and adapted to work       !
    ! box of any size rather than just fftbox.                         !
    ! Modified by Andrea Greco on 13/04/2015 to make it compatible     !
    ! with new FUNCTIONS and FFTBOX_DATA types                         !
    !==================================================================!

    use datatypes, only: COEF, FUNCTIONS, FFTBOX_DATA
    use fft_box, only: FFTBOX_INFO
!$  use rundat, only: pub_threads_per_fftbox
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    ! integer, intent(in) :: box_n1,box_n2,box_n3
    ! real(kind=DP), intent(inout) :: fa_box(box_n1,box_n2,box_n3) ! box data
    type(FFTBOX_DATA), intent(inout) :: fa_box ! box data
    integer, intent(in) :: fa_start1 ! 1-beginning of NGWF tightbox wrt box
    integer, intent(in) :: fa_start2 ! 2-beginning of NGWF tightbox wrt box
    integer, intent(in) :: fa_start3 ! 3-beginning of NGWF tightbox wrt box
    type(FUNCTION_TIGHT_BOX), intent(in) :: fa_tightbox ! NGWF tightbox
    ! real(kind=DP), intent(in), dimension(:) :: fa_on_grid ! NGWF to add
    type(FUNCTIONS), intent(in) :: fa_on_grid ! NGWF to add
    type(SPHERE), intent(in) :: fa_sphere ! sphere of NGWFs PPDs
    ! real(kind=DP), intent(in) :: factor ! ndmh: factor to multiply NGWF by
    type(COEF), intent(in) :: factor ! ndmh: factor to multiply NGWF by
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local Variables
    integer :: fa_ppd      ! PPD counter within sphere
    integer :: fa_point    ! point counter within fa_on_grid
    integer :: a1_pos      ! absolute 1-coordinate of current PPD
    integer :: a2_pos      ! absolute 2-coordinate of current PPD
    integer :: a3_pos      ! absolute 3-coordinate of current PPD
    integer :: start1      ! start of points within tightbox along direction 1 of current PPD
    integer :: start2      ! start of points within tightbox along direction 2 of current PPD
    integer :: start3      ! start of points within tightbox along direction 3 of current PPD
    integer :: finish1     ! end of points within tightbox along direction 1 of current PPD
    integer :: finish2     ! end of points within tightbox along direction 2 of current PPD
    integer :: finish3     ! end of points within tightbox along direction 3 of current PPD
    integer :: box_offset1 ! offset to express 1-start of function in box
    integer :: box_offset2 ! offset to express 2-start of function in box
    integer :: box_offset3 ! offset to express 3-start of function in box
    integer :: box_pt1     ! point counter wrt box along direction-1
    integer :: box_pt2     ! point counter wrt box along direction-2
    integer :: box_pt3     ! point counter wrt box along direction-3
    integer :: box_start1  ! start point of box points along direction-1
    integer :: box_end1    ! end point of box points along direction-1
    integer :: point2      ! point counter wrt current PPD along direction-2
    integer :: point3      ! point counter wrt current PPD along direction-3

    ! -------------------------------------------------------------------------

    call timer_clock('basis_add_function_to_box',1)

    ! agreco: check we copy real to real or complex to complex
    call utils_assert(fa_on_grid%iscmplx.eqv.fa_box%iscmplx,'Error in &
       &basis_add_function_to_box: incompatible argument types')

    if (cell%n_pts == 1) then
       ! cks: case of only one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               cell%n_ppds_a1, cell%n_ppds_a2, cell%n_ppds_a3, fftbox)

          box_pt1 = a1_pos - fa_tightbox%start_ppds1 + fa_start1
          box_pt2 = a2_pos - fa_tightbox%start_ppds2 + fa_start2
          box_pt3 = a3_pos - fa_tightbox%start_ppds3 + fa_start3

          fa_point = fa_sphere%offset + fa_ppd - 1

          ! ndmh: add factor*(ngwf value) to box
          if (fa_box%iscmplx) then
             if (factor%iscmplx) then
                 fa_box%z(box_pt1, box_pt2, box_pt3) = &
                     fa_box%z(box_pt1, box_pt2, box_pt3) + &
                     fa_on_grid%z(fa_point)*factor%z
             else
                fa_box%z(box_pt1, box_pt2, box_pt3) = &
                     fa_box%z(box_pt1, box_pt2, box_pt3) + &
                     fa_on_grid%z(fa_point)*factor%d
             end if
          else
             if (factor%iscmplx) then
                fa_box%d(box_pt1, box_pt2, box_pt3) = &
                     fa_box%d(box_pt1, box_pt2, box_pt3) + &
                     fa_on_grid%d(fa_point)*real(factor%z)
             else
                fa_box%d(box_pt1, box_pt2, box_pt3) = &
                     fa_box%d(box_pt1, box_pt2, box_pt3) + &
                     fa_on_grid%d(fa_point)*factor%d
             end if
          end if

       enddo

    else
       ! cks: case for more than one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
!$OMP PARALLEL DO NUM_THREADS(pub_threads_per_fftbox) DEFAULT(none) &
!$OMP PRIVATE(fa_ppd, a1_pos, a2_pos, a3_pos, start1, finish1, &
!$OMP      box_offset1, start2, finish2, box_offset2, start3, finish3, &
!$OMP      box_offset3, fa_point, box_start1, box_end1, point3, box_pt3, &
!$OMP      point2, box_pt2) &
!$OMP SHARED(fa_sphere, cell, fa_box, factor, fa_tightbox, fa_start1, &
!$OMP      fa_start2, fa_start3, fa_on_grid, fftbox, pub_threads_per_fftbox)
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               cell%n_ppds_a1, cell%n_ppds_a2, cell%n_ppds_a3, fftbox)

          call basis_lims_1d_in_ppd_in_tight(start1, finish1, box_offset1, &
               a1_pos, fa_tightbox%start_ppds1, fa_tightbox%finish_ppds1, &
               fa_tightbox%start_pts1, fa_tightbox%finish_pts1, cell%n_pt1)

          call basis_lims_1d_in_ppd_in_tight(start2, finish2, box_offset2, &
               a2_pos, fa_tightbox%start_ppds2, fa_tightbox%finish_ppds2, &
               fa_tightbox%start_pts2, fa_tightbox%finish_pts2, cell%n_pt2)

          call basis_lims_1d_in_ppd_in_tight(start3, finish3, box_offset3, &
               a3_pos, fa_tightbox%start_ppds3, fa_tightbox%finish_ppds3, &
               fa_tightbox%start_pts3, fa_tightbox%finish_pts3, cell%n_pt3)

          ! cks: loop over the points of the current ppd that also belong to
          !      the box

          fa_point = fa_sphere%offset + (fa_ppd-1)*cell%n_pts - 1 &
               + (start3-1)*cell%n_pt2*cell%n_pt1 &
               + (start2-1)*cell%n_pt1 + start1

          box_start1 = box_offset1 + fa_start1
          box_end1 = box_start1 + finish1 - start1

          do point3=start3,finish3
             ! cks: coordinates of points with respect to the box,
             !      i.e. the first point inside the box in a certain
             !      dimension starts from 1, etc.
             box_pt3 = point3 - start3 + box_offset3 + fa_start3
             do point2=start2,finish2

                box_pt2 = point2 - start2 + box_offset2 + fa_start2

                ! ndmh: removed inner loop for speed
                ! ndmh: add factor*(ngwf value) to box
                if (fa_box%iscmplx) then
                   if(factor%iscmplx) then
                     fa_box%z(box_start1:box_end1,box_pt2,box_pt3) = &
                     fa_box%z(box_start1:box_end1,box_pt2,box_pt3) + &
                     fa_on_grid%z(fa_point:fa_point+(finish1-start1)) * &
                     factor%z
                   else
                     fa_box%z(box_start1:box_end1,box_pt2,box_pt3) = &
                     fa_box%z(box_start1:box_end1,box_pt2,box_pt3) + &
                     fa_on_grid%z(fa_point:fa_point+(finish1-start1)) * &
                     factor%d
                   end if
                else
                   if (factor%iscmplx) then
                     fa_box%d(box_start1:box_end1,box_pt2,box_pt3) = &
                     fa_box%d(box_start1:box_end1,box_pt2,box_pt3) + &
                     fa_on_grid%d(fa_point:fa_point+(finish1-start1)) * &
                     real(factor%z)
                   else
                     fa_box%d(box_start1:box_end1,box_pt2,box_pt3) = &
                     fa_box%d(box_start1:box_end1,box_pt2,box_pt3) + &
                     fa_on_grid%d(fa_point:fa_point+(finish1-start1))*factor%d
                   end if
                end if

                fa_point = fa_point + cell%n_pt1

             end do

             fa_point = fa_point + &
                  (cell%n_pt2-finish2+start2-1)*cell%n_pt1

          end do

       end do
!$OMP END PARALLEL DO

    endif

    call timer_clock('basis_add_function_to_box',2)

  end subroutine basis_add_function_to_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_extract_function_from_box(fa_on_grid, &          ! in-out
       fa_box, fa_sphere, fa_tightbox, fa_start1,  &                ! input
       fa_start2, fa_start3, ngwf_offset, cell, fftbox)             ! input

    !================================================================!
    ! This subroutines extracts PPDs from an FFTbox.                 !
    ! The PPDs which belong to the NGWF "fa" are extracted           !
    ! and stored in the order dictated by the sphere of "fa".        !
    !----------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris on 11/6/2001 for   !
    ! the ONES code and called "basis_extract_ppds_from_pairbox"     !
    ! Renamed by Arash A. Mostofi in August 2003 and removed the     !
    ! "make_pair_even" increment from the FFTbox dimensions.         !
    ! Modified for speed by Chris-Kriton Skylaris on 08/09/2005.     !
    ! Modified to work faster when the ppds contain only one point   !
    ! by Chris-Kriton Skylaris on 31/12/2006.                        !
    ! Modified by Nicholas Hine on 24/06/2009 to combine ppd_loc     !
    ! with ppd_list for MPI efficiency.                              !
    ! Renamed by Nicholas Hine on 28/02/2011 and adapted to work     !
    ! box of any size rather than just fftbox.                       !
    ! Modified by Andrea Greco on 07/04/2015 to ensure compatibility !
    ! with new FUNCTIONS and FFTBOX_DATA types                       !
    !================================================================!

    use constants, only: cmplx_0
    use datatypes, only: FUNCTIONS, FFTBOX_DATA
    use fft_box, only: FFTBOX_INFO
!$  use rundat, only: pub_threads_per_fftbox
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    ! integer, intent(in) :: box_n1,box_n2,box_n3
    type(SPHERE), intent(in) :: fa_sphere ! sphere of NGWFs PPDs we want
    type(FUNCTION_TIGHT_BOX), intent(in) :: fa_tightbox ! tightbox "
    integer, intent(in) :: fa_start1 ! 1-beginning of extract region wrt box
    integer, intent(in) :: fa_start2 ! 2-beginning of extract region wrt box
    integer, intent(in) :: fa_start3 ! 3-beginning of extract region wrt box
    integer, intent(in) :: ngwf_offset ! beginning of NGWF in PPD rep
    ! real(kind=DP), intent(in) :: fa_box(box_n1,box_n2,box_n3) ! extraction box
    type(FFTBOX_DATA), intent(in) :: fa_box ! extraction box
    ! real(kind=DP), intent(inout), dimension(:) :: fa_on_grid ! extracted NGWF
    type(FUNCTIONS), intent(inout) :: fa_on_grid ! extracted NGWF
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local Variables
    integer :: fa_ppd     ! PPD counter within sphere
    integer :: fa_point   ! point counter within fa_on_grid
    integer :: a1_pos     ! absolute 1-coordinate of current PPD
    integer :: a2_pos     ! absolute 2-coordinate of current PPD
    integer :: a3_pos     ! absolute 3-coordinate of current PPD
    integer :: start1     ! start of points within tightbox along direction 1 of current PPD
    integer :: start2     ! start of points within tightbox along direction 2 of current PPD
    integer :: start3     ! start of points within tightbox along direction 3 of current PPD
    integer :: finish1    ! end of points within tightbox along direction 1 of current PPD
    integer :: finish2    ! end of points within tightbox along direction 2 of current PPD
    integer :: finish3    ! end of points within tightbox along direction 3 of current PPD
    integer :: box_offset1 ! offset to express 1-start of function in box
    integer :: box_offset2 ! offset to express 2-start of function in box
    integer :: box_offset3 ! offset to express 3-start of function in box
    integer :: box_pt1    ! point counter wrt box along direction-1
    integer :: box_pt2    ! point counter wrt box along direction-2
    integer :: box_pt3    ! point counter wrt box along direction-3
    integer :: box_start1 ! start point in box along direction-1
    integer :: box_end1   ! end point in box along direction-1
    integer :: point2     ! point counter wrt current PPD along direction-2
    integer :: point3     ! point counter wrt current PPD along direction-3

    ! -------------------------------------------------------------------------

    call timer_clock('basis_extract_function_from_box', 1)

    ! agreco: check we copy real to real or complex to complex
    call utils_assert(fa_on_grid%iscmplx.eqv.fa_box%iscmplx,'Error in &
       &basis_extract_function_from_box: incompatible argument types')

    ! cks: initialise
    if (fa_on_grid%iscmplx) then
       fa_on_grid%z(ngwf_offset: &
            ngwf_offset - 1 + fa_sphere%n_ppds_sphere*cell%n_pts) = cmplx_0
    else
       fa_on_grid%d(ngwf_offset: &
            ngwf_offset - 1 + fa_sphere%n_ppds_sphere*cell%n_pts) = 0.0_DP
    endif

    if (cell%n_pts == 1) then
       ! cks: case of only one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               cell%n_ppds_a1, cell%n_ppds_a2, cell%n_ppds_a3, fftbox)

          box_pt1 = a1_pos - fa_tightbox%start_ppds1 + fa_start1
          box_pt2 = a2_pos - fa_tightbox%start_ppds2 + fa_start2
          box_pt3 = a3_pos - fa_tightbox%start_ppds3 + fa_start3

          fa_point = ngwf_offset + fa_ppd -1

          if (fa_on_grid%iscmplx) then
             fa_on_grid%z(fa_point) = fa_box%z(box_pt1, box_pt2, box_pt3)
          else
             fa_on_grid%d(fa_point) = fa_box%d(box_pt1, box_pt2, box_pt3)
          endif

       enddo

    else
       ! cks: case for more than one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
!$OMP PARALLEL DO NUM_THREADS(pub_threads_per_fftbox) DEFAULT(none) &
!$OMP PRIVATE(fa_ppd, a1_pos, a2_pos, a3_pos, start1, finish1, &
!$OMP      box_offset1, start2, finish2, box_offset2, start3, finish3, &
!$OMP      box_offset3, fa_point, box_start1, box_end1, point3, box_pt3, &
!$OMP      point2, box_pt2) &
!$OMP SHARED(fa_sphere, cell, ngwf_offset, fa_box, fa_tightbox, fa_start1, &
!$OMP      fa_start2, fa_start3, fa_on_grid, fftbox, pub_threads_per_fftbox)
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               cell%n_ppds_a1, cell%n_ppds_a2, cell%n_ppds_a3, fftbox)

          call basis_lims_1d_in_ppd_in_tight(start1, finish1, box_offset1, &
               a1_pos, fa_tightbox%start_ppds1, fa_tightbox%finish_ppds1, &
               fa_tightbox%start_pts1, fa_tightbox%finish_pts1, cell%n_pt1)

          call basis_lims_1d_in_ppd_in_tight(start2, finish2, box_offset2, &
               a2_pos, fa_tightbox%start_ppds2, fa_tightbox%finish_ppds2, &
               fa_tightbox%start_pts2, fa_tightbox%finish_pts2, cell%n_pt2)

          call basis_lims_1d_in_ppd_in_tight(start3, finish3, box_offset3, &
               a3_pos, fa_tightbox%start_ppds3, fa_tightbox%finish_ppds3, &
               fa_tightbox%start_pts3, fa_tightbox%finish_pts3, cell%n_pt3)

          fa_point = ngwf_offset + (fa_ppd-1)*cell%n_pts - 1 &
               + (start3-1)*cell%n_pt2*cell%n_pt1 &
               + (start2-1)*cell%n_pt1 + start1

          box_start1 = box_offset1 + fa_start1
          box_end1 = box_start1 + finish1 - start1

          do point3=start3,finish3
             ! cks: coordinates of points with respect to the FFT box,
             !      i.e. the first point inside the FFT box in a certain
             !      dimension starts from 1, etc.
             box_pt3 = point3 - start3 + box_offset3 + fa_start3
             do point2=start2,finish2

                box_pt2 = point2 - start2 + box_offset2 + fa_start2

                if (fa_on_grid%iscmplx) then
                   fa_on_grid%z(fa_point:fa_point+(finish1-start1)) = &
                        fa_box%z(box_start1:box_end1,box_pt2, box_pt3)
                else
                   fa_on_grid%d(fa_point:fa_point+(finish1-start1)) = &
                        fa_box%d(box_start1:box_end1,box_pt2, box_pt3)
                endif

                fa_point = fa_point + cell%n_pt1

             end do

             fa_point = fa_point + &
                  (cell%n_pt2-finish2+start2-1)*cell%n_pt1

          end do

       end do
!$OMP END PARALLEL DO

    endif

    call timer_clock('basis_extract_function_from_box', 2)

  end subroutine basis_extract_function_from_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_dot_function_with_box(dot_prod, fa_on_grid, &    ! out/input
       fa_box, fa_sphere, fa_tightbox, fa_start1,  &                ! input
       fa_start2, fa_start3, ngwf_offset, cell, fftbox, &           ! input
       factor_type, allow_mix_types)                                ! input

    !================================================================!
    ! This subroutines takes the dot product of a function stored    !
    ! on PPDs (eg an NGWF "fa") with the corresponding regions of    !
    ! an FFTbox (eg the product V_eff|fb> filtered to the coarse     !
    ! grid.                                                          !
    ! The PPDs which belong to the NGWF "fa" and the box are handled !
    ! in the order dictated by the sphere of "fa".                   !
    !----------------------------------------------------------------!
    ! This function was created on 27/06/2013 as an offshoot of      !
    ! basis_extract_function_from_box tailored to a specific purpose !
    ! by Nicholas Hine.                                              !
    !----------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris on 11/6/2001 for   !
    ! the ONES code and called "basis_extract_ppds_from_pairbox"     !
    ! Renamed by Arash A. Mostofi in August 2003 and removed the     !
    ! "make_pair_even" increment from the FFTbox dimensions.         !
    ! Modified for speed by Chris-Kriton Skylaris on 08/09/2005.     !
    ! Modified to work faster when the ppds contain only one point   !
    ! by Chris-Kriton Skylaris on 31/12/2006.                        !
    ! Modified by Nicholas Hine on 24/06/2009 to combine ppd_loc     !
    ! with ppd_list for MPI efficiency.                              !
    ! Renamed by Nicholas Hine on 28/02/2011 and adapted to work     !
    ! box of any size rather than just fftbox.                       !
    ! Modified by Andrea Greco on 13/04/2015 to make it compatible   !
    ! with new FUNCTIONS and FFTBOX_DATA types                       !
    !================================================================!

    use datatypes, only: COEF, FUNCTIONS, FFTBOX_DATA
    use constants, only: cmplx_0
    use fft_box, only: FFTBOX_INFO
!$  use rundat, only: pub_threads_per_fftbox
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! Arguments
    ! integer, intent(in) :: box_n1,box_n2,box_n3
    type(SPHERE), intent(in) :: fa_sphere ! sphere of NGWFs PPDs we want
    type(FUNCTION_TIGHT_BOX), intent(in) :: fa_tightbox ! tightbox "
    integer, intent(in) :: fa_start1 ! 1-beginning of extract region wrt box
    integer, intent(in) :: fa_start2 ! 2-beginning of extract region wrt box
    integer, intent(in) :: fa_start3 ! 3-beginning of extract region wrt box
    integer, intent(in) :: ngwf_offset ! beginning of NGWF in PPD rep
    ! real(kind=DP), intent(in) :: fa_box(box_n1,box_n2,box_n3) ! extraction box
    type(FFTBOX_DATA), intent(in) :: fa_box ! extraction box
    ! real(kind=DP), intent(in), dimension(:) :: fa_on_grid ! extracted NGWF
    type(FUNCTIONS), intent(in) :: fa_on_grid ! extracted NGWF
    ! real(kind=DP), intent(out) :: dot_prod ! dot product of NGWF and box
    type(COEF), intent(inout) :: dot_prod ! dot product of NGWF and box
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    integer, intent(in), optional :: factor_type
    logical, intent(in), optional :: allow_mix_types

    ! Local Variables
    integer :: fa_ppd     ! PPD counter within sphere
    integer :: fa_point   ! point counter within fa_on_grid
    integer :: a1_pos     ! absolute 1-coordinate of current PPD
    integer :: a2_pos     ! absolute 2-coordinate of current PPD
    integer :: a3_pos     ! absolute 3-coordinate of current PPD
    integer :: start1     ! start of points within tightbox along direction 1 of current PPD
    integer :: start2     ! start of points within tightbox along direction 2 of current PPD
    integer :: start3     ! start of points within tightbox along direction 3 of current PPD
    integer :: finish1    ! end of points within tightbox along direction 1 of current PPD
    integer :: finish2    ! end of points within tightbox along direction 2 of current PPD
    integer :: finish3    ! end of points within tightbox along direction 3 of current PPD
    integer :: box_offset1 ! offset to express 1-start of function in box
    integer :: box_offset2 ! offset to express 2-start of function in box
    integer :: box_offset3 ! offset to express 3-start of function in box
    integer :: box_pt1    ! point counter wrt box along direction-1
    integer :: box_pt2    ! point counter wrt box along direction-2
    integer :: box_pt3    ! point counter wrt box along direction-3
    integer :: box_start1 ! start point in box along direction-1
    integer :: box_end1   ! end point in box along direction-1
    integer :: point2     ! point counter wrt current PPD along direction-2
    integer :: point3     ! point counter wrt current PPD along direction-3
    character(len=*), parameter :: myself = 'basis_dot_function_with_box'
    logical :: loc_allow_mix_types
    real(kind=DP) :: dot_prod_real        ! auxiliary variables for reduction
    complex(kind=DP) :: dot_prod_complex
    real(kind=DP) :: a1_neighbour, a2_neighbour, a3_neighbour, factor
!! Not supported by intel fortran 15 or older.
!!!$OMP DECLARE REDUCTION (+:COEF:omp_out=data_coef_sum(omp_out,omp_in)) &
!!!$OMP      initializer(omp_priv=omp_orig)

    ! -------------------------------------------------------------------------

    call timer_clock('basis_dot_function_with_box', 1)

    loc_allow_mix_types = .false.
    if (present(allow_mix_types)) loc_allow_mix_types = allow_mix_types

    ! jme: check types
    if ( (.not. dot_prod%iscmplx) .and. &
         (fa_box%iscmplx .or. fa_on_grid%iscmplx) ) then
       ! real result but one of the inputs is complex
       call utils_abort('Error in '//myself//': incompatible argument types: &
            &iscmplx components of dot_prod (result), fa_on_grid, fa_box are: ', &
            opt_logical_to_print1=dot_prod%iscmplx, &
            opt_logical_to_print2=fa_on_grid%iscmplx, &
            opt_logical_to_print3=fa_box%iscmplx)
    else if ( ( (dot_prod%iscmplx .neqv. fa_box%iscmplx) .or. &
                (dot_prod%iscmplx .neqv. fa_on_grid%iscmplx) ) .and. &
              (.not. loc_allow_mix_types) ) then
       ! types are mixed but this has not been enabled
       call utils_abort('Error in '//myself//': incompatible argument types: &
            &mixing of types not enabled, but iscmplx components of dot_prod &
            &(result), fa_on_grid, fa_box are: ', &
            opt_logical_to_print1=dot_prod%iscmplx, &
            opt_logical_to_print2=fa_on_grid%iscmplx, &
            opt_logical_to_print3=fa_box%iscmplx)
    end if

    ! ndmh: initialise
!!    call data_set_to_zero(dot_prod)
    dot_prod_real = 0.0_DP
    dot_prod_complex = cmplx_0

    if (cell%n_pts == 1) then
       ! cks: case of only one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          if (present(factor_type)) then
             a1_neighbour = real( nint(real(fa_sphere%ppd_list(2,fa_ppd), &
                  kind=DP) / 9.0_DP), kind=DP)
             a2_neighbour = real( nint(real((fa_sphere%ppd_list(2,fa_ppd) - &
                  9 * int(a1_neighbour)), kind=DP) / 3.0_DP), kind=DP)
             a3_neighbour = real( (fa_sphere%ppd_list(2,fa_ppd) - &
                  9 * int(a1_neighbour) - 3 * int(a2_neighbour)), kind=DP)

             if (factor_type == 1) then
                factor = a1_neighbour * cell%a1%X + &
                     a2_neighbour * cell%a2%X + a3_neighbour * cell%a3%X
             else if (factor_type == 2) then
                factor = a1_neighbour * cell%a1%Y + &
                     a2_neighbour * cell%a2%Y + a3_neighbour * cell%a3%Y
             else if (factor_type == 3) then
                factor = a1_neighbour * cell%a1%Z + &
                     a2_neighbour * cell%a2%Z + a3_neighbour * cell%a3%Z
             end if

             ! gcc32: there is a -1 factor still needed because we calculate
             ! \sum_R R \int phi(r-R) p(r) dr = \sum_R R \int phi(r) p(r+R) dr
             ! =  \sum_R (-R) \int phi(r) p(r-R) dr
             ! and the terms before consist of \sum_R R \int phi(r) p (r-R) dr
             factor = -1.0_DP * factor

          else
             factor = 1.0_DP
          end if

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               cell%n_ppds_a1, cell%n_ppds_a2, cell%n_ppds_a3, fftbox)

          box_pt1 = a1_pos - fa_tightbox%start_ppds1 + fa_start1
          box_pt2 = a2_pos - fa_tightbox%start_ppds2 + fa_start2
          box_pt3 = a3_pos - fa_tightbox%start_ppds3 + fa_start3

          fa_point = ngwf_offset + fa_ppd -1

          if (dot_prod%iscmplx) then
             if (fa_on_grid%iscmplx) then
                if(fa_box%iscmplx) then
                   dot_prod_complex = dot_prod_complex + &
                        cmplx(factor,0.0_DP,kind=DP) * &
                        conjg(fa_on_grid%z(fa_point)) * &
                        fa_box%z(box_pt1, box_pt2, box_pt3)
                else
                   dot_prod_complex = dot_prod_complex + &
                        conjg(fa_on_grid%z(fa_point)) * &
                        cmplx(factor * fa_box%d(box_pt1, box_pt2, box_pt3), &
                              0.0_DP, kind=DP)
                end if
             else ! fa_box is real
                if(fa_box%iscmplx) then
                   dot_prod_complex = dot_prod_complex +  &
                        cmplx(factor * fa_on_grid%d(fa_point),0.0_DP,kind=DP) &
                        * fa_box%z(box_pt1, box_pt2, box_pt3)
                else
                   dot_prod_complex = dot_prod_complex + cmplx(factor * &
                        fa_on_grid%d(fa_point) * &
                        fa_box%d(box_pt1, box_pt2, box_pt3), 0.0_DP, kind=DP)
                end if
             end if
          else
             dot_prod_real = dot_prod_real + factor * &
                fa_on_grid%d(fa_point) * fa_box%d(box_pt1, box_pt2, box_pt3)
          endif

       enddo

    else
       ! cks: case for more than one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa


!$OMP PARALLEL DO NUM_THREADS(pub_threads_per_fftbox) DEFAULT(none) &
!$OMP PRIVATE(fa_ppd, a1_pos, a2_pos, a3_pos, start1, finish1, &
!$OMP      box_offset1, start2, finish2, box_offset2, start3, finish3, &
!$OMP      box_offset3, fa_point, box_start1, box_end1, point3, box_pt3, &
!$OMP      point2, box_pt2, a1_neighbour, a2_neighbour, a3_neighbour, factor) &
!$OMP SHARED(fa_sphere, cell, ngwf_offset, fa_box, fa_tightbox, fa_start1, &
!$OMP      fa_start2, fa_start3, fa_on_grid, fftbox, pub_threads_per_fftbox, &
!$OMP      dot_prod, factor_type) &
!$OMP REDUCTION(+:dot_prod_real, dot_prod_complex)
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          a1_neighbour = real( nint(real(fa_sphere%ppd_list(2,fa_ppd), &
               kind=DP) / 9.0_DP), kind=DP)
          a2_neighbour = real( nint(real((fa_sphere%ppd_list(2,fa_ppd) - &
               9 * int(a1_neighbour)), kind=DP) / 3.0_DP), kind=DP)
          a3_neighbour = real( (fa_sphere%ppd_list(2,fa_ppd) - &
               9 * int(a1_neighbour) - 3 * int(a2_neighbour)), kind=DP)

          if (present(factor_type)) then
             if (factor_type == 1) then
                factor = a1_neighbour * cell%a1%X + &
                     a2_neighbour * cell%a2%X + a3_neighbour * cell%a3%X
             else if (factor_type == 2) then
                factor = a1_neighbour * cell%a1%Y + &
                     a2_neighbour * cell%a2%Y + a3_neighbour * cell%a3%Y
             else if (factor_type == 3) then
                factor = a1_neighbour * cell%a1%Z + &
                     a2_neighbour * cell%a2%Z + a3_neighbour * cell%a3%Z
             end if
          else
             factor = 1.0_DP
          end if

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               cell%n_ppds_a1, cell%n_ppds_a2, cell%n_ppds_a3, fftbox)

          call basis_lims_1d_in_ppd_in_tight(start1, finish1, box_offset1, &
               a1_pos, fa_tightbox%start_ppds1, fa_tightbox%finish_ppds1, &
               fa_tightbox%start_pts1, fa_tightbox%finish_pts1, cell%n_pt1)

          call basis_lims_1d_in_ppd_in_tight(start2, finish2, box_offset2, &
               a2_pos, fa_tightbox%start_ppds2, fa_tightbox%finish_ppds2, &
               fa_tightbox%start_pts2, fa_tightbox%finish_pts2, cell%n_pt2)

          call basis_lims_1d_in_ppd_in_tight(start3, finish3, box_offset3, &
               a3_pos, fa_tightbox%start_ppds3, fa_tightbox%finish_ppds3, &
               fa_tightbox%start_pts3, fa_tightbox%finish_pts3, cell%n_pt3)

          fa_point = ngwf_offset + (fa_ppd-1)*cell%n_pts - 1 &
               + (start3-1)*cell%n_pt2*cell%n_pt1 &
               + (start2-1)*cell%n_pt1 + start1

          box_start1 = box_offset1 + fa_start1
          box_end1 = box_start1 + finish1 - start1

          do point3=start3,finish3
             ! cks: coordinates of points with respect to the FFT box,
             !      i.e. the first point inside the FFT box in a certain
             !      dimension starts from 1, etc.
             box_pt3 = point3 - start3 + box_offset3 + fa_start3
             do point2=start2,finish2

               box_pt2 = point2 - start2 + box_offset2 + fa_start2

               if (dot_prod%iscmplx) then
                 if (fa_on_grid%iscmplx) then
                   if(fa_box%iscmplx) then
                     dot_prod_complex = dot_prod_complex + &
                       cmplx(factor,0.0_DP,kind=DP) * sum( &
                       conjg(fa_on_grid%z(fa_point:fa_point+(finish1-start1))) &
                       * fa_box%z(box_start1:box_end1, box_pt2, box_pt3) )
                   else
                     dot_prod_complex = dot_prod_complex + &
                       cmplx(factor,0.0_DP,kind=DP) * sum( &
                       conjg(fa_on_grid%z(fa_point:fa_point+(finish1-start1))) &
                       * fa_box%d(box_start1:box_end1, box_pt2, box_pt3) )
                   end if
                 else ! fa_box is real
                   if(fa_box%iscmplx) then
                     dot_prod_complex = dot_prod_complex + &
                       cmplx(factor,0.0_DP,kind=DP) * sum( &
                       fa_on_grid%d(fa_point:fa_point+(finish1-start1)) &
                       * fa_box%z(box_start1:box_end1, box_pt2, box_pt3) )
                   else
                     dot_prod_complex = dot_prod_complex + &
                       cmplx(factor * sum( &
                       fa_on_grid%d(fa_point:fa_point+(finish1-start1)) &
                       * fa_box%d(box_start1:box_end1, box_pt2, box_pt3) ), &
                       0.0_DP, kind=DP)
                   end if
                 end if

               else
                  dot_prod_real = dot_prod_real + factor * &
                       sum(fa_on_grid%d(fa_point:fa_point+(finish1-start1))* &
                       fa_box%d(box_start1:box_end1, box_pt2, box_pt3))
               endif

               fa_point = fa_point + cell%n_pt1

             end do

             fa_point = fa_point + &
                  (cell%n_pt2-finish2+start2-1)*cell%n_pt1

          end do

       end do
!$OMP END PARALLEL DO

    endif

    if (dot_prod%iscmplx) then
       dot_prod%z = dot_prod_complex
    else
       dot_prod%d = dot_prod_real
    end if

    call timer_clock('basis_dot_function_with_box', 2)

  end subroutine basis_dot_function_with_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 subroutine basis_multiply_function_by_box(fa_on_grid, &           ! in-out
       fa_box, fa_sphere, fa_tightbox, fa_start1, &                ! input
       fa_start2, fa_start3, ngwf_offset, cell, fftbox)            ! input

    !================================================================!
    ! This subroutines multiplies the values of a function in PPDs   !
    ! by a function held in a box.                                   !
    ! The PPDs which belong to the NGWF "fa" are multiplied          !
    ! and stored in the order dictated by the sphere of "fa".        !
    !----------------------------------------------------------------!
    ! Adapted by Nicholas Hine from basis_extract_function_from_box  !
    ! on 28/02/2011.                                                 !
    ! Modified by Andrea Greco on 13/04/2015 to make it compatible   !
    ! with new FUNCTIONS and FFTBOX_DATA types                       !
    !================================================================!

    use datatypes, only: FUNCTIONS, FFTBOX_DATA
    use fft_box, only: FFTBOX_INFO
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock

    implicit none

    ! Arguments
    ! integer, intent(in) :: box_n1,box_n2,box_n3
    type(SPHERE), intent(in) :: fa_sphere ! sphere of NGWFs PPDs we want
    type(FUNCTION_TIGHT_BOX), intent(in) :: fa_tightbox ! tightbox "
    integer, intent(in) :: fa_start1 ! 1-beginning of region wrt box
    integer, intent(in) :: fa_start2 ! 2-beginning of region wrt box
    integer, intent(in) :: fa_start3 ! 3-beginning of region wrt box
    integer, intent(in) :: ngwf_offset ! beginning of NGWF in PPD rep
    ! real(kind=DP), intent(in) :: fa_box(box_n1,box_n2,box_n3) ! box
    type(FFTBOX_DATA), intent(in) :: fa_box ! box
    ! real(kind=DP), intent(inout), dimension(:) :: fa_on_grid ! multiplied NGWF
    type(FUNCTIONS), intent(inout) :: fa_on_grid ! multiplied NGWF
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox


    ! Local Variables
    integer :: fa_ppd     ! PPD counter within sphere
    integer :: fa_point   ! point counter within fa_on_grid
    integer :: a1_pos     ! absolute 1-coordinate of current PPD
    integer :: a2_pos     ! absolute 2-coordinate of current PPD
    integer :: a3_pos     ! absolute 3-coordinate of current PPD
    integer :: start1     ! start of points within tightbox along direction 1 of current PPD
    integer :: start2     ! start of points within tightbox along direction 2 of current PPD
    integer :: start3     ! start of points within tightbox along direction 3 of current PPD
    integer :: finish1    ! end of points within tightbox along direction 1 of current PPD
    integer :: finish2    ! end of points within tightbox along direction 2 of current PPD
    integer :: finish3    ! end of points within tightbox along direction 3 of current PPD
    integer :: box_offset1 ! offset to express 1-start of function in box
    integer :: box_offset2 ! offset to express 2-start of function in box
    integer :: box_offset3 ! offset to express 3-start of function in box
    integer :: box_pt1    ! point counter wrt box along direction-1
    integer :: box_pt2    ! point counter wrt box along direction-2
    integer :: box_pt3    ! point counter wrt box along direction-3
    integer :: box_start1 ! start point in box along direction-1
    integer :: box_end1   ! end point in box along direction-1
    integer :: point2     ! point counter wrt current PPD along direction-2
    integer :: point3     ! point counter wrt current PPD along direction-3

    call timer_clock('basis_multiply_function_by_box', 1)

    if (cell%n_pts == 1) then
       ! cks: case of only one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               cell%n_ppds_a1, cell%n_ppds_a2, cell%n_ppds_a3, fftbox)

          box_pt1 = a1_pos - fa_tightbox%start_ppds1 + fa_start1
          box_pt2 = a2_pos - fa_tightbox%start_ppds2 + fa_start2
          box_pt3 = a3_pos - fa_tightbox%start_ppds3 + fa_start3

          fa_point = ngwf_offset + fa_ppd -1

          ! complex-complex case
          if (fa_on_grid%iscmplx.and.fa_box%iscmplx) then
             fa_on_grid%z(fa_point) = fa_on_grid%z(fa_point) * &
                  fa_box%z(box_pt1, box_pt2, box_pt3)
          ! complex - real case
          else if (fa_on_grid%iscmplx.and.(.not.fa_box%iscmplx)) then
             fa_on_grid%z(fa_point) = fa_on_grid%z(fa_point) * &
                  fa_box%d(box_pt1, box_pt2, box_pt3)
          ! real - complex case
          else if ((.not.fa_on_grid%iscmplx).and.fa_box%iscmplx) then
             fa_on_grid%d(fa_point) = fa_on_grid%d(fa_point) * &
                  real(fa_box%z(box_pt1, box_pt2, box_pt3))
          ! real - real
          else
             fa_on_grid%d(fa_point) = fa_on_grid%d(fa_point) * &
                  fa_box%d(box_pt1, box_pt2, box_pt3)
          endif

       enddo

    else
       ! cks: case for more than one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               cell%n_ppds_a1, cell%n_ppds_a2, cell%n_ppds_a3, fftbox)

          call basis_lims_1d_in_ppd_in_tight(start1, finish1, box_offset1, &
               a1_pos, fa_tightbox%start_ppds1, fa_tightbox%finish_ppds1, &
               fa_tightbox%start_pts1, fa_tightbox%finish_pts1, cell%n_pt1)

          call basis_lims_1d_in_ppd_in_tight(start2, finish2, box_offset2, &
               a2_pos, fa_tightbox%start_ppds2, fa_tightbox%finish_ppds2, &
               fa_tightbox%start_pts2, fa_tightbox%finish_pts2, cell%n_pt2)

          call basis_lims_1d_in_ppd_in_tight(start3, finish3, box_offset3, &
               a3_pos, fa_tightbox%start_ppds3, fa_tightbox%finish_ppds3, &
               fa_tightbox%start_pts3, fa_tightbox%finish_pts3, cell%n_pt3)

          fa_point = ngwf_offset + (fa_ppd-1)*cell%n_pts - 1 &
               + (start3-1)*cell%n_pt2*cell%n_pt1 &
               + (start2-1)*cell%n_pt1 + start1

          box_start1 = box_offset1 + fa_start1
          box_end1 = box_start1 + finish1 - start1

          do point3=start3,finish3
             ! cks: coordinates of points with respect to the FFT box,
             !      i.e. the first point inside the FFT box in a certain
             !      dimension starts from 1, etc.
             box_pt3 = point3 - start3 + box_offset3 + fa_start3
             do point2=start2,finish2

                box_pt2 = point2 - start2 + box_offset2 + fa_start2

                ! agrecocmplx
                ! complex-complex case
                if (fa_on_grid%iscmplx.and.fa_box%iscmplx) then
                   fa_on_grid%z(fa_point:fa_point+(finish1-start1)) = &
                     fa_on_grid%z(fa_point:fa_point+(finish1-start1)) * &
                     fa_box%z(box_start1:box_end1,box_pt2,box_pt3)

                ! complex - real case
                else if (fa_on_grid%iscmplx.and.(.not.fa_box%iscmplx)) then
                   fa_on_grid%z(fa_point:fa_point+(finish1-start1)) = &
                     fa_on_grid%z(fa_point:fa_point+(finish1-start1)) * &
                     fa_box%d(box_start1:box_end1,box_pt2,box_pt3)

                ! real - complex case
                else if ((.not.fa_on_grid%iscmplx).and.fa_box%iscmplx) then
                   fa_on_grid%d(fa_point:fa_point+(finish1-start1)) = &
                     fa_on_grid%d(fa_point:fa_point+(finish1-start1)) * &
                     real(fa_box%z(box_start1:box_end1,box_pt2,box_pt3),kind=DP)

                ! real - real
                else
                   fa_on_grid%d(fa_point:fa_point+(finish1-start1)) = &
                     fa_on_grid%d(fa_point:fa_point+(finish1-start1)) * &
                     fa_box%d(box_start1:box_end1,box_pt2,box_pt3)

                endif

                fa_point = fa_point + cell%n_pt1

             end do

             fa_point = fa_point + &
                  (cell%n_pt2-finish2+start2-1)*cell%n_pt1

          end do

       end do

    endif

    call timer_clock('basis_multiply_function_by_box', 2)

  end subroutine basis_multiply_function_by_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_clean_function(functions_on_grid,  &    ! inout
       function_sphere, cell, fftbox, opt_offset_override) ! in

    !===============================================================!
    ! This subroutine cleans (shaves) an NGWF in ppd representation !
    ! so that it is zero outside its sphere.                        !
    !---------------------------------------------------------------!
    ! @docme                                                        !
    ! opt_offset_override (in, opt): Optional override of the offset!
    !                                used to set point_counter.     !
    !---------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 5/7/2001.                 !
    ! Modified to use cell by Quintin Hill on 15/10/2008.           !
    ! Tidied up and moved to basis_mod by Nick Hine on 12/07/2009   !
    ! Modified to remove shaving for extended NGWFs by Andrea Greco !
    ! Modified to add optional offset override by Jacek Dziedzic.   !
    ! Modified by Andrea Greco on 07/04/2015 to make it compatible  !
    ! with new FUNCTIONS type                                       !
    !===============================================================!

    use constants, only: DP, PI
    use datatypes, only: FUNCTIONS
    use fft_box, only: FFTBOX_INFO
    use geometry, only: operator(*), point, local_displacement, operator(+), &
         geometry_distance, operator(-), operator(.DOT.)
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    ! integer, intent(in)          :: n_func_ppds
    type(CELL_INFO), intent(in)    :: cell
    type(FFTBOX_INFO), intent(in)  :: fftbox
    ! real(kind=DP), intent(inout) :: functions_on_grid(n_func_ppds*cell%n_pts)
    type(FUNCTIONS), intent(inout) :: functions_on_grid
    type(SPHERE), intent(in)       :: function_sphere
    integer, intent(in), optional  :: opt_offset_override

    ! Local Variables
    integer :: point_counter
    integer :: ppd
    integer :: loc
    integer :: loc_1, loc_2, loc_3
    integer :: a1_neighbour, a2_neighbour, a3_neighbour
    integer :: a1_upper,a2_upper,a3_upper
    integer :: a1_lower,a2_lower,a3_lower
    integer :: a1_cell,a2_cell,a3_cell
    integer :: ppd_count
    logical :: in_region
    real(kind=DP) :: local_1, local_2, local_3
    type(POINT) :: current_point,current_displacement,periodic_centre,ppd_origin
    real(kind=DP) :: recip_twopi
    real(kind=DP) :: fcoord_cp(3)
    real(kind=DP) :: fcoord_centre(3)

    recip_twopi = 0.5_DP / pi

    ! ndmh: to check for ppds included multiple times, when cell=fftbox
    ! agreco: if NGWFs are extended along a given direction, no need
    ! to check for NGWFs in periodic images of the simulation cell
    if ((fftbox%coin1).and.(.not.function_sphere%extended(1))) then
       a1_lower = -1
       a1_upper = 1
    else
       a1_lower = 0
       a1_upper = 0
    end if
    if ((fftbox%coin2).and.(.not.function_sphere%extended(2))) then
       a2_lower = -1
       a2_upper = 1
    else
       a2_lower = 0
       a2_upper = 0
    end if
    if ((fftbox%coin3).and.(.not.function_sphere%extended(3))) then
       a3_lower = -1
       a3_upper = 1
    else
       a3_lower = 0
       a3_upper = 0
    end if

    if(present(opt_offset_override)) then
       point_counter = opt_offset_override-1
    else
       point_counter = function_sphere%offset-1
    end if

    ! cks: loop over the ppds that belong to the sphere
    do ppd_count=1,function_sphere%n_ppds_sphere

       ppd = function_sphere%ppd_list(1,ppd_count)
       ppd_origin = basis_ppd_location(ppd,cell)

       loc = function_sphere%ppd_list(2,ppd_count)

       ! cks: loop over points of this ppd
       do loc_3=0,cell%n_pt3 -1
          local_3 = real(loc_3, DP)*cell%d3

          do loc_2=0,cell%n_pt2-1
             local_2 = real(loc_2, DP)*cell%d2

             do loc_1=0,cell%n_pt1-1
                local_1 = real(loc_1, DP)*cell%d1

                current_displacement = local_displacement(cell%a1_unit,&
                     cell%a2_unit,cell%a3_unit,local_1,local_2,local_3)

                current_point = ppd_origin + current_displacement

                ! agreco: remove contributions from directions
                ! along which NGWFs are extended
                fcoord_cp(1) = current_point .dot. cell%b1
                fcoord_cp(2) = current_point .dot. cell%b2
                fcoord_cp(3) = current_point .dot. cell%b3
                fcoord_cp = fcoord_cp * recip_twopi
                fcoord_cp = modulo(fcoord_cp,1.0_DP)

                if (function_sphere%extended(1)) then
                    current_point = current_point - fcoord_cp(1)*cell%a1
                end if
                if (function_sphere%extended(2)) then
                    current_point = current_point - fcoord_cp(2)*cell%a2
                end if
                if (function_sphere%extended(3)) then
                    current_point = current_point - fcoord_cp(3)*cell%a3
                end if

                point_counter = point_counter + 1

                in_region = .false.

                ! ndmh: loop over possible origins of this NGWF from different cells
                ! ndmh: NB: a1_lower = a1_upper = 0 if fftbox%coin1=.false.
                do a1_cell=a1_lower,a1_upper
                   do a2_cell=a2_lower,a2_upper
                      do a3_cell=a3_lower,a3_upper

                         a1_neighbour=nint(real(loc,kind=DP)/9.0_DP)
                         a2_neighbour=nint(real(loc-9*a1_neighbour,DP)/3.0_DP)
                         a3_neighbour=loc-9*a1_neighbour-3*a2_neighbour

                         ! ndmh: In cases where FFTbox==cell, check over nearest
                         ! ndmh: neighbours regardless of loc value
                         if (fftbox%coin1) a1_neighbour = a1_cell
                         if (fftbox%coin2) a2_neighbour = a2_cell
                         if (fftbox%coin3) a3_neighbour = a3_cell

                         periodic_centre = function_sphere%centre &
                              + real(a1_neighbour,DP)*cell%a1 &
                              + real(a2_neighbour,DP)*cell%a2 &
                              + real(a3_neighbour,DP)*cell%a3

                         ! agreco: remove contributions from directions
                         ! along which NGWFs are extended
                         fcoord_centre(1) = periodic_centre .dot. cell%b1
                         fcoord_centre(2) = periodic_centre .dot. cell%b2
                         fcoord_centre(3) = periodic_centre .dot. cell%b3
                         fcoord_centre = fcoord_centre * recip_twopi
                         fcoord_centre = modulo(fcoord_centre,1.0_DP)

                         if (function_sphere%extended(1)) then
                            periodic_centre = periodic_centre - fcoord_centre(1)*cell%a1
                         end if
                         if (function_sphere%extended(2)) then
                            periodic_centre = periodic_centre - fcoord_centre(2)*cell%a2
                         end if
                         if (function_sphere%extended(3)) then
                            periodic_centre = periodic_centre - fcoord_centre(3)*cell%a3
                         end if

                         if (geometry_distance(periodic_centre,current_point) <= &
                              function_sphere%radius) then
                            in_region = .true.
                         end if

                      enddo
                   enddo
                enddo

                if (.not.in_region) then
                   if (functions_on_grid%iscmplx) then
                      functions_on_grid%z(point_counter) = 0.0_DP
                   else
                      functions_on_grid%d(point_counter) = 0.0_DP
                   endif
                endif

             enddo
          enddo
       enddo

    enddo

  end subroutine basis_clean_function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  integer function basis_count_psincs(function_sphere,cell,fftbox)

    !===============================================================!
    ! This subroutine counts the number of unique psinc values in   !
    ! this function sphere.                                         !
    !---------------------------------------------------------------!
    ! Written by Nicholas Hine on 11/11/2010.                       !
    !===============================================================!

    use constants, only: DP
    use fft_box, only: FFTBOX_INFO
    use geometry, only: operator(*), point, local_displacement, operator(+), &
         geometry_distance
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    type(SPHERE), intent(in) :: function_sphere
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local Variables
    integer :: ppd
    integer :: loc
    integer :: loc_1, loc_2, loc_3
    integer :: a1_neighbour, a2_neighbour, a3_neighbour
    integer :: a1_upper,a2_upper,a3_upper
    integer :: a1_lower,a2_lower,a3_lower
    integer :: a1_cell,a2_cell,a3_cell
    integer :: ppd_count
    logical :: in_region
    real(kind=DP) :: local_1, local_2, local_3
    type(POINT) :: current_point,current_displacement,periodic_centre,ppd_origin

    ! ndmh: to check for ppds included multiple times, when cell=fftbox
    if (fftbox%coin1) then
       a1_lower = -1
       a1_upper = 1
    else
       a1_lower = 0
       a1_upper = 0
    end if
    if (fftbox%coin2) then
       a2_lower = -1
       a2_upper = 1
    else
       a2_lower = 0
       a2_upper = 0
    end if
    if (fftbox%coin3) then
       a3_lower = -1
       a3_upper = 1
    else
       a3_lower = 0
       a3_upper = 0
    end if

    basis_count_psincs = 0

    ! cks: loop over the ppds that belong to the sphere
    do ppd_count=1,function_sphere%n_ppds_sphere

       ppd = function_sphere%ppd_list(1,ppd_count)
       ppd_origin = basis_ppd_location(ppd,cell)

       loc = function_sphere%ppd_list(2,ppd_count)

       ! cks: loop over points of this ppd
       do loc_3=0,cell%n_pt3 -1
          local_3 = real(loc_3, DP)*cell%d3

          do loc_2=0,cell%n_pt2-1
             local_2 = real(loc_2, DP)*cell%d2

             do loc_1=0,cell%n_pt1-1
                local_1 = real(loc_1, DP)*cell%d1

                current_displacement = local_displacement(cell%a1_unit,&
                     cell%a2_unit,cell%a3_unit,local_1,local_2,local_3)

                current_point = ppd_origin + current_displacement

                in_region = .false.

                ! ndmh: loop over possible origins of this NGWF from different cells
                ! ndmh: NB: a1_lower = a1_upper = 0 if fftbox%coin1=.false.
                do a1_cell=a1_lower,a1_upper
                   do a2_cell=a2_lower,a2_upper
                      do a3_cell=a3_lower,a3_upper

                         a1_neighbour=nint(real(loc,kind=DP)/9.0_DP)
                         a2_neighbour=nint(real(loc-9*a1_neighbour,DP)/3.0_DP)
                         a3_neighbour=loc-9*a1_neighbour-3*a2_neighbour
                         if (fftbox%coin1) a1_neighbour = a1_cell
                         if (fftbox%coin2) a2_neighbour = a2_cell
                         if (fftbox%coin3) a3_neighbour = a3_cell

                         periodic_centre = function_sphere%centre &
                              + real(a1_neighbour,DP)*cell%a1 &
                              + real(a2_neighbour,DP)*cell%a2 &
                              + real(a3_neighbour,DP)*cell%a3

                         if (geometry_distance(periodic_centre,current_point) <= &
                              function_sphere%radius) then
                            in_region = .true.
                         end if

                      enddo
                   enddo
                enddo

                if (in_region) basis_count_psincs = basis_count_psincs + 1

             enddo
          enddo
       enddo

    enddo

  end function basis_count_psincs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_put_tightbox_in_fftbox(fftbox_inout, &   ! input/output
       start1, start2, start3, &                            ! input
       tightbox_in, npts1, npts2, npts3, factor)            ! input

    !================================================!
    ! Deposits tightbox data into an FFTbox.         !
    !------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/2/2004. !
    ! Modified by Andrea Greco on 07/04/2015 to make !
    ! compatible with new FFTBOX_DATA type           !
    !================================================!

    use datatypes, only: FFTBOX_DATA
    use fft_box, only: FFTBOX_INFO
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    integer, intent(in)        :: start1, start2, start3
    integer, intent(in)        :: npts1, npts2, npts3
    real(kind=DP), intent(in)    :: factor
    ! pdh: changed to avoid unnecessary copy onto stack
    ! real(kind=DP), intent(in)    :: tightbox_in(:,:,:)
    type(FFTBOX_DATA), intent(in) :: tightbox_in
    ! type(FFTBOX_INFO), intent(in) :: fftbox
    ! real(kind=DP), intent(inout) :: fftbox_inout(fftbox%total_ld1, &
    !      fftbox%total_ld2, fftbox%total_pt3)
    type(FFTBOX_DATA), intent(inout) :: fftbox_inout

    ! cks: << local variables>>
    integer :: end1
    integer :: end2
    integer :: end3

    call timer_clock('basis_put_tightbox_in_fftbox',1)

    ! agreco: check we copy real to real or complex to complex
    call utils_assert(fftbox_inout%iscmplx.eqv.tightbox_in%iscmplx,'Error in &
       &basis_put_tightbox_in_fftbox: incompatible argument types')

    end1 = start1 + npts1 - 1
    end2 = start2 + npts2 - 1
    end3 = start3 + npts3 - 1

    ! cks: deposit tightbox (multiplied by factor) to the appropriate
    ! cks: position in fftbox_inout
    if (fftbox_inout%iscmplx) then
       fftbox_inout%z(start1: end1, start2: end2, start3: end3) = &
          fftbox_inout%z(start1: end1, start2: end2, start3: end3) &
         +factor *tightbox_in%z(1: npts1, 1: npts2, 1: npts3)
    else
       fftbox_inout%d(start1: end1, start2: end2, start3: end3) = &
          fftbox_inout%d(start1: end1, start2: end2, start3: end3) &
         +factor *tightbox_in%d(1: npts1, 1: npts2, 1: npts3)
    endif

    call timer_clock('basis_put_tightbox_in_fftbox',2)

  end subroutine basis_put_tightbox_in_fftbox


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_copy_tightbox_to_fftbox(box_out, &
       start1, start2, start3, &
       data_in, npts1, npts2, npts3)

    !===========================================!
    ! Copies tightbox into FFTbox.              !
    !-------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2001. !
    ! Modified by Andrea Greco on 13/04/2015 to !
    ! make it compatible with new FUNCTIONS and !
    ! FFTBOX_DATA types                         !
    !===========================================!

    use datatypes, only: FFTBOX_DATA, data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    integer, intent(in)        :: start1,start2,start3
    integer, intent(in)        :: npts1,npts2,npts3
    ! real(kind=DP), intent(in)  :: data_in(:,:,:)
    type(FFTBOX_DATA), intent(in) :: data_in
    ! type(FFTBOX_INFO), intent(in) :: fftbox
    ! real(kind=DP), intent(out) :: box_out(fftbox%total_ld1, &
    !     fftbox%total_ld2, fftbox%total_pt3)
    type(FFTBOX_DATA), intent(inout) :: box_out

    call timer_clock('basis_copy_tightbox_to_fftbox',1)

    ! agreco: check we copy real to real or complex to complex
    call utils_assert(data_in%iscmplx.eqv.box_out%iscmplx,'Error in &
       &basis_copy_tightbox_to_fftbox: incompatible argument types')

    ! cks: initialise pair_data_row
    call data_set_to_zero(box_out)

    ! cks: copy all the elements of data_in to the appropriate
    !      position in box_out
    if (box_out%iscmplx) then
       box_out%z(start1:start1+npts1-1,start2:start2+npts2-1, &
         start3:start3+npts3-1) = data_in%z(1:npts1, 1:npts2, 1:npts3)
    else
       box_out%d(start1:start1+npts1-1,start2:start2+npts2-1, &
         start3:start3+npts3-1) = data_in%d(1:npts1, 1:npts2, 1:npts3)
    end if

    call timer_clock('basis_copy_tightbox_to_fftbox',2)

  end subroutine basis_copy_tightbox_to_fftbox


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_lims_1d_in_ppd_in_tight(start, finish, box_offset, &
       pos, start_ppds, finish_ppds, start_pts, finish_pts, n_pt)

    !==================================================================!
    ! This subroutine works in 1d and given a tightbox and a ppd       !
    ! inside it returns the starting and ending grid points of the ppd !
    ! within the tightbox and the number of points of the tightbox     !
    ! before the starting grid point of the current ppd.               !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2001, improved in 2007.      !
    !==================================================================!

    ! Arguments
    implicit none
    integer, intent(out) :: start      ! first point (of ppd) in tightbox wrt start of ppd
    integer, intent(out) :: finish     ! last point (of ppd) in tightbox wrt start of ppd
    integer, intent(out) :: box_offset ! number of points in tightbox belonging to previous ppds

    integer, intent(in) :: pos         ! position of current ppd
    integer, intent(in) :: start_ppds  ! first ppd of tightbox
    integer, intent(in) :: finish_ppds ! last ppd of tightbox
    integer, intent(in) :: start_pts   ! first point in tightbox of current ppd
    integer, intent(in) :: finish_pts  ! last point in tightbox of current ppd
    integer, intent(in) :: n_pt        ! number of points in current ppd

    if (pos == start_ppds .and. pos == finish_ppds) then
       ! cks: the box contains only 1 ppd in this dimension
       start = start_pts
       finish = finish_pts
       box_offset = 0
    else if (pos == start_ppds) then
       ! cks: this is the first ppd in the box in this dimension
       start = start_pts
       finish = n_pt
       box_offset = 0
    else if (pos == finish_ppds) then
       ! cks: this is the last ppd in the box in this dimension
       start = 1
       finish = finish_pts
       box_offset = -start_pts + 1 + (finish_ppds-start_ppds) * n_pt
    else
       ! cks: this is neither the first, nor the last ppd in the fftbox in
       !      this dimension
       start = 1
       finish = n_pt
       box_offset = -start_pts + 1 + (pos-start_ppds) * n_pt
    end if

  end subroutine basis_lims_1d_in_ppd_in_tight


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &  ! output
       global_count, global_loc, n_ppds_a1, n_ppds_a2, n_ppds_a3, & ! input
       fftbox)                                                      ! input

    !==========================================================!
    ! This subroutine returns the integer PPD-grid coordinates !
    ! of the current PPD along each lattice vector direction.  !
    ! These coordinates are "absolute" in the sense that they  !
    ! belong to PPDs which may not even be inside the          !
    ! simulation cell so that they always belong to an NGWF    !
    ! sphere whose centre is inside the simulation cell.       !
    !----------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2001.                !
    ! Modified by Chris-Kriton Skylaris on 21/06/2007 to work  !
    ! also when simulation cell and FFT box coincide in one or !
    ! more dimensions.                                         !
    ! Modified by Nicholas Hine in June 2008 for optimal speed !
    !==========================================================!

    use fft_box, only: FFTBOX_INFO

    implicit none

    ! Arguments
    integer, intent(out) :: a1_pos, a2_pos, a3_pos
    integer, intent(in) :: global_count, global_loc, &
         n_ppds_a1, n_ppds_a2, n_ppds_a3
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! cks: internal variable declarations
    integer ::  map_a1a2, a1_neighbour, a2_neighbour, a3_neighbour
    integer, parameter :: neighbours(-13:13,1:3) = &
         reshape((/-1, -1, -1, -1, -1, -1, -1, -1, -1, &
         0,  0,  0,  0,  0,  0,  0,  0,  0, &
         +1, +1, +1, +1, +1, +1, +1, +1, +1, &
         -1, -1, -1,  0,  0,  0, +1, +1, +1, &
         -1, -1, -1,  0,  0,  0, +1, +1, +1, &
         -1, -1, -1,  0,  0,  0, +1, +1, +1, &
         -1,  0, +1, -1,  0, +1, -1,  0, +1, &
         -1,  0, +1, -1,  0, +1, -1,  0, +1, &
         -1,  0, +1, -1,  0, +1, -1,  0, +1 /),(/27,3/))

    ! These are the coordinates of the current ppd in the simulation cell
    a3_pos = (global_count-1) / (n_ppds_a1*n_ppds_a2) + 1
    map_a1a2 = global_count - (n_ppds_a1*n_ppds_a2) * (a3_pos - 1)

    a2_pos = (map_a1a2 - 1) / n_ppds_a1 + 1

    a1_pos = map_a1a2 - n_ppds_a1 * (a2_pos-1)

    ! If the current ppd is not due to a function in the simulation cell
    ! but due to a function in a neighbouring simulation cell
    ! (periodic image), find the coordinates of the ppd that corresponds
    ! to the actual function inside the simulation cell and therefore
    ! is a ppd outside the simulation cell.

    ! ndmh: use lookup table to avoid floating-point arithmetic since this
    ! ndmh: function is called millions of times per iteration
    a1_neighbour = neighbours(global_loc, 1)
    a2_neighbour = neighbours(global_loc, 2)
    a3_neighbour = neighbours(global_loc, 3)

    ! cks: these are now the positions of ppds that can be even outside the
    ! cks: simulation cell, but they all belong to a single function whose
    ! cks: centre is INSIDE the simulation cell.
    if (.not.fftbox%coin3) a3_pos = a3_pos - a3_neighbour * n_ppds_a3
    if (.not.fftbox%coin2) a2_pos = a2_pos - a2_neighbour * n_ppds_a2
    if (.not.fftbox%coin1) a1_pos = a1_pos - a1_neighbour * n_ppds_a1

  end subroutine basis_find_ppd_in_neighbour

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module basis
