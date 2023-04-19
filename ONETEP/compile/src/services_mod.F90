! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Chris-Kriton Skylaris, Arash A. Mostofi, Peter D. Haynes, Nicholas Hine,
!   Simon Dubois, Robert Bell, Gilberto Teobaldi and James C. Womack
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module services

  use constants, only: DP

  implicit none

  private

  public :: services_symm_fd_sec_der
  public :: services_symm_fd_sec_der_KAX
  public :: services_equally_spaced_numbers
  public :: services_locate_interp
  public :: services_linear_interpolation
  public :: services_1d_interpolation
  public :: services_polak_cg_coeff
  public :: services_line_search_parabola
  public :: services_parabolic_step
  public :: services_polynomial_step
  public :: services_cubic_fit_minimum
!  public :: services_cubic_minimum
  public :: services_cubic_fit_maximum
  public :: services_flush
  public :: services_sbessj
  public :: services_radial_transform
  public :: services_regular_transform
  public :: services_radial_integral
  public :: services_radial_integral_rmax
  public :: services_radial_derivative
  public :: services_regular_integral
  public :: services_analytic_limit
  public :: services_read_xyz
  public :: services_write_xyz
  public :: services_open_pspot_file
  public :: services_rationalise_coords
  public :: services_maxboltzdist
  public :: services_maxboltzdist_parallel
  public :: services_gammadist
  public :: services_rms_fit
  public :: services_incomplete_gamma_prepare_lookup
  public :: services_incomplete_gamma_destroy_lookup
  public :: services_incomplete_gamma_from_lookup
  public :: services_gradient_on_grid
  public :: services_factorial

  interface services_polynomial_step
     module procedure services_polynomial_step_real
     module procedure services_polynomial_step_complex
  end interface services_polynomial_step

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine services_write_xyz(elements, xyz_root, title_line)

    !==================================================================!
    ! This subroutine outputs the cartesian coordinates of the atoms   !
    ! in a .xyz format file in Angstroms.                              !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/03/2007.                  !
    !==================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, ANGSTROM
    use ion, only: ELEMENT
    use rundat, only: pub_debug_on_root
    use utils, only: utils_unit, utils_assert

    implicit none

    type(ELEMENT), intent(in) :: elements(:)   ! elements for all atoms on all procs
    character(len=*), intent(in) :: xyz_root   ! root name of file to output
    character(len=*), intent(in) :: title_line ! title line to output in file

    ! Local Variables
    integer :: xyz_unit                ! unit for output file
    integer :: ierr                    ! error flag
    integer :: row                     ! atom counter
    character(len=128) :: xyz_filename ! name of file to output
    character(len=10)  :: cbuf         ! character buffer
    integer :: nat

    if (pub_debug_on_root) &
         write(stdout,'(a)') 'DEBUG: Entering services_write_xyz'

    if (pub_on_root) then

       write(xyz_filename,'(2a)') trim(xyz_root), ".xyz"

       write(stdout,'(3a)',advance='no') &
       ' Writing xyz coordinates to file "',trim(xyz_filename),'" ...'


       ! cks: find free unit to open for the write operation
       xyz_unit =utils_unit()

       ! cks: open file
       open(unit=xyz_unit, form="formatted", file= trim(xyz_filename), &
            action="write", iostat=ierr, position='append')
       call utils_assert(ierr == 0, 'Error in services_write_xyz(): &
               &opening file "'//trim(xyz_filename)//'" failed with code ',ierr)

       ! cks: write the number of atoms
       nat = size(elements)
       write(cbuf,'(i9)') nat
       write(xyz_unit,'(a)',iostat=ierr) adjustl(cbuf)
       call utils_assert(ierr == 0, 'Error in services_write_xyz(): &
               &writing nat failed.')

       ! cks: write the title line
       write(xyz_unit,'(a)',iostat=ierr) adjustl(title_line)
       call utils_assert(ierr == 0, 'Error in services_write_xyz(): &
               &writing title_line failed.')

       ! cks: write the atomic coordinates
       do row =1, nat
          write(xyz_unit,'(a,3f12.6)',iostat=ierr)adjustl(elements(row)%symbol), &
               elements(row)%centre%x/ANGSTROM, &
               elements(row)%centre%y/ANGSTROM, &
               elements(row)%centre%z/ANGSTROM
          call utils_assert(ierr == 0, 'Error in services_write_xyz(): &
                  &writing elements(row)%centre%x/ANGSTROM failed.')
       end do

       ! close file
       close(unit=xyz_unit, iostat=ierr)
       call utils_assert(ierr == 0, 'Error in services_write_xyz(): &
               &closing file "'//trim(xyz_filename)//'" failed with code ',ierr)

       write(stdout,'(a)') ' done'

    endif


    if (pub_debug_on_root) &
         write(stdout,'(a)') 'DEBUG: Leaving services_write_xyz'

    return
  end subroutine services_write_xyz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine services_read_xyz(xyz_file,nat,positions,species)

    !==========================================================!
    ! This subroutine reads position data stored in the xyz    !
    ! file format.                                             !
    ! There are two branches to this routine:                  !
    ! 1) no elements array is provided in the call, and only   !
    !    the number of atoms stored in the file are counted;   !
    ! 2) posiitons array is provided in the call, and the full !
    !    positional data is read and returned.                 !
    !                                                          !
    ! The routine should always be called twice:               !
    !   call services_read_xyz(xyz_file,nat)                   !
    !   allocate(elements(nat))                                !
    !   call services_read_xyz(xyz_file,nat,nsp,pos,species)   !
    !----------------------------------------------------------!
    ! Written by Robert Bell, 22/07/2013                       !
    !==========================================================!

    use comms, only: pub_on_root
    use constants, only:  ANGSTROM
    use utils, only: utils_unit, utils_open_unit_check, utils_close_unit_check,&
         utils_assert, utils_abort

    implicit none

    ! arguments
    character(len=*), intent(in) :: xyz_file
    integer, intent(inout)       :: nat
    real(kind=DP), intent(inout), optional :: positions(:,:) ! dim: (3,nat)
    character(len=4), intent(out), optional :: species(:)    ! dim: (nat)

    ! internal
    integer :: iunit
    integer, parameter :: llength = 80
    integer :: iat, ierr, dummy_nat
    character(len=llength) :: cbuff

    ! root does all the work
    if (pub_on_root) then
       ! open the file
       iunit = utils_unit()
       open(unit=iunit,file=trim(xyz_file),action='read',iostat=ierr)
       call utils_open_unit_check('services_read_xyz',trim(xyz_file),ierr)


       ! first branch: count the number of atoms in the file
       if ((.not.present(positions)).and. (.not.present(species))) then
         read(iunit,*,err=100) nat
       else ! second branch: fill the elements array with the positional data

         ! read nat...
         read(iunit,*,err=100) dummy_nat
         call utils_assert (nat == dummy_nat, 'Error in services_read_xyz():&
              &Inconsistent number of atoms in file "'//trim(xyz_file)//'"&
              & Read i, expecting j atoms. Values of i and j follow ', dummy_nat, nat)
         ! ...and comment line
         read(iunit,'(a)')

         ! ...and position data
         do iat=1,nat
           read(iunit,*,err=200,iostat=ierr) cbuff, positions(:,iat)
           call utils_assert(ierr >= 0, "End of file "//trim(xyz_file)//&
                " reached but only i atoms of n read. Values of i and n &
                &follow: ", iat-1, nat)
           cbuff = adjustl(cbuff)
           species(iat) = cbuff(1:4)
         enddo

         ! convert positions to Bohr
         positions(:,:) = positions(:,:) * ANGSTROM
       endif

       ! tidy up
       close(iunit,iostat=ierr)
       call utils_close_unit_check('services_read_xyz',trim(xyz_file),ierr)
    endif

    return

 100 call utils_abort("Error in services_read_xyz(): Error reading number of &
          &atoms in file "//trim(xyz_file))
 200 call utils_abort("Error in services_read_xyz(): Error reading line n &
          &in file "//trim(xyz_file)//". The value of n was ", iat+2)
  end subroutine services_read_xyz



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_open_pspot_file(iunit,filename,ierr)

    !========================================================================!
    ! This subroutine tries to open a pseudopotential file: first it tries   !
    ! the current directory. If the file does not exist, it then checks the  !
    ! environment variable PSPOT_DIR and looks for tbe file there (ACCELRYS  !
    ! version only).                                                         !
    !------------------------------------------------------------------------!
    ! Arguments                                                              !
    !   iunit (in): unit number to open file with.                           !
    !   filename (in): name of the pseudopotential file to try to open       !
    !   file_exists (out): can find pseudopotential                          !
    !   ierr (out): error flag                                               !
    !------------------------------------------------------------------------!
    ! Written by Victor Milman in May 2010.                                  !
    ! Cleaned up and formatted for ONETEP by Nicholas Hine in June 2010.     !
    ! Moved to services_mod and reworked by Nicholas Hine in June 2018.      !
    !========================================================================!

    use rundat, only: pub_pseudo_path
    use utils, only: utils_abort, utils_open_unit_check

    implicit none

    ! Arguments
    integer,intent(in) :: iunit
    character(len=*) :: filename
    integer,intent(out) :: ierr

    ! Local Variables
    logical :: file_exists

    ! rab207: initialise ierr to zero, in case file does not exist
    ierr = 0

    ! Check if file exists in current directory
    inquire(file=trim(filename),exist=file_exists)

    ! If so, open it and return
    if (file_exists) then
       open(iunit,file=trim(filename),status='old',position='rewind', &
            iostat=ierr)
    else
       ! Check if file exists once pub_pseudo_path is appended
       inquire(file=trim(pub_pseudo_path)//'/'//trim(filename),exist=file_exists)

       ! If so, open it
       if (file_exists) then
          open(iunit,file=trim(pub_pseudo_path)//'/'//trim(filename),status='old', &
               position='rewind',iostat=ierr)
       else
          call utils_abort('ERROR in services_open_pspot_file: Cannot find &
               &file "'//trim(filename)//'".')
       end if
    endif

    call utils_open_unit_check('services_open_pspot_file',trim(filename),ierr)

  end subroutine services_open_pspot_file


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function services_sbessj(l,x)

    !===============================================!
    ! SPHERICAL BESSEL function OF THE FIRST KIND.  !
    !-----------------------------------------------!
    ! Written by Peter D. Haynes in early 2004.     !
    ! Fixed by Jacek Dziedzic in 2014/07 to avoid   !
    ! division overflow for tiny x.                 !
    !===============================================!

    use constants, only: SAFE_DIV_EPS

    implicit none

    ! Arguments
    integer, intent(in) :: l
    real(kind=DP), intent(in) :: x

    ! Local variables
    integer :: j
    real(kind=DP), parameter :: third = 1.0_DP / 3.0_DP
    real(kind=DP), parameter :: ftnth = 1.0_DP / 14.0_DP
    real(kind=DP) :: x2, sb0, sb1, by, bym, byp, ux
    real(kind=DP) :: sbessj

    x2 = 0.5_DP*x*x
    if (abs(x) > 0.001_DP) then
       sb0 = sin(x)/x
    else
       sb0 = 1.0_DP - third*x2*(1.0_DP - 0.1_DP*x2)
    end if
    if (l == 0) then
       sbessj = sb0
    else
       if (abs(x) > 0.001_DP) then
          sb1 = (sb0 - cos(x)) / x
       else
          sb1 = third*x*(1.0_DP - (0.2_DP*x2)*(1.0_DP - ftnth*x2))
       end if
       if (l == 1) then
          sbessj = sb1
       else if (abs(x) < SAFE_DIV_EPS) then
          sbessj = 0.0_DP
       else
          by = sb1
          bym = sb0
          ux = 1.0_DP / x
          do j=1,l-1
             byp = real(2*J+1,DP)*ux*by - bym
             bym = by
             by = byp
          end do
          sbessj = by
       end if
    end if

    services_sbessj =sbessj

  end function services_sbessj


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! cks: written by Chris-Kriton Skylaris in 2000
  function services_symm_fd_sec_der(val,delta_space,order)

    use utils, only: utils_abort

    implicit none

    real(kind=DP) :: services_symm_fd_sec_der

    ! Arguments
    integer, intent(in) :: order
    real(kind=DP), intent(in) :: val(-order/2:order/2), delta_space

    ! cks: internal variables
    real(kind=DP) :: sd


    if (order.eq.2) then
       sd= val(-1) -2.0_DP*val(0) +val(1)
    elseif (order.eq.4) then
       sd= -val(-2)/12.0_DP + 4.0_DP*val(-1)/3.0_DP  -5.0_DP*val(0)/2.0_DP &
            +4.0_DP*val(1)/3.0_DP -val(2)/12.0_DP
    elseif (order.eq.6) then
       sd= ( val(-3) + val(3) )/90.0_DP -3.0_DP*(val(-2) + val(2) )/20.0_DP &
            +3.0_DP*(val(-1)+val(1))/2.0_DP -49.0_DP*val(0)/18.0_DP
    elseif (order.eq.8) then
       sd= -(val(-4)+val(4))/560.0_DP +8.0_DP*(val(-3)+val(3))/315.0_DP &
            -(val(-2)+val(2))/5.0_DP +8.0_DP*(val(-1)+val(1))/5.0_DP &
            -205.0_DP*val(0)/72.0_DP
    elseif (order.eq.10) then
       sd= (val(-5)+val(5))/3150.0_DP -5.0_DP*(val(-4)+val(4))/1008.0_DP &
            +5.0_DP*(val(-3)+val(3))/126.0_DP -5.0_DP*(val(-2)+val(2))/21.0_DP &
            +5.0_DP*(val(-1)+val(1))/3.0_DP -5269.0_DP*val(0)/1800.0_DP
    elseif (order.eq.12) then
       sd= -(val(-6)+val(6))/16632.0_DP + 2.0_DP*(val(-5)+val(5))/1925.0_DP &
            -(val(-4)+val(4))/112.0_DP +10.0_DP*(val(-3)+val(3))/189.0_DP &
            -15.0_DP*(val(-2)+val(2))/56.0_DP + 12.0_DP*(val(-1)+val(1))/7.0_DP &
            -5369.0_DP*val(0)/1800.0_DP
    elseif (order.eq.16) then
       sd= -(val(-8)+val(8))*2.428127428127428127428127428127428E-6_DP &
            +(val(-7)+val(7))*5.074290788576502862217147931433646E-5_DP &
            -(val(-6)+val(6))*5.180005180005180005180005180005181E-4_DP &
            +(val(-5)+val(5))*3.480963480963480963480963480987078E-3_DP &
            -(val(-4)+val(4))*1.767676767676767676767676767673481E-2_DP &
            +(val(-3)+val(3))*7.542087542087542087542087542076068E-2_DP &
            -(val(-2)+val(2))*0.311111111111111111111111111111111_DP  &
            +(val(-1)+val(1))*1.777777777777777777777777777777777_DP   &
            -val(0)*3.05484410430839005820287374945639_DP
    elseif (order.eq.20) then
       sd= -(val(-10)+val(10))*1.082508822446902942258979410682197E-7_DP &
            +(val(-9)+val(9))*2.672861289992352943849331878227647E-6_DP &
            -(val(-8)+val(8))*3.213698066639243109831345125462773E-5_DP &
            +(val(-7)+val(7))*2.518489913447896641173952098321847E-4_DP &
            -(val(-6)+val(6))*1.456876456876456876456876456876457E-3_DP &
            +(val(-5)+val(5))*6.713286713286713286713286713332222E-3_DP &
            -(val(-4)+val(4))*2.622377622377622377622377622392322E-2_DP &
            +(val(-3)+val(3))*9.324009324009324009324009324009325E-2_DP &
            -(val(-2)+val(2))*0.3409090909090909090909090909090909_DP   &
            +(val(-1)+val(1))*1.8181818181818181818181818181818181_DP   &
            -val(0)*3.09953546233308141622756510748108_DP
    elseif (order.eq.24) then
       sd= -(val(-12)+val(12))*5.136127090629715478281907141780610E-9_DP &
            +(val(-11)+val(11))*1.466979770679032784540683560495355E-7_DP &
            -(val(-10)+val(10))*2.041302350899874119688361174429285E-6_DP &
            +(val(-9)+val(9))*1.848092663366141178318680898655856E-5_DP &
            -(val(-8)+val(8))*1.227970945463205525125029768992618E-4_DP &
            +(val(-7)+val(7))*6.415521674256747233306277976777757E-4_DP &
            -(val(-6)+val(6))*2.765208647561589252051467324785477E-3_DP &
            +(val(-5)+val(5))*1.023917259211376858435681965093730E-2_DP &
            -(val(-4)+val(4))*3.399725274725274725274725274737918E-2_DP &
            +(val(-3)+val(3))*0.107448107448107448107448107448107_DP &
            -(val(-2)+val(2))*0.362637362637362637362637362637_DP &
            +(val(-1)+val(1))*1.846153846153846153846153846153_DP &
            -val(0)*3.1299532768418050158602556492717_DP
    elseif (order.eq.28) then
       sd = -(val(-14)+val(14))*2.543605797264240046387230203548493E-10_DP &
            +(val(-13)+val(13))*8.259945926263993712765159382884103E-9_DP &
            -(val(-12)+val(12))*1.308685182692451503866229939725700E-7_DP &
            +(val(-11)+val(11))*1.349784386777007832086822284940225E-6_DP &
            -(val(-10)+val(10))*1.020774442500112173015659352986046E-5_DP &
            +(val(-9)+val(9))*6.049033733333998062315018388065456E-5_DP &
            -(val(-8)+val(8))*2.934726522187822497420020639834882E-4_DP &
            +(val(-7)+val(7))*1.204692403277100313809734420083823E-3_DP &
            -(val(-6)+val(6))*4.304265565875473796285376112133183E-3_DP &
            +(val(-5)+val(5))*1.377364981080151358789129686981574E-2_DP &
            -(val(-4)+val(4))*4.089052287581699346405228758192856E-2_DP &
            +(val(-3)+val(3))*0.118954248366013083531987019838242_DP &
            -(val(-2)+val(2))*0.379166666666666666666666666666666_DP &
            +(val(-1)+val(1))*1.866666666666666666666666666666666_DP &
            -val(0)*3.15199167800108529601965668779366_DP
    else
       call utils_abort("Error in services_symm_fd_sec_der: Finite difference &
            &second derivative not available for order ",order)
       sd = 0.0_DP ! qoh: Prevent compiler warning
    endif

    services_symm_fd_sec_der=sd/(delta_space**2)

  end function services_symm_fd_sec_der


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! cks: written by Chris-Kriton Skylaris, 9/11/2001
  function services_symm_fd_sec_der_KAX(val,delta_space,order)

    use constants, only: PI
    use utils, only: utils_abort

    implicit none

    real(kind=DP) :: services_symm_fd_sec_der_KAX

    ! Arguments
    integer, intent(in) :: order
    real(kind=DP), intent(in) :: val(-order/2:order/2), delta_space

    ! cks: internal variables
    real(kind=DP) :: sd
    real(kind=DP), parameter :: psq=PI**2

    if (order.eq.12) then
       sd=  (val(-6)+val(6))*(1.0_DP/600.0_DP-psq/4096.0_DP) &
            + (val(-5)+val(5))*( 3.0_DP*psq/1024.0_DP - 31.0_DP/1575.0_DP  )  &
            +(val(-4)+val(4))*( 2647.0_DP/25200.0_DP-33.0_DP*psq/2048.0_DP ) &
            +(val(-3)+val(3))*( 55.0_DP*psq/1024.0_DP - 103.0_DP/315.0_DP ) &
            +(val(-2)+val(2))*( 493.0_DP/840.0_DP - 495.0_DP*psq/4096.0_DP ) &
            +(val(-1)+val(1))*( 26.0_DP/75.0_DP + 99.0_DP*psq/512.0_DP ) &
            -val(0)*(2497.0_DP/1800.0_DP+231.0_DP*psq/1024.0_DP )
    else
       call utils_abort("Error in services_symm_fd_sec_der_KAX: Finite &
            &difference second derivative not available for order ",order)
       sd = 0.0_DP ! qoh: Prevent compiler warning
    endif


    services_symm_fd_sec_der_KAX=sd/(delta_space**2)

  end function services_symm_fd_sec_der_KAX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_equally_spaced_numbers(numbers,min_number,max_number, &
       num_numbers)

    use utils, only: utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: num_numbers
    real(kind=DP), intent(out) :: numbers(num_numbers)
    real(kind=DP), intent(in) :: min_number, max_number

    ! Local Variables
    integer :: row
    real(kind=DP) :: delta_number

    call utils_assert(num_numbers>1, 'Error in services_equally_spaced_&
         &numbers(): num_numbers invalid.')

    delta_number = (max_number-min_number)/real(num_numbers-1,kind=DP)

    do row=1,num_numbers
       numbers(row) = min_number + real(row-1,kind=DP)*delta_number
    end do

  end subroutine services_equally_spaced_numbers


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  integer function services_locate_interp(length,points,num_points)

    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: num_points
    real(kind=DP), intent(in) :: length
    real(kind=DP), intent(in) :: points(num_points)

    ! Local Variables
    integer :: n_low,n_high,n_mid
    real(kind=DP) :: eps_tol = 1D-12

    ! CKS: FIRST USE BISECTION METHOD TO FIND THE POINTS BETWEEN WHICH LENGTH LIES
    n_low=0
    n_high=num_points+1
    do
       if (n_high-n_low==1) exit

       n_mid=(n_low+n_high)/2
       if ( length>=points(n_mid) ) then
          n_low=n_mid
       else
          n_high=n_mid
       endif

    enddo

    ! CKS: TAKE CARE OF THE CASE WHERE LENGTH COINCIDES WITH THE FIRST OR THE
    ! CKS: LAST POINT OR STOP IF IT IS OUT OF BOUNDS.
    if (abs(length-points(1)) < eps_tol) then
       services_locate_interp=1
    else if (abs(length-points(num_points)) < eps_tol) then
       services_locate_interp=num_points-1
    else
       services_locate_interp=n_low
       if ((n_low==0) .or. (n_low==num_points)) then
          call utils_abort('Error in services_locate_interp: length out of &
               &interpolating points range')
       endif

    endif

  end function services_locate_interp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function services_linear_interpolation(x_i,y1,y2,x1,x2)

    ! CKS : 1D REAL LINEAR INTERPOLATION BETWEEN TWO GIVEN POINTS

    use constants, only: SAFE_DIV_EPS
    use utils, only: utils_abort

    implicit none


    ! Arguments
    real(kind=DP), intent(in) :: x_i,y1,y2,x1,x2

    if (abs(x2-x1) < SAFE_DIV_EPS) then
       call utils_abort('Error in services_linear_interpolation: Attempted &
            &division by (almost) zero', opt_real_to_print1=x2-x1)
    end if

    services_linear_interpolation = ( y1*(x2-x_i) + y2*(x_i-x1) )/ (x2-x1)


  end function services_linear_interpolation



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function services_1d_interpolation(values, num_values, xx, l_mom)

    !===============================================================!
    ! Interpolation in one dimension suitable for pseudopotentials. !
    !---------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 10/9/2001.                !
    ! Modified by Chris-Kriton Skylaris on 29/3/2004.               !
    ! Cleanup by Nicholas Hine, May 2011.                           !
    ! Additional documentation by James C. Womack, May 2017.        !
    !===============================================================!

    use utils, only: utils_abort

    implicit none

    ! Arguments
    real(kind=DP) :: services_1d_interpolation
    integer, intent(in) ::  num_values
    integer, intent(in) ::  l_mom
    real(kind=DP), intent(in) :: values(num_values)
    real(kind=DP), intent(in) :: xx

    ! Local Variables
    integer :: start_index
    real(kind=DP) :: x2
    real(kind=DP) :: off, f0, f1, f2, f3, f4, t0, t1, t2, t3, t4

    ! initialisations
    start_index = int(xx) + 1
    x2 = real(start_index-1,DP)
    off = xx - x2

    if (start_index < 1) then
       call utils_abort('Error in services_1d_interpolation(): start_index=', &
            start_index)
       services_1d_interpolation = 0.0_DP !qoh: Prevent compiler warning

    elseif (start_index .eq. 1) then
       ! JCW: t0 .. t4 are the polynomial coefficients when the equation
       ! JCW:   fn(x) = t0 + t1*x + t2*x^2 + t3*x^3 + t4*x^4
       ! JCW: is fitted to 5 equally spaced points
       ! JCW:    fn(-2) = f0 = t0 - 2*t1 + 4*t2 - 8*t3 + 16*t4,
       ! JCW:    fn(-1) = f1 = t0 - t1 + t2 - t3 + t4,
       ! JCW:    fn(0)  = f2 = t0,
       ! JCW:    fn(1)  = f3 = t0 + t1 + t2 + t3 + t4,
       ! JCW:    fn(2)  = f4 = t0 + 2*t1 + 4*t2 + 8*t3 + 16*t4
       ! JCW: Using Sage Math to solve this set of linear equations, i.e.
       ! JCW:   sys = [ f0 == fn(-2), f1 == fn(-1), f2 == fn(0), f3 == fn(1), f4 == fn(2)
       ! JCW:   solve( sys, t0, t1, t2, t3, t4 )
       ! JCW: we obtain
       ! JCW: [[
       ! JCW:   t0 == f2,
       ! JCW:   t1 == 1/12*f0 - 2/3*f1 + 2/3*f3 - 1/12*f4,
       ! JCW:   t2 == -1/24*f0 + 2/3*f1 - 5/4*f2 + 2/3*f3 - 1/24*f4,
       ! JCW:   t3 == -1/12*f0 + 1/6*f1 - 1/6*f3 + 1/12*f4,
       ! JCW:   t4 == 1/24*f0 - 1/6*f1 + 1/4*f2 - 1/6*f3 + 1/24*f4
       ! JCW: ]]
       ! JCW: Since we are on point 1 (fn(0)), and we do not know fn(-2), fn(-1), we
       ! JCW: assume that the function is either symmetric or antisymmetric
       ! JCW: about the origin (based on the value of l_mom)

       ! cks: take into account whether function is even or odd
       f0 = values(3)*( (-1.0_DP)**l_mom)
       f1 = values(2)*( (-1.0_DP)**l_mom)
       f2 = values(1)
       f3 = values(2)
       f4 = values(3)

       ! Quartic interpolation
       t0 = f2
       t1 = 2.0_dp*f0-16.0_dp*f1+16.0_dp*f3-2.0_dp*f4
       t2 = -f0+16.0_dp*f1-30.0_dp*f2+16.0_dp*f3-f4
       t3 = -2.0_dp*f0+4.0_dp*f1-4.0_dp*f3+2.0_dp*f4
       t4 = f0-4.0_dp*f1+6.0_dp*f2-4.0_dp*f3+f4

       services_1d_interpolation = t0+off*(t1+off*(t2+off*(t3+off*t4)))/24.0_dp

    elseif (start_index.gt.(num_values-2)) then

       services_1d_interpolation = 0.0_DP

    else
       ! JCW: t0, t1, t2, t3 are the polynomial coefficients when the equation
       ! JCW:    fn(x) = t0 + t1*x + t2*x^2 + t3*x^3
       ! JCW: is fitted to 4 equally spaced points, x = -1, 0, 1, 2, i.e.
       ! JCW:    fn(-1) = f1 = t0 - t1 + t2 - t3
       ! JCW:    fn(0)  = f2 = t0
       ! JCW:    fn(1)  = f3 = t0 + t1 + t2 + t3
       ! JCW:    fn(2)  = f4 = t0 + 2*t1 + 4*t2 + 8*t3
       ! JCW: Using Sage Math to solve this set of linear equations, i.e.
       ! JCW:   solve( [ f1 == fn(-1), f2 == fn(0), f3 == fn(1), f4 == fn(2) ], t0, t1, t2, t3 )
       ! JCW: we obtain
       ! JCW: [[
       ! JCW:   t0 == f2,
       ! JCW:   t1 == -1/3*f1 - 1/2*f2 + f3 - 1/6*f4,
       ! JCW:   t2 == 1/2*f1 - f2 + 1/2*f3,
       ! JCW:   t3 == -1/6*f1 + 1/2*f2 - 1/2*f3 + 1/6*f4
       ! JCW: ]]
       f1 = values(start_index-1)
       f2 = values(start_index)
       f3 = values(start_index+1)
       f4 = values(start_index+2)

       t0 = f2
       t1 = ((6.0_dp*f3)-(2.0_dp*f1)-(3.0_dp*f2)-f4)/6.0_dp
       t2 = (f1+f3-(2.0_dp*f2))/2.0_dp
       t3 = (f4-f1+(3.0_dp*(f2-f3)))/6.0_dp

       services_1d_interpolation = t0+off*(t1+off*(t2+off*t3))

    endif


  end function services_1d_interpolation



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  function services_polak_cg_coeff(prev_cov_dir,&
       cov_grad,contra_grad,prev_contra_grad, vec_size)

    !===============================================================!
    ! This subroutine calculates the conjugate gradient coefficient !
    ! according to the Polak formula:                               !
    !  b_{r+1}=g_{r+1}*(g_{r+1}-g_r)/( p_r * (g_{r+1}-g_r) )        !
    ! taking into account the covariant and contravariant           !
    ! character of the quantities involved.                         !
    !---------------------------------------------------------------!
    !                                                               !
    !                                                               !
    !---------------------------------------------------------------!
    ! Originaly written by Chris-Kriton Skylaris on 25/6/2001       !
    ! for the ONES code.                                            !
    ! Slightly modified by Arash A. Mostofi in March 2003.          !
    ! Parallelised by Chris-Kriton Skylaris on 30/12/2003.          !
    !===============================================================!

    use datatypes, only: COEF, FUNCTIONS, data_coef_abs, data_set_to_zero
    use comms, only: pub_on_root, comms_reduce
    use constants, only: DP, stdout, NORMAL
    use rundat, only: pub_output_detail
    use utils, only: utils_assert
    implicit none

    type(COEF)                  :: services_polak_cg_coeff
    integer,         intent(in) :: vec_size
    type(FUNCTIONS), intent(in) :: prev_cov_dir
    type(FUNCTIONS), intent(in) :: cov_grad
    type(FUNCTIONS), intent(in) :: contra_grad
    type(FUNCTIONS), intent(in) :: prev_contra_grad

    type(COEF) :: denominator  !jmecmplx
    real(kind=DP) :: eps
    logical :: loc_cmplx

    !jmecmplx BEGIN
    ! agrecocmplx: I suppose we either have all complex or
    ! all real only?
    loc_cmplx = prev_cov_dir%iscmplx
    ! agrecocmplx: check this is the case otherwise the following check
    ! will not work as expected...
    call utils_assert((prev_cov_dir%iscmplx .eqv. cov_grad%iscmplx).and.&
         (cov_grad%iscmplx .eqv. contra_grad%iscmplx).and.&
         (contra_grad%iscmplx .eqv. prev_contra_grad%iscmplx), &
         'Error in services_polak_cg_coeff: incompatible argument types')
    ! agrecocmplx: should be compatible now, but need to check it's ok....
    !call utils_assert(.not. loc_cmplx, 'Error in services_polak_cg_coeff: &
    !     &function only ready for real NGWFs.')
    ! agrecocmplx
    ! is this check implemented correctly? I would expect we want all the arrays
    ! to be of the same size?
    if (loc_cmplx) then
       call utils_assert( (size(prev_cov_dir%z) /= vec_size) .or. &
            (size(cov_grad%z) /= vec_size) .or. (size(contra_grad%z) /= vec_size) &
            .or. (size(prev_contra_grad%z) /= vec_size), &
            'Error in services_polak_cg_coeff: sizes of arrays do not match.')
    else
       call utils_assert( (size(prev_cov_dir%d) /= vec_size) .or. &
            (size(cov_grad%d) /= vec_size) .or. (size(contra_grad%d) /= vec_size) &
            .or. (size(prev_contra_grad%d) /= vec_size), &
            'Error in services_polak_cg_coeff: sizes of arrays do not match.')
    end if
    services_polak_cg_coeff%iscmplx = loc_cmplx
    denominator%iscmplx = loc_cmplx
    !jmecmplx END

    eps =epsilon(1.0_DP)

    ! cks: calculate my proc contribution
    ! agrecocmplx: the formula should be the same even in complex case, but need
    ! to take the complex conjugate of left term when making the product?
    if (loc_cmplx) then
       denominator%z = sum( conjg(prev_cov_dir%z(1: vec_size)) &                           !jmecmplx
            *( contra_grad%z(1: vec_size) - prev_contra_grad%z(1: vec_size) )   )   !jmecmplx
    else
       denominator%d = sum( prev_cov_dir%d(1: vec_size) &                           !jmecmplx
            *( contra_grad%d(1: vec_size) - prev_contra_grad%d(1: vec_size) )   )   !jmecmplx
    end if
    ! cks: add up contributions from all procs
    if (loc_cmplx) then
       call comms_reduce('SUM', denominator%z)
    else
       call comms_reduce('SUM', denominator%d)
    end if


    if ( data_coef_abs(denominator) .gt. eps ) then           !jmecmplx
       ! cks: contribution from my proc
       ! agrecocmplx: need to take complex conjugate of left term?
       if (loc_cmplx) then
          services_polak_cg_coeff%z = &                !jmecmplx
               sum(conjg(cov_grad%z(1: vec_size)) &           !jmecmplx
               *(contra_grad%z(1: vec_size) - prev_contra_grad%z(1: vec_size) ) ) & !jmecmplx
               / denominator%z                         !jmecmplx
       else
          services_polak_cg_coeff%d = &                !jmecmplx
               sum(cov_grad%d(1: vec_size) &           !jmecmplx
               *(contra_grad%d(1: vec_size) - prev_contra_grad%d(1: vec_size) ) ) & !jmecmplx
               / denominator%d                         !jmecmplx
       end if
       ! cks: sum of contributions from all procs
       if (loc_cmplx) then
          call comms_reduce('SUM', services_polak_cg_coeff%z)  !jmecmplx
       else
          call comms_reduce('SUM', services_polak_cg_coeff%d)  !jmecmplx
       end if
    else
       if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
            write(stdout,'(a/a)') &
            'WARNING: zero denominator in services_polak_cg_coeff', &
            '         Setting to zero'
       !services_polak_cg_coeff%d = 0.0_DP                   !jmecmplx
       ! agrecocmplx
       call data_set_to_zero(services_polak_cg_coeff)
    endif

    ! agrecocmplx
    if (loc_cmplx) then
       if ( abs(services_polak_cg_coeff%z).gt.(2.0_DP) ) then  !jmecmplx

          if (pub_on_root .and. pub_output_detail >= NORMAL) &
               write(stdout,'(a,f11.5,a)') &
               'WARNING: services_polak_cg_coeff too large &
               &(',services_polak_cg_coeff%z,'). Setting to zero'  !jmecmplx
          services_polak_cg_coeff%z=(0.0_DP,0.0_DP)                !jmecmplx
       endif
    else
       if ( abs(services_polak_cg_coeff%d).gt.(2.0_DP) ) then  !jmecmplx

          if (pub_on_root .and. pub_output_detail >= NORMAL) &
               write(stdout,'(a,f11.5,a)') &
               'WARNING: services_polak_cg_coeff too large &
               &(',services_polak_cg_coeff%d,'). Setting to zero'  !jmecmplx
          services_polak_cg_coeff%d=0.0_DP                         !jmecmplx
       endif
    end if


  end function services_polak_cg_coeff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! LINE SEARCH SERVICES !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_cubic_fit_minimum( &
       min_length, min_value, success, &         ! output
       val1, val2, val3, grad1, coor2, coor3, &  ! input
       max_step)

    !==================================================================!
    ! Do line search by fitting cubic polynomial to f0, f1, f2 and g0. !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 4/12/2001 based on code      !
    ! written earlier by Peter D. Haynes.                              !
    ! Modified by Chris-Kriton Skylaris on 16/7/2003.                  !
    ! Return value to indicate success added by Nick Hine on 04/10/2008!
    !==================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout, NORMAL
    use rundat, only: pub_output_detail

    implicit none

    ! Arguments
    real(kind=DP), intent(out):: min_length, min_value
    logical, intent(out) :: success
    real(kind=DP), intent(in) :: val1, val2, val3, grad1, coor2, coor3
    real(kind=DP), intent(in) :: max_step
    ! coor1 is zero

    ! Local Variables
    real(kind=DP) :: aa, bb, xx, yy, disc, aon3b

    ! cks: initialise
    min_value =0.0_DP ; xx = coor2 ; yy = coor3

    aa = ((yy*yy*yy*val2 - xx*xx*xx*val3)/(yy-xx) &
         - (xx*xx+xx*yy+yy*yy)*val1 - xx*yy*(xx+yy)*grad1) / (xx*xx*yy*yy)

    bb = ((yy*yy*val2 - xx*xx*val3)/(xx-yy) + (xx+yy)*val1 &
         + xx*yy*grad1) / (xx*xx*yy*yy)

    disc=-1.0_DP

    if (abs(bb*yy/aa) > epsilon(1.0_DP)) then        ! avoid div by zero
       aon3b = aa / (3.0_DP * bb)
       disc = aon3b * aon3b - grad1 / (3.0_DP * bb)  ! discriminant

       if (disc >= 0.0_DP) then

          ! cks: set optimal_step to cubic minimum
          min_length = -aon3b + sign(sqrt(disc), bb)

          ! cks: Value of cubic polynomial at minimum
          min_value = val1 + min_length * (grad1 + min_length * &
               (aa + min_length * bb))

       end if
    end if

    ! cks: see if the cubic fit was successful and choose safe length if not.
    if (disc < 0.0_DP ) then
       if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
            write(stdout,'(a)') 'WARNING: Cubic fit was unsuccessful.'
       min_length=0.15_DP
       success = .false.
    ! cks: protection from crazy line search coefficients
    ! ndmh: now uses max_step input parameter
    else if (min_length > max_step) then
       if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
            write(stdout,'(a/a,f11.5,a)') &
            'WARNING in services_cubic_fit_minimum:','  cubic step (', &
            min_length,') too large - setting to safe value'
       min_length = 0.10_DP
       success = .false.
    else if (min_length < 0.0_DP) then
       if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
            write(stdout,'(a/a,f11.5,a)') &
            'WARNING in services_cubic_fit_minimum:','  cubic step (', &
            min_length,') less than zero - setting to safe value'
       min_length = 0.10_DP
       success = .false.
    else
       success = .true.
    endif

  end subroutine services_cubic_fit_minimum


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_line_search_parabola( &
       min_length, min_value, success, &                       ! output
       grad_at_zero, value_at_zero, value_at_trial_length, &   ! input
       trial_length,max_step)                                  ! input

    !==================================================================!
    ! Do line search by fitting a parabola.                            !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/6/2001.                   !
    ! Improved by Chris-Kriton Skylaris on 20/6/2003 and on 26/4/2004. !
    ! Return value to indicate success added by Nick Hine on 04/10/2008!
    !==================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout, NORMAL
    use rundat, only: pub_output_detail
    implicit none

    real(kind=DP), intent(out) :: min_length, min_value
    logical, intent(out) :: success
    real(kind=DP), intent(in) :: trial_length, grad_at_zero, value_at_zero &
         , value_at_trial_length,max_step

    ! cks: internal declarations
    real(kind=DP) :: linear_term
    real(kind=DP) :: quadratic_term
    real(kind=DP) :: pscale
    real(kind=DP) :: q_diff


    ! cks: determine pscale so that precision is maintained in polynomial fitting
    pscale =abs(value_at_trial_length -value_at_zero)
    if ( pscale > tiny(1.0_DP)) then
       pscale =1.0_DP/pscale
       if (value_at_trial_length > value_at_zero) then
          q_diff =1.0_DP
       else
          q_diff =-1.0_DP
       endif
    else
       pscale =100000.0_DP
       q_diff =pscale*value_at_trial_length -pscale*value_at_zero
    endif


    linear_term =  pscale*grad_at_zero  ! cks: b

    quadratic_term = q_diff -linear_term*trial_length


    if ( abs(  quadratic_term/ (pscale*value_at_zero)  ) > epsilon(1.0_DP)) then

       min_length = -0.5_DP *linear_term /quadratic_term ! -b/2a
       min_length = min_length *(trial_length**2)


       min_value = value_at_zero &
            -(0.25_DP * linear_term*linear_term/ (pscale*quadratic_term) ) &
            *(trial_length**2)  ! c -b^2/4a
       success = .true.
    else
       min_length = 0.1_DP
       success = .false.
    endif

    ! cks: protection from crazy line search coefficients
    ! ndmh: changed to use keyword for max line search step
    if (abs(min_length) > max_step) then
       if (pub_on_root .and. pub_output_detail >= NORMAL) then
          write(stdout,'(a/a,f11.5,a)') &
               'WARNING in services_line_search_parabola:', &
               '  quadratic step (', min_length, &
               ') too large - setting to safe value'
          write(stdout,'(a,f11.5,a)') &
               'value at zero =', value_at_zero
          write(stdout,'(a,f11.5,a)') &
               'linear_term=', linear_term
          write(stdout,'(a,f11.5,a)') &
               'quadratic_term=', quadratic_term
          write(stdout,'(a,f11.5,a)') &
               'trial_length=', trial_length
          write(stdout,'(a,f11.5,a)') &
               'pscale=', pscale
       end if
       min_length = 0.15_DP
       success = .false.
    end if


  end subroutine services_line_search_parabola


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_fit_parabola_ptsgrad(aa,bb,cc &
       ,grad_at_zero, value_at_zero, value_at_point, point)

    !=========================================================!
    ! Fits a parabola using its value and gradient at point 0 !
    ! and its value at point x.                               !
    !---------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/6/2001           !
    !=========================================================!

    implicit none

    real(kind=DP), intent(out) :: aa, bb, cc
    real(kind=DP), intent(in) :: grad_at_zero, value_at_zero &
         ,value_at_point, point

    ! cks: aa= ( f(0+x)-f(0)-f'(0)*x ) / ( x^2 )
    aa=(value_at_point - (value_at_zero + grad_at_zero*point) ) / ( point**2 )

    ! cks: bb= f'(0)
    bb=grad_at_zero

    ! cks: cc=f(0)
    cc=value_at_zero

  end subroutine services_fit_parabola_ptsgrad



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  function services_parabolic_step(Qinit, Q1, Q2, L1, L2)

    !==========================================================!
    ! Fits a parabola given value of function at three points. !
    !----------------------------------------------------------!
    ! Written by Arash A. Mostofi, December 2002.              !
    !==========================================================!

    use comms, only : pub_on_root
    use constants, only: stdout
    implicit none

    real(kind=DP) :: services_parabolic_step
    real(kind=DP), intent(in) :: Qinit  ! initial value of function
    real(kind=DP), intent(in) :: Q1,Q2  ! value at trial steps 1 and 2
    real(kind=DP), intent(in) :: L1,L2  ! trial steps 1 and 2

    real(kind=DP) :: ll,aa,bb,cc,d1,d2,eps
    logical, parameter :: info=.false.

    eps = epsilon(1.0_DP)
    services_parabolic_step=10.0_DP
    cc=Qinit ; ll=L2/L1
    d1=Q1-Qinit ; d2 = Q2-Qinit

    bb = d1*ll - d2/ll ; aa = (d2 - d1*ll)/L2

    if (info .and. pub_on_root) then
       write(stdout,'(a,f22.15)') 'LINE SEARCH =====> quadratic coefficient =', aa
       write(stdout,'(a,f22.15)') 'LINE SEARCH =====> linear    coefficient =', bb
    endif

    if (abs(aa).gt.eps) then
       services_parabolic_step = -bb/(2*aa)
    else
       if (pub_on_root) write(stdout,'(a)') 'WARNING: Quadratic fit unsuccessful'
       services_parabolic_step = 0.5_DP
    endif

    !    if (abs(services_parabolic_step).gt.(1.0_DP)) then
    !       print*,'services_parabolic_step=' &
    !            ,services_parabolic_step
    !       print*,'setting services_parabolic_step equal to 0.07'
    !       services_parabolic_step=0.07_DP
    !    endif

  end function services_parabolic_step


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_polynomial_step_real(poly,poly_order,poly_num,Yi,Xi,Xj)

    !==========================================================!
    ! Fits a series of "poly_num" polynomial given values of   !
    ! "poly_num" functions at "poly_order" points.             !
    !----------------------------------------------------------!
    ! Written by Simon M.-M. Dubois, June 2011.                !
    !==========================================================!

    use utils, only: utils_assert

    implicit none

    integer, intent(in) :: poly_order                    ! polynomial order + 1
    integer, intent(in) :: poly_num                      ! number of polynomials
    real(kind=DP), intent(in) :: Yi(poly_order,poly_num) ! set of ordinates
    real(kind=DP), intent(in) :: Xi(poly_order)          ! set of abcissa
    real(kind=DP), intent(in) :: Xj                      ! interpolated abcissa
    real(kind=DP), intent(out) :: poly(poly_num)

    real(kind=DP) :: Vdm_mat(poly_order,poly_order)
    real(kind=DP) :: Vdm_vec(poly_order)
    !real(kind=DP) :: Yj(poly_num,1)
    integer :: ipiv(poly_order)
    integer :: info, io, jo

    ! LAPACK subroutine
    external :: dgesv

    ! Compose Vandermonde matrix
    do io = 1, poly_order
       do jo = 1, poly_order
          Vdm_mat(jo,io) = Xi(jo)**(poly_order-io)
       enddo
       Vdm_vec(io) = Xj**(poly_order-io)
    enddo

    ! Solve Vandermonde system of equations (i.e. compute polynomial coeffs)
    call dgesv(poly_order,poly_num,Vdm_mat,poly_order,ipiv,Yi,poly_order,info)
    call utils_assert(info == 0, 'Error in services_polynomial_step_real(): &
           &computation of the polynomial coefficients with &
           &lapack_dgesv failed with error, ', info)

    poly(:) = 0.0_dp
    do io = 1, poly_num
       do jo = 1, poly_order
          poly(io) = Yi(jo,io)*Vdm_vec(jo)
       enddo
    enddo
    !! Compute the interpolated ordinates
    !call dgemm('T','N',poly_num,1,poly_order,1.0_dp,Yi,poly_order,&
    !       Vdm_vec,1,0.0_dp,Yj,poly_num)

    !poly(:) = Yj(:,1)

  end subroutine services_polynomial_step_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine services_polynomial_step_complex(poly,poly_order,poly_num,Yi,Xi,Xj)

    !==============================================================!
    ! Fits a series of "poly_num" polynomial given values of       !
    ! "poly_num" functions at "poly_order" points.                 !
    !--------------------------------------------------------------!
    ! Adapted by J.M. Escartin from services_polynomial_step_real  !
    ! (by Simon M.-M. Dubois).                                     !
    !==============================================================!

    use constants, only: DP, cmplx_0, stdout
    use utils, only: utils_assert

    implicit none

    integer, intent(in) :: poly_order                    ! polynomial order + 1
    integer, intent(in) :: poly_num                      ! number of polynomials
    complex(kind=DP), intent(in) :: Yi(poly_order,poly_num) ! set of ordinates
    real(kind=DP), intent(in) :: Xi(poly_order)          ! set of abcissa
    real(kind=DP), intent(in) :: Xj                      ! interpolated abcissa
    complex(kind=DP), intent(out) :: poly(poly_num)

    complex(kind=DP) :: Vdm_mat(poly_order,poly_order)
    real(kind=DP) :: Vdm_vec(poly_order)
    !real(kind=DP) :: Yj(poly_num,1)
    integer :: ipiv(poly_order)
    integer :: info, io, jo

    ! LAPACK subroutine
    external :: zgesv

    ! Compose Vandermonde matrix
    do io = 1, poly_order
       do jo = 1, poly_order
          Vdm_mat(jo,io) = Xi(jo)**(poly_order-io)
       enddo
       Vdm_vec(io) = Xj**(poly_order-io)
    enddo

    ! Solve Vandermonde system of equations (i.e. compute polynomial coeffs)
    call zgesv(poly_order,poly_num,Vdm_mat,poly_order,ipiv,Yi,poly_order,info)
    call utils_assert(info == 0, 'Error in services_polynomial_step_complex(): &
           &computation of the polynomial coefficients with &
           &lapack_zgesv failed with error, ', info)

    poly(:) = cmplx_0
    do io = 1, poly_num
       do jo = 1, poly_order
          poly(io) = Yi(jo,io)*Vdm_vec(jo)
       enddo
    enddo
    !! Compute the interpolated ordinates
    !call dgemm('T','N',poly_num,1,poly_order,1.0_dp,Yi,poly_order,&
    !       Vdm_vec,1,0.0_dp,Yj,poly_num)

    !poly(:) = Yj(:,1)

  end subroutine services_polynomial_step_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#if 0
  ! jd (retreat2014): This function doesn't exactly do what's advertised --
  !                   in case of no minimum it returns 0.4, rather than stopping
  !                   There are direct float comparisons and no protection
  !                   against division overflow. Finally, it's never used
  !                   anywhere. Getting rid of it to reduce the no. of warnings
  function services_cubic_minimum(aa,bb,cc,dd)


    !==================================================================!
    ! Returns the x value at which the cubic polynomial                !
    ! f(x)=aa*x^3+bb*x^2+cc*x+dd has a minimium. If the polynomial has !
    ! no minimum the function stops.                                   !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 19/4/2001.                   !
    !==================================================================!


    use comms, only : pub_on_root
    use constants, only: stdout
    implicit none

    real(kind=DP) :: services_cubic_minimum
    real(kind=DP), intent(in) :: aa, bb, cc, dd

    ! cks: internal declarations
    real(kind=DP) :: x1, x2, fx1, fx2, discriminant

    discriminant=4.0_DP*(bb**2-3.0_DP*aa*cc)
    if (discriminant.lt.(0.0_DP)) then
       if (pub_on_root) write(stdout,'(a)')&
            'negative discriminant in services_cubic_minimum'
       !       print*,'ONES execution stops'
       !       stop
       if (pub_on_root) write(stdout,'(a)')&
            'RETURNING services_cubic_minimum=0.4_DP'
       services_cubic_minimum=0.4_DP
       return
    endif
    discriminant=sqrt(discriminant)

    if (aa.ne.(0.0_DP)) then ! jd: Unsafe
       x1=(-2.0_DP*bb+discriminant)/(6.0_DP*aa)
       x2=(-2.0_DP*bb-discriminant)/(6.0_DP*aa)
    else if (bb.ne.(0.0_DP)) then ! jd: Unsafe
       x1=-cc/(2.0_DP*bb)
       x2=x1
    else
       x1=0.0_DP
       x2=x1
    endif

    fx1=aa*x1**3+bb*x1**2+cc*x1+dd
    fx2=aa*x2**3+bb*x2**2+cc*x2+dd

    if (fx1.lt.fx2) then
       services_cubic_minimum=x1
    else
       services_cubic_minimum=x2
    endif

  end function services_cubic_minimum
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!gibo-start


!  subroutine services_cubic_fit_maximum( &
!       min_length, min_value, success, &         ! output
!       val1, val2, val3, grad1, coor2, coor3, &  ! input
!       max_step)

  subroutine services_cubic_fit_maximum( &
       max_length, max_value, success, &         ! output
       y0, y1, y2, g0, x1, x2, &  ! input
       max_step)

    !==================================================================!
    ! Do line search by fitting cubic polynomial to y0, y1, y2 and g0. !
    !          y(x) = ax^3 + bx^2 + cx + d                             !
    !      dy(x)/dx = 3ax^2 + 2b*x + c                                 !
    ! d^2y(x)/da^2x = 6ax + 2b                                         !
    !                                                                  !
    ! Mind that for computational convenience the x-frame is shifted   !
    ! to have x0= 0, yielding                                          !
    ! *** y(x0) = d                                                    !
    ! *** dy(x0)/dx = 3a(x0)^2 +2b(x0) +c = c                          !
    !------------------------------------------------------------------!
    ! Adapted from services_cubic_fit_minimum                          !
    ! by Gilberto Teobaldi on 9/12/11                                  !
    !==================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout, NORMAL
    use rundat, only: pub_output_detail

    implicit none

    ! Arguments
    !real(kind=DP), intent(out):: min_length, min_value
    real(kind=DP), intent(out):: max_length, max_value
    logical, intent(out) :: success
    !real(kind=DP), intent(in) :: val1, val2, val3, grad1, coor2, coor3
    real(kind=DP), intent(in) :: y0, y1, y2, g0, x1, x2
    real(kind=DP), intent(in) :: max_step
    ! coor1 is zero

    ! Local Variables
    !real(kind=DP) :: aa, bb, xx, yy, disc, aon3b
    real(kind=DP) :: a, b, c, d
    !real(kind=DP) :: x0, x1, x2
    real(kind=DP) :: x1_2, x1_3, x2_2
    real(kind=DP) :: x_stat_1, x_stat_2
    real(kind=DP) :: second_deriv_1, second_deriv_2
    !real(kind=DP) :: y0, y1, y2
    !real(kind=DP) :: g0
    real(kind=DP) :: den, num
    real(kind=DP) :: disc
    !real(kind=DP) :: max_value, max_length

  !initialise a negative discriminant for dy/dx= 3ax^2 +2bx + c = 0
  disc = -1.0_DP

  !====STEP-1:
  !fit cubic polynomial (y = a*x^3 + b*x^2 + c*x + d)
  !to (available) y(x0), y(x1), y(x2) and g0= [dy/dx]@x0
  !MIND: from computational convenience the x-frame is shifted
  !to have x0=0., therefore y(x0)=d

  x1_2 = x1*x1
  x1_3 = x1_2*x1

  x2_2 = x2*x2

  den = x2_2*x1 - 1.0_DP

  num = y2*x1_3 - y1 + y0 + g0*x1 -g0*x2*x1_3 - y0*x1_3

  NUM_STAB: if (abs(den) > epsilon(1.0_DP)) then  ! avoid division by zero
  !NUM_STAB: if (abs(den) > 0.) then  ! avoid division by zero

    a = (y1 - y0 - g0*x1 - num/den) / x1_2
    b = num /(x1_2*den)
    c = g0
    d = y0

    !====STEP-2
    ! Calculate the x-coordinate of the stationary points
    ! dy/dx= 3ax^2 +2bx + c = 0
    ! and check sign of 2nd derivative (6ax+2b) at stat.points
    ! [we are after a maximu i.e. d^2y/dx^2 < 0.

    ! discriminant of dy/dx= 3ax^2 +2bx + c = 0
    disc = 4.0_DP*b*b -12.0_DP*a*c
    !disc = 4.0*b*b -12.0*a*c

    if (disc >= 0.0_DP) then
    !if (disc >= 0.) then

        x_stat_1 = (-2.0_DP*b - disc)/( 6.0_DP*a)
        !x_stat_1 = (-2.0*b - disc)/( 6.0*a)

        second_deriv_1 = 6.0_DP*a*x_stat_1 +2.0_DP*b
        !second_deriv_1 = 6.0*a*x_stat_1 +2.0*b

        x_stat_2 = (-2.0_DP*b + disc)/( 6.0_DP*a)
        !x_stat_2 = (-2.0*b + disc)/( 6.0*a)

        second_deriv_2 = 6.0_DP*a*x_stat_2 +2.0_DP*b
        !second_deriv_2 = 6.0*a*x_stat_2 +2.0*b

        if (second_deriv_1 < 0.0_DP) then
        !if (second_deriv_1 < 0.0) then

           ! set optimum step (from x0) to reach cubic maximum
           max_length = x_stat_1

           ! value of cubic polynomial at cubic maximum
           max_value = a*x_stat_1*x_stat_1*x_stat_1 + b*x_stat_1*x_stat_1 &
                     + c*x_stat_1 + d

        elseif (second_deriv_2 > 0.0_DP) then
        !elseif (second_deriv_2 > 0.0) then

           ! set optimum step (from x0) to reach cubic maximum
           max_length = x_stat_2

           ! value of cubic polynomial at cubic maximum
           max_value = a*x_stat_2*x_stat_2*x_stat_2 + b*x_stat_2*x_stat_2 &
                     + c*x_stat_2 + d
        else

         if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
               write(stdout,'(a/a)') &
               'WARNING in services_cubic_fit_maximum:',&
               &' none of the (two) stationary point is a maximum!!!'

        endif

    endif

  endif NUM_STAB

  ! cks: see if the cubic fit was successful and choose safe length if not.
  if (disc < 0.0_DP ) then
     if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
          write(stdout,'(a)') 'WARNING: Cubic fit was unsuccessful.'
     !min_length=0.15_DP
     max_length=0.15_DP
     success = .false.
  ! cks: protection from crazy line search coefficients
  ! ndmh: now uses max_step input parameter
  !else if (min_length > max_step) then
  else if (max_length > max_step) then
     if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
          write(stdout,'(a/a,f11.5,a)') &
          'WARNING in services_cubic_fit_minimum:','  cubic step (', &
          !min_length,') too large - setting to safe value'
          max_length,') too large - setting to safe value'
     !min_length = 0.10_DP
     max_length = 0.10_DP
     success = .false.
  !else if (min_length < 0.0_DP) then
  else if (max_length < 0.0_DP) then
     if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
          write(stdout,'(a/a,f11.5,a)') &
          'WARNING in services_cubic_fit_minimum:','  cubic step (', &
          !min_length,') less than zero - setting to safe value'
          max_length,') less than zero - setting to safe value'
     !min_length = 0.10_DP
     max_length = 0.10_DP
     success = .false.
  else
     success = .true.
  endif

  end subroutine services_cubic_fit_maximum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gibo-end



  subroutine services_flush(unit)

    !==================================================================!
    ! Flushes output (stdout by default)                               !
    !------------------------------------------------------------------!
    ! Written by Peter Haynes 12 November 2004                         !
    ! Simplified with Fortran2003 intrinsic by JM Escartin 5 June 2015 !
    !==================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout

    implicit none

    ! Argument
    integer, optional, intent(in) :: unit

    ! Local variables
    integer :: unit_local     ! Local copy of argument
    integer :: ierr           ! Error flag

    if (present(unit)) then
       unit_local = unit
    else
       unit_local = stdout
    end if

    if (pub_on_root) then

       flush(unit_local, iostat=ierr)
       if (ierr /= 0) write(stdout,'(a,i6)') 'WARNING in services_flush: &
            &flush failed with code ',ierr

    end if

  end subroutine services_flush


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_radial_transform(l,power,nrpts,rlog,rab,nqpts,qmax, &
       rfunc,qfunc)

    !=========================================================================!
    ! Transforms a radial real space function on a logarithmic grid to a      !
    ! regular reciprocal space grid                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  l      (input)  : angular momentum for bessel functions                !
    !  power  (input)  : power of r to multiply integrand by                  !
    !  nrpts  (input)  : number of real space points                          !
    !  rlog   (input)  : the logarithmic grid                                 !
    !  rab    (input)  : the "spacings" of the logarithmic grid               !
    !  nqpts  (input)  : number of points in the regular reciprocal grid      !
    !  qmax   (input)  : the maximum q-vector for the transform               !
    !  rfunc  (input)  : the real space radial function                       !
    !  qfunc  (output) : the reciprocal space function                        !
    !-------------------------------------------------------------------------!
    ! Adapted for ONETEP by Nicholas Hine from a similar routine in CASTEP,   !
    ! which was originally written by Chris J. Pickard on 08/06/2000.         !
    !=========================================================================!

    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer,intent(in) :: nrpts
    integer,intent(in) :: l
    integer,intent(in) :: power
    real(kind=dp),intent(in) :: rlog(nrpts)
    real(kind=dp),intent(in) :: rab(nrpts)
    integer,intent(in) :: nqpts
    real(kind=dp),intent(in) :: qmax
    real(kind=dp),intent(in) :: rfunc(nrpts)
    real(kind=dp),intent(out) :: qfunc(nqpts)

    ! Local Variables
    integer :: ierr        ! Error flag
    integer :: nq,nr       ! Real and recip point counters
    real(kind=dp) :: q,qr        ! q-vector and product q.r
    real(kind=dp),allocatable :: rwgt(:)     ! Integration weights
    real(kind=dp),allocatable :: lookup(:,:) ! Bessel function lookup table

    ! BLAS subroutine
    external :: dgemv

    ! Allocate temporary variables
    allocate(rwgt(nrpts),stat=ierr)
    call utils_alloc_check('services_radial_transform','rwgt',ierr)
    allocate(lookup(nqpts,nrpts),stat=ierr)
    call utils_alloc_check('services_radial_transform','lookup',ierr)

    ! Check arguments
    if (modulo(nrpts,2)==0) call utils_abort('Error in &
         &services_radial_transform: number of grid points must be odd')

    ! Integration weights
    rwgt(1) = rab(1)*4.0_dp/12.0_dp*rlog(1)**power
    do nr=2,nrpts-1,2
       rwgt(nr)=rab(nr)*16.0_dp/12.0_dp*rlog(nr)**power
    end do
    do nr=3,nrpts-2,2
       rwgt(nr)=rab(nr)*8.0_dp/12.0_dp*rlog(nr)**power
    end do
    rwgt(nrpts) = rab(nr)*4.0_dp/12.0_dp*rlog(nrpts)**power

    ! Loop over the radial reciprocal space grid
    do nq=1,nqpts
       q = real(nq-1,dp)*qmax/real(nqpts-1,dp)
       lookup(nq,1) = services_sbessj(l,0.0_dp)*rwgt(1)
       do nr=2,nrpts
          qr = q*rlog(nr)
          lookup(nq,nr) = services_sbessj(l,qr)*rwgt(nr)
       end do
    end do

    ! Do the matrix-vector multiplication
    call dgemv('N',nqpts,nrpts,1.0_dp,lookup(1,1),nqpts,rfunc,1,0.0_dp, &
         qfunc,1)

    ! Deallocate temporary variables
    deallocate(lookup,stat=ierr)
    call utils_dealloc_check('services_radial_transform','lookup',ierr)
    deallocate(rwgt,stat=ierr)
    call utils_dealloc_check('services_radial_transform','rwgt',ierr)

  end subroutine services_radial_transform


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_regular_transform(l,power,nrpts,rmax,nqpts,qmax, &
       rfunc,qfunc)

    !=========================================================================!
    ! Transforms a radial real space function on a logarithmic grid to a      !
    ! regular reciprocal space grid                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  l      (input)  : angular momentum for bessel functions                !
    !  power  (input)  : power of r to multiply integrand by                  !
    !  nrpts  (input)  : number of real space points                          !
    !  rmax   (input)  : the maximum r-point of the regular real grid         !
    !  nqpts  (input)  : number of points in the regular reciprocal grid      !
    !  qmax   (input)  : the maximum q-vector for the transform               !
    !  rfunc  (input)  : the real space radial function                       !
    !  qfunc  (output) : the reciprocal space function                        !
    !-------------------------------------------------------------------------!
    ! Adapted for ONETEP by Nicholas Hine from a similar routine in CASTEP,   !
    ! which was originally written by Chris J. Pickard on 08/06/2000.         !
    !=========================================================================!

    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer,intent(in) :: l
    integer,intent(in) :: power
    integer,intent(in) :: nrpts
    real(kind=dp),intent(in) :: rmax
    integer,intent(in) :: nqpts
    real(kind=dp),intent(in) :: qmax
    real(kind=dp),intent(in) :: rfunc(nrpts)
    real(kind=dp),intent(out) :: qfunc(nqpts)

    ! Local Variables
    integer :: iq,ir         ! Real and recip point counters
    real(kind=dp) :: q,r     ! q-vector and r-point
    real(kind=dp) :: dr
    real(kind=dp), allocatable :: work(:)
    integer       :: ierr

    allocate(work(nrpts),stat=ierr)
    call utils_alloc_check('services_radial_transform','work',ierr)
    dr = rmax/real(nrpts-1,dp)

    ! Check arguments
    if (modulo(nrpts,2)==0) call utils_abort('Error in &
         &services_radial_transform: number of grid points must be odd')

    ! Loop over the radial reciprocal space grid
    do iq=1,nqpts
       q = real(iq-1,dp)*qmax/real(nqpts-1,dp)
       do ir=1,nrpts
          r = real(ir-1,dp)*rmax/real(nrpts-1,dp)
          work(ir) = rfunc(ir)*r**power*services_sbessj(l,q*r)
       end do
       qfunc(iq) = work(1) + work(nrpts)
       do ir=2,nrpts-1,2
          qfunc(iq) = qfunc(iq) + 4.0_DP*work(ir) + 2.0_DP*work(ir+1)
       end do
       qfunc(iq) = qfunc(iq) * dr / 3.0_DP
    end do

    deallocate(work,stat=ierr)
    call utils_dealloc_check('services_radial_transform','work',ierr)

  end subroutine services_regular_transform


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function services_radial_integral(npts,rab,f,inter)

    !=========================================================================!
    ! Simpson's rule integrator for a function on a radial mesh.              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  npts (input)   : number of points in radial mesh                       !
    !  rab (input)    : "spacings" of the radial grid                         !
    !  func (input)   : function to be integrated                             !
    !  inter (output) : intermediate results (optional)                       !
    !-------------------------------------------------------------------------!
    ! Adapted for ONETEP by Nicholas Hine from a similar routine in CASTEP,   !
    ! which was originally written by Chris J. Pickard on 09/06/2000.         !
    !=========================================================================!

    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP), intent(in) :: rab(:)
    real(kind=DP), intent(in) :: f(:)
    real(kind=DP), optional, intent(out) :: inter(:)

    ! Local Variables
    real(kind=DP), parameter :: c1o12 = 1.0_DP/12.0_DP
    real(kind=DP), parameter :: c5o12 = 5.0_DP/12.0_DP
    real(kind=DP), parameter :: c8o12 = 8.0_DP/12.0_DP
    real(kind=DP) :: fm1,f0,fp1
    integer :: ipt

    ! Check arguments
    if (modulo(npts,2)==0) call utils_abort('Error in &
         &services_radial_integral: number of grid points must be odd')

    ! Initialise
    services_radial_integral = 0.0_DP
    if (present(inter)) inter(1) = 0.0_DP

    ! Loop over alternate grid points, adding up integral and calculating
    ! intermediate values if required
    do ipt=2,npts-1,2

       fm1 = f(ipt-1)*rab(ipt-1)
       f0  = f(ipt  )*rab(ipt  )
       fp1 = f(ipt+1)*rab(ipt+1)

       services_radial_integral = services_radial_integral &
            + c5o12*fm1 + c8o12*f0 - c1o12*fp1
       if (present(inter)) inter(ipt) = services_radial_integral

       services_radial_integral = services_radial_integral &
            - c1o12*fm1 + c8o12*f0 + c5o12*fp1
       if (present(inter)) inter(ipt+1) = services_radial_integral

    end do

  end function services_radial_integral


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function services_radial_integral_rmax(npts,rab,r,rmax,f,inter)

    !=========================================================================!
    ! Wrapper for services_radial_integral to integrate up to a fixed rmax    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  npts (input)   : number of points in radial grid                       !
    !  rab (input)    : "spacings" of the radial grid                         !
    !  r   (input)    : radial grid position                                  !
    !  rmax (input)   : upper limit of integral                               !
    !  func (input)   : function to be integrated                             !
    !  inter (output) : array to hold intermediate results                    !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 02/06/2010.                                 !
    !=========================================================================!

    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP), intent(in) :: rab(:)
    real(kind=DP), intent(in) :: r(:)
    real(kind=DP), intent(in) :: rmax
    real(kind=DP), intent(in) :: f(:)
    real(kind=DP), intent(out) :: inter(:)

    ! Local Variables
    integer :: ir
    integer :: npts_int
    real(kind=DP) :: int_full

    ! Calculate the interpolation point
    ir = services_locate_interp(rmax,r,npts)
    ! Find endpoint: min of a) position of rmax + 4 pts, or b) end of array
    npts_int = min(npts,ir+4)
    npts_int = npts_int + modulo(npts_int,2) - 1
    ! Calculate the integral up to npts_int
    int_full = services_radial_integral(npts_int,rab,f,inter)
    ! Do the interpolation
    services_radial_integral_rmax = services_linear_interpolation(rmax, &
         inter(ir),inter(ir+1),r(ir),r(ir+1))

  end function services_radial_integral_rmax


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function services_regular_integral(npts,dr,f)

    !=========================================================================!
    ! Integrates a function on a regular grid.                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  npts (input)   : number of points in grid                              !
    !  dr (input)    :  spacing of the grid                                   !
    !  r   (input)    : radial grid position                                  !
    !  func (input)   : function to be integrated                             !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 09/09/2010.                                 !
    !=========================================================================!

    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP), intent(in) :: dr
    real(kind=DP), intent(in) :: f(:)

    ! Local Variables
    real(kind=DP) :: a(4)
    real(kind=DP) :: intf

    a(1)=17.0_DP/48.0_DP
    a(2)=59.0_DP/48.0_DP
    a(3)=43.0_DP/48.0_DP
    a(4)=49.0_DP/48.0_DP

    intf = a(1)*(f(1)+f(npts)) + &
         a(2)*(f(2)+f(npts-1)) + &
         a(3)*(f(3)+f(npts-2)) + &
         a(4)*(f(4)+f(npts-3))
    intf = intf + sum(f(5:npts-4))
    services_regular_integral = intf*dr

  end function services_regular_integral


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_radial_derivative(grad,func,npts,xmax)

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP),intent(in) :: xmax
    real(kind=DP),intent(in) :: func(:)
    real(kind=DP),intent(out) :: grad(:)

    ! Local Variables
    integer :: ipt
    real(kind=DP) :: t1,t2,t3

    ! Gradient is zero at origin by symmetry
    grad(1) = 0.0_DP
    t1 = 0.0_DP; t2 = 0.0_DP; t3 = 0.0_DP

    ! Cubic interpolation
    ! JCW: t0, t1, t2, t3 are the polynomial coefficients when the equation
    ! JCW:    fn(x) = t0 + t1*x + t2*x^2 + t3*x^3
    ! JCW: is fitted to 4 equally spaced points, x = -1, 0, 1, 2, i.e.
    ! JCW:    fn(-1) = f1 = t0 - t1 + t2 - t3
    ! JCW:    fn(0)  = f2 = t0
    ! JCW:    fn(1)  = f3 = t0 + t1 + t2 + t3
    ! JCW:    fn(2)  = f4 = t0 + 2*t1 + 4*t2 + 8*t3
    ! JCW: Using Sage Math to solve this set of linear equations, i.e.
    ! JCW:   solve( [ f1 == fn(-1), f2 == fn(0), f3 == fn(1), f4 == fn(2) ], t0, t1, t2, t3 )
    ! JCW: we obtain
    ! JCW: [[
    ! JCW:   t0 == f2,
    ! JCW:   t1 == -1/3*f1 - 1/2*f2 + f3 - 1/6*f4,
    ! JCW:   t2 == 1/2*f1 - f2 + 1/2*f3,
    ! JCW:   t3 == -1/6*f1 + 1/2*f2 - 1/2*f3 + 1/6*f4
    ! JCW: ]]
    ! JCW: Since we only want the derivative at x=0, we have
    ! JCW:   fn'(x)  = t1 + 2*t2*x + 3*t3*x^2
    ! JCW:   fn'(0)  = t1
    do ipt=2,npts-2
       t1 = (6.0_DP*func(ipt+1) - 2.0_DP*func(ipt-1) - 3.0_DP*func(ipt) &
            - 1.0_DP*func(ipt+2)) / 6.0_DP
       t2 = (1.0_DP*func(ipt-1) + 1.0_DP*func(ipt+1) - 2.0_DP*func(ipt)) &
            / 2.0_DP
       t3 = (1.0_DP*func(ipt+2) - 1.0_DP*func(ipt-1) + 3.0_DP*func(ipt) &
            - 3.0_DP*func(ipt+1)) / 6.0_DP
       grad(ipt) = t1
    end do

    ! Last two points
    ! JCW: Using the fit where f0 is f(npts-2)
    ! JCW:   fn'(1) = t1 + 2*t2   + 3*t3
    ! JCW:   fn'(2) = t1 + 2*t2*2 + 3*t3*2^2
    grad(npts-1) = t1 + 2.0_DP*t2 + 3.0_DP*t3
    grad(npts) = t1 + 4.0_DP*t2 + 12.0_DP*t3

    ! Normalise for dr
    ! JCW: We have taken the derivative with respect to x_n where x_n
    ! JCW: is related to the variable we want the derivative wrt
    ! JCW:    x = x_n / (npts-1) * xmax
    ! JCW: so use the chain rule
    ! JCW:   dfn/dx_n * dx_n/dx = dfn/dx
    ! JCW: where
    ! JCW    dx_n/dx = (npts-1)/xmax
    grad(:) = grad(:)*real(npts-1,kind=DP)/xmax

  end subroutine services_radial_derivative


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_analytic_limit(npts,rad,func)

    !=========================================================================!
    ! Finds the r->0 limit of a radial function through a spline fit to the   !
    ! first few points on the radial grid.                                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  npts (input)   : number of points in radial grid                       !
    !  rad (input)    : radial grid position                                  !
    !  func (input)   : function to be integrated                             !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 20/07/2010.                                 !
    !=========================================================================!

    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP), intent(in) :: rad(:)
    real(kind=DP), intent(inout) :: func(:)

    if ((size(rad)<npts).or.(size(func)<npts)) then
       call utils_abort('Error in services_analytic_limit: mismatching array &
            & sizes')
    end if

    if (npts>3) then
       func(1) = func(4)+3*(func(2)-func(3))
    else
       func(1) = func(2)
    end if

  end subroutine services_analytic_limit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine services_rationalise_coords(nr,r_abs,cell)

    !=============================================================!
    ! This subroutine rationalises cartesian coordinates by       !
    ! translating them back into the home simulation cell.        !
    !-------------------------------------------------------------!
    ! Originally, md_update_coordinates subroutine written by     !
    ! Arash A. Mostofi in 2006. Modified and displaced by         !
    ! Simon M.-M. Dubois in Oct. 2010.                            !
    !=============================================================!

    use constants,       only: two_pi
    use ion,             only: element
    use simulation_cell, only: CELL_INFO
    use utils,           only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in) :: cell
    integer, intent(in) :: nr
    real(kind=dp), intent(inout) :: r_abs(3,nr)

    ! Local variables
    integer :: j,iat
    integer :: rshift
    real(kind=dp) :: r_frac(3,nr)

    ! convert Cartesian coordinates to fractional and rationalise
    do iat=1,nr
       r_frac(1,iat) = cell%b1%x*r_abs(1,iat) + cell%b1%y*r_abs(2,iat) + cell%b1%z*r_abs(3,iat)
       r_frac(2,iat) = cell%b2%x*r_abs(1,iat) + cell%b2%y*r_abs(2,iat) + cell%b2%z*r_abs(3,iat)
       r_frac(3,iat) = cell%b3%x*r_abs(1,iat) + cell%b3%y*r_abs(2,iat) + cell%b3%z*r_abs(3,iat)
       r_frac(:,iat) = r_frac(:,iat) / two_pi
       do j=1,3
          if (r_frac(j,iat).lt.0.0_dp .or. r_frac(j,iat).ge.1.0_dp) then
             rshift = floor(r_frac(j,iat))
             r_frac(j,iat) = r_frac(j,iat) - real(rshift)
          endif
       enddo
    enddo

    ! convert rationalised fractional coordinates back to Cartesian
    do iat=1,nr
       r_abs(1,iat) = cell%a1%x*r_frac(1,iat) + cell%a2%x*r_frac(2,iat) + cell%a3%x*r_frac(3,iat)
       r_abs(2,iat) = cell%a1%y*r_frac(1,iat) + cell%a2%y*r_frac(2,iat) + cell%a3%y*r_frac(3,iat)
       r_abs(3,iat) = cell%a1%z*r_frac(1,iat) + cell%a2%z*r_frac(2,iat) + cell%a3%z*r_frac(3,iat)
    enddo

    return

  end subroutine services_rationalise_coords


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function services_maxboltzdist()

    !==============================================================!
    ! This function generates a pseudo-random number normally      !
    ! distributed (zero expectation, unit variance)                !
    !--------------------------------------------------------------!
    ! Originally written by Arash A. Mostofi as md_maxboltzdist    !
    ! Displaced by Simon M.-M. Dubois in July 2011.                !
    !==============================================================!

    implicit none

    real(kind=dp)  :: services_maxboltzdist
    real(kind=dp) :: ran(2),x(2),w,f

    do
       call random_number(ran)
       x(:) = 2.0_dp * ran(:) - 1.0_dp
       w = x(1)*x(1) + x(2)*x(2)
       if (w.lt.1.0_dp) exit
    enddo

    f = sqrt((-2.0_dp*log(w))/w)
    services_maxboltzdist = x(1) * f

    return

  end function services_maxboltzdist

  subroutine services_maxboltzdist_parallel(maxboltzdist,nat)

    ! Modified by V. Vitale to use Zuehlsdorff random number       !
    ! generator, in order to avoid different numbers for different !
    ! compilers

    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast
    use rundat, only: pub_rand_seed
    use utils, only: utils_rand_gen, utils_abort, utils_assert, &
         utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=dp), intent(inout), allocatable :: maxboltzdist(:)
    integer, intent(in) :: nat

    ! Local variable
    real(kind=dp), allocatable :: ran(:)
    real(kind=dp) :: x(2),w,f
    integer :: numran,ir, ierr, ir2

    ir = 1
    ir2 = 1
    numran = 30*nat + 50
    allocate(ran(numran),stat=ierr)
    call utils_alloc_check('services_maxboltzdist_parallel','ran',ierr)
    ! If we want the same random numbers independent of the run and the
    ! compiler, then call the congruental random generator on the root proc
    ! and broadcast the result
    if (pub_on_root) then
       call utils_rand_gen(ran,numran,pub_rand_seed)
       do
          x(:) = 2.0_dp * ran(ir:ir+1) - 1.0_dp
          w    = x(1)*x(1) + x(2)*x(2)
          if ( w .lt. 1.0_dp) then
             maxboltzdist(ir2) = x(1) * sqrt((-2.0_dp*log(w)/w))
             ir2 = ir2 + 1
             if (ir2 > 3*nat) exit
          end if
          ir = ir + 2
          call utils_assert(ir .lt. numran, 'Not enough random numbers have &
               & been found, try with a different seed.')
       end do
    end if
    call comms_bcast(pub_root_proc_id, maxboltzdist)

    deallocate(ran,stat=ierr)
    call utils_dealloc_check('services_maxboltzdist_parallel','ran',ierr)

    return

  end subroutine services_maxboltzdist_parallel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function services_gammadist(gamma_alpha)

    !=============================================================!
    ! Generate a pseudo random number out of a Gamma              !
    ! distribution Pg(alpha) :                                    !
    !   Pg(alpha) = x**(alpha-1) * exp(-x) / Gamma(alpha)         !
    ! with alpha integer                                          !
    !                                                             !
    ! For details, see :                                          !
    !   G. Marsaglia, W.W. Tsang,                                 !
    !   ACM Transactions on Mathematical Software,                !
    !   Vol. 26, 363-372 (2000)                                   !
    !-------------------------------------------------------------!
    ! Originally written by Gilberto Teobaldi (May 2011),         !
    ! moved here by Simon M.-M. Dubois in June 2013.              !
    !=============================================================!

      implicit none

      integer, intent(in) :: gamma_alpha

      real(kind=dp) :: services_gammadist
      real(kind=dp) :: ran
      real(kind=dp) :: d, c, x, v, tmp

      d = gamma_alpha - 1.0_DP/3.0_DP
      c = 1.0_DP/SQRT(9.0_DP*d)

      GAMMA_LOOP: do

         v = -1.0_DP
         do while (v < 0.0_DP)
            x = services_maxboltzdist()
            v = (1.0_DP + c*x)
            v = v**3
         enddo

         call random_number(ran)

         if (ran < (1.0_DP - 0.0331_DP * x**4) ) then
            services_gammadist = d*v
            EXIT GAMMA_LOOP
         endif

         tmp = 0.5_DP * x**2 + d * (1.0_DP - v + log(v))

         if ( log(ran) < tmp ) then
            services_gammadist = d*v
         endif

      enddo GAMMA_LOOP

      return

    end function services_gammadist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine services_rms_fit(alpha,nr,ncol,nrow,r_in,r_out)

    !========================================================================!
    ! Computes the coefficients alpha such that the root mean square         !
    ! deviation between R_1 and R = \Sum_{i=2}^{nr} (alpha_i * Rin_i)        !
    ! is minimized.                                                          !
    !                                                                        !
    ! Arguments :                                                            !
    ! alpha     (output) : computed coefficients                             !
    ! nr        (input)  : the number of elements in the basis set           !
    ! nrow      (input)  : the first dimension of Rin and Rout               !
    ! ncol      (input)  : the second dimension of Rin and Rout              !
    ! Rin       (input)  : basis set [dim(nrow,ncol,nr)]                     !
    ! Rout      (input)  : reference vector [dim(nrow,ncol)]                 !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois 06/04/2011                               !
    !========================================================================!

    use constants, only: stdout
    use rundat, only: pub_debug_on_root
    use utils, only: utils_assert

    implicit none

    ! Argument
    integer, intent(in)          :: nr     ! The number of input elements
    integer, intent(in)          :: ncol   !
    integer, intent(in)          :: nrow                ! # of row in each dataset
    real(kind=DP), intent(in)    :: r_in(nrow,ncol,nr)  ! datasets
    real(kind=DP), intent(in)    :: r_out(nrow,ncol)    !
    real(kind=DP), intent(out)   :: alpha(nr)

    ! Internal variables
    real(kind=DP)   :: overlap(nr,nr)
    real(kind=DP)   :: ivec(nr)
    integer  :: ii, jj
    integer  :: icol, irow
    integer  :: ipiv(nr)
    integer  :: ierr

    ! LAPACK subroutine
    external :: dgesv

    !======================================================================!

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Entering services_rms_fit'

    ! Compose the overlap matrix
    overlap(:,:) = 0.0_dp
    do jj = 1, nr
       do ii = 1, nr
          do irow = 1, nrow
             do icol = 1, ncol
                overlap(ii,jj) = overlap(ii,jj) + &
                   r_in(irow,icol,ii)*r_in(irow,icol,jj)
             enddo
          enddo
       enddo
    enddo

    ! Compose the indep. vector
    ivec(:) = 0.0_dp
    do ii = 1, nr
       do irow = 1, nrow
          do icol = 1, ncol
             ivec(ii) = ivec(ii) + &
                r_in(irow,icol,ii)*r_out(irow,icol)
          enddo
       enddo
    enddo

    ! Solve the system of linear equations
    call dgesv(nr,1,overlap,nr,ipiv,ivec,nr,ierr)
    call utils_assert(ierr == 0, 'Error in services_rms_fit(): &
           &computation of the extrapolation coefficients with &
           &lapack_dgesv failed with error, ', ierr)

    if (ierr .ne. 0) then
       alpha(:) = 0.0_dp
       alpha(1) = 1.0_dp
    else
       alpha(:) = ivec(:)
    endif

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Leave services_rms_fit'

    return

  end subroutine services_rms_fit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !========================================================================!
    ! Simple facility for calculating the incomplete gamma function          !
    ! \Gamma(a,x) for non-negative x.                                        !
    !                                                                        !
    ! The canonical way to do this is to use a number of different function  !
    ! fits on a number of subdomains, this is how it's done within netlib.   !
    ! Porting this from netlib is a pain, because there is an endless chain  !
    ! of dependencies. Instead we follow a very simple approach of evalua-   !
    ! ting \Gamma from the definition and storing it in a look-up table.     !
    ! Queries to the look-up are then combined with linear interpolation to  !
    ! get the values.                                                        !
    !                                                                        !
    ! The lookup table is allocated and returned by prepare_lookup() and it's!
    ! the responsibility of the caller to store it for as long as they plan  !
    ! to query it with from_lookup(), then to call destroy_lookup() on it.   !
    !                                                                        !
    ! Note that the lookup only stores values of \Gamma for a single value   !
    ! of 'a' -- so this approach only makes sense for a fixed 'a', or if     !
    ! you only have a few values of 'a', and many values of 'x'. This is the !
    ! usual scenario. Each new 'a' requires a new lookup table.              !
    !                                                                        !
    !========================================================================!

  subroutine services_incomplete_gamma_prepare_lookup(gamma_lookup, a, npts_x, &
       max_x)
    !========================================================================!
    ! Prepares a lookup table for the incomplete gamma function for a given  !
    ! parameter a.                                                           !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    !   gamma_lookup (output): The lookup is returned here. This array will  !
    !      be allocated, filled and returned. The caller needs to keep it,   !
    !      then to call destroy_lookup() on it when no longer needed.        !
    !   a (input): The parameter 'a' for which the lookup is to be calc'd.   !
    !   npts_x (input, opt): The number of points in the lookup.             !
    !      10000001 is a good value and will be used if not specified.       !
    !   max_x (input, opt): Maximum value of the argument x for which the    !
    !       lookup is filled. All queries for larger x later on will return  !
    !       zero. The value 30.0 will be used if not specified. This is OK   !
    !       for small a (say, below 5). Larger a's might need larger cutoffs.!
    !------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2015.                                !
    ! Improved to remove module-wide state by Jacek Dziedzic in October 2015.!
    !========================================================================!

    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out), allocatable :: gamma_lookup(:)
    real(kind=DP), intent(in)               :: a
    integer, intent(in), optional           :: npts_x
    real(kind=DP), intent(in), optional     :: max_x

    ! jd: Local variables
    real(kind=DP), allocatable :: integrand(:)
    real(kind=DP) :: x, hx
    real(kind=DP) :: local_max_x
    integer :: xi
    integer :: local_npts_x
    integer :: ierr
    integer, parameter          :: default_npts_x = 10000001
    real(kind=DP), parameter    :: default_max_x = 30.0_DP
    character(len=*), parameter :: myself = &
         'services_incomplete_gamma_prepare_lookup'

    ! -------------------------------------------------------------------------

    if(present(npts_x)) then
       local_npts_x = npts_x
    else
       local_npts_x = default_npts_x
    end if

    if(present(max_x)) then
       local_max_x = max_x
    else
       local_max_x = default_max_x
    end if

    call utils_assert(local_npts_x >= 1, myself//': Illegal npts_x.')
    allocate(gamma_lookup(local_npts_x),stat=ierr)
    call utils_alloc_check(myself,'services_incomplete_gamma_lookup',ierr)

    allocate(integrand(local_npts_x),stat=ierr)
    call utils_alloc_check(myself,'integrand',ierr)

    hx = local_max_x/(local_npts_x-1)

    ! jd: Fill in the lookup table.
    !     This goes from xmax down to zero, accumulating the integral from
    !     the right.
    do xi = local_npts_x, 1, -1
       x = (xi-1) * hx
       if(xi > 1) then
          integrand(xi) = x**(a-1.0_DP) * exp(-x)
       else
          ! jd: Some values of 'a' have a singularity in the integrand at 0.
          !     We don't care much about this, since we won't be allowing
          !     lookups below hx anyway. Just set it to the value of next point.
          integrand(xi) = integrand(xi+1)
       end if
       if(xi == local_npts_x) then
          gamma_lookup(xi) = 0D0
       else
          ! jd: integral(x) ~= integral(x+hx) + dx * df
          gamma_lookup(xi) = gamma_lookup(xi+1) + &
               0.5_DP * hx * (integrand(xi) + integrand(xi+1))
       end if
    end do

    deallocate(integrand,stat=ierr)
    call utils_dealloc_check(myself,'integrand',ierr)

  end subroutine services_incomplete_gamma_prepare_lookup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine services_incomplete_gamma_destroy_lookup(gamma_lookup)
    !========================================================================!
    ! Destroys the lookup created by prepare_lookup().                       !
    !------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2015.                                !
    ! Improved to remove module-wide state by Jacek Dziedzic in October 2015.!
    !========================================================================!

    use utils, only: utils_dealloc_check, utils_abort

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(inout), allocatable :: gamma_lookup(:)

    ! jd: Local variables
    integer                     :: ierr
    character(len=*), parameter :: myself = &
         'services_incomplete_gamma_destroy_lookup'

    ! -------------------------------------------------------------------------

    if(allocated(gamma_lookup)) then
       deallocate(gamma_lookup,stat=ierr)
       call utils_dealloc_check(myself,'services_incomplete_gamma_lookup',ierr)
    else
       call utils_abort(myself//': Attempt to destroy gamma_lookup &
            &more than once.')
    end if

  end subroutine services_incomplete_gamma_destroy_lookup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function services_incomplete_gamma_from_lookup(a, x, &
       gamma_lookup, max_x)
    !========================================================================!
    ! Serves a value of \Gamma(a,x) from the lookup, using linear interpo-   !
    ! lation for x.                                                          !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    !   a, x (in): The arguments of the function.                            !
    !   gamma_lookup (in): A lookup table prepared in advance using          !
    !                      prepare_lookup().                                 !
    !   max_x (in, opt): The same value that was given to prepare_lookup().  !
    !                    See prepare_lookup() for info on default value.     !
    !------------------------------------------------------------------------!
    ! Notes:                                                                 !
    !   The lookup must have been prepared for the same 'a'. This is checked.!
    !------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2015.                                !
    ! Improved to remove module-wide state by Jacek Dziedzic in October 2015.!
    !========================================================================!

#ifdef DONT_HAVE_GAMMA
    use constants, only: TWO_THIRDS
    use utils, only: utils_feature_not_supported
#endif
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in)              :: a
    real(kind=DP), intent(in)              :: x
    real(kind=DP), intent(in), allocatable :: gamma_lookup(:)
    real(kind=DP), intent(in), optional    :: max_x

    ! jd: Local variables
    integer :: npts_x
    real(kind=DP) :: hx
    real(kind=DP) :: local_max_x
    real(kind=DP), parameter    :: eps_tol = 1D-12
    real(kind=DP), parameter    :: default_max_x = 30.0_DP
    character(len=*), parameter :: myself = &
         'services_incomplete_gamma_from_lookup'

    ! -------------------------------------------------------------------------

    call utils_assert(allocated(gamma_lookup), &
         myself//': The lookup needs to be prepared first.')
    call utils_assert(x >= 0D0, myself//': x must be >=0',x)

    if(present(max_x)) then
       local_max_x = max_x
    else
       local_max_x = default_max_x
    end if

    ! jd: For large values of argument just return 0.0
    if(x >= local_max_x) then
       services_incomplete_gamma_from_lookup = 0.0_DP
    else
       npts_x = size(gamma_lookup)
       hx = local_max_x/(npts_x-1)
       ! jd: For extremely small values our interpolation can go wrong, as
       !     the integrand can have a singularity at zero. Just return
       !     \Gamma(a,0) there. This assumes hx is decently small.
       !     \Gamma(a,0) = \Gamma(a) is an intrinsic, which is technically
       !     F2008, but is widely implemented.
       if(x <= hx) then
#ifdef DONT_HAVE_GAMMA
          ! jd: Compiler doesn't define the intrinsic. If the argument is 2/3,
          !     return a hard-coded value. If not, give up.
          if(abs(a - TWO_THIRDS) < eps_tol) then
             services_incomplete_gamma_from_lookup = &
                  1.3541179394264004169452880281_DP
          else
             call utils_feature_not_supported('The F2008 intrinsic Gamma')
          end if
#else
          ! jd: Compiler defines the intrinsic, good.
          services_incomplete_gamma_from_lookup = &
               gamma(a)
          ! ***************************************************************
          ! jd: If your compiler complains here that it doesn't know what
          !     gamma is, add '-DDONT_HAVE_GAMMA' to your compiler flags.
          ! ***************************************************************
#endif
       else
          ! jd: Otherwise return value from interpolation using the lookup table
          services_incomplete_gamma_from_lookup = services_1d_interpolation( &
               gamma_lookup, npts_x, x/local_max_x*real(npts_x-1,kind=DP),0)
       end if
    end if

  end function services_incomplete_gamma_from_lookup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine services_gradient_on_grid(func_grad,func,grid)

    !========================================================================!
    ! This subroutine calculates the gradient of a function on a grid by     !
    ! transforming the function to reciprocal space and applying the         !
    ! derivative operator along each Cartesian direction. The gradient in    !
    ! real space is obtained by transforming the result back to real space.  !
    !                                                                        !
    ! If we define the forward Fourier transform as                          !
    !   f(r) = 1/V \sum_{G} f(G) \exp(iG.r)                                  !
    ! then                                                                   !
    !   \nabla f(r) = 1/V \sum_{G} f(G) \nabla \exp(iG.r)                    !
    ! which gives                                                            !
    !   \nabla f(r) = 1/V \sum_{G} ( iG f(G) ) \exp(iG.r)                    !
    !------------------------------------------------------------------------!
    ! Arguments                                                              !
    !   func_fine       (in) : Input function on the grid                    !
    !   density_grad (inout) : Gradient of the function on the grid          !
    !   grid            (in) : Grid definition                               !
    !------------------------------------------------------------------------!
    ! Written by James C. Womack, 05/2017                                    !
    ! Derived from xc_gradients routine written by Peter Haynes in February  !
    ! 2004 and later modified by Nicholas Hine in October 2010.              !
    !========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_proc_id
    use constants, only: cmplx_i
    use fourier, only: fourier_apply_cell_forward, fourier_apply_cell_backward
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO),intent(in)   :: grid
    real(kind=DP), intent(inout) :: func_grad(grid%ld1,grid%ld2,grid%max_slabs12,3)
    real(kind=DP), intent(in)    :: func(grid%ld1,grid%ld2,grid%max_slabs12)

    ! Parameters
    character(len=*), parameter :: myself = "services_gradient_on_grid"

    ! Local Variables
    complex(kind=DP), allocatable :: recip_work(:,:,:,:)
    integer :: i3,i2,islab23           ! Loop counters over grid points
    integer :: k                       ! Loop counter for Cartesian components
    complex(kind=DP) :: f_gtmp         ! Temporary copy of reciprocal function at G
    real(kind=DP) :: gvec(3)
    integer :: ierr


    ! Allocate reciprocal space work array
    allocate(recip_work(grid%ld3,grid%ld2,grid%max_slabs23,3),stat=ierr)
    call utils_alloc_check(myself,'recip_work',ierr)

    ! Fourier transform to reciprocal space
    call fourier_apply_cell_forward(func(:,:,:),recip_work(:,:,:,1),grid)

    ! TODO This set of loops is a good candidate for OpenMP parallelization
    ! Apply grad in reciprocal space (multiply by iG)
    do islab23=1,grid%num_slabs23
       do i2=1,grid%n2
          do i3=1,grid%n3

             ! Apply grad
             f_gtmp = recip_work(i3,i2,islab23,1)

             call cell_grid_recip_pt(gvec,islab23 + &
                  grid%first_slab23(pub_my_proc_id) - 1,i2,i3,grid)

             do k=1,3
                recip_work(i3,i2,islab23,k) = f_gtmp * cmplx_i * gvec(k)
             end do

          end do
       end do
    end do

    ! Transform gradient back to real space
    do k=1,3
       call fourier_apply_cell_backward(func_grad(:,:,:,k),recip_work(:,:,:,k),grid)
    end do

    ! Deallocate reciprocal space work array
    deallocate(recip_work,stat=ierr)
    call utils_dealloc_check(myself,'recip_work',ierr)

  end subroutine services_gradient_on_grid

  elemental function services_factorial(n) result(f)
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! Compute the factorial of integer n, i.e.                                 !
    !   n! = \prod_{k=1}^{n} k                                                 !
    ! where 0! == 1 by definition.                                             !
    !                                                                          !
    ! This function only returns correct results for n >= 0. For speed, this   !
    ! is not checked, so care should be taken to provide valid input.          !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS AND RESULT                                                   !
    ! <name> <in/out/inout> <arg descrption>                                   !
    ! n              in     the integer n for which n! will be computed        !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        James C. Womack                                        !
    ! Date of creation: 10/2018                                                !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!
    implicit none

    ! Parameters
    !character(len=*), parameter :: myself = "services_factorial"

    ! Arguments
    integer, intent(in) :: n

    ! Result
    integer :: f

    ! Local variables
    integer :: k

    f = 1
    do k = 2, n
       f = f * k
    end do

  end function services_factorial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module services
