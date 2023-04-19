! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Nicholas D.M. Hine
!
!   based in many cases on existing versions by Chris-Kriton Skylaris,
!   Arash A. Mostofi, Nicholas D.M. Hine and Peter D. Haynes
!
!   Several extra routines have been added by Laura E. Ratcliff
!
!   Thomas Young Centre
!   Imperial College London
!   Exhibition Road
!   UK
!
!   Subsequent additions and modifications by Ed Tait, Gabriel Constantinescu,
!   Andrea Greco, and Jose M Escartin.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module projectors

  use datatypes, only: FUNCTIONS
  use constants, only: DP, PI, stdout
  use geometry, only: POINT

  implicit none

  private

  ! Derived type to store details of a set of projectors
  type PROJECTOR_SET

     ! ndmh: number of projector species
     integer :: n_proj_species

     ! ndmh: definition of radial reciprocal space grid for each species
     real(kind=DP), allocatable :: gmax(:)
     integer, allocatable :: n_rad_pts(:)

     ! ndmh: projectors on radial reciprocal space grid
     integer, allocatable :: num_shells(:)
     integer, allocatable :: ang_mom(:,:)
     real(kind=DP), allocatable :: rad_proj_recip(:,:,:)

     ! ndmh: centres of the projectors of each atom
     type(POINT), allocatable, dimension(:) :: proj_centre

     ! ndmh: max radius of all the projectors on each atom
     real(kind=DP), allocatable, dimension(:) :: proj_max_radius

     ! ndmh: species of each atom
     integer, allocatable, dimension(:) :: proj_species

     ! ndmh: number of projectors on atoms of each species
     integer, allocatable, dimension(:) :: species_num_proj

     ! ndmh: first entry for this species in fftbox_proj_recip
     integer, allocatable, dimension(:) :: species_first_proj

     ! ndmh: reciprocal space FFTbox containing each type of projector
     complex(kind=DP), allocatable, dimension(:,:,:,:) :: fftbox_proj_recip

     ! ndmh: real-space projectors in PPDs (optional)
     type(FUNCTIONS) :: projs_on_grid

     ! ndmh: whether to normalise the projectors of this set
     logical :: normalise

     ! gcc32 (O_prec)_{i,j}(s) (check precond_matrix_calculate for explanations)
     real(kind=DP), allocatable, dimension(:,:,:) :: O_prec

  end type PROJECTOR_SET

  public :: PROJECTOR_SET

  public :: projectors_allocate_set
  public :: projectors_deallocate_set
  public :: projectors_init_fftbox_recip
  public :: projectors_exit_fftbox_recip
  public :: projector_in_box_recip
  public :: projector_in_box_recip_gradg
  public :: projectors_create_real
  public :: projectors_func_ovlp_box
  public :: projectors_func_grad_ovlp_box
  public :: projectors_func_pos_ovlp_box
  public :: projectors_gradient_batch
  public :: projectors_grad_precond_batch
  public :: projectors_proj_overlap
  public :: projectors_commutator
  public :: projectors_func_pos_ovlp_box_gradg

  ! agrecocmplx: interface to switch between real
  ! and complex projectors
  interface projectors_proj_overlap
     module procedure projectors_proj_overlap_real
     module procedure projectors_proj_overlap_cmplx
  end interface

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_allocate_set(proj_set,max_shells,max_npts,par,is_cmplx)

    !==============================================================!
    ! This subroutine allocates storage for a set of projectors    !
    !--------------------------------------------------------------!
    ! Written by Nicholas Hine on 19/05/2010.                      !
    !==============================================================!

    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(PROJECTOR_SET), intent(inout) :: proj_set
    integer, intent(in) :: max_shells, max_npts
    ! agrecocmplx
    logical, optional, intent(in) :: is_cmplx
    type(PARAL_INFO), intent(in) :: par

    ! Local Variables
    integer :: ierr      ! error flag
    ! agrecocmplx
    logical :: loc_cmplx

    ! agrecocmplx
    if (present(is_cmplx)) then
       loc_cmplx = is_cmplx
    else
       loc_cmplx = .false.
    end if

    ! Allocate internal arrays in proj_set relating to projector species
    allocate(proj_set%species_num_proj(proj_set%n_proj_species),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%species_num_proj',ierr)
    allocate(proj_set%species_first_proj(proj_set%n_proj_species),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%species_first_proj',ierr)
    allocate(proj_set%gmax(proj_set%n_proj_species),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%gmax',ierr)
    allocate(proj_set%n_rad_pts(proj_set%n_proj_species),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%n_rad_pts',ierr)
    allocate(proj_set%num_shells(proj_set%n_proj_species),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%num_shells',ierr)
    allocate(proj_set%ang_mom(max_shells,proj_set%n_proj_species),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%ang_mom',ierr)
    allocate(proj_set%rad_proj_recip(max_npts,max_shells, &
         proj_set%n_proj_species),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%rad_proj_recip',ierr)

    ! Allocate arrays in proj_set relating to actual projectors
    allocate(proj_set%proj_centre(par%nat),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%proj_centre',ierr)
    allocate(proj_set%proj_max_radius(par%nat),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%proj_max_radius',ierr)
    allocate(proj_set%proj_species(par%nat),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%proj_species',ierr)

    ! Default is no normalisation
    proj_set%normalise = .false.

    ! agrecocmplx: set real/complex character of projs_on_grid
    proj_set%projs_on_grid%iscmplx = loc_cmplx

#ifdef POINTERFUN
    nullify(proj_set%projs_on_grid%d)
    nullify(proj_set%projs_on_grid%z)
#endif

  end subroutine projectors_allocate_set


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_deallocate_set(proj_set)

    !==============================================================!
    ! This subroutine deallocates storage for a set of projectors  !
    !--------------------------------------------------------------!
    ! Written by Nicholas Hine on 19/05/2010.                      !
    ! Modified by Andrea Greco on 10/05/2015 to use complex NGWFs. !
    !==============================================================!

    use datatypes, only: data_check_alloc, data_functions_dealloc
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(PROJECTOR_SET), intent(inout) :: proj_set

    ! Local Variables
    integer :: ierr      ! error flag

    ! Deallocate storage of evaluated projectors in real/recip space
    if (data_check_alloc(proj_set%projs_on_grid)) then
       call data_functions_dealloc(proj_set%projs_on_grid)
    end if
    if (allocated(proj_set%fftbox_proj_recip)) then
       deallocate(proj_set%fftbox_proj_recip,stat=ierr)
       call utils_dealloc_check('projectors_deallocate_set', &
            'proj_set%fftbox_proj_recip',ierr)
    end if

    ! Deallocate arrays in proj_set relating to actual projectors
    deallocate(proj_set%proj_species,stat=ierr)
    call utils_dealloc_check('projectors_deallocate_set', &
         'proj_set%proj_species',ierr)
    deallocate(proj_set%proj_max_radius,stat=ierr)
    call utils_dealloc_check('projectors_deallocate_set', &
         'proj_set%proj_max_radius',ierr)
    deallocate(proj_set%proj_centre,stat=ierr)
    call utils_dealloc_check('projectors_deallocate_set', &
         'proj_set%proj_centre',ierr)

    ! Deallocate internal arrays in proj_set relating to projector species
    deallocate(proj_set%rad_proj_recip,stat=ierr)
    call utils_dealloc_check('projectors_allocate_set', &
         'proj_set%rad_proj_recip',ierr)
    deallocate(proj_set%ang_mom,stat=ierr)
    call utils_dealloc_check('projectors_allocate_set', &
         'proj_set%ang_mom',ierr)
    deallocate(proj_set%num_shells,stat=ierr)
    call utils_dealloc_check('projectors_allocate_set', &
         'proj_set%num_shells',ierr)
    deallocate(proj_set%n_rad_pts,stat=ierr)
    call utils_dealloc_check('projectors_allocate_set', &
         'proj_set%n_rad_pts',ierr)
    deallocate(proj_set%gmax,stat=ierr)
    call utils_dealloc_check('projectors_allocate_set', &
         'proj_set%gmax',ierr)
    deallocate(proj_set%species_first_proj,stat=ierr)
    call utils_dealloc_check('projectors_deallocate_set', &
         'proj_set%species_first_proj',ierr)
    deallocate(proj_set%species_num_proj,stat=ierr)
    call utils_dealloc_check('projectors_deallocate_set', &
         'proj_set%species_num_proj',ierr)

  end subroutine projectors_deallocate_set


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_init_fftbox_recip(proj_set,cell,fftbox, &
       kpt,delta,cart,fine)

    !==============================================================!
    ! This subroutine initialises a set of projectors in FFTboxes  !
    !--------------------------------------------------------------!
    ! Written by Nicholas Hine on 12/12/2011.                      !
    !==============================================================!

    use fft_box, only: FFTBOX_INFO
!$  use rundat, only: pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check

    implicit none

    ! Arguments
    type(PROJECTOR_SET), intent(inout) :: proj_set
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(POINT), intent(in), optional :: kpt
    real(kind=dp), intent(in), optional :: delta
    integer, intent(in), optional :: cart
    character(len=*), optional, intent(in) :: fine


    ! Local Variables
    integer :: proj_count, proj_count_on_species
    integer :: iproj,tot_num_proj
    integer :: shell, am, azim
    integer :: isp
    integer :: ierr      ! error flag
    type(POINT) :: kpt_loc
    logical :: loc_grad
    logical :: found
    real(kind=DP) :: proj_norm

    call timer_clock('projectors_init_fftbox_recip',1)

    ! ndmh: check for optional arguments
    if (present(kpt)) then
       kpt_loc = kpt
    else
       kpt_loc%x = 0.0_DP ; kpt_loc%y = 0.0_DP ; kpt_loc%z = 0.0_DP
    end if
    if (present(delta).and.present(cart)) then
       if ((cart<1).or.(cart>3)) then
          call utils_abort('Error in projectors_init_fftbox_recip: Invalid &
               &Cartesian direction supplied')
       end if
       loc_grad = .true.
    else
       loc_grad = .false.
    end if

    ! Count total projectors in this set (including m)
    proj_count = 0
    do isp=1,proj_set%n_proj_species
       proj_count = proj_count + proj_set%species_num_proj(isp)
    end do

    ! Allocate projector FFTboxes
    if (present(fine)) then
       allocate(proj_set%fftbox_proj_recip(fftbox%total_ld1_dbl, &
            fftbox%total_ld2_dbl,2*fftbox%total_pt3,proj_count),stat=ierr)
       call utils_alloc_check('projectors_init_fftbox_recip', &
            'proj_set%fftbox_proj_recip',ierr)
    else
       allocate(proj_set%fftbox_proj_recip(fftbox%total_ld1, &
            fftbox%total_ld2,fftbox%total_pt3,proj_count),stat=ierr)
       call utils_alloc_check('projectors_init_fftbox_recip', &
            'proj_set%fftbox_proj_recip',ierr)
    end if


    ! cks: now initialise projectors in fftbox reciprocal representation
    tot_num_proj = proj_count

!$OMP PARALLEL DO NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(iproj,proj_count,found,isp,shell,proj_count_on_species,am, &
!$OMP      azim,proj_norm) &
!$OMP SHARED(tot_num_proj,proj_set,fftbox,cell,delta,cart,kpt_loc, &
!$OMP      loc_grad,pub_threads_num_fftboxes,fine)
    do iproj=1,tot_num_proj

       ! Find isp, shell, azim for this projector
       proj_count = 0
       found = .false.
       do isp=1,proj_set%n_proj_species
          proj_count_on_species = 0
          do shell=1,proj_set%num_shells(isp)
             am = proj_set%ang_mom(shell,isp)
             do azim=-am,am
                proj_count = proj_count + 1
                proj_count_on_species = proj_count_on_species + 1
                if (proj_count==iproj) then
                   found = .true.
                   exit
                end if
                ! Go to next species if we have set up all the projectors on
                ! this atom (even if there are more according to num_shells)
                if (proj_count_on_species>=proj_set%species_num_proj(isp)) exit
             end do
             if (found) exit
             if (proj_count_on_species>=proj_set%species_num_proj(isp)) exit
          end do
          if (found) exit
       end do
       if (.not.found) then
          call utils_abort('Error in projectors_init_fftbox_recip: Projector &
               &counting mismatch.')
       end if

       ! Put projector in reciprocal space box
       ! pdh: change to avoid copy in/out of fftbox_proj_recip subarray
       ! ndmh: change to allow use of proj_box_rec by other modules
       ! ewt23: Modified to accept the option of creating the projectors
       !        on the fine grid
       if (present(fine)) then
          if (.not.loc_grad) then
             call projector_in_box_recip( &
                  proj_set%fftbox_proj_recip(:,:,:,proj_count), &
                  proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                  shell,isp), am, azim, 2*fftbox%total_pt2, &
                  2*fftbox%total_pt2, 2*fftbox%total_pt3, &
                  proj_set%n_rad_pts(isp), proj_set%gmax(isp), cell, kpt_loc,&
                  fine='Fine')
          else
             ! agrecokpt
             call projector_in_box_recip_gradg( &
                  proj_set%fftbox_proj_recip(:,:,:,proj_count), &
                  proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                  shell,isp), am, azim, 2*fftbox%total_pt1, &
                  2*fftbox%total_pt2, 2*fftbox%total_pt3, &
                  proj_set%n_rad_pts(isp), proj_set%gmax(isp), cell, &
                  delta, cart,fine='Fine', kpt=kpt_loc)
          end if
       else
          if (.not.loc_grad) then
             call projector_in_box_recip( &
                  proj_set%fftbox_proj_recip(:,:,:,proj_count), &
                  proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                  shell,isp), am, azim, fftbox%total_pt1, &
                  fftbox%total_pt2, fftbox%total_pt3, &
                  proj_set%n_rad_pts(isp), proj_set%gmax(isp), cell, kpt_loc)
          else
             ! agrecokpt
             call projector_in_box_recip_gradg( &
                  proj_set%fftbox_proj_recip(:,:,:,proj_count), &
                  proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                  shell,isp), am, azim, fftbox%total_pt1, &
                  fftbox%total_pt2, fftbox%total_pt3, &
                  proj_set%n_rad_pts(isp), proj_set%gmax(isp), cell, &
                  delta, cart, kpt=kpt_loc)
          end if
       end if

       ! ndmh: normalise the projector if required
       ! ndmh: (only used for Hubbard hydrogenic projectors)
       if (proj_set%normalise) then
          if (.not. present(fine)) then
             proj_norm = sqrt(sum(abs(proj_set%fftbox_proj_recip(:,:,:, &
                  proj_count))**2)/real(fftbox%total_pt1 * &
                  fftbox%total_pt2*fftbox%total_pt3,kind=DP) &
                  / fftbox%weight)
             proj_set%fftbox_proj_recip(:,:,:,proj_count) = &
                  proj_set%fftbox_proj_recip(:,:,:,proj_count) / proj_norm
          else
             call utils_abort('Error in projectors_init_fftbox_recip: &
                  &Fine option specified for normalised proj set')
          end if
       end if

    end do
!$OMP END PARALLEL DO

    call timer_clock('projectors_init_fftbox_recip',2)

  end subroutine projectors_init_fftbox_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_exit_fftbox_recip(proj_set)

    !==============================================================!
    ! This subroutine deallocates a set of projectors in FFTboxes  !
    !--------------------------------------------------------------!
    ! Written by Nicholas Hine on 12/12/2011.                      !
    !==============================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(PROJECTOR_SET), intent(inout) :: proj_set

    ! Local Variables
    integer :: ierr      ! error flag

    ! Allocate projector FFTboxes
    deallocate(proj_set%fftbox_proj_recip,stat=ierr)
    call utils_dealloc_check('projectors_exit_fftbox_recip', &
         'proj_set%fftbox_proj_recip',ierr)

  end subroutine projectors_exit_fftbox_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projector_in_box_recip(fftbox_recip, &
       radial_projector,ang_mom,azimuthal_ang_mom,n1,n2,n3, &
       n_rad_pts, g_max, cell, kpt, fine)

    !==============================================================!
    ! This subroutine returns a single projector in an FFTbox      !
    ! in reciprocal space.                                         !
    !--------------------------------------------------------------!
    ! Written by Arash A. Mostofi in 2001 based on a subroutine    !
    ! developed earlier by Chris-Kriton Skylaris.                  !
    ! Modified by Arash A. Mostofi in April 2003 to work with      !
    ! complex-to-complex FFTs.                                     !
    ! Modified by Chris-Kriton Skylaris on 25/1/2004 to work with  !
    ! the pseudo_species type (in the calling subroutines).        !
    ! Extension to f-type projectors by Chris-Kriton Skylaris      !
    ! on 17/08/2007.                                               !
    ! Modified to use cell by Quintin Hill on 15/10/2008.      !
    ! Modified to use an output fftbox rather than only depositing !
    ! to the p_species array by  Nicholas Hine on 15/06/2009.      !
    ! Modified to use sw_real_sph_harm by Quintin Hill on          !
    ! 18/09/2009.                                                  !
    ! Moved to projectors_mod by Nicholas Hine on 19/05/2010.      !
    !==============================================================!

    use geometry, only: POINT, OPERATOR(*), OPERATOR(+), geometry_magnitude
    use rundat, only: pub_threads_per_fftbox
    use services, only: services_1d_interpolation
    use simulation_cell, only: CELL_INFO
    use spherical_wave, only: sw_real_sph_harm
    use utils, only: utils_abort, utils_use_var

    implicit none

    ! Arguments
    integer, intent(in) :: n1,n2,n3
    complex(kind=DP), intent(out) :: fftbox_recip(n1,n2,n3)
    integer, intent(in) :: ang_mom,azimuthal_ang_mom
    integer, intent(in) :: n_rad_pts
    real(kind=DP), intent(in) :: g_max
    real(kind=DP), intent(in)  :: radial_projector(n_rad_pts)
    type(CELL_INFO), intent(in)        :: cell
    type(POINT), optional, intent(in) :: kpt
    character(len=*), optional, intent(in) :: fine

    ! Local Variables
    integer     :: row1,row2,row3
    integer     :: K3,K2,K1
    type(POINT) :: g_vector3,g_vector2,g_vector1
    type(POINT) :: pairvec3,pairvec2,pairvec1
    real(kind=DP) :: g_length
    complex(kind=DP) :: phase_fac
    real(kind=DP) :: projector_value,harmonic_factor
    real(kind=DP), parameter :: four_pi=4.0_DP*PI

    pairvec3=( real(cell%total_pt3,DP)/real(n3,DP) ) * cell%b3
    pairvec2=( real(cell%total_pt2,DP)/real(n2,DP) ) * cell%b2
    pairvec1=( real(cell%total_pt1,DP)/real(n1,DP) ) * cell%b1

    if (present(fine)) then
       pairvec3=2.0_DP*pairvec3
       pairvec2=2.0_DP*pairvec2
       pairvec1=2.0_DP*pairvec1
    endif

    ! calculate the phase factor (-i)^l
    select case (ang_mom)
    case (0)
       phase_fac = cmplx(1.0_DP,0.0_DP,kind=DP)
    case (1)
       phase_fac = cmplx(0.0_DP,-1.0_DP,kind=DP)
    case (2)
       phase_fac = cmplx(-1.0_DP,0.0_DP,kind=DP)
    case (3)
       phase_fac = cmplx(0.0_DP,1.0_DP,kind=DP)
    case default
       call utils_abort('Error in projector_in_box_recip: Angular momentum &
            &too high. Only projectors up to f-type symmetry &
            &are currently supported.')
    end select

    call utils_use_var(pub_threads_per_fftbox)

    ! loop over the points of the reciprocal FFT box
!!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) DEFAULT(none) &
!!$OMP PRIVATE(row1, row2, row3, k1, k2, k3, g_vector3, g_vector2, g_vector1, &
!!$OMP         g_length, projector_value, harmonic_factor) &
!!$OMP SHARED(n1, n2, n3, pairvec3, pairvec2, pairvec1, kpt, &
!!$OMP        radial_projector, n_rad_pts, g_max, ang_mom, azimuthal_ang_mom, &
!!$OMP        fftbox_recip, phase_fac, pub_threads_per_fftbox)
!!$OMP DO
    do row3=0,n3-1

       k3=row3
       if ( row3.gt.(n3/2) ) k3=row3-n3

       g_vector3 = real(k3,kind=DP)*pairvec3
       if (present(kpt)) g_vector3 = g_vector3 + kpt

       do row2=0,n2-1

          k2=row2
          if ( row2.gt.(n2/2) ) k2=row2-n2

          g_vector2 = real(k2,kind=DP)*pairvec2 + g_vector3

          do row1=0,n1-1

             k1=row1
             if ( row1.gt.(n1/2) ) k1=row1-n1

             g_vector1 = real(k1,kind=DP)*pairvec1 + g_vector2

             g_length=geometry_magnitude(g_vector1)

             ! cks: interpolation of radial projector
             projector_value=services_1d_interpolation(  &
                  radial_projector(:), n_rad_pts, &
                  (g_length/g_max)*real(n_rad_pts-1, kind=DP), ang_mom)

             ! multiply value of projector at each g-point by the appropriate
             ! phase factor and the appropriate real spherical harmonic factor
             harmonic_factor = four_pi * sw_real_sph_harm(g_vector1%x, &
                  g_vector1%y, g_vector1%z, g_length, ang_mom, azimuthal_ang_mom)

             fftbox_recip(row1+1,row2+1,row3+1) &
                  = projector_value * harmonic_factor * phase_fac

          enddo
       enddo
    enddo
!!$OMP END PARALLEL

  end subroutine projector_in_box_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projector_in_box_recip_gradg(fftbox_recip, &
       radial_projector,ang_mom,azimuthal_ang_mom,n1,n2,n3, &
       n_rad_pts, g_max, cell, fin_diff_shift, cart, fine, kpt)

    !==============================================================!
    ! This subroutine calculates the gradient of the projectors    !
    ! with respect to G in reciprocal space.                       !
    !--------------------------------------------------------------!
    ! Written by Laura Ratcliff in March 2011 based on             !
    ! projector_in_box_recip.                                      !
    ! Modified by Andrea Greco to allow k-point specification,     !
    ! May 2016.                                                    !
    !==============================================================!

    use geometry, only: POINT, OPERATOR(*), OPERATOR(+), geometry_magnitude
    use rundat, only: pub_threads_per_fftbox
    use services, only: services_1d_interpolation
    use simulation_cell, only: CELL_INFO
    use spherical_wave, only: sw_real_sph_harm, sw_grad_real_sph_harm !lr408
    use utils, only: utils_abort, utils_use_var

    implicit none

    ! Arguments
    integer, intent(in) :: n1,n2,n3
    complex(kind=DP), intent(out) :: fftbox_recip(n1,n2,n3)
    integer, intent(in) :: ang_mom,azimuthal_ang_mom
    integer, intent(in) :: n_rad_pts
    real(kind=DP), intent(in) :: g_max
    real(kind=DP), intent(in)  :: radial_projector(n_rad_pts)
    real(kind=dp), intent(in) :: fin_diff_shift
    type(CELL_INFO), intent(in)        :: cell
    integer, intent(in) :: cart
    character(len=*), optional, intent(in) :: fine
    ! agrecokpt
    type(POINT), optional, intent(in) :: kpt

    ! Local Variables
    integer     :: row1,row2,row3
    integer     :: K3,K2,K1
    type(POINT) :: g_vector3,g_vector2,g_vector1
    type(POINT) :: pairvec3,pairvec2,pairvec1
    real(kind=DP) :: g_length
    complex(kind=DP) :: phase_fac
    real(kind=DP) :: projector_value,harmonic_factor
    real(kind=DP), parameter :: four_pi=4.0_DP*PI

    real(kind=DP) :: grad_sph_harm(3), grad_proj(3)
    real(kind=DP) :: g_length2, projector_value2, grad_projector_value

    ! agrecokpt
    type(POINT) :: loc_kpt

    ! agrecokpt
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    pairvec3=( real(cell%total_pt3,DP)/real(n3,DP) ) * cell%b3
    pairvec2=( real(cell%total_pt2,DP)/real(n2,DP) ) * cell%b2
    pairvec1=( real(cell%total_pt1,DP)/real(n1,DP) ) * cell%b1

    if (present(fine)) then
       pairvec3=2.0_DP*pairvec3
       pairvec2=2.0_DP*pairvec2
       pairvec1=2.0_DP*pairvec1
    endif

    ! calculate the phase factor (-i)^l
    select case (ang_mom)
    case (0)
       phase_fac = cmplx(1.0_DP,0.0_DP,kind=DP)
    case (1)
       phase_fac = cmplx(0.0_DP,-1.0_DP,kind=DP)
    case (2)
       phase_fac = cmplx(-1.0_DP,0.0_DP,kind=DP)
    case (3)
       phase_fac = cmplx(0.0_DP,1.0_DP,kind=DP)
    case default
       call utils_abort('Angular momentum too high too high in &
            &projector_in_box_recip_gradg. Only projectors up to f-type &
            &symmetry are currently supported.')
    end select

    call utils_use_var(pub_threads_per_fftbox)

    ! loop over the points of the reciprocal FFT box
!!$OMP PARALLEL DO NUM_THREADS(pub_threads_per_fftbox) DEFAULT(none) &
!!$OMP PRIVATE(row1, row2, row3, k1, k2, k3, g_vector3, g_vector2, g_vector1, &
!!$OMP         g_length, g_length2, projector_value, projector_value2, &
!!$OMP         grad_projector_value, harmonic_factor, grad_sph_harm, &
!!$OMP         grad_proj) &
!!$OMP SHARED(n1, n2, n3, pairvec3, pairvec2, pairvec1, pub_threads_per_fftbox, &
!!$OMP        radial_projector, n_rad_pts, g_max, ang_mom, azimuthal_ang_mom, &
!!$OMP        fftbox_recip, phase_fac, fin_diff_shift, cart, loc_kpt)
    do row3=0,n3-1

       k3=row3
       if ( row3.gt.(n3/2) ) k3=row3-n3

       ! agrecokpt
       g_vector3 = real(k3,kind=DP)*pairvec3 + loc_kpt

       do row2=0,n2-1

          k2=row2
          if ( row2.gt.(n2/2) ) k2=row2-n2

          g_vector2 = real(k2,kind=DP)*pairvec2 + g_vector3

          do row1=0,n1-1

             k1=row1
             if ( row1.gt.(n1/2) ) k1=row1-n1

             g_vector1 = real(k1,kind=DP)*pairvec1 + g_vector2

             g_length=geometry_magnitude(g_vector1)

             ! cks: interpolation of radial projector
             projector_value=services_1d_interpolation(  &
                  radial_projector(:), n_rad_pts, &
                  (g_length/g_max)*real(n_rad_pts-1, kind=DP), ang_mom)


             ! lr408: interpolate projector at a second point for finite diff
             g_length2 = g_length + fin_diff_shift

             projector_value2=services_1d_interpolation(  &
                  radial_projector(:), n_rad_pts, &
                  (g_length2/g_max)*real(n_rad_pts-1, kind=DP), ang_mom)

             grad_projector_value = (projector_value2 - projector_value) / fin_diff_shift

             harmonic_factor = four_pi * sw_real_sph_harm(g_vector1%x, &
                  g_vector1%y, g_vector1%z, g_length, ang_mom, azimuthal_ang_mom)


             if (g_length > 1d-10) then
                grad_proj(1) = grad_projector_value * harmonic_factor &
                     * g_vector1%x / g_length
                grad_proj(2) = grad_projector_value * harmonic_factor &
                     * g_vector1%y / g_length
                grad_proj(3) = grad_projector_value * harmonic_factor &
                     * g_vector1%z / g_length
             else
                if (ang_mom == 1) then  !special value for L=1: lim_{q->0} g(q)/q
                   projector_value = grad_projector_value
                   grad_proj(:) = 0.0_dp
                else
                   projector_value = 0.0_dp
                   grad_proj(:) = 0.0_dp
                end if
             end if

             ! multiply value of projector at each g-point by the appropriate
             ! phase factor (n.b. real part and imaginary part stored as
             ! separate consecutive elements in the x-direction of the array),
             ! and the appropriate real spherical harmonic factor.
             call sw_grad_real_sph_harm(grad_sph_harm, g_vector1%x, &
                  g_vector1%y, g_vector1%z, g_length, ang_mom, azimuthal_ang_mom)

             ! grad_proj * sph_harm
             fftbox_recip(row1+1,row2+1,row3+1) &
                  = grad_proj(cart) * phase_fac

             ! + proj * grad_sph_harm
             fftbox_recip(row1+1,row2+1,row3+1) &
                  = fftbox_recip(row1+1,row2+1,row3+1) + four_pi &
                  * projector_value * grad_sph_harm(cart) * phase_fac

          enddo
       enddo
    enddo
!!$OMP END PARALLEL DO

  end subroutine projector_in_box_recip_gradg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_create_real(proj_basis,proj_set, &
       fftbox,cell,projs_on_grid, is_cmplx, kpt)

    !===============================================================!
    ! This subroutine allocates space for and calculates a set of   !
    ! projectors in real space on PPDs                              !
    !---------------------------------------------------------------!
    ! Written by Nicholas Hine in March 2011.                       !
    !===============================================================!

    ! agrecocmplx
    use basis, only: basis_phase_on_fftbox_recip, basis_point_wrt_box, &
         basis_start_of_box_wrt_cell, basis_extract_function_from_box, &
         basis_find_function_wrt_box
    use comms, only: pub_my_proc_id
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_functions_alloc, &
         data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_debug_on_root, pub_projectors_precalculate, &
         pub_threads_fftbox, pub_kpoint_method
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_use_var, &
         utils_assert

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(PROJECTOR_SET), intent(inout) :: proj_set
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(CELL_INFO), intent(in) :: cell
    !real(kind=DP), intent(out), optional :: projs_on_grid(proj_basis%n_ppds* &
    !     cell%n_pts)
    ! agrecocmplx: use new FUNCTIONS type to allow complex valued projs_on_grid
    type(FUNCTIONS), intent(inout), optional :: projs_on_grid
    ! agrecocmplx: optional argument to specify if we want real or complex projs
    logical, intent(in), optional :: is_cmplx
    ! agrecokpt: optional argument to specify a single k-point
    real(kind=DP), intent(in), optional :: kpt(3)

    ! Local Variables
    integer :: local_proj
    integer :: global_proj
    integer :: atom_of_proj
    integer :: loc_atom_of_proj
    integer :: species_number
    integer :: proj_count, proj_count_atom
    integer :: ierr             ! Error flag
    integer :: n1,n2,n3,ld1,ld2
    integer :: proj_start(1:3)    ! Start of projector function in FFTbox
    integer :: fftbox_start(1:3)  ! Start of FFTbox in cell
    ! {P}rojector {C}entre wrt fft {B}ox in terms of number of {G}ridpoints
    ! in each lattice vector direction. Not necessarily integer values.
    real(kind=DP) :: pcbg1,pcbg2,pcbg3
    !real(kind=DP), allocatable :: projector_fftbox(:,:,:)
    ! agrecocmplx: use new FFTBOX_DATA type to allow complex values
    type(FFTBOX_DATA) :: projector_fftbox
    complex(kind=DP), allocatable :: proj_complex(:,:,:)   ! Workspace
    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt
    type(POINT) :: loc_kpt

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &projectors_create_real'

    ! agrecokpt: currently only KP method is being implemented
    call utils_assert(pub_kpoint_method == 'KP', &
         'Subroutine projectors_create_real currently supports&
         & only KP method for BZ sampling')

    ! Start timer
    call timer_clock('projectors_create_real',1)

    ! Initialise shorthand variables
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2
    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3

    ! agrecocmplx
    if (present(is_cmplx)) then
        loc_cmplx = is_cmplx
    else
        loc_cmplx = .false.
    end if

    ! agrecokpt: default case is Gamma point
    if (present(kpt)) then
       loc_kpt%x = kpt(1)
       loc_kpt%y = kpt(2)
       loc_kpt%z = kpt(3)
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    ! Allocate internal storage in projector set for ppds of projectors
    ! (only if we did not provide an array to hold the result)
    if (.not.present(projs_on_grid)) then
       !allocate(proj_set%projs_on_grid(proj_basis%size_on_grid), &
       !     stat=ierr)
       !call utils_alloc_check('projectors_create_real', &
       !     'proj_set%projs_on_grid',ierr)
       ! agrecocmplx: allocate using appropriate routine
        call data_functions_alloc(proj_set%projs_on_grid, &
             proj_basis%size_on_grid, loc_cmplx)
    end if

    ! agrecokpt: call it with single kpoint (defaults to Gamma)
    if (pub_projectors_precalculate) &
         call projectors_init_fftbox_recip(proj_set,cell,fftbox, &
              kpt=loc_kpt)

    ! Allocate workspace
    allocate(proj_complex(ld1, ld2, n3), stat=ierr)
    call utils_alloc_check('projectors_create_real','proj_complex',ierr)
    !allocate(projector_fftbox(ld1,ld2,n3),stat=ierr)
    !call utils_alloc_check('projectors_create_real','projector_fftbox',ierr)
    ! agrecocmplx: use appropriate routine to allocate FFTBOX_DATA
    call data_fftbox_alloc(projector_fftbox, ld1, ld2, n3, loc_cmplx)

    do local_proj=1,proj_basis%proc_num

       ! Find information about this projector
       global_proj = local_proj + proj_basis%first_on_proc(pub_my_proc_id) - 1
       atom_of_proj = proj_basis%atom_of_func(global_proj)
       loc_atom_of_proj = atom_of_proj - &
            proj_basis%atom_of_func(proj_basis%first_on_proc(pub_my_proc_id)) + 1
       species_number = proj_set%proj_species(atom_of_proj)
       proj_count_atom = global_proj - &
            proj_basis%first_on_atom(atom_of_proj) + 1
       proj_count = proj_set%species_first_proj(species_number) + &
            proj_count_atom - 1

       ! Centre of projector wrt fftbox in terms of grid spacings
       call basis_point_wrt_box(pcbg1, pcbg2, pcbg3, &  ! output
            proj_set%proj_centre(atom_of_proj), n1, n2, n3, cell, &
            fftbox)

       ! Start of fftbox wrt cell in terms of grid-point number
      call basis_start_of_box_wrt_cell(fftbox_start(1), &
           fftbox_start(2),fftbox_start(3), &
           proj_set%proj_centre(atom_of_proj),pcbg1,pcbg2,pcbg3,cell)

       ! Find position of tightbox start wrt fftbox
       call basis_find_function_wrt_box(proj_start(1), proj_start(2), &
            proj_start(3), fftbox_start(1),fftbox_start(2),fftbox_start(3), &
            proj_basis%tight_boxes(local_proj),cell, fftbox)

       ! Apply phase factor e^-(ik.proj_centre_wrt_box) to projector
       call basis_phase_on_fftbox_recip(proj_complex, &
            n1,n2,n3,ld1,ld2,-pcbg1,-pcbg2,-pcbg3, &
            proj_set%fftbox_proj_recip(:,:,:,proj_count))

       ! g=0 element must be real
       ! agrecokpt: but this is not true at a different kpt,
       ! because of the kshift applied? CHECK THIS!!!
       if (loc_kpt%x==0.0_DP.and.loc_kpt%y==0.0_DP.and.&
           loc_kpt%z==0.0_DP) then
          proj_complex(1,1,1) &
          = cmplx(real(proj_complex(1,1,1),kind=DP),0.0_DP,kind=DP)
       end if
       !proj_complex(1,1,1) &
       !   = cmplx(real(proj_complex(1,1,1),kind=DP),0.0_DP,kind=DP)

       ! Fourier transform to real space:
       call fourier_apply_box('Coarse','Backward',proj_complex, &
            omp=pub_threads_fftbox)

       ! Put real part into projector_box
       !projector_fftbox(1:n1,1:n2,1:n3) &
       !     = real(proj_complex(1:n1,1:n2,1:n3),kind=DP)/cell%weight
       ! agrecocmplx: mantain full complex character if needed
       if (projector_fftbox%iscmplx) then
           projector_fftbox%z(1:n1,1:n2,1:n3) &
               = proj_complex(1:n1,1:n2,1:n3)/cell%weight
       else
           projector_fftbox%d(1:n1,1:n2,1:n3) &
               = real(proj_complex(1:n1,1:n2,1:n3),kind=DP)/cell%weight
       end if

       ! Extract projector to ppds
       !if (present(projs_on_grid)) then
       !   call basis_extract_function_from_box(projs_on_grid(:), &
       !        ld1,ld2,n3,projector_fftbox,proj_basis%spheres(local_proj), &
       !        proj_basis%tight_boxes(local_proj),proj_start(1),proj_start(2), &
       !        proj_start(3),proj_basis%spheres(local_proj)%offset,cell, &
       !        fftbox)
       !else
       !   call basis_extract_function_from_box(proj_set%projs_on_grid(:), &
       !        ld1,ld2,n3,projector_fftbox,proj_basis%spheres(local_proj), &
       !        proj_basis%tight_boxes(local_proj),proj_start(1),proj_start(2), &
       !        proj_start(3),proj_basis%spheres(local_proj)%offset,cell, &
       !        fftbox)
       !end if

       ! agrecocmplx: use new routines in basis_new_mod
       if (present(projs_on_grid)) then
          call basis_extract_function_from_box(projs_on_grid, &
               projector_fftbox,proj_basis%spheres(local_proj), &
               proj_basis%tight_boxes(local_proj),proj_start(1),proj_start(2), &
               proj_start(3),proj_basis%spheres(local_proj)%offset,cell, &
               fftbox)
       else
          call basis_extract_function_from_box(proj_set%projs_on_grid, &
               projector_fftbox,proj_basis%spheres(local_proj), &
               proj_basis%tight_boxes(local_proj),proj_start(1),proj_start(2), &
               proj_start(3),proj_basis%spheres(local_proj)%offset,cell, &
               fftbox)
       end if

    end do

    ! Deallocate workspace
    !deallocate(projector_fftbox,stat=ierr)
    !call utils_dealloc_check('projectors_create_real','projector_fftbox',ierr)
    ! agrecocmplx: use appropriate routine to deallocate FFTBOX_DATA
    call data_fftbox_dealloc(projector_fftbox)
    deallocate(proj_complex,stat=ierr)
    call utils_dealloc_check('projectors_create_real','proj_complex',ierr)

    if (pub_projectors_precalculate) call projectors_exit_fftbox_recip(proj_set)

    ! Stop timer
    call timer_clock('projectors_create_real',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &projectors_create_real'

  end subroutine projectors_create_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_func_ovlp_box(sp_overlap, &
       ngwfs_on_grid,ngwf_basis,proj_basis,proj_set,fftbox,cell, &
       kshift,delta,cart,factor_type)

    !===============================================================!
    ! This subroutine calculates the ngwf-projector overlap matrix  !
    ! and stores it in a SPAM3 matrix.                              !
    !---------------------------------------------------------------!
    ! Written by Nicholas Hine in June 2009.                        !
    ! Some code from reused from the old version of                 !
    ! pseudopotentials_sp_ovlp_box written by Arash A. Mostofi in   !
    ! January 2004 and modified by Chris-Kriton Skylaris, Arash     !
    ! Mostofi and Nicholas Hine in 2004-2008.                       !
    ! Moved to projectors_mod by Nicholas Hine on 19/05/2010.       !
    ! Optional k-vector added as input for calculating overlap at a !
    ! different k-point by Laura Ratcliff 17/2/11.                  !
    ! Modified by Andrea Greco on 23/05/2015 to allow use of        !
    ! complex NGWFs.                                                !
    !===============================================================!

    use basis, only: basis_phase_on_fftbox_recip, basis_point_wrt_box, &
         basis_start_of_box_wrt_cell
    use comms, only: pub_my_proc_id
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketfftbox
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_fftbox_batch_size, pub_projectors_precalculate, &
         pub_threads_fftbox, pub_debug_on_root, pub_kpoint_method
!$  use rundat, only: pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_use_var, &
         utils_assert

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: sp_overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(PROJECTOR_SET), intent(inout) :: proj_set
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(CELL_INFO), intent(in) :: cell
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    type(POINT), optional, intent(in) :: kshift
    real(kind=DP), optional, intent(in) :: delta
    integer, optional, intent(in) :: cart
    integer, optional, intent(in) :: factor_type

    ! Local Variables
    integer :: n_batches
    integer :: batch_size
    integer :: batch_count
    integer :: local_start
    integer :: local_end
    integer :: local_proj
    integer :: global_proj
    integer :: batch_proj
    integer :: atom_of_proj
    integer :: loc_atom_of_proj
    integer :: isp
    integer :: iproj, shell, am, azim
    integer :: proj_count, proj_count_atom
    integer :: ierr             ! Error flag
    integer :: idx_len          ! Length of sparse index
    integer :: n1,n2,n3,ld1,ld2
    logical :: found
    ! {P}rojector {C}entre wrt fft {B}ox in terms of number of {G}ridpoints
    ! in each lattice vector direction. Not necessarily integer values.
    real(kind=DP) :: pcbg1,pcbg2,pcbg3
    type(FFTBOX_DATA), allocatable :: projector_fftbox(:)
    integer, allocatable :: fftbox_start(:,:)              ! Ket start positions
    integer, allocatable :: sp_overlap_idx(:)              ! Sparse indexs
    complex(kind=DP), allocatable :: proj_complex(:,:,:)   ! Workspace
    character(len=*), parameter :: myself = 'projectors_func_ovlp_box'
    logical :: ongpu
    type(POINT) :: loc_kpt
    type(PARAL_INFO), pointer :: par
    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &projectors_func_ovlp_box'

    ! agrecokpt: currently only KP method is being implemented
    call utils_assert(pub_kpoint_method == 'KP', &
         'Subroutine projectors_func_ovlp_box currently supports&
         & only KP method for BZ sampling')

    ! Start timer
    call timer_clock('projectors_func_ovlp_box',1)

    call utils_assert(sp_overlap%iscmplx .or. (.not. ngwfs_on_grid%iscmplx), &
         'Error in '//myself//': complex projectors require real output &
         &sp_overlap')

    ! rc2013: get paralle strategy from sparse matrix
    call sparse_get_par(par, sp_overlap, 'R')

    ! agrecokpt: default case is Gamma point
    if (present(kshift)) then
       loc_kpt%x = kshift%x
       loc_kpt%y = kshift%y
       loc_kpt%z = kshift%z
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    ! Initialise shorthand variables
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2
    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3

    ! Obtain index of sp_overlap
    idx_len = sparse_index_length(sp_overlap)
    allocate(sp_overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('projectors_func_ovlp_box','sp_overlap_idx',ierr)
    call sparse_generate_index(sp_overlap_idx,sp_overlap)

    ! Set batch size
    batch_size = pub_fftbox_batch_size

    ! Allocate workspace
    allocate(projector_fftbox(batch_size), stat=ierr)
    call utils_alloc_check('projectors_func_ovlp_box', &
         'projector_fftbox',ierr)
    do batch_count = 1, batch_size
        call data_fftbox_alloc(projector_fftbox(batch_count), ld1, ld2, &
             n3, iscmplx=sp_overlap%iscmplx)
    end do
    allocate(fftbox_start(3,batch_size),stat=ierr)
    call utils_alloc_check('projectors_func_ovlp_box', &
         'fftbox_start',ierr)

    if (pub_projectors_precalculate) &
         call projectors_init_fftbox_recip(proj_set,cell,fftbox,kshift,delta,cart)

    ! Number of row-steps per row-block
    n_batches = proj_basis%max_on_proc / batch_size
    if (mod(proj_basis%max_on_proc,batch_size) > 0) n_batches = n_batches + 1

    local_start = 1
    do batch_count=1,n_batches
       local_end = min(local_start+batch_size-1,proj_basis%proc_num)

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(local_proj,batch_proj,global_proj,atom_of_proj, &
!$OMP      loc_atom_of_proj,isp,proj_count_atom,proj_count,pcbg1,pcbg2,pcbg3, &
!$OMP      proj_complex,ierr,iproj,found,am,azim,shell) &
!$OMP SHARED(local_start,local_end,proj_set,ld1,ld2,n1,n2,n3,batch_size, &
!$OMP      proj_basis,pub_my_proc_id,par,fftbox, &
!$OMP      cell,pub_projectors_precalculate,kshift, &
!$OMP      fftbox_start,projector_fftbox,pub_threads_fftbox, &
!$OMP      pub_threads_num_fftboxes,loc_kpt)
       allocate(proj_complex(ld1, ld2, n3), stat=ierr)
       call utils_alloc_check('projectors_func_ovlp_box', &
            'proj_complex',ierr)

!$OMP DO
       do local_proj=local_start,local_end

          batch_proj = local_proj - local_start + 1

          ! Find information about this projector
          global_proj = local_proj + proj_basis%first_on_proc(pub_my_proc_id) &
               - 1
          atom_of_proj = proj_basis%atom_of_func(global_proj)
          loc_atom_of_proj = atom_of_proj - &
               proj_basis%atom_of_func(proj_basis%first_on_proc(pub_my_proc_id)) + 1
          isp = proj_set%proj_species(atom_of_proj)
          proj_count_atom = global_proj - &
               proj_basis%first_on_atom(atom_of_proj) + 1
          proj_count = proj_set%species_first_proj(isp) + &
               proj_count_atom - 1

          ! Centre of projector wrt fftbox in terms of grid spacings
          call basis_point_wrt_box(pcbg1, pcbg2, pcbg3, &  ! output
               proj_set%proj_centre(atom_of_proj), n1, n2, n3, cell, &
               fftbox)

          ! Start of fftbox wrt cell in terms of grid-point number
          call basis_start_of_box_wrt_cell(fftbox_start(1,batch_proj), &
               fftbox_start(2,batch_proj),fftbox_start(3,batch_proj), &
               proj_set%proj_centre(atom_of_proj),pcbg1,pcbg2,pcbg3, &
               cell)

          if (pub_projectors_precalculate) then

             ! Apply phase factor e^-(ik.proj_centre_wrt_box) to projector
             ! stored in pre-calculated FFTbox
             call basis_phase_on_fftbox_recip(proj_complex, &
                  n1,n2,n3,ld1,ld2,-pcbg1,-pcbg2,-pcbg3, &
                  proj_set%fftbox_proj_recip(:,:,:,proj_count))
          else

             ! Find angular momentum of this projector
             iproj = 1
             found = .false.
             do shell=1,proj_set%num_shells(isp)
                am = proj_set%ang_mom(shell,isp)
                do azim=-am,am
                   if (iproj==proj_count_atom) found=.true.
                   if (found) exit
                   iproj = iproj + 1
                end do
                if (found) exit
             end do

             ! Generate projector in FFTbox
             call projector_in_box_recip(proj_complex, &
                  proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                  shell,isp), am, azim, fftbox%total_pt1, &
                  fftbox%total_pt2, fftbox%total_pt3, &
                  proj_set%n_rad_pts(isp), proj_set%gmax(isp), cell, &
                  kshift)

             ! Apply phase factor e^-(ik.proj_centre_wrt_box) to projector
             call basis_phase_on_fftbox_recip(proj_complex, &
                  n1,n2,n3,ld1,ld2,-pcbg1,-pcbg2,-pcbg3)

          end if

          ! g=0 element must be real
          ! agrecokpt: but this is not true at a different kpt,
          ! because of the kshift applied? CHECK THIS!!!
          if (loc_kpt%x==0.0_DP.and.loc_kpt%y==0.0_DP.and.&
              loc_kpt%z==0.0_DP) then
             proj_complex(1,1,1) &
             = cmplx(real(proj_complex(1,1,1),kind=DP),0.0_DP,kind=DP)
          end if
          !proj_complex(1,1,1) &
          !     = cmplx(real(proj_complex(1,1,1),kind=DP),0.0_DP,kind=DP)

          ! Fourier transform to real space:
          call fourier_apply_box('Coarse','Backward',proj_complex, &
               omp=pub_threads_fftbox)

          ! Put real part into projector_box
          if (projector_fftbox(batch_proj)%iscmplx) then
              projector_fftbox(batch_proj)%z(1:n1,1:n2,1:n3) &
                 = proj_complex(1:n1,1:n2,1:n3)/cell%weight
          else
              projector_fftbox(batch_proj)%d(1:n1,1:n2,1:n3) &
                 = real(proj_complex(1:n1,1:n2,1:n3),kind=DP)/cell%weight
          end if

       end do
!$OMP END DO

       deallocate(proj_complex,stat=ierr)
       call utils_dealloc_check('projectors_func_ovlp_box', &
            'proj_complex',ierr)
!$OMP END PARALLEL

       !  Calculate overlap integrals
       ongpu=.false.
       call function_ops_brappd_ketfftbox(sp_overlap,  &        ! inout
            ngwfs_on_grid, ngwf_basis, cell, fftbox, &          ! input
            projector_fftbox, fftbox_start, batch_size,  &      ! input
            local_start, local_end, idx_len, sp_overlap_idx, &  ! input
            & 'FULL', ongpu, factor_type)    ! input

       local_start = local_start + batch_size

    end do

    if (pub_projectors_precalculate) call projectors_exit_fftbox_recip(proj_set)

    ! Deallocate workspace
    deallocate(fftbox_start,stat=ierr)
    call utils_dealloc_check('projectors_func_ovlp_box', &
         'fftbox_start',ierr)
    do batch_count = 1, batch_size
        call data_fftbox_dealloc(projector_fftbox(batch_count))
    end do
    deallocate(projector_fftbox,stat=ierr)
    call utils_dealloc_check('projectors_func_ovlp_box', &
         'projector_fftbox',ierr)
    deallocate(sp_overlap_idx,stat=ierr)
    call utils_dealloc_check('projectors_func_ovlp_box', &
         'sp_overlap_idx',ierr)

    ! Stop timer
    call timer_clock('projectors_func_ovlp_box',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &projectors_func_ovlp_box'

  end subroutine projectors_func_ovlp_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_func_grad_ovlp_box(siGp_overlap, &
       ngwfs_on_grid,ngwf_basis,proj_basis,proj_set,fftbox,cell,cart,kpt)

    !========================================================================!
    ! This subroutine calculates the overlap matrix between the NGWFs and    !
    ! the first derivatives of the projectors in each of the Cartesian       !
    ! directions, and stores it in a SPAM3 matrix.                           !
    !------------------------------------------------------------------------!
    ! Arguments                                                              !
    ! 1)  siGp_overlap  : output : <ngwf|iG*proj> overlap matrix             !
    ! 2)  ngwfs_on_grid : input  : ngwf values on the grid                   !
    ! 3)  ngwfs_basis   : input  : function basis type for the NGWFs         !
    !------------------------------------------------------------------------!
    ! Written by Nicholas Hine in July 2009.                                 !
    ! Some code reused from the old version of                               !
    ! pseudopotentials_sp_ovlp_box written by Arash A. Mostofi in            !
    ! January 2004 and modified by Chris-Kriton Skylaris, Arash              !
    ! Mostofi and Nicholas Hine in 2004-2008.                                !
    ! Moved to projectors_mod by Nicholas Hine on 19/05/2010.                !
    !========================================================================!

    use basis, only: basis_phase_on_fftbox_recip, basis_point_wrt_box, &
         basis_start_of_box_wrt_cell
    use comms, only: pub_my_proc_id
    use constants, only: cmplx_i
    use datatypes, only: FFTBOX_DATA, data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketfftbox
    use geometry, only: POINT, operator(*), operator(+)
    use rundat, only: pub_fftbox_batch_size, pub_projectors_precalculate, &
         pub_threads_fftbox, pub_debug_on_root
!$  use rundat, only: pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_use_var

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: siGp_overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(PROJECTOR_SET), intent(inout) :: proj_set
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(CELL_INFO), intent(in) :: cell
    !real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%size_on_grid)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    integer, intent(in) :: cart
    ! agrecokpt
    type(POINT), optional, intent(in) :: kpt

    ! Local Variables
    integer :: n_batches
    integer :: batch_size
    integer :: batch_count
    integer :: local_start
    integer :: local_end
    integer :: local_proj
    integer :: global_proj
    integer :: batch_proj
    integer :: atom_of_proj
    integer :: loc_atom_of_proj
    integer :: isp
    integer :: iproj, shell, am, azim
    integer :: proj_count, proj_count_atom
    integer :: ierr             ! Error flag
    integer :: idx_len          ! Length of sparse index
    integer :: n1,n2,n3,ld1,ld2
    logical :: found
    ! {P}rojector {C}entre wrt fft {B}ox in terms of number of {G}ridpoints
    ! in each lattice vector direction. Not necessarily integer values.
    real(kind=DP) :: pcbg1,pcbg2,pcbg3
    !real(kind=DP), allocatable :: projector_fftbox(:,:,:,:)
    type(FFTBOX_DATA), allocatable :: projector_fftbox(:)
    integer, allocatable :: fftbox_start(:,:)              ! Ket start positions
    integer, allocatable :: sp_overlap_idx(:)              ! Sparse index
    complex(kind=DP), allocatable :: proj_complex(:,:,:)   ! Workspace
    logical :: ongpu
    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt
    type(POINT) :: loc_kpt

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &projectors_func_grad_ovlp_box'

    ! Start timer
    call timer_clock('projectors_func_grad_ovlp_box',1)

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid%iscmplx

    ! agrecokpt
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    ! Initialise shorthand variables
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2
    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3

    ! Obtain index of siGp_overlap
    idx_len = sparse_index_length(siGp_overlap)
    allocate(sp_overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('projectors_func_grad_ovlp_box','sp_overlap_idx',ierr)
    call sparse_generate_index(sp_overlap_idx,siGp_overlap)

    ! Set batch size
    batch_size = pub_fftbox_batch_size

    ! Allocate workspace
    !allocate(projector_fftbox(ld1,ld2,n3,batch_size),stat=ierr)
    allocate(projector_fftbox(batch_size), stat=ierr)
    call utils_alloc_check('projectors_func_grad_ovlp_box', &
         'projector_fftbox',ierr)
    do batch_count = 1, batch_size
        call data_fftbox_alloc(projector_fftbox(batch_count), ld1, ld2, n3, &
             iscmplx=loc_cmplx)
    end do
    allocate(fftbox_start(3,batch_size),stat=ierr)
    call utils_alloc_check('projectors_func_grad_ovlp_box', &
         'fftbox_start',ierr)

    if (pub_projectors_precalculate) &
         ! agrecokpt
         call projectors_init_fftbox_recip(proj_set,cell,fftbox,kpt=loc_kpt)

    ! Number of row-steps per row-block
    n_batches = proj_basis%max_on_proc / batch_size
    if (mod(proj_basis%max_on_proc,batch_size) > 0) n_batches = n_batches + 1

    local_start = 1
    do batch_count=1,n_batches
       local_end = min(local_start+batch_size-1,proj_basis%proc_num)

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(local_proj,batch_proj,global_proj,atom_of_proj, &
!$OMP      loc_atom_of_proj,isp,proj_count_atom,proj_count,pcbg1,pcbg2,pcbg3, &
!$OMP      proj_complex,ierr,iproj,found,am,azim,shell) &
!$OMP SHARED(local_start,local_end,proj_set,ld1,ld2,n1,n2,n3,batch_size, &
!$OMP      proj_basis,pub_my_proc_id,fftbox, &
!$OMP      cell,pub_projectors_precalculate,cart,pub_threads_num_fftboxes, &
!$OMP      fftbox_start,projector_fftbox,pub_threads_fftbox,loc_cmplx,loc_kpt)

       allocate(proj_complex(ld1, ld2, n3), stat=ierr)
       call utils_alloc_check('projectors_func_grad_ovlp_box', &
            'proj_complex',ierr)

!$OMP DO
       do local_proj=local_start,local_end

          batch_proj = local_proj - local_start + 1

          ! Find information about this projector
          global_proj = local_proj + proj_basis%first_on_proc(pub_my_proc_id) - 1
          atom_of_proj = proj_basis%atom_of_func(global_proj)
          loc_atom_of_proj = atom_of_proj - &
               proj_basis%atom_of_func(proj_basis%first_on_proc(pub_my_proc_id)) + 1
          isp = proj_set%proj_species(atom_of_proj)
          proj_count_atom = global_proj - &
               proj_basis%first_on_atom(atom_of_proj) + 1
          proj_count = proj_set%species_first_proj(isp) + &
               proj_count_atom - 1

          ! Centre of projector wrt fftbox in terms of grid spacings
          call basis_point_wrt_box(pcbg1, pcbg2, pcbg3, &  ! output
               proj_set%proj_centre(atom_of_proj), n1, n2, n3, cell, &
               fftbox)

          ! Start of fftbox wrt cell in terms of grid-point number
          call basis_start_of_box_wrt_cell(fftbox_start(1,batch_proj), &
               fftbox_start(2,batch_proj), fftbox_start(3,batch_proj), &
               proj_set%proj_centre(atom_of_proj), pcbg1, pcbg2, &
               pcbg3, cell)

          if (pub_projectors_precalculate) then

             ! Apply phase factor e^-(ik.proj_centre_wrt_box) to projector
             ! stored in pre-calculated FFTbox
             call basis_phase_on_fftbox_recip(proj_complex, &
                  n1,n2,n3,ld1,ld2,-pcbg1,-pcbg2,-pcbg3, &
                  proj_set%fftbox_proj_recip(:,:,:,proj_count))
          else

             ! Find angular momentum of this projector
             iproj = 1
             found = .false.
             do shell=1,proj_set%num_shells(isp)
                am = proj_set%ang_mom(shell,isp)
                do azim=-am,am
                   if (iproj==proj_count_atom) found=.true.
                   if (found) exit
                   iproj = iproj + 1
                end do
                if (found) exit
             end do

             ! Generate projector in FFTbox
             ! agrecokpt
             call projector_in_box_recip(proj_complex, &
                  proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                  shell,isp), am, azim, fftbox%total_pt1, &
                  fftbox%total_pt2, fftbox%total_pt3, &
                  proj_set%n_rad_pts(isp), proj_set%gmax(isp), cell, &
                  kpt=loc_kpt)

             ! Apply phase factor e^-(ik.proj_centre_wrt_box) to projector
             call basis_phase_on_fftbox_recip(proj_complex, &
                  n1,n2,n3,ld1,ld2,-pcbg1,-pcbg2,-pcbg3,proj_complex)

          end if

          ! Multiply projector (workspace) by iG and store as proj_complex
          proj_complex(:,:,:) = &
               proj_complex(:,:,:)*fftbox%recip_grid(cart,:,:,:)*cmplx_i

          ! g=0 element must be real
          ! agrecokpt: but this is not true at a different kpt,
          ! because of the kshift applied? CHECK THIS!!!
          if (loc_kpt%x==0.0_DP.and.loc_kpt%y==0.0_DP.and.&
              loc_kpt%z==0.0_DP) then
             proj_complex(1,1,1) &
                = cmplx(real(proj_complex(1,1,1),kind=DP),0.0_DP,kind=DP)
          end if
          !proj_complex(1,1,1) &
          !   = cmplx(real(proj_complex(1,1,1),kind=DP),0.0_DP,kind=DP)

          ! Fourier transform to real space:
          ! ndmh: explicitly in-place transform
          call fourier_apply_box('Coarse','Backward',proj_complex, &
               omp=pub_threads_fftbox)

          ! Put real part into projector_box
          !projector_fftbox(1:n1,1:n2,1:n3,batch_proj) &
          !     = real(proj_complex(1:n1,1:n2,1:n3),kind=DP)/cell%weight
          if (loc_cmplx) then
              projector_fftbox(batch_proj)%z(1:n1,1:n2,1:n3) &
                 = proj_complex(1:n1,1:n2,1:n3)/cell%weight
          else
              projector_fftbox(batch_proj)%d(1:n1,1:n2,1:n3) &
                 = real(proj_complex(1:n1,1:n2,1:n3),kind=DP)/cell%weight
          end if

       end do
!$OMP END DO

       deallocate(proj_complex,stat=ierr)
       call utils_dealloc_check('projectors_func_grad_ovlp_box', &
            'proj_complex',ierr)

!$OMP END PARALLEL

       ! Calculate overlap integrals
       ongpu=.false.
       call function_ops_brappd_ketfftbox(siGp_overlap, &           ! inout
            ngwfs_on_grid, ngwf_basis, cell, fftbox, &   ! input
            projector_fftbox, fftbox_start, &           ! input
            batch_size, local_start, local_end, idx_len, &       ! input
            sp_overlap_idx, 'FULL',ongpu)

       local_start = local_start + batch_size

    end do

    if (pub_projectors_precalculate) call projectors_exit_fftbox_recip(proj_set)

    ! Deallocate workspace
    deallocate(fftbox_start,stat=ierr)
    call utils_dealloc_check('projectors_func_grad_ovlp_box', &
         'fftbox_start',ierr)
    do batch_count = 1, batch_size
        call data_fftbox_dealloc(projector_fftbox(batch_count))
    end do
    deallocate(projector_fftbox,stat=ierr)
    call utils_dealloc_check('projectors_func_grad_ovlp_box', &
         'projector_fftbox',ierr)
    deallocate(sp_overlap_idx,stat=ierr)
    call utils_dealloc_check('projectors_func_grad_ovlp_box', &
         'sp_overlap_idx',ierr)

    ! Stop timer
    call timer_clock('projectors_func_grad_ovlp_box',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &projectors_func_grad_ovlp_box'

  end subroutine projectors_func_grad_ovlp_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_func_pos_ovlp_box(sp_overlap, &
       ngwfs_on_grid,ngwf_basis,proj_basis,proj_set,fftbox,cell,kpt)

    !===============================================================!
    ! This subroutine calculates the overlap matrix between a set   !
    ! of NGWFs, and a set of projectors, with the position operator !
    ! in between the two, and stores the result in a SPAM3 matrix.  !
    !---------------------------------------------------------------!
    ! Written by Laura Ratcliff in November 2011, based on          !
    ! the routine projectors_func_ovlp_box by Nicholas Hine.        !
    !===============================================================!

    use basis, only: basis_phase_on_fftbox_recip, basis_point_wrt_box, &
         basis_start_of_box_wrt_cell
    use comms, only: pub_my_proc_id
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box, fourier_filter, fourier_interpolate
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketfftbox
    use geometry, only: POINT, operator(*), operator(+)
    use rundat, only: pub_fftbox_batch_size, pub_projectors_precalculate, &
         pub_threads_fftbox, pub_debug_on_root, pub_eels_fine_projectors
!$  use rundat, only: pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_use_var, &
         utils_assert

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: sp_overlap(3)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(PROJECTOR_SET), intent(inout) :: proj_set
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(CELL_INFO), intent(in) :: cell
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    ! agrecokpt
    type(POINT), optional, intent(in) :: kpt

    ! Local Variables
    type(POINT) :: r1,r2,r3 ! lr408
    type(POINT) :: a1,a2,a3 ! lr408
    integer :: i1,i2,i3,count

    integer :: n_batches
    integer :: batch_size
    integer :: batch_count
    integer :: local_start
    integer :: local_end
    integer :: local_proj
    integer :: global_proj
    integer :: batch_proj, old_batch_proj
    integer :: atom_of_proj
    integer :: loc_atom_of_proj
    integer :: isp
    integer :: iproj, shell, am, azim
    integer :: proj_count, proj_count_atom
    integer :: ierr             ! Error flag
    integer :: idx_len          ! Length of sparse index
    integer :: n1,n2,n3,ld1,ld2
    integer :: ld1_dbl,ld2_dbl ! lr408
    ! {P}rojector {C}entre wrt fft {B}ox in terms of number of {G}ridpoints
    ! in each lattice vector direction. Not necessarily integer values.
    integer:: xyz
    logical :: found
    real(kind=DP) :: pcbg1,pcbg2,pcbg3
    real(kind=DP) :: old_pcbg1,old_pcbg2,old_pcbg3
    type(FFTBOX_DATA), allocatable :: projector_fftbox(:)
    integer, allocatable :: fftbox_start(:,:)              ! Ket start positions
    integer, allocatable :: sp_overlap_idx(:)              ! Sparse index
    complex(kind=DP), allocatable :: proj_complex(:,:,:)   ! Workspace
    real(kind=DP), allocatable :: old_proj(:,:,:)   ! Workspace
    real(kind=DP), allocatable :: new_proj(:,:,:)   ! Workspace
    real(kind=DP), dimension(:,:,:,:), allocatable :: r_op
    real(kind=DP), dimension(:,:,:), allocatable :: fftbox1_dbl,fftbox2_dbl
    ! agrecocmplx
    complex(kind=DP), dimension(:,:,:), allocatable :: fftbox_cmplx_dbl
    complex(kind=DP), dimension(:,:,:), allocatable :: proj_cmplx_std
    complex(kind=DP), dimension(:,:,:), allocatable :: coarse_work,fine_work
    integer :: thread_id,nthreads,local_size
    integer :: thread_local_start, thread_local_end
    real(kind=DP) :: bc1,bc2,bc3,bc1_dbl,bc2_dbl,bc3_dbl
    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt
    type(POINT) :: loc_kpt
!$  integer, external :: omp_get_thread_num, omp_get_num_threads
    logical :: ongpu

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &projectors_func_pos_ovlp_box'

    !jmecmplx
    ! agrecocmplx: now made compatible, need to check everything works ok!
    !call utils_assert(.not. ngwfs_on_grid%iscmplx, &
    !     'projectors_func_pos_ovlp_box not ready for complex NGWFs yet.')

    ! Start timer
    call timer_clock('projectors_func_pos_ovlp_box',1)

    ! Initialise shorthand variables
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2
    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3

    ld1_dbl = fftbox%total_ld1_dbl ! lr408
    ld2_dbl = fftbox%total_ld2_dbl ! lr408
    old_batch_proj = -999

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid%iscmplx

    ! agrecokpt
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    ! Obtain index of sp_overlap
    idx_len = sparse_index_length(sp_overlap(1))
    allocate(sp_overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box','sp_overlap_idx',ierr)
    call sparse_generate_index(sp_overlap_idx,sp_overlap(1))

    ! Set batch size
    batch_size = pub_fftbox_batch_size

    ! Allocate workspace
    allocate(projector_fftbox(batch_size),stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box', &
         'projector_fftbox',ierr)
    do batch_count = 1, batch_size
       call data_fftbox_alloc(projector_fftbox(batch_count), ld1, ld2, n3, &
            iscmplx=loc_cmplx)
    end do
    allocate(fftbox_start(3,batch_size),stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box', &
         'fftbox_start',ierr)
    allocate(r_op(3,ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box','r_op',ierr)


    if (pub_projectors_precalculate) then
        if(pub_eels_fine_projectors) then
             ! agrecokpt
             call projectors_init_fftbox_recip(proj_set,&
                 cell,fftbox,kpt=loc_kpt,fine="Fine")
        else
             ! agrecokpt
             call projectors_init_fftbox_recip(proj_set,&
                  cell,fftbox,kpt=loc_kpt)
        end if
    end if

    ! calculate vectors between fine grid points
    a1 = (1.0_DP / fftbox%total_pt1_dbl) * fftbox%a1
    a2 = (1.0_DP / fftbox%total_pt2_dbl) * fftbox%a2
    a3 = (1.0_DP / fftbox%total_pt3_dbl) * fftbox%a3

    ! ddor: construct position operatorin FFTbox
    r_op = 0.0_DP
    do i3=1,2*n3
       r3 = real(i3-1,kind=DP) * a3
       do i2=1,ld2_dbl
          r2 = r3 + real(i2-1,kind=DP) * a2
          do i1=1,ld1_dbl
             r1 = r2 + real(i1-1,kind=DP) * a1
             r_op(1,i1,i2,i3) = (r1%X)
             r_op(2,i1,i2,i3) = (r1%Y)
             r_op(3,i1,i2,i3) = (r1%Z)
          enddo
       enddo
    enddo


    bc1 = real(ld1,kind=DP)/2.0_DP
    bc2 = real(ld2,kind=DP)/2.0_DP
    bc3 = real(n3,kind=DP)/2.0_DP


    bc1_dbl = real(ld1_dbl,kind=DP)/2.0_DP
    bc2_dbl = real(ld2_dbl,kind=DP)/2.0_DP
    bc3_dbl = real(2*n3,kind=DP)/2.0_DP

    ! Number of row-steps per row-block
    n_batches = proj_basis%max_on_proc / batch_size
    if (mod(proj_basis%max_on_proc,batch_size) > 0) n_batches = n_batches + 1

    ! loop over x,y,z
    do xyz=1,3

       local_start = 1
       do batch_count=1,n_batches
          local_end = min(local_start+batch_size-1,proj_basis%proc_num)
          local_size = local_end - local_start + 1

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(local_proj,batch_proj,global_proj,atom_of_proj,count, &
!$OMP      loc_atom_of_proj,isp,proj_count_atom,proj_count,pcbg1,pcbg2,pcbg3, &
!$OMP      proj_complex,ierr,iproj,found,am,azim,shell,old_proj,new_proj, &
!$OMP      thread_local_start,thread_local_end,thread_id,nthreads, &
!$OMP      fftbox1_dbl,fftbox2_dbl,old_batch_proj,coarse_work,fine_work, &
!$OMP      old_pcbg1,old_pcbg2,old_pcbg3,fftbox_cmplx_dbl,proj_cmplx_std) &
!$OMP SHARED(local_start,local_end,proj_set,ld1,ld2,n1,n2,n3,batch_size, &
!$OMP      proj_basis,pub_my_proc_id,fftbox,local_size, &
!$OMP      cell,pub_projectors_precalculate,xyz,r_op,loc_cmplx, &
!$OMP      fftbox_start,projector_fftbox,ld1_dbl,ld2_dbl,pub_threads_fftbox, &
!$OMP      pub_threads_num_fftboxes,pub_eels_fine_projectors,loc_kpt)

          if (pub_eels_fine_projectors) then
             allocate(proj_complex(ld1_dbl, ld2_dbl, 2*n3), stat=ierr)
             call utils_alloc_check('projectors_func_pos_ovlp_box', &
                  'proj_complex',ierr)
          else
             allocate(proj_complex(ld1, ld2, n3), stat=ierr)
             call utils_alloc_check('projectors_func_pos_ovlp_box', &
                  'proj_complex',ierr)
          end if
          ! agrecocmplx: allocate required quantities for real/complex case
          ! real case
          if (.not.loc_cmplx) then
             allocate(fftbox1_dbl(ld1_dbl,ld2_dbl,2*n3), stat=ierr)
             call utils_alloc_check('projectors_func_pos_ovlp_box', &
                  'fftbox1_dbl',ierr)
             allocate(fftbox2_dbl(ld1_dbl,ld2_dbl,2*n3), stat=ierr)
             call utils_alloc_check('projectors_func_pos_ovlp_box', &
                  'fftbox2_dbl',ierr)
             allocate(new_proj(ld1, ld2, n3), stat=ierr)
             call utils_alloc_check('projectors_func_pos_ovlp_box', &
                  'new_proj',ierr)
             allocate(old_proj(ld1, ld2, n3), stat=ierr)
             call utils_alloc_check('projectors_func_pos_ovlp_box', &
                  'old_proj',ierr)
          ! complex case
          else
             allocate(fftbox_cmplx_dbl(ld1_dbl,ld2_dbl,2*n3),stat=ierr)
             call utils_alloc_check('projectors_func_pos_ovlp_box', &
                  'fftbox_cmplx_dbl',ierr)
             allocate(proj_cmplx_std(ld1, ld2, n3), stat=ierr)
             call utils_alloc_check('projectors_func_pos_ovlp_box', &
                  'proj_cmplx_std',ierr)
          end if

          allocate(coarse_work(ld1,ld2,n3),stat=ierr)
          call utils_alloc_check('projectors_func_pos_ovlp_box', &
               'coarse_work',ierr)
          allocate(fine_work(ld1_dbl,ld2_dbl,2*n3),stat=ierr)
          call utils_alloc_check('projectors_func_pos_ovlp_box', &
               'fine_work',ierr)

          ! ndmh: Set iteration range of each thread manually here - could use OMP do
          ! ndmh: but then range of iterations of each loop is hard to predict
          ! ndmh: between architectures. This way allows each thread to decide whether
          ! ndmh: to deposit its density box at the end of each iteration.
          thread_id = 0
          nthreads = 1
!$        nthreads = omp_get_num_threads()
!$        thread_id = omp_get_thread_num()
          thread_local_start = (thread_id * local_size) / nthreads &
               + local_start
          thread_local_end = ((thread_id + 1) * local_size) / nthreads &
               + local_start - 1

          count = 1

          do local_proj=thread_local_start,thread_local_end

             batch_proj = local_proj - local_start + 1

             ! Find information about this projector
             global_proj = local_proj + proj_basis%first_on_proc(pub_my_proc_id) &
                  - 1
             atom_of_proj = proj_basis%atom_of_func(global_proj)
             loc_atom_of_proj = atom_of_proj - &
                  proj_basis%atom_of_func(proj_basis%first_on_proc(pub_my_proc_id)) + 1
             isp = proj_set%proj_species(atom_of_proj)
             proj_count_atom = global_proj - &
                  proj_basis%first_on_atom(atom_of_proj) + 1
             proj_count = proj_set%species_first_proj(isp) + &
                  proj_count_atom - 1

             ! Centre of projector wrt fftbox in terms of grid spacings
             call basis_point_wrt_box(pcbg1, pcbg2, pcbg3, &  ! output
                  proj_set%proj_centre(atom_of_proj), n1, n2, n3, cell, &
                  fftbox)

             ! Start of fftbox wrt cell in terms of grid-point number
             call basis_start_of_box_wrt_cell(fftbox_start(1,batch_proj), &
                  fftbox_start(2,batch_proj),fftbox_start(3,batch_proj), &
                  proj_set%proj_centre(atom_of_proj),pcbg1,pcbg2,pcbg3,cell)

             if (pub_projectors_precalculate) then

                ! Apply phase factor e^-(ik.proj_centre_wrt_box) to projector
                ! stored in pre-calculated FFTbox
                ! ewt23
                if (pub_eels_fine_projectors) then
                   call basis_phase_on_fftbox_recip(proj_complex, &
                        2*n1,2*n2,2*n3,ld1_dbl,ld2_dbl,&
                        -2*pcbg1,-2*pcbg2,-2*pcbg3, &
                        proj_set%fftbox_proj_recip(:,:,:,proj_count))
                else
                   call basis_phase_on_fftbox_recip(proj_complex, &
                        n1,n2,n3,ld1,ld2,&
                        -pcbg1,-pcbg2,-pcbg3, &
                        proj_set%fftbox_proj_recip(:,:,:,proj_count))
                end if
             else

                ! Find angular momentum of this projector
                iproj = 1
                found = .false.
                do shell=1,proj_set%num_shells(isp)
                   am = proj_set%ang_mom(shell,isp)
                   do azim=-am,am
                      if (iproj==proj_count_atom) found=.true.
                      if (found) exit
                      iproj = iproj + 1
                   end do
                   if (found) exit
                end do

                if (pub_eels_fine_projectors) then
                     ! Generate projector in FFTbox
                     ! agrecokpt
                     call projector_in_box_recip(proj_complex, &
                          proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                          shell,isp), am, azim, 2*fftbox%total_pt1, &
                          2*fftbox%total_pt2, 2*fftbox%total_pt3, &
                          proj_set%n_rad_pts(isp), proj_set%gmax(isp),cell,&
                          kpt=loc_kpt,fine='Fine')

                     ! Apply phase factor e^-(ik.proj_centre_wrt_box)
                     ! to projector
                     call basis_phase_on_fftbox_recip(proj_complex, &
                          2*n1,2*n2,2*n3,ld1_dbl,ld2_dbl,&
                          -2*pcbg1,-2*pcbg2,-2*pcbg3)
                else
                     ! Generate projector in FFTbox
                     ! agrecokpt
                     call projector_in_box_recip(proj_complex, &
                          proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                          shell,isp), am, azim, fftbox%total_pt1, &
                          fftbox%total_pt2, fftbox%total_pt3, &
                          proj_set%n_rad_pts(isp), proj_set%gmax(isp), cell, &
                          kpt=loc_kpt)

                     ! Apply phase factor e^-(ik.proj_centre_wrt_box) to projector
                     call basis_phase_on_fftbox_recip(proj_complex, &
                          n1,n2,n3,ld1,ld2,-pcbg1,-pcbg2,-pcbg3)
                end if


             end if

             ! g=0 element must be real
             ! agrecokpt: but this is not true at a different kpt,
             ! because of the kshift applied? CHECK THIS!!!
             if (loc_kpt%x==0.0_DP.and.loc_kpt%y==0.0_DP.and.&
                 loc_kpt%z==0.0_DP) then
                proj_complex(1,1,1) &
                   = cmplx(real(proj_complex(1,1,1),kind=DP),0.0_DP,kind=DP)
             end if
             !proj_complex(1,1,1) &
             !   = cmplx(real(proj_complex(1,1,1),kind=DP),0.0_DP,kind=DP)

             ! Fourier transform to real space:
             if (pub_eels_fine_projectors) then
                call fourier_apply_box('Fine','Backward',proj_complex, &
                     omp=pub_threads_fftbox)
             else
                call fourier_apply_box('Coarse','Backward',proj_complex, &
                     omp=pub_threads_fftbox)
             end if

             ! agrecocmplx: distinguish between real and complex case
             ! in complex case, we must interpolate one at a time
             ! real case
             if (.not.loc_cmplx) then
                ! Now interpolate to fine grid and multiply by position operator
                ! Combine two projectors as real and imag parts of box, to avoid
                ! calling fourier_interpolate more than necessary.
                if (count == 2 .or. local_proj == thread_local_end) then
                   if (pub_eels_fine_projectors) then
                      fftbox2_dbl = real(proj_complex,dp)
                   else
                      new_proj = real(proj_complex,dp)
                   end if

                   ! interpolate fftboxes to fine grid
                   if (.not. pub_eels_fine_projectors) &
                        call fourier_interpolate(coarse_work,fine_work, &
                        old_proj,new_proj,fftbox1_dbl,fftbox2_dbl)

                   ! apply FFTbox position operator
                   call internal_apply_r(fftbox1_dbl,xyz,&
                        old_pcbg1,old_pcbg2,old_pcbg3)
                   call internal_apply_r(fftbox2_dbl,xyz,&
                        pcbg1,pcbg2,pcbg3)

                   ! filter fftboxes to standard grid
                   call fourier_filter(coarse_work,fine_work, &
                        fftbox1_dbl,fftbox2_dbl,old_proj,new_proj)

                else
                   !old_atom_of_proj = atom_of_proj
                   !old_local_proj = local_proj
                   if (pub_eels_fine_projectors) then
                        fftbox1_dbl = real(proj_complex,dp)
                   else
                        old_proj = real(proj_complex,dp)
                   end if
                   old_batch_proj = batch_proj
                   !old_r_fft = r_fft
                   old_pcbg1 = pcbg1
                   old_pcbg2 = pcbg2
                   old_pcbg3 = pcbg3
                endif

                ! Put real part into projector_box
                if (count == 2) then

                   !jmecmplx

                   projector_fftbox(old_batch_proj)%d(1:n1,1:n2,1:n3) &
                        = old_proj(1:n1,1:n2,1:n3)/cell%weight

                   projector_fftbox(batch_proj)%d(1:n1,1:n2,1:n3) &
                        = new_proj(1:n1,1:n2,1:n3)/cell%weight

                   count = 1
                else if (local_proj == thread_local_end) then

                   !jmecmplx
                   projector_fftbox(batch_proj)%d(1:n1,1:n2,1:n3) &
                        = new_proj(1:n1,1:n2,1:n3)/cell%weight
                else
                   count = count + 1
                end if
             ! agrecocmplx: complex case
             else
                ! interpolate 1 complex fftbox to fine grid
                ! use fully complex proj_complex, output is fftbox_cmplx_dbl
                if (.not. pub_eels_fine_projectors) then
                   call fourier_interpolate(coarse_work,fftbox_cmplx_dbl, &
                        proj_complex)
                ! projectors already on fine grid, simply copy to
                ! fftbox_cmplx_dbl
                else
                   fftbox_cmplx_dbl(:,:,:) = proj_complex(:,:,:)
                end if

                ! apply FFTbox position operator to complex box
                ! need to call the complex version of the routine
                call internal_apply_r_cmplx(fftbox_cmplx_dbl,xyz,&
                        pcbg1,pcbg2,pcbg3)

                ! filter fftboxes to standard grid
                call fourier_filter(coarse_work,fine_work, &
                     fftbox_cmplx_dbl,proj_cmplx_std)

                ! insert complex in FFTbox batch
                projector_fftbox(batch_proj)%z(1:n1,1:n2,1:n3) &
                      = proj_cmplx_std(1:n1,1:n2,1:n3)/cell%weight

             end if

          end do

          deallocate(fine_work,stat=ierr)
          call utils_dealloc_check('projectors_func_pos_ovlp_box', &
               'fine_work',ierr)
          deallocate(coarse_work,stat=ierr)
          call utils_dealloc_check('projectors_func_pos_ovlp_box', &
               'coarse_work',ierr)
          ! agrecocmplx: deallocate only allocated quantities
          ! real case
          if (.not.loc_cmplx) then
             deallocate(old_proj,stat=ierr)
             call utils_dealloc_check('projectors_func_pos_ovlp_box', &
                  'old_proj',ierr)
             deallocate(new_proj,stat=ierr)
             call utils_dealloc_check('projectors_func_pos_ovlp_box', &
                  'new_proj',ierr)
             deallocate(fftbox2_dbl,stat=ierr)
             call utils_dealloc_check('projectors_func_pos_ovlp_box', &
                  'fftbox2_dbl',ierr)
             deallocate(fftbox1_dbl,stat=ierr)
             call utils_dealloc_check('projectors_func_pos_ovlp_box', &
                  'fftbox1_dbl',ierr)
          ! complex case
          else
             deallocate(proj_cmplx_std,stat=ierr)
             call utils_dealloc_check('projectors_func_pos_ovlp_box', &
                  'proj_cmplx_std',ierr)
             deallocate(fftbox_cmplx_dbl,stat=ierr)
             call utils_dealloc_check('projectors_func_pos_ovlp_box', &
                  'fftbox_cmplx_dbl',ierr)
          end if
          deallocate(proj_complex,stat=ierr)
          call utils_dealloc_check('projectors_func_pos_ovlp_box', &
               'proj_complex',ierr)

!$OMP END PARALLEL

          ! Calculate overlap integrals
          ongpu=.false.
          call function_ops_brappd_ketfftbox(sp_overlap(xyz),  &  ! inout
               ngwfs_on_grid, ngwf_basis, cell, fftbox, &         ! input
               projector_fftbox, fftbox_start, batch_size,  &       ! input
               local_start, local_end, idx_len, sp_overlap_idx, &
               'FULL',ongpu)! input

          local_start = local_start + batch_size

       end do

    end do !xyz

    if (pub_projectors_precalculate) call projectors_exit_fftbox_recip(proj_set)

    ! Deallocate workspace
    deallocate(r_op,stat=ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box','r_op',ierr)
    deallocate(fftbox_start,stat=ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box', &
         'fftbox_start',ierr)
    do batch_count = 1, batch_size
        call data_fftbox_dealloc(projector_fftbox(batch_count))
    end do
    deallocate(projector_fftbox,stat=ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box', &
         'projector_fftbox',ierr)
    deallocate(sp_overlap_idx,stat=ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box', &
         'sp_overlap_idx',ierr)

    ! Stop timer
    call timer_clock('projectors_func_pos_ovlp_box',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &projectors_func_pos_ovlp_box'

contains



    subroutine internal_apply_r(fftbox_dbl,xyz,&
         pcbg1,pcbg2,pcbg3)

      real(kind=DP), intent(inout) :: fftbox_dbl(ld1_dbl,ld2_dbl,2*n3)
      integer, intent(in) :: xyz
      real(kind=DP), intent(in) :: pcbg1,pcbg2,pcbg3
      integer :: i1,i2,i3
      type(POINT) :: r1
      real(kind=DP) :: src_i1,src_i2,src_i3

      if (.not.(fftbox%coin1.or.fftbox%coin2.or.fftbox%coin3)) then
         ! If we're here then FFTbox=/=cell in any direction
         fftbox_dbl = r_op(xyz,:,:,:)*fftbox_dbl
      else
         do i3=1,2*n3
            do i2=1, ld2_dbl
               do i1=1, ld1_dbl
                  src_i1 = real(i1,kind=DP) - 2*pcbg1
                  src_i2 = real(i2,kind=DP) - 2*pcbg2
                  src_i3 = real(i3,kind=DP) - 2*pcbg3

                  if (src_i1 > ld1_dbl/2.0_DP)  src_i1 = src_i1 - ld1_dbl
                  if (src_i2 > ld2_dbl/2.0_DP)  src_i2 = src_i2 - ld2_dbl
                  if (src_i3 >  n3)             src_i3 = src_i3 - n3*2
                  if (src_i1 < -ld1_dbl/2.0_DP) src_i1 = src_i1 + ld1_dbl
                  if (src_i2 < -ld1_dbl/2.0_DP) src_i2 = src_i2 + ld2_dbl
                  if (src_i3 < -n3)             src_i3 = src_i3 + n3*2

                  r1 = (src_i1)*a1+&
                       (src_i2)*a2+&
                       (src_i3)*a3

                  if (xyz == 1) then
                     fftbox_dbl(i1,i2,i3) = fftbox_dbl(i1,i2,i3)*(r1%X)
                  else if (xyz==2) then
                     fftbox_dbl(i1,i2,i3) = fftbox_dbl(i1,i2,i3)*(r1%Y)
                  else
                     fftbox_dbl(i1,i2,i3) = fftbox_dbl(i1,i2,i3)*(r1%Z)
                  end if
               enddo
            enddo
         enddo
      end if
    end subroutine internal_apply_r

    ! complex version of internal_apply_r
    subroutine internal_apply_r_cmplx(fftbox_dbl,xyz,&
         pcbg1,pcbg2,pcbg3)

      complex(kind=DP), intent(inout) :: fftbox_dbl(ld1_dbl,ld2_dbl,2*n3)
      integer, intent(in) :: xyz
      real(kind=DP), intent(in) :: pcbg1,pcbg2,pcbg3
      integer :: i1,i2,i3
      type(POINT) :: r1
      real(kind=DP) :: src_i1,src_i2,src_i3

      if (.not.(fftbox%coin1.or.fftbox%coin2.or.fftbox%coin3)) then
         ! If we're here then FFTbox=/=cell in any direction
         fftbox_dbl = r_op(xyz,:,:,:)*fftbox_dbl
      else
         do i3=1,2*n3
            do i2=1, ld2_dbl
               do i1=1, ld1_dbl
                  src_i1 = real(i1,kind=DP) - 2*pcbg1
                  src_i2 = real(i2,kind=DP) - 2*pcbg2
                  src_i3 = real(i3,kind=DP) - 2*pcbg3

                  if (src_i1 > ld1_dbl/2.0_DP)  src_i1 = src_i1 - ld1_dbl
                  if (src_i2 > ld2_dbl/2.0_DP)  src_i2 = src_i2 - ld2_dbl
                  if (src_i3 >  n3)             src_i3 = src_i3 - n3*2
                  if (src_i1 < -ld1_dbl/2.0_DP) src_i1 = src_i1 + ld1_dbl
                  if (src_i2 < -ld1_dbl/2.0_DP) src_i2 = src_i2 + ld2_dbl
                  if (src_i3 < -n3)             src_i3 = src_i3 + n3*2

                  r1 = (src_i1)*a1+&
                       (src_i2)*a2+&
                       (src_i3)*a3

                  if (xyz == 1) then
                     fftbox_dbl(i1,i2,i3) = fftbox_dbl(i1,i2,i3)*(r1%X)
                  else if (xyz==2) then
                     fftbox_dbl(i1,i2,i3) = fftbox_dbl(i1,i2,i3)*(r1%Y)
                  else
                     fftbox_dbl(i1,i2,i3) = fftbox_dbl(i1,i2,i3)*(r1%Z)
                  end if
               enddo
            enddo
         enddo
      end if
    end subroutine internal_apply_r_cmplx

  end subroutine projectors_func_pos_ovlp_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine projectors_func_pos_ovlp_box_gradg(sp_overlap, &
       ngwfs_on_grid,ngwf_basis,proj_basis,proj_set,fftbox,cell, &
       kpt) ! agrecokpt

    !===============================================================!
    ! This subroutine calculates the overlap matrix between a set   !
    ! of NGWFs, and a set of projectors, with the position operator !
    ! in between the two, and stores the result in a SPAM3 matrix.  !
    ! It uses the derivative of the projector w.r.t. G in recip.    !
    ! space to compute the \vec{r}.projector product                !
    !---------------------------------------------------------------!
    ! Written by Edward Tait in Feb 2015, based on                  !
    ! the routine projectors_func_pos_ovlp_box by Laura Ratcliff.   !
    ! Modified by Andrea Greco in May 2016, to allow specification  !
    ! of k-point and compatibility with complex NGWFs.              !
    ! Simplified by Jose M Escartin in August 2017.                 !
    !===============================================================!

    use datatypes, only: FUNCTIONS
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_axpy, &
         sparse_scale, sparse_copy
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: sp_overlap(3)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(PROJECTOR_SET), intent(inout) :: proj_set
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(CELL_INFO), intent(in) :: cell
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    ! agrecokpt
    type(POINT), optional, intent(in) :: kpt

    ! Local variables
    integer :: xyz ! Cartesian direction
    real(kind=DP) :: delta
    character(len=*), parameter :: myself = 'projectors_func_pos_ovlp_box_gradg'
    ! agrecokpt
    type(POINT) :: loc_kpt

    delta = 0.1_DP

    call utils_assert(all(sp_overlap(:)%iscmplx), 'Error in '//myself//': &
         &the first argument must be an array of three complex SPAM3 matrices.')

    ! agrecokpt
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    ! agrecokpt: specification of kpt allowed only in complex case
    if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.&
        loc_kpt%z/=0.0_DP) then
       call utils_assert(ngwfs_on_grid%iscmplx, 'Error in '//myself//': &
            &non-Gamma k-point requires complex projectors.')
    end if

    ! From projectors_commutator
    ! Calculate <phi|Del_G(proj(G))> overlap matrix
    ! agrecocmplx: compute directly in complex case
    ! agrecokpt: at specified k-point
    do xyz=1,3
       call projectors_func_ovlp_box(sp_overlap(xyz), &
            ngwfs_on_grid, ngwf_basis, proj_basis, proj_set, fftbox, cell, &
            kshift=loc_kpt, delta=delta, cart=xyz)
    end do

    ! Multiply by i
    do xyz=1,3
       call sparse_scale(sp_overlap(xyz), (0.0_DP,1.0_DP) )
    end do

  end subroutine projectors_func_pos_ovlp_box_gradg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine projectors_gradient_batch(fftbox_batch, &
       mmat,nmat,local_start,local_end,batch_size,fa_box_start, &
       ngwf_basis,proj_basis,cell,fftbox,coeff_mat,ps_overlap,proj_set, &
       kpt, directions, weights) ! agrecokpt

    !====================================================================!
    ! This subroutine calculates the non-local potential contribution    !
    ! to the gradient of the energy with respect to the NGWF expansion   !
    ! coefficients in the psinc basis.                                   !
    !--------------------------------------------------------------------!
    ! Original version written by Arash A. Mostofi in January 2004.      !
    ! Modified for speed by Chris-Kriton Skylaris on 31/1/2004.          !
    ! Bug related to initialisation of overlap_test_fa located and       !
    ! fixed by Chris-Kriton Skylaris on 2/5/2004.                        !
    ! Modified to use parallel SPAM by Peter Haynes, July 2006.          !
    ! Modified to use batch system by Nicholas Hine, July 2008.          !
    ! Modified to use SPAM3 matrices by Nicholas Hine, July 2009.        !
    !--------------------------------------------------------------------!
    ! Re-written by Nicholas Hine in October 2009. Outer loop is now     !
    ! over NGWFs fa, inner loop split into two parts: first find overlaps!
    ! for each fa and and store the coefficients of RK and RSK for each  !
    ! NGWF-projector overlap, plus phase arrays e^iG.R for each overlap. !
    ! Next, loop over the overlaps found, and calculate a structure      !
    ! factor for each projector type, multiply this by the projector in  !
    ! reciprocal space and FFT back to real space to get the gradient.   !
    ! OpenMP parallelised by Nicholas Hine in May 2013.                  !
    ! Data parallel OpenMP code added by karl Wilkinson in July 2013     !
    !--------------------------------------------------------------------!
    ! First-order NGWF capabilities added by Gabriel Constantinescu in   !
    ! 2016                                                               !
    !====================================================================!

    use datatypes, only: FFTBOX_DATA
    use comms, only: pub_my_proc_id
    use fourier, only: fourier_apply_box
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_projectors_precalculate, pub_threads_fftbox, &
         pub_debug_on_root, pub_num_spins
!$  use rundat, only: pub_threads_num_fftboxes, pub_threads_per_fftbox
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index, &
         sparse_get_element, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: local_start,local_end,batch_size
    integer, intent(in) :: mmat,nmat
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FFTBOX_DATA), intent(inout) :: fftbox_batch(pub_num_spins, &
         nmat, batch_size)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    integer, intent(in) :: fa_box_start(3,batch_size)
    type(SPAM3), intent(in) :: coeff_mat(pub_num_spins,nmat)
    type(SPAM3), intent(in) :: ps_overlap
    type(PROJECTOR_SET), intent(inout) :: proj_set
    integer, intent(in), optional :: directions(:)!par%nat)
    real(kind=DP), intent(in), optional :: weights(:)!par%nat)
    ! agrecokpt: optional argument to specify k-point
    type(POINT), optional, intent(in) :: kpt

    ! Local Variables
    integer :: n1,n2,n3,ld1,ld2
    integer :: fa, local_fa
    integer :: loc_fa_atom, loc_fa_atom_prev
    integer :: ovlp_atom
    integer :: isp
    integer :: jproj, shell, am, azim
    integer :: proj_count
    integer :: global_proj
    integer :: ierr
    integer :: idx_len
    integer :: fa_atom_overlaps
    integer :: ovlp_list_len
    integer :: iproj, iovlp, novlp
    integer :: i2,i3
    integer :: is,idx,imat
    integer :: batch_count
    logical :: found
    complex(kind=DP) :: phase23
    complex(kind=DP), allocatable :: factors(:)
    integer, allocatable :: ps_overlap_idx(:)
    integer, allocatable :: ovlp_species(:), ovlp_index(:)
    real(kind=DP), allocatable :: ovlp_proj_coeffs_real(:,:,:)
    ! agrecocmplx: complex version of ovlp_proj_coeffs
    complex(kind=DP), allocatable :: ovlp_proj_coeffs_cmplx(:,:,:)
    complex(kind=DP), allocatable :: out_complex(:,:,:)
    complex(kind=DP), allocatable :: proj_complex(:,:,:)
    complex(kind=DP), allocatable :: struct_fac(:,:,:)
    complex(kind=DP), allocatable :: phases1(:,:),phases2(:,:),phases3(:,:)
    ! agrecocmplx
    logical :: loc_cmplx, loc_first_order
    ! agrecokpt
    type(POINT) :: loc_kpt
    ! rc2013: parallel strategy info
    type(PARAL_INFO), pointer :: par

    ! -------------------------------------------------------------------------

    ! cks: time it
    call timer_clock('projectors_gradient_batch', 1)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering projectors_gradient_batch'

    ! agrecocmplx
    loc_cmplx = ps_overlap%iscmplx

    ! agrecokpt
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    if (present(directions).and.present(weights)) then
       loc_first_order = .true.
    else
       loc_first_order = .false.
    end if

    ! rc2013: obtain parallel strategy from sparse matrix
    call sparse_get_par(par, ps_overlap, 'R')

    ! agrecocmplx: should be compatible now, need to check everything
    ! is working properly
    !jmecmplx
    !call utils_assert(.not. loc_cmplx, 'Error in&
    !     & projectors_gradient_batch: not ready yet for complex NGWFs.')

    ! Initialisations
    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2

    ! Obtain index of ps_overlap
    idx_len = sparse_index_length(ps_overlap)
    allocate(ps_overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','ps_overlap_idx',ierr)
    call sparse_generate_index(ps_overlap_idx,ps_overlap)

    ! Find maximum number of projectors with which any fa in batch overlaps
    ovlp_list_len = 0
    do local_fa=local_start,local_end
       fa = local_fa + ngwf_basis%first_on_proc(pub_my_proc_id) - 1
       loc_fa_atom = ngwf_basis%atom_of_func(fa) - &
            ngwf_basis%atom_of_func(ngwf_basis%first_on_proc(pub_my_proc_id)) + 1
       fa_atom_overlaps = ps_overlap_idx(loc_fa_atom+1) - &
            ps_overlap_idx(loc_fa_atom) + 1
       ovlp_list_len = max(ovlp_list_len,fa_atom_overlaps)
    end do

    ! agrecokpt: initialise projectors for specified k-point
    if (pub_projectors_precalculate) &
         call projectors_init_fftbox_recip(proj_set,cell,fftbox,kpt=loc_kpt)

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(local_fa,batch_count,fa,loc_fa_atom,novlp,factors, &
!$OMP      isp,proj_count,out_complex,struct_fac,ovlp_species, ovlp_index, &
!$OMP      ovlp_proj_coeffs_real,ovlp_proj_coeffs_cmplx,is,imat,phase23,phases1, &
!$OMP      phases2,phases3,proj_complex,ierr,iproj,found,am,azim,ovlp_atom,jproj, &
!$OMP      global_proj,loc_fa_atom_prev,i2,i3,idx,shell,iovlp) &
!$OMP SHARED(local_start,local_end,proj_set,ld1,ld2,n1,n2,n3,batch_size, &
!$OMP      proj_basis,pub_my_proc_id,cell,fftbox, &
!$OMP      pub_projectors_precalculate,ovlp_list_len,mmat,nmat, &
!$OMP      ngwf_basis,ps_overlap_idx,coeff_mat,fa_box_start, directions, &
!$OMP      weights,fftbox_batch,pub_threads_fftbox,pub_num_spins, &
!$OMP      pub_threads_num_fftboxes,pub_threads_per_fftbox, &
!$OMP      loc_cmplx,loc_kpt,loc_first_order, par)

    ! Memory allocation
    allocate(out_complex(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','out_complex',ierr)
    allocate(proj_complex(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','proj_complex',ierr)
    allocate(struct_fac(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','struct_fac',ierr)
    if (loc_first_order) then
       allocate(ovlp_index(ovlp_list_len),stat=ierr)
       call utils_alloc_check('projectors_gradient_batch','ovlp_index',ierr)
    end if
    allocate(factors(n1), stat=ierr)
    call utils_alloc_check('projectors_gradient_batch', 'factors',ierr)
    allocate(ovlp_species(ovlp_list_len),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','ovlp_species',ierr)
    ! agrecocmplx
    if (loc_cmplx) then
       allocate(ovlp_proj_coeffs_cmplx(ovlp_list_len,proj_basis%max_on_atom, &
                mmat:nmat), stat=ierr)
       call utils_alloc_check('projectors_gradient_batch', &
            'ovlp_proj_coeffs_cmplx',ierr)
    else
       allocate(ovlp_proj_coeffs_real(ovlp_list_len,proj_basis%max_on_atom, &
                mmat:nmat), stat=ierr)
       call utils_alloc_check('projectors_gradient_batch', &
            'ovlp_proj_coeffs_real',ierr)
    end if
    allocate(phases1(1:n1,ovlp_list_len),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','phases1',ierr)
    allocate(phases2(1:n2,ovlp_list_len),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','phases2',ierr)
    allocate(phases3(1:n3,ovlp_list_len),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','phases3',ierr)

    loc_fa_atom_prev = -9999

    ! Loop over spins
    do is=1,pub_num_spins

       ! Loop over NGWFs in batch
!$OMP DO
       do local_fa=local_start,local_end

          batch_count = local_fa - local_start + 1
          fa = local_fa + ngwf_basis%first_on_proc(pub_my_proc_id) - 1
          loc_fa_atom = ngwf_basis%atom_of_func(fa) - &
               ngwf_basis%atom_of_func(ngwf_basis%first_on_proc(pub_my_proc_id)) + 1

          ! Reset number of overlaps and projectors for this fa
          novlp = 0

          ! Loop over atoms with projectors overlapping fa
          do idx=ps_overlap_idx(loc_fa_atom),ps_overlap_idx(loc_fa_atom+1)-1
             ovlp_atom = ps_overlap_idx(idx)
             novlp = novlp + 1

             ! If fa is on a new atom, we need to update the entries in the
             ! phases arrays for this overlap
             if (loc_fa_atom/=loc_fa_atom_prev) &
                  call internal_phases(phases1,phases2,phases3, &
                  ovlp_atom,novlp,local_fa,fa_box_start(1,batch_count), &
                  fa_box_start(2,batch_count),fa_box_start(3,batch_count), &
                  proj_set,ngwf_basis, &
                  cell,fftbox)

             ! Store information about the species of this atom
             ovlp_species(novlp) = proj_set%proj_species(ovlp_atom)
             if (loc_first_order) ovlp_index(novlp) = ovlp_atom

             ! Loop over all the projectors on this atom storing coefficients
             do iproj=1,proj_basis%num_on_atom(ovlp_atom)
                global_proj = proj_basis%first_on_atom(ovlp_atom)+iproj-1

                ! Store \sum_bj 4*Kab*D_ij*<pj|fb> in ovlp_proj_coeffs(:,i,3)
                ! and \sum_bj 4*Kac*Scb*D_ij*<pj|fb> in ovlp_proj_coeffs(:,i,4)
                ! agrecocmplx
                ! complex case
                if (loc_cmplx) then
                   do imat=mmat,nmat
                      call sparse_get_element(ovlp_proj_coeffs_cmplx(novlp,iproj,imat), &
                           coeff_mat(is,imat),global_proj,fa)
                   end do
                ! real case
                else
                   do imat=mmat,nmat
                      call sparse_get_element(ovlp_proj_coeffs_real(novlp,iproj,imat), &
                           coeff_mat(is,imat),global_proj,fa)
                   end do
                end if
             end do

          end do

          ! Loop over the types of gradient fftbox we are accumulating
          do imat=mmat,nmat ! 3 is ham_fftbox_batch, 4 is tc_ham_fftbox_batch

             ! Zero reciprocal space nonlocal gradient box
             out_complex(:,:,:) = cmplx(0.0_DP,0.0_DP,kind=DP)

             ! Loop over atom species
             do isp=1,proj_set%n_proj_species

                ! Cycle if fa does not overlap any atoms of this species
                if (.not.any(ovlp_species(1:novlp)==isp)) cycle

                ! Loop over all the projector types for this species
                do iproj=1,proj_set%species_num_proj(isp)

                   ! Find entry in global projector fftboxes array
                   proj_count = proj_set%species_first_proj(isp) + &
                        iproj - 1

                   if (.not.pub_projectors_precalculate) then

                      ! Find angular momentum of this projector
                      jproj = 1
                      found = .false.
                      do shell=1,proj_set%num_shells(isp)
                         am = proj_set%ang_mom(shell,isp)
                         do azim=-am,am
                            if (jproj==iproj) found=.true.
                            if (found) exit
                            jproj = jproj + 1
                         end do
                         if (found) exit
                      end do

                      ! Generate projector in FFTbox
                      ! agrecokpt: call with appropriate k-point
                      call projector_in_box_recip(proj_complex, &
                           proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                           shell,isp), am, azim, fftbox%total_pt1, &
                           fftbox%total_pt2, fftbox%total_pt3, &
                           proj_set%n_rad_pts(isp), proj_set%gmax(isp), cell, &
                           kpt=loc_kpt)

                   end if

                   ! Reset structure factor to zero for new projector type
                   struct_fac(:,:,:) = cmplx(0.0_DP,0.0_DP,kind=DP)

                   ! Loop over slabs of reciprocal space points of fftbox
!!$OMP PARALLEL DO NUM_THREADS(pub_threads_per_fftbox) DEFAULT(none) &
!!$OMP PRIVATE(i3,iovlp,phase23, factors) &
!!$OMP SHARED(n1,n2,n3,novlp,isp,ovlp_species,phases1,phases2,phases3, &
!!$OMP      struct_fac,ovlp_proj_coeffs_real,ovlp_proj_coeffs_cmplx, &
!!$OMP      pub_projectors_precalculate,out_complex,proj_set,proj_count, &
!!$OMP      iproj,imat,proj_complex,pub_threads_per_fftbox,loc_cmplx, &
!!$OMP      directions, weights, loc_first_order)
                   do i3=1,n3

                      ! Loop over atoms with projectors overlapping fa
                      do iovlp=1,novlp

                         ! Cycle if overlapping atom is not of the right species
                         if (.not.ovlp_species(iovlp)==isp) cycle

                         ! Loop over recip points in this slab of fftbox to
                         ! calculate structure factor for this projector type
                         ! agrecocmplx
                         ! complex case
                         if (loc_cmplx) then
                            do i2=1,n2
                               phase23 = phases2(i2,iovlp) * phases3(i3,iovlp)
                               ! gcc32: first-order-case
                               if (loc_first_order.and.(imat == 2)) then
                                  factors(1:n1) = cmplx(0.0_DP, &
                                       -fftbox%recip_grid(directions( &
                                       par%orig_atom(ovlp_index(iovlp))), &
                                       1:n1,i2,i3) * weights(par%orig_atom( &
                                       ovlp_index(iovlp))),kind=DP)
                               else
                                  factors(1:n1) = cmplx(1.0_DP,0.0_DP,kind=DP)
                               end if
                               struct_fac(1:n1,i2,i3) = struct_fac(1:n1,i2,i3)+&
                                    phases1(1:n1,iovlp) * phase23 * &
                                    ovlp_proj_coeffs_cmplx(iovlp,iproj,imat) * &
                                    factors(1:n1)
                            end do
                         ! real case
                         else
                            do i2=1,n2
                               phase23 = phases2(i2,iovlp) * phases3(i3,iovlp)
                               if (loc_first_order.and.(imat == 2)) then
                                  factors(1:n1) = cmplx(0.0_DP, &
                                       -fftbox%recip_grid(directions( &
                                       par%orig_atom(ovlp_index(iovlp))), &
                                       1:n1,i2,i3) * weights(par%orig_atom( &
                                       ovlp_index(iovlp))),kind=DP)
                               else
                                  factors(1:n1) = cmplx(1.0_DP,0.0_DP,kind=DP)
                               end if
                               struct_fac(1:n1,i2,i3) = struct_fac(1:n1,i2,i3)+&
                                    phases1(1:n1,iovlp) * phase23 * &
                                    ovlp_proj_coeffs_real(iovlp,iproj,imat) * &
                                    factors(1:n1)
                            end do
                         end if

                      end do

                      ! Multiply the structure factor for this projector type by
                      ! the projector in reciprocal space, and add it to the
                      ! current reciprocal space gradient
                      if (pub_projectors_precalculate) then
                         out_complex(:,:,i3) = out_complex(:,:,i3) + &
                              struct_fac(:,:,i3) * &
                              proj_set%fftbox_proj_recip(:,:,i3,proj_count)
                      else
                         out_complex(:,:,i3) = out_complex(:,:,i3) + &
                              struct_fac(:,:,i3) * proj_complex(:,:,i3)
                      end if

                   end do
!!$OMP END PARALLEL DO

                end do  ! iproj

             end do  ! species_number

             ! Fourier transform gradient contribution to real space:
             call fourier_apply_box('Coarse', 'Backward', out_complex(:,:,:), &
                  omp=pub_threads_fftbox)

             ! agrecocmplx
             if (loc_cmplx) then
                ! Put complex into relevant gradient fftbox
!!$OMP PARALLEL WORKSHARE NUM_THREADS(pub_threads_per_fftbox) &
!!$OMP SHARED(pub_threads_per_fftbox)
                fftbox_batch(is,imat,batch_count)%z(1:n1,1:n2,1:n3) = &
                     fftbox_batch(is,imat,batch_count)%z(1:n1,1:n2,1:n3) + &
                     out_complex(1:n1,1:n2,1:n3)
!!$OMP END PARALLEL WORKSHARE
             else
                ! Put real part of this into relevant gradient fftbox
!!$OMP PARALLEL WORKSHARE NUM_THREADS(pub_threads_per_fftbox) &
!!$OMP SHARED(pub_threads_per_fftbox)
                fftbox_batch(is,imat,batch_count)%d(1:n1,1:n2,1:n3) = &
                     fftbox_batch(is,imat,batch_count)%d(1:n1,1:n2,1:n3) + &
                     real(out_complex(1:n1,1:n2,1:n3),kind=DP)
!!$OMP END PARALLEL WORKSHARE
             end if
          end do

          ! Store previous local index of fa atom, so that we can tell if
          ! phases arrays need regenerating
          loc_fa_atom_prev = loc_fa_atom

       end do ! loop over NGWFs in batch
!$OMP END DO

    end do ! loop over spins

    ! Deallocate Memory
    deallocate(phases3,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','phases3', &
         ierr)
    deallocate(phases2,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','phases2', &
         ierr)
    deallocate(phases1,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','phases1', &
         ierr)
    if (loc_first_order) then
       deallocate(ovlp_index,stat=ierr)
       call utils_dealloc_check('projectors_gradient_batch','ovlp_index',ierr)
    end if
    deallocate(factors, stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch', 'factors',ierr)
    ! agrecocmplx
    ! complex case
    if (loc_cmplx) then
       deallocate(ovlp_proj_coeffs_cmplx,stat=ierr)
       call utils_dealloc_check('projectors_gradient_batch', &
            'ovlp_proj_coeffs_cmplx',ierr)
    ! real case
    else
       deallocate(ovlp_proj_coeffs_real,stat=ierr)
       call utils_dealloc_check('projectors_gradient_batch', &
            'ovlp_proj_coeffs_real',ierr)
    end if
    deallocate(ovlp_species,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','ovlp_species', &
         ierr)
    deallocate(struct_fac,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','struct_fac',&
         ierr)
    deallocate(proj_complex,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','proj_complex',&
         ierr)
    deallocate(out_complex,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','out_complex',&
         ierr)
!$OMP END PARALLEL

    if (pub_projectors_precalculate) call projectors_exit_fftbox_recip(proj_set)

    deallocate(ps_overlap_idx,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','ps_overlap_idx',&
         ierr)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving projectors_gradient_batch'

    ! cks: time it
    call timer_clock('projectors_gradient_batch', 2)

  end subroutine projectors_gradient_batch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine internal_phases(phases1,phases2,phases3,ovlp_atom,novlp, &
       local_fa,fa_box_start1,fa_box_start2, &
       fa_box_start3,proj_set,ngwf_basis,cell,fftbox)

    !=======================================================================!
    ! This subroutine calculates the phase shifts for each reciprocal-space !
    ! grid point in the FFTbox centred on NGWF fa, given a projector at     !
    ! position proj_centre.                                                 !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine in October 2009 based on bits of             !
    ! basis_phase_on_fftbox_recip.                                          !
    ! Made standalone for thread-safety by Nicholas Hine in May 2013.       !
    !=======================================================================!

    use constants, only: two_pi, cmplx_i
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, geometry_magnitude, OPERATOR(.dot.), &
         OPERATOR(+), OPERATOR(*)
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: phases1(:,:),phases2(:,:),phases3(:,:)
    integer, intent(in) :: ovlp_atom,novlp,local_fa
    integer, intent(in) :: fa_box_start1,fa_box_start2,fa_box_start3
    type(PROJECTOR_SET), intent(in) :: proj_set
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local Variables
    real(kind=DP),parameter :: inv_two_pi = 0.5_DP / PI
    real(kind=DP) :: a1,a2,a3
    integer :: i1,i2,i3
    real(kind=DP) :: box_start1,box_start2,box_start3
    real(kind=DP) :: proj_radius
    complex(kind=DP) :: phase_inc, phase_neg
    type(POINT) :: proj_centre
    integer :: n1,n2,n3

    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3

    proj_centre = proj_set%proj_centre(ovlp_atom)
    proj_radius = proj_set%proj_max_radius(ovlp_atom)

    ! Start of fa fftbox wrt sim cell in terms of numbers of grid points
    call internal_origin_nl_fa(box_start1,box_start2, &
         box_start3,fa_box_start1,fa_box_start2,fa_box_start3,proj_radius, &
         proj_centre,ngwf_basis%spheres(local_fa),cell)

    ! Calculate (box_start - proj_centre) in terms of number of grid
    ! points in each lattice direction.
    a1 = (proj_centre .DOT. (cell%b1))
    box_start1 = box_start1 - a1 * inv_two_pi * &
         real(cell%total_pt1,kind=DP)

    a2 = (proj_centre .DOT. (cell%b2))
    box_start2 = box_start2 - a2 * inv_two_pi * &
         real(cell%total_pt2,kind=DP)

    a3 = (proj_centre .DOT. (cell%b3))
    box_start3 = box_start3 - a3 * inv_two_pi * &
         real(cell%total_pt3,kind=DP)

    ! Phases e^(iG1.R) along '1' direction
    phase_inc = exp(two_pi*cmplx_i/real(n1,kind=DP) * box_start1)
    phase_neg = exp(-two_pi*cmplx_i * box_start1)
    phases1(1,novlp) = cmplx(1.0_DP,0.0_DP,kind=DP)
    do i1=2,n1/2+2
       phases1(i1,novlp) = phases1(i1-1,novlp)*phase_inc
    end do
    i1=n1/2+2; phases1(i1,novlp) = phases1(i1,novlp)*phase_neg
    do i1=n1/2+3,n1
       phases1(i1,novlp) = phases1(i1-1,novlp)*phase_inc
    end do

    ! Phases e^(iG2.R) along '2' direction
    phase_inc = exp(two_pi*cmplx_i/real(n2,kind=DP) * box_start2)
    phase_neg = exp(-two_pi*cmplx_i * box_start2)
    phases2(1,novlp) = cmplx(1.0_DP,0.0_DP,kind=DP)
    do i2=2,n2/2+2
       phases2(i2,novlp) = phases2(i2-1,novlp)*phase_inc
    end do
    i2=n2/2+2; phases2(i2,novlp) = phases2(i2,novlp)*phase_neg
    do i2=n2/2+3,n2
       phases2(i2,novlp) = phases2(i2-1,novlp)*phase_inc
    end do

    ! Phases e^(iG3.R) along '3' direction
    phase_inc = exp(two_pi*cmplx_i/real(n3,kind=DP) * box_start3)
    phase_neg = exp(-two_pi*cmplx_i * box_start3)
    phases3(1,novlp) = cmplx(1.0_DP,0.0_DP,kind=DP)
    do i3=2,n3/2+2
       phases3(i3,novlp) = phases3(i3-1,novlp)*phase_inc
    end do
    i3=n3/2+2; phases3(i3,novlp) = phases3(i3,novlp)*phase_neg
    do i3=n3/2+3,n3
       phases3(i3,novlp) = phases3(i3-1,novlp)*phase_inc
    end do

  end subroutine internal_phases


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine internal_origin_nl_fa( &
       box_start1, box_start2, box_start3, &               ! output
       fa_box_start1,fa_box_start2,fa_box_start3, &        ! input
       proj_radius, proj_centre, ngwf_sphere, cell)        ! input

    !============================================================!
    ! Returns the distance, in terms of gridpoints (real values) !
    ! in each lattice direction, from the simulation cell origin !
    ! to the nearest vertex of the fftbox.                       !
    !------------------------------------------------------------!
    ! Written by Arash A. Mostofi, August 2002, using some parts !
    ! from older subroutines written by Chris-Kriton Skylaris.   !
    ! Modified by Arash A. Mostofi on July 2003 and on           !
    ! January 2004.                                              !
    ! Made standalone for thread-safety by Nicholas Hine in May  !
    ! 2013.                                                      !
    !============================================================!

    use basis, only: SPHERE
    use geometry, only: POINT, geometry_magnitude, &
         OPERATOR(.dot.), OPERATOR(+), OPERATOR(*)
    use simulation_cell, only: CELL_INFO

    implicit none

    type(SPHERE), intent(in)           :: ngwf_sphere
    type(POINT), intent(in)            :: proj_centre
    integer, intent(in) :: fa_box_start1,fa_box_start2,fa_box_start3
    real(kind=DP), intent(in)          :: proj_radius
    real(kind=DP), intent(out)         :: box_start1,box_start2,box_start3
    type(CELL_INFO), intent(in)        :: cell

    real(kind=DP),parameter :: inv_two_pi = 0.5_DP / PI

    ! aam: internal declarations
    real(kind=DP) :: sum_radii
    real(kind=DP) :: dist1,dist2,dist3
    type(POINT)   :: ngwf_proj

    box_start1 = real(fa_box_start1-1,kind=DP)
    box_start2 = real(fa_box_start2-1,kind=DP)
    box_start3 = real(fa_box_start3-1,kind=DP)

    ! cks: sum of projector and nwgf radii
    sum_radii = proj_radius + ngwf_sphere%radius

    ! cks: vector from centre of projector to centre of ngwf
    ngwf_proj = ngwf_sphere%centre+(-1.0_DP)*proj_centre

    ! cks: a1_unit, a2_unit and a3_unit components of ngwf_proj
    dist1 = (ngwf_proj.DOT.cell%b1) * &
         inv_two_pi*geometry_magnitude(cell%a1)
    dist2 = (ngwf_proj.DOT.cell%b2) * &
         inv_two_pi*geometry_magnitude(cell%a2)
    dist3 = (ngwf_proj.DOT.cell%b3) * &
         inv_two_pi*geometry_magnitude(cell%a3)

    ! cks: if necessary, move the origin of the ngwf-fftbox by a lattice
    ! cks: vector in each direction. This need arises when the projector
    ! cks: and the ngwf overlap not inside the simulation cell but through
    ! cks: its faces due to the periodic boundary conditions. Then, the
    ! cks: displacement of the ngwf-fftbox places it in the position of
    ! cks: one of its periodic images so that the origin of the projector
    ! cks: w.r.t the origin of the box will have the correct relationship
    ! cks: with the origin of the ngwf w.r.t. the box.
    if ( abs(dist1) > sum_radii ) box_start1 = box_start1 &
         - SIGN(1.0_DP,dist1)*real(cell%total_pt1,kind=DP)

    if ( abs(dist2) > sum_radii ) box_start2 = box_start2 &
         - SIGN(1.0_DP,dist2)*real(cell%total_pt2,kind=DP)

    if ( abs(dist3) > sum_radii ) box_start3 = box_start3 &
         - SIGN(1.0_DP,dist3)*real(cell%total_pt3,kind=DP)

  end subroutine internal_origin_nl_fa


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_grad_precond_batch(fftbox_batch,nmat, &
       local_start,local_end,batch_size,fa_box_start, &
       cell,fftbox,ngwf_basis,ps_overlap,proj_set,precond_func_recip, &
       kpt) ! agrecokpt

    !====================================================================!
    ! This subroutine calculates the PAW-correction to the kinetic       !
    ! energy preconditining of the gradient of the energy with respect   !
    ! to the NGWF expansion coefficients in the psinc basis.             !
    !--------------------------------------------------------------------!
    ! Original version written by Gabriel Constantinescu in May 2014.    !
    ! Constains reverse-engineered parts of projectors_gradient_batch and!
    ! ngwf_gradient_batch                                                !
    !--------------------------------------------------------------------!
    !====================================================================!

    use datatypes, only: FFTBOX_DATA
    use comms, only: pub_my_proc_id
    use constants, only: DP, stdout
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_projectors_precalculate, &
          pub_threads_fftbox, pub_debug_on_root, pub_num_spins
!$  use rundat, only: pub_threads_num_fftboxes, pub_threads_per_fftbox
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: nmat
    integer, intent(in) :: local_start,local_end,batch_size
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(CELL_INFO), intent(in) :: cell
    integer, intent(in) :: fa_box_start(3,batch_size)
    type(SPAM3), intent(in) :: ps_overlap
    type(PROJECTOR_SET), intent(inout) :: proj_set
    type(FFTBOX_DATA), intent(inout) :: fftbox_batch(pub_num_spins, nmat, &
                                                     batch_size)
    real(kind=DP), intent(in)  :: precond_func_recip(:,:,:)
    ! agrecokpt: optional k-point argument
    type(POINT), optional, intent(in) :: kpt

    ! Local Variables
    integer :: n1,n2,n3,ld1,ld2
    integer :: fa, local_fa
    integer :: loc_fa_atom, loc_fa_atom_prev
    integer :: ovlp_atom
    integer :: isp
    integer :: jproj, shell, am, azim
    integer :: proj_count, proj_count_1
    integer :: ierr
    integer :: idx_len
    integer :: fa_atom_overlaps
    integer :: ovlp_list_len
    integer :: iproj, iproj_1, iovlp, novlp
    integer :: i1,i2,i3
    integer :: is,idx
    integer :: batch_count
    logical :: found
    integer :: max_paw_proj_tot
    complex(kind=DP) :: phase23

    real(kind=DP) :: normalization


    integer, allocatable :: ps_overlap_idx(:)
    integer, allocatable :: ovlp_species(:)

    integer, allocatable :: ovlp_atoms(:)

    complex(kind=DP), allocatable :: out_complex(:,:,:)
    complex(kind=DP), allocatable :: proj_complex(:,:,:)
    complex(kind=DP), allocatable :: proj_complex_1(:,:,:)
    complex(kind=DP), allocatable :: struct_fac(:,:,:)
    complex(kind=DP), allocatable :: coarse_work(:,:,:)

    real(kind=DP), allocatable :: A_mat_real(:,:)
    ! agrecocmplx: complex version for complex case
    complex(kind=DP), allocatable :: A_mat_cmplx(:,:)
    complex(kind=DP), allocatable ::  B_mat(:,:)
    complex(kind=DP), allocatable :: phases1(:,:),phases2(:,:),phases3(:,:)
    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt
    type(POINT) :: loc_kpt

    ! time it
    call timer_clock('projectors_grad_precond_batch', 1)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering projectors_grad_precond_batch'

    ! agrecocmplx
    loc_cmplx = ps_overlap%iscmplx

    ! agrecokpt
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    ! agrecocmplx: should be ok now, need to check what happens to O_prec
    ! in the case of k-points/complex NGWFs
    !jmecmplx
    !call utils_assert(.not. loc_cmplx, 'Error in &
    !     & projectors_grad_precond_batch: not ready yet for complex NGWFs.')

    ! Initialisations
    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2

    normalization = 1.0_DP / (n1*n2*n3)

    max_paw_proj_tot = maxval(proj_set%species_num_proj(:))

    ! Obtain index of ps_overlap
    idx_len = sparse_index_length(ps_overlap)
    allocate(ps_overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('projectors_grad_precond_batch','ps_overlap_idx',ierr)
    call sparse_generate_index(ps_overlap_idx,ps_overlap)


    ! Find maximum number of projectors with which any fa in batch overlaps
    ovlp_list_len = 0
    do local_fa=local_start,local_end
       fa = local_fa + ngwf_basis%first_on_proc(pub_my_proc_id) - 1
       loc_fa_atom = ngwf_basis%atom_of_func(fa) - &
            ngwf_basis%atom_of_func(ngwf_basis%first_on_proc(pub_my_proc_id)) + 1
       fa_atom_overlaps = ps_overlap_idx(loc_fa_atom+1) - &
            ps_overlap_idx(loc_fa_atom) + 1
       ovlp_list_len = max(ovlp_list_len,fa_atom_overlaps)
    end do

    ! agrecokpt: call with specified k-point
    if (pub_projectors_precalculate) then
        call projectors_init_fftbox_recip(proj_set,cell,fftbox,kpt=loc_kpt)
    end if


!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(local_fa,batch_count,fa,loc_fa_atom,novlp, &
!$OMP      isp,proj_count,proj_count_1,out_complex,struct_fac,ovlp_species, &
!$OMP      ovlp_atoms,is,phase23,phases1,phases2,phases3,proj_complex, &
!$OMP      proj_complex_1,ierr,iproj,iproj_1,found,am,azim,ovlp_atom,jproj, &
!$OMP      loc_fa_atom_prev,i1,i2,i3,idx,shell,iovlp,coarse_work,A_mat_real, &
!$OMP      A_mat_cmplx,B_mat) &
!$OMP SHARED(local_start,local_end,proj_set,ld1,ld2,n1,n2,n3,batch_size, &
!$OMP      pub_my_proc_id,cell,fftbox,normalization, &
!$OMP      max_paw_proj_tot,pub_projectors_precalculate,ovlp_list_len, nmat, &
!$OMP      ngwf_basis,ps_overlap_idx,fa_box_start, &
!$OMP      fftbox_batch,pub_threads_fftbox,precond_func_recip, &
!$OMP      pub_num_spins,pub_threads_num_fftboxes,pub_threads_per_fftbox, &
!$OMP      loc_cmplx, loc_kpt)


    ! Memory allocation
    allocate(out_complex(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('projectors_grad_precond_batch','out_complex',ierr)

    allocate(proj_complex(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('projectors_grad_precond_batch','proj_complex',ierr)

    allocate(proj_complex_1(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('projectors_grad_precond_batch','proj_complex_1',ierr)

    allocate(struct_fac(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('projectors_grad_precond_batch','struct_fac',ierr)

    allocate(ovlp_species(ovlp_list_len),stat=ierr)
    call utils_alloc_check('projectors_grad_precond_batch','ovlp_species',ierr)

    allocate(ovlp_atoms(ovlp_list_len),stat=ierr)
    call utils_alloc_check('projectors_grad_precond_batch','ovlp_atoms',ierr)

    allocate(phases1(1:n1,ovlp_list_len),stat=ierr)
    call utils_alloc_check('projectors_grad_precond_batch','phases1',ierr)

    allocate(phases2(1:n2,ovlp_list_len),stat=ierr)
    call utils_alloc_check('projectors_grad_precond_batch','phases2',ierr)

    allocate(phases3(1:n3,ovlp_list_len),stat=ierr)
    call utils_alloc_check('projectors_grad_precond_batch','phases3',ierr)

    allocate(coarse_work(fftbox%total_ld1,fftbox%total_ld2,&
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('projectors_grad_precond_batch','coarse_work',ierr)


    loc_fa_atom_prev = -9999

    ! Loop over spins
    do is=1,pub_num_spins
       ! Loop over NGWFs in batch
!$OMP DO
       do local_fa=local_start,local_end

          batch_count = local_fa - local_start + 1
          fa = local_fa + ngwf_basis%first_on_proc(pub_my_proc_id) - 1
          loc_fa_atom = ngwf_basis%atom_of_func(fa) - &
               ngwf_basis%atom_of_func(ngwf_basis%first_on_proc(pub_my_proc_id)) + 1

          ! You are taking as input the KE preconditioned (as for NCPPs)
          ! covariant gradient in real space (fftbox_batch)

          coarse_work = cmplx(0.0_DP,0.0_DP,kind=DP)

          ! Copy the covariant gradient into complex workspace array
          ! agrecocmplx
          ! complex case
          if (loc_cmplx) then
             coarse_work(1:n1,1:n2,1:n3) = &
                 fftbox_batch(is,nmat,batch_count)%z(1:n1,1:n2,1:n3)
          ! real case
          else
             coarse_work(1:n1,1:n2,1:n3) = &
                 cmplx(fftbox_batch(is,nmat,batch_count)%d(1:n1,1:n2,1:n3), &
                 0.0_DP,kind=DP)
          end if

          ! Forward FFT the covariant gradient to reciprocal space
          call fourier_apply_box('Coarse','Forward',coarse_work, &
                      omp=pub_threads_fftbox)

          coarse_work = coarse_work * precond_func_recip

          ! Reset number of overlaps and projectors for this fa
          novlp = 0

          ovlp_species = 0
          ovlp_atoms = 0

          ! Loop over atoms with projectors overlapping fa
          do idx=ps_overlap_idx(loc_fa_atom),ps_overlap_idx(loc_fa_atom+1)-1
             ovlp_atom = ps_overlap_idx(idx)
             novlp = novlp + 1

             ! If fa is on a new atom, we need to update the entries in the
             ! phases arrays for this overlap
             if (loc_fa_atom/=loc_fa_atom_prev) &
                  call internal_phases(phases1,phases2,phases3, &
                  ovlp_atom,novlp,local_fa,fa_box_start(1,batch_count), &
                  fa_box_start(2,batch_count), &
                  fa_box_start(3,batch_count), proj_set,ngwf_basis, &
                  cell,fftbox)

             ! Store information about the species of this atom
             ovlp_species(novlp) = proj_set%proj_species(ovlp_atom)
             ! Store global index of overlapping atom
             ovlp_atoms(novlp) = ovlp_atom

          end do

          ! agrecocmplx
          ! complex case
          if (loc_cmplx) then
             allocate(A_mat_cmplx(max_paw_proj_tot,novlp),stat=ierr)
             call utils_alloc_check('paw_precond_contrib','A_mat_cmplx',ierr)
          ! real case
          else
             allocate(A_mat_real(max_paw_proj_tot,novlp),stat=ierr)
             call utils_alloc_check('paw_precond_contrib','A_mat_real',ierr)
          end if

          allocate(B_mat(max_paw_proj_tot,novlp),stat=ierr)
          call utils_alloc_check('paw_precond_contrib','B_mat',ierr)

          ! Zero reciprocal space gradient-contribution box
          out_complex = cmplx(0.0_DP,0.0_DP,kind=DP)

          ! Loop over atom species
          do isp=1,proj_set%n_proj_species


             ! Cycle if fa does not overlap any atoms of this species
             if (.not.any(ovlp_species(1:novlp)==isp)) cycle


             B_mat = cmplx(0.0_DP,0.0_DP,kind=DP)

             ! Loop over all the projector types for this species
             do iproj_1=1,proj_set%species_num_proj(isp)


                ! Find entry in global projector fftboxes array
                proj_count_1 = proj_set%species_first_proj(isp) + &
                     iproj_1 - 1

                if (.not.pub_projectors_precalculate) then
                   ! Find angular momentum of this projector
                   jproj = 1
                   found = .false.
                   do shell=1,proj_set%num_shells(isp)
                      am = proj_set%ang_mom(shell,isp)
                      do azim=-am,am
                         if (jproj==iproj_1) found=.true.
                         if (found) exit
                         jproj = jproj + 1
                      end do
                      if (found) exit
                   end do

                   proj_complex_1 = cmplx(0.0_DP,0.0_DP,kind=DP)

                   ! Generate projector in FFTbox
                   ! agrecokpt: call with specified k-point
                   call projector_in_box_recip(proj_complex_1, &
                        proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                        shell,isp), am, azim, fftbox%total_pt1, &
                        fftbox%total_pt2, fftbox%total_pt3, &
                        proj_set%n_rad_pts(isp), proj_set%gmax(isp), cell, &
                        kpt=loc_kpt)
                end if


                ! Loop over atoms with projectors overlapping fa
                do iovlp=1,novlp

                   ! Cycle if overlapping atom is not of the right species
                   if (.not.ovlp_species(iovlp)==isp) cycle

                   ! Loop over slabs of reciprocal space points of fftbox

                   if (pub_projectors_precalculate) then

                      do i3=1,n3
                         ! Loop over recip points in this slab of fftbox to
                         ! calculate structure factor for this projector type
                         do i2=1,n2
                            phase23 = phases2(i2,iovlp) * phases3(i3,iovlp)
                            do i1=1,n1
                               B_mat(iproj_1,iovlp) = B_mat(iproj_1,iovlp) + &
                                    conjg( phases1(i1,iovlp) * phase23 * &
                                    proj_set%fftbox_proj_recip(i1,i2,i3, &
                                    proj_count_1) ) * coarse_work(i1,i2,i3)
                            end do
                         end do
                      end do

                   else

                      do i3=1,n3
                         ! Loop over recip points in this slab of fftbox to
                         ! calculate structure factor for this projector type
                         do i2=1,n2
                            phase23 = phases2(i2,iovlp) * phases3(i3,iovlp)
                            do i1=1,n1
                               B_mat(iproj_1,iovlp) = B_mat(iproj_1,iovlp) + &
                                    conjg( phases1(i1,iovlp) * phase23 * &
                                    proj_complex_1(i1,i2,i3) ) * &
                                    coarse_work(i1,i2,i3)
                            end do
                         end do
                      end do


                   end if



                end do  ! End of loop over atoms with |p^j> overlapping fa


             end do !for iproj_1

             ! Reset A_mat to zero for new projector type
             ! agrecocmplx
             if (loc_cmplx) then
                A_mat_cmplx = (0.0_DP,0.0_DP)
             else
                A_mat_real = 0.0_DP
             end if

             ! Loop over all the projector types for this species
             do iproj=1,proj_set%species_num_proj(isp)

                ! Find entry in global projector fftboxes array
                proj_count = proj_set%species_first_proj(isp) + iproj - 1

                ! we must generate the projector in the fftbox even if
                ! projector_precalculate = T, since we need to apply
                ! the kinetic energy preconditioning to it

                ! Find angular momentum of this projector
                if (.not.pub_projectors_precalculate) then
                   jproj = 1
                   found = .false.
                   do shell=1,proj_set%num_shells(isp)
                      am = proj_set%ang_mom(shell,isp)
                      do azim=-am,am
                         if (jproj==iproj) found=.true.
                         if (found) exit
                         jproj = jproj + 1
                      end do
                      if (found) exit
                   end do

                   proj_complex = cmplx(0.0_DP,0.0_DP,kind=DP)

                   ! Generate projector in FFTbox
                   ! agrecokpt: call with specified k-point
                   call projector_in_box_recip(proj_complex, &
                        proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                        shell,isp), am, azim, n1, n2, n3, &
                        proj_set%n_rad_pts(isp), proj_set%gmax(isp), cell, &
                        kpt=loc_kpt)
                end if


                ! Reset structure factor to zero for new projector type
                struct_fac = cmplx(0.0_DP,0.0_DP,kind=DP)


                ! Loop over all the projector types for this species
                do iproj_1=1,proj_set%species_num_proj(isp)


                   ! Loop over atoms with projectors overlapping fa
                   do iovlp=1,novlp

                      ! Cycle if overlapping atom is not of the right species
                      if (.not.ovlp_species(iovlp)==isp) cycle

                      ! agrecocmplx
                      ! complex case
                      if (loc_cmplx) then
                         A_mat_cmplx(iproj,iovlp) = A_mat_cmplx(iproj,iovlp) + &
                               B_mat(iproj_1,iovlp) * &
                               proj_set%O_prec(iproj,iproj_1,isp) * &
                               normalization / fftbox%weight
                      ! real case
                      else
                         A_mat_real(iproj,iovlp) = A_mat_real(iproj,iovlp) + &
                               real(B_mat(iproj_1,iovlp)) * &
                               proj_set%O_prec(iproj,iproj_1,isp) * &
                               normalization / fftbox%weight
                      end if

                   end do  ! End of loop over atoms with |p^j> overlapping fa


                end do !for iproj_1


                ! Loop2 over atoms with projectors overlapping fa
                do iovlp=1,novlp

                   ! Cycle if overlapping atom is not of the right species
                   if (.not.ovlp_species(iovlp)==isp) cycle

                   do i3=1,n3

                         ! Loop over recip points in this slab of fftbox to
                         ! calculate structure factor for this projector type
                         ! agrecocmplx
                         if (loc_cmplx) then
                            do i2=1,n2
                               phase23 = phases2(i2,iovlp) * phases3(i3,iovlp)
                               struct_fac(1:n1,i2,i3) = struct_fac(1:n1,i2,i3) + &
                                    phases1(1:n1,iovlp) * phase23 * &
                                    A_mat_cmplx(iproj,iovlp)
                            end do
                         else
                            do i2=1,n2
                               phase23 = phases2(i2,iovlp) * phases3(i3,iovlp)
                               struct_fac(1:n1,i2,i3) = struct_fac(1:n1,i2,i3) + &
                                    phases1(1:n1,iovlp) * phase23 * &
                                    A_mat_real(iproj,iovlp)
                            end do
                         end if

                         ! Multiply the structure factor for this projector
                         ! type by the projector in reciprocal space, and add
                         ! it to the current reciprocal space gradient

                   end do



                end do  ! End of loop2 over atoms with |p^i> overlapping fa


                if (pub_projectors_precalculate) then

                   do i3=1,n3
                      ! Loop over recip points in this slab of fftbox to
                      ! calculate structure factor for this projector type
                      do i2=1,n2
                         do i1=1,n1

                            out_complex(i1,i2,i3) =  &
                             out_complex(i1,i2,i3) + &
                              struct_fac(i1,i2,i3) * &
                               proj_set%fftbox_proj_recip( &
                               i1,i2,i3,proj_count)

                         end do
                      end do
                   end do

                else

                   do i3=1,n3
                      ! Loop over recip points in this slab of fftbox to
                      ! calculate structure factor for this projector type
                      do i2=1,n2
                         do i1=1,n1
                            out_complex(i1,i2,i3) =  &
                             out_complex(i1,i2,i3) + &
                              struct_fac(i1,i2,i3) * &
                              proj_complex(i1,i2,i3)
                         end do
                      end do
                   end do

                end if



             end do  ! iproj

          end do  ! species_number


          ! agrecocmplx
          if (loc_cmplx) then
             deallocate(A_mat_cmplx,stat=ierr)
             call utils_dealloc_check('paw_precond_contrib','A_mat_cmplx', ierr)
          else
             deallocate(A_mat_real,stat=ierr)
             call utils_dealloc_check('paw_precond_contrib','A_mat_real', ierr)
          end if

          deallocate(B_mat,stat=ierr)
          call utils_dealloc_check('paw_precond_contrib','B_mat', ierr)

          ! Store previous local index of fa atom, so that we can tell if
          ! phases arrays need regenerating
          loc_fa_atom_prev = loc_fa_atom


          ! Backward FFT the covariant gradient to real space
          call fourier_apply_box('Coarse','Backward',out_complex, &
                      omp=pub_threads_fftbox)

          ! Copy the preconditioned covariant gradient out of the complex
          ! workspace array
          ! agrecocmplx
          if (loc_cmplx) then
             fftbox_batch(is,nmat,batch_count)%z(:,:,:) = &
                  fftbox_batch(is,nmat,batch_count)%z(:,:,:) + &
                  out_complex(:,:,:)
          else
             fftbox_batch(is,nmat,batch_count)%d(:,:,:) = &
                  fftbox_batch(is,nmat,batch_count)%d(:,:,:) + &
                  real(out_complex(:,:,:), kind=DP)
          end if
       end do ! loop over NGWFs in batch
!$OMP END DO
    end do ! loop over spins

    ! Deallocate Memory
    deallocate(coarse_work,stat=ierr)
    call utils_dealloc_check('projectors_grad_precond_batch','coarse_work', &
         ierr)
    deallocate(phases3,stat=ierr)
    call utils_dealloc_check('projectors_grad_precond_batch','phases3', &
         ierr)
    deallocate(phases2,stat=ierr)
    call utils_dealloc_check('projectors_grad_precond_batch','phases2', &
         ierr)
    deallocate(phases1,stat=ierr)
    call utils_dealloc_check('projectors_grad_precond_batch','phases1', &
         ierr)
    deallocate(ovlp_species,stat=ierr)
    call utils_dealloc_check('projectors_grad_precond_batch','ovlp_species', &
         ierr)
    deallocate(ovlp_atoms,stat=ierr)
    call utils_dealloc_check('projectors_grad_precond_batch','ovlp_atoms', &
         ierr)
    deallocate(struct_fac,stat=ierr)
    call utils_dealloc_check('projectors_grad_precond_batch','struct_fac',&
         ierr)
    deallocate(proj_complex,stat=ierr)
    call utils_dealloc_check('projectors_grad_precond_batch','proj_complex',&
         ierr)
    deallocate(proj_complex_1,stat=ierr)
    call utils_dealloc_check('projectors_grad_precond_batch','proj_complex_1',&
         ierr)
    deallocate(out_complex,stat=ierr)
    call utils_dealloc_check('projectors_grad_precond_batch','out_complex',&
         ierr)

!$OMP END PARALLEL

    if (pub_projectors_precalculate) then
       call projectors_exit_fftbox_recip(proj_set)
    end if


    deallocate(ps_overlap_idx,stat=ierr)
    call utils_dealloc_check('projectors_grad_precond_batch','ps_overlap_idx',&
         ierr)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving projectors_grad_precond_batch'

    ! cks: time it
    call timer_clock('projectors_grad_precond_batch', 2)

  end subroutine projectors_grad_precond_batch

  ! agrecokpt: should this routine be kpt dependent?
  ! is it sufficient to call projectors_init_fftbox_recip with kpt?
  subroutine projectors_proj_overlap_real(proj_overlap,proj_set, &
         precond_func_recip, cell,fftbox)

    !==================================================================!
    ! This subroutine returns the PAW-projector overlap matrix in      !
    ! normal matrix format.                                            !
    !------------------------------------------------------------------!
    !------------------------------------------------------------------!
    ! This subroutine was written by G. C. Constantinescu on 10/05/14. !
    !==================================================================!

    use constants, only: DP
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_projectors_precalculate
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(PROJECTOR_SET), intent(inout) :: proj_set
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(CELL_INFO), intent(in) :: cell
    real(kind=DP), intent(inout) :: proj_overlap( &
        maxval(proj_set%species_num_proj(:)), maxval( &
        proj_set%species_num_proj(:)), par%num_pspecies)
    real(kind=DP), intent(in)  :: precond_func_recip(:,:,:)

    ! Local Variables
    integer :: i1,i2,i3
    integer :: isp
    integer :: iproj,iproj_1
    integer :: proj_count_1
    integer :: ierr
    integer :: n1,n2,n3,ld1,ld2
    logical :: found
    integer :: jproj, shell, am, azim
    complex(kind=DP), allocatable :: proj_complex(:,:,:)
    complex(kind=DP), allocatable :: proj_complex_1(:,:,:)
    real(kind=DP) normalization
    complex(kind=DP) :: sum

    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2

    allocate(proj_complex(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('projectors_proj_overlap_real','proj_complex',ierr)

    allocate(proj_complex_1(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('projectors_proj_overlap_real','proj_complex_1',ierr)

    ! calculate the projector overlaps only over the different species,
    ! since <p^i|(1+T)^-1|p^j> is the same for two different atoms of
    ! the same species and the same i,j

    normalization = 1.0_DP / real(n1 * n2 * n3,kind=DP)
    proj_overlap = 0.0_DP

    if (pub_projectors_precalculate) then
        call projectors_init_fftbox_recip(proj_set,cell,fftbox)
    end if

    do isp=1,par%num_pspecies

       ! Double loop over projectors i,j
       do iproj=1,proj_set%species_num_proj(isp)

          ! Find angular momentum of this projector
          jproj = 1
          found = .false.
          do shell=1,proj_set%num_shells(isp)
             am = proj_set%ang_mom(shell,isp)
             do azim=-am,am
                if (jproj==iproj) found=.true.
                if (found) exit
                jproj = jproj + 1
             end do
             if (found) exit
          end do

          proj_complex = cmplx(0.0_DP,0.0_DP,kind=DP)

          ! Generate projector in FFTbox
          call projector_in_box_recip(proj_complex, &
               proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
               shell,isp), am, azim, n1, n2, n3, &
               proj_set%n_rad_pts(isp), proj_set%gmax(isp), cell)

          ! apply preconditining to projector |p^i>
          ! i.e. (1+T)^-1|p^i>
          ! the reciprocal preconditioning function has
          ! already been initialised in ngwf_gradient_batch

          proj_complex = proj_complex * precond_func_recip

          ! we will use the complex conjugate of |p^i> later

          do iproj_1=1,iproj

             proj_count_1 = proj_set%species_first_proj(isp) + iproj_1 - 1

             if (.not.pub_projectors_precalculate) then

                ! Find angular momentum of this projector
                jproj = 1
                found = .false.
                do shell=1,proj_set%num_shells(isp)
                   am = proj_set%ang_mom(shell,isp)
                   do azim=-am,am
                      if (jproj==iproj_1) found=.true.
                      if (found) exit
                      jproj = jproj + 1
                   end do
                   if (found) exit
                end do

                proj_complex_1 = cmplx(0.0_DP,0.0_DP,kind=DP)


                ! Generate projector in FFTbox
                call projector_in_box_recip(proj_complex_1, &
                     proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                     shell,isp), am, azim, n1, n2, n3, &
                     proj_set%n_rad_pts(isp), proj_set%gmax(isp), &
                     cell)

             end if


             ! Find C_ij = <p^i|(1+T)^-1|p^j> by doing the integration as
             ! a sum in reciprocal space, on the projector FFTboxes

             ! No phase factor needed, since it is the same origin

             ! Initialize sum
             sum = cmplx(0.0_DP,0.0_DP,kind=DP)

             if (pub_projectors_precalculate) then

                do i3=1,n3
                  do i2=1,n2
                    do i1=1,n1
                       sum = sum + conjg(proj_complex(i1,i2,i3)) *  &
                               proj_set%fftbox_proj_recip(i1,i2,i3, &
                               proj_count_1) * normalization / fftbox%weight
                    end do
                  end do
                end do

             else

                do i3=1,n3
                  do i2=1,n2
                    do i1=1,n1
                       sum = sum + conjg(proj_complex(i1,i2,i3)) *  &
                             proj_complex_1(i1,i2,i3) * normalization &
                             / fftbox%weight
                    end do
                  end do
                end do

             end if

             proj_overlap(iproj,iproj_1,isp) =  real(sum,kind=DP)
             ! symmetrize (since the result is real)

             proj_overlap(iproj_1,iproj,isp) = proj_overlap(iproj,iproj_1,isp)

          end do
       end do

     end do

     if (pub_projectors_precalculate) then
        call projectors_exit_fftbox_recip(proj_set)
     end if

     deallocate(proj_complex,stat=ierr)
     call utils_dealloc_check('projectors_proj_overlap_real','proj_complex',ierr)

     deallocate(proj_complex_1,stat=ierr)
     call utils_dealloc_check('projectors_proj_overlap_real','proj_complex_1',ierr)

  end subroutine projectors_proj_overlap_real

  ! agrecokpt: should this routine be kpt dependent?
  ! is it sufficient to call projectors_init_fftbox_recip with kpt?
  subroutine projectors_proj_overlap_cmplx(proj_overlap,proj_set,precond_func_recip, &
      cell,fftbox,kpt)

    !==================================================================!
    ! This subroutine returns the PAW-projector overlap matrix in      !
    ! normal matrix format (complex case).                             !
    !------------------------------------------------------------------!
    !------------------------------------------------------------------!
    ! This subroutine was written by Andrea Greco on 09/05/16.         !
    !==================================================================!

    use constants, only: DP
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_projectors_precalculate
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(PROJECTOR_SET), intent(inout) :: proj_set
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(CELL_INFO), intent(in) :: cell
    complex(kind=DP), intent(inout) :: proj_overlap( &
        maxval(proj_set%species_num_proj(:)), maxval( &
        proj_set%species_num_proj(:)), par%num_pspecies)
    real(kind=DP), intent(in)  :: precond_func_recip(:,:,:)
    ! agrecokpt
    type(POINT), optional, intent(in) :: kpt

    ! Local Variables
    integer :: i1,i2,i3
    integer :: isp
    integer :: iproj,iproj_1
    integer :: proj_count_1
    integer :: ierr
    integer :: n1,n2,n3,ld1,ld2
    logical :: found
    integer :: jproj, shell, am, azim
    complex(kind=DP), allocatable :: proj_complex(:,:,:)
    complex(kind=DP), allocatable :: proj_complex_1(:,:,:)
    real(kind=DP) normalization
    complex(kind=DP) :: sum
    ! agrecokpt
    type(POINT) :: loc_kpt

    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2

    allocate(proj_complex(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('projectors_proj_overlap_cmplx','proj_complex',ierr)

    allocate(proj_complex_1(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('projectors_proj_overlap_cmplx','proj_complex_1',ierr)

    ! calculate the projector overlaps only over the different species,
    ! since <p^i|(1+T)^-1|p^j> is the same for two different atoms of
    ! the same species and the same i,j

    normalization = 1.0_DP / real(n1 * n2 * n3,kind=DP)
    proj_overlap = (0.0_DP,0.0_DP)

    if (pub_projectors_precalculate) then
        ! agrecokpt
        call projectors_init_fftbox_recip(proj_set,cell,fftbox,kpt=loc_kpt)
    end if

    do isp=1,par%num_pspecies

       ! Double loop over projectors i,j
       do iproj=1,proj_set%species_num_proj(isp)

          ! Find angular momentum of this projector
          jproj = 1
          found = .false.
          do shell=1,proj_set%num_shells(isp)
             am = proj_set%ang_mom(shell,isp)
             do azim=-am,am
                if (jproj==iproj) found=.true.
                if (found) exit
                jproj = jproj + 1
             end do
             if (found) exit
          end do

          proj_complex = cmplx(0.0_DP,0.0_DP,kind=DP)

          ! Generate projector in FFTbox
          ! agrecokpt
          call projector_in_box_recip(proj_complex, &
               proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
               shell,isp), am, azim, n1, n2, n3, &
               proj_set%n_rad_pts(isp), proj_set%gmax(isp), cell, &
               kpt=loc_kpt)

          ! apply preconditining to projector |p^i>
          ! i.e. (1+T)^-1|p^i>
          ! the reciprocal preconditioning function has
          ! already been initialised in ngwf_gradient_batch

          proj_complex = proj_complex * precond_func_recip

          ! we will use the complex conjugate of |p^i> later

          do iproj_1=1,iproj

             proj_count_1 = proj_set%species_first_proj(isp) + iproj_1 - 1

             if (.not.pub_projectors_precalculate) then

                ! Find angular momentum of this projector
                jproj = 1
                found = .false.
                do shell=1,proj_set%num_shells(isp)
                   am = proj_set%ang_mom(shell,isp)
                   do azim=-am,am
                      if (jproj==iproj_1) found=.true.
                      if (found) exit
                      jproj = jproj + 1
                   end do
                   if (found) exit
                end do

                proj_complex_1 = cmplx(0.0_DP,0.0_DP,kind=DP)


                ! Generate projector in FFTbox
                ! agrecokpt
                call projector_in_box_recip(proj_complex_1, &
                     proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                     shell,isp), am, azim, n1, n2, n3, &
                     proj_set%n_rad_pts(isp), proj_set%gmax(isp), &
                     cell, kpt=loc_kpt)

             end if


             ! Find C_ij = <p^i|(1+T)^-1|p^j> by doing the integration as
             ! a sum in reciprocal space, on the projector FFTboxes

             ! No phase factor needed, since it is the same origin

             ! Initialize sum
             sum = cmplx(0.0_DP,0.0_DP,kind=DP)

             if (pub_projectors_precalculate) then

                do i3=1,n3
                  do i2=1,n2
                    do i1=1,n1
                       sum = sum + conjg(proj_complex(i1,i2,i3)) *  &
                               proj_set%fftbox_proj_recip(i1,i2,i3, &
                               proj_count_1) * normalization / fftbox%weight
                    end do
                  end do
                end do

             else

                do i3=1,n3
                  do i2=1,n2
                    do i1=1,n1
                       sum = sum + conjg(proj_complex(i1,i2,i3)) *  &
                             proj_complex_1(i1,i2,i3) * normalization &
                             / fftbox%weight
                    end do
                  end do
                end do

             end if

             ! agrecocmplx: mantain full complex character
             proj_overlap(iproj,iproj_1,isp) =  sum
             ! take the conjugate and symmetrize

             proj_overlap(iproj_1,iproj,isp) = conjg(proj_overlap(iproj,iproj_1,isp))

          end do
       end do

     end do

     if (pub_projectors_precalculate) then
        call projectors_exit_fftbox_recip(proj_set)
     end if

     deallocate(proj_complex,stat=ierr)
     call utils_dealloc_check('projectors_proj_overlap_cmplx','proj_complex',ierr)

     deallocate(proj_complex_1,stat=ierr)
     call utils_dealloc_check('projectors_proj_overlap_cmplx','proj_complex_1',ierr)

  end subroutine projectors_proj_overlap_cmplx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_commutator(nonlocpot_com, proj_basis, &
       ngwf_basis, ngwfs_on_grid, sp_overlap, proj_set, cell, fftbox, &
       delta, dij, kpt) ! agrecokpt

    !==================================================================!
    ! This subroutine calculates the commutator between the nonlocal   !
    ! potential and the position operator for the 3 Cartesian          !
    ! directions.                                                      !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff February 2011.                         !
    ! Modified by Nicholas Hine in April 2011 to not use NGWF_REP and  !
    ! to simplify and cleanup existing code, and for changes to        !
    ! projector initialisation routines.                               !
    ! Moved to projectors_mod by Nicholas Hine in December 2011.       !
    ! Modified by Laura Ratcliff Dec 2011 to fix bug involving complex !
    ! matrices and phase factors.                                      !
    ! Modified by Andrea Greco May 2016 to add specification of kpt    !
    ! and include compatibility with complex NGWFs.                    !
    ! Simplified by Jose M Escartin in September 2017.                 !
    ! Modified for embedding by Joseph Prentice, July 2018             !
    !==================================================================!

    use datatypes, only: FUNCTIONS
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: CELL_INFO
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product, sparse_embed_transpose, &
         sparse_embed_axpy, sparse_embed_scale, sparse_embed_copy, &
         sparse_embed_transpose_structure
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: nonlocpot_com(3)
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid(:)
    type(SPAM3_EMBED), intent(in) :: sp_overlap
    type(PROJECTOR_SET), intent(inout) :: proj_set(:)
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(in) :: delta
    type(SPAM3_EMBED), intent(in) :: dij
    ! agrecokpt
    type(POINT), optional, intent(in) :: kpt

    ! Local Variables
    type(SPAM3_EMBED) :: sp_overlap_dij, ps_overlap, dij_ps_overlap
    type(SPAM3_EMBED) :: sDGp_overlap, DGps_overlap, sp_overlap_cmplx
    type(SPAM3_EMBED) :: nonlocpot_com2
    integer :: cart, isub, jsub
    logical :: ovlp_iscmplx
    logical :: mix_rc
    ! agrecokpt
    type(POINT) :: loc_kpt

    ! agrecocmplx
    ovlp_iscmplx = sp_overlap%iscmplx

    ! agrecokpt
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    ! agrecokpt: specification of kpt allowed only in complex case
    if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.&
        loc_kpt%z/=0.0_DP) then
       if (.not.ovlp_iscmplx) then
          call utils_abort('Error in projectors_commutator: &
               &non-Gamma k-point requires complex projectors.')
       end if
    end if

    ! NB not allowing for sp_overlap or ps_overlap to be at diff k-points
    call sparse_embed_create(nonlocpot_com2,nonlocpot_com(1),iscmplx=.true.)
    call sparse_embed_create(sp_overlap_dij,sp_overlap, &
         iscmplx=(dij%iscmplx.or.ovlp_iscmplx))
    call sparse_embed_create(sDGp_overlap,sp_overlap,iscmplx=.true.)

    ! ndmh: Create ps_overlap (modifying structure code if required)
    call sparse_embed_transpose_structure(ps_overlap%structure,sp_overlap)
    dij_ps_overlap%structure = ps_overlap%structure
    DGps_overlap%structure = ps_overlap%structure

    call sparse_embed_create(dij_ps_overlap,iscmplx=sp_overlap_dij%iscmplx)
    call sparse_embed_create(DGps_overlap,iscmplx=.true.)
    call sparse_embed_create(ps_overlap,iscmplx=sp_overlap%iscmplx)

    ! Create the transpose of the NGWF-projector overlap matrix
    call sparse_embed_transpose(ps_overlap,sp_overlap)

    mix_rc = (ovlp_iscmplx.neqv.dij%iscmplx)
    call sparse_embed_product(sp_overlap_dij,sp_overlap,dij, &
         allow_mix_types=mix_rc)
    call sparse_embed_product(dij_ps_overlap,dij,ps_overlap, &
         allow_mix_types=mix_rc)

    call sparse_embed_destroy(ps_overlap)

    ! Loop over Cartesian directions
    do cart=1,3

       ! Calculate <phi|Del_G(proj(G))> overlap matrix
       ! compute directly for complex case
       ! agrecokpt: at specified k-point
       ! jcap: loop over regions
       do isub=1,size(ngwf_basis(:))
          do jsub=1,size(proj_basis(:))
             call projectors_func_ovlp_box(sDGp_overlap%m(isub,jsub), &
                  ngwfs_on_grid(isub), ngwf_basis(isub), &
                  proj_basis(jsub), proj_set(jsub), fftbox, cell, &
                  kshift=loc_kpt, delta=delta, cart=cart)
          end do
       end do

       ! Transpose <phi|Del_G(proj(G))> to get <Del_G(proj(G))|phi> for ngwf_basis
       call sparse_embed_transpose(DGps_overlap,sDGp_overlap)

       mix_rc = (.not.sp_overlap_dij%iscmplx)
       ! Calculate the matrix \sum_i (<NGWF_a|Del_G(Proj(G))_i>D_ij<Proj_j|NGWF_b>)
       call sparse_embed_product(nonlocpot_com2,sDGp_overlap,dij_ps_overlap, &
            allow_mix_types=mix_rc)

       ! Calculate the matrix \sum_i (<NGWF_a|Proj_i>D_ij<Del_G(Proj(G))_i|NGWF_b>)
       call sparse_embed_product(nonlocpot_com(cart),sp_overlap_dij, &
            DGps_overlap, allow_mix_types=mix_rc)

       ! Sum the two
       call sparse_embed_axpy(nonlocpot_com(cart),nonlocpot_com2,(1.0_DP,0.0_DP))

       ! Multiply by -i
       call sparse_embed_scale(nonlocpot_com(cart),(0.0_DP,-1.0_DP))

    end do

    call sparse_embed_destroy(DGps_overlap)
    call sparse_embed_destroy(dij_ps_overlap)
    call sparse_embed_destroy(sDGp_overlap)
    call sparse_embed_destroy(sp_overlap_dij)
    call sparse_embed_destroy(nonlocpot_com2)

  end subroutine projectors_commutator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module projectors
