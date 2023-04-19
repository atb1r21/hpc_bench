! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                       Visualisation module                     !
!                                                                !
! This module contains subroutines for producing and outputing   !
! quantities such as the charge density or the NGWFs in various  !
! commonly-used formats for plotting.                            !
!----------------------------------------------------------------!
! Written by Chris-Kriton Skylaris, 01/07/2005                   !
!================================================================!

module visual

  implicit none

  private

  public :: visual_scalarfield
  public :: visual_scalarfield_ngwf_basis
  public :: visual_ngwfs
  public :: visual_ngwfs_radial
  public :: visual_efield_on_grid_pbc
  public :: visual_efield_on_grid_obc
  public :: visual_scalarfield_read

  ! jd: This is set by auto-solvation to temporarily change output filenames
  !     so that in-vacuum scalarfields are not overwritten by in-solvent
  !     scalarfield later on. Not pretty, but the alternative is to pass flags
  !     through a forest of calls.
  character(len=7), public, save :: scalarfield_write_filename_prefix = ''

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine visual_ngwfs_radial( &
       ngwfs_on_grid, ngwf_basis, ngwfs_name, ngwfs_type, fftbox, cell, par) ! input

    !===================================================================!
    ! This subroutine outputs radial distributions of NGWFs to files in !
    ! column format. Output happens in parallel (from all processors)   !
    ! and the NGWFs are numbered in their unique order after the space  !
    ! filling rearrangement of atoms.                                   !
    !-------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 06/04/2011                       !
    ! pub_par global variable removed by Robert Charlton, 12/06/2018.   !
    !===================================================================!

    use basis, only: basis_ppd_location
    use comms, only: pub_on_root, pub_my_proc_id
    use constants, only: DP, stdout, PI
    use datatypes, only: FUNCTIONS
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use geometry, only: point, operator(+), operator(*), local_displacement, &
         geometry_distance
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_rootname, pub_write_radial_step, &
         pub_write_radial_smear
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_abort, utils_assert

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid ! NGWFs
    character(len=*), intent(in) :: ngwfs_name ! beginning of NGWFs names
    character(len=*), intent(in) :: ngwfs_type ! beginning of NGWFs names
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) ::fftbox
    type(PARAL_INFO), intent(in) :: par

    ! Local Variables
    integer :: igrid, rgrid_npts
    integer :: point_counter
    integer :: radial_counter
    integer :: ppd
    integer :: loc_ppd
    integer :: loc_1, loc_2, loc_3
    integer :: a1_neighbour, a2_neighbour, a3_neighbour
    integer :: a1_upper,a2_upper,a3_upper
    integer :: a1_lower,a2_lower,a3_lower
    integer :: a1_cell,a2_cell,a3_cell
    integer :: ppd_count
    real(kind=DP) :: local_1, local_2, local_3
    real(kind=DP) :: radial_distance, rrad, fac
    real(kind=DP), allocatable :: func_on_rgrid(:)
    type(POINT) :: current_point, current_displacement
    type(POINT) :: periodic_centre, temp_centre, ppd_origin
    integer :: ngwf_counter   ! local proc NGWF counter
    integer :: orig_atom      ! counter of atoms in input file order
    integer :: loc_atom       ! local atom counter
    integer :: glob_atom      ! global atom (all processors) counter
    integer :: atom_func      ! counter of NGWFs on an atom
    character(len=256) :: fun_name
    integer :: output_unit, ierr
    real(kind=DP), parameter :: eps_tol = 1D-12
    ! agrecocmplx
    logical :: loc_cmplx

    ! Print info about output operation
    if (pub_on_root) write(stdout,'(a)',advance ='no') &
         'Writing NGWFs radial distributions.'

    ! agrecocmplx: now made compatible, check everything is ok
    !call utils_assert(.not. ngwfs_on_grid%iscmplx, 'Error in&
    !     & visual_ngwfs_radial: subroutine not ready yet for complex NGWFs.')

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid%iscmplx

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

    ! Loop over NGWFs stored on local proc
    ngwf_counter = 0

    LOOP_ON_ATOMS : do loc_atom = 1, par%num_atoms_on_proc(pub_my_proc_id)

       glob_atom = loc_atom + par%first_atom_on_proc(pub_my_proc_id) - 1
       orig_atom = par%orig_atom(glob_atom)

       LOOP_ON_NGWFS : do atom_func = 1, ngwf_basis%num_on_atom(glob_atom)

          ngwf_counter = ngwf_counter + 1

          ! Create radial grid
          rgrid_npts = int(ngwf_basis%spheres(ngwf_counter)%radius &
               * 1.1_dp / pub_write_radial_step)

          ! Allocate space to store the radial data
          allocate(func_on_rgrid(rgrid_npts), stat=ierr)
          call utils_alloc_check('visual_ngwfs_radial','func_on_rgrid',ierr)
          func_on_rgrid(:) = 0.0_DP

          ! Loop over ppds that belong to the NGWF sphere
          point_counter = ngwf_basis%spheres(ngwf_counter)%offset-1
          radial_counter = 0

          LOOP_ON_PPDS : do ppd_count = 1, ngwf_basis%spheres(ngwf_counter)%n_ppds_sphere

             ppd =  ngwf_basis%spheres(ngwf_counter)%ppd_list(1,ppd_count)
             ppd_origin = basis_ppd_location(ppd,cell)

             loc_ppd = ngwf_basis%spheres(ngwf_counter)%ppd_list(2,ppd_count)

             ! cks: loop over points of this ppd
             do loc_3 = 0, cell%n_pt3 -1
                local_3 = real(loc_3, DP)*cell%d3

                do loc_2 = 0, cell%n_pt2-1
                   local_2 = real(loc_2, DP)*cell%d2

                   do loc_1 = 0, cell%n_pt1-1
                      local_1 = real(loc_1, DP)*cell%d1

                      current_displacement = local_displacement(cell%a1_unit,&
                           cell%a2_unit,cell%a3_unit,local_1,local_2,local_3)

                      current_point = ppd_origin + current_displacement

                      point_counter = point_counter + 1
                      radial_counter = radial_counter + 1

                      ! ndmh: loop over possible origins of this NGWF from different cells
                      ! ndmh: NB: a1_lower = a1_upper = 0 if fftbox%coin1=.false.
                      periodic_centre = ngwf_basis%spheres(ngwf_counter)%centre

                      do a1_cell=a1_lower,a1_upper
                         do a2_cell=a2_lower,a2_upper
                            do a3_cell=a3_lower,a3_upper

                               a1_neighbour=nint(real(loc_ppd,kind=DP)/9.0_DP)
                               a2_neighbour=nint(real(loc_ppd-9*a1_neighbour,DP)/3.0_DP)
                               a3_neighbour=loc_ppd-9*a1_neighbour-3*a2_neighbour

                               ! ndmh: In cases where FFTbox==cell, check over nearest
                               ! ndmh: neighbours regardless of loc_ppd value
                               if (fftbox%coin1) a1_neighbour = a1_cell
                               if (fftbox%coin2) a2_neighbour = a2_cell
                               if (fftbox%coin3) a3_neighbour = a3_cell

                               temp_centre = ngwf_basis%spheres(ngwf_counter)%centre &
                                    + real(a1_neighbour,DP)*cell%a1 &
                                    + real(a2_neighbour,DP)*cell%a2 &
                                    + real(a3_neighbour,DP)*cell%a3

                               if (geometry_distance(temp_centre,current_point) .lt. &
                                   geometry_distance(periodic_centre,current_point)) then
                                  periodic_centre = temp_centre
                               endif

                            enddo
                         enddo
                      enddo

                      radial_distance = geometry_distance(periodic_centre,current_point)

                      if (abs(radial_distance) < eps_tol) cycle
                      do igrid = 1, rgrid_npts
                         rrad = (igrid * pub_write_radial_step) - radial_distance
                         ! agrecocmplx: complex case, plot the abs of the value
                         if (loc_cmplx) then
                            func_on_rgrid(igrid) = func_on_rgrid(igrid) + &
                                (abs(ngwfs_on_grid%z(point_counter))/(radial_distance**2))* &
                                exp(-rrad**2/pub_write_radial_smear)
                         ! real case
                         else
                            func_on_rgrid(igrid) = func_on_rgrid(igrid) + &
                                (abs(ngwfs_on_grid%d(point_counter))/(radial_distance**2))* & !jmecmplx
                                exp(-rrad**2/pub_write_radial_smear)
                         end if
                      enddo

                   enddo
                enddo
             enddo

          enddo LOOP_ON_PPDS

          ! Compose output file name
          write(fun_name,'(a,i5.5,a,i2.2,a)') &
               trim(pub_rootname)//'_'//trim(adjustl(ngwfs_name))//trim('_atom'), orig_atom, &
               trim('_')//trim(adjustl(ngwfs_type)), atom_func, trim('.radial')
          fun_name = adjustl(fun_name)

          ! Get a unit number that is free
          output_unit = utils_unit()
          open(unit=output_unit,file=fun_name,err=10)

          ! Output the radial data
          fac = sqrt(pi*pub_write_radial_smear)*4*pi
          do igrid = 1, rgrid_npts
             rrad = igrid * pub_write_radial_step
             write(output_unit,*) rrad, func_on_rgrid(igrid)/fac
          enddo

          ! Close output file
          close(output_unit,err=20)

          ! Deallocate radial grid
          deallocate(func_on_rgrid, stat=ierr)
          call utils_dealloc_check('visual_ngwfs_radial','func_on_rgrid',ierr)

       enddo LOOP_ON_NGWFS

    enddo LOOP_ON_ATOMS

    return

10  call utils_abort('Error during opening of file:'//trim(fun_name)//&
         ' in visual_ngwfs_radial().')

20  call utils_abort('Error during closing of file:'//trim(fun_name)//&
         ' in visual_ngwfs_radial().')

  end subroutine visual_ngwfs_radial


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine visual_scalarfield( &
       scalarfield_fine, grid, cell, file_header, scalar_name, &  ! input
       elements, conversion_factor) ! input, optional

    !=====================================================================!
    ! This subroutine outputs a real scalar field (e.g. density,          !
    ! potential, real molecular orbital, etc ) into files in formats      !
    ! suitable for visualisation.                                         !
    !---------------------------------------------------------------------!
    ! Warning: This subroutine allocates on pub_root_proc_id an array     !
    !          which has the size of the fine grid of the whole           !
    !          simulation cell. If your simulation cell is too large      !
    !          calling this routine will cause ONETEP to crash by         !
    !          running out of memory OR to output files on disk           !
    !          which are too large to store and/or plot with the          !
    !          available plotting software.                               !
    ! Note:    Above no longer applies in .cube format as output is split !
    !          over several steps. Easily extendable to .grd format.      !
    ! Note:    Above warning does not apply in .dx format as well.        !
    !---------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 22/06/2005.                     !
    ! Modified by Chris-Kriton Skylaris on 14/02/2006.                    !
    ! Cube output modified by Nicholas Hine on 22/02/2008.                !
    ! Added dx output, Jacek Dziedzic on 08/03/2010.                      !
    ! Added unit conversion functionality, Jacek Dziedzic on 14/05/2010.  !
    ! Removed dependency of elements on par%nat, Rob Bell 28/01/2014      !
    !=====================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_extract_box
    use comms, only: pub_on_root, comms_barrier
    use constants, only: DP, LONG, ANGSTROM, PI, stdout
    use geometry, only: geometry_magnitude, operator(.DOT.)
    use ion, only: element
    use rundat, only: pub_rootname, pub_cube_format, pub_dx_format, &
         pub_grd_format, pub_debug_on_root
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! cks: subroutine arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    real(kind=DP), intent(in) :: scalarfield_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12) ! scalar field to output
    character(len =*), intent(in) :: file_header ! part of header line of output
    ! file
    character(len =*), intent(in) :: scalar_name ! part of file-name to output
    !(scalarfield name)
    type(ELEMENT), intent(in), optional :: elements(:)
    real(kind=DP), intent(in), optional :: conversion_factor

    ! Local variables
    real(kind=DP), dimension(:, :, :), allocatable :: buffer_fine
    ! comms buffer
    real(kind=DP), dimension(:, :, :), allocatable :: scalarfield_box_fine
    ! simcell buffer
    real(kind=DP) :: aa, bb, cc   ! lattice vector lengths
    real(kind=DP) :: alpha, beta, gamma ! lattice angles
    character(len =256) :: output_file  ! file names
    character(len =256) :: title_line   ! title line in output file
    integer :: ldf1,ldf2,nf1,nf2,nf3    ! shorthand working variables
    integer :: ierr                     ! memory alloc error flag
    integer :: nstep                    ! number of steps to output array over
    integer :: istep       ! step counter for splitting output of arrays
    integer :: nslab1      ! number of slabs to write on this step (cube output)
    integer :: nslab1_full ! leading dim of scalarfield_box_fine (grd output)
    integer :: islab1      ! current slab index in '1' direction
    integer :: nslab3      ! number of slabs to write on this step
    integer :: islab3      ! current slab index in '3' direction
    logical :: writeheader ! whether we still need to write output file header
    logical :: grd_2in1    ! whether to output only half the points along each
                           ! axis when doing grd output

    ! ndmh: variables using long integers to avoid overflows
    ! maximum size of buffer on master proc (125000000 real*8's ~ 1GB)
    integer(kind=LONG),parameter :: max_size = 125000000_LONG
    integer(kind=LONG) :: total_size  ! total size of simulation cell data array
    integer(kind=LONG) :: nf1_long      ! long version of nf1
    integer(kind=LONG) :: nf3_long      ! long version of nf3

    ! jd: variables distinguishing CUBE and DX formats
    integer :: n_format ! Index in the loop over CUBE and DX formats
    character(len=4) :: suffix ! 'cube' or 'dx'

    character(len=*), parameter :: myself = 'visual_scalarfield'

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering '//myself//'.'

    ! cks: initialise shorthand variables
    ldf1  = grid%ld1
    ldf2  = grid%ld2
    nf1   = grid%n1
    nf2   = grid%n2
    nf3   = grid%n3

    ! cks:########### SCALAR FIELD OUTPUT ################################

    write(title_line,*) &
         trim(file_header) // ' "' // trim(pub_rootname) // '"'
    title_line = adjustl(title_line)

    ! cks: *** CUBE format ***
    ! cks: output charge density in .cube format
    ! jd:  modified to work with the DX format as well
    do n_format=1,2 ! jd: Loop over formats, to reuse code

       ! jd: Distinguish between cube and dx
       if(n_format == 1 .and. .not. pub_cube_format) cycle
       if(n_format == 2 .and. .not. pub_dx_format) cycle
       if(n_format == 1) suffix = 'cube'
       if(n_format == 2) suffix = 'dx'

       ! ndmh: find output file name
       if (pub_on_root) then
          write(output_file,'(a256)') trim(pub_rootname)//&
               trim(adjustl(scalarfield_write_filename_prefix))//&
               trim(adjustl(scalar_name))//'.'//suffix
          output_file = adjustl(output_file)

          write(stdout,'(3a)',advance ='no') &
               'Writing "', trim(output_file),'" ...'
       end if

       ! ndmh: calculate how many steps to divide output into
       total_size = nf1 * ldf2 * nf3
       if(total_size < max_size)then
          nstep = 1
          nslab1 = nf1
       else
          nf1_long = nf1
          nslab1 = int(nf1_long * max_size / total_size)
          if (nslab1==0) nslab1 = 1
          nslab1 = nslab1 + mod(nslab1,2) ! ensure it is even
          nstep = nf1 / nslab1
          if (nstep*nslab1 < nf1) nstep = nstep + 1
       end if

       if (pub_debug_on_root) then
          write(stdout,*)
          write(stdout,'(4a,i3,a,i3,a)') 'DEBUG: Writing ',suffix,' format ', &
                'simulation cell data in ', nstep, ' steps of ', nslab1, &
                ' slabs'
       endif

       ! ndmh: allocate simulation cell charge density buffer
       allocate(buffer_fine(nslab1,ldf2,grid%max_slabs12), &
            stat=ierr)
       call utils_alloc_check('visual_scalarfield','buffer_fine',ierr)
       allocate(scalarfield_box_fine(nslab1,ldf2,nf3),stat=ierr)
       call utils_alloc_check('visual_scalarfield','scalarfield_box_fine',ierr)

       ! ndmh: write header before first slab only
       writeheader=.true.
       nslab1_full = nslab1
       islab1 = 1

       do istep = 1, nstep

          ! ndmh: do not write beyond end of grid on last step
          if (islab1 + nslab1 - 1 > nf1) nslab1 = nf1 - islab1 + 1

          ! ndmh: extract a box of nslab1 '23'-slabs of cell data
          call cell_grid_extract_box(scalarfield_box_fine,&
               buffer_fine, scalarfield_fine, grid, &
               nslab1, nf2, nf3, nslab1_full, ldf2, islab1, 1, 1, &
               pub_on_root, .false.)

          ! jd: Apply unit conversion factor, if requested
          if (present(conversion_factor)) then
             scalarfield_box_fine = scalarfield_box_fine * conversion_factor
          end if

          ! ndmh: write out this batch of '23'-slabs
          if (pub_on_root .and. n_format == 1) call visual_output_cube( &
               scalarfield_box_fine, &
               "Output from ONETEP calculation", &
               trim(title_line), trim(output_file), 0.0_DP, 0.0_DP, 0.0_DP, &
               cell%a1, cell%a2, cell%a3, nf1, nf2, nf3, &
               nslab1, writeheader, elements)
          if (pub_on_root .and. n_format == 2) call visual_output_dx( &
               scalarfield_box_fine, &
               "Output from ONETEP calculation", &
               trim(title_line), trim(output_file), 0.0_DP, 0.0_DP, 0.0_DP, &
               cell%a1, cell%a2, cell%a3, nf1, nf2, nf3, &
               nslab1, writeheader, (istep == nstep) )

          writeheader=.false.
          islab1 = islab1 + nslab1

       end do

       ! ndmh: deallocate buffers
       deallocate(scalarfield_box_fine,stat=ierr)
       call utils_dealloc_check('visual_scalarfield','scalarfield_box_fine',ierr)
       deallocate(buffer_fine,stat=ierr)
       call utils_dealloc_check('visual_scalarfield','buffer_fine',ierr)

       if (pub_on_root) write(stdout,'(a)') ' done'

    end do ! jd: Loop over CUBE, DX formats

    ! cks: *** GRD format ***
    ! cks: output charge density in .grd format
    if (pub_grd_format) then

       ! ndmh: write standard grid with normal spacing, but for any other grid,
       ! ndmh: only output every other point, to save space
       grd_2in1 = ((grid%n1>cell%total_pt1) .or. &
            (grid%n2>cell%total_pt2) .or. &
            (grid%n3>cell%total_pt3))

       ! ndmh: calculate how many steps to divide output into
       total_size = ldf1 * ldf2 * nf3
       if(total_size < max_size)then
          nstep = 1
          nslab3 = nf3
       else
          nf3_long = nf3
          nslab3 = int(nf3_long * max_size / total_size)
          if (nslab3==0) nslab3 = 1
          nslab3 = nslab3 + mod(nslab3,2) ! ensure it is even
          nstep = nf3 / nslab3
          if (nstep*nslab3 < nf3) nstep = nstep + 1
       end if

       if (pub_debug_on_root) then
          write(stdout,*)
          write(stdout,'(a,a,i3,a,i3,a)') 'DEBUG: Writing grd format ', &
                'simulation cell data in ', nstep, ' steps of ', nslab3, &
                ' slabs'
       endif

       ! cks: Allocate workspace for slabs buffer
       allocate(buffer_fine(ldf1,ldf2,grid%max_slabs12),stat=ierr)
       call utils_alloc_check('visual_scalarfield','buffer_fine',ierr)
       ! cks: allocate simulation cell charge density buffer
       allocate(scalarfield_box_fine(ldf1,ldf2,nslab3),stat=ierr)
       call utils_alloc_check('visual_scalarfield','scalarfield_box_fine',ierr)

       if (pub_on_root) then

          ! Initialise
          scalarfield_box_fine = 0.0_DP

          ! cks: lattice vector lengths in angstrom
          aa =geometry_magnitude(cell%a1)/ANGSTROM
          bb =geometry_magnitude(cell%a2)/ANGSTROM
          cc =geometry_magnitude(cell%a3)/ANGSTROM

          ! cks: lattice angles in degrees
          alpha =acos( cell%a2_unit .DOT. cell%a3_unit )
          alpha =180.0_DP *alpha/PI
          beta =acos( cell%a1_unit .DOT. cell%a3_unit )
          beta =180.0_DP *beta/PI
          gamma =acos( cell%a1_unit .DOT. cell%a2_unit )
          gamma =180.0_DP *gamma/PI

          ! ndmh: find output file name
          write(output_file,'(a256)') trim(pub_rootname)//trim(&
               adjustl(scalar_name))//'.grd'
          output_file = adjustl(output_file)

          write(stdout,'(3a)',advance ='no') &
               'Writing "', trim(output_file),'" ...'

       end if

       ! ndmh: write header before first slab only
       writeheader=.true.
       islab3 = 1

       do istep = 1, nstep

          ! ndmh: do not write beyond end of grid on last step
          if (islab3 + nslab3 - 1 > nf3) nslab3 = nf3 - islab3 + 1

          ! ndmh: extract a box of nslab3 '12'-slabs of cell data
          call cell_grid_extract_box(scalarfield_box_fine,&
               buffer_fine, scalarfield_fine, grid, &
               nf1, nf2, nslab3, ldf1, ldf2, islab3, 1, 1, &
               pub_on_root, .false.)

          ! jd: Apply unit conversion factor, if requested
          if (present(conversion_factor)) then
             scalarfield_box_fine = scalarfield_box_fine * conversion_factor
          end if

          ! cks: Only one every second point is output to mantain manageable
          ! cks: file sizes
          ! ndmh: write out this batch of '12'-slabs
          if (pub_on_root) call visual_output_grd( &
               scalarfield_box_fine, aa, bb, cc, alpha, beta, gamma, &
               nf1, nf2, nf3, nslab3, 0, 0, 0, nf1, nf2, nf3, &
               trim(title_line), trim(output_file), grd_2in1, writeheader)

          writeheader = .false.
          islab3 = islab3 + nslab3

       end do

       ! ndmh: grd data is periodic and includes first slabs twice.
       ! ndmh: so we need to extract the first '12'-slab again
       call cell_grid_extract_box(scalarfield_box_fine,&
            buffer_fine, scalarfield_fine, grid, &
            nf1, nf2, 1, ldf1, ldf2, 1, 1, 1, pub_on_root, .false.)

       ! jd: Apply unit conversion factor here as well, if requested
       if (present(conversion_factor)) then
          scalarfield_box_fine = scalarfield_box_fine * conversion_factor
       end if

       ! ndmh: write first '12'-slab data to the end of the file
       if (pub_on_root) call visual_output_grd( &
            scalarfield_box_fine, aa, bb, cc, alpha, beta, gamma, &
            nf1, nf2, nf3, 1, 0, 0, 0, nf1, nf2, nf3, &
            trim(title_line), trim(output_file), grd_2in1, .false.)

       if (pub_on_root) write(stdout,*) ' done'

       ! cks: Deallocate memory
       deallocate(scalarfield_box_fine,stat=ierr)
       call utils_dealloc_check('visual_scalarfield','scalarfield_box_fine',ierr)

       deallocate(buffer_fine,stat=ierr)
       call utils_dealloc_check('visual_scalarfield','buffer_fine',ierr)

    endif

    ! cks:####### END SCALAR FIELD OUTPUT ################################

    call comms_barrier

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving '//myself//'.'

  end subroutine visual_scalarfield

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine visual_ngwfs( &
       ngwfs_on_grid, ngwf_basis, ngwfs_head_name, elements, cell, fftbox, par)   ! input

    !===================================================================!
    ! This subroutine outputs NGWFs to files in GRD and CUBE formats    !
    ! for plotting. Output happens in parallel (from all processors)    !
    ! and the NGWFs are numbered in their unique order after the space  !
    ! filling rearrangement of atoms.                                   !
    !-------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 02/06/2005.                   !
    ! Modified by Chris-Kriton Skylaris on 13/03/2006 to write only     !
    ! selected NGWFs                                                    !
    ! pub_par global variable removed by Robert Charlton, 12/06/2018.   !
    !===================================================================!

    use basis, only: basis_copy_function_to_box, &
        basis_location_func_wrt_cell
    use comms, only: comms_barrier, pub_on_root, pub_my_proc_id
    use constants, only: DP, stdout, PI, ANGSTROM
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
        data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use geometry, only: operator(.DOT.), point, operator(*), operator(+)
    use ion, only: element
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_cube_format, pub_dx_format, pub_grd_format, &
         pub_rootname, pub_nnho
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid ! NGWFs
    character(len=*), intent(in) :: ngwfs_head_name ! beginning of NGWFs names
    type(PARAL_INFO), intent(in) :: par             ! parallel strategy for these atoms
    type(ELEMENT), intent(in), optional :: elements(par%nat) ! all elements
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) ::fftbox

    ! cks: local variables
    type(FFTBOX_DATA) :: output_box ! tightbox buffer
    real(kind=DP) :: aa, bb, cc ! tightbox "lattice vector" lenghts
    real(kind=DP) :: alpha, beta, gamma ! tightbox "lattice vector" angles
    integer :: npt1, npt2, npt3 ! current tightbox number of points
    integer :: tight_start1 ! origin of tightbox wrt simcell in psinc points
    integer :: tight_start2 ! origin of tightbox wrt simcell in psinc points
    integer :: tight_start3 ! origin of tightbox wrt simcell in psinc points
    integer :: ngwf_counter ! local proc NGWF counter
    integer :: start1, start2, start3 ! start of tightbox from origin
    integer :: end1, end2, end3 ! end of tightbox from origin
    integer :: orig_atom   ! counter of atoms in input file order
    integer :: latom       ! local atom counter
    integer :: gatom       ! global atom (all processors) counter
    integer :: at_func     ! counter of NGWFs on an atom
    character(len =256) :: temp_name, fun_name ! string buffers
    character(len =256) :: temp_name1          ! string buffer
    character(len =6)   :: fun_number ! string buffer
    character(len =256) :: title_line ! line to appear in output file
    TYPE(POINT) :: a1c,a2c,a3c,origin_cube
    ! agrecocmplx
    logical :: loc_cmplx
    real(kind=DP), allocatable :: output_box_abs(:,:,:)
    integer :: ierr

    !jmecmplx
    ! agrecocmplx: now made compatible, need to check everything is ok
    !call utils_assert(.not. ngwfs_on_grid%iscmplx, 'Error in visual_ngwfs:&
    !     & not ready yet for complex NGWFs.')

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid%iscmplx

    ! pdh: universal size to avoid copies elsewhere in ONETEP
    call data_fftbox_alloc(output_box, ngwf_basis%maxtight_pts1, &
         ngwf_basis%maxtight_pts2, ngwf_basis%maxtight_pts3, &
         iscmplx=loc_cmplx)

    ! agrecocmplx: need real array to store abs of NGWFs on grid
    ! if we are in the complex case
    if (loc_cmplx) then
       allocate(output_box_abs(ngwf_basis%maxtight_pts1, ngwf_basis%maxtight_pts2, &
                ngwf_basis%maxtight_pts3),stat=ierr)
       call utils_alloc_check('visual_ngwfs','output_box_abs',ierr)
    end if

    ! cks: Initialise tightbox "cell angles"
    alpha =acos( cell%a2_unit .DOT. cell%a3_unit )
    alpha =180.0_DP *alpha/PI
    beta =acos( cell%a1_unit .DOT. cell%a3_unit )
    beta =180.0_DP *beta/PI
    gamma =acos( cell%a1_unit .DOT. cell%a2_unit )
    gamma =180.0_DP *gamma/PI

    ! cks: print info about output operation
    if (pub_nnho) then
       if (pub_on_root) write(stdout,'(a)',advance ='no') &
            'Writing NNHO plot files in formats:'
    else
       if (pub_on_root) write(stdout,'(a)',advance ='no') &
            'Writing NGWF plot files in formats:'
    endif
    if (pub_grd_format .and. pub_on_root) &
         write(stdout,'(a)',advance ='no')' GRD'
    if (pub_cube_format .and. pub_on_root) &
         write(stdout,'(a)',advance ='no')' CUBE'
    if (pub_dx_format .and. pub_on_root) &
         write(stdout,'(a)',advance ='no')' DX'
    if (pub_on_root) &
         write(stdout,'(a)',advance ='no')' ...'

    ! cks: Loop over NGWFs of pub_my_proc_id
    ngwf_counter =0

    do latom=1, par%num_atoms_on_proc(pub_my_proc_id)

       gatom = latom + par%first_atom_on_proc(pub_my_proc_id) - 1
       orig_atom =par%orig_atom(gatom)

       do at_func=1,ngwf_basis%num_on_atom(gatom)

          ngwf_counter= ngwf_counter+1

          ! cks: create plot files only if needed for current NGWF
          if (par%elements_on_proc(latom)%ngwf_plot) then


             ! cks: construct the title line for output file
             if (pub_nnho) then
                write(title_line,*)  trim(adjustl(ngwfs_head_name)) // &
                     ' NNHO number '
                write(fun_number,'(i2)') at_func
                write(temp_name,*) trim(adjustl(title_line)), trim(fun_number),&
                     ' of atom '
                write(fun_number,'(i5)') orig_atom
                write(title_line,*) trim(adjustl(temp_name)), trim(fun_number), &
                     ' (', trim(adjustl(par%elements_on_proc(latom)%symbol)) &
                     // ')' // ' for: "' // trim(pub_rootname) //'"'
             else
                write(title_line,*)  trim(adjustl(ngwfs_head_name)) // &
                     ' NGWF number '
                write(fun_number,'(i2)') at_func
                write(temp_name,*) trim(adjustl(title_line)), trim(fun_number),&
                     ' of atom '
                write(fun_number,'(i5)') orig_atom
                write(title_line,*) trim(adjustl(temp_name)), trim(fun_number), &
                     ' (', trim(adjustl(par%elements_on_proc(latom)%symbol)) &
                     // ')' // ' for: "' // trim(pub_rootname) //'"'
             endif

             title_line = adjustl(title_line)

             ! cks: ======= PUT NGWF IN TIGHTBOX ==============================

             ! cks: origin of tight_box wrt origin of the simulation cell
             ! cks: in integer number of grid points
             call basis_location_func_wrt_cell(tight_start1, &
                  tight_start2, tight_start3, &
                  ngwf_basis%tight_boxes(ngwf_counter), cell)

             ! cks: number of points in tight box of current function
             npt1 = ngwf_basis%tight_boxes(ngwf_counter)%tight_pts1
             npt2 = ngwf_basis%tight_boxes(ngwf_counter)%tight_pts2
             npt3 = ngwf_basis%tight_boxes(ngwf_counter)%tight_pts3

             ! cks: copy function from ppd representation to output_box 3-D array
             call basis_copy_function_to_box(output_box, 1, 1, 1, &
                  ngwf_basis%tight_boxes(ngwf_counter), ngwfs_on_grid, &
                  ngwf_basis%spheres(ngwf_counter), cell, fftbox)

             ! agrecocmplx: plot the abs of the NGWFs if they are complex
             if (loc_cmplx) then
                output_box_abs(:,:,:) = abs(output_box%z(:,:,:))
             end if
             !cks: ==== END PUT NGWF IN TIGHTBOX ==============================


             ! cks: name of output file without the extension
             write(temp_name,*) trim(adjustl(ngwfs_head_name))// &
                  trim('_atom')
             write(fun_number,'(i5.5)') orig_atom
             write(temp_name1,*) trim(adjustl(temp_name)) // &
                  trim(adjustl(fun_number)) //  '_ngwf'
             write(fun_number,'(i2.2)') at_func
             write(temp_name,*) trim(temp_name1), trim(adjustl(fun_number))

             ! jd, cks: ########### For CUBE and DX output ############
             if (pub_cube_format .or. pub_dx_format) then

                a1c=(real(npt1,kind=DP)/real(cell%total_pt1,kind=DP))&
                     *cell%a1
                a2c=(real(npt2,kind=DP)/real(cell%total_pt2,kind=DP))&
                     *cell%a2
                a3c=(real(npt3,kind=DP)/real(cell%total_pt3,kind=DP))&
                     *cell%a3
                origin_cube=(real(tight_start1-1, kind=DP)/ &
                          real(cell%total_pt1,kind=DP))*cell%a1 &
                          &+(real(tight_start2-1, kind=DP)/ &
                          real(cell%total_pt2,kind=DP))*cell%a2 &
                          &+(real(tight_start3-1, kind=DP)/ &
                          real(cell%total_pt3,kind=DP))*cell%a3
             end if


             ! cks: ######### CUBE OUTPUT ####################################
             if (pub_cube_format) then

                ! cks: add extension to output file name
                write(fun_name,*) trim(temp_name), '.cube'
                fun_name = adjustl(fun_name)

                ! cks: *** CUBE ***
                ! cks: output NGWF in .cube format
                ! agrecocmplx: complex case
                if (loc_cmplx) then
                   call visual_output_cube( &
                        output_box_abs, "Output from ONETEP calculation", &
                        trim(title_line), trim(fun_name), &
                        origin_cube%X,origin_cube%Y,origin_cube%Z, &
                        a1c, a2c, a3c, npt1, npt2, npt3, npt1, .true., elements)
                ! real case
                else
                   call visual_output_cube( &
                        output_box%d, "Output from ONETEP calculation", & !jmecmplx
                        trim(title_line), trim(fun_name), &
                        origin_cube%X,origin_cube%Y,origin_cube%Z, &
                        a1c, a2c, a3c, npt1, npt2, npt3, npt1, .true., elements)
                end if

             endif
             ! cks: ##### END CUBE OUTPUT ####################################

             ! jd:  ######### DX OUTPUT ######################################
             if (pub_dx_format) then

                ! jd: add extension to output file name
                write(fun_name,*) trim(temp_name), '.dx'
                fun_name = adjustl(fun_name)

                ! jd: *** DX ***
                ! jd: output NGWF in .dx format
                ! agrecocmplx: complex case
                if (loc_cmplx) then
                   call visual_output_dx( &
                        output_box_abs, "Output from ONETEP calculation", &
                        trim(title_line), trim(fun_name), &
                        origin_cube%X,origin_cube%Y,origin_cube%Z, &
                        a1c, a2c, a3c, npt1, npt2, npt3, npt1, .true., .true.)
                ! real case
                else
                   call visual_output_dx( &
                        output_box%d, "Output from ONETEP calculation", & !jmecmplx
                        trim(title_line), trim(fun_name), &
                        origin_cube%X,origin_cube%Y,origin_cube%Z, &
                        a1c, a2c, a3c, npt1, npt2, npt3, npt1, .true., .true.)
                end if

             endif
             ! jd:  ##### END DX OUTPUT ######################################

             ! cks: ######### GRD OUTPUT #####################################
             if (pub_grd_format) then

                ! cks: add extension to output file name
                write(fun_name,*) trim(temp_name), '.grd'
                fun_name = adjustl(fun_name)

                aa =(npt1*cell%d1)/ANGSTROM
                bb =(npt2*cell%d2)/ANGSTROM
                cc =(npt3*cell%d3)/ANGSTROM

                start1 =tight_start1 -1
                start2 =tight_start2 -1
                start3 =tight_start3 -1

                end1 =start1 +npt1
                end2 =start2 +npt2
                end3 =start3 +npt3

                ! cks: *** GRD ***
                ! cks: output NGWF in .grd format
                ! agrecocmplx: complex case
                if (loc_cmplx) then
                   call visual_output_grd( &
                        output_box_abs, aa, bb, cc, alpha, beta, gamma, &
                        npt1, npt2, npt3, npt3, start1, start2, start3, &
                        end1, end2, end3, trim(title_line), trim(fun_name), &
                        .false., .true.)   !input
                ! real case
                else
                   call visual_output_grd( &
                        output_box%d, aa, bb, cc, alpha, beta, gamma, & !jmecmplx
                        npt1, npt2, npt3, npt3, start1, start2, start3, &
                        end1, end2, end3, trim(title_line), trim(fun_name), &
                        .false., .true.)   !input
                end if

             endif
             ! cks: ##### END GRD OUTPUT #####################################

          endif

       enddo

    enddo

    call comms_barrier

    ! agrecocmplx
    if (loc_cmplx) then
       deallocate(output_box_abs,stat=ierr)
       call utils_dealloc_check('visual_ngwfs','output_box_abs', ierr)
    end if

    call data_fftbox_dealloc(output_box)

    ! cks: Anounce end of write operations
    if (pub_on_root) write(stdout,*)' done'

  end subroutine visual_ngwfs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine visual_efield_on_grid_pbc(efield_slabs12, pot_slabs12, grid, cell)
    !=========================================================================!
    ! This subroutine calculates the electric field on the grid, given the    !
    ! potential on the grid. The derivative is done in reciprocal space, so   !
    ! periodicity is assumed regardless of whether the calculation actually   !
    ! uses PBCs or not. Use visual_efield_on_grid_obc() for OBC calculations. !
    !                                                                         !
    ! The units of the field are a0^-1 times whatever the potential used.     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   efield_slabs12 (in/out): The calculated field is returned here, in    !
    !                            usual slab representation, except there are  !
    !                            three Cartesian components in the last index.!
    !                            Caller is responsible for allocation.        !
    !   pot_slabs12 (in): Potential generating the field.                     !
    !   grid (in): The grid on which the potential and field are defined.     !
    !   cell (in): Simulation cell descriptor.                                !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2016, using pseudo_local_calculate_   !
    ! forces() due to AAM as template.                                        !
    ! Fixed by Jacek Dziedzic on 2016.10.18 to correctly calculate Y and Z    !
    ! components, using xc_gradients as guideline.                            !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_proc_id
    use constants, only: DP, stdout, cmplx_i
    use fourier, only: fourier_apply_cell_forward, fourier_apply_cell_backward
    use rundat, only: pub_debug_on_root
!$  use rundat, only: pub_threads_max
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)  :: grid
    type(CELL_INFO), intent(in)  :: cell
    real(kind=DP), intent(inout) :: efield_slabs12(grid%ld1, &
         grid%ld2, grid%max_slabs12, 3)
    real(kind=DP), intent(in)    :: pot_slabs12(grid%ld1, &
         grid%ld2, grid%max_slabs12)

    ! Local Variables
    integer :: ipt, i2, i3, islab23, dir
    real(kind=DP) :: gvec(3)
    complex(kind=DP) :: v_of_g
    complex(kind=DP), allocatable :: recip(:,:,:,:)
    integer :: ierr
    character(len=*), parameter :: myself = 'visual_efield_on_grid_pbc'

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering '//myself

    call timer_clock(myself,1)

    ! Allocate workspace and clear output
    efield_slabs12 = 0D0
    allocate(recip(grid%ld3, grid%ld2, grid%max_slabs23,3),stat=ierr)
    call utils_alloc_check(myself,'recip',ierr)

    ! Fourier transform the potential to reciprocal space
    call fourier_apply_cell_forward(pot_slabs12(:,:,:), recip(:,:,:,1), grid)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(ipt,i2,i3,islab23,gvec,v_of_g) &
!$OMP SHARED (pub_my_proc_id,grid,recip,pub_threads_max)

    ! jd: Loop over all g points
!$OMP DO
    do ipt = 1, grid%num_slabs23 * grid%n2 * grid%n3

       i3 = modulo(ipt-1,grid%n3) + 1
       i2 = modulo((ipt-i3)/grid%n3,grid%n2) + 1
       islab23 = (ipt-(i2-1)*grid%n3-i3) / (grid%n3*grid%n2) + 1

       ! jd: Find g-vector
       call cell_grid_recip_pt(gvec, &
            islab23 + grid%first_slab23(pub_my_proc_id) - 1, i2, i3, grid)

       v_of_g = recip(i3,i2,islab23,1)

       ! jd: The grad operator is I*g, minus comes from E = -gradU.
       !     Results verified against numerical derivative of potential
       !     for all Cartesian components
       recip(i3,i2,islab23,1:3) = -CMPLX_I * gvec(1:3) * v_of_g

    end do
!$OMP END DO
!$OMP END PARALLEL

    ! jd: Fourier transform the electric field to back to real space
    do dir = 1,3
       call fourier_apply_cell_backward(efield_slabs12(:,:,:,dir), &
            recip(:,:,:,dir), grid)
    end do
    ! Deallocate workspace
    deallocate(recip,stat=ierr)
    call utils_dealloc_check(myself,'recip',ierr)

    call timer_clock(myself,2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving '//myself

  end subroutine visual_efield_on_grid_pbc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine visual_efield_on_grid_obc(efield_slabs12, pot_slabs12, grid, cell,&
       fd_order)
    !=========================================================================!
    ! This subroutine calculates the electric field on the grid, given the    !
    ! potential on the grid. The derivative is done using finite differences, !
    ! and OBCs are assumed regardless of whether the calculation actually     !
    ! uses OBCs or not. Use visual_efield_on_grid_pbc() for PBC calculations. !
    !                                                                         !
    ! The units of the field are a0^-1 times whatever the potential used.     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   efield_slabs12 (in/out): The calculated field is returned here, in    !
    !                            usual slab representation, except there are  !
    !                            three Cartesian components in the last index.!
    !                            Caller is responsible for allocation.        !
    !   pot_slabs12 (in): Potential generating the field.                     !
    !   grid (in): The grid on which the potential and field are defined.     !
    !   cell (in): Simulation cell descriptor.                                !
    !   fd_order (in): Order of FDs to use for the gradient calculation.      !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2016.                              !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: DP, stdout
    use finite_differences, only: FD_GRID_INFO, finite_difference_set_geometry,&
         finite_difference_gradient
    use rundat, only: pub_debug_on_root
    use services, only: services_1d_interpolation
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)  :: grid
    type(CELL_INFO), intent(in)  :: cell
    real(kind=DP), intent(inout) :: efield_slabs12(grid%ld1, &
         grid%ld2, grid%max_slabs12, 3)
    real(kind=DP), intent(in)    :: pot_slabs12(grid%ld1, &
         grid%ld2, grid%max_slabs12)
    integer, intent(in)          :: fd_order

    ! Local Variables
    type(FD_GRID_INFO) :: fd_grid
    character(len=*), parameter :: myself = 'visual_efield_on_grid_obc'

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering '//myself

    call timer_clock(myself,1)

    ! Allocate workspace and clear output
    efield_slabs12 = 0D0

    fd_grid = finite_difference_set_geometry(grid,1,(/.false., .false., .false./))
    call finite_difference_gradient(efield_slabs12, pot_slabs12, &
         fd_order, fd_grid, grid)

    ! jd: E = -grad U
    efield_slabs12 = -efield_slabs12

    call timer_clock(myself,2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving '//myself

  end subroutine visual_efield_on_grid_obc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine visual_output_cube( &
       cube_contents, comment_line1, comment_line2, &              ! input
       output_file_name, origin1, origin2, origin3, a1, a2, a3, &  ! input
       npt1, npt2, npt3, npt1_write, writeheader, elements)        ! input

    !=========================================================================!
    ! This subroutine outputs in CUBE format to a file, named according to    !
    ! the output_file_name argument, a 3D scalar function represented on a    !
    ! grid. The CUBE format also contains the geometry of the                 !
    ! molecule and it can be read and plotted by a variety of molecular       !
    ! visualisation programs including gOpenMol, Molekel and xCrysDen         !
    !-------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 12/3/2005 for the ONETEP            !
    ! linear-scaling DFT program.                                             !
    ! Modified by Nicholas Hine on 04/12/2007 to support non-orthogonal axes. !
    ! Removed dependency of elements on par%nat, Rob Bell 28/01/2014          !
    !=========================================================================!

    use constants, only: DP
    use geometry, only: point
    use ion, only: element
    use utils, only: utils_unit, utils_abort

    implicit none

    real(kind=DP), intent(in) :: cube_contents(:,:,:) ! buffer with output data
    character(len=*), intent(in) :: output_file_name ! name of file to write
    character(len=*), intent(in) :: comment_line1 ! first comment line to output
    character(len=*), intent(in) :: comment_line2 ! second comment line to output
    real(kind=DP), intent(in) :: origin1 ! coordinates of cube wrt origin
    real(kind=DP), intent(in) :: origin2 ! coordinates of cube wrt origin
    real(kind=DP), intent(in) :: origin3 ! coordinates of cube wrt origin
    TYPE(POINT), intent(in) :: a1, a2, a3 ! ndmh: cell axes rather than d1,d2,d3
    integer, intent(in) :: npt1, npt2, npt3 ! size of cube in grid points
    integer, intent(in) :: npt1_write ! how many of the slabs in a1 dir to output
    logical, intent(in) :: writeheader
    type(ELEMENT), intent(in), optional :: elements(:) ! all elements
    ! in non-spacefill order

    ! cks: << local variables >>
    integer :: atom, ii, remainder ! counting indices
    integer :: dir1, dir2, dir3    ! output loop indices
    integer :: output_unit         ! output unit number
    integer :: nat_out             ! number of atoms to write
    character(len=80) :: remchar   ! string buffer


    ! cks: get a unit number that is free
    output_unit = utils_unit()

    ! ndmh: write header before first slab only
    if(writeheader)then

       ! cks: open up file for output of function
       open(output_unit, file =output_file_name, err =10 )

       ! cks: Write the first two lines which are comments
       write(output_unit,*) comment_line1
       write(output_unit,*) comment_line2

       ! cks: This is now expressed in Cartesian coordinates and
       ! cks: currently it is ONLY CORRECT FOR A CUBIC GRID!
       ! cks: The fastest running direction is direction3 and
       ! cks: the slowest is direction1

       ! ndmh: extended to any grid type

       ! cks: number of atoms, position of the origin of simulation cell
       ! cks: w.r.t origin of simulation cell in cartesian coordinates
       if (present(elements)) then
          nat_out = size(elements)
       else
          nat_out = 1
       end if
       write(output_unit,'(i4,3F13.5)') nat_out, origin1, origin2, origin3

       ! cks: the number of points of the cube and its grid-spacing
       ! ndmh: modified to work with non-orthogonal cells
       write(output_unit,'(i4,3F13.5)') npt1, a1%X/real(npt1,kind=DP),&
             a1%Y/real(npt1,kind=DP),a1%Z/real(npt1,kind=DP)
       write(output_unit,'(i4,3F13.5)') npt2, a2%X/real(npt2,kind=DP),&
             a2%Y/real(npt2,kind=DP),a2%Z/real(npt2,kind=DP)
       write(output_unit,'(i4,3F13.5)') npt3, a3%X/real(npt3,kind=DP),&
             a3%Y/real(npt3,kind=DP),a3%Z/real(npt3,kind=DP)

       ! cks: atomic number, ps-ion-charge and cartesian coordinates of each atom
       if (present(elements)) then
          do atom=1,size(elements)
             write(output_unit,'(i4,4F13.5)') elements(atom)%atomic_number, &
                  real(elements(atom)%ion_charge, kind=DP),  &
                  elements(atom)%centre
          enddo
       else
          write(output_unit,'(i4,4F13.5)') 1, 1.0_DP, (/0.0_DP,0.0_DP,0.0_DP/)
       end if

    else

       ! ndmh: re-open file for continuation of output
       open(output_unit, file=output_file_name, status='old', &
            position='append',err =10 )

    end if

    remainder = mod(npt3,6)
    write(remchar,*) remainder

    ! cks: values of charge density on grid points of the fine grid sim cell
    ! cks: batches of up to six points are output along
    ! cks: the fastest varying direction (dir3).

    ! ndmh: only write out first 1,npt_write lines of cube_contents
    ! ndmh: to remove need to collate entire array at once
    do dir1=1,npt1_write
       do dir2=1,npt2

          do dir3=1,npt3-remainder,6
             write(output_unit, '(6E13.5)') (cube_contents(dir1, dir2, dir3+ii),&
                  ii =0,5  )
          enddo

          if (remainder.gt.0) then
             write(output_unit, '('//trim(adjustl(remchar))//'E13.5)') &
                  (cube_contents(dir1, dir2, npt3 -remainder +ii), ii =1, &
                  remainder  )
          endif

       enddo
    enddo

    close( output_unit, err =20)

    return

10  call utils_abort('Error during opening of file:'//trim(output_file_name)//&
         ' in visual_output_cube().')
20  call utils_abort('Error during closing of file:'//trim(output_file_name)//&
         ' in visual_output_cube().')

  end subroutine visual_output_cube


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine visual_output_grd( &
       grd_contents, aa, bb, cc, alpha, beta, gamma, &         ! input
       npt1, npt2, npt3, npt3_write, start1, start2, start3, & ! input
       end1, end2, end3, title_line, output_file_name, &       ! input
       grd_2in1, writeheader)                                  ! input

    !=======================================================================!
    ! This subroutine outputs in GRD format to a file, named according to   !
    ! the output_file_name argument, a 3D scalar function represented on a  !
    ! grid. The GRD format does not contain the atomic coordinates but      !
    ! supports arbitrary bravais lattices.                                  !
    !-----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 19/6/2005 for the ONETEP          !
    ! linear-scaling DFT program.                                           !
    !=======================================================================!

    use constants, only: DP
    use utils, only: utils_unit, utils_abort

    implicit none

    real(kind=DP), intent(in) :: grd_contents(:, :, :) ! data to output
    real(kind=DP), intent(in) :: aa, bb, cc ! lattice parameters
    real(kind=DP), intent(in) :: alpha, beta, gamma ! lattice angles
    character(len=*), intent(in) :: output_file_name ! name of file to write
    character(len=*), intent(in) :: title_line ! comment line
    integer, intent(in) :: npt1, npt2, npt3 !grid point num along lattice vectors
    integer, intent(in) :: npt3_write       !number of grid pts along a3 to write
    integer, intent(in) :: start1, start2, start3 ! start of tightbox from origin
    integer, intent(in) :: end1, end2, end3 ! end of tightbox from origin
    logical, intent(in) :: grd_2in1
    logical, intent(in) :: writeheader

    ! cks: << local variables >>
    integer :: dir1, dir2, dir3  ! output loop indices
    integer :: output_unit       ! output unit number
    integer :: mm1  ! output array index
    integer :: mm2  ! output array index
    integer :: mm3  ! output array index

    ! cks: get a unit number that is free
    output_unit = utils_unit()

    ! ndmh: write header if this is first block
    if (writeheader) then

       ! cks: open up file for output of function
       open(output_unit, file =output_file_name, err =30)

       ! cks: lattice info
       write(output_unit,*)trim(adjustl(title_line))
       write(output_unit,*)'(1p,e12.5)'
       write(output_unit,'(6f8.3)') aa, bb, cc, alpha, beta, gamma

       ! cks: grid info
       if (.not. grd_2in1) then
          write(output_unit,'(3i8)') npt1, npt2, npt3
       else
          write(output_unit,'(3i8)') npt1/2, npt2/2, npt3/2
       end if
       write(output_unit,'(7i8)') 1, start1, end1, &
            start2, end2, start3, end3

    else

       ! ndmh: re-open file for continuation of output
       open(output_unit, file =output_file_name, status='old', &
            position='append',err =30 )
    end if

    ! ndmh: combined these two routines into one
    if (.not.grd_2in1) then

       ! cks: output plot data
       do dir3 =1, npt3_write +1
          mm3 =dir3
          if (dir3 == npt3 +1) mm3 =1

          do dir2 = 1, npt2 +1
             mm2 =dir2
             if (dir2 == npt2 +1) mm2 =1

             do dir1 =1, npt1 +1
                mm1 =dir1
                if (dir1 == npt1 +1) mm1 =1

                write(output_unit, '(1p,e12.5)') grd_contents(mm1, mm2, mm3)

             enddo
          enddo
       enddo

    else

       ! cks: output plot data, only every second point in each direction
       do dir3 =1, npt3_write, 2
          do dir2 = 1, npt2 +1, 2
             do dir1 =1, npt1 +1, 2
                 write(output_unit, '(1p,e12.5)') &
                      grd_contents(mod(dir1,npt1), mod(dir2,npt2), mod(dir3,npt3))
             enddo
          enddo
       enddo

    end if

    close( output_unit, err =40)

    return

30  call utils_abort('Error during opening of file:'//trim(output_file_name)//&
         ' in visual_output_grd().')
40  call utils_abort('Error during opening of file:'//trim(output_file_name)//&
         ' in visual_output_grd().')

  end subroutine visual_output_grd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine visual_output_dx( &
       cube_contents, comment_line1, comment_line2, &    ! input
       output_file_name, origin1, origin2, origin3, a1, a2, a3, &  ! input
       npt1, npt2, npt3, npt1_write, writeheader, writetrailer)    ! input

    !=========================================================================!
    ! This subroutine outputs in OpenDX format to a file, named according to  !
    ! the output_file_name argument, a 3D scalar function represented on a    !
    ! grid. Operation is similar to visual_output_cube, except this routine   !
    ! also allows for outputting a coarse representation and to tune the      !
    ! number of significant digits.                                           !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 08/03/2010, using visual_output_cube as    !
    ! a template.                                                             !
    !=========================================================================!

    use constants, only: DP, ANGSTROM
    use geometry, only: point
    use ion, only: element
    use rundat, only: pub_dx_coarse, pub_dx_sig_digits
    use utils, only: utils_unit, utils_abort

    implicit none

    real(kind=DP), intent(in) :: cube_contents(:,:,:) ! buffer with output data
    character(len=*), intent(in) :: output_file_name ! name of file to write
    character(len=*), intent(in) :: comment_line1 ! first comment line
    character(len=*), intent(in) :: comment_line2 ! second comment line
    real(kind=DP), intent(in) :: origin1 ! coordinates of cube wrt origin
    real(kind=DP), intent(in) :: origin2 ! coordinates of cube wrt origin
    real(kind=DP), intent(in) :: origin3 ! coordinates of cube wrt origin
    TYPE(POINT), intent(in) :: a1, a2, a3 ! ndmh: cell axes rather than d1,d2,d3
    integer, intent(in) :: npt1, npt2, npt3 ! size of cube in grid points
    integer, intent(in) :: npt1_write   ! how many of the slabs in a1 dir
    logical, intent(in) :: writeheader  ! should the file header be written?
    logical, intent(in) :: writetrailer ! should the file trailer be written?

    ! cks,jd : << local variables >>
    integer :: dir1, dir2, dir3    ! output loop indices
    integer, save :: output_unit         ! output unit number
    integer :: ndata  ! jd: number of output data elements
    integer :: stride ! jd: 2, if dx_format_coarse, 1 otherwise (the usual case)
    integer, save :: column ! jd: current column of .dx output, must be stored
                            !     across batches (!), but not across files
    real(kind=DP) :: curval ! jd: value at the point currently writtn out
    character(len=80) :: format_str, format_str1, format_str2 ! jd: temporaries
    integer :: npt1_out, npt2_out, npt3_out ! jd: actual number of pts to output


    ! jd: determine the stride (1 for fine grid, 2 for coarse grid)
    stride=1
    if(pub_dx_coarse) stride=2

    ! jd: prepare the formatting string, depending on pub_dx_sig_digits
    write(format_str1,*) pub_dx_sig_digits+7
    write(format_str2,*) pub_dx_sig_digits-1
    format_str = &
         'ES'//trim(adjustl(format_str1))//'.'//trim(adjustl(format_str2))

    ! jd: determine how many points along each direction we output
    npt1_out = (npt1/stride + mod(npt1,stride))
    npt2_out = (npt2/stride + mod(npt2,stride))
    npt3_out = (npt3/stride + mod(npt3,stride))

    ! ndmh: write header before first slab only
    if(writeheader) then

       ! jd: Since this is the first batch, reset column
       column = 1

       ! cks: get a unit number that is free
       output_unit = utils_unit()

       ! cks: open up file for output of function
       open(output_unit, file =output_file_name, err =10 )

       ! cks: Write the first two lines which are comments
       ! jd: Truncate comment line at char 80, or else it overflows onto the
       !     next line. Can't keep it on the same line either, as long comment
       !     lines crash VMD.
       write(output_unit,'(a)',advance='no') '# '
       write(output_unit,*) comment_line1(:min(len_trim(comment_line1),77))
       write(output_unit,'(a)',advance='no') '# '
       write(output_unit,*) comment_line2(:min(len_trim(comment_line2),77))

       ! jd: Write the number of ticks along each axis
       !     We don't want padding or risk having asterisks,
       !     hence the formatting mess
       write(output_unit,'(a,i0,a,i0,a,i0)') &
            'object 1 class gridpositions counts ', &
            npt1_out,' ',npt2_out,' ',npt3_out

       ! jd: Write the box origin
       write(output_unit,'(a,3'//trim(format_str)//')') 'origin', &
            origin1/ANGSTROM, origin2/ANGSTROM, origin3/ANGSTROM

       ! jd: Write the box axes deltas
       write(output_unit,'(a,3'//trim(format_str)//')') 'delta', &
            a1%X*real(stride,kind=DP)/real(npt1,kind=DP)/ANGSTROM, &
            a1%Y*real(stride,kind=DP)/real(npt1,kind=DP)/ANGSTROM, &
            a1%Z*real(stride,kind=DP)/real(npt1,kind=DP)/ANGSTROM
       write(output_unit,'(a,3'//trim(format_str)//')') 'delta', &
            a2%X*real(stride,kind=DP)/real(npt2,kind=DP)/ANGSTROM, &
            a2%Y*real(stride,kind=DP)/real(npt2,kind=DP)/ANGSTROM, &
            a2%Z*real(stride,kind=DP)/real(npt2,kind=DP)
       write(output_unit,'(a,3'//trim(format_str)//')') 'delta', &
            a3%X*real(stride,kind=DP)/real(npt3,kind=DP)/ANGSTROM, &
            a3%Y*real(stride,kind=DP)/real(npt3,kind=DP)/ANGSTROM, &
            a3%Z*real(stride,kind=DP)/real(npt3,kind=DP)/ANGSTROM

       ! jd: Write out the number of counts along each axis, again
       write(output_unit,'(a,i0,a,i0,a,i0)') &
            'object 2 class gridconnections counts ', &
            npt1_out,' ',npt2_out,' ',npt3_out

       ! jd: Calculate and write out the number of elements that follows
       ndata = npt1_out * npt2_out * npt3_out
       write(output_unit,'(a,i0,a)') &
            'object 3 class array type double rank 0 items ',ndata, &
            ' data follows'

    end if

    ! jd: Write out in batches, the same way that visual_output_cube does,
    !     except allow for a stride different from 1.
    !     NB: When npt1_write /= npt1 (i.e. for large files in parallel mode)
    !         npt1_write _must_ be even when stride = 2 is used, otherwise
    !         fencepost errors will occur on batch boundaries, and some points
    !         will be skipped. Currently this is not an issue, since the
    !         evenness of npt1_write is ensured in visual_scalarfield, and
    !         visual_ngwfs does not use batching and visual_output_dx is not
    !         called from anywhere else. Tread carefully, however, as this
    !         might become an issue in the future, if stride>2 is ever used
    !         or batching is introduced in visual_ngwfs or visual_scalarfield is
    !         modified or visual_output_dx routine is called from somewhere
    !         else. Thus the sanity check below.
    if(mod(npt1_write,2) /= 0 .and. npt1_write /= npt1) then
       call utils_abort('Sanity check failed in visual_output_dx.')
    end if

    do dir1 =1, npt1_write, stride
       do dir2 = 1, npt2, stride
          do dir3 = 1, npt3, stride
             curval = cube_contents(dir1, dir2, dir3)

             ! jd: Avoid writing out values with exponents that do not fit in
             !     two digits, because then the resultant file cannot be parsed
             !     correctly. Values with exponents smaller than -99 will be
             !     chopped to zero. Values with exponents larger than 98 would
             !     also be chopped to zero, but we do not deal with such large
             !     quantities in onetep.
             if ((curval < 1.0D-99 .and. curval > -1.0D-99) .or. &
                  curval > 1.0D99 .or. curval < -1.0D99) then
                  curval = 0.0_DP
             end if

             write(output_unit, '('//trim(format_str)//')',advance='no') curval

             column = column + 1
             if(column > 3) then ! jd: Output a newline every 3 writes
                write(output_unit,'(a)') ''
                column = 1
             end if
          end do
       end do
    end do


    if(writetrailer) then

       ! jd: Output a newline if last row was not full
       if(column /= 1) write(output_unit,'(a)') ''

       write(output_unit,'(a)') 'attribute "dep" string "positions"'
       write(output_unit,'(a)') 'object "regular positions regular connections"&
            & class field'
       write(output_unit,'(a)') 'component "positions" value 1'
       write(output_unit,'(a)') 'component "connections" value 2'
       write(output_unit,'(a)') 'component "data" value 3'
       write(output_unit,'(a)') ''
       write(output_unit,'(a)') 'end'

       close( output_unit, err =20)

    end if


    return

10  call utils_abort('Error during opening of file:'//trim(output_file_name)//&
         ' in visual_output_dx().')

20  call utils_abort('Error during closing of file:'//trim(output_file_name)//&
         ' in visual_output_dx().')

  end subroutine visual_output_dx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine visual_scalarfield_ngwf_basis(scalarfield_in_ngwf_basis, & ! input
       rep, mdl, ireg, ngwf_basis, file_header, scalar_name, &          ! input
       conversion_factor)                                     ! input, optional
    !=====================================================================!
    ! This subroutine outputs a real scalar field (eg. potential or den-  !
    ! sity) to a file in a format suitable for visualisation. The scalar  !
    ! field V is in the NGWF basis representation. The real-space repre-  !
    ! sentation is reconstructed from |phi_a> | S^-1 V S^-1 | <phi_b|     !
    ! using density_on_grid().                                            !
    !---------------------------------------------------------------------!
    ! Arguments:                                                          !
    !   scalarfield_in_ngwf_basis (in): The scalarfield in matrix rep that!
    !                                   is to be visualised.              !
    !   rep, mdl, ngwf_basis (in): The usual.                             !
    !   file_header, scalar_name, conversion_factor (in): Will be passed  !
    !                                              to visual_scalarfield. !
    !---------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2017/04.                               !
    ! Modified for embedding by Joseph Prentice, September 2018           !
    !=====================================================================!

    use constants, only: DP
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_num_spins
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_product
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! jd: Arguments
    type(SPAM3), intent(in)       :: scalarfield_in_ngwf_basis
    type(NGWF_REP), intent(in) :: rep
    type(MODEL), intent(in)       :: mdl
    type(FUNC_BASIS), intent(in)  :: ngwf_basis
    character(len=*), intent(in)  :: file_header ! part of header line of output
    character(len=*), intent(in)  :: scalar_name ! part of file-name to output
    integer, intent(in)           :: ireg
    real(kind=DP), intent(in), optional :: conversion_factor

    ! jd: Local variables
    integer                          :: pub_num_spins_backup
    type(SPAM3)                      :: VSinv, V_kernel(1)
    real(kind=DP), allocatable       :: V_fine(:,:,:,:)
    integer                          :: ierr
    character(len=*), parameter      :: myself = 'visual_scalarfield_ngwf_basis'

    ! -------------------------------------------------------------------------

    ! jd: Construct V_kernel = S^-1 V S^-1
    call sparse_create(VSinv, rep%inv_overlap%m(ireg,ireg))
    call sparse_create(V_kernel(1), rep%inv_overlap%m(ireg,ireg))
    call sparse_product(VSinv, scalarfield_in_ngwf_basis, rep%inv_overlap%m(ireg,ireg))
    call sparse_product(V_kernel(1), rep%inv_overlap%m(ireg,ireg), VSinv)
    call sparse_destroy(VSinv)

    allocate(V_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12, 1),stat=ierr)
    call utils_alloc_check(myself,'V_fine',ierr)

    ! jd: Leverage density_on_grid to get V(r) = |phi_a> | S^-1 V S^-1 | <phi_b|
    pub_num_spins_backup = pub_num_spins
    pub_num_spins = 1
    call density_on_grid(V_fine, mdl%fine_grid, mdl%dbl_grid, &
         mdl%cell, mdl%fftbox, V_kernel, rep%ngwf_overlap%m(ireg,ireg), &
         rep%ngwfs_on_grid(ireg), ngwf_basis, &
         rep%ngwfs_on_grid(ireg), ngwf_basis)
    pub_num_spins = pub_num_spins_backup

    ! jd: Free V_kernel, output the reconstructed potential
    call sparse_destroy(V_kernel(1))

    call visual_scalarfield(V_fine, mdl%fine_grid, mdl%cell, &
         file_header, scalar_name, conversion_factor = conversion_factor)

    deallocate(V_fine,stat=ierr)
    call utils_dealloc_check(myself,'V_fine',ierr)

  end subroutine visual_scalarfield_ngwf_basis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine visual_scalarfield_read( &
       scalarfield_fine, grid, fd_grid, cell, scalar_name, &  ! input
       conversion_factor) ! input, optional

    !=====================================================================!
    ! This subroutine reads a real scalar field (e.g. density,            !
    ! potential, real molecular orbital, etc ) from a .cube file          !
    !---------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu, adapted from                     !
    ! visual_scalarfield                                                  !
    !=====================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_deposit_box
    use comms, only: pub_on_root, comms_barrier, pub_my_proc_id
    use constants, only: DP, LONG, ANGSTROM, PI, stdout
    use finite_differences, only: FD_GRID_INFO
    use geometry, only: geometry_magnitude, operator(.DOT.), POINT, operator(-)
    use ion, only: element
    use rundat, only: pub_rootname, pub_cube_format, pub_dx_format, &
         pub_grd_format, pub_debug_on_root
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_unit

    implicit none

    ! cks: subroutine arguments
    type(GRID_INFO), intent(in)     :: grid
    type(FD_GRID_INFO), intent(in)  :: fd_grid
    type(CELL_INFO), intent(in)     :: cell
    real(kind=DP), intent(inout)    :: scalarfield_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12) ! scalar field to read in
    character(len =*), intent(in) :: scalar_name ! part of file-name to read
    real(kind=DP), intent(in), optional :: conversion_factor

    ! Local variables
    real(kind=DP), dimension(:, :, :), allocatable :: buffer_fine
    ! comms buffer
    real(kind=DP), dimension(:, :, :), allocatable :: scalarfield_box_fine
    ! simcell buffer
    character(len =256) :: input_file  ! file names
    integer :: ldf1,ldf2,nf1,nf2,nf3    ! shorthand working variables
    integer :: ierr                     ! memory alloc error flag
    integer :: nstep                    ! number of steps to input array over
    integer :: istep       ! step counter for splitting output of arrays
    integer :: nslab1      ! number of slabs to read on this step (cube input)
    integer :: nslab1_full ! leading dim of scalarfield_box_fine (grd input)
    integer :: islab1      ! current slab index in '1' direction
    integer :: nslab3      ! number of slabs to read on this step
    integer :: islab3      ! current slab index in '3' direction

    ! ndmh: variables using long integers to avoid overflows
    ! maximum size of buffer on master proc (125000000 real*8's ~ 1GB)
    integer(kind=LONG),parameter :: max_size = 125000000_LONG
    integer(kind=LONG) :: total_size  ! total size of simulation cell data array
    integer(kind=LONG) :: nf1_long      ! long version of nf1
    integer(kind=LONG) :: nf3_long      ! long version of nf3

    ! jd: variables distinguishing CUBE and DX formats
    character(len=4) :: suffix ! 'cube' or 'dx'

    type(POINT) :: d1, d2, d3, d1_inp, d2_inp, d3_inp
    integer :: npt1_read, npt2_read, npt3_read ! size of cube in grid points
    integer :: nat_in
    real(kind=DP) :: origin1, origin2, origin3
    integer :: atom, ii, remainder ! counting indices
    integer :: dir1, dir2, dir3    ! input loop indices
    integer :: input_unit         ! input unit number
    character(len=80) :: remchar   ! string buffer
    integer :: buf_1
    real(kind=DP) :: buf_2, buf_3(3)
    character(len=256) :: comment_line1 ! first comment line
    character(len=256) :: comment_line2 ! second comment line


    integer :: i1_pq, i2_pq, i3_pq, loc_i3_pq ! jd: Bounds of generic indices

    ! cks: initialise shorthand variables
    ldf1  = fd_grid%ld1f
    ldf2  = fd_grid%ld2f
    nf1   = fd_grid%pq1f
    nf2   = fd_grid%pq2f
    nf3   = fd_grid%pq3f

    ! NOTE that ldf1 and ldf2 are actually equal to grid%ld1 and grid%ld2
    ! it is nf1, nf2 and nf3 that are adapted to granularity

    ! cks: *** CUBE format ***
    suffix = 'cube'

    ! ndmh: find input file name
    if (pub_on_root) then
       write(input_file,'(a256)') trim(pub_rootname)//trim(&
            adjustl(scalar_name))//'.'//suffix
       input_file = adjustl(input_file)

       write(stdout,'(3a)',advance ='no') &
            'Reading "', trim(input_file),'" ...'
    end if

    ! ndmh: calculate how many steps the input is divided into
    total_size = nf1 * nf2 * nf3
    if(total_size < max_size)then
       nstep = 1
       nslab1 = nf1
    else
       nf1_long = nf1
       nslab1 = int(nf1_long * max_size / total_size)
       if (nslab1==0) nslab1 = 1
       nslab1 = nslab1 + mod(nslab1,2) ! ensure it is even
       nstep = nf1 / nslab1
       if (nstep*nslab1 < nf1) nstep = nstep + 1
    end if

    if (pub_debug_on_root) then
       write(stdout,*)
       write(stdout,'(4a,i3,a,i3,a)') 'DEBUG: Reading ',suffix,' format ', &
             'simulation cell data in ', nstep, ' steps of ', nslab1, &
             ' slabs'
    endif

    ! ndmh: allocate simulation cell charge density buffer
    allocate(buffer_fine(nslab1,nf2,grid%max_slabs12), stat=ierr)
    call utils_alloc_check('visual_scalarfield_read','buffer_fine',ierr)
    allocate(scalarfield_box_fine(nslab1,nf2,nf3),stat=ierr)
    call utils_alloc_check('visual_scalarfield_read', &
         'scalarfield_box_fine',ierr)

    ! ndmh: read header before first slab only
    nslab1_full = nslab1
    islab1 = 1

    scalarfield_fine = 0.0_DP

    ! ==== read header =======================================================
    if (pub_on_root) then
       ! get a unit number that is free
       input_unit = utils_unit()

       ! open up file
       open(input_unit, file = trim(input_file), iostat=ierr)
       call utils_assert(ierr == 0, "Error in visual_scalarfield_read: &
            &opening file "//trim(input_file)//" failed with code ", ierr)

       ! cks: read the first two lines which are comments
       read(input_unit,*) comment_line1
       read(input_unit,*) comment_line2

       ! cks: This is now expressed in Cartesian coordinates and
       ! cks: currently it is ONLY CORRECT FOR A CUBIC GRID!
       ! cks: The fastest running direction is direction3 and
       ! cks: the slowest is direction1

       ! number of atoms, position of the origin of simulation cell
       ! w.r.t origin of simulation cell in cartesian coordinates
       read(input_unit,*) nat_in, origin1, origin2, origin3

       ! cks: the number of points of the cube and its grid-spacing
       ! ndmh: modified to work with non-orthogonal cells
       read(input_unit,*) npt1_read, d1%X,d1%Y,d1%Z
       read(input_unit,*) npt2_read, d2%X,d2%Y,d2%Z
       read(input_unit,*) npt3_read, d3%X,d3%Y,d3%Z

       ! see if grid spacing and npts match what we have in ONETEP
       call utils_assert(abs(geometry_magnitude(d1)-fd_grid%d1f) < 0.0001_DP, &
            "Error in visual_scalarfield_read: cube file has wrong spacing along a1")
       call utils_assert(abs(geometry_magnitude(d2)-fd_grid%d2f) < 0.0001_DP, &
            "Error in visual_scalarfield_read: cube file has wrong spacing along a2")
       call utils_assert(abs(geometry_magnitude(d3)-fd_grid%d3f) < 0.0001_DP, &
            "Error in visual_scalarfield_read: cube file has wrong spacing along a3")
       call utils_assert(npt1_read == nf1, "Error in visual_scalarfield_read: &
            &cube file has wrong number of points along a1")
       call utils_assert(npt2_read == nf2, "Error in visual_scalarfield_read: &
            &cube file has wrong number of points along a2")
       call utils_assert(npt3_read == nf3, "Error in visual_scalarfield_read: &
            &cube file has wrong number of points along a3")

       ! there will be no atoms present, just a series of numbers
       ! in one line, such as 1 1.0 0.0 0.0 0.0
       read(input_unit,*) buf_1, buf_2, buf_3(1:3)
    end if ! if root

    ! ========= end of reading header ========================================

    call comms_barrier

    do istep = 1, nstep

       ! ndmh: do not read beyond end of grid on last step
       if (islab1 + nslab1 - 1 > nf1) nslab1 = nf1 - islab1 + 1

       scalarfield_box_fine = 0.0_DP

       ! ndmh: read in batch of '23'-slabs
       if (pub_on_root) then

          remainder = mod(nf3,6)
          write(remchar,*) remainder

          ! values of charge density on grid points of the fine grid sim cell
          ! batches of up to six points are output along
          ! the fastest varying direction (dir3).
          ! only read first 1,npt_write lines of scalarfield_box_fine
          ! to remove need to collate entire array at once

          do dir1=1,nslab1
             do dir2=1,nf2

                do dir3=1,nf3-remainder,6
                   read(input_unit,*) &
                        (scalarfield_box_fine(dir1, dir2, dir3+ii), ii =0,5  )
                enddo

                if (remainder.gt.0) then
                   read(input_unit,*) &
                        (scalarfield_box_fine(dir1, dir2, nf3 -remainder +ii),&
                        ii =1, remainder)
                endif

             enddo
          enddo

          ! jd: Apply unit conversion factor, if requested
          if (present(conversion_factor)) then
             scalarfield_box_fine = scalarfield_box_fine * conversion_factor
          end if

       end if

       ! ndmh: deposit the box of nslab1 '23'-slabs of cell data
       call cell_grid_deposit_box(scalarfield_fine,scalarfield_box_fine, &
            buffer_fine, grid, nslab1, nf2, nf3, nslab1_full, nf2, islab1, &
            1, 1, pub_on_root, .false.)

       islab1 = islab1 + nslab1

    end do

    if (pub_on_root) close(input_unit, iostat=ierr)
    call utils_assert(ierr == 0, "Error in visual_scalarfield_read: &
         &closing file "//trim(input_file)//" failed with code ", ierr)
    call comms_barrier

    ! ndmh: deallocate buffers
    deallocate(scalarfield_box_fine,stat=ierr)
    call utils_dealloc_check('visual_scalarfield_read', &
         'scalarfield_box_fine',ierr)
    deallocate(buffer_fine,stat=ierr)
    call utils_dealloc_check('visual_scalarfield_read', &
         'buffer_fine',ierr)

    if (pub_on_root) write(stdout,'(a)') ' done'

    call comms_barrier

    return

  end subroutine visual_scalarfield_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module visual




