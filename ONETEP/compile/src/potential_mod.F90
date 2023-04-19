! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Chris-Kriton Skylaris
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!   with subsequent optimisations by Nicholas D.M. Hine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module potential

  use rundat, only: pub_debug_on_root

  implicit none

  private

  public :: potential_apply_to_ngwf_batch
  public :: potential_apply_to_fftbox_batch
  public :: potential_apply_to_ppd_funcs
  public :: potential_sawtooth_efield
  public :: potential_add_efield_ion_energy
  public :: potential_input_to_workspace
  public :: potential_confine_ppd_funcs ! gcc32
  public :: potential_confine_ngwf_batch ! gcc32
#ifdef GPU_PGI
  public :: potential_gpu_app2_ngwf_batch
#endif

contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine potential_confine_ngwf_batch(potential_fftbox_batch, &   ! output
       ngwfs_on_grid, ngwf_basis, grid, fftbox, cell, & ! input
       ket_start_in_box, ket_box_start, batch_size, local_start, &     ! input
       local_end, max_current_size)

    !==========================================================================!
    ! This subroutine applies the local confinement operator to a batch        !
    ! of NGWFs in fftboxes.                                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! potential_fftbox_batch (output): Batch off fftboxes with potential-NGWFs !
    ! ngwfs_on_grid (input) : NGWFs on this proc in ppd representation         !
    ! ngwf_basis (input)    : Function basis type for NGWFs                    !
    ! batch_size (input)    : Number of fftboxes in the batch                  !
    ! local_start(input)    : First NGWF of the batch                          !
    ! local_end  (input)    : Last NGWF of the batch                           !
    ! ket_start_in_box (input) :
    !--------------------------------------------------------------------------!
    ! Derived from potential_apply_to_ngwf_batch by Gabriel Constantinescu in  !
    ! July 2015                                                                !
    !==========================================================================!



    use basis, only: basis_copy_function_to_box, basis_box_origin_to_atom
    use cell_grid, only: GRID_INFO, cell_grid_extract_box
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_set_to_zero, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_filter, fourier_interpolate
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, OPERATOR(+), OPERATOR(-), &
         OPERATOR(*), geometry_magnitude
    use rundat, only: pub_threads_num_fftboxes, pub_threads_fftbox, &
         pub_threads_per_fftbox, pub_confined_ngwfs_barrier
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_use_var
    use fourier, only: fourier_filter, fourier_interpolate

    implicit none

    ! Arguments
    integer, intent(in) :: batch_size, local_start, local_end
    integer, intent(in) :: max_current_size ! max size of batch over all procs
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FFTBOX_DATA), intent(inout) :: potential_fftbox_batch(batch_size)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    type(GRID_INFO), intent(in)  :: grid
    integer, intent(in) :: ket_start_in_box(3,batch_size)
    integer, intent(in) :: ket_box_start(3,batch_size)

    ! Local Variables
    type(FFTBOX_DATA) :: row1_box, row2_box
    real(kind=DP), dimension(:,:,:), allocatable :: potential_box_dbl
    type(FFTBOX_DATA) :: row1_box_dbl, row2_box_dbl

    complex(kind=DP), allocatable, dimension(:,:,:) :: coarse_work,fine_work

    integer :: ierr
    integer :: row1, row2
    integer :: n1, n2, n3, ld1, ld2
    integer :: n1_dbl, n2_dbl, n3_dbl, ld1_dbl, ld2_dbl
    integer :: iii
    integer :: batch_count
    integer :: fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl

    logical :: i_need_potential ! whether pub_my_proc_id needs pot'l in fftbox

    integer :: i1,i2,i3
    real(kind=DP) :: magnitude
    type(POINT) :: box_origin,r_pos
    type(POINT) :: box_to_atom
    real(kind=DP) :: width
    ! agrecocmplx
    logical :: loc_cmplx

    ! -------------------------------------------------------------------------


    call utils_use_var(pub_threads_fftbox)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering potential_confine_ngwf_batch'

    ! Start timer
    call timer_clock('potential_confine_ngwf_batch', 1)

    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2
    n1_dbl = fftbox%total_pt1_dbl
    n2_dbl = fftbox%total_pt2_dbl
    n3_dbl = fftbox%total_pt3_dbl
    ld1_dbl = fftbox%total_ld1_dbl
    ld2_dbl = fftbox%total_ld2_dbl

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid%iscmplx

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(potential_box_dbl, iii, batch_count, row1, row2, ierr, &
!$OMP      i_need_potential, fftbox_start1_dbl, fftbox_start2_dbl, &
!$OMP      fftbox_start3_dbl, width, i1,i2,i3,&
!$OMP      coarse_work, fine_work, box_origin,box_to_atom,r_pos,magnitude) &
!$OMP FIRSTPRIVATE(row1_box,row2_box,row1_box_dbl,row2_box_dbl) &
!$OMP SHARED(potential_fftbox_batch, ket_start_in_box, ket_box_start, &
!$OMP      ngwfs_on_grid, n1, n2, n3, ld1, ld2, n1_dbl, &
!$OMP      n2_dbl, n3_dbl, ld1_dbl, ld2_dbl, local_start, local_end, &
!$OMP      cell, fftbox, grid, ngwf_basis, max_current_size, &
!$OMP      pub_threads_per_fftbox,pub_threads_num_fftboxes, &
!$OMP      pub_confined_ngwfs_barrier, loc_cmplx)

    ! ndmh: initialisations to prevent compiler warnings
    fftbox_start1_dbl = -1111111
    fftbox_start2_dbl = -2222222
    fftbox_start3_dbl = -3333333

    ! Allocate workspace arrays
    call data_fftbox_alloc(row1_box, ld1, ld2, n3, &
         iscmplx=loc_cmplx)
    call data_fftbox_alloc(row2_box, ld1, ld2, n3, &
         iscmplx=loc_cmplx)
    allocate(potential_box_dbl(ld1_dbl, ld2_dbl, 2*n3), stat=ierr)
    call utils_alloc_check('potential_confine_to_ngwf_batch',&
         'potential_box_dbl',ierr)
    call data_fftbox_alloc(row1_box_dbl, ld1_dbl, ld2_dbl, 2*n3, &
         iscmplx=loc_cmplx)
    call data_fftbox_alloc(row2_box_dbl, ld1_dbl, ld2_dbl, 2*n3, &
         iscmplx=loc_cmplx)
    allocate(coarse_work(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('potential_confine_ngwf_batch','coarse_work',ierr)
    allocate(fine_work(ld1*2, ld2*2, n3*2),stat=ierr)
    call utils_alloc_check('potential_confine_ngwf_batch','fine_work',ierr)

    ! cks: initialisation
    call data_set_to_zero(row1_box_dbl)
    call data_set_to_zero(row2_box_dbl)
    potential_box_dbl = 0.0_DP


!$OMP DO
    do iii=local_start, local_start+max_current_size-1, 2

       row1 = iii ; row2 = iii + 1

       if (row1 <= local_end) then

          batch_count = iii-local_start+1

          ! ndmh: use new basis function copying straight to fftbox from ppds
          call basis_copy_function_to_box(row1_box, &
               ket_start_in_box(1,batch_count), &
               ket_start_in_box(2,batch_count), &
               ket_start_in_box(3,batch_count), &
               ngwf_basis%tight_boxes(row1), ngwfs_on_grid, &
               ngwf_basis%spheres(row1), cell, fftbox)

       else
          call data_set_to_zero(row1_box)
       end if

       if (row2 <= local_end) then

          batch_count = iii-local_start+2

          ! ndmh: use new basis function copying straight to fftbox from ppds
          call basis_copy_function_to_box(row2_box, &
               ket_start_in_box(1,batch_count), &
               ket_start_in_box(2,batch_count), &
               ket_start_in_box(3,batch_count), &
               ngwf_basis%tight_boxes(row2), ngwfs_on_grid, &
               ngwf_basis%spheres(row2), cell, fftbox)

       else
          call data_set_to_zero(row2_box)
       end if

       if (row1 <= local_end) then
          ! agrecocmplx: if row1_box and row2_box are real,
          ! fourier interpolate them at the same time
          if (.not.loc_cmplx) then
              ! cks: interpolate row1 and row2 ngwfs
              call fourier_interpolate(coarse_work, fine_work, &
                   row1_box%d, row2_box%d, row1_box_dbl%d,row2_box_dbl%d)
          ! agrecocmplx: if row1_box and row2_box are complex,
          ! fourier interpolate one at a time
          else
              call fourier_interpolate(coarse_work, row1_box_dbl%z, &
                   row1_box%z)
              if (row2 <= local_end) then
                  call fourier_interpolate(coarse_work, row2_box_dbl%z, &
                       row2_box%z)
              end if
          end if
       else
          call data_set_to_zero(row1_box_dbl)
          call data_set_to_zero(row2_box_dbl)
       end if

       ! cks: Put potential into FFTbox for first NGWF
       batch_count = iii-local_start+1
       fftbox_start1_dbl = 2*ket_box_start(1,batch_count) - 1
       fftbox_start2_dbl = 2*ket_box_start(2,batch_count) - 1
       fftbox_start3_dbl = 2*ket_box_start(3,batch_count) - 1


       ! gcc32: get distance from box origin to atom of local function
       call basis_box_origin_to_atom(box_to_atom,box_origin, &
            ngwf_basis%spheres(row1)%centre, &
            fftbox_start1_dbl,fftbox_start2_dbl,fftbox_start3_dbl, &
            grid%da1,grid%da2,grid%da3,cell)

       potential_box_dbl=0.0_DP

       width = 4.0_DP * log(100000.0_DP) / &
          ((ngwf_basis%spheres(row1)%radius)**2.0_DP)

       do i1=1,n1_dbl
         do i2=1,n2_dbl
           do i3=1,n3_dbl
               r_pos = box_to_atom - (i1-1)*grid%da1-(i2-1)*grid%da2- &
                       (i3-1)*grid%da3
               magnitude=geometry_magnitude(r_pos)
               potential_box_dbl(i1,i2,i3) = pub_confined_ngwfs_barrier* &
                  exp(-width*((ngwf_basis%spheres(row1)%radius-&
                  magnitude)**2.0_DP))
           end do
         end do
       end do


       ! Apply potential to 'row1' function:
       if (row1 <= local_end) then
!!$OMP PARALLEL WORKSHARE NUM_THREADS(pub_threads_per_fftbox)
          if (loc_cmplx) then
             row1_box_dbl%z(:,:,:) = potential_box_dbl * row1_box_dbl%z
          else
             row1_box_dbl%d(:,:,:) = potential_box_dbl * row1_box_dbl%d
          end if
!!$OMP END PARALLEL WORKSHARE
       end if

       ! cks: Put potential into FFTbox for second NGWF
       batch_count = iii-local_start+2
       fftbox_start1_dbl = 2*ket_box_start(1,batch_count) - 1
       fftbox_start2_dbl = 2*ket_box_start(2,batch_count) - 1
       fftbox_start3_dbl = 2*ket_box_start(3,batch_count) - 1


       ! gcc32: get distance from box origin to atom of local function
       call basis_box_origin_to_atom(box_to_atom,box_origin, &
            ngwf_basis%spheres(row2)%centre, &
            fftbox_start1_dbl,fftbox_start2_dbl,fftbox_start3_dbl, &
            grid%da1,grid%da2,grid%da3,cell)

       potential_box_dbl=0.0_DP

       width = 4.0_DP * log(100000.0_DP) / &
          ((ngwf_basis%spheres(row2)%radius)**2.0_DP)

       do i1=1,n1_dbl
         do i2=1,n2_dbl
           do i3=1,n3_dbl
               r_pos = box_to_atom - (i1-1)*grid%da1-(i2-1)*grid%da2- &
                       (i3-1)*grid%da3
               magnitude=geometry_magnitude(r_pos)
               potential_box_dbl(i1,i2,i3) = pub_confined_ngwfs_barrier* &
                  exp(-width*((ngwf_basis%spheres(row2)%radius-&
                  magnitude)**2.0_DP))
           end do
         end do
       end do


       ! Apply potential to 'row2' function:
       if (row2 <= local_end) then
!!$OMP PARALLEL WORKSHARE NUM_THREADS(pub_threads_per_fftbox)
          if (loc_cmplx) then
              row2_box_dbl%z(:,:,:) = potential_box_dbl * row2_box_dbl%z
          else
              row2_box_dbl%d(:,:,:) = potential_box_dbl * row2_box_dbl%d
          end if
!!$OMP END PARALLEL WORKSHARE
       end if

       ! ndmh: deposit one result of interpolation in potential_fftbox_batch
       ! ndmh: and discard other result
       if ((row1 <= local_end).and.(row2 > local_end)) then

          batch_count = iii-local_start+1

          ! Filter result to standard grid
          ! agrecocmplx: if using real NGWFs, fourier
          ! filter 2 at the same time, otherwise fourier
          ! filter one at a time
          if (.not.loc_cmplx) then
              call fourier_filter(coarse_work, fine_work, &
                   row1_box_dbl%d, row2_box_dbl%d, &
                   potential_fftbox_batch(batch_count)%d, row2_box%d)
          else
              ! only row1 because row2 > local_end
              call fourier_filter(coarse_work, fine_work, &
                   row1_box_dbl%z, &
                   potential_fftbox_batch(batch_count)%z)
          end if

       ! ndmh: deposit both results of interpolation in potential_fftbox_batch
       else if ((row1 <= local_end).and.(row2 <= local_end)) then

          batch_count = iii-local_start+2

          ! Filter result to standard grid
          ! agrecocmplx: if using real NGWFs, fourier
          ! filter 2 at the same time, otherwise fourier
          ! filter one at a time
          if (.not.loc_cmplx) then
              call fourier_filter(coarse_work, fine_work, &
                   row1_box_dbl%d, row2_box_dbl%d, &
                   potential_fftbox_batch(batch_count-1)%d, &
                   potential_fftbox_batch(batch_count)%d)
          else
              ! row1 and row2, one at a time
              call fourier_filter(coarse_work, fine_work, &
                   row1_box_dbl%z, &
                   potential_fftbox_batch(batch_count-1)%z)
              call fourier_filter(coarse_work, fine_work, &
                   row2_box_dbl%z, &
                   potential_fftbox_batch(batch_count)%z)
          end if

       end if
    end do
!$OMP END DO

    deallocate(fine_work,stat=ierr)
    call utils_dealloc_check('potential_confine_ngwf_batch','fine_work',ierr)
    deallocate(coarse_work,stat=ierr)
    call utils_dealloc_check('potential_confine_ngwf_batch','coarse_work',ierr)
    call data_fftbox_dealloc(row2_box_dbl)
    call data_fftbox_dealloc(row1_box_dbl)
    deallocate(potential_box_dbl, stat=ierr)
    call utils_dealloc_check('potential_confine_ngwf_batch',&
         'potential_box_dbl',ierr)
    call data_fftbox_dealloc(row2_box)
    call data_fftbox_dealloc(row1_box)
!$OMP END PARALLEL

    call timer_clock('potential_confine_ngwf_batch', 2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving potential_confine_ngwf_batch'

  end subroutine potential_confine_ngwf_batch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine potential_confine_ppd_funcs(pot_funcs_on_grid, & ! in-output
       fbasis, grid, cell, fftbox)              ! input


    !==========================================================================!
    ! This subroutine applies the confinement potential to a set of functions  !
    ! in PPD storage.                                                          !
    !==========================================================================!
    !  Arguments:                                                              !
    !    pot_funcs_on_grid (inout) : functions in PPD format to multiply by    !
    !                                potential                                 !
    !    fbasis (input)      : The function basis type for the functions       !
    !==========================================================================!
    ! Derived from potential_apply_to_ppd_funcs by Gabriel Constantinescu in   !
    ! July 2015                                                                !
    !==========================================================================!


    use basis, only: basis_location_func_wrt_cell, &
         basis_multiply_function_by_box, basis_box_origin_to_atom
    use cell_grid, only: GRID_INFO, cell_grid_extract_box
    use comms, only: comms_barrier, pub_my_proc_id
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero, &
         data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, OPERATOR(+), OPERATOR(-), &
         OPERATOR(*), geometry_magnitude
    use simulation_cell, only: CELL_INFO
    use rundat, only: pub_confined_ngwfs_barrier
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: fbasis
    type(FUNCTIONS), intent(inout) :: pot_funcs_on_grid
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local Variables
    type(POINT) :: box_origin,r_pos
    type(POINT) :: box_to_atom
    integer :: local_func
    integer :: cell_start1,cell_start2,cell_start3
    integer :: box_n1,box_n2,box_n3
    type(FFTBOX_DATA) :: potential_box
    logical :: i_need_box
    integer :: i1,i2,i3
    real(kind=DP) :: magnitude
    real(kind=DP) :: width

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &potential_confine_ppd_funcs'

    ! Start timer
    call timer_clock('potential_confine_ppd_funcs',1)

    box_n1 = fbasis%maxtight_pts1
    box_n2 = fbasis%maxtight_pts2
    box_n3 = fbasis%maxtight_pts3

    ! Allocate storage for tightbox of potential and buffer
    call data_fftbox_alloc(potential_box, box_n1, box_n2, box_n3, .false.)

    ! Loop over ket functions up to max on any proc
    do local_func=1,fbasis%max_on_proc

       if (local_func <= fbasis%num_on_proc(pub_my_proc_id)) then
          ! Find position of ket function wrt simulation cell:
          call basis_location_func_wrt_cell(cell_start1, &
               cell_start2,cell_start3,fbasis%tight_boxes(local_func), &
               cell)
          i_need_box = .true.
       else
          cell_start1 = -1234
          cell_start2 = -1234
          cell_start3 = -1234
          i_need_box = .false.
       end if

       ! Get distance from box origin to atom of local function
       call basis_box_origin_to_atom(box_to_atom,box_origin, &
            fbasis%spheres(local_func)%centre, &
            cell_start1,cell_start2,cell_start3, &
            grid%da1,grid%da2,grid%da3,cell)

       call data_set_to_zero(potential_box)


       width = 4.0_DP * log(100000.0_DP) / &
          ((fbasis%spheres(local_func)%radius)**2.0_DP)

       do i1=1,box_n1
         do i2=1,box_n2
           do i3=1,box_n3
               r_pos = box_to_atom - (i1-1)*grid%da1-(i2-1)*grid%da2- &
                       (i3-1)*grid%da3
               magnitude=geometry_magnitude(r_pos)
               potential_box%d(i1,i2,i3) = pub_confined_ngwfs_barrier* &
                  exp(-width*((fbasis%spheres(local_func)%radius-&
                  magnitude)**2.0_DP))
           end do
         end do
       end do



       ! Multiply the ket function by the locpot box directly in PPDs
       if (local_func <= fbasis%num_on_proc(pub_my_proc_id)) then
          call basis_multiply_function_by_box(pot_funcs_on_grid, potential_box,&
               fbasis%spheres(local_func), fbasis%tight_boxes(local_func), &
               1, 1, 1, fbasis%spheres(local_func)%offset, cell, fftbox)
       end if

    end do

    ! Deallocate workspace
    call data_fftbox_dealloc(potential_box)

    ! pdh: re-sync procs
    call comms_barrier

    ! Stop timer for total
    call timer_clock('potential_confine_ppd_funcs',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &potential_confine_ppd_funcs'

  end subroutine potential_confine_ppd_funcs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine potential_sawtooth_efield(locpot_on_grid,grid,cell)

    !==========================================================================!
    ! This subroutine adds to the local potential a "sawtooth" contribution    !
    ! due to a uniform external electric field. The electric field is given    !
    ! in the input file in terms of its cartesian coordinates in atomic units. !
    ! 1.0 volt/angstrom = 0.01944655 Eh/(e*a0).                                !
    !==========================================================================!
    ! Arguments:                                                               !
    ! locpot_on_grid (input-output): local potential in simulation cell file   !
    !                             to which the electric field contribution     !
    !                             is to be added.                              !
    !==========================================================================!
    ! Written by Chris-Kriton Skylaris on 9/5/2004 for the ONETEP              !
    ! linear-scaling DFT program.                                              !
    ! Rewritten by Peter Haynes 1/7/2004 for fourier parallelisation.          !
    ! Modified 22/01/2011 by Nicholas Hine to use cell_grid_real_pt routine to !
    ! save on storage.                                                         !
    ! Modified by Andrea Greco on 26/12/2015 to allow e-field origin shift     !
    !==========================================================================!

    use cell_grid,         only: GRID_INFO, cell_grid_real_pt
    use constants,         only: DP, TWO_PI
    use rundat,            only: pub_constant_efield, pub_efield_origin
    use simulation_cell,   only: CELL_INFO
    use geometry,          only: POINT, operator(.DOT.)
    use utils,             only: utils_abort

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)  :: grid
    real(kind=DP), intent(inout) :: locpot_on_grid(grid%ld1,grid%ld2, &
         grid%max_slabs12)
    ! agreco: we need the cell vectors in order to shift
    ! the atomic positions when using a non-zero origin
    ! for the sawtooth potential associated to the e-field
    type(CELL_INFO), intent(in) :: cell

    ! Local variables
    integer :: i1,i2,islab12      ! Loop counters
    real(kind=DP) :: rpt(3)
    ! agreco: POINT types to store e-field origin and current point
    type(POINT) :: origin, pos
    ! agreco: fractional coordinates
    real(kind=DP) :: f1, f2, f3, f1_ion, f2_ion, f3_ion, rx, ry, rz
    logical :: in_a1, in_a2, in_a3, origin_shift
    ! agreco: tolerance to check the specified pub_efield_origin
    ! is within the simulation cell
    real(kind=DP), parameter :: tolerance = 1.0e-10_DP

    ! agreco: initialise origin variable
    origin%x = pub_efield_origin(1)
    origin%y = pub_efield_origin(2)
    origin%z = pub_efield_origin(3)

    f1 = 0.0_DP; f1_ion = 0.0_DP
    f2 = 0.0_DP; f2_ion = 0.0_DP
    f3 = 0.0_DP; f3_ion = 0.0_DP

    origin_shift = .false.

    ! agreco: if the origin of the sawtooth potential is non-zero,
    ! check it is within the simulation cell
    if (any(pub_efield_origin(:) /= 0.0_DP)) then
       f1 = (origin.DOT.cell%b1) / TWO_PI
       f2 = (origin.DOT.cell%b2) / TWO_PI
       f3 = (origin.DOT.cell%b3) / TWO_PI

       in_a1 = (f1 .ge. -tolerance) .and. (f1 .le. 1.0_DP + tolerance)
       in_a2 = (f2 .ge. -tolerance) .and. (f2 .le. 1.0_DP + tolerance)
       in_a3 = (f3 .ge. -tolerance) .and. (f3 .le. 1.0_DP + tolerance)

       ! agreco: abort if origin is outside simulation cell
       if (.not. in_a1 .or. .not. in_a2 .or. .not. in_a3) then
          call utils_abort('Error in potential_add_efield_ion_energy: &
          &the origin of the electric field is outside the simulation cell')
       end if

       origin_shift = .true.
    end if

    ! Loop over real-space grid on this proc and add contribution to locpot
    ! at each point
    do islab12=1,grid%num_my_slabs12
       do i2=1,grid%n2
          do i1=1,grid%n1
             call cell_grid_real_pt(rpt,i1,i2,islab12,grid)
             locpot_on_grid(i1,i2,islab12) = locpot_on_grid(i1,i2,islab12) + &
                  sum(pub_constant_efield(1:3)*rpt(1:3))

             if (origin_shift) then
                rx = 0.0_DP; ry = 0.0_DP; rz = 0.0_DP
                pos%x = rpt(1); pos%y = rpt(2); pos%z = rpt(3)

                ! agreco: fractional coordinates of current atom
                f1_ion = (pos.DOT.cell%b1) / TWO_PI
                f2_ion = (pos.DOT.cell%b2) / TWO_PI
                f3_ion = (pos.DOT.cell%b3) / TWO_PI

                ! shift origin of e-field by a lattice vector in order to get
                ! the correct value of the energy at the position of current atom
                if (f1_ion < f1) then
                   rx = rx - cell%a1%x
                   ry = ry - cell%a1%y
                   rz = rz - cell%a1%z
                end if

                if (f2_ion < f2) then
                   rx = rx - cell%a2%x
                   ry = ry - cell%a2%y
                   rz = rz - cell%a2%z
                end if

                if (f3_ion < f3) then
                   rx = rx - cell%a3%x
                   ry = ry - cell%a3%y
                   rz = rz - cell%a3%z
                end if

                ! compute extra contribution due to shifted origin of e-field
                locpot_on_grid(i1,i2,islab12) = locpot_on_grid(i1,i2,islab12) &
                   - pub_constant_efield(1)*(origin%x+rx) &
                   - pub_constant_efield(2)*(origin%y+ry) &
                   - pub_constant_efield(3)*(origin%z+rz)

             end if

          end do
       end do
    end do

  end subroutine potential_sawtooth_efield


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine potential_add_efield_ion_energy(ewald_energy, cell, par)

    !==========================================================================!
    ! This subroutine adds to the ion-ion energy the contribution from the     !
    ! interaction of a uniform external electric field with the ion charges.   !
    ! The electric field is given in the input file in terms of its cartesian  !
    ! coordinates in atomic units. 1.0 volt/angstrom = 0.01944655 Eh/(e*a0).   !
    !==========================================================================!
    ! Arguments:                                                               !
    ! ewald_energy (inout) : ion-ion energy to which to add efield-ion energy. !
    !==========================================================================!
    ! Written by Nicholas Hine on 10/06/2011.                                  !
    ! Modified by Andrea Greco on 26/12/2015 to allow shift of e-field origin  !
    !==========================================================================!

    use comms,             only: comms_reduce, pub_my_proc_id
    use constants,         only: DP, TWO_PI
    use parallel_strategy, only: PARAL_INFO!par=>pub_par
    use rundat,            only: pub_constant_efield, pub_efield_origin
    use simulation_cell,   only: CELL_INFO
    use geometry,          only: POINT, operator(.DOT.)
    use utils,             only: utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: ewald_energy
    ! agreco: we need the cell vectors in order to shift
    ! the atomic positions when using a non-zero origin
    ! for the sawtooth potential associated to the e-field
    type(CELL_INFO), intent(in) :: cell
    type(PARAL_INFO), intent(in), pointer :: par ! rc2013: parallel strategy

    ! Local variables
    integer :: iat
    real(kind=DP) :: efield_energy
    ! agreco: POINT type to store e-field origin
    type(POINT) :: origin
    ! agreco: fractional coordinates
    real(kind=DP) :: f1, f2, f3, f1_ion, f2_ion, f3_ion, rx, ry, rz
    logical :: in_a1, in_a2, in_a3, origin_shift
    ! agreco: tolerance to check the specified pub_efield_origin
    ! is within the simulation cell
    real(kind=DP), parameter :: tolerance = 1.0e-10_DP

    ! agreco: initialise origin variable
    origin%x = pub_efield_origin(1)
    origin%y = pub_efield_origin(2)
    origin%z = pub_efield_origin(3)

    f1 = 0.0_DP; f1_ion = 0.0_DP
    f2 = 0.0_DP; f2_ion = 0.0_DP
    f3 = 0.0_DP; f3_ion = 0.0_DP

    origin_shift = .false.

    ! agreco: if the origin of the sawtooth potential is non-zero,
    ! check it is within the simulation cell
    if (any(pub_efield_origin(:) /= 0.0_DP)) then
       f1 = (origin.DOT.cell%b1) / TWO_PI
       f2 = (origin.DOT.cell%b2) / TWO_PI
       f3 = (origin.DOT.cell%b3) / TWO_PI

       in_a1 = (f1 .ge. -tolerance) .and. (f1 .le. 1.0_DP + tolerance)
       in_a2 = (f2 .ge. -tolerance) .and. (f2 .le. 1.0_DP + tolerance)
       in_a3 = (f3 .ge. -tolerance) .and. (f3 .le. 1.0_DP + tolerance)

       ! agreco: abort if origin is outside simulation cell
       if (.not. in_a1 .or. .not. in_a2 .or. .not. in_a3) then
          call utils_abort('Error in potential_add_efield_ion_energy: &
          &the origin of the electric field is outside the simulation cell')
       end if

       origin_shift = .true.
    end if

    ! Loop over atoms on this proc and calculate energy of interaction
    ! for each
    efield_energy = 0.0_DP
    do iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       efield_energy = efield_energy - par%elements_on_proc(iat)%ion_charge &
            *(pub_constant_efield(1)*par%elements_on_proc(iat)%centre%x &
            + pub_constant_efield(2)*par%elements_on_proc(iat)%centre%y &
            + pub_constant_efield(3)*par%elements_on_proc(iat)%centre%z)

       if (origin_shift) then
          rx = 0.0_DP; ry = 0.0_DP; rz = 0.0_DP

          ! agreco: fractional coordinates of current atom
          f1_ion = (par%elements_on_proc(iat)%centre.DOT.cell%b1) / TWO_PI
          f2_ion = (par%elements_on_proc(iat)%centre.DOT.cell%b2) / TWO_PI
          f3_ion = (par%elements_on_proc(iat)%centre.DOT.cell%b3) / TWO_PI

          ! shift origin of e-field by a lattice vector in order to get
          ! the correct value of the energy at the position of current atom
          if (f1_ion < f1) then
             rx = rx - cell%a1%x
             ry = ry - cell%a1%y
             rz = rz - cell%a1%z
          end if

          if (f2_ion < f2) then
             rx = rx - cell%a2%x
             ry = ry - cell%a2%y
             rz = rz - cell%a2%z
          end if

          if (f3_ion < f3) then
             rx = rx - cell%a3%x
             ry = ry - cell%a3%y
             rz = rz - cell%a3%z
          end if

          ! compute extra contribution due to shifted origin of e-field
          efield_energy = efield_energy + par%elements_on_proc(iat)%ion_charge &
             *(pub_constant_efield(1)*(origin%x+rx) &
             + pub_constant_efield(2)*(origin%y+ry) &
             + pub_constant_efield(3)*(origin%z+rz))

       end if

    end do

    call comms_reduce('SUM',efield_energy)

    ewald_energy = ewald_energy + efield_energy

  end subroutine potential_add_efield_ion_energy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine potential_apply_to_ngwf_batch(potential_fftbox_batch, &   ! output
       ngwfs_on_grid, ngwf_basis, potential_dbl, grid, fftbox, cell, & ! input
       ket_start_in_box, ket_box_start, batch_size, local_start, &     ! input
       local_end, max_current_size)                                    ! input

    !==========================================================================!
    ! This subroutine applies the local potential operator to a batch          !
    ! of NGWFs in fftboxes.                                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! potential_fftbox_batch (output): Batch off fftboxes with potential-NGWFs !
    ! ngwfs_on_grid (input) : NGWFs on this proc in ppd representation         !
    ! ngwf_basis (input)    : Function basis type for NGWFs                    !
    ! potential_dbl (input) : Local potential on dbl grid in simulation cell   !
    ! batch_size (input)    : Number of fftboxes in the batch                  !
    ! local_start(input)    : First NGWF of the batch                          !
    ! local_end  (input)    : Last NGWF of the batch                           !
    ! ket_start_in_box (input) :
    !--------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 23/11/2003 for the ONETEP program.   !
    ! Modified by Nicholas Hine in 2008 to replace two-stage copy of each row  !
    ! function with single call to basis_copy_function_to_box.                 !
    ! Modified by Nicholas Hine on 21/07/2009 to use function basis type and   !
    ! to deposit result of fourier_filter straight into potential_fftbox array.!
    ! OpenMP parallelised by Karl Wilkinson and Nicholas Hine in May 2013.     !
    ! Modified by Andrea Greco in May 2015 to allow use of complex NGWFs.      !
    !==========================================================================!

    use basis, only: basis_copy_function_to_box
    use cell_grid, only: GRID_INFO, cell_grid_extract_box
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc, data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_filter, fourier_interpolate
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_threads_fftbox
!$  use rundat, only: pub_threads_num_fftboxes, pub_threads_per_fftbox
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_use_var
    use fourier, only: fourier_filter, fourier_interpolate

    implicit none

    ! Arguments
    integer, intent(in) :: batch_size, local_start, local_end
    integer, intent(in) :: max_current_size ! max size of batch over all procs
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FFTBOX_DATA), intent(inout) :: potential_fftbox_batch(batch_size)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    type(GRID_INFO), intent(in)  :: grid
    real(kind=DP), intent(in) :: potential_dbl(grid%ld1, &
         grid%ld2, grid%max_group_slabs12)
    integer, intent(in) :: ket_start_in_box(3,batch_size)
    integer, intent(in) :: ket_box_start(3,batch_size)

    ! Local Variables
    type(FFTBOX_DATA) :: row1_box, row2_box
    real(kind=DP), dimension(:,:,:), allocatable :: potential_box_dbl
    type(FFTBOX_DATA) :: row1_box_dbl, row2_box_dbl
    real(kind=DP), dimension(:,:,:), allocatable :: buffer_dbl

    complex(kind=DP), allocatable, dimension(:,:,:) :: coarse_work,fine_work

    integer :: ierr
    integer :: row1, row2
    integer :: n1, n2, n3, ld1, ld2
    integer :: n1_dbl, n2_dbl, n3_dbl, ld1_dbl, ld2_dbl
    integer :: iii
    integer :: batch_count
    integer :: fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl
    integer :: prev_start1, prev_start2, prev_start3

    logical :: i_need_potential ! whether pub_my_proc_id needs pot'l in fftbox
    ! agrecocmplx
    logical :: loc_cmplx

    ! -------------------------------------------------------------------------

    call utils_use_var(pub_threads_fftbox)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering potential_apply_to_ngwf_batch'

    ! Start timer
    call timer_clock('potential_apply_to_ngwf_batch', 1)

    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2
    n1_dbl = fftbox%total_pt1_dbl
    n2_dbl = fftbox%total_pt2_dbl
    n3_dbl = fftbox%total_pt3_dbl
    ld1_dbl = fftbox%total_ld1_dbl
    ld2_dbl = fftbox%total_ld2_dbl

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid%iscmplx

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(potential_box_dbl, buffer_dbl,iii, batch_count, row1, row2, &
!$OMP      ierr, i_need_potential, fftbox_start1_dbl, fftbox_start2_dbl, &
!$OMP      fftbox_start3_dbl, prev_start1, prev_start2, prev_start3, &
!$OMP      coarse_work, fine_work) &
!$OMP FIRSTPRIVATE(row1_box,row2_box,row1_box_dbl,row2_box_dbl) &
!$OMP SHARED(potential_fftbox_batch, ket_start_in_box, ket_box_start, &
!$OMP      potential_dbl, ngwfs_on_grid, n1, n2, n3, ld1, ld2, n1_dbl, &
!$OMP      n2_dbl, n3_dbl, ld1_dbl, ld2_dbl, local_start, local_end, &
!$OMP      cell, fftbox, grid, ngwf_basis, max_current_size, loc_cmplx, &
!$OMP      pub_threads_per_fftbox,pub_threads_num_fftboxes)

    ! ndmh: initialisations to prevent compiler warnings
    prev_start1 = -1111111
    prev_start2 = -2222222
    prev_start3 = -3333333
    fftbox_start1_dbl = -1111111
    fftbox_start2_dbl = -2222222
    fftbox_start3_dbl = -3333333

    ! Allocate workspace arrays
    call data_fftbox_alloc(row1_box, ld1, ld2, n3, &
         iscmplx=loc_cmplx)
    call data_fftbox_alloc(row2_box, ld1, ld2, n3, &
         iscmplx=loc_cmplx)
    allocate(potential_box_dbl(ld1_dbl, ld2_dbl, 2*n3), stat=ierr)
    call utils_alloc_check('potential_apply_to_ngwf_batch',&
         'potential_box_dbl',ierr)
    call data_fftbox_alloc(row1_box_dbl, ld1_dbl, ld2_dbl, 2*n3, &
         iscmplx=loc_cmplx)
    call data_fftbox_alloc(row2_box_dbl, ld1_dbl, ld2_dbl, 2*n3, &
         iscmplx=loc_cmplx)
    allocate(buffer_dbl(ld1_dbl, ld2_dbl, grid%max_group_slabs12), &
         stat=ierr)
    call utils_alloc_check('potential_apply_to_ngwf_batch','buffer_dbl',&
         ierr)
    allocate(coarse_work(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('potential_apply_to_ngwf_batch','coarse_work',ierr)
    allocate(fine_work(ld1*2, ld2*2, n3*2),stat=ierr)
    call utils_alloc_check('potential_apply_to_ngwf_batch','fine_work',ierr)

    ! cks: initialisation
    buffer_dbl = 0.0_DP
    ! agrecocmplx
    call data_set_to_zero(row1_box_dbl)
    potential_box_dbl = 0.0_DP

!$OMP DO
    do iii=local_start, local_start+max_current_size-1, 2

       row1 = iii ; row2 = iii + 1

       if (row1 <= local_end) then

          batch_count = iii-local_start+1

          ! ndmh: use new basis function copying straight to fftbox from ppds
          ! agrecocmplx: call modified version in basis_new
          call basis_copy_function_to_box(row1_box, &
               ket_start_in_box(1,batch_count), &
               ket_start_in_box(2,batch_count), &
               ket_start_in_box(3,batch_count), &
               ngwf_basis%tight_boxes(row1), ngwfs_on_grid, &
               ngwf_basis%spheres(row1), cell, fftbox)

       else

          !row1_box = 0.0_DP
          ! agrecocmplx
          call data_set_to_zero(row1_box)

       end if

       if (row2 <= local_end) then

          batch_count = iii-local_start+2

          ! ndmh: use new basis function copying straight to fftbox from ppds
          ! agrecocmplx: call modified version in basis_new
          call basis_copy_function_to_box(row2_box, &
               ket_start_in_box(1,batch_count), &
               ket_start_in_box(2,batch_count), &
               ket_start_in_box(3,batch_count), &
               ngwf_basis%tight_boxes(row2), ngwfs_on_grid, &
               ngwf_basis%spheres(row2), cell, fftbox)

       else

          !row2_box = 0.0_DP
          ! agrecocmplx
          call data_set_to_zero(row2_box)

       end if

       if (row1 <= local_end) then
          ! agrecocmplx: if row1_box and row2_box are real,
          ! fourier interpolate them at the same time
          if (.not.loc_cmplx) then
              ! cks: interpolate row1 and row2 ngwfs
              call fourier_interpolate(coarse_work, fine_work, &
                   row1_box%d, row2_box%d, row1_box_dbl%d,row2_box_dbl%d)
          ! agrecocmplx: if row1_box and row2_box are complex,
          ! fourier interpolate one at a time
          else
              call fourier_interpolate(coarse_work, row1_box_dbl%z, &
                   row1_box%z)
              if (row2 <= local_end) then
                  call fourier_interpolate(coarse_work, row2_box_dbl%z, &
                       row2_box%z)
              end if
          end if
       else
          ! agrecocmplx
          call data_set_to_zero(row1_box_dbl)
          call data_set_to_zero(row2_box_dbl)
       end if

       ! cks: Put potential into FFTbox for first NGWF
       batch_count = iii-local_start+1
       fftbox_start1_dbl = 2*ket_box_start(1,batch_count) - 1
       fftbox_start2_dbl = 2*ket_box_start(2,batch_count) - 1
       fftbox_start3_dbl = 2*ket_box_start(3,batch_count) - 1

       ! cks:----- Get potential for FIRST NGWF if needed ------------
       i_need_potential = ((fftbox_start1_dbl /= prev_start1 .or. &
            fftbox_start2_dbl /= prev_start2 .or. &
            fftbox_start3_dbl /= prev_start3) .and. &
            row1 <= local_end)

!$OMP CRITICAL
       call cell_grid_extract_box(potential_box_dbl,&
            buffer_dbl, potential_dbl, grid, n1_dbl, n2_dbl, n3_dbl, &
            ld1_dbl, ld2_dbl, fftbox_start1_dbl, fftbox_start2_dbl, &
            fftbox_start3_dbl, i_need_potential, .true.)
!$OMP END CRITICAL

       if (i_need_potential) then
          prev_start1 = fftbox_start1_dbl
          prev_start2 = fftbox_start2_dbl
          prev_start3 = fftbox_start3_dbl
       end if
       ! cks:-- END Get potential for FIRST NGWF if needed ------------

       ! Apply potential to 'row1' function:
       if (row1 <= local_end) then
!!$OMP PARALLEL WORKSHARE NUM_THREADS(pub_threads_per_fftbox)
          ! agrecocmplx
          if (loc_cmplx) then
              row1_box_dbl%z(:,:,:) = potential_box_dbl * row1_box_dbl%z
          else
              row1_box_dbl%d(:,:,:) = potential_box_dbl * row1_box_dbl%d
          end if
!!$OMP END PARALLEL WORKSHARE
       end if

       ! cks: Put potential into FFTbox for second NGWF
       batch_count = iii-local_start+2
       if (row2 <= local_end) then
          fftbox_start1_dbl = 2*ket_box_start(1,batch_count) - 1
          fftbox_start2_dbl = 2*ket_box_start(2,batch_count) - 1
          fftbox_start3_dbl = 2*ket_box_start(3,batch_count) - 1
       else
          fftbox_start1_dbl = -999
          fftbox_start2_dbl = -999
          fftbox_start3_dbl = -999
       end if

       ! cks:----- Get potential for SECOND NGWF if needed -----------
       i_need_potential = ((fftbox_start1_dbl /= prev_start1 .or. &
            fftbox_start2_dbl /= prev_start2 .or. &
            fftbox_start3_dbl /= prev_start3) .and. &
            row2 <= local_end)

!$OMP CRITICAL
       call cell_grid_extract_box(potential_box_dbl, &
            buffer_dbl, potential_dbl, grid, n1_dbl, n2_dbl, n3_dbl, &
            ld1_dbl, ld2_dbl, fftbox_start1_dbl, fftbox_start2_dbl, &
            fftbox_start3_dbl, i_need_potential, .true.)
!$OMP END CRITICAL

       if (i_need_potential) then
          prev_start1 = fftbox_start1_dbl
          prev_start2 = fftbox_start2_dbl
          prev_start3 = fftbox_start3_dbl
       endif
       ! cks:-- END Get potential for SECOND NGWF if needed -----------

       ! Apply potential to 'row2' function:
       if (row2 <= local_end) then
!!$OMP PARALLEL WORKSHARE NUM_THREADS(pub_threads_per_fftbox)
          ! agrecocmplx
          if (loc_cmplx) then
              row2_box_dbl%z(:,:,:) = potential_box_dbl * row2_box_dbl%z
          else
              row2_box_dbl%d(:,:,:) = potential_box_dbl * row2_box_dbl%d
          end if
!!$OMP END PARALLEL WORKSHARE
       end if

       ! ndmh: deposit one result of interpolation in potential_fftbox_batch
       ! ndmh: and discard other result
       if ((row1 <= local_end).and.(row2 > local_end)) then

          batch_count = iii-local_start+1

          ! Filter result to standard grid
          ! agrecocmplx: if using real NGWFs, fourier
          ! filter 2 at the same time, otherwise fourier
          ! filter one at a time
          if (.not.loc_cmplx) then
              call fourier_filter(coarse_work, fine_work, &
                   row1_box_dbl%d, row2_box_dbl%d, &
                   potential_fftbox_batch(batch_count)%d, row2_box%d)
          else
              ! only row1 because row2 > local_end
              call fourier_filter(coarse_work, fine_work, &
                   row1_box_dbl%z, &
                   potential_fftbox_batch(batch_count)%z)
          end if

       ! ndmh: deposit both results of interpolation in potential_fftbox_batch
       else if ((row1 <= local_end).and.(row2 <= local_end)) then

          batch_count = iii-local_start+2

          ! Filter result to standard grid
          ! agrecocmplx: if using real NGWFs, fourier
          ! filter 2 at the same time, otherwise fourier
          ! filter one at a time
          if (.not.loc_cmplx) then
              call fourier_filter(coarse_work, fine_work, &
                   row1_box_dbl%d, row2_box_dbl%d, &
                   potential_fftbox_batch(batch_count-1)%d, &
                   potential_fftbox_batch(batch_count)%d)
          else
              ! row1 and row2, one at a time
              call fourier_filter(coarse_work, fine_work, &
                   row1_box_dbl%z, &
                   potential_fftbox_batch(batch_count-1)%z)
              call fourier_filter(coarse_work, fine_work, &
                   row2_box_dbl%z, &
                   potential_fftbox_batch(batch_count)%z)
          end if
       end if
    end do
!$OMP END DO

    deallocate(fine_work,stat=ierr)
    call utils_dealloc_check('potential_apply_to_ngwf_batch','fine_work',ierr)
    deallocate(coarse_work,stat=ierr)
    call utils_dealloc_check('potential_apply_to_ngwf_batch','coarse_work',ierr)
    deallocate(buffer_dbl, stat=ierr)
    call utils_dealloc_check('potential_apply_to_ngwf_batch','buffer_dbl',&
         ierr)
    call data_fftbox_dealloc(row2_box_dbl)
    call data_fftbox_dealloc(row1_box_dbl)
    deallocate(potential_box_dbl, stat=ierr)
    call utils_dealloc_check('potential_apply_to_ngwf_batch',&
         'potential_box_dbl',ierr)
    call data_fftbox_dealloc(row2_box)
    call data_fftbox_dealloc(row1_box)
!$OMP END PARALLEL

    call timer_clock('potential_apply_to_ngwf_batch', 2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving potential_apply_to_ngwf_batch'

  end subroutine potential_apply_to_ngwf_batch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine potential_apply_to_fftbox_batch(potential_fftbox_batch, & ! output
       ket_fftbox_batch, iscmplx, potential_dbl, grid, fftbox, & ! input
       ket_start_in_box, ket_box_start, batch_size, local_start, &     ! input
       local_end, max_current_size)                                    ! input

    !==========================================================================!
    ! This subroutine applies the local potential operator to a batch of       !
    ! functions deposited in fftboxes.                                         !
    !--------------------------------------------------------------------------!
    ! TODO Add details for missing arguments                                   !
    ! Arguments:                                                               !
    ! potential_fftbox_batch (output): Batch off fftboxes with potential-NGWFs !
    ! potential_dbl (input) : Local potential on dbl grid in simulation cell   !
    ! batch_size (input)    : Number of fftboxes in the batch                  !
    ! local_start(input)    : First NGWF of the batch                          !
    ! local_end  (input)    : Last NGWF of the batch                           !
    ! ket_start_in_box (input) :                                               !
    !--------------------------------------------------------------------------!
    ! Modified version of potential_apply_to_func_batch so that fftboxes can   !
    ! be directly passed to the routine. J. C. Womack, 2015.                   !
    !--------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 23/11/2003 for the ONETEP program.   !
    ! Modified by Nicholas Hine in 2008 to replace two-stage copy of each row  !
    ! function with single call to basis_copy_function_to_box.                 !
    ! Modified by Nicholas Hine on 21/07/2009 to use function basis type and   !
    ! to deposit result of fourier_filter straight into potential_fftbox array.!
    ! OpenMP parallelised by Karl Wilkinson and Nicholas Hine in May 2013.     !
    ! Modified by Andrea Greco in May 2015 to allow use of complex NGWFs.      !
    !==========================================================================!
    ! TODO JCW: This routine could do with removing some old commented-out code
    ! agrecocmplx
    use basis, only: basis_copy_function_to_box
    use cell_grid, only: GRID_INFO, cell_grid_extract_box
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, &
         data_fftbox_alloc, data_fftbox_dealloc, &
         data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_filter, fourier_interpolate
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_threads_fftbox
!$  use rundat, only: pub_threads_num_fftboxes, pub_threads_per_fftbox
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_use_var, &
         utils_assert
    use fourier, only: fourier_filter, fourier_interpolate

    implicit none

    ! Arguments
    integer, intent(in) :: batch_size, local_start, local_end
    integer, intent(in) :: max_current_size ! max size of batch over all procs
    type(FFTBOX_INFO), intent(in) :: fftbox
    !real(kind=DP), intent(out) :: potential_fftbox_batch(fftbox%total_ld1, &
    !     fftbox%total_ld2, fftbox%total_pt3, batch_size)
    ! agrecocmplx: use new FFTBOX_DATA type
    type(FFTBOX_DATA), intent(inout) :: potential_fftbox_batch(batch_size)
    type(FFTBOX_DATA), intent(in) :: ket_fftbox_batch(batch_size)
    logical, intent(in)          :: iscmplx
    type(GRID_INFO), intent(in)  :: grid
    real(kind=DP), intent(in) :: potential_dbl(grid%ld1, &
         grid%ld2, grid%max_group_slabs12)
    integer, intent(in) :: ket_start_in_box(3,batch_size)
    integer, intent(in) :: ket_box_start(3,batch_size)

    ! Local Variables
    !real(kind=DP), dimension(:,:,:), allocatable :: row1_box, row2_box
    ! agrecocmplx: use new FFTBOX_DATA types
    !type(FFTBOX_DATA) :: row1_box ! JCW: No longer needed
    type(FFTBOX_DATA) :: row2_box
    type(FFTBOX_DATA) :: zero_box
    real(kind=DP), dimension(:,:,:), allocatable :: potential_box_dbl
    !real(kind=DP), dimension(:,:,:), allocatable :: row1_box_dbl
    !real(kind=DP), dimension(:,:,:), allocatable :: row2_box_dbl
    ! agrecocmplx: use new FFTBOX_DATA type
    type(FFTBOX_DATA) :: row1_box_dbl, row2_box_dbl
    real(kind=DP), dimension(:,:,:), allocatable :: buffer_dbl

    complex(kind=DP), allocatable, dimension(:,:,:) :: coarse_work,fine_work

    integer :: ierr
    integer :: row1, row2
    integer :: n1, n2, n3, ld1, ld2
    integer :: n1_dbl, n2_dbl, n3_dbl, ld1_dbl, ld2_dbl
    integer :: iii
    integer :: batch_count1, batch_count2
    integer :: fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl
    integer :: prev_start1, prev_start2, prev_start3

    logical :: i_need_potential ! whether pub_my_proc_id needs pot'l in fftbox

    ! -------------------------------------------------------------------------

    call utils_use_var(pub_threads_fftbox)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering potential_apply_to_fftbox_batch'

    ! Start timer
    call timer_clock('potential_apply_to_fftbox_batch', 1)

    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2
    n1_dbl = fftbox%total_pt1_dbl
    n2_dbl = fftbox%total_pt2_dbl
    n3_dbl = fftbox%total_pt3_dbl
    ld1_dbl = fftbox%total_ld1_dbl
    ld2_dbl = fftbox%total_ld2_dbl

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(potential_box_dbl,buffer_dbl, iii, batch_count1, batch_count2, &
!$OMP      row1, row2, ierr, &
!$OMP      i_need_potential, fftbox_start1_dbl, fftbox_start2_dbl, &
!$OMP      fftbox_start3_dbl, prev_start1, prev_start2, prev_start3, &
!$OMP      coarse_work, fine_work) &
!$OMP FIRSTPRIVATE(row2_box,row1_box_dbl,row2_box_dbl,zero_box) &
!$OMP SHARED(potential_fftbox_batch, ket_start_in_box, ket_box_start, &
!$OMP      ket_fftbox_batch, potential_dbl, n1, n2, n3, ld1, ld2, n1_dbl, &
!$OMP      n2_dbl, n3_dbl, ld1_dbl, ld2_dbl, local_start, local_end, &
!$OMP      iscmplx, fftbox, grid, max_current_size, &
!$OMP      pub_threads_per_fftbox,pub_threads_num_fftboxes)


    ! ndmh: initialisations to prevent compiler warnings
    prev_start1 = -1111111
    prev_start2 = -2222222
    prev_start3 = -3333333
    fftbox_start1_dbl = -1111111
    fftbox_start2_dbl = -2222222
    fftbox_start3_dbl = -3333333

    ! Allocate workspace arrays
    !allocate(row1_box(ld1, ld2, n3), stat=ierr)
    !call utils_alloc_check('potential_apply_to_fftbox_batch','row1_box',ierr)
    !allocate(row2_box(ld1, ld2, n3), stat=ierr)
    !call utils_alloc_check('potential_apply_to_fftbox_batch','row2_box',ierr)
    ! agrecocmplx: allocate using appropriate routine
    !call data_fftbox_alloc(row1_box, ld1, ld2, n3, &
    !     iscmplx=iscmplx) ! JCW: Coarse grid row1_box is no longer needed
    if (.not.iscmplx) then
       ! JCW: Coarse grid row2_box only needed for real NGWFs
       call data_fftbox_alloc(row2_box, ld1, ld2, n3, &
            iscmplx=iscmplx)
    end if
    call data_fftbox_alloc(zero_box, ld1, ld2, n3, &
         iscmplx=iscmplx)
    allocate(potential_box_dbl(ld1_dbl, ld2_dbl, 2*n3), stat=ierr)
    call utils_alloc_check('potential_apply_to_fftbox_batch',&
         'potential_box_dbl',ierr)
    !allocate(row1_box_dbl(ld1_dbl, ld2_dbl, 2*n3), stat=ierr)
    !call utils_alloc_check('potential_apply_to_fftbox_batch','row1_box_dbl',&
    !     ierr)
    !allocate(row2_box_dbl(ld1_dbl, ld2_dbl, 2*n3), stat=ierr)
    !call utils_alloc_check('potential_apply_to_fftbox_batch','row2_box_dbl',&
    !     ierr)
    ! agrecocmplx: allocate using appropriate routine
    call data_fftbox_alloc(row1_box_dbl, ld1_dbl, ld2_dbl, 2*n3, &
         iscmplx=iscmplx)
    call data_fftbox_alloc(row2_box_dbl, ld1_dbl, ld2_dbl, 2*n3, &
         iscmplx=iscmplx)
    allocate(buffer_dbl(ld1_dbl, ld2_dbl, grid%max_group_slabs12), &
         stat=ierr)
    call utils_alloc_check('potential_apply_to_fftbox_batch','buffer_dbl',&
         ierr)
    allocate(coarse_work(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('potential_apply_to_fftbox_batch','coarse_work',ierr)
    allocate(fine_work(ld1*2, ld2*2, n3*2),stat=ierr)
    call utils_alloc_check('potential_apply_to_fftbox_batch','fine_work',ierr)

    ! cks: initialisation
    buffer_dbl = 0.0_DP
    ! agrecocmplx
    ! JCW: TODO zero_box is simply a fftbox filled with zeroes, which allows
    ! JCW: fourier_interpolate to be called with an all-zero real argument
    ! JCW: when the row2 fftbox is not present in the batch.
    ! JCW: Is there a better (less wasteful) way to do this?
    call data_set_to_zero(zero_box)
    call data_set_to_zero(row1_box_dbl)
    potential_box_dbl = 0.0_DP

!$OMP DO
    do iii=local_start, local_start+max_current_size-1, 2

       row1 = iii ; row2 = iii + 1

       batch_count1 = iii-local_start+1
       batch_count2 = iii-local_start+2

       if (row2 <= local_end) then
          ! JCW row1 must be <= local_end, since row1 = row2 - 1
          call utils_assert(.not.iscmplx,&
                   "Error in potential_apply_to_fftbox_batch: &
                   &complex functions not supported currently.")
          ! JCW: If two real functions are available, fourier interpolate
          ! JCW: them at the same time
          ! cks: interpolate row1 and row2 ngwfs
          call fourier_interpolate(coarse_work, fine_work, &
               ket_fftbox_batch(batch_count1)%d, &
               ket_fftbox_batch(batch_count2)%d, &
               row1_box_dbl%d,row2_box_dbl%d)
       else if ( (row2 > local_end).and.(row1 <= local_end) ) then
          ! JCW: fourier interpolate only row1, since no function in row2
          call fourier_interpolate(coarse_work, fine_work, &
               ket_fftbox_batch(batch_count1)%d, &
               zero_box%d, &
               row1_box_dbl%d,row2_box_dbl%d)
       else
         ! JCW: row1 and row2 > local_end: set both row{1,2}_box_dbl = 0.0
          call data_set_to_zero(row1_box_dbl)
          call data_set_to_zero(row2_box_dbl)
       end if

       ! cks: Put potential into FFTbox for first NGWF
       !batch_count1 = iii-local_start+1
       fftbox_start1_dbl = 2*ket_box_start(1,batch_count1) - 1
       fftbox_start2_dbl = 2*ket_box_start(2,batch_count1) - 1
       fftbox_start3_dbl = 2*ket_box_start(3,batch_count1) - 1

       ! cks:----- Get potential for FIRST NGWF if needed ------------
       i_need_potential = ((fftbox_start1_dbl /= prev_start1 .or. &
            fftbox_start2_dbl /= prev_start2 .or. &
            fftbox_start3_dbl /= prev_start3) .and. &
            row1 <= local_end)

!$OMP CRITICAL
       call cell_grid_extract_box(potential_box_dbl,&
            buffer_dbl, potential_dbl, grid, n1_dbl, n2_dbl, n3_dbl, &
            ld1_dbl, ld2_dbl, fftbox_start1_dbl, fftbox_start2_dbl, &
            fftbox_start3_dbl, i_need_potential, .true.)
!$OMP END CRITICAL

       if (i_need_potential) then
          prev_start1 = fftbox_start1_dbl
          prev_start2 = fftbox_start2_dbl
          prev_start3 = fftbox_start3_dbl
       end if
       ! cks:-- END Get potential for FIRST NGWF if needed ------------

       ! Apply potential to 'row1' function:
       if (row1 <= local_end) then
!!$OMP PARALLEL WORKSHARE NUM_THREADS(pub_threads_per_fftbox)
          ! agrecocmplx
          if (iscmplx) then
              row1_box_dbl%z(:,:,:) = potential_box_dbl * row1_box_dbl%z
          else
              row1_box_dbl%d(:,:,:) = potential_box_dbl * row1_box_dbl%d
          end if
!!$OMP END PARALLEL WORKSHARE
       end if

       ! cks: Put potential into FFTbox for second NGWF
       !batch_count2 = iii-local_start+2
       fftbox_start1_dbl = 2*ket_box_start(1,batch_count2) - 1
       fftbox_start2_dbl = 2*ket_box_start(2,batch_count2) - 1
       fftbox_start3_dbl = 2*ket_box_start(3,batch_count2) - 1

       ! cks:----- Get potential for SECOND NGWF if needed -----------
       i_need_potential = ((fftbox_start1_dbl /= prev_start1 .or. &
            fftbox_start2_dbl /= prev_start2 .or. &
            fftbox_start3_dbl /= prev_start3) .and. &
            row2 <= local_end)

!$OMP CRITICAL
       call cell_grid_extract_box(potential_box_dbl, &
            buffer_dbl, potential_dbl, grid, n1_dbl, n2_dbl, n3_dbl, &
            ld1_dbl, ld2_dbl, fftbox_start1_dbl, fftbox_start2_dbl, &
            fftbox_start3_dbl, i_need_potential, .true.)
!$OMP END CRITICAL

       if (i_need_potential) then
          prev_start1 = fftbox_start1_dbl
          prev_start2 = fftbox_start2_dbl
          prev_start3 = fftbox_start3_dbl
       endif
       ! cks:-- END Get potential for SECOND NGWF if needed -----------

       ! Apply potential to 'row2' function:
       if (row2 <= local_end) then
!!$OMP PARALLEL WORKSHARE NUM_THREADS(pub_threads_per_fftbox)
          ! agrecocmplx
          if (iscmplx) then
              row2_box_dbl%z(:,:,:) = potential_box_dbl * row2_box_dbl%z
          else
              row2_box_dbl%d(:,:,:) = potential_box_dbl * row2_box_dbl%d
          end if
!!$OMP END PARALLEL WORKSHARE
       end if

       ! ndmh: deposit one result of interpolation in potential_fftbox_batch
       ! ndmh: and discard other result
       if ((row1 <= local_end).and.(row2 > local_end)) then

          !batch_count1 = iii-local_start+1

          ! Filter result to standard grid
          ! agrecocmplx: if using real NGWFs, fourier
          ! filter 2 at the same time, otherwise fourier
          ! filter one at a time
          if (.not.iscmplx) then
              call fourier_filter(coarse_work, fine_work, &
                   row1_box_dbl%d, row2_box_dbl%d, &
                   potential_fftbox_batch(batch_count1)%d, row2_box%d)
              ! JCW: row2_box is thrown away in this case
          else
              ! only row1 because row2 > local_end
              call fourier_filter(coarse_work, fine_work, &
                   row1_box_dbl%z, &
                   potential_fftbox_batch(batch_count1)%z)
          end if

       ! ndmh: deposit both results of interpolation in potential_fftbox_batch
       else if ((row1 <= local_end).and.(row2 <= local_end)) then

          !batch_count2 = iii-local_start+2

          ! Filter result to standard grid
          ! agrecocmplx: if using real NGWFs, fourier
          ! filter 2 at the same time, otherwise fourier
          ! filter one at a time
          if (.not.iscmplx) then
              call fourier_filter(coarse_work, fine_work, &
                   row1_box_dbl%d, row2_box_dbl%d, &
                   potential_fftbox_batch(batch_count2-1)%d, &
                   potential_fftbox_batch(batch_count2)%d)
          else
              ! row1 and row2, one at a time
              call fourier_filter(coarse_work, fine_work, &
                   row1_box_dbl%z, &
                   potential_fftbox_batch(batch_count2-1)%z)
              call fourier_filter(coarse_work, fine_work, &
                   row2_box_dbl%z, &
                   potential_fftbox_batch(batch_count2)%z)
          end if
       end if
    end do
!$OMP END DO

    deallocate(fine_work,stat=ierr)
    call utils_dealloc_check('potential_apply_to_fftbox_batch','fine_work',ierr)
    deallocate(coarse_work,stat=ierr)
    call utils_dealloc_check('potential_apply_to_fftbox_batch','coarse_work',ierr)
    deallocate(buffer_dbl, stat=ierr)
    call utils_dealloc_check('potential_apply_to_fftbox_batch','buffer_dbl',&
         ierr)
    call data_fftbox_dealloc(zero_box)
    !deallocate(row2_box_dbl, stat=ierr)
    !call utils_dealloc_check('potential_apply_to_fftbox_batch','row2_box_dbl',&
    !    ierr)
    ! agrecocmplx: deallocate using appropriate routine
    call data_fftbox_dealloc(row2_box_dbl)
    !deallocate(row1_box_dbl, stat=ierr)
    !call utils_dealloc_check('potential_apply_to_fftbox_batch','row1_box_dbl',&
    !     ierr)
    ! agrecocmplx: deallocate using appropriate routine
    call data_fftbox_dealloc(row1_box_dbl)
    deallocate(potential_box_dbl, stat=ierr)
    call utils_dealloc_check('potential_apply_to_fftbox_batch',&
         'potential_box_dbl',ierr)
    ! agrecocmplx: deallocate using appropriate routine
    call data_fftbox_dealloc(row2_box)
    !call data_fftbox_dealloc(row1_box) ! JCW: Coarse grid row1_box is no
    !                                    ! JCW: longer needed
!$OMP END PARALLEL

    call timer_clock('potential_apply_to_fftbox_batch', 2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving potential_apply_to_fftbox_batch'

  end subroutine potential_apply_to_fftbox_batch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine potential_apply_to_ppd_funcs(pot_funcs_on_grid, & ! in-output
       fbasis, potential_std, grid, cell, fftbox)              ! input

    !==========================================================================!
    ! This subroutine calculates applies the potential to a set of functions   !
    ! in PPD storage.                                                          !
    !==========================================================================!
    !  Arguments:                                                              !
    !    pot_funcs_on_grid (inout) : functions in PPD format to multiply by    !
    !                                potential                                 !
    !    fbasis (input)      : The function basis type for the functions       !
    !    potential_dbl (in)  : local potential on the "dbl" sim cell grid,     !
    !                          which is actually coarse as dbl_grid_scale=1    !
    !==========================================================================!
    ! Written by Nicholas Hine on 28/02/2011.                                  !
    ! Modifid by Andrea Greco on 30/05/2015 to allow use of complex NGWFs.     !
    !==========================================================================!

    use basis, only: basis_multiply_function_by_box, basis_location_func_wrt_cell
    use cell_grid, only: GRID_INFO, cell_grid_extract_box
    use comms, only: comms_barrier, pub_my_proc_id
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, &
        data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: fbasis
    !real(kind=DP), intent(inout) :: pot_funcs_on_grid(fbasis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(inout) :: pot_funcs_on_grid
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(in) :: potential_std(grid%ld1, &
         grid%ld2, grid%max_group_slabs12)

    ! Local Variables
    integer :: ierr
    integer :: local_func
    integer :: cell_start1,cell_start2,cell_start3
    integer :: box_n1,box_n2,box_n3
    !real(kind=DP), allocatable :: potential_box(:,:,:)
    ! agrecocmplx: use new FFTBOX_DATA type
    type(FFTBOX_DATA) :: potential_box
    real(kind=DP), allocatable :: potential_buffer(:,:,:)
    logical :: i_need_box

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &potential_apply_to_ppd_funcs'

    ! Start timer
    call timer_clock('potential_apply_to_ppd_funcs',1)

    box_n1 = fbasis%maxtight_pts1
    box_n2 = fbasis%maxtight_pts2
    box_n3 = fbasis%maxtight_pts3

    ! Allocate storage for tightbox of potential and buffer
    !allocate(potential_box(box_n1,box_n2,box_n3),stat=ierr)
    !call utils_alloc_check('potential_apply_to_ppd_funcs','potential_box',ierr)
    ! agrecocmplx: allocate using appropriate routine
    ! local potential real in any case, so iscmplx=.false.
    call data_fftbox_alloc(potential_box, box_n1, box_n2, box_n3, .false.)
    allocate(potential_buffer(box_n1,box_n2,grid%max_group_slabs12), &
         stat=ierr)
    call utils_alloc_check('potential_apply_to_ppd_funcs','potential_buffer', &
         ierr)

    ! Loop over ket functions up to max on any proc
    do local_func=1,fbasis%max_on_proc

       if (local_func <= fbasis%num_on_proc(pub_my_proc_id)) then
          ! Find position of ket function wrt simulation cell:
          call basis_location_func_wrt_cell(cell_start1, &
               cell_start2,cell_start3,fbasis%tight_boxes(local_func), &
               cell)
          i_need_box = .true.
       else
          cell_start1 = -1234
          cell_start2 = -1234
          cell_start3 = -1234
          i_need_box = .false.
       end if

       ! Extract tightbox of locpot data from whole-cell grid
       ! agrecocmplx: call with argument potential_box%d
       call cell_grid_extract_box(potential_box%d,&
            potential_buffer, potential_std, grid, &
            box_n1, box_n2, box_n3, box_n1, box_n2, &
            cell_start1, cell_start2, cell_start3, i_need_box, .true.)

       ! Multiply the ket function by the locpot box directly in PPDs
       if (local_func <= fbasis%num_on_proc(pub_my_proc_id)) then
          ! agrecocmplx: call new routine in basis_new
          call basis_multiply_function_by_box(pot_funcs_on_grid, &
               potential_box, &
               fbasis%spheres(local_func), fbasis%tight_boxes(local_func), &
               1, 1, 1, fbasis%spheres(local_func)%offset, cell, fftbox)
       end if

    end do

    ! Deallocate workspace
    deallocate(potential_buffer,stat=ierr)
    call utils_dealloc_check('potential_apply_to_ppd_funcs', &
         'potential_buffer',ierr)
    !deallocate(potential_box,stat=ierr)
    !call utils_dealloc_check('potential_apply_to_ppd_funcs', &
    !     'potential_box',ierr)
    ! agrecocmplx: deallocate using appropriate routine
    call data_fftbox_dealloc(potential_box)

    ! pdh: re-sync procs
    call comms_barrier

    ! Stop timer for total
    call timer_clock('potential_apply_to_ppd_funcs',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &potential_apply_to_ppd_funcs'

  end subroutine potential_apply_to_ppd_funcs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine potential_input_to_workspace(potential_work,potential_in, &
       work_grid,in_grid)

    !==========================================================================!
    ! This subroutine transfers a potential from the input grid to a version   !
    ! on the 'workspace' grid, where all the slabs of the local proc's comms   !
    ! group are duplicated on each proc of the group. The grids may be of      !
    ! different scales, in which case the result is filtered from the input    !
    ! down to the workspace grid.                                              !
    !==========================================================================!
    !  Arguments:                                                              !
    !    potential_work (out) : Workspace array (sized for one group's slabs)  !
    !    potential_in   (in)  : Input array (sized for one proc's slabs)       !
    !    work_grid      (in)  : GRID_INFO defining workspace grid              !
    !    in_grid        (in)  : GRID_INFO defining input grid                  !
    !==========================================================================!
    ! Written by Nicholas Hine on 16/03/2011.                                  !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_allgather, pub_comms_group_size, pub_group_comm, &
         pub_my_rank_in_group
    use constants, only: DP
    use fourier, only: fourier_filter_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: work_grid
    type(GRID_INFO), intent(in) :: in_grid
    real(kind=DP), intent(out) :: potential_work(work_grid%ld1,work_grid%ld2, &
         work_grid%max_group_slabs12)
    real(kind=DP), intent(in) :: potential_in(in_grid%ld1,in_grid%ld2, &
         in_grid%max_slabs12)

    ! Local Variables
    integer :: i3_start, i3_finish
    logical :: in_is_work
    real(kind=DP), allocatable :: disjoint_src_buffer(:,:,:) ! jd: temp
    integer :: ierr

    in_is_work = .false.
    if ((work_grid%n1==in_grid%n1).and.(work_grid%n2==in_grid%n2).and. &
         (work_grid%n3==in_grid%n3)) in_is_work = .true.

    ! ndmh: get data from potential_fine to potential_dbl, depending on scales
    i3_start = work_grid%my_first_slab12_in_group
    i3_finish = work_grid%my_last_slab12_in_group
    if (in_is_work) then
       ! ndmh: copy potential from fine grid to double grid and allgather the
       ! ndmh: other procs' data within the group if required
       if (pub_comms_group_size>1) then
          call comms_allgather(potential_work,potential_in, &
               length_src=work_grid%group_lengths_slabs12(pub_my_rank_in_group), &
               lengths_dest=work_grid%group_lengths_slabs12, &
               displs_dest=work_grid%group_displs_slabs12,comm=pub_group_comm)
       else
          potential_work(:,:,i3_start:i3_finish) = &
               potential_in(:,:,1:work_grid%num_my_slabs12)
       end if
    else
       ! ndmh: filter potential from fine grid to double grid and allgather the
       ! ndmh: other procs' data within the group if required
       i3_finish = i3_start + work_grid%max_slabs12 - 1
       call fourier_filter_cell(potential_in(:,:,:), &
            potential_work(:,:,i3_start:i3_finish),in_grid,work_grid, &
            apply_nyquist=.true.)
       if (pub_comms_group_size>1) then

          ! jd: Copy slice of potential_work to prevent arg aliasing
          allocate(disjoint_src_buffer(work_grid%ld1,work_grid%ld2, &
               i3_start:i3_finish),stat=ierr)
          call utils_alloc_check('potential_input_to_workspace', &
               'disjoint_src_buffer',ierr)
          disjoint_src_buffer = potential_work(:,:,i3_start:i3_finish)

          call comms_allgather(potential_work,disjoint_src_buffer, &
               length_src=work_grid%group_lengths_slabs12(pub_my_rank_in_group), &
               lengths_dest=work_grid%group_lengths_slabs12, &
               displs_dest=work_grid%group_displs_slabs12,comm=pub_group_comm)

          deallocate(disjoint_src_buffer,stat=ierr)
          call utils_dealloc_check('potential_input_to_workspace', &
               'disjoint_src_buffer',ierr)

       end if
    end if

  end subroutine potential_input_to_workspace


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef GPU_PGI
  subroutine potential_gpu_app2_ngwf_batch(potential_fftbox_batch, &   ! output
       ngwfs_on_grid, ngwf_basis, potential_dbl, grid, fftbox, cell, & ! input
       ket_start_in_box, ket_box_start, batch_size, local_start, &     ! input
       local_end, max_current_size)                                    ! input

    !==========================================================================!
    ! This subroutine applies the local potential operator to a batch          !
    ! of NGWFs in fftboxes on GPUs. This version is used to manage the data    !
    ! transfer between device and host                                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! potential_fftbox_batch (output): Batch off fftboxes with potential-NGWFs !
    ! ngwfs_on_grid (input) : NGWFs on this proc in ppd representation         !
    ! ngwf_basis (input)    : Function basis type for NGWFs                    !
    ! potential_dbl (input) : Local potential on dbl grid in simulation cell   !
    ! batch_size (input)    : Number of fftboxes in the batch                  !
    ! local_start(input)    : First NGWF of the batch                          !
    ! local_end  (input)    : Last NGWF of the batch                           !
    !--------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 23/11/2003 for the ONETEP program.   !
    ! Modified by Nicholas Hine in 2008 to replace two-stage copy of each row  !
    ! function with single call to basis_copy_function_to_box.                 !
    ! Modified by Nicholas Hine on 21/07/2009 to use function basis type and   !
    ! to deposit result of fourier_filter straight into potential_fftbox array.!
    ! GPU version created by Karl Wilkinson in 2012. Minor restructuring needed!
    ! and may be reintegrated with the CPU code at a later date. Dynamic       !
    ! precision is currently untested in terms of its effect on a calculation. !
    !==========================================================================!

    use basis, only: basis_copy_function_to_box, basis_location_func_wrt_cell
    use cell_grid, only: GRID_INFO, cell_grid_extract_box
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc, data_set_to_zero
    use fft_box, only: FFTBOX_INFO
#ifdef GPU_SP_TEST_POT
    use fourier, only: fourier_gpu_filter_SP, fourier_gpu_interpolate_SP, &
                       fourier_gpu_filter, fourier_gpu_interpolate
#else
    use fourier, only: fourier_gpu_filter, fourier_gpu_interpolate
#endif
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use fourier_gpu_wrapper_mod
    use accel_lib !! External dependency
    use cudafor   !! External dependency

    implicit none

    ! Arguments
    integer, intent(in) :: batch_size, local_start, local_end
    integer, intent(in) :: max_current_size ! max size of batch over all procs
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FFTBOX_DATA), intent(inout) :: potential_fftbox_batch(batch_size)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    type(GRID_INFO), intent(in)  :: grid
    type(CELL_INFO), intent(in) :: cell
    real(kind=DP), intent(in) :: potential_dbl(grid%ld1, &
         grid%ld2, grid%max_group_slabs12)
    integer, intent(in) :: ket_start_in_box(3,batch_size)
    integer, intent(in) :: ket_box_start(3,batch_size)

    ! Local Variables
    !KAW merge retreat real(kind=DP), dimension(:,:,:), allocatable :: row1_box, row2_box
    type(FFTBOX_DATA) :: row1_box, row2_box
    real(kind=DP), dimension(:,:,:,:), allocatable :: potential_box_dbl_pair
    real(kind=DP), dimension(:,:,:), allocatable :: potential_box_dbl
    real(kind=DP), dimension(:,:,:), allocatable :: buffer_dbl

    integer :: ierr
    integer :: row1, row2
    integer :: n1, n2, n3, ld1, ld2
    integer :: n1_dbl, n2_dbl, n3_dbl, ld1_dbl, ld2_dbl
    integer :: iii
    integer :: batch_count, box_pos
    integer :: fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl
    integer :: prev_start1, prev_start2, prev_start3

    logical :: i_need_potential ! whether pub_my_proc_id needs pot'l in fftbox
    logical :: apply_pot_2

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering potential_gpu_app2_ngwf_batch'

    ! Start timer
    call timer_clock('potential_apply_to_ngwf_batch_gpu', 1)

    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2
    n1_dbl = fftbox%total_pt1_dbl
    n2_dbl = fftbox%total_pt2_dbl
    n3_dbl = fftbox%total_pt3_dbl
    ld1_dbl = fftbox%total_ld1_dbl
    ld2_dbl = fftbox%total_ld2_dbl

    ! ndmh: initialisations to prevent compiler warnings
    prev_start1 = -1111111
    prev_start2 = -2222222
    prev_start3 = -3333333
    fftbox_start1_dbl = -1111111
    fftbox_start2_dbl = -2222222
    fftbox_start3_dbl = -3333333

    ! Allocate workspace arrays
    !KAW retreat mergeallocate(row1_box(ld1, ld2, n3), stat=ierr)
    !KAW retreat mergecall utils_alloc_check('potential_gpu_app2_ngwf_batch','row1_box',ierr)
    !KAW retreat mergeallocate(row2_box(ld1, ld2, n3), stat=ierr)
    !KAW retreat mergecall utils_alloc_check('potential_gpu_app2_ngwf_batch','row2_box',ierr)
    call data_fftbox_alloc(row1_box, ld1, ld2, n3, &
         iscmplx=ngwfs_on_grid%iscmplx)
    call data_fftbox_alloc(row2_box, ld1, ld2, n3, &
         iscmplx=ngwfs_on_grid%iscmplx)



    allocate(potential_box_dbl(ld1_dbl, ld2_dbl, 2*n3), stat=ierr)
    call utils_alloc_check('potential_gpu_app2_ngwf_batch',&
         'potential_box_dbl',ierr)
    allocate(potential_box_dbl_pair(ld1_dbl, ld2_dbl, 2*n3, 2), stat=ierr)
    call utils_alloc_check('potential_gpu_app2_ngwf_batch',&
         'potential_box_dbl_pair',ierr)
    allocate(buffer_dbl(ld1_dbl, ld2_dbl, grid%max_group_slabs12), &
         stat=ierr)
    call utils_alloc_check('potential_gpu_app2_ngwf_batch','buffer_dbl',&
         ierr)

    ! cks: initialisation
    buffer_dbl = 0.0_DP
    potential_box_dbl = 0.0_DP



    batch_count = 0
    do iii=local_start, local_start+max_current_size-1, 2
       row1 = iii ; row2 = iii + 1

       if (row1 <= local_end) then

          ! ndmh: use new basis function copying straight to fftbox from ppds
          call basis_copy_function_to_box(row1_box, &
               ket_start_in_box(1,iii-local_start+1), &
               ket_start_in_box(2,iii-local_start+1), &
               ket_start_in_box(3,iii-local_start+1), &
               ngwf_basis%tight_boxes(row1), ngwfs_on_grid, &
               ngwf_basis%spheres(row1), cell, fftbox)

       else
          ! agrecocmplx
          call data_set_to_zero(row1_box)
          !KAW retreatmerge row1_box = 0.0_DP
       end if

       if (row2 <= local_end) then

          ! ndmh: use new basis function copying straight to fftbox from ppds
          call basis_copy_function_to_box(row2_box, &
               ket_start_in_box(1,iii-local_start+2), &
               ket_start_in_box(2,iii-local_start+2), &
               ket_start_in_box(3,iii-local_start+2), &
               ngwf_basis%tight_boxes(row2), ngwfs_on_grid, &
               ngwf_basis%spheres(row2), cell, fftbox)

       else
          ! agrecocmplx
          call data_set_to_zero(row2_box)
          !KAW retreatmerge row2_box = 0.0_DP

       end if

       ! cks: Put potential into FFTbox for first NGWF
       box_pos = iii-local_start+1
       fftbox_start1_dbl = 2*ket_box_start(1,box_pos) - 1
       fftbox_start2_dbl = 2*ket_box_start(2,box_pos) - 1
       fftbox_start3_dbl = 2*ket_box_start(3,box_pos) - 1

       ! cks:----- Get potential for FIRST NGWF if needed ------------
       i_need_potential = ((fftbox_start1_dbl /= prev_start1 .or. &
            fftbox_start2_dbl /= prev_start2 .or. &
            fftbox_start3_dbl /= prev_start3) .and. &
            row1 <= local_end)

       call cell_grid_extract_box(potential_box_dbl(:,:,:),&
            buffer_dbl, potential_dbl, grid, n1_dbl, n2_dbl, n3_dbl, &
            ld1_dbl, ld2_dbl, fftbox_start1_dbl, fftbox_start2_dbl, &
            fftbox_start3_dbl, i_need_potential, .true.)
       potential_box_dbl_pair(:,:,:,1) = potential_box_dbl(:,:,:)

       if (i_need_potential) then
          prev_start1 = fftbox_start1_dbl
          prev_start2 = fftbox_start2_dbl
          prev_start3 = fftbox_start3_dbl
       end if
       ! cks:-- END Get potential for FIRST NGWF if needed ------------

       ! cks: Put potential into FFTbox for second NGWF
       box_pos = iii-local_start+2
       fftbox_start1_dbl = 2*ket_box_start(1,box_pos) - 1
       fftbox_start2_dbl = 2*ket_box_start(2,box_pos) - 1
       fftbox_start3_dbl = 2*ket_box_start(3,box_pos) - 1

       ! cks:----- Get potential for SECOND NGWF if needed -----------
       i_need_potential = ((fftbox_start1_dbl /= prev_start1 .or. &
            fftbox_start2_dbl /= prev_start2 .or. &
            fftbox_start3_dbl /= prev_start3) .and. &
            row2 <= local_end)

       call cell_grid_extract_box(potential_box_dbl(:,:,:), &
            buffer_dbl, potential_dbl, grid, n1_dbl, n2_dbl, n3_dbl, &
            ld1_dbl, ld2_dbl, fftbox_start1_dbl, fftbox_start2_dbl, &
            fftbox_start3_dbl, i_need_potential, .true.)
       potential_box_dbl_pair(:,:,:,2) = potential_box_dbl(:,:,:)

       if (i_need_potential) then
          prev_start1 = fftbox_start1_dbl
          prev_start2 = fftbox_start2_dbl
          prev_start3 = fftbox_start3_dbl
       endif
       ! cks:-- END Get potential for SECOND NGWF if needed -----------

       apply_pot_2 = .false.
       if (row1 <= local_end) then
          if (row2 <= local_end) apply_pot_2 = .true.

          ! cks: interpolate row1 and row2 ngwfs
#ifdef GPU_SP_TEST_POT
          ! kaw: single precision version of the routines
          if (GPU_SP) then
             call fourier_gpu_interpolate_SP(row1_box, row2_box, &
                     potential_box_dbl_pair, apply_pot_2)
          else
             call fourier_gpu_interpolate(row1_box, row2_box, &
                     potential_box_dbl_pair, apply_pot_2)
          end if
#else


             call fourier_gpu_interpolate(row1_box%d, row2_box%d, &
                     potential_box_dbl_pair, apply_pot_2)
#endif

       else

#ifdef GPU_SP_TEST_POT
          if (GPU_SP) then
             fftbox_fine_1_gpuSP = 0.0
             fftbox_fine_2_gpuSP = 0.0
          else
             fftbox_fine_1_gpu = 0.0_DP
             fftbox_fine_2_gpu = 0.0_DP
          end if
#else
             fftbox_fine_1_gpu = 0.0_DP
             fftbox_fine_2_gpu = 0.0_DP
#endif
       end if

       ! ndmh: deposit one result of interpolation in potential_fftbox_batch
       ! ndmh: and discard other result
       if ((row1 <= local_end).and.(row2 > local_end)) then
          batch_count = batch_count + 1
          ! Filter result to standard grid
#ifdef GPU_SP_TEST_POT
          if (GPU_SP) then
             call fourier_gpu_filter_SP(batch_count, .false.)
          else
             call fourier_gpu_filter(batch_count, .false.)
          end if
#else
             call fourier_gpu_filter(batch_count, .false.)
#endif

       ! ndmh: deposit both results of interpolation in potential_fftbox_batch
       else if ((row1 <= local_end).and.(row2 <= local_end)) then
          batch_count = batch_count + 2
          ! Filter result to standard grid
#ifdef GPU_SP_TEST_POT
          if (GPU_SP) then
             call fourier_gpu_filter_SP(batch_count, .true.)
          else
             call fourier_gpu_filter(batch_count, .true.)
          end if
#else
             call fourier_gpu_filter(batch_count, .true.)
#endif
       end if
    end do


    deallocate(buffer_dbl, stat=ierr)
    call utils_dealloc_check('potential_gpu_app2_ngwf_batch','buffer_dbl',&
         ierr)
    deallocate(potential_box_dbl, stat=ierr)
    call utils_dealloc_check('potential_gpu_app2_ngwf_batch',&
         'potential_box_dbl',ierr)

    !KAW retreat mergedeallocate(row2_box, stat=ierr)
    !KAW retreat mergecall utils_dealloc_check('potential_gpu_app2_ngwf_batch','row2_box', &
    !KAW retreat merge        ierr)
    !KAW retreat mergedeallocate(row1_box, stat=ierr)
    !KAW retreat call utils_dealloc_check('potential_gpu_app2_ngwf_batch','row1_box', &
    !KAW retreat         ierr)
    call data_fftbox_dealloc(row2_box)
    call data_fftbox_dealloc(row1_box)

    call timer_clock('potential_apply_to_ngwf_batch_gpu', 2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving potential_gpu_app2_ngwf_batch'

  end subroutine potential_gpu_app2_ngwf_batch

#else
  subroutine dummy_gpu
  use fourier_gpu_wrapper_mod
  !kaw: This is a dummy use statement to make the check_dependencies script happy.
  end subroutine dummy_gpu
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module potential
