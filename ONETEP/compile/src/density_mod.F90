! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Chris-Kriton Skylaris, Arash A. Mostofi, Peter D. Haynes,
!   Nicholas D.M. Hine and Karl Wilkinson
!
!   between January 2001 and May 2014.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module density

  use constants, only: DP
  use rundat, only: pub_debug_on_root

  implicit none

  private

  public :: density_on_grid
  public :: density_plot_slice
  public :: density_workspace_to_output ! required by ke_density module

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_on_grid(density_fine, &
       fine_grid, dbl_grid, cell, fftbox, denskern, overlap, &
       ngwfs_on_grid_left, ngwf_basis_left, ngwfs_on_grid_right, &
       ngwf_basis_right)

    !==========================================================================!
    ! This subroutine calculates the charge density of a fine grid in the      !
    ! simulation cell. Wrapper for calculating the density on the fine grid    !
    ! (which may be any scale), by calling the appropriate routine to generate !
    ! the density on the coarse or double grids and then interpolating up (if  !
    ! necessary) to the grid scale of the fine grid.                           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! density_fine (output) : The total charge density of the system on the    !
    ! fine grid of the simulation cell                                         !
    ! fine_grid     (input) : GRID_INFO describing fine grid on which density  !
    !   is eventually constructed (may be finer than twice the psinc cutoff).  !
    ! dbl_grid      (input) : GRID INFO structure describing a grid with twice !
    !  the cutoff of the psinc grid, which is enough to represent the soft     !
    !  part of the density obtained from the NGWF products                     !
    ! denskern      (input) : Density kernel sparse matrix                     !
    ! overlap       (input) : Overlap matrix                                   !
    ! ngwfs_on_grid (input) : The NGWFs on the real-space ppd representation   !
    !  that belong to this proc (left and right versions to allow construction !
    !  of transition densities in LRTDDFT in form phi_a K^ab chi_b             !
    ! ngwf_basis            : The function basis for the NGWFs (also allows    !
    !  separate specification of left and right versions)                      !
    !--------------------------------------------------------------------------!
    ! Current form by Nicholas Hine 02/28/11. Old code moved to new routine    !
    ! density_on_dbl_grid (see below).                                         !
    ! Modified by Andrea Greco on 08/05/2015 to allow use of complex NGWFs.    !
    !==========================================================================!

    ! agrecocmplx
    use datatypes, only: FUNCTIONS
    use cell_grid, only: GRID_INFO
    use constants, only: DP
    use function_basis, only: FUNC_BASIS
    use fft_box, only: FFTBOX_INFO
    use rundat, only: pub_dbl_is_std, pub_num_spins, pub_check_density, &
        pub_imag_thr
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    ! JCW: For outputting binary density dumps for input to CASTEP
    use rundat, only: pub_devel_code, pub_rootname
    use utils, only: utils_abort, utils_devel_code

    implicit none

    ! Arguments
    ! agrecokpt: change this to denskern(pub_num_spins, pub_num_kpoints)
    ! or SPAM3_ARRAY denskern
    ! rc2013: EMBED_FIX! problems with deferred shape arrays when
    ! passing embedding structures
    type(SPAM3), intent(in) :: denskern(pub_num_spins,1)
    type(SPAM3), intent(in) :: overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis_left
    type(FUNC_BASIS), intent(in) :: ngwf_basis_right
    type(GRID_INFO), intent(in) :: fine_grid
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(out) :: density_fine(fine_grid%ld1, &
         fine_grid%ld2, fine_grid%max_slabs12, pub_num_spins)
    ! agrecokpt: need to change this to ngwf_on_grid_left(pub_num_kpoints)
    ! and ngwf_on_grid_right(pub_num_kpoints)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid_left
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid_right

    ! Local Variables
    integer :: ierr
    logical :: loc_cmplx
    real(kind=DP), allocatable, dimension(:,:,:,:) :: density_work
    ! agrecocmplx: array to store imaginary part of originally complex
    ! density, to check the density is real when complex NGWFs are used
    real(kind=DP), allocatable, dimension(:,:,:,:) :: density_work_imag

    ! agrecocmplx
    loc_cmplx = overlap%iscmplx

    ! ndmh: if the output is to be on the standard grid, use that as workspace
    if (pub_dbl_is_std) then

       ! ndmh: allocate coarse-grid whole cell workspace array
       allocate(density_work(dbl_grid%ld1,dbl_grid%ld2,&
            dbl_grid%max_group_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('density_on_grid','density_work',ierr)

       ! agrecocmplx: if complex NGWFs and request to check density
       if (loc_cmplx.and.pub_check_density) then
          allocate(density_work_imag(dbl_grid%ld1,dbl_grid%ld2,&
             dbl_grid%max_group_slabs12,pub_num_spins),stat=ierr)
          call utils_alloc_check('density_on_grid','density_work_imag',ierr)

          ! agrecocmplx: calculate the density on the coarse grid
          ! also accumulate imaginary part to check density is real
          call density_on_coarse_grid(density_work, dbl_grid, cell, fftbox, &
               fftbox, denskern, overlap, ngwfs_on_grid_left, &
               ngwf_basis_left, ngwfs_on_grid_right, ngwf_basis_right, &
               density_work_imag)

          ! agrecocmplx: abort if imaginary part of density
          ! is above threshold; this usually means something has
          ! gone wrong..... need to check if density needs to be
          ! real at this stage even in the case of PAW formalism
          if (any(abs(density_work_imag)>pub_imag_thr)) then
             call utils_abort('Error in routine density_on_grid: &
                  &imaginary part of (complex) density is above threshold')
          end if

          deallocate(density_work_imag,stat=ierr)
          call utils_dealloc_check('density_on_grid','density_work_imag',ierr)

       else
          ! ndmh: calculate the density on the coarse grid
          call density_on_coarse_grid(density_work, dbl_grid, cell, fftbox, &
               fftbox, denskern, overlap, ngwfs_on_grid_left, &
               ngwf_basis_left, ngwfs_on_grid_right, ngwf_basis_right)
       end if

       ! ndmh: copy or upscale to output grid as required
       call density_workspace_to_output(density_work,density_fine, &
            dbl_grid,fine_grid)

       deallocate(density_work,stat=ierr)
       call utils_dealloc_check('density_on_grid','density_work',ierr)

    else ! ndmh: otherwise, default is to use dbl grid as workspace

       ! ndmh: allocate double-grid whole cell workspace array
       allocate(density_work(dbl_grid%ld1,dbl_grid%ld2,&
            dbl_grid%max_group_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check('density_on_grid','density_work',ierr)

       ! agrecocmplx: if complex NGWFs and request to check density
       if (loc_cmplx.and.pub_check_density) then
          allocate(density_work_imag(dbl_grid%ld1,dbl_grid%ld2,&
             dbl_grid%max_group_slabs12,pub_num_spins),stat=ierr)
          call utils_alloc_check('density_on_grid','density_work_imag',ierr)

          ! agrecocmplx: calculate the density on the double grid
          ! also accumulate imaginary part to check density is real
          call density_on_dbl_grid(density_work, dbl_grid, cell, fftbox, &
               denskern, overlap, ngwfs_on_grid_left, ngwf_basis_left, &
               ngwfs_on_grid_right, ngwf_basis_right, density_work_imag)

          ! agrecocmplx: abort if imaginary part of density
          ! is above threshold; this usually means something has
          ! gone wrong..... need to check if density needs to be
          ! real at this stage even in the case of PAW formalism
          if (any(abs(density_work_imag)>pub_imag_thr)) then
             call utils_abort('Error in routine density_on_grid: &
                  &imaginary part of (complex) density is above threshold')
          end if

          deallocate(density_work_imag,stat=ierr)
          call utils_dealloc_check('density_on_grid','density_work_imag',ierr)

       else
          ! ndmh: calculate the density on the double grid
          call density_on_dbl_grid(density_work, dbl_grid, cell, fftbox, &
               denskern, overlap, ngwfs_on_grid_left, ngwf_basis_left, &
               ngwfs_on_grid_right, ngwf_basis_right)
       end if

       ! ndmh: copy or upscale to 'fine' grid as required
       call density_workspace_to_output(density_work,density_fine, &
            dbl_grid,fine_grid)


       deallocate(density_work,stat=ierr)
       call utils_dealloc_check('density_on_grid','density_work',ierr)

    end if

    ! JCW: castep_dump is true, then dump the charge density in CASTEP-compatible
    ! JCW: binary format to disk. This is set by DEVEL_CODE: DEN:CASTEP_DUMP=T:DEN
    if (utils_devel_code(.false.,"DEN","CASTEP_DUMP",pub_devel_code)) then
       ! The definition of rho (electronic or total density) depends on what
       ! keywords are used to run ONETEP. See the documentation for this routine
       ! (above) for details.
       ! The density on the full grid (including parts that are truncated when
       ! dimensions are set to multigrid magic numbers) is dumped.

       if (pub_num_spins == 1) then
          call castep_density_dump(density_fine(:,:,:,1),&
               trim(pub_rootname)//"_castep_dump.den",cell,fine_grid)
       else if (pub_num_spins == 2) then
          call castep_density_dump(density_fine(:,:,:,1)+density_fine(:,:,:,2),&
               trim(pub_rootname)//"_castep_dump.den",cell,fine_grid)
       else
          call utils_abort("Error in density_on_grid: pub_num_spins must be &
               &1 or 2 for CASTEP density dump to proceed.")
       end if
    end if

  end subroutine density_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_workspace_to_output(density_work,density_out, &
       work_grid,out_grid)

    !==========================================================================!
    ! This subroutine transfers a density from the workspace grid to a version !
    ! on the output grid. On the workspace grid all the slabs of the local     !
    ! proc's comms group are duplicated on each proc of the group. The grids   !
    ! may be of different scales, in which case the result is interpolated     !
    ! from the workspace up to the output grid.                                !
    !==========================================================================!
    !  Arguments:                                                              !
    !    density_work (in)  : Workspace array (sized for one group's slabs)    !
    !    density_out  (out) : Output array (sized for one proc's slabs)        !
    !    work_grid    (in)  : GRID_INFO defining workspace grid                !
    !    out_grid     (in)  : GRID_INFO defining output grid                   !
    !==========================================================================!
    ! Written by Nicholas Hine on 16/03/2011.                                  !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce, pub_comms_group_size, &
         pub_group_comm, pub_first_proc_in_group
    use constants, only: DP
    use fourier, only: fourier_interpolate_cell
    use rundat, only: pub_num_spins

    ! Arguments
    type(GRID_INFO), intent(in) :: work_grid
    type(GRID_INFO), intent(in) :: out_grid
    real(kind=DP), intent(inout) :: density_work(work_grid%ld1,work_grid%ld2, &
         work_grid%max_group_slabs12,pub_num_spins)
    real(kind=DP), intent(out) :: density_out(out_grid%ld1,out_grid%ld2, &
         out_grid%max_slabs12,pub_num_spins)

    ! Local Variables
    integer :: is
    integer :: proc
    integer :: i3_start, i3_finish
    logical :: work_is_out

    work_is_out = .false.
    if ((work_grid%n1==out_grid%n1).and.(work_grid%n2==out_grid%n2).and. &
         (work_grid%n3==out_grid%n3)) work_is_out = .true.

    ! ndmh: transfer data from density_work to density_out
    do is=1,pub_num_spins
       if (pub_comms_group_size>1) then
          ! ndmh: sum density data over all procs in this group
          do proc=0,pub_comms_group_size-1
             i3_start = work_grid%first_slab12(proc + &
                  pub_first_proc_in_group) - &
                  work_grid%first_slab12(pub_first_proc_in_group) + 1
             i3_finish = work_grid%last_slab12(proc + &
                  pub_first_proc_in_group) - &
                  work_grid%first_slab12(pub_first_proc_in_group) + 1
             ! ndmh: if workspace grid is output grid, sum result
             ! ndmh: directly in the density_out array, otherwise
             ! ndmh: sum it in-place and transfer
             if (work_is_out) then
                call comms_reduce('SUM',density_out(:,:, &
                     1:(i3_finish-i3_start+1),is), &
                     comm=pub_group_comm,root=proc, &
                     d_array_src=density_work(:,:,i3_start:i3_finish,is))
             else
                call comms_reduce('SUM',density_work(:,:, &
                     i3_start:i3_finish,is),comm=pub_group_comm,root=proc)
             end if
          end do
       else ! ndmh: no summation required - copy straight to fine grid
          if (work_is_out) then
             i3_start = work_grid%my_first_slab12_in_group
             i3_finish = work_grid%my_last_slab12_in_group
             density_out(:,:,1:work_grid%num_my_slabs12,is) = &
                  density_work(:,:,i3_start:i3_finish,is)
          end if
       end if
       ! ndmh: now interpolate from workspace grid to output grid if required
       if (.not.work_is_out) then
          i3_start = work_grid%my_first_slab12_in_group
          i3_finish = i3_start + work_grid%max_slabs12 - 1
          call fourier_interpolate_cell(density_work(:,:,i3_start:i3_finish,is), &
               density_out(:,:,:,is),work_grid,out_grid,apply_nyquist=.true.)
       end if
       ! ndmh: zero any padding in density_out
       density_out(:,:,(out_grid%num_my_slabs12+1):out_grid%max_slabs12,is) &
            = 0.0_DP

    end do

  end subroutine density_workspace_to_output


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_on_coarse_grid(density_std, grid, cell, fftbox, &
       uni_tightbox, denskern, overlap, ngwfs_on_grid_left, ngwf_basis_left, &
       ngwfs_on_grid_right, ngwf_basis_right, density_std_imag)

    !==========================================================================!
    ! This subroutine calculates the charge density on the coarse grid in the  !
    ! simulation cell.                                                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! density  (output) : The total charge density of the system on the        !
    ! coarse grid of the simulation cell                                       !
    ! denskern      (input) : Density kernel sparse matrix                     !
    ! overlap       (input) : Overlap matrix                                   !
    ! ngwfs_on_grid (input) : The NGWFs on the real-space ppd representation   !
    !  that belong to this proc                                                !
    ! ngwf_basis            : The function basis for the NGWFs                 !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 28/02/2011.                                  !
    ! Modified by Andrea Greco on 08/05/2015 to allow use of complex NGWFs.    !
    !==========================================================================!

    use basis, only: basis_location_func_wrt_cell, basis_copy_function_to_box
    use cell_grid, only: GRID_INFO, cell_grid_deposit_box
    use comms, only: comms_barrier, pub_my_proc_id
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, &
         data_functions_alloc, data_functions_dealloc, &
         data_fftbox_alloc, data_fftbox_dealloc, &
         data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_sum_ppd_funcs
    use rundat, only: pub_num_spins, pub_imag_thr
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_abort

    implicit none

    ! Arguments
    ! agrecokpt: change to denskern(pub_num_spins, pub_num_kpoints)
    ! or SPAM3_ARRAY denskern
    type(SPAM3), intent(in) :: denskern(pub_num_spins)
    type(SPAM3), intent(in) :: overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis_left
    type(FUNC_BASIS), intent(in) :: ngwf_basis_right
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FFTBOX_INFO), intent(in) :: uni_tightbox
    real(kind=DP), intent(out) :: density_std(grid%ld1, &
         grid%ld2, grid%max_group_slabs12, pub_num_spins)
    ! agrecokpt: need to change this to ngwf_on_grid_left(pub_num_kpoints)
    ! and ngwf_on_grid_right(pub_num_kpoints)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid_left
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid_right
    ! agrecocmplx: accumulate imaginary part of density if required
    real(kind=DP), intent(out), optional :: density_std_imag(grid%ld1, &
         grid%ld2, grid%max_group_slabs12, pub_num_spins)

    ! Local Variables
    integer :: ierr
    integer :: is
    integer :: local_col
    integer :: col_cell_start1,col_cell_start2,col_cell_start3
    integer :: box_n1,box_n2,box_n3
    ! agrecokpt: need to introduce partial sum for
    ! k-point dependence, i.e. row_sum_on_grid_partial(ik),
    ! then row_sum_on_grid is obtained by doing the weighted
    ! sum of all row_sum_on_grid_partial(ik)
    type(FUNCTIONS) :: row_sum_on_grid(pub_num_spins)
    type(FFTBOX_DATA) :: density_box
    real(kind=DP), allocatable :: density_buffer(:,:,:)
    logical :: i_have_box
    ! agrecocmplx: local variable for complex quantities
    logical :: loc_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering density_on_coarse_grid'

    call timer_clock('density_on_coarse_grid', 1)

    box_n1 = uni_tightbox%total_pt1
    box_n2 = uni_tightbox%total_pt2
    box_n3 = uni_tightbox%total_pt3

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid_left%iscmplx

    ! Allocate temporary storage
    !allocate(row_sum_on_grid(ngwf_basis_left%size_on_grid, &
    !     pub_num_spins),stat=ierr)
    !call utils_alloc_check('density_on_coarse_grid','row_sum_on_grid',ierr)
    ! agrecocmplx: allocate using appropriate routine and initialise
    ! agrecokpt: for each kpt, need to allocate row_sum_on_grid_partial(ik)
    ! as well
    do is=1,pub_num_spins
       call data_functions_alloc(row_sum_on_grid(is), &
            ngwf_basis_left%size_on_grid, iscmplx=loc_cmplx)
       call data_set_to_zero(row_sum_on_grid(is))
    end do
    !allocate(density_box(box_n1,box_n2,box_n3),stat=ierr)
    !call utils_alloc_check('density_on_coarse_grid','density_box',ierr)
    ! agrecocmplx: allocate using appropriate routine
    call data_fftbox_alloc(density_box, box_n1, box_n2, &
                            box_n3, iscmplx=loc_cmplx)
    allocate(density_buffer(box_n1,box_n2,grid%max_group_slabs12), &
         stat=ierr)
    call utils_alloc_check('density_on_coarse_grid','density_buffer',ierr)

    !row_sum_on_grid = 0.0_DP
    !density_box = 0.0_DP
    ! agrecocmplx: initialise accordingly
    call data_set_to_zero(density_box)
    density_std = 0.0_DP

    ! agrecocmplx
    if (present(density_std_imag)) then
       density_std_imag = 0.0_DP
    end if

    ! Calculate \sum_b K^ab \phi_b(r)
    ! agrecokpt: this needs to be called for each row_sum_on_grid_partial(ik),
    ! since in general we will have different denskern/overlap/NGWFs for different
    ! k-points; row_sum_on_grid_partial(ik), denskern, overlap, ngwfs_on_grid_right
    ! all depend on the specific kpoint kpt
    call function_ops_sum_ppd_funcs(row_sum_on_grid,ngwf_basis_left,denskern, &
         1,pub_num_spins,overlap,ngwfs_on_grid_right,ngwf_basis_right)

    do is=1,pub_num_spins

       ! Multiply by \phi_a(r)
       !row_sum_on_grid(:,is) = row_sum_on_grid(:,is)*ngwfs_on_grid_left(:)
       ! agrecocmplx: if complex NGWFs, take complex conjugate of bra
       ! agrecokpt: in general this must be done for each k-point,
       ! using row_sum_on_grid_partial(ik) and ngwfs_on_grid_left(ik)
       if (loc_cmplx) then
           row_sum_on_grid(is)%z(:) = &
               row_sum_on_grid(is)%z(:)*conjg(ngwfs_on_grid_left%z(:))
       else
           row_sum_on_grid(is)%d(:) = &
               row_sum_on_grid(is)%d(:)*ngwfs_on_grid_left%d(:)
       end if

       ! Loop over col functions up to max on any proc
       do local_col=1,ngwf_basis_left%max_on_proc

          ! If there are any cols left on this proc...
          if (local_col <= ngwf_basis_left%num_on_proc(pub_my_proc_id)) then

             ! Copy density for this column to box
             !call basis_copy_function_to_box(density_box,box_n1,box_n2,box_n3, &
             !     1,1,1,ngwf_basis_left%tight_boxes(local_col), &
             !     row_sum_on_grid(:,is),ngwf_basis_left%spheres(local_col), &
             !     cell, fftbox)
             ! agrecocmplx: use new version of routine in basis_new_mod
             call basis_copy_function_to_box(density_box, &
                  1,1,1,ngwf_basis_left%tight_boxes(local_col), &
                  row_sum_on_grid(is),ngwf_basis_left%spheres(local_col), &
                  cell, fftbox)

             ! Find position of col function wrt simulation cell:
             call basis_location_func_wrt_cell(col_cell_start1, &
                  col_cell_start2,col_cell_start3, &
                  ngwf_basis_left%tight_boxes(local_col),cell)

             i_have_box = .true.
          else
             col_cell_start1 = -1234
             col_cell_start2 = -1234
             col_cell_start3 = -1234
             i_have_box = .false.
          end if

          ! Deposit box of density to whole-cell grid
          !call cell_grid_deposit_box(density_std(:,:,:,is),&
          !     density_box, density_buffer, grid, &
          !     box_n1, box_n2, box_n3, box_n1, box_n2, &
          !     col_cell_start1, col_cell_start2, col_cell_start3, &
          !     i_have_box, .true.)
          ! agrecocmplx: if using complex NGWFs, take real part of
          ! density_box%z, since the density must be real
          if (loc_cmplx) then
              call cell_grid_deposit_box(density_std(:,:,:,is),&
               real(density_box%z), density_buffer, grid, &
               box_n1, box_n2, box_n3, box_n1, box_n2, &
               col_cell_start1, col_cell_start2, col_cell_start3, &
               i_have_box, .true.)

               ! agrecocmplx: accumulate imaginary part if required
               if (present(density_std_imag)) then
                  call cell_grid_deposit_box(density_std_imag(:,:,:,is),&
                   aimag(density_box%z), density_buffer, grid, &
                   box_n1, box_n2, box_n3, box_n1, box_n2, &
                   col_cell_start1, col_cell_start2, col_cell_start3, &
                   i_have_box, .true.)

               end if
          else
              call cell_grid_deposit_box(density_std(:,:,:,is),&
               density_box%d, density_buffer, grid, &
               box_n1, box_n2, box_n3, box_n1, box_n2, &
               col_cell_start1, col_cell_start2, col_cell_start3, &
               i_have_box, .true.)
          end if

       end do

    end do

    ! Deallocate workspace
    deallocate(density_buffer,stat=ierr)
    call utils_dealloc_check('density_on_coarse_grid','density_buffer',ierr)
    !deallocate(density_box,stat=ierr)
    !call utils_dealloc_check('density_on_coarse_grid','density_box',ierr)
    ! agrecocmplx: deallocate using appropriate routine
    call data_fftbox_dealloc(density_box)
    ! agrecocmplx: deallocate using appropriate routine
    ! agrecokpt: deallocate row_sum_on_grid_partial as well
    do is=1,pub_num_spins
        call data_functions_dealloc(row_sum_on_grid(is))
    end do
    !deallocate(row_sum_on_grid,stat=ierr)
    !call utils_dealloc_check('density_on_coarse_grid','row_sum_on_grid',ierr)

    ! Re-sync procs
    call comms_barrier

    call timer_clock('density_on_coarse_grid', 2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving density_on_coarse_grid'

  end subroutine density_on_coarse_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_on_dbl_grid(density_dbl, dbl_grid, cell, fftbox, &
       denskern, overlap, ngwfs_on_grid_left, ngwf_basis_left, &
       ngwfs_on_grid_right, ngwf_basis_right, density_dbl_imag)

    !==========================================================================!
    ! This subroutine calculates the charge density on the double grid in the  !
    ! simulation cell.                                                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! density_dbl  (output) : The total charge density of the system on the    !
    ! double grid of the simulation cell                                       !
    ! denskern      (input) : Density kernel sparse matrix                     !
    ! overlap       (input) : Overlap matrix                                   !
    ! ngwfs_on_grid (input) : The NGWFs on the real-space ppd representation   !
    !  that belong to this proc                                                !
    ! ngwf_basis            : The function basis for the NGWFs                 !
    !--------------------------------------------------------------------------!
    ! Key internal variables:                                                  !
    !   batch_size: This is set equal to pub_fftbox_batch_size which comes from   !
    !     the rundat module. The value of batch size determines how large is   !
    !     the batch of accumulated fftboxes. Increasing this number decreases  !
    !     the communication per processor but increases the allocated memory   !
    !     per processor.                                                       !
    !--------------------------------------------------------------------------!
    ! Comment on parallel efficiency:                                          !
    ! This code for computing the charge density in parallel has good parallel !
    ! scalability because it adheres to two important conditions:              !
    ! 1) Interpolation, multiplication and deposition to the fine-grid         !
    !    charge density in the simulation cell happens only once per           !
    !    column of the density kernel.                                         !
    ! 2) Only NGWFs in ppd representation are communicated between the         !
    !    processors in one-to-one non-blocking fashion. This is crucial        !
    !    as communication of NGWFs in fftboxes rather than ppd representation  !
    !    is too costly and destroys parallel scalability.                      !
    ! The above goals are achieved by calculating the charge density in two    !
    ! stages: First density_batch_row_sums (involves one-to-one communication) !
    ! is called to accumulate sums of rows for a batch of columns of the       !
    ! density kernel. Then the subroutine density_batch_interp_deposit         !
    ! (involves no communication) is called to interpolate the accumulated     !
    ! sums and their columns and deposit them in the fine grid of the          !
    ! simulation cell to build the charge density.                             !
    !--------------------------------------------------------------------------!
    ! Originaly written by Chris-Kriton Skylaris in January 2001,              !
    ! to use a "pair-box".                                                     !
    ! Improved by Oswaldo Dieguez in 2001 so that it "cshifts" the "pair-box"  !
    ! Rewritten by Arash Mostofi in 2003 so that it uses a "triple-box" and    !
    ! interpolates only once per row of the density kernel.                    !
    ! Rewritten and parallelised by Chris-Kriton Skylaris in November 2003.    !
    ! Improved parallel version written by Chris-Kriton Skylaris on 15/1/2004. !
    ! Modified by Peter Haynes, July 2006 to use parallel SPAM 2.              !
    ! Spin polarised by Peter Haynes, July 2006                                !
    ! Modified in July 2009 by Nicholas Hine to use function_basis type and    !
    ! function sum routines.                                                   !
    ! Created in current form by Nicholas Hine on 28/02/2011 - Code all comes  !
    ! from previous version of density_on_fine_grid.                           !
    ! Modified to incorporate calls to GPU routines and dynamic precision by   !
    ! Karl Wilkinson in 2012                                                   !
    ! Modified by Andrea Greco on 08/05/2015 to allow use of complex NGWFs.    !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero, &
         data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_interpolate_cell
! Dynamic Precision
#ifdef GPU_SP_TEST
    use fourier_gpu_wrapper_mod, only: GPU_SP, GPU_SP_switch, deltaE_gpu
#endif
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_sum_fftbox_batch, &
         function_ops_batch_col_start
    use rundat, only: pub_fftbox_batch_size, pub_num_spins
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    ! agrecokpt: need to use denskern(pub_num_spins,pub_num_kpoints)
    type(SPAM3), intent(in) :: denskern(pub_num_spins,1)
    type(SPAM3), intent(in) :: overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis_left
    type(FUNC_BASIS), intent(in) :: ngwf_basis_right
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(out) :: density_dbl(dbl_grid%ld1, &
         dbl_grid%ld2, dbl_grid%max_group_slabs12, pub_num_spins)
    ! agrecokpt: need to change to ngwfs_on_grid_left(pub_num_kpoints),
    ! ngwfs_on_grid_right(pub_num_kpoints)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid_left
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid_right
    ! agrecocmplx: accumulate imaginary part of density if required
    real(kind=DP), intent(out), optional :: density_dbl_imag(dbl_grid%ld1, &
       dbl_grid%ld2, dbl_grid%max_group_slabs12, pub_num_spins)

    ! Local Variables
    integer :: batch_size
    integer :: max_current_size
    integer :: batch_count
    integer :: n_batches
    integer :: local_start, local_end
    integer :: ierr    ! pdh: error flag
    integer :: idx_len
    integer, allocatable :: overlap_idx(:)            ! Index for overlap matrix
    integer, allocatable :: col_box_start(:,:)
    integer, allocatable :: col_start_in_box(:,:)
    ! agrecocmplx: extra variables to allocate row_fftbox_sum_batch
    integer :: is, batch
    ! agrecocmplx: extra variable to use complex quantities
    logical :: loc_cmplx

#ifdef GPU_PGI
    ! agrecocmplx: these are not yet compatible with complex NGWFs
    ! agrecokpt: need to add partial sums for each kpoint
    !type(FFTBOX_DATA), allocatable, dimension(:,:) :: row_fftbox_sum_batch
    type(FFTBOX_DATA), allocatable, pinned, dimension(:,:) :: row_fftbox_sum_batch
    !KAW merge retreat real(kind=DP), allocatable, pinned, dimension(:,:,:,:,:) :: row_fftbox_sum_batch
    real(kind=DP), allocatable, pinned, dimension(:,:,:,:,:) :: row_fftbox_sum_batch_gpu

    real(kind=DP), allocatable, pinned, dimension(:,:,:,:,:) :: row_fftbox_dbl
    !type(FFTBOX_DATA) :: col_fftbox
    type(FFTBOX_DATA), pinned :: col_fftbox
    real(kind=DP), allocatable, pinned, dimension(:,:,:,:) :: final_row_fftbox_dbl
    real(kind=DP), allocatable, pinned, dimension(:,:,:,:) :: buffer_dbl
#else
    type(FFTBOX_DATA), allocatable, dimension(:,:) :: row_fftbox_sum_batch
#endif
    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering density_on_dbl_grid'

    call timer_clock('density_on_dbl_grid',1)

    ! Obtain index for overlap matrix
    idx_len = sparse_index_length(overlap)
    allocate(overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('density_on_dbl_grid','overlap_idx',ierr)
    call sparse_generate_index(overlap_idx,overlap)

    ! Set batch size according to input file
    batch_size = pub_fftbox_batch_size

    ! agrecocmplx: optional variable for complex NGWFs
    loc_cmplx = ngwfs_on_grid_left%iscmplx

    ! Allocate workspace
    allocate(col_box_start(3, batch_size), stat=ierr)
    call utils_alloc_check('density_on_dbl_grid','col_box_start',ierr)
    allocate(col_start_in_box(3, batch_size), stat=ierr)
    call utils_alloc_check('density_on_dbl_grid','col_start_in_box',ierr)
    !allocate(row_fftbox_sum_batch(fftbox%total_ld1, &
    !     fftbox%total_ld2,fftbox%total_pt3, &
    !     pub_num_spins,batch_size),stat=ierr)
    !call utils_alloc_check('density_on_dbl_grid','row_fftbox_sum_batch',ierr)
    ! agrecocmplx: use appropriate routine to allocate row_fftbox_sum_batch
    ! agrecokpt: need to allocate partial sums row_fftbox_sum_batch_partial
    ! to take into account contributions from different k-points;
    ! alternatively, add extra dimension directly to row_fftbox_sum_batch
    allocate(row_fftbox_sum_batch(pub_num_spins,batch_size),stat=ierr)
    call utils_alloc_check('density_on_dbl_grid','row_fftbox_sum_batch',ierr)
    do is=1, pub_num_spins
        do batch=1, batch_size
            call data_fftbox_alloc(row_fftbox_sum_batch(is,batch), &
                fftbox%total_ld1, fftbox%total_ld2, fftbox%total_pt3, &
                iscmplx=loc_cmplx)
        end do
    end do

#ifdef GPU_PGI
    ! agrecokpt: TO-DO: check/implement k-points here as well
    call utils_assert(loc_cmplx == .false., &
         'GPU version of subroutine density_on_dbl_grid not tested yet for &
         & complex NGWFs.')

    allocate(row_fftbox_sum_batch_gpu(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3, &
         pub_num_spins,batch_size),stat=ierr)
    call utils_alloc_check('density_on_dbl_grid','row_fftbox_sum_batch_gpu',ierr)

    allocate(final_row_fftbox_dbl(fftbox%total_ld1_dbl,&
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl,pub_num_spins),stat=ierr)
    call utils_alloc_check('density_gpu_batch_interp_dep','final_row_fftbox',ierr)

    !allocate(col_fftbox(fftbox%total_ld1, &
    !     fftbox%total_ld2,fftbox%total_pt3),stat=ierr)
    !call utils_alloc_check('density_on_dbl_grid','col_fftbox',ierr)
    ! agrecocmplx
    call data_fftbox_alloc(col_fftbox, fftbox%total_ld1, &
             fftbox%total_ld2,fftbox%total_pt3, iscmplx=loc_cmplx)
    allocate(row_fftbox_dbl(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl,pub_num_spins,2),stat=ierr)
    call utils_alloc_check('density_on_dbl_grid','row_fftbox_dbl',ierr)
    allocate(buffer_dbl(fftbox%total_ld1_dbl, fftbox%total_ld2_dbl, &
         dbl_grid%max_group_slabs12, pub_num_spins),stat=ierr)
    call utils_alloc_check('density_on_dbl_grid','buffer_dbl',ierr)
#endif

    ! cks: zero before accumulation
    density_dbl = 0.0_DP

    ! agrecocmplx
    if (present(density_dbl_imag)) then
       density_dbl_imag = 0.0_DP
    end if

    ! cks: number of row-steps per row-block
    n_batches = ngwf_basis_right%max_on_proc / batch_size
    if (mod(ngwf_basis_right%max_on_proc, batch_size) > 0) n_batches = n_batches + 1

    ! kaw: Feedback to show switch between single and double precision and the
    !      current values of the related variables.
#ifdef GPU_SP_TEST
    if (GPU_SP) then
        if (pub_on_root) write(stdout,*) 'Using SP charge density routine', &
                         deltaE_gpu,GPU_SP_switch
    else
        if (pub_on_root) write(stdout,*) 'Using DP charge density routine', &
                         deltaE_gpu,GPU_SP_switch
    end if
#endif

    ! cks: loop over batches of NGWFs belonging to pub_my_proc_id
    local_start = 1
    do batch_count=1,n_batches

    if (pub_debug_on_root) write(stdout,'(a,2(i4,a))') 'DEBUG: Batch loop ', &
         batch_count,' of ',n_batches,' in density_on_fine_grid'

       ! cks: limits of my current batch in pub_my_proc_id-local NGWF counting
       ! cks: scheme
       local_end = min(local_start+batch_size-1,ngwf_basis_right%proc_num)

       ! cks: maximum size of current batch over all procs
       max_current_size = local_end - local_start + 1
       call comms_reduce('MAX', max_current_size)

       ! ndmh: find start positions of functions in this batch
       call function_ops_batch_col_start(col_box_start,col_start_in_box, &
            batch_size,local_start,local_end,fftbox,cell, &
            ngwf_basis_right)

       ! cks: zero before accumulation
       !row_fftbox_sum_batch = 0.0_DP
       ! agrecocmplx: initialise accordingly
       ! agrecokpt: need to loop over kpoints as well
       do is=1, pub_num_spins
           do batch=1, batch_size
               call data_set_to_zero(row_fftbox_sum_batch(is,batch))
           end do
       end do

       call timer_clock('density_batch_row_sums',1)

       ! ndmh: density-kernel block contributions to fftboxes of current batch
       ! agrecokpt: need to call this for each k-point, using
       ! row_fftbox_sum_batch_partial(ik), ngwfs_on_grid_left(ik),
       ! denskern(ik)
       call function_ops_sum_fftbox_batch(row_fftbox_sum_batch, &
            fftbox, cell, ngwfs_on_grid_left, ngwf_basis_left, &
            col_box_start, batch_size, local_start, local_end, &
            overlap_idx, idx_len, denskern, 1, 1, 1.0_DP)

       call timer_clock('density_batch_row_sums',2)
       call timer_clock('density_batch_interp_deposit',1)

       ! cks: interpolation of fftboxes of current batch, multiplication
       ! cks: and accumulation to charge-density 12-slabs of each proc
       ! agrecokpt: this needs to be called for each k-point, i.e.
       ! using row_fftbox_sum_batch_partial(ik) and ngwfs_on_grid_right(ik);
       ! final result will be the weighted sum for each ik component
#ifdef GPU_PGI
       !KAw: Introduced due to the issues arising from the combination of
       !derived types and pinned memory.
       do is=1, pub_num_spins
           do batch=1, batch_size
              row_fftbox_sum_batch_gpu(:,:,:,is,batch)=row_fftbox_sum_batch(is,batch)%d
           end do
       end do
       call density_gpu_batch_interp_dep(density_dbl, dbl_grid, &
            cell,fftbox,row_fftbox_sum_batch_gpu,ngwfs_on_grid_right, &
            ngwf_basis_right,col_box_start,col_start_in_box,local_start, &
            local_end,batch_size, max_current_size, &
            final_row_fftbox_dbl, col_fftbox, row_fftbox_dbl, buffer_dbl)
#else
       ! agrecocmplx
       if (present(density_dbl_imag)) then
          call density_batch_interp_deposit(density_dbl,dbl_grid, &
               cell,fftbox,row_fftbox_sum_batch,ngwfs_on_grid_right, &
               ngwf_basis_right,col_box_start,col_start_in_box,local_start,&
               local_end,batch_size,max_current_size,density_dbl_imag)
       else
          call density_batch_interp_deposit(density_dbl,dbl_grid, &
               cell,fftbox,row_fftbox_sum_batch,ngwfs_on_grid_right, &
               ngwf_basis_right,col_box_start,col_start_in_box,local_start,&
               local_end,batch_size,max_current_size)
       end if
#endif

       call timer_clock('density_batch_interp_deposit',2)

       ! Update first NGWF of next batch
       local_start = local_start + batch_size

    end do

    ! Deallocate workspace
    ! agrecocmplx: deallocate using appropriate routine
    ! agrecokpt: need to deallocate row_fftbox_sum_batch_partial
    ! as well, for each k-point
    do is=1, pub_num_spins
        do batch=1, batch_size
            call data_fftbox_dealloc(row_fftbox_sum_batch(is,batch))
        end do
    end do
    deallocate(row_fftbox_sum_batch,stat=ierr)
    call utils_dealloc_check('density_on_dbl_grid','row_fftbox_sum_batch',ierr)
    deallocate(col_start_in_box,stat=ierr)
    call utils_dealloc_check('density_on_dbl_grid','col_start_in_box',ierr)
    deallocate(col_box_start,stat=ierr)
    call utils_dealloc_check('density_on_dbl_grid','col_box_start',ierr)
    deallocate(overlap_idx,stat=ierr)
    call utils_dealloc_check('density_on_dbl_grid','overlap_idx',ierr)
#ifdef GPU_PGI
    deallocate(buffer_dbl,stat=ierr)
    call utils_dealloc_check('density_on_dbl_grid','row_fftbox_sum_batch',ierr)
    deallocate(row_fftbox_dbl,stat=ierr)
    call utils_dealloc_check('density_on_dbl_grid','row_fftbox_sum_batch',ierr)
    !deallocate(col_fftbox,stat=ierr)
    !call utils_dealloc_check('density_on_dbl_grid','row_fftbox_sum_batch',ierr)
    call data_fftbox_dealloc(col_fftbox)
    deallocate(final_row_fftbox_dbl,stat=ierr)
    call utils_dealloc_check('density_on_dbl_grid','row_fftbox_sum_batch',ierr)
#endif

    call timer_clock('density_on_dbl_grid',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving density_on_dbl_grid'

  end subroutine density_on_dbl_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_batch_interp_deposit(density_dbl, &
       dbl_grid, cell, fftbox, row_fftbox_sum_batch, ngwfs_on_grid, &
       ngwf_basis, col_box_start, col_start_in_box, local_start, local_end, &
       batch_size, max_current_size, density_dbl_imag)

    !==========================================================================!
    ! This subroutines interpolates each of the sums-in-fftboxes of the        !
    ! batch and its corresponding col NGWF. Then it deposits them in the       !
    ! correct position in the fine grid charge density in the simulation cell. !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !  density_dbl (input/output): Charge density in double-grid real space    !
    !    representation.                                                       !
    !  row_fftbox_sum_batch (input): Batch of fftboxes. Each fftbox            !
    !   contains the quantity phi_b*K^{ba} for the function phi_a.             !
    !  ngwfs_on_grid (input): Current NGWFs for pub_my_proc_id in ppd          !
    !    representation.                                                       !
    !  local_start (input): Number of the first NGWF of the current batch      !
    !    in the counting scheme of all NGWFs of pub_my_proc_id.                !
    !  local_end (input): Number of the last NGWF of the current batch         !
    !    in the counting scheme of all NGWFs of pub_my_proc_id.                !
    !  batch_size (input): Number of fftboxes (phi_a NGWFs) in each batch      !
    !--------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 15/1/2004 for the ONETEP code.       !
    ! Modified by Chris-Kriton Skylaris on 15/6/2004 so that it works          !
    ! with the data-parallel charge density.                                   !
    ! Modified by Nicholas Hine on 21/07/2009 to use function basis type.      !
    ! Modified by Nicholas Hine on 04/09/2009 to combine fourier_interpolate   !
    ! and multiplication - removes need for col_fftbox_dbl.                    !
    ! Modified by Nicholas Hine in November 2009 to accumulate the density for !
    ! each atom before depositing it to cell, for a parallel-efficiency gain.  !
    ! OpenMP parallelised by Karl Wilkinson and Nicholas Hine in May 2013.     !
    ! Modified by Andrea Greco on 08/05/2015 to allow use of complex NGWFs.    !
    !--------------------------------------------------------------------------!

    use basis, only: basis_copy_function_to_box, basis_location_func_wrt_cell
    use cell_grid, only: GRID_INFO, cell_grid_deposit_box
    use comms, only: comms_barrier, pub_my_proc_id, pub_rank_comm
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, &
         data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_interpolate, fourier_interpolate_product
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_num_spins
!$  use rundat, only: pub_threads_num_fftboxes, pub_threads_per_fftbox
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: local_start
    integer, intent(in) :: local_end
    integer, intent(in) :: batch_size
    integer, intent(in) :: max_current_size
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(inout) :: density_dbl(dbl_grid%ld1, &
         dbl_grid%ld2, dbl_grid%max_group_slabs12, pub_num_spins)
    type(FFTBOX_DATA), intent(in) :: row_fftbox_sum_batch(pub_num_spins,batch_size)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    integer, intent(in) :: col_start_in_box(3,batch_size)
    integer, intent(in) :: col_box_start(3,batch_size)
    ! agrecocmplx: accumulate imaginary part of density if required
    real(kind=DP), intent(inout), optional :: density_dbl_imag(dbl_grid%ld1, &
         dbl_grid%ld2, dbl_grid%max_group_slabs12, pub_num_spins)

    ! Local Variables
    integer :: is
    integer :: fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl
    integer :: col, local_col
    integer :: atom_of_col
    integer :: first_on_col_atom, last_on_col_atom
    integer :: batch_count
    integer :: ierr
    logical :: i_have_box
    real(kind=DP), allocatable, dimension(:,:,:) :: buffer_dbl
    real(kind=DP), allocatable, dimension(:,:,:,:) :: row_fftbox_dbl
    ! agreco: container to accumulate imaginary part of complex density
    ! used to check the computed density is real as expected
    real(kind=DP), allocatable, dimension(:,:,:,:) :: row_fftbox_dbl_imag
    type(FFTBOX_DATA) :: col_fftbox
    complex(kind=DP), allocatable, dimension(:,:,:) :: coarse_work,fine_work
    ! agrecocmplx: extra workspaces for complex NGWFs
    complex(kind=DP), allocatable, dimension(:,:,:) :: fine_work_2
    integer :: thread_id,nthreads
    integer :: thread_local_start, thread_local_end
    ! agrecocmplx: extra variable for complex quantities
    logical :: loc_cmplx

!$  integer, external :: omp_get_thread_num, omp_get_num_threads

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering density_batch_interp_deposit'

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid%iscmplx

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(local_col, col, atom_of_col, i_have_box, &
!$OMP      first_on_col_atom,  last_on_col_atom, batch_count, &
!$OMP      is, fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl,&
!$OMP      row_fftbox_dbl, row_fftbox_dbl_imag, buffer_dbl,ierr, &
!$OMP      thread_local_start, thread_local_end, thread_id, nthreads, &
!$OMP      coarse_work, fine_work, fine_work_2) &
!$OMP FIRSTPRIVATE(col_fftbox) &
!$OMP SHARED(row_fftbox_sum_batch, density_dbl, density_dbl_imag, &
!$OMP      ngwfs_on_grid, ngwf_basis, dbl_grid, col_box_start, &
!$OMP      col_start_in_box, fftbox, cell, pub_num_spins, local_start, &
!$OMP      local_end, max_current_size, pub_my_proc_id, pub_rank_comm, &
!$OMP      pub_threads_num_fftboxes, pub_threads_per_fftbox, loc_cmplx)

    ! Allocate Workspace
    call data_fftbox_alloc(col_fftbox, fftbox%total_ld1, &
             fftbox%total_ld2,fftbox%total_pt3, iscmplx=loc_cmplx)
    allocate(row_fftbox_dbl(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl,2),stat=ierr)
    call utils_alloc_check('density_batch_interp_deposit','row_fftbox_dbl',ierr)
    allocate(buffer_dbl(fftbox%total_ld1_dbl, fftbox%total_ld2_dbl, &
         dbl_grid%max_group_slabs12),stat=ierr)
    call utils_alloc_check('density_batch_interp_deposit','buffer_dbl',ierr)
    allocate(coarse_work(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('density_batch_interp_deposit','coarse_work',ierr)
    allocate(fine_work(fftbox%total_ld1*2, &
         fftbox%total_ld2*2,fftbox%total_pt3*2),stat=ierr)
    call utils_alloc_check('density_batch_interp_deposit','fine_work',ierr)
    allocate(fine_work_2(fftbox%total_ld1*2, &
         fftbox%total_ld2*2,fftbox%total_pt3*2),stat=ierr)
    call utils_alloc_check('density_batch_interp_deposit','fine_work_2',ierr)
    ! agreco: extra arrays only if we need to check density is real when
    ! using complex NGWFs
    if (present(density_dbl_imag)) then
       allocate(row_fftbox_dbl_imag(fftbox%total_ld1_dbl, &
          fftbox%total_ld2_dbl,fftbox%total_pt3_dbl,2),stat=ierr)
       call utils_alloc_check('density_batch_interp_deposit', &
            'row_fftbox_dbl_imag', ierr)
       row_fftbox_dbl_imag = 0.0_DP
    end if

    ! pdh: zero workspace
    row_fftbox_dbl = 0.0_DP

    ! ndmh: Set iteration range of each thread manually here - could use OMP do
    ! ndmh: but then range of iterations of each loop is hard to predict
    ! ndmh: between architectures. This way allows each thread to decide whether
    ! ndmh: to deposit its density box at the end of each iteration.
    thread_id = 0
    nthreads = 1
!$  nthreads = omp_get_num_threads()
!$  thread_id = omp_get_thread_num()
    thread_local_start = (thread_id * max_current_size) / nthreads &
         + local_start
    thread_local_end = ((thread_id + 1) * max_current_size) / nthreads &
         + local_start - 1

    do is=1,pub_num_spins
       ! cks: loop over the members of the largest current batch (from all procs)
       ! cks: and interpolate each pair of fftboxes of each proc, multiply them
       ! cks: together and deposit them accordingly to the 12-slabs of all procs
       do local_col=thread_local_start,thread_local_end
          batch_count = local_col - local_start + 1

          call timer_clock('density_fftbox_interpolate_multiply', 1)

          ! cks: interpolate and multiply my col NGWF if pub_my_proc_id
          ! cks: posseses an fftbox for this local_col value
          if ( local_col .le. local_end) then

             ! ndmh: find information about this column NGWF
             col = local_col + ngwf_basis%first_on_proc(pub_my_proc_id) - 1
             atom_of_col = ngwf_basis%atom_of_func(col)
             first_on_col_atom = ngwf_basis%first_on_atom(atom_of_col)
             last_on_col_atom = first_on_col_atom + &
                  ngwf_basis%num_on_atom(atom_of_col) - 1

             ! ndmh: copy column NGWF to column fftbox
             !call basis_copy_function_to_box(col_fftbox, &
             !     fftbox%total_pt1,fftbox%total_pt2,fftbox%total_pt3, &
             !     col_start_in_box(1,batch_count), &
             !     col_start_in_box(2,batch_count), &
             !     col_start_in_box(3,batch_count), &
             !     ngwf_basis%tight_boxes(local_col), ngwfs_on_grid, &
             !     ngwf_basis%spheres(local_col), cell, fftbox)
             ! agrecocmplx: call modified routine in basis_new_mod
             call basis_copy_function_to_box(col_fftbox, &
                  col_start_in_box(1,batch_count), &
                  col_start_in_box(2,batch_count), &
                  col_start_in_box(3,batch_count), &
                  ngwf_basis%tight_boxes(local_col), ngwfs_on_grid, &
                  ngwf_basis%spheres(local_col), cell, fftbox)

             ! ndmh: interpolate col and sum of rows together with c2c ffts
             ! ndmh: then multiply them together to find contribution to the
             ! ndmh: charge density
             !call fourier_interpolate_product(coarse_work, fine_work, &
             !     row_fftbox_sum_batch(:,:,:,is,batch_count), col_fftbox, &
             !     row_fftbox_dbl(:,:,:,1))
             ! agrecocmplx: when using real NGWFs, use fourier_interpolate_product
             ! when using complex NGWFs, interpolate each complex object independently,
             ! then take the real part of the product of the results;
             if (.not.loc_cmplx) then
                 call fourier_interpolate_product(coarse_work, fine_work, &
                      row_fftbox_sum_batch(is,batch_count)%d, col_fftbox%d, &
                      row_fftbox_dbl(:,:,:,1))
             else
                 call fourier_interpolate(coarse_work, fine_work, &
                      row_fftbox_sum_batch(is,batch_count)%z)
                 call fourier_interpolate(coarse_work, fine_work_2, &
                      col_fftbox%z)
                 ! agrecocmplx: need to take conjugate of LHS ngwfs_on_grid
                 ! check this since it's critical
                 ! agreco_critical
                 fine_work(:,:,:) = conjg(fine_work(:,:,:))*fine_work_2(:,:,:)

                 row_fftbox_dbl(:,:,:,1) = real(fine_work(:,:,:))

                 ! agrecocmplx
                 if (present(density_dbl_imag)) then
                    row_fftbox_dbl_imag(:,:,:,1) = aimag(fine_work(:,:,:))
                 end if

             end if

             ! ndmh: if this is a new atom or the start of a new batch, reset the
             ! ndmh: accumulated density for this atom to the density for this NGWF
             if ((col == first_on_col_atom).or.(local_col==thread_local_start)) then
                row_fftbox_dbl(:,:,:,2) = row_fftbox_dbl(:,:,:,1)
                ! agrecocmplx
                if (present(density_dbl_imag)) then
                   row_fftbox_dbl_imag(:,:,:,2) = row_fftbox_dbl_imag(:,:,:,1)
                end if
             else ! ndmh: add the density for this NGWF to the sum for this atom
                row_fftbox_dbl(:,:,:,2) = row_fftbox_dbl(:,:,:,2) + &
                     row_fftbox_dbl(:,:,:,1)
                ! agrecocmplx
                if (present(density_dbl_imag)) then
                   row_fftbox_dbl_imag(:,:,:,2) = row_fftbox_dbl_imag(:,:,:,2) + &
                        row_fftbox_dbl_imag(:,:,:,1)
                end if
             end if

             ! ndmh: if this is the last NGWF on this atom or the end of the batch,
             ! ndmh: we will need to deposit the density to the whole-cell array
             if ((col == last_on_col_atom).or.(local_col==thread_local_end)) then
                i_have_box = .true.
             else
                i_have_box = .false.
             end if

          else
             i_have_box =.false.
          end if

          ! ndmh: synchronise procs in this rank so that load-balancing time is
          ! ndmh: reported in density_fftbox_interpolate_multiply rather than as
          ! ndmh: the wait for the alltoall in cell_grid_deposit_box.
!$        if (.false.) then
          call comms_barrier(pub_rank_comm)
!$        end if
          call timer_clock('density_fftbox_interpolate_multiply', 2)
          call timer_clock('density_fftbox_deposit_to_cell', 1)

          ! cks: get into the depositing fftboxes to 12-slabs of all procs
          ! cks: regardless of if pub_my_proc_id has an fftbox to deposit
          ! ndmh: box_to_cell routine moved to basis_mod
          fftbox_start1_dbl = 2*col_box_start(1,batch_count) - 1
          fftbox_start2_dbl = 2*col_box_start(2,batch_count) - 1
          fftbox_start3_dbl = 2*col_box_start(3,batch_count) - 1
!$OMP CRITICAL
          call cell_grid_deposit_box(density_dbl(:,:,:,is), &
               row_fftbox_dbl(:,:,:,2), buffer_dbl, dbl_grid, &
               fftbox%total_pt1_dbl, fftbox%total_pt2_dbl, &
               fftbox%total_pt3_dbl, fftbox%total_ld1_dbl, &
               fftbox%total_ld2_dbl, fftbox_start1_dbl, &
               fftbox_start2_dbl, fftbox_start3_dbl, i_have_box, .true.)

          if (present(density_dbl_imag)) then
             call cell_grid_deposit_box(density_dbl_imag(:,:,:,is), &
                  row_fftbox_dbl_imag(:,:,:,2), buffer_dbl, dbl_grid, &
                  fftbox%total_pt1_dbl, fftbox%total_pt2_dbl, &
                  fftbox%total_pt3_dbl, fftbox%total_ld1_dbl, &
                  fftbox%total_ld2_dbl, fftbox_start1_dbl, &
                  fftbox_start2_dbl, fftbox_start3_dbl, i_have_box, .true.)
          end if
!$OMP END CRITICAL

          call timer_clock('density_fftbox_deposit_to_cell', 2)

       end do

!$OMP BARRIER

    end do

    deallocate(fine_work,stat=ierr)
    call utils_dealloc_check('density_batch_interp_deposit','fine_work',ierr)
    ! agrecocmplx
    deallocate(fine_work_2,stat=ierr)
    call utils_dealloc_check('density_batch_interp_deposit','fine_work_2',ierr)
    deallocate(coarse_work,stat=ierr)
    call utils_dealloc_check('density_batch_interp_deposit','coarse_work',ierr)

    deallocate(buffer_dbl,stat=ierr)
    call utils_dealloc_check('density_batch_interp_deposit','buffer_dbl',ierr)
    deallocate(row_fftbox_dbl,stat=ierr)
    call utils_dealloc_check('density_batch_interp_deposit','row_fftbox_dbl',ierr)
    call data_fftbox_dealloc(col_fftbox)
!$OMP END PARALLEL

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving density_batch_interp_deposit'

  end subroutine density_batch_interp_deposit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_plot_slice(coord,xvec,yvec,ovec,nx,ny,data,cell,n1,n2,n3, &
       fname,title,xlabel,ylabel)

    !=================================================================!
    ! This routine plots a slice of a 3D array in MTV format          !
    !-----------------------------------------------------------------!
    ! Written by Peter D. Haynes in summer 2004.                      !
    !=================================================================!


    use comms, only: pub_on_root
    use constants, only: DP, stdout, PI
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, &
         utils_assert

    implicit none

    ! Arguments
    character, intent(in) :: coord         ! (C)artesian or (F)ractional
    real(kind=DP), intent(in) :: xvec(:)   ! Vector giving x-axis
    real(kind=DP), intent(in) :: yvec(:)   ! Vector giving y-axis
    real(kind=DP), intent(in) :: ovec(:)   ! Vector giving origin
    integer, intent(in) :: nx              ! Number of points along x
    integer, intent(in) :: ny              ! Number of points along y
    real(kind=DP), intent(in) :: data(:,:,:)  ! The data to plot
    integer, intent(in) :: n1,n2,n3        ! Size of data
    type(CELL_INFO), intent(in) :: cell
    character(len=*), intent(in) :: fname  ! Filename
    character(len=*), optional, intent(in) :: title
    character(len=*), optional, intent(in) :: xlabel
    character(len=*), optional, intent(in) :: ylabel

    ! Local variables
    integer, parameter :: iounit = 20   ! Unit for output file
    integer, parameter :: finesse = 4   ! Accuracy of spline interpolation
    integer :: ierr                     ! Error flag
    integer :: ix,iy                    ! Plot point counters
    integer :: i1,i2,i3                 ! Interpolation counters
    integer :: j1,j2,j3                 ! Interpolation counters
    integer :: i,j,k                    ! Loop counters
    integer :: n(3)                     ! Copy of n1,n2,n3
    integer :: m1,m2,m3                 ! Number of points to interpolate from
    real(kind=DP), parameter :: cutoff = 2.0_DP     ! Interpolation cutoff
    real(kind=DP) :: x,h,delta                      ! Grid spacing/position
    real(kind=DP) :: normfac                        ! PSinc normalisation
    real(kind=DP) :: xaxrel(3),yaxrel(3),orgrel(3)  ! Plot frame definition
    real(kind=DP) :: ptrel(3)                       ! Plot point
    real(kind=DP) :: psinc23,psinc3                 ! Interpolation variables
    real(kind=DP), allocatable :: psinc(:,:)        ! Value of PSinc for spline
    real(kind=DP), allocatable :: d2psinc(:,:)      ! Value of 2nd deriv
    real(kind=DP), allocatable :: psincpt(:,:)      ! Value of PSinc at point
    real(kind=DP), allocatable :: line(:)           ! Line of data for output

    ! Check arguments
    if (coord /= 'C' .and. coord /= 'c' .and. coord /= 'F' .and. &
         coord /= 'f') then
       call utils_abort('Error in density_plot_slice: &
            &invalid coordinate specifier: '//trim(coord))
    end if
    if (size(xvec) < 3) then
       call utils_abort('Error in density_plot_slice: &
            &x-axis requires a 3-component vector')
    end if
    if (size(yvec) < 3) then
       call utils_abort('Error in density_plot_slice: &
            &y-axis requires a 3-component vector')
    end if
    if (size(ovec) < 3) then
       call utils_abort('Error in density_plot_slice: &
            &origin requires a 3-component vector')
    end if
    call utils_assert(nx>=1, 'Error in density_plot_slice: # of x-points < 1')
    call utils_assert(ny>=1, 'Error in density_plot_slice: # of y-points < 1')

    call utils_assert(size(data,1) >= n1, 'Error in density_plot_slice: &
            &data array mismatch in index 1')
    call utils_assert(size(data,2) >= n2, 'Error in density_plot_slice: &
            &data array mismatch in index 2')
    call utils_assert(size(data,3) >= n3, 'Error in density_plot_slice: &
            &data array mismatch in index 3')

    if (pub_on_root) write(stdout,'(2a)') &
         '>>>>>>>>>>>> Writing slice to file: ', trim(fname)

    ! Allocate workspace
    allocate(line(nx),stat=ierr)
    call utils_alloc_check('density_plot_slice','line',ierr)
    allocate(psinc(finesse*max(n1,n2,n3),3),stat=ierr)
    call utils_alloc_check('density_plot_slice','psinc',ierr)
    allocate(d2psinc(finesse*max(n1,n2,n3),3),stat=ierr)
    call utils_alloc_check('density_plot_slice','d2psinc',ierr)
    allocate(psincpt(max(n1,n2,n3),3),stat=ierr)
    call utils_alloc_check('density_plot_slice','psincpt',ierr)

    ! Convert to fractional coordinates if necessary
    if (coord == 'C' .or. coord == 'c') then
       xaxrel(1) = (xvec(1) * cell%b1%x + xvec(2) * cell%b1%y + &
            xvec(3) * cell%b1%z) * 0.5_DP / pi
       xaxrel(2) = (xvec(1) * cell%b2%x + xvec(2) * cell%b2%y + &
            xvec(3) * cell%b2%z) * 0.5_DP / pi
       xaxrel(3) = (xvec(1) * cell%b3%x + xvec(2) * cell%b3%y + &
            xvec(3) * cell%b3%z) * 0.5_DP / pi
       yaxrel(1) = (yvec(1) * cell%b1%x + yvec(2) * cell%b1%y + &
            yvec(3) * cell%b1%z) * 0.5_DP / pi
       yaxrel(2) = (yvec(1) * cell%b2%x + yvec(2) * cell%b2%y + &
            yvec(3) * cell%b2%z) * 0.5_DP / pi
       yaxrel(3) = (yvec(1) * cell%b3%x + yvec(2) * cell%b3%y + &
            yvec(3) * cell%b3%z) * 0.5_DP / pi
       orgrel(1) = (ovec(1) * cell%b1%x + ovec(2) * cell%b1%y + &
            ovec(3) * cell%b1%z) * 0.5_DP / pi
       orgrel(2) = (ovec(1) * cell%b2%x + ovec(2) * cell%b2%y + &
            ovec(3) * cell%b2%z) * 0.5_DP / pi
       orgrel(3) = (ovec(1) * cell%b3%x + ovec(2) * cell%b3%y + &
            ovec(3) * cell%b3%z) * 0.5_DP / pi
    else
       xaxrel = xvec(1:3)
       yaxrel = yvec(1:3)
       orgrel = ovec(1:3)
    end if

    ! Write header to file
    if (pub_on_root) then
       open(unit=iounit,file=trim(fname))
       write(iounit,'(a)') '$DATA=CONTOUR'
       write(iounit,'(a,i4)') '% XMIN=0  XMAX=1  NX=',nx
       write(iounit,'(a,i4)') '% YMIN=0  YMAX=1  NY=',ny
       if (present(title)) write(iounit,'(3a)') '% TOPLABEL="',trim(title),'"'
       write(iounit,'(a)') '% EQUALSCALE'
       if (present(xlabel)) write(iounit,'(3a)') '% XLABEL="',trim(xlabel),'"'
       if (present(ylabel)) write(iounit,'(3a)') '% YLABEL="',trim(ylabel),'"'
       write(iounit,'(a)') '% CONTSTYLE=2'
       write(iounit,'(a)') '% NSTEPS=10'
    end if

    ! Set up spline-fitted PSinc functions
    n(1) = n1 ; n(2) = n2 ; n(3) = n3
    do k=1,3
       h = 1.0_DP / real(finesse*n(k),kind=DP)
       do i=1,finesse*n(k)
          x = (i-1) * h
          psinc(i,k) = 0.5_DP
          do j=1,(n(k)-1)/2   ! omit Nyquist for even grid to keep psinc real
             psinc(i,k) = psinc(i,k) + cos(2.0_DP*pi*j*x)
          end do
       end do
       normfac = 1.0_DP / (0.5_DP + real((n(k)-1)/2,kind=DP))
       psinc(1:finesse*n(k),k) = psinc(1:finesse*n(k),k) * normfac
       call internal_spline_uni_fit(finesse*n(k),h,psinc(:,k),d2psinc(:,k))
    end do

    ! Calculate set of grid points to use in interpolation
    m1 = min(int(2.0_DP * cutoff / cell%d1),(n1-1)/2)
    m2 = min(int(2.0_DP * cutoff / cell%d2),(n2-1)/2)
    m3 = min(int(2.0_DP * cutoff / cell%d3),(n3-1)/2)

    ! Loop over plotting points
    do iy=1,ny
       do ix=1,nx

          ! Get position of point in fractional coordinates
          ptrel = orgrel + (ix-1)*xaxrel/nx + (iy-1)*yaxrel/ny

          ! Interpolate PSinc's for this point
          h = 1.0_DP / real(finesse*n1,kind=DP)
          delta = 1.0_DP / real(n1,kind=DP)
          do j1=-m1,m1
             i1 = modulo(int(ptrel(1)*n1) + j1,n1) + 1
             x = modulo((i1-1) * delta - ptrel(1),1.0_DP)
             psincpt(i1,1) = internal_spline_uni_eval(finesse*n1,h,psinc(:,1), &
                  d2psinc(:,1),x)
          end do
          h = 1.0_DP / real(finesse*n2,kind=DP)
          delta = 1.0_DP / real(n2,kind=DP)
          do j2=-m2,m2
             i2 = modulo(int(ptrel(2)*n2) + j2,n2) + 1
             x = modulo((i2-1) * delta - ptrel(2),1.0_DP)
             psincpt(i2,2) = internal_spline_uni_eval(finesse*n2,h,psinc(:,2), &
                  d2psinc(:,2),x)
          end do
          h = 1.0_DP / real(finesse*n3,kind=DP)
          delta = 1.0_DP / real(n3,kind=DP)
          do j3=-m3,m3
             i3 = modulo(int(ptrel(3)*n3) + j3,n3) + 1
             x = modulo((i3-1) * delta - ptrel(3),1.0_DP)
             psincpt(i3,3) = internal_spline_uni_eval(finesse*n3,h,psinc(:,3), &
                  d2psinc(:,3),x)
          end do

          ! Evaluate function at this point
          line(ix) = 0.0_DP
          do j3=-m3,m3
             i3 = modulo(int(ptrel(3)*n3) + j3,n3) + 1
             psinc3 = psincpt(i3,3)
             do j2=-m2,m2
                i2 = modulo(int(ptrel(2)*n2) + j2,n2) + 1
                psinc23 = psincpt(i2,2) * psinc3
                do j1=-m1,m1
                   i1 = modulo(int(ptrel(1)*n1) + j1,n1) + 1
                   line(ix) = line(ix) + data(i1,i2,i3) * psincpt(i1,1) * &
                        psinc23
                end do
             end do
          end do

       end do

       ! Write this data to file
       if (pub_on_root) write(iounit,'(8(f9.4,1x))') line

    end do

    ! Write footer to file
    if (pub_on_root) then
       write(iounit,'(a)') '$ END'
       close(iounit)
    end if

    ! Deallocate workspace
    deallocate(psincpt,stat=ierr)
    call utils_dealloc_check('density_plot_slice','psincpt',ierr)
    deallocate(d2psinc,stat=ierr)
    call utils_dealloc_check('density_plot_slice','d2psinc',ierr)
    deallocate(psinc,stat=ierr)
    call utils_dealloc_check('density_plot_slice','psinc',ierr)
    deallocate(line,stat=ierr)
    call utils_dealloc_check('density_plot_slice','line',ierr)

  contains

    !==========================================================================
    ! Fit a natural cubic spline to data tabulated on a uniform grid
    !==========================================================================

    subroutine internal_spline_uni_fit(n,h,y,ypp)

      use utils, only: utils_assert

      implicit none
      integer, intent(in) :: n
      real(kind=DP), intent(in) :: h,y(n)
      real(kind=DP), intent(out) :: ypp(n)
      real(kind=DP), allocatable :: d(:),e(:)
      integer :: i, ierr, info

      ! LAPACK subroutine
      external :: dptsv

      allocate(d(n),stat=ierr)
      call utils_alloc_check('internal_spline_uni_fit (density_mod.F90)',&
           'd',ierr)
      allocate(e(n-1),stat=ierr)
      call utils_alloc_check('internal_spline_uni_fit (density_mod.F90)',&
           'e',ierr)

      d(:) = 4.0_DP*h*h
      e(:) = h*h
      ypp(1) = y(2)-2.0_DP*y(1)+y(n)
      do i=2,n-1
         ypp(i) = y(i+1)-2.0_DP*y(i)+y(i-1)
      end do
      ypp(n) = y(1)-2.0_DP*y(n)+y(n-1)
      ypp(:) = 6.0_DP * ypp(:)
      call dptsv(n,1,d,e,ypp,n,info)
      call utils_assert(info == 0, &
           'Error in internal_spline_uni_fit (density_mod.F90): ', info)
      deallocate(e,stat=ierr)
      call utils_dealloc_check('internal_spline_uni_fit (density_mod.F90)',&
           'e',ierr)
      deallocate(d,stat=ierr)
      call utils_dealloc_check('internal_spline_uni_fit (density_mod.F90)',&
           'd',ierr)

    end subroutine internal_spline_uni_fit


    !==========================================================================
    ! Evaluate a natural cubic spline on a uniform grid assumed to start at x=0
    !==========================================================================

    real(kind=DP) function internal_spline_uni_eval(n,h,y,ypp,x)

      implicit none
      integer, intent(in) :: n
      real(kind=DP), intent(in) :: h,x,y(n),ypp(n)
      integer :: ix
      real(kind=DP), parameter :: sixth = 1.0d0 / 6.0d0
      real(kind=DP) :: a,b,pos

      pos = x / h
      ix = int(pos)
      b = pos - real(ix,kind=DP)
      a = 1.0_DP - b
      internal_spline_uni_eval = a*y(ix+1)+b*y(ix+2) + &
           ((a*a*a-a)*ypp(ix+1)+(b*b*b-b)*ypp(ix+2))*h*h*sixth

    end function internal_spline_uni_eval

  end subroutine density_plot_slice

  subroutine castep_density_dump(den,filename,cell,grid)
    !=========================================================================!
    ! Output charge density used by multigrid solver to a binary file, in     !
    ! CASTEP-compatible format.                                               !
    ! All file IO is done by the root proc and the full simulation cell       !
    ! density is temporarily stored on the root proc for output.              !
    !                                                                         !
    ! Notes:                                                                  !
    ! CASTEP stores the electron density as a total charge. Therefore spin    !
    ! components should be summed before passing a total charge density, den, !
    ! as an argument to this routine. The charge density is multiplied by the !
    ! cell volume inside this routine.                                        !
    !                                                                         !
    ! CASTEP appears to have standardized on big endian byte ordering. Since  !
    ! this routine outputs unformatted binary data to a file, it is important !
    ! that the endianness of the data is in agreement with the endianness     !
    ! expected by CASTEP (typically compiled ONETEP with -fconvert=big-endian !
    ! compiler option).                                                       !
    !-------------------------------------------------------------------------!
    ! J. C. Womack, 12/2016                                                   !
    !=========================================================================!
    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root, pub_root_proc_id, pub_total_num_procs, &
         pub_my_proc_id, comms_barrier, &
         comms_recv, comms_reduce, comms_send, comms_wait
    use constants, only: stdout
    use geometry, only: OPERATOR(.dot.), operator(.cross.)
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_assert, &
         utils_dealloc_check, utils_unit
    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: den(:,:,:)  !< charge density for this MPI rank
    character(len=*), intent(in) :: filename !< filename to which density will be dumped
    type(CELL_INFO), intent(in) :: cell
    type(GRID_INFO), intent(in) :: grid

    ! Parameters
    character(len=*), parameter :: myself = "castep_density_dump"

    ! Local variables
    real(kind=DP), allocatable :: den_to_dump(:,:,:) !< Charge density for all procs
    real(kind=DP), allocatable :: den_recv(:,:,:)    !< Temporary array to recv data

    integer :: dump_unit
    integer :: ierr
    integer :: irank
    integer :: num_slabs
    integer :: nx, ny
    real(kind=DP) :: norm
    real(kind=DP) :: charge
    real(kind=DP) :: volume
    integer :: my_request_handle

    ! Calculate volume of simulation cell (so density can be converted to charge)
    volume = abs((cell%a1.cross.cell%a2).dot.cell%a3)

    if (pub_on_root) then
       ! Allocate space for full charge density
       allocate(den_to_dump(grid%n1,grid%n2,grid%n3), stat=ierr)
       call utils_alloc_check(myself,'den_to_dump',ierr)
       ! Allocate space to recv data into
       ! Note that the charge density stores on each notes has the dimensions
       ! (ld1,ld2,max_slabs12) and ld1, ld2 do not in general equal n1, n2
       allocate(den_recv(grid%ld1,grid%ld2,grid%max_slabs12), stat=ierr)
       call utils_alloc_check(myself,'den_recv',ierr)
       call utils_assert(pub_root_proc_id==0,"Error in "//myself//": &
            &pub_root_proc_id==0 assumed but this is not the case. Abort.")
       ! Place root MPI rank data into den_to_dump
       irank = 0
       num_slabs = grid%last_slab12(irank) - grid%first_slab12(irank) + 1
       den_to_dump(1:grid%n1,1:grid%n2,grid%first_slab12(irank):grid%last_slab12(irank)) = &
            den(1:grid%n1,1:grid%n2,1:num_slabs)
       ! Receive data from other MPI ranks and copy into den_to_dump
       do irank = 1, pub_total_num_procs-1
          num_slabs = grid%last_slab12(irank) - grid%first_slab12(irank) + 1
          call comms_recv(irank,den_recv,length=grid%ld1*grid%ld2*num_slabs)
          !call comms_recv(irank,den_recv)
          den_to_dump(1:grid%n1,1:grid%n2,grid%first_slab12(irank):grid%last_slab12(irank)) = &
               den_recv(1:grid%n1,1:grid%n2,1:num_slabs)
       end do
       deallocate(den_recv, stat=ierr)
       call utils_dealloc_check(myself,'den_recv',ierr)

    else
       call comms_send(pub_root_proc_id, den, &
            length=grid%ld1*grid%ld2*grid%num_my_slabs12, &
            return_handle=my_request_handle,&
            add_to_stack = .false.)
       ! Wait for all sending and receiving to complete
       call comms_wait(my_request_handle)
    end if

    ! Barrier before proceeding
    call comms_barrier

    ! Check norm of data gathered to root proc against norm obtained
    ! by a comms_reduce
    charge = sum( den(1:grid%n1,1:grid%n2,1:grid%num_my_slabs12 ) ) * grid%weight
    norm   = sum( den(1:grid%n1,1:grid%n2,1:grid%num_my_slabs12 )**2 )
    call comms_reduce('SUM',norm,root=pub_root_proc_id)
    call comms_reduce('SUM',charge,root=pub_root_proc_id)
    if (pub_on_root) then
       ! Check total norm of collected data
       write(stdout,'(a)') "<< CASTEP density dump started >>"
       write(stdout,'(a,es20.10)') "Volume of cell / bohr**3 =     ", volume
       norm = sqrt( norm )
       write(stdout,'(a,es20.10)') "Norm by comms_reduce =         ", norm
       write(stdout,'(a,es20.10)') "Total charge by comms_reduce = ", charge
       charge = sum( den_to_dump(1:grid%n1,1:grid%n2,1:grid%n3 ) ) * grid%weight
       norm = sqrt( sum( den_to_dump(1:grid%n1,1:grid%n2,1:grid%n3)**2 ) )
       write(stdout,'(a,es20.10)') "Norm by gathering    =         ", norm
       write(stdout,'(a,es20.10)') "Total charge by gathering    = ", charge
    end if

    if (pub_on_root) then
       ! CASTEP expects the charge for each grid point, as if it was uniformly
       ! spread over the volume of the simulation cell. Therefore, multiply
       ! electron density at each grid point by volume of cell
       den_to_dump = den_to_dump * volume


       ! Get free IO unit
       dump_unit = utils_unit()
       ! Open binary file for output
       open(unit=dump_unit,file=trim(filename),action="write", form="unformatted", &
            status="replace", access="sequential" )

       ! Dump density data to file in format readable by CASTEP's density_read_serial
       ! routine


       ! Data has format: nx, ny, charge_col(1:nz)
       ! where nx, ny are the x and y coordinates of a z-column of density
       ! and charge_col is an array of data containing the charge for all
       ! grid points (nx,ny,1:nz).

       ! NOTE: CASTEP expects BIG ENDIAN unformatted data, so ONETEP must be compiled
       ! to output in big endian format, rather than little endian (little endian is
       ! the standard on Intel x86_64 platforms.
       ! NOTE: CASTEP expects complex input, so even though the imaginary part is always
       ! 0.0_DP we write the density data in complex form.
       do nx = 1, grid%n1
          do ny = 1, grid%n2
             ! Dump density in complex form, since CASTEP expects a complex variable
             write(dump_unit) nx, ny, cmplx( den_to_dump(nx,ny,1:grid%n3), 0.0_DP, kind=DP)
          end do
       end do
       close(unit=dump_unit)

       write(stdout,'(a,a)') "Saved charge density to ", filename
       write(stdout,'(a)') "<< CASTEP density dump completed >>"
    end if

    if (pub_on_root) then
    ! Deallocate temporary space for full charge density
       deallocate(den_to_dump, stat=ierr)
       call utils_dealloc_check(myself,'den_to_dump',ierr)
    end if

  end subroutine castep_density_dump


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef GPU_PGI
  subroutine density_gpu_batch_interp_dep(density_dbl, &
       dbl_grid, cell, fftbox, row_fftbox_sum_batch, ngwfs_on_grid, &
       ngwf_basis, col_box_start, col_start_in_box, local_start, &
       local_end, batch_size, max_current_size, &
       final_row_fftbox_dbl, col_fftbox, row_fftbox_dbl, buffer_dbl )
    !==========================================================================!
    ! This subroutines interpolates each of the sums-in-fftboxes of the        !
    ! batch and its corresponding col NGWF. Then it deposits them in the       !
    ! correct position in the fine grid charge density in the simulation cell. !
    !                                                                          !
    ! The routine has been modified for GPU execution by changing the order in !
    ! which the stages of the process occur. This allows asynchronous data     !
    ! transfer and density calculation.                                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !  density_dbl (input/output): Charge density in double-grid real space    !
    !    representation.                                                       !
    !  row_fftbox_sum_batch (input): Batch of fftboxes. Each fftbox            !
    !   contains the quantity phi_b*K^{ba} for the function phi_a.             !
    !  ngwfs_on_grid (input): Current NGWFs for pub_my_proc_id in ppd          !
    !    representation.                                                       !
    !  local_start (input): Number of the first NGWF of the current batch      !
    !    in the counting scheme of all NGWFs of pub_my_proc_id.                !
    !  local_end (input): Number of the last NGWF of the current batch         !
    !    in the counting scheme of all NGWFs of pub_my_proc_id.                !
    !  batch_size (input): Number of fftboxes (phi_a NGWFs) in each batch      !
    !--------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 15/1/2004 for the ONETEP code.       !
    ! Modified by Chris-Kriton Skylaris on 15/6/2004 so that it works          !
    ! with the data-parallel charge density.                                   !
    ! Modified by Nicholas Hine on 21/07/2009 to use function basis type.      !
    ! Modified by Nicholas Hine on 04/09/2009 to combine fourier_interpolate   !
    ! and multiplication - removes need for col_fftbox_dbl.                    !
    ! Modified by Nicholas Hine in November 2009 to accumulate the density for !
    ! each atom before depositing it to cell, for a parallel-efficiency gain.  !
    ! Modified by Karl Wilkinson in August 2012 for GPU implementation.        !
    !--------------------------------------------------------------------------!

    !kaw merge change use basis, only: basis_copy_function_to_box, basis_location_func_wrt_cell
    use basis, only: basis_copy_function_to_box, basis_location_func_wrt_cell
    use cell_grid, only: GRID_INFO, cell_grid_deposit_box
    use comms, only: comms_barrier, pub_my_proc_id, &
         pub_on_root, pub_rank_comm
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, &
         data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
#ifdef GPU_SP_TEST
    use fourier, only: fourier_gpu_interp_prod_SP, fourier_gpu_interp_prod
#else
    use fourier, only: fourier_gpu_interp_prod
#endif
    use fourier_gpu_wrapper_mod
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_num_spins, pub_cmplx_ngwfs
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort,
         utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: local_start
    integer, intent(in) :: local_end
    integer, intent(in) :: batch_size
    integer, intent(in) :: max_current_size
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(inout) :: density_dbl(dbl_grid%ld1, &
         dbl_grid%ld2, dbl_grid%max_group_slabs12, pub_num_spins)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    !kaw merge changereal(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%size_on_grid)
    integer, intent(in) :: col_start_in_box(3,batch_size)
    integer, intent(in) :: col_box_start(3,batch_size)

    ! cks: <<local variables>>
    integer :: is
    integer :: col_start1, col_start2, col_start3
    integer :: fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl
    integer :: col, local_col, col_prev, local_col_prev
    integer :: atom_of_col
    integer :: first_on_col_atom, last_on_col_atom
    integer :: batch_count
    logical :: i_have_box
    integer :: ierr

    integer :: atom_of_col_prev
    integer :: first_on_col_atom_prev, last_on_col_atom_prev
    logical :: gpu_sum, gpu_copy_off_previous, gpu_copy_off_current

    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: row_fftbox_dbl
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: row_fftbox_sum_batch
    !real(kind=DP), dimension(:,:,:) :: row_fftbox_sum_batch
    type(FFTBOX_DATA) :: col_fftbox
    !kaw merge change real(kind=DP), allocatable, dimension(:,:,:) :: col_fftbox
    real(kind=DP), allocatable, dimension(:,:,:,:) :: final_row_fftbox_dbl
    real(kind=DP), allocatable, dimension(:,:,:,:) :: buffer_dbl


    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering density_gpu_batch_interp_dep'

    ! agrecokpt: TO-DO: check/implement k-points here as well
    call utils_assert(pub_cmplx_ngwfs == .false., &
         'Subroutine density_gpu_batch_interp_dep not ready yet for &
         & complex NGWFs.')

    ! cks: loop over the members of the largest current batch (from all procs)
    ! cks: and interpolate each pair of fftboxes of each proc, multiply them
    ! cks: together and deposit them accordingly to the 12-slabs of all procs
    batch_count = 0
    do local_col=local_start,local_start+max_current_size-1
       batch_count = batch_count + 1

       col_start1 = col_start_in_box(1,local_col - local_start + 1)
       col_start2 = col_start_in_box(2,local_col - local_start + 1)
       col_start3 = col_start_in_box(3,local_col - local_start + 1)

       call timer_clock('density_fftbox_interpolate_multiply', 1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  SETUP VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! cks: interpolate and multiply my col NGWF if pub_my_proc_id
       ! cks: posseses an fftbox for this local_col value
       if ( local_col .le. local_end) then

          ! kaw: Gather information and set flags for current density
          col = local_col + ngwf_basis%first_on_proc(pub_my_proc_id) - 1
          atom_of_col = ngwf_basis%atom_of_func(col)
          first_on_col_atom = ngwf_basis%first_on_atom(atom_of_col)
          last_on_col_atom = first_on_col_atom + &
               ngwf_basis%num_on_atom(atom_of_col) - 1

          ! kaw: Toggle for copy off of data from current col ngwf if batch is complete
          gpu_copy_off_current = (local_col == local_end)

          ! kaw: Toggle for summation of density from an atom
          if ((col == first_on_col_atom).or.(local_col==local_start)) then
             gpu_sum = .false.
          else
             gpu_sum = .true.
          end if

          ! kaw: Gather information and set flags for previous density
          col_prev = col - 1
          local_col_prev = local_col - 1
          atom_of_col_prev = ngwf_basis%atom_of_func(col_prev)
          first_on_col_atom_prev = ngwf_basis%first_on_atom(atom_of_col_prev)
          last_on_col_atom_prev = first_on_col_atom_prev + &
               ngwf_basis%num_on_atom(atom_of_col_prev) - 1

          ! kaw: Toggle for copy off of data from previous atom if summation is complete
          gpu_copy_off_previous = (col_prev == last_on_col_atom_prev)

          if (local_col == local_start) gpu_copy_off_previous = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  CALCULATE DENSITY  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! ndmh: copy column NGWF to column fftbox
          call basis_copy_function_to_box(col_fftbox, &
               col_start1, col_start2, col_start3, &
               ngwf_basis%tight_boxes(local_col), ngwfs_on_grid, &
               ngwf_basis%spheres(local_col),cell,fftbox)

          ! ndmh: interpolate col and sum of rows together with c2c ffts
          ! ndmh: then multiply them together to find contribution to the
          ! ndmh: charge density

          do is=1,pub_num_spins
#ifdef GPU_SP_TEST
             if (GPU_SP) then
             call fourier_gpu_interp_prod_SP( &
                  row_fftbox_sum_batch(:,:,:,is,batch_count), col_fftbox, &
                  row_fftbox_dbl(:,:,:,is,2), final_row_fftbox_dbl(:,:,:,is),&
                  gpu_sum, gpu_copy_off_current, gpu_copy_off_previous, is)
             else
#endif
             call fourier_gpu_interp_prod( &
                  !row_fftbox_sum_batch, col_fftbox%d, &     !inputs
                  row_fftbox_sum_batch(:,:,:,is,batch_count), col_fftbox%d, &     !inputs
                  row_fftbox_dbl(:,:,:,is,2), final_row_fftbox_dbl(:,:,:,is),&  !outputs
                  gpu_sum, gpu_copy_off_current, gpu_copy_off_previous, is)
#ifdef GPU_SP_TEST
             end if
#endif
          end do
       else

          ! kaw: Toggle for copy off of data from current col ngwf if batch is complete
          gpu_copy_off_previous = .false.
          gpu_copy_off_current = .false.
       end if

       ! ndmh: synchronise procs in this rank so that load-balancing time is
       ! ndmh: reported in density_fftbox_interpolate_multiply rather than as
       ! ndmh: the wait for the alltoall in cell_grid_deposit_box.
       call comms_barrier(pub_rank_comm)
       call timer_clock('density_fftbox_interpolate_multiply', 2)
       call timer_clock('density_fftbox_deposit_to_cell', 1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  DISTRIBUTION OF PREVIOUS !!!!!!!!!!!!!!!!!!!!!!!!!!!!

       i_have_box = gpu_copy_off_previous

       ! cks: get into the depositing fftboxes to 12-slabs of all procs
       ! cks: regardless of if pub_my_proc_id has an fftbox to deposit
       ! ndmh: box_to_cell routine moved to basis_mod
       fftbox_start1_dbl = 2*col_box_start(1,local_col_prev-local_start+1) - 1
       fftbox_start2_dbl = 2*col_box_start(2,local_col_prev-local_start+1) - 1
       fftbox_start3_dbl = 2*col_box_start(3,local_col_prev-local_start+1) - 1

       do is=1,pub_num_spins
          call cell_grid_deposit_box(density_dbl(:,:,:,is), &
               row_fftbox_dbl(:,:,:,is,2), buffer_dbl, dbl_grid, &
               fftbox%total_pt1_dbl, fftbox%total_pt2_dbl, &
               fftbox%total_pt3_dbl, fftbox%total_ld1_dbl, &
               fftbox%total_ld2_dbl, fftbox_start1_dbl, &
               fftbox_start2_dbl, fftbox_start3_dbl, i_have_box, .true.)
       end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  DISTRIBUTION OF CURRENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!

       i_have_box = gpu_copy_off_current

       fftbox_start1_dbl = 2*col_box_start(1,local_col-local_start+1) - 1
       fftbox_start2_dbl = 2*col_box_start(2,local_col-local_start+1) - 1
       fftbox_start3_dbl = 2*col_box_start(3,local_col-local_start+1) - 1

       do is=1,pub_num_spins
          call cell_grid_deposit_box(density_dbl(:,:,:,is), &
               final_row_fftbox_dbl(:,:,:,is), buffer_dbl, dbl_grid, &
               fftbox%total_pt1_dbl, fftbox%total_pt2_dbl, &
               fftbox%total_pt3_dbl, fftbox%total_ld1_dbl, &
               fftbox%total_ld2_dbl, fftbox_start1_dbl, &
               fftbox_start2_dbl, fftbox_start3_dbl, i_have_box, .true.)
       end do

       call timer_clock('density_fftbox_deposit_to_cell', 2)

    end do

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving density_gpu_batch_interp_dep'
  end subroutine density_gpu_batch_interp_dep

#else

  subroutine dummy_gpu
  use fourier_gpu_wrapper_mod
  !kaw: This is a dummy use statement to make the check_dependencies script happy.
  end subroutine dummy_gpu

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module density

