! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                  Kinetic energy density module                 !
!                                                                !
!                  James C. Womack, 2015 - 2016                  !
!                                                                !
!----------------------------------------------------------------!
! A module providing subroutines for evaluation of the kinetic   !
! energy density (tau), in terms of the density kernel and       !
! NGWFs.                                                         !
!                                                                !
! The kinetic energy density is required for the implementation  !
! of meta-GGA exchange-correlation functionals.                  !
!                                                                !
! ke_density_on_grid and associated routines are modified        !
! versions of the density_on_grid subroutine (and others) from   !
! the density module.                                            !
!================================================================!
module ke_density

  use constants, only: DP, stdout, stderr
  use rundat, only: pub_debug_on_root

  implicit none
  private

  public :: ke_density_on_grid
  public :: weizsacker_ke_density_on_grid
  public :: ke_density_kinetic_energy

contains

  real(kind=DP) function ke_density_kinetic_energy(ke_density_fine,fine_grid) &
      result(kinetic_energy)
    !==========================================================================!
    ! This function evaluates the kinetic energy by integration over the       !
    ! the kinetic energy density, passed in as ke_density_fine.                !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! ke_density_fine (input) : The total kinetic energy density of the system !
    ! on the fine grid of the simulation cell                                  !
    ! fine_grid     (input) : GRID_INFO describing fine grid on which kinetic  !
    ! energy density is represented (may be finer than twice the psinc cutoff).!
    !                                                                          !
    ! Return value:                                                            !
    ! kinetic_energy (output) : Evaluated by integration over ke_density_fine  !
    ! and summation over spin components                                       !
    !--------------------------------------------------------------------------!
    ! Written by James C. Womack (2015)                                        !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: DP
    use integrals, only: integrals_trace_on_grid
    use rundat, only: pub_num_spins

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: fine_grid
    real(kind=DP), intent(in) :: ke_density_fine(fine_grid%ld1, &
         fine_grid%ld2, fine_grid%max_slabs12, pub_num_spins)

    ! Local variables
    integer :: spin_count

    ! kinetic energy calculated by integrating kinetic energy density
    kinetic_energy = 0.0_DP

    ! Integrate ke_density_fine on grid
    do spin_count = 1, pub_num_spins
       kinetic_energy = kinetic_energy + &
            integrals_trace_on_grid(ke_density_fine(:,:,:,spin_count),fine_grid)
    end do

  end function ke_density_kinetic_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ke_density_on_grid(ke_density_fine, &
       fine_grid, dbl_grid, cell, fftbox, denskern, overlap, &
       ngwfs_on_grid, ngwf_basis)
    !==========================================================================!
    ! This subroutine calculates the kinetic energy density on a fine grid in  !
    ! the simulation cell. Calls ke_density_on_dbl_grid to evaluate kinetic    !
    ! energy density on the double grid, and then uses                         !
    ! density_workspace_to_output (from density module) to interpolate to the  !
    ! scale of the fine grid for output.                                       !
    ! The kinetic energy density, \tau, is defined as                          !
    !   \tau = 0.5 * \sum_{\alpha\beta} K^{\alpha\beta}                        !
    !                (\nabla\phi_{\alpha}).(\nabla\phi_{\beta})                !
    ! with the factor of 0.5 included by convention (as in                     !
    ! Ernzerhof, M. & Scuseria, J. Chem. Phys. 111, 911-915 (1999) ).          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! ke_density_fine (output) : The total kinetic energy density on the fine  !
    ! grid of the simulation cell                                              !
    ! fine_grid     (input) : GRID_INFO describing fine grid on which kinetic  !
    !  energy density is eventually constructed                                !
    ! dbl_grid      (input) : GRID_INFO structure describing a grid with twice !
    !  the cutoff of the psinc grid                                            !
    ! cell (input) : CELL_INFO structure describing the simulation cell        !
    ! fftbox (input) : FFTBOX_INFO structure describing the fftbox used for    !
    !  FFT operations for calculating the kinetic energy density               !
    ! denskern      (input) : Density kernel sparse matrix                     !
    ! overlap       (input) : Overlap matrix                                   !
    ! ngwfs_on_grid (input) : The NGWFs on the real-space ppd representation   !
    !  that belong to this proc.                                               !
    ! ngwf_basis            : The function basis for the NGWFs.                !
    !--------------------------------------------------------------------------!
    ! Modified version of density_on_grid (Nicholas Hine 02/28/11), created by !
    ! James C. Womack, 2015.                                                   !
    !==========================================================================!

    use datatypes, only: FUNCTIONS
    use cell_grid, only: GRID_INFO
    use constants, only: DP
    use density, only: density_workspace_to_output
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_num_spins, pub_dbl_grid_scale, pub_dbl_is_std
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: denskern(pub_num_spins,1)
    type(SPAM3), intent(in) :: overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(GRID_INFO), intent(in) :: fine_grid
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(out) :: ke_density_fine(fine_grid%ld1, &
         fine_grid%ld2, fine_grid%max_slabs12, pub_num_spins)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid

    ! Local Variables
    integer :: ierr
    real(kind=DP), allocatable, dimension(:,:,:,:) :: ke_density_work

    ! JCW: Assert that double grid is has double the number of points as
    ! JCW: standard grid, i.e. dbl_grid_scale = 2.0 (default).
    ! JCW: ke_density_on_dbl_grid assumes that this is the case.
    call utils_assert(pub_dbl_grid_scale.eq.2.0_DP.and.(.not.pub_dbl_is_std),&
         "Error in ke_density_on_grid: calculation of KE density requires &
         &dbl_grid_scale = 2.0.")

    ! JCW: Allocate per-proc double grid real space work array.
    allocate(ke_density_work(dbl_grid%ld1,dbl_grid%ld2,&
         dbl_grid%max_group_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('ke_density_on_grid','ke_density_work',ierr)

    ! JCW: Calculate the kinetic energy density on the double grid
    call ke_density_on_dbl_grid(ke_density_work, dbl_grid, cell, fftbox, &
         denskern, overlap, ngwfs_on_grid, ngwf_basis)

    ! ndmh: copy or upscale to 'fine' grid as required
    ! JCW: (use original routine from density module, since the copy/upscale
    ! JCW: is identical to routine used for charge density).
    call density_workspace_to_output(ke_density_work,ke_density_fine, &
         dbl_grid,fine_grid)

    deallocate(ke_density_work,stat=ierr)
    call utils_dealloc_check('ke_density_on_grid','ke_density_work',ierr)

  end subroutine ke_density_on_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ke_density_on_dbl_grid(ke_density_dbl, dbl_grid, cell, fftbox, &
       denskern, overlap, ngwfs_on_grid, ngwf_basis)

    !==========================================================================!
    ! This subroutine calculates the kinetic energy density on the double grid !
    ! in the  simulation cell.                                                 !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! ke_density_dbl  (output) : The total kinetic energy density on the       !
    !  double grid of the simulation cell                                      !
    ! dbl_grid      (input) : GRID_INFO structure describing a grid with twice !
    !  the cutoff of the psinc grid                                            !
    ! cell   (input) : CELL_INFO structure describing the simulation cell      !
    ! fftbox (input) : FFTBOX_INFO structure describing the fftbox used for    !
    !  FFT operations for calculating the kinetic energy density               !
    ! denskern      (input) : Density kernel sparse matrix                     !
    ! overlap       (input) : Overlap matrix                                   !
    ! ngwfs_on_grid (input) : The NGWFs on the real-space ppd representation   !
    !  that belong to this proc                                                !
    ! ngwf_basis    (input) : The function basis for the NGWFs                 !
    !--------------------------------------------------------------------------!
    ! Modified version of density_on_dbl_grid for calculation of kinetic       !
    ! energy density, created by James C. Womack, 2015.                        !
    ! See below for comments from original density_on_dbl_grid subroutine      !
    ! relating to parallel efficiency and internal variables (which should     !
    ! apply to this modified subroutine, too).                                 !
    !==========================================================================!
    ! Comments from original density_on_dbl_grid routine:                      !
    ! JCW: It appears that "density_batch_row_sums" has been replaced by       !
    ! JCW: "function_ops_sum_fftbox_batch".                                    !
    !--------------------------------------------------------------------------!
    ! Key internal variables:                                                  !
    !   batch_size: This is set equal to pub_fftbox_batch_size which comes from    !
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
    ! from previous version of ke_density_on_fine_grid.                        !
    ! Modified to incorporate calls to GPU routines and dynamic precision by   !
    ! Karl Wilkinson in 2012                                                   !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero, &
         data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_interpolate_cell, fourier_apply_box, &
         fourier_apply_box_pair
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_sum_fftbox_batch, &
         function_ops_batch_col_start
    use kinetic, only: kinetic_grad_on_box
    use rundat, only: pub_fftbox_batch_size, pub_num_spins, pub_threads_fftbox, &
                      pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: denskern(pub_num_spins,1)
    type(SPAM3), intent(in) :: overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(out) :: ke_density_dbl(dbl_grid%ld1, &
         dbl_grid%ld2, dbl_grid%max_group_slabs12, pub_num_spins)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid

    ! Local Variables
    integer :: batch_size
    integer :: max_current_size
    integer :: batch_count1, batch_count2
    integer :: batch_number
    integer :: n_batches
    integer :: local_start, local_end
    integer :: iii, row1, row2
    integer :: is
    integer :: dim
    integer :: ierr    ! pdh: error flag
    integer :: idx_len
    integer, allocatable :: overlap_idx(:)            ! Index for overlap matrix
    integer, allocatable :: col_box_start(:,:)
    integer, allocatable :: col_start_in_box(:,:)
    ! agrecocmplx: extra variable to use complex quantities
    logical :: loc_cmplx


    type(FFTBOX_DATA), allocatable, dimension(:,:) :: row_fftbox_sum_batch
    type(FFTBOX_DATA), allocatable, dimension(:,:,:) :: &
         row_fftbox_sum_grad_batch
    complex(kind=DP), allocatable, dimension(:,:,:)  :: coarse_work
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: zwork_grad_box
    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering ke_density_on_dbl_grid'

    call timer_clock('ke_density_on_dbl_grid',1)

    ! Obtain index for overlap matrix
    idx_len = sparse_index_length(overlap)
    allocate(overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('ke_density_on_dbl_grid','overlap_idx',ierr)
    call sparse_generate_index(overlap_idx,overlap)

    ! Set batch size according to input file
    batch_size = pub_fftbox_batch_size

    ! agrecocmplx: optional variable for complex NGWFs
    loc_cmplx = ngwfs_on_grid%iscmplx
    ! JCW: ke_density_on_dbl_grid has not been updated to support complex NGWFs
    ! JCW: yet
    call utils_assert(.not.loc_cmplx, "Error in ke_density_on_dbl_grid: &
         &Complex NGWFs not currently supported for ke_density evaluation.")

    ! Allocate workspace
    allocate(col_box_start(3, batch_size), stat=ierr)
    call utils_alloc_check('ke_density_on_dbl_grid','col_box_start',ierr)
    allocate(col_start_in_box(3, batch_size), stat=ierr)
    call utils_alloc_check('ke_density_on_dbl_grid','col_start_in_box',ierr)
    allocate(row_fftbox_sum_batch(pub_num_spins,batch_size),stat=ierr)
    call utils_alloc_check('ke_density_on_dbl_grid','row_fftbox_sum_batch',ierr)
    do is=1, pub_num_spins
        do batch_count1=1, batch_size
            call data_fftbox_alloc(row_fftbox_sum_batch(is,batch_count1), &
                fftbox%total_ld1, fftbox%total_ld2, fftbox%total_pt3, &
                iscmplx=loc_cmplx)
        end do
    end do


    ! Allocate workspace for taking gradient
    allocate(row_fftbox_sum_grad_batch(3,pub_num_spins,batch_size),stat=ierr)
    call utils_alloc_check('ke_density_on_dbl_grid',&
         'row_fftbox_sum_grad_batch',ierr)
    do is=1, pub_num_spins
        do batch_count1=1, batch_size
           do dim=1,3
              call data_fftbox_alloc(&
                   row_fftbox_sum_grad_batch(dim,is,batch_count1), &
                   fftbox%total_ld1, fftbox%total_ld2, fftbox%total_pt3, &
                   iscmplx=loc_cmplx)
           end do
        end do
    end do

    ! cks: zero before accumulation
    ke_density_dbl = 0.0_DP

    ! cks: number of row-steps per row-block
    n_batches = ngwf_basis%max_on_proc / batch_size
    if (mod(ngwf_basis%max_on_proc, batch_size) > 0) n_batches = n_batches + 1

    ! cks: loop over batches of NGWFs belonging to pub_my_proc_id
    local_start = 1
    do batch_number=1,n_batches

       if (pub_debug_on_root) write(stdout,'(a,2(i4,a))') 'DEBUG: Batch loop ',&
          batch_number,' of ',n_batches,' in ke_density_on_fine_grid'

       ! cks: limits of my current batch in pub_my_proc_id-local NGWF counting
       ! cks: scheme
       local_end = min(local_start+batch_size-1,ngwf_basis%proc_num)

       ! cks: maximum size of current batch over all procs
       max_current_size = local_end - local_start + 1
       call comms_reduce('MAX', max_current_size)

       ! ndmh: find start positions of functions in this batch
       call function_ops_batch_col_start(col_box_start,col_start_in_box, &
            batch_size,local_start,local_end,fftbox,cell, &
            ngwf_basis)

       ! cks: zero before accumulation
       ! agrecocmplx: initialise accordingly
       do is=1, pub_num_spins
           do batch_count1=1, batch_size
               call data_set_to_zero(row_fftbox_sum_batch(is,batch_count1))
           end do
       end do


       call timer_clock('ke_density_batch_row_sums',1)

       ! ndmh: density-kernel block contributions to fftboxes of current batch
       ! JCW: This involves inter-proc communication
       call function_ops_sum_fftbox_batch(row_fftbox_sum_batch, &
            fftbox, cell, ngwfs_on_grid, ngwf_basis, &
            col_box_start, batch_size, local_start, local_end, &
            overlap_idx, idx_len, denskern, 1, 1, 1.0_DP)

       call timer_clock('ke_density_batch_row_sums',2)

       call timer_clock('ke_density_batch_grad_apply',1)
!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(coarse_work,zwork_grad_box, &
!$OMP         is,iii,ierr,row1,row2,batch_count1,batch_count2 ) &
!$OMP SHARED(row_fftbox_sum_batch,row_fftbox_sum_grad_batch, &
!$OMP        fftbox, pub_num_spins,pub_threads_fftbox, &
!$OMP        pub_threads_num_fftboxes, local_start, local_end )

       ! JCW: Allocate per-FFTbox-batch work arrays
       allocate(coarse_work(fftbox%total_ld1, &
            fftbox%total_ld2,fftbox%total_pt3),stat=ierr)
       call utils_alloc_check('ke_density_on_dbl_grid','coarse_work',&
            ierr)
       allocate(zwork_grad_box(fftbox%total_ld1,fftbox%total_ld2, &
            fftbox%total_pt3,3),stat=ierr)
       call utils_alloc_check('ke_density_on_dbl_grid',&
            'zwork_grad_box',ierr)


       ! JCW: Apply kinetic_grad_on_box to each FFTbox in batch
       do is = 1, pub_num_spins
!$OMP DO
          do iii = local_start, local_end, 2
             ! JCW: Loop increments by 2, so that the fourier transform can be
             ! JCW: applied in pairs, via fourier_apply_box_pair
             row1 = iii
             row2 = iii+1
             batch_count1 = row1 - local_start + 1
             batch_count2 = row2 - local_start + 1
             ! JCW: Copy row_fftbox(es) to complex coarse_work array
             if ( row2 <= local_end ) then
                ! JCW: A pair of FFTboxes can be transformed simultaneously
                coarse_work = cmplx(&
                     row_fftbox_sum_batch(is,batch_count1)%d(:,:,:),&
                     row_fftbox_sum_batch(is,batch_count2)%d(:,:,:),&
                     kind=DP)
             else
                ! JCW: There is only one FFTbox left to transform
                coarse_work = cmplx(&
                     row_fftbox_sum_batch(is,batch_count1)%d(:,:,:),&
                     0.0_DP,kind=DP)
             end if

             ! JCW: Fourier transform FFTboxes to reciprocal space
             ! JCW: Two real FFTs can be done simultaneously, with one
             ! JCW: real function in the real part of coarse_work and
             ! JCW: the other real function in the imaginary part of
             ! JCW: coarse_work.
             call fourier_apply_box('Coarse','Forward',coarse_work,&
                  omp=pub_threads_fftbox)

             ! JCW: Apply gradient operator in reciprocal space
             call kinetic_grad_on_box(coarse_work,zwork_grad_box,fftbox)

             do dim=1,3
                ! JCW: Reverse transform FFTboxes to real space following
                ! JCW: application of gradient operator
                call fourier_apply_box('Coarse','Backward',&
                     zwork_grad_box(:,:,:,dim), omp=pub_threads_fftbox)
                ! JCW: Copy real components of zwork_grad_box into
                ! JCW: row_fftbox_sum_grad_batch
                row_fftbox_sum_grad_batch(dim,is,batch_count1)%d(:,:,:) = &
                     real(zwork_grad_box(:,:,:,dim),kind=DP)
                if (row2 <= local_end) then
                   ! JCW: Second fftbox is in imaginary part of zwork_grad_box
                   row_fftbox_sum_grad_batch(dim,is,batch_count2)%d(:,:,:) = &
                        real(aimag(zwork_grad_box(:,:,:,dim)),kind=DP)
                end if
             end do
          end do
!$OMP ENDDO
       end do


       ! JCW: Deallocate per-FFTbox-batch work arrays
       deallocate(zwork_grad_box,stat=ierr)
       call utils_dealloc_check('ke_density_on_dbl_grid',&
            'zwork_grad_box',ierr)
       deallocate(coarse_work,stat=ierr)
       call utils_dealloc_check('ke_density_on_dbl_grid','coarse_work',ierr)

!$OMP END PARALLEL
       call timer_clock('ke_density_batch_grad_apply',2)

       call timer_clock('ke_density_batch_interp_deposit',1)

       ! cks: interpolation of fftboxes of current batch, multiplication
       ! cks: and accumulation to charge-density 12-slabs of each proc
       call ke_density_batch_interp_deposit(ke_density_dbl,dbl_grid, &
            cell,fftbox,row_fftbox_sum_grad_batch,ngwfs_on_grid, &
            ngwf_basis,col_box_start,col_start_in_box,local_start,&
            local_end,batch_size,max_current_size)

       call timer_clock('ke_density_batch_interp_deposit',2)

       ! Update first NGWF of next batch
       local_start = local_start + batch_size

    end do

    ! Deallocate workspace for taking gradient
    ! agrecocmplx: deallocate using appropriate routine
    do is=1, pub_num_spins
        do batch_count1=1, batch_size
           do dim=1,3
              call data_fftbox_dealloc(row_fftbox_sum_grad_batch(dim,is,&
                   batch_count1))
           end do
        end do
    end do
    deallocate(row_fftbox_sum_grad_batch,stat=ierr)
    call utils_dealloc_check('ke_density_on_dbl_grid',&
         'row_fftbox_sum_grad_batch',ierr)

    ! Deallocate workspace
    ! agrecocmplx: deallocate using appropriate routine
    do is=1, pub_num_spins
        do batch_count1=1, batch_size
            call data_fftbox_dealloc(row_fftbox_sum_batch(is,batch_count1))
        end do
    end do
    deallocate(row_fftbox_sum_batch,stat=ierr)
    call utils_dealloc_check('ke_density_on_dbl_grid','row_fftbox_sum_batch',&
         ierr)
    deallocate(col_start_in_box,stat=ierr)
    call utils_dealloc_check('ke_density_on_dbl_grid','col_start_in_box',ierr)
    deallocate(col_box_start,stat=ierr)
    call utils_dealloc_check('ke_density_on_dbl_grid','col_box_start',ierr)
    deallocate(overlap_idx,stat=ierr)
    call utils_dealloc_check('ke_density_on_dbl_grid','overlap_idx',ierr)

    call timer_clock('ke_density_on_dbl_grid',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving ke_density_on_dbl_grid'

  end subroutine ke_density_on_dbl_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ke_density_batch_interp_deposit(ke_density_dbl, &
       dbl_grid, cell, fftbox, row_fftbox_sum_grad_batch, ngwfs_on_grid, &
       ngwf_basis, col_box_start, col_start_in_box, local_start, local_end, &
       batch_size, max_current_size)

    !==========================================================================!
    ! This subroutines interpolates each of the sums-in-fftboxes (with grad    !
    ! operator applied) of the batch and its corresponding col NGWF. Then it   !
    ! takes the dot product of the gradients of the NGWF and gradients of the  !
    ! sums-in-fftboxes to produce the kinetic energy density,                  !
    !   \tau = 0.5 * \sum_{\alpha\beta} K^{\alpha\beta}                        !
    !                (\nabla\phi_{\alpha}).(\nabla\phi_{\beta})                !
    !  which is deposited on the fine grid of the simulation cell.             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !  ke_density_dbl (input/output): Charge density in double-grid real space !
    !    representation.                                                       !
    !  dbl_grid     (input) : GRID_INFO structure describing a grid with twice !
    !  the cutoff of the psinc grid                                            !
    !  cell   (input) : CELL_INFO structure describing the simulation cell     !
    !  fftbox (input) : FFTBOX_INFO structure describing the fftbox used for   !
    !  row_fftbox_sum_grad_batch (input): Batch of fftboxes for each Cartesian !
    !    direction. Each fftbox contains the quantity \nabla_{i} phi_b*K^{ba}  !
    !    for the function phi_a with Cartesian direction i = x, y or z         !
    !  ngwfs_on_grid (input): Current NGWFs for pub_my_proc_id in ppd          !
    !    representation.                                                       !
    !  ngwf_basis (input) : The function basis for the NGWFs                   !
    !  col_box_start(input): list of column function fftbox start positions in !
    !    the simulation cell                                                   !
    !  col_start_in_box (input): list of column function start positions in an !
    !    fftbox, i.e. where each NGWF starts in an fftbox.                     !
    !  local_start (input): Number of the first NGWF of the current batch      !
    !    in the counting scheme of all NGWFs of pub_my_proc_id.                !
    !  local_end (input): Number of the last NGWF of the current batch         !
    !    in the counting scheme of all NGWFs of pub_my_proc_id.                !
    !  batch_size (input): Number of fftboxes (phi_a NGWFs) in each batch      !
    !  max_current_size (input): maximum size of current batch for all procs   !
    !--------------------------------------------------------------------------!
    ! Modified version of density_batch_interp_deposit for calculation of      !
    ! kinetic energy density, created by James C. Womack, 2015.                !
    ! See below for author details from original density_batch_interp_deposit  !
    ! subroutine.                                                              !
    !==========================================================================!
    ! Author details from original density_batch_interp_deposit routine:       !
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
    !==========================================================================!

    use basis, only: basis_copy_function_to_box, basis_location_func_wrt_cell
    use cell_grid, only: GRID_INFO, cell_grid_deposit_box
    use comms, only: comms_barrier, pub_my_proc_id, pub_rank_comm
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, &
         data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box, fourier_interpolate_product
    use function_basis, only: FUNC_BASIS
    use kinetic, only: kinetic_grad_on_box
    use rundat, only: pub_threads_num_fftboxes, pub_num_spins, &
         pub_threads_per_fftbox, pub_threads_fftbox
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: local_start
    integer, intent(in) :: local_end
    integer, intent(in) :: batch_size
    integer, intent(in) :: max_current_size
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(inout) :: ke_density_dbl(dbl_grid%ld1, &
         dbl_grid%ld2, dbl_grid%max_group_slabs12, pub_num_spins)
    type(FFTBOX_DATA), intent(in) :: row_fftbox_sum_grad_batch(&
          3,pub_num_spins,batch_size)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNCTIONS) :: ngwfs_on_grid
    integer, intent(in) :: col_start_in_box(3,batch_size)
    integer, intent(in) :: col_box_start(3,batch_size)

    ! Local Variables
    integer :: is
    integer :: dim
    integer :: fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl
    integer :: col, local_col
    integer :: atom_of_col
    integer :: first_on_col_atom, last_on_col_atom
    integer :: batch_count
    integer :: ierr
    logical :: i_have_box
    real(kind=DP), allocatable, dimension(:,:,:) :: buffer_dbl
    type(FFTBOX_DATA) :: col_fftbox
    real(kind=DP), allocatable, dimension(:,:,:,:) :: row_fftbox_dbl
    real(kind=DP), allocatable, dimension(:,:,:) :: fftbox_dbl_work
    complex(kind=DP), allocatable, dimension(:,:,:) :: coarse_work,fine_work
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: zwork_grad_box
    ! agrecocmplx: extra variable for complex quantities
    logical :: loc_cmplx
    integer :: thread_id,nthreads
    integer :: thread_local_start, thread_local_end
!$  integer, external :: omp_get_thread_num, omp_get_num_threads

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering ke_density_batch_interp_deposit'

    ! agrecocmplx: optional variable for complex NGWFs
    loc_cmplx = ngwfs_on_grid%iscmplx
    ! JCW: ke_density_batch_interp_deposit has not been updated to support
    ! JCW: complex NGWFs yet
    call utils_assert(.not.loc_cmplx, &
         "Error in ke_density_batch_interp_deposit: &
         &Complex NGWFs not currently supported for ke_density evaluation.")

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(local_col, col, atom_of_col, i_have_box, &
!$OMP      first_on_col_atom,  last_on_col_atom, batch_count, &
!$OMP      is, fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl,&
!$OMP      row_fftbox_dbl, buffer_dbl,ierr, &
!$OMP      thread_local_start, thread_local_end, thread_id, nthreads, &
!$OMP      coarse_work, fine_work, zwork_grad_box, fftbox_dbl_work, dim ) &
!$OMP FIRSTPRIVATE(col_fftbox) &
!$OMP SHARED(row_fftbox_sum_grad_batch, ke_density_dbl, ngwfs_on_grid, &
!$OMP      ngwf_basis, dbl_grid, col_box_start, col_start_in_box, fftbox, &
!$OMP      cell, pub_num_spins, local_start, local_end, max_current_size, &
!$OMP      pub_my_proc_id, pub_rank_comm, pub_threads_num_fftboxes, &
!$OMP      pub_threads_fftbox, pub_threads_per_fftbox,loc_cmplx)

    ! Allocate Workspace
    ! agrecocmplx: use approproate routine
    call data_fftbox_alloc(col_fftbox, fftbox%total_ld1, &
             fftbox%total_ld2,fftbox%total_pt3, iscmplx=loc_cmplx)
    allocate(row_fftbox_dbl(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl,2),stat=ierr)
    call utils_alloc_check('ke_density_batch_interp_deposit','row_fftbox_dbl',&
         ierr)
    allocate(fftbox_dbl_work(fftbox%total_ld1_dbl, &
         fftbox%total_ld2_dbl,fftbox%total_pt3_dbl),stat=ierr)
    call utils_alloc_check('ke_density_batch_interp_deposit','fftbox_dbl_work',&
         ierr)
    allocate(buffer_dbl(fftbox%total_ld1_dbl, fftbox%total_ld2_dbl, &
         dbl_grid%max_group_slabs12),stat=ierr)
    call utils_alloc_check('ke_density_batch_interp_deposit','buffer_dbl',ierr)
    allocate(coarse_work(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('ke_density_batch_interp_deposit','coarse_work',ierr)
    allocate(fine_work(fftbox%total_ld1*2, &
         fftbox%total_ld2*2,fftbox%total_pt3*2),stat=ierr)
    call utils_alloc_check('ke_density_batch_interp_deposit','fine_work',ierr)
    allocate(zwork_grad_box(fftbox%total_ld1,fftbox%total_ld2, &
         fftbox%total_pt3,3),stat=ierr)
    call utils_alloc_check('ke_density_batch_interp_deposit','zwork_grad_box',&
         ierr)

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
       ! cks: loop over the members of the largest current batch (from all
       ! cks: procs) and interpolate each pair of fftboxes of each proc,
       ! cks: multiply them together and deposit them accordingly to the
       ! cks: 12-slabs of all procs
       do local_col=thread_local_start,thread_local_end
          batch_count = local_col - local_start + 1

          call timer_clock('ke_density_fftbox_interpolate_multiply', 1)

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
             call basis_copy_function_to_box(col_fftbox, &
                  col_start_in_box(1,batch_count), &
                  col_start_in_box(2,batch_count), &
                  col_start_in_box(3,batch_count), &
                  ngwf_basis%tight_boxes(local_col), ngwfs_on_grid, &
                  ngwf_basis%spheres(local_col), cell, fftbox)

             ! JCW: Copy col_fftbox%d to complex coarse_work array
             coarse_work(:,:,:) = cmplx(col_fftbox%d,0.0_DP,kind=DP)

             ! JCW: Fourier transform a single NGWF in an FFTbox
             call fourier_apply_box('Coarse','Forward',coarse_work,&
                  omp=pub_threads_fftbox)

             ! JCW: Apply kinetic_grad_on_box to a single NGWF (1 batch)
             call kinetic_grad_on_box(coarse_work,zwork_grad_box,fftbox)

             ! JCW: Zero row_fftbox_dbl(:,:,:,1), since we will accumulate into
             ! JCW: it.
             row_fftbox_dbl(:,:,:,1) = 0.0_DP
             do dim=1,3
                call fourier_apply_box('Coarse','Backward',&
                     zwork_grad_box(:,:,:,dim), omp=pub_threads_fftbox)
                ! JCW: If double grid is not the same as standard grid, then
                ! JCW: evaluate product with interpolation.
                ! JCW: We can assume this because it is asserted in
                ! JCW: ke_density_on_grid
                ! ndmh: interpolate col and sum of rows together with c2c ffts
                ! ndmh: then multiply them together to find contribution to the
                ! ndmh: charge density

                ! JCW: FIX TO AVOID OPENMP STACKSIZE ISSUES
                ! JCW: Typecasting zwork_grad_box to a real as a subroutine
                ! JCW: argument causes a crash with >1 OpenMP threads sometimes
                ! JCW: with Intel compiler and OpenMP.
                ! JCW: Increading KMP_STACKSIZE environment variable above
                ! JCW: (small) default fixes this issue. It seems likely that
                ! JCW: the thread stack is used for temporary space for the
                ! JCW: real component of zwork_grad_box.
                ! JCW: For safety, we should perform the function call in a way
                ! JCW: which does not need to temporarily allocate memory in
                ! JCW: this way.
                !call fourier_interpolate_product(coarse_work, fine_work, &
                !     row_fftbox_sum_grad_batch(dim,is,batch_count)%d(:,:,:), &
                !     real(zwork_grad_box(:,:,:,dim),kind=DP), fftbox_dbl_work)
                ! JCW: Copy real component of zwork_grad_box to unused
                ! JCW: unused col_fftbox%d before calling routine.
                col_fftbox%d(:,:,:) = real(zwork_grad_box(:,:,:,dim),kind=DP)
                call fourier_interpolate_product(coarse_work, fine_work, &
                     row_fftbox_sum_grad_batch(dim,is,batch_count)%d, &
                     col_fftbox%d, fftbox_dbl_work)
                ! JCW: This alternative method is robust for small OpenMP stack
                ! JCW: sizes.

                ! JCW: Accumulate product for each dimension
                ! JCW: (multiply by 1/2 for \tau defined in line with
                ! JCW: convention
                row_fftbox_dbl(:,:,:,1) = row_fftbox_dbl(:,:,:,1) + &
                      0.5_DP * fftbox_dbl_work(:,:,:)
             end do


             ! ndmh: if this is a new atom or the start of a new batch, reset
             ! ndmh: the accumulated density for this atom to the density for
             ! ndmh: this NGWF
             if ((col == first_on_col_atom).or.(local_col==thread_local_start))&
                  then
                row_fftbox_dbl(:,:,:,2) = row_fftbox_dbl(:,:,:,1)
             else ! ndmh: add the density for this NGWF to the sum for this atom
                row_fftbox_dbl(:,:,:,2) = row_fftbox_dbl(:,:,:,2) + &
                     row_fftbox_dbl(:,:,:,1)
             end if

             ! ndmh: if this is the last NGWF on this atom or the end of the
             ! ndmh: batch, we will need to deposit the density to the
             ! ndmh: whole-cell array
             if ((col == last_on_col_atom).or.(local_col==thread_local_end)) &
                 then
                i_have_box = .true.
             else
                i_have_box = .false.
             end if

          else
             i_have_box =.false.
          end if


          ! ndmh: synchronise procs in this rank so that load-balancing time is
          ! ndmh: reported in ke_density_fftbox_interpolate_multiply rather than
          ! ndmh: as the wait for the alltoall in cell_grid_deposit_box.
!$        if (.false.) then
          call comms_barrier(pub_rank_comm)
!$        end if
          call timer_clock('ke_density_fftbox_interpolate_multiply', 2)
          call timer_clock('ke_density_fftbox_deposit_to_cell', 1)

          ! cks: get into the depositing fftboxes to 12-slabs of all procs
          ! cks: regardless of if pub_my_proc_id has an fftbox to deposit
          ! ndmh: box_to_cell routine moved to basis_mod
          fftbox_start1_dbl = 2*col_box_start(1,batch_count) - 1
          fftbox_start2_dbl = 2*col_box_start(2,batch_count) - 1
          fftbox_start3_dbl = 2*col_box_start(3,batch_count) - 1
!$OMP CRITICAL
          call cell_grid_deposit_box(ke_density_dbl(:,:,:,is), &
               row_fftbox_dbl(:,:,:,2), buffer_dbl, dbl_grid, &
               fftbox%total_pt1_dbl, fftbox%total_pt2_dbl, &
               fftbox%total_pt3_dbl, fftbox%total_ld1_dbl, &
               fftbox%total_ld2_dbl, fftbox_start1_dbl, &
               fftbox_start2_dbl, fftbox_start3_dbl, i_have_box, .true.)
!$OMP END CRITICAL

          call timer_clock('ke_density_fftbox_deposit_to_cell', 2)

       end do

!$OMP BARRIER

    end do


    deallocate(zwork_grad_box,stat=ierr)
    call utils_dealloc_check('ke_density_batch_interp_deposit',&
         'zwork_grad_box',ierr)
    deallocate(fine_work,stat=ierr)
    call utils_dealloc_check('ke_density_batch_interp_deposit',&
         'fine_work',ierr)
    deallocate(coarse_work,stat=ierr)
    call utils_dealloc_check('ke_density_batch_interp_deposit',&
         'coarse_work',ierr)

    deallocate(buffer_dbl,stat=ierr)
    call utils_dealloc_check('ke_density_batch_interp_deposit',&
         'buffer_dbl',ierr)
    deallocate(fftbox_dbl_work,stat=ierr)
    call utils_dealloc_check('ke_density_batch_interp_deposit',&
         'fftbox_dbl_work',ierr)
    deallocate(row_fftbox_dbl,stat=ierr)
    call utils_dealloc_check('ke_density_batch_interp_deposit',&
         'row_fftbox_dbl',ierr)
    ! agrecocmplx: deallocate using appropriate routine
    call data_fftbox_dealloc(col_fftbox)
!$OMP END PARALLEL

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving ke_density_batch_interp_deposit'

  end subroutine ke_density_batch_interp_deposit

  subroutine weizsacker_ke_density_on_grid(w_ke_density_fine, &
       density_fine, fine_grid )
    !==========================================================================!
    !                  EXPERIMENTAL ROUTINE --- UNTESTED!                      !
    ! This subroutine calculates the von Weizsacker approximation to the       !
    ! kinetic energy density,                                                  !
    !   \tau^{W} = (1/8) * |grad(rho)|^{2} / rho                               !
    ! on a fine grid using the charge density.                                 !
    !                                                                          !
    ! The von Weizsacker kinetic energy functional was developed in the 1930s  !
    ! by Weizsacker, and is described in many more recent papers, e.g.         !
    ! Acharya, P. K., Bartolotti, L. J., Sears, S. B. & Parr, R. G.            !
    !                                               PNAS 77, 6978-6982 (1980). !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! w_ke_density_fine (output) : The total von Weizsacker kinetic energy     !
    !  density on the fine grid of the simulation cell                         !
    ! density_fine  (input) : The charge density on the fine grid of the       !
    !  simulation cell                                                         !
    ! fine_grid     (input) : GRID_INFO describing fine grid on which kinetic  !
    !  energy density is eventually constructed                                !
    !--------------------------------------------------------------------------!
    ! Written by James C. Womack, May 2016.                                    !
    !==========================================================================!

    use constants, only: DP
    use cell_grid, only: GRID_INFO
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_num_spins
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert
    use xc, only: xc_gradients

    type(GRID_INFO), intent(in) :: fine_grid
    real(kind=DP), intent(out) :: w_ke_density_fine(fine_grid%ld1, &
         fine_grid%ld2, fine_grid%max_slabs12, pub_num_spins)
    real(kind=DP), intent(inout) :: density_fine(fine_grid%ld1, &
         fine_grid%ld2, fine_grid%max_slabs12, pub_num_spins)

    real(kind=DP), allocatable :: density_grad(:,:,:,:,:)
    complex(kind=DP), allocatable :: recip_work(:,:,:,:)
    integer :: i1, i2, islab12, is
    integer :: ierr
    real(kind=DP) :: o_o_den

    ! Allocate workspace for calculating gradient of the density
    allocate(density_grad(fine_grid%ld1,fine_grid%ld2,fine_grid%max_slabs12,3, &
         pub_num_spins),stat=ierr)
    call utils_alloc_check('weizsacker_ke_density_fine', &
         'density_grad',ierr)

    ! Allocate reciprocal space workspace
    allocate(recip_work(fine_grid%ld3,fine_grid%ld2,fine_grid%max_slabs23,3),stat=ierr)
    call utils_alloc_check('internal_allocate_workspace (xc)', &
         'recip_work',ierr)

    ! Get density gradient
    call xc_gradients(density_fine,density_grad,recip_work, &
                      fine_grid,0)

    ! Square and sum Cartesian components of density_grad
    density_grad(:,:,:,1,:) = &
         density_grad(:,:,:,1,:)*density_grad(:,:,:,1,:) + &
         density_grad(:,:,:,2,:)*density_grad(:,:,:,2,:) + &
         density_grad(:,:,:,3,:)*density_grad(:,:,:,3,:)
    ! ...density_grad(:,:,:,1,:) now contains |grad(rho)|^{2}

    ! Calculate von Weizacker kinetic energy density
    !   \tau^{W} = (1/8) * |grad(rho)|^{2} / rho
    do is=1,pub_num_spins
       do islab12=1,fine_grid%num_my_slabs12
          do i2=1,fine_grid%n2
             do i1=1,fine_grid%n1
                if (density_fine(i1,i2,islab12,is) > 0.0_DP) then
                   o_o_den = 1.0_DP / density_fine(i1,i2,islab12,is)
                else
                   o_o_den = 0.0_DP
                end if
                w_ke_density_fine(i1,i2,islab12,is) = &
                     0.125_DP * density_grad(i1,i2,islab12,1,is) &
                     /  density_fine(i1,i2,islab12,is)
             end do
          end do
       end do
    end do

    ! Deallocate workspaces
    deallocate(density_grad,stat=ierr)
    call utils_dealloc_check('weizsacker_ke_density_fine', &
         'density_grad',ierr)
    deallocate(recip_work,stat=ierr)
    call utils_dealloc_check('weizsacker_ke_density_fine', &
         'recip_work',ierr)
  end subroutine weizsacker_ke_density_on_grid

end module
