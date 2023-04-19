! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by:
!
!   Chris-Kriton Skylaris, Arash A. Mostofi, Nicholas D.M. Hine,
!   Simon M.M.Dubois and Jacek Dziedzic.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module integrals

  use constants, only: DP
  use rundat, only: pub_debug_on_root

  implicit none

  private

  public :: integrals_confinement_pot ! gcc32
  public :: integrals_locpot
  public :: integrals_div_locpot_grad_functions ! JCW
  public :: test_integrals_div_locpot_grad_functions ! JCW (temporary)
  public :: integrals_kinetic
  public :: integrals_grad
  public :: integrals_pos
  public :: integrals_exr
  public :: integrals_trace_on_grid
  public :: integrals_product_on_grid
  public :: integrals_product_on_grid_vec

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine integrals_confinement_pot_dbl_grid(locpot,  &       ! input-output
       bras_on_grid, bra_basis, kets_on_grid, ket_basis, &
       dbl_grid, cell, fftbox)

    !==========================================================================!
    ! This subroutine applies a local confining potential to each NGWF and     !
    ! calculates the braket:                                                   !
    ! confinement_pot_\alpha\beta = < bra_\alpha | ConfPot | ket_\beta >       !
    ! The non-diagonal terms are removed, as each NGWF only feels its own      !
    ! confinement                                                              !
    !==========================================================================!
    !  Arguments:                                                              !
    !    locpot     (inout)  : diagonal sparse local potential matrix elements !
    !    bras_on_grid  (in)  : bra functions on grid in PPD format             !
    !    bra_basis (input)   : The basis type for the bras                     !
    !    kets_on_grid  (in)  : ket functions on grid in PPD format             !
    !    ket_basis (input)   : The basis type for the kets                     !
    !==========================================================================!
    ! Derived from integrals_locpot_dbl_grid by Gabriel Constantinescu in      !
    ! July 2015                                                                !
    !==========================================================================!


    use cell_grid, only: GRID_INFO
    use comms, only: comms_barrier, comms_reduce, pub_on_root, pub_my_proc_id
    use constants, only: DP, stdout, cmplx_0
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_batch_col_start, &
         function_ops_brappd_ketfftbox
    use parallel_strategy, only: PARAL_INFO
    use potential, only: potential_confine_ngwf_batch
    use rundat, only: pub_fftbox_batch_size, pub_locpot_scheme
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_generate_index, sparse_index_length, &
         sparse_scale, sparse_create, sparse_transpose, sparse_axpy, &
         sparse_destroy, sparse_expand, pattern_lower, pattern_alternate, &
         sparse_put_element, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: locpot
    type(FUNC_BASIS), intent(in) :: bra_basis
    type(FUNCTIONS), intent(in) :: bras_on_grid
    type(FUNC_BASIS), intent(in) :: ket_basis
    type(FUNCTIONS), intent(in) :: kets_on_grid
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local Variables
    type(FFTBOX_DATA), allocatable, dimension(:) :: potential_fftbox_batch
    integer :: batch_size
    integer :: batch_count
    integer :: n_batches
    integer :: local_start, local_end
    integer :: max_current_size ! maximum batch size ovar all procs
    integer :: idx_len  ! pdh: sparse matrix index length
    integer :: ierr ! pdh: error flag
    integer, allocatable, dimension(:) :: locpot_idx ! pdh: sparse index
    integer, allocatable :: ket_box_start(:,:)
    integer, allocatable :: ket_start_in_box(:,:)
    type(SPAM3) :: locpot_transpose  ! ndmh: for symmetrising locpot
    character(len=80) :: locpot_scheme
    logical :: loc_cmplx, ongpu

    integer :: loc_bra,jket,ibra,iat,loc_iat
    type(PARAL_INFO), pointer :: row_par ! rc2013: column parallel strategy


    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering integrals_confinement_pot_dbl_grid'

    ! Start timer
    call timer_clock('integrals_confinement_pot_dbl_grid',1)

    ! agrecocmplx
    loc_cmplx = kets_on_grid%iscmplx.or.bras_on_grid%iscmplx
    ! Set batch size
    batch_size = pub_fftbox_batch_size

    ! rc2013: obtain parallel strategy for this SPAM3
    call sparse_get_par(row_par, locpot, 'R')

    ! pdh: obtain index for locpot
    idx_len = sparse_index_length(locpot)
    allocate(locpot_idx(idx_len),stat=ierr)
    call utils_alloc_check('integrals_confinement_pot_dbl_grid','locpot_idx',ierr)
    call sparse_generate_index(locpot_idx,locpot)

    ! rc2013: check that the bases refer to the same region
    if (bra_basis%name/=ket_basis%name .or. bra_basis%ireg/=ket_basis%ireg) then
       locpot_scheme = 'FULL'
    else
       locpot_scheme = pub_locpot_scheme
    end if

    ! Allocate workspace
    allocate(potential_fftbox_batch(batch_size), stat=ierr)
    call utils_alloc_check('integrals_confinement_pot_dbl_grid', &
           'potential_fftbox_batch',ierr)
    do batch_count =1, batch_size
        ! agrecocmplx: potential_fftbox_batch is complex when NGWFs are complex
        call data_fftbox_alloc(potential_fftbox_batch(batch_count), &
             fftbox%total_ld1, fftbox%total_ld2, fftbox%total_pt3, loc_cmplx) !jmecmplx
    end do

    allocate(ket_box_start(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_confinement_pot_dbl_grid', &
           'ket_box_start',ierr)
    allocate(ket_start_in_box(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_confinement_pot_dbl_grid', &
           'ket_start_in_box',ierr)

    n_batches = ket_basis%max_on_proc / batch_size
    if (mod(ket_basis%max_on_proc,batch_size) > 0) n_batches = n_batches + 1

    local_start = 1
    do batch_count=1,n_batches

       if (pub_debug_on_root) write(stdout,'(a,2(i4,a))') &
            'DEBUG: Batch loop',batch_count, ' of',n_batches, &
            ' in integrals_confinement_pot_dbl_grid'

       local_end = min(local_start+batch_size-1,ket_basis%proc_num)

       max_current_size = local_end - local_start + 1
       call comms_reduce('MAX',max_current_size)

       call function_ops_batch_col_start(ket_box_start,ket_start_in_box, &
            batch_size,local_start,local_end,fftbox,cell,ket_basis)

       ! apply confinement operator on to kets
       call potential_confine_ngwf_batch(potential_fftbox_batch, &
            kets_on_grid, ket_basis, dbl_grid, fftbox, cell, &
            ket_start_in_box, ket_box_start, batch_size, local_start, &
            local_end, max_current_size)

       ! ndmh: calculate locpot integrals in SPAM3 format
       call timer_clock('integrals_confinement_pot_mat_els_batch',1)
       ongpu=.false.
       call function_ops_brappd_ketfftbox(locpot,              & ! inout
            bras_on_grid, bra_basis, cell, fftbox,                   & ! input
            potential_fftbox_batch, ket_box_start, batch_size,       & ! input
            local_start, local_end, idx_len, locpot_idx,             & ! input
            pub_locpot_scheme,ongpu)                                   ! input
       call timer_clock('integrals_confinement_pot_mat_els_batch',2)
       local_start = local_start + batch_size

    end do

    do loc_bra=1,bra_basis%num_on_proc(pub_my_proc_id)
       ibra = loc_bra + bra_basis%first_on_proc(pub_my_proc_id) - 1
       iat = bra_basis%atom_of_func(ibra)
       loc_iat = iat - row_par%first_atom_on_proc(pub_my_proc_id) + 1
       do jket=ket_basis%first_on_atom(iat), &
            ket_basis%first_on_atom(iat)+ket_basis%num_on_atom(iat)-1
           if(jket.ne.ibra) then
              ! agrecocmplx
              if (loc_cmplx) then
                 call sparse_put_element(cmplx_0,locpot,jket,ibra)
              else
                 call sparse_put_element(0.0_DP,locpot,jket,ibra)
              end if
           end if
       end do
    end do


    ! Deallocate workspace
    deallocate(ket_start_in_box,stat=ierr)
    call utils_dealloc_check('integrals_confinement_pot_dbl_grid', &
          'ket_start_in_box',ierr)
    deallocate(ket_box_start,stat=ierr)
    call utils_dealloc_check('integrals_confinement_pot_dbl_grid', &
          'ket_box_start',ierr)
    do batch_count =1, batch_size
        call data_fftbox_dealloc(potential_fftbox_batch(batch_count))
    end do
    deallocate(potential_fftbox_batch,stat=ierr)
    call utils_dealloc_check('integrals_confinement_pot_dbl_grid', &
          'potential_fftbox_batch',ierr)
    deallocate(locpot_idx,stat=ierr)
    call utils_dealloc_check('integrals_confinement_pot_dbl_grid', &
          'locpot_idx',ierr)

    ! ndmh: expand to full matrix
    select case (locpot_scheme)
    case ('LOWER','lower')
          call sparse_expand(locpot,PATTERN_LOWER)
    case ('ALTERNATE','alternate')
          call sparse_expand(locpot,PATTERN_ALTERNATE)
    case ('FULL','full')
          ! ndmh: only symmetrise if bra_basis==ket_basis
          if (bra_basis%name==ket_basis%name) then
             call sparse_scale(locpot,0.5_DP)
             call sparse_create(locpot_transpose, locpot)
             call sparse_transpose(locpot_transpose, locpot)
             call sparse_axpy(locpot, locpot_transpose, 1.0_DP)
             call sparse_destroy(locpot_transpose)
          end if
    case ('ASYM','asym')
       if (pub_on_root) write(stdout,'(a)') &
            'WARNING: Skipping symmetrisation of locpot in &
              &integrals_confinement_pot - '
       if (pub_on_root) write(stdout,'(a)') &
            'Resulting matrix may not be Hermitian.'
    case default
       call utils_abort('Error in integrals_confinement_pot: &
            &calculation pattern"'//&
            trim(pub_locpot_scheme)//'" not recognised')
    end select

    call comms_barrier

    ! Stop timer
    call timer_clock('integrals_confinement_pot_dbl_grid',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &integrals_confinement_pot_dbl_grid'

  end subroutine integrals_confinement_pot_dbl_grid



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine integrals_confinement_pot(locpot,  &                ! input-output
       bras_on_grid, bra_basis, kets_on_grid, ket_basis, &
       dbl_grid, cell, fftbox)

    !==========================================================================!
    ! This subroutine applies a local confining potential to each NGWF and     !
    ! calculates the braket:                                                   !
    ! confinement_pot_\alpha\beta = < bra_\alpha | ConfPot | ket_\beta >       !
    ! The non-diagonal terms are removed, as each NGWF only feels its own      !
    ! confinement                                                              !
    !==========================================================================!
    !  Arguments:                                                              !
    !    locpot     (inout)  : diagonal sparse local potential matrix elements !
    !    bras_on_grid  (in)  : bra functions on grid in PPD format             !
    !    bra_basis (input)   : The basis type for the bras                     !
    !    kets_on_grid  (in)  : ket functions on grid in PPD format             !
    !    ket_basis (input)   : The basis type for the kets                     !
    !==========================================================================!
    ! Derived from integrals_locpot by Gabriel Constantinescu in July 2015     !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_my_proc_id
    use constants, only: DP, stdout, cmplx_0
    use datatypes, only: FUNCTIONS, data_functions_alloc, &
         data_functions_dealloc, data_functions_copy
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketppd
    use parallel_strategy, only: PARAL_INFO
    use potential, only: potential_confine_ppd_funcs, &
         potential_input_to_workspace
    use rundat, only: pub_dbl_is_std
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3,sparse_put_element, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: locpot
    type(FUNC_BASIS), intent(in) :: bra_basis
    type(FUNCTIONS), intent(in) :: bras_on_grid
    type(FUNC_BASIS), intent(in) :: ket_basis
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FUNCTIONS), intent(in) :: kets_on_grid
    type(PARAL_INFO), pointer :: row_par ! rc2013: column parallel strategy

    ! Local Variables
    integer :: loc_bra,jket,ibra,iat,loc_iat
    type(FUNCTIONS) :: locpot_kets_on_grid
    logical :: loc_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &integrals_confinement_pot'

    loc_cmplx = kets_on_grid%iscmplx.or.bras_on_grid%iscmplx

    ! rc2013: obtain parallel strategy for this SPAM3
    call sparse_get_par(row_par, locpot, 'R')

    ! ndmh: if the potential is on the standard grid, use that as workspace
    if (pub_dbl_is_std) then

       ! Start timer
       call timer_clock('integrals_confinement_pot_coarse_grid',1)

       call data_functions_alloc(locpot_kets_on_grid, ket_basis%size_on_grid, &
                                  iscmplx=loc_cmplx)

       ! Copy in the ket functions
       call data_functions_copy(locpot_kets_on_grid, kets_on_grid)

       ! Apply the potential to the kets
       call potential_confine_ppd_funcs(locpot_kets_on_grid,ket_basis,&
            dbl_grid, cell, fftbox)

       ! Now do the overlap matrix of the bras with the locpot times the kets
       call function_ops_brappd_ketppd(locpot,bras_on_grid,bra_basis, &
            locpot_kets_on_grid,ket_basis,cell)

       do loc_bra=1,bra_basis%num_on_proc(pub_my_proc_id)
          ibra = loc_bra + bra_basis%first_on_proc(pub_my_proc_id) - 1
          iat = bra_basis%atom_of_func(ibra)
          loc_iat = iat - row_par%first_atom_on_proc(pub_my_proc_id) + 1
          do jket=ket_basis%first_on_atom(iat), &
               ket_basis%first_on_atom(iat)+ket_basis%num_on_atom(iat)-1
             if(jket.ne.ibra) then
                ! agrecocmplx
                if (loc_cmplx) then
                   call sparse_put_element(cmplx_0,locpot,jket,ibra)
                else
                   call sparse_put_element(0.0_DP,locpot,jket,ibra)
                end if
             end if
          end do
       end do

       ! Deallocate temporary storage
       call data_functions_dealloc(locpot_kets_on_grid)
       ! Stop timer
       call timer_clock('integrals_confinement_pot_coarse_grid',2)

    else

       ! Calculate the locpot integrals on the double grid, with interpolation
       call integrals_confinement_pot_dbl_grid(locpot,bras_on_grid,bra_basis, &
            kets_on_grid,ket_basis,dbl_grid,cell,fftbox)

    end if

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving  &
          &integrals_confinement_pot'

  end subroutine integrals_confinement_pot



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine integrals_locpot(locpot,  &                        ! input-output
       bras_on_grid, bra_basis, kets_on_grid, ket_basis, &
       fine_grid, dbl_grid, cell, fftbox, potential_fine)

    !==========================================================================!
    ! This subroutine calculates local potential integrals between two         !
    ! sets of functions.                                                       !
    ! Result is locpot_\alpha\beta = < bra_\alpha | Vloc | ket_\beta >         !
    !==========================================================================!
    !  Arguments:                                                              !
    !    locpot     (inout)  : sparse local potential matrix elements          !
    !    bras_on_grid  (in)  : bra functions on grid in PPD format             !
    !    bra_basis (input)   : The basis type for the bras                     !
    !    kets_on_grid  (in)  : ket functions on grid in PPD format             !
    !    ket_basis (input)   : The basis type for the kets                     !
    !    potential_fine (in) : local potential on fine sim cell grid           !
    !==========================================================================!
    ! Originaly written by Chris-Kriton Skylaris in January 2001,              !
    ! to use a "pair-box".                                                     !
    ! Improved by Oswaldo Dieguez in 2001 so that it "cshifts" the "pair-box". !
    ! Rewritten by Arash Mostofi in 2003 so that it uses a "triple-box".       !
    ! Rewritten by Chris-Kriton Skylaris on 23/11/2003 so that it runs on      !
    ! parallel computers.                                                      !
    ! Modifications and speed-ups by Peter Haynes 16/03/2005                   !
    ! Modified to use planned sums system by Nicholas Hine, May 2008           !
    ! Symmetrisation options added by Nicholas Hine, July 2008.                !
    ! Modified to remove integer buffers by Nicholas Hine, June 2009           !
    ! Adapted to use SPAM3 matrices by Nicholas Hine, June 2009                !
    ! Modified to take different function sets for bras and kets by            !
    ! Nicholas Hine on 11/11/2009.                                             !
    ! Modified to skip symmetrisation of the locpot matrix for the case of     !
    ! cross-Hamiltonians (bras/=kets) by Nicholas Hine, Oct 2010.              !
    ! Modified for possibility of performing integrals on coarse grid by       !
    ! Nicholas Hine in March 2011.                                             !
    ! Modified by Andrea Greco in May 2015 to allow use of complex NGWFs.      !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, data_functions_alloc, &
        data_functions_dealloc, data_functions_copy
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    ! agrecocmplx
    use function_ops, only: function_ops_brappd_ketppd
    ! agrecocmplx
    use potential, only: potential_apply_to_ppd_funcs, &
         potential_input_to_workspace
    use rundat, only: pub_dbl_is_std
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: locpot
    type(FUNC_BASIS), intent(in) :: bra_basis
    !real(kind=DP), intent(in) :: bras_on_grid(bra_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: bras_on_grid
    type(FUNC_BASIS), intent(in) :: ket_basis
    type(GRID_INFO), intent(in) :: fine_grid
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    !real(kind=DP), intent(in) :: kets_on_grid(ket_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: kets_on_grid
    real(kind=DP), intent(in) :: potential_fine(fine_grid%ld1,fine_grid%ld2, &
         fine_grid%max_slabs12)

    ! Local Variables
    integer :: ierr
    real(kind=DP), allocatable :: potential_work(:,:,:)
    !real(kind=DP), allocatable :: locpot_kets_on_grid(:)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS) :: locpot_kets_on_grid
    ! agrecocmplx: local variable for complex values
    logical :: loc_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering integrals_locpot'

    ! agrecocmplx
    loc_cmplx = kets_on_grid%iscmplx.or.bras_on_grid%iscmplx

    ! ndmh: if the potential is on the standard grid, use that as workspace
    if (pub_dbl_is_std) then

       ! Start timer
       call timer_clock('integrals_locpot_coarse_grid',1)

       ! Allocate coarse-grid workspace
       allocate(potential_work(dbl_grid%ld1,dbl_grid%ld2,&
            dbl_grid%max_group_slabs12),stat=ierr)
       call utils_alloc_check('integrals_locpot','potential_work',ierr)
       !allocate(locpot_kets_on_grid(ket_basis%size_on_grid),stat=ierr)
       !call utils_alloc_check('integrals_locpot','locpot_kets_on_grid',ierr)
       ! agrecocmplx: allocate using appropriate routine
       call data_functions_alloc(locpot_kets_on_grid, ket_basis%size_on_grid, &
                                  iscmplx=loc_cmplx)

       ! Get workspace potential from input potential
       call potential_input_to_workspace(potential_work,potential_fine, &
            dbl_grid,fine_grid)

       ! Copy in the ket functions
       ! agrecocmplx
       call data_functions_copy(locpot_kets_on_grid, kets_on_grid)

       ! Apply the potential to the kets
       call potential_apply_to_ppd_funcs(locpot_kets_on_grid,ket_basis,&
            potential_work, dbl_grid, cell, fftbox)

       ! Now do the overlap matrix of the bras with the locpot times the kets
       call function_ops_brappd_ketppd(locpot,bras_on_grid,bra_basis, &
            locpot_kets_on_grid,ket_basis,cell)

       ! Deallocate temporary storage
       !deallocate(locpot_kets_on_grid,stat=ierr)
       !call utils_dealloc_check('integrals_locpot','locpot_kets_on_grid',ierr)
       ! agrecocmplx: deallocate using appropriate routine
       call data_functions_dealloc(locpot_kets_on_grid)
       deallocate(potential_work,stat=ierr)
       call utils_dealloc_check('integrals_locpot','potential_work',ierr)

       ! Stop timer
       call timer_clock('integrals_locpot_coarse_grid',2)

    else

       ! Allocate double-grid workspace
       allocate(potential_work(dbl_grid%ld1,dbl_grid%ld2,&
            dbl_grid%max_group_slabs12),stat=ierr)
       call utils_alloc_check('integrals_locpot','potential_work',ierr)

       ! Get workspace potential from input potential
       call potential_input_to_workspace(potential_work,potential_fine, &
            dbl_grid,fine_grid)

       ! Calculate the locpot integrals on the double grid, with interpolation
       call integrals_locpot_dbl_grid(locpot,bras_on_grid,bra_basis, &
            kets_on_grid,ket_basis,potential_work,dbl_grid,cell,fftbox)

       ! Deallocate temporary storage
       deallocate(potential_work,stat=ierr)
       call utils_dealloc_check('integrals_locpot','potential_work',ierr)

    end if

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving integrals_locpot'

  end subroutine integrals_locpot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine integrals_div_locpot_grad_functions(mat, &       ! input-output
       bras_on_grid, bra_basis, kets_on_grid, ket_basis, &
       fine_grid, dbl_grid, cell, fftbox, potential_fine)
    !==========================================================================!
    ! This subroutine calculates integrals of the form                         !
    !   int_\alpha\beta = < bra_\alpha | \nabla . ( Vloc \nabla ket_\beta ) >  !
    ! This is similar to the kinetic energy integrals, but with divergence     !
    ! operator applied to the product of the local potential and ket functions !
    ! rather than the ket functions alone.                                     !
    !==========================================================================!
    !  Arguments:                                                              !
    ! TODO Add details for missing arguments
    !    mat      (inout)    : sparse matrix of integrals                      !
    !    bras_on_grid  (in)  : bra functions on grid in PPD format             !
    !    bra_basis (input)   : The basis type for the bras                     !
    !    kets_on_grid  (in)  : ket functions on grid in PPD format             !
    !    ket_basis (in)      : The basis type for the kets                     !
    !    potential_fine (in) : local potential on fine sim cell grid           !
    !==========================================================================!
    ! integrals_div_locpot_grad_functions created by James C. Womack, 2015.    !
    !==========================================================================!
    ! Based on integrals_kinetic and integrals_locpot_dbl_grid (code copied    !
    ! and modified by James C. Womack)                                         !
    !==========================================================================!
    use cell_grid, only: GRID_INFO
    use comms, only: comms_barrier, comms_reduce
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_batch_col_start, &
        function_ops_brappd_ketfftbox
    use kinetic, only: kinetic_apply_to_func_batch, kinetic_grad_to_func_batch,&
         kinetic_grad_on_cart_box
    use potential, only: potential_apply_to_ngwf_batch, &
         potential_input_to_workspace, potential_apply_to_fftbox_batch
    use rundat, only: pub_fftbox_batch_size, pub_dbl_is_std, &
         pub_threads_fftbox
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_generate_index, sparse_index_length, &
         sparse_scale, sparse_create, sparse_transpose, sparse_axpy, &
         sparse_destroy, sparse_expand, pattern_lower, pattern_alternate
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat
    type(FUNC_BASIS), intent(in) :: bra_basis
    type(FUNCTIONS), intent(in) :: bras_on_grid
    type(FUNC_BASIS), intent(in) :: ket_basis
    type(FUNCTIONS), intent(in) :: kets_on_grid
    type(GRID_INFO), intent(in) :: fine_grid
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(in) :: potential_fine(fine_grid%ld1,fine_grid%ld2, &
         fine_grid%max_slabs12)

    ! Local Variables
    type(FFTBOX_DATA), allocatable, dimension(:) :: ket_fftbox_batch
    type(FFTBOX_DATA), allocatable, dimension(:,:,:) :: grad_work_fftbox_batch
    complex(kind=DP), allocatable, dimension(:,:,:,:,:) :: zwork_grad_box
    integer :: batch_size
    integer :: batch_count
    integer :: batch_number
    integer :: n_batches
    integer :: local_start, local_end
    integer :: max_current_size ! maximum batch size ovar all procs
    integer :: idx_len  ! pdh: sparse matrix index length
    integer :: ierr ! pdh: error flag
    integer :: iii
    integer :: row1
    integer, allocatable, dimension(:) :: mat_idx ! pdh: sparse index
    integer, allocatable :: ket_box_start(:,:)
    integer, allocatable :: ket_start_in_box(:,:)
    !type(SPAM3) :: mat_transpose  ! ndmh: for symmetrising mat <-- unused
    character(len=80) :: mat_scheme
    ! agrecocmplx
    logical :: loc_cmplx, ongpu
    integer :: dim
    real(kind=DP), allocatable :: potential_dbl(:,:,:)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering integrals_div_locpot_grad_functions'


    ! Check for unsupported conditions
    if (pub_dbl_is_std) then
       ! JCW: dbl_grid_scale must be 2.0
       call utils_abort("Error in integrals_div_locpot_grad_functions: &
         &using dbl_grid_scale == 1.0 is unsupported.")
    end if

    ! Start timer
    call timer_clock('integrals_div_locpot_grad_functions',1)

    ! Obtain index for matrix
    idx_len = sparse_index_length(mat)
    allocate(mat_idx(idx_len),stat=ierr)
    call utils_alloc_check('integrals_div_locpot_grad_functions','mat_idx',ierr)
    call sparse_generate_index(mat_idx,mat)

    ! Define shorthand variables
    batch_size = pub_fftbox_batch_size

    ! agrecocmplx
    loc_cmplx = kets_on_grid%iscmplx

    ! TODO: Check that these mat_scheme values are appropriate for
    !       dfdtau integrals on a theoretical basis
    ! lr408: Force calculation of the full integral if this is for a cross matrix
    ! rc2013: check that the bases refer to the same region
    if (bra_basis%name/=ket_basis%name .or. bra_basis%ireg/=ket_basis%ireg) then
       mat_scheme = 'FULL'
    else
       mat_scheme = 'ALTERNATE'
    end if

    ! Allocate workspace
    ! grad_work_fftbox_batch has 3 Cartesian components and two indepent
    ! sets of these. Set 1 will initially be used to store gradients of NGWFs
    ! and set 2 will initially be used to store the gradients of NGWFs with
    ! local potential applied.
    allocate(grad_work_fftbox_batch(batch_size,3,2),stat=ierr)
    call utils_alloc_check('integrals_div_locpot_grad_functions',&
         'grad_work_fftbox_batch',ierr)
    ! zwork_grad_box is a complex array, with the dimensions of an FFTbox, and
    ! with 3 Cartesian components
    allocate(zwork_grad_box(fftbox%total_ld1,fftbox%total_ld2, &
         fftbox%total_pt3,3,2))
    call utils_alloc_check('integrals_div_locpot_grad_functions',&
         'zwork_grad_box',ierr)
    allocate(ket_fftbox_batch(batch_size), stat=ierr)
    call utils_alloc_check('integrals_div_locpot_grad_functions',&
         'ket_fftbox_batch',ierr)
    do batch_count = 1, batch_size
        ! JCW: Allocate fftboxes in FFTBOX_DATA objects
        do dim = 1, 3
           call data_fftbox_alloc(grad_work_fftbox_batch(batch_count,dim,1), &
                fftbox%total_ld1, fftbox%total_ld2, fftbox%total_pt3, &
                iscmplx=loc_cmplx)
           call data_fftbox_alloc(grad_work_fftbox_batch(batch_count,dim,2), &
                fftbox%total_ld1, fftbox%total_ld2, fftbox%total_pt3, &
                iscmplx=loc_cmplx)
        end do
        call data_fftbox_alloc(ket_fftbox_batch(batch_count), &
             fftbox%total_ld1, fftbox%total_ld2, fftbox%total_pt3, &
             iscmplx=loc_cmplx)
    end do
    allocate(ket_box_start(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_div_locpot_grad_functions',&
         'ket_box_start',ierr)
    allocate(ket_start_in_box(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_div_locpot_grad_functions',&
         'ket_start_in_box',ierr)

    ! Allocate double-grid workspace
    allocate(potential_dbl(dbl_grid%ld1,dbl_grid%ld2,&
         dbl_grid%max_group_slabs12),stat=ierr)
    call utils_alloc_check('integrals_div_locpot_grad_functions',&
         'potential_dbl',ierr)

    ! Get workspace potential from input potential
    call potential_input_to_workspace(potential_dbl,potential_fine, &
         dbl_grid,fine_grid)

    ! cks: number of row-steps per row-block
    n_batches = ket_basis%max_on_proc / batch_size
    if (mod(ket_basis%max_on_proc,batch_size) > 0) n_batches = n_batches + 1

    local_start = 1
    do batch_number=1,n_batches

       if (pub_debug_on_root) write(stdout,'(a,2(i4,a))') 'DEBUG: Batch loop ',&
            batch_number, ' of ',n_batches,' in integrals_div_locpot_grad_functions'

       local_end = min(local_start+batch_size-1,ket_basis%proc_num)

       ! cks: maximum size of current batch over all procs
       max_current_size = local_end - local_start + 1
       call comms_reduce('MAX',max_current_size)

       call function_ops_batch_col_start(ket_box_start,ket_start_in_box, &
            batch_size,local_start,local_end,fftbox,cell,ket_basis)

       ! JCW: Apply grad operator to kets and place result into
       ! JCW: grad_work_fftbox_batch(:,:,1)
       call kinetic_grad_to_func_batch(grad_work_fftbox_batch(:,:,1), &   ! out
            kets_on_grid, ket_basis, cell, fftbox, ket_start_in_box, &
            batch_size, local_start, local_end) ! in

       ! JCW: Apply potential to gradient of kets and place result into
       ! JCW: grad_work_fftbox_batch(:,:,2)
       do dim = 1, 3
          call potential_apply_to_fftbox_batch(grad_work_fftbox_batch(:,dim,2), &
               grad_work_fftbox_batch(:,dim,1), kets_on_grid%iscmplx, &
               potential_dbl, dbl_grid, fftbox, &
               ket_start_in_box, ket_box_start, batch_size, local_start, &
               local_end, max_current_size)
       end do

       ! JCW: Apply divergence operator to gradient of kets with local potential
       ! JCW: applied and place result back into grad_work_fftbox_batch(:,1,1).

       ! JCW: Get Cartesian components of divergence and place into
       ! JCW: grad_work_fftbox_batch(:,:,1)

       ! TODO: Investigate simulataneous Fourier-transform of two FFTboxes
       ! JCW: Fourier transform FFTboxes containing potential applied to
       ! JCW: gradients of NGWFs.
       do iii = local_start, local_end, 1
          row1 = iii
          batch_count = iii - local_start + 1
          ! JCW: Copy real space FFTboxes to complex array and Fourier
          ! JCW: transform to reciprocal space
          do dim = 1, 3
             zwork_grad_box(:,:,:,dim,1) = cmplx(&
                  grad_work_fftbox_batch(batch_count,dim,2)%d(:,:,:),&
                  0.0_DP,kind=DP)
             call fourier_apply_box('Coarse','Forward',&
                  zwork_grad_box(:,:,:,dim,1),omp=pub_threads_fftbox)
          end do

          ! JCW: Apply gradient operator in reciprocal space
          call kinetic_grad_on_cart_box(zwork_grad_box(:,:,:,:,1),&
               zwork_grad_box(:,:,:,:,2),fftbox)

          ! JCW: Zero ket_fftbox_batch for summing into
          ket_fftbox_batch(batch_count)%d(:,:,:) = 0.0_DP
          do dim=1,3
             ! JCW: Reverse transform FFTboxes to real space following
             ! JCW: application of gradient operator
             call fourier_apply_box('Coarse','Backward',&
                  zwork_grad_box(:,:,:,dim,2), omp=pub_threads_fftbox)
             ! JCW: Sum real components of zwork_grad_box(:,:,:,:,2) into
             ! JCW: ket_fftbox_batch, multiply by -0.5_DP
             ket_fftbox_batch(batch_count)%d(:,:,:) = &
                  ket_fftbox_batch(batch_count)%d(:,:,:) - &
                  0.5_DP * real(zwork_grad_box(:,:,:,dim,2),kind=DP)
          end do
       end do

       if (pub_debug_on_root) write(stdout,'(a,2(a))') 'DEBUG:  ', &
            'after kinetic_apply_to_func_batch'

       ! ndmh: calculate kinetic energy integrals in SPAM3 format
       ! JCW: ongpu is always .false., since integrals_div_locpot_grad_functions
       ! JCW: does not currently support GPUs.
       ongpu=.false.
       call function_ops_brappd_ketfftbox(mat, &
            bras_on_grid, bra_basis, cell, fftbox, &
            ket_fftbox_batch, ket_box_start, batch_size, &
            local_start, local_end, idx_len, mat_idx, mat_scheme,ongpu)

       local_start = local_start + batch_size

    end do

    ! ndmh: expand to full matrix
    if (mat_scheme=='ALTERNATE') then
       call sparse_expand(mat,PATTERN_ALTERNATE)
    end if

    ! deallocate double-grid workspace
    deallocate(potential_dbl,stat=ierr)
    call utils_dealloc_check('integrals_div_locpot_grad_functions',&
         'potential_dbl',ierr)

    ! pdh: deallocate workspace
    deallocate(ket_start_in_box,stat=ierr)
    call utils_dealloc_check('integrals_div_locpot_grad_functions',&
         'ket_start_in_box',ierr)
    deallocate(ket_box_start,stat=ierr)
    call utils_dealloc_check('integrals_div_locpot_grad_functions',&
         'ket_box_start',ierr)
    ! agrecocmplx: deallocate using appropriate routine
    do batch_count = 1, batch_size
        call data_fftbox_dealloc(ket_fftbox_batch(batch_count))
    end do
    deallocate(ket_fftbox_batch,stat=ierr)
    call utils_dealloc_check('integrals_div_locpot_grad_functions',&
         'ket_fftbox_batch',ierr)
    deallocate(zwork_grad_box,stat=ierr)
    call utils_dealloc_check('integrals_div_locpot_grad_functions',&
         'zwork_grad_box',ierr)
    do batch_count = 1, batch_size
        do dim = 1, 3
           call data_fftbox_dealloc(grad_work_fftbox_batch(batch_count,dim,1))
           call data_fftbox_dealloc(grad_work_fftbox_batch(batch_count,dim,2))
        end do
    end do
    deallocate(grad_work_fftbox_batch,stat=ierr)
    call utils_dealloc_check('integrals_div_locpot_grad_functions',&
         'grad_work_fftbox_batch',ierr)
    deallocate(mat_idx,stat=ierr)
    call utils_dealloc_check('integrals_div_locpot_grad_functions','mat_idx',ierr)

    ! pdh: sync procs
    call comms_barrier

    ! Stop timer
    call timer_clock('integrals_div_locpot_grad_functions',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &integrals_div_locpot_grad_functions'

  end subroutine integrals_div_locpot_grad_functions

  subroutine test_integrals_div_locpot_grad_functions(mat, &
       bras_on_grid, bra_basis, kets_on_grid, ket_basis, &
       fine_grid, dbl_grid, cell, fftbox, potential_fine,denskern)
    !==========================================================================!
    ! This subroutine is designed to test integrals_div_locpot_grad_functions. !
    ! This is done through a sequence of test cases                            !
    ! TODO Provide more detailed documentation                                 !
    !==========================================================================!
    use comms, only: pub_on_root, &
         comms_reduce, comms_recv, comms_send, comms_wait
    use cell_grid, only: GRID_INFO
    use constants, only: DP, stdout, file_maxsize
    use datatypes, only: FUNCTIONS, data_functions_alloc, data_functions_dealloc
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use kinetic, only: kinetic_grad_to_funcs_on_grid
    use rundat, only: pub_num_spins, pub_num_kpoints, PUB_1K, pub_rootname
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_axpy, sparse_create, sparse_destroy, &
         sparse_get_col, sparse_first_elem_on_proc, &
         sparse_scale, sparse_trace, sparse_transpose, sparse_transpose_structure
    use sparse_embed, only: SPAM3_EMBED_ARRAY
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, &
         utils_assert, utils_unit, utils_open_unit_check, utils_close_unit_check

    implicit none

    ! Arguments
    ! Previously evaluated SPAM3 matrices of div_locpot_grad integrals
    type(SPAM3), intent(in) :: mat
    type(FUNC_BASIS), intent(in) :: bra_basis
    type(FUNCTIONS), intent(in) :: bras_on_grid
    type(FUNC_BASIS), intent(in) :: ket_basis
    type(FUNCTIONS), intent(in) :: kets_on_grid
    type(GRID_INFO), intent(in) :: fine_grid
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(in) :: potential_fine(fine_grid%ld1,fine_grid%ld2, &
         fine_grid%max_slabs12,pub_num_spins)
    ! rc2013: EMBED_FIX! This assumes 1 subsystem
    type(SPAM3_EMBED_ARRAY), intent(in)   :: denskern

    ! Local variables
    integer :: is
    integer :: dim
    ! mat_local_1 : div_locpot_grad using potential_fine passed into routine
    ! mat_local_2 : div_locpot_grad using potential_fine = 1.0_DP
    type(SPAM3) :: mat_local_1(pub_num_spins), mat_local_1_trans(pub_num_spins)
    type(SPAM3) :: mat_local_2
    type(SPAM3) :: locpot_local(pub_num_spins), locpot_sum(pub_num_spins)
    type(SPAM3) :: kinet_local, kinet_local_trans
    type(SPAM3) :: mat_tmp
    character(len=30) :: struc_trans
    real(kind=DP), allocatable :: potential_local(:,:,:,:)
    type(FUNCTIONS) :: bra_grads_on_grid(3), ket_grads_on_grid(3)
    real(kind=DP) :: max_diff_arr(3) ! 1: max_diff, 2: max_val1, 3: max_val2
    real(kind=DP) :: trace_diff_arr(3) ! 1: abs_diff, 2: trace1, 3: trace2
    real(kind=DP) :: trace1, trace2
    integer :: ierr, iunit
    character(len=file_maxsize) :: output_file
    logical,save :: file_created = .false.
    integer,save :: itest = 0


    ! Check for unsupported situations
    call utils_assert(bras_on_grid%iscmplx.eqv.kets_on_grid%iscmplx,&
         "Error in test_integrals_div_locpot_grad_functions: &
         &bra and ket must be both real or both complex.")
    call utils_assert(.not.bras_on_grid%iscmplx,&
         "Error in test_integrals_div_locpot_grad_functions: &
         &Complex NGWFs not supported.")
    call utils_assert(pub_num_kpoints == PUB_1K, &
         "Error in test_integrals_div_locpot_grad_functions: &
         &More than one k-point not supported.")

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering test_integrals_div_locpot_grad_functions'

    ! Start timer
    call timer_clock('test_integrals_div_locpot_grad_functions',1)

    ! rc2013: allocate space for array
    allocate(potential_local(fine_grid%ld1, fine_grid%ld2, &
         fine_grid%max_slabs12, pub_num_spins), stat=ierr)
    call utils_alloc_check('test_integrals_div_locpot_grad_functions', &
         'potential_local',ierr)

    ! Open unit for writing out test data (rather than using output file)
    if (pub_on_root) then
       iunit = utils_unit()
       output_file = trim(pub_rootname)//"_integral_test.out"
       if (pub_debug_on_root) then
          write(stdout,'(a,a,a,i5)') "DEBUG: Writing test output to file ", &
               trim(output_file),", unit: ", iunit
       end if
       if (.not.file_created) then
          ! Start by overwriting any pre-existing file of the same name
          open(unit=iunit,file=output_file,action="write",status="replace",&
               iostat=ierr)
          file_created = .true.
          ! Write header to file
          write(iunit,*) "JCW test code: integrals_div_locpot_grad_functions"
          ! Write intro to test
          write(iunit,*)
          write(iunit,'(a,i5)') "Test sequence ",itest
          itest = itest + 1
       else
          ! File already exists, so simply append to this
          open(unit=iunit,file=output_file,action="write",status="old",&
               position="append",iostat=ierr)
          ! Write intro for new test
          write(iunit,*)
          write(iunit,'(a,i5)') "Test sequence ",itest
          itest = itest + 1
       end if
       call utils_open_unit_check('test_integrals_div_locpot_grad_functions',&
            output_file,ierr)
    else
      iunit = -100 ! Set to nonsensical value -- should not be used on non-root
                   ! proc
    end if

    ! Set up local SPAM3 objects with sparsity of input mat
    call sparse_create(kinet_local,mat,iscmplx=bras_on_grid%iscmplx)
    call sparse_create(mat_local_2,mat,iscmplx=bras_on_grid%iscmplx)
    do is=1,pub_num_spins
       call sparse_create(mat_local_1(is),mat,iscmplx=bras_on_grid%iscmplx)
    end do

    ! Set up local SPAM3 objects with sparsity of transpose of input mat
    call sparse_transpose_structure(struc_trans,mat)
    mat_tmp%structure=struc_trans
    call sparse_create(mat_tmp)
    call sparse_create(kinet_local_trans,mat_tmp,iscmplx=bras_on_grid%iscmplx)
    do is=1,pub_num_spins
       call sparse_create(mat_local_1_trans(is),mat_tmp,iscmplx=bras_on_grid%iscmplx)
    end do
    call sparse_destroy(mat_tmp)

    ! 1. Check that matrix of div_locpot_grad integrals is hermitian
    if (pub_on_root) then
       write(iunit,*) "Test 1: check if matrix is hermitian"
    end if
    ! Evaluate kinetic energy integrals
    call integrals_kinetic(kinet_local, bras_on_grid, bra_basis, &
         kets_on_grid, ket_basis, cell, fftbox)
    ! Evaluate div_locpot_grad integrals
    do is = 1, pub_num_spins
       call integrals_div_locpot_grad_functions(mat_local_1(is), &
            bras_on_grid,bra_basis,kets_on_grid,ket_basis, &
            fine_grid, dbl_grid, cell, fftbox, &
            potential_fine(:,:,:,is))
    end do
    ! Evaluate matrix transposes
    call sparse_transpose(kinet_local_trans,kinet_local)
    do is = 1, pub_num_spins
       call sparse_transpose(mat_local_1_trans(is),mat_local_1(is))
    end do

    ! Compare integrals with transposed integrals
    ! kinetic energy integrals
    ! Element-by-element of comparison the SPAM3 objects
    max_diff_arr = 0.0_DP
    if (pub_on_root) write(iunit,*) "kinetic energy integrals..."
    call element_by_element_comparison(kinet_local,kinet_local_trans,&
         bra_basis,max_diff_arr)
    if (pub_on_root) then
       write(iunit,'(a,es24.16)') "Maximum per-integral difference = ", &
            max_diff_arr(1)
       write(iunit,'(a,es24.16,a)') "(kinetic integral value         = ", &
            max_diff_arr(2),")"
       write(iunit,'(a,es24.16,a)') "(transposed integral value      = ", &
            max_diff_arr(3),")"
    end if
    ! Compare of trace of SPAM3 matrices with density kernel
    trace_diff_arr = 0.0_DP
    do is =1 , pub_num_spins
       call trace_comparison(kinet_local,kinet_local_trans,&
            denskern%m(is,PUB_1K)%p,trace1,trace2)
          trace_diff_arr(2) = trace_diff_arr(2) + trace1
          trace_diff_arr(3) = trace_diff_arr(3) + trace2
    end do
    trace_diff_arr(1) = abs( trace_diff_arr(2) - trace_diff_arr(3) )
    if (pub_on_root) then
       write(iunit,'(a,es24.16)') "Abs. diff. in trace with K      = ", &
            trace_diff_arr(1)
       write(iunit,'(a,es24.16,a)') "(kinetic integral trace         = ", &
            trace_diff_arr(2),")"
       write(iunit,'(a,es24.16,a)') "(transposed integral trace      = ", &
            trace_diff_arr(3),")"
    end if

    ! div_locpot_grad integrals
    ! Element-by-element of comparison the SPAM3 objects
    max_diff_arr = 0.0_DP
    if (pub_on_root) write(iunit,*) "div_locpot_grad integrals..."
    do is=1,pub_num_spins
      call element_by_element_comparison(mat_local_1(is),mat_local_1_trans(is),&
           bra_basis,max_diff_arr)
      if (pub_on_root) then
         write(iunit,'(a,es24.16)') "Maximum per-integral difference = ", &
              max_diff_arr(1)
         write(iunit,'(a,es24.16,a)') "(div_locpot_grad integral value = ", &
              max_diff_arr(2),")"
         write(iunit,'(a,es24.16,a)') "(transposed integral value      = ", &
              max_diff_arr(3),")"
      end if
    end do
    ! Compare of trace of SPAM3 matrices with density kernel
    trace_diff_arr = 0.0_DP
    do is = 1, pub_num_spins
       call trace_comparison(mat_local_1(is),mat_local_1_trans(is),&
            denskern%m(is,PUB_1K)%p,trace1,trace2)
       trace_diff_arr(2) = trace_diff_arr(2) + trace1
       trace_diff_arr(3) = trace_diff_arr(3) + trace2
    end do
    trace_diff_arr(1) = abs( trace_diff_arr(2) - trace_diff_arr(3) )
    if (pub_on_root) then
       write(iunit,'(a,es24.16)') "Abs. diff. in trace with K      = ", &
            trace_diff_arr(1)
       write(iunit,'(a,es24.16,a)') "(div_locpot_grad integral trace = ", &
            trace_diff_arr(2),")"
       write(iunit,'(a,es24.16,a)') "(transposed integral trace      = ", &
            trace_diff_arr(3),")"
    end if


    ! 2. Evaluate div_locpot_grad integrals with potential = 1.0_DP.
    !    Compare to the kinetic energy integrals. These should be the
    !    equivalent.
    !    (No need to consider spin, since potential_local is equivalent for
    !    all spin components)

    if (pub_on_root) then
       write(iunit,*) "Test 2: compare with kinetic energy integrals using &
            &potential_local = 1.0_DP"
    end if
    potential_local(:,:,:,:) = 1.0_DP
    ! Evaluate new div_locpot_grad integrals with potential_local = 1.0_DP
    call integrals_div_locpot_grad_functions(mat_local_2, &
         bras_on_grid,bra_basis,kets_on_grid,ket_basis, &
         fine_grid, dbl_grid, cell, fftbox, &
         potential_local(:,:,:,1))
    ! (Use kinetic energy integrals from Test 1)

    ! Compare SPAM3 objects, element-by-element
    max_diff_arr = 0.0_DP
    call element_by_element_comparison(mat_local_2,kinet_local,bra_basis,&
         max_diff_arr)
    if (pub_on_root) then
       write(iunit,'(a,es24.16)') "Maximum per-integral difference = ", &
            max_diff_arr(1)
       write(iunit,'(a,es24.16,a)') "(div_locpot_grad integral value = ", &
            max_diff_arr(2),")"
       write(iunit,'(a,es24.16,a)') "(kinetic integral value         = ", &
            max_diff_arr(3),")"
    end if

    ! Compare of trace of SPAM3 matrices with density kernel
    trace_diff_arr = 0.0_DP
    do is = 1, pub_num_spins
       call trace_comparison(mat_local_2,kinet_local,denskern%m(is,PUB_1K)%p,&
            trace1,trace2)
       trace_diff_arr(2) = trace_diff_arr(2) + trace1
       trace_diff_arr(3) = trace_diff_arr(3) + trace2
    end do
    trace_diff_arr(1) = abs( trace_diff_arr(2) - trace_diff_arr(3) )
    if (pub_on_root) then
       write(iunit,'(a,es24.16)') "Abs. diff. in trace with K      = ", &
            trace_diff_arr(1)
       write(iunit,'(a,es24.16,a)') "(div_locpot_grad integral trace = ", &
            trace_diff_arr(2),")"
       write(iunit,'(a,es24.16,a)') "(kinetic integral trace         = ", &
            trace_diff_arr(3),")"
    end if

    if (pub_on_root) then
       write(iunit,*) "Test 2a: compare with kinetic energy integrals using &
            &potential_fine"
    end if
    ! 2a. Compare div_locpot_grad and kinetic integrals to ensure that the
    ! kinetic and div_locpot_grad integrals do not also agree when the
    ! external potential_fine is used to evaluate div_locpot_grad integrals
    ! Compare SPAM3 objects, element-by-element
    max_diff_arr = 0.0_DP
    do is = 1, pub_num_spins
       call element_by_element_comparison(mat_local_1(is),kinet_local,&
            bra_basis,max_diff_arr)
       if (pub_on_root) then
          write(iunit,'(a,es24.16,a,i3)') "Maximum per-integral difference = ", &
               max_diff_arr(1), ", is = ", is

          write(iunit,'(a,es24.16,a,i3,a)') "(div_locpot_grad integral value = ", &
               max_diff_arr(2),", is = ", is,")"
          write(iunit,'(a,es24.16,a,i3,a)') "(kinetic integral value         = ", &
               max_diff_arr(3),", is = ", is,")"
       end if
    end do
    ! Compare of trace of SPAM3 matrices with density kernel
    trace_diff_arr = 0.0_DP
    do is = 1, pub_num_spins
       call trace_comparison(mat_local_1(is),kinet_local,denskern%m(is,PUB_1K)%p,&
            trace1,trace2)
       trace_diff_arr(2) = trace_diff_arr(2) + trace1
       trace_diff_arr(3) = trace_diff_arr(3) + trace2
    end do
    trace_diff_arr(1) = abs( trace_diff_arr(2) - trace_diff_arr(3) )
    if (pub_on_root) then
       write(iunit,'(a,es24.16)') "Abs. diff. in trace with K      = ", &
            trace_diff_arr(1)
       write(iunit,'(a,es24.16,a)') "(div_locpot_grad integral trace = ", &
            trace_diff_arr(2),")"
       write(iunit,'(a,es24.16,a)') "(kinetic integral trace         = ", &
            trace_diff_arr(3),")"
    end if

    ! 3. Evaluate locpot integrals with gradients of NGWFs on grid passed
    !    as arguments (rather than the NGWFs themselves).
    !    Compare to the matrix of div_locpot_grad integrals in mat.
    !    The results should be close, though some difference is expected
    !    because the gradients of NGWFs on grid will be "shaved" to the
    !    localization spheres of the NGWFs.

    if (pub_on_root) then
       write(iunit,*) "Test 3: compare with locpot integrals"
    end if

    ! Allocate attributes of FUNCTIONS objects for storing gradients of NGWFs
    do dim = 1, 3
       call data_functions_alloc(bra_grads_on_grid(dim),&
            bra_basis%size_on_grid,bras_on_grid%iscmplx)
       call data_functions_alloc(ket_grads_on_grid(dim),&
            ket_basis%size_on_grid,kets_on_grid%iscmplx)
    end do

    ! Create SPAM3 object to hold result of integrals_locpot
    do is=1,pub_num_spins
       call sparse_create(locpot_local(is),mat,iscmplx=bras_on_grid%iscmplx)
    end do
    ! Create SPAM3 object for result of integrals_locpot to be summed into
    do is=1,pub_num_spins
       call sparse_create(locpot_sum(is),mat,iscmplx=bras_on_grid%iscmplx)
    end do


    ! Take gradients of bras_on_grid and kets_on_grid in Fourier space and
    ! place on grid (PPDs)
    call kinetic_grad_to_funcs_on_grid(bra_grads_on_grid,bras_on_grid,&
         bra_basis,cell,fftbox)
    call kinetic_grad_to_funcs_on_grid(ket_grads_on_grid,kets_on_grid,&
         ket_basis,cell,fftbox)


    ! Evaluate locpot integrals
    do is=1,pub_num_spins
       ! Ensure locpot_sum set to zero before summing into
       call sparse_scale(locpot_sum(is),0.0_DP)
       do dim = 1, 3
          call integrals_locpot(locpot_local(is), &
               bra_grads_on_grid(dim),bra_basis, &
               ket_grads_on_grid(dim),ket_basis, &
               fine_grid, dbl_grid, cell, fftbox, &
               potential_fine(:,:,:,is))
          call sparse_axpy(locpot_sum(is),locpot_local(is),0.5_DP)
       end do
       ! (Use div_locpot_grad integrals from Test 1)
    end do

    ! Compare div_locpot_grad integrals with locpot integrals
    ! calculated locally, element-by-element of the SPAM3 objects
    max_diff_arr = 0.0_DP
    do is=1,pub_num_spins
      call element_by_element_comparison(mat_local_1(is),locpot_sum(is),bra_basis,&
           max_diff_arr)
      if (pub_on_root) then
         write(iunit,'(a,es24.16,a,i3)') "Maximum per-integral difference = ", &
              max_diff_arr(1),", is = ", is
         write(iunit,'(a,es24.16,a,i3,a)') "(div_locpot_grad integral value = ", &
              max_diff_arr(2),", is = ", is,")"
         write(iunit,'(a,es24.16,a,i3,a)') "(locpot integral value          = ", &
              max_diff_arr(3),", is = ", is,")"
      end if
    end do

    ! Compare of trace of SPAM3 matrices with density kernel
    trace_diff_arr = 0.0_DP
    do is = 1, pub_num_spins
       call trace_comparison(mat_local_1(is),locpot_sum(is),denskern%m(is,PUB_1K)%p,&
            trace1,trace2)
       trace_diff_arr(2) = trace_diff_arr(2) + trace1
       trace_diff_arr(3) = trace_diff_arr(3) + trace2
    end do
    trace_diff_arr(1) = abs( trace_diff_arr(2) - trace_diff_arr(3) )
    if (pub_on_root) then
       write(iunit,'(a,es24.16)') "Abs. diff. in trace with K      = ", &
            trace_diff_arr(1)
       write(iunit,'(a,es24.16,a)') "(div_locpot_grad integral trace = ", &
            trace_diff_arr(2),")"
       write(iunit,'(a,es24.16,a)') "(locpot integral trace          = ", &
            trace_diff_arr(3),")"
    end if

    ! rc2013: deallocate space for array
    deallocate(potential_local, stat=ierr)
    call utils_dealloc_check('test_integrals_div_locpot_grad_functions', &
         'potential_local',ierr)

    ! Destroy mat_local_1 SPAM3 object
    call sparse_destroy(kinet_local_trans)
    call sparse_destroy(kinet_local)
    call sparse_destroy(mat_local_2)
    do is=1,pub_num_spins
       call sparse_destroy(mat_local_1_trans(is))
       call sparse_destroy(mat_local_1(is))
       call sparse_destroy(locpot_local(is))
       call sparse_destroy(locpot_sum(is))
    end do

    ! Deallocate attributes of FUNCTIONS objects for storing gradients of NGWFs
    do dim = 1, 3
       call data_functions_dealloc(bra_grads_on_grid(dim))
    end do

    ! Close
    if (pub_on_root) then
       if (pub_debug_on_root) then
          write(stdout,'(a,a,a,i5)') "DEBUG: Closing file ", &
               trim(output_file),", unit: ", iunit
       end if
       close(iunit,iostat=ierr)
       call utils_close_unit_check('test_integrals_div_locpot_grad_functions',&
            output_file,ierr)
    end if
    ! Stop timer
    call timer_clock('test_integrals_div_locpot_grad_functions',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving test_integrals_div_locpot_grad_functions'


    contains
      subroutine element_by_element_comparison(mat1,mat2,fbasis,&
           max_diff_arr)
        !========================================================================!
        ! Compares elements in SPAM3 matrices mat1 and mat2 and returns the      !
        ! maximum per-element difference, and values in mat1 and mat2            !
        ! corresponding to this difference in max_diff_arr.                      !
        !========================================================================!
        use comms, only: pub_my_proc_id, pub_on_root, pub_total_num_procs, &
             pub_root_proc_id, &
             comms_recv, comms_send, comms_wait
        use constants, only: DP, stdout
        use datatypes, only: FUNCTIONS
        use function_basis, only: FUNC_BASIS
        use sparse, only: SPAM3, sparse_num_rows, &
             sparse_get_col, sparse_num_elems_on_proc, &
             sparse_first_elem_on_proc, sparse_proc_num_element
        use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

        implicit none
        ! Arguments
        type(SPAM3),intent(in) :: mat1, mat2
        type(FUNC_BASIS), intent(in) :: fbasis
        real(kind=DP),intent(out) :: max_diff_arr(3) ! 1: max_diff, 2: max_val1, 3: max_val2

        ! Local variables
        integer :: ierr
        integer :: iproc
        integer :: num_cols, num_rows
        integer :: first_col
        integer :: local_col, local_row
        real(kind=DP) :: val1, val2
        real(kind=DP), allocatable :: col_data1(:), col_data2(:)
        real(kind=DP) :: diff
        real(kind=DP) :: proc_max_diff_arr(3) ! 1: max_diff, 2: max_val1, 3: max_val2
        integer :: ihandle

        ! Element-by-element comparison
        ! Get number of rows (i.e. number of elements in a column)
        num_rows = sparse_num_rows(mat1)
        call utils_assert(&
             num_rows.eq.sparse_num_rows(mat2),&
             "Error in test_integrals_div_locpot_grad_functions: &
             &div_locpot_grad and kinet SPAM3 matrices do not have same number &
             &of column elements (rows).")
        allocate(col_data1(num_rows),stat=ierr)
        call utils_alloc_check('element_by_element_comparison',&
             'col_data1',ierr)
        allocate(col_data2(num_rows),stat=ierr)
        call utils_alloc_check('element_by_element_comparison',&
             'col_data2',ierr)

        ! Get first column on proc
        first_col = sparse_first_elem_on_proc(pub_my_proc_id,mat1,'C')
        num_cols  = sparse_num_elems_on_proc(pub_my_proc_id,mat1,'C')
        call utils_assert(&
             (num_cols.eq.int(&
             sparse_num_elems_on_proc(pub_my_proc_id,mat2,'C')).and.&
             first_col.eq.int(&
             sparse_first_elem_on_proc(pub_my_proc_id,mat2,'C'))),&
             "Error in test_integrals_div_locpot_grad_functions: &
             &div_locpot_grad and kinet SPAM3 matrices are not distributed &
             &across procs in the same way.")
        ! Check data from sparse module agrees with data from NGWF basis
        call utils_assert( num_cols.eq.fbasis%proc_num, &
             "Error in test_integrals_div_locpot_grad_functions: &
             &data from function basis and sparse algegbra module do not agree &
             &( sparse_num_elems_on_proc / fbasis%num_proc ).")
        call utils_assert( first_col.eq.fbasis%first_on_proc(pub_my_proc_id), &
             "Error in test_integrals_div_locpot_grad_functions: &
             &data from function basis and sparse algegbra module do not agree &
             &( sparse_first_elem_on_proc / fbasis%first_on_proc ).")
        call utils_assert( num_rows.eq.fbasis%num, &
             "Error in test_integrals_div_locpot_grad_functions: &
             &data from function basis and sparse algegbra module do not agree &
             &( sparse_num_rows / fbasis%num ).")

        col_data1 = 0.0_DP; col_data2 = 0.0_DP
        proc_max_diff_arr = 0.0_DP
        do local_col = first_col, first_col + num_cols - 1
           call sparse_get_col(col_data1,mat1,local_col)
           call sparse_get_col(col_data2,mat2,local_col)
           do local_row = 1, num_rows
              val1 = col_data1(local_row); val2 = col_data2(local_row)
              diff = abs(val1-val2)
              if (diff.gt.max_diff_arr(1)) then
                proc_max_diff_arr(1) = diff
                proc_max_diff_arr(2) = val1
                proc_max_diff_arr(3) = val2
              end if
           end do
        end do

        !col_data1 = 0.0_DP; col_data2 = 0.0_DP
        !proc_max_diff_arr = 0.0_DP
        !do local_row = 1, sparse_proc_num_element(mat1)
        !   val1 = mat1%dmtx(local_row); val2 = mat2%dmtx(local_row)
        !   diff = abs(val1-val2)
        !   if (diff.gt.max_diff_arr(1)) then
        !     proc_max_diff_arr(1) = diff
        !     proc_max_diff_arr(2) = val1
        !     proc_max_diff_arr(3) = val2
        !   end if
        !end do

        if (.not.pub_on_root) then
           !write(stdout,'(a,i3,a,3es24.16)') &
           !  "Data sent from proc ", pub_my_proc_id, ":", &
           !  proc_max_diff_arr(1:3)
           call comms_send(pub_root_proc_id,proc_max_diff_arr,length=3,return_handle=ihandle)
           ! Use comms_wait to block process until root proc has received
           call comms_wait(ihandle,free_stack=.true.)
        end if
        if (pub_on_root) then
           max_diff_arr(:) = proc_max_diff_arr(:)
           !write(stdout,'(a,3x,a,3es24.16)') &
           !  "Data from root proc     ",":", &
           !  proc_max_diff_arr(1:3)
           if (pub_total_num_procs.gt.1) then
              do iproc = 0, pub_total_num_procs-1
                 if (iproc.ne.pub_root_proc_id) then
                    call comms_recv(iproc,proc_max_diff_arr,length=3)
                    !write(stdout,'(a,i3,a,3es24.16)') &
                    !  "Data received from proc ", iproc, ":", &
                    !  proc_max_diff_arr(1:3)
                    if (proc_max_diff_arr(1).gt.max_diff_arr(1)) then
                      max_diff_arr(:) = proc_max_diff_arr(:)
                    end if
                 end if
              end do
           end if
        end if

        deallocate(col_data1,stat=ierr)
        call utils_dealloc_check('element_by_element_comparison',&
             'col_data1',ierr)
        deallocate(col_data2,stat=ierr)
        call utils_dealloc_check('element_by_element_comparison',&
             'col_data2',ierr)
    end subroutine element_by_element_comparison

    subroutine trace_comparison(mat1,mat2,mat3,trace1,trace2)
      !==========================================================================!
      ! Traces mat1 with mat3, trace mat2 with mat3 and returns the traces in    !
      ! trace1, and trace2, respectively                                         !
      !==========================================================================!
      use sparse, only: SPAM3, sparse_trace

      implicit none

      ! Arguments
      type(SPAM3),intent(in) :: mat1, mat2, mat3
      real(kind=DP),intent(out) :: trace1, trace2

      trace1 = &
          sparse_trace(mat3,mat1)
      trace2 = &
          sparse_trace(mat3,mat2)

    end subroutine trace_comparison

  end subroutine test_integrals_div_locpot_grad_functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine integrals_locpot_dbl_grid(locpot,  &                ! input-output
       bras_on_grid, bra_basis, kets_on_grid, ket_basis, potential_dbl, &
       dbl_grid, cell, fftbox)

    !==========================================================================!
    ! This subroutine calculates local potential integrals between two         !
    ! sets of functions.                                                       !
    ! Result is locpot_\alpha\beta = < bra_\alpha | Vloc | ket_\beta >         !
    !==========================================================================!
    !  Arguments:                                                              !
    !    locpot     (inout)  : sparse local potential matrix elements          !
    !    bras_on_grid  (in)  : bra functions on grid in PPD format             !
    !    bra_basis (input)   : The basis type for the bras                     !
    !    kets_on_grid  (in)  : ket functions on grid in PPD format             !
    !    ket_basis (input)   : The basis type for the kets                     !
    !    potential_fine (in) : local potential on fine sim cell grid           !
    !==========================================================================!
    ! Originaly written by Chris-Kriton Skylaris in January 2001,              !
    ! to use a "pair-box".                                                     !
    ! Improved by Oswaldo Dieguez in 2001 so that it "cshifts" the "pair-box". !
    ! Rewritten by Arash Mostofi in 2003 so that it uses a "triple-box".       !
    ! Rewritten by Chris-Kriton Skylaris on 23/11/2003 so that it runs on      !
    ! parallel computers.                                                      !
    ! Modifications and speed-ups by Peter Haynes 16/03/2005                   !
    ! Modified to use planned sums system by Nicholas Hine, May 2008           !
    ! Symmetrisation options added by Nicholas Hine, July 2008.                !
    ! Modified to remove integer buffers by Nicholas Hine, June 2009           !
    ! Adapted to use SPAM3 matrices by Nicholas Hine, June 2009                !
    ! Modified to take different function sets for bras and kets by            !
    ! Nicholas Hine on 11/11/2009.                                             !
    ! Modified to skip symmetrisation of the locpot matrix for the case of     !
    ! cross-Hamiltonians (bras/=kets) by Nicholas Hine, Oct 2010.              !
    ! Split off from integrals_locpot for clarity of above by Nicholas Hine in !
    ! March 2011.                                                              !
    ! Modified by Andrea Greco in May 2015 to allow use of complex NGWFs.      !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_barrier, comms_reduce, pub_on_root
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    ! agrecocmplx
    use function_ops, only: function_ops_batch_col_start, &
        function_ops_brappd_ketfftbox
#ifdef GPU_PGI
    use potential, only: potential_gpu_app2_ngwf_batch
    use cudafor !! External dependency
    use fourier_gpu_wrapper_mod
#else
    use potential, only: potential_apply_to_ngwf_batch
#endif
    use rundat, only: pub_fftbox_batch_size, pub_locpot_scheme
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_generate_index, sparse_index_length, &
         sparse_scale, sparse_create, sparse_transpose, sparse_axpy, &
         sparse_destroy, sparse_expand, pattern_lower, pattern_alternate
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: locpot
    type(FUNC_BASIS), intent(in) :: bra_basis
    !real(kind=DP), intent(in) :: bras_on_grid(bra_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: bras_on_grid
    type(FUNC_BASIS), intent(in) :: ket_basis
    !real(kind=DP), intent(in) :: kets_on_grid(ket_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: kets_on_grid
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(in) :: potential_dbl(dbl_grid%ld1, &
         dbl_grid%ld2, dbl_grid%max_group_slabs12)

    ! Local Variables
    !real(kind=DP), allocatable, dimension(:,:,:,:) :: potential_fftbox_batch
    ! agrecocmplx: use new FFTBOX_DATA type
    type(FFTBOX_DATA), allocatable, dimension(:) :: potential_fftbox_batch
    integer :: batch_size
    integer :: batch_count
    integer :: n_batches
    integer :: local_start, local_end
    integer :: max_current_size ! maximum batch size ovar all procs
    integer :: idx_len  ! pdh: sparse matrix index length
    integer :: ierr ! pdh: error flag
    integer, allocatable, dimension(:) :: locpot_idx ! pdh: sparse index
    integer, allocatable :: ket_box_start(:,:)
    integer, allocatable :: ket_start_in_box(:,:)
    type(SPAM3) :: locpot_transpose  ! ndmh: for symmetrising locpot
    character(len=80) :: locpot_scheme
    ! agrecocmplx
    logical :: loc_cmplx, ongpu

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering integrals_locpot_dbl_grid'

    ! Start timer
    call timer_clock('integrals_locpot_dbl_grid',1)

    ! agrecocmplx
    loc_cmplx = kets_on_grid%iscmplx.or.bras_on_grid%iscmplx
    ! Set batch size
    batch_size = pub_fftbox_batch_size

    ! pdh: obtain index for locpot
    idx_len = sparse_index_length(locpot)
    allocate(locpot_idx(idx_len),stat=ierr)
    call utils_alloc_check('integrals_locpot_dbl_grid','locpot_idx',ierr)
    call sparse_generate_index(locpot_idx,locpot)

    ! rc2013: check that the bases refer to the same region
    if (bra_basis%name/=ket_basis%name .or. bra_basis%ireg/=ket_basis%ireg) then
       locpot_scheme = 'FULL'
    else
       locpot_scheme = pub_locpot_scheme
    end if

    ! Allocate workspace
    allocate(potential_fftbox_batch(batch_size), stat=ierr)
    call utils_alloc_check('integrals_locpot_dbl_grid','potential_fftbox_batch',ierr)
    ! agrecocmplx: potential_fftbox_batch complex if NGWFs are complex
    do batch_count =1, batch_size
        call data_fftbox_alloc(potential_fftbox_batch(batch_count), &
             fftbox%total_ld1, fftbox%total_ld2, fftbox%total_pt3, iscmplx=loc_cmplx)
    end do

    allocate(ket_box_start(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_locpot_dbl_grid','ket_box_start',ierr)
    allocate(ket_start_in_box(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_locpot_dbl_grid','ket_start_in_box',ierr)

    ! cks: number of row-steps per row-block
    n_batches = ket_basis%max_on_proc / batch_size
    if (mod(ket_basis%max_on_proc,batch_size) > 0) n_batches = n_batches + 1

    local_start = 1
    do batch_count=1,n_batches

    if (pub_debug_on_root) write(stdout,'(a,2(i4,a))') &
         'DEBUG: Batch loop',batch_count, ' of',n_batches, &
         ' in integrals_locpot_dbl_grid'

       local_end = min(local_start+batch_size-1,ket_basis%proc_num)

       ! cks: maximum size of current batch over all procs
       max_current_size = local_end - local_start + 1
       call comms_reduce('MAX',max_current_size)

       call function_ops_batch_col_start(ket_box_start,ket_start_in_box, &
            batch_size,local_start,local_end,fftbox,cell,ket_basis)

       ! cks: apply locpot operator on to kets
#ifdef GPU_PGI
       call potential_gpu_app2_ngwf_batch(potential_fftbox_batch, &
            kets_on_grid, ket_basis, potential_dbl, dbl_grid, fftbox, cell, &
            ket_start_in_box, ket_box_start, batch_size, local_start, &
            local_end, max_current_size)
            ongpu = .true.
#else
       call potential_apply_to_ngwf_batch(potential_fftbox_batch, &
            kets_on_grid, ket_basis, potential_dbl, dbl_grid, fftbox, cell, &
            ket_start_in_box, ket_box_start, batch_size, local_start, &
            local_end, max_current_size)
            ongpu = .false.
#endif

       ! ndmh: calculate locpot integrals in SPAM3 format
       ! rc2013: use local locpot_scheme, rather than public version
       call timer_clock('integrals_locpot_mat_els_batch',1)
       call function_ops_brappd_ketfftbox(locpot,                    & ! inout
            bras_on_grid, bra_basis, cell, fftbox,                   & ! input
            potential_fftbox_batch, ket_box_start, batch_size,       & ! input
            local_start, local_end, idx_len, locpot_idx,             & ! input
            locpot_scheme,ongpu)                                       ! input
       call timer_clock('integrals_locpot_mat_els_batch',2)
       local_start = local_start + batch_size

    end do

    ! Deallocate workspace
    deallocate(ket_start_in_box,stat=ierr)
    call utils_dealloc_check('integrals_locpot_dbl_grid','ket_start_in_box',ierr)
    deallocate(ket_box_start,stat=ierr)
    call utils_dealloc_check('integrals_locpot_dbl_grid','ket_box_start',ierr)
    ! agrecocmplx: deallocate using appropriate routine
    do batch_count =1, batch_size
        call data_fftbox_dealloc(potential_fftbox_batch(batch_count))
    end do
    deallocate(potential_fftbox_batch,stat=ierr)
    call utils_dealloc_check('integrals_locpot_dbl_grid','potential_fftbox_batch',ierr)
    deallocate(locpot_idx,stat=ierr)
    call utils_dealloc_check('integrals_locpot_dbl_grid','locpot_idx',ierr)

    ! ndmh: expand to full matrix
    select case (locpot_scheme)
    case ('LOWER','lower')
          call sparse_expand(locpot,PATTERN_LOWER)
    case ('ALTERNATE','alternate')
          call sparse_expand(locpot,PATTERN_ALTERNATE)
    case ('FULL','full')
          ! ndmh: only symmetrise if bra_basis==ket_basis
          ! rc2013: and if the bases are in the same region
          if (bra_basis%name==ket_basis%name .and. bra_basis%ireg==ket_basis%ireg) then
             call sparse_scale(locpot,0.5_DP)
             call sparse_create(locpot_transpose, locpot)
             call sparse_transpose(locpot_transpose, locpot)
             call sparse_axpy(locpot, locpot_transpose, 1.0_DP)
             call sparse_destroy(locpot_transpose)
          end if
    case ('ASYM','asym')
       if (pub_on_root) write(stdout,'(a)') &
            'WARNING: Skipping symmetrisation of locpot in integrals_locpot - '
       if (pub_on_root) write(stdout,'(a)') &
            'Resulting matrix may not be Hermitian.'
    case default
       call utils_abort('Error in integrals_locpot: calculation pattern"'//&
            trim(pub_locpot_scheme)//'" not recognised')
    end select

    ! pdh: re-sync procs
    call comms_barrier

    ! Stop timer
    call timer_clock('integrals_locpot_dbl_grid',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &integrals_locpot_dbl_grid'

  end subroutine integrals_locpot_dbl_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine integrals_kinetic(kinet, bras_on_grid, bra_basis, kets_on_grid, &
       ket_basis, cell, fftbox)

    !========================================================================!
    ! This subroutine calculates kinetic integrals between two function sets.!
    ! Result is kinet_\alpha\beta = < bra_\alpha | \nabla^2 | ket_\beta >    !
    !========================================================================!
    !  Arguments:                                                            !
    !    kinet      (inout) : sparse kinetic energy matrix elements          !
    !    bras_on_grid  (in) : bra functions on grid in PPD format            !
    !    bra_basis (input)  : The basis type for the bras                    !
    !    kets_on_grid  (in) : ket functions on grid in PPD format            !
    !    ket_basis (input)  : The basis type for the kets                    !
    !========================================================================!
    ! Originaly written by Chris-Kriton Skylaris in January 2001,            !
    ! to use a "pair-box".                                                   !
    ! Rewritten by Arash Mostofi in 2003 so that it uses a "triple-box".     !
    ! Rewritten by Chris-Kriton Skylaris on 20/11/2003 so that it runs on    !
    ! parallel computers.                                                    !
    ! Modifications and speed-ups by Peter Haynes 16/03/2005                 !
    ! Further modification to use parallel SPAM 2, July 2006                 !
    ! Modified to remove integer buffers by Nicholas Hine, June 2009         !
    ! Adapted to use SPAM3 matrices by Nicholas Hine, June 2009              !
    ! Modified to take different function sets for bras and kets by          !
    ! Nicholas Hine on 11/11/2009.                                           !
    ! Modified to calculate full pattern of the kinetic matrix for the case  !
    ! of valence-conduction Hamiltonians by Laura Ratcliff, Oct 2010.        !
    ! Modified by Andrea Greco on 31/05/2015 to allow use of complex NGWFs.  !
    !========================================================================!

    use comms, only: comms_barrier
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_batch_col_start, &
         function_ops_brappd_ketfftbox
#ifdef GPU_PGI_KINETIC
    use kinetic, only: kinetic_gpu_app2_func_batch
#else
    use kinetic, only: kinetic_apply_to_func_batch
#endif
    use rundat, only: pub_fftbox_batch_size
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_generate_index, sparse_index_length, &
         sparse_expand, PATTERN_ALTERNATE
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: kinet
    type(FUNC_BASIS), intent(in) :: bra_basis
    !real(kind=DP), intent(in) :: bras_on_grid(bra_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: bras_on_grid
    type(FUNC_BASIS), intent(in) :: ket_basis
    !real(kind=DP), intent(in) :: kets_on_grid(ket_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: kets_on_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local Variables
    type(FFTBOX_DATA), allocatable :: kinetic_fftbox_batch(:)
    integer :: batch_size
    integer :: batch_count
    integer :: n_batches
    integer :: local_start, local_end
    integer :: ierr                            ! pdh: error flag
    integer :: idx_len                         ! pdh: length of sparse index
    integer, allocatable :: kinet_idx(:)       ! pdh: sparse index
    integer, allocatable :: ket_box_start(:,:)
    integer, allocatable :: ket_start_in_box(:,:)
    character(20) :: kinetic_scheme
    ! agrecocmplx: local variable for complex NGWFs
    logical :: loc_cmplx
    logical :: ongpu

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering integrals_kinetic'

    ! Start timer
    call timer_clock('integrals_kinetic',1)

    ! Obtain index for kinet
    idx_len = sparse_index_length(kinet)
    allocate(kinet_idx(idx_len),stat=ierr)
    call utils_alloc_check('integrals_kinetic','kinet_idx',ierr)
    call sparse_generate_index(kinet_idx,kinet)

    ! Define shorthand variables
    batch_size = pub_fftbox_batch_size

    ! agrecocmplx
    loc_cmplx = kets_on_grid%iscmplx

    ! lr408: Force calculation of the full integral if this is for a cross matrix
    ! rc2013: check that the bases refer to the same region
    if (bra_basis%name/=ket_basis%name .or. bra_basis%ireg/=ket_basis%ireg) then
       kinetic_scheme = 'FULL'
    else
       kinetic_scheme = 'ALTERNATE'
    end if

    ! Allocate workspace
    allocate(kinetic_fftbox_batch(batch_size), stat=ierr)
    call utils_alloc_check('integrals_kinetic','kinetic_fftbox_batch',ierr)
    do batch_count = 1, batch_size
        call data_fftbox_alloc(kinetic_fftbox_batch(batch_count), &
             fftbox%total_ld1, fftbox%total_ld2, fftbox%total_pt3, &
             iscmplx=loc_cmplx)
    end do
    allocate(ket_box_start(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_kinetic','ket_box_start',ierr)
    allocate(ket_start_in_box(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_kinetic','ket_start_in_box',ierr)

    ! cks: number of row-steps per row-block
    n_batches = ket_basis%max_on_proc / batch_size
    if (mod(ket_basis%max_on_proc,batch_size) > 0) n_batches = n_batches + 1

    local_start = 1
    do batch_count=1,n_batches

       if (pub_debug_on_root) write(stdout,'(a,2(i4,a))') 'DEBUG: Batch loop ',&
            batch_count, ' of ',n_batches,' in integrals_kinetic'

       local_end = min(local_start+batch_size-1,ket_basis%proc_num)

       call function_ops_batch_col_start(ket_box_start,ket_start_in_box, &
            batch_size,local_start,local_end,fftbox,cell,ket_basis)

       ! cks: apply kinetic operator on to kets
#ifdef GPU_PGI_KINETIC
       call kinetic_gpu_app2_func_batch(kinetic_fftbox_batch, &          ! out
            kets_on_grid, ket_basis, cell, fftbox, ket_start_in_box, &
            batch_size, local_start, local_end) ! in
            ongpu = .true.
#else
       call kinetic_apply_to_func_batch(kinetic_fftbox_batch, &          ! out
            kets_on_grid, ket_basis, cell, fftbox, ket_start_in_box, &
            batch_size, local_start, local_end) ! in
            ongpu = .false.
#endif

       if (pub_debug_on_root) write(stdout,'(a,2(a))') 'DEBUG:  ', &
            'after kinetic_apply_to_func_batch'

       ! ndmh: calculate kinetic energy integrals in SPAM3 format
       call function_ops_brappd_ketfftbox(kinet, &
            bras_on_grid, bra_basis, cell, fftbox, &
            kinetic_fftbox_batch, ket_box_start, batch_size, &
            local_start, local_end, idx_len, kinet_idx, kinetic_scheme,ongpu)

       local_start = local_start + batch_size

    end do

    ! ndmh: expand to full matrix
    if (kinetic_scheme=='ALTERNATE') then
       call sparse_expand(kinet,PATTERN_ALTERNATE)
    end if

    ! pdh: deallocate workspace
    deallocate(ket_start_in_box,stat=ierr)
    call utils_dealloc_check('integrals_kinetic','ket_start_in_box',ierr)
    deallocate(ket_box_start,stat=ierr)
    call utils_dealloc_check('integrals_kinetic','ket_box_start',ierr)
    ! agrecocmplx: deallocate using appropriate routine
    do batch_count = 1, batch_size
        call data_fftbox_dealloc(kinetic_fftbox_batch(batch_count))
    end do
    deallocate(kinetic_fftbox_batch,stat=ierr)
    call utils_dealloc_check('integrals_kinetic','kinetic_fftbox_batch',ierr)
    deallocate(kinet_idx,stat=ierr)
    call utils_dealloc_check('integrals_kinetic','kinet_idx',ierr)

    ! pdh: sync procs
    call comms_barrier

    ! Stop timer
    call timer_clock('integrals_kinetic',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving integrals_kinetic'

  end subroutine integrals_kinetic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine integrals_grad(grad, bras_on_grid, bra_basis, kets_on_grid, &
       ket_basis, cell, fftbox)

    !=======================================================================!
    ! This subroutine calculates grad integrals between two function sets.  !
    ! Result is grad(i)_\alpha\beta = < bra_\alpha | \nabla_i | ket_\beta > !
    !=======================================================================!
    !  Arguments:                                                           !
    !    grad(3)    (inout) : sparse grad matrix elements                   !
    !    bras_on_grid  (in) : bra functions on grid in PPD format           !
    !    bra_basis (input)  : The basis type for the bras                   !
    !    kets_on_grid  (in) : ket functions on grid in PPD format           !
    !    ket_basis (input)  : The basis type for the kets                   !
    !=======================================================================!
    ! Originally written by Chris-Kriton Skylaris in January 2001,          !
    ! to use a "pair-box".                                                  !
    ! Rewritten by Arash Mostofi in 2003 so that it uses a "triple-box".    !
    ! Rewritten by Chris-Kriton Skylaris on 20/11/2003 so that it runs on   !
    ! parallel computers.                                                   !
    ! Modifications and speed-ups by Peter Haynes 16/03/2005                !
    ! Further modification to use parallel SPAM 2, July 2006                !
    ! Modified to remove integer buffers by Nicholas Hine, June 2009.       !
    ! Adapted for SPAM3 by Nicholas Hine, July 2009.                        !
    ! Modified to take different function sets for bras and kets by         !
    ! Nicholas Hine on 11/11/2009.                                          !
    ! Modified by Andrea Greco on 16/06/2015 to allow use of complex NGWfs. !
    !=======================================================================!

    use comms, only: comms_barrier
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, &
         FFTBOX_DATA, data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    ! agrecocmplx
    use function_ops, only: function_ops_batch_col_start, &
         function_ops_brappd_ketfftbox
    ! agrecocmplx
    use kinetic, only: kinetic_grad_to_func_batch
    use rundat, only: pub_fftbox_batch_size
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_generate_index, sparse_index_length, &
         sparse_expand, pattern_alternate
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: grad(3)
    type(FUNC_BASIS), intent(in) :: bra_basis
    !real(kind=DP), intent(in) :: bras_on_grid(bra_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: bras_on_grid
    type(FUNC_BASIS), intent(in) :: ket_basis
    !real(kind=DP), intent(in) :: kets_on_grid(ket_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: kets_on_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local Variables
    !real(kind=DP), allocatable, dimension(:,:,:,:,:) :: grad_fftbox_batch
    ! agrecocmplx: use new FFTBOX_DATA type
    type(FFTBOX_DATA), allocatable, dimension(:,:) :: grad_fftbox_batch
    integer :: batch_size
    integer :: batch_count
    integer :: n_batches
    integer :: local_start, local_end
    integer :: ierr         ! pdh: error flag
    integer :: dim          ! pdh: cartesian direction
    integer :: idx_len      ! pdh: length of sparse index
    integer, allocatable :: ket_box_start(:,:)
    integer, allocatable :: ket_start_in_box(:,:)
    integer, allocatable, dimension(:) :: grad_idx   ! pdh: sparse index
    character(20) :: pattern
    logical :: loc_cmplx ! agrecocmplx
    logical :: ongpu

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering integrals_grad'

    ! Start timer
    call timer_clock('integrals_grad',1)

    ! agrecocmplx
    loc_cmplx = kets_on_grid%iscmplx.or.bras_on_grid%iscmplx

    ! Obtain index for grad
    idx_len = sparse_index_length(grad(1))
    allocate(grad_idx(idx_len),stat=ierr)
    call utils_alloc_check('integrals_grad','grad_idx',ierr)
    call sparse_generate_index(grad_idx,grad(1))

    ! Define shorthand variables
    batch_size = pub_fftbox_batch_size

    ! ddor: FULL integral is calculated for non-square matrices
    ! lr408: changed to different basis names
    ! rc2013: check that the bases refer to the same region
    if (bra_basis%name .eq. ket_basis%name .or. bra_basis%ireg/=ket_basis%ireg) then
       pattern='ALTERNATE'
    else
       pattern='FULL'
    end if

    ! Allocate workspace
    !allocate(grad_fftbox_batch(fftbox%total_ld1, fftbox%total_ld2, &
    !     fftbox%total_pt3, batch_size, 3), stat=ierr)
    !call utils_alloc_check('integrals_grad','grad_fftbox_batch',ierr)
    ! agrecocmplx: allocate using appropriate routine
    allocate(grad_fftbox_batch(batch_size,3), stat=ierr)
    call utils_alloc_check('integrals_grad','grad_fftbox_batch',ierr)
    do dim = 1,3
        do batch_count = 1,batch_size
            call data_fftbox_alloc(grad_fftbox_batch(batch_count,dim), &
                 fftbox%total_ld1, fftbox%total_ld2, fftbox%total_pt3, &
                 iscmplx=loc_cmplx)
        end do
    end do
    allocate(ket_box_start(3,batch_size),stat=ierr)
    call utils_alloc_check('integrals_grad','ket_box_start',ierr)
    allocate(ket_start_in_box(3,batch_size),stat=ierr)
    call utils_alloc_check('integrals_grad','ket_start_in_box',ierr)

    ! cks: number of row-steps per row-block
    n_batches = ket_basis%max_on_proc / batch_size
    if (mod(ket_basis%max_on_proc,batch_size) > 0) n_batches = n_batches + 1

    local_start = 1
    do batch_count=1,n_batches

       if (pub_debug_on_root) write(stdout,'(a,2(i4,a))') 'DEBUG: Batch loop ', &
            batch_count,' of ',n_batches,' in integrals_grad'

       local_end = min(local_start+batch_size-1,ket_basis%proc_num)

       call function_ops_batch_col_start(ket_box_start,ket_start_in_box, &
            batch_size,local_start,local_end,fftbox,cell,ket_basis)

       ! cks: apply grad operator on to kets
       call kinetic_grad_to_func_batch(grad_fftbox_batch, &
            kets_on_grid, ket_basis, cell, fftbox, ket_start_in_box, &
            batch_size, local_start, local_end)  ! in

       if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: after &
            &kinetic_grad_to_func_batch'

       do dim=1,3

          ! ndmh: calculate grad integrals in SPAM3 format
          ! agrecocmplx
          ! kaw: Leaving this here as these integrals are a future GPU target.
          ongpu = .false.
          call function_ops_brappd_ketfftbox(grad(dim),&              ! inout
               bras_on_grid, bra_basis, cell, fftbox, &            ! in
               grad_fftbox_batch(:,dim), &                   ! in
               ket_box_start, batch_size, local_start, local_end,& ! in
               idx_len, grad_idx, pattern,ongpu)                         ! in
       end do

       local_start = local_start + batch_size

    end do

    ! pdh: expand to full matrix (which is antisymmetric)
    ! ddor: non-square matrix needs to be calculated fully, otherwise expand
    if (pattern=='ALTERNATE') then
       do dim=1,3
          call sparse_expand(grad(dim),PATTERN_ALTERNATE,.false.)
       end do
    endif

    ! pdh: deallocate workspace
    deallocate(ket_start_in_box,stat=ierr)
    call utils_dealloc_check('integrals_grad','ket_start_in_box',ierr)
    deallocate(ket_box_start,stat=ierr)
    call utils_dealloc_check('integrals_grad','ket_box_start',ierr)
    !deallocate(grad_fftbox_batch,stat=ierr)
    !call utils_dealloc_check('integrals_grad','grad_fftbox_batch',ierr)
    ! agrecocmplx: deallocate using appropriate routines
    do dim = 1,3
        do batch_count = 1,batch_size
            call data_fftbox_dealloc(grad_fftbox_batch(batch_count,dim))
        end do
    end do
    deallocate(grad_fftbox_batch,stat=ierr)
    call utils_dealloc_check('integrals_grad','grad_fftbox_batch',ierr)
    deallocate(grad_idx,stat=ierr)
    call utils_dealloc_check('integrals_grad','grad_idx',ierr)

    ! pdh: sync procs
    call comms_barrier

    ! Stop timer
    call timer_clock('integrals_grad',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving integrals_grad'

  end subroutine integrals_grad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !=======================================================================!
    ! This subroutine calculates matrix elements of complex exponentials of !
    ! the position operator (within the FFTbox) between two function sets.  !
    ! Gives exr(i)_\alpha\beta =                                            !
    !         < bra_\alpha | exp(qp.r^FFT_i^pow) | ket_\beta >              !
    ! where r^FFT_i = r_i - R_FFT, with R_FFT being the origin of the beta  !
    ! function's FFTbox.                                                    !
    ! This function does not account for PAW augmentation charges and the   !
    ! result must be post-processed to include them.                        !
    !=======================================================================!
    !  Arguments:                                                           !
    !    exr(3)     (inout) : sparse position operator matrix elements      !
    !    bras_on_grid  (in) : bra functions on grid in PPD format           !
    !    bra_basis  (input) : The basis type for the bras                   !
    !    kets_on_grid  (in) : ket functions on grid in PPD format           !
    !    ket_basis  (input) : The basis type for the kets                   !
    !    qp         (input) : vector in reciprocal space                    !
    !    axis   (opt input) : The axis to calculate, if only doing one.     !
    !=======================================================================!
    ! Written by Simon M.-M. Dubois in October 2013 on the basis of         !
    ! the integrals_pos subroutine                                          !
    ! Modified by Andrea Greco on 18/06/2015 to allow use of complex NGWFs. !
    !=======================================================================!

  subroutine integrals_exr(exr, bras_on_grid, bra_basis, kets_on_grid, &
       ket_basis,dbl_grid, cell, fftbox, qp,axis)

    use basis, only: basis_copy_function_to_box
    use cell_grid, only: GRID_INFO
    use comms, only : comms_reduce
    use constants, only : DP, cmplx_0, cmplx_i
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, &
         data_fftbox_alloc, data_fftbox_dealloc, data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_interpolate, fourier_filter
    use function_basis, only: FUNC_BASIS
    ! agrecocmplx
    use function_ops, only: function_ops_batch_col_start, &
         function_ops_brappd_ketfftbox
    use geometry, only: POINT, operator(*), operator(+), operator(.DOT.) !ddor
    use rundat, only: pub_fftbox_batch_size
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index, & !ddor
         sparse_create, sparse_destroy, sparse_scale, sparse_axpy !smmd
    use utils, only: utils_alloc_check, utils_dealloc_check
    ! Arguments
    type(SPAM3), intent(inout)  :: exr(3)
    type(FUNC_BASIS), intent(in) :: bra_basis
    !real(kind=DP), intent(in) :: bras_on_grid(bra_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: bras_on_grid
    type(FUNC_BASIS), intent(in) :: ket_basis
    !real(kind=DP), intent(in) :: kets_on_grid(ket_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: kets_on_grid
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(POINT), intent(in) :: qp(3)
    ! ddor: Power of position operator for which matrix elements are calculated
    ! ddor: Specifies if we only need one direction
    integer, optional, intent(in)        :: axis

    ! Local Variables
    type(SPAM3) :: exr_imag
    type(POINT) :: r1,r2,r3
    type(POINT) :: a1,a2,a3
    integer :: i1,i2,i3
    integer :: n3,ld1,ld2,ld1_dbl,ld2_dbl
    integer :: max_current_size
    integer :: local_start,local_end
    integer :: batch,n_batches,batch_size,batch_count
    integer, allocatable :: ket_box_start(:,:)
    integer, allocatable :: ket_start_in_box(:,:)
    integer :: ngwf,f1
    integer :: idx_len
    integer :: xyz,ierr
    ! ddor: Integers used to select axis if only one is needed
    integer :: xyz_index,xyz_count

    ! workspace
    !real(kind=DP), allocatable, dimension(:,:,:,:) :: fftbox_batch_real
    !real(kind=DP), allocatable, dimension(:,:,:,:) :: fftbox_batch_imag
    ! agrecocmplx: use new FFTBOX_DATA type
    type(FFTBOX_DATA), allocatable, dimension(:) :: fftbox_batch_real
    type(FFTBOX_DATA), allocatable, dimension(:) :: fftbox_batch_imag
    integer,       allocatable, dimension(:)       :: idx

    ! position operator, FFTbox and tightbox
    complex(kind=DP), dimension(:,:,:,:), allocatable :: exr_op
    ! agrecocmplx: uncommented unused variables, changed fftbox1_real,
    ! fftbox1_imag, fftbox1_dbl_real, fftbox1_dbl_imag  to new type
    !real(kind=DP), dimension(:,:,:), allocatable :: fftbox1_real,fftbox2_real
    !real(kind=DP), dimension(:,:,:), allocatable :: fftbox1_imag,fftbox2_imag
    type(FFTBOX_DATA) :: fftbox1_real, fftbox1_imag
    type(FFTBOX_DATA) :: fftbox1_dbl_real, fftbox1_dbl_imag
    !real(kind=DP), dimension(:,:,:), allocatable :: fftbox1_dbl_real,fftbox2_dbl_real
    !real(kind=DP), dimension(:,:,:), allocatable :: fftbox1_dbl_imag,fftbox2_dbl_imag
    complex(kind=DP), dimension(:,:,:), allocatable :: coarse_work,fine_work
    logical :: loc_cmplx ! agrecocmplx
    logical :: ongpu

    ! local copies of fftbox info
    n3 = fftbox%total_pt3
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2
    ld1_dbl = fftbox%total_ld1_dbl
    ld2_dbl = fftbox%total_ld2_dbl

    ! agrecocmplx
    loc_cmplx = bras_on_grid%iscmplx.or.kets_on_grid%iscmplx

    ! collect batch size and determine number of batches
    batch_size = pub_fftbox_batch_size
    n_batches = ket_basis%max_on_proc / batch_size
    if (mod(ket_basis%max_on_proc,batch_size) > 0) n_batches = n_batches + 1

    ! obtain index for exr (all 3 are the same)
    idx_len = sparse_index_length(exr(1))
    allocate(idx(idx_len),stat=ierr)
    call utils_alloc_check('integrals_exr','idx',ierr)
    call sparse_generate_index(idx,exr(1))

    ! allocate workspace
    !allocate(fftbox_batch_real(ld1,ld2,n3,batch_size),stat=ierr)
    !call utils_alloc_check('integrals_exr','fftbox_batch_real',ierr)
    !allocate(fftbox_batch_imag(ld1,ld2,n3,batch_size),stat=ierr)
    !call utils_alloc_check('integrals_exr','fftbox_batch_imag',ierr)
    ! agrecocmplx: allocate using appropriate routines
    allocate(fftbox_batch_real(batch_size),stat=ierr)
    call utils_alloc_check('integrals_exr','fftbox_batch_real',ierr)
    allocate(fftbox_batch_imag(batch_size),stat=ierr)
    call utils_alloc_check('integrals_exr','fftbox_batch_imag',ierr)
    do batch_count = 1,batch_size
        call data_fftbox_alloc(fftbox_batch_real(batch_count), ld1, &
             ld2, n3, iscmplx=loc_cmplx)
        call data_fftbox_alloc(fftbox_batch_imag(batch_count), ld1, &
             ld2, n3, iscmplx=loc_cmplx)
    end do
    ! allocate position operator and FFTboxes
    allocate(ket_box_start(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_exr','ket_box_start',ierr)
    allocate(ket_start_in_box(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_exr','ket_start_in_box',ierr)
    !allocate(fftbox1_real(ld1,ld2,n3), stat=ierr)
    !call utils_alloc_check('integrals_exr','fftbox1_real',ierr)
    !allocate(fftbox1_imag(ld1,ld2,n3), stat=ierr)
    !call utils_alloc_check('integrals_exr','fftbox1_imag',ierr)
    ! agrecocmplx: allocate using appropriate routines
    call data_fftbox_alloc(fftbox1_real, ld1, ld2, n3, &
         iscmplx=loc_cmplx)
    call data_fftbox_alloc(fftbox1_imag, ld1, ld2, n3, &
         iscmplx=loc_cmplx)
    ! agrecocmplx: not used
    !allocate(fftbox2_real(ld1,ld2,n3), stat=ierr)
    !call utils_alloc_check('integrals_exr','fftbox2_real', ierr)
    !allocate(fftbox2_imag(ld1,ld2,n3), stat=ierr)
    !call utils_alloc_check('integrals_exr','fftbox2_imag', ierr)
    !allocate(fftbox1_dbl_real(ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    !call utils_alloc_check('integrals_exr','fftbox1_dbl_real',ierr)
    !allocate(fftbox1_dbl_imag(ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    !call utils_alloc_check('integrals_exr','fftbox1_dbl_imag',ierr)
    ! agrecocmplx: allocate using appropriate routines
    call data_fftbox_alloc(fftbox1_dbl_real, ld1_dbl, ld2_dbl, 2*n3, &
         iscmplx=loc_cmplx)
    call data_fftbox_alloc(fftbox1_dbl_imag, ld1_dbl, ld2_dbl, 2*n3, &
         iscmplx=loc_cmplx)
    ! agrecocmplx: not used
    !allocate(fftbox2_dbl_real(ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    !call utils_alloc_check('integrals_exr','fftbox2_dbl_real',ierr)
    !allocate(fftbox2_dbl_imag(ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    !call utils_alloc_check('integrals_exr','fftbox2_dbl_imag',ierr)
    allocate(exr_op(3,ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    call utils_alloc_check('integrals_exr','exr_op',ierr)
    allocate(coarse_work(ld1,ld2,n3), stat=ierr)
    call utils_alloc_check('integrals_exr','coarse_work',ierr)
    allocate(fine_work(ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    call utils_alloc_check('integrals_exr','fine_work',ierr)

    ! calculate vectors between fine grid points
    a1 = (1.0_DP / dbl_grid%n1) * cell%a1
    a2 = (1.0_DP / dbl_grid%n2) * cell%a2
    a3 = (1.0_DP / dbl_grid%n3) * cell%a3

    exr_op = cmplx_0
    do i3=1,2*n3
       r3 = real(i3-1,kind=DP) * a3
       do i2=1,ld2_dbl
          r2 = r3 + real(i2-1,kind=DP) * a2
          do i1=1,ld1_dbl
             r1 = r2 + real(i1-1,kind=DP) * a1
             exr_op(1,i1,i2,i3) = exp(cmplx(0.0_DP,qp(1).DOT.r1,kind=DP))
             exr_op(2,i1,i2,i3) = exp(cmplx(0.0_DP,qp(2).DOT.r1,kind=DP))
             exr_op(3,i1,i2,i3) = exp(cmplx(0.0_DP,qp(3).DOT.r1,kind=DP))
          enddo
       enddo
    enddo

    ! ddor: Do one direction only if axis is present
    if (present(axis)) then
       xyz_count = 1
    else
       xyz_count = 3
    endif

    call sparse_create(exr_imag,exr(1),iscmplx=.false.)

    ! loop over x,y,z (or just axis)
    do xyz_index=1,xyz_count

       ! ddor: xyz is the label on the axis direction
       if (present(axis)) then
          xyz = axis
       else
          xyz = xyz_index
       endif

       call sparse_scale(exr(xyz),0.0_DP)
       call sparse_scale(exr_imag,0.0_DP)

       ! loop over batches of NGWFs on this proc to construct
       ! the matrix elements R_ab = < phi_a | exp(i qp.r^pow) | phi_b >
       local_start = 1
       do batch=1,n_batches

          ! find last NGWF in batch
          local_end = min(local_start+batch_size-1,ket_basis%proc_num)
          max_current_size = local_end - local_start
          call comms_reduce('MAX',max_current_size)

          call function_ops_batch_col_start(ket_box_start,ket_start_in_box, &
               batch_size,local_start,local_end,fftbox,cell,ket_basis)

          ! place NGWFs in FFTboxes
          batch_count = 1

          !!do ngwf=local_start,local_start+max_current_size,2
          !!   f1 = ngwf
          !!   f2 = ngwf +1

          do ngwf=local_start,local_start+max_current_size,1
             f1 = ngwf

             ! copy function 1 into first fftbox if required
             if (f1 <= local_end) then
                ! agrecocmplx: use modified routine in basis_new
                call basis_copy_function_to_box(fftbox1_real, &
                     ket_start_in_box(1,batch_count), &
                     ket_start_in_box(2,batch_count), &
                     ket_start_in_box(3,batch_count), &
                     ket_basis%tight_boxes(f1),kets_on_grid, &
                     ket_basis%spheres(f1),cell, fftbox)
                ! agrecocmplx
                call data_set_to_zero(fftbox1_imag)
             else
                ! agrecocmplx
                call data_set_to_zero(fftbox1_real)
                call data_set_to_zero(fftbox1_imag)
             endif

             !!! copy function 2 into second fftbox if required
             !!if (f2 <= local_end) then
             !!   call basis_copy_function_to_box(fftbox2_real,ld1,ld2,n3, &
             !!      ket_start_in_box(1,batch_count), &
             !!      ket_start_in_box(2,batch_count), &
             !!      ket_start_in_box(3,batch_count), &
             !!        ket_basis%tight_boxes(f2),kets_on_grid, &
             !!        ket_basis%spheres(f2))
             !!   fftbox2_imag = 0.0_DP
             !!else
             !!   fftbox2_real = 0.0_DP
             !!   fftbox2_imag = 0.0_DP
             !!endif

             if (f1 <= local_end) then
                ! agrecocmplx: if complex functions, interpolate one at
                ! a time
                if (.not.loc_cmplx) then
                    ! interpolate fftboxes to fine grid
                    !!call fourier_interpolate(coarse_work, fine_work, &
                    !!     fftbox1_real,fftbox2_real,fftbox1_dbl_real, &
                    !!     fftbox2_dbl_real)
                    call fourier_interpolate(coarse_work, fine_work, &
                         fftbox1_real%d,fftbox1_imag%d,fftbox1_dbl_real%d, &
                         fftbox1_dbl_imag%d)

                    ! apply operator to fftbox1
                    fftbox1_dbl_real%d(:,:,:) = real(exr_op(xyz,:,:,:) * fine_work(:,:,:))
                    fftbox1_dbl_imag%d(:,:,:) = aimag(exr_op(xyz,:,:,:) * fine_work(:,:,:))
                    !!fftbox1_dbl_real = real(exr_op(xyz,:,:,:) * fftbox1_dbl_real)
                    !!fftbox1_dbl_imag = aimag(exr_op(xyz,:,:,:) * fftbox1_dbl_real)

                    ! filter fftbox1 to standard grid
                    call fourier_filter(coarse_work, fine_work, &
                         fftbox1_dbl_real%d,fftbox1_dbl_imag%d,fftbox1_real%d,fftbox1_imag%d)

                    !!if (f2 <= local_end) then
                    !!   ! apply operator to fftbox2
                    !!   fftbox2_dbl_real = real(exr_op(xyz,:,:,:) * fftbox2_dbl_real)
                    !!   fftbox2_dbl_imag = aimag(exr_op(xyz,:,:,:) * fftbox2_dbl_real)
                    !!
                    !!   ! filter fftbox1 to standard grid
                    !!   call fourier_filter(coarse_work, fine_work, &
                    !!        fftbox2_dbl_real,fftbox2_dbl_imag,fftbox2_real,fftbox2_imag)
                    !!
                    !!endif
                else
                    ! agrecocmplx: only 1 fftbox is needed in this case
                    call fourier_interpolate(coarse_work, fine_work, fftbox1_real%z)

                    fftbox1_dbl_real%z(:,:,:) = exr_op(xyz,:,:,:) * fine_work(:,:,:)

                    call fourier_filter(coarse_work, fine_work, fftbox1_dbl_real%z, &
                         fftbox1_real%z)

                end if

             endif

             ! agrecocmplx
             if (.not.loc_cmplx) then

                 ! copy functions into batch
                 fftbox_batch_real(batch_count)%d(:,:,:) = fftbox1_real%d(:,:,:)
                 fftbox_batch_imag(batch_count)%d(:,:,:) = fftbox1_imag%d(:,:,:)

             else

                 fftbox_batch_real(batch_count)%z(:,:,:) = fftbox1_real%z(:,:,:)

             end if


             batch_count = batch_count +1
             !!if (f2 <= local_end) then
             !!   fftbox_batch_real(:,:,:,batch_count) = fftbox2_real
             !!   fftbox_batch_imag(:,:,:,batch_count) = fftbox2_imag
             !!   batch_count = batch_count +1
             !!endif
          enddo

          ! ndmh: calculate operator matrix elements in SPAM3 format
          ongpu=.false.
          call function_ops_brappd_ketfftbox(exr(xyz),               & ! inout
               bras_on_grid, bra_basis, cell, fftbox,             & ! input
               fftbox_batch_real, ket_box_start, batch_size,      & ! input
               local_start, local_end, idx_len, idx, 'FULL',ongpu)

          ! agrecocmplx: in the real case, need to add the imaginary
          ! part explicitly
          if (.not.loc_cmplx) then

              ongpu=.false.
              call function_ops_brappd_ketfftbox(exr_imag,               & ! inout
                   bras_on_grid, bra_basis, cell, fftbox,             & ! input
                   fftbox_batch_imag, ket_box_start, batch_size,      & ! input
                   local_start, local_end, idx_len, idx, 'FULL',ongpu)

          end if

          ! increment batch start
          local_start = local_start + batch_size

       enddo  ! end batch loop

       ! agrecocmplx: in the real case, need to add the imaginary
       ! part explicitly
       if (.not.loc_cmplx) then

           call sparse_axpy(exr(xyz),exr_imag,cmplx_i)

       end if

    enddo  ! end loop over xyz

    call sparse_destroy(exr_imag)

    ! deallocate FFTbox and tightbox
    deallocate(fine_work,stat=ierr)
    call utils_dealloc_check('integrals_exr','fine_work',ierr)
    deallocate(coarse_work,stat=ierr)
    call utils_dealloc_check('integrals_exr','coarse_work',ierr)
    deallocate(exr_op,stat=ierr)
    call utils_dealloc_check('integrals_exr','exr_op',ierr)
    ! agrecocmplx: not used
    !deallocate(fftbox2_dbl_imag,stat=ierr)
    !call utils_dealloc_check('integrals_exr','fftbox2_dbl_imag',ierr)
    !deallocate(fftbox2_dbl_real,stat=ierr)
    !call utils_dealloc_check('integrals_exr','fftbox2_dbl_real',ierr)
    !deallocate(fftbox2_imag,stat=ierr)
    !call utils_dealloc_check('integrals_exr','fftbox2_imag',ierr)
    !deallocate(fftbox2_real,stat=ierr)
    !call utils_dealloc_check('integrals_exr','fftbox2_real',ierr)
    ! agrecocmplx: deallocate using appropriate routines
    call data_fftbox_dealloc(fftbox1_dbl_imag)
    call data_fftbox_dealloc(fftbox1_dbl_real)
    call data_fftbox_dealloc(fftbox1_imag)
    call data_fftbox_dealloc(fftbox1_real)
    deallocate(ket_start_in_box,stat=ierr)
    call utils_dealloc_check('integrals_exr','ket_start_in_box',ierr)
    deallocate(ket_box_start,stat=ierr)
    call utils_dealloc_check('integrals_exr','ket_box_start',ierr)

    ! deallocate workspace
    ! agrecocmplx: deallocate using appropriate routines
    do batch_count = 1,batch_size
        call data_fftbox_dealloc(fftbox_batch_imag(batch_count))
        call data_fftbox_dealloc(fftbox_batch_real(batch_count))
    end do

    deallocate(fftbox_batch_imag,stat=ierr)
    call utils_dealloc_check('integrals_exr','fftbox_batch_imag',ierr)
    deallocate(fftbox_batch_real,stat=ierr)
    call utils_dealloc_check('integrals_exr','fftbox_batch_real',ierr)

    ! deallocate index
    deallocate(idx,stat=ierr)
    call utils_dealloc_check('integrals_exr','idx',ierr)

  end subroutine integrals_exr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !=======================================================================!
    ! This subroutine calculates matrix elements of powers of the position  !
    ! operator (within the FFTbox) between two function sets.               !
    ! Gives rmat(i)_\alpha\beta = < bra_\alpha | r^FFT_i^pow | ket_\beta >  !
    ! where r^FFT_i = r_i - R_FFT, with R_FFT being the origin of the beta  !
    ! function's FFTbox.                                                    !
    ! This function does not account for PAW augmentation charges and the   !
    ! result must be post-processed to include them.                        !
    !=======================================================================!
    !  Arguments:                                                           !
    !    rmat(3)    (inout) : sparse position operator matrix elements      !
    !    bras_on_grid  (in) : bra functions on grid in PPD format           !
    !    bra_basis  (input) : The basis type for the bras                   !
    !    kets_on_grid  (in) : ket functions on grid in PPD format           !
    !    ket_basis  (input) : The basis type for the kets                   !
    !    order      (input) : Power of the position operator r_i            !
    !    axis   (opt input) : The axis to calculate, if only doing one.     !
    !=======================================================================!
    ! Originally written by Mark Robinson as internal_matrix_elements in    !
    ! the routine polarisation_calculate.                                   !
    ! Re-worked by Nicholas Hine in 2008 to use SPAM3 and new version of    !
    ! function_ops_brappd_ketfftbox                                         !
    ! Moved to integrals_mod by Nicholas Hine on 11/11/2009.                !
    ! Changed from pattern ALTERNATE to pattern FULL, as with the former it !
    ! was impossible to correctly account for the origin shift (12/04/2010) !
    ! Modified by Andrea Greco on 21/06/2015 to alow use of complex NGWFs.  !
    !=======================================================================!

  subroutine integrals_pos(rmat, bras_on_grid, bra_basis, kets_on_grid, &
       ket_basis,dbl_grid,cell,fftbox,order,axis)

    use basis, only: basis_copy_function_to_box
    use cell_grid, only: GRID_INFO
    use comms, only : comms_reduce
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, &
         data_fftbox_alloc, data_fftbox_dealloc, data_set_to_zero
    use constants, only : DP
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_interpolate, fourier_filter
    use function_basis, only: FUNC_BASIS
    ! agrecocmplx
    use function_ops, only: function_ops_batch_col_start, &
         function_ops_brappd_ketfftbox
    use geometry, only: POINT, operator(*), operator(+) !ddor
    use rundat, only: pub_fftbox_batch_size
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index !ddor

    ! Arguments
    type(SPAM3), intent(inout)  :: rmat(3)
    type(FUNC_BASIS), intent(in) :: bra_basis
    !real(kind=DP), intent(in) :: bras_on_grid(bra_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: bras_on_grid
    type(FUNC_BASIS), intent(in) :: ket_basis
    !real(kind=DP), intent(in) :: kets_on_grid(ket_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: kets_on_grid
    type(GRID_INFO), intent(in) :: dbl_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    ! ddor: Power of position operator for which matrix elements are calculated
    integer, intent(in) :: order
    ! ddor: Specifies if we only need one direction
    integer, optional, intent(in)        :: axis

    ! Local Variables
    type(POINT) :: r1,r2,r3
    type(POINT) :: a1,a2,a3
    integer :: i1,i2,i3
    integer :: n3,ld1,ld2,ld1_dbl,ld2_dbl
    integer :: max_current_size
    integer :: local_start,local_end
    integer :: batch,n_batches,batch_size,batch_count
    integer, allocatable :: ket_box_start(:,:)
    integer, allocatable :: ket_start_in_box(:,:)
    integer :: ngwf,f1,f2
    integer :: idx_len
    integer :: xyz,ierr
    ! ddor: Integers used to select axis if only one is needed
    integer :: xyz_index,xyz_count

    ! workspace
    !real(kind=DP), allocatable, dimension(:,:,:,:) :: fftbox_batch
    ! agrecocmplx: use new FFTBOX_DATA type
    type(FFTBOX_DATA), allocatable, dimension(:) :: fftbox_batch
    integer,       allocatable, dimension(:)       :: idx

    ! position operator, FFTbox and tightbox
    real(kind=DP), dimension(:,:,:,:), allocatable :: r_op
    ! agrecocmplx: use new FFTBOX_DATA type
    !real(kind=DP), dimension(:,:,:), allocatable :: fftbox1,fftbox2
    !real(kind=DP), dimension(:,:,:), allocatable :: fftbox1_dbl,fftbox2_dbl
    type(FFTBOX_DATA) :: fftbox1,fftbox2
    type(FFTBOX_DATA) :: fftbox1_dbl,fftbox2_dbl
    complex(kind=DP), dimension(:,:,:), allocatable :: coarse_work,fine_work
    logical :: ongpu


    ! local copies of fftbox info
    n3 = fftbox%total_pt3
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2
    ld1_dbl = fftbox%total_ld1_dbl
    ld2_dbl = fftbox%total_ld2_dbl

    ! collect batch size and determine number of batches
    batch_size = pub_fftbox_batch_size
    n_batches = ket_basis%max_on_proc / batch_size
    if (mod(ket_basis%max_on_proc,batch_size) > 0) n_batches = n_batches + 1

    ! obtain index for rmat (all 3 are the same)
    idx_len = sparse_index_length(rmat(1))
    allocate(idx(idx_len),stat=ierr)
    call utils_alloc_check('integrals_pos','idx',ierr)
    call sparse_generate_index(idx,rmat(1))

    ! allocate workspace
    ! agrecocmplx: allocate using appropriate routines
    !allocate(fftbox_batch(ld1,ld2,n3,batch_size),stat=ierr)
    !call utils_alloc_check('integrals_pos','fftbox_batch',ierr)
    allocate(fftbox_batch(batch_size),stat=ierr)
    call utils_alloc_check('integrals_pos','fftbox_batch',ierr)
    do batch_count = 1,batch_size
        call data_fftbox_alloc(fftbox_batch(batch_count), &
            ld1,ld2,n3, iscmplx=kets_on_grid%iscmplx)
    end do

    ! allocate position operator and FFTboxes
    allocate(ket_box_start(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_pos','ket_box_start',ierr)
    allocate(ket_start_in_box(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_pos','ket_start_in_box',ierr)
    !allocate(fftbox1(ld1,ld2,n3), stat=ierr)
    !call utils_alloc_check('integrals_pos','fftbox1',ierr)
    !allocate(fftbox2(ld1,ld2,n3), stat=ierr)
    !call utils_alloc_check('integrals_pos','fftbox2', ierr)
    !allocate(fftbox1_dbl(ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    !call utils_alloc_check('integrals_pos','fftbox1_dbl',ierr)
    !allocate(fftbox2_dbl(ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    !call utils_alloc_check('integrals_pos','fftbox2_dbl',ierr)
    ! agrecocmplx: allocate using appropriate routines
    call data_fftbox_alloc(fftbox1, ld1,ld2,n3, &
         iscmplx=kets_on_grid%iscmplx)
    call data_fftbox_alloc(fftbox2, ld1,ld2,n3, &
         iscmplx=kets_on_grid%iscmplx)
    call data_fftbox_alloc(fftbox1_dbl, ld1_dbl,ld2_dbl,2*n3, &
         iscmplx=kets_on_grid%iscmplx)
    call data_fftbox_alloc(fftbox2_dbl, ld1_dbl,ld2_dbl,2*n3, &
         iscmplx=kets_on_grid%iscmplx)
    allocate(r_op(3,ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    call utils_alloc_check('integrals_pos','r_op',ierr)
    allocate(coarse_work(ld1,ld2,n3), stat=ierr)
    call utils_alloc_check('integrals_pos','coarse_work',ierr)
    allocate(fine_work(ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    call utils_alloc_check('integrals_pos','fine_work',ierr)

    ! calculate vectors between fine grid points
    a1 = (1.0_DP / dbl_grid%n1) * cell%a1
    a2 = (1.0_DP / dbl_grid%n2) * cell%a2
    a3 = (1.0_DP / dbl_grid%n3) * cell%a3

    ! ddor: construct position operator to power `order' in FFTbox
    r_op = 0.0_DP
    do i3=1,2*n3
       r3 = real(i3-1,kind=DP) * a3
       do i2=1,ld2_dbl
          r2 = r3 + real(i2-1,kind=DP) * a2
          do i1=1,ld1_dbl
             r1 = r2 + real(i1-1,kind=DP) * a1
             r_op(1,i1,i2,i3) = (r1%X)**order
             r_op(2,i1,i2,i3) = (r1%Y)**order
             r_op(3,i1,i2,i3) = (r1%Z)**order
          enddo
       enddo
    enddo

    ! ddor: Do one direction only if axis is present
    if (present(axis)) then
       xyz_count = 1
    else
       xyz_count = 3
    endif

    ! loop over x,y,z (or just axis)
    do xyz_index=1,xyz_count

       ! ddor: xyz is the label on the axis direction
       if (present(axis)) then
          xyz = axis
       else
          xyz = xyz_index
       endif

       ! loop over batches of NGWFs on this proc to construct
       ! the matrix elements R_ab = < phi_a | r | phi_b >
       local_start = 1
       do batch=1,n_batches

          ! find last NGWF in batch
          local_end = min(local_start+batch_size-1,ket_basis%proc_num)
          max_current_size = local_end - local_start
          call comms_reduce('MAX',max_current_size)

          call function_ops_batch_col_start(ket_box_start,ket_start_in_box, &
               batch_size,local_start,local_end,fftbox,cell,ket_basis)

          ! place NGWFs in FFTboxes
          batch_count = 1

          ! agrecocmplx: one at a time if using complex NGWFs
          if (.not.kets_on_grid%iscmplx) then

              do ngwf=local_start,local_start+max_current_size,2
                 f1 = ngwf
                 f2 = ngwf +1

                 ! copy function 1 into first fftbox if required
                 if (f1 <= local_end) then
                    ! agrecocmplx: use modified routine in basis_new
                    call basis_copy_function_to_box(fftbox1, &
                         ket_start_in_box(1,batch_count), &
                         ket_start_in_box(2,batch_count), &
                         ket_start_in_box(3,batch_count), &
                         ket_basis%tight_boxes(f1),kets_on_grid, &
                         ket_basis%spheres(f1), cell, fftbox)
                 else
                    fftbox1%d(:,:,:) = 0.0_DP
                 endif

                 ! copy function 2 into second fftbox if required
                 if (f2 <= local_end) then
                    ! agrecocmplx: use modified routine in basis_new
                    call basis_copy_function_to_box(fftbox2, &
                         ket_start_in_box(1,batch_count+1), &
                         ket_start_in_box(2,batch_count+1), &
                         ket_start_in_box(3,batch_count+1), &
                         ket_basis%tight_boxes(f2),kets_on_grid, &
                         ket_basis%spheres(f2),cell, fftbox)
                 else
                    fftbox2%d(:,:,:) = 0.0_DP
                 endif

                 if (f1 <= local_end) then
                    ! interpolate fftboxes to fine grid
                    call fourier_interpolate(coarse_work, fine_work, &
                         fftbox1%d,fftbox2%d,fftbox1_dbl%d,fftbox2_dbl%d)

                    ! apply FFTbox position operator and increment batch_count
                    fftbox1_dbl%d(:,:,:) = r_op(xyz,:,:,:) * fftbox1_dbl%d
                    if (f2 <= local_end) fftbox2_dbl%d(:,:,:) = r_op(xyz,:,:,:) &
                         * fftbox2_dbl%d

                    ! filter fftboxes to standard grid
                    call fourier_filter(coarse_work, fine_work, &
                         fftbox1_dbl%d,fftbox2_dbl%d,fftbox1%d,fftbox2%d)
                 endif

                 ! copy functions into batch
                 ! agrecocmplx
                 fftbox_batch(batch_count)%d(:,:,:) = fftbox1%d(:,:,:)
                 batch_count = batch_count +1
                 if (f2 <= local_end) then
                    fftbox_batch(batch_count)%d(:,:,:) = fftbox2%d(:,:,:)
                    batch_count = batch_count +1
                 endif
              enddo

          else

              do ngwf=local_start,local_start+max_current_size,1
                 f1 = ngwf

                 ! copy function 1 into first fftbox if required
                 if (f1 <= local_end) then
                     call basis_copy_function_to_box(fftbox1, &
                       ket_start_in_box(1,batch_count), &
                       ket_start_in_box(2,batch_count), &
                       ket_start_in_box(3,batch_count), &
                       ket_basis%tight_boxes(f1),kets_on_grid, &
                       ket_basis%spheres(f1), cell, fftbox)

                 else
                     call data_set_to_zero(fftbox1)
                 endif

                 if (f1 <= local_end) then
                     ! interpolate fftbox to fine grid
                     call fourier_interpolate(coarse_work, fine_work, &
                          fftbox1%z)

                     ! apply FFTbox position operator and increment batch_count
                     fftbox1_dbl%z(:,:,:) = &
                          r_op(xyz,:,:,:) * fftbox1_dbl%z(:,:,:)

                     ! filter fftbox to standard grid
                     call fourier_filter(coarse_work, fine_work, &
                          fftbox1_dbl%z, fftbox1%z)

                 end if

                 fftbox_batch(batch_count)%z(:,:,:) = fftbox1%z(:,:,:)
                 batch_count = batch_count +1

              end do

          end if

          ! ndmh: calculate position operator matrix elements in SPAM3 format
          ongpu=.false.
          call function_ops_brappd_ketfftbox(rmat(xyz),                   & ! inout
               bras_on_grid, bra_basis, cell, fftbox,                  & ! input
               fftbox_batch, ket_box_start, batch_size,                & ! input
               local_start, local_end, idx_len, idx, 'FULL',ongpu)

          ! increment batch start
          local_start = local_start + batch_size

       enddo  ! end batch loop

    enddo  ! end loop over xyz

    ! deallocate FFTbox and tightbox
    deallocate(fine_work,stat=ierr)
    call utils_dealloc_check('integrals_pos','fine_work',ierr)
    deallocate(coarse_work,stat=ierr)
    call utils_dealloc_check('integrals_pos','coarse_work',ierr)
    deallocate(r_op,stat=ierr)
    call utils_dealloc_check('integrals_pos','r_op',ierr)
    ! agrecocmplx
    call data_fftbox_dealloc(fftbox2_dbl)
    call data_fftbox_dealloc(fftbox1_dbl)
    call data_fftbox_dealloc(fftbox2)
    call data_fftbox_dealloc(fftbox1)
    deallocate(ket_start_in_box,stat=ierr)
    call utils_dealloc_check('integrals_pos','ket_start_in_box',ierr)
    deallocate(ket_box_start,stat=ierr)
    call utils_dealloc_check('integrals_pos','ket_box_start',ierr)

    ! deallocate workspace
    !deallocate(fftbox_batch,stat=ierr)
    !call utils_dealloc_check('integrals_pos','fftbox_batch',ierr)
    ! agrecocmplx: deallocate using appropriate routine
    do batch_count = 1,batch_size
        call data_fftbox_dealloc(fftbox_batch(batch_count))
    end do

    deallocate(fftbox_batch,stat=ierr)
    call utils_dealloc_check('integrals_pos','fftbox_batch',ierr)

    ! deallocate index
    deallocate(idx,stat=ierr)
    call utils_dealloc_check('integrals_pos','idx',ierr)

  end subroutine integrals_pos




  real(kind=DP) function integrals_trace_on_grid(integrand,grid)

    !======================================================================!
    ! Returns the integral of a function on the passed grid.               !
    !----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000. Modified by Peter D. Haynes!
    ! on 1/7/2004 for Fourier parallelisation.                             !
    ! Modified to use cell by Quintin Hill on 15/10/2008.                  !
    !======================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP
!$  use rundat, only: pub_threads_max

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(in) :: integrand(grid%ld1, &
         grid%ld2, grid%max_slabs12)

    ! Local variables
    integer :: i1,i2,islab12,ipt
    real(kind=DP) :: integral

    ! Sum over real-space on this proc
    integral = 0.0_DP
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(ipt,i1,i2,islab12) &
!$OMP SHARED(grid,integrand,pub_threads_max) REDUCTION(+:integral)
    do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
       i1 = modulo(ipt-1,grid%n1) + 1
       i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
       islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1

       integral = integral + integrand(i1,i2,islab12)
    end do
!$OMP END PARALLEL DO

    ! Sum over all procs
    call comms_reduce('SUM',integral)

    ! Multiply by weight per grid point
    integrals_trace_on_grid = integral * grid%weight

  end function integrals_trace_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function integrals_product_on_grid(grid, x1, x2, m1, m2, m3)
    !======================================================================!
    ! Returns the integral of a product of two functions, x1 and x2 on the !
    ! passed grid. If x2 is omitted, integrates x1 only. Optional arguments!
    ! m1, m2 and m3 can be used to narrow down the integration to a subset !
    ! of the grid.                                                         !
    !----------------------------------------------------------------------!
    ! Arguments:                                                           !
    !   grid (input): The grid on which the integration takes place.       !
    !   x1, x2 (input): The quantities to integrate. x2 is optional.       !
    !   m1, m2, m3 (input, optional): Values to override grid%n1, grid%n2  !
    !                                 and grid%num_my_slabs12 with.        !
    ! Return value: Integral of x1(r)*x2(r) on the grid.                   !
    !----------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in June 2010.                              !
    !======================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP
    use utils, only: utils_assert

    ! jd: Arguments
    type(GRID_INFO), intent(in)         :: grid
    real(kind=DP), intent(in)           :: x1(grid%ld1, grid%ld2, &
         grid%max_slabs12)
    real(kind=DP), intent(in), optional :: x2(grid%ld1, grid%ld2, &
         grid%max_slabs12)
    integer, intent(in), optional :: m1, m2, m3

    ! jd: Internal variables
    integer :: i1, i2, i3
    real(kind=DP) :: integral
    integer :: p1, p2, p3

    !------------------------------------------------------------------------

    ! jd: Default upper bounds of integration
    p1 = grid%n1
    p2 = grid%n2
    p3 = grid%num_my_slabs12

    ! jd: Override them, if the user desires
    if(present(m1)) p1 = m1
    if(present(m2)) p2 = m2
    if(present(m3)) p3 = m3

    call utils_assert(p1>0 .and. p2>0 .and. p3>=0,'Bad integration bounds in &
         &integrals_product_on_grid')

    integral = 0.0_DP
    ! jd: The following duplicates code in order not to check the
    !     if(present) condition repeatedly in a tight loop
    if(present(x2)) then
       do i3=1, p3
          do i2=1, p2
             do i1=1, p1
                integral = integral + x1(i1,i2,i3) * x2(i1,i2,i3)
             end do
          end do
       end do
    else
       do i3=1, p3
          do i2=1, p2
             do i1=1, p1
                integral = integral + x1(i1,i2,i3)
             end do
          end do
       end do
    end if

    call comms_reduce('SUM',integral)

    ! jd: NB, the weight remains unchanged even if integration bounds overridden
    integrals_product_on_grid = integral * grid%weight

  end function integrals_product_on_grid
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function integrals_product_on_grid_vec(grid, x1, x2, m1, m2, m3)
    !======================================================================!
    ! Returns the integral of a product of two functions, x1 and x2 on the !
    ! passed grid. x1 is assumed to have three Cartesian components. x2 is !
    ! a scalar function. The integral is thus a Cartesian vector.          !
    ! If x2 is omitted, integrates x1 only. Optional arguments m1, m2 and  !
    ! m3 can be used to narrow down the integration to a subset of the     !
    ! grid.                                                                !
    !----------------------------------------------------------------------!
    ! The rationale for not just using integrals_product_on_grid for every !
    ! component is cache-friendliness and compactness of notation in the   !
    ! caller.                                                              !
    !----------------------------------------------------------------------!
    ! Arguments:                                                           !
    !   grid (input): The grid on which the integration takes place.       !
    !   x1, x2 (input): The quantities to integrate. x2 is optional.       !
    !                   x1 is a vector field, x2 is a scalar field.        !
    !   m1, m2, m3 (input, optional): Values to override grid%n1, grid%n2  !
    !                                 and grid%num_my_slabs12 with.        !
    ! Return value: Integral of x1(r)*x2(r) on the grid.                   !
    !               (a Cartesian vector).                                  !
    !----------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2015.                         !
    !======================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP
    use utils, only: utils_assert

    ! jd: Arguments
    type(GRID_INFO), intent(in)         :: grid
    real(kind=DP), intent(in)           :: x1(3, grid%ld1, grid%ld2, &
         grid%max_slabs12)
    real(kind=DP), intent(in), optional :: x2(grid%ld1, grid%ld2, &
         grid%max_slabs12)
    integer, intent(in), optional :: m1, m2, m3

    ! jd: Return value
    real(kind=DP), dimension(3) :: integrals_product_on_grid_vec

    ! jd: Internal variables
    integer :: i1, i2, i3
    real(kind=DP) :: integral(3)
    integer :: p1, p2, p3

    !------------------------------------------------------------------------

    ! jd: Default upper bounds of integration
    p1 = grid%n1
    p2 = grid%n2
    p3 = grid%num_my_slabs12

    ! jd: Override them, if the user desires
    if(present(m1)) p1 = m1
    if(present(m2)) p2 = m2
    if(present(m3)) p3 = m3

    call utils_assert(p1>0 .and. p2>0 .and. p3>=0,'Bad integration bounds in &
         &integrals_product_on_grid_vec')

    integral(1:3) = 0.0_DP
    ! jd: The following duplicates code in order not to check the
    !     if(present) condition repeatedly in a tight loop
    if(present(x2)) then
       do i3=1, p3
          do i2=1, p2
             do i1=1, p1
                integral(1:3) = integral(1:3) + x1(1:3,i1,i2,i3) * x2(i1,i2,i3)
             end do
          end do
       end do
    else
       do i3=1, p3
          do i2=1, p2
             do i1=1, p1
                integral(1:3) = integral(1:3) + x1(1:3,i1,i2,i3) * x2(i1,i2,i3)
             end do
          end do
       end do
    end if

    call comms_reduce('SUM',integral(1:3))

    ! jd: NB, the weight remains unchanged even if integration bounds overridden
    integrals_product_on_grid_vec(1:3) = integral(1:3) * grid%weight

  end function integrals_product_on_grid_vec
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


end module integrals

