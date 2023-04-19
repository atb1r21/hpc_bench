! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   The subroutines in this file were created by                              !
!                                                                             !
!   Nicholas D.M. Hine                                                        !
!                                                                             !
!   in March 2015.                                                            !
!                                                                             !
!   based in large part on earlier versions in function_basis_mod.F90.        !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module function_ops

  use basis, only: SPHERE
  use constants, only: DP
  use datatypes, only: FUNCTIONS
  use function_basis, only: FUNC_BASIS

  implicit none

  private

  ! ndmh: This structure contains information required for parallel
  ! ndmh: communication of a set of functions
  type FUNC_COMM

    integer :: plan_steps  ! number of steps in plan
    integer :: nbuf        ! number of buffers
    integer :: rbuf        ! next buffer to receive to
    integer :: ubuf        ! next buffer to use

    ! ndmh: buffer to store requested indices
    integer :: probe_count

    integer, allocatable :: req_buffer(:)
    logical, allocatable :: reqs_out(:)
    logical, allocatable :: reqs_in(:)
    integer, allocatable :: plan(:,:)    ! Plan for this proc
    integer, allocatable :: func_requests(:)
    integer, allocatable :: handles(:)

    type(SPHERE), allocatable :: recv_spheres(:)
    !real(kind=DP), allocatable :: recv_func_on_grid(:,:)
    ! agrecocmplx: use new type FUNCTIONS
    type(FUNCTIONS), allocatable :: recv_func_on_grid(:)

  end type FUNC_COMM

  ! ndmh: special function request values and tags
  integer, parameter :: FUNCS_WAITING = -1000
  integer, parameter :: FUNCS_DONE = -2000
  integer, parameter :: req_tag = 10000000
  integer, parameter :: probe_frequency = 4
  integer, parameter :: num_buffers = 1

  ! ndmh: public subroutines

  ! ndmh: row sums / integrals routines
  public :: function_ops_batch_col_start
  public :: function_ops_sum_fftbox_batch
  public :: function_ops_sum_ppd_funcs
  public :: function_ops_brappd_ketfftbox
  public :: function_ops_brappd_ketppd

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine function_ops_sum_fftbox_batch(fftbox_batch, fftbox, cell, &
       funcs_on_grid, row_basis, col_box_start, batch_size, &
       local_start, local_end, mat_idx, idx_len, coeff_mat, mmat, nmat, &
       common_fac)

    !=========================================================================!
    ! This subroutine calculates sums in fftboxes of functions multiplied by  !
    ! coefficients in the form of matrices. The fftbox data is A_\alpha(r)    !
    ! corresponding to a given column function g_\alpha (where the set of     !
    ! column functions is not necessarily the same as the row functions).     !
    ! The fftbox data A_\alpha(r) is constructed from the sum of functions    !
    ! f_\beta multiplied by coefficients M_\beta\alpha, for functions f_\beta !
    ! which overlap \alpha (ie have a nonzero matrix element S_{\beta\alpha}  !
    ! The sum is given by                                                     !
    ! A_\alpha(r) = \sum_beta M_{\beta\alpha} f_\beta(r)                      !
    ! The batch contains a subset of the functions f_\alpha that belong to    !
    ! the proc pub_my_proc_id. Increasing the size of the batch               !
    ! increases parallel efficiency by reducing communication but also        !
    ! increases the memory use per processor.                                 !
    ! Both A_\alpha and M may be arrays, running from mmat to nmat - the sum  !
    ! will be taken over this range of arrays.                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! fftbox_batch (input/output) : FFTboxes containing results A(r)          !
    ! funcs_on_grid (input) : funcs of pub_my_proc_id in ppd representation   !
    ! row_basis (input) : Function basis type for the functions being summed  !
    ! col_box_start (input): Cell start pos of each column function FFTbox    !
    ! batch_size (input) : Number of fftboxes (phi_a functions) in each batch !
    ! local_start (input): Number of the first function of the current batch  !
    !    in the counting scheme of all functions of pub_my_proc_id.           !
    ! local_end (input)  : Number of the last function of the current batch   !
    !    in the counting scheme of all functions of pub_my_proc_id.           !
    ! mat_idx (input) : Overlap matrix index between row funcs and col funcs  !
    ! idx_len (input) : Length of overlap matrix index                        !
    ! coeff_mat (input) : Coefficients matrix M_{\alpha\beta} in SPAM3 format.!
    ! cell (input)           : The CELL_INFO type describing the cell         !
    ! fftbox (input)         : The FFTBOX_INFO type describing the col fftbox !
    ! mmat (input) : First entry in the array of A's and M's to use.          !
    ! nmat (input) : Last entry in the array of A's and M's to use.           !
    ! common_fac (input) : Common factor to multiply each function by         !
    !-------------------------------------------------------------------------!
    ! This subroutine was written by Chris-Kriton Skylaris on 18/9/2003       !
    ! and is capable of running on parallel computers with an                 !
    ! arbitary number of processors.                                          !
    ! Rewritten for SPAM3, moved to function_basis and made able to support   !
    ! different function sets by Nicholas Hine, April-May 2009.               !
    ! Modified to remove integer buffers and consolidate MPI usage            !
    ! by Nicholas Hine, June 2009.                                            !
    ! Modified to allow different function basis for rows and columns by      !
    ! David O'Regan, September 2009.                                          !
    ! Modified to use FUNC_COMM type by Nicholas Hine on 15/06/2012.          !
    ! Moved to function_ops_mod by Nicholas Hine on 03/03/2015.               !
    ! Modified by Andrea Greco on 30/04/2015 to make it compatible with       !
    ! the use of complex-valued NGWFs                                         !
    !=========================================================================!

    use basis, only: SPHERE, basis_add_function_to_box, &
         basis_find_function_wrt_box
    use comms, only: comms_barrier, comms_free, pub_my_proc_id
    use datatypes, only: COEF, FFTBOX_DATA, FUNCTIONS
    use fft_box, only: FFTBOX_INFO
    use rundat, only: pub_num_spins
!$  use rundat, only: pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_get_element, sparse_proc_of_elem, &
         sparse_first_elem_on_proc, sparse_proc_num_element
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    integer, intent(in) :: batch_size
    integer, intent(in) :: nmat
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    !real(kind=DP), intent(inout) :: fftbox_batch(fftbox%total_ld1, &
    !     fftbox%total_ld2, fftbox%total_pt3, pub_num_spins, &
    !     nmat, batch_size)
    ! agrecocmplx: use new type FFTBOX_DATA
    type(FFTBOX_DATA), intent(inout) :: fftbox_batch(pub_num_spins, &
          nmat, batch_size)
    integer, intent(in) :: col_box_start(3, batch_size)
    integer, intent(in) :: local_start, local_end
    integer, intent(in) :: mmat ! ddor: Starting index
    integer, intent(in) :: idx_len
    integer, intent(in) :: mat_idx(idx_len)
    type(FUNC_BASIS), intent(in) :: row_basis
    type(SPAM3), intent(in) :: coeff_mat(pub_num_spins,nmat)
    !real(kind=DP), intent(in)  :: funcs_on_grid(row_basis%size_on_grid)
    ! agrecocmplx: use new type FUNCTIONS
    type(FUNCTIONS), intent(in) :: funcs_on_grid
    real(kind=DP), intent(in) :: common_fac

    ! Local Variables
    integer :: local_row, local_col
    integer :: col
    integer :: recv_row
    integer :: first_row_on_proc
    integer :: recv_proc
    integer :: batch_count
    integer :: row_start1, row_start2, row_start3
    integer :: is
    integer :: imat
    logical :: remote_row, remote_req_row
    !real(kind=DP) :: coeff(pub_num_spins,nmat)
    ! agrecocmplx: use new type COEF
    type(COEF) :: coeff(pub_num_spins,nmat)
    ! ndmh: variables for planned execution
    integer :: req_row
    integer :: plan_step                 ! counter for current step in plan
    integer :: prev_req_row              ! Last row function requested
    integer :: prev_recv_row             ! Last row function received
    integer :: req_proc
    integer :: last_step
    integer :: loop_step
!$  integer, external :: omp_get_thread_num, omp_get_num_threads

    type(FUNC_COMM) :: row_comm

#ifdef ITC_TRACE
    call VTBEGIN(vt_function_ops_sum_fftbox_batch,vt_err)
#endif

    ! ndmh: allocate comms workspace
    row_comm%plan_steps = sparse_proc_num_element(coeff_mat(1,mmat))
    call function_ops_func_comm_init(row_comm,row_basis,2,funcs_on_grid%iscmplx)

    ! ndmh: create a plan from the matrix index
    call function_ops_batch_row_plan(row_comm,idx_len,mat_idx, &
         coeff_mat(1,mmat),local_start,local_end,'FULL')

    ! ndmh: initialisation of send request receive operations
    call function_ops_init_requests(row_comm)

    ! ndmh: initializations
    first_row_on_proc = sparse_first_elem_on_proc( &
         pub_my_proc_id,coeff_mat(1,mmat),'R')

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(loop_step,recv_row,col,imat,coeff,local_col,batch_count, &
!$OMP      row_start1,row_start2,row_start3,plan_step,last_step, &
!$OMP      recv_proc,prev_recv_row,req_row,req_proc,prev_req_row, &
!$OMP      remote_row,remote_req_row,local_row) &
!$OMP SHARED(row_comm,mmat,nmat,coeff_mat,common_fac, &
!$OMP      pub_my_proc_id,local_start,row_basis,cell, &
!$OMP      col_box_start,funcs_on_grid,fftbox_batch, &
!$OMP      first_row_on_proc,fftbox,pub_threads_num_fftboxes,pub_num_spins)

    ! ndmh: initializations
    prev_recv_row = -2
    prev_req_row = -2
    local_row = -1
    remote_row = .false.
    remote_req_row = .false.
    recv_row = -1
    plan_step = -1

    ! ndmh: loop over the steps of the plan
    do

       plan_step = plan_step + 1

       if (plan_step>row_comm%plan_steps) exit

       ! ndmh: skip recv if still pre-start
       if (plan_step > 0) then

          ! ndmh: find the row of the next matrix element to calculate
          recv_row = row_comm%plan(1,plan_step)

          ! ndmh: receive or copy the ppd list and offset into the buffer sphere
          if (prev_recv_row /= recv_row) then

             ! ndmh: find the proc the row function is stored on
             recv_proc = sparse_proc_of_elem(recv_row,coeff_mat(1,mmat),'R')

             ! ndmh: finish the recv operation for upcoming step
             if (recv_proc==pub_my_proc_id) then
                local_row = recv_row - first_row_on_proc + 1
                remote_row = .false.
             else
!$OMP MASTER
                call function_ops_finish_recv(recv_row,row_comm, &
                     row_basis,funcs_on_grid)
!$OMP END MASTER
!$OMP FLUSH(row_comm)
                remote_row = .true.
             end if
             prev_recv_row = recv_row

          end if

       end if

       ! ndmh: find last step in plan which uses current recv_row
       last_step = plan_step

       do
          if (last_step==row_comm%plan_steps) exit
          if (row_comm%plan(1,last_step+1)==recv_row) then
             last_step=last_step + 1
          else
             exit
          end if
       end do

      ! ndmh: check if we need to request a function to be sent to this proc
       ! ndmh: for the next step
       if (last_step+1<=row_comm%plan_steps) then
          req_row = row_comm%plan(1,last_step+1)
          req_proc = sparse_proc_of_elem(req_row,coeff_mat(1,mmat),'R')
          ! ndmh: if this function is not local to this proc and we do not
          ! ndmh: already have it, then send a request for it to req_proc
          if (req_row /= prev_req_row) then
             if (req_proc /= pub_my_proc_id) then
                remote_req_row = .true.
!$OMP MASTER
                call function_ops_start_recv(req_proc,req_row,row_comm, &
                     row_basis)
                call function_ops_request(req_proc,req_row,row_comm)
!$OMP END MASTER
             else
                remote_req_row = .false.
             end if
             prev_req_row = req_row
          end if
       end if

!$OMP MASTER
       ! ndmh: respond to any send requests made of this proc
       call function_ops_respond_to_reqs(row_comm,row_basis,funcs_on_grid)
!$OMP END MASTER

       if (plan_step < 1) cycle

       if (remote_row) then
          ! synchronise threads, so that master has finished recv operation
!$OMP BARRIER
       end if

!!$OMP DO
       do loop_step=plan_step,last_step

          ! ndmh: find relevant row and column information for this pair
          recv_row = row_comm%plan(1,loop_step)
          col = row_comm%plan(2,loop_step)
          local_col = col - sparse_first_elem_on_proc(pub_my_proc_id, &
               coeff_mat(1,mmat),'C') + 1
          batch_count = local_col - local_start + 1
!$        if (modulo(batch_count,omp_get_num_threads())/=omp_get_thread_num()) cycle

          ! ndmh: find coefficients multiplying this function for each fftbox
          do imat=mmat,nmat
             do is=1,pub_num_spins
                ! agrecocmplx: coeff is complex if coeff_mat is complex
                coeff(is,imat)%iscmplx = coeff_mat(is,imat)%iscmplx
                if (coeff(is,imat)%iscmplx) then
                   call sparse_get_element(coeff(is,imat)%z,coeff_mat(is,imat), &
                        recv_row,col)
                   coeff(is,imat)%z = coeff(is,imat)%z * common_fac
                else
                   call sparse_get_element(coeff(is,imat)%d,coeff_mat(is,imat), &
                        recv_row,col)
                   coeff(is,imat)%d = coeff(is,imat)%d * common_fac
                end if
             end do
          end do

          ! Find position of bra tightbox start wrt ket fftbox
          call basis_find_function_wrt_box(row_start1, row_start2, row_start3, &
               col_box_start(1,batch_count), col_box_start(2,batch_count), &
               col_box_start(3,batch_count), row_basis%all_tbs(recv_row), &
               cell, fftbox)

          ! ndmh: add function times coeff into the fftboxes of the batch
          do imat=mmat,nmat
             do is=1,pub_num_spins
                if (remote_row) then
                   ! agrecocmplx: call new function in basis_new_mod
                   call basis_add_function_to_box( &
                        fftbox_batch(is,imat,batch_count), &
                        row_start1,row_start2,row_start3, &
                        row_basis%all_tbs(recv_row), &
                        row_comm%recv_func_on_grid(row_comm%ubuf), &
                        row_comm%recv_spheres(row_comm%ubuf),coeff(is,imat), &
                        cell, fftbox)
                else
                   call basis_add_function_to_box( &
                        fftbox_batch(is,imat,batch_count), &
                        row_start1,row_start2,row_start3, &
                        row_basis%all_tbs(recv_row), &
                        funcs_on_grid, &
                        row_basis%spheres(local_row),coeff(is,imat), &
                        cell, fftbox)
                end if
             end do  ! is
          end do  ! imat

!$OMP MASTER
          ! ndmh: respond to any send requests made of this proc
          call function_ops_respond_to_reqs(row_comm,row_basis,funcs_on_grid)
!$OMP END MASTER

       end do
!!$OMP END DO NOWAIT

       plan_step = last_step

       if (remote_req_row) then
!$OMP BARRIER
       end if

    end do  ! plan_step
!$OMP END PARALLEL

    call function_ops_finalise(row_comm,row_basis,funcs_on_grid) !jmecmplx

    ! ndmh: synchronise again so that non-blocking sends are completed and all
    ! ndmh: handles restored before deallocation
    call comms_free
    call comms_barrier

    call function_ops_func_comm_exit(row_comm)

#ifdef ITC_TRACE
    call VTEND(vt_function_ops_sum_fftbox_batch,vt_err)
#endif

  end subroutine function_ops_sum_fftbox_batch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_ops_sum_ppd_funcs(sum_on_grid, &        ! inout
       sum_basis, coeff_mat, mmat, nmat, index_mat, &           ! input
       funcs_on_grid, fun_basis)                                ! input

    !=========================================================================!
    ! This subroutine calculates sums in ppds of functions in ppds multiplied !
    ! by coefficients in the form of SPAM3 matrices, where only the functions !
    ! for which there is a nonzero element of an overlap matrix are included  !
    ! in the summation.                                                       !
    ! The sum is given by                                                     !
    ! |sum_a> = \sum_b M_ab |f_b>     for all b where S_ab /= 0               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! sum_on_grid (inout) : sum funcs of pub_my_proc_id in ppd representation !
    ! sum_basis (input) : Function basis type for the result functions        !
    ! funcs_on_grid (input) : funcs of pub_my_proc_id in ppd representation   !
    ! fun_basis (input) : Function basis type for the functions being summed  !
    ! coeff_mat (input) : coefficients matrix M_{\alpha\beta} in SPAM3 format.!
    ! mmat (input) : First entry in the array of M's to use.                  !
    ! nmat (input) : Last entry in the array of M's to use.                   !
    ! index_mat (input) : index matrix S_{\alpha\beta} in SPAM3 format.       !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine in May 2010, reusing bits of the routine       !
    ! integrals_brappd_ketppd (originally integrals_brackets).                !
    ! Modified to use FUNC_COMM type by Nicholas Hine on 15/06/2012.          !
    ! Moved to function_ops_mod by Nicholas Hine on 03/03/2015.               !
    ! Modified by Andrea Greco on 30/04/2015 to make it compatible with the   !
    ! use of complex valued NGWFs                                             !
    !=========================================================================!

    use datatypes, only: COEF, FUNCTIONS
    use comms, only: comms_barrier, comms_free, pub_my_proc_id
!$  use rundat, only: pub_threads_max
    use sparse, only: SPAM3, sparse_get_element, sparse_proc_of_elem, &
         sparse_first_elem_on_proc, sparse_num_elems_on_proc, &
         sparse_proc_num_element, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    integer, intent(in) :: mmat, nmat
    type(FUNC_BASIS), intent(in) :: fun_basis
    type(FUNC_BASIS), intent(in) :: sum_basis
    !real(kind=DP), intent(inout) :: sum_on_grid(sum_basis%size_on_grid,nmat)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(inout) :: sum_on_grid(nmat)
    type(SPAM3), intent(in) :: coeff_mat(nmat)
    type(SPAM3), intent(in) :: index_mat
    !real(kind=DP), intent(in) :: funcs_on_grid(fun_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: funcs_on_grid

    ! Local Variables
    integer :: col, recv_row
    integer :: local_col, local_row
    integer :: first_row_on_proc
    integer :: recv_proc
    logical :: remote_row, remote_req_row
    !real(kind=DP) :: coeff_mat_el      ! matrix element
    ! agrecocmplx: use new COEF type
    type(COEF) :: coeff_mat_el          ! matrix element
    integer :: ierr                     ! ndmh: error flag
    integer :: idx_len
    integer :: imat
    integer,allocatable :: index_mat_idx(:)
    ! ndmh: variables for planned execution
    integer :: req_row
    integer :: plan_step                ! counter for current step in plan
    integer :: prev_req_row             ! Last row function requested
    integer :: prev_recv_row            ! Last row function received
    integer :: req_proc
    integer :: last_step
    integer :: loop_step
    type(FUNC_COMM) :: row_comm
!$  integer, external :: omp_get_thread_num, omp_get_num_threads

    call timer_clock('function_ops_sum_ppd_funcs',1)

#ifdef ITC_TRACE
    call VTBEGIN(vt_function_ops_sum_ppd_funcs,vt_err)
#endif

    ! cks: synchronize PEs
    call comms_barrier

    ! ndmh: initialise comms workspace
    row_comm%plan_steps = sparse_proc_num_element(index_mat)
    call function_ops_func_comm_init(row_comm,fun_basis,2,funcs_on_grid%iscmplx)

    ! ndmh: get the matrix index
    idx_len = sparse_index_length(index_mat)
    allocate(index_mat_idx(idx_len),stat=ierr)
    call utils_alloc_check('function_ops_sum_ppd_funcs','index_mat_idx',ierr)
    call sparse_generate_index(index_mat_idx,index_mat)

    ! ndmh: create a plan from the matrix index
    local_col = sparse_num_elems_on_proc(pub_my_proc_id,coeff_mat(mmat),'C')
    call function_ops_batch_row_plan(row_comm,idx_len,index_mat_idx, &
         coeff_mat(mmat),1,local_col,'FULL')

    deallocate(index_mat_idx,stat=ierr)
    call utils_dealloc_check('function_ops_sum_ppd_funcs','index_mat_idx', &
         ierr)

    ! ndmh: initialisation of send request receive operations
    call function_ops_init_requests(row_comm)

    first_row_on_proc = sparse_first_elem_on_proc( &
          pub_my_proc_id,coeff_mat(mmat),'R')

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(plan_step,last_step,loop_step,col,local_col,imat,coeff_mat_el, &
!$OMP      prev_recv_row,recv_row,local_row,remote_row,recv_proc,prev_req_row, &
!$OMP      req_row,req_proc,remote_req_row) &
!$OMP SHARED(row_comm,pub_my_proc_id,coeff_mat, &
!$OMP      fun_basis,sum_basis,mmat,nmat,first_row_on_proc, &
!$OMP      sum_on_grid,funcs_on_grid,pub_threads_max)

    ! ndmh: initializations
    prev_recv_row = -2
    prev_req_row = -2
    local_row = -1
    remote_row = .false.
    remote_req_row = .false.
    plan_step = -1
    recv_row = -1

    ! ndmh: loop over the steps of the plan
    do

       plan_step = plan_step + 1
       if (plan_step>row_comm%plan_steps) exit

       ! ndmh: skip recv if still pre-start
       if (plan_step > 0) then

          ! ndmh: find the row of the next matrix element to calculate
          recv_row = row_comm%plan(1,plan_step)

          ! ndmh: receive or copy the ppd list and offset into the buffer sphere
          if (prev_recv_row /= recv_row) then

             ! ndmh: find the proc the row function is stored on
             recv_proc = sparse_proc_of_elem(recv_row,coeff_mat(mmat),'R')

             ! ndmh: finish the recv operation for upcoming step
             if (recv_proc==pub_my_proc_id) then
                local_row = recv_row - first_row_on_proc + 1
                remote_row = .false.
             else
!$OMP MASTER
                call function_ops_finish_recv(recv_row,row_comm,fun_basis, &
                     funcs_on_grid)
!$OMP END MASTER
!$OMP FLUSH(row_comm)
                remote_row = .true.
             end if
             prev_recv_row = recv_row

          end if

       end if

       ! ndmh: find last step in plan which uses current recv_row
       last_step = plan_step

       do
          if (last_step==row_comm%plan_steps) exit
          if (row_comm%plan(1,last_step+1)==recv_row) then
             last_step=last_step + 1
          else
             exit
          end if
       end do

       ! ndmh: check if we need to request a function to be sent to this proc
       ! ndmh: for the next step
       if (last_step+1<=row_comm%plan_steps) then
          req_row = row_comm%plan(1,last_step+1)
          req_proc = sparse_proc_of_elem(req_row,coeff_mat(mmat),'R')
          ! ndmh: if this bra is not local to this proc and we do not already
          ! ndmh: have it, then send a request for it to req_proc
          if (req_row /= prev_req_row) then
             if (req_proc /= pub_my_proc_id) then
                remote_req_row = .true.
!$OMP MASTER
                call function_ops_start_recv(req_proc,req_row,row_comm, &
                     fun_basis)
                call function_ops_request(req_proc,req_row,row_comm)
!$OMP END MASTER
             else
                remote_req_row = .false.
             end if
             prev_req_row = req_row
          end if
       end if

       ! ndmh: respond to any send requests made of this proc
!$OMP MASTER
       call function_ops_respond_to_reqs(row_comm,fun_basis,funcs_on_grid)
!$OMP END MASTER

       if (plan_step < 1) cycle

       if (remote_row) then
!$OMP BARRIER
       end if

!$OMP DO
       do loop_step=plan_step,last_step

          ! Find local index of col on this proc
          col = row_comm%plan(2,loop_step)
          local_col = col - &
               sparse_first_elem_on_proc(pub_my_proc_id,coeff_mat(mmat),'C') + 1

          do imat=mmat,nmat

             ! ndmh: get matrix element from SPAM3
             ! agrecocmplx: coeff_mat_el is complex if coeff_mat is complex
             coeff_mat_el%iscmplx = coeff_mat(imat)%iscmplx
             if (coeff_mat_el%iscmplx) then
                call sparse_get_element(coeff_mat_el%z,coeff_mat(imat), &
                     recv_row,col)
             else
                call sparse_get_element(coeff_mat_el%d,coeff_mat(imat), &
                     recv_row,col)
             end if

             ! ndmh: add coeff_mat matrix element times function to result
             if (remote_row) then
                call internal_add_ppds(sum_on_grid(imat), &
                     row_comm%recv_func_on_grid(row_comm%ubuf), sum_basis%n_ppds, &
                     row_comm%recv_spheres(row_comm%ubuf)%n_ppds_sphere, &
                     sum_basis%spheres(local_col),&
                     row_comm%recv_spheres(row_comm%ubuf), coeff_mat_el)
             else
                call internal_add_ppds(sum_on_grid(imat), &
                     funcs_on_grid, sum_basis%n_ppds, fun_basis%n_ppds, &
                     sum_basis%spheres(local_col),fun_basis%spheres(local_row), &
                     coeff_mat_el)
             end if

          end do

!$        if (omp_get_thread_num()==0) then
          ! ndmh: respond to any send requests made of this proc
          call function_ops_respond_to_reqs(row_comm,fun_basis,funcs_on_grid)
!$        end if

       end do
!$OMP END DO
       plan_step = last_step

       if (remote_req_row) then
!$OMP BARRIER
       end if

    end do  ! plan_step
!$OMP END PARALLEL

    ! ndmh: wait until completion messages have been received from all
    ! ndmh: other procs with which communication is required.
    call function_ops_finalise(row_comm,fun_basis,funcs_on_grid)

    ! ndmh: synchronise again so that non-blocking sends are completed and all
    ! ndmh: handles restored before deallocation
    call comms_free
    call comms_barrier

    ! ndmh: deallocate workspace
    call function_ops_func_comm_exit(row_comm)

    call timer_clock('function_ops_sum_ppd_funcs',2)

#ifdef ITC_TRACE
    call VTEND(vt_function_ops_sum_fftbox_batch,vt_err)
#endif

  contains

    ! ndmh: copy of integrals_bra_dot_ket_ppds, turned around logically
    subroutine internal_add_ppds(sumfunc, func, n_sum_ppds, n_func_ppds, &
         sum_sphere, func_sphere, coeff)

      use datatypes, only: data_functions_axpy
      implicit none

      ! Arguments
      integer, intent(in) :: n_sum_ppds, n_func_ppds
      !real(kind=DP), intent(inout) :: sumfunc(n_sum_ppds*sum_basis%n_pts)
      ! agrecocmplx: use new FUNCTIONS type
      type(FUNCTIONS), intent(inout) :: sumfunc
      !real(kind=DP), intent(in) :: func(n_func_ppds*fun_basis%n_pts)
      type(FUNCTIONS), intent(in) :: func
      type(SPHERE), intent(in)  :: sum_sphere
      type(SPHERE), intent(in)  :: func_sphere
      !real(kind=DP), intent(in) :: coeff
      ! agrecocmplx: use new COEF type
      type(COEF), intent(in) :: coeff

      ! Local Variables
      integer :: isum, ifunc
      integer :: sum_ppd, func_ppd
      integer :: sum_start, func_start
      !integer :: i                       ! loop counter

      ! Start at the beginning of the sum sphere ppd list
      ifunc = 1
      func_ppd = func_sphere%ppd_list(1,ifunc)

      ! ndmh: calculate sum = sum + alpha*func
      do isum=1,sum_sphere%n_ppds_sphere
         ! ndmh: find ppd number and start position for isum
         sum_ppd = sum_sphere%ppd_list(1,isum)
         sum_start = sum_sphere%offset + (isum-1)*sum_basis%n_pts
         do
            ! ndmh: keep moving on while func_ppd is less than sum_ppd
            if (func_ppd < sum_ppd) then
               ifunc = ifunc + 1
               if (ifunc > func_sphere%n_ppds_sphere) exit
               func_ppd = func_sphere%ppd_list(1,ifunc)
               cycle
            end if

            ! ndmh: ppd numbers match, so add this ppd of func to sum
            if (func_ppd == sum_ppd) then
               func_start = func_sphere%offset + (ifunc-1)*fun_basis%n_pts
               !do i=0,fun_basis%n_pts-1
               !   sumfunc(sum_start+i) = sumfunc(sum_start+i) + &
               !        coeff*func(func_start+i)
               !
               !end do

               ! agreco: perform axpy using new routine in basis_new_mod
               call data_functions_axpy(sumfunc, func, coeff, &
                    sum_start, func_start, fun_basis%n_pts)

            end if

            ! ndmh: move on to next sum ppd as soon as we have found a match
            ! ndmh: or moved beyond sum_ppd
            if (func_ppd >= sum_ppd) exit

         end do
      end do

    end subroutine internal_add_ppds

  end subroutine function_ops_sum_ppd_funcs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_ops_brappd_ketfftbox(bracket,               &    ! inout
       bras_on_grid, bra_basis, cell, fftbox,                     &    ! input
       ket_fftbox_batch, ket_box_start, batch_size,               &    ! input
       local_start, local_end, idx_len, bracket_idx, pattern,     &    ! input
       ongpu, factor_type)                                             ! input

    !=========================================================================!
    ! This subroutine computes a "bracket" <fa|fb> in sparse matrix (SPAM3)   !
    ! storage form. In general, the "bras" fa, which are represented on the   !
    ! coarse grid in ppd-indexed form, are a different set of functions from  !
    ! the "kets" fb, which are in fftboxes. The matrix is divided into blocks !
    ! calculated by different procs. In general the ket functions belong      !
    ! to the local proc while the bra functions are received from other procs.!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! bracket (input/output) : SPAM3 structure for matrix elements            !
    ! bracket_idx (input)    : index for above matrix                         !
    ! idx_len (input)        : Length of index for above matrix               !
    ! bras_on_grid (input)   : Bra functions in ppd representation            !
    ! ket_fftbox_batch (input) : Ket functions in fftbox representation       !
    ! ket_box_start (input)  : Start positions of ket functions in fftbox     !
    ! local_start (input)    : First ket in batch                             !
    ! local_end (input)      : Last ket in batch                              !
    ! batch_size (input)     : Size of batch                                  !
    ! bra_basis (input)      : The basis type for the bra functions           !
    ! cell (input)           : The CELL_INFO type describing the cell         !
    ! fftbox (input)         : The FFTBOX_INFO type describing the ket fftbox !
    ! pattern (input)        : Which terms of the matrix to calculate         !
    !                          pattern = 'FULL' calculates all terms          !
    !-------------------------------------------------------------------------!
    ! Integral subroutines were originally written by Chris-Kriton Skylaris   !
    ! in 2000.                                                                !
    ! Arash Mostofi rewrote them so that they use a "triple box" and          !
    ! complex-to-complex FFTs.                                                !
    ! This subroutine was written by Chris-Kriton Skylaris on 18/9/2003       !
    ! and is capable of running on parallel computers with an                 !
    ! arbitary number of processors.                                          !
    ! Modifications and speed-ups by Peter Haynes 16/03/2005.                 !
    ! Modifications and speed-ups by Chris-Kriton Skylaris 08/09/2005.        !
    ! Further modification to use parallel SPAM 2, July 2006                  !
    ! Modifications to symmetrise the locpot matrix consistently between      !
    ! upper and lower triangles by Nicholas Hine, 28/04/2008                  !
    ! Modified to use planned exection by Nicholas Hine, 14/05/2008.          !
    ! Modified for SPAM3, function basis type and request-based comms by      !
    ! Nicholas Hine in June 2009.                                             !
    ! Modified to use FUNC_COM by Nicholas Hine in on 15/06/2012.             !
    ! Moved to function_ops_mod by Nicholas Hine on 03/03/2015.               !
    ! Modified by Andrea Greco on 30/04/2015 to make it compatible with the   !
    ! use of complex-valued NGWFs                                             !
    !=========================================================================!
#ifdef GPU_PGI
    use fourier_gpu_wrapper_mod
    use accel_lib !! External dependency
    use cudafor !! External dependency
    use constants, only: DP
#endif

    use basis, only: basis_dot_function_with_box, basis_find_function_wrt_box
    use comms, only: comms_barrier, comms_free, pub_my_proc_id
    use datatypes, only: FFTBOX_DATA, FUNCTIONS, COEF, data_coef_scale
    use fft_box, only: FFTBOX_INFO
!$  use rundat, only: pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_put_element, sparse_proc_of_elem, &
         sparse_first_elem_on_proc, sparse_proc_num_element
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    integer, intent(in)        :: batch_size, local_start, local_end
    type(SPAM3), intent(inout) :: bracket
    type(FUNC_BASIS), intent(in) :: bra_basis
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    integer, intent(in)        :: idx_len
    integer, intent(in)        :: bracket_idx(idx_len)
    !real(kind=DP), intent(in)  :: bras_on_grid(bra_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: bras_on_grid
    !real(kind=DP), intent(in)  :: ket_fftbox_batch(fftbox%total_ld1, &
    !     fftbox%total_ld2, fftbox%total_pt3, batch_size)
    ! agrecocmplx: use new FFTBOX_DATA type
    type(FFTBOX_DATA), intent(in) :: ket_fftbox_batch(batch_size)
    integer, intent(in)        :: ket_box_start(3, batch_size)
    character(*), intent(in)   :: pattern
    integer, intent(in), optional :: factor_type ! if 1, 2, or 3, it implies
        ! that one does \sum_R R_i \int bra(r-R) ket(r) dr, where "i"
        ! is the factor type. In normal operation, there is no R_i ! gcc32

    ! internal declarations
    integer :: col, recv_row
    integer :: local_col, local_row
    integer :: first_row_on_proc
    integer :: recv_proc
    integer :: bra_start
    integer :: bra_start1, bra_start2, bra_start3
    integer :: batch_count
    integer :: dot_npts
    logical :: remote_row, remote_req_row
    logical, intent(in) :: ongpu ! kaw: Copies data back from the GPU
                                           !      If called from the GPU version
                                           !      of a routine.
    type(COEF) :: bracket_el ! matrix element

    ! ndmh: variables for planned execution
    integer :: req_row
    integer :: plan_step                 ! counter for current step in plan
    integer :: prev_req_row              ! Last row function requested
    integer :: prev_recv_row             ! Last row function received
    integer :: req_proc
    integer :: last_step
    integer :: loop_step
    type(FUNC_COMM) :: row_comm

!$  integer, external :: omp_get_thread_num, omp_get_num_threads

#ifdef ITC_TRACE
    call VTBEGIN(vt_integrals_brappd_ketfftbox,vt_err)
#endif

    ! cks: synchronize PEs
    call comms_barrier

    if (present(factor_type)) then
       if ((factor_type .gt. 3).or.(factor_type .lt. 1)) then
          call utils_abort('Wrong factor for brappd_ketppd')
       end if
    end if

    ! ndmh: allocate workspace for plan
    row_comm%plan_steps = sparse_proc_num_element(bracket)
    call function_ops_func_comm_init(row_comm,bra_basis,2,bras_on_grid%iscmplx)

    ! ndmh: create a plan from the matrix index
    call function_ops_batch_row_plan(row_comm,idx_len,bracket_idx, &
         bracket,local_start,local_end,pattern)

    ! ndmh: initialisation of send request receive operations
    call function_ops_init_requests(row_comm)

    first_row_on_proc = sparse_first_elem_on_proc(pub_my_proc_id,bracket,'R')

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(plan_step,last_step,loop_step,col,bracket_el,local_col, &
!$OMP      batch_count,bra_start1,bra_start2,bra_start3,dot_npts,bra_start, &
!$OMP      prev_recv_row,prev_req_row,local_row,remote_row,recv_row, &
!$OMP      recv_proc,req_proc,req_row,remote_req_row) &
!$OMP SHARED(row_comm,bracket,pub_my_proc_id,local_start,bra_basis,cell, &
!$OMP      ket_box_start,bras_on_grid,ket_fftbox_batch,fftbox, &
!$OMP      first_row_on_proc,pub_threads_num_fftboxes, factor_type)

    ! ndmh: initializations
    prev_recv_row = -2
    prev_req_row = -2
    local_row = -1
    remote_row = .false.
    remote_req_row = .false.
    plan_step = -1
    recv_row = -1


#ifdef GPU_PGI
       if (ongpu) then
          do batch_count = 1, batch_size
             ket_fftbox_batch(batch_count)%d(:,:,:) = potential_fftbox_batch_gpu(:,:,:,batch_count)
          end do
       end if
#endif

    ! ndmh: loop over the steps of the plan

    do
       plan_step = plan_step + 1
       if (plan_step>row_comm%plan_steps) exit

       ! ndmh: skip recv if still pre-start
       if (plan_step > 0) then

          ! ndmh: find the row and col of the matrix element to calculate
          recv_row = row_comm%plan(1,plan_step)

          ! ndmh: receive or copy the ppd list and offset into the buffer sphere
          if (prev_recv_row /= recv_row) then

             ! ndmh: find the proc the row function is stored on
             recv_proc = sparse_proc_of_elem(recv_row,bracket,'R')

             ! ndmh: finish the recv operation for upcoming step
             if (recv_proc==pub_my_proc_id) then
                local_row = recv_row - first_row_on_proc + 1
                remote_row = .false.
             else
!$OMP MASTER
                call function_ops_finish_recv(recv_row,row_comm,bra_basis, &
                     bras_on_grid)
!$OMP END MASTER
!$OMP FLUSH(row_comm)
                remote_row = .true.
             end if
             prev_recv_row = recv_row

          end if

       end if

       ! ndmh: find last step in plan which uses current recv_row
       last_step = plan_step
       do
          if (last_step==row_comm%plan_steps) exit
          if (row_comm%plan(1,last_step+1)==recv_row) then
             last_step=last_step + 1
          else
             exit
          end if
       end do

       ! ndmh: check if we need to request a function to be sent to this proc
       ! ndmh: for the next step
       if (last_step+1<=row_comm%plan_steps) then
          req_row = row_comm%plan(1,last_step+1)
          req_proc = sparse_proc_of_elem(req_row,bracket,'R')
          ! ndmh: if this function is not local to this proc and we do not already
          ! ndmh: have it, then send a request for it to req_proc
          if (req_row /= prev_req_row) then
             if (req_proc /= pub_my_proc_id) then
                remote_req_row = .true.
!$OMP MASTER
                call function_ops_start_recv(req_proc,req_row,row_comm, &
                     bra_basis)
                call function_ops_request(req_proc,req_row,row_comm)
!$OMP END MASTER
             else
                remote_req_row = .false.
             end if
             prev_req_row = req_row
          end if
       end if

       ! ndmh: respond to any send requests made of this proc
!$OMP MASTER
       call function_ops_respond_to_reqs(row_comm,bra_basis,bras_on_grid)
!$OMP END MASTER

       ! ndmh: cycle if still pre-start
       if (plan_step < 1) cycle

       if (remote_row) then
!$OMP BARRIER
       end if

!!$OMP DO
       do loop_step=plan_step,last_step

          ! ndmh: find relevant column information for this pair
          col = row_comm%plan(2,loop_step)
          local_col = col - &
               sparse_first_elem_on_proc(pub_my_proc_id,bracket,'C') + 1
          batch_count = local_col - local_start + 1
!$        if (modulo(batch_count,omp_get_num_threads())/=omp_get_thread_num()) cycle

          ! Find position of bra tightbox start wrt ket fftbox
          call basis_find_function_wrt_box(bra_start1, bra_start2, bra_start3, &
               ket_box_start(1,batch_count), ket_box_start(2,batch_count), &
               ket_box_start(3,batch_count), bra_basis%all_tbs(recv_row), cell, &
               fftbox)

          ! agrecocmplx: bracket_el complex if bracket is complex
          bracket_el%iscmplx = bracket%iscmplx

          ! cks: extract ppds belonging to bra function from ket fftbox
          ! cks: then ddot ppds - bra and ket representations have the same ppds
          if (remote_row) then
             bra_start = row_comm%recv_spheres(row_comm%ubuf)%offset
             dot_npts = row_comm%recv_spheres(row_comm%ubuf)%n_ppds_sphere * cell%n_pts
             call basis_dot_function_with_box(bracket_el, &
                  row_comm%recv_func_on_grid(row_comm%ubuf), &
                  ket_fftbox_batch(batch_count), &
                  row_comm%recv_spheres(row_comm%ubuf), &
                  bra_basis%all_tbs(recv_row), bra_start1, bra_start2, &
                  bra_start3, bra_start, cell, fftbox, factor_type, &
                  allow_mix_types=.true.)
          else
             bra_start = bra_basis%spheres(local_row)%offset
             dot_npts = bra_basis%spheres(local_row)%n_ppds_sphere * cell%n_pts
             call basis_dot_function_with_box(bracket_el, &
                  bras_on_grid, ket_fftbox_batch(batch_count), &
                  bra_basis%spheres(local_row), bra_basis%all_tbs(recv_row), &
                  bra_start1, bra_start2, bra_start3, bra_start, &
                  cell, fftbox, factor_type, allow_mix_types=.true.)
          end if

          ! cks: scale with grid point weight
          ! agreco: use function in basis_mod
          call data_coef_scale(bracket_el, cell%weight)

          ! ndmh: deposit in matrix
          if (bracket_el%iscmplx) then
              call sparse_put_element(bracket_el%z,bracket,recv_row,col)
          else
              call sparse_put_element(bracket_el%d,bracket,recv_row,col)
          end if

!$OMP MASTER
          ! ndmh: respond to any send requests made of this proc
          call function_ops_respond_to_reqs(row_comm,bra_basis,bras_on_grid)
!$OMP END MASTER

       end do
!!$OMP END DO NOWAIT

       if (remote_req_row) then
!$OMP BARRIER
       end if

       plan_step = last_step

    end do

!$OMP END PARALLEL

    call function_ops_finalise(row_comm,bra_basis,bras_on_grid)

    ! ndmh: synchronise again so that non-blocking sends are completed and all
    ! ndmh: handles restored before deallocation
    call comms_free
    call comms_barrier

    ! ndmh: deallocate comms workspace
    call function_ops_func_comm_exit(row_comm)

#ifdef ITC_TRACE
    call VTEND(vt_integrals_brappd_ketfftbox,vt_err)
#endif

  end subroutine function_ops_brappd_ketfftbox


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine function_ops_brappd_ketppd(bracket,  &                ! inout
       bras_on_grid, bra_basis, kets_on_grid, ket_basis, cell, &   ! input
       factor_type)     ! input

    !=========================================================================!
    ! This subroutine computes a "bracket" <fa|fb> in sparse matrix (SPAM3)   !
    ! storage form. In general, the "bras" fa can be a different set of       !
    ! functions from the "kets" fb, but both are represented on the same grid !
    ! in ppd-indexed form. The matrix is divided into blocks calculated by    !
    ! different procs. In general the ket functions belong to the local proc  !
    ! while the bra functions are received from other procs.                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! bracket (input/output) : SPAM3 structure for matrix elements            !
    ! bras_on_grid (input)   : Bra functions in ppd representation            !
    ! kets_on_grid (input)   : Ket functions in ppd representation            !
    ! bra_basis (input)      : The basis type for the bra functions           !
    ! ket_basis (input)      : The basis type for the ket functions           !
    ! cell (input)           : The CELL_INFO type describing the cell         !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine in July 2009 based on integrals_brackets,      !
    ! originally written by Chris-Kriton Skylaris in 2000 and modified by     !
    ! Oswaldo Dieguez, Arash Mostofi, Peter Haynes and Nicholas Hine in       !
    ! 2003 to 2009.                                                           !
    ! Modified to use FUNC_COM by Nicholas Hine on 15/06/2012.                !
    ! Moved to function_ops_mod by Nicholas Hine on 03/03/2015.               !
    ! Modified by Andrea Greco on 30/04/2015 to make it compatible with       !
    ! complex-valued NGWFs                                                    !
    !=========================================================================!

    use datatypes, only: COEF, FUNCTIONS, data_coef_scale
    use simulation_cell, only: CELL_INFO
    use comms, only: comms_barrier, comms_free, pub_my_proc_id
    use constants, only: stdout
!$  use rundat, only: pub_threads_max
    use sparse, only: SPAM3, sparse_put_element, sparse_proc_of_elem, &
         sparse_first_elem_on_proc, sparse_num_elems_on_proc, &
         sparse_proc_num_element, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: bracket
    type(FUNC_BASIS), intent(in) :: bra_basis
    type(FUNC_BASIS), intent(in) :: ket_basis
    type(CELL_INFO), intent(in) :: cell
    !real(kind=DP), intent(in) :: bras_on_grid(bra_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: bras_on_grid
    !real(kind=DP), intent(in) :: kets_on_grid(ket_basis%size_on_grid)
    type(FUNCTIONS), intent(in) :: kets_on_grid
    integer, intent(in), optional :: factor_type ! if 1, 2, or 3, it implies
        ! that one does \sum_R R_i \int bra(r-R) ket(r) dr, where "i"
        ! is the factor type. In normal operation, there is no R_i ! gcc32

    ! Local Variables
    integer :: col, recv_row
    integer :: local_col, local_row
    integer :: first_row_on_proc
    integer :: recv_proc
    logical :: remote_row, remote_req_row
    type(COEF) :: bracket_el             ! pdh: matrix element
    integer :: ierr                      ! ndmh: error flag
    integer :: idx_len
    integer,allocatable :: bracket_idx(:)
    ! ndmh: variables for planned execution
    integer :: req_row
    integer :: plan_step                 ! counter for current step in plan
    integer :: prev_req_row              ! Last row function requested
    integer :: prev_recv_row             ! Last row function received
    integer :: req_proc
    integer :: last_step
    integer :: loop_step
    type(FUNC_COMM) :: row_comm
!$  integer, external :: omp_get_thread_num

    call timer_clock('function_ops_brappd_ketppd',1)

#ifdef ITC_TRACE
    call VTBEGIN(vt_integrals_brappd_ketppd,vt_err)
#endif

    ! cks: synchronize PEs
    call comms_barrier

    ! ndmh: allocate workspace
    row_comm%plan_steps = sparse_proc_num_element(bracket)
    call function_ops_func_comm_init(row_comm,bra_basis,2,bras_on_grid%iscmplx)

    ! ndmh: get the matrix index
    idx_len = sparse_index_length(bracket)
    allocate(bracket_idx(idx_len),stat=ierr)
    call utils_alloc_check('function_ops_brappd_ketppd','bracket_idx',ierr)
    call sparse_generate_index(bracket_idx,bracket)

    ! ndmh: create a plan from the matrix index
    local_col = sparse_num_elems_on_proc(pub_my_proc_id,bracket,'C')
    call function_ops_batch_row_plan(row_comm,idx_len,bracket_idx, &
         bracket,1,local_col,'FULL')

    deallocate(bracket_idx,stat=ierr)
    call utils_dealloc_check('function_ops_brappd_ketppd','bracket_idx',ierr)


    first_row_on_proc = sparse_first_elem_on_proc(pub_my_proc_id,bracket,'R')

    ! ndmh: initialisation of send request receive operations
    call function_ops_init_requests(row_comm)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(plan_step,last_step,loop_step,col,local_col,bracket_el, &
!$OMP      local_row,remote_row,recv_row,prev_recv_row,recv_proc, &
!$OMP      req_row,req_proc,remote_req_row,prev_req_row) &
!$OMP SHARED(row_comm,pub_my_proc_id,bracket,first_row_on_proc, &
!$OMP      bras_on_grid,kets_on_grid,bra_basis,ket_basis,cell, &
!$OMP      pub_threads_max, factor_type)

    ! ndmh: initializations
    prev_recv_row = -2
    prev_req_row = -2
    local_row = -1
    remote_row = .false.
    remote_req_row = .false.
    plan_step = -1
    recv_row = -1

    ! ndmh: loop over the steps of the plan
    do
       plan_step = plan_step + 1
       if (plan_step>row_comm%plan_steps) exit

       ! ndmh: skip recv if still pre-start
       if (plan_step > 0) then

          ! ndmh: find the row and col of the matrix element to calculate
          recv_row = row_comm%plan(1,plan_step)

          ! ndmh: receive or copy the ppd list and offset into the buffer sphere
          if (prev_recv_row /= recv_row) then

             ! ndmh: find the proc the row function is stored on
             recv_proc = sparse_proc_of_elem(recv_row,bracket,'R')

             ! ndmh: finish the recv operation for upcoming step
             if (recv_proc==pub_my_proc_id) then
                local_row = recv_row - first_row_on_proc + 1
                remote_row = .false.
             else
!$OMP MASTER
                call function_ops_finish_recv(recv_row,row_comm,bra_basis, &
                     bras_on_grid)
!$OMP END MASTER
!$OMP FLUSH(row_comm)
                remote_row = .true.
             end if
             prev_recv_row = recv_row

          end if

       end if

       ! ndmh: find last step in plan which uses current recv_row
       last_step = plan_step
       do
          if (last_step==row_comm%plan_steps) exit
          if (row_comm%plan(1,last_step+1)==recv_row) then
             last_step=last_step + 1
          else
             exit
          end if
       end do

       ! ndmh: check if we need to request a function to be sent to this proc
       ! ndmh: for the next step
       if (last_step+1<=row_comm%plan_steps) then
          req_row = row_comm%plan(1,last_step+1)
          req_proc = sparse_proc_of_elem(req_row,bracket,'R')
          ! ndmh: if this function is not local to this proc and we do not already
          ! ndmh: have it, then send a request for it to req_proc
          if (req_row /= prev_req_row) then
             if (req_proc /= pub_my_proc_id) then
                remote_req_row = .true.
!$OMP MASTER
                call function_ops_start_recv(req_proc,req_row,row_comm, &
                     bra_basis)
                call function_ops_request(req_proc,req_row,row_comm)
!$OMP END MASTER
             else
                remote_req_row = .true.
             end if
             prev_req_row = req_row
          end if
       end if

       ! ndmh: respond to any send requests made of this proc
!$OMP MASTER
       call function_ops_respond_to_reqs(row_comm,bra_basis,bras_on_grid)
!$OMP END MASTER

       ! ndmh: cycle if still pre-start
       if (plan_step < 1) cycle

       if (remote_row) then
!$OMP BARRIER
       end if

!$OMP DO
       do loop_step=plan_step,last_step

          ! Find local index of col on this proc
          col = row_comm%plan(2,loop_step)
          local_col = col - &
               sparse_first_elem_on_proc(pub_my_proc_id,bracket,'C') + 1

          ! ndmh: calculate bracket matrix element
          if (remote_row) then
             bracket_el = &
                  internal_bra_dot_ket_ppds(row_comm%recv_func_on_grid(row_comm%ubuf), &
                  kets_on_grid, row_comm%recv_spheres(row_comm%ubuf), &
                  ket_basis%spheres(local_col), factor_type)
          else
             bracket_el = &
                  internal_bra_dot_ket_ppds(bras_on_grid, &
                  kets_on_grid, bra_basis%spheres(local_row), &
                  ket_basis%spheres(local_col), factor_type)
          end if
          ! bracket_el should be complex if bracket is complex
          call utils_assert(bracket_el%iscmplx .eqv. bracket%iscmplx, &
             'Error in function_ops_brappd_ketppd: &
             &mismatched real/complex brackets.')

          ! cks: scale with grid point weight

          call data_coef_scale(bracket_el, cell%weight)

          if (bracket_el%iscmplx) then
              ! ndmh: deposit in matrix
              call sparse_put_element(bracket_el%z,bracket,recv_row,col)
          else
              ! ndmh: deposit in matrix
              call sparse_put_element(bracket_el%d,bracket,recv_row,col)
          end if

!$        if (omp_get_thread_num()==0) then
          ! ndmh: respond to any send requests made of this proc
          call function_ops_respond_to_reqs(row_comm,bra_basis,bras_on_grid)
!$        end if

       end do
!$OMP END DO NOWAIT
       plan_step = last_step

       if (remote_req_row) then
!$OMP BARRIER
       end if

    end do
!$OMP END PARALLEL

    ! ndmh: wait until completion messages have been received from all
    ! ndmh: other procs with which communication is required.
    call function_ops_finalise(row_comm,bra_basis,bras_on_grid)

    ! ndmh: synchronise again so that non-blocking sends are completed and all
    ! ndmh: handles restored before deallocation
    call comms_free
    call comms_barrier

    ! ndmh: deallocate comms workspace
    call function_ops_func_comm_exit(row_comm)

    call timer_clock('function_ops_brappd_ketppd',2)

#ifdef ITC_TRACE
    call VTEND(vt_integrals_brappd_ketppd,vt_err)
#endif

contains

    type(COEF) function internal_bra_dot_ket_ppds(bra_on_grid, ket_on_grid, &
         bra_sphere, ket_sphere, factor_type)

      !==================================================================!
      ! This function calculates a <fa|fb> matrix element as the dot     !
      ! between the points in the ppds in common between the fa function !
      ! (the "bra") and the fb function (the "ket").                     !
      !------------------------------------------------------------------!
      ! Arguments:                                                       !
      ! bra_on_grid (input): The bra functions in ppd representation     !
      ! ket_on_grid (input): The ket functions in ppd representation     !
      ! n_bra_ppds (input) : Number of ppds of bra_on_grid               !
      ! n_ket_ppds (input) : Number of ppds of ket_on_grid               !
      ! bra_sphere (input) : Sphere for bra function                     !
      ! ket_sphere (input) : Sphere for ket function                     !
      !------------------------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 5/12/2003                    !
      ! Updated by Nicholas Hine in December 2009 to take advantage of   !
      ! the fact that ppds are listed in ascending order, so there is no !
      ! need for a full double loop over all ppds of the bras and kets.  !
      ! OpenMP parallelised by Nicholas Hine in May 2013.                !
      !==================================================================!

      use basis, only: SPHERE
      use constants, only: stdout
      use datatypes, only: COEF, FUNCTIONS, data_functions_dot
      use utils, only: utils_abort
      !use rundat, only: pub_threads_max

      implicit none

      ! Arguments
      !real(kind=DP), intent(in) :: bra_on_grid(:)
      ! agrecocmplx: use new FUNCTIONS type
      type(FUNCTIONS), intent(in) :: bra_on_grid
      !real(kind=DP), intent(in) :: ket_on_grid(:)
      type(FUNCTIONS), intent(in) :: ket_on_grid
      type(SPHERE), intent(in)  :: bra_sphere
      type(SPHERE), intent(in)  :: ket_sphere
      integer, intent(in), optional :: factor_type

      ! Local Variables
      integer :: ibra, iket              ! index in bra and ket ppd lists
      integer :: bra_ppd, ket_ppd        ! ppd index of bra and ket ppds
      integer :: bra_start, ket_start    ! index of first point of ppd
      !integer :: i                       ! point loop counter
      !real(kind=DP) :: ppd_sum
      real(kind=DP) :: a1_neighbour, a2_neighbour, a3_neighbour, factor ! gcc32
      type(COEF) :: temp_sum
      real(KIND=DP) :: ppd_sum_real
      real(KIND=DP) :: ppd_sum_imag

         ppd_sum_real = 0.0_DP
         ppd_sum_imag = 0.0_DP

      if (present(factor_type)) then
         if ((factor_type .gt. 3).or.(factor_type .lt. 1)) then
            call utils_abort('Wrong factor for brappd_ketppd')
         end if
      end if


!!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!!$OMP PRIVATE (ibra,bra_ppd,bra_start,ket_ppd,iket,ket_start,i,temp_sum) &
!!$OMP SHARED(bra_sphere,ket_sphere,cell,bra_on_grid,ket_on_grid, &
!!$OMP      pub_threads_max) &
!!$OMP REDUCTION(+:ppd_sum_real,ppd_sum_imag)

      ! Start at the beginning of the ket ppd list
      iket = 1
      ket_ppd = ket_sphere%ppd_list(1,iket)

!!$OMP DO
      do ibra=1,bra_sphere%n_ppds_sphere
         ! ndmh: find ppd number and start position for ibra
         bra_ppd = bra_sphere%ppd_list(1,ibra)
         bra_start = bra_sphere%offset + (ibra-1)*cell%n_pts
         do
            ! ndmh: keep moving on while ket_ppd is less than bra_ppd
            if (ket_ppd < bra_ppd) then
               iket = iket + 1
               if (iket > ket_sphere%n_ppds_sphere) exit
               ket_ppd = ket_sphere%ppd_list(1,iket)
               cycle
            end if

            ! ndmh: ppd numbers match, so add product of these ppds
            if (ket_ppd == bra_ppd) then
               ket_start = ket_sphere%offset + (iket-1)*cell%n_pts
               !do i=0,cell%n_pts-1
               !   ppd_sum = ppd_sum + bra_on_grid(bra_start+i) * &
               !        ket_on_grid(ket_start+i)
               !end do

               ! agrecocmplx: use new function
               temp_sum = data_functions_dot(bra_on_grid, ket_on_grid, &
                          bra_start, ket_start, cell%n_pts)

               if (present(factor_type)) then
                  a1_neighbour = real( nint(real(ket_sphere%ppd_list(2,iket), &
                       kind=DP) / 9.0_DP), kind=DP)
                  a2_neighbour = real( nint(real((ket_sphere%ppd_list(2,iket) -&
                       9 * int(a1_neighbour)), kind=DP) / 3.0_DP), kind=DP)
                  a3_neighbour = real( (ket_sphere%ppd_list(2,iket) - &
                       9 * int(a1_neighbour) - 3 * int(a2_neighbour)), kind=DP)

                  ! gcc32: we used the KET sphere because we calculate
                  ! \sum_R R \int BRA(r-R) KET(r) dr, i.e. we are looping over
                  ! the ppds of BRA duplicates from other cells that
                  ! overlap the sphere of the KET. One can see from which
                  ! neighbouring cell these BRA ppds come from by using
                  ! ket_sphere%ppd_list(2,iket), since we are only interested
                  ! in the ppds that appear in the KET sphere.


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

               ! agrecocmplx: distinguish between real-real and
               ! complex-complex case
               if (temp_sum%iscmplx) then
                   ppd_sum_real = ppd_sum_real + factor * real(temp_sum%z)
                   ppd_sum_imag = ppd_sum_imag + factor * aimag(temp_sum%z)
               else
                   ppd_sum_real = ppd_sum_real + factor * temp_sum%d
                   ! ppd_sum_imag = ppd_sum_imag + 0.0_DP
               end if

            end if

            ! ndmh: move on to next bra ppd as soon as we have found a match
            ! ndmh: or moved beyond bra_ppd
            if (ket_ppd >= bra_ppd) exit

         end do
      end do
!!$OMP END DO

!!$OMP END PARALLEL

      if (bra_on_grid%iscmplx.or.ket_on_grid%iscmplx) then
          internal_bra_dot_ket_ppds%iscmplx = .true.
          internal_bra_dot_ket_ppds%z = cmplx(ppd_sum_real,ppd_sum_imag,kind=DP)
      else
          internal_bra_dot_ket_ppds%iscmplx = .false.
          internal_bra_dot_ket_ppds%d = ppd_sum_real
      end if

    end function internal_bra_dot_ket_ppds

  end subroutine function_ops_brappd_ketppd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_ops_batch_row_plan(fcomm,idx_len, &
       sparse_idx,mat,local_start,local_end,scheme)

    !=========================================================================!
    ! This subroutine creates a list of row, column pairs of funcs that need  !
    ! to be calculated on the local proc. This can then be used to execute    !
    ! row sums such as density_batch_row_sums and ngwf_gradient_fb_sums in a  !
    ! maximally efficient manner                                              !
    !-------------------------------------------------------------------------!
    !  Arguments:                                                             !
    !    num_plan_steps (inout) : on input: number of steps in my_plan array  !
    !                             on output: number of steps actually needed  !
    !    plan(1,:,:)    (out) : rows of row/col pairs this proc will          !
    !                           need to calculate                             !
    !    plan(2,:,:)    (out) : cols of row/col pairs this proc will          !
    !                           need to calculate                             !
    !    idx_len        (in)  : length of sparse index sparse_idx             !
    !    sparse_idx     (in)  : index of the matrix being calculated          !
    !    local_start    (in)  : start of col functions in this batch          !
    !    local_end      (in)  : end of col functions in this batch            !
    !    scheme         (in)  : string identifying which scheme to use        !
    !                           'FULL' = calculate all elements               !
    !                           'LOWER' = calculate lower triangle only       !
    !                           'ALTERNATE' = calculate alternating elements  !
    !                           from upper and lower triangles                !
    !    reqs           (out) : array of logical flags to store whether this  !
    !                           proc will need to request anything from each  !
    !                           of the other procs                            !
    !-------------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine in May 2008 based on code  !
    ! from the old integrals_brappd_ketfftbox by Chris-Kriton Skylaris        !
    !=========================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs
    use constants, only: stdout
    use parallel_strategy, only: PARAL_INFO
    use sparse, only: SPAM3, sparse_first_elem_on_proc, sparse_atom_of_elem, &
         pattern_full, pattern_alternate, pattern_lower, sparse_get_par, &
         sparse_last_elem_on_proc
    use utils, only: utils_abort

    ! Arguments
    type(FUNC_COMM), intent(inout) :: fcomm
    integer, intent(in) :: idx_len
    integer, intent(in) :: sparse_idx(idx_len)
    type(SPAM3),intent(in) :: mat
    integer, intent(in) :: local_start, local_end
    character(*), intent(in) :: scheme

    ! Locals
    integer :: idx                   ! sparse matrix index counter
    integer :: local_col_atom        ! atom of local_col
    integer :: loc_local_col_atom    ! local index of above atom
    integer :: pattern               ! pattern determining elements calculated
    integer :: step                  ! proc-block loop counter
    integer :: plan_step             ! position counter for plan array
    integer :: recv_proc             ! proc to receive data from
    integer :: recv_row              ! the function being received
    integer :: recv_row_atom         ! the atom of the function being received
    integer :: col                   ! global index of column function
    integer :: local_col             ! local index of column function
    integer :: first_col             ! first column on this proc
    integer :: first_row             ! first row on recvproc
    integer :: last_row              ! last row on recvproc
    type(PARAL_INFO), pointer :: col_par ! rc2013: column parallel strategy

    ! ndmh: initialise plan variables
    fcomm%plan = 0
    fcomm%reqs_in = .false.
    plan_step = 1

    ! rc2013: obtain parallel strategy for this SPAM3
    call sparse_get_par(col_par, mat, 'C')

    ! ndmh: Set calculation pattern
    select case (scheme)
    case ('FULL','full','ASYM','asym')
       pattern = pattern_full
    case ('LOWER','lower')
       pattern = pattern_lower
    case ('ALT','alt','ALTERNATE','alternate')
       pattern = pattern_alternate
    case default
       call utils_abort('Error in function_ops_batch_row_plan: &
            &calculation pattern "'//trim(scheme)//'" not recognised')
    end select

    ! ndmh: Loop over segments of the matrix, starting from
    ! ndmh: the diagonal segments
    first_col = sparse_first_elem_on_proc(pub_my_proc_id,mat,'C')
    do step=0,pub_total_num_procs-1

       ! ndmh: find proc to check overlaps with on this step
       recv_proc = modulo(pub_my_proc_id + step, pub_total_num_procs)

       ! Loop over the rows on recv_proc
       first_row = sparse_first_elem_on_proc(recv_proc,mat,'R')
       last_row = sparse_last_elem_on_proc(recv_proc,mat,'R')
       do recv_row=first_row,last_row
          recv_row_atom = sparse_atom_of_elem(recv_row,mat,'R')

          ! ndmh: loop over cols on this proc to see whether they
          ! ndmh: overlap recv_row
          do local_col=local_start,local_end
             col = local_col + first_col - 1
             local_col_atom = sparse_atom_of_elem(col,mat,'C')
             loc_local_col_atom = local_col_atom - &
                  col_par%first_atom_on_proc(pub_my_proc_id) + 1
             do idx=sparse_idx(loc_local_col_atom), &
                  sparse_idx(loc_local_col_atom+1)-1

                ! ndmh: add to plan only if sparse element exists...
                if (sparse_idx(idx) == recv_row_atom) then

                   ! ndmh: ... and is in the pattern being used
                   if (internal_in_pattern(recv_row,col)) then
                      fcomm%plan(1,plan_step) = recv_row
                      fcomm%plan(2,plan_step) = col
                      fcomm%reqs_in(recv_proc) = .true.
                      plan_step = plan_step + 1
                   end if
                end if

             end do
          end do
       end do

    end do

    ! ndmh: return index of last entry in plan on this proc
    fcomm%plan_steps = plan_step - 1

  contains

    !=========================================================================!
    ! Checks if this proc needs to receive and use the data for this col, row !
    ! combination from another proc, based on the variable pattern.           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! row (input) : row function to consider receiving in full matrix         !
    ! col (input) : column function to consider receiving in full matrix      !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, April 2008                                    !
    !=========================================================================!

    logical function internal_in_pattern(row,col)

      use utils, only: utils_abort

      integer, intent(in) :: col
      integer, intent(in) :: row

      select case (pattern)
      case (pattern_full)
         internal_in_pattern = .true.
      case (pattern_lower)
         if (col<=row) then
            internal_in_pattern = .true.
         else
            internal_in_pattern = .false.
         end if
      case (pattern_alternate)
         if (((col<row).and.(mod(col+row,2)==1)).or. &
              ((col>row).and.(mod(col+row,2)==0)).or. &
              (col==row)) then
            internal_in_pattern = .true.
         else
            internal_in_pattern = .false.
         end if
      case default
         internal_in_pattern = .false.
         call utils_abort('Invalid pattern supplied to &
                 &internal_in_pattern (function_ops_batch_row_plan)')
      end select

    end function internal_in_pattern

  end subroutine function_ops_batch_row_plan


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_ops_batch_col_start(col_box_start,col_start_in_box, &
       batch_size,local_start,local_end,fftbox,cell,col_basis)

    !=========================================================================!
    ! This subroutine creates a list of column function / ket start positions !
    ! for a batch of functions defined by local_start and local_end. It       !
    ! chooses where they start in an FFTbox
    !-------------------------------------------------------------------------!
    !  Arguments:                                                             !
    !    col_box_start    (out) : Start positions of fftboxes for each func   !
    !    col_start_in_box (out) : Start positions of functions within fftboxes!
    !    batch_size        (in) : Number of functions in the batch            !
    !    local_start       (in) : Local index of first function in batch      !
    !    local_end         (in) : Local index of last function in batch       !
    !    fftbox            (in) : FFTBOX_INFO description of fftbox           !
    !    col_basis         (in) : FUNC_BASIS description of the functions     !
    !-------------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine in November 2014           !
    ! incorporating bits of earlier code by Nicholas Hine and Chris-Kriton    !
    ! Skylaris from 2004-2009.                                                !
    !=========================================================================!

    use basis, only: basis_location_func_wrt_cell, basis_ket_start_wrt_fftbox, &
         basis_point_wrt_box, basis_start_of_box_wrt_cell
    use fft_box, only: FFTBOX_INFO
    use simulation_cell, only: CELL_INFO

    ! Arguments
    integer, intent(in) :: batch_size
    integer, intent(in) :: local_start
    integer, intent(in) :: local_end
    integer, intent(out) :: col_box_start(3,batch_size)
    integer, intent(out) :: col_start_in_box(3,batch_size)
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(CELL_INFO), intent(in) :: cell
    type(FUNC_BASIS), intent(in) :: col_basis

    ! Local Variables
    integer :: local_col
    integer :: batch_index
    integer :: col_cell_start(3)
    ! {C}entre wrt fft {B}ox in terms of number of {G}ridpoints
    ! in each lattice vector direction. Not necessarily integer values.
    real(kind=DP) :: cbg1,cbg2,cbg3

    integer, parameter :: nonsense = -123456789

    ! -------------------------------------------------------------------------
    col_start_in_box(:,:) = nonsense !jd: Take care to initialise entire out arg

    ! ndmh: find fftbox start positions for each col/ket function in this batch
    batch_index = 1
    do local_col=local_start,local_end

#ifdef FFTBOX_OLD_POS
       ! ndmh: determine start position of ket within FFTbox
       call basis_ket_start_wrt_fftbox(col_start_in_box(1,batch_index), &
            col_start_in_box(2,batch_index),col_start_in_box(3,batch_index), &
            fftbox%total_pt1, fftbox%total_pt2, fftbox%total_pt3, fftbox)

       ! ndmh: determine position of col function wrt start of simulation cell
       call basis_location_func_wrt_cell(col_cell_start(1), &
            col_cell_start(2), col_cell_start(3), &
            col_basis%tight_boxes(local_col), cell)

       ! ndmh: determine start of FFTbox wrt start of simulation cell
       col_box_start(1,batch_index) = col_cell_start(1) &
            - col_start_in_box(1,batch_index) + 1
       col_box_start(2,batch_index) = col_cell_start(2) &
            - col_start_in_box(2,batch_index) + 1
       col_box_start(3,batch_index) = col_cell_start(3) &
            - col_start_in_box(3,batch_index) + 1
#else
       ! Centre of projector wrt fftbox in terms of grid spacings
       call basis_point_wrt_box(cbg1, cbg2, cbg3, &  ! output
            col_basis%spheres(local_col)%centre, fftbox%total_pt1, &
            fftbox%total_pt2, fftbox%total_pt3, cell, fftbox)

       ! Start of fftbox wrt cell in terms of grid-point number
       call basis_start_of_box_wrt_cell(col_box_start(1,batch_index), &
            col_box_start(2,batch_index),col_box_start(3,batch_index), &
            col_basis%spheres(local_col)%centre,cbg1,cbg2,cbg3,cell)

       ! ndmh: determine position of col function wrt start of simulation cell
       call basis_location_func_wrt_cell(col_cell_start(1), &
            col_cell_start(2), col_cell_start(3), &
            col_basis%tight_boxes(local_col), cell)

       ! ndmh: determine start of col function wrt start of box
       col_start_in_box(1,batch_index) = col_cell_start(1) &
            - col_box_start(1,batch_index) + 1
       col_start_in_box(2,batch_index) = col_cell_start(2) &
            - col_box_start(2,batch_index) + 1
       col_start_in_box(3,batch_index) = col_cell_start(3) &
            - col_box_start(3,batch_index) + 1
#endif

       batch_index = batch_index + 1

    end do

  end subroutine function_ops_batch_col_start


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_ops_func_comm_init(fcomm,fbasis,nbuf,is_cmplx)

    !=========================================================================!
    ! Allocates storage for a FUNC_COMM object, which manages communication   !
    ! of functions between parallel procs.                                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! fcomm (inout) : FUNC_COMM type to allocate                              !
    ! fbasis (input) : corresponding FUNC_BASIS                               !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 15/06/2012                                  !
    !=========================================================================!

    use datatypes, only: data_functions_alloc
    use comms, only: pub_total_num_procs
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(FUNC_COMM), intent(inout) :: fcomm
    type(FUNC_BASIS), intent(in) :: fbasis
    integer, intent(in) :: nbuf
    logical, intent(in), optional :: is_cmplx

    ! Local Variables
    integer :: ibuf
    integer :: ierr
    logical :: loc_cmplx

    if(present(is_cmplx)) then
       loc_cmplx = is_cmplx
    else
       loc_cmplx = .false.
    end if

    fcomm%nbuf = nbuf

    ! ndmh: allocate comms workspace
    allocate(fcomm%plan(2,fcomm%plan_steps),stat=ierr)
    call utils_alloc_check('function_ops_func_comm_init','fcomm%plan',ierr)
    allocate(fcomm%func_requests(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('function_ops_func_comm_init', &
         'fcomm%func_requests',ierr)
    allocate(fcomm%handles(0:pub_total_num_procs+2*fcomm%nbuf-1),stat=ierr)
    call utils_alloc_check('function_ops_func_comm_init', &
         'fcomm%handles',ierr)
    allocate(fcomm%reqs_in(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('function_ops_func_comm_init', &
         'fcomm%reqs_in',ierr)
    allocate(fcomm%reqs_out(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('function_ops_func_comm_init', &
         'fcomm%reqs_out',ierr)

    ! ndmh: set up receiving spheres and buffers
    allocate(fcomm%req_buffer(1:fcomm%nbuf),stat=ierr)
    call utils_alloc_check('function_ops_func_comm_init', &
         'fcomm%req_buffer',ierr)
    !allocate(fcomm%recv_func_on_grid(fbasis%func_on_grid_buffer_size, &
    !     1:fcomm%nbuf),stat=ierr)
    !call utils_alloc_check('function_ops_func_comm_init', &
    !     'fcomm%recv_func_on_grid',ierr)
    allocate(fcomm%recv_func_on_grid(1:fcomm%nbuf),stat=ierr)
    call utils_alloc_check('function_ops_func_comm_init', &
         'fcomm%recv_func_on_grid',ierr)
    allocate(fcomm%recv_spheres(1:fcomm%nbuf),stat=ierr)
    call utils_alloc_check('function_ops_func_comm_init', &
         'fcomm%recv_spheres',ierr)
    do ibuf=1,fcomm%nbuf
       fcomm%recv_spheres(ibuf)%n_ppds_sphere = fbasis%func_on_grid_buffer_size
       fcomm%recv_spheres(ibuf)%offset = 1
       fcomm%recv_spheres(ibuf)%radius = 0.0_DP
       allocate(fcomm%recv_spheres(ibuf)%ppd_list(2, &
            fcomm%recv_spheres(ibuf)%n_ppds_sphere),stat=ierr)
       call utils_alloc_check('function_ops_func_comm_init',&
            'fcomm%recv_spheres(ibuf)%ppd_list',ierr)
       ! agrecocmplx: allocate FUNCTIONS type with the appropriate routine
       call data_functions_alloc(fcomm%recv_func_on_grid(ibuf), &
            fbasis%func_on_grid_buffer_size, loc_cmplx)
    end do

    fcomm%rbuf = 1
    fcomm%ubuf = 0

  end subroutine function_ops_func_comm_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_ops_func_comm_exit(fcomm)

    !=========================================================================!
    ! Deallocates storage for a FUNC_COMM object, which manages communication !
    ! of functions between parallel procs.                                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! fcomm (inout) : FUNC_COMM type to deallocate                            !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 15/06/2012                                  !
    !=========================================================================!

    use datatypes, only: data_functions_dealloc
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_COMM), intent(inout) :: fcomm

    ! Local Variables
    integer :: ibuf
    integer :: ierr

    ! ndmh: deallocate receiving spheres and buffers
    do ibuf=fcomm%nbuf,1,-1
       deallocate(fcomm%recv_spheres(ibuf)%ppd_list,stat=ierr)
       call utils_dealloc_check('function_ops_func_comm_exit',&
            'fcomm%recv_spheres(ibuf)%ppd_list',ierr)
       ! agrecocmplx: deallocate FUNCTIONS type with
       ! appropriate routine
       call data_functions_dealloc(fcomm%recv_func_on_grid(ibuf))
    end do
    deallocate(fcomm%recv_spheres,stat=ierr)
    call utils_dealloc_check('function_ops_func_comm_exit', &
         'fcomm%recv_spheres',ierr)
    deallocate(fcomm%recv_func_on_grid,stat=ierr)
    call utils_dealloc_check('function_ops_func_comm_exit', &
         'fcomm%recv_func_on_grid',ierr)
    deallocate(fcomm%req_buffer,stat=ierr)
    call utils_dealloc_check('function_ops_func_comm_exit', &
         'fcomm%req_buffer',ierr)

    ! ndmh: deallocate comms workspace
    deallocate(fcomm%reqs_out,stat=ierr)
    call utils_dealloc_check('function_ops_func_comm_exit', &
         'fcomm%reqs_out',ierr)
    deallocate(fcomm%reqs_in,stat=ierr)
    call utils_dealloc_check('function_ops_func_comm_exit', &
         'fcomm%reqs_in',ierr)
    deallocate(fcomm%handles,stat=ierr)
    call utils_dealloc_check('function_ops_func_comm_exit', &
         'fcomm%handles',ierr)
    deallocate(fcomm%func_requests,stat=ierr)
    call utils_dealloc_check('function_ops_func_comm_exit', &
         'fcomm%func_requests',ierr)
    deallocate(fcomm%plan,stat=ierr)
    call utils_dealloc_check('function_ops_func_comm_exit', &
         'fcomm%plan',ierr)

  end subroutine function_ops_func_comm_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_ops_init_requests(fcomm)

    use comms, only: comms_alltoall, comms_irecv, pub_my_proc_id, &
         pub_null_handle, pub_total_num_procs

    implicit none

    ! Arguments
    type(FUNC_COMM), intent(inout) :: fcomm

    ! Locals
    integer :: proc

    ! ndmh: initialisation
    fcomm%probe_count = 0
    fcomm%handles = pub_null_handle
    fcomm%reqs_in(pub_my_proc_id) = .false.
    fcomm%reqs_out = fcomm%reqs_in
    call comms_alltoall(fcomm%reqs_out,1)

    ! ndmh: loop over all procs
    do proc=0,pub_total_num_procs-1
       if (fcomm%reqs_out(proc)) then
          ! ndmh: start asynchronous receive for request from this proc
          ! jd: Relying on this value persisting after the call to MPI_IRECV
          !     is flaky:
          !     https://bitbucket.org/onetep/onetep/issues/1890/illegal-peeking-into-mpi_irecvs-buffer
          !     ... but we decided to let it slide for the time being, cf.
          !     ONETEP development slack / repository, 2022.04.06.
          fcomm%func_requests(proc) = FUNCS_WAITING
          call comms_irecv(proc,fcomm%func_requests(proc),1,tag=req_tag+proc, &
               handle=fcomm%handles(proc))
       else
          fcomm%func_requests(proc) = FUNCS_DONE
          fcomm%handles(proc) = pub_null_handle
       end if
    end do

    fcomm%req_buffer(:) = 0

  end subroutine function_ops_init_requests


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_ops_request(req_proc,req_row,fcomm)

    use comms, only: comms_send, pub_my_proc_id

    implicit none

    ! Arguments
    type(FUNC_COMM), intent(inout) :: fcomm
    integer,intent(in) :: req_proc
    integer,intent(in) :: req_row

    ! Store requested row index in buffer to asynchronously send
    fcomm%req_buffer(fcomm%rbuf) = req_row

    ! ndmh: send request for function required
    call comms_send(req_proc,fcomm%req_buffer(fcomm%rbuf),1, &
         tag=req_tag+pub_my_proc_id)

    ! ndmh: advance buffer index and cycle if filled
    fcomm%rbuf = fcomm%rbuf + 1
    if (fcomm%rbuf > fcomm%nbuf) fcomm%rbuf = 1

  end subroutine function_ops_request


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_ops_finalise(fcomm,fbasis,funcs_on_grid)

    use datatypes, only: FUNCTIONS
    use comms, only: comms_waitany, comms_wait, pub_total_num_procs

    implicit none

    ! Arguments
    type(FUNC_COMM), intent(inout) :: fcomm
    type(FUNC_BASIS), intent(in) :: fbasis
    !real(kind=DP), intent(in)  :: funcs_on_grid(fbasis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: funcs_on_grid

    ! Locals
    integer :: proc

    ! ndmh: send signal indicating completion to all procs
    do proc=0,pub_total_num_procs-1
       if (fcomm%reqs_in(proc)) then
          call function_ops_request(proc,-2000,fcomm)
       end if
    end do

    do
       if (all(fcomm%func_requests(:)==FUNCS_DONE)) exit

       ! ndmh: respond to any outstanding requests
       call function_ops_respond_to_reqs(fcomm,fbasis,funcs_on_grid)

       ! ndmh: if all other procs have not yet finished, wait for a new request
       if (any(fcomm%func_requests(:)==FUNCS_WAITING)) then
          call comms_waitany(pub_total_num_procs,fcomm%handles)
       end if
    end do

    ! ndmh: ensure all procs have completed their final receive operation
    do proc=0,pub_total_num_procs-1
       call comms_wait(fcomm%handles(proc))
    end do

  end subroutine function_ops_finalise


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_ops_respond_to_reqs(fcomm,fbasis,funcs_on_grid)

    use basis, only: SPHERE
    use comms, only: comms_irecv, comms_wait, comms_probe, pub_my_proc_id, &
         pub_total_num_procs
    use datatypes, only: FUNCTIONS

    implicit none

    ! Arguments
    type(FUNC_COMM), intent(inout) :: fcomm
    type(FUNC_BASIS), intent(in) :: fbasis
    !real(kind=DP), intent(in)  :: funcs_on_grid(fbasis%size_on_grid)
    ! agreco: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: funcs_on_grid

    ! Locals
    integer :: proc
    logical :: result
    logical :: can_exit

    ! Every probe_frequency times this routine is called, we must enter
    ! the MPI library to ensure any ongoing sends complete even if no
    ! requests are received.
    fcomm%probe_count = fcomm%probe_count + 1
    if (fcomm%probe_count >= probe_frequency) then
       call comms_probe(result,pub_my_proc_id)
       fcomm%probe_count = 0
    end if

    ! ndmh: do not continue until there are no outstanding requests
    can_exit = .true.
    do
       ! ndmh: loop over procs checking for requests
       do proc=0,pub_total_num_procs-1

          if (fcomm%func_requests(proc)<0) cycle
          can_exit = .false.

          ! ndmh: ensure receive has completed
          call comms_wait(fcomm%handles(proc))

          ! ndmh: send requested ngwf to proc
          call function_ops_send(proc,fcomm%func_requests(proc),fbasis, &
               funcs_on_grid)

          ! ndmh: re-initialise the request receive operation
          fcomm%func_requests(proc) = FUNCS_WAITING
          call comms_irecv(proc,fcomm%func_requests(proc),1,tag=req_tag+proc, &
               handle=fcomm%handles(proc))
          call comms_probe(result,proc)

       end do
       if (can_exit) then
          exit
       else
          can_exit = .true.
       end if
    end do

  end subroutine function_ops_respond_to_reqs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_ops_send(send_proc,send_row,fbasis,funcs_on_grid)

    use basis, only: SPHERE
    use comms, only: comms_send, pub_my_proc_id
    use datatypes, only: FUNCTIONS

    implicit none

    ! Arguments
    integer,intent(in) :: send_proc
    integer,intent(in) :: send_row
    type(FUNC_BASIS), intent(in) :: fbasis
    !real(kind=DP), intent(in)  :: funcs_on_grid(fbasis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: funcs_on_grid

    ! Locals
    integer :: local_row
    integer :: n_ppds,n_pts,start_pt

    local_row = send_row - fbasis%first_on_proc(pub_my_proc_id) + 1
    n_ppds = fbasis%n_ppds_sphere(send_row)
    n_pts = n_ppds*fbasis%n_pts
    start_pt = fbasis%spheres(local_row)%offset

    ! ndmh: send sphere information and ppd data for this function
    call comms_send(send_proc,fbasis%spheres(local_row)%ppd_list(:,:), &
         2*n_ppds,tag=send_row*2+0)
    ! agrecocmplx: send real or complex data if dealing with
    ! real or complex functions
    if (funcs_on_grid%iscmplx) then
       call comms_send(send_proc,funcs_on_grid%z(start_pt:start_pt+n_pts-1), &
            n_pts,tag=send_row*2+1)
    else
       call comms_send(send_proc,funcs_on_grid%d(start_pt:start_pt+n_pts-1), &
            n_pts,tag=send_row*2+1)
    end if

  end subroutine function_ops_send


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_ops_start_recv(recv_proc,recv_row,fcomm, &
       fbasis)

    use basis, only: SPHERE
    use comms, only: comms_irecv, pub_total_num_procs

    implicit none

    ! Arguments
    integer,intent(in) :: recv_proc
    integer,intent(in) :: recv_row
    type(FUNC_COMM), intent(inout) :: fcomm
    type(FUNC_BASIS), intent(in) :: fbasis

    ! Locals
    integer :: n_ppds,n_pts
    integer :: ihandle

    ! ndmh: find sizes of buffers
    n_ppds = fbasis%n_ppds_sphere(recv_row)
    n_pts = n_ppds*fbasis%n_pts
    fcomm%recv_spheres(fcomm%rbuf)%n_ppds_sphere = n_ppds
    fcomm%recv_spheres(fcomm%rbuf)%offset = 1
    ihandle = pub_total_num_procs + fcomm%rbuf - 1

    ! ndmh: start asynchronous receive operations for the incoming ngwf
    call comms_irecv(recv_proc,fcomm%recv_spheres(fcomm%rbuf)%ppd_list(:,:), &
         2*n_ppds,tag=recv_row*2+0,handle=fcomm%handles(ihandle))
    ! agrecocmplx: receive real or complex data if using
    ! real or complex valued functions
    if (fcomm%recv_func_on_grid(fcomm%rbuf)%iscmplx) then
       call comms_irecv(recv_proc,fcomm%recv_func_on_grid(fcomm%rbuf)%z(:),n_pts, &
            tag=recv_row*2+1,handle=fcomm%handles(ihandle+1))
    else
       call comms_irecv(recv_proc,fcomm%recv_func_on_grid(fcomm%rbuf)%d(:),n_pts, &
            tag=recv_row*2+1,handle=fcomm%handles(ihandle+1))
    end if

  end subroutine function_ops_start_recv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_ops_finish_recv(recv_row,fcomm,fbasis,funcs_on_grid)

    use datatypes, only: FUNCTIONS
    use comms, only: comms_waitany, comms_wait, &
         pub_null_handle, pub_total_num_procs
    use utils, only: utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: recv_row
    type(FUNC_COMM), intent(inout) :: fcomm
    type(FUNC_BASIS), intent(in) :: fbasis
    !real(kind=DP), intent(in)  :: funcs_on_grid(fbasis%size_on_grid)
    ! agrecocmplx: use new type FUNCTIONS
    type(FUNCTIONS), intent(in) :: funcs_on_grid

    ! Local Variables
    integer :: ibuf, ihandle

    ! ndmh: find the buffer this function was being received to
    do ibuf=1,fcomm%nbuf
       if (fcomm%req_buffer(ibuf)==recv_row) exit
       call utils_assert(ibuf/=fcomm%nbuf, 'Error in function_ops_finish_recv&
               &: row to receive was not found in request buffer: ',recv_row)
    end do
    fcomm%ubuf = ibuf
    ihandle = pub_total_num_procs + ibuf - 1

    ! ndmh: now ensure this receive has completed
    do

       if (fcomm%handles(ihandle+1)==pub_null_handle) then
          ! call comms_wait to ensure receive has completed before exiting
          call comms_wait(fcomm%handles(ihandle))
          call comms_wait(fcomm%handles(ihandle+1))
          exit
       end if

       if (any(fcomm%func_requests(:)>0)) then
          ! ndmh: send any functions required by other procs, if there are
          ! ndmh: outstanding requests
          call function_ops_respond_to_reqs(fcomm,fbasis,funcs_on_grid)
       else
          ! ndmh: wait for either an incoming request, or one of the recv
          ! ndmh: operations to complete
          call comms_waitany(pub_total_num_procs+2*fcomm%nbuf,fcomm%handles)
       end if

    end do

  end subroutine function_ops_finish_recv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module function_ops
