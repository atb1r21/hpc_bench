! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-

module nonsc_forces

  !=====================================================================!
  ! This module contains the subroutines to calculate the correction to !
  ! the ONETEP forces due to non full self-consistency of the NGWF      !
  ! optimisation.                                                       !
  ! If the NGWFs are optimised to a tight threshold (as in a standard   !
  ! calculation), then the contribution of the non self-consistent force!
  ! should be residual. However, if the NGWFs are not totally optimised !
  ! this force should be taken into account.                            !
  !---------------------------------------------------------------------!
  ! *** Please note that this capability is still under development.    !
  !---------------------------------------------------------------------!
  ! This module was created by Alvaro Ruiz Serrano in November 2010.    !
  !=====================================================================!

  implicit none

  public :: nonsc_forces_ngwfs_calc
  public :: nonsc_forces_overlap_deriv


contains

  subroutine nonsc_forces_ngwfs_calc(nonsc_forces,&
       ngwfs_on_grid, contra_grad_on_grid, ngwf_basis, cell, fftbox, par, suffix)

    !======================================================================!
    ! This subroutine calculates the forces due to non self-consistency of !
    ! the outer loop. These forces behave like Pulay forces as if the NGWF !
    ! functions were the ultimate basis set.                               !
    ! See: J. Chem. Phys. vol 121, no 13, October 2004.                    !
    !======================================================================!
    !  Arguments:                                                          !
    !    nonsc_forces (out)  : NGWF non self-consistent forces in          !
    !                              cartesian coordinates array             !
    !    ngwfs_on_grid  (in) : current set of NGWFS on grid in PPD format  !
    !    contra_gradient (in): contravariant NGFW gradient on grid (PPD)   !
    !    ngwf_basis (in)     : the basis type for the NGWFs                !
    !======================================================================!
    ! Written by Alvaro Ruiz Serrano in January 2010.                      !
    ! Clean-up by Alvaro Ruiz Serrano in December 2010.                    !
    ! Modified to work with embedding by Joseph Prentice, May 2018         !
    !======================================================================!

    use datatypes, only: FUNCTIONS
    use bibliography, only: bibliography_cite
    use constants, only: DP, VERBOSE, CRLF
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_grad
    use parallel_strategy, only : PARAL_INFO
    use rundat, only: pub_output_detail
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_create, sparse_destroy
    use timer, only: timer_clock
    use utils, only: utils_flushed_string_output

    implicit none

    ! Arguments
    real(kind=DP),   intent(out) :: nonsc_forces(:,:)
    type(FUNC_BASIS),intent(in) :: ngwf_basis
    type(FFTBOX_INFO),intent(in) :: fftbox
    type(CELL_INFO), intent(in) :: cell
    type(FUNCTIONS),intent(in) :: ngwfs_on_grid
    type(FUNCTIONS),intent(in) :: contra_grad_on_grid
    type(PARAL_INFO),intent(in) :: par
    character(*),intent(in),optional :: suffix

    ! Local Variables
    ! ars: nabla matrices. The diagonal contains the forces
    type(SPAM3) :: nabla(3)
    integer :: ndim
    ! agrecocmplx
    logical :: loc_cmplx


    ! Start timer
    call timer_clock('nonsc_forces_ngwfs_calc',1)
    call bibliography_cite('PULAYFORCES')

    if(pub_output_detail >= VERBOSE) then
       call utils_flushed_string_output(&
            'Calculating NGWF non self-consistent forces...')
    end if

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid%iscmplx

    ! ars: create nablas
    do ndim = 1, 3
       ! jcap: add subregion suffix on if present
       if(present(suffix))then
          nabla(ndim)%structure = 'D'//trim(suffix)
       else
          nabla(ndim)%structure = 'D'
       end if
       ! agrecocmplx
       call sparse_create(nabla(ndim),iscmplx=loc_cmplx)
    end do

    ! ars: calculate nabla = <contra_grad|d(ngwf)/dx>
    call integrals_grad(nabla, contra_grad_on_grid, ngwf_basis, ngwfs_on_grid, &
         ngwf_basis, cell, fftbox)

    ! ars: at this point, the deriv of the NGWFs should be scaled according to
    !      d(ngwf)/dx = - d(ngwf)/dX_alpha
    !      and then, when extracting the forces,
    !      F = - sum_alpha <contra_grad|d(ngwf)/dX_alpha>
    ! ---> avoid multiply twice times -1

    ! ars: extract nonsc_forces forces
    call nonsc_forces_extract_forces(nonsc_forces, nabla, ngwf_basis, cell, par)

    ! ars: destroy nablas
    do ndim = 3, 1, -1
       call sparse_destroy(nabla(ndim))
    end do

    if(pub_output_detail >= VERBOSE) then
       call utils_flushed_string_output('done.'//CRLF)
    end if

    ! Stop timer
    call timer_clock('nonsc_forces_ngwfs_calc',2)


  end subroutine nonsc_forces_ngwfs_calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nonsc_forces_overlap_deriv(ov_deriv, ngwfs_on_grid, ngwf_basis, &
       cell,fftbox,par)

    !======================================================================!
    ! This subroutine calculates the derivative of the overlap matrix wrt  !
    ! the atomic positions in NGWF representation and stores the data in   !
    ! SPAM3 omega matrices corresponding to X, Y and Z coordinates.        !
    !======================================================================!
    !  Arguments:                                                          !
    !    ov_deriv (inout) : deriv of the overlap matrix wrt atpos elements !
    !    ngwfs_on_grid(in): current set of NGWFS on grid in PPD format     !
    !    ngwf_basis (in)  : The basis type for the NGWFs                   !
    !======================================================================!
    ! Written by Alvaro Ruiz Serrano in November 2009                      !
    ! Clean-up by Alvaro Ruiz Serrano in December 2010.                    !
    ! Modified to pass through par by Joseph Prentice, May 2018            !
    !======================================================================!

    use datatypes, only: FUNCTIONS
    use constants, only: DP
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_grad
    use parallel_strategy, only : PARAL_INFO
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_scale
    use timer, only: timer_clock


    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: ov_deriv(3)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FFTBOX_INFO),intent(in) :: fftbox
    type(CELL_INFO), intent(in) :: cell
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    type(PARAL_INFO), intent(in) :: par

    ! Local Variables
    integer :: ndim


    ! Start timer
    call timer_clock('nonsc_forces_overlap_deriv',1)

    ! ars: calculate ov_deriv = d/dx S_{alpha,beta}
    call integrals_grad(ov_deriv, ngwfs_on_grid, ngwf_basis, ngwfs_on_grid, &
         ngwf_basis, cell, fftbox)

    ! ars: scale according to
    !      d(ngwf)/dx = - d(ngwf)/dX_alpha
    do ndim = 1, 3
       call sparse_scale(ov_deriv(ndim), -1.0_DP)
    end do

    ! ars: add transpose to diagonal blocks
    call nonsc_forces_diagonal_blocks(ov_deriv, ngwf_basis, par)

    ! Stop timer
    call timer_clock('nonsc_forces_overlap_deriv',2)


  end subroutine nonsc_forces_overlap_deriv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nonsc_forces_extract_forces(nonsc_forces, nabla, ngwf_basis, cell, par)

    !======================================================================!
    ! This subroutine extracts the nonsc_forces forces from the diagonal   !
    ! terms of nablaXYZ matrices:                                          !
    ! F_atom = -sum(NGWFs_on atom) * nabla_atomblock(NGWF, NGWF)           !
    !======================================================================!
    !  Arguments:                                                          !
    !    nonsc_forces (out)  : nonsc_forces in cartesian coordinates array !
    !    nabla        (in)   : matrices whose diagonal elements contain the!
    !                          forces. The rest of the elements are junk.  !
    !======================================================================!
    ! Written by Alvaro Ruiz Serrano in January 2010.                      !
    ! Clean-up by Alvaro Ruiz Serrano in December 2010.                    !
    ! Modified to remove pub_par by Joseph Prentice, May 2018              !
    !======================================================================!

    use comms, only: pub_my_proc_id, comms_reduce
    use constants, only : DP
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: PARAL_INFO
    ! agrecocmplx
    use rundat, only: pub_imag_thr
    use simulation_cell, only: CELL_INFO
    use sparse, only: sparse_get_block, SPAM3
    use utils, only: utils_alloc_check, utils_dealloc_check, &
            utils_abort

    implicit none


    ! ars: arguments
    type(CELL_INFO), intent(in) :: cell
    real(kind=DP), intent(out) :: nonsc_forces(:,:)
    type(SPAM3), intent(in) :: nabla(3)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(PARAL_INFO), intent(in) :: par


    ! ars: local variables
    real(kind=DP), allocatable :: block_real(:,:)
    ! agrecocmplx
    complex(kind=DP), allocatable :: block_cmplx(:,:)

    integer :: iat, loc_iat, ingwf, ndim, ierr
    real(kind=DP) :: block_tr
    ! agrecocmplx
    logical :: loc_cmplx

    ! agrecocmplx
    loc_cmplx = nabla(1)%iscmplx

    ! ars: initialise
    nonsc_forces(:,:) = 0.0_DP

    ! agrecocmplx
    if (loc_cmplx) then
       allocate(block_cmplx(ngwf_basis%max_on_atom,ngwf_basis%max_on_atom),stat=ierr)
       call utils_alloc_check('nonsc_forces_extract_forces','block_cmplx',ierr)
    else
       allocate(block_real(ngwf_basis%max_on_atom,ngwf_basis%max_on_atom),stat=ierr)
       call utils_alloc_check('nonsc_forces_extract_forces','block_real',ierr)
    end if

    ! ars: extract forces on each atom
    iat = par%first_atom_on_proc(pub_my_proc_id)
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       do ndim = 1, 3
          block_tr = 0.0_DP
          ! agrecocmplx
          if (loc_cmplx) then
             block_cmplx(:,:) = (0.0_DP,0.0_DP)
             call sparse_get_block(block_cmplx,nabla(ndim),iat,iat)

             ! agrecocmplx: forces must be real
             do ingwf=1,ngwf_basis%num_on_atom(iat)
                ! agrecocmplx: abort if imaginary part of forces
                ! is above threshold; this usually means something has
                ! gone wrong..... need to check if forces need to be
                ! real at this stage even in the case of PAW formalism
                if (any(abs(aimag(block_cmplx))>pub_imag_thr)) then
                   call utils_abort('Error in routine nonsc_forces_extract_forces: &
                        &imaginary part of (complex) forces is above threshold')
                end if

                block_tr = block_tr + real(block_cmplx(ingwf,ingwf),kind=DP)
             end do

          else
             block_real(:,:) = 0.0_DP
             call sparse_get_block(block_real,nabla(ndim),iat,iat)

             do ingwf=1,ngwf_basis%num_on_atom(iat)
                block_tr = block_tr + block_real(ingwf,ingwf)
             end do

          end if

          nonsc_forces(ndim,par%orig_atom(iat)) = block_tr
       end do

       iat = iat + 1
    end do

    ! ars: deallocate memory
    ! agrecocmplx
    if (loc_cmplx) then
       deallocate(block_cmplx, stat=ierr)
       call utils_dealloc_check('nonsc_forces_extract_forces','block_cmplx',ierr)
    else
       deallocate(block_real, stat=ierr)
       call utils_dealloc_check('nonsc_forces_extract_forces','block_real',ierr)
    end if

    ! ars: add factors and distribute over the procs
    nonsc_forces(:,:) = nonsc_forces(:,:)/cell%weight
    call comms_reduce('SUM', nonsc_forces)


  end subroutine nonsc_forces_extract_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nonsc_forces_diagonal_blocks(ov_deriv, ngwf_basis, par)

    !======================================================================!
    ! This subroutine extracts the diagonal atom-blocks of ov_deriv and    !
    ! adds the transpose of themselves to calculate the correct derivative !
    ! of the overlap matrix wrt the centre of the NGWF spheres.            !
    !                                                                      !
    ! ***Note: in the vast majority of calculations, this routine will set !
    ! the diagonal blocks of ov_deriv to zero or almost zero.              !
    !======================================================================!
    !  Arguments:                                                          !
    !   ov_deriv (inout): deriv of the overlap wrt the centre of the sphere!
    !======================================================================!
    ! Originally written by Alvaro Ruiz Serrano in November 2009.          !
    ! Clean-up by Alvaro Ruiz Serrano in December 2010.                    !
    ! Modified to remove pub_par by Joseph Prentice, May 2018              !
    !======================================================================!


    use comms, only: pub_my_proc_id
    use constants, only : DP
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: PARAL_INFO
    use sparse, only: sparse_get_block, sparse_put_block, SPAM3
    use utils, only: utils_alloc_check, utils_dealloc_check


    implicit none

    ! ars: <<arguments>>
    type(SPAM3), intent(inout) :: ov_deriv(3)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(PARAL_INFO), intent(in) :: par

    ! ars: <<local variables>>
    real(kind=DP), allocatable :: block(:,:)
    integer :: atom, no_ngwfs, local_atom, ierr, ndim


    ! ars: loop over the atoms on this proc
    do local_atom = 1,par%num_atoms_on_proc(pub_my_proc_id)
       atom = par%first_atom_on_proc(pub_my_proc_id) + local_atom - 1

       ! ars: define block
       no_ngwfs = ngwf_basis%num_on_atom(atom)
       allocate(block(no_ngwfs, no_ngwfs), stat=ierr)
       call utils_alloc_check('nonsc_forces_diagonal_blocks','block',ierr)

       do ndim = 1, 3

          ! ars: init block
          block(:,:) = 0.0_DP

          ! ars: extract block corresponding to this atom
          call sparse_get_block(block, ov_deriv(ndim), atom, atom)

          ! ars: add the transpose
          block(:,:) = block(:,:) + transpose(block(:,:))

          ! ars: put block back in omega SPAM3
          call sparse_put_block(block, ov_deriv(ndim), atom, atom)

       end do

       ! ars: deallocate and start new atom
       deallocate(block, stat=ierr)
       call utils_dealloc_check('nonsc_forces_diagonal_blocks','block',ierr)

    end do


  end subroutine nonsc_forces_diagonal_blocks


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module nonsc_forces
