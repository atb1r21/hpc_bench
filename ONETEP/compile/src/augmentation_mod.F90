! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!      Nicholas D.M. Hine
!
!      Thomas Young Centre
!      Imperial College London
!      Exhibition Road
!      UK
!
!   Subsequent additions and modifications by: Gabriel C Constantinescu,
!      Jose M Escartin, Andrea Greco.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module augmentation

  use constants, only: DP, PI, stdout

  implicit none

  private

  ! Max number of augmentation functions on any atom
  integer :: max_proj_tot

  integer, parameter :: pad_pts = 4

  public :: aug_projector_denskern
  public :: aug_nonlocal_mat
  !public :: aug_nonlocal_commutator_mat
  public :: augmentation_overlap
  public :: augmentation_pos
  public :: augmentation_grad
  public :: augmentation_density_on_grid
  public :: augmentation_box_init
  public :: augmentation_box_exit
  public :: augmentation_density_forces
  public :: augmentation_screen_dij
  public :: aug_nl_calculate_forces
  public :: augmentation_FO_screen_dij ! gcc32
  public :: aug_FO_density_on_grid ! gcc32

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_box_init(aug_box,cell,grid,nl_projectors)

    !==================================================================!
    ! This subroutine calculates the size of the augmentation box used !
    ! in the calculation of nhat and other quantities.                 !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  grid (in) : the whole-cell grid for which the augmentation      !
    !              box is required.                                    !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30/06/10.                            !
    ! Moved to augmentation_mod by Nicholas Hine on 14/02/11.          !
    !==================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce, pub_on_root
    use constants, only: VERBOSE, stdout
    use fft_box, only: FFTBOX_INFO, fftbox_init, fftbox_find_size
    use geometry, only: geometry_magnitude
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_augbox_pref, pub_output_detail
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_banner

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(PROJECTOR_SET), intent(in) :: nl_projectors
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(inout) :: aug_box

    ! Local variables
    real(kind=DP) :: max_rad
    integer :: aug_box_n1, aug_box_n2, aug_box_n3
    logical :: extend(3)
    real(kind=DP) :: halo

    ! Find size of augmentation density box
    max_rad = maxval(nl_projectors%proj_max_radius(:))
    max_proj_tot = maxval(nl_projectors%species_num_proj(:))
    extend = (/.false.,.false.,.false./)
    halo = 0.0_DP

    call fftbox_find_size(aug_box_n1, aug_box_n2, aug_box_n3, &
       max_rad,1,extend,0.0_DP,pub_augbox_pref, &
       cell%a1,cell%a2,cell%a3,cell%b1,cell%b2,cell%b3, &
       geometry_magnitude(grid%da1),geometry_magnitude(grid%da2), &
       geometry_magnitude(grid%da3),grid%n1,grid%n2,grid%n3)

    ! Report size of aug box
    if (pub_on_root) then
       write(stdout,'(a)') ''
       write(stdout,'(a)') utils_banner('=', 'Charge Augmentation')
       write(stdout,'(3(a,i4))') '                           Aug box size: ',&
            aug_box_n1,' x',aug_box_n2,' x',aug_box_n3
       write(stdout,'(a)') repeat('=',80)
    end if

    call fftbox_init(aug_box,aug_box_n1,aug_box_n2,aug_box_n3,cell)

  end subroutine augmentation_box_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_box_exit(aug_box)

    !==================================================================!
    ! This subroutine deallocates storage for the augmentation box     !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 13/08/14.                            !
    !==================================================================!

    use fft_box, only: FFTBOX_INFO, fftbox_exit

    ! Arguments
    type(FFTBOX_INFO), intent(inout) :: aug_box

    call fftbox_exit(aug_box)

  end subroutine augmentation_box_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine aug_projector_denskern(proj_denskern,denskern,sp_overlap, &
       sp_overlap_ket)

    !====================================================================!
    ! This subroutine creates the reduced density matrix for each        !
    ! atomic site, in a diagonal matrix of size nproj x nproj. This      !
    ! is rho_ij in the PAW and USP formalisms, where \rho_{ij} is the    !
    ! occupancy of partial wave i,j.                                     !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  proj_denskern (inout) : The projector density kernel rho^ij       !
    !  denskern (in) : The NGWF density kernel K^ab                      !
    !  sp_overlap (in) : The NGWF-Projector overlap matrix <phi_a|p_i>.  !
    !--------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.             !
    ! Moved to projectors_mod by Nicholas Hine, 11/02/11.                !
    !====================================================================!

    use rundat, only: pub_num_spins, pub_imag_thr
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_product, &
         sparse_transpose, sparse_transpose_structure, sparse_take_real_part

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: proj_denskern(pub_num_spins)
    type(SPAM3), intent(in) :: sp_overlap
    type(SPAM3), optional, intent(in) :: sp_overlap_ket
    type(SPAM3), intent(in) :: denskern(pub_num_spins)

    ! Local Variables
    type(SPAM3) :: ps_overlap
    type(SPAM3) :: ksp
    integer :: is
    ! agrecocmplx
    logical :: loc_cmplx
    ! needed for sparse_product when using complex NGWFs
    type(SPAM3) :: proj_cmplx

    ! agrecocmplx
    loc_cmplx = denskern(1)%iscmplx

    ! Transpose <phi_a|p^i> to get <p^j|phi_b>
    ! tjz07: if there are two different NGWF sets, use the appropriate here
    call sparse_transpose_structure(ps_overlap%structure,sp_overlap)
    call sparse_create(ps_overlap, iscmplx=loc_cmplx)
    call sparse_transpose(ps_overlap,sp_overlap)
    if(present(sp_overlap_ket)) then
       call sparse_create(ksp,denskern(1),sp_overlap_ket)
    else
       call sparse_create(ksp,denskern(1),sp_overlap)
    endif

    ! agrecocmplx
    if (loc_cmplx) then
       call sparse_create(proj_cmplx,proj_denskern(1),iscmplx=loc_cmplx)
    end if

    ! Calculate <p^i|\hat{\rho}|p^j>, the projector density kernel:
    ! rho^ij = <p^i|phi_a>.(K^ab.<phi_b|p^j>)
    ! tjz07: Extend routine so it is valid for different NGWF sets
    ! agrecokpt: need to loop over k-points as well
    do is=1,pub_num_spins
       if(present(sp_overlap_ket)) then
          call sparse_product(ksp,denskern(is),sp_overlap_ket)
       else
          call sparse_product(ksp,denskern(is),sp_overlap)
       endif
       ! agrecocmplx: convert complex product to real proj_denskern
       if (loc_cmplx) then
          call sparse_product(proj_cmplx,ps_overlap,ksp)
          ! convert complex to real: use safe copy routine
          ! to check imaginary part of temp complex matrix
          !call sparse_copy(proj_denskern(is),proj_cmplx)
          call sparse_take_real_part(proj_denskern(is), &
               proj_cmplx, pub_imag_thr)
       ! standard real case
       else
          call sparse_product(proj_denskern(is),ps_overlap,ksp)
       end if
    end do

    ! agrecocmplx
    if (loc_cmplx) then
       call sparse_destroy(proj_cmplx)
    end if

    call sparse_destroy(ksp)
    call sparse_destroy(ps_overlap)

  end subroutine aug_projector_denskern


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_overlap(overlap,pseudo_sp,paw_sp,bra_proj_overlap, &
       ket_proj_overlap)

    !====================================================================!
    ! This subroutine augments an overlap matrix by adding the           !
    ! contribution from the overlap operator S = 1 + |p_i>O_ij<p_j|.     !
    ! The matrix returned is                                             !
    !      <bra_a|S|ket_b> = <bra_a|ket_b> + <bra_a|p_i>O_ij<p_j|ket_b>  !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  overlap (inout) : On entry: a general overlap matrix between two  !
    !      set of functions <bra_a|ket_b>.                               !
    !      On exit: the overlap matrix augmented with the projector part !
    !  bra_proj_overlap (in) : The Bra-Projector overlap matrix.         !
    !  ket_proj_overlap (in) : The Ket-Projector overlap matrix.         !
    !--------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 28/05/10.             !
    ! Moved to augmentation_mod by Nicholas Hine 14/02/11.               !
    ! Adapted to cope with augmenting overlaps between different sets of !
    ! functions in bras and kets.                                        !
    ! Modified to remove pub_par by Robert Charlton, 06/11/2018.         !
    !====================================================================!

    use pseudopotentials, only: PSEUDO_SPECIES, pseudo_aug_Q_matrix
    use paw, only: PAW_SPECIES, paw_projector_overlap
    use rundat, only: pub_paw, pub_usp
    use sparse, only: SPAM3, sparse_axpy, sparse_create, sparse_destroy, &
         sparse_product, sparse_transpose, &
         sparse_transpose_structure
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: overlap
    type(SPAM3), intent(in) :: bra_proj_overlap
    type(PSEUDO_SPECIES), intent(in) :: pseudo_sp(:)
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    type(SPAM3), intent(in), optional :: ket_proj_overlap

    ! Local Variables
    type(SPAM3) :: proj_overlap
    type(SPAM3) :: proj_ket_overlap
    type(SPAM3) :: bra_proj_overlap_O
    ! agrecocmplx
    logical :: loc_cmplx

    ! Check inputs
    if (present(ket_proj_overlap)) then
       if (ket_proj_overlap%iscmplx.neqv.bra_proj_overlap%iscmplx) then
          call utils_abort('Error in augmentation_overlap: bra and ket overlap &
               &matrices must be both real or both complex')
       end if
    end if

    ! agrecocmplx
    loc_cmplx = overlap%iscmplx

    ! Create temporary matrices to store <bra_a|p_i>, <bra_a|p_i> O_ij  and O_ij
    if (present(ket_proj_overlap)) then
       call sparse_transpose_structure(proj_ket_overlap%structure, &
            ket_proj_overlap)
    else
       call sparse_transpose_structure(proj_ket_overlap%structure, &
            bra_proj_overlap)
    end if
    ! agrecocmplx
    call sparse_create(proj_ket_overlap,iscmplx=loc_cmplx)
    call sparse_create(bra_proj_overlap_O,bra_proj_overlap)
    proj_overlap%structure = 'E'
    call sparse_create(proj_overlap,iscmplx=loc_cmplx)

    if (.not.present(ket_proj_overlap)) then
       ! If functions of the bras are the same as the functions of the kets, then just
       ! transpose <bra_a|p_i> to get <p_i|ket_a>
       call sparse_transpose(proj_ket_overlap,bra_proj_overlap)
    else
       ! If they are different, transpose the supplied <ket_a|p_i> to get <p_i|ket_a>
       call sparse_transpose(proj_ket_overlap,ket_proj_overlap)
    end if

    ! Get projector overlap matrix O_ij
    if (pub_paw) then
       call paw_projector_overlap(proj_overlap,paw_sp)
    else if (pub_usp) then
       call pseudo_aug_Q_matrix(proj_overlap,pseudo_sp)
    end if

    ! Create product <bra_a|p_i> O_ij
    call sparse_product(bra_proj_overlap_O,bra_proj_overlap,proj_overlap)

    ! Create NGWF-sized proj_overlap matrix <bra_a|p_i>O_ij<p_j|ket_b>
    call sparse_destroy(proj_overlap)
    call sparse_create(proj_overlap,overlap)
    call sparse_product(proj_overlap,bra_proj_overlap_O,proj_ket_overlap)

    ! Add it to the unaugmented overlap matrix <bra_a|ket_b>
    call sparse_axpy(overlap,proj_overlap,1.0_DP)

    ! Clean up temporary matrices
    call sparse_destroy(proj_overlap)
    call sparse_destroy(bra_proj_overlap_O)
    call sparse_destroy(proj_ket_overlap)

  end subroutine augmentation_overlap


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_pos(r_elements, proj_basis, pseudo_sp, paw_sp, &
       bra_proj_overlap, ket_proj_overlap, axis, first_order, direction, &
       weight)

    !====================================================================!
    ! This subroutine augments the matrix for the position operator in   !
    ! the NGWF representation, by adding to the r_elements matrix the    !
    ! contribution |p_i>R_ij<p_j| from the sphere terms.                 !
    ! where in PAW, R_ij = <phi_i|r|phi_j>-<tphi_i|r|tphi_j>.            !
    ! The matrices returned (one for each cartesian direction) are       !
    !   <bra_a|r|ket_b> = <bra_a|r|ket_b> + <bra_a|p_i>R_ij<p_j|ket_b>   !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  r_elements(3) (inout) : On entry: position operator matrices      !
    !      between two set of functions <bra_a|r|ket_b>.                 !
    !      On exit: the pos matrix augmented with the projector part     !
    !  proj_basis (in) : FUNC_BASIS type describing projectors.          !
    !  bra_proj_overlap (in) : The Bra-Projector overlap matrix.         !
    !  ket_proj_overlap (in) : The Ket-Projector overlap matrix.         !
    !  axis (in,optional)    : axis specifier (all axes if not present)  !
    !--------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 02/12/11.             !
    ! Modified to remove pub_par by Robert Charlton, 06/11/2018.         !
    !====================================================================!

    use comms, only: pub_my_proc_id
    use function_basis, only: FUNC_BASIS
    use pseudopotentials, only: PSEUDO_SPECIES, pseudo_aug_Q_matrix
    use parallel_strategy, only: PARAL_INFO
    use paw, only: PAW_SPECIES, paw_projector_overlap, paw_position_operator
    use rundat, only: pub_paw, pub_usp
    use sparse, only: SPAM3, sparse_axpy, sparse_create, sparse_destroy, &
         sparse_product, sparse_put_element, sparse_copy, sparse_get_par, &
         sparse_transpose, sparse_transpose_structure, sparse_get_element
    use utils, only: utils_abort, utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: r_elements(3)
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(SPAM3), intent(in) :: bra_proj_overlap
    type(PSEUDO_SPECIES), intent(in) :: pseudo_sp(:)
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    type(SPAM3), intent(in), optional :: ket_proj_overlap
    integer, intent(in), optional :: axis
    logical, intent(in), optional :: first_order
    integer, intent(in), optional :: direction(:)
    real(kind=DP), intent(in), optional :: weight(:)

    ! Local Variables
    type(SPAM3) :: proj_overlap
    type(SPAM3) :: r_sphere(3)
    type(SPAM3) :: proj_r_elements
    type(SPAM3) :: proj_ket_overlap
    type(SPAM3) :: bra_proj_overlap_R
    integer :: xyz, axmin, axmax
    integer :: loc_iproj, iproj, jproj
    integer :: iat, loc_iat
    real(kind=DP) :: R_atom(3), r_el, o_el
    ! agrecocmplx
    logical :: loc_cmplx, loc_first_order
    type(PARAL_INFO), pointer :: par

    ! Check inputs
    if (present(ket_proj_overlap)) then
       if (ket_proj_overlap%iscmplx.neqv.bra_proj_overlap%iscmplx) then
          call utils_abort('Error in augmentation_pos: bra and ket overlap &
               &matrices must be both real or both complex')
       end if
    end if
    if (present(axis)) then
       axmin = axis
       axmax = axis
    else
       axmin = 1
       axmax = 3
    end if

    ! rc2013: get parallel strategy, check arguments are compatible
    call sparse_get_par(par, bra_proj_overlap)
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &augmentation_pos: allocated parallel strategy is &
         &incompatible with paw_sp.')
    call utils_assert(par%num_pspecies == size(pseudo_sp), 'Error in &
         &augmentation_pos: allocated parallel strategy is &
         &incompatible with pseudo_sp.')
    ! jcap: only check these if present
    if (present(direction)) then
       call utils_assert(par%nat == size(direction), 'Error in augmentation_pos: &
            &allocated parallel strategy is incompatible with direction.')
    end if
    if (present(weight)) then
       call utils_assert(par%nat == size(weight), 'Error in augmentation_pos: &
            &allocated parallel strategy is incompatible with weight.')
    end if

    if (present(first_order)) then
       loc_first_order = first_order
    else
       loc_first_order = .false.
    end if

    ! agrecocmplx
    loc_cmplx = bra_proj_overlap%iscmplx

    ! Create temporary matrices to store <bra_a|p_i>, <bra_a|p_i> R_ij  and R_ij
    if (present(ket_proj_overlap)) then
       call sparse_transpose_structure(proj_ket_overlap%structure, &
            ket_proj_overlap)
    else
       call sparse_transpose_structure(proj_ket_overlap%structure, &
            bra_proj_overlap)
    end if
    call sparse_create(proj_ket_overlap,iscmplx=loc_cmplx)
    call sparse_create(bra_proj_overlap_R,bra_proj_overlap)
    proj_overlap%structure = 'E'
    ! agrecocmplx: these matrices are real anyway
    call sparse_create(proj_overlap)
    do xyz=axmin,axmax
       r_sphere(xyz)%structure = 'E'
       call sparse_create(r_sphere(xyz))
    end do

    if (.not.present(ket_proj_overlap)) then
       ! If functions of the bras are the same as the functions of the kets, then just
       ! transpose <bra_a|p_i> to get <p_i|ket_a>
       call sparse_transpose(proj_ket_overlap,bra_proj_overlap)
    else
       ! If they are different, transpose the supplied <ket_a|p_i> to get <p_i|ket_a>
       call sparse_transpose(proj_ket_overlap,ket_proj_overlap)
    end if

    ! Get projector position operator R_ij(1:3) and overlap matrix O_ij
    if (pub_paw) then
       call paw_projector_overlap(proj_overlap,paw_sp)
       ! ndmh: calculate sphere part of position operator
       call paw_position_operator(r_sphere,paw_sp)
    else if (pub_usp) then
       call pseudo_aug_Q_matrix(proj_overlap,pseudo_sp)
       call utils_abort('Error in augmentation_pos: USP position operator &
            &does not exist')
    end if

    ! ndmh: cycle over projectors on this proc, applying correction to
    ! ndmh: r_sphere to move it to the atom centre
    do loc_iproj=1,proj_basis%num_on_proc(pub_my_proc_id)
       iproj = loc_iproj + proj_basis%first_on_proc(pub_my_proc_id) - 1
       iat = proj_basis%atom_of_func(iproj)
       loc_iat = iat - par%first_atom_on_proc(pub_my_proc_id) + 1
       do jproj=proj_basis%first_on_atom(iat), &
            proj_basis%first_on_atom(iat)+proj_basis%num_on_atom(iat)-1

          ! Extract overlap element
          call sparse_get_element(o_el,proj_overlap,jproj,iproj)

          R_atom(1) = par%elements_on_proc(loc_iat)%centre%x
          R_atom(2) = par%elements_on_proc(loc_iat)%centre%y
          R_atom(3) = par%elements_on_proc(loc_iat)%centre%z

          ! Extract element from r_sphere and shift by R_atom*o_el
          ! ddor: get elements for one direction only if axis is specified
          do xyz=axmin,axmax
             call sparse_get_element(r_el,r_sphere(xyz),jproj,iproj)
             if (loc_first_order) then
                if (xyz .eq. direction(par%orig_atom(iat))) then
                   r_el = o_el * weight(par%orig_atom(iat))
                else
                   r_el = 0.0_DP
                end if
             else
                r_el = R_atom(xyz) * o_el + r_el
             end if
             call sparse_put_element(r_el,r_sphere(xyz),jproj,iproj)
          end do

       end do
    end do

    call sparse_create(proj_r_elements,r_elements(axmin))

    do xyz=axmin,axmax

       ! Create product <bra_a|p_i> R_ij
       call sparse_product(bra_proj_overlap_R, bra_proj_overlap, &
            r_sphere(xyz), allow_mix_types=loc_cmplx)

       ! Create NGWF-sized r_elements matrix <bra_a|p_i>R_ij<p_j|ket_b>
       call sparse_product(proj_r_elements,bra_proj_overlap_R,proj_ket_overlap)

       ! Add it to the unaugmented overlap matrix <bra_a|ket_b>
       call sparse_axpy(r_elements(xyz),proj_r_elements,1.0_DP)

    end do

    nullify(par)

    ! Clean up temporary matrices
    call sparse_destroy(proj_r_elements)
    do xyz=axmax,axmin,-1
       call sparse_destroy(r_sphere(xyz))
    end do
    call sparse_destroy(proj_overlap)
    call sparse_destroy(bra_proj_overlap_R)
    call sparse_destroy(proj_ket_overlap)

  end subroutine augmentation_pos


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_grad(grad_elements,pseudo_sp,paw_sp,bra_proj_overlap, &
       ket_proj_overlap,axis)

    !====================================================================!
    ! This subroutine augments the matrix for the position operator in   !
    ! the NGWF representation, by adding to the r_elements matrix the    !
    ! contribution |p_i>R_ij<p_j| from the sphere terms.                 !
    ! where in PAW, R_ij = <phi_i|r|phi_j>-<tphi_i|r|tphi_j>.            !
    ! The matrices returned (one for each cartesian direction) are       !
    !   <bra_a|r|ket_b> = <bra_a|r|ket_b> + <bra_a|p_i>R_ij<p_j|ket_b>   !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  grad_elements(3) (inout) : On entry: grad operator matrices       !
    !      between two set of functions <bra_a|nabla|ket_b>.             !
    !      On exit: the grad matrix augmented with the projector part    !
    !  bra_proj_overlap (in) : The Bra-Projector overlap matrix.         !
    !  ket_proj_overlap (in) : The Ket-Projector overlap matrix.         !
    !  axis (in,optional)    : axis specifier (all axes if not present)  !
    !--------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 02/12/11.             !
    ! Modified to remove pub_par by Robert Charlton, 06/11/2018.         !
    !====================================================================!

    use function_basis, only: FUNC_BASIS
    use pseudopotentials, only: PSEUDO_SPECIES, pseudo_aug_Q_matrix
    use paw, only: PAW_SPECIES, paw_projector_overlap, paw_grad_operator
    use rundat, only: pub_paw, pub_usp
    use sparse, only: SPAM3, sparse_axpy, sparse_create, sparse_destroy, &
         sparse_product, sparse_transpose, sparse_transpose_structure, &
         sparse_copy
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: grad_elements(3)
    type(SPAM3), intent(in) :: bra_proj_overlap
    type(PSEUDO_SPECIES), intent(in) :: pseudo_sp(:)
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    type(SPAM3), intent(in), optional :: ket_proj_overlap
    integer, intent(in), optional :: axis

    ! Local Variables
    type(SPAM3) :: proj_overlap
    type(SPAM3) :: grad_sphere(3)
    type(SPAM3) :: proj_grad_elements
    type(SPAM3) :: proj_ket_overlap
    type(SPAM3) :: bra_proj_overlap_grad
    integer :: xyz, axmin, axmax
    ! agrecocmplx
    logical :: loc_cmplx

    ! Check inputs
    if (present(ket_proj_overlap)) then
       if (ket_proj_overlap%iscmplx.neqv.bra_proj_overlap%iscmplx) then
          call utils_abort('Error in augmentation_grad: bra and ket overlap &
               &matrices must be both real or both complex')
       end if
    end if
    if (present(axis)) then
       axmin = axis
       axmax = axis
    else
       axmin = 1
       axmax = 3
    end if

    ! agrecocmplx
    loc_cmplx = bra_proj_overlap%iscmplx

    ! Create temporary matrices to store <bra_a|p_i>, <bra_a|p_i> grad_ij
    if (present(ket_proj_overlap)) then
       call sparse_transpose_structure(proj_ket_overlap%structure, &
            ket_proj_overlap)
    else
       call sparse_transpose_structure(proj_ket_overlap%structure, &
            bra_proj_overlap)
    end if
    call sparse_create(proj_ket_overlap,iscmplx=loc_cmplx)
    call sparse_create(bra_proj_overlap_grad,bra_proj_overlap)
    proj_overlap%structure = 'E'
    ! agrecocmplx: these matrices are real anyway
    call sparse_create(proj_overlap)
    do xyz=axmin,axmax
       grad_sphere(xyz)%structure = 'E'
       call sparse_create(grad_sphere(xyz))
    end do

    if (.not.present(ket_proj_overlap)) then
       ! If functions of the bras are the same as the functions of the kets,
       ! then just transpose <bra_a|p_i> to get <p_i|ket_a>
       call sparse_transpose(proj_ket_overlap,bra_proj_overlap)
    else
       ! If they are different, transpose the supplied <ket_a|p_i>
       ! to get <p_i|ket_a>
       call sparse_transpose(proj_ket_overlap,ket_proj_overlap)
    end if

    ! Get projector position operator grad_ij(1:3) and overlap matrix O_ij
    if (pub_paw) then
       call paw_projector_overlap(proj_overlap,paw_sp)
       ! ndmh: calculate sphere part of grad operator
       call paw_grad_operator(grad_sphere,paw_sp)
    else if (pub_usp) then
       call pseudo_aug_Q_matrix(proj_overlap,pseudo_sp)
       call utils_abort('Error in augmentation_grad: USP grad operator &
            &does not exist')
    end if

    call sparse_create(proj_grad_elements,grad_elements(axmin))

    do xyz=axmin,axmax

       ! Create product <bra_a|p_i> grad_ij
       call sparse_product(bra_proj_overlap_grad, bra_proj_overlap, &
            grad_sphere(xyz), allow_mix_types=loc_cmplx)

       ! Create NGWF-sized r_elements matrix <bra_a|p_i>grad_ij<p_j|ket_b>
       call sparse_product(proj_grad_elements,bra_proj_overlap_grad, &
            proj_ket_overlap)

       ! Add it to the unaugmented overlap matrix <bra_a|ket_b>
       call sparse_axpy(grad_elements(xyz),proj_grad_elements,1.0_DP)

    end do

    ! Clean up temporary matrices
    call sparse_destroy(proj_grad_elements)
    do xyz=axmax,axmin,-1
       call sparse_destroy(grad_sphere(xyz))
    end do
    call sparse_destroy(proj_overlap)
    call sparse_destroy(bra_proj_overlap_grad)
    call sparse_destroy(proj_ket_overlap)

  end subroutine augmentation_grad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine aug_nonlocal_mat(nonlocpot,dijhat,rho_ij,sp_overlap, &
       pseudo_sp,paw_sp,paw_sphere_energies,show_matrices)

    !==================================================================!
    ! This subroutine creates the PAW nonlocal matrix given by         !
    !  V^nl_ab = <phi_a|p_i> D_ij <p_j|phi_b>                          !
    ! where D_ij are the nonlocal energies given by                    !
    !  D_ij = \hat{D}_ij + D^1_ij - \tilde{D}^1_ij                     !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  nonlocpot (inout) : Nonlocal potential matrix V^nl_ab           !
    !  dijhat (in) : Augmentation functions screened with locpot       !
    !  rho_ij (in) : Projector density kernel rho_ij                   !
    !  sp_overlap (in) : NGWF-Projector overlap matrix <phi_a|p_i>     !
    !  paw_sphere_energies (inout) : Sphere energy contribution terms  !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.           !
    ! Modified for embedding by Joseph Prentice, June 2018             !
    !==================================================================!

    use comms, only: pub_on_root
    use constants, only: VERBOSE, paw_en_size, paw_en_dijhat
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: PARAL_INFO
    use paw, only: PAW_SPECIES, paw_nonlocal_energies, paw_show_atomblocks
    use pseudopotentials, only: PSEUDO_SPECIES, pseudo_get_dij
    use rundat, only: pub_debug_on_root, pub_num_spins, pub_paw, pub_usp, &
         pub_paw_output_detail
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_product, &
         sparse_trace, sparse_axpy, sparse_transpose, &
         sparse_transpose_structure, sparse_copy, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: nonlocpot(pub_num_spins)
    type(SPAM3), intent(in) :: dijhat(pub_num_spins)
    type(SPAM3), intent(in) :: rho_ij(pub_num_spins)
    type(SPAM3), intent(in) :: sp_overlap
    type(PSEUDO_SPECIES), intent(in) :: pseudo_sp(:)
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    real(kind=DP), intent(inout), optional :: paw_sphere_energies(paw_en_size)
    logical, intent(in), optional :: show_matrices

    ! Local Variables
    type(SPAM3), allocatable :: dij(:)
    type(SPAM3) :: sp_overlap_dij
    type(SPAM3) :: ps_overlap
    integer :: ierr
    integer :: is
    logical :: iscmplx
    logical :: loc_show_matrices
    type(PARAL_INFO), pointer :: par

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering aug_nonlocal_mat'

    ! Start Timer
    call timer_clock('aug_nonlocal_mat',1)

    ! Optional argument
    loc_show_matrices = .false.
    if (present(show_matrices)) loc_show_matrices = show_matrices

    ! rc2013: get parallel strategy from SPAM3
    call sparse_get_par(par, nonlocpot(1))
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &aug_nonlocal_mat: allocated parallel strategy is &
         &incompatible with paw_sp.')
    call utils_assert(par%num_pspecies == size(pseudo_sp), 'Error in &
         &aug_nonlocal_mat: allocated parallel strategy is &
         &incompatible with pseudo_sp.')

    ! Allocate arrays of matrices
    allocate(dij(pub_num_spins),stat=ierr)
    call utils_alloc_check('aug_nonlocal_mat','dij',ierr)
    do is=1,pub_num_spins
       call sparse_create(dij(is),rho_ij(is))
    end do

    ! Calculate the nonlocal energies D_ij, and the sphere energies if
    ! required
    if (pub_paw) then
       if (present(paw_sphere_energies)) then
          call paw_nonlocal_energies(dij,rho_ij,paw_sp,par, &
               paw_sphere_energies,loc_show_matrices)
       else
          call paw_nonlocal_energies(dij,rho_ij,paw_sp,par,&
               show_matrices=loc_show_matrices)
       end if

       do is=1,pub_num_spins
          if ((loc_show_matrices).and.(pub_paw_output_detail>=VERBOSE)) then
             if (pub_on_root) write(stdout,'(a,i4)') 'dijhat', is
             call paw_show_atomblocks(dijhat(is),paw_sp)
          end if
          call sparse_axpy(dij(is),dijhat(is),1.0_DP)
          if (present(paw_sphere_energies)) then
             paw_sphere_energies(paw_en_dijhat) = &
                  paw_sphere_energies(paw_en_dijhat) &
                  + sparse_trace(rho_ij(is),dijhat(is))
          end if
       end do

       if ((loc_show_matrices).and.(pub_paw_output_detail>=VERBOSE)) then
          do is=1,pub_num_spins
             if (pub_on_root) write(stdout,'(a,i4)') 'dij total', is
             call paw_show_atomblocks(dij(is),paw_sp)
          end do
       end if

    else if (pub_usp) then

       do is=1,pub_num_spins
          call pseudo_get_dij(dij(is),pseudo_sp)
       end do

    end if

    ! Create the matrix structures (set as real or complex depending on
    ! whether sp_overlap is real or complex) for sp and ps matrices
    iscmplx = sp_overlap%iscmplx
    call sparse_create(sp_overlap_dij,sp_overlap,iscmplx=iscmplx)
    call sparse_transpose_structure(ps_overlap%structure,sp_overlap)
    call sparse_create(ps_overlap,iscmplx=iscmplx)

    ! Create the transpose of the NGWF-projector overlap matrix
    call sparse_transpose(ps_overlap,sp_overlap)

    ! Loop over spins
    do is=1,pub_num_spins

       ! Calculate the matrix <NGWF_a|Proj_i> * D_ij
       call sparse_product(sp_overlap_dij, sp_overlap, dij(is), &
            allow_mix_types=iscmplx)

       ! Calculate the matrix \sum_i (<NGWF_a|Proj_i>D_ij<Proj_j|NGWF_b>)
       call sparse_product(nonlocpot(is),sp_overlap_dij,ps_overlap)

    end do

    call sparse_destroy(ps_overlap)
    call sparse_destroy(sp_overlap_dij)
    do is=pub_num_spins,1,-1
       call sparse_destroy(dij(is))
    end do
    deallocate(dij,stat=ierr)
    call utils_dealloc_check('aug_nonlocal_mat','dij',ierr)

    nullify(par)

    ! Start Timer
    call timer_clock('aug_nonlocal_mat',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving aug_nonlocal_mat'

  end subroutine aug_nonlocal_mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! rc2013: this subroutine is never used
!  subroutine aug_nonlocal_commutator_mat(nonlocpot_com, proj_basis, &
!       dijhat, rho_ij,  nl_projectors, ngwf_basis, ngwfs_on_grid, &
!       sp_overlap, fftbox, cell, pseudo_sp, paw_sp, delta_in)
!
!    !==================================================================!
!    ! This subroutine calculates the commutator between the nonlocal   !
!    ! potential and the position operator for the 3 Cartesian          !
!    ! directions.                                                      !
!    !------------------------------------------------------------------!
!    ! Arguments:                                                       !
!    !  nonlocpot (inout) : Nonlocal potential, position op commutator  !
!    !  dijhat (in) : Augmentation functions screened with locpot       !
!    !  rho_ij (in) : Projector density kernel rho_ij                   !
!    !  sp_overlap (in) : NGWF-Projector overlap matrix <phi_a|p_i>     !
!    !------------------------------------------------------------------!
!    ! This subroutine was written by Nicholas Hine 02/12/11.           !
!    !==================================================================!
!
!    use datatypes, only: FUNCTIONS
!    use cell_grid, only: GRID_INFO
!    use fft_box, only: FFTBOX_INFO
!    use function_basis, only: FUNC_BASIS
!    use parallel_strategy, only: PARAL_INFO
!    use paw, only: PAW_SPECIES, paw_nonlocal_energies
!    use projectors, only: PROJECTOR_SET, projectors_commutator
!    use pseudopotentials, only: PSEUDO_SPECIES, pseudo_get_dij
!    use rundat, only: pub_paw, pub_usp, pub_num_spins, pub_debug_on_root
!    use simulation_cell, only: CELL_INFO
!    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_axpy, &
!         sparse_get_par
!    use timer, only: timer_clock
!    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
!         utils_assert
!
!    implicit none
!
!    ! Arguments
!    type(SPAM3), intent(inout) :: nonlocpot_com(pub_num_spins)
!    type(FUNC_BASIS), intent(in) :: proj_basis
!    type(FUNC_BASIS), intent(in) :: ngwf_basis
!    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
!    type(PROJECTOR_SET), intent(inout) :: nl_projectors
!    type(SPAM3), intent(in) :: dijhat(pub_num_spins)
!    type(SPAM3), intent(in) :: rho_ij(pub_num_spins)
!    type(SPAM3), intent(in) :: sp_overlap
!    real(kind=DP), intent(in) :: delta_in ! finite difference shift
!    type(PSEUDO_SPECIES), intent(in) :: pseudo_sp(:)
!    type(PAW_SPECIES), intent(inout) :: paw_sp(:)
!    type(CELL_INFO), intent(in) :: cell
!    type(FFTBOX_INFO), intent(in) :: fftbox
!
!    ! Local Variables
!    type(SPAM3), allocatable :: dij(:)
!    integer :: ierr
!    integer :: is
!    type(PARAL_INFO), pointer :: par
!
!    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
!         &aug_nonlocal_commutator_mat'
!
!    ! Start Timer
!    call timer_clock('aug_nonlocal_commutator_mat',1)
!
!    ! rc2013: get parallel strategy from SPAM3
!    call sparse_get_par(par, nonlocpot_com(1))
!    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
!         &aug_nonlocal_mat: allocated parallel strategy is &
!         &incompatible with paw_sp.')
!    call utils_assert(par%num_pspecies == size(pseudo_sp), 'Error in &
!         &aug_nonlocal_mat: allocated parallel strategy is &
!         &incompatible with pseudo_sp.')
!
!    ! Allocate arrays of matrices
!    allocate(dij(pub_num_spins),stat=ierr)
!    call utils_alloc_check('aug_nonlocal_commutator_mat','dij',ierr)
!    do is=1,pub_num_spins
!       call sparse_create(dij(is),rho_ij(is))
!    end do
!
!    ! Calculate the nonlocal energies D_ij
!    if (pub_paw) then
!       call paw_nonlocal_energies(dij,rho_ij,paw_sp,par)
!       do is=1,pub_num_spins
!          call sparse_axpy(dij(is),dijhat(is),1.0_DP)
!       end do
!    else if (pub_usp) then
!       do is=1,pub_num_spins
!          call pseudo_get_dij(dij(is),pseudo_sp)
!       end do
!    end if
!
!    do is=1,pub_num_spins,1
!       call projectors_commutator(nonlocpot_com(is), proj_basis, &
!            ngwf_basis, ngwfs_on_grid, sp_overlap, nl_projectors, &
!            cell, fftbox, delta_in, dij(is))
!    end do
!
!    ! Deallocate nonlocal energies
!    do is=pub_num_spins,1,-1
!       call sparse_destroy(dij(is))
!    end do
!    deallocate(dij,stat=ierr)
!    call utils_dealloc_check('aug_nonlocal_commutator_mat','dij',ierr)
!
!    nullify(par)
!
!    ! Start Timer
!    call timer_clock('aug_nonlocal_commutator_mat',2)
!
!    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
!         &aug_nonlocal_commutator_mat'
!
!  end subroutine aug_nonlocal_commutator_mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_density_on_grid(nhat_den_grad,grid,cell,pseudo_sp, &
       paw_sp,aug_box,denskern,sp_overlap,sp_overlap_ket)

    !==================================================================!
    ! This subroutine creates the compensation density \hat{n}(r) on   !
    ! the simulation cell fine grid, and also the gradient of the      !
    ! augmentation density in each cartesian direction.                !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  grid (in) : Grid definition for the fine grid.                  !
    !  nhat (inout) : The compensation density \hat{n}(r) on the grid. !
    !  denskern (in) : The NGWF density kernel K^ab.                   !
    !  sp_overlap (in) : The NGWF-Projector overlap matrix <phi_a|p_i>.!
    !  sp_overlap_ket (in) :: NGWF-Projector overlap for second NGWF   !
    !    species if the bra and ket species are different.             !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.           !
    ! Moved to projectors_mod by Nicholas Hine on 11/02/11.            !
    ! Modified to remove pub_par by Robert Charlton, 06/11/2018.       !
    !==================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_deposit_box, &
         cell_grid_box_start_wrt_atom
    use comms, only: comms_reduce, pub_on_root, pub_my_proc_id
    use constants, only: max_spins, VERBOSE
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use geometry, only: geometry_magnitude
    use parallel_strategy, only: PARAL_INFO
    use paw, only: PAW_SPECIES, paw_atom_aug_den
    use pseudopotentials, only: PSEUDO_SPECIES, pseudo_atom_aug_den
    use rundat, only: pub_paw, pub_usp, pub_output_detail, pub_aug_den_dim, &
         pub_aug_funcs_recip, pub_debug_on_root, pub_num_spins
!$  use rundat, only: pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_get_block, &
         sparse_get_par
    use spherical_wave, only: sw_init
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_assert
    use xc, only: pub_xc_gradient_corrected

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in) :: cell
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: nhat_den_grad(grid%ld1, &
         grid%ld2, grid%max_slabs12, pub_num_spins, 0:pub_aug_den_dim)
    ! agrecokpt: change this to full denskern(pub_num_spins, pub_num_kpoints)
    type(SPAM3), intent(in) :: denskern(pub_num_spins)
    type(SPAM3), intent(in) :: sp_overlap
    type(FFTBOX_INFO), intent(in) :: aug_box
    type(PSEUDO_SPECIES), intent(in) :: pseudo_sp(:)
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    type(SPAM3), optional, intent(in) :: sp_overlap_ket ! tjz07: only used if bra
      ! and ket are corresponding to different NGWF species as in lr_tddft

    ! Local Variables
    type(SPAM3), allocatable :: rho_ij(:)
    character(20) :: fmt,tmp
    integer :: iat, loc_iat
    integer :: ierr
    integer :: is
    integer :: isp
    integer :: box_n1,box_n2,box_n3
    integer :: box_start1,box_start2,box_start3
    integer :: cart
    integer :: lmax
    real(kind=DP) :: total_nhat(max_spins)
    real(kind=DP) :: total_nhat_targ(max_spins)
    real(kind=DP),allocatable :: rho_ij_block(:,:,:)
    real(kind=DP), allocatable :: atom_nhat(:,:,:,:)
    real(kind=DP), allocatable :: atom_grad_nhat(:,:,:,:,:)
    real(kind=DP), allocatable :: atom_aug_func(:,:,:)
    real(kind=DP), allocatable :: atom_aug_func_grad(:,:,:,:)
    real(kind=DP), allocatable :: buffer(:,:,:)
    complex(kind=DP), allocatable :: atom_nhat_recip(:,:,:,:)
    complex(kind=DP), allocatable :: atom_grad_nhat_recip(:,:,:,:,:)
    complex(kind=DP), allocatable :: atom_aug_func_recip(:,:,:)
    complex(kind=DP), allocatable :: atom_aug_func_grad_recip(:,:,:,:)
    logical :: i_have_box
    type(PARAL_INFO), pointer :: par

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &augmentation_density_on_grid'

    ! Start Timer
    call timer_clock('augmentation_density_on_grid',1)

    ! rc2013: get parallel strategy, check arguments are compatible
    call sparse_get_par(par, sp_overlap)
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &augmentation_pos: allocated parallel strategy is &
         &incompatible with paw_sp.')
    call utils_assert(par%num_pspecies == size(pseudo_sp), 'Error in &
         &augmentation_pos: allocated parallel strategy is &
         &incompatible with pseudo_sp.')

    ! Find size of box
    box_n1 = aug_box%total_pt1
    box_n2 = aug_box%total_pt2
    box_n3 = aug_box%total_pt3

    ! agrecokpt: TO-DO: check/implement k-points here as well
    ! agrecocmplx: routine should now be compatible with complex NGWFs
    !call utils_assert(loc_cmplx.eqv..false., &
    !     'Subroutine augmentation_density_on_grid not ready yet for &
    !     & complex NGWFs.')

    ! Ensure spherical waves have been initialised to high enough lmax
    ! (must be before OMP region starts to avoid race condition)
    if (pub_paw) then
       lmax = maxval(paw_sp(:)%lmax)
       if (lmax>3) call sw_init(lmax+1,1)
    end if

    ! Create projector density kernel
    allocate(rho_ij(pub_num_spins),stat=ierr)
    call utils_alloc_check('augmentation_density_on_grid','rho_ij',ierr)
    do is=1,pub_num_spins
       ! agrecocmplx: rho_ij is real even when using complex NGWFs
       rho_ij(is)%structure = 'E'
       call sparse_create(rho_ij(is))
    end do

    if(present(sp_overlap_ket)) then
       call aug_projector_denskern(rho_ij,denskern,sp_overlap,sp_overlap_ket)
    else
       call aug_projector_denskern(rho_ij,denskern,sp_overlap)
    endif

    total_nhat = 0.0_DP
    total_nhat_targ = 0.0_DP

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(rho_ij_block,atom_nhat,atom_grad_nhat,buffer, &
!$OMP      atom_aug_func_recip,atom_aug_func_grad_recip,atom_nhat_recip, &
!$OMP      atom_grad_nhat_recip,atom_aug_func,atom_aug_func_grad, &
!$OMP      loc_iat,iat,isp,is,ierr,box_start1,box_start2,box_start3, &
!$OMP      cart,i_have_box) &
!$OMP SHARED(pub_my_proc_id,grid,cell,par,max_proj_tot,pub_num_spins,paw_sp, &
!$OMP      pseudo_sp,aug_box,box_n1,box_n2,box_n3,pub_aug_funcs_recip, &
!$OMP      pub_threads_num_fftboxes,rho_ij,pub_paw,pub_usp,nhat_den_grad, &
!$OMP      pub_xc_gradient_corrected) &
!$OMP REDUCTION (+:total_nhat,total_nhat_targ)

    ! Allocate workspace
    allocate(rho_ij_block(max_proj_tot,max_proj_tot, &
         pub_num_spins),stat=ierr)
    call utils_alloc_check('augmentation_density_on_grid','rho_ij_block',ierr)
    allocate(atom_nhat(box_n1,box_n2,box_n3,pub_num_spins), &
         stat=ierr)
    call utils_alloc_check('augmentation_density_on_grid','atom_nhat',ierr)
    allocate(atom_grad_nhat(box_n1,box_n2,box_n3,pub_num_spins,3), &
         stat=ierr)
    call utils_alloc_check('augmentation_density_on_grid','atom_grad_nhat', &
         ierr)
    allocate(buffer(box_n1,box_n2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('augmentation_density_on_grid','buffer',ierr)

    if (pub_aug_funcs_recip) then
       allocate(atom_aug_func_recip(box_n1,box_n2,box_n3),stat=ierr)
       call utils_alloc_check('augmentation_density_on_grid', &
            'atom_aug_func_recip',ierr)
       allocate(atom_aug_func_grad_recip(box_n1,box_n2,box_n3,3),stat=ierr)
       call utils_alloc_check('augmentation_density_on_grid', &
            'atom_aug_func_grad_recip',ierr)
       allocate(atom_nhat_recip(box_n1,box_n2,box_n3,pub_num_spins), &
            stat=ierr)
       call utils_alloc_check('augmentation_density_on_grid', &
            'atom_nhat_recip',ierr)
       allocate(atom_grad_nhat_recip(box_n1,box_n2,box_n3,pub_num_spins, &
            3),stat=ierr)
       call utils_alloc_check('augmentation_density_on_grid', &
            'atom_grad_nhat_recip',ierr)
    else
       allocate(atom_aug_func(box_n1,box_n2,box_n3),stat=ierr)
       call utils_alloc_check('augmentation_density_on_grid', &
            'atom_aug_func',ierr)
       allocate(atom_aug_func_grad(box_n1,box_n2,box_n3,3),stat=ierr)
       call utils_alloc_check('augmentation_density_on_grid', &
            'atom_aug_func_grad',ierr)
    end if

    !if (pub_on_root) write(stdout,'(a)')' iat lup mup ipt ipw  li  mi jpt &
    !     &jpw  lj  mj            nLij             Gij           rhoij         &
    !     &    qij             qLM'
!$OMP DO
    do loc_iat=1,par%max_atoms_on_proc

       if (loc_iat<=par%num_atoms_on_proc(pub_my_proc_id)) then
          iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
          isp = par%elements_on_proc(loc_iat)%pspecies_number

          ! Find where box for this atom is located in simulation cell
          call cell_grid_box_start_wrt_atom( &
               box_start1, box_start2, box_start3, &
               par%elements_on_proc(loc_iat)%centre, box_n1, box_n2, box_n3, &
               grid,cell)

          ! Get block of \rho_ij for this atom
          ! agrecocmplx
          do is=1,pub_num_spins
             call sparse_get_block(rho_ij_block(:,:,is),rho_ij(is),iat,iat)
          end do

          ! Call appropriate routine to generate the augmentation density for
          ! this atom in the augmentation box using rhoij
          if (pub_paw) then
             if (pub_aug_funcs_recip) then
                call paw_atom_aug_den(atom_nhat,atom_grad_nhat,total_nhat, &
                     total_nhat_targ,rho_ij_block, &
                     isp,par%elements_on_proc(loc_iat)%centre,grid,cell,paw_sp, &
                     aug_box,box_start1,box_start2,box_start3, &
                     atom_aug_func_recip=atom_aug_func_recip, &
                     atom_aug_func_grad_recip=atom_aug_func_grad_recip, &
                     atom_nhat_recip=atom_nhat_recip, &
                     atom_grad_nhat_recip=atom_grad_nhat_recip)
             else
                call paw_atom_aug_den(atom_nhat,atom_grad_nhat,total_nhat, &
                     total_nhat_targ,rho_ij_block, &
                     isp,par%elements_on_proc(loc_iat)%centre,grid,cell,paw_sp, &
                     aug_box,box_start1,box_start2,box_start3, &
                     atom_aug_func_real=atom_aug_func, &
                     atom_aug_func_grad_real=atom_aug_func_grad)
             end if
          else if (pub_usp) then
             call pseudo_atom_aug_den(atom_nhat,atom_grad_nhat,atom_aug_func, &
                  atom_aug_func_grad,total_nhat,total_nhat_targ,rho_ij_block, &
                  isp,par%elements_on_proc(loc_iat)%centre,grid,cell,pseudo_sp, &
                  box_n1,box_n2,box_n3,box_start1,box_start2,box_start3)
          end if

          i_have_box = .true.
       else
          ! Nothing to deposit on this proc
          i_have_box = .false.
       end if

       ! Deposit this box to the simulation cell if present, or just wait for
       ! data from other procs if no box
!$OMP CRITICAL
       do is=1,pub_num_spins
          call cell_grid_deposit_box(nhat_den_grad(:,:,:,is,0), &
               atom_nhat(:,:,:,is), buffer, grid, &
               box_n1, box_n2, box_n3, box_n1, box_n2, &
               box_start1, box_start2, box_start3, i_have_box, .false.)
          ! Deposit gradient to whole-cell arrays
          if (pub_xc_gradient_corrected) then
             do cart=1,3
                call cell_grid_deposit_box(nhat_den_grad(:,:,:,is,cart), &
                     atom_grad_nhat(:,:,:,is,cart), buffer, grid, &
                     box_n1, box_n2, box_n3, box_n1, box_n2, &
                     box_start1, box_start2, box_start3, i_have_box, .false.)
             end do
          end if
       end do
!$OMP END CRITICAL
    end do
!$OMP END DO

    ! Deallocate temporary arrays and matrices
    if (pub_aug_funcs_recip) then
       deallocate(atom_grad_nhat_recip,stat=ierr)
       call utils_dealloc_check('augmentation_density_on_grid', &
            'atom_grad_nhat_recip',ierr)
       deallocate(atom_nhat_recip,stat=ierr)
       call utils_dealloc_check('augmentation_density_on_grid', &
            'atom_nhat_recip',ierr)
       deallocate(atom_aug_func_grad_recip,stat=ierr)
       call utils_dealloc_check('augmentation_density_on_grid', &
            'atom_aug_func_grad_recip',ierr)
       deallocate(atom_aug_func_recip,stat=ierr)
       call utils_dealloc_check('augmentation_density_on_grid', &
            'atom_aug_func_recip',ierr)
    else
       deallocate(atom_aug_func_grad,stat=ierr)
       call utils_dealloc_check('augmentation_density_on_grid', &
            'atom_aug_func_grad',ierr)
       deallocate(atom_aug_func,stat=ierr)
       call utils_dealloc_check('augmentation_density_on_grid', &
            'atom_aug_func',ierr)
    end if

    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('augmentation_density_on_grid','buffer',ierr)
    deallocate(atom_grad_nhat,stat=ierr)
    call utils_dealloc_check('augmentation_density_on_grid','atom_grad_nhat', &
         ierr)
    deallocate(atom_nhat,stat=ierr)
    call utils_dealloc_check('augmentation_density_on_grid','atom_nhat',ierr)
    deallocate(rho_ij_block,stat=ierr)
    call utils_dealloc_check('augmentation_density_on_grid','rho_ij_block',ierr)

!$OMP END PARALLEL

    ! Get sum of compensation density over all procs, compare to target
    call comms_reduce('SUM',total_nhat(:))
    call comms_reduce('SUM',total_nhat_targ(:))
    if (any(abs(total_nhat_targ(1:pub_num_spins) &
         - total_nhat(1:pub_num_spins)) > 1.0d-15)) then
       if (pub_on_root.and.(pub_output_detail>=VERBOSE).and. &
            (.not.pub_aug_funcs_recip)) then
          write(tmp,'(i5)') pub_num_spins
          write(fmt,'(3a)') '(a,',trim(adjustl(tmp)),'f14.8)'
          write(stdout,fmt) 'Total Compensation Charge Target     &
               &    : ',total_nhat_targ(1:pub_num_spins)
          write(stdout,fmt) 'Total Compensation Charge on Regular &
               &Grid: ',total_nhat(1:pub_num_spins)
       end if
    end if

    do is=pub_num_spins,1,-1
       call sparse_destroy(rho_ij(is))
    end do

    nullify(par)

    deallocate(rho_ij,stat=ierr)
    call utils_dealloc_check('augmentation_density_on_grid','rho_ij',ierr)

    ! Stop Timer
    call timer_clock('augmentation_density_on_grid',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving augmentation_density_on_grid'

  end subroutine augmentation_density_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_screen_dij(dij,locpot,aug_box,cell,grid,pseudo_sp, &
       paw_sp)

    !=====================================================================!
    ! This subroutine calculates the contribution from the interaction of !
    ! the augmentation charge with the effective potential                !
    !     \hat{D}_ij = \sum_LM \int \tilde{v}_eff(r) \hat{Q}_ij^LM (r) dr !
    ! to the nonlocal term Dij of the Hamiltonian, where \tilde{v}_eff(r) !
    ! is the effective potential resulting from the smooth part of the    !
    ! density, v_H[\tilde{n}+\hat{n}+\tilde{n}_Zc]                        !
    !         +v_xc[\tilde{n}+\hat{n}+\tilde{n}_c]                        !
    !---------------------------------------------------------------------!
    ! Arguments:                                                          !
    !   dij(inout)    : Compensation density contribution to nonlocal     !
    !                   energy term \hat{D}_ij.                           !
    !   locpot(in)    : Local potential on fine grid.                     !
    !   grid(in)      : Grid definition for fine grid.                    !
    ! @docme: cell, aug_box, pseudo_sp, paw_sp.                           !
    !---------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30 May 2010.                            !
    !=====================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_extract_box, &
         cell_grid_box_start_wrt_atom
    use comms, only: pub_my_proc_id
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use geometry, only: geometry_magnitude
    use parallel_strategy, only: PARAL_INFO
    use paw, only: PAW_SPECIES, paw_atom_aug_integrals
    use pseudopotentials, only: PSEUDO_SPECIES, pseudo_atom_aug_integrals
    use rundat, only: pub_paw, pub_usp, pub_aug_funcs_recip, pub_debug_on_root,&
         pub_debug
!$  use rundat, only: pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_get_block, sparse_put_block, sparse_get_par
    use spherical_wave, only: sw_init
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_sanity_check, utils_assert

    implicit none

    ! Arguments
    type(CELL_INFO),intent(in) :: cell
    type(GRID_INFO),intent(in) :: grid
    type(FFTBOX_INFO),intent(in) :: aug_box
    type(SPAM3),intent(inout) :: dij(:)
    real(kind=DP),intent(in) :: locpot(grid%ld1,grid%ld2, &
         grid%max_slabs12,size(dij))
    type(PSEUDO_SPECIES), intent(in) :: pseudo_sp(:)
    type(PAW_SPECIES), intent(in) :: paw_sp(:)

    ! Local Variables
    integer :: loc_iat, iat
    integer :: isp
    integer :: is
    integer :: ierr
    integer :: box_n1
    integer :: box_n2
    integer :: box_n3
    integer :: box_start1
    integer :: box_start2
    integer :: box_start3
    logical :: i_need_box
    integer :: num_spins
    integer :: lmax
    real(kind=DP), allocatable :: dij_at(:,:,:)
    real(kind=DP), allocatable :: locpot_box(:,:,:,:)
    real(kind=DP), allocatable :: atom_aug_func(:,:,:)
    real(kind=DP), allocatable :: buffer(:,:,:)
    complex(kind=DP), allocatable :: atom_aug_func_recip(:,:,:)
    complex(kind=DP), allocatable :: locpot_box_recip(:,:,:,:)
    type(PARAL_INFO), pointer :: par

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &augmentation_screen_dij'

    ! Start Timer
    call timer_clock('augmentation_screen_dij',1)

    if (size(locpot,4)/=size(dij,1)) then
       call utils_abort('Error in augmentation_screen_dij: Inconsistent &
            &sizes of locpot and dij')
    end if
    num_spins = size(dij,1)
    isp = 0

    do is = 1, num_spins
       call utils_sanity_check(locpot(:,:,:,is),'locpot')
    end do

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, dij(1))

    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &augmentation_screen_dij: allocated parallel strategy is &
         &incompatible with paw_sp.')
    call utils_assert(par%num_pspecies == size(pseudo_sp), 'Error in &
         &augmentation_screen_dij: allocated parallel strategy is &
         &incompatible with pseudo_sp.')

    ! Find size of box
    box_n1 = aug_box%total_pt1
    box_n2 = aug_box%total_pt2
    box_n3 = aug_box%total_pt3

    ! Ensure spherical waves have been initialised to high enough lmax
    ! (must be before OMP region starts to avoid race condition)
    if (pub_paw) then
       lmax = maxval(paw_sp(:)%lmax)
       if (lmax>3) call sw_init(lmax+1,1)
    end if

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(dij_at,locpot_box,buffer,atom_aug_func_recip,locpot_box_recip, &
!$OMP      atom_aug_func,loc_iat,iat,isp,is,ierr,box_start1,box_start2, &
!$OMP      box_start3,i_need_box) &
!$OMP SHARED (pub_my_proc_id,grid,cell,par,aug_box,max_proj_tot,locpot,num_spins, &
!$OMP      box_n1,box_n2,box_n3,pub_aug_funcs_recip,pub_threads_num_fftboxes, &
!$OMP      dij,pseudo_sp,paw_sp,pub_paw,pub_usp)

    ! Allocate temporary arrays
    allocate(dij_at(max_proj_tot,max_proj_tot, &
         num_spins),stat=ierr)
    call utils_alloc_check('augmentation_screen_dij','dij_at',ierr)
    allocate(locpot_box(box_n1,box_n2,box_n3,num_spins),stat=ierr)
    call utils_alloc_check('augmentation_screen_dij','locpot_box',ierr)
    allocate(buffer(box_n1,box_n2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('augmentation_screen_dij','buffer',ierr)
    if (pub_aug_funcs_recip) then
       allocate(atom_aug_func_recip(box_n1,box_n2,box_n3),stat=ierr)
       call utils_alloc_check('augmentation_screen_dij','atom_aug_func_recip',ierr)
       allocate(locpot_box_recip(box_n1,box_n2,box_n3,num_spins),stat=ierr)
       call utils_alloc_check('augmentation_screen_dij','locpot_box_recip',ierr)
    else
       allocate(atom_aug_func(box_n1,box_n2,box_n3),stat=ierr)
       call utils_alloc_check('augmentation_screen_dij','atom_aug_func',ierr)
    end if
    locpot_box(:,:,:,:) = 0.0_DP

    ! Loop over atoms
!$OMP DO
    do loc_iat=1,par%max_atoms_on_proc

       ! Only need to extract if there is an atom left on this proc
       if (loc_iat<=par%num_atoms_on_proc(pub_my_proc_id)) then
          iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
          isp = par%elements_on_proc(loc_iat)%pspecies_number

          ! Find where box for this atom is located in simulation cell
          call cell_grid_box_start_wrt_atom( &
               box_start1, box_start2, box_start3, &
               par%elements_on_proc(loc_iat)%centre, box_n1, box_n2, box_n3, &
               grid,cell)

          i_need_box = .true.
       else
          i_need_box = .false.
       end if

       ! Extract tightbox of data from effective potential over simulation
       ! cell for this atom
!$OMP CRITICAL
       do is=1,num_spins
          call cell_grid_extract_box(locpot_box(:,:,:,is), &
               buffer, locpot(:,:,:,is), grid, &
               box_n1, box_n2, box_n3, box_n1, box_n2, &
               box_start1, box_start2, box_start3, i_need_box, .false.)
       end do  ! is
!$OMP END CRITICAL
       ! Only need to take overlap if there is an atom left on this proc
       if (loc_iat<=par%num_atoms_on_proc(pub_my_proc_id)) then

          do is=1,num_spins
             call utils_sanity_check(locpot_box(:,:,:,is),'locpot_box')
          end do

          ! Get block of dij from SPAM3 matrix
          !do is=1,num_spins
          !   call sparse_get_block(dij_at(:,:,is),dij(is),iat,iat)
          !end do  ! is
          dij_at(:,:,:) = 0.0_DP

          if (pub_aug_funcs_recip) then
             do is=1,num_spins
                locpot_box_recip(:,:,:,is) = locpot_box(:,:,:,is)
                call fourier_apply_box('F','F',locpot_box_recip(:,:,:,is), &
                     box=aug_box)
             end do
          end if

          ! gcc32: we actually need v(-G), which for purely real potentials
          ! is the same as v(G)^*
          locpot_box_recip = conjg(locpot_box_recip) / &
               real(box_n1*box_n2*box_n3,kind=DP)

          ! Call appropriate routine to calculate the integral of the
          ! augmentation function for this atom with the local potential
          ! in the augmentation box
          if (pub_paw) then
             if (pub_aug_funcs_recip) then
                call paw_atom_aug_integrals(dij_at,num_spins,isp, &
                     par%elements_on_proc(loc_iat)%centre,grid,cell,paw_sp, &
                     box_n1,box_n2,box_n3,box_start1,box_start2,box_start3, &
                     locpot_box_recip=locpot_box_recip, &
                     atom_aug_func_recip=atom_aug_func_recip)
             else
                call paw_atom_aug_integrals(dij_at,num_spins,isp, &
                     par%elements_on_proc(loc_iat)%centre,grid,cell,paw_sp, &
                     box_n1,box_n2,box_n3,box_start1,box_start2,box_start3, &
                     locpot_box,atom_aug_func)
             end if
          else if (pub_usp) then
             call pseudo_atom_aug_integrals(locpot_box,atom_aug_func, &
                  dij_at,num_spins,isp,par%elements_on_proc(loc_iat)%centre, &
                  grid,cell,pseudo_sp,box_n1,box_n2,box_n3, &
                  box_start1,box_start2,box_start3)
          end if

          ! Put block of dij into SPAM3 matrix
          do is=1,num_spins
             call sparse_put_block(dij_at(:,:,is),dij(is),iat,iat)
          end do  ! is

       end if  ! loc_iat < loc_nat
    end do  ! loc_iat
!$OMP END DO

    ! Deallocate temporary arrays
    if (pub_aug_funcs_recip) then
       deallocate(locpot_box_recip,stat=ierr)
       call utils_dealloc_check('augmentation_screen_dij','locpot_box_recip',ierr)
       deallocate(atom_aug_func_recip,stat=ierr)
       call utils_dealloc_check('augmentation_screen_dij','atom_aug_func_recip',ierr)
    else
       deallocate(atom_aug_func,stat=ierr)
       call utils_dealloc_check('augmentation_screen_dij','atom_aug_func',ierr)
    end if
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('augmentation_screen_dij','buffer',ierr)
    deallocate(locpot_box,stat=ierr)
    call utils_dealloc_check('augmentation_screen_dij','locpot_box',ierr)
    deallocate(dij_at,stat=ierr)
    call utils_dealloc_check('augmentation_screen_dij','dij_at',ierr)

!$OMP END PARALLEL

    ! Stop Timer
    call timer_clock('augmentation_screen_dij',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving augmentation_screen_dij'

  end subroutine augmentation_screen_dij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine augmentation_FO_screen_dij(dij, locpot, FO_locpot, aug_box, cell, &
       grid, direction1, weight1, paw_sp, only_FO_locpot, only_FO_Qs, &
       only_SO_Qs, direction2, weight2)

    !======================================================================!
    ! This subroutine calculates the contribution from the FO (first-order)!
    ! interaction of the augmentation charge with a local potential        !
    ! \frac{d}{d\epsilon} \hat{D}_ij = \sum_LM \int \frac{d}{d\epsilon}    !
    !    [locpot(r) \hat{Q}_ij^LM (r) dr]                                  !
    ! If locpot == lhxc potential, this term is the first-order \hat{Dij}  !
    !                                                                      !
    ! d/d_eps is the atomic perturbation \epsilon, while lhxc(r) is the    !
    ! effective potential resulting from the smooth part of the density:   !
    ! v_H[\tilde{n}+\hat{n}+\tilde{n}_Zc]                                  !
    !         +v_xc[\tilde{n}+\hat{n}+\tilde{n}_c]                         !
    !----------------------------------------------------------------------!
    ! Arguments:                                                           !
    !   dij(inout)    : first-order (FO) energy term d/d_eps \hat{D}_ij.   !
    !   locpot(in)    : Local potential on fine grid.                      !
    !   FO_locpot(in) : First-order local potential on fine grid           !
    !   grid(in)      : Grid definition for fine grid.                     !
    !   only_FO_locpot (in) : if yes, compute only                         !
    !    \sum_LM\int \hat{Q}_ij^LM(r) \frac{d}{d\epsilon}                  !
    !                                     \tilde{v}_H[\tilde{rho}_Zc](r)   !
    !   only_FO_Qs (in) : if yes, compute only                             !
    !    \sum_LM\int \tilde{v}_eff(r) \frac{d}{d\epsilon}\hat{Q}_ij^{LM}(r)!
    !   only_SO_Qs (in) : if present and yes, and if cart_2 is present,    !
    !   compute only  \sum_LM\int \tilde{v}_eff(r) \frac{d^2}{d\epsilon~   !
    !   d\lambda} \hat{Q}_ij^{LM}(r)!
    !----------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in May 2015, based on              !
    !      augmentation_screen_dij                                         !
    !======================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_extract_box, &
         cell_grid_box_start_wrt_atom
    use comms, only: pub_my_proc_id
    use constants, only: max_spins
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use geometry, only: geometry_magnitude
    use parallel_strategy, only: PARAL_INFO
    use paw, only: paw_atom_aug_integrals, paw_atom_grad_aug_integrals, &
         PAW_SPECIES
    use rundat, only: pub_paw, pub_usp, pub_aug_funcs_recip, pub_debug_on_root,&
         pub_debug, pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_get_block, sparse_put_block, sparse_get_par
    use spherical_wave, only: sw_init
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_sanity_check, utils_assert

    implicit none

    ! Arguments
    type(FFTBOX_INFO),intent(in) :: aug_box
    type(CELL_INFO),intent(in) :: cell
    type(GRID_INFO),intent(in) :: grid
    type(SPAM3),intent(inout) :: dij(:)
    real(kind=DP),intent(in) :: locpot(grid%ld1,grid%ld2, grid%max_slabs12, &
         size(dij))
    real(kind=DP),intent(in) :: FO_locpot(grid%ld1,grid%ld2, grid%max_slabs12, &
         size(dij))
    integer, intent(in) :: direction1(:)
    real(kind=DP), intent(in) :: weight1(:)
    type(PAW_SPECIES), intent(in) :: paw_sp(:)

    logical, intent(in) :: only_FO_locpot
    logical, intent(in) :: only_FO_Qs
    logical, intent(in) :: only_SO_Qs
    integer, intent(in), optional :: direction2(:)
    real(kind=DP), intent(in), optional :: weight2(:)

    ! Local Variables
    integer :: loc_iat, iat
    integer :: isp
    integer :: is
    integer :: ierr
    integer :: iter1,iter2
    integer :: box_n1
    integer :: box_n2
    integer :: box_n3
    integer :: box_start1
    integer :: box_start2
    integer :: box_start3
    logical :: i_need_box
    integer :: num_spins
    integer :: lmax
    real(kind=DP), allocatable :: dij_at(:,:,:)
    real(kind=DP), allocatable :: dij_at_2(:,:,:)
    real(kind=DP), allocatable :: locpot_box(:,:,:,:)
    real(kind=DP), allocatable :: FO_locpot_box(:,:,:,:)
    real(kind=DP), allocatable :: atom_aug_func(:,:,:)
    real(kind=DP), allocatable :: atom_grad_aug_func(:,:,:,:)
    real(kind=DP), allocatable :: buffer(:,:,:)
    complex(kind=DP), allocatable :: atom_aug_func_recip(:,:,:)
    complex(kind=DP), allocatable :: atom_grad_aug_func_recip(:,:,:,:)
    complex(kind=DP), allocatable :: locpot_box_recip(:,:,:,:)
    complex(kind=DP), allocatable :: FO_locpot_box_recip(:,:,:,:)
    type(PARAL_INFO), pointer :: par

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &augmentation_FO_screen_dij'
    call utils_sanity_check(locpot(:,:,:,1),'locpot')

    ! Start Timer
    call timer_clock('augmentation_FO_screen_dij',1)

    if (size(locpot,4)/=size(dij,1)) then
       call utils_abort('Error in augmentation_FO_screen_dij: Inconsistent &
            &sizes of locpot and dij')
    end if

    if (only_FO_locpot.and.only_FO_Qs) then
       call utils_abort('Error: choose between only_FO_locpot and only_FO_Qs, &
            &if any at all')
    end if

    if (only_FO_locpot.and.only_SO_Qs) then
       call utils_abort('Error: choose between only_FO_locpot and only_SO_Qs, &
            &if any at all')
    end if

    if (only_FO_Qs.and.only_SO_Qs) then
       call utils_abort('Error: choose between only_FO_Qs and only_SO_Qs')
    end if

    if (only_SO_Qs.and. (.not.present(direction2))) then
       call utils_abort('Error: need cart2 if only_SO_Qs = TRUE')
    end if

    if ((.not.only_SO_Qs).and.present(direction2)) then
       call utils_abort('Error: need only_SO_Qs = TRUE if cart2 is present')
    end if

    ! rc2013: get parallel strategy, check arguments are compatible
    call sparse_get_par(par, dij(1))
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &augmentation_F0_screen_dij: allocated parallel strategy is &
         &incompatible with paw_sp.')
    call utils_assert(par%nat == size(direction1), 'Error in &
         &augmentation_F0_screen_dij: allocated parallel strategy is &
         &incompatible with direction1.')
    call utils_assert(par%nat == size(weight1), 'Error in &
         &augmentation_F0_screen_dij: allocated parallel strategy is &
         &incompatible with weight1.')

    num_spins = size(dij,1)
    isp = 0

    ! Find size of box
    box_n1 = aug_box%total_pt1
    box_n2 = aug_box%total_pt2
    box_n3 = aug_box%total_pt3

    ! Ensure spherical waves have been initialised to high enough lmax
    ! (must be before OMP region starts to avoid race condition)
    if (pub_paw) then
       lmax = maxval(paw_sp(:)%lmax)
       if (lmax>3) call sw_init(lmax+1,1)
    end if

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(dij_at,dij_at_2,locpot_box,FO_locpot_box,buffer, &
!$OMP      atom_aug_func_recip, atom_grad_aug_func_recip,locpot_box_recip, &
!$OMP      FO_locpot_box_recip, &
!$OMP      atom_aug_func,atom_grad_aug_func,loc_iat,iat,isp,is,ierr, &
!$OMP      box_start1,box_start2, box_start3,i_need_box) &
!$OMP SHARED (pub_my_proc_id,grid,cell,par,max_proj_tot,locpot,num_spins, &
!$OMP      aug_box,box_n1,box_n2,box_n3,pub_aug_funcs_recip, &
!$OMP      pub_threads_num_fftboxes,dij,pub_paw,pub_usp, &
!$OMP      direction1,weight1,direction2,weight2,&
!$OMP      FO_locpot,paw_sp, only_FO_locpot,only_FO_Qs,only_SO_Qs)

    ! Allocate temporary arrays
    allocate(dij_at(max_proj_tot, max_proj_tot, num_spins),stat=ierr)
    call utils_alloc_check('augmentation_FO_screen_dij','dij_at',ierr)
    allocate(dij_at_2(max_proj_tot, max_proj_tot, num_spins),stat=ierr)
    call utils_alloc_check('augmentation_FO_screen_dij','dij_at_2',ierr)
    allocate(locpot_box(box_n1,box_n2,box_n3,num_spins),stat=ierr)
    call utils_alloc_check('augmentation_FO_screen_dij','locpot_box',ierr)
    allocate(FO_locpot_box(box_n1,box_n2,box_n3,num_spins),stat=ierr)
    call utils_alloc_check('augmentation_FO_screen_dij','FO_locpot_box',ierr)
    allocate(buffer(box_n1,box_n2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('augmentation_FO_screen_dij','buffer',ierr)
    if (pub_aug_funcs_recip) then
       allocate(atom_aug_func_recip(box_n1,box_n2,box_n3),stat=ierr)
       call utils_alloc_check('augmentation_FO_screen_dij', &
            'atom_aug_func_recip',ierr)
       allocate(locpot_box_recip(box_n1,box_n2,box_n3,num_spins),stat=ierr)
       call utils_alloc_check('augmentation_FO_screen_dij', &
            'locpot_box_recip',ierr)
       allocate(atom_grad_aug_func_recip(box_n1,box_n2,box_n3,3),stat=ierr)
       call utils_alloc_check('augmentation_FO_screen_dij', &
            'atom_grad_aug_func_recip',ierr)
       allocate(FO_locpot_box_recip(box_n1,box_n2,box_n3,num_spins),stat=ierr)
       call utils_alloc_check('augmentation_FO_screen_dij', &
            'FO_locpot_box_recip', ierr)
    else
       allocate(atom_aug_func(box_n1,box_n2,box_n3),stat=ierr)
       call utils_alloc_check('augmentation_FO_screen_dij', &
            'atom_grad_aug_func', ierr)
       allocate(atom_grad_aug_func(box_n1,box_n2,box_n3,3),stat=ierr)
       call utils_alloc_check('augmentation_FO_screen_dij', &
            'atom_grad_aug_func', ierr)
    end if

    ! Loop over atoms
!$OMP DO
    do loc_iat = 1,par%max_atoms_on_proc

       ! Only need to extract if there is an atom left on this proc
       if (loc_iat<=par%num_atoms_on_proc(pub_my_proc_id)) then
          iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
          isp = par%elements_on_proc(loc_iat)%pspecies_number

          ! Find where box for this atom is located in simulation cell
          call cell_grid_box_start_wrt_atom(box_start1, box_start2, &
               box_start3, par%elements_on_proc(loc_iat)%centre, box_n1, &
               box_n2, box_n3, grid,cell)
          i_need_box = .true.
       else
          i_need_box = .false.
       end if

       ! Extract tightbox of data from effective potential over simulation
       ! cell for this atom

       locpot_box = 0.0_DP
       FO_locpot_box = 0.0_DP

!$OMP CRITICAL
       do is = 1,num_spins
          call cell_grid_extract_box(locpot_box(:,:,:,is), buffer, &
               locpot(:,:,:,is), grid, box_n1, box_n2, box_n3, box_n1, &
               box_n2, box_start1, box_start2, box_start3, i_need_box, .false.)

          call cell_grid_extract_box(FO_locpot_box(:,:,:,is), buffer, &
               FO_locpot(:,:,:,is), grid, box_n1, box_n2, box_n3, box_n1, &
               box_n2, box_start1, box_start2, box_start3, i_need_box, .false.)
       end do  ! over is
!$OMP END CRITICAL

       ! Only need to take overlap if there is an atom left on this proc
       if (loc_iat<=par%num_atoms_on_proc(pub_my_proc_id)) then
          iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
          isp = par%elements_on_proc(loc_iat)%pspecies_number

          if (pub_debug) call utils_sanity_check(locpot_box(:,:,:,1), &
               'locpot_box')
          if (pub_debug) call utils_sanity_check(FO_locpot_box(:,:,:,1), &
               'FO_locpot_box')


          dij_at(:,:,:) = 0.0_DP
          dij_at_2(:,:,:) = 0.0_DP

          if (pub_aug_funcs_recip) then
             do is = 1,num_spins
               FO_locpot_box_recip(:,:,:,is) = FO_locpot_box(:,:,:,is)
               call fourier_apply_box('F','F',FO_locpot_box_recip(:,:,:,is),&
                    box=aug_box)

               locpot_box_recip(:,:,:,is) = locpot_box(:,:,:,is)
               call fourier_apply_box('F','F',locpot_box_recip(:,:,:,is), &
                    box=aug_box)

             end do ! over is
          end if


          ! gcc32: we actually need v(-G), which for purely real potentials
          ! is the same as v(G)^*
          locpot_box_recip = conjg(locpot_box_recip) / &
               real(box_n1*box_n2*box_n3,kind=DP)
          FO_locpot_box_recip = conjg(FO_locpot_box_recip) / &
               real(box_n1*box_n2*box_n3,kind=DP)


          ! Call appropriate routine to calculate the integral of the
          ! augmentation function for this atom with the local potential
          ! in the augmentation box
          if (pub_paw) then
             if (pub_aug_funcs_recip) then
                call paw_atom_aug_integrals(dij_at, num_spins, isp, &
                     par%elements_on_proc(loc_iat)%centre, grid, cell, &
                     paw_sp, box_n1, box_n2, box_n3, box_start1, &
                     box_start2, box_start3, &
                     locpot_box_recip=FO_locpot_box_recip, &
                     atom_aug_func_recip=atom_aug_func_recip)

             else
                call paw_atom_aug_integrals(dij_at, num_spins, isp, &
                     par%elements_on_proc(loc_iat)%centre, grid, cell, &
                     paw_sp, box_n1, box_n2, box_n3, box_start1, &
                     box_start2,box_start3, FO_locpot_box, atom_aug_func)
             end if
          end if


          ! Call appropriate routine to calculate the integral of the
          ! augmentation function for this atom with the local potential
          ! in the augmentation box
          if (pub_paw) then
             if (pub_aug_funcs_recip) then
                if (only_SO_Qs) then
                   call paw_atom_grad_aug_integrals(dij_at_2,num_spins,isp, &
                        par%elements_on_proc(loc_iat)%centre, grid, cell, &
                        box_n1, box_n2, box_n3, box_start1, box_start2, &
                        box_start3, direction1(par%orig_atom(iat)), paw_sp, &
                        locpot_box_recip=locpot_box_recip, &
                        atom_grad_aug_func_recip=atom_grad_aug_func_recip, &
                        cart2 = direction2(par%orig_atom(iat)))

                    dij_at_2 = dij_at_2 * weight1(par%orig_atom(iat)) * &
                         weight2(par%orig_atom(iat))
                else
                   call paw_atom_grad_aug_integrals(dij_at_2,num_spins,isp, &
                        par%elements_on_proc(loc_iat)%centre, grid, cell, &
                        box_n1, box_n2, box_n3, box_start1, box_start2, &
                        box_start3, direction1(par%orig_atom(iat)), paw_sp, &
                        locpot_box_recip=locpot_box_recip, &
                        atom_grad_aug_func_recip=atom_grad_aug_func_recip)

                   dij_at_2 = dij_at_2 * weight1(par%orig_atom(iat))
                end if
             else
                if (only_SO_Qs) then
                   call utils_abort('Error: only_SO_Qs implemented only &
                        &for pub_aug_funcs_recip = TRUE')
                else
                   call paw_atom_grad_aug_integrals(dij_at_2,num_spins,isp, &
                        par%elements_on_proc(loc_iat)%centre, grid, cell, &
                        box_n1, box_n2, box_n3, box_start1, box_start2, &
                        box_start3, direction1(par%orig_atom(iat)), paw_sp, &
                        locpot_box, atom_grad_aug_func)
                   dij_at_2 = dij_at_2 * weight1(par%orig_atom(iat))
                end if
             end if
          end if

          ! GATHER ALL THE TERMS

          do is = 1,num_spins
             do iter1 = 1,max_proj_tot
                do iter2 = 1,max_proj_tot

                   ! The minus appears since in paw_atom_grad_aug_integrals
                   ! we have d/dr, instead of d/dR_I

                   if ((.not.only_FO_locpot) .and.  (.not.(only_FO_Qs.or. &
                        only_SO_Qs))) then
                      dij_at(iter1,iter2,is) = dij_at(iter1,iter2,is) - &
                           dij_at_2(iter1,iter2,is)
                   end if

                   if  ((.not.only_FO_locpot) .and. only_FO_Qs) then
                      dij_at(iter1,iter2,is) = -1.0_DP*dij_at_2(iter1,iter2,is)
                   end if

                   if (only_FO_locpot .and.  (.not.(only_FO_Qs.or. &
                        only_SO_Qs))) then
                    ! this is already true:
                    !  dij_at(iter1,iter2,is) = dij_at(iter1,iter2,is)
                   end if

                   ! here we do not need an extra minus, since for second-
                   ! order, the extra '-' signs cancel
                   if  ((.not.only_FO_locpot) .and. only_SO_Qs) then
                      dij_at(iter1,iter2,is) = dij_at_2(iter1,iter2,is)
                   end if


                end do ! over iter2
             end do ! over iter1
          end do ! over is

          ! Put block of dij into SPAM3 matrix
          do is=1,num_spins
             call sparse_put_block(dij_at(:,:,is), dij(is), iat, iat)
          end do  ! is

       end if  ! loc_iat < loc_nat
    end do  ! loc_iat
!$OMP END DO

    ! Deallocate temporary arrays
    if (pub_aug_funcs_recip) then
       deallocate(locpot_box_recip,stat=ierr)
       call utils_dealloc_check('augmentation_FO_screen_dij', &
            'locpot_box_recip',ierr)
       deallocate(atom_aug_func_recip,stat=ierr)
       call utils_dealloc_check('augmentation_FO_screen_dij', &
            'atom_aug_func_recip',ierr)
       deallocate(FO_locpot_box_recip,stat=ierr)
       call utils_dealloc_check('augmentation_FO_screen_dij', &
            'FO_locpot_box_recip',ierr)
       deallocate(atom_grad_aug_func_recip,stat=ierr)
       call utils_dealloc_check('augmentation_FO_screen_dij', &
            'atom_grad_aug_func_recip',ierr)

    else
       deallocate(atom_aug_func,stat=ierr)
       call utils_dealloc_check('augmentation_FO_screen_dij','atom_aug_func', &
            ierr)
       deallocate(atom_grad_aug_func,stat=ierr)
       call utils_dealloc_check('augmentation_FO_screen_dij', &
            'atom_grad_aug_func',ierr)
    end if
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('augmentation_FO_screen_dij','buffer',ierr)
    deallocate(locpot_box,stat=ierr)
    call utils_dealloc_check('augmentation_FO_screen_dij','locpot_box',ierr)
    deallocate(FO_locpot_box,stat=ierr)
    call utils_dealloc_check('augmentation_FO_screen_dij','FO_locpot_box',ierr)
    deallocate(dij_at,stat=ierr)
    call utils_dealloc_check('augmentation_FO_screen_dij','dij_at',ierr)
    deallocate(dij_at_2,stat=ierr)
    call utils_dealloc_check('augmentation_FO_screen_dij','dij_at_2',ierr)

!$OMP END PARALLEL

    ! Stop Timer
    call timer_clock('augmentation_FO_screen_dij',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving augmentation_FO_screen_dij'

  end subroutine augmentation_FO_screen_dij


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine aug_nl_calculate_forces(nlps_forces, &
       ngwfs_on_grid,ngwf_basis,proj_basis,nl_projectors, &
       sp_overlap,mdl,inv_overlap,pur_denskern,ham, dijhat, kpt) ! agrecokpt

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the nonlocal projectors.                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  nlps_forces     : output : nonlocal forces                             !
    !  sp_overlap      : input  : ngwf-projector overlap matrix               !
    !  ngwfs_on_grid   : input  : NGWF data in ppd format                     !
    !  ngwf_basis      : input  : Function basis description for NGWFs        !
    !  proj_basis      : input  : Function basis description for projectors   !
    !  pur_denskern    : input  : purified density kernel SPAM3               !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 08/06/2010.                                 !
    ! Modified to deal with embedding and to remove pub_par by Joseph         !
    ! Prentice, May 2018                                                      !
    !=========================================================================!

    use datatypes, only: FUNCTIONS
    use cell_grid, only: GRID_INFO
    use comms, only: comms_barrier, comms_reduce, pub_my_proc_id
    use function_basis, only: FUNC_BASIS
    ! agrecokpt
    use geometry, only: POINT
    use model_type, only: MODEL
    use paw, only: paw_projector_overlap, paw_nonlocal_energies
    use projectors, only: PROJECTOR_SET, projectors_func_grad_ovlp_box
    use pseudopotentials, only: pseudo_get_dij, pseudo_aug_Q_matrix
    use rundat, only: pub_num_spins, pub_paw, pub_usp, pub_debug_on_root, &
         pub_imag_thr
    use sparse, only: SPAM3, sparse_get_element, sparse_create, &
         sparse_copy, sparse_destroy, sparse_take_real_part, sparse_axpy, &
         sparse_product, sparse_transpose, sparse_transpose_structure, &
         sparse_scale
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: nlps_forces(:,:)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(PROJECTOR_SET), intent(inout) :: nl_projectors
    type(SPAM3), intent(in) :: sp_overlap
    type(SPAM3), intent(in) :: inv_overlap
    type(SPAM3), intent(in) :: pur_denskern(pub_num_spins)
    type(SPAM3), intent(in) :: ham(pub_num_spins)
    type(SPAM3), intent(in) :: dijhat(pub_num_spins)
    type(MODEL), intent(in) :: mdl
    ! agrecokpt
    type(POINT), optional, intent(in) :: kpt

    ! Local Variables
    type(SPAM3) :: siGp_overlap,iGps_overlap,ps_overlap
    type(SPAM3) :: rkq,dij_rkq,kh,khsq,oij
    type(SPAM3),allocatable :: dij(:)
    type(SPAM3),allocatable :: rho_ij(:)
    type(SPAM3),allocatable :: khs(:)
    integer :: cart, is
    integer :: iat, orig_iat, first_atom, last_atom
    integer :: atom_proj, global_proj
    integer :: ierr
    real(kind=DP) :: proj_force
    ! jme: real/complex objects and pointer to select them
    type(SPAM3),       target  :: nl_force_mat(3)
    type(SPAM3),       target  :: nl_force_mat_cmplx
    type(SPAM3),       pointer :: nl_force_mat_ptr
    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt
    type(POINT) :: loc_kpt

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering aug_nl_calculate_forces'

    ! Start timer
    call timer_clock('aug_nl_calculate_forces',1)

    ! Initialise
    nlps_forces  = 0.0_DP

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

    ! agrecokpt: specification of kpt allowed only in complex case
    if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.&
        loc_kpt%z/=0.0_DP) then
       if (.not.loc_cmplx) then
          call utils_abort('Error in aug_nl_calculate_forces: &
               &non-Gamma k-point requires complex projectors.')
       end if
    end if


    ! Create result matrices to hold nonlocal forces, projector-by-projector
    do cart=1,3
       nl_force_mat(cart)%structure = 'E'
       call sparse_create(nl_force_mat(cart))
    end do

    ! Create matrix arrays and structures
    oij%structure = 'E'
    call sparse_create(oij)
    allocate(dij(pub_num_spins),stat=ierr)
    call utils_alloc_check('aug_nl_calculate_forces','dij',ierr)
    allocate(rho_ij(pub_num_spins),stat=ierr)
    call utils_alloc_check('aug_nl_calculate_forces','rho_ij',ierr)
    do is=1,pub_num_spins
       dij(is)%structure = 'E'
       call sparse_create(dij(is))
    end do
    do is=1,pub_num_spins
       call sparse_create(rho_ij(is),dij(is))
    end do

    if (loc_cmplx) then
       ! Create auxiliary complex matrix
       call sparse_create(nl_force_mat_cmplx, nl_force_mat(1), &
            iscmplx=loc_cmplx)
    end if

    if (pub_paw) then

       ! Get projector overlap matrix
       call paw_projector_overlap(oij,mdl%regions(1)%paw_sp)

       ! Get projector density kernel
       call aug_projector_denskern(rho_ij,pur_denskern,sp_overlap)

       ! Get nonlocal energies in diagonal matrix
       call paw_nonlocal_energies(dij,rho_ij,mdl%regions(1)%paw_sp,mdl%par)
       do is=1,pub_num_spins
          call sparse_axpy(dij(is),dijhat(is),1.0_DP)
       end do

    else if (pub_usp) then

       ! Get projector overlap matrix
       call pseudo_aug_Q_matrix(oij,mdl%regions(1)%pseudo_sp)

       ! Get nonlocal energies in diagonal matrix
       do is=1,pub_num_spins
          call pseudo_get_dij(dij(is),mdl%regions(1)%pseudo_sp)
          call sparse_axpy(dij(is),dijhat(is),1.0_DP)
       end do

    end if

    do is=pub_num_spins,1,-1
       call sparse_destroy(rho_ij(is))
    end do
    deallocate(rho_ij,stat=ierr)
    call utils_dealloc_check('aug_nl_calculate_forces','rho_ij',ierr)

    ! Create matrices to hold <phi|iG*proj> overlap matrix and transpose
    call sparse_transpose_structure(ps_overlap%structure,sp_overlap)
    ! agrecocmplx
    call sparse_create(ps_overlap, iscmplx=loc_cmplx)
    call sparse_transpose(ps_overlap,sp_overlap)
    call sparse_create(siGp_overlap,sp_overlap)
    call sparse_create(iGps_overlap,ps_overlap)

    ! Create temporary matrices rkq and dij_rkq
    ! agrecocmplx
    call sparse_create(rkq,nl_force_mat(1),iscmplx=loc_cmplx)
    call sparse_create(dij_rkq,rkq)

    ! Calculate K H S^-1
    allocate(khs(pub_num_spins),stat=ierr)
    call utils_alloc_check('aug_nl_calculate_forces','khs',ierr)
    call sparse_create(kh,pur_denskern(1),ham(1))
    call sparse_destroy(kh)
    do is=1,pub_num_spins
       call sparse_create(khs(is),kh,inv_overlap)
    end do
    call sparse_create(kh)
    do is=1,pub_num_spins
       call sparse_product(kh,pur_denskern(is),ham(is))
       call sparse_product(khs(is),kh,inv_overlap)
    end do
    call sparse_destroy(kh)
    call sparse_create(khsq,khs(1),sp_overlap)

    ! Loop over Cartesian directions
    do cart=1,3

       ! jme: initialise pointers
       if (loc_cmplx) then
          nl_force_mat_ptr => nl_force_mat_cmplx
          ! agrecocmplx: reset nl_force_mat_cmplx for the
          ! new cartesian component being computed
          call sparse_scale(nl_force_mat_cmplx, 0.0_DP)
       else
          nl_force_mat_ptr => nl_force_mat(cart)
       end if

       ! Calculate <phi|iG*proj> overlap matrix
       ! agrecokpt: at specified k-point
       call projectors_func_grad_ovlp_box(siGp_overlap, &
            ngwfs_on_grid,ngwf_basis,proj_basis,nl_projectors, &
            mdl%fftbox, mdl%cell, cart, kpt=loc_kpt)

       ! Transpose it to get <iG*proj|phi> overlap matrix
       call sparse_transpose(iGps_overlap,siGp_overlap)

       do is=1,pub_num_spins

          ! Calculate D_ij <proj_j|phi_a> K^ab <phi_b|iG.proj_i>
          call sparse_product(khsq,pur_denskern(is),siGp_overlap)
          call sparse_product(rkq,ps_overlap,khsq)
          call sparse_product(dij_rkq, dij(is), rkq, &
               allow_mix_types=loc_cmplx)
          call sparse_axpy(nl_force_mat_ptr, dij_rkq, 1.0_DP)

          ! Calculate <iG.proj_i|phi_a> K^ab <phi_b|proj_j> D_ji
          call sparse_product(khsq,pur_denskern(is),sp_overlap)
          call sparse_product(rkq,iGps_overlap,khsq)
          call sparse_product(dij_rkq, rkq, dij(is), &
               allow_mix_types=loc_cmplx)
          call sparse_axpy(nl_force_mat_ptr, dij_rkq, 1.0_DP)

          ! Calculate O_ij <proj_j|phi_c>K^cd H_da S^ab <phi_b|iG.proj_i>
          call sparse_product(khsq,khs(is),siGp_overlap)
          call sparse_product(rkq,ps_overlap,khsq)
          call sparse_product(dij_rkq, oij, rkq, &
               allow_mix_types=loc_cmplx)
          call sparse_axpy(nl_force_mat_ptr, dij_rkq, -1.0_DP)

          ! Calculate <iG.proj_i|phi_a> K^ag H_gd S^da <phi_b|proj_j> O_ji
          call sparse_product(khsq,khs(is),sp_overlap)
          call sparse_product(rkq,iGps_overlap,khsq)
          call sparse_product(dij_rkq, rkq, oij, &
               allow_mix_types=loc_cmplx)
          call sparse_axpy(nl_force_mat_ptr, dij_rkq, -1.0_DP)

       end do

       if (loc_cmplx) then
          ! agrecocmplx: convert to real matrix since forces are real
          ! use safe routine to check imaginary part of complex matrix
          call sparse_take_real_part(nl_force_mat(cart), &
               nl_force_mat_cmplx, pub_imag_thr)
       end if

    end do

    ! jme: nullify pointer
    nullify(nl_force_mat_ptr)
    ! Destroy temporary matrices and deallocate arrays
    call sparse_destroy(khsq)
    do is=pub_num_spins,1,-1
       call sparse_destroy(khs(is))
    end do
    deallocate(khs,stat=ierr)
    call utils_dealloc_check('aug_nl_calculate_forces','khs',ierr)
    call sparse_destroy(dij_rkq)
    call sparse_destroy(rkq)
    call sparse_destroy(iGps_overlap)
    call sparse_destroy(siGp_overlap)
    call sparse_destroy(ps_overlap)
    call sparse_destroy(oij)
    do is=pub_num_spins,1,-1
       call sparse_destroy(dij(is))
    end do
    deallocate(dij,stat=ierr)
    call utils_dealloc_check('aug_nl_calculate_forces','dij',ierr)
    ! agrecocmplx: destroy temp complex matrices if needed
    if (loc_cmplx) call sparse_destroy(nl_force_mat_cmplx)

    ! Loop over atoms
    first_atom = mdl%par%first_atom_on_proc(pub_my_proc_id)
    last_atom  = first_atom + &
         mdl%par%num_atoms_on_proc(pub_my_proc_id) - 1
    do iat = first_atom, last_atom

       ! Find atom number in input file order
       orig_iat = mdl%par%orig_atom(iat)

       ! Loop over projectors on this atom
       do atom_proj=1,proj_basis%num_on_atom(iat)
          global_proj = proj_basis%first_on_atom(iat) + atom_proj - 1

          ! Loop over Cartesian co-ordinates
          do cart=1,3

             ! Find contribution of this projector to force on this atom
             ! from diagonal elements of nl_force_mat for this coordinate
             call sparse_get_element(proj_force,nl_force_mat(cart), &
                  global_proj, global_proj)
             nlps_forces(cart,orig_iat) = nlps_forces(cart,orig_iat) + proj_force

          end do  ! cart

       end do  ! atom_proj

    end do   ! loc_iat

    ! Reduce result across procs
    call comms_barrier
    call comms_reduce('SUM',nlps_forces,3*mdl%nat)

    ! Destroy temporary matrices
    do cart=1,3
       call sparse_destroy(nl_force_mat(cart))
    end do

    ! Stop timer
    call timer_clock('aug_nl_calculate_forces',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving aug_nl_calculate_forces'

  end subroutine aug_nl_calculate_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_density_forces(nhat_forces,denskern,sp_overlap, &
       locpot,mdl)

    !=========================================================================!
    ! This subroutine calculates the contribution to the forces on each atom  !
    ! from the interaction of the compensation density nhat with the          !
    ! effective potential.                                                    !
    !  F_I = -d/dR_I(\int v_eff(r) n_hat(r) dr)                               !
    !      = \int v_eff(r) \sum_ijLM d/dR_I(g_L(r)S_LM(r)) rho_ij q_ij^LM dr  !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30 May 2010.                                !
    ! Modified to remove pub_par and deal with embedding by Joseph Prentice,  !
    ! May 2018                                                                !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_extract_box, &
         cell_grid_box_start_wrt_atom
    use comms, only: comms_barrier, comms_reduce, pub_my_proc_id
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use geometry, only: geometry_magnitude
    use model_type, only: MODEL
    use paw, only: paw_atom_aug_force
    use pseudopotentials, only: pseudo_atom_aug_force
    use rundat, only: pub_paw, pub_usp, pub_aug_funcs_recip, &
         pub_num_spins, pub_debug_on_root
!$  use rundat, only: pub_threads_num_fftboxes
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_get_block
    use spherical_wave, only: sw_init
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    real(kind=DP),intent(out) :: nhat_forces(3,mdl%nat)
    real(kind=DP),intent(in) :: locpot(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins)
    type(SPAM3),intent(in) :: denskern(pub_num_spins)
    type(SPAM3),intent(in) :: sp_overlap

    ! Local Variables
    type(SPAM3), allocatable :: rho_ij(:)
    integer :: loc_iat, orig_iat, iat
    integer :: isp
    integer :: is
    integer :: ierr
    integer :: box_n1
    integer :: box_n2
    integer :: box_n3
    integer :: box_start1
    integer :: box_start2
    integer :: box_start3
    logical :: i_need_box
    integer :: lmax
    real(kind=DP), allocatable :: buffer(:,:,:)
    real(kind=DP), allocatable :: rhoij_at(:,:,:)
    real(kind=DP), allocatable :: locpot_box(:,:,:,:)
    real(kind=DP), allocatable :: atom_grad_aug_func(:,:,:,:)
    complex(kind=DP), allocatable :: locpot_box_recip(:,:,:,:)
    complex(kind=DP), allocatable :: atom_grad_aug_func_recip(:,:,:,:)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering augmentation_density_forces'

    ! Start Timer
    call timer_clock('augmentation_density_forces',1)

    isp = 0

    ! Find size of box
    box_n1 = mdl%aug_box%total_pt1
    box_n2 = mdl%aug_box%total_pt2
    box_n3 = mdl%aug_box%total_pt3

    ! Ensure spherical waves have been initialised to high enough lmax
    ! (must be before OMP region starts to avoid race condition)
    if (pub_paw) then
       lmax = maxval(mdl%paw_sp(:)%lmax)
       if (lmax>3) call sw_init(lmax+1,1)
    end if

    ! Create projector density kernel
    allocate(rho_ij(pub_num_spins),stat=ierr)
    call utils_alloc_check('augmentation_density_forces','rho_ij',ierr)

    do is=1,pub_num_spins
       rho_ij(is)%structure = 'E'
       call sparse_create(rho_ij(is))
    end do  ! is

    call aug_projector_denskern(rho_ij,denskern,sp_overlap)

    ! Initialise
    nhat_forces = 0.0_DP

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(rhoij_at,locpot_box,buffer,locpot_box_recip,atom_grad_aug_func, &
!$OMP      atom_grad_aug_func_recip,loc_iat,iat,isp,is,ierr,box_start1, &
!$OMP      box_start2,box_start3,orig_iat,i_need_box) &
!$OMP SHARED (pub_my_proc_id,mdl,max_proj_tot,locpot,pub_num_spins, &
!$OMP      box_n1,box_n2,box_n3,pub_aug_funcs_recip,pub_threads_num_fftboxes, &
!$OMP      rho_ij,pub_paw,pub_usp), &
!$OMP REDUCTION (+:nhat_forces)

    ! Allocate temporary arrays
    allocate(rhoij_at(max_proj_tot,max_proj_tot,pub_num_spins),stat=ierr)
    call utils_alloc_check('augmentation_density_forces','rhoij_at',ierr)
    allocate(locpot_box(box_n1,box_n2,box_n3,pub_num_spins),stat=ierr)
    call utils_alloc_check('augmentation_density_forces','locpot_box',ierr)
    allocate(buffer(box_n1,box_n2,mdl%fine_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('augmentation_density_forces','buffer',ierr)
    if (pub_aug_funcs_recip) then
       allocate(locpot_box_recip(box_n1,box_n2,box_n3,pub_num_spins),stat=ierr)
       call utils_alloc_check('augmentation_density_forces', &
            'locpot_box_recip',ierr)
       allocate(atom_grad_aug_func_recip(box_n1,box_n2,box_n3,3),stat=ierr)
       call utils_alloc_check('augmentation_density_forces', &
            'atom_grad_aug_func_recip',ierr)
    else
       allocate(atom_grad_aug_func(box_n1,box_n2,box_n3,3),stat=ierr)
       call utils_alloc_check('augmentation_density_forces', &
            'atom_grad_aug_func',ierr)
    end if

    ! Loop over atoms
!$OMP DO
    do loc_iat=1,mdl%par%max_atoms_on_proc

       ! Only need to extract if there is an atom left on this proc
       if (loc_iat<=mdl%par%num_atoms_on_proc(pub_my_proc_id)) then

          ! Find atom number in input file order and species number
          iat = mdl%par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
          isp = mdl%par%elements_on_proc(loc_iat)%pspecies_number
          orig_iat = mdl%par%orig_atom(iat)

          ! Find where box for this atom is located in simulation cell
          call cell_grid_box_start_wrt_atom( &
               box_start1, box_start2, box_start3, &
               mdl%par%elements_on_proc(loc_iat)%centre, &
               box_n1, box_n2, box_n3, mdl%fine_grid, mdl%cell)

          i_need_box = .true.
       else
          i_need_box = .false.
          orig_iat = -1  ! suppress warning
       end if ! loc_iat < loc_nat

       ! Extract tightbox of data from effective potential over simulation
       ! cell for this atom
!$OMP CRITICAL
       do is=1,pub_num_spins
          call cell_grid_extract_box(locpot_box(:,:,:,is), &
               buffer, locpot(:,:,:,is), mdl%fine_grid, &
               box_n1, box_n2, box_n3, box_n1, box_n2, &
               box_start1, box_start2, box_start3, i_need_box, .false.)
       end do  ! is
!$OMP END CRITICAL

       ! Fourier transform the locpot box if required
       if (pub_aug_funcs_recip) then
          do is=1,pub_num_spins
             locpot_box_recip(:,:,:,is) = locpot_box(:,:,:,is)
             call fourier_apply_box('F','B',locpot_box_recip(:,:,:,is), &
                  box=mdl%aug_box)
          end do
       end if

       ! Only need to calculate force if there is an atom left on this proc
       if (loc_iat<=mdl%par%num_atoms_on_proc(pub_my_proc_id)) then

          ! Get blocks of projector density matrix for this atom
          do is=1,pub_num_spins
             call sparse_get_block(rhoij_at(:,:,is),rho_ij(is),iat,iat)
          end do  ! is

          ! Call appropriate routine to calculate the contribution to the
          ! augmentation density force for this atom using the augmentation box
          if (pub_paw) then
             if (pub_aug_funcs_recip) then
                call paw_atom_aug_force(nhat_forces(:,orig_iat), &
                     rhoij_at,pub_num_spins,isp, &
                     mdl%par%elements_on_proc(loc_iat)%centre, &
                     mdl%fine_grid,mdl%cell,mdl%regions(1)%paw_sp,&
                     box_n1,box_n2,box_n3, &
                     box_start1,box_start2,box_start3, &
                     locpot_box_recip=locpot_box_recip, &
                     atom_grad_aug_func_recip=atom_grad_aug_func_recip)
             else
                call paw_atom_aug_force(nhat_forces(:,orig_iat), &
                     rhoij_at,pub_num_spins,isp, &
                     mdl%par%elements_on_proc(loc_iat)%centre, &
                     mdl%fine_grid,mdl%cell,mdl%regions(1)%paw_sp,&
                     box_n1,box_n2,box_n3, &
                     box_start1,box_start2,box_start3, &
                     locpot_box_real=locpot_box, &
                     atom_grad_aug_func_real=atom_grad_aug_func)
             end if
          else if (pub_usp) then
             call pseudo_atom_aug_force(nhat_forces(:,orig_iat),locpot_box, &
                  atom_grad_aug_func,rhoij_at,pub_num_spins,isp, &
                  mdl%par%elements_on_proc(loc_iat)%centre,&
                  mdl%fine_grid,mdl%cell,mdl%regions(1)%pseudo_sp,&
                  box_n1,box_n2,box_n3,box_start1,box_start2,box_start3)
          end if

       end if  ! loc_iat < loc_nat
    end do  ! loc_iat
!$OMP END DO

    ! Deallocate temporary arrays
    if (pub_aug_funcs_recip) then
       deallocate(atom_grad_aug_func_recip,stat=ierr)
       call utils_dealloc_check('augmentation_density_forces', &
            'atom_grad_aug_func_recip',ierr)
       deallocate(locpot_box_recip,stat=ierr)
       call utils_dealloc_check('augmentation_density_forces', &
            'locpot_box_recip',ierr)
    else
       deallocate(atom_grad_aug_func,stat=ierr)
       call utils_dealloc_check('augmentation_density_forces', &
            'atom_grad_aug_func',ierr)
    end if
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('augmentation_density_forces','buffer',ierr)
    deallocate(locpot_box,stat=ierr)
    call utils_dealloc_check('augmentation_density_forces','locpot_box',ierr)
    deallocate(rhoij_at,stat=ierr)
    call utils_dealloc_check('augmentation_density_forces','rhoij_at',ierr)

!$OMP END PARALLEL

    do is=pub_num_spins,1,-1
       call sparse_destroy(rho_ij(is))
    end do  ! is

    deallocate(rho_ij,stat=ierr)
    call utils_dealloc_check('augmentation_density_forces','rho_ij',ierr)

    ! Reduce result across procs
    call comms_barrier
    call comms_reduce('SUM',nhat_forces,3*mdl%nat)

    ! Stop Timer
    call timer_clock('augmentation_density_forces',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving augmentation_density_forces'

  end subroutine augmentation_density_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aug_FO_density_on_grid(nhat_FO_den_grad, aug_box, grid, cell, &
       denskern, sp_overlap, rho_ij_FO, direction, weight, paw_sp)

    !==================================================================!
    ! This subroutine creates the FO (first-order) compensation density!
    ! \frac{d}{d\epsilon} \hat{\rho}(r), and deposits it  on the       !
    ! simulation cell fine grid.                                       !
    !                                                                  !
    ! \frac{d}{d\epsilon} is a perturbation of atom 'atom' along       !
    ! direction 'cart'                                                 !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  grid (in) : Grid definition for the fine grid.                  !
    !  nhat_FO_den_grad (inout) : The compensation density \hat{n}(r)  !
    !  on the grid, including its gradient in each cartesian direction.!
    !  denskern (in) : The NGWF density kernel K^ab.                   !
    !  sp_overlap (in) : The NGWF-Projector overlap matrix <phi_a|p_i>.!
    !  rho_ij_FO(in) : first-order (FO) projector density-kernel       !
    !  atom : atom considered                                          !
    !  cart : cartesian direction considered (1, 2 or 3)               !
    !------------------------------------------------------------------!

    use cell_grid, only: GRID_INFO, cell_grid_deposit_box, &
         cell_grid_box_start_wrt_atom
    use comms, only: pub_on_root, pub_my_proc_id
    use constants, only: max_spins, VERBOSE
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: PARAL_INFO
    use paw, only: paw_atom_aug_den, PAW_SPECIES
    use rundat, only: pub_paw, pub_output_detail, pub_aug_den_dim, &
         pub_aug_funcs_recip, pub_debug_on_root, pub_num_spins
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_get_block, &
         sparse_get_par
    use spherical_wave, only: sw_init
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_assert

    implicit none

    ! Arguments
    type(FFTBOX_INFO),intent(in) :: aug_box
    type(CELL_INFO), intent(in) :: cell
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: nhat_FO_den_grad(grid%ld1, &
         grid%ld2, grid%max_slabs12, pub_num_spins, 0:pub_aug_den_dim)
    type(SPAM3), intent(in) :: denskern(pub_num_spins)
    type(SPAM3), intent(in) :: sp_overlap
    type(SPAM3), intent(in) :: rho_ij_FO(pub_num_spins)
    integer, intent(in) :: direction(:)
    real(kind=DP), intent(in) :: weight(:)
    type(PAW_SPECIES), intent(in) :: paw_sp(:)

    ! Local Variables
    type(SPAM3), allocatable :: rho_ij(:)
    integer :: iat, loc_iat
    integer :: ierr
    integer :: is, i1,i2,i3, idimen
    integer :: isp
    integer :: box_n1,box_n2,box_n3
    integer :: box_start1,box_start2,box_start3
    integer :: lmax
    real(kind=DP),allocatable :: rho_ij_block(:,:,:)
    real(kind=DP),allocatable :: rho_ij_FO_block(:,:,:)
    real(kind=DP), allocatable :: atom_nhat(:,:,:,:,:)
    real(kind=DP), allocatable :: atom_grad_nhat(:,:,:,:,:,:)
    real(kind=DP), allocatable :: atom_aug_func(:,:,:,:)
    real(kind=DP), allocatable :: atom_aug_func_grad(:,:,:,:,:)
    real(kind=DP), allocatable :: buffer(:,:,:)
    complex(kind=DP), allocatable :: atom_nhat_recip(:,:,:,:,:)
    complex(kind=DP), allocatable :: atom_grad_nhat_recip(:,:,:,:,:,:)
    complex(kind=DP), allocatable :: atom_aug_func_recip(:,:,:,:)
    complex(kind=DP), allocatable :: atom_aug_func_grad_recip(:,:,:,:,:)
    logical :: i_have_box

    real(kind=DP) :: total_nhat(max_spins,4)
    real(kind=DP) :: total_nhat_targ(max_spins,4)
    type(PARAL_INFO), pointer :: par

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &aug_FO_density_on_grid'

    ! Start Timer
    call timer_clock('aug_FO_density_on_grid',1)

    ! rc2013: get parallel strategy, check arguments are compatible
    call sparse_get_par(par, sp_overlap)
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &augmentation_F0_screen_dij: allocated parallel strategy is &
         &incompatible with paw_sp.')
    call utils_assert(par%nat == size(direction), 'Error in &
         &augmentation_F0_screen_dij: allocated parallel strategy is &
         &incompatible with direction.')
    call utils_assert(par%nat == size(weight), 'Error in &
         &augmentation_F0_screen_dij: allocated parallel strategy is &
         &incompatible with weight.')

    ! Find size of box
    box_n1 = aug_box%total_pt1
    box_n2 = aug_box%total_pt2
    box_n3 = aug_box%total_pt3

    ! Ensure spherical waves have been initialised to high enough lmax
    ! (must be before OMP region starts to avoid race condition)
    if (pub_paw) then
       lmax = maxval(paw_sp(:)%lmax)
       if (lmax>3) call sw_init(lmax+1,1)
    end if

    ! Create projector density kernel
    allocate(rho_ij(pub_num_spins),stat=ierr)
    call utils_alloc_check('aug_FO_density_on_grid','rho_ij',ierr)
    do is = 1,pub_num_spins
       rho_ij(is)%structure = 'E'
       call sparse_create(rho_ij(is))
    end do

    call aug_projector_denskern(rho_ij, denskern, sp_overlap)

    ! Allocate workspace
    allocate(rho_ij_block(max_proj_tot,max_proj_tot,pub_num_spins),stat=ierr)
    call utils_alloc_check('aug_FO_density_on_grid','rho_ij_block',ierr)
    allocate(rho_ij_FO_block(max_proj_tot,max_proj_tot,pub_num_spins),stat=ierr)
    call utils_alloc_check('aug_FO_density_on_grid','rho_ij_FO_block',ierr)
    allocate(atom_nhat(box_n1,box_n2,box_n3,pub_num_spins,4),stat=ierr)
    call utils_alloc_check('aug_FO_density_on_grid','atom_nhat',ierr)
    allocate(atom_grad_nhat(box_n1,box_n2,box_n3,pub_num_spins,3,4),stat=ierr)
    call utils_alloc_check('aug_FO_density_on_grid','atom_grad_nhat',ierr)
    allocate(buffer(box_n1,box_n2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('aug_FO_density_on_grid','buffer',ierr)

    if (pub_aug_funcs_recip) then
       allocate(atom_aug_func_recip(box_n1,box_n2,box_n3,4),stat=ierr)
       call utils_alloc_check('aug_FO_density_on_grid', &
            'atom_aug_func_recip',ierr)
       allocate(atom_aug_func_grad_recip(box_n1,box_n2,box_n3,3,4),stat=ierr)
       call utils_alloc_check('aug_FO_density_on_grid', &
            'atom_aug_func_grad_recip',ierr)
       allocate(atom_nhat_recip(box_n1,box_n2,box_n3,pub_num_spins,4),stat=ierr)
       call utils_alloc_check('aug_FO_density_on_grid', &
            'atom_nhat_recip',ierr)
       allocate(atom_grad_nhat_recip(box_n1,box_n2,box_n3,pub_num_spins, &
            3,4),stat=ierr)
       call utils_alloc_check('aug_FO_density_on_grid', &
            'atom_grad_nhat_recip',ierr)
    else
       call utils_abort('Cannot calculate gradient of first-order &
            &augmentation charge in real-space')
    end if

    total_nhat = 0.0_DP
    total_nhat_targ = 0.0_DP

    do loc_iat=1,par%max_atoms_on_proc

       if (loc_iat<=par%num_atoms_on_proc(pub_my_proc_id)) then
          iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
          isp = par%elements_on_proc(loc_iat)%pspecies_number

          atom_nhat = 0.0_DP
          atom_grad_nhat = 0.0_DP

          if (pub_aug_funcs_recip) then
             atom_nhat_recip = 0.0_DP
             atom_grad_nhat_recip = 0.0_DP
          end if


          ! Find where box for this atom is located in simulation cell
          call cell_grid_box_start_wrt_atom(box_start1, box_start2, &
               box_start3, par%elements_on_proc(loc_iat)%centre, box_n1, &
               box_n2, box_n3, grid, cell)

          ! Get block of \rho_ij for this atom
          do is = 1,pub_num_spins
             call sparse_get_block(rho_ij_block(:,:,is), rho_ij(is), iat,iat)
          end do

          ! Get block of \rho_ij_FO for this atom
          rho_ij_FO_block = 0.0_DP
          do is = 1,pub_num_spins
             call sparse_get_block(rho_ij_FO_block(:,:,is), rho_ij_FO(is), &
                  iat, iat)
          end do

          if (pub_aug_funcs_recip) then
             call paw_atom_aug_den(atom_nhat(:,:,:,:,1), &
                  atom_grad_nhat(:,:,:,:,:,1), &
                  total_nhat(:,1), total_nhat_targ(:,1), rho_ij_FO_block,&
                  isp, par%elements_on_proc(loc_iat)%centre, grid, cell, &
                  paw_sp, aug_box, box_start1, box_start2,&
                  box_start3, &
                  atom_aug_func_recip = atom_aug_func_recip(:,:,:,1), &
                  atom_aug_func_grad_recip= &
                       atom_aug_func_grad_recip(:,:,:,:,1), &
                  atom_nhat_recip=atom_nhat_recip(:,:,:,:,1), &
                  atom_grad_nhat_recip=atom_grad_nhat_recip(:,:,:,:,:,1))
           end if


          ! Call appropriate routine to generate the aug density for
          ! this atom in the aug box using rhoij
          if (pub_paw) then
             if (pub_aug_funcs_recip) then

                if (pub_aug_den_dim > 0) then

                   call paw_atom_aug_den(atom_nhat(:,:,:,:,2), &
                        atom_grad_nhat(:,:,:,:,:,2), total_nhat(:,2), &
                        total_nhat_targ(:,2), rho_ij_block, isp, &
                        par%elements_on_proc(loc_iat)%centre, grid, cell, &
                        paw_sp, aug_box, box_start1, box_start2, box_start3, &
                        atom_aug_func_recip = atom_aug_func_recip(:,:,:,2), &
                        atom_aug_func_grad_recip = &
                             atom_aug_func_grad_recip(:,:,:,:,2), &
                        atom_nhat_recip = atom_nhat_recip(:,:,:,:,2), &
                        atom_grad_nhat_recip=atom_grad_nhat_recip(:,:,:,:,:,2))

                   call paw_atom_aug_den(atom_nhat(:,:,:,:,3), &
                        atom_grad_nhat(:,:,:,:,:,3), total_nhat(:,3), &
                        total_nhat_targ(:,3), rho_ij_block, isp, &
                        par%elements_on_proc(loc_iat)%centre, grid, cell, &
                        paw_sp, aug_box, box_start1, box_start2, box_start3, &
                        atom_aug_func_recip=atom_aug_func_recip(:,:,:,3), &
                        atom_aug_func_grad_recip= &
                             atom_aug_func_grad_recip(:,:,:,:,3), &
                        atom_nhat_recip=atom_nhat_recip(:,:,:,:,3), &
                        atom_grad_nhat_recip= &
                             atom_grad_nhat_recip(:,:,:,:,:,3), &
                        grad_cartgrad = direction(par%orig_atom(iat)))
                else
                   call paw_atom_aug_den(atom_nhat(:,:,:,:,2), &
                        atom_grad_nhat(:,:,:,:,:,2), total_nhat(:,2), &
                        total_nhat_targ(:,2), rho_ij_block, isp, &
                        par%elements_on_proc(loc_iat)%centre, grid, cell, &
                        paw_sp, aug_box, box_start1, box_start2, box_start3, &
                        atom_aug_func_recip = atom_aug_func_recip(:,:,:,2), &
                        atom_aug_func_grad_recip = &
                             atom_aug_func_grad_recip(:,:,:,:,2), &
                        atom_nhat_recip = atom_nhat_recip(:,:,:,:,2), &
                        atom_grad_nhat_recip = &
                             atom_grad_nhat_recip(:,:,:,:,:,2), &
                        input_cart = direction(par%orig_atom(iat)))
                end if


             else
                call utils_abort('Cannot calculate gradient of first-order &
                     &augmentation charge in real-space')
             end if
          end if

          if (pub_aug_den_dim > 0) then

             ! copy atom_grad_nhat(:,:,:,is,cart,2) to atom_nhat(:,:,:,is,1)
             ! to get the full FO augmentation charge for atom 'atom'
             ! perturbed along direction 'cart'

             do is = 1,pub_num_spins
                atom_nhat(1:box_n1,1:box_n2,1:box_n3,is,1) = &
                     atom_nhat(1:box_n1,1:box_n2,1:box_n3,is,1) - &
                     weight(par%orig_atom(iat)) * &
                     atom_grad_nhat(1:box_n1,1:box_n2,1:box_n3,is, &
                          direction(par%orig_atom(iat)),2)
                !gcc32: the - appears because d/d\epsilon is -d/dr
             end do

             !copy atom_grad_nhat(:,:,:,is,:,3) to atom_grad_nhat(:,:,:,is,:,1)
             ! to get the gradient of the first-irder augmentation charge for
             ! this atom perturbed along a certain direction with a certain
             ! weight

             do idimen = 1,3
                do is = 1,pub_num_spins
                   atom_grad_nhat(1:box_n1,1:box_n2,1:box_n3,is,idimen,1) =  &
                        atom_grad_nhat(1:box_n1,1:box_n2,1:box_n3,is,idimen,1)-&
                        weight(par%orig_atom(iat)) * &
                        atom_grad_nhat(1:box_n1,1:box_n2,1:box_n3,is,idimen,3)
                   !gcc32: the - appears because
                   ! d/dr d/d_lambda [..] = - d/dr (d/dr)_lambda [..]
                end do
             end do
          else ! if pub_aug_den_dim == 0, i.e. LDA

             do is = 1,pub_num_spins
                atom_nhat(1:box_n1,1:box_n2,1:box_n3,is,1) = &
                     atom_nhat(1:box_n1,1:box_n2,1:box_n3,is,1) - &
                     atom_nhat(1:box_n1,1:box_n2,1:box_n3,is,2) * &
                     weight(par%orig_atom(iat))
                !gcc32: the - appears because d/d\epsilon is -d/dr
             end do

          end if  ! if pub_aug_den_dim > 0

          i_have_box = .true.
       else
          ! Nothing to deposit on this proc
          i_have_box = .false.
       end if

       ! Deposit this box to the simulation cell if present, or just wait for
       ! data from other procs if no box
       do is = 1,pub_num_spins
          call cell_grid_deposit_box(nhat_FO_den_grad(:,:,:,is,0), &
               atom_nhat(:,:,:,is,1), buffer, grid, box_n1, box_n2, box_n3, &
               box_n1, box_n2, box_start1, box_start2, box_start3, &
               i_have_box, .false.)
       end do

       ! now for the gradient of the first-order augmentation charge
       if (pub_aug_den_dim > 0) then

          do idimen = 1,3
             do is = 1,pub_num_spins
                call cell_grid_deposit_box(nhat_FO_den_grad(:,:,:,is,idimen), &
                     atom_grad_nhat(:,:,:,is,idimen,1), buffer, grid, &
                     box_n1, box_n2, box_n3, box_n1, box_n2, box_start1, &
                     box_start2, box_start3, i_have_box, .false.)
             end do
          end do

       end if

    end do ! over atoms

    ! Deallocate temporary arrays and matrices
    if (pub_aug_funcs_recip) then
       deallocate(atom_grad_nhat_recip,stat=ierr)
       call utils_dealloc_check('aug_FO_density_on_grid', &
            'atom_grad_nhat_recip',ierr)
       deallocate(atom_nhat_recip,stat=ierr)
       call utils_dealloc_check('aug_FO_density_on_grid', &
            'atom_nhat_recip',ierr)
       deallocate(atom_aug_func_grad_recip,stat=ierr)
       call utils_dealloc_check('aug_FO_density_on_grid', &
            'atom_aug_func_grad_recip',ierr)
       deallocate(atom_aug_func_recip,stat=ierr)
       call utils_dealloc_check('aug_FO_density_on_grid', &
            'atom_aug_func_recip',ierr)
    end if

    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('aug_FO_density_on_grid','buffer',ierr)
    deallocate(atom_grad_nhat,stat=ierr)
    call utils_dealloc_check('aug_FO_density_on_grid','atom_grad_nhat',ierr)
    deallocate(atom_nhat,stat=ierr)
    call utils_dealloc_check('aug_FO_density_on_grid','atom_nhat',ierr)
    deallocate(rho_ij_block,stat=ierr)
    call utils_dealloc_check('aug_FO_density_on_grid','rho_ij_block',ierr)
    deallocate(rho_ij_FO_block,stat=ierr)
    call utils_dealloc_check('aug_FO_density_on_grid','rho_ij_FO_block',ierr)


    do is = 1,pub_num_spins
       call sparse_destroy(rho_ij(is))
    end do

    deallocate(rho_ij,stat=ierr)
    call utils_dealloc_check('aug_FO_density_on_grid','rho_ij',ierr)

    ! Stop Timer
    call timer_clock('aug_FO_density_on_grid',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving aug_FO_density_on_grid'

  end subroutine aug_FO_density_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module augmentation
