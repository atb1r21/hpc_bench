! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!          S C F   F O R   M O L E C U L A R   I N T E R A C T I O N S        !
!                                M O D U L E                                  !
!=============================================================================!
!                                                                             !
! This module performs the SCF-MI algorithm for molecular interactions of     !
! Stoll, Wagenblast and Preuss:                                               !
!  H. Stoll, G. Wagenblast, and H. Preuss, Theor. Chem. Acc. 57, 169 (1980).  !
!                                                                             !
! This involves projection of the Hamiltonian matrix to localise the MOs to   !
! specific fragments: this state therefore avoids interfragmental charge      !
! delocalisation and produces charge transfer-restricted states.              !
!                                                                             !
! NOTE: This module is currently implemented as a state machine:              !
! - The procedures of this module depend of the state of                      !
!   the private object denskern_R. This object must be updated at each        !
!   SCF iteration.                                                            !
! - denskern_R is the proper (idempotent) density kernel constructed from the !
!   (non-idempotent) density kernels of the fragments in the full overlap     !
!   regime.  This quantity is constructed from the (interfragmentally         !
!   non-orthogonal) MOs obtained from the non-idempotent kernel as,           !
!                  R  =  M ( M^\dagger S M)^-1 M^\dagger                      !
!   where M are the MO coefficients obtained from the kernel, and S is the    !
!   (full) NGWF overlap matrix.                                               !
! - The denskern_R kernel is updated by calling (with the current fragment    !
!   density kernel as an argument) scfmi_construct_nonorth_kernel() .         !
!                                                                             !
!-----------------------------------------------------------------------------!
! Originally written by Max Phipps March 2015                                 !
! Modified to convert pub_R_spam3 (SPAM3_ARRAY)                               !
!      to denskern_R (DKERN) December 2015.                                   !
! Modified to allow diagonalisation-free EDA June 2016.                       !
!-----------------------------------------------------------------------------!

module fragment_scfmi

  use constants, only: DP
  use dense, only: DEM
  use kernel, only: DKERN
  use rundat, only: PUB_1K, pub_eda_nodiag, pub_num_spins
  use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY


  implicit none


  public :: scfmi_allocate_matrices
  public :: scfmi_deallocate_matrices

  public :: scfmi_hamiltonian_projection
  public :: scfmi_construct_nonorth_kernel
  public :: scfmi_construct_nonorth_kernel_diag   ! with diagonalisation
  public :: scfmi_construct_nonorth_kernel_nodiag ! diagonalisation-free

  public :: scfmi_ngwf_gradient


  ! orthogonal representation of the density kernel
  type(DKERN), public, save :: denskern_R

  ! full overlap matrix, dense format
  type(DEM), public         :: overlap_f_dens


  ! The MO overlap matrix and its inverse:
  type(DEM), allocatable, dimension(:), private :: MtSM
  type(DEM), allocatable, dimension(:), private :: inv_MtSM

  ! Intermediate components used for projecting:
  type(DEM), allocatable, dimension(:), private :: M_occ  ! occupied-only MOs
  type(DEM), allocatable, dimension(:), private :: Mt, MtS !, SMiz !, SM


  ! whether the module-global arrays have been allocated
  logical, private :: arrays_allocated

contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine scfmi_allocate_matrices(denskern,rep)

    !======================================================================!
    ! This subroutine is used to allocate the objects used in the SCF-MI   !
    ! module.                                                              !
    !----------------------------------------------------------------------!
    ! Written by Max Phipps, Mid 2015                                      !
    ! Modified for embedding by Joseph Prentice, September 2018            !
    !======================================================================!

    use constants, only: stdout
    use dense, only: dense_create
    use fragment_matrix_ops, only: fmo_allocate_masks
    use kernel, only: kernel_create
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_debug_on_root
    use sparse_embed, only: sparse_embed_num_rows
    use utils, only: utils_alloc_check

    ! Arguments
    type(DKERN), intent(in)       :: denskern
    type(NGWF_REP), intent(in)    :: rep

    ! Local Variables
    integer            :: is, ierr, num, num_occMO

    if (pub_debug_on_root) write(stdout,'(/a)') 'DEBUG: Entering &
        &scfmi_allocate_matrices'

    if (.not. arrays_allocated) then

       ! mjsp: Allocate fragment masks
       call fmo_allocate_masks()

       ! mjsp: Representation of density kernel with full overlap
       call kernel_create(denskern_R, denskern%kern%m(1,PUB_1K)%structure, .false.)

       ! number of NGWFs
       num = sparse_embed_num_rows(denskern_R%kern%m(1,PUB_1K))

       ! dense NGWF overlap storage:
       call dense_create(overlap_f_dens,num,num)

       ! mjsp: Allocate MO coefficient dependent matrices
       if (pub_eda_nodiag .ne. 3) then

          allocate(M_occ(1:pub_num_spins), stat=ierr)
          call utils_alloc_check('scfmi_allocate_matrices', 'M_occ', ierr)
          allocate(Mt(1:pub_num_spins), stat=ierr)
          call utils_alloc_check('scfmi_allocate_matrices', 'Mt', ierr)
          allocate(MtS(1:pub_num_spins), stat=ierr)
          call utils_alloc_check('scfmi_allocate_matrices', 'MtS', ierr)
          allocate(MtSM(1:pub_num_spins), stat=ierr)
          call utils_alloc_check('scfmi_allocate_matrices', 'MtSM', ierr)
          allocate(inv_MtSM(1:pub_num_spins), stat=ierr)
          call utils_alloc_check('scfmi_allocate_matrices', 'inv_MtSM', ierr)

          do is=1,pub_num_spins

             ! number of occupied MOs (spin dependent)
             num_occMO = rep%n_occ(is,PUB_1K)

             ! Eigenvector storage:
             ! Isolated fragment eigenvectors:
             call dense_create(M_occ(is), num, num_occMO)

             ! (nonorthogonal MO) density matrix components:
             call dense_create(Mt(is),num_occMO,num)
             call dense_create(MtS(is),num_occMO,num)

             ! Nonorthogonal MO metrics:
             ! mjsp: MO overlap matrix (sigma)
             call dense_create(MtSM(is),num_occMO,num_occMO)
             ! mjsp: Inverse of MO overlap matrix (sigma^-1):
             call dense_create(inv_MtSM(is),num_occMO,num_occMO)

          end do

       end if

       arrays_allocated = .true.

    end if

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
        &scfmi_allocate_matrices'

  end subroutine scfmi_allocate_matrices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine scfmi_deallocate_matrices()

    !======================================================================!
    ! This subroutine is used to deallocate the objects used in the SCF-MI !
    ! module.                                                              !
    !----------------------------------------------------------------------!
    ! Written by Max Phipps, Mid 2015                                      !
    !======================================================================!

    use dense, only: dense_destroy
    use fragment_matrix_ops, only: fmo_destroy_masks
    use kernel, only: kernel_destroy
    use utils, only: utils_dealloc_check

    ! Local Variables
    integer :: is, ierr

    if (arrays_allocated) then

       ! mjsp: Deallocate MO coefficient dependent matrices
       if (pub_eda_nodiag .ne. 3) then

          do is=1,pub_num_spins

             ! Eigenvector storage:
             ! Isolated fragment eigenvectors:
             call dense_destroy(M_occ(is))

             ! (nonorthogonal MO) density matrix components:
             call dense_destroy(Mt(is))
             call dense_destroy(MtS(is))

             ! Nonorthogonal MO metrics:
             ! mjsp: MO overlap matrix (sigma)
             call dense_destroy(MtSM(is))
             ! mjsp: Inverse of MO overlap matrix (sigma^-1):
             call dense_destroy(inv_MtSM(is))

          end do

          deallocate(M_occ, stat=ierr)
          call utils_dealloc_check('scfmi_deallocate_matrices', 'M_occ', ierr)
          deallocate(Mt, stat=ierr)
          call utils_dealloc_check('scfmi_deallocate_matrices', 'Mt', ierr)
          deallocate(MtS, stat=ierr)
          call utils_dealloc_check('scfmi_deallocate_matrices', 'MtS', ierr)
          deallocate(MtSM, stat=ierr)
          call utils_dealloc_check('scfmi_deallocate_matrices', 'MtSM', ierr)
          deallocate(inv_MtSM, stat=ierr)
          call utils_dealloc_check('scfmi_deallocate_matrices', 'inv_MtSM', ierr)

       end if

       ! mjsp: Deallocate fragment masks
       call fmo_destroy_masks()

       ! mjsp: Representation of density kernel with full overlap
       call kernel_destroy(denskern_R)

       ! dense NGWF overlap storage:
       call dense_destroy(overlap_f_dens)

       arrays_allocated = .false.

    end if

  end subroutine scfmi_deallocate_matrices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine scfmi_hamiltonian_projection(ham, rep)

    !======================================================================!
    ! This subroutine acts as a driver for calling the SCF-MI projection   !
    ! equations in order to locally-project the Hamiltonian, ham, given a  !
    ! density kernel (passed by calling scfmi_construct_nonorth_kernel()). !
    !                                                                      !
    ! Currently the only implemented equations are those of                !
    ! Stoll et.al.[1], however this driver has been included to allow for  !
    ! easy extention of the code for the alternative SCF-MI approaches,    !
    ! i.e. the equations of Gianinetti et.al.[2] or Nagata et.al.[3].      !
    !                                                                      !
    ! [1] H. Stoll, G. Wagenblast, and H. Preuss, Theor. Chem. Acc. 57,    !
    !   169 (1980).                                                        !
    ! [2] E. Gianinetti, M. Raimondi, E. Tornaghi, Int. J. Quantum Chem.   !
    !   60, 157-166 (1996)                                                 !
    ! [3] T. Nagata, O. Takahashi, K. Saito, and S. Iwata, J. Chem. Phys.  !
    !   115, 3553 2001.
    !                                                                      !
    !----------------------------------------------------------------------!
    ! WARNING: The current implementation of this module requires          !
    ! denskern_R to have been updated via scfmi_construct_nonorth_kernel() !
    ! before calling this routine!                                         !
    !----------------------------------------------------------------------!
    ! Written by Max Phipps, April 2015                                    !
    ! Modified for embedding by Joseph Prentice, September 2018            !
    !======================================================================!

    use constants, only: stdout
    use rundat, only: pub_debug_on_root, pub_eda_nodiag
    use ngwf_representation, only: NGWF_HAM, NGWF_REP
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(NGWF_HAM),   intent(inout) :: ham
    type(NGWF_REP),   intent(in   ) :: rep

    ! Start Timer
    call timer_clock('fragment_scfmi_projection',1)

    ! --------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
        &scfmi_hamiltonian_projection'

    ! locally-project the full Hamiltonian:

    ! mjsp: If all goes wrong, then diagonalisation approach
    ! mjsp: may be used as per below:
    if (pub_eda_nodiag .ge. 2) then
       call scfmi_hamiltonian_proj_stoll_nodiag(rep, ham)
    else
       call scfmi_hamiltonian_proj_stoll(rep, ham)
    end if

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
        &scfmi_hamiltonian_projection'

    ! Stop Timer
    call timer_clock('fragment_scfmi_projection',2)


  end subroutine scfmi_hamiltonian_projection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine scfmi_construct_nonorth_kernel(denskern, rep)

    !======================================================================!
    ! This subroutine updates the state of the density kernel denskern_R   !
    ! from the nonorthogonal orbitals described by the denskern passed to  !
    ! this method.                                                         !
    ! This may be achieved via a kernel diagonalisation approach or        !
    ! via purification (diagonalisation-free).                             !
    !----------------------------------------------------------------------!
    ! WARNING: The current implementation requires denskern_R to have      !
    ! been updated via a call to this method before calling the other      !
    ! state-dependent methods of this module!                              !
    !----------------------------------------------------------------------!
    ! Written by Max Phipps, April 2016                                    !
    ! Modified for embedding by Joseph Prentice, September 2018            !
    !======================================================================!

    use constants, only: EDA_FROZENIDEM
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_eda_mode
    use fragment_matrix_ops, only: fmo_freeze_mat_spam_nxn

    ! Arguments
    type(DKERN),intent(inout) :: denskern
    type(NGWF_REP), intent(in)   :: rep

    ! Local Variables
    integer :: is

    ! Ensure the kernel is block diagonal
    do is=1,pub_num_spins
       call fmo_freeze_mat_spam_nxn(denskern%kern%m(is,PUB_1K),0)
    end do

    ! Construct denskern_R from the fragment denskern
    if (((pub_eda_mode .eq. EDA_FROZENIDEM).and. & ! frozen density component
         (pub_eda_nodiag .eq. 1 .or. pub_eda_nodiag .eq. 3)) .or. &
         ((pub_eda_mode .ne. EDA_FROZENIDEM).and. & ! polarisation
         (pub_eda_nodiag .eq. 2 .or. pub_eda_nodiag .eq. 3 ))) then
       call scfmi_construct_nonorth_kernel_nodiag(denskern, rep)
    else
       call scfmi_construct_nonorth_kernel_diag(denskern%kern, rep)
    end if

  end subroutine scfmi_construct_nonorth_kernel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine scfmi_construct_nonorth_kernel_diag(denskern, rep)

    !======================================================================!
    ! This subroutine updates the state of the density kernel denskern_R   !
    ! from the nonorthogonal orbitals described by the denskern passed to  !
    ! this method.  The construction is given by,                          !
    !   R  =  M ( M^\dagger S M)^-1 M^\dagger                              !
    ! where R is the object denskern_R, M are the MO expansion             !
    ! coefficients obtained from denskern, and S is the (full) NGWF        !
    ! overlap matrix.                                                      !
    !                                                                      !
    ! This routine has a diagonalisation bottleneck: a purification        !
    ! approach to constructing denskern_R is coded in the routine          !
    ! scfmi_construct_nonorth_kernel_nodiag                                !
    !----------------------------------------------------------------------!
    ! WARNING: The current implementation requires denskern_R to have      !
    ! been updated via a call to this method before calling the other      !
    ! state-dependent methods of this module!                              !
    !----------------------------------------------------------------------!
    ! Written by Max Phipps, April 2015                                    !
    ! Improved optimisation by Max Phipps, December 2015                   !
    ! Modified for embedding by Joseph Prentice, September 2018            !
    !======================================================================!

    use constants, only: stdout, EDA_CTFRAGLOC_DEVEL
    use dense, only: DEM, dense_create, dense_convert, dense_destroy, &
         dense_eigensolve, dense_product, dense_scale, dense_copy, &
         dense_transpose, &
         dense_get_col, dense_put_col, &
         dense_invert
    use fragment_matrix_ops, only: fmo_get_frag_blk_dens_nxn
    use fragment_data, only: NGWF_INDEX_MAP_TYPE, pub_frag_data, &
         fragment_data_get_ngwf_index_map
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_frag_iatm, pub_eda_mode, &
         pub_frag_counter, pub_frag_counter2, &
         pub_debug_on_root, &
         pub_eigensolver_orfac, pub_eigensolver_abstol
    use sparse_embed, only: sparse_embed_num_rows
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    ! Arguments
    type(SPAM3_EMBED_ARRAY),intent(in) :: denskern
    type(NGWF_REP), intent(in)   :: rep

    ! Local Variables
    integer :: num        ! number of NGWFs
    integer :: num_occMO  ! number of occupied MOs
    integer :: num_frag_ngwf   ! number of MOs on a frag
    integer :: num_frag_occMO  ! number of occupied MOs on a frag
    integer :: ierr
    integer :: ngwf_offset, ingwf_frag, mo_offset ! for MO coefficient extraction
    integer :: ingwf_in_super
    integer :: is                              ! spin loop counter
    integer :: it, it2, it_i                   ! fragment loop counters
    integer :: join_skip_mo, join_skip_ngwf    ! skip for if combining fragments
    real(kind=DP), dimension(:), allocatable :: eigenvalues
    ! for fragment to supermolecule vector redistribution:
    real(kind=DP)              :: col(pub_frag_data(0)%rep%n_occ(1,PUB_1K))
    real(kind=DP), allocatable :: colbuf(:)

    type(DEM) :: frag_dkn_dens ! fragment density kernels
    type(DEM) :: frag_eigs_dens   ! The fragment MO eigenvectors (contravariant)
    type(DEM) :: frag_overlap_dens ! inverse of fragment overlap blocks
    type(DEM) :: frag_eigs_dens_t ! The fragment MO eigenvectors (contravariant)
                                                               ! (transposed)

    type(DEM) :: dkn_dens ! full (block diagonal) density kernel
    type(DEM) :: Miz ! temporary intermediate matrix
    type(DEM) :: R ! the proper representation of the density kernel in the full
                   ! overlap

    logical :: jf_log ! join_frags logical (for fragment pair delocalisations)

    type(NGWF_INDEX_MAP_TYPE) :: ngwf_index_map

    ! Start Timer
    call timer_clock('fragment_scfmi_construct_kern_diag',1)

    ! --------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
        &scfmi_construct_nonorth_kernel_diag'

    ngwf_index_map = fragment_data_get_ngwf_index_map()

    ! number of NGWFs
    num = sparse_embed_num_rows(rep%overlap)

    ! mjsp: DEM buffers:
    ! Input kernel's dense storage:
    call dense_create(dkn_dens,num,num)
    ! Representation of non-orth kernel:
    call dense_create(R,num,num)


    ! Overview:
    ! 1. Calculate the fragment localised MO eigenvectors
    ! (by diagonalisation using the fragment-recast kernels and the
    ! fragment-cast overlap matrix)
    ! 2. Construct KS state overlap and inverse overlap metrics from KS
    ! eigenvectors (using the true full NGWF overlap matrix).
    !
    ! denskern_R now prepared for projecting the Hamiltonian in a fragment-wise
    ! fashion.

    ! mjsp: loop over spins
    do is=1,pub_num_spins

       ! mjsp: expand kernel to dense matrix
       call dense_convert(dkn_dens,denskern%m(is,PUB_1K))

       ! =========== 1. OBTAIN FRAGMENT-LOCALISED MO COEFFICIENTS  ========

       ! mjsp: Loop initialisations
       ngwf_offset = 0
       mo_offset = 0

       do it=1,pub_frag_iatm(0)

          ! Use the pub_frag_data frozen kernel storage to obtain the
          ! dimensions of the matrices and eigenvalue storage:
          if ((pub_eda_mode == EDA_CTFRAGLOC_DEVEL) .and. &
              ((it .eq. pub_frag_counter2) .or. (it .eq. pub_frag_counter))) then
             ! intrafragmental (charge transfer) delocalisations:

             if (it .eq. pub_frag_counter2) then

                ! iterate fragment NGWF offset and MO offset
                ngwf_offset = ngwf_offset + pub_frag_data(it)%ngwf_basis(1)%num
                mo_offset = mo_offset + pub_frag_data(it)%rep%n_occ(is,PUB_1K)

                ! omit blocks that relate to pub_frag_counter2:
                ! (these will be grouped with pub_frag_counter fragment)
                cycle

             else if (it .eq. pub_frag_counter) then

                ! join pub_frag_counter and pub_frag_counter2 fragments:

                ! set optional second fragment ID:
                it2 = pub_frag_counter2
                jf_log = .true.

                ! number of NGWFs on this fragment pair:
                num_frag_ngwf = pub_frag_data(it)%ngwf_basis(1)%num + &
                                pub_frag_data(it2)%ngwf_basis(1)%num

                ! number of occupied MOs on this fragment pair:
                num_frag_occMO = pub_frag_data(it)%rep%n_occ(is,PUB_1K) + &
                                 pub_frag_data(it2)%rep%n_occ(is,PUB_1K)

             end if

          else
             ! standard:

             ! fragment ID2 is redundant:
             it2 = 0
             jf_log = .false.

             ! number of NGWFs on this fragment:
             num_frag_ngwf = pub_frag_data(it)%ngwf_basis(1)%num
             num_frag_occMO = pub_frag_data(it)%rep%n_occ(is,PUB_1K)

          end if


          ! mjsp: allocate eigenvalues
          allocate(eigenvalues(num_frag_ngwf),stat=ierr)
          call utils_alloc_check('scfmi_construct_nonorth_kernel_diag', &
               'eigenvalues',ierr)
          eigenvalues(:) = 0.0_DP

          ! mjsp: Initialise fragment storage.
          call dense_create(frag_dkn_dens, num_frag_ngwf, num_frag_ngwf)
          call dense_create(frag_eigs_dens, num_frag_ngwf, num_frag_ngwf)
          call dense_create(frag_overlap_dens, num_frag_ngwf, num_frag_ngwf)
          call dense_create(frag_eigs_dens_t, num_frag_ngwf, num_frag_ngwf)

          ! mjsp: Isolate the block of the density kernel that
          ! relates to this fragment.
          call fmo_get_frag_blk_dens_nxn(frag_dkn_dens, dkn_dens, it, it2, &
               frag_loc=.true.)

          ! mjsp: multiply by -1 to get eigenvecs in increasing order of occupancy
          call dense_scale(frag_dkn_dens,-1.0_DP)

          ! mjsp: Isolate the overlap block.
          call fmo_get_frag_blk_dens_nxn(frag_overlap_dens, &
               overlap_f_dens,it, it2, frag_loc=.true.)

          ! mjsp: diagonalise to obtain MO coefficients, M, for this fragment
          ! using inverse overlap:
          !    K^ab M_bi = S^ab M_bi f_i
          call timer_clock('fragment_scfmi_eigensolve',1)
          call dense_eigensolve(num_frag_occMO,eigenvalues,frag_dkn_dens, &
               frag_overlap_dens, 2, frag_eigs_dens, &
               pub_eigensolver_orfac, pub_eigensolver_abstol)
          call timer_clock('fragment_scfmi_eigensolve',2)


          ! mjsp: Put frag_eigs_dens(it) into the supermolecular MO matrix,
          ! mjsp: M_occ, via the transposed MO matrix, Mt:

          ! mjsp: Transpose fragment MO matrix so that our columns
          ! mjsp: are now NGWFs rather than MOs:
          ! mjsp: (we will be placing these columns into Mt)
          call dense_transpose(frag_eigs_dens_t,frag_eigs_dens)


          ! mjsp: Initialisations
          allocate(colbuf(num_frag_ngwf),stat=ierr)
          call utils_alloc_check('scfmi_construct_nonorth_kernel_diag','colbuf',ierr)

          ! mjsp: If joining fragments
          if (jf_log) then
             ! mjsp: Calculate skips required to place the vectors into
             ! mjsp: the fragment ordered matrix. (e.g. joining fragment A+C
             ! mjsp: must require skip for fragment B)
             join_skip_mo = 0
             join_skip_ngwf = 0
             do it_i=it+1,it2-1 ! loop fragment(s) between the fragment pair
                join_skip_mo   = join_skip_mo   + &
                     pub_frag_data(it_i)%rep%n_occ(is,PUB_1K)
                join_skip_ngwf = join_skip_ngwf + &
                     pub_frag_data(it_i)%ngwf_basis(1)%num
             end do
          end if

          ! loop NGWFs of this fragment (cols)
          do ingwf_frag=1,num_frag_ngwf

             ! mjsp: (1) obtain the (fragment-local) column for this NGWF and
             ! convert to the full fragment-adjusted ordering:

             ! mjsp: Initialisations
             col(:) = 0.0_DP
             colbuf(:) = 0.0_DP

             ! mjsp: get the fragment NGWF column and place in colbuf
             call dense_get_col(colbuf(:), frag_eigs_dens_t, ingwf_frag)

             ! if not joining two fragments
             if (.not.(jf_log)) then

                ! mjsp: extract occupied MOs
                col(mo_offset+1: &
                    mo_offset+num_frag_occMO) = &
                        colbuf(1:pub_frag_data(it)%rep%n_occ(is,PUB_1K))

                ! mjsp: get NGWF in the supermolecule indexing scheme
                ingwf_in_super = ngwf_index_map%supers_by_cum_frag( &
                     ingwf_frag + ngwf_offset)

             ! if joining two fragments
             else

                ! mjsp: split the colbuf array into the two fragments in the
                ! fragment ordered array, col
                ! FRAGMENT 1:
                col(mo_offset+1: &
                    mo_offset+pub_frag_data(it)%rep%n_occ(is,PUB_1K)) = &
                        colbuf(1:pub_frag_data(it)%rep%n_occ(is,PUB_1K))
                ! FRAGMENT 2:
                col(mo_offset+join_skip_mo+pub_frag_data(it)%rep%n_occ(is,PUB_1K)+1: &
                    mo_offset+join_skip_mo+pub_frag_data(it)%rep%n_occ(is,PUB_1K) &
                                + pub_frag_data(it2)%rep%n_occ(is,PUB_1K)) = &
                        colbuf(1+pub_frag_data(it)%rep%n_occ(is,PUB_1K): &
                               pub_frag_data(it)%rep%n_occ(is,PUB_1K) &
                                  + pub_frag_data(it2)%rep%n_occ(is,PUB_1K))

                ! mjsp: get NGWF in the supermolecule indexing scheme
                ! mjsp: if joining fragments and we are on the second fragment
                if (ingwf_frag .gt. pub_frag_data(it)%ngwf_basis(1)%num) then
                   ! mjsp: then skip is neccesary to ensure correct indexing
                   ! mjsp: of pub_frag2super_ngwf_idx.
                   ingwf_in_super = ngwf_index_map%supers_by_cum_frag( &
                        ingwf_frag + ngwf_offset + join_skip_ngwf)
                else
                   ! mjsp: else on first fragment of pair -> no skip neccesary
                   ingwf_in_super = ngwf_index_map%supers_by_cum_frag( &
                        ingwf_frag + ngwf_offset)
                end if

             end if


             ! mjsp: (2) Place column into the transposed supermolecule MO
             ! mjsp: matrix, Mt, at the supermolecule-scheme NGWF index
             call dense_put_col(col, Mt(is), ingwf_in_super)

          end do ! NGWF

          ! cleanup
          deallocate(colbuf,stat=ierr)
          call utils_dealloc_check('scfmi_construct_nonorth_kernel_diag','colbuf', &
               ierr)


          ! cleanup:
          ! mjsp: deallocate eigenvalues
          deallocate(eigenvalues,stat=ierr)
          call utils_dealloc_check('scfmi_construct_nonorth_kernel_diag', &
               'eigenvalues',ierr)

          ! mjsp: destroy temporary fragment matrices
          call dense_destroy(frag_dkn_dens)
          call dense_destroy(frag_eigs_dens)
          call dense_destroy(frag_overlap_dens)
          call dense_destroy(frag_eigs_dens_t)

          ! iterate fragment NGWF offset and MO offset
          ngwf_offset = ngwf_offset + pub_frag_data(it)%ngwf_basis(1)%num
          mo_offset = mo_offset + pub_frag_data(it)%rep%n_occ(is,PUB_1K)

       end do ! fragment

       ! mjsp: Construct M_occ from the transposed MO matrix:
       call dense_transpose(M_occ(is), Mt(is))



       ! ============ 2. CONSTRUCT KS ORBITAL INVERSE OVERLAP METRIC =========

       ! mjsp: construct M'*S (where S is the full NGWF overlap matrix)
       call dense_product(MtS(is),Mt(is),overlap_f_dens)

       ! mjsp: construct M'*S*M = overlap of KS states
       ! (shorthand='z')
       call dense_product(MtSM(is),MtS(is),M_occ(is))

       ! mjsp: construct inverse of M'*S*M
       ! shorthand='iz'
       call dense_copy(inv_MtSM(is),MtSM(is))
       call dense_invert(inv_MtSM(is))


       ! number of occupied MOs of full system (spin dependent)
       num_occMO = rep%n_occ(is,PUB_1K)

       ! construct temporary intermediate matrix using num_occMO:
       call dense_create(Miz,num,num_occMO)

       ! R = M (M'SM)^-1 M'
       call dense_product(Miz,M_occ(is),inv_MtSM(is))
       call dense_product(R,Miz,Mt(is))

       ! Convert dense R to SPAM3
       call dense_convert(denskern_R%kern%m(is,PUB_1K),R)

       ! cleanup
       call dense_destroy(Miz)

    end do ! spin

    ! cleanup
    call dense_destroy(R)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
        &scfmi_construct_nonorth_kernel_diag'

    ! Stop Timer
    call timer_clock('fragment_scfmi_construct_kern_diag',2)


  end subroutine scfmi_construct_nonorth_kernel_diag

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine scfmi_construct_nonorth_kernel_nodiag(denskern, rep)

    !======================================================================!
    ! This subroutine updates the state of the density kernel denskern_R   !
    ! from the nonorthogonal orbitals described by the denskern passed to  !
    ! this method.  The construction is equivalent to:                     !
    !   R  =  M ( M^\dagger S M)^-1 M^\dagger                              !
    ! where R is the object denskern_R, M are the MO expansion             !
    ! coefficients obtained from denskern, and S is the (full) NGWF        !
    ! overlap matrix.                                                      !
    !                                                                      !
    ! This routine avoids the diagonalisation bottleneck present within    !
    ! scfmi_construct_nonorth_kernel_diag                                  !
    !----------------------------------------------------------------------!
    ! WARNING: The current implementation requires denskern_R to have      !
    ! been updated via a call to this method before calling the other      !
    ! state-dependent methods of this module!                              !
    !----------------------------------------------------------------------!
    ! Written by Max Phipps, April 2016                                    !
    ! Modified for embedding by Joseph Prentice, September 2018            !
    !======================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, NORMAL
    use kernel, only: kernel_create, kernel_destroy
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_debug_on_root, pub_output_detail
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy, &
         sparse_embed_product, sparse_embed_array_copy
    use timer, only: timer_clock
    !use utils, only: utils_banner

    ! Arguments
    type(DKERN),intent(inout) :: denskern
    type(NGWF_REP), intent(in)   :: rep

    ! Local Variables
    integer :: pur_loop, ierr
    type(DKERN)   :: denskern_R_buf
    type(SPAM3_EMBED)   :: tmp, ls, lsl, lsls, lslsl, err
    real(kind=DP), allocatable :: min_occ(:,:), max_occ(:,:)
    real(kind=DP) :: steplen, rmse
    real(kind=DP), parameter :: rmse_conv = 1.0E-8_DP
    logical :: log_adaptive


    ! Start Timer
    call timer_clock('fragment_scfmi_construct_kern_nodiag',1)

    ! --------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
        &scfmi_construct_nonorth_kernel_nodiag'


    ! mjsp: Adaptive purification approach that guarantees
    ! mjsp: occupation number stability

    call kernel_create(denskern_R_buf, denskern%kern%m(1,PUB_1K)%structure, .false.)
    call sparse_embed_copy(denskern_R_buf%kern%m(1,PUB_1K), denskern%kern%m(1,PUB_1K))


    ! Allocate the purification matrices
    call internal_allocate(denskern_R, rep%overlap)

    if (pub_on_root.and.(pub_output_detail>NORMAL)) then
       !write(stdout,'(/a)') utils_banner('-', 'ADAPTIVE PURIFICATION')
       !write(stdout,'(/a)') '--------------------- &
       !      &ADAPTIVE PURIFICATION ---------------------'
       write(stdout,'(/a)') '-------------- &
            &INTER-FRAGMENT ADAPTIVE PURIFICATION -------------'
       write(stdout, '("|",2x,a,3x,a,3x,a,3x,a,1x,"|")') "Iteration", &
            "Extremal occupancies", "Step length", "RMSE(KSK-K)"
    end if

    log_adaptive = .true.

    ! Purify with respect to this overlap:
    pur_loop = 1
    do

       if (log_adaptive) then

          ! mjsp: Calculate the steepest descent step that will ensure the
          ! occupancy numbers are not affected.
          call internal_steplen(steplen, &
               denskern_R_buf, rep%overlap_scfmi_full)

       else
          steplen = 0.5
       end if

       ! mjsp: Do the purification step
       if (log_adaptive) then
          call internal_purify_adaptive(denskern_R%kern%m(1,PUB_1K), &
               denskern_R_buf%kern%m(1,PUB_1K), rep%overlap_scfmi_full, steplen)
       else
          call internal_purify(denskern_R%kern%m(1,PUB_1K), &
               denskern_R_buf%kern%m(1,PUB_1K), rep%overlap_scfmi_full)
       end if


       ! mjsp: Calculate convergence metric
       call sparse_embed_product(ls, denskern_R%kern%m(1,PUB_1K), rep%overlap_scfmi_full)
       call sparse_embed_product(lsl, ls, denskern_R%kern%m(1,PUB_1K))
       rmse = internal_check_conv(denskern_R%kern%m(1,PUB_1K),lsl)

       if (pub_on_root.and.(pub_output_detail>NORMAL)) then
          if (log_adaptive) then
             write(stdout,'("|",4x,i3,6x,"[",f9.6,",",f9.6,"]",3x,f9.6,5x,f10.8,2x,"|")') &
                  pur_loop, minval(min_occ), maxval(max_occ), steplen, rmse
          else
             write(stdout,'("|",4x,i3,6x,21x,3x,f9.6,5x,f10.8,2x,"|")') &
                  pur_loop, steplen, rmse
          end if
       end if

       if (rmse .lt. rmse_conv) exit

       ! mjsp: Copy the step purified kernel for the next iteration
       call sparse_embed_array_copy(denskern_R_buf%kern, denskern_R%kern)

       pur_loop = pur_loop + 1

    end do

    ! if (pub_on_root.and.(pub_output_detail>NORMAL)) utils_banner('-')
    if (pub_on_root.and.(pub_output_detail>NORMAL)) &
       write(stdout,'(a/)') '----------------------&
            &-------------------------------------------'

    call kernel_destroy(denskern_R_buf)

    ! Deallocate the purification matrices
    call internal_deallocate()

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
        &scfmi_construct_nonorth_kernel_nodiag'

    ! Stop Timer
    call timer_clock('fragment_scfmi_construct_kern_nodiag',2)

  contains

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine internal_allocate(denskern, s)

      ! mjsp: Allocates the matrices for the purification

      use sparse_embed, only: sparse_embed_create
      use utils, only: utils_alloc_check

      ! Arguments
      type(DKERN), intent(inout)  :: denskern     ! Density kernel
      type(SPAM3_EMBED), intent(in)     :: s  ! Overlap

      ! Local Variables
      integer :: nspins

      ! Adaptive purification
      call sparse_embed_create(tmp, denskern%kern%m(1,PUB_1K))
      call sparse_embed_create(ls, denskern%kern%m(1,PUB_1K), s)
      call sparse_embed_create(lsl, ls, denskern%kern%m(1,PUB_1K))
      call sparse_embed_create(lsls, ls, ls)
      call sparse_embed_create(lslsl, lsls, denskern%kern%m(1,PUB_1K))

      ! Convergence metric
      call sparse_embed_create(err, denskern%kern%m(1,PUB_1K)) ! error of matrices

      ! Step length calculation
      ! KPOINTS_DANGER
      nspins = denskern%kern%num_spins
      allocate(min_occ(nspins, denskern%kern%num_kpoints), stat=ierr)
      call utils_alloc_check('scfmi_construct_nonorth_kernel_nodiag', 'min_occ', ierr)
      allocate(max_occ(nspins, denskern%kern%num_kpoints), stat=ierr)
      call utils_alloc_check('scfmi_construct_nonorth_kernel_nodiag', 'max_occ', ierr)

    end subroutine internal_allocate

    subroutine internal_deallocate()

      ! mjsp: Deallocates the matrices for the purification

      use sparse_embed, only: sparse_embed_destroy
      use utils, only: utils_dealloc_check

      call sparse_embed_destroy(tmp)

      ! Adaptive purification
      call sparse_embed_destroy(ls)
      call sparse_embed_destroy(lsl)
      call sparse_embed_destroy(lsls)
      call sparse_embed_destroy(lslsl)

      ! Convergence metric
      call sparse_embed_destroy(err)

      ! Step length calculation
      deallocate(min_occ, stat=ierr)
      call utils_dealloc_check('scfmi_construct_nonorth_kernel_nodiag', 'min_occ', ierr)
      deallocate(max_occ, stat=ierr)
      call utils_dealloc_check('scfmi_construct_nonorth_kernel_nodiag', 'max_occ', ierr)

    end subroutine internal_deallocate

    subroutine internal_purify_adaptive(pur_k, l, s, steplen)

      ! mjsp: Adaptive purification

      use constants, only: DP
      use sparse_embed, only: sparse_embed_product, sparse_embed_scale, &
           sparse_embed_copy, sparse_embed_axpy

      ! Arguments
      type(SPAM3_EMBED), intent(inout) :: pur_k
      type(SPAM3_EMBED), intent(in) :: l
      type(SPAM3_EMBED), intent(in) :: s
      real(kind=DP), intent(in) :: steplen

      ! K = L + 2*step*(3LSL-2LSLSL-L)
      call sparse_embed_product(ls, l, s)
      call sparse_embed_product(lsl, ls, l)
      call sparse_embed_product(lslsl, ls, lsl)
      call sparse_embed_copy(pur_k, lsl)             ! pur_k = LSL
      call sparse_embed_scale(pur_k, 3.0_DP)         ! pur_k = 3LSL
      call sparse_embed_axpy(pur_k, lslsl, -2.0_DP)  ! pur_k = 3LSL-2LSLSL
      call sparse_embed_axpy(pur_k, l, -1.0_DP)      ! pur_k = 3LSL-2LSLSL-L
      call sparse_embed_scale(pur_k, 2.0_DP*steplen) ! pur_k = 2*step*(3LSL-2LSLSL-L)
      call sparse_embed_axpy(pur_k, l, 1.0_DP)       ! pur_k = L+2*step*(3LSL-2LSLSL-L)

    end subroutine internal_purify_adaptive

    subroutine internal_purify(pur_k, l, s)

      ! mjsp: McWeeny purification. This is used when occupancies within range
      ! that adaptive purification is no longer required.

      use constants, only: DP
      use sparse_embed, only:sparse_embed_product, sparse_embed_scale, &
           sparse_embed_copy, sparse_embed_axpy

      ! Arguments
      type(SPAM3_EMBED), intent(inout) :: pur_k
      type(SPAM3_EMBED), intent(in) :: l
      type(SPAM3_EMBED), intent(in) :: s

      ! K = 3LSL-2LSLSL
      call sparse_embed_product(ls, l, s)
      call sparse_embed_product(lsl, ls, l)
      call sparse_embed_product(lslsl, ls, lsl)
      call sparse_embed_copy(pur_k, lsl)             ! pur_k = LSL
      call sparse_embed_scale(pur_k, 3.0_DP)         ! pur_k = 3LSL
      call sparse_embed_axpy(pur_k, lslsl, -2.0_DP)  ! pur_k = 3LSL-2LSLSL

    end subroutine internal_purify

    function internal_check_conv(k, ksk) result(rmse)

      ! mjsp: Checks purification convergence as error of idempotency:
      ! Root mean squared error, RMSE = sqrt(av((KSK - K).^2)))

      use constants, only: DP
      use sparse_embed, only: sparse_embed_copy, sparse_embed_axpy, &
           sparse_embed_entrywise_norm, sparse_embed_num_cols

      ! Arguments
      type(SPAM3_EMBED), intent(in) :: k
      type(SPAM3_EMBED), intent(in) :: ksk
      real(kind=DP)           :: rmse

      call sparse_embed_copy(err, k)
      call sparse_embed_axpy(err, ksk, -1.0_DP)
      rmse = sparse_embed_entrywise_norm(err, 2)
      rmse = rmse * (1.0_DP/sparse_embed_num_cols(err)) ! NOTE: Correct for square matrices only

    end function internal_check_conv

    subroutine internal_steplen(steplen, denskern, overlap)

      ! mjsp: Calculates steplength that will ensure no occupancy flipping.
      ! Based on kernel_fix

      use constants, only: DP
      use kernel, only: kernel_occupancy_bounds

      ! Arguments
      real(kind=DP), intent(inout) :: steplen      ! Step length (output)
      type(DKERN), intent(inout)   :: denskern     ! Density kernel
      type(SPAM3_EMBED), intent(in)      :: overlap      ! Overlap matrix

      ! Local Variables
      real(kind=DP) :: steplen0, steplen1
      real(kind=DP), parameter :: minbound = -0.36602540378443864676_DP
      real(kind=DP), parameter :: maxbound = 1.36602540378443864676_DP

      ! Set bounds
      min_occ = 0.0_DP
      max_occ = 1.0_DP

      ! mjsp: Calculate bounds of occupancies in the full overlap
      call kernel_occupancy_bounds(max_occ, min_occ, denskern, overlap)

      ! mjsp: Calculate step length:
      ! mjsp: If occupancy out of interval [minbound,maxval] then
      ! mjsp:   steplen is the minimum of the step length required by the extremal
      ! mjsp:   occupancies over 1.0 and under 0.0 to ensure no occupancy flipping.
      ! mjsp: Else
      ! mjsp:   steplen = 0.5

      if (minval(min_occ) .lt. minbound .or. &
           maxval(max_occ) .gt. maxbound) then

         ! minimal:
         steplen0 = 1E9_DP
         if (minval(min_occ) .lt. minbound) &
            steplen0 = 1.0_DP/(4.0_DP*(minval(min_occ)-1.0_DP)*minval(min_occ))

         ! maximal:
         steplen1 = 1E9_DP
         if (maxval(max_occ) .gt. maxbound) &
            steplen1 = 1.0_DP/(4.0_DP*(maxval(max_occ)-1.0_DP)*maxval(max_occ))

         ! minimum of minimal and maximal:
         steplen = min(steplen0, steplen1)-0.001_DP

      else

         steplen = 0.5_DP

         ! if adaptive purification no longer neccesary,
         ! then use steplen=0.5 from here onwards in the purification
         if (pub_on_root.and.(pub_output_detail>NORMAL)) then
            write(stdout,'("|",a,7x,"|")') &
                 " --> Occupancies in safe range. Fixing step size to 0.5."
         end if
         log_adaptive = .false.

      end if

    end subroutine internal_steplen

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  end subroutine scfmi_construct_nonorth_kernel_nodiag

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine scfmi_hamiltonian_proj_stoll(rep, ham)

    !=========================================================================!
    ! This subroutine is used to calculate the projected Stoll Hamiltonian    !
    ! matrix. This is called internally via the wrapper                       !
    ! scfmi_hamiltonian_projection(). For further details, see:               !
    !     H. Stoll, G. Wagenblast, and H. Preuss, Theor. Chem. Acc. 57,       !
    !     169 (1980).                                                         !
    !-------------------------------------------------------------------------!
    ! rep            (in) : NGWF representation                               !
    ! ham            (in) : Hamiltonian matrices                              !
    !-------------------------------------------------------------------------!
    ! WARNING: The current implementation of this module requires             !
    ! denskern_R to have been updated via scfmi_construct_nonorth_kernel()    !
    ! before calling this routine!                                            !
    !-------------------------------------------------------------------------!
    ! Written by Max Phipps in April 2015.                                    !
    ! Modified for embedding by Joseph Prentice, September 2018               !
    !=========================================================================!

    use constants, only: stdout
    use dense, only: DEM, dense_create, dense_convert, dense_destroy, &
         dense_product, dense_transpose
    use fragment_matrix_ops, only: fmo_masks
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_debug_on_root
    use sparse_embed, only: sparse_embed_num_rows, sparse_embed_scale

    implicit none

    ! Arguments
    type(NGWF_REP), intent(in)     :: rep
    type(NGWF_HAM), intent(inout)  :: ham

    ! Local variables
    integer :: is    ! spin loop counter
    type(DEM) ::  ham_dens, proj_ham, & ! projected and unprojected Hamiltonians
      izMtS, SM, &  ! involved in projecting the Hamiltonian
      HM, HMiz, MtHM, izMtHM, izMtHMiz, SMizMtHMiz, &
      HMizMtS, SMizMtHM, SMizMtHMizMtS, &
      HMiz_MtS_t_blk, SMizMtHMiz_MtS_blk, &
      SMizMtHMiz_MtS_t_blk, SM_izMtHMiz_MtS_blk, &
      HMizMtS_t_blk, HMiz_MtS_blk, &
      H_blk, HMizMtS_blk, & ! fragment matrix blocks
      SMizMtHMizMtS_blk, HMiz_blk, &
      MtS_blk, SMizMtHMiz_blk, SM_blk, izMtHMiz_blk
    integer :: num_ngwf, num_occMO
    type(DEM) :: tmpNM_blk   ! N=NGWFs, M=(occupied) MOs

    ! --------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
        &scfmi_hamiltonian_proj_stoll'

    num_ngwf = sparse_embed_num_rows(ham%ham(1))

    call dense_create(  HMizMtS,num_ngwf,num_ngwf)
    call dense_create(SMizMtHMizMtS,num_ngwf,num_ngwf)
    call dense_create(ham_dens,num_ngwf,num_ngwf)
    call dense_create(proj_ham,num_ngwf,num_ngwf)
    call dense_create(  H_blk,      num_ngwf,num_ngwf)
    call dense_create(  HMizMtS_blk,num_ngwf,num_ngwf)
    call dense_create(SMizMtHMizMtS_blk,num_ngwf,num_ngwf)
    call dense_create(HMiz_MtS_blk,num_ngwf,num_ngwf)
    call dense_create(HMiz_MtS_t_blk,num_ngwf,num_ngwf)
    call dense_create(SMizMtHMiz_MtS_blk,num_ngwf,num_ngwf)
    call dense_create(SM_izMtHMiz_MtS_blk,num_ngwf,num_ngwf)
    call dense_create(HMizMtS_t_blk,num_ngwf,num_ngwf)
    call dense_create(SMizMtHMiz_MtS_t_blk,num_ngwf,num_ngwf)

    ! mjsp: loop over spins
    do is=1,pub_num_spins

       ! number of occupied MOs of full system (spin dependent)
       num_occMO = rep%n_occ(is,PUB_1K)

       ! create the intermediate matrices with their
       ! appropriate dimensions (supermolecule matrices)
       call dense_create(  HM,num_ngwf,num_occMO)
       call dense_create(MtHM,num_occMO,num_occMO)
       call dense_create(  HMiz,num_ngwf,num_occMO)
       call dense_create(    izMtS,num_occMO,num_ngwf)
       call dense_create(  izMtHM,num_occMO,num_occMO)
       call dense_create(  izMtHMiz,num_occMO,num_occMO)
       call dense_create(SM,num_ngwf,num_occMO)
       call dense_create(SMizMtHM,num_ngwf,num_occMO)
       call dense_create(SMizMtHMiz,num_ngwf,num_occMO)

       call dense_create(  HMiz_blk,    num_ngwf,num_occMO)  ! NxM
       call dense_create(  MtS_blk,   num_occMO,num_ngwf)    ! MxN
       call dense_create(izMtHMiz_blk,num_occMO,num_occMO)   ! MxM
       call dense_create(SMizMtHMiz_blk,num_ngwf,num_occMO)  ! NxM
       call dense_create(SM_blk,        num_ngwf,num_occMO)  ! NxM

       ! initialise temporary matrix
       call dense_create(tmpNM_blk,num_ngwf,num_occMO)


       ! expand Hamiltonian to dense matrix
       call dense_convert(ham_dens,ham%ham(is))

       ! clear the full Hamiltonian
       call sparse_embed_scale(ham%ham(is), 0.0_DP)


       ! (1) Calculate components of the sub-blocks of the projected Hamiltonian
       ! NOTE: SMiz, and MtS already constructed

       ! izMtS
       call dense_product(izMtS,inv_MtSM(is),MtS(is))

       ! SM
       call dense_product(SM,overlap_f_dens,M_occ(is))

       ! HM
       call dense_product(HM,ham_dens,M_occ(is))
       ! HMiz =  H M z^-1
       call dense_product(HMiz,HM,inv_MtSM(is))

       ! MtHM
       call dense_product(MtHM,Mt(is),HM)

       !   izMtHM
       call dense_product(izMtHM,inv_MtSM(is),MtHM)
       ! SMizMtHM
       call dense_product(SMizMtHM,SM,izMtHM)


       ! (2) Calculate sub-blocks of the projected Hamiltonian
       call dense_product(HMizMtS,HMiz,MtS(is))
       call dense_product(  izMtHMiz,izMtHM,inv_MtSM(is))
       call dense_product(SMizMtHMiz,SM,izMtHMiz)
       call dense_product(SMizMtHMizMtS,SMizMtHMiz,MtS(is))

       ! mjsp: Use fragment masks to 'fragment' the matrices
       H_blk%dmtx(:,:) = ham_dens%dmtx(:,:) * &
            fmo_masks%nxn%dmtx(:,:)
       HMizMtS_blk%dmtx(:,:) = HMizMtS%dmtx(:,:) * &
            fmo_masks%nxn%dmtx(:,:)
       SMizMtHMizMtS_blk%dmtx(:,:) = SMizMtHMizMtS%dmtx(:,:) * &
            fmo_masks%nxn%dmtx(:,:)
       HMiz_blk%dmtx(:,:) = HMiz%dmtx(:,:) * &
            fmo_masks%nxm(is)%dmtx(:,:)
       MtS_blk%dmtx(:,:) = MtS(is)%dmtx(:,:) * &
            fmo_masks%mxn(is)%dmtx(:,:)
       izMtHMiz_blk%dmtx(:,:) = izMtHMiz%dmtx(:,:) * &
            fmo_masks%mxm(is)%dmtx(:,:)
       SMizMtHMiz_blk%dmtx(:,:) = SMizMtHMiz%dmtx(:,:) * &
            fmo_masks%nxm(is)%dmtx(:,:)
       SM_blk%dmtx(:,:) = SM%dmtx(:,:) * &
            fmo_masks%nxm(is)%dmtx(:,:)

       ! (3) Construct the required block products and transposes
       call dense_transpose(HMizMtS_t_blk,HMizMtS_blk)

       ! product of frag blocks:
       call dense_product(HMiz_MtS_blk,HMiz_blk,MtS_blk)

       call dense_transpose(HMiz_MtS_t_blk,HMiz_MtS_blk)

       ! product of frag blocks:
       call dense_product(SMizMtHMiz_MtS_blk,SMizMtHMiz_blk,MtS_blk)

       call dense_transpose(SMizMtHMiz_MtS_t_blk,SMizMtHMiz_MtS_blk)

       ! product of frag blocks:
       call dense_product(tmpNM_blk,SM_blk,izMtHMiz_blk)
       call dense_product(SM_izMtHMiz_MtS_blk,tmpNM_blk,MtS_blk)


       ! (4) Construct the projected Hamiltonian block from these components
       proj_ham%dmtx(:,:) = 0.0_DP
       proj_ham%dmtx(:,:) = proj_ham%dmtx(:,:) + H_blk%dmtx(:,:)
       proj_ham%dmtx(:,:) = proj_ham%dmtx(:,:) - HMizMtS_blk%dmtx(:,:)
       proj_ham%dmtx(:,:) = proj_ham%dmtx(:,:) - HMizMtS_t_blk%dmtx(:,:)
       proj_ham%dmtx(:,:) = proj_ham%dmtx(:,:) + SMizMtHMizMtS_blk%dmtx(:,:)
       proj_ham%dmtx(:,:) = proj_ham%dmtx(:,:) + HMiz_MtS_blk%dmtx(:,:)
       proj_ham%dmtx(:,:) = proj_ham%dmtx(:,:) + HMiz_MtS_t_blk%dmtx(:,:)
       proj_ham%dmtx(:,:) = proj_ham%dmtx(:,:) - SMizMtHMiz_MtS_blk%dmtx(:,:)
       proj_ham%dmtx(:,:) = proj_ham%dmtx(:,:) - SMizMtHMiz_MtS_t_blk%dmtx(:,:)
       proj_ham%dmtx(:,:) = proj_ham%dmtx(:,:) + SM_izMtHMiz_MtS_blk%dmtx(:,:)


       ! copy the now fully projected Hamiltonian back into the
       ! sparse storage for this spin
       call dense_convert(ham%ham(is),proj_ham)

       ! mjsp: cleanup:
       call dense_destroy(tmpNM_blk)

       call dense_destroy(  HM)
       call dense_destroy(MtHM)
       call dense_destroy(  HMiz)
       call dense_destroy(    izMtS)
       call dense_destroy(  HMizMtS)
       call dense_destroy(  izMtHM)
       call dense_destroy(  izMtHMiz)
       call dense_destroy(SM)
       call dense_destroy(SMizMtHM)
       call dense_destroy(SMizMtHMiz)

    end do

    ! mjsp: cleanup:
    call dense_destroy(  HMizMtS_blk)
    call dense_destroy(SMizMtHMizMtS)
    call dense_destroy(ham_dens)
    call dense_destroy(proj_ham)
    call dense_destroy(  H_blk)
    call dense_destroy(SMizMtHMizMtS_blk)
    call dense_destroy(  HMiz_blk )
    call dense_destroy(  MtS_blk )
    call dense_destroy(izMtHMiz_blk)
    call dense_destroy(SMizMtHMiz_blk)
    call dense_destroy(SM_blk)
    call dense_destroy(HMiz_MtS_blk)
    call dense_destroy(HMiz_MtS_t_blk)
    call dense_destroy(SMizMtHMiz_MtS_blk)
    call dense_destroy(SM_izMtHMiz_MtS_blk)
    call dense_destroy(HMizMtS_t_blk)
    call dense_destroy(SMizMtHMiz_MtS_t_blk)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
        &scfmi_hamiltonian_proj_stoll'

  end subroutine scfmi_hamiltonian_proj_stoll

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine scfmi_hamiltonian_proj_stoll_nodiag(rep, ham)

    !=========================================================================!
    ! This subroutine is used to calculate the projected Stoll Hamiltonian    !
    ! matrix whilst avoiding diagonalisation. This is called internally via   !
    ! the wrapper scfmi_hamiltonian_projection(). For further details, see:   !
    !     H. Stoll, G. Wagenblast, and H. Preuss, Theor. Chem. Acc. 57,       !
    !     169 (1980).                                                         !
    !-------------------------------------------------------------------------!
    ! ham            (in) : Hamiltonian matrices                              !
    !-------------------------------------------------------------------------!
    ! WARNING: The current implementation of this module requires             !
    ! denskern_R to have been updated via scfmi_construct_nonorth_kernel()    !
    ! before calling this routine!                                            !
    !-------------------------------------------------------------------------!
    ! Written by Max Phipps in April 2016.                                    !
    ! Modified for embedding by Joseph Prentice, September 2018               !
    !=========================================================================!

    use constants, only: stdout
    use fragment_matrix_ops, only: fmo_masks
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_debug_on_root, pub_frag_iatm
    use sparse_embed, only: sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_transpose, sparse_embed_copy, sparse_embed_num_rows, &
         sparse_embed_scale, sparse_embed_product, sparse_embed_axpy

    implicit none

    ! Arguments
    type(NGWF_REP), intent(in   ) :: rep
    type(NGWF_HAM), intent(inout)  :: ham

    ! Local variables
    integer :: is    ! spin loop counter
    type(SPAM3_EMBED) ::  proj_ham, & ! unprojected Hamiltonian
      Rr, & ! Row-fragmented density kernel R
      RrS, HRrS, HRrS_blk, HRrS_t_blk, &
      SRHRrS, SRHRrS_blk, SRHRrS_t_blk, &
      Rrt, SRrt, SRrtHRrS, SRrtHRrS_blk, & ! t indicates transpose
      RS, HRS, SR, SRHRS, &
      HRS_t_blk, &
      H_blk, HRS_blk, & ! fragment matrix blocks
      SRHRS_blk
    integer :: num_ngwf, it

    ! --------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
        &scfmi_hamiltonian_proj_stoll_nodiag'

    num_ngwf = sparse_embed_num_rows(ham%ham(1))

    call sparse_embed_create(proj_ham, ham%ham(1))

    Rr%structure = 'K'
    call sparse_embed_create(   Rr)
    RrS%structure = 'KS'
    call sparse_embed_create(   RrS)
    call sparse_embed_create(  HRrS,        ham%ham(1), RrS)
    Rrt%structure = 'K'
    call sparse_embed_create(      Rrt)
    SRrt%structure = 'SK'
    call sparse_embed_create(     SRrt)
    call sparse_embed_create(     SRrtHRrS, SRrt, HRrS)
    SR%structure = 'SK'
    call sparse_embed_create(    SR)
    call sparse_embed_create(SRHRrS,        SR, HRrS)

    call sparse_embed_create(  H_blk,      ham%ham(1))
    RS%structure = 'KS'
    call sparse_embed_create(       RS)
    call sparse_embed_create(      HRS,    ham%ham(1), RS)
    call sparse_embed_create(    SRHRS,    SR, HRS)
    call sparse_embed_create(  HRS_blk,    HRS)
    call sparse_embed_create(  HRS_t_blk,  HRS)
    call sparse_embed_create(  HRrS_blk,   HRS)
    call sparse_embed_create(  HRrS_t_blk, HRS)
    call sparse_embed_create(SRHRS_blk,    SRHRS)
    call sparse_embed_create(SRHRrS_t_blk, SRHRS)
    call sparse_embed_create(SRHRrS_blk,   SRHRS)
    call sparse_embed_create(SRrtHRrS_blk, SRHRS)

    ! mjsp: loop over spins
    do is=1,pub_num_spins

       ! Calculate sub-blocks of the projected Hamiltonian
       call sparse_embed_product(RS,denskern_R%kern%m(is,PUB_1K), rep%overlap_scfmi_full)

       call sparse_embed_product(HRS,ham%ham(is),RS)

       !call sparse_product(SR,rep%overlap,denskern_R)
       call sparse_embed_transpose(SR,RS) ! SR=(RS)'

       call sparse_embed_product(SRHRS,SR,HRS)

       ! mjsp: Clear the Rr dependent buffers
       ! Rr = row-fragmented R
       call sparse_embed_scale(HRrS_blk, 0.0_DP)
       call sparse_embed_scale(SRHRrS_blk, 0.0_DP)
       call sparse_embed_scale(SRrtHRrS_blk, 0.0_DP) ! 't' indicates transpose

       ! mjsp: Construct the matrices that are built using
       ! mjsp: the row-fragmented (column intact) kernel
       do it=1,pub_frag_iatm(0) ! loop fragments

         ! mjsp: Obtain column-fragmented density kernel, R,
         ! mjsp: for this fragment, it.
         Rr%p%dmtx(:) = denskern_R%kern%m(is,PUB_1K)%p%dmtx(:) * &
              fmo_masks%nxn_srK(it)%p%dmtx(:)

         ! Rr dependent matrix components:
         ! HRrS
         call sparse_embed_product(RrS,Rr,rep%overlap_scfmi_full)
         call sparse_embed_product(HRrS,ham%ham(is),RrS)

         ! SRHRrS
         call sparse_embed_product(SRHRrS,SR,HRrS)

         ! SRrtHRrS
         call sparse_embed_transpose(Rrt,Rr)
         call sparse_embed_product(SRrt,rep%overlap_scfmi_full,Rrt)
         call sparse_embed_product(SRrtHRrS,SRrt,HRrS)

         ! mjsp: Extract the fragment
         ! mjsp: matrix block (fragment it)
         HRrS%p%dmtx(:) = HRrS%p%dmtx(:) * &
              fmo_masks%nxn_scHKS(it)%p%dmtx(:)
         HRrS%p%dmtx(:) = HRrS%p%dmtx(:) * &
              fmo_masks%nxn_srHKS(it)%p%dmtx(:)
         SRHRrS%p%dmtx(:) = SRHRrS%p%dmtx(:) * &
              fmo_masks%nxn_scSKHKS(it)%p%dmtx(:)
         SRHRrS%p%dmtx(:) = SRHRrS%p%dmtx(:) * &
              fmo_masks%nxn_srSKHKS(it)%p%dmtx(:)
         SRrtHRrS%p%dmtx(:) = SRrtHRrS%p%dmtx(:) * &
              fmo_masks%nxn_scSKHKS(it)%p%dmtx(:)
         SRrtHRrS%p%dmtx(:) = SRrtHRrS%p%dmtx(:) * &
              fmo_masks%nxn_srSKHKS(it)%p%dmtx(:)

         ! mjsp: Add the fragment block to the supermolecule
         ! mjsp: matrix component
         HRrS_blk%p%dmtx(:) = HRrS_blk%p%dmtx(:) + &
              HRrS%p%dmtx(:)
         SRHRrS_blk%p%dmtx(:) = SRHRrS_blk%p%dmtx(:) + &
              SRHRrS%p%dmtx(:)
         SRrtHRrS_blk%p%dmtx(:) = SRrtHRrS_blk%p%dmtx(:) + &
              SRrtHRrS%p%dmtx(:)

       end do

       ! mjsp: Use fragment masks to 'fragment' the matrices
       H_blk%p%dmtx(:) = ham%ham(is)%p%dmtx(:) * &
            fmo_masks%nxn_sH%p%dmtx(:)
       HRS_blk%p%dmtx(:) = HRS%p%dmtx(:) * &
            fmo_masks%nxn_sHKS%p%dmtx(:)
       SRHRS_blk%p%dmtx(:) = SRHRS%p%dmtx(:) * &
            fmo_masks%nxn_sSKHKS%p%dmtx(:)

       ! (4) Construct the required block products and transposes
       call sparse_embed_transpose(HRS_t_blk, HRS_blk)
       call sparse_embed_transpose(HRrS_t_blk,HRrS_blk)
       call sparse_embed_transpose(SRHRrS_t_blk,SRHRrS_blk)

       ! (5) Construct the projected Hamiltonian block from these components
       call sparse_embed_scale(proj_ham, 0.0_DP)
       call sparse_embed_axpy(proj_ham, H_blk, 1.0_DP)

       call sparse_embed_axpy(proj_ham, HRS_blk, -1.0_DP)
       call sparse_embed_axpy(proj_ham, HRS_t_blk, -1.0_DP)

       call sparse_embed_axpy(proj_ham, SRHRS_blk, +1.0_DP)

       call sparse_embed_axpy(proj_ham, HRrS_blk, +1.0_DP)
       call sparse_embed_axpy(proj_ham, HRrS_t_blk, +1.0_DP)

       call sparse_embed_axpy(proj_ham, SRHRrS_blk, -1.0_DP)
       call sparse_embed_axpy(proj_ham, SRHRrS_t_blk, -1.0_DP)

       call sparse_embed_axpy(proj_ham, SRrtHRrS_blk, 1.0_DP)

       ! copy the now fully projected Hamiltonian back into the
       ! sparse storage for this spin
       call sparse_embed_copy(ham%ham(is),proj_ham)

    end do

    ! mjsp: cleanup:
    call sparse_embed_destroy(proj_ham)

    call sparse_embed_destroy(Rr)
    call sparse_embed_destroy(RrS)
    call sparse_embed_destroy(HRrS)
    call sparse_embed_destroy( Rrt)
    call sparse_embed_destroy(SRrt)
    call sparse_embed_destroy(SRrtHRrS)
    call sparse_embed_destroy(SRHRrS)

    call sparse_embed_destroy(  H_blk)
    call sparse_embed_destroy(  RS)
    call sparse_embed_destroy(  SR)
    call sparse_embed_destroy(SRHRS)
    call sparse_embed_destroy(  HRS)
    call sparse_embed_destroy(  HRS_blk)
    call sparse_embed_destroy(HRS_t_blk)
    call sparse_embed_destroy(  HRrS_blk)
    call sparse_embed_destroy(  HRrS_t_blk)
    call sparse_embed_destroy(SRHRS_blk)
    call sparse_embed_destroy(SRHRrS_t_blk)
    call sparse_embed_destroy(SRHRrS_blk)
    call sparse_embed_destroy(SRrtHRrS_blk)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
        &scfmi_hamiltonian_proj_stoll_nodiag'

  end subroutine scfmi_hamiltonian_proj_stoll_nodiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine scfmi_ngwf_gradient(pmat, qmat, denskern, ham, inv_overlap)

    !==========================================================================!
    ! Calculates the NGWF gradient during fragment SCF-MI calculations.        !
    ! This is currently equivalent to the unconstrained gradient (SCF-MI       !
    ! constraints are currently placed only on the density kernel).            !
    !--------------------------------------------------------------------------!
    ! Written by Max Phipps in March 2015 (modified from edft_ngwf_gradient).  !
    ! Modified for embedding by Joseph Prentice, September 2018                !
    !==========================================================================!

    use constants, only: DP, stdout
    use rundat, only: pub_num_spins, pub_debug_on_root
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product, sparse_embed_copy, &
         sparse_embed_scale

    implicit none

    type(SPAM3_EMBED), intent(inout) :: pmat(1:pub_num_spins)
    type(SPAM3_EMBED), intent(inout) :: qmat(1:pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: denskern(1:pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: ham(1:pub_num_spins)
    type(SPAM3_EMBED), intent(in   ) :: inv_overlap

    integer :: is
    type(SPAM3_EMBED) :: kh

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: Entering scfmi_ngwf_gradient"

    ! mjsp: create temporary structures
    call sparse_embed_create(kh, denskern(1), ham(1))

    do is = 1, pub_num_spins

       ! mjsp: copy denskern to pmat
       call sparse_embed_copy(pmat(is), denskern(is))

       ! mjsp: calculate qmat
       call sparse_embed_product(kh, denskern(is), ham(is))
       call sparse_embed_product(qmat(is), kh, inv_overlap)
       call sparse_embed_scale(qmat(is), -1.0_DP)

    end do

    ! mjsp: destroy structures
    call sparse_embed_destroy(kh)

    if (pub_debug_on_root) write(stdout,'(a)') &
         "DEBUG: Leaving scfmi_ngwf_gradient"

  end subroutine scfmi_ngwf_gradient

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module fragment_scfmi
