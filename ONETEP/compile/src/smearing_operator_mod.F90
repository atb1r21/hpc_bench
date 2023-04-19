! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                                        S M E A R I N G         !
!                                     ____ ___  ____  ____/ /    !
!                                    / __ `__ \/ __ \/ __  /     !
!                                   / / / / / / /_/ / /_/ /      !
!                                  /_/ /_/ /_/\____/\__,_/       !
!                                                                !
!                   Smearing Operator Module                     !
!                         Version 1.0                            !
!                         06 Feb 2014                            !
!                                                                !
!                   Written by Jolyon Aarons                     !
!                                                                !
!----------------------------------------------------------------!
!                                                                !
! This module applies smearing functions to Hamiltonians,        !
! calculating the corresponding density matrix, without          !
! having to diagonalise and avoiding cubically scaling steps.    !
!                                                                !
! Notes:                                                         !
! As of Nov 2014, the extremal eigenvalue bottlenecks have been  !
! mostly fixed. MPI_Bcast are still quite slow on some clusters  !
! Iridis in particular, but Archer is also affected.             !
!                                                                !
! As of Jun 2014, Fermi operator expansions are still the only   !
! operators available, but a better performing option has been   !
! implemented, as has support for optimising the chemical        !
! potential with respect to a number of particles constraint.    !
!                                                                !
! As of Feb 2014, the only operator that has been implemented is !
! a Fermi operator expansion in Chebyshev polynomials.           !
!                                                                !
! The routine has been written to operate on an abstract         !
! "matrix" type, so that interfaces could be written to both     !
! SPAM and DEM matrices, as well as standard, serial (FORTRAN)   !
! matrices, while performing exactly the same high-level         !
! abstract operations.                                           !
!                                                                !
! Modified for embedding by Joseph Prentice, June 2018           !
! SPAM3_EMBED routines added by Robert Charlton, August 2018.    !
!                                                                !
! Modified to include routines for the mermin method by Emiliano !
! Poli, January 2022                                             !
!================================================================!
module smearing_operator
  use constants, only: DP, stdout, pi
  use dense, only: DEM
  use sparse, only: SPAM3
  use sparse_embed, only: SPAM3_EMBED
  implicit none

  public apply_fermi_function
  public calculate_fermi_entropy
  public lowdin_transformation
  public smearing_matrix
  !ep : subroutine to calculate entropy energy/derivatives
  public calculate_fermi_entropy_mermin

  private
  logical :: debug_info_toggle = .false. ! spam the screen with loads of debug info if true.
  integer, parameter :: matrix_type_standard=1
  integer, parameter :: matrix_type_DEM=2
  integer, parameter :: matrix_type_SPAM3=3
  integer, parameter :: matrix_type_SPAM3_EMBED=4 ! jcap: for embedding

  character(len=16), parameter :: SIGNUM_PROJECTION="SIGNUM_PROJECTION"
  character(len=16), parameter :: CONTOUR_PROJECTIO="CONTOUR_PROJECTIO"
  character(len=16), parameter :: CERIOTTI_KUHNE_PA="CERIOTTI_KUHNE_PA"
  character(len=16), parameter :: LIANG_BAER_SARAVA="LIANG_BAER_SARAVA"
  character(len=16), parameter :: LIN_CONTOUR_INTEG="LIN_CONTOUR_INTEG"
  character(len=16), parameter :: ANNEALING        ="ANNEALING_METHOD2"
  character(len=16) :: PUB_FOE_TYPE=ANNEALING !CONTOUR_PROJECTIO

  type smearing_matrix
     integer :: matrix_type=matrix_type_standard
     logical :: standard_is_cmplx=.false.
     logical :: matrix_is_allocated=.false.
     real(kind=dp), dimension(:,:), allocatable :: data
     complex(kind=dp), dimension(:,:), allocatable :: zdata
     type(DEM)   :: dataDEM
     type(SPAM3) :: dataSPAM3
     type(SPAM3_EMBED) :: dataSPAM3_EMBED ! jcap: for embedding
  end type smearing_matrix
  type vector
     integer :: matrix_type=matrix_type_standard
     logical :: standard_is_cmplx=.false.
     real(kind=dp), dimension(:), allocatable :: data
     complex(kind=dp), dimension(:), allocatable :: zdata
     type(DEM)   :: dataDEM
     type(SPAM3) :: dataSPAM3
     type(SPAM3_EMBED) :: dataSPAM3_EMBED ! jcap: for embedding
  end type vector

  interface apply_fermi_function
     module procedure apply_fermi_function_standard
     module procedure apply_fermi_function_dem
     module procedure apply_fermi_function_spam3
     module procedure apply_fermi_function_spam3_embed
     module procedure apply_fermi_function_CCH_spam3                ! contracovariant Hamiltonian interface
     module procedure apply_fermi_function_CCH_spam3_embed          ! contracovariant Hamiltonian interface
  end interface

  interface calculate_fermi_entropy
     module procedure calculate_fermi_entropy_standard
     module procedure calculate_fermi_entropy_dem
     module procedure calculate_fermi_entropy_spam3
     module procedure calculate_fermi_entropy_spam3_embed
  end interface

  interface matrix_allocate
     module procedure matrix_allocate_num
     module procedure matrix_allocate_mat
  end interface
  interface vector_allocate
     module procedure vector_allocate_num
     module procedure vector_allocate_vec
  end interface

  interface jacobi_elliptic
     module procedure jacobi_elliptic_real
     module procedure jacobi_elliptic_complex
  end interface jacobi_elliptic

contains

  !============================================================================!
  ! This subroutine applies the Fermi operator to a given Hamiltonian and      !
  ! overlap matrix, and returns the density matrix. The Fermi-level and the    !
  ! smearing temperature must be given. This form of the routine performs the  !
  ! operation on standard, serial matrices.                                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The Hamiltonian matrix                           !
  !   overlap      (input)  : The overlap matrix                               !
  !   fermi_level  (input)  : The Fermi level                                  !
  !   smearing     (input)  : The smearing temerature                          !
  !   densmat      (output) : The computed density matrix                      !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Feb 2014.                                        !
  !----------------------------------------------------------------------------!
  ! Seperated common operations into apply_fermi_function_common, Feb 2014.    !
  !============================================================================!
  subroutine apply_fermi_function_standard(H,overlap,inv_overlap,fermi_level,&
       & smearing,densmat,correct_chemical_potential, Ne, boundl, boundu, &
       & fermi_tol, entropy, entropy_approx)
    use comms, only : pub_on_root
    use utils, only : utils_assert
    implicit none
    real(kind=dp), dimension(:,:), intent(in)  :: H
    real(kind=dp), dimension(:,:), intent(in)  :: overlap
    real(kind=dp), dimension(:,:), intent(in)  :: inv_overlap
    real(kind=dp),                 intent(inout)  :: fermi_level
    real(kind=dp),                 intent(in)  :: smearing
    real(kind=dp), dimension(:,:), intent(out) :: densmat
    logical,       optional,       intent(in)     :: correct_chemical_potential
    real(kind=dp), optional,       intent(in)     :: Ne
    real(kind=dp), optional,       intent(in)     :: boundl
    real(kind=dp), optional,       intent(in)     :: boundu
    real(kind=dp), optional,       intent(in)     :: fermi_tol
    real(kind=dp), optional,       intent(out)    :: entropy
    integer,       optional,       intent(in)     :: entropy_approx

    integer :: N
    type(smearing_matrix) :: ham
    type(smearing_matrix) :: rho
    type(smearing_matrix) :: M
    type(smearing_matrix) :: invS

    logical       :: loc_correct
    real(kind=dp) :: loc_ubound, loc_lbound, loc_Ne, loc_fermi_tol
    integer       :: loc_entropy_approx
    real(kind=dp) :: traces(2)


    if(present(correct_chemical_potential)) then
       loc_correct=correct_chemical_potential
       if(loc_correct) then
          call utils_assert(present(Ne), "Error in apply_fermi_function_common : Must supply Ne to correct chemical potential!")
          loc_Ne=Ne

          if(present(boundl)) then
             loc_lbound=boundl
          else
             write(stdout,*) "Smearing-> apply_fermi_function_common : lbound not supplied, using -5.0"
             loc_lbound=-5.0_dp
          end if
          if(present(boundu)) then
             loc_ubound=boundu
          else
             write(stdout,*) "Smearing-> apply_fermi_function_common : ubound not supplied, using 5.0"
             loc_ubound=5.0_dp
          end if
          if(present(fermi_tol)) then
             loc_fermi_tol=fermi_tol
          else
             write(stdout,*) "Smearing-> apply_fermi_function_common : fermi_tol not supplied, using 1e-7"
             loc_fermi_tol=1e-7_dp
          end if
       end if
    else
       loc_correct=.false.
    end if



    N=size(H,1)
    ham%matrix_type=matrix_type_standard
    call matrix_allocate(ham,N,N)
    ham%data=H
    call matrix_allocate(rho,ham)
    call matrix_allocate(M,ham)
    M%data=overlap
    call matrix_allocate(invS,ham)
    invS%data=inv_overlap

    traces(1)=matrix_trace(ham)
    traces(2)=matrix_old_trace(ham)
    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing->", traces, "<== SMEAR"
    end if

    if(present(entropy_approx)) then
       loc_entropy_approx=entropy_approx
    else
       loc_entropy_approx=0
    end if

    if(present(entropy)) then
       call apply_fermi_function_common(ham,M,invS,fermi_level,smearing,rho,loc_correct,&
            &loc_Ne,loc_lbound,loc_ubound,loc_fermi_tol,entropy,loc_entropy_approx)
    else
       call apply_fermi_function_common(ham,M,invS,fermi_level,smearing,rho,loc_correct,&
            &loc_Ne,loc_lbound,loc_ubound,loc_fermi_tol)
    end if

    densmat=rho%data
    call matrix_free(rho)
    call matrix_free(ham)
    call matrix_free(M)
    call matrix_free(invS)
  end subroutine apply_fermi_function_standard

  !============================================================================!
  ! This subroutine applies the Fermi operator to a given Hamiltonian and      !
  ! overlap matrix, and returns the density matrix. The Fermi-level and the    !
  ! smearing temperature must be given. This form of the routine performs the  !
  ! operation on ONETEP dense matrices.                                        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The Hamiltonian matrix                           !
  !   overlap      (input)  : The overlap matrix                               !
  !   fermi_level  (input)  : The Fermi level                                  !
  !   smearing     (input)  : The smearing temerature                          !
  !   densmat      (output) : The computed density matrix                      !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Feb 2014.                                        !
  !----------------------------------------------------------------------------!
  ! Seperated common operations into apply_fermi_function_common, Feb 2014.    !
  !============================================================================!
  subroutine apply_fermi_function_dem(H,overlap,inv_overlap,fermi_level, &
       & smearing,densmat, correct_chemical_potential, Ne, boundl, boundu, &
       & fermi_tol, entropy, entropy_approx)
    use utils, only : utils_abort, utils_assert
    use comms, only : pub_on_root
    use dense, only : dense_copy
    implicit none
    type(DEM),        intent(in)     :: H
    type(DEM),        intent(in)     :: overlap
    type(DEM),        intent(in)     :: inv_overlap
    real(kind=dp),    intent(inout)  :: fermi_level
    real(kind=dp),    intent(in)     :: smearing
    type(DEM),        intent(inout)  :: densmat
    logical,       optional,       intent(in)     :: correct_chemical_potential
    real(kind=dp), optional,       intent(in)     :: Ne
    real(kind=dp), optional,       intent(in)     :: boundl
    real(kind=dp), optional,       intent(in)     :: boundu
    real(kind=dp), optional,       intent(in)     :: fermi_tol
    real(kind=dp), optional,       intent(out)    :: entropy
    integer,       optional,       intent(in)     :: entropy_approx
    integer       :: loc_entropy_approx
    integer :: N
    type(smearing_matrix) :: ham
    type(smearing_matrix) :: rho
    type(smearing_matrix) :: M
    type(smearing_matrix) :: invS
    logical       :: loc_correct
    real(kind=dp) :: loc_ubound, loc_lbound, loc_Ne, loc_fermi_tol
    ! agrecocmplx
    logical :: loc_cmplx

    ! agrecocmplx
    loc_cmplx = overlap%iscmplx

    if(present(correct_chemical_potential)) then
       loc_correct=correct_chemical_potential
       if(loc_correct) then
          call utils_assert(present(Ne), "Error in apply_fermi_function_common : Must supply Ne to correct chemical potential!")
          loc_Ne=Ne

          if(present(boundl)) then
             loc_lbound=boundl
          else
             write(stdout,*) "smearing-> apply_fermi_function_common : lbound not supplied, using -5.0"
             loc_lbound=-5.0_dp
          end if
          if(present(boundu)) then
             loc_ubound=boundu
          else
             write(stdout,*) "smearing-> apply_fermi_function_common : ubound not supplied, using 5.0"
             loc_ubound=5.0_dp
          end if
          if(present(fermi_tol)) then
             loc_fermi_tol=fermi_tol
          else
             write(stdout,*) "smearing-> apply_fermi_function_common : fermi_tol not supplied, using 1e-7"
             loc_fermi_tol=1e-7_dp
          end if
       end if
    else
       loc_correct=.false.
    end if

    N=H%nrows
    if(N/=H%mcols) then
       call utils_abort('Error in apply_fermi_function_dem: H matrix not square!')
    end if
    ham%matrix_type=matrix_type_dem
    ! agrecocmplx
    ham%standard_is_cmplx = loc_cmplx

    call matrix_allocate(ham,N,N)

    call dense_copy(ham%dataDEM,H)

    call matrix_allocate(rho,ham)
    call matrix_allocate(M,ham)
    if(N/=M%dataDEM%nrows) then
       call utils_abort('Error in apply_fermi_function_dem: H and overlap dims not consistent!')
    end if
    if(N/=M%dataDEM%mcols) then
       call utils_abort('Error in apply_fermi_function_dem: overlap matrix not square!')
    end if

    call dense_copy(M%dataDEM,overlap)

    call matrix_allocate(invS,ham)
    call dense_copy(invS%dataDEM, inv_overlap)

    if(present(entropy_approx)) then
       loc_entropy_approx=entropy_approx
    else
       loc_entropy_approx=0
    end if

    if(present(entropy)) then
       call apply_fermi_function_common(ham,M,invS,fermi_level,smearing,rho,loc_correct,&
            &loc_Ne,loc_lbound,loc_ubound,loc_fermi_tol,entropy,loc_entropy_approx)
    else
       call apply_fermi_function_common(ham,M,invS,fermi_level,smearing,rho,loc_correct,&
            &loc_Ne,loc_lbound,loc_ubound,loc_fermi_tol)
    end if


    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> enter dense copy", "<== SMEAR"
    end if
    call dense_copy(densmat,rho%dataDEM)
    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> leave dense copy", "<== SMEAR"
    end if

    call matrix_free(rho)
    call matrix_free(ham)
    call matrix_free(M)
    call matrix_free(invS)
  end subroutine apply_fermi_function_dem


  !============================================================================!
  ! This subroutine applies the Fermi operator to a given Hamiltonian and      !
  ! overlap matrix, and returns the density matrix. The Fermi-level and the    !
  ! smearing temperature must be given. This form of the routine performs the  !
  ! operation on ONETEP dense matrices.                                        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The Hamiltonian matrix                           !
  !   overlap      (input)  : The overlap matrix                               !
  !   fermi_level  (input)  : The Fermi level                                  !
  !   smearing     (input)  : The smearing temerature                          !
  !   densmat      (output) : The computed density matrix                      !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Mar 2015.                                        !
  !============================================================================!
  subroutine apply_fermi_function_spam3(H,overlap,inv_overlap,fermi_level, &
       & smearing,densmat, correct_chemical_potential, Ne, boundl, boundu, &
       & fermi_tol, entropy, entropy_approx, lowdin_orthog)
    use comms, only : pub_on_root
    use sparse, only : sparse_copy, sparse_num_cols, sparse_num_rows, sparse_create
    use utils, only : utils_abort, utils_assert
    implicit none
    type(SPAM3),      intent(in)     :: H
    type(SPAM3),      intent(in)     :: overlap
    type(SPAM3),      intent(in)     :: inv_overlap
    real(kind=dp),    intent(inout)  :: fermi_level
    real(kind=dp),    intent(in)     :: smearing
    type(SPAM3),      intent(inout)  :: densmat
    logical,       optional,       intent(in)     :: correct_chemical_potential
    real(kind=dp), optional,       intent(in)     :: Ne
    real(kind=dp), optional,       intent(in)     :: boundl
    real(kind=dp), optional,       intent(in)     :: boundu
    real(kind=dp), optional,       intent(in)     :: fermi_tol
    real(kind=dp), optional,       intent(out)    :: entropy
    integer,       optional,       intent(in)     :: entropy_approx
    logical,       optional,       intent(in)     :: lowdin_orthog
    integer       :: loc_entropy_approx
    integer :: N
    type(smearing_matrix) :: ham
    type(smearing_matrix) :: rho
    type(smearing_matrix) :: M
    type(smearing_matrix) :: invS
    logical       :: loc_correct
    real(kind=dp) :: loc_ubound, loc_lbound, loc_Ne, loc_fermi_tol
    ! agrecocmplx
    logical :: loc_cmplx
    logical :: loc_lowdin_orthog

    ! agrecocmplx
    loc_cmplx = overlap%iscmplx


    loc_lowdin_orthog=.true.
    if(present(lowdin_orthog)) then
       loc_lowdin_orthog=lowdin_orthog
    end if

    if(present(correct_chemical_potential)) then
       loc_correct=correct_chemical_potential
       if(loc_correct) then
          call utils_assert(present(Ne), "Error in apply_fermi_function_spam3 : Must supply Ne to correct chemical potential!")
          loc_Ne=Ne

          if(present(boundl)) then
             loc_lbound=boundl
          else
             write(stdout,*) "smearing-> apply_fermi_function_spam3 : lbound not supplied, using -5.0"
             loc_lbound=-5.0_dp
          end if
          if(present(boundu)) then
             loc_ubound=boundu
          else
             write(stdout,*) "smearing-> apply_fermi_function_spam3 : ubound not supplied, using 5.0"
             loc_ubound=5.0_dp
          end if
          if(present(fermi_tol)) then
             loc_fermi_tol=fermi_tol
          else
             write(stdout,*) "smearing-> apply_fermi_function_spam3 : fermi_tol not supplied, using 1e-7"
             loc_fermi_tol=1e-7_dp
          end if
       end if
    else
       loc_correct=.false.
    end if

    N = sparse_num_rows(H)
    if(N/=sparse_num_cols(H)) then
       call utils_abort('Error in apply_fermi_function_spam3: H matrix not square!')
    end if
    ham%matrix_type=matrix_type_spam3
    ! agrecocmplx
    ham%standard_is_cmplx = loc_cmplx
    !    call matrix_allocate(ham,N,N)
    call sparse_create(ham%dataSPAM3,H)
    ham%matrix_is_allocated=.true.

    call sparse_copy(ham%dataSPAM3,H)

    !call matrix_allocate(rho,ham)
    rho%matrix_type=matrix_type_spam3
    rho%standard_is_cmplx = loc_cmplx
    call sparse_create(rho%dataSPAM3,densmat)
    rho%matrix_is_allocated=.true.


    !call matrix_allocate(M,ham)
    M%matrix_type=matrix_type_spam3
    M%standard_is_cmplx = loc_cmplx
    call sparse_create(M%dataSPAM3,overlap)
    M%matrix_is_allocated=.true.

    if(N/=sparse_num_rows(M%dataSPAM3)) then
       call utils_abort('Error in apply_fermi_function_spam3: H and overlap dims not consistent!')
    end if
    if(N/=sparse_num_cols(M%dataSPAM3)) then
       call utils_abort('Error in apply_fermi_function_spam3: overlap matrix not square!')
    end if

    call sparse_copy(M%dataSPAM3,overlap)

    invS%matrix_type=matrix_type_spam3
    invS%standard_is_cmplx = loc_cmplx
    call sparse_create(invS%dataSPAM3,inv_overlap)
    invS%matrix_is_allocated=.true.
    call sparse_copy(invS%dataSPAM3,inv_overlap)


    if(present(entropy_approx)) then
       loc_entropy_approx=entropy_approx
    else
       loc_entropy_approx=0
    end if

    if(present(entropy)) then
       call apply_fermi_function_common(ham,M,invS,fermi_level,smearing,rho,loc_correct,&
            &loc_Ne,loc_lbound,loc_ubound,loc_fermi_tol,entropy,loc_entropy_approx,&
            &lowdin_on=loc_lowdin_orthog)
    else
       call apply_fermi_function_common(ham,M,invS,fermi_level,smearing,rho,loc_correct,&
            &loc_Ne,loc_lbound,loc_ubound,loc_fermi_tol,lowdin_on=loc_lowdin_orthog)
    end if


    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "enter sparse copy", "<== SMEAR"
    end if
    call sparse_copy(densmat,rho%dataSPAM3)
    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "leave sparse copy", "<== SMEAR"
    end if

    call matrix_free(rho)
    call matrix_free(ham)
    call matrix_free(M)
    call matrix_free(invS)
  end subroutine apply_fermi_function_spam3

  !============================================================================!
  ! This subroutine applies the Fermi operator to a given Hamiltonian and      !
  ! overlap matrix, and returns the density matrix. The Fermi-level and the    !
  ! smearing temperature must be given. This form of the routine performs the  !
  ! operation on ONETEP SPAM3_EMBED matrices.                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The Hamiltonian matrix                           !
  !   overlap      (input)  : The overlap matrix                               !
  !   fermi_level  (input)  : The Fermi level                                  !
  !   smearing     (input)  : The smearing temerature                          !
  !   densmat      (output) : The computed density matrix                      !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, August 2018.                                   !
  !============================================================================!
  subroutine apply_fermi_function_spam3_embed(H,overlap,inv_overlap,fermi_level,&
       & smearing,densmat, correct_chemical_potential, Ne, boundl, boundu, &
       & fermi_tol, entropy, entropy_approx, lowdin_orthog)
    use comms, only : pub_on_root
    use sparse_embed, only : sparse_embed_copy, sparse_embed_num_cols, &
         sparse_embed_num_rows, sparse_embed_create
    use utils, only : utils_abort, utils_assert
    implicit none
    type(SPAM3_EMBED),      intent(in)     :: H
    type(SPAM3_EMBED),      intent(in)     :: overlap
    type(SPAM3_EMBED),      intent(in)     :: inv_overlap
    real(kind=dp),    intent(inout)  :: fermi_level
    real(kind=dp),    intent(in)     :: smearing
    type(SPAM3_EMBED),      intent(inout)  :: densmat
    logical,       optional,       intent(in)     :: correct_chemical_potential
    real(kind=dp), optional,       intent(in)     :: Ne
    real(kind=dp), optional,       intent(in)     :: boundl
    real(kind=dp), optional,       intent(in)     :: boundu
    real(kind=dp), optional,       intent(in)     :: fermi_tol
    real(kind=dp), optional,       intent(out)    :: entropy
    integer,       optional,       intent(in)     :: entropy_approx
    logical,       optional,       intent(in)     :: lowdin_orthog
    integer       :: loc_entropy_approx
    integer :: N
    type(smearing_matrix) :: ham
    type(smearing_matrix) :: rho
    type(smearing_matrix) :: M
    type(smearing_matrix) :: invS
    logical       :: loc_correct
    real(kind=dp) :: loc_ubound, loc_lbound, loc_Ne, loc_fermi_tol
    ! agrecocmplx
    logical :: loc_cmplx
    logical :: loc_lowdin_orthog

    ! agrecocmplx
    loc_cmplx = overlap%iscmplx


    loc_lowdin_orthog=.true.
    if(present(lowdin_orthog)) then
       loc_lowdin_orthog=lowdin_orthog
    end if

    if(present(correct_chemical_potential)) then
       loc_correct=correct_chemical_potential
       if(loc_correct) then
          call utils_assert(present(Ne), "Error in apply_fermi_function_spam3 : Must supply Ne to correct chemical potential!")
          loc_Ne=Ne

          if(present(boundl)) then
             loc_lbound=boundl
          else
             write(stdout,*) "smearing-> apply_fermi_function_spam3 : lbound not supplied, using -5.0"
             loc_lbound=-5.0_dp
          end if
          if(present(boundu)) then
             loc_ubound=boundu
          else
             write(stdout,*) "smearing-> apply_fermi_function_spam3 : ubound not supplied, using 5.0"
             loc_ubound=5.0_dp
          end if
          if(present(fermi_tol)) then
             loc_fermi_tol=fermi_tol
          else
             write(stdout,*) "smearing-> apply_fermi_function_spam3 : fermi_tol not supplied, using 1e-7"
             loc_fermi_tol=1e-7_dp
          end if
       end if
    else
       loc_correct=.false.
    end if

    N = sparse_embed_num_rows(H)
    if(N/=sparse_embed_num_cols(H)) then
       call utils_abort('Error in apply_fermi_function_spam3: H matrix not square!')
    end if
    ham%matrix_type=matrix_type_spam3_embed
    ! agrecocmplx
    ham%standard_is_cmplx = loc_cmplx
    !    call matrix_allocate(ham,N,N)
    call sparse_embed_create(ham%dataSPAM3_EMBED,H)
    ham%matrix_is_allocated=.true.

    call sparse_embed_copy(ham%dataSPAM3_EMBED,H)

    !call matrix_allocate(rho,ham)
    rho%matrix_type=matrix_type_spam3_embed
    rho%standard_is_cmplx = loc_cmplx
    call sparse_embed_create(rho%dataSPAM3_EMBED,densmat)
    rho%matrix_is_allocated=.true.


    !call matrix_allocate(M,ham)
    M%matrix_type=matrix_type_spam3_embed
    M%standard_is_cmplx = loc_cmplx
    call sparse_embed_create(M%dataSPAM3_EMBED,overlap)
    M%matrix_is_allocated=.true.

    if(N/=sparse_embed_num_rows(M%dataSPAM3_EMBED)) then
       call utils_abort('Error in apply_fermi_function_spam3: H and overlap dims not consistent!')
    end if
    if(N/=sparse_embed_num_cols(M%dataSPAM3_EMBED)) then
       call utils_abort('Error in apply_fermi_function_spam3: overlap matrix not square!')
    end if

    call sparse_embed_copy(M%dataSPAM3_EMBED,overlap)

    invS%matrix_type=matrix_type_spam3_embed
    invS%standard_is_cmplx = loc_cmplx
    call sparse_embed_create(invS%dataSPAM3_EMBED,inv_overlap)
    invS%matrix_is_allocated=.true.
    call sparse_embed_copy(invS%dataSPAM3_EMBED,inv_overlap)

    if(present(entropy_approx)) then
       loc_entropy_approx=entropy_approx
    else
       loc_entropy_approx=0
    end if

    if(present(entropy)) then
       call apply_fermi_function_common(ham,M,invS,fermi_level,smearing,rho,loc_correct,&
            &loc_Ne,loc_lbound,loc_ubound,loc_fermi_tol,entropy,loc_entropy_approx,&
            &lowdin_on=loc_lowdin_orthog)
    else
       call apply_fermi_function_common(ham,M,invS,fermi_level,smearing,rho,loc_correct,&
            &loc_Ne,loc_lbound,loc_ubound,loc_fermi_tol,lowdin_on=loc_lowdin_orthog)
    end if


    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "enter sparse copy", "<== SMEAR"
    end if
    call sparse_embed_copy(densmat,rho%dataSPAM3_EMBED)
    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "leave sparse copy", "<== SMEAR"
    end if

    call matrix_free(rho)
    call matrix_free(ham)
    call matrix_free(M)
    call matrix_free(invS)
  end subroutine apply_fermi_function_spam3_embed


  !============================================================================!
  ! This subroutine applies the Fermi operator to a given Hamiltonian and      !
  ! overlap matrix, and returns the density matrix. The Fermi-level and the    !
  ! smearing temperature must be given. This form of the routine performs the  !
  ! operation on the internal, abstracted matrix type.                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The Hamiltonian matrix                           !
  !   overlap      (input)  : The overlap matrix                               !
  !   fermi_level  (input)  : The Fermi level                                  !
  !   smearing     (input)  : The smearing temerature                          !
  !   densmat      (output) : The computed density matrix                      !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Feb 2014.                                        !
  !----------------------------------------------------------------------------!
  ! Seperated common operations into this routine, Feb 2014.                   !
  !============================================================================!
  subroutine apply_fermi_function_common(H,overlap,inv_overlap,fermi_level, &
       & smearing,densmat,correct_chemical_potential, Ne, boundl, boundu, &
       & fermi_tol,entropy, entropy_approx, lowdin_on)
    use comms, only: pub_on_root
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_assert
    implicit none
    type(smearing_matrix), intent(inout)  :: H
    type(smearing_matrix), intent(inout)  :: overlap
    type(smearing_matrix), intent(inout)  :: inv_overlap
    real(kind=dp),           intent(inout)   :: fermi_level
    real(kind=dp),intent(in)     :: smearing
    type(smearing_matrix),            intent(inout)   :: densmat ! Inout to pass allocation data
    logical,       optional, intent(in)      :: correct_chemical_potential
    real(kind=dp), optional, intent(in)      :: Ne
    real(kind=dp), optional, intent(in)      :: boundl
    real(kind=dp), optional, intent(in)      :: boundu
    real(kind=dp), optional, intent(in)      :: fermi_tol
    real(kind=dp), optional, intent(out)     :: entropy
    integer,       optional, intent(in)      :: entropy_approx
    logical,       optional, intent(in)      :: lowdin_on


    integer :: entropy_approx_loc

    type(smearing_matrix) :: inv_sqrt_S
    type(smearing_matrix) :: inv_sqrt_S_trans
    type(smearing_matrix) :: tmp, invH
!    type(smearing_matrix) :: invS
    type(smearing_matrix) :: CCH

    real(kind=dp) :: mineval,maxeval
    integer ::i
    integer :: matmuls=0

    real(kind=dp) :: fermi, rnum, trval
    logical :: loc_correct
    real(kind=dp) :: loc_ubound, loc_lbound, loc_fermi_tol
    real(kind=dp) :: loc_smear
    real(kind=dp) :: dmtrace


    logical, save :: have_invinit_a, have_invinit_b
    type(smearing_matrix), save :: invinit_a, invinit_b
    integer :: matmuls_changebeta

    integer :: proj_matmuls
    real(kind=dp),save :: threshold=10.0_dp

    real(kind=dp) :: max_eig
    logical :: staging
    integer :: stage_fac
    ! agrecocmplx
    logical :: loc_cmplx
    logical :: loc_lowdin_on
    logical :: Nan

    character(len=16) :: foe_type

    call timer_clock('apply_fermi_function_com',1)


    !ja531->sparsefoe:
    loc_lowdin_on=.false.
    if(present(lowdin_on)) then
       loc_lowdin_on=lowdin_on
    end if

    if(loc_lowdin_on) then
       foe_type=SIGNUM_PROJECTION
       foe_type=CERIOTTI_KUHNE_PA
       !       foe_type=CONTOUR_PROJECTIO
       foe_type=ANNEALING
    else
       foe_type=ANNEALING
!       foe_type=CONTOUR_PROJECTIO
    end if

    matmuls=0
    proj_matmuls=0

    ! agrecocmplx
    loc_cmplx = overlap%standard_is_cmplx

    if(present(correct_chemical_potential)) then
       loc_correct=correct_chemical_potential
       if(loc_correct) then
          call utils_assert(present(Ne), "Error in apply_fermi_function_common : Must supply Ne to correct chemical potential!")
          if(present(boundl)) then
             loc_lbound=boundl
          else
             write(stdout,*) "apply_fermi_function_common : lbound not supplied, using -5.0"
             loc_lbound=-5.0_dp
          end if
          if(present(boundu)) then
             loc_ubound=boundu
          else
             write(stdout,*) "apply_fermi_function_common : ubound not supplied, using 5.0"
             loc_ubound=5.0_dp
          end if
          if(present(fermi_tol)) then
             loc_fermi_tol=fermi_tol
          else
             write(stdout,*) "apply_fermi_function_common : fermi_tol not supplied, using 1e-7"
             loc_fermi_tol=1e-7_dp
          end if
       end if
    else
       loc_correct=.false.
    end if

    staging = .false.
    stage_fac = 4

!!!!!
    !    loc_smear=smearing/16.0_dp
    if(staging) then
       loc_smear=smearing*real(stage_fac,dp)
    else
       loc_smear=smearing
    end if
    !    loc_smear=smearing
!!!!!

    if(loc_lowdin_on) then
       ! ja531 -> sparsity fixes
       call matrix_allocate(inv_sqrt_S,densmat)

       ! ja531-> I think a Cholesky factorisation would be better, but unfortunately
       ! ja531-> we don't have a sparse routine.
       call lowdin_transformation(overlap,inv_sqrt_S)

       ! ja531 -> sparsity fixes
       call matrix_allocate(inv_sqrt_S_trans,densmat)
       call matrix_transpose(inv_sqrt_S,inv_sqrt_S_trans)


       ! ja531-> Apply Lowdin:
       call matrix_allocate(tmp,densmat)
       call matrix_scale(tmp,0.0_dp)
       call matrix_multiply(1.0_dp,inv_sqrt_S,H,0.0_dp,tmp)
       call matrix_multiply(1.0_dp,tmp,inv_sqrt_S,0.0_dp,H)
       call matrix_free(tmp)
    else ! natural representation
       !ja531-> assume H --> H^bar = S^-1 H
       call matrix_allocate(CCH, inv_overlap, H)

!       call matrix_info(CCH,"CCH")
!       call matrix_info(inv_overlap,"inv_overlap")
!       call matrix_info(H,"H")

       call matrix_multiply(1.0_dp,inv_overlap,H,0.0_dp,CCH)
    end if

    ! ja531-> put Hamiltonian matrix into units of smearing widths and move chemical potential to zero.
    if(.not.loc_lowdin_on) then
       call matrix_scale(H,1.0_dp/loc_smear)
       call matrix_axpy(H,overlap,-fermi_level/loc_smear)
       call matrix_scale(CCH,1.0_dp/loc_smear,-fermi_level/loc_smear)
    else
       call matrix_scale(H,1.0_dp/loc_smear,-fermi_level/loc_smear)
    end if

    ! ja531-> The problem with smear campaigns is that too often they work.
    ! ja531-> Here we construct the density kernel at the trial chemical potential

    !    call fermi_operator_chebyshev(H,densmat,matmuls)
!    call fermi_operator_projection(H,densmat,proj_matmuls)


    select case(FOE_TYPE)
    case(SIGNUM_PROJECTION)
       call utils_assert(loc_lowdin_on,"Error in Smearing_changemu: Sign proj method only supports orthogonal Hams.")
       call fermi_operator_projection(H,densmat,proj_matmuls)
    case(CONTOUR_PROJECTIO)
       if(loc_lowdin_on) then
          call fermi_operator_projection_contourint(H,H,H,H,densmat,orthogonal=.true.,&
               contravariant_kernel=.false.,avoid_solve=.true.,matmuls_out=proj_matmuls)
       else
          call fermi_operator_projection_contourint(H,overlap,inv_overlap,CCH,densmat,&
               orthogonal=.false.,contravariant_kernel=.true.,avoid_solve=.true.,matmuls_out=proj_matmuls)
       end if
    case(CERIOTTI_KUHNE_PA)
       call utils_assert(loc_lowdin_on,"Error in Smearing_changemu: Ceriotti method only supports orthogonal Hams.")

       call fermi_operator_chebyshev(H,densmat,proj_matmuls)

    case(ANNEALING)
       if(loc_lowdin_on) then
          call fermi_operator_annealing(H,overlap,H,densmat,.true.,.false.,proj_matmuls)
       else
          call fermi_operator_annealing(H,overlap,CCH,densmat,.false.,.false.,proj_matmuls)
       end if
       !       case("HEAD_GORDON")
       !          call utils_assert(orthogonal,"Error in Smearing_changemu: Head-Gordon method only supports orthogonal Hams")
       !          call fermi_operator_head_gordon(H,densmat,proj_matmuls)
       !  !       case("LIN")
       !  !          call fermi_operator_pex(H,densmat,proj_matmuls)
    case default
       call utils_abort("Error: Unrecognised FOE type in smearing_changemu")
    end select


    matmuls=matmuls+proj_matmuls
    ! ja531-> Now we check to see if the FOE succeeded... and try an ad-hoc scheme (until it works...) if not.
    do ! ja531-> ***ISSUE*** This never ends! Is it guaranteed to work?
       if(abs(matrix_trace(densmat))>matrix_dimension(densmat).or.matrix_any_isnan(densmat)) then
          ! Probably hit a pole (sign function directly on top on an eigenvalue...
          ! easiest thing we can do is to randomly perturb until it works)
          ! This would be worth someone thinking about...
          dmtrace=matrix_trace(densmat)
          if(pub_on_root) then
             if(debug_info_toggle) write(stdout,*) "smearing-> pole : ", dmtrace, fermi_level, Ne, "<== SMEAR"
          end if
          call random_number(rnum)
          rnum=(1.0_dp-rnum)*2.0_dp
          !          if(H%matrix_type/=matrix_type_standard) then
          !             call comms_bcast(pub_root_proc_id,rnum)
          !          end if
          ! We probably hit an eigenvalue with the shift... pole.


          if(.not.loc_lowdin_on) then
             call matrix_scale(H,1.0_dp,fermi_level/loc_smear)
             !          fermi_level=fermi_level+0.00001_dp*fermi_level
             fermi_level=fermi_level+(rnum*fermi_level)
             call matrix_scale(H,1.0_dp,-fermi_level/loc_smear)
             call matrix_axpy(CCH,overlap,fermi_level/loc_smear)
             fermi_level=fermi_level+(rnum*fermi_level)
             call matrix_axpy(CCH,overlap,-fermi_level/loc_smear)
          else
             call matrix_scale(H,1.0_dp,fermi_level/loc_smear)
             !          fermi_level=fermi_level+0.00001_dp*fermi_level
             fermi_level=fermi_level+(rnum*fermi_level)
             call matrix_scale(H,1.0_dp,-fermi_level/loc_smear)
          end if

          select case(FOE_TYPE)
          case(SIGNUM_PROJECTION)
             call utils_assert(loc_lowdin_on,"Error in Smearing_changemu: Sign proj method only supports orthogonal Hams.")
             call fermi_operator_projection(H,densmat,proj_matmuls)
          case(CONTOUR_PROJECTIO)
             if(loc_lowdin_on) then
                call fermi_operator_projection_contourint(H,H,H,H,densmat,orthogonal=.true.,&
                     contravariant_kernel=.false.,avoid_solve=.true.,matmuls_out=proj_matmuls)
             else
                call fermi_operator_projection_contourint(H,overlap,inv_overlap,CCH,densmat,&
                     orthogonal=.false.,contravariant_kernel=.true.,avoid_solve=.true.,matmuls_out=proj_matmuls)
             end if

          case(CERIOTTI_KUHNE_PA)
             call utils_assert(loc_lowdin_on,"Error in Smearing_changemu: Ceriotti method only supports orthogonal Hams.")

             call fermi_operator_chebyshev(H,densmat,proj_matmuls)

          case(ANNEALING)
!             call utils_assert(.not.loc_lowdin_on,"Error in Smearing_changemu: Annealing method only supports nonorthogonal Hams.")
!             call fermi_operator_annealing(H,overlap,CCH,densmat,.false.,.false.,proj_matmuls)
             if(loc_lowdin_on) then
                call fermi_operator_annealing(H,overlap,H,densmat,.true.,.false.,proj_matmuls)
             else
                call fermi_operator_annealing(H,overlap,CCH,densmat,.false.,.false.,proj_matmuls)
             end if
             !       case("HEAD_GORDON")
             !          call utils_assert(orthogonal,"Error in Smearing_changemu: Head-Gordon method only supports orthogonal Hams")
             !          call fermi_operator_head_gordon(H,densmat,proj_matmuls)
             !  !       case("LIN")
             !  !          call fermi_operator_pex(H,densmat,proj_matmuls)
          case default
             call utils_abort("Error: Unrecognised FOE type in smearing_changemu")
          end select

          matmuls=matmuls+proj_matmuls
       else
          exit
       end if
    end do

    if(pub_on_root.and.debug_info_toggle) then
       write(stdout,*) "smearing-> proj_matmuls, loc_smear, smearing, threshold : ", &
            & proj_matmuls, loc_smear, smearing, threshold, "<== SMEAR"
!       write(stdout,*) "JOLY->","fermi_level=",fermi_level,"smearing=",smearing,&
!            & "correct_chempot=",correct_chemical_potential,"Ne=", Ne,"boundl=", boundl,&
!            & "bound=",boundu, "fermi_tol=",fermi_tol,"entropy=",entropy,"entropy_approx=",entropy_approx
    end if


    !call changemu(T,mu,Dmu,smearing,orth_ham,invinit,have_invinit,matmuls)
    dmtrace=matrix_trace(densmat)

    nan=matrix_any_isnan(densmat)
    if(pub_on_root) then

       if(debug_info_toggle) write(stdout,*) "smearing-> matmuls = ", matmuls, fermi_level, &
            nan, Ne, dmtrace,&
            & "<== SMEAR"
    end if

    fermi=fermi_level


    ! ja531-> Here we optionally correct the chemical potential to conserve electron number.
    if(loc_correct) then
!       call findfermi_safe_newton(densmat,Ne,fermi,loc_smear,loc_fermi_tol,H, &
!            loc_lbound,loc_ubound,threshold,proj_matmuls,matmuls)

       if(.not.loc_lowdin_on) then
          call findfermi_safe_newton(H,overlap,inv_overlap,CCH,densmat,.false.,.false.,Ne,fermi,loc_smear,&
               loc_fermi_tol,loc_lbound,loc_ubound,foe_type,threshold,proj_matmuls,matmuls,.true.)
       else
          call findfermi_safe_newton(H,H,H,H,densmat,.true.,.true.,Ne,fermi,loc_smear,loc_fermi_tol,loc_lbound,loc_ubound,&
               &foe_type,threshold,proj_matmuls,matmuls,.true.)
       end if

       fermi_level=fermi
       if(pub_on_root) then
          if(debug_info_toggle) write(stdout,*) &
               "smearing-> apply_fermi_function_common -> new chemical potential  : ", fermi_level, "<== SMEARING"
          if(debug_info_toggle) write(stdout,*) &
               "smearing-> apply_fermi_function_common -> target electron number   : ", Ne, "<== SMEARING"
       end if
       if(.not.loc_lowdin_on) then
          trval=matrix_trace(overlap,densmat)
       else
          trval=matrix_trace(densmat)
       end if
       if(pub_on_root) then
          if(debug_info_toggle) write(stdout,*) &
               "smearing-> apply_fermi_function_common -> actual electron number   : ", trval, "<== SMEARING"
          if(debug_info_toggle) write(stdout,*) &
               "smearing-> apply_fermi_function_common -> electon number residual : ", Ne-trval, "<== SMEARING"
       end if
    end if
    if(pub_on_root.and.debug_info_toggle) write(stdout,*) "smearing-> matmuls = ", matmuls, "<== SMEAR"




    if(staging) then
       call matrix_scale(densmat,2.0_dp,-1.0_dp)
       if(.not.invinit_a%matrix_is_allocated) then
          call matrix_allocate(invinit_a,densmat)
          have_invinit_a=.false.
       end if
       matmuls_changebeta=0
       call changebeta(densmat,stage_fac,invinit_a,have_invinit_a,matmuls_changebeta)
       matmuls=matmuls+matmuls_changebeta

       call matrix_scale(densmat,0.5_dp,0.5_dp)
!!!!!

       if(pub_on_root.and.debug_info_toggle) write(stdout,*) "smearing-> matmuls_changebeta = ", matmuls_changebeta, "<== SMEAR"

       if(pub_on_root.and.debug_info_toggle) write(stdout,*) &
            "smearing-> proj_matmuls, smearing, threshold : ", proj_matmuls, smearing, threshold, "<== SMEAR"

       !call changemu(T,mu,Dmu,smearing,orth_ham,invinit,have_invinit,matmuls)
       dmtrace=matrix_trace(densmat)
       nan=matrix_any_isnan(densmat)
       if(pub_on_root) then

          if(debug_info_toggle) write(stdout,*) "smearing-> matmuls = ", matmuls, fermi_level, &
               nan, Ne, dmtrace,&
               & "<== SMEAR"

       end if
       fermi=fermi_level

       if(loc_correct) then
!          call findfermi_safe_newton(densmat,Ne,fermi,smearing,loc_fermi_tol,H, &
!               loc_lbound,loc_ubound,threshold,proj_matmuls,matmuls)
          call findfermi_safe_newton(H,H,H,H,densmat,.true.,.true.,Ne,fermi,smearing,loc_fermi_tol,loc_lbound,loc_ubound,&
               &foe_type,threshold,proj_matmuls,matmuls,.true.)
          fermi_level=fermi
          if(pub_on_root) then
             if(debug_info_toggle) write(stdout,*) &
                  "smearing-> apply_fermi_function_common -> new chemical potential  : ", fermi_level, "<== SMEARING"
             if(debug_info_toggle) write(stdout,*) &
                  "smearing-> apply_fermi_function_common -> target electron number   : ", Ne, "<== SMEARING"
          end if
          trval=matrix_trace(densmat)
          if(pub_on_root) then
             if(debug_info_toggle) write(stdout,*) &
                  "smearing-> apply_fermi_function_common -> actual elctron number   : ", trval, "<== SMEARING"
             if(debug_info_toggle) write(stdout,*) &
                  "smearing-> apply_fermi_function_common -> electron number residual : ", Ne-trval, "<== SMEARING"
          end if
       end if
       if(pub_on_root.and.debug_info_toggle) write(stdout,*) "smearing-> matmuls = ", matmuls, "<== SMEAR"

    end if

    !write(6,*)
    !flush(6)
    if(present(entropy)) then
       if(present(entropy_approx)) then
          entropy_approx_loc=entropy_approx
       else
          entropy_approx_loc=0
       end if
       select case(entropy_approx_loc)
       case(0) ! Exact entropy
          entropy = calculate_fermi_entropy_common(densmat)
       case(1) ! Coarse approximation to entropic term
          entropy = fermi_entropy_approx(densmat)
       case(2)
          entropy = fermi_entropy_approx(densmat,refine=.true.)
       case default
          call utils_abort('This entropy approximation is not implemented!')
       end select
    end if


    ! Un-Lowdin:
    if(loc_lowdin_on) then
       call matrix_allocate(tmp,H)
       call matrix_multiply(1.0_dp,inv_sqrt_S,densmat,0.0_dp,tmp)
       call matrix_multiply(1.0_dp,tmp,inv_sqrt_S_trans,0.0_dp,densmat)
       call matrix_free(tmp)
       call matrix_free(inv_sqrt_S)
       call matrix_free(inv_sqrt_S_trans)
    else
       call matrix_free(CCH)
    end if


    call timer_clock('apply_fermi_function_com',2)

  end subroutine apply_fermi_function_common


  !============================================================================!
  ! This subroutine applies the Fermi operator to a given Hamiltonian and      !
  ! overlap matrix, and returns the density matrix. The Fermi-level and the    !
  ! smearing temperature must be given. This form of the routine performs the  !
  ! operation on SPAM3 matrices, and where no orthogonalisation is performed   !
  ! because the input Hamiltonian matrix is also expected in mixed             !
  ! contra-covariant form as the third argument.
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The Hamiltonian matrix                           !
  !   overlap      (input)  : The overlap matrix                               !
  !   ccham        (input)  : The contra-covariant Hamiltonian matrix          !
  !   chempot      (input)  : The chemical potential                           !
  !   smearing     (input)  : The smearing temerature                          !
  !   Ne_target    (input)  : The target number of electrons                   !
  !   entropy      (output) : The entropy of the system calculated with the    !
  !                           contra-covariant form of the density kernel,     !
  !                           which is freed at the end of this routine.       !
  !   denskern     (output) : The computed density matrix                      !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, May 2017.                                        !
  !============================================================================!
  subroutine apply_fermi_function_CCH_spam3(ham,overlap,ccham,chempot,smearing,Ne_target,entropy,denskern,inverse_overlap)
    use comms, only : pub_on_root
    use sparse, only : sparse_copy, sparse_create
    use rundat, only : pub_foe_mu_tol

    implicit none
    type(SPAM3),      intent(in)     :: ham
    type(SPAM3),      intent(in)     :: overlap
    type(SPAM3),      intent(inout)  :: ccham
    real(kind=dp),    intent(inout)  :: chempot
    real(kind=dp),    intent(in)     :: smearing
    integer,          intent(in)     :: Ne_target
    real(kind=dp),    intent(out)    :: entropy
    type(SPAM3),      intent(inout)  :: denskern
    type(SPAM3), optional, intent(in)     :: inverse_overlap

    logical :: orthogonal = .false.
    integer :: matmuls_tmp, matmuls_proj
    integer :: matmuls
    logical :: loc_cmplx
    integer :: outunit
    logical :: loc_correct

    logical :: use_nsh
    real(kind=dp) :: loc_fermi_tol, loc_lbound, loc_ubound
    real(kind=dp) :: trval
    real(kind=dp) :: threshold

    type(smearing_matrix) :: H
    type(smearing_matrix) :: S
    type(smearing_matrix) :: CCH
    type(smearing_matrix) :: K
    type(smearing_matrix) :: KS
    type(smearing_matrix) :: invS

    integer :: N=8
    real(kind=dp) :: x(8), w(8)

    real(kind=dp) :: Ne

    character(len=16) :: foe_type

    Ne=real(Ne_target,dp)

    matmuls=0

    ! ja531-> loc_correct specifies whether we are going to correct the chemical potential of the system
    ! to have the target number of electrons. This is assumed to be true for now...
    loc_correct=.true.
!    loc_fermi_tol=1e-7_dp
    loc_fermi_tol=pub_foe_mu_tol

    ! ja531-> these bounds on the chemical potential search are arbitrary, but probably big enough that
    ! there will never be a problem.
    loc_lbound=-10.0_dp
    loc_ubound=10.0_dp
    ! ja531-> this tells the routine which changes the chemical potential not to use a Newton-Schulz-Hotelling algorithm.
    use_nsh=.false. ! This is solely used in changemu
    threshold=0.500000_dp

    loc_cmplx = overlap%iscmplx

    foe_type=annealing

    ! ja531-> Allocate memory for local matrices

    H%matrix_type=matrix_type_spam3
    H%standard_is_cmplx = loc_cmplx
    call sparse_create(H%dataSPAM3,ham)
    H%matrix_is_allocated=.true.

    S%matrix_type=matrix_type_spam3
    S%standard_is_cmplx = loc_cmplx
    call sparse_create(S%dataSPAM3,overlap)
    S%matrix_is_allocated=.true.

    CCH%matrix_type=matrix_type_spam3
    CCH%standard_is_cmplx = loc_cmplx
    call sparse_create(CCH%dataSPAM3,ccham)
    CCH%matrix_is_allocated=.true.

    K%matrix_type=matrix_type_spam3
    K%standard_is_cmplx = loc_cmplx
    call sparse_create(K%dataSPAM3,denskern)
    K%matrix_is_allocated=.true.

    if(present(inverse_overlap)) then
       invS%matrix_type=matrix_type_spam3
       invS%standard_is_cmplx = loc_cmplx
       call sparse_create(invS%dataSPAM3,inverse_overlap)
       invS%matrix_is_allocated=.true.
       call sparse_copy(invS%dataSPAM3,inverse_overlap)
    end if

    call sparse_copy(H%dataSPAM3,ham)
    call sparse_copy(CCH%dataSPAM3,ccham)
    call sparse_copy(S%dataSPAM3,overlap)


!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="H.mtx")
!    end if
!    call sparse_show_matrix(H%dataSPAM3,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
    !    end if


    ! ja531-> Shift the Hamiltonian so that the eigenvalue origin lies on the chemical potential

    call matrix_axpy(H,S,-chempot)
    call matrix_scale(H,1.0_dp/smearing)

    call matrix_scale(CCH,1.0_dp,-chempot)
    call matrix_scale(CCH,1.0_dp/smearing)

    ! ja531-> at this point we need to perform the operator expansion on the shifted Hamiltonian matrix
    ! we have a large number of choices as to which algorithm to use, but for now we hard-code the
    ! annealing and quenching algorithm.

!    call fermi_operator_projection_contourint(H,S,CCH,K,orthogonal,.true.,matmuls_tmp)

    if(present(inverse_overlap)) then
       call fermi_operator_annealing(H,S,CCH,K,orthogonal,.true.,matmuls_tmp,invS)
    else
       call fermi_operator_annealing(H,S,CCH,K,orthogonal,.true.,matmuls_tmp)
    end if
    matmuls=matmuls+matmuls_tmp

    ! ja531-> at this point we try to correct the chemical potential to give the target number of electrons
    ! given in the routine argument Ne. This is possibly not something we always want to do, so this
    ! lives in a conditional and "loc_correct" may migrate to being an argument at some point if this is
    ! important.

    if(loc_correct) then

       ! ja531-> The only important thing to note about the findfermi routine beyond the normal interface
       ! is that if we are dealing with non-orthogonal matrices and we have a contra-covariant Hamiltonian
       ! matrix then we need to ensure that use_nsh is false.

       if(present(inverse_overlap)) then
          call findfermi_safe_newton(H,S,invS,CCH,K,orthogonal,use_nsh,Ne,chempot,smearing,loc_fermi_tol,loc_lbound,loc_ubound,&
               &foe_type,threshold,matmuls_proj,matmuls_tmp,avoid_solve=.true.)
       else
          call findfermi_safe_newton(H,S,invS,CCH,K,orthogonal,use_nsh,Ne,chempot,smearing,loc_fermi_tol,loc_lbound,loc_ubound,&
               &foe_type,threshold,matmuls_proj,matmuls_tmp,avoid_solve=.false.)
       end if
!       call findfermi_safe_newton_anneal(H,S,CCH,K,orthogonal,use_nsh,Ne,chempot,smearing,loc_fermi_tol,loc_lbound,loc_ubound,&
!            &threshold,matmuls_proj,matmuls_tmp)

       ! ja531-> This is some info which may be useful for debugging, but is not always printed as it is assumed that the
       ! chemical potential search gets close enough to the right answer that this is all fine.

       matmuls=matmuls+matmuls_tmp
        if(pub_on_root) then
           if(debug_info_toggle) write(stdout,*) &
                "smearing-> apply_fermi_function_common -> new chemical potential  : ", chempot, "<== SMEARING"
           if(debug_info_toggle) write(stdout,*) &
                "smearing-> apply_fermi_function_common -> target electron number   : ", Ne, "<== SMEARING"
        end if
        trval=matrix_trace(K,S)
        if(pub_on_root) then
           if(debug_info_toggle) write(stdout,*) &
                "smearing-> apply_fermi_function_common -> actual elctron number   : ", trval, "<== SMEARING"
           if(debug_info_toggle) write(stdout,*) &
                "smearing-> apply_fermi_function_common -> electron number residual : ", Ne-trval, "<== SMEARING"
        end if
     end if
     if(pub_on_root.and.debug_info_toggle) write(stdout,*) "smearing-> matmuls = ", matmuls, "<== SMEAR"

     ! ja531-> at this point we should have a contravariant density kernel produced by one of the Fermi operator expansions.
     ! We want to use this to get a value for the entropy, we we can do as an expansion, just like the FOE, if we
     ! have access to a contra-covariant density kernel. So this is set up as:
     call matrix_allocate(KS,K)
     call matrix_multiply(1.0_dp,S,K,0.0_dp,KS)

     ! ja531-> and then fed into a routine which computes the entropy as an approximate expansion of the Fermi-Dirac
     ! entropy expression
     entropy=fermi_entropy_approx(KS,S=S,refine=.false.,cocontra_K=.true.)
     if(pub_on_root.and.debug_info_toggle) write(stdout,*) &
                "smearing-> entropy = ", entropy, "<== SMEAR"
     call matrix_free(KS)

     ! ja531-> put the contravariant density kernel in the output argument.
    call sparse_copy(denskern,K%dataSPAM3)

    call matrix_free(H)
    call matrix_free(S)
    call matrix_free(CCH)
    call matrix_free(K)
    if(present(inverse_overlap)) then
       call matrix_free(invS)
    end if

  end subroutine apply_fermi_function_CCH_spam3



  !============================================================================!
  ! This subroutine applies the Fermi operator to a given Hamiltonian and      !
  ! overlap matrix, and returns the density matrix. The Fermi-level and the    !
  ! smearing temperature must be given. This form of the routine performs the  !
  ! operation on SPAM3_EMBED matrices, and where no orthogonalisation is       !
  ! performed because the input Hamiltonian matrix is also expected in mixed   !
  ! contra-covariant form as the third argument.                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The Hamiltonian matrix                           !
  !   overlap      (input)  : The overlap matrix                               !
  !   ccham        (input)  : The contra-covariant Hamiltonian matrix          !
  !   chempot      (input)  : The chemical potential                           !
  !   smearing     (input)  : The smearing temerature                          !
  !   Ne_target    (input)  : The target number of electrons                   !
  !   entropy      (output) : The entropy of the system calculated with the    !
  !                           contra-covariant form of the density kernel,     !
  !                           which is freed at the end of this routine.       !
  !   denskern     (output) : The computed density matrix                      !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, August 2018.                                   !
  !============================================================================!
  subroutine apply_fermi_function_CCH_spam3_embed(ham,overlap,ccham,chempot,smearing,Ne_target,entropy,denskern,inverse_overlap)
    use comms, only : pub_on_root
    use sparse_embed, only: sparse_embed_copy, &
         sparse_embed_create
    use rundat, only : pub_foe_mu_tol

    implicit none
    type(SPAM3_EMBED),      intent(in)     :: ham
    type(SPAM3_EMBED),      intent(in)     :: overlap
    type(SPAM3_EMBED),      intent(inout)  :: ccham
    real(kind=dp),    intent(inout)  :: chempot
    real(kind=dp),    intent(in)     :: smearing
    integer,          intent(in)     :: Ne_target
    real(kind=dp),    intent(out)    :: entropy
    type(SPAM3_EMBED),      intent(inout)  :: denskern
    type(SPAM3_EMBED), optional, intent(in)     :: inverse_overlap

    logical :: orthogonal = .false.
    integer :: matmuls_tmp, matmuls_proj
    integer :: matmuls
    logical :: loc_cmplx
    integer :: outunit
    logical :: loc_correct

    logical :: use_nsh
    real(kind=dp) :: loc_fermi_tol, loc_lbound, loc_ubound
    real(kind=dp) :: trval
    real(kind=dp) :: threshold

    type(smearing_matrix) :: H
    type(smearing_matrix) :: S
    type(smearing_matrix) :: CCH
    type(smearing_matrix) :: K
    type(smearing_matrix) :: KS
    type(smearing_matrix) :: invS

    integer :: N=8
    real(kind=dp) :: x(8), w(8)

    real(kind=dp) :: Ne

    character(len=16) :: foe_type

    Ne=real(Ne_target,dp)

    matmuls=0

    ! ja531-> loc_correct specifies whether we are going to correct the chemical potential of the system
    ! to have the target number of electrons. This is assumed to be true for now...
    loc_correct=.true.
!    loc_fermi_tol=1e-7_dp
    loc_fermi_tol=pub_foe_mu_tol

    ! ja531-> these bounds on the chemical potential search are arbitrary, but probably big enough that
    ! there will never be a problem.
    loc_lbound=-10.0_dp
    loc_ubound=10.0_dp
    ! ja531-> this tells the routine which changes the chemical potential not to use a Newton-Schulz-Hotelling algorithm.
    use_nsh=.false. ! This is solely used in changemu
    threshold=0.500000_dp

    loc_cmplx = overlap%iscmplx

    foe_type=annealing

    ! ja531-> Allocate memory for local matrices

    H%matrix_type=matrix_type_spam3_embed
    H%standard_is_cmplx = loc_cmplx
    call sparse_embed_create(H%dataSPAM3_EMBED,ham)
    H%matrix_is_allocated=.true.

    S%matrix_type=matrix_type_spam3_embed
    S%standard_is_cmplx = loc_cmplx
    call sparse_embed_create(S%dataSPAM3_EMBED,overlap)
    S%matrix_is_allocated=.true.

    CCH%matrix_type=matrix_type_spam3_embed
    CCH%standard_is_cmplx = loc_cmplx
    call sparse_embed_create(CCH%dataSPAM3_EMBED,ccham)
    CCH%matrix_is_allocated=.true.

    K%matrix_type=matrix_type_spam3_embed
    K%standard_is_cmplx = loc_cmplx
    call sparse_embed_create(K%dataSPAM3_EMBED,denskern)
    K%matrix_is_allocated=.true.

    if(present(inverse_overlap)) then
       invS%matrix_type=matrix_type_spam3_embed
       invS%standard_is_cmplx = loc_cmplx
       call sparse_embed_create(invS%dataSPAM3_EMBED,inverse_overlap)
       invS%matrix_is_allocated=.true.
       call sparse_embed_copy(invS%dataSPAM3_EMBED,inverse_overlap)
    end if

    call sparse_embed_copy(H%dataSPAM3_EMBED,ham)
    call sparse_embed_copy(CCH%dataSPAM3_EMBED,ccham)
    call sparse_embed_copy(S%dataSPAM3_EMBED,overlap)


!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="H.mtx")
!    end if
!    call sparse_show_matrix(H%dataSPAM3_EMBED,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
    !    end if


    ! ja531-> Shift the Hamiltonian so that the eigenvalue origin lies on the chemical potential

    call matrix_axpy(H,S,-chempot)
    call matrix_scale(H,1.0_dp/smearing)

    call matrix_scale(CCH,1.0_dp,-chempot)
    call matrix_scale(CCH,1.0_dp/smearing)

    ! ja531-> at this point we need to perform the operator expansion on the shifted Hamiltonian matrix
    ! we have a large number of choices as to which algorithm to use, but for now we hard-code the
    ! annealing and quenching algorithm.

!    call fermi_operator_projection_contourint(H,S,CCH,K,orthogonal,.true.,matmuls_tmp)

    if(present(inverse_overlap)) then
       call fermi_operator_annealing(H,S,CCH,K,orthogonal,.true.,matmuls_tmp,invS)
    else
       call fermi_operator_annealing(H,S,CCH,K,orthogonal,.true.,matmuls_tmp)
    end if
    matmuls=matmuls+matmuls_tmp

    ! ja531-> at this point we try to correct the chemical potential to give the target number of electrons
    ! given in the routine argument Ne. This is possibly not something we always want to do, so this
    ! lives in a conditional and "loc_correct" may migrate to being an argument at some point if this is
    ! important.

    if(loc_correct) then

       ! ja531-> The only important thing to note about the findfermi routine beyond the normal interface
       ! is that if we are dealing with non-orthogonal matrices and we have a contra-covariant Hamiltonian
       ! matrix then we need to ensure that use_nsh is false.

       if(present(inverse_overlap)) then
          call findfermi_safe_newton(H,S,invS,CCH,K,orthogonal,use_nsh,Ne,chempot,smearing,loc_fermi_tol,loc_lbound,loc_ubound,&
               &foe_type,threshold,matmuls_proj,matmuls_tmp,avoid_solve=.true.)
       else
          call findfermi_safe_newton(H,S,invS,CCH,K,orthogonal,use_nsh,Ne,chempot,smearing,loc_fermi_tol,loc_lbound,loc_ubound,&
               &foe_type,threshold,matmuls_proj,matmuls_tmp,avoid_solve=.false.)
       end if
!       call findfermi_safe_newton_anneal(H,S,CCH,K,orthogonal,use_nsh,Ne,chempot,smearing,loc_fermi_tol,loc_lbound,loc_ubound,&
!            &threshold,matmuls_proj,matmuls_tmp)

       ! ja531-> This is some info which may be useful for debugging, but is not always printed as it is assumed that the
       ! chemical potential search gets close enough to the right answer that this is all fine.

       matmuls=matmuls+matmuls_tmp
        if(pub_on_root) then
           if(debug_info_toggle) write(stdout,*) &
                "smearing-> apply_fermi_function_common -> new chemical potential  : ", chempot, "<== SMEARING"
           if(debug_info_toggle) write(stdout,*) &
                "smearing-> apply_fermi_function_common -> target electron number   : ", Ne, "<== SMEARING"
        end if
        trval=matrix_trace(S,K)
        if(pub_on_root) then
           if(debug_info_toggle) write(stdout,*) &
                "smearing-> apply_fermi_function_common -> actual elctron number   : ", trval, "<== SMEARING"
           if(debug_info_toggle) write(stdout,*) &
                "smearing-> apply_fermi_function_common -> electron number residual : ", Ne-trval, "<== SMEARING"
        end if
     end if
     if(pub_on_root.and.debug_info_toggle) write(stdout,*) "smearing-> matmuls = ", matmuls, "<== SMEAR"

     ! ja531-> at this point we should have a contravariant density kernel produced by one of the Fermi operator expansions.
     ! We want to use this to get a value for the entropy, we we can do as an expansion, just like the FOE, if we
     ! have access to a contra-covariant density kernel. So this is set up as:
     call matrix_allocate(KS,K)
     call matrix_multiply(1.0_dp,S,K,0.0_dp,KS)

     ! ja531-> and then fed into a routine which computes the entropy as an approximate expansion of the Fermi-Dirac
     ! entropy expression
     entropy=fermi_entropy_approx(KS,S=S,refine=.false.,cocontra_K=.true.)
     if(pub_on_root.and.debug_info_toggle) write(stdout,*) &
                "smearing-> entropy = ", entropy, "<== SMEAR"
     call matrix_free(KS)

     ! ja531-> put the contravariant density kernel in the output argument.
    call sparse_embed_copy(denskern,K%dataSPAM3_EMBED)

    call matrix_free(H)
    call matrix_free(S)
    call matrix_free(CCH)
    call matrix_free(K)
    if(present(inverse_overlap)) then
       call matrix_free(invS)
    end if

  end subroutine apply_fermi_function_CCH_spam3_embed



  !============================================================================!
  ! The following routine calculates the entropy only as an expansion in the   !
  ! density kernel. This assumes that the density kernel has already been      !
  ! computed as a series, or otherwise and can be passed into these routines.  !
  ! This is the interface for fortran serial matrices.                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   denskern     (input)  : The density kernel                               !
  !   overlap      (input)  : The overlap matrix                               !
  !   S            (output) : The computed entropy                             !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Mar 2015.                                        !
  !============================================================================!
  subroutine calculate_fermi_entropy_standard(denskern,overlap,s)
    !    use comms, only:
    use utils, only: utils_abort
    implicit none
    real(kind=dp), dimension(:,:), intent(in)     :: denskern
    real(kind=dp), dimension(:,:), intent(in)     :: overlap
    real(kind=dp),                 intent(out)    :: s

    integer :: N
    type(smearing_matrix) :: Kmat
    type(smearing_matrix) :: Smat
    integer :: entropy_approx_loc

    N=size(denskern,1)
    if(N/=size(denskern,2)) then
       call utils_abort('Error in calculate_fermi_entropy_standard: denskern matrix not square!')
    end if
    if(N/=size(overlap,1)) then
       call utils_abort('Error in calculate_fermi_entropy_standard: overlap matrix not commensurate with denskern!')
    end if
    if(N/=size(overlap,2)) then
       call utils_abort('Error in calculate_fermi_entropy_standard: overlap matrix not square!')
    end if

    Kmat%matrix_type=matrix_type_standard
    call matrix_allocate(Kmat,N,N)
    Smat%matrix_type=matrix_type_standard
    call matrix_allocate(Smat,Kmat)

    Kmat%data = denskern
    Smat%data = overlap
    !
    !   if(present(entropy_approx)) then
    !      entropy_approx_loc=entropy_approx
    !   else
    !      entropy_approx_loc=0
    !   end if
    entropy_approx_loc=0
    select case(entropy_approx_loc)
    case(0) ! Exact entropy
       s = calculate_fermi_entropy_common(Kmat,Smat)
    case(1) ! Coarse approximation to entropic term
       s = fermi_entropy_approx(Kmat,Smat)
    case(2)
       s = fermi_entropy_approx(Kmat,Smat,refine=.true.)
    case default
       call utils_abort('This entropy approximation is not implemented!')
    end select

    !    s=calculate_fermi_entropy_common(Kmat,Smat)

    call matrix_free(Kmat)
    call matrix_free(Smat)

  end subroutine calculate_fermi_entropy_standard

  !============================================================================!
  ! The following routine calculates the entropy only as an expansion in the   !
  ! density kernel. This assumes that the density kernel has already been      !
  ! computed as a series, or otherwise and can be passed into these routines.  !
  ! This is the interface for DEM, dense matrices.                             !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   denskern     (input)  : The density kernel                               !
  !   overlap      (input)  : The overlap matrix                               !
  !   S            (output) : The computed entropy                             !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Mar 2015.                                        !
  !============================================================================!
  subroutine calculate_fermi_entropy_dem(denskern,overlap,s)
    use comms, only: pub_on_root
    use dense, only: dense_copy
    use utils, only: utils_abort
    implicit none
    type(DEM),        intent(in)     :: denskern
    type(DEM),        intent(in)     :: overlap
    real(kind=dp),    intent(out)    :: s

    integer :: N
    type(smearing_matrix) :: Kmat
    type(smearing_matrix) :: Smat

    integer :: entropy_approx_loc
    ! agrecocmplx
    logical :: loc_cmplx

    ! agrecocmplx
    loc_cmplx = overlap%iscmplx

    N=denskern%nrows
    if(N/=denskern%mcols) then
       call utils_abort('Error in calculate_fermi_entropy_dem: denskern matrix not square!')
    end if
    if(N/=overlap%nrows) then
       call utils_abort('Error in calculate_fermi_entropy_dem: overlap matrix not commensurate with denskern!')
    end if
    if(N/=overlap%mcols) then
       call utils_abort('Error in calculate_fermi_entropy_dem: overlap matrix not square!')
    end if

    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> Calculating entropy <== SMEARING"
    end if

    Kmat%matrix_type=matrix_type_dem
    ! agrecocmplx
    Kmat%standard_is_cmplx = loc_cmplx

    call matrix_allocate(Kmat,N,N)
    Smat%matrix_type=matrix_type_dem
    call matrix_allocate(Smat,Kmat)

    call dense_copy(Kmat%dataDEM,denskern)
    call dense_copy(Smat%dataDEM,overlap)


    entropy_approx_loc=1
    select case(entropy_approx_loc)
    case(0) ! Exact entropy
       s = calculate_fermi_entropy_common(Kmat,Smat)
    case(1) ! Coarse approximation to entropic term
       s = fermi_entropy_approx(Kmat,Smat)
    case(2)
       s = fermi_entropy_approx(Kmat,Smat,refine=.true.)
    case default
       call utils_abort('This entropy approximation is not implemented!')
    end select
    !    s=calculate_fermi_entropy_common(Kmat,Smat)

    call matrix_free(Kmat)
    call matrix_free(Smat)

  end subroutine calculate_fermi_entropy_dem

  !============================================================================!
  ! The following routine calculates the entropy only as an expansion in the   !
  ! density kernel. This assumes that the density kernel has already been      !
  ! computed as a series, or otherwise and can be passed into these routines.  !
  ! This is the interface for SPAM3 matrices.                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   denskern     (input)  : The density kernel                               !
  !   overlap      (input)  : The overlap matrix                               !
  !   S            (output) : The computed entropy                             !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Mar 2015.                                        !
  !============================================================================!
  subroutine calculate_fermi_entropy_spam3(denskern,overlap,s)
    use comms, only: pub_on_root
    use sparse, only: sparse_copy, sparse_num_rows, sparse_num_cols, sparse_create
    use utils, only: utils_abort
    implicit none
    type(SPAM3),      intent(in)     :: denskern
    type(SPAM3),      intent(in)     :: overlap
    real(kind=dp),    intent(out)    :: s

    integer :: N
    type(smearing_matrix) :: Kmat
    type(smearing_matrix) :: Smat

    integer :: entropy_approx_loc
    ! agrecocmplx
    logical :: loc_cmplx

    ! agrecocmplx
    loc_cmplx = overlap%iscmplx

    N = sparse_num_rows(denskern)
    if(N/=sparse_num_cols(denskern)) then
       call utils_abort('Error in calculate_fermi_entropy_spam3: denskern matrix not square!')
    end if
    if(N/=sparse_num_rows(overlap)) then
       call utils_abort('Error in calculate_fermi_entropy_spam3: overlap matrix not commensurate with denskern!')
    end if
    if(N/=sparse_num_cols(overlap)) then
       call utils_abort('Error in calculate_fermi_entropy_spam3: overlap matrix not square!')
    end if

    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> Calculating entropy <== SMEARING"
    end if

    Kmat%matrix_type=matrix_type_spam3
    ! agrecocmplx
    Kmat%standard_is_cmplx = loc_cmplx

    !    call matrix_allocate(Kmat,N,N)
    call sparse_create(Kmat%dataSPAM3,denskern)
    Kmat%matrix_is_allocated=.true.

    Smat%matrix_type=matrix_type_spam3
    call sparse_create(Smat%dataSPAM3,overlap)
    Smat%matrix_is_allocated=.true.
!    call matrix_allocate(Smat,Kmat)

    call sparse_copy(Kmat%dataSPAM3,denskern)
    call sparse_copy(Smat%dataSPAM3,overlap)

    entropy_approx_loc=0
    select case(entropy_approx_loc)
    case(0) ! Exact entropy
       s = calculate_fermi_entropy_common(Kmat,Smat)
    case(1) ! Coarse approximation to entropic term
       s = fermi_entropy_approx(Kmat,Smat)
    case(2)
       s = fermi_entropy_approx(Kmat,Smat,refine=.true.)
    case default
       call utils_abort('This entropy approximation is not implemented!')
    end select

    !    s=calculate_fermi_entropy_common(Kmat,Smat)

    call matrix_free(Kmat)
    call matrix_free(Smat)

  end subroutine calculate_fermi_entropy_spam3

  !============================================================================!
  ! The following routine calculates the entropy only as an expansion in the   !
  ! density kernel. This assumes that the density kernel has already been      !
  ! computed as a series, or otherwise and can be passed into these routines.  !
  ! This is the interface for SPAM3_EMBED matrices.                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   denskern     (input)  : The density kernel                               !
  !   overlap      (input)  : The overlap matrix                               !
  !   S            (output) : The computed entropy                             !
  !----------------------------------------------------------------------------!
  ! SPAM3_EMBED version added by Robert Charlton, 10th August 2018.            !
  ! Written by Jolyon Aarons, Mar 2015.                                        !
  !============================================================================!
  subroutine calculate_fermi_entropy_spam3_embed(denskern,overlap,s)
    use comms, only: pub_on_root
    use sparse_embed, only: sparse_embed_copy, sparse_embed_num_rows, &
         sparse_embed_num_cols, sparse_embed_create
    use utils, only: utils_abort
    implicit none
    type(SPAM3_EMBED),      intent(in)     :: denskern
    type(SPAM3_EMBED),      intent(in)     :: overlap
    real(kind=dp),          intent(out)    :: s

    integer :: N
    type(smearing_matrix) :: Kmat
    type(smearing_matrix) :: Smat

    integer :: entropy_approx_loc
    ! agrecocmplx
    logical :: loc_cmplx

    ! agrecocmplx
    loc_cmplx = overlap%iscmplx

    N = sparse_embed_num_rows(denskern)
    if(N/=sparse_embed_num_cols(denskern)) then
       call utils_abort('Error in calculate_fermi_entropy_spam3: denskern matrix not square!')
    end if
    if(N/=sparse_embed_num_rows(overlap)) then
       call utils_abort('Error in calculate_fermi_entropy_spam3: overlap matrix not commensurate with denskern!')
    end if
    if(N/=sparse_embed_num_cols(overlap)) then
       call utils_abort('Error in calculate_fermi_entropy_spam3: overlap matrix not square!')
    end if

    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> Calculating entropy <== SMEARING"
    end if

    Kmat%matrix_type=matrix_type_spam3_embed
    ! agrecocmplx
    Kmat%standard_is_cmplx = loc_cmplx

    !    call matrix_allocate(Kmat,N,N)
    call sparse_embed_create(Kmat%dataSPAM3_EMBED,denskern)
    Kmat%matrix_is_allocated=.true.

    Smat%matrix_type=matrix_type_spam3_embed
    call sparse_embed_create(Smat%dataSPAM3_EMBED,overlap)
    Smat%matrix_is_allocated=.true.
!    call matrix_allocate(Smat,Kmat)

    call sparse_embed_copy(Kmat%dataSPAM3_EMBED,denskern)
    call sparse_embed_copy(Smat%dataSPAM3_EMBED,overlap)

    entropy_approx_loc=0
    select case(entropy_approx_loc)
    case(0) ! Exact entropy
       s = calculate_fermi_entropy_common(Kmat,Smat)
    case(1) ! Coarse approximation to entropic term
       s = fermi_entropy_approx(Kmat,Smat)
    case(2)
       s = fermi_entropy_approx(Kmat,Smat,refine=.true.)
    case default
       call utils_abort('This entropy approximation is not implemented!')
    end select

    !    s=calculate_fermi_entropy_common(Kmat,Smat)

    call matrix_free(Kmat)
    call matrix_free(Smat)

  end subroutine calculate_fermi_entropy_spam3_embed

  !ep: mermin subroutine to calculate entropy and its derivative
  subroutine calculate_fermi_entropy_mermin(denskern,overlap, &
    entropy_approx_loc,s,expans,deriv,inv_overlap)
    use comms, only: pub_on_root
    use sparse_embed, only: sparse_embed_copy, sparse_embed_num_rows, &
         sparse_embed_num_cols, sparse_embed_create
    use utils, only: utils_abort
    implicit none
    type(SPAM3_EMBED),      intent(in)     :: denskern
    type(SPAM3_EMBED),      intent(in)     :: overlap
    integer, intent(in)                    :: entropy_approx_loc
    type(SPAM3_EMBED), optional, intent(inout)  :: deriv
    type(SPAM3_EMBED), optional, intent(in)  :: inv_overlap
    integer, optional, intent(in)  :: expans
    real(kind=dp),          intent(out)    :: s

    integer :: N
    type(smearing_matrix) :: Kmat
    type(smearing_matrix) :: Smat
    type(smearing_matrix) :: invSmat

    logical :: loc_cmplx

    ! agrecocmplx
    loc_cmplx = overlap%iscmplx

    N = sparse_embed_num_rows(denskern)
    if(N/=sparse_embed_num_cols(denskern)) then
       call utils_abort('Error in calculate_fermi_entropy_mermin: denskern &
           & matrix not square!')
    end if
    if(N/=sparse_embed_num_rows(overlap)) then
       call utils_abort('Error in calculate_fermi_entropy_mermin: overlap &
       & matrix not commensurate with denskern!')
    end if
    if(N/=sparse_embed_num_cols(overlap)) then
       call utils_abort('Error in calculate_fermi_entropy_mermin: overlap &
       & matrix not square!')
    end if

    if( (present(deriv) .and. entropy_approx_loc>3) .or. &
     ( entropy_approx_loc==2 .and. .not. present(deriv) ) .or. &
     ( entropy_approx_loc==3 .and. .not. present(deriv) )) then
       call utils_abort('Error: required entropy deriv but flag or &
     &  matrix not provided !')
    end if

    Kmat%matrix_type=matrix_type_spam3_embed
    Kmat%standard_is_cmplx = loc_cmplx

    call sparse_embed_create(Kmat%dataSPAM3_EMBED,denskern)
    Kmat%matrix_is_allocated=.true.

    Smat%matrix_type=matrix_type_spam3_embed
    call sparse_embed_create(Smat%dataSPAM3_EMBED,overlap)
    Smat%matrix_is_allocated=.true.

    call sparse_embed_copy(Kmat%dataSPAM3_EMBED,denskern)
    call sparse_embed_copy(Smat%dataSPAM3_EMBED,overlap)

    select case(entropy_approx_loc)
    case(0) ! Exact entropy (joly approx)
       s = calculate_fermi_entropy_common(Kmat,Smat)
    case(1)
      s = fermi_entropy_approx_mermin(Kmat,Smat,expans)
    case(2)
       call fermi_deriv_approx(Kmat,deriv,s,Smat,expans)
    case(3)
       invSmat%matrix_type=matrix_type_spam3_embed
       call sparse_embed_create(invSmat%dataSPAM3_EMBED,inv_overlap)
       invSmat%matrix_is_allocated=.true.
       call sparse_embed_copy(invSmat%dataSPAM3_EMBED,inv_overlap)
       call deriv_ngwf_approx(Kmat,deriv,Smat,expans,invSmat)
    case default
       call utils_abort('This entropy approximation is not implemented!')
    end select

    call matrix_free(Kmat)
    call matrix_free(Smat)

  end subroutine calculate_fermi_entropy_mermin
  !ep

  !============================================================================!
  ! The following routine calculates the entropy only as an expansion in the   !
  ! density kernel. This assumes that the density kernel has already been      !
  ! computed as a series, or otherwise and can be passed into these routines.  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   denskern     (input)  : The density kernel                               !
  !   overlap      (input)  : The overlap matrix                               !
  !   S            (output) : The computed entropy                             !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Mar 2015.                                        !
  !============================================================================!
  function calculate_fermi_entropy_common(denskern,overlap) result(s)
    use comms, only: comms_barrier
    use timer, only: timer_clock
    implicit none
    type(smearing_matrix),           intent(inout) :: denskern
    type(smearing_matrix), optional, intent(inout) :: overlap
    real(kind=dp)                         :: s
    type(smearing_matrix) :: K,L,C
    type(smearing_matrix) :: inv_sqrt_S, inv_sqrt_S_trans, tmp_mat
    integer :: lowdin_maxiters
    real(kind=dp) :: log_tol

    call timer_clock('calculate_fermi_entropy_com',1)

    lowdin_maxiters = 100
    log_tol = 1e-12_dp

    call matrix_allocate(K,denskern)
    call matrix_allocate(L,denskern)
    !ja531->sparsefoe:
    !call matrix_allocate(C,K,L)
    call matrix_allocate(C,denskern)

    ! Lowdin transform overlap matrix:
    if(present(overlap)) then
       !ja531->sparsefoe: surely inv_sqrt_S should have denskern sparsity? That's what I'm assuming now...
!       call matrix_allocate(inv_sqrt_S,overlap)
       call matrix_allocate(inv_sqrt_S,denskern)

       call lowdin_transformation(overlap,sqrt_S_out=inv_sqrt_S,maxiters=lowdin_maxiters)
       !ja531->sparsefoe: surely inv_sqrt_S should have denskern sparsity? That's what I'm assuming now...
!       call matrix_allocate(inv_sqrt_S_trans,overlap)
       call matrix_allocate(inv_sqrt_S_trans,denskern)

       call matrix_transpose(inv_sqrt_S,inv_sqrt_S_trans)
       !ja531->sparsefoe:
!       call matrix_allocate(tmp_mat,inv_sqrt_S_trans,denskern)
       call matrix_allocate(tmp_mat,denskern)

       call matrix_multiply(1.0_dp,inv_sqrt_S_trans,denskern,0.0_dp,tmp_mat)
       call matrix_free(inv_sqrt_S_trans)
       call matrix_multiply(1.0_dp,tmp_mat,inv_sqrt_S,0.0_dp,K)
       call matrix_free(inv_sqrt_S)
       call matrix_free(tmp_mat)
    else
       call matrix_copy(denskern,K)
    end if

    ! from here vvv
    call comms_barrier()
    !ja531 --> debug
    if(matrix_any_isnan(K).and.debug_info_toggle) write(stdout,*) "smearing-> EntNaN:1", "<== SMEAR"

    call matrix_safe_logarithm(K,L,lowdin_maxiters=lowdin_maxiters,mintol=log_tol)
    call comms_barrier()
    !ja531 --> debug
    if(matrix_any_isnan(L).and.debug_info_toggle) write(stdout,*) "smearing-> EntNaN:2", "<== SMEAR"

    call matrix_multiply(1.0_dp,K,L,0.0_dp,C)
    !ja531 --> debug

    if(matrix_any_isnan(C).and.debug_info_toggle) write(stdout,*) "smearing-> EntNaN:3", "<== SMEAR"
    call comms_barrier()


    call matrix_scale(K,-1.0_dp,1.0_dp)
    !ja531 --> debug
    if(matrix_any_isnan(K).and.debug_info_toggle) write(stdout,*) "smearing-> EntNaN:4", "<== SMEAR"
    call comms_barrier()

    call matrix_scale(L,0.0_dp)
    call matrix_safe_logarithm(K,L,lowdin_maxiters=lowdin_maxiters,mintol=log_tol)
    !ja531 --> debug
    if(matrix_any_isnan(L).and.debug_info_toggle) write(stdout,*) "smearing-> EntNaN:5", "<== SMEAR"
    call comms_barrier()
    call matrix_multiply(1.0_dp,K,L,1.0_dp,C)

    !ja531 --> debug
    if(matrix_any_isnan(C).and.debug_info_toggle) write(stdout,*) "smearing-> EntNaN:6", "<== SMEAR"
    call comms_barrier()
    s=matrix_trace(C)

    call matrix_free(K)
    call matrix_free(L)
    call matrix_free(C)

    call timer_clock('calculate_fermi_entropy_com',2)

  end function calculate_fermi_entropy_common


  !============================================================================!
  ! This routine returns the approximate eigenvalues of a matrix and their     !
  ! error radii in the complex plane (as we are dealing with real matrices,    !
  ! these are error bars).                                                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   A            (input)  : The input matrix                                 !
  !   eigs         (output) : The eigenvalue.                                  !
  !   radii        (output) : The eigenvector.                                 !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Jan 2014,                                        !
  !----------------------------------------------------------------------------!
  ! Gerschgorin often gives very poor eigenvalues, this is only useful for a   !
  ! maximum bound on the range of eigenvalues.                                 !
  !============================================================================!
  subroutine gerschgorin(A,eigs,radii)
    use dense, only : dense_get_element
    use utils, only : utils_abort
    implicit none
    type(smearing_matrix), intent(in) :: A
    real(kind=dp), dimension(:), intent(out) :: eigs
    real(kind=dp), dimension(:), intent(out) :: radii
    complex(kind=dp) :: zval
    real(kind=dp)    :: rval
    integer :: i,j,N
    N=matrix_dimension(A)
    select case(A%matrix_type)
    case(matrix_type_standard)
       do i=1,N
          eigs(i)=A%data(i,i)
          radii(i)=0.0_dp
          do j=1,i-1
             radii(i)=radii(i)+abs(A%data(j,i))
          end do
          do j=i+1,N
             radii(i)=radii(i)+abs(A%data(j,i))
          end do
       end do
    case(matrix_type_dem)
       if(A%dataDEM%iscmplx) then
          do i=1,N
             call dense_get_element(zval,A%dataDEM,i,i)
             eigs(i)=real(zval,dp)
             radii(i)=0.0_dp
             do j=1,i-1
                call dense_get_element(zval,A%dataDEM,j,i)
                radii(i)=radii(i)+abs(zval)
             end do
             do j=i+1,N
                call dense_get_element(zval,A%dataDEM,j,i)
                radii(i)=radii(i)+abs(zval)
             end do
          end do
       else
          do i=1,N
             call dense_get_element(rval,A%dataDEM,i,i)
             eigs(i)=rval
             radii(i)=0.0_dp
             do j=1,i-1
                call dense_get_element(rval,A%dataDEM,j,i)
                radii(i)=radii(i)+rval
             end do
             do j=i+1,N
                call dense_get_element(rval,A%dataDEM,j,i)
                radii(i)=radii(i)+rval
             end do
          end do
       end if
    case default
       call utils_abort('gerschgorin is not implemented for this type of matrix.')
    end select
  end subroutine gerschgorin

  !============================================================================!
  ! This routine computes the matrix exponential of... a matrix. This is based !
  ! completely upon Nicholas Higham, of Manchester University Maths School's   !
  ! work. In particular his 2008 paper on using order 13 Pade approximants.    !
  ! This implementation firstly decides whether we can get away with doing a   !
  ! lower order approximant, from 3 to 9. If not then it performs scaling and  !
  ! squaring until the order 13 PA is good enough.                             !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   A            (input)  : The input matrix                                 !
  !   R            (output) : The computed matrix exponential                  !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Jan 2014,                                        !
  !     (based on a version written for a generalized langevin equation        !
  !      thermostat in a path integral molecular dynamics code I wrote)        !
  !----------------------------------------------------------------------------!
  ! Fixed NaN bug. Apparently everything has to be explicity zeroed. Who Knew? !
  !   Feb 2014.                                                                !
  ! Scaling and squaring is non-optimal in this version. Does one step of      !
  ! divide and conquer. Could be much better, but we shouldn't need to end up  !
  ! doing it if we have well conditioned matrices.                             !
  !============================================================================!
  subroutine mat_exp_higham(A,R)
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none
    type(smearing_matrix), intent(inout) :: A ! Inout for dense_create_copy (in)
    type(smearing_matrix), intent(inout) :: R
    type(smearing_matrix) :: U, V
    type(smearing_matrix) :: A2, A4, A6, pow_mat
    type(smearing_matrix), dimension(:), allocatable :: tmp_mat
    real(kind=dp) :: norm
    real(kind=dp), dimension(:), allocatable :: c, b
    integer :: i, j, l
    real(kind=dp) :: s
    integer :: sint
    integer :: ierr

    call timer_clock('mat_exp_higham',1)

    ! Don't need balance transformation --> only dealing with Hermitian matrices

    norm = matrix_norm(A,1)
    allocate(c(5),stat=ierr)
    call utils_alloc_check('mat_exp_higham','c',ierr)
    allocate(b(14),stat=ierr)
    call utils_alloc_check('mat_exp_higham','b',ierr)
    allocate(tmp_mat(2),stat=ierr)
    call utils_alloc_check('mat_exp_higham','tmp_mat',ierr)

    c=[0.015_dp, 0.25_dp, 0.95_dp, 2.1_dp, 5.37_dp]

    ! If the norm is small enough, use the Pade-Approximation (PA) directly
    if (norm <= c(4)) then

       do i=1,4
          if(norm<=c(i)) exit
       end do
       l=i

       ! Special magic and northern wizardry went into computing these Pade approximant
       ! coefficients (MATLAB).
       select case(l)
       case(1)
          b=[120.0_dp,60.0_dp,12.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
               & 0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
       case(2)
          b=[30240.0_dp,15120.0_dp,3360.0_dp,420.0_dp,30.0_dp,1.0_dp,0.0_dp,0.0_dp,&
               & 0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
       case(3)
          b=[17297280.0_dp,8648640.0_dp,1995840.0_dp,277200.0_dp,25200.0_dp,1512.0_dp,&
               & 56.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
       case(4)
          b=[17643225600.0_dp,8821612800.0_dp,2075673600.0_dp,302702400.0_dp,30270240.0_dp,&
               & 2162160.0_dp,110880.0_dp,3960.0_dp,90.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
       case default

       end select

       call matrix_allocate(A2,A)
       call matrix_multiply(1.0_dp,A,A,0.0_dp,A2)

       call matrix_allocate(U,A)
       call matrix_scale(U,0.0_dp)
       call matrix_allocate(V,A)
       call matrix_scale(V,0.0_dp)
       call matrix_scale(U,0.0_dp,b(2))
       call matrix_scale(V,0.0_dp,b(1))

       call matrix_allocate(tmp_mat(1),A)
       call matrix_copy(A2,tmp_mat(1))
       call matrix_axpy(U,tmp_mat(1),b(4))
       call matrix_axpy(V,tmp_mat(1),b(3))

       call matrix_allocate(tmp_mat(2),A)
       do i=2,l
          call matrix_multiply(1.0_dp,tmp_mat(1+mod(i,2)),A2,0.0_dp,tmp_mat(1+mod(i-1,2)))
          call matrix_axpy(U,tmp_mat(1+mod(i-1,2)),b((2*i)+2))
          call matrix_axpy(V,tmp_mat(1+mod(i-1,2)),b((2*i)+1))
       end do
       call matrix_free(A2)

       call matrix_multiply(1.0_dp,A,U,0.0_dp,tmp_mat(1))
       call matrix_free(U)
       call matrix_copy(V,tmp_mat(2))
       call matrix_axpy(V,tmp_mat(1),-1.0_dp)
       call matrix_axpy(tmp_mat(2),tmp_mat(1),1.0_dp)
       call matrix_free(tmp_mat(1))

       ! ja531-> A wise man once said "Never invert a matrix!".
       call matrix_invert(V)
       ! ja531-> A Cholesky factorization would be \emph{really} useful here.
       call matrix_multiply(1.0_dp,V,tmp_mat(2),0.0_dp,R)
       call matrix_free(V)
       call matrix_free(tmp_mat(2))

       ! Else, check if norm of A is small enough for m=13.
       ! If not, scale the matrix
    else

       s=log(norm/5.37_dp)
       call matrix_allocate(tmp_mat(1),A)
       call matrix_allocate(tmp_mat(2),A)
       call matrix_copy(A,tmp_mat(2))
       ! Scaling
       if (s > 0.0_dp) then
          sint = ceiling(s)
          call matrix_scale(tmp_mat(2),1.0_dp/(2.0_dp**sint))
       end if
       ! Calculate PA
       b=[64764752532480000.0_dp,32382376266240000.0_dp,7771770303897600.0_dp,1187353796428800.0_dp, &
            & 129060195264000.0_dp,10559470521600.0_dp,670442572800.0_dp,33522128640.0_dp, 1323241920.0_dp, &
            & 40840800.0_dp,960960.0_dp,16380.0_dp, 182.0_dp, 1.0_dp]

       call matrix_allocate(A2,A)
       call matrix_allocate(A4,A)
       call matrix_allocate(A6,A)

       call matrix_allocate(pow_mat,A)
       call matrix_copy(tmp_mat(2),pow_mat)

       call matrix_multiply(1.0_dp,tmp_mat(2),pow_mat,0.0_dp,A2)

       call matrix_copy(A2,pow_mat)
       call matrix_multiply(1.0_dp,A2,pow_mat,0.0_dp,A4)

       call matrix_multiply(1.0_dp,A2,A4,0.0_dp,A6)

       call matrix_free(pow_mat)


       ! Construct U
       call matrix_allocate(U,A)
       call matrix_scale(tmp_mat(1),0.0_dp)
       call matrix_axpy(tmp_mat(1),A6,b(14))
       call matrix_axpy(tmp_mat(1),A4,b(12))
       call matrix_axpy(tmp_mat(1),A2,b(10))
       call matrix_multiply(1.0_dp,A6,tmp_mat(1),0.0_dp,U)
       call matrix_axpy(U,A6,b(8))
       call matrix_axpy(U,A4,b(6))
       call matrix_axpy(U,A2,b(4))
       !       call matrix_shift(U,b(2))
       call matrix_scale(U,1.0_dp,b(2))
       call matrix_copy(U,tmp_mat(1))
       call matrix_multiply(1.0_dp,tmp_mat(2),tmp_mat(1),0.0_dp,U)

       ! Construct V
       call matrix_allocate(V,A)
       call matrix_scale(tmp_mat(1),0.0_dp)
       call matrix_axpy(tmp_mat(1),A6,b(13))
       call matrix_axpy(tmp_mat(1),A4,b(11))
       call matrix_axpy(tmp_mat(1),A2,b(9))
       call matrix_multiply(1.0_dp,A6,tmp_mat(1),0.0_dp,V)
       call matrix_axpy(V,A6,b(7))
       call matrix_axpy(V,A4,b(5))
       call matrix_axpy(V,A2,b(3))
       !      call matrix_shift(V,b(1))
       call matrix_scale(V,1.0_dp,b(1))



       ! Solve for r13
       call matrix_copy(V, tmp_mat(1))
       call matrix_axpy(tmp_mat(1), U, -1.0_dp)
       call matrix_axpy(V,U,1.0_dp)
       ! Eurrrgh. Again...
       call matrix_invert(tmp_mat(1))

       call matrix_free(A2)
       call matrix_free(A4)
       call matrix_free(A6)
       call matrix_free(U)

       if(s > 0.0_dp) then
          call matrix_multiply(1.0_dp,tmp_mat(1),V,0.0_dp,tmp_mat(2))
          do i=1,sint-1
             call matrix_multiply(1.0_dp,tmp_mat(1+mod(i,2)),tmp_mat(1+mod(i,2)),0.0_dp,&
                  &tmp_mat(1+mod(i-1,2)))
          end do
          call matrix_multiply(1.0_dp,tmp_mat(1+mod(sint,2)),tmp_mat(1+mod(sint,2)),0.0_dp,R)
       else
          call matrix_multiply(1.0_dp,tmp_mat(1),V,0.0_dp,R)
       end if
       call matrix_free(tmp_mat(1))
       call matrix_free(tmp_mat(2))
       call matrix_free(V)
    end if

    deallocate(b,stat=ierr)
    call utils_dealloc_check('mat_exp_higham','b',ierr)
    deallocate(c,stat=ierr)
    call utils_dealloc_check('mat_exp_higham','c',ierr)
    deallocate(tmp_mat,stat=ierr)
    call utils_dealloc_check('mat_exp_higham','tmp_mat',ierr)

    call timer_clock('mat_exp_higham',2)

  end subroutine mat_exp_higham


  !============================================================================!
  ! This routine computes the logarithm of a matrix using the algorithm of     !
  ! Higham et al. See:                                                         !
  !                                                                            !
  ! Higham, N.J., 2001. Evaluating Pad\'e approximants of the matrix logarithm.!
  ! SIAM Journal on Matrix Analysis and Applications, 22(4), pp.1126-1135.     !
  !                                                                            !
  ! It requires a lot of matrix products which means that it is expensive to   !
  ! run. I have avoided using it wherever possible, but there is still the     !
  ! option to use it for computing the entropy if the exact (non-series exp)   !
  ! of the Fermi-Dirac entropy expression is desired.                          !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   A            (input)  : The input matrix                                 !
  !   X            (output) : The computed matrix logarithm                    !
  ! lowdin_maxiters(input)  : (Optional) impose a maximum number of iterations !
  !                           in the Lowdin orthogonalisation.
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, May 2015,                                        !
  !============================================================================!
  subroutine matrix_logarithm(A,X,lowdin_maxiters)
    use timer, only : timer_clock

    implicit none
    type(smearing_matrix), intent(inout) :: A
    type(smearing_matrix), intent(inout) :: X
    integer, optional, intent(in) :: lowdin_maxiters

    real(kind=dp), dimension(16) :: xvals
    real(kind=dp), dimension(6)  :: d
    real(kind=dp), dimension(5)  :: alpha
    integer :: nmax
    integer :: n
    integer :: s, its, k, cost, flag
    integer :: p
    real(kind=dp) :: d2, d3, d6, d7
    real(kind=dp) :: alpha2, alpha5, alpha6
    real(kind=dp) :: eta, eta6
    integer :: kp, j1, j2
    integer :: i
    integer :: mmax
    integer :: m

    type(smearing_matrix) :: scaled_A, current_A, tmp_mat
    type(smearing_matrix) :: R, Z0

    logical :: have_lowdin_maxiters

    call timer_clock('matrix_logarithm',1)

    have_lowdin_maxiters = present(lowdin_maxiters)

    ! ja531-> These values are given in the paper (Pade coeffs)
    xvals = [1.586970738772063e-005_dp, 2.313807884242979e-003_dp, 1.938179313533253e-002_dp, 6.209171588994762e-002_dp, &
         & 1.276404810806775e-001_dp, 2.060962623452836e-001_dp, 2.879093714241194e-001_dp, 3.666532675959788e-001_dp, &
         & 4.389227326152340e-001_dp, 5.034050432047666e-001_dp, 5.600071293013720e-001_dp, 6.092525642521717e-001_dp, &
         & 6.519202543720032e-001_dp, 6.888477797186464e-001_dp, 7.208340678820352e-001_dp, 7.485977242539218e-001_dp]

    mmax = 16
    n = matrix_dimension(A,1)

    s = 0
    its = 5
    k = 0
    cost = 0
    flag = 0

    call matrix_allocate(scaled_A,A)
    call matrix_allocate(current_A,A)
    call matrix_copy(A,scaled_A)
    call matrix_copy(A,current_A)
    call matrix_scale(scaled_A,1.0_dp,-1.0_dp)

    call matrix_allocate(tmp_mat,A)
    call matrix_allocate(R,A)
    call matrix_allocate(Z0,A)


    d2 = normAm(scaled_A,2)**(0.5_dp)
    d3 = normAm(scaled_A,3)**(1.0_dp/3.0_dp)

    alpha2 = max(d2,d3)
    if(alpha2 <= xvals(2)) then
       if(alpha2 <= xvals(1)) then
          m=1
       else
          m=2
       end if
       flag = 2
       !    disp('m:'), m
       !    [norm(A-eye(n),1) d2 d3]
       !    A
    end if


    kp=0
    j1=0
    j2=0
    m=0
    k=0

    do
       if(flag==2) exit

       d = 0.0_dp
       alpha = 0.0_dp
       flag = 0
       !    call matrix_copy(A,scaled_A)
       !    call matrix_scale(scaled_A,1.0_dp,-1.0_dp)
       d(3) = normAm(scaled_A,3)**(1.0_dp/3.0_dp)
       alpha(2) = huge(alpha(2))
       ploop: do p=3,5
          !       call matrix_copy(A,scaled_A)
          !       call matrix_scale(scaled_A,1.0_dp,-1.0_dp)
          d(p+1) = normAm(scaled_A,p+1)**(1.0_dp/real(p+1,dp))
          alpha(p) = max(d(p),d(p+1))
          eta = min(alpha(p-1),alpha(p))
          if(eta <= xvals(16)) then

             kp = ceiling(real(p*(p-1)-1,dp)/2.0_dp)
             do i=kp,mmax
                if(eta <= xvals(i)) then
                   j1 = i
                   exit
                end if
             end do
             do i=kp,mmax
                if(eta/2.0_dp <= xvals(i)) then
                   j2 = i
                   exit
                end if
             end do
             kp = ceiling(real(p*(p+1)-1,dp)/2.0_dp)

             if((2.0_dp*real(j1-j2,dp)/3.0_dp) < its .and. j1 <= kp) then
                m = j1
                flag = 2
                exit
             else
                if(2.0_dp*real(j1-j2,dp)/3.0_dp >= its .and. k < 2) then
                   k = k+1
                   flag = 1
                   exit ploop
                end if
             end if
          end if
       end do ploop
       d6 = d(6)
       alpha5 = alpha(5)

       if(flag==0) then
          !       call matrix_copy(A,scaled_A)
          !       call matrix_scale(scaled_A,1.0_dp,-1.0_dp)
          d7 = normAm(scaled_A,7)**(1.0_dp/7.0_dp)
          alpha6 = max(d6,d7)
          eta6 = min(alpha5,alpha6)
          if(eta6 <= xvals(mmax)) then
             if(eta6 <= xvals(15)) then
                m=15
             else
                m=mmax
             end if
          end if
       end if

       if(flag /= 2) then
          call matrix_copy(current_A,scaled_A)
          if(have_lowdin_maxiters) then
             call lowdin_transformation(scaled_A,sqrt_S_out=current_A,niter=its,maxiters=lowdin_maxiters)
          else
             call lowdin_transformation(scaled_A,sqrt_S_out=current_A,niter=its)
          end if
          call matrix_copy(current_A,scaled_A)

          call matrix_scale(scaled_A,1.0_dp,-1.0_dp)
          !        [A,~,its] = sqrtm_dbp(A,1)
          s = s+1
          cost = cost + its
          if(s == 1) then
             call matrix_copy(current_A,Z0)
             call matrix_scale(Z0,1.0_dp,-1.0_dp)
             !           Z0 = A - eye(n)
          end if
          if(s == 2) then
             call matrix_copy(current_A,R)
             call matrix_scale(R,1.0_dp,1.0_dp)
             !           R = eye(n) + A
          end if
          if(s > 2) then
             call matrix_copy(R,tmp_mat)
             call matrix_multiply(1.0_dp,tmp_mat,current_A,1.0_dp,R)
             !           R = R*(eye(n) + A)
          end if
       end if



    end do



    if(s >= 2) then
       call matrix_copy(R,tmp_mat)

       call matrix_invert(tmp_mat)

       call matrix_multiply(1.0_dp,Z0,tmp_mat,0.0_dp,R)
       !     R = Z0/R
    else
       call matrix_copy(current_A,R)
       call matrix_scale(R,1.0_dp,-1.0_dp)
       !     R = A - eye(n)
    end if
    call logm_pf(R,m,X)
    call matrix_scale(X,2.0_dp**s)
    !  X = 2^s*logm_pf(R,m)
    !  if (isreal(A)) then
    !     X = real(X)
    !  end if

    !  cost = 4*cost*n^3 + 8*m*n^3/3
    !  if(s >= 2) then
    !     cost = cost + 2*(s-2/3)*n^3
    !  end if

    call matrix_free(scaled_A)
    call matrix_free(current_A)
    call matrix_free(tmp_mat)
    call matrix_free(R)
    call matrix_free(Z0)


    call timer_clock('matrix_logarithm',2)

  contains
    subroutine logm_pf(A,m,S)
      use utils,only: utils_alloc_check, utils_dealloc_check
      !LOGM_PF   Pade approximation to matrix log by partial fraction expansion.
      !   LOGM_PF(A,m) is an [m/m] Pade approximant to LOG(EYE(SIZE(A))+A).
      implicit none
      type(smearing_matrix), intent(inout) :: A
      integer, intent(in) :: m
      type(smearing_matrix), intent(inout) :: S

      integer :: j
      real(kind=dp), dimension(:), allocatable :: nodes
      real(kind=dp), dimension(:), allocatable :: wts

      real(kind=dp) :: mineval,maxeval
      type(smearing_matrix) :: testmat,testmat2!, tmp_mat
      integer :: inverr
      integer :: ierr

      call timer_clock('matrix_logarithm_logm_pf',1)

      !call matrix_allocate(tmp_mat,A)
      !call matrix_allocate(testmat,tmp_mat)
      !call matrix_allocate(testmat2,tmp_mat)

      allocate(nodes(m),stat=ierr)
      call utils_alloc_check('logm_pf','nodes',ierr)
      allocate(wts(m),stat=ierr)
      call utils_alloc_check('logm_pf','wts',ierr)

      call gauss_legendre(m,nodes,wts)

      ! Convert from [-1,1] to [0,1].
      nodes = (nodes + 1.0_dp)/2.0_dp
      wts = wts/2.0_dp

      call matrix_scale(S,0.0_dp)

      do j=1,m
         call matrix_copy(A,tmp_mat)
         call matrix_scale(tmp_mat,nodes(j),1.0_dp)

         ! ja531 -> could we use an inverse initialisation here?

         call matrix_invert(tmp_mat)!,inverr)

         call matrix_multiply(wts(j),A,tmp_mat,1.0_dp,S)
         !       S = S + wts(j)*(A/(eye(n) + nodes(j)*A))
      end do

      deallocate(nodes,stat=ierr)
      call utils_dealloc_check('logm_pf','nodes',ierr)
      deallocate(wts,stat=ierr)
      call utils_dealloc_check('logm_pf','wts',ierr)

      !call matrix_free(testmat)
      !call matrix_free(testmat2)
      !call matrix_free(tmp_mat)
      call timer_clock('matrix_logarithm_logm_pf',2)

    end subroutine logm_pf


    function normAm(A,m) result(c)
      !NORMAM   Estimate of 1-norm of power of matrix (not as cheap as I'd like).
      implicit none
      type(smearing_matrix), intent(inout) :: A
      integer, intent(in) :: m
      real(kind=dp) :: c
      integer :: n
      type(smearing_matrix), dimension(1:2) :: X
      integer :: i

      call timer_clock('matrix_logarithm_normam',1)

      n = matrix_dimension(A,1)
      !if isequal(A,abs(A))
      !    e = ones(n,1)
      !    do j=1,m
      !        e = A'*e;
      !    end
      !    c = norm(e,inf);
      !    mv = m;
      !else
      call matrix_allocate(X(1),A,A)
      call matrix_allocate(X(2),X(1),A)
      call matrix_scale(X(1),0.0_dp)
      call matrix_scale(X(1),0.0_dp,1.0_dp)
      do i=1,m
         call matrix_multiply(1.0_dp,X(mod(i-1,2)+1),A,0.0_dp,X(mod(i,2)+1))
      end do

      c = matrix_induced_norm(X(mod(m,2)+1))

      call matrix_free(X(1))
      call matrix_free(X(2))

      call timer_clock('matrix_logarithm_normam',2)

    end function normAm

  end subroutine matrix_logarithm


    subroutine gauss_legendre(n,x,w)
      !GAUSS_LEGENDRE  Nodes and weights for Gauss-Legendre quadrature.
      !   [X,W] = GAUSS_LEGENDRE(N) computes the nodes X and weights W
      !   for N-point Gauss-Legendre quadrature.

      ! Reference:
      ! G. H. Golub and J. H. Welsch, Calculation of Gauss quadrature
      ! rules, Math. Comp., 23(106):221-230, 1969.
      use utils, only: utils_alloc_check, utils_dealloc_check, utils_safe_nint, utils_assert
      use timer, only: timer_clock
      implicit none

      ! Arguments
      integer, intent(in) :: n
      real(kind=dp), dimension(1:n), intent(out) :: x
      real(kind=dp), dimension(1:n), intent(out) :: w

      ! LAPACK subroutine
      external :: dsyev

      ! Local variables
      real(kind=dp), dimension(:,:), allocatable :: vmat
      real(kind=dp), dimension(:),   allocatable :: work
      integer :: work_size, info
      integer :: i
      integer :: ierr

      call timer_clock('matrix_logarithm_gaussleg',1)

      allocate(vmat(n,n),stat=ierr)
      call utils_alloc_check('gauss_legendre','vmat',ierr)

      vmat=0.0_dp
      do i = 1,n-1
         vmat(i,i+1) = real(i,dp)/sqrt(real(((2*i)**2)-1,dp))
         vmat(i+1,i) = vmat(i,i+1)
      end do

      work_size=-1
      allocate(work(1),stat=ierr)
      call utils_alloc_check('gauss_legendre','work',ierr)
      if(debug_info_toggle) write(stdout,*) "smearing-> dsyev : ", n, size(vmat,1), size(x), size(work), work_size, "<== SMEAR"
      call dsyev('V','U',n,vmat,size(vmat,1),x,work,work_size,info)
      call utils_assert(info==0, "Error in gauss legendre : dsyev workspace query failed")
      work_size = utils_safe_nint(work(1))
      deallocate(work,stat=ierr)
      call utils_dealloc_check('gauss_legendre','work',ierr)
      allocate(work(work_size),stat=ierr)
      call utils_alloc_check('gauss_legendre','work',ierr)
      call dsyev('V','U',n,vmat,size(vmat,1),x,work,work_size,info)
      call utils_assert(info==0, "Error in gauss legendre : dsyev failed")
      w = 2.0_dp*(vmat(1,:)**2)

      deallocate(work,stat=ierr)
      call utils_dealloc_check('gauss_legendre','work',ierr)
      deallocate(vmat,stat=ierr)
      call utils_dealloc_check('gauss_legendre','vmat',ierr)

      call timer_clock('matrix_logarithm_gaussleg',2)
    end subroutine gauss_legendre


  !============================================================================!
  ! The problem with the Fermi-Dirac entropy function is that we want to       !
  ! evaluate arguments close to 0 and 1. These end up being 0 arguments to the !
  ! ln() function. To avoid this, we can impose a cutoff and set the results   !
  ! of those arguments to zero as this algorithm does. This is achieved for    !
  ! matrices by setting up a projector which projects these values to 1.       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   A            (input)  : The input matrix                                 !
  !   X            (output) : The computed matrix logarithm                    !
  !   mintol       (input)  : The cutoff (1e-11 if not set) (Optional)
  ! lowdin_maxiters(input)  : (Optional) impose a maximum number of iterations !
  !                           in the Lowdin orthogonalisation.
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, May 2015,                                        !
  !============================================================================!
  subroutine matrix_safe_logarithm(A,X,mintol,lowdin_maxiters)
    use comms, only: comms_reduce, pub_on_root
    use timer, only: timer_clock
    implicit none
    type(smearing_matrix), intent(inout) :: A
    type(smearing_matrix), intent(inout) :: X
    real(kind=dp), optional, intent(in) :: mintol
    integer, optional, intent(in) :: lowdin_maxiters

    type(smearing_matrix) :: tmp_mat, tmp_mat2
    real(kind=dp) :: mintol_loc
    real(kind=dp) :: sign_tol
    integer       :: matmuls_loc
    integer :: status
    logical :: have_lowdin_maxiters
    real(kind=dp) :: shift
    integer :: i
    logical :: logfail=.false.
    real(kind=dp) :: tmpvec(4)

    real(kind=dp) :: abs_emax

    call timer_clock('matrix_safe_logarithm',1)

    have_lowdin_maxiters=present(lowdin_maxiters)

    if(present(mintol)) then
       mintol_loc=abs(mintol)
    else
       mintol_loc=1e-11_dp
    end if

    sign_tol=1e-13_dp

    call matrix_allocate(tmp_mat,A)

    call matrix_copy(A,tmp_mat)

    call matrix_allocate(tmp_mat2,tmp_mat)
    shift=0.0_dp
    i=0
    do
       i=i+1
       call matrix_scale(tmp_mat,1.0_dp,-(mintol_loc+shift))
       matmuls_loc=0
       abs_emax=-1.0_dp
       call sign_NSH(tmp_mat,sign_tol,tmp_mat2,abs_emax,matmuls_loc,status)
       tmpvec(1)=matrix_norm(tmp_mat,2)
       tmpvec(2)=matrix_trace(tmp_mat)
       tmpvec(3)=matrix_norm(tmp_mat2,2)
       tmpvec(4)=matrix_trace(tmp_mat2)
       if(pub_on_root) then
          if(debug_info_toggle) write(stdout,*) "smearing-> sign nsh in safe log", tmpvec, mintol_loc, shift, "<== SMEAR"
       end if
       if(status/=0) then
          call random_number(shift)
          shift=shift*10.0_dp**(log10(mintol_loc)-2+i)
          cycle
       end if
       call matrix_scale(tmp_mat2,-0.5_dp,0.5_dp)

       call matrix_axpy(tmp_mat2,A,1.0_dp)

       if(have_lowdin_maxiters) then
          call matrix_logarithm(tmp_mat2,X,lowdin_maxiters=lowdin_maxiters)
       else
          call matrix_logarithm(tmp_mat2,X)
       end if

       logfail=.not.matrix_any_isnan(X)
       if(A%matrix_type/=matrix_type_standard) then
          call comms_reduce("OR",logfail)
       end if
       if(logfail) then
          exit
       else
          call matrix_scale(tmp_mat2,0.0_dp)
          call matrix_scale(X,0.0_dp)
          call random_number(shift)
          shift=shift*10.0_dp**(log10(mintol_loc)-2+i)
       end if

    end do

    call matrix_free(tmp_mat)
    call matrix_free(tmp_mat2)

    call timer_clock('matrix_safe_logarithm',2)

  end subroutine matrix_safe_logarithm


  !============================================================================!
  ! This routine uses a locally preconditioned conjugate gradients algorithm   !
  ! to find extremal eigenvalues of an (possibly generalized) eigenproblem     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   A            (input)  : The input matrix                                 !
  !   M            (input)  : The input metric (I if standard e-problem)       !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, May 2015,                                        !
  !============================================================================!
  subroutine local_precond_cg(A,M,uplo,precond,tol,eval,const,normal)!,evec)
    use comms, only: comms_bcast, pub_root_proc_id, pub_on_root

    use rundat, only: pub_debug_on_root
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, &
         utils_assert, utils_safe_nint
    implicit none
    type(smearing_matrix),                             intent(inout)  :: A
    type(smearing_matrix),                             intent(inout)  :: M
    character(len=*),                         intent(in)  :: uplo
    character(len=*),                         intent(in)  :: precond
    real(kind=dp),                            intent(in)  :: tol
    real(kind=dp),                            intent(out) :: eval
    real(kind=dp),                  optional, intent(in)  :: const
    logical,                        optional, intent(in)  :: normal ! Normal eigenproblem... i.e. ignore M
    !    real(kind=dp), dimension(:), optional, intent(out) :: evec

    ! LAPACK subroutine
    external :: dsygv

    type(smearing_matrix) :: B

    real(kind=dp) :: rnums(2)
    real(kind=dp) :: vnorm,unorm
    type(vector) :: x
    type(vector) :: u
    type(vector) :: v
    type(vector) :: g
    type(vector) :: p
    type(vector) :: tmpvec, tmpvec2
    real(kind=dp), dimension(:),   allocatable :: locx, locu
    real(kind=dp), dimension(:,:), allocatable :: tmpmat1
    real(kind=dp), dimension(:,:), allocatable :: tmpmat2
    real(kind=dp), dimension(:,:), allocatable :: mm
    real(kind=dp), dimension(:,:), allocatable :: aa
    real(kind=dp) :: q
    real(kind=dp) :: gnorm
    real(kind=dp) :: eig
    real(kind=dp) :: loctol, outvec(6)
    integer :: i, k
    integer :: N
    logical :: lower, upper

    integer :: info
    integer :: work_size
    real(kind=dp), dimension(:), allocatable :: work
    real(kind=dp), dimension(:), allocatable :: rwork
    real(kind=dp), dimension(:), allocatable :: delta
    real(kind=dp), dimension(:), allocatable :: eigs

    logical :: converged

    real(kind=dp) :: alloc_time,free_time,matvec_time,scale_time,bcast_time,dot_time,axpy_time,eig_time,copy_time, &
         & total_time,norm_time,bef_time
    real(kind=dp),save :: free_time_tot=0.0_dp, matvec_time_tot=0.0_dp, scale_time_tot=0.0_dp, bcast_time_tot=0.0_dp, &
         & dot_time_tot=0.0_dp, axpy_time_tot=0.0_dp, eig_time_tot=0.0_dp, copy_time_tot=0.0_dp, alloc_time_tot=0.0_dp, &
         & total_time_tot=0.0_dp, norm_time_tot=0.0_dp
    integer,save :: entries=0

    logical :: loc_normal
    integer :: ierr
    ! agrecocmplx
    logical :: loc_cmplx

    call timer_clock('local_precond_cg',1)

    ! jmecmplx
    loc_cmplx = ( A%standard_is_cmplx .or. M%standard_is_cmplx )

    ! agrecocmplx
    call utils_assert(.not. loc_cmplx, &
         'Subroutine local_precond_cg not ready yet &
         & for complex NGWFs.')

    loc_normal=.true.
    if(present(normal)) loc_normal=normal

    entries=entries+1

    if(uplo=="upper") then
       lower=.false.
       upper=.true.
    else if(uplo=="lower") then
       lower=.true.
       upper=.false.
    else
       call utils_abort("Error in local_precond_cg : unrecognised uplo.")
    end if
    N=matrix_dimension(A)
    !    allocate(x(N))
    !    allocate(u(N))
    !    allocate(v(N))
    !    allocate(g(N))
    !    allocate(p(N))
    !    allocate(tmpvec(N))


    free_time   = 0.0_dp
    matvec_time = 0.0_dp
    scale_time  = 0.0_dp
    bcast_time  = 0.0_dp
    dot_time    = 0.0_dp
    norm_time   = 0.0_dp
    axpy_time   = 0.0_dp
    eig_time    = 0.0_dp
    copy_time   = 0.0_dp
    alloc_time  = 0.0_dp
    total_time  = 0.0_dp

    x%matrix_type=A%matrix_type
    call vector_allocate(x,N,rand=.true.)
    call vector_allocate(u,x)
    call vector_allocate(v,x)
    call vector_allocate(g,x)
    call vector_allocate(p,x)
    call vector_allocate(tmpvec,x)
    call vector_allocate(tmpvec2,x)
    !    allocate(tmpmat1(3,N))
    !    allocate(tmpmat2(N,3))

    allocate(aa(3,3),stat=ierr)
    call utils_alloc_check('local_precond_cg','aa',ierr)
    allocate(mm(3,3),stat=ierr)
    call utils_alloc_check('local_precond_cg','mm',ierr)
    allocate(delta(3),stat=ierr)
    call utils_alloc_check('local_precond_cg','delta',ierr)
    allocate(eigs(3),stat=ierr)
    call utils_alloc_check('local_precond_cg','eigs',ierr)


    vnorm=1.0_dp/vector_norm(x)
    if(A%matrix_type/=matrix_type_standard) then
       call comms_bcast(pub_root_proc_id,vnorm)
    end if

    call vector_scale(x,vnorm)

    if(.not.loc_normal) then
       call matrix_multiply_vec(M,x,u)

       q=sqrt(vector_dot(x,u))

       if(A%matrix_type/=matrix_type_standard) then
          call comms_bcast(pub_root_proc_id,q)
       end if

       call vector_scale(x,1.0_dp/q)
       call vector_scale(u,1.0_dp/q)

    else
       q=vector_norm(x)
       if(A%matrix_type/=matrix_type_standard) then
          call comms_bcast(pub_root_proc_id,q)
       end if
       call vector_scale(x,1.0_dp/q)
    end if
    !    x=x/q
    !    u=u/q

    call matrix_multiply_vec(A,x,v)
    !    if(lower) then
    !       v=-v
    !    end if
    loctol=tol
    !    eig=dot_product(conjg(x),v)
    eig=vector_dot(x,v)

    k=0
    gnorm=1.0_dp
    eigsearch: do
       if(A%matrix_type/=matrix_type_standard) then
          call comms_bcast(pub_root_proc_id,eig)
       end if

       if(present(const)) then
          loctol=10.0_dp**(log10(tol)-abs(nint(log10(abs(eig))-log10(abs(eig+const)))))

       end if
       converged=gnorm <= loctol .or. k>100

       if(A%matrix_type/=matrix_type_standard) then
          call comms_bcast(pub_root_proc_id,converged)
       end if

       if(converged) exit eigsearch
       k=k+1

       call vector_copy(v,g)
       if(.not.loc_normal) then
          call vector_axpy(g,-eig,u)
       else
          call vector_axpy(g,-eig,x)
       end if

       !       g=v-eig*u

       gnorm=vector_norm(g)
       if(A%matrix_type/=matrix_type_standard) then
          call comms_bcast(pub_root_proc_id,gnorm)
       end if

       select case(precond)
       case("inverse")
          if(k==1) then
             call matrix_allocate(B,A)
             B=A

             call matrix_invert(B)
          end if
          call matrix_multiply_vec(B,g,tmpvec)
       case("gaussseidel")
          if(pub_on_root) then
             write(stdout,*) "Stub: Gauss-Seidel preconditioner not implemented yet. Doing nothing!"
          end if
          ! Should be :
          !           : tmpvec=matmul(inverse(upper_tri(A)),matmul(diag(A),inverse(lower_tri(A)*g)))
          tmpvec=g
       case("none")
          !          tmpvec=g
          !bef_time=mpi_wtime()
          !          call vector_copy(g,tmpvec)
          !copy_time=copy_time+(mpi_wtime()-bef_time)
       case default
          !          tmpvec=g
          !bef_time=mpi_wtime()
          !          call vector_copy(g,tmpvec)
          !copy_time=copy_time+(mpi_wtime()-bef_time)
       end select

       if(upper) then
          if(precond=="inverse".or.precond=="gaussseidel") then
             call vector_copy(tmpvec,g)
          end if
          !          g=-tmpvec
       else if(lower) then
          if(precond=="inverse".or.precond=="gaussseidel") then
             call vector_copy(tmpvec,g) ! g=tmpvec
          end if
          call vector_scale(g,-1.0_dp)
       end if

       if(k==1) then
          call vector_scale(p,0.0_dp) !p=0.0_dp
       end if

       call matrix_multiply_vec(A,g,tmpvec)
       call matrix_multiply_vec(A,p,tmpvec2)

       aa(1,1)=vector_dot(x,v)
       aa(2,1)=vector_dot(g,v)
       aa(3,1)=vector_dot(p,v)
       aa(1,2)=vector_dot(x,tmpvec)
       aa(2,2)=vector_dot(g,tmpvec)
       aa(3,2)=vector_dot(p,tmpvec)
       aa(1,3)=vector_dot(x,tmpvec2)
       aa(2,3)=vector_dot(g,tmpvec2)
       aa(3,3)=vector_dot(p,tmpvec2)

       !       aa=(aa+conjg(transpose(aa)))/2.0_dp
       aa=(aa+(transpose(aa)))/2.0_dp

       if(.not.loc_normal) then
          call matrix_multiply_vec(M,g,tmpvec)
          call matrix_multiply_vec(M,p,tmpvec2)

          mm(1,1)=vector_dot(x,u)
          mm(2,1)=vector_dot(g,u)
          mm(3,1)=vector_dot(p,u)
          mm(1,2)=vector_dot(x,tmpvec)
          mm(2,2)=vector_dot(g,tmpvec)
          mm(3,2)=vector_dot(p,tmpvec)
          mm(1,3)=vector_dot(x,tmpvec2)
          mm(2,3)=vector_dot(g,tmpvec2)
          mm(3,3)=vector_dot(p,tmpvec2)


          !       mm=(mm+conjg(transpose(mm)))/2.0_dp
          mm=(mm+(transpose(mm)))/2.0_dp

       else
          mm(1,1)=1.0_dp            !vector_dot(x,x)
          mm(3,2)=vector_dot(p,g)
          mm(2,1)=vector_dot(g,x)
          mm(3,1)=vector_dot(p,x)
          mm(2,2)=vector_norm(g)**2 !vector_dot(g,g)
          mm(3,3)=vector_norm(p)**2 !vector_dot(p,p)
          mm(1,2)=mm(2,1)           !vector_dot(x,g)
          mm(1,3)=mm(3,1)           !vector_dot(x,p)
          mm(2,3)=mm(3,2)           !vector_dot(g,p)
       end if

       work_size=-1
       eigs=0.0_dp
       if(k==1) then
          i=2
       else
          i=3
       end if
       allocate(work(1),stat=ierr)
       call utils_alloc_check('local_precond_cg','work',ierr)
       call dsygv(1,'V','U',i,aa,size(aa,1),mm,size(mm,1),eigs,work,work_size,info)
       work_size = utils_safe_nint(work(1))
       deallocate(work,stat=ierr)
       call utils_dealloc_check('local_precond_cg','work',ierr)
       allocate(work(work_size),stat=ierr)
       call utils_alloc_check('local_precond_cg','work',ierr)
       call dsygv(1,'V','U',i,aa,size(aa,1),mm,size(mm,1),eigs,work,work_size,info)
       deallocate(work,stat=ierr)
       call utils_dealloc_check('local_precond_cg','work',ierr)

       ! If info > 0 and we are not on the first iteration, then tolerance was probably too optimistic
       ! Precision fell off a cliff and made our subspace matrix non SPD. This could do with better
       ! checks, but for now, just returning the eigenvalue from the previous iteration.

       if(A%matrix_type/=matrix_type_standard) then
          call comms_bcast(pub_root_proc_id,info)
       end if

       if(info > 0) then
          exit eigsearch
       end if

       if(upper) then
          i=maxloc(eigs(1:i),1)
          eig=maxval(eigs(1:i))
       else if(lower) then
          i=minloc(eigs(1:i),1)
          eig=minval(eigs(1:i))
       end if
       delta=aa(:,i)

       if(A%matrix_type/=matrix_type_standard) then
          call comms_bcast(pub_root_proc_id,delta)
       end if

       call vector_scale(p,delta(3))
       call vector_axpy(p,delta(2),g)

       call vector_scale(x,delta(1))
       call vector_axpy(x,1.0_dp,p)


       if(.not.loc_normal) then
          call matrix_multiply_vec(M,x,u)
          q=1.0_dp/sqrt(vector_dot(x,u))

          if(A%matrix_type/=matrix_type_standard) then
             call comms_bcast(pub_root_proc_id,q)
          end if

          !       x=x/q
          !       u=u/q
          call vector_scale(x,q)
          call vector_scale(u,q)
       else
          q=vector_norm(x)
          if(A%matrix_type/=matrix_type_standard) then
             call comms_bcast(pub_root_proc_id,q)
          end if
          call vector_scale(x,1.0_dp/q)
       end if

       call matrix_multiply_vec(A,x,v)
       !       if(lower) then
       !         v=-v
       !       end if


    end do eigsearch

    eval=eig
    call matrix_free(B)
    call vector_free(x)
    call vector_free(u)
    call vector_free(v)
    call vector_free(g)
    call vector_free(tmpvec)
    call vector_free(tmpvec2)

    deallocate(aa,stat=ierr)
    call utils_dealloc_check('local_precond_cg','aa',ierr)
    deallocate(mm,stat=ierr)
    call utils_dealloc_check('local_precond_cg','mm',ierr)
    deallocate(delta,stat=ierr)
    call utils_dealloc_check('local_precond_cg','delta',ierr)
    deallocate(eigs,stat=ierr)
    call utils_dealloc_check('local_precond_cg','eigs',ierr)

    total_time=alloc_time+free_time+matvec_time+scale_time+bcast_time+dot_time+axpy_time+eig_time+copy_time+norm_time

    alloc_time_tot=alloc_time_tot+alloc_time
    matvec_time_tot=matvec_time_tot+matvec_time
    scale_time_tot=scale_time_tot+scale_time
    bcast_time_tot=bcast_time_tot+bcast_time
    dot_time_tot=dot_time_tot+dot_time
    norm_time_tot=norm_time_tot+norm_time
    axpy_time_tot=axpy_time_tot+axpy_time
    eig_time_tot=eig_time_tot+eig_time
    copy_time_tot=copy_time_tot+copy_time
    free_time_tot=free_time_tot+free_time
    total_time_tot=total_time_tot+total_time

    if(pub_debug_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> alloc time  = ", alloc_time  , alloc_time_tot , "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> matvec time = ", matvec_time , matvec_time_tot, "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> scale time  = ", scale_time  , scale_time_tot , "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> bcast time  = ", bcast_time  , bcast_time_tot , "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> dot time    = ", dot_time    , dot_time_tot, "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> norm time   = ", norm_time   , norm_time_tot, "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> axpy time   = ", axpy_time   , axpy_time_tot  , "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> eig time    = ", eig_time    , eig_time_tot   , "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> copy time   = ", copy_time   , copy_time_tot  , "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> free time   = ", free_time   , free_time_tot  , "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> total time  = ", total_time  , total_time_tot , "<== TIME"
    end if


    call timer_clock('local_precond_cg',2)
  contains

    ! This routine computes the nth order norm of
    ! a vector.
    function internal_vecnorm(x,ord) result(norm)
      implicit none
      real(kind=dp) :: norm
      real(kind=dp), dimension(:), intent(in) :: x
      integer, intent(in) :: ord
      integer :: i
      norm=0.0_dp
      do i=1,size(x)
         norm=norm+(abs(x(i))**ord)
      end do
      norm=norm**(1.0_dp/real(ord,dp))
    end function internal_vecnorm

  end subroutine local_precond_cg


  !============================================================================!
  ! This routine computes the highest and lowest eigenvalues of a matrix,      !
  ! using the locally preconditioned cg algorithm (see above). Unlike the uplo !
  ! parameter in that routine, this routine computes the true lowest eigenval  !
  ! by using the Gerschgorin circle theorem to put an upper bound on the range !
  ! of eigenvalues, and then scale and shift, to calculate the lower eig of    !
  ! the modified system and back-transform.                                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   A            (input)  : The input matrix                                 !
  !   M            (input)  : The rotation matrix in the gen-eig problem       !
  !   tol          (input)  : Desired accuracy of eigenvalue                   !
  !   mineval      (output) : The minimum eigenvalue                           !
  !   maxeval      (output) : The maximum eigenvalue                           !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Jan 2014,                                        !
  ! Modified to handle SPAM3_EMBED by Robert Charlton, Aug/Sept 2018.          !
  !----------------------------------------------------------------------------!
  ! Really ought to sort out better conditioning.                              !
  !============================================================================!
  subroutine extremum_eigs(A,M,tol,mineval,maxeval)
    use comms, only : pub_root_proc_id, comms_bcast
    use sparse,only : sparse_extremal_eigenvalue
    use sparse_embed, only : sparse_embed_extremal_eigenvalue
    use timer, only : timer_clock
    implicit none
    type(smearing_matrix),  intent(inout)  :: A ! Inout for dense_create_copy in matrix_allocate (in)
    type(smearing_matrix),  intent(inout)  :: M ! Inout for dense_create_copy in local_precond_cg (in)
    real(kind=dp), intent(in)  :: tol
    real(kind=dp), optional, intent(out) :: mineval
    real(kind=dp), optional, intent(out) :: maxeval

    type(smearing_matrix) :: B
    real(kind=dp) :: twonorm,fac
    real(kind=dp), dimension(:), allocatable :: geigs
    real(kind=dp),    dimension(:), allocatable :: grads
    integer :: i, N

    call timer_clock('smearing_extremum_eigs',1)

    if(.not.present(mineval).and..not.present(maxeval)) then
       call timer_clock('smearing_extremum_eigs',2)
       return
    end if

    N=matrix_dimension(A)

    call matrix_allocate(B,A)
    call matrix_copy(A,B)

    twonorm=matrix_norm(B,2)
    if(A%matrix_type/=matrix_type_standard) then
       call comms_bcast(pub_root_proc_id,twonorm)
    end if

    !    call gerschgorin(B,geigs,grads)
    !    fac=abs(minval(real(geigs,dp)-real(grads,dp)))

    fac=matrix_induced_norm(B)

    if(A%matrix_type/=matrix_type_standard) then
       call comms_bcast(pub_root_proc_id,fac)
    end if

    if(present(maxeval)) then
       ! rc2013: allow SPAM3_EMBED option too
       select case(A%matrix_type)
       case(matrix_type_SPAM3)
          call sparse_extremal_eigenvalue(B%dataSPAM3,M%dataSPAM3,maxeval,tol)
       case(matrix_type_SPAM3_EMBED)
          call sparse_embed_extremal_eigenvalue(B%dataSPAM3_EMBED,&
               M%dataSPAM3_EMBED,maxeval,tol)
       case default
          call local_precond_cg(B,M,"upper","none",tol,maxeval)
       end select
    end if

    if(present(mineval)) then
       call matrix_scale(B,-1.0_dp,fac)
       call matrix_scale(B,1.0_dp/twonorm)
       ! rc2013: allow SPAM3_EMBED option too
       select case(A%matrix_type)
       case(matrix_type_SPAM3)
          call sparse_extremal_eigenvalue(B%dataSPAM3,M%dataSPAM3,maxeval,tol)
       case(matrix_type_SPAM3_EMBED)
          call sparse_embed_extremal_eigenvalue(B%dataSPAM3_EMBED,&
               M%dataSPAM3_EMBED,maxeval,tol)
       case default
          call local_precond_cg(B,M,"upper","none",tol,mineval,const=-fac)
       end select

       mineval=(fac-(twonorm*mineval))


    end if
    call matrix_free(B)


    call timer_clock('smearing_extremum_eigs',2)

  end subroutine extremum_eigs




  !============================================================================!
  ! Since a Cholesky decomposition is currently unavailable in the Sparse mod, !
  ! instead a Lowdin transformation of the overlap matrix is used to           !
  ! orthogonalise the Hamiltonian. This routine calculates this with a         !
  ! modified Newton-Schultz-Hotelling algorithm, of the first order.           !
  ! By using this approach, the L^(-1/2) matrix can be obtained in linear      !
  ! scaling time (for sparse matrices). Where S = L^(-1/2)*L^(-1/2)'           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   A            (input)  : The input matrix                                 !
  !   M            (input)  : The rotation matrix in the gen-eig problem       !
  !   uplo         (input)  : Whether we want highest or lowest (trunc) eigval !
  !   precond      (input)  : Which preconditioner to use i.e. none/inverse... !
  !   tol          (input)  : Desired accuracy of eigenvalue                   !
  !   eval         (output) : The eigenvalue.                                  !
  !   evec (opt)   (output) : The eigenvector.                                 !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Jan 2014,                                        !
  !----------------------------------------------------------------------------!
  ! Lanczos?                                                                   !
  ! Sparsity patterns of intermediate matrices may need more thought? (NB1)    !
  !============================================================================!
  subroutine lowdin_transformation(overlap,inv_sqrt_S_out,tol,sqrt_S_out,niter,maxiters,status)
    use comms, only : pub_on_root
    use timer, only : timer_clock
    use utils, only : utils_alloc_check, utils_dealloc_check

    implicit none
    type(smearing_matrix),           intent(inout)  :: overlap ! Inout for dense_create_copy in matrix_allocate (in)
    type(smearing_matrix), optional, intent(inout) :: inv_sqrt_S_out
    real(kind=dp),optional, intent(in)  :: tol
    type(smearing_matrix), optional, intent(inout) :: sqrt_S_out
    integer,      optional, intent(out)   :: niter
    integer,      optional, intent(in)    :: maxiters
    integer,      optional, intent(out)   :: status

    type(smearing_matrix)  :: sqrt_S, inv_sqrt_S, Xk, tmp_mat, tmp_mat2
    integer       :: i, j
    real(kind=dp) :: epsmax, epsmin, lambda
    logical       :: finaliter
    real(kind=dp) :: loctol
    real(kind=dp) :: trace_powers(0:4)
    real(kind=dp) :: fake_vec(1:5)
    real(kind=dp) :: opt_tol
    integer       :: loc_stat
    integer       :: method

    integer, parameter :: powermethod=0, frobfunc=1, gersfunc=2
    logical :: have_maxiters
    real(kind=dp) :: old_norm, new_norm

    integer :: recorded_hist_points
    integer, parameter :: max_hist_points=3
    integer, parameter :: num_positive_to_stop=5
    integer, parameter :: num_positive_to_reset=3
    integer, parameter :: miniters=3
    integer :: num_positive_grads
    real(kind=dp), allocatable :: conv_history(:), coeffs(:)

    real(kind=dp) :: total_time, alloc_time, prod_time, copy_time, set_time, scale_time, free_time, bef_time, trace_time

    integer :: ierr

    call timer_clock('lowdin_transformation',1)

    total_time=0.0_dp
    alloc_time=0.0_dp
    prod_time=0.0_dp
    copy_time=0.0_dp
    set_time=0.0_dp
    scale_time=0.0_dp
    free_time=0.0_dp
    bef_time=0.0_dp
    trace_time=0.0_dp

    loc_stat=0

    if(present(tol)) then
       loctol=tol
    else
       loctol=10e-9_dp
    end if
    if(.not.present(inv_sqrt_S_out).and..not.present(sqrt_S_out)) then
       loc_stat=-1
       if(present(status)) then
          status=loc_stat
       end if
       call timer_clock('lowdin_transformation',2)
       return
    end if

    !    bef_time=mpi_wtime()
    ! ja531 -> NB1
    if(present(sqrt_S_out)) then
       call matrix_allocate(sqrt_S,sqrt_S_out)
    else
       if(present(inv_sqrt_S_out)) then
          call matrix_allocate(sqrt_S,inv_sqrt_S_out)
       else
          call matrix_allocate(sqrt_S,overlap,overlap)
       end if
    end if
    if(present(inv_sqrt_S_out)) then
       call matrix_allocate(inv_sqrt_S,inv_sqrt_S_out)
    else
       if(present(sqrt_S_out)) then
          call matrix_allocate(inv_sqrt_S,sqrt_S_out,sqrt_S_out)
       else
          call matrix_allocate(inv_sqrt_S,overlap,overlap)
       end if
    end if

    call matrix_copy(overlap,sqrt_S)

    call matrix_scale(inv_sqrt_S,0.0_dp)

    call matrix_scale(inv_sqrt_S,0.0_dp,1.0_dp)

    ! ja531 -> NB1
    if(present(inv_sqrt_S_out)) then
       call matrix_allocate(tmp_mat,inv_sqrt_S_out)
    elseif(present(sqrt_S_out)) then
       call matrix_allocate(tmp_mat,sqrt_S_out,sqrt_S_out)
    end if

    method=frobfunc
    if(method==powermethod) then ! Iterative method
       call extremum_eigs(overlap,inv_sqrt_S,loctol,epsmin,epsmax)
       lambda=2.0_dp/(epsmax+epsmin)
    else if(method==frobfunc) then ! Frobenius approximation
       opt_tol=1e-4_dp
       call matrix_allocate(tmp_mat2,overlap)
       call matrix_scale(tmp_mat,0.0_dp)
       call matrix_scale(tmp_mat2,0.0_dp)
       call matrix_copy(overlap,tmp_mat2)
       call matrix_multiply(1.0_dp,overlap,tmp_mat2,0.0_dp,tmp_mat)

       trace_powers(0)=matrix_trace(inv_sqrt_S)
       trace_powers(1)=matrix_trace(overlap)
       trace_powers(2)=matrix_trace(tmp_mat)
       trace_powers(3)=matrix_trace(tmp_mat,overlap)
       trace_powers(4)=matrix_trace(tmp_mat,tmp_mat)

       fake_vec(1)=0.0_dp ! Start
       fake_vec(2)=1.0_dp ! Middle
       fake_vec(3)=1.0_dp ! End

       call internal_bracket_minimum(frobfunc, fake_vec(1:1), fake_vec(3:3), loc_stat, middle=fake_vec(2:2))
       call internal_brents_method(frobfunc, fake_vec(1:1), fake_vec(2:2), fake_vec(3:3), &
            & opt_tol, fake_vec(4:4), fake_vec(5), loc_stat)

       if(loc_stat==-2) then
          call extremum_eigs(overlap,inv_sqrt_S,loctol,epsmin,epsmax)
          lambda=2.0_dp/(epsmax+epsmin)
          if(pub_on_root.and.debug_info_toggle) then
             write(stdout,*) "smearing-> eig lambda = ",lambda,epsmax,epsmin
          end if

       else

          lambda=fake_vec(4)
          if(lambda<0.0_dp) then
             loc_stat=-2
             if(present(status)) then
                status=loc_stat
             end if
             call matrix_free(sqrt_S)
             call matrix_free(inv_sqrt_S)
             call matrix_free(tmp_mat)
             call matrix_free(tmp_mat2)
             call timer_clock('lowdin_transformation',2)
             return
          end if
       end if

       if(pub_on_root.and.debug_info_toggle) then
          write(stdout,*) "smearing-> est lambda = ",lambda, internal_frobenius_est([lambda]), "<== SMEAR"
       end if

       !lambda=10.0_dp
       j=0
       !       do
       !          j=j+1
       !       bef_time=mpi_wtime()
       call matrix_copy(overlap,tmp_mat2) ! Is this needed? Scale should take care of it.
       !       copy_time=copy_time+(mpi_wtime()-bef_time)

       !       bef_time=mpi_wtime()
       old_norm=sqrt(matrix_trace(tmp_mat2,tmp_mat2))
       !       trace_time=trace_time+(mpi_wtime()-bef_time)

       !       bef_time=mpi_wtime()
       call matrix_scale(tmp_mat2,lambda,-1.0_dp)
       !       scale_time=scale_time+(mpi_wtime()-bef_time)

       !       bef_time=mpi_wtime()
       new_norm=sqrt(matrix_trace(tmp_mat2,tmp_mat2))
       !       trace_time=trace_time+(mpi_wtime()-bef_time)
       if(pub_on_root.and.debug_info_toggle) then
          write(stdout,*) "smearing-> lambda = ",lambda, old_norm,new_norm, "<== SMEAR"
       end if

       if(new_norm>old_norm) then
          !          if(sqrt(matrix_trace(tmp_mat2,tmp_mat2))>1.0_dp) then
          lambda=0.9_dp*lambda
          !          else
          !             exit
       end if
       !          if(j>3) then
       !             loc_stat=-3
       !             if(present(status)) then
       !                status=loc_stat
       !             end if
       !             call matrix_free(sqrt_S)
       !             call matrix_free(inv_sqrt_S)
       !             call matrix_free(tmp_mat)
       !             call matrix_free(tmp_mat2)
       !             call timer_clock('lowdin_transformation',2)
       !             return
       !          end if
       !       end do
       !          bef_time=mpi_wtime()
       call matrix_free(tmp_mat2)
       !          free_time=free_time+(mpi_wtime()-bef_time)
    else if(method==gersfunc) then ! Gerschgorin approximation

    end if

    recorded_hist_points=0
    allocate(conv_history(max_hist_points),stat=ierr)
    call utils_alloc_check('lowdin_transformation','conv_history',ierr)
    allocate(coeffs(2),stat=ierr)
    call utils_alloc_check('lowdin_transformation','coeffs',ierr)
    num_positive_grads=0


    have_maxiters=present(maxiters)

    finaliter=.false.
    i=0
    !    bef_time=mpi_wtime()
    ! ja531 -> NB1
    call matrix_allocate(Xk,tmp_mat)


    !    alloc_time=alloc_time+(mpi_wtime()-bef_time)

    do
       i=i+1

       old_norm=new_norm



       if(num_positive_grads>num_positive_to_reset) then
          if(pub_on_root.and.debug_info_toggle) then
             write(stdout,*) "smearing-> low reset = ",i,lambda, "<== SMEAR"
          end if
          !          bef_time=mpi_wtime()
          call matrix_copy(overlap,sqrt_S)
          !          copy_time=copy_time+(mpi_wtime()-bef_time)
          !          bef_time=mpi_wtime()
          call matrix_scale(inv_sqrt_S,0.0_dp)
          !          set_time=set_time+(mpi_wtime()-bef_time)
          !          bef_time=mpi_wtime()
          call matrix_scale(inv_sqrt_S,0.0_dp,1.0_dp)
          !          scale_time=scale_time+(mpi_wtime()-bef_time)
          lambda=0.9_dp*lambda
       end if

       !       bef_time=mpi_wtime()
       call matrix_multiply(lambda,sqrt_S,inv_sqrt_S,0.0_dp,Xk)
       !       prod_time=prod_time+(mpi_wtime()-bef_time)
       !       bef_time=mpi_wtime()
       call matrix_scale(Xk,1.0_dp,-1.0_dp)
       !       scale_time=scale_time+(mpi_wtime()-bef_time)


       !       bef_time=mpi_wtime()
       new_norm=sqrt(matrix_trace(Xk,Xk))
       !       trace_time=trace_time+(mpi_wtime()-bef_time)

       recorded_hist_points=recorded_hist_points+1
       if(recorded_hist_points<=max_hist_points) then
          conv_history(recorded_hist_points)=new_norm
          call poly_fit(1, conv_history, coeffs, recorded_hist_points)
       else
          conv_history(1:max_hist_points-1)=conv_history(2:max_hist_points)
          conv_history(max_hist_points)=new_norm
          call poly_fit(1, conv_history, coeffs)
       end if

       if(pub_on_root.and.debug_info_toggle) then
          write(stdout,*) "smearing-> low conv = ",i,lambda, new_norm, coeffs, "<== SMEAR"
       end if

       if(i>=miniters) then

          if(coeffs(2)>0.0_dp) then
             num_positive_grads=num_positive_grads+1
          else
             num_positive_grads=0
          end if

          if(num_positive_grads>num_positive_to_stop) then
             if(pub_on_root) then
                if(debug_info_toggle) then
                   write(stdout,*) "smearing-> stopping lowdin ==> not converging any more : ", i, " : ", coeffs,"<== SMEAR"
                end if
             end if
             loc_stat=-4
             if(present(status)) then
                status=loc_stat
             end if
             call matrix_free(sqrt_S)
             call matrix_free(inv_sqrt_S)
             call matrix_free(tmp_mat)
             call matrix_free(Xk)
             deallocate(conv_history,stat=ierr)
             call utils_dealloc_check('lowdin_transformation','conv_history',ierr)
             deallocate(coeffs,stat=ierr)
             call utils_dealloc_check('lowdin_transformation','coeffs',ierr)
             call timer_clock('lowdin_transformation',2)
             return
          end if

          if(new_norm<loctol) then
             finaliter=.true.
          end if
          if(have_maxiters) then
             if(i>=maxiters) then
                finaliter=.true.
             end if
          end if

       end if

       !       bef_time=mpi_wtime()
       call matrix_scale(Xk,-0.5_dp,1.0_dp)
       !       scale_time=scale_time+(mpi_wtime()-bef_time)
       !       call matrix_shift(Xk,1.0_dp)
       !       T2=0.5*(3*eye(117)-Xk)
       !       inv_sqrt_S=inv_sqrt_S*T2
       !       sqrt_S=T2*sqrt_S
       !       bef_time=mpi_wtime()
       call matrix_multiply(1.0_dp,inv_sqrt_S,Xk,0.0_dp,tmp_mat)
       !       prod_time=prod_time+(mpi_wtime()-bef_time)
       !       bef_time=mpi_wtime()
       call matrix_copy(tmp_mat,inv_sqrt_S)
       !       copy_time=copy_time+(mpi_wtime()-bef_time)
       !       bef_time=mpi_wtime()
       call matrix_multiply(1.0_dp,Xk,sqrt_S,0.0_dp,tmp_mat)
       !       prod_time=prod_time+(mpi_wtime()-bef_time)
       !       bef_time=mpi_wtime()
       call matrix_copy(tmp_mat,sqrt_S)
       !       copy_time=copy_time+(mpi_wtime()-bef_time)
       if(finaliter) then
          if(present(niter)) niter=i
          exit
       end if
    end do
    !    bef_time=mpi_wtime()
    call matrix_free(Xk)
    call matrix_free(tmp_mat)
    !    free_time=free_time+(mpi_wtime()-bef_time)
    !i
    !    bef_time=mpi_wtime()
    call matrix_scale(inv_sqrt_S,sqrt(lambda))
    call matrix_scale(sqrt_S,sqrt(lambda))
    !    scale_time=scale_time+(mpi_wtime()-bef_time)

    if(present(inv_sqrt_S_out)) call matrix_copy(inv_sqrt_S,inv_sqrt_S_out)
    if(present(sqrt_S_out)) call matrix_copy(sqrt_S,sqrt_S_out)

    !    bef_time=mpi_wtime()
    call matrix_free(sqrt_S)
    call matrix_free(inv_sqrt_S)
    !    free_time=free_time+(mpi_wtime()-bef_time)

    deallocate(conv_history,stat=ierr)
    call utils_dealloc_check('lowdin_transformation','conv_history',ierr)
    deallocate(coeffs,stat=ierr)
    call utils_dealloc_check('lowdin_transformation','coeffs',ierr)


    total_time=alloc_time+prod_time+copy_time+set_time+scale_time+free_time+trace_time

    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> alloc_time = ", alloc_time, "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> prod_time  = ", prod_time , "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> copy_time  = ", copy_time , "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> set_time   = ", set_time  , "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> scale_time = ", scale_time, "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> free_time  = ", free_time , "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> trace_time = ", trace_time, "<== TIME"
       if(debug_info_toggle) write(stdout,*) "smearing-> total_time = ", total_time, "<== TIME"
    end if

    if(present(status)) then
       status=loc_stat
    end if

    call timer_clock('lowdin_transformation',2)

    return

  contains

    function internal_frobenius_est(pos) result(two_norm_est)
      real(kind=dp), dimension(:),intent(in) :: pos
      real(kind=dp) :: two_norm_est
      real(kind=dp) :: lambda
      lambda=pos(1)
      two_norm_est = sqrt((((lambda**4)*trace_powers(4)) - (4.0_dp*(lambda**3)*trace_powers(3)) + &
           &(6.0_dp*(lambda**2)*trace_powers(2)) - (4.0_dp*lambda*trace_powers(1)) + trace_powers(0)) / &
           &(((lambda**2)*trace_powers(2)) - (2.0_dp*lambda*trace_powers(1)) + trace_powers(0)))
    end function internal_frobenius_est

    function internal_gerschgorin_est(pos) result(two_norm_est)
      real(kind=dp), dimension(:),intent(in) :: pos
      real(kind=dp) :: two_norm_est
      real(kind=dp) :: lambda
      two_norm_est = 0.0_DP ! @fixme
    end function internal_gerschgorin_est

    ! The following routines were written for other projects by me (Jolyon Aarons, ja531) but I allow for them
    ! to be used in ONETEP.

    subroutine internal_brents_method(whichfunc, start, middle, finish, tol, optimum, root, status)
      use comms, only: comms_bcast, pub_root_proc_id
      use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert
      implicit none

      !      interface
      !         real(kind=dp) function func(position)
      !           import dp
      !           real(kind=dp), dimension(:), intent(in) :: position
      !         end function func
      !      end interface
      integer, intent(in) :: whichfunc
      real(kind=dp), dimension(:), intent(in) :: start, middle, finish
      real(kind=dp), intent(in) :: tol
      real(kind=dp), dimension(:), intent(out) :: optimum
      real(kind=dp), intent(out) :: root
      integer, intent(out) :: status

      real(kind=dp), dimension(7) :: xvalues, fvalues

      integer, parameter :: startpoint=1,prev=2,prevprev=3,best=4,endpoint=5,latest=6,midpoint=7

      real(kind=dp) :: golden_ratio
      integer :: maxiters
      real(kind=dp) :: maxgrad
      real(kind=dp) :: curtol
      real(kind=dp) :: lim
      real(kind=dp) :: tmp(2), frac(2), diff(2), fdiff(2)
      real(kind=dp) :: steplength, quadstep, prevsteplength
      real(kind=dp) :: phi
      integer :: i
      logical :: reverse
      logical :: dogolden

      real(kind=dp), parameter :: numtol=0.0000001_dp

      real(kind=dp), dimension(:), allocatable :: direction
      real(kind=dp), dimension(:), allocatable :: initial

      integer :: ierr

      allocate(direction(size(start)),stat=ierr)
      call utils_alloc_check('brents_method','direction',ierr)

      allocate(initial(size(start)),stat=ierr)
      call utils_alloc_check('brents_method','initial',ierr)

      fvalues(:) = 0.0_dp
      xvalues(:) = 0.0_dp

      status=-1

      ! Calculate angle between the two halves of the space.
      phi=dot_product((middle-start)/sqrt(dot_product(middle-start,middle-start)),&
           &(middle-finish)/sqrt(dot_product(middle-finish,middle-finish)))
      call comms_bcast(pub_root_proc_id,phi)
      !      call utils_assert(.not.((abs(phi)-1.0_dp.ge.numtol).or.(abs(phi)-1.0_dp.le.-numtol)),&
      !           &"Error : The three initial points do not lie on a line to within specified numerical precision = ", numtol)
      if((abs(phi)-1.0_dp.ge.numtol).or.(abs(phi)-1.0_dp.le.-numtol)) then
         status=-2
         return
      end if


      maxgrad=100.0_dp
      maxiters=1000
      !    golden_ratio=(1.0_dp + sqrt(5.0_dp))/2.0_dp
      golden_ratio=1.0_dp-(2.0_dp/(1.0_dp + sqrt(5.0_dp))) !Conjugate GR

      xvalues(startpoint)=0.0_dp
      xvalues(endpoint)=1.0_dp


      if(xvalues(endpoint)<xvalues(startpoint)) then
         initial=finish
         direction=(start-finish)
         reverse=.true.
      else
         initial=start
         direction=-(start-finish)
         reverse=.false.
      end if

      call utils_assert(all(abs(direction) > numtol), 'Error in internal_brents_method:&
           & interval has to be finite.')

      xvalues(best)=sqrt(dot_product(middle-initial,middle-initial))/sqrt(dot_product(direction,direction))
      xvalues(prev)=xvalues(best)
      xvalues(prevprev)=xvalues(prev)


      if(whichfunc==frobfunc) then
         fvalues(best)=internal_frobenius_est(initial+direction*xvalues(best))
      else if(whichfunc==gersfunc) then
         fvalues(best)=internal_gerschgorin_est(initial+direction*xvalues(best))
      end if

      fvalues(prev)=fvalues(best)
      fvalues(prevprev)=fvalues(prev)

      quadstep=0.0_dp
      steplength=0.0_dp

      linmin: do i=1,maxiters
         xvalues(midpoint) = 0.5_dp * (xvalues(startpoint) + xvalues(endpoint))
         curtol = tol*abs(xvalues(best)) + epsilon(xvalues(best))

         ! Stopping condition in Brent's method is that startpoint and endpoint are 2*x*tol apart
         if(abs((xvalues(best)-xvalues(midpoint))) <= (2.0_dp*curtol-0.5_dp*(xvalues(endpoint)-xvalues(startpoint)))) then
            status=0
            exit linmin
         end if

         dogolden=.true.
         if(abs(steplength)>curtol) then ! Not on first iteration
            ! Take quadratic step.
            diff(1)=xvalues(best)-xvalues(prev)
            diff(2)=xvalues(best)-xvalues(prevprev)
            fdiff(1)=fvalues(best)-fvalues(prev)
            fdiff(2)=fvalues(best)-fvalues(prevprev)
            tmp(1)=diff(1)*fdiff(2)
            tmp(2)=diff(2)*fdiff(1)

            frac(1)=tmp(1)*diff(1) - tmp(2)*diff(2)
            frac(2)=tmp(1)-tmp(2)


            if(frac(2)>0.0_dp) frac(1)=-frac(1)
            frac(2)=abs(frac(2))
            prevsteplength=steplength
            steplength=quadstep
            if((abs(frac(1)) >= abs(frac(2)*prevsteplength)).or.(0.5_dp*frac(1) >= (xvalues(startpoint)-xvalues(best))*frac(2)).or.&
                 &(0.5_dp*frac(1) <= (xvalues(endpoint)-xvalues(best))*frac(2))) then
               dogolden=.true.
            else
               dogolden=.false.
               quadstep=0.5_dp*(frac(1)/frac(2))
               xvalues(latest)=xvalues(best)+quadstep
               if((xvalues(latest)-xvalues(startpoint)<2.0_dp*curtol).or.(xvalues(endpoint)-xvalues(latest)<2.0_dp*curtol)) &
                    &quadstep=sign(curtol,xvalues(midpoint)-xvalues(best))
            end if
         end if

         if(dogolden) then
            if(xvalues(best)>=xvalues(midpoint)) then
               steplength=golden_ratio*(xvalues(startpoint)-xvalues(best))
            else
               steplength=golden_ratio*(xvalues(endpoint)-xvalues(best))
            end if
         end if

         if(abs(steplength)>=curtol) then
            xvalues(latest)=xvalues(best)+steplength
         else
            xvalues(latest)=xvalues(best)+sign(curtol,steplength)
         end if


         if(whichfunc==frobfunc) then
            fvalues(latest)=internal_frobenius_est(initial+direction*xvalues(latest))
         else if(whichfunc==gersfunc) then
            fvalues(latest)=internal_gerschgorin_est(initial+direction*xvalues(latest))
         end if

         if(fvalues(latest)<=fvalues(best)) then
            if(xvalues(latest)>=xvalues(best)) then
               xvalues(startpoint)=xvalues(best)
            else
               xvalues(endpoint)=xvalues(best)
            end if
            xvalues(prevprev)=xvalues(prev)
            fvalues(prevprev)=fvalues(prev)
            xvalues(prev)=xvalues(best)
            fvalues(prev)=fvalues(best)
            xvalues(best)=xvalues(latest)
            fvalues(best)=fvalues(latest)
         else
            if(xvalues(latest)<xvalues(best)) then
               xvalues(startpoint)=xvalues(latest)
            else
               xvalues(endpoint)=xvalues(latest)
            end if
            if((fvalues(latest)<=fvalues(prev)).or.(xvalues(prev)==xvalues(best))) then
               xvalues(prevprev)=xvalues(prev)
               fvalues(prevprev)=fvalues(prev)
               xvalues(prev)=xvalues(latest)
               fvalues(prev)=fvalues(latest)
            else if((fvalues(latest)<=fvalues(prevprev)).or.(xvalues(prevprev)==xvalues(best)).or.&
                 &(xvalues(prevprev)==xvalues(prev))) then
               xvalues(prevprev)=xvalues(latest)
               fvalues(prevprev)=fvalues(latest)
            end if
         end if
      end do linmin

      optimum = initial+direction*xvalues(best)
      root = fvalues(best)

      deallocate(initial,stat=ierr)
      call utils_dealloc_check('brents_method','initial',ierr)
      deallocate(direction,stat=ierr)
      call utils_dealloc_check('brents_method','direction',ierr)

    end subroutine internal_brents_method


    subroutine internal_bracket_minimum(whichfunc, start, finish, status, xvals, fvals, middle)
      use utils, only: utils_alloc_check, utils_dealloc_check
      implicit none
      !      interface
      !         real(kind=dp) function func(position)
      !           import dp
      !           real(kind=dp), dimension(:), intent(in) :: position
      !         end function func
      !      end interface
      integer :: whichfunc
      real(kind=dp), dimension(:), intent(inout) :: start, finish
      integer, intent(out) :: status

      real(kind=dp), dimension(3), optional, intent(out) :: xvals, fvals
      real(kind=dp), dimension(:), optional, intent(out) :: middle

      real(kind=dp), dimension(4) :: xvalues, fvalues

      real(kind=dp) :: golden_ratio
      integer :: maxiters
      real(kind=dp) :: maxgrad

      real(kind=dp) :: lim
      real(kind=dp) :: tmp(2)
      integer :: i
      integer :: ierr
      logical :: reverse

      real(kind=dp), dimension(:), allocatable :: direction
      real(kind=dp), dimension(:), allocatable :: initial

      allocate(direction(size(start)),stat=ierr)
      call utils_alloc_check('bracket_minimum','direction',ierr)

      allocate(initial(size(start)),stat=ierr)
      call utils_alloc_check('bracket_minimum','initial',ierr)

      status=-1

      if(whichfunc==frobfunc) then
         fvalues(1)=internal_frobenius_est(start)
         fvalues(2)=internal_frobenius_est(finish)
      else if(whichfunc==gersfunc) then
         fvalues(1)=internal_gerschgorin_est(start)
         fvalues(2)=internal_gerschgorin_est(finish)
      end if

      maxgrad=100.0_dp
      maxiters=1000
      golden_ratio=(1.0_dp + sqrt(5.0_dp))/2.0_dp

      if(fvalues(2)>fvalues(1)) then
         fvalues(4)=fvalues(2)
         fvalues(2)=fvalues(1)
         fvalues(1)=fvalues(4)
         initial=finish
         direction=start-finish
         reverse=.true.
      else
         initial=start
         direction=finish-start
         reverse=.false.
      end if

      xvalues(1)=0.0_dp
      xvalues(2)=1.0_dp
      xvalues(3)=1.0_dp+golden_ratio

      if(whichfunc==frobfunc) then
         fvalues(3)=internal_frobenius_est(initial+direction*xvalues(3))
      else if(whichfunc==gersfunc) then
         fvalues(3)=internal_gerschgorin_est(initial+direction*xvalues(3))
      end if

      golden_section_loop: do i=1,maxiters
         if(fvalues(2)<fvalues(3)) then
            status=0
            exit golden_section_loop
         end if
         tmp(1)=(xvalues(2)-xvalues(1))*(fvalues(2)-fvalues(3))
         tmp(2)=(xvalues(2)-xvalues(3))*(fvalues(2)-fvalues(1))
         xvalues(4)=xvalues(2)-((xvalues(2)-xvalues(3))*tmp(2)-(xvalues(2)-xvalues(1))*tmp(1))&
              &/(2.0_dp*sign(max(abs(tmp(2)-tmp(1)),tiny(1.0_dp)),tmp(2)-tmp(1)))
         lim=xvalues(2)+maxgrad*(xvalues(3)-xvalues(2))

         if((xvalues(2)-xvalues(4))*(xvalues(4)-xvalues(3)) > 0.0_dp) then
            if(whichfunc==frobfunc) then
               fvalues(4)=internal_frobenius_est(initial+direction*xvalues(4))
            else if(whichfunc==gersfunc) then
               fvalues(4)=internal_gerschgorin_est(initial+direction*xvalues(4))
            end if
            if(fvalues(4)<fvalues(3)) then
               xvalues(1)=xvalues(2)
               fvalues(1)=fvalues(2)
               xvalues(2)=xvalues(4)
               fvalues(2)=fvalues(4)
               status=0
               exit golden_section_loop
            else if(fvalues(4)>fvalues(2)) then
               xvalues(3)=xvalues(4)
               fvalues(3)=fvalues(4)
               status=0
               exit golden_section_loop
            end if
            xvalues(4)=xvalues(3)+golden_ratio*(xvalues(3)-xvalues(2))
            if(whichfunc==frobfunc) then
               fvalues(4)=internal_frobenius_est(initial+direction*xvalues(4))
            else if(whichfunc==gersfunc) then
               fvalues(4)=internal_gerschgorin_est(initial+direction*xvalues(4))
            end if
         else if((xvalues(3)-xvalues(4))*(xvalues(4)-lim) > 0.0_dp) then
            if(whichfunc==frobfunc) then
               fvalues(4)=internal_frobenius_est(initial+direction*xvalues(4))
            else if(whichfunc==gersfunc) then
               fvalues(4)=internal_gerschgorin_est(initial+direction*xvalues(4))
            end if
            if(fvalues(4)<fvalues(3)) then
               xvalues(2:4)=eoshift(xvalues(2:4),1)
               fvalues(2:4)=eoshift(fvalues(2:4),1)
               xvalues(4)=xvalues(3)+golden_ratio*(xvalues(3)-xvalues(2))
               if(whichfunc==frobfunc) then
                  fvalues(4)=internal_frobenius_est(initial+direction*xvalues(4))
               else if(whichfunc==gersfunc) then
                  fvalues(4)=internal_gerschgorin_est(initial+direction*xvalues(4))
               end if
            end if
         else if((xvalues(4)-lim)*(lim-xvalues(3)) >= 0.0_dp) then
            xvalues(4)=lim
            if(whichfunc==frobfunc) then
               fvalues(4)=internal_frobenius_est(initial+direction*xvalues(4))
            else if(whichfunc==gersfunc) then
               fvalues(4)=internal_gerschgorin_est(initial+direction*xvalues(4))
            end if
         else
            xvalues(4)=xvalues(3)+golden_ratio*(xvalues(3)-xvalues(2))
            if(whichfunc==frobfunc) then
               fvalues(4)=internal_frobenius_est(initial+direction*xvalues(4))
            else if(whichfunc==gersfunc) then
               fvalues(4)=internal_gerschgorin_est(initial+direction*xvalues(4))
            end if
         end if
         xvalues=eoshift(xvalues,1)
         fvalues=eoshift(fvalues,1)
      end do golden_section_loop

      if(reverse) then
         if(pub_on_root) then
            if(debug_info_toggle) write(stdout,*) "smearing-> reverse",i, "<== SMEAR"
         end if
         start=initial+direction*xvalues(3)
         finish=initial+direction*xvalues(1)
      else
         start=initial+direction*xvalues(1)
         finish=initial+direction*xvalues(3)
         if(pub_on_root) then
            if(debug_info_toggle) write(stdout,*) "smearing-> not reverse",i, "<== SMEAR"
         end if
      end if

      if(present(xvals)) then
         xvals=xvalues(1:3)
      end if
      if(present(xvals)) then
         fvals=fvalues(1:3)
      end if
      if(present(middle)) then
         middle=initial+direction*xvalues(2)
      end if

      deallocate(initial,stat=ierr)
      call utils_dealloc_check('bracket_minimum','initial',ierr)
      deallocate(direction,stat=ierr)
      call utils_dealloc_check('bracket_minimum','direction',ierr)

    end subroutine internal_bracket_minimum
  end subroutine lowdin_transformation

  !============================================================================!
  ! This routine applies the fermi operator to a modified Hamiltonian, such    !
  ! that H -> B*(L'*H*L - I*u), where B is thermodynamic beta... of the        !
  ! electrons, u is the Fermi energy and L is the inverse of a symmetric       !
  ! matrix decomposition of the overlap matrix (such as Cholesky or Lowdin).   !
  ! The result is a density matrix in the orthogonal space. To achieve the     !
  ! density matrix in the original space, the result must be transformed back  !
  ! using L. This routine performs quite well conditioned operations, provided !
  ! that certain parameters are set sensibly. For instance D is the precision  !
  ! of the eventual result, i.e. an error of 10^-D, eigtol is the tolerance of !
  ! the eigenvalue error in H, but only rough estimates of the eigenrange are  !
  ! required, so it is safe to set this fairly coarsly. P and lbar are the     !
  ! most important parameters. P is the total number of terms to factor the    !
  ! operator expansion into. Formally this is exact for any P, but as a        !
  ! Chebyshev approximation is used to the inverse of the subspace matrices,   !
  ! it is important that P is large (>100), as this is only OK in the limit of !
  ! large P. Secondly, lbar controls the lowest l (i.e. the matrix             !
  ! contributions, from 1 to P are labeled l) matrix for which the Chebyshev   !
  ! approximation is used, beneath this l, a Newton-Schultz-Hotelling approach !
  ! is used, with extrapolation from the next highest l matrix. for the most   !
  ! part (i.e. lbar > 2) this will give the same answers for any lbar, but     !
  ! this parameter most strongly defines how many matrix-matrix products are   !
  ! performed and hence the performance of the routine. With P=1000, this has  !
  ! been empirically determined to be usually in the range of 6 to 15.
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The input scaled Hamiltonian matrix.             !
  !   rho          (output) : The corresponding density matrix.                !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Jan 2014,                                        !
  !----------------------------------------------------------------------------!
  ! This routine scales as O(N^(1/3)) in matrix-matrix multiplications. This   !
  ! Means that if the bandwidth of the Hamiltonian scales linearly with system !
  ! size, as when a density kernel cut-off or thresholding is used, that the   !
  ! whole order of the operation is O(N^(4/3)), which compares favourably with !
  ! diagonalisation -> O(N^3)... but not quite linear scaling.                 !
  !============================================================================!
  subroutine fermi_operator_chebyshev(H,rho,matmuls)
    use comms, only : pub_on_root,comms_barrier
    use utils, only : utils_alloc_check, utils_dealloc_check
    implicit none
    type(smearing_matrix), intent(inout)  :: H ! Inout for dense_create_copy in matrix_allocate (in)
    type(smearing_matrix), intent(inout) :: rho
    integer, optional, intent(out) :: matmuls

    integer :: N
    real(kind=dp) :: z0,chi
    type(smearing_matrix) :: H_p, H_2p, X
    integer :: P, D, lbar
    integer :: ml
    real(kind=dp) :: min_eig, max_eig
    real(kind=dp), dimension(:), allocatable :: coeffs
    real(kind=dp), dimension(:,:), allocatable :: dccoeffs
    integer :: li
    complex(kind=dp) :: b
    real(kind=dp) :: sgnb
    integer :: di0, i
    real(kind=dp), dimension(:), allocatable :: c
    integer :: layers
    integer :: jh,jl, jlold, jhold
    integer :: j,k
    integer, dimension(:), allocatable :: num, parent, start
    integer :: calclayerstart
    type(smearing_matrix), dimension(:), allocatable :: exes
    type(smearing_matrix), dimension(:), allocatable :: accumat
    type(smearing_matrix), dimension(:), allocatable :: NLcheb_acc
    integer :: l
    logical :: nans
    integer :: matmuls_divconq, matmuls_shn
    type(smearing_matrix) :: NLcheb, NLprop
    type(smearing_matrix) :: T0, T1, T
    type(smearing_matrix) :: tmp, tmp2
    type(smearing_matrix) :: rhobcheb, invMcheb
    type(smearing_matrix) :: compcheb
    real(kind=dp) :: conv
    real(kind=dp) :: conv_test, conv_test_old
    real(kind=dp) :: infone
    integer :: mn
    integer :: nmats
    real(kind=dp) :: nu
    integer :: storing
    real(kind=dp) :: eigtol
    integer :: ierr

    N=matrix_dimension(H)

    P=1000
    eigtol=1e-5_dp
    D=9
    lbar=15

    call matrix_allocate(tmp,H)
    call matrix_scale(tmp,0.0_dp)
    call matrix_scale(tmp,0.0_dp,1.0_dp)

    call extremum_eigs(H,tmp,eigtol,min_eig,max_eig)

    z0=(exp(max_eig/(2.0_dp*real(P,dp)))+exp(min_eig/(2.0_dp*real(P,dp))))/2.0_dp
    chi=(exp(max_eig/(2.0_dp*real(P,dp)))-exp(min_eig/(2.0_dp*real(P,dp))))/2.0_dp

    call matrix_allocate(H_p,H)
    call matrix_allocate(H_2p,H)

    call matrix_copy(H,tmp)
    call matrix_scale(tmp,1.0_dp/real(2*P,dp))
    call mat_exp_higham(tmp,H_2p)

    !    call matrix_copy(H_2p,tmp)
    !    call matrix_multiply(1.0_dp, H_2p, tmp, 0.0_dp, H_p)

    call matrix_copy(H,tmp)
    call matrix_scale(tmp,1.0_dp/real(P,dp))
    !    call matrix_allocate(H_2p,H)
    call mat_exp_higham(tmp,H_p)

    call matrix_allocate(X,H)

    call matrix_copy(H_2p,X)
    call matrix_scale(X,1.0_dp/chi,-z0/chi)

    ml=floor(((1.0_dp/2.0_dp)+(((max_eig)*real(D,dp)*log(10.0_dp))/(pi*real(2*lbar-1,dp)))))
    layers=floor(log10(real(ml,dp))/log10(2.0_dp))-2

    layers=max(layers,1)


    allocate(coeffs(1:(ml+1)*(layers+1)),stat=ierr)
    call utils_alloc_check('fermi_operator_chebyshev','coeffs',ierr)
    allocate(c(1:ml+1),stat=ierr)
    call utils_alloc_check('fermi_operator_chebyshev','c',ierr)
    allocate(dccoeffs(P-lbar+1,size(coeffs)),stat=ierr)
    call utils_alloc_check('fermi_operator_chebyshev','dccoeffs',ierr)

    allocate(num(1:2**(layers+1)-1),stat=ierr)
    call utils_alloc_check('fermi_operator_chebyshev','num',ierr)
    allocate(parent(1:2**(layers+1)-1),stat=ierr)
    call utils_alloc_check('fermi_operator_chebyshev','parent',ierr)
    allocate(start(1:2**(layers+1)-1),stat=ierr)
    call utils_alloc_check('fermi_operator_chebyshev','start',ierr)

    do li=lbar,P
       b=(exp(cmplx(0.0_dp,pi*(real(2*li-1,dp)/real(2*P,dp)),dp))-z0)/chi
       sgnb=sign(1.0_dp,real(b,dp))

       c=0.0_dp
       do i=0,ml
          if(i==0) then
             di0=1
          else
             di0=0
          end if
          c(i+1)=(real(di0-2,dp)/((chi**2)*aimag(b)))*aimag((sgnb/sqrt((b**2)-1.0_dp))*((b-(sqrt((b**2)-1.0_dp)*sgnb))**i))
       end do

       ! Construct binary tree
       jh=0
       jl=0
       jlold=1
       jhold=1

       num=0
       num(1)=(ml+1)*2

       parent=0

       start=0
       start(1)=1
       do i=0,layers
          jl=jl+2**(max(0,i-1))
          jh=jh+2**(i)
          do j=jl,jh
             parent(j)=jlold+(floor(real((j-jl+2),dp)/2.0_dp)-1)
             if(mod(j,2)==0) then
                num(j)=floor(real(num(parent(j)),dp)/2.0_dp)
             else
                num(j)=floor(real((num(parent(j))-1),dp)/2.0_dp)
             end if
          end do
          jlold=jl
          jhold=jh
       end do

       do i=2,size(num)
          start(i)=start(i-1)+num(i-1)+1
       end do

       coeffs=0
       coeffs(1:ml+1)=c
       jl=0
       jh=0
       jlold=1
       jhold=ml+1
       do i=2,size(num)
          jl=jhold+1
          jh=jl+num(i)
          if(mod(i,2)==0) then
             do j=jh,jl,-1
                k=j-jl
                coeffs(j)=coeffs(2*k+start(parent(i)))
             end do
          else
             coeffs(jh)=2*coeffs(2*num(i)+1+start(parent(i)))
             do j=jh-1,jl+1,-1
                k=j-jl
                coeffs(j)=2*coeffs(2*k+1+start(parent(i)))-coeffs(k+1+jl)
             end do
             coeffs(jl)=coeffs(start(parent(i))+1)-(coeffs(jl+1)/2)
          end if
          jlold=jl
          jhold=jh

       end do
       !dccoeffs=0.0_dp
       dccoeffs(li-lbar+1,1:1+size(coeffs)-start(size(num)-2**layers+1)) = coeffs(start(size(num)-2**layers+1):size(coeffs))
       !       stop
    end do
    deallocate(c,stat=ierr)
    call utils_dealloc_check('fermi_operator_chebyshev','c',ierr)
    deallocate(coeffs,stat=ierr)
    call utils_dealloc_check('fermi_operator_chebyshev','coeffs',ierr)
    deallocate(parent,stat=ierr)
    call utils_dealloc_check('fermi_operator_chebyshev','parent',ierr)

    calclayerstart=start(size(num)-2**layers+1)
    allocate(exes(1:layers+1),stat=ierr)
    call utils_alloc_check('fermi_operator_chebyshev','exes',ierr)
    do i = 1,layers+1
       call matrix_allocate(exes(i),H)
       call matrix_scale(exes(i),0.0_dp)
    end do
    matmuls_divconq=0


    call comms_barrier()
    call matrix_copy(X,exes(1))
    call matrix_free(X)

    do i=2,layers+1
       matmuls_divconq=matmuls_divconq+1
       call matrix_multiply(1.0_dp,exes(i-1),exes(i-1),0.0_dp,exes(i))
       call matrix_scale(exes(i),2.0_dp,-1.0_dp)
    end do

    allocate(accumat(1:2**layers),stat=ierr)
    call utils_alloc_check('fermi_operator_chebyshev','accumat',ierr)
    do i=1,2**layers
       call matrix_allocate(accumat(i),H)
       call matrix_scale(accumat(i),0.0_dp)
    end do

    allocate(NLcheb_acc(1:2**layers),stat=ierr)
    call utils_alloc_check('fermi_operator_chebyshev','NLcheb_acc',ierr)
    do i=1,2**layers
       call matrix_allocate(NLcheb_acc(i),H)
       call matrix_scale(NLcheb_acc(i),0.0_dp)
    end do

    call matrix_allocate(T0,H)
    call matrix_allocate(T1,H)
    call matrix_allocate(T,H)

    do i=0,maxval(num(size(num)-2**layers+1:size(num)))
       if(i==0) then
          call matrix_scale(T0,0.0_dp)
          call matrix_scale(T0,0.0_dp,1.0_dp)

          call matrix_copy(T0,T)
       else if(i==1) then

          call matrix_copy(exes(layers+1),T1)

          call matrix_copy(T1,T)
       else

          call matrix_copy(T0,T)
          call matrix_multiply(2.0_dp,exes(layers+1),T1,-1.0_dp,T)
          matmuls_divconq=matmuls_divconq+1

          call matrix_copy(T1,T0)

          call matrix_copy(T,T1)
       end if
       do j=1,(2**layers)
          if(i<=num(size(num)-2**layers+j)) then

             do li=lbar,P
                call matrix_axpy(accumat(j),T,dccoeffs(li-lbar+1,start(size(num)-2**layers+j)-start(size(num)-2**layers+1)+1+i))
             end do
             call matrix_axpy(NLcheb_acc(j),T,dccoeffs(1,start(size(num)-2**layers+j)-start(size(num)-2**layers+1)+1+i))

          end if
       end do
    end do

    deallocate(dccoeffs,stat=ierr)
    call utils_dealloc_check('fermi_operator_chebyshev','decoeffs',ierr)
    call matrix_free(T0)
    call matrix_free(T1)
    call matrix_free(T)

    !    call matrix_allocate(compcheb,H)
    !    call matrix_scale(compcheb,0.0_dp)
    k=0
    do i=1,layers
       nmats=2**i
       do j=1,2**(layers-i)
          k=k+1
          matmuls_divconq=matmuls_divconq+1
          call matrix_multiply(1.0_dp,exes(layers-i+1),accumat(nmats*(j-1)+1+(2**(i-1))),1.0_dp,accumat(nmats*(j-1)+1))
          call matrix_multiply(1.0_dp,exes(layers-i+1),NLcheb_acc(nmats*(j-1)+1+(2**(i-1))),1.0_dp,NLcheb_acc(nmats*(j-1)+1))
          matmuls_divconq=matmuls_divconq+1
       end do
    end do
    call matrix_allocate(NLcheb,H)


    call matrix_copy(accumat(1),NLcheb)

    do i = 1,layers+1
       call matrix_free(exes(i))
    end do
    deallocate(exes,stat=ierr)
    call utils_dealloc_check('fermi_operator_chebyshev','exes',ierr)

    do i=1,2**layers
       call matrix_free(accumat(i))
    end do
    deallocate(accumat,stat=ierr)
    call utils_dealloc_check('fermi_operator_chebyshev','accumat',ierr)

    call matrix_allocate(invMcheb,H)

    call matrix_copy(H_p,invMcheb)
    call matrix_scale(invMcheb,0.5_dp,-0.5_dp)
    call matrix_allocate(rhobcheb,H)
    call matrix_multiply(-1.0_dp,invMcheb,NLcheb,0.0_dp,rhobcheb)
    call matrix_scale(rhobcheb,1.0_dp,0.5*real(P-lbar+1,dp))

    call matrix_copy(NLcheb_acc(1),NLcheb)

    do i=1,2**layers
       call matrix_free(NLcheb_acc(i))
    end do
    deallocate(NLcheb_acc,stat=ierr)
    call utils_dealloc_check('fermi_operator_chebyshev','NLcheb_acc',ierr)

    call matrix_allocate(tmp2,H)
    call matrix_allocate(NLprop,H)

    storing=2**layers + layers+1
    matmuls_shn=0
    do l=lbar-1,1,-1
       nu=cos(pi/real(P,dp))+sin(pi/real(P,dp))*tan(pi*real(2*l-1,dp)/(2.0_dp*real(P,dp)))
       !call matrix_copy(NLcheb,tmp)
       !call matrix_multiply((1.0_dp-(1.0_dp/nu)),NLcheb,H_p,(1.0_dp-(1.0_dp/nu)),tmp)
       !NLcheb=(NLcheb/nu)*(eye(N)+tmp+mpower(tmp,2))
       call matrix_scale(NLcheb,1.0_dp/nu)
       conv=(8.0_dp*real(l,dp))/real((1+2*l)**2,dp)
       mn=nint(((1.0_dp)/(log(2.0_dp)))*log((log(real(1-conv,dp))-real(D,dp)*log(10.0_dp))/log(conv)))

       call matrix_copy(H_p,NLprop)

       call matrix_axpy(NLprop,H_2p,-2.0_dp*cos(pi*real(2*l-1,dp)/(2.0_dp*real(P,dp))))
       call matrix_scale(NLprop,1.0_dp,1.0_dp)


       nans=.false.

       i=0
       conv_test_old=huge(1.0_dp)
       do
          i=i+1
          !        tmp=(2.0_dp*NLcheb - (NLcheb*NLprop*NLcheb))

          call matrix_multiply(1.0_dp,Nlprop,NLcheb,0.0_dp,tmp)

          call matrix_copy(Nlcheb,tmp2)
          call matrix_multiply(-1.0_dp,Nlcheb,tmp,2.0_dp,tmp2)


          matmuls_shn=matmuls_shn+2

          call matrix_multiply(1.0_dp,NLprop,tmp2,0.0_dp,tmp)
          call matrix_scale(tmp,1.0_dp,-1.0_dp)

          conv_test=abs(matrix_trace(tmp))

          if(conv_test>conv_test_old.or.conv_test/=conv_test) then ! NaN is not equal to itself
             nans=.true.
             exit
          end if

          call matrix_copy(tmp2,NLcheb)
          if(i>=mn) then
             if(conv_test<=10.0_dp**(-D)) then
                exit
             end if
          end if

          conv_test_old=conv_test
       end do
       if(nans) then
          if(pub_on_root) then
             if(debug_info_toggle) write(stdout,*) "smearing-> Warning: reinitialising inverse in Newton-&
                  &Shultz-Hotelling, low-l Chebyshev matrix contribution.", "<== SMEAR"
          end if
          !return
          mn=nint(log(max_eig-min_eig)+log(real(D,dp)*log(10.0_dp))+log(real(N,dp))/((pi**2)*real(2*l-1,dp)**2))
          !        NLcheb=transpose(NLprop)*(1.0_dp/(norm(NLprop,1)*norm(NLprop,inf)))

          call matrix_transpose(NLprop,NLcheb)

          infone=1.0_dp/(matrix_induced_norm(NLprop)*matrix_induced_norm(NLcheb)) ! one over the 1 norm * inf norm

          call matrix_scale(NLcheb,infone)

          !for i=1:mn+10
          i=0
          do
             i=i+1
             !           tmp=(2.0_dp*NLcheb - (NLcheb*NLprop*NLcheb))
             call matrix_multiply(1.0_dp,Nlprop,NLcheb,0.0_dp,tmp)

             call matrix_copy(Nlcheb,tmp2)
             call matrix_multiply(-1.0_dp,Nlcheb,tmp,2.0_dp,tmp2)
             matmuls_shn=matmuls_shn+2

             call matrix_copy(tmp2,NLcheb)
             if(i>=mn) then
                call matrix_multiply(1.0_dp,NLprop,tmp2,0.0_dp,tmp)
                call matrix_scale(tmp,1.0_dp,-1.0_dp)
                if(abs(matrix_trace(tmp))<=10.0_dp**(-D)) then
                   exit
                end if
             end if
          end do
       end if

       call matrix_copy(NLcheb,invMcheb)

       call matrix_multiply(1.0_dp,H_p,NLcheb,-1.0_dp,invMcheb)
       call matrix_scale(invMcheb,-0.5_dp,0.5_dp)
       call matrix_axpy(rhobcheb,invMcheb,1.0_dp)
    end do

    deallocate(num,stat=ierr)
    call utils_dealloc_check('fermi_operator_chebyshev','num',ierr)
    deallocate(start,stat=ierr)
    call utils_dealloc_check('fermi_operator_chebyshev','start',ierr)

    call matrix_free(tmp)
    call matrix_free(tmp2)
    call matrix_free(NLcheb)
    call matrix_free(invMcheb)
    call matrix_free(NLprop)
    call matrix_free(H_p)
    call matrix_free(H_2p)

    call matrix_scale(rhobcheb,1.0_dp/real(P,dp))

    call matrix_copy(rhobcheb,rho)
    call matrix_free(rhobcheb)

    if(present(matmuls)) then
       matmuls=matmuls_shn+matmuls_divconq
    end if
    !    matmuls_shn
    !    matmuls_divconq
    !    dkrc=inv(chol(overlap))*rhobcheb*transpose(inv(chol(overlap)))
    !    diff=abs(trace(dkrc)-trace(denskern))

  end subroutine fermi_operator_chebyshev

  !============================================================================!
  ! Condition number of a matrix!                                              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   A            (input)  : The input matrix.                                !
  !   n            (input)  : The norm order of the condition number
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, May 2014,                                        !
  !============================================================================!
  function condition_number(A,n) result(cond)
    implicit none
    type(smearing_matrix), intent(inout) :: A
    integer, optional, intent(in) :: n ! Optional norm order
    integer :: m
    real(kind=dp) :: cond
    type(smearing_matrix) :: A_inv
    if(present(n)) then
       m=n
    else
       m=2
    end if
    call matrix_allocate(A_inv,A)
    call matrix_copy(A,A_inv)
    call matrix_invert(A_inv)
    cond = matrix_norm(A_inv,m)*matrix_norm(A,m)
    call matrix_free(A_inv)
  end function condition_number


  !============================================================================!
  ! This routine applies the fermi operator to a modified Hamiltonian, such    !
  ! that H -> B*(L'*H*L - I*u), where B is thermodynamic beta... of the        !
  ! electrons, u is the Fermi energy and L is the inverse of a symmetric       !
  ! matrix decomposition of the overlap matrix (such as Cholesky or Lowdin).   !
  ! The result is a density matrix in the orthogonal space. To achieve the     !
  ! density matrix in the original space, the result must be transformed back  !
  ! using L.                                                                   !
  ! This routine constructs two projection operators using a shifted matrix    !
  ! sign function and uses them to construct three eigenvalue subspaces; the   !
  ! first where the eigenvalues are under some -threshold and give 0 to        !
  ! numerical precision, the second, where the eigenvalues are above a         !
  ! +threshold and give 1 to numerical precision (via logistic function) and   !
  ! those that are between -threshold and +threshold, which are dealt with     !
  ! separately, via a series expansion.
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The input scaled Hamiltonian matrix.             !
  !   rho          (output) : The corresponding density matrix.                !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, May 2014,                                        !
  !----------------------------------------------------------------------------!
  ! This routine scales as O(log(N)) in matrix-matrix multiplications.         !
  !============================================================================!
  subroutine fermi_operator_projection_contourint(H,S,invS,CCH,rho,orthogonal,&
       contravariant_kernel, avoid_solve, matmuls_out)
    use comms, only : pub_on_root


    use timer, only : timer_clock
    use utils, only : utils_dealloc_check
    implicit none
    type(smearing_matrix), intent(inout) :: H
    type(smearing_matrix), intent(inout) :: S
    type(smearing_matrix), intent(inout) :: invS
    type(smearing_matrix), intent(inout) :: CCH
    type(smearing_matrix), intent(inout) :: rho ! Inout to make operations at the end of the routine easier.
    logical, optional, intent(in)  :: orthogonal ! Ignore S
    logical, optional, intent(in)  :: contravariant_kernel
    logical, optional, intent(in)  :: avoid_solve
    integer, optional, intent(out) :: matmuls_out

    type(smearing_matrix) :: A
    type(smearing_matrix) :: inside_T
    type(smearing_matrix) :: inside_A, outside_A, outside_T
    type(smearing_matrix) :: CC_inside_A, CC_outside_A, CC_outside_A2, outside_A2

    type(smearing_matrix) :: CC_inside_T, CC_outside_T, CC_righthand_T, CC_T, CC_outside_T_inv
    type(smearing_matrix) :: CC_T_trans, CC_inside_TP
    type(smearing_matrix) :: T_trans
    type(smearing_matrix) :: S_trans
    type(smearing_matrix) :: projector
    type(smearing_matrix) :: eye

    real(kind=dp) :: tmp
    type(smearing_matrix) :: tmp1, tmp2
    real(kind=dp) :: cutoff, eigtol, signtol, skew
    integer :: ncontour
    real(kind=dp) :: min_abs_eig, min_eig, max_eig
    real(kind=dp), dimension(:), allocatable :: cheby_coeffs, cheby_coeffs_full, coeffs
    integer :: mmsi, matmuls_cheby, matmuls
    real(kind=dp) :: specwidth,frobnorm

    integer :: status
    real(kind=dp) :: shift
    integer :: i
    logical :: no_extremal_eigs
    logical :: use_gerschgorin
    integer :: ierr
    real(kind=dp) :: res
    real(kind=dp) :: cheby_threshold
    logical :: ignore_S
    integer :: outunit
    logical :: use_nsh
    logical :: contravariant_output
    logical :: no_solve

    real(kind=dp) :: abs_emax, maxeval

    call timer_clock('fermi_operator_projection_contourint',1)

    matmuls = 0

    ignore_S=.false.
    if(present(orthogonal)) then
       ignore_S=orthogonal
    end if

    contravariant_output=.true.
    if(present(contravariant_kernel)) then
       contravariant_output=contravariant_kernel
    end if

    no_solve=.true.
    if(present(avoid_solve)) then
       no_solve=avoid_solve
    end if

    cutoff=15.0_dp
    ncontour=32
    skew=0.5_dp

    cheby_threshold=1e-12_dp

    signtol=1e-7_dp

    ! ja531-> In the mapping from Fermi-Dirac to tanh, the argument A-> -H/2
    call matrix_allocate(A,H)
    call matrix_copy(H,A)
    call matrix_scale(A,-0.5_dp)


    call matrix_allocate(inside_A,A)


!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="A.mtx")
!    end if
!    call sparse_show_matrix(A%dataSPAM3,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
!    end if
!
!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="S.mtx")
!    end if
!    call sparse_show_matrix(S%dataSPAM3,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
!    end if

    call matrix_allocate(projector,A)

    call boxcar_projection_contourint(A,S,ncontour,0.0_dp,cutoff,skew,ignore_S,ignore_S,inside_A,projector,no_solve)


!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="inside_A.mtx")
!    end if
!    call sparse_show_matrix(inside_A%dataSPAM3,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
!    end if
!
!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="proj.mtx")
!    end if
!    call sparse_show_matrix(projector%dataSPAM3,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
!    end if


    call matrix_allocate(outside_A,A)

    call matrix_copy(A,outside_A)
    call matrix_axpy(outside_A,inside_A,-1.0_dp)



!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="outside_A.mtx")
!    end if
!    call sparse_show_matrix(outside_A%dataSPAM3,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
!    end if



    call matrix_allocate(inside_T,rho)
    if(.not.ignore_S) then
       call matrix_allocate(CC_inside_T,rho)
    end if

    if(.not.ignore_S) then
       call matrix_allocate(CC_inside_A,CCH)
       if(no_solve) then
          call matrix_multiply(1.0_dp,invS,inside_A,0.0_dp,CC_inside_A)
       else
          call matrix_solve(inside_A,S,CC_inside_A)
       end if
    end if


!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="CC_inside_A.mtx")
!    end if
!    call sparse_show_matrix(CC_inside_A%dataSPAM3,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
!    end if



    call matrix_scale(CC_inside_A,1.0_dp/cutoff)
    call get_tanh_chebyshev_coeffs(-cutoff,cutoff,cheby_threshold,coeffs)

    if(ignore_S) then
       call fast_cheby_resum(inside_A,coeffs,inside_T,matmuls_cheby)
    else
       call fast_cheby_resum(CC_inside_A,coeffs,CC_inside_T,matmuls_cheby)

       call matrix_free(CC_inside_A)
    end if

    deallocate(coeffs,stat=ierr)
    call utils_dealloc_check('fermi_operator_projection_contourint','cheby_coeffs',ierr)

    matmuls=matmuls+matmuls_cheby

    call matrix_allocate(CC_inside_TP,CC_inside_T)
    call matrix_scale(CC_inside_T,0.5_dp,0.5_dp)
    call matrix_multiply(1.0_dp,CC_inside_T,projector,0.0_dp,CC_inside_TP)
    call matrix_copy(CC_inside_TP,CC_inside_T)
    call matrix_scale(CC_inside_T,2.0_dp,-1.0_dp)
    call matrix_free(CC_inside_TP)


!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="CC_inside_T.mtx")
!    end if
!    call sparse_show_matrix(CC_inside_T%dataSPAM3,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
!    end if



    call matrix_allocate(outside_T,rho)

    if(.not.ignore_S) then
       call matrix_allocate(CC_outside_A,CCH)
       call matrix_allocate(CC_outside_T,rho)
    else
       call matrix_allocate(CC_outside_A,H)
       call matrix_allocate(CC_outside_T,rho)
    end if



    use_nsh=.false.

    if(use_nsh) then

       !vvv might need this vvv after all! vvv
       ! set so that all conduction eigs and projected out eigs -> 0.
       ! We are effectively setting the origin to be half-way between the
       ! non-projected valence eigs and the zero eigs, by doing this.
       if(ignore_S) then
          call matrix_scale(outside_A,1.0_dp,cutoff/2.0_dp)
       else
          call matrix_axpy(outside_A,S,cutoff/2.0_dp)
       end if

       shift=cutoff/2.0_dp


!       if(pub_on_root) then
!          outunit=utils_unit()
!          open(unit=outunit,file="outside_A_scale.mtx")
!       end if
!       call sparse_show_matrix(outside_A%dataSPAM3,outunit,matlab_format=.true.)
!       if(pub_on_root) then
!          close(outunit)
!       end if


       i=0
       do
          i=i+1
          if(ignore_S) then
             abs_emax=-1.0_dp
             call sign_nsh(outside_A,signtol,outside_T,abs_emax,mmsi,status)
          else
             if(no_solve) then
                call matrix_multiply(1.0_dp,invS,outside_A,0.0_dp,CC_outside_A)
             else
                call matrix_solve(outside_A,S,CC_outside_A)
             end if
!             if(pub_on_root) then
!                outunit=utils_unit()
!                open(unit=outunit,file="CC_outside_A.mtx")
!             end if
!             call sparse_show_matrix(CC_outside_A%dataSPAM3,outunit,matlab_format=.true.)
!             if(pub_on_root) then
!                close(outunit)
!             end if

             abs_emax=-1.0_dp
             call sign_nsh(CC_outside_A,signtol,CC_outside_T,abs_emax,mmsi,status)
          end if

          if(pub_on_root) then
             if(debug_info_toggle) then
                write(stdout,*) "smearing-> 1st Sign: ",mmsi
             end if
          end if
          matmuls=matmuls+mmsi
          if(status==0) then
             exit
          else
             !undo old shift
             if(ignore_S) then
                call matrix_scale(outside_A,1.0_dp,-shift)
             else
                call matrix_axpy(outside_A,S,-shift)
             end if

             call random_number(shift)
             shift=(shift-0.5_dp)/2.0_dp !*2.0_dp
             !          shift=shift*10.0_dp**(i-4)
             !          call matrix_scale(outside_A,1.0_dp,shift)

             ! shift: -0.25 -- 0.25

             shift=(shift*cutoff) + cutoff/2.0_dp

             if(ignore_S) then
                call matrix_scale(outside_A,1.0_dp,shift)
             else
                call matrix_axpy(outside_A,S,shift)
             end if


          end if
       end do

    ! Back to Fermi-Dirac from tanh...
    call matrix_scale(CC_inside_T,0.5_dp,0.5_dp)
    call matrix_scale(CC_outside_T,0.5_dp,0.5_dp)

    call matrix_allocate(CC_righthand_T,rho)
    call matrix_multiply(1.0_dp,CC_inside_T,CC_outside_T,1.0_dp,CC_righthand_T)

    else ! using contour integral

       if(ignore_S) then
          call matrix_copy(outside_A,CC_outside_A)
       else
          if(no_solve) then
             call matrix_multiply(1.0_dp,invS,outside_A,0.0_dp,CC_outside_A)
          else
             call matrix_solve(outside_A,S,CC_outside_A)
          end if
       end if

!       if(pub_on_root) then
!          outunit=utils_unit()
!          open(unit=outunit,file="CC_outside_A.mtx")
!       end if
!       call sparse_show_matrix(CC_outside_A%dataSPAM3,outunit,matlab_format=.true.)
!       if(pub_on_root) then
!          close(outunit)
       !       end if


       call matrix_allocate(CC_outside_A2,CC_outside_A)


       specwidth=matrix_norm(CC_outside_A2,2)

       specwidth=sqrt(specwidth)

       frobnorm=sqrt(matrix_trace(CC_outside_A2))/sqrt(real(matrix_dimension(CC_outside_A2),dp))


       call matrix_allocate(outside_A2,CCH)
       ! WARNING: THIS NUKES CC_OUTSIDE_A !!! MAKE SURE WE DON'T NEED IT!

       if(ignore_S) then
          call matrix_copy(A,CC_outside_A)
       else
          if(no_solve) then
             call matrix_multiply(1.0_dp,invS,A,0.0_dp,CC_outside_A)
          else
             call matrix_solve(A,S,CC_outside_A)
          end if
       end if


       call matrix_multiply(1.0_dp,CC_outside_A,CC_outside_A,0.0_dp,CC_outside_A2)

       if(ignore_S) then
          call matrix_allocate(eye,H)
          call matrix_scale(eye,0.0_dp,1.0_dp)
          call matrix_copy(CC_outside_A2,outside_A2)
          call extremum_eigs(outside_A2,eye,tol=1.0_dp,maxeval=maxeval)
          call matrix_free(eye)
       else
          call matrix_multiply(1.0_dp,S,CC_outside_A2,0.0_dp,outside_A2)
          call extremum_eigs(outside_A2,S,tol=1.0_dp,maxeval=maxeval)
       end if


       call matrix_free(outside_A2)

       maxeval=sqrt(maxeval)

       if(pub_on_root) then
          if(debug_info_toggle) then
             write(stdout,*) "Spectral width:",specwidth,frobnorm,maxeval
          end if
       end if

       specwidth=maxeval!frobnorm

       call matrix_free(CC_outside_A2)

       call boxcar_projection_contourint(outside_A,S,ncontour,cutoff + (specwidth/2.0_dp),cutoff + specwidth/2.0_dp, &
            & skew,ignore_S,ignore_S,inside_A,CC_outside_T,no_solve)


       call matrix_scale(projector,-1.0_dp,1.0_dp)
       call matrix_allocate(CC_inside_TP,CC_outside_T)
       call matrix_multiply(1.0_dp,CC_outside_T,projector,0.0_dp,CC_inside_TP)
       call matrix_copy(CC_inside_TP,CC_outside_T)
       call matrix_free(CC_inside_TP)


       ! Back to Fermi-Dirac from tanh...
       call matrix_scale(CC_inside_T,0.5_dp,0.5_dp)

       call matrix_allocate(CC_righthand_T,rho)
!       call matrix_multiply(1.0_dp,CC_inside_T,CC_outside_T,1.0_dp,CC_righthand_T)

    end if

    call matrix_free(projector)
    call matrix_free(A)


!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="CC_outside_T.mtx")
!    end if
!    call sparse_show_matrix(CC_outside_T%dataSPAM3,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
!    end if

    if(use_nsh) then

!       if(pub_on_root) then
!          outunit=utils_unit()
!          open(unit=outunit,file="CC_righthand_T.mtx")
!       end if
!       call sparse_show_matrix(CC_righthand_T%dataSPAM3,outunit,matlab_format=.true.)
!       if(pub_on_root) then
!          close(outunit)
!       end if



       if(ignore_S) then
          call matrix_scale(outside_A,1.0_dp,-cutoff)
       else
          call matrix_axpy(outside_A,S,-cutoff)
       end if


!       if(pub_on_root) then
!          outunit=utils_unit()
!          open(unit=outunit,file="outside_A_scale2.mtx")
!       end if
!       call sparse_show_matrix(outside_A%dataSPAM3,outunit,matlab_format=.true.)
!       if(pub_on_root) then
!          close(outunit)
!       end if

       shift=0.0_dp
       i=0
       do
          i=i+1
          if(ignore_S) then
             abs_emax=-1.0_dp
             call sign_nsh(outside_A,signtol,outside_T,abs_emax,mmsi,status)
          else

             if(no_solve) then
                call matrix_multiply(1.0_dp,invS,outside_A,0.0_dp,CC_outside_A)
             else
                call matrix_solve(outside_A,S,CC_outside_A)
             end if

!             if(pub_on_root) then
!                outunit=utils_unit()
!                open(unit=outunit,file="CC_outside_A2.mtx")
!             end if
!             call sparse_show_matrix(CC_outside_A%dataSPAM3,outunit,matlab_format=.true.)
!             if(pub_on_root) then
!                close(outunit)
!             end if

             abs_emax=-1.0_dp
             call sign_nsh(CC_outside_A,signtol,CC_outside_T,abs_emax,mmsi,status)
          end if

          if(pub_on_root) then
             if(debug_info_toggle) then
                write(stdout,*) "smearing-> 1st Sign: ",mmsi
             end if
          end if
          matmuls=matmuls+mmsi
          if(status==0) then
             exit
          else
             call random_number(shift)
             shift=(shift-1.0_dp)*2.0_dp
             shift=shift*10.0_dp**(i-4)
             call matrix_scale(outside_A,1.0_dp,shift)
          end if
       end do

       if(.not.ignore_S) then
          call matrix_free(CC_outside_A)
       end if


       ! Back to Fermi-Dirac from tanh...
       call matrix_scale(CC_outside_T,0.5_dp,0.5_dp)

!       if(pub_on_root) then
!          outunit=utils_unit()
!          open(unit=outunit,file="CC_outside_T2.mtx")
!       end if
!       call sparse_show_matrix(CC_outside_T%dataSPAM3,outunit,matlab_format=.true.)
!       if(pub_on_root) then
!          close(outunit)
!       end if


       call matrix_allocate(CC_T,rho)

       call matrix_allocate(CC_outside_T_inv,CC_outside_T)
       call matrix_copy(CC_outside_T,CC_outside_T_inv)
       call matrix_scale(CC_outside_T_inv,-1.0_dp,1.0_dp)

       call matrix_multiply(1.0_dp,CC_righthand_T,CC_outside_T_inv,1.0_dp,CC_T)
       call matrix_free(CC_outside_T_inv)
       call matrix_axpy(CC_T,CC_outside_T,1.0_dp)

       if(.not.ignore_S) then
          call matrix_free(CC_righthand_T)
       end if

!       if(pub_on_root) then
!          outunit=utils_unit()
!          open(unit=outunit,file="CC_T.mtx")
!       end if
!       call sparse_show_matrix(CC_T%dataSPAM3,outunit,matlab_format=.true.)
!       if(pub_on_root) then
!          close(outunit)
!       end if


    else
       if(.not.ignore_S) then
          call matrix_free(CC_outside_A)
          call matrix_free(CC_righthand_T)

          call matrix_allocate(CC_T,rho)
          call matrix_copy(CC_outside_T,CC_T)
          call matrix_axpy(CC_T,CC_inside_T,1.0_dp)
       end if

    end if

    call matrix_free(CC_inside_T)
    call matrix_free(CC_outside_T)


    if(.not.ignore_S) then
       if(contravariant_output) then

          if(no_solve) then
             call matrix_multiply(1.0_dp,CC_T,invS,0.0_dp,rho)
          else
             call matrix_allocate(CC_T_trans,rho)

             ! transpose CC_outside_T
             call matrix_transpose(CC_T,CC_T_trans)

             call matrix_free(CC_T)

             ! transpose S
             call matrix_allocate(S_trans,S)
             call matrix_transpose(S,S_trans)

             ! solve S'outside_T' = CC_outside_T'
             call matrix_allocate(T_trans,rho)

             call matrix_solve(CC_T_trans,S_trans,T_trans)

             call matrix_free(S_trans)
             call matrix_free(CC_T_trans)

             ! transpose outside_T
             call matrix_transpose(T_trans,rho)

             call matrix_free(T_trans)
          end if

       else ! contra-covariant kernel

          call matrix_copy(CC_T,rho)
          call matrix_free(CC_T)

       end if


    end if

!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="outside_T.mtx")
!    end if
!    call sparse_show_matrix(outside_T%dataSPAM3,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
!    end if


!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="rho.mtx")
!    end if
!    call sparse_show_matrix(rho%dataSPAM3,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
!    end if


    call matrix_free(outside_A)
    call matrix_free(inside_A)
    call matrix_free(inside_T)
    call matrix_free(outside_T)


    call timer_clock('fermi_operator_projection_contourint',2)

  end subroutine fermi_operator_projection_contourint



  !============================================================================!
  ! This routine finds a density matrix for a contra-covariant Hamiltonian     !
  ! matrix CCH. It either outputs a contravariant or contra-covariant density- !
  ! kernel. It uses the so called AQuA-FOE method which we have submitted to   !
  ! JCP for review.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The input covariant Hamiltonian matrix.          !
  !   S            (input)  : The overlap matrix / metric.                     !
  !   CCH          (input)  : The contra-covariant Hamiltonian matrix          !
  !   rho          (output) : The density kernel (sparsity pattern input)      !
  !   orthogonal   (input)  : whether the Hamiltonian is orthogonal            !
  ! contravariant_kernel (input) : make the density kernel contravariant?      !
  ! matmuls_out    (output) : Number of matmuls we needed                      !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, June 2017,                                       !
  !============================================================================!
 subroutine fermi_operator_annealing(H,S,CCH,rho,orthogonal,contravariant_kernel, matmuls_out, invS)


    use rundat, only : pub_FOE_cheby_thres


    use utils, only : utils_dealloc_check
    implicit none
    type(smearing_matrix), intent(inout) :: H
    type(smearing_matrix), intent(inout) :: S
    type(smearing_matrix), intent(inout) :: CCH
    type(smearing_matrix), intent(inout) :: rho ! Inout to make operations at the end of the routine easier.
    logical, optional, intent(in)  :: orthogonal ! Ignore S
    logical, optional, intent(in)  :: contravariant_kernel
    integer, optional, intent(out) :: matmuls_out
    type(smearing_matrix), optional, intent(in) :: invS

    logical :: contravariant_output
    integer :: ierr
    integer :: l

    type(smearing_matrix) :: CCH_loc, CCH_bkp, CCH2, H2
    type(smearing_matrix) :: CCT, CCT_bkp, CCT2
    type(smearing_matrix) :: CCT_trans, S_trans, T_trans

    real(kind=dp), dimension(:), allocatable :: coeffs
    real(kind=dp) :: cheby_threshold
    real(kind=dp) :: cutoff
    real(kind=dp) :: specwidth
    integer :: nlevels
    integer :: matmuls_loc, matmuls_tmp
    integer :: quench_method
    integer,parameter :: quench_solve=111, quench_chebyshev=222, quench_invert=333
    real(kind=dp) :: ttrace
    real(kind=dp), allocatable :: eigvals(:)
    real(kind=dp), allocatable :: eigvecs(:,:)
    integer :: nume


    matmuls_loc=0

    ! ja531-> Instead of doing the quenching with an inversion we can expand the quench formula as a chebyshev
    ! expansion and use only matrix products, which is hard-coded here.
    quench_method=quench_chebyshev
!    quench_method=quench_solve
    ! ja531-> Set max error of Chebyshev expansion of tanh

    cheby_threshold=1e-9_dp
    cheby_threshold=pub_FOE_cheby_thres

    ! ja531-> cutoff is the maximum argument of tanh that is accurate
    cutoff=15.0_dp

    contravariant_output=.true.
    if(present(contravariant_kernel)) then
       contravariant_output=contravariant_kernel
    end if

    if(orthogonal) then
       call matrix_allocate(CCH_loc,H)
       call matrix_copy(H,CCH_loc)
    else
       call matrix_allocate(CCH_loc,CCH)
       call matrix_copy(CCH,CCH_loc)
    end if

    ! convert from argument of FD to argument of tanh: H -> -H/2
    call matrix_scale(CCH_loc,-0.5_dp)

    ! work out spectral width of Hamiltonian = specwidth
    call matrix_allocate(CCH_bkp,CCH_loc)
    call matrix_copy(CCH_loc,CCH_bkp)
    call matrix_allocate(CCH2,CCH_loc,CCH_bkp)
    call matrix_multiply(1.0_dp,CCH_loc,CCH_bkp,0.0_dp,CCH2)
    call matrix_free(CCH_bkp)

!    call matrix_allocate(H2,S,CCH2)
!    call matrix_multiply(1.0_dp,S,CCH2,0.0_dp,H2)


    !We'll use the Frob. norm to give an upper limit to the spectral width
    ! rather than using extremal eigs.
    specwidth=sqrt(matrix_trace(CCH2))

    call matrix_free(CCH2)


!    call extremum_eigs(H2,S,tol=0.00001_dp,maxeval=specwidth)
!    call matrix_free(H2)
!    specwidth=sqrt(specwidth)
!    if(pub_debug_on_root) then
!       write(stdout,*) "Spectral width:",specwidth
!    end if



    ! work out number of annealing levels = log2(specwidth/cutoff)
    ! max to avoid dividing by 0!
!    nlevels = max(ceiling(log(specwidth/cutoff)/log(2.0_dp)),1)+1
    nlevels = max(ceiling(log(specwidth/cutoff)/log(2.0_dp)),0)
!    if(pub_debug_on_root) then
!       write(stdout,*) "Annealing Nlevels=",nlevels
!    end if

    ! H' -> H/2^nlevels
    call matrix_scale(CCH_loc,1.0_dp/real(2**nlevels,dp))

!    ttrace=matrix_trace(CCH_loc)
!    if(pub_debug_on_root) then
!       write(stdout,*) "Annealing CCH_loc: Ktrace=",ttrace
!    end if

    ! calculate tanh(H') -> high electron temperature
    call matrix_scale(CCH_loc,1.0_dp/cutoff)

    call get_tanh_chebyshev_coeffs(-cutoff,cutoff,cheby_threshold,coeffs)

    call matrix_allocate(CCT,rho)
    call fast_cheby_resum(CCH_loc,coeffs,CCT,matmuls_loc)
!    call slow_cheby_resum(CCH_loc,coeffs,CCT)

    call matrix_free(CCH_loc)
    deallocate(coeffs,stat=ierr)
    call utils_dealloc_check('fermi_operator_projection_contourint','cheby_coeffs',ierr)

    matmuls_loc=matmuls_loc+2 ! For the two in the spectral width calculation

    ! Now quench it nlevels times until we have tanh(H)

    call matrix_allocate(CCT2,CCT)
    call matrix_allocate(CCT_bkp,CCT)


    if(quench_method==quench_chebyshev) then
       call get_hypmultid_chebyshev_coeffs(-1.0_dp,1.0_dp,cheby_threshold,coeffs)
    end if



    do l=1,nlevels
       ! tanh(2*x0)=2*tanh(x0)/(tanh(x0)^2 + 1)

       select case(quench_method)
       case(quench_solve)
          call matrix_scale(CCT2,0.0_dp,1.0_dp)                ! T2 = I
          call matrix_copy(CCT,CCT_bkp)
          call matrix_multiply(1.0_dp,CCT,CCT_bkp,1.0_dp,CCT2) ! T2 = T^2 + I


!          call matrix_solve(CCT_bkp,CCT2,CCT)                  ! T -> T/(T^2 + I)

          call matrix_invert(CCT2)
          call matrix_multiply(1.0_dp,CCT_bkp,CCT2,0.0_dp,CCT)

          call matrix_scale(CCT,2.0_dp)                        ! T -> 2T/(T^2 + I)
       case(quench_chebyshev)
          call fast_cheby_resum(CCT,coeffs,CCT_bkp,matmuls_tmp)
          call matrix_copy(CCT_bkp,CCT)
!          ttrace=matrix_trace(CCT)
!          if(pub_debug_on_root) then
!             write(stdout,*) "Quenching ",l," CCT: Ktrace=",ttrace
!          end if
          matmuls_loc=matmuls_loc+matmuls_tmp
       end select

       ! one of these "matmuls" could be a solve...
       matmuls_loc=matmuls_loc+2
    end do

    if(quench_method==quench_chebyshev) then
       deallocate(coeffs,stat=ierr)
       call utils_dealloc_check('fermi_operator_projection_contourint','cheby_coeffs',ierr)
    end if

    call matrix_free(CCT_bkp)
    call matrix_free(CCT2)

    ! Back to Fermi-Dirac from tanh...
    call matrix_scale(CCT,0.5_dp,0.5_dp)


    if(contravariant_output) then

       if(.not.present(invS)) then

          call matrix_allocate(CCT_trans,rho)

          ! transpose CC_outside_T
          call matrix_transpose(CCT,CCT_trans)
          call matrix_free(CCT)

          ! transpose S
          call matrix_allocate(S_trans,S)
          call matrix_transpose(S,S_trans)

          ! solve S'outside_T' = CC_outside_T'
          call matrix_allocate(T_trans,rho)

          call matrix_solve(CCT_trans,S_trans,T_trans)
          matmuls_loc=matmuls_loc+1 ! count it as a matmul for now...

          call matrix_free(S_trans)
          call matrix_free(CCT_trans)

          ! transpose outside_T
          call matrix_transpose(T_trans,rho)

          call matrix_free(T_trans)

       else

          ! rho=CCT*invS
          call matrix_multiply(1.0_dp,CCT,invS,0.0_dp,rho)
          matmuls_loc=matmuls_loc+1
          call matrix_free(CCT)

       end if

    else ! contra-covariant kernel

       call matrix_copy(CCT,rho)
       call matrix_free(CCT)

    end if

    if(present(matmuls_out)) then
       matmuls_out=matmuls_loc
    end if


  end subroutine fermi_operator_annealing

  !============================================================================!
  ! This routine applies the fermi operator to a modified Hamiltonian, such    !
  ! that H -> B*(L'*H*L - I*u), where B is thermodynamic beta... of the        !
  ! electrons, u is the Fermi energy and L is the inverse of a symmetric       !
  ! matrix decomposition of the overlap matrix (such as Cholesky or Lowdin).   !
  ! The result is a density matrix in the orthogonal space. To achieve the     !
  ! density matrix in the original space, the result must be transformed back  !
  ! using L.                                                                   !
  ! This routine constructs two projection operators using a shifted matrix    !
  ! sign function and uses them to construct three eigenvalue subspaces; the   !
  ! first where the eigenvalues are under some -threshold and give 0 to        !
  ! numerical precision, the second, where the eigenvalues are above a         !
  ! +threshold and give 1 to numerical precision (via logistic function) and   !
  ! those that are between -threshold and +threshold, which are dealt with     !
  ! separately, via a series expansion.
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The input scaled Hamiltonian matrix.             !
  !   rho          (output) : The corresponding density matrix.                !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, May 2014,                                        !
  !----------------------------------------------------------------------------!
  ! This routine scales as O(log(N)) in matrix-matrix multiplications.         !
  !============================================================================!
  subroutine fermi_operator_projection(H,rho,matmuls_out)
    use comms, only : pub_on_root
    use rundat, only : pub_edft_smearing_width
    use timer, only : timer_clock
    use utils, only : utils_alloc_check, utils_dealloc_check
    implicit none
    type(smearing_matrix), intent(inout) :: H
    type(smearing_matrix), intent(inout) :: rho ! Inout to make operations at the end of the routine easier.
    integer, optional, intent(out) :: matmuls_out

    type(smearing_matrix) :: A
    type(smearing_matrix) :: low_T, medium_T, high_T, sign_cutoff
    type(smearing_matrix) :: low_A, medium_A
    real(kind=dp) :: tmp
    type(smearing_matrix) :: tmp1, tmp2
    real(kind=dp) :: eigtol, signtol
    real(kind=dp) :: min_abs_eig, min_eig, max_eig
    real(kind=dp), dimension(:), allocatable :: cheby_coeffs, cheby_coeffs_full, coeffs
    integer :: mmsi, matmuls_cheby, matmuls

    integer :: status
    real(kind=dp) :: shift
    integer :: i
    logical :: no_extremal_eigs
    logical :: use_gerschgorin
    integer :: ierr
    real(kind=dp) :: res
    real(kind=dp) :: cheby_threshold

    real(kind=dp) :: abs_emax

    real(kind=dp), save :: cutoff=15.000_dp

    integer :: outunit

    call timer_clock('fermi_operator_projection',1)

    ! ja531-> In the mapping from Fermi-Dirac to tanh, the argument A-> -H/2
    call matrix_allocate(A,H)
    call matrix_copy(H,A)
    call matrix_scale(A,-0.5_dp)


!    call matrix_allocate(medium_A,A)
!
!
!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="A.mtx")
!    end if
!    call sparse_show_matrix(A%dataSPAM3,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
!    end if
!
!
!    call boxcar_projection_contourint(A,32,0.0_dp,15.0_dp,0.5_dp,medium_A)
!
!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="medium_A.mtx")
!    end if
!    call sparse_show_matrix(medium_A%dataSPAM3,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
!    end if
!
!    call matrix_free(medium_A)



    matmuls=0

    eigtol=1e-6_dp
    signtol=1e-7_dp
    res=1e-7_dp
    cheby_threshold=res
    !commented the next line... because of convergence issues with sparse matrices with too small a cutoff....
!    cutoff=min(log(((1.0_dp-res)+1.0_dp)/(-(1.0_dp-res)+1.0_dp))/2.0_dp,15.0_dp) !(that's atanh!)

    no_extremal_eigs=.true.
    use_gerschgorin=.false.
    ! Do the method without extremal eigs?
    if(no_extremal_eigs) then
       if(use_gerschgorin) then
          write(*,*) "max_eig:::"
          max_eig=matrix_induced_norm(A)
          write(*,*) "max_eig was:",max_eig,pub_edft_smearing_width

!          call matrix_allocate(tmp1,A,A)
!          call matrix_multiply(1.0_dp,A,A,0.0_dp,tmp1)
!          max_eig=sqrt(abs(matrix_trace(tmp1)))
!          write(*,*) "max_eig was:",max_eig
!          call matrix_free(tmp1)
!          max_eig=5.0_dp


!          cutoff=log(((max_eig*1.1_dp)+1.0_dp)/(-(max_eig*1.1_dp)+1.0_dp))/2.0_dp !(that's atanh!)
          cutoff=1.1_dp*max_eig
          write(*,*) "cutoff was:",cutoff,pub_edft_smearing_width

          min_eig=-max_eig
!          if(1.1_dp*max_eig>15.0_dp) then
!             cutoff=15.0_dp
!          else
!             cutoff=1.1_dp*max_eig
!          end if
!
         min_abs_eig=0.1_dp*cutoff

       else
       min_eig=-1.1_dp*cutoff
       max_eig=1.1_dp*cutoff
       min_abs_eig=0.9_dp*cutoff
       end if
    else
       ! ja531 --> 18/09/14
       call matrix_allocate(tmp1,H)
       call matrix_allocate(tmp2,H)
       ! ja531-> tmp2 is just an identity matrix...
       call matrix_scale(tmp2,0.0_dp)
       call matrix_scale(tmp2,0.0_dp,1.0_dp)
       call matrix_multiply(-0.5_dp,H,A,0.0_dp,tmp1)
       matmuls=matmuls+1
       ! ja531 --> 18/09/14
       call extremum_eigs(tmp1,tmp2,eigtol,mineval=min_abs_eig)
       min_abs_eig=sqrt(min_abs_eig)

       call extremum_eigs(A,tmp2,eigtol,min_eig,max_eig)
       call matrix_free(tmp1)
       call matrix_free(tmp2)
    end if


    ! ja531 --> 18/09/14
    !min_abs_eig=abs(max_eig)/condition_number(A)

    if(min_abs_eig>cutoff) then
       !ja531->sparsefoe:
!       call matrix_allocate(medium_T,H)
       call matrix_allocate(medium_T,rho)

       shift=0.0_dp
       i=0
       do
          i=i+1
          abs_emax=-1.0_dp
          call sign_nsh(A,signtol,medium_T,abs_emax,mmsi,status)
          if(pub_on_root) then
             if(debug_info_toggle) then
                write(stdout,*) "smearing-> 1st Sign: ",mmsi
             end if
          end if
          matmuls=matmuls+mmsi
          if(status==0) then
             exit
          else
             call random_number(shift)
             shift=(shift-1.0_dp)*2.0_dp
             shift=shift*10.0_dp**(i-8)
             call matrix_scale(A,1.0_dp,shift)
          end if
       end do

    else
       if(max_eig>cutoff) then
          call matrix_allocate(tmp1,H)
          call matrix_copy(A,tmp1)
          call matrix_scale(tmp1,1.0_dp,-cutoff)
          !ja531->sparsefoe:
!          call matrix_allocate(sign_cutoff,H)
          call matrix_allocate(sign_cutoff,rho)


          shift=0.0_dp
          i=0
          do
             i=i+1
             if(pub_on_root) then
                if(debug_info_toggle) then
                   write(stdout,*) "smearing-> 2nd Sign: cutoff=",cutoff
                end if
             end if

             abs_emax=-1.0_dp
             call sign_nsh(tmp1,signtol,sign_cutoff,abs_emax,mmsi,status)
             if(pub_on_root) then
                if(debug_info_toggle) then
                   write(stdout,*) "smearing-> 2nd Sign: ",mmsi
                end if
             end if
             matmuls=matmuls+mmsi
             if(status==0) then
                exit
             else
                call random_number(shift)
                shift=(shift)*2.0_dp
                shift=shift*10.0_dp**(i-4) !14 is the smallest value that actually results in any change at all...
                cutoff=cutoff+shift
                call matrix_scale(tmp1,1.0_dp,-shift)
             end if
          end do

          !ja531->sparsefoe:
!          call matrix_allocate(high_T,H)
          call matrix_allocate(high_T,rho)
          call matrix_copy(sign_cutoff,high_T)
          call matrix_scale(high_T,0.5_dp,0.5_dp)
          !ja531->sparsefoe:
          call matrix_free(tmp1)
          call matrix_allocate(tmp1,rho)

          call matrix_multiply(1.0_dp,A,sign_cutoff,0.0_dp,tmp1)
          call matrix_free(sign_cutoff)
          matmuls=matmuls+1

          !ja531->sparsefoe:
!          call matrix_allocate(low_A,H)
          call matrix_allocate(low_A,rho)

          call matrix_copy(A,low_A)
          call matrix_axpy(low_A,tmp1,-1.0_dp)
          call matrix_scale(low_A,0.5_dp)
          call matrix_free(tmp1)
       else
          !ja531->sparsefoe:
!          call matrix_allocate(low_A,H)
          call matrix_allocate(low_A,rho)
          call matrix_copy(A,low_A)
       end if

       if(min_eig<(-cutoff)) then
          !ja531->sparsefoe:
!          call matrix_allocate(tmp1,H)
          call matrix_allocate(tmp1,rho)

          call matrix_copy(low_A,tmp1)
          call matrix_scale(tmp1,1.0_dp,cutoff)
          !ja531->sparsefoe:
!          call matrix_allocate(sign_cutoff,H)
          call matrix_allocate(sign_cutoff,rho)

          shift=0.0_dp
          i=0
          do
             i=i+1
             if(pub_on_root) then
                if(debug_info_toggle) then
                   write(stdout,*) "smearing-> 3rd Sign: cutoff=",cutoff
                end if
             end if

             abs_emax=-1.0_dp
             call sign_nsh(tmp1,signtol,sign_cutoff,abs_emax,mmsi,status)
             if(pub_on_root) then
                if(debug_info_toggle) then
                   write(stdout,*) "smearing-> 3rd Sign: ",mmsi
                end if
             end if
             matmuls=matmuls+mmsi
             if(status==0) then
                exit
             else
                call random_number(shift)
                shift=(shift)*2.0_dp
                shift=shift*10.0_dp**(i-4)
                cutoff=cutoff+shift
                call matrix_scale(tmp1,1.0_dp,shift)
             end if
          end do


          call matrix_multiply(1.0_dp,low_A,sign_cutoff,0.0_dp,tmp1)
          matmuls=matmuls+1
          !ja531->sparsefoe:
!          call matrix_allocate(medium_A,H)
          call matrix_allocate(medium_A,rho)

          call matrix_copy(low_A,medium_A)
          call matrix_axpy(medium_A,tmp1,1.0_dp)
          call matrix_scale(medium_A,0.5_dp)
          call matrix_free(tmp1)
          !ja531->sparsefoe:
!          call matrix_allocate(low_T,H)
          call matrix_allocate(low_T,rho)

          call matrix_copy(sign_cutoff,low_T)
          call matrix_scale(low_T,0.5_dp,-0.5_dp)
          call matrix_free(sign_cutoff)
       else
          !ja531->sparsefoe:
!          call matrix_allocate(medium_A,H)
          call matrix_allocate(medium_A,rho)

          call matrix_copy(A,medium_A)
       end if
    end if

    ! Do the method without extremal eigs?
    if(no_extremal_eigs) then
       !ja531->sparsefoe:
!       call matrix_allocate(tmp1,medium_A,medium_A)
       call matrix_allocate(tmp1,rho)

       call matrix_multiply(1.0_dp,medium_A,medium_A,0.0_dp,tmp1)
       if(sqrt(abs(matrix_trace(tmp1)))<eigtol) then
          min_abs_eig=1.1_dp*cutoff
       else
          min_abs_eig=0.9_dp*cutoff
       end if
       call matrix_free(tmp1)
    end if

    if(min_abs_eig<cutoff) then
       !ja531->sparsefoe:
!       call matrix_allocate(medium_T,H)
       call matrix_allocate(medium_T,rho)

       call matrix_scale(medium_A,1.0_dp/cutoff)
       call get_tanh_chebyshev_coeffs(-cutoff,cutoff,cheby_threshold,coeffs)
       ! Could calculate Chebyshev coeffs for different cutoffs... at the moment it's just hard-coded for a cutoff of 15.
       allocate(cheby_coeffs_full(337),stat=ierr)
       call utils_alloc_check('fermi_operator_projection','cheby_coeffs_full',ierr)
       cheby_coeffs_full=(/-4.9563527885051629e-19_dp, 1.2709033457233794e+00_dp, 8.9937151641416602e-18_dp,&
            & -4.1747689700886825e-01_dp,&
            & 1.6108146562641781e-18_dp, 2.4332169991823077e-01_dp, -9.4170702981598104e-18_dp, -1.6650421357427858e-01_dp,&
            & 1.1151793774136617e-17_dp, 1.2244297878905804e-01_dp, -1.5695117163599683e-17_dp, -9.3559089908144502e-02_dp,&
            & -2.8323041779230504e-18_dp, 7.3096077046332114e-02_dp, -4.5440540164708237e-18_dp, -5.7887575992751496e-02_dp,&
            & -2.1533828539238099e-17_dp, 4.6234862082162362e-02_dp, 7.6968141943436615e-18_dp, -3.7129231260117988e-02_dp,&
            & -1.2490039732122380e-17_dp, 2.9922093849227099e-02_dp, -5.5308203234549153e-18_dp, -2.4169440040929274e-02_dp,&
            & 7.6298479910485601e-18_dp, 1.9552246139928151e-02_dp, 4.9241770249532959e-18_dp, -1.5832829825013831e-02_dp,&
            & 6.8426626842610182e-18_dp, 1.2829378889803532e-02_dp, -1.8112274081546521e-18_dp, -1.0400193108220981e-02_dp,&
            & -9.1900707643032933e-18_dp, 8.4333893838576448e-03_dp, -1.0374131478332715e-17_dp, -6.8398369646402845e-03_dp,&
            & -2.3175098592939044e-19_dp, 5.5481002094077197e-03_dp, 1.2142141463121059e-17_dp, -4.5006921033925347e-03_dp,&
            & 1.0265106964206103e-17_dp, 3.6512242609227159e-03_dp, -7.3894890018990132e-19_dp, -2.9621959000831455e-03_dp,&
            & -2.8326586485436946e-18_dp, 2.4032541929801306e-03_dp, -1.4693550585431695e-17_dp, -1.9498119212930805e-03_dp,&
            & -3.0064339889014400e-18_dp, 1.5819415586354833e-03_dp, -7.4770170201150985e-18_dp, -1.2834864000734550e-03_dp,&
            & -5.8906809586705316e-18_dp, 1.0413439063346928e-03_dp, -2.1794253960530038e-18_dp, -8.4488668391553400e-04_dp,&
            & 1.2188528430963106e-17_dp, 6.8549401096665199e-04_dp, 7.1372575505099958e-19_dp, -5.5617244870159862e-04_dp,&
            & -2.7003332757397526e-18_dp, 4.5124840605533049e-04_dp, -9.5466193298613829e-18_dp, -3.6611890782022445e-04_dp,&
            & -9.6628826705081059e-18_dp, 2.9704949143643605e-04_dp, 9.6937847515425776e-19_dp, -2.4101029054447373e-04_dp,&
            & -1.2250343148327254e-18_dp, 1.9554307382404305e-04_dp, -8.1637305697375222e-20_dp, -1.5865338454094423e-04_dp,&
            & -2.5574551858756606e-18_dp, 1.2872303752695780e-04_dp, 3.2893500384959963e-18_dp, -1.0443912882177937e-04_dp,&
            & 8.8011737346884447e-18_dp, 8.4736440512242213e-05_dp, -2.2922088616697569e-18_dp, -6.8750712487534213e-05_dp,&
            & 1.1088788921554525e-18_dp, 5.5780730369880899e-05_dp, 2.4679855119358769e-18_dp, -4.5257566067381714e-05_dp,&
            & -3.4702813760646623e-19_dp, 3.6719621395211069e-05_dp, -4.5519423139845883e-18_dp, -2.9792379897539745e-05_dp,&
            & 1.0730904586118607e-18_dp, 2.4171978601862854e-05_dp, -2.2037323627361310e-18_dp, -1.9611879034348674e-05_dp,&
            & 9.9963186958215794e-18_dp, 1.5912052799005866e-05_dp, 7.1845748667449528e-18_dp, -1.2910207331881769e-05_dp,&
            & -1.5311446429506578e-18_dp, 1.0474666943055726e-05_dp, 6.5246561608025744e-18_dp, -8.4985968691403879e-06_dp,&
            & -2.0041089592077674e-17_dp, 6.8953169733607542e-06_dp, -8.0958225257451980e-19_dp, -5.5944995277822019e-06_dp,&
            & -3.9509794498474893e-18_dp, 4.5390842931327300e-06_dp, -1.0816963975765377e-17_dp, -3.6827755761471381e-06_dp,&
            & -6.1915441555418331e-18_dp, 2.9880114730440191e-06_dp, 7.9561981956189958e-18_dp, -2.4243162198494409e-06_dp,&
            & 7.6670914399214032e-19_dp, 1.9669633758083432e-06_dp, 3.3096272010303460e-18_dp, -1.5958912002326034e-06_dp,&
            & -2.2038415912061615e-18_dp, 1.2948226460611237e-06_dp, 2.2660388331764904e-17_dp, -1.0505513687538957e-06_dp,&
            & 2.1914978377380792e-18_dp, 8.5236243106621825e-07_dp, -4.9494224726115693e-19_dp, -6.9156229337547085e-07_dp,&
            & 3.7610619479572694e-20_dp, 5.6109747236687427e-07_dp, 2.3917541390661184e-18_dp, -4.5524514066136444e-07_dp,&
            & -6.9837497234454390e-18_dp, 3.6936209531155301e-07_dp, 2.0450723930024455e-18_dp, -2.9968108450221034e-07_dp,&
            & -1.0057806152030866e-18_dp, 2.4314555701544287e-07_dp, -4.0447303954275294e-18_dp, -1.9727558712741883e-07_dp,&
            & 4.0394847245935221e-18_dp, 1.6005909281500775e-07_dp, 1.5629194061188890e-18_dp, -1.2986357595446523e-07_dp,&
            & -1.1920526555932951e-17_dp, 1.0536451296337744e-07_dp, -1.9468965721730738e-17_dp, -8.5487254667945299e-08_dp,&
            & 1.3449585684849298e-17_dp, 6.9359887007541155e-08_dp, -2.1195460920643224e-17_dp, -5.6274984408289297e-08_dp,&
            & 1.7258014423064345e-17_dp, 4.5658578865406015e-08_dp, 2.2193166858357647e-18_dp, -3.7044982701198159e-08_dp,&
            & 2.7876381708010217e-18_dp, 3.0056361303590750e-08_dp, -6.6864029180296343e-19_dp, -2.4386159475186283e-08_dp,&
            & 1.9417271395704195e-18_dp, 1.9785654263780395e-08_dp, 5.7824115865893567e-19_dp, -1.6053044973386456e-08_dp,&
            & 1.6603781841492296e-17_dp, 1.3024601046370782e-08_dp, -4.8737469086967434e-18_dp, -1.0567480010550504e-08_dp,&
            & 2.4162219843962669e-18_dp, 8.5739005455572527e-09_dp, -8.0024446064406280e-18_dp, -6.9564144282903877e-09_dp,&
            & 6.6084703846735505e-19_dp, 5.6440708399930181e-09_dp, -7.2176887482606436e-18_dp, -4.5793038498070131e-09_dp,&
            & 6.4019556851525027e-19_dp, 3.7154076407323811e-09_dp, 2.4781763942525814e-19_dp, -3.0144873983332600e-09_dp,&
            & -6.0302292260146153e-18_dp, 2.4457974118621298e-09_dp, 6.3606527452482924e-18_dp, -1.9843920983932842e-09_dp,&
            & 5.6264804260744297e-19_dp, 1.6100319541533819e-09_dp, -2.5082324116898487e-18_dp, -1.3062957134641578e-09_dp,&
            & 1.8824830074395601e-17_dp, 1.0598600126665846e-09_dp, -6.2567039748132132e-18_dp, -8.5991499693061608e-10_dp,&
            & -5.4386347588675154e-18_dp, 6.9769004637636701e-10_dp, 4.2578230893332992e-18_dp, -5.6606920118591162e-10_dp,&
            & 8.6379849752428556e-18_dp, 4.5927893622822404e-10_dp, -1.7632256044822513e-18_dp, -3.7263491730108035e-10_dp,&
            & -2.3491362856089391e-18_dp, 3.0233646885638853e-10_dp, 1.4734107304753761e-18_dp, -2.4530000980196202e-10_dp,&
            & 1.9287721959967893e-18_dp, 1.9902360265785672e-10_dp, 4.8092341001429454e-18_dp, -1.6147733047360334e-10_dp,&
            & -1.1737570455030390e-18_dp, 1.3101426972489262e-10_dp, 6.5732777246935424e-18_dp, -1.0629811705126488e-10_dp,&
            & -5.0991660705137990e-18_dp, 8.6244700653384062e-11_dp, 4.1198239429776594e-19_dp, -6.9974511350818263e-11_dp,&
            & 6.8524945638801110e-18_dp, 5.6773649934721822e-11_dp, 6.5375027384734948e-18_dp, -4.6063142532281901e-11_dp,&
            & -1.9143271026867869e-18_dp, 3.7373259125347438e-11_dp, -5.4743363756023469e-18_dp, -3.0322698551812250e-11_dp,&
            & -4.2981232251026185e-18_dp, 2.4602248355276050e-11_dp, -1.6432062045786659e-18_dp, -1.9960966323746963e-11_dp,&
            & -8.0846944316387981e-19_dp, 1.6195291124421792e-11_dp, 6.1725351460013987e-18_dp, -1.3140021687400827e-11_dp,&
            & 7.9839973687529438e-18_dp, 1.0661130255361907e-11_dp, 1.8427084162685988e-17_dp, -8.6498770185878782e-12_dp,&
            & 2.5869686866723071e-17_dp, 7.0180509043370498e-12_dp, 2.1140508992040653e-18_dp, -5.6940900588152695e-12_dp,&
            & -1.7291067932484984e-17_dp, 4.6198864787888725e-12_dp, -4.6843958152511333e-18_dp, -3.7483475163445817e-12_dp,&
            & 8.2909257095095754e-19_dp, 3.0412072083535731e-12_dp, -7.1845748667449528e-18_dp, -2.4674780293155810e-12_dp,&
            & 3.2206220735255234e-18_dp, 2.0019860387108209e-12_dp, 7.1682652618458181e-19_dp, -1.6243067107533596e-12_dp,&
            & -8.2527281918660270e-19_dp, 1.3178899277370921e-12_dp, 1.0313702430621841e-17_dp, -1.0692422044131845e-12_dp,&
            & -3.6180540931976645e-18_dp, 8.6756278714733004e-13_dp, 1.3536903700945378e-17_dp, -7.0385623250549572e-13_dp,&
            & 2.7735974588402591e-18_dp, 5.7108908479358037e-13_dp, 1.7965735828192408e-18_dp, -4.6334036041625765e-13_dp,&
            & 1.3006778534734273e-17_dp, 3.7595116648606174e-13_dp, -1.9676559615612862e-18_dp, -3.0502658930421841e-13_dp,&
            & -2.1475303344843273e-18_dp, 2.4747183056091018e-13_dp, -1.0476773682678570e-17_dp, -2.0079824872948613e-13_dp,&
            & 3.2310555621261372e-18_dp, 1.6291876495364673e-13_dp, 1.2239456867930135e-18_dp, -1.3219749369260589e-13_dp,&
            & -5.7358919068462931e-19_dp, 1.0724799851112907e-13_dp, -2.2133423852104709e-17_dp, -8.7029424672160242e-14_dp,&
            & 3.6039309736331378e-18_dp, 7.0609213747072202e-14_dp, -6.6342805035750825e-19_dp, -5.7287136344198945e-14_dp,&
            & 6.2378677145683460e-18_dp, 4.6480697102071371e-14_dp, -3.1050496935963413e-18_dp, -3.7714044850053106e-14_dp,&
            & 1.1070951698044603e-17_dp, 3.0594574492884669e-14_dp, -2.4397175129065546e-18_dp, -2.4838761999593626e-14_dp,&
            & -2.5575130538886930e-18_dp, 2.0132374603388738e-14_dp, -4.7311624347445844e-18_dp, -1.6343986349494816e-14_dp,&
            & 1.6965749876426747e-17_dp, 1.3261362081214078e-14_dp, -2.0647876331739830e-18_dp, -1.0740602355920258e-14_dp,&
            & -6.1906462459523470e-19_dp, 8.7272285958796989e-15_dp, 4.1864841213256094e-18_dp, -7.0766805114276721e-15_dp,&
            & 1.0895945697008061e-18_dp, 5.7558124932910459e-15_dp, 3.6620044274111045e-18_dp, -4.6616150093487231e-15_dp,&
            & -8.6734328486994116e-18_dp, 3.7800450600332715e-15_dp, -1.5721308660460798e-18_dp, -3.0663302584885277e-15_dp,&
            & -1.6244411996817680e-17_dp, 2.4913933350219286e-15_dp, -3.1577930527930302e-18_dp, -2.0182268554793023e-15_dp,&
            & 1.9505764380572537e-18_dp, 1.6309704909374323e-15_dp, -3.8050039774011130e-17_dp, -1.3164073006269714e-15_dp,&
            & -1.1011575487757301e-17_dp, 1.0705722023171152e-15_dp, 7.7294474766523072e-19_dp, -8.5777945593062691e-16_dp,&
            & -1.7599375828611617e-18_dp, 6.9124600223685343e-16_dp, 1.1152301055049299e-17_dp, -5.4453795969710063e-16_dp,&
            & -8.7350379656867455e-19_dp, 4.2822888092684610e-16_dp, -4.9563527885051629e-19_dp, -3.1720657846433045e-16_dp,&
            & 8.2605879808419381e-20_dp, 2.2204460492503131e-16_dp, 1.4869058365515489e-18_dp, -1.6917684184764289e-16_dp,&
            & -1.5075573065036537e-17_dp, 8.4588420923821446e-17_dp, -2.4058962494202145e-18_dp, -1.6917684184764289e-16_dp,&
            & 1.1564823173178713e-18_dp/)

       call fast_cheby_resum(medium_A,coeffs,medium_T,matmuls_cheby)


       deallocate(coeffs,stat=ierr)
       call utils_dealloc_check('fermi_operator_projection','cheby_coeffs',ierr)

       matmuls=matmuls+matmuls_cheby
    end if

    call matrix_scale(rho,0.0_dp)
    if(low_T%matrix_is_allocated) then
       call matrix_axpy(rho,low_T,1.0_dp)
       call matrix_free(low_T)
       call matrix_free(low_A)
    end if
    if(medium_T%matrix_is_allocated) then
       call matrix_axpy(rho,medium_T,1.0_dp)
       call matrix_free(medium_T)
       call matrix_free(medium_A)
    end if
    if(high_T%matrix_is_allocated) then
       call matrix_axpy(rho,high_T,1.0_dp)
       call matrix_free(high_T)
    end if

    !ja531->sparsefoe: (avoiding this)
! !    call matrix_allocate(tmp1,H)
!     call matrix_allocate(tmp1,rho)
!
!     call matrix_transpose(rho,tmp1)
!     call matrix_axpy(rho,tmp1,1.0_dp)
!     call matrix_free(tmp1)
!    call matrix_scale(rho,0.5_dp)

    ! ja531-> Convert back to F-D function from tanh
    call matrix_scale(rho,1.0_dp,1.0_dp)
    call matrix_scale(rho,0.5_dp)

    call matrix_free(A)

    if(present(matmuls_out)) then
       matmuls_out=matmuls
    end if

    call timer_clock('fermi_operator_projection',2)

  end subroutine fermi_operator_projection



  !============================================================================!
  ! This routine applies a Fermi-Operator expansion using the algorithm of     !
  ! Liang and Head Gordon:                                                     !
  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The input scaled Hamiltonian matrix.             !
  !   rho          (output) : The corresponding density matrix.                !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Mar 2017,                                        !
  !============================================================================!
  subroutine fermi_operator_head_gordon(H,rho,matmuls_out)



    use timer, only : timer_clock
    use utils, only : utils_dealloc_check
    implicit none
    type(smearing_matrix), intent(inout) :: H
    type(smearing_matrix), intent(inout) :: rho ! Inout to make operations at the end of the routine easier.
    integer, optional, intent(out) :: matmuls_out

    type(smearing_matrix) :: A
    type(smearing_matrix) :: T
    real(kind=dp) :: cutoff
    real(kind=dp), dimension(:), allocatable :: coeffs
    integer :: mmsi, matmuls_cheby, matmuls

    integer :: ierr
    integer :: status
    integer :: i
    real(kind=dp) :: cheby_threshold

    integer :: outunit

    call timer_clock('fermi_operator_head_gordon',1)

    ! ja531-> In the mapping from Fermi-Dirac to tanh, the argument A-> -H/2
    call matrix_allocate(A,H)
    call matrix_copy(H,A)
    call matrix_scale(A,-0.5_dp)

    matmuls=0

    cutoff=1.1_dp*matrix_induced_norm(A)

    cheby_threshold=1e-7_dp

    !ja531->sparsefoe:
    !       call matrix_allocate(medium_T,H)
    call matrix_allocate(T,rho)

    call matrix_scale(A,1.0_dp/cutoff)
    call get_tanh_chebyshev_coeffs(-cutoff,cutoff,cheby_threshold,coeffs)

    call fast_cheby_resum(A,coeffs,T,matmuls_cheby)

    deallocate(coeffs,stat=ierr)
    call utils_dealloc_check('fermi_operator_head_gordon','cheby_coeffs',ierr)

    matmuls=matmuls+matmuls_cheby

    call matrix_scale(rho,0.0_dp)

    call matrix_axpy(rho,T,1.0_dp)
    call matrix_free(T)
    call matrix_free(A)

    ! ja531-> Convert back to F-D function from tanh
    call matrix_scale(rho,1.0_dp,1.0_dp)
    call matrix_scale(rho,0.5_dp)

    if(present(matmuls_out)) then
       matmuls_out=matmuls
    end if

    call timer_clock('fermi_operator_head_gordon',2)

  end subroutine fermi_operator_head_gordon



  !============================================================================!
  ! This routine gets the tanh Chebyshev coefficients to approximate the       !
  ! function between two interval limits to a given threshold. The algorithm   !
  ! Uses a discrete cosine transform to give the coefficients (which uses the  !
  ! 1D Fourier transform of the linalg module).
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   lower        (input)  : The lower limit of the interval                  !
  !   upper        (input)  : The upper limit of the interval                  !
  !   threshold    (input)  : The accuracy threshold of the Chebyshev          !
  !                            approximation on the interval                   !
  !   coeffs      (output)  : The Chebyshev coefficients                       !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Apr 2017,                                        !
  !============================================================================!

  subroutine get_tanh_chebyshev_coeffs(lower,upper,threshold,coeffs)
    use comms, only: pub_on_root
    use constants, only: dp
    use utils, only : utils_alloc_check, utils_dealloc_check
    implicit none
    real(kind=dp), intent(in) :: lower
    real(kind=dp), intent(in) :: upper
    real(kind=dp), intent(in) :: threshold
    real(kind=dp), dimension(:), allocatable, intent(inout) :: coeffs

    integer :: n
    integer :: ierr
    real(kind=dp) :: halfwidth, midpoint
    real(kind=dp) :: shift
    integer :: i


    if(allocated(coeffs)) then
       deallocate(coeffs,stat=ierr)
       call utils_dealloc_check('get_tanh_chebyshev_coeffs','coeffs',ierr)
    end if

    halfwidth=(upper-lower)/2.0_dp
    midpoint=(lower+upper)/2.0_dp

    n = howmany(lower,upper,threshold)
    !n=337
    if(pub_on_root) then
       if(debug_info_toggle) then
          write(*,*) "n=",n," lower=",lower," upper=",upper
       end if
    end if

    if(n>0) then
       allocate(coeffs(n),stat=ierr)
       call utils_alloc_check('get_tanh_chebyshev_coeffs','coeffs',ierr)

       call tanh_chebyshev_coeffs(n,lower,upper,coeffs,shift)

    else
       allocate(coeffs(-n),stat=ierr)
       call utils_alloc_check('get_tanh_chebyshev_coeffs','coeffs',ierr)

       call tanh_chebyshev_coeffs(-n,lower,upper,coeffs,shift)

    end if

  contains

    ! A function to give the number of Chebyshev coefficients needed to represent
    ! the tanh function to a given accuracy threshold on the interval [lower, upper].
    function howmany(lower,upper,threshold) result(nterms)
      implicit none
      integer :: nterms
      real(kind=dp), intent(in) :: lower
      real(kind=dp), intent(in) :: upper
      real(kind=dp), intent(in) :: threshold

      real(kind=dp) :: halfwidth, midpoint
      real(kind=dp), dimension(:), allocatable :: coeffs
      integer :: i
      integer :: maxterms
      real(kind=dp) :: error, olderror
      integer :: wentup
      real(kind=dp) :: shift
      integer :: ierr

      halfwidth=(upper-lower)/2.0_dp
      midpoint=(lower+upper)/2.0_dp

      ! This is arbitrary... For a purely Head-Gordon scheme, we wouldn't want this limit.
      ! It is set to this because at 341, we achieve values of tanh(argument) which are
      ! indistinguishable from -1 and 1 to machine precision and can hence project eigs
      ! outside of this range.
!      maxterms=341
      maxterms=10000

      error=huge(1.0_dp)
      olderror=huge(1.0_dp)

      nterms=-1

      do i=3,maxterms,1
         allocate(coeffs(i),stat=ierr)
         call utils_alloc_check('get_tanh_chebyshev_coeffs_howmany','coeffs',ierr)

         call tanh_chebyshev_coeffs(i,lower,upper,coeffs,shift)
         olderror=error
         error=(abs(chebpolval(coeffs, (lower - midpoint)/halfwidth) - tanh(lower) + shift)+ &
              & abs(chebpolval(coeffs, (upper - midpoint)/halfwidth) - tanh(upper) + shift))/2.0_dp
!         write(*,*) i,": delta=",error, (chebpolval(coeffs, (lower - midpoint)/halfwidth) - tanh(lower) + shift), &
!              & (chebpolval(coeffs, (upper - midpoint)/halfwidth) - tanh(upper) + shift)
         if(error>olderror) then
            wentup=wentup+1
         else
            wentup=0
         end if
         if(wentup>=3) then
            ! failed to converge
            nterms=-i
            exit
         end if

         if(error<threshold) then
            nterms=i
            deallocate(coeffs,stat=ierr)
            call utils_dealloc_check('get_tanh_chebyshev_coeffs_howmany','coeffs',ierr)
            exit
         end if
         deallocate(coeffs,stat=ierr)
         call utils_dealloc_check('get_tanh_chebyshev_coeffs_howmany','coeffs',ierr)
      end do

    end function howmany

    ! This function evaluates a Chebyshev polynomial with coeffs c, for a given
    ! argument x (between -1 and 1).
    function chebpolval(c,x) result(u)
      implicit none
      real(kind=dp) :: u
      real(kind=dp), dimension(:), intent(in) :: c
      real(kind=dp),               intent(in) :: x
      integer :: n
      real(kind=dp) :: ujp1, ujp2
      integer :: j

      n = size(c)
      u = c(n)
      if(n > 1) then
         ujp1 = u
         u = c(n-1) + 2.0_dp*x*c(n)
         do j = n-2,1,-1
            ujp2 = ujp1
            ujp1 = u
            u = c(j) + 2.0_dp*x*ujp1 - ujp2
         end do
         u = u - x*ujp1
      end if
    end function chebpolval


    ! This routine gives Chebyshev coefficients which approximate tanh on an
    ! interval [start, finish]. The number of coeffs is requested using n.
    ! the shift  gives the constant shift between the approximation and the
    ! tanh function.
    subroutine tanh_chebyshev_coeffs(n,start,finish,alpha, shift)
      implicit none
      integer, intent(in) :: n
      real(kind=dp), intent(in) :: start
      real(kind=dp), intent(in) :: finish
      real(kind=dp), dimension(n), intent(out) :: alpha ! function coefficients
      real(kind=dp), intent(out) :: shift

      real(kind=dp),  dimension(n) :: y
      !    real(kind=dp), dimension(n) :: alpha ! function coefficients
      real(kind=dp), dimension(n) :: theta, costheta
      integer(kind=dp) :: plandct
      integer :: i
      real(kind=dp) :: halfwidth, midpoint

      halfwidth=(finish-start)/2.0_dp
      midpoint=(start+finish)/2.0_dp

      do i=1,n
         theta(i)=(real(i-1,dp)+0.5_dp)*pi/real(n,dp)
         costheta(i)=cos(theta(i))
         y(i)=tanh(midpoint + costheta(i)*halfwidth)
      end do

      ! DCT
      ! We could use a DCT directly by doing this:
      !    call dfftw_plan_r2r_1d(plandct, n, y, alpha, FFTW_REDFT10, FFTW_ESTIMATE)
      !    call dfftw_execute_r2r(plandct, y, alpha)
      !    call dfftw_destroy_plan(plandct)

      call cos_trans(y,alpha)

      alpha=(alpha/real(n,dp))
      shift=(tanh(midpoint)-chebpolval(alpha,0.0_dp))

    end subroutine tanh_chebyshev_coeffs


    ! An inefficient cosine transform using ONETEPs built-in 1d FFT.
    subroutine cos_trans(invec, outvec)
      use comms, only: pub_on_root, pub_root_proc_id,comms_bcast
      use linalg, only: linalg_1d_fft
      implicit none
      real(kind=dp),    dimension(:), intent(in)    :: invec
      real(kind=dp),    dimension(:), intent(out)   :: outvec
      integer :: n
      complex(kind=dp), dimension(:), allocatable   :: dftvec
      complex(kind=dp), dimension(:), allocatable   :: tmpcmplx
      real(kind=dp),    dimension(:), allocatable   :: tmpvec
      integer :: ierr
      integer :: i

      n=size(invec)
      allocate(dftvec(2*n),stat=ierr)
      call utils_alloc_check('get_tanh_chebyshev_coeffs_costran','dftvec',ierr)
      allocate(tmpvec(2*n),stat=ierr)
      call utils_alloc_check('get_tanh_chebyshev_coeffs_costran','tmpvec',ierr)
      allocate(tmpcmplx(2*n),stat=ierr)
      call utils_alloc_check('get_tanh_chebyshev_coeffs_costran','tmpcmplx',ierr)

      tmpvec(1:n)=invec
      tmpvec(n+1:2*n)=invec(n:1:-1)
      tmpcmplx=cmplx(tmpvec,0.0_dp)
      dftvec=cmplx(0.0_dp,0.0_dp,dp)

      if(pub_on_root) then
         call linalg_1d_fft('F', 2*n, tmpcmplx, dftvec)
      end if
      call comms_bcast(pub_root_proc_id,dftvec)

      do i=2,n
         outvec(i)=(1.0_dp/sqrt(2.0_dp*real(n,dp)))*real(exp(cmplx(0.0_dp,-real(i-1,dp)*pi)/real(2*n,dp))*dftvec(i),dp)
      end do
      outvec(1)=real(dftvec(1),dp)/sqrt(4.0_dp/real(n,dp))

      outvec=outvec*(sqrt(2.0_dp*real(n,dp)))
      outvec(1)=outvec(1)*(sqrt(2.0_dp)/real(n,dp))

      deallocate(dftvec,stat=ierr)
      call utils_dealloc_check('get_tanh_chebyshev_coeffs_costran','dftvec',ierr)
      deallocate(tmpvec,stat=ierr)
      call utils_dealloc_check('get_tanh_chebyshev_coeffs_costran','tmpvec',ierr)
      deallocate(tmpcmplx,stat=ierr)
      call utils_dealloc_check('get_tanh_chebyshev_coeffs_costran','tmpcmplx',ierr)

    end subroutine cos_trans


  end subroutine get_tanh_chebyshev_coeffs


   !============================================================================!
  ! This routine gets the tanh Chebyshev coefficients to approximate the       !
  ! function between two interval limits to a given threshold. The algorithm   !
  ! Uses a discrete cosine transform to give the coefficients (which uses the  !
  ! 1D Fourier transform of the linalg module).
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   lower        (input)  : The lower limit of the interval                  !
  !   upper        (input)  : The upper limit of the interval                  !
  !   threshold    (input)  : The accuracy threshold of the Chebyshev          !
  !                            approximation on the interval                   !
  !   coeffs      (output)  : The Chebyshev coefficients                       !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Apr 2017,                                        !
  !============================================================================!

  subroutine get_hypmultid_chebyshev_coeffs(lower,upper,threshold,coeffs)
    use comms, only: pub_on_root
    use constants, only: dp
    use utils, only : utils_alloc_check, utils_dealloc_check
    implicit none
    real(kind=dp), intent(in) :: lower
    real(kind=dp), intent(in) :: upper
    real(kind=dp), intent(in) :: threshold
    real(kind=dp), dimension(:), allocatable, intent(inout) :: coeffs

    integer :: n
    integer :: ierr
    real(kind=dp) :: halfwidth, midpoint
    real(kind=dp) :: shift
    integer :: i


    if(allocated(coeffs)) then
       deallocate(coeffs,stat=ierr)
       call utils_dealloc_check('get_hypmultid_chebyshev_coeffs','coeffs',ierr)
    end if

    halfwidth=(upper-lower)/2.0_dp
    midpoint=(lower+upper)/2.0_dp

    n = howmany(lower,upper,threshold)
!n=337
    if(pub_on_root) then
       if(debug_info_toggle) then
          write(*,*) "n=",n," lower=",lower," upper=",upper
       end if
    end if
    if(n>0) then
       allocate(coeffs(n),stat=ierr)
       call utils_alloc_check('get_hypmultid_chebyshev_coeffs','coeffs',ierr)

       call hypmultid_chebyshev_coeffs(n,lower,upper,coeffs,shift)

    else
       allocate(coeffs(-n),stat=ierr)
       call utils_alloc_check('get_hypmultid_chebyshev_coeffs','coeffs',ierr)

       call hypmultid_chebyshev_coeffs(-n,lower,upper,coeffs,shift)

    end if

  contains

    ! A function to give the number of Chebyshev coefficients needed to represent
    ! the hypmultid function to a given accuracy threshold on the interval [lower, upper].
    function howmany(lower,upper,threshold) result(nterms)
      implicit none
      integer :: nterms
      real(kind=dp), intent(in) :: lower
      real(kind=dp), intent(in) :: upper
      real(kind=dp), intent(in) :: threshold

      real(kind=dp) :: halfwidth, midpoint
      real(kind=dp), dimension(:), allocatable :: coeffs
      integer :: i
      integer :: maxterms
      real(kind=dp) :: error, olderror
      integer :: wentup
      real(kind=dp) :: shift
      integer :: ierr

      halfwidth=(upper-lower)/2.0_dp
      midpoint=(lower+upper)/2.0_dp

      ! This is arbitrary... For a purely Head-Gordon scheme, we wouldn't want this limit.
      ! It is set to this because at 341, we achieve values of hypmultid(argument) which are
      ! indistinguishable from -1 and 1 to machine precision and can hence project eigs
      ! outside of this range.
!      maxterms=341
      maxterms=10000

      error=huge(1.0_dp)
      olderror=huge(1.0_dp)

      nterms=-1

      do i=3,maxterms,1
         allocate(coeffs(i),stat=ierr)
         call utils_alloc_check('get_hypmultid_chebyshev_coeffs_howmany','coeffs',ierr)

         call hypmultid_chebyshev_coeffs(i,lower,upper,coeffs,shift)
         olderror=error
         error=(abs(chebpolval(coeffs, (lower - midpoint)/halfwidth) - hypmultid(lower) + shift)+ &
              & abs(chebpolval(coeffs, (upper - midpoint)/halfwidth) - hypmultid(upper) + shift))/2.0_dp
!         write(*,*) i,": delta=",error, (chebpolval(coeffs, (lower - midpoint)/halfwidth) - hypmultid(lower) + shift), &
!              & (chebpolval(coeffs, (upper - midpoint)/halfwidth) - hypmultid(upper) + shift)
         if(error>olderror) then
            wentup=wentup+1
         else
            wentup=0
         end if
         if(wentup>=3) then
            ! failed to converge
            nterms=-i
            exit
         end if

         if(error<threshold) then
            nterms=i
            deallocate(coeffs,stat=ierr)
            call utils_dealloc_check('get_hypmultid_chebyshev_coeffs_howmany','coeffs',ierr)
            exit
         end if
         deallocate(coeffs,stat=ierr)
         call utils_dealloc_check('get_hypmultid_chebyshev_coeffs_howmany','coeffs',ierr)
      end do

    end function howmany

    ! This function evaluates a Chebyshev polynomial with coeffs c, for a given
    ! argument x (between -1 and 1).
    function chebpolval(c,x) result(u)
      implicit none
      real(kind=dp) :: u
      real(kind=dp), dimension(:), intent(in) :: c
      real(kind=dp),               intent(in) :: x
      integer :: n
      real(kind=dp) :: ujp1, ujp2
      integer :: j

      n = size(c)
      u = c(n)
      if(n > 1) then
         ujp1 = u
         u = c(n-1) + 2.0_dp*x*c(n)
         do j = n-2,1,-1
            ujp2 = ujp1
            ujp1 = u
            u = c(j) + 2.0_dp*x*ujp1 - ujp2
         end do
         u = u - x*ujp1
      end if
    end function chebpolval


    ! This routine gives Chebyshev coefficients which approximate hypmultid on an
    ! interval [start, finish]. The number of coeffs is requested using n.
    ! the shift  gives the constant shift between the approximation and the
    ! hypmultid function.
    subroutine hypmultid_chebyshev_coeffs(n,start,finish,alpha, shift)
      implicit none
      integer, intent(in) :: n
      real(kind=dp), intent(in) :: start
      real(kind=dp), intent(in) :: finish
      real(kind=dp), dimension(n), intent(out) :: alpha ! function coefficients
      real(kind=dp), intent(out) :: shift

      real(kind=dp),  dimension(n) :: y
      !    real(kind=dp), dimension(n) :: alpha ! function coefficients
      real(kind=dp), dimension(n) :: theta, costheta
      integer(kind=dp) :: plandct
      integer :: i
      real(kind=dp) :: halfwidth, midpoint

      halfwidth=(finish-start)/2.0_dp
      midpoint=(start+finish)/2.0_dp

      do i=1,n
         theta(i)=(real(i-1,dp)+0.5_dp)*pi/real(n,dp)
         costheta(i)=cos(theta(i))
         y(i)=hypmultid(midpoint + costheta(i)*halfwidth)
      end do

      ! DCT
      ! We could use a DCT directly by doing this:
      !    call dfftw_plan_r2r_1d(plandct, n, y, alpha, FFTW_REDFT10, FFTW_ESTIMATE)
      !    call dfftw_execute_r2r(plandct, y, alpha)
      !    call dfftw_destroy_plan(plandct)

      call cos_trans(y,alpha)

      alpha=(alpha/real(n,dp))
      shift=(hypmultid(midpoint)-chebpolval(alpha,0.0_dp))

    end subroutine hypmultid_chebyshev_coeffs


    ! An inefficient cosine transform using ONETEPs built-in 1d FFT.
    subroutine cos_trans(invec, outvec)
      use comms, only: pub_on_root, pub_root_proc_id,comms_bcast
      use linalg, only: linalg_1d_fft
      implicit none
      real(kind=dp),    dimension(:), intent(in)    :: invec
      real(kind=dp),    dimension(:), intent(out)   :: outvec
      integer :: n
      complex(kind=dp), dimension(:), allocatable   :: dftvec
      complex(kind=dp), dimension(:), allocatable   :: tmpcmplx
      real(kind=dp),    dimension(:), allocatable   :: tmpvec
      integer :: ierr
      integer :: i

      n=size(invec)
      allocate(dftvec(2*n),stat=ierr)
      call utils_alloc_check('get_hypmultid_chebyshev_coeffs_costran','dftvec',ierr)
      allocate(tmpvec(2*n),stat=ierr)
      call utils_alloc_check('get_hypmultid_chebyshev_coeffs_costran','tmpvec',ierr)
      allocate(tmpcmplx(2*n),stat=ierr)
      call utils_alloc_check('get_hypmultid_chebyshev_coeffs_costran','tmpcmplx',ierr)

      tmpvec(1:n)=invec
      tmpvec(n+1:2*n)=invec(n:1:-1)
      tmpcmplx=cmplx(tmpvec,0.0_dp)
      dftvec=cmplx(0.0_dp,0.0_dp,dp)

      if(pub_on_root) then
         call linalg_1d_fft('F', 2*n, tmpcmplx, dftvec)
      end if
      call comms_bcast(pub_root_proc_id,dftvec)

      do i=2,n
         outvec(i)=(1.0_dp/sqrt(2.0_dp*real(n,dp)))*real(exp(cmplx(0.0_dp,-real(i-1,dp)*pi)/real(2*n,dp))*dftvec(i),dp)
      end do
      outvec(1)=real(dftvec(1),dp)/sqrt(4.0_dp/real(n,dp))

      outvec=outvec*(sqrt(2.0_dp*real(n,dp)))
      outvec(1)=outvec(1)*(sqrt(2.0_dp)/real(n,dp))

      deallocate(dftvec,stat=ierr)
      call utils_dealloc_check('get_hypmultid_chebyshev_coeffs_costran','dftvec',ierr)
      deallocate(tmpvec,stat=ierr)
      call utils_dealloc_check('get_hypmultid_chebyshev_coeffs_costran','tmpvec',ierr)
      deallocate(tmpcmplx,stat=ierr)
      call utils_dealloc_check('get_hypmultid_chebyshev_coeffs_costran','tmpcmplx',ierr)

    end subroutine cos_trans

    real(kind=dp) function hypmultid(x)
      implicit none
      real(kind=dp) :: x

      hypmultid = (2.0_dp*x)/(x**2 + 1.0_dp)

    end function hypmultid


  end subroutine get_hypmultid_chebyshev_coeffs



     !============================================================================!
  ! This routine gets the tanh Chebyshev coefficients to approximate the       !
  ! function between two interval limits to a given threshold. The algorithm   !
  ! Uses a discrete cosine transform to give the coefficients (which uses the  !
  ! 1D Fourier transform of the linalg module).
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   lower        (input)  : The lower limit of the interval                  !
  !   upper        (input)  : The upper limit of the interval                  !
  !   threshold    (input)  : The accuracy threshold of the Chebyshev          !
  !                            approximation on the interval                   !
  !   coeffs      (output)  : The Chebyshev coefficients                       !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Apr 2017,                                        !
  !============================================================================!

  subroutine get_hypsumid_chebyshev_coeffs(lower,upper,threshold,delta,coeffs)
    use comms, only: pub_on_root
    use constants, only: dp
    use utils, only : utils_alloc_check, utils_dealloc_check
    implicit none
    real(kind=dp), intent(in) :: lower
    real(kind=dp), intent(in) :: upper
    real(kind=dp), intent(in) :: threshold
    real(kind=dp), intent(in) :: delta
    real(kind=dp), dimension(:), allocatable, intent(inout) :: coeffs

    integer :: n
    integer :: ierr
    real(kind=dp) :: halfwidth, midpoint
    real(kind=dp) :: shift
    integer :: i


    if(allocated(coeffs)) then
       deallocate(coeffs,stat=ierr)
       call utils_dealloc_check('get_hypsumid_chebyshev_coeffs','coeffs',ierr)
    end if

    halfwidth=(upper-lower)/2.0_dp
    midpoint=(lower+upper)/2.0_dp

    n = howmany(lower,upper,threshold,delta)
    !n=337,
    if(pub_on_root) then
       if(debug_info_toggle) then
          write(*,*) "n=",n," lower=",lower," upper=",upper
       end if
    end if

    if(n>0) then
       allocate(coeffs(n),stat=ierr)
       call utils_alloc_check('get_hypsumid_chebyshev_coeffs','coeffs',ierr)

       call hypsumid_chebyshev_coeffs(n,lower,upper,delta,coeffs,shift)

    else
       allocate(coeffs(-n),stat=ierr)
       call utils_alloc_check('get_hypsumid_chebyshev_coeffs','coeffs',ierr)

       call hypsumid_chebyshev_coeffs(-n,lower,upper,delta,coeffs,shift)

    end if

    if(abs(n)>0) then
       coeffs(1)=coeffs(1)+shift
    end if

  contains

    ! A function to give the number of Chebyshev coefficients needed to represent
    ! the hypsumid function to a given accuracy threshold on the interval [lower, upper].
    function howmany(lower,upper,threshold,delta) result(nterms)
      implicit none
      integer :: nterms
      real(kind=dp), intent(in) :: lower
      real(kind=dp), intent(in) :: upper
      real(kind=dp), intent(in) :: threshold
      real(kind=dp), intent(in) :: delta

      real(kind=dp) :: halfwidth, midpoint
      real(kind=dp), dimension(:), allocatable :: coeffs
      integer :: i
      integer :: maxterms
      real(kind=dp) :: error, olderror
      integer :: wentup
      real(kind=dp) :: shift
      integer :: ierr

      halfwidth=(upper-lower)/2.0_dp
      midpoint=(lower+upper)/2.0_dp

      ! This is arbitrary... For a purely Head-Gordon scheme, we wouldn't want this limit.
      ! It is set to this because at 341, we achieve values of hypsumid(argument) which are
      ! indistinguishable from -1 and 1 to machine precision and can hence project eigs
      ! outside of this range.
!      maxterms=341
      maxterms=80

      error=huge(1.0_dp)
      olderror=huge(1.0_dp)

      nterms=-1

      do i=3,maxterms,1
         allocate(coeffs(i),stat=ierr)
         call utils_alloc_check('get_hypsumid_chebyshev_coeffs_howmany','coeffs',ierr)

         call hypsumid_chebyshev_coeffs(i,lower,upper,delta,coeffs,shift)
         olderror=error
         error=(abs(chebpolval(coeffs, (lower - midpoint)/halfwidth) - hypsumid(lower,delta) + shift)+ &
              & abs(chebpolval(coeffs, (upper - midpoint)/halfwidth) - hypsumid(upper,delta) + shift))/2.0_dp
!         write(*,*) i,": delta=",error, (chebpolval(coeffs, (lower - midpoint)/halfwidth) - hypsumid(lower,delta) + shift), &
!              & (chebpolval(coeffs, (upper - midpoint)/halfwidth) - hypsumid(upper,delta) + shift)
         if(error>olderror) then
            wentup=wentup+1
         else
            wentup=0
         end if
         if(wentup>=3) then
            ! failed to converge
            nterms=-i
            exit
         end if

         if(error<threshold) then
            nterms=i
            deallocate(coeffs,stat=ierr)
            call utils_dealloc_check('get_hypsumid_chebyshev_coeffs_howmany','coeffs',ierr)
            exit
         end if
         deallocate(coeffs,stat=ierr)
         call utils_dealloc_check('get_hypsumid_chebyshev_coeffs_howmany','coeffs',ierr)
      end do

      if(nterms==-1) then
         nterms=maxterms
      end if

    end function howmany

    ! This function evaluates a Chebyshev polynomial with coeffs c, for a given
    ! argument x (between -1 and 1).
    function chebpolval(c,x) result(u)
      implicit none
      real(kind=dp) :: u
      real(kind=dp), dimension(:), intent(in) :: c
      real(kind=dp),               intent(in) :: x
      integer :: n
      real(kind=dp) :: ujp1, ujp2
      integer :: j

      n = size(c)
      u = c(n)
      if(n > 1) then
         ujp1 = u
         u = c(n-1) + 2.0_dp*x*c(n)
         do j = n-2,1,-1
            ujp2 = ujp1
            ujp1 = u
            u = c(j) + 2.0_dp*x*ujp1 - ujp2
         end do
         u = u - x*ujp1
      end if
    end function chebpolval


    ! This routine gives Chebyshev coefficients which approximate hypsumid on an
    ! interval [start, finish]. The number of coeffs is requested using n.
    ! the shift  gives the constant shift between the approximation and the
    ! hypsumid function.
    subroutine hypsumid_chebyshev_coeffs(n,start,finish,delta,alpha, shift)
      implicit none
      integer, intent(in) :: n
      real(kind=dp), intent(in) :: start
      real(kind=dp), intent(in) :: finish
      real(kind=dp), intent(in) :: delta
      real(kind=dp), dimension(n), intent(out) :: alpha ! function coefficients
      real(kind=dp), intent(out) :: shift

      real(kind=dp),  dimension(n) :: y
      !    real(kind=dp), dimension(n) :: alpha ! function coefficients
      real(kind=dp), dimension(n) :: theta, costheta
      integer(kind=dp) :: plandct
      integer :: i
      real(kind=dp) :: halfwidth, midpoint

      halfwidth=(finish-start)/2.0_dp
      midpoint=(start+finish)/2.0_dp

      do i=1,n
         theta(i)=(real(i-1,dp)+0.5_dp)*pi/real(n,dp)
         costheta(i)=cos(theta(i))
         y(i)=hypsumid(midpoint + costheta(i)*halfwidth,delta)
      end do

      ! DCT
      ! We could use a DCT directly by doing this:
      !    call dfftw_plan_r2r_1d(plandct, n, y, alpha, FFTW_REDFT10, FFTW_ESTIMATE)
      !    call dfftw_execute_r2r(plandct, y, alpha)
      !    call dfftw_destroy_plan(plandct)

      call cos_trans(y,alpha)

      alpha=(alpha/real(n,dp))
      shift=(hypsumid(midpoint,delta)-chebpolval(alpha,0.0_dp))

    end subroutine hypsumid_chebyshev_coeffs


    ! An inefficient cosine transform using ONETEPs built-in 1d FFT.
    subroutine cos_trans(invec, outvec)
      use comms, only: pub_on_root, pub_root_proc_id,comms_bcast
      use linalg, only: linalg_1d_fft
      implicit none
      real(kind=dp),    dimension(:), intent(in)    :: invec
      real(kind=dp),    dimension(:), intent(out)   :: outvec
      integer :: n
      complex(kind=dp), dimension(:), allocatable   :: dftvec
      complex(kind=dp), dimension(:), allocatable   :: tmpcmplx
      real(kind=dp),    dimension(:), allocatable   :: tmpvec
      integer :: ierr
      integer :: i

      n=size(invec)
      allocate(dftvec(2*n),stat=ierr)
      call utils_alloc_check('get_hypsumid_chebyshev_coeffs_costran','dftvec',ierr)
      allocate(tmpvec(2*n),stat=ierr)
      call utils_alloc_check('get_hypsumid_chebyshev_coeffs_costran','tmpvec',ierr)
      allocate(tmpcmplx(2*n),stat=ierr)
      call utils_alloc_check('get_hypsumid_chebyshev_coeffs_costran','tmpcmplx',ierr)

      tmpvec(1:n)=invec
      tmpvec(n+1:2*n)=invec(n:1:-1)
      tmpcmplx=cmplx(tmpvec,0.0_dp)
      dftvec=cmplx(0.0_dp,0.0_dp,dp)

      if(pub_on_root) then
         call linalg_1d_fft('F', 2*n, tmpcmplx, dftvec)
      end if
      call comms_bcast(pub_root_proc_id,dftvec)

      do i=2,n
         outvec(i)=(1.0_dp/sqrt(2.0_dp*real(n,dp)))*real(exp(cmplx(0.0_dp,-real(i-1,dp)*pi)/real(2*n,dp))*dftvec(i),dp)
      end do
      outvec(1)=real(dftvec(1),dp)/sqrt(4.0_dp/real(n,dp))

      outvec=outvec*(sqrt(2.0_dp*real(n,dp)))
      outvec(1)=outvec(1)*(sqrt(2.0_dp)/real(n,dp))

      deallocate(dftvec,stat=ierr)
      call utils_dealloc_check('get_hypsumid_chebyshev_coeffs_costran','dftvec',ierr)
      deallocate(tmpvec,stat=ierr)
      call utils_dealloc_check('get_hypsumid_chebyshev_coeffs_costran','tmpvec',ierr)
      deallocate(tmpcmplx,stat=ierr)
      call utils_dealloc_check('get_hypsumid_chebyshev_coeffs_costran','tmpcmplx',ierr)

    end subroutine cos_trans

    real(kind=dp) function hypsumid(x,y)
      implicit none
      real(kind=dp) :: x,y

      hypsumid = (x+tanh(y))/(1.0_dp + (x*tanh(y)))

    end function hypsumid


  end subroutine get_hypsumid_chebyshev_coeffs





   !============================================================================!
  ! This routine gets the tanh Chebyshev coefficients to approximate the       !
  ! function between two interval limits to a given threshold. The algorithm   !
  ! Uses a discrete cosine transform to give the coefficients (which uses the  !
  ! 1D Fourier transform of the linalg module).
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   lower        (input)  : The lower limit of the interval                  !
  !   upper        (input)  : The upper limit of the interval                  !
  !   threshold    (input)  : The accuracy threshold of the Chebyshev          !
  !                            approximation on the interval                   !
  !   coeffs      (output)  : The Chebyshev coefficients                       !
  !   expans      (input)   : Chebyshev expansion used if defined              !
  !   dcoeffs     (output)  : The Chebyshev coefficients  for the derivative   !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Apr 2017,                                        !
  ! Modified by Emiliano Poli in 2022 to work with the Mermin Method           !
  !============================================================================!

  subroutine get_approxentropy_chebyshev_coeffs(lower,upper,threshold,coeffs, &
    expans,dcoeffs)
    use comms, only: pub_on_root
    use constants, only: dp, PI
    use rundat, only: pub_mermin
    use utils, only : utils_alloc_check, utils_dealloc_check
    implicit none
    real(kind=dp), intent(in) :: lower
    real(kind=dp), intent(in) :: upper
    real(kind=dp), intent(in) :: threshold
    real(kind=dp), dimension(:), allocatable, intent(inout) :: coeffs
    !ep : mermin variables
    integer      , optional, intent(in) :: expans
    real(kind=dp), optional, dimension(:), allocatable, intent(inout) :: dcoeffs


    integer :: n
    integer :: ierr
    real(kind=dp) :: halfwidth, midpoint
    real(kind=dp) :: shift, shiftb
    integer :: i


    if(allocated(coeffs)) then
       deallocate(coeffs,stat=ierr)
       call utils_dealloc_check('get_approxentropy_chebyshev_coeffs','coeffs',ierr)
    end if

    halfwidth=(upper-lower)/2.0_dp
    midpoint=(lower+upper)/2.0_dp

    !ep
    if (pub_mermin) then
       if (present(expans)) then
          n=expans
       else
          !ep: default expansion in Mermin mode
          n=11
       end if
    else
       n = howmany(lower,upper,threshold)
    end if

    if(pub_on_root) then
       if(debug_info_toggle) then
          write(*,*) "n=",n," lower=",lower," upper=",upper
       end if
    end if
    if(n>0) then
       allocate(coeffs(n),stat=ierr)
       call utils_alloc_check('get_approxentropy_chebyshev_coeffs','coeffs',ierr)

       call approxentropy_chebyshev_coeffs(n,lower,upper,coeffs,shift)

    else
       allocate(coeffs(-n),stat=ierr)
       call utils_alloc_check('get_approxentropy_chebyshev_coeffs','coeffs',ierr)
       call approxentropy_chebyshev_coeffs(-n,lower,upper,coeffs,shift)
    end if

    coeffs(1)=coeffs(1)+shift

    ! ep: calculate Chebyshev derivatove coeffs if requested
    if (pub_mermin) then
       if(present(dcoeffs)) then
          if(allocated(dcoeffs)) then
             deallocate(dcoeffs,stat=ierr)
             call utils_dealloc_check('get_approxderiv_chebyshev_coeffs', &
                  'coeffs',ierr)
          end if
          allocate(dcoeffs(n),stat=ierr)
          call utils_alloc_check('get_approxentropy_chebyshev_coeffs','dcoeffs',ierr)
          call approxderiv_chebyshev_coeffs(n,lower,upper,coeffs,dcoeffs,shiftb)
          dcoeffs(1)=dcoeffs(1)+shiftb
       end if
    end if
    ! ep: calculate Chebyshev derivatove coeffs if requested

  contains

    ! A function to give the number of Chebyshev coefficients needed to represent
    ! the approxentropy function to a given accuracy threshold on the interval [lower, upper].
    function howmany(lower,upper,threshold) result(nterms)
      implicit none
      integer :: nterms
      real(kind=dp), intent(in) :: lower
      real(kind=dp), intent(in) :: upper
      real(kind=dp), intent(in) :: threshold

      real(kind=dp) :: halfwidth, midpoint
      real(kind=dp), dimension(:), allocatable :: coeffs
      integer :: i
      integer :: maxterms
      real(kind=dp) :: error, olderror
      integer :: wentup
      real(kind=dp) :: shift
      integer :: ierr

      halfwidth=(upper-lower)/2.0_dp
      midpoint=(lower+upper)/2.0_dp

      ! This is arbitrary... For a purely Head-Gordon scheme, we wouldn't want this limit.
      ! It is set to this because at 341, we achieve values of approxentropy(argument) which are
      ! indistinguishable from -1 and 1 to machine precision and can hence project eigs
      ! outside of this range.
!      maxterms=341
      maxterms=1000

      error=huge(1.0_dp)
      olderror=huge(1.0_dp)

      nterms=-1

      do i=3,maxterms,1
         allocate(coeffs(i),stat=ierr)
         call utils_alloc_check('get_approxentropy_chebyshev_coeffs_howmany','coeffs',ierr)

         call approxentropy_chebyshev_coeffs(i,lower,upper,coeffs,shift)
         olderror=error
         error=(abs(chebpolval(coeffs, (lower - midpoint)/halfwidth) - approxentropy(lower) + shift)+ &
              & abs(chebpolval(coeffs, (upper - midpoint)/halfwidth) - approxentropy(upper) + shift))/2.0_dp
!         write(*,*) i,": delta=",error, (chebpolval(coeffs, (lower - midpoint)/halfwidth) - approxentropy(lower) + shift), &
!              & (chebpolval(coeffs, (upper - midpoint)/halfwidth) - approxentropy(upper) + shift)
         if(error>olderror) then
            wentup=wentup+1
         else
            wentup=0
         end if
         if(wentup>=3) then
            ! failed to converge
            nterms=-i
            exit
         end if

         if(error<threshold) then
            nterms=i
            deallocate(coeffs,stat=ierr)
            call utils_dealloc_check('get_approxentropy_chebyshev_coeffs_howmany','coeffs',ierr)
            exit
         end if
         deallocate(coeffs,stat=ierr)
         call utils_dealloc_check('get_approxentropy_chebyshev_coeffs_howmany','coeffs',ierr)
      end do

    end function howmany

    ! This function evaluates a Chebyshev polynomial with coeffs c, for a given
    ! argument x (between -1 and 1).
    function chebpolval(c,x) result(u)
      implicit none
      real(kind=dp) :: u
      real(kind=dp), dimension(:), intent(in) :: c
      real(kind=dp),               intent(in) :: x
      integer :: n
      real(kind=dp) :: ujp1, ujp2
      integer :: j

      n = size(c)
      u = c(n)
      if(n > 1) then
         ujp1 = u
         u = c(n-1) + 2.0_dp*x*c(n)
         do j = n-2,1,-1
            ujp2 = ujp1
            ujp1 = u
            u = c(j) + 2.0_dp*x*ujp1 - ujp2
         end do
         u = u - x*ujp1
      end if
    end function chebpolval


    ! This routine gives Chebyshev coefficients which approximate approxentropy on an
    ! interval [start, finish]. The number of coeffs is requested using n.
    ! the shift  gives the constant shift between the approximation and the
    ! approxentropy function.
    subroutine approxentropy_chebyshev_coeffs(n,start,finish,alpha, shift)
      implicit none
      integer, intent(in) :: n
      real(kind=dp), intent(in) :: start
      real(kind=dp), intent(in) :: finish
      real(kind=dp), dimension(n), intent(out) :: alpha ! function coefficients
      real(kind=dp), intent(out) :: shift

      real(kind=dp),  dimension(n) :: y
      !    real(kind=dp), dimension(n) :: alpha ! function coefficients
      real(kind=dp), dimension(n) :: theta, costheta
      integer(kind=dp) :: plandct
      integer :: i
      real(kind=dp) :: halfwidth, midpoint

      halfwidth=(finish-start)/2.0_dp
      midpoint=(start+finish)/2.0_dp

      do i=1,n
         theta(i)=(real(i-1,dp)+0.5_dp)*pi/real(n,dp)
         costheta(i)=cos(theta(i))
         y(i)=approxentropy(midpoint + costheta(i)*halfwidth)
      end do

      ! DCT
      ! We could use a DCT directly by doing this:
      !    call dfftw_plan_r2r_1d(plandct, n, y, alpha, FFTW_REDFT10, FFTW_ESTIMATE)
      !    call dfftw_execute_r2r(plandct, y, alpha)
      !    call dfftw_destroy_plan(plandct)

      call cos_trans(y,alpha)

      alpha=(alpha/real(n,dp))
      shift=(approxentropy(midpoint)-chebpolval(alpha,0.0_dp))

    end subroutine approxentropy_chebyshev_coeffs

    !ep: This routine gives Chebyshev coefficients which approximate the derivative
    ! of the approxiamte entropy on an interval [start, finish]. i
    ! The number of coeffs is requested using n.
    ! The shift  gives the constant shift between the approximation and the
    ! approxentropy function.
    subroutine approxderiv_chebyshev_coeffs(n,start,finish,alpha,beta,shift)
      implicit none
      integer, intent(in) :: n
      real(kind=dp), intent(in) :: start
      real(kind=dp), intent(in) :: finish
      real(kind=dp), dimension(n), intent(in) :: alpha ! cheb function coefficients
      real(kind=dp), dimension(n), intent(out) :: beta ! cheb deriv function coefficients
      integer :: i, j
      real(kind=dp) :: con
      real(kind=dp) :: halfwidth, midpoint
      real(kind=dp), intent(out) :: shift

      halfwidth=(finish-start)/2.0_dp
      midpoint=(start+finish)/2.0_dp

      beta(n)=0
      beta(n-1)=2*(n-1)*alpha(n)
      do j=n-2,1,-1
         beta(j)=beta(j+2)+2*j*alpha(j+1)
      enddo
      con=2./(finish-start)
      do j=1,n
         beta(j)=beta(j)*con
      enddo
      shift=(approxderiv(0.5_dp)-chebpolval(beta,0.0_dp))

    end subroutine approxderiv_chebyshev_coeffs
    !ep

    ! An inefficient cosine transform using ONETEPs built-in 1d FFT.
    subroutine cos_trans(invec, outvec)
      use comms, only: pub_on_root, pub_root_proc_id,comms_bcast
      use linalg, only: linalg_1d_fft
      implicit none
      real(kind=dp),    dimension(:), intent(in)    :: invec
      real(kind=dp),    dimension(:), intent(out)   :: outvec
      integer :: n
      complex(kind=dp), dimension(:), allocatable   :: dftvec
      complex(kind=dp), dimension(:), allocatable   :: tmpcmplx
      real(kind=dp),    dimension(:), allocatable   :: tmpvec
      integer :: ierr
      integer :: i

      n=size(invec)
      allocate(dftvec(2*n),stat=ierr)
      call utils_alloc_check('get_approxentropy_chebyshev_coeffs_costran','dftvec',ierr)
      allocate(tmpvec(2*n),stat=ierr)
      call utils_alloc_check('get_approxentropy_chebyshev_coeffs_costran','tmpvec',ierr)
      allocate(tmpcmplx(2*n),stat=ierr)
      call utils_alloc_check('get_approxentropy_chebyshev_coeffs_costran','tmpcmplx',ierr)

      tmpvec(1:n)=invec
      tmpvec(n+1:2*n)=invec(n:1:-1)
      tmpcmplx=cmplx(tmpvec,0.0_dp)
      dftvec=cmplx(0.0_dp,0.0_dp,dp)

      if(pub_on_root) then
         call linalg_1d_fft('F', 2*n, tmpcmplx, dftvec)
      end if
      call comms_bcast(pub_root_proc_id,dftvec)

      do i=2,n
         outvec(i)=(1.0_dp/sqrt(2.0_dp*real(n,dp)))*real(exp(cmplx(0.0_dp,-real(i-1,dp)*pi)/real(2*n,dp))*dftvec(i),dp)
      end do
      outvec(1)=real(dftvec(1),dp)/sqrt(4.0_dp/real(n,dp))

      outvec=outvec*(sqrt(2.0_dp*real(n,dp)))
      outvec(1)=outvec(1)*(sqrt(2.0_dp)/real(n,dp))

      deallocate(dftvec,stat=ierr)
      call utils_dealloc_check('get_approxentropy_chebyshev_coeffs_costran','dftvec',ierr)
      deallocate(tmpvec,stat=ierr)
      call utils_dealloc_check('get_approxentropy_chebyshev_coeffs_costran','tmpvec',ierr)
      deallocate(tmpcmplx,stat=ierr)
      call utils_dealloc_check('get_approxentropy_chebyshev_coeffs_costran','tmpcmplx',ierr)

    end subroutine cos_trans

    real(kind=dp) function approxentropy(x)
      use utils, only : utils_assert
      implicit none
      real(kind=dp) :: x,z

      !ep: z shift not used in mermin method
      !ep: we use the real function
      if (.not. pub_mermin) then
         z=(x/2.0_dp)+0.5_dp
      end if
      !ep

      !ep
      if (pub_mermin) then
         call utils_assert(x>0.0, 'Error in approxentropy function:&
              & x has to be positive')
         call utils_assert((1-x)>0.0, 'Error in approxentropy function:&
              & 1-x has to be positive')
         approxentropy = x*log(x)+(1-x)*log(1-x)
      else
!        approxentropy = 1.96056_dp*(x**2) + 0.0286723_dp/(0.114753_dp + 1.9888_dp*x -1.9888_dp*(x**2)) - 0.24986_dp -1.96056_dp*x
         approxentropy = 1.96056_dp*(z**2) + 0.0286723_dp/(0.114753_dp + 1.9888_dp*z -1.9888_dp*(z**2)) - 0.24986_dp -1.96056_dp*z
!        approxentropy = (2.0_dp*x)/(x**2 + 1.0_dp)
      end if

    end function approxentropy

    !ep
    real(kind=dp) function approxderiv(x)
      use utils, only : utils_assert
      implicit none
      real(kind=dp) :: x,z

      call utils_assert(x>0.0, 'Error in approxderiv function:&
           & x has to be positive')
      call utils_assert((1-x)>0.0, 'Error in approxderiv function:&
           & 1-x has to be positive')
      approxderiv = log(x)-log(1-x)

    end function approxderiv
    !ep

  end subroutine get_approxentropy_chebyshev_coeffs


  !============================================================================!
  ! To caclulate the sign function of a matrix, iteratively and without        !
  ! diagonalisation, a modified Newton's scheme can be employed much like the  !
  ! one in this module used for calculating the matrix square-root iteratively !
  ! for use in the Lowdin transformation. In particular, this routine makes    !
  ! use of the accelerated, stabilised update of [Chen and Chow].              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   A            (input)  : The matrix for which the sign function will be   !
  !                           computed.                                        !
  !   tol          (input)  : The tolerance on the 2-norm stopping criterion.  !
  !   X            (output) : The computed sign-function matrix                !
  !   matmuls      (output) : The number of matrix-matrix operations required  !
  !                           to apply the function.                           !
  !   status       (output) : Success = 0.
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, May 2014,                                        !
  !----------------------------------------------------------------------------!
  ! This routine scales as O(log(N)) in matrix-matrix multiplications.         !
  !============================================================================!
  subroutine sign_nsh(A,tol,X,abs_emax,matmuls,status)
    use comms, only : pub_on_root
    use timer, only : timer_clock
    use utils, only : utils_alloc_check, utils_dealloc_check
    implicit none
    type(smearing_matrix),  intent(inout) :: A
    real(kind=dp), intent(in) :: tol
    type(smearing_matrix),  intent(inout) :: X ! Inout to preserve allocation
    real(kind=dp), intent(inout) :: abs_emax
    integer,       intent(out) :: matmuls
    integer,       intent(out) :: status

    type(smearing_matrix) :: tmp1, tmp2
    type(smearing_matrix) :: X0, X2, X0_bkp
    real(kind=dp) :: emin, emax, y, alpha, alphahat
    real(kind=dp) :: eigtol
    real(kind=dp) :: mnorm
    real(kind=dp), dimension(:), allocatable :: mnorm_history
    integer :: mnorm_history_length, mnorm_history_num, mnorm_history_maxlength
    integer :: i
    integer :: low_stat
    integer :: low_maxiters, low_iter
    integer :: ierr

    real(kind=dp) :: coeffs(2)
    integer :: num_positive_grads, num_positive_to_stop
    real(kind=dp) :: t

    call timer_clock('sign_nsh',1)

    status=0

    matmuls=0

    t=1.0_dp
    call internal_cubic_roots([-t, 3.0_dp/2.0_dp, 0.0_dp, -1.0_dp/2.0_dp],alphahat)
!    write(*,*) alphahat

!    alphahat = 1.697702485255767507654391_dp

    mnorm_history_maxlength=5
    mnorm_history_length=0
    mnorm_history_num=0
    num_positive_grads=0
    num_positive_to_stop=4

    allocate(mnorm_history(mnorm_history_maxlength), stat=ierr)
    call utils_alloc_check('sign_nsh','mnorm_history',ierr)
    mnorm_history=0.0_dp

    !    if(matrix_dimension(A)<10000.and..false.) then ! Always true until the matrix vector multiplication can be improved.
    if(.false.) then
       call matrix_allocate(tmp1,A)
       low_maxiters=100
       call lowdin_transformation(A,inv_sqrt_S_out=tmp1,niter=low_iter,maxiters=low_maxiters,status=low_stat)
       if(low_stat/=0) then
          status=-10
          call matrix_free(tmp1)
          call timer_clock('sign_nsh',2)
          return
       end if
       matmuls=matmuls+(low_iter*3)
       call matrix_allocate(tmp2,tmp1,tmp1)
       call matrix_multiply(1.0_dp,tmp1,tmp1,0.0_dp,tmp2)
       matmuls=matmuls+1
       call matrix_free(tmp1)
       call matrix_multiply(1.0_dp,A,tmp2,0.0_dp,X)
       matmuls=matmuls+1
       call matrix_free(tmp2)

    else

       eigtol=1e-5_dp

       !ja531-> tmp2 is just used here as a temporary copy of A.
       call matrix_allocate(tmp2,A)
       call matrix_copy(A,tmp2)
       !ja531-> sparsefoe:
!       call matrix_allocate(tmp1,A)
       call matrix_allocate(tmp1,X)

       call matrix_multiply(1.0_dp,A,tmp2,0.0_dp,tmp1)
       matmuls=matmuls+1
       call matrix_free(tmp2)
       ! Use the Gerschgorin approximation (via 1-norm) to estimate maximum eigenvalue.
       if(abs_emax>0.0_dp) then
          emax=abs_emax
       else
          abs_emax=sqrt(matrix_induced_norm(tmp1))
          emax=abs_emax
       end if

       !ja531-> sparsefoe:
!       call matrix_allocate(X0,A)
       call matrix_allocate(X0,X,X)
       call matrix_allocate(X0_bkp,X,X)

       call matrix_copy(A,X0)
       call matrix_scale(X0,1.0_dp/emax)
       !ja531-> sparsefoe:
!       call matrix_allocate(X2,A)
       call matrix_allocate(X2,X,X)

       call matrix_copy(X0,tmp1)
       call matrix_multiply(1.0_dp,X0,tmp1,0.0_dp,X2)


       call matrix_free(tmp1)
       call matrix_allocate(tmp1,X)

       matmuls=matmuls+1
       !       y=emin/emax

       y=1e-6_dp

       i=0
       do
          i=i+1

          call matrix_copy(X0,X0_bkp)

          alpha=min(sqrt(3.0_dp/(1.0_dp + y + y**2)),alphahat)
          call matrix_copy(X2,tmp1)
          call matrix_scale(tmp1,-(alpha**2),3.0_dp)
          call matrix_multiply(0.5_dp*alpha,X0,tmp1,0.0_dp,X0)
          matmuls=matmuls+1
          y=0.5_dp*alpha*y*(3.0_dp-(alpha**2)*y**2)
          call matrix_copy(X0,tmp1)
          call matrix_multiply(1.0_dp,X0,tmp1,0.0_dp,X2)
          matmuls=matmuls+1
          call matrix_copy(X2,tmp1)
          call matrix_scale(tmp1,1.0_dp,-1.0_dp)
          mnorm=matrix_norm(tmp1,2)

!          mnorm_history_num=mnorm_history_num+1
!          if(mnorm_history_num>mnorm_history_maxlength) mnorm_history_num=1
!          if(mnorm_history_length<mnorm_history_maxlength) mnorm_history_length=mnorm_history_length+1

!          mnorm_history(1:mnorm_history_maxlength-1)=mnorm_history(2:mnorm_history_maxlength)
!          mnorm_history(mnorm_history_length)=mnorm


          mnorm_history_length=mnorm_history_length+1
          if(mnorm_history_length<=mnorm_history_maxlength) then
             mnorm_history(mnorm_history_length)=mnorm
             call poly_fit(1, mnorm_history, coeffs, mnorm_history_length)
          else
             mnorm_history(1:mnorm_history_maxlength-1)=mnorm_history(2:mnorm_history_maxlength)
             mnorm_history(mnorm_history_maxlength)=mnorm
             call poly_fit(1, mnorm_history, coeffs)
          end if

          if(coeffs(2)>0.0_dp) then
             num_positive_grads=num_positive_grads+1
          else
             num_positive_grads=0
          end if

          if(mnorm_history_length>1) then
             if(mnorm_history(mnorm_history_maxlength)>mnorm_history(mnorm_history_maxlength-1)) then
                num_positive_grads=num_positive_grads+1

                if(pub_on_root) then
                   if(mod(i,1)==0 .and. debug_info_toggle) write(stdout,*) "smearing-> ",mnorm, i, tol, "<== SMEAR"
                end if
                call matrix_copy(X0_bkp,X)
                exit
             end if
          end if

          if(pub_on_root) then
             if(mod(i,1)==0 .and. debug_info_toggle) write(stdout,*) "smearing-> ",mnorm, i, tol, "<== SMEAR"
          end if
          if(mnorm<tol) then
             call matrix_copy(X0,X)
             exit
          elseif(mnorm/=mnorm) then ! something's gone very wrong...
             status=-2
             if(pub_on_root) then
                if(debug_info_toggle) write(stdout,*) "smearing-> NaNreturn :::: ", emin,emax,mnorm, "<== SMEAR"
             end if
             exit
          elseif(num_positive_grads>num_positive_to_stop) then
             status=-3
             if(pub_on_root) then
                if(debug_info_toggle) write(stdout,*) "smearing-> Gradient went up, so stopping NSH. Final norm:",mnorm,"<== SMEAR"
             end if
             exit
          end if
       end do
       call matrix_free(X0)
       call matrix_free(X0_bkp)
       call matrix_free(X2)
       call matrix_free(tmp1)

!       if(pub_on_root) then
!          do i=1,mnorm_history_maxlength
!             write(stdout,*) "mnorm_history:",i,"->",mnorm_history(i)
!          end do
!       end if

    end if


    deallocate(mnorm_history, stat=ierr)
    call utils_dealloc_check('sign_nsh','mnorm_history',ierr)


    call timer_clock('sign_nsh',2)
contains

  ! Adapted from internal_quartic_min in kernel_mod
    subroutine internal_cubic_roots(a,x)

      use constants, only: DP, PI
      use utils, only: utils_assert

      implicit none

      ! Arguments
      real(kind=DP), intent(in) :: a(4)
      real(kind=DP), intent(out) :: x

      ! Local variables
      integer :: i
      logical :: foundmin
      real(kind=DP), parameter :: THIRD = 1.0_DP / 3.0_DP
      real(kind=DP) :: aa,bb,cc
      real(kind=DP) :: q,r,theta,rtq,y,z
      real(kind=DP) :: x1(3),q1(3),qmin, delta

      ! If a4 vanishes then the gradient must be zero in which case choose
      ! step length to correspond to purification transformation
!      if (abs(a(4)) < epsilon(1.0_DP)) then
!         x = -0.5_DP
!         return
!      end if

      ! write result as  x**3 + a*x**2 + b*x + c = 0
      aa =  a(3) / a(4)
      bb = a(2) / a(4)
      cc =  a(1) / a(4)


      ! Solve cubic equation
      q = (aa*aa - 3.0_DP*bb)/9.0_DP !delta0
      r = (2*aa*aa*aa - 9.0_DP*aa*bb + 27.0_DP*cc)/54.0_DP !delta1
      if (r*r < q*q*q) then   ! three real roots
         rtq = sqrt(q)
         theta = acos(r/(rtq*q))
         x1(1) = -2.0_DP * rtq * cos(theta*THIRD) - aa*THIRD
         x1(2) = -2.0_DP * rtq * cos((theta+2.0_DP*PI)*THIRD) - aa*THIRD
         x1(3) = -2.0_DP * rtq * cos((theta-2.0_DP*PI)*THIRD) - aa*THIRD

! write(*,*) x1
         foundmin = .false.
         qmin = 0.0_DP ! qoh: Initialise to prevent compiler warning
         do i=1,3
            if((x1(i)<=sqrt(3.0_dp)).and.(x1(i)>=1.0_dp)) then
               x=x1(i)
               foundmin = .true.
            end if
         end do
         call utils_assert(foundmin, 'Error in internal_quartic_min &
              &(kernel_mod.F90): no minimum found')
      else   ! only one root -> must be a minimum
         delta=((r*54.0_dp)**2 - 4.0_dp*(q*9.0_dp)**3)/(-27.0_dp)
         if(abs(delta)<epsilon(1.0_dp).and.abs(q*9.0_dp)<epsilon(1.0_dp)) then
            y = -sign(1.0_DP,r)*(abs(r)+sqrt(r*r-q*q*q))**THIRD
            if (abs(y) < epsilon(1.0_DP)) then
               z = 0.0_DP
            else
               z = q / y
            end if
            x = y + z - aa*THIRD
         else
            y=(9.0_dp*cc - aa*bb)/(2.0_dp*q*9.0_dp)
            if(y>sqrt(3.0_dp).or.y<1.0_dp) then
               y=(4.0_dp*aa*bb - 9.0_dp*cc - aa**3)/(q*9.0_dp)
            end if
            x=y
         end if
! write(*,*) x
      end if



    end subroutine internal_cubic_roots

  end subroutine sign_nsh



  !============================================================================!
  ! To caclulate the sign function of a matrix, iteratively and without        !
  ! diagonalisation, a modified Newton's scheme can be employed much like the  !
  ! one in this module used for calculating the matrix square-root iteratively !
  ! for use in the Lowdin transformation. In particular, this routine makes    !
  ! use of the accelerated, stabilised update of [Chen and Chow].              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   A            (input)  : The matrix for which the sign function will be   !
  !                           computed.                                        !
  !   tol          (input)  : The tolerance on the 2-norm stopping criterion.  !
  !   X            (output) : The computed sign-function matrix                !
  !   matmuls      (output) : The number of matrix-matrix operations required  !
  !                           to apply the function.                           !
  !   status       (output) : Success = 0.
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, May 2014,                                        !
  !----------------------------------------------------------------------------!
  ! This routine scales as O(log(N)) in matrix-matrix multiplications.         !
  !============================================================================!
  subroutine sign_newton(A,tol,X,matmuls,status)
    use comms, only : pub_on_root
    use timer, only : timer_clock
    use utils, only : utils_alloc_check, utils_dealloc_check
    implicit none
    type(smearing_matrix),  intent(inout) :: A
    real(kind=dp), intent(in) :: tol
    type(smearing_matrix),  intent(inout) :: X ! Inout to preserve allocation
    integer,       intent(out) :: matmuls
    integer,       intent(out) :: status

    type(smearing_matrix) :: Xinv, Xdiff
    type(smearing_matrix) :: X0
    real(kind=dp) :: emin, emax, y, alpha, alphahat
    real(kind=dp) :: eigtol
    real(kind=dp) :: mnorm
    real(kind=dp), dimension(:), allocatable :: mnorm_history
    integer :: mnorm_history_length, mnorm_history_num, mnorm_history_maxlength
    integer :: i
    integer :: low_stat
    integer :: low_maxiters, low_iter
    integer :: ierr
    logical :: scale
    real(kind=dp) :: mu
    integer :: D
    real(kind=dp) :: coeffs(2)
    integer :: num_positive_grads, num_positive_to_stop
    real(kind=dp) :: tol_cgce, tol_scale

    call timer_clock('sign_newton',1)

    status=0

    matmuls=0

    alphahat = 1.697702485255767507654391_dp

    mnorm_history_maxlength=5
    mnorm_history_length=0
    mnorm_history_num=0
    num_positive_grads=0
    num_positive_to_stop=3

    tol_scale=1e-5_dp
    tol_cgce=1e-7_dp

    allocate(mnorm_history(mnorm_history_maxlength), stat=ierr)
    call utils_alloc_check('sign_newton','mnorm_history',ierr)
    mnorm_history=0.0_dp

    !    if(matrix_dimension(A)<10000.and..false.) then ! Always true until the matrix vector multiplication can be improved.
    eigtol=1e-5_dp

    !ja531-> sparsefoe:
    !       call matrix_allocate(X0,A)
    call matrix_allocate(X0,X)

    call matrix_copy(A,X0)
    !ja531-> sparsefoe:
    !       call matrix_allocate(X2,A)
    call matrix_allocate(Xinv,X)

    call matrix_allocate(Xdiff,X)
    !       y=emin/emax

    y=1e-6_dp


    call matrix_transpose(X0,Xinv)
    call matrix_scale(Xinv,1.0_dp/(matrix_induced_norm(X0)*matrix_induced_norm(Xinv)))
    D=8

    i=0
    do
       i=i+1

       call nsh_inverse(X0,Xinv,D,.true.,matmuls)

       if(scale) then
          mu=sqrt((matrix_norm(Xinv,2))/(matrix_norm(X0,2)))
       else
          mu=1.0_dp
       end if

       call matrix_copy(X0,Xdiff)
       call matrix_scale(Xdiff,-1.0_dp)

       call matrix_scale(X0,mu)
       call matrix_axpy(X0,Xinv,1.0_dp/mu)
       call matrix_scale(X0,0.5_dp)

       call matrix_axpy(Xdiff,X0,1.0_dp)

       mnorm=sqrt(matrix_trace(Xdiff,Xdiff))/sqrt(matrix_trace(X0,X0))

       if(pub_on_root) then
          if(debug_info_toggle) write(stdout,*) "smearing-> delta=",mnorm,"<== SMEAR"
       end if


       mnorm_history_length=mnorm_history_length+1
       if(mnorm_history_length<=mnorm_history_maxlength) then
          mnorm_history(mnorm_history_length)=mnorm
          mnorm_history_num=mnorm_history_length
          call poly_fit(1, mnorm_history, coeffs, mnorm_history_length)
       else
          mnorm_history(1:mnorm_history_maxlength-1)=mnorm_history(2:mnorm_history_maxlength)
          mnorm_history(mnorm_history_maxlength)=mnorm
          mnorm_history_num=mnorm_history_maxlength
          call poly_fit(1, mnorm_history, coeffs)
       end if

       if(coeffs(2)>0.0_dp) then
          num_positive_grads=num_positive_grads+1
       else
          num_positive_grads=0
       end if


       if(scale.and.mnorm<tol_scale) scale=.false.

       write(*,*) sqrt(matrix_trace(Xdiff,Xdiff)),sqrt(tol_cgce*sqrt(matrix_trace(X0,X0))/sqrt(matrix_trace(Xinv,Xinv))),&
            &sqrt(matrix_trace(X0,X0)),sqrt(matrix_trace(Xinv,Xinv)),mnorm,mnorm_history(max(mnorm_history_num-1,1)),scale


       if(mnorm_history_length>1) then
          if((sqrt(matrix_trace(Xdiff,Xdiff))<=sqrt(tol_cgce*sqrt(matrix_trace(X0,X0))/sqrt(matrix_trace(Xinv,Xinv)))).or.&
               & (.not.scale.and.mnorm>mnorm_history(mnorm_history_num-1)/2.0_dp)) then
             call matrix_copy(X0,X)
             if(pub_on_root) then
                if(debug_info_toggle) write(stdout,*) "smearing-> exit1:<== SMEAR"
             end if
             exit
          end if
       else
          if((sqrt(matrix_trace(Xdiff,Xdiff))<=sqrt(tol_cgce*sqrt(matrix_trace(X0,X0))/sqrt(matrix_trace(Xinv,Xinv))))) then
             call matrix_copy(X0,X)
             if(pub_on_root) then
                if(debug_info_toggle) write(stdout,*) "smearing-> exit2:<== SMEAR"
             end if

             exit
          end if
       end if


       if(pub_on_root) then
          if(mod(i,1)==0 .and. debug_info_toggle) write(stdout,*) "smearing-> ",mnorm, i, tol, "<== SMEAR"
       end if
    end do
    call matrix_free(X0)
    call matrix_free(Xinv)
    call matrix_free(Xdiff)

    if(pub_on_root) then
       do i=1,mnorm_history_maxlength
          write(stdout,*) "mnorm_history:",i,"->",mnorm_history(i)
       end do
    end if


    deallocate(mnorm_history, stat=ierr)
    call utils_dealloc_check('sign_newton','mnorm_history',ierr)


    call timer_clock('sign_newton',2)
  end subroutine sign_newton



  !============================================================================!
  ! This routine implements a Fermi Operator Expansion where the fractional    !
  ! occupancy eigenvalues are treated differently. The 1 and 0 eigenvalues are !
  ! treated with a purification method, the fractional eigenvalues get a       !
  ! Fermi-Dirac Chebyshev expansion. To project the regions out, a boxcar      !
  ! function is implemented as a contour expansion.                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   A            (input)  : The input matrix / Hamiltonian                   !
  !   S            (input)  : Metric / Overlap matrix                          !
  !   npoints      (input)  : number of points on the contour integration      !
  !   origin       (input)  : Origin of the boxcar                             !
  !   halfwidth    (input)  : The halfwidth of the boxcar                      !
  !   skew         (input)  : The skew of the boxcar                           !
  !   ignore_S     (input)  : Treat the Hamiltonian as if S=identity           !
  !   ignore_proj  (input)  : Don't return the projector matrix in proj        !
  !   no_solve     (input)  : Use inversion or other means rather than solves  !
  !   Ap           (output) : The result / density kernel                      !
  !   proj         (output) : The projector matrix used to make Ap             !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, sometime in 2019.                                !
  !============================================================================!
  subroutine boxcar_projection_contourint(A,S,npoints,origin,halfwidth,skew,&
       ignore_S,ignore_proj,Ap,proj,no_solve)
    use comms, only: pub_on_root, pub_my_proc_id
    use constants, only: cmplx_0, pi, stdout

    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none
    type(smearing_matrix), intent(inout) :: A
    type(smearing_matrix), intent(inout) :: S
    integer,               intent(in)    :: npoints
    real(kind=dp),         intent(in)    :: origin
    real(kind=dp),         intent(in)    :: halfwidth
    real(kind=dp),         intent(in)    :: skew
    logical,               intent(in)    :: ignore_S
    logical,               intent(in)    :: ignore_proj
    type(smearing_matrix), intent(inout) :: Ap
    type(smearing_matrix), intent(inout) :: proj
    logical,               intent(in)    :: no_solve

    integer :: N
    real(kind=dp) :: len
    real(kind=dp) :: alpha
    type(smearing_matrix) :: tmpA
    type(smearing_matrix) :: tmpB
    type(smearing_matrix) :: tmpC
    type(smearing_matrix) :: tmpAp
    type(smearing_matrix) :: tmpAp_real
    type(smearing_matrix) :: tmpproj_real
    type(smearing_matrix) :: tmpproj
    type(smearing_matrix) :: tmpS
    integer :: i
    real(kind=dp), dimension(:), allocatable :: theta
    complex(kind=dp), dimension(:), allocatable :: z
    complex(kind=dp), dimension(:), allocatable :: omega
    integer :: ierr

    N=npoints
    len=halfwidth
    alpha=skew

    if(pub_on_root) then
       if(debug_info_toggle) then
          write(stdout,*) "entering boxcar_projection_contourint"
       end if
    end if

    allocate(theta(N),stat=ierr)
    call utils_alloc_check('boxcar_projection_contourint','theta',ierr)
    allocate(z(N),stat=ierr)
    call utils_alloc_check('boxcar_projection_contourint','z',ierr)
    allocate(omega(N),stat=ierr)
    call utils_alloc_check('boxcar_projection_contourint','omega',ierr)


    call matrix_allocate(tmpA,A,iscmplx=.true.)
    call matrix_copy(A,tmpA)

    call matrix_allocate(tmpB,A,iscmplx=.true.)

    call matrix_scale_cmplx(tmpB,cmplx_0)

    call matrix_allocate(tmpC,A,iscmplx=.true.)
    call matrix_scale_cmplx(tmpC,cmplx_0)

    call matrix_allocate(tmpAp,Ap,iscmplx=.true.)
    call matrix_scale_cmplx(tmpAp,cmplx_0)

    if(.not.ignore_S) then
       call matrix_allocate(tmpAp_real,Ap)
       call matrix_scale(tmpAp_real,0.0_dp)
    end if

    if(.not.ignore_proj) then
       call matrix_allocate(tmpS,S,iscmplx=.true.)
       call matrix_copy(S,tmpS)

       call matrix_allocate(tmpproj,proj,iscmplx=.true.)
       call matrix_scale_cmplx(tmpproj,cmplx_0)
    end if

    call matrix_scale(Ap,0.0_dp)

    if(pub_on_root) then
       if(debug_info_toggle) then
          write(stdout,*) "beginning of contour loop. N=",N
          write(stdout,*) "PROC:",pub_my_proc_id, Ap%standard_is_cmplx, S%standard_is_cmplx, tmpAp_real%standard_is_cmplx
       end if
    end if



    if(mod(N,2)==0) then
       do i=1,(N/2)


          theta(i)=(2.0_dp*pi/real(N,dp))*(real(i,dp)-0.5_dp)
          z(i)=cmplx(origin,0.0_dp,dp) + cmplx(len,0.0_dp,dp)*cmplx(cos(theta(i)),alpha*sin(theta(i)),dp)
          omega(i)=cmplx(len,0.0_dp,dp)*cmplx(alpha*cos(theta(i)),sin(theta(i)))

          if(pub_on_root) then
             if(debug_info_toggle) then
                write(stdout,*) "iteration=",i, theta(i),z(i),omega(i)
             end if
          end if


          call matrix_copy(tmpA,tmpB)
          if(ignore_S) then
             call matrix_scale_cmplx(tmpB,cmplx(-1.0_dp,0.0_dp,dp),z(i))
          else

             call matrix_scale_cmplx(tmpB,cmplx(-1.0_dp,0.0_dp,dp))

             call matrix_axpy_cmplx(tmpB,S,z(i))
          end if


          if(no_solve) then
             call matrix_invert(tmpB)
             call matrix_multiply(1.0_dp,tmpB,tmpA,0.0_dp,tmpC)
          else
             call matrix_solve(tmpA,tmpB,tmpC)
          end if

          call matrix_axpy_cmplx(tmpAp,tmpC,omega(i))

          if(.not.ignore_proj) then
             call matrix_copy(tmpA,tmpB)
             if(ignore_S) then
                call matrix_scale_cmplx(tmpB,cmplx(-1.0_dp,0.0_dp,dp),z(i))
             else
                call matrix_scale_cmplx(tmpB,cmplx(-1.0_dp,0.0_dp,dp))
                call matrix_axpy_cmplx(tmpB,S,z(i))
             end if

             if(no_solve) then
                call matrix_invert(tmpB)
                call matrix_multiply(1.0_dp,tmpB,tmpS,0.0_dp,tmpC)
             else
                call matrix_solve(tmpS,tmpB,tmpC)
             end if

             call matrix_axpy_cmplx(tmpproj,tmpC,omega(i))
          end if

       end do

       if(ignore_S) then
          call matrix_copy(tmpAp,Ap)
          call matrix_scale(Ap,1.0_dp/real(N/2,dp))
       else
          call matrix_copy(tmpAp,tmpAp_real)
          if(pub_on_root) then
             if(debug_info_toggle) then
                write(stdout,*) "PROC:",pub_my_proc_id, Ap%standard_is_cmplx, S%standard_is_cmplx, tmpAp_real%standard_is_cmplx
             end if
          end if

          call matrix_multiply(1.0_dp/real(N/2,dp),S,tmpAp_real,0.0_dp,Ap)

          call matrix_free(tmpAp_real)


          if(.not.ignore_proj) then
             call matrix_allocate(tmpproj_real,proj)

             call matrix_copy(tmpproj,tmpproj_real)

             call matrix_scale(tmpproj_real,1.0_dp/real(N/2,dp))

!             call matrix_multiply(1.0_dp/real(N/2,dp),S,tmpproj_real,0.0_dp,proj)

             call matrix_copy(tmpproj_real,proj)

             call matrix_free(tmpproj_real)
          end if

       end if

    else
       do i=1,N
          theta(i)=(2.0_dp*pi/real(N,dp))*(real(i,dp)-0.5_dp)
          z(i)=cmplx(origin,0.0_dp,dp) + cmplx(len,0.0_dp,dp)*cmplx(cos(theta(i)),alpha*sin(theta(i)),dp)
          omega(i)=cmplx(len,0.0_dp,dp)*cmplx(alpha*cos(theta(i)),sin(theta(i)))

          call matrix_copy(tmpA,tmpB)
          if(ignore_S) then
             call matrix_scale_cmplx(tmpB,cmplx(-1.0_dp,0.0_dp,dp),z(i))
          else
             call matrix_scale_cmplx(tmpB,cmplx(-1.0_dp,0.0_dp,dp))
             call matrix_axpy_cmplx(tmpB,S,z(i))
          end if

          if(no_solve) then
             call matrix_invert(tmpB)
             call matrix_multiply(1.0_dp,tmpB,tmpA,0.0_dp,tmpC)
          else
             call matrix_solve(tmpA,tmpB,tmpC)
          end if

          call matrix_axpy_cmplx(tmpAp,tmpC,omega(i))


          if(.not.ignore_proj) then
             call matrix_copy(tmpA,tmpB)
             if(ignore_S) then
                call matrix_scale_cmplx(tmpB,cmplx(-1.0_dp,0.0_dp,dp),z(i))
             else
                call matrix_scale_cmplx(tmpB,cmplx(-1.0_dp,0.0_dp,dp))
                call matrix_axpy_cmplx(tmpB,S,z(i))
             end if

             if(no_solve) then
                call matrix_invert(tmpB)
                call matrix_multiply(1.0_dp,tmpB,tmpS,0.0_dp,tmpC)
             else
                call matrix_solve(tmpS,tmpB,tmpC)
             end if

             call matrix_axpy_cmplx(tmpproj,tmpC,omega(i))
          end if


       end do

       if(ignore_S) then
          call matrix_copy(tmpAp,Ap)
          call matrix_scale(Ap,1.0_dp/real(N,dp))
       else
          call matrix_copy(tmpAp,tmpAp_real)
          if(pub_on_root) then
             if(debug_info_toggle) then
                write(stdout,*) "PROC:",pub_my_proc_id, Ap%standard_is_cmplx, S%standard_is_cmplx, tmpAp_real%standard_is_cmplx
             end if
          end if

          call matrix_multiply(1.0_dp/real(N,dp),S,tmpAp_real,0.0_dp,Ap)
          call matrix_free(tmpAp_real)

          if(.not.ignore_proj) then
             call matrix_allocate(tmpproj_real,proj)

             call matrix_copy(tmpproj,tmpproj_real)

             call matrix_scale(tmpproj_real,1.0_dp/real(N,dp))

!             call matrix_multiply(1.0_dp/real(N,dp),S,tmpproj_real,0.0_dp,proj)

             call matrix_copy(tmpproj_real,proj)

             call matrix_free(tmpproj_real)
          end if

       end if

    end if


    if(.not.ignore_proj) then
       call matrix_free(tmpS)
       call matrix_free(tmpproj)
    end if

    call matrix_free(tmpA)
    call matrix_free(tmpB)
    call matrix_free(tmpC)
    call matrix_free(tmpAp)



    deallocate(theta,stat=ierr)
    call utils_dealloc_check('boxcar_projection_contourint','theta',ierr)
    deallocate(z,stat=ierr)
    call utils_dealloc_check('boxcar_projection_contourint','z',ierr)
    deallocate(omega,stat=ierr)
    call utils_dealloc_check('boxcar_projection_contourint','omega',ierr)


  end subroutine boxcar_projection_contourint


  !============================================================================!
  ! This routine calculates a Chebyshev expansion of a matrix.                 !
  ! It needs N matrix multiplications, where N is the length of the expansion. !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The input matrix.                                !
  !   chebycoeffs  (input)  : The coefficients of the expansion.               !
  !   Ap           (output) : The result / density kernel                      !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, May 2020.                                        !
  !============================================================================!
  subroutine slow_cheby_resum(H,cheby_coeffs,rho)
    implicit none
    type(smearing_matrix), intent(inout)  :: H
    real(kind=dp), dimension(:), intent(in) :: cheby_coeffs
    type(smearing_matrix), intent(inout) :: rho

    type(smearing_matrix) :: T0, T1, T
    integer :: N, i

    N=size(cheby_coeffs)

    call matrix_scale(rho,0.0_dp)

    if(N>0) then
       call matrix_allocate(T0,rho)
       call matrix_scale(T0,0.0_dp,1.0_dp)
       call matrix_axpy(rho,T0,cheby_coeffs(1))
    end if

    if(N>1) then
       call matrix_allocate(T1,rho)
       call matrix_copy(H,T1)
       call matrix_axpy(rho,T1,cheby_coeffs(2))
    end if

    if(N>2) then
       call matrix_allocate(T,rho)
       do i=3,N
          call matrix_copy(T0,T)
          call matrix_multiply(2.0_dp,H,T1,-1.0_dp,T)
          call matrix_axpy(rho,T,cheby_coeffs(i))
          call matrix_copy(T1,T0)
          call matrix_copy(T,T1)
       end do
       call matrix_free(T)
    end if

    if(N>1) then
       call matrix_free(T1)
    end if
    if(N>0) then
       call matrix_free(T0)
    end if

  end subroutine slow_cheby_resum

  !============================================================================!
  ! This routine resums a Chebyshev expansion using the algorithm of           !
  ! Liang and Head Gordon:                                                     !
  ! It needs 2*sqrt(N) matrix multiplications, where N is the length of the    !
  ! expansion.                                                                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The input matrix.                                !
  !   chebycoeffs  (input)  : The coefficients of the expansion.               !
  !   rho          (output) : The result / density kernel                      !
  !   matmuls      (output) : The number of matmuls required.                  !
  !   forcelayers  (input)  : Force a certain number of layers in the method.  !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Mar 2017,                                        !
  !============================================================================!
  subroutine fast_cheby_resum(H,cheby_coeffs,rho,matmuls,force_layers)


    use timer, only : timer_clock
    use utils, only : utils_alloc_check, utils_dealloc_check
    implicit none
    type(smearing_matrix), intent(inout)  :: H
    real(kind=dp), dimension(:), intent(in) :: cheby_coeffs
    type(smearing_matrix), intent(inout) :: rho
    integer, intent(out) :: matmuls
    integer, optional, intent(in) :: force_layers

    integer :: N
    integer :: ml
    real(kind=dp), dimension(:), allocatable :: coeffs
    real(kind=dp), dimension(:,:), allocatable :: dccoeffs
    integer :: i
    integer :: layers
    integer :: jh,jl, jlold, jhold
    integer :: j,k
    integer, dimension(:), allocatable :: num, parent, start
    type(smearing_matrix), dimension(:), allocatable :: exes
    type(smearing_matrix), dimension(:), allocatable :: accumat
    !    integer :: l
    integer :: matmuls_divconq
    type(smearing_matrix) :: T0, T1, T
    integer :: nmats

    integer :: ierr

    call timer_clock('fast_cheby_resum',1)

    N=matrix_dimension(H)

    ml=size(cheby_coeffs)-1

    layers=max(floor(log10(real(ml,dp))/log10(2.0_dp))-1,1)
    layers=min(layers,4)
    if(present(force_layers)) then
       layers=force_layers
!       write(*,*) "we have ",layers, " layers"
    end if
    allocate(coeffs(1:(ml+1)*(layers+1)),stat=ierr)
    call utils_alloc_check('fast_cheby_resum','coeffs',ierr)
    coeffs(1:ml+1)=cheby_coeffs
    ! set coeffs here

    allocate(dccoeffs(1,size(cheby_coeffs)),stat=ierr)
    call utils_alloc_check('fast_cheby_resum','dccoeffs',ierr)

!if(pub_on_root) then
!    write(*,*) "ml=",ml
!    write(*,*) "we have ",layers, " layers"
!    write(*,*) "num should be: ", 2**(layers+1)-1
!end if

    allocate(num(1:2**(layers+1)-1),stat=ierr)
    call utils_alloc_check('fast_cheby_resum','num',ierr)
    allocate(parent(1:2**(layers+1)-1),stat=ierr)
    call utils_alloc_check('fast_cheby_resum','parent',ierr)
    allocate(start(1:2**(layers+1)-1),stat=ierr)
    call utils_alloc_check('fast_cheby_resum','start',ierr)

    ! Construct binary tree
    jh=0
    jl=0
    jlold=1
    jhold=1

    num=0
    num(1)=(ml+1)*2

    parent=0

    start=0
    start(1)=1
    do i=0,layers
       jl=jl+2**(max(0,i-1))
       jh=jh+2**(i)
       do j=jl,jh
          parent(j)=jlold+(floor(real((j-jl+2),dp)/2.0_dp)-1)
          if(mod(j,2)==0) then
             num(j)=floor(real(num(parent(j)),dp)/2.0_dp)
          else
             num(j)=floor(real((num(parent(j))-1),dp)/2.0_dp)
          end if
       end do
       jlold=jl
       jhold=jh
    end do

    do i=2,size(num)
       start(i)=start(i-1)+num(i-1)+1
    end do

    !    coeffs=0
    !    coeffs(1:ml+1)=c
    jl=0
    jh=0
    jlold=1
    jhold=ml+1
    do i=2,size(num)
       jl=jhold+1
       jh=jl+num(i)
       if(mod(i,2)==0) then
          do j=jh,jl,-1
             k=j-jl
             coeffs(j)=coeffs(2*k+start(parent(i)))
          end do
       else
          coeffs(jh)=2*coeffs(2*num(i)+1+start(parent(i)))
          do j=jh-1,jl+1,-1
             k=j-jl
             coeffs(j)=2*coeffs(2*k+1+start(parent(i)))-coeffs(k+1+jl)
          end do
          coeffs(jl)=coeffs(start(parent(i))+1)-(coeffs(jl+1)/2)
       end if
       jlold=jl
       jhold=jh

    end do
    !dccoeffs=0.0_dp
    dccoeffs(1,1:1+size(coeffs)-start(size(num)-2**layers+1)) = coeffs(start(size(num)-2**layers+1):size(coeffs))

    deallocate(coeffs,stat=ierr)
    call utils_dealloc_check('fast_cheby_resum','coeffs',ierr)
    deallocate(parent,stat=ierr)
    call utils_dealloc_check('fast_cheby_resum','parent',ierr)

    allocate(exes(1:layers+1),stat=ierr)
    call utils_alloc_check('fast_cheby_resum','exes',ierr)
    do i = 1,layers+1
       !ja531->sparsefoe:
!       call matrix_allocate(exes(i),H)
       call matrix_allocate(exes(i),rho)
       call matrix_scale(exes(i),0.0_dp)
    end do
    matmuls_divconq=0
    call matrix_copy(H,exes(1))

    do i=2,layers+1
       matmuls_divconq=matmuls_divconq+1
       call matrix_multiply(1.0_dp,exes(i-1),exes(i-1),0.0_dp,exes(i))
       call matrix_scale(exes(i),2.0_dp,-1.0_dp)
    end do

    allocate(accumat(1:2**layers),stat=ierr)
    call utils_alloc_check('fast_cheby_resum','accumat',ierr)
    do i=1,2**layers
       !ja531->sparsefoe:
!       call matrix_allocate(accumat(i),H)
       call matrix_allocate(accumat(i),rho)
       call matrix_scale(accumat(i),0.0_dp)
    end do

    !ja531->sparsefoe:
!    call matrix_allocate(T0,H)
    call matrix_allocate(T0,rho)
!    call matrix_allocate(T1,H)
    call matrix_allocate(T1,rho)
!    call matrix_allocate(T,H)
    call matrix_allocate(T,rho)

    do i=0,maxval(num(size(num)-2**layers+1:size(num)))
       if(i==0) then
          call matrix_scale(T0,0.0_dp)
          call matrix_scale(T0,0.0_dp,1.0_dp)
          call matrix_copy(T0,T)
       else if(i==1) then
          call matrix_copy(exes(layers+1),T1)
          call matrix_copy(T1,T)
       else
          call matrix_copy(T0,T)
          call matrix_multiply(2.0_dp,exes(layers+1),T1,-1.0_dp,T)
          matmuls_divconq=matmuls_divconq+1
          call matrix_copy(T1,T0)
          call matrix_copy(T,T1)
       end if
       !       do j=1,(2**layers)
       !ja531-> ODDFIX: EXPERIMENTAL!
       ! old version
!       do j=(2**(layers-1))+1,(2**layers)
!          if(i<=num(size(num)-2**layers+j)) then
!             !             do li=lbar,P
!             call matrix_axpy(accumat(j),T,dccoeffs(1,start(size(num)-2**layers+j)-start(size(num)-2**layers+1)+1+i))
!             !             end do
!             !             call matrix_axpy(NLcheb_acc(j),T,dccoeffs(1,start(size(num)-2**layers+j)&
!             !                  -start(size(num)-2**layers+1)+1+i))
!
!          end if
!       end do
       ! new version
       do j=1,(2**layers)
          if(i<=num(size(num)-2**layers+j)) then
             !             do li=lbar,P
             call matrix_axpy(accumat(j),T,dccoeffs(1,start(size(num)-2**layers+j)-start(size(num)-2**layers+1)+1+i))
             !             end do
             !             call matrix_axpy(NLcheb_acc(j),T,dccoeffs(1,start(size(num)-2**layers+j)&
             !                  -start(size(num)-2**layers+1)+1+i))

          end if
       end do

    end do

    deallocate(dccoeffs,stat=ierr)
    call utils_dealloc_check('fast_cheby_resum','dccoeffs',ierr)
    call matrix_free(T0)
    call matrix_free(T1)
    call matrix_free(T)

    k=0
!ja531-> ODDFIX: EXPERIMENTAL!
    ! old version
!    do i=1,layers-1
!       nmats=2**i
!
!       do j=(2**(layers-i-1))+1,2**(layers-i)
!          k=k+1
!          matmuls_divconq=matmuls_divconq+1
!
!          call matrix_multiply(1.0_dp,exes(layers-i+1),accumat(nmats*(j-1)+1+(2**(i-1))),1.0_dp,accumat(nmats*(j-1)+1))
!          matmuls_divconq=matmuls_divconq+1
!       end do
!    end do
    ! new version
    do i=1,layers-1
       nmats=2**i

       do j=1,2**(layers-i)
          k=k+1
          matmuls_divconq=matmuls_divconq+1

          call matrix_multiply(1.0_dp,exes(layers-i+1),accumat(nmats*(j-1)+1+(2**(i-1))),1.0_dp,accumat(nmats*(j-1)+1))
          matmuls_divconq=matmuls_divconq+1
       end do
    end do

    call matrix_multiply(1.0_dp,exes(1),accumat(1+(2**(layers-1))),1.0_dp,accumat(1))
    matmuls_divconq=matmuls_divconq+1
    call matrix_copy(accumat(1),rho)

    do i = 1,layers+1
       call matrix_free(exes(i))
    end do
    deallocate(exes,stat=ierr)
    call utils_dealloc_check('fast_cheby_resum','exes',ierr)

    do i=1,2**layers
       call matrix_free(accumat(i))
    end do
    deallocate(accumat,stat=ierr)
    call utils_dealloc_check('fast_cheby_resum','accumat',ierr)

    deallocate(num,stat=ierr)
    call utils_dealloc_check('fast_cheby_resum','num',ierr)
    deallocate(start,stat=ierr)
    call utils_dealloc_check('fast_cheby_resum','start',ierr)

    matmuls=matmuls_divconq

    call timer_clock('fast_cheby_resum',2)

  end subroutine fast_cheby_resum



  !============================================================================!
  ! This routine searches for the electron number conserving chemical          !
  ! potential. It either uses hyperbolic trigonometry to modify the chemical   !
  ! potential of the input density kernel, or recalculates the density kernel  !
  ! at the new chemical potential, depending on the condition number of the    !
  ! matrix.
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The Hamiltonian matrix.                          !
  !   S            (input)  : The overlap matrix.                              !
  !   CCH          (input)  : The contra-covariant Hamiltonian (invS.H)        !
  !   rho          (inout)  : The density kernel                               !
  !   orthogonal   (input)  : Ignore S (treat it as identity)                  !
  !   use_nsh      (input)  : Use Hotelling method in change-mu routine        !
  !   Ne           (input)  : Number of electrons                              !
  !   fermi        (inout)  : chemical potential                               !
  !   smearing     (input)  : The electron temperature                         !
  !   tol          (input)  : Tolerance of error in number of electrons        !
  !   lowb         (input)  : Low limit of chemical potential search window    !
  !   highb        (input)  : High limit of chemical potential search window   !
  !   foe_type     (input)  : The FOE method to use, if recalculating          !
  !   threshold    (input)  : Magnitude of change in chemical potential        !
  !                           that will result in recalculating density kernel !
  !   proj_matmuls (output) : number of matmuls used to project                !
  !   totalmatmuls (output) : number of matmuls used to do everything          !
  !   avoid_solve  (input)  : Avoid using the solve routine & use inverse      !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons... in 2017?                                       !
  !============================================================================!
  subroutine findfermi_safe_newton(H,S,invS,CCH,rho,orthogonal,use_nsh,Ne,fermi,&
       smearing,tol,lowb,highb,foe_type, threshold,proj_matmuls,total_matmuls,&
       avoid_solve)
    use timer, only : timer_clock
    use comms, only: pub_on_root



    implicit none
    type(smearing_matrix),  intent(inout) :: H
    type(smearing_matrix),  intent(inout) :: S
    type(smearing_matrix),  intent(inout) :: invS
    type(smearing_matrix),  intent(inout) :: CCH
    type(smearing_matrix),  intent(inout) :: rho
    logical, intent(in) :: orthogonal
    logical, intent(in) :: use_nsh
    real(kind=dp), intent(in)    :: Ne  ! I can't imagine a time when this would be non-integer... but just in case.
    real(kind=dp), intent(inout) :: fermi
    real(kind=dp), intent(in)    :: smearing
    real(kind=dp), intent(in)    :: tol
!    type(smearing_matrix),  intent(inout) :: orth_ham
    real(kind=dp), intent(in)    :: lowb
    real(kind=dp), intent(in)    :: highb
    character(len=16), intent(in) :: foe_type
    real(kind=dp), intent(inout) :: threshold
    integer,       intent(inout) :: proj_matmuls
    integer,       intent(out)   :: total_matmuls
    logical,       intent(in)    :: avoid_solve

    integer :: maxit
    integer :: N
    integer :: matmuls=0
    type(smearing_matrix) :: T
    type(smearing_matrix) :: Tinit
    type(smearing_matrix) :: invdenom
    type(smearing_matrix) :: CC_T_trans, S_trans, T_trans
    real(kind=dp) :: Dmu
    real(kind=dp) :: fermi_tmp
    real(kind=dp) :: fermi_new
    logical :: have_invdenom
    real(kind=dp) :: fl, fh
    real(kind=dp) :: xl, xh
    real(kind=dp) :: swap
    real(kind=dp) :: dxold, dx
    real(kind=dp) :: df
    real(kind=dp) :: f
    integer :: j,i

    real(kind=dp) :: rnum

    real(kind=dp) :: ttrace

    !    real(kind=dp) :: threshold
    character     :: method_used

    integer :: outunit

    call timer_clock('find_fermi_safe_newton',1)

    maxit=100
    N=matrix_dimension(rho)


    call matrix_allocate(T,rho)
    call matrix_allocate(Tinit,rho)
    call matrix_allocate(invdenom,rho)

    if(orthogonal) then
       call matrix_copy(rho,Tinit)
    else
       ! Set up Tinit to be contra-covariant
       call matrix_multiply(1.0_dp,rho,S,0.0_dp,Tinit)
    end if

    ttrace=matrix_trace(Tinit)
!    if(pub_on_root) then
!       write(stdout,*) "FFSN START: Ktrace=",ttrace
!    end if

    call matrix_scale(Tinit,2.0_dp,-1.0_dp)

    total_matmuls=0
    Dmu=lowb-fermi
    fermi_tmp=fermi
    fermi_new=0.0_dp
    have_invdenom=.false.
    call matrix_copy(Tinit,T)
    call changemu(H,S,invS,CCH,T,orthogonal,use_nsh,fermi_tmp,Dmu,smearing,threshold,&
         foe_type,avoid_solve,invdenom,have_invdenom,matmuls,method_used)



!           if(pub_on_root) then
!              write(stdout,*) "FIRST"
!              write(stdout,*) "fermi=",fermi_tmp
!              write(stdout,*) "Dmu=",Dmu
!              write(stdout,*) "smearing=",smearing
!              write(stdout,*) "threshold=",threshold
!              outunit=utils_unit()
!              open(unit=outunit,file="first_H.mtx")
!           end if
!           call sparse_show_matrix(H%dataSPAM3,outunit,matlab_format=.true.)
!           if(pub_on_root) then
!              close(outunit)
!           end if
!           if(pub_on_root) then
!              outunit=utils_unit()
!              open(unit=outunit,file="first_S.mtx")
!           end if
!           call sparse_show_matrix(S%dataSPAM3,outunit,matlab_format=.true.)
!           if(pub_on_root) then
!              close(outunit)
!           end if
!           if(pub_on_root) then
!              outunit=utils_unit()
!              open(unit=outunit,file="first_T.mtx")
!           end if
!           call sparse_show_matrix(T%dataSPAM3,outunit,matlab_format=.true.)
!           if(pub_on_root) then
!              close(outunit)
!           end if
!           if(pub_on_root) then
!              outunit=utils_unit()
!              open(unit=outunit,file="first_CCH.mtx")
!           end if
!           call sparse_show_matrix(CCH%dataSPAM3,outunit,matlab_format=.true.)
!           if(pub_on_root) then
!              close(outunit)
!           end if
!
    i=0
    ! This loop is to catch bad chemical potentials... most likely this meant that a sign function fell
    ! directly on an eigenvalue.
    do
       i=i+1
       !       f=0.5_dp*(matrix_trace(T)+real(N,dp))-Ne
       fl=0.5_dp*(matrix_trace(T)+real(N,dp))-Ne
       ttrace=matrix_trace(T)
       if(pub_on_root) then
          if(debug_info_toggle) write(stdout,*) "smearing-> ", i,matmuls, "1=+=> ", &
               0.5_dp*(ttrace+real(N,dp)), Ne,Dmu, fermi_tmp, fermi_new, "<== SMEAR"
       end if
       if(fl/=fl.or.abs(0.5_dp*(matrix_trace(T)))>matrix_dimension(T)) then

          !          Dmu=fermi_new-fermi_tmp-(fermi_new*0.001_dp)
          call random_number(rnum)
          rnum=(1.0_dp-rnum)*2.0_dp
          ! We probably hit an eigenvalue with the shift... pole.

          if(abs(fermi_new)>epsilon(fermi_new)) then
             rnum=rnum*fermi_new
          else
             rnum=rnum*(10.0_dp**(i-8))
          end if

          !          if(T%matrix_type/=matrix_type_standard) then
          !             call comms_bcast(pub_root_proc_id,rnum)
          !          end if

          Dmu=fermi_new-fermi_tmp+(rnum)

          fermi_new=fermi_tmp+Dmu
          call changemu(H,S,invS,CCH,T,orthogonal,use_nsh,fermi_tmp,Dmu,smearing,threshold,&
               &foe_type,avoid_solve,invdenom,have_invdenom,matmuls,method_used)
          !          call changemu(T,fermi_tmp,Dmu,smearing,threshold,orth_ham, &
          !               invdenom,have_invdenom,matmuls,method_used)
          total_matmuls=total_matmuls+matmuls

       else
          exit
       end if
    end do
    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> ---2", matmuls,"<== SMEAR"
    end if

    ! Undo the shift
    if(orthogonal) then
       call matrix_scale(H,1.0_dp,+Dmu/smearing)
    else
       call matrix_axpy(H,S,+Dmu/smearing)
       call matrix_scale(CCH,1.0_dp,+Dmu/smearing)
    end if

    call matrix_copy(Tinit,T)

    Dmu=highb-fermi
    fermi_tmp=fermi
!    call changemu(T,fermi_tmp,Dmu,smearing,threshold,orth_ham,invdenom,have_invdenom,matmuls,method_used)
    call changemu(H,S,invS,CCH,T,orthogonal,use_nsh,fermi_tmp,Dmu,smearing,threshold,&
         foe_type,avoid_solve,invdenom,have_invdenom,matmuls,method_used)



!           if(pub_on_root) then
!              write(stdout,*) "SECOND"
!              write(stdout,*) "fermi=",fermi_tmp
!              write(stdout,*) "Dmu=",Dmu
!              write(stdout,*) "smearing=",smearing
!              write(stdout,*) "threshold=",threshold
!              outunit=utils_unit()
!              open(unit=outunit,file="second_H.mtx")
!           end if
!           call sparse_show_matrix(H%dataSPAM3,outunit,matlab_format=.true.)
!           if(pub_on_root) then
!              close(outunit)
!           end if
!           if(pub_on_root) then
!              outunit=utils_unit()
!              open(unit=outunit,file="second_S.mtx")
!           end if
!           call sparse_show_matrix(S%dataSPAM3,outunit,matlab_format=.true.)
!           if(pub_on_root) then
!              close(outunit)
!           end if
!           if(pub_on_root) then
!              outunit=utils_unit()
!              open(unit=outunit,file="second_T.mtx")
!           end if
!           call sparse_show_matrix(T%dataSPAM3,outunit,matlab_format=.true.)
!           if(pub_on_root) then
!              close(outunit)
!           end if
!           if(pub_on_root) then
!              outunit=utils_unit()
!              open(unit=outunit,file="second_CCH.mtx")
!           end if
!           call sparse_show_matrix(CCH%dataSPAM3,outunit,matlab_format=.true.)
!           if(pub_on_root) then
!              close(outunit)
!           end if
!

    i=0
    ! This loop is to catch bad chemical potentials... most likely this meant that a sign function fell
    ! directly on an eigenvalue.
    do
       i=i+1
       !       f=0.5_dp*(matrix_trace(T)+real(N,dp))-Ne
       fh=0.5_dp*(matrix_trace(T)+real(N,dp))-Ne
       ttrace=matrix_trace(T)
       if(pub_on_root) then
          if(debug_info_toggle) write(stdout,*) "smearing-> ",i, matmuls,"2=+=> ", &
               0.5_dp*(ttrace+real(N,dp)), Ne,Dmu, fermi_tmp, fermi_new, "<== SMEAR"
       end if
       if(fh/=fh.or.abs(0.5_dp*(matrix_trace(T)))>matrix_dimension(T)) then
          !          Dmu=fermi_new-fermi_tmp-(fermi_new*0.001_dp)

          call random_number(rnum)
          rnum=(1.0_dp-rnum)*2.0_dp
          ! We probably hit an eigenvalue with the shift... pole.

          if(abs(fermi_new)>epsilon(fermi_new)) then
             rnum=rnum*fermi_new
          else
             rnum=rnum*(10.0_dp**(i-8))
          end if

          !          if(T%matrix_type/=matrix_type_standard) then
          !             call comms_bcast(pub_root_proc_id,rnum)
          !          end if

          Dmu=fermi_new-fermi_tmp+(rnum)

          fermi_new=fermi_tmp+Dmu
!          call changemu(T,fermi_tmp,Dmu,smearing,threshold,orth_ham,invdenom, &
          !               have_invdenom,matmuls,method_used)
          call changemu(H,S,invS,CCH,T,orthogonal,use_nsh,fermi_tmp,Dmu,smearing,threshold,&
               &foe_type, avoid_solve,invdenom,have_invdenom,matmuls,method_used)
          total_matmuls=total_matmuls+matmuls
       else
          exit
       end if

    end do




    ttrace=matrix_trace(T)
    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> ---4", matmuls, 0.5_dp*(ttrace+real(N,dp)), "<== SMEAR"
    end if
    !fermi_new=fermi_tmp+Dmu

    if(fl<0.0_dp) then
       xl=lowb
       xh=highb
    else
       xh=lowb
       xl=highb
       swap=fl
       fl=fh
       fh=swap
    end if


    ! Check if the Chemical potential bounds were wide enough
    ! Don't think this is needed...
!     swap=xh
!     do i=1,10
!        if(fh<Ne) then
!
!           Dmu=(highb-lowb)/2.0_dp
!
!
!           if(pub_on_root) then
!              if(debug_info_toggle) write(stdout,*) "smearing-> bounds not wide enough!",fh , "<== SMEAR"
!           end if
!
!           call changemu(H,S,CCH,T,orthogonal,use_nsh,fermi_tmp,Dmu,smearing,threshold,&
!                &invdenom,have_invdenom,matmuls,method_used)
!           total_matmuls=total_matmuls+matmuls
!
!           fh=0.5_dp*(matrix_trace(T)+real(N,dp))-Ne
!           ttrace=matrix_trace(T)
!           if(pub_on_root) then
!              if(debug_info_toggle) write(stdout,*) "smearing-> ",i, matmuls,"69=+=> ", &
!                   0.5_dp*(ttrace+real(N,dp)), Ne,Dmu, fermi_tmp, fermi_new, "<== SMEAR"
!           end if
!           swap=swap+Dmu
!        else
!
!           xh=swap
!
!        end if
!     end do


    !    fermi_tmp=highb
    !    fermi_new=highb!+Dmu
    dxold=abs(xh-xl)
    dx=dxold

!    call matrix_scale(H,1.0_dp,+Dmu/smearing)
    if(orthogonal) then
       call matrix_scale(H,1.0_dp,+Dmu/smearing)
    else
       call matrix_axpy(H,S,+Dmu/smearing)
       call matrix_scale(CCH,1.0_dp,+Dmu/smearing)
    end if


    fermi_tmp=fermi_new
    Dmu=-Dmu
    fermi_new=fermi

    call matrix_copy(Tinit,T)
    f=0.5_dp*(matrix_trace(T)+real(N,dp))-Ne

    if(f < 0.0_dp) then
       xl=fermi_new
       fl=f
    else
       xh=fermi_new
       fh=f
    end if
    df=deriv_chemicalpot(T,smearing)
    total_matmuls=total_matmuls+1

    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> DERIV:", df,xl,fl,fh,fermi_new, "<== SMEAR"
    end if

    !    call dense_write(exes(layers+1)%dataDEM,"exes.mat")

    do j=1,maxit
       !    call dense_write(T%dataDEM,"T")

       if(((fermi_new-xh)*df-f)*((fermi_new-xl)*df-f) >= 0.0_dp .or. abs(2.0_dp*f) > abs(dxold*df) .or. df/=df) then
          dxold=dx
          dx=0.5_dp*(xh-xl)
          fermi_tmp=fermi_new
          fermi_new=xl+dx
          if(xl == fermi_new) then
             exit
          end if
       else
          dxold=dx
          dx=f/df
          fermi_tmp=fermi_new
          fermi_new=fermi_new-dx
          if(fermi_tmp == fermi_new) then
             if(pub_on_root) then
                if(debug_info_toggle) write(stdout,*) "smearing-> B",fermi_tmp,fermi_new, Dmu, "<== SMEAR"
             end if
             if(j==1) then
                Dmu=0.0_dp
                call matrix_copy(Tinit,T)
             end if

             exit
          end if
       end if
       if(abs(dx) < tol) then
          exit
       end if
       Dmu=fermi_new-fermi_tmp

!       call changemu(T,fermi_tmp,Dmu,smearing,threshold,orth_ham,invdenom,have_invdenom,matmuls,method_used)
       call changemu(H,S,invS,CCH,T,orthogonal,use_nsh,fermi_tmp,Dmu,smearing,threshold,&
            foe_type,avoid_solve,invdenom,have_invdenom,matmuls,method_used)
       if(method_used=="P") then
          proj_matmuls=matmuls
       else if(method_used=="H".and.matmuls>((proj_matmuls*110)/100)) then
          !          threshold=(Dmu/smearing)-threshold
          threshold=-((3.0_dp*Dmu/(4.0_dp*smearing)) - 0.5*threshold)
       end if

       total_matmuls=total_matmuls+matmuls
       i=0
       ! This loop is to catch bad chemical potentials... most likely this meant that a sign function fell
       ! directly on an eigenvalue.
       do
          i=i+1
          f=0.5_dp*(matrix_trace(T)+real(N,dp))-Ne
          ttrace=matrix_trace(T)
          if(pub_on_root) then
             if(debug_info_toggle) write(stdout,*) "smearing-> ",i, matmuls, "3=+=> ", &
                  0.5_dp*(ttrace+real(N,dp)), Ne,Dmu, fermi_tmp, fermi_new, "<== SMEAR"
          end if
          if(f/=f.or.abs(f)>Ne*100.0_dp) then
             fermi_new=fermi_new-(fermi_new*0.001_dp)
             Dmu=fermi_new-fermi_tmp
!             call changemu(T,fermi_tmp,Dmu,smearing,threshold,orth_ham,invdenom,have_invdenom,matmuls,method_used)
             call changemu(H,S,invS,CCH,T,orthogonal,use_nsh,fermi_tmp,Dmu,smearing,threshold,&
                  foe_type,avoid_solve,invdenom,have_invdenom,matmuls,method_used)
             if(method_used=="P") then
                proj_matmuls=matmuls
             else if(method_used=="H".and.matmuls>((proj_matmuls*110)/100)) then
                !    threshold=(Dmu/smearing)-threshold
                threshold=-((3.0_dp*Dmu/(4.0_dp*smearing)) - 0.5*threshold)
             end if
             total_matmuls=total_matmuls+matmuls
          else
             exit
          end if
       end do

       df=deriv_chemicalpot(T,smearing)
       total_matmuls=total_matmuls+1
       if(f < 0.0_dp) then
          xl=fermi_new
          fl=f
       else
          xh=fermi_new
          fh=f
       end if

    end do

    call matrix_scale(T,0.5_dp,0.5_dp)



    if(orthogonal) then

       call matrix_copy(T,rho)

    else ! Need to make rho contravariant

       if(.not.avoid_solve) then

          call matrix_allocate(CC_T_trans,rho)

          ! transpose CC_outside_T
          call matrix_transpose(T,CC_T_trans)

          !       call matrix_free(T)

          ! transpose S
          call matrix_allocate(S_trans,S)
          call matrix_transpose(S,S_trans)

          ! solve S'outside_T' = CC_outside_T'
          call matrix_allocate(T_trans,rho)
          !       if(pub_on_root) then
          !          write(stdout,*) "Entering FindFermi Contra output Solve"
          !       end if
          call matrix_solve(CC_T_trans,S_trans,T_trans)

          call matrix_free(S_trans)
          call matrix_free(CC_T_trans)

          ! transpose outside_T
          call matrix_transpose(T_trans,rho)

          call matrix_free(T_trans)

       else

          ! rho=CCT*invS
          call matrix_multiply(1.0_dp,T,invS,0.0_dp,rho)
          total_matmuls=total_matmuls+1



       end if


    end if


    fermi=fermi_new
    call matrix_free(Tinit)
    call matrix_free(T)
    call matrix_free(invdenom)

    call timer_clock('find_fermi_safe_newton',2)

  end subroutine findfermi_safe_newton


  !============================================================================!
  ! Same as above, but this version corrects the chemical potential at each    !
  ! annealing step rather than just at the end. (experimental)                 !
  ! This routine searches for the electron number conserving chemical          !
  ! potential. It either uses hyperbolic trigonometry to modify the chemical   !
  ! potential of the input density kernel, or recalculates the density kernel  !
  ! at the new chemical potential, depending on the condition number of the    !
  ! matrix.                                                                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The Hamiltonian matrix.                          !
  !   S            (input)  : The overlap matrix.                              !
  !   CCH          (input)  : The contra-covariant Hamiltonian (invS.H)        !
  !   rho          (inout)  : The density kernel                               !
  !   orthogonal   (input)  : Ignore S (treat it as identity)                  !
  !   use_nsh      (input)  : Use Hotelling method in change-mu routine        !
  !   Ne           (input)  : Number of electrons                              !
  !   fermi        (inout)  : chemical potential                               !
  !   smearing     (input)  : The electron temperature                         !
  !   tol          (input)  : Tolerance of error in number of electrons        !
  !   lowb         (input)  : Low limit of chemical potential search window    !
  !   highb        (input)  : High limit of chemical potential search window   !
  !   foe_type     (input)  : The FOE method to use, if recalculating          !
  !   threshold    (input)  : Magnitude of change in chemical potential        !
  !                           that will result in recalculating density kernel !
  !   proj_matmuls (output) : number of matmuls used to project                !
  !   totalmatmuls (output) : number of matmuls used to do everything          !
  !   avoid_solve  (input)  : Avoid using the solve routine & use inverse      !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons... in 2017?                                       !
  !============================================================================!
  subroutine findfermi_safe_newton_anneal(H,S,invS,CCH,rho,orthogonal,use_nsh,Ne,fermi,smearing,tol,lowb,highb,&
       &foe_type,threshold,proj_matmuls,total_matmuls,avoid_solve)
    use timer, only : timer_clock
    use comms, only: pub_on_root


    use utils, only: utils_dealloc_check
    implicit none
    type(smearing_matrix),  intent(inout) :: H
    type(smearing_matrix),  intent(inout) :: S
    type(smearing_matrix),  intent(inout) :: invS
    type(smearing_matrix),  intent(inout) :: CCH
    type(smearing_matrix),  intent(inout) :: rho
    logical, intent(in) :: orthogonal
    logical, intent(in) :: use_nsh
    real(kind=dp), intent(in)    :: Ne  ! I can't imagine a time when this would be non-integer... but just in case.
    real(kind=dp), intent(inout) :: fermi
    real(kind=dp), intent(in)    :: smearing
    real(kind=dp), intent(in)    :: tol
!    type(smearing_matrix),  intent(inout) :: orth_ham
    real(kind=dp), intent(in)    :: lowb
    real(kind=dp), intent(in)    :: highb
    character(len=16), intent(in) :: foe_type
    real(kind=dp), intent(inout) :: threshold
    integer,       intent(inout) :: proj_matmuls
    integer,       intent(out)   :: total_matmuls
    logical,       intent(in)    :: avoid_solve

    integer :: maxit
    integer :: N
    integer :: matmuls=0
    type(smearing_matrix) :: T, Thigh, Tlow
    type(smearing_matrix) :: CCH_loc, H_loc
    type(smearing_matrix) :: CCH2, CCH_bkp, H2
    type(smearing_matrix) :: Tinit
    type(smearing_matrix) :: invdenom
    type(smearing_matrix) :: CC_T_trans, S_trans, T_trans
    real(kind=dp) :: Dmu
    real(kind=dp) :: fermi_tmp
    real(kind=dp) :: fermi_new
    logical :: have_invdenom
    real(kind=dp) :: fl, fh
    real(kind=dp) :: xl, xh
    real(kind=dp) :: swap
    real(kind=dp) :: dxold, dx
    real(kind=dp) :: df
    real(kind=dp) :: f
    integer :: j,i

    real(kind=dp) :: rnum

    real(kind=dp) :: ttrace

    real(kind=dp) :: cutoff, smearing_loc, specwidth, cheby_threshold, tol_loc
    integer :: ierr
    real(kind=dp), dimension(:), allocatable :: coeffs
    integer :: ilevel, nlevels

    !    real(kind=dp) :: threshold
    character     :: method_used

    integer :: outunit

    call timer_clock('find_fermi_safe_newton_anneal',1)

    maxit=100
    N=matrix_dimension(rho)



    call matrix_allocate(T,rho)
    call matrix_allocate(Tinit,rho)
    call matrix_allocate(Tlow,rho)
    call matrix_allocate(Thigh,rho)
    call matrix_allocate(invdenom,rho)

    ! Set max error of Chebyshev expansion of tanh
    cheby_threshold=1e-9_dp

    ! cutoff is the maximum argument of tanh that is accurate
    cutoff=15.0_dp

    call matrix_allocate(CCH_loc,CCH)
    call matrix_copy(CCH,CCH_loc)

    call matrix_allocate(H_loc,H)
    call matrix_copy(H,H_loc)

    ! convert from argument of FD to argument of tanh: H -> -H/2
    call matrix_scale(CCH_loc,-0.5_dp)

    ! work out spectral width of Hamiltonian = specwidth
    call matrix_allocate(CCH_bkp,CCH_loc)
    call matrix_copy(CCH_loc,CCH_bkp)
    call matrix_allocate(CCH2,CCH_loc,CCH_bkp)
    call matrix_multiply(1.0_dp,CCH_loc,CCH_bkp,0.0_dp,CCH2)
    call matrix_free(CCH_bkp)
    call matrix_allocate(H2,H,H)
    call matrix_multiply(1.0_dp,S,CCH2,0.0_dp,H2)
    call matrix_free(CCH2)
    call extremum_eigs(H2,S,tol=0.00001_dp,maxeval=specwidth)
    call matrix_free(H2)
    specwidth=sqrt(specwidth)
    if(pub_on_root) then
       write(stdout,*) "FF Spectral width:",specwidth
    end if

    ! work out number of annealing levels = log2(specwidth/cutoff)
    ! max to avoid dividing by 0!
    nlevels = max(ceiling(log(specwidth/cutoff)/log(2.0_dp)),1)+1
    if(pub_on_root) then
       write(stdout,*) "FF Annealing Nlevels=",nlevels
    end if

    ! H' -> H/2^nlevels
    call matrix_scale(CCH_loc,1.0_dp/real(2**nlevels,dp))
    call matrix_scale(H_loc,1.0_dp/real(2**nlevels,dp))

    ttrace=matrix_trace(CCH_loc)
    if(pub_on_root) then
       write(stdout,*) "FF Annealing CCH_loc: Ktrace=",ttrace
    end if

    smearing_loc=smearing*real(2**nlevels,dp)


    call fermi_operator_annealing(H_loc,S,CCH_loc,Tinit,orthogonal,.false.,matmuls)
    total_matmuls=total_matmuls+matmuls


    ttrace=matrix_trace(Tinit)
    if(pub_on_root) then
       write(stdout,*) "FFSN START: Ktrace=",ttrace
    end if

    call matrix_scale(Tinit,2.0_dp,-1.0_dp)

    total_matmuls=0
    Dmu=lowb-fermi
    fermi_tmp=fermi
    fermi_new=0.0_dp
    have_invdenom=.false.
    call matrix_copy(Tinit,Tlow)
    call changemu(H_loc,S,invS,CCH_loc,Tlow,orthogonal,use_nsh,fermi_tmp,Dmu,smearing_loc,threshold,&
         foe_type, avoid_solve,invdenom, have_invdenom,matmuls,method_used)
    total_matmuls=total_matmuls+matmuls

!
    i=0
    ! This loop is to catch bad chemical potentials... most likely this meant that a sign function fell
    ! directly on an eigenvalue.
    do
       i=i+1
       !       f=0.5_dp*(matrix_trace(T)+real(N,dp))-Ne
       fl=0.5_dp*(matrix_trace(Tlow)+real(N,dp))-Ne
       ttrace=matrix_trace(Tlow)
       if(pub_on_root) then
          if(debug_info_toggle) write(stdout,*) "smearing-> ", i,matmuls, "1=+=> ", &
               0.5_dp*(ttrace+real(N,dp)), Ne,Dmu, fermi_tmp, fermi_new, "<== SMEAR"
       end if
       if(fl/=fl.or.abs(0.5_dp*(matrix_trace(Tlow)))>matrix_dimension(Tlow)) then

          !          Dmu=fermi_new-fermi_tmp-(fermi_new*0.001_dp)
          call random_number(rnum)
          rnum=(1.0_dp-rnum)*2.0_dp
          ! We probably hit an eigenvalue with the shift... pole.

          if(abs(fermi_new)>epsilon(fermi_new)) then
             rnum=rnum*fermi_new
          else
             rnum=rnum*(10.0_dp**(i-8))
          end if

          !          if(T%matrix_type/=matrix_type_standard) then
          !             call comms_bcast(pub_root_proc_id,rnum)
          !          end if

          Dmu=fermi_new-fermi_tmp+(rnum)

          fermi_new=fermi_tmp+Dmu
          call changemu(H_loc,S,invS,CCH_loc,Tlow,orthogonal,use_nsh,fermi_tmp,Dmu,smearing_loc,threshold,&
               &foe_type,avoid_solve,invdenom,have_invdenom,matmuls,method_used)
          !          call changemu(T,fermi_tmp,Dmu,smearing,threshold,orth_ham, &
          !               invdenom,have_invdenom,matmuls,method_used)
          total_matmuls=total_matmuls+matmuls

       else
          exit
       end if
    end do
    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> ---2", matmuls,"<== SMEAR"
    end if

    ! Undo the shift
    if(orthogonal) then
       call matrix_scale(H_loc,1.0_dp,+Dmu/smearing_loc)
    else
       call matrix_axpy(H_loc,S,+Dmu/smearing_loc)
       call matrix_scale(CCH_loc,1.0_dp,+Dmu/smearing_loc)
    end if

    call matrix_copy(Tinit,Thigh)

    Dmu=highb-fermi
    fermi_tmp=fermi
!    call changemu(T,fermi_tmp,Dmu,smearing,threshold,orth_ham,invdenom,have_invdenom,matmuls,method_used)
    call changemu(H_loc,S,invS,CCH_loc,Thigh,orthogonal,use_nsh,fermi_tmp,Dmu,smearing_loc, &
         & threshold,foe_type,avoid_solve,invdenom,have_invdenom,matmuls,method_used)



    i=0
    ! This loop is to catch bad chemical potentials... most likely this meant that a sign function fell
    ! directly on an eigenvalue.
    do
       i=i+1
       !       f=0.5_dp*(matrix_trace(T)+real(N,dp))-Ne
       fh=0.5_dp*(matrix_trace(Thigh)+real(N,dp))-Ne
       ttrace=matrix_trace(Thigh)
       if(pub_on_root) then
          if(debug_info_toggle) write(stdout,*) "smearing-> ",i, matmuls,"2=+=> ", &
               0.5_dp*(ttrace+real(N,dp)), Ne,Dmu, fermi_tmp, fermi_new, "<== SMEAR"
       end if
       if(fh/=fh.or.abs(0.5_dp*(matrix_trace(Thigh)))>matrix_dimension(Thigh)) then
          !          Dmu=fermi_new-fermi_tmp-(fermi_new*0.001_dp)

          call random_number(rnum)
          rnum=(1.0_dp-rnum)*2.0_dp
          ! We probably hit an eigenvalue with the shift... pole.

          if(abs(fermi_new)>epsilon(fermi_new)) then
             rnum=rnum*fermi_new
          else
             rnum=rnum*(10.0_dp**(i-8))
          end if

          !          if(T%matrix_type/=matrix_type_standard) then
          !             call comms_bcast(pub_root_proc_id,rnum)
          !          end if

          Dmu=fermi_new-fermi_tmp+(rnum)

          fermi_new=fermi_tmp+Dmu
!          call changemu(T,fermi_tmp,Dmu,smearing,threshold,orth_ham,invdenom, &
          !               have_invdenom,matmuls,method_used)
          call changemu(H_loc,S,invS,CCH_loc,Thigh,orthogonal,use_nsh,fermi_tmp,Dmu,smearing_loc,threshold,&
               &foe_type,avoid_solve,invdenom,have_invdenom,matmuls,method_used)
          total_matmuls=total_matmuls+matmuls
       else
          exit
       end if

    end do




    ttrace=matrix_trace(Thigh)
    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> ---4", matmuls, 0.5_dp*(ttrace+real(N,dp)), "<== SMEAR"
    end if
    !fermi_new=fermi_tmp+Dmu

    if(fl<0.0_dp) then
       xl=lowb
       xh=highb
    else
       xh=lowb
       xl=highb
       swap=fl
       fl=fh
       fh=swap
    end if



    !    fermi_tmp=highb
    !    fermi_new=highb!+Dmu
    dxold=abs(xh-xl)
    dx=dxold

!    call matrix_scale(H,1.0_dp,+Dmu/smearing)
    if(orthogonal) then
       call matrix_scale(H_loc,1.0_dp,+Dmu/smearing_loc)
    else
       call matrix_axpy(H_loc,S,+Dmu/smearing_loc)
       call matrix_scale(CCH_loc,1.0_dp,+Dmu/smearing_loc)
    end if


    fermi_tmp=fermi_new
    Dmu=-Dmu
    fermi_new=fermi

    call matrix_copy(Tinit,T)
    f=0.5_dp*(matrix_trace(T)+real(N,dp))-Ne

    if(f < 0.0_dp) then
       xl=fermi_new
       fl=f
    else
       xh=fermi_new
       fh=f
    end if
    df=deriv_chemicalpot(T,smearing_loc)
    total_matmuls=total_matmuls+1

    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> DERIV:", df,xl,fl,fh,fermi_new, "<== SMEAR"
    end if

    !    call dense_write(exes(layers+1)%dataDEM,"exes.mat")


    call get_hypmultid_chebyshev_coeffs(-1.0_dp,1.0_dp,cheby_threshold,coeffs)

    tol_loc=tol*(10.0_dp**real(nlevels,dp))

    do ilevel=1,nlevels+1

       if(ilevel>1) then



          smearing_loc=smearing_loc/2.0_dp

          if(orthogonal) then
             call matrix_scale(H_loc,2.0_dp)
          else
             call matrix_scale(H_loc,2.0_dp)
             call matrix_scale(CCH_loc,2.0_dp)
          end if

!          call fermi_operator_annealing(H_loc,S,CCH_loc,T,orthogonal,.false.,matmuls)
!          total_matmuls=total_matmuls+matmuls

          call fast_cheby_resum(T,coeffs,Tinit,matmuls)
          call matrix_copy(Tinit,T)
          ttrace=matrix_trace(T)
          if(pub_on_root) then
             write(stdout,*) "FF Quenching ",ilevel," CCT: Ttrace=",ttrace
          end if
          total_matmuls=total_matmuls+matmuls

          f=0.5_dp*(matrix_trace(T)+real(N,dp))-Ne

          Dmu=-10.0_dp*tol_loc
          call matrix_copy(T,Tlow)
          call changemu(H_loc,S,invS,CCH_loc,Tlow,orthogonal,use_nsh,fermi_tmp,Dmu,smearing_loc,threshold,&
               &foe_type,avoid_solve,invdenom,have_invdenom,matmuls,method_used)
!          call fast_cheby_resum(Tlow,coeffs,Tinit,matmuls)
!          call matrix_copy(Tinit,Tlow)
          ttrace=matrix_trace(Tlow)
          if(pub_on_root) then
             write(stdout,*) "FF Quenching ",ilevel," CCT: Tlowtrace=",ttrace
          end if
          total_matmuls=total_matmuls+matmuls

          xl=fermi_tmp+Dmu
          fl=0.5_dp*(matrix_trace(Tlow)+real(N,dp))-Ne

          if(orthogonal) then
             call matrix_scale(H_loc,1.0_dp,+Dmu/smearing_loc)
          else
             call matrix_axpy(H_loc,S,+Dmu/smearing_loc)
             call matrix_scale(CCH_loc,1.0_dp,+Dmu/smearing_loc)
          end if

          Dmu=20.0_dp*tol_loc
          call matrix_copy(T,Thigh)
          call changemu(H_loc,S,invS,CCH_loc,Thigh,orthogonal,use_nsh,fermi_tmp,Dmu,smearing_loc,threshold,&
               &foe_type,avoid_solve,invdenom,have_invdenom,matmuls,method_used)
!          call fast_cheby_resum(Thigh,coeffs,Tinit,matmuls)
!          call matrix_copy(Tinit,Thigh)
          ttrace=matrix_trace(Thigh)
          if(pub_on_root) then
             write(stdout,*) "FF Quenching ",ilevel," CCT: Thightrace=",ttrace
          end if
          total_matmuls=total_matmuls+matmuls
          !
          xh=fermi_tmp+Dmu
          fh=0.5_dp*(matrix_trace(Thigh)+real(N,dp))-Ne

          if(orthogonal) then
             call matrix_scale(H_loc,1.0_dp,+Dmu/smearing_loc)
          else
             call matrix_axpy(H_loc,S,+Dmu/smearing_loc)
             call matrix_scale(CCH_loc,1.0_dp,+Dmu/smearing_loc)
          end if


!          xl=lowb
!          xh=highb

          dxold=abs(xh-xl)
          dx=dxold

          if(f < 0.0_dp) then
             xl=fermi_new
             fl=f
          else
             xh=fermi_new
             fh=f
          end if

          df=deriv_chemicalpot(T,smearing_loc)
          total_matmuls=total_matmuls+1

          if(pub_on_root) then
             if(debug_info_toggle) write(stdout,*) "smearing-> DERIV:", df,xl,fl,fh,fermi_new, "<== SMEAR"
          end if

          call matrix_copy(T,Tinit)

          tol_loc=tol_loc/10.0_dp

       end if

       do j=1,maxit
          !    call dense_write(T%dataDEM,"T")
          if(pub_on_root) then
             if(debug_info_toggle) write(stdout,*) "smearing-> INFOZ:",(((fermi_new-xh)*df-f)*((fermi_new-xl)*df-f)), &
                  & abs(2.0_dp*f), abs(dxold*df), "dxold=",dxold,"df=", df,"fermi_new=", fermi_new, "xl=",xl,"xh=",xh,&
                  & "dx=",dx,"f=",f, "<== SMEAR"
          end if

          if(((fermi_new-xh)*df-f)*((fermi_new-xl)*df-f) >= 0.0_dp .or. abs(2.0_dp*f) > abs(dxold*df) .or. df/=df) then
             dxold=dx
             dx=0.5_dp*(xh-xl)
             fermi_tmp=fermi_new
             fermi_new=xl+dx
             if(xl == fermi_new) then
                exit
             end if
          else
             dxold=dx
             dx=f/df
             fermi_tmp=fermi_new
             fermi_new=fermi_new-dx
             if(fermi_tmp == fermi_new) then
                if(pub_on_root) then
                   if(debug_info_toggle) write(stdout,*) "smearing-> B",fermi_tmp,fermi_new, Dmu, "<== SMEAR"
                end if
                if(j==1) then
                   Dmu=0.0_dp
                   call matrix_copy(Tinit,T)
                end if

                exit
             end if
          end if
          if(abs(dx) < tol_loc) then
             exit
          end if
          Dmu=fermi_new-fermi_tmp

          if(pub_on_root) then
             if(debug_info_toggle) write(stdout,*) "smearing-> Fermi level:",fermi_tmp,Dmu, "<== SMEAR"
          end if

          !       call changemu(T,fermi_tmp,Dmu,smearing,threshold,orth_ham,invdenom,have_invdenom,matmuls,method_used)
          call changemu(H_loc,S,invS,CCH_loc,T,orthogonal,use_nsh,fermi_tmp,Dmu,smearing_loc,threshold,&
               &foe_type, avoid_solve, invdenom,have_invdenom,matmuls,method_used)
          if(method_used=="P") then
             proj_matmuls=matmuls
          else if(method_used=="H".and.matmuls>((proj_matmuls*110)/100)) then
             !          threshold=(Dmu/smearing)-threshold
             threshold=-((3.0_dp*Dmu/(4.0_dp*smearing_loc)) - 0.5*threshold)
          end if

          total_matmuls=total_matmuls+matmuls
          i=0
          ! This loop is to catch bad chemical potentials... most likely this meant that a sign function fell
          ! directly on an eigenvalue.
          do
             i=i+1
             f=0.5_dp*(matrix_trace(T)+real(N,dp))-Ne
             ttrace=matrix_trace(T)
             if(pub_on_root) then
                if(debug_info_toggle) write(stdout,*) "smearing-> ",i, matmuls, "3=+=> ", &
                     0.5_dp*(ttrace+real(N,dp)), Ne,Dmu, fermi_tmp, fermi_new, "<== SMEAR"
             end if
             if(f/=f.or.abs(f)>Ne*100.0_dp) then
                fermi_new=fermi_new-(fermi_new*0.001_dp)
                Dmu=fermi_new-fermi_tmp
                !             call changemu(T,fermi_tmp,Dmu,smearing,threshold,orth_ham,invdenom,have_invdenom,matmuls,method_used)
                call changemu(H_loc,S,invS,CCH_loc,T,orthogonal,use_nsh,fermi_tmp,Dmu,smearing_loc,&
                     &threshold,foe_type,avoid_solve,invdenom,have_invdenom,matmuls,method_used)
                if(method_used=="P") then
                   proj_matmuls=matmuls
                else if(method_used=="H".and.matmuls>((proj_matmuls*110)/100)) then
                   !    threshold=(Dmu/smearing)-threshold
                   threshold=-((3.0_dp*Dmu/(4.0_dp*smearing_loc)) - 0.5*threshold)
                end if
                total_matmuls=total_matmuls+matmuls
             else
                exit
             end if
          end do

          df=deriv_chemicalpot(T,smearing_loc)
          total_matmuls=total_matmuls+1
          if(f < 0.0_dp) then
             xl=fermi_new
             fl=f
          else
             xh=fermi_new
             fh=f
          end if

       end do

    end do

    call matrix_scale(T,0.5_dp,0.5_dp)

    if(orthogonal) then

       call matrix_copy(T,rho)

    else ! Need to make rho contravariant

       if(.not.avoid_solve) then

          call matrix_allocate(CC_T_trans,rho)

          ! transpose CC_outside_T
          call matrix_transpose(T,CC_T_trans)

          !       call matrix_free(T)

          ! transpose S
          call matrix_allocate(S_trans,S)
          call matrix_transpose(S,S_trans)

          ! solve S'outside_T' = CC_outside_T'
          call matrix_allocate(T_trans,rho)
          if(pub_on_root) then
             write(stdout,*) "Entering FindFermi Contra output Solve"
          end if
          call matrix_solve(CC_T_trans,S_trans,T_trans)

          call matrix_free(S_trans)
          call matrix_free(CC_T_trans)

          ! transpose outside_T
          call matrix_transpose(T_trans,rho)

          call matrix_free(T_trans)

       else

          ! rho=CCT*invS
          call matrix_multiply(1.0_dp,T,invS,0.0_dp,rho)
          total_matmuls=total_matmuls+1

       end if

    end if

    deallocate(coeffs,stat=ierr)
    call utils_dealloc_check('fermi_operator_projection_contourint','cheby_coeffs',ierr)

    fermi=fermi_new
    call matrix_free(CCH_loc)
    call matrix_free(H_loc)

    call matrix_free(Tinit)
    call matrix_free(T)
    call matrix_free(Tlow)
    call matrix_free(Thigh)
    call matrix_free(invdenom)

    call timer_clock('find_fermi_safe_newton_anneal',2)

  end subroutine findfermi_safe_newton_anneal




  !============================================================================!
  ! This routine gives the derivative of the trace of the density kernel, with !
  ! respect to chemical potential. It uses the scaled and shifted tanh matrix  !
  ! to do it.                                                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   T            (input)  : The tanh matrix.                                 !
  !   smearing     (input)  : The electron temperature.                        !
  !   df           (output) : The change in trace(tanh)                        !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons... in 2017?                                       !
  !============================================================================!
  function deriv_chemicalpot(T,smearing) result(df)
    use timer, only : timer_clock
    implicit none
    type(smearing_matrix),  intent(inout) :: T
    real(kind=dp), intent(in)    :: smearing
    type(smearing_matrix) :: T2, tmp
    real(kind=dp) :: df
    call timer_clock('smearing_deriv_chemicalpot',1)

    !ja531->sparsefoe: assuming that tanh sparsity (K) is good enough for product.
    call matrix_allocate(T2,T)
    call matrix_allocate(tmp,T)
    call matrix_copy(T,tmp) ! Do we really need a copy?
    call matrix_multiply(1.0_dp,T,tmp,0.0_dp,T2)
    df=-(1.0_dp/4.0_dp)*(1.0_dp-matrix_trace(T2))
    !    df=-(1.0_dp/(4.0_dp*smearing))*(1.0_dp-matrix_trace(T2))
    call matrix_free(T2)
    call matrix_free(tmp)

    call timer_clock('smearing_deriv_chemicalpot',2)
  end function deriv_chemicalpot


  !============================================================================!
  ! This routine gives the derivative of trace(HK) with respect to H. This     !
  ! routine has never been tried or tested. It is not used in this module (at  !
  ! the moment).                                                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   T            (input)  : The tanh matrix.                                 !
  !   K            (input)  : The density kernel.                              !
  !   H            (input)  : The Hamiltonian.                                 !
  !   smearing     (input)  : The electron temperature.                        !
  !   deriv        (output) : The change in trace(KH)                          !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons... in 2017?                                       !
  !============================================================================!
  function deriv_E_H(T,K,H,smearing) result(deriv)
    use timer, only : timer_clock
    implicit none
    type(smearing_matrix),  intent(inout) :: T
    type(smearing_matrix),  intent(inout) :: K
    type(smearing_matrix),  intent(inout) :: H
    real(kind=dp), intent(in)    :: smearing
    type(smearing_matrix) :: T2, tmp
    type(smearing_matrix) :: deriv
    call timer_clock('smearing_deriv_E_H',1)

    call matrix_allocate(T2,T)
    call matrix_allocate(tmp,T)
    call matrix_copy(T,tmp) ! Do we really need a copy?
    call matrix_multiply(1.0_dp,T,tmp,0.0_dp,T2)

    call matrix_scale(T2,1.0_dp/(4.0_dp*smearing),-1.0_dp/(4.0_dp*smearing))

    call matrix_copy(K,deriv)
    call matrix_multiply(1.0_dp,H,T2,1.0_dp,deriv)

    call matrix_free(T2)
    call matrix_free(tmp)

    call timer_clock('smearing_deriv_E_H',2)
  end function deriv_E_H


  !============================================================================!
  ! This routine modifies the implicit chemical potential in a density kernel  !
  ! if possible, or creates a new density kernel at the new chemical potential !
  ! It either uses hyperbolic trigonometry to modify the chemical              !
  ! potential of the input density kernel, or recalculates the density kernel  !
  ! at the new chemical potential, depending on the condition number of the    !
  ! matrix.                                                                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The Hamiltonian matrix.                          !
  !   S            (input)  : The overlap matrix.                              !
  !   invS         (input)  : The inverse overlap matrix.                      !
  !   CCH          (input)  : The contra-covariant Hamiltonian (invS.H)        !
  !   T            (inout)  : Scaled & shifted density kernel (tanh matrix)    !
  !   orthogonal   (input)  : Ignore S (treat it as identity)                  !
  !   use_nsh      (input)  : Use Hotelling method in this routine             !
  !   mu           (inout)  : chemical potential                               !
  !   Dmu          (input)  : Requested change in chemical potential           !
  !   smearing     (input)  : The electron temperature                         !
  !   threshold    (input)  : Magnitude of change in chemical potential        !
  !                           that will result in recalculating density kernel !
  !   foe_type     (input)  : The FOE method to use, if recalculating          !
  !   avoid_solve  (input)  : Avoid using the solve routine & use inverse      !
  !   invinit      (input)  : Initial inverse matrix                           !
  !   have_invinit (input)  : Whether we have invinit or not                   !
  !   matmuls      (output) : number of matmuls used to do everything          !
  !   method_used  (output) : Whether the denskern was modified or rebuilt     !
  !   force_reinit (input)  : Always start inverse from scratch                !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons... in 2017?                                       !
  !============================================================================!
  ! Assuming that if orthogonal==.false., then T is contra-covariant
  subroutine changemu(H,S,invS,CCH,T,orthogonal,use_nsh,mu,Dmu,smearing,threshold,&
       foe_type,avoid_solve,invinit,have_invinit,matmuls,method_used,force_reinit)
    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast
    use rundat, only: pub_FOE_cheby_thres
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_assert, utils_dealloc_check
    implicit none
    type(smearing_matrix),  intent(inout) :: H
    type(smearing_matrix),  intent(inout) :: S
    type(smearing_matrix),  intent(inout) :: invS
    type(smearing_matrix),  intent(inout) :: CCH
    type(smearing_matrix),  intent(inout) :: T
    logical,       intent(in)    :: orthogonal
    logical,       intent(in)    :: use_nsh
    real(kind=dp), intent(in)    :: mu
    real(kind=dp), intent(in)    :: Dmu
    real(kind=dp), intent(in)    :: smearing
    real(kind=dp), intent(in)    :: threshold
    character(len=16), intent(in) :: foe_type
    logical,       intent(in)    :: avoid_solve
    type(smearing_matrix),  intent(inout) :: invinit
    logical,       intent(inout) :: have_invinit
    integer,       intent(out)   :: matmuls
    character,     intent(out)   :: method_used
    logical,       intent(in), optional :: force_reinit


    real(kind=dp), dimension(:), allocatable :: coeffs
    real(kind=dp) :: cheby_threshold
    integer :: ierr

    integer :: N
    type(smearing_matrix) :: temp1
    type(smearing_matrix) :: temp2
    type(smearing_matrix) :: denom
    type(smearing_matrix) :: CC_T_trans, denom_trans, T_trans
    real(kind=dp) :: eigtol
    real(kind=dp) :: min_abs_eig, max_abs_eig
    real(kind=dp) :: condi
    integer :: matmuls_reinit=0
    integer :: matmuls_inv=0
    real(kind=dp) :: newmu
    integer :: D
    logical :: reinit
    real(kind=dp) :: norm, condest, ttrace
    logical :: condestnan

    logical :: do_proj

    call timer_clock('smearing_changemu',1)

    reinit=.false.
    if(present(force_reinit)) then
       reinit=force_reinit
    end if



    matmuls=0

    N=matrix_dimension(T)

    eigtol=1e-6_dp
    cheby_threshold=1e-8_dp

    cheby_threshold=pub_FOE_cheby_thres

    call matrix_allocate(denom,T)

    ! Here we use a locally preconditioned-CG algorithm on the square of a matrix
    ! to find the infinum and supremum absolute eigenvalues of the matrix.
    ! There must be a better way than this?
    ! Perhaps the better way is just not to bother and to use matmul counts from
    ! previous iterations to determine a threshold...
    call matrix_scale(denom,0.0_dp,1.0_dp)
    call matrix_axpy(denom,T,tanh(Dmu/(2.0_dp*smearing)))
    call matrix_allocate(temp1,denom)


    matmuls_reinit=0
    matmuls_inv=0
    newmu=Dmu + mu

    norm=matrix_induced_norm(T)
    ttrace=matrix_trace(T)

    condestnan=.false.
    if(abs(1.0_dp-abs(tanh(Dmu/(2.0_dp*smearing))))<epsilon(condest)) then
       condestnan=.true.
       condest=huge(condest)
    else
       condest=((1.0_dp+abs(tanh(Dmu/(2.0_dp*smearing))))/(1.0_dp-abs(tanh(Dmu/(2.0_dp*smearing)))))
    end if

    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> changemu-> do_proj->", &
            condest,abs(Dmu/(2.0_dp*smearing)),norm
    end if


    do_proj=(condest>1e4_dp.or.condest<1.0_dp.or.condest/=condest.or.abs(Dmu/(2.0_dp*smearing))>threshold.or. &
         & norm/=norm.or.ttrace/=ttrace.or.abs(0.5_dp*(ttrace+real(N,dp)))>real(N,dp))


    call comms_bcast(pub_root_proc_id,do_proj)

    if(do_proj) then
       method_used="P"!rojection
       if(pub_on_root.and.debug_info_toggle) write(stdout,*) &
            "smearing-> Recomputing Fermi operator using projection method. Delta-mu : ",Dmu,Dmu/(2.0_dp*smearing), "<== SMEAR"
       call matrix_scale(temp1,0.0_dp)


       select case(FOE_TYPE)
       case(SIGNUM_PROJECTION)
          call utils_assert(orthogonal,"Error in Smearing_changemu: Sign proj method only supports orthogonal Hams.")
          call matrix_scale(H,1.0_dp,-Dmu/smearing)
          call fermi_operator_projection(H,temp1,matmuls)
       case(CONTOUR_PROJECTIO)
          if(orthogonal) then
             call matrix_scale(H,1.0_dp,-Dmu/smearing)
             call fermi_operator_projection_contourint(H,H,H,H,temp1,orthogonal,&
                  contravariant_kernel=.true.,avoid_solve=avoid_solve, matmuls_out=matmuls)
          else
             call matrix_axpy(H,S,-Dmu/smearing)
             call matrix_scale(CCH,1.0_dp,-Dmu/smearing)
             call fermi_operator_projection_contourint(H,S,invS,CCH,temp1,orthogonal,&
                  contravariant_kernel=.false.,avoid_solve=avoid_solve, matmuls_out=matmuls)
          end if
       case(CERIOTTI_KUHNE_PA)
          call utils_assert(orthogonal,"Error in Smearing_changemu: Ceriotti method only supports orthogonal Hams.")
          call fermi_operator_chebyshev(H,temp1,matmuls)

       case(ANNEALING)
!          call utils_assert(.not.orthogonal,"Error in Smearing_changemu: Annealing method only supports nonorthogonal Hams.")

!          call fermi_operator_annealing(H,S,CCH,temp1,orthogonal,.false.,matmuls)
          if(orthogonal) then

             call matrix_scale(H,1.0_dp,-Dmu/smearing)

             call fermi_operator_annealing(H,S,H,temp1,.true.,.false.,matmuls)
          else
             call matrix_axpy(H,S,-Dmu/smearing)
             call matrix_scale(CCH,1.0_dp,-Dmu/smearing)

             call fermi_operator_annealing(H,S,CCH,temp1,.false.,.false.,matmuls)
          end if
          !       case("HEAD_GORDON")
          !          call utils_assert(orthogonal,"Error in Smearing_changemu: Head-Gordon method only supports orthogonal Hams")
          !          call matrix_scale(H,1.0_dp,-Dmu/smearing)
          !          call fermi_operator_head_gordon(H,temp1,matmuls)
          !  !       case("LIN")
          !             call matrix_scale(H,1.0_dp,-Dmu/smearing)
          !  !          call fermi_operator_pex(H,temp1,matmuls)
       case default
          call utils_abort("Error: Unrecognised FOE type in smearing_changemu")
       end select

       matmuls_reinit=matmuls_reinit+matmuls
       call matrix_scale(temp1,2.0_dp,-1.0_dp)
       call matrix_copy(temp1,T)
       have_invinit=.false.

    else
       method_used="H"!yperbolic identity
       if(pub_on_root.and.debug_info_toggle) write(stdout,*) &
            "smearing-> Modifying Fermi operator to new chemical potential. Delta-mu : ",Dmu,Dmu/(2.0_dp*smearing), "<== SMEAR"

       if(orthogonal) then
          call matrix_scale(H,1.0_dp,-Dmu/smearing)
       else
          call matrix_axpy(H,S,-Dmu/smearing)
          call matrix_scale(CCH,1.0_dp,-Dmu/smearing)
       end if


       if(use_nsh) then
          if(.not.have_invinit) then

             call matrix_transpose(denom,invinit)
             call matrix_scale(invinit,1.0_dp/(matrix_induced_norm(denom)*matrix_induced_norm(invinit)))
          end if
          D=8
          call nsh_inverse(denom,invinit,D,.true.,matmuls)
          have_invinit=.true.

          call matrix_scale(T,1.0_dp,tanh(Dmu/(2.0_dp*smearing)))
          call matrix_multiply(1.0_dp,T,invinit,0.0_dp,temp1)
          matmuls_inv=matmuls_inv+matmuls+1
          call matrix_copy(temp1,T)

       else

!          call matrix_scale(T,1.0_dp,tanh(Dmu/(2.0_dp*smearing)))
!          call matrix_solve(T,denom,temp1)
!          call matrix_copy(temp1,T)

          if(pub_on_root) then
             if(debug_info_toggle) write(stdout,*) "Starting hypsumid:"
          end if
          call get_hypsumid_chebyshev_coeffs(-1.0_dp,1.0_dp,cheby_threshold,Dmu/(2.0_dp*smearing),coeffs)
          if(pub_on_root) then
             !if(debug_info_toggle)
!             write(stdout,*) "NUM COEFFS:",size(coeffs)

          end if

!    if(pub_on_root) then
!       outunit=utils_unit()
!       open(unit=outunit,file="medium_A.mtx")
!    end if
!    call sparse_show_matrix(medium_A%dataSPAM3,outunit,matlab_format=.true.)
!    if(pub_on_root) then
!       close(outunit)
!    end if

          call fast_cheby_resum(T,coeffs,temp1,matmuls)
!          call slow_cheby_resum(T,coeffs,temp1)
          matmuls_inv=matmuls_inv+matmuls
          call matrix_copy(temp1,T)
          deallocate(coeffs,stat=ierr)
          call utils_dealloc_check('fermi_operator_projection_contourint','cheby_coeffs',ierr)

       end if

    end if
    matmuls=matmuls_inv+matmuls_reinit
    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> changemu==>", &
            matmuls_inv,matmuls_reinit,abs(Dmu/(2.0_dp*smearing)),Dmu/(2.0_dp*smearing),&
            & (1.0_dp+abs(tanh(Dmu/(2.0_dp*smearing)))),(1.0_dp-abs(tanh(Dmu/(2.0_dp*smearing)))),condest,norm, &
            & threshold, "<== SMEAR"

    end if

    call matrix_free(temp1)

    call matrix_free(denom)

    call timer_clock('smearing_changemu',2)

  end subroutine changemu


  !============================================================================!
  ! This routine modifies the implicit electron temperature in a density       !
  ! kernel without re-calculating it.                                          !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   T            (inout)  : Scaled & shifted density kernel (tanh matrix)    !
  !   multiply_by  (input)  : Factor to multiply the electron temperature by   !
  !   invinit      (input)  : Initial inverse matrix                           !
  !   have_invinit (input)  : Whether we have invinit or not                   !
  !   matmuls      (output) : number of matmuls used to do everything          !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons... in 2017?                                       !
  !============================================================================!
  subroutine changebeta(T,multiply_by,invinit,have_invinit,matmuls)

    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none
    type(smearing_matrix),  intent(inout) :: T
    integer,       intent(in)    :: multiply_by

    type(smearing_matrix),  intent(inout) :: invinit
    logical,       intent(inout) :: have_invinit
    integer,       intent(out)   :: matmuls

    integer :: N
    type(smearing_matrix) :: T1, T2
    type(smearing_matrix) :: x2,x3,x4,x6,x8,x12,x16
    integer, dimension(:), allocatable :: ci
    real(kind=dp), dimension(:), allocatable :: c
    integer :: inv_matmuls
    integer :: ierr

    call timer_clock('smearing_changebeta',1)

    matmuls=0
    N=matrix_dimension(T)

    allocate(ci(multiply_by+1),stat=ierr)
    call utils_alloc_check('changebeta','ci',ierr)

    call pascal_tri_row(multiply_by,ci)
    allocate(c(multiply_by),stat=ierr)
    call utils_alloc_check('changebeta','c',ierr)
    c=real(ci(2:),dp)
    deallocate(ci,stat=ierr)
    call utils_dealloc_check('changebeta','ci',ierr)


    call matrix_allocate(T1,T)
    call matrix_allocate(T2,T)

    select case(multiply_by)
    case(16)
       ! ---------- 16 ---------
       !
       ! X2 = x*x
       ! x4 = x2*x2
       ! x8 = x4*x4
       ! x12= x4*x8
       ! x16= x4*x12
       !
       ! T1 = I + c8*x8 + c12*x12 + c16*x16 + x2*(c2*I + c4*x2 + c6*x4 + c10*x8 + c14*x12)
       ! T2 = x *(c1*I + c3*x2 + c5*x4 + c9*x8 + c13*x12 + x2*(c7*x4 + c11*x8 + c15*x12))
       !
       ! ans=T1*inv(T2)

       call matrix_allocate(x2,T,T)
       call matrix_allocate(x4,x2,x2)
       call matrix_allocate(x8,x4,x4)
       call matrix_allocate(x12,x4,x8)
       call matrix_allocate(x16,x4,x12)

       call matrix_multiply(1.0_dp,T,T,0.0_dp,x2)
       call matrix_multiply(1.0_dp,x2,x2,0.0_dp,x4)
       call matrix_multiply(1.0_dp,x4,x4,0.0_dp,x8)
       call matrix_multiply(1.0_dp,x4,x8,0.0_dp,x12)
       call matrix_multiply(1.0_dp,x4,x12,0.0_dp,x16)

       call matrix_copy(x4,T1)
       call matrix_scale(T1,c(7))
       call matrix_axpy(T1,x8,c(11))
       call matrix_axpy(T1,x12,c(15))
       call matrix_multiply(1.0_dp,x2,T1,0.0_dp,T2)

       call matrix_scale(T2,1.0_dp,c(1))
       call matrix_axpy(T2,x2,c(3))
       call matrix_axpy(T2,x4,c(5))
       call matrix_axpy(T2,x8,c(9))
       call matrix_axpy(T2,x12,c(13))

       call matrix_multiply(1.0_dp,T,T2,0.0_dp,T1)
       call matrix_copy(T1,T2)

       call matrix_scale(T1,0.0_dp,c(2))
       call matrix_axpy(T1,x2,c(4))
       call matrix_axpy(T1,x4,c(6))
       call matrix_axpy(T1,x8,c(10))
       call matrix_axpy(T1,x12,c(14))
       call matrix_multiply(1.0_dp,x2,T1,0.0_dp,x4)
       call matrix_copy(x4,T1)

       call matrix_scale(T1,1.0_dp,1.0_dp)
       call matrix_axpy(T1,x8,c(8))
       call matrix_axpy(T1,x12,c(12))
       call matrix_axpy(T1,x16,c(16))

       call matrix_free(x2)
       call matrix_free(x4)
       call matrix_free(x8)
       call matrix_free(x12)
       call matrix_free(x16)

       matmuls=matmuls+8

    case(8)
       !
       ! 8
       !
       ! ---------- 8 ---------
       !
       ! x2 = x*x
       ! x4 = x2*x2
       ! x6 = x4*x2
       ! x8 = x4*x4
       !
       ! T1 = I + c2*x2 + c4*x4 + c6*x6 + c8*x8
       ! T2 = x*(c1*I + c3*x2 + c5*x4 + c7*x6)
       !
       ! ans=T1*inv(T2)
       !

       call matrix_allocate(x2,T,T)
       call matrix_allocate(x4,x2,x2)
       call matrix_allocate(x6,x2,x4)
       call matrix_allocate(x8,x4,x4)

       call matrix_multiply(1.0_dp,T,T,0.0_dp,x2)
       call matrix_multiply(1.0_dp,x2,x2,0.0_dp,x4)
       call matrix_multiply(1.0_dp,x2,x4,0.0_dp,x6)
       call matrix_multiply(1.0_dp,x4,x4,0.0_dp,x8)

       call matrix_scale(T1,0.0_dp,c(1))
       call matrix_axpy(T1,x2,c(3))
       call matrix_axpy(T1,x4,c(5))
       call matrix_axpy(T1,x6,c(7))
       call matrix_multiply(1.0_dp,T,T1,0.0_dp,T2)

       call matrix_copy(x2,T1)
       call matrix_scale(T1,c(2),1.0_dp)
       call matrix_axpy(T1,x4,c(4))
       call matrix_axpy(T1,x6,c(6))
       call matrix_axpy(T1,x8,c(8))

       call matrix_free(x2)
       call matrix_free(x4)
       call matrix_free(x6)
       call matrix_free(x8)

       matmuls=matmuls+5

    case(4)
       !
       ! 5
       !
       ! ---------- 4 ---------
       !
       ! x2 = x*x
       ! x3 = x*x2
       ! x4 = x2*x2
       !
       ! T1 = I + c2*x2 + c4*x4
       ! T2 = c1*x + c3*x3
       !
       ! ans=T1*inv(T2)
       !

       call matrix_allocate(x2,T,T)
       call matrix_allocate(x3,T,x2)
       call matrix_allocate(x4,x2,x2)

       call matrix_multiply(1.0_dp,T,T,0.0_dp,x2)
       call matrix_multiply(1.0_dp,T,x2,0.0_dp,x3)
       call matrix_multiply(1.0_dp,x2,x2,0.0_dp,x4)

       call matrix_copy(T,T2)
       call matrix_scale(T2,c(1))
       call matrix_axpy(T2,x3,c(3))

       call matrix_copy(x2,T1)
       call matrix_scale(T1,c(2),1.0_dp)
       call matrix_axpy(T1,x4,c(4))

       call matrix_free(x2)
       call matrix_free(x3)
       call matrix_free(x4)

       matmuls=matmuls+3

    case(2)
       !
       ! 3
       !
       ! ---------- 2 ---------
       !
       ! x2 = x*x
       !
       ! T1 = I + c2*x2
       ! T2 = c1*x
       !
       ! ans=T1*inv(T2)
       !

       call matrix_allocate(x2,T,T)

       call matrix_multiply(1.0_dp,T,T,0.0_dp,x2)
       call matrix_copy(T,T2)
       call matrix_scale(T2,c(1))

       call matrix_copy(x2,T1)
       call matrix_scale(T1,c(2),1.0_dp)

       call matrix_free(x2)

       matmuls=matmuls+1
       ! 1
    end select

    inv_matmuls=0
    if(have_invinit) then
       call matrix_invert(T1,invinit,matmuls=inv_matmuls)
    else
       call matrix_invert(T1,matmuls=inv_matmuls)
       call matrix_copy(T1,invinit)
       have_invinit=.true.
    end if
    matmuls=matmuls+inv_matmuls
    call matrix_multiply(1.0_dp, T1,T2, 0.0_dp, T)
    matmuls=matmuls+1

    call matrix_free(T1)
    call matrix_free(T2)

    deallocate(c,stat=ierr)
    call utils_dealloc_check('changebeta','c',ierr)

    call timer_clock('smearing_changebeta',2)

  end subroutine changebeta


  !============================================================================!
  ! This routine calculates a row from Pascal's triangle.                      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   n            (input)  : which row                                        !
  !   row          (output) : The row                                          !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons... in 2017?                                       !
  !============================================================================!
  subroutine pascal_tri_row(n,row)
    implicit none
    integer, intent(in) :: n
    integer, dimension(:), intent(out) :: row
    integer :: i, j

    row(1)=1
    if(n>0) then
       do i=1,n
          row(2:i+1)=row(1:i)
          row(1)=1
          if(i>1) then
             do j=2,i
                row(j)=sum(row(j:j+1))
             end do
          end if
       end do
    end if

  end subroutine pascal_tri_row


  !============================================================================!
  ! This routine inverts a matrix using the Newton-Schulz-Hotelling method.    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   denom        (input)  : the matrix to invert                             !
  !   invdenom     (inout)  : matrix inverse + initialization                  !
  !   D            (input)  : Precision                                        !
  !   haveinvdenom (input)  : Whether we have inv init or not                  !
  !   matmuls      (output) : number of matmuls used to do everything          !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons... in 2017?                                       !
  !============================================================================!
  ! Checked, works.
  subroutine nsh_inverse(denom,invdenom,D,have_invdenom,matmuls)
    use comms, only : pub_on_root
    use timer, only : timer_clock
    use utils, only : utils_alloc_check, utils_dealloc_check
    implicit none
    type(smearing_matrix), intent(inout) :: denom
    type(smearing_matrix), intent(inout) :: invdenom
    integer,      intent(in)    :: D
    logical,      intent(in)    :: have_invdenom
    integer,      intent(out)   :: matmuls

    logical :: nans
    integer :: miniters!, maxiters
    real(kind=dp) :: conv_test_old, conv_test
    type(smearing_matrix) :: tmp, tmp2
    integer :: i
    integer :: N
    real(kind=dp) :: scale, trace_test
    integer :: recorded_hist_points
    integer, parameter :: max_hist_points=3
    integer, parameter :: num_positive_to_stop=3
    integer :: num_positive_grads
    real(kind=dp), allocatable :: conv_history(:), coeffs(:)
    integer :: ierr

    call timer_clock('smearing_nsh_inverse',1)

    matmuls=0
    nans=.false.

    i=0
    N=matrix_dimension(denom)
    miniters=1
    !    maxiters=1000
    conv_test_old=huge(1.0_dp)

    if(.not.have_invdenom) then
       call matrix_transpose(denom,invdenom)
       scale=1.0_dp/(matrix_induced_norm(denom)*matrix_induced_norm(invdenom))
       call matrix_scale(invdenom,scale)
    end if

    call matrix_allocate(tmp, denom)
    call matrix_allocate(tmp2, denom)

    do
       i=i+1

       ! NLcheb=invdenom
       call matrix_multiply(1.0_dp,denom,invdenom,0.0_dp,tmp)
       call matrix_copy(invdenom,tmp2)
       call matrix_multiply(-1.0_dp,invdenom,tmp,2.0_dp,tmp2)


       matmuls=matmuls+2

       call matrix_multiply(1.0_dp,denom,tmp2,0.0_dp,tmp)
       call matrix_scale(tmp,1.0_dp,-1.0_dp)

       conv_test=abs(matrix_trace(tmp))

       if(conv_test>conv_test_old.or.conv_test/=conv_test) then ! NaN is not equal to itself
          nans=.true.
          exit
       end if

       call matrix_copy(tmp2,invdenom)
       if(i>=miniters) then
          if(conv_test<=10.0_dp**(-D)) then
             exit
          end if
       end if

       conv_test_old=conv_test
    end do


    if(nans) then
       if(pub_on_root) then
          if(debug_info_toggle) write(stdout,*) &
               &"smearing-> Warning: reinitialising inverse in Newton-Shultz-Hotelling.", "<== SMEAR"
       end if

       call matrix_transpose(denom,invdenom)

       scale=1.0_dp/(matrix_induced_norm(denom)*matrix_induced_norm(invdenom)) ! one over the 1 norm * inf norm

       call matrix_scale(invdenom,scale)

       i=0 ! What exactly does maxiters refer to?

       recorded_hist_points=0
       allocate(conv_history(max_hist_points),stat=ierr)
       call utils_alloc_check('smearing_nsh_inverse','conv_history',ierr)
       allocate(coeffs(2),stat=ierr)
       call utils_alloc_check('smearing_nsh_inverse','coeffs',ierr)

       num_positive_grads=0
       do
          i=i+1
          !           tmp=(2.0_dp*NLcheb - (NLcheb*NLprop*NLcheb))
          call matrix_multiply(1.0_dp,denom,invdenom,0.0_dp,tmp)
          call matrix_copy(invdenom,tmp2)
          call matrix_multiply(-1.0_dp,invdenom,tmp,2.0_dp,tmp2)
          matmuls=matmuls+2

          call matrix_copy(tmp2,invdenom)
          if(i>=miniters) then
             call matrix_multiply(1.0_dp,denom,tmp2,0.0_dp,tmp)
             call matrix_scale(tmp,1.0_dp,-1.0_dp)

             trace_test=abs(matrix_trace(tmp))

             recorded_hist_points=recorded_hist_points+1
             if(recorded_hist_points<=max_hist_points) then
                conv_history(recorded_hist_points)=trace_test
                call poly_fit(1, conv_history, coeffs, recorded_hist_points)
             else
                conv_history(1:max_hist_points-1)=conv_history(2:max_hist_points)
                conv_history(max_hist_points)=trace_test
                call poly_fit(1, conv_history, coeffs)
             end if

             if(coeffs(2)>0.0_dp) then
                num_positive_grads=num_positive_grads+1
             else
                num_positive_grads=0
             end if

             if(num_positive_grads>num_positive_to_stop) then
                if(pub_on_root) then
                   if(debug_info_toggle) write(stdout,*) "smearing-> stopping nsh_inverse ==> not converging any more &
                        &: ", i, matmuls, " : ", coeffs,"<== SMEAR"
                end if
                exit
             end if

             if(pub_on_root) then
                if(debug_info_toggle) write(stdout,*) "smearing-> nsh_inverse ==> ", i, matmuls, " : ", coeffs,"<== SMEAR"
             end if

             if(trace_test<=10.0_dp**(-D)) then
                exit
             end if
          end if
       end do
       deallocate(conv_history,stat=ierr)
       call utils_dealloc_check('smearing_nsh_inverse','conv_history',ierr)
       deallocate(coeffs,stat=ierr)
       call utils_dealloc_check('smearing_nsh_inverse','coeffs',ierr)
    end if

    call matrix_free(tmp)
    call matrix_free(tmp2)

    call timer_clock('smearing_nsh_inverse',2)

  end subroutine nsh_inverse

  !============================================================================!
  ! This routine fits a polynomial to data.                                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   fit_order    (input)  : order of polynomial                              !
  !   conv_history (input)  : data to fit to                                   !
  !   coeffs       (output) : polynomial coeffs                                !
  !   hist_points  (input)  : Limit the data to this number of points          !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons... in 2017?                                       !
  !============================================================================!
  subroutine poly_fit(fit_order, conv_history, coeffs, hist_points)
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_safe_nint
    implicit none
    integer,                 intent(in)    :: fit_order
    real(kind=dp),           intent(in)    :: conv_history(:)
    real(kind=dp),           intent(out)   :: coeffs(:)
    integer,       optional, intent(in)    :: hist_points

    ! LAPACK subroutines
    external :: dgetrf, dgetri

    integer :: fi,hj
    integer :: num_hist_points
    real(kind=dp), allocatable :: fitmat(:,:), tmpmat(:,:), fitwork(:)
    integer, allocatable :: ipiv(:)
    integer :: info
    integer :: work_size
    integer :: ierr

    if(present(hist_points)) then
       num_hist_points=hist_points
    else
       num_hist_points=size(conv_history,1)
    end if
    allocate(fitmat(min(num_hist_points,num_hist_points),fit_order+1),stat=ierr)
    call utils_alloc_check('poly_fit','fitmat',ierr)
    allocate(tmpmat(size(fitmat,2),size(fitmat,2)),stat=ierr)
    call utils_alloc_check('poly_fit','tmpmat',ierr)
    allocate(ipiv(size(tmpmat,1)),stat=ierr)
    call utils_alloc_check('poly_fit','ipiv',ierr)

    do fi = 0, fit_order
       do hj = 1, min(num_hist_points,num_hist_points)
          fitmat(hj,fi+1) = real(hj**fi,dp)
       end do
    end do
    tmpmat=matmul(transpose(fitmat),fitmat)
    call dgetrf(fit_order+1, fit_order+1, tmpmat, fit_order+1, ipiv, info)
    work_size=-1
    allocate(fitwork(1),stat=ierr)
    call utils_alloc_check('poly_fit','fitwork',ierr)
    call dgetri(fit_order+1, tmpmat, fit_order+1, ipiv, fitwork, work_size, info)
    work_size = utils_safe_nint(fitwork(1))
    deallocate(fitwork,stat=ierr)
    call utils_dealloc_check('poly_fit','fitwork',ierr)
    allocate(fitwork(work_size),stat=ierr)
    call utils_alloc_check('poly_fit','fitwork',ierr)
    call dgetri(fit_order+1, tmpmat, fit_order+1, ipiv, fitwork, work_size, info)
    deallocate(fitwork,stat=ierr)
    call utils_dealloc_check('poly_fit','fitwork',ierr)

    coeffs(1:fit_order+1)=matmul(matmul(tmpmat, transpose(fitmat)), conv_history(1:num_hist_points))

    deallocate(fitmat,stat=ierr)
    call utils_dealloc_check('poly_fit','fitmat',ierr)
    deallocate(tmpmat,stat=ierr)
    call utils_dealloc_check('poly_fit','tmpmat',ierr)
    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('poly_fit','ipiv',ierr)

  end subroutine poly_fit


  !============================================================================!
  ! A rough approximation to the Fermi-Dirac entropy of a Hamiltonian of       !
  ! finite-temperature, non-interacting particles                              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   K             (inout)  : Density kernel                                  !
  !   S             (input)  : The overlap matrix.                             !
  !   inv_init      (input)  : Initial inverse                                 !
  !   invinitrefine (input)  : Initial inverse (for refined approximation)     !
  !   refine        (input)  : perform refined approximation                   !
  !   cocontra_K    (input)  : The co-contravariant kernel                     !
  !   entropy       (output) : The entropy                                     !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons... in 2017?                                       !
  !============================================================================!
  function fermi_entropy_approx(K,S,inv_init,inv_init_refine,refine,cocontra_K) result(entropy)
    use comms,only:pub_on_root
    use rundat,only:pub_foe_check_entropy
    use utils,only:utils_dealloc_check
    type(smearing_matrix),           intent(inout) :: K
    type(smearing_matrix), optional, intent(inout) :: S
    type(smearing_matrix), optional, intent(inout) :: inv_init
    type(smearing_matrix), optional, intent(inout) :: inv_init_refine
    logical,      optional, intent(in)    :: refine
    logical,      optional, intent(in)    :: cocontra_K
    real(kind=dp)                         :: entropy

    real(kind=dp) :: quad_approx_entropy

    type(smearing_matrix) :: K2, numer, accumu
    integer :: matmuls, inv_matmuls, matmuls_refine

    type(smearing_matrix) :: inv_sqrt_S, inv_sqrt_S_trans, tmp_mat, K_loc
    integer :: lowdin_maxiters

    logical :: cocontra_K_local

    real(kind=dp), dimension(:), allocatable :: coeffs
    integer :: ierr
    real(kind=dp) :: threshold

    threshold=1e-8_dp

    cocontra_K_local=.false.
    if(present(cocontra_K)) then
       cocontra_K_local=cocontra_K
    end if

    call matrix_allocate(K_loc, K)

    if(present(S).and.(.not.cocontra_K_local)) then
       lowdin_maxiters = 100 ! Why (not)?
       call matrix_allocate(inv_sqrt_S,S)
       call lowdin_transformation(S,sqrt_S_out=inv_sqrt_S,maxiters=lowdin_maxiters)
       call matrix_allocate(inv_sqrt_S_trans,S)
       call matrix_transpose(inv_sqrt_S,inv_sqrt_S_trans)

       call matrix_allocate(tmp_mat,inv_sqrt_S_trans,K)
       call matrix_multiply(1.0_dp,inv_sqrt_S_trans,K,0.0_dp,tmp_mat)
       call matrix_free(inv_sqrt_S_trans)
       call matrix_multiply(1.0_dp,tmp_mat,inv_sqrt_S,0.0_dp,K_loc)
       call matrix_free(inv_sqrt_S)
       call matrix_free(tmp_mat)
    else
       call matrix_copy(K,K_loc)
    end if

 quad_approx_entropy = matrix_trace(K_loc,K_loc) ! KSKS
 quad_approx_entropy = quad_approx_entropy - matrix_trace(K_loc) ! KSKS - KS
 quad_approx_entropy = 3.0_dp*quad_approx_entropy


    !    A simple model for the Fermi-dirac entropy of non-interacting, finite-temperature
    !    particles is given as:
    !    y = 1.96056*x^2 + 0.0286723/(0.114753 + 1.9888*x -1.9888*x^2) - 0.24986 -1.96056*x

    !ja531->sparsefoe:
!    call matrix_allocate(K2,K_loc,K_loc)

if(.not.cocontra_K_local) then

    call matrix_allocate(K2,K_loc)

    call matrix_allocate(numer,K2)
    call matrix_allocate(accumu,K2)


    call matrix_multiply(1.0_dp,K_loc,K_loc,0.0_dp,K2)
    matmuls=1

    call matrix_copy(K2,numer)

    call matrix_scale(numer,-1.9888_dp,0.114753_dp)

    call matrix_axpy(numer,K_loc,1.9888_dp)

    if(.not.cocontra_K_local) then
       if(present(inv_init)) then
          call matrix_invert(numer, inv_init, matmuls=inv_matmuls)
          call matrix_copy(numer, inv_init)
       else
          call matrix_invert(numer,matmuls=inv_matmuls)
       end if
       matmuls=matmuls+inv_matmuls
    else
       if(pub_on_root) then
          write(stdout,*) "Entering Entropy approx Solve"
       end if
       call matrix_solve(S,numer,accumu)

       if(pub_on_root) then
          write(stdout,*) "Entering Entropy approx contra Solve"
       end if
       call matrix_solve(accumu,S,numer)
    end if


    call matrix_copy(numer,accumu)

    call matrix_scale(accumu,0.0286723_dp,-0.24986_dp)

    call matrix_axpy(accumu,K_loc,-1.96056_dp)

    call matrix_axpy(accumu,K2,1.96056_dp)

    call matrix_free(K2)
    call matrix_free(numer)

 else
    call matrix_allocate(accumu,K_loc)

    call get_approxentropy_chebyshev_coeffs(-1.0_dp,1.0_dp,threshold,coeffs)

    if(pub_on_root) then
       if(debug_info_toggle) then
          write(stdout,*) "smearing-> Entropy coeffs:", coeffs
       end if
    end if

    call matrix_scale(K_loc,2.0_dp,-1.0_dp)
    call fast_cheby_resum(K_loc,coeffs,accumu,matmuls,force_layers=1)

    deallocate(coeffs,stat=ierr)
    call utils_dealloc_check('fermi_operator_projection_contourint','cheby_coeffs',ierr)

 end if


    entropy = matrix_trace(accumu)

    if(present(refine)) then
       if(refine) then
          if(present(inv_init_refine)) then
             entropy = fermi_entropy_refineapprox(accumu,inv_init_refine,matmuls_ref=matmuls_refine)
          else
             entropy = fermi_entropy_refineapprox(accumu,matmuls_ref=matmuls_refine)
          end if
          matmuls=matmuls+matmuls_refine
       end if
    end if


 !ja531-> check entropy:
 ! the quadratic f(x) = ((x-0.5)^2 - 0.25)*3.0
 ! (x-0.5)(x-0.5) - 0.25 = x^2 - x + 0.25 - 0.25
 ! f(x) = 3x^2 - 3x
 ! the coefficients of 3 minimise the integral of x*Ln(x) + (1-x)*Ln(1-x) - cx^2 - cx over [0:1]
 ! The maximum error of this approximation is at 0.5, where the error is 0.05685 or 8.2% of the
 ! correct value at 0.5.
 ! This indicates that in the very worst case where every eigenvalue of KS is 0.5, then the error
 ! on the entropy from this approximation is 8.2%
 ! The value of this is that KSKS is easy to calculate reliably, compared with a high order series
 ! expansion... i.e. if there is an eigenvalue of K which is a little above 1 because of noise due
 ! to truncation, the corresponding eigenvalues in the entropy matrix can be very large. This can
 ! accumulate in the entropy value, if there are many of these greater than 1 eigenvalues, to give
 ! a completely wrong value.
 ! We can, therefore, check the value of entropy we have calculated by comparing it with this
 ! quadratic approximation. If the absolute difference between the two values is less than 8.2% of
 ! the high order value, then we assume that things have gone well and do nothing. If this is not
 ! the case, then there was a high amount of error accumulated and so we'll print a warning and use
 ! the quadratic approximation.

    if(pub_FOE_check_entropy) then
       if(abs(quad_approx_entropy-entropy)/abs(entropy) > 0.082_dp) then
          if(pub_on_root) then
             write(stdout,*) "WARNING: entropy calculation accumulated a lot of error, using quad. approx."
             if(debug_info_toggle) write(stdout,*) "WARNING: bad entropy=",entropy, " quad approx=",quad_approx_entropy
          end if
          entropy=quad_approx_entropy
       else
          if(pub_on_root) then
             if(debug_info_toggle) write(stdout,*) "INFO: good entropy=",entropy, " quad approx=",quad_approx_entropy
          end if
       end if
    end if

    if(pub_on_root) then
       if(debug_info_toggle) write(stdout,*) "smearing-> Fermi_entropy_approx matmuls : ", matmuls
    end if

    call matrix_free(accumu)
    call matrix_free(K_loc)
  end function fermi_entropy_approx

  !ep
  !============================================================================!
  ! A Mermin method compliant approximation to the Fermi-Dirac entropy of a    !
  !  Hamiltonian of finite-temperature, non-interacting particles              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   K             (inout)  : Density kernel                                  !
  !   S             (input)  : The overlap matrix.                             !
  !   inv_init      (input)  : Initial inverse                                 !
  !   invinitrefine (input)  : Initial inverse (for refined approximation)     !
  !   refine        (input)  : perform refined approximation                   !
  !   cocontra_K    (input)  : The co-contravariant kernel                     !
  !   entropy       (output) : The entropy                                     !
  !----------------------------------------------------------------------------!
  ! Adapted by Emliano Poli from Jolyon Aarons version to work with the        !
  ! Mermin method                                                             !
  !============================================================================!
  function fermi_entropy_approx_mermin(K,S, expans, &
    inv_init,inv_init_refine,refine,cocontra_K) result(entropy)
    use constants, only: k_B
    use utils, only: utils_unit
    use comms,only:pub_on_root
    use rundat,only:pub_foe_check_entropy, pub_mermin_smearing_width
    use utils,only:utils_dealloc_check
    type(smearing_matrix),           intent(inout) :: K
    type(smearing_matrix), optional, intent(inout) :: S
    type(smearing_matrix), optional, intent(inout) :: inv_init
    type(smearing_matrix), optional, intent(inout) :: inv_init_refine
    logical,      optional, intent(in)    :: refine
    logical,      optional, intent(in)    :: cocontra_K
    integer,      optional, intent(in)    :: expans
    real(kind=dp)                         :: entropy
    real(kind=dp)                         :: tmp_pop
    real(kind=dp)                         :: tmp_pop2
    real(kind=dp) :: quad_approx_entropy
    type(smearing_matrix) ::  numer, accumu, accumu2
    integer :: matmuls, inv_matmuls, matmuls_refine
    type(smearing_matrix) :: inv_sqrt_S, inv_sqrt_S_trans, tmp_mat, K_loc
    integer :: lowdin_maxiters , outunit
    logical :: cocontra_K_local
    real(kind=dp), dimension(:), allocatable :: coeffs
    integer :: ierr
    real(kind=dp) :: threshold

    threshold=1e-8_dp
    cocontra_K_local=.true.
    if(present(cocontra_K)) then
       cocontra_K_local=cocontra_K
    end if

    call matrix_allocate(K_loc, K)

    if(present(S).and.(.not.cocontra_K_local)) then
       lowdin_maxiters = 100 ! Why (not)?
       call matrix_allocate(inv_sqrt_S,S)
       call lowdin_transformation(S,sqrt_S_out=inv_sqrt_S, &
            maxiters=lowdin_maxiters)
       call matrix_allocate(inv_sqrt_S_trans,S)
       call matrix_transpose(inv_sqrt_S,inv_sqrt_S_trans)

       call matrix_allocate(tmp_mat,inv_sqrt_S_trans,K)
       call matrix_multiply(1.0_dp,inv_sqrt_S_trans,K,0.0_dp,tmp_mat)
       call matrix_free(inv_sqrt_S_trans)
       call matrix_multiply(1.0_dp,tmp_mat,inv_sqrt_S,0.0_dp,K_loc)
       call matrix_free(inv_sqrt_S)
       call matrix_free(tmp_mat)
    else
       call matrix_multiply(1.0_dp,K,S,0.0_dp,K_loc)
    end if

    call matrix_allocate(accumu,K_loc)
    call get_approxentropy_chebyshev_coeffs(0.0_dp,1.0_dp,threshold, &
         coeffs,expans)
    call matrix_scale(K_loc,2.0_dp,-1.0_dp)
    call fast_cheby_resum(K_loc,coeffs,accumu,matmuls,force_layers=1)
    deallocate(coeffs,stat=ierr)
    call utils_dealloc_check('fermi_entropy_approx_mermin', &
         'coeffs',ierr)
    !end if

    entropy = matrix_trace(accumu)
    entropy = k_B * entropy

    if(present(refine)) then
       if(refine) then
          if(present(inv_init_refine)) then
             entropy = fermi_entropy_refineapprox(accumu,inv_init_refine,matmuls_ref=matmuls_refine)
          else
             entropy = fermi_entropy_refineapprox(accumu,matmuls_ref=matmuls_refine)
          end if
          matmuls=matmuls+matmuls_refine
       end if
    end if

    call matrix_free(accumu)
    call matrix_free(K_loc)
  end function fermi_entropy_approx_mermin
  !ep

  !ep
  !============================================================================!
  ! An  approximation to the Derivative Fermi-Dirac entropy                    !
  ! of a Hamiltonian of finite-temperature, non-interacting particles w.r.t    !
  ! the Density Kernel                                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   K             (inout)  : Density kernel                                  !
  !   S             (input)  : The overlap matrix.                             !
  !   inv_init      (input)  : Initial inverse                                 !
  !   invinitrefine (input)  : Initial inverse (for refined approximation)     !
  !   refine        (input)  : perform refined approximation                   !
  !   cocontra_K    (input)  : The co-contravariant kernel                     !
  !   derivative    (output) : Derivative  of entropy                          !
  !----------------------------------------------------------------------------!
  ! Written by Poli Emiliano                                                   !
  !============================================================================!
  subroutine fermi_deriv_approx(K,derivative,entropy,S,expans, &
    inv_init,inv_init_refine,refine,cocontra_K)
    use constants, only: k_B
    use utils, only: utils_unit
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy, &
        sparse_embed_create, sparse_embed_scale
    use sparse, only: SPAM3, sparse_show_matrix
    use comms,only:pub_on_root
    use rundat,only:pub_foe_check_entropy
    use utils,only:utils_dealloc_check

    type(smearing_matrix),           intent(inout) :: K
    type(SPAM3_EMBED),               intent(inout) :: derivative
    real(kind=dp) ,                  intent(inout) :: entropy
    type(smearing_matrix), optional, intent(inout) :: S
    type(smearing_matrix), optional, intent(inout) :: inv_init
    type(smearing_matrix), optional, intent(inout) :: inv_init_refine
    integer,      optional, intent(in)    :: expans
    logical,      optional, intent(in)    :: refine
    logical,      optional, intent(in)    :: cocontra_K
    type(smearing_matrix) :: K2, numer, accumu, accumu2, accumu3
    integer :: matmuls, inv_matmuls, matmuls_refine
    type(smearing_matrix) :: inv_sqrt_S, inv_sqrt_S_trans, &
        tmp_mat, K_loc, tmp_invo
    integer :: lowdin_maxiters , outunit
    logical :: cocontra_K_local
    real(kind=dp), dimension(:), allocatable :: coeffs
    real(kind=dp), dimension(:), allocatable :: dcoeffs
    integer :: ierr
    real(kind=dp) :: threshold

    threshold=1e-8_dp
    cocontra_K_local=.true.
    if(present(cocontra_K)) then
       cocontra_K_local=cocontra_K
    end if

    call matrix_allocate(K_loc, K)

    if(present(S).and.(.not.cocontra_K_local)) then
       lowdin_maxiters = 100 ! Why (not)?
       call matrix_allocate(inv_sqrt_S,S)
       call lowdin_transformation(S,sqrt_S_out=inv_sqrt_S, &
            maxiters=lowdin_maxiters)
       call matrix_allocate(inv_sqrt_S_trans,S)
       call matrix_transpose(inv_sqrt_S,inv_sqrt_S_trans)

       call matrix_allocate(tmp_mat,inv_sqrt_S_trans,K)
       call matrix_multiply(1.0_dp,inv_sqrt_S_trans,K,0.0_dp,tmp_mat)
       call matrix_free(inv_sqrt_S_trans)
       call matrix_multiply(1.0_dp,tmp_mat,inv_sqrt_S,0.0_dp,K_loc)
       call matrix_free(inv_sqrt_S)
       call matrix_free(tmp_mat)
    else
       call matrix_multiply(1.0_dp,K,S,0.0_dp,K_loc)
    end if

    call matrix_allocate(accumu,K_loc)
    call matrix_allocate(accumu2,K_loc)

    call get_approxentropy_chebyshev_coeffs(0.0_dp,1.0_dp,threshold, &
         coeffs, expans, dcoeffs)
    call matrix_scale(K_loc,2.0_dp,-1.0_dp)
    call fast_cheby_resum(K_loc,coeffs,accumu,matmuls,force_layers=1)
    call fast_cheby_resum(K_loc,dcoeffs,accumu2,matmuls,force_layers=1)
    call sparse_embed_scale(accumu2%dataSPAM3_EMBED,k_B,0.0_DP)
    call sparse_embed_copy(derivative,accumu2%dataSPAM3_EMBED)

    deallocate(coeffs,stat=ierr)
    deallocate(dcoeffs,stat=ierr)

    entropy = matrix_trace(accumu)
    entropy = k_B * entropy

    call matrix_free(accumu)
    call matrix_free(accumu2)
    call matrix_free(K_loc)

  end subroutine fermi_deriv_approx
  !ep

  !ep
  !============================================================================!
  ! An  approximation to the Derivative Fermi-Dirac entropy                    !
  ! of a Hamiltonian of finite-temperature, non-interacting particles  w.r.t.  !
  ! to the NGWFs                                                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   K             (inout)  : Density kernel                                  !
  !   derivative    (output) : Derivative of entropy w.r.t. of NGWFs           !
  !   S             (input)  : The overlap matrix.                             !
  !   expans        (input)  : Chebyshev expansion of the approximation        !
  !   Sinv          (input)  : Inverse of the overlap matrix                   !
  !----------------------------------------------------------------------------!
  ! Written by Poli Emiliano                                                   !
  !============================================================================!
  subroutine deriv_ngwf_approx(K,derivative,S,expans,Sinv)
    use constants, only: k_B
    use utils, only: utils_unit
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy, &
        sparse_embed_create, sparse_embed_scale
    use sparse, only: SPAM3, sparse_show_matrix
    use comms,only:pub_on_root
    use rundat,only:pub_foe_check_entropy
    use utils,only:utils_dealloc_check

    type(smearing_matrix),           intent(inout) :: K
    type(SPAM3_EMBED),               intent(inout) :: derivative
    integer,      optional, intent(in)    :: expans
    type(smearing_matrix), intent(inout) :: S
    type(smearing_matrix), intent(inout) :: Sinv
    type(smearing_matrix) :: K2, numer, accumu2
    integer :: matmuls, inv_matmuls, matmuls_refine
    type(smearing_matrix) :: inv_S, inv_sqrt_S_trans, tmp_mat, K_loc, Sloc
    integer :: lowdin_maxiters , outunit
    logical :: cocontra_K_local
    real(kind=dp), dimension(:), allocatable :: coeffs
    real(kind=dp), dimension(:), allocatable :: dcoeffs
    integer :: ierr
    real(kind=dp) :: threshold , dentrop

    threshold=1e-8_dp
    call matrix_allocate(K_loc, K)
    call matrix_allocate(K2, K)
    call matrix_allocate(Sloc, S)
    call matrix_allocate(inv_S, Sinv)
    call matrix_copy(K,K2)
    call matrix_copy(S,Sloc)
    call matrix_copy(Sinv, inv_S)
    call matrix_multiply(1.0_dp,K,S,0.0_dp,K_loc)

    call matrix_allocate(accumu2,K)
    call get_approxentropy_chebyshev_coeffs(0.0_dp,1.0_dp,threshold, &
        coeffs, expans,dcoeffs)
    !ep: test routines
    !call slow_cheby_deriv(K_loc,K2,coeffs,accumu2, Sinv)
    !call matrix_scale(K_loc,2.0_dp,-1.0_dp)
    !call hardcode_deriv(K_loc,K2,coeffs,accumu2, Sinv)
    !ep: test routines
    call slow_deriv(K_loc,K2,coeffs,accumu2)
    call sparse_embed_scale(accumu2%dataSPAM3_EMBED,k_B,0.0_DP)
    call sparse_embed_copy(derivative,accumu2%dataSPAM3_EMBED)

    call matrix_free(accumu2)
    call matrix_free(K_loc)
    call matrix_free(K2)

  end subroutine deriv_ngwf_approx
  !ep

  !ep
  !============================================================================!
  ! This routine calculates a Chebyshev expansion of a matrix modified to      !
  ! comply with Entropy derivative formula.                                    !
  ! It needs N matrix multiplications, where N is the length of the expansion. !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   H            (input)  : The input matrix(KS).                            !
  !   Ker          (input)  : The Kernel matrix(K).                            !
  !   cheby_coeffs (input)  : The coefficients of the expansion.               !
  !   rho          (output) : The result                                       !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons and modified by Emiliano Poli to calculated the   !
  ! explicit derivative of the Entropy w.r.t the NGWFs in 2022                 !
  !============================================================================!


  subroutine slow_deriv(H,Ker,cheby_coeffs,rho)
    use utils, only: utils_unit
    use comms,only:pub_on_root
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy, sparse_embed_create
    use sparse, only: SPAM3, sparse_show_matrix

    implicit none

    type(smearing_matrix), intent(inout)  :: H
    type(smearing_matrix), intent(inout)  :: Ker
    real(kind=dp), dimension(:), intent(in) :: cheby_coeffs
    type(smearing_matrix), intent(inout) :: rho
    character(len=100) :: outfile, outfile2, outfile3, outfile4
    character(len=10)  :: cycle_char
    type(smearing_matrix) :: T0, T1, T, HS, T0p, T1p, Tp, rhop
    type(smearing_matrix) :: tmp, tmp2, tmp3, tmp4
    type(smearing_matrix) :: x, dx, kerloc
    integer :: N, i
    integer       :: outunit
    real(kind=dp) :: Ntmp

    N=size(cheby_coeffs)

    call matrix_allocate(kerloc,Ker)
    call matrix_copy(Ker,kerloc)
    call matrix_scale(rho,0.0_dp)

    if(N>0) then
       call matrix_allocate(T0,rho)
       call matrix_allocate(T0p,rho)
       call matrix_scale(T0,0.0_dp,1.0_dp)
       call matrix_scale(T0p,0.0_dp,0.0_dp)
       call matrix_axpy(rho,T0p,cheby_coeffs(1))
    end if

    if(N>1) then
       call matrix_allocate(T1p,rho)
       call matrix_allocate(T1,rho)
       call matrix_allocate(x,H)
       call matrix_allocate(dx,Ker)
       call matrix_copy(H,x)
       call matrix_scale(x,4.0_dp,-2.0_dp) ! x = 4KS-2I
       call matrix_copy(kerloc, dx)
       call matrix_scale(dx,4.0_dp,0.0_dp) ! dx = 4K
       call matrix_copy(kerloc,T1p)
       call matrix_scale(T1p,2.0_dp,0.0_dp) ! T1p = 2K
       call matrix_copy(H,T1)
       call matrix_scale(T1,2.0_dp,-1.0_dp) ! Tp = 2KS-I
       call matrix_axpy(rho,T1p,cheby_coeffs(2))
    end if

    if(N>2) then
       call matrix_allocate(T,rho)
       call matrix_allocate(tmp,rho)
       call matrix_allocate(tmp2,rho)
       call matrix_allocate(Tp,rho)
       do i=3, N
          !POLYN
          call matrix_scale(tmp,0.0_dp,0.0_dp)
          call matrix_scale(tmp2,0.0_dp,0.0_dp)
          call matrix_copy(T0,T)
          call matrix_multiply(1.0_dp,x,T1,-1.0_dp,T)
          !POLYN
          !DERIV
          !call matrix_copy(T0p,Tp) ! T0' ===> T'
          call matrix_multiply(4.0_dp,T1,kerloc,1.0_dp,tmp) ! 4KS - 2 * T1' ===> tmp
          call matrix_multiply(4.0_dp,H,T1p,1.0_dp,tmp2) ! 4KS - 2 * T1' ===> tmp
          call matrix_axpy(tmp,tmp2,1.0_dp) !  tmp - T0' ===> tmp
          call matrix_axpy(tmp,T1p,-2.0_dp) !  tmp - T0' ===> tmp
          call matrix_axpy(tmp,T0p,-1.0_dp) !  tmp - T0' ===> tmp
          call matrix_copy(tmp,Tp) ! tmp ===> T'
          call matrix_axpy(rho,Tp,cheby_coeffs(i)) ! rho ===> T' * coeff
          !DERIV
          !COPY
          call matrix_copy(T1,T0) ! T1 ===> T0
          call matrix_copy(T,T1) ! T ===> T1
          call matrix_copy(T1p,T0p) ! T1' ===> T0'
          call matrix_copy(Tp,T1p) ! T' ===> T1'
       end do
    end if
    call matrix_free(T)
    call matrix_free(Tp)
    call matrix_free(tmp)
    call matrix_free(tmp2)


    if(N>1) then
       call matrix_free(T1)
       call matrix_free(T1p)
       call matrix_free(x)
       call matrix_free(dx)
    end if
    if(N>0) then
       call matrix_free(T0)
       call matrix_free(T0p)
    end if

    call matrix_free(rhop)
    call matrix_free(kerloc)

  end subroutine slow_deriv
  !ep


  !============================================================================!
  ! A refined approximation to the Fermi-Dirac entropy of a Hamiltonian of     !
  ! finite-temperature, non-interacting particles                              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   entmat        (input)  : matrix from rough approximation                 !
  !   inv_init      (input)  : inverse initialization                          !
  !   matmuls       (output) : number of matmuls required                      !
  !   entropy       (output) : the entropy                                     !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons... in 2017?                                       !
  !============================================================================!
  function fermi_entropy_refineapprox(entmat,inv_init,matmuls_ref) result(entropy)
    use comms,only:pub_on_root
    type(smearing_matrix),           intent(inout) :: entmat
    type(smearing_matrix), optional, intent(inout) :: inv_init
    integer,      optional, intent(out)   :: matmuls_ref
    real(kind=dp)                         :: entropy

    type(smearing_matrix) :: K2, K3, K4, numer, denom, accumu
    integer :: matmuls, inv_matmuls

    type(smearing_matrix) :: inv_sqrt_S, inv_sqrt_S_trans, tmp_mat
    integer :: lowdin_maxiters

    matmuls=0

    !    y = x + (0.0154879202579394*x + x^4 + 1.15423486644271*x^3 + 0.341889372337663*x^2 &
    !         - 4.96444747147679e-6)/(0.0444653992157101 - 1.80660650711283*x - 1.15423486644271*x^2 &
    !         - 5.06114049570097*x^3)

    !ja531->sparsefoe:
!    call matrix_allocate(K2,entmat,entmat)
    call matrix_allocate(K2,entmat)

    call matrix_multiply(1.0_dp,entmat,entmat,0.0_dp,K2)
    matmuls=matmuls+1
    !ja531->sparsefoe:
!    call matrix_allocate(K3,K2,entmat)
    call matrix_allocate(K3,entmat)

    call matrix_multiply(1.0_dp,K2,entmat,0.0_dp,K3)
    matmuls=matmuls+1
    !ja531->sparsefoe:
!    call matrix_allocate(K4,K3,entmat)
    call matrix_allocate(K4,entmat)

    call matrix_multiply(1.0_dp,K3,entmat,0.0_dp,K4)
    matmuls=matmuls+1

    call matrix_allocate(numer,K4)
    call matrix_allocate(denom,K3)

    call matrix_copy(entmat, numer)
    call matrix_copy(entmat, denom)

    call matrix_scale(numer,0.0154879202579394_dp,-4.96444747147679E-6_dp)
    call matrix_scale(denom,-1.80660650711283_dp,0.0444653992157101_dp)

    call matrix_axpy(numer,K2,0.341889372337663_dp)
    call matrix_axpy(denom,K2,-1.15423486644271_dp)

    call matrix_free(K2)
    call matrix_axpy(numer,K3,1.15423486644271_dp)
    call matrix_axpy(denom,K3,-5.06114049570097_dp)

    call matrix_free(K3)
    call matrix_axpy(numer,K4,1.0_dp)
    call matrix_free(K4)

    if(present(inv_init)) then
       call matrix_invert(denom, inv_init, matmuls=inv_matmuls)
       call matrix_copy(denom, inv_init)
    else
       call matrix_invert(denom,matmuls=inv_matmuls)
    end if
    matmuls=matmuls+inv_matmuls

    !ja531->sparsefoe:
!    call matrix_allocate(accumu,numer,denom)
    call matrix_allocate(accumu,entmat)

    call matrix_multiply(1.0_dp,numer,denom,0.0_dp,accumu)
    matmuls=matmuls+1

    call matrix_axpy(accumu,entmat,1.0_dp)

    entropy = matrix_trace(accumu)


    if(.not.present(matmuls_ref)) then
       if(pub_on_root) then
          write(stdout,*) "Entropy refined to : ", entropy
          if(debug_info_toggle) write(stdout,*) "Fermi_entropy_refine matmuls : ", matmuls
       end if
    else
       matmuls_ref=matmuls
    end if


    call matrix_free(accumu)
    call matrix_free(denom)
    call matrix_free(numer)
  end function fermi_entropy_refineapprox


  !============================================================================!
  ! An implementation of the contour from Lin Lin's PEXSI method.              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons... in 2017?                                       !
  !============================================================================!
  subroutine linlin_contour(beta,em,Q)
    implicit none
    real(kind=dp), intent(in) :: beta
    real(kind=dp), intent(in) :: em   ! eigenvalue spectrum width i.e. -4....4
    integer, intent(in) :: Q          ! number of points
    real(kind=dp) :: m, MM, k, kp, KK, KKp
    complex(kind=dp) :: sn,cn,dn, t, z
    complex(kind=dp) :: eps_p, eps_m, g_p, g_m

    integer :: j, Qh

    Qh=Q/2

    m=(pi/beta)**2
    MM=(em**2)+m
    k=(sqrt(MM/m)-1.0_dp)/(sqrt(MM/m)+1.0_dp)
    kp=sqrt(1.0_dp - k**2)

    KK=complete_elliptic_integral(k)
    KKp=complete_elliptic_integral(kp)

    write(*,*) m,MM,k,kp,KK,KKp

!    KKp=KKp/2.0_dp
!    write(*,*) "KK->",KK, KKp

    do j=1,Qh
       t=cmplx(-KK + 2.0_dp*(real(j,dp)-0.5_dp)*KK/real(Qh,dp),KKp/2.0_dp,dp)

!       if(aimag(t)>0.5_dp*KKp) then
!          call jacobi_elliptic(KKp-t,1.0_dp-(k**2),sn,cn,dn)
!       else
          call jacobi_elliptic(t,1.0_dp-(kp**2),sn,cn,dn)
!       end if

       call get_z(z,k,sn,m,MM)
       eps_p=sqrt(z-m)
       write(*,*) eps_p
       eps_m=-eps_p
!       write(69,*) real(eps_p,dp), aimag(eps_p)
!       write(69,*) real(eps_m,dp), aimag(eps_m)


       ! Intel Fortran...
!       g_p=tanh(beta*eps_p/2.0_dp)
!       g_m=tanh(beta*eps_m/2.0_dp)

       g_p=(exp(beta*eps_p/2.0_dp) - exp(-beta*eps_p/2.0_dp))/(exp(beta*eps_p/2.0_dp) + exp(-beta*eps_p/2.0_dp))
       g_m=(exp(beta*eps_m/2.0_dp) - exp(-beta*eps_m/2.0_dp))/(exp(beta*eps_m/2.0_dp) + exp(-beta*eps_m/2.0_dp))


    end do

    do j=Qh,1,-1
       t=cmplx(-KK + 2.0_dp*(real(j,dp)-0.5_dp)*KK/real(Qh,dp),KKp/2.0_dp,dp)
       call jacobi_elliptic(t,1.0_dp-(kp**2),sn,cn,dn)

       call get_z(z,k,sn,m,MM)
       eps_p=sqrt(z-m)
!       eps_m=-eps_p
       write(69,*) -real(eps_p,dp), aimag(eps_p)
    end do

    do j=1,Qh
       t=cmplx(-KK + 2.0_dp*(real(j,dp)-0.5_dp)*KK/real(Qh,dp),KKp/2.0_dp,dp)
       call jacobi_elliptic(t,1.0_dp-(kp**2),sn,cn,dn)

       call get_z(z,k,sn,m,MM)
       eps_p=sqrt(z-m)
!       eps_m=-eps_p
       write(69,*) real(eps_p,dp), aimag(eps_p)
    end do



  end subroutine linlin_contour

  subroutine get_z(z,k,sn,m,MM)
    implicit none
    complex(kind=dp), intent(out) :: z
    real(kind=dp), intent(in) :: k
    complex(kind=dp), intent(in) :: sn
    real(kind=dp), intent(in) :: m
    real(kind=dp), intent(in) :: MM
    real(kind=dp)             :: km1

    km1=1.0_dp/k
    z=sqrt(m*MM)*(km1+sn)/(km1-sn)

  end subroutine get_z

  subroutine jacobi_elliptic_complex(uu,emmc,sn,cn,dn)
    implicit none
    complex(kind=dp), intent(in)  :: uu
    real(kind=dp),    intent(in)  :: emmc
    complex(kind=dp), intent(out) :: sn,cn,dn
    real(kind=dp) :: emcp, denom, cn1,dn1,sn1,cn2,dn2,sn2

    emcp = 1.0_dp - emmc

    call jacobi_elliptic_real(real(uu,dp),emcp,sn1,cn1,dn1)
    call jacobi_elliptic_real(aimag(uu),emmc,sn2,cn2,dn2)

    denom = 1.0_dp - (dn1**2)*(sn2**2)

    sn=cmplx(sn1*dn2,cn1*dn1*sn2*cn2,dp)/denom
    cn=cmplx(cn1*cn2,-sn1*dn1*sn2*dn2,dp)/denom
    dn=cmplx(dn1*cn2*dn2,-sn1*cn1*sn2*emmc,dp)/denom

  end subroutine jacobi_elliptic_complex

  subroutine jacobi_elliptic_real(uu,emmc,sn,cn,dn)
    implicit none
    real(kind=dp), intent(in)  :: uu, emmc
    real(kind=dp), intent(out) :: cn,dn,sn
    real(kind=dp), parameter :: CA=0.0000003_dp ! The accuracy is the square of CA.
    !Returns the Jacobian elliptic functions sn(u,kc), cn(u,kc), and dn(u,kc). Here uu = u,
    !while emmc = kc2 = 1-k^2
    integer :: i,ii,l
    real(kind=dp) :: a,b,c,d,emc,u,em(13),en(13)
    logical :: bo

    emc=emmc
    u=uu
    if(abs(emc)>epsilon(1.0_dp))then
       bo=(emc.lt.0.0_dp)
       if(bo) then
          d=1.0_dp-emc
          emc=-emc/d
          d=sqrt(d)
          u=d*u
       end if
       a=1.0_dp
       dn=1.0_dp

       do i=1,13
!          write(*,*) i
          l=i
          em(i)=a
          emc=sqrt(emc)
          en(i)=emc
          c=0.5_dp*(a+emc)
          if(abs(a-emc).le.CA*a) exit
          emc=a*emc
          a=c
       end do
       u=c*u
       sn=sin(u)
       cn=cos(u)
       if(abs(sn)>epsilon(1.0_dp)) then
          a=cn/sn
          c=a*c
          do ii=l,1,-1
!             write(*,*) ii
             b=em(ii)
             a=c*a
             c=dn*c
             dn=(en(ii)+a)/(b+a)
             a=c/b
          end do
          a=1.0_dp/sqrt(c**2+1.0_dp)
          if(sn.lt.0.0_dp)then
             sn=-a
          else
             sn=a
          end if
          cn=c*sn
       end if
       if(bo)then
          a=dn
          dn=cn
          cn=a
          sn=sn/d
       end if
    else
       cn=1.0_dp/cosh(u)
       dn=cn
       sn=tanh(u)
    end if
    return
  end subroutine jacobi_elliptic_real

  function incomplete_elliptic_integral(theta,k) result(FF)
    implicit none
    real(kind=dp) :: FF
    real(kind=dp), intent(in) :: theta, k
    FF=sin(theta)*carlson_elliptic_first(cos(theta)**2,1.0_dp-(k**2)*(sin(theta)**2),1.0_dp)
  end function incomplete_elliptic_integral

  function complete_elliptic_integral(k) result(FF)
    implicit none
    real(kind=dp) :: FF
    real(kind=dp), intent(in) :: k
    FF=carlson_elliptic_first(0.0_dp,1.0_dp-(k**2),1.0_dp)
  end function complete_elliptic_integral

  function carlson_elliptic_first(x,y,z) result(rf)
    use utils, only: utils_abort
    implicit none
    real(kind=dp), intent(in) :: x,y,z
    real(kind=dp) :: rf
    real(kind=dp),parameter :: ERRTOL=0.08_dp,TINY=1.5e-38_dp,BIG=3.0e37_dp,THIRD=1.0_dp/3.0_dp,&
         & C1=1.0_dp/24.0_dp,C2=0.1_dp,C3=3.0_dp/44.0_dp,C4=1.0_dp/14.0_dp
    ! Computes Carlson's elliptic integral of the 1st  kind, R_F(x,y,z).
    ! x,y and z must be nonnegative, and at most one can be zero.
    ! TINY must be at least 5 times the machine underflow limit, BIG at most one fifth the machine overflow limit.
    real(kind=dp) ::  alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt

    if(min(x,y,z).lt.0.0_dp.or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,z).gt.BIG) then
       call utils_abort("Error in carlson_elliptic_first : invalid arguments.")
    end if
    xt=x
    yt=y
    zt=z
    do
       sqrtx=sqrt(xt)
       sqrty=sqrt(yt)
       sqrtz=sqrt(zt)
       alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
       xt=0.25_dp*(xt+alamb)
       yt=0.25_dp*(yt+alamb)
       zt=0.25_dp*(zt+alamb)
       ave=THIRD*(xt+yt+zt)
       delx=(ave-xt)/ave
       dely=(ave-yt)/ave
       delz=(ave-zt)/ave
       if(.not.(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)) exit
    end do
    e2=delx*dely-delz**2
    e3=delx*dely*delz
    rf=(1.0_dp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
  end function carlson_elliptic_first


  subroutine matrix_info(mat,name)
    use comms, only: pub_on_root
    use sparse, only: sparse_num_cols, sparse_num_rows
    use sparse_embed, only: sparse_embed_num_cols, sparse_embed_num_rows
    implicit none
    type(smearing_matrix), intent(in) :: mat
    character(len=*) :: name

    integer :: nrows, mcols

    if(pub_on_root) then
       write(stdout,*) "Matrix: ",trim(name)
    end if
    if(mat%matrix_type==matrix_type_standard) then
       if(pub_on_root) then
          write(stdout,*) "STANDARD"
          if(mat%standard_is_cmplx) then
             write(stdout,*) "cmplx"
             if(mat%matrix_is_allocated) then
                write(stdout,*) "size=",size(mat%zdata,1),size(mat%zdata,2)
             else
                write(stdout,*) "not allocated"
             end if
          else
             write(stdout,*) "real"
             if(mat%matrix_is_allocated) then
                write(stdout,*) "size=",size(mat%data,1),size(mat%data,2)
             else
                write(stdout,*) "not allocated"
             end if
          end if

       end if
    elseif(mat%matrix_type==matrix_type_dem) then
       if(pub_on_root) then
          write(stdout,*) "DEM"
          if(mat%datadem%iscmplx) then
             write(stdout,*) "cmplx"
          else
             write(stdout,*) "real"
          end if
          write(stdout,*) "size=", mat%datadem%nrows, mat%datadem%mcols
       end if

    elseif(mat%matrix_type==matrix_type_spam3) then
       if(pub_on_root) then
          write(stdout,*) "SPAM3"
          if(mat%dataspam3%iscmplx) then
             write(stdout,*) "cmplx"
          else
             write(stdout,*) "real"
          end if
       end if
       nrows=sparse_num_rows(mat%dataspam3)
       mcols=sparse_num_cols(mat%dataspam3)
       if(pub_on_root) then
          write(stdout,*) "size=",nrows,mcols
       end if

    elseif(mat%matrix_type==matrix_type_spam3_embed) then
       if(pub_on_root) then
          write(stdout,*) "SPAM3_EMBED:",mat%dataspam3_embed%mrows,mat%dataspam3_embed%ncols
          if(mat%dataspam3_embed%iscmplx) then
             write(stdout,*) "cmplx"
          else
             write(stdout,*) "real"
          end if
       end if
       nrows=sparse_embed_num_rows(mat%dataspam3_embed)
       mcols=sparse_embed_num_cols(mat%dataspam3_embed)
       if(pub_on_root) then
          write(stdout,*) "size=",nrows,mcols
       end if

    end if

  end subroutine matrix_info


    !============================================================================!
  ! These are the internal wrapper routines to the high-level abstract matrix  !
  ! type.                                                                      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   A,B,C        (input)  : Matrices                                         !
  !   alpha,beta   (input)  : Scalars                                          !
  !   p,q          (input)  : Vectors                                          !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Feb 2014.                                        !
  !----------------------------------------------------------------------------!
  ! Put matrix_type case statement inside each to call the appropriate routine !
  ! for this type when routine is called. Feb 2014                             !
  ! This section ends with the next header!                                    !
  !============================================================================!


  ! Finds the dimension of a vector.
  function vector_dimension(v) result(N)
    use sparse, only: sparse_num_rows
    use sparse_embed, only: sparse_embed_num_rows
    use utils, only: utils_abort
    implicit none
    type(vector), intent(in) :: v
    integer :: N

    select case(v%matrix_type)
    case(matrix_type_standard)
       if(v%standard_is_cmplx) then
          N=size(v%zdata)
       else
          N=size(v%data)
       end if
    case(matrix_type_dem)
       N=v%dataDEM%nrows
    case(matrix_type_spam3)
       N = sparse_num_rows(v%dataSPAM3)
    case(matrix_type_spam3_embed)
       N = sparse_embed_num_rows(v%dataSPAM3_EMBED)
    case default
       call utils_abort('vector_dimension is not implemented for this matrix type')
    end select
  end function vector_dimension

  ! Finds the leading dimension of a matrix. They should be square, by the way!
  function matrix_dimension(A,dim) result(N)
    use sparse, only: sparse_num_rows, sparse_num_cols
    use sparse_embed, only: sparse_embed_num_rows, sparse_embed_num_cols
    use utils, only: utils_abort
    implicit none
    type(smearing_matrix), intent(in) :: A
    integer, optional, intent(in) :: dim
    integer :: N
    integer :: locdim
    if(present(dim)) then
       if(dim>2.or.dim<1) then
          call utils_abort("Error in Matrix_dimension : dimension other than 1 or 2 specified!")
       end if
       locdim=dim
    else
       locdim=1
    end if

    select case(A%matrix_type)
    case(matrix_type_standard)
       if(A%standard_is_cmplx) then
          N=size(A%zdata,locdim)
       else
          N=size(A%data,locdim)
       end if
    case(matrix_type_dem)
       if(locdim==1) then
          N=A%dataDEM%nrows
       else if(locdim==2) then
          N=A%dataDEM%mcols
       end if
    case(matrix_type_spam3)
       if(locdim==1) then
          N = sparse_num_rows(A%dataSPAM3)
       else if(locdim==2) then
          N = sparse_num_cols(A%dataSPAM3)
       end if
    case(matrix_type_spam3_embed)
       if(locdim==1) then
          N = sparse_embed_num_rows(A%dataSPAM3_EMBED)
       else if(locdim==2) then
          N = sparse_embed_num_cols(A%dataSPAM3_EMBED)
       end if
    case default
       call utils_abort('matrix_dimension is not implemented for this matrix type')
    end select
  end function matrix_dimension

  ! Computes the entry-wise norm of a matrix, of order ord.
  function matrix_norm(A,ord) result(B)
    use comms, only : comms_reduce

    use sparse, only : sparse_entrywise_norm
    use sparse_embed, only : sparse_embed_entrywise_norm
    use utils, only : utils_abort
    implicit none
    type(smearing_matrix), intent(in) :: A
    integer, intent(in) :: ord
    real(kind=dp) :: B
    real(kind=dp), dimension(:), allocatable :: col
    integer :: i,j,N

    N=matrix_dimension(A)
    B=0.0_dp
    select case(A%matrix_type)
    case(matrix_type_standard)
       ! agrecocmplx
       ! complex case
       if (A%standard_is_cmplx) then
          do i=1,size(A%zdata,2)
             do j=1,size(A%zdata,1)
                B=B+(abs(A%zdata(j,i))**ord)
             end do
          end do
       ! real case
       else
          do i=1,size(A%data,2)
             do j=1,size(A%data,1)
                B=B+(abs(A%data(j,i))**ord)
             end do
          end do
       end if
       !       B=B**(1.0_dp/real(ord,dp))
    case(matrix_type_dem)
       !       allocate(col(N))
       !       do i=1,N
       !          call dense_get_col(col,A%dataDEM,i)
       !          B=B+sum(abs(col)**ord)
       !       end do
       !       deallocate(col)
       if(A%dataDEM%iscmplx) then
          B=sum(abs(A%dataDEM%zmtx)**ord)
       else
          B=sum(abs(A%dataDEM%dmtx)**ord)
       end if
       call comms_reduce("SUM",B)
    case(matrix_type_spam3)
       B=sparse_entrywise_norm(A%dataSPAM3,ord)
    case(matrix_type_spam3_embed)
       B=sparse_embed_entrywise_norm(A%dataSPAM3_EMBED,ord)
    case default
       call utils_abort('matrix_norm is not implemented for this matrix type')
    end select
    ! March 2015 -> moved out of select
    if(A%matrix_type/=matrix_type_spam3) then
       B=B**(1.0_dp/real(ord,dp))
    end if
  end function matrix_norm

  ! Computes the induced 1-norm of a matrix. If you want the infinity
  ! induced norm, then transpose first!
  function matrix_induced_norm(A) result(B)

    use dense, only : dense_norm
    use sparse,only : sparse_1norm
    use sparse_embed,only : sparse_embed_1norm
    use utils, only : utils_abort
    implicit none
    type(smearing_matrix), intent(inout) :: A
    !    integer, intent(in) :: ord
    real(kind=dp) :: B
    real(kind=dp) :: s
    real(kind=dp), dimension(:), allocatable :: col
    integer :: i,j,N

    N=matrix_dimension(A)
    B=0.0_dp
    select case(A%matrix_type)
    case(matrix_type_standard)
       ! agrecocmplx
       ! complex case
       if (A%standard_is_cmplx) then
          do i=1,size(A%zdata,2)
             s=0.0_dp
             do j=1,size(A%zdata,1)
                s=s+abs(A%zdata(j,i))
             end do
             if(s>B) B=s
          end do
       ! real case
       else
          do i=1,size(A%data,2)
             s=0.0_dp
             do j=1,size(A%data,1)
                s=s+abs(A%data(j,i))
             end do
             if(s>B) B=s
          end do
       end if
    case(matrix_type_dem)
       B=dense_norm("1",A%dataDEM)
    case(matrix_type_spam3)
       call sparse_1norm(B,A%dataSPAM3)
    case(matrix_type_spam3_embed)
       ! rc2013: this functionality has not been set up properly with embedding
       if(A%dataSPAM3_EMBED%mrows .gt. 1 .or. A%dataSPAM3_EMBED%ncols .gt. 1) &
            call utils_abort('matrix_induced_norm not &
            &implemented/tested with embedding')
       call sparse_embed_1norm(B,A%dataSPAM3_EMBED)
    case default
       call utils_abort('matrix_induced_norm is not implemented for this matrix type')
    end select

  end function matrix_induced_norm

  ! Allocate a vector of dimension (n).
  subroutine vector_allocate_num(v,n,rand)

    use dense, only : dense_create
    use utils, only : utils_abort, utils_alloc_check, utils_dealloc_check
    implicit none
    type(vector), intent(inout) :: v
    integer, intent(in) :: n
    logical,      optional, intent(in) :: rand
    real(kind=dp), allocatable :: locvec(:)
    integer :: ierr

    select case(v%matrix_type)
    case(matrix_type_standard)
       if(allocated(v%data)) then
          deallocate(v%data,stat=ierr)
          call utils_dealloc_check('vector_allocate_num','v%data',ierr)
       end if
       allocate(v%data(n),stat=ierr)
       call utils_alloc_check('vector_allocate_num','v%data',ierr)
       if(present(rand)) then
          if(rand) then
             call random_number(v%data)
          else
             v%data=0.0_dp
          end if
       else
          v%data=0.0_dp
       end if

    case(matrix_type_dem)
       call dense_create(v%dataDEM,N,1)
       if(present(rand)) then
          if(rand) then
             call random_number(v%dataDEM%dmtx)
             !             allocate(locvec(N))
             !             call random_number(locvec)
             !             call dense_put_col(locvec,v%dataDEM,1)
             !             deallocate(locvec)
          else
             v%dataDEM%dmtx=0.0_dp
          end if
       else
          v%dataDEM%dmtx=0.0_dp
       end if
    case default
       call utils_abort('vector_allocate_num is not implemented for this matrix type')
    end select
  end subroutine vector_allocate_num

  ! Allocate a vector v, based on u, or u'A.
  subroutine vector_allocate_vec(v,u,A,iscmplx,rand)
    use dense, only : dense_create
    use utils, only : utils_abort, utils_alloc_check, utils_dealloc_check
    implicit none
    type(vector), intent(inout) :: v
    type(vector), intent(inout) :: u
    type(smearing_matrix), optional, intent(inout) :: A
    logical,      optional, intent(in) :: iscmplx
    logical,      optional, intent(in) :: rand
    logical :: iscmplx_loc
    integer :: n
    integer :: ierr

    iscmplx_loc=.false.
    if(present(iscmplx)) iscmplx_loc=iscmplx

    v%matrix_type=u%matrix_type
    n=vector_dimension(u)

    select case(v%matrix_type)
    case(matrix_type_standard)
       if(allocated(v%data)) then
          deallocate(v%data,stat=ierr)
          call utils_dealloc_check('vector_allocate_vec','v%data',ierr)
       end if
       if(.not.present(A)) then
          allocate(v%data(vector_dimension(u)),stat=ierr)
          call utils_alloc_check('vector_allocate_vec','v%data',ierr)
       else
          allocate(v%data(matrix_dimension(A,2)),stat=ierr)
          call utils_alloc_check('vector_allocate_vec','v%data',ierr)
       end if
       if(present(rand)) then
          if(rand) then
             call random_number(v%data)
          else
             v%data=0.0_dp
          end if
       else
          v%data=0.0_dp
       end if
    case(matrix_type_dem)
       if(.not.present(A)) then
          !          call dense_create(v%dataDEM,u%dataDEM)
          call dense_create(v%dataDEM,N,1)
       else
          call dense_create(v%dataDEM,matrix_dimension(A,2),1)
       end if

       if(present(rand)) then
          if(rand) then
             call random_number(v%dataDEM%dmtx)
          else
             v%dataDEM%dmtx=0.0_dp
          end if
       else
          v%dataDEM%dmtx=0.0_dp
       end if
    case default
       call utils_abort('vector_allocate_vec is not implemented for this matrix type')
    end select
  end subroutine vector_allocate_vec


  ! Take the vector (2) norm.
  function vector_norm(v) result(norm)
    use dense,only: dense_norm_vec
    use utils, only : utils_abort
    implicit none
    type(vector), intent(inout) :: v
    real(kind=dp) :: norm
    integer :: n

    n=vector_dimension(v)
    select case(v%matrix_type)
    case(matrix_type_standard)
       norm=sqrt(sum(v%data**2))
    case(matrix_type_dem)
       norm=dense_norm_vec(v%dataDEM)
    case default
       call utils_abort('vector_norm is not implemented for this matrix type')
    end select
  end function vector_norm

  ! Scale a vector.
  subroutine vector_scale(v,a)
    use utils, only : utils_abort
    implicit none
    type(vector),  intent(inout) :: v
    real(kind=dp), intent(in)    :: a
    integer                      :: n

#ifdef SCALAPACK
    ! ScaLAPACK subroutine
    external :: pdscal
#endif

    n=vector_dimension(v)
    select case(v%matrix_type)
    case(matrix_type_standard)
       v%data=v%data*a
    case(matrix_type_dem)
       ! Move to dense_mod!!!
#ifdef SCALAPACK
       call pdscal(n, a, v%dataDEM%dmtx, 1, 1, v%dataDEM%blacs_desc, 1)
#else
       v%dataDEM%dmtx=v%dataDEM%dmtx*a
#endif
    case default
       call utils_abort('vector_scale is not implemented for this matrix type')
    end select
  end subroutine vector_scale

  ! Vector axpy. u=u+a*v
  subroutine vector_axpy(u,a,v)

    use dense, only: dense_vec_axpy
    use utils, only : utils_abort
    implicit none
    type(vector),  intent(inout) :: u
    real(kind=dp), intent(in)    :: a
    type(vector),  intent(inout) :: v
    integer                      :: n

    n=vector_dimension(v)
    select case(v%matrix_type)
    case(matrix_type_standard)
       u%data=u%data + a*v%data
    case(matrix_type_dem)
       call dense_vec_axpy(u%dataDEM,a,v%dataDEM)
    case default
       call utils_abort('vector_axpy is not implemented for this matrix type')
    end select
  end subroutine vector_axpy

  ! Vector copy
  subroutine vector_copy(u,v)

    use utils, only : utils_abort
    implicit none

    ! Arguments
    type(vector),  intent(inout) :: u
    type(vector),  intent(inout) :: v

#ifdef SCALAPACK
    ! ScaLAPACK subroutine
    external :: pdcopy
#endif

    ! Local variable
    integer                      :: n

    n=vector_dimension(v)
    select case(v%matrix_type)
    case(matrix_type_standard)
       v%data=u%data
    case(matrix_type_dem)
       ! Move to dense_mod!!!
#ifdef SCALAPACK
       call pdcopy(n, u%dataDEM%dmtx, 1, 1, u%dataDEM%blacs_desc, 1, &
            & v%dataDEM%dmtx, 1, 1, v%dataDEM%blacs_desc, 1)
#else
       v%dataDEM%dmtx=u%dataDEM%dmtx
#endif

    case default
       call utils_abort('vector_copy is not implemented for this matrix type')
    end select
  end subroutine vector_copy


  ! Allocate a matrix of dimension (n,m).
  subroutine matrix_allocate_num(A,n,m)
    use dense, only : dense_create
    use utils, only : utils_abort, utils_alloc_check, utils_dealloc_check
    implicit none
    type(smearing_matrix), intent(inout) :: A
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer :: ierr

    select case(A%matrix_type)
    case(matrix_type_standard)
       ! agrecocmplx
       ! complex case
       if (A%standard_is_cmplx) then
          if(allocated(A%zdata)) then
             deallocate(A%zdata,stat=ierr)
             call utils_dealloc_check('matrix_allocate_num','A%zdata',ierr)
          end if
          allocate(A%zdata(n,m),stat=ierr)
          call utils_alloc_check('matrix_allocate_num','A%zdata',ierr)
          A%zdata=(0.0_DP,0.0_DP)
       ! real case
       else
          if(allocated(A%data)) then
             deallocate(A%data,stat=ierr)
             call utils_dealloc_check('matrix_allocate_num','A%data',ierr)
          end if
          allocate(A%data(n,m),stat=ierr)
          call utils_alloc_check('matrix_allocate_num','A%data',ierr)
          A%data=0.0_dp
       end if
       A%matrix_is_allocated=.true.
    case(matrix_type_dem)
       ! agrecocmplx
       call dense_create(A%dataDEM,n,m,iscmplx=A%standard_is_cmplx)
       A%matrix_is_allocated=.true.
    case default
       call utils_abort('matrix_allocate_num is not implemented for this matrix type')
    end select
  end subroutine matrix_allocate_num

  ! Allocate a matrix A, based on B, or B*C.
  subroutine matrix_allocate_mat(A,B,C,iscmplx)
    use dense, only : dense_create
    use sparse,only : sparse_create
    use sparse_embed, only : sparse_embed_create
    use utils, only : utils_abort, utils_alloc_check, utils_dealloc_check, utils_assert
    implicit none
    type(smearing_matrix), intent(inout) :: A
    type(smearing_matrix), intent(inout) :: B
    type(smearing_matrix), optional, intent(inout) :: C
    logical,      optional, intent(in) :: iscmplx
    logical :: iscmplx_loc
    integer :: ierr

    ! agrecocmplx: change this so that by default A is complex if B is complex,
    ! A is real if B is real
    !iscmplx_loc=.false.
    iscmplx_loc = B%standard_is_cmplx
    if(present(iscmplx)) iscmplx_loc=iscmplx

    A%matrix_type=B%matrix_type
    ! agrecocmplx
    A%standard_is_cmplx = iscmplx_loc

    call utils_assert(.not.A%matrix_is_allocated,"Error in smearing-> matrix_allocate_mat: matrix_already allocated.")

    select case(B%matrix_type)
    case(matrix_type_standard)
       ! agrecocmplx
       ! complex case
       if (iscmplx_loc) then
          if(allocated(A%zdata)) then
             deallocate(A%zdata,stat=ierr)
             call utils_dealloc_check('matrix_allocate_mat','A%zdata',ierr)
          end if
          if(.not.present(C)) then
             allocate(A%zdata(matrix_dimension(B,1),matrix_dimension(B,2)),stat=ierr)
             call utils_alloc_check('matrix_allocate_mat','A%zdata',ierr)
          else
             allocate(A%zdata(matrix_dimension(B,1),matrix_dimension(C,2)),stat=ierr)
             call utils_alloc_check('matrix_allocate_mat','A%zdata',ierr)
          end if
          A%zdata=(0.0_dp,0.0_dp)
       ! real case
       else
          if(allocated(A%data)) then
             deallocate(A%data,stat=ierr)
             call utils_dealloc_check('matrix_allocate_mat','A%data',ierr)
          end if
          if(.not.present(C)) then
             allocate(A%data(matrix_dimension(B,1),matrix_dimension(B,2)),stat=ierr)
             call utils_alloc_check('matrix_allocate_mat','A%data',ierr)
          else
             allocate(A%data(matrix_dimension(B,1),matrix_dimension(C,2)),stat=ierr)
             call utils_alloc_check('matrix_allocate_mat','A%data',ierr)
          end if
          A%data=0.0_dp
       end if
       A%matrix_is_allocated=.true.
    case(matrix_type_dem)
       if(.not.present(C)) then
          ! agrecocmplx: allow for matrix A to have iscmplx different from B
          call dense_create(A%dataDEM,B%dataDEM,iscmplx=iscmplx_loc)
       else
          ! agrecocmplx
          call dense_create(A%dataDEM,matrix_dimension(B,1),matrix_dimension(C,2),iscmplx=iscmplx_loc)
       end if
       A%matrix_is_allocated=.true.
    case(matrix_type_spam3)
       if(.not.present(C)) then
          ! agrecocmplx: allow for matrix A to have iscmplx different from B
          call sparse_create(A%dataSPAM3,B%dataSPAM3,iscmplx=iscmplx_loc)
       else
          ! agrecocmplx
          call sparse_create(A%dataSPAM3,B%dataSPAM3,C%dataSPAM3,iscmplx=iscmplx_loc)
       end if
       A%matrix_is_allocated=.true.
    case(matrix_type_spam3_embed)
       if(.not.present(C)) then
          ! agrecocmplx: allow for matrix A to have iscmplx different from B
          call sparse_embed_create(A%dataSPAM3_EMBED,B%dataSPAM3_EMBED,&
               iscmplx=iscmplx_loc)
       else
          ! agrecocmplx
          call sparse_embed_create(A%dataSPAM3_EMBED,B%dataSPAM3_EMBED,&
               C%dataSPAM3_EMBED,iscmplx=iscmplx_loc)
       end if
       A%matrix_is_allocated=.true.
    case default
       call utils_abort('matrix_allocate_mat is not implemented for this matrix type')
    end select
  end subroutine matrix_allocate_mat

  ! Just like BLAS gemm. C <- beta*C + alpha*A*B. No
  ! built-in transpose support.
  subroutine matrix_multiply(alpha,A,B,beta,C)
    use dense, only : dense_create, dense_product, dense_scale, dense_axpy, dense_destroy
    use sparse,only : sparse_create, sparse_product, sparse_scale, sparse_axpy, sparse_destroy
    use sparse_embed,only : sparse_embed_create, sparse_embed_product, &
         sparse_embed_scale, sparse_embed_axpy, sparse_embed_destroy
    use timer, only : timer_clock
    use utils, only : utils_abort, utils_assert
    implicit none
    real(kind=dp), intent(in) :: alpha
    type(smearing_matrix), intent(inout) :: A !Inout for dense_create_copy (in)
    type(smearing_matrix), intent(in) :: B
    real(kind=dp), intent(in) :: beta
    type(smearing_matrix), intent(inout) :: C
    type(DEM) :: tmpdem
    type(SPAM3) :: tmpspam
    type(SPAM3_EMBED) :: tmpspam_embed

    call timer_clock('matrix_multiply',1)

    ! agrecocmplx: all matrices either complex or real
    call utils_assert((A%standard_is_cmplx .eqv. B%standard_is_cmplx) .and. &
         (C%standard_is_cmplx .eqv. B%standard_is_cmplx), &
         'Error in matrix_multiply: incompatible argument types.')



    select case(A%matrix_type)
    case(matrix_type_standard)
       ! agrecocmplx
       ! complex case
       if (A%standard_is_cmplx) then
          C%zdata=beta*C%zdata + alpha*matmul(A%zdata,B%zdata)
       ! real case
       else
          C%data=beta*C%data + alpha*matmul(A%data,B%data)
       end if
    case(matrix_type_dem)
       call dense_create(tmpdem,C%dataDEM,iscmplx=C%dataDEM%iscmplx)
       call dense_product(tmpdem,A%dataDEM,B%dataDEM)
       call dense_scale(C%dataDEM,beta) !tick

       call dense_axpy(C%dataDEM,tmpdem,alpha)
       call dense_destroy(tmpdem)
    case(matrix_type_spam3)
       if(abs(beta)>epsilon(beta)) then
          call sparse_create(tmpspam,C%dataSPAM3,iscmplx=C%dataSPAM3%iscmplx)
          call sparse_product(tmpspam, A%dataSPAM3, B%dataSPAM3)
          call sparse_scale(C%dataSPAM3,beta)
          call sparse_axpy(C%dataSPAM3,tmpspam,alpha)
          call sparse_destroy(tmpspam)
       else
          call sparse_product(C%dataSPAM3, A%dataSPAM3, B%dataSPAM3)
          call sparse_scale(C%dataSPAM3,alpha)
       end if
    case(matrix_type_spam3_embed)
       if(abs(beta)>epsilon(beta)) then
          call sparse_embed_create(tmpspam_embed,C%dataSPAM3_EMBED,iscmplx=C%dataSPAM3_EMBED%iscmplx)
          call sparse_embed_product(tmpspam_embed, A%dataSPAM3_EMBED, &
               B%dataSPAM3_EMBED)
          call sparse_embed_scale(C%dataSPAM3_EMBED,beta)
          call sparse_embed_axpy(C%dataSPAM3_EMBED,tmpspam_embed,alpha)
          call sparse_embed_destroy(tmpspam_embed)
       else
          call sparse_embed_product(C%dataSPAM3_EMBED, A%dataSPAM3_EMBED, &
               B%dataSPAM3_EMBED)
          call sparse_embed_scale(C%dataSPAM3_EMBED,alpha)
       end if
    case default
       call utils_abort('matrix_multiply is not implemented for this matrix type')
    end select

    call timer_clock('matrix_multiply',2)

  end subroutine matrix_multiply

  ! This routine solves BC=A for C.
  subroutine matrix_solve(A,B,C)


    use sparse, only: sparse_solve2
    use timer, only: timer_clock
    use utils, only: utils_abort
    implicit none
    type(smearing_matrix), intent(inout) :: A
    type(smearing_matrix), intent(in) :: B
    type(smearing_matrix), intent(inout) :: C

    integer :: N

    call timer_clock('matrix_solve',1)

!    if(pub_on_root) then
!       write(stdout,*) "entering matrix_solve"
!    end if


    N=matrix_dimension(A)
    select case(A%matrix_type)
    case(matrix_type_spam3)
!       if(pub_on_root) then
!          write(stdout,*) "calling sparse_solve2"
!       end if

       call sparse_solve2(C%dataSPAM3,A%dataSPAM3,B%dataSPAM3,aplusb_sparsity=.false.)
    case default
       call utils_abort('matrix_solve is not implemented for this matrix type')
    end select

    call timer_clock('matrix_solve',2)

  end subroutine matrix_solve

  ! Matrix-vector multiply. q=Ap.
  subroutine matrix_multiply_vec(A,p,q)
    use dense, only : dense_vector_product

    use timer, only : timer_clock
    use utils, only : utils_abort
    implicit none
    type(smearing_matrix),                intent(inout) :: A !Inout for dense_create_copy (in)
    type(vector),                intent(inout) :: p
    type(vector),                intent(inout) :: q
    !    real(kind=dp), dimension(:), allocatable :: locq
    !    type(dem) :: tmpdem1, tmpdem2, demA_loc
    integer :: N
    integer :: status
    integer, save :: count=0

    call timer_clock('matrix_multiply_vec',1)

    count=count+1

    N=matrix_dimension(A)
    select case(A%matrix_type)
    case(matrix_type_standard)
       q%data=matmul(A%data,p%data)
    case(matrix_type_dem)
       call dense_vector_product(q%dataDEM,A%dataDEM,p%dataDEM)
    case default
       call utils_abort('matrix_multiply_vec is not implemented for this matrix type')
    end select

    call timer_clock('matrix_multiply_vec',2)
  end subroutine matrix_multiply_vec

  ! Find the dot-product of two vectors
  function vector_dot(v,u) result(dot)
    use dense, only: dense_dot_product
    use utils, only: utils_assert, utils_abort
    implicit none
    type(vector), intent(in) :: u,v
    real(kind=dp)            :: dot
    integer                  :: N

    call utils_assert(v%matrix_type==u%matrix_type,&
         &'Error in vector_dot: u and v must be the same type.')
    N=vector_dimension(u)
    call utils_assert(vector_dimension(v)==N,&
         &'Error in vector_dot: u and v must be the same length.')
    select case(v%matrix_type)
    case(matrix_type_standard)
       dot=dot_product(v%data,u%data)
    case(matrix_type_dem)
       dot=dense_dot_product(v%dataDEM,u%dataDEM)
    case default
       call utils_abort('vector_dot is not implemented for this matrix type')
    end select

  end function vector_dot


  ! Multiply a matrix by a scalar and then add a scaled identity
  ! matrix to the result, i.e. A <- (alpha*A) + beta*I.
  ! Scaled identity is optional.
  subroutine matrix_scale(A,alpha,beta)
    use dense, only : dense_scale
    use sparse,only : sparse_scale
    use sparse_embed, only : sparse_embed_scale
    use utils, only : utils_abort
    implicit none
    type(smearing_matrix), intent(inout) :: A
    real(kind=dp), intent(in) :: alpha
    real(kind=dp), optional, intent(in) :: beta
    integer :: i
    select case(A%matrix_type)
    case(matrix_type_standard)
       ! agrecocmplx
       ! complex case
       if (A%standard_is_cmplx) then
          A%zdata=alpha*A%zdata
          if(present(beta)) then
             do i=1,min(size(A%zdata,1),size(A%zdata,2))
                A%zdata(i,i)=A%zdata(i,i)+cmplx(beta,0.0_dp,kind=DP)
             end do
          end if
       ! real case
       else
          A%data=alpha*A%data
          if(present(beta)) then
             do i=1,min(size(A%data,1),size(A%data,2))
                A%data(i,i)=A%data(i,i)+beta
             end do
          end if
       end if
    case(matrix_type_dem)
       if(present(beta)) then
          call dense_scale(A%dataDEM,alpha,beta)
       else
          call dense_scale(A%dataDEM,alpha)
       end if
    case(matrix_type_spam3)
       if(present(beta)) then
          call sparse_scale(A%dataSPAM3,alpha,beta)
       else
          call sparse_scale(A%dataSPAM3,alpha)
       end if
    case(matrix_type_spam3_embed)
       if(present(beta)) then
          call sparse_embed_scale(A%dataSPAM3_EMBED,alpha,beta)
       else
          call sparse_embed_scale(A%dataSPAM3_EMBED,alpha)
       end if
    case default
       call utils_abort('matrix_scale is not implemented for this matrix type')
    end select
  end subroutine matrix_scale


  ! Multiply a matrix by a scalar and then add a scaled identity
  ! matrix to the result, i.e. A <- (alpha*A) + beta*I.
  ! Scaled identity is optional.
  ! Complex version
  subroutine matrix_scale_cmplx(A,alpha,beta)
    use dense, only : dense_scale
    use sparse,only : sparse_scale
    use sparse_embed,only: sparse_embed_scale
    use utils, only : utils_abort
    implicit none
    type(smearing_matrix), intent(inout) :: A
    complex(kind=dp), intent(in) :: alpha
    complex(kind=dp), optional, intent(in) :: beta
    integer :: i
    select case(A%matrix_type)
    case(matrix_type_standard)
       ! agrecocmplx
       ! complex case
       if (A%standard_is_cmplx) then
          A%zdata=alpha*A%zdata
          if(present(beta)) then
             do i=1,min(size(A%zdata,1),size(A%zdata,2))
                A%zdata(i,i)=A%zdata(i,i)+beta
             end do
          end if
       ! real case
       else
          call utils_abort('matrix_scale_cmplx cannot be used with real matrices!')
       end if
    case(matrix_type_dem)
       if(present(beta)) then
          call dense_scale(A%dataDEM,alpha,beta)
       else
          call dense_scale(A%dataDEM,alpha)
       end if
    case(matrix_type_spam3)
       if(present(beta)) then
          call sparse_scale(A%dataSPAM3,alpha,beta)
       else
          call sparse_scale(A%dataSPAM3,alpha)
       end if
    case(matrix_type_spam3_embed)
       if(present(beta)) then
          call sparse_embed_scale(A%dataSPAM3_embed,alpha,beta)
       else
          call sparse_embed_scale(A%dataSPAM3_embed,alpha)
       end if
    case default
       call utils_abort('matrix_scale_cmplx is not implemented for this matrix type')
    end select
  end subroutine matrix_scale_cmplx


  ! Copy a matrix from A to B.
  subroutine matrix_copy(A,B)
    use dense, only : dense_copy
    use sparse,only : sparse_copy
    use sparse_embed, only : sparse_embed_copy
    use utils, only : utils_abort
    implicit none
    type(smearing_matrix), intent(in) :: A
    type(smearing_matrix), intent(inout) :: B


    select case(A%matrix_type)
    case(matrix_type_standard)
!       ! agrecocmplx
!       B%standard_is_cmplx = A%standard_is_cmplx

       ! agrecocmplx
       ! complex case
       if (A%standard_is_cmplx) then
          if(B%standard_is_cmplx) then
             B%zdata = A%zdata
          else
             B%data = real(A%zdata,dp)
          end if
       ! real case
       else
          if(B%standard_is_cmplx) then
             B%zdata = cmplx(A%data,0.0_dp)
          else
             B%data = A%data
          end if
       end if
    case(matrix_type_dem)
       call dense_copy(B%dataDEM,A%dataDEM)
    case(matrix_type_spam3)
       if(.not.B%dataSPAM3%iscmplx.and.A%dataSPAM3%iscmplx) then
          call sparse_copy(B%dataSPAM3,A%dataSPAM3,cmplx_to_real=.true.)
       else
          call sparse_copy(B%dataSPAM3,A%dataSPAM3)
       end if
    case(matrix_type_spam3_embed)
       ! rc2013: SPAM3_EMBED functionality with cmplx not fully implemented
       if(.not.B%dataSPAM3_EMBED%iscmplx.and.A%dataSPAM3_EMBED%iscmplx) then
!         call utils_abort('Complex SPAM3_EMBED not ready yet.')
          call sparse_embed_copy(B%dataSPAM3_EMBED,A%dataSPAM3_EMBED, &
               cmplx_to_real=.true.)
       else
          call sparse_embed_copy(B%dataSPAM3_EMBED,A%dataSPAM3_EMBED)
       end if
    case default
       call utils_abort('matrix_copy is not implemented for this matrix type')
    end select
  end subroutine matrix_copy

  ! Perform a BLAS axpy operation.
  ! A <- A+alpha*B
  subroutine matrix_axpy(A,B,alpha)
    use dense, only : dense_axpy
    use sparse,only : sparse_axpy
    use sparse_embed, only : sparse_embed_axpy
    use utils, only : utils_abort
    implicit none
    type(smearing_matrix), intent(inout) :: A
    type(smearing_matrix), intent(in) :: B
    real(kind=dp), intent(in) :: alpha
    select case(A%matrix_type)
    case(matrix_type_standard)
       ! agrecocmplx
       ! B complex
       if (B%standard_is_cmplx) then
          ! A complex
          if (A%standard_is_cmplx) then
             A%zdata=A%zdata+(alpha*B%zdata)
          else
             ! A real
             A%data=A%data+(alpha*real(B%zdata,kind=DP))
          end if
       ! B real
       else
          ! A complex
          if (A%standard_is_cmplx) then
             A%zdata=A%zdata+(alpha*cmplx(B%data,0.0_DP,kind=DP))
          ! A real
          else
             A%data=A%data+(alpha*B%data)
          end if
       end if
    case(matrix_type_dem)
       call dense_axpy(A%dataDEM,B%dataDEM,alpha)
    case(matrix_type_spam3)
       call sparse_axpy(A%dataSPAM3,B%dataSPAM3,alpha)
    case(matrix_type_spam3_embed)
       call sparse_embed_axpy(A%dataSPAM3_EMBED,B%dataSPAM3_EMBED,alpha)
    case default
       call utils_abort('matrix_axpy is not implemented for this matrix type')
    end select
  end subroutine matrix_axpy

  ! Perform a complex BLAS axpy operation.
  ! A <- A+alpha*B
  subroutine matrix_axpy_cmplx(A,B,alpha)
    use dense, only : dense_axpy
    use sparse,only : sparse_axpy
    use sparse_embed,only: sparse_embed_axpy
    use utils, only : utils_abort
    implicit none
    type(smearing_matrix), intent(inout) :: A
    type(smearing_matrix), intent(in) :: B
    complex(kind=dp), intent(in) :: alpha
    select case(A%matrix_type)
    case(matrix_type_standard)
       ! agrecocmplx
       ! B complex
       if (B%standard_is_cmplx) then
          ! A complex
          if (A%standard_is_cmplx) then
             A%zdata=A%zdata+(alpha*B%zdata)
          else
             call utils_abort('matrix_axpy_cmplx is only available for complex matrices (A real)')
             ! A real
!             A%data=A%data+(alpha*real(B%zdata,kind=DP))
          end if
       ! B real
       else
          ! A complex
          if (A%standard_is_cmplx) then
!             call utils_abort('matrix_axpy_cmplx is only available for complex matrices (B real)')
             A%zdata=A%zdata+(alpha*cmplx(B%data,0.0_DP,kind=DP))
          ! A real
          else
             call utils_abort('matrix_axpy_cmplx is only available for complex matrices (A&B real)')
!             A%data=A%data+(alpha*B%data)
          end if
       end if
    case(matrix_type_dem)
       call dense_axpy(A%dataDEM,B%dataDEM,alpha)
    case(matrix_type_spam3)
       call sparse_axpy(A%dataSPAM3,B%dataSPAM3,alpha)
    case(matrix_type_spam3_embed)
       call sparse_embed_axpy(A%dataSPAM3_embed,B%dataSPAM3_embed,alpha)
    case default
       call utils_abort('matrix_axpy is not implemented for this matrix type')
    end select
  end subroutine matrix_axpy_cmplx



  ! Deallocate vector memory.
  subroutine vector_free(v)
    use dense, only : dense_destroy
    use utils, only : utils_abort, utils_dealloc_check
    implicit none
    type(vector), intent(inout) :: v
    integer :: ierr

    select case(v%matrix_type)
    case(matrix_type_standard)
       if(allocated(v%data)) then
          deallocate(v%data,stat=ierr)
          call utils_dealloc_check('vector_free','v%data',ierr)
       end if
    case(matrix_type_dem)
       call dense_destroy(v%dataDEM)
    case default
       call utils_abort('vector_free is not implemented for this matrix type')
    end select
  end subroutine vector_free

  ! Deallocate matrix memory.
  subroutine matrix_free(A)
    use dense, only : dense_destroy
    use sparse,only : sparse_destroy
    use sparse_embed, only : sparse_embed_destroy
    use utils, only : utils_abort, utils_dealloc_check
    implicit none
    type(smearing_matrix), intent(inout) :: A
    integer :: ierr
    select case(A%matrix_type)
    case(matrix_type_standard)
       ! agrecocmplx
       ! complex case
       if (A%standard_is_cmplx) then
          if(allocated(A%zdata)) then
             deallocate(A%zdata,stat=ierr)
             call utils_dealloc_check('matrix_free','A%zdata',ierr)
          end if
       ! real case
       else
          if(allocated(A%data)) then
             deallocate(A%data,stat=ierr)
             call utils_dealloc_check('matrix_free','A%data',ierr)
          end if
       end if
       A%matrix_is_allocated=.false.
    case(matrix_type_dem)
       call dense_destroy(A%dataDEM)
       A%matrix_is_allocated=.false.
    case(matrix_type_spam3)
       call sparse_destroy(A%dataSPAM3)
       A%matrix_is_allocated=.false.
    case(matrix_type_spam3_embed)
       call sparse_embed_destroy(A%dataSPAM3_EMBED)
       A%matrix_is_allocated=.false.
    case default
       call utils_abort('matrix_free is not implemented for this matrix type')
    end select
  end subroutine matrix_free

  ! Invert matrix. Hopefully we don't need to use this in the long term.
  subroutine matrix_invert(A,Ainv_appr,prec,method,matmuls,error)
    use dense, only : dense_invert
    use timer, only : timer_clock
    use sparse,only : sparse_hotelling_invert, sparse_hotelling_init
    use utils, only : utils_alloc_check, utils_dealloc_check, &
         utils_assert, utils_safe_nint
    implicit none

    ! Arguments
    type(smearing_matrix), intent(inout) :: A
    type(smearing_matrix), optional, intent(inout) :: Ainv_appr
    integer, optional, intent(in) :: prec
    character(len=*), optional, intent(in) :: method
    integer, optional, intent(out) :: matmuls
    integer, optional, intent(out) :: error

    ! LAPACK subroutines
    external :: zgetrf, zgetri, dgetrf, dgetri

    ! Local variables
    integer, dimension(:), allocatable :: pivot
    real(kind=dp), dimension(:), allocatable :: work
    integer :: info
    integer :: work_size
    integer :: D
    type(smearing_matrix) :: loc_invinit
    integer :: loc_matmuls
    logical :: have_invinit
    character :: inversion_method
    integer :: ierr
    ! agrecocmplx
    complex(kind=dp), dimension(:), allocatable :: zwork

    call timer_clock('matrix_invert',1)

    have_invinit=present(Ainv_appr)

    if (have_invinit) then
       ! agrecocmplx: matrices either complex or real
       call utils_assert(A%standard_is_cmplx.eqv.Ainv_appr%standard_is_cmplx, &
            'Error in matrix_invert: incompatible argument types.')
    end if

    if(present(method)) then
       inversion_method=method(1:1)
    else
       inversion_method="I"
    end if

    select case(inversion_method)
    case("d","D")
       select case(A%matrix_type)
       case(matrix_type_standard)
          ! agrecocmplx
          ! complex case
          if (A%standard_is_cmplx) then
             allocate(pivot(min(size(A%zdata,1),size(A%zdata,2))),stat=ierr)
             call utils_alloc_check('matrix_invert','pivot',ierr)
             call zgetrf(size(A%zdata,1),size(A%zdata,2),A%zdata, size(A%zdata,1), pivot, info)
             work_size=-1
             allocate(zwork(1),stat=ierr)
             call utils_alloc_check('matrix_invert','zwork',ierr)
             call zgetri(size(A%zdata,1),A%zdata, size(A%zdata,1), pivot, zwork, work_size, info)
             work_size = utils_safe_nint(real(zwork(1)))
             deallocate(zwork,stat=ierr)
             call utils_dealloc_check('matrix_invert','zwork',ierr)
             allocate(zwork(work_size),stat=ierr)
             call utils_alloc_check('matrix_invert','zwork',ierr)
             call zgetri(size(A%zdata,1),A%zdata, size(A%zdata,1), pivot, zwork, work_size, info)
             deallocate(zwork,stat=ierr)
             call utils_dealloc_check('matrix_invert','zwork',ierr)
             deallocate(pivot,stat=ierr)
             call utils_dealloc_check('matrix_invert','pivot',ierr)
          ! real case
          else
             allocate(pivot(min(size(A%data,1),size(A%data,2))),stat=ierr)
             call utils_alloc_check('matrix_invert','pivot',ierr)
             call dgetrf(size(A%data,1),size(A%data,2),A%data, size(A%data,1), pivot, info)
             work_size=-1
             allocate(work(1),stat=ierr)
             call utils_alloc_check('matrix_invert','work',ierr)
             call dgetri(size(A%data,1),A%data, size(A%data,1), pivot, work, work_size, info)
             work_size = utils_safe_nint(work(1))
             deallocate(work,stat=ierr)
             call utils_dealloc_check('matrix_invert','work',ierr)
             allocate(work(work_size),stat=ierr)
             call utils_alloc_check('matrix_invert','work',ierr)
             call dgetri(size(A%data,1),A%data, size(A%data,1), pivot, work, work_size, info)
             deallocate(work,stat=ierr)
             call utils_dealloc_check('matrix_invert','work',ierr)
             deallocate(pivot,stat=ierr)
             call utils_dealloc_check('matrix_invert','pivot',ierr)
          end if
       case(matrix_type_dem)
          call dense_invert(A%dataDEM,info)
       case(matrix_type_spam3)
          write(stdout,*) "smearing-> Warning: Direct inversion not yet implemented for sparse matrices, using"&
               &"         iterative method instead!"
          D=8
          info=0
          if(present(prec)) D=prec
          call matrix_allocate(loc_invinit,A)
          if(.not.have_invinit) then
             call sparse_hotelling_init(loc_invinit%dataSPAM3,A%dataSPAM3)
          else
             call matrix_copy(Ainv_appr,loc_invinit)
          end if

          !          call nsh_inverse(A,loc_invinit,D,.true.,loc_matmuls)
          call sparse_hotelling_invert(loc_invinit%dataSPAM3,A%dataSPAM3,.false., 10.0_dp**(-D), 100)
          ! Maxiters hard-coded to 100...

          call matrix_copy(loc_invinit,A)
          call matrix_free(loc_invinit)
       case default
       end select
    case("i","I")
       D=8
       info=0
       if(present(prec)) D=prec
       call matrix_allocate(loc_invinit,A)

       if(.not.have_invinit) then
          call nsh_inverse(A,loc_invinit,D,.false.,loc_matmuls)
       else
          call matrix_copy(Ainv_appr,loc_invinit)
          call nsh_inverse(A,loc_invinit,D,.true.,loc_matmuls)
       end if
       call matrix_copy(loc_invinit,A)
       call matrix_free(loc_invinit)
    end select

    if(present(matmuls)) then
       matmuls=loc_matmuls
    end if

    if(present(error)) then
       error=info
    end if

    call timer_clock('matrix_invert',2)

  end subroutine matrix_invert

  ! Transpose A out of place -> B.
  subroutine matrix_transpose(A,B)
    use dense, only : dense_transpose
    use sparse,only : sparse_transpose
    use sparse_embed, only : sparse_embed_transpose
    use utils, only : utils_abort, utils_assert
    implicit none
    type(smearing_matrix), intent(in) :: A
    type(smearing_matrix), intent(inout) :: B

    ! agrecocmplx: either both real or both complex
    call utils_assert(A%standard_is_cmplx .eqv. B%standard_is_cmplx, &
         'Error in matrix_transpose: mixed real/complex argument types.')

    select case(A%matrix_type)
    case(matrix_type_standard)
       ! agrecocmplx
       ! complex case: take the hermitian
       if (A%standard_is_cmplx) then
          B%zdata=transpose(conjg(A%zdata))
       ! real case
       else
          B%data=transpose(A%data)
       end if
    case(matrix_type_dem)
       call dense_transpose(B%dataDEM,A%dataDEM)
    case(matrix_type_spam3)
       call sparse_transpose(B%dataSPAM3,A%dataSPAM3)
    case(matrix_type_spam3_embed)
       call sparse_embed_transpose(B%dataSPAM3_EMBED,A%dataSPAM3_EMBED)
    case default
       call utils_abort('matrix_transpose is not implemented for this matrix type')
    end select
  end subroutine matrix_transpose

  ! Work out the matrix trace of A.
  function matrix_old_trace(A) result(trace)
    use utils, only : utils_abort
    implicit none
    type(smearing_matrix), intent(in) :: A
    real(kind=dp) :: trace
    integer :: i
    trace=0.0_dp
    select case(A%matrix_type)
    case(matrix_type_standard)
       ! agrecocmplx
       ! complex case
       if (A%standard_is_cmplx) then
          do i=1,size(A%zdata,1)
             trace=trace+real(A%zdata(i,i),kind=DP)
          end do
       ! real case
       else
          do i=1,size(A%data,1)
             trace=trace+A%data(i,i)
          end do
       end if
    case default
       call utils_abort('matrix_old_trace is not implemented for this matrix type')
    end select
  end function matrix_old_trace

  ! Work out the matrix trace of A, or product AB
  function matrix_trace(A,B) result(trace)

    use dense, only : dense_trace
    use sparse,only : sparse_trace
    use sparse_embed, only : sparse_embed_trace
    use utils, only : utils_abort, utils_assert
    implicit none
    type(smearing_matrix), intent(inout) :: A
    type(smearing_matrix), optional, intent(inout) :: B
    type(smearing_matrix) :: tmp_mat, tmp_mat2
    real(kind=dp) :: trace
    real(kind=dp) :: tmp
    real(kind=dp), dimension(:), allocatable :: col
    integer :: i,j,N,M
    integer :: stat
    trace=0.0_dp
    if(present(B)) then
       ! agrecocmplx: all matrices either complex or real
       call utils_assert(A%standard_is_cmplx .eqv. B%standard_is_cmplx, &
            'Error in matrix_trace: incompatible argument types.')
       N=min(matrix_dimension(A,1),matrix_dimension(B,2))
       M=matrix_dimension(A,2)
       if(M/=matrix_dimension(B,1)) then
          call utils_abort("Error in matrix_trace: dim 2 of A and dim 1 of B different!")
       end if
    else
       N=min(matrix_dimension(A,1),matrix_dimension(A,2))
    end if
    select case(A%matrix_type)
    case(matrix_type_standard)
       if(.not.present(B)) then
          ! agrecocmplx
          ! complex case
          if (A%standard_is_cmplx) then
             do i=1,N
                trace=trace+real(A%zdata(i,i),kind=DP)
             end do
          ! real case
          else
             do i=1,N
                trace=trace+A%data(i,i)
             end do
          end if
       else
          ! agrecocmplx
          ! complex case
          if (A%standard_is_cmplx) then
             do i=1,N
                do j=1,M
                   trace=trace+real(A%zdata(i,j)*B%zdata(j,i),kind=DP)
                end do
             end do
          ! real case
          else
             do i=1,N
                trace=trace+dot_product(A%data(i,:),B%data(:,i))
             end do
          end if
       end if
    case(matrix_type_dem)
       if(.not.present(B)) then
          trace = dense_trace(A%dataDEM)
       else
          trace = dense_trace(A%dataDEM,B%dataDEM)
       end if
    case(matrix_type_spam3)
       if(.not.present(B)) then
          trace = sparse_trace(A%dataSPAM3)
       else
          trace = sparse_trace(A%dataSPAM3,B%dataSPAM3)
       end if
    case(matrix_type_spam3_embed)
       if(.not.present(B)) then
          call sparse_embed_trace(trace,A%dataSPAM3_EMBED)
       else
          call sparse_embed_trace(trace,A%dataSPAM3_EMBED,B%dataSPAM3_EMBED)
       end if
    case default
       call utils_abort('matrix_trace is not implemented for this matrix type')
    end select
  end function matrix_trace

  ! Debugging routine for writing a (standard matrix)
  ! to screen. This isn't called in the production version.
  ! Consider deleting it?
  subroutine matrix_write(A)
    use dense, only : dense_get_element
    use utils, only : utils_abort, utils_alloc_check, utils_dealloc_check
    implicit none
    type(smearing_matrix), intent(in) :: A
    integer :: i,j,N
    real(kind=dp),dimension(:),allocatable :: vec
    integer :: ierr

    N=matrix_dimension(A)
    write(stdout,*) "smearing-> matrix ---------------------------------------------------------------"
    select case(A%matrix_type)
    case(matrix_type_standard)
       do i=1,N
          write(stdout,*) A%data(i,:)
       end do
    case(matrix_type_dem)
       allocate(vec(N),stat=ierr)
       call utils_alloc_check('matrix_write','vec',ierr)
       do i=1,N
          do j=1,N
             call dense_get_element(vec(j),A%dataDEM,i,j)
          end do
          write(stdout,*) vec
       end do
       deallocate(vec,stat=ierr)
       call utils_dealloc_check('matrix_write','vec',ierr)
    case default
       call utils_abort('matrix_write is not implemented for this matrix type')
    end select
    write(stdout,*) "smearing-> matrix ---------------------------------------------------------------"
  end subroutine matrix_write

  ! Check whether any element of a matrix is NaN.
  function matrix_any_isnan(A) result(isnan)
    use comms,  only: comms_reduce
    use sparse, only: sparse_any_isnan
    use sparse_embed, only: sparse_embed_any_isnan
    use utils,  only: utils_isnan, utils_abort
    implicit none
    type(smearing_matrix), intent(in) :: A
    logical :: isnan

    logical :: loc_isnan
    integer :: i,j

    loc_isnan=.false.

    select case(A%matrix_type)
    case(matrix_type_standard)
       if(A%standard_is_cmplx) then
          do i=1,matrix_dimension(A,2)
             do j=1,matrix_dimension(A,1)
                loc_isnan=loc_isnan.or.utils_isnan(real(A%zdata(j,i),dp))
                loc_isnan=loc_isnan.or.utils_isnan(aimag(A%zdata(j,i)))
             end do
          end do
       else
          do i=1,matrix_dimension(A,2)
             do j=1,matrix_dimension(A,1)
                loc_isnan=loc_isnan.or.utils_isnan(A%data(j,i))
             end do
          end do
       end if
    case(matrix_type_dem)
       if(A%dataDEM%iscmplx) then
          do i=1,size(A%dataDEM%dmtx,2)
             do j=1,size(A%dataDEM%dmtx,1)
                loc_isnan=loc_isnan.or.utils_isnan(real(A%dataDEM%zmtx(j,i),dp))
                loc_isnan=loc_isnan.or.utils_isnan(aimag(A%dataDEM%zmtx(j,i)))
             end do
          end do
       else
          do i=1,size(A%dataDEM%dmtx,2)
             do j=1,size(A%dataDEM%dmtx,1)
                loc_isnan=loc_isnan.or.utils_isnan(A%dataDEM%dmtx(j,i))
             end do
          end do
       end if
    case(matrix_type_spam3)
       loc_isnan=sparse_any_isnan(A%dataSPAM3)
    case(matrix_type_spam3_embed)
       loc_isnan=sparse_embed_any_isnan(A%dataSPAM3_embed)
    case default
       call utils_abort('matrix_any_isnan is not implemented for this matrix type')
    end select

    call comms_reduce('OR',loc_isnan)

    isnan=loc_isnan

  end function matrix_any_isnan

  !============================================================================!
  !                     End of matrix wrapper routines!                        !
  !============================================================================!


end module smearing_operator
