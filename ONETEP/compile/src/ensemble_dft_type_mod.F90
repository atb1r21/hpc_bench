!================================================================!
!                                                                !
!                    Ensemble DFT model container                !
!                                                                !
!----------------------------------------------------------------!
! Defines a container for ensemble DFT-associated information.   !
! This makes it easier to pass EDFT data across the call chain.  !
! The type itself cannot reside in ensemble_dft_mod.F90, because !
! it is needed in hamiltonian_mod, and hamiltonian_mod cannot    !
! depend on ensemble_dft_mod, because it is called from there.   !
!----------------------------------------------------------------!
! Consider moving edft_create() and edft_destroy() here too.     !
!----------------------------------------------------------------!
! Split off from ensemble_dft_mod by Jacek Dziedzic in June 2022.!
!================================================================!

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

module ensemble_dft_type

  use constants, only: DP
  use dense, only: DEM
  use sparse_embed, only: SPAM3_EMBED

  implicit none

  private

  public :: EDFT_MODEL

  type EDFT_MODEL

     ! ars: flag to indicate allocation
     logical :: allocd

     ! ars: flag to indicate initialisation
     logical :: initialised

     ! ars: number of molecular orbitals in the model
     integer :: num

     ! ars: number of MOs with non-zero occupancy
     integer :: nbands

     ! kkbd: total number of electrons
     real(kind=DP) :: nelec

     ! ars: fractional occupancies
     real(kind=DP), allocatable :: occ(:,:)

     ! ars: eigenvalues of the Hamiltonian
     real(kind=DP), allocatable :: h_evals(:,:)

     ! ars: Fermi level
     real(kind=DP), allocatable :: fermi(:)

     ! ars: Integrated number of electrons per spin channel
     real(kind=DP), allocatable :: integrated_ne(:)

     ! kkbd: Spin-occupancy of each spin channel. This is the target for fixed-
     !       spin runs, essentially a real equivalent of n_occ
     real(kind=DP), allocatable :: s_occ(:)

     ! ars: Residual contribution to the Lagrangian due to non-orthogonality
     real(kind=DP), allocatable :: orth_residue(:)

     ! ars: rotation from NGWFs to molecular orbitals basis set (M^\alpha_i)
     type(DEM), allocatable :: mo(:)

     ! ja531 --> NGWF update transformation matrix.
     type(SPAM3_EMBED) :: rot

     ! ja531 --> backup Hamiltonian for FOE rotation
     type(SPAM3_EMBED), allocatable :: old_ham(:)

     ! gab: history of ham matrices
     type(SPAM3_EMBED), allocatable :: ham_in_his(:,:)

     ! gab: history of residual ham matrices
     type(SPAM3_EMBED), allocatable :: ham_err_his(:,:)

     ! ars: and last but not least, the free energy, the energy and the entropy
     real(kind=DP) :: free_energy, energy, entropy

  end type EDFT_MODEL

end module ensemble_dft_type
