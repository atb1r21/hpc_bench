!================================================================!
!                                                                !
!               Density Functional Tight Binding module          !
!                                                                !
! This module contains subroutines for performing DFTB calcula-  !
! tions. So far, we are focusing on xTB methods by the Grimme    !
! group.                                                         !
!----------------------------------------------------------------!
! Work started in 2021/06. Written by Jacek Dziedzic and Arihant !
! Bhandari under the supervision of Chris-Kriton Skylaris.       !
!================================================================!
! References:                                                    !
! [1] P. Pracht, E. Caldeweyher, S. Ehlert, and S. Grimme,       !
!     A robust non-self-consistent tight-binding quantum         !
!     chemistry method for large molecules.                      !
!     ChemRxiv. 2019: 1-19.                                      !
!     https://doi.org/10.26434/chemrxiv.8326202.                 !
! [2] S. Grimme, J. Antony, S. Ehrlich, and H. Krieg,            !
!     A consistent and accurate ab initio parametrization of     !
!     density functional dispersion correction (DFT-D) for the   !
!     94 elements H-Pu.                                          !
!     J. Chem. Phys. 2010 (132), 154104.                         !
! [3] C. Bannwarth, S. Ehlert, and S. Grimme,                    !
!     GFN2-xTB - An Accurate and Broadly Parametrized            !
!     Self-Consistent Tight-Binding Quantum Chemical Method      !
!     with Multipole Electrostatics and Density-Dependent        !
!     Dispersion Contributions,                                  !
!     J. Chem. Theory Comput. 2019 (15), 1652.                   !
!     https://doi.org/10.1021/acs.jctc.8b01176.                  !
! [4] S. Grimme, C. Bannwarth, and P. Shushkov,                  !
!     A Robust and Accurate Tight-Binding Quantum Chemical       !
!     Method for Structures, Vibrational Frequencies, and        !
!     Noncovalent Interactions of Large Molecular Systems        !
!     Parametrized for All Spd-Block Elements (Z = 1-86),        !
!     J. Chem. Theory Comput. 2017 (13), 1989.                   !
!     https://doi.org/10.1021/acs.jctc.7b00118.                  !
! [5] C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen,        !
!     P. Pracht, J. Seibert, S. Spicher, and S. Grimme,          !
!     Extended tight-binding quantum chemistry methods,          !
!     Wiley Interdiscip. Rev. Comput. Mol. Sci. 2021 (11), 1     !
!     https://doi.org/10.1002/wcms.1493.                         !
! [6] The Computational Modelling of Heavy Atom Chemistry,       !
!     Chris-Kriton Skylaris (PhD thesis), Cambridge, 1999.       !
!================================================================!


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

module dftb

  use constants, only: stdout, DP, garbage_int, garbage_real
  use dense, only: DEM
  use geometry, only: POINT
  use neighbour_list, only: NL_NEIGHBOUR_LIST
  use rundat, only: pub_debug_on_root

  implicit none

  private

  public :: dftb_init_stage_1
  public :: dftb_init_stage_2
  public :: dftb_exit
  public :: dftb_free
  public :: dftb_rep
  public :: dftb_srb
  public :: dftb_ies
  public :: dftb_calculate_ham
  public :: dftb_eht_energy
  public :: dftb_electronic_force
  public :: dftb_analytical_gaussian_overlap
  public :: dftb_print_energy_components
  public :: dftb_print_force
  public :: dftb_qc_print

  abstract interface
  ! ab: so that we can pass a function as argument to another function
  pure function functions(dr,rc)

    use constants, only: DP
    use geometry, only: POINT

    implicit none

    ! ab: arguments
    type(POINT), intent(in)   :: dr ! ab: a vector
    real(kind=DP), intent(in) :: rc ! ab: some parameter

    ! ab: local variable
    real(kind=DP) :: functions(1:3)

  end function functions
  end interface
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  ! ---------------------------------------------------------------------------
  ! ----- P U B L I C   S U B R O U T I N E S   A N D   F U N C T I O N S -----
  ! ---------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_init_stage_1(mdl)
    !==========================================================================!
    ! This subroutine performs the first stage of DFTB initialisation.         !
    !   - Parameter files are read.                                            !
    !   - The number of NGWFs is inferred from the DFTB AO string for each     !
    !     ONETEP species. The number of NGWFs is then corrected in the species !
    !     and elements arrays on all MPI ranks.                                !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   mdl (inout): %dftb_par is modified, as well as the species and elements!
    !                arrays.                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 2021/06/03.                                 !
    !==========================================================================!

    use model_type, only: MODEL
    use rundat, only: pub_dftb_method_param_file, pub_dftb_common_param_file, &
         pub_dftb
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(MODEL), intent(inout) :: mdl

    ! jd: Local variables
    character(len=*), parameter :: myself = 'dftb_init_stage_1'

    ! -------------------------------------------------------------------------

    call utils_assert(pub_dftb, &
         myself//': Must only be invoked if DFTB is in use.')

    call dftb_read_params(mdl, pub_dftb_method_param_file, &
         pub_dftb_common_param_file)

    call dftb_adjust_num_ngwfs(mdl)

    ! ab: create inferred parameters from read parameters
    call dftb_infer_params(mdl)

  end subroutine dftb_init_stage_1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_init_stage_2(mdl, rep, ngwf_basis, ham)
    !==========================================================================!
    ! This subroutine initialises the DFTB module.                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 2021/06/03.                                 !
    !==========================================================================!

    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use neighbour_list, only: neighbour_list_init_from_sparse
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_dftb
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(MODEL), intent(inout)     :: mdl
    type(NGWF_REP), intent(inout)  :: rep
    type(FUNC_BASIS), intent(in)   :: ngwf_basis(rep%nsub)
    type(NGWF_HAM), intent(inout)  :: ham

    ! jd: Local variables
    character(len=*), parameter :: myself = 'dftb_init_stage_2'
    ! -------------------------------------------------------------------------

    call utils_assert(pub_dftb, &
         myself//': Must only be invoked if DFTB is in use')

    call neighbour_list_init_from_sparse(mdl%dftb_var%nl, 'dftb_nl', rep%overlap%m(1,1))

    call dftb_calculate_overlap(mdl, rep, ngwf_basis)

    ! ab: initialize coordination numbers
    call dftb_init_coordination(mdl)

   end subroutine dftb_init_stage_2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_free(mdl)
    !==========================================================================!
    ! Cleans up after the DFTB module.                                         !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 2021/06/03.                                 !
    !==========================================================================!

    use dense, only: dense_destroy
    use neighbour_list, only: neighbour_list_free
    use model_type, only: MODEL
    use rundat, only: pub_dftb_calculate_force

    implicit none

    ! jd: Arguments
    type(MODEL), intent(inout) :: mdl

    ! jd: Local variables
    character(len=*), parameter :: myself = 'dftb_free'
    integer :: ic

    ! -------------------------------------------------------------------------

    ! ab: deallocate matrices
    if (pub_dftb_calculate_force) then
       do ic = 1, 3
          call dense_destroy(mdl%dftb_var%dqdr(ic))
          call dense_destroy(mdl%dftb_var%dcdr(ic))
          call dense_destroy(mdl%dftb_var%dxdr(ic))
       end do
    end if

    call neighbour_list_free(mdl%dftb_var%nl)

  end subroutine dftb_free

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_exit(mdl)
    !==========================================================================!
    ! Cleans up after the DFTB module.                                         !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari on 2022/06/30.                               !
    !==========================================================================!

    use dftb_parameters, only: dftb_num_species
    use model_type, only: MODEL
    use rundat, only: pub_dftb
    use utils, only: utils_assert, utils_dealloc_check

    implicit none

    ! jd: Arguments
    type(MODEL), intent(inout) :: mdl

    ! jd: Local variables
    character(len=*), parameter :: myself = 'dftb_free'
    integer :: ierr, iat, isp, zsp, ishell, nshells

    ! -------------------------------------------------------------------------

    call utils_assert(pub_dftb, &
         myself//': Must only be invoked if DFTB is in use')

    do iat=1,dftb_num_species
       nshells = mdl%dftb_par%z(iat)%nshells
       do ishell = 1,nshells
          deallocate(mdl%dftb_par%z(iat)%sto(ishell)%alpha, stat=ierr)
          call utils_dealloc_check(myself, &
               'mdl%dftb_par%z(iat)%sto(ishell)%alpha', ierr)
          deallocate(mdl%dftb_par%z(iat)%sto(ishell)%coeff, stat=ierr)
          call utils_dealloc_check(myself, &
               'mdl%dftb_par%z(iat)%sto(ishell)%coeff', ierr)
       end do
       deallocate(mdl%dftb_par%z(iat)%sto, stat=ierr)
       call utils_dealloc_check(myself, &
            'mdl%dftb_par%z(iat)%sto', ierr)
    end do

    do isp = 1, mdl%regions(1)%par%num_species
       zsp = mdl%regions(1)%species(isp)%atomic_number
       deallocate(mdl%dftb_par%z(zsp)%ang_mom, stat=ierr)
       call utils_dealloc_check(myself, &
            'mdl%dftb_par%z(zsp)%ang_mom', ierr)
       deallocate(mdl%dftb_par%z(zsp)%ang, stat=ierr)
       call utils_dealloc_check(myself, &
            'mdl%dftb_par%z(zsp)%ang', ierr)
       deallocate(mdl%dftb_par%z(zsp)%ngwf_type, stat=ierr)
       call utils_dealloc_check(myself, &
            'mdl%dftb_par%z(zsp)%ngwf_type', ierr)
    end do

  end subroutine dftb_exit

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_rep(dftb_par, elements, par, cell, e_rep, f_rep)
    !==========================================================================!
    ! Calculates the repulsive energy term, E_rep, cf. [1]:(3).                !
    ! Hybrid parallelism is in effect, the calculated value is the global      !
    ! energy.                                                                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   dftb_par (in): DFTB model parameters.                                  !
    !   elements (in): Needed for ionic positions.                             !
    !   par (in): Parallel distribution info.                                  !
    !   cell (in): Needed for periodicity, if in effect.                       !
    !   e_rep (out): Total repulsion energy.                                   !
    !   f_rep (out): Total repulsion force.                                    !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 2021/06/15.                                 !
    ! Corrected by Arihant Bhandari in Feb 2022.                               !
    ! Forces added by Arihant Bhandari in Mar 2022.                            !
    !==========================================================================!

    use comms, only: pub_my_proc_id, comms_reduce
    use dftb_parameters, only: DFTB_PARAMS, dftb_max_species
    use geometry, only: POINT, OPERATOR(.DOT.), OPERATOR(-), OPERATOR(+)
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_threads_max, pub_dftb_calculate_force, pub_print_qc, &
         pub_dftb_rep_cutoff, pub_dftb_bc_is_periodic
    use simulation_cell, only: CELL_INFO, &
         simulation_cell_lattice_points_init,simulation_cell_lattice_points_exit
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_abort, utils_int_to_str

    implicit none

    ! jd: Arguments
    type(DFTB_PARAMS), intent(in) :: dftb_par
    type(ELEMENT), intent(in)     :: elements(:)
    type(PARAL_INFO), intent(in)  :: par
    type(CELL_INFO), intent(in)   :: cell
    real(kind=DP), intent(out)    :: e_rep
    real(kind=DP), intent(inout)  :: f_rep(1:3,1:par%nat)

    ! jd: Local variables
    integer :: local_iat, itr
    integer :: global_iat, global_jat
    integer :: orig_iat, orig_jat
    integer :: zi, zj ! jd: atomic numbers of atoms I and J
    real(kind=DP) :: z_eff_i, z_eff_j ! jd: Effective core charges of atoms I, J
    real(kind=DP) :: alpha_i, alpha_j ! jd: Element-specific parameters
    real(kind=DP) :: eni, enj  ! ab: element specific electronegativity
    real(kind=DP) :: R1, R2    ! jd: distance between atoms I, J
    real(kind=DP) :: prefactor ! jd: Fraction in [1]:(3)
    real(kind=DP) :: expterm   ! jd: Exponential term in [1]:(3)
    real(kind=DP) :: d2, d4, rs, alpha, cutoff
    type(POINT)   :: R_i, R_j  ! jd: Positions of atoms I, J
    type(POINT)   :: R_ij
    real(kind=DP) :: Rij(1:3)
    real(kind=DP) :: dedr, expo
    type(POINT)   :: lattice(3) ! ab: lattice vectors
    type(POINT), allocatable :: points(:) ! ab: translational lattice points
    real(kind=DP), parameter :: lcut = 1D-6 ! ab: for excluding I==J
    character(len=*), parameter :: myself = 'dftb_rep'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    call utils_assert(par%nat == size(elements), myself//': &
         &allocated parallel strategy is incompatible with elements.')

    ! ab: choose cutoff for periodic translations as per boundary conditions
    if(.not. any(pub_dftb_bc_is_periodic)) then
       cutoff = 0.0_DP
    else if(all(pub_dftb_bc_is_periodic)) then
       cutoff = pub_dftb_rep_cutoff
    else
       call utils_abort(myself//': Mixed BCs are not currently supported.')
    end if

    ! ab: lattice vectors:
    lattice(1) = cell%a1
    lattice(2) = cell%a2
    lattice(3) = cell%a3

    ! ab: generate lattice points
    call simulation_cell_lattice_points_init(lattice, cutoff, points)

    ! ab: initialize
    e_rep = 0.0_DP
    f_rep = 0.0_DP

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(STATIC)  &
!$OMP PRIVATE(local_iat,global_iat,zi,R_i,z_eff_i,alpha_i,itr,eni,enj,d2,d4,rs,&
!$OMP      global_jat,zj,R_j,z_eff_j,alpha_j,R_ij,Rij,R1,R2,prefactor,alpha,   &
!$OMP      expterm,orig_iat,orig_jat,dedr,expo)                                &
!$OMP SHARED(par,elements,dftb_par,cell,pub_my_proc_id,points,&
!$OMP      pub_dftb_calculate_force,pub_dftb_rep_cutoff) &
!$OMP REDUCTION(+:e_rep,f_rep)
    ! jd: This rank's atoms I
    loop_atom_I:                                                               &
    do local_iat = 1, par%num_atoms_on_proc(pub_my_proc_id)
       global_iat = par%first_atom_on_proc(pub_my_proc_id) + local_iat-1
       orig_iat = par%orig_atom(global_iat)
       zi = elements(orig_iat)%atomic_number
       R_i = elements(orig_iat)%centre
       if(zi > dftb_max_species .or. zi < 1) then
          call utils_abort(myself//': atomic number zi out of range: '//&
               trim(utils_int_to_str(zi))//'.', avoid_mpi_calls = .true.)
       end if
       alpha_i = dftb_par%z(zi)%repa
       z_eff_i = dftb_par%z(zi)%repb
       if (dftb_par%z(zi)%en==garbage_real) then
          eni = dftb_par%z(zi)%paulingen
       else
          eni = dftb_par%z(zi)%en
       end if
       ! jd: All atoms J < I
       loop_atom_J:                                                            &
       do global_jat = 1, global_iat
          orig_jat = par%orig_atom(global_jat)
          zj = elements(orig_jat)%atomic_number
          R_j = elements(orig_jat)%centre
          if (dftb_par%z(zj)%en==garbage_real) then
             enj = dftb_par%z(zj)%paulingen
          else
             enj = dftb_par%z(zj)%en
          end if
          if(zj > dftb_max_species .or. zj < 1) then
             call utils_abort(myself//': atomic number zj out of range: '//&
                  trim(utils_int_to_str(zj))//'.', avoid_mpi_calls = .true.)
          end if
          alpha_j = dftb_par%z(zj)%repa
          z_eff_j = dftb_par%z(zj)%repb
          ! ab: loop over all translational lattice points
          do itr = 1, size(points)
             R_ij = R_I - R_J + points(itr)
             R2 = R_ij .DOT. R_ij

             if (R2 > pub_dftb_rep_cutoff**2 .or. R2 < lcut ) cycle

             Rij = (/ R_ij%X, R_ij%Y, R_ij%Z /)
             R1 = sqrt(R2)
             d2 = (eni-enj)**2
             d4 = d2**2
             prefactor = z_eff_i * z_eff_j
             rs = dftb_par%renscale
             alpha = (1.0_DP + rs * (d2+d4) / 100.0_DP) * sqrt(alpha_i*alpha_j)
             expo = alpha * R1**1.5_DP
             expterm = prefactor * exp( -expo )
             e_rep = e_rep + expterm / R1
             if (pub_dftb_calculate_force) then
                dedr = expterm * (1.0_DP + 1.5_DP * expo) / R1**3
                f_rep(:,orig_iat) = f_rep(:,orig_iat) + dedr * Rij(:)
                f_rep(:,orig_jat) = f_rep(:,orig_jat) - dedr * Rij(:)
             end if
          end do
       end do loop_atom_J
    end do loop_atom_I
!$OMP END PARALLEL DO

    call comms_reduce('SUM', e_rep)
    if (pub_dftb_calculate_force) then
       call comms_reduce('SUM', f_rep)
    end if

    call simulation_cell_lattice_points_exit(points)

    if (pub_print_qc) then
       call dftb_qc_print(myself,e_rep,size(elements),f_rep)
    end if

    call timer_clock(myself,2)

  end subroutine dftb_rep

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_srb(dftb_par, elements, par, cell, dcdr, e_srb, f_srb)
    !==========================================================================!
    ! Calculates the short-range basis/bond correction term E_srb, [1]:(4).    !
    ! Hybrid parallelism is in effect, the calculated value is the global      !
    ! energy.                                                                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   dftb_par (in): DFTB model parameters.                                  !
    !   elements (in): Needed for ionic positions.                             !
    !   par (in): Parallel distribution info.                                  !
    !   cell (in): Needed for periodicity, if in effect.                       !
    !   dcdr(in): Coordination number derivative for calculating forces.       !
    !   e_srb (out): Total srb energy.                                         !
    !   f_srb (out): Total srb force.                                          !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 2021/06/16.                                 !
    ! Corrected by Arihant Bhandari in Feb 2022.                               !
    ! Forces added by Arihant Bhandari in Mar 2022.                            !
    !==========================================================================!

    use comms, only: pub_my_proc_id, comms_reduce
    use constants, only: DFTB_GFN0
    use dftb_parameters, only: DFTB_PARAMS, dftb_max_species
    use dense, only: DEM, dense_create, dense_product, dense_destroy, &
         dense_get_col, dense_put_col
    use geometry, only: POINT, magnitude, OPERATOR(-), OPERATOR(+)
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_threads_max, pub_dftb_method, pub_dftb_calculate_force, &
         pub_print_qc, pub_dftb_srb_cutoff, pub_dftb_bc_is_periodic
    use simulation_cell, only: CELL_INFO, minimum_image_distance, &
         simulation_cell_lattice_points_init,simulation_cell_lattice_points_exit
    use utils, only: utils_assert, utils_abort, utils_int_to_str

    implicit none

    ! jd: Arguments
    type(DFTB_PARAMS), intent(in) :: dftb_par
    type(ELEMENT), intent(in)     :: elements(:)
    type(PARAL_INFO), intent(in)  :: par
    type(CELL_INFO), intent(in)   :: cell
    type(DEM), intent(in)         :: dcdr(3)
    real(kind=DP), intent(out)    :: e_srb
    real(kind=DP), intent(out)    :: f_srb(1:3,1:par%nat)

    ! jd: Local variables
    integer :: local_iat, itr, ntr, ic
    integer :: global_iat, global_jat
    integer :: orig_iat, orig_jat
    integer :: zi, zj       ! jd: Atomic numbers of atoms I and J
    type(POINT) :: R_i, R_j ! jd: Positions of atoms I, J
    real(kind=DP) :: term1  ! jd: -\eta_srb * (1+g_scal*(\Delta EN)^2) in [1]:(4)
    real(kind=DP) :: dterm, expterm, dr0dci, dr0dcj, dfdc(par%nat+1)
    real(kind=DP) :: r_ij   ! jd: minimum image distance between atoms I, J
    real(kind=DP) :: lat_cut! ab: cutoff for lattice translations
    type(POINT) :: rvec     ! ab: translational posiion vector
    real(kind=DP) :: rv(1:3)
    real(kind=DP) :: rij    ! ab: translational distance between atoms I, J
    real(kind=DP) :: delta_r! jd: r_ij - r0_ij
    real(kind=DP) :: r0_ij  ! jd: CN-scaled damping distance, cf. [1]:(5)
    real(kind=DP) :: delta_en ! jd: Difference of electronegativities
    type(POINT)   :: lattice(3) ! ab: lattice vectors
    integer, parameter :: z_min = 5 ! jd: Minimum z for which SRB is applied
    integer, parameter :: z_max = 9 ! jd: Maximum z for which SRB is applied
    type(POINT), allocatable :: points(:) ! ab: translational lattice points
    type(DEM) :: dfc, fc
    character(len=*), parameter :: myself = 'dftb_srb'

    ! -------------------------------------------------------------------------

    ! ab: initialize
    e_srb = 0.0_DP
    f_srb = 0.0_DP
    dfdc  = 0.0_DP

    ! jd: This term is only present in GFN0
    if (pub_dftb_method /= DFTB_GFN0) then
       return
    end if

    call utils_assert(par%nat == size(elements), myself//': &
         &allocated parallel strategy is incompatible with elements.')

    ! ab: choose cutoff for periodic translations as per boundary conditions
    if (.not. any(pub_dftb_bc_is_periodic)) then
       lat_cut = 0.0_DP
    else if(all(pub_dftb_bc_is_periodic)) then
       lat_cut = pub_dftb_srb_cutoff
    else
       call utils_abort(myself//': Mixed BCs are not currently supported.')
    end if

    ! ab: lattice vectors:
    lattice(1) = cell%a1
    lattice(2) = cell%a2
    lattice(3) = cell%a3

    ! ab: generate lattice points
    call simulation_cell_lattice_points_init(lattice, lat_cut, points)
    ntr = size(points)

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(STATIC)  &
!$OMP PRIVATE(local_iat,global_iat,orig_iat,zi,rvec,rv,R_i,expterm,dr0dci,itr, &
!$OMP dr0dcj,global_jat,orig_jat,zj,R_j,delta_en,term1,rij,r_ij,delta_r,r0_ij, &
!$OMP      dterm) &
!$OMP SHARED(par,elements,dftb_par,cell,pub_my_proc_id,points,ntr, &
!$OMP      pub_dftb_calculate_force,pub_dftb_srb_cutoff) &
!$OMP REDUCTION(+:e_srb,f_srb,dfdc)
    ! jd: This rank's atoms I
    loop_atom_I:                                                               &
    do local_iat = 1, par%num_atoms_on_proc(pub_my_proc_id)
       global_iat = par%first_atom_on_proc(pub_my_proc_id) + local_iat-1
       orig_iat = par%orig_atom(global_iat)
       zi = elements(orig_iat)%atomic_number

       ! jd: If atomic number not in [z_min,z_max], SRB term is absent
       if(zi < z_min .or. zi > z_max) cycle

       R_i = elements(orig_iat)%centre
       if(zi > dftb_max_species .or. zi < 1) then
          call utils_abort(myself//': atomic number zi out of range: '//&
               trim(utils_int_to_str(zi))//'.', avoid_mpi_calls = .true.)
       end if

       ! jd: All atoms J
       loop_atom_J:                                                            &
       do global_jat = 1, global_iat - 1

          orig_jat = par%orig_atom(global_jat)
          zj = elements(orig_jat)%atomic_number
          ! jd: If atomic number not in [z_min,z_max], SRB term is absent
          if(zj < z_min .or. zj > z_max) cycle

          ! jd: Only heteroatomic bonds
          if(zi == zj) cycle

          R_j = elements(orig_jat)%centre
          if(zj > dftb_max_species .or. zj < 1) then
             call utils_abort(myself//': atomic number zj out of range: '//&
                  trim(utils_int_to_str(zj))//'.', avoid_mpi_calls = .true.)
          end if

          r_ij = minimum_image_distance(R_j,R_i,cell)
          ! ab: distance-based cutoff
          if (r_ij > pub_dftb_srb_cutoff) cycle

          ! ab: delta electronegativity for R0_ij
          delta_en = abs(dftb_par%z(zi)%elecneg - dftb_par%z(zj)%elecneg)
          call dftb_calculate_r0_ij(orig_iat, orig_jat, zi, zj, &
               delta_en, elements, dftb_par, r0_ij, dr0dci, dr0dcj)

          ! ab: delta electronegativity for term1
          delta_en = dftb_par%z(zi)%paulingen - dftb_par%z(zj)%paulingen
          ! jd: -\eta_srb * (1+g_scal*(\Delta EN)^2)
          term1 = dftb_par%srbexp*(1.0_DP+dftb_par%srbken*delta_en*delta_en)

          ! ab: loop over all translational lattice points
          do itr = 1, ntr
             rvec= R_i - R_j + points(itr)
             rij = magnitude(rvec)
             rv = (/ rvec%X, rvec%Y, rvec%Z /)
             delta_r = rij - r0_ij
             expterm = exp(-term1 * delta_r * delta_r)
             e_srb = e_srb + expterm
             dterm = term1 * expterm * delta_r
             if (pub_dftb_calculate_force) then
                f_srb(:,orig_iat) = f_srb(:,orig_iat) + dterm*rv(:)/rij
                f_srb(:,orig_jat) = f_srb(:,orig_jat) - dterm*rv(:)/rij
                dfdc(orig_iat) = dfdc(orig_iat) - dterm * dr0dci
                dfdc(orig_jat) = dfdc(orig_jat) - dterm * dr0dcj
             end if
          end do

       end do loop_atom_J
    end do loop_atom_I
!$OMP END PARALLEL DO

    call comms_reduce('SUM', e_srb)
    ! jd: Apply prefactor k_SRB in front of the sum in [1]:(4)
    e_srb = dftb_par%srbpre * e_srb

    if (pub_dftb_calculate_force) then
       ! ab: two body forces
       call comms_reduce('SUM', f_srb)
       ! ab: three body forces
       call comms_reduce('SUM', dfdc)
       call dense_create(fc, par%nat+1, 1, iscmplx=.false.)
       call dense_create(dfc, par%nat+1, 1, iscmplx=.false.)
       call dense_put_col(dfdc, dfc, 1)
       do ic = 1, 3
          call dense_product(fc, dcdr(ic), dfc, opA='T', opB='N')
          call dense_get_col(dfdc, fc, 1)
          f_srb(ic,:) = f_srb(ic,:) + dfdc(1:par%nat)
       end do
       f_srb = 2.0_DP * dftb_par%srbpre * f_srb
       call dense_destroy(fc)
       call dense_destroy(dfc)
    end if

    ! ab: deallocate lattice points
    call simulation_cell_lattice_points_exit(points)

    if (pub_print_qc) then
       call dftb_qc_print(myself,e_srb,size(elements),f_srb)
    end if

  end subroutine dftb_srb

  ! ---------------------------------------------------------------------------
  ! ---- P R I V A T E   S U B R O U T I N E S   A N D   F U N C T I O N S ----
  ! ---------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine dftb_calculate_r0_ij(global_iat, global_jat, zi, zj, delta_en, &
       elements, dftb_par, r0_ij, dr0dci, dr0dcj)
    !========================================================================!
    ! Calculates rab0 ([1]:(5)) between two atoms.                           !
    ! This is coordination-number-dependent, so it depends on more than the  !
    ! species of the two atoms.                                              !
    !------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 2021/06/20.                               !
    ! Corrected by Arihant Bhandari in Feb 2022.                             !
    ! Derivatives added by Arihant Bhandari in Mar 2022.                     !
    !========================================================================!

    use dftb_parameters, only: DFTB_PARAMS, dftb_max_species
    use ion, only: ELEMENT

    implicit none

    ! jd: Arguments
    integer, intent(in)         :: global_iat, global_jat
    integer, intent(in)         :: zi, zj
    real(kind=DP), intent(in)   :: delta_en
    type(ELEMENT), intent(in)   :: elements(:)
    type(DFTB_PARAMS), intent(in) :: dftb_par
    real(kind=DP), intent(out)  :: r0_ij
    real(kind=DP), intent(out)  :: dr0dci
    real(kind=DP), intent(out)  :: dr0dcj

    ! jd: Local variables
    integer :: row_i, row_j     ! jd: Row in periodic table of atoms I and J
    real(kind=DP) :: r0_i, r0_j ! jd: Per-species radii before CN-dependence
    real(kind=DP) :: r_i, r_j   ! jd: Per-species radii after CN-dependence
    real(kind=DP) :: c1, c2     ! jd: Constants in [1]:(5)
    real(kind=DP) :: c_i, c_j   ! ab: Coordination number (CN)
    real(kind=DP) :: cf_i, cf_j, enf

    integer :: dummy
    integer, parameter :: row_lookup(dftb_max_species) = [ &
         1,1, &
         2,2,2,2,2,2,2,2, &
         3,3,3,3,3,3,3,3, &
         (4,dummy=19,dftb_max_species) &
         ]

    ! -----------------------------------------------------------------------

    row_i = row_lookup(zi)
    row_j = row_lookup(zj)

    r0_i = dftb_par%z(zi)%r0
    r0_j = dftb_par%z(zj)%r0

    cf_i = dftb_par%z(zi)%cnfactor
    cf_j = dftb_par%z(zj)%cnfactor

    ! ab: coordination number
    c_i = elements(global_iat)%dftb_coord
    c_j = elements(global_jat)%dftb_coord

    r_i = r0_i + cf_i*c_i + dftb_par%srbshift
    r_j = r0_j + cf_j*c_j + dftb_par%srbshift

    c1 = 0.005_DP * (dftb_par%p(row_i,1) + dftb_par%p(row_j,1))
    c2 = 0.005_DP * (dftb_par%p(row_i,2) + dftb_par%p(row_j,2))

    enf = (1.0_DP - c1 * delta_en - c2 * delta_en * delta_en)

    r0_ij = (r_i + r_j) * enf

    dr0dci = cf_i * enf
    dr0dcj = cf_j * enf

  end subroutine dftb_calculate_r0_ij

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function dftb_gfn0_dcndr(ri,rj,rc_i,rc_j,c_i,points,elements,dftb_par) result(dcndr)
    !========================================================================!
    ! Calculates the derivative of coordination number d c_i / d rj.         !
    !------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Apr 2022.                               !
    !========================================================================!

    use dftb_parameters, only: DFTB_PARAMS
    use geometry, only: POINT, magnitude, OPERATOR(-)
    use ion, only: ELEMENT
    use timer, only: timer_clock

    implicit none

    ! ab: Arguments
    type(POINT), intent(in)    :: ri, rj       ! ab: atomic coordinates
    real(kind=DP), intent(in)  :: rc_i, rc_j   ! ab: covalent radii
    real(kind=DP), intent(in)  :: c_i          ! ab: coordination number
    type(POINT), intent(in)    :: points(:)    ! ab: translational points
    type(ELEMENT), intent(in)   :: elements(:) ! ab: elements array
    type(DFTB_PARAMS), intent(in) :: dftb_par  ! ab: dftb parameters

    ! ab: Local variables
    integer :: kat, zk
    type(POINT) :: dr, rk
    real(kind=DP) :: rc_k, addition(1:3), dcndr(1:3)
    real(kind=DP), parameter :: low_cutoff = 1D-6 ! ab: For I==J
    character(len=*), parameter :: myself = 'dftb_gfn0_dcndr'

    ! -----------------------------------------------------------------------

    call timer_clock(myself,1)

    addition = 0.0_DP
    dcndr = 0.0_DP
    dr = ri - rj
    if (magnitude(dr) < low_cutoff) then
       do kat = 1, size(elements)
          zk = elements(kat)%atomic_number
          rk = elements(kat)%centre
          rc_k = dftb_par%z(zk)%covrad
          addition = addition + dftb_nsum(ri-rk,rc_i+rc_k,dftb_gfn0_dcijn,points)
       end do
    end if

    dcndr = ( addition-dftb_nsum(ri-rj,rc_i+rc_j,dftb_gfn0_dcijn,points) ) &
         *dftb_gfn0_dcn(c_i)

    call timer_clock(myself,2)

  end function dftb_gfn0_dcndr

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function dftb_gfn0_dcn(c_i) result(dcn)
    !========================================================================!
    ! Used in Calculation of the derivative of coordination number           !
    !------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in May 2022.                               !
    !========================================================================!

    implicit none

    ! ab: Arguments
    real(kind=DP), intent(in)  :: c_i          ! ab: coordination number

    ! ab: Local variables
    real(kind=DP) :: dcn
    real(kind=DP), parameter :: cm = 8.0_DP ! ab: max coordination number

    ! -----------------------------------------------------------------------

    dcn = 1.0_DP - exp(c_i)/(1.0_DP+exp(cm))

  end function dftb_gfn0_dcn

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function dftb_gfn0_dx(ci,ki) result(dx)
    !========================================================================!
    ! Used in Calculation of the derivative of elextronegativity (X).        !
    !------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in May 2022.                               !
    !========================================================================!

    use constants, only: SAFE_DIV_EPS

    implicit none

    ! ab: Arguments
    real(kind=DP), intent(in)  :: ci          ! ab: coordination number
    real(kind=DP), intent(in)  :: ki          ! ab: kappa

    ! ab: Local variables
    real(kind=DP) :: dx

    ! -----------------------------------------------------------------------

    dx = ki / 2.0_DP / (sqrt(ci + SAFE_DIV_EPS))

  end function dftb_gfn0_dx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function dftb_nsum(rij, rc_ij, func, points) result(nsum)
    !========================================================================!
    ! Calculates the sum of func(rij+n,rc_ij) over all translational points. !
    !------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Apr 2022.                               !
    !========================================================================!

    use geometry, only: POINT, OPERATOR(+)

    implicit none

    ! ab: Arguments
    type(POINT), intent(in)    :: rij   ! ab: position vector, ri-rj
    real(kind=DP), intent(in)  :: rc_ij ! ab: a parameter
    procedure(functions)       :: func  ! ab: a function
    type(POINT), intent(in)    :: points(:) ! ab: translational lattice points

    ! ab: Local variables
    integer :: itr
    type(POINT) :: dr
    real(kind=DP) :: nsum(1:3)

    ! -----------------------------------------------------------------------

    nsum = 0.0_DP
    do itr = 1, size(points)
       dr = rij + points(itr)
       nsum = nsum +  func(dr, rc_ij)
    end do

  end function dftb_nsum

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function dftb_gfn0_dcijn(dr, rc_ij) result(dcijn)
    !========================================================================!
    ! Used in the calculation of the derivative of coordination number.      !
    !------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Apr 2022.                               !
    !========================================================================!

    use constants, only: PI
    use geometry, only: POINT, magnitude
    use rundat, only: pub_dftb_coord_cutoff

    implicit none

    ! ab: Arguments
    type(POINT), intent(in)    :: dr    ! ab: dr vector =  ri - rj + n
    real(kind=DP), intent(in)  :: rc_ij ! ab: sum of covalent radii, rc_i+rc_j

    ! ab: Local variables
    real(kind=DP) :: dcijn(1:3), rijn
    real(kind=DP), parameter :: c0 = 7.5_DP ! ab: Hard-coded parameter
    real(kind=DP), parameter :: low_cutoff = 1D-6 ! ab: For excluding I==J

    ! -----------------------------------------------------------------------

    rijn = magnitude(dr)

    if (rijn < low_cutoff .or. rijn > pub_dftb_coord_cutoff) then
       dcijn = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
    else
       dcijn = -c0/rc_ij/sqrt(PI) / rijn &
            * exp(-(c0*(1.0_DP-rijn/rc_ij))**2) * (/ dr%X, dr%Y, dr%Z /)
    end if

  end function dftb_gfn0_dcijn

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_init_coordination(mdl)
    !========================================================================!
    ! Initializes the fractional coordination number of an atom following:   !
    ! GFN0: [1]:(10), GFN1: [2]:(15), GFN2: [3]:(18)                         !
    !------------------------------------------------------------------------!
    ! mdl(inout) : mdl%elements%dftb_coord is calculated.                         !
    !------------------------------------------------------------------------!
    ! Modified by Arihant Bhandari in Aug 2021 from dftb_calculate_coord()   !
    ! originally written by Jacek Dziedzic.                                  !
    !========================================================================!

    use comms, only: pub_my_proc_id, comms_reduce
    use constants, only: DFTB_GFN0, DFTB_GFN1, DFTB_GFN2
    use dense, only: dense_create, dense_put_element
    use geometry, only: POINT, magnitude, OPERATOR(+), OPERATOR(-), &
         OPERATOR(*)
    use model_type, only: MODEL
    use rundat, only: pub_dftb_method, pub_dftb_bc_is_periodic, &
         pub_dftb_coord_cutoff, pub_dftb_calculate_force, pub_threads_max
    use simulation_cell, only: simulation_cell_lattice_points_init, &
         simulation_cell_lattice_points_exit
    use timer, only: timer_clock
    use utils, only: utils_erf, utils_assert, utils_abort

    implicit none

    ! ab: Arguments
    type(MODEL), intent(inout) :: mdl

    ! ab: Local variables
    type(POINT) :: lattice(3) ! ab: lattice vectors
    integer :: iat, jat, itr, ntr, ic, nat, local_iat, global_iat, local_jat, &
         global_jat
    integer :: zi, zj      ! jd: Atomic numbers of atoms I and J
    type(POINT) :: R_i, R_j  ! jd: Positions of atoms I, J
    type(POINT) :: vec       ! ab: translational vector between I and J
    real(kind=DP) :: R_ij    ! jd: Distance between atoms I, J
    real(kind=DP) :: R_0     ! jd: Sum of covalent radii
    real(kind=DP) :: Rc_i, Rc_j, ci, dcnidrj(3), dxidrj(3), dxi
    real(kind=DP) :: cn      ! jd: pairwise contribution
    real(kind=DP) :: cutoff  ! ab: cutoff for periodic translations
    real(kind=DP) :: ex0, ex1, ex2, coord(mdl%nat)
    real(kind=DP), parameter :: cmax = 8.0_DP ! ab: max coordination number
    real(kind=DP), parameter :: c0 = 7.5_DP ! jd: Hard-coded parameter
    real(kind=DP), parameter :: k1 = 16.0_DP ! ab: Hard-coded parameter
    real(kind=DP), parameter :: k2 = 10.0_DP ! ab: Hard-coded parameter
    real(kind=DP), parameter :: p4 = 2.0_DP ! ab: Hard-coded parameter
    real(kind=DP), parameter :: low_cutoff = 1D-6 ! jd: For excluding I==J
    character(len=*), parameter :: myself = 'dftb_init_coordination'
    type(POINT), allocatable :: coord_points(:)
#ifdef SCALAPACK
    integer :: nprow, npcol, myrow, mycol, lj, li, nb
    ! ab: BLACS subroutines
    external :: blacs_gridinfo
#endif

    ! -----------------------------------------------------------------------

    call timer_clock(myself,1)

    call utils_assert(mdl%nsub == 1, &
       myself//': DFTB does not support embedding regions')

    ! ab: lattice vectors:
    lattice(1) = mdl%cell%a1
    lattice(2) = mdl%cell%a2
    lattice(3) = mdl%cell%a3

    ! ab: check boundary conditions and set cutoff for periodic translations
    if(.not. any(pub_dftb_bc_is_periodic)) then
       cutoff = 0.0_DP
    else if (all(pub_dftb_bc_is_periodic)) then
       cutoff = pub_dftb_coord_cutoff
    else
       call utils_abort(myself//': Mixed BCs are not currently supported.')
    end if

    ! ab: initialize translational lattice points
    call simulation_cell_lattice_points_init(lattice, cutoff, coord_points)

    ! ab: initialize coordination number
    coord = 0.0_DP
    nat = size(mdl%elements)

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(STATIC)  &
!$OMP PRIVATE(local_iat,global_iat,zi,R_i,cn,itr,ex0,ex1,ex2,                  &
!$OMP      zj,R_j,R_ij,R_0,iat,jat,vec)                                        &
!$OMP SHARED(mdl,coord_points,pub_my_proc_id,pub_dftb_coord_cutoff,            &
!$OMP      pub_dftb_method) &
!$OMP REDUCTION(+:coord)
    do local_iat = 1, mdl%regions(1)%par%num_atoms_on_proc(pub_my_proc_id)
       global_iat = mdl%regions(1)%par%first_atom_on_proc(pub_my_proc_id) + local_iat-1
       iat = mdl%regions(1)%par%orig_atom(global_iat)
       R_i = mdl%elements(iat)%centre
       zi = mdl%elements(iat)%atomic_number
       do jat = 1, iat
          R_j = mdl%elements(jat)%centre
          zj = mdl%elements(jat)%atomic_number
          do itr = 1, size(coord_points)
             vec =  R_i - R_j + coord_points(itr)
             R_ij = magnitude(vec)

             ! jd: Exclude i==j, but not i==image(j)
             if(R_ij < low_cutoff .or. R_ij > pub_dftb_coord_cutoff) cycle

             R_0 = mdl%dftb_par%z(zi)%covrad + mdl%dftb_par%z(zj)%covrad

             if (pub_dftb_method == DFTB_GFN0) then
                cn = (1.0_DP + utils_erf(-c0 * (R_ij-R_0)/R_0)) * 0.5_DP
                ! ab: commenting out the derivative as it is calculated on the
                ! fly now, but retaining the formula here temporarily
                !dcn=-c0/sqrt(PI)/R_0*&
                !     exp(-c0**2*(R_ij-R_0)**2/R_0**2)
             else if (pub_dftb_method == DFTB_GFN1) then
                ex0= exp(-k1*(R_0/R_ij-1.0_DP))
                cn = 1.0_DP / (1.0_DP + ex0)
                !dcn=-k1*R_0*ex0/(R_ij**2*((ex0+1.0_DP)**2))
             else if (pub_dftb_method == DFTB_GFN2) then
                ex1= exp(-k2*(R_0/R_ij-1.0_DP))
                ex2= exp(-p4*k2*((R_0+p4)/R_ij-1.0_DP))
                cn = 1.0_DP / (1.0_DP + ex1) * 1.0_DP / (1.0_DP + ex2)
                !dcn=-k2*R_0*ex1/(R_ij**2*((ex1+1.0_DP)**2))*1.0_DP/(1.0_DP+ex2)&
                !     -1.0_DP/(1.0_DP+ex1)*p4*k2*(R_0+p4)*ex2/&
                !     (R_ij**2*((ex2+1.0_DP)**2))
             else
                call utils_abort(myself//': Unsupported DFTB method. &
                     &Only GFN0, GFN1 and GFN2 are supported.')
             end if

             coord(iat) = coord(iat) + cn

             if (iat /= jat) then ! ab: avoid double counting
                coord(jat) = coord(jat) + cn
             end if

          end do
       end do
    end do
!$OMP END PARALLEL DO

    call comms_reduce('SUM',coord)

    ! ab: cut coordination number in GFN0
    if (pub_dftb_method == DFTB_GFN0) then
       do iat = 1, nat
          coord(iat) = log(1.0_DP + exp(cmax)) - &
               log(1.0_DP + exp(cmax - coord(iat)))
       end do
    end if

    mdl%elements(1:nat)%dftb_coord = coord(1:nat)
    mdl%regions(1)%elements(1:nat)%dftb_coord = coord(1:nat)

    ! ab: also calculate the derivative in forces calculation
    if (pub_dftb_calculate_force .and. pub_dftb_method == DFTB_GFN0) then

       do ic = 1, 3
          call dense_create(mdl%dftb_var%dcdr(ic), nat+1, nat+1, iscmplx=.false.)
          call dense_create(mdl%dftb_var%dxdr(ic), nat+1, nat+1, iscmplx=.false.)
       end do

#ifdef SCALAPACK

       ! ab: block-cyclically distributed matrix
       call blacs_gridinfo(mdl%dftb_var%dcdr(1)%blacs_desc(2), nprow, npcol, &
            myrow, mycol)
       nb = mdl%dftb_var%dcdr(1)%blacs_nb

       do lj = 0, mdl%dftb_var%dcdr(1)%blacs_ncol - 1
          jat=((lj/nb)*npcol+mycol)*nb + mod(lj,nb) + 1

          do li = 0,  mdl%dftb_var%dcdr(1)%blacs_ld - 1
             iat=((li/nb)*nprow+myrow)*nb + mod(li,nb) + 1

             if (jat==nat+1 .or. iat==nat+1) cycle

             zj = mdl%elements(jat)%atomic_number
             r_j = mdl%elements(jat)%centre
             rc_j = mdl%dftb_par%z(zj)%covrad

             zi = mdl%elements(iat)%atomic_number
             ci = mdl%elements(iat)%dftb_coord
             r_i = mdl%elements(iat)%centre
             rc_i = mdl%dftb_par%z(zi)%covrad
             dxi = dftb_gfn0_dx(ci,mdl%dftb_par%z(zi)%kappa)
             dcnidrj = dftb_gfn0_dcndr(r_i, r_j, rc_i, rc_j, ci, &
                  coord_points, mdl%elements, mdl%dftb_par)
             dxidrj = dxi * dcnidrj
             do ic = 1, 3
                call dense_put_element(dcnidrj(ic), mdl%dftb_var%dcdr(ic), iat, jat)
                call dense_put_element(dxidrj(ic), mdl%dftb_var%dxdr(ic), iat, jat)
             end do

          end do

       end do

#else

       ! ab: undistributed matrices, (allocated on all procs)
       do local_jat = 1, mdl%regions(1)%par%num_atoms_on_proc(pub_my_proc_id)
          global_jat = mdl%regions(1)%par%first_atom_on_proc(pub_my_proc_id)+local_jat-1
          jat = mdl%regions(1)%par%orig_atom(global_jat)
          zj = mdl%elements(jat)%atomic_number
          r_j = mdl%elements(jat)%centre
          rc_j = mdl%dftb_par%z(zj)%covrad

          do iat = 1 , nat
             zi = mdl%elements(iat)%atomic_number
             ci = mdl%elements(iat)%dftb_coord
             r_i = mdl%elements(iat)%centre
             rc_i = mdl%dftb_par%z(zi)%covrad
             dxi = dftb_gfn0_dx(ci,mdl%dftb_par%z(zi)%kappa)
             dcnidrj = dftb_gfn0_dcndr(r_i, r_j, rc_i, rc_j, ci, &
                  coord_points, mdl%elements, mdl%dftb_par)
             dxidrj = dxi * dcnidrj

             do ic = 1, 3
                mdl%dftb_var%dcdr(ic)%dmtx(iat,jat) = dcnidrj(ic)
                mdl%dftb_var%dxdr(ic)%dmtx(iat,jat) = dxidrj(ic)
             end do

          end do

       end do

       do ic = 1, 3
          call comms_reduce('SUM', mdl%dftb_var%dcdr(ic)%dmtx)
          call comms_reduce('SUM', mdl%dftb_var%dxdr(ic)%dmtx)
       end do

#endif

    end if

    ! ab: deallocate translational lattice points
    call simulation_cell_lattice_points_exit(coord_points)

    call timer_clock(myself,2)

  end subroutine dftb_init_coordination

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_read_params(mdl, par_filename_method, par_filename_common)
    !==========================================================================!
    ! Reads DFTB parameters in the GFN-xTB format from par_filename_method,    !
    ! and common parameters from par_filename_common.                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !    mdl (inout): %dftb_par is modified.                                   !
    !    par_filename_method (in): The filename for the method-specific DFTB   !
    !                              parameters.                                 !
    !    par_filename_common (in): The filename for the method-agnostic DFTB   !
    !                              parameters.                                 !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 2021/06/03.                                 !
    ! Extended with the common parameters by Jacek Dziedzic on 2021/06/23.     !
    !==========================================================================!

    use comms, only: pub_on_root
    use dftb_parameters, only: dftb_max_species, dftb_num_species, dftb_strlen
    use model_type, only: MODEL
    use utils, only: utils_unit, utils_abort, utils_open_unit_check, &
         utils_close_unit_check, utils_int_to_str, utils_assert, &
         utils_alloc_check

    implicit none

    ! jd: Arguments
    type(MODEL), intent(inout)   :: mdl
    character(len=*), intent(in) :: par_filename_method
    character(len=*), intent(in) :: par_filename_common

    ! jd: Local variables
    real(kind=DP) :: cur_z_real
    real(kind=DP) :: pairparval
    logical :: match
    logical :: in_pairpar
    integer :: iunit
    integer :: cur_z
    integer :: z1, z2
    integer :: ierr
    integer :: ci
    integer :: nshells, ishell, cur_shell, nprimitive
    character :: c1
    character(len=dftb_strlen) :: line ! jd: line read from param file
    character(len=dftb_strlen) :: line_nls ! jd: with leading spaces removed
    character(len=dftb_strlen) :: word ! jd: first word in line
    character(len=*), parameter :: myself = 'dftb_read_params'

    ! -------------------------------------------------------------------------

    ! ----------------------------------------------
    ! Parameters of the DFTB method (GFN0/GFN1/GFN2)
    ! ----------------------------------------------

    if(pub_on_root) then
       write(stdout,'(a)') 'Reading DFTB-method parameters from "'//&
            trim(par_filename_method)//'".'
    end if

    ! jd: Read on *all* ranks, as bcasting a DFTB_PARAMS structure
    !     would be clunky.
    iunit = utils_unit()
    open(unit=iunit,file=trim(par_filename_method),action='read',iostat=ierr)
    call utils_open_unit_check(myself,trim(par_filename_method),ierr)

    ! jd: Read all the lines in the param file
    cur_shell = 0
    cur_z = 0
    in_pairpar = .false.
    do
       read(iunit,'(a)',err=100,iostat=ierr) line
       if(ierr /= 0) exit

       line_nls = adjustl(line)

       ! jd: Parse line, skpping empty lines and section delimiters.
       !     Pay attention to whether we are inside a $pairpar section, as
       !     this one has no keys.
       if(len(trim(line_nls)) == 0) cycle
       if(line_nls(1:8) == "$pairpar") in_pairpar = .true.
       if(line_nls(1:4) == "$end") in_pairpar = .false.
       if(line_nls(1:1) == "$" .and. line_nls(2:2) /= "Z") cycle
       read(line_nls,*) word

       ! jd: Global information about the parameter set
       call internal_parse_string(line_nls,'level', mdl%dftb_par%level)
       call internal_parse_string(line_nls,'name', mdl%dftb_par%name)
       call internal_parse_string(line_nls,'doi', mdl%dftb_par%doi)

       ! jd: Global parameters common to GFN0, GFN1, GFN2
       call internal_parse_real(line_nls,'ks', mdl%dftb_par%ks)
       call internal_parse_real(line_nls,'kp', mdl%dftb_par%kp)
       call internal_parse_real(line_nls,'kd', mdl%dftb_par%kd)
       call internal_parse_real(line_nls,'kdiff', mdl%dftb_par%kdiff)
       call internal_parse_real(line_nls,'ipeashift', mdl%dftb_par%ipeashift)
       call internal_parse_real(line_nls,'a1', mdl%dftb_par%a1)
       call internal_parse_real(line_nls,'a2', mdl%dftb_par%a2)
       call internal_parse_real(line_nls,'s8', mdl%dftb_par%s8)
       call internal_parse_real(line_nls,'s9', mdl%dftb_par%s9)
       call internal_parse_real(line_nls,'kexp', mdl%dftb_par%kexp)
       call internal_parse_real(line_nls,'kexplight', mdl%dftb_par%kexplight)

       call internal_parse_real(line_nls,'ens', mdl%dftb_par%ens)
       call internal_parse_real(line_nls,'enp', mdl%dftb_par%enp)
       call internal_parse_real(line_nls,'end', mdl%dftb_par%end)
       call internal_parse_real(line_nls,'enscale4', mdl%dftb_par%enscale4)
       call internal_parse_real(line_nls,'srbshift', mdl%dftb_par%srbshift)
       call internal_parse_real(line_nls,'srbpre', mdl%dftb_par%srbpre)
       call internal_parse_real(line_nls,'srbexp', mdl%dftb_par%srbexp)
       call internal_parse_real(line_nls,'srbken', mdl%dftb_par%srbken)
       call internal_parse_real(line_nls,'renscale', mdl%dftb_par%renscale)

       ! jd: Global parameters for GFN1 only
       call internal_parse_real(line_nls,'ksp', mdl%dftb_par%ksp)
       call internal_parse_real(line_nls,'cns', mdl%dftb_par%cns)
       call internal_parse_real(line_nls,'cnp', mdl%dftb_par%cnp)
       call internal_parse_real(line_nls,'cnd1', mdl%dftb_par%cnd1)
       call internal_parse_real(line_nls,'cnd2', mdl%dftb_par%cnd2)
       call internal_parse_real(line_nls,'xbdamp', mdl%dftb_par%xbdamp)
       call internal_parse_real(line_nls,'xbrad', mdl%dftb_par%xbrad)

       ! jd: Global parameters for GFN2 only
       call internal_parse_real(line_nls,'ksd', mdl%dftb_par%ksd)
       call internal_parse_real(line_nls,'kpd', mdl%dftb_par%kpd)
       call internal_parse_real(line_nls,'gam3s', mdl%dftb_par%gam3s)
       call internal_parse_real(line_nls,'gam3p', mdl%dftb_par%gam3p)
       call internal_parse_real(line_nls,'gam3d1', mdl%dftb_par%gam3d1)
       call internal_parse_real(line_nls,'gam3d2', mdl%dftb_par%gam3d2)
       call internal_parse_real(line_nls,'aesshift', mdl%dftb_par%aesshift)
       call internal_parse_real(line_nls,'aesexp', mdl%dftb_par%aesexp)
       call internal_parse_real(line_nls,'aesrmax', mdl%dftb_par%aesrmax)
       call internal_parse_real(line_nls,'aesdmp3', mdl%dftb_par%aesdmp3)
       call internal_parse_real(line_nls,'aesdmp5', mdl%dftb_par%aesdmp5)

       ! jd: Global parameters for GFN1 and GFN2 only
       call internal_parse_real(line_nls,'enscale', mdl%dftb_par%enscale)
       call internal_parse_real(line_nls,'alphaj', mdl%dftb_par%alphaj)

       ! jd: Per-species tag "$Z="
       ! jd: For double-digit Z's there is no space between "$Z=" and the number.
       !     Make sure to put it there so the parser can dinstinguish the fields.
       if(line_nls(1:2) == '$Z') line_nls(3:3) = ' '
       cur_z_real = garbage_real ! jd: Needed or parse_real() will panic
       call internal_parse_real(line_nls,'$Z', cur_z_real, match)
       if(match) then
          cur_z = nint(cur_z_real)
          call utils_assert(cur_z <= dftb_max_species .and. cur_z > 0, myself//&
               ': Z value out of range for a per-species parameter. Z='//&
               trim(utils_int_to_str(cur_z))//'.')
       end if

       ! ----------------------------------------------------------------------
       ! jd: Per-species parameters
       ! ----------------------------------------------------------------------

       ! jd: There is no space between "ao=" and the value.
       !     Make sure to put it there so the parser can dinstinguish the fields.
       if(line_nls(1:2) == 'ao') line_nls(3:3) = ' '
       call internal_parse_string(line_nls,'ao', mdl%dftb_par%z(cur_z)%ao, match)

       ! ab: store the number of shells and allocate STO-nG datastructure
       if (match) then
          nshells = 0
          do ci = 1, len(mdl%dftb_par%z(cur_z)%ao)
             c1 = mdl%dftb_par%z(cur_z)%ao(ci:ci)
             select case(c1)
             case('1','2','3','4','5','6','7','8',' ')
             case('s')
                nshells = nshells + 1
             case('p')
                nshells = nshells + 1
             case('d')
                nshells = nshells + 1
             case default
                call utils_abort(myself//': Unrecognised subshell "'//c1//&
                     '" in DFTB AO string "'&
                     //trim(mdl%dftb_par%z(cur_z)%ao)//'".')
             end select
          end do
          mdl%dftb_par%z(cur_z)%nshells = nshells
          allocate(mdl%dftb_par%z(cur_z)%ang(nshells), stat=ierr)
          call utils_alloc_check(myself, 'mdl%dftb_par%z(cur_z)%ang(nshells)', ierr)
          allocate(mdl%dftb_par%z(cur_z)%sto(nshells), stat=ierr)
          call utils_alloc_check(myself, 'mdl%dftb_par%z(cur_z)%sto(nshells)', ierr)

          ! ab: store the name of individual ao and angular momentum
          do ishell = 1, nshells
             mdl%dftb_par%z(cur_z)%sto(ishell)%ao = &
                  mdl%dftb_par%z(cur_z)%ao(2*ishell-1:2*ishell)
             c1 = mdl%dftb_par%z(cur_z)%sto(ishell)%ao(2:2)
             select case(c1)
             case('s')
                mdl%dftb_par%z(cur_z)%ang(ishell) = 0
             case('p')
                mdl%dftb_par%z(cur_z)%ang(ishell) = 1
             case('d')
                mdl%dftb_par%z(cur_z)%ang(ishell) = 2
             case default
                call utils_abort(myself//': Unrecognised subshell "'//c1//&
                     '" in DFTB AO string "'&
                     //trim(mdl%dftb_par%z(cur_z)%sto(ishell)%ao)//'".')
             end select
          end do
       end if


       ! ab: store the number of primitives
       if(line_nls(3:7) == ' STO-') then
          ! ab: check which shell
          do ishell = 1, mdl%dftb_par%z(cur_z)%nshells
             if (mdl%dftb_par%z(cur_z)%sto(ishell)%ao == line_nls(1:2)) then
                cur_shell=ishell
             end if
          end do

          read(line_nls(8:8),'(i1)',err=100) nprimitive
          mdl%dftb_par%z(cur_z)%sto(cur_shell)%nprim= nprimitive

          ! ab: allocate space for STO-nG exponents and contraction coefficients
          allocate(mdl%dftb_par%z(cur_z)%sto(cur_shell)%alpha(nprimitive), &
               stat=ierr)
          call utils_alloc_check(myself, &
               'mdl%dftb_par%z(cur_z)%sto(cur_shell)%alpha(nprimitive)', ierr)
          mdl%dftb_par%z(cur_z)%sto(cur_shell)%alpha(:)=garbage_real

          allocate(mdl%dftb_par%z(cur_z)%sto(cur_shell)%coeff(nprimitive), &
               stat=ierr)
          call utils_alloc_check(myself, &
               'mdl%dftb_par%z(cur_z)%sto(cur_shell)%coeff(nprimitive)', ierr)
          mdl%dftb_par%z(cur_z)%sto(cur_shell)%coeff(:)=garbage_real
       end if

       ! ab: store STO-nG exponents and contraction coefficients
       if (line_nls(1:3)=='sto') then
          call internal_parse_reals(mdl%dftb_par%z(cur_z)%sto(cur_shell)%nprim, &
               line_nls,'sto_alpha=', mdl%dftb_par%z(cur_z)%sto(cur_shell)%alpha)
          call internal_parse_reals(mdl%dftb_par%z(cur_z)%sto(cur_shell)%nprim, &
               line_nls,'sto_coeff=', mdl%dftb_par%z(cur_z)%sto(cur_shell)%coeff)
       end if

       if (line_nls(1:3)=='lev') then
          call internal_parse_reals(mdl%dftb_par%z(cur_z)%nshells, &
               line_nls,'lev=', mdl%dftb_par%z(cur_z)%lev)
       end if

       if (line_nls(1:3)=='exp') then
          call internal_parse_reals(mdl%dftb_par%z(cur_z)%nshells, &
               line_nls,'exp=', mdl%dftb_par%z(cur_z)%exp)
       end if

       ! jd: Per-species parameters common to GFN0, GFN1, GFN2
       call internal_parse_real(line_nls,'GAM=', mdl%dftb_par%z(cur_z)%gam)

       call internal_parse_real(line_nls,'POLYS=', mdl%dftb_par%z(cur_z)%polys)
       call internal_parse_real(line_nls,'POLYP=', mdl%dftb_par%z(cur_z)%polyp)
       call internal_parse_real(line_nls,'POLYD=', mdl%dftb_par%z(cur_z)%polyd)
       call internal_parse_real(line_nls,'REPA=', mdl%dftb_par%z(cur_z)%repa)
       call internal_parse_real(line_nls,'REPB=', mdl%dftb_par%z(cur_z)%repb)

       ! jd: Per-species parameters for GFN0 only
       call internal_parse_real(line_nls,'ALPG=', mdl%dftb_par%z(cur_z)%alpg)
       call internal_parse_real(line_nls,'EN=', mdl%dftb_par%z(cur_z)%en)
       call internal_parse_real(line_nls,'KAPPA=', mdl%dftb_par%z(cur_z)%kappa)
       call internal_parse_real(line_nls,'KQAT2=', mdl%dftb_par%z(cur_z)%kqat2)
       call internal_parse_real(line_nls,'KQS=', mdl%dftb_par%z(cur_z)%kqs)
       call internal_parse_real(line_nls,'KQP=', mdl%dftb_par%z(cur_z)%kqp)
       call internal_parse_real(line_nls,'KQD=', mdl%dftb_par%z(cur_z)%kqd)
       call internal_parse_real(line_nls,'XI=', mdl%dftb_par%z(cur_z)%xi)

       ! jd: Per-species parameters for GFN1 only
       call internal_parse_real(line_nls,'CXB=', mdl%dftb_par%z(cur_z)%cxb)

       ! jd: Per-species parameters for GFN2 only
       call internal_parse_real(line_nls,'DPOL=', mdl%dftb_par%z(cur_z)%dpol)
       call internal_parse_real(line_nls,'QPOL=', mdl%dftb_par%z(cur_z)%qpol)

       ! jd: Per-species parameters for GFN0 and GFN2 only
       call internal_parse_real(line_nls,'KCNS=', mdl%dftb_par%z(cur_z)%kcns)
       call internal_parse_real(line_nls,'KCNP=', mdl%dftb_par%z(cur_z)%kcnp)
       call internal_parse_real(line_nls,'KCND=', mdl%dftb_par%z(cur_z)%kcnd)

       ! jd: Per-species parameters for GFN1 and GFN2 only
       call internal_parse_real(line_nls,'GAM3=', mdl%dftb_par%z(cur_z)%gam3)

       ! ----------------------------------------------------------------------
       ! jd: Per-species-pair parameters
       ! ----------------------------------------------------------------------

       ! jd: These are handled differently, because there is no key in these
       !     lines. Instead 'in_pairpar' tells us if we are inside of a
       !     per-species-pair block.
       if(in_pairpar) then
          read(line_nls,*) z1, z2, pairparval
          call utils_assert(z1 <= dftb_max_species .and. z1 > 0 .and. &
               z2 <= dftb_max_species .and. z2 > 0, myself//&
               ': Z value out of range for a per-species-pair parameter. Z1='//&
               trim(utils_int_to_str(z1))//', Z2='//trim(utils_int_to_str(z1))&
               //'.')
          mdl%dftb_par%zz(z1,z2)%x = pairparval
          mdl%dftb_par%zz(z2,z1)%x = pairparval
       end if

    end do

    dftb_num_species = cur_z

    close(iunit,iostat=ierr)
    call utils_close_unit_check(myself,trim(par_filename_method),ierr)

    if(pub_debug_on_root) then
       write(stdout,'(a)') 'DFTB-method parameters for '//&
            trim(mdl%dftb_par%name)//' read correctly.'
    end if

    ! ---------------------------------------------------------------
    ! Parameters that are common to all methods (electronegativities,
    ! CN factors and such). These are read from a separate file.
    ! ---------------------------------------------------------------

    if(pub_on_root) then
       write(stdout,'(a)') 'Reading DFTB-common parameters from "'//&
            trim(par_filename_common)//'".'
    end if

    ! jd: Read on *all* ranks, as bcasting a mdl%dftb_parAMS structure
    !     would be clunky.
    iunit = utils_unit()
    open(unit=iunit,file=trim(par_filename_common),action='read',iostat=ierr)
    call utils_open_unit_check(myself,trim(par_filename_common),ierr)

    ! jd: Read all the lines in the param file
    cur_z = 0
    do
       read(iunit,'(a)',err=100,iostat=ierr) line
       if(ierr /= 0) exit

       line_nls = adjustl(line)

       ! jd: Parse line, skpping empty lines and section delimiters.
       !     Pay attention to whether we are inside a $pairpar section, as
       !     this one has no keys.
       if(len(trim(line_nls)) == 0) cycle
       if(line_nls(1:1) == "$" .and. line_nls(2:2) /= "Z") cycle
       read(line_nls,*) word

       ! jd: Global parameters -- only 'p' so far
       call internal_parse_real(line_nls,'p11', mdl%dftb_par%p(1,1))
       call internal_parse_real(line_nls,'p21', mdl%dftb_par%p(2,1))
       call internal_parse_real(line_nls,'p31', mdl%dftb_par%p(3,1))
       call internal_parse_real(line_nls,'p41', mdl%dftb_par%p(4,1))
       call internal_parse_real(line_nls,'p12', mdl%dftb_par%p(1,2))
       call internal_parse_real(line_nls,'p22', mdl%dftb_par%p(2,2))
       call internal_parse_real(line_nls,'p32', mdl%dftb_par%p(3,2))
       call internal_parse_real(line_nls,'p42', mdl%dftb_par%p(4,2))

       ! jd: For double-digit Z's there is no space between "$Z=" and the number.
       !     Make sure to put it there so the parser can dinstinguish the fields.
       if(line_nls(1:2) == '$Z') line_nls(3:3) = ' '
       cur_z_real = garbage_real ! jd: Needed or parse_real() will panic
       call internal_parse_real(line_nls,'$Z', cur_z_real, match)
       if(match) then
          cur_z = nint(cur_z_real)
          call utils_assert(cur_z <= dftb_max_species .and. cur_z > 0, myself//&
               ': Z value out of range for a per-species parameter. Z='//&
               trim(utils_int_to_str(cur_z))//'.')
       end if

       ! ----------------------------------------------------------------------
       ! jd: Per-species parameters
       ! ----------------------------------------------------------------------

       call internal_parse_real(line_nls,'ELECNEG=', mdl%dftb_par%z(cur_z)%elecneg)
       call internal_parse_real(line_nls,'R0=', mdl%dftb_par%z(cur_z)%r0)
       call internal_parse_real(line_nls,'CNFACTOR=', mdl%dftb_par%z(cur_z)%cnfactor)
       call internal_parse_real(line_nls,'COVRAD=', mdl%dftb_par%z(cur_z)%covrad)
       call internal_parse_real(line_nls,'ATOMICRAD=',mdl%dftb_par%z(cur_z)%atomicrad)
       call internal_parse_real(line_nls,'PAULINGEN=',mdl%dftb_par%z(cur_z)%paulingen)

    end do

    close(iunit,iostat=ierr)
    call utils_close_unit_check(myself,trim(par_filename_common),ierr)

    if(pub_debug_on_root) then
       write(stdout,'(a)') 'DFTB-common parameters read correctly.'
    end if

    return

100   call utils_abort(myself//&
           ': Error reading "'//trim(par_filename_method)//'".')
200   call utils_abort(myself//&
           ': Error reading "'//trim(par_filename_common)//'".')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains

    subroutine internal_parse_string(line, key, store, match)
      !========================================================================!
      ! Parses a single line, and if the first word matches 'key', the second  !
      ! word is stored in 'store' and 'match', if present, is set to .true.    !
      ! If the first word does not match 'key', 'store' remains unchanged, and !
      ! 'match', if present, is set to .false..                                !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic on 2021/06/04.                               !
      !========================================================================!

      implicit none

      ! jd: Arguments
      character(len=*), intent(in)              :: line
      character(len=*), intent(in)              :: key
      character(len=dftb_strlen), intent(inout) :: store
      logical, intent(out), optional            :: match

      ! jd: Local variables
      character(len=dftb_strlen) :: dummy
      character(len=*), parameter :: myself = &
           'dftb_read_params::internal_parse_string'

      ! -----------------------------------------------------------------------

      ! jd: Check for a match.
      if(line(1:len(key)) == trim(key))then
         read(line,*,err=100) dummy, store
         if(present(match)) match = .true.
      else
         if(present(match)) match = .false.
      end if

      return

100   call utils_abort(myself//': Error reading or parsing "'//&
           trim(par_filename_method)//'" or "'//&
           trim(par_filename_common)//'" (string key: "'//trim(key)//'").')

    end subroutine internal_parse_string

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_parse_real(line, key, store, match)
      !========================================================================!
      ! Parses a single line, and if the first word matches 'key', the second  !
      ! word (which must be a real value) is stored in 'store' and 'match',    !
      ! if present, is set to .true. If the first word does not match 'key',   !
      ! 'store' remains unchanged, and 'match', if present, is set to .false.. !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic on 2021/06/04.                               !
      !========================================================================!

      use utils, only: utils_real_to_str, utils_assert

      implicit none

      ! jd: Arguments
      character(len=*), intent(in)   :: line
      character(len=*), intent(in)   :: key
      real(kind=DP), intent(inout)   :: store
      logical, intent(out), optional :: match

      ! jd: Local variables
      character(len=dftb_strlen) :: dummy
      real(kind=DP) :: to_store
      character(len=*), parameter :: myself = &
           'dftb_read_params::internal_parse_real'

      ! -----------------------------------------------------------------------

      ! jd: Check for a match. The comparison against " " ensures keys like "ks"
      !     do not match lines like "ksp".
      if(line(1:len(key)) == trim(key) .and. &
           (line(len(key)+1:len(key)+1) == " ")) then
         read(line,*,err=100) dummy, to_store
         call utils_assert(store == garbage_real, myself//&
              ': Trying to store a parameter that has already been stored. &
              &This could indicate a malformed parameter file &
              &(real-valued key: "'//trim(key)//'"). '//&
              'Previous value: '//trim(utils_real_to_str(store))//&
              '. Trying to store: '//trim(utils_real_to_str(to_store)))
         store = to_store
         if(present(match)) match = .true.
      else
         if(present(match)) match = .false.
      end if

      return

100   call utils_abort(myself//': Error reading or parsing "'//&
           trim(par_filename_method)//'" or "'//&
           trim(par_filename_common)//'" (real-valued key: "'//trim(key)//'").')

    end subroutine internal_parse_real

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_parse_reals(num, line, key, store, match)
      !========================================================================!
      ! Parses a single line, and if the first word matches 'key', subsequent  !
      ! words, up to 'num' of them, which must all be real values, are stored  !
      ! in 'store' and 'match', if present, is set to .true. If the first word !
      ! does not match 'key', 'store' remains unchanged, and 'match',          !
      ! if present, is set to .false..                                         !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic on 2021/06/05.                               !
      !========================================================================!

      use utils, only: utils_int_to_str, utils_real_to_str, utils_assert

      implicit none

      ! jd: Arguments
      integer, intent(in)            :: num
      character(len=*), intent(in)   :: line
      character(len=*), intent(in)   :: key
      real(kind=DP), intent(inout)   :: store(num)
      logical, intent(out), optional :: match

      ! jd: Local variables
      character(len=dftb_strlen) :: dummy
      real(kind=DP) :: to_store(num)
      integer :: i
      integer :: ios
      character(len=*), parameter :: myself = &
           'dftb_read_params::internal_parse_reals'

      ! -----------------------------------------------------------------------

      ! jd: Check for a match. The comparison against " " ensures keys like "ks"
      !     do not match lines like "ksp".
      if(line(1:len(key)) == trim(key) .and. &
           (line(len(key)+1:len(key)+1) == " ")) then
         to_store(:) = garbage_real
         read(line,*,err=100,iostat=ios) dummy, to_store(1:num)
         if(ios /= 0) goto 100 ! error condition
         do i = 1, num
            call utils_assert(store(i) == garbage_real, myself//&
                 ': Trying to store a parameter that has already been stored. &
                 &This could indicate a malformed parameter file &
                 &(real'//trim(utils_int_to_str(num))//'-valued key: "'//trim(key)//&
                 '"). '//'Previous value: '//trim(utils_real_to_str(store(i)))&
                 //'. Trying to store: '//trim(utils_real_to_str(to_store(i))))
            store(i) = to_store(i)
         end do
         if(present(match)) match = .true.
      else
         if(present(match)) match = .false.
      end if

      return

100   call utils_abort(myself//': Error reading or parsing "'//&
           trim(par_filename_method)//'" or "'//&
           trim(par_filename_common)//'" (real('//trim(utils_int_to_str(num))&
           //')-valued key: "'//trim(key)//'").')

    end subroutine internal_parse_reals

  end subroutine dftb_read_params

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_adjust_num_ngwfs(mdl)
    !==========================================================================!
    ! Infers the number of NGWFs from the DFTB parameter file for every ONETEP !
    ! species. Adjusts the number of NGWFs in the elements(:) and species(:)   !
    ! arrays -- both in mdl and in mdl%regions(1) on all MPI ranks.            !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   mdl (inout): %dftb_par is needed for the per-species AO strings.       !
    !                %elements(:) and %species(:) is adjusted, along with the  !
    !                regions(1) versions.                                      !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Does not support embedding regions, an assert guards against that.     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 2021/06/26.                                 !
    !==========================================================================!

    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast
    use dftb_parameters, only: dftb_max_species
    use model_type, only: MODEL
    use utils, only: utils_int_to_str, utils_assert

    implicit none

    ! jd: Arguments
    type(MODEL), intent(inout) :: mdl

    ! jd: Local variables
    integer :: isp
    integer :: iat
    integer :: z
    integer :: nfunctions
    character(len=*), parameter :: myself = 'dftb_adjust_num_ngwfs'

    ! -------------------------------------------------------------------------

    call utils_assert(mdl%nsub == 1, &
       myself//': DFTB does not support embedding regions')

    if(pub_on_root) then
       write(stdout,'(a)') '+-----------------------------------------------&
            &-----------------------------+'
       write(stdout,'(a)') '|                  DFTB basis information for ON&
            &ETEP species                 |'
       write(stdout,'(a)') '+-----------------------------------------------&
            &-----------------------------+'
       write(stdout,'(a)') '| species | symbol | atomic number | DFTB valenc&
            &e orbitals | Inferred #NGWFs |'

       do isp = 1, mdl%regions(1)%par%num_species
          z = mdl%regions(1)%species(isp)%atomic_number
          call utils_assert(z > 0 .and. z < dftb_max_species, &
               myself//': Species atomic number out of range for DFTB: '//&
               trim(utils_int_to_str(z))//'.')

          nfunctions = internal_ao_string_to_nfunctions(mdl%dftb_par%z(z)%ao)

          ! jd: Store the new #NGWFs in species
          mdl%regions(1)%species(isp)%nfunctions = nfunctions
          mdl%species(isp)%nfunctions = nfunctions

          ! jd: Store the new #NGWFs in all elements that match this species
          do iat = 1, mdl%regions(1)%par%nat
             if(mdl%regions(1)%elements(iat)%species_id == &
                  mdl%regions(1)%species(isp)%species_id) then
                mdl%regions(1)%elements(iat)%nfunctions = nfunctions
                mdl%elements(iat)%nfunctions = nfunctions
             end if
          end do

          write(stdout,'(a,i7,a,a6,a,i13,a,a21,a,i15, a)') '| ', isp, ' | ', &
               adjustr(mdl%regions(1)%species(isp)%symbol), ' | ', z, ' | ', &
               adjustr(trim(mdl%dftb_par%z(z)%ao)), ' | ', nfunctions, ' |'

       end do

       write(stdout,'(a)') '+-----------------------------------------------&
            &-----------------------------+'
    end if

    ! jd: Broadcast the changes to non-root nodes
    do isp = 1, mdl%regions(1)%par%num_species
       call comms_bcast(pub_root_proc_id, &
            mdl%regions(1)%species(isp)%nfunctions)
       call comms_bcast(pub_root_proc_id, mdl%species(isp)%nfunctions)

       do iat = 1, mdl%regions(1)%par%nat
          if(mdl%regions(1)%elements(iat)%species_id == &
               mdl%regions(1)%species(isp)%species_id) then
             call comms_bcast(pub_root_proc_id, &
                  mdl%regions(1)%elements(iat)%nfunctions)
             call comms_bcast(pub_root_proc_id, mdl%elements(iat)%nfunctions)
          end if
       end do

    end do

    return

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains

    integer function internal_ao_string_to_nfunctions(ao_string)
      !========================================================================!
      ! Returns the number of NGWFs suitable for an AO basis string used by    !
      ! xTB, like "5s5p4d".                                                    !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   ao_string (in): Per-species parameter of xTB methods describing the  !
      !                   AO's used in the basis set.                          !
      !------------------------------------------------------------------------!
      ! Return value:                                                          !
      !   Number of NGWFs in a minimal set corresponding to ao_string.         !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic on 2021/06/26.                               !
      !========================================================================!

      use dftb_parameters, only: dftb_strlen
      use rundat, only: pub_dftb_cartesian_ngwfs
      use utils, only: utils_abort

      implicit none

      ! jd: Arguments
      character(len=dftb_strlen), intent(in) :: ao_string

      ! jd: Local variables
      integer :: ci
      integer :: nfunctions
      character :: c
      character(len=*), parameter :: myself = &
           'dftb_adjust_num_ngwfs::internal_ao_string_to_nfunctions'

      ! -----------------------------------------------------------------------

      nfunctions = 0

      do ci = 1, dftb_strlen
         c = ao_string(ci:ci) !ab: correct c
         select case(c)
         case('1','2','3','4','5','6','7','8', ' ')
            ! jd: no-op for shell numbers and trailing space
         case('s')
            nfunctions = nfunctions + 1
         case('p')
            nfunctions = nfunctions + 3
         case('d')
            if (pub_dftb_cartesian_ngwfs) then
               nfunctions = nfunctions + 6
            else
               nfunctions = nfunctions + 5
            end if
         case default
            call utils_abort(myself//': Unrecognised subshell "'//c//&
                 '" in DFTB AO string "'//trim(ao_string)//'".')
         end select
      end do

      internal_ao_string_to_nfunctions = nfunctions

    end function internal_ao_string_to_nfunctions

  end subroutine dftb_adjust_num_ngwfs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_infer_params(mdl)
    !========================================================================!
    ! Infers angular momentum array and ngwf_type array for NGWFs from an AO !
    ! basis string read from the input file like "5s5p4d".                   !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    ! mdl(inout) : mdl%dftb_par is modified to include inferred parameters.  !
    !------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in July 2021.                              !
    !========================================================================!

    use dftb_parameters, only: dftb_max_species, dftb_strlen
    use model_type, only: MODEL
    use utils, only: utils_alloc_check, utils_abort, utils_int_to_str, &
         utils_assert
    use rundat, only: pub_dftb_cartesian_ngwfs

    implicit none

    ! ab: Arguments
    type(MODEL), intent(inout) :: mdl

    ! ab: Local variables
    integer :: isp, zsp
    integer :: ifun, nfun, istr, ierr, ngwf_type
    character :: string
    character(len=dftb_strlen) :: ao_string
    character(len=*), parameter :: myself = 'dftb_infer_params'

    ! -----------------------------------------------------------------------

    call utils_assert(mdl%nsub == 1, &
       myself//': DFTB does not support embedding regions')

    do isp = 1, mdl%regions(1)%par%num_species
       zsp = mdl%regions(1)%species(isp)%atomic_number
       call utils_assert(zsp > 0 .and. zsp < dftb_max_species, &
            myself//': Species atomic number out of range for DFTB: '//&
            trim(utils_int_to_str(zsp))//'.')

       ao_string = mdl%dftb_par%z(zsp)%ao
       ifun = 1
       nfun = mdl%regions(1)%species(isp)%nfunctions
       allocate(mdl%dftb_par%z(zsp)%ang_mom(nfun), stat=ierr)
       call utils_alloc_check(myself, 'mdl%dftb_par%z(zsp)%ang_mom', ierr)

       ngwf_type = 1
       allocate(mdl%dftb_par%z(zsp)%ngwf_type(nfun), stat=ierr)
       call utils_alloc_check(myself, 'mdl%dftb_par%z(zsp)%ngwf_type', ierr)

       do istr = 1, dftb_strlen
          string = ao_string(istr:istr)
          select case(string)
          case('1','2','3','4','5','6','7','8', ' ')
             ! jd: no-op for shell numbers and trailing space
          case('s')
             mdl%dftb_par%z(zsp)%ang_mom(ifun)=0
             mdl%dftb_par%z(zsp)%ngwf_type(ifun)=ngwf_type
             ifun=ifun+1
             ngwf_type=ngwf_type+1
          case('p')
             mdl%dftb_par%z(zsp)%ang_mom(ifun:ifun+2)=1
             mdl%dftb_par%z(zsp)%ngwf_type(ifun:ifun+2)=ngwf_type
             ifun=ifun+3
             ngwf_type=ngwf_type+1
          case('d')
             if (pub_dftb_cartesian_ngwfs) then
                mdl%dftb_par%z(zsp)%ang_mom(ifun:ifun+5)=2
                mdl%dftb_par%z(zsp)%ngwf_type(ifun:ifun+5)=ngwf_type
                ifun=ifun+6
             else
                mdl%dftb_par%z(zsp)%ang_mom(ifun:ifun+4)=2
                mdl%dftb_par%z(zsp)%ngwf_type(ifun:ifun+4)=ngwf_type
                ifun=ifun+5
             end if
             ngwf_type=ngwf_type+1
          case default
             call utils_abort(myself//': Unrecognised subshell "'//string//&
                  '" in DFTB AO string "'//trim(ao_string)//'".')
          end select
          if (ifun > nfun) exit
       end do
    end do

  end subroutine dftb_infer_params

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_calculate_overlap(mdl, rep, ngwf_basis)
    !==========================================================================!
    ! This subroutine calculates the overlap matrix from NGWFs.                !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   mdl(in)    : mdl%cell is needed for calculating the overlap matrix.    !
    !   rep(inout) : rep%ngwfs_on_grid are the input ngwfs and rep%overlap is  !
    !                the resultant output overlap matrix.                      !
    !   ngwf_basis : NGWF basis set. Both PAOs and STO-3G are allowed.         !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in June 2021                                 !
    !==========================================================================!

    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketppd
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_dftb_overlap_analytical
    use sparse, only: sparse_show_matrix

    implicit none

    ! ab: Arguments
    type(MODEL), intent(in)        :: mdl
    type(NGWF_REP), intent(inout)  :: rep
    type(FUNC_BASIS), intent(in)   :: ngwf_basis(rep%nsub)

    ! ab: Local variables
    integer :: isub, jsub, ncols, mrows
    character(len=*), parameter :: myself = 'dftb_calculate_overlap'

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering '//myself

    ncols = rep%overlap%ncols
    mrows = rep%overlap%mrows

    ! ab: calculate bracket of ngwfs
    do isub=1,ncols
       do jsub=1,mrows
          if(pub_dftb_overlap_analytical) then
             call dftb_analytical_gaussian_overlap(rep%overlap%m(jsub,isub), &
                  mdl, ngwf_basis(1))
          else
             call function_ops_brappd_ketppd(rep%overlap%m(jsub,isub), &
                  rep%ngwfs_on_grid(jsub), ngwf_basis(jsub), &
                  rep%ngwfs_on_grid(isub), ngwf_basis(isub), mdl%cell)
          end if
       end do
    end do

! jd: For debugging only
#if 0
    write (stdout, '(a)') "Overlap matrix"
    call sparse_show_matrix(rep%overlap%m(1,1), show_elems=.true.)
#endif

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving '//myself

  end subroutine dftb_calculate_overlap

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_calculate_ham(mdl, ham, rep, ngwf_basis)
    !==========================================================================!
    ! This subroutine calculates the hamiltonian matrix.                       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   mdl(in)        : mdl is used for obtaining dftb_par, elements and cell !
    !   ham(inout)     : hamiltonian matrix                                    !
    !   rep(in)        : rep%overlap matrix is used in the calculation of h0ij !
    !   ngwf_basis(in) : NGWF basis set. Both PAOs and STO-3G are allowed.     !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Jan 2022.                                 !
    !==========================================================================!

    use constants, only: DFTB_GFN0
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_HAM, NGWF_REP
    use rundat, only: pub_dftb_method
    use utils, only: utils_abort

    implicit none

    ! ab: Arguments
    type(MODEL), intent(in)        :: mdl
    type(NGWF_HAM), intent(inout)  :: ham
    type(NGWF_REP), intent(in)     :: rep
    type(FUNC_BASIS), intent(in)   :: ngwf_basis(rep%nsub)

    ! ab: Local variables
    character(len=*), parameter    :: myself='dftb_calculate_ham'

    ! -------------------------------------------------------------------------

    select case (pub_dftb_method)
    case(DFTB_GFN0)
       call dftb_calculate_h0ij(mdl, ham, rep, ngwf_basis)
    case default
       call utils_abort(myself//': Unsupported DFTB method. &
            &Only GFN0 is supported.')
    end select

  end subroutine dftb_calculate_ham

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_calculate_h0ij(mdl, ham, rep, ngwf_basis)
    !==========================================================================!
    ! This subroutine calculates the zeroth order hamiltonian h0 from [5]:(23) !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   mdl(in)        : mdl is used for obtaining dftb_par, elements and cell !
    !   ham(inout)     : hamiltonian matrix                                    !
    !   rep(in)        : rep%overlap matrix is used in the calculation of h0ij !
    !   ngwf_basis(in) : NGWF basis set. Both PAOs and STO-3G are allowed.     !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in July 2021                                 !
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, magnitude, OPERATOR(-)
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_HAM, NGWF_REP
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_num_spins, pub_dftb_bc_is_periodic
    use simulation_cell, only: minimum_image_displacement
    use sparse, only: sparse_put_element, sparse_get_element, &
         sparse_get_par, sparse_show_matrix
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! ab: Arguments
    type(MODEL), intent(in)        :: mdl
    type(NGWF_HAM), intent(inout)  :: ham
    type(NGWF_REP), intent(in)     :: rep
    type(FUNC_BASIS), intent(in)   :: ngwf_basis(rep%nsub)

    ! ab: Local variables
    character(len=*), parameter    :: myself='dftb_calculate_h0ij'
    type(PARAL_INFO), pointer      :: par

    ! Atom and NGWF counters
    integer           :: local_i, global_i, first_ngwf_of_i, nfun_i, orig_i
    integer           :: local_j, global_j, first_ngwf_of_j, nfun_j, orig_j
    integer           :: zi, zj, li, lj, ngwf_i, ngwf_j, ni, nj
    integer           :: global_ii_ngwf, global_jj_ngwf
    real(kind=DP)     :: aij, kij, pij, xij, yij
    real(kind=DP)     :: sij, h0ij, d_ij
    type(POINT)       :: ri, rj, r_ij
    ! --------------------------------------------------------------------------

    call timer_clock(myself,1)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering '//myself

    call sparse_get_par(par, rep%overlap%m(1,1))
    ! --------------------------------------------------------------------------
    ! ab: Loop over all atoms on this core                               III
    ! --------------------------------------------------------------------------
    loop_I:                                                                   &
    do local_i = 1, par%num_atoms_on_proc(pub_my_proc_id)
       global_i = par%first_atom_on_proc(pub_my_proc_id) + local_i-1
       first_ngwf_of_i = ngwf_basis(1)%first_on_atom(global_i)
       orig_i = par%orig_atom(global_i)
       zi = mdl%elements(orig_i)%atomic_number
       ri = mdl%elements(orig_i)%centre
       nfun_i= mdl%elements(orig_i)%nfunctions
       ! -----------------------------------------------------------------------
       ! ab: Loop over all atoms J that are s-neighbours with I          JJJ
       ! -----------------------------------------------------------------------
       loop_J:                                                                &
       do local_j = mdl%dftb_var%nl%first_idx(global_i),mdl%dftb_var%nl%last_idx(global_i)
          global_j = mdl%dftb_var%nl%neighbours(local_j)
          first_ngwf_of_j = ngwf_basis(1)%first_on_atom(global_j)
          orig_j = par%orig_atom(global_j)
          zj = mdl%elements(orig_j)%atomic_number
          rj = mdl%elements(orig_j)%centre
          nfun_j = mdl%elements(orig_j)%nfunctions
          ! ab: Take care of boundary conditions
          if (.not. any(pub_dftb_bc_is_periodic)) then
             ! ab: in OBCs
             r_ij = ri - rj
             d_ij = magnitude(r_ij)
          else if (all(pub_dftb_bc_is_periodic)) then
             ! ab: find the minimum image of atom I
             call minimum_image_displacement(d_ij, r_ij, ri, rj, mdl%cell)
          else
             call utils_abort(myself//': Mixed BCs are not currently supported.')
          end if
          ! --------------------------------------------------------------------
          ! ab: loop over all ngwfs i on atom I                          iii
          ! --------------------------------------------------------------------
          loop_ngwf_i:                                                        &
          do ngwf_i = 1, nfun_i
             global_ii_ngwf = first_ngwf_of_i + ngwf_i - 1
             li = mdl%dftb_par%z(zi)%ang_mom(ngwf_i)
             ni = mdl%dftb_par%z(zi)%ngwf_type(ngwf_i)
             ! -----------------------------------------------------------------
             ! ab: loop over all ngwfs j on J                            jjj
             ! -----------------------------------------------------------------
             loop_ngwf_j:                                                     &
             do ngwf_j = 1, nfun_j
                global_jj_ngwf = first_ngwf_of_j + ngwf_j - 1
                lj = mdl%dftb_par%z(zj)%ang_mom(ngwf_j)
                nj = mdl%dftb_par%z(zj)%ngwf_type(ngwf_j)

                ! --------------------------------------------------------------
                ! ab: Diagonal values obtained from dftb_calculate_h0ii
                ! off diagonal terms calculated as:
                ! h0ij = aij * kij * pij * xij * yij * sij
                ! --------------------------------------------------------------
                if (global_i == global_j) then

                   if (ngwf_i == ngwf_j) then
                      h0ij = dftb_calculate_h0ii(orig_i, li, ni, mdl%dftb_par, &
                           mdl%elements)
                   else
                      h0ij = 0.0_DP
                   end if

                else

                   aij = dftb_calculate_aij(orig_i, orig_j, li, lj, ni, nj, &
                        mdl%dftb_par, mdl%elements)

                   kij = dftb_calculate_kij(zi, zj, li, lj, mdl%dftb_par, &
                        vi = ni.eq. mdl%dftb_par%z(zi)%ngwf_type(nfun_i), &
                        vj = nj.eq. mdl%dftb_par%z(zj)%ngwf_type(nfun_j))

                   pij = dftb_calculate_pij(zi, zj, li, lj, d_ij, &
                        mdl%dftb_par)

                   xij = dftb_calculate_xij(zi, zj, li, lj, mdl%dftb_par, &
                        vi = ni.eq. mdl%dftb_par%z(zi)%ngwf_type(nfun_i), &
                        vj = nj.eq. mdl%dftb_par%z(zj)%ngwf_type(nfun_j))

                   yij = dftb_calculate_yij(zi, zj, ni, nj, mdl%dftb_par)

                   call sparse_get_element(sij, rep%overlap%m(1,1), &
                        global_jj_ngwf, global_ii_ngwf)

                   h0ij = aij * kij * xij * yij * pij * sij

                end if

                ! ab: put the h0ij element in the hamiltonian matrix
                call sparse_put_element(h0ij, ham%ham(1)%m(1,1), &
                     global_jj_ngwf, global_ii_ngwf)

             end do loop_ngwf_j
          end do loop_ngwf_i
       end do loop_J
    end do loop_I

    ! ab: spin degenerate matrix
    if(pub_num_spins==2) ham%ham(2) = ham%ham(1)

! ab: For debugging only
#if 0
    write (stdout, '(a)') "Hamiltonian matrix"
    call sparse_show_matrix(ham%ham(1)%m(1,1), show_elems=.true.)
#endif

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving '//myself

    call timer_clock(myself,2)

  end subroutine dftb_calculate_h0ij

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_electronic_force(ngwf_basis,dk,mdl,rep,ev,oc,mo,num,F_dftb_el)
    !==========================================================================!
    ! This subroutine calculates the electronic forces.                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   ngwf_basis(in) : NGWF basis set. Both PAOs and STO-3G are allowed.     !
    !   dk(in)         : density kernel                                        !
    !   mdl(in)        : mdl is used for obtaining dftb_par, elements and cell !
    !   rep(in)        : rep%overlap matrix is used in the calculation.        !
    !   ev(in)         : eigenvalues of the hamiltonian matrix.                !
    !   oc(in)         : fractional occupancies of orbitals.                   !
    !   mo(in): rotation from NGWFs to molecular orbitals basis set(M^\alpha_i)!
    !   num(in)        : number of molecular orbitals.                         !
    !   F_dftb_el(out) : resultant electronic force.                           !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in July 2022.                                !
    !==========================================================================!

    use dense, only: DEM, dense_create, dense_product, dense_destroy, &
         dense_convert, dense_put_element, dense_get_col, dense_put_col
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_num_spins, PUB_1K, pub_spin_fac
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
         sparse_embed_create, sparse_embed_destroy
    use timer, only: timer_clock

    implicit none

    ! ab: Arguments
    type(FUNC_BASIS), intent(in)     :: ngwf_basis
    type(SPAM3_EMBED_ARRAY),intent(in) :: dk
    type(MODEL), intent(in)          :: mdl
    type(NGWF_REP), intent(in)       :: rep
    real(kind=DP), intent(in)        :: ev(:,:)
    real(kind=DP), intent(in)        :: oc(:,:)
    type(DEM),     intent(in)        :: mo(:)
    integer, intent(in)              :: num
    real(kind=DP), intent(out)       :: F_dftb_el(1:3,1:mdl%nat)

    ! ab: Local variables
    integer       :: is, iorb, ic, nat
    type(DEM)     :: ev_dens, mo_eval, wk_dens, dcn(3), dfc, dfq, fq, fc
    real(kind=DP) :: wo, dfdc(mdl%nat+1), dfdq(mdl%nat+1)
    type(SPAM3_EMBED) :: wk(1:pub_num_spins)
    character(len=*), parameter    :: myself='dftb_electronic_force'
    ! --------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering '//myself

    call timer_clock(myself,1)

    ! ab: initialize the forces and create matrices
    F_dftb_el = 0.0_DP
    dfdc = 0.0_DP
    dfdq = 0.0_DP
    nat = mdl%nat
    call dense_create(fc, nat+1, 1, iscmplx=.false.)
    call dense_create(fq, nat+1, 1, iscmplx=.false.)
    call dense_create(dfc, nat+1, 1, iscmplx=.false.)
    call dense_create(dfq, nat+1, 1, iscmplx=.false.)
    call dense_create(mo_eval, num, num, iscmplx=.false.)
    call dense_create(ev_dens, num, num, iscmplx=.false.)
    call dense_create(wk_dens, num, num, iscmplx=.false.)

    ! ab: create energy weighted density kernel
    do is = 1, pub_num_spins
       call sparse_embed_create(wk(is),dk%m(is,PUB_1K))
       do iorb = 1, num
          wo = oc(iorb,is) * ev(iorb,is)
          call dense_put_element(wo, ev_dens, iorb, iorb)
       end do
       call dense_product(mo_eval, mo(is), ev_dens, &
            opA='N',opB='N',first_k=1,last_k=num)
       call dense_product(wk_dens, mo_eval, mo(is), &
            opA='N',opB='T',first_k=1,last_k=num)
       call dense_convert(wk(is), wk_dens)
    end do

    ! ab: Calculate forces on each atom by breaking into two-body & three-body
    ! ab: two-body forces
    call dftb_deldr(dk,wk,ngwf_basis,mdl,rep,F_dftb_el,dfdc(1:nat),dfdq(1:nat))
    ! ab: three body forces
    call dense_put_col(dfdc, dfc, 1)
    call dense_put_col(dfdq, dfq, 1)
    do ic = 1, 3
       call dense_product(fq, mdl%dftb_var%dqdr(ic), dfq, opA='T', opB='N')
       call dense_product(fc, mdl%dftb_var%dcdr(ic), dfc, opA='T', opB='N')
       call dense_get_col(dfdq, fq, 1)
       call dense_get_col(dfdc, fc, 1)
       F_dftb_el(ic,:) = F_dftb_el(ic,:) + dfdq(1:nat) + dfdc(1:nat)
    end do

    F_dftb_el = F_dftb_el * pub_spin_fac

    ! ab: deallocate the matrices
    call dense_destroy(mo_eval)
    call dense_destroy(ev_dens)
    call dense_destroy(wk_dens)
    call dense_destroy(fc)
    call dense_destroy(fq)
    call dense_destroy(dfc)
    call dense_destroy(dfq)

    do is = 1, pub_num_spins
       call sparse_embed_destroy(wk(is))
    end do

    call timer_clock(myself,2)

  end subroutine dftb_electronic_force

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_deldr(dk,wk,ngwf_basis,mdl,rep,force,dfdc,dfdq)
    !==========================================================================!
    ! This subroutines returns the two-body force and its partial derviatives. !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   dk(in)         : density kernel.                                       !
    !   wk(in)         : energy-weighted density kernel.                       !
    !   ngwf_basis(in) : NGWF basis set. Both PAOs and STO-3G are allowed.     !
    !   mdl(in)        : mdl is used for obtaining dftb_par, elements and cell !
    !   rep(in)        : rep%overlap matrix is used in the calculation.        !
    !   force(out)     : two body force.                                       !
    !   dfdc(out) : derivative of the force with respect to coordination number!
    !   dfdq(out) : derivative of the force with respect to charges.           !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in July 2022.                                !
    !==========================================================================!

    use comms, only: pub_my_proc_id, comms_reduce
    use constants, only: HARTREE_IN_EVS
    use function_basis, only: FUNC_BASIS
    use gaussian_ops, only: gaussian_overlap_integral
    use geometry, only: POINT, magnitude, OPERATOR(-), OPERATOR(+)
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_dftb_bc_is_periodic, pub_dftb_overlap_cutoff, &
         PUB_1K, pub_num_spins, pub_threads_max
    use simulation_cell, only: minimum_image_displacement
    use sparse, only: sparse_get_element, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_abort
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY

    implicit none

    ! ab: Arguments
    type(SPAM3_EMBED_ARRAY), intent(in) :: dk
    type(SPAM3_EMBED), intent(in)    :: wk(1:pub_num_spins)
    type(FUNC_BASIS), intent(in)     :: ngwf_basis
    type(MODEL), intent(in)          :: mdl
    type(NGWF_REP), intent(in)       :: rep
    real(kind=DP), intent(inout)     :: force(1:3,mdl%nat)
    real(kind=DP), intent(inout)     :: dfdc(mdl%nat), dfdq(mdl%nat)

    ! ab: Local variables
    character(len=*), parameter    :: myself='dftb_deldr'

    ! Atom and NGWF counters
    integer           :: local_i, global_i, first_ngwf_of_i, nfun_i, orig_i
    integer           :: local_j, global_j, first_ngwf_of_j, nfun_j, orig_j
    integer           :: zi, zj, li, lj, ngwf_i, ngwf_j, ni, nj, is
    integer           :: global_ii_ngwf, global_jj_ngwf, nshells_i, nshells_j
    integer           :: s_i, s_j, iorb
    integer           :: nprim_s_i, nprim_s_j, prim_i, prim_j
    integer           :: angmom_s_i, angmom_s_j
    integer           :: ix, iy, iz, jx, jy, jz
    real(kind=DP)     :: aij, kij, pij, xij, yij, a_i, a_j, c_i, c_j, cij
    real(kind=DP)     :: sij, d_ij, mk, pk, qi, qj, kqi, kqj, dh(1:3), bij
    real(kind=DP)     :: dsx, dsy, dsz, dsx1, dsx2, dsy1, dsy2, dsz1, dsz2
    real(kind=DP)     :: kcni, kcnj, kq2i, kq2j, df(1:3), dpij(1:3), ds(1:3)
    type(POINT)       :: ri, rj, r_ij
    type(PARAL_INFO), pointer :: par
    ! --------------------------------------------------------------------------

    call timer_clock(myself,1)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering '//myself

    call sparse_get_par(par, rep%overlap%m(1,1))

    ! -----------------------------------------------------------------------
    ! ab: Loop over all atoms on this core                               III
    ! -----------------------------------------------------------------------
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(local_i, global_i, first_ngwf_of_i, orig_i, zi, ri, nfun_i, &
!$OMP      nshells_i, ngwf_i, s_i, nprim_s_i, angmom_s_i, ix, iy, iz, &
!$OMP      global_ii_ngwf, li, ni, local_j, global_j, first_ngwf_of_j, orig_j, &
!$OMP      zj, rj, nfun_j, nshells_j, ngwf_j, r_ij, d_ij, s_j, nprim_s_j, &
!$OMP      angmom_s_j, jx, jy, jz, sij, global_jj_ngwf, lj, nj, dsx, dsy, dsz, &
!$OMP      a_i, c_i, a_j, c_j, dsx1, dsy1, dsz1, dsx2, ds, &
!$OMP      kcni, kcnj, kqi, kqj, dh, dpij, df, bij, cij, qi, kq2i, qj, kq2j, &
!$OMP      dsy2, dsz2, kij, xij, yij, pij, aij, pk, mk) &
!$OMP SHARED(par, mdl, pub_my_proc_id, ngwf_basis, dk, rep, &
!$OMP      pub_dftb_bc_is_periodic, pub_dftb_overlap_cutoff, pub_num_spins, wk)&
!$OMP REDUCTION(+:force, dfdc, dfdq)
    loop_I:                                                                   &
    do local_i = 1, par%num_atoms_on_proc(pub_my_proc_id)
       global_i = par%first_atom_on_proc(pub_my_proc_id) + local_i-1
       first_ngwf_of_i = ngwf_basis%first_on_atom(global_i)
       orig_i = par%orig_atom(global_i)
       zi = mdl%elements(orig_i)%atomic_number
       ri = mdl%elements(orig_i)%centre
       qi = mdl%elements(orig_i)%dftb_charge
       kq2i = mdl%dftb_par%z(zi)%kqat2 / HARTREE_IN_EVS
       nfun_i= mdl%elements(orig_i)%nfunctions
       nshells_i = mdl%dftb_par%z(zi)%nshells

       ngwf_i = 0
       ! --------------------------------------------------------------------
       ! ab: Loop over all shells s_i on atom I                          s_i
       ! --------------------------------------------------------------------
       loop_si:                                                             &
       do s_i = 1, nshells_i
          nprim_s_i = mdl%dftb_par%z(zi)%sto(s_i)%nprim
          angmom_s_i = mdl%dftb_par%z(zi)%ang(s_i)
          !-----------------------------------------------------------------
          ! ab: Loop over all orbitals in s_i with angular moment ix, iy, iz
          !-----------------------------------------------------------------
          loop_ix:                                                         &
          do ix = 0, angmom_s_i
             loop_iy:                                                      &
             do iy = 0, angmom_s_i
                loop_iz:                                                   &
                do iz = 0, angmom_s_i

                   if (ix+iy+iz .ne. angmom_s_i) cycle

                   ! ab: ngwf counter
                   ngwf_i = ngwf_i + 1
                   global_ii_ngwf = first_ngwf_of_i + ngwf_i - 1
                   li = mdl%dftb_par%z(zi)%ang_mom(ngwf_i)
                   ni = mdl%dftb_par%z(zi)%ngwf_type(ngwf_i)
                   kcni = dftb_kcn(zi, li, ni, mdl%dftb_par)
                   kqi = dftb_kq(zi,li,mdl%dftb_par)+2.0_DP*qi*kq2i

                   ! -------------------------------------------------------------------
                   ! ab: Loop over all atoms J that are s-neighbours with I          JJJ
                   ! -------------------------------------------------------------------
                   loop_J:                                                                &
                   do local_j = mdl%dftb_var%nl%first_idx(global_i), mdl%dftb_var%nl%last_idx(global_i)
                      global_j = mdl%dftb_var%nl%neighbours(local_j)
                      first_ngwf_of_j = ngwf_basis%first_on_atom(global_j)
                      orig_j = par%orig_atom(global_j)
                      zj = mdl%elements(orig_j)%atomic_number
                      rj = mdl%elements(orig_j)%centre
                      qj = mdl%elements(orig_j)%dftb_charge
                      kq2j = mdl%dftb_par%z(zj)%kqat2 / HARTREE_IN_EVS
                      nfun_j = mdl%elements(orig_j)%nfunctions
                      nshells_j = mdl%dftb_par%z(zj)%nshells

                      ! ab: Take care of boundary conditions
                      if (.not. any(pub_dftb_bc_is_periodic)) then
                         ! ab: in OBCs
                         r_ij = ri - rj
                         d_ij = magnitude(r_ij)
                      else if (all(pub_dftb_bc_is_periodic)) then
                         ! ab: find the minimum image of atom I
                         call minimum_image_displacement(d_ij, r_ij, ri, rj, mdl%cell)
                         ri = rj + r_ij
                      else
                         call utils_abort(myself//': Mixed BCs are not currently supported.')
                      end if

                      ngwf_j = 0
                      ! -------------------------------------------------------------
                      ! jd: Loop over all shells s_j on atom J                    s_j
                      ! -------------------------------------------------------------
                      loop_sj:                                                       &
                      do s_j = 1, nshells_j
                         nprim_s_j = mdl%dftb_par%z(zj)%sto(s_j)%nprim
                         angmom_s_j = mdl%dftb_par%z(zj)%ang(s_j)
                         !---------------------------------------------------
                         ! ab: Loop over all orbitals in s_j with  jx, jy, jz
                         !---------------------------------------------------
                         loop_jx:                                           &
                         do jx = 0, angmom_s_j
                            loop_jy:                                        &
                            do jy = 0, angmom_s_j
                               loop_jz:                                     &
                               do jz = 0, angmom_s_j

                                  if (jx+jy+jz .ne. angmom_s_j) cycle

                                  ! ab: ngwf counter
                                  ngwf_j = ngwf_j + 1
                                  global_jj_ngwf = first_ngwf_of_j + ngwf_j - 1
                                  lj = mdl%dftb_par%z(zj)%ang_mom(ngwf_j)
                                  nj = mdl%dftb_par%z(zj)%ngwf_type(ngwf_j)
                                  kcnj = dftb_kcn(zj, lj, nj, mdl%dftb_par)
                                  kqj = dftb_kq(zj,lj,mdl%dftb_par)+2.0_DP*qj*kq2j

                                  ! ab: overlap element
                                  call sparse_get_element(sij, rep%overlap%m(1,1), &
                                       global_jj_ngwf, global_ii_ngwf)

                                  ! ab: initialize overlap derivatives to be zero
                                  dsx = 0.0_DP
                                  dsy = 0.0_DP
                                  dsz = 0.0_DP

                                  ! ab: calculate only if distance less than cutoff
                                  if (d_ij .le. pub_dftb_overlap_cutoff) then
                                     !---------------------------------------------
                                     ! ab: Loop over all primitives          prim_i
                                     !---------------------------------------------
                                     loop_prim_i:                                 &
                                     do prim_i = 1, nprim_s_i

                                        ! ab: exponent
                                        a_i= mdl%dftb_par%z(zi)%sto(s_i)%&
                                             alpha(prim_i)
                                        ! ab: contraction coefficient
                                        c_i= mdl%dftb_par%z(zi)%sto(s_i)%&
                                             coeff(prim_i)

                                        !------------------------------------------
                                        ! ab: Loop over all primitives       prim_j
                                        !------------------------------------------
                                        loop_prim_j:                              &
                                        do prim_j = 1, nprim_s_j

                                           ! ab: exponent
                                           a_j= mdl%dftb_par%z(zj)%sto(s_j)%&
                                                alpha(prim_j)
                                           ! ab: contraction coefficient
                                           c_j= mdl%dftb_par%z(zj)%sto(s_j)%&
                                                coeff(prim_j)

                                           ! ab: derivatives
                                           dsx1 = gaussian_overlap_integral( &
                                             a_i,a_j,ri,rj,ix+1,jx,iy,jy,iz,jz)
                                           dsy1 = gaussian_overlap_integral( &
                                             a_i,a_j,ri,rj,ix,jx,iy+1,jy,iz,jz)
                                           dsz1 = gaussian_overlap_integral( &
                                             a_i,a_j,ri,rj,ix,jx,iy,jy,iz+1,jz)
                                           if (ix .ge. 1) then
                                              dsx2 = gaussian_overlap_integral( &
                                                a_i,a_j,ri,rj,ix-1,jx,iy,jy,iz,jz)
                                           else
                                              dsx2 = 0.0_DP
                                           end if
                                           if (iy .ge. 1) then
                                              dsy2 = gaussian_overlap_integral( &
                                                a_i,a_j,ri,rj,ix,jx,iy-1,jy,iz,jz)
                                           else
                                              dsy2 = 0.0_DP
                                           end if
                                           if (iz .ge. 1) then
                                              dsz2 = gaussian_overlap_integral( &
                                                a_i,a_j,ri,rj,ix,jx,iy,jy,iz-1,jz)
                                           else
                                              dsz2 = 0.0_DP
                                           end if

                                           dsx=dsx+c_i*c_j*(2.0_DP*a_i*dsx1-ix*dsx2)
                                           dsy=dsy+c_i*c_j*(2.0_DP*a_i*dsy1-iy*dsy2)
                                           dsz=dsz+c_i*c_j*(2.0_DP*a_i*dsz1-iz*dsz2)

                                        end do loop_prim_j

                                     end do loop_prim_i

                                  end if
                                  ds = (/ dsx, dsy, dsz/)

                                  ! ab: orbitals on the same atom
                                  if (orig_i == orig_j) then
                                     ! ab: same orbital
                                     if (ngwf_i == ngwf_j) then
                                        loop_spin:                                   &
                                        do is = 1, pub_num_spins

                                           call sparse_get_element(pk,dk%m(is,PUB_1K)%m(1,1), &
                                                global_jj_ngwf, global_ii_ngwf)
                                           dfdc(orig_i) = dfdc(orig_i) + pk*kcni
                                           dfdq(orig_i) = dfdq(orig_i) + pk*kqi

                                        end do loop_spin
                                     end if

                                  ! ab: orbitals on different atoms
                                  else

                                     aij = dftb_calculate_aij(orig_i, orig_j, li, lj,&
                                           ni, nj, mdl%dftb_par, mdl%elements)

                                     kij = dftb_calculate_kij(zi, zj, li, lj, mdl%dftb_par, &
                                          vi = ni.eq. mdl%dftb_par%z(zi)%ngwf_type(nfun_i), &
                                          vj = nj.eq. mdl%dftb_par%z(zj)%ngwf_type(nfun_j))

                                     xij = dftb_calculate_xij(zi, zj, li, lj, mdl%dftb_par, &
                                          vi = ni.eq. mdl%dftb_par%z(zi)%ngwf_type(nfun_i), &
                                          vj = nj.eq. mdl%dftb_par%z(zj)%ngwf_type(nfun_j))

                                     yij = dftb_calculate_yij(zi, zj, ni, nj, mdl%dftb_par)

                                     ! ab: shell polynomial
                                     pij = dftb_calculate_pij(zi,zj,li,lj,d_ij,&
                                          mdl%dftb_par)
                                     dpij = dftb_calculate_dpij(zi,zj,li,lj,r_ij, &
                                          mdl%dftb_par)

                                     bij = kij*xij*yij
                                     dh = aij*bij*(dpij*sij + pij*ds)

                                     loop_spins:                                   &
                                     do is = 1, pub_num_spins

                                        call sparse_get_element(pk,dk%m(is,PUB_1K)%m(1,1), &
                                             global_jj_ngwf, global_ii_ngwf)
                                        call sparse_get_element(mk, wk(is)%m(1,1), &
                                             global_jj_ngwf, global_ii_ngwf)

                                        df = mk*ds - pk*dh
                                        force(:,orig_i) = force(:,orig_i) + df
                                        force(:,orig_j) = force(:,orig_j) - df

                                        cij = 0.5_DP*pk*bij*pij*sij
                                        dfdc(orig_i) = dfdc(orig_i) + cij*kcni
                                        dfdc(orig_j) = dfdc(orig_j) + cij*kcnj
                                        dfdq(orig_i) = dfdq(orig_i) + cij*kqi
                                        dfdq(orig_j) = dfdq(orig_j) + cij*kqj

                                     end do loop_spins

                                  end if

                               end do loop_jz

                            end do loop_jy

                         end do loop_jx

                      end do loop_sj

                   end do loop_J

                end do loop_iz

             end do loop_iy

          end do loop_ix

       end do loop_si

    end do loop_I

    call comms_reduce('SUM', force)
    call comms_reduce('SUM', dfdc)
    call comms_reduce('SUM', dfdq)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving '//myself

    call timer_clock(myself,2)

  end subroutine dftb_deldr

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_eht_energy(denskern, ham, eht)
    !==========================================================================!
    ! This subroutine calculates the energy via contracting the density kernel !
    ! with hamiltonian matrix.                                                 !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   denskern(in)   : density kernel                                        !
    !   ham(in)        : hamiltonian matrix                                    !
    !   eht(out)       : energy = sum_\alpha\beta K^\alpha\beta H_\beta\alpha. !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Jan 2022.                                 !
    !==========================================================================!

    use ngwf_representation, only: NGWF_HAM
    use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_trace
    use rundat, only: pub_num_spins, PUB_1K, pub_spin_fac

    implicit none

    ! ab: Arguments
    type(SPAM3_EMBED_ARRAY),intent(in) :: denskern
    type(NGWF_HAM), intent(in)         :: ham
    real(kind=DP), intent(out)         :: eht

    ! ab: Local variables
    integer :: is
    real(kind=DP) :: trace
    character(len=*), parameter :: myself='dftb_eht_energy'

    ! --------------------------------------------------------------------------

    eht = 0.0_DP
    do is=1, pub_num_spins
       call sparse_embed_trace(trace,denskern%m(is,PUB_1K),ham%ham(is))
       eht=eht+pub_spin_fac*trace
    end do

  end subroutine dftb_eht_energy

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function dftb_calculate_aij(global_iat, global_jat, li, lj, &
       ni, nj, dftb_par, elements)
    !==========================================================================!
    ! This subroutine calculates the contribution from atomic energy levels as !
    ! average of the diagonal elements of hamiltonian: aij = 1/2 (h0ii + h0jj).!
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! global_iat, global_jat (in) : global ids for corresponding atoms I and J.!
    ! li, lj (in)                 : angular momentum for functions i and j.    !
    ! ni, nj (in)                 : ngwf_type for functions i and j.           !
    ! dftb_par (in)               : dftb parameter set.                        !
    ! elements(in)                : elements array.                            !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in July 2021                                 !
    !==========================================================================!

    use dftb_parameters, only: DFTB_PARAMS
    use ion, only: ELEMENT

    implicit none

    ! ab: Arguments
    integer, intent(in)           :: global_iat, global_jat
    integer, intent(in)           :: li, lj
    integer, intent(in)           :: ni, nj
    type(DFTB_PARAMS), intent(in) :: dftb_par
    type(ELEMENT), intent(in)     :: elements(:)

    ! ab: Local variables
    character(len=*), parameter :: myself='dftb_calculate_aij'

    ! -------------------------------------------------------------------------

    dftb_calculate_aij = 0.50_DP * &
         (dftb_calculate_h0ii(global_iat, li, ni, dftb_par, elements) &
         +dftb_calculate_h0ii(global_jat, lj, nj, dftb_par, elements))

  end function dftb_calculate_aij

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function dftb_calculate_h0ii(global_iat, li, ni, dftb_par, &
       elements)
    !==========================================================================!
    ! This subroutine calculates a diagonal entry of the hamiltonian (h0ii).   !
    ! GFN0: [5]:(51), GFN1: [5]:(29), GFN2: [5]:(38)                           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! global_iat(in): global ids for atom I.                                   !
    ! li(in)        : angular momentum for function i.                         !
    ! ni(in)        : ngwf_type for functions i.                               !
    ! dftb_par (in) : dftb parameter set.                                      !
    ! elements(in)  : elements array.                                          !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in July 2021                                 !
    !==========================================================================!

    use constants, only: DFTB_GFN0, HARTREE_IN_EVS
    use dftb_parameters, only: DFTB_PARAMS
    use ion, only: ELEMENT
    use rundat, only: pub_dftb_method

    implicit none

    ! ab: Arguments
    integer, intent(in)           :: global_iat
    integer, intent(in)           :: li
    integer, intent(in)           :: ni
    type(DFTB_PARAMS), intent(in) :: dftb_par
    type(ELEMENT), intent(in)     :: elements(:)

    ! ab: Local variables
    real(kind=DP) :: hii, kcni, coordi, qi, kqi, kq2i
    integer :: zi
    character(len=*), parameter :: myself='dftb_calculate_h0ii'

    ! -------------------------------------------------------------------------

    zi = elements(global_iat)%atomic_number
    hii = dftb_par%z(zi)%lev(ni) / HARTREE_IN_EVS
    kcni = dftb_kcn(zi, li, ni, dftb_par)
    coordi = elements(global_iat)%dftb_coord
    dftb_calculate_h0ii = hii - kcni * coordi

    ! ab: include charge dependence in case of GFN0
    if (pub_dftb_method==DFTB_GFN0) then
       qi = elements(global_iat)%dftb_charge
       kqi = dftb_kq(zi, li, dftb_par)
       kq2i = dftb_par%z(zi)%kqat2 / HARTREE_IN_EVS
       dftb_calculate_h0ii = dftb_calculate_h0ii - kqi * qi - kq2i*qi**2
    end if

  end function dftb_calculate_h0ii

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function dftb_kcn(zi, li, ni, dftb_par) result(kcni)
    !==========================================================================!
    ! This function returns kcn (coefficient of coordination number).          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! zi(in)        : atomic number for atom I.                                !
    ! li(in)        : angular momentum for function i.                         !
    ! dftb_par (in) : dftb parameter set.                                      !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in May 2022.                                 !
    !==========================================================================!

    use constants, only: DFTB_GFN0, DFTB_GFN1, DFTB_GFN2, HARTREE_IN_EVS
    use dftb_parameters, only: DFTB_PARAMS
    use rundat, only: pub_dftb_method
    use utils, only: utils_abort

    implicit none

    ! ab: Arguments
    integer, intent(in)           :: zi
    integer, intent(in)           :: li
    integer, intent(in)           :: ni
    type(DFTB_PARAMS), intent(in) :: dftb_par

    ! ab: Local variables
    real(kind=DP) :: kcni, hii
    integer :: kcn_type
    character(len=*), parameter :: myself='dftb_kcn'

    integer, parameter :: gfn1kcn(118) = [& ! hard-coded element-specific data
    &  1,                                                 1, &! H-He
    &  0, 0,                               0, 1, 1, 1, 1, 1, &! Li-Ne
    &  0, 0,                               0, 1, 1, 1, 1, 1, &! Na-Ar
    &  0, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, &! K-Kr
    &  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, &! Rb-Xe
    &  0, 0, &! Cs/Ba
    &        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &!La-Lu
    &        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, &! Lu-Rn
    &  0, 0, &
    &        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &!Fr-
    &        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1 ]! -Og
    ! -------------------------------------------------------------------------

    kcni = 0.0_DP
    hii = dftb_par%z(zi)%lev(ni) / HARTREE_IN_EVS

    if(pub_dftb_method==DFTB_GFN1) then
       ! ab: values read from the parameter file are divided by 100 as per
       ! values given in Table 2 of ref. [4] and similar rescaling done in
       ! line 837 of xtb/src/xtb/gfn1.f90
       kcn_type = gfn1kcn(zi)
       if(kcn_type > 0) then
          if(li==0) then
             kcni = -hii * dftb_par%cns / 100.0_DP
          else if (li==1) then
             kcni = -hii * dftb_par%cnp / 100.0_DP
          else if (li==2) then
             if (kcn_type == 1) then
                kcni = -hii * dftb_par%cnd1 / 100.0_DP
             else if (kcn_type == 2) then
                kcni = -hii * dftb_par%cnd2 / 100.0_DP
             else
                call utils_abort(myself//': invalid kcn_type')
             end if
          else
             call utils_abort(myself//': dftb works for s, p, d elements only')
          end if
       end if
    else if(pub_dftb_method==DFTB_GFN0 .or. pub_dftb_method==DFTB_GFN2) then
       ! ab: values read from the parameter file are divided by 10 as per
       ! values given in Table S52 of ref. [3] and similar rescaling done in
       ! line 563-565 of xtb/src/read_gfn_param.f90
       if(li==0) then
          kcni = dftb_par%z(zi)%kcns / 10.0_DP / HARTREE_IN_EVS
       else if (li==1) then
          kcni = dftb_par%z(zi)%kcnp / 10.0_DP / HARTREE_IN_EVS
       else if (li==2) then
          kcni = dftb_par%z(zi)%kcnd / 10.0_DP / HARTREE_IN_EVS
       else
          call utils_abort(myself//': DFTB works for s, p, d elements only')
       end if
    else
       call utils_abort(myself//': Unsupported DFTB method. &
            &Only GFN0, GFN1 and GFN2 are supported.')
    end if

  end function dftb_kcn

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function dftb_kq(zi, li, dftb_par) result(kqi)
    !==========================================================================!
    ! This function returns kq (coefficient of charge).                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! zi(in)        : atomic number for atom I.                                !
    ! li(in)        : angular momentum for function i.                         !
    ! dftb_par (in) : dftb parameter set.                                      !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in May 2022.                                 !
    !==========================================================================!

    use constants, only: DFTB_GFN0, DFTB_GFN1, DFTB_GFN2, HARTREE_IN_EVS
    use dftb_parameters, only: DFTB_PARAMS
    use rundat, only: pub_dftb_method
    use utils, only: utils_abort

    implicit none

    ! ab: Arguments
    integer, intent(in)           :: zi
    integer, intent(in)           :: li
    type(DFTB_PARAMS), intent(in) :: dftb_par

    ! ab: Local variables
    real(kind=DP) :: kqi
    character(len=*), parameter :: myself='dftb_kq'

    ! -----------------------------------------------------------------------

    if (pub_dftb_method==DFTB_GFN0) then
       if (li==0) then
          kqi = dftb_par%z(zi)%kqs / HARTREE_IN_EVS
       else if (li==1) then
          kqi = dftb_par%z(zi)%kqp / HARTREE_IN_EVS
       else if (li==2) then
          kqi = dftb_par%z(zi)%kqd / HARTREE_IN_EVS
       else
          call utils_abort(myself//': DFTB works for s, p, d elements only')
       end if
    else if(pub_dftb_method==DFTB_GFN1 .or. pub_dftb_method==DFTB_GFN2) then
       kqi = 0.0_DP
    else
       call utils_abort(myself//': Unsupported DFTB method. &
            &Only GFN0, GFN1 and GFN2 are supported.')
    end if

  end function dftb_kq

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function dftb_calculate_kij(zi, zj, li, lj, dftb_par, vi, vj)
    !==========================================================================!
    ! This subroutine calculates the scaling factor k(zi,zj,li,lj) for         !
    ! GFN0: [1]:(19), GFN1: [4]:(10), GFN2: [3]:(16)                           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! zi, zj (in)                 : atomic no. for corresponding atoms I and J.!
    ! li, lj (in)                 : angular momentum for functions i and j.    !
    ! dftb_par (in)               : dftb parameter set.                        !
    ! vi, vj (in)                 : logical indicating a valence orbital.      !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in July 2021                                 !
    !==========================================================================!

    use constants, only: DFTB_GFN0, DFTB_GFN1, DFTB_GFN2
    use dftb_parameters, only: DFTB_PARAMS
    use rundat, only: pub_dftb_method
    use utils, only: utils_abort

    implicit none

    ! ab: Arguments
    integer, intent(in)           :: zi, zj
    integer, intent(in)           :: li, lj
    type(DFTB_PARAMS), intent(in) :: dftb_par
    logical, intent(in)           :: vi, vj

    ! ab: Local variables
    integer :: ali, blj
    real(kind=DP) :: k_l(0:2), k_li_lj(0:2,0:2), k_zi_zj, kv
    character(len=*), parameter :: myself='dftb_calculate_kij'

    ! -------------------------------------------------------------------------

    k_l(0) = dftb_par%ks
    k_l(1) = dftb_par%kp
    k_l(2) = dftb_par%kd

    if (pub_dftb_method==DFTB_GFN0) then
       if (zi==1 .and. vi) then
          if (zj==1 .and. vj) then
             dftb_calculate_kij = 0.0_DP
             return
          else
             kv = dftb_par%kdiff
          end if
       else
          if (zj==1 .and. vj) then
             kv = dftb_par%kdiff
          else
             kv = 1.0_DP
          end if
       end if
    else if (pub_dftb_method==DFTB_GFN1) then
       if (zi==1 .and. vi) then
          if (zj==1 .and. vj) then
             dftb_calculate_kij = dftb_par%kdiff
             return
          else
             dftb_calculate_kij = 0.5_DP * ( dftb_par%kdiff + k_l(lj) )
             return
          end if
       else
          if (zj==1 .and. vj) then
             dftb_calculate_kij = 0.5_DP * ( dftb_par%kdiff + k_l(li) )
             return
          else
             kv = 1.0_DP
          end if
       end if
    else if (pub_dftb_method==DFTB_GFN2) then
       kv=1.0_DP
    else
       call utils_abort(myself//': Unsupported DFTB method. &
            &Only GFN0, GFN1 and GFN2 are supported.')
    end if

    ! ab: k_li_lj is calculated as average of k_li and k_lj
    do blj = 0,2
       do ali = 0,2
          k_li_lj(ali,blj) = 0.5_DP*(k_l(ali) + k_l(blj))
       end do
    end do

    ! ab: except the following cases:
    if(pub_dftb_method==DFTB_GFN0) then
    else if (pub_dftb_method==DFTB_GFN1) then
       k_li_lj(0,1) = dftb_par%ksp
       k_li_lj(1,0) = dftb_par%ksp
    else if (pub_dftb_method==DFTB_GFN2) then
       k_li_lj(0,2) = dftb_par%ksd
       k_li_lj(2,0) = dftb_par%ksd
       k_li_lj(2,1) = dftb_par%kpd
       k_li_lj(1,2) = dftb_par%kpd
    else
       call utils_abort(myself//': Unsupported DFTB method. &
            &Only GFN0, GFN1 and GFN2 are supported.')
    end if

    k_zi_zj = 1.0_DP ! ab: generally, except following cases:

    if(pub_dftb_method==DFTB_GFN0) then

       if(  (zi == 29 .or. zi == 47 .or. zi == 79) .and. &
            (zj == 29 .or. zj == 47 .or. zj == 79)  ) then

          k_zi_zj = 0.90_DP ! for Cu(29), Ag(47) and Au(79)

       else if ( (((zi>20).and.(zi<31)) &       ! 3d elements
            .or.  ((zi>38).and.(zi<49)) &       ! 4d elements
            .or.  ((zi>71).and.(zi<81))) &      ! 5d elements
            .and.(((zj>20).and.(zj<31)) &       ! 3d elements
            .or.  ((zj>38).and.(zj<49)) &       ! 4d elements
            .or.  ((zj>71).and.(zj<81))) ) then ! 5d elements

          k_zi_zj = 1.10_DP ! for transition elements pair

       end if

    else if(pub_dftb_method==DFTB_GFN1) then

       ! ab: specific pair scaling factor for Hydrides and Borides
       if(  ((zi == 1) .and. (zj == 1))  .or. & ! HH
            ((zi == 1) .and. (zj == 5))  .or. & ! HB
            ((zi == 5) .and. (zj == 1))  .or. & ! BH
            ((zi == 1) .and. (zj == 7))  .or. & ! HN
            ((zi == 7) .and. (zj == 1))  .or. & ! NH
            ((zi == 1) .and. (zj ==28))  .or. & ! HNi
            ((zi ==28) .and. (zj == 1))  .or. & ! NiH
            ((zi == 1) .and. (zj ==75))  .or. & ! HRe
            ((zi ==75) .and. (zj == 1))  .or. & ! ReH
            ((zi == 1) .and. (zj ==78))  .or. & ! HPt
            ((zi ==78) .and. (zj == 1))  .or. & ! PtH
            ((zi ==15) .and. (zj == 5))  .or. & ! PB
            ((zi == 5) .and. (zj ==15))  .or. & ! BP
            ((zi ==14) .and. (zj == 5))  .or. & ! SiB
            ((zi == 5) .and. (zj ==14))  ) then ! BSi

          k_zi_zj = dftb_par%zz(zi,zj)%x      ! for above pairs

       ! ab: specific pair scaling factor for transition elements
       else if ( (((zi>20).and.(zi<31)) &       ! 3d elements
            .or.  ((zi>38).and.(zi<49)) &       ! 4d elements
            .or.  ((zi>71).and.(zi<81))) &      ! 5d elements
            .and.(((zj>20).and.(zj<31)) &       ! 3d elements
            .or.  ((zj>38).and.(zj<49)) &       ! 4d elements
            .or.  ((zj>71).and.(zj<81))) ) then ! 5d elements


          ! ab: case 1: both are 3d elements
          if ( ((zi>20).and.(zi<31)) .and. &    ! 3d elements
               ((zj>20).and.(zj<31)) ) then     ! 3d elements

             k_zi_zj = 1.10_DP

          ! ab: case 2: both are 4d/5d elements
          else if ( (((zi>38).and.(zi<49))  &      ! 4d elements
               .or.  ((zi>71).and.(zi<81))) &      ! 5d elements
               .and.(((zj>38).and.(zj<49))  &      ! 4d elements
               .or.  ((zj>71).and.(zj<81))) ) then ! 5d elements

             k_zi_zj = 1.20_DP

          ! ab: case 3: mixed pair 3d - 4d/5d
          else

             k_zi_zj = 1.15_DP

          end if

       end if

    end if

    dftb_calculate_kij = kv * k_zi_zj * k_li_lj(li,lj)

  end function dftb_calculate_kij

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function dftb_calculate_pij(zi, zj, li, lj, R_ij, dftb_par)
    !==========================================================================!
    ! This subroutine calculates the function P(Rij,li,lj) from [5]:(24)       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! zi, zj (in)                 : atomic no. for corresponding atoms I and J.!
    ! li, lj (in)                 : angular momentum for functions i and j.    !
    ! dftb_par (in)               : dftb parameter set.                        !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in July 2021                                 !
    !==========================================================================!

    use constants, only: SAFE_DIV_EPS
    use dftb_parameters, only: DFTB_PARAMS
    use utils, only: utils_abort

    implicit none

    ! ab: Arguments
    integer, intent(in)           :: zi, zj
    integer, intent(in)           :: li, lj
    real(kind=DP), intent(in)     :: R_ij      ! ab: distance between atoms I, J
    type(DFTB_PARAMS), intent(in) :: dftb_par

    ! ab: Local variables
    real(kind=DP) :: R_atomicrad, k_poly_i, k_poly_j
    character(len=*), parameter :: myself='dftb_calculate_pij'

    ! -------------------------------------------------------------------------

    R_atomicrad = dftb_par%z(zi)%atomicrad + dftb_par%z(zj)%atomicrad

    k_poly_i = dftb_kpoly(zi, li, dftb_par)
    k_poly_j = dftb_kpoly(zj, lj, dftb_par)

    if(abs(R_atomicrad) >= SAFE_DIV_EPS) then
       dftb_calculate_pij =  &
            ( 1.0_DP + k_poly_i * sqrt(R_ij/R_atomicrad) ) * &
            ( 1.0_DP + k_poly_j * sqrt(R_ij/R_atomicrad) )
    else
       call utils_abort(myself//': Covalent radius cannot be zero.')
    end if

  end function dftb_calculate_pij

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function dftb_calculate_dpij(zi, zj, li, lj, dR, dftb_par) result (dpij)
    !==========================================================================!
    ! This subroutine calculates the derivative of function P(Rij,li,lj).      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! zi, zj (in)                 : atomic no. for corresponding atoms I and J.!
    ! li, lj (in)                 : angular momentum for functions i and j.    !
    ! dftb_par (in)               : dftb parameter set.                        !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in May 2022.                                 !
    !==========================================================================!

    use constants, only: SAFE_DIV_EPS
    use dftb_parameters, only: DFTB_PARAMS
    use geometry, only: POINT, magnitude
    use utils, only: utils_abort

    implicit none

    ! ab: Arguments
    integer, intent(in)           :: zi, zj
    integer, intent(in)           :: li, lj
    type(POINT), intent(in)       :: dR
    type(DFTB_PARAMS), intent(in) :: dftb_par

    ! ab: Local variables
    real(kind=DP) :: R_ij      ! ab: distance between atoms I, J
    real(kind=DP) :: R_atomicrad, ki, kj, dpij(1:3), rf
    character(len=*), parameter :: myself='dftb_calculate_dpij'

    ! -------------------------------------------------------------------------

    R_ij = magnitude(dR)
    R_atomicrad = dftb_par%z(zi)%atomicrad + dftb_par%z(zj)%atomicrad

    ki = dftb_kpoly(zi, li, dftb_par)
    kj = dftb_kpoly(zj, lj, dftb_par)

    if(abs(R_atomicrad) >= SAFE_DIV_EPS) then
       rf = sqrt(R_ij/R_atomicrad)
    else
       call utils_abort(myself//': Covalent radius cannot be zero.')
    end if

    dpij =  0.5_DP*(ki+kj+2.0_DP*ki*kj*rf)*rf/R_ij**2*(/ dR%X,dR%Y,dR%Z /)

  end function dftb_calculate_dpij

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function dftb_kpoly(zi, li, dftb_par) result(kpoly)
    !==========================================================================!
    ! This function returns kq (coefficient of charge).                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! zi(in)        : atomic number for atom I.                                !
    ! li(in)        : angular momentum for function i.                         !
    ! dftb_par (in) : dftb parameter set.                                      !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in May 2022.                                 !
    !==========================================================================!

    use dftb_parameters, only: DFTB_PARAMS
    use utils, only: utils_abort

    implicit none

    ! ab: Arguments
    integer, intent(in)           :: zi
    integer, intent(in)           :: li
    type(DFTB_PARAMS), intent(in) :: dftb_par

    ! ab: Local variables
    real(kind=DP) :: kpoly
    character(len=*), parameter :: myself='dftb_kpoly'

    ! -----------------------------------------------------------------------

    ! ab: k_poly is divided by 100 as per similar rescaling done in xTB
    ! cf. line 63 - 64 of xtb/src/grad_core.f90 &
    ! cf. line 719 - 720 of xtb/src/scc_core.f90
    ! also the values of k_poly given in Table S52 of ref.[3] are smaller
    ! than that in the input file by factor of 100.
    if(li==0) then
       kpoly = dftb_par%z(zi)%polys / 100.0_DP
    else if (li==1) then
       kpoly = dftb_par%z(zi)%polyp / 100.0_DP
    else if (li==2) then
       kpoly = dftb_par%z(zi)%polyd / 100.0_DP
    else
       call utils_abort(myself//': DFTB works for s, p, d elements only')
    end if

  end function dftb_kpoly

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function dftb_calculate_xij(zi, zj, li, lj, dftb_par, vi, vj)
    !==========================================================================!
    ! This subroutine calculates the function X(zi,zj,li,lj) from Pauling data !
    ! of electronegativity of elements of periodic table following:            !
    ! GFN0: [1]:(20), GFN1: [4]:(10) & [5]:(28), GFN2: [3]:(16)                !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! zi, zj (in)                 : atomic no. for corresponding atoms I and J.!
    ! li, lj (in)                 : angular momentum for functions i and j.    !
    ! dftb_par (in)               : dftb parameter set.                        !
    ! vi, vj (in)                 : logical indicating a valence orbital.      !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in July 2021                                 !
    !==========================================================================!

    use constants, only: DFTB_GFN0, DFTB_GFN1, DFTB_GFN2
    use dftb_parameters, only: DFTB_PARAMS
    use rundat, only: pub_dftb_method
    use utils, only: utils_abort

    implicit none

    ! ab: Arguments
    integer, intent(in)           :: zi, zj
    integer, intent(in)           :: li, lj
    type(DFTB_PARAMS), intent(in) :: dftb_par
    logical, intent(in)           :: vi, vj

    ! ab: Local variables
    integer :: ali, blj
    real(kind=DP) :: eni, enj
    real(kind=DP) :: delta_en, k_en, b_en, k_en_l(0:2), k_en_li_lj(0:2,0:2), xij
    character(len=*), parameter :: myself='dftb_calculate_xij'

    ! -------------------------------------------------------------------------

    if (pub_dftb_method==DFTB_GFN1) then
       if (zi==1 .and. vi) then
          dftb_calculate_xij = 1.0_DP
          return
       else
          if (zj==1 .and. vj) then
             dftb_calculate_xij = 1.0_DP
             return
          end if
       end if
    end if

    if (dftb_par%z(zi)%en==garbage_real) then
       eni = dftb_par%z(zi)%paulingen
    else
       eni = dftb_par%z(zi)%en
    end if
    if (dftb_par%z(zj)%en==garbage_real) then
       enj = dftb_par%z(zj)%paulingen
    else
       enj = dftb_par%z(zj)%en
    end if

    delta_en = abs(eni - enj)

    if(pub_dftb_method==DFTB_GFN0) then
    ! ab: parameters read from the input file are divided by 100 as per values
    ! given in Table 2 of ref. [1] and similar rescaling done inside xTB
    ! cf. line 800 of xtb/src/xtb/gfn0.f90
       k_en_l(0) = 0.01_DP * dftb_par%ens
       k_en_l(1) = 0.01_DP * dftb_par%enp
       k_en_l(2) = 0.01_DP * dftb_par%end
       b_en      = dftb_par%enscale4

       do blj = 0,2
          do ali = 0,2
             k_en_li_lj(ali,blj) = 0.5_DP*(k_en_l(ali) + k_en_l(blj))
          end do
       end do

       xij = 1.0_DP + k_en_li_lj(li,lj) * delta_en**2.0_DP * &
            (1.0_DP + b_en              * delta_en**2.0_DP)

    else if(pub_dftb_method==DFTB_GFN1 .or. pub_dftb_method==DFTB_GFN2) then
    ! ab: parameters read from the input file are divided by 100 as per values
    ! given in Table 2 of ref. [3 and 4] and similar rescaling done inside xTB
    ! cf. line 653 of xtb/src/xtb/gfn1.f90 & line 837 of xtb/src/xtb/gfn2.f90
       k_en = 0.01_DP * dftb_par%enscale

       xij = 1.0_DP + k_en * delta_en**2.0_DP

    else
       call utils_abort(myself//': Unsupported DFTB method. &
            &Only GFN0, GFN1 and GFN2 are supported.')

    end if

    dftb_calculate_xij = xij

  end function dftb_calculate_xij

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function dftb_calculate_yij(zi, zj, ni, nj, dftb_par)
    !==========================================================================!
    ! This subroutine calculates the function y of slater exponents for        !
    ! GFN0: [1]:(18), GFN1: [4]:(10) & [5]:(above 28), GFN2: [5]:(37)          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! zi, zj (in)                 : atomic no. for corresponding atoms I and J.!
    ! ni, nj (in)                 : ngwf_type for functions i and j.           !
    ! dftb_par (in)               : dftb parameter set.                        !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in July 2021                                 !
    !==========================================================================!

    use constants, only: DFTB_GFN0, DFTB_GFN1, DFTB_GFN2, SAFE_DIV_EPS
    use dftb_parameters, only: DFTB_PARAMS
    use rundat, only: pub_dftb_method
    use utils, only: utils_abort

    implicit none

    ! ab: Arguments
    integer, intent(in)           :: zi, zj
    integer, intent(in)           :: ni, nj
    type(DFTB_PARAMS), intent(in) :: dftb_par

    ! ab: Local variables
    real(kind=DP) :: expi, expj, yij
    character(len=*), parameter :: myself='dftb_calculate_yij'

    ! -------------------------------------------------------------------------

    if (zi==zj .and. ni==nj) then
       dftb_calculate_yij = 1.0_DP
       return
    end if

    expi = dftb_par%z(zi)%exp(ni)
    expj = dftb_par%z(zj)%exp(nj)

    if(pub_dftb_method==DFTB_GFN0) then

       if( abs(expi+expj) >= SAFE_DIV_EPS ) then
          yij = 2.0_DP * sqrt(expi*expj) / (expi + expj)
       else
          call utils_abort(myself//': Slater exponents sum up to be zero.')
       end if

    else if(pub_dftb_method==DFTB_GFN1) then

       yij = 1.0_DP

    else if(pub_dftb_method==DFTB_GFN2) then

       if( abs(expi+expj) >= SAFE_DIV_EPS ) then
          yij = sqrt( 2.0_DP * sqrt(expi*expj) / (expi + expj) )
       else
          call utils_abort(myself//': Slater exponents sum up to be zero.')
       end if

    else
       call utils_abort(myself//': Unsupported DFTB method. &
            &Only GFN0, GFN1 and GFN2 are supported.')

    end if

    dftb_calculate_yij = yij

  end function dftb_calculate_yij

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_analytical_gaussian_overlap(overlap, mdl, ngwf_basis)
    !==========================================================================!
    ! Calculates the overlap matrix through analytical overlap integrals       !
    ! between contracted Gaussians.                                            !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   overlap (inout): The overlap SPAM3 to fill in. Assumed to have correct !
    !                    sparsity when passed in.                              !
    !   mdl (in): MODEL_TYPE required for dftb_par and elements.               !
    !   ngwf_basis (in): Required for counting NGWFs.                          !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2021                                !
    ! Re-written by Arihant Bhandari in Nov-Dec 2021.                          !
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use function_basis, only: FUNC_BASIS
    use gaussian_ops, only: gaussian_overlap_integral
    use geometry, only: POINT, OPERATOR(+), OPERATOR(-), magnitude
    use model_type, only: MODEL
    use parallel_strategy, only: PARAL_INFO
    use simulation_cell, only: minimum_image_displacement
    use sparse, only: SPAM3, sparse_put_element, sparse_get_par
    use rundat, only: pub_dftb_bc_is_periodic, pub_dftb_overlap_cutoff
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    type(SPAM3), intent(inout)   :: overlap
    type(MODEL), intent(in)      :: mdl
    type(FUNC_BASIS), intent(in) :: ngwf_basis

    ! jd: Local variables
    type(PARAL_INFO), pointer :: par
    type(POINT) :: r_i, r_j, r_ij    ! ab: Positions of atoms I, J & displacement
    integer :: local_i, local_j
    integer :: global_i, global_j
    integer :: orig_i, orig_j
    integer :: first_ngwf_of_i, first_ngwf_of_j
    integer :: ngwf_i, ngwf_j
    integer :: z_i, z_j
    integer :: s_i, s_j
    integer :: nshells_i, nshells_j
    integer :: nprim_s_i, nprim_s_j, prim_i, prim_j
    integer :: angmom_s_i, angmom_s_j
    integer :: ix, iy, iz, jx, jy, jz
    real(kind=DP) :: a_i, a_j
    real(kind=DP) :: c_i, c_j
    real(kind=DP) :: oij, int_ij, d_ij
    character(len=*), parameter :: myself = 'dftb_analytical_gaussian_overlap'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering '//myself

    call sparse_get_par(par, overlap)

    ! --------------------------------------------------------------------------
    ! jd: Loop over all atoms on this rank                                  III
    ! --------------------------------------------------------------------------
    loop_I:                                                                    &
    do local_i = 1, par%num_atoms_on_proc(pub_my_proc_id)
       global_i = par%first_atom_on_proc(pub_my_proc_id) + local_i-1
       first_ngwf_of_i = ngwf_basis%first_on_atom(global_i)
       orig_i = par%orig_atom(global_i)
       z_i = mdl%elements(orig_i)%atomic_number
       r_i = mdl%elements(orig_i)%centre
       nshells_i = mdl%dftb_par%z(z_i)%nshells

       ngwf_i = first_ngwf_of_i - 1
       ! --------------------------------------------------------------------
       ! jd: Loop over all shells s_i on atom I                          s_i
       ! --------------------------------------------------------------------
       loop_si:                                                             &
       do s_i = 1, nshells_i
          nprim_s_i = mdl%dftb_par%z(z_i)%sto(s_i)%nprim
          angmom_s_i = mdl%dftb_par%z(z_i)%ang(s_i)

          !-----------------------------------------------------------------
          ! ab: Loop over all orbitals in s_i with angular moment ix, iy, iz
          !-----------------------------------------------------------------
          loop_ix:                                                         &
          do ix = 0, angmom_s_i

             loop_iy:                                                      &
             do iy = 0, angmom_s_i

                loop_iz:                                                   &
                do iz = 0, angmom_s_i

                   if (ix+iy+iz .ne. angmom_s_i) cycle

                   ! ab: ngwf counter
                   ngwf_i = ngwf_i + 1

                   ! -----------------------------------------------------------------------
                   ! jd: Loop over all atoms J that are s-neighbours with I             JJJ
                   ! -----------------------------------------------------------------------
                   loop_J:                                                                &
                   do local_j = mdl%dftb_var%nl%first_idx(global_i),mdl%dftb_var%nl%last_idx(global_i)
                      global_j = mdl%dftb_var%nl%neighbours(local_j)
                      first_ngwf_of_j = ngwf_basis%first_on_atom(global_j)
                      orig_j = par%orig_atom(global_j)
                      z_j = mdl%elements(orig_j)%atomic_number
                      r_j = mdl%elements(orig_j)%centre
                      nshells_j = mdl%dftb_par%z(z_j)%nshells

                      ngwf_j = first_ngwf_of_j - 1

                      ! ab: Take care of boundary conditions
                      if (.not. any(pub_dftb_bc_is_periodic)) then
                         ! ab: do nothing in OBCs
                         r_ij = r_i - r_j
                         d_ij = magnitude(r_ij)
                      else if (all(pub_dftb_bc_is_periodic)) then
                         ! ab: find the minimum image of atom I
                         call minimum_image_displacement(d_ij, r_ij, r_i, r_j, mdl%cell)
                         r_i = r_j + r_ij
                      else
                         call utils_abort(myself//': Mixed BCs are not currently supported.')
                      end if
                      ! -------------------------------------------------------------
                      ! jd: Loop over all shells s_j on atom J                    s_j
                      ! -------------------------------------------------------------
                      loop_sj:                                                       &
                      do s_j = 1, nshells_j
                         nprim_s_j = mdl%dftb_par%z(z_j)%sto(s_j)%nprim
                         angmom_s_j = mdl%dftb_par%z(z_j)%ang(s_j)

                         !---------------------------------------------------
                         ! ab: Loop over all orbitals in s_j with  jx, jy, jz
                         !---------------------------------------------------
                         loop_jx:                                           &
                         do jx = 0, angmom_s_j

                            loop_jy:                                        &
                            do jy = 0, angmom_s_j

                               loop_jz:                                     &
                               do jz = 0, angmom_s_j

                                  if (jx+jy+jz .ne. angmom_s_j) cycle

                                  ! ab: ngwf counter
                                  ngwf_j = ngwf_j + 1

                                  ! ab: initialize overlap integral to be zero
                                  oij = 0.0_DP

                                  ! ab: cutoff for overlap matrix
                                  if (d_ij .le. pub_dftb_overlap_cutoff) then
                                     !---------------------------------------------
                                     ! ab: Loop over all primitives          prim_i
                                     !---------------------------------------------
                                     loop_prim_i:                                 &
                                     do prim_i = 1, nprim_s_i

                                        ! ab: exponent
                                        a_i= mdl%dftb_par%z(z_i)%sto(s_i)%&
                                             alpha(prim_i)
                                        ! ab: contraction coefficient
                                        c_i= mdl%dftb_par%z(z_i)%sto(s_i)%&
                                             coeff(prim_i)

                                        !------------------------------------------
                                        ! ab: Loop over all primitives       prim_j
                                        !------------------------------------------
                                        loop_prim_j:                              &
                                        do prim_j = 1, nprim_s_j

                                           ! ab: exponent
                                           a_j= mdl%dftb_par%z(z_j)%sto(s_j)%&
                                                alpha(prim_j)
                                           ! ab: contraction coefficient
                                           c_j= mdl%dftb_par%z(z_j)%sto(s_j)%&
                                                coeff(prim_j)

                                           ! ab: get overlap integral
                                           int_ij = gaussian_overlap_integral( &
                                                a_i,a_j,r_i,r_j,ix,jx,iy,jy,iz,jz)

                                           ! ab: sum over all contractions
                                           oij = oij + c_i*c_j*int_ij

                                        end do loop_prim_j

                                     end do loop_prim_i

                                  end if

                                  ! ab: put element in the overlap matrix
                                  call sparse_put_element(oij, overlap, &
                                       ngwf_j, ngwf_i)

                               end do loop_jz

                            end do loop_jy

                         end do loop_jx

                      end do loop_sj

                   end do loop_J

                end do loop_iz

             end do loop_iy

          end do loop_ix

       end do loop_si

    end do loop_I

    call timer_clock(myself,2)

  end subroutine dftb_analytical_gaussian_overlap

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_ies(mdl, energy, force)
    !==========================================================================!
    ! Calculates the isotropic electrostatic (ies) energy and force.           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   mdl(inout):  mdl%elements%dftb_charge is modified.                     !
    !   energy(out): total energy.                                             !
    !   force(out):  force on each atom.                                       !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Nov 2021.                                 !
    ! Forces added in May 2022.                                                !
    !==========================================================================!

    use comms, only: pub_my_proc_id, comms_reduce
    use constants, only: DFTB_GFN0
    use dense, only: DEM, dense_create, dense_product, dense_destroy, &
         dense_get_col, dense_put_col, dense_get_element, dense_put_element, &
         dense_invert, dense_vector_product, dense_axpy, dense_dot_product
    use model_type, only: MODEL
    use rundat, only: pub_dftb_method, pub_dftb_bc_is_periodic, pub_print_qc, &
         pub_dftb_calculate_force, pub_charge
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_assert

    implicit none

    ! ab: Arguments
    type(MODEL), intent(inout)    :: mdl
    real(kind=DP), intent(out)    :: energy
    real(kind=DP), intent(out)    :: force(1:3,size(mdl%elements))

    ! ab: Local variables
    character(len=*), parameter :: myself = 'dftb_ies'
    character(len=3) :: err_str
    integer      :: ierr, ic, zi
    integer      :: nat       ! ab: total number of atoms
    integer      :: iat, local_iat, global_iat, kat  ! ab: atom id
    real(kind=DP):: qk, ak(mdl%nat+1), bk(mdl%nat+1), daqk, ci
    type(DEM) :: aa, bb, qq, daij(3), es
#ifdef SCALAPACK
    real(kind=DP):: bi
    integer :: nprow, npcol, myrow, mycol, li
    ! ab: BLACS subroutines
    external :: blacs_gridinfo
#endif

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    select case (pub_dftb_method)
    case(DFTB_GFN0)
    case default
       call utils_abort(myself//': Unsupported DFTB method. &
            &Only GFN0 is supported.')
    end select

    nat = mdl%nat

    ! ab: allocate matrices
    call dense_create(aa, nat+1, nat+1, iscmplx=.false.)
    call dense_create(bb, nat+1, 1, iscmplx=.false.)
    call dense_create(qq, nat+1, 1, iscmplx=.false.)
    call dense_create(es, 1, 1, iscmplx=.false.)

    if (pub_dftb_calculate_force) then
       do ic = 1, 3
          call dense_create(daij(ic), nat+1, nat+1, iscmplx=.false.)
          call dense_create(mdl%dftb_var%dqdr(ic), nat+1, nat+1, iscmplx=.false.)
       end do
    end if

    ! ab: create matrices and their derivatives
    if (.not. any(pub_dftb_bc_is_periodic)) then
       call dftb_gfn0_ies_obc(mdl, aa, daij)
    else if(all(pub_dftb_bc_is_periodic)) then
       call dftb_gfn0_ies_pbc(mdl, aa, daij)
    else
       call utils_abort(myself//': Mixed BCs are not currently supported.')
    end if

    ! ab: create B matrix
#ifdef SCALAPACK
    ! ab: block-cyclically distributed matrix
    call blacs_gridinfo(bb%blacs_desc(2), nprow, npcol, myrow, mycol)
    if (bb%blacs_ncol==1) then
       do li = 0, bb%blacs_ld - 1
          iat=((li/bb%blacs_nb)*nprow+myrow)*bb%blacs_nb+mod(li,bb%blacs_nb)+1
          if (iat==nat+1) then
             call dense_put_element(pub_charge, bb, iat, 1)
             cycle
          end if
          zi = mdl%elements(iat)%atomic_number
          ci = mdl%elements(iat)%dftb_coord
          bi = mdl%dftb_par%z(zi)%kappa*sqrt(ci) - mdl%dftb_par%z(zi)%xi
          call dense_put_element(bi, bb, iat, 1)
       end do
    end if
#else
    ! ab: undistributed matrix (allocated on all procs)
    do local_iat = 1, mdl%regions(1)%par%num_atoms_on_proc(pub_my_proc_id)
       global_iat = mdl%regions(1)%par%first_atom_on_proc(pub_my_proc_id) + local_iat-1
       iat = mdl%regions(1)%par%orig_atom(global_iat)
       zi = mdl%elements(iat)%atomic_number
       ci = mdl%elements(iat)%dftb_coord
       bb%dmtx(iat,1) = mdl%dftb_par%z(zi)%kappa*sqrt(ci) - mdl%dftb_par%z(zi)%xi
    end do
    call comms_reduce('SUM', bb%dmtx)
    bb%dmtx(nat+1,1) = pub_charge
#endif

    ! ab: calculate inverse of A.
    call dense_invert(aa, ierr)
    write(err_str, '(i3)') ierr
    call utils_assert(ierr == 0, myself//'matrix inversion failed with exit code ' &
         // adjustl(err_str))

    ! ab: qq = aa^-1 * bb
    call dense_vector_product(qq, aa, bb)

    ! ab: calculate isotropic electrostatic energy
    ! ab: e_ies = q^T ( 0.5*A*q - B ), eq. [1]:(7)
    ! ab:       = -1/2 q^T * B
#ifdef SCALAPACK
    call dense_product(es, qq, bb, opA='T', opB='N')
    call dense_get_element(energy, es, 1, 1)
    energy = -0.5_DP * energy
#else
    energy = -0.5_DP * dense_dot_product(qq, bb)
#endif

    ! ab: store the charges in mdl
    do iat = 1, nat
       call dense_get_element(mdl%elements(iat)%dftb_charge, qq, iat, 1)
    end do

    ! ab: forces
    if (pub_dftb_calculate_force) then
       do ic = 1, 3
          call dense_vector_product(bb, daij(ic), qq)
          do kat = 1, nat
             call dense_get_element(qk, qq, kat, 1)
             call dense_get_col(ak, daij(ic), kat)
             call dense_get_col(bk, mdl%dftb_var%dxdr(ic), kat)
             bk = bk + qk * ak
             call dense_put_col(bk, mdl%dftb_var%dxdr(ic), kat)
             call dense_get_col(ak, aa, kat)
             call dense_get_element(daqk, bb, kat, 1)
             bk = -ak * daqk
             call dense_put_col(bk, mdl%dftb_var%dqdr(ic), kat)
          end do
          call dense_product(bb, mdl%dftb_var%dxdr(ic), qq, opA='T', opB='N')
          call dense_get_col(bk, bb, 1)
          force(ic,:) = bk(1:nat)
          call dense_product(daij(ic), aa, mdl%dftb_var%dxdr(ic))
          call dense_axpy(mdl%dftb_var%dqdr(ic), daij(ic), 1.0_DP)
          call dense_destroy(daij(ic))
       end do
    end if

    ! ab: deallocate matrices
    call dense_destroy(aa)
    call dense_destroy(bb)
    call dense_destroy(qq)
    call dense_destroy(es)

    if (pub_print_qc) then
       call dftb_qc_print(myself,energy,size(mdl%elements),force)
    end if

    call timer_clock(myself,2)

  end subroutine dftb_ies

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_gfn0_ies_obc(mdl, aa, daij)
    !==========================================================================!
    ! This subroutine calculates the matrix 'A' of coulombic interactions of   !
    ! smeared Gaussians and its gradient.                                      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   mdl(in):  model type.                                                  !
    !   aa(inout):  The 'A' matrix.                                            !
    !   daij(3)(inout): Gradient of the A matrix.                              !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Nov 2021                                  !
    !==========================================================================!

    use comms, only: pub_my_proc_id, comms_reduce
    use constants, only: PI, SAFE_DIV_EPS
    use dense, only: DEM, dense_put_element
    use ewald, only: ewald_real_derivative
    use geometry, only: POINT, magnitude, OPERATOR(-)
    use model_type, only: MODEL
    use rundat, only: pub_dftb_calculate_force
    use utils, only: utils_abort, utils_erf

    implicit none

    ! ab: Arguments
    type(MODEL), intent(in)            :: mdl
    type(DEM), intent(inout)           :: aa
    type(DEM), intent(inout)           :: daij(3)

    ! ab: Local variables
    character(len=*), parameter  :: myself='dftb_gfn0_ies_obc'
    integer      :: ic, local_jat, global_jat
    integer      :: nat       ! ab: total number of atoms
    integer      :: iat, jat  ! ab: atom id
    integer      :: zi, zj    ! ab: atomic number
    type(POINT)  :: ri, rj    ! ab: Positions of atoms I, J
    type(POINT)  :: dr        ! ab: displacement vector
    real(kind=DP):: rij       ! ab: distance between atoms I and J
    real(kind=DP):: gij       ! ab: gamma in eq. [1]:(9)
    real(kind=DP):: ad(3)
    real(kind=DP):: aij
#ifdef SCALAPACK
    integer :: nprow, npcol, myrow, mycol, lj, li
    ! ab: BLACS subroutines
    external :: blacs_gridinfo
#endif

    ! --------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering '//myself

    nat = mdl%nat

#ifdef SCALAPACK
    ! ab: block-cyclically distributed matrices
    call blacs_gridinfo(aa%blacs_desc(2), nprow, npcol, myrow, mycol)

    do lj = 0, aa%blacs_ncol - 1
       jat=((lj/aa%blacs_nb)*npcol+mycol)*aa%blacs_nb+mod(lj,aa%blacs_nb)+1

       do li = 0, aa%blacs_ld - 1
          iat=((li/aa%blacs_nb)*nprow+myrow)*aa%blacs_nb+mod(li,aa%blacs_nb)+1

          if (jat==nat+1) then
             if (iat/=nat+1) then
                call dense_put_element(1.0_DP, aa, iat, jat)
             end if
             cycle
          else
             if (iat==nat+1) then
                call dense_put_element(1.0_DP, aa, iat, jat)
                cycle
             end if
          end if

          zj = mdl%elements(jat)%atomic_number
          rj = mdl%elements(jat)%centre

          zi = mdl%elements(iat)%atomic_number
          ri = mdl%elements(iat)%centre

          gij = sqrt(mdl%dftb_par%z(zi)%alpg**2+mdl%dftb_par%z(zj)%alpg**2)
          if (gij > SAFE_DIV_EPS) then
             gij = 1.0_DP/gij
          else
             call utils_abort(myself//': Division by zero in calculating gij.')
          end if

          dr = ri - rj
          rij = magnitude(dr)

          ! ab: A matrix in eq. [1]:(9)
          if (iat == jat) then
             aij = mdl%dftb_par%z(zi)%gam + 2.0_DP*gij/sqrt(PI)
          else
             aij = utils_erf(gij*rij)/rij
          end if
          call dense_put_element(aij, aa, iat, jat)

          if (pub_dftb_calculate_force) then
             if (iat /= jat) then
                ad = ewald_real_derivative(dr,gij)
                do ic = 1, 3
                   call dense_put_element(ad(ic), daij(ic), iat, jat)
                end do
             end if
          end if
       end do
    end do

#else
    ! ab: undistributed matrix (allocated on all procs)
    aa%dmtx(nat+1,:) = 1.0_DP
    aa%dmtx(:,nat+1) = 1.0_DP
    aa%dmtx(nat+1, nat+1) = 0.0_DP

    do local_jat = 1, mdl%regions(1)%par%num_atoms_on_proc(pub_my_proc_id)
       global_jat = mdl%regions(1)%par%first_atom_on_proc(pub_my_proc_id) + local_jat-1
       jat = mdl%regions(1)%par%orig_atom(global_jat)
       zj = mdl%elements(jat)%atomic_number
       rj = mdl%elements(jat)%centre

       do iat = 1, nat
          zi = mdl%elements(iat)%atomic_number
          ri = mdl%elements(iat)%centre

          gij = sqrt(mdl%dftb_par%z(zi)%alpg**2+mdl%dftb_par%z(zj)%alpg**2)
          if (gij > SAFE_DIV_EPS) then
             gij = 1.0_DP/gij
          else
             call utils_abort(myself//': Division by zero in calculating gij.')
          end if

          dr = ri - rj
          rij = magnitude(dr)

          ! ab: A matrix in eq. [1]:(9)
          if (iat == jat) then
             aa%dmtx(iat,jat) = mdl%dftb_par%z(zi)%gam + 2.0_DP*gij/sqrt(PI)
          else
             ! ab: calculate distance between atoms
             aa%dmtx(iat,jat) = utils_erf(gij*rij)/rij
          end if

          if (pub_dftb_calculate_force) then
             if (iat == jat) then
                ad = 0.0_DP
             else
                ad = ewald_real_derivative(dr,gij)
             end if
             do ic = 1, 3
                daij(ic)%dmtx(iat,jat) = ad(ic)
             end do
          end if

       end do
    end do

    call comms_reduce('SUM', aa%dmtx)
    if (pub_dftb_calculate_force) then
       do ic = 1, 3
          call comms_reduce('SUM', daij(ic)%dmtx)
       end do
    end if

#endif

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving '//myself

  end subroutine dftb_gfn0_ies_obc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_gfn0_ies_pbc(mdl, aa, daij)
    !==========================================================================!
    ! This subroutine calculates the matrix 'A' of coulombic interactions of   !
    ! smeared Gaussians and its gradient in PBCs via Ewald sum.                !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   mdl(in):  model type.                                                  !
    !   aa(inout):  The 'A' matrix.                                            !
    !   daij(3)(inout): Gradient of the A matrix.                              !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Nov 2021                                  !
    !==========================================================================!

    use comms, only: pub_my_proc_id, comms_reduce
    use constants, only: PI, SAFE_DIV_EPS
    use dense, only: DEM, dense_put_element
    use ewald, only: ewald_real_sum, ewald_recp_sum, ewald_get_eta, &
         ewald_real_cutoff, ewald_recp_cutoff, ewald_real_derivative, &
         ewald_recp_sum_dr
    use geometry, only: POINT, magnitude, OPERATOR(-), OPERATOR(.CROSS.), &
         OPERATOR(.DOT.), OPERATOR(+)
    use model_type, only: MODEL
    use rundat, only: pub_dftb_ewald_parameter, &
         pub_dftb_ewald_replicate_xtb, pub_dftb_calculate_force
    use simulation_cell, only: simulation_cell_lattice_points_init, &
         simulation_cell_lattice_points_exit
    use utils, only: utils_abort

    implicit none

    ! ab: Arguments
    type(MODEL), intent(in)            :: mdl
    type(DEM), intent(inout)           :: aa
    type(DEM), intent(inout)           :: daij(3)

    ! ab: Local variables
    character(len=*), parameter  :: myself='dftb_gfn0_ies_obc'
    integer      :: ic, local_jat, global_jat
    integer      :: nat       ! ab: total number of atoms
    integer      :: iat, jat  ! ab: atom id
    integer      :: zi, zj    ! ab: atomic number
    type(POINT)  :: ri, rj    ! ab: Positions of atoms I, J
    type(POINT)  :: dr        ! ab: displacement vector
    real(kind=DP):: rij       ! ab: distance between atoms I and J
    real(kind=DP):: gij       ! ab: gamma in eq. [1]:(9)
    real(kind=DP):: ad(3)
    real(kind=DP):: aij
    real(kind=dp), parameter :: tol  = 1D-8   ! ab: tolerance
    type(POINT), allocatable :: real_points(:)
    type(POINT), allocatable :: recp_points(:)
    type(POINT)  :: real_lattice(3) ! ab: lattice vectors
    type(POINT)  :: recp_lattice(3) ! ab: reciprocal lattice vectors
    real(kind=DP):: real_cell_length(3) ! ab: magnitude of lattice vectors
    real(kind=DP):: recp_cell_length(3) ! ab: magnitude of rec. lattice vectors
    real(kind=DP):: real_cutoff ! ab: cutoff for real space sum
    real(kind=DP):: recp_cutoff ! ab: cutoff for reciprocal space sum
    real(kind=dp):: gterm, rterm
    real(kind=DP):: eta, volume
#ifdef SCALAPACK
    integer :: nprow, npcol, myrow, mycol, lj, li
    ! ab: BLACS subroutines
    external :: blacs_gridinfo
#endif

    ! --------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering '//myself

    nat = mdl%nat

    ! ab: lattice vectors:
    real_lattice(1) = mdl%cell%a1
    real_lattice(2) = mdl%cell%a2
    real_lattice(3) = mdl%cell%a3
    real_cell_length(1) = magnitude(real_lattice(1))
    real_cell_length(2) = magnitude(real_lattice(2))
    real_cell_length(3) = magnitude(real_lattice(3))
    volume = abs((real_lattice(1).CROSS.real_lattice(2)).DOT.real_lattice(3))

    ! ab: reciprocal lattice vectors:
    recp_lattice(1) = mdl%cell%b1
    recp_lattice(2) = mdl%cell%b2
    recp_lattice(3) = mdl%cell%b3
    recp_cell_length(1) = magnitude(recp_lattice(1))
    recp_cell_length(2) = magnitude(recp_lattice(2))
    recp_cell_length(3) = magnitude(recp_lattice(3))

    ! ab: get Ewald convergence parameter
    if (pub_dftb_ewald_parameter > 0.0_DP) then
       eta = pub_dftb_ewald_parameter
    else
       call ewald_get_eta(real_lattice, recp_lattice, tol, eta)
    end if

    ! ab: get cutoff for real and reciprocal space sums
    call ewald_real_cutoff(eta, tol, real_cutoff)
    call ewald_recp_cutoff(eta, volume, tol, recp_cutoff)

    ! ab: generate lattice points
    call simulation_cell_lattice_points_init(real_lattice, real_cutoff, real_points)
    call simulation_cell_lattice_points_init(recp_lattice, recp_cutoff, recp_points)

#ifdef SCALAPACK
    ! ab: block-cyclically distributed matrices
    call blacs_gridinfo(aa%blacs_desc(2), nprow, npcol, myrow, mycol)

    do lj = 0, aa%blacs_ncol - 1
       jat=((lj/aa%blacs_nb)*npcol+mycol)*aa%blacs_nb+mod(lj,aa%blacs_nb)+1

       do li = 0, aa%blacs_ld - 1
          iat=((li/aa%blacs_nb)*nprow+myrow)*aa%blacs_nb+mod(li,aa%blacs_nb)+1

          if (jat==nat+1) then
             if (iat/=nat+1) then
                call dense_put_element(1.0_DP, aa, iat, jat)
             end if
             cycle
          else
             if (iat==nat+1) then
                call dense_put_element(1.0_DP, aa, iat, jat)
                cycle
             end if
          end if

          zj = mdl%elements(jat)%atomic_number
          rj = mdl%elements(jat)%centre

          zi = mdl%elements(iat)%atomic_number
          ri = mdl%elements(iat)%centre

          gij = sqrt(mdl%dftb_par%z(zi)%alpg**2+mdl%dftb_par%z(zj)%alpg**2)
          if (gij > SAFE_DIV_EPS) then
             gij = 1.0_DP/gij
          else
             call utils_abort(myself//': Division by zero in calculating gij.')
          end if

          dr = ri - rj
          rij = magnitude(dr)

          gterm = ewald_recp_sum(dr,recp_points(2:),volume,eta)

          ! ab: temporary fix around to replicate the bug in xTB
          if (pub_dftb_ewald_replicate_xtb) then
             gterm = gterm - PI/volume/eta**2
             if (iat == jat) then
                dr = real_points(2)
             else
                rij = 0.0_DP
             end if
          end if
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          rterm = ewald_real_sum(dr,real_points,gij) &
                - ewald_real_sum(dr,real_points,eta)

          if (iat == jat) then
             aij = gterm + rterm + &
                  mdl%dftb_par%z(zi)%gam + 2.0_DP*(gij-eta)/sqrt(PI)
          else
             aij = gterm + rterm + &
                  2.0_DP*PI/3.0_DP/volume * rij**2
          end if
          call dense_put_element(aij, aa, iat, jat)

          if (pub_dftb_calculate_force) then
             if (iat /= jat) then
                ad = dftb_nsum(dr,gij,ewald_real_derivative,real_points) &
                     -dftb_nsum(dr,eta,ewald_real_derivative,real_points) &
                     +ewald_recp_sum_dr(dr,recp_points(2:),volume,eta) &
                     +4.0_DP*PI/3.0_DP/volume*(/dr%X,dr%Y,dr%Z/)
                do ic = 1, 3
                   call dense_put_element(ad(ic), daij(ic), iat, jat)
                end do
             end if
          end if
       end do
    end do

#else
    ! ab: undistributed matrix (allocated on all procs)
    aa%dmtx(nat+1,:) = 1.0_DP
    aa%dmtx(:,nat+1) = 1.0_DP
    aa%dmtx(nat+1, nat+1) = 0.0_DP

    do local_jat = 1, mdl%regions(1)%par%num_atoms_on_proc(pub_my_proc_id)
       global_jat = mdl%regions(1)%par%first_atom_on_proc(pub_my_proc_id) + local_jat-1
       jat = mdl%regions(1)%par%orig_atom(global_jat)
       zj = mdl%elements(jat)%atomic_number
       rj = mdl%elements(jat)%centre

       do iat = 1, nat
          zi = mdl%elements(iat)%atomic_number
          ri = mdl%elements(iat)%centre

          gij = sqrt(mdl%dftb_par%z(zi)%alpg**2+mdl%dftb_par%z(zj)%alpg**2)
          if (gij > SAFE_DIV_EPS) then
             gij = 1.0_DP/gij
          else
             call utils_abort(myself//': Division by zero in calculating gij.')
          end if

          dr = ri - rj
          rij = magnitude(dr)

          gterm = ewald_recp_sum(dr,recp_points(2:),volume,eta)

          ! ab: temporary fix around to replicate the bug in xTB
          if (pub_dftb_ewald_replicate_xtb) then
             gterm = gterm - PI/volume/eta**2
             if (iat == jat) then
                dr = real_points(2)
             else
                rij = 0.0_DP
             end if
          end if
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          rterm = ewald_real_sum(dr,real_points,gij) &
                - ewald_real_sum(dr,real_points,eta)

          if (iat == jat) then
             aa%dmtx(iat, jat) = gterm + rterm + &
                  mdl%dftb_par%z(zi)%gam + 2.0_DP*(gij-eta)/sqrt(PI)
          else
             aa%dmtx(iat, jat) = gterm + rterm + &
                  2.0_DP*PI/3.0_DP/volume * rij**2
          end if

          if (pub_dftb_calculate_force) then
             if (iat /= jat) then
                ad = dftb_nsum(dr,gij,ewald_real_derivative,real_points) &
                     -dftb_nsum(dr,eta,ewald_real_derivative,real_points) &
                     +ewald_recp_sum_dr(dr,recp_points(2:),volume,eta) &
                     +4.0_DP*PI/3.0_DP/volume*(/dr%X,dr%Y,dr%Z/)
                do ic = 1, 3
                   daij(ic)%dmtx(iat,jat) = ad(ic)
                end do
             end if
          end if
       end do
    end do

    call comms_reduce('SUM', aa%dmtx)
    if (pub_dftb_calculate_force) then
       do ic = 1, 3
          call comms_reduce('SUM', daij(ic)%dmtx)
       end do
    end if

#endif

    ! ab: deallocate lattice points
    call simulation_cell_lattice_points_exit(real_points)
    call simulation_cell_lattice_points_exit(recp_points)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving '//myself

  end subroutine dftb_gfn0_ies_pbc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_qc_print(tag,energy,nat,force)
    !==========================================================================!
    ! Wrapper for printing <QC> output.                                        !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Nov 2021.                                 !
    !==========================================================================!

    use comms, only: pub_on_root
    use utils, only: utils_qc_print

    implicit none

    ! ab: Arguments
    character(len=*), intent(in) :: tag
    real(kind=DP), intent(in)    :: energy
    integer, intent(in)          :: nat
    real(kind=DP), intent(in)    :: force(1:3,nat)

    ! ab: Local variables
    integer :: icomp, iatom
    character(len=64) :: icomp_iatom_string

    ! -------------------------------------------------------------------------

    if(pub_on_root) then
       call utils_qc_print(tag//'_energy', energy)
       do icomp = 1, 3
          do iatom = 1, nat
             write(icomp_iatom_string,'(a1,i0,a1,i0,a1)') '(', icomp , ',', &
                  iatom, ')'
             call utils_qc_print(tag//'_force'//trim(icomp_iatom_string), &
                  force(icomp,iatom))
          end do
       end do
    end if


  end subroutine dftb_qc_print

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dftb_print_energy_components(dftb_eht_energy, dftb_rep_energy, &
         dftb_srb_energy, dftb_ies_energy, dftb_disp_energy, dftb_total_energy)
    !==========================================================================!
    ! For printing energy components.                                          !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in May 2022.                                 !
    !==========================================================================!

    use comms, only: pub_on_root

    implicit none

    ! ab: Arguments
    real(kind=DP), intent(in) :: dftb_eht_energy, dftb_rep_energy, &
         dftb_srb_energy, dftb_ies_energy, dftb_disp_energy, dftb_total_energy

    ! -------------------------------------------------------------------------
    if (pub_on_root) then
       write(stdout,'(/a)') '===================================================&
            &============================='
       write(stdout,'(11x,a)') '---------------- ENERGY COMPONENTS (Eh) &
            &----------------'
       write(stdout,'(11x,a,f24.14,a)') &
            '| DFTB Electronic Energy     :', dftb_eht_energy,   ' |'
       write(stdout,'(11x,a,f24.14,a)') &
            '| DFTB Repulsion Energy      :', dftb_rep_energy,   ' |'
       write(stdout,'(11x,a,f24.14,a)') &
            '| DFTB SRB Correction Energy :', dftb_srb_energy,   ' |'
       write(stdout,'(11x,a,f24.14,a)') &
            '| DFTB Electrostatic (IES)   :', dftb_ies_energy,   ' |'
       write(stdout,'(11x,a,f24.14,a)') &
            '| Dispersion Correction      :', dftb_disp_energy,  ' |'
       write(stdout,'(11x,a,f24.14,a)') &
            '| Total                      :', dftb_total_energy, ' |'

       write(stdout,'(11x,a)') '----------------------------------------&
            &----------------'
       write(stdout,'(a/)') '===================================================&
            &============================='
   end if

   end subroutine dftb_print_energy_components

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine dftb_print_force(forces,component_name,title,mdl)
    !==========================================================================!
    ! Adapted from internal_print_forces() from forces_mod for printing forces.!
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in May 2022.                                 !
    !==========================================================================!

     use comms, only: pub_on_root
     use model_type, only: MODEL

     ! Arguments
     type(MODEL), intent(in) :: mdl
     real(kind=DP), intent(in) :: forces(1:3, 1:mdl%nat)
     character(*), intent(in)  :: component_name, title

     ! Local variables
     integer :: bef, aft, w, l1, l2, atom

     !--------------------------------------------------------------------------
     w = 58
     l1 = len(title)
     l2 = len(component_name)

     if (pub_on_root) then
        write(stdout,'(a)') ' '
        bef = (w-l1-2)/2
        aft = (w-l1-2)/2 + modulo(w-l1-2,2)
        write(stdout,'(a)') repeat('*',bef)//' '//title//' '//repeat('*',aft)
        bef = (w-l2-2)/2
        aft = (w-l2-2)/2 + modulo(w-l2-2,2)
        write(stdout,'(a)') repeat('*',bef)//' '//component_name// &
             ' '//repeat('*',aft)
        write(stdout,'(a)') '*'//repeat(' ',w-2)//'*'
        write(stdout,'(a)') '* Element  Atom         &
             &Cartesian components (Eh/a)      *'
        write(stdout,'(a)') '* '//repeat('-',w-4)//' *'
        write(stdout,'(a)') '*                       x            &
             &y            z      *'
        write(stdout,'(a)') '*'//repeat(' ',w-2)//'*'
        do atom=1,mdl%nat
           write(stdout,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
                mdl%elements(atom)%symbol,' ',&
                atom,'   ',forces(:,atom),' *'
        end do
        write(stdout,'(a)') '*'//repeat(' ',w-2)//'*'
        write(stdout,'(a)') repeat('*',w)
        write(stdout,'(a,3f13.8)') 'TOTAL:           ',&
             sum(forces(1,1:mdl%nat)),  &
             sum(forces(2,1:mdl%nat)),  &
             sum(forces(3,1:mdl%nat))
     end if

   end subroutine dftb_print_force

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module dftb
