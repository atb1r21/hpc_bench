! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                  Energy Eigenstate module                      !
!                                                                !
! This module analyses the energy eigenstates obtained by        !
! diagonalising the hamiltonian in the NGWF basis, after a       !
! converged ONETEP calculation has completed.                    !
!----------------------------------------------------------------!
! Module created from parts of properties_mod.F90 by Nicholas    !
! Hine on 16/11/2011.                                            !
! The original routines were written by Chris-Kriton Skylaris,   !
! Peter Haynes, Nicholas Hine, and others, 2002-2011.            !
! Reorganised by Nicholas Hine in August 2015.                   !
!================================================================!

module eigenstates

  use constants, only: DP
  use function_basis, only: FUNC_BASIS
  use rundat, only: pub_debug_on_root
  use simulation_cell, only: CELL_INFO
  implicit none

  private

  public :: eigenstates_calculate
  public :: eigenstates_create_storage
  public :: eigenstates_destroy_storage
  public :: eigenstates_plot_mo
  public :: eigenstates_orb_title_and_filename
  public :: eigenstates_init_pdos
  public :: eigenstates_cleanup_pdos
  public :: eigenstates_lumo_search

  ! JA: sw basis for pDOS.
  type(func_basis), allocatable, save :: lcao_basis(:)
  type(CELL_INFO), save :: fine_cell

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !======================================================================!
  ! Helper routine to create appropriately-sized storage for the results !
  ! of a diagonalisation of the hamiltonian in a given NGWF basis        !
  !----------------------------------------------------------------------!
  ! Written by Nicholas Hine in August 2015.                             !
  ! Embedding version by Robert Charlton, 24/06/2018.                    !
  !======================================================================!

  subroutine eigenstates_create_storage(ngwf_basis,eigen_en,eigen_occ, &
       eigs_dens,is_cmplx)

    use dense, only: DEM, dense_create
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_num_spins
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(DEM), allocatable, intent(inout) :: eigs_dens(:)
    real(kind=DP), allocatable, intent(inout) :: eigen_en(:,:)
    real(kind=DP), allocatable, intent(inout) :: eigen_occ(:,:)
    ! agrecocmplx
    logical, optional, intent(in) :: is_cmplx

    ! Local Variables
    integer :: ierr
    integer :: is
    integer :: num
    ! agrecocmplx
    logical :: loc_cmplx

    ! agrecocmplx
    if (present(is_cmplx)) then
       loc_cmplx = is_cmplx
    else
       loc_cmplx = .false.
    end if
    num = sum(ngwf_basis(:)%num)

    ! Create dense matrix storage for eigenvectors
    allocate(eigs_dens(pub_num_spins),stat=ierr)
    call utils_alloc_check('properties_calculate','eigs_dens',ierr)
    do is=1,pub_num_spins
       ! agrecocmplx: complex eigenstates for complex matrices
       call dense_create(eigs_dens(is),num,num,iscmplx=loc_cmplx)
    end do

    ! Create arrays to store eigenvalues and kernel occupations
    ! agrecocmplx: these should stay real anyway
    allocate(eigen_en(num,pub_num_spins),stat=ierr)
    call utils_alloc_check('eigenstates_create_storage','eigen_en',ierr)
    allocate(eigen_occ(num,pub_num_spins),stat=ierr)
    call utils_alloc_check('eigenstates_create_storage','eigen_occ',ierr)

  end subroutine eigenstates_create_storage


  !======================================================================!
  ! Helper routine to destroy storage previously created for the results !
  ! of a diagonalisation of the hamiltonian in a given NGWF basis        !
  !----------------------------------------------------------------------!
  ! Written by Nicholas Hine in August 2015.                             !
  !======================================================================!

  subroutine eigenstates_destroy_storage(eigen_en,eigen_occ, &
       eigs_dens)

    use dense, only: DEM, dense_destroy
    use rundat, only: pub_num_spins
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(DEM), allocatable, intent(inout) :: eigs_dens(:)
    real(kind=DP), allocatable, intent(inout) :: eigen_en(:,:)
    real(kind=DP), allocatable, intent(inout) :: eigen_occ(:,:)

    ! Local Variables
    integer :: ierr
    integer :: is

    ! Destroy storage for eigenvalues and kernel occupations
    deallocate(eigen_occ,stat=ierr)
    call utils_dealloc_check('eigenstates_destroy_storage','eigen_occ',ierr)
    deallocate(eigen_en,stat=ierr)
    call utils_dealloc_check('eigenstates_destroy_storage','eigen_en',ierr)

    ! Destroy dense storage for eigenvectors
    do is=pub_num_spins,1,-1
       call dense_destroy(eigs_dens(is))
    end do
    deallocate(eigs_dens,stat=ierr)
    call utils_dealloc_check('eigenstates_destroy_storage','eigs_dens',ierr)

  end subroutine eigenstates_destroy_storage

  subroutine eigenstates_calculate(denskern, ham, rep, &
       ngwf_basis, mdl, hfxstate, ham_type, eigs_dens, eigen_en, eigen_occ, &
       nl_projectors, proj_basis, core_wvfns, core_basis, cond_n_occ, &
       cross_overlap_cj, cond_denskern)

    !======================================================================!
    ! This subroutine generates and diagonalises a Hamiltonian and         !
    ! performs various kinds of analysis with the resulting eigenvectors   !
    ! and eigenenergies.                                                   !
    !----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 10/05/2006 as part of the routine!
    ! properties_calculate.                                                !
    ! Modified by Nicholas Hine in July 2009 to use function_basis and     !
    ! SPAM3 routines.                                                      !
    ! Modified by Nicholas Hine in October 2009 to make use of ScaLAPACK   !
    ! for diagonalisation of large systems.                                !
    ! Modified by Nicholas Hine in February 2010 to remove #ifdef's by     !
    ! moving code relating to eigensolving to dense_mod.                   !
    ! Modified by Nicholas Hine in October 2010 to use NGWF_REP.           !
    ! Modified by Laura Ratcliff in October 2010 for conduction            !
    ! calculations.                                                        !
    ! Simplified by Nicholas Hine in April 2011 to reduce branching of     !
    ! code in spectra calculations.                                        !
    ! Moved to eigenstates_mod by Nicholas Hine in November 2011.          !
    ! Modified for embedding by Robert Charlton, 25/06/2018.               !
    ! Modified to shift LDOS outputs setting the Fermi energy as the zero  !
    ! in the energy scale by Arihant Bhandari in Jan, 2021.                !
    !======================================================================!

    use augmentation, only: augmentation_overlap
    use comms, only: pub_on_root
    use constants, only: stdout, HARTREE_IN_EVS
    use datatypes, only: FUNCTIONS, data_functions_alloc, data_functions_dealloc
    use dense, only: DEM, dense_create, dense_destroy, dense_eigensolve, &
         dense_product, dense_convert, dense_get_element, dense_get_col
    use eels, only: eels_calculate
    use ensemble_dft, only: edft_zerok_fermi
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE, hf_exchange_calculate
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use parallel_strategy, only: PARAL_INFO
    use population, only: population_analysis_mulliken
    use projectors, only: PROJECTOR_SET
    use optics, only: optics_calculate_spectra
    use rundat, only: pub_homo_dens_plot, pub_lumo_dens_plot, pub_homo_plot, &
         pub_lumo_plot, pub_dos_smear, pub_ldos_ngroups, &
         pub_num_eigenvalues, pub_aug, pub_cond_calculate, &
         pub_spectra_calculate, pub_hubbard, &
         pub_eels_calculate, pub_devel_code, pub_eigensolver_orfac, &
         pub_eigensolver_abstol, pub_num_spins, pub_pdos_max_l, pub_spin, &
         pub_popn_mulliken_partial, pub_num_kpoints, PUB_1K, pub_edft, pub_use_hfx, &
         pub_active_region, pub_use_activehfx
    use sparse, only: SPAM3, sparse_create, sparse_destroy, &
         sparse_outer_product, sparse_trace, sparse_scale
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_transpose, sparse_embed_transpose_structure, &
         sparse_embed_outer_product, SPAM3_EMBED_ARRAY, sparse_embed_trace
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_banner
    use xc, only: pub_hfxfraction, xc_embed_swap_functional

    implicit none

    ! Arguments
    type(NGWF_REP), intent(in) :: rep
    type(NGWF_HAM), intent(inout) :: ham
    type(SPAM3_EMBED_ARRAY), intent(in) :: denskern
    type(FUNC_BASIS), intent(in) :: ngwf_basis(rep%nsub)
    type(FUNC_BASIS), optional, intent(in) :: proj_basis(rep%nsub)
    type(PROJECTOR_SET), optional, intent(inout) :: nl_projectors(rep%nsub)
    type(MODEL), intent(inout) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate

    ! lr408: String defining which type of Hamiltonian this is for
    character(len=*), intent(in) :: ham_type
    ! smmd: dense matrix for eigenvectors
    type(DEM), intent(inout) :: eigs_dens(pub_num_spins)
    ! smmd: hamiltonian eigenvalues
    real(kind=DP), intent(inout) :: eigen_en(sum(ngwf_basis(:)%num),pub_num_spins)
    ! smmd: denskern occupancies
    real(kind=DP), intent(inout) :: eigen_occ(sum(ngwf_basis(:)%num),pub_num_spins)
    ! lr408: Optional integer defining the total number of states to include for optical spectra
    integer, optional, intent(in) :: cond_n_occ(:,:)
    type(FUNC_BASIS), intent(in), optional :: core_basis(rep%nsub)
    type(PROJECTOR_SET), intent(inout), optional :: core_wvfns(rep%nsub)
    ! ndmh: Optional matrix arg for joint basis calculations, overlap of joint and cross NGWFs
    type(SPAM3_EMBED), optional, intent(in) :: cross_overlap_cj
    ! ndmh: Optional DKERN arg for joint basis calculations, giving joint denskern
    type(SPAM3_EMBED), optional, intent(in) :: cond_denskern(pub_num_spins)

    ! Local Variables
    type(SPAM3_EMBED)  :: mo_kern(1)   ! density kernel of one molecular orbital
    type(SPAM3_EMBED)  :: hub_proj     ! Hubbard projection in SPAM3 format
    type(SPAM3_EMBED)  :: aug_overlap  ! Augmentation part of overlap matrix
    type(SPAM3_EMBED)  :: cross_olap_trans        ! temporary cross overlap transpose
    type(SPAM3_EMBED)  :: kt,ukt                  ! Density kernel projected into current rep
    type(SPAM3_EMBED)  :: cross_olap_jc_trans     ! temporary cross overlap transpose
    type(SPAM3_EMBED)  :: kc_tcj,ujc_kc_tcj       ! Density kernel projected into current rep
    type(DEM) :: buffer_dens     ! dense matrix for hamiltonian and density kernel
    type(DEM) :: overlap_dens    ! dense matrix for overlap
    type(DEM) :: hub_proj_dens   ! dense matrix for Hubbard projection
    ! agrecocmplx: use FUNCTIONS type to switch between real and complex case
    type(FUNCTIONS) :: mo_coeff    ! coefficients of MO's
    type(PARAL_INFO), pointer :: par   ! parallel strategy
    real(kind=DP), allocatable, dimension(:,:) :: eigen_cond_occ ! occupancies from cond_denskern
    integer :: num                 ! number of NGWFs to diagonalise over (may be joint)
    integer :: ierr                ! memory allocation error flag
    integer :: is                  ! spin loop counter
    integer :: orbital_index       ! orbital counting index
    integer :: first_orbital,last_orbital
    integer :: num_opt_states      ! number of states to calculate optical properties for
    integer :: homo_max            ! maximum index of occupied orbital to plot
    integer :: homo_min            ! minimum index of occupied orbital to plot
    integer :: lumo_max            ! maximum index of virtual orbital to plot
    integer :: lumo_min            ! minimum index of virtual orbital to plot
    character(len=256) :: output_file  ! file names
    character(len=256) :: title_line   ! title line in output file
    ! agrecocmplx
    logical :: loc_cmplx
    complex(kind=DP) :: temp_cmplx
    ! kkbd ! agrecocmplx
    ! The second dimension should be kpoints
    real(kind=dp), dimension(pub_num_spins,1) :: s_occ ! Real version of rep%n_occ
    integer, dimension(pub_num_spins,1) :: OK_homo_orb, OK_lumo_orb
    real(kind=dp), dimension(pub_num_spins) :: zerok_fermi
    integer :: iorb, jorb
    real(kind=dp) :: e_split
    type(SPAM3)  :: lumo_kern(pub_num_spins)   ! density kernel of one molecular orbital
    type(SPAM3)  :: homo_kern(pub_num_spins)   ! density kernel of one molecular orbital
    type(FUNCTIONS) :: homo_coeff, lumo_coeff    ! coefficients of MO's
    real(kind=DP), allocatable, dimension(:,:) :: energies ! buffer for shifted energies
    real(kind=DP) :: fermi_e  ! Fermi energy

    !    integer :: pub_pdos_max_l = 1

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering eigenstates_calculate.'

    ! Start timer
    call timer_clock("eigenstates_calculate",1)

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine eigenstates_calculate not ready yet for more&
         & than one k-point.')

    ! agrecocmplx
    loc_cmplx = rep%overlap%iscmplx

    ! lr408: Print diagonalisation headers
    if (pub_cond_calculate) then
       if (pub_on_root) then
          if (ham_type == 'valence') then
             write(stdout,'(/a)') ''
             write(stdout,'(/a)') utils_banner('=', 'Valence diagonalisation')
          else if (ham_type == 'proj') then
             write(stdout,'(/a)') ''
             write(stdout,'(/a)') utils_banner('=', &
                  'Projected conduction diagonalisation')
          else if (ham_type == 'cond') then
             write(stdout,'(/a)') ''
             write(stdout,'(/a)') utils_banner('=', &
                  'Unprojected conduction diagonalisation')
          else if (ham_type == 'joint') then
             write(stdout,'(/a)') ''
             write(stdout,'(/a)') utils_banner('=', &
                  'Joint valence and conduction diagonalisation')
          end if
       end if
    end if

    ! rc2013: take sum of all NGWF bases
    num = sum(ngwf_basis(:)%num)

    ! ndmh: Allocate storage for coefficient vector
    ! agrecocmplx: allocate using appropriate routine
    call data_functions_alloc(mo_coeff,num,iscmplx=loc_cmplx)
    call data_functions_alloc(homo_coeff,num,iscmplx=loc_cmplx)
    call data_functions_alloc(lumo_coeff,num,iscmplx=loc_cmplx)

    ! Create dense matrices for Hamiltonian and overlap
    call dense_create(buffer_dens,num,num,iscmplx=loc_cmplx)
    call dense_create(overlap_dens,num,num,iscmplx=loc_cmplx)


    if ((ham_type=='valence').or.(ham_type=='cond').or. &
         (ham_type == 'proj').or.(ham_type=='aux')) then

       ! ndmh: Create matrix for <chi_c|phi_b> K^ba <phi_a|chi_d>
       call sparse_embed_create(kt,denskern%m(1,1),rep%overlap)
       call sparse_embed_create(ukt,rep%overlap,kt)

    else if (ham_type=='joint') then

       ! ndmh: Create matrix for <xi_c|phi_b> K^ba <phi_a|xi_d>
       call sparse_embed_create(kt,denskern%m(1,1),rep%cross_overlap)
       call sparse_embed_transpose_structure(cross_olap_trans%structure, &
            rep%cross_overlap)
       call sparse_embed_create(cross_olap_trans,iscmplx=loc_cmplx)
       call sparse_embed_transpose(cross_olap_trans,rep%cross_overlap)
       call sparse_embed_create(ukt,cross_olap_trans,kt)

       ! create matrix for <xi_c|chi_b> K_c^ba <chi_a|xi_d>
       call sparse_embed_create(kc_tcj,cond_denskern(1),cross_overlap_cj)
       call sparse_embed_transpose_structure(cross_olap_jc_trans%structure, &
            cross_overlap_cj)
       call sparse_embed_create(cross_olap_jc_trans,iscmplx=loc_cmplx)
       call sparse_embed_transpose(cross_olap_jc_trans,cross_overlap_cj)
       call sparse_embed_create(ujc_kc_tcj,cross_olap_jc_trans,kc_tcj)

       ! ndmh: Create storage for occupancies of joint orbitals according
       ! ndmh: to conduction kernel
       allocate(eigen_cond_occ(num,pub_num_spins),stat=ierr)
       call utils_alloc_check('eigenstates_calculate','eigen_cond_occ',ierr)

    end if

    ! ddor: Create matrices for DFT+U projection
    if (pub_hubbard) then
       call dense_create(hub_proj_dens,num,num,iscmplx=loc_cmplx)
       call sparse_embed_create(hub_proj,rep%hub_overlap,rep%hub_overlap_t)
    endif

    ! ndmh: Create matrix to hold MO kernel
    if ((pub_aug.and.(pub_homo_plot >= 0 .or. pub_lumo_plot >= 0)).or. &
         (pub_homo_dens_plot >= 0 .or. pub_lumo_dens_plot >= 0).or. &
         pub_popn_mulliken_partial .or. pub_use_hfx .or. pub_use_activehfx) then
       mo_kern(1)%structure = 'K'//rep%postfix
       ! agrecocmplx
       call sparse_embed_create(mo_kern(1),iscmplx=loc_cmplx)
    end if

    ! ndmh: Create matrix for aug part of overlap
    if (pub_aug.and.(pub_homo_plot >= 0 .or. pub_lumo_plot >= 0)) then
       call sparse_embed_create(aug_overlap,rep%overlap)
       ! rc2013: EMBED_FIX! not compatible multiple subsytems
       call augmentation_overlap(aug_overlap%p,mdl%pseudo_sp,mdl%paw_sp, &
            rep%sp_overlap%p)
    end if

    ! pdh: loop over spins
    do is=1,pub_num_spins

       ! ndmh: diagonalise hamiltonian for eigenvalues
       call dense_convert(overlap_dens,rep%overlap)
       ! lr408: Choose correct Hamiltonian for conduction calculation
       if (ham_type=='cond') then
          call dense_convert(buffer_dens,ham%cond_non_proj_ham(is))
       else
          call dense_convert(buffer_dens,ham%ham(is))
       end if
       call dense_eigensolve(num,eigen_en(:,is),buffer_dens, &
            overlap_dens,1,eigs_dens(is),pub_eigensolver_orfac, &
            pub_eigensolver_abstol)

       ! ndmh: now calculate occupancy of each state by projection onto
       ! ndmh: valence density matrix for this spin
       if ((ham_type=='valence').or.(ham_type=='cond').or. &
            (ham_type == 'proj').or.(ham_type=='aux')) then
          call sparse_embed_product(kt,denskern%m(is,PUB_1K),rep%overlap)
          call sparse_embed_product(ukt,rep%overlap,kt)
       else if (ham_type=='joint') then
          call sparse_embed_product(kt,denskern%m(is,PUB_1K),rep%cross_overlap)
          call sparse_embed_product(ukt,cross_olap_trans,kt)
          call sparse_embed_product(kc_tcj,cond_denskern(is),cross_overlap_cj)
          call sparse_embed_product(ujc_kc_tcj,cross_olap_jc_trans,kc_tcj)
       end if
       call dense_convert(overlap_dens,ukt)
       call dense_product(buffer_dens,eigs_dens(is),overlap_dens,opA='T')
       call dense_product(overlap_dens,buffer_dens,eigs_dens(is))
       do orbital_index=1,num
          ! agrecocmplx: need temporary complex to get complex element
          if (loc_cmplx) then
             call dense_get_element(temp_cmplx,overlap_dens, &
                  orbital_index,orbital_index)
             ! agrecocmplx: occupancies are real anyway
             eigen_occ(orbital_index,is) = real(temp_cmplx,kind=DP)
          else
             call dense_get_element(eigen_occ(orbital_index,is),overlap_dens, &
                  orbital_index,orbital_index)
          end if
       end do

       ! ndmh: calculate 'occupancy' of joint basis states via projection
       ! ndmh: onto conduction density kernel, to print as extra column
       if (ham_type=='joint') then
          call dense_convert(overlap_dens,ujc_kc_tcj)
          call dense_product(buffer_dens,eigs_dens(is),overlap_dens,opA='T')
          call dense_product(overlap_dens,buffer_dens,eigs_dens(is))
          do orbital_index=1,num
             ! agrecocmplx: need temporary complex to get complex element
             if (loc_cmplx) then
                call dense_get_element(temp_cmplx,overlap_dens, &
                     orbital_index,orbital_index)
                eigen_cond_occ(orbital_index,is) = real(temp_cmplx,kind=DP)
             else
                ! agrecocmplx: occupancies are real anyway
                call dense_get_element(eigen_cond_occ(orbital_index,is),overlap_dens, &
                     orbital_index,orbital_index)
             end if
          end do
       end if

    end do

    if (pub_edft) then
       ! kkbd: If we've done an EDFT run we might have a non-integer net spin
       !       that also doesn't equal its original value!
       do is=1,pub_num_spins
          call sparse_embed_trace(s_occ(is,PUB_1K),denskern%m(is,PUB_1K),rep%overlap)
       end do

       call edft_zerok_fermi(zerok_fermi, s_occ(:,PUB_1K), eigen_en)
    else
       s_occ = real(rep%n_occ,kind=dp)
    end if

    ! cks: print out orbital energies and denskern occupancies
    if (pub_num_eigenvalues > 0 .or. pub_cond_calculate) then
       if (ham_type/='joint') then
          call eigenstates_print_ens_occs(eigen_en, eigen_occ, &
               s_occ(:,PUB_1K), num)
       else
          call eigenstates_print_ens_occs(eigen_en, eigen_occ, &
               s_occ(:,PUB_1K), &
               num,eigen_cond_occ,cond_n_occ(:,PUB_1K))
       end if
    end if

    do is=1,pub_num_spins

       if (pub_edft) then
          do iorb=1,num
             if (eigen_en(iorb,is) > zerok_fermi(is)) then
                ! kkbd: previous orbital was HOMO of this spin channel
                !       and this one is LUMO
                lumo_min = iorb
                homo_max = iorb - 1
                exit
             end if
          end do
       else
          lumo_min = nint(s_occ(is,PUB_1K)) + 1
          homo_max = nint(s_occ(is,PUB_1K))
       end if

       ! kkbd: Keep 0K homo and lumo for later
       OK_lumo_orb(is,:) = lumo_min
       OK_homo_orb(is,:) = homo_max

       ! ndmh: find range of states to plot
       if (pub_cond_calculate .and. ham_type == 'proj') then
          lumo_min = 1
          lumo_max = 1  + pub_lumo_plot
       else
          lumo_max = lumo_min + pub_lumo_plot
       end if
       if (lumo_max > num) lumo_max = num
       homo_min = homo_max - pub_homo_plot
       if (homo_min <= 0) homo_min = 1

       ! cks: Create canonical orbital plot files

       if(mdl%nsub == 1) then
          call internal_plot_orbitals
       else
           if(pub_on_root) write(stdout,'(/a)') &
               'WARNING: plotting molecular orbitals requested with embedding, &
               & but this combination is not implemented/tested. Skipping MO plot.'
       end if

    end do

    ! ab: calculate Fermi level for shifting DOS outputs
    call eigenstates_fermi_level(eigen_en, s_occ(:,PUB_1K), num, fermi_e)

    allocate(energies(num,pub_num_spins),stat=ierr)
    call utils_alloc_check('eigenstates_calculate','energies',ierr)

    ! cks: set fermi energy as the zero of energy
    ! ndmh: or don't, in the case of COND calculations
    if (.not.pub_cond_calculate) then
       energies = eigen_en - fermi_e
    else
       energies = eigen_en
    end if

    ! ndmh: Calculate and output LDOS for this spin
    if ((pub_dos_smear > 0.0_DP).and.(pub_ldos_ngroups>0)) then

       do is=1,pub_num_spins

          ! ndmh: Standard LDOS, projected onto atoms
          call dense_convert(overlap_dens,rep%overlap)
          call dense_product(buffer_dens,overlap_dens,eigs_dens(is))
          call eigenstates_ldos_gp(energies(:,is),eigs_dens(is),buffer_dens, &
               is,mdl,ngwf_basis,.false.,ham_type)

          ! ddor: Compute projection of LDOS onto Hubbard manifold
          if (pub_hubbard) then
             call sparse_embed_product(hub_proj,&
                  rep%hub_overlap,rep%hub_overlap_t)
             call dense_convert(hub_proj_dens,hub_proj)
             call dense_product(buffer_dens,hub_proj_dens,eigs_dens(is))
             call eigenstates_ldos_gp(energies(:,is),eigs_dens(is),buffer_dens, &
                  is,mdl,ngwf_basis,.true.,ham_type)
          endif
       end do
    end if

    ! lr408: Calculate matrix elements for spectra - only in conduction
    ! lr408: calculation and for valence and joint basis sets only
    if ((pub_spectra_calculate .and. pub_cond_calculate).and. &
         present(nl_projectors).and.present(proj_basis)) then
       if ((ham_type=='valence'.or.ham_type=='joint') .and. mdl%nsub==1) then
          ! kkbd: TODO: Is this supposed to be 0K HOMO?
          num_opt_states = maxval(rep%n_occ(:,PUB_1K))
          if (present(cond_n_occ)) num_opt_states = &
               maxval(rep%n_occ(:,PUB_1K)+cond_n_occ(:,PUB_1K))
          if (num_opt_states>num) num_opt_states = num
          do is=1,pub_num_spins
             call optics_calculate_spectra(is,eigs_dens(is),eigen_en(:,is), &
                  rep,ngwf_basis,nl_projectors,proj_basis,mdl, &
                  num_opt_states)
          end do
       else if (ham_type=='valence'.or.ham_type=='joint' &
            .and. mdl%nsub .gt. 1) then
          ! rc2013: not compatible with subsystem calculations yet
          if(pub_on_root) write(stdout,'(/a)') &
               'WARNING: Optical spectra calculation requested with embedding, &
               &but this combination is not implemented/tested. Skipping spectra.'
       end if
    end if

    if (pub_eels_calculate .and. pub_cond_calculate .and. &
         present(core_basis) .and. present(core_wvfns)) then
       if (ham_type=='valence'.or.ham_type=='joint' .and. mdl%nsub==1) then
          do is=1,pub_num_spins
             call eels_calculate(is,eigs_dens(is),eigen_en(:,is), &
                  rep,ngwf_basis(1),proj_basis(1),core_basis(1),core_wvfns(1),mdl)
          end do
       ! rc2013: not compatible with embedding yet
       else if (ham_type=='valence'.or.ham_type=='joint' &
            .and. mdl%nsub .gt. 1) then
          ! rc2013: not compatible with subsystem calculations yet
          if(pub_on_root) write(stdout,'(/a)') &
               'WARNING: EELS calculation requested with embedding, but this &
               &combination is not implemented/tested. Skipping EELS.'
       end if
    end if

    if ((index(pub_devel_code,'TRANSITION_DENPOT')>0).and. &
         (ham_type=='joint')) then
       do is=1,pub_num_spins
          call eigenstates_transition_denpot(eigs_dens(is), &
               ngwf_basis,rep,mdl)
       end do
    end if

    ! cks: output gamma-point Gaussian-smeared density of states
    ! cks: and CASTEP format .bands file
    if (pub_dos_smear > 0.0_DP) then
       call eigenstates_dos_gp(energies, s_occ(:,PUB_1K), &
            num, ham_type)
    end if
    if (pub_dos_smear > 0.0_DP) then
       call eigenstates_write_bands(ngwf_basis, mdl%cell, eigen_en, &
            s_occ(:,PUB_1K), ham_type)
    end if

    ! JA: write pdos weights for optados and/or Gaussian smeared
    ! JA: angular momentum projected density of states.
    if(pub_dos_smear > 0.0_dp .and. pub_pdos_max_l >= 0 .and. mdl%nsub == 1) then
!       do is=1,pub_num_spins
          call eigenstates_pdos(rep%ngwfs_on_grid, ngwf_basis, &
               proj_basis, nl_projectors, energies, eigen_occ, eigs_dens, &
               s_occ(:,PUB_1K), rep, mdl, ham_type=ham_type)
!       end do
    ! rc2013: EMBED_FIX!
    else if(pub_dos_smear > 0.0_dp .and. pub_pdos_max_l >= 0) then
       if(pub_on_root) write(stdout,'(/a)') &
            'WARNING: PDOS requested with embedding, but this &
            & combination is not implemented/tested. Skipping PDOS calculation.'
    end if

    ! ndmh: Mulliken analysis of partial charges of each orbital
    if (pub_popn_mulliken_partial) then

       ! ndmh: find range to print partial charges for
       do is=1,pub_num_spins
          first_orbital = max(OK_homo_orb(is,PUB_1K)-pub_num_eigenvalues,0)+1
          last_orbital = min(OK_homo_orb(is,PUB_1K)+pub_num_eigenvalues,num)

          do orbital_index=first_orbital,last_orbital

             ! agrecocmplx: distinguish between real and complex case
             if (loc_cmplx) then
                call dense_get_col(mo_coeff%z,eigs_dens(is),orbital_index)
                call sparse_embed_outer_product(mo_kern(1),mo_coeff%z,mo_coeff%z)
             else
                ! get M.O. eigenvector from dense matrix
                call dense_get_col(mo_coeff%d,eigs_dens(is),orbital_index)
                ! form sparse kernel for M.O. as SPAM3 matrix
                call sparse_embed_outer_product(mo_kern(1),mo_coeff%d,mo_coeff%d)
             end if

             call population_analysis_mulliken(rep%overlap, mo_kern, &
                  mdl, ngwf_basis, partial_charge_orbital=orbital_index)
          end do
       end do

    end if

    ! rc2013: calculate the singlet-triplet splitting for triplet/Hartree-Fock
    if ((pub_use_hfx .or. pub_use_activehfx) .and. pub_spin == 2) then

       ! rc2013: for now this only makes sense for 1 spin
       do is=1,pub_num_spins

          ! jcap: swap functionals for EMFT, to get pub_hfxfraction right
          call xc_embed_swap_functional(.true.)

          ! rc2013: HOMO index
          iorb = nint(s_occ(is,PUB_1K))
          if(pub_num_spins == 1) then
             iorb = nint(s_occ(is,PUB_1K))
          else if(is == 1) then
             iorb = nint(s_occ(is,PUB_1K)) - 1
          else
             iorb = nint(s_occ(is,PUB_1K)) + 1
          end if
          ! rc2013: LUMO index
          jorb = iorb + 1

          ! agrecocmplx: distinguish between real and complex case
          call sparse_create(homo_kern(is), mo_kern(1)%p,iscmplx=loc_cmplx)
          call sparse_create(lumo_kern(is), mo_kern(1)%p,iscmplx=loc_cmplx)
          if (loc_cmplx) then
             call dense_get_col(homo_coeff%z,eigs_dens(is),iorb)
             call dense_get_col(lumo_coeff%z,eigs_dens(is),jorb)
             call sparse_outer_product(homo_kern(is),homo_coeff%z,homo_coeff%z)
             call sparse_outer_product(lumo_kern(is),lumo_coeff%z,lumo_coeff%z)
          else
             ! get M.O. eigenvector from dense matrix
             call dense_get_col(homo_coeff%d,eigs_dens(is),iorb)
             call dense_get_col(lumo_coeff%d,eigs_dens(is),jorb)
             ! form sparse kernel for M.O. as SPAM3 matrix
             call sparse_outer_product(homo_kern(is),homo_coeff%d,homo_coeff%d)
             call sparse_outer_product(lumo_kern(is),lumo_coeff%d,lumo_coeff%d)
          end if
          call sparse_scale(ham%hfexchange(is)%p, 0.0_DP)

       end do
       call hf_exchange_calculate(hfxstate, ham%hfexchange, &
            rep, pub_active_region, &
            homo_kern, lumo_kern, &
            ngwf_basis(1), mdl%fftbox, mdl%cell, &
            mdl%regions(pub_active_region)%elements, .false.)

       do is=1,pub_num_spins
          if(pub_hfxfraction .ne. 0.0_DP) then
             e_split =  sparse_trace(homo_kern(is), ham%hfexchange(is)%p)/pub_hfxfraction
          else
             write(*,*) 'WARNING: pub_hfxfraction = 0!'
             e_split =  sparse_trace(homo_kern(is), ham%hfexchange(is)%p)
          endif
          if(pub_on_root) write(stdout,'(a28,f15.9,a12,i4)') &
               'Singlet-triplet splitting = ', e_split*HARTREE_IN_EVS, ' eV for spin ', is

          call sparse_destroy(lumo_kern(is))
          call sparse_destroy(homo_kern(is))
       end do

    end if

    ! ab: deallocate energies
    deallocate(energies,stat=ierr)
    call utils_dealloc_check('eigenstates_calculate','energies',ierr)

    ! ndmh: Destroy matrix for aug part of overlap
    if (pub_aug.and.(pub_homo_plot >= 0 .or. pub_lumo_plot >= 0)) then
       call sparse_embed_destroy(aug_overlap)
    end if

    ! ndmh: Destroy matrix for MO kernel
    if ((pub_aug.and.(pub_homo_plot >= 0 .or. pub_lumo_plot >= 0)).or. &
         (pub_homo_dens_plot >= 0 .or. pub_lumo_dens_plot >= 0).or. &
         pub_popn_mulliken_partial) then
       call sparse_embed_destroy(mo_kern(1))
    end if

    ! ndmh: Destroy matrices for kernel projection
    if (ham_type=='joint') then
       deallocate(eigen_cond_occ,stat=ierr)
       call utils_dealloc_check('eigenstates_calculate','eigen_cond_occ',ierr)
    end if
    call sparse_embed_destroy(ukt)
    if (ham_type=='joint') then
       call sparse_embed_destroy(cross_olap_jc_trans)
       call sparse_embed_destroy(ujc_kc_tcj)
       call sparse_embed_destroy(kc_tcj)
       call sparse_embed_destroy(cross_olap_trans)
    end if
    call sparse_embed_destroy(kt)

    ! Destroy dense matrices for Hamiltonian, denskern and overlap
    if (pub_hubbard) then
       call sparse_embed_destroy(hub_proj)
       call dense_destroy(hub_proj_dens)
    endif
    call dense_destroy(overlap_dens)
    call dense_destroy(buffer_dens)

    ! agrecocmplx: deallocate using appropriate routine
    call data_functions_dealloc(mo_coeff)

    ! cks: stop timer
    call timer_clock("eigenstates_calculate",2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving eigenstates_calculate.'

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_plot_orbitals

      use dense, only: dense_get_col

      character(len=1) :: Hermitian_conjg_op ! Hermitian conjugate operator

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering internal_plot_orbitals (eigenstates_calculate).'

      if (loc_cmplx) then
         Hermitian_conjg_op = 'C'
      else
         Hermitian_conjg_op = 'T'
      end if

      if ((pub_homo_plot >= 0 .or. pub_lumo_plot >= 0)) then

         if (pub_on_root) then
            write(stdout,'(/a)') utils_banner('=', &
                 'Writing canonical molecular orbitals')
            if (pub_num_spins == 2) write(stdout,'(a,i1)') &
                 '    Spin ',is
            write(stdout,'(a)')'| Orbital  Energy (eV)  Norm     |'
         end if

         ! ndmh: output homo and below for plotting
         do orbital_index=homo_min,homo_max
            call eigenstates_orb_title_and_filename('HOMO',is, &
                 OK_homo_orb(is,PUB_1K), &
                 orbital_index,title_line,output_file, &
                 eigen_en(orbital_index,is),ham_type)

            ! ndmh: if augmentation is active, we need to calculate the
            ! ndmh: aug part of the norm, so we need the MO 'kernel'
            if (pub_aug) then
               ! ndmh: create dense 'kernel' for MO from selected eigenvector
               call dense_product(buffer_dens,eigs_dens(is),eigs_dens(is), &
                    opA='N',opB=Hermitian_conjg_op,first_k=orbital_index, &
                    last_k=orbital_index)
               ! ndmh: convert to SPAM3
               call dense_convert(mo_kern(1),buffer_dens)
            end if

            ! agrecocmplx
            if (loc_cmplx) then
               call dense_get_col(mo_coeff%z,eigs_dens(is),orbital_index)
            else
               call dense_get_col(mo_coeff%d,eigs_dens(is),orbital_index)
            end if

            call eigenstates_plot_mo(mo_coeff, rep%ngwfs_on_grid(1), &
                 ngwf_basis(1), mdl, title_line, output_file, &
                 aug_overlap%p, mo_kern(1)%p)
         end do

         ! ndmh: output lumo and above for plotting
         do orbital_index=lumo_min,lumo_max
            if (ham_type == 'proj') then
               call eigenstates_orb_title_and_filename('LUMO',is,1, &
                    orbital_index,title_line,output_file, &
                    eigen_en(orbital_index,is),ham_type)
            else
               call eigenstates_orb_title_and_filename('LUMO',is, &
                    OK_lumo_orb(is,PUB_1K),orbital_index,title_line, &
                    output_file,eigen_en(orbital_index,is),ham_type)
            end if

            ! ndmh: if augmentation is active, we need to calculate the
            ! ndmh: aug part of the norm, so we need the MO 'kernel'
            if (pub_aug) then
               ! ndmh: create dense 'kernel' for MO from selected eigenvector
               call dense_product(buffer_dens,eigs_dens(is),eigs_dens(is), &
                    opA='N',opB=Hermitian_conjg_op,first_k=orbital_index, &
                    last_k=orbital_index)
               ! ndmh: convert to SPAM3
               call dense_convert(mo_kern(1),buffer_dens)
            end if

            ! agrecocmplx
            if (loc_cmplx) then
               call dense_get_col(mo_coeff%z,eigs_dens(is),orbital_index)
            else
               call dense_get_col(mo_coeff%d,eigs_dens(is),orbital_index)
            end if

            call eigenstates_plot_mo(mo_coeff, rep%ngwfs_on_grid(1), &
                 ngwf_basis(1), mdl, title_line, output_file, &
                 aug_overlap%p, mo_kern(1)%p)
         end do

         if (pub_on_root) write(stdout,'(a)') repeat('=',80)
      end if

      ! ndmh: find range of states to plot
      if (pub_cond_calculate .and. ham_type=='proj') then
         lumo_min = 1
         lumo_max = 1 + pub_lumo_dens_plot
      else
         lumo_min = OK_lumo_orb(is,PUB_1K)
         lumo_max = lumo_min + pub_lumo_dens_plot
      end if
      if (lumo_max > num) lumo_max = num
      homo_max = OK_homo_orb(is,PUB_1K)
      homo_min = homo_max - pub_homo_dens_plot
      if (homo_min <= 0) homo_min = 1

      if ((pub_homo_dens_plot >= 0 .or. pub_lumo_dens_plot >= 0)) then

         if (pub_on_root) then
            write(stdout,'(/a)') utils_banner('=', &
                 'Writing densities of canonical molecular orbitals')
            if (pub_num_spins == 2) write(stdout,'(a,i1)') &
                 '    Spin ',is
            write(stdout,'(a)')'| Orbital  Energy (eV)  Norm     |'
         end if

         ! cks: output homo and below for plotting
         do orbital_index=homo_min,homo_max
            ! ndmh: create dense 'kernel' for MO from selected eigenvector
            call dense_product(buffer_dens,eigs_dens(is),eigs_dens(is), &
                 opA='N',opB=Hermitian_conjg_op,first_k=orbital_index, &
                 last_k=orbital_index)
            ! ndmh: convert to SPAM3
            call dense_convert(mo_kern(1),buffer_dens)
            ! ndmh: plot orbital density scalar field
            call eigenstates_orb_title_and_filename('HOMO_density',is, &
                 OK_homo_orb(is,PUB_1K),orbital_index,title_line, &
                 output_file,eigen_en(orbital_index,is),ham_type)
            ! rc2013: this enforces 1 subsystem
            call eigenstates_plot_squared_mo( &
                 mo_kern(1)%m(1,1), rep%ngwf_overlap%p, rep%ngwfs_on_grid(1), &
                 ngwf_basis(1), mdl, title_line, output_file)
         end do

         ! cks: output lumo and above for plotting
         if (lumo_max > num) lumo_max = num
         do orbital_index=lumo_min,lumo_max
            ! ndmh: create dense 'kernel' for MO from selected eigenvector
            call dense_product(buffer_dens,eigs_dens(is),eigs_dens(is), &
                 opA='N',opB=Hermitian_conjg_op,first_k=orbital_index, &
                 last_k=orbital_index)
            ! ndmh: convert to SPAM3
            call dense_convert(mo_kern(1),buffer_dens)
            ! ndmh: plot orbital density scalar field
            if (ham_type == 'proj') then
               call eigenstates_orb_title_and_filename('LUMO_density',is,1, &
                    orbital_index,title_line,output_file, &
                    eigen_en(orbital_index,is),ham_type)
            else
               call eigenstates_orb_title_and_filename('LUMO_density',is, &
                    OK_lumo_orb(is,PUB_1K),orbital_index,title_line, &
                    output_file,eigen_en(orbital_index,is),ham_type)
            end if
            call eigenstates_plot_squared_mo( &
                 mo_kern(1)%m(1,1), rep%ngwf_overlap%p, rep%ngwfs_on_grid(1), &
                 ngwf_basis(1), mdl, title_line, output_file)
         end do

         if (pub_on_root) write(stdout,'(a)') repeat('=',80)
      end if

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving internal_plot_orbitals (eigenstates_calculate).'

    end subroutine internal_plot_orbitals

  end subroutine eigenstates_calculate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_orb_title_and_filename(id,is,origin,orb_idx, &
       title_line,output_file,eigen_en,ham_type)

    use comms, only: pub_on_root
    use constants, only: stdout, UP, HARTREE_IN_EVS
    use rundat, only: pub_cond_calculate, pub_num_spins

    implicit none

    ! Arguments
    integer, intent(in) :: orb_idx
    integer, intent(in) :: is
    character(*),intent(in) :: id
    integer,intent(in) :: origin
    character(len=256), intent(out) :: output_file
    character(len=256), intent(out) :: title_line
    real(kind=DP), intent(in) :: eigen_en
    character(len=*), intent(in) :: ham_type

    ! Local Variables
    character(len=50)  :: txt_buffer   ! text buffer
    character(len=50)  :: orb_buffer   ! text buffer
    character(len=256) :: title_buffer ! text buffer

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering eigenstates_orb_title_and_filename.'

    write(title_buffer,'(a30,i4)')'Density of canonical orbital:', &
         orb_idx

    if (orb_idx == origin) then
       if (pub_num_spins == 1) then
          write(output_file,*) '_',id
       else
          if (is == UP) then
             write(output_file,*) '_',id,'_UP'
          else
             write(output_file,*) '_',id,'_DN'
          end if
       end if
       if (pub_cond_calculate) write(output_file,*) &
            '_'//trim(ham_type)//trim(adjustl(output_file))
       if (pub_on_root) write(stdout,'(a8,f14.8,a)',advance ='no') &
            adjustl('| '//id//'  '),&
            eigen_en*HARTREE_IN_EVS,' '
       write(title_line,*)trim(adjustl(title_buffer)),' (',id,')'
    else
       if ((orb_idx-origin)>0) then
          write(txt_buffer,*) orb_idx-origin
          write(orb_buffer,*) '+',trim(adjustl(txt_buffer))
       else
          write(orb_buffer,*) orb_idx-origin
       end if
       if (pub_num_spins == 1) then
          write(output_file,*)'_',id,trim(adjustl(orb_buffer))
       else
          if (is == UP) then
             write(output_file,*)'_',id, &
                  trim(adjustl(orb_buffer)),'_UP'
          else
             write(output_file,*)'_',id, &
                  trim(adjustl(orb_buffer)),'_DN'
          end if
       end if
       if (pub_cond_calculate) write(output_file,*) &
            '_'//trim(ham_type)//trim(adjustl(output_file))
       write(title_line,*)trim(adjustl(title_buffer)),' (',id, &
            trim(adjustl(orb_buffer)),')'
       write(txt_buffer,*) '| ',id,trim(adjustl(orb_buffer))
       if (pub_on_root) write(stdout,'(a9,f13.8,a)',advance ='no')&
            adjustl(txt_buffer),eigen_en*HARTREE_IN_EVS,' '
    end if

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving eigenstates_orb_title_and_filename.'

  end subroutine eigenstates_orb_title_and_filename


  subroutine eigenstates_plot_mo(mo_coeff, &
       ngwfs_on_grid, ngwf_basis, mdl, file_header, scalar_name, &
       aug_overlap, mo_kern,print_norm)

    !==================================================================!
    ! This subroutine prepares and outputs a plotfile for the value of !
    ! a canonical molecular orbital.                                   !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 13/01/2011.                          !
    ! rab207: Added option to not print the state norm, 28/04/2014     !
    ! Modified for embedding by Joseph Prentice, September 2018        !
    !==================================================================!

    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc, data_set_to_zero
    use comms, only: pub_on_root
    use constants, only: stdout
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_product_on_grid
    use model_type, only: MODEL
    use rundat, only: pub_cube_format, pub_dx_format, pub_grd_format, pub_aug
    use sparse, only: SPAM3, sparse_trace
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, utils_assert
    use visual, only: visual_scalarfield

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNCTIONS), intent(in) :: mo_coeff         ! NGWF coeffs of orb.
    type(MODEL), intent(in) :: mdl
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid    ! ngwfs on this proc
    character(len=*), intent(in) :: file_header ! part of header line of output file
    character(len=*), intent(in) :: scalar_name ! part of file-name to output (scalafield name)
    type(SPAM3), intent(in), optional :: aug_overlap    ! sp overlaps
    type(SPAM3), intent(in), optional :: mo_kern       ! for aug norm
    logical, intent(in), optional :: print_norm

    ! Local Variables
    integer :: ierr
    real(kind=DP) :: norm, aug_norm
    type(FFTBOX_DATA) ::  orbital_std
    real(kind=DP), allocatable :: abs_orbital_std(:,:,:)
    logical :: print_norm_loc
    logical :: loc_cmplx

    if (pub_debug_on_root) &
         write(stdout,'(a)') 'DEBUG: Entering eigenstates_plot_mo.'

    ! jcap: check if we have more than one subsystem
    call utils_assert(mdl%nsub.eq.1,'Subroutine eigenstates_plot_mo is not &
         &yet ready for more than one subsystem')

    loc_cmplx = .false.
    if(ngwfs_on_grid%iscmplx) loc_cmplx = .true.

    print_norm_loc = .true.
    if (present(print_norm)) print_norm_loc = print_norm

    ! ndmh: allocate storage for orbital
    call data_fftbox_alloc(orbital_std, mdl%std_grid%ld1, mdl%std_grid%ld2, &
         mdl%std_grid%max_slabs12, loc_cmplx)

    ! lr408: set to zero here instead to allow joint val-cond basis calculations
    call data_set_to_zero(orbital_std)

    ! ndmh: calculate the M.O. on the standard grid
    call eigenstates_calculate_mo(orbital_std, mo_coeff, &
         ngwfs_on_grid, ngwf_basis, mdl%std_grid, mdl%cell, mdl%fftbox)

    ! jme: in the complex case, compute the absolute value of the M.O.
    if (loc_cmplx) then
       allocate(abs_orbital_std(mdl%std_grid%ld1, mdl%std_grid%ld2, &
            mdl%std_grid%max_slabs12), stat=ierr)
       call utils_alloc_check('eigenstates_plot_mo', 'abs_orbital_std', ierr)
       abs_orbital_std(:,:,:) = abs(orbital_std%z)
    end if

    ! ndmh: get norm of M.O.
    if (print_norm_loc) then
       if (loc_cmplx) then
          norm = integrals_product_on_grid(mdl%std_grid, &
               abs_orbital_std, abs_orbital_std)
       else
          norm = integrals_product_on_grid(mdl%std_grid, &
               orbital_std%d, orbital_std%d)
       end if
    else
       norm = 0.0_DP
    endif

    ! ndmh: augment norm with augmentation part in PAW/USP formalism
    if (pub_aug) then
       if ((.not.present(mo_kern)).or.(.not.present(aug_overlap))) then
          call utils_abort('Error in eigenstates_plot_mo: mo_kern or &
               &aug_overlap was not provided, yet augmentation is active')
       end if
       aug_norm = sparse_trace(mo_kern,aug_overlap,'C')
       norm = norm + aug_norm
    end if

    if (pub_on_root.and.print_norm_loc) &
         write(stdout,'(f9.6,a)',advance='no') norm, ' | '

    ! ndmh: plot this orbital on standard grid
    if (loc_cmplx) then
       ! jme: in the complex case, for now, we only plot
       ! the absolute value of the orbital
       call visual_scalarfield(abs_orbital_std, mdl%std_grid, mdl%cell, &
            file_header, scalar_name, mdl%elements, 1.0_DP)  ! input
    else
       call visual_scalarfield(orbital_std%d, mdl%std_grid, mdl%cell, &
            file_header, scalar_name, mdl%elements, 1.0_DP)  ! input
    end if

    ! ndmh: print end of line if no message was written
    if ((.not.pub_grd_format).and.(.not.pub_cube_format).and. &
         (.not.pub_dx_format).and.pub_on_root) write(stdout,*)

    ! jme: deallocate storage for absolute value
    if (loc_cmplx) then
       deallocate(abs_orbital_std, stat=ierr)
       call utils_dealloc_check('eigenstates_plot_mo', 'abs_orbital_std', ierr)
    end if

    ! ndmh: deallocate storage for orbital
    call data_fftbox_dealloc(orbital_std)

    if (pub_debug_on_root) &
         write(stdout,'(a)') 'DEBUG: Leaving eigenstates_plot_mo.'

  end subroutine eigenstates_plot_mo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_calculate_mo(orbital_std, mo_coeff, ngwfs_on_grid, &
       ngwf_basis, grid, cell, fftbox)

    !==================================================================!
    ! This subroutine calculates a canonical molecular orbital from    !
    ! coefficients of each NGWF, returned by a diagonalisation of the  !
    ! Hamiltonian.                                                     !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 13/01/2011.                          !
    ! Modified for embedding by Joseph Prentice, September 2018        !
    !==================================================================!

    use basis, only: basis_copy_function_to_box, &
         basis_location_func_wrt_cell
    use cell_grid, only: GRID_INFO, cell_grid_deposit_box
    use comms, only: pub_my_proc_id
    use constants, only: stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc, data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: CELL_INFO

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    ! agrecocmplx
    type(FUNCTIONS), intent(in) :: mo_coeff              ! NGWF coeffs of orb.
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid         ! NGWFs on this proc
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FFTBOX_DATA), intent(inout) :: orbital_std

    ! Local Variables
    integer :: ierr
    integer :: ingwf, loc_ingwf
    integer :: ngwf_cell_start1, ngwf_cell_start2, ngwf_cell_start3
    logical :: i_have_box
    integer :: n1, n2, n3, nsub, prev_ngwfs
    type(FFTBOX_DATA) :: ngwf_box_std
    type(FFTBOX_DATA) ::  buffer_std
    ! agrecocmplx
    logical :: loc_cmplx

    if (pub_debug_on_root) &
         write(stdout,'(a)') 'DEBUG: Entering eigenstates_calculate_mo.'

    !jmecmplx: compatibility with complext NGWFs still needs to be checked
    !call utils_assert(.not. ngwfs_on_grid%iscmplx, 'Error in&
    !     & eigenstates_calculate_mo: not ready for complex NGWFs yet.')

    ! agrecocmplx: mo_coeff and ngwfs_on_grid either both real or both complex
    loc_cmplx = .false.
    if (ngwfs_on_grid%iscmplx) loc_cmplx=.true.

    ! ndmh: initialisations
    n1 = ngwf_basis%maxtight_pts1
    n2 = ngwf_basis%maxtight_pts2
    n3 = ngwf_basis%maxtight_pts3

    ! ndmh: allocate storage for buffer and NGWF tightbox
    call data_fftbox_alloc(buffer_std, n1, n2, grid%max_slabs12, iscmplx=loc_cmplx)
    call data_fftbox_alloc(ngwf_box_std, n1, n2, n3, iscmplx=loc_cmplx)
    call data_set_to_zero(ngwf_box_std)
    call data_set_to_zero(buffer_std)

    ! ndmh: Loop over NGWFs on each proc
    do loc_ingwf=1,ngwf_basis%max_on_proc
       ingwf = loc_ingwf + ngwf_basis%first_on_proc(pub_my_proc_id) - 1

       if (loc_ingwf<=ngwf_basis%num_on_proc(pub_my_proc_id)) then

          ! ndmh: copy NGWF to tightbox
          call basis_copy_function_to_box(ngwf_box_std, &
               1,1,1,ngwf_basis%tight_boxes(loc_ingwf), &
               ngwfs_on_grid,ngwf_basis%spheres(loc_ingwf),&
               cell,fftbox)

          ! agrecocmplx
          if (loc_cmplx) then
             ngwf_box_std%z(:,:,:) = ngwf_box_std%z * mo_coeff%z(ingwf)
          else
             ngwf_box_std%d(:,:,:) = ngwf_box_std%d * mo_coeff%d(ingwf)
          end if

          ! ndmh: find start of tightbox of loc_ingwf in simulation cell
          call basis_location_func_wrt_cell( &
               ngwf_cell_start1, ngwf_cell_start2, ngwf_cell_start3, &
               ngwf_basis%tight_boxes(loc_ingwf),cell)

          i_have_box = .true.
       else
          ngwf_cell_start1 = -1234
          ngwf_cell_start2 = -1234
          ngwf_cell_start3 = -1234
          i_have_box = .false.
       end if

       ! ndmh: deposit NGWF scaled by orbital coefficient to standard grid
       !jmecmplx
       ! jme: select whether complex or real case
       ! jcap: cell_grid_deposit_box automatically sums up
       ! orbital_std as we loop over regions
       if (loc_cmplx) then
          call cell_grid_deposit_box(orbital_std%z,ngwf_box_std%z,buffer_std%z, &
               grid,n1,n2,n3,n1,n2,ngwf_cell_start1,ngwf_cell_start2, &
               ngwf_cell_start3,i_have_box,.false.)
       else
          call cell_grid_deposit_box(orbital_std%d,ngwf_box_std%d,buffer_std%d, &
               grid,n1,n2,n3,n1,n2,ngwf_cell_start1,ngwf_cell_start2, &
               ngwf_cell_start3,i_have_box,.false.)
       end if
    end do

    ! ndmh: deallocate universal tightbox array
    call data_fftbox_dealloc(ngwf_box_std)
    call data_fftbox_dealloc(buffer_std)

    if (pub_debug_on_root) &
         write(stdout,'(a)') 'DEBUG: Leaving eigenstates_calculate_mo.'

  end subroutine eigenstates_calculate_mo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_plot_squared_mo( mo_kern, &
       overlap, ngwfs_on_grid, ngwf_basis, mdl, &
       file_header, scalar_name)

    !==================================================================!
    ! This subroutine prepares and outputs a plotfile for the square   !
    ! of a selected canonical molecular orbital.                       !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 10/05/2006.                  !
    ! Modified by Nicholas Hine on 04/03/2010 to take in mo_kern       !
    ! already formed, rather than create it from the orbitals.         !
    ! Trivially modified by Jacek Dziedzic on 14/05/2010 to delegate   !
    ! the unit conversion to visual_scalarfield.                       !
    !==================================================================!

    use datatypes, only: FUNCTIONS
    use comms, only: pub_on_root
    use constants, only: DP, ANGSTROM, stdout
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_trace_on_grid
    use model_type, only: MODEL
    use rundat, only: pub_cube_format, pub_dx_format, pub_grd_format,&
         pub_num_spins
    use sparse, only: SPAM3
    use utils, only: utils_alloc_check, utils_dealloc_check
    use visual, only: visual_scalarfield

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(SPAM3), intent(in) :: mo_kern(1) ! denskern of single mo as SPAM3
    type(SPAM3), intent(in)    :: overlap    ! overlap as SPAM3
    type(MODEL), intent(in) :: mdl
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid ! ngwfs on this proc
    character(len =*), intent(in) :: file_header ! part of header line of output file
    character(len =*), intent(in) :: scalar_name ! part of file-name to output (scalafield name)

    ! Local Variables
    real(kind=DP), allocatable, dimension(:,:,:,:) :: mo_density_fine
    real(kind=DP) :: norm
    integer :: ierr        ! memory allocation error flag
    integer :: nspins

    if (pub_debug_on_root) &
         write(stdout,'(a)') 'DEBUG: Entering eigenstates_plot_squared_mo.'

    ! cks: allocate memory for orbital charge density
    allocate(mo_density_fine(mdl%fine_grid%ld1, &
         mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, 1),stat=ierr)
    call utils_alloc_check('eigenstates_plot_squared_mo','mo_density_fine',ierr)

    ! cks: generate MO "density"
    ! pdh: nasty hack
    nspins = pub_num_spins
    pub_num_spins = 1
    call density_on_grid(mo_density_fine, mdl%fine_grid, mdl%dbl_grid, &
         mdl%cell, mdl%fftbox, mo_kern(1), overlap, ngwfs_on_grid, ngwf_basis, &
         ngwfs_on_grid, ngwf_basis)
    pub_num_spins = nspins

    ! ndmh: get sum of M.O. density
    norm = integrals_trace_on_grid(mo_density_fine(:,:,:,1),mdl%fine_grid)
    if (pub_on_root) write(stdout,'(f9.6,a)',advance='no') norm, ' | '

    ! cks: output MO density to plotfile formats
    ! cks: output orbital density in Angstrom^-3 rather than in Bohr^-3
    ! jd:  ... but leave the unit conversion to visual_scalarfield
    call visual_scalarfield( &
         mo_density_fine(:,:,:,1), mdl%fine_grid, mdl%cell, file_header, &
         scalar_name, mdl%elements, ANGSTROM**3)  ! input

    ! ndmh: print end of line if no message was written
    if ((.not.pub_grd_format).and.(.not.pub_cube_format).and. &
         (.not.pub_dx_format).and.pub_on_root) write(stdout,*)

    ! cks: deallocate memory for orbital charge density
    deallocate(mo_density_fine,stat=ierr)
    call utils_dealloc_check('eigenstates_plot_squared_mo', &
         'mo_density_fine',ierr)

    if (pub_debug_on_root) &
         write(stdout,'(a)') 'DEBUG: Leaving eigenstates_plot_squared_mo.'

  end subroutine eigenstates_plot_squared_mo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_dos_gp(eigen_en, s_occ, num, ham_type,&
       name,weight)

    !==================================================================!
    ! This subroutine prepares and outputs in simple txt format a      !
    ! density of states (DOS) plot which has been generated with the   !
    ! Gaussian smearing method using gamma point only energies. The    !
    ! DOS is calculated as a histogram of the sum of all Gaussian-     !
    ! broadened energies. Energies in the output are in eV and so is   !
    ! the half-width of the smearing Gaussians which is set by the     !
    ! input parameter dos_smear.                                       !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 10/05/2006.                  !
    ! Moved to eigenstates_mod by Nicholas Hine in November 2011.      !
    ! Modified for embedding by Joseph Prentice, September 2018        !
    !==================================================================!

    use comms, only: comms_barrier, pub_on_root
    use constants, only: DP, UP, DN, stdout, HARTREE_IN_EVS, PI
    use rundat, only: pub_rootname, pub_dos_smear, pub_cond_calculate,&
         pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_open_unit_check, utils_close_unit_check, utils_abort, &
         utils_banner

    implicit none

    real(kind=DP), intent(in) :: eigen_en(:,:) ! molecular orbital energies
    real(kind=dp), intent(in) :: s_occ(:)      ! total occupancy of spin channel
    integer, intent(in) :: num                 ! number of NGWFs
    character(len=*), intent(in) :: ham_type
    character(len=*), optional, intent(in) :: name
    real(kind=dp), optional, intent(in) :: weight

    ! Local Variables
    real(kind=DP) :: alpha     ! Gaussian exponent
    real(kind=DP) :: delta_e   ! energy increment
    real(kind=DP) :: e_point   ! energy point
    real(kind=DP) :: histo_val ! histogram value
    real(kind=DP) :: dist_sq   ! squared distance
    real(kind=DP) :: gnorm     ! Gaussian normalisation factor
    real(kind=DP), allocatable, dimension(:,:) :: energies ! buffer for energies in eV
    real(kind=DP), parameter :: en_offset =3.0_DP ! left and right offset for energy scale
    integer, parameter :: histonum=2001 ! number of histogram points
    integer :: row  ! row counting index
    integer :: is   ! spin loop counter
    integer :: output_unit ! fortran output unit number
    integer :: ierr ! memory allocation error flag
    integer :: io_status ! file access error flag
    integer :: orbital_index ! orbital counting index
    character(len=256) :: output_file  ! output file name buffer
    character(len=256) :: rootname

    real(kind=dp) :: loc_weight


    if (pub_debug_on_root) &
         write(stdout,'(a)') 'DEBUG: Entering eigenstates_dos_gp.'

    if (pub_on_root) then

       loc_weight=1.0_dp
       if(present(weight)) then
          loc_weight = weight
       end if


       allocate(energies(num,pub_num_spins),stat=ierr)
       call utils_alloc_check('eigenstates_dos_gp','energies',ierr)

       ! cks: store eigen energies in buffer and convert to eV
       energies = eigen_en*HARTREE_IN_EVS

       ! cks: smearing Gaussian exponent
       ! cks: pub_dos_smear is the half-width of the Gaussian
       alpha =log(2.0_DP) /((pub_dos_smear*HARTREE_IN_EVS)**2)

       ! cks: normalisation factor for Gaussians
       gnorm = sqrt(alpha/PI)*loc_weight

       ! cks: generate DOS for each spin
       spin_loop: do is=1,pub_num_spins

          if (pub_on_root) write(stdout,'(a)')' '

          if(present(name)) then
             rootname=trim(adjustl(name))
          else
             rootname=trim(adjustl(pub_rootname))
          end if

          write(output_file,*)trim(adjustl(rootname))//'_'
          if (pub_cond_calculate) &
               write(output_file,*)trim(output_file)//trim(ham_type)//'_'
          if ((pub_num_spins==2).and.(is == UP)) &
               write(output_file,*) trim(output_file)//'UP'//'_'
          if ((pub_num_spins==2).and.(is == DN)) &
               write(output_file,*) trim(output_file)//'DN'//'_'
          write(output_file,*) trim(output_file)//'DOS.txt'

          if (pub_num_spins == 1) then
             write(stdout,'(a)') &
                  utils_banner('=', 'Density of States (DOS) calculation')
          elseif (is == UP) then
             write(stdout,'(a)') &
                  utils_banner('=', 'Density of States (DOS) calculation for UP spin')
          elseif (is == DN) then
             write(stdout,'(a)') &
                  utils_banner('=', 'Density of States (DOS) calculation for DOWN spin')
          endif

          output_file = adjustl(output_file)

          ! cks: print output warning
          write(stdout,'(3a)',advance ='no') &
               'Writing "', trim(output_file),'" ...'

          ! cks: get a unit number that is free
          output_unit = utils_unit()

          open(unit=output_unit, form="formatted" ,file=trim(output_file), &
               action="write",iostat=io_status)
          call utils_open_unit_check('eigenstates_dos_gp','output_file', &
               io_status)

          ! cks: write first line
          write(output_unit,'(a)',err =100)'#  Energy (eV) |  DOS (states/eV)'

          delta_e =(2.0_DP*en_offset + maxval(energies(num,:)) &
               - minval(energies(1,:)))/real(histonum-1, kind=DP)

          e_point = minval(energies(1,:)) - en_offset

          ! cks: Loop over DOS histogram points
          do row=1, histonum

             histo_val = 0.0_DP
             ! cks: accumulate value at current histogram point
             do orbital_index=1,num
                dist_sq =(e_point -energies(orbital_index, is))**2
                ! CW: protect exponential function
                if (abs(alpha*dist_sq)<500.0_DP) then
                   histo_val = histo_val +gnorm*exp(-alpha *dist_sq)
                endif
             end do

             ! cks: write to file current histo-point
             write(output_unit,'(f14.8,f16.10)',err=100) e_point, histo_val

             ! cks: update next histogram coordinate
             e_point = e_point +delta_e
          end do

          close(unit=output_unit,iostat=io_status)
          call utils_close_unit_check('eigenstates_dos_gp','output_unit', &
               io_status)

          ! cks: notify of end of output
          write(stdout,*)' done'

          write(stdout,'(a)') repeat('=',80)

       enddo spin_loop



       deallocate(energies,stat=ierr)
       call utils_dealloc_check('eigenstates_dos_gp','energies',ierr)

    endif
    call comms_barrier

    if (pub_debug_on_root) &
         write(stdout,'(a)') 'DEBUG: Leaving eigenstates_dos_gp.'

    return


100 call utils_abort('Problem writing to file in eigenstates_dos_gp.')

  end subroutine eigenstates_dos_gp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_ldos_gp(eigen_en, eigenvecs_dens, s_eigenvecs_dens, &
       is, mdl, ngwf_basis, hubbard_switch, ham_type)

    !===================================================================!
    ! This subroutine prepares and outputs in simple txt format a       !
    ! Local Density of States (LDOS) plot which has been generated      !
    ! with the Gaussian smearing method, using gamma point only         !
    ! energies. The groups over which to sum the LDOS are defined by    !
    ! pub_ldos_groups, which is set through the species_ldos_groups     !
    ! block in the input file. LDOS is calculated as a histogram of the !
    ! sum of all Gaussian-broadened energies multiplied by a coefficient!
    ! expressing the contribution of the NGWFs on each atom in the LDOS !
    ! group to that orbital. Energies in the output are in eV and so is !
    ! the half-width of the smearing Gaussians which is set by the      !
    ! input parameter ldos_smear.                                       !
    !-------------------------------------------------------------------!
    ! Written by Nicholas Hine on 15th February 2010, based on parts of !
    ! eigenstates_dos_gp which was written by Chris-Kriton Skylaris on  !
    ! 10/05/2006.                                                       !
    ! Modified to include Hubbard subspaces by David O'Regan, 10/2/2011.!
    ! Moved to eigenstates_mod by Nicholas Hine in November 2011.       !
    ! Modified for embedding structures by Robert Charlton, 26.06/2018. !
    !===================================================================!

    use bibliography, only: bibliography_cite
    use comms, only: comms_barrier, comms_reduce, &
         pub_my_proc_id, pub_on_root
    use constants, only: DP, UP, DN, stdout, HARTREE_IN_EVS, PI
    ! agrecocmplx
    use datatypes, only: FUNCTIONS, data_functions_alloc, data_functions_dealloc
    use dense, only: DEM, dense_get_col
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_rootname, pub_dos_smear, pub_ldos_ngroups, &
         pub_ldos_group_nsp, pub_ldos_groups, pub_cond_calculate, pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_open_unit_check, utils_close_unit_check, utils_abort, &
         utils_banner

    implicit none

    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)   ! Function basis for NGWFs
    real(kind=DP), intent(in) :: eigen_en(sum(ngwf_basis(:)%num))
    type(DEM), intent(in) :: eigenvecs_dens
    type(DEM), intent(in) :: s_eigenvecs_dens
    type(MODEL), intent(in), target :: mdl
    integer, intent(in) :: is                    ! Current spin
    logical, intent(in) :: hubbard_switch        ! ddor: DFT+U projected LDOS
    character(len=*), intent(in) :: ham_type

    ! Local Variables
    real(kind=DP) :: alpha     ! Gaussian exponent
    real(kind=DP) :: delta_e   ! energy increment
    real(kind=DP) :: e_point   ! energy point
    real(kind=DP) :: gnorm     ! Gaussian normalisation factor
    real(kind=DP), allocatable, dimension(:) :: orb_coeff
    real(kind=DP), allocatable, dimension(:) :: energies ! buffer for energies in eV
    real(kind=DP), allocatable, dimension(:,:) :: histo_val ! histogram values
    real(kind=DP), parameter :: en_offset =3.0_DP ! left and right offsets
    ! agrecocmplx: use FUNCTIONS type to switch between real and complex case
    type(FUNCTIONS) :: eig_col, s_eig_col
    integer, parameter :: histonum=2001 ! number of histogram points
    integer :: row  ! row counting index
    integer :: output_unit ! fortran output unit number
    integer :: ierr ! memory allocation error flag
    integer :: igroup
    integer :: iat, loc_iat
    integer :: isp
    integer :: ingwf, isub, last_ngwf
    integer :: io_status ! file access error flag
    integer :: orbital_index ! orbital counting index
    logical :: iat_in_group
    character(len=256) :: output_file  ! output file name buffer
    character(len=4) :: iat_orig_id
    ! agrecocmplx
    logical :: loc_cmplx
    real(kind=DP) :: temp_real
    integer :: num ! rc2013: total number of NGWFs in all bases
    type(PARAL_INFO), pointer :: par ! rc2013: regional parallel strategies

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering eigenstates_ldos_gp'

    call bibliography_cite('LDOS')

    ! agrecocmplx
    loc_cmplx = eigenvecs_dens%iscmplx

    ! rc2013: get total number of NGWFs
    num = sum(ngwf_basis(:)%num)

    ! ndmh: Allocate temporary arrays
    allocate(energies(num),stat=ierr)
    call utils_alloc_check('eigenstates_ldos_gp','energies',ierr)
    allocate(histo_val(histonum,0:pub_ldos_ngroups+1),stat=ierr)
    call utils_alloc_check('eigenstates_ldos_gp','histo_val',ierr)
    ! agrecocmplx: allocate using appropriate routines
    call data_functions_alloc(eig_col,num,iscmplx=loc_cmplx)
    call data_functions_alloc(s_eig_col,num,iscmplx=loc_cmplx)
    allocate(orb_coeff(pub_ldos_ngroups+1),stat=ierr)
    call utils_alloc_check('eigenstates_ldos_gp','orb_coeff',ierr)

    ! cks: store eigen energies in buffer and convert to eV
    energies = eigen_en
    energies = energies*HARTREE_IN_EVS
    ! qoh: Initialise to avoid compiler warning
    output_unit = -1

    if (pub_on_root) then

       write(stdout,'(/a)') &
            utils_banner('=', 'Local Density of States (LDOS) calculation')

       ! ndmh: find name for output files
       ! ddor: Include possibility of DFT+U projected LDOS
       ! ndmh: Include possibility of various COND task ham_types
       write(output_file,*) trim(pub_rootname)
       if (pub_cond_calculate) write(output_file,*) &
            trim(output_file)//'_'//trim(ham_type)
       if (hubbard_switch) write(output_file,*) &
            trim(output_file)//'_'//'HUB'
       if ((pub_num_spins==2).and.(is == UP)) &
            write(output_file,*) trim(output_file)//'_'//'UP'
       if ((pub_num_spins==2).and.(is == DN)) &
            write(output_file,*) trim(output_file)//'_'//'DN'
       write(output_file,*) trim(output_file)//'_LDOS.txt'

       output_file = adjustl(output_file)

       ! cks: print output warning
       write(stdout,'(3a)',advance ='no') 'Writing "', trim(output_file),'" ...'

       ! cks: get a unit number that is free
       output_unit = utils_unit()

       open(unit=output_unit, form="formatted" ,file=trim(output_file), &
            action="write",iostat=io_status)
       call utils_open_unit_check('eigenstates_ldos_gp','output_file', &
            io_status)

       ! ndmh: write first line
       write(output_unit,'(a)',advance='no',err =100)'#  Energy (eV) |'
       do igroup=1,pub_ldos_ngroups
          write(output_unit,'(a)',advance='no',err=100) '  LDOS ('
          do isp=1,pub_ldos_group_nsp(igroup)
             write(output_unit,'(a)',advance='no',err=100) &
                  pub_ldos_groups(isp,igroup)
             if (isp<pub_ldos_group_nsp(igroup)) &
                  write(output_unit,'(a)',advance='no',err=100) ', '
          end do
          write(output_unit,'(a)',advance='no',err=100) ') |'
       end do
       write(output_unit,'(a)',err=100) '  Total DOS    |'

    end if

    ! ndmh: initialisations
    histo_val(:,:) = 0.0_DP
    delta_e = (2.0_DP*en_offset + energies(num) - energies(1)) / &
         real(histonum-1, kind=DP)

    ! ndmh: smearing Gaussian exponent
    ! ndmh: pub_ldos_smear is the half-width in eVs
    alpha = log(2.0_DP) / ((pub_dos_smear*HARTREE_IN_EVS)**2)
    gnorm = sqrt(alpha/PI)

    ! ndmh: first calculate the total DOS, as last column of histogram
    ! ddor: However this is not needed for DFT+U
    if (.not. hubbard_switch) then
       do orbital_index=1,num

          e_point = energies(1) - en_offset
          do row=1,histonum
             ! CW: protect exponential function
             if (abs(alpha * (e_point - energies(orbital_index))**2) <  500.0_DP) then
                histo_val(row,pub_ldos_ngroups+1) = &
                     histo_val(row,pub_ldos_ngroups+1) &
                     + gnorm*exp(-alpha*(e_point - energies(orbital_index))**2)
             endif
             e_point = e_point + delta_e
          end do

       end do
    endif

    ! ndmh: accumulate contributions to LDOS for each orbital
    do orbital_index=1,num

       ! ndmh: Loop over LDOS histogram points to set up exp(-alpha*(e-ei)^2)
       ! ndmh: array for each point, where ei is the current eigenvalue
       e_point = energies(1) - en_offset
       do row=1,histonum
          ! CW: protect exponential function
          if (abs(alpha * (e_point - energies(orbital_index))**2) <  500.0_DP) then
             histo_val(row,0) = &
                  gnorm*exp(-alpha*(e_point - energies(orbital_index))**2)
          else
             histo_val(row,0) = 0.0_DP
          endif
          e_point = e_point + delta_e
       end do

       ! agrecocmplx
       if (loc_cmplx) then
          call dense_get_col(eig_col%z,eigenvecs_dens,orbital_index)
          call dense_get_col(s_eig_col%z,s_eigenvecs_dens,orbital_index)
       else
          call dense_get_col(eig_col%d,eigenvecs_dens,orbital_index)
          call dense_get_col(s_eig_col%d,s_eigenvecs_dens,orbital_index)
       end if

       orb_coeff(:) = 0.0_DP
       last_ngwf = 0
       ingwf = 1

       ! rc2013: loop over regions of atoms
       do isub=1,mdl%nsub
          par=>mdl%regions(isub)%par

          ! ndmh: loop over atoms to add up contributions to the local DOS
          do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
             iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
             iat_orig_id = &
                  mdl%regions(isub)%elements(par%orig_atom(iat))%species_id

             ! ndmh: go to next atom if this atom does not occur in any LDOS groups
             iat_in_group = .false.
             do igroup=1,pub_ldos_ngroups
                if (any(pub_ldos_groups(1:pub_ldos_group_nsp(igroup),igroup) &
                     ==iat_orig_id)) iat_in_group=.true.
             end do
             ! ddor: Full summation for DFT+U total LDOS
             if ((.not.iat_in_group) .and. (.not. hubbard_switch)) cycle

             ! ndmh: loop over NGWFs on this atom
             ! rc2013: start from the last NGWF in the preceding basis
             do ingwf = last_ngwf + ngwf_basis(isub)%first_on_atom(iat), &
                  last_ngwf + ngwf_basis(isub)%first_on_atom(iat) +  &
                  ngwf_basis(isub)%num_on_atom(iat)-1

                ! ndmh: evaluate orbital coefficient for NGWF a
                ! ndmh: given by M(dagger)_an * sum_b(S_ab M_nb)

                ! ndmh: add contribution of this eigenvalue to any of the LDOS
                ! ndmh: groups in which this atom is listed
                ! ddor: iat_in_group check for DFT+U
                ! agrecocmplx: the contribution should be real anyway, use temp_real
                ! to store the appropriate quantity for real/complex case
                if (loc_cmplx) then
                   temp_real = real(conjg(eig_col%z(ingwf)) * s_eig_col%z(ingwf),kind=DP)
                else
                   temp_real = eig_col%d(ingwf) * s_eig_col%d(ingwf)
                end if

                if (iat_in_group) then
                   do igroup=1,pub_ldos_ngroups

                      if (any(pub_ldos_groups(1:pub_ldos_group_nsp(igroup),igroup) &
                           ==iat_orig_id)) orb_coeff(igroup) = &
                           orb_coeff(igroup) + temp_real

                   end do  ! igroup
                endif

                ! ddor: Total DFT+U LDOS
                if (hubbard_switch) then
                   orb_coeff(pub_ldos_ngroups+1) = orb_coeff(pub_ldos_ngroups+1) + &
                        temp_real
                endif

             end do  ! ingwf

          end do  ! iat
          ! rc2013: store the number of NGWFs we've included so far so that we &
          ! get the contributions from each region right
          last_ngwf = ingwf - 1

       end do  ! isub

       ! ndmh: Add up contributions to orb_coeff for each group from all procs
       call comms_reduce('SUM',orb_coeff)

       ! ndmh: Add up contributions for each group from this orbital from all NGWFs
       do igroup=1,pub_ldos_ngroups
          histo_val(:,igroup) = histo_val(:,igroup) + &
               orb_coeff(igroup) * histo_val(:,0)
       end do

       ! ddor: Total DFT+U LDOS
       if (hubbard_switch) then
          histo_val(:,pub_ldos_ngroups+1) = histo_val(:,pub_ldos_ngroups+1) + &
               orb_coeff(pub_ldos_ngroups+1) * histo_val(:,0)
       endif

    end do  ! orbital_index

    if (pub_on_root) then

       ! ndmh: loop over histogram points writing output file
       e_point = energies(1) - en_offset
       do row=1,histonum

          ! ndmh: write to file current histo-point for each group, and total
          write(output_unit,'(f14.8,f16.10)',advance='no',err=100) e_point
          do igroup=1,pub_ldos_ngroups
             write(output_unit,'(f16.10)',advance='no',err=100) &
                  histo_val(row,igroup)
          end do
          write(output_unit,'(f16.10)') histo_val(row,pub_ldos_ngroups+1)

          e_point = e_point + delta_e

       end do

       ! ndmh: close output file
       close(unit=output_unit,iostat=io_status)
       call utils_close_unit_check('eigenstates_ldos_gp','output_unit', &
            io_status)

       ! ndmh: finish writing message to stdout
       write(stdout,*)' done'
       write(stdout,'(a)') repeat('=',80)

    endif

    ! ndmh: deallocate temporary arrays
    deallocate(orb_coeff,stat=ierr)
    call utils_dealloc_check('eigenstates_ldos_gp','orb_coeff',ierr)
    ! agrecocmplx: deallocate using appropriate arrays
    call data_functions_dealloc(s_eig_col)
    call data_functions_dealloc(eig_col)
    deallocate(histo_val,stat=ierr)
    call utils_dealloc_check('eigenstates_ldos_gp','histo_val',ierr)
    deallocate(energies,stat=ierr)
    call utils_dealloc_check('eigenstates_ldos_gp','energies',ierr)

    call comms_barrier

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving eigenstates_ldos_gp'

    return


100 call utils_abort('Problem writing to file in eigenstates_ldos_gp.')

  end subroutine eigenstates_ldos_gp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_write_bands(ngwf_basis,cell,eigen_en, s_occ, ham_type)

    !============================================================!
    ! This subroutine writes to a file the orbitals in CASTEP    !
    ! .bands format.                                             !
    !------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 18/02/2007.            !
    ! Bug in Fermi level choice for fully-spin polarised systems !
    ! fixed by Nicholas Hine on 08/06/2010. Also cleaned up the  !
    ! code somewhat.                                             !
    ! Moved to eigenstates_mod by Nicholas Hine in November 2011.!
    ! Modified for embedding by Joseph Prentice, September 2018  !
    !============================================================!

    use comms, only: pub_on_root, comms_barrier
    use constants, only: DP, max_spins, stdout
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_rootname, pub_num_spins, pub_edft
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort, utils_unit, utils_open_unit_check, &
         utils_close_unit_check
    use ensemble_dft, only: edft_zerok_fermi

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(CELL_INFO), intent(in) :: cell
    real(kind=DP), intent(in) :: eigen_en(:,:) ! band eigen-energies
    real(kind=dp), intent(in) :: s_occ(:)      ! spin channel occupancy
    ! lr408: Optional argument to set different filenames for different diagonalisations in
    ! lr408: a conduction calculation
    character(*), optional, intent(in) :: ham_type

    ! Local Variables
    integer :: output_unit  ! unit name of output file
    integer :: ierr         ! unit open/close error flag
    integer :: is           ! spin counter
    integer :: neigen       ! eigenvalues counter
    character(len=256) :: filename  ! output file name buffer
    real(kind=DP)     :: efermi(max_spins) ! "Fermi level" for each spin
    integer :: num

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &eigenstates_write_bands'

    ! jcap: set total number of NGWFs
    num=sum(ngwf_basis(:)%num)

    if (pub_on_root) then

       ! cks: Initialise "Fermi levels" to HOMO energies
       ! ndmh: modified to check for invalid n_occ
       if (pub_edft) then
          ! kkbd: Use the 0K Fermi level if appropriate
          call edft_zerok_fermi(efermi,s_occ,eigen_en)
       else
          do is=1,pub_num_spins
             if ((s_occ(is) > 0).and.(s_occ(is) <= num)) then
                efermi(is) = eigen_en(nint(s_occ(is)),is)
             else if (nint(s_occ(is)) == 0) then
                efermi(is) = eigen_en(1,is)
             else
                call utils_abort('Error in eigenstates_write_bands: n_occ(is) < 0 &
                     &or n_occ(is) > N_NGWFs')
             end if
          end do
       end if

       ! cks: name of output file
       if (present(ham_type)) then
          if (ham_type=='valence') then
             write(filename,*)trim(pub_rootname)//'.val_bands'
          else if (ham_type=='proj') then
             write(filename,*)trim(pub_rootname)//'.proj_bands'
          else if (ham_type=='cond') then
             write(filename,*)trim(pub_rootname)//'.cond_bands'
          else if (ham_type=='joint') then
             write(filename,*)trim(pub_rootname)//'.joint_bands'
          else if (ham_type=='aux') then !jhl52: aux
             write(filename,*)trim(pub_rootname)//'.aux_bands'
          end if
       else
          write(filename,*)trim(pub_rootname)//'.bands'
       end if

       filename = adjustl(filename)

       write(stdout,*)''
       write(stdout,'(3a)', advance='no')&
            'Writing bands to file "',trim(adjustl(filename)),'" ...'

       ! cks: find free unit to open
       output_unit =utils_unit()

       ! cks: open output unit
       open(unit=output_unit, file=trim(filename), form='formatted', &
            action='write', iostat=ierr, status='replace')
       call utils_open_unit_check('eigenstates_write_bands','output_unit',ierr)

       ! Write the number of k-points, spin components and eigenvalues to file
       write(output_unit,'(a,i4)') 'Number of k-points        ', 1
       write(output_unit,'(a,i4)') 'Number of spin components ', &
            pub_num_spins

       if (pub_num_spins == 1) then
          write(output_unit,'(a,g11.4)') 'Number of electrons       ', &
               2.0*sum(s_occ)
       else
          write(output_unit,'(a,2g11.4)') 'Number of electrons       ', &
               s_occ(1), s_occ(2)
       end if

       if (pub_num_spins == 1) then
          write(output_unit,'(a,i6)') 'Number of eigenvalues     ', &
               num
       else
          write(output_unit,'(a,2i6)') 'Number of eigenvalues     ', &
               num, num
       end if

       if( pub_num_spins == 1 ) then
          write(output_unit,'(a,f12.6)') '0K Fermi energy (in atomic units) ', &
               efermi(1)
       else
          write(output_unit,'(a,2f12.6)') '0K Fermi energy (in atomic units) ', &
               efermi(1:pub_num_spins)
       end if

       write(output_unit,'(a)') 'Unit cell vectors'
       write(output_unit,'(3f12.6)') cell%a1%x, cell%a1%y, cell%a1%z
       write(output_unit,'(3f12.6)') cell%a2%x, cell%a2%y, cell%a2%z
       write(output_unit,'(3f12.6)') cell%a3%x, cell%a3%y, cell%a3%z

       ! cks: write the actual band energies for each spin component
       write(output_unit,'(a,i4,4f12.8)') 'K-point ', &
            1, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP
       do is=1,pub_num_spins
          write(output_unit,'(a,i1)') 'Spin component ',is
          do neigen=1,num
             write(output_unit,'(f14.8)') eigen_en(neigen,is)
          end do
       end do

       ! cks: close the output unit (file)
       close(unit=output_unit,iostat=ierr)
       call utils_close_unit_check('eigenstates_write_bands','output_unit',ierr)
       write(stdout,'(a)') ' done'

    endif
    call comms_barrier

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &eigenstates_write_bands'

    return

  end subroutine eigenstates_write_bands

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eigenstates_init_pdos(mdl,ngwf_basis,nmax,par)
    !==========================================================================!
    ! This routine sets up the function basis that we will need for a PDOS     !
    ! calculation
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! mdl (inout)              : To modify the pdos nfunctions in %elements    !
    ! ngwf_basis (input)       : The function basis for the NGWFs              !
    ! nmax (input)             : Max size of AM-resolved basis set             !
    !--------------------------------------------------------------------------!
    ! Written by Jolyon Aarons. Modified most recently in Jan 2019             !
    ! Modified for embedding by Joseph Prentice, September 2018                !
    ! Time travel is possible...                                               !
    !==========================================================================!

    use comms, only: comms_reduce
    use constants, only: stdout
    use function_basis, only: FUNC_BASIS, function_basis_allocate, &
         function_basis_distribute, function_basis_init_spheres, &
         function_basis_init_tight_boxes, function_basis_gath_all_tbs
    use geometry, only: POINT
    use ion, only: ELEMENT
    use model_type, only: MODEL

    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_pdos_max_l, pub_pdos_reduce_sws, pub_pdos_max_n, pub_pdos_stride_n,&
         pub_pdos_construct_basis, pub_pdos_pseudoatomic

    use sparse, only: BLKS_PDOS, sparse_init_blocking_scheme
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    type(MODEL), intent(inout) :: mdl !inout because I needed to modify %elements
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    integer, optional, intent(in) ::nmax
    type(PARAL_INFO), intent(in), pointer :: par

    integer, allocatable       :: num_orbitals_on_atom(:)
    type(element), allocatable :: lcao_elements(:)
    integer :: num_orbitals
    integer :: atom, ll
    integer :: ierr, ireg
    integer :: max_n, num_n

    type(POINT) :: a1, a2, a3
    real(kind=DP) :: d1, d2, d3
    integer       :: n_pt1, n_pt2, n_pt3
    integer       :: n1,n2,n3

    if (pub_debug_on_root) &
         write(stdout,'(a)') 'DEBUG: Entering eigenstates_init_pdos.'

    ! jcap: get region from par
    ireg=par%par_index

    ! JA: work out how many angular momentum projected functions we will have
    if(pub_pdos_reduce_sws) then
       if(pub_pdos_construct_basis) then
          num_orbitals = par%num_ngwfs*((pub_pdos_max_l+1)**2)
       else
          num_orbitals=0
          do atom=1,par%nat
             do ll=0,pub_pdos_max_l
                num_orbitals = num_orbitals + (2*ll + 1)
             end do
          end do
       end if
       num_n=1

    else if(pub_pdos_pseudoatomic) then
       num_orbitals = par%num_ngwfs
    else
       ! check optional arguments
       if (present (nmax)) then
          max_n = nmax
       elseif(pub_pdos_max_n>0) then
          max_n=pub_pdos_max_n
       else
          ! default is number of psincs across largest NGWF
          max_n = ceiling(maxval(ngwf_basis%spheres(:)%radius) / &
               min(mdl%cell%d1,mdl%cell%d2,mdl%cell%d3))
          call comms_reduce('MAX',max_n)

       endif

       num_n=((max_n-1)/pub_pdos_stride_n)+1

       num_orbitals=0
       do atom=1,par%nat
          do ll=0,pub_pdos_max_l
             num_orbitals = num_orbitals + (2*ll + 1)
          end do
       end do



       num_orbitals=num_orbitals*num_n
    end if

    ! jcap: allocate lcao_basis, if it hasn't been allocated yet
    if (.not.allocated(lcao_basis)) then
       allocate(lcao_basis(mdl%nsub),stat=ierr)
       call utils_alloc_check('eigenstates','lcao_basis',ierr)
    end if

    allocate(lcao_elements(size(mdl%regions(ireg)%elements)),stat=ierr)
    call utils_alloc_check('eigenstates_init_pdos','lcao_elements',ierr)
    lcao_elements=mdl%regions(ireg)%elements

    ! Number of orbitals on atoms
    allocate(num_orbitals_on_atom(1:par%nat),stat=ierr)
    call utils_alloc_check('eigenstates_init_pdos','num_orbitals_on_atom',ierr)

    if(pub_pdos_reduce_sws) then
       if(pub_pdos_construct_basis) then
          num_orbitals_on_atom = mdl%regions(ireg)%elements(:)%nfunctions*((pub_pdos_max_l+1)**2)
       else
          num_orbitals_on_atom=(pub_pdos_max_l+1)**2
       end if
    else if(pub_pdos_pseudoatomic) then
       num_orbitals_on_atom = mdl%regions(ireg)%elements(:)%nfunctions
    else
       num_orbitals_on_atom=((pub_pdos_max_l+1)**2)*num_n
    end if

    ! set up functions

    mdl%elements(:)%nfunctions_pdos = num_orbitals_on_atom
    lcao_elements(:)%nfunctions = num_orbitals_on_atom

    call function_basis_allocate(lcao_basis(ireg),num_orbitals,"pdos_sws", par=par)

!    call function_basis_distribute(lcao_basis(ireg),lcao_elements,par)
    call function_basis_distribute(lcao_basis(ireg),mdl%elements,par)

    call function_basis_init_spheres(lcao_basis(ireg),mdl%cell,mdl%fftbox, &
         par%elements_on_proc, par=par)

    call function_basis_init_tight_boxes(lcao_basis(ireg),mdl%fftbox,mdl%cell, par=par)

    call function_basis_gath_all_tbs(lcao_basis(ireg), par=par)

    ! Haven't moved this to Energy and Force because that would mean passing lcao_basis around... when it is only used here.
    ! rc2013: pass par explicitly for sparse initialisation
    call sparse_init_blocking_scheme(BLKS_PDOS,num_orbitals,lcao_basis(ireg)%num_on_proc,&
         lcao_basis(ireg)%num_on_atom,lcao_basis(ireg)%first_on_proc,&
         lcao_basis(ireg)%first_on_atom,lcao_basis(ireg)%atom_of_func,&
         lcao_basis(ireg)%proc_of_func, par)

    deallocate(num_orbitals_on_atom,stat=ierr)
    call utils_dealloc_check('eigenstates_init_pdos','num_orbitals_on_atom',ierr)
    deallocate(lcao_elements,stat=ierr)
    call utils_dealloc_check('eigenstates_init_pdos','lcao_elements',ierr)

    if (pub_debug_on_root) &
         write(stdout,'(a)') 'DEBUG: Leaving eigenstates_init_pdos.'

  end subroutine eigenstates_init_pdos

  subroutine eigenstates_cleanup_pdos
    use function_basis, only: function_basis_deallocate
    use utils, only: utils_dealloc_check
    implicit none

    integer :: isub,ierr

    do isub=1,size(lcao_basis)
       call function_basis_deallocate(lcao_basis(isub))
    end do

    deallocate(lcao_basis,stat=ierr)
    call utils_dealloc_check('eigenstates','lcao_basis',ierr)

  end subroutine eigenstates_cleanup_pdos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eigenstates_pdos(ngwfs_on_grid,ngwf_basis,proj_basis,nl_projectors, &
       & eigen_en,eigen_occ,eigenvecs_dens,s_occ,rep,mdl,ham_type,nmax,output)
    !==========================================================================!
    ! This routine builds the angular momentum projected density of states     !
    ! (p-DOS) weights matrix and uses it to build a Gaussian smeared p-DOS     !
    ! histogram, which is outputted to file. Optionally, p-DOS weights can be  !
    ! written to Castep-style *seed*.pdos_weights file, suitable for use in    !
    ! OptaDOS, along with band gradients and a Castep-style 'output cell'      !
    ! file.                                                                    !
    !                                                                          !
    ! The expression used for the weights is:                                  !
    ! M_i^{\dagger\alpha}\left<\varphi_\alpha|\chi^{(A)}\right>                !
    !         \Lambda_{AB}^{-1}\left<\chi^{(B)}\varphi_\beta\right>M_i^{\beta} !
    !                                                                          !
    ! Where M_i are the eigenvectors of the Hamiltonian, \varphi are the NGWFs !
    ! \chi are the spherical waves and \Lambda is the SW-SW overlap matrix.    !
    ! The spherical waves are indexed by l,m,n quantum numbers, so to produce  !
    ! the pDOS weights matrix, they reduced over n, leaving a hilbert space of !
    ! natoms*l*m. This 'full' weights matrix is written out for OptaDOS, and   !
    ! also further reduced over m and the atoms to give l+1 vectors of length  !
    ! number of eigen-energies, which which to weight the outputted histograms !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! ngwfs_on_grid (input)    : NGWFs in ppd representation on this proc      !
    ! ngwf_basis (input)       : The function basis for the NGWFs              !
    ! elements (input)         : Element array for atoms                       !
    ! lmax (input)             : Largest angular momentum SW to include        !
    ! nmax (input)             : Size of spherical wave basis set              !
    ! output (input)           : Should we write out SW coefficients to file   !
    !--------------------------------------------------------------------------!
    ! Written by Jolyon Aarons. 23/07/2014                                     !
    ! Modified by Jolyon Aarons to do Lowdin orthogonalization and support the !
    ! pseudoatomic functions as an AM-resolved basis in 2019.
    ! Modified for embedding by Joseph Prentice, September 2018                !
    !==========================================================================!

    use augmentation, only: augmentation_overlap
    use basis, only: basis_extract_function_from_box, &
         basis_location_func_wrt_cell, basis_add_function_to_box, &
         basis_clean_function, basis_dot_function_with_box
    use comms, only: pub_on_root, pub_total_num_procs, pub_my_proc_id, &
         comms_barrier, comms_reduce, comms_allgather
    use constants, only: DP, stdout, two_pi, cmplx_0, UP, DN, NORMAL
    use datatypes, only: COEF, FUNCTIONS, FFTBOX_DATA, &
         data_fftbox_alloc, data_fftbox_dealloc, &
         data_functions_alloc, data_functions_dealloc, data_set_to_zero, &
         data_fftbox_scale, data_functions_copy
    use dense, only: DEM, dense_create, &
         dense_destroy,dense_convert, dense_get_element, dense_product, &
         dense_invert, dense_copy, dense_trace, dense_normal_eigensolve, &
         dense_scale
    use electronic, only: electronic_energy
    use electronic_init, only: electronic_init_denskern
    use ensemble_dft, only: edft_create,edft_destroy, edft_calculate, &
         edft_matrix_from_eigs
    use ensemble_dft_type, only: EDFT_MODEL
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketppd, function_ops_sum_ppd_funcs
    use geometry, only: point, operator(.DOT.), operator(*), operator(+)
    use hubbard_build, only: hubbard_model
    use ion, only: element
    use kernel, only: DKERN, kernel_create, kernel_destroy
    use model_type, only: MODEL, model_write_castep_cell_file

    use ngwf_representation, only: NGWF_REP, ngwf_rep_create, ngwf_rep_destroy, &
         NGWF_HAM, ngwf_ham_create, ngwf_ham_destroy
    use ngwfs, only: ngwfs_initialise_ngwfs

    use projectors, only: projector_set, projectors_func_ovlp_box
    use restart, only: restart_ngwfs_tightbox_output, restart_hamiltonian_write, restart_kernel_write
    use rundat, only: pub_rootname, pub_pdos_max_l, pub_any_nl_proj, pub_paw, pub_aug, pub_debug, &
         pub_pdos_optados_output, pub_print_qc, pub_num_spins, pub_pdos_reduce_sws, pub_realspace_projectors, &
         pub_pdos_max_n, pub_pdos_stride_n, pub_pdos_construct_basis, pub_num_kpoints, lnv_threshold_orig, &
         pub_firstit_lnv, pub_charge, pub_spin, pub_1K, pub_eigensolver_abstol, pub_eigensolver_orfac, &
         pub_edft, pub_pdos_lowdin, pub_pdos_lcao_optimize, pub_pdos_orth_atom_blocks, pub_spin_fac,&
         pub_pdos_output_basis, pub_pdos_output_swopt_kernham, pub_pdos_pseudoatomic, pub_output_detail
    use sparse, only: SPAM3, sparse_create,&
         sparse_hotelling_init, sparse_scale, &
         sparse_hotelling_invert, sparse_product, sparse_transpose, sparse_num_cols, &
         sparse_num_rows, sparse_destroy, sparse_copy, sparse_trace, sparse_show_matrix,&
         sparse_get_block, sparse_put_block, sparse_get_element, sparse_1norm
    use sparse_array, only: SPAM3_ARRAY
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_scale, sparse_embed_product, sparse_embed_transpose, &
         sparse_embed_num_cols, sparse_embed_num_rows, sparse_embed_destroy, &
         sparse_embed_copy, sparse_embed_trace, sparse_embed_1norm
    use spherical_wave, only: sw_bessel_zeros, sw_init, sw_recp_generate

    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_open_unit_check, utils_close_unit_check, utils_abort, utils_qc_print, &
         utils_assert, utils_banner, utils_feature_not_supported
    use visual, only: visual_ngwfs

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(FUNC_BASIS),              intent(in)    :: proj_basis(:)
    type(PROJECTOR_SET),           intent(inout) :: nl_projectors(:)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid(:)
    real(kind=DP), intent(in) :: eigen_en(:,:) ! band eigen-energies
    real(kind=DP), intent(in) :: eigen_occ(:,:) ! band eigen-occupancies
    type(DEM), intent(inout)  :: eigenvecs_dens(:)
    real(kind=dp), intent(in)        :: s_occ(:)      ! number of occupied bands
    type(NGWF_REP), intent(in) :: rep
    type(MODEL), intent(inout) :: mdl
    character(*), optional, intent(in) :: ham_type
    integer, optional, intent(in) :: nmax
    logical, optional, intent(in) :: output

    integer :: max_n ! Internal version of nmax
    integer :: num_n

    ! counters
    integer       :: atom,spin,neigen,kpoint
    integer       :: atom_ind
    integer       :: ll,l2,mm,nn,proc        ! loop counters
    integer       :: icol

    ! variables
    type(POINT)   :: fbcentre,centre,frac,disp
    real(kind=DP) :: radius
    real(kind=DP) :: spil_param,spil_param_occ,spil_elec  !spilling parameter
    integer       :: ierr
    integer       :: sx,sy,sz
    integer       :: num_orbitals
    integer       :: num_orbitals_on_proc
    integer       :: max_orbitals_on_proc
    integer       :: old_atom
    !    integer       :: old_orb_index_on_proc
    !    integer       :: orb_index_on_proc

    ! The number of p-DOS orbitals on each proc.
    integer, allocatable       :: num_orbitals_on_procs(:)
    integer, allocatable       :: num_orbitals_on_atom(:)
    !    integer, allocatable       :: cum_orbitals_on_atom(:)

    integer, allocatable       :: atom_in_species(:)
    integer, allocatable       :: count_in_species(:)
    integer :: num_species

    ! spherical bessel function zeros array
    real(kind=DP), allocatable :: qnl(:)

    ! summation arrays
    real(kind=DP), allocatable :: pdos_weights(:,:,:)
    !    real(kind=DP), allocatable :: pdos_weights_trans(:,:)

    ! output arrays
    integer, allocatable :: orbital_order(:)
    integer, allocatable :: orbital_species(:)
    integer, allocatable :: orbital_ion(:)
    integer, allocatable :: orbital_l(:)

    ! file data (when printing out pDOS weights)
    integer            :: fileunit
    character(len=256) :: filename

    character(len=8)  :: date
    character(len=80) :: file_header
    real(kind=dp)     :: file_version

    ! JA: New
    !    type(element), allocatable :: lcao_elements(:)
    !    type(func_basis) :: lcao_basis  ! Moved to global scope
    integer :: func_on_atom
    integer :: lcao, loc_lcao
    integer :: loc_ngwf
    integer :: funcs_before_atom
    type(FFTBOX_DATA) :: lcao_fftbox

    type(FUNCTIONS) :: lcao_on_grid(mdl%nsub)
    type(FUNCTIONS) :: orth_lcao_on_grid(mdl%nsub,1)
    complex(kind=DP), allocatable :: lcao_work(:,:,:)
    type(COEF) :: temp_coef, ngwf_coef
    type(SPAM3_EMBED) :: lcao_overlap
    type(spam3_embed) :: lcao_overlap_block_diag
    type(SPAM3_EMBED) :: ngwf_lcao_overlap
    type(SPAM3_EMBED) :: lcao_ngwf_overlap
    type(SPAM3_EMBED) :: invoverlap

    !integer     :: max_hotelling_iters
    !real(kind=dp) :: max_resid_hotelling

    integer :: iat, iat2, iatproc, iblk, jrow
    integer,allocatable :: seg_nzb(:)  ! Number of nonzero blocks on each seg
    integer,allocatable :: seg_nze(:)  ! Number of nonzero elements on each seg
    type(ELEMENT), allocatable :: elems_sfc(:)   ! List of atoms in SFC order
    !type(OVERLAP_LIST) :: overlaps

    type(DEM) :: pdos_weights_DEM
    type(DEM) :: pdos_weights_DEM2
    !    type(spam3) :: ngwf_overlap
    real(kind=dp) :: rval, rval2
    real(kind=dp) :: lcao_trace, inv_lcao_trace

    integer :: matfile

    type(SPAM3_EMBED) :: sp_overlap

    type(dem) :: lcao_ngwf_overlap_DEM
    !    type(dem) :: dem_prod_mat,dem_right_trans
    type(dem) :: lcao_overlap_dem

    real(kind=dp) :: rmaxsq, gmaxsq
    integer :: n1half,n2half,n3half

    integer :: jrow2
    ! agrecocmplx
    logical :: loc_cmplx
    ! temp variable to store real SW to copy to complex box
    real(kind=DP), allocatable :: lcao_fftbox_real(:,:,:)

    integer :: isub, jsub, lcao_num, num, ireg

    character :: suffix

    type(NGWF_REP)    :: sw_rep
    type(NGWF_HAM)    :: sw_ham
    type(EDFT_MODEL)  :: sw_edft
    type(DKERN)       :: sw_denskern
    type(DKERN)       :: sw_pur_denskern
    real(kind=dp), allocatable :: sw_lhxc_fine(:,:,:,:)
    type(FUNC_BASIS)  :: hub_proj_basis
    type(HUBBARD_MODEL) :: hub
    integer :: edft_maxiter
    real(kind=dp)     :: sw_mu(pub_num_spins,pub_num_kpoints)
    real(kind=dp)     :: sw_energy

    integer :: prev_atom, cur_atom
    real(kind=dp), allocatable ::  temp_block(:,:)
    integer :: i,j
    real(kind=dp) :: nelecs
    integer :: is
    type(DEM) :: sw_overlap_dens
    type(DEM), allocatable :: sw_eigs_dens(:)
    type(DEM) :: sw_ham_dens
    real(kind=dp), allocatable :: sw_eigen_en(:,:)
    real(kind=dp), allocatable :: sw_eigen_occ(:,:)
    type(DEM) :: sw_buffer_dens
    type(SPAM3) :: SW_KS
    type(SPAM3) :: SW_SKS
    complex(kind=dp) :: sw_temp_cmplx
    integer :: orbital_index
    type(DEM) :: sw_lowdin
    type(DEM) :: sw_S_evecs
    type(DEM) :: kernel_dens
    real(kind=dp), allocatable :: sw_S_evals(:)
    real(kind=dp), allocatable :: sw_s_occ(:,:)
    logical :: non_pos_def_overlap
    real(kind=dp) :: diag_el
    integer :: work_length
    real(kind=dp), allocatable :: overlap_eigenvalues(:)
    real(kind=dp), allocatable :: dsyev_work(:)

    integer :: iset, orig_atom, em, lcao_on_atom, lcao_in_species, shell
    integer :: distr_atom

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering eigenstates_pdos'

    ! JA: start timer
    call timer_clock('eigenstates_pdos',1)

    if(pub_on_root) then
       write(stdout,*)
       write(stdout,'(a)') utils_banner('=', &
            'Projected Density of States (pDOS) calculation')
    end if

    call date_and_time(date)
    file_header='Generated with ONETEP on '//date//'.'
    file_version=1.0_dp

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid(1)%iscmplx

    ! JA: check optional arguments
    !jmecmplx
    ! agrecocmplx: now made compatible, but need to check everything is ok!
    !call utils_assert(.not. ngwfs_on_grid%iscmplx, 'Error in eigenstates_pdos:&
    !     & subroutine not yet checked for complex NGWFs.')
    ! agrecocmplx
    temp_coef%iscmplx = loc_cmplx
    if (loc_cmplx) then
       temp_coef%z = (1.0_DP,0.0_DP)
    else
       temp_coef%d = 1.0_DP
    end if
    ngwf_coef%iscmplx = loc_cmplx
    if (loc_cmplx) then
       ngwf_coef%z = (0.0_DP,0.0_DP)
    else
       ngwf_coef%d = 0.0_DP
    end if



    if (present (nmax)) then
       max_n = nmax

    elseif(pub_pdos_max_n>0) then
       max_n=pub_pdos_max_n
    else
       ! JA: default is number of psincs across largest NGWF
       ! jcap: loop over regions to get the biggest overall
       max_n=0
       do isub=1,mdl%nsub
          max_n = max(max_n, ceiling(maxval(lcao_basis(isub)%spheres(:)%radius) / &
               min(mdl%cell%d1,mdl%cell%d2,mdl%cell%d3)) )
       end do
       call comms_reduce('MAX',max_n)
    endif

    num_n=((max_n-1)/pub_pdos_stride_n)+1

    ! JA: allocate qnl array. This is the wave-vectors array required for calculating the S-Bessel functions.
    allocate(qnl(max_n),stat=ierr)
    call utils_alloc_check('eigenstates_pdos','qnl',ierr)

    ! JA: init SW module if not already - with max kr
    n1half = mdl%fftbox%total_pt1/2+1
    n2half = mdl%fftbox%total_pt2/2+1
    n3half = mdl%fftbox%total_pt3/2+1

    gmaxsq = 2.0_DP*mdl%fftbox%recip_grid(5,n1half,n2half,n3half)
    ! jcap: loop over regions to get this maximum
    rmaxsq=0.0_DP
    do isub=1,mdl%nsub
       rmaxsq = max(rmaxsq, gmaxsq*(maxval(lcao_basis(isub)%spheres(:)%radius)**2) )
    end do

    call comms_reduce('MAX',rmaxsq)

    call sw_init(pub_pdos_max_l,max_n,rmaxsq)

    ! JA: calculate shift required to centre SW in FFTbox on a grid point
    fbcentre = real(mdl%fftbox%total_pt1-1,kind=DP)/2/mdl%fftbox%total_pt1 * mdl%fftbox%a1 + &
         real(mdl%fftbox%total_pt2-1,kind=DP)/2/mdl%fftbox%total_pt2 * mdl%fftbox%a2 + &
         real(mdl%fftbox%total_pt3-1,kind=DP)/2/mdl%fftbox%total_pt3 * mdl%fftbox%a3

    ! JA: decide how many electrons we will have in each spin channel at each k-point
    ! bearing in mind that this routine is presently limited to Gamma point only.
    if(pub_pdos_reduce_sws) then
       if(pub_pdos_construct_basis) then
          num_orbitals_on_proc = 0
          do isub=1,mdl%nsub
             num_orbitals_on_proc = num_orbitals_on_proc + &
                  ngwf_basis(isub)%num_on_proc(pub_my_proc_id)*((pub_pdos_max_l+1)**2)
          end do
       else
          num_orbitals_on_proc=0
          do isub=1,mdl%nsub
             do atom=1,mdl%regions(isub)%par%num_atoms_on_proc(pub_my_proc_id)
                do ll=0,pub_pdos_max_l
                   num_orbitals_on_proc = num_orbitals_on_proc + (2*ll + 1)
                end do
             end do
          end do
       end if
    else if(pub_pdos_pseudoatomic) then
       num_orbitals_on_proc=0
       ! jcap: loop over regions
       do isub=1,mdl%nsub
          do atom=1,mdl%regions(isub)%par%num_atoms_on_proc(pub_my_proc_id)
             orig_atom=mdl%regions(isub)%par%orig_atom(atom+mdl%regions(isub)%par%first_atom_on_proc(pub_my_proc_id)-1)
             num_orbitals_on_proc=num_orbitals_on_proc + &
                  mdl%regions(isub)%elements(orig_atom)%nfunctions
          end do
       end do

    else
       num_orbitals_on_proc=0
       ! jcap: loop over regions
       do isub=1,mdl%nsub
          do atom=1,mdl%regions(isub)%par%num_atoms_on_proc(pub_my_proc_id)
             do ll=0,pub_pdos_max_l
                num_orbitals_on_proc = num_orbitals_on_proc + (2*ll + 1)
             end do
          end do
       end do
       num_orbitals_on_proc=num_orbitals_on_proc*num_n
    end if

    ! JA: determine where all the orbitals are
    num_orbitals=num_orbitals_on_proc
    max_orbitals_on_proc=num_orbitals_on_proc
    call comms_reduce('MAX',max_orbitals_on_proc)
    call comms_reduce('SUM',num_orbitals)
    allocate(num_orbitals_on_procs(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('eigenstates_pdos','num_orbitals_on_procs',ierr)
    call comms_allgather(num_orbitals_on_procs,num_orbitals_on_proc)

    if(pub_pdos_reduce_sws) then
       if(pub_pdos_construct_basis) then
          num_orbitals=sum(mdl%regions(:)%par%nat*((pub_pdos_max_l+1)**2))
       end if
    end if

    if((.not.pub_pdos_reduce_sws).and.(.not.pub_pdos_pseudoatomic)) then
       num_orbitals=num_orbitals/num_n
    end if


    ! Create some work space
    ! agrecocmplx
    call data_fftbox_alloc(lcao_fftbox, mdl%fftbox%total_ld1, &
         mdl%fftbox%total_ld2, mdl%fftbox%total_pt3, &
         iscmplx=loc_cmplx)
    allocate(lcao_work(mdl%fftbox%total_ld1,mdl%fftbox%total_ld2,&
         mdl%fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('eigenstates_pdos','lcao_work',ierr)
    ! allocate sw_on_grid array
    ! agrecocmplx
    ! jcap: loop over regions
    do isub=1,mdl%nsub
       call data_functions_alloc(lcao_on_grid(isub), lcao_basis(isub)%size_on_grid, &
            iscmplx=loc_cmplx)
    end do

    ! JA: Transfer atoms to space-filling curve order
    allocate(elems_sfc(mdl%nat),stat=ierr)
    call utils_alloc_check('sparse_init','elems_sfc',ierr)
    allocate(seg_nzb(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('sparse_init','seg_nzb',ierr)
    allocate(seg_nze(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('sparse_init','seg_nze',ierr)
    iblk = 1
    do isub=1,mdl%nsub
       do proc=0,pub_total_num_procs-1
          do iatproc=1,mdl%regions(isub)%par%num_atoms_on_proc(proc)
             iat = iatproc + mdl%regions(isub)%par%first_atom_on_proc(proc) - 1
             elems_sfc(iblk) = mdl%regions(isub)%elements(iat)
             iblk = iblk + 1
          end do
       end do
    end do


    if(present(ham_type)) then
       select case(trim(adjustl(ham_type)))
          case('valence')
             suffix=''
          case('proj')
             suffix=''
          case('cond')
             suffix='c'
          case('joint')
             suffix='j'
          case('aux')
             suffix='a'
          case default
             call utils_abort('Error in eigenstates_pdos: &
                  &unrecognised ham_type: '//trim(adjustl(ham_type))//'.')
       end select
    else
          suffix=''
    end if

    suffix=rep%postfix(1:1)

    ! JA: Create the various overlap and inverse overlap matrices
    ! agrecocmplx: lcao_overlap will be real anyway, but when we need
    ! to compute SW-NGWF overlap, we need lcao_on_grid to be complex
    ! for the routines to work properly; currently having lcao_overlap
    ! complex instead of defining a complex copy
    lcao_overlap%structure = 'L$'
    call sparse_embed_create(lcao_overlap,iscmplx=loc_cmplx)
    lcao_ngwf_overlap%structure = 'P$'//trim(adjustl(suffix))
    call sparse_embed_create(lcao_ngwf_overlap,iscmplx=loc_cmplx)
    call sparse_embed_create(invoverlap, lcao_overlap, lcao_overlap)
    !    call sparse_create(invoverlap_lcao_ngwf, invoverlap, lcao_ngwf_overlap)
    ngwf_lcao_overlap%structure = 'P%'//trim(adjustl(suffix))
    call sparse_embed_create(ngwf_lcao_overlap,iscmplx=loc_cmplx)


    ! agrecocmplx
    call dense_create(pdos_weights_DEM,int(sparse_embed_num_rows(lcao_overlap)),&
         int(sparse_embed_num_cols(lcao_ngwf_overlap)),iscmplx=loc_cmplx)
    call dense_create(pdos_weights_DEM2,pdos_weights_DEM)
    ! agrecocmplx: allocate temp variables for complex case
    if (loc_cmplx) then
       allocate(lcao_fftbox_real(mdl%fftbox%total_ld1,mdl%fftbox%total_ld2,&
            mdl%fftbox%total_pt3),stat=ierr)
       call utils_alloc_check('eigenstates_pdos','lcao_fftbox_real',ierr)
    end if

    ! agrecocmplx
    call dense_create(lcao_ngwf_overlap_DEM,pdos_weights_DEM)

    old_atom = -999
    call data_set_to_zero(lcao_fftbox)
    do isub=1,mdl%nsub
       call data_set_to_zero(lcao_on_grid(isub))
    end do

    ! JA: Construct Angular Momentum resolved functions
    if(pub_on_root) then
       write(stdout,*)''
       write(stdout,'(a)', advance='no') 'Constructing AM resolved functions  ...'
    end if


    ! jcap: loop over regions
    do isub=1,mdl%nsub
       ! If we are using the pseudoatomic set, we don't need to set up SWs!
       if(.not.pub_pdos_pseudoatomic) then

       loc_lcao = 0
       do lcao=lcao_basis(isub)%first_on_proc(pub_my_proc_id), &
            lcao_basis(isub)%first_on_proc(pub_my_proc_id+1)-1

          loc_lcao = loc_lcao + 1
          centre = lcao_basis(isub)%spheres(loc_lcao)%centre
          radius = lcao_basis(isub)%spheres(loc_lcao)%radius

          ! JA: calculate fractional coordinates of sw centre
          frac%X = (centre .DOT. mdl%cell%b1) / two_pi * real(mdl%cell%total_pt1,kind=DP)
          frac%Y = (centre .DOT. mdl%cell%b2) / two_pi * real(mdl%cell%total_pt2,kind=DP)
          frac%Z = (centre .DOT. mdl%cell%b3) / two_pi * real(mdl%cell%total_pt3,kind=DP)

          ! JA: calculate displacement of centre from nearest grid point
          disp = mod(frac%X,1.0_DP) / mdl%fftbox%total_pt1 * mdl%fftbox%a1 + &
               mod(frac%Y,1.0_DP) / mdl%fftbox%total_pt2 * mdl%fftbox%a2 + &
               mod(frac%Z,1.0_DP) / mdl%fftbox%total_pt3 * mdl%fftbox%a3

          ! JA: calculate first grid point of tightbox wrt cell

          call basis_location_func_wrt_cell(sx,sy,sz,&
               lcao_basis(isub)%tight_boxes(loc_lcao),mdl%cell)

          ! JA: calculate location of tightbox wrt FFTbox and displacement of centre wrt
          ! FFTbox depending whether FFTbox side coincides with sim cell
          if(mdl%fftbox%coin1) then
             sx = 1
             disp%X = centre%X
          else
             sx = (mdl%fftbox%total_pt1-1)/2 - floor(frac%X) + sx
             disp%X = fbcentre%X + disp%X
          endif
          if(mdl%fftbox%coin2) then
             sy = 1
             disp%Y = centre%Y
          else
             sy = (mdl%fftbox%total_pt2-1)/2 - floor(frac%Y) + sy
             disp%Y = fbcentre%Y + disp%Y
          endif
          if(mdl%fftbox%coin3) then
             sz = 1
             disp%Z = centre%Z
          else
             sz = (mdl%fftbox%total_pt3-1)/2 - floor(frac%Z) + sz
             disp%Z = fbcentre%Z + disp%Z
          endif


          func_on_atom=lcao-lcao_basis(isub)%first_on_atom(lcao_basis(isub)%atom_of_func(lcao))+1

          if(pub_pdos_reduce_sws) then

             if(pub_pdos_construct_basis) then

                !JA: loop over embedding regions... I guess this is how to do it(!)
                do jsub=1,mdl%nsub

                   loc_ngwf=((func_on_atom-1)/((pub_pdos_max_l+1)**2)) &
                        + ngwf_basis(jsub)%first_on_atom(lcao_basis(isub)%atom_of_func(lcao)) &
                        - ngwf_basis(jsub)%first_on_proc(pub_my_proc_id) + 1

                   ! JA: Set azimuthal quantum number
                   ! In this case mod(func_on_atom-1,(pub_pdos_max_l+1)**2)+1
                   ! represents the index of the new func on an ngwf.
                   ll=ceiling(sqrt(real(mod(func_on_atom-1,(pub_pdos_max_l+1)**2)+1)))-1
                   ! magenetic quantum number
                   mm=(mod(func_on_atom-1,(pub_pdos_max_l+1)**2)+1)-(ll**2)-ll-1

                   ! JA: generate bessel function roots for current NGWF radius and l
                   qnl(:) = sw_bessel_zeros(1:max_n,ll) / radius

                   ! JA: Sum over principal quantum number
                   do nn=1,max_n,pub_pdos_stride_n!1

                      ! JA: Generate a spherical wave in the fft box
                      ! agrecocmplx: use temp copy to store real SW before copying to complex box
                      if (loc_cmplx) then
                         call sw_recp_generate(lcao_fftbox_real,lcao_work,ll,mm,qnl(nn),radius, &
                              disp,mdl%fftbox)
                         ! copy to complex
                         lcao_fftbox%z(:,:,:) = cmplx(lcao_fftbox_real(:,:,:),kind=DP)
                         ! real case: use real box directly
                      else
                         call sw_recp_generate(lcao_fftbox%d,lcao_work,ll,mm,qnl(nn),radius, & !jmecmplx
                              disp,mdl%fftbox)
                      end if

                      !JA: get fitting coeffs by dotting with an NGWF
                      call basis_dot_function_with_box(ngwf_coef, &
                           ngwfs_on_grid(jsub), lcao_fftbox, &
                           lcao_basis(isub)%spheres(loc_lcao), lcao_basis(isub)%tight_boxes(loc_lcao), &
                           sx, sy, sz, ngwf_basis(jsub)%spheres(loc_ngwf)%offset, &
                           mdl%cell, mdl%fftbox)

                      call data_fftbox_scale(lcao_fftbox, ngwf_coef)

                      ! JA: Sum PPDs with what is in box
                      call basis_add_function_to_box(lcao_fftbox, &
                           sx, sy, sz, &
                           lcao_basis(isub)%tight_boxes(loc_lcao), lcao_on_grid(isub), &
                           lcao_basis(isub)%spheres(loc_lcao), &
                           temp_coef, mdl%cell, mdl%fftbox)

                      ! JA: Put everything into PPDs
                      call basis_extract_function_from_box(lcao_on_grid(isub), &
                           lcao_fftbox, &
                           lcao_basis(isub)%spheres(loc_lcao),&
                           lcao_basis(isub)%tight_boxes(loc_lcao), sx,sy,sz, &
                           lcao_basis(isub)%spheres(loc_lcao)%offset, mdl%cell, mdl%fftbox)

                   end do

                end do

             else
                ! JA: Set azimuthal quantum number
                ll=ceiling(sqrt(real(func_on_atom)))-1
                ! magenetic quantum number
                mm=func_on_atom-(ll**2)-ll-1

                ! JA: generate bessel function roots for current NGWF radius and l
                qnl(:) = sw_bessel_zeros(1:max_n,ll) / radius

                ! JA: Sum over principal quantum number
                do nn=1,max_n,pub_pdos_stride_n!1

                   ! JA: Generate a spherical wave in the fft box
                   ! agrecocmplx: use temp copy to store real SW before copying to complex box
                   if (loc_cmplx) then
                      call sw_recp_generate(lcao_fftbox_real,lcao_work,ll,mm,qnl(nn),radius, &
                           disp,mdl%fftbox)
                      ! copy to complex
                      lcao_fftbox%z(:,:,:) = cmplx(lcao_fftbox_real(:,:,:),kind=DP)
                      ! real case: use real box directly
                   else
                      call sw_recp_generate(lcao_fftbox%d,lcao_work,ll,mm,qnl(nn),radius, & !jmecmplx
                           disp,mdl%fftbox)
                   end if

                   ! JA: Sum PPDs with what is in box
                   call basis_add_function_to_box(lcao_fftbox, &
                        sx, sy, sz, &
                        lcao_basis(isub)%tight_boxes(loc_lcao), lcao_on_grid(isub), &
                        lcao_basis(isub)%spheres(loc_lcao), &
                        temp_coef, mdl%cell, mdl%fftbox)

                   ! JA: Put everything into PPDs
                   call basis_extract_function_from_box(lcao_on_grid(isub), &
                        lcao_fftbox, &
                        lcao_basis(isub)%spheres(loc_lcao),&
                        lcao_basis(isub)%tight_boxes(loc_lcao), sx,sy,sz, &
                        lcao_basis(isub)%spheres(loc_lcao)%offset, mdl%cell, mdl%fftbox)

                end do
             end if

             ! JA: Clean up
             call basis_clean_function(lcao_on_grid(isub), lcao_basis(isub)%spheres(loc_lcao), &
                  mdl%cell, mdl%fftbox)

          else

             ! JA: set quantum numbers
             nn=ceiling(real(func_on_atom)/real((pub_pdos_max_l+1)**2))
             ll=ceiling(sqrt(real(mod(func_on_atom-1,(pub_pdos_max_l+1)**2))+1.0))-1
             mm=func_on_atom-((nn-1)*((pub_pdos_max_l+1)**2))-(ll**2)-ll-1

             ! JA: Don't sum over principal quantum number!

             ! JA: Generate a spherical wave in the fft box
             ! agrecocmplx
             if (loc_cmplx) then
                call sw_recp_generate(lcao_fftbox_real,lcao_work,ll,mm,&
                     sw_bessel_zeros(((nn-1)*pub_pdos_stride_n)+1,ll)/radius, &
                     radius, disp,mdl%fftbox)
                ! copy to complex
                lcao_fftbox%z(:,:,:) = cmplx(lcao_fftbox_real(:,:,:),kind=DP)
             else

                call sw_recp_generate(lcao_fftbox%d,lcao_work,ll,mm,&
                     sw_bessel_zeros(((nn-1)*pub_pdos_stride_n)+1,ll)/radius,radius, &
                     disp,mdl%fftbox)
             end if

             ! JA: Put everything into PPDs
             call basis_extract_function_from_box(lcao_on_grid(isub), &
                  lcao_fftbox, &
                  lcao_basis(isub)%spheres(loc_lcao),&
                  lcao_basis(isub)%tight_boxes(loc_lcao), sx,sy,sz, &
                  lcao_basis(isub)%spheres(loc_lcao)%offset, mdl%cell, mdl%fftbox)

             ! JA: Clean up
             call basis_clean_function(lcao_on_grid(isub), &
                  lcao_basis(isub)%spheres(loc_lcao), &
                  mdl%cell, mdl%fftbox)

          end if

       end do

    else

       call ngwfs_initialise_ngwfs(lcao_on_grid(isub),lcao_basis(isub),mdl,isub)

       if (pub_on_root.and.(pub_output_detail<NORMAL)) &
            write(stdout,'(a)') '... done'

    end if

    deallocate(qnl,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','qnl',ierr)


! JA: Uncomment to make visualisation of SW functions for debugging:
       call visual_ngwfs(lcao_on_grid(isub), lcao_basis(isub), "lcao", &
            mdl%regions(isub)%elements, mdl%cell, mdl%fftbox, mdl%regions(isub)%par)

end do

    if(pub_on_root.and..not.pub_pdos_pseudoatomic) write(stdout,*) ' done'

    if(lcao_overlap%iscmplx) then
       call sparse_embed_scale(lcao_overlap,cmplx_0)
    else
       call sparse_embed_scale(lcao_overlap,0.0_DP)
    end if

    if(pub_on_root) then
       write(stdout,*)''
       write(stdout,'(a)', advance='no') 'Performing overlap integrals ...'
    end if

    ! jcap: loop over regions
    do isub=1,mdl%nsub
       do jsub=1,mdl%nsub
          ! JA: SW-SW overlap ("\Lambda")
          call function_ops_brappd_ketppd(lcao_overlap%m(isub,jsub), &
               lcao_on_grid(isub), lcao_basis(isub), &
               lcao_on_grid(jsub), lcao_basis(jsub), mdl%cell)

       end do
    end do


    if (pub_aug) then
       sp_overlap%structure='Pr$'//trim(adjustl(suffix))

       call sparse_embed_create(sp_overlap,iscmplx=loc_cmplx)
       ! JA: calculate the ngwf-projector overlap matrix
       ! jcap: loop over regions
       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub

             if (pub_any_nl_proj.or.pub_paw) then
                if ((.not.pub_realspace_projectors).or.pub_paw) then
                   ! agrecokpt: need to make this k-dependent?
                   call projectors_func_ovlp_box(sp_overlap%m(isub,jsub), &
                        lcao_on_grid(isub), lcao_basis(isub), &
                        proj_basis(jsub), nl_projectors(jsub), mdl%fftbox,mdl%cell)
                else
                   call function_ops_brappd_ketppd(sp_overlap%m(isub,jsub), &
                        lcao_on_grid(isub), lcao_basis(isub), &
                        nl_projectors(jsub)%projs_on_grid, proj_basis(jsub), mdl%cell)
                end if
             end if
             call augmentation_overlap(lcao_ngwf_overlap%m(isub,jsub),&
                  mdl%regions(isub)%pseudo_sp,mdl%regions(jsub)%paw_sp,&
                  sp_overlap%m(isub,jsub),rep%sp_overlap%m(isub,jsub))

          end do
       end do


    end if

    ! JA: Check for zeros on diagonal... indicating a function with zero weight...
    !          need to deal with this or we'll have a non-positive-definite overlap matrix.
    non_pos_def_overlap=.false.
    do isub=1,mdl%nsub
       do lcao=lcao_basis(isub)%first_on_proc(pub_my_proc_id), &
            lcao_basis(isub)%first_on_proc(pub_my_proc_id+1)-1

          call sparse_get_element(diag_el, lcao_overlap%m(isub,isub), lcao,lcao)
          non_pos_def_overlap=non_pos_def_overlap.or.(diag_el<epsilon(diag_el))

       end do
    end do

    ! JA: Just going to consider it to be non-positive definite, since we were suffering from numerical problems
    non_pos_def_overlap=.true.

    call comms_reduce('OR',non_pos_def_overlap)

    ! JA: This will orthonormalize the SW functions on an atom, if enabled.
    if(pub_pdos_orth_atom_blocks) then

       ! Copied this from below (electronic minimization section) for this experiment
       ! so that we have access to the "Dp" structure.
       call ngwf_rep_create(sw_rep,'p',mdl,loc_cmplx)

       lcao_overlap_block_diag%structure = "Dp"
       call sparse_embed_create(lcao_overlap_block_diag)
       call sparse_embed_copy(lcao_overlap_block_diag,lcao_overlap)


       ! Number of orbitals on atoms
       allocate(num_orbitals_on_atom(1:sum(mdl%regions(:)%par%nat)),stat=ierr)
       call utils_alloc_check('eigenstates_pdos','num_orbitals_on_atom',ierr)

       iat2=0
       do isub=1,mdl%nsub
          do iat=1,mdl%regions(isub)%par%nat
             iat2=iat2+1
             num_orbitals_on_atom(iat2) = lcao_basis(isub)%num_on_atom(iat)
          end do
       end do

       allocate(temp_block(1:maxval(num_orbitals_on_atom),1:maxval(num_orbitals_on_atom)),stat=ierr)
       call utils_alloc_check('eigenstates_pdos','temp_block',ierr)

       if(non_pos_def_overlap) then
          work_length=-1
          allocate(overlap_eigenvalues(size(temp_block,1)),stat=ierr)
          call utils_alloc_check('eigenstates_pdos','overlap_eigenvalues',ierr)
          allocate(dsyev_work(1),stat=ierr)
          call utils_alloc_check('eigenstates_pdos','dsyev_work',ierr)
          call dsyev('V', 'U', size(temp_block,1), temp_block, size(temp_block,1), overlap_eigenvalues,&
               dsyev_work, work_length, ierr)
          call utils_assert(ierr==0,'Error in eigenstates_pdos : dsyev reported a nonzero status in workspace query')
          work_length=int(dsyev_work(1))
          deallocate(dsyev_work,stat=ierr)
          call utils_dealloc_check('eigenstates_pdos','dsyev_work',ierr)
          allocate(dsyev_work(work_length),stat=ierr)
          call utils_alloc_check('eigenstates_pdos','dsyev_work',ierr)
       end if

       loc_lcao = 0
       prev_atom = -1
       do isub=1,mdl%nsub
          do lcao=lcao_basis(isub)%first_on_proc(pub_my_proc_id), &
               lcao_basis(isub)%first_on_proc(pub_my_proc_id+1)-1

             loc_lcao = loc_lcao + 1

             cur_atom = lcao_basis(isub)%atom_of_func(lcao)

             if(cur_atom /= prev_atom) then

                temp_block=0.0_dp

                call sparse_get_block(temp_block(1:num_orbitals_on_atom(cur_atom),1:num_orbitals_on_atom(cur_atom)),&
                     lcao_overlap_block_diag%m(isub,isub),cur_atom,cur_atom)

                ! JA: If we know there's a zero on the diagonal, we have to do this differently...
                if(.not.non_pos_def_overlap) then
                   call dpotrf('L', num_orbitals_on_atom(cur_atom), temp_block, size(temp_block,1), ierr)
                   call utils_assert(ierr==0,'Error in eigenstates_pdos : dpotrf reported a nonzero status')

                   do i=1,num_orbitals_on_atom(cur_atom)
                      do j=1,i-1
                         temp_block(j,i)=0.0_dp
                      end do
                   end do

                   call dtrtri('L','N', num_orbitals_on_atom(cur_atom), temp_block, size(temp_block,1), ierr)
                   call utils_assert(ierr==0,'Error in eigenstates_pdos : dtrtri reported a nonzero status')
                else ! ... using an eigensolver.

                   call dsyev('V', 'U', num_orbitals_on_atom(cur_atom), temp_block, size(temp_block,1),&
                        overlap_eigenvalues,dsyev_work, work_length, ierr)

                   call utils_assert(ierr==0,'Error in eigenstates_pdos : dsyev reported a nonzero status')

                   do i=1,num_orbitals_on_atom(cur_atom)
                      if(overlap_eigenvalues(i)>epsilon(overlap_eigenvalues(i))) then
                         overlap_eigenvalues(i)=1.0_dp/overlap_eigenvalues(i)
                      else
                         overlap_eigenvalues(i)=0.0_dp
                      end if
                      temp_block(1:num_orbitals_on_atom(cur_atom),i)=&
                           temp_block(1:num_orbitals_on_atom(cur_atom),i)*sqrt(overlap_eigenvalues(i))
                   end do

                   temp_block(1:num_orbitals_on_atom(cur_atom),1:num_orbitals_on_atom(cur_atom))=&
                        matmul(temp_block(1:num_orbitals_on_atom(cur_atom),1:num_orbitals_on_atom(cur_atom)),&
                        transpose(temp_block(1:num_orbitals_on_atom(cur_atom),1:num_orbitals_on_atom(cur_atom))))

                end if


                temp_block(1:num_orbitals_on_atom(cur_atom),1:num_orbitals_on_atom(cur_atom)) = &
                     transpose(temp_block(1:num_orbitals_on_atom(cur_atom),1:num_orbitals_on_atom(cur_atom)))

                call sparse_put_block(temp_block(1:num_orbitals_on_atom(cur_atom),1:num_orbitals_on_atom(cur_atom)),&
                     lcao_overlap_block_diag%m(isub,isub),cur_atom,cur_atom)

             end if

             prev_atom = cur_atom

          end do
       end do

       if(non_pos_def_overlap) then
          deallocate(overlap_eigenvalues,stat=ierr)
          call utils_dealloc_check('eigenstates_pdos','overlap_eigenvalues',ierr)
          deallocate(dsyev_work,stat=ierr)
          call utils_dealloc_check('eigenstates_pdos','dsyev_work',ierr)
       end if

       deallocate(num_orbitals_on_atom,stat=ierr)
       call utils_dealloc_check('eigenstates_pdos','num_orbitals_on_atom',ierr)

       deallocate(temp_block,stat=ierr)
       call utils_dealloc_check('eigenstates_pdos','temp_block',ierr)


       do isub=1,mdl%nsub
          call data_functions_alloc(orth_lcao_on_grid(isub,1), lcao_basis(isub)%size_on_grid, &
               iscmplx=loc_cmplx)
          call data_set_to_zero(orth_lcao_on_grid(isub,1))

          call function_ops_sum_ppd_funcs(orth_lcao_on_grid(isub,1:1), &
               lcao_basis(isub), [lcao_overlap_block_diag%m(isub,isub)], 1, 1, lcao_overlap_block_diag%m(isub,isub), &
               lcao_on_grid(isub), lcao_basis(isub))

          call data_functions_copy(lcao_on_grid(isub),orth_lcao_on_grid(isub,1))

          call data_functions_dealloc(orth_lcao_on_grid(isub,1))

          call sparse_embed_destroy(lcao_overlap_block_diag)

       end do

       call ngwf_rep_destroy(sw_rep)

    end if

       ! JA: SW-SW overlap ("\Lambda")

       do isub=1,mdl%nsub
          do jsub=1,mdl%nsub
             call function_ops_brappd_ketppd(lcao_overlap%m(isub,jsub), &
                  lcao_on_grid(isub), lcao_basis(isub), &
                  lcao_on_grid(jsub), lcao_basis(jsub), mdl%cell)
          end do
       end do

       if (pub_aug) then
          do isub=1,mdl%nsub
             do jsub=1,mdl%nsub
                call augmentation_overlap(lcao_overlap%m(isub,jsub),&
                     mdl%regions(isub)%pseudo_sp,mdl%regions(jsub)%paw_sp,&
                     sp_overlap%m(isub,jsub))
             end do
          end do
       end if


    if(lcao_ngwf_overlap%iscmplx) then
       call sparse_embed_scale(lcao_ngwf_overlap,cmplx_0)
    else
       call sparse_embed_scale(lcao_ngwf_overlap,0.0_DP)
    end if

    ! jcap: loop over regions

    do isub=1,mdl%nsub
       do jsub=1,mdl%nsub
          ! JA: SW-NGWF overlap
          call function_ops_brappd_ketppd(lcao_ngwf_overlap%m(isub,jsub), &
               lcao_on_grid(isub), lcao_basis(isub), &
               ngwfs_on_grid(jsub), ngwf_basis(jsub), mdl%cell)
       end do
    end do
    if(invoverlap%iscmplx) then
       call sparse_embed_scale(invoverlap,cmplx_0)
    else
       call sparse_embed_scale(invoverlap,0.0_DP)
    end if

    ! JA: Transpose into NGWF-SW overlap before (possible) augmentation.
    call sparse_embed_transpose(ngwf_lcao_overlap,lcao_ngwf_overlap)

    ! JA: Now deal with Augmentation region contribution to overlaps
    !sp_overlap%structure='Pr$'//trim(adjustl(suffix))
    ! agrecocmplx
       if (pub_aug) then
          do isub=1,mdl%nsub
             do jsub=1,mdl%nsub
                call augmentation_overlap(lcao_ngwf_overlap%m(isub,jsub),&
                     mdl%regions(isub)%pseudo_sp,mdl%regions(jsub)%paw_sp,&
                     sp_overlap%m(isub,jsub),rep%sp_overlap%m(isub,jsub))
             end do
          end do
          call sparse_embed_destroy(sp_overlap)
       end if










    ! JA: Opting to do inverse overlap with dense_invert because of convergence problems...
    ! agrecocmplx
    call dense_create(lcao_overlap_dem,int(sparse_embed_num_rows(lcao_overlap)), &
         int(sparse_embed_num_cols(lcao_overlap)),iscmplx=loc_cmplx)
    call dense_convert(lcao_overlap_dem,lcao_overlap)


    if(.not.non_pos_def_overlap) then
       call dense_invert(lcao_overlap_dem)
    end if



!    rval=sparse_trace(invoverlap)
    if(pub_on_root) write(stdout,*) ' done'


    if(pub_pdos_output_basis) then
       do isub=1,mdl%nsub
          call restart_ngwfs_tightbox_output(lcao_on_grid(isub), lcao_basis(isub), &
               mdl%cell, mdl%fftbox, mdl%elements, "sw_basis.tightbox_ngwfs",mdl%regions(isub))
       end do
    end if

    ! jcap: get total number of lcao
    lcao_num=sum(lcao_basis(:)%num)
    ! jcap: get total number of NGWFs
    num=sum(ngwf_basis(:)%num)






    !JA: setup arrays with metadata about lcao basis functions : need this for optados output and
    !        to make the histogram later. Doing it now in case we're doing pdos_lcao_optimize.
    ! orbital_order gives the original position of the atom of this function in the input file
    ! orbital_species gives the species of the atom of this function
    ! orbital_l gives the l quantum number of this function
    ! orbital_ion gives the index of the atom of this function IN A SPECIES!
    !
    if (pub_on_root) then
       allocate(orbital_order(num_orbitals),stat=ierr)
       call utils_alloc_check('eigenstates_pdos','orbital_order',ierr)
       allocate(orbital_species(num_orbitals),stat=ierr)
       call utils_alloc_check('eigenstates_pdos','orbital_species',ierr)
       allocate(orbital_ion(num_orbitals),stat=ierr)
       call utils_alloc_check('eigenstates_pdos','orbital_ion',ierr)
       allocate(orbital_l(num_orbitals),stat=ierr)
       call utils_alloc_check('eigenstates_pdos','orbital_l',ierr)
       orbital_species=-1
       orbital_ion=-1
       orbital_l=-1
       orbital_order=-1

       num_species=maxval((/(mdl%elements(atom)%species_number,atom=1,sum(mdl%regions(:)%par%nat))/))
       allocate(atom_in_species(sum(mdl%regions(:)%par%nat)),stat=ierr)
       call utils_alloc_check('eigenstates_pdos','atom_in_species',ierr)
       allocate(count_in_species(num_species),stat=ierr)
       call utils_alloc_check('eigenstates_pdos','count_in_species',ierr)

       atom_in_species=0
       count_in_species=0
       do isub=1,mdl%nsub
          do atom=1,mdl%regions(isub)%par%nat
             orig_atom = mdl%regions(isub)%par%orig_atom(atom)
             count_in_species(mdl%elements(orig_atom)%species_number)=count_in_species(mdl%elements(orig_atom)%species_number)+1
             atom_in_species(orig_atom)=count_in_species(mdl%elements(orig_atom)%species_number)
          end do
       end do

       if(pub_pdos_reduce_sws) then
          lcao=0
          do isub=1,mdl%nsub
             do atom=1,mdl%regions(isub)%par%nat
                distr_atom=mdl%regions(isub)%par%distr_atom(atom)
                do jrow=(distr_atom-1)*((pub_pdos_max_l+1)**2)+1, distr_atom*((pub_pdos_max_l+1)**2)
                   lcao=lcao+1
                   orbital_order(lcao)=jrow
                end do
                atom_ind=mdl%regions(isub)%elements(atom)%global_atom_number
                orbital_species((atom_ind-1)*((pub_pdos_max_l+1)**2)+1:atom_ind*((pub_pdos_max_l+1)**2))=&
                     mdl%regions(isub)%elements(atom)%species_number
                orbital_ion((atom_ind-1)*((pub_pdos_max_l+1)**2)+1:atom_ind*((pub_pdos_max_l+1)**2))=atom_in_species(atom_ind)!atom
                orbital_l((atom_ind-1)*((pub_pdos_max_l+1)**2)+1:atom_ind*((pub_pdos_max_l+1)**2)) = &
                     (/((/(ll,l2=1,((ll*2)+1))/),ll=0,pub_pdos_max_l)/)
             end do

          end do
          !          end if

       else if(pub_pdos_pseudoatomic) then

          iset=1 ! JA: normal NGWFs -> maybe we should allow this to change for conduction / TDDFT? Would this ever be needed?

          lcao=0
          do jsub=1,mdl%nsub
             do orig_atom=1,mdl%regions(jsub)%par%nat
                distr_atom=mdl%regions(jsub)%par%distr_atom(orig_atom)
                do isub=1,mdl%nsub
                   do jrow=lcao_basis(isub)%first_on_atom(distr_atom), &
                        lcao_basis(isub)%num_on_atom(distr_atom)+lcao_basis(isub)%first_on_atom(distr_atom)-1
                      lcao_on_atom = jrow - lcao_basis(isub)%first_on_atom(distr_atom) + 1
                      lcao=lcao+1
                      orbital_species(lcao)=mdl%regions(isub)%elements(orig_atom)%species_number
                      orbital_order(lcao)=jrow
                      lcao_in_species = 0
                      do shell=1,mdl%regions(isub)%radial_ngwfs(orbital_species(lcao),iset)%nshells
                         do em=-mdl%regions(isub)%radial_ngwfs(orbital_species(lcao),iset)%angmom(shell),&
                              mdl%regions(isub)%radial_ngwfs(orbital_species(lcao),iset)%angmom(shell)
                            lcao_in_species = lcao_in_species + 1
                            if (lcao_in_species.eq.lcao_on_atom) exit
                         enddo
                         if (lcao_in_species.eq.lcao_on_atom) exit
                      enddo

                      orbital_ion(lcao)=atom_in_species(orig_atom)
                      orbital_l(lcao)=mdl%regions(isub)%radial_ngwfs(orbital_species(lcao),iset)%angmom(shell)
                   end do
                end do
             end do
          end do

       else

          lcao=0
          do isub=1,mdl%nsub
             do atom=1,mdl%regions(isub)%par%nat
                distr_atom=mdl%regions(isub)%par%distr_atom(atom)

                do jrow=(distr_atom-1)*((pub_pdos_max_l+1)**2)+1, distr_atom*((pub_pdos_max_l+1)**2)
                   lcao=lcao+1
                   orbital_order(lcao)=jrow
                end do

                orbital_species((atom-1)*((pub_pdos_max_l+1)**2)+1:atom*((pub_pdos_max_l+1)**2))=&
                     &mdl%elements(atom)%species_number
                orbital_ion((atom-1)*((pub_pdos_max_l+1)**2)+1:atom*((pub_pdos_max_l+1)**2))=atom_in_species(atom)!atom
                orbital_l((atom-1)*((pub_pdos_max_l+1)**2)+1:atom*((pub_pdos_max_l+1)**2)) = &
                     (/((/((/(ll,l2=1,((ll*2)+1))/),ll=0,pub_pdos_max_l)/),nn=1,1)/)

             end do
          end do
       end if

    end if




!     ! JA: allocate pdos weights on root
!     if(pub_on_root) then
!        if(pub_pdos_reduce_sws) then
!           allocate(pdos_weights(lcao_num,num,pub_num_spins),stat=ierr)
!           call utils_alloc_check('eigenstates_pdos','pdos_weights',ierr)
!           pdos_weights=0.0_dp
!        else
!           allocate(pdos_weights(lcao_num/num_n,num,pub_num_spins),stat=ierr)
!           call utils_alloc_check('eigenstates_pdos','pdos_weights',ierr)
!           pdos_weights=0.0_dp
!        end if
!     end if


    ! JA: This will perform an inner loop optimization with the SW basis, if requested
    ! and then perform a PDOS.
    ! This is experimental functionality and may be removed in a later revision of ONETEP, depending how
    ! useful it is.
    ! search for
    !                                       end if lcao_optimize
    ! to find the end of this section

!     lcao_optimize: if(pub_pdos_lcao_optimize) then
!
!
!        call ngwf_rep_create(sw_rep,'p',mdl%elements,mdl%cell,loc_cmplx, par)
!
!        ! JA: Work out number of electrons. I stole this from energy_and_force_mod...
!        ! ... no idea what the significance of the 4 is!
!        nelecs = 0.0_dp
!        do atom=1,par%nat
!           nelecs = nelecs + &
!                real(nint(4*(mdl%elements(atom)%ion_charge)), kind=DP)
!        end do
!        nelecs = nelecs/4.0_dp - pub_charge
!        sw_rep%n_occ(UP,PUB_1K) = (nelecs + pub_spin) / 2
!        if (pub_num_spins > 1) then
!           sw_rep%n_occ(DN,PUB_1K) = (nelecs - pub_spin) / 2
!        end if
!
!
!        ! moved call to create sw_rep to above ( experiment in overlap)
!
!        call ngwf_ham_create(sw_ham,sw_rep)
!
!        if (pub_edft) then
!           call edft_create(sw_edft, sw_ham, lcao_num, nelecs)
!        end if
!
!        call kernel_create(sw_denskern, 'K'//sw_rep%postfix, loc_cmplx)
!        call kernel_create(sw_pur_denskern, 'K'//sw_rep%postfix, loc_cmplx)
!
!        allocate(sw_lhxc_fine(mdl%fine_grid%ld1,&
!             mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12, pub_num_spins),stat=ierr)
!        call utils_alloc_check('eigenstates_pdos','sw_lhxc_fine',ierr)
!
!        call data_functions_alloc(sw_rep%ngwfs_on_grid, lcao_basis%size_on_grid, &
!             iscmplx=loc_cmplx)
!        call data_functions_copy(sw_rep%ngwfs_on_grid,lcao_on_grid)
!
!        call electronic_init_denskern(sw_rep, sw_ham, sw_denskern, &
!             sw_lhxc_fine, lcao_basis, proj_basis, nl_projectors, &
!             hub_proj_basis, hub, mdl, force_no_IO=.true., force_init_pot=.true., &
!             allow_pseuinvS=.true.) ! JA: allow pseudoinverse S
!
!        sw_mu=0.0_dp
!        sw_energy = electronic_energy(sw_denskern, sw_pur_denskern, sw_ham, &
!             sw_lhxc_fine, sw_mu, sw_edft, sw_rep, lcao_basis, hub_proj_basis, hub, &
!             mdl, lnv_threshold_orig, pub_firstit_lnv, kernel_update=.true., &
!             force_no_IO=.true.)
!
!        if(pub_pdos_output_swopt_kernham) then
!           call restart_kernel_write(sw_denskern%kern,file_extension='sw_basis.dkn')
!           call restart_hamiltonian_write(sw_ham%ham,file_extension='sw_basis.ham')
!        end if
!
!
!        ! JA: construct Lowdin factor
!        !          S_lowdin*evecs
!
!        call dense_create(sw_overlap_dens,lcao_num,lcao_num, &
!             iscmplx=loc_cmplx)
!
!        call dense_convert(sw_overlap_dens,sw_rep%overlap)
!
!        call dense_create(sw_S_evecs,lcao_num,lcao_num, &
!             iscmplx=loc_cmplx)
!        call dense_create(sw_lowdin,lcao_num,lcao_num, &
!             iscmplx=loc_cmplx)
!
!        allocate(sw_S_evals(lcao_num),stat=ierr)
!        call utils_alloc_check('eigenstates_pdos','sw_S_evals',ierr)
!
!        call dense_normal_eigensolve(lcao_num,sw_S_evals,sw_overlap_dens,sw_S_evecs,&
!             arg_orfac=pub_eigensolver_orfac, arg_abstol=pub_eigensolver_abstol, rtype=1)
!
!
!        ! JA: construct S^(1/2)
!        call dense_scale(sw_lowdin, 0.0_dp)
!        do orbital_index= 1,lcao_num
! !          call utils_assert(sw_S_evals(orbital_index)>0.0_dp,&
! !               'Error in eigenstates_pdos: SW overlap matrix is not positive definite! Index: ',orbital_index)
!           if(sw_S_evals(orbital_index)>epsilon(0.0_dp)) then
!              sw_S_evals(orbital_index) = sqrt(sw_S_evals(orbital_index))
!           else
!              sw_S_evals(orbital_index) = 0.0_dp
!           end if
!        end do
!        call edft_matrix_from_eigs(sw_lowdin, sw_S_evecs, sw_S_evals, lcao_num)
!
!        call dense_create(sw_buffer_dens,lcao_num,lcao_num, &
!             iscmplx=loc_cmplx)
!
!        do orbital_index= 1,lcao_num
!           if(sw_S_evals(orbital_index)>epsilon(0.0_dp)) then
!              sw_S_evals(orbital_index) = 1.0_dp/sw_S_evals(orbital_index)
!           else
!              sw_S_evals(orbital_index) = 0.0_dp
!           end if
!        end do
!
! !       call dense_copy(sw_buffer_dens,sw_lowdin)
! !       call dense_invert(sw_buffer_dens)
!
!        call edft_matrix_from_eigs(sw_buffer_dens, sw_S_evecs, sw_S_evals, lcao_num)
!
!
!        call dense_create(sw_ham_dens,lcao_num,lcao_num, &
!             iscmplx=loc_cmplx)
!
!        allocate(sw_eigs_dens(pub_num_spins),stat=ierr)
!        call utils_alloc_check('eigenstates_pdos','sw_eigs_dens',ierr)
!        do is=1,pub_num_spins
!           ! agrecocmplx: complex eigenstates for complex matrices
!           call dense_create(sw_eigs_dens(is),lcao_num,lcao_num, &
!                iscmplx=loc_cmplx)
!        end do
!
!        ! Create arrays to store eigenvalues and kernel occupations
!        ! agrecocmplx: these should stay real anyway
!        allocate(sw_eigen_en(lcao_num,pub_num_spins),stat=ierr)
!        call utils_alloc_check('eigenstates_pdos','sw_eigen_en',ierr)
!        allocate(sw_eigen_occ(lcao_num,pub_num_spins),stat=ierr)
!        call utils_alloc_check('eigenstates_pdos','sw_eigen_occ',ierr)
!
!
!        do is=1,pub_num_spins
!           ! lr408: Choose correct Hamiltonian for conduction calculation
!           if(present(ham_type)) then
!              if (ham_type=='cond') then
!                 call dense_convert(sw_ham_dens,sw_ham%cond_non_proj_ham(is))
!              else
!                 call dense_convert(sw_ham_dens,sw_ham%ham(is))
!              end if
!           else
!              call dense_convert(sw_ham_dens,sw_ham%ham(is))
!           end if
!
!           call dense_scale(sw_S_evecs,0.0_dp)
!           call dense_product(sw_S_evecs, sw_buffer_dens,sw_ham_dens)
!           call dense_scale(sw_ham_dens,0.0_dp)
!           call dense_product(sw_ham_dens,sw_S_evecs,sw_buffer_dens)
!
!           call dense_normal_eigensolve(lcao_num,sw_eigen_en(:,is),sw_ham_dens,sw_eigs_dens(is),&
!                arg_orfac=pub_eigensolver_orfac, arg_abstol=pub_eigensolver_abstol, rtype=1)
!        end do
!
!        call dense_destroy(sw_S_evecs)
!
!        deallocate(sw_S_evals,stat=ierr)
!        call utils_dealloc_check('eigenstates_pdos','sw_S_evals',ierr)
!
!        call sparse_create(sw_KS, sw_denskern%kern%m(1,1),sw_rep%overlap)
!        call sparse_create(sw_SKS, sw_rep%overlap,sw_KS)
!
!
!        do is=1,pub_num_spins
!
!
!           ! agrecocmplx: complex eigenstates for complex matrices
!           call dense_create(kernel_dens,lcao_num,lcao_num, &
!                iscmplx=loc_cmplx)
!
!           call dense_convert(kernel_dens,sw_denskern%kern%m(is,PUB_1K))
!
!           ! ndmh: now calculate occupancy of each state by projection onto
!           ! ndmh: valence density matrix for this spin
!           if(present(ham_type)) then
!              if ((ham_type=='valence').or.(ham_type=='cond').or. &
!                   (ham_type == 'proj').or.(ham_type=='aux')) then
!                 call sparse_product(sw_KS,sw_denskern%kern%m(is,PUB_1K),sw_rep%overlap)
!                 call sparse_product(sw_SKS,sw_rep%overlap,sw_KS)
!              else if (ham_type=='joint') then
!                 call utils_abort('joint basis not supported with optimized PDOS at the moment.')
!              end if
!           else
!              call sparse_product(sw_KS,sw_denskern%kern%m(is,PUB_1K),sw_rep%overlap)
!              call sparse_product(sw_SKS,sw_rep%overlap,sw_KS)
!           end if
!
!
!           call dense_product(sw_buffer_dens,kernel_dens,sw_lowdin)
!           call dense_product(kernel_dens,sw_lowdin,sw_buffer_dens)
!
!           call dense_product(sw_buffer_dens,sw_eigs_dens(is),kernel_dens,opA='T')
!           call dense_product(sw_overlap_dens,sw_buffer_dens,sw_eigs_dens(is))
!           do orbital_index=1,lcao_num
!              ! agrecocmplx: need temporary complex to get complex element
!              if (loc_cmplx) then
!                 call dense_get_element(sw_temp_cmplx,sw_overlap_dens, &
!                      orbital_index,orbital_index)
!                 ! agrecocmplx: occupancies are real anyway
!                 sw_eigen_occ(orbital_index,is) = real(sw_temp_cmplx,kind=DP)
!              else
!                 call dense_get_element(sw_eigen_occ(orbital_index,is),sw_overlap_dens, &
!                      orbital_index,orbital_index)
!              end if
!           end do
!
!           call dense_destroy(kernel_dens)
!        end do
!
!
!
!        do is=1,pub_num_spins
!           ! JA: sw_buffer_dens will be the weights
! !          call dense_product(sw_buffer_dens, sw_lowdin, sw_eigs_dens(is))
!
!           call dense_copy(sw_buffer_dens, sw_eigs_dens(is))
!
!           do icol=1,lcao_num
!              do jrow=1,lcao_num
!                 call dense_get_element(rval,sw_buffer_dens,jrow,icol)
!                 if(pub_on_root) then
!                    pdos_weights(jrow,icol,is) = rval**2
!                 end if
!              end do
!           end do
!
!           if(pub_on_root) then
!               if(pub_num_spins==1) then
!                   open(newunit=matfile, file="pdos_weights_swopt.dat",recl=lcao_num*50)
!               else
!                 if(spin==UP) then
!                   open(newunit=matfile, file="pdos_weights_swopt_UP.dat",recl=lcao_num*50)
!                 elseif(spin==DN) then
!                    open(newunit=matfile, file="pdos_weights_swopt_DOWN.dat",recl=lcao_num*50)
!                 end if
!               end if
!                 do icol=1,lcao_num
!                    do jrow=1, lcao_num
!                       write(matfile,'(f18.10,2x)',advance="NO") pdos_weights(jrow,icol,is)
!                    end do
!                    write(matfile,*)
!                 end do
!                 close(matfile)
!           end if
!
!        end do
!
!
!        allocate(sw_s_occ(pub_num_spins,pub_num_kpoints),stat=ierr)
!        call utils_alloc_check('eigenstates_pdos','sw_s_occ',ierr)
!        sw_s_occ=0.0_dp
!
!        if (pub_edft) then
!           ! kkbd: If we've done an EDFT run we might have a non-integer net spin
!           !       that also doesn't equal its original value!
!           do is=1,pub_num_spins
!              sw_s_occ(is,PUB_1K) = sparse_trace(sw_denskern%kern%m(is,PUB_1K),sw_rep%overlap)
!           end do
!        else
!           sw_s_occ = real(sw_rep%n_occ,kind=dp)
!        end if
!
!        call eigenstates_dos_gp(sw_eigen_en, sw_s_occ(:,PUB_1K), lcao_basis, ham_type,&
!             name=trim(pub_rootname)//"_swopt")!, weight=real(ngwf_basis%num,dp)/real(lcao_basis%num,dp))
!
!
!        ! JA: Not weighted by occupancy
!        call internal_pdos_sum(sw_eigen_en, pdos_weights, &
!             orbital_species, orbital_ion, orbital_l,orbital_order,&
!             sw_s_occ(:,pub_1k), lcao_basis, lcao_basis, &
!             ham_type, .true., .false., .false., trim(pub_rootname)//"_swopt")!,&
!
!        ! JA: Weighted by occupancy for calculation of band centres
!        call internal_pdos_sum(sw_eigen_en, pdos_weights, &
!             orbital_species, orbital_ion, orbital_l,orbital_order,&
!             sw_s_occ(:,pub_1k), lcao_basis, lcao_basis, &
!             ham_type, .true., print_qc,.true.,trim(pub_rootname)//"_swopt_occ",sw_eigen_occ)!,&
!
!        deallocate(sw_s_occ,stat=ierr)
!        call utils_dealloc_check('eigenstates_pdos','sw_s_occ',ierr)
!
!
!        if(pub_on_root) then
!           deallocate(pdos_weights,stat=ierr)
!           call utils_dealloc_check('eigenstates_pdos','pdos_weights',ierr)
!        end if
!
!        call dense_destroy(sw_lowdin)
!        call dense_destroy(sw_ham_dens)
!        call dense_destroy(sw_overlap_dens)
!        call sparse_destroy(sw_KS)
!        call sparse_destroy(sw_SKS)
!        call dense_destroy(sw_buffer_dens)
!
!        deallocate(sw_eigen_occ,stat=ierr)
!        call utils_dealloc_check('eigenstates_pdos','sw_eigen_occ',ierr)
!        deallocate(sw_eigen_en,stat=ierr)
!        call utils_dealloc_check('eigenstates_pdos','sw_eigen_en',ierr)
!
!        ! Destroy dense storage for eigenvectors
!        do is=1,pub_num_spins
!           call dense_destroy(sw_eigs_dens(is))
!        end do
!        deallocate(sw_eigs_dens,stat=ierr)
!        call utils_dealloc_check('eigenstates_pdos','sw_eigs_dens',ierr)
!
!        call data_functions_dealloc(sw_rep%ngwfs_on_grid)
!
!        deallocate(sw_lhxc_fine,stat=ierr)
!        call utils_dealloc_check('eigenstates_pdos','sw_lhxc_fine',ierr)
!
!
!        call kernel_destroy(sw_denskern)
!
!        call kernel_destroy(sw_pur_denskern)
!
!        if(pub_edft) then
!           call edft_destroy(sw_edft)
!        end if
!
!        call ngwf_ham_destroy(sw_ham)
!
!        call ngwf_rep_destroy(sw_rep)
!
!     end if lcao_optimize



    if(pub_on_root) then
       write(stdout,*)''
       write(stdout,'(a)', advance='no') 'Computing pDOS weights ...'
    end if


    call dense_convert(lcao_ngwf_overlap_DEM,lcao_ngwf_overlap)

    ! JA: allocate pdos weights on root
    if(pub_on_root) then
       if(pub_pdos_reduce_sws) then
          if(pub_pdos_construct_basis) then
             allocate(pdos_weights(sum(mdl%regions(:)%par%nat)*((pub_pdos_max_l+1)**2),num,pub_num_spins),stat=ierr)
             call utils_alloc_check('eigenstates_pdos','pdos_weights',ierr)
          else
             allocate(pdos_weights(lcao_num,num,pub_num_spins),stat=ierr)
             call utils_alloc_check('eigenstates_pdos','pdos_weights',ierr)
          end if
          pdos_weights=0.0_dp
       else if(pub_pdos_pseudoatomic) then
          allocate(pdos_weights(lcao_num,num,pub_num_spins),stat=ierr)
          call utils_alloc_check('eigenstates_pdos','pdos_weights',ierr)
          pdos_weights=0.0_dp
       else
          allocate(pdos_weights(lcao_num/num_n,num,pub_num_spins),stat=ierr)
          call utils_alloc_check('eigenstates_pdos','pdos_weights',ierr)
          pdos_weights=0.0_dp
       end if
    end if

    spil_elec = 0.0_dp
    spil_param = 0.0_dp
    spil_param_occ = 0.0_dp

    do spin=1,pub_num_spins

       ! agrecocmplx: weights will be real, but need complex matrix
       ! to actually make the product
       ! JA: Form the right half of the resolution of identity
       call dense_product(pdos_weights_DEM2,lcao_ngwf_overlap_DEM,eigenvecs_dens(spin))

       if(spin==1) then
          if(pub_pdos_lowdin.or.non_pos_def_overlap) then

             call dense_create(sw_S_evecs,lcao_num,lcao_num, &
                  iscmplx=loc_cmplx)

             call dense_scale(sw_S_evecs, 0.0_dp)

             allocate(sw_S_evals(lcao_num),stat=ierr)
             call utils_alloc_check('eigenstates_pdos','sw_S_evals',ierr)

             sw_S_evals=0.0_dp

             call dense_normal_eigensolve(lcao_num,sw_S_evals,lcao_overlap_dem,sw_S_evecs,&
                  arg_orfac=pub_eigensolver_orfac, arg_abstol=pub_eigensolver_abstol, rtype=1)

             if(pub_on_root) then
                open(newunit=fileunit, file="overlap_eigvals.txt")
                do i=1,lcao_num
                   write(fileunit,*) sw_S_evals(i)
                end do
                close(fileunit)
             end if

             call dense_scale(lcao_overlap_dem, 0.0_dp)

             ! JA: construct S^(-1/2) or S^(-1)
             if(pub_pdos_lowdin) then
                if(non_pos_def_overlap) then
                   do orbital_index= 1,lcao_num
                      if(sw_S_evals(orbital_index)>epsilon(0.0_dp)) then
                         sw_S_evals(orbital_index) = 1.0_dp/sqrt(sw_S_evals(orbital_index))
                      else
                         sw_S_evals(orbital_index) = 0.0_dp
                      end if
                   end do
                else
                   do orbital_index= 1,lcao_num
                      call utils_assert(sw_S_evals(orbital_index)>0.0_dp,&
                           'Error in eigenstates_pdos: SW overlap matrix is not positive definite! Index: ',orbital_index)
                      sw_S_evals(orbital_index) = sqrt(sw_S_evals(orbital_index))
                   end do
                end if
             else
                do orbital_index= 1,lcao_num
                   if(sw_S_evals(orbital_index)>epsilon(0.0_dp)) then
                      sw_S_evals(orbital_index) = 1.0_dp/sw_S_evals(orbital_index)
                   else
                      sw_S_evals(orbital_index) = 0.0_dp
                   end if
                end do
             end if

             call edft_matrix_from_eigs(lcao_overlap_dem, sw_S_evecs, sw_S_evals, lcao_num)

             call dense_destroy(sw_S_evecs)

             deallocate(sw_S_evals,stat=ierr)
             call utils_dealloc_check('eigenstates_pdos','sw_S_evals',ierr)

          end if
       end if

       if(spin==1) then

          ! JA: Print some QC test values.
          if (pub_print_qc) then
             call sparse_embed_trace(lcao_trace,lcao_overlap)
          end if

          if (pub_print_qc) then
             inv_lcao_trace=dense_trace(lcao_overlap_dem)
          end if

       end if


       ! JA: And the left half
       call dense_product(pdos_weights_DEM,lcao_overlap_dem,pdos_weights_DEM2)

       if(pub_pdos_lowdin) then

          call dense_copy(pdos_weights_DEM2, pdos_weights_DEM)

       end if


       ! JA: build projected density

!        if(pub_pdos_projected_density) then
!
!           if(pub_pdos_lowdin) then
!              call ngwf_rep_create(sw_rep,'p',mdl%elements,mdl%cell,loc_cmplx, par)
!
!              lowdin_factor%structure = "Kp"
!              call sparse_create(lowdin_factor)
!              call sparse_convert(lowdin_factor,lcao_overlap_dem)
!
!              lcao_overlap_block_diag%structure = "Dp"
!              call sparse_create(lcao_overlap_block_diag)
!              call sparse_scale(lcao_overlap_block_diag,0.0_dp)
!
!
!              call internal_projected_densities(lowdin_factor, lcao_ngwf_overlap, ngwf_lcao_overlap,&
!                   lcao_overlap_block_diag, ngwf_basis, ngwfs_on_grid, denskern)
!
!              call sparse_destroy(lcao_overlap_block_diag)
!              call sparse_destroy(lowdin_factor)
!              call ngwf_rep_destroy(sw_rep)
!           end if
!
!        end if


       ! JA: The final pdos_weights to be produced... probably only need this on the root proc, actually.
       if(pub_pdos_reduce_sws) then


          if(pub_pdos_construct_basis) then

             do icol=1,num
                do jrow=1,lcao_num
                   ! agrecocmplx: need to check if rval and rval2 are independently real
                   if (loc_cmplx) then
                      call dense_get_element(temp_coef%z,pdos_weights_DEM,jrow,icol)
                      rval = real(temp_coef%z,kind=DP)
                      call dense_get_element(temp_coef%z,pdos_weights_DEM2,jrow,icol)
                      rval2 = real(temp_coef%z,kind=DP)
                      ! real case
                   else
                      call dense_get_element(rval,pdos_weights_DEM,jrow,icol)
                      call dense_get_element(rval2,pdos_weights_DEM2,jrow,icol)
                   end if

                   do isub=1,mdl%nsub
                      do jsub=1,mdl%nsub
                         funcs_before_atom = &
                              (ngwf_basis(jsub)%first_on_atom(lcao_basis(isub)%atom_of_func(jrow))-1)*((pub_pdos_max_l+1)**2)

                         func_on_atom = jrow - funcs_before_atom

                         jrow2=(mod(func_on_atom-1,&
                              (pub_pdos_max_l+1)**2)+1+(lcao_basis(isub)%atom_of_func(jrow)-1)*((pub_pdos_max_l+1)**2))

                         if(pub_on_root) then
                            pdos_weights(jrow2,icol,spin)=pdos_weights(jrow2,icol,spin)+(rval*rval2)
                         end if
                      end do
                   end do
                end do
             end do


             if(pub_on_root) then
                if(pub_num_spins==1) then
                   open(newunit=matfile, file="pdos_weights.dat", recl=num*150)
                else
                   if(spin==UP) then
                      open(newunit=matfile, file="pdos_weights_UP.dat", recl=num*150)
                   elseif(spin==DN) then
                      open(newunit=matfile, file="pdos_weights_DOWN.dat", recl=num*150)
                   end if
                end if
                do icol=1,num
                   do isub=1,mdl%nsub
                      do atom=1,mdl%regions(isub)%par%nat
                         distr_atom=mdl%regions(isub)%par%distr_atom(atom)
                         do jrow=(((pub_pdos_max_l+1)**2)*(distr_atom-1))+1,((pub_pdos_max_l+1)**2)*distr_atom
                            write(matfile,'(f18.10,2x)',advance="NO") pdos_weights(jrow,icol,spin)
                            spil_param = spil_param + (pdos_weights(jrow,icol,spin)/(num*pub_num_spins))
                            spil_param_occ = spil_param_occ +(eigen_occ(icol,spin)*pdos_weights(jrow,icol,spin))
                         end do
                      end do
                   end do
                   write(matfile,*)
                   spil_elec = spil_elec + eigen_occ(icol,spin)
                end do
                close(matfile)
             end if


          else

             ! JA: Stick into pdos weights array
             ! agrecocmplx: need to differentiate between real and complex case
             ! need to take real part to get weights in the complex case
             do icol=1,num
                do jrow=1,lcao_num
                   ! agrecocmplx: rval and rval2 should be independently real, must check!
                   if (loc_cmplx) then
                      call dense_get_element(temp_coef%z,pdos_weights_DEM,jrow,icol)
                      rval = real(temp_coef%z,kind=DP)
                      call dense_get_element(temp_coef%z,pdos_weights_DEM2,jrow,icol)
                      rval2 = real(temp_coef%z,kind=DP)
                      ! real case
                   else
                      call dense_get_element(rval,pdos_weights_DEM,jrow,icol)
                      call dense_get_element(rval2,pdos_weights_DEM2,jrow,icol)
                   end if
                   if(pub_on_root) then
                      pdos_weights(jrow,icol,spin)=rval*rval2
                   end if
                end do
             end do

             if(pub_on_root) then
                if(pub_num_spins==1) then
                   open(newunit=matfile, file="pdos_weights.dat", recl=lcao_num*150)
                else
                   if(spin==UP) then
                      open(newunit=matfile, file="pdos_weights_UP.dat", recl=lcao_num*150)
                   elseif(spin==DN) then
                      open(newunit=matfile, file="pdos_weights_DOWN.dat", recl=lcao_num*150)
                   end if
                end if
                do icol=1,num

                   do jsub=1,mdl%nsub
                      do atom=1,mdl%regions(jsub)%par%nat
                         distr_atom=mdl%regions(jsub)%par%distr_atom(atom)
                         do isub=1,mdl%nsub
                            do jrow=lcao_basis(isub)%first_on_atom(distr_atom), &
                                 lcao_basis(isub)%num_on_atom(distr_atom)+lcao_basis(isub)%first_on_atom(distr_atom)-1

                               write(matfile,'(f18.10,2x)',advance="NO") pdos_weights(jrow,icol,spin)
                               spil_param = spil_param + (pdos_weights(jrow,icol,spin)/(num*pub_num_spins))
                               spil_param_occ = spil_param_occ +(eigen_occ(icol,spin)*pdos_weights(jrow,icol,spin))
                            end do
                         end do
                      end do
                   end do

                   write(matfile,*)
                   spil_elec = spil_elec + eigen_occ(icol,spin)
                end do
                close(matfile)
             end if
          end if

       else if(pub_pdos_pseudoatomic) then

          ! JA: Stick into pdos weights array
          ! agrecocmplx: need to differentiate between real and complex case
          ! need to take real part to get weights in the complex case
          do icol=1,num
             do jrow=1,lcao_num
                ! agrecocmplx: rval and rval2 should be independently real, must check!
                if (loc_cmplx) then
                   call dense_get_element(temp_coef%z,pdos_weights_DEM,jrow,icol)
                   rval = real(temp_coef%z,kind=DP)
                   call dense_get_element(temp_coef%z,pdos_weights_DEM2,jrow,icol)
                   rval2 = real(temp_coef%z,kind=DP)
                   ! real case
                else
                   call dense_get_element(rval,pdos_weights_DEM,jrow,icol)
                   call dense_get_element(rval2,pdos_weights_DEM2,jrow,icol)
                end if
                if(pub_on_root) then
                   pdos_weights(jrow,icol,spin)=rval*rval2
                end if
             end do
          end do

          if(pub_on_root) then
             if(pub_num_spins==1) then
                open(newunit=matfile, file="pdos_weights.dat", recl=lcao_num*150)
             else
                if(spin==UP) then
                   open(newunit=matfile, file="pdos_weights_UP.dat", recl=lcao_num*150)
                elseif(spin==DN) then
                   open(newunit=matfile, file="pdos_weights_DOWN.dat", recl=lcao_num*150)
                end if
             end if
             do icol=1,num
                do jsub=1,mdl%nsub
                   do atom=1,mdl%regions(jsub)%par%nat
                      distr_atom=mdl%regions(jsub)%par%distr_atom(atom)
                      do isub=1,mdl%nsub
                         do jrow=lcao_basis(isub)%first_on_atom(distr_atom), lcao_basis(isub)%num_on_atom(distr_atom)+&
                              lcao_basis(isub)%first_on_atom(distr_atom)-1

                            write(matfile,'(f18.10,2x)',advance="NO") pdos_weights(jrow,icol,spin)
                            spil_param = spil_param + (pdos_weights(jrow,icol,spin)/(num*pub_num_spins))
                            spil_param_occ = spil_param_occ +(eigen_occ(icol,spin)*pdos_weights(jrow,icol,spin))
                         end do
                      end do
                   end do
                end do
                write(matfile,*)
                spil_elec = spil_elec + eigen_occ(icol,spin)
             end do
             close(matfile)
          end if


       else

          do icol=1,num
             do jrow=1,lcao_num
             ! agrecocmplx: need to check if rval and rval2 are independently real
                if (loc_cmplx) then
                   call dense_get_element(temp_coef%z,pdos_weights_DEM,jrow,icol)
                   rval = real(temp_coef%z,kind=DP)
                   call dense_get_element(temp_coef%z,pdos_weights_DEM2,jrow,icol)
                   rval2 = real(temp_coef%z,kind=DP)
                ! real case
                else
                   call dense_get_element(rval,pdos_weights_DEM,jrow,icol)
                   call dense_get_element(rval2,pdos_weights_DEM2,jrow,icol)
                end if
                if(pub_on_root) then
                   jrow2=((mod(jrow-1,((pub_pdos_max_l+1)**2))+1)+&
                        &(((jrow-1)/(((pub_pdos_max_l+1)**2)*num_n))*((pub_pdos_max_l+1)**2)))
                   pdos_weights(jrow2,icol,spin)=pdos_weights(jrow2,icol,spin)+(rval*rval2)
                end if
             end do
          end do

          if(pub_on_root) then
             if(pub_num_spins==1) then
                open(newunit=matfile, file="pdos_weights.dat", recl=num*150)
             else
                if(spin==UP) then
                   open(newunit=matfile, file="pdos_weights_UP.dat", recl=num*150)
                elseif(spin==DN) then
                  open(newunit=matfile, file="pdos_weights_DOWN.dat", recl=num*150)
                end if
             end if
             do icol=1,num
                do jsub=1,mdl%nsub
                   do atom=1,mdl%regions(jsub)%par%nat
                      distr_atom=mdl%regions(jsub)%par%distr_atom(atom)
                      do isub=1,mdl%nsub
                         do jrow=((lcao_basis(isub)%first_on_atom(distr_atom)-1)/num_n)+1, &
                              ((lcao_basis(isub)%num_on_atom(distr_atom)+lcao_basis(isub)%first_on_atom(distr_atom)-2)/num_n)+1
                            write(matfile,'(f18.10,2x)',advance="NO") pdos_weights(jrow,icol,spin)
                            spil_param = spil_param + (pdos_weights(jrow,icol,spin)/(num*pub_num_spins))
                            spil_param_occ = spil_param_occ +(eigen_occ(icol,spin)*pdos_weights(jrow,icol,spin))
                         end do
                      end do
                   end do
                end do
                write(matfile,*)
                spil_elec = spil_elec + eigen_occ(icol,spin)
             end do
             close(matfile)
          end if

       end if

    end do


    if(pub_on_root) write(stdout,*) ' done'


    ! JA: Print some QC test values.
    if (pub_print_qc) then
       if (pub_on_root) then
          call utils_qc_print('lcao_overlap_trace',lcao_trace)
       end if
    end if

    if (pub_print_qc) then
       if (pub_on_root) then
          call utils_qc_print('inv_lcao_overlap_trace',inv_lcao_trace)
       end if
    end if

    !JA: we only need this on root & components exist solely on root anyway.
    if(pub_on_root) then
       spil_param = 1.0_dp - spil_param
       !JA: be extra careful not to divide by zero...
       if(abs(spil_elec)>epsilon(spil_elec)) then
          spil_param_occ = 1.0_dp - (spil_param_occ/spil_elec)
       else
          spil_param_occ = spil_param
       end if
    end if

    if(pub_on_root) then
       write(stdout,*)''
!    write(stdout,'(a)')'Computed Spilling Parameters '
       write(stdout,'(a31,f6.2,1x,a1)')'All bands spilling parameter = ', (100*spil_param),'%'
       write(stdout,'(a40,f6.2,1x,a1)')'Occupancy-weighted spilling parameter = ',(100*spil_param_occ),'%'
    endif


    call comms_barrier()


    if (pub_print_qc) then
       call sparse_embed_trace(rval,lcao_ngwf_overlap,ngwf_lcao_overlap)
       if (pub_on_root) then
          call utils_qc_print('lcao_ngwf_lcao_trace',rval)
       end if
    end if


    ! JA: Clean up
    ! JA: Don't need the SW overlap anymore...
    call dense_destroy(lcao_overlap_dem)


    !    call sparse_destroy(ngwf_overlap)
    call dense_destroy(pdos_weights_DEM)
    call dense_destroy(pdos_weights_DEM2)
    call dense_destroy(lcao_ngwf_overlap_DEM)


    call data_fftbox_dealloc(lcao_fftbox)
    ! agrecocmplx
    if (loc_cmplx) then
       deallocate(lcao_fftbox_real,stat=ierr)
       call utils_dealloc_check('eigenstates_pdos','lcao_fftbox_real',ierr)
    end if


    deallocate(lcao_work,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','lcao_work',ierr)

    ! jcap: loop over regions
    do isub=1,mdl%nsub
       call data_functions_dealloc(lcao_on_grid(isub))
    end do

    deallocate(elems_sfc,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','elems_sfc',ierr)
    deallocate(seg_nzb,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','seg_nzb',ierr)
    deallocate(seg_nze,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','seg_nze',ierr)


    deallocate(num_orbitals_on_procs,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','num_orbitals_on_procs',ierr)


    call sparse_embed_destroy(lcao_overlap)
    call sparse_embed_destroy(lcao_ngwf_overlap)
    call sparse_embed_destroy(invoverlap)

    call sparse_embed_destroy(ngwf_lcao_overlap)

    call comms_barrier

    ! JA: write out pdos weights
    if (pub_on_root.and.pub_pdos_optados_output) then
       write(stdout,*)''
       write(stdout,'(a)') ' => Outputting data for OptaDOS <='

       ! JA: Decide on a suffix
       if (present(ham_type)) then
          if (ham_type=='valence') then
             write(filename,*)trim(pub_rootname)//'.val_pdos_bin'
          else if (ham_type=='proj') then
             write(filename,*)trim(pub_rootname)//'.proj_pdos_bin'
          else if (ham_type=='cond') then
             write(filename,*)trim(pub_rootname)//'.cond_pdos_bin'
          else if (ham_type=='joint') then
             write(filename,*)trim(pub_rootname)//'.joint_pdos_bin'
          end if
       else
          write(filename,*)trim(pub_rootname)//'.pdos_bin'
       end if

       filename = adjustl(filename)


       write(stdout,*)''
       write(stdout,'(3a)', advance='no')&
            'Writing pDOS weights to file "',trim(adjustl(filename)),'" ...'

       fileunit = utils_unit()
#ifndef F2008
       ! jd: 'convert' is non-standard Fortran.
       open(unit=fileunit, file=trim(filename), form='unformatted', &
            action='write', iostat=ierr, status='replace', convert='big_endian')
#else
       call utils_feature_not_supported('eigenstates_pdos: '//&
            'OptaDos big-endian output has not been compiled in because it &
            &needs "convert", which is not part of the F2008 standard, &
            &and you asked for an F2008-compliant binary with -D2008.')
#endif
       ! JA: NB Convert is required to maintain compatibility.
       call utils_open_unit_check('eigenstates_pdos','output_unit',ierr)

       write(fileunit) file_version
       write(fileunit) file_header

       kpoint=1
       write(fileunit) kpoint ! JA: Single k-point hard coded...
       write(fileunit) pub_num_spins
       write(fileunit) num_orbitals
       write(fileunit) num

!       allocate(orbital_species(num_orbitals),stat=ierr)
!       call utils_alloc_check('eigenstates_pdos','orbital_species',ierr)
!       allocate(orbital_ion(num_orbitals),stat=ierr)
!       call utils_alloc_check('eigenstates_pdos','orbital_ion',ierr)
!       allocate(orbital_l(num_orbitals),stat=ierr)
!       call utils_alloc_check('eigenstates_pdos','orbital_l',ierr)
!       orbital_species=-1
!       orbital_ion=-1
!       orbital_l=-1

!       num_species=maxval((/(mdl%elements(atom)%species_number,atom=1,mdl%nat)/))
!       allocate(atom_in_species(mdl%nat),stat=ierr)
!       call utils_alloc_check('eigenstates_pdos','atom_in_species',ierr)
!       allocate(count_in_species(num_species),stat=ierr)
!       call utils_alloc_check('eigenstates_pdos','count_in_species',ierr)
!
!       atom_in_species=0
!       count_in_species=0
!       do atom=1,mdl%nat
!          count_in_species(mdl%elements(atom)%species_number)=count_in_species(mdl%elements(atom)%species_number)+1
!          atom_in_species(atom)=count_in_species(mdl%elements(atom)%species_number)
!       end do
!
!       if(pub_pdos_reduce_sws) then
!          do atom=1,mdl%nat
!             orbital_species((atom-1)*((pub_pdos_max_l+1)**2)+1:atom*((pub_pdos_max_l+1)**2))=mdl%elements(atom)%species_number
!             orbital_ion((atom-1)*((pub_pdos_max_l+1)**2)+1:atom*((pub_pdos_max_l+1)**2))=atom_in_species(atom)!atom
!             orbital_l((atom-1)*((pub_pdos_max_l+1)**2)+1:atom*((pub_pdos_max_l+1)**2)) = &
!                  (/((/(ll,l2=1,((ll*2)+1))/),ll=0,pub_pdos_max_l)/)
!
!          end do
!       else
!          do atom=1,mdl%nat
!             orbital_species((atom-1)*num_n*((pub_pdos_max_l+1)**2)+1:atom*num_n*((pub_pdos_max_l+1)**2))=&
!                  &mdl%elements(atom)%species_number
!             orbital_ion((atom-1)*num_n*((pub_pdos_max_l+1)**2)+1:atom*num_n*((pub_pdos_max_l+1)**2))=atom_in_species(atom)!atom
!             orbital_l((atom-1)*num_n*((pub_pdos_max_l+1)**2)+1:atom*num_n*((pub_pdos_max_l+1)**2)) = &
!                  (/((/((/(ll,l2=1,((ll*2)+1))/),ll=0,pub_pdos_max_l)/),nn=1,num_n)/)
!          end do
!
!       end if

       write(fileunit) orbital_species
       write(fileunit) orbital_ion
       write(fileunit) orbital_l


       kpoint=1
       write(fileunit) kpoint, 0.0_dp, 0.0_dp, 0.0_dp ! JA: first k-point of 1: Gamma Point.
       ! JA: NGWFs are not spin-polarised, so we just output twice.
       ! JA: Edit... but the eigenstates of the Hamiltonian are... so it's actually spin polarised
       if(pub_pdos_reduce_sws) then
          if(pub_pdos_construct_basis) then
             do spin = 1, pub_num_spins
                write(fileunit) spin
                write(fileunit) num
                do neigen = 1,num
                   do jsub=1,mdl%nsub
                      do atom=1,mdl%regions(jsub)%par%nat
                         distr_atom=mdl%regions(jsub)%par%distr_atom(atom)
                         do jrow=(((pub_pdos_max_l+1)**2)*(distr_atom-1))+1,((pub_pdos_max_l+1)**2)*distr_atom
                            write(fileunit) real(pdos_weights(jrow,neigen,spin),dp)
                         end do
                      end do
                   end do
!                   write(fileunit) real(pdos_weights(1:num_orbitals,neigen,spin),dp)
                end do
             end do

          else
             do spin = 1, pub_num_spins
                write(fileunit) spin
                write(fileunit) num
                do neigen = 1,num
                   do jsub=1,mdl%nsub
                      do atom=1,mdl%regions(jsub)%par%nat
                         distr_atom=mdl%regions(jsub)%par%distr_atom(atom)
                         do isub=1,mdl%nsub
                            do jrow=lcao_basis(isub)%first_on_atom(distr_atom), lcao_basis(isub)%num_on_atom(distr_atom)+&

                                 lcao_basis(isub)%first_on_atom(distr_atom)-1
                               write(fileunit) real(pdos_weights(jrow,neigen,spin),dp)
                            end do
                         end do
                      end do
                   end do
!                   write(fileunit) real(pdos_weights(1:num_orbitals,neigen,spin),dp)
                end do
             end do
          end if
       else if(pub_pdos_pseudoatomic) then
          do spin = 1, pub_num_spins
             write(fileunit) spin
             write(fileunit) num
             do neigen = 1,num
                do jsub=1,mdl%nsub
                   do atom=1,mdl%regions(jsub)%par%nat
                      distr_atom=mdl%regions(jsub)%par%distr_atom(atom)
                      do isub=1,mdl%nsub
                         do jrow=lcao_basis(isub)%first_on_atom(distr_atom), lcao_basis(isub)%num_on_atom(distr_atom)+&

                              lcao_basis(isub)%first_on_atom(distr_atom)-1
                            write(fileunit) real(pdos_weights(jrow,neigen,spin),dp)

                         end do
                      end do
                   end do
                end do
             end do
          end do
       else
          do spin = 1, pub_num_spins
             write(fileunit) spin

             write(fileunit) num
             do neigen = 1,num
                do jsub=1,mdl%nsub
                   do atom=1,mdl%regions(jsub)%par%nat
                      distr_atom=mdl%regions(jsub)%par%distr_atom(atom)
                      do isub=1,mdl%nsub
                         do jrow=((lcao_basis(isub)%first_on_atom(distr_atom)-1)/num_n)+1, &
                              ((lcao_basis(isub)%num_on_atom(distr_atom)+lcao_basis(isub)%first_on_atom(distr_atom)-2)/num_n)+1
                            write(fileunit) real(pdos_weights(jrow,neigen,spin),dp)
                         end do
                      end do
                   end do
                end do
                !                write(fileunit) real(pdos_weights(1:num_orbitals/num_n,neigen,spin),dp)
             end do
          end do
       end if

       close(unit=fileunit,iostat=ierr)
       call utils_close_unit_check('eigenstates_pdos','output_unit',ierr)
       write(stdout,'(a)') ' done'

       ! JA: Spoof band gradients for OptaDOS.
       if (present(ham_type)) then
          if (ham_type=='valence') then
             write(filename,*)trim(pub_rootname)//'.val_dome_bin'
          else if (ham_type=='proj') then
             write(filename,*)trim(pub_rootname)//'.proj_dome_bin'
          else if (ham_type=='cond') then
             write(filename,*)trim(pub_rootname)//'.cond_dome_bin'
          else if (ham_type=='joint') then
             write(filename,*)trim(pub_rootname)//'.joint_dome_bin'
          end if
       else
          write(filename,*)trim(pub_rootname)//'.dome_bin'
       end if
       filename = adjustl(filename)
       write(stdout,*)''
       write(stdout,'(3a)', advance='no')&
            'Writing band gradients to file "',trim(adjustl(filename)),'" ...'
       fileunit = utils_unit()

#ifndef F2008
       ! jd: 'convert' is non-standard Fortran.
       open(unit=fileunit, file=trim(filename), form='unformatted', &
            action='write', iostat=ierr, status='replace', convert='big_endian')
       ! JA: NB Convert is required to maintain compatibility.
#else
       call utils_abort('eigenstates_pdos: big-endian override not &
            &compiled in because "convert" is a GCC extension.')
#endif


       call utils_open_unit_check('eigenstates_pdos','output_unit',ierr)

       write(fileunit) file_version
       write(fileunit) file_header

       ! JA: ONETEP is Gamma point only... so we can safely assume that the band gradients are zero... right?
       do spin = 1, pub_num_spins
          write(fileunit) (0.0_dp,neigen=1,3*num)
          ! (((band_gradient_in_group(na,indx,1,spin)),na=1,num_eigenvalues(spin)),indx=1,3)
       end do

       close(unit=fileunit,iostat=ierr)
       call utils_close_unit_check('eigenstates_pdos','output_unit',ierr)
       write(stdout,'(a)') ' done'

    endif
    call comms_barrier

    if(pub_pdos_optados_output) then
       write(filename,*) trim(adjustl(pub_rootname))//'-out.cell'
       call model_write_castep_cell_file(filename,mdl)
    end if


    !JA: Made internal_pdos_sum called only on root - to avoid compiler
    !    complaining about non allocated arguments.
    if(pub_on_root) then

       ! JA: Not weighted by occupancy
       call internal_pdos_sum(eigen_en, pdos_weights, orbital_species, orbital_ion, orbital_l,&
            orbital_order, s_occ, ngwf_basis, &
            ham_type, .true., .false., .false., trim(pub_rootname))
       ! JA: Weighted by occupancy for calculation of band centres
       call internal_pdos_sum(eigen_en, pdos_weights, orbital_species, orbital_ion, orbital_l,&
            orbital_order, s_occ, ngwf_basis, &
            ham_type, .true., pub_print_qc, .true.,trim(pub_rootname)//"_occ",eigen_occ)


       deallocate(orbital_species,stat=ierr)
       call utils_dealloc_check('eigenstates_pdos','orbital_species',ierr)
       deallocate(orbital_ion,stat=ierr)
       call utils_dealloc_check('eigenstates_pdos','orbital_ion',ierr)
       deallocate(orbital_l,stat=ierr)
       call utils_dealloc_check('eigenstates_pdos','orbital_l',ierr)
       deallocate(orbital_order,stat=ierr)
       call utils_dealloc_check('eigenstates_pdos','orbital_order',ierr)
       deallocate(atom_in_species,stat=ierr)
       call utils_dealloc_check('eigenstates_pdos','atom_in_species',ierr)
       deallocate(count_in_species,stat=ierr)
       call utils_dealloc_check('eigenstates_pdos','count_in_species',ierr)


       deallocate(pdos_weights,stat=ierr)
       call utils_dealloc_check('eigenstates_pdos','pdos_weights',ierr)
    end if

    if(pub_on_root) write(stdout,'(a)') repeat('=',80)

    call comms_barrier

    ! JA: stop timer
    call timer_clock('eigenstates_pdos',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving eigenstates_pdos'

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_pdos_sum(eigen_en, pdos_weights, orbital_species,&
         orbital_ion, orbital_l, orbital_order, s_occ, ngwf_basis, &
         ham_type, write_out, write_qc,write_band_centres, &
         filename,eigen_occ, weight)
    !==================================================================!
    ! This subroutine sums over columns of the PDOS weights matrix to  !
    ! produce a histogram column matrix. It can optionally sum over    !
    ! all atoms, or groups of atoms, over all angular momentum states  !
    ! or some and over magnetic quantum numbers, or not. It then feeds !
    ! this matrix into eigenstates_histogram to produce the Gaussian   !
    ! weighted PDOS.                                                   !
    ! To be called only on root!                                       !
    !------------------------------------------------------------------!
    ! Written by Jolyon Aarons, 2019.                                  !
    !==================================================================!
      use constants, only: DP, stdout
      use function_basis, only: FUNC_BASIS
      use rundat, only: pub_pdos_groups,pub_pdos_group_nsp, pub_pdos_ngroups, &
           pub_pdos_sum_mag
      use utils, only: utils_alloc_check, utils_dealloc_check, &
           utils_qc_print

      implicit none

      real(kind=DP), intent(in) :: eigen_en(:,:)       ! molecular orbital energies
      real(kind=dp), intent(in) :: pdos_weights(:,:,:) ! weights from the ang mom projection
      integer,       intent(in) :: orbital_species(:)
      integer,       intent(in) :: orbital_ion(:)
      integer,       intent(in) :: orbital_l(:)
      integer,       intent(in) :: orbital_order(:)
      real(kind=dp), intent(in) :: s_occ(:)            ! spin-channel occupancy
      type(FUNC_BASIS), intent(in) :: ngwf_basis(:)       ! Function basis for NGWFs

      character(len=*), intent(in) :: ham_type
      logical,          intent(in) :: write_out
      logical,          intent(in) :: write_qc
      logical,          intent(in) :: write_band_centres
      character(len=*), intent(in) :: filename
      real(kind=DP), optional, intent(in) :: eigen_occ(:,:)     ! molecular orbital occupancies
      real(kind=dp), optional, intent(in) :: weight

      ! Local Variables
      real(kind=DP), allocatable, dimension(:,:,:) :: histo_weights ! histogram weight
      integer, allocatable, dimension(:,:) :: histo_groups
      logical, allocatable, dimension(:,:,:) :: ang_mom_in_group
      integer, allocatable, dimension(:,:,:) :: num_on_atom_l
      integer, allocatable, dimension(:,:,:) :: n_in_group
      integer, allocatable, dimension(:,:) :: ang_mom_of_histo
      integer, allocatable, dimension(:,:) :: ang_mom_of_pdos
      integer, allocatable, dimension(:,:) :: mag
      real(kind=DP), allocatable, dimension(:,:) :: occs
      integer, allocatable, dimension(:) :: nat_in_species
      integer, allocatable, dimension(:) :: orbital_number

      integer :: max_ang_mom


      integer :: is   ! spin loop counter
      integer :: ierr ! memory allocation error flag
      integer :: orbital_index ! orbital counting index
      integer :: mm

      integer :: ia1,ia2, iat, iat_orig
      character(len=4) :: nchar
      character(len=2) :: lchar

      logical :: iat_in_group
      character(len=4) :: iat_orig_id
      integer :: igroup, ingwf
      real(kind=DP) :: rval
      real(kind=dp) :: loc_weight
      integer :: num_histos
      integer :: nat

      if (pub_debug_on_root) &
           write(stdout,'(a)') 'DEBUG: Entering internal_pdos_sum (eigenstates_pdos).'


      loc_weight=1.0_dp
      if(present(weight)) then
         loc_weight=weight
      end if


      max_ang_mom=maxval(orbital_l)

      allocate(occs(sum(ngwf_basis(:)%num),pub_num_spins),stat=ierr)
      call utils_alloc_check('eigenstates_pdos','occs',ierr)

      occs=1.0_dp
      if(present(eigen_occ)) then
         occs=eigen_occ
      end if


      ! JA: n_in_group represents how many histograms there will be in each angular momentum channel (l) of each group
      !         at the moment, the only meaning of this is how many magnetic quantum number histograms there will be
      !         in each angular momentum channel (l) of each group.
      if(pub_pdos_ngroups>0) then
         allocate(n_in_group(pub_pdos_ngroups,max_ang_mom+1,pub_num_spins),stat=ierr)
         call utils_alloc_check('eigenstates_pdos','n_in_group',ierr)
      else
         allocate(n_in_group(1,max_ang_mom+1,pub_num_spins),stat=ierr)
         call utils_alloc_check('eigenstates_pdos','n_in_group',ierr)
      end if
      n_in_group=0

      ! JA: work out how many atoms there are in each species
      allocate(nat_in_species(maxval(orbital_species)),stat=ierr)
      call utils_alloc_check('eigenstates_pdos','nat_in_species',ierr)
      nat_in_species=0
      do orbital_index=1,size(orbital_order)
         if(orbital_ion(orbital_index)>nat_in_species(orbital_species(orbital_index))) then
            nat_in_species(orbital_species(orbital_index))=orbital_ion(orbital_index)
         end if
      end do

      ! JA: work out how many atoms there are per orbital
      allocate(orbital_number(size(orbital_species)),stat=ierr)
      call utils_alloc_check('eigenstates_pdos','orbital_number',ierr)
      orbital_number=0
      do orbital_index=1,size(orbital_order)
         orbital_number(orbital_index)=sum(nat_in_species(1:orbital_species(orbital_index)-1))+orbital_ion(orbital_index)
      end do
      deallocate(nat_in_species,stat=ierr)
      call utils_dealloc_check('eigenstates_pdos','nat_in_species',ierr)


      ! total number of atoms
      nat=maxval(orbital_number)


      allocate(num_on_atom_l(maxval(orbital_l)+1,nat,pub_num_spins),stat=ierr)
      call utils_alloc_check('eigenstates_pdos','num_on_atom_l',ierr)

      allocate(mag(size(orbital_order),pub_num_spins),stat=ierr)
      call utils_alloc_check('eigenstates_pdos','mag',ierr)

      ! JA: num_on_atom_l gives the number of magnetic quantum numbers (m) in each angular momentum channel (l) of each atom
      !         mag is the index of each of these magnetic quantum numbers in each orbital
      !         setting them all to 1 will sum over m.
      num_on_atom_l=0
      mag=0
      if(pub_pdos_sum_mag) then
         num_on_atom_l=1
         mag=1
      else
         do is=1,pub_num_spins
            do orbital_index=1,size(orbital_order)
               ia1=orbital_l(orbital_index)+1
               iat_orig_id = mdl%species(orbital_species(orbital_index))%species_id

               ia2=orbital_number(orbital_index)

               num_on_atom_l(ia1,ia2,is)=num_on_atom_l(ia1,ia2,is)+1
               mag(orbital_index,is)=num_on_atom_l(ia1,ia2,is)
            end do
         end do
      end if

      do is=1,pub_num_spins
         do orbital_index=1,size(orbital_order)
            iat_orig_id = mdl%species(orbital_species(orbital_index))%species_id
            ia2=orbital_number(orbital_index)
            ia1=orbital_l(orbital_index)+1
            if(pub_pdos_ngroups>0) then
               do igroup=1,pub_pdos_ngroups
                  if(any(pub_pdos_groups(1:pub_pdos_group_nsp(igroup),igroup)==iat_orig_id)) then
                     n_in_group(igroup,ia1,is)=max(num_on_atom_l(ia1,ia2,is),n_in_group(igroup,ia1,is))
                  end if
               end do
            else
               n_in_group(1,ia1,is)=max(num_on_atom_l(ia1,ia2,is),n_in_group(1,ia1,is))
            end if
         end do
      end do
      num_histos=sum(n_in_group)/pub_num_spins


      deallocate(num_on_atom_l,stat=ierr)
      call utils_dealloc_check('eigenstates_pdos','num_on_atom_l',ierr)


      allocate(ang_mom_of_histo(num_histos,pub_num_spins),stat=ierr)
      call utils_alloc_check('eigenstates_pdos','ang_mom_of_histo',ierr)

      ! ang_mom_of_histo give the angular momentum (l) of each histogram to be made

      do is=1,pub_num_spins
         ia2=0
         do igroup=1,max(pub_pdos_ngroups,1)
            do ia1=1,max_ang_mom+1
               ang_mom_of_histo(ia2+1:ia2+n_in_group(igroup,ia1,is),is) = ia1-1
               ia2=ia2+n_in_group(igroup,ia1,is)
            end do
         end do
      end do

      allocate(histo_weights(num_histos,size(pdos_weights,2),size(pdos_weights,3)),stat=ierr)
      call utils_alloc_check('eigenstates_pdos','histo_weights',ierr)

      allocate(histo_groups(num_histos,pub_num_spins),stat=ierr)
      call utils_alloc_check('eigenstates_pdos','histo_groups',ierr)
      histo_groups=-1

      ! now we sum over the atoms in the groups and optionally over m.

      histo_weights=0.0_dp
      do is=1,pub_num_spins
         do orbital_index=1,size(orbital_order)
            ia1=orbital_l(orbital_index)+1
            ia2=orbital_order(orbital_index)
            iat_orig_id = mdl%species(orbital_species(orbital_index))%species_id
            mm=mag(orbital_index,is)
            !                  orbital_l orbital_ion orbital_species

            iat_in_group = .false.
            if(pub_pdos_ngroups>0) then
               do igroup=1,pub_pdos_ngroups
                  if(any(pub_pdos_groups(1:pub_pdos_group_nsp(igroup),igroup)==iat_orig_id)) then
                     iat_in_group=.true.
                  end if
               end do
            else
               iat_in_group=.true.
            end if
            if (.not.iat_in_group) cycle

            if(pub_pdos_ngroups>0) then
               do igroup=1,pub_pdos_ngroups
                  if(any(pub_pdos_groups(1:pub_pdos_group_nsp(igroup),igroup)==iat_orig_id)) then
                     do ingwf=1,sum(ngwf_basis(:)%num)

                        histo_weights(sum(n_in_group(1:igroup-1,:,is))+sum(n_in_group(igroup,1:ia1-1,is))+mm,ingwf,is) = &
                             & histo_weights(sum(n_in_group(1:igroup-1,:,is))+sum(n_in_group(igroup,1:ia1-1,is))+mm,ingwf,is) +&
                             & pdos_weights(ia2,ingwf,is)
                     end do
                     histo_groups(sum(n_in_group(1:igroup-1,:,is))+sum(n_in_group(igroup,1:ia1-1,is))+mm,is)=igroup
                  end if
               end do
            else
               do ingwf=1,sum(ngwf_basis(:)%num)
                  igroup=1
                  histo_weights(sum(n_in_group(1:igroup-1,:,is))+sum(n_in_group(igroup,1:ia1-1,is))+mm,ingwf,is) = &
                       & histo_weights(sum(n_in_group(1:igroup-1,:,is))+sum(n_in_group(igroup,1:ia1-1,is))+mm,ingwf,is) +&
                       & pdos_weights(ia2,ingwf,is)
               end do
               histo_groups(sum(n_in_group(1:igroup-1,:,is))+sum(n_in_group(igroup,1:ia1-1,is))+mm,is)=1
            end if
         end do
      end do


      deallocate(orbital_number,stat=ierr)
      call utils_dealloc_check('eigenstates_pdos','orbital_number',ierr)

      deallocate(mag,stat=ierr)
      call utils_dealloc_check('eigenstates_pdos','mag',ierr)

      deallocate(n_in_group,stat=ierr)
      call utils_dealloc_check('eigenstates_pdos','n_in_group',ierr)


      if (write_qc.and.write_out) then
         do is=1,pub_num_spins
            do ia1=1,size(histo_weights,2)
               do ia2=1,size(histo_weights,1)
                  write(lchar,'(i2)') ia2-1
                  write(nchar,'(i4)') ia1
                  call utils_qc_print('pdoswts(l='//trim(adjustl(lchar))//',eig_i='//trim(adjustl(nchar))//')', &
                       histo_weights(ia2,ia1,is))
               end do
            end do
         end do
      end if


      call eigenstates_histogram(eigen_en, occs, present(eigen_occ), s_occ, ngwf_basis, num_histos, &
           ang_mom_of_histo, histo_weights, histo_groups, loc_weight, ham_type, filename, write_out, &
           write_band_centres, write_qc)

      deallocate(histo_groups,stat=ierr)
      call utils_dealloc_check('eigenstates_pdos','histo_groups',ierr)

      deallocate(occs,stat=ierr)
      call utils_dealloc_check('eigenstates_pdos','occs',ierr)

      deallocate(histo_weights,stat=ierr)
      call utils_dealloc_check('eigenstates_pdos','histo_weights',ierr)

      deallocate(ang_mom_of_histo,stat=ierr)
      call utils_dealloc_check('eigenstates_pdos','ang_mom_of_histo',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving internal_pdos_sum (eigenstates_pdos).'

      return

    end subroutine internal_pdos_sum


  end subroutine eigenstates_pdos


  subroutine eigenstates_histogram(eigen_en, occs, occ_weighted, s_occ, ngwf_basis, num_histos, &
       ang_mom_of_histo, histo_weights, histo_groups, gaussian_weight, ham_type, filename, write_out, &
       write_band_centres, write_qc)
    !==================================================================!
    ! This subroutine prepares and outputs in simple txt format a      !
    ! projected density of states (pDOS) plot which has been generated !
    ! with the Gaussian smearing method using gamma point only         !
    ! energies. The pDOS is calculated as a histogram of the sum of    !
    ! all Gaussian-broadened energies. Energies in the output are in   !
    ! eV and so is the half-width of the smearing Gaussians which is   !
    ! set by the input parameter pdos_smear.                           !
    !------------------------------------------------------------------!
    ! Adapted from DOS version by Jolyon Aarons on 25/10/2014.         !
    ! Which was:                                                       !
    ! Written by Chris-Kriton Skylaris on 10/05/2006.                  !
    ! Moved to eigenstates_mod by Nicholas Hine in November 2011.      !
    !==================================================================!
    use constants, only: hartree_in_evs, stdout, UP, DN, pi
    use comms, only: pub_on_root
    use rundat, only: pub_num_spins, pub_cond_calculate,pub_dos_smear, &
         pub_pdos_ngroups, pub_pdos_max_l, pub_spin_fac
    use utils, only: utils_alloc_check, utils_unit, utils_open_unit_check,&
         utils_dealloc_check, utils_qc_print, utils_close_unit_check, utils_abort
    implicit none
    real(kind=DP), intent(in) :: eigen_en(:,:)       ! molecular orbital energies
    real(kind=DP), intent(in) :: occs(:,:)     ! molecular orbital occupancies
    logical,       intent(in) :: occ_weighted
    real(kind=dp), intent(in) :: s_occ(:)            ! spin-channel occupancy
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)       ! Function basis for NGWFs
    integer,       intent(in) :: num_histos
    integer,       intent(in) :: ang_mom_of_histo(:,:)
    real(kind=dp), intent(in) :: histo_weights(:,:,:)
    integer,       intent(in) :: histo_groups(:,:)
    real(kind=dp), intent(in) :: gaussian_weight
    character(len=*), intent(in) :: ham_type
    character(len=*), intent(in) :: filename
    logical,          intent(in) :: write_out
    logical,          intent(in) :: write_band_centres
    logical,          intent(in) :: write_qc

    real(kind=DP) :: alpha     ! Gaussian exponent
    real(kind=DP) :: delta_e   ! energy increment
    real(kind=DP) :: e_point   ! energy point

    real(kind=DP) :: dist_sq   ! squared distance
    real(kind=DP) :: gnorm     ! Gaussian normalisation factor

    character(len=256) :: output_file  ! output file name buffer
    real(kind=DP), allocatable, dimension(:,:) :: energies ! buffer for energies in eV

    real(kind=DP), allocatable, dimension(:) :: histo_val ! histogram value
    real(kind=DP), allocatable, dimension(:) :: histo_val_old ! histogram value
    real(kind=DP), allocatable, dimension(:) :: histo_val_oldold ! histogram value
    real(kind=DP), allocatable, dimension(:) :: histo_pro ! histogram value*x
    real(kind=DP), allocatable, dimension(:) :: histo_pro_old ! histogram value*x
    real(kind=DP), allocatable, dimension(:) :: histo_pro_oldold ! histogram value*x
    real(kind=dp),  allocatable, dimension(:) :: integral_num, integral_den

    character(len=:), allocatable :: histo_format

    integer :: output_unit ! fortran output unit number

    integer, parameter :: histonum=2001 ! number of histogram points
    integer :: row  ! row counting index
    real(kind=DP), parameter :: en_offset =3.0_DP ! left and right offset for energy scale

    integer :: io_status ! file access error flag
    integer :: igroup, ia1, ia2
    integer :: ierr
    integer :: is
    integer :: orbital_index
    integer :: recl


    ! JA: Histogram construction from this point onwards Vvvvvv


    allocate(energies(sum(ngwf_basis(:)%num),pub_num_spins),stat=ierr)
    call utils_alloc_check('eigenstates_pdos','energies',ierr)


    allocate(histo_val(num_histos),stat=ierr)
    call utils_alloc_check('eigenstates_pdos','histo_val',ierr)
    allocate(histo_val_old(num_histos),stat=ierr)
    call utils_alloc_check('eigenstates_pdos','histo_val_old',ierr)
    allocate(histo_val_oldold(num_histos),stat=ierr)
    call utils_alloc_check('eigenstates_pdos','histo_val_oldold',ierr)
    allocate(histo_pro(num_histos),stat=ierr)
    call utils_alloc_check('eigenstates_pdos','histo_pro',ierr)
    allocate(histo_pro_old(num_histos),stat=ierr)
    call utils_alloc_check('eigenstates_pdos','histo_pro_old',ierr)
    allocate(histo_pro_oldold(num_histos),stat=ierr)
    call utils_alloc_check('eigenstates_pdos','histo_pro_oldold',ierr)
    allocate(integral_num(num_histos),stat=ierr)
    call utils_alloc_check('eigenstates_pdos','integral_num',ierr)
    allocate(integral_den(num_histos),stat=ierr)
    call utils_alloc_check('eigenstates_pdos','integral_den',ierr)

    ! cks: store eigen energies in buffer and convert to eV
    energies = eigen_en*HARTREE_IN_EVS

    ! cks: smearing Gaussian exponent
    ! cks: pub_pdos_smear is the half-width of the Gaussian
    alpha =log(2.0_DP) /((pub_dos_smear*HARTREE_IN_EVS)**2)

    ! cks: normalisation factor for Gaussians
    gnorm = gaussian_weight*sqrt(alpha/PI)




    ! cks: generate DOS for each spin
    spin_loop: do is=1,pub_num_spins


       if(write_out) then
          if (pub_on_root) write(stdout,'(a)')' '

          write(output_file,*)trim(filename)//'_'
          if (pub_cond_calculate) &
               write(output_file,*)trim(output_file)//trim(ham_type)//'_'
          if ((pub_num_spins==2).and.(is == UP)) &
               write(output_file,*) trim(output_file)//'UP'//'_'
          if ((pub_num_spins==2).and.(is == DN)) &
               write(output_file,*) trim(output_file)//'DN'//'_'
          write(output_file,*) trim(output_file)//'PDOS.txt'

          if(occ_weighted) then
             if (pub_num_spins == 1) then
                write(stdout,'(a)') ' => Computing Occupancy-weighted Gaussian smeared pDOS <= '
             elseif (is == UP) then
                write(stdout,'(a)') ' => Computing Occupancy-weighted Gaussian smeared pDOS for UP spin<= '
             elseif (is == DN) then
                write(stdout,'(a)') ' => Computing Occupancy-weighted Gaussian smeared pDOS for DOWN spin<= '
             endif
          else
             if (pub_num_spins == 1) then
                write(stdout,'(a)') ' => Computing Gaussian smeared pDOS <= '
             elseif (is == UP) then
                write(stdout,'(a)') ' => Computing Gaussian smeared pDOS for UP spin<= '
             elseif (is == DN) then
                write(stdout,'(a)') ' => Computing Gaussian smeared pDOS for DOWN spin<= '
             endif
          end if

          output_file = adjustl(output_file)

          ! cks: print output warning
          write(stdout,'(3a)',advance ='no') &
               'Writing "', trim(output_file),'" ...'

          ! cks: get a unit number that is free
          output_unit = utils_unit()

          recl=130 + max(pub_pdos_ngroups,1)*(pub_pdos_max_l+1)*25*num_histos

          open(unit=output_unit, form="formatted" ,file=trim(output_file), &
               action="write",iostat=io_status,recl=recl)
          call utils_open_unit_check('eigenstates_pdos','output_file', &
               io_status)

          ! cks: write first line
          if(pub_pdos_ngroups>0) then

             write(output_unit,'(a)',err =100)'#                   pDOS (states/eV)'
             write(output_unit,'(a)',advance='no',err =100)'#  Energy (eV)   | '

             do ia1=1,num_histos

                igroup = histo_groups(ia1,is)

                if(ang_mom_of_histo(ia1,is)+1==1) then
                   write(output_unit,'(a7,i3.1,a6)',advance='no',err =100) ' group ',igroup,': S   '
                else if(ang_mom_of_histo(ia1,is)+1==2) then
                   write(output_unit,'(a7,i3.1,a6)',advance='no',err =100) ' group ',igroup,': P   '
                else if(ang_mom_of_histo(ia1,is)+1==3) then
                   write(output_unit,'(a7,i3.1,a6)',advance='no',err =100) ' group ',igroup,': D   '
                else if(ang_mom_of_histo(ia1,is)+1==4) then
                   write(output_unit,'(a7,i3.1,a6)',advance='no',err =100) ' group ',igroup,': F   '
                else
                   write(output_unit,'(a7,i3.1,a6)',advance='no',err =100) ' group ',igroup,&
                        ': '//achar(66+ang_mom_of_histo(ia1,is)+1)//'   '
                end if

             end do
             write(output_unit,'(a16)') '      SUM       '

          else

             write(output_unit,'(a)',err =100)'#                   pDOS (states/eV)'
             write(output_unit,'(a)',advance='no',err =100)'#  Energy (eV)   | '

             do ia1=1,size(integral_num)

                if(ang_mom_of_histo(ia1,is)+1==1) then
                   write(output_unit,'(a16)',advance='no',err =100) '       S        '
                else if(ang_mom_of_histo(ia1,is)+1==2) then
                   write(output_unit,'(a16)',advance='no',err =100) '       P        '
                else if(ang_mom_of_histo(ia1,is)+1==3) then
                   write(output_unit,'(a16)',advance='no',err =100) '       D        '
                else if(ang_mom_of_histo(ia1,is)+1==4) then
                   write(output_unit,'(a16)',advance='no',err =100) '       F        '
                else
                   write(output_unit,'(a16)',advance='no',err =100) &
                        '       '//achar(66+ang_mom_of_histo(ia1,is)+1)//'        '
                end if

             end do
             write(output_unit,'(a16)') '      SUM       '

          end if

       end if

       delta_e =(2.0_DP*en_offset + maxval(energies(ngwf_basis%num,:)) &
            - minval(energies(1,:)))/real(histonum-1, kind=DP)

       e_point = minval(energies(1,:)) - en_offset

       ! cks: Loop over DOS histogram points
       histo_val_old=0.0_dp
       histo_val_oldold=0.0_dp
       histo_pro_old=0.0_dp
       histo_pro_oldold=0.0_dp
       integral_num=0.0_dp
       integral_den=0.0_dp

       histo_format=repeat('F16.8,',size(histo_val)+1)//'F16.8'

       do row=1, histonum

          histo_val = 0.0_DP
          ! cks: accumulate value at current histogram point

          do orbital_index=1,sum(ngwf_basis(:)%num)
             dist_sq =(e_point -energies(orbital_index, is))**2
             histo_val = histo_val + occs(orbital_index,is)*histo_weights(:,orbital_index,is)*gnorm*exp(-alpha *dist_sq)
          end do
          histo_pro = histo_val*e_point

          ! JA: Integrals for band centres : numerator -> occ weighted band centres, denom -> integral of band
          if(row >= 3 .and. mod(row,2)>0) then
             integral_num = integral_num + (2.0_dp*delta_e / 6.0_dp) * (histo_pro_oldold + histo_pro + 4.0_dp*histo_pro_old)
             integral_den = integral_den + (2.0_dp*delta_e / 6.0_dp) * (histo_val_oldold + histo_val + 4.0_dp*histo_val_old)
          end if
          histo_val_oldold=histo_val_old
          histo_val_old=histo_val
          histo_pro_oldold=histo_pro_old
          histo_pro_old=histo_pro

          ! cks: write to file current histo-point

          if(write_out) write(output_unit,'('//histo_format//')',err=100) e_point, histo_val, sum(histo_val)

          ! cks: update next histogram coordinate
          e_point = e_point +delta_e
       end do

       if(write_out) then
          close(unit=output_unit,iostat=io_status)
          call utils_close_unit_check('eigenstates_pdos','output_unit', &
               io_status)
          write(stdout,*)' done'
       end if
       ! cks: notify of end of output

       ! JA:! end of histogram ^^^^^^


       ! JA: write band centres (Energy weighted integrals of DOS)
       if(write_band_centres) then
          if(pub_num_spins==1) then
             write(stdout,*) ' => Band centres: '
          else
             if(is==UP) then
                write(stdout,*) ' => Band centres for UP spin channel: '
             else if(is==DN) then
                write(stdout,*) ' => Band centres for DOWN spin channel: '
             end if
          end if
          !               ia2=1
          do ia1=1,size(integral_num)
             if (write_qc.and.write_out) then
                call utils_qc_print('band_centres',integral_num(ia1)/integral_den(ia1))
             end if

             igroup = histo_groups(ia1,is)

             ! JA: Could run into problems if we reach the Qth orbital(/s)!
             if(ang_mom_of_histo(ia1,is)+1==1) then
                write(stdout,'(a,i0,a,f12.6,a)') " S band centre of group ",igroup,":", &
                     integral_num(ia1)/integral_den(ia1)," eV"
             else if(ang_mom_of_histo(ia1,is)+1==2) then
                write(stdout,'(a,i0,a,f12.6,a)') " P band centre of group ",igroup,":", &
                     integral_num(ia1)/integral_den(ia1)," eV"
             else if(ang_mom_of_histo(ia1,is)+1==3) then
                write(stdout,'(a,i0,a,f12.6,a)') " D band centre of group ",igroup,":", &
                     integral_num(ia1)/integral_den(ia1)," eV"
             else if(ang_mom_of_histo(ia1,is)+1==4) then
                write(stdout,'(a,i0,a,f12.6,a)') " F band centre of group ",igroup,":", &
                     integral_num(ia1)/integral_den(ia1),"eV "
             else
                write(stdout,'(a,i0,a,f12.6,a)') " "//achar(66+ang_mom_of_histo(ia1,is)+1)//" band centre of group ",&
                     igroup,":", integral_num(ia1)/integral_den(ia1)," eV"
             end if

          end do
          if(pub_num_spins==1) then
             write(stdout,*) ' Band centres done. <='
          else
             if(is==UP) then
                write(stdout,*) ' Band centres for UP spin channel done. <= '
             else if(is==DN) then
                write(stdout,*) ' Band centres for DOWN spin channel done. <='
             end if
          end if

          ! JA: write integrated number of electrons in each band.
          if(pub_num_spins==1) then
             write(stdout,*) ' => Integrated number of electrons in each AM band: '
          else
             if(is==UP) then
                write(stdout,*) ' => Integrated number of electrons in each AM band for UP spin channel: '
             else if(is==DN) then
                write(stdout,*) ' => Integrated number of electrons in each AM band for DOWN spin channel: '
             end if
          end if

          do ia1=1,size(integral_num)

             igroup = histo_groups(ia1,is)

             ! JA: Could run into problems if we reach the Qth orbital(/s)!
             if(ang_mom_of_histo(ia1,is)+1==1) then
                write(stdout,'(a,i0,a,f12.6,a)') " S num electrons of group ",igroup,":", &
                     integral_den(ia1)*pub_spin_fac
             else if(ang_mom_of_histo(ia1,is)+1==2) then
                write(stdout,'(a,i0,a,f12.6,a)') " P num electrons of group ",igroup,":", &
                     integral_den(ia1)*pub_spin_fac
             else if(ang_mom_of_histo(ia1,is)+1==3) then
                write(stdout,'(a,i0,a,f12.6,a)') " D num electrons of group ",igroup,":", &
                     integral_den(ia1)*pub_spin_fac
             else if(ang_mom_of_histo(ia1,is)+1==4) then
                write(stdout,'(a,i0,a,f12.6,a)') " F num electrons of group ",igroup,":", &
                     integral_den(ia1)*pub_spin_fac
             else
                write(stdout,'(a,i0,a,f12.6,a)') " "//achar(66+ang_mom_of_histo(ia1,is)+1)//" num electrons of group ",igroup,":", &
                     integral_den(ia1)*pub_spin_fac
             end if

          end do
          if(pub_num_spins==1) then
             write(stdout,*) ' Integrated number of electrons done. <='
          else
             if(is==UP) then
                write(stdout,*) ' Integrated number of electrons for UP spin channel done. <= '
             else if(is==DN) then
                write(stdout,*) ' Integrated number of electrons for DOWN spin channel done. <='
             end if
          end if
       end if


    enddo spin_loop

    deallocate(energies,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','energies',ierr)


    deallocate(histo_val,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','histo_val',ierr)
    deallocate(histo_val_old,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','histo_val_old',ierr)
    deallocate(histo_val_oldold,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','histo_val_oldold',ierr)
    deallocate(histo_pro,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','histo_pro',ierr)
    deallocate(histo_pro_old,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','histo_pro_old',ierr)
    deallocate(histo_pro_oldold,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','histo_pro_oldold',ierr)
    deallocate(integral_num,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','integral_num',ierr)
    deallocate(integral_den,stat=ierr)
    call utils_dealloc_check('eigenstates_pdos','integral_den',ierr)

    return

100 call utils_abort('Problem writing to file in eigenstates_pdos.')

  end subroutine eigenstates_histogram




  !   subroutine internal_projected_densities(lowdin_factor, lcao_ngwf_overlap)
  !     implicit none
  !     type(SPAM3), intent(in) :: lowdin_factor
  !     type(SPAM3), intent(in) :: lcao_ngwf_overlap
  !
  !     call sparse_create(DP, denskern%kern%m(1,pub_1k))
  !     call sparse_create(PDP, denskern%kern%m(1,pub_1k)
  !
  !
  !     do igroup=1, num_groups
  !
  !        call sparse_scale(diag_matrix,1.0_dp)
  !
  !        do istate=first_elem_on_proc(pub_my_proc_id), first_elem_on_proc(pub_my_proc_id+1)-1
  !
  !           ang_mom =
  !
  !           atom =
  !
  !           atom_in_group=.false.
  !           atom_index_in_group=-1
  !           do i=1,num_in_group(igroup)
  !              if(atom_in_group(i,igroup)==atom) then
  !                 atom_in_group=.true.
  !                 atom_index_in_group=i
  !              end if
  !           end do
  !
  !           if(.not.atom_in_group) then
  !              call sparse_put_element(0.0_dp,diag_matrix,istate,istate)
  !           end if
  !
  !        end do
  !
  !        call sparse_product(tmp_square, diag_matrix, lowdin_factor
  !
  !        call sparse_product(tmp_rect, tmp_square, lcao_ngwf_overlap)
  !
  !        call sparse_product(tmp_rect2, lowdin_factor, lcao_ngwf_overlap)
  !
  !        call sparse_product(projection_mat, ngwf_lcao_overlap, tmp_rect2)
  !
  !     call sparse_product(DP, denskern, projection_mat)
  !     call sparse_product(PDP, projection_mat, DP)
  !
  !
  !        call function_ops_sum_ppd_funcs(proj_ngwfs_on_grid(1:1), &
  !             ngwf_basis, [projection_mat], 1, 1, overlap, &
  !             ngwfs_on_grid, ngwf_basis)
  !
  !        call density_on_grid(density_fine, &
  !             mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, denskern, &
  !             rep%ngwf_overlap, proj_ngwfs_on_grid, ngwf_basis, proj_ngwfs_on_grid, &
  !             ngwf_basis)
  !
  !     end do
  !
  !     call sparse_destroy(DP)
  !     call sparse_destroy(PDP)
  !
  !   end subroutine internal_projected_densities

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_lumo_search(rep, ham, denskern, ngwf_basis, mdl, &
       hfxstate)

    !========================================================================!
    ! This subroutine                                                        !
    !------------------------------------------------------------------------!
    !  denskern  (inout) : The ground state density kernel                   !
    !  rep          (in) : Valence NGWF representation                       !
    !  ham          (in) : Valence NGWF Hamiltonian                          !
    !  ngwf_basis   (in) : Function basis type for the valence NGWFs         !
    !------------------------------------------------------------------------!
    ! Written by Nicholas Hine in February 2016.                             !
    ! Modified for embedding structures by Robert Charlton, October 2018.    !
    !========================================================================!

    use datatypes, only: FUNCTIONS, data_functions_alloc, data_functions_dealloc
    use comms, only: pub_on_root
    use constants, only: stdout, VERBOSE
    use dense, only: DEM
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE
    use kernel, only: DKERN, kernel_bandgap, &
         kernel_init_core_ham
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM, ngwf_ham_create, &
         ngwf_ham_destroy
    use palser_mano, only: palser_mano_kernel_optimise
    use population, only: population_analysis_mulliken
    use restart, only: restart_hamiltonian_read, restart_overlap_read
    use rundat, only: pub_output_detail, pub_debug_on_root, &
         pub_num_spins, pub_max_resid_hotelling, pub_maxit_hotelling, &
         pub_maxit_palser_mano, pub_num_kpoints, pub_popn_mulliken_partial, &
         pub_num_eigenvalues, pub_read_hamiltonian, PUB_1K, &
         pub_devel_code
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_copy, &
         sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_hotelling_init, sparse_embed_hotelling_invert, &
         sparse_embed_product, sparse_embed_outer_product, sparse_embed_axpy, &
         sparse_embed_extremal_eigenvalue, sparse_embed_trace
    use utils, only: utils_abort, utils_assert, utils_alloc_check, &
         utils_dealloc_check, utils_isnan, utils_devel_code, &
         utils_banner

    implicit none

    ! Arguments
    type(MODEL), intent(inout) :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    type(FUNC_BASIS), intent(in) :: ngwf_basis(mdl%nsub)
    type(DKERN), intent(inout)   :: denskern
    type(NGWF_REP), intent(inout)   :: rep
    type(NGWF_HAM), intent(inout):: ham

    ! Local Variables
    integer :: ilumo, nlumo
    integer :: num_ngwfs
    integer :: is
    integer :: ierr
    logical :: diagonalise
    type(SPAM3_EMBED), allocatable :: S_invH(:), P(:), PS(:)
    type(SPAM3_EMBED) :: SP, SPS, SPH, SPSPS, SPHPS
    real(kind=DP), allocatable :: eval(:,:)
    type(FUNCTIONS), allocatable, dimension(:,:) :: evec
    real(kind=DP) :: final_max_resid
    real(kind=DP) :: temp_val(2)
    type(NGWF_HAM) :: proj_ham
    integer :: ndopant, zdopant, qdopant
    logical :: loc_cmplx

    ! dense matrix for eigenvectors
    type(DEM), allocatable :: eigs_dens(:)
    ! hamiltonian eigenvalues
    real(kind=DP), allocatable :: eigen_en(:,:)
    ! denskern occupancies
    real(kind=DP), allocatable :: eigen_occ(:,:)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &conduction_ngwf_optimise'

    ! jme: SPINS_DANGER: this subroutine mixes spins
    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine eigenstates_lumo_search (conduction) not ready yet for &
         &more than one k-point.')

    if (pub_on_root) write(stdout,'(a)') 'Beginning LUMO Search'

    ! jmecmplx
    loc_cmplx = rep%overlap%iscmplx

    ! Allocate storage for eigenvectors and eigenvalues
    nlumo = pub_num_eigenvalues
    num_ngwfs = sum(ngwf_basis(:)%num)
    allocate(evec(nlumo,pub_num_spins),stat=ierr)
    call utils_alloc_check('eigenstates_lumo_search','evec',ierr)
    do is=1, pub_num_spins
       do ilumo = 1, nlumo
          call data_functions_alloc(evec(ilumo,is),num_ngwfs,iscmplx=loc_cmplx)
       end do
    end do
    allocate(eval(nlumo,pub_num_spins),stat=ierr)
    call utils_alloc_check('eigenstates_lumo_search','evec',ierr)

    ! Allocate matrix arrays
    allocate(S_invH(pub_num_spins),stat=ierr)
    call utils_alloc_check('eigenstates_lumo_search','S_invH',ierr)
    allocate(P(pub_num_spins),stat=ierr)
    call utils_alloc_check('eigenstates_lumo_search','P',ierr)
    allocate(PS(pub_num_spins),stat=ierr)
    call utils_alloc_check('eigenstates_lumo_search','PS',ierr)

    ! Create matrices
    do is=1,pub_num_spins
       call sparse_embed_create(S_invH(is),rep%inv_overlap,ham%ham(1))
       call sparse_embed_create(P(is),denskern%kern%m(1,PUB_1K))
       call sparse_embed_create(PS(is),denskern%kern%m(1,PUB_1K),rep%overlap)
    end do
    call sparse_embed_create(SP,rep%overlap,denskern%kern%m(1,PUB_1K))
    call sparse_embed_create(SPH,SP,ham%ham(1))
    call sparse_embed_create(SPS,rep%overlap,PS(1))
    call sparse_embed_create(SPSPS,SPS,PS(1))
    call sparse_embed_create(SPHPS,SPH,PS(1))

    ! ndmh: create projected hamiltonian
    call ngwf_ham_create(proj_ham,rep)

    ! Load in Hamiltonian and Overlap from file, if required
    if (pub_read_hamiltonian) then
       call restart_hamiltonian_read(ham%ham)
       call restart_overlap_read(rep%overlap)
    end if

    ! Find atomic number and charge of dopant atoms (defaults to P, 1)
    zdopant = utils_devel_code(15,'LUMO_SEARCH', &
         'ZDOPANT',pub_devel_code)
    qdopant = utils_devel_code(1,'LUMO_SEARCH', &
         'QDOPANT',pub_devel_code)

    ! Adjust number of electrons in n_occ by subtracting number
    ! of dopant atoms
    ndopant = count(mdl%elements(:)%atomic_number==zdopant)
    if (pub_num_spins==1) then
       rep%n_occ(1,PUB_1K) = rep%n_occ(1,PUB_1K) - qdopant * ndopant / 2
    else
       rep%n_occ(1,PUB_1K) = rep%n_occ(1,PUB_1K) - qdopant * ndopant / 2 - modulo(ndopant,2)
       rep%n_occ(2,PUB_1K) = rep%n_occ(2,PUB_1K) - qdopant * ndopant / 2
    end if
    if (pub_on_root) print *,rep%n_occ(:,:)

    ! Construct Inverse overlap matrix with hotelling algorithm
    call sparse_embed_hotelling_init(rep%inv_overlap,rep%overlap)

    if (pub_output_detail>=VERBOSE .and. pub_on_root) then
       write(stdout,'(a)') utils_banner('=', &
            'Calculation of NGWF S^-1 using Hotelling algorithm')
       write(stdout,'(a)')'         Iteration  |  Frobenius norm and Maximum &
            &abs value    (of I-S*S_n^-1)   '
    end if

    ! ndmh: call sparse_mod version of hotelling's algorithm
    call sparse_embed_hotelling_invert(rep%inv_overlap, & !inout
         rep%overlap,show_output=(pub_output_detail>=VERBOSE), &
         max_resid_converged=pub_max_resid_hotelling*0.001_DP, &
         num_iter=pub_maxit_hotelling,final_max_resid=final_max_resid)

    ! ndmh: test if Hotelling failed to invert the matrix
    if (final_max_resid > 0.95_DP) then
          call utils_abort('Error in eigenstates_lumo_search: &
               &Inversion of overlap matrix failed')
    end if

    if (pub_output_detail>=VERBOSE.and.pub_on_root) then
       write(stdout,'(a)') repeat('=',80)
    end if

    ! Call Palser Manolopoulos to generate kernel for sub-gap states
    do is=1,pub_num_spins
       if (pub_maxit_palser_mano < 1) then
          ! Diagonalise the initial Hamiltonian
          call kernel_init_core_ham(denskern%kern%m(is,PUB_1K), &
               ham%ham(is), rep%overlap, rep%n_occ(is,PUB_1K))
       else
          ! Palser-Manolopoulos Canonical Purification
          call palser_mano_kernel_optimise(denskern%kern%m(is,PUB_1K), &
               ham%ham(is), rep%overlap, rep%inv_overlap, &
               rep%n_occ(is,PUB_1K), &
               num_iter=pub_maxit_palser_mano)
       endif
    end do

    temp_val(:) = kernel_bandgap(denskern,ham%ham,rep%overlap,rep%inv_overlap)
    if (pub_on_root) write(stdout,*) 'bandgap: ', temp_val(1:pub_num_spins)

    do is=1,pub_num_spins
       call sparse_embed_trace(temp_val(is),denskern%kern%m(is,PUB_1K),rep%overlap)
       if (pub_on_root) write(stdout,*) 'bandgap: ', temp_val(is)
    end do

    temp_val(:) = 0.0_DP
    do is=1,pub_num_spins
       ! start by computing the largest eigenvalue of the Hamiltonian
       ! this fixes the shift to which states are projected
       call sparse_embed_product(S_invH(is),rep%inv_overlap,ham%ham(is))
       call sparse_embed_extremal_eigenvalue(S_invH(is),rep%overlap,temp_val(is), &
            tol=0.000001_DP)
       temp_val(is) = temp_val(is) + 0.2_DP
    end do
    proj_ham%cond_shift = maxval(temp_val(:))

    if (pub_on_root) write(stdout,'(a,f16.12)') &
        'Shift to apply to valence states: ',proj_ham%cond_shift

    diagonalise = utils_devel_code(.false.,'LUMO_SEARCH', &
         'DIAGONALISE',pub_devel_code)

    ! Show eigenvalues of valence hamiltonian
    if (diagonalise) then
       call eigenstates_create_storage(ngwf_basis, eigen_en, eigen_occ, &
            eigs_dens)
       call eigenstates_calculate(denskern%kern, ham, &
            rep, ngwf_basis, mdl, hfxstate, 'valence', &
            eigs_dens, eigen_en, eigen_occ)
       call eigenstates_destroy_storage(eigen_en, eigen_occ, &
            eigs_dens)
    end if

    ! Loop over spins
    do is=1,pub_num_spins

       ! create a projected Hamiltonian, with all valence states
       ! projected to the shift value.
       call sparse_embed_product(SP,rep%overlap,denskern%kern%m(is,PUB_1K))
       call sparse_embed_product(PS(is),denskern%kern%m(is,PUB_1K),rep%overlap)
       call sparse_embed_product(SPH,SP,ham%ham(is))
       call sparse_embed_product(SPS,SP,rep%overlap)
       call sparse_embed_axpy(SPH,SPS,-proj_ham%cond_shift)
       call sparse_embed_product(SPHPS,SPH,PS(is))
       call sparse_embed_copy(proj_ham%ham(is),ham%ham(is))
       call sparse_embed_axpy(proj_ham%ham(is),SPHPS,-1.0_DP)

    end do

    ! Show eigenvalues obtained from diagonalisation of proj_ham
    if (diagonalise) then
       rep%n_occ(:, PUB_1K) = nlumo
       call eigenstates_create_storage(ngwf_basis, eigen_en, eigen_occ, &
            eigs_dens)
       call eigenstates_calculate(denskern%kern, proj_ham, &
            rep, ngwf_basis, mdl, hfxstate, 'valence', eigs_dens, &
            eigen_en, eigen_occ)
       call eigenstates_destroy_storage(eigen_en, eigen_occ, &
            eigs_dens)
    end if

    ! Loop over spins
    do is=1,pub_num_spins

       ! Create S^{-1} H
       call sparse_embed_product(S_invH(is),rep%inv_overlap,proj_ham%ham(is))

       ! Loop over number of states required
       do ilumo=1,nlumo

          ! Call sparse_extremal_eigenvalue to find kernel for lowest state
          call sparse_embed_extremal_eigenvalue(S_invH(is),rep%overlap, &
               eval(ilumo,is),tol=0.000002_DP, min_val=.true., &
               evec=evec(ilumo,is))

          ! Repeat with lower threshold if failed to converge
          if (utils_isnan(eval(ilumo,is))) then
             if (pub_on_root) write(stdout, '(a,i3,a,i3,a,f16.10)') &
                  'Spin ',is,' Eigenvalue ',ilumo,' search failed. &
                  &Repeating for tol=0.0001'
             call sparse_embed_extremal_eigenvalue(S_invH(is),rep%overlap, &
                  eval(ilumo,is),tol=0.0001_DP, min_val=.true., &
                  evec=evec(ilumo,is))
          end if

          ! Print eigenvalue and save eigenvector
          ! create kernel as the outer product of the eigenvector
          if (loc_cmplx) then
             call sparse_embed_outer_product(P(is), &
                  evec(ilumo,is)%z,evec(ilumo,is)%z)  !jmecmplx
          else
             call sparse_embed_outer_product(P(is), &
                  evec(ilumo,is)%d,evec(ilumo,is)%d)
          end if
          if (pub_popn_mulliken_partial) then
             call population_analysis_mulliken(rep%overlap, P(is:is), &
                  mdl, ngwf_basis, ilumo)
          end if

          ! Subtract outer product of eigenvector from kernel
          ! project out this value from the hamiltonian
          temp_val(is) = -(proj_ham%cond_shift-eval(ilumo,is))
          call sparse_embed_product(PS(is),P(is),rep%overlap)
          call sparse_embed_axpy(S_invH(is),PS(is),-1.0_DP*temp_val(is))

          if (pub_on_root) write(stdout, '(a,i3,a,i3,a,1f16.10)') &
               'Spin ',is,' Eigenvalue ',ilumo,':',eval(ilumo,is)

       ! End loop over number of states
       end do

    end do

    if (pub_on_root) write(stdout,'(a)') 'Completed LUMO Search'

    ! ndmh: create projected hamiltonian
    call ngwf_ham_destroy(proj_ham)

    ! Destroy matrices
    call sparse_embed_destroy(SPHPS)
    call sparse_embed_destroy(SPSPS)
    call sparse_embed_destroy(SPS)
    call sparse_embed_destroy(SPH)
    call sparse_embed_destroy(SP)
    do is=pub_num_spins,1,-1
       call sparse_embed_destroy(PS(is))
       call sparse_embed_destroy(P(is))
       call sparse_embed_destroy(S_invH(is))
    end do

    ! Deallocate matrix arrays
    deallocate(PS,stat=ierr)
    call utils_dealloc_check('eigenstates_lumo_search','PS',ierr)
    deallocate(P,stat=ierr)
    call utils_dealloc_check('eigenstates_lumo_search','P',ierr)
    deallocate(S_invH,stat=ierr)
    call utils_dealloc_check('eigenstates_lumo_search','S_invH',ierr)

    ! Deallocate storage for eigenvectors and eigenvalues
    deallocate(eval,stat=ierr)
    call utils_dealloc_check('eigenstates_lumo_search','evec',ierr)
    do is = pub_num_spins, 1, -1
       do ilumo = nlumo, 1, -1
          call data_functions_dealloc(evec(ilumo,is))
       end do
    end do
    deallocate(evec, stat=ierr)
    call utils_dealloc_check('eigenstates_lumo_search','evec',ierr)

  end subroutine eigenstates_lumo_search


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_print_ens_occs(eigen_en, eigen_occ, s_occ, neig, &
       eigen_cond_occ, cond_n_occ)

    !=================================================================!
    ! This subroutine prints the eigenvalues of the Hamiltonian and   !
    ! the density kernel and prints maximum and minimum values as     !
    ! well as a range of them around the band gap.                    !
    !-----------------------------------------------------------------!
    ! NOTE:                                                           !
    !   The occupancies printed by this subroutine are the projection !
    !   of the orbitals |phi_c> M^c_n onto the density operator       !
    !   given by |phi_a>K^ab<phi_b|.                                  !
    !-----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 21/05/2006 by modifying     !
    ! internal_find_eigen_info which was written by Chris-Kriton      !
    ! Skylaris and Peter Haynes in 2004-2005.                         !
    ! Modified format strings, and enabled printing of wider range    !
    ! of unoccupied states if required, by Nicholas Hine on 19/10/09. !
    ! Moved to eigenstates_mod by Nicholas Hine in November 2011.     !
    !=================================================================!

    use comms, only: pub_on_root, comms_barrier
    use constants, only: DP, stdout
    use rundat, only: pub_num_eigenvalues, pub_print_qc, pub_num_spins, &
                      pub_edft
    use utils, only: utils_qc_print, utils_banner
    use ensemble_dft, only: edft_zerok_fermi

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: eigen_en(:,:)  ! orbital energies from Hamiltonian
    real(kind=DP), intent(in) :: eigen_occ(:,:) ! orbital occupancies from density kernel
    real(kind=DP), intent(in) :: s_occ(:)       ! number of occupied orbitals of each spin
    integer, intent(in)        :: neig           ! number of NGWFs
    ! optional argument for orbital occupancies from conduction density kernel
    real(kind=DP), optional, intent(in) :: eigen_cond_occ(:,:)
    integer, optional, intent(in)        :: cond_n_occ(:)   ! number of occupied cond orbitals of each spin

    ! Local variables
    integer :: is         ! Spin counter
    integer :: ieig       ! Loop counter
    integer :: num_print  ! Number of eigenvalues to actually print
    logical :: cond_occs
    real(kind=dp), dimension(pub_num_spins) :: zerok_fermi
    integer :: iorb, lumo_orb, homo_orb
    character(len=57) :: gap_string

    cond_occs = present(eigen_cond_occ).and.present(cond_n_occ)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &eigenstates_print_ens_occs'

    ! Diagonalisation is not parallelised, so only work on root
    if (pub_on_root) then

       ! kkbd: Work out 0K Fermi energy
       if (pub_edft) then
          call edft_zerok_fermi(zerok_fermi, s_occ(:), eigen_en)
       end if

       ! pdh: loop over spins
       do is=1,pub_num_spins

          ! kkbd: Find 0K HOMO and LUMO in EDFT runs
          if (pub_edft) then
             homo_orb = 0
             lumo_orb = 0
             do iorb=1,neig
                if (eigen_en(iorb,is) > zerok_fermi(is)) then
                   ! kkbd: previous orbital was HOMO of this spin channel
                   !       and this one is LUMO
                   lumo_orb = iorb
                   homo_orb = iorb - 1
                   exit
                end if
             end do
          else
             lumo_orb = nint(s_occ(is)) + 1
             homo_orb = nint(s_occ(is))
          end if

          ! Print out the required results
          if (pub_num_spins == 1) then
             write(stdout,'(/a/)') utils_banner('=', &
                  'Orbital energy and occupancy information')
          else
             write(stdout,'(/a,i1,a/)') '============== &
                  &Orbital energy and occupancy information for spin ',is, &
                  ' ============='
          end if
          write(stdout,'(a,i8)')   '      Total number of orbitals:', neig
          !write(stdout,'(a,i8)') '   Number of occupied orbitals:', n_occ(is)
          if (.not.pub_edft) then
             write(stdout,'(a,f18.9)') '   Number of occupied orbitals:', s_occ(is)
          else
             write(stdout,'(a,f18.9)') '        Occupancy sum (kernel):', s_occ(is)
          end if
          write(stdout,'(a,f18.9)') '   Occupancy sum (eigenvalues):', &
               sum(eigen_occ(:,is))
          if (cond_occs) then
             write(stdout,'(a,i8)')    ' Num of occupied cond orbitals:', cond_n_occ(is)
             write(stdout,'(a,f18.9)') '            Cond Occupancy sum:', &
                  sum(eigen_cond_occ(:,is))
          end if

          ! pdh: check for no states occupied as well as all
          if (homo_orb<neig .and. homo_orb > 0) then
             if( eigen_occ(homo_orb,is)-eigen_occ(lumo_orb,is) < 0.0_DP ) then
                write(stdout,'(a,f18.9,a)') '                 HOMO-LUMO gap:', &
                     eigen_en(homo_orb,is) - eigen_en(lumo_orb,is),' Eh'
             else
                write(stdout,'(a,f18.9,a)') '                 HOMO-LUMO gap:', &
                     eigen_en(lumo_orb,is) - eigen_en(homo_orb,is),' Eh'
             end if
             write(stdout,'(a,f18.9,a)') '                 Mid-gap level:', &
                  0.5_DP*(eigen_en(homo_orb,is) + eigen_en(lumo_orb,is)),' Eh'
          else if (lumo_orb==0) then
             write(stdout,'(a,a)') '                 HOMO-LUMO gap:', &
                  '  Unavailable - no orbitals occupied'
             write(stdout,'(a,a)') '                 Mid-gap level:', &
                  '  Unavailable - no orbitals occupied'
          else
             write(stdout,'(a,a)') '                 HOMO-LUMO gap:', &
                  '  Unavailable - all orbitals occupied'
             write(stdout,'(a,a)') '                 Mid-gap level:', &
                  '  Unavailable - all orbitals occupied'
          end if
          write(stdout,'(/a)',advance='no') '                        &
               &Orbital | Energy (Eh) | Occupancy'
          if (cond_occs) then
             write(stdout,'(a)') ' | Cond Occupancy'
          else
             write(stdout,*)
          end if

          ! pdh: don't print lowest state if no states are occupied
          if (lumo_orb > 0) then
             write(stdout,'(a,i6,f16.9,f12.7)',advance='no') &
                  '                       ', &
                  1,eigen_en(1,is), eigen_occ(1,is)
             if (cond_occs) then
                write(stdout,'(f12.7)') eigen_cond_occ(1,is)
             else
                write(stdout,*)
             end if
          end if
          ! pdh: avoid repeating lowest state
          num_print = max(homo_orb-pub_num_eigenvalues,1)+1

          if (num_print > 2) write(stdout,'(a)') &
               '                        .......   ...........   .........'

          do ieig=num_print,homo_orb
             write(stdout,'(a,i6,f16.9,f12.7)',advance='no') &
                  '                       ', &
                  ieig, eigen_en(ieig,is), eigen_occ(ieig,is)
             if (cond_occs) then
                write(stdout,'(f12.7)') eigen_cond_occ(ieig,is)
             else
                write(stdout,*)
             end if
          end do

          ! pdh: don't print gap if no states are occupied
          if (pub_edft) then
             gap_string = '                        &
                          &....... -- gap at 0K -- .........'
          else
             gap_string = '                        &
                          &.......   --- gap ---   .........'
          end if
          if (lumo_orb > 0) write(stdout,'(a)') gap_string

          ! Don't print more virtual orbitals than are available
          num_print = min(homo_orb+pub_num_eigenvalues,neig-1)

          do ieig=lumo_orb,num_print
             write(stdout,'(a,i6,f16.9,f12.7)',advance='no') &
                  '                       ', &
                  ieig, eigen_en(ieig,is), eigen_occ(ieig,is)
             if (cond_occs) then
                write(stdout,'(f12.7)') eigen_cond_occ(ieig,is)
             else
                write(stdout,*)
             end if
          end do
          if (num_print < neig-1) write(stdout,'(a)') '                        &
               &.......   ...........   .........'

          ! pdh: don't print highest state if all states are occupied
          if (lumo_orb < neig) then
             write(stdout,'(a,i6,f16.9,f12.7)',advance='no') &
                  '                       ', &
                  neig, eigen_en(neig,is),eigen_occ(neig,is)
             if (cond_occs) then
                write(stdout,'(f12.7)') eigen_cond_occ(neig,is)
             else
                write(stdout,*)
             end if
          end if

          if (pub_print_qc) then
             write(stdout,*)
             call utils_qc_print('energy_eigenvalue_1',eigen_en(1,is))
             call utils_qc_print('energy_eigenvalue_n',eigen_en(neig,is))
             call utils_qc_print('occupancy_1',eigen_occ(1,is))
             call utils_qc_print('occupancy_n',eigen_occ(neig,is))
             if (cond_occs) then
                call utils_qc_print('cond_occupancy_1',eigen_cond_occ(1,is))
                call utils_qc_print('cond_occupancy_n',eigen_cond_occ(neig,is))
             end if
          end if

          write(stdout,'(/a)') repeat('=',80)

       end do

    end if

    call comms_barrier

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &eigenstates_print_ens_occs'

  end subroutine eigenstates_print_ens_occs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! jd (retreat2014): Removed s_eigs_dens -- was unreferenced
  subroutine eigenstates_transition_denpot(eigs_dens, &
       ngwf_basis,rep,mdl)

    !==================================================================!
    ! This subroutine calculates the interactions of a pair of         !
    ! transition densities.                                            !
    !------------------------------------------------------------------!
    ! Routine created by Nicholas Hine in July 2013                    !
    ! Modified for embedding by Joseph Prentice, September 2018        !
    !==================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, HARTREE_IN_EVS, ANGSTROM
    use cutoff_coulomb, only: cutoff_coulomb_hartree
    use dense, only: DEM, dense_convert, dense_rank1_update, &
         dense_create, dense_destroy, dense_scale
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use hartree, only: hartree_on_grid, hartree_via_multigrid
    use model_type, only: MODEL
    use is_solvation, only: have_rho_ion
    use is_smeared_ions, only: rho_ion
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_coulomb_cutoff, pub_num_spins, &
         pub_multigrid_hartree, pub_devel_code, pub_is_bulk_permittivity,&
         pub_turn_off_hartree
    use sparse, only: SPAM3
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_devel_code
    use visual, only: visual_scalarfield

    implicit none

    ! Arguments
    type(DEM), intent(in) :: eigs_dens
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(NGWF_REP), intent(in) :: rep
    type(MODEL), intent(in) :: mdl

    ! Local variables
    integer :: is, ieig, jeig
    integer :: ierr, num
    real(kind=DP) :: eps_trans, eps_normal
    type(SPAM3_EMBED) :: trans_kern(1:pub_num_spins)
    type(SPAM3) :: kern_array(1:pub_num_spins)
    type(DEM) :: trans_kern_dens
    real(kind=DP), allocatable :: trans_dens_fine(:,:,:,:)
    real(kind=DP), allocatable :: buffer_fine(:,:,:,:)
    real(kind=DP), allocatable :: trans_hart_fine(:,:,:,:)
    character(len=256) :: output_file  ! file names
    character(len=256) :: title_line   ! title line in output file
    character(len=16) :: tmp1,tmp2
    ! agrecocmplx
    logical :: loc_cmplx

    integer :: isub, jsub

    ! agrecocmplx: for the case of complex NGWFs,
    ! also eigs_dens is complex
    loc_cmplx = rep%ngwfs_on_grid(1)%iscmplx

    ! jcap: find total number of NGWFs
    num=sum(ngwf_basis(:)%num)

    do is=1,pub_num_spins
       trans_kern(is)%structure = 'K'//rep%postfix
       ! agrecocmplx
       call sparse_embed_create(trans_kern(is),iscmplx=loc_cmplx)
    end do

    ! agrecocmplx
    call dense_create(trans_kern_dens,num,num,iscmplx=loc_cmplx)
    call dense_scale(trans_kern_dens,0.0_DP)

    ieig = utils_devel_code(-1,'TRANSITION_DENPOT','IEIG',pub_devel_code)
    jeig = utils_devel_code(-1,'TRANSITION_DENPOT','JEIG',pub_devel_code)
    jeig = utils_devel_code(-1,'TRANSITION_DENPOT','JEIG',pub_devel_code)

    call dense_rank1_update(trans_kern_dens,eigs_dens,eigs_dens, &
         ieig,jeig,1.00_DP)

    call dense_convert(trans_kern(1),trans_kern_dens)

    ! ndmh: Allocate charge density slabs
    allocate(trans_dens_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('eigenstates_transition_denpot &
         &(properties_calculate)','trans_dens_fine',ierr)
    allocate(buffer_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('eigenstates_transition_denpot &
         &(properties_calculate)','buffer_fine',ierr)
    allocate(trans_hart_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('eigenstates_transition_denpot &
         &(properties_calculate)','trans_hart_fine',ierr)

    if (pub_on_root) write(stdout,*) 'Calculating transition density for i=', &
         ieig,' j=',jeig

    ! ndmh: Calculate charge associated with transition density
    ! jcap: loop over regions
    trans_dens_fine=0.0_DP
    do isub=1,mdl%nsub
       do jsub=1,mdl%nsub
          ! jcap: create temporary kernel
          call sparse_embed_extract_from_array(kern_array,trans_kern,isub,jsub)
          buffer_fine=0.0_DP
          call density_on_grid(buffer_fine,mdl%fine_grid,mdl%dbl_grid, &
               mdl%cell,mdl%fftbox,kern_array,rep%ngwf_overlap%m(isub,jsub),&
               rep%ngwfs_on_grid(isub),ngwf_basis(isub),&
               rep%ngwfs_on_grid(jsub),ngwf_basis(jsub))
          trans_dens_fine=trans_dens_fine+buffer_fine
          call sparse_embed_destroy_extracted_array(kern_array)
       end do
    end do

    ! ndmh: write i,j to temp strings
    write(tmp1,*) ieig
    write(tmp2,*) jeig

    ! ndmh: output transition density to plotfile formats
    ! ndmh: output orbital density in Angstrom^-3 rather than in Bohr^-3
    write(output_file,'(4a)') '_trans_dens_i',trim(adjustl(tmp1)), &
         'j',trim(adjustl(tmp2))
    write(title_line,'(4a)') 'Transition density for i=',trim(adjustl(tmp1)), &
         'and j=',trim(adjustl(tmp2))
    call visual_scalarfield( &
         trans_dens_fine(:,:,:,1), mdl%fine_grid, mdl%cell, title_line, &
         output_file, mdl%elements, ANGSTROM**3)  ! input

    if (pub_on_root) write(stdout,*) 'Calculating Hartree potential of &
         &transition density for i=',ieig,' j=',jeig

    if (.not. pub_turn_off_hartree) then
       if (pub_coulomb_cutoff) then
          call cutoff_coulomb_hartree(trans_hart_fine, trans_dens_fine, &
               mdl)
       else if (pub_multigrid_hartree) then
          eps_normal = pub_is_bulk_permittivity
          eps_trans = utils_devel_code(eps_normal,'TRANSITION_DENPOT','EPS', &
               pub_devel_code)
          pub_is_bulk_permittivity = eps_trans
          have_rho_ion = .false.
          rho_ion = 0.0_DP
          call hartree_via_multigrid(trans_hart_fine, trans_dens_fine, &
               mdl%fine_grid,mdl%cell,elements=mdl%elements)
          pub_is_bulk_permittivity = eps_normal
       else
          call hartree_on_grid(trans_hart_fine, trans_dens_fine, &
               mdl%fine_grid, mdl%cell)
       end if
    else
       trans_hart_fine = 0.0_DP
    endif

    ! ndmh: output transition density to plotfile formats
    ! ndmh: output orbital density in Angstrom^-3 rather than in Bohr^-3
    write(output_file,'(4a)') '_trans_hart_i',trim(adjustl(tmp1)), &
         'j',trim(adjustl(tmp2))
    write(title_line,'(4a)') 'Hartree Pot of transition density for i=', &
         trim(adjustl(tmp1)),'and j=',trim(adjustl(tmp2))
    call visual_scalarfield( &
         trans_hart_fine(:,:,:,1), mdl%fine_grid, mdl%cell, title_line, &
         output_file,mdl%elements, HARTREE_IN_EVS)  ! input

    ! ndmh: Dellocate charge density slabs
    deallocate(trans_hart_fine,stat=ierr)
    call utils_dealloc_check('eigenstates_transition_denpot &
         &(properties_calculate)','trans_hart_fine',ierr)
    deallocate(trans_dens_fine,stat=ierr)
    call utils_dealloc_check('eigenstates_transition_denpot &
         &(properties_calculate)','trans_dens_fine',ierr)
    deallocate(buffer_fine,stat=ierr)
    call utils_dealloc_check('eigenstates_transition_denpot &
         &(properties_calculate)','buffer_fine',ierr)


    call dense_destroy(trans_kern_dens)

    do is=pub_num_spins,1,-1
       call sparse_embed_destroy(trans_kern(is))
    end do

  end subroutine eigenstates_transition_denpot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eigenstates_fermi_level(eigen_en, s_occ, num, fermi_e)

    !==================================================================!
    ! This subroutine calculates the Fermi level for shifting DOS.     !
    !==================================================================!
    use constants, only: DP
    use rundat, only: pub_num_spins, pub_edft
    use ensemble_dft, only: edft_zerok_fermi

    implicit none

    real(kind=DP), intent(in) :: eigen_en(:,:) ! molecular orbital energies
    real(kind=dp), intent(in) :: s_occ(:)      ! total occupancy of spin channel
    integer, intent(in)       :: num           ! number of ngwfs
    real(kind=DP), intent(out) :: fermi_e      ! Fermi level for shifting DOS

    integer, dimension(pub_num_spins) :: n_occ
    real(kind=dp), dimension(pub_num_spins) :: zerok_fermi
    integer :: is

    fermi_e = 0.0_DP

    if (pub_edft) then
       ! kkbd: use the 0K Fermi level based on calculation context
       call edft_zerok_fermi(zerok_fermi,s_occ,eigen_en)
       ! kkbd: Convention was to average spin channel fermi levels...
       fermi_e = sum(zerok_fermi) / real(pub_num_spins,kind=dp)
    else
       n_occ = nint(s_occ)
       do is=1,pub_num_spins
          ! cks: assume Fermi energy as average of HOMO and LUMO and spins
          if (n_occ(is)<num.and.n_occ(is)>0) then
             fermi_e = fermi_e + 0.5_DP*(eigen_en(n_occ(is)+1,is)+ &
                  eigen_en(n_occ(is),is))
             ! ndmh: unless we have all orbitals occupied
          else if (n_occ(is) > 0) then
             fermi_e = fermi_e + eigen_en(n_occ(is),is) + tiny(1.0_DP)
          else
             fermi_e = fermi_e + eigen_en(1,is) + tiny(1.0_DP)
          end if
       end do
       fermi_e = fermi_e / pub_num_spins
    end if

  end subroutine eigenstates_fermi_level

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  subroutine eigenstates_print_eigenvecs_comp(eigs_dens, elements, ngwf_basis, n_occ)
!
!    !=================================================================!
!    ! This subroutine prints the relative contributions of NGWFs to   !
!    ! the Kohn-Sham orbitals from regions specified by the user.      !
!    !                                                                 !
!    !-----------------------------------------------------------------!
!    ! Arguments:                                                      !
!    ! eigs_dens  :  Matrix holding the eigenvectors c (Hc=ESc).       !
!    ! elements   :  Elements list containing atom region information. !
!    ! ngwf_basis :  Needed to identify the relevant atom for each NGWF!
!    !                                                                 !
!    !-----------------------------------------------------------------!
!    ! Written by Robert Charlton, 30/11/2016.                         !
!    !=================================================================!
!
!    use comms, only: pub_on_root, comms_barrier
!    use constants, only: DP, stdout
!    use dense, only: DEM
!    use ion, only: ELEMENT
!    use rundat, only: pub_num_spins
!
!    implicit none
!
!    ! Arguments
!    type(DEM), intent(inout) :: eigs_dens(pub_num_spins) ! orbital energies from Hamiltonian
!    type(ELEMENT), intent(in) :: elements(:)             ! elements list
!    type(FUNC_BASIS), intent(in) :: ngwf_basis           ! Function basis
!    ! rc2013: n_occ not used at the moment
!    integer, intent(in)        :: n_occ(:)       ! number of occupied orbitals of each spin
!
!    ! Local variables
!    integer :: is         ! Spin counter
!    integer :: num_print  ! Number of eigenvalues to actually print
!    integer :: orb_i      ! Orbital counter
!    integer :: orb_j      ! NGWF counter
!    integer :: num        ! No. of NGWFs
!    real    :: ngwf_sum   ! Normalisation constant
!    real    :: ngwf_tot1  ! Total NGWFs in MOL1 (replace with array at some point)
!    real    :: ngwf_tot2
!    real(kind=DP), dimension(ngwf_basis%num) :: ngwf_cont ! NGWF contributions to each orbital
!
!    !cond_occs = present(eigen_cond_occ).and.present(cond_n_occ)
!    num = ngwf_basis%num
!
!    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
!         &eigenstates_print_eigenvecs_comp'
!
!    ! Diagonalisation is not parallelised, so only work on root
!    if (pub_on_root) then
!       ! pdh: loop over spins
!       do is=1,pub_num_spins
!          ! Print out the required results
!          if (pub_num_spins == 1) then
!             write(stdout,'(/a/)') '=================== &
!                  &Orbital composition information ==================='
!          else
!             write(stdout,'(/a,i1,a/)') '============== &
!                  &Orbital composition information for spin ',is, &
!                  ' ============='
!          end if
!          write(stdout,'(/a)') '                        &
!               &Orbital | MOL1 (%) | MOL2 (%)'
!             ! rc2013: should adapt to look more like regular eigenvalue printing
!             ! at some point...
!             do orb_i=1,num
!             !do orb_i=num_print,n_occ(is)
!                ngwf_tot1 = 0
!                ngwf_tot2 = 0
!                ! Loop over NGWFs
!                do orb_j=1,num
!                   ! rc2013: not entirely sure why but we should take the trace of the
!                   ! eigenvectors matrix. This gives "reasonable" results.
!                   ngwf_cont(orb_j) = eigs_dens(is)%dmtx(orb_j,orb_i)*&
!                        eigs_dens(is)%dmtx(orb_j,orb_i)
!                   ! MPI warning! This may need to be adapted for calcs on multiple cpus.
!                   select case (elements(ngwf_basis%atom_of_func(orb_j))%ngwf_constraint)
!                   case('MOL1')
!                      ngwf_tot1 = ngwf_tot1 + ngwf_cont(orb_j)
!                   case('MOL2')
!                     ngwf_tot2 = ngwf_tot2 + ngwf_cont(orb_j)
!                   end select
!                end do
!                ! Total NGWF contributions to this KS orbital
!                ngwf_sum  = ngwf_tot1+ngwf_tot2
!                ! Normalise the array of NGWF contributions
!                ngwf_cont = ngwf_cont/ngwf_sum
!                ngwf_tot1 = 100.0_DP*ngwf_tot1/ngwf_sum
!                ngwf_tot2 = 100.0_DP*ngwf_tot2/ngwf_sum
!                write(stdout,'(a23,i6,f13.3,f11.3)') '                        ',&
!                   orb_i, ngwf_tot1, ngwf_tot2
!             end do
!       end do
!       write(stdout,'(/a)') '=========================================&
!          &======================================='
!    end if
!
!    call comms_barrier
!
!    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
!         &eigenstates_print_eigenvecs_comp'
!
!  end subroutine eigenstates_print_eigenvecs_comp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module eigenstates

