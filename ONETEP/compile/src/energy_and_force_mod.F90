! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                      E N E R G Y    A N D    F O R C E                      !
!=============================================================================!
!                                                                             !
! This module calculates the total energy and forces, as well as converged    !
! ngwfs, density kernel and related quantities, for a given ionic             !
! configuration.                                                              !
!                                                                             !
!-----------------------------------------------------------------------------!
!                                                                             !
!                 The ONETEP code is written and maintained by                !
! Chris-Kriton Skylaris, Arash A. Mostofi, Nicholas Hine and Peter D. Haynes  !
!                                                                             !
!-----------------------------------------------------------------------------!
! Originally written by Chris-Kriton Skylaris in 2000.                        !
! Improved and parallelised by Chris-Kriton Skylaris in November 2003.        !
! Improvements by Peter D. Haynes in 2004.                                    !
! Turned into a module by Arash A. Mostofi, Version 0.01, 01/11/2004.         !
! Modified by Chris-Kriton Skylaris on 08/11/2004 to work with ppd_strategy.  !
! Modified by Nicholas Hine on 11/08/2008 to remove legacy SPAM denskern      !
! Modified by Nicholas Hine in July 2009 for SPAM3 and function_basis.        !
! Modified by Nicholas Hine in November 2009 to remove use of workspace_mod.  !
! PAW functionality by Nicholas Hine, May-July 2010.                          !
! Conduction NGWF optimisation added by Laura Ratcliff in October 2010.       !
! Reorganisation to split off energy_and_force_init/exit_cell as separate     !
! routines, by Nicholas Hine in May 2012.                                     !
! Restructured for embedding by Robert Charlton, August 2018.                 !
!-----------------------------------------------------------------------------!

module energy_and_force

  implicit none

  private

  public :: energy_and_force_calculate
  public :: energy_and_force_init_cell
  public :: energy_and_force_exit_cell

contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine energy_and_force_init_cell(mdl)

    use cell_grid, only: cell_grid_distribute
    use comms, only: pub_on_root
    use constants, only: CRLF, DP, EDA_INIT, stdout
    use cutoff_coulomb, only: cutoff_coulomb_init
    use dense, only: dense_init
    use density_init, only: density_init_radial_init
    use fft_box, only: fftbox_init, fftbox_find_size
    use fourier, only: fourier_init_cell, fourier_init_threads
    use k_points, only: kpoints_generate_unique_list
    use model_type, only: MODEL
    use multigrid_methods, only: multigrid_initialise
    use ngwfs_radial, only: ngwfs_radial_set_nfunctions, ngwfs_radial_init
    use openbc_locps, only: openbc_locps_init
    use paw, only: paw_read_species, paw_exc_core_init
    use ppd_strategy, only: ppd_strategy_check_and_print
    use pseudopotentials, only: pseudopotentials_read_species
    use rundat, only: pub_coulomb_cutoff, pub_initial_dens_realspace, &
         pub_devel_code, pub_fine_grid_scale, pub_dbl_grid_scale, &
         pub_paw, pub_fine_is_dbl, pub_dbl_is_std, pub_aug, &
         pub_aug_den_dim, pub_is_implicit_solvent, pub_ngwf_halo, &
         pub_is_smeared_ion_rep, pub_use_remote_ngwfs, pub_use_swx, &
         pub_pspot_bc_is_periodic, pub_fftbox_pref, pub_extend_ngwf, &
         pub_use_hfx, pub_dma_calculate, rundat_set_pub_use_swx, &
         pub_cond_calculate_any_task, pub_use_aux_ngwfs, &
         pub_hubbard_atomsolve, pub_eda, pub_eda_mode, pub_any_nl_proj, &
         pub_kpoints_specified_list, pub_kpoint_list, &
         pub_multigrid_in_use, pub_multigrid_bc_is_periodic, &
         pub_kp_grid_size, pub_kp_grid_shift, &
         pub_use_activehfx, rundat_set_pub_use_activeswx, &
         pub_active_region, pub_emft, pub_emft_follow, pub_do_fandt, &
         pub_cdft, pub_dmft_points, pub_foe, pub_hubbard, &
         pub_pol_emb_qmstar, pub_write_nbo, pub_use_activeswx, &
         pub_do_tddft, pub_lr_phonons_calculate, pub_xc_ke_density_required, &
         pub_block_orthogonalise, pub_build_bo, pub_usp, pub_nlcc
    use services, only: services_flush
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_devel_code, utils_assert, &
         utils_banner
    use xc, only: xc_hfxinit, xc_init, pub_xc_gradient_corrected, &
         xc_embed_swap_functional

    implicit none

    ! Arguments
    type(MODEL), target, intent(inout)   :: mdl

    ! Local Variables
    integer :: n1, n2, n3
    integer :: ierr
    integer :: iat, glob_iat, isub
    logical :: constrain_for_mg
    character(2) :: isub_str
    real(kind=DP) :: max_rad
    ! jcap: temporary global variables
    logical :: orig_pub_aug, orig_pub_nlcc, orig_pub_usp
    character(len=*), parameter :: myself = 'energy_and_force_init_cell'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    ! ndmh: initialise FFT threading
    call fourier_init_threads

    ! cks: initialise the PSEUDO_SPECIES type array and set up
    ! cks: projector numbers, and parts of elements array
    nullify(mdl%pseudo_sp)
    nullify(mdl%paw_sp)
    if (.not.pub_paw) then
       ! rc2013: PAR_WARNING!
       if(pub_on_root) write(stdout,'(a)') utils_banner('<', 'Full system')
       call pseudopotentials_read_species(mdl%pseudo_sp,mdl%elements,mdl%cell, &
            mdl%species, mdl%par, mdl%any_nl_proj) !input/output
       allocate(mdl%paw_sp(size(mdl%pseudo_sp)),stat=ierr)
       call utils_alloc_check(myself,'mdl%paw_sp',ierr)
    end if
    if (pub_paw) then
       ! rc2013: PAR_WARNING!
       if(pub_on_root) write(stdout,'(a)') utils_banner('<', 'Full system')
       call paw_read_species(mdl%paw_sp,mdl%elements,mdl%species,mdl%par)
       allocate(mdl%pseudo_sp(size(mdl%paw_sp)),stat=ierr)
       call utils_alloc_check(myself,'mdl%pseudo_sp',ierr)
    end if

    ! rc2013: cycle through embedding regions
    ! There must be a better way of avoiding doing the regions and mdl
    ! separately?
    ! jcap: if we only have one region, just point to full system arrays
    if (mdl%nsub.eq.1) then
       nullify(mdl%regions(1)%pseudo_sp)
       nullify(mdl%regions(1)%paw_sp)
       mdl%regions(1)%pseudo_sp => mdl%pseudo_sp
       mdl%regions(1)%paw_sp => mdl%paw_sp
       mdl%regions(1)%any_nl_proj = mdl%any_nl_proj
       ! jcap: elements information needs to be set up individually to
       ! avoid copying over previously set up data that is different
       ! for the region and overall elements arrays
       mdl%regions(1)%elements(:)%pspecies_number = &
            mdl%elements(:)%pspecies_number
       mdl%regions(1)%elements(:)%atomic_number = &
            mdl%elements(:)%atomic_number
       mdl%regions(1)%elements(:)%nprojectors = &
            mdl%elements(:)%nprojectors
       mdl%regions(1)%elements(:)%npawpws = mdl%elements(:)%npawpws
       mdl%regions(1)%elements(:)%ncorewfs = mdl%elements(:)%ncorewfs
       mdl%regions(1)%elements(:)%max_core_radius = &
            mdl%elements(:)%max_core_radius
       mdl%regions(1)%elements(:)%max_core_wf_radius = &
            mdl%elements(:)%max_core_wf_radius
       mdl%regions(1)%elements(:)%ion_charge = mdl%elements(:)%ion_charge
       if (mdl%any_nl_proj) pub_any_nl_proj = .true.
       ! ndmh: determine number of NGWFs on each atom if this is set to automatic
       if (.not.(pub_eda) .or. (pub_eda .and. pub_eda_mode == EDA_INIT)) &
            call ngwfs_radial_set_nfunctions(mdl%regions(1))
    else
       ! jcap: need to save the full system values of these global
       ! variables to prevented them being overwritten
       orig_pub_usp = pub_usp
       orig_pub_aug = pub_aug
       orig_pub_nlcc = pub_nlcc
       do isub=1,mdl%nsub
          write(isub_str,'(i0)')isub
          nullify(mdl%regions(isub)%pseudo_sp)
          nullify(mdl%regions(isub)%paw_sp)
          if (.not.pub_paw) then
             if(pub_on_root) write(stdout,'(a)') &
                  utils_banner('<', 'Region '//isub_str)
             call pseudopotentials_read_species(mdl%regions(isub)%pseudo_sp, &
                  mdl%regions(isub)%elements, mdl%cell, mdl%regions(isub)%species,&
                  mdl%regions(isub)%par, mdl%regions(isub)%any_nl_proj) !in/out
             allocate(mdl%regions(isub)%paw_sp( &
                  mdl%regions(isub)%par%num_pspecies),stat=ierr)
             call utils_alloc_check(myself,'mdl%regions(isub)%paw_sp',ierr)
             ! rc2013: only allocate pub_any_nl_proj when we've checked all regions
             ! rc2013: CHECK!
             if (mdl%regions(isub)%any_nl_proj) pub_any_nl_proj = .true.
          end if
          if (pub_paw) then
             if(pub_on_root) write(stdout,'(a)') &
                  utils_banner('<', 'Region '//isub_str)
             call paw_read_species(mdl%regions(isub)%paw_sp, &
                  mdl%regions(isub)%elements, mdl%regions(isub)%species, &
                  mdl%regions(isub)%par) !input/output
             allocate(mdl%regions(isub)%pseudo_sp( &
                  mdl%regions(isub)%par%num_pspecies),stat=ierr)
             call utils_alloc_check(myself,'mdl%regions(isub)%pseudo_sp',ierr)
          end if

          ! ndmh: determine number of NGWFs on each atom if this is set to automatic
          if (.not.(pub_eda) .or. (pub_eda .and. pub_eda_mode == EDA_INIT)) &
               call ngwfs_radial_set_nfunctions(mdl%regions(isub))

       end do
       ! jcap: revert to full system values for these global variables
       pub_usp = orig_pub_usp
       pub_aug = orig_pub_aug
       pub_nlcc = orig_pub_nlcc
    end if

    ! Print atom counting information
    if (pub_on_root) call internal_print_num_species(mdl)
    call services_flush

    ! jcap: make checks that we only have one subregion - several
    ! functionalities only work for one subregion at the moment
    if (mdl%nsub.gt.1) then
       call utils_assert(.not.pub_dma_calculate, &
         'DMA is not yet compatible with embedding calculations!'//CRLF//&
         'Please only use one subregion when doing a DMA calculation.')
       call utils_assert(.not.pub_hubbard, &
         'Hubbard is not yet compatible with embedding calculations!'//CRLF//&
         'Please only use one subregion when doing a Hubbard calculation.')
       call utils_assert(.not.pub_paw, &
         'PAW is not yet compatible with embedding calculations!'//CRLF//&
         'Please only use one subregion when doing a PAW calculation.')
       call utils_assert(.not.pub_eda, &
         'EDA is not yet compatible with embedding calculations!'//CRLF//&
         'Please only use one subregion when doing a EDA calculation.')
       call utils_assert(.not.pub_write_nbo, &
         'Writing NBO is not yet compatible with embedding calculations!'//CRLF//&
         'Please only use one subregion when writing out a NBO analysis.')
       call utils_assert(.not.pub_cdft, &
         'CDFT is not yet compatible with embedding calculations!'//CRLF//&
         'Please only use one subregion when doing a CDFT calculation.')
       call utils_assert(pub_dmft_points==0, &
         'DMFT is not yet compatible with embedding calculations!'//CRLF//&
         'Please only use one subregion when doing a DMFT calculation.')
       call utils_assert(.not.pub_pol_emb_qmstar, &
         'QM/MM is not yet compatible with quantum embedding calculations!'//CRLF//&
         'Please only use one subregion when doing a QM/MM calculation.')
       call utils_assert(.not.pub_foe, &
         'FOE is not yet compatible with embedding calculations!'//CRLF//&
         'Please only use one subregion when doing a FOE calculation.')
       call utils_assert(.not.pub_do_tddft, &
         'TDDFT is not yet compatible with embedding calculations!'//CRLF//&
         'Please only use one subregion when doing a TDDFT calculation.')
       call utils_assert(.not.pub_lr_phonons_calculate, &
         'LR_PHONONS is not yet compatible with embedding calculations!'//CRLF//&
         'Please only use one subregion when doing a LR_PHONONS calculation.')
       call utils_assert(.not.pub_xc_ke_density_required, &
         'Meta-GGAs are not yet compatible with embedding calculations!'//CRLF//&
         'Please only use one subregion when doing a meta-GGA calculation.')
    end if

    ! rc2013: if needed print the information for EMFT
    if(pub_on_root .and. (pub_emft .or. pub_emft_follow .or. pub_do_fandt)) &
         call internal_print_embedding_info

    ! qoh: Initialise pub_usehfx and pub_hfxfraction: this needs to be done before
    ! qoh: mdl%fftbox is initialised, but cannot be done in rundat.
    call xc_hfxinit

    ! jd: Spherical wave expansion module in use if HFx or DMA in use.
    call rundat_set_pub_use_swx(pub_use_hfx .or. pub_dma_calculate)

    ! jcap: Same for active subregion
    call rundat_set_pub_use_activeswx(pub_use_activehfx)

    ! jd: Remote NGWF facility currently only used if SWx is in use
    ! rc2013: or if using SWx in active region
    pub_use_remote_ngwfs = pub_use_swx .or. pub_use_activeswx

    ! ndmh: find biggest NGWF radius, including conduction NGWFs if relevant
    max_rad = maxval(mdl%elements(:)%radius)
    if (pub_cond_calculate_any_task) then
       max_rad = max(max_rad, maxval(mdl%elements(:)%radius_cond))
    end if

    ! cks: determine FFTbox dimensions in gridpoints
    call fftbox_find_size(n1,n2,n3,max_rad,3,pub_extend_ngwf,pub_ngwf_halo, &
         pub_fftbox_pref,mdl%cell%a1,mdl%cell%a2,mdl%cell%a3,mdl%cell%b1, &
         mdl%cell%b2,mdl%cell%b3,mdl%cell%d1,mdl%cell%d2,mdl%cell%d3, &
         mdl%cell%total_pt1,mdl%cell%total_pt2,mdl%cell%total_pt3)
    ! cks: initialise mdl%fftbox structure
    call fftbox_init(mdl%fftbox,n1,n2,n3,mdl%cell)

    call services_flush

    ! cks: print psinc grid info and check fftbox and ppd sizes
    call ppd_strategy_check_and_print(mdl%fftbox,mdl%cell)

    ! ndmh: initialise standard grid and standard grid FFTs
    call cell_grid_distribute(mdl%std_grid,mdl%cell,1.0_DP,'Coarse grid',.true.)
    call fourier_init_cell(mdl%std_grid)

    ! Determine whether we need to constrain the fine grid to satisfy the
    ! multigrid solver
    if (pub_multigrid_in_use) then
       constrain_for_mg = pub_multigrid_in_use.and.any(pub_multigrid_bc_is_periodic)
    else
       constrain_for_mg = .false.
    end if

    ! ndmh: set up double grid
    ! ndmh: first, see if dbl grid must be initialised separately:
    ! ndmh: if dbl_grid_scale is 1.0, then double grid is same as fine grid
    ! ndmh: otherwise it is a different grid.
    if (pub_dbl_grid_scale > 1.0_DP) then

       ! ndmh: initialise dbl grid and dbl grid FFTs
       ! jd: pass .true. as the last param if double is fine,
       !     to (potentially) tweak the distribution for multigrid
       ! JCW: If we are constraining the fine grid for multigrid, then
       ! JCW: is_fine must be .false., even if pub_fine_grid_scale ==
       ! JCW: pub_dbl_grid_scale, since the fine grid scale must be
       ! JCW: modified to satisfy the constraints of the multigrid solver
       call cell_grid_distribute(mdl%dbl_grid,mdl%cell,&
            pub_dbl_grid_scale,'Double grid',.true., &
            is_fine=(pub_fine_grid_scale==pub_dbl_grid_scale)&
            .and.(.not.constrain_for_mg))
       call fourier_init_cell(mdl%dbl_grid)
       pub_dbl_is_std = .false.

    else

       ! ndmh: copy dbl and fine grids from std grid
       mdl%dbl_grid = mdl%std_grid
       pub_dbl_is_std = .true.

    end if

    ! ndmh: set up fine grid (optionally)
    ! ndmh: first, see if fine grid must be initialised separately:
    ! ndmh: if fine_grid_scale differs from dbl_grid scale, it will
    ! ndmh: need to be counted as a separate grid.
    ! JCW:  The fine grid should also be initialised separately if the
    ! JCW:  multigrid solver is in use with periodic BCs, since we need
    ! JCW:  to ensure that it satisfies the grid constraints.
    if ( pub_fine_grid_scale/=pub_dbl_grid_scale.or.&
         constrain_for_mg ) then

       ! Only constrain the fine grid for multigrid solver if
       ! we are using the multigrid solver and if the fine grid is
       ! periodic

       ! ndmh: initialise fine grid
       call cell_grid_distribute(mdl%fine_grid,mdl%cell, &
            pub_fine_grid_scale,'  Fine grid',.true.,&
            is_fine=.true.,&
            constrain_for_mg=constrain_for_mg)
       call fourier_init_cell(mdl%fine_grid)
       pub_fine_is_dbl = .false.

    else if ( utils_devel_code(.false.,"MG","FORCE_PBC_GRID_CONSTRAINT",&
         pub_devel_code) ) then
       ! JCW: Force constraint of the fine grid for full PBCs, so we can
       ! JCW: obtain identical fine grid sizes to those used in multigrid
       ! JCW: calculations when not using the solver.
       call cell_grid_distribute(mdl%fine_grid,mdl%cell, &
            pub_fine_grid_scale,'  Fine grid',.true.,&
            is_fine=.true.,&
            constrain_for_mg=.true.,&
            mg_is_periodic = [.true.,.true.,.true.] )
       call fourier_init_cell(mdl%fine_grid)
       pub_fine_is_dbl = .false.

    else

       ! ndmh: copy dbl and fine grids from std grid
       mdl%fine_grid = mdl%dbl_grid
       pub_fine_is_dbl = .true.

    end if

    ! jd: Initialise multigrid module, if required
    if (pub_is_implicit_solvent .or. pub_is_smeared_ion_rep) then
       call multigrid_initialise(mdl%fine_grid)
    end if

    ! jd: If local openbc pseudo is in use, initialise it
    ! JCW: This is the case when pub_pspot_bc_is_periodic is .false. along
    ! JCW: each lattice direction
    if (.not.pub_coulomb_cutoff) then
       ! JCW: Only use openbc_locps when not performing a cutoff Coulomb
       ! JCW: calculation. The cutoff_coulomb module deals with the local
       ! JCW: pseudopotential itself via cutoff_coulomb_localpseudo.
       if (.not.any(pub_pspot_bc_is_periodic)) then
          ! rc2013: do this for all regions
          do isub=1,mdl%nsub
             call openbc_locps_init(mdl%cell, mdl%regions(isub)%pseudo_sp, &
                  mdl%regions(isub)%paw_sp,mdl%regions(isub)%par%num_pspecies)
          end do
       end if
    end if

    ! ndmh: set up parallel strategy and Fourier transforms for padded cell
    if (pub_coulomb_cutoff) call cutoff_coulomb_init(mdl)

    ! rc2013: update mdl data with what we've initialised in regions so far
    ! This is needed here so that spin can be setup correctly
    do isub=1,mdl%nsub
       do iat=1,mdl%regions(isub)%par%nat
          glob_iat = mdl%regions(isub)%elements(iat)%global_atom_number
          mdl%elements(glob_iat) = mdl%regions(isub)%elements(iat)
       end do
    end do

    ! ndmh: determine spin polarisation strategy
    call internal_setup_spin

    ! ndmh: Initialise xc module
    call xc_init

    ! ndmh: Initialise PAW core xc terms
    !if (pub_paw) call paw_exc_core_init(mdl%paw_sp)
    ! rc2013: should translate this data to regions somehow...
    do isub=1,mdl%nsub
       if (pub_paw) call paw_exc_core_init(mdl%regions(isub)%paw_sp)
    end do

    ! ndmh: Allocate storage for radial densities
    ! rc2013: initialise subsystem density storage too
    !if (pub_initial_dens_realspace) call density_init_radial_init(mdl)
    if (pub_initial_dens_realspace) then
       do isub=1,mdl%nsub
          ! rc2013: also should check mdl%nsp...
          call density_init_radial_init(mdl%regions(isub))
       end do
    end if

    ! rc2013: if doing EMFT as a follow-up to a regular optimisation, then
    ! initialise everything at the lower level of theory
    if(pub_emft_follow) pub_emft = .false.
    ! Similarly if doing EMFT with block orthogonalisation then we don't want
    ! to use it during the low-level optimisation
    if(pub_block_orthogonalise .and. .not. pub_emft_follow) pub_build_bo=.true.

    ! rc2013: initialise radial NGWF functions for all subsystems
    do isub=1,mdl%nsub
       ! rc2013: if we're in the active region then we need to initialise
       ! the correct XC functional (in principle -- this may be unnecessary)
       if(pub_emft .and. (isub == pub_active_region)) &
            call xc_embed_swap_functional(.true.)

       ! ndmh: initialise radial functions for NGWFs for each species
       ! jcap: state which region's NGWFs we are initialising first
       if (mdl%nsub.gt.1) then
          write(isub_str,'(i0)')isub
          if(pub_on_root) write(stdout,'(a)') &
               utils_banner('<', 'Region '//isub_str)
       end if
       call ngwfs_radial_init('ngwfs               ', &
            mdl%regions(isub),mdl%dftb_par)
       if (pub_cond_calculate_any_task) call ngwfs_radial_init(&
            'ngwfs_cond          ', mdl%regions(isub),mdl%dftb_par)
       if (pub_use_aux_ngwfs) call ngwfs_radial_init('ngwfs_aux           ', &
            mdl%regions(isub),mdl%dftb_par)
       if (pub_hubbard_atomsolve) call ngwfs_radial_init('hub_ngwfs           ', &
            mdl%regions(isub),mdl%dftb_par)

       ! rc2013: copy info from regions%elements to global list
       do iat=1,mdl%regions(isub)%par%nat
          glob_iat = mdl%regions(isub)%elements(iat)%global_atom_number
          mdl%elements(glob_iat) = mdl%regions(isub)%elements(iat)
       end do

       ! rc2013: now let's swap the functionals back
       if(pub_emft .and. isub == pub_active_region) &
            call xc_embed_swap_functional(.false.)
    end do

    ! ndmh: initialise dense matrix module
    call dense_init

    if (pub_aug) then
       ! ndmh: set up num of components of PAW compensation density & gradient
       pub_aug_den_dim = 0
       if (pub_xc_gradient_corrected) pub_aug_den_dim = 3
    end if

    ! agrecokpt: generate k-points list from pub_kpoint_list
    ! only if k-points are specified in the input file
    if (pub_kpoints_specified_list) then
       call kpoints_generate_unique_list(mdl%unique_kpoints, &
            mdl%nkpoints, pub_kpoint_list, mdl%cell)
    ! agrecokpt: initialise from k-point mesh
    ! checks to avoid specification of both mesh and list are in rundat module
    ! currently using is_time_reversal=1 only
    else if (any(pub_kp_grid_size /= 1).or.any(pub_kp_grid_shift /= 0)) then
       call kpoints_generate_unique_list(mdl%unique_kpoints, &
            mdl%nkpoints, pub_kp_grid_size, pub_kp_grid_shift, &
            1, mdl%cell, mdl%elements, mdl%nat)

    ! agrecokpt: standard Gamma-point only sampling
    else
       mdl%nkpoints = 1
    end if

    call timer_clock(myself,2)

    return

contains

    subroutine internal_setup_spin

      !======================================================================!
      ! This suboutine counts the number of electrons in the cell and uses   !
      ! this to check the spin polarisation strategy.                        !
      !----------------------------------------------------------------------!
      ! Written by Peter Haynes, July 2006, based on some code from          !
      ! internal_ppd_sph_box_kb_init by Chris-Kriton Skylaris.               !
      ! Modified by Phil Avraam, November 2008 for fractional electronic     !
      ! charges (tidied up by Nicholas Hine)                                 !
      !======================================================================!

      use constants, only: stdout, garbage_int
      use rundat, only: pub_spin, pub_spin_polarised, pub_num_spins, &
           pub_spin_fac, pub_lr_tddft_calculate, pub_edft, pub_real_spin

      implicit none

      ! Local variables
      integer :: nelecs
      integer :: nspins

      ! ndmh: count number of electrons
      nelecs = nint(energy_and_force_count_elecs(mdl%elements, mdl%nat))

      ! kkbd: If the user does not specify a spin then apply the default:
      !         even nelecs -> net spin 0
      !         odd nelecs  -> net spin 1
      !       This is because in EDFT, it is perfectly valid to have
      !       no net spin when we have an odd number of electrons but
      !       this is not the behavior a user would expect by default.
      !       There must therefore be a different behaviour if the user
      !       specifies spin 0 vs not specifying a spin.
      if (pub_spin == garbage_int) then
         if (mod(nelecs, 2) == 0) then
            pub_spin = 0
         else
            pub_spin = 1
         end if
         pub_real_spin = real(pub_spin, kind=dp)
      end if

      ! Make sure spin polarization is enabled if we need it
      if ((.not.pub_edft) .and. (abs(pub_spin) > 0)) then
         pub_spin_polarised = .true.
      else if (pub_edft .and. (abs(pub_real_spin) > 0.000000001)) then
         pub_spin_polarised = .true.
      end if

      ! pdh: check for inconsistencies between charge and spin
      ! pdh: and set number of occupied states accordingly
      if ((.not. pub_edft) .and. (mod(nelecs,2) /= 0)) then
         if (.not. pub_spin_polarised) then
            if (pub_on_root) write(stdout,'(a/a)') myself//&
                 ': WARNING in internal_setup_spin:',&
                 '        odd number of electrons in cell - &
                 &spin-polarised calculation will be performed'
            pub_spin_polarised = .true.
         end if
         if (pub_spin == 0) pub_spin = 1
         if (mod(pub_spin,2) == 0) then
            if (pub_on_root) write(stdout,'(a/a)') myself//&
                 ': WARNING in internal_setup_spin:',&
                 '        odd number of electrons in cell - &
                 &even spin has been over-ridden to unity'
            pub_spin = 1
         end if
      else if (.not. pub_edft) then
         if (pub_spin_polarised .and. pub_spin == 0) then
            if (pub_on_root) write(stdout,'(a/a/a)') myself//&
                 ': WARNING in internal_setup_spin:',&
                 '        even number of electrons in cell and zero spin', &
                 '        but continuing with spin-polarised calculation'
         end if
         if (mod(pub_spin,2) /= 0) then
            if (pub_on_root) write(stdout,'(a/a)') myself//&
                 ': WARNING in internal_setup_spin:',&
                 '        even number of electrons in cell - &
                 &odd spin has been over-ridden to zero'
            pub_spin = 0
         end if
      end if

      ! pdh: set number of spins
      nspins = 1
      if (pub_spin_polarised) nspins = 2
      pub_num_spins = nspins
      pub_spin_fac = 2.0_DP / real(pub_num_spins,kind=DP)

    end subroutine internal_setup_spin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine internal_print_num_species(mdl)

    !=========================================================-!
    ! This subroutine prints out the number of atoms and NGWFs !
    ! for each species.                                        !
    !----------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/03/2007.          !
    ! Conduction NGWFs added by Nicholas Hine on 12/10/2011.   !
    ! Subsystem info added by Robert Charlton, 07/08/2018.     !
    !==========================================================!

    use constants, only: stdout
    use model_type, only: MODEL
    use rundat, only: pub_cond_calculate_any_task, pub_paw
    use utils, only: utils_banner

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl

    ! Local Variables
    integer :: n_atoms
    integer :: n_ngwfs
    integer :: n_ngwfs_cond
    integer :: n_projs
    integer :: species_counter
    integer :: row
    character(len=2) :: el_symbol
    integer :: isub
    ! rc2013: totals for all systems
    integer :: n_tot_atoms, n_tot_ngwfs, n_tot_ngwfs_cond, n_tot_projs

    write(stdout,'(a)') utils_banner('-','Atom counting information')

    ! rc2013: print info for all subsystems, then totals at end
    n_tot_atoms = 0
    n_tot_ngwfs = 0
    n_tot_ngwfs_cond = 0
    n_tot_projs = 0
    regions_loop: do isub=1,mdl%nsub
       if(mdl%nsub .gt. 1) then
          write(stdout,'(a)') repeat('=',50)
          write(stdout,'(a28, i3, a3)') 'Printing info for subsystem ', isub, '...'
          write(stdout,'(a)') repeat('=',50)
       end if
       if (pub_cond_calculate_any_task) then
          write(stdout,"(a)")"Symbol    Natoms  Nvalngwfs Ncondngwfs  Nprojs"
       else
          write(stdout,"(a)")"Symbol    Natoms    Nngwfs    Nprojs"
       end if

       do species_counter=1,mdl%regions(isub)%par%num_pspecies

          n_atoms = 0
          n_ngwfs = 0
          n_ngwfs_cond = 0
          n_projs = 0
          do row=1,mdl%regions(isub)%par%nat

             if (mdl%regions(isub)%elements(row)%pspecies_number &
                  == species_counter) then
                n_atoms = n_atoms + 1
                n_ngwfs = n_ngwfs + mdl%regions(isub)%elements(row)%nfunctions
                n_ngwfs_cond = n_ngwfs_cond + &
                     mdl%regions(isub)%elements(row)%nfunctions_cond
                if (pub_paw) then
                   n_projs = n_projs + mdl%regions(isub)%elements(row)%npawpws
                else
                   n_projs = n_projs + &
                        mdl%regions(isub)%elements(row)%nprojectors
                end if
                el_symbol = mdl%regions(isub)%elements(row)%symbol
             endif
          enddo

          n_tot_atoms = n_tot_atoms + n_atoms
          n_tot_ngwfs = n_tot_ngwfs + n_ngwfs
          n_tot_ngwfs_cond = n_tot_ngwfs_cond + n_ngwfs_cond
          n_tot_projs = n_tot_projs + n_projs

          if (pub_cond_calculate_any_task) then
              write(stdout,"(a4,4i10)")el_symbol,n_atoms,n_ngwfs,n_ngwfs_cond,n_projs
          else
              write(stdout,"(a4,3i10)")el_symbol,n_atoms,n_ngwfs,n_projs
          end if

       enddo

       if (pub_paw) then
           n_projs = mdl%regions(isub)%par%num_pawpws
       else
           n_projs = mdl%regions(isub)%par%num_projectors
       end if
       if (pub_cond_calculate_any_task) then
          write(stdout,"(a)")".......   ......    ......    ......    ......"
          write(stdout,"(a7,i7,3i10)")"Totals:",mdl%regions(isub)%par%nat, &
               mdl%regions(isub)%par%num_ngwfs, &
               mdl%regions(isub)%par%num_ngwfs_cond, n_projs
       else
          write(stdout,"(a)")".......   ......    ......    ......"
          write(stdout,"(a7,i7,2i10)")"Totals:",mdl%regions(isub)%par%nat, &
               mdl%regions(isub)%par%num_ngwfs, n_projs
       end if
    enddo regions_loop
    if(mdl%nsub .gt. 1) then
       write(stdout,'(a)') repeat('=',50)
       write(stdout,'(a28)') 'Results for all subsystems: '
       if (pub_cond_calculate_any_task) then
          write(stdout,"(a7,i7,3i10)")"Totals:", n_tot_atoms, n_tot_ngwfs, &
               n_tot_ngwfs_cond, n_tot_projs
       else
          write(stdout,"(a7,i7,2i10)")"Totals:", n_tot_atoms, n_tot_ngwfs, &
               n_tot_projs
       endif
    end if
    write(stdout,'(a/)') repeat('-',80)

  end subroutine internal_print_num_species

    subroutine internal_print_embedding_info

      !============================================================!
      ! This subroutine prints out information for embedding runs. !
      !------------------------------------------------------------!
      ! Written by Robert Charlton, 16/11/2018.                    !
      !============================================================!

      use rundat, only: pub_active_region, pub_ngwf_regions_ngroups, pub_do_fandt, &
           pub_xc_functional, pub_active_xc_functional, pub_freeze_envir_ngwfs, &
           pub_emft, pub_emft_follow, pub_emft_lnv_only, &
           pub_block_orthogonalise
      use utils, only: utils_banner

      implicit none

      ! Local variables
      integer :: isub, isp

      write(stdout,'(a)') utils_banner('=','Embedding information')
      write(stdout,'(i1,a)') pub_ngwf_regions_ngroups, ' subsystems specified...'
      if(pub_do_fandt) then
         write(stdout,'(a)') '... freeze-and-thaw (F+T) requested'
      else if(pub_freeze_envir_ngwfs) then
         write(stdout,'(a)') '... frozen density embedding (FDE) requested'
      else
         write(stdout,'(a)') '... all NGWFs will be optimised'
      end if
      if(pub_block_orthogonalise) then
         write(stdout,'(a)') '... block orthogonalisation requested'
      end if
      if(pub_emft) then
         write(stdout,'(a)') '... Embedded mean-field theory (EMFT) requested'
         ! rc2013: state how the optimisation will be done
         if (pub_emft_follow .and. pub_emft_lnv_only) then
            write(stdout,'(a)') '... EMFT optimisation of density kernel only &
                 &in active region will be done after ' //CRLF// &
                 '    a low-level optimisation of the full system'
         else if (pub_emft_follow) then
            write(stdout,'(a)') '... EMFT optimisation of NGWFs and kernel &
                 &in active region will be done after ' //CRLF// &
                 '    a low-level optimisation of the full system'
         else
            write(stdout,'(a)') '... full EMFT optimisation of NGWFs and kernel'
         end if
      else
         write(stdout,'(a)') '... 1 level of theory will be used for all &
              &regions'
      end if
      do isub=1,pub_ngwf_regions_ngroups
         write(stdout,'(a)') repeat('-',50)
         write(stdout,'(a11,i1,a,i4,a)',advance='no') '    Region ', isub, &
              ' contains ', mdl%regions(isub)%par%nat, ' atoms from the &
              &following species: '
         do isp=1,size(mdl%regions(isub)%species)
            write(stdout,'(a)',advance='no') &
                 mdl%regions(isub)%species(isp)%species_id
         end do
         if(isub == pub_active_region) then
            write(stdout,'(/3a)') '    These will be evaluated at the ',&
                 adjustl(trim(pub_active_xc_functional)), ' level.'
         else
            write(stdout,'(/3a)') '    These will be evaluated at the ',&
                 adjustl(trim(pub_xc_functional)), ' level.'
         end if
      end do
      write(stdout,'(a/)') repeat('=',80)

    end subroutine internal_print_embedding_info

  end subroutine energy_and_force_init_cell

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  real(kind=dp) function energy_and_force_count_elecs(elements, nat)

    use constants, only: DP
    use ion, only: ELEMENT
    use rundat, only: pub_charge, pub_edft
    use utils, only: utils_abort

    ! Arguments
    integer,       intent(in)   :: nat
    type(ELEMENT), intent(in)   :: elements(nat)

    ! Local variables
    integer,parameter :: cfac=4
    real(kind=DP),parameter :: ctol=0.00001_DP
    integer :: atom
    real(kind=DP) :: cfac_nelecs !, frac_nelecs ! ndmh: removed duplicate

    ! Count number of electrons for charge neutral cell
    cfac_nelecs = 0.0_DP
    do atom=1,nat
       ! pa: number of electrons times common factor
       ! pa: used to round the ionic charge to the nearest 1/cfac
       cfac_nelecs = cfac_nelecs + &
            real(nint(cfac*(elements(atom)%ion_charge)), kind=DP)
    end do
    ! pa: check that the total number of electrons is integer
    ! kkbd: only do this if we're NOT doing an edft run
    if (.not.pub_edft .and. (abs(cfac_nelecs/real(cfac,DP) - pub_charge - &
         nint(cfac_nelecs/real(cfac,DP) - pub_charge)) > ctol)) then
       call utils_abort('Error in energy_and_force_count_elecs: total &
            &electronic charge is non-integer.')
    endif

    ! cks: add the total charge of the system (from the input file)
    ! pa: allow for fractional ionic charge but not fractional total charge
    ! kkbd: Allow fractional total charge, preserve previous behavior
    energy_and_force_count_elecs = cfac_nelecs/cfac - pub_charge

  end function energy_and_force_count_elecs

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine energy_and_force_exit_cell(mdl)

    use cell_grid, only: cell_grid_exit
    use cutoff_coulomb, only: cutoff_coulomb_exit
    use dense, only: dense_exit
    use density_init, only: density_init_radial_exit
    use fft_box, only: fftbox_exit
    use fourier, only: fourier_exit_cell, fourier_exit_threads
    use is_solvation, only: implicit_solvent_exit
    use is_solvation_boltzmann, only: implicit_solvent_boltzmann_exit_ions
    use model_type, only: MODEL
    use multigrid_methods, only: multigrid_exit
    use ngwfs_radial, only: ngwfs_radial_exit
    use openbc_locps, only: openbc_locps_exit
    use paw, only: paw_all_species_exit
    use pbc_corrections, only: pbc_corr_exit
    use pseudopotentials, only: pseudopotentials_species_exit
    use rundat, only: pub_coulomb_cutoff, pub_is_implicit_solvent, &
         pub_fine_is_dbl, pub_dbl_is_std, pub_paw, &
         pub_pspot_bc_is_periodic, &
         pub_is_smeared_ion_rep, pub_is_pbe, pub_initial_dens_realspace
    use utils, only: utils_dealloc_check
    use xc, only: xc_exit

    implicit none

    ! Arguments
    type(MODEL), intent(inout) :: mdl

    ! Local Variables
    integer :: ierr
    integer :: isub
    character(len=*), parameter :: myself = 'energy_and_force_exit_cell'

    ! ndmh: deallocate storage for radial NGWFs
    ! rc2013: deallocate each region
    do isub=1,mdl%nsub
       call ngwfs_radial_exit(mdl%regions(isub))

       ! ndmh: deallocate storage for radial densities
       if (pub_initial_dens_realspace) &
            call density_init_radial_exit(mdl%regions(isub))
    enddo

    ! clean up dense module
    call dense_exit

    ! clean up XC routines
    call xc_exit

    ! ndmh: clean up cutoff coulomb module
    if (pub_coulomb_cutoff) call cutoff_coulomb_exit(mdl)

    ! jd: clean up pbc_correction module
    call pbc_corr_exit

    ! jd: clean up Boltzmann solvation (global state)
    if (pub_is_pbe /= 'NONE') then
       call implicit_solvent_boltzmann_exit_ions
    end if

    ! jd: clean up is_solvation module
    if (pub_is_implicit_solvent) call implicit_solvent_exit

    ! jd: clean up multigrid module
    if (pub_is_implicit_solvent .or. pub_is_smeared_ion_rep) then
       call multigrid_exit
    end if

    ! jd: clean up openbc local pseudo, if in use
    ! JCW: This is the case when pub_pspot_bc_is_periodic is .false. along
    ! JCW: each lattice direction
    if (.not.pub_coulomb_cutoff) then
       ! JCW: Only use openbc_locps when not performing a cutoff Coulomb
       ! JCW: calculation. The cutoff_coulomb module deals with the local
       ! JCW: pseudopotential itself via cutoff_coulomb_localpseudo.
       if (.not.any(pub_pspot_bc_is_periodic)) then
          do isub=1,mdl%nsub
             call openbc_locps_exit(mdl%regions(isub)%pseudo_sp, &
                  mdl%regions(isub)%paw_sp, mdl%regions(isub)%par%num_pspecies)
          enddo
       end if
    end if

    ! clean up whole-cell FFT routines
    call fourier_exit_cell

    ! clean up cell grids
    if (.not.pub_fine_is_dbl) call cell_grid_exit(mdl%fine_grid)
    if (.not.pub_dbl_is_std) call cell_grid_exit(mdl%dbl_grid)
    call cell_grid_exit(mdl%std_grid)

    ! deallocate arrays in cell
    call fftbox_exit(mdl%fftbox)

    ! clean up arrays in the PAW/pseudopotentials modules
    if (.not.pub_paw) then
       call pseudopotentials_species_exit(mdl%pseudo_sp)
       deallocate(mdl%paw_sp,stat=ierr)
       call utils_dealloc_check(myself,'mdl%paw_sp',ierr)
    else
       call paw_all_species_exit(mdl%paw_sp)
       deallocate(mdl%pseudo_sp,stat=ierr)
       call utils_dealloc_check(myself, 'mdl%pseudo_sp',ierr)
    end if

    ! rc2013: repeat for the subsystem quantities
    ! jcap: if we are in a multiple region calculation, as otherwise
    ! these will be pointers to the overall system arrays
    if (mdl%nsub.gt.1) then
       do isub=1,mdl%nsub
          if (.not.pub_paw) then
             call pseudopotentials_species_exit(mdl%regions(isub)%pseudo_sp)
             deallocate(mdl%regions(isub)%paw_sp,stat=ierr)
             call utils_dealloc_check(myself,'mdl%regions(isub)%paw_sp',ierr)
          else
             call paw_all_species_exit(mdl%regions(isub)%paw_sp)
             deallocate(mdl%regions(isub)%pseudo_sp,stat=ierr)
             call utils_dealloc_check(myself,'mdl%regions(isub)%pseudo_sp',ierr)
          end if
       enddo
    end if

    ! agrecokpt: deallocate k-point array if required
    if (allocated(mdl%unique_kpoints)) then
       deallocate(mdl%unique_kpoints,stat=ierr)
        call utils_dealloc_check('generate_unique_kpoints_list', &
             'mdl%unique_kpoints', ierr)
    end if

    ! close down FFT threading
    call fourier_exit_threads

  end subroutine energy_and_force_exit_cell

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine energy_and_force_calculate(total_energy,total_forces,mdl, &
       properties_only,return_converged, is_cmplx, kpt, custom_restart_name)
    !==========================================================================!
    ! Wrapper to an energy and force calculation.                              !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Sep 2022.                                 !
    !==========================================================================!

    use constants, only: DP
    use model_type, only: MODEL
    use rundat, only: pub_dftb, pub_is_auto_solvation, pub_cond_calculate, &
         pub_lr_tddft_calculate

    implicit none

    ! ab: arguments
    type(MODEL), target, intent(inout)  :: mdl
    real(kind=DP), intent(out)  :: total_energy
    real(kind=DP), intent(out)  :: total_forces(1:3,1:mdl%nat)
    logical,intent(in),optional :: properties_only
    logical,intent(out),optional :: return_converged
    ! agrecocmplx: optional argument for complex NGWFs
    logical, intent(in), optional :: is_cmplx
    ! agrecokpt: optional argument for single kpt, in fractional coordinates
    real(kind=DP), intent(in), optional :: kpt(3)
    character(len=*), optional, intent(in) :: custom_restart_name

    ! -----------------------------------------------------------------------

    ! jd: Should a vacuum calculation automatically precede solvation?
    ! ndmh: deactivate vacuum calculation if doing COND or LR-TDDFT, as
    ! ndmh: it will have already been done
    if ((pub_is_auto_solvation).and. &
         (.not.pub_cond_calculate).and.(.not.pub_lr_tddft_calculate)) then
       call energy_and_forces_auto_solvation(total_energy,total_forces,mdl, &
       properties_only,return_converged, is_cmplx, kpt, custom_restart_name)
    else
       call energy_and_forces(total_energy,total_forces,mdl, &
       properties_only,return_converged, is_cmplx, kpt, custom_restart_name)
    end if

  end subroutine energy_and_force_calculate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine energy_and_forces_auto_solvation(total_energy,total_forces,mdl, &
       properties_only,return_converged, is_cmplx, kpt, custom_restart_name)
    !==========================================================================!
    ! Wrapper to an energy and force calculation for auto solvation.           !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Sep 2022.                                 !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout, CRLF, HARTREE_TO_KCAL_PER_MOL
    use model_type, only: MODEL
    use is_solvation, only: implicit_solvent_energy_terms
    use is_solvation_boltzmann, only: implicit_solvent_boltzmann_energy_contrib
    use restart, only: restart_read_filename_prefix, restart_write_filename_prefix
    use rundat, only: pub_is_auto_solvation, pub_is_implicit_solvent,&
         pub_is_include_apolar, pub_is_pbe, pub_is_pbe_bc_debye_screening, &
         pub_is_pbe_neutralisation_scheme, pub_edft, &
         pub_multigrid_bc_is_periodic, pub_write_denskern, pub_task, &
         pub_write_tightbox_ngwfs, pub_write_hamiltonian, pub_read_hamiltonian, &
         pub_is_restart_vac_from_vac, pub_read_denskern, pub_maxit_pen, &
         pub_read_tightbox_ngwfs, pub_is_apolar_method

    implicit none

    ! ab: arguments
    type(MODEL), target, intent(inout)  :: mdl
    real(kind=DP), intent(out)  :: total_energy
    real(kind=DP), intent(out)  :: total_forces(1:3,1:mdl%nat)
    logical,intent(in),optional :: properties_only
    logical,intent(out),optional :: return_converged
    ! agrecocmplx: optional argument for complex NGWFs
    logical, intent(in), optional :: is_cmplx
    ! agrecokpt: optional argument for single kpt, in fractional coordinates
    real(kind=DP), intent(in), optional :: kpt(3)
    character(len=*), optional, intent(in) :: custom_restart_name

    ! -----------------------------------------------------------------------

    ! jd: Stage of auto-solvation:
    !     V - in vacuum, S - in solvent, X - no auto solvation
    character :: auto_solvation_state = 'X'
    character(len=50) :: change_string
    character(len=50) :: empty_string=repeat(' ',50)
    character(len=80) :: orig_is_pbe, orig_is_pbe_neutralisation_scheme
    integer :: orig_maxit_pen
    logical :: ngwfs_converged, orig_is_include_apolar, orig_write_tightbox_ngwfs, &
         orig_write_denskern, orig_write_hamiltonian, orig_is_implicit_solvent, &
         orig_is_pbe_bc_debye_screening
    real(kind=DP) :: total_energy_in_vacuum ! jd: Needed in auto solvation
    real(kind=DP) :: nonelectrostatic_dft_energy
    real(kind=DP) :: nonelectrostatic_dft_energy_in_vacuum
    real(kind=DP) :: electrostatic_energy
    real(kind=DP) :: electrostatic_energy_in_vacuum
    real(kind=DP) :: apolar_energy
    real(kind=DP) :: polar_energy
    real(kind=DP) :: solvation_energy
    real(kind=DP) :: solvation_energy_check

    ! ab: calculation in vacuum
    auto_solvation_state = 'V'
    orig_is_implicit_solvent = pub_is_implicit_solvent
    pub_is_implicit_solvent = .false.
    orig_is_include_apolar = pub_is_include_apolar
    pub_is_include_apolar = .false.
    orig_is_pbe = pub_is_pbe     ! jd: We do not want Boltzmann ions
    pub_is_pbe = 'NONE'          !     in the vacuum calculation.
    orig_is_pbe_bc_debye_screening = pub_is_pbe_bc_debye_screening
    pub_is_pbe_bc_debye_screening = .false. ! jd: ... or Debye screening
    orig_is_pbe_neutralisation_scheme = pub_is_pbe_neutralisation_scheme
    if(all(pub_multigrid_bc_is_periodic(:))) then
      ! jd: For PBCs
      pub_is_pbe_neutralisation_scheme = 'JELLIUM'
    else
      ! jd: for OBCs
      pub_is_pbe_neutralisation_scheme = 'NONE'
    end if
    orig_write_denskern = pub_write_denskern
    pub_write_denskern = .true.
    orig_write_tightbox_ngwfs = pub_write_tightbox_ngwfs
    pub_write_tightbox_ngwfs = .true.
    orig_write_hamiltonian = pub_write_hamiltonian
    if(pub_edft) pub_write_hamiltonian = .true.
    ! ebl: altering geometry optimisation to read from previous vacuum run, not solvated run
    ! jcap: or if asked to in the input file
    if (pub_task=='GEOMETRYOPTIMIZATION'.or.pub_is_restart_vac_from_vac) then
       restart_read_filename_prefix = 'vacuum_'
    else
       restart_read_filename_prefix = ''
    end if
    restart_write_filename_prefix = 'vacuum_'
    if(pub_on_root) then
       write(stdout,'(a)') CRLF//'================&
            &========== Calculation in vacuum =============================='
       write(stdout,'(a)') '= Adjusting the following parameters:'
       write(change_string,'(a,l1,a)') &
            ' (changed from ', orig_is_implicit_solvent, ' set already in &
            &rundat)'
       write(stdout,'(a,l1,a)') '= is_implicit_solvent: ', &
            pub_is_implicit_solvent, &
            merge(change_string, empty_string, &
            pub_is_implicit_solvent .neqv. orig_is_implicit_solvent)
       write(change_string,'(a,l1,a)') &
            ' (changed from ', orig_is_include_apolar, ')'
       write(stdout,'(a,l1,a)') '= is_include_apolar: ', &
            pub_is_include_apolar, &
            merge(change_string, empty_string, &
            pub_is_include_apolar .neqv. orig_is_include_apolar)
       write(change_string,'(a,a,a)') &
            ' (changed from ', trim(orig_is_pbe), ')'
       write(stdout,'(a,a,a)') '= is_pbe: ', trim(pub_is_pbe), &
            merge(change_string, empty_string, &
            trim(pub_is_pbe) /= trim(orig_is_pbe))
       write(change_string,'(a,l1,a)') &
            ' (changed from ', orig_is_pbe_bc_debye_screening, ')'
       write(stdout,'(a,l1,a)') '= is_pbe_bc_debye_screening: ', &
            pub_is_pbe_bc_debye_screening, &
            merge(change_string, empty_string, &
            pub_is_pbe_bc_debye_screening .neqv. orig_is_pbe_bc_debye_screening)
       write(change_string,'(a,a,a)') &
            ' (changed from ', trim(orig_is_pbe_neutralisation_scheme), ')'
       write(stdout,'(a,a,a)') '= is_pbe_neutralisation_scheme: ', &
            trim(pub_is_pbe_neutralisation_scheme), &
            merge(change_string, empty_string, &
            trim(pub_is_pbe) /= trim(orig_is_pbe))
       write(change_string,'(a,l1,a)') &
            ' (changed from ', orig_write_denskern, ')'
       write(stdout,'(a,l1,a)') '= write_denskern: ', pub_write_denskern, &
            merge(change_string, empty_string, &
            pub_write_denskern .neqv. orig_write_denskern)
       write(change_string,'(a,l1,a)') &
            ' (changed from ', orig_write_tightbox_ngwfs, ')'
       write(stdout,'(a,l1,a)') '= write_tightbox_ngwfs: ', &
            pub_write_tightbox_ngwfs, &
            merge(change_string, empty_string, &
            pub_write_tightbox_ngwfs .neqv. orig_write_tightbox_ngwfs)
       if(pub_edft .and. orig_write_hamiltonian .neqv. pub_write_hamiltonian) then
          write(stdout,'(a,l1)') '= write_hamiltonian: ', &
               pub_write_hamiltonian
       end if
       write(stdout,'(3a)') '= restart read filename prefix: "',&
            trim(restart_read_filename_prefix),'"'
       write(stdout,'(3a)') '= restart write filename prefix: "',&
            trim(restart_write_filename_prefix),'"'
       write(stdout,'(a)') '================&
            &==============================================================='
    end if
    call energy_and_forces(total_energy,total_forces,mdl, &
    properties_only,ngwfs_converged, is_cmplx, kpt, custom_restart_name)
    total_energy_in_vacuum = total_energy
    if(pub_on_root) then
       write(stdout,'(a)') CRLF//'===================================&
            &============================================='
       write(stdout,'(a,a)') 'IS: Calculation in vacuum completed', &
            merge(' (NGWFs converged).      ', ' (*NO* NGWF convergence).', &
            ngwfs_converged)
       write(stdout,'(a)') '=========================================&
            &======================================='
    end if

    ! ab: calculation in solvent
    auto_solvation_state = 'S'
    pub_is_auto_solvation = .false.
    pub_is_implicit_solvent = .true.
    pub_is_pbe = orig_is_pbe
    pub_is_pbe_bc_debye_screening = orig_is_pbe_bc_debye_screening
    pub_is_pbe_neutralisation_scheme = orig_is_pbe_neutralisation_scheme
    pub_is_include_apolar = orig_is_include_apolar
    pub_write_denskern = orig_write_denskern
    pub_write_tightbox_ngwfs = orig_write_tightbox_ngwfs
    if(pub_edft) pub_write_hamiltonian = orig_write_hamiltonian
    pub_read_denskern = .true.
    pub_read_tightbox_ngwfs = .true.
    if(pub_edft) pub_read_hamiltonian = .true.
    orig_maxit_pen = pub_maxit_pen
    pub_maxit_pen = 0
    restart_write_filename_prefix = ''
    restart_read_filename_prefix = 'vacuum_'
    if(pub_on_root) then
       write(stdout,'(a)') CRLF//'================&
            &========== Calculation in solvent =============================='
       write(stdout,'(a)') '= Adjusting the following parameters: '
       write(stdout,'(a,l1)') '= is_implicit_solvent: ', &
            pub_is_implicit_solvent
       write(stdout,'(a,l1,a)') '= is_include_apolar: ', &
            pub_is_include_apolar,' (restored to user-specified value)'
       if(trim(orig_is_pbe) /= 'NONE') then
          write(stdout,'(a,a,a)') '= is_pbe: ', trim(pub_is_pbe), &
               ' (restored to user-specified value)'
       end if
       if(orig_is_pbe_bc_debye_screening) then
          write(stdout,'(a,l1,a)') '= is_pbe_bc_debye_screening: ', &
               orig_is_pbe_bc_debye_screening, &
               ' (restored to user-specified value)'
       end if
       if(trim(orig_is_pbe_neutralisation_scheme) /= &
            trim(pub_is_pbe_neutralisation_scheme)) then
          write(stdout,'(a,a,a)') '= is_pbe_neutralisation_scheme: ', &
               trim(pub_is_pbe_neutralisation_scheme), &
               ' (restored to user-specified value)'
       end if
       write(stdout,'(a,l1,a)') '= write_denskern: ', &
            pub_write_denskern,' (restored to user-specified value)'
       write(stdout,'(a,l1,a)') '= write_tightbox_ngwfs: ', &
            pub_write_tightbox_ngwfs,' (restored to user-specified value)'
       if(pub_edft) then
          write(stdout,'(a,l1,a)') '= write_hamiltonian: ', &
               pub_write_hamiltonian, ' (restored to user-specified value)'
       end if
       write(stdout,'(a,l1)') '= read_denskern: ', pub_read_denskern
       write(stdout,'(a,l1)') '= read_tightbox_ngwfs: ', pub_read_tightbox_ngwfs
       if(pub_edft) then
          write(stdout,'(a,l1)') '= read_hamiltonian: ', pub_read_hamiltonian
       end if
       write(stdout,'(a,i0)') '= maxit_pen: ', pub_maxit_pen
       write(stdout,'(3a)') '= restart read filename prefix: "',&
            trim(restart_read_filename_prefix), &
            '"'
       write(stdout,'(3a)') '= restart write filename prefix: "',&
            trim(restart_write_filename_prefix), &
            '"'
       write(stdout,'(a)') '========================&
            &========================================================'//CRLF
    end if
    call energy_and_forces(total_energy,total_forces,mdl, &
         properties_only,ngwfs_converged, is_cmplx, kpt)
    pub_maxit_pen = orig_maxit_pen
    ! jd: Needed if more vacuum-solvent cycles will be performed
    !     (e.g. in force test or MD)
    pub_is_auto_solvation = .true.
    ! jd: Needed not to have the next vacuum calculation read restart files
    !     (in forcetest solvation only makes sense with READ_TB_NGWFS=F,
    !      and READ_DKN=F). Also pub_maxit_pen needs to be restored.
    ! ebl: the exception to the above being geometry optimisation, where
    !      reading of denskern and ngwfs is determined by the geometry_optimise
    !      module, rather than here
    if (pub_task /= 'GEOMETRYOPTIMIZATION') then
       pub_read_denskern = .false.
       pub_read_tightbox_ngwfs = .false.
       pub_read_hamiltonian = .false.
    end if
    restart_read_filename_prefix = ''
    restart_write_filename_prefix = ''
    if(pub_on_root) then
       write(stdout,'(a)') CRLF//'================&
            &=== Finished the calculation in solvent ========================'
       write(stdout,'(a)') '= Adjusting the following parameters: '
       write(stdout,'(a,l1)') '= read_denskern: ', pub_read_denskern
       write(stdout,'(a,l1)') '= read_tightbox_ngwfs: ', pub_read_tightbox_ngwfs
       write(stdout,'(a,l1)') '= read_hamiltonian: ', pub_read_hamiltonian
       write(stdout,'(a,i0,a)') '= maxit_pen: ', &
            pub_maxit_pen,' (restored to user-specified value)'
       write(stdout,'(3a)') '= restart read filename prefix: "',&
            trim(restart_read_filename_prefix),'"'
       write(stdout,'(3a)') '= restart write filename prefix: "',&
            trim(restart_write_filename_prefix),'"'
       write(stdout,'(a)') '================&
            &================================================================'

       ! jd: Figure out the non-solvent contribution to energy, i.e. the
       !     usual DFT terms except molecular Hartree.
       nonelectrostatic_dft_energy = total_energy - &
            implicit_solvent_energy_terms%E_elec_fixed - &
            implicit_solvent_boltzmann_energy_contrib(&
            implicit_solvent_energy_terms)
       if(pub_is_include_apolar) then
          nonelectrostatic_dft_energy = nonelectrostatic_dft_energy - &
               implicit_solvent_energy_terms%E_apolar_cavitation - &
               implicit_solvent_energy_terms%E_apolar_disrep
       end if
       nonelectrostatic_dft_energy_in_vacuum = total_energy_in_vacuum - &
            implicit_solvent_energy_terms%E_elec_fixed_vac

       write(stdout,'(a)') CRLF//'==========================================&
            &======================================'
       write(stdout,'(a,a)') 'IS: Calculation in solvent completed', &
            merge(' (NGWFs converged).      ', ' (*NO* NGWF convergence).', &
            ngwfs_converged)
       write(stdout,'(a)') CRLF//'Individual components of total energy in &
            &solvent:     hartree           kcal/mol'
       write(stdout,'(a,f24.14,a,f16.6,a)') &
            '- Usual non-electrostatic DFT terms: ', &
            nonelectrostatic_dft_energy, '   ', &
            nonelectrostatic_dft_energy * HARTREE_TO_KCAL_PER_MOL
       write(stdout,'(a,f22.14,a,f16.6,a)') &
            '- Electrostatic fixed charge energy:   ', &
            implicit_solvent_energy_terms%E_elec_fixed, '   ', &
            implicit_solvent_energy_terms%E_elec_fixed * &
            HARTREE_TO_KCAL_PER_MOL
       if(pub_is_pbe /= "NONE") then
          write(stdout,'(a,f22.14,a,f16.6,a)') &
               '- Electrostatic mobile charge energy:  ', &
               implicit_solvent_energy_terms%E_elec_mob, '   ', &
               implicit_solvent_energy_terms%E_elec_mob * &
               HARTREE_TO_KCAL_PER_MOL
          write(stdout,'(a,f22.14,a,f16.6,a)') &
               '- Accessibility (steric) correction:   ', &
               implicit_solvent_energy_terms%E_acc, '   ', &
               implicit_solvent_energy_terms%E_acc * &
               HARTREE_TO_KCAL_PER_MOL
          write(stdout,'(a,f22.14,a,f16.6,a)') &
               '- Osmotic pressure contribution:       ', &
               implicit_solvent_energy_terms%E_osmo, '   ', &
               implicit_solvent_energy_terms%E_osmo * &
               HARTREE_TO_KCAL_PER_MOL
          write(stdout,'(a,f22.14,a,f16.6,a)') &
               '- Ionic atmosph. rearrangement entropy:', &
               implicit_solvent_energy_terms%E_atmo, '   ', &
               implicit_solvent_energy_terms%E_atmo * &
               HARTREE_TO_KCAL_PER_MOL
          write(stdout,'(a,f22.14,a,f16.6,a)') &
               '- Chemical potential contribution:     ', &
               implicit_solvent_energy_terms%E_chempot, '   ', &
               implicit_solvent_energy_terms%E_chempot * &
               HARTREE_TO_KCAL_PER_MOL
       end if
       if(pub_is_include_apolar) then
          if(pub_is_apolar_method=='SAV') then
             write(stdout,'(a,f22.14,a,f16.6,a)') &
                  '- Apolar SASA cavitation energy:       ', &
                  implicit_solvent_energy_terms%E_apolar_sasa, '   ', &
                  implicit_solvent_energy_terms%E_apolar_sasa * &
                  HARTREE_TO_KCAL_PER_MOL
             write(stdout,'(a,f22.14,a,f16.6,a)') &
                  '- Apolar SAV cavitation energy:        ', &
                  implicit_solvent_energy_terms%E_apolar_sav, '   ', &
                  implicit_solvent_energy_terms%E_apolar_sav * &
                  HARTREE_TO_KCAL_PER_MOL
             write(stdout,'(a,f22.14,a,f16.6,a)') &
                  '- Apolar dispersion-repulsion energy:  ', &
                  implicit_solvent_energy_terms%E_apolar_disrep, '   ', &
                  implicit_solvent_energy_terms%E_apolar_disrep * &
                  HARTREE_TO_KCAL_PER_MOL
          else
             write(stdout,'(a,f22.14,a,f16.6,a)') &
                  '- Apolar cavitation energy:            ', &
                  implicit_solvent_energy_terms%E_apolar_cavitation, '   ', &
                  implicit_solvent_energy_terms%E_apolar_cavitation * &
                  HARTREE_TO_KCAL_PER_MOL
             write(stdout,'(a,f22.14,a,f16.6,a)') &
                  '- Apolar dispersion-repulsion energy:  ', &
                  implicit_solvent_energy_terms%E_apolar_disrep, '   ', &
                  implicit_solvent_energy_terms%E_apolar_disrep * &
                  HARTREE_TO_KCAL_PER_MOL
          end if
       end if
       write(stdout,'(a)') ' -----------------------------------------------&
            &--------------------------------'
       write(stdout,'(a,f22.14,a,f16.6,a)') &
            '- Total energy in solvent:             ', &
            total_energy, '   ', &
            total_energy * HARTREE_TO_KCAL_PER_MOL

       ! gab: print the SASA and SAV values.
       write(stdout,'(a)') CRLF//'Surface accessible surface area and volume&
            &:    SASA / bohr**2      SAV / bohr**3'
       write(stdout,'(a,f16.5,a,f16.5,a)') &
            '- Quantum surface area/volume:               ', &
            implicit_solvent_energy_terms%cavity_surface_area, '   ', &
            implicit_solvent_energy_terms%cavity_volume

       ! jd: Electrostatic energy (1/2 \int \rho_tot \phi dr).
       electrostatic_energy = implicit_solvent_energy_terms%E_elec_fixed
       if(pub_is_pbe /= "NONE") then
          electrostatic_energy = electrostatic_energy + &
               implicit_solvent_energy_terms%E_elec_mob
       end if
       electrostatic_energy_in_vacuum = &
            implicit_solvent_energy_terms%E_elec_fixed_vac

       ! jd: Apolar energy (cavitation, dispersion, repulsion)
       apolar_energy = implicit_solvent_energy_terms%E_apolar_cavitation + &
            implicit_solvent_energy_terms%E_apolar_disrep

       ! jd: Polar energy (electrostatic & change in non-es DFT terms)
       polar_energy = electrostatic_energy + nonelectrostatic_dft_energy &
            - electrostatic_energy_in_vacuum &
            - nonelectrostatic_dft_energy_in_vacuum

       write(stdout,'(a)') CRLF//'Components of total energy in &
            &solvent:                hartree           kcal/mol'
       write(stdout,'(a,f24.14,a,f16.6,a)') &
            '- Usual non-electrostatic DFT terms: ', &
            nonelectrostatic_dft_energy, '   ', &
            nonelectrostatic_dft_energy * HARTREE_TO_KCAL_PER_MOL
       write(stdout,'(a,f22.14,a,f16.6,a)') &
            '- Electrostatic energy:                ', &
            electrostatic_energy, '   ', &
            electrostatic_energy * HARTREE_TO_KCAL_PER_MOL
       if(pub_is_pbe /= "NONE") then
          write(stdout,'(a,f22.14,a,f16.6,a)') &
               '- Non-electrostatic electrolyte terms: ', &
            implicit_solvent_boltzmann_energy_contrib(&
            implicit_solvent_energy_terms), '   ', &
            implicit_solvent_boltzmann_energy_contrib(&
            implicit_solvent_energy_terms) * HARTREE_TO_KCAL_PER_MOL
       end if
       if(pub_is_include_apolar) then
          write(stdout,'(a,f22.14,a,f16.6,a)') &
               '- Apolar energy terms:                 ', &
               implicit_solvent_energy_terms%E_apolar_cavitation + &
               implicit_solvent_energy_terms%E_apolar_disrep, '   ', &
               (implicit_solvent_energy_terms%E_apolar_cavitation + &
               implicit_solvent_energy_terms%E_apolar_disrep) * &
               HARTREE_TO_KCAL_PER_MOL
       end if
       write(stdout,'(a)') ' -----------------------------------------------&
            &--------------------------------'
       write(stdout,'(a,f22.14,a,f16.6,a)') &
            '- Total energy in solvent:             ', &
            total_energy, '   ', &
            total_energy * HARTREE_TO_KCAL_PER_MOL

       write(stdout,'(a)') CRLF//'Calculation of free energy of solvation:&
            &              hartree           kcal/mol'
       write(stdout,'(a,f22.14,a,f16.6,a)') &
            '- Total energy in solvent:         (+) ', &
            total_energy, '   ', &
            total_energy * HARTREE_TO_KCAL_PER_MOL
       write(stdout,'(a,f22.14,a,f16.6,a)') &
            '- Total energy in vacuum:          (-) ', &
            total_energy_in_vacuum, '   ', &
            total_energy_in_vacuum * HARTREE_TO_KCAL_PER_MOL
       if(pub_is_pbe /= "NONE") then
          write(stdout,'(a,f22.14,a,f16.6,a)') &
               '- Energy of pure electrolyte:      (-) ', &
               implicit_solvent_energy_terms%E_pure, '   ', &
               implicit_solvent_energy_terms%E_pure * HARTREE_TO_KCAL_PER_MOL
       end if

       ! jd: Calculate final solvation energy as the difference between
       !     in-solvent and in-vacuum
       solvation_energy = total_energy - total_energy_in_vacuum

       ! jd: Mind the pure electrolyte term in Boltzmann solvation
       if(pub_is_pbe /= "NONE") then
          solvation_energy = solvation_energy - &
               implicit_solvent_energy_terms%E_pure
       end if

       write(stdout,'(a)') ' -----------------------------------------------&
            &--------------------------------'
       write(stdout,'(a,f22.14,a,f16.6,a)') &
            '- Total free energy of solvation:      ', solvation_energy, &
            '   ', solvation_energy * HARTREE_TO_KCAL_PER_MOL

       if(any(pub_multigrid_bc_is_periodic)) then
          write(stdout,'(/a)') 'WARNING: Non-OBC calculation. Due to electro&
               &static interactions between periodic'//CRLF//'         image&
               &s being mostly screened in solvent and unscreened in vacuum,&
               &the'//CRLF//'         calculated solvation energy is meaning&
               &less unless you subsequently'//CRLF//'         correct the i&
               &n-vacuum energy to remove periodic interactions.'
       end if

       write(stdout,'(a)') CRLF//'Components of polar term in f.e. of &
            &solvation:        hartree           kcal/mol'
       write(stdout,'(a,f22.14,a,f16.6,a)') &
            '- Electrostatic:                       ', &
            electrostatic_energy - electrostatic_energy_in_vacuum, &
            '   ', (electrostatic_energy - electrostatic_energy_in_vacuum) *&
            HARTREE_TO_KCAL_PER_MOL
       write(stdout,'(a,f22.14,a,f16.6,a)') &
            '- Change in nonelectrostatic DFT terms:', &
            (nonelectrostatic_dft_energy - &
            nonelectrostatic_dft_energy_in_vacuum), '   ', &
            (nonelectrostatic_dft_energy - &
            nonelectrostatic_dft_energy_in_vacuum) * HARTREE_TO_KCAL_PER_MOL
       write(stdout,'(a)') ' -----------------------------------------------&
            &--------------------------------'
       write(stdout,'(a,f22.14,a,f16.6,a)') &
            '- Polar term in f.e. of solvation:     ', &
            polar_energy, '   ', polar_energy * HARTREE_TO_KCAL_PER_MOL

       ! jd: Check for the free energy of solvation
       solvation_energy_check = polar_energy + apolar_energy

       if(pub_is_pbe /= "NONE") then
          ! jd: Mind the pure electrolyte term in Boltzmann solvation
          solvation_energy_check = solvation_energy_check - &
               implicit_solvent_energy_terms%E_pure
          ! jd: ... and the non-electrostatic terms from the electrolyte
          solvation_energy_check = solvation_energy_check + &
               implicit_solvent_boltzmann_energy_contrib(&  ! <- all
               implicit_solvent_energy_terms) &             ! <- all
               - implicit_solvent_energy_terms%E_elec_mob   ! minus elec
       end if

       write(stdout,'(a)') CRLF//'Components of free energy of solvation: &
            &              hartree           kcal/mol'
       write(stdout,'(a,f22.14,a,f16.6,a)') &
            '- Polar term in f.e. of solvation: (+) ', &
            polar_energy, '   ', polar_energy * HARTREE_TO_KCAL_PER_MOL
       write(stdout,'(a,f22.14,a,f16.6,a)') &
            '- Apolar (cavitation, dis., rep.): (+) ', &
            apolar_energy, '   ', apolar_energy * HARTREE_TO_KCAL_PER_MOL
       if(pub_is_pbe /= "NONE") then
          write(stdout,'(a,f22.14,a,f16.6,a)') &
               '- Non-es. electrolyte terms:       (+) ', &
               implicit_solvent_boltzmann_energy_contrib(&
               implicit_solvent_energy_terms), '   ', &
               implicit_solvent_boltzmann_energy_contrib(&
               implicit_solvent_energy_terms) * HARTREE_TO_KCAL_PER_MOL
          write(stdout,'(a,f22.14,a,f16.6,a)') &
               '- Energy of pure electrolyte:      (-) ', &
               implicit_solvent_energy_terms%E_pure, '   ', &
               implicit_solvent_energy_terms%E_pure * HARTREE_TO_KCAL_PER_MOL
       end if
       write(stdout,'(a)') ' -----------------------------------------------&
            &--------------------------------'
       write(stdout,'(a,f22.14,a,f16.6,a)') &
            '- Total free energy of solvation:      ', &
            solvation_energy_check, '   ', solvation_energy_check * &
            HARTREE_TO_KCAL_PER_MOL

       if(any(pub_multigrid_bc_is_periodic)) then
          write(stdout,'(/a)') 'WARNING: Non-OBC calculation. Due to electro&
               &static interactions between periodic'//CRLF//'         image&
               &s being mostly screened in solvent and unscreened in vacuum,&
               &the'//CRLF//'         calculated solvation energy is meaning&
               &less unless you subsequently'//CRLF//'         correct the i&
               &n-vacuum energy to remove periodic interactions.'
       end if

       write(stdout,'(a)') '===============================================&
            &================================='//CRLF
    end if

    if (present(return_converged)) return_converged = ngwfs_converged

  end subroutine energy_and_forces_auto_solvation

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine energy_and_forces(total_energy,total_forces,mdl, &
       properties_only,return_converged, is_cmplx, kpt, custom_restart_name)

    use augmentation, only: augmentation_box_init, augmentation_box_exit, &
         augmentation_density_on_grid
    use cdft_intermediate_cg, only: cdft_intermediate_restart_read, &
         cdft_intermediate_update_matrices,dft_nu_intermediate_restart_read,&
         dft_nu_intermediate_update_matrices
    use cell_grid, only: cell_grid_distribute, cell_grid_exit
    use classical_pot, only: classical_pot_ii_energy
    use comms, only: comms_barrier, comms_reduce, pub_on_root, &
         pub_total_num_procs
    use constants, only: stdout, DP, NORMAL, VERBOSE, UP, DN, CRLF, &
         EDA_ISOLATED, EDA_FROZENNONIDEM, EDA_FULL, paw_en_size, &
         REP_SWEX_PROPERTIES_DMA_1, REP_SWEX_PROPERTIES_DMA_2, &
         HARTREE_TO_KCAL_PER_MOL, LONG
    use conduction, only: conduction_ngwf_optimise
    use couplings, only: couplings_calculate
    use cutoff_coulomb, only: cutoff_coulomb_ii_energy, &
         cutoff_coulomb_init, cutoff_coulomb_exit, &
         cutoff_coulomb_check_boundaries
    use datatypes, only: data_functions_alloc, data_functions_dealloc
    use density, only: density_on_grid
    use dense, only: dense_create, dense_convert
    use dftb, only: dftb_init_stage_2, dftb_free, dftb_rep, dftb_electronic_force, &
         dftb_srb, dftb_ies, dftb_print_energy_components, dftb_print_force, &
         dftb_qc_print
    use dma, only: dma_expand_ngwf_pairs, dma_free_swexes
    use eda_driver_supermol, only: eda_driver_supermol_run, &
         eda_driver_prep_supermol
    use eigenstates, only: eigenstates_init_pdos, eigenstates_cleanup_pdos, &
         eigenstates_lumo_search
    use electronic_history, only: elec_history_restore_history_parm
    use electronic_init, only: electronic_init_denskern
    use ensemble_dft, only: edft_create, edft_destroy
    use ensemble_dft_type, only: EDFT_MODEL
    use ewald, only: ewald_calculate_energy
    use fourier, only: fourier_init_fftbox, fourier_init_box, &
         fourier_exit_fftbox, fourier_exit_box
    use forces, only: forces_calculate, forces_calculate_correction
    use function_basis, only: FUNC_BASIS, function_basis_allocate, &
         function_basis_deallocate, function_basis_distribute, &
         function_basis_init_tight_boxes, function_basis_init_spheres, &
         function_basis_gath_all_tbs, function_basis_init_uni_tb, &
         function_basis_exit_uni_tb, function_basis_estimate_size, &
         function_basis_copy_spheres
    use fragment_data, only: pub_frag_data
    use fragment_matrix_ops, only: fmo_apply_intrafrag_partitions
    use fragment_scfmi, only: scfmi_allocate_matrices, scfmi_deallocate_matrices
    use hamiltonian, only : hamiltonian_dens_indep_matrices
    use hf_exchange, only: HFX_STATE, hf_exchange_init, hf_exchange_free
    use hubbard_build, only: HUBBARD_MODEL, hubbard_model_init, &
         hubbard_model_exit, hubbard_species_proj, &
         hubbard_build_consist, &
         hubbard_projector_update, hubbard_test_convergence, &
         hubbard_species_exit_proj, &
         hubbard_build_matrices, hubbard_build_matrices_exit
    use function_ops, only: function_ops_brappd_ketppd
    use is_poisson, only: POISSON_PROBLEM, deallocate_poisson_problem
    use is_smeared_ions, only: smeared_ion_exit, smeared_ion_initialise
    use is_solvation, only: initialize_solvation_problem, have_rho_ion, &
         have_initial_eps, implicit_solvent_energy_terms
    use is_solvation_boltzmann, only: implicit_solvent_boltzmann_init, &
         implicit_solvent_boltzmann_exit
    use kernel, only: DKERN, kernel_create, kernel_destroy
    use lnv, only: lnv_reset_trial_step_to_initial
    use lr_phonons, only: lr_phonons_calculate
    use lr_tddft, only: lr_tddft_calc_excitations
    use model_type, only: MODEL
    use ngwf_cg, only: ngwf_cg_optimise, ngwf_cg_emft_reset
    use ngwf_representation, only: NGWF_REP, ngwf_rep_create, ngwf_rep_destroy, &
         NGWF_HAM, ngwf_ham_create, ngwf_ham_destroy, ngwf_rep_register_change, &
         ngwf_ham_register_change
    use ngwfs, only: ngwfs_initialise, ngwfs_initialise
    use parallel_strategy, only: parallel_strategy_distr_atoms, &
         parallel_strategy_check_atoms, parallel_strategy_exit, &
         parallel_strategy_distr_kpoints, PARAL_INFO
    use paw, only: paw_read_species, paw_species_init_proj, &
         paw_species_init_core_wvfns, paw_species_exit_core_wvfns
    use pbc_corrections, only: pbc_corr_check_cell_size, pbc_corr_exit
    use penalty, only: pen_reset_trial_step_to_initial
    use potential, only: potential_add_efield_ion_energy
    use projectors, only: PROJECTOR_SET, projectors_create_real, &
         projectors_deallocate_set, projectors_func_ovlp_box
    use properties, only: properties_calculate
    use polarisable_embedding, only: polarisable_embedding_init, &
         polarisable_embedding_exit
    use pseudopotentials, only: pseudo_species_init_proj
    use remote, only: remote_init, remote_exit
    use restart, only: restart_kernel_read, restart_ngwfs_tightbox_input, &
         restart_read_filename_prefix, restart_write_filename_prefix
    use rundat, only: pub_maxit_pen, pub_coulomb_cutoff, &
         pub_coulomb_cutoff_type, pub_dispersion, pub_do_properties, &
         pub_do_tddft, pub_md_properties, pub_hubbard, &
         pub_hubbard_restart, pub_hubbard_atomsolve, pub_hub_on_the_fly, &
         pub_hub_max_iter, pub_cond_calculate, pub_lr_tddft_calculate,&
         pub_mt_cutoff, pub_lr_phonons_calculate,&
         pub_output_detail, pub_any_nl_proj, &
         pub_paw, pub_write_forces, pub_task, pub_nlcc, pub_spin, &
         pub_maxit_ngwf_cg, pub_aug_den_dim, &
         pub_nonsc_forces, pub_aug, pub_realspace_projectors, &
         pub_constant_efield, pub_cdft, pub_edft_init_maxit, &
         pub_cdft_continuation, pub_maxit_cdft_u_cg, pub_spin_fac, &
         pub_use_aux_ngwfs, pub_use_hfx, pub_eels_calculate, pub_edft, &
         pub_is_smeared_ion_rep, pub_is_separate_restart_files, &
         pub_is_auto_solvation, pub_is_implicit_solvent, pub_num_spins, &
         pub_is_include_apolar, pub_write_denskern, pub_write_tightbox_ngwfs, &
         pub_read_denskern, pub_read_tightbox_ngwfs, pub_use_swx, pub_is_pbe, &
         pub_dos_smear, pub_pdos_max_l, pub_dft_nu, pub_dft_nu_continuation, &
         md_global_restart, pub_num_kpoints, PUB_1K, &
         pub_eda, pub_frag_counter, pub_tightbox_fft_coarse, &
         pub_eda_mode, pub_eda_scfmi, pub_eda_preptool, pub_eda_split_atoms, &
         pub_dmft_points, pub_kpoint_method, &
         pub_turn_off_ewald, pub_print_qc, &
         pub_xc_ke_density_required, pub_pol_emb_pot, pub_hhf_nstates, &
         pub_hhf_factor, pub_dma_calculate, pub_rootname, pub_debug, &
         pub_ion_ion_bc_is_periodic, pub_smeared_ion_bc_is_periodic, &
         pub_dmft_read, pub_do_fandt, &
         pub_quit_region, pub_freeze_switch_steps, pub_parallel_scheme, &
         pub_use_activehfx, pub_use_activeswx, pub_active_region, &
         pub_emft_follow, pub_is_pbe_bc_debye_screening, &
         pub_multigrid_bc_is_periodic, pub_ngwf_gradient_needed, &
         pub_is_pbe_neutralisation_scheme, pub_write_hamiltonian, &
         pub_read_hamiltonian, pub_is_apolar_method, &
         pub_is_solvent_surf_tension, pub_is_apolar_scaling_factor, &
         pub_is_solvent_pressure, pub_dftb, pub_dftb_bc_is_periodic, &
         pub_is_restart_vac_from_vac, pub_is_emft_cavity, pub_emft, &
         pub_block_orthogonalise, pub_build_bo

    use services, only: services_flush
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_destroy, &
         sparse_exit, sparse_mod_init, sparse_init_blocking_scheme, &
         BLKS_NGWF, BLKS_PROJ, BLKS_HUB_PROJ, &
         BLKS_COND, BLKS_JOINT, BLKS_CORE, BLKS_AUX
    use sparse_embed, only: sparse_embed_extract_from_array, &
         sparse_embed_destroy_extracted_array
    use sph_harm_rotation, only: sph_harm_rot_initialise, sph_harm_rot_exit
    use sw_expansion, only: swx_init_module, swx_cleanup_module
    use sw_resolution_of_identity, only: swri_init_module_stage_1, &
         swri_init_module_stage_2, swri_cleanup_module
    use tddft, only: tddft_calculate
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, &
         utils_assert, utils_banner, utils_qc_print
    use vdwcorrection, only: vdwcorrection_calculate_energy,&
         pub_dispersion_energy, vdwcorrection_calculate_forces
    use xc, only: xc_hfxinit, xc_init, xc_exit

    implicit none

    ! Arguments
    type(MODEL), target, intent(inout)  :: mdl
    real(kind=DP), intent(out)  :: total_energy
    real(kind=DP), intent(out)  :: total_forces(1:3,1:mdl%nat)
    logical,intent(in),optional :: properties_only
    logical,intent(out),optional :: return_converged
    ! agrecocmplx: optional argument for complex NGWFs
    logical, intent(in), optional :: is_cmplx
    ! agrecokpt: optional argument for single kpt,
    ! in fractional coordinates
    real(kind=DP), intent(in), optional :: kpt(3)
    character(len=*), optional, intent(in) :: custom_restart_name

    ! Local Variables
    ! rc2013: allow structures to include components from multiple systems
    type(FUNC_BASIS) :: ngwf_basis(mdl%nsub)
    type(FUNC_BASIS) :: aux_ngwf_basis(mdl%nsub)
    type(FUNC_BASIS) :: cond_ngwf_basis(mdl%nsub)
    type(FUNC_BASIS) :: joint_ngwf_basis(mdl%nsub)
    type(FUNC_BASIS) :: proj_basis(mdl%nsub)
    type(FUNC_BASIS) :: hub_proj_basis(mdl%nsub)
    type(FUNC_BASIS) :: core_basis(mdl%nsub)
    type(NGWF_REP) :: rep
    type(NGWF_HAM) :: ham
    type(HUBBARD_MODEL) :: hub
    type(PROJECTOR_SET) :: nl_projectors(mdl%nsub)
    type(PROJECTOR_SET) :: core_wvfns(mdl%nsub)
    type(DKERN) :: denskern
    type(HFX_STATE), target :: hfxstate
    real(kind=DP), allocatable, dimension(:,:) :: ngwf_nonsc_forces
    real(kind=DP), dimension(:,:,:,:), allocatable :: temp_density_fine
    real(kind=DP), dimension(:,:,:,:,:), allocatable :: nhat_den_grad
    ! rc2013: lhxc_fine can vary across regions -- include extra index
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: lhxc_fine
    ! JCW: Gradient of XC energy per unit volume wrt kinetic energy
    ! JCW: density. Required to evaluate lhxc matrix elements for meta-GGAs.
    real(kind=DP), allocatable, dimension(:,:,:,:) :: dfdtau_fine
    real(kind=DP) :: throwaway
    real(kind=DP) :: proj_sphere_padding(mdl%nsub)
    integer :: ierr
    integer :: local_size
    integer(kind=LONG) :: global_size
    integer :: nelecs
    integer :: is
    logical :: local_properties
    logical :: ngwfs_converged
    logical :: out_of_runtime
    integer :: hub_proj_iteration ! ddor: The DFT+U projector-consistency
                                  ! ddor: iteration we are on
    logical :: check_dmft_store

    type(EDFT_MODEL) :: edft  ! ars: ensemble DFT container
    type(POISSON_PROBLEM) :: dummy_solv_problem ! jd: for solvent init

    ! jhl52: aux
    type(NGWF_REP) :: aux_rep
    type(NGWF_HAM) :: aux_ham
    type(DKERN) :: aux_denskern

    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecokpt: single k-point, in Cartesian coordinates
    real(kind=DP) :: kpt_cart(3)

    ! rc2013: par pointer
    type(PARAL_INFO), pointer :: par
    integer :: isub, jsub
    logical :: first_sparse_call  ! rc2013: sparse mod boolean
    character(len=1) :: isub_str

    ! jcap: temporary variables for embedding loops
    real(kind=DP),allocatable :: sub_temp_density_fine(:,:,:,:)
    integer :: tot_num
    type(SPAM3) :: kern_array(pub_num_spins)
    integer :: ii
    character(200) :: vac_ngwfs_ext
    logical :: temp_emft, temp_emft_follow, temp_build_bo
    character(len=*), parameter :: myself = 'energy_and_forces'
    ! jd: Storage for when we temporarily overwrite pub_rootname
    character(len=80) :: old_pub_rootname

    ! jd: DFTB parameters and energy terms
    real(kind=DP) :: dftb_rep_energy
    real(kind=DP) :: dftb_srb_energy
    real(kind=DP) :: dftb_ies_energy
    real(kind=DP) :: dftb_eht_energy
    real(kind=DP) :: dftb_disp_energy
    real(kind=DP) :: F_dftb_rep(1:3,1:mdl%nat)
    real(kind=DP) :: F_dftb_srb(1:3,1:mdl%nat)
    real(kind=DP) :: F_dftb_ies(1:3,1:mdl%nat)
    real(kind=DP) :: F_dftb_dis(1:3,1:mdl%nat)
    real(kind=DP) :: F_dftb_eht(1:3,1:mdl%nat)

    ! -------------------------------------------------------------------------

    ! Initialisations

    call timer_clock(myself,1)

    ! agrecocmplx
    if (present(is_cmplx)) then
        loc_cmplx = is_cmplx
    else
        loc_cmplx = .false.
    end if

    if (loc_cmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Complex-valued NGWFs will be used.'
    end if

    ! jme: initialise k-points
    call internal_kpoints_init(mdl, kpt_cart, kpt)

    ! set local variable properties_only: this overrides global setting of
    ! whether or not to do a properties calculation on this run
    local_properties = .false.
    if (present(properties_only)) local_properties = properties_only
    ngwfs_converged = .true. ! set to true if no optimisation requested
    out_of_runtime = .false.

    ! jd: Fresh calculation, reset lnv and penalty trial steps to initial
    lnv_reset_trial_step_to_initial = .true.
    pen_reset_trial_step_to_initial = .true.

    ! rc2013: Boolean to check if this is the first call of sparse_mod_init
    first_sparse_call = .true.

    ! rc2013: decide how many procs each subsystem should have
    call internal_initialise_par

    ! Adopt parallel strategy for atom distribution
    if (pub_on_root .and. (pub_output_detail >= NORMAL) .and. &
         .not.(pub_eda) .and. (mdl%nsub .eq. 1)) &
         write(stdout,'(/a)',advance='no') 'Determining parallel strategy ...'

    ! rc2013: loop over regions for atom distibution and allocations
    do isub=1,mdl%nsub
       par => mdl%regions(isub)%par
       write(isub_str,'(i1)') isub
       if (pub_on_root .and. (pub_output_detail >= NORMAL) .and. &
            .not.(pub_eda) .and. (mdl%nsub .gt. 1)) &
            write(stdout,'(a,i1,a3)',advance='no') &
            'Determining parallel strategy for subsystem ', isub, '...'
       call parallel_strategy_distr_atoms(par,mdl%regions(isub)%elements, &
            mdl%cell)

       if (pub_on_root .and. (pub_output_detail >= NORMAL) .and. &
            .not.(pub_eda) .and. (mdl%nsub .gt. 1) .and. &
            (pub_parallel_scheme .ne. 'NONE')) &
            write(stdout,'(a4,i3,a)') '... ', par%num_procs, &
            ' processors allocated to this region'

       ! CHECK CENTRES
       ! cks: make sure that no atom is outside the simulation cell
       call parallel_strategy_check_atoms(par,mdl%regions(isub)%elements, &
            mdl%cell)

       ! jd: make sure that the cell is large enough if MT correction is in effect
       if (pub_on_root .and. pub_mt_cutoff /= 0.0_DP) &
            call pbc_corr_check_cell_size(mdl%regions(isub)%elements,mdl%cell)

       ! Distribute the NGWFs
       if (.not.(pub_eda)) then
          call function_basis_allocate(ngwf_basis(isub),par%num_ngwfs, &
               'ngwfs', par)
          call function_basis_distribute(ngwf_basis(isub), &
               mdl%regions(isub)%elements, par)
       else
          ! mjsp: Copy the precomputed fragment/supermolecule NGWFs
          ! mjsp: object to this module's ngwf_basis object
          ngwf_basis(isub) = pub_frag_data(pub_frag_counter)%ngwf_basis(isub)
       end if


       ! Exit if there is a proc with no NGWFs
       ! rc2013: unless we're using a regional parallel strategy
       if(pub_parallel_scheme == 'NONE') then
          call utils_assert (minval(ngwf_basis(isub)%num_on_proc) >= 1, &
               'Too many processors (MPI ranks) for this system size. Sorry!'//CRLF//&
               'Your options include: '//CRLF//&
               '1) Reducing the number of processors (MPI ranks), '//CRLF//&
               '2) Increasing the number of atoms, '//CRLF//&
               '3) Using OpenMP parallelisation to your advantage, by reducing'//CRLF//&
               '   the number of MPI ranks, and simultaneously using threads_max'//CRLF//&
               '   and other threads_* keywords in your input file to have each'//CRLF//&
               '   MPI rank spawn multiple threads.'//CRLF//CRLF//&
               'Occasionally, turning off the space-filling curve by specifying'//CRLF//&
               '"use_space_filling_curve F" in the input file might help --'//CRLF//&
               'but only if the resultant atom ordering better lends itself to '//CRLF//&
               'distributing across processors. If you do this, be aware that'//CRLF//&
               'if the ions move (geometry optimisation, TS search, MD), the'//CRLF//&
               'ordering may once again become unfavourable and this error might'//CRLF//&
               'come back in the course of the calculation.')
       end if

       ! Distribute the nonlocal projectors if required
       if (pub_any_nl_proj) then
          call function_basis_allocate(proj_basis(isub),par%num_projectors, &
               'projs',par)
          call function_basis_distribute(proj_basis(isub), &
               mdl%regions(isub)%elements, par)
       end if

       ! Distribute PAW partial waves if required
       if (pub_paw) then
          call function_basis_allocate(proj_basis(isub),par%num_pawpws, &
               'pawpws', par)
          call function_basis_distribute(proj_basis(isub), &
               mdl%regions(isub)%elements, par)

          if (pub_eels_calculate) then
             call function_basis_allocate(core_basis(isub),par%num_corewfs, &
                  'corewfs', par)
             call function_basis_distribute(core_basis(isub), &
                  mdl%regions(isub)%elements, par)
          end if

       end if

       ! Distribute Hubbard projectors if required
       if (pub_hubbard) then
          call function_basis_allocate(hub_proj_basis(isub),par%num_hub_proj, &
               'hub_projs', par)
          call function_basis_distribute(hub_proj_basis(isub), &
               mdl%regions(isub)%elements, par)
       end if

       ! Distribute Conduction NGWFs if required (for COND, TDDFT or LR_PHONONS)
       if (pub_cond_calculate .or. pub_lr_tddft_calculate) then
          call function_basis_allocate(cond_ngwf_basis(isub), &
               par%num_ngwfs_cond, 'ngwfs_cond', par)
          call function_basis_distribute(cond_ngwf_basis(isub), &
               mdl%regions(isub)%elements, par)

          call function_basis_allocate(joint_ngwf_basis(isub), par%num_ngwfs &
               + par%num_ngwfs_cond,'ngwfs_joint', par)
          call function_basis_distribute(joint_ngwf_basis(isub), &
               mdl%regions(isub)%elements, par)
       end if

       ! Distribute Auxiliary NGWFs if required (not currently used)
       if (pub_use_aux_ngwfs) then
          call function_basis_allocate(aux_ngwf_basis(isub),par%num_ngwfs_aux, &
               'ngwfs_aux', par)
          call function_basis_distribute(aux_ngwf_basis(isub), &
               mdl%regions(isub)%elements, par)
       end if

       ! agrecokpt: distribute k-points if required
       if (mdl%nkpoints>1) then
          call utils_assert(mdl%nsub .eq. 1, myself//': subsystem calculations &
               &are not compatible with multiple k-points. Consider placing &
               &all atoms within the same region in %block species_ngwf_regions.')
          par%num_kpoints = mdl%nkpoints
          call parallel_strategy_distr_kpoints(par, mdl%unique_kpoints)
       end if

       ! Report results of distribution
       if (pub_on_root .and. pub_output_detail >= VERBOSE) then
          write(stdout,'(3(a,i6))') '... NGWF load balancing: max ', &
               maxval(ngwf_basis(isub)%num_on_proc),';  min ', &
               minval(ngwf_basis(isub)%num_on_proc),';  average ', &
               nint(ngwf_basis(isub)%num/real(pub_total_num_procs,kind=DP))
          if(pub_debug .and. mdl%nsub .gt. 1) then
             do ii=0,pub_total_num_procs-1
                write(stdout,'(2(a,i6))') 'DEBUG: NGWFs on proc:', ii, &
                     ' = ', ngwf_basis(isub)%num_on_proc(ii)
             end do
          end if
          if (pub_any_nl_proj.or.pub_paw) then
             write(stdout,'(3(a,i6))') '... Projector load balancing: max ', &
                  maxval(proj_basis(isub)%num_on_proc),';  min ', &
                  minval(proj_basis(isub)%num_on_proc),';  average ', &
                  nint(proj_basis(isub)%num/real(pub_total_num_procs,kind=DP))
          end if
          if (pub_hubbard) then
             write(stdout,'(3(a,i6))') '... Hubbard Projector load balancing: max ', &
                  maxval(hub_proj_basis(isub)%num_on_proc),';  min ', &
                  minval(hub_proj_basis(isub)%num_on_proc),';  average ', &
                  nint(hub_proj_basis(isub)%num/real(pub_total_num_procs,kind=DP))
          end if
          if (pub_cond_calculate) then
             write(stdout,'(3(a,i6))') '... Conduction NGWF load balancing: max ', &
                  maxval(cond_ngwf_basis(isub)%num_on_proc),';  min ', &
                  minval(cond_ngwf_basis(isub)%num_on_proc),';  average ', &
                  nint(cond_ngwf_basis(isub)%num/real(pub_total_num_procs,kind=DP))
             write(stdout,'(3(a,i6))') '... Joint NGWF load balancing: max ', &
                  maxval(joint_ngwf_basis(isub)%num_on_proc),';  min ', &
                  minval(joint_ngwf_basis(isub)%num_on_proc),';  average ', &
                  nint(joint_ngwf_basis(isub)%num/real(pub_total_num_procs,kind=DP))
          end if
          if (pub_use_aux_ngwfs) then
             write(stdout,'(3(a,i6))') '... Auxiliary NGWF load balancing: max ', &
                  maxval(aux_ngwf_basis(isub)%num_on_proc),';  min ', &
                  minval(aux_ngwf_basis(isub)%num_on_proc),';  average ', &
                  nint(aux_ngwf_basis(isub)%num/real(pub_total_num_procs,kind=DP))
          end if
          ! agrecokpt: print information about k-points distribution
          if (mdl%nkpoints>1) then
             write(stdout,'(3(a,i6))') '... k-points load balancing: max ', &
                  maxval(par%num_kpoints_on_proc),';  min ', &
                  minval(par%num_kpoints_on_proc),';  average ', &
                  nint(par%num_kpoints/real(pub_total_num_procs,kind=DP))
          end if
       end if

       call comms_barrier
       if (pub_on_root .and. (pub_output_detail >= NORMAL)) &
            write(stdout,'(a)') '... done'
       call services_flush
    end do

    ! ndmh: for cutoff coulomb, check no NGWFs extend beyond cell boundaries
    if (pub_coulomb_cutoff) call cutoff_coulomb_check_boundaries(mdl)

    ! ebl: only recalculate if DMFT store files do not exist and pub_dmft_read
    !      is true
    check_dmft_store = .false.
    if(pub_dmft_points>0 .and. pub_dmft_read) then
        inquire(file=trim(pub_rootname)//'.ham1',exist=check_dmft_store)
    end if
    if (.not. check_dmft_store) then

       ! Calculate the Ion-Ion energy
       if ( (.not.pub_cond_calculate) .and. (.not.pub_lr_tddft_calculate) .and. &
            (.not.pub_lr_phonons_calculate) .and. (.not.pub_turn_off_ewald) &
            .and. (.not.pub_dftb) ) then
          ! jd: Direct summation requested instead of Ewald
          ! JCW: If not fully periodic...
          if (.not.all(pub_ion_ion_bc_is_periodic(1:3)) ) then
             if (pub_on_root) then
                write(stdout,'(a)',advance='no') &
                     'Calculating direct Ion-Ion energy ...'
             end if

             if (pub_coulomb_cutoff) then
                ! JCW: For cutoff Coulomb, we allow pub_coulomb_cutoff_type to
                ! JCW: determine the ion-ion boundary conditions

                ! cks: Different approach depending on whether there are classical atoms
                if (mdl%nat_classical > 0) then

                   select case (pub_coulomb_cutoff_type)
                      ! 0D periodicity
                   case('SPHERE','sphere','CYLINDER','cylinder')
                      call classical_pot_ii_energy(mdl%ewald_energy, & ! output
                           mdl%classical_elements,mdl%elements, &
                           mdl%nat,mdl%nat_classical)
                   case default
                      call utils_abort('Illegal cutoff type for cutoff Coulomb &
                           &with embedding')
                   end select

                else
                   ! jd: Usual cutoff Coulomb calculation
                   call cutoff_coulomb_ii_energy(mdl%regions,mdl%cell, &
                        mdl%ewald_energy,pub_coulomb_cutoff_type)
                endif

             else
                ! jd: Reuse the CC code for direct summation for other calculations
                !     that use direct summation, like MT or smeared-ions
                if (.not.any(pub_ion_ion_bc_is_periodic(1:3)) ) then
                   ! Fully open, use sphere
                   call cutoff_coulomb_ii_energy(mdl%regions,mdl%cell, &
                        mdl%ewald_energy, 'SPHERE')
                else
                   ! Mixed BC, not yet supported
                   call utils_abort("Mixed BCs for non-cutoff-Coulomb calculations &
                        &are not yet supported.")
                end if
             end if
          ! jd: Ewald by default
          ! JCW: If fully periodic...
          else
             if (pub_on_root) then
                write(stdout,'(a)',advance='no') &
                     'Calculating Ewald energy ...'
             end if

             if (mdl%nat_classical > 0) then
                ! cks: include classical atom charges in Ewald energy
                call ewald_calculate_energy(mdl%elements,mdl%nat,mdl%cell,&
                     mdl%ewald_energy,.false.,mdl%nat_classical, &
                     mdl%classical_elements)
                call ewald_calculate_energy(mdl%elements,mdl%nat,mdl%cell, &
                     throwaway,.true.,mdl%nat_classical)
             else
                call ewald_calculate_energy(mdl%elements,mdl%nat,mdl%cell, &
                     mdl%ewald_energy,.false.,mdl%nat_classical)
                call ewald_calculate_energy(mdl%elements,mdl%nat,mdl%cell, &
                     throwaway,.true.,mdl%nat_classical)
             endif

          end if
       else
          ! No need to calculate Ewald Energy during a COND or TDDFT Calculation
          ! jd: or if turn_off_ewald is specified.
          ! ab: or if doing a dftb calculation.
          mdl%ewald_energy = 0.0_DP
       end if

       ! ndmh: calculate energy of ions in external field
       if (any(abs(pub_constant_efield) > 0.0_DP)) &
            call potential_add_efield_ion_energy(mdl%ewald_energy, mdl%cell,par)

       call services_flush

       ! qoh: Calculate dispersion energy
       if (.not.(pub_dispersion == '0' .or. pub_dispersion == 'NONE')) then
          call vdwcorrection_calculate_energy(mdl%elements,mdl%cell,mdl%par)
       else
          pub_dispersion_energy = 0.0_DP
       end if
    end if

    ! INITIALISE PPDS, SPHERES, BOXES, KB_DENOMINATORS
    if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
         write(stdout,'(a)',advance='no') 'Basis initialisation ...'

    ! rc2013: do this for all subsystems
    do isub=1,mdl%nsub
       par => mdl%regions(isub)%par
       ! jcap: state which region we are initialising first
       if (pub_on_root .and. mdl%nsub.gt.1) then
          write(isub_str,'(i0)')isub
          write(stdout,'(/a)') utils_banner('<', 'Region '//isub_str)
       end if

       ! NGWF SPHERES
       call function_basis_init_spheres(ngwf_basis(isub), mdl%cell, mdl%fftbox, &
            par%elements_on_proc, par=par)

       ! FRAGMENT NGWF SPHERES
       ! mjsp: Also allocate the fragment NGWF storage if neccesary.
       ! mjsp: This is used to construct the frozen density state later.
       if ((pub_eda) .and. (pub_eda_mode == EDA_ISOLATED)) call &
          function_basis_init_spheres(pub_frag_data(pub_frag_counter)%ngwf_basis(isub), &
          pub_frag_data(pub_frag_counter)%mdl%cell, &
          pub_frag_data(pub_frag_counter)%mdl%fftbox, par%elements_on_proc, par=par)

       if (pub_on_root .and. pub_output_detail >= VERBOSE) &
            write(stdout,'(a)') '... NGWF spheres initialised'

       ! PROJECTOR SPHERES
       if (pub_any_nl_proj.or.pub_paw) then
          if (.not.pub_realspace_projectors) then
             call function_basis_init_spheres(proj_basis(isub), mdl%cell, &
                  mdl%fftbox, par%elements_on_proc, par=par)
          else
             proj_sphere_padding(isub) = &
                  2.0_DP*maxval(ngwf_basis(isub)%spheres(:)%radius)
             ! rc2013: need to be careful with padding for all NGWFs
             ! rc2013: currently this assigns different padding to each region
             call comms_reduce('MAX',proj_sphere_padding(isub))
             call function_basis_init_spheres(proj_basis(isub), mdl%cell, &
                  mdl%fftbox, par%elements_on_proc, proj_sphere_padding(isub), &
                  par=par)
          end if
          if (pub_on_root .and. pub_output_detail >= VERBOSE) &
               write(stdout,'(a)') '... Projector spheres initialised'

          if (pub_paw.and.pub_eels_calculate) then
             if (.not.pub_realspace_projectors) then
                call function_basis_init_spheres(core_basis(isub), mdl%cell, &
                     mdl%fftbox, par%elements_on_proc, par=par)
             else
                ! rc2013: currently this assigns different padding to each region
                proj_sphere_padding(isub) = &
                     2.0_DP*maxval(ngwf_basis(isub)%spheres(:)%radius)
                call comms_reduce('MAX',proj_sphere_padding(isub))
                call function_basis_init_spheres(core_basis(isub), mdl%cell, &
                     mdl%fftbox, par%elements_on_proc,proj_sphere_padding(isub), par=par)
             end if
             if (pub_on_root .and. pub_output_detail >= VERBOSE) &
                  write(stdout,'(a)') '... Core spheres initialised'

          end if
       end if

       ! HUBBARD PROJECTOR SPHERES
       if (pub_hubbard) then
          call function_basis_init_spheres(hub_proj_basis(isub), mdl%cell, &
               mdl%fftbox, par%elements_on_proc, par=par)
          if (pub_on_root .and. pub_output_detail >= VERBOSE) &
               write(stdout,'(a)') '... Hubbard Projector spheres initialised'
          ! ndmh: Hubbard Model allocation and setup
          call hubbard_model_init(hub,mdl%regions(isub)%elements)
       endif

       ! CONDUCTION NGWF SPHERES
       if (pub_cond_calculate .or. pub_lr_tddft_calculate) then
          call function_basis_init_spheres(cond_ngwf_basis(isub), mdl%cell, &
               mdl%fftbox, par%elements_on_proc, par=par)
          if (pub_on_root .and. pub_output_detail >= VERBOSE) &
               write(stdout,'(a)') '... Conduction NGWF spheres initialised'
          call function_basis_copy_spheres(joint_ngwf_basis(isub),ngwf_basis(isub),&
               cond_ngwf_basis(isub), par=par)!mdl%regions(isub)%par)
       end if

       ! AUXILIARY NGWF SPHERES (not currently used)
       if (pub_use_aux_ngwfs) then
          call function_basis_init_spheres(aux_ngwf_basis(isub), mdl%cell, &
               mdl%fftbox, par%elements_on_proc, par=par)
          if (pub_on_root .and. pub_output_detail >= VERBOSE) &
               write(stdout,'(a)') '... Auxiliary NGWF spheres initialised'
       end if

       ! TIGHT BOXES
       call function_basis_init_tight_boxes(ngwf_basis(isub),mdl%fftbox,mdl%cell, &
            par=par)
       if (pub_on_root .and. pub_output_detail >= VERBOSE) &
            write(stdout,'(a)') '... NGWF tight boxes initialised'

       ! PAW SPHERE TIGHT BOXES
       if (pub_paw.or.pub_realspace_projectors) then
          call function_basis_init_tight_boxes(proj_basis(isub),mdl%fftbox,mdl%cell,&
               par=par)
          if (pub_on_root .and. pub_output_detail >= VERBOSE) &
               write(stdout,'(a)') '... Projector tight boxes initialised'
          if (pub_paw.and.pub_eels_calculate) then
             call function_basis_init_tight_boxes(core_basis(isub),mdl%fftbox, &
                  mdl%cell, par=par)
             if (pub_on_root .and. pub_output_detail >= VERBOSE) &
                  write(stdout,'(a)') '... Core tight boxes initialised'
          end if
       end if

       ! CONDUCTION NGWF TIGHT BOXES
       if (pub_cond_calculate .or. pub_lr_tddft_calculate) then
          call function_basis_init_tight_boxes(cond_ngwf_basis(isub), &
               mdl%fftbox, mdl%cell, par=par)
          if (pub_on_root .and. pub_output_detail >= VERBOSE) &
               write(stdout,'(a)') '... Conduction NGWF tight boxes initialised'
          call function_basis_init_tight_boxes(joint_ngwf_basis(isub), &
               mdl%fftbox, mdl%cell, par=par)
          if (pub_on_root .and. pub_output_detail >= VERBOSE) &
               write(stdout,'(a)') '... Joint NGWF tight boxes initialised'
       end if

       ! AUXILIARY NGWF TIGHT BOXES (not currently used)
       if (pub_use_aux_ngwfs) then
          call function_basis_init_tight_boxes(aux_ngwf_basis(isub),mdl%fftbox, &
               mdl%cell, par=par)
          if (pub_on_root .and. pub_output_detail >= VERBOSE) &
               write(stdout,'(a)') '... Auxiliary NGWF tight boxes initialised'
       end if

       call comms_barrier

       ! cks: gather all tightboxes from procs and join them in a large array
       call function_basis_gath_all_tbs(ngwf_basis(isub), par=par)

       ! lr408: Gather all conduction tightboxes from procs and join them in a large array
       ! lr408: Pass both NGWF sets to ensure universal tightbox has correct size
       ! tjz07: Add lr_tddft option
       ! rc2013: Don't initialise uni_tightboxes until we've got all regions!
       if (pub_cond_calculate .or. pub_lr_tddft_calculate) then
          call function_basis_gath_all_tbs(cond_ngwf_basis(isub), par=par)
          call function_basis_gath_all_tbs(joint_ngwf_basis(isub), par=par)
       end if

       if (pub_use_aux_ngwfs) call function_basis_gath_all_tbs(aux_ngwf_basis(isub), &
            par=par)

       call comms_barrier
       if (pub_on_root .and. pub_output_detail >= VERBOSE) write(stdout,'(a)') &
            '... Tight-boxes gathered'

       ! Print out an estimate of the memory usage of the NGWF basis
       call function_basis_estimate_size(ngwf_basis(isub),local_size,global_size)
       if (pub_on_root .and. pub_output_detail >= VERBOSE) write(stdout,'(a,2i13)') &
            '... NGWF basis size (local,global): ', local_size, global_size

       ! lr408: Estimate conduction basis size as well
       if (pub_cond_calculate .or. pub_lr_tddft_calculate) then
          call function_basis_estimate_size(cond_ngwf_basis(isub),local_size, &
               global_size)
          if (pub_on_root .and. pub_output_detail >= VERBOSE) &
               write(stdout,'(a,2i13)') '... Conduction NGWF basis size &
               &(local,global): ', local_size, global_size
          call function_basis_estimate_size(joint_ngwf_basis(isub),local_size, &
               global_size)
          if (pub_on_root .and. pub_output_detail >= VERBOSE) &
               write(stdout,'(a,2i13)') '... Joint NGWF basis size &
               &(local,global): ', local_size, global_size
       end if

       ! cks: initialise nonlocal pseudpotential projectors in
       ! cks: complex fftbox reciprocal representation
       if (pub_any_nl_proj) then
          call pseudo_species_init_proj(nl_projectors(isub), &
               mdl%regions(isub)%pseudo_sp, mdl%regions(isub)%elements, par, &
               is_cmplx=loc_cmplx)
       end if
       if (pub_paw) then
          call paw_species_init_proj(nl_projectors(isub), &
               mdl%regions(isub)%paw_sp, mdl%regions(isub)%elements, par, &
               is_cmplx=loc_cmplx)
       end if

       ! Load and initialise core wavefunctions if required
       if (pub_paw.and.pub_eels_calculate) then
          call paw_species_init_core_wvfns(core_wvfns(isub), &
               mdl%regions(isub)%elements,mdl%regions(isub)%paw_sp,par)
       end if

       call comms_barrier
       if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
            write(stdout,'(a)') '... done'
       call services_flush

       if (pub_aug) then
          ! ndmh: Initialise PAW augmentation density box
          call augmentation_box_init(mdl%aug_box,mdl%cell,mdl%fine_grid, &
               nl_projectors(isub))
       end if
    end do

    ! rc2013: pass all NGWF bases to ensure universal tightbox has correct size
    ! EMBED_FIX: consider doing this with other function_basis routines to remove
    ! subsystem loops above
    if (pub_cond_calculate .or. pub_lr_tddft_calculate) then
       call function_basis_init_uni_tb(mdl%uni_tightbox,mdl%cell, &
            ngwf_basis, cond_ngwf_basis)
    else
       call function_basis_init_uni_tb(mdl%uni_tightbox,mdl%cell,ngwf_basis)
    end if

    ! ndmh: Initialise FFTbox, tightbox and aug box Fourier transforms
    call fourier_init_fftbox(mdl%fftbox)
    if (pub_aug) call fourier_init_box(mdl%aug_box)
    if (pub_tightbox_fft_coarse) call fourier_init_box(mdl%uni_tightbox)
    call services_flush

    ! jcap: before we enter the regions loop again, initialise
    ! remote_mod if we need to
    if (pub_use_swx) then
       ! rc2013: otherwise there's only 1 subsystem
       call remote_init(ngwf_basis(1), loc_cmplx)
    else if (pub_use_activeswx) then
       ! rc2013: this is only initialising for the active subsystem
       call remote_init(ngwf_basis(pub_active_region), loc_cmplx)
    end if

    do isub=1,mdl%nsub
       par => mdl%regions(isub)%par
       ! Convert projectors to real space (if required)
       ! agrecokpt: initialize projectors on grid with complex values
       ! if using complex NGWFs in the calculation
       ! we can probably always use real projectors on grid and then create
       ! complex copies to store the k-point dependent phase factor
       if (pub_any_nl_proj.and.pub_realspace_projectors) then
          ! agrecokpt: in TB method for k-points, the k-point dependent terms appear
          ! as phases that multiply the entire hamiltonian/overlap matrices; no
          ! need to evaluate projectors as k-dependent quantities at this stage
          select case (pub_kpoint_method)
             ! KP method
             case('KP')
                call projectors_create_real(proj_basis(isub),nl_projectors(isub), &
                     mdl%fftbox,mdl%cell, is_cmplx=loc_cmplx, kpt=kpt_cart)
             ! TB method
             case('TB')
                call projectors_create_real(proj_basis(isub),nl_projectors(isub), &
                     mdl%fftbox,mdl%cell, is_cmplx=loc_cmplx)

             case default
                call utils_abort('Illegal k-point method specified')
          end select
       end if

       if (pub_paw.and.pub_realspace_projectors) then
          ! agrecokpt: in TB method for k-points, the k-point dependent terms appear
          ! as phases that multiply the entire hamiltonian/overlap matrices; no
          ! need to evaluate projectors as k-dependent quantities at this stage
          select case (pub_kpoint_method)
             ! KP method
             case('KP')
                call projectors_create_real(proj_basis(isub),nl_projectors(isub), &
                     mdl%fftbox,mdl%cell, is_cmplx=loc_cmplx, kpt=kpt_cart)

             ! TB method
             case('TB')
                call projectors_create_real(proj_basis(isub),nl_projectors(isub), &
                     mdl%fftbox,mdl%cell, is_cmplx=loc_cmplx)

             case default
                call utils_abort('Illegal k-point method specified')
          end select
       end if

       !if (pub_paw.and.pub_realspace_projectors) &
       !     call projectors_create_real(core_basis,core_wvfns,mdl%fftbox,mdl%cell)

       if (pub_hubbard) then
          ! ddor: initialise hubbard projectors in
          !      complex fftbox reciprocal representation
          call hubbard_species_proj(hub,mdl%regions(isub)%elements,mdl%fftbox, &
               is_cmplx=loc_cmplx)

          if (pub_task == 'HUBBARDSCF' .or. pub_hubbard_restart &
               & .or. pub_hubbard_atomsolve .or. pub_hub_on_the_fly) then
             call function_basis_init_tight_boxes(hub_proj_basis(isub), &
                  mdl%fftbox,mdl%cell,ngwf_basis(isub), par=par)
             call function_basis_gath_all_tbs(hub_proj_basis(isub), &
                  ngwf_basis(isub), par=par)
             call hubbard_build_consist(hub, hub_proj_basis(isub), &
                  ngwf_basis(isub), loc_cmplx)
          endif
       endif

       ! jd: First stage of spherical wave resolution of identity init.
       !     This *must* precede the SPAM3 module initialisation below.
       ! jcap: if we only need spherical waves in the active
       ! subregion for an embedding calculation, only do this if we
       ! are in the active region
       if((pub_use_activeswx.and.(isub==pub_active_region)).or.pub_use_swx) then
          call swri_init_module_stage_1(mdl%regions(isub)%elements, mdl%cell, &
               mdl%uni_tightbox, mdl%regions(isub)%par)
       end if

       ! ja531: pdos init. This also *must* precede SPAM3 stuff.
       if(pub_dos_smear > 0.0_dp .and. pub_pdos_max_l >= 0) then
           par%max_pdos_n = ceiling(maxval(ngwf_basis(isub)%spheres(:)%radius) / &
               & min(mdl%cell%d1,mdl%cell%d2,mdl%cell%d3))
           call comms_reduce('MAX',par%max_pdos_n)
       end if

       ! ndmh: initialise SPAM3 module
       ! rc2013: pass par explicitly for initialisation
       ! rc2013: require initialisation of cross overlap terms
       ! rc2013: need a way of distinguishing subsystem sparse structures
       call sparse_init_blocking_scheme(BLKS_NGWF,ngwf_basis(isub)%num, &
            ngwf_basis(isub)%num_on_proc, ngwf_basis(isub)%num_on_atom, &
            ngwf_basis(isub)%first_on_proc, ngwf_basis(isub)%first_on_atom, &
            ngwf_basis(isub)%atom_of_func, ngwf_basis(isub)%proc_of_func, par)
       if (pub_any_nl_proj.or.pub_paw) then
          call sparse_init_blocking_scheme(BLKS_PROJ,proj_basis(isub)%num, &
               proj_basis(isub)%num_on_proc, proj_basis(isub)%num_on_atom, &
               proj_basis(isub)%first_on_proc, proj_basis(isub)%first_on_atom, &
               proj_basis(isub)%atom_of_func, proj_basis(isub)%proc_of_func,par)
       end if
       if (pub_hubbard) then
          call sparse_init_blocking_scheme(BLKS_HUB_PROJ, &
               hub_proj_basis(isub)%num, hub_proj_basis(isub)%num_on_proc, &
               hub_proj_basis(isub)%num_on_atom, &
               hub_proj_basis(isub)%first_on_proc, &
               hub_proj_basis(isub)%first_on_atom, &
               hub_proj_basis(isub)%atom_of_func, &
               hub_proj_basis(isub)%proc_of_func, par)
       end if
       if (pub_cond_calculate .or. pub_lr_tddft_calculate) then
          call sparse_init_blocking_scheme(BLKS_COND,cond_ngwf_basis(isub)%num, &
               cond_ngwf_basis(isub)%num_on_proc, cond_ngwf_basis(isub)%num_on_atom, &
               cond_ngwf_basis(isub)%first_on_proc,  &
               cond_ngwf_basis(isub)%first_on_atom, &
               cond_ngwf_basis(isub)%atom_of_func, &
               cond_ngwf_basis(isub)%proc_of_func, par)
          call sparse_init_blocking_scheme(BLKS_JOINT,joint_ngwf_basis(isub)%num, &
               joint_ngwf_basis(isub)%num_on_proc,  &
               joint_ngwf_basis(isub)%num_on_atom, &
               joint_ngwf_basis(isub)%first_on_proc, &
               joint_ngwf_basis(isub)%first_on_atom, &
               joint_ngwf_basis(isub)%atom_of_func, &
               joint_ngwf_basis(isub)%proc_of_func, par)
       end if
       if (pub_paw.and.pub_eels_calculate) then
          call sparse_init_blocking_scheme(BLKS_CORE,core_basis(isub)%num, &
               core_basis(isub)%num_on_proc, core_basis(isub)%num_on_atom, &
               core_basis(isub)%first_on_proc, core_basis(isub)%first_on_atom, &
               core_basis(isub)%atom_of_func, core_basis(isub)%proc_of_func, par)
       end if
       if (pub_use_aux_ngwfs) then
          call sparse_init_blocking_scheme(BLKS_AUX,aux_ngwf_basis(isub)%num, &
               aux_ngwf_basis(isub)%num_on_proc, aux_ngwf_basis(isub)%num_on_atom, &
               aux_ngwf_basis(isub)%first_on_proc, aux_ngwf_basis(isub)%first_on_atom, &
               aux_ngwf_basis(isub)%atom_of_func, aux_ngwf_basis(isub)%proc_of_func, par)
       end if

       ! ja531 --> Set up SW function basis / blocking scheme in eigenstates.
       if(pub_dos_smear > 0.0_dp .and. pub_pdos_max_l >= 0) then
          call eigenstates_init_pdos(mdl,ngwf_basis(isub), par=par)
       end if

       ! rc2013: initialise sparse and related features for this par
       call sparse_mod_init(par, first_sparse_call)
       call services_flush
    end do

    ! ndmh: allocate storage for local pseudopotential on fine grid
    allocate(mdl%localpseudo_fine(mdl%fine_grid%ld1, &
        mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12),stat=ierr)
    call utils_alloc_check(myself,'mdl%localpseudo_fine', ierr)

    ! ndmh: Allocate SPAM3 structures for overlap matrix, kinetic energy,
    ! ndmh: density kernel, inverse overlap
    call ngwf_rep_create(rep,'',mdl, is_cmplx=loc_cmplx)

    do isub=1,mdl%nsub
       ! ddor: initialise SPAM3 matrices for DFT+U if necessary
       if (pub_hubbard) then
           call hubbard_build_matrices(hub,hub_proj_basis(isub),is_cmplx=loc_cmplx)
       endif

       ! For constrained DFT calculations only, read the .cdft file and
       ! use it to set the constraining U potentials. Also update the
       ! cDFT overlap matrices
       ! NB: we are implicitly assuming that a single-point cDFT run has been
       !     already performed and the .cdft file is present]
       !gom
       if ((pub_cdft.or.pub_dft_nu).AND.(pub_maxit_cdft_u_cg>0).AND. &
                          ( (pub_cdft_continuation)       .OR. &
                          (pub_dft_nu_continuation)       .OR. &
                          (pub_task=='GEOMETRYOPTIMIZATION')  .OR. &
                          (pub_task=='MOLECULARDYNAMICS')     .OR. &
                          (pub_task=='TRANSITIONSTATESEARCH') .OR. &
                          (pub_task=='PROPERTIES')            .OR. &
                          (pub_task=='COND')                  .OR. &
                          (pub_task=='PROPERTIES_COND') )) then
                          ! pub_task=TDDFT to be potentially added...
          if(pub_cdft) then
              call cdft_intermediate_restart_read(hub,mdl)
              call cdft_intermediate_update_matrices(hub, hub_proj_basis(isub))
          elseif(pub_dft_nu) then
             call dft_nu_intermediate_restart_read(hub,mdl)
             call dft_nu_intermediate_update_matrices(hub, hub_proj_basis(isub))
          endif
       endif

       if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
            write(stdout,'(a)') ' done'

       ! ndmh: allocate storage for local pseudopotential on fine grid
       ! jcap: the regional localpseudo_fine arrays are now
       ! initialised in electronic_init, because that's the only place
       ! they are used, to save memory

    enddo
    ! rc2013: END OF LOOP: allocate grid data

    ! ndmh: allocate storage for core charge if NLCC is present
    if (pub_nlcc) then
       allocate(mdl%core_density_fine(mdl%fine_grid%ld1, &
            mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12),stat=ierr)
       call utils_alloc_check(myself, 'mdl%core_density_fine',ierr)
       ! jcap: do same for each subsystem
       do isub=1,mdl%nsub
          allocate(mdl%regions(isub)%core_density_fine(mdl%fine_grid%ld1, &
               mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12),stat=ierr)
          call utils_alloc_check(myself, 'mdl%regions%core_density_fine',ierr)
       end do ! isub
    else
       ! ndmh: no core charge required - allocate dummy storage
       allocate(mdl%core_density_fine(1,1,1),stat=ierr)
       call utils_alloc_check(myself, 'mdl%core_density_fine',ierr)
       ! jcap: do same for each subsystem
       do isub=1,mdl%nsub
          allocate(mdl%regions(isub)%core_density_fine(1,1,1),stat=ierr)
          call utils_alloc_check(myself, 'mdl%regions%core_density_fine',ierr)
       end do ! isub
    end if

    ! lr408: Allocate storage for slabs of local, hartree and XC pot
    allocate(lhxc_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,pub_num_spins,mdl%nsub),stat=ierr)
    call utils_alloc_check(myself,'lhxc_fine',ierr)

    if (pub_xc_ke_density_required) then
       ! JCW: If a tau-dependent XC functional is requested, then allocate space
       ! JCW: to store the gradient of the XC energy per unit volume wrt tau.
       allocate(dfdtau_fine(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check(myself,'dfdtau_fine',ierr)
    else
       ! JCW: If no tau-dependent XC functional is requested, then allocate a
       ! JCW: 1 element dummy array (this avoids the need for conditional
       ! JCW: subroutine calls and optional arguments)
       allocate(dfdtau_fine(1,1,1,1),stat=ierr)
       call utils_alloc_check(myself,'dfdtau_fine',ierr)
    end if

    ! ndmh: Setup and create storage for Hamiltonian
    call ngwf_ham_create(ham,rep)

    ! Electron counting & apply spin polarisation to get n_occ
    ! agrecokpt: occupation numbers are independent of k-point,
    ! hence set same value for all k-points with same spin channel
    ! need to consider if we rescale by kpoint weight here, or
    ! we do it later
    !do ik=1,pub_num_kpoints
    !    rep%n_occ(UP,ik) = ((nelecs+pub_spin)/2)*mdl%unique_kpoints(ik)%weight
    !    if (pub_num_spins > 1) then
    !       rep%n_occ(DN,:) = ((nelecs-pub_spin)/2)*mdl%unique_kpoints(ik)%weight
    !    end if
    !end do
    nelecs = energy_and_force_count_elecs(mdl%elements, mdl%nat)
    rep%n_occ(UP,PUB_1K) = (nelecs + pub_spin) / 2
    if(pub_num_spins > 1) then
       rep%n_occ(DN,PUB_1K) = (nelecs - pub_spin) / 2
    end if

    ! cks: Increase number of occupied states if doing hyper Hartree Fock,
    ! cks: and compute HHF scaling factor
    if (pub_hhf_nstates > 0) then
       pub_hhf_factor = 0.0_DP
       do is=1,pub_num_spins
          pub_hhf_factor = pub_hhf_factor + &
               real(rep%n_occ(is,PUB_1K),kind=DP)
       enddo
       if (pub_num_spins == 1) then
          pub_hhf_factor = 2.0_DP* pub_hhf_factor
       endif
       pub_hhf_factor = (pub_hhf_factor-1.0_DP) / &
            (pub_hhf_factor+2.0_DP*real(pub_hhf_nstates,kind=DP)-1.0_DP)
       rep%n_occ(:,:) = rep%n_occ(:,:) + pub_hhf_nstates
    endif

    ! mjsp : Store fragment electron counts (used in EDA
    ! supermolecule kernel polarisation methods)
    if (pub_eda .and. pub_eda_mode == EDA_ISOLATED) then

       ! Allocate n_occ according to public numbers of spins and k-points.
       allocate(pub_frag_data(pub_frag_counter)%rep%n_occ(pub_num_spins, &
            pub_num_kpoints), stat=ierr)
       call utils_alloc_check(myself, 'pfd(it)%rep%n_occ', &
            ierr)

       pub_frag_data(pub_frag_counter)%rep%n_occ = rep%n_occ

    end if

    ! rc2013: create full system density kernel
    call kernel_create(denskern, 'K'//trim(rep%postfix), is_cmplx=loc_cmplx)

    ! jd: Read denskern file
    ! jcap: This needs to be done for the whole system at once,
    ! so do it before we enter the loop over regions
    if (pub_is_separate_restart_files) then
       call restart_kernel_read(denskern%kern, read_solv = .true., &
            read_emft = pub_is_emft_cavity)
    end if

    ! JCW: Initialise spherical harmonic rotation module (should be done prior
    ! JCW: to initialization of evaluation of metric matrices in
    ! JCW: swri_init_module_stage_2). The initialization routine returns
    ! JCW: without doing anything if spherical harmonic rotations will not be
    ! JCW: needed.
    call sph_harm_rot_initialise

    ! rc2013: NEW LOOP
    do isub=1,mdl%nsub
       ! rc2013: get total n_occ
       !rep%n_occ(UP,PUB_1K,isub) = rep%n_occ(UP,PUB_1K,isub)
       par => mdl%regions(isub)%par
       ! rc2013: use new rep structure instead
       call data_functions_alloc(rep%ngwfs_on_grid(isub), &
             ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)

       ! jd: If spherical wave expansion in use, complete swri initialisation.
       !     This fills the metric matrices and does additional setup (NLs, etc).
       ! jcap: If we only need spherical waves in the active
       ! subregion, only do it there
       !     Initialise swex module.
       if((pub_use_activeswx.and.(isub==pub_active_region)).or.pub_use_swx) then
          call swri_init_module_stage_2(mdl%regions(isub)%elements, &
               rep%overlap%m(isub,isub), mdl%cell)
          call swx_init_module(rep)
       end if
    ! rc2013: EMBED_FIX!
    !! ndmh: Allocate storage for NGWFs
    !!allocate(rep%ngwfs_on_grid(ngwf_basis%size_on_grid),stat=ierr)
    !!call utils_alloc_check(myself,'rep%ngwfs_on_grid', &
    !!     ierr)
    !! agrecocmplx: use appropriate routine to allocate ngwfs_on_grid
    !call data_functions_alloc(rep%ngwfs_on_grid, ngwf_basis%size_on_grid, &
    !     iscmplx=loc_cmplx)

    !! JCW: Initialise spherical harmonic rotation module (should be done prior
    !! JCW: to initialization of evaluation of metric matrices in
    !! JCW: swri_init_module_stage_2). The initialization routine returns
    !! JCW: without doing anything if spherical harmonic rotations will not be
    !! JCW: needed.
    !call sph_harm_rot_initialise

    !! jd: If spherical wave expansion in use, complete swri initialisation.
    !!     This fills the metric matrices and does additional setup (NLs, etc).
    !!     Initialise swex module.
    !if (pub_use_swx) then
    !   call swri_init_module_stage_2(mdl%elements, rep%overlap)
    !   call swx_init_module(rep)
    !end if

       ! jd: Set up polarisable embedding, if needed.
       if (pub_pol_emb_pot) call polarisable_embedding_init( &
            ngwf_basis, mdl, rep%n_occ, rep%overlap)!, denskern)

       ! jd: Initialise HFx, if needed.
       ! jcap: again, only do the active subregion if required
       ! jcap: this also works if we are doing a full system calculation
       if (pub_use_hfx.or.(pub_use_activehfx.and.(isub==pub_active_region))) then
          call hf_exchange_init(hfxstate, rep, ham%hfexchange(1)%m(isub,isub),&
               mdl%cell, size(mdl%regions(isub)%elements))
       end if

    end do

    ! jd: Establish if NGWFs are going to be optimised. Currently this informa-
    !     tion is only needed to inform the HFx engine on how to optimise SWOP
    !     caching (whether to optimise for the calculation of NGWF gradient).
    !     This needs to be done before the call to ngwfs_initialise at [**].
    !     The logic below reproduces the logic later on, each section of the
    !     sections like [A1*] below has a matching section like [B1*] later on.
    pub_ngwf_gradient_needed = .false.
    if ((pub_task/='PROPERTIES').and.(pub_task/='COND').and.(pub_task/='HUBBARDSCF').and. &
         (pub_task/='TDDFT').and.(pub_task/='PROPERTIES_COND').and. &
         (pub_task/='LR_TDDFT').and.(pub_task/='LR_PHONONS').and. &
         (pub_task/='LUMOSEARCH').and. (pub_task/='COUPLINGS').and. &
         .not.(((pub_task=='EDA') .or. (pub_task=='EDA_PREP')) &
         .and. pub_eda_mode/=EDA_ISOLATED)) then                      ! [A1*]
       if ((pub_maxit_ngwf_cg .ge. 0) .or. &
            (pub_nonsc_forces .and. pub_write_forces)) then
          pub_ngwf_gradient_needed = .true.
       end if
    end if
    if (pub_task == 'HUBBARDSCF' .and. pub_hub_max_iter >=1) then     ! [A2*]
       pub_ngwf_gradient_needed = .true.
    end if
    if (pub_task == 'COND') then ! [A3*] cf. cond_properites_only in conduction
       pub_ngwf_gradient_needed = .true.
    end if
    if (pub_eda .and. pub_eda_mode/=EDA_ISOLATED) then ! [A4*]
       if ((pub_maxit_ngwf_cg .gt. 0).or.pub_nonsc_forces) then
          pub_ngwf_gradient_needed = .true.
       end if
    end if

    ! jcap: do these bits on a system-wide basis
    ! jd: This is a good place to initialise smeared ions, if in use
    if (pub_is_smeared_ion_rep) call smeared_ion_initialise(mdl, &
         pub_smeared_ion_bc_is_periodic)

    ! jd: Initialise Boltzmann solvation's ion-position-dependent state
    !     (reports on Boltzmann ions, calculates steric potential, calls
    !     dl_mg_init_nonlin). This is done here, as when ions move, the steric
    !     potential must be recalculated and dl_mg must be made aware.
    if (pub_is_pbe /= 'NONE') then
       call implicit_solvent_boltzmann_init(mdl)
    end if

    ! jd: Also, solvation must know they were re-initialised. Annoying
    !     vestige of Hatem's approach ('poisson_problem' structure).
    if (pub_is_smeared_ion_rep) have_rho_ion = .false.

    ! jd: If separate restart files are needed to construct the solvent cavity,
    !     deal with these now.
    if (pub_is_separate_restart_files) then
       allocate(temp_density_fine(mdl%fine_grid%ld1, &
            mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check(myself,'temp_density_fine', ierr)
       allocate(sub_temp_density_fine(mdl%fine_grid%ld1, &
            mdl%fine_grid%ld2, mdl%fine_grid%max_slabs12,pub_num_spins),stat=ierr)
       call utils_alloc_check(myself,'sub_temp_density_fine', ierr)
       ! jcap: loop over regions to read in NGWFs
       do isub=1,mdl%nsub

          ! jcap: add region string to name if required
          if (mdl%nsub.gt.1) then
             write(vac_ngwfs_ext,'(a,i1)') 'vacuum_tightbox_ngwfs_',isub
          else
             vac_ngwfs_ext = 'vacuum_tightbox_ngwfs'
          end if

          ! jd: Read tightbox NGWFs
          call restart_ngwfs_tightbox_input(rep%ngwfs_on_grid(isub), &
               ngwf_basis(isub), mdl%cell, mdl%fftbox, mdl%regions(isub)%elements, &
               trim(vac_ngwfs_ext),mdl%regions(isub))
       end do
       call ngwf_rep_register_change(rep,'restart from vacuum tightbox NGWFs')
       call ngwf_ham_register_change(ham,rep)
       call restart_kernel_read(denskern%kern, read_solv = .true., &
            read_emft = pub_is_emft_cavity)
       ! jcap: if using PAW, allocate nhat_den_grad
       if (pub_aug) then
          allocate(nhat_den_grad(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
               mdl%fine_grid%max_slabs12, pub_num_spins, 0:pub_aug_den_dim), &
               stat=ierr)
          call utils_alloc_check(myself,'nhat_den_grad', ierr)
       end if
       ! jcap: now loop over regions again
       temp_density_fine=0.d0
       ! jcap: need to calculate the overlap matrix, but if we're
       ! using the EMFT kernel and we need to block orthogonalise,
       ! need to do this by calling hamiltonian_dens_indep_matrices
       if (pub_block_orthogonalise.and.pub_is_emft_cavity) then
          temp_build_bo=pub_build_bo
          pub_build_bo=.true.
          call hamiltonian_dens_indep_matrices(rep,ngwf_basis,proj_basis,nl_projectors,&
               hub_proj_basis,hub,mdl)
          pub_build_bo=temp_build_bo
       end if
       do isub=1,mdl%nsub
          ! jd: Calculate overlap matrix
          ! jcap: This needs to be done between all possible pairs of regions
          ! jcap: if it hasn't already been done because we've block orthogonalised
          do jsub=1,mdl%nsub
             if (.not.(pub_block_orthogonalise.and.pub_is_emft_cavity)) then
                call function_ops_brappd_ketppd(rep%overlap%m(isub,jsub), &
                     rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                     rep%ngwfs_on_grid(jsub), ngwf_basis(jsub), mdl%cell)
             end if
             ! jd: Calculate density on fine grid
             call sparse_embed_extract_from_array(kern_array,&
                  denskern%kern%m(:,PUB_1K),isub,jsub)
             call density_on_grid(sub_temp_density_fine, &
                  mdl%fine_grid, mdl%dbl_grid, mdl%cell, mdl%fftbox, &
                  kern_array, rep%ngwf_overlap%m(isub,jsub), &
                  rep%ngwfs_on_grid(isub), ngwf_basis(isub), &
                  rep%ngwfs_on_grid(jsub), ngwf_basis(jsub))
             temp_density_fine=temp_density_fine+sub_temp_density_fine
             ! ndmh: for PAW, get augmentation density nhat
             if (pub_aug) then
                nhat_den_grad = 0.0_DP
                call projectors_func_ovlp_box(rep%sp_overlap%m(isub,jsub), &
                     rep%ngwfs_on_grid(isub),ngwf_basis(isub), &
                     proj_basis(jsub),nl_projectors(jsub), mdl%fftbox,mdl%cell)
                ! jcap: PAW still only compatible with one subregion
                call augmentation_density_on_grid(nhat_den_grad, &
                     mdl%fine_grid,mdl%cell,mdl%regions(isub)%pseudo_sp,&
                     mdl%regions(isub)%paw_sp,mdl%aug_box, &
                     kern_array,rep%sp_overlap%m(isub,jsub))
                temp_density_fine = temp_density_fine + nhat_den_grad(:,:,:,:,0)
             end if
             call sparse_embed_destroy_extracted_array(kern_array)
          end do
       end do
       ! jcap: deallocate nhat_den_grad
       if (pub_aug) then
          deallocate(nhat_den_grad,stat=ierr)
          call utils_dealloc_check(myself, 'nhat_den_grad', ierr)
       end if

       ! jd: Take spin into account
       temp_density_fine = temp_density_fine * pub_spin_fac
       ! jd: Initialise and destroy a dummy solvation problem, this is enough (!)
       !     for the initial eps and eps_half to be calculated.
       have_rho_ion = .true.
       call initialize_solvation_problem(temp_density_fine, &
            dummy_solv_problem,mdl%fine_grid,mdl%cell,mdl%elements)
       call deallocate_poisson_problem(dummy_solv_problem,'DIRECT')
       have_rho_ion = .false.
       deallocate(temp_density_fine,stat=ierr)
       call utils_dealloc_check(myself,'temp_density_fine', ierr)
       deallocate(sub_temp_density_fine,stat=ierr)
       call utils_dealloc_check(myself,'sub_temp_density_fine', ierr)
    end if


    ! vv: Restore electronic_history now
    if (md_global_restart) then
       call elec_history_restore_history_parm(mdl%elements, ngwf_basis, &
            rep%ngwfs_on_grid(1)%iscmplx, mdl%cell, mdl%fftbox)
       md_global_restart = .false.
    end if

    if (pub_use_aux_ngwfs) call internal_create_aux_rep

    ! mjsp: EDA frozen density calculation
    if (pub_eda_scfmi .or. (pub_eda .and. (pub_eda_mode == EDA_FROZENNONIDEM))) then

       ! mjsp: Initialise supermolecule NGWFs and the EDA driver
       ! mjsp: supermolecule-fragment data arrays
       call eda_driver_prep_supermol(mdl, hfxstate, rep, ngwf_basis)

       ! mjsp: If intrafragment partitioning involved
       ! then apply the modifications to the fragment data to enable this:
       if (pub_eda_split_atoms) &
          call fmo_apply_intrafrag_partitions()

       ! mjsp: allocate the Stoll SCF MI matrices
       call scfmi_allocate_matrices(denskern,rep)

    else

       ! ndmh: Now actually initialise the NGWF values (either from PAOs, STO-3G,
       ! ndmh: or by loading from disk (or memory)
       if (pub_task/='LUMOSEARCH') then
          call ngwfs_initialise(rep,ngwf_basis,mdl,hfxstate,&
               restart_rootname=custom_restart_name) ! [**]
       end if

    end if

    tot_num=0
    do isub=1,mdl%nsub
       tot_num=tot_num+ngwf_basis(isub)%num
    end do

    call comms_barrier
    call services_flush

    call comms_barrier

    ! ars: create EDFT model if required
    ! jcap: for the full system
    if ((pub_edft).or.(pub_edft_init_maxit>0)) then
       call edft_create(edft, ham, tot_num, &
            energy_and_force_count_elecs(mdl%elements, mdl%nat))
    end if
    ! jd: Initialise DFTB, if needed
    if(pub_dftb) then
       call utils_assert(mdl%nsub == 1, &
            'DFTB does not support embedding regions')
       dftb_disp_energy = pub_dispersion_energy
       call dftb_init_stage_2(mdl, rep, ngwf_basis, ham)
       call dftb_rep(mdl%dftb_par, mdl%regions(1)%elements, par, &
            mdl%cell, dftb_rep_energy, F_dftb_rep)
       call dftb_srb(mdl%dftb_par, mdl%regions(1)%elements, par, &
            mdl%cell, mdl%dftb_var%dcdr(1:3), dftb_srb_energy, F_dftb_srb)
       call dftb_ies(mdl, dftb_ies_energy, F_dftb_ies)
    end if

    ! LUMOSEARCH: Enter driver for finding LUMO in linear-scaling effort
    if (pub_task=='LUMOSEARCH') then
       call eigenstates_lumo_search(rep, ham, denskern, ngwf_basis, mdl, &
            hfxstate)
    end if

    ! jcap: In an TD-EMFT calculation, we need to make sure pub_emft
    ! is true and emft_follow is false, as we want to include the
    ! active functional at all points
    temp_emft_follow=.false.
    temp_emft=.false.
    if (pub_task=='LR_TDDFT') then
       if (pub_emft.or.pub_emft_follow) then
          if (pub_on_root) write(stdout, '(a)') &
               'Forcing pub_emft = T, pub_emft_follow = F for LR-TDDFT calculation'
          temp_emft_follow=pub_emft_follow
          temp_emft=pub_emft
          pub_emft_follow=.false.
          pub_emft=.true.
          ! jcap: We also want to block orthogonalise our matrices if required
          if (pub_block_orthogonalise) then
             temp_build_bo=pub_build_bo
             pub_build_bo=.true.
          end if
       end if
    end if

    ! cks: initialise density kernel, local pseudopotential on fine grid
    ! cks: and density-independent matrices
    ! ebl: only initialise density kernel if DMFT store file not present,
    !      and pub_dmft_read is true.
    check_dmft_store=.false.
    if(pub_dmft_points>0 .and. pub_dmft_read) then
       inquire(file=trim(pub_rootname)//'.ham1',exist=check_dmft_store)
    end if
    if (.not. check_dmft_store .and. pub_task/='LUMOSEARCH') then
       call electronic_init_denskern(rep, ham, denskern, &
            lhxc_fine, ngwf_basis, proj_basis, nl_projectors, &
            hub_proj_basis, hub, mdl, hfxstate, edft_for_init=edft, &
            kpt=kpt_cart, dfdtau_fine = dfdtau_fine, &
            restart_rootname=custom_restart_name)
       ! JCW: At present (01/2016), the initial guess for KE density is 0.0_DP
       ! JCW: for all grid points, so dfdtau_fine should also be 0.0_DP at all
       ! JCW: grid points on exiting electronic_init_denskern.
    end if

    ! ndmh: if this is not an EDFT run, but we used EDFT for the
    ! ndmh: initialisation, clean up the EDFT container now.
    if ((.not.pub_edft).and.(pub_edft_init_maxit>0)) then
       call edft_destroy(edft)
    end if

    call comms_barrier

    ! ars: allocate NGWF non self-consistent forces only if required
    if (pub_nonsc_forces) then
       allocate(ngwf_nonsc_forces(1:3,mdl%nat),stat=ierr)
       call utils_alloc_check(myself, 'ngwf_nonsc_forces', ierr)
       ngwf_nonsc_forces(:,:) = -999_DP ! ndmh: init to flagged value
    else
       allocate(ngwf_nonsc_forces(0,0),stat=ierr)
       call utils_alloc_check(myself,'ngwf_nonsc_forces', ierr)
    end if

    ! rc2013: set embedding optimisation parameters if needed
    if (mdl%nsub == 1) then
       pub_do_fandt = .false.
       pub_quit_region = .false.
       pub_freeze_switch_steps = -1
    else if (pub_do_fandt) then
       ! rc2013: don't quit region for F+T tests
       ! quitting is the "right" way to do this (since we shouldn't take any
       ! further steps in this region once we're below the convergence threshold)
       ! but would be inconsistent with old F+T runs
       pub_quit_region = .false.
    else if (pub_emft_follow) then
       ! rc2013: do 2 NGWF optimisations: 1 at low level, one at EMFT level
       pub_quit_region = .true.
    else
       pub_quit_region = .false.
       pub_freeze_switch_steps = -1
    endif

    ! rc2013: print header for embedding calculation
    !if((mdl%nsub) .gt. 1) then
    !   if(pub_on_root) write(stdout,'(a)') &
    !        utils_banner('=', 'Begin embedding optimisation')
    !endif

    ! ndmh: standard valence NGWF optimisation (used for most tasks)
    if ((pub_task/='PROPERTIES').and.(pub_task/='COND').and.(pub_task/='HUBBARDSCF').and. &
         (pub_task/='TDDFT').and.(pub_task/='PROPERTIES_COND').and.&
         (pub_task/='LR_TDDFT').and.(pub_task/='LR_PHONONS').and. &
         (pub_task/='LUMOSEARCH').and. (pub_task/='COUPLINGS').and. &
         .not.(((pub_task=='EDA') .or. (pub_task=='EDA_PREP')) &
         .and. pub_eda_mode/=EDA_ISOLATED)) then ! [B1*]
         ! jd: ^^ PLEASE REFLECT ANY CHANGES TO THE CONDITIONS HERE IN AN
         !     ANALOGOUS BLOCK [A1*] SOME 250 LINES ABOVE
       ! optimise density kernel and NGWF coefficients
       if ((pub_maxit_ngwf_cg .gt. 0).or.pub_nonsc_forces.or.pub_dftb) then
          call ngwf_cg_optimise(total_energy, ngwfs_converged, &
               out_of_runtime, ham, denskern, edft, rep, ngwf_nonsc_forces, &
               lhxc_fine, ngwf_basis, proj_basis,nl_projectors, &
               hub_proj_basis, hub, mdl, hfxstate, kpt=kpt_cart, &
               dfdtau_fine = dfdtau_fine)
       end if
    end if

    if (pub_task == 'HUBBARDSCF') then
       ! jd: ^^ PLEASE REFLECT ANY CHANGES TO THE CONDITION HERE IN AN ANALOGOUS
       !        BLOCK [A2*] SOME 250 LINES ABOVE

       ! ddor: DFT+U energy minimisation with self-consistent NGWF projectors.
       ! JCW: Abort if tau-dependent XC functional requested
       call utils_assert(.not.pub_xc_ke_density_required, myself//": &
            &pub_xc_ke_density_required is true &
            &but a Hubbard SCF calculation has been requested. This &
            &combination has not been implemented.")
       ! Loop over Hubbard projector optimisation counter
       do hub_proj_iteration=1,pub_hub_max_iter
          if (pub_on_root) then
             write(stdout,'(/a/)') repeat('#',80)
             write(stdout,'(a,i5)') 'HUBBARD PROJECTOR SCF ITERATION ', &
                  hub_proj_iteration
             write(stdout,'(/a/)') repeat('#',80)
          endif
          hub%consistency_iteration = hub_proj_iteration ! public
          ! ddor: Depending on the projector iteration we are on,
          !       carry out some projector and metric mixing,
          !       and generate the new projector-NGWF overlap matrix
          ! rc2013: not ready for embedding
          call hubbard_projector_update( &
               rep%ngwfs_on_grid(1), ngwf_basis(1), nl_projectors(1), &
               proj_basis(1), hub_proj_basis(1), hub, rep%inv_overlap%p, &
               rep%hub_overlap%p, rep%hub_overlap_t%p, &
               rep%sp_overlap%p, rep%hub_proj_paw_overlap%p, &
               rep%hub_ngwf_paw_overlap%p, mdl)

          ! optimise density kernel and NGWF coefficients
          ! agrecokpt: call the routine with the appropriate k-point
          call ngwf_cg_optimise(total_energy, ngwfs_converged, out_of_runtime, &
               ham, denskern, edft, rep, ngwf_nonsc_forces, lhxc_fine, &
               ngwf_basis, proj_basis, nl_projectors, hub_proj_basis, hub, mdl,&
               hfxstate, kpt=kpt_cart)

          ! ddor: Call to Hubbard projector optimisation
          !       energy convergence criterion
          if (hubbard_test_convergence(hub,total_energy)) exit
       end do
    end if

    if ((pub_task == 'COND').or.(pub_task == 'PROPERTIES_COND')) then
       ! jd: ^^ PLEASE REFLECT ANY CHANGES TO THE CONDITIONS HERE IN AN
       !     ANALOGOUS BLOCK [A3*] SOME 300 LINES ABOVE

       ! lr408: Call conduction driver routine
       ! agrecokpt: at specified k-point?
       call conduction_ngwf_optimise(total_energy, denskern, rep, ham, &
            ngwf_basis, proj_basis, nl_projectors, hub_proj_basis, hub, edft, &
            cond_ngwf_basis, joint_ngwf_basis, core_basis, core_wvfns, mdl, &
            hfxstate, lhxc_fine, cond_properties_only=(pub_task=='PROPERTIES_COND'), &
            kpt=kpt_cart)
    end if

    if (((pub_task=='PROPERTIES').or.(local_properties).or.(pub_do_properties) &
         .or.(pub_md_properties)).and.(pub_task/='COND').and. &
         (pub_task/='PROPERTIES_COND').and.(pub_task/='LR_TDDFT').and. &
         (pub_task/='LR_PHONONS').and. (pub_task/='COUPLINGS').and. &
         (.not.out_of_runtime)) then
       if (pub_use_aux_ngwfs) then
          ! JCW: Abort if KE-density dependent XC functional, since
          ! JCW: compatibility with auxiliary NGWFs is untested
          call utils_assert(.not.pub_xc_ke_density_required,myself//": &
               &pub_xc_ke_density_required is true &
               &but auxiliary NGWFs have been requested (block species_aux). &
               &This combination has not been tested.")
          ! jhl52: Properties, e.g. energy and occupancy, calculated using val
          ! Hamiltonian that is diagonalised in the AUX-NGWFs set
          call properties_calculate(aux_denskern, aux_ham, &
               aux_rep, aux_ngwf_basis, proj_basis, nl_projectors, &
               hub_proj_basis, hub, core_basis, core_wvfns, mdl, hfxstate, &
               lhxc_fine, 'aux')
       else
          ! jd: If properties-DMA requested, perform an SW expansion now,
          !     while 'rep' is still writable.
          if(pub_dma_calculate) then
             ! jcap: loop over regions
             do isub=1,mdl%nsub
                call dma_expand_ngwf_pairs(rep, ngwf_basis(isub), mdl, isub, &
                     REP_SWEX_PROPERTIES_DMA_1)
             end do
          end if
          ! jd: [*] If auto-solvation is in effect and we are in the vacuum
          !     stage, modify pub_rootname temporarily so that files output
          !     during properties are not overwritten later in the solvent stage.
          if (pub_is_auto_solvation .and. .not. pub_is_implicit_solvent) then
             old_pub_rootname = pub_rootname
             pub_rootname = trim(pub_rootname)//'_vacuum'
          end if
          ! cks: calculate electronic properties
          call properties_calculate(denskern, ham, &
               rep, ngwf_basis, proj_basis, nl_projectors, hub_proj_basis, hub, &
               core_basis, core_wvfns, mdl, hfxstate, lhxc_fine, 'valence',&
               dfdtau_fine=dfdtau_fine)
          ! jd: Undo potential change from [*]
          if (pub_is_auto_solvation .and. .not. pub_is_implicit_solvent) then
             pub_rootname = old_pub_rootname
          end if
          ! jd: Clean up after properties-DMA's SW expansion
          if(pub_dma_calculate) then
             call dma_free_swexes(rep%swexes(&
                  REP_SWEX_PROPERTIES_DMA_1:REP_SWEX_PROPERTIES_DMA_2))
          end if
       end if
    end if

    if ((pub_task == 'TDDFT').or.(pub_do_tddft)) then
       ! ddor: calculate dielectric response from TDDFT
       if(mdl%nsub .gt. 1) call utils_abort('ERROR: "TDDFT" &
            &with quantum embedding has not been implemented/tested &
            &yet. Aborting calculation.')
       call tddft_calculate(denskern%kern, rep, ngwf_basis(1), proj_basis(1), mdl)
    end if

    if (((pub_task == 'LR_TDDFT') .or. (pub_lr_tddft_calculate)) .and. &
         (pub_task /= 'COUPLINGS')) then
       ! tjz07: Call main linear response TDDFT routine
       if (pub_use_aux_ngwfs) then ! jhl52: aux option for QNTO analysis
          ! jcap: auxiliary NGWFs not set up with embedding yet
          ! jcap: ensure idempotency of kernels before putting it in
          call lr_tddft_calc_excitations(denskern%kern, &
               rep, ngwf_basis, cond_ngwf_basis, &
               mdl, hfxstate, lhxc_fine, proj_basis, hub_proj_basis, &
               total_energy, hub, nl_projectors, joint_ngwf_basis, &
               aux_rep, aux_ngwf_basis, aux_denskern%kern)
       else
          ! jcap: ensure idempotency of val kernel before putting it in
          call lr_tddft_calc_excitations(denskern%kern, &
               rep, ngwf_basis, cond_ngwf_basis,&
               mdl, hfxstate, lhxc_fine, proj_basis, hub_proj_basis, &
               total_energy, hub, nl_projectors, joint_ngwf_basis)
       endif ! aux
    end if

    if ((pub_task == 'LR_PHONONS') .or. (pub_lr_phonons_calculate)) then
       call lr_phonons_calculate(denskern%kern, rep, ngwf_basis, &
            hub_proj_basis, proj_basis, nl_projectors, mdl, hfxstate, &
            lhxc_fine(:,:,:,:,1), hub)
    end if

    ! cks: calculate forces if required
    if ((pub_write_forces .or. pub_task == 'GEOMETRYOPTIMIZATION' .or. &
          pub_task == 'TRANSITIONSTATESEARCH' .or. pub_task == 'MOLECULARDYNAMICS' .or. &
          pub_task == 'PHONON').and..not.out_of_runtime) then
       ! ab: forces in DFTB:
       if (pub_dftb) then
          ! ab: for dispersion use vdw D2 model
          if (pub_dispersion /= '0') then
             call vdwcorrection_calculate_forces(F_dftb_dis,mdl%elements,mdl%cell,mdl%par)
          else
             F_dftb_dis = 0.0_DP
          end if
          ! ab: calculate electronic forces
          call dftb_electronic_force(ngwf_basis(1),denskern%kern,mdl,rep, &
               edft%h_evals(:,:),edft%occ(:,:),edft%mo(:),edft%num,F_dftb_eht)

          total_forces = F_dftb_rep + F_dftb_srb + F_dftb_ies + F_dftb_dis + &
               F_dftb_eht
       else
          ! ndmh: Calculate ionic forces
          call forces_calculate(total_forces, &
               denskern, ham, lhxc_fine, rep, &
               ngwf_basis,proj_basis,nl_projectors,hub_proj_basis,hub, &
               mdl, hfxstate, ngwf_nonsc_forces,dfdtau_fine=dfdtau_fine)
       end if
    end if

    if (pub_dftb) then

       dftb_eht_energy = total_energy
       total_energy = dftb_eht_energy + dftb_rep_energy + &
            dftb_srb_energy + dftb_ies_energy + dftb_disp_energy

       ! ab: print qc info
       if (pub_print_qc) then
          call dftb_qc_print('dftb_electronic',dftb_eht_energy,mdl%nat,F_dftb_eht)
       end if

    end if


    ! dhpt: calculate couplings if required
    if (pub_task == 'COUPLINGS') then
       call utils_assert(mdl%nsub == 1, 'Error: "COUPLINGS" task requested &
            &but this is incompatible with quantum embedding. Please &
            &Please only use one subregion when doing a Couplings calculation.')
       call couplings_calculate(ngwf_basis(1), rep, hub, mdl, &
            hub_proj_basis(1), cond_ngwf_basis(1))
    end if

    ! mjsp: If not an isolated fragment EDA calculation:
    if (pub_eda .and. pub_eda_mode/=EDA_ISOLATED) then
       ! jd: ^^ PLEASE REFLECT ANY CHANGES TO THE CONDITIONS HERE IN AN
       !     ANALOGOUS BLOCK [A4*] SOME 400 LINES ABOVE
       ! mjsp: EDA supermolecule driver
       call eda_driver_supermol_run(mdl, hfxstate, ngwf_basis, aux_ngwf_basis, &
            cond_ngwf_basis, joint_ngwf_basis, proj_basis, &
            hub_proj_basis, hub, core_basis, core_wvfns, rep, ham, denskern, &
            edft, nl_projectors, ngwf_nonsc_forces, lhxc_fine(:,:,:,:,1))
    end if  ! if pub_eda

    ! jcap: if we have been doing an EMFT follow calculation, we need
    ! to reset the parameters in case further NGWF optimisations are
    ! required (e.g. conduction)
    if (pub_emft_follow) then
       call ngwf_cg_emft_reset
    end if

    ! jcap: if we have been doing a TD-EMFT calculation, reset
    ! pub_emft and emft_follow back to what they were
    if (temp_emft.or.temp_emft_follow) then
       if (pub_on_root) write(stdout, '(a)') &
            'Resetting pub_emft, pub_emft_follow to previous value'
       pub_emft_follow=temp_emft_follow
       pub_emft=temp_emft
       ! jcap: similar for build_bo
       if (pub_block_orthogonalise) pub_build_bo = temp_build_bo
    end if

    ! -------------------------------------------------------------------------
    ! ------------------------------- CLEANUP ---------------------------------
    ! -------------------------------------------------------------------------

    ! Deallocate
    if (pub_use_aux_ngwfs) call internal_destroy_aux_rep

    ! ars: deallocate EDFT container
    if (pub_edft) call edft_destroy(edft)

    ! rc2013: now let's destroy all the subsystem structures
    do isub=1,mdl%nsub
       ! ndmh: deallocate storage for NGWFs
       call data_functions_dealloc(rep%ngwfs_on_grid(isub))
    end do

    ! ndmh: deallocate storage for Hamiltonian, and whole-cell arrays
    call ngwf_ham_destroy(ham)
    deallocate(ngwf_nonsc_forces,stat=ierr)
    call utils_dealloc_check('energy_and_force_calculate','ngwf_nonsc_forces', &
         ierr)
    ! JCW: Deallocate storage for dfdtau_fine --- always allocated, since if a
    ! JCW: non-tau-dependent XC functional is used, then a single element
    ! JCW: dummy array is allocated.
    deallocate(dfdtau_fine,stat=ierr)
    call utils_dealloc_check(myself,'dfdtau_fine',ierr)
    deallocate(lhxc_fine,stat=ierr)
    call utils_dealloc_check(myself,'lhxc_fine',ierr)
    deallocate(mdl%core_density_fine,stat=ierr)
    call utils_dealloc_check(myself,'mdl%core_density_fine',ierr)
    ! jcap: do the same for individual regions
    do isub=1,mdl%nsub
       ! jcap: need to check if this is still allocated for this
       ! region, as most will be deallocated immediately after use in
       ! electronic_init_pot, as they aren't used elsewhere (exception
       ! for the active region in embedding)
       if (allocated(mdl%regions(isub)%core_density_fine)) then
          deallocate(mdl%regions(isub)%core_density_fine,stat=ierr)
          call utils_dealloc_check(myself,'mdl%regions%core_density_fine',ierr)
       end if
    end do

    deallocate(mdl%localpseudo_fine,stat=ierr)
    call utils_dealloc_check(myself,'mdl%localpseudo_fine', ierr)

    ! cks: memory deallocation for sparse matrices
    call kernel_destroy(denskern)
    call ngwf_rep_destroy(rep)

    ! rc2013: now let's destroy all the subsystem structures
    do isub=1,mdl%nsub

       ! ndmh: deallocation for core wavefunctions if required
       if (pub_paw.and.pub_eels_calculate) then
          call paw_species_exit_core_wvfns(core_wvfns(isub), &
               mdl%regions(isub)%paw_sp)
       end if
    end do

    ! ddor: DFT+U deallocation if necessary
    if (pub_hubbard) then
       call hubbard_build_matrices_exit(hub)
       ! ddor: deallocate DFT+U projector information
       call hubbard_species_exit_proj(hub)
       ! ddor: deallocate hubbard_atoms type array in DFT+U
       call hubbard_model_exit(hub)
    endif

    ! mjsp: SCF-MI polarisation or EDA supermolecule calculation:
    if ((pub_eda_scfmi) .or. ((pub_eda) .and. &
         ((pub_eda_mode == EDA_FULL .or. &
          (pub_eda_mode == EDA_FROZENNONIDEM .and. (pub_eda_preptool)))))) &
       call scfmi_deallocate_matrices()

    do isub=1,mdl%nsub
       ! pdh: clean up sparse module
       call sparse_exit(mdl%regions(isub)%par)
    end do

    ! jd: Clean up Hartree-Fock exchange, if in use
    ! jcap: including if it's just in the active subregion
    if (pub_use_hfx .or. pub_use_activehfx) call hf_exchange_free(hfxstate)

    ! jd: Clean up polarisable embedding module
    if (pub_pol_emb_pot) call polarisable_embedding_exit

    ! jd: Clean up spherical wave expansion and resolution of identity modules
    ! jcap: including if it's just in the active subregion
    if (pub_use_swx.or.pub_use_activeswx) then
       call remote_exit
       call swx_cleanup_module
       call swri_cleanup_module
    end if

    ! JCW: Clean up spherical harmonic rotation module (routine returns
    ! JCW: immediately if there is nothing to be done)
    call sph_harm_rot_exit

    ! jd: Clean up Boltzmann solvation (ion-position-dependent state)
    if (pub_is_pbe /= 'NONE') then
       call implicit_solvent_boltzmann_exit
    end if

    ! jd: Clean up smeared ions, if in use
    if(pub_is_smeared_ion_rep) call smeared_ion_exit

    ! jd: Force a recalculation of IS dielectric
    have_initial_eps = .false.

    do isub=1,mdl%nsub
       ! ndmh: Deallocate function basis objects
       if (pub_hubbard) call function_basis_deallocate(hub_proj_basis(isub))
       if (pub_any_nl_proj.or.pub_paw) &
            call function_basis_deallocate(proj_basis(isub))
       if (pub_paw.and.pub_eels_calculate) &
            call function_basis_deallocate(core_basis(isub))
       call function_basis_deallocate(ngwf_basis(isub))
       if (pub_use_aux_ngwfs) then
          call function_basis_deallocate(aux_ngwf_basis(isub))
       end if
       if (pub_cond_calculate) then
          call function_basis_deallocate(cond_ngwf_basis(isub))
          call function_basis_deallocate(joint_ngwf_basis(isub))
       end if
    end do

    ! ja531 --> Deallocate SW function basis in eigenstates
    if(pub_dos_smear > 0.0_dp .and. pub_pdos_max_l >= 0) then
       call eigenstates_cleanup_pdos()
    end if

    ! deallocate universal tightbox array
    call function_basis_exit_uni_tb(mdl%uni_tightbox)

    ! clean up parallel strategy
    do isub=1,mdl%nsub
       call parallel_strategy_exit(mdl%regions(isub)%par)
    end do

    ! ndmh: deallocate box FFT arrays
    if (pub_aug) call fourier_exit_box(mdl%aug_box)
    if (pub_tightbox_fft_coarse) call fourier_exit_box(mdl%uni_tightbox)
    call fourier_exit_fftbox

    ! ndmh: deallocate augmentation box
    call augmentation_box_exit(mdl%aug_box)

    ! Deallocate projector set
    if (pub_paw.or.pub_any_nl_proj) then
       do isub=1,mdl%nsub
          call projectors_deallocate_set(nl_projectors(isub))
       enddo
    end if

    ! jd: If this was a Boltzmann in-solvent calculation but we're not in
    !     auto-solvation mode, let the user know what the pure electrolyte
    !     energy was. Otherwise they won't be able to calculate dG_solvation.
    if(pub_is_implicit_solvent .and. .not. pub_is_auto_solvation .and. &
         pub_is_pbe /= "NONE") then
       if(pub_on_root) then
          write(stdout,'(/a)') '                                         &
               &             hartree           kcal/mol'
          write(stdout,'(a,f22.14,a,f16.6,a)') &
               'IS: NB: Energy of pure electrolyte:    ', &
               implicit_solvent_energy_terms%E_pure, '   ', &
               implicit_solvent_energy_terms%E_pure * HARTREE_TO_KCAL_PER_MOL
          if(pub_print_qc) then
             call utils_qc_print('E_pure_electrolyte', &
                  implicit_solvent_energy_terms%E_pure * HARTREE_TO_KCAL_PER_MOL)
          end if
       end if
    end if

    ! ndmh: if return value indicating convergence was requested, set it now
    if (present(return_converged)) return_converged = ngwfs_converged

    ! ab: print energy and force components in dftb
    if (pub_dftb) then
       call dftb_print_energy_components(dftb_eht_energy, dftb_rep_energy, &
            dftb_srb_energy, dftb_ies_energy, dftb_disp_energy, total_energy)
       if ((pub_write_forces .or. pub_task == 'GEOMETRYOPTIMIZATION' .or. &
             pub_task == 'TRANSITIONSTATESEARCH' .or. pub_task == 'MOLECULARDYNAMICS' .or. &
             pub_task == 'PHONON').and..not.out_of_runtime) then
          call dftb_print_force(F_dftb_eht,'DFTB Electronic Forces','Unconstrained',mdl)
          call dftb_print_force(F_dftb_rep,'DFTB Repulsion Forces','Unconstrained',mdl)
          call dftb_print_force(F_dftb_srb,'DFTB SRB Forces','Unconstrained',mdl)
          call dftb_print_force(F_dftb_ies,'DFTB Electrostatic Forces','Unconstrained',mdl)
          call dftb_print_force(F_dftb_dis,'Dispersion Forces','Unconstrained',mdl)
          call forces_calculate_correction(total_forces, mdl)
       end if
       ! jd: Clean up DFTB
       call dftb_free(mdl)
    end if

    call timer_clock(myself,2)

    return

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  contains

!-------------------------------------------------------------------------------

  subroutine internal_kpoints_init(mdl, kpt_cart, kpt)

    !==========================================================================!
    ! k-points initialisations and checks for energy_and_force_calculate       !
    !                                                                          !
    ! Created by JM Escartin, from existing code in energy_and_force_calculate !
    ! August 2016.                                                             !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout
    use model_type, only: MODEL
    use rundat, only: pub_num_kpoints, &
         pub_kpoint_method, PUB_1K
    use utils, only: utils_assert
    implicit none

    ! Arguments
    type(MODEL), intent(in)  :: mdl
    real(kind=DP),intent(out) :: kpt_cart(3) ! single k-points, in Cartesian coordinates
    real(kind=DP), intent(in), optional :: kpt(3) ! optional input

    ! Local variables
    integer :: ik

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, myself//&
         ': not ready yet for more than one k-point.')

    ! agrecokpt: currently only KP method is being implemented
    ! agrecokpt: now TB method should work for the basic case,
    ! I added warnings in the individual components if a given
    ! term or scheme has not been checked for TB method
    !call utils_assert(pub_kpoint_method == 'KP', myself//&
    !     ': currently supports only KP method for BZ sampling')

    ! agrecokpt: default case is Gamma point
    if (present(kpt)) then
       kpt_cart(:) = kpt(:)
    else
       kpt_cart = 0.0_DP
    end if

    if (any(kpt_cart .ne. 0.0_DP)) then
       if (pub_on_root) then
          ! agrecokpt: print the method it will be used
          write(stdout,'(a,a,a)') 'Non Gamma point calculation with ', &
               trim(pub_kpoint_method), ' method will be performed.'
          ! agrecokpt: test conversion from fractional to cartesian
          ! coordinates
          write(stdout,'(a,3f16.5)') 'Cartesian coordinates of k-point are ', &
               kpt_cart(:)
       end if
    end if

  end subroutine internal_kpoints_init

!-------------------------------------------------------------------------------

    subroutine internal_create_aux_rep

      use datatypes, only: data_functions_alloc
      use hamiltonian, only: hamiltonian_dens_indep_matrices, & ! jhl52
           hamiltonian_dens_dep_nonsc ! jhl52
      use kernel, only: kernel_create, kernel_destroy
      use sparse_embed, only: sparse_embed_product, sparse_embed_transpose_structure, &
           sparse_embed_transpose, SPAM3_EMBED, sparse_embed_create, &
           sparse_embed_destroy, SPAM3_EMBED_ARRAY, sparse_embed_array_create, &
           sparse_embed_array_product, sparse_embed_array_destroy

      ! Local Variables
      type(SPAM3_EMBED) :: Tt, SiTt, TSi
      type(SPAM3_EMBED_ARRAY) :: KTSi

      ! Set occ
      aux_rep%n_occ(UP,PUB_1K) = rep%n_occ(UP,PUB_1K)
      if (pub_num_spins > 1) then
         aux_rep%n_occ(DN,PUB_1K) = rep%n_occ(DN,PUB_1K)
      end if

      ! Allocate SPAM3_EMBED structures for overlap matrix, kinetic energy,
      ! density kernel, inverse overlap
      call ngwf_rep_create(aux_rep,'a',mdl, is_cmplx=loc_cmplx)
      ! agrecocmplx
      call kernel_create(aux_denskern, 'Ka', is_cmplx=loc_cmplx)

      ! Allocate storage for NGWFs
      do isub=1,mdl%nsub
         call data_functions_alloc(aux_rep%ngwfs_on_grid(isub), &
              aux_ngwf_basis(isub)%size_on_grid, iscmplx=loc_cmplx)
      end do

      ! Create the aux_ham
      call ngwf_ham_create(aux_ham,aux_rep)

      ! Now actually initialise the NGWF values (either from PAOs, STO-3G,
      ! or by loading from disk (or memory)
      call ngwfs_initialise(aux_rep,aux_ngwf_basis,mdl,hfxstate)

      ! Build up the density-indep part Hamiltonian and aux_rep%overlap,
      ! cross_overlap, and inv_overlap
      ! agrecokpt: at specified k-point?
      call hamiltonian_dens_indep_matrices(aux_rep, aux_ngwf_basis, &
           proj_basis, nl_projectors, hub_proj_basis, hub, mdl, rep, ngwf_basis, &
           kpt=kpt_cart)

      ! Construct the aux_denskern
      ! agrecocmplx
      call sparse_embed_transpose_structure(Tt%structure,aux_rep%cross_overlap)
      call sparse_embed_create(Tt,iscmplx=loc_cmplx)
      call sparse_embed_transpose(Tt,aux_rep%cross_overlap)
      call sparse_embed_create(SiTt,aux_rep%inv_overlap,Tt)
      call sparse_embed_product(SiTt,aux_rep%inv_overlap,Tt)

      call sparse_embed_create(TSi,aux_rep%cross_overlap,aux_rep%inv_overlap)
      call sparse_embed_product(TSi,aux_rep%cross_overlap,aux_rep%inv_overlap)
      call sparse_embed_array_create(KTSi, denskern%kern, TSi)

      call sparse_embed_array_product(KTSi, denskern%kern, TSi)
      call sparse_embed_array_product(aux_denskern%kern, SiTt, KTSi)

      ! Build up aux_ham
      call hamiltonian_dens_dep_nonsc(aux_ham, aux_rep, aux_ngwf_basis, &
           mdl,hfxstate,lhxc_fine,hub,rep,ham,denskern%kern%m(:,PUB_1K))

      call sparse_embed_array_destroy(KTSi)
      call sparse_embed_destroy(TSi)
      call sparse_embed_destroy(SiTt)
      call sparse_embed_destroy(Tt)

    end subroutine internal_create_aux_rep


    subroutine internal_destroy_aux_rep
      use datatypes, only: data_functions_dealloc
      implicit none

      ! ndmh: deallocate storage for NGWFs
      do isub=1,mdl%nsub
         call data_functions_dealloc(aux_rep%ngwfs_on_grid(isub))
      end do

      ! ndmh: deallocate storage for Hamiltonian, and whole-cell arrays
      call ngwf_ham_destroy(aux_ham)

      ! ndmh: memory deallocation for sparse matrices
      call kernel_destroy(aux_denskern)

      ! ndmh: deallocate storage for overlaps etc
      call ngwf_rep_destroy(aux_rep)

    end subroutine internal_destroy_aux_rep

!-------------------------------------------------------------------------------

    !==========================================================================!
    ! Allocate the number of procs to be included for each parallel strategy.  !
    !--------------------------------------------------------------------------!
    ! Written by Robert Charlton, 11/05/2018.                                  !
    !==========================================================================!

    subroutine internal_initialise_par

      use rundat, only: pub_parallel_scheme

      implicit none

      ! local variables
      integer :: procs_per_par(mdl%nsub)
      integer :: proc_to_give(1), ii
      real(kind=DP) :: remainder(mdl%nsub)

      do isub=1,mdl%nsub
         ! rc2013: set the region index
         mdl%regions(isub)%par%par_index = isub
      end do

      ! rc2013: decide how to distribute the procs
      select case(pub_parallel_scheme)
      ! Spread atoms across all processors regardless of region
      case('NONE')
         procs_per_par = pub_total_num_procs

         if(mdl%nsub .gt. 1 .and. pub_on_root) write(stdout,'(a/a)') &
              myself//': WARNING in internal_initialise_par. &
              &pub_parallel_scheme == NONE, which is not an efficient &
              &way of distributing resources; consider using HOUSE instead.'

      ! Split processore equally between regions
      case('SENATE')
         ! Make sure we have enough resources for this
         call utils_assert(pub_total_num_procs .ge. mdl%nsub, myself//&
              ' : pub_parallel_scheme: '//CRLF// &
              '   SENATE requested, but the number of regions in '//CRLF// &
              '   %block species_ngwf_regions exceeds the number of '//CRLF// &
              '   processors. Please either allocate additional '//CRLF// &
              '   processors,  or set pub_parallel_scheme: NONE in input file.')
         procs_per_par = pub_total_num_procs/mdl%nsub

      ! Default: Distribute processors proportionally between regions
      case('HOUSE')
         ! Make sure we have enough resources for this
         call utils_assert(pub_total_num_procs .ge. mdl%nsub, myself//&
              ' : pub_parallel_scheme: '//CRLF// &
              '   HOUSE requested, but the number of regions in '//CRLF// &
              '   %block species_ngwf_regions exceeds the number of '//CRLF// &
              '   processors. Please either allocate additional '//CRLF// &
              '   processors,  or set pub_parallel_scheme: NONE in input file.')
         procs_per_par = 1

         ! Calculate the number of procs per par with Huntingdon-Hill algorithm
         ! Each region already has 1 proc, so subtract that from total
         do ii=1,pub_total_num_procs-mdl%nsub
            do isub=1,mdl%nsub
               remainder(isub) = mdl%regions(isub)%par%nat/ &
                    sqrt(real(procs_per_par(isub)*(procs_per_par(isub)+1)))
            end do
            ! Give the next proc to the region with the highest remainder
            proc_to_give = maxloc(remainder)
            procs_per_par(proc_to_give(1)) = procs_per_par(proc_to_give(1)) + 1
         end do
      end select

      mdl%regions(1)%par%first_proc = 0
      mdl%regions(1)%par%num_procs = procs_per_par(1)
      do isub=2,mdl%nsub
         if(pub_parallel_scheme == 'NONE') then
            mdl%regions(isub)%par%first_proc = 0
         else
            mdl%regions(isub)%par%first_proc = &
                 mdl%regions(isub-1)%par%first_proc + procs_per_par(isub-1)
         end if
         mdl%regions(isub)%par%num_procs = procs_per_par(isub)
      end do

      if(mdl%nsub == 1) then
         call utils_assert(mdl%regions(1)%par%first_proc == 0, &
              'Error in internal_initialise_par: failure to allocate &
              &parallel strategies.')
         call utils_assert(procs_per_par(1) == pub_total_num_procs, &
              'Error in internal_initialise_par: there is only one parallel &
              &strategy, but the total number of procs allocated to this par &
              &does not match the number of procs available.')
      endif

    end subroutine internal_initialise_par

  end subroutine energy_and_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module energy_and_force
