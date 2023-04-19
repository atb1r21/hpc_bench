! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:

module rundat_blocks

  !=====================================================================!
  ! This module contains the subroutines that read the simulation cell  !
  ! and species parameters from the input file in form of data blocks.  !
  !                                                                     !
  ! The routines in this module are adapted from the former routines    !
  ! services_read_system_input and services_read_system_input_old.      !
  !---------------------------------------------------------------------!
  ! This module was created by Alvaro Ruiz Serrano in November 2010.    !
  ! based on code by Chris-Kriton Skylaris, Nicholas Hine, Peter Haynes !
  ! and Arash Mostofi.                                                  !
  ! Subsequent modifications by Robert Bell, Jacek Dziedzic, Gilberto   !
  ! Teobaldi, Quintin Hill, Louis Lee and Robert Charlton.              !
  !=====================================================================!

  use constants, only: DP
  use rundat, only: pub_debug_on_root, pub_debug

  implicit none

  public :: rundat_blocks_exec
  public :: rundat_blocks_exit

contains



  subroutine rundat_blocks_exec(mdl,tag)

    !==============================================================!
    ! Driver routine that reads esdf_block from the input file     !
    ! using new input file syntax and initialises various          !
    ! parameters of the simulation cell and the ppd grid.          !
    !--------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano on 17/11/2010.                !
    !==============================================================!

    use dftb, only: dftb_init_stage_1
    use geometry, only: POINT
    use model_type, only: MODEL
    use rundat, only: pub_is_implicit_solvent, pub_ngwf_regions_ngroups, &
         pub_dftb
    use utils, only: utils_assert, utils_dealloc_check

    implicit none

    ! ars: arguments
    type(MODEL),                intent(inout) :: mdl
    character(len=*), optional, intent(in   ) :: tag

    ! ars: local variables
    type(POINT) :: a1,a2,a3
    type(POINT) :: a1_pad,a2_pad,a3_pad
    integer, allocatable      :: nat_reg(:), nsp_reg(:)
    integer     :: isub, ierr

    ! ars: read blocks from input file
    call rundat_blocks_read(a1, a2, a3, a1_pad, a2_pad, a3_pad, &
         mdl%nat, mdl%nat_classical, mdl%nsp, &
         mdl%elements, mdl%classical_elements, mdl%species, tag, &
         mdl%regions, nat_reg, nsp_reg)

    ! rc2013: multiple regions not compatible with TSSearch
    if (present(tag)) then
       mdl%nsub = 1
       mdl%regions(1)%par%nat = mdl%nat
       mdl%regions(1)%par%num_species = mdl%nsp
       mdl%regions(1)%reg_count = 1
       mdl%regions(1)%par%nat_classical = mdl%nat_classical
    else
       mdl%nsub                = pub_ngwf_regions_ngroups
       ! rc2013: assign subsystem structures
       do isub=1,mdl%nsub
          mdl%regions(isub)%par%nat           = nat_reg(isub)
          mdl%regions(isub)%par%num_species   = nsp_reg(isub)
          mdl%regions(isub)%reg_count         = isub
          mdl%regions(isub)%par%nat_classical = 0
       enddo
    end if

    if(pub_debug) then
       ! ars: check data if #DEBUG
       call internal_debug()
    end if

    ! ars: initialise grid
    call rundat_blocks_init_grid(a1, a2, a3, a1_pad, a2_pad, a3_pad, mdl)

    ! ars: check compatibility
    call utils_assert(.not. ((mdl%nat_classical .gt. 0) &
         .and. pub_is_implicit_solvent),'is_implicit_solvent cannot be used &
         &simultaneously with classical atoms, sorry.')

    ! rc2013: deallocate embedding arrays
    if(allocated(nat_reg)) then
       deallocate(nat_reg, stat=ierr)
       call utils_dealloc_check('internal_species_ngwf_regions','nat_reg',ierr)
    end if
    if(allocated(nsp_reg)) then
       deallocate(nsp_reg, stat=ierr)
       call utils_dealloc_check('internal_species_ngwf_regions','nsp_reg',ierr)
    end if

    ! jd: Stage-1 of DFTB initialisation -- we want to adjust #NGWFs here
    if(pub_dftb) then
       call dftb_init_stage_1(mdl)
    end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains

      subroutine internal_debug()

      use comms, only: comms_barrier, pub_on_root
      use constants, only: stdout
      use hubbard_init, only: h_species
      use rundat, only: pub_hubbard, pub_ngwf_regions_ngroups

      implicit none

      integer :: ii, row, ireg

      call comms_barrier

      if (pub_on_root) then
         write(stdout,'(/a,i6)') 'DEBUG: No of Atoms:',mdl%nat
         write(stdout,'(a,i5)') 'DEBUG: No of Species: ',mdl%nsp
         if (pub_hubbard) then
            ! jcap: Changed this to explicitly reference
            ! regions(1)%par, as mdl%par hasn't been set to point to
            ! mdl%regions(1)%par yet
            write(stdout,'(/a,i6)') 'DEBUG: No of Hubbard Atoms:',mdl%regions(1)%par%nat_hub
            write(stdout,'(a,i5)') 'DEBUG: No of Hubbard Species: ',mdl%regions(1)%par%num_hub_species
         endif
         write(stdout,*) 'DEBUG: =======================================&
              &============================'
         do ii=1,mdl%nat
            write(stdout,*) 'DEBUG: Atom       :',ii
            write(stdout,*) 'DEBUG: Symbol     : ',mdl%elements(ii)%symbol
            write(stdout,*) 'DEBUG: Z          :',mdl%elements(ii)%atomic_number
            write(stdout,*) 'DEBUG: NGWFs      :',mdl%elements(ii)%nfunctions
            write(stdout,*) 'DEBUG: Radius     :',mdl%elements(ii)%radius
            write(stdout,*) 'DEBUG: Pseudopot  : ',&
                 trim(mdl%species(mdl%elements(ii)%species_number)%pseudo_name)
            write(stdout,*) 'DEBUG: Fireball   : ',&
                 trim(mdl%species(mdl%elements(ii)%species_number)%ngwf_set)
            write(stdout,'(a,3f10.5)') ' DEBUG: Centre     :', &
                 mdl%elements(ii)%centre%x,&
                 mdl%elements(ii)%centre%y,mdl%elements(ii)%centre%z
            write(stdout,*) 'DEBUG: Constraint : ', &
                 mdl%elements(ii)%ion_constraint_type
            write(stdout,'(a,3f10.5)') ' DEBUG: Constraint :', &
                 mdl%elements(ii)%ion_constraint(1),&
                 mdl%elements(ii)%ion_constraint(2),&
                 mdl%elements(ii)%ion_constraint(3)
            if (pub_hubbard) then
               ! jcap: As above, explicitly reference
               ! mdl%regions(1)%par rather than mdl%par
               do row=1,mdl%regions(1)%par%num_hub_species
                  ! rc2013: EMBED_FIX!
                  !if (regions(ireg)%elements(ii)%species_id == h_species(row)%hub_species) then
                     write(stdout,*) 'DEBUG: Hubbard projector angular momentum :',&
                          h_species(row)%hub_ang_mom
                     write(stdout,*) 'DEBUG: Hubbard U parameter (Ha)           :',&
                          h_species(row)%hub_u
                     write(stdout,*) 'DEBUG: Hubbard projector effective charge :',&
                          h_species(row)%hub_charge
                     write(stdout,*) 'DEBUG: Hubbard alpha parameter (Ha)       :',&
                          h_species(row)%hub_alpha
                     write(stdout,*) 'DEBUG: Hubbard spin splitting (Ha)        :',&
                          h_species(row)%hub_spin_splitting
                  !endif
               enddo
            endif
            if (any(mdl%elements(:)%radius_cond>0.0_DP)) then
               write(stdout,*) 'DEBUG: Conduction NGWFs   :',mdl%elements(ii)%nfunctions_cond
               write(stdout,*) 'DEBUG: Conduction radius  :',mdl%elements(ii)%radius_cond
            end if
            write(stdout,*) 'DEBUG: ====================================&
                 &==============================='
         end do
      end if

      ! rc2013: subsystem info if there's more than 1 system specified
      if (pub_on_root .and. (pub_ngwf_regions_ngroups .gt. 1)) then
         write(stdout,'(/a,i6)') 'DEBUG: number of regions: ',pub_ngwf_regions_ngroups
         do ireg=1,pub_ngwf_regions_ngroups
            write(stdout,'(/a,i6)') 'DEBUG: No of Atoms:  ',nat_reg(ireg)
            write(stdout,'(a,i5)') 'DEBUG: No of Species: ',nsp_reg(ireg)
            if (pub_hubbard) then
               ! jcap: As above, explicitly reference
               ! mdl%regions(1)%par rather than mdl%par
                write(stdout,'(/a,i6)') 'DEBUG: No of Hubbard Atoms:',mdl%regions(1)%par%nat_hub
                write(stdout,'(a,i5)') 'DEBUG: No of Hubbard Species: ',mdl%regions(1)%par%num_hub_species
            endif
            write(stdout,*) 'DEBUG: =======================================&
                &============================'
            do ii=1,nat_reg(ireg)
               write(stdout,*) 'DEBUG: Atom       : ',ii
               write(stdout,*) 'DEBUG: Symbol     : ',mdl%regions(ireg)%elements(ii)%symbol
               write(stdout,*) 'DEBUG: Z          : ',mdl%regions(ireg)%elements(ii)%atomic_number
               write(stdout,*) 'DEBUG: NGWFs      : ',mdl%regions(ireg)%elements(ii)%nfunctions
               write(stdout,*) 'DEBUG: Radius     : ',mdl%regions(ireg)%elements(ii)%radius
               write(stdout,*) 'DEBUG: Pseudopot  : ',&
                   trim(mdl%regions(ireg)%species(mdl%regions(ireg)%elements(ii)%species_number)%pseudo_name)
               write(stdout,*) 'DEBUG: Fireball   : ',&
                   trim(mdl%regions(ireg)%species(mdl%regions(ireg)%elements(ii)%species_number)%ngwf_set)
               write(stdout,'(a,3f10.5)') ' DEBUG: Centre     : ', &
                    mdl%regions(ireg)%elements(ii)%centre%x,&
                    mdl%regions(ireg)%elements(ii)%centre%y,&
                    mdl%regions(ireg)%elements(ii)%centre%z
               write(stdout,*) 'DEBUG: Constraint : ', &
                   mdl%regions(ireg)%elements(ii)%ion_constraint_type
               write(stdout,'(a,3f10.5)') ' DEBUG: Constraint : ', &
                   mdl%regions(ireg)%elements(ii)%ion_constraint(1),&
                   mdl%regions(ireg)%elements(ii)%ion_constraint(2),&
                   mdl%regions(ireg)%elements(ii)%ion_constraint(3)
               if (pub_hubbard) then
                  ! jcap: As above, explicitly reference
                  ! mdl%regions(1)%par rather than mdl%par
                  do row=1,mdl%regions(1)%par%num_hub_species
                     if (mdl%regions(ireg)%elements(ii)%species_id == h_species(row)%hub_species) then
                        write(stdout,*) 'DEBUG: Hubbard projector angular momentum :',&
                             h_species(row)%hub_ang_mom
                        write(stdout,*) 'DEBUG: Hubbard U parameter (Ha)           :',&
                             h_species(row)%hub_u
                        write(stdout,*) 'DEBUG: Hubbard projector effective charge :',&
                             h_species(row)%hub_charge
                        write(stdout,*) 'DEBUG: Hubbard alpha parameter (Ha)       :',&
                             h_species(row)%hub_alpha
                        write(stdout,*) 'DEBUG: Hubbard spin splitting (Ha)        :',&
                             h_species(row)%hub_spin_splitting
                     endif
                  enddo
               endif
               if (any(mdl%regions(ireg)%elements(:)%radius_cond>0.0_DP)) then
               write(stdout,*) 'DEBUG: Conduction NGWFs   :',mdl%regions(ireg)%elements(ii)%nfunctions_cond
               write(stdout,*) 'DEBUG: Conduction radius  :',mdl%regions(ireg)%elements(ii)%radius_cond
               end if
               write(stdout,*) 'DEBUG: ====================================&
                   &==============================='
            end do
         end do
      end if

    end subroutine internal_debug

  end subroutine rundat_blocks_exec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine rundat_blocks_read(a1, a2, a3, a1_pad, a2_pad, a3_pad, &
       nat, class_nat, nsp, elements, classical_elements, species, tag, &
       regions, nat_reg, nsp_reg)

    !=============================================================!
    ! This subroutine reads the input file using the esdf module  !
    !-------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2001.                   !
    ! Modified by Peter D. Haynes in 2004 so that the root proc   !
    ! reads and communicates across all procs.                    !
    ! Modified by Chris-Kriton Skylaris on 09/11/2004 so that it  !
    ! determines grids by calling ppd_strategy rather than        !
    ! reading them from input.                                    !
    !=============================================================!
    ! Rewritten by Arash A Mostofi, July 2005                     !
    ! NGWF plot block added by Chris-Kriton Skylaris on 13/03/2006!
    ! Hubbard DFT+U block added by D.D. O'Regan on 15/04/2009     !
    !=============================================================!
    ! Adapted and cleaned up from services_read_system_input by   !
    ! Alvaro Ruiz Serrano, 16/11/2010.                            !
    !=============================================================!


    use comms, only: pub_on_root
    use constants, only: stdout
    use geometry, only: POINT
    use ion, only: ELEMENT, SPECIE
    use model_type, only: REGION
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_ngwf_regions_ngroups
    use sw_resolution_of_identity, only: swri_library_size
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! ars: arguments
    type(POINT),                intent(  out) :: a1,a2,a3
    type(POINT),                intent(  out) :: a1_pad,a2_pad,a3_pad
    integer,                    intent(  out) :: nat
    integer,                    intent(  out) :: class_nat
    integer,                    intent(  out) :: nsp
    type(REGION), allocatable, target, intent(inout) :: regions(:)
    type(ELEMENT), allocatable, intent(inout) :: elements(:)
    type(ELEMENT), allocatable, intent(inout) :: classical_elements(:)
    type(SPECIE), allocatable,  intent(inout) :: species(:)
    ! rc2013: number of atoms in each "block"
    integer, allocatable,       intent(  out) :: nat_reg(:)
    integer, allocatable,       intent(  out) :: nsp_reg(:)
    character(len=*), optional, intent(in   ) :: tag

    ! ars: ---- local variables ----
    integer :: hub_nsp, cdft_nsp, ierr
    integer :: num_hub_proj
    character(len=4), pointer, dimension(:) :: id_list
    character(len=4), allocatable :: all_species_id(:)
    integer :: swri_entry
    type(PARAL_INFO), pointer     :: par

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering rundat_blocks_read'


    ! ars: ============= READ LATTICE BLOCKS =============
    call internal_lattice()
    call internal_padded_lattice()
    ! ars: =========== END READ LATTICE BLOCKS ===========




    ! ars: ============ READ ATOMIC POSITIONS ============
    call internal_positions()

    if (present(tag)) then
       return
    end if

    ! ars: allocate global arrays
    if (pub_on_root) then
       allocate(id_list(1:nsp),stat=ierr)
       call utils_alloc_check('rundat_blocks_read','=>id_list',ierr)
    end if
    allocate(all_species_id(1:nsp), stat=ierr)
    call utils_alloc_check('rundat_blocks_read','all_species_id',ierr)
    ! ars: ========== END READ ATOMIC POSITIONS ==========

    ! rc2013: ===== READ SPECIES - NGWF REGIONS ====
    call internal_species_ngwf_regions()



    ! smmd: ============ READ ATOMIC VELOCITIES ==========
    call internal_velocities()
    ! smmd: ========== END READ ATOMIC VELOCITIES ========

    ! smmd: ======== READ THERMOSTAT DEFINITION ==========
    call internal_thermostat()
    ! smmd: ====== END READ THERMOSTAT DEFINITION ========


    ! ars: ======= READ SPECIES - MANDATORY BLOCKS =======
    call internal_species()
    call internal_species_pot()
    ! ars: ===== END READ SPECIES - MANDATORY BLOCKS =====

    ! jd: ============ READ SWRI BLOCK ==========
    !     This needs to come before species_swri.
    call internal_swri()
    ! jd: ========== END READ SWRI BLOCK ========

    ! ars: ======= READ SPECIES - OPTIONAL BLOCKS ========
    call internal_species_cond()
    call internal_species_aux()
    call internal_species_core_wf()

    species(:)%ngwf_set = "UNSET"
    species(:)%cond_ngwf_set = "UNSET"
    species(:)%aux_ngwf_set = "UNSET"

    call internal_species_atomic_set('')
    call internal_species_atomic_set('_cond')
    call internal_species_atomic_set('_aux')
    call internal_species_constraints()
    call internal_species_ngwf_plot()
    call internal_species_ldos_groups()
    call internal_species_bsunfld_groups()
    call internal_species_bsunfld_projatoms()
    call internal_species_pdos_groups()
    call internal_species_locdipole_groups()
    call internal_species_scissor()
    call internal_species_solvent_radius()
    do swri_entry = 1, swri_library_size
       call internal_species_swri(swri_entry)
    end do

    ! tjz07: ====== READ SPECIES - OPTIONAL BLOCK TDDFT ==
    call internal_species_tddft_kernel()
    call internal_species_tddft_ct()
    ! these define the sparsity pattern of the TDDFT kernel

    ! lpl: NBO block
    call internal_species_write_nbo()

    ! ars: ===== END READ SPECIES - OPTIONAL BLOCKS ======

    ! smmd: =========== READ TRANSPORT SETUP =============
    call internal_etrans_setup()
    call internal_etrans_leads()
    call internal_etrans_plot_eigchan_energies()
    ! smmd: ========= END READ TRANSPORT SETUP ===========

    ! lpl: =========== READ DDEC RCOMP BLOCK =============
    call internal_ddec_read_rcomp()
    ! lpl: This block is related to the still-dodgy 'ddec_rmse' routine
    !call internal_ddec_rmse_vdW()
    ! lpl: ========= END READ DDEC RCOMP BLOCK ===========

    ! ars: ============= READ HUBBARD BLOCK ==============
    call internal_hubbard()
    ! ars: =========== END READ HUBBARD BLOCK ============

    ! gibo: ======== READ CONSTRAINED_DFT BLOCK ======
    call internal_cdft()
    ! gibo:  ======== READ CONSTRAINED_DFT BLOCK ======

    ! gom: ======== READ DFT_NU_BLOCK ======
    call internal_dft_nu()
    ! gom:  ======== READ DFT_NU BLOCK ======


    ! ars: ======== READ CLASSICAL CHARGES BLOCK =========
    call internal_classical_atoms()
    ! ars: ====== END READ CLASSICAL CHARGES BLOCK =======

    ! jd: ============ READ SOL_IONS BLOCK ==========
    call internal_sol_ions()
    ! jd: ========== END READ SOL_IONS BLOCK ========

    ! jd: ============ READ IS+DIELECTRIC_EXCLUSIONS BLOCK ==========
    call internal_is_dielectric_exclusions()
    ! jd: ========== END READ IS_DIELECTRIC_EXCLUSIONS BLOCK ========

    ! jd: ============ READ THOLE_POLARISABILITIES BLOCK ==========
    call internal_thole_polarisabilities()
    ! jd: ========== END READ THOLE_POLARISABILITIES BLOCK ========

    ! jd: ============ READ MM_REP_PARAMS BLOCK ==========
    call internal_mm_rep_params()
    ! jd: ========== END READ MM_REP_PARAMS BLOCK ========

    ! gab: ============ READ IS_SOFT_SPHERE_RADII BLOCK =====================
    call internal_is_soft_sphere_radii()
    ! gab: ========== END READ IS_SOFT_SPHERE_RADII BLOCK ===================

    ! dhpt: ======== READ COUPLINGS BLOCK ======
    call internal_couplings()
    ! dhpt:  ======== READ COUPLINGS BLOCK ======

    ! mjsp: ============ READ EDA_IATM BLOCK ==========
    call internal_eda()
    ! mjsp: ========== END READ EDA_IATM BLOCK ========

    !qoh: VDW parameter override block
    call internal_vdwparam_override()
    ! rc2013: subsystem elements + species allocation
    call internal_subsystem_allocate()
    ! ars: deallocate global arrays
    if (pub_on_root) then
       deallocate(id_list,stat=ierr)
       call utils_dealloc_check('rundat_blocks_read','=>id_list',ierr)
    endif
    deallocate(all_species_id, stat=ierr)
    call utils_dealloc_check('rundat_blocks_read','all_species_id',ierr)

    if(pub_debug_on_root) then
       write(stdout,'(/a)') 'DEBUG: Leaving rundat_blocks_read'
    end if

    return

  contains


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !----------------------------    DATA BLOCKS    ----------------------------!
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!




    !-----------------------------
    !--------- %block lattice_cart and lattice_abc
    !-----------------------------

    subroutine internal_lattice()

      use comms, only: pub_root_proc_id, comms_bcast
      use constants, only: TWO_PI, PI
      use esdf, only: block_data, esdf_block, esdf_convfac
      use utils, only: utils_assert, utils_abort

      implicit none

      !-------------
      ! ars: dummies
      character(len=80) :: dummy_len_unit
      !-------------
      ! ars: buffers
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: real
      real(kind=dp) :: lenconfac
      ! ars: integer
      integer :: num_lines

      ! rab207: for lattice_abc
      logical :: lattice_defined
      integer, parameter :: LATTICE_NDEF=0, LATTICE_CART = 1, LATTICE_ABC = 2
      integer :: block_type
      real(kind=DP) :: a,b,c,alph,beta,gamm
      character(len=*), parameter :: myself = 'internal_lattice'

      ! -----------------------------------------------------------------------

      ! pdh: allocate buffers
      ! ddor: Altered for DFT+U
      allocate(dbuf(9),stat=ierr)
      call utils_alloc_check(myself,'dbuf',ierr)
      dbuf = 0.0_DP

      block_type = LATTICE_NDEF
      lattice_defined = .false.

      ! aam: Root proc only reads information
      if (pub_on_root) then

         ! <lattice_cart> block of input file
         if (esdf_block('lattice_cart',num_lines)) then

            lattice_defined = .true.
            block_type = LATTICE_CART

            ! aam: 3/6/09 added ability to read optional  unit string
            !      ("bohr" or "ang") in first line of lattice_cart block
            !      (default is bohr)
            lenconfac=0.0_dp ! jd: kills compiler warning
            if (num_lines == 4) then
               read(block_data(1),'(a)') dummy_len_unit
               read(block_data(2),*) dbuf(1:3)
               read(block_data(3),*) dbuf(4:6)
               read(block_data(4),*) dbuf(7:9)
               lenconfac = esdf_convfac(dummy_len_unit,'bohr')
            else if (num_lines == 3) then
               read(block_data(1),*) dbuf(1:3)
               read(block_data(2),*) dbuf(4:6)
               read(block_data(3),*) dbuf(7:9)
               lenconfac=1.0_dp
            else
               call utils_abort('Error in '//myself//': &
                    &malformed %block lattice_cart specification')
            end if
            dbuf(1:9)=lenconfac*dbuf(1:9)
         end if

         ! <lattice_abc> block of input file
         if (esdf_block('lattice_abc',num_lines)) then
            call utils_assert(.not. lattice_defined, &
                 'Error in '//myself//': Cannot define both %block &
                 &lattice_abc and %block lattice_cart.')

            lattice_defined = .true.
            block_type = LATTICE_ABC

            lenconfac=0.0_dp ! jd: kills compiler warning
            if (num_lines == 3) then
               read(block_data(1),'(a)') dummy_len_unit
               read(block_data(2),*) dbuf(1:3)
               read(block_data(3),*) dbuf(4:6)
               lenconfac = esdf_convfac(dummy_len_unit,'bohr')
            else if (num_lines == 2) then
               read(block_data(1),*) dbuf(1:3)
               read(block_data(2),*) dbuf(4:6)
               lenconfac=1.0_dp
            else
               call utils_abort('Error in '//myself//': &
                    &malformed %block lattice_abc specification.')
            end if
            ! convert lengths
            dbuf(1:3)=lenconfac*dbuf(1:3)
            ! convert angles to radians
            dbuf(4:6)=PI/180.0_DP*dbuf(4:6)
         endif
      end if

      ! pdh: broadcast and unpack
      call comms_bcast(pub_root_proc_id,dbuf,9)
      call comms_bcast(pub_root_proc_id,block_type)

      if (block_type == LATTICE_CART) then
         a1%x = dbuf(1) ; a1%y = dbuf(2) ; a1%z = dbuf(3)
         a2%x = dbuf(4) ; a2%y = dbuf(5) ; a2%z = dbuf(6)
         a3%x = dbuf(7) ; a3%y = dbuf(8) ; a3%z = dbuf(9)
      elseif (block_type == LATTICE_ABC) then
         ! rab207: use same convention as castep:
         !   1) a lies along x axis
         !   2) b lies in xy plane
         !   3) c chosen as positive square root
         a = dbuf(1) ; b = dbuf(2)   ; c = dbuf(3)
         alph=dbuf(4); beta = dbuf(5); gamm = dbuf(6)

         ! check for errors
         call utils_assert(a > 0.0_DP, &
              & 'Error in '//myself//': requires a > 0')
         call utils_assert(b > 0.0_DP, &
              & 'Error in '//myself//': requires b > 0')
         call utils_assert(c > 0.0_DP, &
              & 'Error in '//myself//': requires c > 0')
         call utils_assert(alph > 0.0_DP .and. alph < PI, &
              & 'Error in '//myself//': requires 0 < alpha < 180')
         call utils_assert(beta > 0.0_DP .and. beta < PI, &
              & 'Error in '//myself//': requires 0 < beta < 180')
         call utils_assert(gamm > 0.0_DP .and. gamm < PI, &
              & 'Error in '//myself//': requires 0 < gamma < 180')
         call utils_assert(alph+beta+gamm < TWO_PI, &
              & 'Error in '//myself//': requires alpha+beta+gamma < 360')
         call utils_assert(abs(alph-beta) < gamm, &
              & 'Error in '//myself//': requires abs(alpha-beta) < gamma')
         call utils_assert(abs(beta-gamm) < alph, &
              & 'Error in '//myself//': requires abs(beta-gamma) < alpha')
         call utils_assert(abs(gamm-alph) < beta, &
              & 'Error in '//myself//': requires abs(gamma-alpha) < beta')

         a1%x = a ; a1%y = 0.0_DP; a1%z = 0.0_DP
         a2%x = b*cos(gamm); a2%y = b*sin(gamm); a2%z = 0.0_DP
         a3%x = c*cos(beta)
         a3%y = c*(cos(alph)-cos(beta)*cos(gamm))/sin(gamm)
         a3%z = sqrt(c**2 - a3%x**2 - a3%y**2)

      else
         call utils_abort('Error in '//myself//'): Neither &
              &%block lattice_cart nor %block lattice_abc specified.')
      endif

      ! pdh: deallocate buffers
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check(myself,'dbuf',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <lattice_???> block'

    end subroutine internal_lattice



    !------------------------------------
    !--------- %block padded_lattice_cart and padded_lattice_abc
    !------------------------------------

    subroutine internal_padded_lattice

      use comms, only: pub_root_proc_id, comms_bcast
      use constants, only: PI, TWO_PI
      use esdf, only: block_data, esdf_block, esdf_convfac
      use geometry, only: magnitude, operator(.DOT.)
      use rundat, only: pub_coulomb_cutoff
      use utils, only: utils_abort, utils_assert

      implicit none


      !-------------
      ! ars: dummies
      character(len=80) :: dummy_len_unit
      !-------------
      ! ars: buffers
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: real
      real(kind=dp) :: lenconfac
      ! ars: integer
      integer :: num_lines

      ! rab207: for lattice_abc
      logical :: lattice_defined
      integer, parameter :: LATTICE_NDEF=0, LATTICE_CART = 1, LATTICE_ABC = 2
      integer :: block_type
      real(kind=DP) :: a,b,c,alph,beta,gamm
      real(kind=DP) :: cell_consistency_tolerance = 1.0d-7
      logical :: cell_consistent
      character(len=*), parameter :: myself = 'internal_padded_lattice'

      ! -----------------------------------------------------------------------

      ! ndmh: Root proc only: read information from
      !      <padded_lattice_cart> block of input file
      if (pub_coulomb_cutoff) then

         ! pdh: allocate buffers
         allocate(dbuf(9),stat=ierr)
         call utils_alloc_check(myself,'dbuf',ierr)

         block_type = LATTICE_NDEF
         lattice_defined = .false.

         if (pub_on_root) then

            ! <lattice_cart> block of input file
            if (esdf_block('padded_lattice_cart',num_lines)) then

               lattice_defined = .true.
               block_type = LATTICE_CART

               ! aam: 3/6/09 added ability to read optional  unit string
               !      ("bohr" or "ang") in first line of lattice_cart block
               !      (default is bohr)

               lenconfac=0.0_dp ! jd: kills compiler warning
               if (num_lines == 4) then
                  read(block_data(1),'(a)') dummy_len_unit
                  read(block_data(2),*) dbuf(1:3)
                  read(block_data(3),*) dbuf(4:6)
                  read(block_data(4),*) dbuf(7:9)
                  lenconfac = esdf_convfac(dummy_len_unit,'bohr')
               else if (num_lines == 3) then
                  read(block_data(1),*) dbuf(1:3)
                  read(block_data(2),*) dbuf(4:6)
                  read(block_data(3),*) dbuf(7:9)
                  lenconfac=1.0_dp
               else
                  call utils_abort('Error in '//myself//': &
                       &malformed %block padded_lattice_cart specification')
               end if
               dbuf(1:9)=lenconfac*dbuf(1:9)
            end if

            ! <lattice_abc> block of input file
            if (esdf_block('padded_lattice_abc',num_lines)) then
               call utils_assert(.not. lattice_defined, &
                    'Error in '//myself//': Cannot define both %block &
                    &padded_lattice_abc and %padded_block lattice_cart.')

               lattice_defined = .true.
               block_type = LATTICE_ABC

               lenconfac=0.0_dp ! jd: kills compiler warning
               if (num_lines == 3) then
                  read(block_data(1),'(a)') dummy_len_unit
                  read(block_data(2),*) dbuf(1:3)
                  read(block_data(3),*) dbuf(4:6)
                  lenconfac = esdf_convfac(dummy_len_unit,'bohr')
               else if (num_lines == 2) then
                  read(block_data(1),*) dbuf(1:3)
                  read(block_data(2),*) dbuf(4:6)
                  lenconfac=1.0_dp
               else
               call utils_abort('Error in '//myself//': &
                    &malformed %block padded_lattice_abc specification.')
               end if
               ! convert lengths
               dbuf(1:3)=lenconfac*dbuf(1:3)
               ! convert angles to radians
               dbuf(4:6)=PI/180.0_DP*dbuf(4:6)
            endif

            if (.not.lattice_defined) then
               block_type = LATTICE_CART
               dbuf(1) = -1.0_DP; dbuf(2) =  0.0_DP; dbuf(3) =  0.0_DP
               dbuf(4) =  0.0_DP; dbuf(5) = -1.0_DP; dbuf(6) =  0.0_DP
               dbuf(7) =  0.0_DP; dbuf(8) =  0.0_DP; dbuf(9) = -1.0_DP
            endif
         end if

         ! ndmh: broadcast and unpack
         call comms_bcast(pub_root_proc_id,dbuf,9)
         call comms_bcast(pub_root_proc_id,block_type)

         if (block_type == LATTICE_CART) then
            a1_pad%x = dbuf(1) ; a1_pad%y = dbuf(2) ; a1_pad%z = dbuf(3)
            a2_pad%x = dbuf(4) ; a2_pad%y = dbuf(5) ; a2_pad%z = dbuf(6)
            a3_pad%x = dbuf(7) ; a3_pad%y = dbuf(8) ; a3_pad%z = dbuf(9)
         elseif (block_type == LATTICE_ABC) then
            ! rab207: use same convention as castep:
            !   1) a lies along x axis
            !   2) b lies in xy plane
            !   3) c chosen as positive square root
            a = dbuf(1) ; b = dbuf(2)   ; c = dbuf(3)
            alph=dbuf(4); beta = dbuf(5); gamm = dbuf(6)

            ! check for errors
            call utils_assert(a > 0.0_DP, &
                 & 'Error in '//myself//': requires a > 0')
            call utils_assert(b > 0.0_DP, &
                 & 'Error in '//myself//': requires b > 0')
            call utils_assert(c > 0.0_DP, &
                 & 'Error in '//myself//': requires c > 0')
            call utils_assert(alph > 0.0_DP .and. alph < PI, &
                 & 'Error in '//myself//': requires 0 < alpha < 180')
            call utils_assert(beta > 0.0_DP .and. beta < PI, &
                 & 'Error in '//myself//': requires 0 < beta < 180')
            call utils_assert(gamm > 0.0_DP .and. gamm < PI, &
                 & 'Error in '//myself//': requires 0 < gamma < 180')
            call utils_assert(alph+beta+gamm < TWO_PI, &
                 & 'Error in '//myself//': requires alpha+beta+gamma < 360')
            call utils_assert(abs(alph-beta) < gamm, &
                 & 'Error in '//myself//': requires abs(alpha-beta) < gamma')
            call utils_assert(abs(beta-gamm) < alph, &
                 & 'Error in '//myself//': requires abs(beta-gamma) < alpha')
            call utils_assert(abs(gamm-alph) < beta, &
                 & 'Error in '//myself//': requires abs(gamma-alpha) < beta')

            a1_pad%x = a ; a1_pad%y = 0.0_DP; a1_pad%z = 0.0_DP
            a2_pad%x = b*cos(gamm); a2_pad%y = b*sin(gamm); a2_pad%z = 0.0_DP
            a3_pad%x = c*cos(beta)
            a3_pad%y = c*(cos(alph)-cos(beta)*cos(gamm))/sin(gamm)
            a3_pad%z = sqrt(c**2 - a3_pad%x**2 - a3_pad%y**2)

            ! rab207: check that padded cell definition is consistent with
            !         simulation cell definition
            cell_consistent = &
                 & abs((a1.dot.a1_pad)/magnitude(a1)/magnitude(a1_pad)-1.0_DP) < cell_consistency_tolerance
            cell_consistent = cell_consistent .and. &
                 & (abs((a2.dot.a2_pad)/magnitude(a2)/magnitude(a2_pad)-1.0_DP) < cell_consistency_tolerance)
            cell_consistent = cell_consistent .and. &
                 & (abs((a3.dot.a3_pad)/magnitude(a3)/magnitude(a3_pad)-1.0_DP) < cell_consistency_tolerance)
            call utils_assert(cell_consistent,'Error in '//myself//': padded&
                 & cell vectors are not aligned with simulation cell vectors')

         else
            call utils_abort('Error in '//myself//'): Neither %block padded_lat&
                 &tice_cart nor %block padded_lattice_abc specified.')
         endif

         ! pdh: allocate buffers
         ! ddor: Altered for DFT+U
         deallocate(dbuf,stat=ierr)
         call utils_dealloc_check(myself,'dbuf',ierr)

         if (pub_debug_on_root) write(stdout,'(a)') &
              'DEBUG: Read <padded_lattice_cart> block'

      end if

    end subroutine internal_padded_lattice


    !------------------------------
    !--------- %block positions
    !------------------------------

    subroutine internal_positions()

      !=============================================================!
      ! Modified 09/13 by Robert Bell                               !
      ! This subroutine now looks automatically for positional data !
      ! from absolute and fractional data blocks. If new positions  !
      ! blocks are added, an esdf key must be made for **both** the !
      ! abs and frac blocks. The code can also check for an .xyz    !
      ! file, however this needs to be coded in, and a corresponding!
      ! tag to rundat_mod added.                                    !
      !=============================================================!

      use comms, only: pub_root_proc_id, comms_bcast, comms_barrier
      use constants, only: CRLF
      use esdf, only: block_data, esdf_block, &
           esdf_convfac, esdf_line_divide, esdf_reduce
      use geometry, only: point, operator(+), operator(*)
      use rundat, only: pub_input_xyz_file, pub_input_xyz_file_product, &
           pub_input_xyz_file_intermediate
      use services, only: services_read_xyz
      use utils, only: utils_abort, utils_nth_word_of

      implicit none

      !--------------
      ! ars: dummies
      integer :: dummy_nat
      character(len=80) :: dummy_len_unit
      !--------------
      ! ars: buffers
      integer, allocatable :: ibuf(:)
      real(kind=DP), allocatable :: posbuf(:,:)
      character(len=4),allocatable :: species_id(:)
      !--------------
      ! ars: integers
      integer           :: ishift, row, ii
      ! ars: real
      real(kind=dp)     :: lenconfac
      !--------------
      ! smmd: allow for specification of groups of atoms
      integer, allocatable :: group_buf(:)
      integer :: maxp
      character(len=80) :: cbuf(4)

      ! rab207: allow positions defined in _abs _frac and .xyz file
      integer :: positions_defined
      integer :: positions_type
      character(len=80) :: tag_loc
      character(len=80) :: block_tag
      character(len=80) :: loc_input_xyz_file
      integer, parameter :: ABS_BLOCK = 1, FRAC_BLOCK = 2, XYZ_FILE = 3
      logical :: log_dummy ! required to take the return value of esdf_block, but never used
      ! jd:
      logical :: has_four_columns
      character(len=*), parameter :: myself = 'internal_positions'

      ! -----------------------------------------------------------------------

      ! check whether tag is provided
      if (present(tag)) then
         tag_loc = '_'//trim(tag)
      else
         tag_loc = ''
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read positions_abs or positions_frac or positions_xyz
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ishift=0 ! qoh: initialise on all processors

      ! aam: Root proc only: read absolute cartesian positions from
      !      <positions_???> block of input file
      if (pub_on_root) then
         lenconfac=1.0_dp
         nsp = 0
         ! loop over all different ways of defining positions in input
         positions_defined = 0
         if (esdf_block('positions_abs'//trim(tag_loc),dummy_nat)) then
            positions_defined = positions_defined + 1
            positions_type = ABS_BLOCK
            block_tag = 'positions_abs'//trim(tag_loc)
         endif
         if (esdf_block('positions_frac'//trim(tag_loc),dummy_nat)) then
            positions_defined = positions_defined + 1
            positions_type = FRAC_BLOCK
            block_tag = 'positions_frac'//trim(tag_loc)
         endif
         ! check the correct .xyz file for this tag
         if ( (trim(tag_loc) == '' .and. trim(pub_input_xyz_file)/='') .or. &
              (trim(tag_loc) == '_PRODUCT' .and. trim(pub_input_xyz_file_product)/='') .or. &
              (trim(tag_loc) == '_INTERMEDIATE' .and. trim(pub_input_xyz_file_intermediate)/='')) then
            positions_defined = positions_defined + 1
            positions_type = XYZ_FILE
            ! find the correct .xyz file for this tag
            select case (tag_loc)
            case ('')
               loc_input_xyz_file = pub_input_xyz_file
            case ('_PRODUCT')
               loc_input_xyz_file = pub_input_xyz_file_product
            case ('_INTERMEDIATE')
               loc_input_xyz_file = pub_input_xyz_file_intermediate
            case default
               call utils_abort('Error in '//myself//': reading from &
                    &xyz file not implemented for positions tag "'//trim(tag_loc)//'"')
            end select
            call services_read_xyz(loc_input_xyz_file, dummy_nat)
         endif

         ! check positions are defined once and only once
         if (positions_defined > 1) then
            call utils_abort('Error in '//myself//': &
                 &multiple positions blocks found for %block POSITIONS'//trim(tag_loc))
         endif
         if (positions_defined == 0) then
            ! should not reach here as check already done in onetep.F90
            call utils_abort('Error in '//myself//': &
                 &no positions data found for %block POSITIONS'//trim(tag_loc))
         endif


         ! sort out inputs contained within the .dat file
         if (positions_type == ABS_BLOCK .or. positions_type == FRAC_BLOCK) then

            log_dummy = esdf_block(trim(block_tag),dummy_nat)

            ! aam: 3/6/09 added ability to read optional  unit string
            !      ("bohr" or "ang") in first line of positions_abs block
            !      (default is bohr)
            if (positions_type == ABS_BLOCK) then
               read(block_data(1),'(a)') dummy_len_unit
               dummy_len_unit=esdf_reduce(dummy_len_unit)
               if ((index(dummy_len_unit,"ang").ne.0) .or. &
                    (index(dummy_len_unit,"bohr").ne.0)) then
                  lenconfac = esdf_convfac(dummy_len_unit,'bohr')
                  ishift=1
                  dummy_nat = dummy_nat - 1
               else
                  ! jd: Cover corner case of something specified that does not
                  !     look like 'ang', 'bohr' or a position of an atom,
                  !     which we distinguish by having four values columns.
                  read(block_data(1),'(a)',iostat=ierr) dummy_len_unit
                  if(ierr /= 0) dummy_len_unit = ''
                  has_four_columns = &
                       trim(utils_nth_word_of(dummy_len_unit,4)) /= ''
                  if(ierr /= 0 .or. .not. has_four_columns) then
                     call utils_abort('Error in '//myself//': &
                          &Unrecognized unit "'//trim(dummy_len_unit)//'" in &
                          &%block positions_abs'//trim(tag_loc)//' or &
                          &block positions_frac'//trim(tag_loc)//'.'//&
                          CRLF//'Alternatively, a malformatted position of &
                          &first atom in this block.')
                  end if
               end if

            end if
         end if

         nat = dummy_nat
      end if

      ! broadcast number of atoms found to all procs
      call comms_bcast(pub_root_proc_id,nat)

      ! allocate elements array if it is not already allocated
      if (.not.allocated(elements)) then
         allocate(elements(nat),stat=ierr)
         call utils_alloc_check(myself,'elements',ierr)
      else if (size(elements)/=nat) then
         call utils_abort('Error in '//myself//' &
              &size of elements array does not match number of atoms')
      end if

      ! pdh: allocate buffers
      allocate(posbuf(3,nat),stat=ierr)
      call utils_alloc_check(myself,'posbuf',ierr)
      allocate(ibuf(nat*4),stat=ierr)
      call utils_alloc_check(myself,'ibuf',ierr)
      allocate(group_buf(nat),stat=ierr)
      call utils_alloc_check(myself,'group_buf',ierr)
      allocate(species_id(nat),stat=ierr)
      call utils_alloc_check(myself,'species_id',ierr)

      if (pub_on_root) then

         if (positions_type == ABS_BLOCK .or. positions_type == FRAC_BLOCK) then

            if (nat > 0) nsp = 1
            do row=1,nat
               ! smmd: modified in order to allow for groups definition
               maxp = 5
               call esdf_line_divide(maxp,cbuf,block_data(row+ishift))
               if (maxp == 4) then
                  read(block_data(row+ishift),*) species_id(row),posbuf(:,row)
                  group_buf(row) = 1
               elseif (maxp == 5) then
                  read(block_data(row+ishift),*) species_id(row),posbuf(:,row),group_buf(row)
               else
                  call utils_abort('Error in '//myself//&
                       &': malformed positions block for tag "'//&
                       trim(tag_loc)//'" or broken logic in esdf_line_divide.')
               endif

            end do
            posbuf(:,:) = lenconfac*posbuf(:,:)
         elseif (positions_type == XYZ_FILE) then
            ! find the correct .xyz file for this tag
            select case (tag_loc)
            case ('')
               loc_input_xyz_file = pub_input_xyz_file
            case ('_PRODUCT')
               loc_input_xyz_file = pub_input_xyz_file_product
            case ('_INTERMEDIATE')
               loc_input_xyz_file = pub_input_xyz_file_intermediate
            case default
               call utils_abort('Error in '//myself//': reading from &
                    &xyz file not implemented for positions tag "'//trim(tag_loc)//'"')
            end select

            dummy_nat = nat
            write(stdout,'(/2x,2a)') "Reading POSITIONS"//trim(tag_loc)// &
                 & " geometry from ", trim(loc_input_xyz_file)
            call services_read_xyz(loc_input_xyz_file,dummy_nat,posbuf,species_id)
            group_buf(:) = 1
         end if

         if (nat > 0) nsp = 1
         do row=1,nat
            do ii=1,4
               ibuf((row-1)*4+ii) = iachar(species_id(row)(ii:ii))
            enddo

            ! aam: count number of distinct species
            sp_count: do ii=1,row-1
               if (species_id(row) == species_id(ii)) then
                  exit sp_count
               else
                  if (ii == row-1) nsp = nsp + 1
               end if
            end do sp_count

         enddo

      end if

      call comms_barrier

      ! pdh: broadcast and unpack positions
      call comms_bcast(pub_root_proc_id,posbuf)
      call comms_bcast(pub_root_proc_id,group_buf,nat)
      call comms_bcast(pub_root_proc_id,positions_type,1)

      if (positions_type == ABS_BLOCK .or. positions_type == XYZ_FILE) then
         do row=1,nat
            elements(row)%centre%x = posbuf(1,row)
            elements(row)%centre%y = posbuf(2,row)
            elements(row)%centre%z = posbuf(3,row)
            elements(row)%group_id = group_buf(row)
         end do
      elseif (positions_type == FRAC_BLOCK) then
         do row=1,nat
            elements(row)%centre = posbuf(1,row)*a1+posbuf(2,row)*a2+posbuf(3,row)*a3
            elements(row)%group_id = group_buf(row)
         enddo
      endif


      ! pdh: broadcast and unpack species_id's
      call comms_bcast(pub_root_proc_id,ibuf)
      do row=1,nat
         do ii=1,4
            elements(row)%species_id(ii:ii) = achar(ibuf((row-1)*4+ii))
         end do
      end do

      ! pdh: deallocate buffers
      deallocate(posbuf,stat=ierr)
      call utils_dealloc_check(myself,'posbuf',ierr)
      deallocate(ibuf,stat=ierr)
      call utils_dealloc_check(myself,'ibuf',ierr)
      deallocate(group_buf,stat=ierr)
      call utils_dealloc_check(myself,'group_buf',ierr)
      deallocate(species_id,stat=ierr)
      call utils_dealloc_check(myself,'species_id',ierr)

      ! pdh: broadcast number of species
      call comms_bcast(pub_root_proc_id,nsp)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <positions'//trim(tag_loc)//'> block'

    end subroutine internal_positions


    !------------------------------
    !--------- %block velocities
    !------------------------------

    subroutine internal_velocities()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_convfac
      use rundat, only: md_init_velocities
      use utils, only: utils_abort

      implicit none

      !--------------
      ! smmd: dummies
      integer :: dummy_nat
      character(len=80) :: dummy_vel_unit
      !--------------
      ! smmd: buffers
      real(kind=DP), allocatable :: dbuf(:)
      !--------------
      ! smmd: integers
      integer           :: ishift, row
      ! smmd: real
      real(kind=dp)     :: velconfac

      ! pdh: allocate buffers
      allocate(dbuf(nat*4),stat=ierr)
      call utils_alloc_check('internal_velocities','dbuf',ierr)
      dbuf(:) = 0.d0

      ishift=0 ! qoh: initialise on all processors
      ! aam: Root proc only: read absolute cartesian positions from
      !      <positions_abs> block of input file
      if (pub_on_root) then
         velconfac=1.0_dp
         if (esdf_block('velocities',dummy_nat)) then
            if (dummy_nat == nat+1) then
               read(block_data(1),'(a)') dummy_vel_unit
               velconfac = esdf_convfac(dummy_vel_unit,'auv')
               ishift=1
            elseif (dummy_nat.ne.nat) then
               call utils_abort('Error in internal_velocities(): mismatching &
                    &numbers of atoms in velocities specification.')
            endif
            do row=1,nat
               read(block_data(row+ishift),*) dbuf(row*3-2:row*3)
            end do
            dbuf(1:3*nat) = velconfac*dbuf(1:3*nat)
            md_init_velocities = .false.
         else
            md_init_velocities = .true.
         end if
      end if

      call comms_bcast(pub_root_proc_id,md_init_velocities)
      call comms_bcast(pub_root_proc_id,dbuf,nat*3)
      do row=1,nat
         elements(row)%ion_velocity(1:3) = dbuf(row*3-2:row*3)
      end do

      ! smmd: deallocate buffers
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_velocities','dbuf',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Read <velocities> block'

    end subroutine internal_velocities



    !------------------------------
    !--------- %block etrans_setup
    !------------------------------

    subroutine internal_etrans_setup()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_line_divide
      use rundat, only: pub_etrans_lcr
      use etrans, only: pub_lcr
      use utils, only: utils_abort, utils_assert

      implicit none

      !--------------
      character(len=80)    :: line, cjunk(2)
      integer           :: ibuf(2), maxp, nrow, il
      logical           :: use_etrans

      ! smmd: allocate buffers
      if (pub_on_root) then
         use_etrans = esdf_block('etrans_setup',nrow)

         if (use_etrans) then
            call utils_assert(nrow == 1, 'Error in internal_etrans_setup: &
                 &incorrect number of lines in etrans_setup block.')

            line = block_data(1)
            maxp = 2
            call esdf_line_divide(maxp,cjunk,line)
            do il=1,2
               read(cjunk(il),*) ibuf(il)
            enddo
         endif

      endif
      call comms_bcast(pub_root_proc_id,use_etrans)

      if (use_etrans) then

         call comms_bcast(pub_root_proc_id,ibuf,2)
         pub_lcr%atms(1) = ibuf(1)
         pub_lcr%atms(2) = ibuf(2)
         pub_lcr%atms(3) = ibuf(1)
         pub_lcr%atms(4) = ibuf(2)
         pub_lcr%ioread = .false.
         pub_lcr%iofile = 'none'
         pub_lcr%type   = 'lcr'

         ! check for errors
         call utils_assert(all(pub_lcr%atms(:) > 0), &
              & 'Transport LCR indices must be positive and non-zero!')
         call utils_assert(all(pub_lcr%atms(:) <= nat), &
              & 'Transport LCR indices must be smaller than number of atoms!')
         call utils_assert(pub_lcr%atms(2) > pub_lcr%atms(1), &
              & 'Transport LCR indices must be in ascending order!')
      else
         if (pub_etrans_lcr) then
            call utils_abort('etrans_setup block required for etrans_lcr calculation')
         endif
      endif

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Read <etrans_setup> block'

    end subroutine internal_etrans_setup



    !------------------------------
    !--------- %block etrans_leads
    !------------------------------

    subroutine internal_etrans_leads()

      use comms, only: pub_root_proc_id, comms_bcast
      use constants, only: DP, CRLF, ANGSTROM
      use esdf, only: block_data, esdf_block, esdf_line_divide, esdf_reduce
      use geometry, only: POINT, magnitude, operator(-)
      use rundat, only: pub_etrans_lcr, pub_etrans_bulk, &
           pub_etrans_lead_disp_tol, pub_etrans_lead_lcheck, &
           pub_etrans_seed_lead, pub_cond_calculate_any_task
      use etrans, only: etrans_init_leads, pub_nleads, pub_leads, pub_lcr
      use utils, only: utils_abort, utils_assert

      implicit none

      !--------------
      integer, allocatable :: ibuf(:,:), ibuf2(:)
      character(len=80)    :: line, cjunk(8), leadsname, msg, msg2, msg3
      character(len=80), allocatable :: cbuf(:)
      integer           :: ilead, nrow, il, maxp, iat, jat, idx, lead_nat
      logical           :: use_etrans
      type(POINT)       :: r0, r1
      real(kind=DP)     :: lead_length ! periodic length of lead
      real(kind=DP)     :: max_ngwf_rad
!      real(kind=DP), parameter :: lead_symm_disp_tol = 1.0d-3

      ! smmd: allocate buffers
      if (pub_on_root) then
         use_etrans = esdf_block('etrans_leads',nrow)
         if (.not. use_etrans) then
            nrow = 0
         elseif(use_etrans .and. nrow .lt. 1) then
            write(stdout,'(/a)') 'Error in internal_etrans_leads: &
                 &incorrect number of lines in'
            write(stdout,'(a)') 'leads specification'
            use_etrans = .not. use_etrans
         endif
      endif
      call comms_bcast(pub_root_proc_id,use_etrans)
      call comms_bcast(pub_root_proc_id,nrow)

      if (use_etrans) then
         allocate(ibuf(nrow,4),stat=ierr)
         call utils_alloc_check('internal_etrans_leads','ibuf',ierr)
         allocate(ibuf2(nrow),stat=ierr)
         call utils_alloc_check('internal_etrans_leads','ibuf2',ierr)
         allocate(cbuf(nrow),stat=ierr)
         call utils_alloc_check('internal_etrans_leads','cbuf',ierr)

         if (pub_on_root) then
            do ilead = 1, nrow
               line = block_data(ilead)
               maxp = 8

               ! initialise defaults
               cbuf(ilead) = 'none' ! File name
               ibuf2(ilead) = 1     ! lead unit cells

               cjunk(:) = ''

               call esdf_line_divide(maxp,cjunk,line)
               ! must have 4 entries
               if (maxp < 4) then
                  call utils_abort('Error in internal_etrans_leads: malformed entry&
                       & in etrans_lead block "'//block_data(ilead)//'": 4 block indices&
                       & required. Problem with lead',ilead)
               endif
               ! read block indices
               do il=1,4
                  read(cjunk(il),*) ibuf(ilead,il)
               enddo

               ! check remaining slots
               do il=5,7,2

                  ! check for a file name
                  idx = index(esdf_reduce(cjunk(il)),'file')
                  if (idx > 0) cbuf(ilead) = cjunk(il+1)

                  ! check for number of unit cells
                  idx = index(esdf_reduce(cjunk(il)),'unitcells')
                  if (idx > 0) read(cjunk(il+1),*) ibuf2(ilead)
               enddo
            end do
         end if


         call comms_bcast(pub_root_proc_id,ibuf)
         call comms_bcast(pub_root_proc_id,ibuf2)
         do ilead = 1, nrow
            call comms_bcast(pub_root_proc_id,cbuf(ilead))
         enddo

         pub_nleads = nrow
         call etrans_init_leads

         do ilead = 1, nrow
            pub_leads%type = 'lead'
            pub_leads(ilead)%atms(1:4) = ibuf(ilead,1:4)
            pub_leads(ilead)%num_unit_cells = ibuf2(ilead)
            leadsname = trim(adjustl(cbuf(ilead)))
            if (leadsname.ne.'none') then
               pub_leads(ilead)%ioread = .true.
               pub_leads(ilead)%iofile = leadsname
            else
               pub_leads(ilead)%ioread = .false.
               pub_leads(ilead)%iofile = 'none'
            endif
         end do

         ! rab207: check for errors
         do ilead=1,pub_nleads
            write(msg,'(a,i3,a)') 'Problem with lead', ilead,':'
            call utils_assert(pub_leads(ilead)%atms(2) > pub_leads(ilead)%atms(1), &
                 & trim(msg), CRLF//'First principle layer contains no atoms!')
            call utils_assert(pub_leads(ilead)%atms(4) > pub_leads(ilead)%atms(3), &
                 & trim(msg), CRLF//'Second principle layer contains no atoms!')
            call utils_assert(pub_leads(ilead)%atms(3) > pub_leads(ilead)%atms(2) .or. &
                 & pub_leads(ilead)%atms(1) > pub_leads(ilead)%atms(4), trim(msg), &
                 & CRLF//'The two principle layers must use different atoms!')
            call utils_assert(all(pub_leads(ilead)%atms(:) > 0), trim(msg), &
                 & CRLF//'Invalid atom indices: all indices must be > 0.')
            call utils_assert(all(pub_leads(ilead)%atms(:) <= nat), trim(msg), &
                 & CRLF//'Invalid atom indices: all indices must be < number of atoms.')
            call utils_assert(pub_leads(ilead)%atms(2)-pub_leads(ilead)%atms(1) == &
                 & pub_leads(ilead)%atms(4)-pub_leads(ilead)%atms(3), trim(msg), &
                 & CRLF//'The two principle layers must have the same number of atoms!')
            if (pub_etrans_lcr) then
               call utils_assert(all(pub_leads(ilead)%atms(:) <= maxval(pub_lcr%atms(:))) &
                           .and. all(pub_leads(ilead)%atms(:) >= minval(pub_lcr%atms(:))), &
                    & trim(msg), CRLF//'The leads principle layers must be contained within&
                    & the device region.')
            endif

            lead_nat = pub_leads(ilead)%atms(2) - pub_leads(ilead)%atms(1) + 1
            call utils_assert(modulo(lead_nat,pub_leads(ilead)%num_unit_cells)==0, &
                 trim(msg), CRLF//'The number of atoms do not divide evenly into &
                 &the requested number of unit cells!')

            ! check atoms in lead block and 1st principle layer are the same structure,
            ! in the same atomic order
            ! NOTE: this ignores the periodic boundary conditions, and currently spits
            ! out an error if the leads cross the cell boundary.
            iat = pub_leads(ilead)%atms(1)
            jat = pub_leads(ilead)%atms(3)
            ! a translation vector from lead block to 1st principle layer...
            r0 = elements(jat)%centre - elements(iat)%centre
            lead_length = magnitude(r0)

            ! ...check it is the same as all other translation vectors
            do il = 1, pub_leads(ilead)%atms(2) - pub_leads(ilead)%atms(1)
               r1 = elements(jat+il)%centre - elements(iat+il)%centre
               lead_length = max(lead_length,magnitude(r1))

               call utils_assert( magnitude(r1-r0) < pub_etrans_lead_disp_tol, trim(msg), &
                      CRLF//'Atoms in the two lead principles layers must have the same geometry'// &
                    & CRLF//'(i.e. be periodic repeats of each other) and be in the same atomic'//&
                    & CRLF//'ordering. Note that the internal routine that checks the geometry'// &
                    & CRLF//'and the ordering cannot deal with leads that cross the cell boundary.'// &
                    & CRLF//'Increase etrans_lead_disp_tol if you wish to ignore this warning.')
            enddo

            ! check that the lead size is sufficient
            jat = pub_leads(ilead)%atms(2)
            max_ngwf_rad = maxval(elements(iat:jat)%radius)
            if (pub_cond_calculate_any_task) then
               max_ngwf_rad = max(max_ngwf_rad,maxval(elements(iat:jat)%radius_cond))
            endif

            write(msg2,'(f8.2,a)') 2.0_DP*max_ngwf_rad*ANGSTROM, ' Angstrom'
            write(msg3,'(f8.2,a)') lead_length*ANGSTROM, ' Angstrom'
            if (pub_etrans_lead_lcheck) then
               call utils_assert(lead_length > 2.0_DP*max_ngwf_rad, trim(msg)// &
                    CRLF//'The principle layer periodic length is shorter than the maximum NGWF'//&
                    CRLF//'diameter and so does not form a complete principle layer. Increase the'//&
                    CRLF//'number of primitive unit cells to avoid truncation errors, or override'//&
                    CRLF//'with etrans_lead_size_check = .false.'//&
                    CRLF//&
                    CRLF//'The minimum principle layer size for the NGWF radii is '//trim(adjustl(msg2))//&
                    CRLF//'The current principle layer size is '//trim(adjustl(msg3)))
            endif
            if (lead_length < 2.0_DP*max_ngwf_rad .and. pub_on_root) then
               write(stdout,'(/2a/)') 'WARNING: ', trim(msg)// &
                    CRLF//'The principle layer periodic length is shorter than the maximum NGWF'//&
                    CRLF//'diameter and so does not form a complete principle layer. Increase the'//&
                    CRLF//'number of primitive unit cells to avoid truncation errors.'//&
                    CRLF//&
                    CRLF//'The minimum principle layer size for the NGWF radii is '//trim(adjustl(msg2))//&
                    CRLF//'The current principle layer size is '//trim(adjustl(msg3))
            endif

            ! check seed lead is valid
            call utils_assert(pub_etrans_seed_lead <= pub_nleads .and. &
                 pub_etrans_seed_lead > 0, 'Error: etrans_seed_lead must be &
                 &between 1 and the number of leads')

            ! check that the lead unit cells have the same structure
            if (pub_leads(ilead)%num_unit_cells > 1) then
               lead_nat = pub_leads(ilead)%atms(2) - pub_leads(ilead)%atms(1) + 1
               lead_nat = lead_nat / pub_leads(ilead)%num_unit_cells
               iat = pub_leads(ilead)%atms(1)
               jat = pub_leads(ilead)%atms(2)
               r0 = elements(iat+lead_nat)%centre - elements(iat)%centre

               do il=iat+lead_nat,jat
                  r1 = elements(il)%centre - elements(il-lead_nat)%centre
                  call utils_assert( magnitude(r1-r0) < pub_etrans_lead_disp_tol, trim(msg)// &
                       CRLF//'To symmetrise the lead unit cells, the geometry of all primitive unit'//&
                       CRLF//'cells must be identical with the atoms in the same atomic ordering.'//&
                       CRLF//'Note that the internal routine that checks the atomic ordering cannot'//&
                       CRLF//'deal with principle layers that cross the cell boundary.'//&
                       CRLF//'Increase etrans_lead_disp_tol if you wish to ignore this warning.')
               enddo
            endif

         enddo


         ! smmd: deallocate buffers
         deallocate(ibuf,stat=ierr)
         call utils_dealloc_check('internal_etrans_leads','ibuf',ierr)
         deallocate(ibuf2,stat=ierr)
         call utils_dealloc_check('internal_etrans_leads','ibuf2',ierr)
         deallocate(cbuf,stat=ierr)
         call utils_dealloc_check('internal_etrans_leads','cbuf',ierr)
      else if(pub_etrans_lcr .or. pub_etrans_bulk) then
         call utils_abort('etrans_leads block required for etrans_lcr or etrans_bulk&
              & calculation')
      endif

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <etrans_leads> block'

    end subroutine internal_etrans_leads

    !------------------------------------
    !--------- %block etrans_eigenchannel_energies
    !------------------------------------

    subroutine internal_etrans_plot_eigchan_energies()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: esdf_convfac, block_data, esdf_block, esdf_reduce
      use etrans, only: pub_etrans_eigchan_en, pub_etrans_plot_num

      implicit none

      !---------------
      logical :: use_etrans
      integer :: nrow, irow, ishift
      integer :: ierr
      real(kind=DP) :: enconvfac
      character(len=80) :: dummy_en_unit

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <etrans_plot_eigchan_energies> block'

      if (pub_on_root) then
         use_etrans = esdf_block('etrans_plot_eigchan_energies',nrow)


         enconvfac = 1.0_DP
         ishift = 0
         ! check if units on first line
         if (nrow > 0 .and. (index(esdf_reduce(block_data(1)),'ha')>0 .or. &
              index(esdf_reduce(block_data(1)),'ev')>0) ) then
            nrow = nrow - 1
            ishift = 1
            read(block_data(1),'(a)') dummy_en_unit
            enconvfac = esdf_convfac(dummy_en_unit,'ha')
         endif
         if (nrow == 0) use_etrans = .false.
      endif
      call comms_bcast(pub_root_proc_id,use_etrans)
      call comms_bcast(pub_root_proc_id,nrow)
      pub_etrans_plot_num = nrow

      if (use_etrans) then
         ! allocate space for the energies
         allocate(pub_etrans_eigchan_en(nrow),stat=ierr)
         call utils_alloc_check('internal_etrans_ldos','pub_etrans_eigchan_en',ierr)

         ! root reads the energies
         if (pub_on_root) then
            do irow=1,nrow
               read(block_data(irow+ishift),*) pub_etrans_eigchan_en(irow)
            enddo
            ! convert units
            pub_etrans_eigchan_en(:) = pub_etrans_eigchan_en(:) * enconvfac
         endif

         ! broadcast
         call comms_bcast(pub_root_proc_id,pub_etrans_eigchan_en)
      endif

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <etrans_plot_eigchan_energies> block'

    end subroutine internal_etrans_plot_eigchan_energies

    !------------------------------------
    !--------- %block thermostat
    !------------------------------------

    subroutine internal_thermostat()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_convfac
      use md_thermostat, only: md_thermo, md_thermo_num
      use rundat, only: md_delta_t
      use utils, only: utils_abort

      implicit none

      !-------------
      real(kind=DP)    :: io_temp, io_tgrad
      real(kind=DP)    :: io_tau, io_mix, io_damp
      integer          :: io_kind, io_start, io_stop
      integer          :: io_group, io_nhc, io_nhi
      logical          :: io_upd
      integer          :: nthermo
      !-------------
      integer, allocatable :: ibuf(:)
      logical, allocatable :: lbuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      integer :: row, nrow
      integer :: ith, ierr
      character (len=80) :: cjunk
      !character (len=4) :: thermotype(4)
      logical :: thermo
      !-------------
      !data thermotype /'none','andersen','langevin','nosehoover'/

      ! smmd: allocate buffers
      if (pub_on_root) thermo = esdf_block('thermostat',nrow)
      call comms_bcast(pub_root_proc_id,nrow)

      if (nrow .ge. 1) then
         allocate(lbuf(nrow),stat=ierr)
         call utils_alloc_check('internal_thermostat','lbuf',ierr)
         allocate(dbuf(nrow*5),stat=ierr)
         call utils_alloc_check('internal_thermostat','dbuf',ierr)
         allocate(ibuf(nrow*6),stat=ierr)
         call utils_alloc_check('internal_thermostat','ibuf',ierr)
         dbuf(:) = 0.0_dp
         ibuf(:) = 0
      endif

      ! smmd: Root proc only: read information from <thermostat>
      ! block of input file
      ierr = 0
      if (pub_on_root .and. nrow .ge. 1) then
         if (esdf_block('thermostat',nrow)) then
            nthermo = 0
            row = 1
            do while (row .le. nrow)

               ! Parse each user defined thermostat
               cjunk = block_data(row)

               if (index(cjunk,'none').gt.0 .or. index(cjunk,'andersen').gt.0 .or.&
                    index(cjunk,'langevin').gt.0 .or. index(cjunk,'nosehoover').gt.0 .or.&
                    index(cjunk,'berendsen').gt.0 .or. index(cjunk,'bussi').gt.0) then

                  ! Read mandatory parameters
                  call thermostat_read_params(cjunk, &
                       io_temp,io_kind,io_start,io_stop)
                  nthermo = nthermo + 1
                  row = row + 1

                  ! Default parameters
                  io_tgrad  = 0.0_dp
                  io_tau    = 10.0_dp*md_delta_t
                  io_mix    = 1.0_dp
                  io_damp   = 0.2_dp
                  io_group  = 0
                  io_nhc    = 4
                  io_nhi    = 10
                  io_upd    = .false.

                  ! User defined parameters
                  options_loop : do
                     cjunk = block_data(row)
                     if (index(cjunk,'none').gt.0 .or. index(cjunk,'andersen').gt.0 .or.&
                          index(cjunk,'langevin').gt.0 .or. index(cjunk,'nosehoover').gt.0 .or.&
                          index(cjunk,'berendsen').gt.0 .or. index(cjunk,'bussi').gt.0 .or.&
                          row .gt. nrow) then
                        exit options_loop
                     else
                        ! Look for optional parameters
                        call thermostat_read_options(cjunk,io_tgrad, &
                             io_tau,io_mix,io_damp,io_group,io_nhc,io_nhi,io_upd)
                        row = row + 1
                     endif
                  enddo options_loop

               else
                  call utils_abort('Error in internal_thermostat(): &
                       &wrong thermostat definition in input file.')
               endif

               ! Load buffer arrays
               dbuf((nthermo-1)*5+1:nthermo*5) = &
                    (/io_temp,io_tgrad,io_tau,io_mix,io_damp/)
               ibuf((nthermo-1)*6+1:nthermo*6) = &
                    (/io_kind,io_start,io_stop,io_group,io_nhc,io_nhi/)
               lbuf(nthermo) = io_upd

            enddo
         endif
      endif

      if (nrow .ge. 1) then
         ! smmd: Broadcast buffers
         call comms_bcast(pub_root_proc_id,nthermo)
         call comms_bcast(pub_root_proc_id,dbuf,nrow*5)
         call comms_bcast(pub_root_proc_id,ibuf,nrow*6)
         call comms_bcast(pub_root_proc_id,lbuf,nrow)

         ! smmd: allocate public arrays in thermostat_mod
         md_thermo_num = nthermo
         allocate(md_thermo(md_thermo_num),stat=ierr)
         call utils_alloc_check('internal_thermostat','md_thermo',ierr)

         ! smmd: process thermostat parameters
         do ith=1,md_thermo_num

            md_thermo(ith)%type   = ibuf((ith-1)*6+1)
            md_thermo(ith)%start  = ibuf((ith-1)*6+2)
            md_thermo(ith)%stop   = ibuf((ith-1)*6+3)
            md_thermo(ith)%group  = ibuf((ith-1)*6+4)
            md_thermo(ith)%nhc_length = ibuf((ith-1)*6+5)
            md_thermo(ith)%nhc_integ_nstep = ibuf((ith-1)*6+6)

            md_thermo(ith)%Tinit = dbuf((ith-1)*5+1)
            md_thermo(ith)%Tgrad = dbuf((ith-1)*5+2)
            md_thermo(ith)%tau  = dbuf((ith-1)*5+3)
            md_thermo(ith)%mix   = dbuf((ith-1)*5+4)
            md_thermo(ith)%damp  = dbuf((ith-1)*5+5)

            md_thermo(ith)%nhc_upd = lbuf(ith)

         enddo

         ! smmd: deallocate buffers
         deallocate(ibuf,stat=ierr)
         call utils_dealloc_check('internal_thermostat','ibuf',ierr)
         deallocate(dbuf,stat=ierr)
         call utils_dealloc_check('internal_thermostat','dbuf',ierr)
         deallocate(lbuf,stat=ierr)
         call utils_dealloc_check('internal_thermostat','lbuf',ierr)
      endif

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <thermostat> block'

    end subroutine internal_thermostat


    !------------------------
    !--------- %block_species
    !------------------------

    subroutine internal_species

      use comms, only: pub_root_proc_id, comms_bcast
      use constants, only: CRLF, kernel_cutoff_sanity_threshold
      use esdf, only: block_data, esdf_block, esdf_convfac
      use rundat, only: pub_ngwf_halo, pub_kernel_cutoff
      use utils, only: utils_abort, utils_assert, utils_flushed_string_output

      implicit none

      !-------------
      ! ars: dummies
      character(len=80) :: dummy_len_unit
      character(len=4)  :: dummy_id
      character(len=2)  :: dummy_symb
      integer :: dummy_nsp, dummy_nfunc, dummy_Z
      real(kind=dp) :: dummy_rad
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      integer, allocatable :: ibuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      character(len=256) :: message
      character(len=10) :: valuestr, valuestr2
      !-------------
      ! ars: integers
      integer :: ishift, row, ii
      integer :: iostat
      ! ars: logical
      logical :: element_found
      ! ars: real
      real(kind=dp) :: lenconfac
      real(kind=DP) :: max_ngwf_rad
      character(len=*), parameter :: myself = 'internal_species'

      ! -----------------------------------------------------------------------

      ! Allocate array if not already allocated, even if size is zero
      if (.not.allocated(species)) then
         allocate(species(nsp),stat=ierr)
         call utils_alloc_check('internal_species','species',ierr)
         ! jd: Initialise species member hc_steric_cutoff to zero
         species(1:nsp)%hc_steric_cutoff = 0D0
      else if (size(species)/=nsp) then
         call utils_abort('Error in internal_species &
              &(rundat_blocks_read) size of species array does &
              &not match number of species')
      end if

      ! pdh: allocate buffers
      ! ddor: Altered for DFT+U
      allocate(dbuf(max(nat*4,9)),stat=ierr)
      call utils_alloc_check('internal_species','dbuf',ierr)
      allocate(ibuf(nat*4),stat=ierr)
      call utils_alloc_check('internal_species','ibuf',ierr)
      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_species','cbuf',ierr)

      ! aam: Root proc only: read information from <species> block of input file
      if (pub_on_root) then
         lenconfac=1.0_dp
         ishift=0
         if (esdf_block('species',dummy_nsp)) then
            ! aam: 3/6/09 added ability to read optional  unit string
            !      ("bohr" or "ang") in first line of species block
            !      (default is bohr)
            if (dummy_nsp == nsp+1) then
               read(block_data(1),'(a)') dummy_len_unit
               lenconfac = esdf_convfac(dummy_len_unit,'bohr')
               ishift=1
            endif
            if ( (dummy_nsp.ne.nsp) .and. (dummy_nsp.ne.nsp+1) ) then
               call utils_abort('Error in internal_species: mismatching &
                    &numbers of species in <species> specification.'//CRLF//&
                    'Perhaps your positions block, or input .xyz file (if &
                    &using one) contains species that are not listed in the &
                    &species block.')
            endif
            ! aam: check for multiply defined species
            call internal_check_species(nsp,'<species>') !ddor
            do row=1,nsp
               read(block_data(row+ishift),*,iostat=ierr) dummy_id,dummy_symb, &
                    ibuf(row*2-1:row*2),dbuf(row)
               ! jd: Catch the corner case where an extra species in the
               !     positions block is initially missed because there's a
               !     unit specification in the species block.
               !     If unchecked, this leads to a crash with an end-of-file
               !     error reported by the RTL.
               call utils_assert(ierr==0, 'Error in '//myself//': &
                    &detected reading the <species> specification. '//CRLF//&
                    'Perhaps your positions block, or input .xyz file (if &
                    &using one) contains species that are not listed in the &
                    &species block.')
               cbuf(row)(1:4) = dummy_id
               cbuf(row)(5:6) = dummy_symb
               ! ars: save species_id's in buffer
               all_species_id(row)(1:4) = cbuf(row)(1:4)
            end do
            dbuf(1:nsp)=lenconfac*dbuf(1:nsp)
         else
            call utils_abort('Error in internal_species: &
                 &<species> block not found in input file.')
         end if
      end if

      ! pdh: broadcast species information
      call comms_bcast(pub_root_proc_id,ibuf,nsp*2)
      call comms_bcast(pub_root_proc_id,dbuf,nsp)
      do row=1,nsp
         call comms_bcast(pub_root_proc_id,cbuf(row),6)
      end do

      ! extract element information from species
      element_found = .false. ! qoh: Initialise to avoid compiler warnings
      do row=1,nsp
         element_found = .false. ! aam: error check flag
         dummy_id = cbuf(row)(1:4)
         dummy_symb = cbuf(row)(5:6)
         dummy_Z = ibuf(row*2-1)
         dummy_nfunc = ibuf(row*2)
         dummy_rad = dbuf(row)
         species(row)%species_id = dummy_id
         species(row)%symbol = dummy_symb
         species(row)%nfunctions = dummy_nfunc
         species(row)%atomic_number = dummy_Z
         do ii=1,nat
            if (elements(ii)%species_id == dummy_id) then
               elements(ii)%symbol         = dummy_symb
               elements(ii)%atomic_number  = dummy_Z
               elements(ii)%nfunctions     = dummy_nfunc
               elements(ii)%species_number = row
               ! rc2013: label for species numbers across all subsystems
               elements(ii)%global_species_number = row
               elements(ii)%radius         = dummy_rad
               species(row)%region = elements(ii)%region
               element_found = .true.
            end if
         end do
         if (.not. element_found) then
            call utils_abort('Error in internal_species: &
                 &mismatching species in <species> and <positions_abs> &
                 &blocks of input file.')
         end if
      end do

      ! cks: 21/04/2006 Warning! If there is a halo, the element radius
      ! cks: is the NGWF radius plus the halo length and so it is not
      ! cks: equal to the NGWF-sphere radius (which does not include the halo)
      if (pub_ngwf_halo > 0.0_DP) then
         ! cks: increase element radius by halo length
         do ii=1,nat
            elements(ii)%radius = elements(ii)%radius + pub_ngwf_halo
         end do
      end if

      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_species','cbuf',ierr)
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_species','dbuf',ierr)
      deallocate(ibuf,stat=ierr)
      call utils_dealloc_check('internal_species','ibuf',ierr)

      ! jd: Check against unreasonably low kernel cutoffs
      max_ngwf_rad = maxval(elements(:)%radius)
      if(pub_kernel_cutoff <= kernel_cutoff_sanity_threshold) then
         write(valuestr,'(f10.3)') pub_kernel_cutoff
         write(message,'(a,a,a)') &
              'WARNING: The kernel cutoff you specified (', trim(adjustl(valuestr)), &
              ' bohr) seems unreasonably low.'
         call utils_flushed_string_output(message//CRLF)
      end if
      if(pub_kernel_cutoff <= 2.0_DP*max_ngwf_rad) then
         write(valuestr,'(f10.3)') pub_kernel_cutoff
         write(valuestr2,'(f10.3)') max_ngwf_rad
         write(message,'(a,a,a,a,a)') &
              'The kernel cutoff (you specified ',trim(adjustl(valuestr)), &
              ' bohr) must be larger than twice the maximum NGWF radius&
              & (which you specified to be ',trim(adjustl(valuestr2)),' bohr).'
         call utils_abort(trim(message))
      end if

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Read <species> block'

    end subroutine internal_species


    !------------------------
    !--------- %block_sol_ions
    !------------------------

    ! jd: by Jacek Dziedzic in 07/2014, largely cloned from internal_species.
    subroutine internal_sol_ions

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use is_solvation_boltzmann, only: &
           implicit_solvent_boltzmann_init_ions, &
           implicit_solvent_boltzmann_max_n_ions
      use rundat, only: pub_is_pbe
      use utils, only: utils_abort, utils_assert, utils_int_to_str

      implicit none

      integer            :: n_sol_ions
      ! jd: Names of ions
      character(len=10)  :: cbuf(implicit_solvent_boltzmann_max_n_ions)
      ! jd: Charges, concentrations, optional neutralisation ratios
      real(kind=DP)      :: dbuf(implicit_solvent_boltzmann_max_n_ions*3)
      integer            :: row
      integer            :: io
      character(len=4)   :: dummy_id

      ! ----------------------------------------------------------------------

      if(pub_is_pbe == 'NONE') return

      ! jd: Take care of padding between n_sol_ions and max_sol_ions
      dbuf(:) = 0D0; cbuf(:) = '          '

      ! jd: On root read the block and parse it
      if (pub_on_root) then
         if (esdf_block('sol_ions',n_sol_ions)) then
            call utils_assert(&
                 n_sol_ions <= implicit_solvent_boltzmann_max_n_ions, &
                 'Too many types of ions in %block sol_ions')
            do row=1, n_sol_ions
               ! jd: Read ID, charge, concentration. These must be present.
               read(block_data(row),*,iostat=io) dummy_id, dbuf(row*3-2:row*3-1)
               call utils_assert(io == 0, 'Error parsing row #'//&
                    trim(utils_int_to_str(row))//' of %block sol_ions. Make sur&
                    &e at least three columns are present (name, charge, concen&
                    &tration) and that the 2nd and 3rd are valid numbers.')
               read(block_data(row),*,iostat=io) dummy_id, dbuf(row*3-2:row*3)
               if(io < 0) then
                  ! EOF on 3rd field -- optional neutralisation ratio absent.
                  dbuf(row*3) = 0.0_DP
               end if
               call utils_assert(io <= 0, 'Error parsing row #'//&
                    trim(utils_int_to_str(row))//' of %block sol_ions. Looks li&
                    &ke the optional 4th column (neutralisation ratio) is not a&
                    & number.')
               cbuf(row)(1:10) = trim(dummy_id)
            end do
         else
            call utils_abort('Error in internal_sol_ions: &
                 &when is_pbe is used, "%block sol_ions" must be specified.')
         end if
      end if

      ! jd: Broadcast parsed block contents
      call comms_bcast(pub_root_proc_id,n_sol_ions)
      call comms_bcast(pub_root_proc_id,dbuf,n_sol_ions*3)
      do row=1,n_sol_ions
         call comms_bcast(pub_root_proc_id,cbuf(row),10)
      end do

      call implicit_solvent_boltzmann_init_ions(cbuf,dbuf,n_sol_ions)

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Read <sol_ions> block'

    end subroutine internal_sol_ions


    !------------------------------------
    !--------- %block species_solvent_radius
    !------------------------------------

    !ab: adapted from other species specific blocks

    subroutine internal_species_solvent_radius()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_convfac
      use utils, only: utils_abort, utils_assert
      use rundat, only: pub_is_pbe
      use constants, only: CRLF

      implicit none

      !-------------
      ! dummies
      character(len=80) :: dummy_rad_unit
      real(kind=dp)     :: dummy_rad
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      !-------------
      ! buffers
      character(len=68), allocatable :: cbuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      integer :: ishift, row, ii, isp, ierr
      logical :: species_found, element_found
      real(kind=dp) :: lenconfac
      character(len=*), parameter :: myself='internal_species_solvent_radius'

      ! ab: no role of solvent radius when there are no electrolyte ions
      if(pub_is_pbe == 'NONE') return

      ! allocate buffers
      allocate(dbuf(nsp),stat=ierr)
      call utils_alloc_check(myself,'dbuf',ierr)
      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check(myself,'cbuf',ierr)

      ! read information from <species_solvent_radius> block of input file
      ! on root proc only
      if (pub_on_root) then
         lenconfac=1.0_dp
         ishift=0
         if (esdf_block('species_solvent_radius',dummy_nsp)) then
            if (dummy_nsp == nsp+1) then
                read(block_data(1),'(a)') dummy_rad_unit
                lenconfac = esdf_convfac(dummy_rad_unit,'bohr')
                ishift=1
            endif
            call utils_assert((dummy_nsp.eq.nsp).or.(dummy_nsp.eq.nsp+1),&
                 'Error in internal_species_solvent_radius: mismatching &
                 &numbers of species in <species_solvent_radius> &
                 &specification.'//CRLF//&
                 'Perhaps your positions block, or input .xyz file (if &
                 &using one) contains species that are not listed in the &
                 &species_solvent_radius block.')

            ! check for multiply defined species
            call internal_check_species(nsp,'<species_solvent_radius>')
            do row=1,nsp
               read(block_data(row+ishift),*,iostat=ierr) dummy_id,dbuf(row)
               ! jd: Catch the corner case where an extra species in the
               ! positions block is initially missed because there's a
               ! unit specification in the block.
               ! If unchecked, this leads to a crash with an end-of-file
               ! error reported by the RTL.
               call utils_assert(ierr==0, 'Error in '//myself//': &
                    &reading the <species_solvent_radius> specification. '//&
                    CRLF//&
                    'Perhaps your positions block, or input .xyz file (if &
                    &using one) contains species that are not listed in the &
                    &species block.')
               cbuf(row)(1:4) = dummy_id
            end do
            dbuf(1:nsp)=lenconfac*dbuf(1:nsp)
         end if
      end if

      ! broadcast
      do row=1,nsp
         call comms_bcast(pub_root_proc_id,cbuf(row),4)
      end do
      call comms_bcast(pub_root_proc_id,dbuf,nsp)

     ! extract solvent_radius
      do row=1,nsp
         species_found = .false. ! error check flag
         element_found = .false. ! error check flag
         dummy_id = cbuf(row)(1:4)
         dummy_rad = dbuf(row)
         do isp=1,nsp
            if (species(isp)%species_id == dummy_id) then
               species(isp)%solvent_radius = dummy_rad
               species_found=.true.
            end if
         end do
         if (.not. species_found) then
            call utils_abort('Error in '//myself//': &
                 &mismatching species in <species_solvent_radius> and &
                 &<positions_abs> blocks of input file.')
         end if
         do ii=1,nat
            if (elements(ii)%species_id == dummy_id) then
               elements(ii)%solvent_radius = dummy_rad
               element_found = .true.
            end if
         end do
         if (.not. element_found) then
            call utils_abort('Error in '//myself//': &
                 &mismatching species in <species_solvent_radius> and &
                 &<positions_abs> blocks of input file.')
         end if
      end do

      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check(myself,'cbuf',ierr)
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check(myself,'dbuf',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Read <species_solvent_radius> block'

    end subroutine internal_species_solvent_radius


    !------------------------
    !--------- %block_is_dielectric_exclusions
    !------------------------

    ! jd: by Jacek Dziedzic in 08/2017, largely cloned from internal_thole_polarisabilities.
    subroutine internal_is_dielectric_exclusions

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_is_implicit_solvent
      use is_smeared_ions, only: dielectric_exclusions, n_max_exclusions
      use utils, only: utils_assert

      implicit none

      integer            :: n_exclusions
      integer            :: row
      logical            :: is_dielectric_exclusions_exist
      character(len=20)  :: exclusion_type, dummy

      ! ----------------------------------------------------------------------

      dielectric_exclusions(:,:) = 0D0
      if(.not. pub_is_implicit_solvent) return

      if(pub_on_root) then
         write(stdout,'(a)') 'IS: Dielectric excluded regions:'
         is_dielectric_exclusions_exist = esdf_block('is_dielectric_exclusions', &
              n_exclusions)
         call utils_assert(n_exclusions <= n_max_exclusions, &
              'Error in internal_is_dielectric_exclusions: too many exclusions')
      end if

      call comms_bcast(pub_root_proc_id, is_dielectric_exclusions_exist)
      call comms_bcast(pub_root_proc_id, n_exclusions)

      if(.not. is_dielectric_exclusions_exist .or. n_exclusions == 0) then
         if(pub_on_root) write(stdout, '(a)') '    - NONE.'
         return
      end if

      ! jd: On root read the block and parse it
      if (pub_on_root) then
         do row=1, n_exclusions
            read(block_data(row),*) exclusion_type
            if(trim(exclusion_type) == 'sphere' .or. &
                 trim(exclusion_type) == 'SPHERE') then
               dielectric_exclusions(row,1) = 1
               read(block_data(row),*) dummy, dielectric_exclusions(row,2:5)
               write(stdout,'(a,f9.3,a,f9.3,a,f9.3,a,f9.3)') &
                    '    - SPHERE @ ', dielectric_exclusions(row,2), &
                    ', ', dielectric_exclusions(row,3), &
                    ', ', dielectric_exclusions(row,4), &
                    ' r = ', dielectric_exclusions(row,5)
               call utils_assert(dielectric_exclusions(row,5) >= 0D0, &
                    'Illegal sphere dielectric exclusion: radius < 0.')
            else if(trim(exclusion_type) == 'box' .or. &
                 trim(exclusion_type) == 'BOX') then
               dielectric_exclusions(row,1) = 2
               read(block_data(row),*) dummy, dielectric_exclusions(row,2:7)
               write(stdout,'(a,f9.3,a,f9.3,a,f9.3,a,f9.3,a,f9.3,a,f9.3,a)') &
                    '    - BOX    @ ', dielectric_exclusions(row,2), &
                    ', ', dielectric_exclusions(row,4), &
                    ', ', dielectric_exclusions(row,6), &
                    ' - ', dielectric_exclusions(row,3), &
                    ', ', dielectric_exclusions(row,5), &
                    ', ', dielectric_exclusions(row,7),'.'
               call utils_assert(dielectric_exclusions(row,2) <= dielectric_exclusions(row,3), &
                    'Illegal box dielectric exclusion: xmin larger than xmax.')
               call utils_assert(dielectric_exclusions(row,4) <= dielectric_exclusions(row,5), &
                    'Illegal box dielectric exclusion: ymin larger than ymax.')
               call utils_assert(dielectric_exclusions(row,6) <= dielectric_exclusions(row,7), &
                    'Illegal box dielectric exclusion: zmin larger than zmax.')
            else if(trim(exclusion_type) == 'xcyl' .or. &
                 trim(exclusion_type) == 'XCYL') then
               dielectric_exclusions(row,1) = 3
               read(block_data(row),*) dummy, dielectric_exclusions(row,2:4)
               write(stdout,'(a,f9.3,a,f9.3,a,f9.3)') &
                    '    - XCYL @ Y= ', dielectric_exclusions(row,2), &
                    ', Z= ', dielectric_exclusions(row,3), &
                    ' r = ', dielectric_exclusions(row,4)
               call utils_assert(dielectric_exclusions(row,4) >= 0D0, &
                    'Illegal xcyl dielectric exclusion: radius < 0.')
            else if(trim(exclusion_type) == 'ycyl' .or. &
                 trim(exclusion_type) == 'YCYL') then
               dielectric_exclusions(row,1) = 4
               read(block_data(row),*) dummy, dielectric_exclusions(row,2:4)
               write(stdout,'(a,f9.3,a,f9.3,a,f9.3)') &
                    '    - XCYL @ X= ', dielectric_exclusions(row,2), &
                    ', Z= ', dielectric_exclusions(row,3), &
                    ' r = ', dielectric_exclusions(row,4)
               call utils_assert(dielectric_exclusions(row,4) >= 0D0, &
                    'Illegal ycyl dielectric exclusion: radius < 0.')
            else if(trim(exclusion_type) == 'zcyl' .or. &
                 trim(exclusion_type) == 'ZCYL') then
               dielectric_exclusions(row,1) = 5
               read(block_data(row),*) dummy, dielectric_exclusions(row,2:4)
               write(stdout,'(a,f9.3,a,f9.3,a,f9.3)') &
                    '    - ZCYL @ X= ', dielectric_exclusions(row,2), &
                    ', Y= ', dielectric_exclusions(row,3), &
                    ' r = ', dielectric_exclusions(row,4)
               call utils_assert(dielectric_exclusions(row,4) >= 0D0, &
                    'Illegal zcyl dielectric exclusion: radius < 0.')
            else
               call utils_assert(.false., 'Unrecognized dielectric exclusion &
                    &type: "'//trim(exclusion_type)//'".')
            end if
         end do
      end if

      ! jd: Broadcast parsed block contents.
      call comms_bcast(pub_root_proc_id,dielectric_exclusions)

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Read <is_dielectric_exclusions> block'

    end subroutine internal_is_dielectric_exclusions

    !------------------------
    !--------- %block_thole_polarisabilities
    !------------------------

    ! jd: by Jacek Dziedzic in 11/2015, largely cloned from internal_species.
    subroutine internal_thole_polarisabilities

      use comms, only: pub_root_proc_id, comms_bcast
      use constants, only: garbage_real
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_pol_emb_qmstar, pub_pol_emb_fixed_charge, &
           pub_pol_emb_polscal
      use utils, only: utils_assert, utils_int_to_str

      implicit none

      integer            :: n_polarisabilities
      logical            :: thole_polarisabilities_exist
      ! jd: Names of ions, will be ignored, except in sanity check
      character(len=10), allocatable  :: cbuf(:)
      ! jd: Polarisabilities
      real(kind=DP), allocatable      :: dbuf(:)
      integer            :: row
      character(len=4)   :: dummy_id

      ! ----------------------------------------------------------------------

      elements(:)%thole_polarisability = -garbage_real
      if(.not. pub_pol_emb_qmstar) return
      if(pub_pol_emb_fixed_charge) return ! Thole block absent in fixed charge FFs

      if(pub_on_root) then
         thole_polarisabilities_exist = esdf_block('thole_polarisabilities', &
              n_polarisabilities)

         call utils_assert(thole_polarisabilities_exist, &
              'Error in internal_thole_polarisabilities: when polarisable &
              &embedding is used (pol_emb_pot_filename), &
              &"%block thole_polarisabilities" must be specified, unless the &
              &force field is not polarisable (pol_emb_fixed_charge T) or &
              &the QM* description is never used (pol_emb_qmstar F)')
      end if

      call comms_bcast(pub_root_proc_id,n_polarisabilities)

      call utils_assert(n_polarisabilities == size(elements), &
           'Error in internal_thole_polarisabilities: The number of specified &
           &Thole polarisabilities ('//&
           trim(utils_int_to_str(n_polarisabilities))//&
           ') is different from the number of atoms in the system ('//&
           trim(utils_int_to_str(size(elements)))//')')

      allocate(cbuf(n_polarisabilities),stat=ierr)
      call utils_alloc_check('internal_thole_polarisabilities','cbuf',ierr)
      allocate(dbuf(n_polarisabilities),stat=ierr)
      call utils_alloc_check('internal_thole_polarisabilities','dbuf',ierr)

      ! jd: On root read the block and parse it
      if (pub_on_root) then
         do row=1, n_polarisabilities
            read(block_data(row),*) dummy_id, dbuf(row)
            cbuf(row)(1:10) = trim(dummy_id)

            ! jd: Most rudimentary sanity check -- ensure that the first letter of
            !     the symbol in elements matches the first letter specified in the
            !     polarisability block. This will trap obvious mismatches, while
            !     still allowing N1 to match N and so on.
            call utils_assert(cbuf(row)(1:1) == elements(row)%symbol(1:1), &
                 'Error in internal_thole_polarisabilities: Mismatch between &
                 &the symbol of the polarisability of atom '//&
                 trim(utils_int_to_str(row))//': "'//trim(cbuf(row))//&
                 '" and the symbol inferred from the positions block: "'//&
                 trim(elements(row)%symbol)//".")

         end do
      end if

      ! jd: Broadcast parsed block contents.
      !     Ignore cbuf. Store dbuf in elements.
      !     Enforce polscal, if any.
      call comms_bcast(pub_root_proc_id,dbuf,n_polarisabilities)
      do row = 1, n_polarisabilities
         elements(row)%thole_polarisability = dbuf(row) * pub_pol_emb_polscal
      end do

      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_thole_polarisabilities','dbuf',ierr)
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_thole_polarisabilities','cbuf',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Read <thole_polarisabilities> block'

    end subroutine internal_thole_polarisabilities

    !------------------------
    !--------- %block_mm_rep_params
    !------------------------

    ! jd: by Jacek Dziedzic in 03/2018,
    !     largely cloned from internal_thole_polarisabilities.
    subroutine internal_mm_rep_params

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use polarisable_embedding, only: polarisable_embedding_init_mm_rep_params
      use rundat, only: pub_pol_emb_qmstar

      implicit none

      integer            :: n_mm_species
      logical            :: mm_rep_params_exist
      ! jd: Names of ions
      character(len=10), allocatable  :: cbuf(:)
      ! jd: MM params
      real(kind=DP), allocatable      :: dbuf(:,:)
      integer            :: row
      character(len=4)   :: dummy_id

      ! ----------------------------------------------------------------------

      if(.not. pub_pol_emb_qmstar) return

      if(pub_on_root) then
         mm_rep_params_exist = esdf_block('mm_rep_params', n_mm_species)
      end if

      call comms_bcast(pub_root_proc_id,mm_rep_params_exist)

      ! jd: If block is absent, leave it, polemb's internal structure will
      !     remain unallocated, this is fine.
      if(.not. mm_rep_params_exist) return

      call comms_bcast(pub_root_proc_id,n_mm_species)

      allocate(cbuf(n_mm_species),stat=ierr)
      call utils_alloc_check('internal_mm_rep_params','cbuf',ierr)
      allocate(dbuf(2,n_mm_species),stat=ierr)
      call utils_alloc_check('internal_mm_rep_params','dbuf',ierr)

      ! jd: On root read the block and parse it
      if (pub_on_root) then
         do row = 1, n_mm_species
            read(block_data(row),*) dummy_id, dbuf(1,row), dbuf(2,row)
            cbuf(row)(1:10) = trim(dummy_id)
         end do
      end if

      ! jd: Broadcast parsed block contents.
      do row = 1, n_mm_species
         call comms_bcast(pub_root_proc_id,dbuf(:,row))
         call comms_bcast(pub_root_proc_id,cbuf(row))
      end do

      ! jd: Have polemb interpret the data, allocate internal structures and
      !     store the data there.
      call polarisable_embedding_init_mm_rep_params(n_mm_species, cbuf, dbuf)

      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_mm_rep_params','dbuf',ierr)
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_mm_rep_params','cbuf',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Read <mm_rep_params> block'

    end subroutine internal_mm_rep_params

    !------------------------
    !--------- %block_is_soft_sphere_radii
    !------------------------

    ! gab: by Gabriel Bramley, 01/2019.
    subroutine internal_is_soft_sphere_radii

      use comms, only: pub_root_proc_id, comms_bcast
      use constants, only: DP, periodic_table_ss_radii, CRLF
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_is_implicit_solvent, pub_is_dielectric_function, &
           pub_is_soft_sphere_scale
      use utils, only: utils_abort

      implicit none

!     gab: Dummy variables
      real(kind=DP) :: dummy_ss_radius
      character(len=4) :: dummy_spec

!     gab: Local variables
      logical :: is_ss_exists, species_found
      integer :: n_ss_radii, iat, nrows, atomic_no, sat
      real(kind=DP), allocatable :: cbuf(:)
      character(len=4), allocatable :: dbuf(:)

      ! gab: If you're interested in using soft_sphere_radii for something
      ! other than soft spheres please feel free to remove these conditions/add
      ! to them. But for now, these will be in place to prevent errors
      ! occuring in standard calculations.

      if (.not. pub_is_implicit_solvent) return
      if (pub_is_dielectric_function /= 'SOFT_SPHERE') return

      ! gab: Sets up all soft sphere radii with the default Alvarez values,
      ! scaled by the input value.

      if (pub_on_root) then

         is_ss_exists = esdf_block('is_soft_sphere_radii', &
              n_ss_radii)
         write(stdout,'(a)') 'IS: Assigning soft sphere radii for soft&
              & sphere solvation...'

      end if

      do iat = 1,nat
         atomic_no = elements(iat)%atomic_number
         elements(iat)%soft_sphere_radius = pub_is_soft_sphere_scale * &
              periodic_table_ss_radii(atomic_no)
      end do

      call comms_bcast(pub_root_proc_id, is_ss_exists)
      call comms_bcast(pub_root_proc_id, n_ss_radii)

      if (.not. is_ss_exists .or. n_ss_radii == 0) then

         if (pub_on_root) then
            write(stdout,'(a)') 'IS: No soft sphere radii specified. &
                 &Default Alvarez radii used for all atoms.'
         end if

      else

         if (pub_on_root) then

            write(stdout,'(/a/)') 'IS: Reading soft sphere radii from &
                 & is_soft_sphere_radii...'

            if (n_ss_radii .lt. nsp) then
               write(stdout,'(/a)') 'Warning! The number of species in &
                    &is_soft_sphere_radii does not match'//CRLF//'the species&
                    & block. Using default Alvarez radii for unspecified &
                    &species.'
            end if

         end if

         ! gab: Read the inputs into the cbuf array

         allocate(dbuf(n_ss_radii),stat=ierr)
         call utils_alloc_check('internal_soft_sphere_radii','dbuf',ierr)
         allocate(cbuf(n_ss_radii),stat=ierr)
         call utils_alloc_check('internal_soft_sphere_radii','cbuf',ierr)

         if (pub_on_root) then

            do nrows = 1, n_ss_radii
               read(block_data(nrows),*) dbuf(nrows), cbuf(nrows)
            end do

         end if

         ! gab: Assign soft_sphere_radii to all species, regardless of order
         ! in block.

         do iat = 1,n_ss_radii
            call comms_bcast(pub_root_proc_id,cbuf(iat))
            call comms_bcast(pub_root_proc_id,dbuf(iat))
         end do

         do nrows = 1,n_ss_radii

            species_found = .false.
            dummy_spec = dbuf(nrows)
            dummy_ss_radius = cbuf(nrows)

            do iat = 1,nat
               if (dummy_spec == elements(iat)%species_id) then
                  elements(iat)%soft_sphere_radius = dummy_ss_radius
                  species_found = .true.
               end if
            end do

            if (.not. species_found .and. pub_on_root) then
               call utils_abort('Error in internal_is_soft_sphere_radii: &
                    &species specified in the is_soft_sphere_radii block &
                    &mismatch with those in the species block.')
            end if

         end do

         deallocate(cbuf,stat=ierr)
         call utils_dealloc_check('internal_is_soft_sphere_radii','cbuf',ierr)
         deallocate(dbuf,stat=ierr)
         call utils_dealloc_check('internal_is_soft_sphere_radii','dbuf',ierr)

      end if

      ! gab: Print out soft spheres assigned to each species type.

      if (pub_on_root) then
         write (stdout, '(/a)') '----------- SOFT SPHERE RADII -----------'

         do sat = 1,nsp
            do iat = 1,nat

               if (species(sat)%species_id == elements(iat)%species_id) then
                  write (stdout, '(a,f5.3)') species(sat)%species_id, &
                       elements(iat)%soft_sphere_radius
                  exit
               end if

            end do
         end do
         write (stdout, '(a/)') '-----------------------------------------'
      end if

    end subroutine internal_is_soft_sphere_radii

    !-------------------------------------------------------------------------

    !------------------------
    !--------- %block_swri
    !------------------------

    ! jd: by Jacek Dziedzic in 02/2015, largely cloned from internal_species.
    ! JCW: Increased rwbuf from 5 to 10 characters to provide headroom for new
    ! JCW: flags required for invoking new integration scheme
    subroutine internal_swri

      use comms, only: pub_root_proc_id, comms_bcast
      use constants, only: max_rundat_swri_block_entries
      use esdf, only: block_data, esdf_block
      use sw_resolution_of_identity, only: swri_init_module_stage_0
      use utils, only: utils_assert

      implicit none

      ! Parameters
      integer, parameter :: nflags = 10

      ! Local variables
      integer            :: n_entries
      character(len=32)  :: nbuf(max_rundat_swri_block_entries) ! Names of entries
      integer            :: lbuf(max_rundat_swri_block_entries) ! max_l
      integer            :: qbuf(max_rundat_swri_block_entries) ! max_q
      character(len=2)   :: mbuf(max_rundat_swri_block_entries) ! 1 or 2 characters for metric(s)
      integer            :: cibuf(max_rundat_swri_block_entries) ! cheb_intervals
      integer            :: cobuf(max_rundat_swri_block_entries) ! cheb_order
      character(len=nflags) :: rwbuf(max_rundat_swri_block_entries) ! 1 to 10 characters for read/write/assemble/disassemble/etc
      integer            :: row
      character(len=32)  :: dummy_id
      character(len=2)   :: dummy_m
      character(len=5)   :: dummy_rw

      ! ----------------------------------------------------------------------

      ! jd: Take care of padding between n_entires and max_rundat_block_entries
      nbuf(:) = '                                '
      lbuf(:) = 0; qbuf(:) = 0; mbuf(:) = '  '; cibuf(:) = 0; cobuf(:) = 0;
      rwbuf(:) = repeat(' ',nflags)

      ! jd: On root read the block and parse it
      if (pub_on_root) then
         if (esdf_block('swri',n_entries)) then
            call utils_assert(&
                 n_entries <= max_rundat_swri_block_entries, &
                 'Too many SWRIs in %block swri. Current limit is 1 due to &
                 &limit of one blocking scheme for SWs currently supported in &
                 &sparse_initialise and par%num_sw being a scalar.')
            do row=1, n_entries
               read(block_data(row),*) dummy_id, lbuf(row), qbuf(row), &
                    dummy_m, cibuf(row), cobuf(row), dummy_rw
               nbuf(row)(1:32) = trim(dummy_id)
               mbuf(row)(1:2)  = trim(dummy_m)
               rwbuf(row)(1:5) = trim(dummy_rw)
            end do
         else
            n_entries = 0
         end if
      end if

      ! jd: Broadcast parsed block contents
      call comms_bcast(pub_root_proc_id,n_entries)
      call comms_bcast(pub_root_proc_id,lbuf,n_entries)
      call comms_bcast(pub_root_proc_id,qbuf,n_entries)
      call comms_bcast(pub_root_proc_id,cibuf,n_entries)
      call comms_bcast(pub_root_proc_id,cobuf,n_entries)
      do row=1, n_entries
         call comms_bcast(pub_root_proc_id,nbuf(row),32)
         call comms_bcast(pub_root_proc_id,mbuf(row),2)
         call comms_bcast(pub_root_proc_id,rwbuf(row),5)
      end do

      call swri_init_module_stage_0(nbuf, lbuf, qbuf, mbuf, cibuf, cobuf, &
           rwbuf, n_entries)

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Read <swri> block'

    end subroutine internal_swri
!==========================================================================================  START

    !-------------------------------------------------------------------!
    !--------- %internal_eda
    !-------------------------------------------------------------------!
    ! mjsp: EDA fragment atoms list (29/07/2014)                        !
    ! mjsp: EDA fragment loading from disk (21/09/2015)                 !
    !       Modified from internal_species_write_nbo in this module     !
    !-------------------------------------------------------------------!


    subroutine internal_eda()

      use comms,  only: pub_root_proc_id, comms_bcast
      use constants, only: CRLF, EDA_INIT
      use esdf,   only: block_data, esdf_block
      use rundat, only: pub_eda, pub_frag_iatm, pub_frag_charge, &
         pub_eda_split_atoms, pub_frag_atoms, pub_frag_deloc, &
         pub_eda_frag_isol_ct, &
         pub_eda_frag_isol_pol, &
         pub_eda_deltadens, pub_eda_continuation, &
         pub_eda_read_frags, pub_frag_file_prefix, &
         pub_eda_read_super, pub_super_file_prefix, pub_maxit_palser_mano, &
         pub_eda_mode, pub_eda_nodiag, pub_spin_polarised
      use utils, only: utils_abort, utils_assert

      implicit none

      ! mjsp: Local variables
      integer :: nfrag, nfrag_test, frag_nat, frag_charge_int, &
         nsuper, natoms, npairs
      real(kind=DP) :: frag_charge
      character(len=80) :: file_prefix

      ! mjsp: Misc. scalars
      integer ::  row, ii, ierr

      ! mjsp: Checking flags
      logical :: fraglist_atoms_exists, list_files_exists, skip

      ! mjsp: Catch illogical reading of both fragments and supermolecule
      ! in same calculation
      if((pub_eda) .and. (pub_eda_read_super) .and. (pub_eda_read_frags)) &
         call utils_abort('Error in internal_eda(): &
           &cannot read supermolecule and fragments in the same calculation.')

     ! mjsp: continuation calculation and delta density calculation incompatibility
     if (pub_eda_deltadens .and. pub_eda_continuation) &
         call utils_abort('Error in internal_eda(): &
           &cannot calculate electron density difference &
           &with'//CRLF//'EDA continuation calculation. &
           &Electron density differences may be externally '//CRLF//&
           'calculated from electron densities (see utils).')

      ! mjsp: Catch illogical continuation settings
      if((pub_eda) .and. (.not.(pub_eda_read_super)) .and. (pub_eda_continuation)) &
         call utils_abort('Error in internal_eda(): &
           &must specify supermolecule EDA data in a continuation calculation.')
      if((pub_eda) .and. (pub_eda_read_frags) .and. (pub_eda_continuation)) &
         call utils_abort('Error in internal_eda(): &
           &fragment data cannot be read in an EDA continuation calculation.')

      ! mjsp: Spin-polarised case is in development
      if((pub_eda) .and. (pub_spin_polarised) .and. (pub_eda_read_frags)) &
         call utils_abort('Error in internal_eda(): &
           &spin-polarised calculations are in development.')

      ! mjsp: Palser-Manolopoulos canonical purification not coded for
      ! EDA.  I think this would be relatively easy to implement by adding calls to
      ! fragment_matrix_ops routines in palser_mano_mod.F90 - but strong chance
      ! that this might be rather slow to run considering the current
      ! bulkiness of fragment_matrix_ops routines.
      call utils_assert(.not.(pub_eda .and. pub_maxit_palser_mano .ge. 1), &
          'Error in internal_eda(): Please ensure Palser-Manolopoulos'//CRLF//'&
            &canonical purification is not used with SCF-MI by setting'//CRLF//' &
            &maxit_palser_mano = 0.')

      if (pub_eda .and. pub_eda_nodiag /= 3) then
         call utils_abort( &
              "'eda_nodiag' keyword must be set to 3 (default). " // &
              "Other values are disabled due to issues in using " // &
              "diagonalization to construct the antisymmetrized " // &
              "frozen density kernel and possibly later steps.")
      end if

      if((pub_eda) .and. (pub_eda_nodiag .gt. 0) .and. ((pub_eda_frag_isol_pol) .or. &
           (pub_eda_frag_isol_ct))) &
         call utils_abort('Error in internal_eda(): &
              &diagonalisation-free EDA not yet available for'//CRLF//' &
              &fragment-specific EDA calculations')


      ! mjsp: Allocate necessary arrays
      if (allocated(pub_frag_iatm)) then
         deallocate(pub_frag_iatm,stat=ierr)
         call utils_dealloc_check('internal_eda', &
              'pub_frag_iatm',ierr)
      end if

      if (allocated(pub_frag_atoms)) then
         deallocate(pub_frag_atoms,stat=ierr)
         call utils_dealloc_check('internal_eda', &
              'pub_frag_atoms',ierr)
      end if

      if (allocated(pub_frag_deloc)) then
         deallocate(pub_frag_deloc,stat=ierr)
         call utils_dealloc_check('internal_eda', &
              'pub_frag_deloc',ierr)
      end if

      if (allocated(pub_frag_charge)) then
         deallocate(pub_frag_charge,stat=ierr)
         call utils_dealloc_check('internal_eda', &
              'pub_frag_charge',ierr)
      end if

      nfrag = 0
      if (pub_eda) then

         ! mjsp: Initialise global EDA flag
         pub_eda_mode = EDA_INIT

         ! mjsp: Get EDA frag atoms list block if pub_eda = T
         if(pub_on_root) fraglist_atoms_exists &
              = esdf_block('eda_iatm',nfrag)

         call comms_bcast(pub_root_proc_id,fraglist_atoms_exists)
         call comms_bcast(pub_root_proc_id,nfrag)

         if(fraglist_atoms_exists .and. nfrag > 0) then
            allocate(pub_frag_iatm(0:nfrag),stat=ierr)
            call utils_alloc_check('internal_eda', &
                 'pub_frag_iatm',ierr)

            allocate(pub_frag_charge(0:nfrag),stat=ierr)
            call utils_alloc_check('internal_eda', &
                 'pub_frag_charge',ierr)

            pub_frag_iatm(0) = 0 ! idx 0 stores # fragments
            pub_frag_iatm(1:nfrag) = -1 ! mjsp: -1 for error checking
            pub_frag_charge(0:nfrag) = 0.0_DP ! idx 0 stores total charge
                                              ! on supermolecule
         else
            call utils_abort('Error in internal_eda() reading &
                 &block eda_iatm.')
         end if

         if(pub_on_root) then ! mjsp: Read block

            ! mjsp: Read eda_iatm
            do row=1,nfrag

               read(block_data(row),*) frag_nat, frag_charge_int
               frag_charge = real(frag_charge_int,kind=DP)  ! int->real

               if(frag_nat < 0) then
                  call utils_abort('Error in internal_eda(): &
                       &frag_nat < 0')
               end if
               if(frag_charge < -50) then
                  write(stdout,'(a,I4)') 'WARNING in internal_eda(): &
                       &frag_charge < -50 for fragment ID:', row
               end if
               if(frag_charge > 50) then
                  write(stdout,'(a,I4)') 'WARNING in internal_eda(): &
                       &frag_charge > 50 for fragment ID:', row
               end if

               ! mjsp: total number of fragments:
               pub_frag_iatm(0) = &
                    pub_frag_iatm(0) + 1

               ! mjsp: atoms in each fragment:
               pub_frag_iatm(pub_frag_iatm(0)) = &
                    frag_nat

               ! mjsp: charge on the full supermolecule:
               pub_frag_charge(0) = pub_frag_charge(0) + frag_charge

               ! mjsp: charge on each fragment:
               pub_frag_charge(row) = frag_charge

            end do ! END do row=1,nfrag

         end if ! END if(pub_on_root)

         ! mjsp: Get EDA frag dkn and NGWFs list block if pub_eda = T and
         ! pub_eda_read_frags = T
         if (pub_eda_read_frags .eqv. .true.) then

           if(pub_on_root) list_files_exists &
                = esdf_block('eda_frags',nfrag_test)

           call comms_bcast(pub_root_proc_id,list_files_exists)
           call comms_bcast(pub_root_proc_id,nfrag_test)

           if(list_files_exists .and. nfrag_test > 0) then
              allocate(pub_frag_file_prefix(1:nfrag_test),stat=ierr)
              call utils_alloc_check('internal_eda', &
                   'pub_frag_charge',ierr)

              !pub_frag_file_prefix(0:nfrag_test) = ''

           else
              call utils_abort('Error in internal_eda() reading &
                   &block eda_frags.')
           end if

           ! ensure number of lines in this block equals number of fragments
           ! (i.e. that the eda_iatm and eda_frags blocks agree)
           call utils_assert(nfrag .eq. nfrag_test, 'Error in internal_eda(): &
                &number of fragments given in eda_frags &
                &does not equal'//CRLF//'number of fragments given &
                &in eda_iatm')

           if(pub_on_root) then ! mjsp: Read block

             ! mjsp: Read eda_frags
             do row=1,nfrag

                read(block_data(row),*) file_prefix

                pub_frag_file_prefix(row) = file_prefix

                if(len(pub_frag_file_prefix(row)) < 0) then
                   call utils_abort('Error in internal_eda(): &
                        &fragment filename prefix must not be blank.')
                end if

             end do ! END do row=1,nfrag

           end if ! END if(pub_on_root)

           do row=1,nfrag
             call comms_bcast(pub_root_proc_id,pub_frag_file_prefix(row))
           end do

         end if ! END if(pub_eda_read_frags)


         ! mjsp: Get EDA supermolecule dkn and NGWFs list block if pub_eda = T and
         ! pub_eda_read_super = T
         if (pub_eda_read_super .eqv. .true.) then

           if(pub_on_root) list_files_exists = esdf_block('eda_super',nsuper)

           call comms_bcast(pub_root_proc_id,list_files_exists)
           call comms_bcast(pub_root_proc_id,nsuper)

           if(.not.(list_files_exists .and. nsuper == 1)) then
              call utils_abort('Error in internal_eda(): EDA_SUPER block &
                  &must contain only one supermolecule filename &
                  &prefix')
           end if

           if(pub_on_root) then ! mjsp: Read block

             ! mjsp: Read eda_super
             read(block_data(1),*) pub_super_file_prefix

             if(len(pub_super_file_prefix) < 0) then
                call utils_abort('Error in internal_eda(): &
                     &supermolecule filename prefix must not be blank.')
             end if

           end if ! END if(pub_on_root)

           call comms_bcast(pub_root_proc_id,pub_super_file_prefix)

         end if ! END if(pub_eda_read_super)


         ! mjsp: Get fragment pairs to calculate delocalisations between
         if (pub_eda_frag_isol_ct) then

            ! mjsp: Fragment delocalisations are currently under development and should not
            ! be used.
            call utils_abort('Error in internal_eda(): &
                 &Interfragment delocalisations are currently' &
                 //CRLF//'under development.')

            if(pub_on_root) list_files_exists = esdf_block('eda_deloc',npairs)

            call comms_bcast(pub_root_proc_id,list_files_exists)
            call comms_bcast(pub_root_proc_id,npairs)

            if(.not.(list_files_exists)) then
               call utils_abort('Error in internal_eda(): eda_frag_isol_ct &
                  &parameter set to TRUE: &
                  &expected EDA_DELOC block but no block found! ')
            end if

            ! mjsp: allocate list of fragment pairs (frag1#, frag2#)
            allocate(pub_frag_deloc(2,npairs),stat=ierr)
            call utils_alloc_check('internal_eda', &
             'pub_frag_deloc',ierr)

            if(pub_on_root) then ! mjsp: Read block

              ! mjsp: Read eda_iatm
              do row=1,npairs

                 ! mjsp: Read fragment pairs
                 read(block_data(row),*) pub_frag_deloc(1,row), pub_frag_deloc(2,row)

                 ! mjsp: frag2 must be greater than frag1
                 call utils_assert(pub_frag_deloc(1,row) .lt. pub_frag_deloc(2,row), &
                   'Error in internal_eda(): fragment index #2 must be greater than &
                   &fragment index #1 in frag_deloc block')

              end do

            end if ! END if(pub_on_root)

            call comms_bcast(pub_root_proc_id,pub_frag_deloc)

         end if ! END if (pub_eda_frag_isol_ct)


         ! mjsp: Get EDA atoms to split (block) if pub_eda = T and
         ! pub_eda_split_atoms = T
         if (pub_eda_split_atoms .eqv. .true.) then

            ! mjsp: Fragment splitting is in development and should not be used
            ! PLACEHOLDER
            call utils_abort('Error in internal_eda(): &
                 &Atom splitting is currently under development.')

            if(pub_on_root) list_files_exists = esdf_block('eda_frag_atoms',natoms)

            call comms_bcast(pub_root_proc_id,list_files_exists)
            call comms_bcast(pub_root_proc_id,natoms)

            if(.not.(list_files_exists)) then
               call utils_abort('Error in internal_eda(): eda_split_atoms parameter set to TRUE: &
                  &expected EDA_FRAG_ATOMS block but no block found! ')
            end if

            ! mjsp: allocate list of atoms (frag#, atom#, NGWF#)
            allocate(pub_frag_atoms(3,natoms),stat=ierr)
            call utils_alloc_check('internal_eda', &
             'pub_frag_atoms',ierr)

            if(pub_on_root) then ! mjsp: Read block

              ! mjsp: Read eda_iatm
              do row=1,natoms

                 ! mjsp: Read atoms (frag#, atom#, NGWF#)
                 read(block_data(row),*) pub_frag_atoms(1,row), pub_frag_atoms(2,row), pub_frag_atoms(3,row)

              end do

            end if ! END if(pub_on_root)

            call comms_bcast(pub_root_proc_id,pub_frag_atoms)

         end if ! END if(pub_eda_read_super)

         call comms_bcast(pub_root_proc_id,pub_frag_iatm)
         call comms_bcast(pub_root_proc_id,pub_frag_charge)

      end if ! END if(pub_eda)

    end subroutine internal_eda

!==========================================================================================  END

    !----------------------------
    !--------- %block species_pot
    !----------------------------

    subroutine internal_species_pot()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_pseudo_path
      use utils, only: utils_abort

      implicit none

      !-------------
      ! ars: dummies
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      character(len=64) :: dummy_pseudoname
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      !-------------
      ! ars: integers
      integer :: row, ii
      ! ars: logical
      logical :: element_found, species_found


      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_species_pot','cbuf',ierr)

      ! aam: Root proc only: read information from
      !      <species_pot> block of input file
      if (pub_on_root) then
         if (esdf_block('species_pot',dummy_nsp)) then
            if (nsp == dummy_nsp)  then
               ! aam: check for multiply defined species
               call internal_check_species(nsp,'<species_pot>')
               do row=1,nsp
                  read(block_data(row),*) dummy_id, dummy_pseudoname
                  cbuf(row)(1:4) = dummy_id
                  cbuf(row)(5:68) = dummy_pseudoname
               end do
            else
               call utils_abort('Error in internal_species_pot(): &
                    &mismatching numbers of species in <species_pot> block of &
                    &input file')
            end if
         else
            call utils_abort('Error in internal_species_pot(): &
                 &<species_pot> block not found in input file')
         end if
      end if

      ! pdh: broadcast species information
      do row=1,nsp
         call comms_bcast(pub_root_proc_id,cbuf(row))
      end do

      ! extract element information from species
      do row=1,nsp
         element_found = .false. ! aam: error check flag
         dummy_id = cbuf(row)(1:4)
         dummy_pseudoname = cbuf(row)(5:68)

         ! set species version
         species_found = .false.
         do ii=1,nsp
            if (species(ii)%species_id == dummy_id) then
               species(ii)%pseudo_name = dummy_pseudoname
               species_found = .true.
            end if
         end do

         ! set elements version
         do ii=1,nat
            if (elements(ii)%species_id == dummy_id) then
               element_found = .true.
            end if
         end do
         if (.not. element_found) then
            call utils_abort('Error in internal_species_pot: &
                 &mismatching species in <species_pot> and <positions_abs> &
                 &blocks of input file')
         end if
      end do


      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_species_pot','cbuf',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_pot> block'

    end subroutine internal_species_pot


    !----------------------------
    !--------- %block species_core_wf
    !----------------------------

    subroutine internal_species_core_wf()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use utils, only: utils_abort

      implicit none

      ! Local Variables
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      character(len=64) :: dummy_corewfname
      character(len=68), allocatable :: cbuf(:)
      integer :: row, ii
      logical :: element_found, species_found

      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_species_core_wf','cbuf',ierr)

      ! ndmh: Root proc only: read information from
      !       <species_core_wf> block of input file
      if (pub_on_root) then
         if (esdf_block('species_core_wf',dummy_nsp)) then
            if (nsp == dummy_nsp)  then
               ! aam: check for multiply defined species
               call internal_check_species(nsp,'<species_core_wf>')
               do row=1,nsp
                  read(block_data(row),*) dummy_id, dummy_corewfname
                  cbuf(row)(1:4) = dummy_id
                  cbuf(row)(5:68) = dummy_corewfname
               end do
            else
               call utils_abort('Error in internal_species_core_wf: &
                    &mismatching numbers of species in <species_core_wf> &
                    &block of input file.')
            end if
         else
            !write(stdout,'(/a/)',advance='no') &
            !     '  <species_core_wf> block not found'
            do row=1,nsp
               cbuf(row)(1:4) = all_species_id(row)(1:4)
               cbuf(row)(5:68) = 'NONE'//repeat(' ',60)
            end do
         end if
      end if

      ! ndmh: broadcast species information
      do row=1,nsp
         call comms_bcast(pub_root_proc_id,cbuf(row))
      end do

      ! extract element information from species
      do row=1,nsp
         element_found = .false. ! aam: error check flag
         dummy_id = cbuf(row)(1:4)
         dummy_corewfname = cbuf(row)(5:68)

         ! set species version
         species_found = .false.
         do ii=1,nsp
            if (species(ii)%species_id == dummy_id) then
               species(ii)%core_wf_name = dummy_corewfname
               species_found = .true.
            end if
         end do

         ! set elements version
         do ii=1,nat
            if (elements(ii)%species_id == dummy_id) then
               element_found = .true.
            end if
         end do
         if (.not. element_found) then
            call utils_abort('Error in internal_species_core_wf: &
                 &mismatching species in <species_core_wf> and <positions_abs> &
                 &blocks of input file.')
         end if
      end do


      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_species_core_wf','cbuf',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_core_wf> block'

    end subroutine internal_species_core_wf


    !-----------------------------
    !--------- %block species_cond
    !-----------------------------

    subroutine internal_species_cond()

      use comms, only: pub_root_proc_id, comms_bcast
      use constants, only: CRLF, kernel_cutoff_sanity_threshold
      use esdf, only: block_data, esdf_block, esdf_convfac
      use rundat, only: pub_ngwf_halo, pub_cond_calculate_any_task, &
           cond_kernel_cutoff
      use utils, only: utils_abort, utils_assert, utils_flushed_string_output

      implicit none

      !-------------
      ! ars: dummies
      character(len=80) :: dummy_len_unit
      character(len=4)  :: dummy_id
      character(len=2)  :: dummy_symb
      integer :: dummy_nsp, dummy_nfunc, dummy_Z
      real(kind=dp) :: dummy_rad
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      integer, allocatable :: ibuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      character(len=256) :: message
      character(len=10) :: valuestr, valuestr2
      !-------------
      ! ars: integers
      integer :: ishift, row, ii
      integer :: ierr
      ! ars: logical
      logical :: element_found, species_cond_block_found
      ! ars: real
      real(kind=dp) :: lenconfac
      real(kind=DP) :: max_ngwf_rad_cond
      character(len=*), parameter :: myself = 'internal_species_cond'

      ! -----------------------------------------------------------------------

      ! pdh: allocate buffers
      ! ddor: Altered for DFT+U
      allocate(dbuf(max(nat*4,9)),stat=ierr)
      call utils_alloc_check(myself,'dbuf',ierr)
      allocate(ibuf(nat*4),stat=ierr)
      call utils_alloc_check(myself,'ibuf',ierr)
      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check(myself,'cbuf',ierr)


      ! ndmh: Root proc only: read information from <species_cond> block of input
      ! ndmh: file for conduction band NGWFs
      ! lr408: If this isn't a conduction calculation, don't bother checking
      ! lr408: for conduction species block
      if (pub_on_root) then
         lenconfac=1.0_dp
         ishift=0
         species_cond_block_found = esdf_block('species_cond',dummy_nsp)
         if (species_cond_block_found) then
            ! aam: 3/6/09 added ability to read optional  unit string
            !      ("bohr" or "ang") in first line of species block
            !      (default is bohr)
            if (dummy_nsp == nsp+1) then
               read(block_data(1),'(a)') dummy_len_unit
               lenconfac = esdf_convfac(dummy_len_unit,'bohr')
               ishift=1
            endif
            call utils_assert( (dummy_nsp.eq.nsp) .or. (dummy_nsp.eq.nsp+1), &
                 'Error in internal_species_cond: mismatching numbers of &
                 &species in <species_cond> specification.'//CRLF//&
                 'Perhaps your positions block, or input .xyz file (if &
                 &using one) contains species that are not listed in the &
                 &species_cond block.')

            ! aam: check for multiply defined species
            call internal_check_species(nsp,'<species_cond>') !ddor
            do row=1,nsp
               read(block_data(row+ishift),*,iostat=ierr) dummy_id,dummy_symb, &
                    ibuf(row*2-1:row*2),dbuf(row)
               ! jd: Catch the corner case where an extra species in the
               !     positions block is initially missed because there's a
               !     unit specification in the species_cond block.
               !     If unchecked, this leads to a crash with an end-of-file
               !     error reported by the RTL.
               call utils_assert(ierr==0, 'Error in '//myself//': &
                    &detected reading the <species_cond> specification. '//&
                    CRLF//&
                    'Perhaps your positions block, or input .xyz file (if &
                    &using one) contains species that are not listed in the &
                    &species block.')
               cbuf(row)(1:4) = dummy_id
               cbuf(row)(5:6) = dummy_symb
            end do
            dbuf(1:nsp)=lenconfac*dbuf(1:nsp)
         else
            call utils_assert(.not. pub_cond_calculate_any_task, &
                 'Error in '//myself//': conduction calculation speci&
                 &fied in task, but no species_cond block has been specified.')
         end if
      end if

      call comms_bcast(pub_root_proc_id,species_cond_block_found)

      if (species_cond_block_found) then

         ! ndmh: broadcast species information
         call comms_bcast(pub_root_proc_id,ibuf,nsp*2)
         call comms_bcast(pub_root_proc_id,dbuf,nsp)
         do row=1,nsp
            call comms_bcast(pub_root_proc_id,cbuf(row),6)
         end do

         ! ndmh: extract element information from species
         element_found = .false.
         do row=1,nsp
            element_found = .false.
            dummy_id = cbuf(row)(1:4)
            dummy_symb = cbuf(row)(5:6)
            dummy_Z = ibuf(row*2-1)
            dummy_nfunc = ibuf(row*2)
            dummy_rad = dbuf(row)
            species(row)%nfunctions_cond = dummy_nfunc
            do ii=1,nat
               if (elements(ii)%species_id == dummy_id) then
                  ! ndmh: check Z and symbol match
                  if (elements(ii)%symbol/=dummy_symb) then
                     call utils_abort('Error in '//myself//': &
                          &mismatching element symbols in <species_cond> and &
                          &<species> blocks.')
                  end if
                  if (elements(ii)%atomic_number /= dummy_Z) then
                     call utils_abort('Error in '//myself//': &
                          &mismatching atomic numbers in <species_cond> and &
                          &<species> blocks.')
                  end if
                  ! ndmh: copy info to elements array
                  elements(ii)%nfunctions_cond = dummy_nfunc
                  elements(ii)%radius_cond     = dummy_rad
                  element_found = .true.
               end if
            end do
            if (.not. element_found) then
               call utils_abort('Error in '//myself//': &
                    &mismatching species in <species_cond> and <positions_abs> &
                    &blocks of input file.')
            end if
         end do

         ! cks: 21/04/2006 Warning! If there is a halo, the element radius
         ! cks: is the NGWF radius plus the halo length and so it is not
         ! cks: equal to the NGWF-sphere radius (which does not include the halo)
         if (pub_ngwf_halo > 0.0_DP) then
            ! cks: increase element radius by halo length
            do ii=1,nat
               elements(ii)%radius_cond = elements(ii)%radius_cond + pub_ngwf_halo
            end do
         end if

         ! jd: Check against unreasonably low kernel cutoffs
         max_ngwf_rad_cond = maxval(elements(:)%radius_cond)
         if(cond_kernel_cutoff <= kernel_cutoff_sanity_threshold) then
            write(valuestr,'(f10.3)') cond_kernel_cutoff
            write(message,'(a,a,a)') &
                 'WARNING: The conduction kernel cutoff you specified (', &
                 trim(adjustl(valuestr)), &
                 ' bohr) seems unreasonably low.'
            call utils_flushed_string_output(message//CRLF)
         end if
         if(cond_kernel_cutoff <= 2.0_DP*max_ngwf_rad_cond) then
            write(valuestr,'(f10.3)') cond_kernel_cutoff
            write(valuestr2,'(f10.3)') max_ngwf_rad_cond
            write(message,'(a,a,a,a,a)') &
                 'The conduction kernel cutoff (you specified ', &
                 trim(adjustl(valuestr)), &
                 ' bohr) must be larger than twice the maximum conduction NGWF&
                 & radius (which you specified to be ', &
                 trim(adjustl(valuestr2)),' bohr).'
            call utils_abort(trim(message))
         end if

         if (pub_debug_on_root) write(stdout,'(a)') &
              'DEBUG: Read <species_cond> block'

      else
         ! ndmh: no conduction NGWF species block is present, so set all related
         ! ndmh: entries to zero
         species(:)%nfunctions_cond = 0
         do ii=1,nat
            elements(ii)%radius_cond = 0.0_DP
            elements(ii)%nfunctions_cond = 0
         end do
      end if

      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check(myself,'cbuf',ierr)
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check(myself,'dbuf',ierr)
      deallocate(ibuf,stat=ierr)
      call utils_dealloc_check(myself,'ibuf',ierr)

    end subroutine internal_species_cond


    !-----------------------------
    !--------- %block species_aux
    !-----------------------------

    subroutine internal_species_aux()

      use comms, only: pub_root_proc_id, comms_bcast
      use constants, only: CRLF
      use esdf, only: block_data, esdf_block, esdf_convfac
      use rundat, only: pub_use_aux_ngwfs
      use utils, only: utils_abort, utils_assert

      implicit none

      !-------------
      ! ars: dummies
      character(len=80) :: dummy_len_unit
      character(len=4)  :: dummy_id
      character(len=2)  :: dummy_symb
      integer :: dummy_nsp, dummy_nfunc, dummy_Z
      integer :: ierr
      real(kind=dp) :: dummy_rad
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      integer, allocatable :: ibuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: integers
      integer :: ishift, row, ii
      ! ars: logical
      logical :: element_found
      ! ars: real
      real(kind=dp) :: lenconfac
      character(len=*), parameter :: myself = 'internal_species_aux'

      ! -----------------------------------------------------------------------

      if (pub_on_root) then
         ! ndmh: activate auxiliary NGWFs if block is present
         pub_use_aux_ngwfs = esdf_block('species_aux',dummy_nsp)
      end if
      call comms_bcast(pub_root_proc_id,pub_use_aux_ngwfs)

      if (pub_use_aux_ngwfs) then

         ! pdh: allocate buffers
         allocate(dbuf(max(nat*4,9)),stat=ierr)
         call utils_alloc_check(myself,'dbuf',ierr)
         allocate(ibuf(nat*4),stat=ierr)
         call utils_alloc_check(myself,'ibuf',ierr)
         allocate(cbuf(nsp),stat=ierr)
         call utils_alloc_check(myself,'cbuf',ierr)


         ! ndmh: Root proc only: read information from <species_aux> block of input
         ! ndmh: file for auxiliary NGWFs
         if (pub_on_root) then
            lenconfac=1.0_dp
            ishift=0
            if (esdf_block('species_aux',dummy_nsp)) then
               ! aam: 3/6/09 added ability to read optional  unit string
               !      ("bohr" or "ang") in first line of species block
               !      (default is bohr)
               if (dummy_nsp == nsp+1) then
                  read(block_data(1),'(a)') dummy_len_unit
                  lenconfac = esdf_convfac(dummy_len_unit,'bohr')
                  ishift=1
               endif
               if ( (dummy_nsp.ne.nsp) .and. (dummy_nsp.ne.nsp+1) ) then
                  call utils_abort('Error in internal_species_aux(): &
                       &mismatching numbers of species in <species_aux> &
                       &specification.'//CRLF//'Perhaps your positions block, &
                       &or input .xyz file (if using one) contains species &
                       &that are not listed in the species_aux block.')
               endif
               ! aam: check for multiply defined species
               call internal_check_species(nsp,'<species_aux>') !ddor
               do row=1,nsp
                  read(block_data(row+ishift),*,iostat=ierr) dummy_id, &
                       dummy_symb, ibuf(row*2-1:row*2),dbuf(row)
                  ! jd: Catch the corner case where an extra species in the
                  !     positions block is initially missed because there's a
                  !     unit specification in the species_aux block.
                  !     If unchecked, this leads to a crash with an end-of-file
                  !     error reported by the RTL.
                  call utils_assert(ierr==0, 'Error in '//myself//': &
                       &detected reading the <species_aux> specification. '//&
                       CRLF//&
                       'Perhaps your positions block, or input .xyz file (if &
                       &using one) contains species that are not listed in the &
                       &species block.')
                  cbuf(row)(1:4) = dummy_id
                  cbuf(row)(5:6) = dummy_symb
               end do
               dbuf(1:nsp)=lenconfac*dbuf(1:nsp)
            else
               call utils_abort('Error in '//myself//': &
                    &auxiliary NGWFs required due to task, but no species_aux &
                    &block has been specified.')
            end if
         end if

         ! ndmh: broadcast species information
         call comms_bcast(pub_root_proc_id,ibuf,nsp*2)
         call comms_bcast(pub_root_proc_id,dbuf,nsp)
         do row=1,nsp
            call comms_bcast(pub_root_proc_id,cbuf(row),6)
         end do

         ! ndmh: extract element information from species
         element_found = .false.
         do row=1,nsp
            element_found = .false.
            dummy_id = cbuf(row)(1:4)
            dummy_symb = cbuf(row)(5:6)
            dummy_Z = ibuf(row*2-1)
            dummy_nfunc = ibuf(row*2)
            dummy_rad = dbuf(row)
            species(row)%nfunctions_aux = dummy_nfunc
            do ii=1,nat
               if (elements(ii)%species_id == dummy_id) then
                  ! ndmh: check Z and symbol match
                  if (elements(ii)%symbol/=dummy_symb) then
                     call utils_abort('Error in '//myself//': &
                          &mismatching element symbols in <species_aux> &
                          &and <species> blocks.')
                  end if
                  if (elements(ii)%atomic_number /= dummy_Z) then
                     call utils_abort('Error in '//myself//': &
                          &mismatching atomic numbers in <species_aux> &
                          &and <species> blocks.')
                  end if
                  ! ndmh: copy info to elements array
                  elements(ii)%nfunctions_aux = dummy_nfunc
                  element_found = .true.
               end if
            end do
            if (.not. element_found) then
               call utils_abort('Error in '//myself//': &
                    &mismatching species in <species_aux> &
                    &and <species> blocks.')
            end if
         end do

         ! pdh: deallocate buffers
         deallocate(cbuf,stat=ierr)
         call utils_dealloc_check(myself,'cbuf',ierr)
         deallocate(dbuf,stat=ierr)
         call utils_dealloc_check(myself,'dbuf',ierr)
         deallocate(ibuf,stat=ierr)
         call utils_dealloc_check(myself,'ibuf',ierr)

         if (pub_debug_on_root) write(stdout,'(a)') &
              'DEBUG: Read <species_aux> block'

      else
         ! ndmh: no auxiliary NGWF species block is present, so set all related
         ! ndmh: entries to zero/false
         pub_use_aux_ngwfs = .false.
         species(:)%nfunctions_aux = 0
         do ii=1,nat
            elements(ii)%nfunctions_aux = 0
         end do
      end if

    end subroutine internal_species_aux





    !-----------------------------------
    !--------- %block species_atomic_set
    !-----------------------------------

    subroutine internal_species_atomic_set(which_type)

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_initial_dens_realspace
      use utils, only: utils_abort, utils_assert

      implicit none

      ! Arguments
      character(*), intent(in) :: which_type

      !-------------
      ! ars: dummies
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      character(len=128) :: dummy_ngwf_set
      !-------------
      ! ars: buffers
      character(len=132), allocatable :: cbuf(:)
      !-------------
      ! ars: integers
      integer :: row, ii
      ! ars: logical
      logical :: element_found, species_found
      logical :: copy_val_ngwfs

      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_species_atomic_set','cbuf',ierr)

      ! aam: Root proc only: read information from
      !      <species_atomic_set> block of input file
      if (pub_on_root) then
         if (esdf_block('species_atomic_set'//trim(which_type),dummy_nsp)) then
            if (nsp == dummy_nsp)  then
               ! aam: check for multiply defined species
               call internal_check_species(nsp,'<species_atomic_set'// &
                    trim(which_type)//'>')
               do row=1,nsp
                  read(block_data(row),*) dummy_id,dummy_ngwf_set
                  cbuf(row)(1:4) = dummy_id
                  cbuf(row)(5:132) = dummy_ngwf_set
               end do
            else
               call utils_abort('Error in internal_species_atomic_set(): &
                    &mismatching numbers of species in <species_atomic_set'//&
                    trim(which_type)//'> block of input file.')
            end if
            copy_val_ngwfs = .false.
         else
            copy_val_ngwfs = .true.
            if (trim(which_type)=='') then
               write(stdout,'(/a/)',advance='no') &
                    '  <species_atomic_set> block not found:&
                    & NGWF initialisation set to SOLVE'
               ! cks: Use automatic NGWF initialisation if species_atomic_set
               ! cks: is missing
               do row=1,nsp
                  cbuf(row)(1:4) = all_species_id(row)(1:4)
#ifdef OLD_DEFAULTS
                  cbuf(row)(5:132) = 'AUTO'//repeat(' ',124)
#else
                  cbuf(row)(5:132) = 'SOLVE'//repeat(' ',123)
#endif
               end do
               copy_val_ngwfs = .false.
            end if
         end if
      end if

      ! pdh: broadcast species information
      call comms_bcast(pub_root_proc_id,copy_val_ngwfs)
      do row=1,nsp
         call comms_bcast(pub_root_proc_id,cbuf(row))
      end do

      do row=1,nsp
         element_found = .false. ! aam: error check flag
         ! ndmh: for cond/aux NGWFs, if we are copying the NGWF sets from
         ! ndmh: the valence NGWF sets, do so now and exit
         if (copy_val_ngwfs) then
            ! set species version
            if (trim(which_type)=='_cond') then
               species(row)%cond_ngwf_set = species(row)%ngwf_set
            else if (trim(which_type)=='_aux') then
               species(row)%aux_ngwf_set = species(row)%ngwf_set
            end if
         else
            dummy_id = cbuf(row)(1:4)
            dummy_ngwf_set = cbuf(row)(5:132)

            ! set species version
            species_found = .false.
            do ii=1,nsp
               if (species(ii)%species_id == dummy_id) then
                  if (trim(which_type)=='') then
                     species(ii)%ngwf_set = dummy_ngwf_set
                  else if (trim(which_type)=='_cond') then
                     species(ii)%cond_ngwf_set = dummy_ngwf_set
                  else if (trim(which_type)=='_aux') then
                     species(ii)%aux_ngwf_set = dummy_ngwf_set
                  else
                     write(stdout,'(/a)') 'Error in internal_species_&
                          &atomic_set: unexpected type of NGWF set'
                  end if
                  species_found = .true.
               end if
            end do
            if (.not. species_found) then
               call utils_abort('Error in internal_species_atomic_set():&
                    & mismatching species in <species_atomic_set> and &
                    &<positions_abs> blocks of input file.')
            end if

            ! JCW: Abort for "AUTO" type initialization, if
            ! JCW: pub_initial_dens_realspace keyword is not F, since the
            ! JCW: AUTO initialization of NGWFs (via ngwfs_generate_sto3g)
            ! JCW: is not compatible with realspace density initialization
            ! JCW: (via density_init_guess_real), and requires reciprocal
            ! JCW: space density initialization
            if (any(species(1:nsp)%ngwf_set ==  'AUTO').or.&
                any(species(1:nsp)%cond_ngwf_set == 'AUTO').or.&
                any(species(1:nsp)%aux_ngwf_set == 'AUTO')) then
               call utils_assert(.not.pub_initial_dens_realspace,&
                    'Error in internal_species_atomic_set(): &
                    &AUTO-type initialization requires the keyword &
                    &initial_dens_realspace is set to F. Please modify your &
                    &input file and re-run the calculation.')
            end if

            do ii=1,nat
               if (elements(ii)%species_id == dummy_id) then
                  element_found = .true.
               end if
            end do
            if (.not. element_found) then
               call utils_abort('Error in internal_species_atomic_set():&
                    & mismatching species in <species_atomic_set> and &
                    &<positions_abs> blocks of input file.')
            end if
         end if
      end do

      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_species_atomic_set','cbuf',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_atomic_set> block'

    end subroutine internal_species_atomic_set






    !------------------------------------
    !--------- %block species_constraints
    !------------------------------------

    subroutine internal_species_constraints()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use utils, only: utils_abort, utils_assert

      implicit none


      !-------------
      ! ars: dummies
      character(len=5)  :: dummy_contype
      real(kind=dp)     :: dummy_convec(3)
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: integers
      integer :: row, ii, ierr
      ! ars: logical
      logical :: element_found



      ! pdh: allocate buffers
      ! ddor: Altered for DFT+U
      allocate(dbuf(max(nat*4,9)),stat=ierr)
      call utils_alloc_check('internal_species_constraints','dbuf',ierr)
      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_species_constraints','cbuf',ierr)

      ! aam: Root proc only: read information from
      !      <species_constraints> block of input file
      if (pub_on_root) then
         if (esdf_block('species_constraints',dummy_nsp)) then
            call utils_assert(dummy_nsp <= nsp, 'Error in internal_species_&
                 &constraints(): too many species in <species_constraints> &
                 &block of input file.')

            ! aam: check for multiply defined species
            call internal_check_species(dummy_nsp,'<species_constraints>')
            do row=1,dummy_nsp ! aam: dummy_nsp is not necessarily equal to nsp
               if (index(block_data(row),'FIXED') /= 0 .or. &
                    index(block_data(row),'NONE') /= 0) then
                  ! aam: only species id and constraint type required here
                  read(block_data(row),*) dummy_id, dummy_contype
                  dbuf(row*3-2:row*3) = 0.0_DP
               else if (index(block_data(row),'LINE') /= 0 .or. &
                    index(block_data(row),'PLANE') /= 0) then
                  ! aam: also require the constraint vector here
                  read(block_data(row),*) dummy_id, dummy_contype, &
                       dbuf(row*3-2:row*3)
               else
                  call utils_abort('Error in internal_species_constraints(): &
                       &unknown constraint identifier in <species_constraints> &
                       &block of input file.')
               end if
               cbuf(row)(1:4) = dummy_id
               cbuf(row)(5:9) = dummy_contype
            end do
         end if
      end if

      ! pdh: broadcast
      call comms_bcast(pub_root_proc_id,dummy_nsp)
      do row=1,dummy_nsp
         call comms_bcast(pub_root_proc_id,cbuf(row),9)
      end do
      call comms_bcast(pub_root_proc_id,dbuf,dummy_nsp*3)

      ! aam: initialise to default values
      do ii=1,nat
         elements(ii)%ion_constraint_type = 'NONE'
         elements(ii)%ion_constraint      = 0.0_DP
      end do

      ! pdh: extract constraint information
      do row=1,dummy_nsp
         element_found = .false. ! aam: error check flag
         dummy_id = cbuf(row)(1:4)
         dummy_contype = cbuf(row)(5:9)
         dummy_convec = dbuf(row*3-2:row*3)
         do ii=1,nat
            if (elements(ii)%species_id == dummy_id) then
               elements(ii)%ion_constraint_type = dummy_contype
               elements(ii)%ion_constraint      = dummy_convec
               element_found = .true.
            end if
         end do
         if (.not. element_found) then
            call utils_abort('Error in internal_species_constraints(): &
                 &mismatching species in <species_constraints> and &
                 &<positions_abs> blocks of input file.')
         end if
      end do

      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_species_constraints','cbuf',ierr)
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_species_constraints','dbuf',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_constraints> block'

    end subroutine internal_species_constraints



    !----------------------------------
    !--------- %block_species_ngwf_plot
    !----------------------------------

    subroutine internal_species_ngwf_plot()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use utils, only: utils_abort, utils_assert

      implicit none

      !-------------
      ! ars: dummies
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      !-------------
      ! ars: integers
      integer ::  row, ii, ierr
      ! ars: logical
      logical :: element_found



      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_species_ngwf_plot','cbuf',ierr)

      ! cks: Root proc only: read information from
      !      <species_ngwf_plot> block of input file
      if (pub_on_root) then
         if (esdf_block('species_ngwf_plot',dummy_nsp)) then
            call utils_assert(dummy_nsp <= nsp, &
                 'Error in internal_species_ngwf_plot(): &
                 &too many species in <species_ngwf_plot> block of input file.')
            ! aam: check for multiply defined species
            call internal_check_species(dummy_nsp,'<species_ngwf_plot>')
            do row=1,dummy_nsp ! aam: dummy_nsp is not necessarily equal to nsp
               read(block_data(row),*) dummy_id
               cbuf(row)(1:4) = dummy_id
            end do
         end if
      end if

      ! pdh: broadcast
      call comms_bcast(pub_root_proc_id,dummy_nsp)
      do row=1,dummy_nsp
         call comms_bcast(pub_root_proc_id,cbuf(row),4)
      end do

      ! cks: initialise to default values
      do ii=1,nat
         elements(ii)%ngwf_plot = .false.
      end do

      ! pdh: extract information
      do row=1,dummy_nsp
         element_found = .false. ! aam: error check flag
         dummy_id = cbuf(row)(1:4)
         do ii=1,nat
            if (elements(ii)%species_id == dummy_id) then
               elements(ii)%ngwf_plot = .true.
               element_found = .true.
            end if
         end do
         if (.not. element_found) then
            call utils_abort('Error in internal_species_ngwf_plot: &
                 &mismatching species in <species_ngwf_plot> and &
                 &<positions_abs> blocks of input file.')
         endif
      end do

      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_species_ngwf_plot','cbuf',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_ngwf_plot> block'

    end subroutine internal_species_ngwf_plot


    !-----------------------------------------
    !--------- %block species_locdipole_groups
    !-----------------------------------------

    subroutine internal_species_locdipole_groups()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use utils, only: utils_abort, utils_assert

      implicit none

      ! vv: dummy variables
      character(len=4) :: dummy_id
      integer :: dummy_nsp
      ! vv: dummy buffers
      character(len=68), allocatable :: cbuf(:)
      ! vv: integers and locals
      integer :: row,ii,ierr
      logical :: element_found



      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_species_locdipole_groups','cbuf', ierr)

      ! cks: Root proc only: read information from
      !      <species_ngwf_plot> block of input file
      if (pub_on_root) then
         if (esdf_block('species_locdipole_groups',dummy_nsp)) then
            call utils_assert(dummy_nsp <= nsp, &
               'Error in internal_species_locdipole_groups(): &
               &too many species in <species_locdipole_groups> &
               &block of input file.')
            ! aam: check for multiply defined species
            call internal_check_species(dummy_nsp,'<species_locdipole_groups>')
            do row=1,dummy_nsp ! aam: dummy_nsp is not necessarily equal to nsp
               read(block_data(row),*) dummy_id
               cbuf(row)(1:4) = dummy_id
            end do
         end if
      end if

      ! pdh: broadcast
      call comms_bcast(pub_root_proc_id,dummy_nsp)
      do row=1,dummy_nsp
         call comms_bcast(pub_root_proc_id,cbuf(row),4)
      end do

         do ii=1,nat
            elements(ii)%loc_dipole = .false.
         end do

         do row=1,dummy_nsp
            element_found = .false.
            dummy_id = cbuf(row)(1:4)
            do ii = 1, nat
               if ( elements(ii)%species_id == dummy_id) then
                  elements(ii)%loc_dipole = .true.
                  element_found = .true.
               end if
            end do
            if (.not.element_found) then
            call utils_abort('Error in internal_species_locdipole_groups: &
                 &mismatching species in <species_locdipole_groups> and &
                 &<positions_abs> blocks of input file.')
         endif
         end do

         deallocate(cbuf,stat=ierr)
         call utils_dealloc_check('internal_species_locdipole_groups','cbuf', &
              ierr)

         if (pub_debug_on_root) write(stdout,'(a)') &
              'DEBUG: Read <species_locdipole_groups> block'

    end subroutine internal_species_locdipole_groups

    !-----------------------------------------
    !--------- %block species_swri_*
    !-----------------------------------------

    subroutine internal_species_swri(swri_entry)
      ! Cloned and extended on 2015/04/01 by Jacek Dziedzic from
      ! internal_species_locdipole_groups by Valerio Vitale

      use comms, only: pub_root_proc_id, comms_bcast
      use constants, only: max_rundat_swri_block_entries
      use esdf, only: block_data, esdf_block
      use sw_resolution_of_identity, only: swri_rundat_block_data
      use utils, only: utils_abort, utils_assert

      implicit none

      integer, intent(in) :: swri_entry

      ! vv: dummy variables
      character(len=4) :: dummy_id
      integer :: dummy_nsp
      ! vv: dummy buffers
      character(len=68), allocatable :: cbuf(:)
      ! vv: integers and locals
      integer :: row,ii,ierr
      logical :: element_found
      character(len=128) :: blockname

      ! -----------------------------------------------------------------------

      call utils_assert(swri_entry <= max_rundat_swri_block_entries, &
           'Bounds check in internal_species_swri().')

      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_species_swri','cbuf', ierr)

      ! cks: Root proc only: read information from
      !      <species_swri_*> block of input file
      blockname = 'species_swri-'//&
           trim(swri_rundat_block_data(swri_entry)%entry_name)
      if (pub_on_root) then
         if (esdf_block(trim(blockname),dummy_nsp,.true.)) then
            call utils_assert(dummy_nsp <= nsp, &
                 'Error in internal_species_swri(): &
                 &too many species in '//trim(blockname))
            ! aam: check for multiply defined species
            call internal_check_species(dummy_nsp,trim(blockname))
            do row=1,dummy_nsp ! aam: dummy_nsp is not necessarily equal to nsp
               read(block_data(row),*) dummy_id
               cbuf(row)(1:4) = dummy_id
            end do
         end if
      end if

      ! pdh: broadcast
      call comms_bcast(pub_root_proc_id,dummy_nsp)
      do row=1,dummy_nsp
         call comms_bcast(pub_root_proc_id,cbuf(row),4)
      end do

      do ii=1,nat
         elements(ii)%in_swri(swri_entry) = .false.
      end do

      do row=1,dummy_nsp
         element_found = .false.
         dummy_id = cbuf(row)(1:4)
         do ii = 1, nat
            if ( elements(ii)%species_id == dummy_id) then
               elements(ii)%in_swri(swri_entry) = .true.
               element_found = .true.
            end if
         end do
         if (.not.element_found) then
            call utils_abort('Error in internal_species_swri: &
                 &mismatching species in '//trim(blockname)//' and &
                 &<positions_abs> (or equivalent) blocks of input file.')
         endif
      end do

      call utils_assert(any(elements(:)%in_swri(swri_entry)), &
           'SWRI ['//trim(swri_rundat_block_data(swri_entry)%entry_name)//&
           '] has no atoms assigned to it. You either forgot to specify a &
           &"%block '//trim(blockname)//'", or this block is empty, or it &
           &is effectively empty (specifies species none of which have &
           &corresponding atoms). If you simply want to apply the SWRI to &
           &*all* atoms in your system, list all your species names in this &
           &block.')

      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_species_swri','cbuf',ierr)

      if (pub_debug_on_root) write(stdout,'(a,i0)') &
           'DEBUG: Read <species_swri_*> block #', swri_entry

    end subroutine internal_species_swri

    !------------------------------------
    !--------- %block species_ldos_groups
    !------------------------------------

    subroutine internal_species_ldos_groups()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_ldos_ngroups, pub_ldos_group_nsp, pub_ldos_groups
      use utils, only: utils_abort

      implicit none

      !-------------
      ! ars: dummies
      character(len=80) :: dummy_ldos_group
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      !-------------
      ! ars: integers
      integer :: row, ii, jj
      ! ars: logical
      logical :: element_found



      ! ndmh: Root proc only: read information from
      !      <species_ldos_groups> block of input file
      ! ndmh: check if the block exists
      if (pub_on_root) then
         if (esdf_block('species_ldos_groups',pub_ldos_ngroups)) then
            if (pub_ldos_ngroups > nat) then
               call utils_abort('Error in internal_species_ldos_groups(): &
                    &more groups than atoms defined in <species_ldos_groups>.')
            end if
         else
            pub_ldos_ngroups = 0
            dummy_nsp = 0
         end if
      end if

      ! ndmh: allocate array for numbers of species in each group
      call comms_bcast(pub_root_proc_id,pub_ldos_ngroups)
      if (allocated(pub_ldos_group_nsp)) then
         deallocate(pub_ldos_group_nsp,stat=ierr)
         call utils_dealloc_check('internal_species_ldos_groups', &
              'pub_ldos_group_nsp',ierr)
      end if
      allocate(pub_ldos_group_nsp(pub_ldos_ngroups),stat=ierr)
      call utils_alloc_check('internal_species_ldos_groups', &
           'pub_ldos_group_nsp',ierr)

      ! ndmh: count number of species in each group
      if (pub_on_root.and.(pub_ldos_ngroups > 0)) then
         dummy_nsp = 0
         ! loop over lines in block (LDOS groups)
         do row=1,pub_ldos_ngroups
            dummy_ldos_group = trim(adjustl(block_data(row)))
            pub_ldos_group_nsp(row) = 1
            ! ndmh: count number of spaces (distinct species) in this group
            jj = 0
            do ii=1,len_trim(dummy_ldos_group)
               if ((dummy_ldos_group(ii:ii)==' ').and. &
                    (dummy_ldos_group(ii+1:ii+1)/=' ')) then
                  if (ii > jj + 1) then
                     pub_ldos_group_nsp(row) = pub_ldos_group_nsp(row) + 1
                     jj = ii
                  else
                     jj = ii
                  end if
               end if
            end do
         end do
      end if

      ! ndmh: allocate array for storing LDOS group species ID's
      call comms_bcast(pub_root_proc_id,pub_ldos_group_nsp)
      dummy_nsp = maxval(pub_ldos_group_nsp)
      call comms_bcast(pub_root_proc_id,dummy_nsp)
      if (allocated(pub_ldos_groups)) then
         deallocate(pub_ldos_groups,stat=ierr)
         call utils_dealloc_check('internal_species_ldos_groups','pub_ldos_groups',&
              ierr)
      end if
      allocate(pub_ldos_groups(dummy_nsp,pub_ldos_ngroups),stat=ierr)
      call utils_alloc_check('internal_species_ldos_groups','pub_ldos_groups',ierr)
      pub_ldos_groups(:,:) = ''

      ! ndmh: now fill the pub_ldos_groups array with these strings
      if (pub_on_root.and.(pub_ldos_ngroups > 0)) then
         do row=1,pub_ldos_ngroups
            dummy_ldos_group = trim(adjustl(block_data(row)))
            ii = 1
            jj = len_trim(dummy_ldos_group)
            do dummy_nsp=1,pub_ldos_group_nsp(row)
               if (index(dummy_ldos_group(ii:jj),' ') > 0) then
                  jj = ii + index(dummy_ldos_group(ii:jj),' ') - 1
                  do
                     if (dummy_ldos_group(jj+1:jj+1)/=' ') exit
                     jj = jj + 1
                  end do
               end if
               pub_ldos_groups(dummy_nsp,row) = &
                    adjustl(trim(dummy_ldos_group(ii:jj)))
               ii = jj + 1
               jj = len_trim(dummy_ldos_group)
            end do
         end do
      end if

      ! ndmh: broadcast LDOS groups from root proc
      do row=1,pub_ldos_ngroups
         do dummy_nsp=1,pub_ldos_group_nsp(row)
            call comms_bcast(pub_root_proc_id,pub_ldos_groups(dummy_nsp,row))
         end do
      end do

      ! ndmh: check LDOS group species match elements array species id's
      do row=1,pub_ldos_ngroups
         do ii=1,pub_ldos_group_nsp(row)
            element_found = .false.
            dummy_id = pub_ldos_groups(ii,row)
            do jj=1,nat
               if (elements(jj)%species_id == dummy_id) then
                  element_found = .true.
               end if
            end do
            if (.not. element_found) then
               call utils_abort('Error in internal_species_ldos_groups: &
                    &mismatching species in <species_ldos_groups> and &
                    &<positions_abs> blocks of input file.')
            endif
         end do
      end do

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_ldos_groups> block'

    end subroutine internal_species_ldos_groups

    !------------------------------------
    !--------- %block species_bsunfld_groups
    !------------------------------------

    subroutine internal_species_bsunfld_groups()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_bsunfld_ngroups, pub_bsunfld_group_nsp, &
           pub_bsunfld_groups
      use utils, only: utils_abort

      implicit none

      !-------------
      ! gcc: dummies
      character(len=80) :: dummy_bsunfld_group
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      !-------------
      ! gcc: integers
      integer :: row, ii, jj
      ! gcc: logical
      logical :: element_found



      ! ndmh: Root proc only: read information from
      !      <species_bsunfld_groups> block of input file
      ! ndmh: check if the block exists
      if (pub_on_root) then
         if (esdf_block('species_bsunfld_groups',pub_bsunfld_ngroups)) then
            ! gcc32: allow at most a single group FOR THE MOMENT
            if (pub_bsunfld_ngroups > 1) then
               call utils_abort('Error in internal_species_bsunfld_groups(): &
                    &more groups than atoms defined in &
                    &<species_bsunfld_groups>.')
            end if
         else
            pub_bsunfld_ngroups = 0
            dummy_nsp = 0
         end if
      end if

      ! ndmh: allocate array for numbers of species in each group
      call comms_bcast(pub_root_proc_id,pub_bsunfld_ngroups)
      if (allocated(pub_bsunfld_group_nsp)) then
         deallocate(pub_bsunfld_group_nsp,stat=ierr)
         call utils_dealloc_check('internal_species_bsunfld_groups', &
              'pub_bsunfld_group_nsp',ierr)
      end if
      allocate(pub_bsunfld_group_nsp(pub_bsunfld_ngroups),stat=ierr)
      call utils_alloc_check('internal_species_bsunfld_groups', &
           'pub_bsunfld_group_nsp',ierr)

      ! ndmh: count number of species in each group
      if (pub_on_root.and.(pub_bsunfld_ngroups > 0)) then
         dummy_nsp = 0
         ! loop over lines in block (BSUNFLD groups)
         do row=1,pub_bsunfld_ngroups
            dummy_bsunfld_group = trim(adjustl(block_data(row)))
            pub_bsunfld_group_nsp(row) = 1
            ! ndmh: count number of spaces (distinct species) in this group
            jj = 0
            do ii=1,len_trim(dummy_bsunfld_group)
               if ((dummy_bsunfld_group(ii:ii)==' ').and. &
                    (dummy_bsunfld_group(ii+1:ii+1)/=' ')) then
                  if (ii > jj + 1) then
                     pub_bsunfld_group_nsp(row) = pub_bsunfld_group_nsp(row) + 1
                     jj = ii
                  else
                     jj = ii
                  end if
               end if
            end do
         end do
      end if

      ! ndmh: allocate array for storing BSUNFLD group species ID's
      call comms_bcast(pub_root_proc_id,pub_bsunfld_group_nsp)
      dummy_nsp = maxval(pub_bsunfld_group_nsp)
      call comms_bcast(pub_root_proc_id,dummy_nsp)
      if (allocated(pub_bsunfld_groups)) then
         deallocate(pub_bsunfld_groups,stat=ierr)
         call utils_dealloc_check('internal_species_bsunfld_groups','pub_bsunfld_groups',&
              ierr)
      end if
      allocate(pub_bsunfld_groups(dummy_nsp,pub_bsunfld_ngroups),stat=ierr)
      call utils_alloc_check('internal_species_bsunfld_groups','pub_bsunfld_groups',ierr)
      pub_bsunfld_groups(:,:) = ''

      ! ndmh: now fill the pub_bsunfld_groups array with these strings
      if (pub_on_root.and.(pub_bsunfld_ngroups > 0)) then
         do row=1,pub_bsunfld_ngroups
            dummy_bsunfld_group = trim(adjustl(block_data(row)))
            ii = 1
            jj = len_trim(dummy_bsunfld_group)
            do dummy_nsp=1,pub_bsunfld_group_nsp(row)
               if (index(dummy_bsunfld_group(ii:jj),' ') > 0) then
                  jj = ii + index(dummy_bsunfld_group(ii:jj),' ') - 1
                  do
                     if (dummy_bsunfld_group(jj+1:jj+1)/=' ') exit
                     jj = jj + 1
                  end do
               end if
               pub_bsunfld_groups(dummy_nsp,row) = &
                    adjustl(trim(dummy_bsunfld_group(ii:jj)))
               ii = jj + 1
               jj = len_trim(dummy_bsunfld_group)
            end do
         end do
      end if

      ! ndmh: broadcast BSUNFLD groups from root proc
      do row=1,pub_bsunfld_ngroups
         do dummy_nsp=1,pub_bsunfld_group_nsp(row)
            call comms_bcast(pub_root_proc_id,pub_bsunfld_groups(dummy_nsp,row))
         end do
      end do

      ! ndmh: check BSUNFLD group species match elements array species id's
      do row=1,pub_bsunfld_ngroups
         do ii=1,pub_bsunfld_group_nsp(row)
            element_found = .false.
            dummy_id = pub_bsunfld_groups(ii,row)
            do jj=1,nat
               if (elements(jj)%species_id == dummy_id) then
                  element_found = .true.
               end if
            end do
            if (.not. element_found) then
               call utils_abort('Error in internal_species_bsunfld_groups: &
                    &mismatching species in <species_bsunfld_groups> and &
                    &<positions_abs> blocks of input file.')
            endif
         end do
      end do

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_bsunfld_groups> block'

    end subroutine internal_species_bsunfld_groups

    !------------------------------------
    !--------- %block species_bsunfld_projatoms
    !------------------------------------

    subroutine internal_species_bsunfld_projatoms()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_bsunfld_nprojatoms, &
           pub_bsunfld_projatoms_nsp, pub_bsunfld_projatoms
      use utils, only: utils_abort

      implicit none

      !-------------
      ! gcc: dummies
      character(len=80) :: dummy_bsunfld_projatoms
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      !-------------
      ! gcc: integers
      integer :: row, ii, jj
      ! gcc: logical
      logical :: element_found



      ! ndmh: Root proc only: read information from
      !      <species_bsunfld_projatoms> block of input file
      ! ndmh: check if the block exists
      if (pub_on_root) then
         if (esdf_block('species_bsunfld_projatoms', &
             pub_bsunfld_nprojatoms)) then
         else
            pub_bsunfld_nprojatoms = 0
            dummy_nsp = 0
         end if
      end if

      ! ndmh: allocate array for numbers of species in each group
      call comms_bcast(pub_root_proc_id,pub_bsunfld_nprojatoms)
      if (allocated(pub_bsunfld_projatoms_nsp)) then
         deallocate(pub_bsunfld_projatoms_nsp,stat=ierr)
         call utils_dealloc_check('internal_species_bsunfld_projatoms', &
              'pub_bsunfld_projatoms_nsp',ierr)
      end if
      allocate(pub_bsunfld_projatoms_nsp(pub_bsunfld_nprojatoms), &
          stat=ierr)
      call utils_alloc_check('internal_species_bsunfld_projatoms', &
           'pub_bsunfld_projatoms_nsp',ierr)

      ! ndmh: count number of species in each group
      if (pub_on_root.and.(pub_bsunfld_nprojatoms > 0)) then
         dummy_nsp = 0
         ! loop over lines in block (BSUNFLD groups)
         do row=1,pub_bsunfld_nprojatoms
            dummy_bsunfld_projatoms = trim(adjustl(block_data(row)))
            pub_bsunfld_projatoms_nsp(row) = 1
            ! ndmh: count number of spaces (distinct species) in this group
            jj = 0
            do ii=1,len_trim(dummy_bsunfld_projatoms)
               if ((dummy_bsunfld_projatoms(ii:ii)==' ').and. &
                    (dummy_bsunfld_projatoms(ii+1:ii+1)/=' ')) then
                  if (ii > jj + 1) then
                     pub_bsunfld_projatoms_nsp(row) = &
                          pub_bsunfld_projatoms_nsp(row) + 1
                     jj = ii
                  else
                     jj = ii
                  end if
               end if
            end do
         end do
      end if

      ! ndmh: allocate array for storing BSUNFLD group species ID's
      call comms_bcast(pub_root_proc_id,pub_bsunfld_projatoms_nsp)
      dummy_nsp = maxval(pub_bsunfld_projatoms_nsp)
      call comms_bcast(pub_root_proc_id,dummy_nsp)
      if (allocated(pub_bsunfld_projatoms)) then
         deallocate(pub_bsunfld_projatoms,stat=ierr)
         call utils_dealloc_check('internal_species_bsunfld_projatoms', &
              'pub_bsunfld_projatoms', ierr)
      end if
      allocate(pub_bsunfld_projatoms(dummy_nsp, &
           pub_bsunfld_nprojatoms),stat=ierr)
      call utils_alloc_check('internal_species_bsunfld_projatoms', &
           'pub_bsunfld_projatoms',ierr)
      pub_bsunfld_projatoms(:,:) = ''

      ! ndmh: now fill the pub_bsunfld_projatoms array with these strings
      if (pub_on_root.and.(pub_bsunfld_nprojatoms > 0)) then
         do row=1,pub_bsunfld_nprojatoms
            dummy_bsunfld_projatoms = trim(adjustl(block_data(row)))
            ii = 1
            jj = len_trim(dummy_bsunfld_projatoms)
            do dummy_nsp=1,pub_bsunfld_projatoms_nsp(row)
               if (index(dummy_bsunfld_projatoms(ii:jj),' ') > 0) then
                  jj = ii + index(dummy_bsunfld_projatoms(ii:jj),' ') - 1
                  do
                     if (dummy_bsunfld_projatoms(jj+1:jj+1)/=' ') exit
                     jj = jj + 1
                  end do
               end if
               pub_bsunfld_projatoms(dummy_nsp,row) = &
                    adjustl(trim(dummy_bsunfld_projatoms(ii:jj)))
               ii = jj + 1
               jj = len_trim(dummy_bsunfld_projatoms)
            end do
         end do
      end if

      ! ndmh: broadcast BSUNFLD groups from root proc
      do row=1,pub_bsunfld_nprojatoms
         do dummy_nsp=1,pub_bsunfld_projatoms_nsp(row)
            call comms_bcast(pub_root_proc_id,pub_bsunfld_projatoms(dummy_nsp, &
                 row))
         end do
      end do

      ! ndmh: check BSUNFLD group species match elements array species id's
      do row=1,pub_bsunfld_nprojatoms
         do ii=1,pub_bsunfld_projatoms_nsp(row)
            element_found = .false.
            dummy_id = pub_bsunfld_projatoms(ii,row)
            do jj=1,nat
               if (elements(jj)%species_id == dummy_id) then
                  element_found = .true.
               end if
            end do
            if (.not. element_found) then
               call utils_abort('Error in internal_species_bsunfld_projatoms: &
                    &mismatching species in <species_bsunfld_projatoms> and &
                    &<positions_abs> blocks of input file.')
            endif
         end do
      end do

      ! jcap: deallocate to avoid memory leakage
      deallocate(pub_bsunfld_projatoms_nsp,stat=ierr)
      call utils_dealloc_check('internal_species_bsunfld_projatoms', &
           'pub_bsunfld_projatoms_nsp',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_bsunfld_projatoms> block'

    end subroutine internal_species_bsunfld_projatoms


    !------------------------------------
    !--------- %block species_pdos_groups
    !------------------------------------

    subroutine internal_species_pdos_groups()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_pdos_ngroups, pub_pdos_group_nsp, pub_pdos_groups
      use utils, only: utils_abort

      implicit none

      !-------------
      ! ars: dummies
      character(len=80) :: dummy_pdos_group
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      !-------------
      ! ars: integers
      integer :: row, ii, jj
      ! ars: logical
      logical :: element_found



      ! ndmh: Root proc only: read information from
      !      <species_pdos_groups> block of input file
      ! ndmh: check if the block exists
      if (pub_on_root) then
         if (esdf_block('species_pdos_groups',pub_pdos_ngroups)) then
            if (pub_pdos_ngroups > nat) then
               call utils_abort('Error in internal_species_pdos_groups(): &
                    &more groups than atoms defined in <species_pdos_groups>.')
            end if
         else
            pub_pdos_ngroups = 0
            dummy_nsp = 0
         end if
      end if

      ! ndmh: allocate array for numbers of species in each group
      call comms_bcast(pub_root_proc_id,pub_pdos_ngroups)
      if (allocated(pub_pdos_group_nsp)) then
         deallocate(pub_pdos_group_nsp,stat=ierr)
         call utils_dealloc_check('internal_species_pdos_groups', &
              'pub_pdos_group_nsp',ierr)
      end if
      allocate(pub_pdos_group_nsp(pub_pdos_ngroups),stat=ierr)
      call utils_alloc_check('internal_species_pdos_groups', &
           'pub_pdos_group_nsp',ierr)

      ! ndmh: count number of species in each group
      if (pub_on_root.and.(pub_pdos_ngroups > 0)) then
         dummy_nsp = 0
         ! loop over lines in block (PDOS groups)
         do row=1,pub_pdos_ngroups
            dummy_pdos_group = trim(adjustl(block_data(row)))
            pub_pdos_group_nsp(row) = 1
            ! ndmh: count number of spaces (distinct species) in this group
            jj = 0
            do ii=1,len_trim(dummy_pdos_group)
               if ((dummy_pdos_group(ii:ii)==' ').and. &
                    (dummy_pdos_group(ii+1:ii+1)/=' ')) then
                  if (ii > jj + 1) then
                     pub_pdos_group_nsp(row) = pub_pdos_group_nsp(row) + 1
                     jj = ii
                  else
                     jj = ii
                  end if
               end if
            end do
         end do
      end if

      ! ndmh: allocate array for storing PDOS group species ID's
      call comms_bcast(pub_root_proc_id,pub_pdos_group_nsp)
      dummy_nsp = maxval(pub_pdos_group_nsp)
      call comms_bcast(pub_root_proc_id,dummy_nsp)
      if (allocated(pub_pdos_groups)) then
         deallocate(pub_pdos_groups,stat=ierr)
         call utils_dealloc_check('internal_species_pdos_groups','pub_pdos_groups',&
              ierr)
      end if
      allocate(pub_pdos_groups(dummy_nsp,pub_pdos_ngroups),stat=ierr)
      call utils_alloc_check('internal_species_pdos_groups','pub_pdos_groups',ierr)
      pub_pdos_groups(:,:) = ''

      ! ndmh: now fill the pub_pdos_groups array with these strings
      if (pub_on_root.and.(pub_pdos_ngroups > 0)) then
         do row=1,pub_pdos_ngroups
            dummy_pdos_group = trim(adjustl(block_data(row)))
            ii = 1
            jj = len_trim(dummy_pdos_group)
            do dummy_nsp=1,pub_pdos_group_nsp(row)
               if (index(dummy_pdos_group(ii:jj),' ') > 0) then
                  jj = ii + index(dummy_pdos_group(ii:jj),' ') - 1
                  do
                     if (dummy_pdos_group(jj+1:jj+1)/=' ') exit
                     jj = jj + 1
                  end do
               end if
               pub_pdos_groups(dummy_nsp,row) = &
                    adjustl(trim(dummy_pdos_group(ii:jj)))
               ii = jj + 1
               jj = len_trim(dummy_pdos_group)
            end do
         end do
      end if

      ! ndmh: broadcast PDOS groups from root proc
      do row=1,pub_pdos_ngroups
         do dummy_nsp=1,pub_pdos_group_nsp(row)
            call comms_bcast(pub_root_proc_id,pub_pdos_groups(dummy_nsp,row))
         end do
      end do

      ! ndmh: check PDOS group species match elements array species id's
      do row=1,pub_pdos_ngroups
         do ii=1,pub_pdos_group_nsp(row)
            element_found = .false.
            dummy_id = pub_pdos_groups(ii,row)
            do jj=1,nat
               if (elements(jj)%species_id == dummy_id) then
                  element_found = .true.
               end if
            end do
            if (.not. element_found) then
               call utils_abort('Error in internal_species_pdos_groups: &
                    &mismatching species in <species_pdos_groups> and &
                    &<positions_abs> blocks of input file.')
            endif
         end do
      end do

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_pdos_groups> block'

    end subroutine internal_species_pdos_groups

    !------------------------------------
    !--------- %block species_scissor
    !------------------------------------

    subroutine internal_species_scissor()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_convfac
      use rundat, only: pub_scissor_ngroups, pub_scissor_group_nsp, &
           pub_scissor_groups, pub_scissor_group_shifts
      use utils, only: utils_abort, utils_assert

      implicit none

      character(len=4)  :: dummy_id
      character(len=80) :: scissor_group
      character(len=80) :: dummy_shift
      integer :: ngroups, nsp, nshifts
      integer :: iat
      integer :: read_success
      integer :: row, ii, jj
      logical :: element_found
      logical :: shift_missing
      real(kind=dp) :: shift, confac

      ! Root proc only: read information from
      !    <species_scissor> block of input file
      ! check if the block exists
      if (pub_on_root) then
         if (esdf_block('species_scissor',ngroups)) then
            call utils_assert(ngroups <= nat, 'Error in internal_species_&
                 &scissor(): too many groups in <species_scissor> &
                 &block of input file.')
         else
            ngroups = 0
         end if
      end if

      call comms_bcast(pub_root_proc_id,ngroups)

      ! allocate array for number of species in each group
      if (allocated(pub_scissor_group_nsp)) then
         deallocate(pub_scissor_group_nsp,stat=ierr)
         call utils_dealloc_check('internal_species_scissor', &
              'pub_scissor_group_nsp',ierr)
      end if
      allocate(pub_scissor_group_nsp(ngroups),stat=ierr)
      call utils_alloc_check('internal_species_scissor', &
           'pub_scissor_group_nsp',ierr)

      ! allocate array for valence and conduction shifts in each group
      if (allocated(pub_scissor_group_shifts)) then
         deallocate(pub_scissor_group_shifts,stat=ierr)
         call utils_dealloc_check('internal_species_scissor', &
              'pub_scissor_group_shifts',ierr)
      end if
      allocate(pub_scissor_group_shifts(2,ngroups),stat=ierr)
      call utils_alloc_check('internal_species_scissor', &
           'pub_scissor_group_shifts',ierr)

      shift_missing=.false.
      nsp = 0
      ! count number of species and shifts in each group
      if (pub_on_root.and.(ngroups > 0)) then
         ! loop over lines in block (scissor groups)
         do row=1,ngroups
            scissor_group = trim(adjustl(block_data(row)))
            pub_scissor_group_nsp(row) = 1
            nshifts = 0
            jj = len_trim(scissor_group)

            ! loop through each character until it reaches spaces
            ! this starts from the end because we need to find two shift
            ! parameters first
            do ii=len_trim(scissor_group),2,-1
               if ((scissor_group(ii:ii)==' ') .and. &
                    (scissor_group(ii-1:ii-1)/=' ') .and. ii<jj-1) then
                  ! if it's numeric then it must be a shift
                  ! otherwise it's a specie
                  dummy_shift = trim(adjustl(scissor_group(ii:jj)))
                  read(dummy_shift,*,iostat=read_success) shift
                  if (read_success==0) then
                     nshifts = nshifts + 1
                     confac = esdf_convfac('eV','hartree')
                     shift = shift*confac
                     pub_scissor_group_shifts(nshifts,row) = shift
                  else
                     pub_scissor_group_nsp(row) = pub_scissor_group_nsp(row) + 1
                  end if
                  jj = ii
               end if
            end do

            ! there should be exactly two shifts
            if(nshifts/=2) then
               shift_missing=.true.
               exit
            end if
         end do
      end if

      call comms_bcast(pub_root_proc_id,pub_scissor_group_shifts)
      call comms_bcast(pub_root_proc_id, shift_missing)
      call utils_assert(.not.shift_missing, &
           'Error in rundat_blocks->internal_species_scissor: &
           &you must include exactly 2 shift parameters for each group')

      ! allocate array for storing scissor-group species ID's
      call comms_bcast(pub_root_proc_id,pub_scissor_group_nsp)
      nsp = maxval(pub_scissor_group_nsp)
      call comms_bcast(pub_root_proc_id,nsp)
      if (allocated(pub_scissor_groups)) then
         deallocate(pub_scissor_groups,stat=ierr)
         call utils_dealloc_check('internal_species_scissor', &
              'pub_scissor_groups', ierr)
      end if
      allocate(pub_scissor_groups(nsp,ngroups),stat=ierr)
      call utils_alloc_check('internal_species_scissor','pub_scissor_groups',ierr)
      pub_scissor_groups(:,:) = ''

      ! now fill the pub_scissor array with these strings
      if (pub_on_root.and.(ngroups > 0)) then
         do row=1,ngroups
            scissor_group = trim(adjustl(block_data(row)))
            ii = 1
            jj = len_trim(scissor_group)
            do nsp=1,pub_scissor_group_nsp(row)
               if (index(scissor_group(ii:jj),' ') > 0) then
                  jj = ii + index(scissor_group(ii:jj),' ') - 1
                  do
                     if (scissor_group(jj+1:jj+1)/=' ') exit
                     jj = jj + 1
                  end do
               end if
               pub_scissor_groups(nsp,row) = &
                    adjustl(trim(scissor_group(ii:jj)))
               ii = jj + 1
               jj = len_trim(scissor_group)
            end do
         end do
      end if

      ! broadcast scissor groups from root proc
      do row=1,ngroups
         do nsp=1,pub_scissor_group_nsp(row)
            call comms_bcast(pub_root_proc_id,pub_scissor_groups(nsp,row))
         end do
      end do

      ! check scissor group species match elements array species id's
      do row=1,pub_scissor_ngroups
         do ii=1,pub_scissor_group_nsp(row)
            element_found = .false.
            dummy_id = pub_scissor_groups(ii,row)
            do jj=1,nat
               if (elements(jj)%species_id == dummy_id) then
                  element_found = .true.
               end if
            end do
            if (.not. element_found) then
               call utils_abort('Error in internal_species_scissor: &
                    &mismatching species in <species_scissor> and &
                    &<positions_abs> blocks of input file.')
            endif
         end do
      end do

      pub_scissor_ngroups = ngroups
      call comms_bcast(pub_root_proc_id, pub_scissor_ngroups)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_scissor> block'

    end subroutine internal_species_scissor

    !--------------------------------
    !--------- %block tddft_kernel_ct
    !-------------------------------
    ! while the normal tddft_kernel block defines the
    ! general regions tddft_kernel_ct defines charge transfer
    ! regions between the general block layout
    ! tjz21: Mainly copied from LDOS block
    subroutine internal_species_tddft_ct()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_tddft_ct_ngroups, pub_tddft_ct_group_nsp,&
           pub_tddft_ct_groups,pub_lr_tddft_sparse_region, &
           pub_lr_tddft_subsystem_coupling
      use utils, only: utils_abort

      implicit none

      !-------------
      ! ars: dummies
      character(len=80) :: dummy_tddft_ct_group
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      !-------------
      ! ars: integers
      integer :: row, ii, jj
      ! ars: logical
      logical :: element_found

      ! ndmh: Root proc only: read information from
      !      <species_ldos_groups> block of input file
      ! ndmh: check if the block exists
      if (pub_on_root) then
         if (esdf_block('species_tddft_ct',pub_tddft_ct_ngroups)) then
            if (pub_tddft_ct_ngroups > nat) then
               call utils_abort('Error in internal_species_tddft_ct(): &
                    &more groups than atoms defined in <species_tddft_ct>.')
            end if
         else
            pub_tddft_ct_ngroups = 0
            dummy_nsp = 0
         end if
      end if
      ! ndmh: allocate array for numbers of species in each group
      call comms_bcast(pub_root_proc_id,pub_tddft_ct_ngroups)
      if (allocated(pub_tddft_ct_group_nsp)) then
         deallocate(pub_tddft_ct_group_nsp,stat=ierr)
         call utils_dealloc_check('internal_species_tddft_ct', &
              'pub_tddft_ct_group_nsp',ierr)
      end if
      allocate(pub_tddft_ct_group_nsp(pub_tddft_ct_ngroups),stat=ierr)
      call utils_alloc_check('internal_species_tddft_ct', &
           'pub_tddft_ct_group_nsp',ierr)

      ! ndmh: count number of species in each group
      if (pub_on_root.and.(pub_tddft_ct_ngroups > 0)) then
         dummy_nsp = 0
         ! loop over lines in block (LDOS groups)
         do row=1,pub_tddft_ct_ngroups
            dummy_tddft_ct_group = trim(adjustl(block_data(row)))
            pub_tddft_ct_group_nsp(row) = 1
            ! ndmh: count number of spaces (distinct species) in this group
            jj = 0
            do ii=1,len_trim(dummy_tddft_ct_group)
               if ((dummy_tddft_ct_group(ii:ii)==' ').and. &
                    (dummy_tddft_ct_group(ii+1:ii+1)/=' ')) then
                  if (ii > jj + 1) then
                     pub_tddft_ct_group_nsp(row) = pub_tddft_ct_group_nsp(row) + 1
                     jj = ii
                  else
                     jj = ii
                  end if
               end if
            end do
         end do
      end if

      ! ndmh: allocate array for storing LDOS group species ID's
      call comms_bcast(pub_root_proc_id,pub_tddft_ct_group_nsp)
      dummy_nsp = maxval(pub_tddft_ct_group_nsp)
      call comms_bcast(pub_root_proc_id,dummy_nsp)
      if (allocated(pub_tddft_ct_groups)) then
         deallocate(pub_tddft_ct_groups,stat=ierr)
         call utils_dealloc_check('internal_species_tddft_ct','pub_tddft_ct_groups',&
              ierr)
      end if
      allocate(pub_tddft_ct_groups(dummy_nsp,pub_tddft_ct_ngroups),stat=ierr)
      call utils_alloc_check('internal_species_tddft_ct','pub_tddft_ct_groups',ierr)
      pub_tddft_ct_groups(:,:) = ''

      ! ndmh: now fill the pub_ldos_groups array with these strings
      if (pub_on_root.and.(pub_tddft_ct_ngroups > 0)) then
         do row=1,pub_tddft_ct_ngroups
            dummy_tddft_ct_group = trim(adjustl(block_data(row)))
            ii = 1
            jj = len_trim(dummy_tddft_ct_group)
            do dummy_nsp=1,pub_tddft_ct_group_nsp(row)
               if (index(dummy_tddft_ct_group(ii:jj),' ') > 0) then
                  jj = ii + index(dummy_tddft_ct_group(ii:jj),' ') - 1
                  do
                     if (dummy_tddft_ct_group(jj+1:jj+1)/=' ') exit
                     jj = jj + 1
                  end do
               end if
               pub_tddft_ct_groups(dummy_nsp,row) = &
                    adjustl(trim(dummy_tddft_ct_group(ii:jj)))
               ii = jj + 1
               jj = len_trim(dummy_tddft_ct_group)
            end do
         end do
      end if

      ! ndmh: broadcast LDOS groups from root proc
      do row=1,pub_tddft_ct_ngroups
         do dummy_nsp=1,pub_tddft_ct_group_nsp(row)
            call comms_bcast(pub_root_proc_id,pub_tddft_ct_groups(dummy_nsp,row))
         end do
      end do

      ! ndmh: check LDOS group species match elements array species id's
      do row=1,pub_tddft_ct_ngroups
         do ii=1,pub_tddft_ct_group_nsp(row)
            element_found = .false.
            dummy_id = pub_tddft_ct_groups(ii,row)
            do jj=1,nat
               if (elements(jj)%species_id == dummy_id) then
                  element_found = .true.
               end if
            end do
            if (.not. element_found) then
               call utils_abort('Error in internal_species_tddft_ct: &
                    &mismatching species in <species_tddft_ct> and &
                    &<positions_abs> blocks of input file.')
            endif
         end do
      end do

      ! sanity check. If this block is present, the standard TDDFT kernel block
      ! should also be present and not empty
      if(pub_tddft_ct_ngroups>0 .and. .not. pub_lr_tddft_subsystem_coupling) then
         if(.not. pub_lr_tddft_sparse_region) then
            call utils_abort('Error in internal_species_tddft_ct: &
                &species_tddft_ct block must be used in conjunction with &
                &species_tddft_kernel block')
         endif
      endif

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_tddft_ct> block'

    end subroutine internal_species_tddft_ct

    !-----------------------------
    !--------- %block tddft_kernel
    !-----------------------------
    ! tjz21: Mainly copied from LDOS block
    subroutine internal_species_tddft_kernel()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_tddft_kernel_ngroups, pub_tddft_kernel_group_nsp,&
           pub_tddft_kernel_groups,pub_lr_tddft_sparse_region, &
           pub_lr_tddft_subsystem_coupling
      use utils, only: utils_abort

      implicit none

      !-------------
      ! ars: dummies
      character(len=80) :: dummy_tddft_kernel_group
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      !-------------
      ! ars: integers
      integer :: row, ii, jj
      ! ars: logical
      logical :: element_found



      ! ndmh: Root proc only: read information from
      !      <species_ldos_groups> block of input file
      ! ndmh: check if the block exists
      if (pub_on_root) then
         if (esdf_block('species_tddft_kernel',pub_tddft_kernel_ngroups)) then
            if (pub_tddft_kernel_ngroups > nat) then
               call utils_abort('Error in internal_species_tddft_kernel(): &
                    &more groups than atoms defined in <species_tddft_kernel>.')
            end if
         else
            pub_tddft_kernel_ngroups = 0
            dummy_nsp = 0
         end if
      end if
      ! ndmh: allocate array for numbers of species in each group
      call comms_bcast(pub_root_proc_id,pub_tddft_kernel_ngroups)
      if (allocated(pub_tddft_kernel_group_nsp)) then
         deallocate(pub_tddft_kernel_group_nsp,stat=ierr)
         call utils_dealloc_check('internal_species_tddft_kernel', &
              'pub_tddft_kernel_group_nsp',ierr)
      end if
      allocate(pub_tddft_kernel_group_nsp(pub_tddft_kernel_ngroups),stat=ierr)
      call utils_alloc_check('internal_species_tddft_kernel', &
           'pub_tddft_kernel_group_nsp',ierr)

      ! ndmh: count number of species in each group
      if (pub_on_root.and.(pub_tddft_kernel_ngroups > 0)) then
         dummy_nsp = 0
         ! loop over lines in block (LDOS groups)
         do row=1,pub_tddft_kernel_ngroups
            dummy_tddft_kernel_group = trim(adjustl(block_data(row)))
            pub_tddft_kernel_group_nsp(row) = 1
            ! ndmh: count number of spaces (distinct species) in this group
            jj = 0
            do ii=1,len_trim(dummy_tddft_kernel_group)
               if ((dummy_tddft_kernel_group(ii:ii)==' ').and. &
                    (dummy_tddft_kernel_group(ii+1:ii+1)/=' ')) then
                  if (ii > jj + 1) then
                     pub_tddft_kernel_group_nsp(row) = pub_tddft_kernel_group_nsp(row) + 1
                     jj = ii
                  else
                     jj = ii
                  end if
               end if
            end do
         end do
      end if

      ! ndmh: allocate array for storing LDOS group species ID's
      call comms_bcast(pub_root_proc_id,pub_tddft_kernel_group_nsp)
      dummy_nsp = maxval(pub_tddft_kernel_group_nsp)
      call comms_bcast(pub_root_proc_id,dummy_nsp)
      if (allocated(pub_tddft_kernel_groups)) then
         deallocate(pub_tddft_kernel_groups,stat=ierr)
         call utils_dealloc_check('internal_species_tddft_kernel','pub_tddft_kernel_groups',&
              ierr)
      end if
      allocate(pub_tddft_kernel_groups(dummy_nsp,pub_tddft_kernel_ngroups),stat=ierr)
      call utils_alloc_check('internal_species_tddft_kernel','pub_tddft_kernel_groups',ierr)
      pub_tddft_kernel_groups(:,:) = ''

      ! ndmh: now fill the pub_ldos_groups array with these strings
      if (pub_on_root.and.(pub_tddft_kernel_ngroups > 0)) then
         do row=1,pub_tddft_kernel_ngroups
            dummy_tddft_kernel_group = trim(adjustl(block_data(row)))
            ii = 1
            jj = len_trim(dummy_tddft_kernel_group)
            do dummy_nsp=1,pub_tddft_kernel_group_nsp(row)
               if (index(dummy_tddft_kernel_group(ii:jj),' ') > 0) then
                  jj = ii + index(dummy_tddft_kernel_group(ii:jj),' ') - 1
                  do
                     if (dummy_tddft_kernel_group(jj+1:jj+1)/=' ') exit
                     jj = jj + 1
                  end do
               end if
               pub_tddft_kernel_groups(dummy_nsp,row) = &
                    adjustl(trim(dummy_tddft_kernel_group(ii:jj)))
               ii = jj + 1
               jj = len_trim(dummy_tddft_kernel_group)
            end do
         end do
      end if

      ! ndmh: broadcast LDOS groups from root proc
      do row=1,pub_tddft_kernel_ngroups
         do dummy_nsp=1,pub_tddft_kernel_group_nsp(row)
            call comms_bcast(pub_root_proc_id,pub_tddft_kernel_groups(dummy_nsp,row))
         end do
      end do

      ! ndmh: check LDOS group species match elements array species id's
      do row=1,pub_tddft_kernel_ngroups
         do ii=1,pub_tddft_kernel_group_nsp(row)
            element_found = .false.
            dummy_id = pub_tddft_kernel_groups(ii,row)
            do jj=1,nat
               if (elements(jj)%species_id == dummy_id) then
                  element_found = .true.
               end if
            end do
            if (.not. element_found) then
               call utils_abort('Error in internal_species_tddft_kernel: &
                    &mismatching species in <species_tddft_kernel> and &
                    &<positions_abs> blocks of input file.')
            endif
         end do
      end do

      ! only set sparse_region to .true. if this is not a subsystem couling
      ! calculation. In a coupling calculation the kernels read in are from
      ! different subsystem regions and should therefore be treated without
      ! a region specific truncation
      if(pub_tddft_kernel_ngroups>0 .and. .not. pub_lr_tddft_subsystem_coupling) then
         pub_lr_tddft_sparse_region=.true.
      else
         pub_lr_tddft_sparse_region=.false.
      endif

      call comms_bcast(pub_root_proc_id,pub_lr_tddft_sparse_region)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_tddft_kernel> block'

    end subroutine internal_species_tddft_kernel


    !------------------------
    !--------- %block hubbard
    !------------------------

    subroutine internal_hubbard()

      use comms, only: pub_root_proc_id, comms_bcast
      use constants, only: HARTREE_IN_EVS
      use esdf, only: block_data, esdf_block
      use hubbard_init, only: hubbard_init_species, h_species
      use rundat, only: pub_hubbard, pub_hubbard_atomsolve, pub_hubbard_restart,&
           pub_task, pub_hub_read_projectors, pub_hub_proj_mixing, &
           pub_cdft_multi_proj
      use utils, only: utils_abort

      implicit none

      !-------------
      ! ars: dummies
      character(len=4)  :: dummy_id
      integer           :: dummy_l, dummy_h
      real(kind=dp)     :: dummy_u, dummy_j, dummy_c, dummy_a, dummy_s !ddor
      integer           :: dummy_hub_nsp !ddor
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      integer, allocatable :: ibuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: integers
      integer :: row, ii
      ! ars: logical
      logical :: element_found



      ! pdh: allocate buffers
      ! ddor: Altered for DFT+U
      allocate(dbuf(max(nat*5,9)),stat=ierr)
      call utils_alloc_check('internal_hubbard','dbuf',ierr)
      allocate(ibuf(nat*4),stat=ierr)
      call utils_alloc_check('internal_hubbard','ibuf',ierr)
      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_hubbard','cbuf',ierr)

      ! ddor: Root proc only: read DFT+U information from <hubbard> block of
      ! ddor: input file
      hub_nsp = 0 ! Number of Hubbard species in cell
      pub_hubbard = .false.
      element_found = .false. ! error check flag
      if (pub_on_root) then
         if (esdf_block('hubbard',dummy_hub_nsp)) then
            if ( (dummy_hub_nsp .le. nsp) .and. ( dummy_hub_nsp > 0 ) ) then
               ! aam: check for multiply defined species
               call internal_check_species(dummy_hub_nsp,'<hubbard>')
               hub_nsp = dummy_hub_nsp ! The number of Hubbard species
               pub_hubbard = .true. ! Turn on DFT+U if Hubbard atom found
               do row=1,dummy_hub_nsp
                  read(block_data(row),*) dummy_id, &
                       ibuf(row),dbuf(row*5-4:row*5)
                  cbuf(row)(1:4) = dummy_id
               end do
            else
               call utils_abort('Error in internal_hubbard(): mismatching &
                    &number of species in <hubbard> block of input file.')
            end if
         else
            if (pub_hubbard) then
               call utils_abort('Error in internal_hubbard(): &
                    &<hubbard> block not found in input file, yet the DFT+U &
                    &feature is activated.')
            endif
         end if
      end if

      ! ddor: broadcast flag for DFT+U calculation
      call comms_bcast(pub_root_proc_id,pub_hubbard)
      ! ddor: broadcast number of Hubbard species
      call comms_bcast(pub_root_proc_id,hub_nsp)
      par%num_hub_species = hub_nsp

      if (pub_hubbard) then

         ! ddor: broadcast species information
         call comms_bcast(pub_root_proc_id,ibuf,hub_nsp)
         call comms_bcast(pub_root_proc_id,dbuf,hub_nsp*5)
         do row=1,hub_nsp
            call comms_bcast(pub_root_proc_id,cbuf(row),4)
         end do

         ! ddor: allocate public DFT+U h_species type
         call hubbard_init_species(par%num_hub_species)

         ! ddor: extract Hubbard element information to h_species
         dummy_h = 0 ! initialise counter of Hubbard atoms
         num_hub_proj = 0 ! initialise counter of Hubbard projectors
         ! Set some defaults for non-Hubbard atoms
         do row=1,hub_nsp
            dummy_id = cbuf(row)(1:4)
            dummy_l = ibuf(row)
            dummy_u = dbuf(row*5-4)
            dummy_j = dbuf(row*5-3)
            dummy_c = dbuf(row*5-2)
            dummy_a = dbuf(row*5-1)
            dummy_s = dbuf(row*5)
            do ii=1,nat ! The total number of atoms in the cell
               if (elements(ii)%species_id == dummy_id) then
                  ! If it's a Hubbard atom then...
                  ! Store namelist of Hubbard atoms
                  h_species(row)%hub_species = dummy_id(1:4)
                  ! Angular momentum channel of projector
                  h_species(row)%hub_ang_mom = dummy_l
                  ! Hubbard U parameter (eV)
                  h_species(row)%hub_u       = dummy_u / HARTREE_IN_EVS
                  ! Hund's exchange J parameter (eV)
                  h_species(row)%hub_j       = dummy_j / HARTREE_IN_EVS
                  ! Effective charge for radial function
                  h_species(row)%hub_charge  = ABS(dummy_c)
                  ! Hubbard alpha parameter (eV)
                  h_species(row)%hub_alpha   = dummy_a / HARTREE_IN_EVS
                  ! Hubbard spin Zeeman splitting (eV)
                  h_species(row)%hub_spin_splitting = dummy_s / HARTREE_IN_EVS
                  h_species(row)%orig_hub_spin_splitting = dummy_s / HARTREE_IN_EVS
                  ! ddor: Trigger fireball NGWFs by using negative Z parameter
                  if (( dummy_c < 0.0_DP ).AND.&
                       & (.not.pub_hub_read_projectors)) then
                     pub_hubbard_atomsolve = .true.
                  endif

                  !gom: for cdft_multi_proj set number cDFT-proj=number
                  !val_NGWFs
                  if (pub_cdft_multi_proj) &
                       h_species(row)%cdft_num_proj  = elements(ii)%nfunctions

                  if ( ( h_species(row)%hub_ang_mom < 0 ) .or. &
                       &( h_species(row)%hub_ang_mom > 3 ) ) then
                     call utils_abort('Error in internal_hubbard(): unphysical &
                          &value of angular momentum for projector in &
                          &<hubbard> block.')
                  endif

                  ! ddor: count total number of Hubbard atoms
                  dummy_h = dummy_h + 1
                  !gom: change projector counting scheme for cdft_multi_proj
                  if (pub_cdft_multi_proj) then
                     num_hub_proj = num_hub_proj + h_species(row)%cdft_num_proj
                  else
                     num_hub_proj = num_hub_proj + 2 * h_species(row)%hub_ang_mom +1
                  endif
                  element_found = .true.

               end if
            end do

            ! gibo: ============= DFT+U PROPERTIES-friendly
            ! gibo: spare the user the need to set (i) pub_task=PROPERTIES,
            ! (ii) hubbard_proj_mixing <0, (iii) Z (constrained_dft block) >0
            ! for properties run
            if ( (pub_task=='PROPERTIES') .OR. (pub_task=='COND') .OR.      &
                 (pub_task=='PROPERTIES_COND') .OR. (pub_task=='TDDFT')) then
               pub_hubbard_atomsolve = .FALSE.
               pub_hubbard_restart   = .TRUE.
               pub_hub_proj_mixing   = -1.0_DP
               if (pub_on_root) write(stdout,'(a)') 'WARNING in internal_hubbard:&
                    &pub_hubbard_restart overriden to T'
               if (pub_on_root) write(stdout,'(a)') 'WARNING in internal_hubbard:&
                    &pub_hub_atomsolve overriden to F'
               if (pub_on_root) write(stdout,'(a)') 'WARNING in internal_hubbard:&
                    &pub_hub_proj_mixing overriden to "-1.0_DP"'
            endif
            ! gibo: ============= DFT+U PROPERTIES-friendly

            if (pub_hubbard_restart .and. pub_hubbard_atomsolve) then
               call utils_abort('Error in internal_hubbard(): Cannot have &
                    &negative values for DFT+U charge Z with PAW or &
                    &negative projector mixing parameter in &
                    &<hubbard> block.')
            endif

            ! gibo: ============= DFT+U PROPERTIES-friendly
            call comms_bcast(pub_root_proc_id, pub_hubbard_restart)
            call comms_bcast(pub_root_proc_id, pub_hub_proj_mixing)
            ! gibo: ============= DFT+U PROPERTIES-friendly
            call comms_bcast(pub_root_proc_id, pub_hubbard_atomsolve)
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_species)
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_ang_mom)
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_u)
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_j)
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_charge)
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_alpha)
            call comms_bcast(pub_root_proc_id, &
                 & h_species(row)%hub_spin_splitting)
            call comms_bcast(pub_root_proc_id, &
                 & h_species(row)%orig_hub_spin_splitting)
         end do

         if (.not. element_found) then
            call utils_abort('Error in internal_hubbard(): &
                 &mismatching species in <hubbard> and <positions_abs> &
                 &blocks of input file.')
         else
            par%nat_hub = dummy_h ! The number of Hubbard atoms in the cell
            par%num_hub_proj = num_hub_proj
         end if

         ! ddor: broadcast number of Hubbard atoms in the cell
         call comms_bcast(pub_root_proc_id,par%nat_hub)
         call comms_bcast(pub_root_proc_id,par%num_hub_proj)

         if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Read <hubbard> block'

      end if

      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_hubbard','cbuf',ierr)
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_hubbard','dbuf',ierr)
      deallocate(ibuf,stat=ierr)
      call utils_dealloc_check('internal_hubbard','ibuf',ierr)

    end subroutine internal_hubbard


    !------------------------
    !--------- %block constrained_dft
    !------------------------

    ! IDEA: to (internally) parse the constrained_dft block as a modified version of
    ! hubbard_init (written by ddor)
    !
    subroutine internal_cdft()

      use comms, only: pub_root_proc_id, comms_bcast
      use constants, only: HARTREE_IN_EVS, CRLF
      use esdf, only: block_data, esdf_block
      use hubbard_init, only: hubbard_init_species, h_species
      use rundat, only: pub_hubbard, pub_hubbard_atomsolve, &
           pub_cdft,pub_cdft_atom_charge, pub_cdft_atom_spin, &
           pub_cdft_group_charge_acceptor, pub_cdft_group_charge_donor, &
           pub_cdft_group_spin_acceptor, pub_cdft_group_spin_donor, &
           pub_cdft_group_charge_acceptor_u, pub_cdft_group_charge_donor_u, &
           pub_cdft_group_spin_acceptor_u, pub_cdft_group_spin_donor_u, &
           pub_cdft_group_charge_diff, pub_cdft_group_spin_diff, &
           pub_cdft_group_charge_diff_u, pub_cdft_group_spin_diff_u, &
           pub_spin_polarised, pub_cdft_guru, pub_cdft_hubbard, &
           pub_hubbard_restart, pub_hub_proj_mixing, pub_task, &
           pub_cdft_read_projectors, pub_cdft_multi_proj
      use utils, only: utils_abort

      implicit none

      !-------------
      ! ars: dummies
      character(len=4)  :: dummy_id
      integer           :: dummy_l, dummy_h
      real(kind=dp)     :: dummy_u, dummy_j, dummy_c
      real(kind=dp)     :: dummy_u_s, dummy_s_target
      real(kind=dp)     :: dummy_u_q_up, dummy_u_q_down, dummy_n_up, dummy_n_down
      integer           :: dummy_cdft_nsp
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      integer, allocatable :: ibuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: integers
      integer :: row, ii
      ! ars: logical
      logical :: element_found
      !gibo: smallest representable real
      real(KIND=DP) :: eps

      ! pdh: allocate buffers
      ! gibo: modified for cDFT 04.10.11
      allocate(dbuf(max(nat*9,9)),stat=ierr)
      call utils_alloc_check('internal_cdft','dbuf',ierr)
      ! gibo: modified for cDFT 04.10.11
      allocate(ibuf(nat*8),stat=ierr)
      call utils_alloc_check('internal_cdft','ibuf',ierr)
      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_cdft','cbuf',ierr)

      ! Root proc only: read DFT+U information from <constrained_dft> block of input file
      element_found = .false. ! error check flag

      cdft_nsp = 0 ! Number of cDFT(=Hubbard) species in cell
      pub_cdft = .false.

      if (pub_on_root) then
         if (esdf_block('constrained_dft',dummy_cdft_nsp)) then
            if ( (dummy_cdft_nsp .le. nsp) .and. ( dummy_cdft_nsp > 0 ) ) then
               ! aam: check for multiply defined species
               call internal_check_species(dummy_cdft_nsp,'<constrained_dft>')
               cdft_nsp = dummy_cdft_nsp ! The number of constrained species
               pub_cdft = .true. ! Turn on cDFT

               !gibo: make sure the <hubbard> and <constrained_cdft> blocks
               !      are *NOT* simultanously defined in input file
               !MIND: internal_hubbard is called before internal_cdft therefore
               !      pub_hubbard=T if <hubbard> block in .dat file...
               if ( (pub_cdft) .AND. (pub_hubbard) ) then
                  call utils_abort('Error in internal_cdft(): <hubbard> and &
                       &<constrained_dft> blocks cannot be simultanously &
                       &present in input file. Re-run with either <hubbard> or &
                       &<constrained_dft> block in input file.')
               endif

               ! *** CRUCIAL **** if (pub_cdft) then turn on also DFT+U...
               hub_nsp = cdft_nsp ! number of Hub species=number of cDFT species
               pub_hubbard = .true. ! Turn on also DFT+U
               ! *** CRUCIAL **** if (pub_cdft) then turn on also DFT+U...

               do row=1,dummy_cdft_nsp
                  read(block_data(row),*,end=100) dummy_id, &
                       ibuf(row),dbuf(row*9-8:row*9)
                  cbuf(row)(1:4) = dummy_id
               end do
            else
               call utils_abort('Error in internal_cdft(): mismatching number &
                    &of species in <constrained_dft> block of input file.')
            end if
         else
            if (pub_cdft) then
               call utils_abort('Error in internal_cdft(): &
                    &<constrained_dft> block not found in input file, yet the &
                    &cDFT feature is activated [pub_cdft=.TRUE.].')
            endif
         end if
      end if

      !gibo: broadcast flag for cDFT calculation
      call comms_bcast(pub_root_proc_id,pub_cdft)

      !gibo: to avoid messing with pure Hubbard calculations...
      if (pub_cdft) then
         ! ddor: broadcast flag for DFT+U calculation
         call comms_bcast(pub_root_proc_id,pub_hubbard)
         ! ddor: broadcast number of Hubbard species
         call comms_bcast(pub_root_proc_id,hub_nsp)
         par%num_hub_species  = hub_nsp
         call comms_bcast(pub_root_proc_id,cdft_nsp)
         par%num_cdft_species = cdft_nsp
      endif

      ! gibo: follow ddor's DFT+U framework as much as possible...
      CHECK_CDFT: if (pub_cdft) then

         ! If required, prepare to (hopefully) help unexperienced users
         ! [cdft_guru = .FALSE.]
         if (.not.pub_cdft_guru) then
            eps = epsilon(1.0_DP)
            if (pub_on_root) then
               write(stdout,'(/a)') 'WARNING in &
                    &internal_cdft: '
               write(stdout,'(a)') 'Fail-safe initialisation of &
                    & cDFT (non-zero) |U-potentials| to 1 eV'
               write(stdout,'(a)') 'If unhappy about this, &
                    &and you really know what you are doing, &
                    &set "cdft_guru: T" in input file'
            end if
         endif

         ! gibo: turn spin-polarisation on for constrained-DFT runs
         pub_spin_polarised = .true.
         call comms_bcast(pub_root_proc_id, pub_spin_polarised)

         if (pub_on_root) write(stdout,'(a)') &
              'WARNING in internal_cdft: turning spin-polarisation on for &
              &constrained-DFT run'

         ! ddor: broadcast species information
         call comms_bcast(pub_root_proc_id,ibuf,hub_nsp)
         call comms_bcast(pub_root_proc_id,dbuf,hub_nsp*9)  !gibo: modified for cDFT 04.10.11
         do row=1,hub_nsp
            call comms_bcast(pub_root_proc_id,cbuf(row),8)  !gibo: modified for cDFT 04.10.11
         end do

         ! ddor: allocate public DFT+U h_species type
         call hubbard_init_species(par%num_hub_species)

         ! ddor: extract Hubbard element information to h_species
         dummy_h = 0 ! initialise counter of Hubbard atoms
         num_hub_proj = 0 ! initialise counter of Hubbard projectors

         !gibo: initialise counters of cDFT_charge/spin acceptor/donor
         par%nat_hub_charge_acceptor = 0
         par%nat_hub_charge_donor    = 0
         par%nat_hub_spin_acceptor   = 0
         par%nat_hub_spin_donor      = 0

         !gibo: initialise to zero number of cdft_active atoms (needed for cdft_hubbard:T)
         par%nat_cdft_active     = 0

         ! Set some defaults for non-Hubbard atoms
         do row=1,hub_nsp
            dummy_id = cbuf(row)(1:4)
            dummy_l = ibuf(row)

            ! different U_q for alpha and beta spin-channels 04.10.11
            dummy_c         = dbuf(row*9-8)
            dummy_u         = dbuf(row*9-7)
            dummy_j         = dbuf(row*9-6)
            dummy_u_q_up    = dbuf(row*9-5)
            dummy_u_q_down  = dbuf(row*9-4)
            dummy_u_s       = dbuf(row*9-3)
            dummy_n_up      = dbuf(row*9-2)
            dummy_n_down    = dbuf(row*9-1)
            dummy_s_target  = dbuf(row*9)
            ! different U_q for alpha and beta spin-channels 04.10.11

            ! Unless the user claims to be a cDFT guru [pub_cdft_guru = .TRUE.]
            ! [who should know about the best cDFT U-potentials initialisation]
            ! initialise the (non-zero) cDFT U-potentials to 1 eV
            if (.not.pub_cdft_guru) then
               if (ABS(dummy_u_q_up) > eps) &
                    dummy_u_q_up = dummy_u_q_up/ABS(dummy_u_q_up)
               if (ABS(dummy_u_q_down) > eps) &
                    dummy_u_q_down = dummy_u_q_down/ABS(dummy_u_q_down)
               if (ABS(dummy_u_s) > eps) &
                    dummy_u_s = dummy_u_s/ABS(dummy_u_s)
            endif

            ATOMS_CELL: do ii=1,nat ! The total number of atoms in the cell
               if (elements(ii)%species_id == dummy_id) then
                  ! gibo: If it is a cDFT (=Hubbard) atom then...
                  ! Store namelist of Hubbard atoms
                  h_species(row)%hub_species = dummy_id(1:4)

                  ! modified ==== start
                  ! Angular momentum channel of projector
                  h_species(row)%cdft_ang_mom        = dummy_l
                  ! Copy cdft angular momentum into DFT-U counterpart
                  h_species(row)%hub_ang_mom         = dummy_l

                  ! Effective charge for radial function
                  h_species(row)%hub_charge          = ABS(dummy_c)

                  ! Hubbard U parameter (eV)
                  h_species(row)%hub_u               = dummy_u / HARTREE_IN_EVS

                  ! Hund's exchange J parameter (eV)
                  h_species(row)%hub_j               = dummy_j / HARTREE_IN_EVS

                  ! different U_q for alpha and beta spin-channels 04.10.11
                  h_species(row)%cdft_u_charge_up    = dummy_u_q_up / HARTREE_IN_EVS

                  h_species(row)%cdft_u_charge_down  = dummy_u_q_down / HARTREE_IN_EVS
                  ! different U_q for alpha and beta spin-channels 04.10.11

                  ! SPIN constrain (U) parameter (eV)
                  h_species(row)%cdft_u_spin         = dummy_u_s / HARTREE_IN_EVS

                  ! TARGET: number of electrons_UP
                  h_species(row)%cdft_target_up      = dummy_n_up

                  ! TARGET: number of electrons_DOWN
                  h_species(row)%cdft_target_down    = dummy_n_down

                  ! TARGET: spin population (N_alpha-N_beta)
                  h_species(row)%cdft_target_spin    = dummy_s_target

                  ! pub_hubbard_atomsolve = .true.
                  ! ddor: Trigger fireball NGWFs by using negative Z parameter
                  if (( dummy_c < 0.0_DP ).AND.&
                       & (.not.pub_cdft_read_projectors)) then
                     pub_hubbard_atomsolve = .true.
                  endif

                  !gibo: for cdft_multi_proj set number cDFT-proj=number val_NGWFs
                  if (pub_cdft_multi_proj) &
                       h_species(row)%cdft_num_proj  = elements(ii)%nfunctions
                  ! modified ==== end

                  if ( ( h_species(row)%hub_ang_mom < 0 ) .and. &
                       &( h_species(row)%hub_ang_mom > 3 ) ) then
                     call utils_abort('Error in internal_cdft(): unphysical &
                          &value of angular momentum for projector in &
                          &<constrained_dft> block.')
                  endif

                  ! make sure that ABS(h_species(row)%hub_u) > 0. for cdft_hubbard=.T.
                  if ((.not.pub_cdft_hubbard) .AND. &
                       (ABS(h_species(row)%hub_u) > 0.0_DP)) then
                     call utils_abort('Error in internal_cdft(): set "cdft_&
                          &"hubbard : T" for a constrained-DFT+U simulation.')
                  endif

                  ! make sure that ABS(h_species(row)%hub_j) > 0. for cdft_hubbard=.T.
                  if ((.not.pub_cdft_hubbard) .AND. &
                       (ABS(h_species(row)%hub_j) > 0.0_DP)) then
                     call utils_abort('Error in internal_cdft(): set "cdft_&
                          &"hubbard : T" for a constrained-DFT+U+J simulation.')
                  endif

                  ! increase number of cDFT-active atoms (U_hub >= 0 )
                  if ((pub_cdft_hubbard).AND.(h_species(row)%hub_u >= 0.0_DP)) then
                     par%nat_cdft_active =  par%nat_cdft_active + 1
                     h_species(row)%cdft_active = .TRUE.
                  elseif ((pub_cdft_hubbard).AND.(h_species(row)%hub_u < 0.0_DP)) then
                     ! U_hub < 0: this is a DFT+U_only atom
                     h_species(row)%cdft_active = .FALSE.
                  else
                     ! for cDFT-only runs, all species (in the cDFT-block) are cDFT-active
                     h_species(row)%cdft_active = .TRUE.
                     par%nat_cdft_active =  par%nat_cdft_active + 1
                  endif

                  ! make sure DFT+U_only atoms have zero cDFT-potentials/targets
                  if ((pub_cdft_hubbard).AND.(.not.h_species(row)%cdft_active)) then
                     h_species(row)%cdft_u_charge_up   = 0.0_DP
                     h_species(row)%cdft_u_charge_down = 0.0_DP
                     h_species(row)%cdft_u_spin        = 0.0_DP
                     h_species(row)%cdft_target_up     = 0.0_DP
                     h_species(row)%cdft_target_down   = 0.0_DP
                     h_species(row)%cdft_target_spin   = 0.0_DP
                  endif

                  !gibo: for cDFT+U runs, avoid checks on groups for cdft-INactive atoms
                  SKIP_cDFT_INACTIVE: if (h_species(row)%cdft_active) then

                     ! If constraining group CHARGE difference,
                     ! decide whether a given atom is an acceptor or a donor...
                     ! dhpt: initialise defaults
                     h_species(row)%cdft_charge_donor    = .FALSE.
                     h_species(row)%cdft_charge_acceptor = .FALSE.

                     if (pub_cdft_group_charge_acceptor .OR. &
                          pub_cdft_group_charge_donor    .OR. &
                          pub_cdft_group_charge_diff) then

                        if (dummy_u_q_up > 0._DP) then
                           h_species(row)%cdft_charge_donor    = .TRUE.
                           h_species(row)%cdft_charge_acceptor = .FALSE.
                           par%nat_hub_charge_donor = &
                                par%nat_hub_charge_donor + 1

                        elseif (dummy_u_q_up < 0._DP) then
                           h_species(row)%cdft_charge_donor    = .FALSE.
                           h_species(row)%cdft_charge_acceptor = .TRUE.
                           par%nat_hub_charge_acceptor = &
                                par%nat_hub_charge_acceptor + 1

                        endif

                        if(pub_debug_on_root) write(stdout,'(a,f10.5,a,l2)') &
                             'dummy_u_q_up = ',&
                             &dummy_u_q_up, ', cdft_charge_donor = ',&
                             &h_species(row)%cdft_charge_donor

                     endif

                     ! If constraining group CHARGE (difference),
                     ! avoid silly inputs like (U_up > 0, U_down <0 and
                     ! cdft_group_charge/spin_[diff] active)
                     if (pub_cdft_group_charge_acceptor .OR. &
                          pub_cdft_group_charge_donor    .OR. &
                          pub_cdft_group_charge_diff) then

                        if ( dummy_u_q_up*dummy_u_q_down < 0._DP) then
                           call utils_abort('Error in internal_cdft(): &
                                &cdft_group_charge_acceptor/donor active, &
                                &yet different signs for U_up and U_down of the &
                                &same atom in <constrained_dft> block.')
                        endif
                     endif

                     ! If constraining group SPIN difference,
                     ! decide whether a given atom is an acceptor or a donor...
                     ! dhpt: initialise defaults
                     h_species(row)%cdft_spin_donor    = .FALSE.
                     h_species(row)%cdft_spin_acceptor = .FALSE.

                     if (pub_cdft_group_spin_acceptor .OR. &
                          pub_cdft_group_spin_donor    .OR. &
                          pub_cdft_group_spin_diff) then

                        if (dummy_u_s > 0._DP) then
                           h_species(row)%cdft_spin_donor    = .TRUE.
                           h_species(row)%cdft_spin_acceptor = .FALSE.
                           par%nat_hub_spin_donor = &
                                par%nat_hub_spin_donor + 1

                        elseif (dummy_u_s < 0._DP) then
                           h_species(row)%cdft_spin_donor    = .FALSE.
                           h_species(row)%cdft_spin_acceptor = .TRUE.
                           par%nat_hub_spin_acceptor = &
                                par%nat_hub_spin_acceptor + 1

                        endif

                     endif

                  endif SKIP_cDFT_INACTIVE

                  ! ddor: count total number of Hubbard atoms
                  ! gibo: number of Hubbard atoms = number of cDFT atoms
                  dummy_h = dummy_h + 1

                  !gibo: change projector counting scheme for cdft_multi_proj
                  if (pub_cdft_multi_proj) then
                     num_hub_proj = num_hub_proj + h_species(row)%cdft_num_proj
                  else
                     num_hub_proj = num_hub_proj + 2 * h_species(row)%hub_ang_mom + 1
                  endif

                  element_found = .true.

               end if

            enddo ATOMS_CELL

            ! gibo: ============= cDFT-PROPERTIES-friendly
            ! gibo: spare the user the need to set (i) pub_task=PROPERTIES,
            ! (ii) hubbard_proj_mixing <0, (iii) Z (constrained_dft block) >0
            ! for properties run
            if ( (pub_task=='PROPERTIES') .OR. (pub_task=='COND') .OR.       &
                 (pub_task=='PROPERTIES_COND') .OR. (pub_task=='TDDFT')) then
               pub_hubbard_atomsolve = .FALSE.
               pub_hubbard_restart   = .TRUE.
               pub_hub_proj_mixing   = -1.0_DP
               if (pub_on_root) write(stdout,'(a)') 'WARNING in internal_cdft: &
                    &pub_hubbard_restart overriden to T'
               if (pub_on_root) write(stdout,'(a)') 'WARNING in internal_cdft: &
                    &pub_hub_atomsolve overriden to F'
               if (pub_on_root) write(stdout,'(a)') 'WARNING in internal_cdft: &
                    &pub_hub_proj_mixing overriden to "-1.0_DP"'
            endif
            ! gibo: ============= cDFT-PROPERTIES-friendly

            ! 4 GIBO: TO BE CHECKED AFTER PAW GO AHEAD FROM NICK ==== START
            !: cDFT within PAW is something we fancy...
            if (pub_hubbard_restart .and. pub_hubbard_atomsolve) then
               call utils_abort('Error in internal_cdft(): Cannot have &
                    &negative values for DFT+U charge Z with PAW or &
                    &negative projector mixing parameter in &
                    &<hubbard> block.')
            endif
            ! 4 GIBO: TO BE CHECKED AFTER PAW GO AHEAD FROM NICK ==== END

            ! gibo: cDFT-PROPERTIES-friendly
            call comms_bcast(pub_root_proc_id, pub_hubbard_restart)
            call comms_bcast(pub_root_proc_id, pub_hub_proj_mixing)
            ! gibo: cDFT-PROPERTIES-friendly
            call comms_bcast(pub_root_proc_id, pub_hubbard_atomsolve)
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_species)
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_ang_mom)
            !gibo: make sure U_hub > 0 for (c)DFT-U_only species
            if ((pub_cdft_hubbard).AND.(h_species(row)%hub_u < 0.)) then
               h_species(row)%hub_u = -h_species(row)%hub_u
            endif
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_u)
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_j)
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_charge)
            call comms_bcast(pub_root_proc_id, h_species(row)%cdft_num_proj)
            call comms_bcast(pub_root_proc_id, h_species(row)%cdft_u_charge_up)
            call comms_bcast(pub_root_proc_id, h_species(row)%cdft_u_charge_down)
            call comms_bcast(pub_root_proc_id, h_species(row)%cdft_u_spin)
            call comms_bcast(pub_root_proc_id, h_species(row)%cdft_target_up)
            call comms_bcast(pub_root_proc_id, h_species(row)%cdft_target_down)
            call comms_bcast(pub_root_proc_id, h_species(row)%cdft_target_spin)
            call comms_bcast(pub_root_proc_id, h_species(row)%cdft_charge_donor)
            call comms_bcast(pub_root_proc_id, h_species(row)%cdft_charge_acceptor)
            call comms_bcast(pub_root_proc_id, h_species(row)%cdft_spin_donor)
            call comms_bcast(pub_root_proc_id, h_species(row)%cdft_spin_acceptor)

            !gibo: to avoid problems, initialise (unused) Hubbard variables
            !      and send them around
            h_species(row)%hub_alpha = 0._DP
            h_species(row)%hub_spin_splitting = 0._DP
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_alpha)
            call comms_bcast(pub_root_proc_id, &
                 & h_species(row)%hub_spin_splitting)

         end do

         ! gibo: if <constrained_dft> block is present
         ![we are inside the CHECK_CDFT IF-THEN construct]
         ! make sure at least one cDFT-mode is activated
         if (.not.(pub_cdft_atom_charge .OR. pub_cdft_atom_spin .OR. &
              pub_cdft_group_charge_acceptor .OR. &
              pub_cdft_group_charge_donor    .OR. &
              pub_cdft_group_spin_acceptor   .OR. &
              pub_cdft_group_spin_donor      .OR. &
              pub_cdft_group_charge_diff     .OR. &
              pub_cdft_group_spin_diff) ) then
            call utils_abort('Error in internal_cdft(): no cDFT mode selected, &
                 &despite the presence of  <constrained_dft> block. &
                 &Select *ONE* from available cDFT-modes: &
                 &[1] cdft_atom_charge           : T, &
                 &[2] cdft_atom_spin             : T, &
                 &[3] cdft_group_charge_acceptor : T, &
                 &[4] cdft_group_charge_donor    : T, &
                 &[5] cdft_group_spin_acceptor   : T, &
                 &[6] cdft_group_spin_donor      : T, &
                 &[7] cdft_group_charge_diff     : T, &
                 &[8] cdft_group_spin_diff       : T.')
         endif

         !gibo: for group-charge_acceptor/donor runs, make sure all the acceptor and donor
         !      atoms have the same constraining potentials
         if (pub_cdft_group_charge_acceptor .OR. &
              pub_cdft_group_charge_donor    .OR. &
              pub_cdft_group_charge_diff) then

            LOOP_OUT_1: do ii=1,hub_nsp

               !gibo: avoid checking for cDFT-INactive atoms
               if ((pub_cdft_hubbard).AND.(.not.h_species(ii)%cdft_active)) CYCLE LOOP_OUT_1

               ! allow single point runs cDFT with 1-atom groups...
               LOOP_IN_1: do row = ii, hub_nsp

                  !gibo: avoid checking for cDFT-INactive atoms
                  if ((pub_cdft_hubbard).AND.(.not.h_species(row)%cdft_active)) CYCLE LOOP_IN_1

                  if (h_species(ii)%cdft_charge_donor .AND. h_species(row)%cdft_charge_donor ) then

                     ! to avoid another loop...
                     if (.not.pub_cdft_group_charge_diff) &
                          pub_cdft_group_charge_donor_u = h_species(ii)%cdft_u_charge_up

                     if ( (h_species(ii)%cdft_u_charge_up   .NE. &
                          h_species(row)%cdft_u_charge_up) .OR. &
                          (h_species(ii)%cdft_u_charge_down .NE. &
                          h_species(row)%cdft_u_charge_down) ) then

                        call utils_abort('Error in internal_cdft(): &
                             &mismatching U_charge values in <constrained_dft> block&
                             & for pub_cdft_group_charge_donor simulation [DONOR].')
                     endif

                  endif

                  if (h_species(ii)%cdft_charge_acceptor .AND. h_species(row)%cdft_charge_acceptor ) then

                     ! to avoid another loop...
                     if (.not.pub_cdft_group_charge_diff) &
                          pub_cdft_group_charge_acceptor_u =  h_species(ii)%cdft_u_charge_up

                     if ( (h_species(ii)%cdft_u_charge_up   .NE. &
                          h_species(row)%cdft_u_charge_up) .OR. &
                          (h_species(ii)%cdft_u_charge_down .NE. &
                          h_species(row)%cdft_u_charge_down) ) then

                        call utils_abort('Error in internal_cdft(): &
                             &mismatching U_charge values in <constrained_dft> block &
                             &for pub_cdft_group_charge_acceptor simulation [ACCEPTOR].')
                     endif

                  endif

               enddo LOOP_IN_1
            enddo LOOP_OUT_1

            !gibo: for group_charge_diff runs, make sure the absolute value of
            !      the cDFT-potential of donor- and acceptor- atoms is the same
            if (pub_cdft_group_charge_diff) then

               LOOP_OUT_2: do ii=1,hub_nsp

                  !gibo: avoid checking for cDFT-INactive atoms
                  if ((pub_cdft_hubbard).AND.(.not.h_species(ii)%cdft_active)) CYCLE LOOP_OUT_2

                  LOOP_IN_2: do row = ii+1, hub_nsp

                     !gibo: avoid checking for cDFT-INactive atoms
                     if ((pub_cdft_hubbard).AND.(.not.h_species(row)%cdft_active)) CYCLE LOOP_IN_2

                     if (h_species(ii)%cdft_charge_acceptor .AND. &
                          h_species(row)%cdft_charge_donor ) then

                        if ( ( ABS(h_species(ii)%cdft_u_charge_up)   .NE. &
                             ABS(h_species(row)%cdft_u_charge_up)) .OR. &
                             ( ABS(h_species(ii)%cdft_u_charge_down) .NE. &
                             ABS(h_species(row)%cdft_u_charge_down) ) ) then

                           call utils_abort('Error in internal_cdft(): &
                                &mismatching U_charge values between acceptor- and &
                                &donor-atoms in <constrained_dft> block &
                                &for cdft_group_charge_diff simulation].')
                        endif
                     endif
                  enddo LOOP_IN_2
               enddo LOOP_OUT_2
            endif

            ! cdft_group_charge_diff mode
            if (pub_cdft_group_charge_diff) then
               !gibo: since all acceptor and donor atoms have the same constraining
               !      potentials, set pub_cdft_charge_group_diff_u
               !      [from 1st hub_species], and broadcast as needed in
               !      sbrtne cdft_energy_total [hubbard_build_mod.F90]
               !gibo: use ABS and act accordingly in 'sbrtne cdft_energy_total'
               ![alternative: re-loop over hub_nsp and single out acc. and don. species]
               pub_cdft_group_charge_diff_u = ABS( h_species(1)%cdft_u_charge_up )
               call comms_bcast(pub_root_proc_id,pub_cdft_group_charge_diff_u)

            else
               call comms_bcast(pub_root_proc_id,pub_cdft_group_charge_donor_u)
               call comms_bcast(pub_root_proc_id,pub_cdft_group_charge_acceptor_u)
            endif
         endif


         !gibo: for group-spin_diff runs,
         !      make sure that all the acceptor and donor atoms have the same
         !      constraining potentials
         if (pub_cdft_group_spin_acceptor .OR. &
              pub_cdft_group_spin_donor    .OR. &
              pub_cdft_group_spin_diff) then

            LOOP_OUT_3: do ii=1,hub_nsp

               !gibo: avoid checking for cDFT-INactive atoms
               if ((pub_cdft_hubbard).AND.(.not.h_species(ii)%cdft_active)) CYCLE LOOP_OUT_3

               ! allow single point runs cDFT with 1-atom groups...
               LOOP_IN_3: do row = ii, hub_nsp

                  !gibo: avoid checking for cDFT-INactive atoms
                  if ((pub_cdft_hubbard).AND.(.not.h_species(row)%cdft_active)) CYCLE LOOP_IN_3

                  if (h_species(ii)%cdft_spin_donor .AND. h_species(row)%cdft_spin_donor ) then

                     ! to avoid another loop...
                     if (.not.pub_cdft_group_spin_diff) &
                          pub_cdft_group_spin_donor_u =  h_species(ii)%cdft_u_spin

                     if (h_species(ii)%cdft_u_spin .NE. h_species(row)%cdft_u_spin) then

                        call utils_abort('Error in internal_cdft(): &
                             &mismatching U_spin values in <constrained_dft> block for &
                             &pub_cdft_group_spin_diff simulation [DONOR].')
                     endif

                  endif

                  if (h_species(ii)%cdft_spin_acceptor .AND. h_species(row)%cdft_spin_acceptor ) then

                     ! to avoid another loop...
                     if (.not.pub_cdft_group_spin_diff) &
                          pub_cdft_group_spin_acceptor_u =  h_species(ii)%cdft_u_spin

                     if (h_species(ii)%cdft_u_spin .NE. h_species(row)%cdft_u_spin) then

                        call utils_abort('Error in internal_cdft(): &
                             &mismatching U_spin values in <constrained_dft> block for &
                             &pub_cdft_group_spin_diff simulation [ACCEPTOR].')

                     endif

                  endif

               enddo LOOP_IN_3
            enddo LOOP_OUT_3

            !gibo: for group_spin_diff runs, make sure the absolute value of
            !      the cDFT-potential of donor- and acceptor- atoms is the same
            if (pub_cdft_group_spin_diff) then

               LOOP_OUT_4: do ii=1,hub_nsp

                  !gibo: avoid checking for cDFT-INactive atoms
                  if ((pub_cdft_hubbard).AND.(.not.h_species(ii)%cdft_active)) CYCLE LOOP_OUT_4

                  LOOP_IN_4: do row = ii+1, hub_nsp

                     !gibo: avoid checking for cDFT-INactive atoms
                     if ((pub_cdft_hubbard).AND.(.not.h_species(row)%cdft_active)) CYCLE LOOP_IN_4

                     if (h_species(ii)%cdft_spin_acceptor .AND. &
                          h_species(row)%cdft_spin_donor ) then

                        if ( ABS(h_species(ii)%cdft_u_spin)   .NE. &
                             ABS(h_species(row)%cdft_u_spin)) then

                           call utils_abort('Error in internal_cdft(): &
                                &mismatching U_spin values between acceptor- and &
                                &donor-atoms in <constrained_dft> block &
                                &for cdft_group_spin_diff simulation].')
                        endif
                     endif
                  enddo LOOP_IN_4
               enddo LOOP_OUT_4
            endif

            ! cdft_group_spin_diff mode
            if (pub_cdft_group_spin_diff) then
               !gibo: since all acceptor and donor atoms have the same constraining
               !      potentials, set pub_cdft_charge_group_diff_u
               !      [from 1st hub_species], and broadcast as needed in
               !      sbrtne cdft_energy_total [hubbard_build_mod.F90]
               !gibo: use ABS and act accordingly in 'sbrtne cdft_energy_total'
               ![alternative: re-loop over hub_nsp and single out acc. and don. species]
               pub_cdft_group_spin_diff_u = ABS( h_species(1)%cdft_u_spin )
               call comms_bcast(pub_root_proc_id,pub_cdft_group_spin_diff_u)
            else
               call comms_bcast(pub_root_proc_id,pub_cdft_group_spin_donor_u)
               call comms_bcast(pub_root_proc_id,pub_cdft_group_spin_acceptor_u)
            endif
         endif


         if (.not. element_found) then
            call utils_abort('Error in internal_cdft(): &
                 &mismatching species in <constrained_dft> and <positions_abs> &
                 &blocks of input file.')
         else
            par%nat_hub = dummy_h ! The number of Hubbard atoms in the cell
            par%num_hub_proj = num_hub_proj
         end if

         if ((pub_cdft_hubbard) .AND. &
              (par%nat_cdft_active > par%nat_hub)) then
            call utils_abort('Error in internal_cdft(): &
                 &number of cDFT-active atoms LARGER than number of &
                 &U_hub-corrected ones.')
         endif

         ! ddor: broadcast number of Hubbard atoms in the cell
         call comms_bcast(pub_root_proc_id,par%nat_hub)
         call comms_bcast(pub_root_proc_id,par%num_hub_proj)

         !gibo: and number of cDFT charge/spin acceptor/donor atoms in the cell
         call comms_bcast(pub_root_proc_id,par%nat_cdft_active)
         call comms_bcast(pub_root_proc_id,par%nat_hub_charge_acceptor)
         call comms_bcast(pub_root_proc_id,par%nat_hub_charge_donor)
         call comms_bcast(pub_root_proc_id,par%nat_hub_spin_acceptor)
         call comms_bcast(pub_root_proc_id,par%nat_hub_spin_donor)


         ! Check safety of relevant modes vs numbers of atoms
         if (pub_cdft_group_spin_diff) then
            if (par%nat_hub_spin_acceptor<1) call utils_abort( &
                 'Error in internal_cdft(): cdft_group_spin_diff specified, &
                 &yet no spin'//CRLF//'acceptors were specified.')
            if (par%nat_hub_spin_donor<1) call utils_abort( &
                 'Error in internal_cdft(): cdft_group_spin_diff specified, &
                 &yet no spin'//CRLF//'donors were specified.')
         end if

         if (pub_cdft_group_charge_diff) then
            if (par%nat_hub_charge_acceptor<1) call utils_abort( &
                 'Error in internal_cdft(): cdft_group_charge_diff specified, &
                 &yet no charge'//CRLF//'acceptors were specified.')
            if (par%nat_hub_charge_donor<1) call utils_abort( &
                 'Error in internal_cdft(): cdft_group_charge_diff specified, &
                 &yet no charge'//CRLF//'donors were specified.')
         end if

         !##############################################################################
         !
         !   4 Nick: Would it be possible to keep the (now) commented TWEAK_OPT in devel?
         !           I am to define fail-safe initialisation as we get more experienced
         !           with cDFT....
         !
         !##############################################################################
         !         !***Tweak pub_cdft_cg_max and cdft_trial_length for unexperienced user
         !         TWEAK_OPT: if (.not.pub_cdft_guru) then
         !!            if ((pub_cdft_atom_spin.OR.                    &
         !!                 pub_cdft_group_spin_diff.OR.              &
         !!                 pub_cdft_group_spin_acceptor.OR.          &
         !!                 pub_cdft_group_spin_donor).AND.           &
         !!                 (.not.(pub_cdft_atom_charge.OR.           &
         !!                        pub_cdft_group_charge_diff.OR.     &
         !!                        pub_cdft_group_charge_acceptor.OR. &
         !!                        pub_cdft_group_charge_donor)) ) then
         !!
         !!                  !spin-only cDFT simulation: start easy...
         !!                  cdft_trial_length = 0.1_DP
         !!                  call comms_bcast(pub_root_proc_id,cdft_trial_length)
         !
         !            if ((pub_cdft_atom_charge.OR.             &
         !                 pub_cdft_group_charge_diff.OR.       &
         !                 pub_cdft_group_charge_acceptor.OR.   &
         !                 pub_cdft_group_charge_donor).AND.    &
         !                 (.not.(pub_cdft_atom_spin.OR.        &
         !                        pub_cdft_group_spin_diff.OR.  &
         !                        pub_cdft_group_spin_acceptor.OR.   &
         !                        pub_cdft_group_spin_donor)) ) then
         !                  !charge-only cDFT simulation: default (0.1_DP) is too low
         !                  cdft_trial_length = 5.0_DP
         !                  call comms_bcast(pub_root_proc_id,cdft_trial_length)
         !            endif
         !
         !            !***Refine pub_cdft_cg_max according to cDFT-mode and cDFT-atoms
         !            !atom-charge: twice as many U-potentials as cDFT(=Hubbard) atoms
         !            if (pub_cdft_atom_charge)  then
         !               pub_cdft_cg_max = 2*par%nat_hub
         !
         !            !atom-spin: as many U-potentials as cDFT(=Hubbard) atoms
         !            elseif (pub_cdft_atom_spin) then
         !                pub_cdft_cg_max = par%nat_hub
         !
         !            ! group_charge/spin_diff: one U-potential, set steepest descents
         !            elseif (pub_cdft_group_charge_diff.OR.pub_cdft_group_spin_diff) then
         !                pub_cdft_cg_max = 1
         !
         !            ! group_charge/spin: one U-potential for each group/mode
         !            else
         !              pub_cdft_cg_max = 0
         !              if (pub_cdft_group_charge_acceptor) &
         !                  pub_cdft_cg_max=pub_cdft_cg_max+1
         !
         !              if (pub_cdft_group_charge_donor) &
         !                  pub_cdft_cg_max=pub_cdft_cg_max+1
         !
         !              if (pub_cdft_group_spin_acceptor) &
         !                  pub_cdft_cg_max=pub_cdft_cg_max+1
         !
         !              if (pub_cdft_group_spin_donor) &
         !                  pub_cdft_cg_max=pub_cdft_cg_max+1
         !
         !            endif
         !
         !            call comms_bcast(pub_root_proc_id,pub_cdft_cg_max)
         !
         !         endif TWEAK_OPT
         !gibo-end-20.03.12

         if (pub_debug_on_root) write(stdout,'(a)') &
              'DEBUG: Read <constrained_dft> block'

      endif CHECK_CDFT

      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_cdft','cbuf',ierr)
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_cdft','dbuf',ierr)
      deallocate(ibuf,stat=ierr)
      call utils_dealloc_check('internal_cdft','ibuf',ierr)

      return

100   call utils_abort('Error in internal_cdft(): <constrained_dft> block &
           &does not appear '//CRLF//'to contain correct number of columns. &
           &Required format is: '//CRLF//'# L Z Uh J Uq(UP) Uq(DOWN) Us N(UP) &
           &N(DOWN) [N(UP)-N(DOWN)]')

    end subroutine internal_cdft


    !------------------------
    !--------- %block dft_nu
    !------------------------
    subroutine internal_dft_nu()

      use comms, only: pub_root_proc_id, comms_bcast
      use constants, only: HARTREE_IN_EVS
      use esdf, only: block_data, esdf_block
      use hubbard_init, only: hubbard_init_species, h_species
      use rundat, only: pub_hubbard, pub_hubbard_atomsolve, &
           pub_spin_polarised, &
           pub_hubbard_restart, pub_hub_proj_mixing, pub_task, &
           pub_cdft_read_projectors, pub_cdft_multi_proj, &
           pub_dft_nu
      use utils, only: utils_abort

      implicit none

      !-------------
      ! ars: dummies
      character(len=4)  :: dummy_id
      integer           :: dummy_l, dummy_h
      real(kind=dp)     :: dummy_c
      integer           :: dummy_group !gom
      real(kind=dp)     :: dummy_u1_up, dummy_u1_dn, dummy_u2_up, dummy_u2_dn
      real(kind=dp)     :: dummy_n_up, dummy_nu_up, dummy_n_dn, dummy_nu_dn
      integer           :: dummy_cdft_nsp
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      integer, allocatable :: ibuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: integers
      integer :: row, ii
      ! ars: logical
      logical :: element_found

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering <dft_nu>'

      ! gom: allocate buffers
      allocate(dbuf(max(nat*7,7)),stat=ierr)
      call utils_alloc_check('internal_dft_nu','dbuf',ierr)
      allocate(ibuf(max(nat*2,2)),stat=ierr)
      call utils_alloc_check('internal_dft_nu','ibuf',ierr)
      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_dft_nu','cbuf',ierr)

      ! Root proc only: read DFT+U information from <constrained_dft> block of input file
      element_found = .false. ! error check flag

      cdft_nsp = 0 ! Number of cDFT(=Hubbard) species in cell
      pub_dft_nu = .false.

      if (pub_on_root) then
         if (esdf_block('dft_nu',dummy_cdft_nsp)) then
            if ( (dummy_cdft_nsp .le. nsp) .and. ( dummy_cdft_nsp > 0 ) ) then

               ! aam: check for multiply defined species
               call internal_check_species(dummy_cdft_nsp,'<dft_nu>')
               cdft_nsp = dummy_cdft_nsp ! The number of constrained species
               pub_dft_nu = .true.

               !gom: make sure the <hubbard> and <dft_nu> blocks
               !      are *NOT* simultanously defined in input file
               !MIND: internal_hubbard is called before internal_dft_nu therefore
               !      pub_hubbard=T if <hubbard> block in .dat file...
               if ( (pub_dft_nu) .AND. (pub_hubbard) ) then
                  call utils_abort('Error in internal_dft_nu(): <hubbard> &
                       &and <dft_nu> blocks cannot be simultanously &
                       &present in input file. Re-run with either <hubbard> or &
                       &<constrained_dft> block in input file.')
               endif

               hub_nsp = cdft_nsp ! number of Hub species=number of cDFT species
               pub_hubbard = .true. ! Turn on also DFT+U

               !gom
               do row=1,dummy_cdft_nsp
                  read(block_data(row),*) dummy_id, &
                       ibuf(row*2-1),dbuf(row*7-6:row*7),ibuf(row*2)
                  cbuf(row)(1:4) = dummy_id
               end do

            else
               call utils_abort('Error in internal_dft_nu(): mismatching &
                    &number of species in <constrained_dft> block of input file.')
            end if
         else
            if (pub_dft_nu) then
               call utils_abort('Error in internal_dft_nu(): &
                    &<dft_nu> block not found in input file, yet the &
                    &DFT+nu feature is activated [pub_dft_nu=.TRUE.].')
            endif
         end if
      end if

      ! gom : broadcast flag for DFT+nu calculation
      call comms_bcast(pub_root_proc_id,pub_dft_nu)

      if (pub_dft_nu) then
         ! gom: broadcast flag for DFT+U calculation
         call comms_bcast(pub_root_proc_id,pub_hubbard)
         ! gom: broadcast number of Hubbard species
         call comms_bcast(pub_root_proc_id,hub_nsp)
         par%num_hub_species  = hub_nsp
         call comms_bcast(pub_root_proc_id,cdft_nsp)
         par%num_cdft_species = cdft_nsp
      endif

      CHECK_CDFT: if (pub_dft_nu) then

         ! gom: turn spin-polarisation on for constrained-DFT runs
         pub_spin_polarised = .true.
         call comms_bcast(pub_root_proc_id, pub_spin_polarised)

         if (pub_on_root) write(stdout,'(a)') &
              'WARNING in internal_cdft: turning spin-polarisation on for &
              &constrained-DFT run'

         ! gom: broadcast species information
         call comms_bcast(pub_root_proc_id,ibuf,hub_nsp)
         call comms_bcast(pub_root_proc_id,dbuf,hub_nsp*7)
         do row=1,hub_nsp
            call comms_bcast(pub_root_proc_id,cbuf(row),7)
         end do

         ! gom: allocate public DFT+U h_species type
         call hubbard_init_species(par%num_hub_species)

         ! gom: extract Hubbard element information to h_species
         dummy_h = 0 ! initialise counter of Hubbard atoms
         num_hub_proj = 0 ! initialise counter of Hubbard projectors

         do row=1,hub_nsp
            dummy_id = cbuf(row)(1:4)
            dummy_l = ibuf(row*2-1)

            !gom: read DFT+nu block elements
            dummy_c      = dbuf(row*7-6)
            dummy_u1_up  = dbuf(row*7-5)
            dummy_u1_dn  = dbuf(row*7-4)
            dummy_u2_up  = dbuf(row*7-3)
            dummy_u2_dn  = dbuf(row*7-2)
            dummy_n_up   = dbuf(row*7-1)
            dummy_n_dn   = dbuf(row*7)
            dummy_group  = ibuf(row*2)

            if(.not. pub_spin_polarised.and. &
               (dummy_u1_up .ne. dummy_u1_dn).or.&
               (dummy_u2_up .ne. dummy_u2_dn).or.&
               (dummy_n_up .ne.  dummy_n_dn)) then
                 call utils_abort('Error in internal_dft_nu(): &
                     &mismatching values in <dft_nu> block for &
                     &spin 0 system. All values must match.')
             endif

            ! nu is not read from file but determined by either
            ! N (if positive) or U1, U2 if N is negative

            ! Check N_up - if positive then nu takes fixed value
            if (dummy_n_up > 0.0_DP) then
               dummy_nu_up = dummy_n_up - floor(dummy_n_up)
               dummy_n_up = floor(dummy_n_up)
               ! if N_up negative then nu is taken from U1/U2
            else
               ! first check if 0/0 case
               if (dummy_u2_up .eq. 0.0_DP .and. dummy_u1_up .ne. 0.0_DP) then
                  call utils_abort('Error in internal_dft_nu(): &
                       &Hubbard U2_up = 0 in dft_nu block while &
                       &attempting nu=U1/U2!')
               elseif (dummy_u2_up .eq. 0.0_DP .and. &
                    & dummy_u1_up .eq. 0.0_DP) then
                  dummy_n_up = floor((-1.0_DP)*dummy_n_up)
                  dummy_nu_up  = 0.0_DP
               else
                  dummy_n_up = floor((-1.0_DP)*dummy_n_up)
                  dummy_nu_up  = 0.5*(dummy_u1_up/dummy_u2_up)
               endif
            endif

            ! Check N_dn - if positive then nu takes fixed value
            if (dummy_n_dn > 0.0_DP) then
               dummy_nu_dn = dummy_n_dn - floor(dummy_n_dn)
               dummy_n_dn = floor(dummy_n_dn)
               ! if N_dn negative then nu is taken from U1/U2
            else
               ! first check if 0/0 case
               if (dummy_u2_dn .eq. 0.0_DP .and. dummy_u1_dn .ne. 0.0_DP) then
                  call utils_abort('Error in internal_dft_nu(): &
                       &Hubbard U2_dn = 0 in dft_nu block while &
                       &attempting nu=U1/U2!')
               elseif (dummy_u2_dn .eq. 0.0_DP .and. &
                    & dummy_u1_dn .eq. 0.0_DP) then
                  dummy_n_dn = floor((-1.0_DP)*dummy_n_dn)
                  dummy_nu_dn  = 0.0_DP
               else
                  dummy_n_dn = floor((-1.0_DP)*dummy_n_dn)
                  dummy_nu_dn  = 0.5*(dummy_u1_dn/dummy_u2_dn)
               endif
            endif

            ATOMS_CELL: do ii=1,nat ! Loop over all atoms in cell
               if (elements(ii)%species_id == dummy_id) then

                  ! Store namelist of Hubbard atoms
                  h_species(row)%hub_species = dummy_id(1:4)

                  ! Angular momentum channel of projector
                  h_species(row)%cdft_ang_mom        = dummy_l
                  ! Copy cdft angular momentum into DFT-U counterpart
                  h_species(row)%hub_ang_mom         = dummy_l

                  ! Effective charge for radial function
                  h_species(row)%hub_charge          = ABS(dummy_c)

                  ! Hubbard U1 up
                  h_species(row)%nu_u1_up       = dummy_u1_up / HARTREE_IN_EVS

                  ! Hubbard U1 down
                  h_species(row)%nu_u1_down     = dummy_u1_dn / HARTREE_IN_EVS

                  ! Hubbard U2 up
                  h_species(row)%nu_u2_up       = dummy_u2_up / HARTREE_IN_EVS

                  ! Hubbard U2 down
                  h_species(row)%nu_u2_down     = dummy_u2_dn / HARTREE_IN_EVS

                  ! Target N up
                  h_species(row)%nu_target_n_up      = dummy_n_up

                  ! Target nu up
                  h_species(row)%nu_nu_up     = dummy_nu_up

                  ! Target N down
                  h_species(row)%nu_target_n_down    = dummy_n_dn

                  ! Target nu down
                  h_species(row)%nu_nu_down   = dummy_nu_dn

                  ! Species Group
                  h_species(row)%nu_group       = dummy_group

                  ! pub_hubbard_atomsolve = .true.
                  ! ddor: Trigger fireball NGWFs by using negative Z parameter
                  if (( dummy_c < 0.0_DP ).AND.&
                       & (.not.pub_cdft_read_projectors)) then
                     pub_hubbard_atomsolve = .true.
                  endif

                  ! gom: for cdft_multi_proj set number cDFT-proj=number val_NGWFs
                  if (pub_cdft_multi_proj) &
                       h_species(row)%cdft_num_proj  = elements(ii)%nfunctions

                  if ( ( h_species(row)%hub_ang_mom < 0 ) .and. &
                       &( h_species(row)%hub_ang_mom > 3 ) ) then
                     call utils_abort('Error in internal_dft_nu(): &
                          &unphysical value of angular momentum for &
                          &projector in <constrained_dft> block.')
                  endif


                  ! gom: count total number of Hubbard atoms
                  dummy_h = dummy_h + 1

                  ! gom: change projector counting scheme for cdft_multi_proj
                  if (pub_cdft_multi_proj) then
                     num_hub_proj = num_hub_proj + h_species(row)%cdft_num_proj
                  else
                     num_hub_proj = num_hub_proj + 2 * h_species(row)%hub_ang_mom + 1
                  endif

                  element_found = .true.

               end if

            enddo ATOMS_CELL

            ! gom: code taken from gibo's cdft module
            ! gibo: spare the user the need to set (i) pub_task=PROPERTIES,
            ! (ii) hubbard_proj_mixing <0, (iii) Z (constrained_dft block) >0
            ! for properties run
            if ( (pub_task=='PROPERTIES') .OR. (pub_task=='COND') .OR.       &
                 (pub_task=='PROPERTIES_COND') .OR. (pub_task=='TDDFT')) then
               pub_hubbard_atomsolve = .FALSE.
               pub_hubbard_restart   = .TRUE.
               pub_hub_proj_mixing   = -1.0_DP
               if (pub_on_root) write(stdout,'(a)') 'WARNING in internal_dft_nu&
                    &: pub_hubbard_restart overriden to T'
               if (pub_on_root) write(stdout,'(a)') 'WARNING in internal_dft_nu&
                    &: pub_hub_atomsolve overriden to F'
               if (pub_on_root) write(stdout,'(a)') 'WARNING in internal_dft_nu&
                    &: pub_hub_proj_mixing overriden to "-1.0_DP"'
            endif

            ! 4 GIBO: TO BE CHECKED AFTER PAW GO AHEAD FROM NICK ==== START
            !: cDFT within PAW is something we fancy...
            if (pub_hubbard_restart .and. pub_hubbard_atomsolve) then
               call utils_abort('Error in internal_dft_nu(): &
                    &Cannot have negative values for DFT+U charge Z &
                    &with PAW or negative projector mixing parameter &
                    &in <hubbard> block.')
            endif
            ! 4 GIBO: TO BE CHECKED AFTER PAW GO AHEAD FROM NICK ==== END

            ! gibo: cDFT-PROPERTIES-friendly
            call comms_bcast(pub_root_proc_id, pub_hubbard_restart)
            call comms_bcast(pub_root_proc_id, pub_hub_proj_mixing)
            ! gibo: cDFT-PROPERTIES-friendly
            call comms_bcast(pub_root_proc_id, pub_hubbard_atomsolve)
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_species)
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_ang_mom)
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_u)
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_charge)

            ! gom : Broadcast DFT+nu parameters
            call comms_bcast(pub_root_proc_id, h_species(row)%cdft_num_proj)
            call comms_bcast(pub_root_proc_id, h_species(row)%nu_u1_up)
            call comms_bcast(pub_root_proc_id, h_species(row)%nu_u1_down)
            call comms_bcast(pub_root_proc_id, h_species(row)%nu_u2_up)
            call comms_bcast(pub_root_proc_id, h_species(row)%nu_u2_down)
            call comms_bcast(pub_root_proc_id, h_species(row)%nu_target_n_up)
            call comms_bcast(pub_root_proc_id, h_species(row)%nu_nu_up)
            call comms_bcast(pub_root_proc_id, h_species(row)%nu_target_n_down)
            call comms_bcast(pub_root_proc_id, h_species(row)%nu_nu_down)
            call comms_bcast(pub_root_proc_id, h_species(row)%nu_group)

            !gibo: to avoid problems, initialise (unused) Hubbard variables
            !      and send them around
            h_species(row)%hub_alpha = 0._DP
            h_species(row)%hub_spin_splitting = 0._DP
            call comms_bcast(pub_root_proc_id, h_species(row)%hub_alpha)
            call comms_bcast(pub_root_proc_id, &
                 & h_species(row)%hub_spin_splitting)

         end do


         if (.not. element_found) then
            call utils_abort('Error in internal_dft_nu(): &
                 &mismatching species in <dft_nu> and <positions_abs> &
                 &blocks of input file.')
         else
            par%nat_hub = dummy_h ! The number of Hubbard atoms in the cell
            par%num_hub_proj = num_hub_proj
         end if


         ! ddor: broadcast number of Hubbard atoms in the cell
         call comms_bcast(pub_root_proc_id,par%nat_hub)
         call comms_bcast(pub_root_proc_id,par%num_hub_proj)


         if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Read <dft_nu> block'
      endif CHECK_CDFT

      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_dft_nu','cbuf',ierr)
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_dft_nu','dbuf',ierr)
      deallocate(ibuf,stat=ierr)
      call utils_dealloc_check('internal_dft_nu','ibuf',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving <internal_dft_nu>'

    end subroutine internal_dft_nu

    !-------------------------------
    !--------- %block classical_info
    !-------------------------------

    subroutine internal_classical_atoms()

      use classical_pot, only: classical_pot_init
      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_reduce, esdf_convfac
      use utils, only: utils_abort

      implicit none

      !-------------
      ! ars: buffers
      integer, allocatable :: ibuf_class(:)
      real(kind=DP), allocatable :: dbuf_class(:)
      ! ars: integer
      integer :: row, ii
      ! ars: character
      character(len=4),allocatable  :: species_id_classical(:)
      ! ars: real
      real(kind=DP), allocatable    :: classical_charge(:)
      ! ars: logical
      logical :: class_nats_exist

      ! rab207: units on first line
      character(len=80) :: dummy_string
      real(kind=DP) :: lenconfac
      integer :: block_offset

      ! cks: initialise number of classical atoms
      class_nat =0
      lenconfac =1.0_DP
      if (pub_on_root) then
         class_nats_exist =esdf_block('classical_info',class_nat)
         if(class_nats_exist) then ! jd: needed, or else class_nat can become -1
            ! rab207: check if 'bohr' or 'ang' units on first line
            read(block_data(1),*) dummy_string
            dummy_string = esdf_reduce(dummy_string)
            if ((index(dummy_string,"ang").ne.0).or. &
                 & (index(dummy_string,"bohr").ne.0)) then
               class_nat = class_nat - 1
               lenconfac = esdf_convfac(dummy_string,'bohr')
               block_offset = 1
            else
               lenconfac = 1.0_DP
               block_offset = 0
            endif
        endif
      endif

      ! kaw: Check if the arrays have already been allocated
      if (allocated(classical_elements)) class_nats_exist = .false.

      ! cks: make sure all cores have info about classical nats
      call comms_bcast(pub_root_proc_id,class_nat)
      call comms_bcast(pub_root_proc_id,class_nats_exist)
      call comms_bcast(pub_root_proc_id,lenconfac)

      ! Allocate array if not already allocated, even if size is zero
      if (.not.allocated(classical_elements)) then
         call classical_pot_init(classical_elements,class_nat)
      else if (size(classical_elements)/=class_nat) then
         call utils_abort('Error in internal_classical_atoms &
              &(rundat_blocks_read) size of classical_elements array does &
              &not match number of classical atoms')
      end if

      ! cks: allocate classical arrays only if there are classical atoms
      if (class_nats_exist) then

         allocate(species_id_classical(class_nat),stat=ierr)
         call utils_alloc_check('internal_classical_atoms',&
              'species_id_classical',ierr)
         allocate(classical_charge(class_nat),stat=ierr)
         call utils_alloc_check('internal_classical_atoms','classical_charge',&
              ierr)

         allocate(dbuf_class(3*class_nat),stat=ierr)
         call utils_alloc_check('internal_classical_atoms','dbuf_class',&
              ierr)

         allocate(ibuf_class(4*class_nat),stat=ierr)
         call utils_alloc_check('internal_classical_atoms','ibuf_class',&
              ierr)

         ! ars: Root proc only: read absolute cartesian positions from
         !      <classical_info> block of input file
         if (pub_on_root) then

            do row=1,class_nat
               read(block_data(row+block_offset),*) species_id_classical(row),&
                    dbuf_class(row*3-2:row*3),classical_charge(row)
               do ii=1,4
                  ibuf_class((row-1)*4+ii) = iachar(species_id_classical(row)(ii:ii))
               end do

               ! rab207: convert units
               dbuf_class(row*3-2:row*3) = lenconfac*dbuf_class(row*3-2:row*3)

            end do
         end if


         ! cks: broadcast positions and charges
         call comms_bcast(pub_root_proc_id,dbuf_class,class_nat*3)
         call comms_bcast(pub_root_proc_id,classical_charge,class_nat)

         ! cks: broadcast and unpack species_id_classical
         call comms_bcast(pub_root_proc_id,ibuf_class)
         do row=1,class_nat
            do ii=1,4
               species_id_classical(row)(ii:ii) = achar(ibuf_class((row-1)*4+ii))
            end do
         end do

         ! cks: unpack and store ids, coordinates and charges
         do row=1,class_nat
            classical_elements(row)%symbol = species_id_classical(row)(1:2)
            classical_elements(row)%centre%x = dbuf_class(row*3-2)
            classical_elements(row)%centre%y = dbuf_class(row*3-1)
            classical_elements(row)%centre%z = dbuf_class(row*3)
            classical_elements(row)%ion_charge = classical_charge(row)
         end do

         ! cks: deallocate classical arrays
         deallocate(ibuf_class,stat=ierr)
         call utils_dealloc_check('internal_classical_atoms',&
              'ibuf_class',ierr)

         deallocate(dbuf_class, stat=ierr)
         call utils_dealloc_check('internal_classical_atoms',&
              'dbuf_class',ierr)

         deallocate(classical_charge,stat=ierr)
         call utils_dealloc_check('internal_classical_atoms',&
              'classical_charge',ierr)

         deallocate(species_id_classical,stat=ierr)
         call utils_dealloc_check('internal_classical_atoms',&
              'species_id_classical',ierr)
      endif

    end subroutine internal_classical_atoms


    !--------------------------------------------------------------------------!
    ! lpl: NBO Partial matrix atom list & NGWF labelling (16/06/2011)          !
    !      Modified from %block_species_ngwf_plot in this module               !
    !      Reads: %block nbo_write_species                                     !
    !             %block nbo_species_ngwflabel                                 !
    !             %block nbo_list_plotnbo                                      !
    !--------------------------------------------------------------------------!

    subroutine internal_species_write_nbo()

      use comms,  only: pub_root_proc_id, comms_bcast
      use esdf,   only: block_data, esdf_block
      use rundat, only: pub_nbo_write_species, pub_nbo_ngwf_label, &
           pub_nbo_list_plotnbo, pub_plot_nbo
      use utils,  only: utils_abort

      implicit none

      ! lpl: Local variables
      character(len=4) :: dummy_id
      character(len=252) :: dummy_ngwf_label
      integer :: dummy_nsp_cbuf1   ! lpl: For <nbo_write_species>
      integer :: dummy_nsp_cbuf2   ! lpl: For <nbo_species_ngwflabel>
      integer :: nnbo, nbo_it

      ! lpl: Local list buffers
      ! lpl: For <nbo_write_species>
      character(len=68), allocatable :: cbuf1(:)
      ! lpl: For <nbo_species_ngwflabel>
      character(len=256),  allocatable :: cbuf2(:)

      ! lpl: Misc. scalars
      integer ::  row, ii, ierr

      ! lpl: Checking flags
      logical :: element_found
      logical :: skip_species_list   ! lpl: For <nbo_write_species>
      logical :: nbolist_exists, skip

      allocate(cbuf1(nsp),stat=ierr)
      call utils_alloc_check('internal_species_write_nbo','cbuf1',ierr)
      allocate(cbuf2(nsp),stat=ierr)
      call utils_alloc_check('internal_species_write_nbo','cbuf2',ierr)

      ! lpl: Allocate necessary arrays
      if (allocated(pub_nbo_write_species)) then
         deallocate(pub_nbo_write_species,stat=ierr)
         call utils_dealloc_check('internal_species_write_nbo', &
              'pub_nbo_write_species',ierr)
      end if
      allocate(pub_nbo_write_species(nat),stat=ierr)
      call utils_alloc_check('internal_species_write_nbo', &
           'pub_nbo_write_species',ierr)

      if (allocated(pub_nbo_ngwf_label)) then
         deallocate(pub_nbo_ngwf_label,stat=ierr)
         call utils_dealloc_check('internal_species_write_nbo', &
              'pub_nbo_ngwf_label',ierr)
      end if
      allocate(pub_nbo_ngwf_label(nsp),stat=ierr)
      call utils_alloc_check('internal_species_write_nbo', &
           'pub_nbo_ngwf_label',ierr)

      if (allocated(pub_nbo_list_plotnbo)) then
         deallocate(pub_nbo_list_plotnbo,stat=ierr)
         call utils_dealloc_check('internal_species_write_nbo', &
              'pub_nbo_list_plotnbo',ierr)
      end if

      ! lpl: Get NBO list block if plot_nbo = T
      nnbo = 0
      if(pub_plot_nbo) then

         if(pub_on_root) nbolist_exists &
              = esdf_block('nbo_list_plotnbo',nnbo)

         call comms_bcast(pub_root_proc_id,nbolist_exists)
         call comms_bcast(pub_root_proc_id,nnbo)

         if(nbolist_exists .and. nnbo > 0) then
            allocate(pub_nbo_list_plotnbo(0:nnbo),stat=ierr)
            call utils_alloc_check('internal_species_write_nbo', &
                 'pub_nbo_list_plotnbo',ierr)

            pub_nbo_list_plotnbo(0) = 0 ! idx 0 stores #NBOs
            pub_nbo_list_plotnbo(1:nnbo) = -1 ! lpl: -1 for error checking
         else
            call utils_abort('Error in internal_species_write_nbo() reading &
                 &block nbo_list_plotnbo.')
         end if

         if(pub_on_root) then ! lpl: Read block

            do row=1,nnbo
               read(block_data(row),*) nbo_it
               if(nbo_it < 0) then
                  call utils_abort('Error in internal_species_write_nbo(): &
                       &nbo_it < 0')
               end if

               ! lpl: Check for multiply defined NBOs
               if(row == 1) then ! Read 1st element
                  pub_nbo_list_plotnbo(0) = 1
                  pub_nbo_list_plotnbo(1) = nbo_it
               else ! Check against read elements

                  skip = .false.
                  do ii=1,pub_nbo_list_plotnbo(0)-1
                     if(nbo_it == pub_nbo_list_plotnbo(ii)) then
                        write(stdout,'(a,I5,a)') 'WARNING: NBO # ', &
                             nbo_it,' multiply defined.'
                        skip = .true.
                        exit
                     end if
                  end do ! END  do col=1,row-1
                  if(.not. skip) then
                     pub_nbo_list_plotnbo(0) = & ! idx 0 stores # NBOs
                          pub_nbo_list_plotnbo(0) + 1
                     pub_nbo_list_plotnbo(pub_nbo_list_plotnbo(0)) &
                          = nbo_it
                  end if

               end if
            end do ! END do row=1,nnbo

         end if ! END if(pub_on_root)

      end if ! END if(pub_plot_nbo)

      ! lpl: Broadcast pub_nbo_list_plotnbo if present & allocated
      if (allocated(pub_nbo_list_plotnbo)) &
           call comms_bcast(pub_root_proc_id,pub_nbo_list_plotnbo)

      ! lpl: initialise info arrays to default
      pub_nbo_write_species = .false.
      do row=1,nsp
         pub_nbo_ngwf_label(row)(1:4) = all_species_id(row)(1:4)
         pub_nbo_ngwf_label(row)(5:)  = 'AUTO'
      end do

      ! lpl: Read info from <nbo_write_species> and <nbo_species_ngwflabel>
      !      blocks in input file on root proc only. If blocks are
      !      not present default values are assumed
      if (pub_on_root) then

         ! lpl: <nbo_write_species>
         skip_species_list = .false.
         if (esdf_block('nbo_write_species',dummy_nsp_cbuf1)) then
            if (dummy_nsp_cbuf1 > nsp) then ! aam: error check
               call utils_abort('Error in internal_species_write_nbo(): too&
                    & many species in <nbo_write_species> block of input file.')
            else
               ! aam: check for multiply defined species
               call internal_check_species(dummy_nsp_cbuf1, &
                    '<nbo_write_species>')
               ! aam: dummy_nsp is not necessarily equal to nsp
               do row=1,dummy_nsp_cbuf1
                  read(block_data(row),*) dummy_id
                  cbuf1(row)(1:4) = dummy_id
               end do
            end if
         else
            ! lpl: If block isn't specified, skip reading species list
            skip_species_list = .true.
         end if

         ! lpl: <nbo_species_ngwflabel>
         if (esdf_block('nbo_species_ngwflabel',dummy_nsp_cbuf2)) then
            if (dummy_nsp_cbuf2 > nsp) then ! aam: error check
               call utils_abort('Error in internal_species_write_nbo: too many &
                    &species in <nbo_species_ngwflabel> block of input file.')
            else
               ! aam: check for multiply defined species
               call internal_check_species(dummy_nsp_cbuf2, &
                    '<nbo_species_ngwflabel>')
               ! aam: dummy_nsp is not necessarily equal to nsp
               do row=1,dummy_nsp_cbuf2
                  read(block_data(row),*) dummy_id,dummy_ngwf_label
                  cbuf2(row)(1:4) = dummy_id
                  cbuf2(row)(5:)  = dummy_ngwf_label
               end do
            end if
         end if

      end if

      ! lpl: broadcast string buffers & indicators
      call comms_bcast(pub_root_proc_id,dummy_nsp_cbuf1)
      call comms_bcast(pub_root_proc_id,dummy_nsp_cbuf2)
      do row=1,dummy_nsp_cbuf1
         call comms_bcast(pub_root_proc_id,cbuf1(row),4)
      end do
      do row=1,dummy_nsp_cbuf2
         call comms_bcast(pub_root_proc_id,cbuf2(row))
      end do
      call comms_bcast(pub_root_proc_id,skip_species_list)

      ! lpl: extract information from <nbo_write_species>
      if(.not. skip_species_list) then
         do row=1,dummy_nsp_cbuf1

            element_found = .false. ! aam: error check flag
            dummy_id = cbuf1(row)(1:4)

            ! lpl: <nbo_write_species> loop
            do ii=1,nat
               if (elements(ii)%species_id == dummy_id) then
                  pub_nbo_write_species(ii) = .true.
                  element_found = .true.
               end if
            end do

            if (.not. element_found) then
               call utils_abort('Error in internal_nbo_write_species(): &
                    &mismatching species in <internal_species_write_nbo> and &
                    &<positions_abs> blocks of input file.')
            end if
         end do
      else
         ! lpl: If no nbo_write_species block is specified, include all atoms
         !      in output
         pub_nbo_write_species = .true.
      end if

      ! lpl: extract information from <nbo_species_ngwflabel>
      do row=1,dummy_nsp_cbuf2

         element_found = .false. ! aam: error check flag
         dummy_id = cbuf2(row)(1:4)
         dummy_ngwf_label = cbuf2(row)(5:)

         ! lpl: Check that this species exists
         do ii=1,nat
            if (elements(ii)%species_id == dummy_id) element_found = .true.
         end do

         if (.not. element_found) then
            call utils_abort('Error in internal_species_write_nbo: mismatching &
                 &species in <nbo_species_ngwflabel> and <positions_abs> blocks&
                 & of input file.')
         end if

         ! lpl: Overwrite default ('AUTO') with use-specified configuration
         do ii=1,nsp
            if(pub_nbo_ngwf_label(ii)(1:4) == dummy_id) then
               pub_nbo_ngwf_label(ii)(5:)  = dummy_ngwf_label
            end if
         end do

      end do

      ! lpl: Deallocate buffers
      deallocate(cbuf1,stat=ierr)
      call utils_dealloc_check('internal_species_write_nbo','cbuf1',ierr)
      deallocate(cbuf2,stat=ierr)
      call utils_dealloc_check('internal_species_write_nbo','cbuf2',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <nbo_write_species> block'

    end subroutine internal_species_write_nbo

    !--------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !--------------------------------------------------------------------------!

    subroutine internal_ddec_read_rcomp()

      !-----------------------------------------------------------------------!
      ! lpl: Reads %block ddec_rcomp                                          !
      !      Only checks for unidentified species, not unspecified ones       !
      !-----------------------------------------------------------------------!

      use comms, only: pub_root_proc_id, comms_barrier, comms_bcast
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_ddec_rcomp
      use utils, only: utils_abort

      implicit none

      character(len=512)  :: dummy_id
      integer :: num_lines
      integer :: row, ii, isp
      logical :: isp_found, block_exists

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Start reading <ddec_rcomp> block'

      block_exists = .false.
      ! lpl: Allocate necessary arrays
      if (allocated(pub_ddec_rcomp)) then
         deallocate(pub_ddec_rcomp,stat=ierr)
         call utils_dealloc_check('internal_ddec_read_rcomp', &
              'pub_ddec_rcomp',ierr)
      end if

      ! Root proc only: read information
      if (pub_on_root) then
         block_exists = esdf_block('ddec_rcomp',num_lines)
      end if
      call comms_bcast(pub_root_proc_id,block_exists)
      if ( .not. block_exists) return

      call comms_bcast(pub_root_proc_id,num_lines)

      ! Allocate num_lines of data
      allocate(pub_ddec_rcomp(num_lines),stat=ierr)
      call utils_alloc_check('internal_ddec_read_rcomp', &
           'pub_ddec_rcomp',ierr)
      pub_ddec_rcomp = ''

      call comms_barrier

      ! Resume reading data on root proc and bcast to all
      if (pub_on_root) then
         do row=1,num_lines
            read(block_data(row),'(a)') pub_ddec_rcomp(row)
         end do
      end if
      do row=1,num_lines
         call comms_bcast(pub_root_proc_id,pub_ddec_rcomp(row))
      end do

      do row=1,num_lines
         isp_found = .false. ! error check flag
         read(pub_ddec_rcomp(row),*) dummy_id

         ! lpl: Find real isp for this dummy_id
         do ii=1,nat
            if ( trim(adjustl(elements(ii)%species_id)) &
                 == trim(adjustl(dummy_id)) ) then
               isp = elements(ii)%species_number
               isp_found = .true.
               exit
            end if
         end do

         ! lpl: If species specified in block does not match any
         if (.not. isp_found) call utils_abort('Error in &
              &internal_ddec_read_rcomp: unidentified &
              &species in <ddec_rcomp> block of input file')

         if(pub_debug_on_root) write(stdout,'(a,1x,a,1x,a)') &
              'DEBUG: <ddec_rcomp>: ', &
              elements(isp)%species_id, pub_ddec_rcomp(row)

      end do

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Finished reading <ddec_rcomp> block'

    end subroutine internal_ddec_read_rcomp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if 0
    subroutine internal_ddec_rmse_vdW()

      !-----------------------------------------------------------------------!
      ! lpl: Reads %block ddec_rmse_vdW                                       !
      !      Used by the still-dodgy subroutine 'ddec_rmse'. Presence of this !
      !      block will initiate the said subroutine                          !
      !-----------------------------------------------------------------------!

      use comms, only: comms_barrier, comms_bcast, pub_root_proc_id
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_ddec_rmse_vdW
      use utils, only: utils_abort

      implicit none

      ! dummies
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      character(len=64) :: dummy_ddec_vdW
      ! buffers
      character(len=68), allocatable :: cbuf(:)
      ! integers
      integer :: row, ii, isp
      ! logical
      logical :: isp_found

      ! lpl: Allocate necessary arrays
      if (allocated(pub_ddec_rmse_vdW)) then
         deallocate(pub_ddec_rmse_vdW,stat=ierr)
         call utils_dealloc_check('internal_ddec_rmse_vdW', &
              'pub_ddec_rmse_vdW',ierr)
      end if
      allocate(pub_ddec_rmse_vdW(nsp),stat=ierr)
      call utils_alloc_check('internal_ddec_rmse_vdW', &
           'pub_ddec_rmse_vdW',ierr)

      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_ddec_rmse_vdW','cbuf',ierr)

      ! Root proc only: read information from
      ! <species_atomic_set> block of input file
      if (pub_on_root) then
         ! lpl: Verify block for all species
         if (esdf_block('ddec_rmse_vdW',dummy_nsp)) then
            if (nsp == dummy_nsp)  then
               ! check for multiply defined species
               call internal_check_species(nsp,'<ddec_rmse_vdW>')
               do row=1,nsp
                  read(block_data(row),*) dummy_id, dummy_ddec_vdW
                  cbuf(row)(1:4) = dummy_id
                  cbuf(row)(5:68) = dummy_ddec_vdW
               end do
            else
               call utils_abort('Error in internal_ddec_rmse_vdW: &
                    &mismatching numbers of species in <ddec_rmse_vdW> block &
                    &of input file')
            end if

            ! lpl: Transfer verified block into pub_ddec_rmse_vdW
            do row=1,nsp
               isp_found = .false. ! error check flag
               dummy_id = cbuf(row)(1:4)
               dummy_ddec_vdW = cbuf(row)(5:68)

               ! lpl: Find real isp for this dummy_id
               do ii=1,nat
                  if (elements(ii)%species_id == dummy_id) then
                     isp = elements(ii)%species_number
                     exit
                  end if
               end do
               pub_ddec_rmse_vdW(isp) = trim(adjustl(dummy_id))//' '//&
                    trim(adjustl(dummy_ddec_vdW))
               isp_found = .true.
               if (.not. isp_found) call utils_abort('Error in &
                    &internal_ddec_rmse_vdW: mismatching species in &
                    &<ddec_rmse_vdW> and <positions_abs> blocks of input file')
            end do
         end if
      end if
      call comms_bcast(pub_root_proc_id,pub_ddec_rmse_vdW)

      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_ddec_rmse_vdW','cbuf',ierr)

      call comms_barrier

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <ddec_rmse_vdW> block'

    end subroutine internal_ddec_rmse_vdW
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

    !----------------------------
    !--------- %block couplings
    !----------------------------

    subroutine internal_couplings()

      use comms, only: pub_root_proc_id, comms_bcast, &
           comms_barrier
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_lr_tddft_calculate, pub_task, &
           pub_cond_calculate, pub_cond_calculate_any_task
      use utils, only: utils_abort
      use couplings, only: couplings_state_ids, couplings_state_names, &
           couplings_state_energies, couplings_state_types, couplings_nstates, &
           couplings_state_exc_idxs

      implicit none


      ! dhpt: dummys
      character(len=4)  :: dummy_id
      character(len=64) :: dummy_name
      real(kind=DP)     :: dummy_energy
      character(len=5)  :: dummy_type
      integer           :: dummy_exc_idx
      ! ars: integers
      integer :: row
      integer :: nstates
      logical :: states_exist

      ! dhpt: Root proc only: read information from
      !       <couplings> block of input file
      if (pub_on_root) then
         states_exist = esdf_block('couplings_states',nstates)
         if (.not.states_exist) nstates = 0
      end if

      call comms_bcast(pub_root_proc_id,nstates)
      call comms_barrier

      ! allocate couplings arrays
      if (.not.allocated(couplings_state_ids)) then
         allocate(couplings_state_ids(nstates),stat=ierr)
         call utils_alloc_check('internal_couplings','couplings_state_ids',ierr)
         allocate(couplings_state_names(nstates),stat=ierr)
         call utils_alloc_check('internal_couplings','couplings_state_names',ierr)
         allocate(couplings_state_energies(nstates),stat=ierr)
         call utils_alloc_check('internal_couplings','couplings_state_energies',ierr)
         allocate(couplings_state_types(nstates),stat=ierr)
         call utils_alloc_check('internal_couplings','couplings_state_types',ierr)
         allocate(couplings_state_exc_idxs(nstates),stat=ierr)
         call utils_alloc_check('internal_couplings','couplings_state_exc_idxs',ierr)
      end if

      couplings_nstates = nstates

      if (pub_on_root) then
         do row=1,nstates
            read(block_data(row),*) dummy_id, dummy_type, &
                                   &dummy_name, dummy_exc_idx, &
                                   &dummy_energy
            couplings_state_ids(row)      = dummy_id
            couplings_state_names(row)    = dummy_name
            couplings_state_types(row)    = dummy_type
            couplings_state_energies(row) = dummy_energy
            couplings_state_exc_idxs(row) = dummy_exc_idx
            if (trim(dummy_type) == 'tddft') then
               if (pub_task=='COUPLINGS') then
                  if (pub_debug_on_root) write(stdout,'(a)') 'TDDFT type detected; &
                      &setting pub_lr_tddft_calculate and pub_cond_calculate &
                      &to true.'
                  pub_lr_tddft_calculate = .true.
                  pub_cond_calculate     = .true.
               end if
               pub_cond_calculate_any_task = .true.
            else if (trim(dummy_type) /= 'cdft') then
               call utils_abort('Error in internal_couplings(): &
                    &Unknown type '//trim(dummy_type))
            end if
         end do
      end if

      ! broadcast
      call comms_bcast(pub_root_proc_id,pub_lr_tddft_calculate)
      call comms_bcast(pub_root_proc_id,pub_cond_calculate)
      do row=1,nstates
         call comms_bcast(pub_root_proc_id,couplings_state_ids(row))
         call comms_bcast(pub_root_proc_id,couplings_state_names(row))
         call comms_bcast(pub_root_proc_id,couplings_state_types(row))
         call comms_bcast(pub_root_proc_id,couplings_state_energies(row))
         call comms_bcast(pub_root_proc_id,couplings_state_exc_idxs(row))
      end do
      call comms_barrier

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <couplings> block'

    end subroutine internal_couplings


    !-----------------------------
    !--------- %block species_ngwf_regions
    !-----------------------------
    ! rc2013: New block to allocate atoms to each region.
    subroutine internal_species_ngwf_regions()

      use comms, only: pub_root_proc_id, comms_bcast
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_ngwf_regions_ngroups, pub_ngwf_regions_group_nsp,&
           pub_ngwf_regions_groups
      use utils, only: utils_abort

      implicit none

      !-------------
      ! ars: dummies
      character(len=80) :: dummy_ngwf_regions_group
      character(len=4)  :: dummy_id
      character(len=8)       :: row_string
      integer :: dummy_nsp
      !-------------
      ! ars: integers
      integer :: row, ii, jj, iat, isub
      ! ars: logical
      logical :: element_found
      integer :: true_ngroups

      character(len=*), parameter :: myself = 'internal_species_ngwf_regions'


      ! ndmh: Root proc only: read information from
      !      <species_ngwf_regions> block of input file
      ! ndmh: check if the block exists
      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering species_ngwf_regions'
      if (pub_on_root) then
         if (esdf_block('species_ngwf_regions',pub_ngwf_regions_ngroups)) then
            if (pub_ngwf_regions_ngroups > nat) then
               call utils_abort('Error in internal_species_ngwf_regions(): &
                    &more groups than atoms defined in <species_ngwf_regions>.')
            end if
         else
            pub_ngwf_regions_ngroups = 0
            dummy_nsp = 0
         end if
      end if
      ! ndmh: allocate array for numbers of species in each group
      call comms_bcast(pub_root_proc_id,pub_ngwf_regions_ngroups)
      if (allocated(pub_ngwf_regions_group_nsp)) then
         deallocate(pub_ngwf_regions_group_nsp,stat=ierr)
         call utils_dealloc_check('internal_species_ngwf_regions', &
              'pub_ngwf_regions_group_nsp',ierr)
      end if
      allocate(pub_ngwf_regions_group_nsp(pub_ngwf_regions_ngroups),stat=ierr)
      call utils_alloc_check('internal_species_ngwf_regions', &
           'pub_ngwf_regions_group_nsp',ierr)

      ! ndmh: count number of species in each group
      if (pub_on_root.and.(pub_ngwf_regions_ngroups > 0)) then
         dummy_nsp = 0
         ! loop over lines in block (LDOS groups)
         do row=1,pub_ngwf_regions_ngroups
            dummy_ngwf_regions_group = trim(adjustl(block_data(row)))
            pub_ngwf_regions_group_nsp(row) = 1
            ! ndmh: count number of spaces (distinct species) in this group
            jj = 0
            do ii=1,len_trim(dummy_ngwf_regions_group)
               if ((dummy_ngwf_regions_group(ii:ii)==' ').and. &
                    (dummy_ngwf_regions_group(ii+1:ii+1)/=' ')) then
                  if (ii > jj + 1) then
                     pub_ngwf_regions_group_nsp(row) = pub_ngwf_regions_group_nsp(row) + 1
                     jj = ii
                  else
                     jj = ii
                  end if
               end if
            end do
         end do
      end if

      ! ndmh: allocate array for storing LDOS group species ID's
      call comms_bcast(pub_root_proc_id,pub_ngwf_regions_group_nsp)
      dummy_nsp = maxval(pub_ngwf_regions_group_nsp)
      call comms_bcast(pub_root_proc_id,dummy_nsp)
      if (allocated(pub_ngwf_regions_groups)) then
         deallocate(pub_ngwf_regions_groups,stat=ierr)
         call utils_dealloc_check('internal_species_ngwf_regions','pub_ngwf_regions_groups',&
              ierr)
      end if
      allocate(pub_ngwf_regions_groups(dummy_nsp,pub_ngwf_regions_ngroups),stat=ierr)
      call utils_alloc_check('internal_species_ngwf_regions','pub_ngwf_regions_groups',ierr)
      pub_ngwf_regions_groups(:,:) = ''

      ! ndmh: now fill the pub_ldos_groups array with these strings
      if (pub_on_root.and.(pub_ngwf_regions_ngroups > 0)) then
         do row=1,pub_ngwf_regions_ngroups
            dummy_ngwf_regions_group = trim(adjustl(block_data(row)))
            ii = 1
            jj = len_trim(dummy_ngwf_regions_group)
            do dummy_nsp=1,pub_ngwf_regions_group_nsp(row)
               if (index(dummy_ngwf_regions_group(ii:jj),' ') > 0) then
                  jj = ii + index(dummy_ngwf_regions_group(ii:jj),' ') - 1
                  do
                     if (dummy_ngwf_regions_group(jj+1:jj+1)/=' ') exit
                     jj = jj + 1
                  end do
               end if
               pub_ngwf_regions_groups(dummy_nsp,row) = &
                    adjustl(trim(dummy_ngwf_regions_group(ii:jj)))
               ii = jj + 1
               jj = len_trim(dummy_ngwf_regions_group)
            end do
         end do
      end if

      ! ndmh: broadcast LDOS groups from root proc
      do row=1,pub_ngwf_regions_ngroups
         do dummy_nsp=1,pub_ngwf_regions_group_nsp(row)
            call comms_bcast(pub_root_proc_id,pub_ngwf_regions_groups(dummy_nsp,row))
         end do
      end do

      ! rc2013: if no block was provided let's treat the atoms
      !         like they're all in 1 region
      true_ngroups = pub_ngwf_regions_ngroups
      if(pub_ngwf_regions_ngroups == 0) pub_ngwf_regions_ngroups = 1

      ! rc2013: allocate regions and species structures here
      if (.not.allocated(nat_reg)) then
          allocate(nat_reg(pub_ngwf_regions_ngroups),stat=ierr)
         call utils_alloc_check('internal_species_ngwf_regions','nat_reg',ierr)
      end if
      if (.not.allocated(nsp_reg)) then
          allocate(nsp_reg(pub_ngwf_regions_ngroups),stat=ierr)
         call utils_alloc_check('internal_species_ngwf_regions','nsp_reg',ierr)
      end if
      !if (.not.allocated(num)) then
      !    allocate(num(pub_ngwf_regions_ngroups),stat=ierr)
      !   call utils_alloc_check('internal_species_ngwf_regions','num',ierr)
      !end if
      !if (.not.allocated(num_cond)) then
      !    allocate(num_cond(pub_ngwf_regions_ngroups),stat=ierr)
      !   call utils_alloc_check('internal_species_ngwf_regions','num_cond',ierr)
      !end if
      !if (.not.allocated(num_aux)) then
      !    allocate(num_aux(pub_ngwf_regions_ngroups),stat=ierr)
      !   call utils_alloc_check('internal_species_ngwf_regions','num_aux',ierr)
      !end if

      ! ndmh: check LDOS group species match elements array species id's
      ! rc2013: divide info into subsystems
      nat_reg = 0
      if(true_ngroups .gt. 0) nsp_reg = pub_ngwf_regions_group_nsp
      if (pub_ngwf_regions_ngroups > 1) then
         do row=1,pub_ngwf_regions_ngroups
            do ii=1,pub_ngwf_regions_group_nsp(row)
               element_found = .false.
               dummy_id = pub_ngwf_regions_groups(ii,row)
               ! rc2013: write the integer to a string
               write(row_string,'(i5)') row
               do iat=1,nat
                  if (elements(iat)%species_id == dummy_id) then
                     element_found = .true.
                     elements(iat)%region          = row
                     ! rc2013: count the number of atoms in this region
                     nat_reg(row) = nat_reg(row) + 1
                  end if
               end do
               if (.not. element_found) then
                  call utils_abort('Error in internal_species_ngwf_regions: &
                       &mismatching species in <species_ngwf_regions> and &
                       &<positions_abs> blocks of input file.')
               endif
            end do
         end do
      else
         ! rc2013: set the regions label even if there are none specified
         do iat=1,nat
            elements(iat)%region = 1
            nat_reg(1) = nat
            nsp_reg(1) = nsp
         end do
      end if

            ! rc2013: allocate regions lists
      if(.not.(allocated(regions))) then
         allocate(regions(pub_ngwf_regions_ngroups),stat=ierr)
         call utils_alloc_check(myself,'regions',ierr)
      end if

      do isub=1,pub_ngwf_regions_ngroups
         ! rc2013: Allocate the regions%elements lists
         if (.not.allocated(regions(isub)%elements)) then
            allocate(regions(isub)%elements(nat_reg(isub)),stat=ierr)
            call utils_alloc_check(myself,'elements',ierr)
         else if (size(regions(isub)%elements)/=nat_reg(isub)) then
            call utils_abort('Error in '//myself//' &
                 &size of elements array does not match number of atoms')
         end if
         ! Do the same with species
         if (.not.allocated(regions(isub)%species)) then
            allocate(regions(isub)%species(nsp_reg(isub)),stat=ierr)
            call utils_alloc_check(myself,'species',ierr)
         else if (size(regions(isub)%species)/=nsp_reg(isub)) then
            call utils_abort('Error in '//myself//' &
                 &(rundat_blocks_read) size of species array does &
                 &not match number of species')
         end if
      end do

      par=>regions(1)%par

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_ngwf_regions> block'

    end subroutine internal_species_ngwf_regions

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!---------------------------   END DATA BLOCKS    --------------------------!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!





    !-----------------------
    !--------- check_species
    !-----------------------
    subroutine internal_check_species(nrecords,blockname)

      use esdf, only: block_data
      use utils, only: utils_abort

      implicit none

      integer, intent(in) :: nrecords
      character(len=*), intent(in) :: blockname

      ! <<<local variables>>>
      integer :: pp,qq

      do pp=1,nrecords
         read(block_data(pp),*) id_list(pp)
      enddo
      do pp=1,nrecords
         do qq=pp+1,nrecords
            if (id_list(pp).eq.id_list(qq)) then
               call utils_abort('Error in internal_check_species(): &
                    &multiply defined species in '//trim(blockname)//' block &
                    &of input file.')
            endif
         enddo
      enddo

    end subroutine internal_check_species

    subroutine internal_vdwparam_override()

      !==================================================================!
      ! This subroutine reads the VDW parameter override block from the  !
      ! input file and dumps it into vdwparam_override in the            !
      ! vdwcorrection module.                                            !
      !------------------------------------------------------------------!
      ! Written by Quintin Hill in December 2010.                        !
      !==================================================================!

      use comms, only: pub_on_root
      use esdf, only: block_data, esdf_block
      use vdwcorrection, only: vdwcorrection_override_alloc

      implicit none

      integer :: num_rows ! number of rows

      if (esdf_block('vdw_params',num_rows)) then
         call vdwcorrection_override_alloc(block_data(1:num_rows),num_rows)
      end if

    end subroutine internal_vdwparam_override

    subroutine internal_subsystem_allocate()

      !==================================================================!
      ! This subroutine takes the mdl%elements and mdl%species arrays    !
      ! and allocates their data to the relevant subsystems.             !
      !------------------------------------------------------------------!
      ! Created by Robert Charlton, 12/12/2016.                          !
      !==================================================================!

      use utils, only: utils_abort, utils_assert, utils_flushed_string_output
      use rundat, only: pub_ngwf_regions_ngroups

      implicit none

      integer :: ierr
      integer :: ii, jj, kk, ll ! rc2013: loop iterators
      integer :: iat, loc_isp1, loc_isp2, iregion
      ! rc2013: arrays to count local parameters
      integer :: jat(pub_ngwf_regions_ngroups)
      integer :: jsp(pub_ngwf_regions_ngroups)
      integer :: loc_isp(pub_ngwf_regions_ngroups)
      logical :: new_species
      character(len=*), parameter :: myself = 'internal_subsystem_allocate'

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering internal_subsystem_allocate'


      ! rc2013: Allocate data from mdl%elements to regions%elements arrays.
      jj=1
      kk=1
      ! Initialise arrays
      jat = 0
      jsp = 0
      loc_isp = 1

      do ii=1,nat
         ! Get the region counter
         iregion = elements(ii)%region
         ! Count the no. of atoms allocated to each region
         jat(iregion) = jat(iregion) + 1
         ! Copy all elements info to the regions structure
         regions(iregion)%elements(jat(iregion)) = elements(ii)
         ! Mark each atom in regions with corresponding label in mdl%elements
         regions(iregion)%elements(jat(iregion))%global_atom_number = ii
         ! Count the NGWFs in each region
         !num(iregion) = num(iregion) + elements(ii)%nfunctions
         !num_aux(iregion) = num_aux(iregion) + elements(ii)%nfunctions_aux
         !num_cond(iregion) = num_cond(iregion) + elements(ii)%nfunctions_cond
      end do
      ! Repeat for species structure
      do ii=1,nsp
         ! Get the region counter
         iregion = species(ii)%region
         ! Count the no. of species allocated to each region
         jsp(iregion) = jsp(iregion) + 1
         ! Copy all species info to the regions structure
         regions(iregion)%species(jsp(iregion)) = species(ii)

         ! rc2013: set local labels for each elements list's species number
         new_species = .false.
         do iat=1,nat_reg(iregion)
            if (regions(iregion)%elements(iat)%global_species_number == ii) then
               ! rc2013: localise species references for this subsystem
               regions(iregion)%elements(iat)%species_number = loc_isp(iregion)
               new_species = .true.
            end if
         end do
         if (new_species) loc_isp(iregion) = loc_isp(iregion) + 1
      end do
      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving internal_subsystem_allocate'

    end subroutine internal_subsystem_allocate

  end subroutine rundat_blocks_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine rundat_blocks_exit(mdl)

    !=============================================================!
    ! This subroutine deallocates arrays read in this module      !
    !=============================================================!
    ! Written by Nicholas Hine in April 2016.                     !
    ! Modified for multiple subsystems by Robert Charlton,        !
    ! 19/12/2017.                                                 !
    !=============================================================!

    use dftb, only: dftb_exit
    use model_type, only: MODEL
    use rundat, only: pub_dftb
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(MODEL), intent(inout) :: mdl

    ! Local Variables
    integer :: ierr, isub

    if (pub_dftb) call dftb_exit(mdl)

    ! rc2013: loop over all subsystems
    do isub=1,mdl%nsub
       if(allocated(mdl%regions(isub)%species)) then
          deallocate(mdl%regions(isub)%species,stat=ierr)
          call utils_dealloc_check('internal_species_ngwf_regions','species', ierr)
       endif
       if(allocated(mdl%regions(isub)%classical_elements)) then
          deallocate(mdl%regions(isub)%classical_elements,stat=ierr)
          call utils_dealloc_check('internal_species_ngwf_regions','classical_elements', ierr)
       end if
       if(allocated(mdl%regions(isub)%elements)) then
          deallocate(mdl%regions(isub)%elements,stat=ierr)
          call utils_dealloc_check('internal_species_ngwf_regions','elements', ierr)
       endif
    enddo

    ! rc2013: now deallocate the mdl structures and regions too
    deallocate(mdl%species, stat=ierr)
    call utils_dealloc_check('internal_species','species', ierr)
    deallocate(mdl%classical_elements, stat=ierr)
    call utils_dealloc_check('classical_pot_init','classical_elements', ierr)
    deallocate(mdl%elements, stat=ierr)
    call utils_dealloc_check('internal_positions','elements', ierr)
    deallocate(mdl%regions, stat=ierr)
    call utils_dealloc_check('internal_species_ngwf_regions','regions', ierr)

  end subroutine rundat_blocks_exit

  subroutine thermostat_read_params(line,temp,tkind,tstart,tstop)

    use constants, only: dp,stdout
    use esdf,      only: esdf_line_divide, esdf_convfac, esdf_reduce
    use utils,     only: utils_abort

    implicit none

    ! Arguments
    character(len=80), intent(in) :: line
    real(kind=dp), intent(out)  :: temp
    integer, intent(out)        :: tkind, tstart, tstop

    ! Local variables
    character(len=80) :: cbuf(5), cjunk
    character(len=20) :: runit
    real(kind=dp)     :: confac
    integer           :: maxp

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Enter thermostat_read_params'

    temp = 0.0_dp
    tkind = 0
    tstart = 0
    tstop = 0

    maxp = 5
    call esdf_line_divide(maxp,cbuf,line)

    if (maxp .ne. 5) then
       call utils_abort('Error in thermostat_read_params: &
            &mandatory parameters not found in input file.')
    endif

    ! Read thermostat range
    read(cbuf(1),*) tstart
    read(cbuf(2),*) tstop

    ! Process thermostat type
    cjunk = esdf_reduce(cbuf(3))
    select case (cjunk)
    case ('none')
       tkind = 0
    case ('andersen')
       tkind = 1
    case ('langevin')
       tkind = 2
    case ('nosehoover')
       tkind = 3
    case ('berendsen')
       tkind = 4
    case ('bussi')
       tkind = 5
    end select

    ! Read thermostat temperature
    read(cbuf(4),*) temp
    read(cbuf(5),*) runit
    confac = esdf_convfac(runit,'hartree')
    temp = temp*confac


  end subroutine thermostat_read_params

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine thermostat_read_options(line,tgrad,tau,mix,damp,group,nhc,nhi,upd)

    use constants, only: dp
    use esdf, only: esdf_line_divide, esdf_convfac, esdf_reduce
    use rundat, only: md_delta_t

    implicit none

    ! Arguments
    character(len=80), intent(in) :: line
    real(kind=dp), intent(inout) :: tgrad, tau, mix, damp
    integer, intent(inout) :: group, nhc, nhi
    logical, intent(inout) :: upd

    ! Local variables
    character(len=80) :: cbuf(16), cjunk
    character(len=20) :: runit
    real(kind=dp)     :: confac
    integer           :: ip, maxp

    maxp = 16
    call esdf_line_divide(maxp,cbuf,line)

    ip = 1
    do while (ip.lt.maxp)
       cjunk = esdf_reduce(cbuf(ip))
       select case (cjunk)
       case ('tgrad')
          read(cbuf(ip+1),*) tgrad
          read(cbuf(ip+2),*) runit
          confac = esdf_convfac(runit,'hartree')
          tgrad = tgrad*confac/md_delta_t
          ip = ip + 3
       case ('tau')
          read(cbuf(ip+1),*) tau
          read(cbuf(ip+2),*) runit
          confac = esdf_convfac(runit,'aut')
          tau = tau*confac
          ip = ip + 3
       case ('mix')
          read(cbuf(ip+1),*) mix
          ip = ip + 2
       case ('damp')
          read(cbuf(ip+1),*) damp
          ip = ip + 2
       case ('group')
          read(cbuf(ip+1),*) group
          ip = ip + 2
       case ('nchain')
          read(cbuf(ip+1),*) nhc
          ip = ip + 2
       case ('nstep')
          read(cbuf(ip+1),*) nhi
          ip = ip + 2
       case ('update')
          read(cbuf(ip+1),*) upd
          ip = ip + 2
       case DEFAULT
          ip = ip + 1
       end select
    enddo

  end subroutine thermostat_read_options

  !------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine rundat_blocks_init_grid(a1, a2, a3, a1_pad, a2_pad, a3_pad, &
       mdl)

    use constants, only: DP, stdout
    use geometry, only: point
    use ion, only: element
    use model_type, only: MODEL
    use ppd_strategy, only: ppd_strategy_determine
    use rundat, only: pub_coulomb_cutoff, pub_psinc_spacing
    use simulation_cell, only: simulation_cell_initialise, &
         simulation_cell_add_padding
    use services, only: services_rationalise_coords
    use rundat, only: pub_eda, pub_frag_iatm
    use fragment_data, only: pub_frag_data
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(POINT),   intent(in   ) :: a1, a2, a3
    type(POINT),   intent(inout) :: a1_pad, a2_pad, a3_pad
    type(MODEL), intent(inout) :: mdl

    ! Local Variables
    real(kind=DP) :: d1, d2, d3
    real(kind=DP) :: r_abs(3,mdl%nat)
    integer       :: iat, loc_isp, ii
    integer       :: n_pt1, n_pt2, n_pt3
    integer       :: row, iatspec, rowspec
    logical       :: new_species
    integer       :: it, ierr, isub, iregion
    integer       :: jat(mdl%nsub)

    ! cks: determine grids and ppds based on KE-cutoff input parameter
    call ppd_strategy_determine(d1,d2,d3,n_pt1,n_pt2,n_pt3, &
         a1,a2,a3,mdl%elements)

    ! Initialise mdl%cell
    call simulation_cell_initialise(mdl%cell,a1,a2,a3, &
         n_pt1,n_pt2,n_pt3,d1,d2,d3)

    ! ndmh: initialise mdl%padded_cell
    if (pub_coulomb_cutoff) then

       ! ndmh: apply override to psinc spacing so that spacing in padded
       ! ndmh: cell matches that of original cell
       pub_psinc_spacing(1) = d1
       pub_psinc_spacing(2) = d2
       pub_psinc_spacing(3) = d3

       ! ndmh: check padding
       call simulation_cell_add_padding(d1,d2,d3,n_pt1,n_pt2,n_pt3, &
            a1,a2,a3,a1_pad,a2_pad,a3_pad)

       ! ndmh: determine grids and ppds of mdl%padded_cell
       call ppd_strategy_determine(d1, d2, d3, n_pt1, n_pt2, n_pt3, &
            a1_pad, a2_pad, a3_pad, mdl%elements)

       ! ndmh: initialise mdl%padded_cell
       call simulation_cell_initialise(mdl%padded_cell,a1_pad,a2_pad,a3_pad, &
            n_pt1,n_pt2,n_pt3,d1,d2,d3)

    end if

    ! mjsp: If this calculation uses fragments:
    if (pub_eda) then

       ! mjsp: Allocate fragment and supermolecule data workspace
       allocate(pub_frag_data(0:pub_frag_iatm(0)),stat=ierr) ! mjsp: allocate number
       call utils_alloc_check('rundat_blocks_init_grid','pub_frag_data',ierr)

       ! mjsp: Initialise the fragment cells
       do it=0,pub_frag_iatm(0)

          ! mjsp: Use same cell for all fragments:
          pub_frag_data(it)%mdl%cell = mdl%cell

          ! mjsp: initialise mdl%padded_cell if used
          if (pub_coulomb_cutoff) then

             ! mjsp: initialise mdl%padded_cell
             call simulation_cell_initialise(pub_frag_data(it)%mdl%padded_cell,&
                  a1_pad,a2_pad,a3_pad, &
                  n_pt1,n_pt2,n_pt3,d1,d2,d3)

          end if

       end do

    end if

    ! rc2013: cycle over all regions
    do isub=1,mdl%nsub

       ! cks: initialise total number of different atomic pseudopotential species
       ! cks, 25/1/2004: A distinct "species" is determined by the
       ! cks: pseudopotential name rather than the atomic number
       ! ndmh: also allow species to differ based on core wvfn name for EELS
       mdl%regions(isub)%par%num_pspecies = 1
       do iat=2,mdl%regions(isub)%par%nat

          new_species = .true.
          do row=1,iat-1
             iatspec  = mdl%regions(isub)%elements(iat)%species_number
             rowspec  = mdl%regions(isub)%elements(row)%species_number
             if ((mdl%regions(isub)%species(iatspec)&
                  %pseudo_name == &
                  mdl%regions(isub)%species(rowspec)&
                  %pseudo_name) .and. &
                  (mdl%regions(isub)%species(iatspec)&
                  %core_wf_name == &
                  mdl%regions(isub)%species(rowspec)&
                  %core_wf_name)) &
                  new_species = .false.
          end do

          if (new_species) mdl%regions(isub)%par%num_pspecies = &
               mdl%regions(isub)%par%num_pspecies + 1

       end do

       ! qoh: initialise total number of different atomic species
       ! qoh: (as listed in % block species).

       mdl%regions(isub)%par%num_species = 1
       do iat=2,mdl%regions(isub)%par%nat

          new_species = .true.
          do row=1,iat-1
             if (mdl%regions(isub)%elements(iat)%species_number == &
                  mdl%regions(isub)%elements(row)%species_number) &
                  new_species = .false.
          end do

          if (new_species) mdl%regions(isub)%par%num_species = &
               mdl%regions(isub)%par%num_species + 1

       end do

    enddo

    ! jcap: do the same for the whole system
    ! rc2013: place num_pspecies directly into mdl
    mdl%num_pspecies = 1
    do iat=2,mdl%nat

       new_species = .true.
       do row=1,iat-1
          if ((mdl%species(mdl%elements(iat)%species_number)&
               %pseudo_name == &
               mdl%species(mdl%elements(row)%species_number)&
               %pseudo_name) .and. &
              (mdl%species(mdl%elements(iat)%species_number)&
               %core_wf_name == &
               mdl%species(mdl%elements(row)%species_number)&
               %core_wf_name)) &
               new_species = .false.
       end do

       if (new_species) mdl%num_pspecies = mdl%num_pspecies + 1

    end do

    mdl%num_species = 1
    do iat=2,mdl%nat

       new_species = .true.
       do row=1,iat-1
          if (mdl%elements(iat)%species_number == mdl%elements(row)%species_number) &
               new_species = .false.
       end do

       if (new_species) mdl%num_species = mdl%num_species + 1

    end do

    ! smmd: ensure all ions are contained within the cell
    do iat=1,mdl%nat
       r_abs(1,iat) = mdl%elements(iat)%centre%x
       r_abs(2,iat) = mdl%elements(iat)%centre%y
       r_abs(3,iat) = mdl%elements(iat)%centre%z
    enddo

    call services_rationalise_coords(mdl%nat,r_abs,mdl%cell)

    jat = 0
    do iat=1,mdl%nat
       mdl%elements(iat)%centre%x =  r_abs(1,iat)
       mdl%elements(iat)%centre%y =  r_abs(2,iat)
       mdl%elements(iat)%centre%z =  r_abs(3,iat)
       ! rc2013: update the local system counters too
       iregion = mdl%elements(iat)%region
       ! rc2013: count the number of atoms in this region
       jat(iregion) = jat(iregion) + 1
       if (mdl%regions(iregion)%elements(jat(iregion))%global_atom_number &
            == iat) then
          mdl%regions(iregion)%elements(jat(iregion))%centre%x =  r_abs(1,iat)
          mdl%regions(iregion)%elements(jat(iregion))%centre%y =  r_abs(2,iat)
          mdl%regions(iregion)%elements(jat(iregion))%centre%z =  r_abs(3,iat)
       endif
    enddo

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: exit rundat_blocks_init_grid'

  end subroutine rundat_blocks_init_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module rundat_blocks
