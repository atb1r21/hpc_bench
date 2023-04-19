! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                    T S S E A R C H                                          !
!=============================================================================!
!                                                                             !
! $Id: tssearch_mod.F90,v 1.11 2009/09/23 16:56:58 cks22 Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! This module performs a transition state search starting from a reactant,    !
! product and intermediate (if required) structures                           !
!                                                                             !
!-----------------------------------------------------------------------------!
! Written from  "Module Specification for New Plane Wave Code",  M. Segall,   !
!    P. Lindan, M. Probert, C. Pickard, P. Hasnip, S. Clark and M. Payne      !
!                           Copyright 1999/2000                               !
!-----------------------------------------------------------------------------!
! Written by Niri Govind, v0.1, 19/11/2001                                    !
!-----------------------------------------------------------------------------!
! modification information                                                    !
!=============================================================================!

module tssearch
  use simulation_cell, only: castep_cell

  implicit none                                 !Impose strong typing

  private                                       !Everything is private ...

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: tssearch_run
  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!


  integer, save :: num_ionic_constraints

  ! Cell data in CASTEP format
  type(castep_cell)  :: current_cell

contains

  subroutine tssearch_run(mdl,output_file)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   cmdl, intent=inout, the model to be optimised                         !
    !   output_file, intent=in, the filename to which results will be written.!
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   cmdl has already been properly read and initialised                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, 20/11/2001                                !
    ! Modified to handle subsystem structures by Robert Charlton, 22/08/2019. !
    !=========================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout
    use model_type, only: MODEL
    use rundat, only: tssearch_method
    use simulation_cell, only: copy_to_castep_cell, &
         castep_model_dealloc, castep_cell_dealloc, &
         castep_model_alloc, castep_cell_copy, castep_cell_2_elements
    use tssearch_utils, only: tssearch_utils_io_abort, constrain_ions, &
         tssearch_utils_initialize
    use utils, only: utils_assert
    use services, only: services_flush

    implicit none

    ! Arguments
    type(MODEL), intent(inout) :: mdl
    character(len=*),  intent(in)    :: output_file

    integer :: ndim  !the number of dimensions in the search space
                   !so for the augmented-Hessian BFGS method used here,
                   !we have ndim=9+3*num_ions as we seek to optimize the
                   !cell (9 vars) and the ions (3*num_ions d.of.f) simultaneously.
    integer       :: ii,ireg,jat(mdl%nsub)

    !kkbd: NEB does not use any tssearch code
    if(tssearch_method == "NEB")then
       if (pub_on_root) write(stdout,*) "Starting NEB"
       call services_flush()
       call tssearch_neb_driver(mdl, output_file)
       return
    end if

    ! qoh: Actually initialise ndim!
    ndim = 9 + 3*mdl%nat

    ! fill CASTEP-style cell data structures from mdl%cell and "elements"
    ! kaw: Directly analogous to cks changes for energy calculation
    if (mdl%nat_classical > 0) then
       call copy_to_castep_cell(current_cell,mdl%cell,mdl%nat,mdl%elements, &
            mdl%nat_classical,mdl%classical_elements)
       write(stdout,*) 'Allocating for classical atoms',mdl%nat_classical
    else
       call copy_to_castep_cell(current_cell,mdl%cell,mdl%nat,mdl%elements, &
            mdl%nat_classical)
    end if

    ! Consistency check
    call utils_assert(current_cell%num_ions > 0, 'Error in tssearch_run(): &
         &current_cell uninitialised.')

    ! Allocate model data
    call castep_model_alloc(mdl%cmdl,mdl%nat+mdl%nat_classical,ndim, &
         constrain_ions,.false.)

    call castep_cell_copy(current_cell,mdl%cmdl%cell)

    !set the cmdl_ptr and output_file_ptr
    call tssearch_utils_initialize(mdl%cmdl,mdl)

    !Analyse constraints to see what is being asked to relax
    constrain_ions=.false.

    num_ionic_constraints=min(3*mdl%cmdl%cell%num_ions,num_ionic_constraints)  !catch any sillies
    num_ionic_constraints=max(0,num_ionic_constraints)                    !catch any sillies
    if (num_ionic_constraints>0) constrain_ions=.true.
    if (num_ionic_constraints==3*mdl%cmdl%cell%num_ions) call tssearch_utils_io_abort &
       & ('ERROR: Cannot perform transition state search run when all atoms are fixed')

    !Write main header
    if (pub_on_root) then
      write(stdout,1)
      write(stdout,2) 'Starting ONETEP Transition State Search'
      write(stdout,1)
      write(stdout,'(/a,a,a)') &
            ' Transition State Search: output file is "',trim(output_file),'"'
    end if


    !Now select transition state search method
    select case (tssearch_method)
    case ('LSTQST')
       call tssearch_LSTQST(mdl,output_file)
    case default
       call tssearch_utils_io_abort('Error in tssearch_run - unrecognised transition state search request='//tssearch_method//'.')
    end select

    ! Make sure that the elements ONETEP array is updated
    call castep_cell_2_elements(mdl%cmdl%cell,mdl%elements,mdl%nat)

    ! rc2013: update ionic coordinates in regions lists
    jat = 0
    do ii=1,mdl%nat
       ! Get the region counter
       ireg = mdl%elements(ii)%region
       ! Count the no. of atoms allocated to each region
       jat(ireg) = jat(ireg) + 1
       ! Copy the co-ordinates over
       mdl%regions(ireg)%elements(jat(ireg))%centre = mdl%elements(ii)%centre
    end do

    ! Deallocate memory
    call castep_model_dealloc(mdl%cmdl,.false.)
    call castep_cell_dealloc(current_cell)

1   format(80('='))
2   format(19('<'),1x,a,1x,20('>'))

    return
  end subroutine tssearch_run

  subroutine tssearch_neb_driver(mdl,output_file)
    use constants, only: dp, stdout
    use rundat, only: pub_devel_code
    use model_type, only: model
    use neb, only: neb_path, neb_init, neb_energy_and_force, neb_optimize
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_devel_code
    use geometry_optimiser, only: geometry_optimise

    implicit none

    integer :: ierr
    logical :: use_geomopt

    type(NEB_PATH) :: neb_mep

    type(model),intent(inout) :: mdl
    character(len=*),  intent(in) :: output_file

    real(kind=dp) :: image_energy
    real(kind=dp),allocatable,dimension(:,:) :: image_force

    allocate(image_force(1:3,mdl%nat),stat=ierr)
    call utils_alloc_check('neb_driver', 'image force', ierr)

    call neb_init(mdl,neb_mep)

    ! kkbd: We can use BFGS or the NEB internal optimizer
    use_geomopt = utils_devel_code(.false.,"NEB","USE_GEOMOPT", pub_devel_code)
    if (.not.use_geomopt) then
      call neb_optimize(mdl,neb_mep)
    else
      call geometry_optimise(image_energy, image_force, &
        mdl, output_file, path=neb_mep)
    end if

    deallocate(image_force,stat=ierr)
    call utils_dealloc_check('neb_driver', 'image force', ierr)

  end subroutine tssearch_neb_driver

  subroutine tssearch_LSTQST(mdl,output_file)
    !=========================================================================!
    ! Initializes the reactant, product and intermediate cells                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   I/O                                                                   !
    !   Wavefunction                                                          !
    !   Density                                                               !
    !   Cell                                                                  !
    !   Parameters                                                            !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, version 0.1, 21/11/01                           !
    !=========================================================================!

    use comms,    only : pub_on_root,comms_bcast,pub_root_proc_id
    use constants,only : DP
    use esdf,     only : esdf_init, esdf_close
    use ion,      only : element
    use model_type, only: MODEL
    use rundat,   only : pub_rootname, tssearch_lstqst_protocol, &
         tssearch_cg_max_iter, tssearch_force_tol, &
         tssearch_disp_tol, tssearch_qst_max_iter
    use rundat_blocks, only: rundat_blocks_exec
    use simulation_cell, only: castep_cell_dealloc, &
         castep_cell_alloc, castep_cell_copy, castep_cell_frac_to_cart, &
         copy_to_castep_cell
    use tssearch_algor,  only: tssearch_algor_QST, tssearch_algor_lst
    use tssearch_utils,  only: tssearch_utils_io_abort, &
         tssearch_utils_get_coordinates, tssearch_utils_fix_coordinates, &
         tssearch_utils_check_coord, tssearch_utils_deallocate
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments

    type(MODEL), intent(inout) :: mdl ! The model to be initialised
    character(len=*),  intent(in)     :: output_file

    ! Local Variables
    type(ELEMENT),allocatable,dimension(:) :: elements_tmp
    integer :: ndim_atoms
    integer :: nat
    integer :: nat_reac
    integer :: nat_prod
    integer :: nat_diff
    integer :: ierr
    real (kind=dp), allocatable, dimension(:,:,:) :: xbuff
    real (kind=dp), allocatable, dimension(:)     :: x_reac,x_prod,x_intm,xclean1,xclean2
    real (kind=dp), allocatable, dimension(:)     :: g_reac,g_prod,g_intm
    real (kind=dp)                                :: e_reac,e_prod,e_intm
    real (kind=dp)                                :: p12,p13,p23
    integer,  allocatable, dimension(:)           :: idfix
    logical                                       :: tsfound,period,impose
    integer                                       :: neg,nfix
    integer                                       :: mtsmethod,nstrucs
    integer                                       :: icontinue

    ! First initialize the cells

    ! Now, allocate the original, model cell, reactant, product and intermediate cells
    ! and copy some dummy contents
!    call castep_cell_alloc(mdl%cmdl%orig_cell,mdl%nat+mdl%nat_classical)
!    call castep_cell_alloc(mdl%cmdl%cell,mdl%nat+mdl%nat_classical)
!    call castep_cell_alloc(mdl%cmdl%reac_cell,mdl%nat+mdl%nat_classical)
!    call castep_cell_alloc(mdl%cmdl%prod_cell,mdl%nat+mdl%nat_classical)
!    call castep_cell_alloc(mdl%cmdl%intm_cell,mdl%nat+mdl%nat_classical)

    call castep_cell_copy(current_cell,mdl%cmdl%orig_cell)
    call castep_cell_copy(current_cell,mdl%cmdl%cell)
    call castep_cell_copy(current_cell,mdl%cmdl%reac_cell)
    call castep_cell_copy(current_cell,mdl%cmdl%prod_cell)
    call castep_cell_copy(current_cell,mdl%cmdl%intm_cell)

    ! Copy relevant coordinates into appropriate locations
    nstrucs = 0

    !check model has all the bits needed
    if (.not.allocated(mdl%cmdl%forces)) then
       allocate(mdl%cmdl%forces(1:3,1:mdl%cmdl%cell%max_ions_in_species,1:mdl%cmdl%cell%num_species),stat=ierr)
       call utils_alloc_check('tssearch_LSTQST','mdl%cmdl%forces',ierr)
       mdl%cmdl%forces=0.0_dp
    end if

    ! Reactant structure - already copied (take current one)
    nstrucs = nstrucs + 1
    nat_reac  = mdl%cmdl%reac_cell%num_ions

    !Calculate number of degrees of freedom
    nat = nat_reac
    ndim_atoms = 3*nat

    ! prepare a spare array for reading configurations
    allocate(elements_tmp(nat),stat=ierr)
    call utils_alloc_check('tssearch_LSTQST','elements_tmp',ierr)

    if (pub_on_root) call esdf_init(trim(pub_rootname)//'.dat',ierr)

    ! store initial elements array
    elements_tmp = mdl%elements

    ! Product structure
    nstrucs = nstrucs + 1
    call rundat_blocks_exec(mdl,'PRODUCT')
    call copy_to_castep_cell(mdl%cmdl%prod_cell,mdl%cell,mdl%nat,mdl%elements, &
         mdl%nat_classical,mdl%classical_elements)
    nat_prod  = mdl%cmdl%prod_cell%num_ions

    ! Intermediate structure
    select case (tssearch_lstqst_protocol)
    case ('QST/OPTIMIZATION')
        call rundat_blocks_exec(mdl,'INTERMEDIATE')
        call copy_to_castep_cell(mdl%cmdl%intm_cell,mdl%cell,mdl%nat, &
             mdl%elements,mdl%nat_classical,mdl%classical_elements)
        nstrucs = nstrucs + 1
    end select

    ! restore initial elements array from copy
    mdl%elements = elements_tmp

    if (pub_on_root) call esdf_close()

    !Check number of atoms in reactant and product models
    nat_diff = nat_reac - nat_prod
    if (nat_diff /=0) ierr = 1
    if (ierr/=0) call tssearch_utils_io_abort('Error in tssearch_LSTQST - Unequal number of atoms &
            & in reactant and product')

    !Allocate reactant, product, intermediate arrays
    !Reactant
    ndim_atoms = (3*nat) - (3*mdl%nat_classical)
    allocate(x_reac(1:ndim_atoms),stat=ierr)
    call utils_alloc_check('tssearch_LSTQST','x_reac',ierr)
    x_reac=0.0_dp
    allocate (g_reac(1:ndim_atoms),stat=ierr)
    call utils_alloc_check('tssearch_LSTQST','g_reac',ierr)
    g_reac=0.0_dp

    !Product
    allocate (x_prod(1:ndim_atoms),stat=ierr)
    call utils_alloc_check('tssearch_LSTQST','x_prod',ierr)
    x_prod=0.0_dp
    allocate (g_prod(1:ndim_atoms),stat=ierr)
    call utils_alloc_check('tssearch_LSTQST','g_prod',ierr)
    g_prod=0.0_dp

    !Intermediate
    allocate (x_intm(1:ndim_atoms),stat=ierr)
    call utils_alloc_check('tssearch_LSTQST','x_intm',ierr)
    x_intm=0.0_dp
    allocate (g_intm(1:ndim_atoms),stat=ierr)
    call utils_alloc_check('tssearch_LSTQST','g_intm',ierr)
    g_intm=0.0_dp

    !Arrays for clean structures
    allocate (xclean1(1:ndim_atoms),stat=ierr)
    call utils_alloc_check('tssearch_LSTQST','xclean1',ierr)
    xclean1=0.0_dp
    allocate (xclean2(1:ndim_atoms),stat=ierr)
    call utils_alloc_check('tssearch_LSTQST','xclean2',ierr)
    xclean2=0.0_dp

    !Buffer variable coordinates
    allocate(xbuff(1:3,1:mdl%cmdl%cell%max_ions_in_species,1:mdl%cmdl%cell%num_species),stat=ierr)
    call utils_alloc_check('tssearch_LSTQST','xbuff',ierr)
    xbuff=0.0_dp

    !Fixed coordinate information
    allocate (idfix(1:ndim_atoms),stat=ierr)
    call utils_alloc_check('tssearch_LSTQST','idfix',ierr)
    idfix=0
    nfix=0

    !Reactant
    call castep_cell_frac_to_cart(mdl%cmdl%reac_cell,xbuff)
    call tssearch_utils_get_coordinates(mdl%cmdl,xbuff,x_reac)
    xclean1 = x_reac      !preserve a clean set of reactant coordinates

    !Product
    call castep_cell_frac_to_cart(mdl%cmdl%prod_cell,xbuff)
    call tssearch_utils_get_coordinates(mdl%cmdl,xbuff,x_prod)
    xclean2 = x_prod      !preserve a clean set of product coordinates

    !Intermediate
    call castep_cell_frac_to_cart(mdl%cmdl%intm_cell,xbuff)
    call tssearch_utils_get_coordinates(mdl%cmdl,xbuff,x_intm)

    !get fixed coordinate information
    call tssearch_utils_fix_coordinates(nfix,idfix)

    !check the consistency of the structures with fixed atom information and in relation to each other
    call tssearch_utils_check_coord(ndim_atoms,x_reac,x_prod,p12,idfix)

    !set necessary flags before calling main drivers
    impose = .FALSE.        !no superposition for periodic models
    period = .TRUE.         !periodic boundary condition flag
    tsfound = .FALSE.       !transition state flag
    neg = 0                 !energy-gradient counter

    !set the method flags
    select case (tssearch_lstqst_protocol)
    case ('LSTMAXIMUM')
            mtsmethod = 1
            if (nstrucs < 2) call tssearch_utils_io_abort('Error in tssearch_LSTQST - &
            & Insufficient number of structures provided. Please check the seed.dat file')
    case ('HALGREN-LIPSCOMB')
            mtsmethod = 12
            if (nstrucs < 2) call tssearch_utils_io_abort('Error in tssearch_LSTQST - &
            & Insufficient number of structures provided. Please check the seed.dat file')
    case ('LST/OPTIMIZATION')
            mtsmethod = 12
            if (nstrucs < 2) call tssearch_utils_io_abort('Error in tssearch_LSTQST - &
            & Insufficient number of structures provided. Please check the seed.dat file')
    case ('COMPLETELSTQST')
            mtsmethod = 123
            if (nstrucs < 2) call tssearch_utils_io_abort('Error in tssearch_LSTQST - &
            & Insufficient number of structures provided. Please check the seed.dat file')
    case ('QST/OPTIMIZATION')
            mtsmethod = 3
            if (nstrucs /= 3) call tssearch_utils_io_abort('Error in tssearch_LSTQST - &
            & Insufficient number of structures provided. Please check the seed.dat file')
    !check the consistency of the structures with fixed atom information and in relation to each other
            call tssearch_utils_check_coord(ndim_atoms,x_reac,x_intm,p13,idfix)
            call tssearch_utils_check_coord(ndim_atoms,x_prod,x_intm,p23,idfix)
    end select

    !select driver
    icontinue=0
    if ( mtsmethod .eq. 1  .or. mtsmethod .eq. 12 .or.  mtsmethod .eq. 123 ) then
        if(pub_on_root)then
           ! Do algorithimic things only on the root proc
           ! all other procs are "on alert" in the  tssearch_dummy
           ! where the are waiting for the signal to join electronic calculations
           call tssearch_algor_LST(mtsmethod, nat - mdl%nat_classical, ndim_atoms, &
                                period, impose, &
                                x_reac, g_reac, e_reac, x_prod, g_prod, e_prod, &
                                x_intm, g_intm, e_intm, xclean1,xclean2, nfix, idfix,&
                                tsfound,neg,output_file,&
                                tssearch_cg_max_iter,tssearch_force_tol,&
                                tssearch_disp_tol)
           ! signal to let other procs to run:
           call comms_bcast(pub_root_proc_id,icontinue)

        else
           call tssearch_dummy(ndim_atoms,x_reac, g_reac)
        endif
        ! let's inform all the procs about current status:
        call comms_bcast(pub_root_proc_id,tsfound)

    end if

    if (.not.tsfound ) then
       if ( mtsmethod .eq. 123 .or.  mtsmethod .eq. 3   ) then
        if(pub_on_root)then
           ! Do algorithimic things only on the root proc
           ! all other procs are "on alert" in the  tssearch_dummy
           ! where the are waiting for the signal to join electronic calculations
           call tssearch_algor_QST(mtsmethod, nat, ndim_atoms, period, impose, &
                                x_reac, g_reac, e_reac, x_prod, g_prod, e_prod, &
                                x_intm, g_intm, e_intm, xclean1,xclean2, nfix, idfix,&
                                tsfound,neg,output_file,tssearch_qst_max_iter,&
                                tssearch_cg_max_iter,tssearch_force_tol,&
                                tssearch_disp_tol)
           ! signal to let other procs to run:
           call comms_bcast(pub_root_proc_id,icontinue)
        else
           call tssearch_dummy(ndim_atoms,x_reac, g_reac)
        endif

       end if
    end if

    !Deallocate local arrays
    deallocate (x_reac,stat=ierr)
    call utils_dealloc_check('tssearch_LSTQST','x_reac',ierr)
    deallocate (g_reac,stat=ierr)
    call utils_dealloc_check('tssearch_LSTQST','g_reac',ierr)
    deallocate (x_prod,stat=ierr)
    call utils_dealloc_check('tssearch_LSTQST','x_prod',ierr)
    deallocate (g_prod,stat=ierr)
    call utils_dealloc_check('tssearch_LSTQST','g_prod',ierr)
    deallocate (x_intm,stat=ierr)
    call utils_dealloc_check('tssearch_LSTQST','x_intm',ierr)
    deallocate (g_intm,stat=ierr)
    call utils_dealloc_check('tssearch_LSTQST','g_intm',ierr)
    deallocate (xbuff,stat=ierr)
    call utils_dealloc_check('tssearch_LSTQST','xbuff',ierr)
    deallocate (idfix,stat=ierr)
    call utils_dealloc_check('tssearch_LSTQST','idfix',ierr)
    deallocate (xclean1,stat=ierr)
    call utils_dealloc_check('tssearch_LSTQST','xclean1',ierr)
    deallocate (xclean2,stat=ierr)
    call utils_dealloc_check('tssearch_LSTQST','xclean2',ierr)
    deallocate(elements_tmp,stat=ierr)
    call utils_dealloc_check('tssearch_LSTQST','elements_tmp',ierr)

    !all done
    !Deallocate public arrays in tssearch_utils
    call tssearch_utils_deallocate()

    return
  end subroutine tssearch_LSTQST

  subroutine tssearch_dummy(ndim,x, g)
    !=========================================================================!
    ! This routine is called on each non root proc. An execution stops here   !
    ! and program wait for the data from the root proc to start force         !
    ! calculations for the geometry defined by the TSSEARCH algorithm (root   !
    ! proc). When data are ready variable icontinue=1 is sent to each proc.   !
    ! After everything is done icontinue=0 is sent and routine return control !
    ! the main flow                                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:(see below)                                                   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used: constants,comms,tssearch_utils                            !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   cmdl has already been properly read and initialised                    !
    !-------------------------------------------------------------------------!
    ! Written by Alexander Perlov 20/Dec/2006                                 !
    !=========================================================================!
    use constants, only: DP
    use comms, only : comms_bcast,pub_root_proc_id
    use tssearch_utils, only : tssearch_utils_energy_gradient
    implicit none
    ! arguments
    integer, intent(in)           ::ndim         ! number of atom dimensions
    real(kind=dp), intent(inout)  ::x(:)         ! Cartesian coords array
    real(kind=dp), intent(inout)  ::g(:)         ! Gradients array


    !local variables
    real(kind=dp)                 ::energy   ! Total energy
    integer                       :: icontinue
    call comms_bcast(pub_root_proc_id,icontinue)
    do while (icontinue == 1)
       call tssearch_utils_energy_gradient(ndim,x,g,energy)
       call comms_bcast(pub_root_proc_id,icontinue)
    end do
  end  subroutine tssearch_dummy

end module tssearch
