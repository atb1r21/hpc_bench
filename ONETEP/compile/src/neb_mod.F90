! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!-----------------------------------------------------------------------!
!                 N U D G E D   E L A S T I C   B A N D                 !
!-----------------------------------------------------------------------!
!                                                                       !
! This module implements the nudged elastic band method for transition  !
! state searching. It makes use of the image comms functionality such   !
! that each NEB bead is a separate ONETEP image. Chain optimization is  !
! threaded through the geometry optimizer which makes calls to NEB      !
! force nudging routines.                                               !
!                                                                       !
! To use the functionality, include the below in your input file,       !
! where N is the number of beads you want and the number of MPI         !
! processes you want must be divisible by N (unless you use advanced    !
! image specification functionality). You must also specify a product   !
! positions block using, eg, %block positions_abs_product, where the    !
! order of the product atoms will define the linear interpolation used  !
! between reactant and product as well as which atoms are chained       !
! together. Climbing image can be switched on after M geometry          !
! optimization steps, and M=-1 (default) disables climbing image.       !
! Convergence tolerances can be specified using the tssearch tol        !
! keywords (default atomic units)                                       !
!                                                                       !
! Each image has its own output files, which can be used for            !
! continuation (e.g. "geom_continuation : T" continues an existing      !
! NEB run).                                                             !
!                                                                       !
! Currently, reactant and product energies are provided through devel   !
! codes. If they're not provided a simple approach to make the tangent  !
! well-behaved will be performed.                                       !
!                                                                       !
!-----------------------------------------------------------------------!
!                                                                       !
! To enable a NEB run:                                                  !
!   task : transitionstatesearch                                        !
!   tssearch_method : NEB                                               !
!   num_images : N                                                      !
!   %block positions_{abs, frac}_product                                !
!     [product specification]                                           !
!   %endblock positions_{abs, frac}_product                             !
! Optional:                                                             !
!   neb_ci_delay : M                                                    !
!   tssearch_energy_tol : q                                             !
!   tssearch_disp_tol   : r                                             !
!   tssearch_force_tol  : s                                             !
! Devel Codes:                                                          !
!   Block: NEB                                                          !
!   Real keywords: Reactant_Energy, Product_Energy (both in Ha)         !
!                                                                       !
!-----------------------------------------------------------------------!
!                                                                       !
! Written by Kevin Duff mostly in 2017 and 2018.                        !
!                                                                       !
!-----------------------------------------------------------------------!
! To-do list:                                                           !
! - Add a references section here and proper citations to the function  !
!   calls                                                               !
! - Investigate/develop dynamic spring constant methods                 !
! - Allow user to specify multiple beads per ONETEP image               !
! - Use linear interpolation to allow users to change number of beads   !
!   on continuation.                                                    !
!-----------------------------------------------------------------------!

module neb

  use constants,  only: DP, stdout
  use model_type, only: MODEL

  implicit none

  private

  ! Types

  type NEB_PATH

    ! Positions and nudged forces of product, images, and reactant
    real(kind=DP), allocatable, dimension(:,:) ::  im_positions,  im_forces, &
                                                  old_positions, old_forces

    ! Energies of product, images, and reactant
    real(kind=DP), allocatable, dimension(:) :: im_energies
    ! NOTE: we could wrap these into im_energies similar to im_positions
    real(kind=DP) :: reactant_energy, product_energy

    ! Spring constant
    ! A potential extension of this would be to allow each spring to vary its
    ! spring constant to change sampling across the path.
    real(kind=DP) :: spring_k

    ! Image index of climbing image - less than 0 if none
    integer :: climbing_image

    ! Information used by each process, not just the root process
    real(kind=dp) :: prev_energy, next_energy
    real(kind=dp),allocatable,dimension(:) :: prev_pos, next_pos, &
                                              prev_to_cur, cur_to_next, &
                                              my_force
    integer :: nat
    integer :: ci_delay

    ! My image's total energy
    real(kind=dp) :: my_energy

    ! Optimization methods variables
    character(len=80) :: update_method

    ! FIRE-related
    integer :: cut, n_min
    real(kind=dp) :: dt_max, fire_adec, fire_dec, fire_inc, fire_alpha, orig_alpha, fire_dt

    ! LBFGS-related
    real(kind=dp), allocatable, dimension(:)   :: lbfgs_q, lbfgs_z, old_q, old_z ! (1:3*nat*num_images)
    real(kind=dp), allocatable, dimension(:)   :: lbfgs_rho, lbfgs_alpha         ! (1:history_size)
    real(kind=dp), allocatable, dimension(:,:) :: lbfgs_s, lbfgs_y               ! (1:3*nat*num_images, 1:history_size)

  end type NEB_PATH

  public :: neb_path

  ! Private Variables

  ! Public Stuff

  public :: neb_init
  public :: neb_energy_and_force
  public :: neb_sync_path
  public :: neb_converged
  public :: neb_optimize

  logical, public :: neb_calculated_endpoint

  contains

  subroutine neb_init(mdl, mep)
    !====================================================================!
    ! Set up the NEB chain.                                              !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  mdl (inout) - Model corresponding to the reactant                 !
    !  mep (  out) - Structure containing information for the whole MEP  !
    !--------------------------------------------------------------------!
    ! Necessary conditions:                                              !
    !  - Image comms correctly set up                                    !
    !  - Product structure defined                                       !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff                                              !
    !====================================================================!
    use comms, only: pub_on_root, pub_image_comm, pub_imroots_comm, &
         comms_send, comms_recv, comms_bcast, comms_barrier, pub_world_comm, &
         pub_root_proc_id
    use energy_and_force, only: energy_and_force_calculate
    use esdf, only: esdf_init, esdf_close, esdf_block
    use image_comms, only: image_comms_init, pub_my_rank_in_image, &
         pub_my_image, pub_my_rank_in_imroots, pub_on_imroots_root, &
         pub_my_rank_in_world, pub_imroots_root_id, pub_base_rootname, &
         orig_stdout
    use rundat, only:  pub_num_images, pub_write_xyz, pub_rootname, &
         pub_devel_code, pub_debug_on_root, pub_reactant_energy,  &
         pub_reactant_rootname, pub_product_energy, pub_product_rootname, &
         pub_input_xyz_file_intermediate, pub_neb_spring_constant, pub_neb_read_xyz, &
         pub_neb_ci_delay, pub_neb_update_method, pub_print_qc, pub_neb_glbfgs_history_size
    use rundat_blocks, only: rundat_blocks_exec
    use services, only: services_flush, services_write_xyz, services_read_xyz
    use utils, only: utils_alloc_check, utils_abort, utils_dealloc_check, &
         utils_devel_code, utils_assert
    implicit none

    ! Arguments
    type(MODEL),               intent(inout) :: mdl
    type(NEB_PATH),            intent(  out) :: mep

    ! Local Variables
    real(kind=DP), allocatable, dimension(:)   :: temp_pos, reac_pos, prod_pos,&
                                                  interm_pos
    real(kind=DP), allocatable, dimension(:,:) :: temp_force, posbuf
    character(len=4), allocatable, dimension(:) :: dummy_species
    real(kind=DP) :: ri_dist, ip_dist, total_dist
    integer :: ierr, i, product_image, interm_bead, dummy_nat, row
    logical :: interm_specified
    character(len=4) :: img_num
    character(len=24) :: outtag

    ! NOTE: Maybe we should output a warning if <6 DoF are fixed
    !   (e.g. fix_all_cell is 9 DoF
    !         2 fixed atoms
    !         fixed cell vector & fixed atom
    !         3 atoms with point, line, and plane constraints)
    ! For now we're just using fixed cells... Don't bother with constraints as > 6 DoF are fixed.

    mep%climbing_image = -1 !Start with no climbing image
    mep%nat = mdl%nat
    mep%ci_delay = pub_neb_ci_delay
    mep%update_method = pub_neb_update_method

    ! If on root start setting up the neb_path
    if (pub_on_imroots_root) then
      ! N.B. first and last are product, reactant
      allocate(mep%im_positions(1:3*mdl%nat,1:pub_num_images+2), stat=ierr)
      call utils_alloc_check('neb_init', 'mep%im_positions', ierr)
      allocate(mep%old_positions(1:3*mdl%nat,1:pub_num_images), stat=ierr)
      call utils_alloc_check('neb_init', 'mep%old_positions', ierr)
      allocate(mep%im_forces(1:3*mdl%nat,1:pub_num_images), stat=ierr)
      call utils_alloc_check('neb_init', 'mep%im_forces', ierr)
      allocate(mep%old_forces(1:3*mdl%nat,1:pub_num_images), stat=ierr)
      call utils_alloc_check('neb_init', 'mep%old_forces', ierr)
      allocate(mep%im_energies(1:pub_num_images), stat=ierr)
      call utils_alloc_check('neb_init', 'mep%im_energies', ierr)

      mep%im_positions = 0.0_dp
      mep%old_positions = 0.0_dp
      mep%im_forces = 0.0_dp
      mep%old_forces = 0.0_dp
      mep%im_energies = 0.0_dp
      mep%reactant_energy = 0.0_dp
      mep%product_energy = 0.0_dp
    end if

    ! LBFGS matrices
    ! We could get away with only allocating these matrices where they're used
    ! but we can switch optimization method on the fly so I just allocate them anyway.
    allocate(mep%lbfgs_s(3*mep%nat*pub_num_images, pub_neb_glbfgs_history_size), &
             mep%lbfgs_y(3*mep%nat*pub_num_images, pub_neb_glbfgs_history_size), &
             mep%lbfgs_q(3*mep%nat*pub_num_images), &
             mep%lbfgs_z(3*mep%nat*pub_num_images), &
             mep%old_q(3*mep%nat*pub_num_images), &
             mep%old_z(3*mep%nat*pub_num_images), &
             mep%lbfgs_rho(pub_neb_glbfgs_history_size), &
             mep%lbfgs_alpha(pub_neb_glbfgs_history_size), stat=ierr)
    call utils_alloc_check("neb_init","G-LBFGS arrays",ierr)

    call utils_assert(pub_num_images > 1, "Error in neb_init: &
            &NEB currently only supports 2 or more beads. &
            &Currently running with num_images = ", &
            pub_num_images)

    neb_calculated_endpoint = .false.

    allocate(mep%prev_to_cur(1:3*mdl%nat), &
             mep%cur_to_next(1:3*mdl%nat), stat=ierr)
    call utils_alloc_check('neb_driver', 'neighbor vectors', ierr)
    allocate(mep%prev_pos(1:3*mdl%nat), mep%next_pos(1:3*mdl%nat), stat=ierr)
    call utils_alloc_check('neb_init', 'previous and next positions',ierr)
    allocate(mep%my_force(1:3*mdl%nat), stat=ierr)
    call utils_alloc_check('neb_init', "image's force",ierr)

    mep%prev_to_cur = 0.0_dp
    mep%cur_to_next = 0.0_dp
    mep%prev_pos = 0.0_dp
    mep%next_pos = 0.0_dp

    ! GL-BFGS parameters
    mep%lbfgs_q = 0.0_dp
    mep%lbfgs_z = 0.0_dp
    mep%lbfgs_rho = 0.0_dp
    mep%lbfgs_alpha = 0.0_dp
    mep%lbfgs_s = 0.0_dp
    mep%lbfgs_y = 0.0_dp

    ! FIRE parameters
    ! values taken from "Structural Relaxation Made Simple (FIRE reference)" or estimated
    mep%cut = 1
    mep%fire_alpha = 0.1_dp
    mep%orig_alpha = mep%fire_alpha
    mep%n_min = 5
    mep%fire_inc = 1.1_dp
    mep%fire_dec = 0.5_dp
    mep%fire_adec = 0.99_dp
    mep%fire_dt = 0.1_dp
    mep%dt_max = 1.0_dp

    ! Start setting up this image's model and info
    allocate(temp_pos(1:3*mdl%nat), temp_force(1:3, 1:mdl%nat), stat=ierr)
    call utils_alloc_check('neb_init', 'temp arrays', ierr)
    allocate(reac_pos(1:3*mdl%nat), prod_pos(1:3*mdl%nat), interm_pos(1:3*mdl%nat), stat=ierr)
    call utils_alloc_check('neb_init', 'reactant, product, and intermediate arrays', ierr)

    if (pub_num_images > 1) then
      product_image = pub_num_images-1
    else
      product_image = 0
    end if

    if (pub_on_imroots_root) then
      ! Reactant
      ! Handle current mdl, throw its info into the neb_path
      if (.not.pub_on_root) then
        call utils_abort("ERROR in image comms setup.")
      end if
      ! Extract reactant information
      call neb_extract_coords_from_mdl(mdl, mep%im_positions(:,1))
      reac_pos = mep%im_positions(:,1)

      mep%old_positions = 0.0_dp

      ! Product
      call esdf_init(trim(pub_base_rootname)//'.dat',ierr)
      if (ierr /= 0) then
        call utils_abort("Error reading input file "&
              //trim(pub_base_rootname)//".dat:",ierr)
      end if
    end if ! pub_on_imroots_root

    ! This must be done by the whole image
    if (pub_my_image == 0) then
      call comms_barrier
      call rundat_blocks_exec(mdl, 'PRODUCT')
      call comms_barrier
    end if

    if (pub_on_imroots_root) then
      ! Handle product mdl, throw its info into the neb_path
      call neb_extract_coords_from_mdl(mdl, mep%im_positions(:,pub_num_images+2))
    end if

    ! Send the product position to the product image
    if (product_image == 0) then
      ! We're doing both the reactant & product on just one image!
      prod_pos = mep%im_positions(:,pub_num_images+2)
    else
      ! Send the product info to the product image
      if (pub_on_imroots_root) then
        call comms_send(product_image,mep%im_positions(:,pub_num_images+2), &
                 length=3*mdl%nat, tag=1337, comm=pub_imroots_comm)
      else if ((pub_my_image == product_image) .and. pub_on_root) then
        call comms_recv(0,prod_pos,length=3*mdl%nat, &
                 tag=1337,comm=pub_imroots_comm)
      end if
    end if

    if (pub_my_image == 0) then
      call comms_barrier(pub_image_comm)
      call comms_bcast(pub_root_proc_id,reac_pos,length=3*mdl%nat,comm=pub_image_comm)
    end if

    if (pub_my_image == product_image) then
      call comms_barrier(pub_image_comm)
      call comms_bcast(pub_root_proc_id,prod_pos,length=3*mdl%nat,comm=pub_image_comm)
    end if

    ! Do a singlepoint on the product and reactant.
    if (pub_my_image == 0) then
      if (trim(adjustl(pub_reactant_rootname)) /= "NONE") then
        if (pub_on_root) write(stdout,*) &
             "NEB: Reading reactant energy from rootname: ", &
             pub_reactant_rootname
        ! Read in reactant tighbox NGWF & DKN and calculate total energy from that
        call neb_load_coords_into_mdl(mdl, reac_pos) !mep%im_positions(:,1))
        !!! COMPLEX, KPOINT
        call energy_and_force_calculate(mep%reactant_energy,temp_force,mdl, &
             custom_restart_name=trim(pub_reactant_rootname))

      else if (pub_reactant_energy < 0.0_DP) then
        if (pub_on_root) write(stdout,*) &
             "NEB: Reactant energy specified in input: ", pub_reactant_energy
        ! User specified reactant total energy directly
        mep%reactant_energy = pub_reactant_energy
      else
        ! Calculate reactant energy from scratch
        if (pub_on_root) write(stdout,*) &
             "NEB: Performing total energy calculation on reactant."

        call neb_load_coords_into_mdl(mdl, reac_pos) !mep%im_positions(:,1))
        !!! COMPLEX, KPOINT
        call energy_and_force_calculate(mep%reactant_energy,temp_force,mdl)
      end if

      if (pub_on_root) write(orig_stdout,*) "NEB: Reactant Energy - ", mep%reactant_energy
    end if

    ! Work out product energy
    if (pub_my_image == product_image) then
      if (trim(adjustl(pub_product_rootname)) /= "NONE") then
        if (pub_on_root) write(stdout,*) &
             "NEB: Reading product energy from rootname: ", pub_product_rootname
        ! Read in product tighbox NGWF & DKN and calculate total energy from that
        call neb_load_coords_into_mdl(mdl,prod_pos)
        !!! COMPLEX, KPOINT
        call energy_and_force_calculate(mep%product_energy,temp_force,mdl, &
             custom_restart_name=trim(pub_product_rootname))
      else if (pub_product_energy < 0.0_DP) then
        if (pub_on_root) write(stdout,*) &
             "NEB: Product energy specified in input: ", pub_product_energy
        ! User specified product total energy directly
        mep%product_energy = pub_product_energy
      else
        ! Calculate product energy from scratch
        if (pub_on_root) write(stdout,*) &
             "NEB: Performing total energy calculation on product."

        call neb_load_coords_into_mdl(mdl,prod_pos)
        !!! COMPLEX, KPOINT
        call energy_and_force_calculate(mep%product_energy,temp_force,mdl)
      end if
      if (pub_on_root) write(orig_stdout,*) "NEB: Product Energy - ", mep%product_energy
    end if

    if (pub_on_root) call comms_bcast(0, mep%reactant_energy, comm=pub_imroots_comm)
    call comms_bcast(pub_root_proc_id, mep%reactant_energy, comm=pub_image_comm)
    if (pub_on_root) call comms_bcast(product_image, mep%product_energy, comm=pub_imroots_comm)
    call comms_bcast(pub_root_proc_id, mep%product_energy, comm=pub_image_comm)

    ! Done working out reactant and product total energies

    ! Setup spring constant. For now, we just either read it in or calculate a sensible one.
    ! I think later on some work can be put into `smart' spring constants.
    if (pub_neb_spring_constant < 0.0_DP) then
       !mep%spring_k = neb_auto_spring_const()
       mep%spring_k = 0.02
    else
       ! User specified one in the input file, use that
       mep%spring_k = pub_neb_spring_constant
    end if

    if (pub_debug_on_root) write(orig_stdout,*)"NEB: Spring Constant - ",mep%spring_k

    if (pub_my_image == 0) then
       ! Initialize neb images in the neb_path type

       ! Linear interpolation of product and reactant.
       ! If intermediate image is supplied do 2 interpolations.
       ! Intermediate

       interm_specified = .false.
       if (pub_on_root) then
          interm_specified = (esdf_block('positions_abs_intermediate',i) .or. &
                              esdf_block('positions_frac_intermediate',i) .or. &
                              (trim(pub_input_xyz_file_intermediate)/=''))
       end if

       call comms_bcast(0, interm_specified, comm=pub_image_comm)

      if (pub_neb_read_xyz) then
         ! We're restarting from xyzs of a previous NEB run or just a user-specified path
         if (pub_on_root) write(orig_stdout,*)"NEB: Reading initial NEB beads from xyz files"

         if (pub_on_imroots_root) then
            allocate(posbuf(3,mdl%nat),dummy_species(mdl%nat),stat=ierr)
            call utils_alloc_check("neb_init","xyz read arrays",ierr)

            do i = 0, pub_num_images-1
              write(img_num,'(I4)')i
              call services_read_xyz(trim(pub_base_rootname)//trim(adjustl(img_num))//".xyz",dummy_nat)

              call services_read_xyz(trim(pub_base_rootname)//trim(adjustl(img_num))//".xyz",dummy_nat,posbuf,dummy_species)

              do row=0,mdl%nat-1
                mep%im_positions(row*3 + 1,i+2) = posbuf(1,row+1)
                mep%im_positions(row*3 + 2,i+2) = posbuf(2,row+1)
                mep%im_positions(row*3 + 3,i+2) = posbuf(3,row+1)
                !elements(row)%group_id = group_buf(row)
              end do
            end do

            deallocate(posbuf,dummy_species,stat=ierr)
            call utils_dealloc_check("neb_init","xyz read arrays",ierr)
         end if
      else if (interm_specified) then
         ! Intermediate configuration exists. Do linear interpolation R->I->P
         !call rundat_blocks_exec(mdl, 'INTERMEDIATE')
         call rundat_blocks_exec(mdl, 'INTERMEDIATE')

         if (pub_on_imroots_root) then
           call neb_extract_coords_from_mdl(mdl, interm_pos)
           ! Work out the R-I distance and T-P distance
           temp_pos(:) = (interm_pos - mep%im_positions(:,1))
           ri_dist = sqrt(dot_product(temp_pos,temp_pos))

           temp_pos(:) = (mep%im_positions(:,pub_num_images+2) - interm_pos)
           ip_dist = sqrt(dot_product(temp_pos,temp_pos))

           total_dist = ri_dist + ip_dist

           ! kkbd: The intermediate bead should be placed so beads on either side are
           !       as equally spaced as possible (it shouldnt naively be the middle
           !       image)
           !       Also make sure we always have a bead explicitly placed at the
           !       intermediate guess.
           interm_bead = nint(ri_dist * real(pub_num_images,dp) / total_dist) + 1
           interm_bead = max(interm_bead,2)
           interm_bead = min(interm_bead,pub_num_images+1)

           write(orig_stdout,*)"NEB: Guessed intermediate will be &
                          &initialized at bead number ",interm_bead-1

           do i = 2,pub_num_images+1
              if (i < interm_bead) then
                call neb_subtract(interm_pos, mep%im_positions(:,1),mdl, temp_pos)
                temp_pos = (temp_pos / real(interm_bead-1,dp))
                mep%im_positions(:,i) = mep%im_positions(:,1) &
                                       + (temp_pos(:) * real(i-1,dp))
                call neb_wrap(mep%im_positions(:,i), mdl)
              else if (i > interm_bead) then
                call neb_subtract(mep%im_positions(:,pub_num_images+2), interm_pos, mdl, temp_pos)
                temp_pos = (temp_pos / real(pub_num_images-interm_bead+2,dp))
                mep%im_positions(:,i) = interm_pos &
                                       + (temp_pos(:) * real(i-interm_bead,dp))
                call neb_wrap(mep%im_positions(:,i), mdl)
              else if (i == interm_bead) then
                mep%im_positions(:,i) = interm_pos
              end if
           end do
         end if
       else
         if (pub_on_imroots_root) then
           ! No intermediate structure specified. Linear interpolation R->P
           call neb_subtract(mep%im_positions(:,pub_num_images+2), &
                         mep%im_positions(:,1),mdl, temp_pos)
           temp_pos = (temp_pos / real(pub_num_images+1,dp))
           do i=2,pub_num_images+1
             mep%im_positions(:,i) = mep%im_positions(:,1) &
                                  + (temp_pos(:) * real(i-1,dp))
             call neb_wrap(mep%im_positions(:,i), mdl)
           end do
         end if
       end if

      if (pub_on_imroots_root) call esdf_close

    end if !imroots_root

    ! QC Test Output. Unfortunately as we're not using stdout we can't use the utils routine.
    if (pub_on_imroots_root .and. pub_print_qc) then
      outtag = '[reac_en]'
      write(orig_stdout,'(a30,f23.12)') '<QC> '//adjustr(outtag)//':', mep%reactant_energy
      outtag = '[prod_en]'
      write(orig_stdout,'(a30,f23.12)') '<QC> '//adjustr(outtag)//':', mep%product_energy
    end if

    call comms_barrier

    call neb_scatter_mep_info(mep, mdl)

    call neb_sync_path(mdl, mep, print_summary=.false.)

    if (pub_on_imroots_root) write(orig_stdout,'(a,a)')"NEB: Starting with ",trim(adjustl(mep%update_method))

    if (pub_write_xyz) call services_write_xyz(mdl%elements, pub_rootname, &
                                               "Initial NEB coordinates")

    deallocate(temp_pos, temp_force,stat=ierr)
    call utils_dealloc_check('neb_init','temp arrays',ierr)
    deallocate(reac_pos, prod_pos, interm_pos, stat=ierr)
    call utils_dealloc_check('neb_init', 'reactant, product, and intermediate positions',ierr)

  end subroutine neb_init

  subroutine neb_optimize(mdl, neb_mep)
    !====================================================================!
    ! Calculate energies and nudged forces of the beads on the NEB chain !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  energy  (  out) - Total energy of this NEB bead                   !
    !  forces  (  out) - Nudged force on this NEB bead                   !
    !  mdl     (inout) - mdl corresponding to this NEB bead              !
    !  neb_mep (inout) - path corresponding to the NEB chain             !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Aug 2018                                    !
    !====================================================================!

    use image_comms, only: pub_imroots_root_id, pub_my_image, pub_on_imroots_root, &
           orig_stdout
    use comms, only: pub_on_root, comms_reduce, comms_bcast, pub_image_comm, &
           pub_imroots_comm, pub_root_proc_id, pub_imroots_comm, comms_allgather
    use rundat, only: pub_num_images, pub_edft, pub_rootname, &
           pub_read_denskern, pub_read_tightbox_ngwfs, pub_read_hamiltonian, &
           pub_write_denskern, pub_write_tightbox_ngwfs, pub_write_hamiltonian, &
           tssearch_energy_tol, tssearch_force_tol, tssearch_disp_tol, &
           pub_edft_spin_fix, pub_edft_spin_fix_orig, pub_neb_converge_all, &
           pub_debug_on_root, pub_neb_max_iter, pub_neb_glbfgs_history_size, &
           pub_neb_print_summary, pub_neb_continuation
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort
    use services, only: services_write_xyz
    use forces, only: forces_apply_constraints

    implicit none

    ! Arguments
    type(MODEL),    intent(inout) :: mdl
    type(NEB_PATH), intent(inout) :: neb_mep

    ! Local Variables
    real(kind=dp), parameter :: neb_max_ionic_step = 0.5 !bohr

    integer :: i, j, ierr, iteration, old_ci
    character(len=80) :: xyz_header
    real(kind=dp), dimension(3) :: temp_dr, temp_f, temp_r
    real(kind=dp), allocatable, dimension(:)   :: v, f, r, new_r, dr, old_f, old_r, ci_v, temp_q, temp_z
    real(kind=dp), allocatable, dimension(:,:) :: force_2d
    real(kind=dp) :: old_e, energy, vf, dE, Fmax, Fcur, dRmax, dRcur
    real(kind=dp) :: dt = 0.2_dp
    logical :: conv

    integer :: my_start

    allocate(v(3*neb_mep%nat), r(3*neb_mep%nat), dr(3*neb_mep%nat), &
             new_r(3*neb_mep%nat), f(3*neb_mep%nat), old_r(3*neb_mep%nat), &
             old_f(3*neb_mep%nat), ci_v(3*neb_mep%nat), stat=ierr)
    call utils_alloc_check("neb_optimize","vecs",ierr)

    allocate(temp_q(3*neb_mep%nat*pub_num_images), &
             temp_z(3*neb_mep%nat*pub_num_images), stat=ierr)
    call utils_alloc_check("neb_optimize","temp arrays",ierr)

    allocate(force_2d(1:3,mdl%nat), stat=ierr)
    call utils_alloc_check("neb_optimize","2D force",ierr)

    v = 0.0_dp
    energy = 0.0_dp
    old_e = 0.0_dp
    old_r = 0.0_dp
    old_f = 0.0_dp

    old_ci = neb_mep%climbing_image

    pub_write_tightbox_ngwfs = .true.
    if (pub_edft) then
      pub_write_hamiltonian = .true.
    else
      pub_write_denskern = .true.
    end if

    iteration = 0

    if (pub_neb_continuation) then
      call neb_cont_read(mdl, neb_mep, v, iteration)

      ! Periodically reset NGWFs etc
      if ((iteration > 0) .and. mod(iteration, 6) == 0) then
        ! Reset NGWFs etc
        if(pub_on_imroots_root) write(orig_stdout,*) "NEB: Resetting NGWFs etc on NEB iteration",iteration
        pub_read_tightbox_ngwfs = .false.
        if (pub_edft) then
          pub_read_hamiltonian = .false.
        else
          pub_read_denskern = .false.
        end if
        if (pub_edft .and. (pub_edft_spin_fix_orig > 0)) then
          ! kkbd: Reset the fixed spin duration if we're doing a free-
          !       spin run with a fixed-spin start
          pub_edft_spin_fix = pub_edft_spin_fix_orig
          if (pub_on_root) write(orig_stdout,*) &
                  "NEB: Resetting: edft_spin_fix &
                  &parameter set to",pub_edft_spin_fix
        end if
      else
        pub_read_tightbox_ngwfs = .true.
        if (pub_edft) then
          pub_read_hamiltonian = .true.
        else
          pub_read_denskern = .true.
        end if
      end if
    end if

    if(pub_on_root) then
      xyz_header = "Initial NEB Coordinates"
      call services_write_xyz(mdl%elements, pub_rootname, trim(xyz_header))
    end if

    main_loop : do

      iteration = iteration + 1

      if ((pub_neb_max_iter > 0) .and. (iteration > pub_neb_max_iter)) then
        if (pub_on_imroots_root) write(orig_stdout,*) "Reached max NEB iterations, exiting"
        exit
      end if

      call neb_energy_and_force(energy, force_2d, mdl, neb_mep, is_cmplx=.false., en_sync=.true.)

      call forces_apply_constraints(force_2d, mdl)

      do i=1,neb_mep%nat
        f((i*3)-2) = force_2d(1,i)
        f((i*3)-1) = force_2d(2,i)
        f((i*3))   = force_2d(3,i)
      end do

      neb_mep%my_force = f

      call neb_extract_coords_from_mdl(mdl, r)

      if ((neb_mep%update_method == "GLBFGS") .and. (neb_mep%ci_delay == 0)) then
        if (pub_on_imroots_root) write(orig_stdout,*) "NEB: Climbing image enabled on a GL-BFGS run. This is not stable."
        if (pub_on_imroots_root) write(orig_stdout,*) "NEB: Automatically switching to FIRE."
        neb_mep%update_method = "FIRE"
      end if

      ! Now determine a new step
      ! These optimization methods were adapted from "Optimization Methods for finding minimum energy paths"
      ! https://doi.org/10.1063/1.2841941
      if (neb_mep%update_method == "QUICKMIN") then
        ! kkbd: This isn't worth ever actually using. It's just for testing.
        vf = dot_product(v,f)

        if (vf < 0.0_dp) then
          v = 0.0_dp
        end if

        v = v + dt * f
        dr = dt * v
      else if (neb_mep%update_method == "FIRE") then
        vf = dot_product(v,f)
        v = ((1.0_dp - neb_mep%fire_alpha) * v) &
            + (neb_mep%fire_alpha * f * sqrt(dot_product(v,v))/sqrt(dot_product(f,f)))

        if (vf < 0.0_dp) then
          v = 0.0_dp
          neb_mep%cut = iteration
          neb_mep%fire_dt = neb_mep%fire_dt * neb_mep%fire_dec
          neb_mep%fire_alpha = neb_mep%orig_alpha
        else if ((vf >= 0.0_dp) .and. (iteration - neb_mep%cut > neb_mep%n_min)) then
          dt = min(dt * neb_mep%fire_inc, neb_mep%dt_max)
          neb_mep%fire_alpha = neb_mep%fire_alpha * neb_mep%fire_adec
        end if

        v = v + dt * f

        dr = dt * v
      else if (neb_mep%update_method == "GLBFGS") then
        ! Gather to imroots_root for optimization there
        if (pub_on_root) then
          call comms_allgather(temp_q,f,comm=pub_imroots_comm,gather_not_allgather=.true.)
          neb_mep%lbfgs_q = temp_q
        end if

        if (iteration == 1) then
          ! Do a simple MD step
          dr = f * dt * dt
          neb_mep%old_q = -neb_mep%lbfgs_q
          if (pub_on_root) then
            call comms_allgather(temp_z,dr, &
                 comm=pub_imroots_comm,gather_not_allgather=.true.)
            neb_mep%old_z = temp_z
          end if

        else ! not iteration 1
          if (pub_on_imroots_root) then
            ! Do the LBFGS step
            neb_mep%lbfgs_q = -neb_mep%lbfgs_q

            j = modulo(iteration-1, pub_neb_glbfgs_history_size) + 1

            temp_q = neb_mep%lbfgs_q
            temp_z = neb_mep%old_q
            neb_mep%lbfgs_y(:,j) =  temp_q - temp_z!cur_y
            neb_mep%lbfgs_s(:,j) = neb_mep%old_z !cur_s
            neb_mep%lbfgs_rho(j) = 1.0_dp / dot_product(neb_mep%lbfgs_y(:,j),neb_mep%lbfgs_s(:,j))

            neb_mep%old_q = neb_mep%lbfgs_q

            do i = iteration - 1, iteration - pub_neb_glbfgs_history_size, -1
              j = modulo(i,pub_neb_glbfgs_history_size) + 1
              neb_mep%lbfgs_alpha(j) = neb_mep%lbfgs_rho(j) &
                                       * dot_product(neb_mep%lbfgs_s(:,j),neb_mep%lbfgs_q)
              neb_mep%lbfgs_q = neb_mep%lbfgs_q - neb_mep%lbfgs_y(:,j) * neb_mep%lbfgs_alpha(j)
            end do

            j = modulo(iteration-1, pub_neb_glbfgs_history_size) + 1
            neb_mep%lbfgs_z = neb_mep%lbfgs_q &
                              * (dot_product(neb_mep%lbfgs_y(:,j),neb_mep%lbfgs_s(:,j))&
                              / dot_product(neb_mep%lbfgs_y(:,j),neb_mep%lbfgs_y(:,j)))

            do i = iteration - pub_neb_glbfgs_history_size, iteration - 1
              j = modulo(i,pub_neb_glbfgs_history_size) + 1
              neb_mep%lbfgs_z = neb_mep%lbfgs_z + neb_mep%lbfgs_s(:,j) &
                                * (neb_mep%lbfgs_alpha(j) - neb_mep%lbfgs_rho(j) &
                                  * dot_product(neb_mep%lbfgs_y(:,j),neb_mep%lbfgs_z))
            end do

            neb_mep%old_z = neb_mep%lbfgs_z

          end if ! pub_on_imroots_root

          ! Distibute results
          if (pub_on_root) then
            write(stdout,*)"lbfgs_z",neb_mep%lbfgs_z
            temp_z = neb_mep%lbfgs_z
            call comms_bcast(pub_imroots_root_id, temp_z, comm=pub_imroots_comm)

            ! separate out my part of the array
            my_start = (3*neb_mep%nat*pub_my_image) + 1
            dr = -temp_z(my_start:my_start + 3*neb_mep%nat - 1)
          end if

          call comms_bcast(pub_root_proc_id, dr, comm=pub_image_comm)
        end if ! Not first iteration

      else
        call utils_abort("NEB: Bad Update Method specified: "//trim(neb_mep%update_method))
      end if ! mep%update_method

      if (pub_debug_on_root) write(stdout,*) "NEB: dr ",dr

      old_f = f
      old_r = r
      old_ci = neb_mep%climbing_image

      ! Limit the size of dr and find dRmax for position convergence
      dRmax = 0.0_dp
      do i=1,neb_mep%nat
        temp_dr(1) = dr((i*3)-2)
        temp_dr(2) = dr((i*3)-1)
        temp_dr(3) = dr((i*3)-0)
        dRcur = sqrt(dot_product(temp_dr,temp_dr))
        if (dRcur > dRmax) dRmax = dRcur
      end do

      if (dRmax > neb_max_ionic_step) then
        ! Scale all ionic steps down so max is neb_max_ionic_step
        do i=1,neb_mep%nat
          temp_dr(1) = dr((i*3)-2)
          temp_dr(2) = dr((i*3)-1)
          temp_dr(3) = dr((i*3)-0)

          temp_dr = temp_dr * neb_max_ionic_step / dRmax

          dr((i*3)-2) = temp_dr(1)
          dr((i*3)-1) = temp_dr(2)
          dr((i*3)-0) = temp_dr(3)
        end do

        dRmax = neb_max_ionic_step
      end if

      ! I think this can be cleaned up similarly to the new geomopt routines.
      if (pub_on_root) then
        call comms_allgather(temp_z,dr,comm=pub_imroots_comm,gather_not_allgather=.true.)
        neb_mep%old_z = temp_z
      end if

      ! Force and energy convergence
      Fmax = 0.0_dp
      do i=1,neb_mep%nat
        temp_f(1) = f((i*3)-2)
        temp_f(2) = f((i*3)-1)
        temp_f(3) = f((i*3)-0)

        Fcur = sqrt(dot_product(temp_f,temp_f))
        if (Fcur > Fmax) Fmax = Fcur
      end do

      dE = abs(energy - old_e)
      old_e = energy

      conv = .false.
      ! Do full convergence if we're doing so for all beads, OR we're the climbing image.
      if (pub_neb_converge_all .or. (neb_mep%climbing_image == pub_my_image)) then
        if ((dRmax < tssearch_disp_tol) &
          .and. (Fmax < tssearch_force_tol) &
          .and. (dE < tssearch_energy_tol)) conv = .true.
      else
        ! Just converge forces. NOTE: Eventually we could also converge projected displacement
        if (Fmax < tssearch_force_tol) conv = .true.
      end if

      if (pub_on_root) then
        call comms_reduce('AND',conv,pub_imroots_comm)
      end if
      call comms_bcast(pub_root_proc_id,conv,1,pub_image_comm)

      new_r = r + dr
      do i=1, neb_mep%nat
        ! Wrap coords back into mdl
        temp_r(1) = new_r((i*3)-2)
        temp_r(2) = new_r((i*3)-1)
        temp_r(3) = new_r((i*3)-0)

        call neb_wrap(temp_r, mdl)

        new_r((i*3)-2) = temp_r(1)
        new_r((i*3)-1) = temp_r(2)
        new_r((i*3)-0) = temp_r(3)
      end do

      call neb_load_coords_into_mdl(mdl, new_r)

      call neb_sync_path(mdl, neb_mep)

      if(pub_on_root) then
        write(xyz_header,'(a, i5)') "NEB Coordinates After Iteration",iteration
        call services_write_xyz(mdl%elements, pub_rootname, trim(xyz_header))
      end if

      if (pub_on_imroots_root) then
        ! Write out path XYZ
        call neb_path_xyz(mdl, neb_mep)
      end if

      if (pub_neb_print_summary) call neb_print_conv(conv,iteration, dE, Fmax, dRmax, neb_mep)

      if (conv) then
        if (neb_mep%ci_delay > 0) then
          if (pub_on_imroots_root) write(orig_stdout,*) "NEB: Path converged but the climbing image wasn't switched on."
          if (pub_on_imroots_root) write(orig_stdout,*) "NEB: Skipping convergence and continuing with climbing image."
          if (pub_on_imroots_root) write(orig_stdout,*) "NEB: Enabling climbing image at iteration: ",iteration
          neb_mep%ci_delay = 0
        else
          if (pub_on_imroots_root) write(orig_stdout,*) "NEB: Converged! Exiting."
          call neb_print_path(neb_mep)
          exit
        end if
      end if

      ! Decrease CI delay counter
      if (neb_mep%ci_delay > 0) then
        neb_mep%ci_delay = neb_mep%ci_delay - 1
        if (neb_mep%ci_delay == 0) then
          if (pub_on_imroots_root) write(orig_stdout,*)"NEB: Enabling climbing image at iteration: ",iteration
        else
          if (pub_on_imroots_root) write(orig_stdout,*)"NEB: Reducing climbing image delay to: ",neb_mep%ci_delay
        end if
      end if

      ! Write out continuation
      if (pub_on_root) then
        call neb_cont_write(mdl, neb_mep, v, iteration)
      end if

      ! Periodically reset NGWFs etc
      if (mod(iteration, 6) == 0) then
        ! Reset NGWFs etc
        if(pub_on_imroots_root) write(orig_stdout,*) "NEB: Resetting NGWFs etc on NEB iteration",iteration
        pub_read_tightbox_ngwfs = .false.
        if (pub_edft) then
          pub_read_hamiltonian = .false.
        else
          pub_read_denskern = .false.
        end if
        if (pub_edft .and. (pub_edft_spin_fix_orig > 0)) then
          ! kkbd: Reset the fixed spin duration if we're doing a free-
          !       spin run with a fixed-spin start
          pub_edft_spin_fix = pub_edft_spin_fix_orig
          if (pub_on_imroots_root) write(orig_stdout,*) &
                  "NEB: Resetting: edft_spin_fix &
                  &parameter set to",pub_edft_spin_fix
        end if
      else
        pub_read_tightbox_ngwfs = .true.
        if (pub_edft) then
          pub_read_hamiltonian = .true.
        else
          pub_read_denskern = .true.
        end if
      end if

    end do main_loop

    deallocate(v, r, dr, new_r, f, old_r, old_f, ci_v, stat=ierr)
    call utils_dealloc_check("neb_optimize","vecs",ierr)

    deallocate(force_2d, stat=ierr)
    call utils_dealloc_check("neb_optimize","2D force",ierr)

  end subroutine neb_optimize

  subroutine neb_energy_and_force(energy, forces, mdl, neb_mep, is_cmplx, en_sync)
    !====================================================================!
    ! Calculate energies and nudged forces of the beads on the NEB chain !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  energy  (  out) - Total energy of this NEB bead                   !
    !  forces  (  out) - Nudged force on this NEB bead                   !
    !  mdl     (inout) - mdl corresponding to this NEB bead              !
    !  neb_mep (inout) - path corresponding to the NEB chain             !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Jan 2018                                    !
    !====================================================================!

    use comms, only: pub_on_root
    use energy_and_force, only: energy_and_force_calculate
    use rundat, only: pub_debug_on_root
    use services, only: services_flush
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    ! Arguments
    real(kind=dp),               intent(  out) :: energy
    real(kind=dp),dimension(:,:),intent(  out) :: forces
    type(MODEL),                 intent(inout) :: mdl
    type(NEB_PATH),              intent(inout) :: neb_mep
    logical, optional,           intent(in   ) :: is_cmplx
    logical, optional,           intent(in   ) :: en_sync

    ! Local Variables
    integer :: ierr
    real(kind=dp),allocatable,dimension(:,:) :: raw_forces

    allocate(raw_forces(1:3,mdl%nat),stat=ierr)
    call utils_alloc_check('neb_energy_and_force','raw_forces',ierr)

    ! First calculate actual total energy and raw forces
    ! agrecocmplx
    call energy_and_force_calculate(energy, raw_forces, mdl,  &
         is_cmplx=is_cmplx)

    neb_mep%my_energy = energy

    if (present(en_sync)) then
      if (pub_debug_on_root) write(stdout,*)"NEB: Syncing during e&f calc"
      if (en_sync) call neb_sync_path(mdl,neb_mep, print_summary=.false.)
    end if

    call neb_nudge_forces(neb_mep, energy, raw_forces, forces)

    if (pub_debug_on_root) then
      write(stdout,*)"NEB: Done nudging forces"
      call services_flush
    end if

    deallocate(raw_forces,stat=ierr)
    call utils_dealloc_check('neb_energy_and_force','raw_forces',ierr)

  end subroutine neb_energy_and_force

  subroutine neb_sync_path(mdl, neb_mep, print_summary)
    !====================================================================!
    ! Distribute bead energy and position information to neighbors       !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Jan 2018                                    !
    !====================================================================!

    use comms, only: pub_on_root
    use energy_and_force, only: energy_and_force_calculate
    use rundat, only: pub_neb_print_summary
    use services, only: services_flush
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    ! Arguments
    type(MODEL),                 intent(inout) :: mdl
    type(NEB_PATH),              intent(inout) :: neb_mep
    logical, optional,           intent(in   ) :: print_summary

    ! Local Variables
    integer :: ierr
    real(kind=dp),allocatable,dimension(:)   :: pos, forces_1d
    real(kind=dp),allocatable,dimension(:,:) :: forces
    real(kind=dp) :: energy
    logical :: loc_print_summary

    allocate(pos(1:3 * mdl%nat),stat=ierr)
    call utils_alloc_check('neb_energy_and_force','positions',ierr)

    allocate(forces(3,mdl%nat),stat=ierr)
    call utils_alloc_check('neb_energy_and_force','raw_forces',ierr)

    allocate(forces_1d(1:3 * mdl%nat),stat=ierr)
    call utils_alloc_check('neb_energy_and_force','forces_1d',ierr)

    loc_print_summary = pub_neb_print_summary

    if(present(print_summary)) loc_print_summary = print_summary

    energy = neb_mep%my_energy

    forces_1d = neb_mep%my_force

    ! Make sure stdout is up-to-date while we're waiting
    call services_flush

    ! Get positions
    call neb_extract_coords_from_mdl(mdl, pos)

    ! Now sync NEB path info
    call neb_gather_mep_info(energy, forces_1d, pos, neb_mep)
    call neb_scatter_mep_info(neb_mep, mdl)

    if (loc_print_summary) call neb_print_path(neb_mep)

    deallocate(pos,stat=ierr)
    call utils_dealloc_check('neb_sync_path','helper arrays',ierr)

    deallocate(forces,stat=ierr)
    call utils_dealloc_check('neb_sync_path','raw_forces',ierr)

    deallocate(forces_1d,stat=ierr)
    call utils_dealloc_check('neb_sync_path','forces_1d',ierr)

  end subroutine neb_sync_path

  subroutine neb_tangent(mep, image_energy, tangent)
    !====================================================================!
    ! Approximates the unit tangent to the path at our point. There are  !
    ! several ways to go about doing this. For stability it is reported  !
    ! that an upwinding tangent is favourable. This tangent points to    !
    ! the higher energy neighbor, or in the direction of a linear        !
    ! interpolation of the two neighbors at extrema.                     !
    ! See G. Henkelman and H. J\'onsson, J. Chem. Phys. 113, 9978 (2000) !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Jan 2018                                    !
    !====================================================================!

    use comms, only: pub_on_root
    use image_comms, only: pub_my_image
    use rundat, only: pub_num_images, pub_debug_on_root, pub_devel_code
    use services, only: services_flush
    use utils, only: utils_devel_code
    implicit none

    ! Arguments
    type(NEB_PATH),            intent(in   ) :: mep
    real(kind=dp),             intent(in   ) :: image_energy
    real(kind=dp),dimension(:),intent(  out) :: tangent

    ! Local Variables
    real(kind=dp) :: prev_e, next_e, v_max, v_min
    logical :: loc_write_debug_info

    loc_write_debug_info = utils_devel_code(.false.,"NEB","DEBUG_OUT", pub_devel_code)
    loc_write_debug_info = pub_debug_on_root .or. (pub_on_root.and.loc_write_debug_info)

    prev_e = mep%prev_energy
    next_e = mep%next_energy

    if (pub_my_image == 0) then
      ! Try to read reactant energy from devel_code
      prev_e = mep%reactant_energy
      !prev_e = image_energy - abs(next_e - image_energy)
    else if (pub_my_image == pub_num_images-1) then
      next_e = mep%product_energy
      ! Try to read product energy from devel_code
      !next_e = image_energy - abs(prev_e - image_energy)
    end if

    if ((prev_e > image_energy) .and. (image_energy > next_e)) then
      tangent = mep%prev_to_cur
    else if ((prev_e < image_energy) .and. (image_energy < next_e)) then
      tangent = mep%cur_to_next
    else ! We're at an extreme
      v_max = max(abs(next_e-image_energy),abs(prev_e-image_energy))
      v_min = min(abs(next_e-image_energy),abs(prev_e-image_energy))
      if (next_e > prev_e) then
        tangent = v_max * (mep%cur_to_next) + v_min * (mep%prev_to_cur)
      else
        tangent = v_min * (mep%cur_to_next) + v_max * (mep%prev_to_cur)
      end if
    end if

    if (loc_write_debug_info) then
      write(stdout,*) "NEB: --NEB TANGENT INFO--"
      write(stdout,*) "NEB: tangent",tangent
      write(stdout,*) "NEB: prev_e", prev_e
      write(stdout,*) "NEB: next_e", next_e
      write(stdout,*) "NEB: prev_to_cur", mep%prev_to_cur
      write(stdout,*) "NEB: cur_to_next", mep%cur_to_next
      write(stdout,*) "NEB: |tangent|",sqrt(dot_product(tangent,tangent))
    end if

    tangent = tangent / sqrt(dot_product(tangent,tangent))

    if (loc_write_debug_info) then
      write(stdout,*) "NEB: unit tangent",tangent
      call services_flush
    end if

  end subroutine neb_tangent

  subroutine neb_nudge_forces(mep, energy, force, nudged_force)
    !====================================================================!
    ! Calculated the nudged force on this NEB bead.                      !
    ! For most beads this is the sum of:                                 !
    !  Real force with the component parallel to the path projected out  !
    !  Spring force with the component perpendicular to the path         !
    !    projected out                                                   !
    ! The climbing image acts differently. It feels no spring force and  !
    ! instead moves in a transformed potential surface such that the     !
    ! minimum lies at the transition state.                              !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Jan 2018                                    !
    !====================================================================!

    use comms, only: pub_on_root
    use geometry, only: geometry_magnitude
    use image_comms, only: pub_my_image
    use model_type, only: MODEL
    use rundat, only: pub_debug_on_root, pub_devel_code
    use services, only: services_flush
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_devel_code
    implicit none

    ! Arguments
    type(NEB_PATH),              intent(in   ) :: mep
    real(kind=dp),               intent(in   ) :: energy
    real(kind=dp),dimension(:,:),intent(in   ) :: force
    real(kind=dp),dimension(:,:),intent(  out) :: nudged_force

    ! Local Variables
    integer :: i
    integer :: ierr
    real(kind=dp),allocatable,dimension(:) :: tangent
    real(kind=dp),allocatable,dimension(:) :: force_1d, nudged_real, &
                                              nudged_spring, full_spring
    logical :: loc_write_debug_info

    allocate(tangent(1:3*mep%nat), stat=ierr)
    call utils_alloc_check('neb_nudge_forces', 'path tangent', ierr)

    allocate(force_1d(1:3*mep%nat), stat=ierr)
    call utils_alloc_check('neb_nudge_forces', 'force_1d', ierr)

    allocate(full_spring(1:3*mep%nat), stat=ierr)
    call utils_alloc_check('neb_full_forces', 'full force components', ierr)

    allocate(nudged_real(1:3*mep%nat), nudged_spring(1:3*mep%nat), stat=ierr)
    call utils_alloc_check('neb_nudge_forces', 'nudged force components', ierr)

    ! Write out debug info if we're in debug mode or NEB debug mode
    loc_write_debug_info = utils_devel_code(.false.,"NEB","DEBUG_OUT", pub_devel_code)
    loc_write_debug_info = pub_debug_on_root .or. (pub_on_root .and. loc_write_debug_info)

    do i=1,mep%nat
      force_1d((i*3)-2) = force(1,i)
      force_1d((i*3)-1) = force(2,i)
      force_1d((i*3))   = force(3,i)
    end do

    call neb_tangent(mep, energy, tangent)

    ! Calculate nudged image force
    !   = nudged real force
    !     + nudged spring force
    if (pub_my_image /= mep%climbing_image) then

      if (loc_write_debug_info) then
        write(stdout,*)"NEB: --NEB FORCE INFO--"
        write(stdout,*)"NEB: full real force",force_1d
        call services_flush
      end if

      ! Kept separate for debug outputs
      nudged_real = force_1d - dot_product(force_1d,tangent) * tangent
      if (loc_write_debug_info) then
        write(stdout,*)"NEB: nudged real force",nudged_real
        call services_flush
      end if

      ! Rewrote for improved tangent
      full_spring = (mep%spring_k * mep%cur_to_next)
      full_spring(:) = full_spring(:) - (mep%spring_k * mep%prev_to_cur(:))
      !do i = 1, 3*mep%nat
      !  full_spring(i) = full_spring(i) - (spring_k * mep%prev_to_cur(i))
      !end do

      if (loc_write_debug_info) then
        write(stdout,*)"NEB: full spring force",full_spring
        call services_flush
      end if

      !nudged_spring = dot_product(full_spring,tangent) * tangent
      ! kkbd: This is a superior nudged spring force, see DOI: 10.1063/1.1323224
      nudged_spring = mep%spring_k * (sqrt(dot_product(mep%cur_to_next,mep%cur_to_next)) - &
                                      sqrt(dot_product(mep%prev_to_cur,mep%prev_to_cur))) * tangent

      if (loc_write_debug_info) then
        write(stdout,*)"NEB: nudged spring force",nudged_spring
        call services_flush
      end if

      force_1d  = nudged_real + nudged_spring

    else ! If climbing image we climb up the path and feel no spring force
      if (pub_on_root) then
        write(stdout,*)"NEB: This image is our climbing image"
      end if

      if (loc_write_debug_info) then
        write(stdout,*)"NEB: --NEB FORCE INFO--"
        write(stdout,*)"NEB: full real force",force_1d
        call services_flush
      end if

      force_1d = force_1d - (2.0_dp * dot_product(force_1d,tangent) * tangent)
    end if

    if (loc_write_debug_info) then
      write(stdout,*)"NEB: Nudged force",force_1d
      call services_flush
    end if

    do i=1,mep%nat
      nudged_force(1,i) = force_1d((i*3)-2)
      nudged_force(2,i) = force_1d((i*3)-1)
      nudged_force(3,i) = force_1d((i*3))
    end do

    deallocate(tangent, stat=ierr)
    call utils_dealloc_check('neb_nudge_forces', 'path tangent', ierr)

    deallocate(force_1d, stat=ierr)
    call utils_dealloc_check('neb_nudge_forces', 'force_1d', ierr)

    deallocate(full_spring,stat=ierr)
    call utils_dealloc_check('neb_nudge_forces','force components',ierr)

    deallocate(nudged_real, nudged_spring, stat=ierr)
    call utils_dealloc_check('neb_nudge_forces', 'nudged force components', ierr)

  end subroutine neb_nudge_forces

!-------------------------------------------------------------------------------
! Convergence Checking Routine
!-------------------------------------------------------------------------------

  logical function neb_converged(image_conv, iteration, dE, Fmax, dRmax, neb_mep)
    !====================================================================!
    ! Currently the path is converged if each image is converged         !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Jan 2018                                    !
    !====================================================================!

    use comms, only: comms_barrier, pub_image_comm,pub_on_root, &
         pub_imroots_comm, pub_root_proc_id, &
         comms_reduce, comms_bcast, comms_send, comms_recv
    use image_comms, only: pub_on_imroots_root, &
         pub_imroots_root_id, pub_my_rank_in_imroots
    use rundat, only: pub_num_images, pub_neb_print_summary
    use utils, only: utils_flush
    implicit none

    ! Arguments
    logical,       intent(in   ) :: image_conv
    real(kind=dp), intent(in   ) :: dE, Fmax, dRmax
    integer,       intent(in   ) :: iteration
    type(NEB_PATH),intent(in   ) :: neb_mep

    ! Local Variables
    integer :: i, im, ierr, string_index
    logical :: convd

    neb_converged = image_conv

    if (pub_on_root) then
      call comms_reduce('AND',neb_converged,pub_imroots_comm)
    end if

    call comms_bcast(pub_root_proc_id,neb_converged,1,pub_image_comm)

    if (pub_neb_print_summary) then
      call neb_print_conv(neb_converged, iteration, dE, Fmax, dRmax, neb_mep)
    end if

    return

  end function neb_converged

!-------------------------------------------------------------------------------
! MEP Distribution Routines
!-------------------------------------------------------------------------------

  subroutine neb_gather_mep_info(im_energy, im_forces, im_pos, mep)
    !====================================================================!
    ! Gather together information on the positions and energies of each  !
    ! bead in the NEB chain, so they can be scattered back to neighbors  !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Jan 2018                                    !
    !====================================================================!

    use comms, only: pub_imroots_comm, pub_image_comm, comms_barrier, &
                     comms_send, comms_recv, pub_on_root
    use image_comms, only: pub_my_rank_in_imroots, pub_on_imroots_root, &
                           pub_imroots_root_id
    use rundat, only: pub_num_images, pub_debug_on_root
    use services, only: services_flush
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    ! Arguments
    real(kind=DP),              intent(in   ) :: im_energy
    real(kind=DP), dimension(:),intent(in   ) :: im_forces
    real(kind=DP), dimension(:),intent(in   ) :: im_pos
    type(NEB_PATH),             intent(inout) :: mep

    ! Local Variables
    integer :: i, ierr, temp_climb(1)
    real(kind=dp) :: temp_energy
    real(kind=dp),allocatable,dimension(:) :: temp_forces
    real(kind=dp),allocatable,dimension(:) :: temp_pos

    allocate(temp_forces(size(im_forces)),temp_pos(size(im_pos)), stat=ierr)
    call utils_alloc_check("neb_gather_mep_info","temp arrays",ierr)

    if (pub_on_root) call comms_barrier(pub_imroots_comm)
    call comms_barrier(pub_image_comm)
    if (pub_on_imroots_root) then
      mep%im_forces(:,pub_my_rank_in_imroots+1) = im_forces
      mep%im_energies(pub_my_rank_in_imroots+1) = im_energy
      mep%im_positions(:,pub_my_rank_in_imroots+2) = im_pos
      do i=0,pub_num_images-1
        if (i /= pub_my_rank_in_imroots) then
          call comms_recv(i, temp_forces, tag=i, &
            length=size(mep%im_forces(:,i)), comm=pub_imroots_comm) !recv forces
          call comms_recv(i, temp_pos, tag=i+10*pub_num_images+2, &
            length=size(mep%im_positions(:,i)), comm=pub_imroots_comm) !recv positions
          call comms_recv(i, temp_energy, tag=20*pub_num_images+1+i, &
            comm=pub_imroots_comm) !recv energies
          mep%im_energies(i+1) = temp_energy
          mep%im_forces(:,i+1) = temp_forces
          mep%im_positions(:,i+2) = temp_pos
        end if
      end do
    else if (pub_on_root) then ! .not. pub_on_imroots_root
      call comms_send(pub_imroots_root_id, im_forces, length=size(im_forces), &
        & tag=pub_my_rank_in_imroots, comm=pub_imroots_comm) !send forces

      call comms_send(pub_imroots_root_id, im_pos, length=size(im_pos), &
        & tag=pub_my_rank_in_imroots+10*pub_num_images+2, comm=pub_imroots_comm) !send pos

      call comms_send(pub_imroots_root_id, im_energy, &
        & tag=20*pub_num_images+1+pub_my_rank_in_imroots, comm=pub_imroots_comm) !send energies

    end if ! pub_on_imroots_root

    if (pub_on_imroots_root .and. pub_debug_on_root) then
      write(stdout,*)"NEB: mep energies",mep%im_energies
      write(stdout,*)"NEB: forces",mep%im_forces
      call services_flush
    end if

    if (pub_on_imroots_root .and. (mep%ci_delay == 0)) then
      ! Slightly roundabout way because maxloc returns an array of size 1....
      temp_climb = maxloc(mep%im_energies) - 1
      mep%climbing_image = temp_climb(1)
      if (pub_debug_on_root) then
        write(stdout,*)"NEB: New climbing image:",mep%climbing_image
      end if
    end if

    if (pub_debug_on_root) then
      write(stdout,*)"NEB: my energy",im_energy
      write(stdout,*)"NEB: my force",im_forces
    end if

    call comms_barrier(pub_image_comm)

    deallocate(temp_forces,temp_pos, stat=ierr)
    call utils_dealloc_check("neb_gather_mep_info","temp arrays",ierr)

  end subroutine neb_gather_mep_info

  subroutine neb_scatter_mep_info(mep, mdl)
    !====================================================================!
    ! Scatter bead energies and positions to neighboring beads           !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Jan 2018                                    !
    !====================================================================!

    use comms, only: pub_imroots_comm, pub_image_comm, comms_barrier, &
                     comms_send, comms_recv, pub_on_root, pub_root_proc_id, &
                     comms_bcast
    use image_comms, only: pub_my_rank_in_imroots, pub_on_imroots_root, &
                           pub_imroots_root_id, pub_my_image
    use rundat, only: pub_num_images
    use services, only: services_flush
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    ! Arguments
    type(NEB_PATH),intent(inout) :: mep
    type(MODEL),   intent(inout) :: mdl

    ! Local Variables
    integer :: i, ierr, climb_im
    real(kind=dp) :: e_prev, e_next
    real(kind=dp), allocatable, dimension(:) :: im_cur, im_prev, im_next

    allocate(im_cur(1:3*mdl%nat), im_prev(1:3*mdl%nat), &
             im_next(1:3*mdl%nat), stat=ierr)
    call utils_alloc_check("neb_scatter_mep_info","image positions",ierr)

    if (pub_on_imroots_root) then
      do i=0, pub_num_images-1
        if (i /= pub_my_rank_in_imroots)then
          call comms_send(i,mep%im_positions(:,i+2), &
            length=size(mep%im_positions(:,i)), tag=i, comm=pub_imroots_comm)
        else
          im_cur = mep%im_positions(:,i+2)
        end if
      end do
    else if (pub_on_root) then ! Root but not imroots root
      call comms_recv(pub_imroots_root_id, im_cur, length=size(im_cur), &
        & tag=pub_my_rank_in_imroots, comm=pub_imroots_comm)
    end if

    call comms_barrier(pub_image_comm)
    call comms_bcast(pub_root_proc_id,im_cur,length=3*mdl%nat,comm=pub_image_comm)

    ! Update mdl
    call neb_load_coords_into_mdl(mdl,im_cur)

    ! Now send prev neighbor information
    if (pub_on_imroots_root) then
      do i=0,pub_num_images-1
        if (i /= pub_my_rank_in_imroots)then
          call comms_send(i,mep%im_positions(:,i+1), &
            length=size(mep%im_positions(:,i+1)), tag=i, comm=pub_imroots_comm)
        else
          im_prev = mep%im_positions(:,i+1)
        end if
      end do
    else if (pub_on_root) then
      call comms_recv(pub_imroots_root_id, im_prev, length=size(im_prev), &
        & tag=pub_my_rank_in_imroots, comm=pub_imroots_comm)
    end if

    call comms_barrier(pub_image_comm)
    call comms_bcast(pub_root_proc_id,im_prev,length=3*mdl%nat,comm=pub_image_comm)

    ! Now send next neighbor information
    if (pub_on_imroots_root) then
      do i=0,pub_num_images-1
        if (i /= pub_my_rank_in_imroots)then
          call comms_send(i,mep%im_positions(:,i+3), &
            length=size(mep%im_positions(:,i+3)), tag=i, comm=pub_imroots_comm)
        else
          im_next = mep%im_positions(:,i+3)
        end if
      end do
    else if (pub_on_root) then
      call comms_recv(pub_imroots_root_id, im_next, length=size(im_next), &
        & tag=pub_my_rank_in_imroots, comm=pub_imroots_comm)
    end if

    call comms_barrier(pub_image_comm)
    call comms_bcast(pub_root_proc_id,im_next,length=3*mdl%nat,comm=pub_image_comm)

    !Now send previous energy
    if (pub_on_imroots_root) then
      do i=1,pub_num_images-1
        if (i /= pub_my_rank_in_imroots)then
          call comms_send(i,mep%im_energies(i), tag=i, comm=pub_imroots_comm)
        else
          e_prev = mep%im_energies(i)
        end if
      end do
    else if (pub_on_root) then
      call comms_recv(pub_imroots_root_id, e_prev, tag=pub_my_rank_in_imroots, &
        comm=pub_imroots_comm)
    end if

    call comms_barrier(pub_image_comm)
    call comms_bcast(pub_root_proc_id,e_prev,comm=pub_image_comm)

    !Now send next energy
    if (pub_on_imroots_root) then
      do i=0,pub_num_images-2
        if (i /= pub_my_rank_in_imroots) then
          call comms_send(i,mep%im_energies(i+2), tag=i, comm=pub_imroots_comm)
        else
          e_next = mep%im_energies(i+2)
        end if
      end do
    else if (pub_on_root .and. (pub_my_rank_in_imroots < pub_num_images-1)) then
      call comms_recv(pub_imroots_root_id, e_next, tag=pub_my_rank_in_imroots, &
        comm=pub_imroots_comm)
    end if

    call comms_barrier(pub_image_comm)
    call comms_bcast(pub_root_proc_id,e_next,comm=pub_image_comm)

    ! Send climbing image info
    if (pub_on_imroots_root) then
      do i=0,pub_num_images-1
        if (i /= pub_my_rank_in_imroots) then
          call comms_send(i,mep%climbing_image, tag=i, comm=pub_imroots_comm)
        else
          climb_im = mep%climbing_image
        end if
      end do
    else if (pub_on_root) then
      call comms_recv(pub_imroots_root_id, climb_im, tag=pub_my_rank_in_imroots,&
        comm=pub_imroots_comm)
    end if

    call comms_barrier(pub_image_comm)
    call comms_bcast(pub_root_proc_id,climb_im,comm=pub_image_comm)

!    if (pub_on_root) then
    mep%prev_pos = im_prev
    call neb_subtract(im_cur, im_prev, mdl, mep%prev_to_cur)
    mep%next_pos = im_next
    call neb_subtract(im_next, im_cur, mdl, mep%cur_to_next)
    mep%climbing_image = climb_im
    if (pub_my_image > 0)mep%prev_energy = e_prev
    if (pub_my_image < pub_num_images-1)mep%next_energy = e_next
!    end if

    deallocate(im_cur,im_prev,im_next,stat=ierr)
    call utils_dealloc_check("neb_scatter_mep_info","image positions",ierr)

    call comms_barrier(pub_image_comm)

  end subroutine neb_scatter_mep_info

!---------------------------------------------------------------------------------------
! Model-to-Coords-to-Model
!---------------------------------------------------------------------------------------

  subroutine neb_extract_coords_from_mdl(mdl, positions)
    implicit none

    ! Arguments
    type(MODEL),               intent(in   ) :: mdl
    real(kind=dp),dimension(:),intent(  out) :: positions

    ! Local Variables
    integer :: i

    do i=1,mdl%nat
      positions(3*(i-1)+1) = mdl%elements(i)%centre%x
      positions(3*(i-1)+2) = mdl%elements(i)%centre%y
      positions(3*(i-1)+3) = mdl%elements(i)%centre%z
    end do

  end subroutine neb_extract_coords_from_mdl

  subroutine neb_load_coords_into_mdl(mdl, positions)
    implicit none

    ! Arguments
    type(MODEL),               intent(inout) :: mdl
    real(kind=dp),dimension(:),intent(in   ) :: positions

    ! Local Variables
    integer :: i

    do i=1,mdl%nat
      mdl%elements(i)%centre%x = positions(3*(i-1)+1)
      mdl%elements(i)%centre%y = positions(3*(i-1)+2)
      mdl%elements(i)%centre%z = positions(3*(i-1)+3)

      ! rc2013: also load into regions list (only 1 subsystem allowed with NEB)
      mdl%regions(1)%elements(i)%centre%x = positions(3*(i-1)+1)
      mdl%regions(1)%elements(i)%centre%y = positions(3*(i-1)+2)
      mdl%regions(1)%elements(i)%centre%z = positions(3*(i-1)+3)
    end do
  end subroutine neb_load_coords_into_mdl

!---------------------------------------------------------------------------------------
! Misc. Helper Routines
!---------------------------------------------------------------------------------------

  subroutine neb_subtract(a,b,mdl,result)
    !====================================================================!
    ! Given full configurations a and b, returns a-b with the minimum    !
    ! image convention in the cell contained in mdl.                     !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Apr 2018                                    !
    !====================================================================!

    use comms, only: pub_on_root
    use rundat, only: pub_debug_on_root

    ! Arguments
    real(kind=dp), dimension(:), intent(in   ) :: a, b
    type(MODEL),                 intent(in   ) :: mdl
    real(kind=dp), dimension(:), intent(inout) :: result

    ! Local variables
    integer :: i, j, k, l, ierr
    real(kind=dp) :: distance, min_dist
    real(kind=dp), dimension(3) :: temp_vec, best_vec

    result = a - b
    min_dist = sqrt(dot_product(result,result))

    do i=1,size(a),3
      ! For each atom
      temp_vec = result(i:i+2)
      min_dist = sqrt(dot_product(temp_vec,temp_vec))
      best_vec = temp_vec
      do j=-1,1
        do k=-1,1
          do l=-1,1
            ! For each translation
            ! Add each unit cell component
            temp_vec(1) = result(i+0) + real(j,kind=dp) * mdl%cell%a1%X &
                                      + real(k,kind=dp) * mdl%cell%a2%X &
                                      + real(l,kind=dp) * mdl%cell%a3%X
            temp_vec(2) = result(i+1) + real(j,kind=dp) * mdl%cell%a1%Y &
                                      + real(k,kind=dp) * mdl%cell%a2%Y &
                                      + real(l,kind=dp) * mdl%cell%a3%Y
            temp_vec(3) = result(i+2) + real(j,kind=dp) * mdl%cell%a1%Z &
                                      + real(k,kind=dp) * mdl%cell%a2%Z &
                                      + real(l,kind=dp) * mdl%cell%a3%Z
            distance = sqrt(dot_product(temp_vec,temp_vec))
            if (distance < min_dist) then
              min_dist = distance
              best_vec = temp_vec
            end if
          end do
        end do
      end do
      result(i:i+2) = best_vec
    end do

  end subroutine neb_subtract

  subroutine neb_wrap(pos, mdl)
    !====================================================================!
    ! Wraps the given set of positions back into the simulation cell.    !
    ! Currently it's pretty lazily implemented, creating a castep model  !
    !  just to call cart_to_frac etc...                                  !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Apr 2018                                    !
    !====================================================================!

    use comms, only: pub_on_root
    use model_type, only: MODEL
    use simulation_cell, only: castep_cell_frac_to_cart, &
         castep_cell_cart_to_frac, castep_cell, copy_to_castep_cell

    ! Arguments
    real(kind=dp), dimension(:), intent(inout) :: pos
    type(MODEL),                 intent(in   ) :: mdl

    ! Local variables
    type(CASTEP_CELL) :: temp_cell
    real(kind=dp), dimension(3) :: frac_pos
    integer :: i, j, ierr

    call copy_to_castep_cell(temp_cell,mdl%cell,mdl%nat,mdl%elements, &
             mdl%nat_classical,mdl%classical_elements)

    do i=1,size(pos),3
      call castep_cell_cart_to_frac(temp_cell, pos(i:i+2), frac_pos)
      do j=1,3
        frac_pos(j) = MODULO(frac_pos(j),1.0_DP)
      end do
      call castep_cell_frac_to_cart(temp_cell, frac_pos, pos(i:i+2))
    end do

  end subroutine neb_wrap

!---------------------------------------------------------------------------------------
! File I/O Routines
!---------------------------------------------------------------------------------------

  subroutine neb_cont_write(mdl, neb_mep, v, iteration)
    !====================================================================!
    ! Write out a continuation file for the NEB path.                    !
    ! This includes:                                                     !
    ! - Update Method                                                    !
    ! - Iteration number                                                 !
    ! - Bead positions                                                   !
    ! - Bead velocities                                                  !
    ! - Internal optimizer variables                                     !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Oct 2018                                    !
    !====================================================================!
    use image_comms, only: pub_my_image, pub_on_imroots_root
    use utils, only: utils_unit, utils_open_unit_check, &
                     utils_close_unit_check, utils_write_check, &
                     utils_alloc_check, utils_dealloc_check
    use rundat, only: pub_num_images

    ! Arguments
    type(MODEL),                 intent(in   ) :: mdl
    type(NEB_PATH),              intent(in   ) :: neb_mep
    real(kind=dp), dimension(:), intent(in   ) :: v
    integer,                     intent(in   ) :: iteration

    ! Local variables
    integer :: iunit, ios, i, ierr
    real(kind=dp), dimension(:), allocatable :: r
    character(len=80) :: bead_string, filename

    allocate(r(3*neb_mep%nat),stat=ierr)
    call utils_alloc_check("neb_cont_write","Bead Position",ierr)

    write(bead_string,*)pub_my_image
    filename = "neb_"//trim(adjustl(bead_string))//".neb_cont"

    iunit = utils_unit()

    open(unit=iunit,file=trim(filename),iostat=ios, &
            form='UNFORMATTED',action='WRITE')
    call utils_open_unit_check('geom_opt_continuation_write',filename,ios)
    rewind(iunit)

    call neb_extract_coords_from_mdl(mdl, r)

    ! Chain Properties
    write(iunit, iostat=ios) neb_mep%update_method
    call utils_write_check('neb_cont_write','update_method',ios)
    write(iunit, iostat=ios) iteration
    call utils_write_check('neb_cont_write','iteration number',ios)
    write(iunit, iostat=ios) r
    call utils_write_check('neb_cont_write','bead position',ios)
    write(iunit, iostat=ios) v
    call utils_write_check('neb_cont_write','bead velocity',ios)

    ! FIRE Related
    write(iunit, iostat=ios) neb_mep%cut
    call utils_write_check('neb_cont_write','cut',ios)
    write(iunit, iostat=ios) neb_mep%fire_alpha
    call utils_write_check('neb_cont_write','fire_alpha',ios)
    write(iunit, iostat=ios) neb_mep%fire_dt
    call utils_write_check('neb_cont_write','fire_dt',ios)

    if (pub_on_imroots_root) then
      ! GL-BFGS Related - only computed on root proc
      write(iunit, iostat=ios) neb_mep%lbfgs_q, neb_mep%lbfgs_z, &
                                neb_mep%old_q, neb_mep%old_z, &
                                neb_mep%lbfgs_rho, neb_mep%lbfgs_alpha, &
                                neb_mep%lbfgs_s, neb_mep%lbfgs_y
      call utils_write_check('neb_cont_write','GL-BFGS Arrays',ios)
    end if

    close(iunit, iostat=ios)
    call utils_close_unit_check('neb_cont_write','iunit',ios)

    deallocate(r,stat=ierr)
    call utils_dealloc_check("neb_cont_write","Bead Position",ierr)

  end subroutine neb_cont_write

  subroutine neb_cont_read(mdl, neb_mep, v, iteration)
    !====================================================================!
    ! Read in a continuation file for the NEB path.                      !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Oct 2018                                    !
    !====================================================================!
    use image_comms, only: pub_my_image, pub_on_imroots_root
    use utils, only: utils_unit, utils_open_unit_check, &
                     utils_close_unit_check, utils_read_check, &
                     utils_alloc_check, utils_dealloc_check
    use rundat, only: pub_num_images

    ! Arguments
    type(MODEL),                 intent(inout) :: mdl
    type(NEB_PATH),              intent(inout) :: neb_mep
    real(kind=dp), dimension(:), intent(  out) :: v
    integer,                     intent(  out) :: iteration

    ! Local variables
    integer :: iunit, ios, i, ierr
    real(kind=dp), dimension(:), allocatable :: r
    character(len=80) :: bead_string, filename

    allocate(r(3*neb_mep%nat),stat=ierr)
    call utils_alloc_check("neb_cont_read","Bead Position",ierr)

    ! Find available unit specifier for reading continuation file
    iunit = utils_unit()

    write(bead_string,*)pub_my_image
    filename = "neb_"//trim(adjustl(bead_string))//".neb_cont"

    ! Open continuation file
    open(unit=iunit,file=filename,iostat=ios,&
         form='UNFORMATTED',action='READ',status='OLD')
    call utils_open_unit_check('geom_opt_continuation_read',filename,ios)
    rewind(iunit)

    ! Chain Properties
    read(iunit, iostat=ios) neb_mep%update_method
    call utils_read_check('neb_cont_read','update_method',ios)
    read(iunit, iostat=ios) iteration
    call utils_read_check('neb_cont_read','iteration number',ios)
    read(iunit, iostat=ios) r
    call utils_read_check('neb_cont_read','bead position',ios)
    read(iunit, iostat=ios) v
    call utils_read_check('neb_cont_read','bead velocity',ios)

    ! FIRE Related
    read(iunit, iostat=ios) neb_mep%cut
    call utils_read_check('neb_cont_read','cut',ios)
    read(iunit, iostat=ios) neb_mep%fire_alpha
    call utils_read_check('neb_cont_read','fire_alpha',ios)
    read(iunit, iostat=ios) neb_mep%fire_dt
    call utils_read_check('neb_cont_read','fire_dt',ios)

    if (pub_on_imroots_root) then
      ! GL-BFGS Related - only computed on root proc
      read(iunit, iostat=ios) neb_mep%lbfgs_q, neb_mep%lbfgs_z, &
                                neb_mep%old_q, neb_mep%old_z, &
                                neb_mep%lbfgs_rho, neb_mep%lbfgs_alpha, &
                                neb_mep%lbfgs_s, neb_mep%lbfgs_y
      call utils_read_check('neb_cont_read','GL-BFGS Arrays',ios)
    end if

    call neb_load_coords_into_mdl(mdl, r)

    if (neb_mep%ci_delay >= 0) then
      neb_mep%ci_delay = max(0,neb_mep%ci_delay - iteration)
    else ! No climbing image
      neb_mep%ci_delay = neb_mep%ci_delay - iteration
    end if

    close(iunit, iostat=ios)
    call utils_close_unit_check('neb_cont_read','iunit',ios)

    deallocate(r,stat=ierr)
    call utils_dealloc_check("neb_cont_write","Bead Position",ierr)
  end subroutine neb_cont_read

  subroutine neb_print_path(neb_mep)
    !====================================================================!
    ! Print out normalized reaction coordinates and image energies       !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Mar 2018                                    !
    !====================================================================!

    use comms, only: pub_on_root, comms_send, comms_recv, pub_imroots_comm, &
         comms_barrier
    use image_comms, only: orig_stdout, pub_imroots_root_id, &
         pub_my_rank_in_imroots, pub_on_imroots_root
    use rundat, only: pub_num_images, pub_print_qc
    use services, only: services_flush
    implicit none

    ! Arguments
    type(NEB_PATH), intent(in   ) :: neb_mep

    ! Local variables
    character(len=80) :: status_divider_string, status_label_string, &
                        status_announce_string, data_string
    character(len=3)  :: is_climb
    real(kind=DP) :: mep_len, coord_frag, ts_energy
    real(kind=DP), dimension(1:pub_num_images) :: im_coords
    integer :: image_label_len, dble_output_len, units_len, bool_len, &
              coord_output_len, string_index, i, ts_loc
    integer, dimension(1) :: temp_loc
    character(len=24) :: outtag

    status_announce_string = '------------------ Reaction Pathway ------------------'
    !                         1234567 234567890 23456789012345678 234567890123 23456
    status_divider_string  = '+-----+---------+-----------------+------------+-----+'
    status_label_string    = '| Img |  Coord  |   rel. energy   |    Units   | CI? |'

    image_label_len = 7
    coord_output_len = 10
    dble_output_len = 18
    units_len = 13
    bool_len = 6

    data_string = repeat(' ',len(data_string))

    if (pub_on_imroots_root) then
      ! NOTE: this might not work if a 1-bead run becomes supported.
      do i = 0, pub_num_images - 1
        if (i /= pub_my_rank_in_imroots) then
          call comms_recv(i, im_coords(i+1), tag=400+i, comm=pub_imroots_comm)
        else
          im_coords(i+1) = sqrt(dot_product(neb_mep%prev_to_cur,neb_mep%prev_to_cur))
        end if
      end do

      call comms_recv(pub_num_images-1, coord_frag, tag=500+pub_num_images-1, comm=pub_imroots_comm)

      mep_len = coord_frag + SUM(im_coords)
      do i = 2, pub_num_images
        im_coords(i) = im_coords(i-1) + im_coords(i)
      end do
    else if (pub_on_root) then
      coord_frag = sqrt(dot_product(neb_mep%prev_to_cur,neb_mep%prev_to_cur))
      call comms_send(pub_imroots_root_id, coord_frag, tag=400+pub_my_rank_in_imroots, comm=pub_imroots_comm)

      if (pub_my_rank_in_imroots == pub_num_images - 1) then
        coord_frag = sqrt(dot_product(neb_mep%cur_to_next,neb_mep%cur_to_next))
        call comms_send(pub_imroots_root_id, coord_frag, tag=500+pub_my_rank_in_imroots, comm=pub_imroots_comm)
      end if
    end if

    if (pub_on_imroots_root) then
      ts_energy = maxval(neb_mep%im_energies(:)) - neb_mep%reactant_energy
      temp_loc  = maxloc(neb_mep%im_energies(:))
      ts_loc = temp_loc(1) - 1

      write(orig_stdout,*)""
      write(orig_stdout,*)""
      write(orig_stdout,'(a)')status_announce_string
      write(orig_stdout,'(a)')status_divider_string
      write(orig_stdout,'(a)')status_label_string
      write(orig_stdout,'(a)')status_divider_string
      data_string = repeat(' ',len(data_string))
      ! Construct the data string for the reactant
      string_index = 1
      write(data_string(string_index:string_index+image_label_len),'(a,a4,1x,a)')'|','   R','|'
      string_index=string_index+image_label_len
      write(data_string(string_index:string_index+coord_output_len),'(1x,f7.5,1x,a)')0.0_dp,'|'
      string_index = string_index + coord_output_len
      write(data_string(string_index:string_index+dble_output_len),'(1x,es15.6e3,1x,a)')0.0_dp,'|'
      string_index = string_index + dble_output_len
      write(data_string(string_index:string_index+units_len),'(1x,a10,1x,a)')"Ha",'|'
      string_index = string_index + units_len
      write(data_string(string_index:string_index+bool_len),'(1x,a3,1x,a)')"No",'|'
      write(orig_stdout,'(a)')data_string

      do i=0,pub_num_images-1

        if (i == neb_mep%climbing_image) then
          is_climb = "Yes"
        else
          is_climb = " No"
        end if

        data_string = repeat(' ',len(data_string))
        ! Construct the data string for the product
        string_index = 1
        write(data_string(string_index:string_index+image_label_len),'(a,i4,1x,a)')'|',i,'|'
        string_index=string_index+image_label_len
        write(data_string(string_index:string_index+coord_output_len),'(1x,f7.5,1x,a)') &
                                        im_coords(i+1)/mep_len,'|'
        string_index = string_index + coord_output_len
        write(data_string(string_index:string_index+dble_output_len),'(1x,es15.6e3,1x,a)') &
                                        neb_mep%im_energies(i+1) - neb_mep%reactant_energy,'|'
        string_index = string_index + dble_output_len
        write(data_string(string_index:string_index+units_len),'(1x,a10,1x,a)')"Ha",'|'
        string_index = string_index + units_len
        write(data_string(string_index:string_index+bool_len),'(1x,a3,1x,a)')is_climb,'|'
        write(orig_stdout,'(a)')data_string
      end do

      data_string = repeat(' ',len(data_string))
      ! Construct the data string for the product
      string_index = 1
      write(data_string(string_index:string_index+image_label_len),'(a,a4,1x,a)')'|','   P','|'
      string_index=string_index+image_label_len
      write(data_string(string_index:string_index+coord_output_len),'(1x,f7.5,1x,a)')1.0_dp,'|'
      string_index = string_index + coord_output_len
      write(data_string(string_index:string_index+dble_output_len),'(1x,es15.6e3,1x,a)')&
                                      neb_mep%product_energy - neb_mep%reactant_energy,'|'
      string_index = string_index + dble_output_len
      write(data_string(string_index:string_index+units_len),'(1x,a10,1x,a)')"Ha",'|'
      string_index = string_index + units_len
      write(data_string(string_index:string_index+bool_len),'(1x,a3,1x,a)')"No",'|'
      write(orig_stdout,'(a)')data_string

      write(orig_stdout,'(a)')status_divider_string

      ! Also write out the current TS info
      ! We should handle the unfortunate cases that the reactant or product are our highest energy image
      write(orig_stdout,'(1x,a,1x,i4,1x,a,1x,f7.5,1x,a)')"TS is bead",ts_loc,"with relative energy",ts_energy,"Ha."
      call services_flush(orig_stdout)

      ! QC Test Output. Unfortunately as we're not using stdout we can't use the utils routine.
      if (pub_print_qc) then
        outtag = '[first_im_en]'
        write(orig_stdout,'(a30,f23.12)') '<QC> '//adjustr(outtag)//':', neb_mep%im_energies(1) &
                                                                 - neb_mep%reactant_energy
        outtag = '[last_im_en]'
        write(orig_stdout,'(a30,f23.12)') '<QC> '//adjustr(outtag)//':', neb_mep%im_energies(pub_num_images) &
                                                                - neb_mep%reactant_energy
      end if
    end if

  end subroutine neb_print_path

  subroutine neb_print_conv(conv, iteration, dE, Fmax, dRmax, neb_mep)
    !====================================================================!
    ! Print out normalized reaction coordinates and image energies       !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Mar 2018                                    !
    !====================================================================!

    use comms, only: pub_on_root, comms_send, comms_recv, pub_imroots_comm
    use image_comms, only: pub_on_imroots_root, orig_stdout, &
         pub_imroots_root_id, pub_my_rank_in_imroots
    use rundat, only: pub_num_images, tssearch_energy_tol, tssearch_force_tol, &
         tssearch_disp_tol, pub_neb_converge_all
    use services, only: services_flush
    implicit none

    ! Arguments
    logical,        intent(in   ) :: conv
    integer,        intent(in   ) :: iteration
    real(kind=DP),  intent(in   ) :: dE, Fmax, dRmax
    type(NEB_PATH), intent(in   ) :: neb_mep

    ! Local variables
    real(kind=dp) :: tde, tfmax, tdrmax, dE_gmax, F_gmax, dR_gmax
    character(len=80) :: conv_divider_string, conv_label_string, &
                        data_string, conv_string, &
                        conv_announce_string
    integer :: image_label_len, dble_output_len, units_len, bool_len, &
              coord_output_len, i, string_index
    logical :: loc_write_disp_energy_tols

    loc_write_disp_energy_tols = pub_neb_converge_all &
                                 .or. (neb_mep%climbing_image .ge. 0)

    de_gmax = 0.0_dp
    F_gmax = 0.0_dp
    dR_gmax = 0.0_dp

    data_string     = repeat(' ',len(data_string))
    conv_divider_string  = repeat(' ',len(conv_divider_string))
    conv_label_string    = repeat(' ',len(conv_label_string))

    ! Formatting stuff - partly taken from geomopt

    conv_announce_string   = '------------------------ Convergence ------------------------'
    !                         1234567 23456789012345678 23456789012345678 23456789012345678 23456
    conv_divider_string    = '+-----+-----------------+-----------------+-----------------+'
    conv_label_string      = '| Img |       dE        |      Fmax       |      dRmax      |'

    image_label_len = 7
    coord_output_len = 17
    dble_output_len = 18
    units_len = 13
    bool_len = 6

    if (conv) then
      conv_string = "converged"
    else
      !              1234567890123
      conv_string = "not converged"
    end if

    if (pub_on_imroots_root) then
      write(orig_stdout,*)""
      write(orig_stdout,*)""
      write(orig_stdout,'(a,a,a,i4,a,a13)')" NEB: Finished ",trim(adjustl(neb_mep%update_method)),&
                                          " iteration ",iteration,', ',conv_string
      write(orig_stdout,'(a)')conv_announce_string
      write(orig_stdout,'(a)')conv_divider_string
      write(orig_stdout,'(a)')conv_label_string
      write(orig_stdout,'(a)')conv_divider_string
      do i=0,pub_num_images-1
        if (i /= pub_my_rank_in_imroots) then
          call comms_recv(i, tde, tag=100+i, comm=pub_imroots_comm)
          call comms_recv(i, tfmax, tag=200+i, comm=pub_imroots_comm)
          call comms_recv(i, tdrmax, tag=300+i, comm=pub_imroots_comm)
        else
          tde = dE
          tfmax = Fmax
          tdrmax = dRmax
        end if
        if ((.not.pub_neb_converge_all) .and. (neb_mep%climbing_image .gt. 0)) then
          ! If we're only converging dE and dRmax of the climbing image
          if (i == neb_mep%climbing_image) then
            dE_gmax = tde
            dR_gmax = tdrmax
          end if
        else
          if (tde > dE_gmax) dE_gmax = tde
          if (tdrmax > dR_gmax) dR_gmax = tdrmax
        end if

        if (tfmax > F_gmax) F_gmax = tfmax

        data_string = repeat(' ',len(data_string))
        ! Construct the data string
        string_index = 1
        write(data_string(string_index:string_index+image_label_len),'(a,i4,1x,a)')'|',i,'|'
        string_index=string_index+image_label_len
        write(data_string(string_index:string_index+dble_output_len),'(1x,es15.6e3,1x,a)')tde,'|'
        string_index=string_index+dble_output_len
        write(data_string(string_index:string_index+dble_output_len),'(1x,es15.6e3,1x,a)')tfmax,'|'
        string_index=string_index+dble_output_len
        write(data_string(string_index:string_index+dble_output_len),'(1x,es15.6e3,1x,a)')tdrmax,'|'
        write(orig_stdout,'(a)')data_string
      end do
      write(orig_stdout,'(a)')conv_divider_string

      ! Write out the max value across all images
      data_string = repeat(' ',len(data_string))
      string_index = 1
      write(data_string(string_index:string_index+image_label_len),'(a,a,1x,a)')'|',' max','|'
      string_index=string_index+image_label_len
      write(data_string(string_index:string_index+dble_output_len),'(1x,es15.6e3,1x,a)')dE_gmax,'|'
      string_index=string_index+dble_output_len
      write(data_string(string_index:string_index+dble_output_len),'(1x,es15.6e3,1x,a)')F_gmax,'|'
      string_index=string_index+dble_output_len
      write(data_string(string_index:string_index+dble_output_len),'(1x,es15.6e3,1x,a)')dR_gmax,'|'
      write(orig_stdout,'(a)')data_string
      write(orig_stdout,'(a)')conv_divider_string

      ! Write out the tolerances
      data_string = repeat(' ',len(data_string))
      string_index = 1
      write(data_string(string_index:string_index+image_label_len),'(a,a,1x,a)')'|',' tol','|'
      string_index=string_index+image_label_len
      if (loc_write_disp_energy_tols) then
        write(data_string(string_index:string_index+dble_output_len),'(1x,es15.6e3,1x,a)')tssearch_energy_tol,'|'
      else
        write(data_string(string_index:string_index+dble_output_len),'(1x,a,1x,a)')'      ---      ','|'
      end if
      string_index=string_index+dble_output_len
      write(data_string(string_index:string_index+dble_output_len),'(1x,es15.6e3,1x,a)')tssearch_force_tol,'|'
      string_index=string_index+dble_output_len
      if (loc_write_disp_energy_tols) then
        write(data_string(string_index:string_index+dble_output_len),'(1x,es15.6e3,1x,a)')tssearch_disp_tol,'|'
      else
        write(data_string(string_index:string_index+dble_output_len),'(1x,a,1x,a)')'      ---      ','|'
      end if
      write(orig_stdout,'(a)')data_string
      write(orig_stdout,'(a)')conv_divider_string

      if ((.not.pub_neb_converge_all) .and. (neb_mep%climbing_image .gt. 0)) then
        ! We don't usually write out energy and disp tol, but because we've got a CI, write them
         !                                                      123456789012345678901
        write(orig_stdout,'(a)')"(Energy and displacement only converged for climbing image)"
      end if
      call services_flush(orig_stdout)
    else if (pub_on_root) then
      call comms_send(pub_imroots_root_id, dE, tag=100+pub_my_rank_in_imroots, comm=pub_imroots_comm)
      call comms_send(pub_imroots_root_id, Fmax, tag=200+pub_my_rank_in_imroots, comm=pub_imroots_comm)
      call comms_send(pub_imroots_root_id, dRmax, tag=300+pub_my_rank_in_imroots, comm=pub_imroots_comm)
    end if

  end subroutine neb_print_conv

  subroutine neb_path_xyz(mdl, neb_mep)
    !====================================================================!
    ! Write outh neb_path.xyz, rewriting any previously written file.    !
    !--------------------------------------------------------------------!
    ! Written by Kevin Duff, Aug 2018                                    !
    !====================================================================!

    use rundat, only: pub_num_images
    use services, only: services_write_xyz
    use utils, only: utils_unit, utils_alloc_check, utils_dealloc_check, &
                     utils_assert
    use image_comms, only: pub_on_imroots_root
    implicit none

    ! Arguments
    type(MODEL),    intent(inout) :: mdl
    type(NEB_PATH), intent(inout) :: neb_mep

    ! Local variables
    integer :: i, iunit, ierr
    real(kind=dp), dimension(:), allocatable :: orig_coords, tmp_coords
    character(len=80) :: header_string
    character(len=12) :: bead_string

    ! Do nothing if we're not the master proc
    if (.not.pub_on_imroots_root) return

    allocate(orig_coords(3*neb_mep%nat), tmp_coords(3*neb_mep%nat), stat=ierr)
    call utils_alloc_check("neb_path_xyz", "coordinate arrays", ierr)

    ! Open and clear the file
    iunit = utils_unit()
    open(unit=iunit, file="neb_path.xyz", status="replace",iostat=ierr)
    call utils_assert(ierr == 0, "Error in neb_path_xyz: Failed to clear neb_path.xyz")
    close(iunit)

    ! Back up the image's model, as we'll be using it to write the path xyzs
    call neb_extract_coords_from_mdl(mdl,orig_coords)

    do i=1,pub_num_images+2
      if (i == 1) then
        header_string = "Reactant"
      else if (i == pub_num_images + 2) then
        header_string = "Product"
      else
        write(bead_string,*)i-2
        header_string = "Bead "//trim(adjustl(bead_string))
      end if
      tmp_coords = neb_mep%im_positions(:,i)
      call neb_load_coords_into_mdl(mdl, tmp_coords)
      call services_write_xyz(mdl%elements, "neb_path", trim(header_string))
    end do

    ! Return the image coordinates to the mdl
    call neb_load_coords_into_mdl(mdl,orig_coords)

    deallocate(orig_coords, tmp_coords, stat=ierr)
    call utils_dealloc_check("neb_path_xyz", "coordinate arrays", ierr)

  end subroutine neb_path_xyz

  subroutine neb_cleanup(mep)
    use image_comms, only: pub_on_imroots_root
    use utils, only: utils_dealloc_check
    implicit none

    ! Arguments
    type(NEB_PATH),intent(inout) :: mep

    ! Local Variables
    integer :: ierr

    if (pub_on_imroots_root) then
      deallocate(mep%im_positions, stat=ierr)
      call utils_dealloc_check('neb_init', 'mep%im_positions', ierr)
      deallocate(mep%old_positions, stat=ierr)
      call utils_dealloc_check('neb_init', 'mep%old_positions', ierr)
      deallocate(mep%im_forces, stat=ierr)
      call utils_dealloc_check('neb_init', 'mep%im_forces', ierr)
      deallocate(mep%im_energies, stat=ierr)
      call utils_dealloc_check('neb_init', 'mep%im_energies', ierr)
    end if

    deallocate(mep%prev_to_cur, mep%cur_to_next, stat=ierr)
    call utils_dealloc_check('neb_driver', 'neighbor vectors', ierr)
    deallocate(mep%prev_pos, mep%next_pos, stat=ierr)
    call utils_dealloc_check('neb_init', 'previous and next positions',ierr)

    deallocate(mep%lbfgs_s, mep%lbfgs_y, mep%lbfgs_q, mep%lbfgs_z, &
             mep%old_q, mep%old_z, mep%lbfgs_rho, mep%lbfgs_alpha, stat=ierr)
    call utils_dealloc_check("neb_optimize","G-LBFGS arrays",ierr)
  end subroutine

end module neb
