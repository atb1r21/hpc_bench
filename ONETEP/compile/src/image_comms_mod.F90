! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!==================================================================!
!                                                                  !
!               Image-Parallel Communications Module               !
!                                                                  !
!------------------------------------------------------------------!
! This module contains comms infrastructure that allows multiple   !
! ONETEP simulations to be running in the same MPI comm world.     !
! This allows for communications between simulations, for example  !
! for chain-of-states methods like NEB where the force on one bead !
! is affected by properties of connected beads.                    !
!                                                                  !
! When initialized, communication within an image is handled with  !
! pub_image_comm, which becomes the default communicator across    !
! ONETEP. Communication between root procs of images is through    !
! pub_imroots_comm.                                                !
!------------------------------------------------------------------!
! Written by Kevin Duff                                            !
!==================================================================!

module image_comms
  use comms, only: pub_image_comm, pub_imroots_comm
  use constants, only : DP, stdout, LONG, stderr
#if defined(MPI) && !defined(USE_INCLUDE_MPIF)
  use mpi !! External dependency
#endif
  implicit none

#ifdef MPI
#ifdef USE_INCLUDE_MPIF
#include "mpif.h"
#endif
#else
  integer, parameter :: MPI_COMM_NULL = 0
  integer, parameter :: MPI_COMM_WORLD = -1
  integer, parameter :: MPI_ANY_TAG = -1
#endif

  private

  ! Public Variables

  integer, public :: orig_stdout
  integer, public :: orig_stderr

  !Information about separate images
  logical, public :: pub_comms_images_initialized
  integer, public :: pub_image_size
  integer, public :: pub_my_image
  integer, public :: pub_my_rank_in_image

  !Information about images intercommunicator
  logical, public :: pub_on_imroots_root
  integer, public :: pub_my_rank_in_imroots
  integer, public :: pub_imroots_root_id

  !Comms World information viewed from one image
  integer, public :: pub_comms_world_size
  integer, public :: pub_my_rank_in_world

  integer,allocatable,dimension(:),public :: pub_roots_in_world

  !Misc. information
  character(len=80), public :: pub_base_rootname

  ! Public Subroutines
  public :: image_comms_init
  public :: image_comms_finalize

  contains

  subroutine image_comms_init()

    !=========================================================================!
    ! This subroutine creates equally-sized groups of processors, each of     !
    ! which act as a separate ONETEP simulation, complete with distinct root, !
    ! run parameters, task, output files, and networking via pub_image_comm.  !
    ! It also creates pub_imroots_comm on which the root processes of each    !
    ! ONETEP image can communicate. This essentially allows ONETEP to task    !
    ! farm and to complete special tasks that require communicating systems   !
    ! such as path-integral molecular dynamics or nudged elastic band in a    !
    ! scalable manor.                                                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables used:                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   comms_init                                                            !
    !-------------------------------------------------------------------------!
    ! Written by Kevin Duff                                                   !
    !=========================================================================!

    use comms, only: comms_set_default_comm, pub_total_num_procs, &
         pub_my_proc_id, pub_root_proc_id, pub_on_root, comms_barrier
    use rundat, only: pub_num_images, pub_rootname, pub_debug, &
         pub_debug_on_root, pub_image_sizes
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_unit, utils_update_rootname

    ! Local Variables
    integer :: ierr
    integer :: i,j
    integer :: world_group, imroots_group
    integer,dimension(:),allocatable :: image_sizes
    character(len=30) :: imagechar
    integer :: myzero

    ! Leave if we aren't doing an image-parallel run
    if (pub_num_images < 1) then
       ! Prevent strange interaction later on...
       call utils_abort("Error in image_comms_init: &
            &Please don't request less than one ONETEP image.", &
            pub_num_images)
    end if

#ifndef MPI
    ! Set variables to avoid warnings
    pub_comms_images_initialized = .false.
    pub_image_size = 1
    pub_my_image = 0
    pub_my_rank_in_image = pub_my_proc_id
    pub_on_imroots_root = pub_on_root
    pub_imroots_root_id = pub_root_proc_id
    pub_comms_world_size = -1
    pub_my_rank_in_world = pub_my_proc_id
    pub_base_rootname = pub_rootname
    allocate(pub_roots_in_world(1))
    pub_roots_in_world = pub_my_proc_id

    return
#else

    if (pub_num_images == 1) then
       pub_comms_images_initialized = .false.
       return
    end if

    call comms_barrier

    allocate(image_sizes(pub_num_images),stat=ierr)
    call utils_alloc_check('image_comms_init','image_sizes',ierr)

    if (trim(adjustl(pub_image_sizes)) /= "DEFAULT") then
       ! Parse image size string
       call parse_image_sizes(image_sizes,pub_image_sizes)

    else
       if (mod(pub_total_num_procs, pub_num_images) .ne. 0) then
          call utils_abort("Error in image_comms_init: &
               &Total number of MPI processes must be divisible by num_images:",&
               pub_total_num_procs, pub_num_images)
       end if

       image_sizes = pub_total_num_procs / pub_num_images

    end if

    ! Note down all the image root procs in mpi_comm_world
    ! TODO make this work for arbitrary root proc and image distribution
    allocate(pub_roots_in_world(pub_num_images),stat=ierr)
    call utils_alloc_check('image_comms_init','pub_roots_in_world',ierr)

    myzero = 0
    do i=1,pub_num_images
       pub_roots_in_world(i) = myzero
       myzero = myzero + image_sizes(i)
    end do
    !if(pub_on_root)write(stdout,*)pub_roots_in_world

    if (pub_on_root) then
       write(stdout,'(a,i0,a)')"Initializing image-parallel run with ", &
            pub_num_images," images."
    end if

    ! Now create the image communicators
    j = 0
    i = pub_my_proc_id
    do
       i = i - image_sizes(j+1)
       if (i < 0) then
          pub_my_image = j
          exit
       end if
       j = j + 1
    end do

    call mpi_comm_split(mpi_comm_world, pub_my_image, &
         pub_my_proc_id, pub_image_comm, ierr)
    if (ierr .ne. 0) then
       call utils_abort('Error in image_comms_init: &
            &Could not create image communicators')
    end if

    call mpi_comm_rank(pub_image_comm, pub_my_rank_in_image, ierr)
    if (ierr .ne. 0) then
       call utils_abort('Error in image_comms_init: &
            &Could not get image communicator ranks')
    end if

    pub_image_size = image_sizes(pub_my_image+1)

    call mpi_comm_group(mpi_comm_world, world_group, ierr)
    if (ierr .ne. 0) then
       call utils_abort('Error in image_comms_init:&
            & Could not create world group')
    end if

    ! Finally create the group and its communicator
    pub_imroots_comm = MPI_COMM_NULL
    call mpi_group_incl(world_group, pub_num_images, pub_roots_in_world, &
         imroots_group, ierr)
    if (ierr .ne. 0) then
       call utils_abort('Error in image_comms_init:&
            & Could not create imroots group')
    end if
    call mpi_comm_create(mpi_comm_world, imroots_group, pub_imroots_comm, ierr)
    if (ierr .ne. 0) then
       call utils_abort('Error in image_comms_init:&
            & Could not create imroots communicator')
    end if

    pub_my_rank_in_imroots = -1
    if (pub_imroots_comm.ne.MPI_COMM_NULL) then
       call mpi_comm_rank(pub_imroots_comm,pub_my_rank_in_imroots,ierr)
       if (ierr .ne. 0) then
          call utils_abort('Error in image_comms_init: &
               &Could not get imroots ranks')
       end if
    end if

    pub_imroots_root_id = 0
    pub_on_imroots_root = .false.
    if (pub_my_rank_in_imroots == pub_imroots_root_id) then
       pub_on_imroots_root = .true.
    end if

    ! Change default communicator
    call comms_set_default_comm(pub_image_comm)

    ! Give each image a root process
    pub_on_root = .false.
    if (pub_my_rank_in_image == 0) then
       ! TODO get rid of this test
       if(pub_my_rank_in_imroots == -1)write(stdout,*)"Problem - not in imroots"
       pub_on_root = .true.
       pub_debug_on_root = pub_debug
    end if

    ! Change other comms-related variables to insulate images
    pub_comms_world_size = pub_total_num_procs
    pub_my_rank_in_world = pub_my_proc_id
    pub_total_num_procs  = pub_image_size
    pub_my_proc_id       = pub_my_rank_in_image
    pub_root_proc_id     = 0

    ! kkbd: Preserve the ability to write to stdout
    orig_stdout = stdout
    orig_stderr = stderr

    stdout = utils_unit()
    stderr = stdout

    ! Give each image its own set of output files
    write(imagechar,*)pub_my_image
    pub_base_rootname = pub_rootname
    pub_rootname = trim(adjustl(pub_rootname))//trim(adjustl(imagechar))
    ! ebl: make sure utils_mod also knows that this has been updated
    call utils_update_rootname(pub_rootname, delete_err_file = .true.)

    if (pub_on_root) then
       open(file=trim(pub_rootname)//trim(".onetep"),unit=stdout)
    end if

    if (pub_on_root) then
       write(stdout,*)"Initialized image",pub_my_image," with ",pub_image_size,"processes."
    end if

    deallocate(image_sizes,stat=ierr)
    call utils_dealloc_check('image_comms_init','image_sizes',ierr)

    call comms_barrier(MPI_COMM_WORLD)
    pub_comms_images_initialized = .true.
#endif

  end subroutine image_comms_init

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine image_comms_finalize()

    !=========================================================================!
    ! This subroutine cleans up anything to do with images just before ONETEP !
    ! shuts down. Only really does anything if images were initalized.        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables used:                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Kevin Duff                                                   !
    !=========================================================================!

    use comms, only: pub_on_root, comms_barrier
    use utils, only: utils_dealloc_check

    implicit none

    ! Local Variables
    integer :: ierr

    call comms_barrier(MPI_COMM_WORLD)

    if (pub_comms_images_initialized)then
       deallocate(pub_roots_in_world,stat=ierr)
       call utils_dealloc_check('image_comms_finalize','pub_roots_in_world',ierr)
       if (pub_on_root) close(stdout)
       pub_comms_images_initialized = .false.
    end if

  end subroutine image_comms_finalize

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine parse_image_sizes(image_array, size_string)

    !=========================================================================!
    use comms, only: pub_total_num_procs
    use rundat, only: pub_num_images
    use utils, only: utils_abort

    implicit none

    ! Local Variables
    character(len=*),    intent(in   ) :: size_string
    integer,dimension(:),intent(inout) :: image_array
    character(len=80) :: string
    character(len=80) :: temp_string
    integer :: i, delim_index

    image_array = 0

    string = adjustl(trim(size_string))

    do i=1,pub_num_images
       delim_index = index(string, "|")
       if (delim_index == 0) then
          ! Potentially do something if we end too early? kkbd
          temp_string = string
       else
          temp_string = string(1:delim_index-1)
       end if

       read(temp_string,*) image_array(i)

       string = string(delim_index+1:)

    end do

    if (sum(image_array) /= pub_total_num_procs) then
       call utils_abort("Error in image_comms_init:&
            & Image sizes must add up to number of MPI procs!")
    end if

  end subroutine parse_image_sizes

end module image_comms
