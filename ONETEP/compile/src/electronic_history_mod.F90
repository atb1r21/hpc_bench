! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!============================================================================!
!                                                                            !
!               Electronic history module                                    !
!                                                                            !
!   Subroutines for dealing with the history of electronic degrees of        !
!   freedom (edofs) that are used to create initial guesses for NGWFs and    !
!   density kernel during a molecular dynamic calculation or a geometry      !
!   optimization.                                                            !
!                                                                            !
!----------------------------------------------------------------------------!
!   Written by Simon M.-M. Dubois, March 2013                                !
!                                                                            !
!   TCM Group, Cavendish laboratory                                          !
!   Madingley Road                                                           !
!   Cambridge CB3 0HE                                                        !
!   UK                                                                       !
!============================================================================!


module electronic_history

  use datatypes, only: FUNCTIONS, FFTBOX_DATA
  use constants, only: dp, NORMAL, BRIEF
  use rundat, only: pub_md_output_detail, pub_debug_on_root

  implicit none

  private

  !=== Subroutines ==========================================================!

  ! smmd: New subroutines to keep electronic degrees of freedom
  ! smmd: in memory.
  public :: elec_history_initialise_method
  public :: elec_history_update_method
  public :: elec_history_reset_method
  public :: elec_history_use_method

  public :: elec_history_reset

  public :: elec_history_create_storage
  public :: elec_history_destroy_storage

  public :: elec_history_create_entry
  public :: elec_history_remove_entry

  ! jd (retreat2014): The subroutines that are private need not be listed here,
  !                   since they can't be accessed from outside this module
  !                   anyway, and everything is private by default owing to the
  !                   'private' clause at line 28. Normally it would not matter,
  !                   but gfortran bug 54224 flags such functions as unused and
  !                   generates spurious warnings. So for now I commented them out.
  !                   cf. https://gcc.gnu.org/bugzilla/show_bug.cgi?id=54224

  public :: elec_history_store_ngwfs
!  private :: elec_history_read_ngwfs
  public :: elec_history_compose_ngwfs

  public :: elec_history_store_dkn
!  private :: elec_history_read_dkn
  public :: elec_history_compose_dkn

  public :: elec_history_backup_history
  public :: elec_history_restore_history_info
  public :: elec_history_restore_history_parm
! public :: elec_history_write_ngwfs
! public :: elec_history_write_dkn
!  private :: methods_index
!  private :: history_index
  public :: elec_history_update_dkn
  public :: elec_history_compute_temp

  !=== Type definitions =====================================================!

  ! NGWFs storage format
  type, public  :: ngwfs_storage_format
     integer                         :: tb_size(3)
     integer, allocatable            :: procs(:)
     real(kind=DP), allocatable      :: orig(:,:)
     type(FFTBOX_DATA), allocatable :: first(:)
     type(FFTBOX_DATA), allocatable :: init(:)
     type(FFTBOX_DATA), allocatable :: scf(:)
  end type ngwfs_storage_format

  ! Density kernel storage format
  type, public  :: dkn_storage_format
     integer, allocatable            :: idx(:,:)
     integer, allocatable            :: idxinit(:,:)
     integer, allocatable            :: idxvel(:,:)
     type(FUNCTIONS), allocatable :: init(:)
     type(FUNCTIONS), allocatable :: scf(:)
     type(FUNCTIONS), allocatable :: vel(:)
  end type

  ! Storage format for electronic degrees of freedom
  type, public ::  edf_storage_format
     ! Label of the item
     integer                         :: label
     ! Atomic coordinates
     real(kind=DP), allocatable      :: coord(:,:)
     ! NGWFs
     type(ngwfs_storage_format), allocatable :: ngwfs(:)  ! jcap: make this an array with an element for each region
     ! Density kernel
     type(dkn_storage_format)        :: dkn
     ! Lowdin matrix
  end type

  ! Electronic degrees of freedom composition variables
  type, public :: edf_compo_param
     ! Name of the composition scheme
     character(len=10)                :: name
     ! Number of items to be used in the mixing
     integer                          :: ngwfs_mix_size
     ! Type of mixing
     character(len=10)                :: ngwfs_mix_type
     ! Number of mixing step before reset
     integer                          :: ngwfs_mix_reset
     ! Number of initialization steps
     integer                          :: ngwfs_init_num
     ! Type of initialization steps
     character(len=10)                :: ngwfs_init_type
     ! Mixing coefficient
     real(kind=DP)                    :: ngwfs_coeff
     ! Internal counter
     integer                          :: ngwfs_count
     ! Number of items to be used in the mixing
     integer                          :: dkn_mix_size
     ! Type of mixing
     character(len=10)                :: dkn_mix_type
     ! Number of mixing step before reset
     integer                          :: dkn_mix_reset
     ! Number of initialization steps
     integer                          :: dkn_init_num
     ! Type of initialization steps
     character(len=10)                :: dkn_init_type
     ! Internal counter
     integer                          :: dkn_count
     ! List of items linked to the composition scheme
     integer, allocatable             :: ptr(:)
     ! Maximum number of items required by the scheme
     integer                          :: maxit
     ! Current number of items available
     integer                          :: size
     ! Max radius for local mixing schems
     real(kind=DP)                    :: mix_radius
     ! Smearing length for local mixing schems
     real(kind=DP)                    :: mix_smear
  end type edf_compo_param

  !=== Module variables =====================================================!

  type(edf_storage_format), allocatable, save  :: history(:)
  integer, save  :: history_max = 0
  integer, save  :: history_size = 0
  integer, save  :: history_count = 0
  integer, save  :: history_io = 0

  type(edf_compo_param), allocatable, save  :: methods(:)
  integer, save  :: methods_num = 0
  integer, save  :: methods_max = 4
  integer, save  :: methods_io = 0



contains

!============================================================================!
!============================================================================!
!============================================================================!

  subroutine elec_history_initialise_method(label, dkn_mix_type, &
     dkn_mix_num, dkn_mix_reset, dkn_mix_init_type, dkn_mix_init_num, &
     ngwfs_mix_type, ngwfs_mix_num, ngwfs_mix_reset, ngwfs_mix_init_type, &
     ngwfs_mix_init_num, mix_coeff, mix_radius, mix_smear,&
     denskern_count,ngwf_count)

    !========================================================================!
    ! This subroutine initializes the mixing method as defined by the input  !
    ! variables mix_dkn_num, mix_dkn_type, ....                              !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 27/03/2013                            !
    !========================================================================!

    use constants, only: stdout
    use utils, only: utils_alloc_check,  utils_assert

    implicit none

    ! Parameters
    character(len=*), intent(in) :: label

    integer, intent(in) :: dkn_mix_num
    integer, intent(in) :: ngwfs_mix_num
    integer, intent(in) :: dkn_mix_init_num
    integer, intent(in) :: ngwfs_mix_init_num
    integer, intent(in) :: dkn_mix_reset
    integer, intent(in) :: ngwfs_mix_reset
    real(kind=DP), intent(in) :: mix_radius
    real(kind=DP), intent(in) :: mix_smear
    real(kind=DP), intent(in) :: mix_coeff
    character(len=*), intent(in) :: dkn_mix_type
    character(len=*), intent(in) :: dkn_mix_init_type
    character(len=*), intent(in) :: ngwfs_mix_type
    character(len=*), intent(in) :: ngwfs_mix_init_type
    integer, optional, intent(in) :: denskern_count, ngwf_count

    ! Local variables
    integer :: maxsize
    integer :: idx, ierr

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering elec_history_initialise_method'

    ! If required allocate methods
    if (.not. allocated(methods)) then
       allocate(methods(methods_max),stat=ierr)
       call utils_alloc_check('elec_history_initialise_method','methods',ierr)
    endif

    ! Find an index for the new mixing method
    idx = methods_index()
    call utils_assert(idx <= methods_max, 'Error in elec_history_initialise_method: &
                 &method idx is out of range!')

    ! Initialise composition parameters
    methods(idx)%name = label
    methods(idx)%dkn_count = 0
    methods(idx)%ngwfs_count = 0

    ! Allocate ptr
    maxsize = max(dkn_mix_num, ngwfs_mix_num)
    if ((dkn_mix_init_num > 0 .and. dkn_mix_init_type == 'REUSE') .or. &
        (ngwfs_mix_init_num > 0 .and. ngwfs_mix_init_type == 'REUSE')) then
       maxsize = max(maxsize, 1)
    endif
    methods(idx)%size = 0
    methods(idx)%maxit = maxsize
    if (methods(idx)%maxit > 0) then
       allocate(methods(idx)%ptr(maxsize),stat=ierr)
       call utils_alloc_check('elec_history_initialise_method','methods(idx)%ptr',ierr)
       methods(idx)%ptr(:) = 0
    endif

    ! Fill in mixing parameters
    methods(idx)%dkn_mix_type  = dkn_mix_type
    methods(idx)%dkn_mix_size  = dkn_mix_num
    methods(idx)%dkn_mix_reset = dkn_mix_reset
    methods(idx)%dkn_init_num  = dkn_mix_init_num
    methods(idx)%dkn_init_type = dkn_mix_init_type

    methods(idx)%ngwfs_mix_type  = ngwfs_mix_type
    methods(idx)%ngwfs_mix_size  = ngwfs_mix_num
    methods(idx)%ngwfs_mix_reset = ngwfs_mix_reset
    methods(idx)%ngwfs_init_num  = ngwfs_mix_init_num
    methods(idx)%ngwfs_init_type = ngwfs_mix_init_type
    methods(idx)%ngwfs_coeff = mix_coeff

    methods(idx)%mix_radius = mix_radius
    methods(idx)%mix_smear  = mix_smear

    ! Update elec_mix_num
    methods_num = methods_num + 1

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving elec_history_initialise_method'

    return

  end subroutine elec_history_initialise_method

!============================================================================!
!============================================================================!
!============================================================================!

  subroutine elec_history_use_method(label)

    !========================================================================!
    ! This subroutine sets the methods_io public variable to the method      !
    ! associated with label                                                  !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 27/03/2013                            !
    !========================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout
    use utils, only: utils_assert

    implicit none

    ! Parameters
    character(len=*), intent(in) :: label

    ! Local variables
    integer :: idx

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering elec_history_use_method'

    ! Preliminary checks
    call utils_assert(allocated(methods), 'Error in elec_history_use_method: &
            &no method defined!')

    ! Find an index for the new mixing method
    idx = methods_index(label)
    call utils_assert(idx <= methods_max,'Error in elec_history_use_method: &
                &method idx is out of range!')
    methods_io = idx

    if (pub_on_root .and. pub_md_output_detail>=NORMAL ) &
          write(stdout,'(a,a)') 'Electronic history : &
          &composition method set to ', label

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving elec_history_use_method'

    return

  end subroutine elec_history_use_method

!============================================================================!
!============================================================================!
!============================================================================!

  subroutine elec_history_update_method(label)

    !========================================================================!
    ! This subroutine update the pointers and counters associated with the   !
    ! composition method                                                     !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 27/03/2013                            !
    !========================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout
    use utils, only: utils_assert

    implicit none

    ! Parameters
    character(len=*), intent(in) :: label

    ! Local variables
    integer :: iptr, idx
    logical :: dkn_reset, ngwfs_reset

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Entering elec_history_update_method'

    ! If required allocate elec_method
    call utils_assert(allocated(methods),'Error in elec_history_update_method: &
          &no method defined!')

    ! Find the index of the mixing method
    idx = methods_index(label)
    call utils_assert(idx <= methods_max,'Error in elec_history_update_method: &
          &method idx is out of range!')

    ! If method is empty
    if (methods(idx)%maxit == 0) return

    ! Update the ptr array
    if (methods(idx)%maxit.gt.1) then
       do iptr = methods(idx)%maxit, 2, -1
         methods(idx)%ptr(iptr) = methods(idx)%ptr(iptr-1)
       enddo
    endif
    methods(idx)%ptr(1) = history(history_io)%label

    if (pub_on_root .and. pub_md_output_detail>NORMAL) then
       write(stdout,'(a,a)') 'Electronic history : update composition method ', label
    endif

    ! Update the size and counter variables
    methods(idx)%size = min(methods(idx)%size+1,methods(idx)%maxit)
    methods(idx)%ngwfs_count = methods(idx)%ngwfs_count + 1
    methods(idx)%dkn_count = methods(idx)%dkn_count + 1

    ! If required, reset the method
    dkn_reset = .false.
    ngwfs_reset = .false.
    if (methods(idx)%ngwfs_count == methods(idx)%ngwfs_mix_reset) then
       ngwfs_reset = .true.
    endif
    if (methods(idx)%dkn_count == methods(idx)%dkn_mix_reset) then
       dkn_reset = .true.
    endif
    call elec_history_reset_method(ngwfs_reset,dkn_reset,label)

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Leaving elec_history_update_method'

    return

  end subroutine elec_history_update_method

!============================================================================!
!============================================================================!
!============================================================================!

  subroutine elec_history_reset()

    !========================================================================!
    ! This subroutine fully reset the history by resetting the pointers and  !
    ! counters associated with all the composition methods                   !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 27/05/2013                            !
    !========================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout

    implicit none

    ! Local variables
    integer :: idx

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Entering elec_history_reset'

    ! Reset the composition methods
    if (pub_on_root .and. pub_md_output_detail>=BRIEF) &
       write(stdout,'(a)') 'Electronic history : reset all methods'

    if (allocated(methods)) then
       if (methods_max.gt.0) then
          do idx = 1, methods_num
             call elec_history_reset_method(.true.,.true.,methods(idx)%name)
          enddo
       endif
    endif

    ! Reset the electronic history
    if (pub_on_root .and. pub_md_output_detail>=BRIEF) &
       write(stdout,'(a)') 'Electronic history : remove all entries'

    if (allocated(history)) then
       if (history_max.gt.0) then
          do idx = 1, history_size
             call elec_history_remove_entry(idx)
          enddo
       endif
    endif

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving elec_history_reset'

    return

  end subroutine elec_history_reset

!============================================================================!
!============================================================================!
!============================================================================!

  subroutine elec_history_reset_method(ngwfs,dkn,label)

    !========================================================================!
    ! This subroutine resets the pointers and counters associated with the   !
    ! composition method                                                     !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 27/03/2013                            !
    !========================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout
    use utils, only: utils_assert

    implicit none

    ! Parameters
    logical, intent(in)           :: ngwfs
    logical, intent(in)           :: dkn
    character(len=*), intent(in)  :: label

    ! Local variables
    integer :: idx

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Entering elec_history_reset_method'

    ! Preliminary checks
    call utils_assert(allocated(methods),'Error in elec_history_reset_method: &
          &no method defined!')

    ! Find the index of the mixing method
    idx = methods_index(label)
    call utils_assert(idx <= methods_max,'Error in elec_history_reset_method: &
          &method idx is out of range!')

    if (methods(idx)%maxit == 0) return

    if (ngwfs) then
       methods(idx)%ngwfs_count = 0
       if (pub_on_root .and. pub_md_output_detail >= BRIEF) &
          write(stdout,'(a,i2)') 'Electronic history : &
          &reset ngwfs in method ',idx
    endif

    if (dkn) then
       methods(idx)%dkn_count = 0
       if (pub_on_root .and. pub_md_output_detail >= BRIEF) &
          write(stdout,'(a,i2)') 'Electronic history : &
          &reset dkn in method ',idx
    endif

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Leaving elec_history_reset_method'

    return

  end subroutine elec_history_reset_method

!============================================================================!
!============================================================================!
!============================================================================!


  function methods_index(name)

    !========================================================================!
    ! If a method with corresponding name is found in methods, then          !
    ! methods_index return its index, otherwise, methods_index return the    !
    ! index of the first empty slot.                                         !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 27/03/2013                            !
    !========================================================================!

    use constants, only: stdout
    use utils, only: utils_assert

    implicit none

    ! Parameters
    character(len=*), intent(in), optional  :: name
    integer    :: methods_index

    ! Local variable
    integer    :: idx
    logical    :: found

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering methods_index'

    if (present(name)) then
       found = .false.
       do idx = 1, methods_num
          if (methods(idx)%name == name) then
             methods_index = idx
             found = .true.
          endif
       enddo
       call utils_assert(found,'Error in methods_index : &
          &not entry corresponding to label!')
    else
       call utils_assert(methods_num < methods_max,'Error in methods_index : &
             &index is out of range!')
       methods_index = methods_num + 1
    endif

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving methods_index'

  end function methods_index

!============================================================================!
!============================================================================!
!============================================================================!


  function history_index(label)

    !========================================================================!
    ! Returns the index of a given history entry.                            !
    !                                                                        !
    ! - If an entry with corresponding label is found in history, its index  !
    !   is returned,                                                         !
    ! - Else if an isolated entry is found (i.e. an entry with no pointer    !
    !   associated), its index is returned,                                  !
    ! - Else, the index of the oldest item is returned                       !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 27/03/2013                            !
    !========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout
    use utils, only: utils_assert

    implicit none

    ! Parameters
    integer, intent(in), optional  :: label
    integer    :: history_index

    ! Local variable
    integer    :: oldest_label
    integer    :: idx, im, iptr
    logical    :: found, inext

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering history_index'

    if (present(label)) then
       if (pub_on_root .and. pub_md_output_detail > NORMAL) &
          write(stdout,'(a,i4,a)') &
          'Electronic history : look for label ',label,'...'
       found = .false.
       do idx = 1, history_size
          if (history(idx)%label == label) then
             history_index = idx
             found = .true.
          endif
       enddo

       call utils_assert(found,'Error in history_index : &
             &not entry corresponding to label!')

    else
       found = .false.
       if (history_size .lt. history_max) then
          history_index = history_size + 1
          found = .true.
          if (pub_on_root .and. pub_md_output_detail > NORMAL) &
             write(stdout,'(a,i2,a,i2,a)') &
             'Electronic history : history_size =',history_size,&
             '< history_max (',history_max,')'

       else
          do idx = 1, history_max
             if (idx == history_io) cycle
             inext = .false.
             do im = 1, methods_num
                do iptr = 1, methods(im)%maxit
                   if (history(idx)%label == methods(im)%ptr(iptr)) then
                      inext = .true.
                      exit
                   endif
                enddo
                if (inext) exit
             enddo
             if (.not. inext) then
                history_index = idx
                found = .true.
                if (pub_on_root .and. pub_md_output_detail > NORMAL) &
                   write(stdout,'(a,i2)') 'Electronic history : &
                   &no ptr associated with idx ',history_index
                exit
             endif
          enddo
       endif

       if (.not. found) then
          oldest_label = history(1)%label
          history_index = 1
          do idx = 1, history_max
             if (history(idx)%label < oldest_label) then
                oldest_label = history(idx)%label
                history_index = idx
             endif
          enddo
          if (pub_on_root .and. pub_md_output_detail > NORMAL) &
             write(stdout,'(a,i2,a,i2,a)') 'Electronic history : &
             &oldest idx is ',history_index,&
             ' (label =',history(history_index)%label ,')'
       endif

       if (pub_on_root .and. pub_md_output_detail > NORMAL) &
          write(stdout,'(a,i2,a)') &
          'Electronic history : return idx ',history_index,')!'

    endif

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving history_index'

  end function history_index

!============================================================================!
!============================================================================!
!============================================================================!

  subroutine elec_history_create_storage(nsub)

    !========================================================================!
    ! This subroutine allocates space to store electronic degrees of freedom !
    ! in memory.                                                             !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 27/03/2013                            !
    !========================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout
    use utils, only: utils_alloc_check, utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: nsub

    ! Local variables
    integer :: idx, ierr

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering elec_history_create_storage'

    ! Check that history has not been allocated yet
    call utils_assert(.not. allocated(history),'Error in elec_history_create_storage : &
          &attempt to allocate history twice!')

    ! Determine the number of steps to be kept in memory
    history_max = 0
    if (methods_num .ge. 1) then
       do idx = 1, methods_num
          history_max = history_max + methods(idx)%maxit
       enddo
    endif
    if (history_max .ge. 1) history_max = history_max + 1

    ! Allocate history
    if (history_max .gt. 1) then
       allocate(history(history_max),stat=ierr)
       call utils_alloc_check('elec_history_create_storage','history',ierr)
    endif

    ! jcap: allocate ngwfs in history
    do idx = 1, history_max
       allocate(history(idx)%ngwfs(nsub),stat=ierr)
       call utils_alloc_check('elec_history_create_storage','history(idx)%ngwfs',ierr)
    end do

    ! Report info on stdout
    if (methods_num .ge. 1 .and. pub_on_root) then
       do idx = 1, methods_num
          write(stdout,*)
          write(stdout,'(a,i2.2,2a)') 'Electronic history : ', &
               methods(idx)%maxit,' steps required for ', &
               methods(idx)%name
       enddo
       write(stdout,'(a,i2.2)') &
          'Electronic history: Total number of steps to be kept in memory = ', history_max
    endif

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving elec_history_create_storage'

    return

  end subroutine elec_history_create_storage

!============================================================================!
!============================================================================!
!============================================================================!

  subroutine elec_history_destroy_storage()

    !========================================================================!
    ! This subroutine allocates space to store electronic degrees of freedom !
    ! in memory.                                                             !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 27/03/2013                            !
    !========================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout
    use utils, only: utils_dealloc_check

    implicit none

    ! Local variables
    integer :: idx, ierr

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering elec_history_destroy_storage'

    if (allocated(history)) then

       if (history_max .gt. 0) then
          do idx = 1, history_max
             call elec_history_remove_entry(idx)
          enddo
       endif
       deallocate(history,stat=ierr)
       call utils_dealloc_check('elec_history_destroy_storage','history',ierr)

       if (pub_on_root .and. pub_md_output_detail>BRIEF) &
          write(stdout,'(a,i2.2)') &
          'Electronic history : clear storage (no item left) '

    endif

    ! Update history parameters
    history_max = 0
    history_size = 0
    history_count = 0
    history_io = 0

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving elec_history_destroy_storage'

    return

  end subroutine elec_history_destroy_storage

!============================================================================!
!============================================================================!
!============================================================================!

  subroutine elec_history_create_entry(coordinates,label)

    !========================================================================!
    ! This subroutine creates a new item in history.                         !
    ! If history_size is smaller than history_max, then the new item is      !
    ! placed in a free slot, otherwise, the new item replaces the oldest     !
    ! item in history.                                                       !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 30/03/2013                            !
    ! Explicit references to par removed by Robert Charlton, 13/12/2016.     !
    !========================================================================!

    use comms, only : pub_on_root
    use constants, only: DP, stdout
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(in)                :: coordinates(:,:) ! size (3,par%nat)
    integer, optional, intent(in)            :: label

    ! Local Variables
    integer    :: new_label, nsub
    integer    :: ii, idx, ierr

    !======================================================================!

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Entering elec_history_create_entry'

    ! If history_max == 0, then no item can be kept in memory
    if (history_max .le. 1) return

    ! Determine where to store the new item,
    ! if required, erase old item in order to free memory space
    if (history_size .lt. history_max) then
       history_size = history_size + 1
       idx = history_size
    else
       idx = history_index()
       nsub = size(history(idx)%ngwfs)  ! jcap: this shouldn't change
       call elec_history_remove_entry(idx)
       ! jcap: reallocate arrays that depend on subsystems (just ngwfs)
       allocate(history(idx)%ngwfs(nsub),stat=ierr)
       call utils_alloc_check('elec_history_create_entry','history(idx)%ngwfs',ierr)
    endif

    ! Determine label given to new item
    if (present(label)) then
       history(idx)%label = label
    else
       new_label = 1
       do ii = 1, history_size
          if (history(ii)%label .ge. new_label) then
             new_label = history(ii)%label + 1
          endif
       enddo
       history(idx)%label = new_label
    endif

    ! Update coordinates
    ! rc2013: extract the length from the coordinates array instead of hardcoding it
    allocate(history(idx)%coord(3,size(coordinates,dim=2)),stat=ierr)
    call utils_alloc_check('elec_history_create_entry','history(idx)%coord',ierr)
    history(idx)%coord(:,:) = coordinates(:,:)

    ! Update counter
    history_count = history_count + 1

    ! Update io buffer
    history_io = idx

    if (pub_on_root .and. pub_md_output_detail>= NORMAL ) &
       write(stdout,'(a,i2,a,i0,a)') &
       'Electronic history : new entry (idx = ',idx,&
       ') (label = ',history(idx)%label,')'

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Leaving elec_history_create_entry'

    return

  end subroutine elec_history_create_entry

!============================================================================!
!============================================================================!
!============================================================================!

  subroutine elec_history_remove_entry(idx)

    !========================================================================!
    ! This subroutine erases an edf item from history.                       !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 30/03/2013                            !
    ! Modified for embedding by Joseph Prentice, June 2018                   !
    !========================================================================!

    use datatypes, only: data_functions_dealloc, data_fftbox_dealloc
    use comms, only : pub_on_root
    use constants, only: stdout
    use utils, only: utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    integer, intent(in)            :: idx

    ! Local Variables
    integer    :: ierr, ii, isub

    !======================================================================!

    If (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Entering elec_history_remove_entry'

    call utils_assert(idx <= history_max,'Error in elec_history_remove_entry : &
          &idx is out of range!')

    ! Free memory space associated with the item
    if (allocated(history(idx)%coord)) then
       deallocate(history(idx)%coord,stat=ierr)
       call utils_dealloc_check('elec_history_remove_entry','history(idx)%coord',ierr)
    endif

    if (allocated(history(idx)%ngwfs)) then
       do isub = 1,size(history(idx)%ngwfs)
          if (allocated(history(idx)%ngwfs(isub)%procs)) then
             deallocate(history(idx)%ngwfs(isub)%procs,stat=ierr)
             call utils_dealloc_check('elec_history_remove_entry',&
                  'history(idx)%procs(isub)',ierr)
          endif
          if (allocated(history(idx)%ngwfs(isub)%orig)) then
             deallocate(history(idx)%ngwfs(isub)%orig,stat=ierr)
             call utils_dealloc_check('elec_history_remove_entry',&
                  'history(idx)%ngwfs(isub)%orig',ierr)
          endif
          if (allocated(history(idx)%ngwfs(isub)%init)) then
             do ii = lbound(history(idx)%ngwfs(isub)%init, 1), &
                  ubound(history(idx)%ngwfs(isub)%init, 1)
                call data_fftbox_dealloc(history(idx)%ngwfs(isub)%init(ii))
             end do
             deallocate(history(idx)%ngwfs(isub)%init, stat=ierr)
             call utils_dealloc_check('elec_history_remove_entry',&
                  'history(idx)%ngwfs(isub)%init',ierr)
          endif
          if (allocated(history(idx)%ngwfs(isub)%scf)) then
             do ii = lbound(history(idx)%ngwfs(isub)%scf, 1), &
                  ubound(history(idx)%ngwfs(isub)%scf, 1)
                call data_fftbox_dealloc(history(idx)%ngwfs(isub)%scf(ii))
             end do
             deallocate(history(idx)%ngwfs(isub)%scf,stat=ierr)
             call utils_dealloc_check('elec_history_remove_entry',&
                  'history(idx)%ngwfs(isub)%scf',ierr)
          endif
       end do
       deallocate(history(idx)%ngwfs,stat=ierr)
       call utils_dealloc_check('elec_history_remove_entry','history(idx)%ngwfs',ierr)
    endif
    if (allocated(history(idx)%dkn%idx)) then
       deallocate(history(idx)%dkn%idx,stat=ierr)
       call utils_dealloc_check('elec_history_remove_entry','history(idx)%dkn%idx',ierr)
    endif
    if (allocated(history(idx)%dkn%idxinit)) then
       deallocate(history(idx)%dkn%idxinit,stat=ierr)
       call utils_dealloc_check('elec_history_remove_entry','history(idx)%dkn%idxinit',ierr)
    endif
    if (allocated(history(idx)%dkn%idxvel)) then
       deallocate(history(idx)%dkn%idxvel,stat=ierr)
       call utils_dealloc_check('elec_history_remove_entry','history(idx)%dkn%idxvel',ierr)
    endif

    if (allocated(history(idx)%dkn%scf)) then
       do ii = lbound(history(idx)%dkn%scf, 1), ubound(history(idx)%dkn%scf, 1)
          call data_functions_dealloc(history(idx)%dkn%scf(ii))
       end do
       deallocate(history(idx)%dkn%scf, stat=ierr)
       call utils_dealloc_check('elec_history_remove_entry','history(idx)%dkn%scf',ierr)
    end if

    if (allocated(history(idx)%dkn%init)) then
       do ii = lbound(history(idx)%dkn%init, 1), ubound(history(idx)%dkn%init, 1)
          call data_functions_dealloc(history(idx)%dkn%init(ii))
       end do
       deallocate(history(idx)%dkn%init, stat=ierr)
       call utils_dealloc_check('elec_history_remove_entry','history(idx)%dkn%init',ierr)
    end if

    if (allocated(history(idx)%dkn%vel)) then
       do ii = lbound(history(idx)%dkn%vel, 1), ubound(history(idx)%dkn%vel, 1)
          call data_functions_dealloc(history(idx)%dkn%vel(ii))
       end do
       deallocate(history(idx)%dkn%vel, stat=ierr)
       call utils_dealloc_check('elec_history_remove_entry','history(idx)%dkn%vel',ierr)
    end if

    if (pub_on_root .and. pub_md_output_detail>=NORMAL) &
       write(stdout,'(a,i2,a,i8,a)') &
       'Electronic history : remove entry (idx = ',idx, &
       ') (label = ',history(idx)%label,')'

    history(idx)%label = 0

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Leaving elec_history_remove_entry'

    return

  end subroutine elec_history_remove_entry

!=============================================================================!
!=============================================================================!
!=============================================================================!

  subroutine elec_history_store_dkn(denskern,denskern_type,label)

    !======================================================================!
    ! This subroutine stores the density kernel in history. Unsegemented   !
    ! format is used so as not to be affected by changes in atomic         !
    ! positions.                                                           !
    !----------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois 30/05/2013                             !
    ! Modified by Andrea Greco in May 2015 to allow use of complex matrices!
    ! Modified for embedding by Joseph Prentice, June 2018                 !
    !======================================================================!


    ! agrecocmplx
    use comms, only : pub_on_root
    use constants, only: stdout
    use datatypes, only: data_functions_alloc, data_functions_dealloc, &
         data_check_unalloc
    use dense, only: dense_convert
    use rundat, only: md_write_out
    use sparse_embed, only: SPAM3_EMBED_ARRAY
    use sparse, only: SPAM3, sparse_convert_unsegment_real, &
         sparse_convert_unsegment_complex, sparse_create, sparse_destroy
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Argument
    type(SPAM3_EMBED_ARRAY), intent(in) :: denskern
    character(len=*), intent(in) :: denskern_type
    integer, optional, intent(in) :: label

    ! Local variables
    integer  :: ii
    integer  :: idxlen, dklen
    integer  :: ierr
    integer  :: ik, is
    logical  :: loc_iscmplx
    logical  :: loc_allocated
    type(SPAM3) :: tmp_kern(denskern%num_spins)

    !======================================================================!

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Entering elec_history_store_dkn'

    ! Return if history has not been allocated
    if (.not. allocated(history)) return

    ! Look for the index of the storage entry
    if (present(label)) then
       ii = history_index(label=label)
    else
       ii = history_io
    endif
    call utils_assert(ii <= history_max,'Error in elec_history_store_dkn : &
          &index is out of range!')

    ! Return if the corresponding entry has not been initialized
    if (history(ii)%label == 0) return

    ! jcap: initialise temporary SPAM3 matrix
    do is = 1, denskern%num_spins
       ! jcap: HACK!
       tmp_kern(is)%structure='K'
       call sparse_create(tmp_kern(is))
    end do

    loc_iscmplx = denskern%m(1,1)%p%iscmplx

    ! Allocate the storage arrays if required,
    ! and store the density kernel in unsegmented format
    if (denskern_type == 'init') then

       ! jcap: this is a very nasty hack, but the
       ! sparse_convert_segment routines can't deal with non-square
       ! matrices, so it can't be done block-by-block within the
       ! SPAM3_EMBED type. Instead, this routine converts the
       ! SPAM3_EMBED matrices to SPAM3 matrices via a DEM matrix
       do is = 1, denskern%num_spins
          call dense_convert(tmp_kern(is),denskern%m(is,1))
       end do
       if (loc_iscmplx) then
          call sparse_convert_unsegment_complex(tmp_kern, idxlen, dklen)
       else
          call sparse_convert_unsegment_real(tmp_kern, idxlen, dklen)
       end if

       if (.not. allocated(history(ii)%dkn%idxinit)) then
          allocate(history(ii)%dkn%idxinit(idxlen,denskern%num_kpoints),stat=ierr)
          call utils_alloc_check('elec_history_store_dkn','history(ii)%dkn%idxinit',ierr)
       endif

       if (.not. allocated(history(ii)%dkn%init)) then
          allocate(history(ii)%dkn%init(denskern%num_kpoints), stat=ierr)
          call utils_alloc_check('elec_history_store_dkn', 'init', ierr)
          do ik = 1, denskern%num_kpoints
             call data_functions_alloc(history(ii)%dkn%init(ik), dklen, &
                  iscmplx=loc_iscmplx)
          end do
       else
          do ik = 1, denskern%num_kpoints
             if (data_check_unalloc(history(ii)%dkn%init(ik))) then
                call data_functions_alloc(history(ii)%dkn%init(ik), dklen, &
                     iscmplx=loc_iscmplx)
             end if
          end do
       end if

       do ik = 1, denskern%num_kpoints
          ! jcap: nasty hack from above again
          do is = 1, denskern%num_spins
             call dense_convert(tmp_kern(is),denskern%m(is,ik))
          end do
          if (loc_iscmplx) then
             call sparse_convert_unsegment_complex(tmp_kern, idxlen, &
                  dklen, history(ii)%dkn%idxinit(:,ik), &
                  history(ii)%dkn%init(ik)%z)
          else
             call sparse_convert_unsegment_real(tmp_kern, idxlen, &
                  dklen, history(ii)%dkn%idxinit(:,ik), &
                  history(ii)%dkn%init(ik)%d)
          end if

       end do

       if (pub_on_root .and. pub_md_output_detail>=normal .and. &
           (md_write_out)) &
          write(stdout,'(a,i2,a,i0,a)') &
          'Electronic history : store auxiliary dkn in entry (idx = ',&
          ii,') (label = ',history(ii)%label,')'

    else if (denskern_type == 'scf') then

       !jme: KPOINTS_DANGER: these calls should fail (arrays not allocated yet)
       ! jcap: nasty hack from above again
       do is = 1, denskern%num_spins
          call dense_convert(tmp_kern(is),denskern%m(is,1))
       end do
       if (loc_iscmplx) then
          call sparse_convert_unsegment_complex(tmp_kern, idxlen, dklen)
       else
          call sparse_convert_unsegment_real(tmp_kern, idxlen, dklen)
       end if

       if (.not. allocated(history(ii)%dkn%idx)) then
          allocate(history(ii)%dkn%idx(idxlen,denskern%num_kpoints),stat=ierr)
          call utils_alloc_check('elec_history_store_dkn',&
               'history(ii)%dkn%idx',ierr)
       endif

       if (.not. allocated(history(ii)%dkn%scf)) then
          allocate(history(ii)%dkn%scf(denskern%num_kpoints), stat=ierr)
          call utils_alloc_check('elec_history_store_dkn',&
               'history(ii)%dkn%scf',ierr)
          do ik = 1, denskern%num_kpoints
             call data_functions_alloc(history(ii)%dkn%scf(ik), dklen, &
                  iscmplx=loc_iscmplx)
          end do
       else
          do ik = 1, denskern%num_kpoints
             if (data_check_unalloc(history(ii)%dkn%scf(ik))) then
                call data_functions_alloc(history(ii)%dkn%scf(ik), dklen, &
                     iscmplx=loc_iscmplx)
             end if
          end do
       end if

       do ik = 1, denskern%num_kpoints
          ! jcap: nasty hack from above again
          do is = 1, denskern%num_spins
             call dense_convert(tmp_kern(is),denskern%m(is,ik))
          end do
          if (loc_iscmplx) then
             call sparse_convert_unsegment_complex(tmp_kern, &
                  idxlen, dklen, history(ii)%dkn%idx(:,ik), &
                  history(ii)%dkn%scf(ik)%z)
          else
             call sparse_convert_unsegment_real(tmp_kern, &
                  idxlen, dklen, history(ii)%dkn%idx(:,ik), &
                  history(ii)%dkn%scf(ik)%d)
          end if

       end do

       if (pub_on_root .and. pub_md_output_detail>=NORMAL .and. &
           (md_write_out)) &
          write(stdout,'(a,i2,a,i0,a)') &
          'Electronic history : store dkn in entry (idx = ',&
          ii,'), (label = ',history(ii)%label,')'

    else if (denskern_type == 'vel') then

       !jme: KPOINTS_DANGER: these calls should fail (arrays not allocated yet)
       ! jcap: nasty hack from above again
       do is = 1, denskern%num_spins
          call dense_convert(tmp_kern(is),denskern%m(is,1))
       end do
       if (loc_iscmplx) then
          call sparse_convert_unsegment_complex(tmp_kern, idxlen, dklen)
       else
          call sparse_convert_unsegment_real(tmp_kern, idxlen, dklen)
       end if

       if (.not. allocated(history(ii)%dkn%idxvel)) then
          allocate(history(ii)%dkn%idxvel(idxlen,denskern%num_kpoints),stat=ierr)
          call utils_alloc_check('elec_history_store_dkn','history(ii)%dkn%idxvel',ierr)
       endif

       if (.not. allocated(history(ii)%dkn%vel)) then
          allocate(history(ii)%dkn%vel(denskern%num_kpoints), stat=ierr)
          call utils_alloc_check('elec_history_store_dkn','history(ii)%dkn%vel',ierr)
          do ik = 1, denskern%num_kpoints
             call data_functions_alloc(history(ii)%dkn%vel(ik), dklen, &
                  iscmplx=loc_iscmplx)
          end do
       else
          do ik = 1, denskern%num_kpoints
             if (data_check_unalloc(history(ii)%dkn%vel(ik))) then
                call data_functions_alloc(history(ii)%dkn%vel(ik), dklen, &
                     iscmplx=loc_iscmplx)
             end if
          end do
       end if

       do ik = 1, denskern%num_kpoints
          ! jcap: nasty hack from above again
          do is = 1, denskern%num_spins
             call dense_convert(tmp_kern(is),denskern%m(is,ik))
          end do
          if (loc_iscmplx) then
             call sparse_convert_unsegment_complex(tmp_kern, &
                  idxlen, dklen, history(ii)%dkn%idxvel(:,ik), &
                  history(ii)%dkn%vel(ik)%z)
          else
             call sparse_convert_unsegment_real(tmp_kern, &
                  idxlen, dklen, history(ii)%dkn%idxvel(:,ik), &
                  history(ii)%dkn%vel(ik)%d)
          end if

       end do

       if (pub_on_root .and. pub_md_output_detail>=NORMAL .and. &
           (md_write_out)) &
          write(stdout,'(a,i2,a,i0,a)') &
          'Electronic history : store vel dkn in entry (idx = ',&
          ii,'), (label = ',history(ii)%label,')'

    endif

    ! jcap: destroy temporary SPAM3 matrix
    do is = 1, denskern%num_spins
       call sparse_destroy(tmp_kern(is))
    end do

    if (pub_debug_on_root) write(stdout,'(a)') &
        'DEBUG: Leaving elec_history_store_dkn'

  end subroutine elec_history_store_dkn

!============================================================================!
!============================================================================!
!============================================================================!

  subroutine elec_history_read_dkn(denskern, denskern_type, label)

    !======================================================================!
    ! This subroutine returns a density kernel stored in history.          !
    ! This operation is not affected by changes in atomic positions.       !
    !----------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois 30/05/2013                             !
    ! Modified for embedding by Joseph Prentice, June 2018                 !
    !======================================================================!


    use comms, only : pub_on_root
    use constants, only: stdout
    use dense, only: dense_convert
    use rundat, only: md_write_out
    use sparse_embed, only: SPAM3_EMBED_ARRAY
    use sparse, only: SPAM3, sparse_convert_segment_real, &
         sparse_convert_segment_complex, sparse_create, sparse_destroy
    use utils, only: utils_alloc_check, utils_assert

    implicit none

    ! Argument
    type(SPAM3_EMBED_ARRAY), intent(inout) :: denskern
    character(len=*), intent(in) :: denskern_type
    integer, intent(in) :: label

    ! Internal variables
    integer  :: ii, ik, is
    type(SPAM3) :: tmp_kern(denskern%num_spins)

    !======================================================================!

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Entering elec_history_read_dkn'

    ! Exit if history has not been allocated
    call utils_assert(allocated(history),'Error in elec_history_read_dkn : &
          &history has not been allocated!')

    ii = history_index(label=label)

    ! Exit if index is out of range
    call utils_assert(ii <= history_max,'Error in elec_history_read_dkn : &
          &index is out of range!')

    ! Exit if the corresponding entry has not been initialized
    call utils_assert(history(ii)%label /= 0, &
          'Error in elec_history_read_dkn : &
          &corresponding entry has not been initialized!')

    ! jcap: initialise temporary SPAM3 matrix
    do is = 1, denskern%num_spins
       tmp_kern(is)%structure='K'
       call sparse_create(tmp_kern(is))
    end do

    ! Read the density kernel in unsegmented format
    if (denskern_type == 'init') then

       if (pub_on_root .and. pub_md_output_detail>=NORMAL .and. &
           (md_write_out)) &
          write(stdout,'(a,i2,a,i0,a)') 'Electronic history : &
          &read auxiliary dkn from entry (idx = ',ii,')&
          &(label = ',history(ii)%label,')'

       do ik = 1, denskern%num_kpoints
          if (history(ii)%dkn%init(ik)%iscmplx) then
             call sparse_convert_segment_complex(tmp_kern, &
                  history(ii)%dkn%idx(:,ik), history(ii)%dkn%init(ik)%z)
          else
             call sparse_convert_segment_real(tmp_kern, &
                  history(ii)%dkn%idx(:,ik), history(ii)%dkn%init(ik)%d)
          end if
          ! jcap: this is a very nasty hack, but the
          ! sparse_convert_segment routines can't deal with non-square
          ! matrices, so it can't be done block-by-block within the
          ! SPAM3_EMBED type. Instead, this routine converts the
          ! SPAM3 matrices to SPAM3_EMBED matrices via a DEM matrix
          do is = 1, denskern%num_spins
             call dense_convert(denskern%m(is,ik),tmp_kern(is))
          end do
       end do

    else if (denskern_type == 'scf') then

       if (pub_on_root .and. pub_md_output_detail>=NORMAL .and. &
           (md_write_out)) &
           write(stdout,'(a,i2,a,i0,a)') 'Electronic history : &
           &read dkn from entry (idx = ',ii,') (label = ', &
           history(ii)%label,')'

       do ik = 1, denskern%num_kpoints
          if (history(ii)%dkn%scf(ik)%iscmplx) then
             call sparse_convert_segment_complex(tmp_kern, &
                  history(ii)%dkn%idx(:,ik), history(ii)%dkn%scf(ik)%z)
          else
             call sparse_convert_segment_real(tmp_kern, &
                  history(ii)%dkn%idx(:,ik), history(ii)%dkn%scf(ik)%d)
          end if
          ! jcap: nasty hack from above again
          do is = 1, denskern%num_spins
             call dense_convert(denskern%m(is,ik),tmp_kern(is))
          end do
       end do

    else if (denskern_type == 'vel') then

       if (pub_on_root .and. pub_md_output_detail>=NORMAL .and. &
           (md_write_out)) &
           write(stdout,'(a,i2,a,i0,a)') 'Electronic history : &
           &read KS vel from entry (idx = ',ii,') (label = ', &
           history(ii)%label,')'

       do ik = 1, denskern%num_kpoints
          if (history(ii)%dkn%vel(ik)%iscmplx) then
             call sparse_convert_segment_complex(tmp_kern, &
                  history(ii)%dkn%idx(:,ik), history(ii)%dkn%vel(ik)%z)
          else
             call sparse_convert_segment_real(tmp_kern, &
                  history(ii)%dkn%idx(:,ik), history(ii)%dkn%vel(ik)%d)
          end if
          ! jcap: nasty hack from above again
          do is = 1, denskern%num_spins
             call dense_convert(denskern%m(is,ik),tmp_kern(is))
          end do
       end do

    endif

    ! jcap: destroy temporary SPAM3 matrix
    do is = 1, denskern%num_spins
       call sparse_destroy(tmp_kern(is))
    end do

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving elec_history_read_dkn'

  end subroutine elec_history_read_dkn


!============================================================================!
!============================================================================!
!============================================================================!


  subroutine elec_history_compose_dkn(denskern, rep, ngwf_basis, &
       nl_projectors, proj_basis, mdl, cstatus, kpt)

    !======================================================================!
    ! This subroutine build a new density kernel according to              !
    ! mix_dkn_type and mix_dkn_num                                         !
    !----------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois 25/2/2010, modified (02/12/2010).      !
    ! Adjusted by Nicholas D.M. Hine 18/05/2012 to take nl_projectors as   !
    ! an input.                                                            !
    ! Modified for embedding by Joseph Prentice, June 2018                 !
    !======================================================================!

    use comms, only : pub_on_root
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use kernel, only: DKERN, kernel_normalise
    use model_type, only : MODEL
    use ngwf_representation, only: NGWF_REP
    use parallel_strategy, only: PARAL_INFO
    use projectors, only: PROJECTOR_SET
    use rundat, only: lnv_threshold_orig, ngwf_threshold_orig, &
         md_lnv_threshold, md_ngwf_threshold, &
         mermin_threshold_orig, pub_mermin
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Argument
    type(DKERN), intent(inout)     :: denskern
    type(FUNC_BASIS), intent(in)   :: ngwf_basis(:)
    type(FUNC_BASIS), intent(in)   :: proj_basis(:)
    type(MODEL), target, intent(inout)   :: mdl
    type(NGWF_REP), intent(in)     :: rep
    type(PROJECTOR_SET), intent(inout) :: nl_projectors(:)
    integer, intent(out)           :: cstatus
    ! agrecokpt: for kpoint dependence of projectors in internal_projection()
    real(kind=DP), optional, intent(in) :: kpt(3)

    ! Internal variables
    integer       :: label, iom, nsub
    real(kind=DP), save :: tmp_ngwf_or, tmp_lnv_or
    type(PARAL_INFO), pointer :: par

    !======================================================================!

    if (pub_debug_on_root) write(stdout,'(a,i4,a,i4,a,i4,a)') &
        'DEBUG: Entering elec_history_compose_dkn (hio =', history_io, &
        '), (mio =',methods_io,'), (hcount =', history_count,')'

    ! Exit if stores haven't been initialized
    if (.not.allocated(history) .or. .not.allocated(methods) .or. &
        history_io==0 .or. methods_io==0 .or. history_count==0) then
       cstatus = -1
       return
    endif

    nsub = size(ngwf_basis)

    ! Determine composition method
    iom = methods_io
    cstatus = 1

    if (methods(iom)%dkn_count .ge. methods(iom)%dkn_init_num + 1) then
       ngwf_threshold_orig = md_ngwf_threshold
       if (pub_mermin) then
           mermin_threshold_orig = md_lnv_threshold
       else
           lnv_threshold_orig = md_lnv_threshold
       end if
    end if


    !=============================================================!
    ! New entry
    !=============================================================!
    if (methods(iom)%dkn_mix_type == 'NONE' &
            .or.   methods(iom)%dkn_count == 0) then

       cstatus = 1

    !=============================================================!
    ! Initialisation phase
    !=============================================================!
    elseif (methods(iom)%dkn_count .le. methods(iom)%dkn_init_num) then

       call internal_init()

    !=============================================================!
    ! No mixing : new denskern = last denskern stored
    !=============================================================!
    elseif (methods(iom)%dkn_mix_type == 'REUSE') then

       call internal_retrieve_last()

    !=============================================================!
    ! Linear or polynomial extrapolations,
    ! transformation of last denskern stored
    !=============================================================!
    elseif (methods(iom)%dkn_mix_type == 'LINEAR') then
       call internal_linear_xtpol()

    elseif (methods(iom)%dkn_mix_type == 'MULTID') then
       call internal_arias_xtpol()

    elseif (methods(iom)%dkn_mix_type == 'POLY') then
       call internal_poly_xtpol()

    elseif (methods(iom)%dkn_mix_type == 'PROJ') then
       call internal_projection()

    elseif (methods(iom)%dkn_mix_type == 'TENSOR') then
       call internal_christoffel()

    !=============================================================!
    ! Propagation of auxiliary kernel :
    ! algorithm based on A. Niklasson et al,
    !           J. Chem. Phys. 130, 214109 (2009)
    !=============================================================!
    elseif (methods(iom)%dkn_mix_type == 'TRPROP') then
       call internal_tr_propagation()

    elseif (methods(iom)%dkn_mix_type == 'XLD') then
       call internal_xl_dissip_propagation()

    elseif (methods(iom)%dkn_mix_type == 'NAIVE') then
       call internal_xl_berendsen_propagation()

    elseif (methods(iom)%dkn_mix_type == 'XLI') then
       call internal_xl_berendsen_propagation()

    elseif (methods(iom)%dkn_mix_type == 'XLIS') then
       call internal_xl_berendsen_scuseria_propagation()
    endif

    !=============================================================!
    ! Check the normalization of new density kernel
    !=============================================================!

    if (cstatus == 0) then
       call kernel_normalise(denskern, rep%overlap, rep%inv_overlap, rep%n_occ)
    endif

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving elec_history_compose_dkn'

  contains

    !################################################################!

    subroutine internal_init()

       use kernel, only: kernel_workspace_create,kernel_workspace_destroy
       use rundat, only: pub_num_spins, md_aux_rep, PUB_1K
       use smearing_operator, only: smearing_matrix, lowdin_transformation
       use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_product, &
            sparse_embed_create, sparse_embed_destroy, sparse_embed_copy, &
            sparse_embed_array_create, sparse_embed_array_destroy, &
            sparse_embed_array_scale, sparse_embed_array_product, &
            sparse_embed_array_axpy, sparse_embed_array_transpose, &
            sparse_embed_array_copy

       implicit none

       type(SPAM3_EMBED_ARRAY) :: tmp_x, dkn_trans, tmp_dkn
       type(smearing_matrix) :: inv_sqrt_overlap, tmp_overlap
       character(len=2) :: structure_type
       integer :: is, nkpts, nspins

       nkpts = denskern%kern%num_kpoints
       nspins = denskern%kern%num_spins

       if (methods(iom)%dkn_mix_type == 'TRPROP') then

          if( denskern%workspace_created) then
             label = methods(iom)%ptr(1)
             call elec_history_read_dkn(denskern%ks,'scf',label)
             call elec_history_store_dkn(denskern%ks,'init')
          else
             call kernel_workspace_create(denskern,rep%overlap)
             label = methods(iom)%ptr(1)
             call elec_history_read_dkn(denskern%ks,'scf',label)
             call elec_history_store_dkn(denskern%ks,'init')
             call kernel_workspace_destroy(denskern)
          end if

          if (methods(iom)%dkn_init_type == 'NONE') then
             cstatus = 1

          elseif (methods(iom)%dkn_init_type == 'REUSE') then
             cstatus = 0
             if (pub_on_root) write(stdout,'(a)') &
                'Electronic history : new density kernel = &
                &last kernel stored '

          endif

       else if ( methods(iom)%dkn_mix_type == 'XLD'    &
            .or. methods(iom)%dkn_mix_type == 'NAIVE'     &
            .or. methods(iom)%dkn_mix_type == 'XLI' &
            .or. methods(iom)%dkn_mix_type == 'XLIS' ) then

          ! vv: Define structures for SPAM3 and SPAM3_ARRAY matrices
          if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
             structure_type = 'K'
          else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
             structure_type = 'KS'
          end if
          call sparse_embed_array_create(tmp_x,n_kpoints=nkpts,n_spins=nspins, &
               structure=structure_type)

          label = methods(iom)%ptr(1)
          call elec_history_read_dkn(tmp_x,'scf',label)
          call elec_history_store_dkn(tmp_x,'init')

          if (methods(iom)%dkn_init_type == 'NONE') then
             cstatus = 1

             call sparse_embed_array_scale(tmp_x,0.0_DP)
             call elec_history_store_dkn(tmp_x,'vel')

          elseif (methods(iom)%dkn_init_type == 'REUSE') then
             ! vv: Go back to the non-orthogonal representation
             if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
                ! vv: Define structures for SPAM3 and SPAM3_ARRAY matrices
                inv_sqrt_overlap%matrix_type = 3
                tmp_overlap%matrix_type = 3
                inv_sqrt_overlap%dataSPAM3_EMBED%structure = 'K'
                tmp_overlap%dataSPAM3_EMBED%structure = 'S'

                ! vv: Allocate workspace
                call sparse_embed_create(inv_sqrt_overlap%dataSPAM3_EMBED)
                call sparse_embed_create(tmp_overlap%dataSPAM3_EMBED)
                call sparse_embed_copy(tmp_overlap%dataSPAM3_EMBED,rep%overlap)


                call sparse_embed_array_copy(denskern%kern,tmp_x)
                call sparse_embed_array_scale(tmp_x,0.0_DP)

                ! vv: Compute K = S.^(1/2).X.S^(1/2)
                call lowdin_transformation(tmp_overlap,inv_sqrt_overlap)

                do is = 1,pub_num_spins
                   call sparse_embed_product(tmp_x%m(is,PUB_1K),&
                        denskern%kern%m(is,PUB_1K),&
                        inv_sqrt_overlap%dataSPAM3_EMBED)
                end do

                call sparse_embed_array_scale(denskern%kern,0.0_DP)
                call elec_history_store_dkn(denskern%kern,'vel')

                do is = 1,pub_num_spins
                   call sparse_embed_product(denskern%kern%m(is,PUB_1K),&
                        inv_sqrt_overlap%dataSPAM3_EMBED,tmp_x%m(is,PUB_1K))
                end do

                call sparse_embed_destroy(inv_sqrt_overlap%dataSPAM3_EMBED)
                call sparse_embed_destroy(tmp_overlap%dataSPAM3_EMBED)
                call sparse_embed_array_destroy(tmp_x)

             else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
                call sparse_embed_array_create(tmp_dkn,n_kpoints=nkpts,&
                     n_spins=pub_num_spins, structure=structure_type)
                call sparse_embed_array_create(dkn_trans,n_kpoints=nkpts,&
                     n_spins=pub_num_spins, structure=structure_type)
                call sparse_embed_array_copy(denskern%kern,tmp_x)
                ! vv: Compute K = X.S^(-1)
                call sparse_embed_array_product(tmp_dkn,denskern%kern,&
                     rep%inv_overlap)

                ! vv: At this stage K and S do not commute, i.e. [K,S] \= 0.
                !     So we need to symmetrize the density kernel
                call sparse_embed_array_transpose(dkn_trans,tmp_dkn)

                call sparse_embed_array_scale(denskern%kern,0.0_DP)
                call elec_history_store_dkn(denskern%kern,'vel')

                call sparse_embed_array_axpy(denskern%kern,tmp_dkn,1.0_dp)
                call sparse_embed_array_axpy(denskern%kern,dkn_trans,1.0_dp)
                call sparse_embed_array_scale(denskern%kern,0.5_dp)

                call sparse_embed_array_destroy(tmp_dkn)
                call sparse_embed_array_destroy(dkn_trans)
             end if
             cstatus = 0
             if (pub_on_root) write(stdout,'(a)') &
                'Electronic history : new density kernel = &
                &last kernel stored '
          end if

          call sparse_embed_array_destroy(tmp_x)

       else

          if (methods(iom)%dkn_init_type == 'NONE') then
             cstatus = 1

          elseif (methods(iom)%dkn_init_type == 'REUSE') then
             label = methods(iom)%ptr(1)
             call elec_history_read_dkn(denskern%kern,'scf',label)
             cstatus = 0
             if (pub_on_root) write(stdout,'(a)') &
                'Electronic history : new density kernel = &
                &last kernel stored '

          endif
       endif

    end subroutine internal_init

    !################################################################!

    subroutine internal_retrieve_last()

       !=============================================================!
       ! No mixing : new density kernel = last density kernel stored
       !=============================================================!

       implicit none

       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : new density kernel = &
          &last kernel stored '
       label = methods(iom)%ptr(1)
       call elec_history_read_dkn(denskern%kern,'scf',label)
       cstatus = 0

    end subroutine internal_retrieve_last

    !################################################################!

    subroutine internal_linear_xtpol()

       !=============================================================!
       ! Linear extrapolation : dkn(new) = 2*dkn(-1) - dkn(-2)
       !=============================================================!

       use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_array_create, &
            sparse_embed_array_scale, sparse_embed_array_axpy, &
            sparse_embed_array_destroy
       implicit none

       type(SPAM3_EMBED_ARRAY) :: old_denskern

       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : linear extrapolation of &
          &density kernels '

       call sparse_embed_array_create(old_denskern, denskern%kern)
       call sparse_embed_array_scale(denskern%kern, 0.0_DP)

       label = methods(iom)%ptr(1)
       call elec_history_read_dkn(old_denskern, 'scf', label)
       call sparse_embed_array_axpy(denskern%kern, old_denskern, 2.0_dp)

       label = methods(iom)%ptr(2)
       call elec_history_read_dkn(old_denskern, 'scf', label)
       call sparse_embed_array_axpy(denskern%kern, old_denskern, -1.0_dp)
       cstatus = 0

       call sparse_embed_array_destroy(old_denskern)

    end subroutine internal_linear_xtpol

    !################################################################!

    subroutine internal_arias_xtpol()

       !=============================================================!
       ! Linear extrapolation :
       ! algorithm based on Arias et al, PRB 45, 1538 (1992)
       !=============================================================!

       use services, only: services_rms_fit
       use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_array_create, &
            sparse_embed_array_scale, sparse_embed_array_axpy, &
            sparse_embed_array_destroy

       implicit none

       real(kind=DP), allocatable  :: linear_coeffs(:)
       real(kind=DP), allocatable  :: old_coords(:,:,:)
       real(kind=DP), allocatable  :: new_coords(:,:)
       type(SPAM3_EMBED_ARRAY) :: old_denskern
       integer  :: ii, istep
       integer  :: ierr
       integer  :: nat

       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : multi-dimensional extrapolation &
          &of density kernels '

       ! rc2013: get no. of atoms from size of elements array
       nat = size(mdl%elements)

       ! Mixing coefficients
       allocate(linear_coeffs(methods(iom)%dkn_mix_size),stat=ierr)
       call utils_alloc_check('elec_history_compose_dkn','linear_coeffs',ierr)

       ! Set of coordinates to be interpolated
       allocate(new_coords(3,nat),stat=ierr)
       call utils_alloc_check('elec_history_compose_dkn','new_coords',ierr)
       new_coords = history(history_io)%coord

       allocate(old_coords(3,nat,methods(iom)%dkn_mix_size),stat=ierr)
       call utils_alloc_check('elec_history_compose_dkn','old_coords',ierr)
       do istep = 1, methods(iom)%dkn_mix_size
          ! Debug
          ii = history_index(label=methods(iom)%ptr(istep))
          old_coords(:,:,istep) = history(ii)%coord
          !old_coords(:,:,istep) = history(methods(iom)%ptr(istep))%coord
       enddo

       ! Interpolation coefficients
       call services_rms_fit(linear_coeffs,methods(iom)%dkn_mix_size,3,nat,&
          old_coords,new_coords)

       deallocate(old_coords,stat=ierr)
       call utils_dealloc_check('elec_history_compose_dkn','old_coords',ierr)
       deallocate(new_coords,stat=ierr)
       call utils_dealloc_check('elec_history_compose_dkn','new_coords',ierr)

       ! Linear combination of density kernels
       call sparse_embed_array_create(old_denskern, denskern%kern)
       call sparse_embed_array_scale(denskern%kern, 0.0_DP)

       do istep = 1, methods(iom)%dkn_mix_size
          label = methods(iom)%ptr(istep)
          call elec_history_read_dkn(old_denskern, 'scf',label)
          call sparse_embed_array_axpy(denskern%kern, old_denskern, linear_coeffs(istep))
       enddo
       cstatus = 0

       call sparse_embed_array_destroy(old_denskern)
       deallocate(linear_coeffs,stat=ierr)
       call utils_dealloc_check('elec_history_compose_dkn','linear_coeffs',ierr)

    end subroutine internal_arias_xtpol

    !################################################################!

    subroutine internal_poly_xtpol()

       !=============================================================!
       ! Polynomial extrapolation
       ! Modified by Andrea Greco to allow use of complex matrices
       !=============================================================!

       use services, only : services_polynomial_step
       use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_array_create, &
            sparse_embed_array_destroy

       implicit none

       type(SPAM3_EMBED_ARRAY) :: old_denskern
       real(kind=DP), allocatable  :: time_steps(:)
       real(kind=DP), allocatable  :: dmtx_history(:,:,:)
       ! agrecocmplx: new variable for complex matrices
       complex(kind=DP), allocatable :: zmtx_history(:,:,:)
       real(kind=DP) :: xtpol_step
       integer, allocatable :: locnum(:,:,:,:)
       integer  :: istep, locsize
       integer  :: nkpoints, nspins, nkns
       integer  :: is, ik, iks, ierr, isub, jsub

       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : polynomial extrapolation of &
          &density kernels '

       allocate(time_steps(methods(iom)%dkn_mix_size), stat=ierr)
       call utils_alloc_check('elec_history_compose_dkn','time_steps',ierr)

       nkpoints = denskern%kern%num_kpoints
       nspins = denskern%kern%num_spins
       nkns = nkpoints * nspins

       call sparse_embed_array_create(old_denskern, denskern%kern)

       allocate(locnum(nspins,nkpoints,nsub,nsub), stat=ierr)
       call utils_alloc_check('elec_history_compose_dkn','locnum',ierr)

       do ik = 1, nkpoints
          do is = 1, nspins
             ! jcap: Loop over all pairs of regions
             do isub=1, nsub
                do jsub=1, nsub
                   if (denskern%kern%m(is,ik)%m(isub,jsub)%iscmplx) then
                      locnum(is,ik,isub,jsub) = &
                           size(denskern%kern%m(is,ik)%m(isub,jsub)%zmtx, 1)
                   else
                      locnum(is,ik,isub,jsub) = &
                           size(denskern%kern%m(is,ik)%m(isub,jsub)%dmtx, 1)
                   end if
                end do
             end do
          end do
       end do

       locsize = maxval(locnum(:,:,:,:))

       ! agrecocmplx: assume all matrices are either real or complex
       ! jmecmplx: add check
       if (denskern%kern%m(1,1)%p%iscmplx) then
           allocate(zmtx_history(methods(iom)%dkn_mix_size, &
                    locsize, nkns), stat=ierr)
           call utils_alloc_check('elec_history_compose_dkn','zmtx_history',ierr)
           zmtx_history(:,:,:) = 0.0_dp
       else
           allocate(dmtx_history(methods(iom)%dkn_mix_size, &
                    locsize, nkns), stat=ierr)
           call utils_alloc_check('elec_history_compose_dkn','dmtx_history',ierr)
           dmtx_history(:,:,:) = 0.0_dp
       end if

       do istep = 1,methods(iom)%dkn_mix_size
          label = methods(iom)%ptr(istep)
          call elec_history_read_dkn(old_denskern, 'scf', label)
          do ik = 1, nkpoints
             do is = 1, nspins
                iks = is + nspins * (ik-1)
                ! jcap: Loop over regions
                do isub = 1, nsub
                   do jsub = 1, nsub
                      if (old_denskern%m(is,ik)%m(isub,jsub)%iscmplx) then
                         zmtx_history(istep,1:locnum(is,ik,isub,jsub),iks) = &
                              old_denskern%m(is,ik)%m(isub,jsub)%zmtx
                      else
                         dmtx_history(istep,1:locnum(is,ik,isub,jsub),iks) = &
                              old_denskern%m(is,ik)%m(isub,jsub)%dmtx
                      end if
                   end do
                end do
             end do
          end do
          time_steps(istep) = real(istep)
       enddo

       call sparse_embed_array_destroy(old_denskern)

       xtpol_step = real(methods(iom)%dkn_mix_size+1)
       do ik = 1, nkpoints
          do is = 1, nspins
             iks = is + nspins * (ik-1)
             ! jcap: Loop over regions
             do isub = 1, nsub
                do jsub = 1, nsub
                   if (denskern%kern%m(is,ik)%m(isub,jsub)%iscmplx) then
                      call services_polynomial_step(&
                           denskern%kern%m(is,ik)%m(isub,jsub)%zmtx, &
                           methods(iom)%dkn_mix_size, locnum(is,ik,isub,jsub), &
                           zmtx_history(1:methods(iom)%dkn_mix_size,&
                           1:locnum(is,ik,isub,jsub),iks), &
                           time_steps, xtpol_step)
                   else
                      call services_polynomial_step(&
                           denskern%kern%m(is,ik)%m(isub,jsub)%dmtx, &
                           methods(iom)%dkn_mix_size, locnum(is,ik,isub,jsub), &
                           dmtx_history(1:methods(iom)%dkn_mix_size,&
                           1:locnum(is,ik,isub,jsub),iks), &
                           time_steps, xtpol_step)
                   end if
                end do
             end do
          end do
       end do
       cstatus = 0

       !deallocate(dmtx_history, stat=ierr)
       !call utils_dealloc_check('elec_history_compose_dkn','dmtx_history',ierr)
       ! agrecocmplx
       if (allocated(zmtx_history)) then
           deallocate(zmtx_history, stat=ierr)
           call utils_dealloc_check('elec_history_compose_dkn',&
                'zmtx_history',ierr)
       else if (allocated(dmtx_history)) then
           deallocate(dmtx_history, stat=ierr)
           call utils_dealloc_check('elec_history_compose_dkn',&
                'dmtx_history',ierr)
       end if
       deallocate(locnum, stat=ierr)
       call utils_dealloc_check('elec_history_compose_dkn','locnum',ierr)
       deallocate(time_steps, stat=ierr)
       call utils_dealloc_check('elec_history_compose_dkn',&
            'time_steps',ierr)

    end subroutine internal_poly_xtpol

    !################################################################!

    subroutine internal_projection()

       !=============================================================!
       ! Project last kernel onto new NGWFs                          !
       ! Modified by Andrea Greco to allow use of complex matrices   !
       !=============================================================!

       use datatypes, only: FUNCTIONS, data_functions_alloc, &
            data_functions_dealloc
       use function_ops, only: function_ops_brappd_ketppd
       use projectors, only: projectors_func_ovlp_box
       use rundat, only: pub_aug, pub_realspace_projectors, &
            pub_kpoint_method
       use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
            sparse_embed_create, sparse_embed_destroy, &
            sparse_embed_array_create, sparse_embed_array_destroy
       use utils, only: utils_abort, utils_assert
       ! agrecokpt: needed for testing in case of single k-point
       use geometry, only: POINT

       implicit none

       ! Local Variables
       type(FUNCTIONS) :: old_ngwfs_on_grid(nsub)
       type(SPAM3_EMBED_ARRAY) :: old_denskern
       type(SPAM3_EMBED) :: old_sp_overlap
       integer  :: isub, jsub
       ! agrecokpt
       type(POINT) :: loc_kpt

       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : projection of density kernel &
          &onto new NGWFs'

       ! agrecokpt: currently only KP method is being implemented
       call utils_assert(pub_kpoint_method == 'KP', &
            'Subroutine internal_projection in elec_history_compose_dkn &
            &currently supports only KP method for BZ sampling')

       ! agrecokpt: single k-point for now, mainly for testing,
       ! in the future take k-points directly from denskern and
       ! NGWFs
       if (present(kpt)) then
          loc_kpt%x = kpt(1)
          loc_kpt%y = kpt(2)
          loc_kpt%z = kpt(3)
       else
          loc_kpt%x = 0.0_DP
          loc_kpt%y = 0.0_DP
          loc_kpt%z = 0.0_DP
       end if

       ! Retrieve last kernel from store
       call sparse_embed_array_create(old_denskern, denskern%kern)
       label = methods(iom)%ptr(1)
       call elec_history_read_dkn(old_denskern, 'scf',label)

       ! Retrieve last ngwfs from store
       label = methods(iom)%ptr(1)
       ! jcap: Loop over regions
       do isub = 1, nsub
          par => mdl%regions(isub)%par
          call data_functions_alloc(old_ngwfs_on_grid(isub), &
               ngwf_basis(isub)%size_on_grid, &
               iscmplx=denskern%kern%m(1,1)%m(isub,isub)%iscmplx)
          ! agrecokpt: need to modify this to read NGWFs for each k-point
          call elec_history_read_ngwfs(old_ngwfs_on_grid(isub),&
               ngwf_basis(isub),mdl%cell,mdl%fftbox,&
               mdl%regions(isub)%elements,isub,'scf',label,par=par)
       end do

       call sparse_embed_create(old_sp_overlap,rep%sp_overlap)

       ! If required, calculate overlap of old NGWFs with projectors
       ! agrecokpt: need to include kpoint dependence here for KP method
       ! as well? currently coded in, since kpoint dependence is already
       ! included if using pub_realspace_projectors
       if (pub_aug) then
          ! jcap: Loop over all pairs of regions
          do isub = 1, nsub
             do jsub = 1,nsub
                if (.not.pub_realspace_projectors) then
                   call projectors_func_ovlp_box(old_sp_overlap%m(isub,jsub), &
                        old_ngwfs_on_grid(isub),ngwf_basis(isub),&
                        proj_basis(jsub),nl_projectors(jsub),mdl%fftbox,&
                        mdl%cell,kshift=loc_kpt)
                   ! agrecokpt: kpoint dependence included in definition of projs_on_grid
                else
                   call function_ops_brappd_ketppd(old_sp_overlap%m(isub,jsub),&
                        old_ngwfs_on_grid(isub),ngwf_basis(isub),&
                        nl_projectors(jsub)%projs_on_grid,proj_basis(jsub),&
                        mdl%cell)
                end if
             end do
          end do
       endif

       ! Apply basis transformation to density kernel
       ! agrecokpt: need to make this compatible with new structures
       !call kernel_basis_transform(ngwf_basis, rep%ngwfs_on_grid, &
       !     old_ngwfs_on_grid, denskern%kern, old_denskern, rep%inv_overlap, &
       !     rep%sp_overlap, old_sp_overlap)
       call utils_abort('Error in internal_projection &
            &(elec_history_compose_dkn): Method not supported')
       cstatus = 0

       ! Deallocate temp arrays
       call sparse_embed_destroy(old_sp_overlap)
       do isub = 1, nsub
          call data_functions_dealloc(old_ngwfs_on_grid(isub))
       end do
       call sparse_embed_array_destroy(old_denskern)

    end subroutine internal_projection

    !################################################################!

    subroutine internal_christoffel()

       !=============================================================!
       ! Restore tensorial integrity to last density kernel
       ! Modified by Andrea Greco to use complex matrices
       !=============================================================!

       use datatypes, only: FUNCTIONS, data_functions_alloc, &
            data_functions_dealloc
       use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
            dense_invert
       use function_ops, only: function_ops_brappd_ketppd
       use projectors, only: projectors_func_ovlp_box
       use rundat, only: pub_aug, pub_realspace_projectors, &
            pub_kpoint_method
       use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
            sparse_embed_create, sparse_embed_destroy, &
            sparse_embed_array_create, sparse_embed_array_destroy
       use utils, only: utils_abort, utils_assert
       ! agrecokpt: needed for testing in case of single k-point
       use geometry, only: POINT

       implicit none

       ! Local Variables
       type(FUNCTIONS) :: old_ngwfs_on_grid(nsub)
       type(SPAM3_EMBED_ARRAY) :: old_denskern
       type(SPAM3_EMBED) :: old_overlap
       type(SPAM3_EMBED) :: old_sp_overlap
       type(SPAM3_EMBED) :: old_inv_overlap
       type(DEM)   :: old_inv_overlap_dens
       integer  :: is, ierr
       ! agrecokpt
       type(POINT) :: loc_kpt
       integer :: tot_num, isub, jsub

       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : add tensorial corrections to &
          &density kernel'

       ! agrecokpt: currently only KP method is being implemented
       call utils_assert(pub_kpoint_method == 'KP', &
            'Subroutine internal_christoffel in elec_history_compose_dkn &
            &currently supports only KP method for BZ sampling')

       ! agrecokpt: single k-point for now, mainly for testing,
       ! in the future take k-points directly from denskern and
       ! NGWFs
       if (present(kpt)) then
          loc_kpt%x = kpt(1)
          loc_kpt%y = kpt(2)
          loc_kpt%z = kpt(3)
       else
          loc_kpt%x = 0.0_DP
          loc_kpt%y = 0.0_DP
          loc_kpt%z = 0.0_DP
       end if

       ! Retrieve density kernel from store
       call sparse_embed_array_create(old_denskern, denskern%kern)
       label = methods(iom)%ptr(1)
       call elec_history_read_dkn(old_denskern,'scf',label)

       ! Retrieve old ngwfs from store
       label = methods(iom)%ptr(1)
       tot_num=0
       ! jcap : Loop over regions
       do isub = 1, nsub
          par => mdl%regions(isub)%par
          tot_num=tot_num+ngwf_basis(isub)%num
          call data_functions_alloc(old_ngwfs_on_grid(isub), &
               ngwf_basis(isub)%size_on_grid, &
                iscmplx=denskern%kern%m(1,1)%m(isub,isub)%iscmplx)
          call elec_history_read_ngwfs(old_ngwfs_on_grid(isub),ngwf_basis(isub), &
               mdl%cell,mdl%fftbox,mdl%regions(isub)%elements,isub,&
               'scf',label,par=par)
       end do

       ! Create sparse arrays for christoffel corrections to density kernel
       call sparse_embed_create(old_sp_overlap,rep%sp_overlap)
       call sparse_embed_create(old_overlap,rep%overlap)
       call sparse_embed_create(old_inv_overlap,rep%inv_overlap)
       call dense_create(old_inv_overlap_dens,tot_num,tot_num, &
                         iscmplx=old_inv_overlap%p%iscmplx)

       ! If required, calculate overlap of old NGWFs with projectors
       if (pub_aug) then
          ! jcap: Loop over all pairs of regions
          do isub = 1, nsub
             do jsub = 1,nsub
                if (.not.pub_realspace_projectors) then
                   ! agrecokpt: include k-point dependence;
                   ! need to check this!
                   call projectors_func_ovlp_box(old_sp_overlap%m(isub,jsub), &
                        old_ngwfs_on_grid(isub),ngwf_basis(isub),&
                        proj_basis(jsub),nl_projectors(jsub),&
                        mdl%fftbox,mdl%cell,kshift=loc_kpt)
                else
                   ! agrecokpt: kpoint dependence included in definition of projs_on_grid
                   call function_ops_brappd_ketppd(old_sp_overlap%m(isub,jsub),&
                        old_ngwfs_on_grid(isub),ngwf_basis(isub),&
                        nl_projectors(jsub)%projs_on_grid,proj_basis(jsub),&
                        mdl%cell)
                end if
             end do
          end do
       endif

       ! Calculate overlap and inv_overlap of old NGWFs
       ! jcap: Loop over all pairs of regions
       do isub = 1, nsub
          do jsub = 1,nsub
             call function_ops_brappd_ketppd(old_overlap%m(isub,jsub), &
                  old_ngwfs_on_grid(isub), ngwf_basis(isub), &
                  old_ngwfs_on_grid(jsub), ngwf_basis(jsub), mdl%cell)
          end do
       end do
       call dense_convert(old_inv_overlap_dens, old_overlap)
       call dense_invert(old_inv_overlap_dens)
       call dense_convert(old_inv_overlap, old_inv_overlap_dens)

       ! Apply christoffel corrections to density kernel
       !call kernel_christoffel(ngwf_basis, rep%ngwfs_on_grid, &
       !     old_ngwfs_on_grid, denskern%kern, old_denskern, rep%inv_overlap, &
       !     old_inv_overlap, rep%sp_overlap, old_sp_overlap)
       call utils_abort('Error in internal_christoffel &
            &(elec_history_compose_dkn): Method not supported')
       cstatus = 0

       ! Deallocate temp arrays
       call dense_destroy(old_inv_overlap_dens)
       call sparse_embed_destroy(old_inv_overlap)
       call sparse_embed_destroy(old_overlap)
       call sparse_embed_destroy(old_sp_overlap)
       call sparse_embed_array_destroy(old_denskern)
       do isub = 1, nsub
          call data_functions_dealloc(old_ngwfs_on_grid(isub))
       end do

    end subroutine internal_christoffel

    !################################################################!


    subroutine internal_tr_propagation()

       !=============================================================!
       ! Propagation of auxiliary kernel :
       ! algorithm based on A. Niklasson et al,
       !           J. Chem. Phys. 130, 214109 (2009)
       !           section IIC, case K = 0
       !=============================================================!

       use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_array_create, &
            sparse_embed_array_scale, sparse_embed_array_axpy, &
            sparse_embed_array_destroy

       implicit none

       type(SPAM3_EMBED_ARRAY) :: old_denskern
       integer  :: is, ierr

       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : time-reversible propagation &
          &of auxiliary kernel'

       call sparse_embed_array_create(old_denskern, denskern%kern)
       call sparse_embed_array_scale(denskern%kern, 0.0_DP)

       label = methods(iom)%ptr(1)
       call elec_history_read_dkn(old_denskern, 'scf', label)
       call sparse_embed_array_axpy(denskern%kern, old_denskern, 2.0_dp)

       label = methods(iom)%ptr(2)
       call elec_history_read_dkn(old_denskern, 'init', label)
       call sparse_embed_array_axpy(denskern%kern, old_denskern, -1.0_dp)
       call elec_history_store_dkn(denskern%kern, 'init')
       cstatus = 0

       call sparse_embed_array_destroy(old_denskern)

    end subroutine internal_tr_propagation

    !################################################################!

    subroutine internal_xl_dissip_propagation()

       !=============================================================!
       ! Propagation of auxiliary kernel in the orthogonal represen- !
       ! tation :                                                    !
       ! algorithms based on A. Niklasson et al,                     !
       !           J. Chem. Phys. 130, 214109 (2009)                 !
       !           section IV, case 3 <= K <= 9                      !
       ! and Bowler et al                                            !
       !=============================================================!

       use comms,  only: pub_on_root
       use constants, only: DP, stdout
       use kernel, only: kernel_workspace_invalidate
       use rundat, only: pub_num_spins, md_aux_rep, PUB_1K, md_autocorr
       use smearing_operator, only: smearing_matrix, lowdin_transformation
       use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
            sparse_embed_create, sparse_embed_destroy, sparse_embed_copy, &
            sparse_embed_product, sparse_embed_array_create, &
            sparse_embed_array_scale, sparse_embed_array_axpy, &
            sparse_embed_array_destroy, sparse_embed_array_product, &
            sparse_embed_array_transpose

       implicit none

       ! Local variables
       type(SPAM3_EMBED_ARRAY)  :: tmp_dkn, old_x, new_x, dkn_trans, vel_x, tmp_x
       type(smearing_matrix) :: inv_sqrt_overlap, tmp_overlap
       character(len=2) :: structure_type
       integer  :: icc, istep, is
       integer  :: nkpts,nspins
       real(kind=DP) :: vel
       integer       :: xl_max = 10
       logical       :: loc_iscmplx
       real(kind=DP), parameter, dimension(7)   :: xl_k = (/  &
           1.690_dp, 1.750_dp, 1.820_dp, 1.840_dp, 1.860_dp, 1.880_dp, 1.890_dp /)
       real(kind=DP), parameter, dimension(70)  :: xl_coef = (/  &
       ! vv : Set of correct coefficients
             0.01000_dp, -0.55000_dp,  0.00000_dp, -0.15000_dp,  0.00000_dp,  0.00000_dp,&
             0.00000_dp,  0.00000_dp,  0.00000_dp,  0.00000_dp,&
             0.07900_dp, -0.65800_dp, -0.11400_dp, -0.11400_dp,  0.05700_dp,  0.00000_dp,&
             0.00000_dp,  0.00000_dp,  0.00000_dp,  0.00000_dp,&
             0.07200_dp, -0.74800_dp, -0.14400_dp, -0.05400_dp,  0.07200_dp, -0.01800_dp,&
             0.00000_dp,  0.00000_dp,  0.00000_dp,  0.00000_dp,&
             0.08300_dp, -0.80200_dp, -0.14850_dp, -0.01100_dp,  0.06600_dp, -0.03300_dp,&
             0.00550_dp,  0.00000_dp,  0.00000_dp,  0.00000_dp,&
             0.08240_dp, -0.84160_dp, -0.14080_dp,  0.01760_dp,  0.05120_dp, -0.04000_dp,&
             0.01280_dp, -0.00160_dp,  0.00000_dp,  0.00000_dp,&
             0.07644_dp, -0.87416_dp, -0.12584_dp,  0.03432_dp,  0.03432_dp, -0.03960_dp,&
             0.01848_dp, -0.00440_dp,  0.00044_dp,  0.00000_dp,&
             0.07568_dp, -0.89704_dp, -0.11232_dp,  0.04368_dp,  0.02016_dp, -0.03600_dp,&
             0.02208_dp, -0.00756_dp,  0.00144_dp, -0.00012_dp /)

       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : extended-Lagrangina propagation of &
          &auxiliary kernel with dissipation'

       nkpts = denskern%kern%num_kpoints
       nspins = denskern%kern%num_spins
       if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
          structure_type = 'K'
          inv_sqrt_overlap%matrix_type = 4
          tmp_overlap%matrix_type = 4
          inv_sqrt_overlap%dataSPAM3_EMBED%structure = 'K'
          tmp_overlap%dataSPAM3_EMBED%structure = 'S'
          ! vv: Allocate workspace
          call sparse_embed_create(inv_sqrt_overlap%dataSPAM3_EMBED)
          call sparse_embed_create(tmp_overlap%dataSPAM3_EMBED)
          call sparse_embed_copy(tmp_overlap%dataSPAM3_EMBED,rep%overlap)
          call sparse_embed_array_create(tmp_x,n_kpoints=nkpts,n_spins=nspins, &
               structure=structure_type)
       else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
          structure_type = 'KS'
          call sparse_embed_array_create(tmp_dkn,n_kpoints=nkpts,n_spins=pub_num_spins,&
               structure=structure_type)
          call sparse_embed_array_create(dkn_trans,n_kpoints=nkpts,n_spins=nspins,&
               structure=structure_type)
       end if

       call sparse_embed_array_create(old_x,n_kpoints=nkpts,n_spins=nspins,&
            structure=structure_type)
       call sparse_embed_array_create(new_x,n_kpoints=nkpts,n_spins=nspins,&
            structure=structure_type)
       call sparse_embed_array_create(vel_x,n_kpoints=nkpts,n_spins=nspins,&
            structure=structure_type)
       call sparse_embed_array_scale(denskern%kern,0.0_dp)

       if (md_autocorr) call elec_history_print_kernel(new_x,'real')
       ! Propagation of the auxiliary K.S with dissipation
       label = methods(iom)%ptr(1)
       call elec_history_read_dkn(old_x,'scf',label)

       icc = methods(iom)%dkn_mix_size-3
       call sparse_embed_array_axpy(new_x,old_x,xl_k(icc))

       do istep = 1, methods(iom)%dkn_mix_size
          label=methods(iom)%ptr(istep)
          call elec_history_read_dkn(old_x,'init',label)
          call sparse_embed_array_axpy(new_x,old_x, &
                  xl_coef((icc-1)*xl_max+istep))
       enddo
       call elec_history_store_dkn(new_x,'init')

       if (md_autocorr) call elec_history_print_kernel(new_x,'aux')

       ! vv: Go back to the non-orthogonal representation
       if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
          ! vv: Compute K = S.^(1/2).X.S^(1/2)
          call lowdin_transformation(tmp_overlap,inv_sqrt_overlap)

          do is = 1,pub_num_spins
             call sparse_embed_product(tmp_x%m(is,PUB_1K),new_x%m(is,PUB_1K),&
                  inv_sqrt_overlap%dataSPAM3_EMBED)
             call sparse_embed_product(denskern%kern%m(is,PUB_1K),&
                  inv_sqrt_overlap%dataSPAM3_EMBED,tmp_x%m(is,PUB_1K))
          end do

       else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
          ! vv: Compute K = X.S^(-1)
          call sparse_embed_array_product(tmp_dkn,new_x,rep%inv_overlap)

          ! vv: At this stage K and S do not commute, i.e. [K,S] \= 0.
          !     So we need to symmetrize the density kernel
          call sparse_embed_array_transpose(dkn_trans,tmp_dkn)

          call sparse_embed_array_axpy(denskern%kern,tmp_dkn,1.0_dp)

          call sparse_embed_array_axpy(denskern%kern,dkn_trans,1.0_dp)

          call sparse_embed_array_scale(denskern%kern,0.5_dp)

       end if

       call kernel_workspace_invalidate(denskern)

       cstatus = 0

       ! vv: Sick and destroy
       call sparse_embed_array_destroy(vel_x)
       call sparse_embed_array_destroy(old_x)
       call sparse_embed_array_destroy(new_x)
       if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
          call sparse_embed_destroy(inv_sqrt_overlap%dataSPAM3_EMBED)
          call sparse_embed_destroy(tmp_overlap%dataSPAM3_EMBED)
          call sparse_embed_array_destroy(tmp_x)
       else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
          call sparse_embed_array_destroy(tmp_dkn)
          call sparse_embed_array_destroy(dkn_trans)
       end if

    end subroutine internal_xl_dissip_propagation

    !################################################################!

    subroutine internal_xl_berendsen_propagation()

       !=============================================================!
       ! Propagation of auxiliary kernel in the orthogonal representa!
       ! tion K.S. :                                                 !
       ! algorithm based on                                          !
       !-------------------------------------------------------------!
       ! Implemented by Valerio Vitale                               !
       !=============================================================!

       use comms,  only: pub_on_root
       use constants, only: DP, stdout
       use kernel, only: kernel_workspace_invalidate
       use rundat, only: pub_num_spins, PUB_1K, md_aux_rep, &
           md_autocorr
       use smearing_operator, only: smearing_matrix,lowdin_transformation
       use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
            sparse_embed_create, sparse_embed_destroy, sparse_embed_copy, &
            sparse_embed_product, sparse_embed_array_axpy, &
            sparse_embed_array_create, sparse_embed_array_scale, &
            sparse_embed_array_destroy, sparse_embed_array_product, &
            sparse_embed_array_transpose

       implicit none

       type(SPAM3_EMBED_ARRAY) :: old_x, vel_x, new_x, tmp_x, tmp_dkn, dkn_trans
       type(smearing_matrix) :: inv_sqrt_overlap, tmp_overlap
       real(kind=DP) :: normS, kappa
       character(len=2) :: structure_type
       integer  :: is, ierr
       integer  :: nspins, nkpts

       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : extended-Lagrangian propagation of &
          &auxiliary kernel with Berendsen thermostat'

       nspins = denskern%kern%num_spins
       nkpts  = denskern%kern%num_kpoints

       if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
          structure_type = 'K'
          inv_sqrt_overlap%matrix_type = 4
          tmp_overlap%matrix_type = 4
          inv_sqrt_overlap%dataSPAM3_EMBED%structure = 'K'
          tmp_overlap%dataSPAM3_EMBED%structure = 'S'
          ! vv: Allocate workspace
          call sparse_embed_create(inv_sqrt_overlap%dataSPAM3_EMBED)
          call sparse_embed_create(tmp_overlap%dataSPAM3_EMBED)
          call sparse_embed_copy(tmp_overlap%dataSPAM3_EMBED,rep%overlap)
          call sparse_embed_array_create(tmp_x,n_kpoints=nkpts,n_spins=nspins, &
               structure=structure_type)

       else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
          structure_type = 'KS'
          call sparse_embed_array_create(tmp_dkn,n_kpoints=nkpts,n_spins=pub_num_spins,&
               structure='K')
          call sparse_embed_array_create(dkn_trans,n_kpoints=nkpts,n_spins=nspins,&
               structure='K')
       end if

       call sparse_embed_array_create(vel_x,n_kpoints=nkpts,n_spins=nspins, &
            structure=structure_type)
       call sparse_embed_array_create(new_x,n_kpoints=nkpts,n_spins=nspins, &
            structure=structure_type)
       call sparse_embed_array_create(old_x,n_kpoints=nkpts,n_spins=nspins, &
            structure=structure_type)
       call sparse_embed_array_scale(denskern%kern,0.0_DP)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! vv: velocity-Verlet update of auxiliary degrees of freedom for XLBOMD !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! vv: Set value for kappa
       kappa = 1.00_DP

       ! vv: Half-step velocity
       ! vv: Read variables from previous step
       label = methods(iom)%ptr(1)
       call elec_history_read_dkn(vel_x,'vel',label)
       call elec_history_read_dkn(old_x,'scf',label)
       call elec_history_read_dkn(new_x,'init',label)

       ! vv: d/dt(d/dt X(n)) = kappa(f(K,S)(n) - X(n))
       call sparse_embed_array_axpy(old_x,new_x,-1.0_DP)
       call sparse_embed_array_axpy(vel_x,old_x,kappa)

       ! print info to external files
       if (md_autocorr) call elec_history_print_kernel(new_x,'aux')

       ! vv: X(n+1) = X(n) + d/dt X(n) * Dt + kappa*(f(K,S)(n) - X(n))
       call sparse_embed_array_axpy(new_x,vel_x,1.0_DP)
       ! Store auxiliary degrees of freedom
       call elec_history_store_dkn(new_x,'init')
       call elec_history_store_dkn(vel_x,'vel')

       ! vv: Go back to the non-orthogonal representation K
       if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
          ! vv: Compute K = S.^(1/2).X.S^(1/2)
          call lowdin_transformation(tmp_overlap,inv_sqrt_overlap)

          do is = 1,pub_num_spins
             call sparse_embed_product(tmp_x%m(is,PUB_1K),new_x%m(is,PUB_1K),&
                  inv_sqrt_overlap%dataSPAM3_EMBED)
             call sparse_embed_product(denskern%kern%m(is,PUB_1K),&
                  inv_sqrt_overlap%dataSPAM3_EMBED,tmp_x%m(is,PUB_1K))
          end do

       else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
          ! vv: Compute K = X.S^(-1)
          call sparse_embed_array_product(tmp_dkn,new_x,rep%inv_overlap)

          ! vv: At this stage K and S do not commute, i.e. [K,S] \= 0.
          !     So we need to symmetrize the density kernel
          call sparse_embed_array_transpose(dkn_trans,tmp_dkn)

          call sparse_embed_array_axpy(denskern%kern,tmp_dkn,1.0_dp)

          call sparse_embed_array_axpy(denskern%kern,dkn_trans,1.0_dp)

          call sparse_embed_array_scale(denskern%kern,0.5_dp)

       end if

       call kernel_workspace_invalidate(denskern)

       ! Sick and destroy
       call sparse_embed_array_destroy(vel_x)
       call sparse_embed_array_destroy(new_x)
       call sparse_embed_array_destroy(old_x)
       if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
          call sparse_embed_destroy(inv_sqrt_overlap%dataSPAM3_EMBED)
          call sparse_embed_destroy(tmp_overlap%dataSPAM3_EMBED)
          call sparse_embed_array_destroy(tmp_x)
       else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
          call sparse_embed_array_destroy(tmp_dkn)
          call sparse_embed_array_destroy(dkn_trans)
       end if

       cstatus = 0

    end subroutine internal_xl_berendsen_propagation

    subroutine internal_xl_berendsen_scuseria_propagation()

       !=============================================================!
       ! Propagation of auxiliary kernel in the orthogonal representa!
       ! tion K.S. :                                                 !
       ! algorithm based on                                          !
       !-------------------------------------------------------------!
       ! Implemented by Valerio Vitale                               !
       !=============================================================!

       use comms,  only: pub_on_root
       use constants, only: DP, stdout
       use kernel, only: kernel_workspace_invalidate
       use rundat, only: pub_num_spins, PUB_1K, md_aux_rep, &
           md_autocorr
       use smearing_operator, only: smearing_matrix,lowdin_transformation
       use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
            sparse_embed_create, sparse_embed_destroy, sparse_embed_trace, &
            sparse_embed_copy, sparse_embed_product, sparse_embed_array_axpy, &
            sparse_embed_array_create, sparse_embed_array_scale, &
            sparse_embed_array_destroy, sparse_embed_array_trace, &
            sparse_embed_array_product, sparse_embed_array_transpose, &
            sparse_embed_array_num_rows, sparse_embed_array_copy
       use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

       implicit none

       type(SPAM3_EMBED_ARRAY) :: old_x, vel_x, new_x, tmp_x, tmp_dkn, dkn_trans, &
          pot, t, tx, xtx, q, tq, qtq
       type(smearing_matrix) :: inv_sqrt_overlap, tmp_overlap
       real(kind=DP) :: kappa, trace, loc_threshold
       real(kind=DP), allocatable :: trace_array(:,:)
       character(len=4) :: structure_type_x, structure_type_xx
       integer  :: is, ierr, counter
       integer  :: nspins, nkpts

       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : extended-Lagrangian propagation of &
          &auxiliary kernel with Berendsen thermostat and idempotency &
          &correction'

       nspins = denskern%kern%num_spins
       nkpts  = denskern%kern%num_kpoints

       if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
          structure_type_x = 'K'
          structure_type_xx = 'K'
          inv_sqrt_overlap%matrix_type = 4
          tmp_overlap%matrix_type = 4
          inv_sqrt_overlap%dataSPAM3_EMBED%structure = 'K'
          tmp_overlap%dataSPAM3_EMBED%structure = 'S'
          ! vv: Allocate workspace
          call sparse_embed_create(inv_sqrt_overlap%dataSPAM3_EMBED)
          call sparse_embed_create(tmp_overlap%dataSPAM3_EMBED)
          call sparse_embed_copy(tmp_overlap%dataSPAM3_EMBED,rep%overlap)

       else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
          structure_type_x = 'KS'
          structure_type_xx = 'KSKS'
          call sparse_embed_array_create(tmp_dkn,n_kpoints=nkpts,n_spins=pub_num_spins,&
               structure='K')
          call sparse_embed_array_create(dkn_trans,n_kpoints=nkpts,n_spins=nspins,&
               structure='K')
       end if

       call sparse_embed_array_create(vel_x,n_kpoints=nkpts,n_spins=nspins, &
            structure=structure_type_x)
       call sparse_embed_array_create(new_x,n_kpoints=nkpts,n_spins=nspins, &
            structure=structure_type_x)
       call sparse_embed_array_create(old_x,n_kpoints=nkpts,n_spins=nspins, &
            structure=structure_type_x)
       call sparse_embed_array_create(tmp_x,n_kpoints=nkpts,n_spins=nspins, &
            structure=structure_type_x)
       call sparse_embed_array_create(t,n_kpoints=nkpts,n_spins=nspins,     &
            structure=structure_type_x)
       call sparse_embed_array_create(q,n_kpoints=nkpts,n_spins=nspins,     &
            structure=structure_type_x)
       call sparse_embed_array_create(tx,n_kpoints=nkpts,n_spins=nspins,    &
            structure=structure_type_x)
       call sparse_embed_array_create(xtx,n_kpoints=nkpts,n_spins=nspins,   &
            structure=structure_type_x)
       call sparse_embed_array_create(tq,n_kpoints=nkpts,n_spins=nspins,    &
            structure=structure_type_x)
       call sparse_embed_array_create(qtq,n_kpoints=nkpts,n_spins=nspins,   &
            structure=structure_type_x)
       call sparse_embed_array_create(pot,n_kpoints=nkpts,n_spins=nspins,   &
            structure=structure_type_x)
       call sparse_embed_array_scale(denskern%kern,0.0_DP)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! vv: velocity-Verlet update of auxiliary degrees of freedom for XLBOMD !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! vv: Set value for kappa
       kappa = 1.00_DP

       ! vv: Half-step velocity
       ! vv: Read variables from previous step
       label = methods(iom)%ptr(1)
       call elec_history_read_dkn(vel_x,'vel',label)
       call elec_history_read_dkn(old_x,'scf',label)
       call elec_history_read_dkn(new_x,'init',label)
       call sparse_embed_array_copy(q,new_x)

       ! vv: d/dt(d/dt X(n)) = kappa(f(K,S)(n) - X(n))
       call sparse_embed_array_axpy(old_x,new_x,-1.0_DP)
       call sparse_embed_array_axpy(vel_x,old_x,kappa)
       call elec_history_store_dkn(vel_x,'vel')

       ! vv: X_i+1 = X_i + W_i + kappa*(f(K,S) - X)_i
       call sparse_embed_array_axpy(new_x,vel_x,1.0_DP)

       call sparse_embed_array_copy(old_x,q)
       call sparse_embed_array_scale(q,-1.0_DP,1.0_DP) ! Q=I - KS

       allocate(trace_array(pot%num_spins,pot%num_kpoints),stat=ierr)
       call utils_alloc_check('internal_xl_berendsen_scuseria_propagation',&
            'trace_array',ierr)

       trace = 1E5_DP
       loc_threshold = 1E-12_DP
       counter = 1
       do while (sqrt(trace)/sparse_embed_array_num_rows(new_x) >= loc_threshold)
          ! vv: XTX matrix, T = 3XX - 2XXX - X
          ! vv: T = -2X + 3
          call sparse_embed_array_copy(t,new_x)
          call sparse_embed_array_scale(t,-2.0_DP,3.0_DP)
          ! vv: T = T.X - I
          call sparse_embed_array_product(tmp_x,t,new_x)
          call sparse_embed_array_scale(tmp_x,1.0_DP,-1.0_DP)
          ! vv: T = T.X
          call sparse_embed_array_product(t,tmp_x,new_x)

          ! vv: X_i.T.X_i
          call sparse_embed_array_product(tx,t,old_x)
          call sparse_embed_array_product(xtx,old_x,tx)
          ! vv: Q_i.T.Q_i matrix
          call sparse_embed_array_product(tq,t,q)
          call sparse_embed_array_product(qtq,q,tq)

          ! vv: X_i+1 = X_i+1 + QTQ + XTX
          call sparse_embed_array_axpy(new_x,xtx,1.0_DP)
          call sparse_embed_array_axpy(new_x,qtq,1.0_DP)

          ! vv: Tr[(X^2 - X)^2]
          call sparse_embed_array_copy(tmp_x,new_x)
          call sparse_embed_array_scale(tmp_x,1.0_DP,-1.0_DP) ! X - I
          call sparse_embed_array_product(pot,tmp_x,new_x) ! (X - I)X = X^2 - X

          call sparse_embed_array_trace(trace_array,pot,pot)
          trace = sum(trace_array) ! Tr[ (X^2 - X)^2]
          counter = counter + 1
          call utils_assert(counter < 20,'Error in electronic_history:&
               & Iterative scheme for minimising Tr[(X.^2 - X).^2] has failed to&
               & convergence after 20 iterations! Aborting.')
       end do
       if (pub_on_root) write(stdout,'(a,i0,a)') &
          'Electronic history : Iterative scheme for minimising Tr[(X.^2 - X).^2]&
           & converged in ',counter,' iterations.'

       deallocate(trace_array,stat=ierr)
       call utils_dealloc_check('internal_xl_berendsen_scuseria_propagation',&
            'trace_array',ierr)

       ! vv: Store auxilary kernel
       call elec_history_store_dkn(new_x,'init')

       ! vv: Go back to the non-orthogonal representation K
       if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
          ! vv: Compute K = S.^(1/2).X.S^(1/2)
          call lowdin_transformation(tmp_overlap,inv_sqrt_overlap)

          do is = 1,pub_num_spins
             call sparse_embed_product(tmp_x%m(is,PUB_1K),new_x%m(is,PUB_1K),&
                  inv_sqrt_overlap%dataSPAM3_EMBED)
             call sparse_embed_product(denskern%kern%m(is,PUB_1K),&
                  inv_sqrt_overlap%dataSPAM3_EMBED,tmp_x%m(is,PUB_1K))
          end do

       else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
          ! vv: Compute K = X.S^(-1)
          call sparse_embed_array_product(tmp_dkn,new_x,rep%inv_overlap)

          ! vv: At this stage K and S do not commute, i.e. [K,S] \= 0.
          !     So we need to symmetrize the density kernel
          call sparse_embed_array_transpose(dkn_trans,tmp_dkn)

          call sparse_embed_array_axpy(denskern%kern,tmp_dkn,1.0_dp)

          call sparse_embed_array_axpy(denskern%kern,dkn_trans,1.0_dp)

          call sparse_embed_array_scale(denskern%kern,0.5_dp)

       end if

       call kernel_workspace_invalidate(denskern)

       ! print info to external files
       if (md_autocorr) call elec_history_print_kernel(new_x,'aux')
       ! Sick and destroy
       call sparse_embed_array_destroy(vel_x)
       call sparse_embed_array_destroy(new_x)
       call sparse_embed_array_destroy(old_x)
       call sparse_embed_array_destroy(tmp_x)
       call sparse_embed_array_destroy(t)
       call sparse_embed_array_destroy(q)
       call sparse_embed_array_destroy(tx)
       call sparse_embed_array_destroy(xtx)
       call sparse_embed_array_destroy(tq)
       call sparse_embed_array_destroy(qtq)
       call sparse_embed_array_destroy(pot)
       if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
          call sparse_embed_destroy(inv_sqrt_overlap%dataSPAM3_EMBED)
          call sparse_embed_destroy(tmp_overlap%dataSPAM3_EMBED)
       else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
          call sparse_embed_array_destroy(tmp_dkn)
          call sparse_embed_array_destroy(dkn_trans)
       end if

       cstatus = 0

    end subroutine internal_xl_berendsen_scuseria_propagation

  end subroutine elec_history_compose_dkn

!============================================================================!
!============================================================================!
!============================================================================!

  subroutine elec_history_store_ngwfs(ngwfs_on_grid, ngwf_basis, &
            cell, fftbox, elements, ireg, ngwfs_type, label, par)

    !========================================================================!
    ! This subroutine store the current NGWFs in a "universal tightbox"      !
    ! representation in list_ngwfs. The order of the NGWFs in list_ngwfs     !
    ! is the one in which they appear in the input file and is not affected  !
    ! by the use of the space-filling curve.                                 !
    !                                                                        !
    ! ndmh 18/12/2009: In cases where the FFTbox coincides with the cell     !
    ! but the individual tightboxes do not need to, a new temporary NGWF     !
    ! basis is constructed, for which the FFTbox only coincides with the     !
    ! cell along direction i if L_i < 2*R_NGWF + \delta_i.                   !
    !                                                                        !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 18/2/2010 based on the                !
    ! restart_ngwfs_tightbox_write subroutine originally written             !
    ! by Chris-Kriton Skylaris on 12/3/2004.                                 !
    ! Modified by Andrea Greco on 16/05/2015 to allow use of complex NGWFs.  !
    !========================================================================!

    use comms, only: comms_reduce, comms_barrier, comms_free, &
         pub_my_proc_id, pub_on_root
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS, function_basis_allocate, &
         function_basis_deallocate, function_basis_distribute, &
         function_basis_copy_spheres, function_basis_init_tight_boxes, &
         function_basis_ppds_to_tightbox
    use geometry, only: magnitude
    use ion, only : ELEMENT
    use parallel_strategy, only: PARAL_INFO!par=>pub_par
    use rundat, only: md_write_out
    use simulation_cell, only: CELL_INFO
    use utils, only : utils_alloc_check, utils_dealloc_check, &
         utils_assert

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in)     :: ngwf_basis
    type(CELL_INFO), intent(in)      :: cell
    type(FFTBOX_INFO), intent(inout) :: fftbox
    type(FUNCTIONS), intent(in)      :: ngwfs_on_grid
    type(PARAL_INFO), intent(in), pointer     :: par
    type(ELEMENT), intent(in)        :: elements(par%nat)
    character(len=*), intent(in)     :: ngwfs_type
    integer, optional, intent(inout) :: label
    integer, intent(in)              :: ireg

    ! Local variables
    type(FFTBOX_DATA) :: uni_tbox
    real(kind=DP) :: tb_orig1   ! a1-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: tb_orig2   ! a2-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: tb_orig3   ! a3-position of NGWF wrt tightbox in (real number of) grid points
    integer :: maxt_n1     ! a1-maximum tightbox points
    integer :: maxt_n2     ! a2-maximum tightbox points
    integer :: maxt_n3     ! a3-maximum tightbox points
    integer :: orig_atom   ! global atom counter in input file order
    integer :: atom_ngwf   ! ngwf of current atom counter
    integer :: distr_ngwf  ! global ngwf counter in spacefil order
    real(kind=DP) :: sendbuf(3) ! comms send pack buffer
    real(kind=DP) :: recvbuf(3) ! comms recv pack buffer
    real(kind=DP) :: max_radius ! Maximum NGWF radius
    integer :: glob_ngwf, loc_ngwf  ! global/local NGWFs counter
    type(FUNC_BASIS) :: tmp_basis ! ndmh: for dealing with FFTbox coinciding with cell
    logical :: orig_coin1, orig_coin2, orig_coin3 ! ndmh: whether the fftbox
                        ! originally coincided with cell in 1,2,3 directions
    integer :: ii
    integer :: storage_proc
    integer :: ierr

    ! agrecocmplx: variable for complex NGWFs
    logical :: loc_cmplx
    integer :: j

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Entering elec_history_store_ngwfs'

    ! Return if history has not been allocated
    if (.not. allocated(history)) return

    ! Look for the index of the storage entry
    if (present(label)) then
       ii = history_index(label=label)
    else
       ii = history_io
    endif

    ! Exit if index is out of range
    call utils_assert(ii <= history_max,'Error in elec_history_store_ngwfs : &
          &index is out of range!')

    ! Return if the corresponding entry has not been initialized
    if (history(ii)%label == 0) return

    ! cks: find maximum number of points for universal tightbox that has odd
    ! cks: number of points in each dimension
    maxt_n1 = ngwf_basis%maxtight_pts1
    maxt_n2 = ngwf_basis%maxtight_pts2
    maxt_n3 = ngwf_basis%maxtight_pts3

    ! ndmh: if the FFTbox coincides with the cell in any direction, but the
    ! ndmh: tightbox does not necessarily need to, we can define a
    ! ndmh: temporary NGWF basis where the tightboxes are as small as they
    ! ndmh: can be, so that we do not get simulation-cell-sized tightboxes
    ! ndmh: written out to files unnecessarily
    orig_coin1 = fftbox%coin1
    orig_coin2 = fftbox%coin2
    orig_coin3 = fftbox%coin3
    if (orig_coin1.or.orig_coin2.or.orig_coin3) then

    max_radius = maxval(ngwf_basis%spheres(:)%radius)
    call comms_reduce('MAX',max_radius)

    ! ndmh: temporarily override coin1,2,3 variables if possible
    if (magnitude(fftbox%a1) > 2.0_DP*max_radius + fftbox%d1) &
        fftbox%coin1 = .false.
    if (magnitude(fftbox%a2) > 2.0_DP*max_radius + fftbox%d2) &
        fftbox%coin2 = .false.
    if (magnitude(fftbox%a3) > 2.0_DP*max_radius + fftbox%d3) &
        fftbox%coin3 = .false.

    end if

    ! ndmh: if any of the coin1,2,3 variables have been changed...
    if ((orig_coin1.neqv.fftbox%coin1).or. &
        (orig_coin2.neqv.fftbox%coin2).or. &
        (orig_coin3.neqv.fftbox%coin3)) then

       ! ndmh: create a temporary function basis where tightboxes
       ! ndmh: do not coincide with cell if possible
       call function_basis_allocate(tmp_basis,ngwf_basis%num, &
           'tmp_'//ngwf_basis%name, par=par)
       call function_basis_distribute(tmp_basis,elements,par)
       call function_basis_copy_spheres(tmp_basis,ngwf_basis, par=par)
       call function_basis_init_tight_boxes(tmp_basis,fftbox,cell, par=par)

       ! ndmh: set a universal tightbox size which depends on the actual size
       ! ndmh: of the temporary NGWF tightboxes
       maxt_n1 = maxval(tmp_basis%tight_boxes(:)%tight_pts1)
       maxt_n2 = maxval(tmp_basis%tight_boxes(:)%tight_pts2)
       maxt_n3 = maxval(tmp_basis%tight_boxes(:)%tight_pts3)
       call comms_reduce('MAX',maxt_n1)
       call comms_reduce('MAX',maxt_n2)
       call comms_reduce('MAX',maxt_n3)

    end if

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid%iscmplx

    ! cks: allocate universal tightbox buffer
    !allocate(uni_tbox(maxt_n1,maxt_n2,maxt_n3),stat=ierr)
    !call utils_alloc_check('restart_ngwfs_tightbox_store','uni_tbox',ierr)
    ! agrecocmplx: allocate using appropriate routine
    call data_fftbox_alloc(uni_tbox, maxt_n1, maxt_n2, maxt_n3, &
         iscmplx=loc_cmplx)

    ! Allocate the storage arrays
    if (.not. allocated(history(ii)%ngwfs(ireg)%orig)) then
       allocate(history(ii)%ngwfs(ireg)%orig(3,ngwf_basis%proc_num),stat=ierr)
       call utils_alloc_check('elec_history_store_ngwfs','history(ii)%ngwfs%orig',ierr)
    endif

    if (.not. allocated(history(ii)%ngwfs(ireg)%procs)) then
       allocate(history(ii)%ngwfs(ireg)%procs(ngwf_basis%num),stat=ierr)
       call utils_alloc_check('elec_history_store_ngwfs','history(ii)%ngwfs%procs',ierr)
    endif

    history(ii)%ngwfs(ireg)%tb_size(:) = (/ maxt_n1, maxt_n2, maxt_n3 /)

    if (ngwfs_type == 'init') then
       if (.not. allocated(history(ii)%ngwfs(ireg)%init)) then

          !allocate(history(ii)%ngwfs%init(maxt_n1,maxt_n2,maxt_n3, &
          !   ngwf_basis%proc_num),stat=ierr)
          !call utils_alloc_check('elec_history_store_ngwfs','init',ierr)
          ! agrecocmplx: allocate using appropriate routines
          allocate(history(ii)%ngwfs(ireg)%init(ngwf_basis%proc_num), stat=ierr)
          call utils_alloc_check('elec_history_store_ngwfs','history(ii)%ngwfs%init',ierr)
          do j = 1, ngwf_basis%proc_num
             call data_fftbox_alloc(history(ii)%ngwfs(ireg)%init(j),maxt_n1, &
                  maxt_n2, maxt_n3, iscmplx=loc_cmplx)
          end do

       endif

       if (pub_on_root .and. pub_md_output_detail>=NORMAL .and. &
           (md_write_out)) &
          write(stdout,'(a,i4,a,i4,a)') &
          'Electronic history : store auxiliary ngwfs in entry &
          &(idx = ',ii,'), (label = ',history(ii)%label,')'

    else if (ngwfs_type.eq.'scf') then
       if (.not. allocated(history(ii)%ngwfs(ireg)%scf)) then

          !allocate(history(ii)%ngwfs%scf(maxt_n1,maxt_n2,maxt_n3, &
          !   ngwf_basis%proc_num),stat=ierr)
          !call utils_alloc_check('elec_history_store_ngwfs','scf',ierr)
          ! agrecocmplx: allocate using appropriate routines
          allocate(history(ii)%ngwfs(ireg)%scf(ngwf_basis%proc_num), stat=ierr)
          call utils_alloc_check('elec_history_store_ngwfs','history(ii)%ngwfs%scf',ierr)
          do j = 1, ngwf_basis%proc_num
              call data_fftbox_alloc(history(ii)%ngwfs(ireg)%scf(j),maxt_n1, &
                  maxt_n2, maxt_n3, iscmplx=loc_cmplx)
          end do

       endif

       if (pub_on_root .and. pub_md_output_detail>=NORMAL .and. &
           (md_write_out)) &
          write(stdout,'(a,i4,a,i4,a)') &
          'Electronic history : store ngwfs in entry (idx = ',&
          ii,'), (label = ',history(ii)%label,')'

    else if (ngwfs_type.eq.'first') then
       if (.not. allocated(history(ii)%ngwfs(ireg)%first)) then

          !allocate(history(ii)%ngwfs%scf(maxt_n1,maxt_n2,maxt_n3, &
          !   ngwf_basis%proc_num),stat=ierr)
          !call utils_alloc_check('elec_history_store_ngwfs','scf',ierr)
          ! agrecocmplx: allocate using appropriate routines
          allocate(history(ii)%ngwfs(ireg)%first(ngwf_basis%proc_num), stat=ierr)
          call utils_alloc_check('elec_history_store_ngwfs','history(ii)%ngwfs%first',ierr)
          do j = 1, ngwf_basis%proc_num
             call data_fftbox_alloc(history(ii)%ngwfs(ireg)%first(j),maxt_n1, &
                  maxt_n2, maxt_n3, iscmplx=loc_cmplx)
          end do

       endif

       if (pub_on_root .and. pub_md_output_detail>=NORMAL .and. &
           (md_write_out)) &
          write(stdout,'(a,i4,a,i4,a)') &
          'Electronic history : store original ngwfs in entry (idx = ',&
          ii,'), (label = ',history(ii)%label,')'

    endif

    ! cks: Loop over all atoms in the order they appear in the input file
    glob_ngwf = 0
    loc_ngwf = 0

    atom_loop: do orig_atom =1, par%nat

       ngwfs_on_atom_loop: do atom_ngwf= &
           1,ngwf_basis%num_on_atom(par%distr_atom(orig_atom))

         ! smmd: global ngwf index (before distribution)
         glob_ngwf = glob_ngwf + 1

         ! smmd: global ngwf index (after distribution)
         distr_ngwf = ngwf_basis%first_on_atom(par%distr_atom(orig_atom)) + &
              atom_ngwf - 1

         ! smmd: storage proc
         storage_proc = ngwf_basis%proc_of_func(distr_ngwf)

         ! cks: this essentially converts the non-blocking sends of
         ! cks: comms_send to blocking sends. Since parallel performance
         ! cks: is not the issue here, this makes the code simpler.
         call comms_barrier

         if ((orig_coin1.neqv.fftbox%coin1).or. &
              (orig_coin2.neqv.fftbox%coin2).or. &
              (orig_coin3.neqv.fftbox%coin3)) then

            call function_basis_ppds_to_tightbox(uni_tbox, &
                 tb_orig1,tb_orig2,tb_orig3,distr_ngwf,storage_proc, &
                 ngwfs_on_grid,tmp_basis,fftbox,cell,sendbuf,recvbuf)
         else
            call function_basis_ppds_to_tightbox(uni_tbox, &
                 tb_orig1,tb_orig2,tb_orig3,distr_ngwf,storage_proc, &
                 ngwfs_on_grid,ngwf_basis,fftbox,cell,sendbuf,recvbuf)
         end if

         ! Store current NGWF on storage_proc
         if (pub_my_proc_id==storage_proc) then

            loc_ngwf = loc_ngwf + 1
            history(ii)%ngwfs(ireg)%orig(:,loc_ngwf) = &
               (/ tb_orig1, tb_orig2, tb_orig3 /)
            if (ngwfs_type.eq.'init') then
               !history(ii)%ngwfs%init(:,:,:,loc_ngwf) = &
               !   uni_tbox(1:maxt_n1,1:maxt_n2,1:maxt_n3)
               ! agrecocmplx
               if (loc_cmplx) then
                   history(ii)%ngwfs(ireg)%init(loc_ngwf)%z(:,:,:) = &
                       uni_tbox%z(1:maxt_n1,1:maxt_n2,1:maxt_n3)
               else
                   history(ii)%ngwfs(ireg)%init(loc_ngwf)%d(:,:,:) = &
                       uni_tbox%d(1:maxt_n1,1:maxt_n2,1:maxt_n3)
               end if
            else if (ngwfs_type.eq.'scf') then
               !history(ii)%ngwfs%scf(:,:,:,loc_ngwf) = &
               !   uni_tbox(1:maxt_n1,1:maxt_n2,1:maxt_n3)
               ! agrecocmplx
               if (loc_cmplx) then
                   history(ii)%ngwfs(ireg)%scf(loc_ngwf)%z(:,:,:) = &
                       uni_tbox%z(1:maxt_n1,1:maxt_n2,1:maxt_n3)
               else
                   history(ii)%ngwfs(ireg)%scf(loc_ngwf)%d(:,:,:) = &
                       uni_tbox%d(1:maxt_n1,1:maxt_n2,1:maxt_n3)
               end if

            else if (ngwfs_type.eq.'first') then
               !history(ii)%ngwfs%scf(:,:,:,loc_ngwf) = &
               !   uni_tbox(1:maxt_n1,1:maxt_n2,1:maxt_n3)
               ! agrecocmplx
               if (loc_cmplx) then
                   history(ii)%ngwfs(ireg)%first(loc_ngwf)%z(:,:,:) = &
                       uni_tbox%z(1:maxt_n1,1:maxt_n2,1:maxt_n3)
               else
                   history(ii)%ngwfs(ireg)%first(loc_ngwf)%d(:,:,:) = &
                       uni_tbox%d(1:maxt_n1,1:maxt_n2,1:maxt_n3)
               end if

            endif
         end if

         ! Keep the identity of storage_proc in memory
         history(ii)%ngwfs(ireg)%procs(glob_ngwf) = storage_proc

       enddo ngwfs_on_atom_loop
    enddo atom_loop

    ! cks: don't go away just yet! Wait for all communications to finish!
    ! cks: also wait before deallocating uni_tbox which is used in non-blocking
    ! cks: send!
    call comms_barrier
    ! cks: clean up send stack too
    call comms_free

    ! cks: deallocate universal tightbox buffer
    !deallocate(uni_tbox,stat=ierr)
    !call utils_dealloc_check('elec_history_store_ngwfs','uni_tbox',ierr)
    ! agrecocmplx: deallocate using appropriate routine
    call data_fftbox_dealloc(uni_tbox)

    ! ndmh: Clean up temporary basis if we created one, and reset coin1,2,3
    if ((orig_coin1.neqv.fftbox%coin1).or. &
        (orig_coin2.neqv.fftbox%coin2).or. &
        (orig_coin3.neqv.fftbox%coin3)) then

       call function_basis_deallocate(tmp_basis)
       fftbox%coin1 = orig_coin1
       fftbox%coin2 = orig_coin2
       fftbox%coin3 = orig_coin3

    end if

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Leave elec_history_store_ngwfs'

    return

  end subroutine elec_history_store_ngwfs


!============================================================================!
!============================================================================!
!============================================================================!

  subroutine elec_history_read_ngwfs(ngwfs_on_grid, ngwf_basis, &
       cell, fftbox, elements, ireg, ngwfs_type, label, coeffs, par)

    !========================================================================!
    ! This subroutine retrieve the  NGWFs in "universal tightbox"            !
    ! representation from list_ngwfs. The order of the NGWFs in list_ngwfs   !
    ! is the one in which they appear in the input file and is not affected  !
    ! by the use of the space-filling curve.                                 !
    !                                                                        !
    ! ndmh 18/12/2009: In cases where the FFTbox coincides with the cell     !
    ! but the individual tightboxes do not need to, a new temporary NGWF     !
    ! basis is constructed, for which the FFTbox only coincides with the     !
    ! cell along direction i if L_i < 2*R_NGWF + \delta_i.                   !
    !                                                                        !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 18/2/2010 based on the                !
    ! restart_ngwfs_tightbox_read subroutine originally written              !
    ! by Chris-Kriton Skylaris on 12/3/2004.                                 !
    ! Modified by Andrea Greco on 16/05/2015 to allow use of complex NGWFs.  !
    !========================================================================!

    use basis, only: basis_ket_start_wrt_fftbox
    use comms, only: comms_barrier, comms_bcast, comms_free, &
         comms_reduce, pub_my_proc_id, pub_on_root, pub_root_proc_id
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
        data_fftbox_dealloc, data_functions_alloc, data_functions_dealloc, &
        data_functions_axpy, data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS, function_basis_allocate, &
         function_basis_deallocate, function_basis_distribute, &
         function_basis_copy_spheres, function_basis_init_tight_boxes, &
         function_basis_tightbox_to_ppds
    use geometry, only: magnitude
    use ion, only : ELEMENT
    use parallel_strategy, only: PARAL_INFO!par=>pub_par
    use rundat, only: md_write_out
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_assert

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(inout) :: fftbox
    type(FUNCTIONS), intent(inout) :: ngwfs_on_grid ! NGWFs on this proc
    type(PARAL_INFO), intent(in), pointer :: par ! parallel strategy
    type(ELEMENT), intent(in) :: elements(par%nat) ! elements of all proc (in input file order)
    character(len=*), intent(in) :: ngwfs_type
    integer, intent(in)  :: label
    real(kind=DP), intent(in), optional :: coeffs(par%nat)
    integer, intent(in) :: ireg ! jcap: region counter


    ! Local Variables
    type(FFTBOX_DATA) :: uni_tbox  ! universal tightbox
    complex(kind=DP), allocatable, dimension(:,:,:) :: fftbox_complex ! complex fftbox space
    complex(kind=DP), allocatable, dimension(:,:,:) :: fftbox_complex_shifted ! complex fftbox space
    type(FFTBOX_DATA) :: fftbox_buffer ! real fftbox space
    real(kind=DP) :: max_radius          ! globaly maximum ngwf radius
    real(kind=DP) :: read_tb_orig1 ! read a1-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: read_tb_orig2 ! read a2-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: read_tb_orig3 ! read a3-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: read_weight   ! read grid point weight
    integer :: n1, n2, n3, ld1, ld2 ! fftbox dimensions
    integer :: orig_atom    ! global atom counter in input file order
    integer :: atom_ngwf    ! ngwf of current atom counter
    integer :: maxt_n1      ! a1-maximum tightbox points
    integer :: maxt_n2      ! a2-maximum tightbox points
    integer :: maxt_n3      ! a3-maximum tightbox points
    integer :: read_maxt_n1 ! a1-maximum tightbox points read from file
    integer :: read_maxt_n2 ! a2-maximum tightbox points read from file
    integer :: read_maxt_n3 ! a3-maximum tightbox points read from file
    integer :: tb_start1    ! a1-grid point start of universal tightbox wrt fftbox
    integer :: tb_start2    ! a2-grid point start of universal tightbox wrt fftbox
    integer :: tb_start3    ! a3-grid point start of universal tightbox wrt fftbox
    integer :: distr_ngwf   ! global ngwf counter in spacefill order
    real(kind=DP) :: sendbuf(3) ! comms send pack buffer
    real(kind=DP) :: recvbuf(3) ! comms recv pack buffer
    integer :: ierr         ! memory allocation error flag

    type(FUNCTIONS) :: buffer_on_grid
    integer :: glob_ngwf, loc_ngwf  ! global/local NGWFs counter
    integer :: ii, storage_proc

    ! ndmh:
    type(FUNC_BASIS) :: tmp_basis
    logical :: orig_coin1,orig_coin2,orig_coin3 ! ndmh: whether the fftbox
                           ! originally coincided with cell in 1,2,3 directions

    ! agrecocmplx: local variable to specify complex NGWFs
    logical :: loc_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Entering elec_history_read_ngwfs'

    ! Exit if history has not been allocated
    call utils_assert(allocated(history),'Error in elec_history_read_ngwfs : &
          &history has not been allocated!')
    ii = history_index(label=label)

    ! Exit if index is out of range
    call utils_assert(ii <= history_max,'Error in elec_history_read_ngwfs : &
          &index is out of range!')

    ! Exit if the corresponding entry has not been initialized
    call utils_assert(history(ii)%label /= 0, &
          'Error in elec_history_read_ngwfs : &
          &corresponding entry has not been initialized!')

    if (pub_on_root .and. pub_md_output_detail>=NORMAL .and. &
        (md_write_out))  then
       if (ngwfs_type.eq.'init') then
          write(stdout,'(a,i2,a,i4,a)') 'Electronic history : &
          &read auxiliary ngwfs from entry (idx = ',ii,')&
          &(label = ',history(ii)%label,')'
       else if (ngwfs_type .eq. 'scf') then
          write(stdout,'(a,i2,a,i4,a)') 'Electronic history : &
          &read ngwfs from entry (idx = ',ii,')&
          &(label = ',history(ii)%label,')'
       else if (ngwfs_type .eq. 'first') then
          write(stdout,'(a,i2,a,i4,a)') 'Electronic history : &
          &read original ngwfs from entry (idx = ',ii,')&
          &(label = ',history(ii)%label,')'
       endif
    endif

    ! cks: ********** INITIALISATIONS ***************************************
    !ngwfs_on_grid = 0.0_DP
    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid%iscmplx

    call data_set_to_zero(ngwfs_on_grid)

    ! cks: find maximum number of (odd numbers) points for universal tightbox
    max_radius = maxval(ngwf_basis%spheres(:)%radius)
    call comms_reduce('MAX',max_radius)
    maxt_n1 = int(2.0_DP*max_radius / cell%d1) + 1
    maxt_n1 = maxt_n1 + 1 - mod(maxt_n1,2)
    maxt_n2 = int(2.0_DP*max_radius / cell%d2) + 1
    maxt_n2 = maxt_n2 + 1 - mod(maxt_n2,2)
    maxt_n3 = int(2.0_DP*max_radius / cell%d3) + 1
    maxt_n3 = maxt_n3 + 1 - mod(maxt_n3,2)

    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2

    read_maxt_n1 = history(ii)%ngwfs(ireg)%tb_size(1)
    read_maxt_n2 = history(ii)%ngwfs(ireg)%tb_size(2)
    read_maxt_n3 = history(ii)%ngwfs(ireg)%tb_size(3)
    ! cks: ****** END INITIALISATIONS ***************************************


    ! ndmh: if the FFTbox coincides with the cell in any direction, but the
    ! ndmh: tightbox does not necessarily need to, we may need to define a
    ! ndmh: temporary NGWF basis where the tightboxes are as small as they
    ! ndmh: can be, unless the file we are reading in was written with
    ! ndmh: tightboxes coinciding with cell
    orig_coin1 = fftbox%coin1
    orig_coin2 = fftbox%coin2
    orig_coin3 = fftbox%coin3
    if (orig_coin1.or.orig_coin2.or.orig_coin3) then

       ! ndmh: temporarily override coin1,2,3 variables if possible (as long
       ! ndmh: as the written values are not for tightbox==fftbox)
       if ((magnitude(fftbox%a1) > 2.0_DP*max_radius + &
            fftbox%d1*cell%n_pt1).and. &
            (read_maxt_n1 /= fftbox%total_pt1)) fftbox%coin1 = .false.
       if ((magnitude(fftbox%a2) > 2.0_DP*max_radius + &
            fftbox%d2*cell%n_pt1).and. &
            (read_maxt_n2 /= fftbox%total_pt2)) fftbox%coin2 = .false.
       if ((magnitude(fftbox%a3) > 2.0_DP*max_radius + &
            fftbox%d3*cell%n_pt1).and. &
            (read_maxt_n3 /= fftbox%total_pt3)) fftbox%coin3 = .false.

       if(pub_debug_on_root) then
          if ((orig_coin1.neqv.fftbox%coin1)) write(stdout,*) &
               'DEBUG: Temporarily overriding fftbox%coin1'
          if ((orig_coin2.neqv.fftbox%coin2)) write(stdout,*) &
               'DEBUG: Temporarily overriding fftbox%coin2'
          if ((orig_coin3.neqv.fftbox%coin3)) write(stdout,*) &
               'DEBUG: Temporarily overriding fftbox%coin3'
       end if
    end if

    ! ndmh: if any of the coin1,2,3 variables have been changed...
    if ((orig_coin1.neqv.fftbox%coin1).or. &
         (orig_coin2.neqv.fftbox%coin2).or. &
         (orig_coin3.neqv.fftbox%coin3)) then

       ! ndmh: create a temporary function basis where tightboxes
       ! ndmh: do not coincide with cell if possible
       call function_basis_allocate(tmp_basis,ngwf_basis%num, &
            'tmp_'//ngwf_basis%name, par=par)
       call function_basis_distribute(tmp_basis,elements,par)
       call function_basis_copy_spheres(tmp_basis, ngwf_basis, par=par)
       call function_basis_init_tight_boxes(tmp_basis, fftbox, cell, par=par)

       ! ndmh: set a universal tightbox size which depends on the actual size
       ! ndmh: of the temporary NGWF tightboxes
       maxt_n1 = maxval(tmp_basis%tight_boxes(:)%tight_pts1)
       maxt_n2 = maxval(tmp_basis%tight_boxes(:)%tight_pts2)
       maxt_n3 = maxval(tmp_basis%tight_boxes(:)%tight_pts3)
       call comms_reduce('MAX',maxt_n1)
       call comms_reduce('MAX',maxt_n2)
       call comms_reduce('MAX',maxt_n3)

    end if

    ! cks: NGWFs will be placed in the centre of the fftbox
    ! cks: or left where they are in the simulation cell
    ! cks: (when the simulation cell and FFT-box coincide)
    call basis_ket_start_wrt_fftbox(tb_start1, tb_start2, tb_start3, &
         n1, n2, n3, fftbox)

    ! ndmh: check tightbox will not extend out of FFTbox if the coin flag has
    ! ndmh: been overridden. If so, try to adjust tb_start to compensate. If
    ! ndmh: that fails, print an error and quit.
    if ((tb_start1 + maxt_n1 > n1).and.(.not.fftbox%coin1)) then

       if (pub_debug_on_root) write(stdout,*) &
            'DEBUG: Attempting to adjust tb_start1'

       call utils_assert(n1-maxt_n1 > 0,'Error in elec_history_read_ngwfs: &
               &Tightbox size along 1-axis is not compatible with FFT-box size.&
               & The values were: ',maxt_n1, n1)
       tb_start1 = (n1 - maxt_n1)/2
    end if
    if ((tb_start2 + maxt_n2 > n2).and.(.not.fftbox%coin2)) then

       if (pub_debug_on_root) write(stdout,*) &
            'DEBUG: Attempting to adjust tb_start2'

       call utils_assert(n2-maxt_n2 > 0,'Error in elec_history_read_ngwfs: &
               &Tightbox size along 2-axis is not compatible with FFT-box size.&
               & The values were: ',maxt_n2, n2)
       tb_start2 = (n2 - maxt_n2)/2
    end if
    if ((tb_start3 + maxt_n3 > n3).and.(.not.fftbox%coin3)) then

       if (pub_debug_on_root) write(stdout,*) &
            'DEBUG: Attempting to adjust tb_start3'

       call utils_assert(n3-maxt_n3 > 0,'Error in elec_history_read_ngwfs: &
               &Tightbox size along 3-axis is not compatible with FFT-box size.&
               & The values were: ',maxt_n3, n3)
       tb_start3 = (n3 - maxt_n3)/2
    end if

    ! cks: allocate universal tightbox buffer
    maxt_n1 =max(min(maxt_n1,n1), read_maxt_n1)
    maxt_n2 =max(min(maxt_n2,n2), read_maxt_n2)
    maxt_n3 =max(min(maxt_n3,n3), read_maxt_n3)

    ! cks: send the new maxt's from the root to all other procs
    call comms_bcast(pub_root_proc_id, maxt_n1)
    call comms_bcast(pub_root_proc_id, maxt_n2)
    call comms_bcast(pub_root_proc_id, maxt_n3)

    ! cks: send read_weight to all procs
    call comms_bcast(pub_root_proc_id, read_weight)

    ! ndmh: allocate temporary storage
    !allocate(uni_tbox(maxt_n1,maxt_n2,maxt_n3),stat=ierr)
    !call utils_alloc_check('elec_history_read_ngwfs', &
    !     'uni_tbox',ierr)
    ! agrecocmplx: allocate using appropriate routine
    call data_fftbox_alloc(uni_tbox, maxt_n1, maxt_n2, maxt_n3, &
        iscmplx=loc_cmplx)
    allocate(fftbox_complex(ld1,ld2,n3),stat=ierr)
    call utils_alloc_check('elec_history_read_ngwfs', &
         'fftbox_complex',ierr)
    allocate(fftbox_complex_shifted(ld1,ld2,n3),stat=ierr)
    call utils_alloc_check('elec_history_read_ngwfs', &
         'fftbox_complex_shifted',ierr)
    !allocate(fftbox_buffer(ld1,ld2,n3),stat=ierr)
    !call utils_alloc_check('elec_history_read_ngwfs', &
    !     'fftbox_buffer',ierr)
    ! agrecocmplx: allocate using appropriate routine
    call data_fftbox_alloc(fftbox_buffer, ld1, ld2, n3, &
        iscmplx=loc_cmplx)
    !allocate(buffer_on_grid(ngwf_basis%size_on_grid), &
    !     stat=ierr)
    !call utils_alloc_check('elec_history_read_ngwfs', &
    !     'buffer_on_grid',ierr)
    ! agrecocmplx: allocate using appropriate routine
    call data_functions_alloc(buffer_on_grid, ngwf_basis%size_on_grid, &
        iscmplx=loc_cmplx)

    ! Initialisation
    loc_ngwf = 0
    glob_ngwf = 0
    !ngwfs_on_grid(:) = 0.d0
    ! agrecocmplx
    call data_set_to_zero(ngwfs_on_grid)

    ! cks: Loop over all atoms in the order they appear in the input file
    atom_loop: do orig_atom =1, par%nat

       ! Initialise buffer
       !buffer_on_grid(:) = 0.0d0
       ! agrecocmplx
       call data_set_to_zero(buffer_on_grid)

       ngwfs_on_atom_loop: do atom_ngwf=1, &
            ngwf_basis%num_on_atom(par%distr_atom(orig_atom))

          ! smmd: global ngwfs index before distribution
          glob_ngwf = glob_ngwf + 1

          ! smmd: global ngwfs index after distribution
          distr_ngwf = ngwf_basis%first_on_atom(par%distr_atom(orig_atom)) + &
               atom_ngwf - 1

          ! smmd: storage proc
          storage_proc = history(ii)%ngwfs(ireg)%procs(glob_ngwf)

          ! cks: this essentially converts the non-blocking sends of
          ! cks: comms_send to blocking sends. Since parallel performance
          ! cks: is not the issue here, this makes the code simpler.
          call comms_barrier

          ! cks: initialise
          ! uni_tbox = 0.0_DP
          ! agrecocmplx
          call data_set_to_zero(uni_tbox)

          ! smmd: if storage proc, read data into buffer
          if (pub_my_proc_id == storage_proc) then

             loc_ngwf = loc_ngwf + 1

             read_tb_orig1 = history(ii)%ngwfs(ireg)%orig(1,loc_ngwf)
             read_tb_orig2 = history(ii)%ngwfs(ireg)%orig(2,loc_ngwf)
             read_tb_orig3 = history(ii)%ngwfs(ireg)%orig(3,loc_ngwf)

             if (ngwfs_type.eq.'init') then

                !uni_tbox(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                !    history(ii)%ngwfs%init(1:read_maxt_n1,1:read_maxt_n2, &
                !    1:read_maxt_n3,loc_ngwf)
                ! agrecocmplx
                if (loc_cmplx) then
                    uni_tbox%z(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                         history(ii)%ngwfs(ireg)%init(loc_ngwf)%z(1:read_maxt_n1,&
                         1:read_maxt_n2,1:read_maxt_n3)
                else
                    uni_tbox%d(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                         history(ii)%ngwfs(ireg)%init(loc_ngwf)%d(1:read_maxt_n1,&
                         1:read_maxt_n2,1:read_maxt_n3)
                endif
             else if (ngwfs_type.eq.'scf') then

                !uni_tbox(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                !    history(ii)%ngwfs%scf(1:read_maxt_n1,1:read_maxt_n2, &
                !    1:read_maxt_n3,loc_ngwf)
                ! agrecocmplx
                if (loc_cmplx) then
                    uni_tbox%z(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                         history(ii)%ngwfs(ireg)%scf(loc_ngwf)%z(1:read_maxt_n1,&
                         1:read_maxt_n2,1:read_maxt_n3)
                else
                    uni_tbox%d(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                         history(ii)%ngwfs(ireg)%scf(loc_ngwf)%d(1:read_maxt_n1,&
                         1:read_maxt_n2,1:read_maxt_n3)
                endif
             else if (ngwfs_type.eq.'first') then

                !uni_tbox(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                !    history(ii)%ngwfs%scf(1:read_maxt_n1,1:read_maxt_n2, &
                !    1:read_maxt_n3,loc_ngwf)
                ! agrecocmplx
                if (loc_cmplx) then
                    uni_tbox%z(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                         history(ii)%ngwfs(ireg)%first(loc_ngwf)%z(1:read_maxt_n1,&
                         1:read_maxt_n2,1:read_maxt_n3)
                else
                    uni_tbox%d(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                         history(ii)%ngwfs(ireg)%first(loc_ngwf)%d(1:read_maxt_n1,&
                         1:read_maxt_n2,1:read_maxt_n3)
                endif
             endif

          endif


          if ((orig_coin1.neqv.fftbox%coin1).or. &
            (orig_coin2.neqv.fftbox%coin2).or. &
            (orig_coin3.neqv.fftbox%coin3)) then

             call function_basis_tightbox_to_ppds(uni_tbox, &
                  read_tb_orig1,read_tb_orig2,read_tb_orig3, &
                  distr_ngwf,storage_proc,maxt_n1,maxt_n2,maxt_n3, &
                  buffer_on_grid,tmp_basis,fftbox,cell,tb_start1,tb_start2,tb_start3, &
                  fftbox_complex,fftbox_complex_shifted,fftbox_buffer, &
                  sendbuf,recvbuf)

          else

             call function_basis_tightbox_to_ppds(uni_tbox, &
                  read_tb_orig1,read_tb_orig2,read_tb_orig3, &
                  distr_ngwf,storage_proc,maxt_n1,maxt_n2,maxt_n3, &
                  buffer_on_grid,ngwf_basis,fftbox,cell,tb_start1,tb_start2,tb_start3, &
                  fftbox_complex,fftbox_complex_shifted,fftbox_buffer, &
                  sendbuf,recvbuf)

          end if

       enddo ngwfs_on_atom_loop

       if (present(coeffs)) then
          !ngwfs_on_grid(:) = ngwfs_on_grid(:) + &
          !          buffer_on_grid(:) * coeffs(orig_atom)
          ! agrecocmplx: use new routine in basis_new_mod
          call data_functions_axpy(ngwfs_on_grid, buffer_on_grid, &
                                   coeffs(orig_atom))
       else
          !ngwfs_on_grid(:) = ngwfs_on_grid(:) + buffer_on_grid(:)
          ! agrecocmplx: use new routine in basis_new_mod
          call data_functions_axpy(ngwfs_on_grid, buffer_on_grid, 1.0_DP)
       endif

    enddo atom_loop

    ! cks: don't go away just yet! Wait for all communications to finish!
    call comms_barrier
    ! cks: clean up send stack too
    call comms_free

    !deallocate(buffer_on_grid, stat=ierr)
    !call utils_dealloc_check('elec_history_read_ngwfs', &
    !     'buffer_on_grid',ierr)
    ! agrecocmplx: deallocate using appropriate routine
    call data_functions_dealloc(buffer_on_grid)
    !deallocate(fftbox_buffer,stat=ierr)
    !call utils_dealloc_check('elec_history_read_ngwfs', &
    !     'fftbox_buffer',ierr)
    ! agrecocmplx: deallocate using appropriate routine
    call data_fftbox_dealloc(fftbox_buffer)
    deallocate(fftbox_complex_shifted,stat=ierr)
    call utils_dealloc_check('elec_history_read_ngwfs', &
         'fftbox_complex_shifted',ierr)
    deallocate(fftbox_complex,stat=ierr)
    call utils_dealloc_check('elec_history_read_ngwfs', &
         'fftbox_complex',ierr)
    !deallocate(uni_tbox,stat=ierr)
    !call utils_dealloc_check('elec_history_read_ngwfs', &
    !     'uni_tbox',ierr)
    ! agrecocmplx: deallocate using appropriate routine
    call data_fftbox_dealloc(uni_tbox)

    ! ndmh: Clean up temporary basis if we created one, and reset coin1,2,3
    if ((orig_coin1.neqv.fftbox%coin1).or. &
         (orig_coin2.neqv.fftbox%coin2).or. &
         (orig_coin3.neqv.fftbox%coin3)) then

       call function_basis_deallocate(tmp_basis)
       fftbox%coin1 = orig_coin1
       fftbox%coin2 = orig_coin2
       fftbox%coin3 = orig_coin3

    end if

    if (pub_debug_on_root) write(stdout,'(a)') &
       'DEBUG: Leaving elec_history_read_ngwfs'

    return

  end subroutine elec_history_read_ngwfs


!============================================================================!
!============================================================================!
!============================================================================!


  subroutine elec_history_compose_ngwfs(ngwfs_on_grid, ngwf_basis, &
             cell, fftbox, elements, ireg, cstatus, par)

    !========================================================================!
    ! This subroutine creates new NGWFs from data kept in store. New ngwfs   !
    ! are in "universal tightbox" representation from previous stored.       !
    !                                                                        !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 20/06/2011                            !
    ! Modified by Andrea Greco on 05/16/2015 to allow use of complex NGWFs.  !
    !========================================================================!

    ! agrecocmplx
    use comms, only: pub_on_root, pub_my_proc_id
    use constants, only: DP, stdout
    use datatypes, only: FUNCTIONS, COEF, data_functions_alloc, &
         data_functions_dealloc, data_functions_axpy, data_functions_dot, &
         data_set_to_zero, data_functions_mix, data_functions_copy
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use ion, only : ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use services, only : services_polynomial_step, services_rms_fit
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    !type(PARAL_INFO), intent(in)  :: par
    type(ELEMENT), intent(in)        :: elements(:) ! size par%nat)
    type(FUNC_BASIS), intent(in)     :: ngwf_basis
    type(CELL_INFO), intent(in)      :: cell
    type(FFTBOX_INFO), intent(inout) :: fftbox
    type(FUNCTIONS), intent(inout)   :: ngwfs_on_grid
    integer, intent(out)             :: cstatus
    type(PARAL_INFO), intent(in), pointer     :: par
    integer, intent(in)              :: ireg  ! jcap: region counter

    ! Local variables
    integer  :: label
    integer  :: ierr, iom

    ! Normalisation variables
    integer       :: row, start, finish
    real(kind=DP) :: norm_constant

    ! Mixing variables
    type(FUNCTIONS) :: mix_ngwfs_on_grid
    ! agrecocmplx: local variable for complex NGWFs
    logical :: loc_cmplx
    ! agrecocmplx
    type(COEF) :: norm_temp

    !======================================================================!

    if (pub_debug_on_root) write(stdout,'(a,i4,a,i4,a,i4,a)') &
        'DEBUG: Entering elec_history_compose_ngwfs (hio =', history_io, &
        '), (mio =',methods_io,'), (hcount =', history_count,')'

    ! Exit if stores haven't been initialized
    if (.not.allocated(history) .or. .not.allocated(methods) .or. &
        history_io==0 .or. methods_io==0 .or. history_count==0) then
       cstatus = -1
       if (pub_debug_on_root) write(stdout,'(a)') &
            'DEBUG: Leaving elec_history_compose_ngwfs'
       return
    endif

    ! Determine composition method
    iom = methods_io
    cstatus = 1

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid%iscmplx

    ! If required keep the input NGWFs into buffer
    if (methods(iom)%ngwfs_coeff .gt. 0.0) then
       !allocate(mix_ngwfs_on_grid(ngwf_basis%size_on_grid), stat=ierr)
       !call utils_alloc_check('elec_history_compose_ngwfs','mix_ngwfs_on_grid',ierr)
       !mix_ngwfs_on_grid = ngwfs_on_grid
       ! agrecocmplx: use appropriate routine to allocate functions on grid
       call data_functions_alloc(mix_ngwfs_on_grid, ngwf_basis%size_on_grid, &
           iscmplx=loc_cmplx)
       call data_functions_copy(mix_ngwfs_on_grid, ngwfs_on_grid)
    endif

    !=============================================================!
    ! New entry
    !=============================================================!
    if (methods(iom)%ngwfs_mix_type == 'NONE' &
            .or.   methods(iom)%ngwfs_count == 0) then

       cstatus = 1

    !=============================================================!
    ! Initialisation phase
    !=============================================================!
    elseif (methods(iom)%ngwfs_count .le. methods(iom)%ngwfs_init_num) then

       call internal_init()

    !=============================================================!
    ! No mixing : new NGWFs = last NGWFs stored
    !=============================================================!
    elseif (methods(iom)%ngwfs_mix_type == 'REUSE') then

       call internal_retrieve_last()

    !=============================================================!
    ! Linear and polynomial extrapolations
    !=============================================================!
    elseif (methods(iom)%ngwfs_mix_type == 'LINEAR') then
       call internal_linear_xtpol()

    elseif (methods(iom)%ngwfs_mix_type == 'MULTID') then
       call internal_arias_xtpol()

    elseif (methods(iom)%ngwfs_mix_type == 'POLY') then
       call internal_poly_xtpol()

    elseif (methods(iom)%ngwfs_mix_type == 'LOCAL') then
       call internal_generalized_arias_xtpol()

    !=============================================================!
    ! Propagation of auxiliary support functions :
    ! algorithm based on A. Niklasson et al,
    !           J. Chem. Phys. 130, 214109 (2009)
    !=============================================================!
    elseif (methods(iom)%ngwfs_mix_type == 'TRPROP') then
       call internal_tr_propagation()

    endif

    !=============================================================!
    ! Precondition new NGWFs
    !=============================================================!

    if (methods(iom)%ngwfs_coeff .gt. 0.0) then

       if (cstatus == 0 ) then
          if (pub_on_root) write(stdout,'(a)') &
             'Electronic history : preconditioning of new NGWFs '
          !ngwfs_on_grid(:) = (1.0-methods(iom)%ngwfs_coeff)*ngwfs_on_grid(:) &
          !    + methods(iom)%ngwfs_coeff*mix_ngwfs_on_grid(:)
          ! agrecocmplx
          call data_functions_mix(ngwfs_on_grid, mix_ngwfs_on_grid, &
               (1.0_DP-methods(iom)%ngwfs_coeff))
       endif

       !deallocate(mix_ngwfs_on_grid, stat=ierr)
       !call utils_dealloc_check('elec_history_compose_ngwfs','mix_ngwfs_on_grid',ierr)
       ! agrecocmplx: deallocate with appropriate routine
       call data_functions_dealloc(mix_ngwfs_on_grid)
    endif


    !=============================================================!
    ! Check the normalization of new NGWFs
    !=============================================================!

    if (cstatus == 0) then
       ! loop over ngwfs on this proc
       do row=1,ngwf_basis%num_on_proc(pub_my_proc_id)

          start=ngwf_basis%spheres(row)%offset
          finish=start+(ngwf_basis%spheres(row)%n_ppds_sphere)*cell%n_pts -1

          !norm_constant=sum(ngwfs_on_grid(start:finish)*ngwfs_on_grid(start:finish) )
          ! agrecocmplx: use data_functions_dot
          norm_temp = data_functions_dot(ngwfs_on_grid, ngwfs_on_grid, &
                          start, start, finish - start + 1)
          !norm_constant=sqrt(norm_constant*cell%weight)
          ! agrecocmplx
          if (norm_temp%iscmplx) then
              norm_constant = sqrt(real(norm_temp%z, kind=DP)) !jmecmplx
          else
              norm_constant = sqrt(norm_temp%d)
          end if
          if (norm_constant.gt.0) then
             !ngwfs_on_grid(start:finish)=ngwfs_on_grid(start:finish)/norm_constant
             ! agrecocmplx
             if (ngwfs_on_grid%iscmplx) then
                 ngwfs_on_grid%z(start:finish)=ngwfs_on_grid%z(start:finish)/norm_constant
             else
                 ngwfs_on_grid%d(start:finish)=ngwfs_on_grid%d(start:finish)/norm_constant
             end if
          endif

       enddo
    endif


    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving elec_history_compose_ngwfs'

  contains

    !################################################################!

    subroutine internal_init()

       use sparse, only: SPAM3, sparse_create, sparse_destroy

       implicit none

       ! Local arguments
       type(SPAM3) :: lwd_overlap

       if (methods(iom)%ngwfs_mix_type == 'TRPROP') then

          lwd_overlap%structure = 'K'
          call sparse_create(lwd_overlap)

          label = methods(iom)%ptr(1)
          call elec_history_read_ngwfs(ngwfs_on_grid,ngwf_basis, &
               cell,fftbox,elements,ireg,'scf',label,par=par)
          call elec_history_store_ngwfs(ngwfs_on_grid,ngwf_basis, &
               cell,fftbox,elements,ireg,'init',par=par)

          if (methods(iom)%ngwfs_init_type == 'NONE') then
             !ngwfs_on_grid(:) = mix_ngwfs_on_grid(:)
             ! agrecocmplx
             call data_functions_copy(ngwfs_on_grid, mix_ngwfs_on_grid)
             cstatus = 1

          elseif (methods(iom)%ngwfs_init_type == 'REUSE') then
             cstatus = 0
             if (pub_on_root) write(stdout,'(a)') &
                'Electronic history : new NGWFs = last NGWFs stored '
          endif

          call sparse_destroy(lwd_overlap)

       else
          if (methods(iom)%ngwfs_init_type == 'NONE') then
             !ngwfs_on_grid(:) = mix_ngwfs_on_grid(:)
             ! agrecocmplx
             call data_functions_copy(ngwfs_on_grid, mix_ngwfs_on_grid)
             cstatus = 1

          elseif (methods(iom)%ngwfs_init_type == 'REUSE') then
             label = methods(iom)%ptr(1)
             call elec_history_read_ngwfs(ngwfs_on_grid,ngwf_basis, &
                  cell,fftbox,elements,ireg,'scf',label,par=par)
             cstatus = 0
             if (pub_on_root) write(stdout,'(a)') &
                'Electronic history : new NGWFs = last NGWFs stored '

          endif
       endif

    end subroutine internal_init

    !################################################################!

    subroutine internal_retrieve_last()

       !=============================================================!
       ! No mixing : new density kernel = last density kernel stored
       !=============================================================!

       implicit none

       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : new NGWFs = last NGWFs stored '
       label = methods(iom)%ptr(1)
       call elec_history_read_ngwfs(ngwfs_on_grid,ngwf_basis, &
            cell,fftbox,elements,ireg,'scf',label,par=par)
       cstatus = 0

    end subroutine internal_retrieve_last

    !################################################################!

    subroutine internal_linear_xtpol()

       !=============================================================!
       ! Linear extrapolation : NGWFs(new) = 2*NGWFs(-1) - NGWFs(-2) !
       ! Modified by Andrea Greco on 16/05/2015 for complex NGWFs.   !
       !=============================================================!

       implicit none

       type(FUNCTIONS) :: old_ngwfs_on_grid
       type(FUNCTIONS) :: tmp_ngwfs_on_grid

       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : linear extrapolation of NGWFs '

       !allocate(tmp_ngwfs_on_grid(ngwf_basis%size_on_grid), stat=ierr)
       !call utils_alloc_check('elec_history_compose_ngwfs','tmp_ngwfs_on_grid',ierr)
       ! agrecocmplx: allocate using appropriate routine
       call data_functions_alloc(tmp_ngwfs_on_grid, ngwf_basis%size_on_grid, &
                                  iscmplx=ngwfs_on_grid%iscmplx)
       !allocate(old_ngwfs_on_grid(ngwf_basis%size_on_grid), stat=ierr)
       !call utils_alloc_check('elec_history_compose_ngwfs','old_ngwfs_on_grid',ierr)
       ! agrecocmplx: allocate using appropriate routine
       call data_functions_alloc(old_ngwfs_on_grid, ngwf_basis%size_on_grid, &
                                  iscmplx=ngwfs_on_grid%iscmplx)

       label = methods(iom)%ptr(1)
       call elec_history_read_ngwfs(old_ngwfs_on_grid,ngwf_basis, &
            cell,fftbox,elements,ireg,'scf',label,par=par)
       !tmp_ngwfs_on_grid = old_ngwfs_on_grid * 2.0_dp
       ! agrecocmplx
       if (ngwfs_on_grid%iscmplx) then
           tmp_ngwfs_on_grid%z(:) = old_ngwfs_on_grid%z * 2.0_DP
       else
           tmp_ngwfs_on_grid%d(:) = old_ngwfs_on_grid%d * 2.0_DP
       end if

       label = methods(iom)%ptr(2)
       call elec_history_read_ngwfs(old_ngwfs_on_grid,ngwf_basis, &
            cell,fftbox,elements,ireg,'scf',label,par=par)
       !ngwfs_on_grid = tmp_ngwfs_on_grid - old_ngwfs_on_grid
       ! agrecocmplx
       if (ngwfs_on_grid%iscmplx) then
           ngwfs_on_grid%z(:) = tmp_ngwfs_on_grid%z - old_ngwfs_on_grid%z
       else
           ngwfs_on_grid%d(:) = tmp_ngwfs_on_grid%d - old_ngwfs_on_grid%d
       end if
       cstatus = 0

       !deallocate(old_ngwfs_on_grid, stat=ierr)
       !call utils_dealloc_check('elec_history_compose_ngwfs','old_ngwfs_on_grid',ierr)
       ! agrecocmplx
       call data_functions_dealloc(old_ngwfs_on_grid)
       !deallocate(tmp_ngwfs_on_grid, stat=ierr)
       !call utils_dealloc_check('elec_history_compose_ngwfs','tmp_ngwfs_on_grid',ierr)
       ! agrecocmplx
       call data_functions_dealloc(tmp_ngwfs_on_grid)

    end subroutine internal_linear_xtpol

    !################################################################!

    subroutine internal_arias_xtpol()

       !=============================================================!
       ! Linear extrapolation :
       ! algorithm based on Arias et al, PRB 45, 1538 (1992)         !
       ! Modified by Andrea Greco to allow use of complex NGWFs.     !
       !=============================================================!

       implicit none

       real(kind=DP), allocatable  :: linear_coeffs(:)
       real(kind=DP), allocatable  :: old_coords(:,:,:)
       real(kind=DP), allocatable  :: new_coords(:,:)
       type(FUNCTIONS) :: old_ngwfs_on_grid
       integer  :: ii, istep


       if (pub_on_root) write(stdout,'(a,i2,a)') &
          'Electronic history : multi-dimensional extrapolation of NGWFs '

       ! Mixing coefficients
       allocate(linear_coeffs(methods(iom)%ngwfs_mix_size),stat=ierr)
       call utils_alloc_check('elec_history_compose_ngwfs','linear_coeffs',ierr)

       ! Set of coordinates to be interpolated
       allocate(new_coords(3,size(elements)),stat=ierr)
       call utils_alloc_check('elec_history_compose_ngwfs','new_coords',ierr)
       new_coords = history(history_io)%coord

       allocate(old_coords(3,size(elements),methods(iom)%ngwfs_mix_size),stat=ierr)
       call utils_alloc_check('elec_history_compose_ngwfs','old_coords',ierr)
       do istep = 1, methods(iom)%ngwfs_mix_size
          ! Debug
          ii = history_index(label=methods(iom)%ptr(istep))
          old_coords(:,:,istep) = history(ii)%coord
          !old_coords(:,:,istep) = history(methods(iom)%ptr(istep))%coord
       enddo

       ! Interpolation coefficients
       call services_rms_fit(linear_coeffs,methods(iom)%ngwfs_mix_size,3, &
               size(elements), old_coords,new_coords)

       deallocate(old_coords,stat=ierr)
       call utils_dealloc_check('elec_history_compose_ngwfs','old_coords',ierr)
       deallocate(new_coords,stat=ierr)
       call utils_dealloc_check('elec_history_compose_ngwfs','new_coords',ierr)

       ! Linear combination of NGWFs
       !allocate(old_ngwfs_on_grid(ngwf_basis%size_on_grid), stat=ierr)
       !call utils_alloc_check('elec_history_compose_ngwfs','old_ngwfs_on_grid',ierr)
       ! agrecocmplx: allocate using appropriate routine
       call data_functions_alloc(old_ngwfs_on_grid, ngwf_basis%size_on_grid, &
                                  iscmplx=ngwfs_on_grid%iscmplx)
       !ngwfs_on_grid(:) = 0.0_dp
       ! agrecocmplx
       call data_set_to_zero(ngwfs_on_grid)

       do istep = 1, methods(iom)%ngwfs_mix_size
          label = methods(iom)%ptr(istep)
          !old_ngwfs_on_grid(:) = 0.0_dp
          ! agrecocmplx
          call data_set_to_zero(old_ngwfs_on_grid)
          call elec_history_read_ngwfs(old_ngwfs_on_grid,ngwf_basis, &
               cell,fftbox,elements,ireg,'scf',label,par=par)
          !ngwfs_on_grid(:) = ngwfs_on_grid(:) + &
          !   linear_coeffs(istep)*old_ngwfs_on_grid(:)
          ! agrecocmplx: use new data_functions_axpy routine
          call data_functions_axpy(ngwfs_on_grid, old_ngwfs_on_grid, &
                                    linear_coeffs(istep))
       enddo
       cstatus = 0

       deallocate(linear_coeffs, stat=ierr)
       call utils_dealloc_check('elec_history_compose_ngwfs','linear_coeffs',ierr)
       !deallocate(old_ngwfs_on_grid, stat=ierr)
       !call utils_dealloc_check('elec_history_compose_ngwfs','old_ngwfs_on_grid',ierr)
       ! agrecocmplx: deallocate using appropriate routine
       call data_functions_dealloc(old_ngwfs_on_grid)
    end subroutine internal_arias_xtpol

    !################################################################!

    subroutine internal_generalized_arias_xtpol()

       !=============================================================!
       ! Generalized linear mixing :
       ! for each atom, a different set of coefficients is computed
       ! as a function of its local environemet
       ! Modified by Andrea Greco to use complex NGWFs.              !
       !=============================================================!

       implicit none

       integer  :: orig_atom
       integer  :: istep, ii
       real(kind=DP), allocatable  :: old_coords(:,:,:)
       real(kind=DP), allocatable  :: new_coords(:,:)
       real(kind=DP), allocatable  :: generalized_coeffs(:,:)
       type(FUNCTIONS) :: old_ngwfs_on_grid

       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : local multi-dimensional extrapolation of NGWFs '

       ! Mixing coefficients
       allocate(generalized_coeffs(methods(iom)%ngwfs_mix_size,size(elements)),stat=ierr)
       call utils_alloc_check('elec_history_compose_ngwfs','generalized_coeffs',ierr)

       ! Set of coordinates to be interpolated
       allocate(new_coords(3,size(elements)),stat=ierr)
       call utils_alloc_check('elec_history_compose_ngwfs','new_coords',ierr)
       new_coords = history(history_io)%coord

       allocate(old_coords(3,size(elements),methods(iom)%ngwfs_mix_size),stat=ierr)
       call utils_alloc_check('elec_history_compose_ngwfs','old_coords',ierr)
       do istep = 1, methods(iom)%ngwfs_mix_size
          ! Debug
          ii = history_index(label=methods(iom)%ptr(istep))
          old_coords(:,:,istep) = history(ii)%coord
          !old_coords(:,:,istep) = history(methods(iom)%ptr(istep))%coord
       enddo

       ! Interpolation coefficients
       do orig_atom = 1, size(elements)

          ! Compute weigthed atomic coordinates
          call generalized_coordinates(methods(iom)%ngwfs_mix_size, &
               cell,old_coords,new_coords, orig_atom)

          ! Interpolation coefficients
          call services_rms_fit(generalized_coeffs(:,orig_atom), &
                  methods(iom)%ngwfs_mix_size, &
                  3,size(elements),old_coords,new_coords)

       enddo

       deallocate(old_coords, stat=ierr)
       call utils_dealloc_check('elec_history_compose_ngwfs','old_coords',ierr)
       deallocate(new_coords, stat=ierr)
       call utils_dealloc_check('elec_history_compose_ngwfs','coord',ierr)

       ! Linear combination of NGWFs
       !allocate(old_ngwfs_on_grid(ngwf_basis%size_on_grid), stat=ierr)
       !call utils_alloc_check('elec_history_compose_ngwfs','old_ngwfs_on_grid',ierr)
       ! agrecocmplx: allocate using appropriate routines
       call data_functions_alloc(old_ngwfs_on_grid, ngwf_basis%size_on_grid, &
                                  iscmplx=ngwfs_on_grid%iscmplx)
       !ngwfs_on_grid(:) = 0.0_dp
       ! agrecocmplx
       call data_set_to_zero(ngwfs_on_grid)
       do istep = 1, methods(iom)%ngwfs_mix_size
          label = methods(iom)%ptr(istep)
          !old_ngwfs_on_grid(:) = 0.0_dp
          ! agrecocmplx
          call data_set_to_zero(old_ngwfs_on_grid)
          call elec_history_read_ngwfs(old_ngwfs_on_grid,ngwf_basis, &
               cell,fftbox,elements,ireg,'scf',label, &
               coeffs=generalized_coeffs(istep,:),par=par)
          !ngwfs_on_grid(:) = ngwfs_on_grid(:) + old_ngwfs_on_grid(:)
          ! agrecocmplx
          call data_functions_axpy(ngwfs_on_grid, old_ngwfs_on_grid, 1.0_DP)
       end do
       cstatus = 0

       !deallocate(old_ngwfs_on_grid, stat=ierr)
       !call utils_dealloc_check('elec_history_compose_ngwfs','buffer_on_grid',ierr)
       ! agrecocmplx: deallocate with appropriate routines
       call data_functions_dealloc(old_ngwfs_on_grid)

       deallocate(generalized_coeffs, stat=ierr)
       call utils_dealloc_check('elec_history_compose_ngwfs','generalized_coeffs',ierr)

    end subroutine internal_generalized_arias_xtpol

    !################################################################!

    subroutine internal_poly_xtpol()

       !=============================================================!
       ! Polynomial extrapolation in real space                      !
       ! Modified by Andrea Greco to allow use of complex NGWFs.     !
       !=============================================================!

       use constants, only: cmplx_0
       use services, only : services_polynomial_step

       implicit none

       real(kind=DP) :: xtpol_step
       integer  :: istep
       real(kind=DP), allocatable  :: time_steps(:)
       type(FUNCTIONS) :: old_ngwfs_on_grid
       real(kind=DP), allocatable  :: ngwfs_history_d(:,:)
       complex(kind=DP), allocatable  :: ngwfs_history_z(:,:)

       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : polynomial extrapolation of NGWFs '

       allocate(time_steps(methods(iom)%ngwfs_mix_size), stat=ierr)
       call utils_alloc_check('elec_history_compose_ngwfs','time_steps',ierr)
       !allocate(old_ngwfs_on_grid(ngwf_basis%size_on_grid), stat=ierr)
       !call utils_alloc_check('elec_history_compose_ngwfs','buffer_on_grid',ierr)
       ! agrecocmplx: allocate using appropriate routine
       call data_functions_alloc(old_ngwfs_on_grid, ngwf_basis%size_on_grid, &
                                  iscmplx=ngwfs_on_grid%iscmplx)
       !ngwfs_on_grid(:) = 0.0_dp
       ! agrecocmplx
       call data_set_to_zero(ngwfs_on_grid)

       if (ngwfs_on_grid%iscmplx) then
          allocate(ngwfs_history_z(methods(iom)%ngwfs_mix_size, &
               ngwf_basis%size_on_grid), stat=ierr)
          call utils_alloc_check('elec_history_compose_ngwfs', &
               'ngwfs_history_z', ierr)
          ngwfs_history_z(:,:) = cmplx_0
       else
          allocate(ngwfs_history_d(methods(iom)%ngwfs_mix_size, &
               ngwf_basis%size_on_grid), stat=ierr)
          call utils_alloc_check('elec_history_compose_ngwfs', &
               'ngwfs_history_d', ierr)
          ngwfs_history_d(:,:) = 0.0_DP
       end if

       ! Loop over the ngwfs to be mixed
       do istep = 1,methods(iom)%ngwfs_mix_size
          label = methods(iom)%ptr(istep)
          !old_ngwfs_on_grid(:) = 0.0_dp
          ! agrecocmplx
          call data_set_to_zero(old_ngwfs_on_grid)
          call elec_history_read_ngwfs(old_ngwfs_on_grid,ngwf_basis, &
               cell,fftbox,elements,ireg,'scf',label,par=par)
          ! agrecocmplx
          if (ngwfs_on_grid%iscmplx) then
              ngwfs_history_z(istep,:) = old_ngwfs_on_grid%z
          else
              ngwfs_history_d(istep,:) = old_ngwfs_on_grid%d
          end if
          time_steps(istep) = real(istep)
       enddo

       xtpol_step = real(methods(iom)%ngwfs_mix_size+1)
       if (ngwfs_on_grid%iscmplx) then
          call services_polynomial_step(ngwfs_on_grid%z, &
               methods(iom)%ngwfs_mix_size, ngwf_basis%size_on_grid,&
               ngwfs_history_z,time_steps,xtpol_step)
       else
          call services_polynomial_step(ngwfs_on_grid%d, &
               methods(iom)%ngwfs_mix_size, ngwf_basis%size_on_grid,&
               ngwfs_history_d,time_steps,xtpol_step)
       end if
       cstatus = 0

       ! agrecocmplx: deallocate with appropriate routines

       if (ngwfs_on_grid%iscmplx) then
          deallocate(ngwfs_history_z, stat=ierr)
          call utils_dealloc_check('elec_history_compose_ngwfs', &
               'ngwfs_history_z', ierr)
       else
          deallocate(ngwfs_history_d, stat=ierr)
          call utils_dealloc_check('elec_history_compose_ngwfs', &
               'ngwfs_history_d', ierr)
       end if
       !deallocate(old_ngwfs_on_grid, stat=ierr)
       !call utils_dealloc_check('elec_history_compose_ngwfs','old_ngwfs_on_grid',ierr)
       ! agrecocmplx: deallocate with appropriate routines
       call data_functions_dealloc(old_ngwfs_on_grid)
       deallocate(time_steps, stat=ierr)
       call utils_dealloc_check('elec_history_compose_ngwfs','time_steps',ierr)

    end subroutine internal_poly_xtpol

    !################################################################!

    subroutine internal_tr_propagation()

       !=============================================================!
       ! Propagation of auxiliary support functions :
       ! algorithm based on A. Niklasson et al,
       !           J. Chem. Phys. 130, 214109 (2009)
       !           section IIC, case K = 0
       ! Modified by Andrea Greco to use complex NGWFs.              !
       !=============================================================!

       ! kappa_mod
       use rundat, only: pub_devel_code
       use comms, only: comms_bcast
       use sparse, only: SPAM3, sparse_create, sparse_destroy
       use utils, only: utils_devel_code

       implicit none

       type(FUNCTIONS) :: aux_ngwfs_on_grid
       type(FUNCTIONS) :: old_aux_ngwfs_on_grid
       ! kappa_mod
       real(kind=DP)      :: kappa
!       character(len=200) :: devel_code
!       integer            :: start_pos, stop_pos, test_pos


       if (pub_on_root) write(stdout,'(a)') &
          'Electronic history : time-reversible propagation of auxiliary NGWFs'

       !allocate(aux_ngwfs_on_grid(ngwf_basis%size_on_grid), stat=ierr)
       !call utils_alloc_check('elec_history_compose_ngwfs','aux_ngwfs_on_grid',ierr)
       !allocate(old_aux_ngwfs_on_grid(ngwf_basis%size_on_grid), stat=ierr)
       !call utils_alloc_check('elec_history_compose_ngwfs','old_aux_ngwfs_on_grid',ierr)
       ! agrecocmplx: allocate with appropriate routines
       call data_functions_alloc(aux_ngwfs_on_grid, ngwf_basis%size_on_grid, &
                                  iscmplx=ngwfs_on_grid%iscmplx)
       call data_functions_alloc(old_aux_ngwfs_on_grid, ngwf_basis%size_on_grid, &
                                  iscmplx=ngwfs_on_grid%iscmplx)

       !- kappa_mod -------------------------------------------------!

       kappa = 2.0_dp

       ! Check local copy of devel_code string and act accordingly
       !if (pub_on_root) then
       !   devel_code=pub_devel_code
       !   if (pub_on_root) write(stdout,*) 'ehist: devel_code ='
       !   if (len_trim(devel_code)>0) then
       !      start_pos=index(devel_code,'ehist:')
       !      stop_pos=index(devel_code,':ehist')
       !      if (stop_pos<=0) stop_pos=len_trim(devel_code)
       !      if (start_pos>0) then
       !         test_pos=index(devel_code,'kappa=')
       !         if (test_pos>start_pos.and.test_pos<stop_pos) then
       !            test_pos=test_pos+len('kappa=')
       !            read(devel_code(test_pos:test_pos+ &
       !                 & index(devel_code(test_pos:stop_pos),':')-2),*) kappa
       !         end if
       !         if (pub_on_root) then
       !            write(stdout,'(/a)') 'ehist: Processing of the devel_code string...'
       !            write(stdout,'(a,i2)') 'ehist: kappa = ', kappa
       !         end if
       !     end if
       !   end if
       !end if
       !call comms_bcast(pub_root_proc_id,kappa)
       kappa = utils_devel_code(2.0_dp,'EHIST','KAPPA',pub_devel_code)
       if (pub_on_root) then
          write(stdout,'(a,i2)') 'ehist: kappa = ', kappa
       endif
       !-------------------------------------------------------------!

       ! Time reversible propagation of NGWFs_aux
       label = methods(iom)%ptr(1)
       call elec_history_read_ngwfs(ngwfs_on_grid,ngwf_basis, &
            cell,fftbox,elements,ireg,'scf',label,par=par)

       label = methods(iom)%ptr(1)
       call elec_history_read_ngwfs(aux_ngwfs_on_grid,ngwf_basis, &
            cell,fftbox,elements,ireg,'init',label,par=par)

       label = methods(iom)%ptr(2)
       call elec_history_read_ngwfs(old_aux_ngwfs_on_grid,ngwf_basis, &
            cell,fftbox,elements,ireg,'init',label,par=par)

       !ngwfs_on_grid = kappa*ngwfs_on_grid + (2.0_dp-kappa)*aux_ngwfs_on_grid &
       !                - old_aux_ngwfs_on_grid
       ! agrecocmplx
       if (ngwfs_on_grid%iscmplx) then
           ngwfs_on_grid%z(:) = kappa*ngwfs_on_grid%z + (2.0_dp-kappa)*aux_ngwfs_on_grid%z &
                           - old_aux_ngwfs_on_grid%z
       else
           ngwfs_on_grid%d(:) = kappa*ngwfs_on_grid%d + (2.0_dp-kappa)*aux_ngwfs_on_grid%d &
                           - old_aux_ngwfs_on_grid%d
       end if
       call elec_history_store_ngwfs(ngwfs_on_grid, ngwf_basis, &
            cell,fftbox,elements,ireg,'init', par=par)
       cstatus = 0

       !deallocate(aux_ngwfs_on_grid, stat=ierr)
       !call utils_dealloc_check('elec_history_compose_ngwfs','aux_ngwfs_on_grid',ierr)
       !deallocate(old_aux_ngwfs_on_grid, stat=ierr)
       !call utils_dealloc_check('elec_history_compose_ngwfs','old_aux_ngwfs_on_grid',ierr)
       ! agrecocmplx: deallocate with appropriate routines
       call data_functions_dealloc(aux_ngwfs_on_grid)
       call data_functions_dealloc(old_aux_ngwfs_on_grid)

    end subroutine internal_tr_propagation

    !################################################################!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! vv: Extended-Lagrangian methods for the propagation of NGWFs are not yet
! available in ONETEP. Further developments need to be done before these schemes
! can be part of the devel functionalities. I keep the implementation
! I have made for the orthogonal (Lowdin transformed) WFs and the
! conventional NGWFs below for future reference.
! However, I have deleted the calls to these subroutines from other modules,
! i.e. ngwf_cg_mod and I have also deleted the related lines within this module.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!   subroutine internal_dissip_ortho_propagation()
!!
!!      !=============================================================!
!!      ! Propagation of auxiliary support functions :
!!      ! algorithm based on A. Niklasson et al,
!!      !           J. Chem. Phys. 130, 214109 (2009)
!!      !           section IV, case 3 <= K <= 9
!!      ! Modified by Andrea Greco to allow use of complex NGWFs
!!      !=============================================================!
!!
!!      use basis, only: basis_clean_function
!!      use datatypes, only: data_set_to_zero, data_functions_copy
!!      use function_ops, only: function_ops_sum_ppd_funcs
!!      use sparse, only: SPAM3, sparse_create, sparse_destroy
!!
!!      implicit none
!!
!!      ! vv: Local variables
!!      type(FUNCTIONS) :: old_ogwfs_on_grid, ogwfs_on_grid ,new_ogwfs_on_grid(1)
!!      integer     :: icc
!!      integer  :: istep, ingwf
!!      integer       :: xl_max = 10
!!      real(kind=DP), parameter, dimension(7)   :: xl_k = (/  &
!!          1.690_dp, 1.750_dp, 1.820_dp, 1.840_dp, 1.860_dp, 1.880_dp, 1.890_dp /)
!!      real(kind=DP), parameter, dimension(70)  :: xl_coef = (/  &
!!      ! vv : Set of correct coefficients
!!            0.01000_dp, -0.55000_dp,  0.00000_dp, -0.15000_dp,  0.00000_dp,  0.00000_dp,&
!!            0.00000_dp,  0.00000_dp,  0.00000_dp,  0.00000_dp,&
!!            0.07900_dp, -0.65800_dp, -0.11400_dp, -0.11400_dp,  0.05700_dp,  0.00000_dp,&
!!            0.00000_dp,  0.00000_dp,  0.00000_dp,  0.00000_dp,&
!!            0.07200_dp, -0.74800_dp, -0.14400_dp, -0.05400_dp,  0.07200_dp, -0.01800_dp,&
!!            0.00000_dp,  0.00000_dp,  0.00000_dp,  0.00000_dp,&
!!            0.08300_dp, -0.80200_dp, -0.14850_dp, -0.01100_dp,  0.06600_dp, -0.03300_dp,&
!!            0.00550_dp,  0.00000_dp,  0.00000_dp,  0.00000_dp,&
!!            0.08240_dp, -0.84160_dp, -0.14080_dp,  0.01760_dp,  0.05120_dp, -0.04000_dp,&
!!            0.01280_dp, -0.00160_dp,  0.00000_dp,  0.00000_dp,&
!!            0.07644_dp, -0.87416_dp, -0.12584_dp,  0.03432_dp,  0.03432_dp, -0.03960_dp,&
!!            0.01848_dp, -0.00440_dp,  0.00044_dp,  0.00000_dp,&
!!            0.07568_dp, -0.89704_dp, -0.11232_dp,  0.04368_dp,  0.02016_dp, -0.03600_dp,&
!!            0.02208_dp, -0.00756_dp,  0.00144_dp, -0.00012_dp /)
!!      type(SPAM3)  :: Umat(1)
!!
!!      if (pub_on_root) write(stdout,'(a)') &
!!         'Electronic history : dissipative propagation of auxiliary NGWFs'
!!
!!      !allocate(old_ngwfs_on_grid(ngwf_basis%size_on_grid), stat=ierr)
!!      !call utils_alloc_check('elec_history_compose_ngwfs','old_ngwfs_on_grid',ierr)
!!      ! agrecocmplx: allocate with appropriate routines
!!      call data_functions_alloc(old_ogwfs_on_grid, ngwf_basis%size_on_grid, &
!!                                 iscmplx=ngwfs_on_grid%iscmplx)
!!      call data_functions_alloc(new_ogwfs_on_grid(1), ngwf_basis%size_on_grid, &
!!                                 iscmplx=ngwfs_on_grid%iscmplx)
!!      call data_functions_alloc(ogwfs_on_grid, ngwf_basis%size_on_grid, &
!!                                 iscmplx=ngwfs_on_grid%iscmplx)
!!
!!      ! Dissipative propagation of NGWFs_aux
!!      label = methods(iom)%ptr(1)
!!      call elec_history_read_ngwfs(old_ogwfs_on_grid,ngwf_basis, &
!!           cell,fftbox,elements,'scf',label)
!!
!!      ! vv: Phase alignement
!!      !=======================================================================
!!       Umat(1)%structure = 'K'
!!       call sparse_create(Umat(1))
!!       call elec_history_read_ngwfs(ogwfs_on_grid,ngwf_basis,&
!!            cell,fftbox,elements,'init',label)
!!       call data_set_to_zero(new_ogwfs_on_grid(1))
!!       call internal_unit_transf(Umat(1),old_ogwfs_on_grid,ogwfs_on_grid)
!!       call function_ops_sum_ppd_funcs(new_ogwfs_on_grid(1),ngwf_basis,&
!!            Umat(1), 1, 1, Umat(1), old_ogwfs_on_grid, ngwf_basis)
!!       call data_set_to_zero(old_ogwfs_on_grid)
!!       call data_functions_copy(old_ogwfs_on_grid,new_ogwfs_on_grid(1))
!!       call sparse_destroy(Umat(1))
!!      !=======================================================================
!!
!!      call data_set_to_zero(ogwfs_on_grid)
!!      icc = methods(iom)%ngwfs_mix_size - 3
!!      ! agrecocmplx
!!      call data_functions_axpy(ogwfs_on_grid,old_ogwfs_on_grid,&
!!           xl_k(icc))
!!
!!      do istep = 1, methods(iom)%ngwfs_mix_size
!!         label = methods(iom)%ptr(istep)
!!         ! vv: Call ogwfs from memory
!!         call elec_history_read_ngwfs(old_ogwfs_on_grid,ngwf_basis, &
!!              cell,fftbox,elements,'init',label)
!!
!!         ! vv: Propagate the orthogonal representation
!!         call data_functions_axpy(ogwfs_on_grid, old_ogwfs_on_grid, &
!!              xl_coef((icc-1)*xl_max+istep))
!!      end do
!!
!!      call data_set_to_zero(ngwfs_on_grid)
!!      call data_functions_copy(ngwfs_on_grid,ogwfs_on_grid)
!!
!!      ! vv: Shave ngwfs
!!      do ingwf = 1,ngwf_basis%proc_num
!!         call basis_clean_function(ngwfs_on_grid,ngwf_basis%spheres(ingwf),&
!!              cell,fftbox)
!!      end do
!!
!!      ! Store the propagated NGWFs_aux and Lowdin overlap
!!      call elec_history_store_ngwfs(ngwfs_on_grid, ngwf_basis, &
!!           cell,fftbox,elements, 'init')
!!      cstatus = 0
!!
!!      call data_functions_dealloc(old_ogwfs_on_grid)
!!      call data_functions_dealloc(ogwfs_on_grid)
!!      call data_functions_dealloc(new_ogwfs_on_grid(1))
!!
!!  end subroutine internal_dissip_ortho_propagation
!!
!!-----------------------------------------------------------------------------!
!!
!!  subroutine internal_unit_transf(Umat,scf_ogwfs,aux_ogwfs)
!!
!!    use function_ops, only : function_ops_brappd_ketppd
!!    use smearing_operator, only : smearing_matrix, lowdin_transformation
!!    use sparse, only : SPAM3, sparse_create, sparse_destroy, &
!!         sparse_transpose, sparse_transpose_structure, sparse_product, &
!!         sparse_trace, sparse_show_matrix
!!
!!    implicit none
!!    ! Arguments
!!    type(SPAM3), intent(inout) :: Umat
!!    type(FUNCTIONS), intent(in) :: scf_ogwfs
!!    type(FUNCTIONS), intent(in) :: aux_ogwfs
!!
!!    ! Local variables
!!    type(smearing_matrix) :: Omat, inv_sqrt_Omat
!!    type(SPAM3) :: overlap, trans_overlap, TUmat
!!    real(kind=DP) :: trace
!!    integer :: num_ngwfs
!!
!!    if (pub_debug_on_root) write(stdout,'(/a)') &
!!        'DEBUG: Entering elec_history_internal_unit_transf'
!!
!!    Omat%matrix_type = 3
!!    inv_sqrt_Omat%matrix_type = 3
!!    Omat%dataSPAM3%structure = 'K'
!!    inv_sqrt_Omat%dataSPAM3%structure = 'K'
!!    call sparse_create(Omat%dataSPAM3)
!!    call sparse_create(inv_sqrt_Omat%dataSPAM3)
!!
!!    overlap%structure = 'K'
!!    call sparse_create(overlap)
!!    call sparse_transpose_structure(trans_overlap%structure,overlap)
!!    call sparse_create(trans_overlap)
!!
!!    ! vv: Calculate the overlap between converged and auxiliary orthogonal
!!    ! functions, i.e. O = <Pshi|Phi>
!!    call function_ops_brappd_ketppd(overlap,scf_ogwfs,ngwf_basis,aux_ogwfs,&
!!         ngwf_basis,cell)
!!    ! vv: Calculate adjoint, i.e. transpose for real matrices
!!    call sparse_transpose(trans_overlap,overlap)
!!    ! vv: Calculate the product O*O^T
!!    call sparse_product(Omat%dataSPAM3,overlap,trans_overlap)
!!
!!    ! vv: Perform Lowdin transformation, i.e. (O*O^T)^(-1/2)
!!    call lowdin_transformation(Omat,inv_sqrt_Omat, tol=10e-12_dp)
!!
!!    ! vv: Generate Unitary transformation for phase alignement,
!!    !     i.e. U=[(O*O^T)^(-1/2)]*O
!!    call sparse_product(Umat,inv_sqrt_Omat%dataSPAM3,overlap)
!!
!!    ! vv: Clean workspace
!!    call sparse_destroy(Omat%dataSPAM3)
!!    call sparse_destroy(inv_sqrt_Omat%dataSPAM3)
!!    call sparse_destroy(overlap)
!!    call sparse_destroy(trans_overlap)
!!    call sparse_destroy(TUmat)
!!
!!    if (pub_debug_on_root) write(stdout,'(/a)') &
!!        'DEBUG: Leaving elec_history_internal_unit_transf'
!!
!!  end subroutine internal_unit_transf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  subroutine internal_dissip_propagation()
!!
!!     !=============================================================!
!!     ! Propagation of auxiliary support functions :
!!     ! algorithm based on A. Niklasson et al,
!!     !           J. Chem. Phys. 130, 214109 (2009)
!!     !           section IV, case 3 <= K <= 9
!!     ! Modified by Andrea Greco to allow use of complex NGWFs
!!     !=============================================================!
!!
!!     use basis, only: basis_clean_function
!!     use datatypes, only: data_set_to_zero
!!
!!     implicit none
!!
!!     type(FUNCTIONS) :: old_ngwfs_on_grid
!!     integer     :: icc
!!     integer  :: istep, ingwf
!!     integer       :: xl_max = 10
!!     real(kind=DP), parameter, dimension(7)   :: xl_k = (/  &
!!         1.690_dp, 1.750_dp, 1.820_dp, 1.840_dp, 1.860_dp, 1.880_dp, 1.890_dp /)
!!     real(kind=DP), parameter, dimension(70)  :: xl_coef = (/  &
!!     ! vv : Set of correct coefficients
!!           0.01000_dp, -0.55000_dp,  0.00000_dp, -0.15000_dp,  0.00000_dp,  0.00000_dp,&
!!           0.00000_dp,  0.00000_dp,  0.00000_dp,  0.00000_dp,&
!!           0.07900_dp, -0.65800_dp, -0.11400_dp, -0.11400_dp,  0.05700_dp,  0.00000_dp,&
!!           0.00000_dp,  0.00000_dp,  0.00000_dp,  0.00000_dp,&
!!           0.07200_dp, -0.74800_dp, -0.14400_dp, -0.05400_dp,  0.07200_dp, -0.01800_dp,&
!!           0.00000_dp,  0.00000_dp,  0.00000_dp,  0.00000_dp,&
!!           0.08300_dp, -0.80200_dp, -0.14850_dp, -0.01100_dp,  0.06600_dp, -0.03300_dp,&
!!           0.00550_dp,  0.00000_dp,  0.00000_dp,  0.00000_dp,&
!!           0.08240_dp, -0.84160_dp, -0.14080_dp,  0.01760_dp,  0.05120_dp, -0.04000_dp,&
!!           0.01280_dp, -0.00160_dp,  0.00000_dp,  0.00000_dp,&
!!           0.07644_dp, -0.87416_dp, -0.12584_dp,  0.03432_dp,  0.03432_dp, -0.03960_dp,&
!!           0.01848_dp, -0.00440_dp,  0.00044_dp,  0.00000_dp,&
!!           0.07568_dp, -0.89704_dp, -0.11232_dp,  0.04368_dp,  0.02016_dp, -0.03600_dp,&
!!           0.02208_dp, -0.00756_dp,  0.00144_dp, -0.00012_dp /)
!!
!!     if (pub_on_root) write(stdout,'(a)') &
!!        'Electronic history : dissipative propagation of auxiliary NGWFs'
!!
!!     !allocate(old_ngwfs_on_grid(ngwf_basis%size_on_grid), stat=ierr)
!!     !call utils_alloc_check('elec_history_compose_ngwfs','old_ngwfs_on_grid',ierr)
!!     ! agrecocmplx: allocate with appropriate routines
!!     call data_functions_alloc(old_ngwfs_on_grid, ngwf_basis%size_on_grid, &
!!                                iscmplx=ngwfs_on_grid%iscmplx)
!!
!!     ! Dissipative propagation of NGWFs_aux
!!     label = methods(iom)%ptr(1)
!!     call elec_history_read_ngwfs(old_ngwfs_on_grid,ngwf_basis, &
!!          cell,fftbox,elements,'scf',label)
!!
!!     icc = methods(iom)%ngwfs_mix_size-3
!!     !ngwfs_on_grid(:) = xl_k(icc) * old_ngwfs_on_grid(:)
!!     ! agrecocmplx
!!     call data_functions_axpy(ngwfs_on_grid, old_ngwfs_on_grid, &
!!          xl_k(icc))
!!
!!     do istep = 1, methods(iom)%ngwfs_mix_size
!!        label = methods(iom)%ptr(istep)
!!        call elec_history_read_ngwfs(old_ngwfs_on_grid,ngwf_basis, &
!!             cell,fftbox,elements,'init',label)
!!        ! agrecocmplx
!!        call data_functions_axpy(ngwfs_on_grid, old_ngwfs_on_grid, &
!!                                  xl_coef((icc-1)*xl_max+istep))
!!     enddo
!!
!!     ! Store the propagated NGWFs_aux
!!     call elec_history_store_ngwfs(ngwfs_on_grid, ngwf_basis, &
!!          cell,fftbox,elements, 'init')
!!     cstatus = 0
!!
!!     !deallocate(old_ngwfs_on_grid, stat=ierr)
!!     !call utils_dealloc_check('elec_history_compose_ngwfs','old_ngwfs_on_grid',ierr)
!!     ! agrecocmplx: deallocate with appropriate routine
!!     call data_functions_dealloc(old_ngwfs_on_grid)
!!
!!  end subroutine internal_dissip_propagation

    !################################################################!

  end subroutine elec_history_compose_ngwfs

!============================================================================!
!============================================================================!
!============================================================================!

  subroutine generalized_coordinates(ncoords,cell,old_coords,coord,iref)

    use constants, only: DP
    !use parallel_strategy, only: par=>pub_par
    use simulation_cell, only: CELL_INFO
    use services, only: services_rationalise_coords
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in)  :: ncoords
    integer, intent(in)  :: iref
    type(CELL_INFO), intent(in)    :: cell
    real(kind=DP), intent(inout) :: coord(:,:)! size (3,par%nat)
    real(kind=DP), intent(inout) :: old_coords(:,:,:)! size (3,par%nat,ncoords)


    ! Local variables
    real(kind=DP) :: ref_point(3)
    real(kind=DP) :: center_shift(3), dist
    real(kind=DP) :: atom_weight(size(coord,dim=2)) ! size (par%nat)
    integer       :: ic, iat, iom
    integer       :: nat ! rc2013: number of atoms

    ! rc2013: size of coord corresponds to no. of atoms in this structure
    nat = size(coord,dim=2)
    ! Current mixing method
    iom = methods_io

    ! Geometric center of simulation cell
    center_shift(1) = (cell%a1%x + cell%a2%x + cell%a3%x)/2.0d0
    center_shift(2) = (cell%a1%y + cell%a2%y + cell%a3%y)/2.0d0
    center_shift(3) = (cell%a1%z + cell%a2%z + cell%a3%z)/2.0d0

    ! Center coordinates around ref_point and rationalise
    do ic = 1, ncoords
       ref_point=old_coords(:,iref,ic)
       do iat = 1, nat
          old_coords(:,iat,ic) = old_coords(:,iat,ic) - ref_point(:) &
                                   + center_shift(:)
       enddo
       call services_rationalise_coords(nat,old_coords(:,:,ic), &
            cell)
    enddo
    ref_point=coord(:,iref)
    do iat = 1, nat
       coord(:,iat) = coord(:,iat) - ref_point(:) &
                                + center_shift(:)
    enddo
    call services_rationalise_coords(nat,coord,cell)

    ! Compute atomic weight
    atom_weight(:) = 0.0d0
    do iat = 1, nat
       dist = sqrt(sum((coord(:,iat)-ref_point(:))**2))
       if (methods(iom)%mix_smear.le.0.0d0) then
          if (dist.lt.methods(iom)%mix_radius) atom_weight(iat) = 1.0d0
       else
          atom_weight(iat) = (exp((dist-methods(iom)%mix_radius)/methods(iom)%mix_smear)+1)**(-1)
       endif
    enddo

    ! Apply atomic weight to coordinates
    do iat = 1, nat
       if (atom_weight(iat).gt.0.0d0) then
          coord(:,iat) = coord(:,iat)*atom_weight(iat)
       endif
    enddo
    do ic = 1, ncoords
       do iat = 1, nat
          if (atom_weight(iat).gt.0.0d0) then
             old_coords(:,iat,ic) = old_coords(:,iat,ic)*atom_weight(iat)
          endif
       enddo
    enddo

  end subroutine generalized_coordinates

!============================================================================!
!============================================================================!
!============================================================================!

  subroutine elec_history_backup_history(nat, mrows, ncols, iscmplx, &
       cell, fftbox)

    use datatypes, only: FUNCTIONS, data_functions_alloc, data_functions_dealloc
    use comms, only: pub_on_root, comms_bcast, comms_barrier
    use constants, only: DP, stdout
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use ion, only : ELEMENT
    use restart, only: restart_kernel_write, restart_ngwfs_tightbox_output
    use rundat, only: pub_num_spins, pub_num_kpoints, PUB_1K, &
         pub_rootname, md_write_out
    use simulation_cell, only: CELL_INFO
    use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_array_create, &
         sparse_embed_array_destroy
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort,&
         utils_unit, utils_open_unit_check, utils_close_unit_check, utils_assert


    implicit none

    ! Arguments
    integer, intent(in) :: nat
    integer, intent(in) :: mrows
    integer, intent(in) :: ncols
    logical, intent(in) :: iscmplx
    type(CELL_INFO), intent(in)    :: cell
    type(FFTBOX_INFO), intent(inout) :: fftbox

    ! Local variables
    type(SPAM3_EMBED_ARRAY) :: tmp_dkn, tmp_dkn_vel
    integer :: ierr, iom, ii, idir, iat, idx, iptr
    integer :: io_unit, io_stat
    integer, allocatable :: loc_ptr(:)
    character(len=80) :: filename
    character(len=2)  :: loc_struc, locvel_struc

    !type(FUNCTIONS) :: tmp_ngwfs_on_grid
    !character(len=22) :: f_extension

    if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering elec_history_backup_history'

    ! jme: KPOINTS_DANGER
    ! Some parts of this function enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine elec_history_backup_history not ready yet &
         &for more than one k-point.')

    ! Find an index for this method
    iom = methods_io

    if(methods(iom)%dkn_mix_type == 'XLD' .or. &
       methods(iom)%dkn_mix_type == 'XLI' .or. &
       methods(iom)%dkn_mix_type == 'XLIS' .or. &
       methods(iom)%dkn_mix_type == 'NAIVE') then
       loc_struc = 'KS'
       locvel_struc = 'K'
    else
       loc_struc = 'K'
       locvel_struc = 'K'
    end if

    ! agrecocmplx
    call sparse_embed_array_create(tmp_dkn, iscmplx=iscmplx, &
         n_spins=pub_num_spins, n_kpoints=pub_num_kpoints, &
         structure=trim(loc_struc))
    ! agrecocmplx
    call sparse_embed_array_create(tmp_dkn_vel, iscmplx=iscmplx, &
         n_spins=pub_num_spins, n_kpoints=pub_num_kpoints, &
         structure=trim(locvel_struc))

    ! vv: NGWFS_PROPAGATION
    ! agrecocmplx
    ! call data_functions_alloc(tmp_ngwfs_on_grid, ngwf_basis%size_on_grid, &
    !      iscmplx=iscmplx)

    ! Abort if stores haven't been initialized
    if (.not. allocated(history) .or. .not.allocated(methods) .or. &
        history_io==0 .or. methods_io==0 .or. history_count==0) then
       call utils_abort('Error in elec_hist_restore_history:&
            & store have not been initialized')
    endif

    if ( history_size < history_max ) then
       call utils_abort('Error in elec_hist_restore_history:&
            & electronic history store has not been filled entirely.')
    end if

    ! vv: Here perform an on-line call to elec_history_update_method
    if(methods(iom)%maxit .gt. 0) then
       allocate(loc_ptr(methods(iom)%maxit), stat=ierr)
       call utils_alloc_check('elec_history_backup_history','loc_ptr',ierr)
       loc_ptr(:) = methods(iom)%ptr(:)
       do iptr = methods(iom)%maxit, 2, -1
         loc_ptr(iptr) = loc_ptr(iptr-1)
       enddo
       loc_ptr(1) = history(history_io)%label
    end if

    ! Open the external file with extension .history.info
    if (pub_on_root) then
       io_unit = utils_unit()
       write(filename,'(a,a)') trim(pub_rootname), '.history.info'
       open(unit=io_unit,iostat=io_stat,file=filename,&
            access='SEQUENTIAL',form='UNFORMATTED',position='REWIND',&
            action='WRITE',status='UNKNOWN')
       call utils_open_unit_check('backup_history',filename,io_stat)

       ! vv: Save composition parameters for continuation runs
       write(io_unit) methods(iom)%dkn_count + 1
       write(io_unit) methods(iom)%ngwfs_count + 1
       write(io_unit) methods(iom)%dkn_mix_type
       write(io_unit) methods(iom)%dkn_mix_size
       write(io_unit) methods(iom)%dkn_mix_reset
       write(io_unit) methods(iom)%dkn_init_num
       write(io_unit) methods(iom)%dkn_init_type
       write(io_unit) methods(iom)%ngwfs_mix_type
       write(io_unit) methods(iom)%ngwfs_mix_size
       write(io_unit) methods(iom)%ngwfs_mix_reset
       write(io_unit) methods(iom)%ngwfs_init_num
       write(io_unit) methods(iom)%ngwfs_init_type
       write(io_unit) methods(iom)%ngwfs_coeff
       write(io_unit) methods(iom)%maxit
       write(io_unit) methods(iom)%mix_radius
       write(io_unit) methods(iom)%mix_smear
       write(io_unit) methods(iom)%name
       write(io_unit) methods(iom)%size
       write(io_unit) methods_io
       write(io_unit) methods_num
       write(io_unit) methods_max
       write(io_unit) history_io
       write(io_unit) history_size
       write(io_unit) history_count
       write(io_unit) history_max
       do ii =1,methods(iom)%maxit
          write(io_unit) loc_ptr(ii)
       end do
    end if

    if(pub_on_root) then
       close(unit=io_unit,iostat=io_stat)
       call utils_close_unit_check('backup_history',filename,io_stat)
    end if

    ! Open the external file with extension .history.var
    if (pub_on_root) then
       io_unit = utils_unit()
       write(filename,'(a,a)') trim(pub_rootname), '.history.var'
       open(unit=io_unit,iostat=io_stat,file=filename,&
            access='SEQUENTIAL',form='UNFORMATTED',position='REWIND',&
            action='WRITE',status='UNKNOWN')
       call utils_open_unit_check('backup_history',filename,io_stat)
    end if

    ! Start the main loop
    if (pub_on_root) then
       do ii = 1, history_max
          ! Start writing
          write(io_unit) history_index(history(ii)%label)
          write(io_unit) history(ii)%label

          do iat = 1,nat
             do idir = 1,3
                write(io_unit) history(ii)%coord(idir,iat)
             end do
          end do

         !if (methods(iom)%dkn_mix_type == 'XLD' .or. &
         !    methods(iom)%dkn_mix_type == 'XLI' .or. &
         !    methods(iom)%dkn_mix_type == 'NAIVE' ) then

         !   write(io_unit) size(history(ii)%dkn%idx,1)
         !   write(io_unit) size(history(ii)%dkn%idxinit,1)
         !end if
         !if (methods(iom)%dkn_mix_type == 'XLI' .or. &
         !    methods(iom)%dkn_mix_type == 'NAIVE' ) then

         !   write(io_unit) size(history(ii)%dkn%idxvel,1)
         !end if

         !if (methods(iom)%dkn_mix_type == 'XLD' .or. &
         !    methods(iom)%dkn_mix_type == 'XLI' .or. &
         !    methods(iom)%dkn_mix_type == 'NAIVE' ) then

         !   do idx = 1,size(history(ii)%dkn%idx,1)
         !      write(io_unit) history(ii)%dkn%idx(idx,PUB_1K)
         !   end do

         !   do idx = 1,size(history(ii)%dkn%idxinit,1)
         !      write(io_unit) history(ii)%dkn%idxinit(idx,PUB_1K)
         !   end do

         !end if

         !if (methods(iom)%dkn_mix_type == 'XLI' .or. &
         !    methods(iom)%dkn_mix_type == 'NAIVE' ) then

         !   do idx = 1, size(history(ii)%dkn%idxvel,1)
         !      write(io_unit) history(ii)%dkn%idxvel(idx,PUB_1K)
         !   end do

         !end if

       ! vv: NGWFS_PROPAGATION
       !  do idx = 1,ngwf_basis%num
       !     write(io_unit) history(ii)%ngwfs%procs(idx)
       !  end do
        end do
     end if

     if(pub_on_root) then
        close(unit=io_unit,iostat=io_stat)
        call utils_close_unit_check('backup_history',filename,io_stat)
     end if

    md_write_out = .false.
    do ii = 1, history_max
     ! vv: NGWFS_PROPAGATION
     ! f_extension = 'tightbox_ngwfs_history'
     ! f_extension = trim(f_extension)
     ! if (methods(iom)%ngwfs_mix_type == 'TRPROP' &
     !     .or. methods(iom)%ngwfs_mix_type == 'XLD') then
     !    call elec_history_read_ngwfs(tmp_ngwfs_on_grid,ngwf_basis, &
     !         cell,fftbox,elements,'init',history(ii)%label)
     !    call restart_ngwfs_tightbox_output(tmp_ngwfs_on_grid, ngwf_basis,&
     !         cell, fftbox, elements, f_extension, history(ii)%label,'init')
     ! end if
     !
     ! call elec_history_read_ngwfs(tmp_ngwfs_on_grid,ngwf_basis, &
     !      cell,fftbox,elements,'scf',history(ii)%label)
     !
     ! call restart_ngwfs_tightbox_output(tmp_ngwfs_on_grid, ngwf_basis,&
     !      cell, fftbox, elements, f_extension, history(ii)%label,'scf')

       if (methods(iom)%dkn_mix_type == 'TRPROP') then
          call elec_history_read_dkn(tmp_dkn,'init',history(ii)%label)
          call restart_kernel_write(tmp_dkn,write_cond=.false.,&
               label=history_index(history(ii)%label),dkn_type='init')
       else if (methods(iom)%dkn_mix_type == 'XLD') then
          call elec_history_read_dkn(tmp_dkn,'init',history(ii)%label)
          call restart_kernel_write(tmp_dkn,write_cond=.false.,&
               label=history_index(history(ii)%label),dkn_type='init')
          call elec_history_read_dkn(tmp_dkn,'scf',history(ii)%label)
          call restart_kernel_write(tmp_dkn,write_cond=.false.,&
               label=history_index(history(ii)%label),dkn_type='scf')
       else if (methods(iom)%dkn_mix_type == 'XLI'  .or. &
                methods(iom)%dkn_mix_type == 'XLIS' .or. &
                methods(iom)%dkn_mix_type == 'NAIVE' ) then
          call elec_history_read_dkn(tmp_dkn,'init',history(ii)%label)
          call restart_kernel_write(tmp_dkn,write_cond=.false.,&
               label=history_index(history(ii)%label),dkn_type='init')
          call elec_history_read_dkn(tmp_dkn,'scf',history(ii)%label)
          call restart_kernel_write(tmp_dkn,write_cond=.false.,&
               label=history_index(history(ii)%label),dkn_type='scf')
          call elec_history_read_dkn(tmp_dkn_vel,'vel',history(ii)%label)
          call restart_kernel_write(tmp_dkn_vel,write_cond=.false.,&
               label=history_index(history(ii)%label),dkn_type='vel')
       else
          call elec_history_read_dkn(tmp_dkn,'scf',history(ii)%label)
          call restart_kernel_write(tmp_dkn,write_cond=.false.,&
               label=history_index(history(ii)%label),dkn_type='scf')
       end if
    end do
    md_write_out = .true.

    if(allocated(loc_ptr)) then
      deallocate(loc_ptr,stat=ierr)
      call utils_dealloc_check('elec_history_backup_history','loc_ptr',ierr)
    end if

    ! Deallocate workspace, i.e. sick and destroy
    ! vv: NGWFS_PROPAGATION
    ! call data_functions_dealloc(tmp_ngwfs_on_grid)
    call sparse_embed_array_destroy(tmp_dkn)
    call sparse_embed_array_destroy(tmp_dkn_vel)

  end subroutine elec_history_backup_history

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine elec_history_restore_history_info(nat,nsub)

    ! vv: This subroutine takes in input all the infos
    use constants, only: DP, stdout
    use comms, only: pub_on_root
    use utils, only: utils_alloc_check, utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: nat, nsub

    ! Local variables
    integer :: dkn_count, ngwfs_count
    integer :: dkn_mix_num
    integer :: ngwfs_mix_num
    integer :: dkn_mix_init_num
    integer :: ngwfs_mix_init_num
    integer :: dkn_mix_reset
    integer :: ngwfs_mix_reset
    integer :: maxit
    real(kind=DP) :: mix_radius
    real(kind=DP) :: mix_smear
    real(kind=DP) :: mix_coeff
    character(len=10) :: dkn_mix_type
    character(len=10) :: dkn_mix_init_type
    character(len=10) :: ngwfs_mix_type
    character(len=10) :: ngwfs_mix_init_type

    ! Local variables
    integer :: ierr, idx

    if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering elec_history_restore_history_info'

    ! Check that history has not been allocated yet
    call utils_assert(.not. allocated(history),'Error in elec_history_create_storage : &
          &attempt to allocate history twice!')

    ! vv: read the info from the external file .history.info
    if (pub_on_root) write(stdout,'(a)',advance='no') &
        '- MD -: Re-filling main electronic history storage...'
    call internal_read_info
    if (pub_on_root) write(stdout,*) 'done'

    ! Report info on stdout
    if (methods_num .ge. 1 .and. pub_on_root) then
       do idx = 1, methods_num
          write(stdout,*)
          write(stdout,'(a,i2.2,2a)') 'Electronic history : ', &
               methods(idx)%maxit,' steps required for ', &
               methods(idx)%name
       enddo
       write(stdout,'(a,i2.2)') &
          'Electronic history: Total number of steps to be kept in memory = ', history_max
    endif

    ! Allocate history
    if (history_max .gt. 1) then
       allocate(history(history_max),stat=ierr)
       call utils_alloc_check('elec_history_restore_history_info','history',ierr)
    endif

    ! vv: allocate the coordinates in the history structure
    ! jcap: also allocate the ngwfs
    do idx = 1, history_max
       allocate(history(idx)%coord(3,nat),stat=ierr)
       call utils_alloc_check('elec_history_restore_history_info','history(idx)%coord',ierr)
       allocate(history(idx)%ngwfs(nsub),stat=ierr)
       call utils_alloc_check('elec_history_restore_history_info','history(idx)%ngwfs',ierr)
    end do

    if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving elec_history_restore_history_info'

contains

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

    subroutine internal_read_info

      use comms, only: pub_on_root, comms_bcast, pub_root_proc_id
      use rundat, only: pub_rootname
      use utils, only: utils_open_unit_check, utils_close_unit_check,&
           utils_abort, utils_unit, utils_alloc_check, utils_dealloc_check

      implicit none

      ! Local variables
      character(len=80) :: filename
      character(len=10) :: methods_name
      integer :: methods_size
      logical :: fileexists
      integer :: io_unit, io_stat, ii, idx
      integer, allocatable :: loc_ptr(:)

      ! Open the external file with extension .history
      if (pub_on_root) then
         write(filename,'(a,a)') trim(pub_rootname), '.history.info'

         inquire(file=filename,exist=fileexists)

         if (fileexists) then
            io_unit = utils_unit()

            open(unit=io_unit,iostat=io_stat,file=filename,&
                 access='SEQUENTIAL',form='UNFORMATTED',position='REWIND',&
                 action='READ')
            call utils_open_unit_check('internal_read_info',filename,io_stat)
         else
            call utils_abort('No history file found!')
         end if
      end if

      if(pub_on_root) then
         ! vv: Restore composition parameters for continuation runs
         read(io_unit)  dkn_count
         read(io_unit)  ngwfs_count
         read(io_unit)  dkn_mix_type
         read(io_unit)  dkn_mix_num
         read(io_unit)  dkn_mix_reset
         read(io_unit)  dkn_mix_init_num
         read(io_unit)  dkn_mix_init_type
         read(io_unit)  ngwfs_mix_type
         read(io_unit)  ngwfs_mix_num
         read(io_unit)  ngwfs_mix_reset
         read(io_unit)  ngwfs_mix_init_num
         read(io_unit)  ngwfs_mix_init_type
         read(io_unit)  mix_coeff
         read(io_unit)  maxit
         read(io_unit)  mix_radius
         read(io_unit)  mix_smear
         read(io_unit)  methods_name
         read(io_unit)  methods_size
         read(io_unit)  methods_io
         read(io_unit)  methods_num
         read(io_unit)  methods_max
         read(io_unit)  history_io
         read(io_unit)  history_size
         read(io_unit)  history_count
         read(io_unit)  history_max
      end if

      ! vv: Let all cores know the values
      call comms_bcast(pub_root_proc_id, dkn_count          )
      call comms_bcast(pub_root_proc_id, ngwfs_count        )
      call comms_bcast(pub_root_proc_id, dkn_mix_type       )
      call comms_bcast(pub_root_proc_id, dkn_mix_num        )
      call comms_bcast(pub_root_proc_id, dkn_mix_reset      )
      call comms_bcast(pub_root_proc_id, dkn_mix_init_num   )
      call comms_bcast(pub_root_proc_id, dkn_mix_init_type  )
      call comms_bcast(pub_root_proc_id, ngwfs_mix_type     )
      call comms_bcast(pub_root_proc_id, ngwfs_mix_num      )
      call comms_bcast(pub_root_proc_id, ngwfs_mix_reset    )
      call comms_bcast(pub_root_proc_id, ngwfs_mix_init_num )
      call comms_bcast(pub_root_proc_id, ngwfs_mix_init_type)
      call comms_bcast(pub_root_proc_id, mix_coeff          )
      call comms_bcast(pub_root_proc_id, maxit              )
      call comms_bcast(pub_root_proc_id, mix_radius         )
      call comms_bcast(pub_root_proc_id, mix_smear          )
      call comms_bcast(pub_root_proc_id, methods_name       )
      call comms_bcast(pub_root_proc_id, methods_size       )
      call comms_bcast(pub_root_proc_id, methods_io         )
      call comms_bcast(pub_root_proc_id, methods_num        )
      call comms_bcast(pub_root_proc_id, methods_max        )
      call comms_bcast(pub_root_proc_id, history_io         )
      call comms_bcast(pub_root_proc_id, history_size       )
      call comms_bcast(pub_root_proc_id, history_count      )
      call comms_bcast(pub_root_proc_id, history_max        )

      allocate(loc_ptr(maxit),stat = ierr)
      call utils_alloc_check('elec_history_restore_history_info','loc_ptr', ierr)

      if(pub_on_root) then
        do ii = 1,maxit
           read(io_unit) loc_ptr(ii)
        end do
      end if
      call comms_bcast(pub_root_proc_id, loc_ptr, maxit         )

      ! vv: Close the unit
      if (pub_on_root) then
         close(io_unit,iostat=io_stat)
         call utils_close_unit_check('backup_history',filename,io_stat)
      end if

      ! vv: Allocate e re-fill the methods structure
      if (.not. allocated(methods)) then
         allocate(methods(methods_max), stat=ierr)
         call utils_alloc_check('elec_history_restore_history_info','methods', &
              ierr)
      end if

      ! vv: Initialise methods
      idx = methods_io
      methods(idx)%name = methods_name
      methods(idx)%dkn_count = dkn_count
      methods(idx)%ngwfs_count = ngwfs_count
      methods(idx)%size = min(methods_size+1,maxit)
      methods(idx)%maxit = maxit
      if (methods(idx)%maxit > 0) then
         allocate(methods(idx)%ptr(maxit),stat=ierr)
         call utils_alloc_check('elec_history_initialise_method','methods(idx)%ptr',ierr)
         methods(idx)%ptr(:) = loc_ptr(:)
      endif

      ! Fill in mixing parameters
      methods(idx)%dkn_mix_type  = dkn_mix_type
      methods(idx)%dkn_mix_size  = dkn_mix_num
      methods(idx)%dkn_mix_reset = dkn_mix_reset
      methods(idx)%dkn_init_num  = dkn_mix_init_num
      methods(idx)%dkn_init_type = dkn_mix_init_type

      methods(idx)%ngwfs_mix_type  = ngwfs_mix_type
      methods(idx)%ngwfs_mix_size  = ngwfs_mix_num
      methods(idx)%ngwfs_mix_reset = ngwfs_mix_reset
      methods(idx)%ngwfs_init_num  = ngwfs_mix_init_num
      methods(idx)%ngwfs_init_type = ngwfs_mix_init_type
      methods(idx)%ngwfs_coeff = mix_coeff

      methods(idx)%mix_radius = mix_radius
      methods(idx)%mix_smear  = mix_smear

      ! vv: Deallocate workspace
      deallocate(loc_ptr, stat=ierr)
      call utils_dealloc_check('elec_history_restore_history_info','loc_ptr', ierr)

    end subroutine internal_read_info

  end subroutine elec_history_restore_history_info

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!


  subroutine elec_history_restore_history_parm(elements, ngwf_basis, &
       iscmplx, cell, fftbox)

    use datatypes, only: FUNCTIONS, data_functions_alloc, data_functions_dealloc
    use comms, only: comms_bcast, pub_root_proc_id, pub_on_root, comms_barrier
    use constants, only: stdout
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use ion, only : ELEMENT
    use md_ions, only: r_t
    use restart, only: restart_kernel_read, restart_ngwfs_tightbox_input
    use rundat, only: pub_num_spins, pub_num_kpoints, PUB_1K, &
         pub_rootname, md_write_out
    use simulation_cell, only: CELL_INFO
    use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_array_create, &
         sparse_embed_array_destroy
    use utils, only: utils_open_unit_check, utils_close_unit_check,&
         utils_assert, utils_alloc_check, utils_dealloc_check, &
         utils_unit, utils_abort, utils_banner

    implicit none

    ! Arguments
    type(ELEMENT), intent(in)        :: elements(:) ! size=(par%nat)
    type(FUNC_BASIS), intent(in)     :: ngwf_basis(:)
    logical, intent(in)              :: iscmplx
    type(CELL_INFO), intent(in)      :: cell
    type(FFTBOX_INFO), intent(inout) :: fftbox

    ! Local variables
    type(SPAM3_EMBED_ARRAY) :: tmp_dkn, tmp_dkn_vel
    type(FUNCTIONS), allocatable :: tmp_ngwfs_on_grid(:)
    real(kind=DP), allocatable :: coord(:,:)
    character(len=80) :: filename
    character(len=22) :: f_extension
    logical :: fileexists
    integer :: io_unit, io_stat, ierr, hist_label, hist_idx, max_label, ihist
    integer :: ii, idx, iat, idir, iom, iptr, idxsize, idxinitsize, idxvelsize
    integer, allocatable :: hist_index(:)
    integer :: nat, nsub

    if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering elec_history_restore_history_parm'

    ! rc2013: extract no. of atoms from size of elements array
    nat = size(elements)
    ! jcap: extract no. of subregions from size of ngwf_basis array
    nsub = size(ngwf_basis)
    ! jcap: allocate tmp_ngwfs_on_grid
!    allocate(tmp_ngwfs_on_grid(nsub),stat=ierr)
!    call utils_alloc_check('elec_history_restore_history_parm',&
!         'tmp_ngwfs_on_grid',ierr)

    ! jme: KPOINTS_DANGER
    ! Some parts of this function enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine elec_history_restore_history_parm not ready yet &
         &for more than one k-point.')

    ! vv: Allocate workspace
    ! Set the index and create the sparse matrices
    iom = methods_io

    if (methods(iom)%dkn_mix_type == 'TRPROP' .or. &
         methods(iom)%dkn_mix_type == 'XLD'   .or. &
         methods(iom)%dkn_mix_type == 'XLI'   .or. &
         methods(iom)%dkn_mix_type == 'XLIS'  .or. &
         methods(iom)%dkn_mix_type == 'NAIVE' ) then
       call sparse_embed_array_create(tmp_dkn, n_spins=pub_num_spins, &
            n_kpoints=pub_num_kpoints, structure='KS')
    else
       call sparse_embed_array_create(tmp_dkn, n_spins=pub_num_spins, &
            n_kpoints=pub_num_kpoints, structure='K')
    end if

    call sparse_embed_array_create(tmp_dkn_vel, n_spins=pub_num_spins, &
         n_kpoints=pub_num_kpoints, structure='K')


    allocate(coord(3,nat), stat=ierr)
    call utils_alloc_check('elec_history_restore_history_parm','coord',ierr)

    allocate(hist_index(history_max), stat=ierr)
    call utils_alloc_check('elec_history_restore_history_parm','hist_index',ierr)

    ! Open the external file with extension <filename>.history.var
    if (pub_on_root) then
       write(filename,'(a,a)') trim(pub_rootname), '.history.var'
       inquire(file=filename,exist=fileexists)
       if (fileexists) then
       ! vv: Find a free unit
       io_unit = utils_unit()
          open(unit=io_unit,iostat=io_stat,file=filename,&
               access='SEQUENTIAL',form='UNFORMATTED',position='REWIND',&
               action='READ')
          call utils_open_unit_check('elec_history_restore_history_parm',&
               filename,io_stat)
       else
          call utils_abort('File .history.var not found')
       end if
    end if

    ! vv: Banner
    if (pub_on_root) then
       write(stdout,'(/a/)') utils_banner('=', &
            'Re-filling the electronic structure history')
    end if

    ! vv: Main loop
    do ii = 1, history_max
       ! vv: Read in the label and the coordinates
       if (pub_on_root) then
          read(io_unit) hist_index(ii)
          read(io_unit) hist_label
          do iat = 1,nat
             do idir = 1,3
                read(io_unit) coord(idir,iat)
             end do
          end do
       end if
       ! vv: Broadcast the label and the coordinates
       call comms_bcast(pub_root_proc_id,hist_index(ii))
       call comms_bcast(pub_root_proc_id,hist_label)
       call comms_bcast(pub_root_proc_id,coord)
       ! vv: Create an entry in the elec_history container and print something
       !     to the stdout
       history(ii)%label = hist_label
       history(ii)%coord(:,:) = coord(:,:)

       if (pub_on_root .and. pub_md_output_detail>= NORMAL) &
          write(stdout,'(a,i2,a,i4,a)') &
          'Electronic history : new entry (idx = ',hist_index(ii),&
          ') (label = ',history(ii)%label,')'

    end do

    ! vv: Close the external file
    if(pub_on_root) then
       close(unit=io_unit,iostat=io_stat)
       call utils_close_unit_check('backup_history',filename,io_stat)
    end if


    ! vv: Now read the ngwfs and density kernels from external files
    !     and repopulate the electronic history
    !     N.B. we can avoid to fill the last entry in the history as
    !     this is not used in the propagation schemes, and it is
    !     removed at the start of a new calculation anyway.
    do ii = 1, methods(iom)%maxit
       hist_idx = history_index(methods(iom)%ptr(ii))
       hist_label = methods(iom)%ptr(ii)

       ! vv: NGWFS_PROPAGATION
       !f_extension = 'tightbox_ngwfs_history'
       !! vv: Read in the tightboxes with label=hist_label
       !!     N.B. We need the auxilary electronic degrees of freedom when we
       !!     use Niklasson's schemes
       !if (methods(iom)%ngwfs_mix_type == 'TRPROP' &
       !    .or. methods(iom)%ngwfs_mix_type == 'XLD') then
       !   call restart_ngwfs_tightbox_input(tmp_ngwfs_on_grid, ngwf_basis,&
       !        cell, fftbox, elements, f_extension, '.', hist_label, 'init')
       !   call elec_history_store_ngwfs(tmp_ngwfs_on_grid,ngwf_basis, &
       !     cell,fftbox,elements,'init', hist_label)
       !end if
       !call restart_ngwfs_tightbox_input(tmp_ngwfs_on_grid, ngwf_basis,&
       !     cell, fftbox, elements, f_extension, '.', hist_label, 'scf')
       !call elec_history_store_ngwfs(tmp_ngwfs_on_grid,ngwf_basis, &
       !     cell,fftbox,elements,'scf', hist_label)
       ! vv: Read in the density kernels with label=hist_label
       !     N.B. We need to read K.S instead of K when we use Niklasson's
       !     schemes

       if (methods(iom)%dkn_mix_type == 'TRPROP') then
          call restart_kernel_read(tmp_dkn, .false., .false.,&
                .false.,.false.,.false.,hist_idx, 'init')
          call elec_history_store_dkn(tmp_dkn,'init', hist_label)
       else if (methods(iom)%dkn_mix_type == 'XLD') then
          call restart_kernel_read(tmp_dkn, .false., .false.,&
                .false.,.false.,.false.,hist_idx, 'init')
          call elec_history_store_dkn(tmp_dkn,'init', hist_label)
          call restart_kernel_read(tmp_dkn, .false., .false., &
               .false.,.false.,.false.,hist_idx, 'scf')
          call elec_history_store_dkn(tmp_dkn,'scf', hist_label)
       else if (methods(iom)%dkn_mix_type == 'XLI' .or. &
                methods(iom)%dkn_mix_type == 'XLIS' .or.&
                methods(iom)%dkn_mix_type == 'NAIVE' ) then
          call restart_kernel_read(tmp_dkn, .false., .false.,&
                .false.,.false.,.false.,hist_idx, 'init')
          call elec_history_store_dkn(tmp_dkn,'init', hist_label)
          call restart_kernel_read(tmp_dkn, .false., .false., &
               .false.,.false.,.false.,hist_idx, 'scf')
          call elec_history_store_dkn(tmp_dkn,'scf', hist_label)
          call restart_kernel_read(tmp_dkn_vel, .false., .false., &
               .false.,.false.,.false.,hist_idx, 'vel')
          call elec_history_store_dkn(tmp_dkn_vel,'vel', hist_label)
       else
          call restart_kernel_read(tmp_dkn, .false., .false., &
               .false.,.false.,.false.,hist_idx, 'scf')
          call elec_history_store_dkn(tmp_dkn,'scf', hist_label)
       end if
    end do
    if(pub_on_root) write(stdout,'(/a/)') utils_banner('=')

    ! vv: Create a new entry
    call elec_history_create_entry(r_t)

    ! vv: Restore the option for printing out
    md_write_out = .true.

    ! vv: Deallocate workspace, i.e. sick and destroy
    ! vv: NGWFS_PROPAGATION
    ! call data_functions_dealloc(tmp_ngwfs_on_grid)
    call sparse_embed_array_destroy(tmp_dkn)
    call sparse_embed_array_destroy(tmp_dkn_vel)

    deallocate(coord, stat=ierr)
    call utils_dealloc_check('elec_history_restore_history_parm','coord',ierr)

    deallocate(hist_index, stat=ierr)
    call utils_dealloc_check('elec_history_restore_history_parm','hist_index',ierr)

    if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving elec_history_restore_history_parm'

  end subroutine elec_history_restore_history_parm

!#############################################################################!
!#############################################################################!
!#############################################################################!

  subroutine elec_history_update_dkn()

    use comms,  only: pub_on_root
    use constants, only: DP, stdout
    use rundat, only: md_aux_rep

    implicit none

    ! Local variables
    integer :: iom

    if (pub_on_root) write(stdout,'(a)') &
       'Electronic history : Entering velocity-Verlet update of &
       &auxiliary kernel velocity with Berendsen thermostat'

    ! vv: Set the index for this method
    iom = methods_io

    ! vv: If we are in the initial regime do nothing
    if (methods(iom)%dkn_count < methods(iom)%dkn_init_num ) then

       return

    else
       if (methods(iom)%dkn_mix_type == 'XLI' .or. methods(iom)%dkn_mix_type == 'NAIVE') then
          if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
             call internal_xl_berendsen_ortho_update
          else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
             call internal_xl_berendsen_asym_update
          end if
       else if (methods(iom)%dkn_mix_type == 'XLIS' ) then
          if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
             call internal_xl_berendsen_scuseria_ortho_update
          else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
             call internal_xl_berendsen_scuseria_asym_update
          end if
       end if

    end if

    if (pub_on_root) write(stdout,'(a)') &
       'Electronic history : Leaving velocity-Verlet update of &
       &auxiliary kernel velocity with Berendsen thermostat'

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine internal_xl_berendsen_ortho_update

    use rundat, only: md_delta_t,pub_num_spins, md_aux_dkn_t, md_aux_beren_tc,&
         pub_rootname, pub_num_kpoints, PUB_1K, md_autocorr
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product, sparse_embed_trace, &
         sparse_embed_array_create, sparse_embed_array_product, &
         sparse_embed_array_destroy, sparse_embed_array_axpy, &
         sparse_embed_array_scale, sparse_embed_array_num_rows
    use utils, only: utils_assert

    implicit none

    type(SPAM3_EMBED_ARRAY) :: old_x1, vel_x, x
    real(kind=DP) :: vel, alpha, tau
    real(kind=DP), dimension(pub_num_spins, PUB_1K) :: kappa
    integer :: label
    character(len=2) :: structure_type
    ! agrecocmplx
    logical :: loc_cmplx

    ! jmecmplx
    loc_cmplx = any(history(history_io)%dkn%scf(:)%iscmplx)

    ! vv: Define sparsity pattern for the orthogonal case
    structure_type = 'K'
    call sparse_embed_array_create(old_x1,n_spins=pub_num_spins, &
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(vel_x,n_spins=pub_num_spins, &
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(x,n_spins=pub_num_spins, &
         n_kpoints=pub_num_kpoints, structure=structure_type)

    if (md_autocorr) call elec_history_print_kernel(x,'real')

    ! vv: Set initial parameters
    alpha = 1.0_dp
    vel = 0.0_dp
    kappa = 1.00_DP
    tau = md_aux_beren_tc

    ! vv: Update velocities.
    ! vv: Calculate accelerations a_n+1
    ! kappa(S^1/2.K.S^1/2 - X)_i+1
    label = history(history_io)%label
    call elec_history_read_dkn(vel_x,'vel',label)
    call elec_history_read_dkn(x,'scf',label)
    call elec_history_read_dkn(old_x1,'init',label)
    call sparse_embed_array_axpy(x,old_x1,-1.0_DP)
    ! vv: d/dt (X)_i+1 = d/dt(X)_i + kappa(S^1/2.K.S^1/2 - X)_i+1
    call sparse_embed_array_axpy(vel_x,x,kappa)

    ! vv: Only Gamma-point for the moment
    call sparse_embed_trace(vel,vel_x%m(1,1),vel_x%m(1,1))
    vel = vel/(sparse_embed_array_num_rows(vel_x)*md_delta_t**2)
    call utils_assert(vel >= 0,'The temperature of the auxiliary degrees of freedom has &
         &become negative. This means that d(KS)/dt has complex eigenvalues. &
         &Try to increase the value of md_aux_beren_tc or md_aux_dkn_t', vel)

    ! vv: Compute the scaling parameter alpha
    if (methods(iom)%dkn_mix_type == 'XLI' .and. &
        (methods(iom)%dkn_count .ne. methods(iom)%dkn_init_num + 1)) then
       ! vv: Compute alpha
       alpha = sqrt(1.0_dp + md_delta_t/tau * (md_aux_dkn_t/vel - 1.0_dp))
       call utils_assert(alpha > 0,'The coupling constant alpha for the &
            &berendsen scaling of the auxiliary velocities has &
            &become negative. Try to increase the value of md_aux_beren_tc or md_aux_dkn_t',alpha)
    end if

    ! vv: Scale velocities according to dkn_mix_type
    call sparse_embed_array_scale(vel_x,alpha)

    ! vv: Store new velocities for next md step
    call elec_history_store_dkn(vel_x,'vel')
    ! print info to external files
    call elec_history_print_log(vel,'aux',alpha)

    ! vv: Free workspace
    call sparse_embed_array_destroy(old_x1)
    call sparse_embed_array_destroy(vel_x)
    call sparse_embed_array_destroy(x)

  end subroutine internal_xl_berendsen_ortho_update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine internal_xl_berendsen_asym_update

    use rundat, only: md_delta_t,pub_num_spins, md_aux_dkn_t, md_aux_beren_tc,&
         pub_rootname, pub_num_kpoints, PUB_1K, md_autocorr
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product, sparse_embed_trace, &
         sparse_embed_array_create, sparse_embed_array_product, &
         sparse_embed_array_destroy, sparse_embed_array_axpy, &
         sparse_embed_array_scale, sparse_embed_array_num_rows, &
         sparse_embed_array_transpose
    use utils, only: utils_assert

    implicit none

    type(SPAM3_EMBED_ARRAY) :: old_x1, vel_x, vel_x_trans, x
    real(kind=DP) :: vel, alpha, tau
    real(kind=DP), dimension(pub_num_spins, PUB_1K) :: kappa
    integer :: label
    character(len=2) :: structure_type
    ! agrecocmplx
    logical :: loc_cmplx

    ! jmecmplx
    loc_cmplx = any(history(history_io)%dkn%scf(:)%iscmplx)

    ! vv: Define sparsity pattern for the orthogonal case
    structure_type = 'KS'
    call sparse_embed_array_create(old_x1,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(vel_x,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(vel_x_trans,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(x,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)

    if (md_autocorr) call elec_history_print_kernel(x,'real')

    ! vv: Set initial parameters
    alpha = 1.0_dp
    vel = 0.0_dp
    kappa = 1.00_DP
    tau = md_aux_beren_tc

    ! vv: Update velocities.
    ! vv: Calculate accelerations a_n+1
    ! kappa(S^1/2.K.S^1/2 - X)_i+1
    label = history(history_io)%label
    call elec_history_read_dkn(vel_x,'vel',label)
    call elec_history_read_dkn(x,'scf',label)
    call elec_history_read_dkn(old_x1,'init',label)
    call sparse_embed_array_axpy(x,old_x1,-1.0_DP)
    ! vv: d/dt (X)_i+1 = d/dt(X)_i + kappa(S^1/2.K.S^1/2 - X)_i+1
    call sparse_embed_array_axpy(vel_x,x,kappa)

    ! vv: Only Gamma-point for the moment
    call sparse_embed_array_transpose(vel_x_trans,vel_x)
    call sparse_embed_trace(vel,vel_x%m(1,1),vel_x_trans%m(1,1))
    vel = vel/(sparse_embed_array_num_rows(vel_x)*md_delta_t**2)
    call utils_assert(vel >= 0,'The temperature of the auxiliary degrees of freedom has &
         &become negative. This means that d(KS)/dt has complex eigenvalues. &
         &Try to increase the value of md_aux_beren_tc or md_aux_dkn_t', vel)

    ! vv: Compute the scaling parameter alpha
    if (methods(iom)%dkn_mix_type == 'XLI' .and. &
        (methods(iom)%dkn_count .ne. methods(iom)%dkn_init_num + 1)) then
       ! vv: Compute alpha
       alpha = sqrt(1.0_dp + md_delta_t/tau * (md_aux_dkn_t/vel - 1.0_dp))
       call utils_assert(alpha > 0,'The coupling constant alpha for the &
            &berendsen scaling of the auxiliary velocities has &
            &become negative. Try to increase the value of md_aux_beren_tc or md_aux_dkn_t',alpha)
    end if

    ! vv: Scale velocities according to dkn_mix_type
    call sparse_embed_array_scale(vel_x,alpha)

    ! vv: Store new velocities for next md step
    call elec_history_store_dkn(vel_x,'vel')
    ! print info to external files
    call elec_history_print_log(vel,'aux',alpha)

    ! vv: Free workspace
    call sparse_embed_array_destroy(old_x1)
    call sparse_embed_array_destroy(vel_x)
    call sparse_embed_array_destroy(vel_x_trans)
    call sparse_embed_array_destroy(x)

  end subroutine internal_xl_berendsen_asym_update


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine internal_xl_berendsen_scuseria_ortho_update

    use rundat, only: md_delta_t,pub_num_spins, md_aux_dkn_t, md_aux_beren_tc,&
         pub_rootname, pub_num_kpoints, PUB_1K, md_autocorr
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
         sparse_embed_destroy, sparse_embed_product, sparse_embed_trace, &
         sparse_embed_array_create, sparse_embed_array_product, &
         sparse_embed_array_destroy, sparse_embed_array_axpy, &
         sparse_embed_array_scale, sparse_embed_array_num_rows, &
         sparse_embed_array_copy
    use utils, only: utils_assert

    implicit none

    type(SPAM3_EMBED_ARRAY) :: old_x1, vel_x, x, tmp_vel, xwx, qwq, q
    real(kind=DP) :: vel, alpha, tau, trace
    real(kind=DP), dimension(pub_num_spins, PUB_1K) :: kappa
    character(len=2) :: structure_type
    integer :: label
    ! agrecocmplx
    logical :: loc_cmplx

    ! jmecmplx
    loc_cmplx = any(history(history_io)%dkn%scf(:)%iscmplx)

    ! vv: Define sparsity pattern for the orthogonal case
    structure_type = 'K'
    call sparse_embed_array_create(old_x1,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(vel_x,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(x,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(tmp_vel,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(q,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(xwx,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(qwq,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)

    ! vv: Print out density kernel matrices if required
    if (md_autocorr) call elec_history_print_kernel(x,'real')

    ! vv: Set initial parameters
    alpha = 0.0_dp
    vel = 0.0_dp
    kappa = 1.0_DP
    tau = md_aux_beren_tc

    ! vv: (d X / dt)_n+1 = (d X / dt)_n + (dt / 2) * ( a_n+1 + a_n)
    ! vv: Call velocities from previous step

    ! vv: Calculate accelerations a_n+1
    ! w * KS_n+1
    label = history(history_io)%label
    call elec_history_read_dkn(vel_x,'vel',label)
    call elec_history_read_dkn(x,'scf',label)
    call elec_history_read_dkn(old_x1,'init',label)
    ! -w * X_n+1
    call sparse_embed_array_axpy(x,old_x1,-1.0_DP)
    call sparse_embed_array_axpy(vel_x,x,kappa)

    call sparse_embed_array_product(tmp_vel,vel_x,old_x1)
    call sparse_embed_array_product(xwx,old_x1,tmp_vel)

    call sparse_embed_array_copy(q,old_x1)
    call sparse_embed_array_scale(q,-1.0_DP,1.0_DP) ! Q = I - KS
    call sparse_embed_array_product(tmp_vel,vel_x,q)
    call sparse_embed_array_product(qwq,q,tmp_vel)

    call sparse_embed_array_axpy(xwx,qwq,1.0_DP)
    call sparse_embed_array_axpy(vel_x,xwx,-1.0_DP)

    ! vv: Only Gamma-point for the moment
    call sparse_embed_trace(trace,vel_x%m(1,1),vel_x%m(1,1))
    vel = vel + trace
    vel = vel/(sparse_embed_array_num_rows(vel_x)*md_delta_t**2)
    call utils_assert(vel >= 0,'The temperature of the auxiliary degrees of freedom has &
         &become negative. This means that d(KS)/dt has complex eigenvalues. &
         &Try to increase the value of md_aux_beren_tc or md_aux_dkn_t', vel)

    ! vv: Scale velocities according to the mix type
    if (methods(iom)%dkn_mix_type == 'XLIS' .and. &
         (methods(iom)%dkn_count .ne. methods(iom)%dkn_init_num + 1)) then
       ! vv: Compute alpha
       alpha = alpha + sqrt(1.0_dp + md_delta_t/tau * (md_aux_dkn_t/vel - 1.0_dp))
       call utils_assert(alpha > 0,'The coupling constant alpha for the &
            &berendsen scaling of the auxiliary velocities has &
            &become negative. Try to increase the value of md_aux_beren_tc or md_aux_dkn_t',alpha)
    else
       alpha = 1.0_dp
    end if

    call sparse_embed_array_scale(vel_x,alpha)
    ! vv: Store new velocities for next md step
    call elec_history_store_dkn(vel_x,'vel')
    ! print info to external files
    call elec_history_print_log(vel,'aux',alpha)


    ! vv: Free workspace
    call sparse_embed_array_destroy(old_x1)
    call sparse_embed_array_destroy(vel_x)
    call sparse_embed_array_destroy(x)
    call sparse_embed_array_destroy(tmp_vel)
    call sparse_embed_array_destroy(q)
    call sparse_embed_array_destroy(xwx)
    call sparse_embed_array_destroy(qwq)

  end subroutine internal_xl_berendsen_scuseria_ortho_update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine internal_xl_berendsen_scuseria_asym_update

    use rundat, only: md_delta_t,pub_num_spins, md_aux_dkn_t, md_aux_beren_tc,&
         pub_rootname, pub_num_kpoints, PUB_1K, md_autocorr
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_product, &
         sparse_trace
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product, sparse_embed_trace, &
         sparse_embed_array_create, sparse_embed_array_product, &
         sparse_embed_array_destroy, sparse_embed_array_axpy, &
         sparse_embed_array_scale, sparse_embed_array_num_rows, &
         sparse_embed_array_transpose, sparse_embed_array_copy
    use utils, only: utils_assert

    implicit none

    type(SPAM3_EMBED_ARRAY) :: old_x1, vel_x, vel_x_trans, x, tmp_vel, xwx, qwq, q
    real(kind=DP) :: vel, alpha, tau, trace
    real(kind=DP), dimension(pub_num_spins, PUB_1K) :: kappa
    character(len=2) :: structure_type
    integer :: label
    ! agrecocmplx
    logical :: loc_cmplx

    ! jmecmplx
    loc_cmplx = any(history(history_io)%dkn%scf(:)%iscmplx)

    structure_type = 'KS'

    ! vv: Update velocities.
    call sparse_embed_array_create(old_x1,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(vel_x,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(vel_x_trans,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(x,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(tmp_vel,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(q,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(xwx,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)
    call sparse_embed_array_create(qwq,n_spins=pub_num_spins,&
         n_kpoints=pub_num_kpoints, structure=structure_type)

    if (md_autocorr) call elec_history_print_kernel(x,'real')

    ! vv: Set initial parameters
    alpha = 0.0_dp
    vel = 0.0_dp
    kappa = 1.0_DP
    tau = md_aux_beren_tc

    ! vv: (d X / dt)_n+1 = (d X / dt)_n + (dt / 2) * ( a_n+1 + a_n)
    ! vv: Call velocities from previous step

    ! vv: Calculate accelerations a_n+1
    ! w * KS_n+1
    label = history(history_io)%label
    call elec_history_read_dkn(vel_x,'vel',label)
    call elec_history_read_dkn(x,'scf',label)
    ! -w * X_n+1
    call elec_history_read_dkn(old_x1,'init',label)
    call sparse_embed_array_axpy(x,old_x1,-1.0_DP)
    call sparse_embed_array_axpy(vel_x,x,kappa)

    call sparse_embed_array_product(tmp_vel,vel_x,old_x1)
    call sparse_embed_array_product(xwx,old_x1,tmp_vel)

    call sparse_embed_array_copy(q,old_x1)
    call sparse_embed_array_scale(q,-1.0_DP,1.0_DP) ! Q = I - KS
    call sparse_embed_array_product(tmp_vel,vel_x,q)
    call sparse_embed_array_product(qwq,q,tmp_vel)

    call sparse_embed_array_axpy(xwx,qwq,1.0_DP)
    call sparse_embed_array_axpy(vel_x,xwx,-1.0_DP)

    ! vv: Only Gamma-point for the moment
    call sparse_embed_array_transpose(vel_x_trans,vel_x)
    call sparse_embed_trace(trace,vel_x%m(1,1),vel_x_trans%m(1,1))
    vel = vel + trace
    vel = vel/(sparse_embed_array_num_rows(vel_x)*md_delta_t**2)
    call utils_assert(vel >= 0,'The temperature of the auxiliary degrees of freedom has &
         &become negative. This means that d(KS)/dt has complex eigenvalues. &
         &Try to increase the value of md_aux_beren_tc or md_aux_dkn_t', vel)

    ! vv: Scale velocities according to the mix type
    if (methods(iom)%dkn_mix_type == 'XLIS' .and. &
       (methods(iom)%dkn_count .ne. methods(iom)%dkn_init_num + 1)) then
       ! vv: Compute alpha
       alpha = alpha + sqrt(1.0_dp + md_delta_t/tau * (md_aux_dkn_t/vel - 1.0_dp))
       call utils_assert(alpha > 0,'The coupling constant alpha for the &
            &berendsen scaling of the auxiliary velocities has &
            &become negative. Try to increase the value of md_aux_beren_tc or md_aux_dkn_t',alpha)
    else
       alpha = 1.0_dp
    end if

    call sparse_embed_array_scale(vel_x,alpha)
    ! vv: Store new velocities for next md step
    call elec_history_store_dkn(vel_x,'vel')
    ! print info to external files
    call elec_history_print_log(vel,'aux',alpha)


    ! vv: Free workspace
    call sparse_embed_array_destroy(old_x1)
    call sparse_embed_array_destroy(vel_x)
    call sparse_embed_array_destroy(vel_x_trans)
    call sparse_embed_array_destroy(x)
    call sparse_embed_array_destroy(tmp_vel)
    call sparse_embed_array_destroy(q)
    call sparse_embed_array_destroy(xwx)
    call sparse_embed_array_destroy(qwq)

  end subroutine internal_xl_berendsen_scuseria_asym_update

  end subroutine elec_history_update_dkn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine elec_history_print_kernel(denskern,label)

    ! Wrapper to print density kernel matrices with sparse_show_matrix

    use comms, only: comms_bcast, pub_root_proc_id, pub_on_root
    use rundat, only: pub_rootname, md_iter_global
    use sparse_embed, only: SPAM3_EMBED_ARRAY
    use sparse, only: sparse_show_matrix
    use utils, only: utils_unit, utils_int_to_str, utils_open_unit_check, &
         utils_close_unit_check

    implicit none

    ! Arguments
    type(SPAM3_EMBED_ARRAY), intent(in) :: denskern
    character(len=*), intent(in) :: label

    ! Local variables
    character(len=32) :: iteration
    character(len=80) :: filename, autocorr_filename
    integer :: io_unit, io_stat

    iteration =  utils_int_to_str(md_iter_global)
    write(autocorr_filename,'(a,a,a,a)') trim(pub_rootname),&
         '.',trim(adjustl(label)),trim(adjustl(iteration))
    if (pub_on_root ) then
       io_unit = utils_unit()
       open(unit=io_unit, file=autocorr_filename, status='UNKNOWN', &
            access='SEQUENTIAL',form='FORMATTED',position='REWIND', &
            action='WRITE',iostat=io_stat)
       call utils_open_unit_check('elec_history_print_kernel',&
            autocorr_filename,io_stat)
    end if
    call comms_bcast(pub_root_proc_id, io_unit)
    call sparse_show_matrix(denskern%m(1,1)%p,outunit=io_unit,&
         matlab_format=.false.)
    if (pub_on_root)  then
       close(io_unit)
       call utils_close_unit_check('elec_history_print_kernel',&
            autocorr_filename,io_stat)
    end if

  end subroutine elec_history_print_kernel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine elec_history_print_log(vel,label,alpha)

    use comms, only: pub_on_root
    use rundat, only: md_iter_global, pub_rootname
    use utils, only: utils_unit, utils_open_unit_check, utils_close_unit_check,&
         utils_assert

    ! Arguments
    real(kind=DP), intent(in) :: vel
    character(len=*), intent(in) :: label
    real(kind=DP), optional, intent(in) :: alpha

    ! Local variables
    character(len=80) :: filename
    character(len=12) :: posfile,statfile
    integer :: io_unit, io_stat, iom

    iom = methods_io

    if (present(alpha)) then
       call utils_assert(methods(iom)%dkn_mix_type == 'XLI' .or. &
            methods(iom)%dkn_mix_type == 'XLIS', &
            'ALPHA is only computed when using mix_dkn_type = XLI or XLIS')
    end if

    if(md_iter_global <= methods(iom)%dkn_init_num + 2) then
       posfile    ='REWIND'
       statfile   ='REPLACE'   ! Create a new output file
    else
       posfile    ='APPEND'    ! Append to existing file
       statfile   ='UNKNOWN'
    end if

    if (pub_on_root .and. label=='real') then
       write(filename,'(a,a)') trim(pub_rootname),'.aux.log'
       io_unit = utils_unit()
       open(unit=io_unit,iostat=io_stat,file=filename,status=statfile,&
            access='SEQUENTIAL',form='FORMATTED',position=posfile,&
            action='WRITE')
       call utils_open_unit_check('elec_history_print_log',filename,&
            io_stat)

       write(io_unit,1) md_iter_global
       write(io_unit,3) vel

       close(io_unit,iostat=io_stat)
       call utils_close_unit_check('elec_history_print_log',filename,&
            io_stat)
    end if

    if (pub_on_root .and. label=='aux') then
       write(filename,'(a,a)') trim(pub_rootname),'.aux.log'
       io_unit = utils_unit()
       open(unit=io_unit,iostat=io_stat,file=filename,status=statfile,&
            access='SEQUENTIAL',form='FORMATTED',position=posfile,&
            action='WRITE')
       call utils_open_unit_check('elec_history_print_log',filename,&
            io_stat)

       if (present(alpha)) then
          write(io_unit,4) alpha
          write(io_unit,2) vel
          write(io_unit,*) ''
       else
          write(io_unit,2) vel
          write(io_unit,*) ''
       end if

       close(io_unit,iostat=io_stat)
       call utils_close_unit_check('elec_history_print_log',filename,&
            io_stat)
    end if


1  format(4x,I0.1,T25,' <--   MD_STEP')
2  format(4x,E17.8,T25,' <--   TAUX')
3  format(4x,E17.8,T25,' <--   TREAL')
4  format(4x,E17.8,T25,' <--   ALPHA')


  end subroutine elec_history_print_log

  subroutine elec_history_compute_temp(method,dkn_type,nkpts,nspins)

    use constants, only : DP
    use rundat, only : md_delta_t, md_aux_rep
    use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_array_create, &
         sparse_embed_array_destroy, sparse_embed_array_axpy, &
         sparse_embed_array_scale, sparse_embed_array_transpose, &
         sparse_embed_array_num_rows, sparse_embed_array_copy, &
         sparse_embed_trace
    use utils, only : utils_assert

    implicit none
    ! Arguments
    character(len=*), intent(in) :: method
    character(len=*), intent(in) :: dkn_type
    integer, intent(in) :: nkpts, nspins

    ! Local variables
    type(SPAM3_EMBED_ARRAY) :: new_dkn, old_dkn, vel_dkn, vel_dkn_trans
    character(len=2) :: structure_type
    real(kind=DP) :: vel, trace
    integer :: iom, label

    iom = methods_io
    vel = 0.0_dp

    if(methods(iom)%dkn_count < methods(iom)%dkn_init_num) return

    if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
       structure_type = 'K'
    else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
       structure_type = 'KS'
    end if
    call sparse_embed_array_create(old_dkn,n_kpoints=nkpts,&
         n_spins=nspins,structure=structure_type)
    call sparse_embed_array_create(new_dkn,n_kpoints=nkpts,&
         n_spins=nspins,structure=structure_type)
    call sparse_embed_array_create(vel_dkn,n_kpoints=nkpts,&
         n_spins=nspins,structure=structure_type)
    call sparse_embed_array_create(vel_dkn_trans,n_kpoints=nkpts, &
         n_spins=nspins,structure=structure_type)


    select case(dkn_type)
    CASE('real')
       ! vv: Calculate velocities of converged kernel
       !label = methods(iom)%ptr(1)
       label = history(history_io)%label
       call elec_history_read_dkn(new_dkn,'scf',label)
       label = methods(iom)%ptr(1)
       call elec_history_read_dkn(old_dkn,'scf',label)
       call sparse_embed_array_axpy(vel_dkn,new_dkn,1.0_dp)
       call sparse_embed_array_axpy(vel_dkn,old_dkn,-1.0_dp)
       call sparse_embed_array_scale(vel_dkn,1.0_dp)


       if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
          call sparse_embed_array_copy(vel_dkn_trans,vel_dkn)
       else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
          call sparse_embed_array_transpose(vel_dkn_trans,vel_dkn)
       end if
       call sparse_embed_trace(trace,vel_dkn%m(1,1),vel_dkn_trans%m(1,1))
       vel = vel + trace
       vel = vel/(sparse_embed_array_num_rows(vel_dkn)*md_delta_t**2)
       call utils_assert(vel >= 0,'The temperature of the auxiliary degrees of freedom has &
            &become negative. This means that d(KS)/dt has complex eigenvalues. &
            &Try to increase the value of md_aux_beren_tc or md_aux_dkn_t', vel)

       ! print info to externel files
       call elec_history_print_log(vel,label=dkn_type)

    CASE('aux')
       if(method == 'XLI' .or. method == 'XLIS') then
          ! vv: Calculate velocities of converged kernel
          label = history(history_io)%label
          call elec_history_read_dkn(vel_dkn,'vel',label)
       elseif (method == 'XLD') then
          ! vv: Calculate velocities of converged kernel
          !label = methods(iom)%ptr(1)
          label = history(history_io)%label
          call elec_history_read_dkn(new_dkn,'init',label)
          label = methods(iom)%ptr(1)
          call elec_history_read_dkn(old_dkn,'init',label)
          call sparse_embed_array_axpy(vel_dkn,new_dkn,1.0_dp)
          call sparse_embed_array_axpy(vel_dkn,old_dkn,-1.0_dp)
       end if

       if (md_aux_rep=='ORTHO' .or. md_aux_rep=='ortho') then
          call sparse_embed_array_copy(vel_dkn_trans,vel_dkn)
       else if (md_aux_rep=='ASYM' .or. md_aux_rep=='asym') then
          call sparse_embed_array_transpose(vel_dkn_trans,vel_dkn)
       end if
       call sparse_embed_trace(trace,vel_dkn%m(1,1),vel_dkn_trans%m(1,1))
       vel = vel + trace
       vel = vel/(sparse_embed_array_num_rows(vel_dkn)*md_delta_t**2)
       call utils_assert(vel >= 0,'The temperature of the auxiliary degrees of freedom has &
            &become negative.', vel)
       ! print info to external files
       call elec_history_print_log(vel,label=dkn_type)

    end select

    call sparse_embed_array_destroy(old_dkn)
    call sparse_embed_array_destroy(new_dkn)
    call sparse_embed_array_destroy(vel_dkn)
    call sparse_embed_array_destroy(vel_dkn_trans)

  end subroutine elec_history_compute_temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module electronic_history
