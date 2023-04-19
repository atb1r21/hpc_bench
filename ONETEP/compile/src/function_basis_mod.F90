! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   The subroutines in this file were written by                              !
!                                                                             !
!   Nicholas D.M. Hine                                                        !
!                                                                             !
!   in July 2009.                                                             !
!                                                                             !
!   Based on in part on previous code by                                      !
!                                                                             !
!   Chris-Kriton Skylaris, Arash A. Mostofi, Peter D. Haynes                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module function_basis

   use basis, only: FUNCTION_TIGHT_BOX, SPHERE
   use constants, only: DP, stdout


  implicit none

  private

  ! ndmh: This structure contains information describing the size, spheres,
  ! ndmh: tightboxes and parallel distribution of a given set of functions,
  ! ndmh: such as: the set of NGWFs, the nonlocal pseudopotential
  ! ndmh: projectors, the PAW partial waves, a DMFT correlated subspace, ...
  type FUNC_BASIS

     ! ndmh: the total number of functions in this set
     integer :: num

     ! ndmh: the number of functions on this proc
     integer :: proc_num

     ! ndmh: the number of functions on each proc (0:pub_total_num_procs-1)
     integer, allocatable :: num_on_proc(:)

     ! ndmh: the first function on each proc (0:pub_total_num_procs)
     integer, allocatable :: first_on_proc(:)

     ! ndmh: the number of functions on each atom (1:par%nat)
     integer, allocatable :: num_on_atom(:)

     ! ndmh: the first function on each atom (1:par%nat)
     integer, allocatable :: first_on_atom(:)

     ! ndmh: the proc of each function (1:num)
     integer, allocatable :: proc_of_func(:)

     ! ndmh: the atom of each function (1:num)
     integer, allocatable :: atom_of_func(:)

     ! ndmh: the maximum number of functions on any proc
     integer :: max_on_proc

     ! ndmh: the maximum number of functions on any atom
     integer :: max_on_atom

     ! ndmh: the number of ppds describing the functions on this proc
     integer :: n_ppds

     ! ndmh: the number of points per ppd (copied from CELL_INFO)
     integer :: n_pts

     ! ndmh: the number of ppds describing the functions on this proc
     ! JCW:  (multiplied by the number of points per ppd)
     integer :: size_on_grid

     ! ndmh: the number of ppds in all spheres over all procs (1:num)
     integer, allocatable :: n_ppds_sphere(:)

     ! ndmh: the maximum number of ppds in any sphere of these functions
     integer :: max_n_ppds_sphere

     ! ndmh: required size of receive buffer for functions from other procs
     integer :: func_on_grid_buffer_size

     ! ndmh: maximum sizes of tightboxes (if allocated)
     integer :: maxtight_pts1
     integer :: maxtight_pts2
     integer :: maxtight_pts3

     ! ndmh: the spheres describing the functions on this proc (1:proc_num)
     type(SPHERE), allocatable :: spheres(:)

     ! ndmh: the tightboxes describing the functions on this proc (1:proc_num)
     ! ndmh: NOT ALWAYS ALLOCATED FOR ALL FUNCTION TYPES
     type(FUNCTION_TIGHT_BOX), allocatable :: tight_boxes(:)

     ! ndmh: the tightboxes describing all the functions (1:num)
     ! ndmh: NOT ALWAYS ALLOCATED FOR ALL FUNCTION TYPES
     type(FUNCTION_TIGHT_BOX), allocatable :: all_tbs(:)

     ! ndmh: string identifying this function basis (eg 'ngwfs')
     character(len=20) :: name

     ! rc2013: integer identifying the region this basis corresponds to
     integer :: ireg

  end type FUNC_BASIS

  public :: FUNC_BASIS

  ! ndmh: special function request values and tags
  integer, parameter :: FUNCS_WAITING = -1000
  integer, parameter :: FUNCS_DONE = -2000
  integer, parameter :: req_tag = 10000000
  integer, parameter :: probe_frequency = 4
  integer, parameter :: num_buffers = 1

  ! ndmh: public subroutines

  ! ndmh: init/exit routines
  public :: function_basis_allocate
  public :: function_basis_distribute
  public :: function_basis_deallocate
  public :: function_basis_init_spheres
  public :: function_basis_copy_spheres
  public :: function_basis_init_tight_boxes
  public :: function_basis_gath_all_tbs
  public :: function_basis_init_uni_tb
  public :: function_basis_exit_uni_tb
  public :: function_basis_est_num_psincs
  public :: function_basis_estimate_size

  ! ndmh: ppd_to_tightbox (parallelised)
  public :: function_basis_ppds_to_tightbox
  public :: function_basis_tightbox_to_ppds
  public :: function_basis_ppds_to_tightbox_basic
  public :: function_basis_tightbox_to_ppds_basic
  public :: function_basis_ppds_to_sph_waves
  public :: function_basis_sph_waves_to_ppds


  interface function_basis_est_num_psincs
      module procedure function_basis_est_num_psincs
      module procedure function_basis_est_num_psincs_array
  end interface function_basis_est_num_psincs

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_allocate(fbasis,num,name,par)

    !=========================================================================!
    ! This subroutine allocates memory for the arrays in the function basis   !
    ! type whose size does not depend on parallel strategy determinations     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! fbasis (inout)            : function basis type describing these funcs  !
    ! num (input)               : total number of functions                   !
    ! name (input)              : identifying string for this function basis  !
    ! par (input)               : parallel strategy for these basis functions !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 10/07/2009                                  !
    ! pub_par global variable removed by Robert Charlton, 20/06/2018.         !
    !=========================================================================!

    use comms, only: pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: fbasis
    integer, intent(in) :: num
    character(len=*) :: name
    type(PARAL_INFO), intent(in), pointer :: par

    ! Local Variables
    integer :: ierr
    character(len=512) :: myself ! jd: Avoids using '//' in alloc_checks

    ! -------------------------------------------------------------------------

    fbasis%num = num
    fbasis%name = name
    fbasis%ireg = par%par_index

    myself = 'function_basis_allocate_'//trim(fbasis%name)

    allocate(fbasis%num_on_proc(0:pub_total_num_procs-1), stat=ierr)
    call utils_alloc_check(myself,'fbasis%num_on_proc',ierr)

    allocate(fbasis%first_on_proc(0:pub_total_num_procs),stat=ierr)
    call utils_alloc_check(myself,'fbasis%first_on_proc',ierr)

    allocate(fbasis%num_on_atom(1:par%nat),stat=ierr)
    call utils_alloc_check(myself,'fbasis%num_on_atom',ierr)

    allocate(fbasis%first_on_atom(1:par%nat),stat=ierr)
    call utils_alloc_check(myself,'fbasis%first_on_atom',ierr)

    allocate(fbasis%proc_of_func(1:num),stat=ierr)
    call utils_alloc_check(myself,'fbasis%proc_of_func',ierr)

    allocate(fbasis%atom_of_func(1:num),stat=ierr)
    call utils_alloc_check(myself,'fbasis%atom_of_func',ierr)

  end subroutine function_basis_allocate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_distribute(fbasis,elements,par)

    !=========================================================================!
    ! This subroutine allocates memory for the arrays in the function basis   !
    ! type whose size does not depend on parallel strategy determinations     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! fbasis (inout)            : function basis type describing these funcs  !
    ! num (input)               : total number of functions                   !
    ! name (input)              : identifying string for this function basis  !
    ! par (input)               : parallel strategy for these basis functions !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 10/07/2009                                  !
    ! pub_par global variable removed by Robert Charlton, 20/06/2018.         !
    !=========================================================================!

    use comms, only: pub_my_proc_id
    use hubbard_init, only: h_species
    use ion, only: ELEMENT
    use parallel_strategy, only: parallel_strategy_distr_funcs, PARAL_INFO
    use rundat, only: pub_cdft_multi_proj
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: fbasis
    type(PARAL_INFO), intent(in), pointer :: par
    type(ELEMENT), intent(in) :: elements(par%nat)

    ! Local Variables
    integer :: ierr, iat, hsp, nat, nhs
    integer, allocatable :: nfuncs(:) ! Number of funcs per atom (input order)
    character(len=512) :: myself ! jd: To avoid using '//' in (de)alloc_checks

    ! ------------------------------------------------------------------------

    myself='function_basis_distribute_'//trim(fbasis%name)

    allocate(nfuncs(par%nat),stat=ierr)
    call utils_alloc_check(myself,'nfuncs',ierr)

    ! ndmh: find number of functions on each atom
    do iat=1,par%nat
       if (fbasis%name=='ngwfs'.or.fbasis%name=='tmp_ngwfs') then
          nfuncs(iat) = elements(iat)%nfunctions
       else if (fbasis%name=='projs') then
          nfuncs(iat) = elements(iat)%nprojectors
       else if (fbasis%name=='pawpws') then
          nfuncs(iat) = elements(iat)%npawpws
       else if (fbasis%name=='hub_projs') then
          nfuncs(iat) = 0
          do hsp=1,par%num_hub_species
             if (h_species(hsp)%hub_species==elements(iat)%species_id) then
                ! gibo: modified for multi-projector cDFT (28.12.12)
                if (pub_cdft_multi_proj) then
                  nfuncs(iat) = h_species(hsp)%cdft_num_proj

                ! gibo: standard hubbard and single-ang.mom per cDFT-site run
                else
                   nfuncs(iat) = 2 * h_species(hsp)%hub_ang_mom + 1
                endif
             end if
          end do
       else if (fbasis%name=='ngwfs_cond'.or.fbasis%name=='tmp_ngwfs_cond') then
          nfuncs(iat) = elements(iat)%nfunctions_cond
       else if (fbasis%name=='ngwfs_joint') then
          nfuncs(iat) = elements(iat)%nfunctions + elements(iat)%nfunctions_cond
       else if (fbasis%name=='ngwfs_FOngwf'.or.fbasis%name=='tmp_ngwfs_FOngwf') then
          nfuncs(iat) = elements(iat)%nfunctions
       else if (fbasis%name=='ngwfs_jointph') then
          nfuncs(iat) = 2*elements(iat)%nfunctions
       else if (fbasis%name=='ngwfs_aux') then
          nfuncs(iat) = elements(iat)%nfunctions_aux
       else if (fbasis%name=='corewfs') then
          nfuncs(iat) = elements(iat)%ncorewfs
       else if (fbasis%name=='pdos_sws'.or.fbasis%name == 'tmp_pdos_sws') then
          nfuncs(iat) = elements(iat)%nfunctions_pdos
       else
          call utils_abort('Error in function_basis_distribute: &
               &unrecognised function basis identifier:'//trim(fbasis%name))
       end if
    end do

    call parallel_strategy_distr_funcs(fbasis%num, nfuncs, par%orig_atom, &
         par, fbasis%first_on_proc, fbasis%num_on_proc, &
         fbasis%first_on_atom, fbasis%num_on_atom, &
         fbasis%proc_of_func, fbasis%atom_of_func, fbasis%max_on_proc, &
         fbasis%max_on_atom)


    fbasis%proc_num = fbasis%num_on_proc(pub_my_proc_id)


    deallocate(nfuncs,stat=ierr)
    call utils_dealloc_check(myself,'nfuncs',ierr)

  end subroutine function_basis_distribute


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_deallocate(fbasis)

    !=========================================================================!
    ! This subroutine deallocates memory for the arrays in the function basis !
    ! type whose size does not depend on parallel strategy determinations     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! fbasis (inout)            : function basis type describing these funcs  !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 10/07/2009                                  !
    !=========================================================================!

    use basis, only: basis_sphere_deallocate
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: fbasis

    ! Local Variables
    integer :: ierr, ifunc
    character(len=512) :: myself ! jd: Avoids using '//' in dealloc_checks

    ! -------------------------------------------------------------------------

    myself = 'function_basis_deallocate_'//trim(fbasis%name)

    if (allocated(fbasis%all_tbs)) then
       deallocate(fbasis%all_tbs,stat=ierr)
       call utils_dealloc_check(myself,'fbasis%all_tbs',ierr)
    end if

    if (allocated(fbasis%tight_boxes)) then
       deallocate(fbasis%tight_boxes,stat=ierr)
       call utils_dealloc_check(myself,'fbasis%tight_boxes',ierr)
    end if

    if (allocated(fbasis%n_ppds_sphere)) then
       deallocate(fbasis%n_ppds_sphere,stat=ierr)
       call utils_dealloc_check(myself,'fbasis%n_ppds_sphere',ierr)
    end if

    if (allocated(fbasis%spheres)) then
       do ifunc=1,fbasis%proc_num
          call basis_sphere_deallocate(fbasis%spheres(ifunc))
       end do
       deallocate(fbasis%spheres,stat=ierr)
       call utils_dealloc_check(myself,'fbasis%spheres',ierr)
    end if

    deallocate(fbasis%atom_of_func,stat=ierr)
    call utils_dealloc_check(myself,'fbasis%atom_of_func',ierr)

    deallocate(fbasis%proc_of_func,stat=ierr)
    call utils_dealloc_check(myself,'fbasis%proc_of_func',ierr)

    deallocate(fbasis%first_on_atom,stat=ierr)
    call utils_dealloc_check(myself,'fbasis%first_on_atom',ierr)

    deallocate(fbasis%num_on_atom,stat=ierr)
    call utils_dealloc_check(myself,'fbasis%num_on_atom',ierr)

    deallocate(fbasis%first_on_proc,stat=ierr)
    call utils_dealloc_check(myself,'fbasis%first_on_proc',ierr)

    deallocate(fbasis%num_on_proc, stat=ierr)
    call utils_dealloc_check(myself,'fbasis%num_on_proc',ierr)

  end subroutine function_basis_deallocate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_copy_spheres(fbasis,fbasis_src1,fbasis_src2,par)

    !========================================================================!
    ! This subroutine copies the spheres array describing the functions of a !
    ! function basis to the spheres array of another function basis.         !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    ! fbasis (inout)            : Function basis to which the spheres are    !
    !                             to be copied.                              !
    ! fbasis_src (inout)        : Function basis from which the spheres are  !
    !                             to be copied.                              !
    ! par (input)               : parallel strategy for these basis functions!
    !------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 08/07/10.                                  !
    ! Modified to accept par as argument by Robert Charlton, 20/06/2018.     !
    !========================================================================!

    use basis, only: basis_copy_sphere
    use comms, only: comms_bcast, pub_my_proc_id, pub_total_num_procs
    use ion, only: element
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: fbasis
    type(FUNC_BASIS), intent(in) :: fbasis_src1
    type(FUNC_BASIS), intent(in), optional :: fbasis_src2
    type(PARAL_INFO), intent(in), pointer  :: par

    ! Local Variables
    integer :: ifunc_src1, ifunc_src2, ifunc
    integer :: global_ifunc, proc
    integer :: count, count_src1, count_src2
    integer :: current_offset
    integer :: iat, loc_iat
    integer :: ierr
    character(len=10) :: iat_str
    character(len=512) :: myself ! jd: Avoids using '//' in alloc_checks

    ! -------------------------------------------------------------------------

    fbasis%n_pts = fbasis_src1%n_pts

    myself = 'function_basis_copy_spheres_'//trim(fbasis%name)

    allocate(fbasis%spheres(fbasis%proc_num),stat=ierr)
    call utils_alloc_check(myself,'fbasis%spheres',ierr)

    allocate(fbasis%n_ppds_sphere(1:fbasis%num),stat=ierr)
    call utils_alloc_check(myself,'fbasis%n_ppds_sphere',ierr)

    ! ndmh: loop over all atoms on proc and copy the spheres of each one
    count = 0
    count_src1 = 0
    count_src2 = 0
    current_offset = 1
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1

       if (.not.present(fbasis_src2)) then
          if (fbasis_src1%num_on_atom(iat)/=fbasis%num_on_atom(iat)) then
             write(iat_str,'(i10)') iat
             iat_str = adjustl(trim(iat_str))
             call utils_abort('Error in function_basis_copy_spheres: &
                  &Mismatching sphere counts on atom '//iat_str)
          end if
       else
          if ((fbasis_src1%num_on_atom(iat)+fbasis_src2%num_on_atom(iat)) &
               /=fbasis%num_on_atom(iat)) then
             write(iat_str,'(i10)') iat
             iat_str = adjustl(trim(iat_str))
             call utils_abort('Error in function_basis_copy_spheres: &
                  &Mismatching sphere counts on atom '//iat_str)
          end if
       end if

       do ifunc_src1=1,fbasis_src1%num_on_atom(iat)
          count_src1 = count_src1 + 1
          count = count + 1

          ! ndmh: copy the sphere from one basis to the other
          call basis_copy_sphere(fbasis%spheres(count),&
               fbasis_src1%spheres(count_src1),current_offset)
          current_offset = current_offset + &
               fbasis%n_pts*fbasis%spheres(count)%n_ppds_sphere
       end do

       if (present(fbasis_src2)) then
          do ifunc_src2=1,fbasis_src2%num_on_atom(iat)
             count_src2 = count_src2 + 1
             count = count + 1

             ! ndmh: copy the sphere from one basis to the other
             call basis_copy_sphere(fbasis%spheres(count),&
                  fbasis_src2%spheres(count_src2),current_offset)
             current_offset = current_offset + &
                  fbasis%n_pts*fbasis%spheres(count)%n_ppds_sphere
          end do
       end if

    end do

    ! ndmh: determine the total number of ppds belonging to the spheres on this
    ! ndmh: proc
    fbasis%n_ppds = 0
    do ifunc=1,fbasis%num_on_proc(pub_my_proc_id)
       fbasis%n_ppds = fbasis%n_ppds + fbasis%spheres(ifunc)%n_ppds_sphere
    end do

    fbasis%size_on_grid = fbasis%n_ppds * fbasis%n_pts

    ! ndmh: fill this proc's entries in global n_ppds_sphere list
    do ifunc=1,fbasis%proc_num
       global_ifunc = ifunc + fbasis%first_on_proc(pub_my_proc_id) - 1
       fbasis%n_ppds_sphere(global_ifunc) = fbasis%spheres(ifunc)%n_ppds_sphere
    end do

    ! ndmh: collate numbers of ppds in all spheres from all procs
    do proc=0,pub_total_num_procs-1
       if (fbasis%num_on_proc(proc) > 0) then
          call comms_bcast(proc,fbasis%n_ppds_sphere( &
               fbasis%first_on_proc(proc)),fbasis%num_on_proc(proc))
       end if
    end do

    fbasis%max_n_ppds_sphere = maxval(fbasis%n_ppds_sphere)
    fbasis%func_on_grid_buffer_size = fbasis%max_n_ppds_sphere * &
         fbasis%n_pts * num_buffers

  end subroutine function_basis_copy_spheres


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_init_spheres(fbasis,cell,fftbox,elements_on_proc, &
       extra_radius, par)

    !========================================================================!
    ! This subroutine initialises the spheres describing the functions on    !
    ! atoms on this proc. It also returns the number of ppds needed for the  !
    ! storage of the functions of the current proc.                          !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    ! fbasis (inout)            : Function basis for which the spheres are   !
    !                             to be initialised.                         !
    ! elements_on_proc (input)  : Array of element structures of my_proc_id. !
    ! par (input)               : parallel strategy for these basis functions!
    !------------------------------------------------------------------------!
    ! Written as basis_init_ngwf_spheres by Chris-Kriton Skylaris in 2000.   !
    ! Modified by Chris-Kriton Skylaris on 26/8/2003 so that it works with   !
    ! the parallel version (ONETEP).                                         !
    ! Modified by Nicholas Hine on 15/09/2008 to copy initialisation of      !
    ! previous sphere when initialising many projectors per atom with same   !
    ! radius                                                                 !
    ! Generalised for function_basis module by Nicholas Hine on 10/07/2009   !
    ! Modified to accept par as argument by Robert Charlton, 20/06/2018.     !
    !========================================================================!

    use basis, only: basis_copy_sphere, basis_initialise_sphere
    use comms, only: comms_bcast, pub_my_proc_id, &
         pub_total_num_procs
    use fft_box, only: FFTBOX_INFO
    use ion, only: element
    use parallel_strategy, only: PARAL_INFO
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort
    ! agreco: currently setting sphere%extended = pub_extend_ngwf
    use rundat, only: pub_extend_ngwf

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: fbasis
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(PARAL_INFO), intent(in), pointer :: par
    type(ELEMENT), intent(in) :: elements_on_proc( &
         par%num_atoms_on_proc(pub_my_proc_id))
    real(kind=DP), intent(in), optional :: extra_radius

    ! Local Variables
    integer, allocatable :: ppd_list(:) ! temporary ppd index list
    integer, allocatable :: ppd_loc(:)  ! temporary ppd loc list
    real(kind=DP) :: radius
    integer :: ifunc, global_ifunc, proc
    integer :: count, current_offset, iat, loc_iat
    integer :: ierr
    character(len=512) :: myself ! jd: Avoids using '//' in (de)alloc_checks

    ! -------------------------------------------------------------------------

    myself = 'function_basis_init_spheres_'//trim(fbasis%name)

    ! Copy PPD size into fbasis structure
    fbasis%n_pts = cell%n_pts

    allocate(fbasis%spheres(fbasis%proc_num),stat=ierr)
    call utils_alloc_check(myself,'fbasis%spheres',ierr)

    allocate(fbasis%n_ppds_sphere(1:fbasis%num),stat=ierr)
    call utils_alloc_check(myself,'fbasis%n_ppds_sphere',ierr)

    ! ndmh: allocate temporary ppd location and number arrays
    allocate(ppd_loc(cell%n_ppds),stat=ierr)
    call utils_alloc_check(myself,'ppd_loc',ierr)
    allocate(ppd_list(cell%n_ppds),stat=ierr)
    call utils_alloc_check(myself,'ppd_list',ierr)

    ! ndmh: loop over all atoms and initialise the spheres of each one
    count = 0
    current_offset = 1
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1

       if (fbasis%name=='ngwfs'.or.fbasis%name=='tmp_ngwfs') then
          radius = elements_on_proc(loc_iat)%radius
       else if (fbasis%name=='projs') then
          radius = elements_on_proc(loc_iat)%max_core_radius
       else if (fbasis%name=='pawpws') then
          radius = elements_on_proc(loc_iat)%max_core_radius
       else if (fbasis%name=='corewfs') then !lr408
          radius = elements_on_proc(loc_iat)%max_core_wf_radius
       else if (fbasis%name=='hub_projs') then
          radius = elements_on_proc(loc_iat)%radius
       else if (fbasis%name=='ngwfs_cond'.or.fbasis%name=='tmp_ngwfs_cond') then
          radius = elements_on_proc(loc_iat)%radius_cond
       else if (fbasis%name=='ngwfs_aux') then
          radius = elements_on_proc(loc_iat)%radius
       else if (fbasis%name=='pdos_sws'.or.fbasis%name=='tmp_pdos_sws') then
          radius = elements_on_proc(loc_iat)%radius
       else
          call utils_abort('Error in function_basis_init_spheres: &
               &unrecognised function basis identifier: '//fbasis%name)
       end if

       if (present(extra_radius)) radius = radius + extra_radius


       do ifunc=1,fbasis%num_on_atom(iat)

          count = count + 1

          ! ndmh: copy the sphere if we are initialising more than 1 per atom
          if (ifunc > 1) then
             call basis_copy_sphere(fbasis%spheres(count),&
                  fbasis%spheres(count-1),current_offset)
          else
             ! agreco: currently setting extended property of a sphere
             ! equal to pub_extend_ngwf; in the future may want to have different
             ! localisations for different functions on same atom
             call basis_initialise_sphere(fbasis%spheres(count), &
                  elements_on_proc(loc_iat)%centre, radius, current_offset, &
                  ppd_list, ppd_loc, cell, fftbox, extended=pub_extend_ngwf)
          end if

          current_offset = current_offset + &
               fbasis%n_pts*fbasis%spheres(count)%n_ppds_sphere

          if (present(extra_radius)) fbasis%spheres(count)%radius = &
               fbasis%spheres(count)%radius - extra_radius

       end do

    end do

    ! ndmh: deallocate temporary ppd location and number arrays
    deallocate(ppd_loc,stat=ierr)
    call utils_dealloc_check(myself,'ppd_loc',ierr)
    deallocate(ppd_list,stat=ierr)
    call utils_dealloc_check(myself,'ppd_list',ierr)

    ! ndmh: determine the total number of ppds belonging to the spheres on this
    ! ndmh: proc
    fbasis%n_ppds = 0
    do ifunc=1,fbasis%num_on_proc(pub_my_proc_id)
       fbasis%n_ppds = fbasis%n_ppds + fbasis%spheres(ifunc)%n_ppds_sphere
    end do

    ! ndmh: set total size of basis
    fbasis%size_on_grid = fbasis%n_ppds * fbasis%n_pts

    ! ndmh: fill this proc's entries in global n_ppds_sphere list
    do ifunc=1,fbasis%proc_num
       global_ifunc = ifunc + fbasis%first_on_proc(pub_my_proc_id) - 1
       fbasis%n_ppds_sphere(global_ifunc) = fbasis%spheres(ifunc)%n_ppds_sphere
    end do

    ! ndmh: collate numbers of ppds in all spheres from all procs
    do proc=0,pub_total_num_procs-1
       if (fbasis%num_on_proc(proc) > 0) then
          call comms_bcast(proc,fbasis%n_ppds_sphere( &
               fbasis%first_on_proc(proc)),fbasis%num_on_proc(proc))
       end if
    end do

    fbasis%max_n_ppds_sphere = maxval(fbasis%n_ppds_sphere)
    fbasis%func_on_grid_buffer_size = fbasis%max_n_ppds_sphere * &
         fbasis%n_pts * num_buffers

  end subroutine function_basis_init_spheres


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_init_tight_boxes(fbasis,fftbox,cell,ngwf_basis,par)

    !========================================================================!
    ! This subroutine initialises the TIGHT BOXES for the functions on atoms !
    ! that belong to the current proc.                                       !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    ! fbasis (inout) : Function basis for which the tightboxes are to be     !
    !                  initialised.                                          !
    ! par (input)    : parallel strategy for these basis functions           !
    !------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000 and revised on 30/5/2001.     !
    ! Modified by Chris-Kriton Skylaris on 21/8/2003 so that it works with   !
    ! the parallel version (ONETEP).                                         !
    ! Modified by Chris-Kriton Skylaris on 17/08/2007 to work when the       !
    ! FFT box coincides with the simulation cell.                            !
    ! Adapted from ngwfs_init_tight_boxes by Nicholas Hine on 10/07/2009.    !
    ! Modified to accept par as argument by Robert Charlton, 20/06/2018.     !
    !========================================================================!

    use comms, only: comms_reduce, pub_my_proc_id
    use fft_box, only: FFTBOX_INFO
    use geometry, only: POINT, geometry_distance, local_displacement
    use parallel_strategy, only: PARAL_INFO
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_abort
    ! agreco: remove this dependency since we read
    ! directly from sphere%extended
    !use rundat, only: pub_extend_ngwf

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FUNC_BASIS), intent(inout) :: fbasis
    type(FUNC_BASIS), optional, intent(in) :: ngwf_basis
    type(PARAL_INFO), intent(in), pointer :: par

    ! Local Variables
    real(kind=DP) :: local_1,local_2,local_3,local_distance

    integer  :: loc_1,loc_2,loc_3,ppd_loc
    integer  :: a1_neighbour,a2_neighbour,a3_neighbour
    integer  :: start1,start2,start3
    integer  :: finish1,finish2,finish3
    integer  :: ppd_a1,ppd_a2,ppd_a3,map_a1a2
    integer  :: ppd_loc1 ! ppd a1-location is always in simcell zero
    integer  :: ppd_loc2
    integer  :: ppd_loc3
    integer  :: points1
    integer  :: points2
    integer  :: points3
    integer  :: ppd_diff,ifunc,ppd,ppd_count
    integer  :: ierr
    integer  :: hubbard_atom, hubbard_proj !ddor

    type(POINT) :: current_point
    character(len=512) :: myself ! jd: Avoids using '//' in alloc_checks

    ! -------------------------------------------------------------------------

    myself = 'function_basis_init_tight_boxes_'//trim(fbasis%name)

    allocate(fbasis%tight_boxes(fbasis%proc_num),stat=ierr)
    call utils_alloc_check(myself,'fbasis%tight_boxes',ierr)

    ! ddor: Use ngwf_basis%tight_boxes on this proc for the Hubbard projectors
    if (fbasis%name(1:9)=='hub_projs') then

       if (present(ngwf_basis)) then

          do hubbard_atom = par%first_atom_on_proc(pub_my_proc_id), &
               par%first_atom_on_proc(pub_my_proc_id+1) - 1

             do hubbard_proj = fbasis%first_on_atom(hubbard_atom) - &
                  fbasis%first_on_proc(pub_my_proc_id) + 1, &
                  fbasis%first_on_atom(hubbard_atom) + &
                  fbasis%num_on_atom(hubbard_atom) - &
                  fbasis%first_on_proc(pub_my_proc_id)
                ! ddor: Use tightbox of first NGWF on a Hubbard atom
                ifunc = ngwf_basis%first_on_atom(hubbard_atom) - &
                     ngwf_basis%first_on_proc(pub_my_proc_id) + 1
                fbasis%tight_boxes(hubbard_proj) = ngwf_basis%tight_boxes(ifunc)

             enddo

          enddo

       else
          call utils_abort('Error in '//myself//': fbasis%name=="hub_projs" &
               &requires optional argument ngwf_basis.')
       end if

       return

    end if

    ! ndmh: normal version
    do ifunc=1,fbasis%proc_num

       ! cks: initialise the padding of the box
       fbasis%tight_boxes(ifunc)%pad1= 0
       fbasis%tight_boxes(ifunc)%pad2= 0
       fbasis%tight_boxes(ifunc)%pad3= 0

       ! cks: initialise border ppd coordinates
       fbasis%tight_boxes(ifunc)%start_ppds1 = 2*cell%n_ppds_a1
       fbasis%tight_boxes(ifunc)%finish_ppds1 = -cell%n_ppds_a1

       fbasis%tight_boxes(ifunc)%start_ppds2 = 2*cell%n_ppds_a2
       fbasis%tight_boxes(ifunc)%finish_ppds2 = -cell%n_ppds_a2

       fbasis%tight_boxes(ifunc)%start_ppds3 = 2*cell%n_ppds_a3
       fbasis%tight_boxes(ifunc)%finish_ppds3 = -cell%n_ppds_a3

       ! cks: initialise the end points of the border ppds
       start1  =cell%n_pt1
       finish1 =1
       start2  =cell%n_pt2
       finish2 =1
       start3  =cell%n_pt3
       finish3 =1

       ! ndmh: loop over all ppds in the sphere
       do ppd_count=1,fbasis%spheres(ifunc)%n_ppds_sphere

          ppd = fbasis%spheres(ifunc)%ppd_list(1,ppd_count)
          ppd_loc = fbasis%spheres(ifunc)%ppd_list(2,ppd_count)

          ! ndmh: find ppd coords of this ppd
          ppd_a3 = (ppd-1) / (cell%n_ppds_a1*cell%n_ppds_a2) + 1
          map_a1a2 = ppd - (cell%n_ppds_a1*cell%n_ppds_a2) * (ppd_a3 - 1)
          ppd_a2 = (map_a1a2 - 1) / cell%n_ppds_a1 + 1
          ppd_a1 = map_a1a2 - cell%n_ppds_a1 * (ppd_a2-1)

          ! cks: find on which periodic cell the ppd is located
          a1_neighbour =nint(real(ppd_loc,DP)/9.0_DP)
          a2_neighbour =nint(real(ppd_loc-9*a1_neighbour,DP)/3.0_DP )
          a3_neighbour =ppd_loc-9*a1_neighbour-3*a2_neighbour

          ! cks: express the position of the ppd in terms of its periodic
          !      image in another cell (in integer ppd coordinates)
          !      leaving the function centre where it is.
          ppd_loc1 =ppd_a1 -a1_neighbour*cell%n_ppds_a1
          ppd_loc2 =ppd_a2 -a2_neighbour*cell%n_ppds_a2
          ppd_loc3 =ppd_a3 -a3_neighbour*cell%n_ppds_a3

          ! cks: redefine border ppd integer coordinates as necessary,
          !      re-initialising their end points whenever a border ppd
          !      is redefined
          if (ppd_loc1 < fbasis%tight_boxes(ifunc)%start_ppds1) then
             fbasis%tight_boxes(ifunc)%start_ppds1 =ppd_loc1
          endif
          if (ppd_loc2 < fbasis%tight_boxes(ifunc)%start_ppds2) then
             fbasis%tight_boxes(ifunc)%start_ppds2 =ppd_loc2
          endif
          if (ppd_loc3 < fbasis%tight_boxes(ifunc)%start_ppds3) then
             fbasis%tight_boxes(ifunc)%start_ppds3 =ppd_loc3
          endif
          if (ppd_loc1 > fbasis%tight_boxes(ifunc)%finish_ppds1) then
             fbasis%tight_boxes(ifunc)%finish_ppds1 =ppd_loc1
          endif
          if (ppd_loc2 > fbasis%tight_boxes(ifunc)%finish_ppds2) then
             fbasis%tight_boxes(ifunc)%finish_ppds2 =ppd_loc2
          endif
          if (ppd_loc3 > fbasis%tight_boxes(ifunc)%finish_ppds3) then
             fbasis%tight_boxes(ifunc)%finish_ppds3 =ppd_loc3
          endif

       end do

       ! cks: loop over all ppds in the tightbox
       ppd=0
       do ppd_a3=fbasis%tight_boxes(ifunc)%start_ppds3, &
            fbasis%tight_boxes(ifunc)%finish_ppds3
          points3=(ppd_a3-1)*cell%n_pt3

          do ppd_a2=fbasis%tight_boxes(ifunc)%start_ppds2, &
               fbasis%tight_boxes(ifunc)%finish_ppds2
             points2 =(ppd_a2-1)*cell%n_pt2

             do ppd_a1=fbasis%tight_boxes(ifunc)%start_ppds1, &
                  fbasis%tight_boxes(ifunc)%finish_ppds1
                points1 =(ppd_a1-1)*cell%n_pt1


                ! cks: loop over all the points of ppd
                do loc_3=0,cell%n_pt3-1
                   local_3 = real(loc_3+points3, DP)*cell%d3

                   do loc_2=0,cell%n_pt2-1
                      local_2 = real(loc_2+points2, DP)*cell%d2

                      do loc_1=0,cell%n_pt1-1
                         local_1 = real(loc_1+points1, DP)*cell%d1

                         current_point = local_displacement( &
                              cell%a1_unit, cell%a2_unit, &
                              cell%a3_unit, local_1, local_2, local_3)

                         local_distance =geometry_distance( &
                              fbasis%spheres(ifunc)%centre, current_point)

                         ! cks: set the limits, in the points of the current ppd,
                         ! of the parellelepiped that exactly inscribes the
                         ! function region into it by testing every single point
                         ! of the ppd that belongs to the function region.
                         if (  local_distance.le. &
                              ( fbasis%spheres(ifunc)%radius )  ) then

                            if ( (loc_1+1.lt.start1).and. &
                                 (ppd_a1.eq.fbasis%tight_boxes(ifunc)%start_ppds1)) &
                                 start1=loc_1+1
                            if ( (loc_1+1.gt.finish1).and.&
                                 (ppd_a1.eq.fbasis%tight_boxes(ifunc)%finish_ppds1) ) &
                                 finish1=loc_1+1
                            if ( (loc_2+1.lt.start2).and.&
                                 (ppd_a2.eq.fbasis%tight_boxes(ifunc)%start_ppds2) ) &
                                 start2=loc_2+1
                            if ( (loc_2+1.gt.finish2).and.&
                                 (ppd_a2.eq.fbasis%tight_boxes(ifunc)%finish_ppds2) ) &
                                 finish2=loc_2+1
                            if ( (loc_3+1.lt.start3).and.&
                                 (ppd_a3.eq.fbasis%tight_boxes(ifunc)%start_ppds3) ) &
                                 start3=loc_3+1
                            if ( (loc_3+1.gt.finish3).and.&
                                 (ppd_a3.eq.fbasis%tight_boxes(ifunc)%finish_ppds3) ) &
                                 finish3=loc_3+1


                         endif

                      enddo  ! loc_1
                   enddo  ! loc_2
                enddo  ! loc_3

             enddo  ! ppd_a1
          enddo  ! ppd_a2
       enddo  ! ppd_a3

       ! agreco: remove "shaving" for border PPDs if extended NGWFs are used
       ! read extended property directly from sphere%extended, so we can specify
       ! different localisations for different functions in same simulation cell
       if (fbasis%spheres(ifunc)%extended(1)) then
            start1  =1
            finish1 =cell%n_pt1
       endif

       if (fbasis%spheres(ifunc)%extended(2)) then
            start2  =1
            finish2 =cell%n_pt2
       endif

       if (fbasis%spheres(ifunc)%extended(3)) then
            start3  =1
            finish3 =cell%n_pt3
       endif

       fbasis%tight_boxes(ifunc)%START_PTS1  = start1
       fbasis%tight_boxes(ifunc)%FINISH_PTS1 = finish1
       fbasis%tight_boxes(ifunc)%START_PTS2  = start2
       fbasis%tight_boxes(ifunc)%FINISH_PTS2 = finish2
       fbasis%tight_boxes(ifunc)%START_PTS3  = start3
       fbasis%tight_boxes(ifunc)%FINISH_PTS3 = finish3

       ! cks: now finish box specification by calculating the number of (grid)
       !      points in every lattice vector direction.
       if (fftbox%coin1) then
          fbasis%tight_boxes(ifunc)%start_ppds1  = 1
          fbasis%tight_boxes(ifunc)%finish_ppds1 = cell%n_ppds_a1
          fbasis%tight_boxes(ifunc)%start_pts1   = 1
          fbasis%tight_boxes(ifunc)%finish_pts1  = cell%n_pt1
          fbasis%tight_boxes(ifunc)%tight_pts1   = cell%total_pt1
       else
          ppd_diff = fbasis%tight_boxes(ifunc)%finish_ppds1 &
               - fbasis%tight_boxes(ifunc)%start_ppds1
          fbasis%tight_boxes(ifunc)%tight_pts1 = &
               fbasis%tight_boxes(ifunc)%finish_pts1 &
               - fbasis%tight_boxes(ifunc)%start_pts1 &
               + 1 + ppd_diff*cell%n_pt1 + 2*fbasis%tight_boxes(ifunc)%pad1
       endif


       if (fftbox%coin2) then
          fbasis%tight_boxes(ifunc)%start_ppds2  = 1
          fbasis%tight_boxes(ifunc)%finish_ppds2 = cell%n_ppds_a2
          fbasis%tight_boxes(ifunc)%start_pts2   = 1
          fbasis%tight_boxes(ifunc)%finish_pts2  = cell%n_pt2
          fbasis%tight_boxes(ifunc)%tight_pts2   = cell%total_pt2
       else
          ppd_diff = fbasis%tight_boxes(ifunc)%finish_ppds2 &
               - fbasis%tight_boxes(ifunc)%start_ppds2
          fbasis%tight_boxes(ifunc)%tight_pts2 = &
               fbasis%tight_boxes(ifunc)%finish_pts2 &
               -fbasis%tight_boxes(ifunc)%start_pts2 &
               + 1 + ppd_diff*cell%n_pt2 + 2*fbasis%tight_boxes(ifunc)%pad2
       endif


       if (fftbox%coin3) then
          fbasis%tight_boxes(ifunc)%start_ppds3  = 1
          fbasis%tight_boxes(ifunc)%finish_ppds3 = cell%n_ppds_a3
          fbasis%tight_boxes(ifunc)%start_pts3   = 1
          fbasis%tight_boxes(ifunc)%finish_pts3  = cell%n_pt3
          fbasis%tight_boxes(ifunc)%tight_pts3   = cell%total_pt3
       else
          ppd_diff = fbasis%tight_boxes(ifunc)%finish_ppds3 &
               - fbasis%tight_boxes(ifunc)%start_ppds3
          fbasis%tight_boxes(ifunc)%tight_pts3 = &
               fbasis%tight_boxes(ifunc)%finish_pts3 &
               - fbasis%tight_boxes(ifunc)%start_pts3 &
               + 1 + ppd_diff*cell%n_pt3 + 2*fbasis%tight_boxes(ifunc)%pad3
       endif

    enddo

    ! ndmh: find maximum size of tightboxes across all functions in set
    fbasis%maxtight_pts1 = maxval(fbasis%tight_boxes(:)%tight_pts1)
    fbasis%maxtight_pts2 = maxval(fbasis%tight_boxes(:)%tight_pts2)
    fbasis%maxtight_pts3 = maxval(fbasis%tight_boxes(:)%tight_pts3)
    call comms_reduce('MAX',fbasis%maxtight_pts1)
    call comms_reduce('MAX',fbasis%maxtight_pts2)
    call comms_reduce('MAX',fbasis%maxtight_pts3)

  end subroutine function_basis_init_tight_boxes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_gath_all_tbs(fbasis, ngwf_basis, par)   ! input

    !=================================================================!
    ! This subroutine initialises the fbasis%all_tbs array that holds !
    ! together all tightboxes of all procs for these functions.       !
    !-----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 24/03/2004.                 !
    ! Modified by Peter Haynes on 26/5/2005.                          !
    ! Moved to function_basis_mod by Nicholas Hine, and adapted       !
    ! to take FUNC_BASIS argument 15/07/2009.                         !
    ! pub_par global variable removed by Robert Charlton, 20/06/2018. !
    !=================================================================!

    use comms, only : comms_bcast, pub_my_proc_id, pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_dealloc_check, utils_alloc_check, utils_abort

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: fbasis
    type(FUNC_BASIS), optional, intent(in) :: ngwf_basis
    type(PARAL_INFO), intent(in), pointer :: par

    ! Local variables
    integer :: local_func             ! Index of function on local proc
    integer :: global_func            ! Global index of function
    integer :: proc                   ! Proc counter
    integer :: ierr                   ! Error flag
    integer, allocatable :: buf(:,:)  ! Buffer for comms
    integer :: atom, proj             ! ddor: Used for Hubbard projector tbs
    character(len=512) :: myself ! jd: Avoids using '//' in (de)alloc_checks

    ! -------------------------------------------------------------------------

    myself = 'function_basis_gath_all_tbs_'//trim(fbasis%name)

    ! pdh: allocate memory for all tight boxes
    if (allocated(fbasis%all_tbs)) then
       if (size(fbasis%all_tbs) /= fbasis%num) then
          deallocate(fbasis%all_tbs, stat=ierr)
          call utils_dealloc_check(myself,'fbasis%all_tbs', ierr)
       end if
    end if

    if (.not. allocated(fbasis%all_tbs)) then
       allocate(fbasis%all_tbs(fbasis%num),stat=ierr)
       call utils_alloc_check(myself,'fbasis%all_tbs',ierr)
    end if

    !ddor: Use ngwf_basis%all_tbs for the Hubbard projectors
    if (fbasis%name(1:9)=='hub_projs') then

       if (present(ngwf_basis)) then

          ! ddor: loop over procs
          do proc=0,pub_total_num_procs-1

             do atom = par%first_atom_on_proc(proc), &
                  par%first_atom_on_proc(proc) + par%num_atoms_on_proc(proc) - 1

                do proj = fbasis%first_on_atom(atom), &
                     fbasis%first_on_atom(atom) + fbasis%num_on_atom(atom) - 1

                   ! ddor: Use tightbox of first NGWF on a Hubbard atom
                   fbasis%all_tbs(proj) = &
                        &ngwf_basis%all_tbs( ngwf_basis%first_on_atom(atom) )

                enddo

             enddo

          enddo

       else
          call utils_abort('Error in '//myself//': fbasis%name=="hub_projs" &
               &requires optional argument ngwf_basis.')
       end if

    else

       ! pdh: allocate comms buffer
       allocate(buf(18,maxval(fbasis%num_on_proc)),stat=ierr)
       call utils_alloc_check(myself,'buf',ierr)

       ! pdh: loop over procs
       do proc=0,pub_total_num_procs-1

          ! pdh: pack up tight boxes on local_proc into buf
          if (proc == pub_my_proc_id) then
             do local_func=1,fbasis%num_on_proc(proc)
                buf( 1,local_func) = fbasis%tight_boxes(local_func)%pad1
                buf( 2,local_func) = fbasis%tight_boxes(local_func)%pad2
                buf( 3,local_func) = fbasis%tight_boxes(local_func)%pad3
                buf( 4,local_func) = fbasis%tight_boxes(local_func)%start_ppds1
                buf( 5,local_func) = fbasis%tight_boxes(local_func)%start_ppds2
                buf( 6,local_func) = fbasis%tight_boxes(local_func)%start_ppds3
                buf( 7,local_func) = fbasis%tight_boxes(local_func)%finish_ppds1
                buf( 8,local_func) = fbasis%tight_boxes(local_func)%finish_ppds2
                buf( 9,local_func) = fbasis%tight_boxes(local_func)%finish_ppds3
                buf(10,local_func) = fbasis%tight_boxes(local_func)%start_pts1
                buf(11,local_func) = fbasis%tight_boxes(local_func)%start_pts2
                buf(12,local_func) = fbasis%tight_boxes(local_func)%start_pts3
                buf(13,local_func) = fbasis%tight_boxes(local_func)%finish_pts1
                buf(14,local_func) = fbasis%tight_boxes(local_func)%finish_pts2
                buf(15,local_func) = fbasis%tight_boxes(local_func)%finish_pts3
                buf(16,local_func) = fbasis%tight_boxes(local_func)%tight_pts1
                buf(17,local_func) = fbasis%tight_boxes(local_func)%tight_pts2
                buf(18,local_func) = fbasis%tight_boxes(local_func)%tight_pts3
             end do
          end if

          ! pdh: broadcast tight-boxes on local_proc
          call comms_bcast(proc,buf,fbasis%num_on_proc(proc)*18)

          ! pdh: loop over funcs
          local_func = 0
          do global_func=1,fbasis%num

             ! pdh: copy into global position
             if (proc == fbasis%proc_of_func(global_func)) then
                local_func = local_func + 1

                fbasis%all_tbs(global_func)%pad1         = buf( 1,local_func)
                fbasis%all_tbs(global_func)%pad2         = buf( 2,local_func)
                fbasis%all_tbs(global_func)%pad3         = buf( 3,local_func)
                fbasis%all_tbs(global_func)%start_ppds1  = buf( 4,local_func)
                fbasis%all_tbs(global_func)%start_ppds2  = buf( 5,local_func)
                fbasis%all_tbs(global_func)%start_ppds3  = buf( 6,local_func)
                fbasis%all_tbs(global_func)%finish_ppds1 = buf( 7,local_func)
                fbasis%all_tbs(global_func)%finish_ppds2 = buf( 8,local_func)
                fbasis%all_tbs(global_func)%finish_ppds3 = buf( 9,local_func)
                fbasis%all_tbs(global_func)%start_pts1   = buf(10,local_func)
                fbasis%all_tbs(global_func)%start_pts2   = buf(11,local_func)
                fbasis%all_tbs(global_func)%start_pts3   = buf(12,local_func)
                fbasis%all_tbs(global_func)%finish_pts1  = buf(13,local_func)
                fbasis%all_tbs(global_func)%finish_pts2  = buf(14,local_func)
                fbasis%all_tbs(global_func)%finish_pts3  = buf(15,local_func)
                fbasis%all_tbs(global_func)%tight_pts1   = buf(16,local_func)
                fbasis%all_tbs(global_func)%tight_pts2   = buf(17,local_func)
                fbasis%all_tbs(global_func)%tight_pts3   = buf(18,local_func)
             end if
          end do

       end do   ! loop over procs

       ! Deallocate comms buffer
       deallocate(buf,stat=ierr)
       call utils_dealloc_check(myself,'buf',ierr)

    endif

  end subroutine function_basis_gath_all_tbs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_init_uni_tb(uni_tightbox,cell, &
       ngwf_basis,ngwf_basis2) ! input

    !=================================================================!
    ! This subroutine finds the size of the universal NGWF tightbox   !
    ! and initialises its reciprocal space grid if required.          !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 19/07/2009 from bits of the         !
    ! routine basis_gath_all_tbs written by Chris-Kriton Skylaris,    !
    ! and modified by Alvaro Ruiz-Serrano.                            !
    ! Modified for multiple NGWF bases by Robert Charlton, 09/10/2017.!
    !=================================================================!

    use fft_box, only: FFTBOX_INFO, fftbox_init
    use geometry, only: point, operator(*)
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    ! lr408: Optional argument for the case of conduction states to ensure universal
    ! lr408: NGWF tightbox is big enough for both valence and conduction NGWFs
    type(FUNC_BASIS), optional, intent(in) :: ngwf_basis2(:)
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(inout) :: uni_tightbox

    ! Local Variables
    integer :: maxtight_pts1, maxtight_pts2, maxtight_pts3
    integer :: isub

    ! ars : the universal tightbox will have an odd number of points
    ! ars : so it can be used to perform FFT without problems
    maxtight_pts1 = ngwf_basis(1)%maxtight_pts1
    maxtight_pts2 = ngwf_basis(1)%maxtight_pts2
    maxtight_pts3 = ngwf_basis(1)%maxtight_pts3

    ! rc2013: find max size from all bases
    if(size(ngwf_basis) .gt. 1) then
        do isub=2,size(ngwf_basis)
            maxtight_pts1 = max(ngwf_basis(isub)%maxtight_pts1,maxtight_pts1)
            maxtight_pts2 = max(ngwf_basis(isub)%maxtight_pts2,maxtight_pts2)
            maxtight_pts3 = max(ngwf_basis(isub)%maxtight_pts3,maxtight_pts3)
        end do
    end if

    ! lr408: Check maximum size needed to fit both conduction and valence NGWFs
    if (present(ngwf_basis2)) then

       ! rc2013: loop over all bases
       do isub=1,size(ngwf_basis)
          maxtight_pts1 = max(ngwf_basis2(isub)%maxtight_pts1,maxtight_pts1)
          maxtight_pts2 = max(ngwf_basis2(isub)%maxtight_pts2,maxtight_pts2)
          maxtight_pts3 = max(ngwf_basis2(isub)%maxtight_pts3,maxtight_pts3)
       end do

    end if

    if (mod(maxtight_pts1,2).eq.0) &
         maxtight_pts1 = maxtight_pts1 + 1

    if (mod(maxtight_pts2,2).eq.0) &
         maxtight_pts2 = maxtight_pts2 + 1

    if (mod(maxtight_pts3,2).eq.0) &
         maxtight_pts3 = maxtight_pts3 + 1

    call fftbox_init(uni_tightbox,maxtight_pts1,maxtight_pts2,maxtight_pts3,cell)

  end subroutine function_basis_init_uni_tb


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_exit_uni_tb(uni_tightbox)

    !=================================================================!
    ! This subroutine deallocates arrays relating to the universal    !
    ! NGWF tightbox.                                                  !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 19/07/2009 from bits of the         !
    ! routine basis_exit by Chris-Kriton Skylaris and Quintin Hill.   !
    !=================================================================!

    use fft_box, only: FFTBOX_INFO, fftbox_exit
    use geometry, only: point, operator(*)
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(inout) :: uni_tightbox

    call fftbox_exit(uni_tightbox)

  end subroutine function_basis_exit_uni_tb


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  integer(kind=LONG) function function_basis_est_num_psincs(fbasis,cell)

    !========================================================================!
    ! This function returns an estimate of the total number of psincs within !
    ! all the spheres of a given function basis, based on the sum of their   !
    ! volumes.                                                               !
    !------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 16/12/2004.                        !
    ! Moved to function_basis by Nicholas Hine on 28/04/2010.                !
    ! Modified to return a long by Jacek Dziedzic in 2022.01 to prevent      !
    ! integer overflows for large systems (~56k atoms).                      !
    !========================================================================!

    use comms, only: comms_reduce, pub_my_proc_id
    use constants, only: DP, PI, LONG
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in) :: cell
    type(FUNC_BASIS) :: fbasis

    ! cks: Local Variables
    integer :: row ! loop counter
    real(kind=DP) :: vol  ! total volume of all NGWF spheres

    ! cks: calculate the sum of all spheres' volumes
    ! agreco: need to modify this for partially extended NGWFs?
    vol = 0.0_DP
    do row=1,fbasis%proc_num
       vol = vol +(4.0_DP/3.0_DP)*PI*(fbasis%spheres(row)%radius**3)
    enddo
    call comms_reduce('SUM', vol)

    function_basis_est_num_psincs = nint(vol/cell%weight,kind=LONG)

  end function function_basis_est_num_psincs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_estimate_size(fbasis,local_size,global_size)

    !=========================================================================!
    ! This subroutine estimates the total size of an allocated function_basis.!
    ! Test/Prototype for eventual widespread memory-usage-estimation code.    !
    !-------------------------------------------------------------------------!
    !  Arguments:                                                             !
    !    fbasis       (in): Function basis type whose size is to be estimated.!
    !    local_size  (out): Size on this proc.                                !
    !    global_size (out): Size on all procs.                                !
    !-------------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine in November 2009.          !
    ! Modified to return a long for global_size by Jacek Dziedzic in 2022.01  !
    ! to prevent integer overflows for large systems (~56k atoms).            !
    !=========================================================================!

    use comms, only: comms_allgather, pub_total_num_procs
    use constants, only: int_size, real_size, logical_size, LONG
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: fbasis
    integer, intent(out) :: local_size
    integer(kind=LONG), intent(out) :: global_size

    ! Local Variables
    integer :: distr_array_size
    integer :: spheres_array_size
    integer :: tightbox_array_size
    integer :: ifunc
    integer, allocatable :: ibuf(:)
    integer :: ierr

    ! ndmh: size of basic information (ints and characters)
    distr_array_size = 7*int_size + 20

    ! ndmh: size of allocated arrays
    distr_array_size = distr_array_size + int_size * &
         ( size(fbasis%num_on_proc) &
         + size(fbasis%first_on_proc) &
         + size(fbasis%num_on_atom) &
         + size(fbasis%first_on_atom) &
         + size(fbasis%proc_of_func) &
         + size(fbasis%atom_of_func) &
         + size(fbasis%n_ppds_sphere))

    ! ndmh: size of spheres array
    ! agreco: add contribution from extended property
    spheres_array_size = 0
    do ifunc=1,fbasis%proc_num
       spheres_array_size = spheres_array_size + real_size + 2*int_size + &
            int_size*2*fbasis%spheres(ifunc)%n_ppds_sphere + 3*logical_size
    end do

    ! ndmh: size of tightboxes arrays
    tightbox_array_size = 0
    if (allocated(fbasis%tight_boxes)) then
       tightbox_array_size = tightbox_array_size + 18*int_size*fbasis%proc_num
    end if
    if (allocated(fbasis%all_tbs)) then
       tightbox_array_size = tightbox_array_size + 18*int_size*fbasis%num
    end if

    local_size = distr_array_size + spheres_array_size + tightbox_array_size

    ! ndmh: find global size summed over all procs
    ! jd: In the absence of a comms_reduce for LONG integers, we avoid integer
    !     overflow by allgathering to an array of integers on each proc, then
    !     locally summing into global_size, which is integer(kind=LONG)
    allocate(ibuf(pub_total_num_procs),stat=ierr)
    call utils_alloc_check('function_basis_estimate_size','ibuf',ierr)
    ibuf(:) = 0
    call comms_allgather(ibuf, local_size)
    global_size = sum(int(ibuf(:),kind=LONG))
    deallocate(ibuf,stat=ierr)
    call utils_dealloc_check('function_basis_estimate_size','ibuf',ierr)

  end subroutine function_basis_estimate_size


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_tightbox_to_ppds(tb_in, &
       read_tb_orig1,read_tb_orig2,read_tb_orig3, &
       ifunc,tb_proc_id,tb_n1,tb_n2,tb_n3,funcs_on_grid,fbasis, &
       fftbox,cell,tb_start1,tb_start2,tb_start3, &
       fftbox_complex,fftbox_complex_shifted, &
       fftbox_real,sendbuf,recvbuf)

    !=========================================================================!
    ! This subroutine fills the ppds of a given function (which can be on     !
    ! any proc) with a tightbox from (potentially) another proc. The          !
    ! arguments specify the original location of the function within the      !
    ! input tightbox, which may differ from the new location of the function  !
    ! with respect to the grid.                                               !
    ! In such cases, the FFTbox is shifted by the Fourier Transform/Phase     !
    ! Shift/Fourier Transform. This routine is intended as part of the NGWF   !
    ! IO system, and as such it effectively serialises the comms. It should   !
    ! NOT be used as part of the main code, as it would be very slow and      !
    ! scale poorly with processor count as currently implemented.             !
    !-------------------------------------------------------------------------!
    !  Arguments:                                                             !
    !  tb_n1,tb_n2,tb_n3 (input)         : Tightbox sizes in 1,2,3-directions !
    !  ifunc (input)                     : Index of the function required     !
    !  tb_proc_id (input)                : Proc with the tightbox to extract  !
    !  fbasis (input)                    : Function basis describing funcs    !
    !  funcs_on_grid (input)             : PPD data of functions              !
    !  sendbuf, recvbuf (input/output)   : Buffers for send/recv of tb origin !
    !  tb_in (input/output)              : Origin tightbox array              !
    !  tb_start1,2,3 (in)                : Origin of a tightbox within FFTbox !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine in February 2010, reusing bits of restart_mod  !
    ! routines (eg restart_ngwfs_tightbox_input) originally written by        !
    ! Chris-Kriton Skylaris in March 2004.                                    !
    ! Modified by Andrea Greco in May 2015 to allow use of complex NGWFs.     !
    !=========================================================================!

    use basis, only: basis_clean_function, basis_extract_function_from_box, &
         basis_function_origin_wrt_tb, basis_phase_on_fftbox_recip
    use comms, only: comms_recv, comms_send, pub_my_proc_id
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    integer, intent(in)           :: tb_n1,tb_n2,tb_n3
    integer, intent(in)           :: tb_start1,tb_start2,tb_start3
    integer, intent(in)           :: ifunc
    integer, intent(in)           :: tb_proc_id
    type(FUNC_BASIS), intent(in)  :: fbasis
    type(FUNCTIONS), intent(inout)  :: funcs_on_grid
    real(kind=DP), intent(inout)  :: sendbuf(3), recvbuf(3)
    type(FFTBOX_DATA), intent(inout)  :: tb_in
    real(kind=DP), intent(inout)  :: read_tb_orig1, read_tb_orig2, read_tb_orig3
    type(FFTBOX_INFO), intent(in)   :: fftbox
    type(CELL_INFO), intent(in)   :: cell
    type(FFTBOX_DATA), intent(inout)  :: fftbox_real
    complex(kind=DP), intent(out) :: fftbox_complex(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    complex(kind=DP), intent(out) :: fftbox_complex_shifted( &
         fftbox%total_ld1,fftbox%total_ld2,fftbox%total_pt3)

    ! Locals
    integer :: loc_ifunc
    integer :: fft_n1, fft_n2, fft_ld1, fft_ld2, fft_n3
    real(kind=DP) :: tb_orig1, tb_orig2, tb_orig3

    ! agrecocmplx: check consistency between fftbox_real,
    ! tb_in, funcs_on_grid
    call utils_assert((tb_in%iscmplx.eqv.fftbox_real%iscmplx) .and. &
         (tb_in%iscmplx.eqv.funcs_on_grid%iscmplx), &
         'Error in function_basis_tightbox_to_ppds: incompatible argument types.')

    fft_ld1 = fftbox%total_ld1
    fft_ld2 = fftbox%total_ld2
    fft_n1 = fftbox%total_pt1
    fft_n2 = fftbox%total_pt2
    fft_n3 = fftbox%total_pt3

    ! ndmh: local proc needs the ppds of this function, so must receive them
    ! ndmh: into its tightbox buffer (and receive the read offset from the
    ! ndmh: origin) and extract the ppds.
    if (pub_my_proc_id==fbasis%proc_of_func(ifunc)) then

       loc_ifunc = ifunc - fbasis%first_on_proc(pub_my_proc_id) + 1

       ! ndmh: recv tightbox data, and origin of atom in read tightbox,
       ! ndmh: if proc holding tb is not local proc
       if (.not. (pub_my_proc_id==tb_proc_id)) then
          call comms_recv(tb_proc_id,recvbuf,tag=ifunc)
          read_tb_orig1 = recvbuf(1)
          read_tb_orig2 = recvbuf(2)
          read_tb_orig3 = recvbuf(3)
          if (tb_in%iscmplx) then
             call comms_recv(tb_proc_id, tb_in%z, tag=ifunc+fbasis%num)
          else
             call comms_recv(tb_proc_id, tb_in%d, tag=ifunc+fbasis%num)
          end if
       end if

       ! cks: find origin of NGWF wrt to tightbox in terms of number of
       ! cks: grid points in each lattice vector direction
       call basis_function_origin_wrt_tb(tb_orig1, tb_orig2, tb_orig3, &
            fbasis%spheres(loc_ifunc)%centre,fbasis%tight_boxes(loc_ifunc), &
            cell)

       ! ndmh: check if the origin has changed since last time
       if ((read_tb_orig1/=tb_orig1) .or. &
            (read_tb_orig2/=tb_orig2) .or. &
            (read_tb_orig3/=tb_orig3)) then

          ! cks: apply phase factor to translate NGWF from beginning of
          ! cks: coordinates in universal tightbox to the atomic centre
          ! cks: of NGWF
          !============== PHASE SHIFT ======================================

          ! cks: initialise complex fftbox
          fftbox_complex = (0.0_DP,0.0_DP)

          ! ndmh: deposit tightbox at (tb_start1,tb_start2,tb_start3)
          ! ndmh: within complex FFTbox
          if (tb_in%iscmplx) then
             fftbox_complex(tb_start1:tb_start1+tb_n1-1, &
                  tb_start2:tb_start2+tb_n2-1,tb_start3:tb_start3+tb_n3-1) = &
                  tb_in%z(1:tb_n1,1:tb_n2,1:tb_n3)
          else
             fftbox_complex(tb_start1:tb_start1+tb_n1-1, &
                  tb_start2:tb_start2+tb_n2-1,tb_start3:tb_start3+tb_n3-1) = &
                  cmplx(tb_in%d(1:tb_n1,1:tb_n2,1:tb_n3),0.0_DP,kind=DP)
          end if

          ! ndmh: in-place FFT to reciprocal space
          call fourier_apply_box('Coarse','Forward',fftbox_complex)

          ! ndmh: apply phase factor to move function by required shift
          call basis_phase_on_fftbox_recip(fftbox_complex_shifted, &
               fft_n1, fft_n2, fft_n3, fft_ld1, fft_ld2, &
               read_tb_orig1-tb_orig1, read_tb_orig2-tb_orig2, &
               read_tb_orig3-tb_orig3, fftbox_complex)

          ! cks: g=0 element must be real
          fftbox_complex_shifted(1,1,1) = &
               cmplx(real(fftbox_complex_shifted(1,1,1),kind=DP),0.0_DP,kind=DP)

          ! ndmh: in-place FFT back to real space
          call fourier_apply_box('Coarse','Backward',fftbox_complex_shifted)

          ! ndmh: transfer to real-valued FFTbox for extraction
          !fftbox_real = real(fftbox_complex_shifted,kind=DP)
          ! agrecocmplx: retain the full complex result if using complex NGWFs
          if (fftbox_real%iscmplx) then
              fftbox_real%z(:,:,:) = fftbox_complex_shifted
          else
              fftbox_real%d(:,:,:) = real(fftbox_complex_shifted,kind=DP)
          end if

          !=========== END PHASE SHIFT ======================================

       else ! ndmh: no need to shift function

          !fftbox_real = 0.0_DP
          ! agrecocmplx: initialise according to it being real or complex
          call data_set_to_zero(fftbox_real)

          ! ndmh: deposit tightbox at (tb_start1,tb_start2,tb_start3)
          ! ndmh: within real FFTbox
          !fftbox_real(tb_start1:tb_start1+tb_n1-1, &
          !     tb_start2:tb_start2+tb_n2-1,tb_start3:tb_start3+tb_n3-1) = &
          !     tb_in(1:tb_n1,1:tb_n2,1:tb_n3)
          ! agrecocmplx: disthinguish between real and complex case
          if (fftbox_real%iscmplx) then
              fftbox_real%z(tb_start1:tb_start1+tb_n1-1, &
                  tb_start2:tb_start2+tb_n2-1,tb_start3:tb_start3+tb_n3-1) = &
                  tb_in%z(1:tb_n1,1:tb_n2,1:tb_n3)
          else
              fftbox_real%d(tb_start1:tb_start1+tb_n1-1, &
                  tb_start2:tb_start2+tb_n2-1,tb_start3:tb_start3+tb_n3-1) = &
                  tb_in%d(1:tb_n1,1:tb_n2,1:tb_n3)
          end if

       end if

       ! ndmh: extract ppds of function from FFTbox
       !call basis_extract_function_from_box(funcs_on_grid,fft_ld1,fft_ld2, &
       !     fft_n3, fftbox_real, fbasis%spheres(loc_ifunc), &
       !     fbasis%tight_boxes(loc_ifunc), tb_start1, tb_start2, tb_start3, &
       !     fbasis%spheres(loc_ifunc)%offset, cell, fftbox)
       ! agrecocmplx: call new routine in basis_new_mod
       call basis_extract_function_from_box(funcs_on_grid, &
            fftbox_real, fbasis%spheres(loc_ifunc), &
            fbasis%tight_boxes(loc_ifunc), tb_start1, tb_start2, tb_start3, &
            fbasis%spheres(loc_ifunc)%offset, cell, fftbox)

       ! ndmh: zero points outside localisation sphere
       !call basis_clean_function(funcs_on_grid, fbasis%spheres(loc_ifunc), &
       !     fbasis%n_ppds, cell, fftbox)
       ! agrecocmplx: call modified routine in basis_new_mod
       call basis_clean_function(funcs_on_grid, fbasis%spheres(loc_ifunc), &
            cell, fftbox)

       ! ndmh: local proc does not have the ppds of this function. Do nothing
       ! ndmh: unless this is the sending proc, in which case, send the
       ! ndmh: tightbox and its position within the tightbox as written.
    else

       if (pub_my_proc_id==tb_proc_id) then
          ! ndmh: send the read function to its proc
          sendbuf(1) = read_tb_orig1
          sendbuf(2) = read_tb_orig2
          sendbuf(3) = read_tb_orig3
          call comms_send(fbasis%proc_of_func(ifunc),sendbuf,tag=ifunc)
          if (tb_in%iscmplx) then
             call comms_send(fbasis%proc_of_func(ifunc), tb_in%z, &
                  tag=ifunc+fbasis%num)
          else
             call comms_send(fbasis%proc_of_func(ifunc), tb_in%d, &
                  tag=ifunc+fbasis%num)
          end if
       endif

    endif

  end subroutine function_basis_tightbox_to_ppds


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_ppds_to_tightbox(tb_out,tb_orig1,tb_orig2,tb_orig3, &
       ifunc,tb_proc_id,funcs_on_grid,fbasis,fftbox,cell, &
       sendbuf,recvbuf)

    !=========================================================================!
    ! This subroutine places the ppds of a given function (which can be on    !
    ! any proc) into a tightbox on another proc. This routine is intended as  !
    ! part of the NGWF IO system, and as such it effectively serialises the   !
    ! comms. It should NOT be used as part of the main code, as that would    !
    ! be very slow in such a context.                                         !
    !-------------------------------------------------------------------------!
    !  Arguments:                                                             !
    !  tb_n1,tb_n2,tb_n3 (input)         : Tightbox sizes in 1,2,3-directions !
    !  ifunc (input)                     : Index of the function required     !
    !  tb_proc_id (input)                : Proc of the destination tightbox   !
    !  fbasis (input)                    : Function basis describing funcs    !
    !  funcs_on_grid (input)             : PPD data of functions              !
    !  sendbuf, recvbuf (input/output)   : Buffers for send/recv of tb origin !
    !  tb_out (output)                   : Destination tightbox array         !
    !  tb_orig1, tb_orig2, tb_orig3 (out): Origin of atom within tightbox     !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine in February 2010, reusing bits of restart_mod  !
    ! routines (eg restart_ngwfs_tightbox_output) by Chris-Kriton Skylaris.   !
    ! Modified by Andrea Greco in May 2015 to allow use of complex NGWFs.     !
    !=========================================================================!

    use basis, only: basis_copy_function_to_box, basis_function_origin_wrt_tb
    use comms, only: comms_recv, comms_send, pub_my_proc_id
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    !integer, intent(in)          :: tb_n1,tb_n2,tb_n3
    integer, intent(in)          :: ifunc
    integer, intent(in)          :: tb_proc_id
    type(FUNC_BASIS), intent(in) :: fbasis
    type(CELL_INFO), intent(in)  :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FUNCTIONS), intent(in) :: funcs_on_grid
    real(kind=DP), intent(inout) :: sendbuf(3), recvbuf(3)
    type(FFTBOX_DATA), intent(inout) :: tb_out
    real(kind=DP), intent(out)   :: tb_orig1, tb_orig2, tb_orig3

    ! Locals
    integer :: loc_ifunc

    ! agrecocmplx: check consistency between tb_in and funcs_on_grid
    call utils_assert(tb_out%iscmplx.eqv.funcs_on_grid%iscmplx, &
         'Error in function_basis_ppds_to_tightbox: incompatible argument types.')

    ! ndmh: local proc has the ppds of this function, so must copy them to its
    ! ndmh: tightbox buffer (and calculate the offset from the origin) and
    ! ndmh: send it (unless it is also the proc that is receiving it).
    if (pub_my_proc_id==fbasis%proc_of_func(ifunc)) then

       loc_ifunc = ifunc - fbasis%first_on_proc(pub_my_proc_id) + 1

       ! Initialise to zero
       !tb_out = 0.0_DP
       ! agrecocmplx: initialise according to NGWFs real or complex
       call data_set_to_zero(tb_out)

       ! Put function in tightbox
       !call basis_copy_function_to_box(tb_out,tb_n1,tb_n2,tb_n3,1,1,1, &
       !     fbasis%tight_boxes(loc_ifunc), funcs_on_grid, &
       !     fbasis%spheres(loc_ifunc), cell, fftbox)
       ! agrecocmplx: call new routine in basis_new_mod
       call basis_copy_function_to_box(tb_out,1,1,1, &
            fbasis%tight_boxes(loc_ifunc), funcs_on_grid, &
            fbasis%spheres(loc_ifunc), cell, fftbox)

       ! cks: find origin of NGWF wrt to tightbox in terms of number of
       ! cks: grid points in each lattice vector direction
       call basis_function_origin_wrt_tb(tb_orig1, tb_orig2, tb_orig3,  &
            fbasis%spheres(loc_ifunc)%centre, &
            fbasis%tight_boxes(loc_ifunc), cell)

       ! ndmh: send, if destination proc is not local
       if (.not. (pub_my_proc_id==tb_proc_id)) then
          sendbuf(1) = tb_orig1
          sendbuf(2) = tb_orig2
          sendbuf(3) = tb_orig3
          call comms_send(tb_proc_id,sendbuf,tag=ifunc)
          if (tb_out%iscmplx) then
             call comms_send(tb_proc_id, tb_out%z, tag=ifunc+fbasis%num)
          else
             call comms_send(tb_proc_id, tb_out%d, tag=ifunc+fbasis%num)
          end if
       end if

       ! ndmh: local proc does not have the ppds of this function. Do nothing
       ! ndmh: unless this is the receiving proc (in which case, receive the
       ! ndmh: tightbox.
    else

       if (pub_my_proc_id==tb_proc_id) then
          ! ndmh: receive the function from its proc
          call comms_recv(fbasis%proc_of_func(ifunc),recvbuf,tag=ifunc)
          tb_orig1 = recvbuf(1)
          tb_orig2 = recvbuf(2)
          tb_orig3 = recvbuf(3)
          if (tb_out%iscmplx) then
             call comms_recv(fbasis%proc_of_func(ifunc), tb_out%z, &
                  tag=ifunc+fbasis%num)
          else
             call comms_recv(fbasis%proc_of_func(ifunc), tb_out%d, &
                  tag=ifunc+fbasis%num)
          end if
       endif

    endif

  end subroutine function_basis_ppds_to_tightbox


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_tightbox_to_ppds_basic(tb_in, &
       read_tb_orig1,read_tb_orig2,read_tb_orig3, &
       ifunc,tb_n1,tb_n2,tb_n3,funcs_on_grid,fbasis,fftbox,cell, &
       tb_start1,tb_start2,tb_start3,fftbox_complex,fftbox_complex_shifted, &
       fftbox_real)

    !=========================================================================!
    ! This subroutine fills the ppds of a given function with a tightbox. The !
    ! arguments specify the original location of the function within the      !
    ! input tightbox, which may differ from the new location of the function  !
    ! with respect to the grid.                                               !
    ! In such cases, the FFTbox is shifted by the Fourier Transform/Phase     !
    ! Shift/Fourier Transform. This routine is intended as part of the NGWF   !
    ! IO system, and as such it effectively serialises the comms. It should   !
    ! NOT be used as part of the main code, as it would be very slow and      !
    ! scale poorly with processor count as currently implemented.             !
    !-------------------------------------------------------------------------!
    !  Arguments:                                                             !
    !  tb_n1,tb_n2,tb_n3 (input)         : Tightbox sizes in 1,2,3-directions !
    !  ifunc (input)                     : Index of the function required     !
    !  fbasis (input)                    : Function basis describing funcs    !
    !  funcs_on_grid (input)             : PPD data of functions              !
    !  tb_in (input/output)              : Origin tightbox array              !
    !  tb_start1,2,3 (in)                : Origin of a tightbox within FFTbox !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine in February 2010, reusing bits of restart_mod  !
    ! routines (eg restart_ngwfs_tightbox_input) originally written by        !
    ! Chris-Kriton Skylaris in March 2004.                                    !
    ! Modified by Gabriel Constantinescu such that the communication part is  !
    ! moved in restart_ngwfs_tightbox_input (now using MPI-IO) and this sub-  !
    ! routine only handles the actual tightbox-ppd transformation. (June 2014)!
    ! Modified by Andrea Greco in May 2015 to allow use of complex NGWFs.     !
    !=========================================================================!

    use basis, only: basis_clean_function, basis_phase_on_fftbox_recip, &
         basis_extract_function_from_box, basis_function_origin_wrt_tb
    use comms, only: pub_my_proc_id
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    integer, intent(in)           :: tb_n1,tb_n2,tb_n3
    integer, intent(in)           :: tb_start1,tb_start2,tb_start3
    integer, intent(in)           :: ifunc
    type(FUNC_BASIS), intent(in)  :: fbasis
    type(CELL_INFO), intent(in)   :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FUNCTIONS), intent(inout) :: funcs_on_grid
    type(FFTBOX_DATA), intent(inout) :: tb_in
    real(kind=DP), intent(inout)  :: read_tb_orig1, read_tb_orig2, read_tb_orig3
    type(FFTBOX_DATA), intent(inout) :: fftbox_real
    complex(kind=DP), intent(out) :: fftbox_complex(fftbox%total_ld1, &
         fftbox%total_ld2,fftbox%total_pt3)
    complex(kind=DP), intent(out) :: fftbox_complex_shifted( &
         fftbox%total_ld1,fftbox%total_ld2,fftbox%total_pt3)

    ! Locals
    integer :: loc_ifunc
    integer :: fft_n1, fft_n2, fft_ld1, fft_ld2, fft_n3
    real(kind=DP) :: tb_orig1, tb_orig2, tb_orig3

    ! agrecocmplx: check consistency between fftbox_real,
    ! tb_in and funcs_on_grid
    call utils_assert((tb_in%iscmplx.eqv.fftbox_real%iscmplx) .and. &
         (tb_in%iscmplx.eqv.funcs_on_grid%iscmplx), &
         'Error in function_basis_tightbox_to_ppds_basic: incompatible argument types.')

    fft_ld1 = fftbox%total_ld1
    fft_ld2 = fftbox%total_ld2
    fft_n1 = fftbox%total_pt1
    fft_n2 = fftbox%total_pt2
    fft_n3 = fftbox%total_pt3

    loc_ifunc = ifunc - fbasis%first_on_proc(pub_my_proc_id) + 1

    ! cks: find origin of NGWF wrt to tightbox in terms of number of
    ! cks: grid points in each lattice vector direction
    call basis_function_origin_wrt_tb(tb_orig1, tb_orig2, tb_orig3, &
         fbasis%spheres(loc_ifunc)%centre,fbasis%tight_boxes(loc_ifunc), &
         cell)

    ! ndmh: check if the origin has changed since last time
    if ((read_tb_orig1/=tb_orig1) .or. &
         (read_tb_orig2/=tb_orig2) .or. &
         (read_tb_orig3/=tb_orig3)) then

       ! cks: apply phase factor to translate NGWF from beginning of
       ! cks: coordinates in universal tightbox to the atomic centre
       ! cks: of NGWF
       !============== PHASE SHIFT ======================================

       ! cks: initialise complex fftbox
       fftbox_complex = (0.0_DP,0.0_DP)

       ! ndmh: deposit tightbox at (tb_start1,tb_start2,tb_start3)
       ! ndmh: within complex FFTbox
       if (tb_in%iscmplx) then
          fftbox_complex(tb_start1:tb_start1+tb_n1-1, &
               tb_start2:tb_start2+tb_n2-1,tb_start3:tb_start3+tb_n3-1) = &
               tb_in%z(1:tb_n1,1:tb_n2,1:tb_n3)
       else
          fftbox_complex(tb_start1:tb_start1+tb_n1-1, &
               tb_start2:tb_start2+tb_n2-1,tb_start3:tb_start3+tb_n3-1) = &
               cmplx(tb_in%d(1:tb_n1,1:tb_n2,1:tb_n3),0.0_DP,kind=DP)
       end if

       ! ndmh: in-place FFT to reciprocal space
       call fourier_apply_box('Coarse','Forward',fftbox_complex)

       ! ndmh: apply phase factor to move function by required shift
       call basis_phase_on_fftbox_recip(fftbox_complex_shifted, &
            fft_n1, fft_n2, fft_n3, fft_ld1, fft_ld2, &
            read_tb_orig1-tb_orig1, read_tb_orig2-tb_orig2, &
            read_tb_orig3-tb_orig3, fftbox_complex)

       ! cks: g=0 element must be real
       fftbox_complex_shifted(1,1,1) = &
            cmplx(real(fftbox_complex_shifted(1,1,1),kind=DP),0.0_DP,kind=DP)

       ! ndmh: in-place FFT back to real space
       call fourier_apply_box('Coarse','Backward',fftbox_complex_shifted)

       ! ndmh: transfer to real-valued FFTbox for extraction
       !fftbox_real = real(fftbox_complex_shifted,kind=DP)
       ! agrecocmplx: retain the full complex result if using complex NGWFs
       if (fftbox_real%iscmplx) then
           fftbox_real%z(:,:,:) = fftbox_complex_shifted
       else
           fftbox_real%d(:,:,:) = real(fftbox_complex_shifted,kind=DP)
       end if

       !=========== END PHASE SHIFT ======================================

    else ! ndmh: no need to shift function

       !fftbox_real = 0.0_DP
       ! agrecocmplx: initialise according to it being real or complex
       call data_set_to_zero(fftbox_real)

       ! ndmh: deposit tightbox at (tb_start1,tb_start2,tb_start3)
       ! ndmh: within real FFTbox
       !fftbox_real(tb_start1:tb_start1+tb_n1-1, &
       !     tb_start2:tb_start2+tb_n2-1,tb_start3:tb_start3+tb_n3-1) = &
       !     tb_in(1:tb_n1,1:tb_n2,1:tb_n3)
       ! agrecocmplx: disthinguish between real and complex case
       if (fftbox_real%iscmplx) then
           fftbox_real%z(tb_start1:tb_start1+tb_n1-1, &
               tb_start2:tb_start2+tb_n2-1,tb_start3:tb_start3+tb_n3-1) = &
               tb_in%z(1:tb_n1,1:tb_n2,1:tb_n3)
       else
           fftbox_real%d(tb_start1:tb_start1+tb_n1-1, &
               tb_start2:tb_start2+tb_n2-1,tb_start3:tb_start3+tb_n3-1) = &
               tb_in%d(1:tb_n1,1:tb_n2,1:tb_n3)
       end if

    end if

    ! ndmh: extract ppds of function from FFTbox
    !call basis_extract_function_from_box(funcs_on_grid,fft_ld1,fft_ld2, &
    !     fft_n3, fftbox_real, fbasis%spheres(loc_ifunc), &
    !     fbasis%tight_boxes(loc_ifunc), tb_start1, tb_start2, tb_start3, &
    !     fbasis%spheres(loc_ifunc)%offset, cell, fftbox)
    ! agrecocmplx: call new routine in basis_new_mod
    call basis_extract_function_from_box(funcs_on_grid, &
         fftbox_real, fbasis%spheres(loc_ifunc), &
         fbasis%tight_boxes(loc_ifunc), tb_start1, tb_start2, tb_start3, &
         fbasis%spheres(loc_ifunc)%offset, cell, fftbox)

    ! ndmh: zero points outside localisation sphere
    !call basis_clean_function(funcs_on_grid, fbasis%spheres(loc_ifunc), &
    !     fbasis%n_ppds, cell, fftbox)
    ! agrecocmplx: call modified routine in basis_new_mod
    call basis_clean_function(funcs_on_grid, fbasis%spheres(loc_ifunc), &
         cell, fftbox)

  end subroutine function_basis_tightbox_to_ppds_basic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_ppds_to_tightbox_basic(tb_out,tb_orig1,tb_orig2, &
       tb_orig3,ifunc,funcs_on_grid,fbasis,fftbox,cell)

    !=========================================================================!
    ! This subroutine places the ppds of a given function (which can be on    !
    ! any proc) into a tightbox on another proc. This routine is intended as  !
    ! part of the NGWF IO system, and as such it effectively serialises the   !
    ! comms. It should NOT be used as part of the main code, as that would    !
    ! be very slow in such a context.                                         !
    !-------------------------------------------------------------------------!
    !  Arguments:                                                             !
    !  tb_n1,tb_n2,tb_n3 (input)         : Tightbox sizes in 1,2,3-directions !
    !  ifunc (input)                     : Index of the function required     !
    !  fbasis (input)                    : Function basis describing funcs    !
    !  funcs_on_grid (input)             : PPD data of functions              !
    !  tb_out (output)                   : Destination tightbox array         !
    !  tb_orig1, tb_orig2, tb_orig3 (out): Origin of atom within tightbox     !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine in February 2010, reusing bits of restart_mod  !
    ! routines (eg restart_ngwfs_tightbox_output) by Chris-Kriton Skylaris.   !
    !                                                                         !
    ! Modified by Gabriel Constantinescu such that the communication part is  !
    ! moved in restart_ngwfs_tightbox_output (now using MPI-IO) and this sub- !
    ! routine only handles the actual ppd-tightbox transformation. (June 2014)!
    ! Modified by Andrea Greco in May 2015 to use complex NGWFs.              !
    !=========================================================================!

    use basis, only: basis_copy_function_to_box, basis_function_origin_wrt_tb
    use comms, only: pub_my_proc_id
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero, &
         data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    !integer, intent(in)          :: tb_n1,tb_n2,tb_n3
    integer, intent(in)          :: ifunc
    type(FUNC_BASIS), intent(in) :: fbasis
    type(CELL_INFO), intent(in)  :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FUNCTIONS), intent(in) :: funcs_on_grid
    type(FFTBOX_DATA), intent(inout) :: tb_out
    real(kind=DP), intent(out)   :: tb_orig1, tb_orig2, tb_orig3

    ! Locals
    integer :: loc_ifunc

    ! agrecocmplx: check consistency between tb_in and funcs_on_grid
    call utils_assert(tb_out%iscmplx.eqv.funcs_on_grid%iscmplx, &
         'Error in function_basis_ppds_to_tightbox_basic: incompatible argument types.')

    ! ndmh: local proc has the ppds of this function, so must copy them to its
    ! ndmh: tightbox buffer (and calculate the offset from the origin)

    loc_ifunc = ifunc - fbasis%first_on_proc(pub_my_proc_id) + 1

    ! Initialise to zero
    !tb_out = 0.0_DP
    ! agrecocmplx: initialise according to NGWFs real or complex
    call data_set_to_zero(tb_out)

    ! Put function in tightbox
    !call basis_copy_function_to_box(tb_out,tb_n1,tb_n2,tb_n3,1,1,1, &
    !     fbasis%tight_boxes(loc_ifunc), funcs_on_grid, &
    !     fbasis%spheres(loc_ifunc), cell, fftbox)
    ! agrecocmplx: call new routine in basis_new_mod
    call basis_copy_function_to_box(tb_out,1,1,1, &
         fbasis%tight_boxes(loc_ifunc), funcs_on_grid, &
         fbasis%spheres(loc_ifunc), cell, fftbox)

    ! cks: find origin of NGWF wrt to tightbox in terms of number of
    ! cks: grid points in each lattice vector direction
    call basis_function_origin_wrt_tb(tb_orig1, tb_orig2, tb_orig3,  &
         fbasis%spheres(loc_ifunc)%centre, &
         fbasis%tight_boxes(loc_ifunc), cell)


  end subroutine function_basis_ppds_to_tightbox_basic



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_ppds_to_sph_waves(sw_coeffs, maxln, maxn, &
       ifunc,sw_proc_id,funcs_on_grid,fbasis,fftbox,uni_tightbox,cell)

    !=======================================================================!
    ! Converts a function in PPD representation to a linear combination of  !
    ! spherical waves. The result is sent to the root proc for writing on a !
    ! file                                                                  !
    !-----------------------------------------------------------------------!
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Originally written by Alvaro Ruiz Serrano in January 2009 based on    !
    ! previous code by Mark Robinson.                                       !
    ! Revisited by Alvaro Ruiz Serrano in November 2011 for re-structuring  !
    ! and bug-fixing inspired on function_basis_ppds_to_tightbox.           !
    !=======================================================================!

    use basis, only: SPHERE, basis_copy_function_to_box, &
         basis_location_func_wrt_cell, basis_func_centre_wrt_fftbox
    use comms, only: comms_recv, comms_send, pub_my_proc_id
    use datatypes, only: FFTBOX_DATA, FUNCTIONS, data_fftbox_alloc, &
         data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use geometry, only: POINT
    use ion, only: element
    use rundat, only: pub_write_max_l
    use simulation_cell, only: CELL_INFO
    use spherical_wave, only: sw_recp_generate_in_tb, sw_bessel_zeros
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! ars: arguments
    integer,       intent(in   ) :: maxn
    real(kind=DP), intent(inout) :: sw_coeffs(0:pub_write_max_l,-pub_write_max_l:pub_write_max_l,1:maxn)
    integer,       intent(in   ) :: maxln(0:pub_write_max_l)
    integer, intent(in)          :: ifunc
    integer, intent(in)          :: sw_proc_id
    type(FUNC_BASIS), intent(in) :: fbasis
    type(FUNCTIONS), intent(in) :: funcs_on_grid
    type(CELL_INFO), intent(in)  :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FFTBOX_INFO), intent(in) :: uni_tightbox

    ! ars: local variables
    integer :: loc_ifunc, tb_n1, tb_n2, tb_n3, start1, start2, start3
    integer :: ll, mm, nn, xx, yy, zz, ierr
    real(kind=DP) :: radius, qnl, sw_dot_func, sw_dot_sw
    real(kind=DP), allocatable :: sw_tb(:,:,:)
    type(FFTBOX_DATA) ::  func_tb
    complex(kind=DP), allocatable :: sw_work(:,:,:)
    type(POINT) :: centre, func_centre
    ! agrecocmplx
    logical :: loc_cmplx

    ! ------------------------------------------------------------------------

    ! agrecocmplx
    loc_cmplx = funcs_on_grid%iscmplx

    ! agrecocmplx
    call utils_assert(.not. loc_cmplx, 'Error in&
         & function_basis_ppds_to_sph_waves: not ready yet for complex NGWFs.')

    ! ars: initialize swcoeff
    sw_coeffs(:,:,:)=0.0_DP


    ! ndmh: local proc has the ppds of this function, so must copy them to its
    ! ndmh: tightbox buffer (and calculate the offset from the origin) and
    ! ndmh: send it (unless it is also the proc that is receiving it).
    if (pub_my_proc_id==fbasis%proc_of_func(ifunc)) then

       ! ars: find local function, its radius, centre and tighbox
       loc_ifunc = ifunc - fbasis%first_on_proc(pub_my_proc_id) + 1
       centre = fbasis%spheres(loc_ifunc)%centre
       radius = fbasis%spheres(loc_ifunc)%radius
       tb_n1 = uni_tightbox%total_pt1
       tb_n2 = uni_tightbox%total_pt2
       tb_n3 = uni_tightbox%total_pt3

       ! ars: allocate workspace
       allocate(sw_tb(tb_n1,tb_n2,tb_n3), stat=ierr)
       call utils_alloc_check('function_basis_ppds_to_sph_waves','sw_tb',ierr)
       allocate(sw_work(tb_n1,tb_n2,tb_n3), stat=ierr)
       call utils_alloc_check('function_basis_ppds_to_sph_waves','sw_work',ierr)
       call data_fftbox_alloc(func_tb, tb_n1, tb_n2, tb_n3, .False.) !jmecmplx

       ! ars : put functions in tightbox
       call basis_copy_function_to_box(func_tb, &
            1, 1, 1, fbasis%tight_boxes(loc_ifunc), funcs_on_grid, &
            fbasis%spheres(loc_ifunc), cell, fftbox)

       ! ars : find the centre of the function wrt the tightbox
       call basis_location_func_wrt_cell(start1,start2,start3, &
            fbasis%tight_boxes(loc_ifunc), cell)

       func_centre=basis_func_centre_wrt_fftbox(centre, &
            1,1,1,start1, start2, start3, cell)


       ! loop over angular momentum l
       do ll=0,pub_write_max_l
          ! loop over azimuthal angular momentum m
          do mm=-ll,+ll
             ! loop over n
             do nn=1,maxln(ll)

                ! ars : calculate q_nl
                qnl = sw_bessel_zeros(nn,ll)/radius

                ! ars : generate the SW in the tightbox
                ! ars : the centre of the function and the spherical wave
                ! ars : must coincide in order to calculate the overlap
                call sw_recp_generate_in_tb(sw_tb,sw_work,ll,mm,&
                     qnl,radius, func_centre, uni_tightbox)


                ! ars: calculate coefficients
                sw_dot_func=0.0_DP
                sw_dot_sw=0.0_DP
                do zz=1, tb_n3
                   do yy=1, tb_n2
                      do xx=1, tb_n1
                         sw_dot_func = sw_dot_func + sw_tb(xx,yy,zz)*func_tb%d(xx,yy,zz)
                         sw_dot_sw = sw_dot_sw + sw_tb(xx,yy,zz)*sw_tb(xx,yy,zz)
                      enddo
                   enddo
                enddo
                sw_coeffs(ll,mm,nn)=sw_dot_func/sw_dot_sw

             enddo  ! end n loop
          enddo     ! end m loop
       enddo        ! end l loop

       ! ndmh: send, if destination proc is not local
       if (.not. (pub_my_proc_id==sw_proc_id)) then
          call comms_send(sw_proc_id,sw_coeffs,tag=ifunc+fbasis%num)
       end if


       ! ars: deallocate workspace
       deallocate(sw_tb,stat=ierr)
       call utils_dealloc_check('function_basis_ppds_to_sph_waves','sw_tb',ierr)
       deallocate(sw_work,stat=ierr)
       call utils_dealloc_check('function_basis_ppds_to_sph_waves','sw_work',ierr)
       call data_fftbox_dealloc(func_tb)


       ! ndmh: local proc does not have the ppds of this function. Do nothing
       ! ndmh: unless this is the receiving proc (in which case, receive the
       ! ndmh: tightbox.
    else

       if (pub_my_proc_id==sw_proc_id) then
          ! ndmh: receive the function from its proc
          call comms_recv(fbasis%proc_of_func(ifunc),sw_coeffs,tag=ifunc+fbasis%num)
       endif

    end if

  end subroutine function_basis_ppds_to_sph_waves


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine function_basis_sph_waves_to_ppds(fbasis,funcs_on_grid,fftbox, &
       uni_tightbox,cell,sw_coeffs,maxln,maxl,maxn,ifunc,sw_proc_id)

    !=======================================================================!
    ! Converts a function in spherical wave representation to PPDs.         !
    ! The result is sent to the proc that owns the function for continuing  !
    ! parallel calculation.                                                 !
    !-----------------------------------------------------------------------!
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Originally written by Alvaro Ruiz Serrano in January 2009.            !
    ! Revisited by Alvaro Ruiz Serrano in November 2011 for re-structuring  !
    ! and bug-fixing inspired on function_basis_tightbox_to_ppds.           !
    !=======================================================================!

    use basis, only: SPHERE, basis_clean_function, &
         basis_location_func_wrt_cell, basis_copy_function_to_box, &
         basis_func_centre_wrt_fftbox, basis_put_tightbox_in_fftbox, &
         basis_extract_function_from_box
    use comms, only: comms_recv, comms_send, pub_my_proc_id
    use datatypes, only: FFTBOX_DATA, FUNCTIONS, &
         data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use geometry, only: POINT
    use ion, only: element
    use simulation_cell, only: CELL_INFO
    use spherical_wave, only: sw_recp_generate_in_tb, sw_bessel_zeros
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! ars: arguments
    type(FUNC_BASIS), intent(in) :: fbasis
    type(FFTBOX_INFO), intent(in):: fftbox
    type(FFTBOX_INFO), intent(in):: uni_tightbox
    type(CELL_INFO), intent(in)  :: cell
    integer,       intent(in   ) :: maxl
    integer,       intent(in   ) :: maxn
    type(FUNCTIONS), intent(inout) :: funcs_on_grid
    real(kind=DP), intent(inout) :: sw_coeffs(0:maxl,-maxl:maxl,1:maxn)
    integer,       intent(in   ) :: maxln(0:maxl)
    integer, intent(in)          :: ifunc
    integer, intent(in)          :: sw_proc_id


    ! ars: local variables
    integer :: loc_ifunc, tb_n1, tb_n2, tb_n3, start1, start2, start3
    integer :: ll, mm, nn, ierr
    integer :: offset, npoints !, lastp
    real(kind=DP) :: radius, qnl
    ! ars: << SW FFTbox and tightbox >>
    type(FFTBOX_DATA) :: func_fftbox, func_tb
    complex(kind=DP), allocatable :: sw_work(:,:,:)
    real(kind=DP), allocatable    :: sw_tb(:,:,:)
    type(POINT) :: centre, func_centre

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~ Generate NGWFs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !################## Generate NGWFs in reciprocal space #####################!


    ! ars : 1st : The SW are generated in the universal tightbox (sw_tb)
    ! ars : 2nd : Generate the NGWFs also in the universal tightbox (ngwfs_tb)
    ! ars : 3rd : Put the NGWFs in a FFTbox in an arbitrary position. In this
    ! ars :       case, the universal tightbox is placed in the bottom-left
    ! ars :       corner of the FFTbox ([1,1,1]). The FFTbox is only used to
    ! ars :       extract ppds, so the position of the universal tightbox
    ! ars :       wrt to the FFTbox is irrelevant.
    ! ars : 4th : Extract ppds from the FFTbox to generate ngwfs_on_grid


    ! ndmh: local proc needs the ppds of this function, so must receive them
    ! ndmh: into its tightbox buffer (and receive the read offset from the
    ! ndmh: origin) and extract the ppds.
    if (pub_my_proc_id==fbasis%proc_of_func(ifunc)) then

       ! ars: find local function, its radius, centre and tighbox
       loc_ifunc = ifunc - fbasis%first_on_proc(pub_my_proc_id) + 1
       centre = fbasis%spheres(loc_ifunc)%centre
       radius = fbasis%spheres(loc_ifunc)%radius
       offset = fbasis%spheres(loc_ifunc)%offset
       npoints = fbasis%spheres(loc_ifunc)%n_ppds_sphere * fbasis%n_pts
       !lastp = offset + npoints -1
       tb_n1 = uni_tightbox%total_pt1
       tb_n2 = uni_tightbox%total_pt2
       tb_n3 = uni_tightbox%total_pt3

       ! ndmh: recv tightbox data, and origin of atom in read tightbox,
       ! ndmh: if proc holding tb is not local proc
       if (.not.(pub_my_proc_id==sw_proc_id)) then
          call comms_recv(sw_proc_id,sw_coeffs,tag=ifunc+fbasis%num)
       end if

       ! ars : allocate workspace
       call data_fftbox_alloc(func_fftbox, fftbox%total_ld1, fftbox%total_ld2, &
            fftbox%total_pt3, .False.)                                !jmecmplx
       call data_fftbox_alloc(func_tb, tb_n1, tb_n2, tb_n3, .False.) !jmecmplx
       allocate(sw_work(tb_n1, tb_n2, tb_n3),stat=ierr)
       call utils_alloc_check('restart_ngwfs_swcoeff_input','sw_work',ierr)
       allocate(sw_tb(tb_n1, tb_n2, tb_n3),stat=ierr)
       call utils_alloc_check('restart_ngwfs_swcoeff_input','sw_tb',ierr)

       ! ars : initialisation of boxes
       func_fftbox%d(:,:,:) = 0.0_DP
       func_tb%d(:,:,:) = 0.0_DP

       ! ars : the SW will be created in the NGWFs position wrt to the
       ! ars : universal tightbox
       call basis_location_func_wrt_cell( &
            start1,start2,start3,fbasis%tight_boxes(loc_ifunc),cell)

       ! ars : vector that points from the origin of the universal tightbox
       ! ars : to the centre of the localization sphere
       func_centre=basis_func_centre_wrt_fftbox( &
            fbasis%spheres(loc_ifunc)%centre,&
            1,1,1,start1, start2, start3, cell)

       ! ars : loop over angular momentum l
       do ll=0,maxl

          ! ars : loop over azimuthal angular momentum m
          do mm=-ll,+ll

             ! ars : loop over n
             do nn=1,maxln(ll)

                qnl=sw_bessel_zeros(nn,ll)/radius

                ! ars : generate the SW in the universal tightbox
                call sw_recp_generate_in_tb(sw_tb,sw_work,ll,mm,qnl,radius,&
                     func_centre,uni_tightbox)

                ! ars : generate the ngwfs in the universal tightbox
                ! ars : by a linear combination of spherical waves
                func_tb%d(:,:,:) = func_tb%d + sw_tb * sw_coeffs(ll,mm,nn)

             enddo  ! end n loop
          enddo     ! end m loop
       enddo        ! end l loop


       ! ars : put the NGWFs in the FFTbox
       call basis_put_tightbox_in_fftbox(func_fftbox, &   ! input/output
            1, 1, 1, & ! ars : position of the tightbox wrt the fftbox
            func_tb, tb_n1, tb_n2, tb_n3, 1.0_DP)! ars : factor = 1. No scaling.

       ! ars : extract ppds from FFTbox and fill ngwfs_on_grid
       call basis_extract_function_from_box(funcs_on_grid, &
            func_fftbox,fbasis%spheres(loc_ifunc),&
            fbasis%tight_boxes(loc_ifunc),&
            1,1,1, offset, cell, fftbox)



       ! ars : deallocate workspace
       deallocate(sw_tb,stat=ierr)
       call utils_dealloc_check('restart_ngwfs_swcoeff_input','sw_tb',ierr)
       deallocate(sw_work,stat=ierr)
       call utils_dealloc_check('restart_ngwfs_swcoeff_input','sw_work',ierr)
       call data_fftbox_dealloc(func_tb)
       call data_fftbox_dealloc(func_fftbox)

       ! ars : shave functions - localise within sphere
       call basis_clean_function(funcs_on_grid,&
            fbasis%spheres(loc_ifunc), cell, fftbox)


       ! ndmh: local proc does not have the ppds of this function. Do nothing
       ! ndmh: unless this is the sending proc, in which case, send the
       ! ndmh: tightbox and its position within the tightbox as written.
    else

       if (pub_my_proc_id==sw_proc_id) then
          ! ndmh: send the read function to its proc
          call comms_send(fbasis%proc_of_func(ifunc),sw_coeffs,tag=ifunc+fbasis%num)
       endif

    endif

    !#################### End generate NGWFs in reciprocal space ###############!
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End generate NGWFs ~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  end subroutine function_basis_sph_waves_to_ppds

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function function_basis_est_num_psincs_array(fbasis,cell) &
      result(est_num_psincs)

    !========================================================================!
    ! This function returns an estimate of the total number of psincs within !
    ! all the spheres of a given function basis, based on the sum of their   !
    ! volumes.                                                               !
    !------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 16/12/2004.                        !
    ! Moved to function_basis by Nicholas Hine on 28/04/2010.                !
    !========================================================================!

    use comms, only: comms_reduce
    use constants, only: DP, PI, LONG
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in) :: cell
    type(FUNC_BASIS) :: fbasis(:)

    ! Return value
    integer(kind=LONG) :: est_num_psincs

    ! cks: Local Variables
    integer :: isub

    est_num_psincs = 0_LONG
    do isub=1,size(fbasis)
        est_num_psincs = est_num_psincs + &
             function_basis_est_num_psincs(fbasis(isub), cell)
    end do

  end function function_basis_est_num_psincs_array


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end module function_basis
