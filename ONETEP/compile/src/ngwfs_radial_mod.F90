! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                   Radial NGWFs module                          !
!                                                                !
! This module initialises the initial radial functions of NGWFs. !
!----------------------------------------------------------------!
! This module was spun off from ngwfs_mod.F90 by Nicholas Hine   !
! in July 2015. It re-uses routines originally written by        !
! Chris-Kriton Skylaris and Nicholas Hine.                       !
!================================================================!
! Reference for Cartesian NGWFs:                                 !
! [1] "The Computational Modelling of Heavy Atom Chemistry",     !
!     Chris-Kriton Skylaris (PhD thesis), Cambridge, 1999.       !
! url: www.southampton.ac.uk/assets/centresresearch/documents/   !
!      compchem/skylaris_phd.pdf                                 !
!================================================================!


module ngwfs_radial

  use constants, only: DP
  use ion, only: RADIAL_NGWF_TYPE
  use ngwf_data, only: GTO_SET

  implicit none
  private

  ! cks: public functions and subroutines
  public :: ngwfs_radial_init
  public :: ngwfs_radial_exit
  public :: ngwfs_radial_set_nfunctions

  ! cks: GTO_SET for each supported element of the periodic table
  type(GTO_SET), allocatable, dimension(:) :: gbasis

  integer, parameter :: report_length = 500

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_radial_set_nfunctions(the_region)

    !=======================================================================!
    ! This subroutine                                                       !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! the_region (in/out) : contains all info relating to the subsystem        !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30/07/2012.                               !
    ! Edited for multiple subsystems by Robert Charlton, 23/01/2017.        !
    !=======================================================================!

    use atom, only: atom_set_nfunc, atom_set_nfunc_nl_cdft
    use model_type, only: REGION
    use rundat, only: pub_cdft_multi_proj

    ! Arguments
    type(REGION), intent(inout) :: the_region

    ! Local Variables
    integer :: isp,iat

    ! ndmh: reset totals of each type
    ! rc2013: par no longer global so refer to this structure specifically
    the_region%par%num_ngwfs      = 0
    the_region%par%num_ngwfs_cond = 0
    the_region%par%num_ngwfs_aux  = 0

    ! loop over atoms
    do isp=1,the_region%par%num_species

       ! ndmh: Find an example of this species
       do iat=1,the_region%par%nat
          if (the_region%elements(iat)%species_number==isp) exit
       end do

       ! check if any ion needs its number of functions set
       if (the_region%species(isp)%nfunctions==-1) &
            call atom_set_nfunc(the_region%species(isp)%nfunctions, &
                 the_region%species(isp)%ngwf_set,the_region%elements(iat))
       if (the_region%elements(iat)%nfunctions_cond==-1) &
            call atom_set_nfunc(the_region%species(isp)%nfunctions_cond, &
                 the_region%species(isp)%cond_ngwf_set, the_region%elements(iat),cond=.true.)
       if (the_region%elements(iat)%nfunctions_aux==-1) &
            call atom_set_nfunc(the_region%species(isp)%nfunctions_aux, &
                 the_region%species(isp)%aux_ngwf_set, the_region%elements(iat))

       ! ndmh: Find an example of this species
       do iat=1,the_region%par%nat
          if (the_region%elements(iat)%species_number==isp) then

             ! update elements array from species
             the_region%elements(iat)%nfunctions = the_region%species(isp)%nfunctions
             the_region%elements(iat)%nfunctions_cond = &
                 the_region%species(isp)%nfunctions_cond
             the_region%elements(iat)%nfunctions_aux = &
                 the_region%species(isp)%nfunctions_aux

             ! count number of NGWFs of each type
             the_region%par%num_ngwfs = &
                 the_region%par%num_ngwfs + the_region%elements(iat)%nfunctions
             the_region%par%num_ngwfs_cond = &
                 the_region%par%num_ngwfs_cond + the_region%elements(iat)%nfunctions_cond
             the_region%par%num_ngwfs_aux = &
                 the_region%par%num_ngwfs_aux + the_region%elements(iat)%nfunctions_aux

             ! gibo: for cdft_multi_proj, work out number of ang.mom. shells for NGWFs
             if (pub_cdft_multi_proj) call &
                 atom_set_nfunc_nl_cdft(the_region%species(isp)%nfunctions, &
                     the_region%species(isp)%ngwf_set, the_region%elements(iat))
          end if
       end do

    end do

  end subroutine ngwfs_radial_set_nfunctions


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_radial_init(ngwf_basis_name, the_region, dftb_par)

    !=======================================================================!
    ! This subroutine initialises the NGWFs of my_proc_id ppd-storage form. !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! ngwf_basis_name (input) : function basis name describing NGWFs        !
    ! the_region (input)      : contains par, elements, species for region  !
    ! dftb_par (input)        : optional parameters for DFTB                !
    !-----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000 and revised on 30/5/2001     !
    ! Modified by Chris-Kriton Skylaris on 21/8/2003 so that it works with  !
    ! the parallel version (ONETEP).                                        !
    ! Modified by Chris-Kriton Skylaris on 2/1/2004 so that it can          !
    ! initialise NGWFs from restart info from file.                         !
    ! Modified by Chris-Kriton Skylaris on 15/03/2005 so that it reads      !
    ! fireball sets only from the root proc.                                !
    ! Modified by Alvaro Ruiz Serrano on March 2009 so that it can          !
    ! initialise the NGWFs from its spherical waves expansion saved in a    !
    ! restart file                                                          !
    !-----------------------------------------------------------------------!
    ! Reorganised and mostly re-written by Nicholas Hine to use the new     !
    ! RADIAL_NGWF_TYPE in September 2010.                                   !
    ! Also now warns about unfilled angular momentum shells, September 2010.!
    ! Reorganised again in October 2011 to allow for NGWF radii to be       !
    ! overridden by choices in the pseudoatomic solver, by Nicholas Hine    !
    ! Adapted for embedding by Robert Charlton, 23/01/2017.                 !
    ! Extended to include the initialization of NGWFs in DFTB by Arihant    !
    ! Bhandari in Oct 2021.                                                 !
    !=======================================================================!

    use comms, only: comms_bcast, pub_my_proc_id, pub_on_root, &
         pub_total_num_procs
    use constants, only: DP, stdout
    use density_init, only: density_init_radial_bcast
    use dftb_parameters, only: DFTB_PARAMS
    use model_type, only: REGION
    use rundat, only: pub_permit_unusual_ngwf_count, pub_dftb
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_abort

    implicit none

    ! Arguments
    character(20) :: ngwf_basis_name
    type(REGION), intent(inout) :: the_region
    type (DFTB_PARAMS), intent(in), optional :: dftb_par

    ! Local Variables
    character(128) :: ngwf_set
    character(10) :: fmt, tmp
    integer :: orig_iat
    integer :: shell,em      ! shell counter for atomic set
    integer :: isp
    integer :: iset, nsets
    integer :: ierr
    integer :: ingwf_in_species
    integer :: num_on_atom
    integer :: proc
    integer :: ireport
    logical :: loc_cond
    logical :: shell_filled
    logical :: store_density
    real(kind=DP) :: radius
    integer, allocatable :: nreport(:)
    character(len=128), allocatable :: report(:,:)

    ! ndmh: check for cond NGWFs
    if (ngwf_basis_name==trim(adjustl('ngwfs_cond'))) then
       loc_cond = .true.
    else
       loc_cond = .false.
    end if

    ! ab: In the case of DFTB, the STO-mG data is stored in dftb_par,
    ! ab: else STO-3G data is initialized in ngwfs_sto3g_create.
    if (pub_dftb) then
       call utils_assert(present(dftb_par), 'Error: DFTB is true, however &
            & DFTB_PAR is not passed to ngwfs_radial_init')
    ! cks: allocate and initialise the gbasis array, in case it is needed
    else
       call ngwfs_sto3g_create
    end if

    nsets = 4

    ! ndmh: Create storage for the radial functions for each l
    if (.not.allocated(the_region%radial_ngwfs)) then
       allocate(the_region%radial_ngwfs(the_region%par%num_species,nsets),stat=ierr)
       call utils_alloc_check('ngwfs_radial_init','the_region%radial_ngwfs',ierr)
    end if
    allocate(report(report_length,the_region%par%num_species),stat=ierr)
    call utils_alloc_check('ngwfs_radial_init','report',ierr)
    allocate(nreport(the_region%par%num_species),stat=ierr)
    call utils_alloc_check('ngwfs_radial_init','nreport',ierr)

    ! ndmh: Evaluate or load the radial functions
    do isp=1,the_region%par%num_species

       ! ndmh: Find an example of this species
       do orig_iat=1,the_region%par%nat
          if (the_region%elements(orig_iat)%species_number==isp) exit
       end do

       ! ndmh: get NGWF string / radius / num NGWFs for this species
       if (ngwf_basis_name==trim(adjustl('ngwfs'))) then
          iset = 1
          num_on_atom = the_region%species(isp)%nfunctions
          radius = the_region%elements(orig_iat)%radius
          ngwf_set = the_region%species(isp)%ngwf_set
       else if (ngwf_basis_name=='ngwfs_cond') then
          iset = 2
          num_on_atom = the_region%species(isp)%nfunctions_cond
          radius = the_region%elements(orig_iat)%radius_cond
          ngwf_set = the_region%species(isp)%cond_ngwf_set
       else if (ngwf_basis_name=='ngwfs_aux') then
          iset = 3
          num_on_atom = the_region%species(isp)%nfunctions_aux
          radius = the_region%elements(orig_iat)%radius
          ngwf_set = the_region%species(isp)%aux_ngwf_set
       else if (ngwf_basis_name=='hub_ngwfs') then
          iset = 4
          num_on_atom = the_region%species(isp)%nfunctions
          radius = the_region%elements(orig_iat)%radius
          ngwf_set = the_region%species(isp)%ngwf_set
       end if

       ! ndmh: allocate storage for NGWF radial functions
       call ngwfs_radial_create(the_region%radial_ngwfs(isp,iset))
       the_region%radial_ngwfs(isp,iset)%nfuncs = num_on_atom

       ! ndmh: now choose how to fill the radial_ngwfs arrays for this species
       proc = modulo(isp,pub_total_num_procs)
       if (proc==pub_my_proc_id) then

          nreport(isp) = 0
          if (index(ngwf_set,'AUTO')>0) then

             ! ndmh: STO-3G orbital
             ! ab: Extended to STO-mG in case of DFTB
             call ngwfs_generate_STO_mG(the_region%radial_ngwfs(isp,iset), &
                  the_region%elements(orig_iat), dftb_par)

          else if (index(ngwf_set,'SOLVE')>0) then

             ! ndmh: pseudoatomic solver orbital
             store_density = .false.
             if ((ngwf_basis_name=='ngwfs')) store_density=.true.
             call ngwfs_solve_atom(the_region%radial_ngwfs(isp,iset), &
                  the_region%elements(orig_iat), &
                  the_region,ngwf_set,radius,num_on_atom, &
                  store_density,report(:,isp),nreport(isp))

          else

             ! ndmh: any other string is assumed to be the name of a fireball
             call ngwfs_read_fireball(the_region%radial_ngwfs(isp,iset), &
                  the_region%species(isp),radius)

          end if
       end if
    end do

    ! ndmh: share between procs, print reports, test for problems
    do isp=1,the_region%par%num_species

       ! ndmh: share the radial function data to all procs
       proc = modulo(isp,pub_total_num_procs)
       call comms_bcast(proc,the_region%radial_ngwfs(isp,iset)%nshells)
       call comms_bcast(proc,the_region%radial_ngwfs(isp,iset)%npts)
       call comms_bcast(proc,the_region%radial_ngwfs(isp,iset)%angmom)
       call comms_bcast(proc,the_region%radial_ngwfs(isp,iset)%rc)
       call comms_bcast(proc,the_region%radial_ngwfs(isp,iset)%rad)
       call comms_bcast(proc,the_region%radial_ngwfs(isp,iset)%func_real)
       call comms_bcast(proc,nreport(isp))
       call comms_bcast(proc,report(1:nreport(isp),isp),128*nreport(isp))

       ! jd: share this radial density, if any, with other procs
       if(allocated(the_region%radial_densities)) then
          call density_init_radial_bcast(the_region%radial_densities(isp),proc)
       end if

       ! Print report
       if (pub_on_root) then
          do ireport=1,nreport(isp)
             write(stdout,'(a)') trim(report(ireport,isp))
          end do
       end if

       ! ndmh: Find an example of this species
       do orig_iat=1,the_region%par%nat
          if (the_region%elements(orig_iat)%species_number==isp) exit
       end do

       ! ndmh: test that this species has no shells of NGWFs where only
       ! ndmh: part of the shell is included in the NGWF set
       ingwf_in_species = 0
       shell_filled = .false.
       do shell=1,the_region%radial_ngwfs(isp,iset)%nshells
          shell_filled = .false.
          do em=-the_region%radial_ngwfs(isp,iset)%angmom(shell), &
               the_region%radial_ngwfs(isp,iset)%angmom(shell)
             ingwf_in_species = ingwf_in_species + 1
             if (em==the_region%radial_ngwfs(isp,iset)%angmom(shell)) shell_filled = .true.
             if (ingwf_in_species.eq.the_region%radial_ngwfs(isp,iset)%nfuncs) exit
          enddo
          if (ingwf_in_species.eq.the_region%radial_ngwfs(isp,iset)%nfuncs) exit
       enddo

       ! ndmh: for final NGWF on atom, check the whole current shell is included
       ! ndmh: in the set of NGWFs. If not, print warning.
       if (.not.shell_filled.and.pub_on_root) then
          write(stdout,*)
          write(stdout,'(a,i2,3a)') 'WARNING in ngwfs_radial_init: Setting &
               &radial NGWF functions for species ',isp,' (', &
               trim(the_region%elements(orig_iat)%species_id),'):'
          write(stdout,'(a,2(i2,a))') 'Last NGWF on atom has m = ',em, &
               ' and does not complete set of m-values for l =', &
               the_region%radial_ngwfs(isp,iset)%angmom(shell),'.'
          write(tmp,'(i10)') the_region%radial_ngwfs(isp,iset)%nshells
          write(fmt,'(a,a,a)') '(a,',trim(adjustl(tmp)),'i5)'
          write(stdout,'(a)') 'This breaks symmetry and may cause unstable &
               &NGWF optimisation and poor '
          write(stdout,'(a,2(i3,a))') 'convergence. Suggest increasing number &
               &of NGWFs per atom from ',the_region%radial_ngwfs(isp,iset)%nfuncs,' to ', &
               the_region%radial_ngwfs(isp,iset)%nfuncs-em+the_region%radial_ngwfs(isp,iset)%angmom(shell), &
               '.'
          write(stdout,fmt) 'NB: All NGWF radial function angular momenta: ', &
               the_region%radial_ngwfs(isp,iset)%angmom(1:the_region%radial_ngwfs(isp,iset)%nshells)
          write(stdout,*)

          if(.not. pub_permit_unusual_ngwf_count) then
             call utils_abort('This WARNING almost certainly means you set up &
                  &your calculation wrong. Check the last lines of your output &
                  &file and make sure you understand the warning. If you are &
                  &absolutely sure you know what you are doing, add &
                  &"permit_unusual_ngwf_count T" to your input file -- this &
                  &will let you run the calculation as is, at your own risk.')
          end if
       end if

    end do

    deallocate(nreport,stat=ierr)
    call utils_dealloc_check('ngwfs_radial_init','nreport',ierr)
    deallocate(report,stat=ierr)
    call utils_dealloc_check('ngwfs_radial_init','report',ierr)

  end subroutine ngwfs_radial_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_radial_exit(the_region)

    !=======================================================================!
    ! This subroutine deallocates storage for radial NGWFs.                 !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! the_region (input)        : REGION container containing NGWF definitions !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine in July 2015.                                !
    ! Edited for multiple subsystems by Robert Charlton, 23/01/2017.        !
    !=======================================================================!

    use model_type, only: REGION
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(REGION), intent(inout) :: the_region

    ! Local Variables
    integer :: isp
    integer :: ierr
    integer :: iset, nsets

    ! ndmh: Deallocate the arrays inside the radial functions type
    nsets = size(the_region%radial_ngwfs,2)
    do iset=1,nsets
       do isp=the_region%par%num_species,1,-1
          call ngwfs_radial_destroy(the_region%radial_ngwfs(isp,iset))
       end do
    end do

    ! ndmh: Deallocate storage for the radial functions
    deallocate(the_region%radial_ngwfs,stat=ierr)
    call utils_dealloc_check('ngwfs_radial_exit','the_region%radial_ngwfs',ierr)

    ! cks: deallocate the gbasis array
    call ngwfs_sto3g_destroy

  end subroutine ngwfs_radial_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_read_fireball(radial_ngwf,current_species,radius)

    !==================================================================!
    ! This function reads from a file a set of NGWFs for a single atom.!
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.                        !
    ! Improved by Peter D. Haynes in 2004.                             !
    ! Modified by Nicholas Hine on 03/09/2010 to use RADIAL_NGWF_TYPE. !
    ! Modified by Ed Tait on 5/9/16 to use SPECIE instead of ELEMENT   !
    !==================================================================!

    use ion, only: SPECIE
    use utils, only: utils_unit, utils_open_unit_check, utils_close_unit_check,&
         utils_assert, utils_abort

    implicit none

    ! Arguments
    type(RADIAL_NGWF_TYPE), intent(inout) :: radial_ngwf
    type(SPECIE), intent(in) :: current_species
    real(kind=DP), intent(in) :: radius

    ! CKS: INTERNAL VARIABLES
    character(len=80) :: shell_name,line
    character(len=3) :: shell_string, mom_string
    character(len=128) :: ngwf_set
    character(len=80) :: current_job
    integer :: row, shell, num_shells
    integer :: ierr
    integer :: iunit

    iunit = utils_unit()
    ngwf_set = current_species%ngwf_set

    open(iunit,file=trim(ngwf_set),status='old',position='rewind', &
         iostat=ierr)
    call utils_open_unit_check('ngwfs_read_fireball',ngwf_set,ierr)

    ! cks : read angular momenta
    current_job = 'seeking marker "ANGULAR_MOMENTA"'
    do
       read(iunit,'(a80)',err=100,end=200) line
       if (index(line,'ANGULAR_MOMENTA') > 0) exit
    end do

    current_job = 'reading block "ANGULAR_MOMENTA"'
    num_shells = 0
    do
       read(iunit,'(a80)',err=100,end=200) line
       if (index(line,'END_ANGULAR_MOMENTA') > 0) exit

       num_shells = num_shells + 1
       if (num_shells <= radial_ngwf%nshells) then
          read(line,*) radial_ngwf%angmom(num_shells)
       end if
    end do

    radial_ngwf%nshells = num_shells

    ! cks : read radial positions
    rewind(iunit,iostat=ierr,err=400)
400 call utils_assert(ierr == 0, 'Error in ngwfs_read_set: rewinding file "'//&
            trim(ngwf_set)//'" failed with code ',ierr)

    current_job = 'seeking marker "RADIAL_POSITIONS"'
    do
       read(iunit,'(a80)',err=100,end=200) line
       if (index(line,'RADIAL_POSITIONS') > 0) exit
    end do

    current_job = 'reading block "RADIAL_POSITIONS"'
    do row=1,radial_ngwf%npts
       read(iunit,'(a80)',err=100,end=200) line
       if (index(line,'END RADIAL_POSITIONS') > 0) exit
       if (index(line,'END_RADIAL_POSITIONS') > 0) exit
       read(line,*,err=100,end=200) radial_ngwf%rad(row)
    end do
    radial_ngwf%npts = row - 1
    radial_ngwf%rc(:) = max(radial_ngwf%rad(radial_ngwf%npts), &
         radius)

    ! cks : read radial functions
    do shell=1,num_shells

       write(shell_string,'(i3)') shell
       write(mom_string,'(i3)') radial_ngwf%angmom(shell)
       write(shell_name,'(a80)') 'SHELL_'//trim(adjustl(shell_string))// &
            '_ANGMOM_'//adjustl(mom_string)
       shell_name = adjustl(shell_name)

       rewind(iunit,iostat=ierr,err=500)
500    call utils_assert(ierr == 0, &
            'Error in ngwfs_read_set: rewinding file "'//trim(ngwf_set)//&
            '" failed with code ',ierr)

       write(current_job,'(a80)') 'seeking marker "'//trim(shell_name)//'"'
       current_job = adjustl(current_job)
       do
          read(iunit,'(a80)',err=100,end=200) line
          if (trim(adjustl(line)) == shell_name) exit
       end do

       write(current_job,'(a80)') 'reading block "'//trim(shell_name)//'"'
       current_job = adjustl(current_job)
       do row=1,radial_ngwf%npts
          read(iunit,*,err=100,end=200) radial_ngwf%func_real(row,shell)
       end do

    end do

    close(iunit,iostat=ierr)
    call utils_close_unit_check('ngwfs_read_fireball',ngwf_set,ierr)

    return

100 call utils_abort('Error in ngwfs_read_fireball: reading file "'//&
         trim(ngwf_set)//'" failed while '//trim(current_job))

200 call utils_abort('Error in ngwfs_read_fireball: file "'//&
         trim(ngwf_set)//'" ended unexpectedly while '//trim(current_job))

  end subroutine ngwfs_read_fireball


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwfs_generate_STO_mG(radial_ngwf,current_element, dftb_par)

    !=========================================================!
    ! This subroutine initialises a radial fireball set for a !
    ! single atom to a set of STO-3G contracted Gaussian (CG) !
    ! functions taking into account also that first CG should !
    ! correspond to the first valence atomic function of the  !
    ! atom.                                                   !
    !---------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 31/03/2006.         !
    ! Modified by Nicholas Hine on 03/09/2010 to use          !
    ! the new RADIAL_NGWF_TYPE.                               !
    ! Extended by Arihant Bhandari to use STO-mG for DFTB in  !
    ! Oct 2021.                                               !
    !=========================================================!

    use constants, only: DP, PI, SAFE_DIV_EPS
    use ion, only: element
    use ngwf_data, only: pub_num_gtatoms
    use services, only: services_equally_spaced_numbers
    use utils, only: utils_assert
    use dftb_parameters, only: DFTB_PARAMS, dftb_num_species
    use rundat, only: pub_dftb, pub_dftb_cartesian_ngwfs

    implicit none

    ! Arguments
    type(RADIAL_NGWF_TYPE), intent(inout) :: radial_ngwf
    type(ELEMENT), intent(in) :: current_element
    type (DFTB_PARAMS), intent(in), optional :: dftb_par

    ! Local Variables
    real(kind=DP), parameter :: max_rad = 10.0_DP ! max extent of Gaussians
    real(kind=DP)  :: rval   ! distance in radial grid
    real(kind=DP)  :: rvsq   ! squared distance in radial grid
    integer :: shell          ! shell counter
    integer :: start_shell    ! starting STO-3G valence shell to use
    real(kind=DP) :: shell_charge   ! electron charge of core shells
    real(kind=DP) :: ion_charge     ! ionic charge of currrent atom
    integer :: fbl_shells     ! number of fireball shells constructed from STO-3G
    integer :: atnum          ! atomic number of current atom
    integer :: max_atnum      ! maximum atomic number of current atom
    integer :: row            ! radial grid point counter
    integer :: num_shells, total_shells
    integer :: iprim, jprim, nprim   ! number of primitive Gaussians
    real(kind=DP) :: lang, summ, norm, aij, cij, overlap_factor
    character(len=*), parameter :: myself = "ngwfs_generate_STO_mG"

    ! ab: In case of DFTB, the STO-mG data is stored in dftb_par,
    if (pub_dftb) then
       call utils_assert(present(dftb_par), 'Error in '//trim(myself)//&
            ': DFTB is true, however DFTB_PAR is not passed.')
    ! ab: else STO-3G data is stored in gbasis.
    else
       call utils_assert(allocated(gbasis), 'Error in '//trim(myself)//&
            ': STO-3G data not initialized.')
    end if

    atnum = current_element%atomic_number
    ion_charge = current_element%ion_charge

    ! ab: find maximum possible atomic number
    if (pub_dftb) then
       max_atnum = dftb_num_species
    else
       max_atnum = pub_num_gtatoms
    end if

    ! cks: stop if set not available
    call utils_assert(atnum <= max_atnum, 'Error: Automatic NGWF &
         &initialisation not supported for element with atomic number ', atnum)

    ! ab: total number of shells
    if (pub_dftb) then
       total_shells = dftb_par%z(atnum)%nshells
    else
       total_shells = gbasis(atnum)%nshells
    end if

    ! cks: find starting shell of STO-3G basis from which to start fireball set
    if (pub_dftb) then
       start_shell = 1
    else
       ! cks: (depends on ion charge of pseudopotential)
       start_shell = -1
       shell_charge = 0
       find_start_shell_loop: do shell=1,total_shells
          shell_charge = shell_charge + 2*(2*gbasis(atnum)%angmom(shell) + 1)
          if (shell_charge >= ( atnum-ion_charge ) ) then
             if (shell_charge == ( atnum-ion_charge ) ) start_shell = shell + 1
             if (shell_charge > ( atnum-ion_charge ) ) start_shell = shell
             exit find_start_shell_loop
          endif
       enddo find_start_shell_loop

       ! ndmh: total number of shells in STO3G set
       fbl_shells = total_shells - start_shell + 1
       call utils_assert(fbl_shells <= radial_ngwf%nshells, &
             'Error in internal_generate_sto3g_set: &
             &number of shells in sto3g set for atomic number A exceeds &
             &maximum N. The values of A and N follow. ',fbl_shells, &
             radial_ngwf%nshells)
    end if

    ! cks: now initialise the radial grid point positions
    call services_equally_spaced_numbers(radial_ngwf%rad(:), &
         0.0_DP,max_rad,radial_ngwf%npts)

    num_shells = 0
    do shell=start_shell,total_shells

       num_shells = num_shells + 1

       ! cks: initialise angular momentum for current shell
       if (pub_dftb) then
          radial_ngwf%angmom(num_shells) = dftb_par%z(atnum)%ang(shell)
       else
          radial_ngwf%angmom(num_shells) = gbasis(atnum)%angmom(shell)
       end if
       radial_ngwf%rc(num_shells) = radial_ngwf%rad(radial_ngwf%npts)

       ! cks: generate current radial contracted Gaussian
       do row=1,radial_ngwf%npts

          rval = radial_ngwf%rad(row)
          rvsq = rval**2

          ! ab: get the number of primitives
          if (pub_dftb) then
             nprim = dftb_par%z(atnum)%sto(shell)%nprim
          else
             nprim = 3
          end if

          ! ab: loop over all the primitive Gaussians
          do iprim=1,nprim
             if (pub_dftb) then
                radial_ngwf%func_real(row,num_shells) = &
                     radial_ngwf%func_real(row,num_shells) + &
                     dftb_par%z(atnum)%sto(shell)%coeff(iprim) * &
                     exp(-dftb_par%z(atnum)%sto(shell)%alpha(iprim)*rvsq)
             else
                ! cks: contract the Gaussian
                radial_ngwf%func_real(row,num_shells) = &
                     radial_ngwf%func_real(row,num_shells) + &
                     gbasis(atnum)%coco(shell, iprim) * &
                     exp(-gbasis(atnum)%expo(shell,iprim)*rvsq)
             end if
          end do

          ! cks: multiply contraction with r^ang_mom for spherical NGWFs
          ! ab:  not in the case of Cartesian NGWFs, cf. Ref[1]:(1.21)
          if (.not. pub_dftb_cartesian_ngwfs) then
             if (radial_ngwf%angmom(num_shells) > 0) then
                radial_ngwf%func_real(row,num_shells) = &
                     radial_ngwf%func_real(row,num_shells) * &
                     rval**radial_ngwf%angmom(num_shells)
             endif
          end if

       enddo

       ! ab: normalize the function in DFTB, but not for Cartesian NGWFs as
       ! ab: STO-mG parameters in GFN input files are already normalized
       if (pub_dftb .and. .not. pub_dftb_cartesian_ngwfs) then
          summ = 0.0_DP
          lang = (2.0_DP * radial_ngwf%angmom(num_shells) + 3.0_DP) / 2.0_DP
          do iprim=1,nprim
             do jprim=1,nprim
                cij= dftb_par%z(atnum)%sto(shell)%coeff(iprim) * &
                     dftb_par%z(atnum)%sto(shell)%coeff(jprim)
                aij= dftb_par%z(atnum)%sto(shell)%alpha(iprim) + &
                     dftb_par%z(atnum)%sto(shell)%alpha(jprim)
                if (aij**lang > SAFE_DIV_EPS) summ = summ + cij / aij**lang
             end do
          end do
          ! ab: norm of the Gaussian
          norm = sqrt (2.0_DP * PI * gamma(lang) * summ)

          ! ab: extra factor missing while calculating the overlap of
          ! ab: ngwfs_on_grid
          overlap_factor = 4.0_DP * PI

          radial_ngwf%func_real(:,num_shells)= &
               radial_ngwf%func_real(:,num_shells)*sqrt(overlap_factor)/norm
       end if
    enddo

    radial_ngwf%nshells = num_shells

  end subroutine ngwfs_generate_STO_mG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwfs_solve_atom(radial_ngwf,current_element,the_region, &
       elem_string,target_radius,target_nfuncs,store_density, &
       report,nreport)

    !==================================================================!
    ! This routine solves the Schrodinger equation for a single ion,   !
    ! and creates the corresponding starting radial functions for the  !
    ! NGWFs.                                                           !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 03/09/2010.                          !
    !==================================================================!

    use atom, only: ATOM_TYPE, BASIS_TYPE, atom_create, atom_destroy, &
         atom_create_basis, atom_destroy_basis, atom_solve, &
         atom_write_orbitals, atom_get_lmax, atom_split_orbital, &
         atom_polarise_orbital, atom_set_density_guess
    use constants, only: DP
    use density_init, only: density_init_radial_store
    use ion, only: ELEMENT
    use model_type, only: REGION
    use rundat, only: pub_cutoff_energy, pub_rootname, &
         pub_initial_dens_realspace, pub_write_initial_radial_ngwfs, &
         pub_num_spins
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(RADIAL_NGWF_TYPE), intent(inout) :: radial_ngwf
    type(ELEMENT), intent(in) :: current_element
    character(128), intent(in) :: elem_string
    real(kind=DP), intent(in) :: target_radius
    integer, intent(inout) :: target_nfuncs
    logical, intent(in) :: store_density
    character(len=128), intent(out) :: report(report_length)
    integer, intent(out) :: nreport
    type(REGION) :: the_region

    ! Local Variables
    type(ATOM_TYPE) :: tatom
    type(BASIS_TYPE) :: tbasis
    character(128) :: config, init
    character(80) :: filename
    character(20) :: tmp
    real(kind=DP) :: cwidth, cwidth_l(0:5), cscale
    real(kind=DP) :: rmax, rmax_l(0:5)
    real(kind=DP) :: init_charge, init_spin
    integer :: isplitnorm, nsplitnorm
    integer :: iorb, jorb, oorb
    integer :: ishell
    integer :: ipol
    integer :: pos
    integer :: lmax
    logical :: fix_derivs

    ! Parse the string following "SOLVE"
    pos = index(elem_string,'conf=')
    if (pos>0) then
       config = elem_string(pos+5:)
    else
       config = ''
    end if

    ! Read the confining potential width
    cwidth_l(:) = 3.0_DP
    pos = index(elem_string,'w=')
    if (pos>0) then
       tmp = elem_string(pos+2:)
       read(tmp,*) cwidth
       cwidth_l(:) = cwidth
    end if
    pos = index(elem_string,'ws=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) cwidth_l(0)
    end if
    pos = index(elem_string,'wp=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) cwidth_l(1)
    end if
    pos = index(elem_string,'wd=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) cwidth_l(2)
    end if
    pos = index(elem_string,'wf=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) cwidth_l(3)
    end if
    pos = index(elem_string,'wg=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) cwidth_l(4)
    end if
    pos = index(elem_string,'wh=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) cwidth_l(5)
    end if

    ! Read the cutoff radius/radii
    rmax_l(:) = target_radius
    pos = index(elem_string,'R=')
    if (pos>0) then
       tmp = elem_string(pos+2:)
       read(tmp,*) rmax
       rmax_l(:) = rmax
    end if
    pos = index(elem_string,'Rs=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) rmax_l(0)
    end if
    pos = index(elem_string,'Rp=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) rmax_l(1)
    end if
    pos = index(elem_string,'Rd=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) rmax_l(2)
    end if
    pos = index(elem_string,'Rf=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) rmax_l(3)
    end if
    pos = index(elem_string,'Rg=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) rmax_l(4)
    end if
    pos = index(elem_string,'Rh=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) rmax_l(5)
    end if

    ! Read the confining potential height
    pos = index(elem_string,'S=')
    if (pos>0) then
       tmp = elem_string(pos+2:)
       read(tmp,*) cscale
    else
       cscale = 100.0_DP
    end if

    ! Parse the string following "INIT"
    pos = index(elem_string,'INIT ')
    if (pos>0) then
       init = elem_string(pos+5:)
    else
       init = ''
    end if
    pos = index(init,'CHARGE=')
    if (pos>0) then
       tmp = init(pos+7:)
       read(tmp,*) init_charge
    else
       init_charge = 0.0_DP
    end if
    pos = index(init,'SPIN=')
    if (pos>0) then
       tmp = init(pos+5:)
       read(tmp,*) init_spin
    else
       init_spin = 0.0_DP
    end if

    ! Parse the string following "FD"
    pos = index(elem_string,'FD')
    if (pos>0) then
       fix_derivs = .true.
    else
       fix_derivs = .false.
    end if

    ! jcap: set the region counter in tatom
    tatom%ireg=the_region%reg_count

    ! Set up the basis and the atom
    call atom_get_lmax(lmax,current_element,target_nfuncs,config, &
         the_region%pseudo_sp,the_region%paw_sp)
    call atom_create_basis(tbasis,rmax_l,lmax,6.0_DP*pub_cutoff_energy,fix_derivs)
    call atom_create(tatom,tbasis,current_element,target_nfuncs, &
         config,cwidth_l,cscale,the_region%pseudo_sp,the_region%paw_sp)

    ! Solve the atom
    call atom_solve(tatom,tbasis,report,nreport,config)

    ! Copy the atom solution into the radial NGWF
    if (tbasis%npts > size(radial_ngwf%rad)) then
       call utils_abort('Error in ngwfs_solve_atom: Atom grid is larger than &
            &RADIAL_NGWF grid')
    end if
    radial_ngwf%npts = tbasis%npts
    radial_ngwf%rad(1:radial_ngwf%npts) = tbasis%rad(1:tbasis%npts)
    radial_ngwf%nshells = 0

    ! Find orbitals which need polarising, and generate polarisation orbitals
    jorb = tatom%norbs - sum(tatom%polarise) + 1
    do iorb=1,tatom%norbs - sum(tatom%polarise)
       ! If splitnorm == -1, ignore this orbital
       if (tatom%splitnorm(1,iorb)==-1.0_DP) cycle

       ! Generate the polarisation orbital
       do ipol=1,tatom%polarise(iorb)
          if (tatom%orb_ang_mom(iorb)+ipol/=tatom%orb_ang_mom(jorb)) then
             call utils_abort('Error polarising orbital: angular momenta do not&
                 & match. Orbital: ', iorb)
          end if
          if (ipol==1) then
             oorb = iorb
          else
             oorb = jorb - 1
          end if
          write(report(nreport+1),'(a,i3,a,i2,a,i2,a)') &
              'Polarising orbital',oorb,' to generate l=', &
               tatom%orb_ang_mom(oorb) + 1,' function (orbital',jorb,')'
          nreport = nreport + 1
          call atom_polarise_orbital(tatom%psi_r(1:tbasis%npts,jorb), &
               tatom%psi_r(1:tbasis%npts,oorb),tatom%eigs(iorb), &
               tatom%orb_ang_mom(oorb),tatom,tbasis)
          jorb = jorb + 1
       end do
    end do

    ! Find orbitals which need splitting, and transfer them to radial NGWF type
    ishell = 0
    do iorb=1,tatom%norbs
       ! If splitnorm == -1, ignore this orbital
       if (tatom%splitnorm(1,iorb)==-1.0_DP) cycle

       ! Find number of target functions to split out of this orbital
       nsplitnorm = count(tatom%splitnorm(:,iorb)>0.0_DP) + 1
       ! Copy the orbital into the radial NGWF nsplitnorm times over
       do isplitnorm=1,nsplitnorm
          ishell = ishell + 1
          radial_ngwf%nshells = radial_ngwf%nshells + 1
          radial_ngwf%angmom(ishell) = tatom%orb_ang_mom(iorb)
          radial_ngwf%rc(ishell) = tbasis%rmax_l(tatom%orb_ang_mom(iorb))
          radial_ngwf%func_real(1:radial_ngwf%npts,ishell) = &
               tatom%psi_r(1:tbasis%npts,iorb)
          if (isplitnorm>1) then
             write(report(nreport+1),'(a,i3,a,f12.9)') &
                 'Splitting orbital ',iorb,', splitnorm=', &
                  tatom%splitnorm(isplitnorm-1,iorb)
             nreport = nreport +1
          end if
       end do
       ! Now do the splitting
       call atom_split_orbital( &
            radial_ngwf%func_real(:,ishell-nsplitnorm+1:ishell), &
            tatom%psi_r(1:tbasis%npts,iorb),nsplitnorm, &
            tatom%splitnorm(:,iorb),tatom%orb_ang_mom(iorb), &
            tatom%work,tbasis)

    end do

    ! Copy the atom density into the radial density storage
    if (pub_initial_dens_realspace.and.store_density) then
       call atom_set_density_guess(tatom,tbasis,init_charge,init_spin,report,nreport)
       call density_init_radial_store(the_region,current_element%species_number, &
            tbasis%npts,tbasis%rad(1:tbasis%npts), &
            tatom%guess_den(1:tbasis%npts,1:pub_num_spins), &
            tatom%guess_comp_den(1:tbasis%npts,1:pub_num_spins), &
            tatom%guess_rhoij(1:tatom%nshells,1:tatom%nshells,1:pub_num_spins))
    end if

    ! Write the orbitals out, if required
    if (pub_write_initial_radial_ngwfs) then
       write(tmp,'(i3.3)') current_element%species_number
       if (any(tatom%splitnorm(:,:)>0.0_DP).or.(any(tatom%polarise(:)>0))) then
          filename = trim(pub_rootname)//'_atom_orbs_'//tmp
          call atom_write_orbitals(tatom,tbasis,filename)
       end if
       filename = trim(pub_rootname)//'_initial_rad_ngwf_'//tmp
       call ngwfs_radial_write(radial_ngwf,filename)
    end if

    ! Clean up the atom and basis
    call atom_destroy(tatom)
    call atom_destroy_basis(tbasis)

  end subroutine ngwfs_solve_atom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_radial_write(radial_ngwf,filename)

    !=====================================================================!
    ! This subroutine allocates storage for a RADIAL_NGWF_TYPE structure. !
    !---------------------------------------------------------------------!
    ! Written by Nicholas Hine on 03/09/2010.                             !
    !=====================================================================!

    use utils, only: utils_open_unit_check, utils_close_unit_check, utils_unit

    implicit none

    ! Arguments
    type(RADIAL_NGWF_TYPE), intent(inout) :: radial_ngwf
    character(len=80), intent(in) :: filename

    ! Local Variables
    integer :: ierr
    integer :: iunit
    integer :: ipt
    integer :: ishell

    ! ndmh: open the file
    iunit = utils_unit()
    open(iunit,file=trim(filename),iostat=ierr)
    call utils_open_unit_check('ngwfs_radial_write',trim(filename),ierr)

    ! ndmh: write the orbitals
    do ipt=1,radial_ngwf%npts
       write(iunit,'(f22.15)',advance='no') radial_ngwf%rad(ipt)
       do ishell=1,radial_ngwf%nshells
          write(iunit,'(f22.15)',advance='no') radial_ngwf%func_real(ipt,ishell)
       end do
       write(iunit,*)
   end do

    ! ndmh: close the file
    close(iunit,iostat=ierr)
    call utils_close_unit_check('ngwfs_radial_write',trim(filename),ierr)

  end subroutine ngwfs_radial_write


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_radial_create(radial_ngwf)

    !=====================================================================!
    ! This subroutine allocates storage for a RADIAL_NGWF_TYPE structure. !
    !---------------------------------------------------------------------!
    ! Written by Nicholas Hine on 03/09/2010.                             !
    !=====================================================================!

    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(RADIAL_NGWF_TYPE), intent(inout) :: radial_ngwf

    ! Local Variables
    integer :: ierr

    ! ndmh: Default maximum sizes (not necessarily all used)
    radial_ngwf%nshells = 100
    radial_ngwf%npts = 4001

    allocate(radial_ngwf%angmom(radial_ngwf%nshells),stat=ierr)
    call utils_alloc_check('ngwfs_radial_create','radial_ngwf%angmom',ierr)
    allocate(radial_ngwf%rc(radial_ngwf%nshells),stat=ierr)
    call utils_alloc_check('ngwfs_radial_create','radial_ngwf%rc',ierr)
    allocate(radial_ngwf%rad(radial_ngwf%npts),stat=ierr)
    call utils_alloc_check('ngwfs_radial_create','radial_ngwf%rad',ierr)
    allocate(radial_ngwf%func_real(radial_ngwf%npts,radial_ngwf%nshells),stat=ierr)
    call utils_alloc_check('ngwfs_radial_create','radial_ngwf%func_real',ierr)

    ! ndmh: Default maximum size for a fireball
    radial_ngwf%npts = 2001
    radial_ngwf%angmom = 0
    radial_ngwf%rc = 0.0_DP
    radial_ngwf%rad = 0.0_DP
    radial_ngwf%func_real = 0.0_DP

  end subroutine ngwfs_radial_create


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_radial_destroy(radial_ngwf)

    !=======================================================================!
    ! This subroutine deallocates storage for a RADIAL_NGWF_TYPE structure. !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine on 03/09/2010.                               !
    !=======================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(RADIAL_NGWF_TYPE), intent(inout) :: radial_ngwf

    ! Local Variables
    integer :: ierr

    if (.not.allocated(radial_ngwf%func_real)) return

    deallocate(radial_ngwf%func_real,stat=ierr)
    call utils_dealloc_check('ngwfs_radial_destroy','radial_ngwf%func_real',ierr)
    deallocate(radial_ngwf%rad,stat=ierr)
    call utils_dealloc_check('ngwfs_radial_destroy','radial_ngwf%rad',ierr)
    deallocate(radial_ngwf%rc,stat=ierr)
    call utils_dealloc_check('ngwfs_radial_destroy','radial_ngwf%rc',ierr)
    deallocate(radial_ngwf%angmom,stat=ierr)
    call utils_dealloc_check('ngwfs_radial_destroy','radial_ngwf%angmom',ierr)

  end subroutine ngwfs_radial_destroy



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_sto3g_create

    !===================================================================!
    ! Allocate appropriate memory and initialise it to hold Gaussian    !
    ! basis set parameters (all STO-3G and only the polarisation        !
    ! functions of 6-31G*) for each element up to atomic number         !
    ! pub_num_gtatoms
    !-------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 31/03/2006.                   !
    !===================================================================!

    use constants, only: DP
    use ngwf_data, only: pub_num_gtatoms, &
       ngwf_data_cocos_and_expos01_20, &
       ngwf_data_cocos_and_expos21_40, &
       ngwf_data_cocos_and_expos41_60, &
       ngwf_data_cocos_and_expos61_80, &
       ngwf_data_cocos_and_expos81_103
    use utils, only: utils_alloc_check
    implicit none

    ! cks: <<local variables>>
    integer :: atom     ! atom counter
    integer :: n_shells ! number of shells
    integer :: n_prim   ! number of primitives in given shell
    integer :: ierr     ! allocation error flag

    if (allocated(gbasis)) return

    ! cks: allocate gbasis array for all supported elements
    allocate(gbasis(pub_num_gtatoms), stat=ierr)
    call utils_alloc_check('ngwfs_sto3g_create','gbasis',ierr)

    ! cks: Number of shells per element, including one polarisation shell
    ! cks: This data was typed by Shyong Chen
    gbasis(1)%nshells=2
    gbasis(2)%nshells=2
    gbasis(3)%nshells=4
    gbasis(4)%nshells=4
    gbasis(5)%nshells=4
    gbasis(6)%nshells=4
    gbasis(7)%nshells=4
    gbasis(8)%nshells=4
    gbasis(9)%nshells=4
    gbasis(10)%nshells=4
    gbasis(11)%nshells=6
    gbasis(12)%nshells=6
    gbasis(13)%nshells=6
    gbasis(14)%nshells=6
    gbasis(15)%nshells=6
    gbasis(16)%nshells=6
    gbasis(17)%nshells=6
    gbasis(18)%nshells=6
    gbasis(19)%nshells=8
    gbasis(20)%nshells=8
    gbasis(21)%nshells=9
    gbasis(22)%nshells=9
    gbasis(23)%nshells=9
    gbasis(24)%nshells=9
    gbasis(25)%nshells=9
    gbasis(26)%nshells=9
    gbasis(27)%nshells=9
    gbasis(28)%nshells=9
    gbasis(29)%nshells=9
    gbasis(30)%nshells=9
    gbasis(31)%nshells=9
    gbasis(32)%nshells=9
    gbasis(33)%nshells=9
    gbasis(34)%nshells=9
    gbasis(35)%nshells=9
    gbasis(36)%nshells=9
    gbasis(37)%nshells=11
    gbasis(38)%nshells=11
    gbasis(39)%nshells=12
    gbasis(40)%nshells=12
    gbasis(41)%nshells=12
    gbasis(42)%nshells=12
    gbasis(43)%nshells=12
    gbasis(44)%nshells=12
    gbasis(45)%nshells=12
    gbasis(46)%nshells=12
    gbasis(47)%nshells=12
    gbasis(48)%nshells=12
    gbasis(49)%nshells=12
    gbasis(50)%nshells=12
    gbasis(51)%nshells=12
    gbasis(52)%nshells=12
    gbasis(53)%nshells=12
    gbasis(54)%nshells=13
    gbasis(55)%nshells=14
    gbasis(56)%nshells=14
    gbasis(57)%nshells=15
    gbasis(58)%nshells=15
    gbasis(59)%nshells=15
    gbasis(60)%nshells=15
    gbasis(61)%nshells=15
    gbasis(62)%nshells=15
    gbasis(63)%nshells=15
    gbasis(64)%nshells=15
    gbasis(65)%nshells=15
    gbasis(66)%nshells=15
    gbasis(67)%nshells=15
    gbasis(68)%nshells=15
    gbasis(69)%nshells=15
    gbasis(70)%nshells=15
    gbasis(71)%nshells=15
    gbasis(72)%nshells=15
    gbasis(73)%nshells=15
    gbasis(74)%nshells=15
    gbasis(75)%nshells=15
    gbasis(76)%nshells=15
    gbasis(77)%nshells=15
    gbasis(78)%nshells=15
    gbasis(79)%nshells=15
    gbasis(80)%nshells=15
    gbasis(81)%nshells=15
    gbasis(82)%nshells=15
    gbasis(83)%nshells=15
    gbasis(84)%nshells=15
    gbasis(85)%nshells=15
    gbasis(86)%nshells=16
    gbasis(87)%nshells=16
    gbasis(88)%nshells=16
    gbasis(89)%nshells=16
    gbasis(90)%nshells=17
    gbasis(91)%nshells=18
    gbasis(92)%nshells=18
    gbasis(93)%nshells=18
    gbasis(94)%nshells=18
    gbasis(95)%nshells=17
    gbasis(96)%nshells=17
    gbasis(97)%nshells=17
    gbasis(98)%nshells=17
    gbasis(99)%nshells=17
    gbasis(100)%nshells=17
    gbasis(101)%nshells=17
    gbasis(102)%nshells=17
    gbasis(103)%nshells=17


    ! cks: allocate memory for the components of each shell
    do atom=1, pub_num_gtatoms

       n_shells = gbasis(atom)%nshells
       n_prim =3

       ! cks: allocate angular momentum memory
       allocate(gbasis(atom)%angmom(n_shells), stat=ierr)
       call utils_alloc_check('ngwfs_sto3g_create','=>gbasis(atom)%angmom',ierr)
       ! cks: initialise
       gbasis(atom)%angmom =-1


       ! cks: allocate exponents memory
       allocate(gbasis(atom)%expo(n_shells, n_prim), stat=ierr)
       call utils_alloc_check('ngwfs_sto3g_create','=>gbasis(atom)%expo',ierr)
       ! cks: initialise
       gbasis(atom)%expo =0.0_DP

       ! cks: allocate contraction coefficients memory
       allocate(gbasis(atom)%coco(n_shells, n_prim), stat=ierr)
       call utils_alloc_check('ngwfs_sto3g_create','=>gbasis(atom)%coco',ierr)
       ! cks: initialise
       gbasis(atom)%coco =0.0_DP

    enddo



    ! cks: load cocos and expos in the allocated memory
    call ngwf_data_cocos_and_expos01_20(gbasis(:))
    call ngwf_data_cocos_and_expos21_40(gbasis(:))
    call ngwf_data_cocos_and_expos41_60(gbasis(:))
    call ngwf_data_cocos_and_expos61_80(gbasis(:))
    call ngwf_data_cocos_and_expos81_103(gbasis(:))


  end subroutine ngwfs_sto3g_create


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_sto3g_destroy

    !=======================================================!
    ! Deallocate all memory of gbasis array.                !
    !---------------------------------------- --------------!
    ! Written by Chris-Kriton Skylaris on 31/03/2006.       !
    !=======================================================!

    use ngwf_data, only: pub_num_gtatoms
    use utils, only: utils_dealloc_check
    implicit none

    ! cks: <<local variables>>
    integer :: atom     ! atom counter
    integer :: ierr     ! allocation error flag

    if (.not.allocated(gbasis)) return

    ! cks: deallocate memory for the components of each shell
    do atom=1, pub_num_gtatoms

       ! cks: deallocate angular momentum memory
       deallocate(gbasis(atom)%angmom, stat=ierr)
       call utils_dealloc_check('ngwfs_sto3g_destroy','=>gbasis(atom)%angmom',ierr)

       ! cks: deallocate exponents memory
       deallocate(gbasis(atom)%expo, stat=ierr)
       call utils_dealloc_check('ngwfs_sto3g_destroy','=>gbasis(atom)%expo',ierr)

       ! cks: deallocate contraction coefficients memory
       deallocate(gbasis(atom)%coco, stat=ierr)
       call utils_dealloc_check('ngwfs_sto3g_destroy','=>gbasis(atom)%coco',ierr)

    enddo

    ! cks: deallocate the gbasis array
    deallocate(gbasis, stat=ierr)
    call utils_dealloc_check('ngwfs_sto3g_destroy','gbasis',ierr)


  end subroutine ngwfs_sto3g_destroy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module ngwfs_radial



