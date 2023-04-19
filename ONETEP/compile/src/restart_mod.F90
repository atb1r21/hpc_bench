! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!============================================================================!
!                                                                            !
!               Electronic optimisation restart module                       !
!                                                                            !
!   The subroutines in this file were written by Chris-Kriton Skylaris       !
!   They are mainly subroutines for storing/retrieving from disk information !
!   in binary format that is used for continuining/restarting a ONETEP       !
!   calculation.                                                             !
!----------------------------------------------------------------------------!
!   Written by Chris-Kriton Skylaris, February 2004                          !
!                                                                            !
!   TCM Group, Cavendish laboratory                                          !
!   Madingley Road                                                           !
!   Cambridge CB3 0HE                                                        !
!   UK                                                                       !
!----------------------------------------------------------------------------!
!   Subsequent modifications by Nicholas Hine, Peter D. Haynes,              !
!   Alvaro Ruiz Serrano, Tim Zuehlsdorff, David O'Regan,                     !
!   Gabriel Constantinescu, Andrea Greco and Jose M. Escartin.               !
!============================================================================!

#if defined(GPU_PGI) || defined(PGI_BUG)
#define NOMPIIO
#endif

module restart

  use constants, only: dp
#ifdef HDF5
  use hdf5
#endif
#if defined(MPI) && defined(HDF5) && !defined(NOMPIIO) && !defined(USE_INCLUDE_MPIF)
  use mpi !! External dependency
#endif
  use rundat, only: pub_debug_on_root

  implicit none
#if defined(MPI) && defined(HDF5) && !defined(NOMPIIO) && defined(USE_INCLUDE_MPIF)
#include "mpif.h"
#endif

  private

  !=== Subroutines ==========================================================!

  ! ndmh: routines to read/write NGWFs
  public :: restart_ngwfs_tightbox_input
  public :: restart_ngwfs_tightbox_output
  public :: restart_sph_waves_input
  public :: restart_sph_waves_output

  ! ndmh: routines to read/write density kernels
  public :: restart_kernel_write
  public :: restart_kernel_read

  ! ndmh: routines to read/write hamiltonian and overlap matrices
  public :: restart_hamiltonian_write
  public :: restart_hamiltonian_read
  public :: restart_overlap_write
  public :: restart_overlap_read

  ! tjz07: TDDFT restart subroutines
  public :: restart_response_kernel_batch_write
  public :: restart_response_kernel_batch_read
  public :: restart_RPA_kernel_read
  public :: restart_RPA_kernel_write

  ! jd: This is set by auto-solvation to temporarily change input or output
  !     filenames for restart. Not pretty, but the alternative is to pass flags
  !     through a forest of calls.
  character(len=7), public, save :: restart_read_filename_prefix = ''
  character(len=7), public, save :: restart_write_filename_prefix = ''

  ! jme: Tightbox file format versions.
  integer, parameter :: OLD_TIGHTBOX_FORMAT  = 100
  integer, parameter :: TIGHTBOX_FORMAT_2017 = 200
  integer, parameter :: CURRENT_TIGHTBOX_FORMAT = TIGHTBOX_FORMAT_2017

  ! jme: Integer values used to store true/false in our tightbox HDF files.
#ifdef HDF5
  integer, parameter :: HDF5_FALSE_INT = 151
  integer, parameter :: HDF5_TRUE_INT  = 446
#endif

contains

!============================================================================!
!============================================================================!

  subroutine restart_sph_waves_output(funcs_on_grid, fbasis, fftbox, uni_tightbox, &
       cell, elements, file_extension, regions)

    !=======================================================================!
    ! Write functions on grid in spherical waves representation on a file so!
    ! calculations can be restarted.                                        !
    !-----------------------------------------------------------------------!
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Originally written by Alvaro Ruiz Serrano in January 2009.            !
    ! Revisited by Alvaro Ruiz Serrano in November 2011 for re-structuring  !
    ! and bug-fixing.                                                       !
    ! Modified to remove pub_par by Joseph Prentice, May 2018               !
    !=======================================================================!

    use datatypes, only: FUNCTIONS
    use constants, only: DP, stdout,NORMAL
    use comms, only: pub_on_root,pub_root_proc_id, &
         comms_reduce, comms_barrier, comms_free
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS, function_basis_ppds_to_sph_waves
    use ion, only : ELEMENT
    use model_type, only: REGION
    use rundat, only : pub_rootname, pub_output_detail, &
         pub_write_max_l,pub_extra_n_sw
    use simulation_cell, only: CELL_INFO
    use spherical_wave, only: sw_init,sw_bessel_zeros, sw_bessel_zeros_init
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_assert, utils_abort

    implicit none

    ! ars: << arguments >>
    type(REGION), intent(in) :: regions
    type(FUNC_BASIS), intent(in) :: fbasis
    type(FUNCTIONS),  intent(in) :: funcs_on_grid
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FFTBOX_INFO), intent(in) :: uni_tightbox
    type(ELEMENT),    intent(in) :: elements(regions%par%nat)
    character(len=*), intent(in) :: file_extension ! ddor: Extension of filename to write

    ! ars: << local variables>>
    character(len=256) :: output_file ! output file name buffer
    integer :: output_unit            ! fortran output unit number
    integer :: distr_func, atom_func  ! ngwf counters
    integer :: orig_atom              ! global atom counter in input file order
    integer :: ll, mm, nn             ! angular momentum counters
    integer :: ierr                   ! memory allocation error flag
    integer :: maxn                   ! ars: maximum number of Bessel zeros into the biggest sphere on each proc
    integer :: num_swcoeff_on_file    ! ars; Number of SW coeffcients (size of the SW basis set)
    real(kind=DP), allocatable :: sw_coeffs (:,:,:) ! ars: spherical waves coefficients
    integer, allocatable :: maxln(:,:) ! ars: maximum n for a given atom and l number
    real(kind=DP) :: alpha_init,alpha,radius ! ars: calculation of maxn
    ! agrecocmplx
    logical :: loc_cmplx


    call timer_clock("restart_sph_waves_output", 1)

    loc_cmplx = funcs_on_grid%iscmplx

    ! agrecocmplx
    call utils_assert(.not. loc_cmplx, 'Error in&
         & restart_sph_waves_output: not ready yet for complex NGWFs.')

    !---------------------------------------------------------------------------
    ! ars: header
    !---------------------------------------------------------------------------
    if (pub_on_root) then
       write(output_file,'(a80)') trim(pub_rootname)//"."//trim(file_extension)
       output_file = adjustl(output_file)

       ! ddor: Writeout depends on the filename extension.
       if (trim(file_extension) .eq. 'sw_ngwfs') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Writing NGWFs to file "', trim(output_file),'"....'
       elseif  (trim(file_extension) .eq. 'sw_hub_projs') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Writing HUBBARD PROJs to file "', trim(output_file),'"....'
       else if (trim(file_extension) .eq. 'sw_ngwfs_cond') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Writing Conduction NGWFs to file "', trim(output_file),'"....'
       else
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Writing functions to file "', trim(output_file),'"....'
       endif
    end if

    !---------------------------------------------------------------------------
    ! ars: initialisation of SW structures
    !---------------------------------------------------------------------------


    ! ars : Initialise sw module using approximate method.
    maxn = ceiling(maxval(fbasis%spheres(:)%radius) / &
         max(cell%d1,cell%d2,cell%d3))
    call comms_reduce('MAX',maxn)
    call sw_bessel_zeros_init(2*maxn,pub_write_max_l)

    allocate(maxln(0:pub_write_max_l,regions%par%nat), stat=ierr)
    call utils_alloc_check('restart_sph_waves_output','maxn',ierr)

    ! ars: calculate maxn for each l
    maxln = 1
    do orig_atom=1,regions%par%nat
       alpha_init = max(cell%d1,cell%d2,cell%d3)/&
                        elements(orig_atom)%radius

       do ll=0, pub_write_max_l
          nn = 1
          alpha = alpha_init

          do while(alpha.ge.alpha_init)
             nn = nn + 1
             ! ars : Security statement that re-initialises sw_bessel_zeros
             ! ars : if maxn_proc increases in more than 2 times its
             ! ars :  previous value
             if (nn.gt.maxn+5) then
                maxn=nn
                call sw_bessel_zeros_init(2*maxn,pub_write_max_l)
             endif
             alpha = 1 - sw_bessel_zeros(nn-1,ll)/sw_bessel_zeros(nn,ll)

          enddo
          maxln(ll,orig_atom)=(nn-1)+pub_extra_n_sw
       enddo
    enddo

    ! ars : Maximum value of n within all procs.
    maxn=maxval(maxln)

    ! ars: init SW module
    radius = maxval(fbasis%spheres(:)%radius)**2 * 2 &
         * maxval(uni_tightbox%recip_grid(5,:,:,:))
    call comms_reduce('MAX',radius)
    call sw_init(pub_write_max_l,maxn,radius)

    ! ars : calculate the total number of SW coefficients
    num_swcoeff_on_file = 0
    do orig_atom =1, regions%par%nat
       do atom_func=1,fbasis%num_on_atom(regions%par%distr_atom(orig_atom))
          do ll =0, pub_write_max_l
             num_swcoeff_on_file = num_swcoeff_on_file + (2*ll+1)*maxln(ll,orig_atom)
          end do
       enddo
    end do

    ! ars : Check if any value of max_nl_matrix is less than or equal to zero as
    ! ars : pub_extra_n_sw can be negative
    call utils_assert(.not. any(maxln.le.0), 'Negative n number for spherical &
         &waves. The values for pub_extra_n_sw and minimum n number were ', &
         pub_extra_n_sw, minval(maxln))

    !---------------------------------------------------------------------------
    ! ars: end initialisation
    !---------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    ! ars: write on disk
    !---------------------------------------------------------------------------

    output_unit =utils_unit()
    if (pub_on_root) then

       open(unit=output_unit, form="unformatted" ,file=trim(output_file), &
            action="write" )

       write(output_unit,err=140) fbasis%num          ! Number of NGWFs
       write(output_unit,err=140) num_swcoeff_on_file ! Number of SW coeffs
       write(output_unit,err=140) pub_write_max_l     ! Max_l
       write(output_unit,err=140) maxn                ! Max_n (global)
       do orig_atom=1,regions%par%nat                           ! Element, Atom,l,max_n
          do ll=0, pub_write_max_l
             write(output_unit, err=140) &
                elements(orig_atom)%species_id,elements(orig_atom)%symbol,&
                orig_atom,ll,maxln(ll,orig_atom)
          enddo
       enddo
    endif

    ! ars : allocate sw coeffficient
    allocate(sw_coeffs(0:pub_write_max_l,-pub_write_max_l:pub_write_max_l,1:maxn),stat=ierr)
    call utils_alloc_check('restart_sph_waves_output','sw_coeffs',ierr)

    ! cks: Loop over all atoms in the order they appear in the input file
    do orig_atom =1, regions%par%nat

       ! ars: loop over functions on atom
       do atom_func= 1,fbasis%num_on_atom(regions%par%distr_atom(orig_atom))

          distr_func = fbasis%first_on_atom(regions%par%distr_atom(orig_atom)) + &
               atom_func - 1

          ! cks: this essentially converts the non-blocking sends of
          ! cks: comms_send to blocking sends. Since parallel performance
          ! cks: is not the issue here, this makes the code simpler.
          call comms_barrier

          call function_basis_ppds_to_sph_waves(sw_coeffs, &
               maxln(:,orig_atom), maxn, distr_func, pub_root_proc_id, &
               funcs_on_grid, fbasis, fftbox, uni_tightbox, cell)


          ! &&&&&&&&&&& WRITE THE BUFFER &&&&&&&&&&&&&&&&&&&&&&&
          if (pub_on_root) then
             do ll = 0, pub_write_max_l
                do mm = -ll,ll
                   do nn = 1, maxln(ll,orig_atom)

                      ! Atom,NGWF,l,m,n,swcoeff
                      write(output_unit,err=140) &
                           elements(orig_atom)%species_id,&
                           elements(orig_atom)%symbol,&
                           orig_atom, atom_func,ll,mm,nn,sw_coeffs(ll,mm,nn)

                   end do
                end do
             end do
          end if
          ! &&&&&&& END WRITE THE BUFFER &&&&&&&&&&&&&&&&&&&&&&&

       end do
    end do

    ! cks: root proc closes output file
    if (pub_on_root) then
       close(unit =output_unit)

       if (pub_output_detail >= NORMAL ) write(stdout,*)' done'

    endif


    ! cks: don't go away just yet! Wait for all communications to finish!
    ! cks: also wait before deallocating uni_tbox which is used in non-blocking
    ! cks: send!
    call comms_barrier
    ! cks: clean up send stack too
    call comms_free

    ! ars: deallocate and finish
    deallocate(sw_coeffs, stat=ierr)
    call utils_dealloc_check('restart_sph_waves_output','sw_coeffs',ierr)
    deallocate(maxln, stat=ierr)
    call utils_dealloc_check('restart_sph_waves_output','max_n',ierr)


    call timer_clock("restart_sph_waves_output", 2)

    return

140 call utils_abort('Problem writing to file in&
                      & restart_sph_waves_output. ONETEP stops')

  end subroutine restart_sph_waves_output


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine restart_sph_waves_input(funcs_on_grid,fbasis,cell,fftbox, &
       uni_tightbox,file_extension,regions)


    !=======================================================================!
    ! Read functions on grid in spherical waves representation on a file so !
    ! calculations can be restarted.                                        !
    !-----------------------------------------------------------------------!
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Originally written by Alvaro Ruiz Serrano in January 2009.            !
    ! Revisited by Alvaro Ruiz Serrano in November 2011 for re-structuring  !
    ! and bug-fixing.                                                       !
    ! Modified to remove pub_par by Joseph Prentice, May 2018               !
    !=======================================================================!

    use datatypes, only: FUNCTIONS, data_set_to_zero
    use constants, only: DP, stdout, NORMAL
    use comms, only: comms_bcast, comms_barrier, comms_free, &
         pub_root_proc_id, pub_on_root
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS,function_basis_sph_waves_to_ppds
    use model_type, only: REGION
    use rundat, only : pub_rootname, pub_read_max_l, pub_output_detail
    use simulation_cell, only: CELL_INFO
    use spherical_wave, only: sw_init
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_assert, utils_abort


    implicit none

    ! ars: Arguments
    type(REGION), intent(in) :: regions
    type(FUNC_BASIS), intent(in   ) :: fbasis
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FFTBOX_INFO), intent(in) :: uni_tightbox
    type(FUNCTIONS),    intent(inout) :: funcs_on_grid
    character(len=*), intent(in   ) :: file_extension ! ddor: Extension of filename to write

    ! ars: Variables
    character(len=256) :: read_file  ! read file name buffer
    integer :: input_unit   ! fortran input unit number
    integer :: ierr         ! memory allocation error flag
    integer :: orig_atom    ! global atom counter in input file order
    integer :: atom_ngwf,distr_ngwf ! ngwf of current atom counter
    integer :: ll,mm,nn,ii
    integer :: num_funcs_on_file    ! total number of ngwfs read from file
    integer :: maxl                 ! maximum value of l within all funcs
    integer :: maxn                 ! maximum value of n within all funcs
    character(len=2) :: read_atom_symbol
    character(len=4) :: read_species_id
    integer :: read_atom,read_ngwf,read_ll,read_mm,read_nn

    integer :: num_swcoeffs, num_swcoeffs_perfunc,num_swcoeffs_on_file
    integer :: skip_lines ! skip a number of lines if required
    real(kind=DP) :: max_radius
    integer, allocatable :: maxln(:,:)    ! ars: n zeros array in serial
    real(kind=DP), allocatable :: sw_coeffs(:,:,:) ! ars: spherical waves coefficients


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Initialisations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    call timer_clock("restart_sph_waves_input", 1)

    ! jme
    call utils_assert(.not. funcs_on_grid%iscmplx, &
         'Error in restart_sph_waves_input: this subroutine is not ready yet &
         &for complex NGWFs.')

    !---------------------------------------------------------------------------
    ! ars: header
    !---------------------------------------------------------------------------
    if (pub_on_root) then
       write(read_file,'(a80)') trim(pub_rootname)//"."//trim(file_extension)
       read_file = adjustl(trim(read_file))

       ! ddor: Writeout depends on the filename extension.
       if (trim(file_extension) .eq. 'sw_ngwfs') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Reading NGWFs from file "', trim(read_file),'"....'
       elseif  (trim(file_extension) .eq. 'sw_hub_projs') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Reading HUBBARD PROJs from file "', trim(read_file),'"....'
       else if (trim(file_extension) .eq. 'sw_ngwfs_cond') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Reading Conduction NGWFs from file "', trim(read_file),'"....'
       else
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Reading functions from file "', trim(read_file),'"....'
       endif
    end if

    ! ars : radius to initialise sw module.
    ! ars : this must be the same than the used in restart_ngwfs_swcoeff_output

    max_radius = maxval(fbasis%spheres(:)%radius)**2 * 2 &
         * maxval(uni_tightbox%recip_grid(5,:,:,:))

    ! ars : initialise checking variables
    num_swcoeffs = 0
    num_swcoeffs_perfunc = 0

    ! ars : this is the array that comes out from this subroutine
    call data_set_to_zero(funcs_on_grid)

    !---------------------------------------------------------------------------
    ! ars: read headers from disk and safe-check
    !---------------------------------------------------------------------------

    ! cks: get a unit number that is free
    input_unit =utils_unit()

    ! cks: root proc opens input file, reads basic info
    if (pub_on_root) then

       open(unit =input_unit, form="unformatted" ,file=trim(read_file), &
            action="read", status='old')


       ! ars : read information on file and broadcast
       read(input_unit,err =150) num_funcs_on_file    ! Number of funcs saved on file
       read(input_unit,err =150) num_swcoeffs_on_file ! Number of SW coeffs
       read(input_unit,err =150) maxl                 ! Max_l
       read(input_unit,err =150) maxn                 ! Max_n

    endif
    call comms_bcast(pub_root_proc_id, num_funcs_on_file)
    call comms_bcast(pub_root_proc_id, num_swcoeffs_on_file)
    call comms_bcast(pub_root_proc_id, maxl)
    call comms_bcast(pub_root_proc_id, maxn)


    ! allocate max_n, read and broadcast
    allocate(maxln(0:maxl,regions%par%nat),stat=ierr)
    call utils_alloc_check('restart_sph_waves_input','maxln',ierr)
    maxln = 0
    if (pub_on_root) then
       ! ars: read from file and fill max_nl_matrix_serial
       do orig_atom=1,regions%par%nat
          do ll=0, maxl
             read(input_unit,err=150) read_species_id, &
                read_atom_symbol,read_atom,read_ll,maxln(read_ll,read_atom)
          enddo
       enddo
    endif
    call comms_bcast(pub_root_proc_id, maxln)

    ! ars: check if the information has been read correctly
    call utils_assert(num_funcs_on_file == fbasis%num, &
         'The number of basis set functions is different from the number of &
         &functions in the file. The values in question were: ', &
         fbasis%num, num_funcs_on_file)

    if (pub_read_max_l.gt.maxl) then
       if (pub_on_root) then
          write(stdout,*)'max l on file = ', maxl
          write(stdout,*)'read_max_l = ', pub_read_max_l
          write(stdout,*)'The angular momentum in the input file is higher &
               &than in the file'
          write(stdout,*)'Maximum angular momentum has been reset to: ',&
               maxl
       endif
       pub_read_max_l = maxl
    endif

    ! ars : initialise spherical_wave module
    call sw_init(maxl,maxn+1,max_radius)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End initialisations ~~~~~~~~~~~~~~~~~~~~~~~~~~!



    !---------------------------------------------------------------------------
    ! ars: init SW module and read coeffs from disk
    !---------------------------------------------------------------------------

    ! ars : allocate arrays
    allocate(sw_coeffs(0:pub_read_max_l,-pub_read_max_l:pub_read_max_l,maxn),stat=ierr)
    call utils_alloc_check ('restart_sph_waves_input','sw_coeffs',ierr)
    sw_coeffs = 0.0_DP

    ! cks: Loop over all atoms in the order they appear in the input file
    do orig_atom =1, regions%par%nat

       ! ars: loop over functions on atom
       do atom_ngwf =1,fbasis%num_on_atom(regions%par%distr_atom(orig_atom))

          distr_ngwf = fbasis%first_on_atom(regions%par%distr_atom(orig_atom)) + &
               atom_ngwf - 1


          ! cks: this essentially converts the non-blocking sends of
          ! cks: comms_send to blocking sends. Since parallel performance
          ! cks: is not the issue here, this makes the code simpler.
          call comms_barrier

          ! ars : n_swcoeffs_ngwf becomes zero after the ngwf loop
          num_swcoeffs_perfunc=0
          ! ars : if pub_read_max_l < read_max_l then the program "skip_liness" extra lines
          skip_lines=0
          do ll = 0,maxl
             skip_lines=skip_lines+maxln(ll,orig_atom)*(2*ll+1)
          enddo

          ! ars : loops over l, m, n
          do ll = 0,pub_read_max_l
             do mm = -ll,ll
                do nn = 1, maxln(ll,orig_atom)


                   ! &&&&&&&&&&& READ INTO BUFFER &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

                   if (pub_on_root) then
                      read(input_unit, err=150)&
                           read_species_id,read_atom_symbol,read_atom,read_ngwf,&
                           read_ll,read_mm,read_nn,sw_coeffs(read_ll,read_mm,read_nn)

                      num_swcoeffs = num_swcoeffs + 1
                      num_swcoeffs_perfunc = num_swcoeffs_perfunc + 1
                   endif
                   ! &&&&&&& END READ INTO BUFFER &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


                enddo  !!! End Loop 5 !!!
             enddo  !!! End Loop 4 !!!
          enddo  !!! End Loop 3 !!!

          if (pub_on_root) then
             do ii=1,skip_lines-num_swcoeffs_perfunc
                read(input_unit, err=150) ! ars: read nothing
             enddo
          end if

          ! ars: create functions in PPD representation
          call function_basis_sph_waves_to_ppds(fbasis,funcs_on_grid, &
               fftbox, uni_tightbox, cell, sw_coeffs, &
               maxln(:,orig_atom),pub_read_max_l,maxn,distr_ngwf,pub_root_proc_id)

       enddo  !!! End Loop 2 !!!
    enddo  !!! End Loop 1 !!!


    ! cks: root proc closes output file
    if (pub_on_root) then
       close(unit=input_unit)

       ! ars : print info about the reading process
       if (pub_output_detail.gt.NORMAL) then
          write(stdout,'(a/)')' done'
          write(stdout,'(a,13x,i0)') "Maximum l of the SW read = ", pub_read_max_l
          write(stdout,'(a,17x,i0)') "Number of SW coefficients read = ",&
               num_swcoeffs
          write(stdout,'(a,i0)') "Number of SW coefficients in the file = ",&
               num_swcoeffs_on_file
       end if
    endif

    call utils_assert(num_swcoeffs <= num_swcoeffs_on_file, &
       'The number of coefficients read is bigger than&
            & the number of coefficients in the file.')


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Deallocate arrays ~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    deallocate(maxln, stat=ierr)
    call utils_dealloc_check('restart_sph_waves_input','maxln',ierr)
    deallocate(sw_coeffs,stat=ierr)
    call utils_dealloc_check('restart_sph_waves_input','sw_coeffs',ierr)
    !~~~~~~~~~~~~~~~~~~~~~~~~~ End deallocate arrays ~~~~~~~~~~~~~~~~~~~~~~~~~~~!


    ! cks: close comms
    call comms_barrier
    call comms_free


    call timer_clock("restart_sph_waves_input", 2)

    return

150 call utils_abort('Problem reading from file in restart_sph_waves_input.')

  end subroutine restart_sph_waves_input


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine restart_ngwfs_tightbox_output(ngwfs_on_grid, ngwf_basis, &
       cell, fftbox, elements, file_extension, regions, label, ngwfs_type)

    !=======================================================================!
    ! This subroutine outputs the current NGWFs in "universal tightbox"     !
    ! representation to a .tightbox_ngwfs file. The order of the NGWFs is   !
    ! the one in which they appear in the input file and is not affected    !
    ! by the use of the space-filling curve. Also the .tightbox_ngwfs       !
    ! remains usable if the positions of the atoms change slightly          !
    ! after the restart.                                                    !
    ! ndmh 18/12/2009: In cases where the FFTbox coincides with the cell    !
    ! but the individual tightboxes do not need to, a new temporary NGWF    !
    ! basis is constructed, for which the FFTbox only coincides with the    !
    ! cell along direction i if L_i < 2*R_NGWF + delta_i.                   !
    !-----------------------------------------------------------------------!
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 12/3/2004.                        !
    ! Modified by Chris-Kriton Skylaris on 13/4/2004 so that it also outputs!
    ! the coordinates of the centre of each NGWF in (real) numbers of grid  !
    ! points from the origin of its tightbox.                               !
    ! Modified by D. O'Regan in November 2009 for different file extensions.!
    ! Modified by Nicholas Hine in December 2009, to allow reduced file     !
    ! sizes by only reading in tightboxes even if FFTbox coincides with     !
    ! cell in some directions.                                              !
    ! Modified by Gabriel Constantinescu in June 2014 to use MPI-IO to write!
    ! the tightboxes, instead of sending them from the root proc            !
    ! Modified by Andrea Greco in May 2015 to allow use of complex NGWFs.   !
    ! Modified by Jose M Escartin in January 2017 to incorporate new format !
    ! for tightbox files.                                                   !
    ! Const-correctness for 'fftbox' fixed by Jacek Dziedzic in June 2017.  !
    ! Modified to remove pub_par by Joseph Prentice, May 2018               !
    !=======================================================================!

    use comms, only: comms_reduce, comms_barrier, comms_free, &
#if defined(MPI) && !defined(NOMPIIO)
         comms_open_file, comms_close_file, comms_write, &
#endif
         comms_bcast,pub_on_root, pub_root_proc_id, pub_my_proc_id
    use constants, only: DP, stdout, NORMAL, LONG
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS, function_basis_allocate, &
         function_basis_deallocate, function_basis_distribute, &
         function_basis_copy_spheres, function_basis_init_tight_boxes
    use geometry, only: magnitude
    use ion, only : ELEMENT
    use model_type, only: REGION
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_rootname, pub_output_detail, pub_debug_on_root, &
         pub_read_real_ngwfs, pub_devel_code
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only : utils_unit, &
         utils_abort, utils_int_to_str, utils_devel_code

    implicit none

    ! Arguments
    type(REGION), target, intent(in) :: regions
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    !real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%size_on_grid) ! NGWFs on this proc
    ! agrecocmplx: use FUNCTIONS type instead
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid ! NGWFs on this proc
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(ELEMENT), intent(in) :: elements(regions%par%nat) ! elements of all proc (in input file order)
    character(len=*), intent(in) :: file_extension ! ddor: Extension of filename to write
    integer, intent(in), optional :: label
    character(len=*), optional, intent(in) :: ngwfs_type

    ! Local variables
    character(len=256) :: output_file  ! outputr file name buffer
    !real(kind=DP), allocatable, dimension(:,:,:) :: uni_tbox  ! universal tightbox
    ! agrecocmplx: use FFTBOX_DATA type instead
    type(FFTBOX_DATA) :: uni_tbox  ! universal tightbox
    real(kind=DP) :: tb_orig1   ! a1-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: tb_orig2   ! a2-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: tb_orig3   ! a3-position of NGWF wrt tightbox in (real number of) grid points
    integer :: maxt_n1     ! a1-maximum tightbox points
    integer :: maxt_n2     ! a2-maximum tightbox points
    integer :: maxt_n3     ! a3-maximum tightbox points
    integer :: orig_atom   ! global atom counter in input file order
    integer :: atom_ngwf   ! ngwf of current atom counter
    integer :: output_unit ! fortran output unit number
    integer :: distr_ngwf  ! global ngwf counter in spacefil order
    real(kind=DP) :: sendbuf(3) ! comms send pack buffer
    real(kind=DP) :: recvbuf(3) ! comms recv pack buffer
    real(kind=DP) :: max_radius ! Maximum NGWF radius
    type(FUNC_BASIS) :: tmp_basis ! ndmh: for dealing with FFTbox coinciding with cell
    logical :: orig_coin1, orig_coin2, orig_coin3 ! ndmh: whether the fftbox
                           ! originally coincided with cell in 1,2,3 directions
    type(FFTBOX_INFO) :: tmp_fftbox ! jd: To preserve const-correctness on fftbox arg

    ! gcc32:
    real(kind=DP) :: tb_orig_buf(3) ! MPI-IO buffer
    logical :: use_tmp_basis

    logical :: loc_cmplx ! agrecocmplx: local variable for complex NGWFs

    ! jcap: par pointer
    type(PARAL_INFO), pointer :: par

#ifdef HDF5
    integer :: ierr        ! memory allocation error flag

    character(50) :: dsetname, str1

#if defined(MPI) && !defined(NOMPIIO)
    integer(hid_t) :: plist_id, plist2_id, plist3_id
#endif
    integer(hid_t) :: filesystem_id
    integer(hid_t) :: file_id
    integer(hid_t) :: group_id(regions%par%nat)
    integer(hid_t) :: attr_id
    integer(hid_t) :: aspace_id1, aspace_id3
    integer(hid_t) :: dset_id, dset2_id, dset3_id
    integer(hid_t) :: dspace_id
    integer(hsize_t) :: dims(3), adims(1)

    integer :: maxt_n(3)
    integer :: iscmplx_int
#else
    ! gcc32:
    integer :: fh
    integer(kind=LONG) :: offset
    integer :: baseint
    integer :: tempsize

    integer :: tightbox_version    ! jme: allowed values are module parameters
    logical :: use_old_tightbox_format ! jme: are we using the old tightbox format?
#endif

    ! -------------------------------------------------------------------------

    if(pub_debug_on_root) write(stdout,*) 'DEBUG: Entering restart_ngwfs_tightbox_output'

    call timer_clock("restart_ngwfs_tightbox_output", 1)

    par => regions%par

    ! agrecocmplx: set loc_cmplx value
    loc_cmplx = ngwfs_on_grid%iscmplx

    ! jd: Local copy to preserve const-correctness of fftbox arg
    tmp_fftbox = fftbox

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
       ! jd: Do this on a copy so that fftbox can be intent(in)
       if (magnitude(fftbox%a1) > (2.0_DP*max_radius + &
            2.0_DP * fftbox%d1 * cell%n_pt1)) &
            tmp_fftbox%coin1 = .false.
       if (magnitude(fftbox%a2) > (2.0_DP*max_radius + &
            2.0_DP * fftbox%d2 * cell%n_pt2)) &
            tmp_fftbox%coin2 = .false.
       if (magnitude(fftbox%a3) > (2.0_DP*max_radius + &
            2.0_DP * fftbox%d3 * cell%n_pt3)) &
            tmp_fftbox%coin3 = .false.

    end if

    ! ndmh: if any of the coin1,2,3 variables have been changed...
    if ((orig_coin1.neqv.tmp_fftbox%coin1).or. &
         (orig_coin2.neqv.tmp_fftbox%coin2).or. &
         (orig_coin3.neqv.tmp_fftbox%coin3)) then

       ! ndmh: create a temporary function basis where tightboxes
       ! ndmh: do not coincide with cell if possible
       call function_basis_allocate(tmp_basis,ngwf_basis%num, &
            'tmp_'//ngwf_basis%name,par)
       call function_basis_distribute(tmp_basis,elements,par)
       call function_basis_copy_spheres(tmp_basis,ngwf_basis,par=par)
       call function_basis_init_tight_boxes(tmp_basis,tmp_fftbox,cell, par=par)

       ! ndmh: set a universal tightbox size which depends on the actual size
       ! ndmh: of the temporary NGWF tightboxes
       maxt_n1 = maxval(tmp_basis%tight_boxes(:)%tight_pts1)
       maxt_n2 = maxval(tmp_basis%tight_boxes(:)%tight_pts2)
       maxt_n3 = maxval(tmp_basis%tight_boxes(:)%tight_pts3)
       call comms_reduce('MAX',maxt_n1)
       call comms_reduce('MAX',maxt_n2)
       call comms_reduce('MAX',maxt_n3)

       use_tmp_basis = .true.

    else
       use_tmp_basis = .false.
    end if

    ! cks: allocate universal tightbox buffer
    ! agrecocmplx: allocate uni_tbox with appropriate routine
    call data_fftbox_alloc(uni_tbox, maxt_n1, maxt_n2, maxt_n3, &
         iscmplx=loc_cmplx)

    if (.not. present(label)) then
       write(output_file,'(a256)') trim(pub_rootname)//"."//&
            trim(restart_write_filename_prefix)//trim(file_extension)
       output_file = adjustl(output_file)
    else
       if (present(ngwfs_type)) then
          write(output_file,'(a256)') trim(pub_rootname)//&
           &".history.tightbox_ngwfs"//trim(utils_int_to_str(label))//&
           "."//trim(ngwfs_type)
       end if
       output_file = adjustl(output_file)
    end if
#ifdef HDF5
    output_file = trim(output_file)//'.h5'
#endif

    ! cks: root proc allocates buffer, writes basic info
    if (pub_on_root) then

       ! cks: get a unit number that is free
       output_unit = utils_unit()
       !open(unit=output_unit, form="unformatted" ,file=trim(output_file), &
       !     action="write", position='rewind', access='sequential' )

       ! ddor: Writeout depends on the filename extension.
       ! rc2013: allow for possibility of embedding subscript in filename
       if (trim(file_extension(1:14)) .eq. 'tightbox_ngwfs') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(3a)',advance ='no') &
               'Writing NGWFs to file "', trim(output_file),'"....'
       elseif  (trim(file_extension(1:18)) .eq. 'tightbox_hub_projs') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(3a)',advance ='no') &
               'Writing HUBBARD PROJs to file "', trim(output_file),'"....'
       else if (trim(file_extension(1:19)) .eq. 'tightbox_ngwfs_cond') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(3a)',advance ='no') &
               'Writing Conduction NGWFs to file "', trim(output_file),'"....'
       else if (trim(file_extension(1:22)) .eq. 'tightbox_ngwfs_history') then
          if (pub_output_detail >= NORMAL ) then
              write(stdout,'(/a)') 'Writing NGWFs history to file '
              write(stdout,'(3a)',advance='no') '"',trim(output_file),'"....'
          end if
       endif

#ifdef HDF5
       call h5open_f(ierr)
       call h5fcreate_f(trim(output_file),h5f_acc_trunc_f,filesystem_id,ierr)
       call h5gcreate_f(filesystem_id,'tightbox_ngwfs',file_id,ierr)

       ! Add version tag
       adims=(/1/)
       call h5screate_simple_f(1,  adims, aspace_id1, ierr)
       call h5acreate_f(file_id, 'VERSION', h5t_native_integer, &
            aspace_id1, attr_id, ierr)
       call h5awrite_f(attr_id, h5t_native_integer, CURRENT_TIGHTBOX_FORMAT, &
            adims, ierr)
       call h5aclose_f(attr_id, ierr)

       adims=(/1/)
       call h5screate_simple_f(1,adims,aspace_id1,ierr)
       call h5acreate_f(file_id,'num_atoms',h5t_native_integer, &
            aspace_id1,attr_id,ierr)
       call h5awrite_f(attr_id,h5t_native_integer,regions%par%nat,adims,ierr)
       call h5aclose_f(attr_id,ierr)
       call h5acreate_f(file_id,'num_ngwfs',h5t_native_integer, &
            aspace_id1,attr_id,ierr)
       call h5awrite_f(attr_id,h5t_native_integer,ngwf_basis%num,adims,ierr)
       call h5aclose_f(attr_id,ierr)
       adims=(/3/)
       call h5screate_simple_f(1,adims,aspace_id3,ierr)
       call h5acreate_f(file_id,'tightbox_dims',h5t_native_integer, &
            aspace_id3,attr_id,ierr)
       maxt_n=(/maxt_n1,maxt_n2,maxt_n3/)
       call h5awrite_f(attr_id,h5t_native_integer,maxt_n,adims,ierr)
       call h5aclose_f(attr_id,ierr)

       ! jme: store in an integer info on real/complex type
       adims=(/1/)
       if (loc_cmplx) then
          iscmplx_int = HDF5_TRUE_INT
       else
          iscmplx_int = HDF5_FALSE_INT
       end if
       call h5acreate_f(file_id,'iscmplx',h5t_native_integer, &
            aspace_id1,attr_id,ierr)
       call h5awrite_f(attr_id,h5t_native_integer,iscmplx_int,adims,ierr)
       call h5aclose_f(attr_id,ierr)

       dims=(/maxt_n1,maxt_n2,maxt_n3/)
       call h5screate_simple_f(3,dims,dspace_id,ierr)
#else
       open(unit=output_unit, form="unformatted" ,file=trim(output_file), &
            action="write" )

       ! jme: The new tightbox format is enforced for complex tightboxes.
       ! For real tightboxes we also use the new format, unless the appropriate
       ! devel code is used.
       use_old_tightbox_format = utils_devel_code(.false., 'RESTART', &
            'OLD_TIGHTBOX_FORMAT', pub_devel_code, no_bcast=.true.)
       if ((.not.loc_cmplx).and.use_old_tightbox_format) then
          tightbox_version = OLD_TIGHTBOX_FORMAT
       else
          tightbox_version = TIGHTBOX_FORMAT_2017
       end if
       select case (tightbox_version)
          case (TIGHTBOX_FORMAT_2017)
             ! Headers added in 2017 tightbox format.
             write(output_unit, err=140) 'VERSION_' ! jme: 8 characters in order
                                                    ! to avoid alignment issues
             write(output_unit, err=140) tightbox_version
             write(output_unit, err=140) loc_cmplx
          case (OLD_TIGHTBOX_FORMAT)
          case default
             call utils_abort('Error in restart_ngwfs_tightbox_output: &
                  &unknown format version.')
       end select
       ! Variables present from initial ('OLD') tightbox format.
       write(output_unit, err=140) ngwf_basis%num
       write(output_unit, err=140) maxt_n1
       write(output_unit, err=140) maxt_n2
       write(output_unit, err=140) maxt_n3
#endif

!for MPIIO close unformatted file here
#if defined(MPI) && !defined(NOMPIIO)
#ifdef HDF5
       call h5sclose_f(dspace_id,ierr)
       call h5sclose_f(aspace_id3,ierr)
       call h5sclose_f(aspace_id1,ierr)
       call h5gclose_f(file_id,ierr)
       call h5fclose_f(filesystem_id,ierr)
       call h5close_f(ierr)
#else
       close(unit =output_unit)
#endif
#endif

    end if ! (pub_on_root)

#if defined(MPI) && !defined(NOMPIIO)
#ifdef HDF5
    call h5open_f(ierr)
    call h5pcreate_f(h5p_file_access_f,plist_id,ierr)
    call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,ierr)
    call h5fopen_f(trim(output_file),h5f_acc_rdwr_f,filesystem_id,ierr, &
         access_prp=plist_id)
    call h5pclose_f(plist_id,ierr)
    call h5pcreate_f(h5p_group_access_f,plist_id,ierr)
    call h5gopen_f(filesystem_id,'tightbox_ngwfs',file_id,ierr, &
         gapl_id=plist_id)
    call h5pclose_f(plist_id,ierr)

    call h5pcreate_f(h5p_attribute_create_f,plist2_id,ierr)
    adims=(/1/)
    call h5screate_simple_f(1,adims,aspace_id1,ierr)
    adims=(/3/)
    call h5screate_simple_f(1,adims,aspace_id3,ierr)
    dims=(/maxt_n1,maxt_n2,maxt_n3/)
    call h5screate_simple_f(3,dims,dspace_id,ierr)
#else
    ! gcc32: open file using MPI-IO
    call comms_open_file(output_file,fh,2)

    ! Initial offset for MPI-IO
    ! jme: the initial offset depends on the format version,
    ! which so far is only known to the root proc.
    call comms_bcast(pub_root_proc_id, tightbox_version)
    ! All formats: variables + 8 length records
    offset = (storage_size(ngwf_basis%num) + storage_size(maxt_n1) + &
         storage_size(maxt_n2) + storage_size(maxt_n3) + &
         8 * storage_size(baseint)) / 8
    select case (tightbox_version)
       case (TIGHTBOX_FORMAT_2017)
          ! Add length of initial tags and complex variable + 6 length records.
          offset = offset + (storage_size('VERSION_') &
               + storage_size(tightbox_version) + storage_size(loc_cmplx) &
               + 6 * storage_size(baseint)) / 8
       case (OLD_TIGHTBOX_FORMAT) ! No additional offset needed.
       case default
          call utils_abort('Error in restart_ngwfs_tightbox_output: &
               &unknown format version.')
    end select

#endif
#endif

    ! cks: Loop over all atoms in the order they appear in the input file
    atom_loop: do orig_atom =1, regions%par%nat

       ngwfs_on_atom_loop: do atom_ngwf= &
            1,ngwf_basis%num_on_atom(regions%par%distr_atom(orig_atom))

          distr_ngwf = ngwf_basis%first_on_atom(regions%par%distr_atom(orig_atom)) + &
               atom_ngwf - 1

#if defined(MPI) && !defined(NOMPIIO)
          call internal_ngwf_write_mpiio
#else
          call internal_ngwf_write
#endif

       enddo ngwfs_on_atom_loop
    enddo atom_loop



    ! cks: don't go away just yet! Wait for all communications to finish!
    ! cks: also wait before deallocating uni_tbox which is used in non-blocking
    ! cks: send!
    call comms_barrier
    ! cks: clean up send stack too
    call comms_free


    ! gcc32: close file using MPI-IO
#if defined(MPI) && !defined(NOMPIIO)
#ifdef HDF5
    call h5sclose_f(dspace_id,ierr)
    call h5sclose_f(aspace_id3,ierr)
    call h5sclose_f(aspace_id1,ierr)
    call h5pclose_f(plist2_id,ierr)
    call h5pclose_f(plist_id,ierr)
    call h5gclose_f(file_id,ierr)
    call h5fclose_f(filesystem_id,ierr)
    call h5close_f(ierr)
#else
    call comms_close_file(fh)
#endif
#endif


    ! cks: root proc closes output file
    if (pub_on_root) then
       !for non-MPIIO close unformatted file here
#if !defined(MPI) || defined(NOMPIIO)
#ifdef HDF5
       call h5sclose_f(dspace_id,ierr)
       call h5sclose_f(aspace_id3,ierr)
       call h5gclose_f(file_id,ierr)
       call h5fclose_f(filesystem_id,ierr)
       call h5close_f(ierr)
#else
       close(unit =output_unit)
#endif
#endif
       if (pub_output_detail >= NORMAL ) write(stdout,*)' done'

    endif

    if (loc_cmplx .and. pub_read_real_ngwfs) then
       pub_read_real_ngwfs = .false.
       if (pub_on_root) then
          write(stdout,'(a)') 'A complex tightbox file has been written. &
               &Disabling the effects of your initial setting &
               &"READ_REAL_NGWFS : T" (from now on : F).'
       end if
    end if

    ! cks: deallocate universal tightbox buffer
    ! agrecocmplx: deallocate using appropriate routine
    call data_fftbox_dealloc(uni_tbox)

    ! ndmh: Clean up temporary basis if we created one
    if (use_tmp_basis) then
       call function_basis_deallocate(tmp_basis)
    end if

    call timer_clock("restart_ngwfs_tightbox_output", 2)

    if(pub_debug_on_root) write(stdout,*) 'DEBUG: Leaving restart_ngwfs_tightbox_output'

    return

140 call utils_abort('Problem writing to file in restart_ngwfs_tightbox_output.')

  contains

#if defined(MPI) && !defined(NOMPIIO)
    ! gcc32: for mpiio use
    subroutine internal_ngwf_write_mpiio
      use datatypes, only: data_set_to_zero
      use function_basis, only: function_basis_ppds_to_tightbox_basic

      implicit none
      logical :: local_ngwf
      integer :: unisize

      if (use_tmp_basis) then
         local_ngwf = (pub_my_proc_id==tmp_basis%proc_of_func(distr_ngwf))
      else
         local_ngwf = (pub_my_proc_id==ngwf_basis%proc_of_func(distr_ngwf))
      end if

      if (uni_tbox%iscmplx) then
         ! jd: Wrapped with int() to work around gfortran's warning about
         !     assignment of an integer(8) to integer(4), due to a bug
         !     in the kind of the storage_size intrinsic (gfortran 4.9.x).
         unisize = (int(storage_size(uni_tbox%z))/8) * size(uni_tbox%z)
      else
         unisize = (int(storage_size(uni_tbox%d))/8) * size(uni_tbox%d)
      end if

#ifdef HDF5
      call h5pcreate_f(h5p_group_create_f,plist_id,ierr)
      if (atom_ngwf==1) then
         write(str1,*) orig_atom
         dsetname="atom_"//trim(adjustl(str1))
         call h5gcreate_f(file_id,adjustl(dsetname),group_id(orig_atom),ierr, &
              gcpl_id=plist_id)
         call h5acreate_f(group_id(orig_atom),'num_ngwfs',h5t_native_integer, &
              aspace_id1,attr_id,ierr,acpl_id=plist2_id)
         adims=(/1/)
         call h5awrite_f(attr_id,h5t_native_integer, &
              ngwf_basis%num_on_atom(regions%par%distr_atom(orig_atom)),adims,ierr)
         call h5aclose_f(attr_id,ierr)
         call h5dcreate_f(group_id(orig_atom),'origin',h5t_native_double, &
              aspace_id3,dset2_id,ierr)
      end if
      write(str1,*) atom_ngwf
      if (uni_tbox%iscmplx) then
         dsetname="ngwf_"//trim(adjustl(str1))//"_real"
         call h5dcreate_f(group_id(orig_atom),adjustl(dsetname), &
              h5t_native_double,dspace_id,dset_id,ierr)
         dsetname="ngwf_"//trim(adjustl(str1))//"_imag"
         call h5dcreate_f(group_id(orig_atom),adjustl(dsetname), &
              h5t_native_double,dspace_id,dset3_id,ierr)
      else
         dsetname="ngwf_"//trim(adjustl(str1))
         call h5dcreate_f(group_id(orig_atom),adjustl(dsetname), &
              h5t_native_double,dspace_id,dset_id,ierr)
      end if
#endif

      if (local_ngwf) then

         ! cks: initialise
         ! agrecocmplx: initialise according to uni_tbox real or complex
         call data_set_to_zero(uni_tbox)
         tb_orig_buf =0.0_DP

         if (use_tmp_basis) then
            call function_basis_ppds_to_tightbox_basic(uni_tbox, &
                 tb_orig1,tb_orig2,tb_orig3,distr_ngwf, &
                 ngwfs_on_grid,tmp_basis,tmp_fftbox,cell)
         else
            call function_basis_ppds_to_tightbox_basic(uni_tbox, &
                 tb_orig1,tb_orig2,tb_orig3,distr_ngwf, &
                 ngwfs_on_grid,ngwf_basis,tmp_fftbox, cell)
         end if

         tb_orig_buf(1) = tb_orig1
         tb_orig_buf(2) = tb_orig2
         tb_orig_buf(3) = tb_orig3

#ifdef HDF5
         call h5pcreate_f(h5p_dataset_xfer_f,plist3_id,ierr)
         call h5pset_dxpl_mpio_f(plist3_id,h5fd_mpio_independent_f,ierr)
         if (atom_ngwf==ngwf_basis%num_on_atom(regions%par%distr_atom(orig_atom))) then
            adims=(/3/)
            call h5dwrite_f(dset2_id,h5t_native_double,tb_orig_buf,adims,ierr, &
                 xfer_prp=plist3_id)
         end if
         if (uni_tbox%iscmplx) then
            call h5dwrite_f(dset_id,h5t_native_double,real(uni_tbox%z,DP),dims,ierr, &
                 xfer_prp=plist3_id)
            call h5dwrite_f(dset3_id,h5t_native_double,aimag(uni_tbox%z),dims,ierr, &
                 xfer_prp=plist3_id)
         else
            call h5dwrite_f(dset_id,h5t_native_double,uni_tbox%d,dims,ierr, &
                 xfer_prp=plist3_id)
         end if
         call h5pclose_f(plist3_id,ierr)
#else
         ! gcc32: get record length
         tempsize = (storage_size(tb_orig_buf)/8) * size(tb_orig_buf)
         call comms_write(fh,offset,tempsize)

         ! gcc32: adjust offset
         offset = offset + storage_size(baseint)/8

         call comms_write(fh,offset,tb_orig_buf)

         ! gcc32: adjust offset
         offset = offset + (storage_size(tb_orig_buf)/8) * size(tb_orig_buf)

         call comms_write(fh,offset,tempsize)

         ! gcc32: adjust offset
         offset = offset + storage_size(baseint)/8

         ! gcc32: get record length
         call comms_write(fh,offset,unisize)

         ! gcc32: adjust offset
         offset = offset + storage_size(baseint)/8

         if (uni_tbox%iscmplx) then
            call comms_write(fh,offset,uni_tbox%z)
         else
            call comms_write(fh,offset,uni_tbox%d)
         end if

         ! gcc32: adjust offset
         offset = offset + unisize

         call comms_write(fh,offset,unisize)

         ! gcc32: adjust offset
         offset = offset + storage_size(baseint)/8

      else

         ! gcc32: adjust offset
         offset = offset + unisize + &
              (storage_size(tb_orig_buf)/8)*size(tb_orig_buf) + &
              4 * storage_size(baseint)/8

#endif

      end if ! (local_ngwf)

#ifdef HDF5
      call h5dclose_f(dset_id,ierr)
      if (uni_tbox%iscmplx) call h5dclose_f(dset3_id,ierr)
      if (atom_ngwf==ngwf_basis%num_on_atom(regions%par%distr_atom(orig_atom))) then
         call h5dclose_f(dset2_id,ierr)
         call h5gclose_f(group_id(orig_atom),ierr)
      end if
#endif

    end subroutine internal_ngwf_write_mpiio

#endif

    ! gcc32: for non-mpiio use
    subroutine internal_ngwf_write
      use function_basis, only: function_basis_ppds_to_tightbox

      implicit none

      ! cks: this essentially converts the non-blocking sends of
      ! cks: comms_send to blocking sends. Since parallel performance
      ! cks: is not the issue here, this makes the code simpler.
      call comms_barrier

      ! ndmh: transfer function to buffer
      if (use_tmp_basis) then
         call function_basis_ppds_to_tightbox(uni_tbox, &
              tb_orig1,tb_orig2,tb_orig3,distr_ngwf,pub_root_proc_id, &
              ngwfs_on_grid,tmp_basis,tmp_fftbox,cell,sendbuf,recvbuf)
      else
         call function_basis_ppds_to_tightbox(uni_tbox, &
              tb_orig1,tb_orig2,tb_orig3,distr_ngwf,pub_root_proc_id, &
              ngwfs_on_grid,ngwf_basis,tmp_fftbox,cell,sendbuf,recvbuf)
      end if

      ! ndmh: write buffer to disk
      if (pub_on_root) then
#ifdef HDF5
         if (atom_ngwf==1) then
            write(str1,*) orig_atom
            dsetname="atom_"//trim(adjustl(str1))
            call h5gcreate_f(file_id,adjustl(dsetname),group_id(orig_atom),ierr)
            call h5acreate_f(group_id(orig_atom),'num_ngwfs', &
                 h5t_native_integer,aspace_id1,attr_id,ierr)
            adims=(/1/)
            call h5awrite_f(attr_id,h5t_native_integer, &
                 ngwf_basis%num_on_atom(regions%par%distr_atom(orig_atom)),adims,ierr)
            call h5aclose_f(attr_id,ierr)
            call h5dcreate_f(group_id(orig_atom),'origin',h5t_native_double, &
                 aspace_id3,dset_id,ierr)
            tb_orig_buf=(/tb_orig1,tb_orig2,tb_orig3/)
            adims=(/3/)
            call h5dwrite_f(dset_id,h5t_native_double,tb_orig_buf,adims,ierr)
            call h5dclose_f(dset_id,ierr)
         end if
         write(str1,*) atom_ngwf
         if (uni_tbox%iscmplx) then
            dsetname="ngwf_"//trim(adjustl(str1))//"_real"
            call h5dcreate_f(group_id(orig_atom),adjustl(dsetname), &
                 h5t_native_double,dspace_id,dset_id,ierr)
            dsetname="ngwf_"//trim(adjustl(str1))//"_imag"
            call h5dcreate_f(group_id(orig_atom),adjustl(dsetname), &
                 h5t_native_double,dspace_id,dset3_id,ierr)
            call h5dwrite_f(dset_id,h5t_native_double,real(uni_tbox%z,DP),dims,ierr)
            call h5dwrite_f(dset3_id,h5t_native_double,aimag(uni_tbox%z),dims,ierr)
         else
            dsetname="ngwf_"//trim(adjustl(str1))
            call h5dcreate_f(group_id(orig_atom),adjustl(dsetname), &
                 h5t_native_double,dspace_id,dset_id,ierr)
            call h5dwrite_f(dset_id,h5t_native_double,uni_tbox%d,dims,ierr)
         end if
         call h5dclose_f(dset_id,ierr)
         if (uni_tbox%iscmplx) call h5dclose_f(dset3_id,ierr)
         if (atom_ngwf==ngwf_basis%num_on_atom(regions%par%distr_atom(orig_atom))) &
              call h5gclose_f(group_id(orig_atom),ierr)
#else
         write(output_unit, err=140) tb_orig1, tb_orig2, tb_orig3

         ! agrecocmplx: distinguish between real and complex cases
         if (uni_tbox%iscmplx) then
             write(output_unit, err=140) uni_tbox%z(1:maxt_n1,1:maxt_n2,1:maxt_n3)
         else
             write(output_unit, err=140) uni_tbox%d(1:maxt_n1,1:maxt_n2,1:maxt_n3)
         end if
#endif
      end if


      return

140   call utils_abort('Problem writing to file in restart_ngwfs_tightbox_output.')

    end subroutine internal_ngwf_write

  end subroutine restart_ngwfs_tightbox_output


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine restart_ngwfs_tightbox_input(ngwfs_on_grid, ngwf_basis, &
       cell, fftbox, elements, file_extension, regions, ref_dir, label, &
       ngwf_type, ref_is_filename, restart_rootname)

    !=======================================================================!
    ! This subroutine reads the current NGWFs in "universal tightbox"       !
    ! representation from a .tightbox_ngwfs file. The order of the NGWFs    !
    ! in the .tightbox_ngwfs file is the one in which they appear in the    !
    ! input file and is not affected by the use of the space-filling curve. !
    ! Most importantly the "universal tightbox" storage allows              !
    ! the .tightbox_ngwfs file of a previous atomic                         !
    ! configuration to be read, i.e. one can initialise the NGWFs with      !
    ! the converged NGWFs of a previous calculation where the atoms         !
    ! were in different positions.                                          !
    ! ndmh 18/12/2009: In cases where the FFTbox coincides with the cell    !
    ! but the individual tightboxes do not need to, a new temporary NGWF    !
    ! basis is constructed, for which the FFTbox only coincides with the    !
    ! cell along direction i if L_i < 2*R_NGWF + delta_i.                   !
    !-----------------------------------------------------------------------!
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 12/3/2004                         !
    ! Modified by Chris-Kriton Skylaris on 13/4/2004 so that the read NGWFs !
    ! are shifted to the centres of their corresponding atoms.              !
    ! Fix by Chris-Kriton Skylaris on 16/02/2005 of bug related to          !
    ! initialisation of tb_start1/2/3 reported by Victor Milman.            !
    ! Modified by D. O'Regan in November 2009 for different file extensions.!
    ! Modified by Nicholas Hine in December 2009, to allow reduced file     !
    ! sizes by only reading in tightboxes even if FFTbox coincides with     !
    ! cell in some directions.                                              !
    ! Modified by Gabriel Constantinescu in June 2014 to use MPI-IO to read !
    ! the tightboxes, instead of sending them from the root proc.           !
    ! Modified by Andrea Greco in May 2015 to allow use of complex NGWFs.   !
    ! Modified by Andrea Greco and Jose M Escartin in January/February 2017 !
    ! to enable reading real NGWFs into complex NGWFs.                      !
    ! Modified by Jose M Escartin in January 2017 to incorporate new format !
    ! for tightbox files.                                                   !
    ! Const-correctness for 'fftbox' fixed by Jacek Dziedzic in June 2017.  !
    ! Modified to remove pub_par by Joseph Prentice, May 2018               !
    !=======================================================================!

    use basis, only: basis_ket_start_wrt_fftbox
    use comms, only: comms_barrier, comms_bcast, comms_free, &
#if defined(MPI) && !defined(NOMPIIO)
         comms_open_file, comms_close_file, comms_read, &
#endif
         comms_reduce, pub_on_root, pub_root_proc_id, pub_my_proc_id
    use constants, only: DP, stdout, LONG
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc, data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS, function_basis_allocate, &
         function_basis_deallocate, function_basis_distribute, &
         function_basis_copy_spheres, function_basis_init_tight_boxes
    use geometry, only: magnitude
    use ion, only : ELEMENT
    use model_type, only: REGION
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_rootname, pub_devel_code, pub_read_real_ngwfs, &
         pub_qnto_analysis
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_assert, utils_abort, utils_devel_code, utils_int_to_str, &
         utils_feature_not_supported

    implicit none

    ! Arguments
    type(REGION), target, intent(in) :: regions
    type(FUNC_BASIS), intent(in)   :: ngwf_basis
    !real(kind=DP), intent(out) :: ngwfs_on_grid(ngwf_basis%size_on_grid) ! NGWFs on this proc
    ! agrecocmplx: use FUNCTIONS type instead
    type(FUNCTIONS), intent(inout) :: ngwfs_on_grid ! NGWFs on this proc
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(ELEMENT), intent(in) :: elements(regions%par%nat) ! elements of all proc (in input file order)
    character(len=*), intent(in) :: file_extension ! ddor: Extension of filename to read
    character(len=*), optional, intent(in) :: ref_dir ! jhl52
    integer, optional, intent(in) :: label
    character(len=*), optional, intent(in) :: ngwf_type
    logical, optional, intent(in) :: ref_is_filename ! mjsp: Modify ref_dir to
                                               ! refer to a filename, rather than a directory
    character(len=*),optional,intent(in) :: restart_rootname

    ! Local Variables
    character(len=256) :: read_file  ! read file name buffer
    ! agrecocmplx: use FFTBOX_DATA type instead
    type(FFTBOX_DATA) :: uni_tbox  ! universal tightbox
    complex(kind=DP), allocatable, dimension(:,:,:) :: fftbox_complex ! complex fftbox space
    complex(kind=DP), allocatable, dimension(:,:,:) :: fftbox_complex_shifted ! complex fftbox space
    !real(kind=DP), allocatable, dimension(:,:,:) :: fftbox_buffer     ! real fftbox space
    ! agrecocmplx: use FFTBOX_DATA type instead
    type(FFTBOX_DATA) :: fftbox_buffer     ! real fftbox space
    type(FFTBOX_INFO) :: tmp_fftbox ! jd: For preserving intent(in) on fftbox arg
    real(kind=DP) :: max_radius          ! globaly maximum ngwf radius
    real(kind=DP) :: read_tb_orig1 ! read a1-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: read_tb_orig2 ! read a2-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: read_tb_orig3 ! read a3-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: read_weight   ! read grid point weight
    integer :: n1, n2, n3, ld1, ld2 ! fftbox dimensions
    integer :: orig_atom    ! global atom counter in input file order
    integer :: atom_ngwf    ! ngwf of current atom counter
    integer :: input_unit   ! fortran input unit number
    integer :: read_num     ! total number of ngwfs read from file
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
    logical :: read_bigendian
    logical :: file_iscmplx ! is the data in the file real or complex?

    ! ndmh:
    type(FUNC_BASIS) :: tmp_basis
    logical :: use_tmp_basis
    logical :: orig_coin1,orig_coin2,orig_coin3 ! ndmh: whether the fftbox
                           ! originally coincided with cell in 1,2,3 directions
    ! gcc32:
    integer :: fh
    real(kind=DP) :: read_tb_buf(3) ! MPI-IO buffer
#ifdef HDF5
    character(50) :: dsetname, str1

    real(kind=DP), allocatable :: buff_real(:,:,:), buff_imag(:,:,:)

#if defined(MPI) && !defined(NOMPIIO)
    integer(hid_t) :: plist_id, plist3_id
#endif
    integer(hid_t) :: filesystem_id
    integer(hid_t) :: file_id
    integer(hid_t) :: group_id(regions%par%nat)
    integer(hid_t) :: attr_id
    integer(hid_t) :: dset_id, dset2_id, dset3_id
    integer(hsize_t) :: dims(3), adims(1)

    integer :: read_maxt_n(3)
#else
    type(FFTBOX_DATA) :: temp_uni_tbox ! temporary buffer
    integer(kind=LONG) :: offset
    integer :: baseint
#endif

    logical :: loc_cmplx        ! flag for reading into complex/real structures
    logical :: read_cmplx       ! jme: flag for reading from a complex/real file
    integer :: tightbox_version ! jme: allowed values are module parameters
#ifdef HDF5
    integer :: iscmplx_int
#else
    character(len=len('VERSION_')) :: tag
    logical :: use_old_tightbox_format ! jme: are we using the old tightbox format?
#endif

    ! -------------------------------------------------------------------------

    ! jcap: par pointer
    type(PARAL_INFO), pointer :: par

    call timer_clock("restart_ngwfs_tightbox_input", 1)

    par => regions%par

    ! jme: set complex character of subroutine output
    loc_cmplx = ngwfs_on_grid%iscmplx
    ! jme: do we want to read a real or a complex file?
    read_cmplx = ( loc_cmplx .and. (.not. pub_read_real_ngwfs) )

    ! cks: ********** INITIALISATIONS ***************************************
    ! agrecocmplx: initialise FUNCTIONS type in the appropriate way
    call data_set_to_zero(ngwfs_on_grid)

    ! jd: Local copy to preserve const-correctness of fftbox arg
    tmp_fftbox = fftbox

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

    read_maxt_n1 =0
    read_maxt_n2 =0
    read_maxt_n3 =0
    ! cks: ****** END INITIALISATIONS ***************************************

    read_file = repeat(' ',len(read_file))

    if (.not. present(label)) then

       if (present(restart_rootname)) then
          ! kkbd: Read in from a custom rootname
          write(read_file,'(a256)') trim(restart_rootname)//"."//trim(file_extension)
       else
          write(read_file,'(a256)') trim(pub_rootname)//"."//&
               trim(restart_read_filename_prefix)//trim(file_extension)
       end if

       if(present(ref_dir)) then ! jhl52

          write(read_file,'(a)') trim(ref_dir)//'/'//trim(pub_rootname)//&
          "."//trim(restart_read_filename_prefix)//trim(file_extension)

          ! mjsp: If ref_dir refers to a filename rather than a directory,
          ! i.e. does not use pub_rootname in the filename.
          ! (usage: for loading in fragment NGWFs, possibly calculated within
          ! a separate calculation folder with different filenames)
          if (present(ref_is_filename)) then
             if (ref_is_filename) then
               write(read_file,'(a)') trim(ref_dir)//"."//trim(file_extension)
             end if
          end if

       end if ! ref_dir

    else

       if (present(ngwf_type)) then
          write(read_file,'(a256)') &
             trim(pub_rootname)//".history.tightbox_ngwfs"//&
             &trim(utils_int_to_str(label))//"."//trim(ngwf_type)
          read_file = trim(read_file)
       end if

    end if

    read_file = adjustl(read_file)
#ifdef HDF5
    read_file = trim(read_file)//'.h5'
#else
    ! cks: get a unit number that is free
    input_unit =utils_unit()
#endif

    ! Check if we want to read big-endian formatted data
    read_bigendian = utils_devel_code(.false., 'RESTART','BIGENDIAN', &
       pub_devel_code)

    ! cks: root proc opens input file, reads basic info
    if (pub_on_root) then

       ! ddor: Writeout depends on the filename extension.
       ! rc2013: allow for possibility of embedding subscript in filename
       if (trim(file_extension(1:14)) .eq. 'tightbox_ngwfs') then
          write(stdout,'(/3a)',advance='no') &
               'Reading NGWFs from file "', trim(read_file),'" ...'
       elseif (trim(file_extension(1:18)) .eq. 'tightbox_hub_projs') then
          write(stdout,'(/3a)',advance='no') &
               'Reading HUBBARD PROJs from file "', trim(read_file),'" ...'
       elseif (trim(file_extension(1:19)) .eq. 'tightbox_ngwfs_cond') then
          write(stdout,'(/3a)',advance='no') &
               'Reading conduction NGWFs from file "', trim(read_file),'" ...'
       elseif (trim(file_extension(1:21)) .eq. 'vacuum_tightbox_ngwfs') then
          write(stdout,'(/a/3a)',advance='no') &
               'Reading NGWFs for initialising the solvent cavity from', &
               'file "',trim(read_file),'" ...'
       elseif (trim(file_extension(1:22)) .eq. 'tightbox_ngwfs_history') then
          write(stdout,'(/a)') 'Reading NGWFs history from file '
          write(stdout,'(3a)',advance ='no') '"',trim(read_file),'" ...'
       endif

#ifdef HDF5
       call h5open_f(ierr)
       call h5fopen_f(trim(read_file),h5f_acc_rdonly_f,filesystem_id,ierr)
       call h5gopen_f(filesystem_id,'tightbox_ngwfs',file_id,ierr)

       ! Obtain version tag.
       adims=(/1/)
       call h5aopen_f(file_id,'VERSION', attr_id, ierr)
       call h5aread_f(attr_id, h5t_native_integer, tightbox_version, adims, &
            ierr)
       select case (ierr)
          case (0)
             if (pub_debug_on_root) then
                write(stdout,'(a,f4.2,a)') 'Opening ONETEP HDF5 tightbox file, &
                     &format v', real(tightbox_version,DP)/100.0_DP, ' .'
             end if
          case (-1)
             tightbox_version = OLD_TIGHTBOX_FORMAT ! = 100
             write(stdout,'(a,f4.2)') 'Opening ONETEP HDF5 tightbox file, &
                  &format not declared so assuming format v1.00 .'
          case default
             call utils_abort('Error in restart_ngwfs_tightbox_input: &
                  &unknown error in HDF5 h5aread_f call.')
       end select
       call h5aclose_f(attr_id, ierr)

       adims=(/1/)
       call h5aopen_f(file_id,'num_ngwfs',attr_id,ierr)
       call h5aread_f(attr_id,h5t_native_integer,read_num,adims,ierr)
       call h5aclose_f(attr_id,ierr)
       adims=(/3/)
       call h5aopen_f(file_id,'tightbox_dims',attr_id,ierr)
       call h5aread_f(attr_id,h5t_native_integer,read_maxt_n,adims,ierr)
       call h5aclose_f(attr_id,ierr)
       read_maxt_n1=read_maxt_n(1)
       read_maxt_n2=read_maxt_n(2)
       read_maxt_n3=read_maxt_n(3)

       ! jme: read info on real/complex type from integer if available.
       if (tightbox_version .ge. TIGHTBOX_FORMAT_2017) then
          adims=(/1/)
          call h5aopen_f(file_id,'iscmplx', attr_id, ierr)
          call h5aread_f(attr_id, h5t_native_integer, iscmplx_int, adims, ierr)
          call utils_assert(ierr==0, 'Error in restart_ngwfs_tightbox_input: &
               &error reading real/complex tag.')
          call h5aclose_f(attr_id, ierr)
          select case (iscmplx_int)
             case(HDF5_FALSE_INT)
                file_iscmplx = .false.
             case(HDF5_TRUE_INT)
                file_iscmplx = .true.
             case default
                call utils_abort('Error in restart_ngwfs_tightbox_input: &
                     &invalid value of iscmplx_int: ', iscmplx_int)
          end select
       else
          file_iscmplx = .false.
       end if

       dims=read_maxt_n
#else
       if (.not.read_bigendian) then
          open(unit =input_unit, form="unformatted" ,file=trim(read_file), &
               action="read", status='old' )
       else
#ifndef F2008
          ! jd: 'convert' is non-standard Fortran.
          write(stdout,*) 'WARNING: Over-riding data format to big-endian &
               &when reading NGWFs'
          write(stdout,*) 'WARNING: and skipping reading number of NGWFs and &
               &tightbox size'
          open(unit =input_unit, form="unformatted" ,file=trim(read_file), &
               action="read", status='old',convert='big_endian')
#else
          call utils_feature_not_supported('restart_ngwfs_tightbox_input:'//&
               'Big-endian input has not been compiled in because it &
               &needs "convert", which is not part of the F2008 standard, &
               &and you asked for an F2008-compliant binary with -D2008.')
#endif

       end if

       ! jme: The new tightbox file format is enforced for complex tightboxes.
       ! For real tightboxes we also use the new format, unless the appropriate
       ! devel code is used.
       ierr = 0
       use_old_tightbox_format = utils_devel_code(.false., 'RESTART', &
            'OLD_TIGHTBOX_FORMAT', pub_devel_code, no_bcast=.true.)
       if (read_cmplx .or. (.not.use_old_tightbox_format)) then
          if (pub_debug_on_root) then
             write(stdout,*) 'Attempting to read new tightbox format file.'
          end if
          read(input_unit, iostat=ierr) tag
          if ((ierr==0) .and. (tag == 'VERSION_')) then
             read(input_unit, err=150) tightbox_version
             if (pub_debug_on_root) then
                write(stdout,'(a,f4.2,a)') &
                     ' ONETEP tightbox format detected v', &
                     real(tightbox_version,DP)/100.0_DP, ' .'
             end if
             read(input_unit, err=150) file_iscmplx
          else
             write(*,*) 'Tightbox format version > v1.00 expected but not &
                  &recognised.'
             ierr = -1
             rewind input_unit
          end if
       end if
       if ( ((.not.read_cmplx).and.use_old_tightbox_format) .or. (ierr.ne.0) ) then
          write(stdout,*) 'Attempting to read old tightbox format file.'
          ! jme: Complex NGWFs were consolidated at the same time that the new
          ! tightbox file format was introduced (and enforced for them), so
          ! nobody should be using files with complex tightboxes
          ! in the old format.
          file_iscmplx = .false.
          tightbox_version = OLD_TIGHTBOX_FORMAT
       end if
       ! jme: Read headers present in all tightbox formats.
       read(input_unit, err=150) read_num
       read(input_unit, err=150) read_maxt_n1
       read(input_unit, err=150) read_maxt_n2
       read(input_unit, err=150) read_maxt_n3

#endif

!for MPIIO close unformatted file here
#if defined(MPI) && !defined(NOMPIIO)
#ifdef HDF5
       call h5gclose_f(file_id,ierr)
       call h5fclose_f(filesystem_id,ierr)
       call h5close_f(ierr)
#else
       close(unit =input_unit)
#endif
#endif

       ! jme: check that the file actually contains the type of data
       ! (real/complex) that we want to read
       call utils_assert(read_cmplx .eqv. file_iscmplx, 'Error in &
            &restart_ngwfs_tightbox_input: Inconsistency in '// &
            trim(file_extension) //', complex character of data &
            &expected to be read differs from that detected in the &
            &tightbox file. Their respective complex attributes are: ', &
            read_cmplx, file_iscmplx)

       ! cks: check that the correct number of NGWFs will be read from file
       if (.not.pub_qnto_analysis) then
          call utils_assert(read_num == ngwf_basis%num, &
               'Error in restart_ngwfs_tightbox_input(): Inconsistency in '//&
               trim(file_extension)//', read_num and ngwfs_basis%num were ', &
               read_num, ngwf_basis%num)
       else
         ! jhl52: In QNTO, the number of REF-NGWFs can be different from that of SYS-NGWFs
         call utils_assert(read_num >= ngwf_basis%num, &
              'Error in restart_ngwfs_tightbox_input(): Inconsistency in '//&
              trim(file_extension)//', read_num and ngwfs_basis%num were ', &
              read_num, ngwf_basis%num)
       endif
    endif

    ! ndmh: broadcast the read_maxt's from the root to all other procs, since
    ! ndmh: they are needed to determine whether to initialise temporary basis
    call comms_bcast(pub_root_proc_id, read_maxt_n1)
    call comms_bcast(pub_root_proc_id, read_maxt_n2)
    call comms_bcast(pub_root_proc_id, read_maxt_n3)

    ! ndmh: if the FFTbox coincides with the cell in any direction, but the
    ! ndmh: tightbox does not necessarily need to, we may need to define a
    ! ndmh: temporary NGWF basis where the tightboxes are as small as they
    ! ndmh: can be, unless the file we are reading in was written with
    ! ndmh: tightboxes coinciding with cell
    ! jd: Modified to maintain const-correctness of fftbox arg
    orig_coin1 = fftbox%coin1
    orig_coin2 = fftbox%coin2
    orig_coin3 = fftbox%coin3
    if (orig_coin1.or.orig_coin2.or.orig_coin3) then

       ! ndmh: temporarily override coin1,2,3 variables if possible (as long
       ! ndmh: as the written values are not for tightbox==fftbox)
       if ((magnitude(fftbox%a1) > 2.0_DP*max_radius + &
            2.0_DP * fftbox%d1*cell%n_pt1).and. &
            (read_maxt_n1 /= fftbox%total_pt1)) tmp_fftbox%coin1 = .false.
       if ((magnitude(fftbox%a2) > 2.0_DP*max_radius + &
            2.0_DP * fftbox%d2*cell%n_pt2).and. &
            (read_maxt_n2 /= fftbox%total_pt2)) tmp_fftbox%coin2 = .false.
       if ((magnitude(fftbox%a3) > 2.0_DP*max_radius + &
            2.0_DP * fftbox%d3*cell%n_pt3).and. &
            (read_maxt_n3 /= fftbox%total_pt3)) tmp_fftbox%coin3 = .false.

       if(pub_debug_on_root) then
          if ((orig_coin1.neqv.tmp_fftbox%coin1)) write(stdout,*) &
               'DEBUG: Temporarily overriding fftbox%coin1'
          if ((orig_coin2.neqv.tmp_fftbox%coin2)) write(stdout,*) &
               'DEBUG: Temporarily overriding fftbox%coin2'
          if ((orig_coin3.neqv.tmp_fftbox%coin3)) write(stdout,*) &
               'DEBUG: Temporarily overriding fftbox%coin3'
       end if

    end if

    ! ndmh: if any of the coin1,2,3 variables have been changed...
    if ((orig_coin1.neqv.tmp_fftbox%coin1).or. &
         (orig_coin2.neqv.tmp_fftbox%coin2).or. &
         (orig_coin3.neqv.tmp_fftbox%coin3)) then

       ! ndmh: create a temporary function basis where tightboxes
       ! ndmh: do not coincide with cell if possible
       call function_basis_allocate(tmp_basis,ngwf_basis%num, &
            'tmp_'//ngwf_basis%name,par)
       call function_basis_distribute(tmp_basis,elements,par)
       call function_basis_copy_spheres(tmp_basis, ngwf_basis, par=par)
       call function_basis_init_tight_boxes(tmp_basis, tmp_fftbox, cell, par=par)

       ! ndmh: set a universal tightbox size which depends on the actual size
       ! ndmh: of the temporary NGWF tightboxes
       maxt_n1 = maxval(tmp_basis%tight_boxes(:)%tight_pts1)
       maxt_n2 = maxval(tmp_basis%tight_boxes(:)%tight_pts2)
       maxt_n3 = maxval(tmp_basis%tight_boxes(:)%tight_pts3)
       call comms_reduce('MAX',maxt_n1)
       call comms_reduce('MAX',maxt_n2)
       call comms_reduce('MAX',maxt_n3)

       use_tmp_basis = .true.
    else
       use_tmp_basis = .false.
    end if

    ! cks: NGWFs will be placed in the centre of the fftbox
    ! cks: or left where they are in the simulation cell
    ! cks: (when the simulation cell and FFT-box coincide)
    call basis_ket_start_wrt_fftbox(tb_start1, tb_start2, tb_start3, &
         n1, n2, n3, tmp_fftbox)

    ! ndmh: check tightbox will not extend out of FFTbox if the coin flag has
    ! ndmh: been overridden. If so, try to adjust tb_start to compensate. If
    ! ndmh: that fails, print an error and quit.
    if ((tb_start1 + maxt_n1 > n1).and.(.not.tmp_fftbox%coin1)) then

       if (pub_debug_on_root) write(stdout,*) &
            'DEBUG: Attempting to adjust tb_start1'

       if (n1 - maxt_n1 > 0) then
          tb_start1 = (n1 - maxt_n1)/2
       else
          call utils_abort('Error in restart_ngwfs_tightbox_input: &
               &Tightbox size along 1-axis is not compatible with FFT-box size.&
               & The values were: ',maxt_n1, n1)
       end if
    end if
    if ((tb_start2 + maxt_n2 > n2).and.(.not.tmp_fftbox%coin2)) then

       if (pub_debug_on_root) write(stdout,*) &
            'DEBUG: Attempting to adjust tb_start2'

       if (n2 - maxt_n2 > 0) then
          tb_start2 = (n2 - maxt_n2)/2
       else
          call utils_abort('Error in restart_ngwfs_tightbox_input: &
               &Tightbox size along 2-axis is not compatible with FFT-box size.&
               & The values were: ',maxt_n2, n2)
       end if
    end if
    if ((tb_start3 + maxt_n3 > n3).and.(.not.tmp_fftbox%coin3)) then

       if (pub_debug_on_root) write(stdout,*) &
            'DEBUG: Attempting to adjust tb_start3'

       if (n3 - maxt_n3 > 0) then
          tb_start3 = (n3 - maxt_n3)/2
       else
          call utils_abort('Error in restart_ngwfs_tightbox_input: &
               &Tightbox size along 3-axis is not compatible with FFT-box size.&
               & The values were: ',maxt_n3, n3)
       end if
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
    ! agrecocmplx: allocate uni_tbox with appropriate routine
    call data_fftbox_alloc(uni_tbox, maxt_n1, maxt_n2, maxt_n3, &
         iscmplx=loc_cmplx)
    allocate(fftbox_complex(ld1,ld2,n3),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_input','fftbox_complex',ierr)
    allocate(fftbox_complex_shifted(ld1,ld2,n3),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_input', &
         'fftbox_complex_shifted',ierr)
    ! agrecocmplx: allocate fftbox_buffer with appropriate routine
    call data_fftbox_alloc(fftbox_buffer, ld1, ld2, n3, iscmplx=loc_cmplx)

#ifdef HDF5
    allocate(buff_real(read_maxt_n1,read_maxt_n2,read_maxt_n3),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_input', 'buff_real', ierr)
    if (read_cmplx) then
       allocate(buff_imag(read_maxt_n1,read_maxt_n2,read_maxt_n3),stat=ierr)
       call utils_alloc_check('restart_ngwfs_tightbox_input', 'buff_imag', ierr)
    end if
#else
    ! gcc32: allocate temporary buffer in case maxt_... != read_maxt_...
    ! agrecocmplx: allocate temp_uni_tbox with appropriate routine
    call data_fftbox_alloc(temp_uni_tbox, read_maxt_n1, read_maxt_n2, &
             read_maxt_n3, iscmplx=read_cmplx)
#endif

#if defined(MPI) && !defined(NOMPIIO)
#ifdef HDF5
    call h5open_f(ierr)
    call h5pcreate_f(h5p_file_access_f,plist_id,ierr)
    call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,ierr)
    call h5fopen_f(trim(read_file),h5f_acc_rdonly_f,filesystem_id,ierr, &
         access_prp=plist_id)
    call h5pclose_f(plist_id,ierr)
    call h5pcreate_f(h5p_group_access_f,plist_id,ierr)
    call h5gopen_f(filesystem_id,'tightbox_ngwfs',file_id,ierr, &
         gapl_id=plist_id)

    dims=read_maxt_n
#else
    ! gcc32: open file using MPI-IO
    call comms_open_file(read_file,fh,1) ! 1 = open for read-only

    ! Initial offset for MPI-IO
    ! jme: the initial offset depends on the format version,
    ! which so far is only known to the root proc.
    call comms_bcast(pub_root_proc_id, tightbox_version)
    ! All formats: variables + 8 length records
    offset = (storage_size(ngwf_basis%num) + storage_size(maxt_n1) + &
         storage_size(maxt_n2) + storage_size(maxt_n3) + &
         8 * storage_size(baseint)) / 8
    select case (tightbox_version)
       case (TIGHTBOX_FORMAT_2017)
          ! Add length of initial tags and complex variable + 6 length records.
          offset = offset + (storage_size('VERSION_') &
               + storage_size(tightbox_version) + storage_size(loc_cmplx) &
               + 6 * storage_size(baseint)) / 8
       case (OLD_TIGHTBOX_FORMAT) ! No additional offset needed.
       case default
          call utils_abort('Error in restart_ngwfs_tightbox_input: &
               &unknown format version.')
    end select
#endif
#endif

    ! cks: Loop over all atoms in the order they appear in the input file
    atom_loop: do orig_atom =1, regions%par%nat

       ngwfs_on_atom_loop: do atom_ngwf=1, &
            ngwf_basis%num_on_atom(regions%par%distr_atom(orig_atom))

          distr_ngwf = ngwf_basis%first_on_atom(regions%par%distr_atom(orig_atom)) + &
               atom_ngwf - 1

#if defined(MPI) && !defined(NOMPIIO)
         call internal_ngwf_read_mpiio
#else
         call internal_ngwf_read
#endif

       enddo ngwfs_on_atom_loop
    enddo atom_loop

    ! cks: don't go away just yet! Wait for all communications to finish!
    call comms_barrier
    ! cks: clean up send stack too
    call comms_free

#if defined(MPI) && !defined(NOMPIIO)
    ! gcc32: close file using MPIIO
#ifdef HDF5
    call h5pclose_f(plist_id,ierr)
    call h5gclose_f(file_id,ierr)
    call h5fclose_f(filesystem_id,ierr)
    call h5close_f(ierr)
#else
    call comms_close_file(fh)
#endif
#endif

    ! cks: root proc closes output file
    if (pub_on_root) then
!for non-MPIIO close unformatted file here
#if !defined(MPI) || defined(NOMPIIO)
#ifdef HDF5
       call h5gclose_f(file_id,ierr)
       call h5fclose_f(filesystem_id,ierr)
       call h5close_f(ierr)
#else
       close(unit =input_unit)
#endif
#endif
       write(stdout,*)' done'
    endif

    ! agrecocmplx: deallocate using appropriate routine
    call data_fftbox_dealloc(fftbox_buffer)
    deallocate(fftbox_complex_shifted,stat=ierr)
    call utils_dealloc_check('restart_ngwfs_tightbox_input', &
         'fftbox_complex_shifted',ierr)
    deallocate(fftbox_complex,stat=ierr)
    call utils_dealloc_check('restart_ngwfs_tightbox_input', &
         'fftbox_complex',ierr)
#ifdef HDF5
    if (read_cmplx) then
       deallocate(buff_imag, stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_input', 'buff_imag', &
            ierr)
    end if
    deallocate(buff_real, stat=ierr)
    call utils_dealloc_check('restart_ngwfs_tightbox_input', 'buff_real', ierr)
#else
    call data_fftbox_dealloc(temp_uni_tbox)
#endif
    ! agrecocmplx: deallocate using appropriate routine
    call data_fftbox_dealloc(uni_tbox)

    ! ndmh: Clean up temporary basis if we created one, and reset coin1,2,3
    if (use_tmp_basis) then
       call function_basis_deallocate(tmp_basis)
    end if


    call timer_clock("restart_ngwfs_tightbox_input", 2)

    return

150 call utils_abort('Problem reading from file in restart_ngwfs_tightbox_&
         &input. Possible causes include: missing file, corrupted/truncated &
         &file, excessively long filename, trying to read a file created on a &
         &different CPU architecture or with a binary using different compiler &
         &settings (big-endian vs. little-endian).')

  ! gcc32: INTERNAL SUBROUTINES:
  contains

#if defined(MPI) && !defined(NOMPIIO)

    ! gcc32: for mpiio use
    subroutine internal_ngwf_read_mpiio
      use function_basis, only: function_basis_tightbox_to_ppds_basic

      implicit none

      logical :: local_ngwf
      integer :: tu_size

      if (use_tmp_basis) then
         local_ngwf = (pub_my_proc_id==tmp_basis%proc_of_func(distr_ngwf))
      else
         local_ngwf = (pub_my_proc_id==ngwf_basis%proc_of_func(distr_ngwf))
      end if

#ifdef HDF5
      if (atom_ngwf==1) then
         write(str1,*) orig_atom
         dsetname="atom_"//trim(adjustl(str1))
         call h5gopen_f(file_id,adjustl(dsetname),group_id(orig_atom),ierr, &
              gapl_id=plist_id)
      end if
      call h5dopen_f(group_id(orig_atom),'origin',dset2_id,ierr)
      write(str1,*) atom_ngwf
      ! agrecocmplx: only need real set if we are reading real NGWFs from
      ! file into complex ones
      if (read_cmplx) then
         dsetname="ngwf_"//trim(adjustl(str1))//"_real"
         call h5dopen_f(group_id(orig_atom),adjustl(dsetname),dset_id,ierr)
         dsetname="ngwf_"//trim(adjustl(str1))//"_imag"
         call h5dopen_f(group_id(orig_atom),adjustl(dsetname),dset3_id,ierr)
      else
         dsetname="ngwf_"//trim(adjustl(str1))
         call h5dopen_f(group_id(orig_atom),adjustl(dsetname),dset_id,ierr)
      end if
#else
      if (temp_uni_tbox%iscmplx) then
         ! jd: Wrapped with int() to work around gfortran's warning about
         !     assignment of an integer(8) to integer(4), due to a bug
         !     in the kind of the storage_size intrinsic (gfortran 4.9.x).
         tu_size = (int(storage_size(temp_uni_tbox%z))/8) *size(temp_uni_tbox%z)
      else
         tu_size = (int(storage_size(temp_uni_tbox%d))/8) *size(temp_uni_tbox%d)
      end if
#endif

      if (local_ngwf) then

         ! cks: initialise
         ! agrecocmplx: initialise according to uni_tbox real or complex
         call data_set_to_zero(uni_tbox)
#ifdef HDF5
         call h5pcreate_f(h5p_dataset_xfer_f,plist3_id,ierr)
         call h5pset_dxpl_mpio_f(plist3_id,h5fd_mpio_independent_f,ierr)
         adims=(/3/)
         call h5dread_f(dset2_id,h5t_native_double,read_tb_buf,adims,ierr, &
              xfer_prp=plist3_id)
         read_tb_orig1=read_tb_buf(1)
         read_tb_orig2=read_tb_buf(2)
         read_tb_orig3=read_tb_buf(3)
         ! Read real part into buff_real
         call h5dread_f(dset_id,h5t_native_double,buff_real,dims,ierr, &
              xfer_prp=plist3_id)
         if (uni_tbox%iscmplx) then
            ! agrecocmplx: if reading real NGWFs from file into complex ones,
            ! set imaginary part to zero instead of reading it from file
            if (pub_read_real_ngwfs) then
            uni_tbox%z(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                 cmplx(buff_real(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3),&
                 0.0_DP, kind=DP)
            else
               ! Read imaginary part into buff_imag
               call h5dread_f(dset3_id,h5t_native_double,buff_imag,dims,ierr, &
                    xfer_prp=plist3_id)
               uni_tbox%z(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                    cmplx(buff_real(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3), &
                    buff_imag(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3), &
                    kind=DP)
            end if
         else
            uni_tbox%d(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                 buff_real(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3)
         end if
         call h5pclose_f(plist3_id,ierr)
#else
         ! agrecocmplx: initialise according to temp_uni_tbox real or complex
         call data_set_to_zero(temp_uni_tbox)

         read_tb_buf =0.0_DP

         ! gcc32: adjust offset
         offset = offset + storage_size(baseint)/8

         call comms_read(fh,offset,read_tb_buf)

         read_tb_orig1 = read_tb_buf(1)
         read_tb_orig2 = read_tb_buf(2)
         read_tb_orig3 = read_tb_buf(3)

         ! gcc32: adjust offset
         offset = offset + (storage_size(read_tb_buf)/8)*size(read_tb_buf)+ &
              2*storage_size(baseint)/8

         ! jme: Real/complex run
         if (uni_tbox%iscmplx) then
             ! agrecocmplx: if reading real NGWFs from file into complex ones,
             ! need temporary real array to store them
             if (pub_read_real_ngwfs) then
                call comms_read(fh,offset,temp_uni_tbox%d)
                uni_tbox%z(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                     cmplx(temp_uni_tbox%d &
                               (1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3), &
                          0.0_DP, kind=DP)
             else
                call comms_read(fh,offset,temp_uni_tbox%z)
                ! gcc32: transfer temporary buffer to (possibly larger) uni_tbox
                uni_tbox%z(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                     temp_uni_tbox%z(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3)
             end if

         else
             call comms_read(fh,offset,temp_uni_tbox%d)
             uni_tbox%d(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                  temp_uni_tbox%d(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3)
         end if

         ! gcc32: adjust offset
         offset = offset + tu_size + storage_size(baseint)/8
#endif

         if (use_tmp_basis) then
            call function_basis_tightbox_to_ppds_basic(uni_tbox, &
                 read_tb_orig1,read_tb_orig2,read_tb_orig3, &
                 distr_ngwf,maxt_n1,maxt_n2,maxt_n3, ngwfs_on_grid, &
                 tmp_basis,tmp_fftbox,cell,tb_start1,tb_start2,tb_start3, &
                 fftbox_complex,fftbox_complex_shifted,fftbox_buffer)
         else
            call function_basis_tightbox_to_ppds_basic(uni_tbox, &
                 read_tb_orig1,read_tb_orig2,read_tb_orig3, &
                 distr_ngwf,maxt_n1,maxt_n2,maxt_n3, ngwfs_on_grid, &
                 ngwf_basis,tmp_fftbox,cell,tb_start1,tb_start2,tb_start3, &
                 fftbox_complex,fftbox_complex_shifted,fftbox_buffer)
         end if

      else

#ifndef HDF5
         ! gcc32: adjust offset
         offset = offset + tu_size + &
              (storage_size(read_tb_buf)/8) * size(read_tb_buf) + &
              4 * storage_size(baseint)/8
#endif

      end if
#ifdef HDF5
      call h5dclose_f(dset_id,ierr)
      ! agrecocmplx
      if (uni_tbox%iscmplx.and.(.not.pub_read_real_ngwfs)) then
         call h5dclose_f(dset3_id,ierr)
      end if
      call h5dclose_f(dset2_id,ierr)
      if (atom_ngwf==ngwf_basis%num_on_atom(regions%par%distr_atom(orig_atom))) &
           call h5gclose_f(group_id(orig_atom),ierr)
#endif

    end subroutine internal_ngwf_read_mpiio

#endif

    ! gcc32: for non-mpiio use
    subroutine internal_ngwf_read
      use function_basis, only: function_basis_tightbox_to_ppds

      implicit none

      ! cks: this essentially converts the non-blocking sends of
      ! cks: comms_send to blocking sends. Since parallel performance
      ! cks: is not the issue here, this makes the code simpler.
      call comms_barrier

      ! cks: initialise
      ! agrecocmplx: initialise according to uni_tbox real or complex
      call data_set_to_zero(uni_tbox)

      ! &&&&&&&&&&& READ INTO BUFFER &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      if (pub_on_root) then

#ifdef HDF5
         if (atom_ngwf==1) then
            write(str1,*) orig_atom
            dsetname="atom_"//trim(adjustl(str1))
            call h5gopen_f(file_id,adjustl(dsetname),group_id(orig_atom),ierr)
            call h5dopen_f(group_id(orig_atom),'origin',dset2_id,ierr)
            adims=(/3/)
            call h5dread_f(dset2_id,h5t_native_double,read_tb_buf,adims,ierr)
            read_tb_orig1=read_tb_buf(1)
            read_tb_orig2=read_tb_buf(2)
            read_tb_orig3=read_tb_buf(3)
            call h5dclose_f(dset2_id,ierr)
         end if
         write(str1,*) atom_ngwf
         if (uni_tbox%iscmplx) then
            ! agrecocmplx
            if (pub_read_real_ngwfs) then
               dsetname="ngwf_"//trim(adjustl(str1))
               call h5dopen_f(group_id(orig_atom),adjustl(dsetname),dset_id,ierr)
               call h5dread_f(dset_id,h5t_native_double,buff_real,dims,ierr)
               uni_tbox%z(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3)= &
                    cmplx(buff_real,0.0_DP,DP)
            else
               dsetname="ngwf_"//trim(adjustl(str1))//"_real"
               call h5dopen_f(group_id(orig_atom),adjustl(dsetname),dset_id,ierr)
               dsetname="ngwf_"//trim(adjustl(str1))//"_imag"
               call h5dopen_f(group_id(orig_atom),adjustl(dsetname),dset3_id,ierr)
               call h5dread_f(dset_id,h5t_native_double,buff_real,dims,ierr)
               call h5dread_f(dset3_id,h5t_native_double,buff_imag,dims,ierr)
               uni_tbox%z(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3)= &
                    cmplx(buff_real,buff_imag,DP)
            end if
         else
            dsetname="ngwf_"//trim(adjustl(str1))
            call h5dopen_f(group_id(orig_atom),adjustl(dsetname),dset_id,ierr)
            call h5dread_f(dset_id,h5t_native_double,uni_tbox%d(1:read_maxt_n1, &
                 1:read_maxt_n2,1:read_maxt_n3),dims,ierr)
         end if
         call h5dclose_f(dset_id,ierr)
         if (read_cmplx) call h5dclose_f(dset3_id,ierr)
         if (atom_ngwf==ngwf_basis%num_on_atom(regions%par%distr_atom(orig_atom))) &
              call h5gclose_f(group_id(orig_atom),ierr)
#else
         read(input_unit, err=150) read_tb_orig1, read_tb_orig2, &
              read_tb_orig3

         ! agrecocmplx: read real or complex values accordingly
         if (uni_tbox%iscmplx) then
             ! agrecocmplx
             if (pub_read_real_ngwfs) then
                read(input_unit, err=150) temp_uni_tbox%d(1:read_maxt_n1, &
                     1:read_maxt_n2, 1:read_maxt_n3)
                uni_tbox%z(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                     cmplx(temp_uni_tbox%d,0.0_DP,kind=DP)
             else
                read(input_unit, err=150) uni_tbox%z(1:read_maxt_n1, &
                     1:read_maxt_n2, 1:read_maxt_n3)
             end if
         else
             read(input_unit, err=150) uni_tbox%d(1:read_maxt_n1, &
                  1:read_maxt_n2, 1:read_maxt_n3)
         end if
#endif

      endif
      ! &&&&&&& END READ INTO BUFFER &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      if (use_tmp_basis) then

         call function_basis_tightbox_to_ppds(uni_tbox, &
              read_tb_orig1,read_tb_orig2,read_tb_orig3, &
              distr_ngwf,pub_root_proc_id,maxt_n1,maxt_n2,maxt_n3, &
              ngwfs_on_grid,tmp_basis,tmp_fftbox,cell,tb_start1,tb_start2,tb_start3, &
              fftbox_complex,fftbox_complex_shifted,fftbox_buffer, &
              sendbuf,recvbuf)

      else

         call function_basis_tightbox_to_ppds(uni_tbox, &
              read_tb_orig1,read_tb_orig2,read_tb_orig3, &
              distr_ngwf,pub_root_proc_id,maxt_n1,maxt_n2,maxt_n3, &
              ngwfs_on_grid,ngwf_basis,tmp_fftbox,cell,tb_start1,tb_start2,tb_start3, &
              fftbox_complex,fftbox_complex_shifted,fftbox_buffer, &
              sendbuf,recvbuf)

      end if

      return

150 call utils_abort('Problem reading from file in restart_ngwfs_tightbox_&
         &input/internal_ngwf_read. Possible causes include: missing file, &
         &corrupted/truncated &
         &file, excessively long filename, trying to read a file created on a &
         &different CPU architecture or with a binary using different compiler &
         &settings (big-endian vs. little-endian).')

    end subroutine internal_ngwf_read

  end subroutine restart_ngwfs_tightbox_input


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine restart_hamiltonian_read(ham,read_cond,restart_rootname)

    !======================================================================!
    ! This subroutine reads in the Hamiltonian from an external binary     !
    ! file.                                                                !
    !----------------------------------------------------------------------!
    ! Created by Alvaro Ruiz Serrano in November 2012 as a variation of    !
    ! restart_kernel_read.                                                 !
    ! Modified for embedding by Joseph Prentice, May 2018                  !
    !======================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout
    use rundat, only : pub_rootname
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_read
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! Argument
    type(SPAM3_EMBED), intent(inout) :: ham(:)
    logical, optional, intent(in) :: read_cond
    character(len=*),optional,intent(in) :: restart_rootname

    ! Local variables
    character(len=80) :: filename
    logical :: fileexists

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering restart_hamiltonian_read'

    call timer_clock("restart_hamiltonian_read", 1)

    ! KKBD: allow a custom rootname to be written
    if (present(restart_rootname)) then
       filename = restart_rootname
    else
       write(filename,'(4a)') trim(pub_rootname),'.', &
            trim(restart_read_filename_prefix),'ham'
    end if

    ! Establish filename from input file
    ! lr408: Check if we are reading in a conduction hamiltonian
    if (present(read_cond)) then
       if (read_cond) then
          write(filename,'(2a)') trim(filename),'.ham_cond'
       end if
    end if
#ifdef HDF5
    filename=trim(filename)//'.h5'
#endif

    if (pub_on_root) write(stdout,'(/3a)',advance='no') &
         'Reading Hamiltonian from file "', trim(filename),'" ...'

    ! Check that the file exists
    ! ndmh: only root proc needs to be able to see the file
    if (pub_on_root) then
       inquire(file=filename,exist=fileexists)
    else
       fileexists = .true.
    end if

    if (fileexists) then

       ! Read Hamiltonian from this file
       call sparse_embed_read(ham,trim(filename),'Hamiltonian')

    else

       call utils_abort('File '''//trim(filename)//''' not found in &
            &restart_hamiltonian_read().')

    end if

    if (pub_on_root) write(stdout,'(a)') ' done'

    call timer_clock("restart_hamiltonian_read", 2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving restart_hamiltonian_read'

  end subroutine restart_hamiltonian_read



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine restart_hamiltonian_write(ham,write_cond,write_soc,&
       file_extension)

    !======================================================================!
    ! This subroutine writes the hamiltonian matrix to an external binary  !
    ! file.                                                                !
    !======================================================================!
    ! Created by Alvaro Ruiz Serrano in November 2012 as a variation of    !
    ! restart_kernel_write.                                                !
    ! Modified for embedding by Joseph Prentice, May 2018                  !
    !======================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout, normal
    use rundat, only : pub_rootname, pub_output_detail
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_write
    use timer, only: timer_clock

    implicit none

    ! Argument
    type(SPAM3_EMBED), intent(inout) :: ham(:)
    logical, optional, intent(in) :: write_cond
    logical, optional, intent(in) :: write_soc
    character(len=*), optional, intent(in) :: file_extension

    ! Local variables
    character(len=256) :: filename
    character(len=256) :: extension

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering restart_hamiltonian_write'

    call timer_clock("restart_hamiltonian_write", 1)

    extension='ham'
    if(present(file_extension)) then
       extension=trim(adjustl(file_extension))
    end if

    ! Establish filename from input file
    write(filename,'(4a)') trim(pub_rootname),'.', &
         trim(restart_write_filename_prefix), trim(extension)

    if (present(write_soc)) then
       if (write_soc) then
          write(filename,'(2a)') trim(filename),'_soc'
       end if
    end if
    if (present(write_cond)) then
       if (write_cond) then
          write(filename,'(2a)') trim(filename),'_cond'
       end if
    end if
#ifdef HDF5
    filename=trim(filename)//'.h5'
#endif
    !if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
    !   write(stdout,'(/3a)',advance='no') &
    !        'Writing Hamiltonian to file "', trim(filename),'" ...'
    !end if

    ! Write Hamiltonian to this file
    call sparse_embed_write(ham,trim(filename),'Hamiltonian')

    !if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
    !   write(stdout,'(a,/)') ' done'
    !end if

    call timer_clock("restart_hamiltonian_write", 2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving restart_hamiltonian_write'

  end subroutine restart_hamiltonian_write


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine restart_overlap_write(ovlp,write_cond,write_inverse)

    !======================================================================!
    ! This subroutine writes the overlap matrix to an external binary      !
    ! file.                                                                !
    !======================================================================!
    ! Created by Nicholas Hine in March 2015 as a variation of             !
    ! restart_kernel_write.                                                !
    ! Modified for embedding by Joseph Prentice, May 2018                  !
    !======================================================================!


    use comms, only : pub_on_root
    use constants, only:stdout, normal
    use rundat, only : pub_rootname, pub_output_detail
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_write
    use timer, only: timer_clock

    implicit none

    ! Argument
    type(SPAM3_EMBED), intent(inout) :: ovlp
    logical, optional, intent(in) :: write_cond
    logical, optional, intent(in) :: write_inverse

    ! Local variables
    character(len=80) :: filename
    character(len=20) :: str_cond, str_inv

    call timer_clock("restart_overlap_write", 1)

    ! Check if we are writing out a conduction overlap
    str_cond = ''
    if (present(write_cond)) then
       if (write_cond) str_cond = '_cond'
    end if

    ! Check if we are writing out the inverse overlap
    str_inv = ''
    if (present(write_inverse)) then
       if (write_inverse) str_inv = 'inv_'
    end if

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering restart_overlap_write'

    ! Establish filename from input file
    write(filename,'(5a)') trim(pub_rootname),'.', &
         trim(str_inv),'ovlp',trim(str_cond)

#ifdef HDF5
    filename=trim(filename)//'.h5'
#endif
    if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
       write(stdout,'(/3a)',advance='no') &
            'Writing Overlap matrix to file "', trim(filename),'" ...'
    end if

    ! Write overlap matrix to this file
    call sparse_embed_write(ovlp,trim(filename),'Overlap matrix')

    if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
       write(stdout,'(a,/)') ' done'
    end if

    call timer_clock("restart_overlap_write", 2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving restart_overlap_write'

  end subroutine restart_overlap_write


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine restart_overlap_read(ovlp,read_cond,read_inverse)

    !======================================================================!
    ! This subroutine reads the overlap matrix from an external binary     !
    ! file.                                                                !
    !======================================================================!
    ! Created by Nicholas Hine in January 2016 as a variation of           !
    ! restart_hamiltonian_read.                                            !
    ! Modified for embedding by Joseph Prentice, May 2018                  !
    !======================================================================!


    use comms, only : pub_on_root
    use constants, only:stdout, normal
    use rundat, only : pub_rootname, pub_output_detail
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_read
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! Argument
    type(SPAM3_EMBED), intent(inout) :: ovlp
    logical, optional, intent(in) :: read_cond
    logical, optional, intent(in) :: read_inverse

    ! Local variables
    character(len=80) :: filename
    character(len=20) :: str_cond, str_inv
    logical :: fileexists

    call timer_clock("restart_overlap_read", 1)

    ! Check if we are reading the conduction overlap
    str_cond = ''
    if (present(read_cond)) then
       if (read_cond) str_cond = '_cond'
    end if

    ! Check if we are reading the inverse overlap
    str_inv = ''
    if (present(read_inverse)) then
       if (read_inverse) str_inv = 'inv_'
    end if

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering restart_overlap_read'

    ! Establish filename from input file
    write(filename,'(5a)') trim(pub_rootname),'.', &
         trim(str_inv),'ovlp',trim(str_cond)

#ifdef HDF5
    filename=trim(filename)//'.h5'
#endif
    if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
       write(stdout,'(/3a)',advance='no') &
            'Reading Overlap matrix from file "', trim(filename),'" ...'
    end if

    ! Check that the file exists
    ! ndmh: only root proc needs to be able to see the file
    if (pub_on_root) then
       inquire(file=filename,exist=fileexists)
    else
       fileexists = .true.
    end if

    if (fileexists) then

       ! Read overlap matrix from this file
       call sparse_embed_read(ovlp,trim(filename),'Overlap matrix')

    else

       call utils_abort('File '''//trim(filename)//''' not found in &
            &restart_overlap_read().')

    end if

    if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
       write(stdout,'(a,/)') ' done'
    end if

    call timer_clock("restart_overlap_read", 2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving restart_overlap_read'

  end subroutine restart_overlap_read


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine restart_kernel_read(denskern,read_cond,read_solv,read_frag, &
         read_super,read_emft,label,dkn_type,read_dmft,restart_rootname)

    !======================================================================!
    ! This subroutine reads in the density kernel from an external binary  !
    ! file.                                                                !
    ! Its operation is not affected by changes in:                         !
    !  1) internal ordering of atoms                                       !
    !  2) number of processors                                             !
    !  3) density kernel cutoff                                            !
    !  4) atomic positions                                                 !
    !----------------------------------------------------------------------!
    ! Written by Peter Haynes 15/2/2007                                    !
    ! Modified by Nicholas Hine 11/08/2008                                 !
    ! Modified by Jacek Dziedzic 19/09/2013                                !
    ! Modified by Max Phipps 21/09/2015                                    !
    ! Modified for embedding by Joseph Prentice, May 2018                  !
    ! Modified further for embedding by Joseph Prentice, Sep 2020          !
    !======================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout
    use rundat, only : pub_rootname, pub_frag_file_prefix, &
         pub_super_file_prefix, pub_frag_counter, &
         pub_eda_continuation, pub_eda_continuation_loaded_dkn, &
         pub_read_sub_denskern
    use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_read
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_assert, utils_int_to_str

    implicit none

    ! Argument
    type(SPAM3_EMBED_ARRAY), intent(inout) :: denskern
    logical, optional, intent(in) :: read_cond
    logical, optional, intent(in) :: read_solv
    logical, optional, intent(in) :: read_frag
    logical, optional, intent(in) :: read_super
    logical, optional, intent(in) :: read_emft
    integer, optional, intent(in) :: label
    character(len=*), optional, intent(in) :: dkn_type
    logical,optional :: read_dmft
    character(len=*),optional,intent(in   ) :: restart_rootname

    ! Local variables
    character(len=256) :: filename
    logical :: fileexists
    logical :: kern_conduction
    logical :: kern_solvent
    logical :: kern_frag
    logical :: kern_super


    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering restart_kernel_read'

    call timer_clock("restart_kernel_read", 1)

    ! Establish filename from input file
    if (present(restart_rootname)) then
       ! kkbd: Read in from a custom rootname
       write(filename,'(3a)') trim(restart_rootname),'.','dkn'
    else if (.not. present(label)) then
       write(filename,'(4a)') trim(pub_rootname),'.', &
            trim(restart_read_filename_prefix),'dkn'
    else
       if (present(dkn_type)) then
          write(filename,'(a256)') trim(pub_rootname)//".history.dkn"//&
             &trim(utils_int_to_str(label))//"."//trim(dkn_type)
       end if
    end if
    filename = adjustl(filename)


    kern_conduction = .false.
    kern_solvent = .false.
    kern_frag = .false.
    kern_super = .false.

    ! lr408: Check if we are reading in a conduction kernel
    if (present(read_cond)) then
       if (read_cond) kern_conduction = .true.
    end if

    if (present(read_solv)) then
       if (read_solv) kern_solvent = .true.
    end if

    if (present(read_frag)) then
       if (read_frag) kern_frag = .true.
    end if
    if (present(read_super)) then
       if (read_super) kern_super = .true.
    end if

    call utils_assert(.not. (kern_conduction .and. kern_solvent), &
         'Logic error in restart_kernel_read().')

    if (kern_conduction) then
       write(filename,'(2a)') trim(pub_rootname),'.dkn_cond'
    end if

    if (kern_solvent) then
       write(filename,'(2a)') trim(pub_rootname),'.vacuum_dkn'
    end if

    if (present(read_dmft)) then
       if (read_dmft) then
          write(filename,'(2a)') trim(pub_rootname),'.dkn_dmft'
       end if
    end if

    if (present(read_emft)) then
       if (read_emft) then
          write(filename,'(2a)') trim(filename),'_emft'
       end if
    end if

#ifdef HDF5
    filename=trim(filename)//'.h5'
#endif

    if (kern_frag) then
       ! mjsp: fragment kernel loading.
       write(filename,'(2a)') trim(pub_frag_file_prefix(pub_frag_counter)),'.dkn'
    end if
    if (kern_super) then
       ! mjsp: supermolecule kernel loading.
       if (pub_eda_continuation .and. &
            .not.(pub_eda_continuation_loaded_dkn)) then
          ! mjsp: From continuation
          write(filename,'(2a)') trim(pub_rootname),'.dkn'
          pub_eda_continuation_loaded_dkn = .true.
       else
          ! mjsp: From EDA_PREP calculation
          write(filename,'(2a)') trim(pub_super_file_prefix),'.dkn'
       end if
    end if


    if (pub_on_root .and. (denskern%mrows==1 .and. denskern%ncols==1)) then
       if(kern_conduction) then
          write(stdout,'(/3a)',advance='no') 'Reading conduction density kernel&
               & from file "', trim(filename),'" ...'
       else if(kern_solvent) then
          write(stdout,'(/a/3a)',advance='no') 'Reading density kernel for initia&
               &lising the solvent cavity from','file "', trim(filename),'" ...'
       else if(kern_frag) then
          write(stdout,'(/3a)',advance='no') 'Reading density kernel for initia&
               &lising the fragment from file "', trim(filename),'" ...'
       else if(kern_super) then
          write(stdout,'(/3a)',advance='no') 'Reading density kernel for initia&
               &lising the supermolecule from file "', trim(filename),'" ...'
       else
          write(stdout,'(/a)',advance='no') 'Reading density kernel from file '
          write(stdout,'(3a)',advance='no')'"',trim(filename),'" ...'
       end if
    end if

    ! Check that the file exists
    ! ndmh: only root proc needs to be able to see the file
    if (pub_on_root) then
       inquire(file=filename,exist=fileexists)
    else
       fileexists = .true.
    end if

    ! rc2013: don't check if this file exists if it's an embedding run
    if(denskern%mrows.gt.1 .or. denskern%ncols.gt.1) then
       call sparse_embed_read(denskern%m,trim(filename),'density kernel', &
            read_diagonal_blocks = pub_read_sub_denskern)
    else if (fileexists) then
       ! Read density kernel from this file
       !call sparse_array_read(denskern,trim(filename))
       call sparse_embed_read(denskern%m,trim(filename))
    else

       call utils_abort('File '''//trim(filename)//''' not found in &
            &restart_kernel_read().')

    end if

    if (pub_on_root) write(stdout,'(a)') ' done'

    call timer_clock("restart_kernel_read", 2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving restart_kernel_read'

  end subroutine restart_kernel_read


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine restart_RPA_kernel_read(denskern,num)

    !======================================================================!
    ! This subroutine reads in a given RPA density kernel pair p,q and     !
    ! stores it in denskern. Num specifies the number of the eigenvector   !
    ! pair that is read in.                                                !
    !----------------------------------------------------------------------!
    ! Written by Tim Zuehlsdorff in March 2015                             !
    ! Modified for embedding by Joseph Prentice, July 2018                 !
    !======================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout
    use rundat, only : pub_rootname, pub_lr_tddft_restart_from_TDA,&
        pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_read, sparse_embed_copy
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! Argument
    integer, intent(in) :: num
    type(SPAM3_EMBED), intent(inout) :: denskern(2,pub_num_spins)

    ! Local variables
    character(len=256) :: filename, loc_filename
    character(len=1) :: isub_str, jsub_str
    logical :: fileexists, loc_fileexists
    integer :: is, isub, jsub

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering restart_RPA_&
       &kernel_read'

    call timer_clock("restart_RPA_kernel_read", 1)

    ! check if restart happens from Tamm-Dancoff result
    if(.not. pub_lr_tddft_restart_from_TDA) then
       !first make sure file you are trying to read actually exists
       ! start with p kernel
       write(filename,'(2a,i0,a)') trim(pub_rootname), '_response_denskern_', &
            num, '.dkn_p'
#ifdef HDF5
       filename=trim(filename)//'.h5'
#endif

       if (pub_on_root) then
          ! jcap: check for the correct files if we have more than one
          ! region
          if (denskern(1,1)%mrows == 1 .and. &
               denskern(1,1)%ncols == 1) then
             inquire(file=filename,exist=fileexists)
          else
             fileexists=.true.
             do isub=1,denskern(1,1)%ncols
                do jsub=1,denskern(1,1)%mrows
                   write(isub_str,'(i1)') isub
                   write(jsub_str,'(i1)') jsub
                   loc_filename = trim(filename)//'_'//isub_str//jsub_str
                   inquire(file=loc_filename,exist=loc_fileexists)
                   fileexists=fileexists.and.loc_fileexists
                end do
             end do
          end if
       else
          fileexists = .true.
       end if

       if (fileexists) then
          ! Read density kernel from this file
          call sparse_embed_read(denskern(1,:),trim(filename),'density kernel')
       else
          call utils_abort('File '''//trim(filename)//''' not found in &
              &restart_RPA_kernel_read().' // &
              'If the file that ONETEP was trying to read had been &
              &written by a version of ONETEP < v5.0, this error &
              &may be due to a discrepancy in the name of the file, &
              &the actual file name having too many spaces for the &
              &current version of ONETEP; if this is the case, &
              &please manually rename the file appropriately and &
              &run ONETEP again.')
       endif

       !first make sure file you are trying to read actually exists
       ! now do q kernel
       write(filename,'(2a,i0,a)') trim(pub_rootname), '_response_denskern_', &
            num, '.dkn_q'
#ifdef HDF5
       filename=trim(filename)//'.h5'
#endif

       if (pub_on_root) then
          ! jcap: check for the correct files if we have more than one
          ! region
          if (denskern(2,1)%mrows == 1 .and. &
               denskern(2,1)%ncols == 1) then
             inquire(file=filename,exist=fileexists)
          else
             fileexists=.true.
             do isub=1,denskern(2,1)%ncols
                do jsub=1,denskern(2,1)%mrows
                   write(isub_str,'(i1)') isub
                   write(jsub_str,'(i1)') jsub
                   loc_filename = trim(filename)//'_'//isub_str//jsub_str
                   inquire(file=loc_filename,exist=loc_fileexists)
                   fileexists=fileexists.and.loc_fileexists
                end do
             end do
          end if
       else
          fileexists = .true.
       end if

       if (fileexists) then
          ! Read density kernel from this file
          call sparse_embed_read(denskern(2,:),trim(filename),'density kernel')
       else
          call utils_abort('File '''//trim(filename)//''' not found in &
               &restart_RPA_kernel_read().' // &
               'If the file that ONETEP was trying to read had been &
               &written by a version of ONETEP < v5.0, this error &
               &may be due to a discrepancy in the name of the file, &
               &the actual file name having too many spaces for the &
               &current version of ONETEP; if this is the case, &
               &please manually rename the file appropriately and &
               &run ONETEP again.')
       endif
    else
       ! restart from TDA results:
       write(filename,'(2a,i0,a)') trim(pub_rootname), '_response_denskern_', &
            num, '.dkn'
#ifdef HDF5
       filename=trim(filename)//'.h5'
#endif

       if (pub_on_root) then
          ! jcap: check for the correct files if we have more than one
          ! region
          if (denskern(1,1)%mrows == 1 .and. &
               denskern(1,1)%ncols == 1) then
             inquire(file=filename,exist=fileexists)
          else
             fileexists=.true.
             do isub=1,denskern(1,1)%ncols
                do jsub=1,denskern(1,1)%mrows
                   write(isub_str,'(i1)') isub
                   write(jsub_str,'(i1)') jsub
                   loc_filename = trim(filename)//'_'//isub_str//jsub_str
                   inquire(file=loc_filename,exist=loc_fileexists)
                   fileexists=fileexists.and.loc_fileexists
                end do
             end do
          end if
       else
          fileexists = .true.
       end if

       if(fileexists) then
          call sparse_embed_read(denskern(1,:),trim(filename),'density kernel')
       else
          call utils_abort('File '''//trim(filename)//''' not found in &
               &restart_RPA_kernel_read().' // &
               'If the file that ONETEP was trying to read had been &
               &written by a version of ONETEP < v5.0, this error &
               &may be due to a discrepancy in the name of the file, &
               &the actual file name having too many spaces for the &
               &current version of ONETEP; if this is the case, &
               &please manually rename the file appropriately and &
               &run ONETEP again.')
       endif

       ! set p component equal to q component when starting from TDA
       do is=1, pub_num_spins
          call sparse_embed_copy(denskern(2,is),denskern(1,is))
       enddo
    endif

    call timer_clock("restart_RPA_kernel_read", 2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving restart_RPA_&
        &kernel_read'


  end subroutine restart_RPA_kernel_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine restart_response_kernel_batch_read(denskern_batch, num_batch)

    !======================================================================!
    ! This subroutine reads in a batch of response density kernels and     !
    ! stores them in the appropriate data structure                        !
    ! Modified for embedding by Joseph Prentice, July 2018                 !
    !======================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout
    use rundat, only : pub_rootname, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_read
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! Argument
    integer, intent(in) :: num_batch
    type(SPAM3_EMBED), intent(inout) :: denskern_batch(num_batch,pub_num_spins)

    ! Local variables
    character(len=256) :: filename, loc_filename
    character(len=1) :: isub_str, jsub_str
    logical :: fileexists, loc_fileexists
    integer :: icount, isub, jsub

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering restart_response_&
       &kernel_batch_read'

    call timer_clock("restart_response_kernel_batch_read", 1)

    do icount=1, num_batch
       !first make sure file you are trying to read actually exists

       write(filename,'(2a,i0,a)') trim(pub_rootname), '_response_denskern_', &
            icount, '.dkn'
#ifdef HDF5
       filename=trim(filename)//'.h5'
#endif

       if (pub_on_root) then
          ! jcap: check for the correct files if we have more than one
          ! region
          if (denskern_batch(icount,1)%mrows == 1 .and. &
               denskern_batch(icount,1)%ncols == 1) then
             inquire(file=filename,exist=fileexists)
          else
             fileexists=.true.
             do isub=1,denskern_batch(icount,1)%ncols
                do jsub=1,denskern_batch(icount,1)%mrows
                   write(isub_str,'(i1)') isub
                   write(jsub_str,'(i1)') jsub
                   loc_filename = trim(filename)//'_'//isub_str//jsub_str
                   inquire(file=loc_filename,exist=loc_fileexists)
                   fileexists=fileexists.and.loc_fileexists
                end do
             end do
          end if
       else
          fileexists = .true.
       end if

       if (fileexists) then

          ! Read density kernel from this file
          call sparse_embed_read(denskern_batch(icount,:),trim(filename),'density kernel')

       else

          call utils_abort('File '''//trim(filename)//''' not found in &
               &restart_response_kernel_batch_read().' // &
               'If the file that ONETEP was trying to read had been &
               &written by a version of ONETEP < v5.0, this error &
               &may be due to a discrepancy in the name of the file, &
               &the actual file name having too many spaces for the &
               &current version of ONETEP; if this is the case, &
               &please manually rename the file appropriately and &
               &run ONETEP again.')

       endif
    enddo

    call timer_clock("restart_response_kernel_batch_read", 2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving restart_response_&
        &kernel_batch_read'

  end subroutine restart_response_kernel_batch_read


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine restart_kernel_write(denskern,write_cond,label,dkn_type,&
       write_dmft,write_emft,file_extension)

    !======================================================================!
    ! This subroutine writes the density kernel to an external binary      !
    ! file.                                                                !
    ! Its operation is not affected by changes in:                         !
    !  1) internal ordering of atoms                                       !
    !  2) number of processors                                             !
    !  3) density kernel cutoff                                            !
    !  4) atomic positions                                                 !
    !----------------------------------------------------------------------!
    ! Written by Peter Haynes 15/2/2007                                    !
    ! Modified for embedding by Joseph Prentice, May 2018                  !
    ! Modified to deal with EMFT follow calculations by Joseph Prentice,   !
    ! April 2019                                                           !
    !======================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, normal
    use rundat, only: pub_rootname, pub_output_detail
    use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_write
    use timer, only: timer_clock
    use utils, only: utils_int_to_str, utils_assert

    implicit none

    ! Argument
    type(SPAM3_EMBED_ARRAY), intent(inout) :: denskern
    logical, optional, intent(in) :: write_cond
    integer, optional, intent(in) :: label
    character(len=*), optional, intent(in) :: dkn_type
    logical, optional, intent(in) :: write_dmft
    logical, optional, intent(in) :: write_emft
    character(len=*), optional, intent(in) :: file_extension

    ! Local variables
    character(len=256) :: filename
    character(len=256) :: extension
    character(len=*), parameter :: myself = 'restart_kernel_write'

    ! -------------------------------------------------------------------------
    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering '//myself

    call timer_clock(myself, 1)

    call utils_assert(present(label) .eqv. present(dkn_type), myself//&
         ": Optional arguments 'label' and 'dkn_type' must either both be &
         &absent or both be present")

    extension='dkn'
    if(present(file_extension)) then
       extension=trim(adjustl(file_extension))
    end if

    ! Establish filename from input file
    ! lr408: Check if we are writing out a conduction kernel
    if (present(write_cond)) then
       if (write_cond) then
          write(filename,'(4a)') trim(pub_rootname),'.',trim(adjustl(extension)),'_cond'
       else
          if (present(label)) then
             write(filename,'(a256)') trim(pub_rootname)//&
                  ".history."//trim(adjustl(extension))//trim(utils_int_to_str(label))//&
                  &"."//trim(dkn_type)
          else
             write(filename,'(4a)') trim(pub_rootname),'.',&
                  trim(restart_write_filename_prefix), trim(adjustl(extension))
          end if
       end if
    else if (present(write_dmft)) then
       if (write_dmft) then
          write(filename,'(4a)') trim(pub_rootname),'.',trim(adjustl(extension)),'_dmft'
       else
          if (present(label)) then
             write(filename,'(a256)') trim(pub_rootname)//&
                  ".history."//trim(adjustl(extension))//trim(utils_int_to_str(label))//&
                  &"."//trim(dkn_type)
          else
             write(filename,'(4a)') trim(pub_rootname),'.',&
                  trim(restart_write_filename_prefix), trim(adjustl(extension))
          end if
       end if
    else
       write(filename,'(4a)') trim(pub_rootname),'.',&
            trim(restart_write_filename_prefix), trim(adjustl(extension))
    end if
    if (present(write_emft)) then
       if (write_emft) then
          write(filename,'(2a)') trim(filename),'_emft'
       end if
    end if

    filename = adjustl(filename)

#ifdef HDF5
    filename=trim(filename)//'.h5'
#endif

    ! rc2013: printing moved to sparse_embed since we could include
    ! regional indices
    !if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
    !      write(stdout,'(/3a)',advance='no') &
    !           'Writing density kernel to file "', trim(filename),'" ...'
    !end if

    ! Write density kernel to this file
    !call sparse_array_write(denskern,trim(filename))
    call sparse_embed_write(denskern%m,trim(filename),'density kernel')

    !if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
    !   write(stdout,'(a)') ' done'
    !end if

    call timer_clock("restart_kernel_write", 2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving restart_kernel_write'

  end subroutine restart_kernel_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine restart_RPA_kernel_write(denskern, num)

    !======================================================================!
    ! This subroutine writes the p and q components of a given RPA eigen-  !
    ! vector into two separate files. An optional parameter num specifies  !
    ! the the number of this RPA eigenvector pair.                         !
    !----------------------------------------------------------------------!
    ! Written by Tim Zuehlsdorff 18/3/2015                                 !
    ! Modified for embedding by Robert Charlton, October 2018.             !
    !======================================================================!

    use comms, only : pub_on_root
    use constants, only:stdout, normal
    use rundat, only : pub_rootname, pub_output_detail, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_write
    use timer, only: timer_clock

    implicit none

    ! Argument
    integer, intent(in) :: num
    type(SPAM3_EMBED), intent(inout) :: denskern(2,pub_num_spins)

    ! Local variables
    character(len=80) :: filename

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering restart_RPA_&
       &kernel_write'

    call timer_clock("restart_RPA_kernel_write", 1)

    ! first specify file name under which this density kernel is
    ! stored. Start with p component
    write(filename,'(2a,i0,a)') trim(pub_rootname), '_response_denskern_', &
         num, '.dkn_p'
#ifdef HDF5
    filename=trim(filename)//'.h5'
#endif

    if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
       write(stdout,'(/3a)',advance='no') &
            'Writing response kernel to file "', trim(filename),'" ...'
    end if

    ! Write density kernel to this file
    call sparse_embed_write(denskern(1,:),trim(filename))

    if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
       write(stdout,'(a)') ' done'
    end if

    ! first specify file name under which this density kernel is
    ! stored. Now do q component
    write(filename,'(2a,i0,a)') trim(pub_rootname), '_response_denskern_', &
         num, '.dkn_q'
#ifdef HDF5
    filename=trim(filename)//'.h5'
#endif

    if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
       write(stdout,'(/3a)',advance='no') &
            'Writing response kernel to file "', trim(filename),'" ...'
    end if

    ! Write density kernel to this file
    call sparse_embed_write(denskern(2,:),trim(filename))

    if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
       write(stdout,'(a)') ' done'
    end if

    call timer_clock("restart_RPA_kernel_write", 2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving restart_RPA_kernel_write'

  end subroutine restart_RPA_kernel_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine restart_response_kernel_batch_write(denskern_batch, num_batch)

    !======================================================================!
    ! This subroutine writes a batch of response density kernels into      !
    ! external binary files. The operation should not be affected by       !
    ! Internal ordering of atoms, number of processors, density kernel     !
    ! cutoff or atomic positions.                                          !
    !----------------------------------------------------------------------!
    ! Written by Tim Zuehlsdorff 4/3/2013                                  !
    ! Modified for embedding by Robert Charlton, October 2018.             !
    !======================================================================!

    use comms, only : pub_on_root
    use constants, only:stdout, normal
    use rundat, only : pub_rootname, pub_output_detail,pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_write
    use timer, only: timer_clock

    implicit none

    ! Argument
    integer, intent(in) :: num_batch
    type(SPAM3_EMBED), intent(inout) :: denskern_batch(num_batch,pub_num_spins)

    ! Local variables
    character(len=80) :: filename
    integer :: icount

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering restart_response_&
       &kernel_batch_write'

    call timer_clock("restart_response_kernel_batch_write", 1)

    ! loop over all response kernels in batch
    do icount=1, num_batch
       ! first specify file name under which this density kernel is
       ! stored.
       write(filename,'(2a,i0,a)') trim(pub_rootname), '_response_denskern_', &
            icount, '.dkn'
#ifdef HDF5
       filename=trim(filename)//'.h5'
#endif

       if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
          write(stdout,'(/3a)',advance='no') &
               'Writing response kernel to file "', trim(filename),'" ...'
       end if

       ! Write density kernel to this file
       call sparse_embed_write(denskern_batch(icount,:),trim(filename))

    enddo

    if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
       write(stdout,'(a)') ' done'
    end if

    call timer_clock("restart_response_kernel_batch_write", 2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving restart_kernel_write'

  end subroutine restart_response_kernel_batch_write

!============================================================================!
!============================================================================!
!============================================================================!

end module restart
