! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:!
!================================================================!
!                                                                !
!                       Anharmonic Module                        !
!                                                                !
! This module calculates the IR frequencies from ab initio MD    !
! simulations. Here the Fourier transform of the dipole moment   !
! autocorrelation function formalism is used, in order to        !
! calculate the IR-spectrum.                                     !
!----------------------------------------------------------------!
! Written by Valerio Vitale (November 2013)                      !
!================================================================!


module anharmonic

  use constants, only: DP

  implicit none

  private

type ATOM_INFO

     character(len=64) :: atom_label
     real(kind=DP)     :: atom_mass
     real(kind=DP),allocatable :: velocity_matrix(:,:)
     real(kind=DP),allocatable :: force_matrix(:,:)
     complex(kind=DP),allocatable :: power_spectrum(:)
     real(kind=DP),allocatable :: nm_velocity(:,:)

end type ATOM_INFO

! vv: Soubroutine and functions:
public :: anharmonic_main

! vv: Global variables
real(kind=DP), public, save :: cell_volume


contains

!==============================================================================!

#ifdef PGI_BUG

  ! ndmh_2014_08_15: Temporary hack to disable compilation of anharmonic
  ! module with PGI compiler so that it does not trigger the compiler bug in the
  ! current version of PGI compiler (will be removed once that is fixed)

  subroutine anharmonic_main(elements,cell)

    !===================================================================!
    ! This subroutine computes and outputs both the IR and the          !
    ! vibratonal spectrum of molecules in gas phase                     !
    ! through the Fourier transform of the dipole moment and velocities !
    ! autocorrelation function, respectively.                           !
    ! Several quantum correlation factors can be used to restore the    !
    ! asymmetry embodied in the detailed balance condition for the IR   !
    ! spectrum.                                                         !
    !-------------------------------------------------------------------!
    ! Written by Valerio Vitale (november 2013)                         !
    !===================================================================!

    use comms,           only: pub_on_root
    use constants,       only: DP
    use dense,           only: DEM
    use geometry,  only: operator(.CROSS.), operator(.DOT.)
    use ion,             only: element
    use linalg,          only: linalg_dsygv_lt
    use rundat,          only: pub_anh_type,pub_rootname,pub_print_qc
    use simulation_cell, only: CELL_INFO
    use utils,           only: utils_alloc_check, utils_dealloc_check, utils_feature_not_supported

    implicit none

    ! vv: Arguments
    type(ELEMENT), intent(in) :: elements(:)
    type(CELL_INFO), intent(in) :: cell

    call utils_feature_not_supported('You appear to be using the PGI_BUG compiler flag, indicating &
         &that you are using a bugged version of the compiler which does not work with the &
         &anharmonic module or your config file needs updating.')

   end subroutine anharmonic_main

#else

  subroutine anharmonic_main(elements,cell)

    !===================================================================!
    ! This subroutine computes and outputs both the IR and the          !
    ! vibratonal spectrum of molecules in gas phase                     !
    ! through the Fourier transform of the dipole moment and velocities !
    ! autocorrelation function, respectively.                           !
    ! Several quantum correlation factors can be used to restore the    !
    ! asymmetry embodied in the detailed balance condition for the IR   !
    ! spectrum.                                                         !
    !-------------------------------------------------------------------!
    ! Written by Valerio Vitale (november 2013)                         !
    !===================================================================!

    use bibliography,    only: bibliography_cite
    use comms,           only: pub_on_root
    use constants,       only: stdout,DP
    use ion,             only: element
    use rundat,          only: pub_anh_type,pub_rootname,pub_print_qc
    use simulation_cell, only: CELL_INFO
    use utils,           only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! vv: Arguments
    type(ELEMENT), intent(in)   :: elements(:)
    type(CELL_INFO), intent(in) :: cell

    ! vv: Local variables
    real(kind=DP),allocatable :: w_mat(:,:)
    real(kind=DP),allocatable :: T_mat(:)
    real(kind=DP)             :: avg_temp,avg_dip
    integer                   :: nat
    integer                   :: ierr
    type(ATOM_INFO),allocatable :: w_infos(:)
    integer  :: Ntot,Npad
    integer  :: itemp

    call bibliography_cite('ANHARM_AND_DMA')

    ! vv: Display Banner
    if (pub_on_root) write(stdout,1) ' Beginning of ONETEP Anharmonic '

    ! vv: Initialise workspace and allocate memory

    ! ndmh: avoid automatic allocation
    nat = size(elements)
    allocate(w_infos(nat),stat=ierr)
    call utils_alloc_check('anharmonic_main','w_infos',ierr)

    call internal_allocate_structure()
    if(pub_on_root) then
       write(stdout,'(a)') ' Initialising workspace ...'
    end if

    ! vv: Get info from the .md file
    if(pub_on_root) then
       write(stdout,'(a,a,a)') 'Opening ',trim(pub_rootname),'.md file and &
            &extracting the info ...'
    end if
    call internal_extract_info()

    ! vv: Calculate the average temperature and tot dipole moment
    if(pub_print_qc) then
    avg_temp = 0.0_DP
    avg_dip  = 0.0_DP
    do itemp = 1,Ntot
       avg_temp = avg_temp + T_mat(itemp)
       avg_dip  = avg_dip + sqrt(w_mat(itemp,1)**2 + w_mat(itemp,2)**2 +&
                   w_mat(itemp,3)**2)
    end do
    avg_temp = avg_temp/Ntot
    avg_dip  = avg_dip/Ntot

    call anharmonic_qc_output('average temperature',avg_temp)
    call anharmonic_qc_output('average dipole moment',avg_dip)

    end if

    ! vv: Select which calculation to perform
    select case(pub_anh_type)

    case('ir_calculation','IR_CALCULATION')
       call ir_spectrum_main(Ntot,Npad,w_mat,cell)

    case('vib_calculation','VIB_CALCULATION')
       call power_spectrum_main(Ntot,Npad,nat,w_infos)

    case('nm_extraction','NM_EXTRACTION')
       call nm_spectrum_main(Ntot,nat,Npad,w_infos,elements)

    end select

   ! Display Banner
   if (pub_on_root) write(stdout,1) ' End of ONETEP Anharmonic '

   deallocate(w_mat,stat=ierr)
   call utils_dealloc_check('anharmonic_main','w_mat',ierr)
   deallocate(w_infos,stat=ierr)
   call utils_dealloc_check('anharmonic_main','w_infos',ierr)

   return

1   format(/,80('x'),/,11('x'),a47,11x,11('x'),/,80('x'),/)

  contains

!=============================================================================!

    subroutine internal_allocate_structure()

      !===================================================================!
      ! This subroutine checks the input, calculates the total number     !
      ! of iterations, and allocates the memory for the workspace.        !
      ! A zero padding region is added in order to get a faster fft       !
      !===================================================================!

      use rundat,    only: pub_anh_first_iter,pub_anh_last_iter,md_num_iter,&
           pub_anh_type
      use utils,     only: utils_assert, utils_alloc_check

      implicit none

      ! Local variables
      integer :: ierr,iat


      Ntot = pub_anh_last_iter + 1 - pub_anh_first_iter

      ! vv: Check if first_iter is > 0
      call utils_assert(pub_anh_first_iter > 0, 'Error in allocate_dip_mat()&
           &The value of anh_first_iter must be greater than:', 0)
      ! vv: Check if last_iter is > 0
      call utils_assert(pub_anh_last_iter > 0, 'Error in allocate_dip_mat()&
           &The value of anh_last_iter must be greater than:', pub_anh_last_iter)
      ! vv: Check if last_iter is i<= md_num_iter
      call utils_assert(pub_anh_last_iter <= md_num_iter, 'Error in allocate_dip_mat()&
           &The value of anh_last_iter must be lower or equal to:', md_num_iter)
      ! vv: Check if first_iter is < Ntot
      call utils_assert(pub_anh_first_iter < pub_anh_last_iter, 'Error in allocate_dip_mat()&
           &The value of anh_first_iter must be lower than:', pub_anh_last_iter )
      call utils_assert(Ntot <= md_num_iter + 1, 'Error in allocate_dip_mat()&
           &The number of iterations must be smaller or equal to:', md_num_iter)

      ! vv: Find the next higher power of 2, for zero-padding the dipole moment matrix
      Npad = nextpow2(Ntot)

      allocate(T_mat(Ntot),stat=ierr)
      call utils_alloc_check('allocate_structure','T_mat',ierr)
      ! vv: Allocate the memory for the zero-padded dipole moment matrix
      allocate(w_mat(Npad,3),stat=ierr)
      call utils_alloc_check('allocate_structure','w_mat',ierr)

      select case(pub_anh_type)

      case('vib_calculation','VIB_CALCULATION')
         do iat=1,nat
            ! vv: Allocate the memory for the zero-padded velocity matrices
            allocate(w_infos(iat)%velocity_matrix(Npad,3),stat=ierr)
            call utils_alloc_check('allocate_structure','w_infos(iat)%velocity_matrix',ierr)
         end do

      case('nm_extraction','NM_EXTRACTION')
         do iat=1,nat
            ! vv: Allocate the memory for the zero-padded velocity matrices
            allocate(w_infos(iat)%velocity_matrix(Npad,3),stat=ierr)
            call utils_alloc_check('allocate_structure','w_infos(iat)%velocity_matrix',ierr)
            ! vv: Allocate the memory for the zero-padded force matrices
            allocate(w_infos(iat)%force_matrix(Npad,3),stat=ierr)
            call utils_alloc_check('allocate_structure','w_infos(iat)%force_matrix',ierr)
         end do

      end select

    end subroutine internal_allocate_structure

!==============================================================================!

    integer function  nextpow2(old_value)

      use constants, only: DP

      ! Arguments
      integer,intent(in)       :: old_value

      ! Local variables
      real(kind=DP)            :: base2_log

      ! Calculate log_2(x) as log(x)/(log(2)), where x = 2N-1
      base2_log = log(2.0_DP*real(old_value,kind=DP)-1)/log(2.0_DP)
      ! 2^next higher power of two
      nextpow2 = 2**floor(base2_log+1)

    end function nextpow2

!==============================================================================!

    subroutine internal_extract_info()

      !===================================================================!
      ! This subroutine reads the .md file and extract the total dipole   !
      ! moment and the velocities of each atom at each iteration,         !
      ! considering only the production part of the md simulation if      !
      ! specified                                                         !
      !===================================================================!

      use comms,          only: pub_on_root,pub_root_proc_id,comms_bcast
      use constants,      only: DP, stdout
      use rundat,         only: md_num_iter,pub_rootname,&
           pub_anh_first_iter,pub_anh_last_iter, pub_anh_type, pub_debug_on_root
      use utils,          only: utils_unit, utils_open_unit_check, &
           utils_close_unit_check, utils_alloc_check, &
           utils_dealloc_check, utils_read_check, utils_abort, &
           utils_sanity_check, utils_assert

      implicit none

      ! vv: Internal variables
      character(len=256) :: md_filename
      character(len=80)  :: linebuffer
      character(len=12)  :: action,form,stat,position,access
      integer            :: md_unit,md_stat
      integer            :: i,il,ios,iat,iter
      integer            :: iline1,iline2,ilinet,ilines
      logical            :: fileexists,fileopened

      ! -----------------------------------------------------------------------

      if (pub_debug_on_root) write(stdout,'(/a)') &
           'DEBUG: Entering extract_info'

      action = 'WRITE'
      form   = 'FORMATTED'
      stat   = 'OLD'
      position = 'APPEND'
      access   = 'SEQUENTIAL'

      ! vv: Open the .md file
      if (pub_on_root) then
        ! vv: Find available unit specifier
         md_unit = utils_unit()
         write(md_filename,'(a,a)') trim(pub_rootname),'.md'
         inquire(file=md_filename,exist=fileexists)
         if (fileexists) then
            open(unit=md_unit,iostat=md_stat,file=md_filename,status=stat,&
                 access=access,form=form,position=position,err=10)
             endfile (md_unit,iostat=ios)
             close(md_unit,err=20)
         else
            call utils_abort('File '''//trim(md_filename)//''' not found in &
                 &extract_info().')
         end if
      end if

      ! vv: Reopen the .md file to extract the info
      if(pub_on_root) then
11       position = 'REWIND'
         inquire(file=md_filename,opened=fileopened)
         if(.not. fileopened) then
            md_unit = utils_unit()
            open(unit=md_unit,iostat=md_stat,file=md_filename,status=stat,&
                 access=access,form=form,position=position,err=10)
         else
            close(md_unit,err=20)
            go to 11
         end if
      end if

      if(pub_on_root) then
         ilinet = 0
         ilines=0
         linebuffer = repeat(' ',80)
         do
            read(md_unit,'(a80)',iostat=ios) linebuffer
            if (is_iostat_end(ios)) exit
            ilines=ilines+1
            ! call utils_read_check('extract_info','linebuffer',ios)
            if(index(linebuffer,"<-- T")>0) then
               ilinet = ilinet + 1
               if (.not.(ilinet > pub_anh_first_iter .and. &
                    ilinet <= pub_anh_last_iter +1)) cycle
               linebuffer = adjustl(linebuffer)
               linebuffer = trim(linebuffer)
               linebuffer = linebuffer(1:index(linebuffer,' <-- T'))
               read(linebuffer,*) T_mat(ilinet-pub_anh_first_iter)
               call utils_read_check('extract_info','linebuffer',ios)
            end if
         end do
         if(all(abs(T_mat)==0.0_DP)) then
            call utils_abort('Error in extract_info() no physical data in the'&
                 //trim(pub_rootname)//'.md file.')
         end if
         call utils_assert(ilinet == md_num_iter + 1, 'Error in extract_info() &
              &the number of md_num_iter must be equal to:', ilinet )
         close(md_unit,err=20)
      end if
      call comms_bcast(pub_root_proc_id,T_mat)

      if(pub_on_root) then
12       position = 'REWIND'
         inquire(file=md_filename,opened=fileopened)
         if(.not. fileopened) then
            md_unit = utils_unit()
            open(unit=md_unit,iostat=md_stat,file=md_filename,status=stat,&
                 access=access,form=form,position=position,err=10)
         else
            close(md_unit,err=20)
            go to 12
         end if
      end if

      ! vv: Extract the total dipole moment matrix from the .md file
      if( pub_anh_type == 'ir_calculation' .or. pub_anh_type == 'IR_CALCULATION') then
         if(pub_on_root) then
            iline1 = 0
            iline2 = 0
            linebuffer = repeat(' ',80)
            ! vv: Read the .md file into linebuffer until the end of file
            do il=1,ilines
               read(md_unit,'(a80)',iostat=ios) linebuffer
               if (is_iostat_end(ios)) exit
               call utils_read_check('extract_info','linebuffer',ios)
               ! vv: If you find <-- P
               if(index(linebuffer,"<-- P")>0) then
                  iline1 = iline1 + 1
                  ! vv: and it is the total dipole moment
                  if(mod(iline1,3)==0)then
                     iline2 = iline2 + 1
                     ! vv: if a production part is specified
                     if (.not.(iline2 > pub_anh_first_iter .and. &
                          iline2 <= pub_anh_last_iter+1)) cycle
                     ! vv: Preprocessing linebuffer in order to have the right
                     ! vv: format
                     linebuffer = adjustl(linebuffer)
                     linebuffer = trim(linebuffer)
                     linebuffer = linebuffer(1:index(linebuffer,' <-- P'))
                     ! vv: Read linebuffer into w_mat
                     read(linebuffer,*) (w_mat(iline2-pub_anh_first_iter,i),i=1,3)
                     call utils_read_check('extract_info','linebuffer',ios)
                  end if
                end if
            end do
            ! vv: Check if all the elements of w_mat are not zeros
            if(all(abs(w_mat)==0.0_DP)) then
               call utils_abort('Error in extract_info() no physical data in the'&
                    //trim(pub_rootname)//'.md file.')
            end if
            ! Check if the number of actual iterations is equal to md_num_iter
            call utils_assert(iline2 == md_num_iter + 1, 'Error in extract_info() &
                 &the number of md_num_iter must be equal to:', iline2 )
            ! Zero padding
            w_mat(Ntot+1:Npad,:) = 0.0_DP
         end if
         call comms_bcast(pub_root_proc_id,w_mat)
         ! vv: Extract the total velocities matrices from the .md file

      elseif(pub_anh_type == 'vib_calculation' .or. pub_anh_type == 'VIB_CALCULATION') then
         if(pub_on_root) then
            iat=1
            iter=1
            iline2 = 0
            linebuffer=repeat(' ',80)
            ! vv: Read the .md file into linebuffer until the end of file
            do il=1,ilines
               read(md_unit,'(a80)',iostat=ios) linebuffer
!               if (is_iostat_end(ios)) exit
               call utils_read_check('extract_info','linebuffer',ios)
               ! vv: If you find <-- V
               if(index(linebuffer,'<-- V')>0) then
                  iline2=iline2 + 1
                  ! vv: count the number of iterations
                  if(iline2 > iter*nat) then
                     iter = iter + 1
                  end if
                  ! vv: if a production part is specified
                  if (.not.(iter > pub_anh_first_iter .and. &
                       iter <= pub_anh_last_iter+1)) cycle
                  ! vv: Preprocessing linebuffer in order to have the right
                  ! vv: format
                  linebuffer = adjustl(linebuffer)
                  linebuffer = trim(linebuffer)
                  ! ja: modified to avoid compiler error with gfortran 9.1
                  linebuffer = &
                       linebuffer(12:index(linebuffer,' <-- V'))//repeat(' ',11)

                  ! vv: Count the atoms
                  if(iat > nat) iat = 1
                  ! vv: Read linebuffer into velocity_matrix
                  read(linebuffer,*)(w_infos(iat)%velocity_matrix(iter-pub_anh_first_iter,&
                       i),i=1,3)
                  call utils_read_check('extract_info','linebuffer',ios)
                  ! vv: Update the number of atoms
                  iat = iat + 1
                end if
            end do
            ! vv: Zero padding
           do iat=1,nat
               w_infos(iat)%velocity_matrix(Ntot+1:Npad,:) = 0.0_DP
           end do
         end if
         ! vv: Broadcast the result
         do iat=1,nat
            call comms_bcast(pub_root_proc_id,w_infos(iat)%velocity_matrix)
         end do

      elseif(pub_anh_type == 'nm_extraction' .or. pub_anh_type == 'NM_EXTRACTION') then
         if(pub_on_root) then
            iat=1
            iter=1
            iline2 = 0
            linebuffer=repeat(' ',80)
            ! vv: Read the .md file into linebuffer until the end of file
            do il=1,ilines
               read(md_unit,'(a80)',iostat=ios) linebuffer
               !if (is_iostat_end(ios)) exit
               call utils_read_check('extract_info','linebuffer',ios)
               ! vv: If you find <-- V
               if(index(linebuffer,'<-- V')>0) then
                  iline2=iline2 + 1
                  ! vv: count the number of iterations
                  if(iline2 > iter*nat) then
                     iter = iter + 1
                  end if
                  ! vv: if a production part is specified
                  if (.not.(iter > pub_anh_first_iter .and. &
                       iter <= pub_anh_last_iter+1)) cycle
                  ! vv: Preprocessing linebuffer in order to have the right
                  ! vv: format
                  linebuffer = adjustl(linebuffer)
                  linebuffer = trim(linebuffer)
                  ! ja: modified to avoid compiler error with gfortran 9.1
                  linebuffer = &
                       linebuffer(12:index(linebuffer,' <-- V'))//repeat(' ',11)


                  ! vv: Count the atoms
                  if(iat > nat) iat = 1
                  ! vv: Read linebuffer into velocity_matrix
                  read(linebuffer,*)(w_infos(iat)%velocity_matrix(iter-pub_anh_first_iter,&
                       i),i=1,3)
                  call utils_read_check('extract_info','linebuffer',ios)
                  ! vv: Update the number of atoms
                  iat = iat + 1
               end if
             end do
           do iat=1,nat
               w_infos(iat)%velocity_matrix(Ntot+1:Npad,:) = 0.0_DP
           end do
         end if

         do iat=1,nat
            call comms_bcast(pub_root_proc_id,w_infos(iat)%velocity_matrix)
         end do

         if(pub_on_root) then
            position = 'REWIND'
               md_unit = utils_unit()
               open(unit=md_unit,iostat=md_stat,file=md_filename,status=stat,&
                 access=access,form=form,position=position,err=10)
         end if

         if(pub_on_root) then
            iat=1
            iter=1
            iline2 = 0
            linebuffer=repeat(' ',80)
            ! vv: Read the .md file into linebuffer until the end of file
            do il=1,ilines
               read(md_unit,'(a80)',iostat=ios) linebuffer
               !if (is_iostat_end(ios)) exit
               call utils_read_check('extract_info','linebuffer',ios)
               ! vv: If you find <-- V
               if(index(linebuffer,'<-- F')>0) then
                  iline2=iline2 + 1
                  ! vv: count the number of iterations
                  if(iline2 > iter*nat) then
                     iter = iter + 1
                  end if
                  ! vv: if a production part is specified
                  if (.not.(iter > pub_anh_first_iter .and. &
                       iter <= pub_anh_last_iter+1)) cycle
                  ! vv: Preprocessing linebuffer in order to have the right
                  ! vv: format
                  linebuffer = adjustl(linebuffer)
                  linebuffer = trim(linebuffer)
                  ! ja: modified to avoid compiler error with gfortran 9.1
                  linebuffer = &
                       linebuffer(12:index(linebuffer,' <-- F'))//repeat(' ',11)

                  ! vv: Count the atoms
                  if(iat > nat) iat = 1
                  ! vv: Read linebuffer into velocity_matrix
                  read(linebuffer,*)(w_infos(iat)%force_matrix(iter-pub_anh_first_iter,&
                       i),i=1,3)
                  call utils_read_check('extract_info','linebuffer',ios)
                  ! vv: Update the number of atoms
                  iat = iat + 1
               end if
            end do
            do iat=1,nat
                w_infos(iat)%force_matrix(Ntot+1:Npad,:) = 0.0_DP
            end do
         end if
         ! vv: Broadcast the result
         do iat=1,nat
            call comms_bcast(pub_root_proc_id,w_infos(iat)%force_matrix)
         end do
      endif

      if(pub_on_root)then
      close(md_unit,err=20)
      end if

      if (pub_debug_on_root) write(stdout,'(/a)') &
           'DEBUG: Leaving extract_info'

      return

10    call utils_abort('Error during opening of file:'//trim(md_filename)//&
           ' in dip_mom_into_matrix().')

20    call utils_abort('Error during closing of file:'//trim(md_filename)//&
           ' in dip_mom_into_matrix().')

    end subroutine internal_extract_info

  end subroutine anharmonic_main
!==============================================================================!

  subroutine power_spectrum_main(tot_iter,tot_iter_pad,nat,w_atom_info)

    !===================================================================!
    ! This subroutine calculates the Fourier transform of the velocities!
    ! autocorrelation function and output the results to an external    !
    ! file.                                                             !
    !-------------------------------------------------------------------!
    ! Written by Valerio Vitale (November 2013)                         !
    !===================================================================!


    use comms,     only: pub_on_root
    use constants, only: DP,stdout
    use rundat,    only: pub_anh_acf_factor,pub_rootname, pub_debug_on_root
    use utils,     only: utils_alloc_check, utils_dealloc_check


    implicit none

    ! vv: Arguments
    integer, intent(in)           :: nat
    type(ATOM_INFO), intent(inout) :: w_atom_info(nat)
    integer, intent(in)  :: tot_iter,tot_iter_pad

    ! vv: Internal variables
    integer :: iat
    real(kind=DP),allocatable    :: sample_f(:)
    complex(kind=DP),allocatable :: crosssection_fun(:)
    integer                      :: ierr

    if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering vib_spectrum_main'

    ! vv: Calculate the frequency sample array
    call anharmonic_freq_sample(sample_f,tot_iter)

    ! vv: Banner
    if (pub_on_root) then
       write(stdout,'(a)') 'Calculating the Fourier transform of the velocities'
       write(stdout,'(a)') 'by mean of the Weiner-Khinchin theorem ...'
    end if

    ! vv: Fourier transform of the velocities matrices
    do iat=1,nat
       call anharmonic_ft_dmacf(tot_iter,tot_iter_pad,w_atom_info(iat)%velocity_matrix,&
            w_atom_info(iat)%power_spectrum,pub_anh_acf_factor)
    end do

    ! vv: Calculate the sum of the autocorrelation functions
    allocate(crosssection_fun(2*tot_iter-1),stat=ierr)
    call utils_alloc_check('power_spectrum_main','crosssection_fun',ierr)

    do iat=1,nat
       crosssection_fun(:) = crosssection_fun(:) + &
            w_atom_info(iat)%power_spectrum(:)
    end do

    ! vv: Banner
    if (pub_on_root) then
       write(stdout,'(a,a,a)') 'Outputing the result into',trim(pub_rootname),'.anh ...'
    end if


    ! vv: Write info into the output file
    call anharmonic_write_freq_info(tot_iter,sample_f,crosssection_fun)

    do iat=1,nat
    deallocate(w_atom_info(iat)%velocity_matrix,stat=ierr)
    call utils_dealloc_check('power_spectrum_main','w_atom_info%velocity_matrix',ierr)
    deallocate(w_atom_info(iat)%power_spectrum,stat=ierr)
    call utils_dealloc_check('power_spectrum_main','w_atom_info%power_matrix',ierr)
    end do

    if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving vib_spectrum_main'

  end subroutine power_spectrum_main
!==============================================================================!

  subroutine ir_spectrum_main(tot_iter,tot_iter_pad,dip_mat,cell)

    !===================================================================!
    ! This subroutine computes and outputs the IR spectrum of gas phase !
    ! molecules through the Fourier transform of the dipole moment      !
    ! autocorrelation function.                                         !
    ! Several quantum correlation factors can be used to restore the    !
    ! asymmetry embodied in the detailed balance condition              !
    !-------------------------------------------------------------------!
    ! Written by Valerio Vitale (november 2013)                         !
    !===================================================================!

    use constants, only: DP,stdout
    use comms,     only: pub_on_root
    use geometry,  only: operator(.CROSS.), operator(.DOT.)
    use simulation_cell, only : CELL_INFO
    use rundat,    only: pub_anh_qc_factor, pub_anh_acf_factor, pub_rootname, &
         pub_debug_on_root
    use utils,     only: utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in)  :: tot_iter,tot_iter_pad
    real(kind=DP),allocatable, intent(in)    :: dip_mat(:,:)
    type(CELL_INFO), intent(in) :: cell

    ! Internal variables
    real(kind=DP),allocatable    :: sample_f(:)
    complex(kind=DP),allocatable :: crosssection_fun(:)

    if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering ir_spectrum_main'

    ! vv: Calculate the volume of the cell
    cell_volume = abs((cell%a1 .CROSS. cell%a2) .DOT. cell%a3)

    ! vv: Calculate the frequency sample array
    call anharmonic_freq_sample(sample_f,tot_iter)

    ! vv: Banner
    if (pub_on_root) then
       write(stdout,'(a)') 'Calculating the Fourier transform of the dipole'
       write(stdout,'(a)') 'moment matrix by mean of the Weiner-Khinchin theorem ...'
    end if

    ! vv: Fourier transform of the Dipole Moment AutoCorrelation Function-DMACF
    call anharmonic_ft_dmacf(tot_iter,tot_iter_pad,dip_mat,crosssection_fun,&
         pub_anh_acf_factor)

    ! vv: Banner
    if (pub_on_root) then
       write(stdout,'(a)') 'Applying the quantum correction factor ...'
    end if

    ! vv: Apply the quantum correction factor
    select case (pub_anh_qc_factor)

    case ('none','NONE')

       call anharmonic_qcf_none(tot_iter,sample_f,crosssection_fun)

    case ('harmonic','HARMONIC')

       call anharmonic_qcf_harmonic(tot_iter,sample_f,crosssection_fun)

    case ('standard', 'STANDARD')

       call anharmonic_qcf_standard(tot_iter,sample_f,crosssection_fun)

    case ('schofield', 'SCHOFIELD')

       call anharmonic_qcf_schofield(tot_iter,sample_f,crosssection_fun)

    end select

    ! vv: Banner
    if (pub_on_root) then
       write(stdout,'(a,a,a)') 'Outputing the result into',trim(pub_rootname),&
            '.anh ...'
    end if

    ! vv: Write info into the output file
    call anharmonic_write_freq_info(tot_iter,sample_f,crosssection_fun)

    if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving ir_spectrum_main'

    return

  end subroutine ir_spectrum_main

!==============================================================================!

  subroutine nm_spectrum_main(tot_iter,nat,tot_iter_pad,w_atom_info,elements)

    use constants,       only: DP,stdout
    use dense,           only: DEM,dense_normal_eigensolve,dense_create,&
              dense_destroy,dense_product,dense_invert
    use ion,             only: element
    use linalg,          only: linalg_dsygv_lt
    use rundat,          only: pub_anh_acf_factor, pub_debug_on_root
    use utils,           only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! vv: Arguments
     integer, intent(in)           :: nat
    type(ATOM_INFO), intent(inout) :: w_atom_info(:)
    type(ELEMENT), intent(in) :: elements(:)
    integer, intent(in) :: tot_iter,tot_iter_pad

    ! vv: Internal variables
    real(kind=DP), allocatable   :: sample_f(:)
    type(DEM)  :: eq_force_matrix,eq_momentum_matrix,&
                  trans_force_matrix,transformation_mat
    real(kind=DP),allocatable :: eigenvalues(:)
    complex(kind=DP),allocatable :: crosssection_fun(:)
    integer                      :: ierr,iat,num
    character(len=64) :: opt

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering nm_spectrum_main'

    num = 3*nat
    allocate(eigenvalues(num),stat=ierr)
    call utils_alloc_check('nm_spectrun_main','eigenvalues',ierr)

    call dense_create(eq_force_matrix,num,num,.false.)
    call dense_create(trans_force_matrix,num,num,.false.)
    call dense_create(eq_momentum_matrix,num,num,.false.)
    call dense_create(transformation_mat,num,num,.false.)

    opt ='forces'
    call anharmonic_equilibrium_matrix(tot_iter,nat,w_atom_info,&
         eq_force_matrix,elements,opt)

    opt ='momenta'
    call anharmonic_equilibrium_matrix(tot_iter,nat,w_atom_info,&
         eq_momentum_matrix,elements,opt)

!    call wrappers_cyclic_jacobi(eq_momentum_matrix,transformation_mat)

    call dense_product(trans_force_matrix,eq_force_matrix,transformation_mat)
    call dense_product(eq_force_matrix,transformation_mat,trans_force_matrix,opA='T')

    call linalg_dsygv_lt(eq_force_matrix%dmtx,eigenvalues,&
      eq_momentum_matrix%dmtx,num)

   ! vv: Change of basis
    call anharmonic_change_of_basis(tot_iter,nat,w_atom_info,eq_force_matrix,elements)

    ! vv: Calculate the power spectrum in the new basis
    do iat=1,nat
       call anharmonic_ft_dmacf(tot_iter,tot_iter_pad,w_atom_info(iat)%velocity_matrix,&
            w_atom_info(iat)%power_spectrum,pub_anh_acf_factor)
    end do

    ! vv: Calculate the sum of the autocorrelation functions
    allocate(crosssection_fun(2*tot_iter-1),stat=ierr)
    call utils_alloc_check('nm_spectrum_main','crosssection_fun',ierr)

    do iat=1,nat
       crosssection_fun(:) = crosssection_fun(:) + &
            w_atom_info(iat)%power_spectrum(:)
    end do

    call anharmonic_freq_sample(sample_f,tot_iter)

    ! vv: Write info into the output file
    call anharmonic_write_freq_info(tot_iter,sample_f,crosssection_fun)

    do iat=1,nat
      deallocate(w_atom_info(iat)%velocity_matrix,stat=ierr)
      call utils_dealloc_check('nm_spectrum_main','w_atom_info%velocity_matrix',ierr)
      deallocate(w_atom_info(iat)%force_matrix,stat=ierr)
      call utils_dealloc_check('nm_spectrum_main','w_atom_info%force_matrix',ierr)
    end do

    call dense_destroy(eq_force_matrix)
    call dense_destroy(eq_momentum_matrix)
    call dense_destroy(transformation_mat)
    call dense_destroy(trans_force_matrix)
    deallocate(eigenvalues,stat=ierr)
    call utils_dealloc_check('nm_spectrum_main','eigenvalues',ierr)

    if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving nm_spectrum_main'

  end subroutine nm_spectrum_main

!==============================================================================!

  subroutine anharmonic_ft_dmacf(Ntot,Npad,w_mat,lineshape_fun,prefactor)

    !===================================================================!
    ! This subroutine calculates the autocorrelation function of time   !
    ! series of data and then the Fourier transform of the resulting    !
    ! function. The former is calculated by mean the Weiner-Khinchin    !
    ! theorem (O(nlogn) rather than O(n^2)).                            !
    ! The autocorrelation function must be real and even, so in order   !
    ! to maintain this symmetry, we need to have a cyclic real vector   !
    !                                                                   !
    !-------------------------------------------------------------------!
    ! Weiner-Khinchin theorem in practice                               !
    ! <mu(0)*mu(t)> = F^-1[F[|mu(t)|^2]]                                !
    ! where F is the Fourier transform and F^-1 is the inverse Fourier  !
    ! transform                                                         !
    !-------------------------------------------------------------------!
    ! Written by Valerio Vitale (November 2013)                         !
    !===================================================================!

    use constants,       only : DP,stdout
    use linalg,          only : linalg_1d_fft
    use rundat,          only : pub_anh_apply_filter, pub_debug_on_root
    use utils,           only : utils_alloc_check, utils_dealloc_check, utils_sanity_check

    implicit none

    ! vv: Argument
    integer, intent(in)       :: Ntot,Npad
    real(kind=DP), intent(in) :: w_mat(Npad,3)
    complex(kind=DP), allocatable, intent(inout) :: lineshape_fun(:)
    character(10),intent(in)  :: prefactor

    ! vv: Internal variables
    integer  :: ierr,i
    ! vv: FFT arrays
    complex(kind=DP), allocatable :: complex_w_mat(:,:),ft_w_mat(:,:),&
         acorr_fun(:,:),cyc_acorr_fun(:,:),ft_lineshape(:)

    if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering anharmonic_ft_dmacf'

    ! vv: FFT arrays allocation
    allocate(ft_w_mat(Npad,3),stat=ierr)
    call utils_alloc_check('anharmonic_ft_dmacf','ft_w_mat',ierr)

    allocate(complex_w_mat(Npad,3),stat=ierr)
    call utils_alloc_check('anharmonic_ft_dmacf','complex_w_mat',ierr)

    ! vv: complexify the real matrix in input
    complex_w_mat = cmplx(w_mat,0.0_dp,kind=DP)

    ! vv: 1 dimensional FFT in x,y,z-direction
    call linalg_1d_fft('F',Npad,complex_w_mat(:,1),ft_w_mat(:,1))
    call linalg_1d_fft('F',Npad,complex_w_mat(:,2),ft_w_mat(:,2))
    call linalg_1d_fft('F',Npad,complex_w_mat(:,3),ft_w_mat(:,3))

    ! vv: Calculate the power spectrum |dip_mat(w)|^2 = dip_mat(w)* x dip_mat(w)
    ft_w_mat = conjg(ft_w_mat)*ft_w_mat
    ! vv: Scale the result for the inverse fft
    ft_w_mat = ft_w_mat/Npad

    ! vv: Calculate the autocorrelation function through the inverse Fourier transform
    allocate(acorr_fun(Npad,3),stat=ierr)
    call utils_alloc_check('ft_dmacf','acorr_fun',ierr)

    ! vv: 1 dimensional IFFT in x,y,z_direction
    call linalg_1d_fft('B',Npad,ft_w_mat(:,1),acorr_fun(:,1))
    call linalg_1d_fft('B',Npad,ft_w_mat(:,2),acorr_fun(:,2))
    call linalg_1d_fft('B',Npad,ft_w_mat(:,3),acorr_fun(:,3))

    ! vv: Strategy:
    ! vv:         1) Double the size of the a.c.f.,i.e 2N -1
    ! vv:         2) make it cyclic
    ! vv:         3) remove the extra points due to zero padding

    allocate(ft_lineshape(2*Ntot-1),stat=ierr)
    call utils_alloc_check('ft_dmacf','ft_lineshape',ierr)
    allocate(cyc_acorr_fun(2*Ntot-1,3),stat=ierr)
    call utils_alloc_check('ft_dmacf','cyc_acorr_fun',ierr)

    ! vv: Double and resize the matrix
    cyc_acorr_fun(1:Ntot-1,:) = acorr_fun(Npad - Ntot + 2: Npad,:)
    cyc_acorr_fun(Ntot :2*Ntot-1,:) = acorr_fun(1:Ntot,:)

    ! vv: <mu(0)*mu(t)> = 1/3 sum <mu_i(0)*mu_i(t)> , where i = x,y,z
    ft_lineshape = (cyc_acorr_fun(:,1) + cyc_acorr_fun(:,2) + cyc_acorr_fun(:,3))/3.0_DP

    ! vv: Apply the gaussian filter
    if(pub_anh_apply_filter) then
       do i=1,2*Ntot-1
          ft_lineshape(i) = cmplx(exp(-10.0_DP*(((-2.0_DP*Ntot-1)/2.0_DP-1+i)&
               /(2*Ntot-1))**2)*real(ft_lineshape(i),kind=DP),0.0_DP,kind=DP)
       end do
    end if

    ! vv: Force the function to be even, i.e. cyclic matrix
    call anharmonic_shift_array(ft_lineshape,2*Ntot-1,'F')

    ! vv: Select the prefactor to scale the auto correlation function
    select case (prefactor)

    case ('none','NONE')

       ! vv: Do nothing
       ft_lineshape = ft_lineshape

    case ('biased','BIASED')

       ! vv: Divide by N = total number of points
       ft_lineshape = ft_lineshape/(2*Ntot-1)

    case ('unbiased','UNBIASED')

       ! vv: Divide the i-th element of the matrix by N-|i|
       do i=1,Ntot
          ft_lineshape(i) =  ft_lineshape(i)/(Ntot + 1 - i)
       end do

       do i=1,Ntot-1
          ft_lineshape(Ntot + i) =  ft_lineshape(Ntot + i)/i
       end do

    case ('normalized','NORMALIZED')

       ! vv: Divide by <mu(0)*mu(0)>
       ft_lineshape =  ft_lineshape/ ft_lineshape(1)

    end select

    ! vv: Calculate the Fourier transform of the dipole moment autocorrelation
    ! vv: function.
    ! vv: Note: this function must be real and even
    allocate(lineshape_fun(2*Ntot-1),stat=ierr)
    call utils_alloc_check('ft_dmacf','lineshape_fun',ierr)

    ! vv: 1 dimensional FFT
    call linalg_1d_fft('F',2*Ntot-1,ft_lineshape,lineshape_fun)

    call anharmonic_shift_array(lineshape_fun,2*Ntot-1,'B')

    ! vv: Free the allocated memory
    deallocate(acorr_fun,stat=ierr)
    call utils_dealloc_check('ft_dmacf','acorr_fun',ierr)
    deallocate(cyc_acorr_fun,stat=ierr)
    call utils_dealloc_check('ft_dmacf','cyc_acorr_fun',ierr)
    deallocate(ft_lineshape,stat=ierr)
    call utils_dealloc_check('ft_dmacf','ft_lineshape',ierr)
    deallocate(ft_w_mat,stat=ierr)
    call utils_dealloc_check('ft_dmacf','ft_w_mat',ierr)
    deallocate(complex_w_mat,stat=ierr)
    call utils_dealloc_check('ft_dmacf','complex_w_mat',ierr)

    if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving anharmonic_ft_dmacf'

    return

  end subroutine anharmonic_ft_dmacf
!==============================================================================!

  subroutine anharmonic_change_of_basis(Ntot,nat,w_mat,eig_mat,elements)

     use constants,       only : DP, periodic_table_mass
     use dense,           only : DEM,dense_invert
     use ion,             only : ELEMENT
     use utils,           only : utils_alloc_check, utils_dealloc_check

     implicit none

     ! vv: Arguments
     integer, intent(in)  :: Ntot
     integer, intent(in)  :: nat
     type(ATOM_INFO), intent(inout) :: w_mat(:)
     type(DEM),intent(inout) :: eig_mat
     type(ELEMENT), intent(in) :: elements(:)

     ! vv: Local variables
     real(kind=DP), allocatable :: qdot_vec(:)
     real(kind=DP), allocatable :: p_vec(:)
     integer :: ierr,iat,icol,jcol,iter,inat
     integer :: atom_I
     real(kind=DP) :: atom_mass_I
!     real(kind=dp), parameter :: electron_mass_si = 9.1093826e-31_dp
!     real(kind=dp), parameter :: avogadro_si = 6.0221415e23_dp

     allocate(p_vec(3*nat),stat=ierr)
     call utils_alloc_check('anharmonic_change_of_basis','p_vec',ierr)

     allocate(qdot_vec(3*nat),stat=ierr)
     call utils_alloc_check('anharmonic_change_of_basis','qdot_vec',ierr)

     do iter=1,Ntot
        iat = 1
        icol = 1
        ! vv: update p_vec
        do inat = 1,3*nat
           if(icol == 4) then
              icol = 1
              iat = iat + 1
           end if
           atom_I = elements(iat)%atomic_number
           atom_mass_I = periodic_table_mass(atom_I)!*1e-3_dp/&
           ! avogadro_si/electron_mass_si
           p_vec(inat) = atom_mass_I*w_mat(iat)%velocity_matrix(iter,icol)
           icol = icol + 1
        end do

        eig_mat%dmtx = transpose(eig_mat%dmtx)

        ! vv: calculate the qdot_vec
        do jcol =1,3*nat
           do icol = 1,3*nat
              qdot_vec(icol) = p_vec(icol)*eig_mat%dmtx(icol,jcol)
           end do
        end do

        icol = 1
        iat = 1
        do inat = 1,3*nat
           if(icol==4) then
             icol = 1
             iat = iat + 1
           end if
           w_mat(iat)%velocity_matrix(iter,icol) = qdot_vec(inat)
           icol = icol + 1
        end do
     end do

     deallocate(p_vec,stat=ierr)
     call utils_dealloc_check('anharmonic_change_of_basis','p_vec',ierr)
     deallocate(qdot_vec,stat=ierr)
     call utils_dealloc_check('anharmonic_change_of_basis','qdot_vec',ierr)

 end subroutine anharmonic_change_of_basis
!==============================================================================!

  subroutine anharmonic_equilibrium_matrix(Ntot,nat,w_mat,eq_dense_matrix,elements,parm)

     use comms,           only : pub_on_root
     use constants,       only : DP,stdout,periodic_table_mass
     use dense,           only : DEM,dense_create
     use ion,             only : element
     use rundat,          only : pub_debug_on_root
     use utils,           only : utils_alloc_check, utils_dealloc_check, utils_sanity_check, utils_abort


     implicit none

     ! vv: Argument
     integer, intent(in)           :: Ntot
     integer, intent(in)           :: nat
     type(ATOM_INFO), intent(inout)   :: w_mat(nat)
     type(DEM), intent(inout)     :: eq_dense_matrix
     type(ELEMENT), intent(in) :: elements(nat)
     character(len=64), intent(in) :: parm

     ! vv: Internal variables
     real(kind=DP),allocatable  :: tmp(:)
     real(kind=DP),allocatable  :: eq_matrix(:,:)
     integer  :: ierr,iter,iat1,iat2,ifm,kfm,icount1,icount2
     integer  :: atom_I,atom_J
     real(kind=DP)  :: atom_mass_I, atom_mass_J
!     real(kind=dp), parameter :: electron_mass_si = 9.1093826e-31_dp
!     real(kind=dp), parameter :: avogadro_si = 6.0221415e23_dp

     if (pub_debug_on_root) write(stdout,'(/a)') &
          'DEBUG: Entering anharmonic_equilibrium_matrix'

     if (eq_dense_matrix%iscmplx) then
        call utils_abort('Error in anharmonic_equilibrium_matrix:&
              &eq_dense_matrix must be real.')
     endif

     allocate(tmp(3*nat),stat=ierr)
     call utils_alloc_check('anharmonic_equilibrium_matrix','tmp',ierr)
     allocate(eq_matrix(3*nat,3*nat),stat=ierr)
     call utils_alloc_check('anharmonic_equilibrium_matrix','eq_matrix',ierr)

     eq_dense_matrix%dmtx(:,:) = 0.0_DP
     eq_matrix = eq_dense_matrix%dmtx

     if(pub_on_root) then
     tmp(:) = 0.0_DP
     icount1 = 0
     icount2 = 0
     do iat1 = 1,nat
      atom_I = elements(iat1)%atomic_number
      atom_mass_I = periodic_table_mass(atom_I)!*1e-3_dp/avogadro_si/electron_mass_si
      do ifm = 1,3
       if(parm == 'forces') then
         icount1 = icount1 + 1
          do iat2 = 1,nat
            do kfm = 1,3
               icount2 = icount2 + 1
                do iter = 1,Ntot
                tmp(icount2) = tmp(icount2) + w_mat(iat1)%force_matrix(iter,ifm)*&
                  w_mat(iat2)%force_matrix(iter,kfm)/Ntot
                end do
            end do
          end do
              eq_matrix(icount1,:) = &
              tmp(:)
          tmp(:) = 0.0_DP
          icount2 = 0
       elseif(parm == 'momenta') then
         icount1 = icount1 + 1
          do iat2 = 1,nat
          atom_J = elements(iat2)%atomic_number
          atom_mass_J = periodic_table_mass(atom_J)!*1e-3_dp/avogadro_si/electron_mass_si
            do kfm = 1,3
               icount2 = icount2 + 1
                  do iter = 1,Ntot
                  tmp(icount2) =  tmp(icount2) + w_mat(iat1)%velocity_matrix(iter,ifm)*&
                  atom_mass_I*w_mat(iat2)%velocity_matrix(iter,kfm)*&
                  atom_mass_J/Ntot
                  end do
            end do
          end do
              eq_matrix(icount1,:) = &
              tmp(:)
          tmp(:) = 0.0_DP
          icount2 = 0
       end if
      end do
     end do

  end if

  eq_dense_matrix%dmtx = eq_matrix

     deallocate(tmp,stat=ierr)
     call utils_dealloc_check('anharmonic_equilibrium_matrix','tmp',ierr)

     if (pub_debug_on_root) write(stdout,'(/a)') &
          'DEBUG: Leaving anharmonic_equilibrium_matrix'

     return

  end subroutine anharmonic_equilibrium_matrix
!==============================================================================!

  subroutine anharmonic_freq_sample(f_sample,Ntot)

    !===================================================================!
    ! This subroutine calculates the time and frequency arrays          !
    !===================================================================!

    use constants,only: DP, TWO_PI,LONG
    use rundat,   only: md_delta_t
    use utils,    only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! vv: Arguments
    real(kind=DP),intent(inout),allocatable  :: f_sample(:)
    integer,intent(in)                     :: Ntot

    ! vv: Internal variables
    integer                    :: ierr
    real(kind=DP), parameter   :: fs2s = 1.0E-15 ! atomic unit of time --> fs
    real(kind=DP), parameter   :: aut2fs = 0.0241888468_DP ! atomic unit of time --> fs

    ! vv: Time parameters
    integer(kind=LONG)         :: i,max_time
    real(kind=DP),allocatable  :: t_sample(:)
    real(kind=DP)              :: delta_t
    ! vv: Frequency parameters
    real(kind=DP)              :: delta_f,max_f,min_f

    ! vv: Allocate the memory for the time and frequency sample arrays
    allocate(t_sample(2*Ntot-1),stat=ierr)
    call utils_alloc_check('freq_sample','t_sample',ierr)
    allocate(f_sample(2*Ntot-1),stat=ierr)
    call utils_alloc_check('freq_sample','f_sample',ierr)

    ! vv: Preparing the axes
    ! vv: Extreme values in the time domain
    max_time =  2*Ntot - 1
    ! min_time = -max_time
    ! vv: time step
    delta_t = md_delta_t*aut2fs*fs2s
    ! vv: sampling time vector
    ! t_sample = (/((min_time*delta_t/2.0_DP-delta_t)+(i*delta_t),i=1,max_time )/)
    ! vv: Extreme values in the frequency domain
    max_f = TWO_PI/delta_t
    min_f = -max_f
    ! vv: frequency step
    delta_f = max_f/max_time
    ! vv: sampling frequency vector
    f_sample = (/((min_f/2.0_DP-delta_f)+(i*delta_f),i=1,max_time)/)
    ! vv: Free the occupied space
    deallocate(t_sample,stat=ierr)
    call utils_dealloc_check('freq_sample','t_sample',ierr)

  end subroutine anharmonic_freq_sample

!==============================================================================!

  subroutine anharmonic_shift_array(array,array_dim,value)

    !===================================================================!
    ! This subroutine swaps the left and the right halves of an array   !
    ! The backward function 'B' undoes the effects of the forward 'F'   !
    !=============================    ======================================!

    use constants, only: DP
    use utils,     only: utils_alloc_check

    implicit none

    ! vv: Arguments
    character,intent(in)            :: value
    integer,intent(in)              :: array_dim
    complex(kind=DP), intent(inout) :: array(array_dim)

    ! vv: Local variables
    integer                         :: new_dim
    integer                         :: ierr
    complex(kind=DP), allocatable   :: temp(:)

    new_dim = floor((array_dim+1)/2.0_DP)

    allocate(temp(array_dim),stat=ierr)
    call utils_alloc_check('shift_array','temp',ierr)

    ! vv: ifftshift
    if (value=='f' .or. value=='F') then
       temp(1:new_dim) = array(new_dim:array_dim)
       temp(new_dim+1:array_dim) = array(1:new_dim-1)
       array(:) = temp(:)
    end if

    ! vv: fftshift
    if (value=='b' .or. value=='B') then
       temp(1:new_dim-1) = array(new_dim+1:array_dim)
       temp(new_dim:array_dim) = array(1:new_dim)
       array(:) = temp(:)
    end if

  end subroutine anharmonic_shift_array

!==============================================================================!

  subroutine anharmonic_qcf_none(Ntot,f_sample,lineshape_fun)

    !===================================================================!
    ! The following four subroutines apply the quantum correction (QC)  !
    ! factors to the absorption cross-section                           !
    !                                                                   !
    ! | QC label   |     QC factor    |      List of subroutines:       !
    ! |------------------------------ |                                 !
    ! | None       |        1         |      anharmonic_qcf_none        !
    ! |            |                  |                                 !
    ! | Standard   | 2/(1+exp(-hbw))  |      anharmonic_qcf_standard    !
    ! |            |                  |                                 !
    ! | Harmonic   | bhw/(1-exp(-hbw))|      anharmonic_qcf_harmonic    !
    ! |            |                  |                                 !
    ! | Schofield  | exp(bhw/2)       |      anharmonic_qcf_schofield   !
    !===================================================================!

    use constants, only: TWO_PI,DP,k_B,FINE_STRUCTURE
    use rundat,    only: pub_anh_md_temp,md_delta_t

    implicit none

    ! vv: Arguments
    integer,intent(in)             :: Ntot
    real(kind=DP),intent(inout)       :: f_sample(2*Ntot-1)
    complex(kind=DP),intent(inout) :: lineshape_fun(2*Ntot-1)

    ! vv: Local variables
    real(kind=DP), parameter :: fs2s = 1.0E-15 ! vv: atomic unit of time --> fs
    real(kind=DP), parameter :: inv_fs2aut = 1.0_DP/0.0241888468_DP
    real(kind=DP), parameter :: au_hbar= 1.0_DP
    real(kind=DP) :: zwork(2*Ntot-1)
    real(kind=DP) :: none_factor(2*Ntot-1)

    ! vv: Calculate the absorption cross-section without any quantum correction
    ! vv: factor

    ! vv: Frequency array in a.u.
    f_sample = f_sample*fs2s*inv_fs2aut

    none_factor = FINE_STRUCTURE*TWO_PI*f_sample*&
         (1.0_DP-exp(-au_hbar*f_sample/(k_B*pub_anh_md_temp)))&
         /(3.0_DP*cell_volume*au_hbar)
    ! vv: absorption cross-section function in a.u.
    zwork = none_factor*real(lineshape_fun,kind=DP)*md_delta_t

    lineshape_fun = cmplx(zwork,0.0_DP,kind=DP)

  end subroutine anharmonic_qcf_none

!==============================================================================!

  subroutine anharmonic_qcf_harmonic(Ntot,f_sample,lineshape_fun)

    use constants, only: TWO_PI,DP,k_B
    use rundat,    only: pub_anh_md_temp,md_delta_t

    implicit none

    ! vv: Arguments
    integer,intent(in)             :: Ntot
    real(kind=DP),intent(inout)    :: f_sample(2*Ntot-1)
    complex(kind=DP),intent(inout) :: lineshape_fun(2*Ntot-1)

    ! vv: Local variables
    real(kind=DP), parameter :: fs2s = 1.0E-15 ! atomic unit of time --> fs
    real(kind=DP), parameter :: inv_fs2aut = 1.0_DP/0.0241888468_DP
    real(kind=DP) :: zwork(2*Ntot-1)
    real(kind=DP) :: harmonic_factor(2*Ntot-1)

    ! vv: Frequency array in a.u.
    f_sample = f_sample*fs2s*inv_fs2aut
    ! vv: Calculate the absorption cross-section with the harmonic quantum
    ! vv: correction factor, i.e. Bhw/(1 - exp(-Bhw)
    harmonic_factor = TWO_PI*f_sample**2/(3.0_DP*cell_volume*&
         k_B*pub_anh_md_temp)
    ! vv: absorption cross-section in a.u.
    zwork = harmonic_factor*real(lineshape_fun,kind=DP)*md_delta_t

    lineshape_fun = cmplx(zwork,0.0_DP,kind=DP)

  end subroutine anharmonic_qcf_harmonic

!==============================================================================!

  subroutine anharmonic_qcf_standard(Ntot,f_sample,lineshape_fun)

    use constants, only: TWO_PI,DP,k_B,FINE_STRUCTURE
    use rundat,    only: pub_anh_md_temp,md_delta_t

    implicit none

    ! vv: Arguments
    integer,intent(in)             :: Ntot
    real(kind=DP), intent(inout)   :: f_sample(2*Ntot-1)
    complex(kind=DP),intent(inout) :: lineshape_fun(2*Ntot-1)

    ! vv: Local variables
    real(kind=DP), parameter :: fs2s = 1.0E-15 ! atomic unit of time --> fs
    real(kind=DP), parameter :: inv_fs2aut = 1.0_DP/0.0241888468_DP
    real(kind=DP), parameter :: au_hbar= 1.0_DP
    real(kind=DP) :: zwork(2*Ntot-1)
    real(kind=DP) :: standard_factor(2*Ntot-1)


    ! vv: Frequency array in a.u.
    f_sample = f_sample*fs2s*inv_fs2aut
    ! vv: Calculate the absorption cross-section with the standard quantum
    ! vv: correction factor, i.e. 2/(1 + exp(-Bhw))
    standard_factor = FINE_STRUCTURE*2*TWO_PI*f_sample*&
         (1.0_DP-exp(-au_hbar*f_sample/(k_B*pub_anh_md_temp)))&
         /(3.0_DP*cell_volume*au_hbar*(1.0_DP + &
         exp(-au_hbar*f_sample/(k_B*pub_anh_md_temp))))
    ! vv: absorption cross-section in a.u.
    zwork = standard_factor*real(lineshape_fun,kind=DP)*md_delta_t

    lineshape_fun = cmplx(zwork,0.0_DP,kind=DP)

  end subroutine anharmonic_qcf_standard

!==============================================================================!

  subroutine anharmonic_qcf_schofield(Ntot,f_sample,lineshape_fun)

    use constants, only: TWO_PI,DP,k_B,FINE_STRUCTURE
    use rundat,    only: pub_anh_md_temp,md_delta_t

    implicit none

    ! vv: Arguments
    integer,intent(in)             :: Ntot
    real(kind=DP),intent(inout)    :: f_sample(2*Ntot-1)
    complex(kind=DP),intent(inout) :: lineshape_fun(2*Ntot-1)

    ! vv: Local variables
    real(kind=DP), parameter :: fs2s = 1.0E-15 ! atomic unit of time --> fs
    real(kind=DP), parameter :: inv_fs2aut = 1.0_DP/0.0241888468_DP
    real(kind=DP), parameter :: au_hbar= 1.0_DP
    real(kind=DP) :: zwork(2*Ntot-1)
    real(kind=DP) :: schofield_factor(2*Ntot-1)
    real(kind=DP) :: max_f_value

    ! vv: frequency in a.u.
    f_sample = f_sample*fs2s*inv_fs2aut
    max_f_value = maxval(f_sample)
    f_sample = f_sample/max_f_value
    ! vv: Calculate the absorption cross-section with the schofield quantum
    ! vv: correction factor, i.e. exp[bhw/2]
    schofield_factor = FINE_STRUCTURE*TWO_PI*f_sample*&
         (1.0_DP-exp(-au_hbar*f_sample/(k_B*pub_anh_md_temp)))&
         /(3.0_DP*cell_volume*au_hbar)*exp(au_hbar*f_sample/&
         (2.0_DP*k_B*pub_anh_md_temp))
    zwork = schofield_factor*real(lineshape_fun,kind=DP)*md_delta_t

    lineshape_fun = cmplx(zwork,0.0_DP,kind=DP)

    f_sample = f_sample*max_f_value

  end subroutine anharmonic_qcf_schofield

!==============================================================================!

  subroutine anharmonic_write_freq_info(Ntot,f_sample,lineshape_fun)

    !===================================================================!
    ! Outputs the result into the .anh_ir and anh_vib file              !
    !===================================================================!

    use comms,     only: pub_on_root
    use constants, only: DP,TWO_PI,SPEED_OF_LIGHT_SI
    use rundat,    only: pub_anh_plot_firstfreq,pub_anh_plot_lastfreq,&
         pub_rootname, pub_anh_plot_all, pub_anh_type
    use utils,     only: utils_abort, utils_unit, utils_assert

    implicit none

    ! vv: Arguments
    integer,intent(in)           :: Ntot
    real(kind=DP),intent(inout)  :: f_sample(2*Ntot-1)
    complex(kind=DP),intent(in)  :: lineshape_fun(2*Ntot-1)

    ! vv: Local variables
    real(kind=DP), parameter :: fs2s = 1.0E-15
    real(kind=DP), parameter   :: aut2fs = 0.0241888468_DP ! vv: atomic unit of time --> fs
    real(kind=DP), parameter :: inv_fs2aut = 1.0_DP/aut2fs

    integer :: i
    character(len=256) :: output_name
    character(len=12)  :: action,form,stat,position,access
    integer            :: output_unit,out_stat
    logical            :: fileexists

    action = 'WRITE'
    form   = 'FORMATTED'
    stat   = 'NEW'
    position = 'REWIND'
    access   = 'SEQUENTIAL'

    ! vv: Check if anh_plot_firstfreq is > anh_first_iter
    call utils_assert(pub_anh_plot_firstfreq >= 0, 'Error in write_freq_info(). &
         &The value of anh_plot_firstfreq must be greater or equal to:', 0)

    select case (pub_anh_type)

    case ('ir_calculation','IR_CALCULATION')
       ! vv: From si units to cm^-1
       f_sample = f_sample*aut2fs/(TWO_PI*SPEED_OF_LIGHT_SI*100.0_DP*fs2s)
       ! vv: Check if plot_firstfreq is < max_fsample
       call utils_assert(pub_anh_plot_firstfreq <= maxval(f_sample), &
            'Error in write_freq_info(). The value of pub_anh_plot_lastfreq must be &
            &lower or equal to:',maxval(f_sample))
       if(pub_on_root) then
          write(output_name,'(a,a)') trim(pub_rootname),'.anh_ir'
          inquire(file=output_name,exist=fileexists)
          if(fileexists) then
             stat='OLD'
          end if
          ! vv: Open the .anh_ir file
          output_unit= utils_unit()
          open(unit=output_unit,iostat=out_stat,file=output_name,status=stat,&
               access=access,form=form,position=position,action=action,err=10)
          write(output_unit,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
               &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
          write(output_unit,'(a)') '  IR frequencies(cm^-1)     &
               &           absorption cross-section(arbitrary units)  '
          do i=1,2*Ntot-1
             if(pub_anh_plot_all) then
                write(output_unit,*)&
                     f_sample(i),real(lineshape_fun(i),kind=DP)
             elseif(f_sample(i) >= pub_anh_plot_firstfreq .and. &
                  f_sample(i) <= pub_anh_plot_lastfreq) then
                write(output_unit,'(f16.8,21x,e19.8)')&
                     f_sample(i),real(lineshape_fun(i),kind=DP)
             end if
          end do
          write(output_unit,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
               &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
          close(output_unit,err=20)
       end if

    case('vib_calculation','VIB_CALCULATION')
       f_sample = f_sample*fs2s*inv_fs2aut
       f_sample = f_sample*aut2fs/(TWO_PI*SPEED_OF_LIGHT_SI*100.0_DP*fs2s)
       ! vv: Check if plot_firstfreq is < max_fsample
       call utils_assert(pub_anh_plot_firstfreq <= maxval(f_sample), &
            'Error in write_freq_info(). The value of pub_anh_plot_lastfreq must be &
            &lower or equal to:',maxval(f_sample))
       if(pub_on_root) then
          write(output_name,'(a,a)') trim(pub_rootname),'.anh_vib'
          inquire(file=output_name,exist=fileexists)
          if(fileexists) then
             stat='OLD'
          end if
          ! vv: Open the .anh_vib file
          output_unit= utils_unit()
          open(unit=output_unit,iostat=out_stat,file=output_name,status=stat,&
               access=access,form=form,position=position,action=action,err=10)
          write(output_unit,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
               &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
          write(output_unit,'(a)') '  frequencies(cm^-1)     &
               &           power spectrum(arbitrary units)  '
          do i=1,2*Ntot-1
             if(pub_anh_plot_all) then
                write(output_unit,*)&
                     f_sample(i),real(lineshape_fun(i),kind=DP)
             elseif(f_sample(i) >= pub_anh_plot_firstfreq .and. &
                  f_sample(i) <= pub_anh_plot_lastfreq) then
                write(output_unit,'(f16.8,21x,e19.8)')&
                     f_sample(i),real(lineshape_fun(i),kind=DP)
             end if
          end do
          write(output_unit,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
               &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
          close(output_unit,err=20)
       end if

    case('nm_extraction','NM_EXTRACTION')
       f_sample = f_sample*fs2s*inv_fs2aut
       f_sample = f_sample*aut2fs/(TWO_PI*SPEED_OF_LIGHT_SI*100.0_DP*fs2s)
       ! vv: Check if plot_firstfreq is < max_fsample
       call utils_assert(pub_anh_plot_firstfreq <= maxval(f_sample), &
            'Error in write_freq_info(). The value of pub_anh_plot_lastfreq must be &
            &lower or equal to:',maxval(f_sample))
       if(pub_on_root) then
          write(output_name,'(a,a)') trim(pub_rootname),'.anh_nme'
          inquire(file=output_name,exist=fileexists)
          if(fileexists) then
             stat='OLD'
          end if
          ! vv: Open the .anh_vib file
          output_unit= utils_unit()
          open(unit=output_unit,iostat=out_stat,file=output_name,status=stat,&
               access=access,form=form,position=position,action=action,err=10)
          write(output_unit,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
               &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
          write(output_unit,'(a)') '  frequencies(cm^-1)     &
               &           localised power spectrum(arbitrary units)  '
          do i=1,2*Ntot-1
             if(pub_anh_plot_all) then
                write(output_unit,*)&
                     f_sample(i),real(lineshape_fun(i),kind=DP)
             elseif(f_sample(i) >= pub_anh_plot_firstfreq .and. &
                  f_sample(i) <= pub_anh_plot_lastfreq) then
                write(output_unit,'(f16.8,21x,e19.8)')&
                     f_sample(i),real(lineshape_fun(i),kind=DP)
             end if
          end do
          write(output_unit,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
               &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
          close(output_unit,err=20)
       end if

    end select

    return


10  call utils_abort('Error during opening of file:'//trim(output_name)//&
         ' in write_freq_info().')

20  call utils_abort('Error during closing of file:'//trim(output_name)//&
         ' in write_freq_info().')

  end subroutine anharmonic_write_freq_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine anharmonic_qc_output(name_data,input_data)

    !===========================================================================!
    ! Prints out quality-control lines <QC> for qc-testing the anharmonic module!
    !---------------------------------------------------------------------------!
    ! Written by Valerio Vitale, 18/12/2013.                                    !
    !===========================================================================!

    use comms, only: pub_on_root
    use constants, only: DP
    use utils, only: utils_qc_print

    implicit none

    ! vv: Arguments
    real(kind=DP), intent(in) :: input_data
    character(len=*),intent(in) :: name_data

    if(pub_on_root) then
             call utils_qc_print(name_data,input_data)
    end if

  end subroutine anharmonic_qc_output

#endif

end module anharmonic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
