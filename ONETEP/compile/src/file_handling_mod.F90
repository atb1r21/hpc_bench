! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                           File Handling Module                              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   Module that allows execution of several UNIX commands within the scope    !
!   of Fortran                                                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   The subroutines in this file were written by Cedric Weber and Edward      !
!   Linscott                                                                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module file_handling

  use utils, only: utils_unit

  implicit none

  private

  public :: file_handling_remove_file
  public :: file_handling_remove_numbered_files
  public :: file_handling_copy_file
  public :: file_handling_move_file

  contains

  subroutine file_handling_remove_file(filename)

    implicit none

    character(*), intent(in)         :: filename
    integer                          :: file_id

    file_id = utils_unit()
    open(unit=file_id, file=filename)
    close(file_id, status='delete')

  end subroutine

  subroutine file_handling_remove_numbered_files(filename)

    implicit none

    character(*), intent(in)         :: filename
    integer                          :: i
    logical                          :: file_exists
    character(len=4)                 :: istr

    i = 1
    file_exists = .true.
    do while (file_exists)
       write(istr,'(i4)') i
       inquire(file=filename//trim(adjustl(istr)), exist=file_exists)
       if (file_exists) call file_handling_remove_file(filename//trim(adjustl(istr)))
    end do

  end subroutine

  subroutine file_handling_copy_file(filename1, filename2, formatted)

  implicit none

  character(*), intent(in) :: filename1, filename2
  integer                  :: unit1, unit2, ierr1, ierr2, irec
  character                :: ch_unform
  character(132)           :: line_form
  logical, optional        :: formatted


  unit1 = utils_unit()
  unit2 = utils_unit()

  if (.not. present(formatted)) formatted = .false.

  if (formatted) then

     open(unit=unit1, file=trim(adjustl(filename1)), action='read', status='old', &
          iostat=ierr1, access='sequential')
     open(unit=unit2, file=trim(adjustl(filename2)), action='write', status='replace', &
          iostat=ierr2, access='sequential')

     do
        read(unit1, '(a)', iostat=ierr1) line_form
        if (ierr1 /= 0) exit
        write(unit2, '(a)') trim(line_form)
     end do

  else

     open(unit=unit1, file=trim(adjustl(filename1)), action='read', status='old', &
          recl=1, iostat=ierr1, access='direct')
     open(unit=unit2, file=trim(adjustl(filename2)), action='write', status='replace', &
          recl=1, iostat=ierr2, access='direct')

     irec = 1
     do
        read(unit=unit1, rec=irec, iostat=ierr1) ch_unform
        if (ierr1 /= 0) exit
        write(unit=unit2, rec=irec) ch_unform
        irec = irec+1
     end do

  end if

  close(unit1)
  close(unit2)

  end subroutine

  subroutine file_handling_move_file(filename1, filename2, formatted)

  implicit none

  character(*), intent(in) :: filename1, filename2
  logical, optional        :: formatted

  if (.not. present(formatted)) formatted = .false.

  call file_handling_copy_file(filename1, filename2, formatted)
  call file_handling_remove_file(filename1)

  end subroutine
end module
