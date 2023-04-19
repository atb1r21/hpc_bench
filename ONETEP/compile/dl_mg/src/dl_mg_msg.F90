!> array with log messages
module dl_mg_msg
  implicit none

  integer, parameters :: nmax =100
  integer, save :: last_loaded=0, last_written=0
  character(len=100), save :: m(nmax)

contains

  subroutine load_msg(s)
    implicit none
    character(len=*) intent(in) :: s

    last_loaded = min(last_loaded+1, nmax)
    m(last_loaded) = s
    
  end subroutine load_msg


  subroutine write_msg
    implicit none

    do i=last_written+1, last_loaded
       write(*,'a') m(i)
    enddo

    last_written = last_loaded
  end subroutine write_msg
  
end module dl_mg_msg
