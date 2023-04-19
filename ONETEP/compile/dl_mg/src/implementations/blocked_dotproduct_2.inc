  subroutine loops_ab(a,b,xloc)
    implicit none
    real(wp), intent(in)           :: a(ba(1):,ba(2):,ba(3):)
    real(wp), intent(in), optional :: b(bb(1):,bb(2):,bb(3):)
    real(wp), intent(inout)        :: xloc

    integer, pointer :: s(:), e(:)
    integer iblk

    if (present(b)) then
       !$OMP DO SCHEDULE(STATIC,1) reduction(+:xloc)
       do iblk = 1, size(blk(t)%start, dim=2)
          s => blk(t)%start(:,iblk)
          e => blk(t)%end(:, iblk)
          xloc = xloc + scl * sum(a(s(1):e(1),s(2):e(2),s(3):e(3)) &
               * b(s(1):e(1),s(2):e(2),s(3):e(3)))
       enddo
       !$OMP ENDDO
    else
       !$OMP DO SCHEDULE(STATIC,1) reduction(+:xloc)
       do iblk = 1, size(blk(t)%start, dim=2)
          s => blk(t)%start(:,iblk)
          e => blk(t)%end(:, iblk)
          xloc = xloc + scl * sum(a(s(1):e(1),s(2):e(2),s(3):e(3)))
       enddo
       !$OMP ENDDO
    endif
  end subroutine loops_ab
