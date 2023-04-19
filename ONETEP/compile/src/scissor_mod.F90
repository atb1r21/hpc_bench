! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                             Scissor Module                                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   The subroutines in this file were written by:                             !
!                                                                             !
!   Nelson Yeung                                                              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module scissor

  use constants, only: DP

  private

  public :: scissor_shift_denskern

contains

  subroutine scissor_shift_denskern(shifted_denskern,denskern,ngwf_basis, &
       mdl,inv_overlap)

    !========================================================================!
    ! This subroutine gets the sum of the shifted valence and conduction
    ! denskern: shifted_denskern = vK + c(S^-1-K)
    !------------------------------------------------------------------------!
    ! shifted_denskern   (inout) : Shifted Density kernel
    ! denskern              (in) : Density kernel
    ! ngwf_basis            (in) : NGWF basis
    ! mdl                   (in) : Model
    !------------------------------------------------------------------------!
    ! Written by Nelson Yeung in January 2019.                               !
    !========================================================================!

    use comms, only: pub_my_proc_id
    use function_basis, only: FUNC_BASIS
    use model_type, only: REGION
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_scissor_ngroups, pub_scissor_groups, &
         pub_scissor_group_nsp, pub_scissor_group_shifts
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_axpy, &
         sparse_index_length, sparse_generate_index, sparse_get_block, &
         sparse_put_block, sparse_destroy
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Argument
    type(SPAM3), intent(inout)     :: shifted_denskern       ! output density kernel
    type(SPAM3), intent(in)        :: denskern               ! input density kernel
    type(FUNC_BASIS), intent(in)   :: ngwf_basis             ! NGWF basis type
    type(REGION), intent(in)       :: mdl                    ! Model region
    type(SPAM3), intent(in)        :: inv_overlap

    ! Local variables
    integer                        :: iat, jat               ! atom counters
    integer                        :: idxlen                 ! index length
    integer                        :: ierr
    integer                        :: ii
    integer                        :: jdx
    integer                        :: loc_iat                ! atom counter on local proc
    integer                        :: max_ngwfs_on_atom      ! max NGWFs on any atom
    integer                        :: orig_iat, orig_jat     ! atom counters input file order
    integer, allocatable           :: idx(:)                 ! sparse index
    complex(kind=DP), allocatable  :: zblock(:,:)            ! complex sparse matrix block
    real(kind=DP), allocatable     :: dblock(:,:)            ! real sparse matrix block
    real(kind=DP)                  :: val_shift, cond_shift  ! scissor shifts
    logical                        :: iscmplx
    logical                        :: region_only            ! shift only the specified regions
    type(SPAM3)                    :: cond_denskern

    iscmplx = denskern%iscmplx
    max_ngwfs_on_atom = maxval(ngwf_basis%num_on_atom)

    call sparse_create(cond_denskern,denskern)

    ! create conduction density kernel: (S^-1-K)
    call sparse_copy(cond_denskern,inv_overlap)
    call sparse_axpy(cond_denskern,denskern,-1.0_DP)

    ! allocate complex or real block
    if (iscmplx) then
       allocate(zblock(max_ngwfs_on_atom,max_ngwfs_on_atom),stat=ierr)
       call utils_alloc_check('scissor_shift_denskern', &
            'zblock',ierr)
       zblock = (0.0_DP,0.0_DP)
    else
       allocate(dblock(max_ngwfs_on_atom,max_ngwfs_on_atom),stat=ierr)
       call utils_alloc_check('scissor_shift_denskern', &
            'dblock',ierr)
       dblock = (0.0_DP,0.0_DP)
    end if

    ! generate index of sparse matrix
    idxlen = sparse_index_length(denskern)
    allocate(idx(idxlen),stat=ierr)
    call utils_alloc_check('scissor_shift_denskern','idx',ierr)
    call sparse_generate_index(idx, denskern)

    ! Check if applying scissor shift to specified regions only
    region_only = .false.
    do ii=1,pub_scissor_ngroups
       if (pub_scissor_group_nsp(ii) > 1) then
          region_only = .true.
       endif
    end do

    ! loop over first atom
    do loc_iat=1,mdl%par%num_atoms_on_proc(pub_my_proc_id)
       iat = loc_iat + mdl%par%first_atom_on_proc(pub_my_proc_id) - 1
       orig_iat = par%orig_atom(iat)

       ! loop over second atom
       do jdx=idx(loc_iat),idx(loc_iat+1)-1
          jat = idx(jdx)
          orig_jat = par%orig_atom(jat)
          val_shift = 0.0_DP
          cond_shift = 0.0_DP

          ! Get valence and conduction scissor shift
          if (region_only) then
             do ii=1,pub_scissor_ngroups
                if (any(pub_scissor_groups(1:pub_scissor_group_nsp(ii),ii) == &
                     mdl%elements(orig_iat)%species_id) .and. &
                     any(pub_scissor_groups(1:pub_scissor_group_nsp(ii),ii) == &
                     mdl%elements(orig_jat)%species_id)) then
                   val_shift = pub_scissor_group_shifts(2,ii)
                   cond_shift = pub_scissor_group_shifts(1,ii)
                   exit
                endif
             end do
          else
             ! get average shift
             do ii=1,pub_scissor_ngroups
                if (mdl%elements(orig_iat)%species_id == &
                     pub_scissor_groups(1,ii) .or. &
                     mdl%elements(orig_jat)%species_id == &
                     pub_scissor_groups(1,ii)) then
                   val_shift = val_shift + pub_scissor_group_shifts(2,ii)
                   cond_shift = cond_shift + pub_scissor_group_shifts(1,ii)

                   ! if same atom just double the shifts
                   if (mdl%elements(orig_iat)%species_id == &
                        mdl%elements(orig_jat)%species_id) then
                      val_shift = val_shift + pub_scissor_group_shifts(2,ii)
                      cond_shift = cond_shift + pub_scissor_group_shifts(1,ii)
                      exit
                   end if
                end if
             end do
             val_shift = val_shift * 0.5_DP
             cond_shift = cond_shift * 0.5_DP
          end if

          if(iscmplx) then
             call sparse_get_block(zblock,denskern,jat,iat)
             zblock = zblock * val_shift
             call sparse_put_block(zblock,shifted_denskern,jat,iat)

             call sparse_get_block(zblock,cond_denskern,jat,iat)
             zblock = zblock * cond_shift
             call sparse_put_block(zblock,cond_denskern,jat,iat)
          else
             call sparse_get_block(dblock,denskern,jat,iat)
             dblock = dblock * val_shift
             call sparse_put_block(dblock,shifted_denskern,jat,iat)

             call sparse_get_block(dblock,cond_denskern,jat,iat)
             dblock = dblock * cond_shift
             call sparse_put_block(dblock,cond_denskern,jat,iat)
          end if
       end do
    end do

    call sparse_axpy(shifted_denskern,cond_denskern,1.0_DP)

    if (iscmplx) then
       deallocate(zblock,stat=ierr)
       call utils_dealloc_check('scissor_shift_denskern', &
            'zblock',ierr)
    else
       deallocate(dblock,stat=ierr)
       call utils_dealloc_check('scissor_shift_denskern', &
            'dblock',ierr)
    end if

    deallocate(idx,stat=ierr)
    call utils_dealloc_check('scissor_shift_denskern', &
         'idx',ierr)

    call sparse_destroy(cond_denskern)

  end subroutine scissor_shift_denskern

end module scissor
