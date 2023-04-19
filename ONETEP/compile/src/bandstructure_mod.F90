! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=================================================================!
!                                                                 !
!                    Bandstructure Module                         !
!                                                                 !
! This module calculates the bandstructure using the NGWFs as a   !
! basis.                                                          !
!-----------------------------------------------------------------!
! Written by Peter Haynes on 10/10/2006 - "k.p" style             !
! Added "tight-binding" style on 16/1/2010                        !
! This file created by Nicholas Hine 09/03/2010 by moving routine !
! from properties_mod.                                            !
! Modified to use NGWF_REP and NGWF_HAM by Nicholas Hine in       !
! October 2010.                                                   !
!=================================================================!

module bandstructure

  implicit none

  private

  public :: bandstructure_calculate
  public :: bsunfold_calculate

contains

  subroutine bandstructure_calculate(ham, rep, ngwf_basis, proj_basis, &
       nl_projectors, rhoij, mdl, ham_type)

    !======================================================================!
    ! This subroutine calculates the bandstructure using the NGWFs as a    !
    ! basis.                                                               !
    !----------------------------------------------------------------------!
    ! Written by Peter Haynes on 10/10/2006 - "k.p" style                  !
    ! Added "tight-binding" style on 16/1/2010                             !
    ! Perturbative spin-orbit couplings added by JM Escartin, Summer 2016. !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,     !
    ! August 2018                                                          !
    !======================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(.DOT.), operator(-), operator(*), &
         geometry_magnitude, operator(+)
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_bs_unfold, pub_bs_kpoint_path_end, &
         pub_bs_num_eigenvalues, pub_bs_kpoint_path_spacing, pub_bs_method, &
         pub_bs_kpoint_path_length, pub_bs_kpoint_path_start, pub_rootname, &
         pub_cond_calculate, pub_debug_on_root, pub_num_spins, &
         pub_perturbative_soc, PUB_1K
    use sparse_embed, only: SPAM3_EMBED
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_abort, utils_banner

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(in) :: ham              ! Hamiltonian matrices
    type(NGWF_REP), intent(in) :: rep              ! NGWF Representation
    type(FUNC_BASIS), intent(in) :: ngwf_basis     ! NGWF basis type
    type(FUNC_BASIS), intent(in) :: proj_basis     ! Projector basis type
    type(PROJECTOR_SET), intent(inout) :: nl_projectors
    type(SPAM3_EMBED), intent(in) :: rhoij(pub_num_spins) ! Projector denskern
    type(MODEL), intent(in) :: mdl
    ! lr408: optional string denoting which basis for conduction calculations
    character(len=*), intent(in) :: ham_type

    ! Local variables
    type(POINT) :: kstart                 ! k-point at start of path segment
    type(POINT) :: kend                   ! k-point at end of path segment
    type(POINT) :: kseg                   ! path segment
    type(POINT) :: kptdiff                ! difference between k-points
    type(POINT), allocatable :: kpt(:)    ! k-point
    integer :: ierr                       ! Error flag
    integer :: is                         ! Spin counters
    integer :: ipath                      ! Bandstructure path counter
    integer :: nsegpts                    ! Number of points along path segment
    integer :: isegpt                     ! Path segment point counter
    integer :: nkpts                      ! Total number of kpts
    integer :: ikpt,jkpt                  ! k-point counter
    integer :: ingwf                      ! NGWF counters
    integer :: nlower,nupper              ! Limits on NGWF bounds for output
    integer :: nkpts_red                  ! Reduced list of unique k-points
    integer :: bands_unit                 ! I/O unit for .bands file
    integer :: agr_unit                   ! I/O unit for .agr file
    integer :: nelec                      ! Number of electrons
    integer :: spin_dim                   ! Number of independent spin dims
    integer :: ngwf_num                   ! Number of NGWFs
    integer, allocatable :: uni_entry(:)  ! Pointer to unique entry
    logical :: unique                     ! Flag for reducing k-points
    logical :: basic_unfold               ! Flag for simple diagonal unfolding
                                          !     for the Brillouin zone
    real(kind=DP) :: seglen               ! Path segment length
    real(kind=DP) :: pathlen              ! Total path length
    real(kind=DP) :: kfrac(3)             ! k-point in fractional coordinates
    real(kind=DP) :: efermi               ! Fermi energies
    real(kind=DP), allocatable :: eigvl_kpts(:,:,:)  ! All eigenvalues (red)
    complex(kind=DP), allocatable :: eigst_kpts(:,:,:,:) ! All eigenstates (red)
    real(kind=DP), allocatable :: kunfold(:,:,:,:)  ! Unfolded k-point
    character(len=90) :: filename                ! Output filename

    ! jme: KPOINTS_DANGER
    ! First k-point enforced in all uses of n_occ via PUB_1K.

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering bandstructure_calculate'

    if (pub_on_root) write(stdout,'(/a)') utils_banner('=', &
         'Bandstructure Calculation')

    ! Check details of unfolding
    basic_unfold = (pub_bs_unfold(1) > 0 .and. pub_bs_unfold(2) > 0 .and. &
         pub_bs_unfold(3) > 0)
    if (basic_unfold) pub_bs_unfold = max(pub_bs_unfold,1)

    ! Count the number of k-points needed
    nkpts = 0
    pathlen = 0.0_DP
    do ipath=1,pub_bs_kpoint_path_length
       kstart = pub_bs_kpoint_path_start(1,ipath) * mdl%cell%b1 + &
            pub_bs_kpoint_path_start(2,ipath) * mdl%cell%b2 + &
            pub_bs_kpoint_path_start(3,ipath) * mdl%cell%b3
       kend = pub_bs_kpoint_path_end(1,ipath) * mdl%cell%b1 + &
            pub_bs_kpoint_path_end(2,ipath) * mdl%cell%b2 + &
            pub_bs_kpoint_path_end(3,ipath) * mdl%cell%b3
       kseg = kend - kstart
       seglen = geometry_magnitude(kseg)
       pathlen = pathlen + seglen
       nsegpts = max(int(seglen / pub_bs_kpoint_path_spacing),1)
       nkpts = nkpts + nsegpts + 1
    end do

    ! Compile a complete list of k-points
    allocate(kpt(nkpts),stat=ierr)
    call utils_alloc_check('bandstructure_calculate','kpt',ierr)
    ikpt = 0
    do ipath=1,pub_bs_kpoint_path_length
       kstart = pub_bs_kpoint_path_start(1,ipath) * mdl%cell%b1 + &
            pub_bs_kpoint_path_start(2,ipath) * mdl%cell%b2 + &
            pub_bs_kpoint_path_start(3,ipath) * mdl%cell%b3
       kend = pub_bs_kpoint_path_end(1,ipath) * mdl%cell%b1 + &
            pub_bs_kpoint_path_end(2,ipath) * mdl%cell%b2 + &
            pub_bs_kpoint_path_end(3,ipath) * mdl%cell%b3
       kseg = kend - kstart
       seglen = geometry_magnitude(kseg)
       nsegpts = max(int(seglen / pub_bs_kpoint_path_spacing),1)
       do isegpt=0,nsegpts
          ikpt = ikpt + 1
          kpt(ikpt) = kstart + (isegpt / real(nsegpts,kind=DP)) * kseg
       end do
    end do

    ! Now compile a reduced list of unique k-points avoiding repetition
    ! (e.g. at start and end of consecutive segments)
    allocate(uni_entry(nkpts),stat=ierr)
    call utils_alloc_check('bandstructure_calculate','uni_entry',ierr)
    nkpts_red = 0
    do ikpt=1,nkpts
       unique = .true.
       do jkpt=1,ikpt-1
          kptdiff = kpt(ikpt) - kpt(jkpt)
          if (geometry_magnitude(kptdiff) < 1.0e-6_DP) then
             unique = .false.
             uni_entry(ikpt) = uni_entry(jkpt)
             exit
          end if
       end do
       if (unique) then
          nkpts_red = nkpts_red + 1
          uni_entry(ikpt) = nkpts_red
       end if
    end do

    ! Set spin dimensions.
    if (pub_perturbative_soc) then
       spin_dim = 1
       ngwf_num = 2 * ngwf_basis%num
    else
       spin_dim = pub_num_spins
       ngwf_num = ngwf_basis%num
    end if

    ! Allocate workspace for diagonalisation and unfolding
    allocate(eigst_kpts(ngwf_num,ngwf_num,nkpts_red,spin_dim), stat=ierr)
    call utils_alloc_check('bandstructure_calculate','eigst_kpts',ierr)
    allocate(eigvl_kpts(ngwf_num,nkpts_red,spin_dim),stat=ierr)
    call utils_alloc_check('bandstructure_calculate','eigvl_kpts',ierr)
    if (basic_unfold) then
       allocate(kunfold(ngwf_num,3,nkpts,spin_dim),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','kunfold',ierr)
    end if

    eigst_kpts = cmplx(0.0_DP, 0.0_DP, kind=DP)
    eigvl_kpts = 0.0_DP

    call bandstructure_core(ham, rep, ngwf_basis, proj_basis, &
         nl_projectors, rhoij, mdl, pub_bs_method, spin_dim, &
         ngwf_num, nkpts_red, nkpts, kpt, uni_entry, basic_unfold, eigvl_kpts, &
         eigst_kpts = eigst_kpts, kunfold = kunfold)

    ! Write results to output file
    if (pub_perturbative_soc) then
       nelec = rep%n_occ(1,PUB_1K)
       select case (pub_num_spins)
          case (1)
             nelec = 2 * nelec
          case (2)
             nelec = nelec + rep%n_occ(2,PUB_1K)
       end select
    end if
    ! ndmh: write all bands if bs_num_eigenvalues is left at default value
    if (pub_bs_num_eigenvalues < 0) then
       if (pub_perturbative_soc) then
          pub_bs_num_eigenvalues = nelec
       else
          pub_bs_num_eigenvalues = maxval(rep%n_occ(:,PUB_1K))
       end if
    end if
    if (pub_on_root .and. pub_bs_num_eigenvalues > 0) then
       write(stdout,'(/a)') &
            ' + ============================================================= +'
       write(stdout,'(a)') &
            ' +                     Electronic energies                       +'
       write(stdout,'(a)') &
            ' +                     -------------------                       +'
       write(stdout,'(a)') &
            ' +                                                               +'
       if (basic_unfold) then
          write(stdout,'(a)') ' +   Band number   Energy in Ha          &
               &Unfolded k-point        +'
       else
          write(stdout,'(a)') ' +   Band number   Energy in Ha          &
               &                        +'
       end if
       write(stdout,'(a)') &
            ' + ============================================================= +'
       write(stdout,'(a)') &
            ' +                                                               +'
       do is = 1, spin_dim
          if (.not. pub_perturbative_soc) nelec = rep%n_occ(is,PUB_1K)
          nlower = max(nelec-pub_bs_num_eigenvalues+1,1)
          nupper = min(nelec+pub_bs_num_eigenvalues,ngwf_num)
          ikpt = 0
          do ipath=1,pub_bs_kpoint_path_length
             kstart = pub_bs_kpoint_path_start(1,ipath) * mdl%cell%b1 + &
                  pub_bs_kpoint_path_start(2,ipath) * mdl%cell%b2 + &
                  pub_bs_kpoint_path_start(3,ipath) * mdl%cell%b3
             kend = pub_bs_kpoint_path_end(1,ipath) * mdl%cell%b1 + &
                  pub_bs_kpoint_path_end(2,ipath) * mdl%cell%b2 + &
                  pub_bs_kpoint_path_end(3,ipath) * mdl%cell%b3
             kseg = kend - kstart
             seglen = geometry_magnitude(kseg)
             nsegpts = max(int(seglen / pub_bs_kpoint_path_spacing),1)
             do isegpt=0,nsegpts
                ikpt = ikpt + 1
                kfrac = pub_bs_kpoint_path_start(:,ipath) + &
                     (isegpt / real(nsegpts,kind=DP)) * &
                     (pub_bs_kpoint_path_end(:,ipath) - &
                     pub_bs_kpoint_path_start(:,ipath))
                write(stdout,'(a)') ' + --------------------------------------&
                     &----------------------- +'
                write(stdout,'(a,i1,a,i5,a,f8.5,1x,f8.5,1x,f8.5,a)') &
                     ' +  Spin=',is,' kpt=',ikpt,' (',kfrac,')                +'
                write(stdout,'(a)') ' + --------------------------------------&
                     &----------------------- +'
                write(stdout,'(a)') ' +                                       &
                     &                        +'
                do ingwf=nlower,nupper
                   if (basic_unfold) then
                      kfrac = kunfold(ingwf,:,ikpt,is)
                      write(stdout,'(a,i12,1x,f15.6,5x,a,3f8.4,a,4x,a)') &
                           ' +',ingwf,eigvl_kpts(ingwf,uni_entry(ikpt),is),'(',&
                           kfrac,')','+'
                   else
                      write(stdout,'(a,i12,1x,f15.6,35x,a)') ' +',ingwf, &
                           eigvl_kpts(ingwf,uni_entry(ikpt),is),'+'
                   end if
                   if (ingwf == nelec) write(stdout,'(a)') &
                     ' +         ... Fermi level ...                          &
                     &         +'
                end do
                write(stdout,'(a)') ' +                                       &
                     &                        +'
             end do
          end do
       end do
       write(stdout,'(a)') &
            ' + ============================================================= +'
    end if

    ! Write out .bands file
    if (pub_on_root) then

       ! Open output file
       bands_unit = utils_unit()
       if (.not.pub_cond_calculate) then
          write(filename,'(2a)') trim(pub_rootname),'_BS.bands'
       else
          write(filename,'(4a)') trim(pub_rootname),'_',trim(ham_type),'_BS.bands'
       end if
       open(bands_unit,file=trim(filename),iostat=ierr)
       if (ierr /= 0) then
          call utils_abort('Error in bandstructure_calculate: opening "'//&
               trim(filename)//'" failed with code ',ierr)
       end if

       write(bands_unit,'(a,i6)') 'Number of k-points',nkpts
       write(bands_unit,'(a,i2)') 'Number of spin components',pub_num_spins
       write(bands_unit,'(a,i2)') 'Number of independent spin components', spin_dim
       write(bands_unit,'(a,2i8)') 'Number of electrons ',rep%n_occ(:,PUB_1K)
       if (pub_perturbative_soc) then
          if (nelec > 0) then
             efermi = maxval(eigvl_kpts(nelec,:,1))
          else
             efermi = 0_DP
          endif
          write(bands_unit,'(a,i8)') 'Number of eigenvalues ',ngwf_num
          write(bands_unit,'(a,f12.6)') 'Fermi energy (in atomic units) ', &
               efermi
       else
          if (rep%n_occ(1,PUB_1K)>0) then
             efermi = maxval(eigvl_kpts(rep%n_occ(1,PUB_1K),:,1))
          else
             efermi = 0_DP
          endif
          if (pub_num_spins == 1) then
             write(bands_unit,'(a,i8)') 'Number of eigenvalues ',ngwf_num
             write(bands_unit,'(a,f12.6)') 'Fermi energy (in atomic units) ', &
                  efermi
          else
             write(bands_unit,'(a,2i8)') 'Number of eigenvalues ',ngwf_num, &
                  ngwf_num
             if (rep%n_occ(2,PUB_1K)>0) then
                efermi = max(efermi,maxval(eigvl_kpts(rep%n_occ(2,PUB_1K),:,2)))
             else
                efermi = 0_DP
             endif
             write(bands_unit,'(a,2f12.6)') 'Fermi energies (in atomic units) ', &
                  efermi,efermi
          end if
       end if
       write(bands_unit,'(a)') 'Unit cell vectors'
       write(bands_unit,'(3f12.6)') mdl%cell%a1%x,mdl%cell%a1%y,mdl%cell%a1%z
       write(bands_unit,'(3f12.6)') mdl%cell%a2%x,mdl%cell%a2%y,mdl%cell%a2%z
       write(bands_unit,'(3f12.6)') mdl%cell%a3%x,mdl%cell%a3%y,mdl%cell%a3%z
       ikpt = 0
       do ipath=1,pub_bs_kpoint_path_length
          kstart = pub_bs_kpoint_path_start(1,ipath) * mdl%cell%b1 + &
               pub_bs_kpoint_path_start(2,ipath) * mdl%cell%b2 + &
               pub_bs_kpoint_path_start(3,ipath) * mdl%cell%b3
          kend = pub_bs_kpoint_path_end(1,ipath) * mdl%cell%b1 + &
               pub_bs_kpoint_path_end(2,ipath) * mdl%cell%b2 + &
               pub_bs_kpoint_path_end(3,ipath) * mdl%cell%b3
          kseg = kend - kstart
          seglen = geometry_magnitude(kseg)
          nsegpts = max(int(seglen / pub_bs_kpoint_path_spacing),1)
          do isegpt=0,nsegpts
             ikpt = ikpt + 1
             kfrac = pub_bs_kpoint_path_start(:,ipath) + &
                  (isegpt / real(nsegpts,kind=DP)) * &
                  (pub_bs_kpoint_path_end(:,ipath) - &
                  pub_bs_kpoint_path_start(:,ipath))
             write(bands_unit,'(a,i6,4f12.8)') 'K-point',ikpt,kfrac,1.0_DP/nkpts
             do is = 1, spin_dim
                write(bands_unit,'(a,i2)') 'Spin component',is
                do ingwf=1,ngwf_num
                   write(bands_unit,'(f14.8)') &
                        eigvl_kpts(ingwf,uni_entry(ikpt),is)
                end do
             end do
          end do
       end do

       ! Close output file
       close(bands_unit,iostat=ierr)
       if (ierr /= 0) then
          call utils_abort('Error in bandstructure_calculate: closing "'//&
               trim(filename)//'" failed with code ',ierr)
       end if

    end if

    ! Write out .agr file
    if (pub_on_root) then

       ! Open output file
       agr_unit = utils_unit()
       write(filename,'(2a)') trim(pub_rootname),'_BS.agr'
       open(agr_unit,file=trim(filename),iostat=ierr)
       if (ierr /= 0) then
          call utils_abort('Error in bandstructure_calculate: opening "'//&
               trim(filename)//'" failed with code ',ierr)
       end if

       write(agr_unit,'(a)') '@version 50112'
       write(agr_unit,'(a)') '@with g0'
       write(agr_unit,'(a,f12.6)') '@    world xmin ',0.0_DP
       write(agr_unit,'(a,f12.6)') '@    world xmax ',pathlen
       write(agr_unit,'(a,f12.6)') '@    world ymin ', &
            real(floor(minval(eigvl_kpts)),kind=DP)
       write(agr_unit,'(a,f12.6)') '@    world ymax ', &
            real(ceiling(maxval(eigvl_kpts)),kind=DP)
       write(agr_unit,'(a)') '@    yaxis  tick on'
       write(agr_unit,'(a,f12.6)') '@    yaxis  tick major ',1.0_DP
       write(agr_unit,'(a)') '@    yaxis  tick minor off'
       write(agr_unit,'(a)') '@    xaxis  tick major off'
       write(agr_unit,'(a)') '@    xaxis  tick minor off'
       write(agr_unit,'(a)') '@    xaxis  ticklabel on'
       write(agr_unit,'(a)') '@    xaxis  ticklabel type spec'
       write(agr_unit,'(a)') '@    xaxis  tick type spec'
       write(agr_unit,'(a,i2,a,f12.6)') '@    xaxis  tick ',0,', ',0.0_DP
       write(agr_unit,'(a,2(i2,a))') '@    xaxis  ticklabel ',0,', "\6k\4\s', &
            0,'"'
       seglen = 0.0_DP
       do ipath=1,pub_bs_kpoint_path_length
          kstart = pub_bs_kpoint_path_start(1,ipath) * mdl%cell%b1 + &
               pub_bs_kpoint_path_start(2,ipath) * mdl%cell%b2 + &
               pub_bs_kpoint_path_start(3,ipath) * mdl%cell%b3
          kend = pub_bs_kpoint_path_end(1,ipath) * mdl%cell%b1 + &
               pub_bs_kpoint_path_end(2,ipath) * mdl%cell%b2 + &
               pub_bs_kpoint_path_end(3,ipath) * mdl%cell%b3
          kseg = kend - kstart
          seglen = seglen + geometry_magnitude(kseg)
          write(agr_unit,'(a,i2,a,f12.6)') '@    xaxis  tick ',ipath,', ',seglen
          write(agr_unit,'(a,2(i2,a))') '@    xaxis  ticklabel ',ipath, &
               ', "\6k\4\s',ipath,'"'
       end do
       write(agr_unit,'(a,i3)') '@    xaxis  tick spec ', &
            pub_bs_kpoint_path_length+1

       do is = 1, spin_dim
          do ingwf = 1, ngwf_num
             write(agr_unit,'(a,i3)') '@TYPE xy'
             ikpt = 0
             pathlen = 0.0_DP
             do ipath=1,pub_bs_kpoint_path_length
                kstart = pub_bs_kpoint_path_start(1,ipath) * mdl%cell%b1 + &
                     pub_bs_kpoint_path_start(2,ipath) * mdl%cell%b2 + &
                     pub_bs_kpoint_path_start(3,ipath) * mdl%cell%b3
                kend = pub_bs_kpoint_path_end(1,ipath) * mdl%cell%b1 + &
                     pub_bs_kpoint_path_end(2,ipath) * mdl%cell%b2 + &
                     pub_bs_kpoint_path_end(3,ipath) * mdl%cell%b3
                kseg = kend - kstart
                seglen = geometry_magnitude(kseg)
                nsegpts = max(int(seglen / pub_bs_kpoint_path_spacing),1)
                do isegpt=0,nsegpts
                   ikpt = ikpt + 1
                   if (isegpt > 0) pathlen = pathlen + seglen / nsegpts
                   write(agr_unit,'(2f12.6)') pathlen, &
                        eigvl_kpts(ingwf,uni_entry(ikpt),is)
                end do
             end do
             write(agr_unit,'(a,i3)') '&'
          end do
       end do

       ! Close output file
       close(agr_unit,iostat=ierr)
       if (ierr /= 0) then
          call utils_abort('Error in bandstructure_calculate: closing "'//&
               trim(filename)//'" failed with code ',ierr)
       end if

    end if

    ! Destroy workspace
    if (basic_unfold) then
       deallocate(kunfold,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','kunfold',ierr)
    end if
    deallocate(eigst_kpts,stat=ierr)
    call utils_dealloc_check('bandstructure_calculate','eigst_kpts',ierr)
    deallocate(eigvl_kpts,stat=ierr)
    call utils_dealloc_check('bandstructure_calculate','eigvl_kpts',ierr)
    deallocate(uni_entry,stat=ierr)
    call utils_dealloc_check('bandstructure_calculate','uni_entry',ierr)
    deallocate(kpt,stat=ierr)
    call utils_dealloc_check('bandstructure_calculate','kpt',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving bandstructure_calculate'

  end subroutine bandstructure_calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine bandstructure_core(ham, rep, ngwf_basis, proj_basis, &
       nl_projectors, rhoij, mdl, bs_method, spin_dim, ngwf_num,&
       nkpts_red, nkpts_total, kpt, uni_entry, basic_unfold, eigvl_kpts, &
       eigst_kpts, dense_eigst, kunfold)

    !======================================================================!
    ! This subroutine calculates the bandstructure using the NGWFs as a    !
    ! basis.                                                               !
    !----------------------------------------------------------------------!
    ! Written by Peter Haynes on 10/10/2006 - "k.p" style                  !
    ! Added "tight-binding" style on 16/1/2010                             !
    ! Perturbative spin-orbit couplings added by JM Escartin, Summer 2016. !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,     !
    ! August 2018                                                          !
    !======================================================================!

    use augmentation, only: aug_nonlocal_mat
    use comms, only: pub_total_num_procs, pub_my_proc_id, &
         comms_bcast, comms_barrier, comms_reduce, pub_on_root
    use constants, only: DP, TWO_PI, cmplx_1, stdout, real_size, cmplx_size, &
         LONG
    use dense, only: DEM, dense_put_col
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(.DOT.), operator(-), operator(*), &
         geometry_magnitude, operator(+)
    use hamiltonian, only: hamiltonian_soc_matrices
    use integrals, only: integrals_grad
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use paw, only: paw_dij_so
    use projectors, only: PROJECTOR_SET, projectors_func_ovlp_box
    use pseudopotentials, only: pseudopotentials_nonlocal_mat, pseudo_get_dij
    use rundat, only: pub_any_nl_proj, pub_aug, pub_bs_unfold, &
         pub_debug_on_root, pub_num_spins, pub_perturbative_soc
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_copy, &
         sparse_axpy, sparse_convert, sparse_get_block, sparse_put_block, &
         sparse_index_length, sparse_generate_index, sparse_product, &
         sparse_transpose, sparse_transpose_structure
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_copy, &
         sparse_embed_destroy, sparse_embed_extract_from_array, &
         sparse_embed_destroy_extracted_array
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, &
         utils_assert

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(in) :: ham              ! Hamiltonian matrices
    type(NGWF_REP), intent(in) :: rep              ! NGWF Representation
    type(FUNC_BASIS), intent(in) :: ngwf_basis     ! NGWF basis type
    type(FUNC_BASIS), intent(in) :: proj_basis     ! Projector basis type
    type(PROJECTOR_SET), intent(inout) :: nl_projectors
    type(SPAM3_EMBED), intent(in) :: rhoij(pub_num_spins) ! Projector denskern
    type(MODEL), intent(in) :: mdl
    character(len=2), intent(in) :: bs_method
    integer, intent(in) :: spin_dim, ngwf_num, nkpts_red, nkpts_total
    type(POINT), intent(in) :: kpt(nkpts_total)    ! k-point
    integer, intent(in) :: uni_entry(nkpts_total)
    logical, intent(in) :: basic_unfold ! a basic diagonal version for unfolding
    ! the following only make sense if basic_unfold is true
    real(kind=DP), intent(inout) :: eigvl_kpts(ngwf_num,nkpts_red,spin_dim)
    complex(kind=DP), intent(inout), optional :: eigst_kpts(ngwf_num,ngwf_num, &
         nkpts_red,spin_dim)
    type(DEM), intent(inout), optional :: dense_eigst(nkpts_red,spin_dim)
    real(kind=DP), intent(inout), optional :: kunfold(ngwf_num,3, &
         nkpts_total,spin_dim)  ! Unfolded k-point



    ! LAPACK subroutines
    external :: ztrsm, zheev, zhegst, dpotrf, zhegv

    ! Local variables
    type(SPAM3), allocatable :: ham0(:)   ! k-independent part of Hamiltonian
    type(SPAM3_EMBED), allocatable :: nl_k(:)   ! Nonlocal matrix at k-point
    type(SPAM3), allocatable :: ham_ks(:) ! Hamiltonian for k-point and spin
    type(SPAM3), allocatable :: over_ks(:) ! Overlap for k-point and spin
    type(SPAM3_EMBED), allocatable :: hamso_k(:,:) ! Hamiltonian (with SOC) for k-point
    type(SPAM3) :: grad(3)                ! Grad matrix
    type(POINT) :: atdisp                 ! atom displacement
    integer :: ierr                       ! Error flag
    integer :: is, isp                    ! Spin counters
    integer :: dim                        ! Dimension counter
    integer :: ikpt,jkpt                  ! k-point counter
    integer :: ingwf,jngwf                ! NGWF counters
    integer :: max_ngwfs_on_atom          ! Max number of NGWFs on any atom
    integer :: proc                       ! Proc counter
    integer :: lwork                      ! Workspace length
    integer :: iat,jat                    ! Atom counters
    integer :: orig_iat, orig_jat         ! Atom counters in input file order
    integer :: loc_iat                    ! Atom counter on local proc
    integer :: jdx                        ! Sparse matrix counters
    integer :: ispec                      ! Species of atom iat
    integer :: wrapped                    ! Flag for wrapping
    integer :: hamidxlen,overidxlen       ! Sparse index lengths
    integer :: max_degeneracy             ! Maximum degeneracy
    integer, allocatable :: translate(:,:)! Translations for unfolding
    integer, allocatable :: hamidx(:)     ! Sparse Hamiltonian index
    integer, allocatable :: overidx(:)    ! Sparse overlap index
    logical :: found                      ! Flag for finding equivalent atom
    logical, allocatable :: flag(:)       ! Flag when finding equivalent atoms
    real(kind=DP), parameter :: tol = 1.0e-8_DP
    real(kind=DP), parameter :: recip_twopi = 1.0_DP / TWO_PI
    real(kind=DP) :: distsq               ! Squared distance
    real(kind=DP) :: phase                ! Translation phase
    real(kind=DP) :: kcart(3)             ! k-point in Cartesians
    real(kind=DP) :: kfrac(3)             ! k-point in fractional coordinates
    real(kind=DP) :: trvec(3)             ! Translation vector in real-space
    real(kind=DP) :: trpos(3)             ! Translated atom position
    real(kind=DP) :: dispfrac(3)          ! Displacement vector (fractional)
    real(kind=DP) :: dispabs(3)           ! Displacement vector (absolute)
    real(kind=DP), allocatable :: ov_sq(:,:)  ! Overlap in dense format
    real(kind=DP), allocatable :: dwork(:,:)  ! Diagonalization workspace
    real(kind=DP), allocatable :: en_ks(:)    ! Eigenvalues for k-point/spin
    real(kind=DP), allocatable :: kpt_ks(:,:) ! Unfolded k-point for k-point/sp
    real(kind=DP), allocatable :: posfrac(:,:)   ! Fractional atomic positions
    complex(kind=DP) :: i2pi                     ! 2i*pi
    complex(kind=DP) :: zphase                   ! Complex translation phase
    complex(kind=DP), allocatable :: hm_sq(:,:)  ! Hamiltonian in dense format
    complex(kind=DP), allocatable :: zov_sq(:,:) ! Overlap in dense format
    complex(kind=DP), allocatable :: bmat(:,:)   ! Diagonalization workspace
    complex(kind=DP), allocatable :: zwork(:)    ! Diagonalization workspace
    complex(kind=DP), allocatable :: trwfn(:)    ! Translated wavefunction
    complex(kind=DP), allocatable :: block(:,:)  ! Sparse matrix block
    complex(kind=DP), allocatable :: psitrpsi(:,:)  ! Inner products
    complex(kind=DP), allocatable :: zlambda(:)  ! Complex eigenvalues
    complex(kind=DP), allocatable :: aux_array(:,:) ! Workspace for SOC merge
    character(len=1024) :: error_message
    integer :: icol
    integer(kind=LONG) :: memtotal_real, memtotal_cmplx


    ! jme: KPOINTS_DANGER
    ! First k-point enforced in all uses of n_occ via PUB_1K.

    ! ndmh: SPAM3 approach to nonlocal pseudopotentials
    type(SPAM3_EMBED) :: sp_overlap

    ! jcap: needed as input to pseudopotentials_nonlocal_mat
    type(SPAM3_EMBED) :: dij
    ! jcap: temporary arrays for workaround for SPAM3 spin arrays for
    ! use with aug_nonlocal_mat
    type(SPAM3) :: aug_array(pub_num_spins), kern_array(pub_num_spins), &
         ham_array(1)
    type(SPAM3_EMBED) :: ham_embed_array(1)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering bandstructure_core'

    !pdh: define constant to avoid problem with gfortran
    i2pi = cmplx(0.0_DP,-TWO_PI,kind=DP)

    ! Allocate workspace for diagonalisation and unfolding
    ! Memory usage check
    memtotal_real = 0
    memtotal_cmplx = 0
    memtotal_real = memtotal_real + ngwf_num
    memtotal_cmplx = memtotal_cmplx + ngwf_num*ngwf_num
    memtotal_real = memtotal_real + ngwf_num
    memtotal_cmplx = memtotal_cmplx + ngwf_num*2
    memtotal_cmplx = memtotal_cmplx + ngwf_num*ngwf_num
    memtotal_cmplx = memtotal_cmplx + ngwf_num*ngwf_num
    if (pub_perturbative_soc) then
       memtotal_cmplx = memtotal_cmplx + ngwf_basis%num*ngwf_basis%num
    end if
    if (pub_on_root) then
       write(stdout,'(a)') 'Estimated total memory usage:'
       write(stdout,'(a,i20)') 'Reals: ',memtotal_real*real_size
       write(stdout,'(a,i20)') 'Cmplx: ',memtotal_cmplx*cmplx_size
    end if

    allocate(en_ks(ngwf_num),stat=ierr)
    call utils_alloc_check('bandstructure_core','en_ks',ierr)
    allocate(hm_sq(ngwf_num,ngwf_num),stat=ierr)
    call utils_alloc_check('bandstructure_core','hm_sq',ierr)
    allocate(dwork(ngwf_num,3),stat=ierr)
    call utils_alloc_check('bandstructure_core','dwork',ierr)
    allocate(zwork(ngwf_num*2),stat=ierr)
    call utils_alloc_check('bandstructure_core','zwork',ierr)
    allocate(bmat(ngwf_num,ngwf_num),stat=ierr)
    call utils_alloc_check('bandstructure_core','bmat',ierr)
    allocate(zov_sq(ngwf_num,ngwf_num),stat=ierr)
    call utils_alloc_check('bandstructure_core','zov_sq',ierr)
    if (basic_unfold) then
       allocate(kpt_ks(ngwf_num,3),stat=ierr)
       call utils_alloc_check('bandstructure_core','kpt_ks',ierr)
       allocate(translate(mdl%nat,3),stat=ierr)
       call utils_alloc_check('bandstructure_core','translate',ierr)
       allocate(flag(mdl%nat),stat=ierr)
       call utils_alloc_check('bandstructure_core','flag',ierr)
       allocate(posfrac(3,mdl%nat),stat=ierr)
       call utils_alloc_check('bandstructure_core','posfrac',ierr)
       allocate(trwfn(ngwf_num),stat=ierr)
       call utils_alloc_check('bandstructure_core','trwfn',ierr)
       max_degeneracy = 9
       allocate(psitrpsi(max_degeneracy,max_degeneracy),stat=ierr)
       call utils_alloc_check('bandstructure_core','psitrpsi',ierr)
       allocate(zlambda(max_degeneracy),stat=ierr)
       call utils_alloc_check('bandstructure_core','zlambda',ierr)
    end if
    if (pub_perturbative_soc) then
       allocate(aux_array(ngwf_basis%num,ngwf_basis%num),stat=ierr)
       call utils_alloc_check('bandstructure_core','aux_array',ierr)
    end if

    ! Prepare translations etc for unfolding
    if (basic_unfold) then
       ! rc2013: assume there is only 1 parallel strategy for the whole system
       do iat=1,mdl%par%nat
          ! ndmh: elements is in input-file order, so use mdl%par%orig_atom
          orig_iat = mdl%par%orig_atom(iat)
          posfrac(1,iat) = mdl%elements(orig_iat)%centre .DOT. mdl%cell%b1
          posfrac(2,iat) = mdl%elements(orig_iat)%centre .DOT. mdl%cell%b2
          posfrac(3,iat) = mdl%elements(orig_iat)%centre .DOT. mdl%cell%b3
       end do
       posfrac = posfrac * recip_twopi
       posfrac = modulo(posfrac,1.0_DP)
       do dim=1,3
          if (pub_bs_unfold(dim) == 1) cycle
          trvec = 0.0_DP
          trvec(dim) = 1.0_DP / pub_bs_unfold(dim)
          flag = .false.
          do iat=1,mdl%par%nat
             orig_iat = mdl%par%orig_atom(iat)
             ispec = mdl%elements(orig_iat)%species_number
             trpos = posfrac(:,iat) + trvec
             wrapped = 1
             if (trpos(dim) >= 1.0_DP) then
                trpos(dim) = trpos(dim) - 1.0_DP
                wrapped = -1
             end if
             found = .false.
             do jat=1,mdl%par%nat
                ! ndmh: elements is in input-file order
                orig_jat = mdl%par%orig_atom(jat)
                if (flag(jat) .or. mdl%elements(orig_jat)%species_number /= ispec) cycle
                dispfrac = posfrac(:,jat) - trpos
                dispabs(1) = dispfrac(1)*mdl%cell%a1%x + &
                     dispfrac(2)*mdl%cell%a2%x + dispfrac(3)*mdl%cell%a3%x
                dispabs(2) = dispfrac(1)*mdl%cell%a1%y + &
                     dispfrac(2)*mdl%cell%a2%y + dispfrac(3)*mdl%cell%a3%y
                dispabs(3) = dispfrac(1)*mdl%cell%a1%z + &
                     dispfrac(2)*mdl%cell%a2%z + dispfrac(3)*mdl%cell%a3%z
                distsq = dispabs(1)*dispabs(1) + dispabs(2)*dispabs(2) + &
                     dispabs(3)*dispabs(3)
                if (distsq < tol) then
                   flag(jat) = .true.
                   translate(iat,dim) = jat * wrapped
                   found = .true.
                   exit
                end if
             end do
             if (.not. found) then
                write(error_message,'(a,i5,a,3(f6.3,a))') &
                     'Error in bandstructure_core: &
                     &no atom equivalent to atom',iat, &
                     ' found when translating by (',trvec(1),',',trvec(2), &
                     ',',trvec(3),')'
                call utils_abort(error_message)
             end if
          end do
       end do
    end if

    ! Prepare perturbative SOC
    if (pub_perturbative_soc) then
       allocate(hamso_k(4, 0:pub_total_num_procs-1), stat=ierr)
       call utils_alloc_check('bandstructure_core','hamso_k',ierr)
       do proc = 0, pub_total_num_procs-1
          do isp = 1, 4
             call sparse_embed_create(hamso_k(isp, proc), ham%ham(1), iscmplx = .true.)
          end do
       end do
    end if

    ! Choose method
    select case (bs_method)

    ! pdh: there is almost certainly more scope for code-sharing below
    ! pdh: and also for modularity that could assist future developments
    ! pdh: but the different approaches can be optimised rather
    ! pdh: differently and one may well be deprecated so for the time
    ! pdh: being separate the code for clarity

    !==================================!
    !                                  !
    !    TIGHT-BINDING STYLE METHOD    !
    !                                  !
    !==================================!

    case ('TB')

       !gcc32: safety check
       if (any( (4.0_DP*ngwf_basis%spheres(:)%radius) .ge. &
          min(geometry_magnitude(mdl%cell%a1), &
          geometry_magnitude(mdl%cell%a2), &
          geometry_magnitude(mdl%cell%a3)) ) ) then
          call utils_abort('ERROR in bandstructure_core: For the TB method to &
               &be correct, the NGWF diameter must be smaller than half of &
               &any lattice vector')
       end if

       ! Diagonalization workspace
       lwork = -1
       if (basic_unfold) then
          call zhegv(1,'V','L',ngwf_num,hm_sq,ngwf_num,zov_sq, &
               ngwf_num,en_ks,zwork,lwork,dwork(1,1),ierr)
       else
          call zhegv(1,'V','L',ngwf_num,hm_sq,ngwf_num,zov_sq, &
               ngwf_num,en_ks,zwork,lwork,dwork(1,1),ierr)
       end if
       if (ierr /= 0) then
          lwork = 2*ngwf_num
       else
          lwork = nint(real(zwork(1),kind=DP))
          deallocate(zwork,stat=ierr)
          call utils_dealloc_check('bandstructure_core','zwork',ierr)
          allocate(zwork(lwork),stat=ierr)
          call utils_alloc_check('bandstructure_core','zwork',ierr)
       end if

       ! Allocate arrays for this method
       ! rc2013: %p is a pointer to the SPAM3 within ham(1)
       hamidxlen = sparse_index_length(ham%ham(1)%p)
       allocate(hamidx(hamidxlen),stat=ierr)
       call utils_alloc_check('bandstructure_core','hamidx',ierr)
       call sparse_generate_index(hamidx,ham%ham(1)%p)
       overidxlen = sparse_index_length(rep%overlap%p)
       allocate(overidx(overidxlen),stat=ierr)
       call utils_alloc_check('bandstructure_core','overidx',ierr)
       call sparse_generate_index(overidx,rep%overlap%p)
       if (.not. pub_perturbative_soc) then
          allocate(ham_ks(0:pub_total_num_procs-1),stat=ierr)
          call utils_alloc_check('bandstructure_core','ham_ks',ierr)
       end if
       allocate(over_ks(0:pub_total_num_procs-1),stat=ierr)
       call utils_alloc_check('bandstructure_core','over_ks',ierr)
       do proc=0,pub_total_num_procs-1
          if (.not. pub_perturbative_soc) then
             call sparse_create(ham_ks(proc),ham%ham(1)%p,iscmplx=.true.)
          end if
          call sparse_create(over_ks(proc),rep%overlap%p,iscmplx=.true.)
       end do
       max_ngwfs_on_atom = maxval(ngwf_basis%num_on_atom)
       allocate(block(max_ngwfs_on_atom,max_ngwfs_on_atom),stat=ierr)
       call utils_alloc_check('bandstructure_core','block',ierr)
       block = (0.0_DP,0.0_DP)

       ! Loop over spin dimensions
       do is = 1, spin_dim

          ! Main loop over unique k-points
          do ikpt=1,nkpts_red,pub_total_num_procs
             do proc=0,pub_total_num_procs-1
                if (ikpt+proc > nkpts_red) exit

                ! Find k-point in original list that matches
                do jkpt=1,nkpts_total
                   if (uni_entry(jkpt) == ikpt+proc) exit
                end do
                if (jkpt > nkpts_total) then
                   call utils_abort('Error in bandstructure_core: &
                        &malformed list of unique k-points')
                end if

                kfrac(1) = kpt(jkpt) .DOT. mdl%cell%a1
                kfrac(2) = kpt(jkpt) .DOT. mdl%cell%a2
                kfrac(3) = kpt(jkpt) .DOT. mdl%cell%a3
                kfrac = kfrac * recip_twopi

                ! Create overlap matrix for this k-point and spin
                call sparse_copy(over_ks(proc),rep%overlap%p)
                do loc_iat=1,mdl%par%num_atoms_on_proc(pub_my_proc_id)
                   iat = loc_iat + mdl%par%first_atom_on_proc(pub_my_proc_id) - 1
                   ! ndmh: elements is in input-file order
                   orig_iat = mdl%par%orig_atom(iat)
                   do jdx=overidx(loc_iat),overidx(loc_iat+1)-1
                      jat = overidx(jdx)
                      orig_jat = mdl%par%orig_atom(jat)
                      atdisp = mdl%elements(orig_iat)%centre - &
                           mdl%elements(orig_jat)%centre
                      dispfrac(1) = atdisp .DOT. mdl%cell%b1
                      dispfrac(2) = atdisp .DOT. mdl%cell%b2
                      dispfrac(3) = atdisp .DOT. mdl%cell%b3
                      dispfrac = dispfrac * recip_twopi
                      zphase = (1.0_DP,0.0_DP)
                      do dim=1,3
                         if (dispfrac(dim) > 0.5_DP) &
                              zphase = zphase * exp(i2pi * kfrac(dim))
                         if (dispfrac(dim) < -0.5_DP) &
                              zphase = zphase * exp(-i2pi * kfrac(dim))
                      end do

                      if (zphase /= (1.0_DP,0.0_DP)) then
                         call sparse_get_block(block,over_ks(proc),jat,iat)
                         block = block * zphase
                         call sparse_put_block(block,over_ks(proc),jat,iat)
                      end if
                   end do
                end do
                if (pub_perturbative_soc) then
                   call sparse_convert(aux_array(:,:),over_ks(proc))
                   if (proc == pub_my_proc_id) then
                      zov_sq = 0.0_DP
                      zov_sq(1:ngwf_basis%num,1:ngwf_basis%num) = &
                                                          aux_array(:,:)
                      zov_sq(ngwf_basis%num+1:ngwf_num, &
                             ngwf_basis%num+1:ngwf_num) = aux_array(:,:)
                   end if
                else
                   call sparse_convert(bmat,over_ks(proc))
                   if (proc == pub_my_proc_id) zov_sq = bmat
                end if

                ! Create Hamiltonian matrix for this k-point and spin
                if (pub_perturbative_soc) then
                   call hamiltonian_soc_matrices(hamso_k(:,proc),ham,rhoij,rep,mdl)
                else
                   call sparse_copy(ham_ks(proc),ham%ham(is)%p)
                end if

                do loc_iat=1,mdl%par%num_atoms_on_proc(pub_my_proc_id)
                   iat = loc_iat + mdl%par%first_atom_on_proc(pub_my_proc_id) - 1
                   ! ndmh: elements is in input-file order
                   orig_iat = mdl%par%orig_atom(iat)
                   do jdx=hamidx(loc_iat),hamidx(loc_iat+1)-1
                      jat = hamidx(jdx)
                      orig_jat = mdl%par%orig_atom(jat)
                      atdisp = mdl%elements(orig_iat)%centre - &
                           mdl%elements(orig_jat)%centre
                      dispfrac(1) = atdisp .DOT. mdl%cell%b1
                      dispfrac(2) = atdisp .DOT. mdl%cell%b2
                      dispfrac(3) = atdisp .DOT. mdl%cell%b3
                      dispfrac = dispfrac * recip_twopi
                      zphase = (1.0_DP,0.0_DP)
                      do dim=1,3
                         if (dispfrac(dim) > 0.5_DP) &
                              zphase = zphase * exp(i2pi * kfrac(dim))
                         if (dispfrac(dim) < -0.5_DP) &
                              zphase = zphase * exp(-i2pi * kfrac(dim))
                      end do
                      if (zphase /= (1.0_DP,0.0_DP)) then
                         if (pub_perturbative_soc) then
                            do isp = 1, 4
                               call sparse_get_block(block,hamso_k(isp,proc)%p,jat,iat)
                               block = block * zphase
                               call sparse_put_block(block,hamso_k(isp,proc)%p,jat,iat)
                            end do
                         else
                            call sparse_get_block(block,ham_ks(proc),jat,iat)
                            block = block * zphase
                            call sparse_put_block(block,ham_ks(proc),jat,iat)
                         end if
                      end if
                   end do
                end do
                if (pub_perturbative_soc) then
                      !     ( |up>  <up|  |up>  <down| )
                      ! H ~ (                          ) ~
                      !     ( |down><up|  |down><down| )
                      !     ( L_z |up><up|  L_- S_+          )   ( #1 #4 )
                      !   ~ (                                ) ~ (       )
                      !     ( L_+ S_-      -L_z |down><down| )   ( #3 #2 )
                   call sparse_convert(aux_array(:,:),hamso_k(1,proc)%p)
                   if (proc == pub_my_proc_id) then
                      hm_sq(1:ngwf_basis%num, &
                            1:ngwf_basis%num) = aux_array(:,:)
                   end if
                   call sparse_convert(aux_array(:,:),hamso_k(3,proc)%p)
                   if (proc == pub_my_proc_id) then
                      hm_sq(ngwf_basis%num+1:ngwf_num, &
                            1:ngwf_basis%num) = aux_array(:,:)
                   end if
                   call sparse_convert(aux_array(:,:),hamso_k(4,proc)%p)
                   if (proc == pub_my_proc_id) then
                      hm_sq(1:ngwf_basis%num, &
                            ngwf_basis%num+1:ngwf_num) = aux_array(:,:)
                   end if
                   call sparse_convert(aux_array(:,:),hamso_k(2,proc)%p)
                   if (proc == pub_my_proc_id) then
                      hm_sq(ngwf_basis%num+1:ngwf_num, &
                            ngwf_basis%num+1:ngwf_num) = aux_array(:,:)
                   end if
                else
                   call sparse_convert(bmat,ham_ks(proc))
                   if (proc == pub_my_proc_id) hm_sq = bmat
                end if

             end do

             ! Diagonalise different k-points on different procs
             if (ikpt+pub_my_proc_id <= nkpts_red) then

                ! Diagonalise Hamiltonian
                if (basic_unfold) then
                   bmat = zov_sq
                   call zhegv(1,'V','L',ngwf_num,hm_sq,ngwf_num, &
                        bmat,ngwf_num,en_ks,zwork,lwork,dwork(1,1),ierr)
                else
                   call zhegv(1,'V','L',ngwf_num,hm_sq,ngwf_num, &
                        zov_sq,ngwf_num,en_ks,zwork,lwork,dwork(1,1),ierr)
                end if
                if (ierr /= 0) then
                   call utils_abort('Error in bandstructure_core: &
                        &zhegv failed with code ',ierr)
                end if

                if (basic_unfold) then

                   ! Recalculate fractional coordinates of relevant k-point
                   do jkpt=1,nkpts_total
                      if (uni_entry(jkpt) == ikpt+pub_my_proc_id) exit
                   end do
                   kfrac(1) = kpt(jkpt) .DOT. mdl%cell%a1
                   kfrac(2) = kpt(jkpt) .DOT. mdl%cell%a2
                   kfrac(3) = kpt(jkpt) .DOT. mdl%cell%a3
                   kfrac = kfrac * recip_twopi

                   call internal_unfold

                end if

             end if

             ! Communicate results between procs
             call comms_barrier

             do proc=0,pub_total_num_procs-1
                if (ikpt+proc > nkpts_red) exit

                if (basic_unfold) then
                   if (proc == pub_my_proc_id) dwork = kpt_ks
                   call comms_bcast(proc,dwork(1,1),3*ngwf_num)
                   do jkpt=1,nkpts_total
                      if (uni_entry(jkpt) == ikpt+proc) &
                           kunfold(:,:,jkpt,is) = dwork
                   end do
                end if

                ! Save eigenvectors to distributed dense matrices,
                ! if requested by passing in dense_eigst
                if (present(dense_eigst)) then
                   if (proc == pub_my_proc_id) bmat(:,:) = hm_sq(:,:)
                   call comms_bcast(proc,bmat,ngwf_num*ngwf_num)
                   do icol = 1, ngwf_num
                      call dense_put_col(bmat(:,icol), &
                           dense_eigst(ikpt+proc,is), icol)
                   end do
                end if

             end do

             ! gcc32: save eigenstates/eigenvalues on this proc/kpt
             if (ikpt+pub_my_proc_id <= nkpts_red) then
                ! Save eigenvectors to full-square array
                if (present(eigst_kpts)) then
                    eigst_kpts(:,:,ikpt+pub_my_proc_id,is) = hm_sq(:,:)
                end if
                eigvl_kpts(:,ikpt+pub_my_proc_id,is) = en_ks(:)
             end if

             call comms_barrier

          end do    ! unique k-points

          call comms_barrier

          ! gcc32: reduce eigenstates/eigenvals over procs for each spin
          if (present(eigst_kpts)) then
             call comms_reduce('SUM', eigst_kpts(:,:,:,is))
          end if
          call comms_reduce('SUM', eigvl_kpts(:,:,is))

       end do   ! spins

       call comms_barrier

       ! Deallocate arrays for this method
       deallocate(block,stat=ierr)
       call utils_dealloc_check('bandstructure_core','block',ierr)
       do proc=pub_total_num_procs-1,0,-1
           call sparse_destroy(over_ks(proc))
           if (.not. pub_perturbative_soc) call sparse_destroy(ham_ks(proc))
       end do
       deallocate(over_ks,stat=ierr)
       call utils_dealloc_check('bandstructure_core','over_ks',ierr)
       if (.not. pub_perturbative_soc) then
          deallocate(ham_ks,stat=ierr)
          call utils_dealloc_check('bandstructure_core','ham_ks',ierr)
       end if
       deallocate(overidx,stat=ierr)
       call utils_dealloc_check('bandstructure_core','overidx',ierr)
       deallocate(hamidx,stat=ierr)
       call utils_dealloc_check('bandstructure_core','hamidx',ierr)

    !==================================!
    !                                  !
    !         K.P STYLE METHOD         !
    !                                  !
    !==================================!

    case ('KP')

       call utils_assert(.not. pub_perturbative_soc, 'ERROR: K.P method&
            & not yet available when using perturbative Spin-Orbit Couplings.')

       ! Diagonalization workspace
       lwork = -1
       if (basic_unfold) then
          call zheev('V','L',ngwf_num,hm_sq,ngwf_num,en_ks, &
               zwork,lwork,dwork(1,1),ierr)
       else
          call zheev('N','L',ngwf_num,hm_sq,ngwf_num,en_ks, &
               zwork,lwork,dwork(1,1),ierr)
       end if
       if (ierr /= 0) then
          lwork = 2*ngwf_num
       else
          lwork = nint(real(zwork(1),kind=DP))
          deallocate(zwork,stat=ierr)
          call utils_dealloc_check('bandstructure_core','zwork',ierr)
          allocate(zwork(lwork),stat=ierr)
          call utils_alloc_check('bandstructure_core','zwork',ierr)
       end if

       ! Prepare k-independent part of Hamiltonian
       allocate(ham0(pub_num_spins),stat=ierr)
       call utils_alloc_check('bandstructure_core','ham0',ierr)
       do is=1,pub_num_spins
          call sparse_create(ham0(is),ham%ham(is)%p,iscmplx=.true.)
          call sparse_copy(ham0(is),ham%ham(is)%p)
          if (pub_any_nl_proj) call sparse_axpy(ham0(is),rep%nonlocpot%p,-1.0_DP)
       end do

       ! Calculate grad matrix elements
       do dim=1,3
          call sparse_create(grad(dim),rep%overlap%p)
       end do
       call integrals_grad(grad, rep%ngwfs_on_grid(1), ngwf_basis, &
            rep%ngwfs_on_grid(1), ngwf_basis, mdl%cell, mdl%fftbox)

       ! Initialise nonlocal matrix and Hamiltonian
       if (pub_any_nl_proj) then
          allocate(nl_k(0:pub_total_num_procs-1),stat=ierr)
          call utils_alloc_check('bandstructure_core','nl_k',ierr)
       end if
       allocate(ham_ks(0:pub_total_num_procs-1),stat=ierr)
       call utils_alloc_check('bandstructure_core','ham_ks',ierr)
       do proc=0,pub_total_num_procs-1
          if (pub_any_nl_proj) &
               call sparse_embed_create(nl_k(proc),rep%nonlocpot,iscmplx=.true.)
          call sparse_create(ham_ks(proc),ham%ham(1)%p,iscmplx=.true.)
       end do

       ! Convert overlap matrix to dense format and obtain
       ! Cholesky decomposition
       allocate(ov_sq(ngwf_basis%num,ngwf_basis%num),stat=ierr)
       call utils_alloc_check('bandstructure_core','ov_sq',ierr)
       call sparse_convert(ov_sq,rep%overlap%p)
       zov_sq = ov_sq

       call dpotrf('L',ngwf_basis%num,ov_sq,ngwf_basis%num,ierr)
       if (ierr /= 0) then
          call utils_abort('Error in bandstructure_core: &
               &dpotrf failed with code ',ierr)
       end if

       ! ndmh: create temporary matrices for real/complex sp overlaps
       if (pub_any_nl_proj) then
          call sparse_embed_create(sp_overlap,rep%sp_overlap,iscmplx=.true.)
       end if

       ! Loop over spins
       do is=1,spin_dim

          ! Main loop over unique k-points
          do ikpt=1,nkpts_red,pub_total_num_procs
             do proc=0,pub_total_num_procs-1
                if (ikpt+proc > nkpts_red) exit

                ! Find k-point in original list that matches
                do jkpt=1,nkpts_total
                   if (uni_entry(jkpt) == ikpt+proc) exit
                end do
                if (jkpt > nkpts_total) then
                   call utils_abort('Error in bandstructure_core: &
                        &malformed list of unique k-points')
                end if

                kcart(1) = kpt(jkpt)%x ; kcart(2) = kpt(jkpt)%y
                kcart(3) = kpt(jkpt)%z

                if (proc == pub_my_proc_id) then
                   kfrac(1) = kpt(jkpt) .DOT. mdl%cell%a1
                   kfrac(2) = kpt(jkpt) .DOT. mdl%cell%a2
                   kfrac(3) = kpt(jkpt) .DOT. mdl%cell%a3
                   kfrac = kfrac * recip_twopi
                end if

                ! Create Hamiltonian at this k-point and spin
                call sparse_copy(ham_ks(proc),ham0(is))

                ! Calculate the nonlocal matrix for this k-point
                if (pub_any_nl_proj) then

                   ! calculate <ngwf|projector> overlap matrix at this kpt
                   call projectors_func_ovlp_box(sp_overlap%p, &
                        rep%ngwfs_on_grid(1),ngwf_basis,proj_basis, nl_projectors, &
                        mdl%fftbox,mdl%cell,kshift=kpt(jkpt))

                   ! ndmh: calculate nonlocal pseudopotential matrix at this kpt
                   if (pub_any_nl_proj) then
                      ! rc2013: create dij to hold nonlocal terms
                      dij%structure = 'E'//trim(rep%postfix)
                      call sparse_embed_create(dij)
                      call pseudo_get_dij(dij%p, mdl%pseudo_sp)
                      call pseudopotentials_nonlocal_mat(nl_k(proc),sp_overlap,&
                           dij)
                      call sparse_embed_destroy(dij)
                   else if (pub_aug) then
                      call sparse_embed_create(ham_embed_array(1),nl_k(proc))
                      call sparse_embed_copy(ham_embed_array(1),nl_k(proc))
                      call sparse_embed_extract_from_array(ham_array,&
                           ham_embed_array)
                      call sparse_embed_extract_from_array(aug_array,ham%dijhat)
                      call sparse_embed_extract_from_array(kern_array,rhoij)
                      call aug_nonlocal_mat(ham_array,aug_array,kern_array,&
                           sp_overlap%p,mdl%regions(1)%pseudo_sp, &
                           mdl%regions(1)%paw_sp)
                      call sparse_embed_destroy_extracted_array(ham_array,&
                           ham_embed_array,.true.)
                      call sparse_embed_copy(nl_k(proc),ham_embed_array(1))
                      call sparse_embed_destroy(ham_embed_array(1))
                      call sparse_embed_destroy_extracted_array(aug_array)
                      call sparse_embed_destroy_extracted_array(kern_array)
                   end if

                   call sparse_axpy(ham_ks(proc),nl_k(proc)%p,1.0_DP)

                end if
                do dim=1,3
                   call sparse_axpy(ham_ks(proc),grad(dim), &
                        cmplx(0.0_DP,-kcart(dim),kind=DP))
                end do
                call sparse_axpy(ham_ks(proc),rep%overlap%p, &
                     0.5_DP*(kcart(1)*kcart(1) + kcart(2)*kcart(2) + &
                     kcart(3)*kcart(3)))
                call sparse_convert(bmat,ham_ks(proc))

                if (proc == pub_my_proc_id) hm_sq = bmat

             end do

             ! Diagonalise different k-points on different procs
             if (ikpt+pub_my_proc_id <= nkpts_red) then

                ! Copy overlap Cholesky factorization
                do ingwf=1,ngwf_basis%num
                   do jngwf=ingwf,ngwf_basis%num
                      bmat(jngwf,ingwf) = &
                           cmplx(ov_sq(jngwf,ingwf),0.0_DP,kind=DP)
                   end do
                end do

                ! Convert Hamiltonian to orthogonal form
                call zhegst(1,'L',ngwf_basis%num,hm_sq,ngwf_basis%num,bmat, &
                     ngwf_basis%num,ierr)
                if (ierr /= 0) then
                   call utils_abort('Error in bandstructure_core: &
                        &zhegst failed with code ',ierr)
                end if

                ! Diagonalise Hamiltonian
                if (basic_unfold) then
                   call zheev('V','L',ngwf_basis%num,hm_sq,ngwf_basis%num,en_ks, &
                        zwork,lwork,dwork(1,1),ierr)
                else
                   call zheev('N','L',ngwf_basis%num,hm_sq,ngwf_basis%num,en_ks, &
                        zwork,lwork,dwork(1,1),ierr)
                end if
                if (ierr /= 0) then
                   call utils_abort('Error in bandstructure_core: &
                        &zheev failed with code ',ierr)
                end if

                if (basic_unfold) then

                   ! Recalculate fractional coordinates of relevant k-point
                   do jkpt=1,nkpts_total
                      if (uni_entry(jkpt) == ikpt+pub_my_proc_id) exit
                   end do
                   kfrac(1) = kpt(jkpt) .DOT. mdl%cell%a1
                   kfrac(2) = kpt(jkpt) .DOT. mdl%cell%a2
                   kfrac(3) = kpt(jkpt) .DOT. mdl%cell%a3
                   kfrac = kfrac * recip_twopi

                   ! Recover eigenvectors
                   call ztrsm('L','L','T','N', ngwf_basis%num, ngwf_basis%num, &
                        cmplx_1, bmat, ngwf_basis%num, hm_sq, ngwf_basis%num)

                   call internal_unfold

                end if
             end if

             ! Communicate results between procs
             call comms_barrier

             do proc=0,pub_total_num_procs-1
                if (ikpt+proc > nkpts_red) exit

                if (basic_unfold) then
                   if (proc == pub_my_proc_id) dwork = kpt_ks
                   call comms_bcast(proc,dwork)
                   do jkpt=1,nkpts_total
                      if (uni_entry(jkpt) == ikpt+proc) &
                           kunfold(:,:,jkpt,is) = dwork
                   end do
                end if

             end do

             ! gcc32: save eigenvalues and eigenstates on this proc/kpt
             if (present(eigst_kpts)) then
                 eigst_kpts(:,:,ikpt+pub_my_proc_id,is) = hm_sq(:,:)
             end if
             if (present(dense_eigst)) then
                do icol = 1, ngwf_num
                   call dense_put_col(hm_sq(:,icol), &
                        dense_eigst(ikpt+pub_my_proc_id,is), icol)
                end do
             end if
             eigvl_kpts(:,ikpt+pub_my_proc_id,is) = en_ks(:)

          end do    ! unique k-points

          ! gcc32: reduce eigenstates/eigenvals over procs for each spin
          if (present(eigst_kpts)) then
             call comms_reduce('SUM', eigst_kpts(:,:,:,is))
          end if
          call comms_reduce('SUM', eigvl_kpts(:,:,is))

       end do   ! spins

       ! Deallocate arrays for this method
       if (pub_any_nl_proj) then
          call sparse_embed_destroy(sp_overlap)
       end if
       deallocate(ov_sq,stat=ierr)
       call utils_dealloc_check('bandstructure_core','ov_sq',ierr)
       do proc=pub_total_num_procs-1,0,-1
          call sparse_destroy(ham_ks(proc))
          if (pub_any_nl_proj) call sparse_embed_destroy(nl_k(proc))
       end do
       deallocate(ham_ks,stat=ierr)
       call utils_dealloc_check('bandstructure_core','ham_ks',ierr)
       if (pub_any_nl_proj) then
          deallocate(nl_k,stat=ierr)
          call utils_dealloc_check('bandstructure_core','nl_k',ierr)
       end if
       do dim=3,1,-1
          call sparse_destroy(grad(dim))
       end do
       do is = pub_num_spins, 1, -1
          call sparse_destroy(ham0(is))
       end do
       deallocate(ham0,stat=ierr)
       call utils_dealloc_check('bandstructure_core','ham0',ierr)

    end select ! bandstructure method

    call comms_barrier

    ! Cleanup perturbative SOC
    if (pub_perturbative_soc) then
       do proc = pub_total_num_procs-1, 0, -1
          do isp = 4, 1, -1
             call sparse_embed_destroy(hamso_k(isp, proc))
          end do
       end do
    end if

    ! Destroy workspace
    if (pub_perturbative_soc) then
       deallocate(aux_array, stat=ierr)
       call utils_dealloc_check('bandstructure_core', 'aux_array', ierr)
    end if

    if (basic_unfold) then
       deallocate(zlambda,stat=ierr)
       call utils_dealloc_check('bandstructure_core','zlambda',ierr)
       deallocate(psitrpsi,stat=ierr)
       call utils_dealloc_check('bandstructure_core','psitrpsi',ierr)
       deallocate(trwfn,stat=ierr)
       call utils_dealloc_check('bandstructure_core','trwfn',ierr)
       deallocate(posfrac,stat=ierr)
       call utils_dealloc_check('bandstructure_core','posfrac',ierr)
       deallocate(flag,stat=ierr)
       call utils_dealloc_check('bandstructure_core','flag',ierr)
       deallocate(translate,stat=ierr)
       call utils_dealloc_check('bandstructure_core','translate',ierr)
       deallocate(kpt_ks,stat=ierr)
       call utils_dealloc_check('bandstructure_core','kpt_ks',ierr)
    end if

    deallocate(zov_sq,stat=ierr)
    call utils_dealloc_check('bandstructure_core','zov_sq',ierr)
    deallocate(bmat,stat=ierr)
    call utils_dealloc_check('bandstructure_core','bmat',ierr)
    deallocate(zwork,stat=ierr)
    call utils_dealloc_check('bandstructure_core','zwork',ierr)
    deallocate(dwork,stat=ierr)
    call utils_dealloc_check('bandstructure_core','dwork',ierr)
    deallocate(hm_sq,stat=ierr)
    call utils_dealloc_check('bandstructure_core','hm_sq',ierr)
    deallocate(en_ks,stat=ierr)
    call utils_dealloc_check('bandstructure_core','en_ks',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving bandstructure_core'

  contains

    subroutine internal_unfold

      implicit none

      ! Local variables
      integer :: iwfn,jwfn,kwfn             ! Wavefunction counters
      integer :: degen                      ! Degeneracy
      real(kind=DP), parameter :: DEGEN_TOL = 1.0e-6_DP
      complex(kind=DP) :: dotp              ! Dot product
      complex(kind=DP) :: vl(1), vr(1)

      ! LAPACK subroutines
      external :: zgeev

      call utils_assert(.not. pub_perturbative_soc, 'ERROR: internal_unfold&
           & not yet available when using perturbative Spin-Orbit Couplings.')

      call utils_assert(.not. pub_perturbative_soc, 'ERROR: internal_unfold&
           & not yet available when using perturbative Spin-Orbit Couplings.')

      ! Loop over eigenstates
      iwfn = 0
      degen = 1
      do
         iwfn = iwfn + degen
         if (iwfn > ngwf_num) exit

         ! Count degeneracy of this state
         do jwfn=iwfn+1,ngwf_num
            if (en_ks(jwfn) - en_ks(iwfn) > DEGEN_TOL) exit
         end do
         degen = jwfn - iwfn

         ! Reallocate memory if necessary
         if (degen > max_degeneracy) then
            deallocate(zlambda,stat=ierr)
            call utils_dealloc_check( &
                 'internal_unfold (bandstructure_calculate)','zlambda',ierr)
            deallocate(psitrpsi,stat=ierr)
            call utils_dealloc_check( &
                 'internal_unfold (bandstructure_calculate)','psitrpsi',ierr)
            allocate(psitrpsi(degen,degen),stat=ierr)
            call utils_alloc_check( &
                 'internal_unfold (bandstructure_calculate)','psitrpsi',ierr)
            allocate(zlambda(degen),stat=ierr)
            call utils_alloc_check( &
                 'internal_unfold (bandstructure_calculate)','zlambda',ierr)
            max_degeneracy = degen
         end if

         ! Attempt displacements along each direction
         do dim=1,3

            if (pub_bs_unfold(dim) == 1) then
               kpt_ks(iwfn:iwfn+degen-1,dim) = 0.0_DP
            else

               ! Loop over degenerate states
               do jwfn=1,degen

                  ! Apply translation operator
                  do ingwf=1,ngwf_num
                     iat = ngwf_basis%atom_of_func(mod(ingwf-1,ngwf_basis%num)+1)
                     jat = translate(iat,dim)
                     wrapped = sign(1,jat)
                     jat = abs(jat)
                     jngwf = ngwf_basis%first_on_atom(jat) + ingwf - &
                          ngwf_basis%first_on_atom(iat)
                     trwfn(jngwf) = hm_sq(ingwf,iwfn+jwfn-1)
                     if (wrapped == -1) trwfn(jngwf) = trwfn(jngwf) * &
                          exp(i2pi*kfrac(dim))
                  end do

                  ! Calculate expectation value (inner product)
                  do kwfn=1,degen
                     psitrpsi(kwfn,jwfn) = (0.0_DP,0.0_DP)
                     do jngwf=1,ngwf_num
                        dotp = (0.0_DP,0.0_DP)
                        do ingwf=1,ngwf_num
                           dotp = dotp + conjg(hm_sq(ingwf,iwfn+kwfn-1)) * &
                                zov_sq(ingwf,jngwf)
                        end do
                        psitrpsi(kwfn,jwfn) = psitrpsi(kwfn,jwfn) + &
                             dotp * trwfn(jngwf)
                     end do
                  end do

               end do

               ! Calculate complex phases
               if (degen == 1) then
                  zlambda(1) = psitrpsi(1,1)
               else
                  call zgeev('N','N',degen,psitrpsi(1,1),max_degeneracy, &
                       zlambda, vl, 1, vr, 1, zwork, lwork, dwork(1,1), ierr)
                  if (ierr /= 0) then
                     call utils_abort('Error in internal_unfold &
                          &(bandstructure_calculate): zgeev failed with code ',&
                          ierr)
                  end if
               end if

               ! Calculate unfolded k-points
               do jwfn=1,degen
                  phase = atan2(aimag(zlambda(jwfn)), &
                       real(zlambda(jwfn),kind=DP))
                  kpt_ks(iwfn+jwfn-1,dim) = -phase * recip_twopi
               end do

            end if

         end do

      end do

    end subroutine internal_unfold

  end subroutine bandstructure_core

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine bsunfold_calculate(ham, rep, ngwf_basis, proj_basis, &
       nl_projectors, rhoij, mdl)

    !======================================================================!
    ! This subroutine unfolds the bandstructure from the supercell         !
    ! representation to the primitive-cell representation                  !
    !----------------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in the Summer of 2017              !
    !----------------------------------------------------------------------!
    ! Perturbative spin-orbit couplings added by JM Escartin, Summer 2016. !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,     !
    ! August 2018                                                          !
    !======================================================================!

    use comms, only: pub_on_root, pub_total_num_procs, pub_my_proc_id, &
         comms_bcast, comms_barrier, pub_root_proc_id, comms_reduce
    use constants, only: DP, cmplx_i, stdout, ANGSTROM, &
         HARTREE_IN_EVS
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
         dense_get_element, dense_get_col
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(.DOT.), operator(-), operator(*), &
         geometry_magnitude, operator(+)
    use linalg, only: linalg_invert_serial
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_bsunfld_kpoint_path_end, pub_bsunfld_num_eigenvalues,&
         pub_bsunfld_kpoint_path_length, pub_bsunfld_kpoint_path_start, &
         pub_debug_on_root, pub_num_spins, pub_rootname, &
         pub_perturbative_soc, PUB_1K, pub_bsunfld_num_atoms_prim, &
         pub_bsunfld_num_kpts_path, pub_bsunfld_transfo, pub_bsunfld_ngroups, &
         pub_bsunfld_groups, pub_bsunfld_nprojatoms, pub_bsunfld_projatoms, &
         pub_bsunfld_restart, pub_print_qc
!$  use rundat, only: pub_threads_max
    use sparse_embed, only: SPAM3_EMBED
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_abort, utils_banner, utils_qc_print

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(in) :: ham              ! Hamiltonian matrices
    type(NGWF_REP), intent(in) :: rep              ! NGWF Representation
    type(FUNC_BASIS), intent(in) :: ngwf_basis     ! NGWF basis type
    type(FUNC_BASIS), intent(in) :: proj_basis     ! Projector basis type
    type(PROJECTOR_SET), intent(inout) :: nl_projectors
    type(SPAM3_EMBED), intent(in) :: rhoij(pub_num_spins) ! Projector denskern
    type(MODEL), intent(in) :: mdl

    ! Local variables
    type(POINT) :: kstart                 ! k-point at start of path segment
    type(POINT) :: kend                   ! k-point at end of path segment
    type(POINT) :: kseg                   ! path segment
    type(POINT) :: kptdiff                ! difference between k-points
    type(POINT), allocatable :: kpt(:)    ! primitive-cell k-point
    type(POINT), allocatable :: SCkpt(:)  ! supercell k-point
    integer :: ierr                       ! Error flag
    integer :: ipath                      ! Bandstructure path counter
    integer :: nsegpts                    ! Number of points along path segment
    integer :: isegpt                     ! Path segment point counter
    integer :: nkpts                      ! Total number of kpts
    integer :: ikpt,jkpt,ikpt_red,jkpt_red  ! k-point counter
    integer :: ival, icell, icell_inuse
    integer :: ingwf, ingwf_mod, jngwf, jngwf_mod, loc_ingwf  ! NGWF counters
    integer :: nlower,nupper              ! Limits on NGWF bounds for output
    integer :: neigvals                   ! Actual number of eigenvalues used
    integer :: nkpts_red                  ! Reduced list of unique k-points
    integer :: nelec                      ! Number of electrons
    integer :: spin_dim                   ! Number of independent spin dims
    integer :: ngwf_num                   ! Number of NGWFs
    integer, allocatable :: uni_entry(:)  ! Pointer to unique entry (prim cell)
    integer, allocatable :: SCuni_entry(:)  ! Pointer to unique entry
    integer, allocatable :: revuni_entry(:,:) ! inverse of uni_entry
    logical :: unique                     ! Flag for reducing k-points
    complex(kind=DP), allocatable :: spect_func_red(:,:,:)
    real(kind=DP), allocatable :: eigvals_func_red(:,:)
    real(kind=DP), allocatable :: eigvl_kpts(:,:,:)  ! SC eigenvalues (red)
    complex(kind=DP), allocatable :: eigst_kpts(:,:) ! SC eigenstates (red)
    complex(kind=DP), allocatable :: buffer_col(:)
    type(DEM), allocatable :: dense_eigst(:,:)
    real(kind=DP) :: transfo_mat(3,3), inv_transfo_mat(3,3), tr_transfo_mat(3,3)
    integer :: iter1, iter2, iter3
    type(POINT) :: a1_prim, a2_prim, a3_prim, b1_prim, b2_prim, b3_prim

    real(kind=DP) :: temp(8), dotprod
    complex(kind=DP), allocatable :: exp_idotprod(:)
    integer, allocatable :: ingwf_eqv(:)
    real(kind=DP), allocatable :: coeff_SC(:,:)
    integer, allocatable :: ngwfs_idx_use(:), atoms_idx_use(:)
    integer, allocatable :: ngwfs_idx_use_cells(:) ! prim cells of useful ngwfs
    integer, allocatable :: ngwfidx_wrt_cell(:)
    integer :: num_atoms_use, num_ngwfs_use_init, num_ngwfs_use_final, &
         ngwf_counter, ngwf_use_counter, atom_use_counter, num_temp
    integer :: idx, jdx

    integer :: max_unientries, num_unientries, ceil
    type(POINT), allocatable :: pos_cells_use(:)
    integer :: iat_ext, iter, dir, num_ngwfs_cell
    character(len=256) :: specfunc_file(pub_bsunfld_nprojatoms + 1), &
        specfuncred_file(pub_bsunfld_nprojatoms + 1), name_group
    integer :: specfunc_unit, specfuncred_unit
    integer :: SCnkpts_red
    complex(kind=DP) :: element_0, element_1, element_2
    compleX(kind=DP), allocatable :: row_spectral(:,:)
    type(DEM) :: dense_ovlp
    real(kind=DP), allocatable :: overlap_array(:,:)

    integer, allocatable :: list_kptred_done(:)
    integer :: num_kptred_done, ikptred_done
    logical :: kpt_done
    integer, allocatable :: to_include(:,:)

    ! jme: KPOINTS_DANGER
    ! First k-point enforced in all uses of n_occ via PUB_1K.

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering bsunfold_calculate'

    if (pub_on_root) write(stdout,'(/a)') utils_banner('=', &
         'Bandstructure Unfolding Calculation')

    ! gcc32: for now, abort if spin-polarised
    if (pub_num_spins > 1) call utils_abort('Cannot deal with spin-polarised &
         &calculations yet')

    call dense_create(dense_ovlp,ngwf_basis%num,ngwf_basis%num, &
         iscmplx=rep%ngwfs_on_grid(1)%iscmplx)
    call dense_convert(dense_ovlp, rep%overlap)

    allocate(overlap_array(ngwf_basis%num,ngwf_basis%num),stat=ierr)
    call utils_alloc_check('bsunfold_calculate','overlap_array',ierr)

    do ingwf = 1, ngwf_basis%num
       do jngwf = 1, ngwf_basis%num
          call dense_get_element(overlap_array(ingwf,jngwf), dense_ovlp, &
               ingwf, jngwf)
       end do
    end do

    call dense_destroy(dense_ovlp)


    ! recover the transformation matrix
    do iter1 = 1,3
       do iter2 = 1,3
          iter3 = (iter1-1)*3 + iter2
          transfo_mat(iter1, iter2) = pub_bsunfld_transfo(iter3)
       end do
    end do

    inv_transfo_mat = transfo_mat
    call linalg_invert_serial(inv_transfo_mat, ierr)
    if (ierr.ne.0) call utils_abort('Inversion of transformation matrix failed')

    do iter1 = 1,3
       do iter2 = 1,3
          tr_transfo_mat(iter1, iter2) = transfo_mat(iter2, iter1)
       end do
    end do

    ! get the primitive cell real-space and reciprocal lattice vectors
    a1_prim = inv_transfo_mat(1,1) * mdl%cell%a1 + &
         inv_transfo_mat(1,2) * mdl%cell%a2  + &
         inv_transfo_mat(1,3) * mdl%cell%a3
    a2_prim = inv_transfo_mat(2,1) * mdl%cell%a1 + &
         inv_transfo_mat(2,2) * mdl%cell%a2  + &
         inv_transfo_mat(2,3) * mdl%cell%a3
    a3_prim = inv_transfo_mat(3,1) * mdl%cell%a1 + &
         inv_transfo_mat(3,2) * mdl%cell%a2  + &
         inv_transfo_mat(3,3) * mdl%cell%a3

    b1_prim = tr_transfo_mat(1,1) * mdl%cell%b1 + &
         tr_transfo_mat(1,2) * mdl%cell%b2  + &
         tr_transfo_mat(1,3) * mdl%cell%b3
    b2_prim = tr_transfo_mat(2,1) * mdl%cell%b1 + &
         tr_transfo_mat(2,2) * mdl%cell%b2  + &
         tr_transfo_mat(2,3) * mdl%cell%b3
    b3_prim = tr_transfo_mat(3,1) * mdl%cell%b1 + &
         tr_transfo_mat(3,2) * mdl%cell%b2  + &
         tr_transfo_mat(3,3) * mdl%cell%b3

    if (pub_on_root) then
       write(stdout,*) '    Real Lattice                                 Reciprocal Lattice '
       write(stdout,'(2(3f14.8,3x))') a1_prim%x, a1_prim%y, a1_prim%z, b1_prim%x, b1_prim%y, b1_prim%z
       write(stdout,'(2(3f14.8,3x))') a2_prim%x, a2_prim%y, a2_prim%z, b2_prim%x, b2_prim%y, b2_prim%z
       write(stdout,'(2(3f14.8,3x))') a3_prim%x, a3_prim%y, a3_prim%z, b3_prim%x, b3_prim%y, b3_prim%z
    end if

    ! Count the number of primitive-cell k-points needed
    nkpts = pub_bsunfld_num_kpts_path * pub_bsunfld_kpoint_path_length


    if (pub_perturbative_soc) then
       nelec = rep%n_occ(1,PUB_1K)
       select case (pub_num_spins)
          case (1)
             nelec = 2 * nelec
          case (2)
             nelec = nelec + rep%n_occ(2,PUB_1K)
       end select
    else
       nelec = maxval(rep%n_occ(:,PUB_1K))
    end if

    ! ndmh: write all bands if bs_num_eigenvalues is left at default value
    if (pub_bsunfld_num_eigenvalues <= 0) then
       pub_bsunfld_num_eigenvalues = nelec
    else
       ! the number of eigenvalues considered cannot exceed the number of
       ! (nondegenerate) occupied states
       pub_bsunfld_num_eigenvalues = min(pub_bsunfld_num_eigenvalues, nelec)
    end if

    ! Set spin dimensions.
    if (pub_perturbative_soc) then
       spin_dim = 1
       ngwf_num = 2 * ngwf_basis%num
    else
       spin_dim = pub_num_spins
       ngwf_num = ngwf_basis%num
    end if

    nlower = max(nelec-pub_bsunfld_num_eigenvalues+1,1)
    nupper = min(nelec+pub_bsunfld_num_eigenvalues,ngwf_num)

    if (pub_on_root)  write(stdout,'(a,2i8)') 'Eigenvalue range: ',nlower,nupper
    if(pub_print_qc) then
       call utils_qc_print('bsunfold_nlower',nlower)
       call utils_qc_print('bsunfold_nupper',nupper)
    end if
    neigvals = nupper - nlower + 1

    allocate(row_spectral(2*neigvals, &
         1 + pub_bsunfld_nprojatoms), stat=ierr)
    call utils_alloc_check('bsunfold_calculate','row_spectral',ierr)


    ! Compile a complete list of primitive-cell k-points
    allocate(kpt(nkpts),stat=ierr)
    call utils_alloc_check('bsunfold_calculate','kpt',ierr)
    ikpt = 0
    do ipath=1,pub_bsunfld_kpoint_path_length
       kstart = pub_bsunfld_kpoint_path_start(1,ipath) * b1_prim + &
            pub_bsunfld_kpoint_path_start(2,ipath) * b2_prim + &
            pub_bsunfld_kpoint_path_start(3,ipath) * b3_prim
       kend = pub_bsunfld_kpoint_path_end(1,ipath) * b1_prim + &
            pub_bsunfld_kpoint_path_end(2,ipath) * b2_prim + &
            pub_bsunfld_kpoint_path_end(3,ipath) * b3_prim
       kseg = kend - kstart
       nsegpts = max(pub_bsunfld_num_kpts_path-1,1)
       do isegpt=0,nsegpts
          ikpt = ikpt + 1
          kpt(ikpt) = kstart + (isegpt / real(nsegpts,kind=DP)) * kseg
       end do
    end do

    ! Now compile a reduced list of unique prim-cell k-points avoiding
    ! repetition (e.g. at start and end of consecutive segments)
    allocate(uni_entry(nkpts),stat=ierr)
    call utils_alloc_check('bsunfold_calculate','uni_entry',ierr)
    nkpts_red = 0
    do ikpt=1,nkpts
       unique = .true.
       do jkpt=1,ikpt-1
          kptdiff = kpt(ikpt) - kpt(jkpt)
          if (geometry_magnitude(kptdiff) < 1.0e-6_DP) then
             unique = .false.
             uni_entry(ikpt) = uni_entry(jkpt)
             exit
          end if
       end do
       if (unique) then
          nkpts_red = nkpts_red + 1
          uni_entry(ikpt) = nkpts_red
       end if
    end do

    ! gcc32: get maximum number of time a unique entry appears
    max_unientries = 0
    do ikpt = 1, nkpts_red
        num_unientries = 0
        do jkpt = 1, nkpts
           if (uni_entry(jkpt).eq.ikpt) then
              num_unientries = num_unientries + 1
           end if
        end do

        if (num_unientries > max_unientries) max_unientries = num_unientries
    end do

    allocate(spect_func_red(nkpts_red,neigvals, &
         1 + pub_bsunfld_nprojatoms), stat=ierr)
    call utils_alloc_check('bsunfold_calculate','spect_func_red',ierr)
    allocate(eigvals_func_red(nkpts_red,neigvals),stat=ierr)
    call utils_alloc_check('bsunfold_calculate','eigvals_func_red',ierr)

    ! transform primitive cell unique k-points to supercell (SC) kpts
    ! use the transformation matrix
    allocate(SCkpt(nkpts),stat=ierr)
    call utils_alloc_check('bsunfold_calculate','SCkpt',ierr)

    nsegpts = max(pub_bsunfld_num_kpts_path-1,1)
    allocate(coeff_SC(nsegpts+1,3),stat=ierr)
    call utils_alloc_check('bsunfold_calculate','coeff_SC',ierr)

    ! first get the fractional coordinates for the BZ of the SC for the
    ! full list of primitive-cell kpts:
    ikpt = 1
    do ipath=1,pub_bsunfld_kpoint_path_length
       coeff_SC(1,:) = pub_bsunfld_kpoint_path_start(1,ipath) * &
            tr_transfo_mat(1,:) + pub_bsunfld_kpoint_path_start(2,ipath) * &
            tr_transfo_mat(2,:) +  pub_bsunfld_kpoint_path_start(3,ipath) * &
            tr_transfo_mat(3,:)

       coeff_SC(nsegpts+1,:) = pub_bsunfld_kpoint_path_end(1,ipath) * &
            tr_transfo_mat(1,:) + pub_bsunfld_kpoint_path_end(2,ipath) * &
            tr_transfo_mat(2,:) +  pub_bsunfld_kpoint_path_end(3,ipath) * &
            tr_transfo_mat(3,:)

       do iter = 2, nsegpts
           coeff_SC(iter,:) = coeff_SC(1,:) + (iter - 1) * &
                (coeff_SC(nsegpts+1,:) - coeff_SC(1,:)) / nsegpts
       end do

       ! return all points to the first BZ
       do iter = 1, nsegpts+1
          do dir = 1,3
             if (coeff_SC(iter,dir) < -0.5_DP) coeff_SC(iter,dir) = &
                  coeff_SC(iter,dir) + 1.0_DP
             if (coeff_SC(iter,dir) >= 0.5_DP) coeff_SC(iter,dir) = &
                  coeff_SC(iter,dir) - 1.0_DP
          end do

          SCkpt(ikpt) = coeff_SC(iter,1) * mdl%cell%b1 + &
               coeff_SC(iter,2) * mdl%cell%b2 + &
               coeff_SC(iter,3) * mdl%cell%b3

          ikpt = ikpt + 1
       end do

    end do

    deallocate(coeff_SC,stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','coeff_SC',ierr)

    ! Now compile a reduced list of unique SUPERCELL k-points avoiding
    ! repetition (e.g. at start and end of consecutive segments)
    allocate(SCuni_entry(nkpts),stat=ierr)
    call utils_alloc_check('bsunfold_calculate','SCuni_entry',ierr)
    SCnkpts_red = 0
    do ikpt=1,nkpts
       unique = .true.
       do jkpt=1,ikpt-1
          kptdiff = SCkpt(ikpt) - SCkpt(jkpt)
          if (geometry_magnitude(kptdiff) < 1.0e-6_DP) then
             unique = .false.
             SCuni_entry(ikpt) = SCuni_entry(jkpt)
             exit
          end if
       end do
       if (unique) then
          SCnkpts_red = SCnkpts_red + 1
          SCuni_entry(ikpt) = SCnkpts_red
       end if
    end do

    ! Allocate workspace for diagonalisation and unfolding
    allocate(dense_eigst(SCnkpts_red,spin_dim), stat=ierr)
    call utils_alloc_check('bsunfold_calculate','dense_eigst', ierr)
    do ikpt = 1, SCnkpts_red
       do idx = 1, spin_dim
          call dense_create(dense_eigst(ikpt,idx), ngwf_num, ngwf_num, &
               iscmplx = .true.)
       end do
    end do

    allocate(eigst_kpts(ngwf_num,ngwf_num), stat=ierr)
    call utils_alloc_check('bsunfold_calculate','eigst_kpts',ierr)
    allocate(buffer_col(ngwf_num), stat=ierr)
    call utils_alloc_check('bsunfold_calculate','buffer_col',ierr)
    allocate(eigvl_kpts(ngwf_num,SCnkpts_red,spin_dim),stat=ierr)
    call utils_alloc_check('bsunfold_calculate','eigvl_kpts',ierr)
    ! allocate reverse uni entries
    allocate(revuni_entry(nkpts_red, max_unientries), stat=ierr)
    call utils_alloc_check('bsunfold_calculate','revuni_entry',ierr)

    ! actually calculate the reverse uni_entries for the prim and SC:
    revuni_entry = -1 ! initialise it to -1 to avoid false entries

    do ikpt = 1, nkpts_red
       do jkpt = 1, nkpts
          iter = 1
          if (uni_entry(jkpt).eq.ikpt) then
             revuni_entry(ikpt,iter) = jkpt
             iter = iter + 1
          end if
       end do
    end do

    eigst_kpts = cmplx(0.0_DP, 0.0_DP, kind=DP)
    eigvl_kpts = 0.0_DP

    if (pub_on_root) write(stdout,'(/a)') 'Hamiltonian Diagonalisation ...'

    ! get eigenvalues and eigenvectors for selected supercell k-points
    call bandstructure_core(ham, rep, ngwf_basis, proj_basis, nl_projectors, &
         rhoij, mdl, 'TB', spin_dim, ngwf_num, SCnkpts_red, &
         nkpts, SCkpt, SCuni_entry, .false., eigvl_kpts, &
         dense_eigst = dense_eigst)

    if (pub_on_root) write(stdout,'(/a)') 'Hamiltonian Diagonalisation complete'

    call comms_barrier

    if (pub_num_spins.ne.1) call utils_abort('ERROR in bsunfold_calculate: Can &
         &only deal with 1 spin at present')

    ! see how many atoms/ngwfs are in the projected groups:
    if (pub_bsunfld_ngroups.ne.1) call utils_abort('ERROR in bsunfold_calculate: &
         &Can only deal with one group of atoms at the moment!')

    num_atoms_use = 0
    num_ngwfs_use_init = 0
    num_ngwfs_use_final = 0 ! twice num_ngwfs_use_init if SOC is used

    do iat_ext = 1, mdl%nat ! in EXTERNAL (input) ORDER
       if (any(pub_bsunfld_groups(:,1) .eq. &
            mdl%elements(iat_ext)%species_id)) then

          num_atoms_use = num_atoms_use + 1
          num_ngwfs_use_init = num_ngwfs_use_init + &
               ngwf_basis%num_on_atom(mdl%par%distr_atom(iat_ext))
       end if
    end do ! over atoms


    if (pub_perturbative_soc) then
       num_ngwfs_use_final = 2 * num_ngwfs_use_init
    else
       num_ngwfs_use_final = num_ngwfs_use_init
    end if

    if (num_atoms_use.le.0) call utils_abort('No atoms matched the &
         &group description!')
    if (mod(num_atoms_use,pub_bsunfld_num_atoms_prim).ne.0) then
       call utils_abort('Number of projected atoms and unit cells do not match')
    end if

    ! indices of useful atoms (in EXTERNAL ORDER)
    allocate(atoms_idx_use(num_atoms_use), stat=ierr)
    call utils_alloc_check('bsunfold_calculate','atoms_idx_use',ierr)
    ! indices of useful ngwfs (in INTERNAL ORDER)
    allocate(ngwfs_idx_use(num_ngwfs_use_final), stat=ierr)
    call utils_alloc_check('bsunfold_calculate','ngwfs_idx_use',ierr)
    allocate(ngwfs_idx_use_cells(num_ngwfs_use_final), stat=ierr)
    call utils_alloc_check('bsunfold_calculate','ngwfs_idx_use_cells',ierr)
    allocate(ngwfidx_wrt_cell(num_ngwfs_use_final), stat=ierr)
    call utils_alloc_check('bsunfold_calculate','ngwfidx_wrt_cell',ierr)

    allocate(to_include(num_ngwfs_use_final, 1 + pub_bsunfld_nprojatoms), &
         stat = ierr)
    call utils_alloc_check('bsunfold_calculate', 'to_include', ierr)

    ! indices of primitive cell positions (in EXTERNAL ORDER)
    allocate(pos_cells_use(num_atoms_use/pub_bsunfld_num_atoms_prim), stat=ierr)
    call utils_alloc_check('bsunfold_calculate','pos_cells_use',ierr)

    ngwf_use_counter = 1
    atom_use_counter = 1
    to_include(:,:) = 0
    to_include(:,1) = 1 ! the case where you consider all the atoms in
                              ! the projected subsystem

    do iat_ext = 1, mdl%nat ! atom order in EXTERNAL order
       ! if atoms belong to the projection group
       if (any(pub_bsunfld_groups(:,1) .eq. &
            mdl%elements(iat_ext)%species_id)) then

          atoms_idx_use(atom_use_counter) = iat_ext ! External order

          ngwf_counter = ngwf_basis%first_on_atom(mdl%par%distr_atom(iat_ext))

          do iter = 1, ngwf_basis%num_on_atom(mdl%par%distr_atom(iat_ext))
             ngwfs_idx_use(ngwf_use_counter+iter-1) = ngwf_counter+iter-1
             do idx = 2, 1 + pub_bsunfld_nprojatoms
                if (any(pub_bsunfld_projatoms(:,idx-1) .eq. &
                     mdl%elements(iat_ext)%species_id)) then
                   to_include(ngwf_counter+iter-1, idx) = 1
                end if
             end do

             ngwfs_idx_use_cells(ngwf_use_counter+iter-1) = &
                  1 + (atom_use_counter - 1) / pub_bsunfld_num_atoms_prim
          end do

          ngwf_use_counter = ngwf_use_counter + &
               ngwf_basis%num_on_atom(mdl%par%distr_atom(iat_ext))
          atom_use_counter = atom_use_counter + 1
       end if ! if atom in projection group
    end do ! over atoms

    ! if SOC is in use, duplicate the list of ngwfs
    if (pub_perturbative_soc) then
       ngwfs_idx_use(num_ngwfs_use_init+1:num_ngwfs_use_final) = &
            ngwfs_idx_use(1:num_ngwfs_use_init)
       ngwfs_idx_use_cells(num_ngwfs_use_init+1:num_ngwfs_use_final) = &
            ngwfs_idx_use_cells(1:num_ngwfs_use_init)
       do idx = 2, 1+ pub_bsunfld_nprojatoms
          to_include(num_ngwfs_use_init+1:num_ngwfs_use_final,idx) = &
               to_include(1:num_ngwfs_use_init,idx)
       end do
    end if

    ! get positions of unit cells, assuming that the atoms on which one
    ! projects have been grouped by unit cells in the same order
    do iter = 1, num_atoms_use / pub_bsunfld_num_atoms_prim
       iat_ext = 1 + (iter-1) * pub_bsunfld_num_atoms_prim
       pos_cells_use(iter) = mdl%elements(atoms_idx_use(iat_ext))%centre
    end do


    ! get num_ngwfs per cell (wihtout factor any factor of 2,
    ! even when pub_perturbative_soc is involved)
    num_ngwfs_cell = 0
    do atom_use_counter = 1, pub_bsunfld_num_atoms_prim
       iat_ext = atoms_idx_use(atom_use_counter)
       num_ngwfs_cell = num_ngwfs_cell + &
            ngwf_basis%num_on_atom(mdl%par%distr_atom(iat_ext))
    end do

    ngwfidx_wrt_cell = -1 ! initialise
    ngwf_use_counter = 1
    ! go over the used atoms
    do atom_use_counter = 1, num_atoms_use
       iat_ext = atoms_idx_use(atom_use_counter) ! atom in full external order
       ! how many ngwfs this atom has
       num_temp = ngwf_basis%num_on_atom(mdl%par%distr_atom(iat_ext))

       ! reset ngwf counter if we are going into a new cell
       if (mod((atom_use_counter), pub_bsunfld_num_atoms_prim) == 1) then
          ngwf_counter = 1
       end if

       do iter = 1, num_temp
          ngwfidx_wrt_cell(ngwf_use_counter) = ngwf_counter + iter - 1
          ngwf_use_counter = ngwf_use_counter + 1
       end do
       ngwf_counter = ngwf_counter + num_temp
    end do

    if (pub_perturbative_soc) then
       ngwfidx_wrt_cell(num_ngwfs_use_init+1:num_ngwfs_use_final) = &
            ngwfidx_wrt_cell(1:num_ngwfs_use_init)
    end if

    write(specfuncred_file(1),'(2a)') trim(pub_rootname), &
         '_unfolded_specfunc_red.txt'
    write(specfunc_file(1),'(2a)') trim(pub_rootname), &
         '_unfolded_specfunc.txt'

    do idx = 2, pub_bsunfld_nprojatoms + 1
       write(name_group,'(i2)') idx - 1
       name_group = trim(adjustl(name_group))
       write(specfuncred_file(idx),'(2a)') trim(pub_rootname), &
            trim(adjustl('_unfolded_specfunc_red_G'//trim(name_group)//'.txt'))
       write(specfunc_file(idx),'(2a)') trim(pub_rootname), &
            trim(adjustl('_unfolded_specfunc_G'//trim(name_group)//'.txt'))
    end do

    ! the restarts are always written in the "unfolded_specfunc_red.txt", with
    ! "unfolded_specfunc_red.txt" being written only at the end after all
    ! calculations are performed

    spect_func_red = 0.0_DP
    eigvals_func_red = 0.0_DP

    if (pub_on_root) then
       if (pub_bsunfld_restart) then
          do idx = 1, 1 + pub_bsunfld_nprojatoms
             write(stdout,*) 'RESTARTING FROM FILE ', specfuncred_file(idx)
             ! see if one can open the file
             specfuncred_unit = utils_unit()
             open(specfuncred_unit,file=trim(adjustl(specfuncred_file(idx))), &
                  form='FORMATTED',action='READ',iostat=ierr)
             if (ierr /= 0) then
                call utils_abort('Error in bsunfold_calculate: opening "'//&
                     trim(adjustl(specfuncred_file(idx)))//'" failed with&
                     & code ',ierr)
             end if

             ! read stored values until end of file
             do
                read(specfuncred_unit,*,iostat=ierr) temp(1:8)
                if (ierr.ne.0) exit ! exit if error in reading line
                ival = int(temp(4))
                ikpt = int(temp(8))
                ! read eigenvalues in eV, transform in Ha internally
                if (idx .eq. 1) then
                   eigvals_func_red(ikpt,ival) = temp(1) / HARTREE_IN_EVS
                end if
                spect_func_red(ikpt,ival,idx) = cmplx(temp(2),temp(3),kind=DP)
             end do

             close(specfuncred_unit)
          end do
       end if

       ! gcc32: count how many reduced kpts have already been calculated
       num_kptred_done = 0
       do ikpt_red = 1, nkpts_red
          ! if any eigenvalue is non-zero, then it has been read from the
          ! restart file, meaning that the respective kpt should be skipped
          if (any(abs(eigvals_func_red(ikpt_red,:)) > 1.0e-6_DP)) &
              num_kptred_done = num_kptred_done + 1
       end do

    end if

    ! gcc32: there is no need to send the read spectral function values or
    ! eigenvalues to other procs, since only the root deals with printing

    ! gcc32: However, each proc should have a list of the reduced kpts that
    ! should be skipped
    call comms_bcast(pub_root_proc_id, num_kptred_done)
    call comms_barrier

    allocate(list_kptred_done(num_kptred_done), stat=ierr)
    call utils_alloc_check('bsunfold_calculate','list_kptred_done',ierr)

    ! see what kpts should be skipped
    list_kptred_done = 0

    if (pub_on_root.and.pub_bsunfld_restart) then
       ikptred_done = 1
       do ikpt_red = 1, nkpts_red
          if (any(abs(eigvals_func_red(ikpt_red,:)) > 1.0e-6_DP)) then
             list_kptred_done(ikptred_done) = ikpt_red
             ikptred_done = ikptred_done + 1
          end if
       end do
    end if

    ! send list from root to all procs
    call comms_bcast(pub_root_proc_id, list_kptred_done)
    call comms_barrier

    ! get upper limit of reduced kpt index for each proc
    ceil = ceiling(real(nkpts_red,kind=DP)/real(pub_total_num_procs,kind=DP)) *&
         pub_total_num_procs

    ! go over primitive-cell reduced kpts, parallelise over them
    do ikpt_red = 1+pub_my_proc_id, ceil, pub_total_num_procs

       ! if reduced kpt has been dealt with before, skip it
       kpt_done = .false.
       if (any(list_kptred_done(:) .eq. ikpt_red)) kpt_done = .true.

       ! get eigenstates for this kpt and spin:
       do jdx = 1, SCnkpts_red
          do idx = 1, ngwf_num
             buffer_col = 0.0_DP
             call dense_get_col(buffer_col,dense_eigst(jdx,1),idx)

             ! this weirdly placed if statement appears here because
             ! dense_get_col needs to get called on all cores
             if ((ikpt_red <= nkpts_red).and.(.not.kpt_done)) then
                ! get equivalent supercell reduced k-point
                ikpt = revuni_entry(ikpt_red,1) ! we can just chose one of the
                                  ! duplicated values, it does not matter which
                jkpt_red = SCuni_entry(ikpt)
                if (jdx == jkpt_red) eigst_kpts(:,idx) = buffer_col
             end if

          end do
       end do

       if ((ikpt_red <= nkpts_red).and.(.not.kpt_done)) then

          write(stdout,*) 'Working on reduced KPT ', ikpt_red, ' on proc ', &
              pub_my_proc_id

          ! use ikpt and jkpt_red from before

          ! fill corresponding eigenvalue table
          eigvals_func_red(ikpt_red,1:neigvals) = &
               eigvl_kpts(nlower:nupper,jkpt_red,1)

          ! initialise
          row_spectral = 0.0_DP

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(ival,loc_ingwf,ingwf,ingwf_mod,ingwf_eqv,jngwf,jngwf_mod,&
!$OMP      element_0,element_1,element_2, dotprod,exp_idotprod,icell, &
!$OMP      ierr,icell_inuse) &
!$OMP SHARED(neigvals,ngwf_basis,pub_threads_max,&
!$OMP      pos_cells_use,ngwfs_idx_use_cells,ngwfs_idx_use,ikpt_red, &
!$OMP      num_ngwfs_use_final,num_ngwfs_use_init,rep, pub_on_root, &
!$OMP      ngwf_num,ikpt,jkpt_red,num_atoms_use,pub_bsunfld_num_atoms_prim,&
!$OMP      nlower,kpt, overlap_array,ngwfidx_wrt_cell,num_ngwfs_cell, &
!$OMP      to_include, eigst_kpts) &
!$OMP REDUCTION (+:row_spectral)

       allocate(exp_idotprod(num_atoms_use/pub_bsunfld_num_atoms_prim), &
            stat=ierr)
       allocate(ingwf_eqv(num_atoms_use/pub_bsunfld_num_atoms_prim), &
            stat=ierr)

       ! go over eigenvalues
!$OMP DO
          do ival = 1,neigvals

             ! go over useful ngwfs
             do loc_ingwf = 1,num_ngwfs_use_final
                ! get actual ngwf index
                ingwf = ngwfs_idx_use(loc_ingwf)
                if (loc_ingwf.gt.num_ngwfs_use_init) &
                     ingwf = ingwf + ngwf_basis%num

                ingwf_mod = 1 + mod((ingwf-1), ngwf_basis%num)
                icell_inuse = ngwfs_idx_use_cells(loc_ingwf)
                element_0 = eigst_kpts(ingwf,nlower+ival-1)

                ! precalculate exp_idotprod and ingwf_eqv for each
                ! icell for this loc_ingwf
                do icell=1,num_atoms_use/pub_bsunfld_num_atoms_prim
                   dotprod = kpt(ikpt).DOT.(pos_cells_use(icell) - &
                        pos_cells_use(icell_inuse))
                   exp_idotprod(icell) = exp(cmplx_i*dotprod)
                   ingwf_eqv(icell) = ngwfs_idx_use(num_ngwfs_cell * (icell-1) + &
                        ngwfidx_wrt_cell(loc_ingwf))
                end do

                ! go over all ngwfs (internal order)
                do jngwf = 1, ngwf_num

                   jngwf_mod = 1 + mod((jngwf-1), ngwf_basis%num)

                   element_1 = element_0 * &
                        conjg(eigst_kpts(jngwf,nlower+ival-1))

                   element_2 = cmplx(0.0_DP,0.0_DP,kind=DP)

                   do icell = 1, num_atoms_use/pub_bsunfld_num_atoms_prim
                      element_2 = element_2 + exp_idotprod(icell) * &
                           overlap_array(jngwf_mod,ingwf_eqv(icell))
                   end do ! over cell

                   row_spectral(ival,:) = row_spectral(ival,:) + &
                        element_2 * element_1 * to_include(loc_ingwf,:)

                end do ! all ngwfs
             end do ! useful ngwfs
          end do ! eigenvals
!$OMP END DO
       deallocate(exp_idotprod,stat=ierr)
       deallocate(ingwf_eqv,stat=ierr)

!$OMP END PARALLEL

          ! normalise by number of cells and number of spins
          do idx = 1, 1 + pub_bsunfld_nprojatoms
             spect_func_red(ikpt_red,:,idx) = row_spectral(:,idx) / &
                  real((num_atoms_use / pub_bsunfld_num_atoms_prim), kind=DP)
          end do
       end if ! if ikpt_red <= nkpts_red

       call comms_barrier
       ! reduce values across procs and store them on the root proc
       call comms_reduce('SUM', spect_func_red, root=pub_root_proc_id)
       call comms_reduce('SUM', eigvals_func_red, root=pub_root_proc_id)
       call comms_barrier

       ! reset spectral function / eigenvalues on procs other than the root
       if (.not.pub_on_root) then
          spect_func_red = cmplx(0.0_DP, 0.0_DP, kind=DP)
          eigvals_func_red = 0.0_DP
       end if

       if (pub_on_root) then
          do idx = 1, 1 + pub_bsunfld_nprojatoms
             specfuncred_unit = utils_unit()
             open(specfuncred_unit,file=trim(adjustl(specfuncred_file(idx))), &
                  form='FORMATTED',action='WRITE',iostat=ierr)
             if (ierr.ne.0) then
                call utils_abort('Error in opening '//&
                     &trim(adjustl(specfuncred_file(idx))))
             end if

             do ikpt = 1,nkpts_red

                do ival = 1,neigvals
                   write(specfuncred_unit,'(3F20.10,I8,3F12.6,I8)') &
                        HARTREE_IN_EVS * eigvals_func_red(ikpt,ival), &
                        real(spect_func_red(ikpt,ival,idx),kind=DP), &
                        aimag(spect_func_red(ikpt,ival,idx)), ival, &
                        kpt(revuni_entry(ikpt,1))%X * ANGSTROM, &
                        kpt(revuni_entry(ikpt,1))%Y * ANGSTROM, &
                        kpt(revuni_entry(ikpt,1))%Z * ANGSTROM, ikpt
                end do
             end do

             close(specfuncred_unit)
          end do
       end if ! if on root

    end do ! over reduced kpts

    call comms_barrier
    if (pub_on_root .and. pub_print_qc) then
       call utils_qc_print('bsunfold_spect_func_red_real:',real(spect_func_red(1,1,1),kind=DP))
       call utils_qc_print('bsunfold_spect_func_red_imag:',aimag(spect_func_red(1,1,1)))
    end if

    if (pub_on_root) then
       do idx = 1, 1 + pub_bsunfld_nprojatoms
          specfunc_unit = utils_unit()
          open(specfunc_unit,file=trim(adjustl(specfunc_file(idx))), &
               form='FORMATTED',action='WRITE',iostat=ierr)
          if (ierr.ne.0) then
             call utils_abort('Error in opening '//&
                  &trim(adjustl(specfunc_file(idx))))
          end if

          do ikpt = 1, nkpts
             ikpt_red = uni_entry(ikpt)
             do ival = 1,neigvals
                write(specfunc_unit,'(3F20.10,I8,3F20.10,I8)') &
                     HARTREE_IN_EVS * eigvals_func_red(ikpt_red,ival), &
                     real(spect_func_red(ikpt_red,ival,idx),kind=DP), &
                     aimag(spect_func_red(ikpt_red,ival,idx)), ival, &
                     kpt(ikpt)%X * ANGSTROM, kpt(ikpt)%Y * ANGSTROM, &
                     kpt(ikpt)%Z * ANGSTROM, ikpt
             end do
          end do

          close(specfunc_unit)
       end do
    end if ! if on root

    ! Destroy workspace
    ! Allocate workspace for diagonalisation and unfolding
    do ikpt = 1, SCnkpts_red
       do idx = 1, spin_dim
          call dense_destroy(dense_eigst(ikpt,idx))
       end do
    end do
    deallocate(dense_eigst, stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','dense_eigst', ierr)
    deallocate(buffer_col, stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','buffer_col',ierr)
    deallocate(list_kptred_done, stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','list_kptred_done',ierr)
    deallocate(overlap_array,stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','overlap_array',ierr)
    deallocate(eigst_kpts,stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','eigst_kpts',ierr)
    deallocate(eigvl_kpts,stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','eigvl_kpts',ierr)
    deallocate(uni_entry,stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','uni_entry',ierr)
    deallocate(SCuni_entry,stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','SCuni_entry',ierr)
    deallocate(revuni_entry, stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','revuni_entry',ierr)
    deallocate(kpt,stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','kpt',ierr)
    deallocate(SCkpt,stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','SCkpt',ierr)
    deallocate(spect_func_red, stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','spect_func_red',ierr)
    deallocate(eigvals_func_red, stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','eigvals_func_red',ierr)
    deallocate(ngwfs_idx_use_cells, stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','ngwfs_idx_use_cells',ierr)
    deallocate(ngwfs_idx_use, stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','ngwfs_idx_use',ierr)
    deallocate(atoms_idx_use, stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','atoms_idx_use',ierr)
    deallocate(pos_cells_use, stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','pos_cells_use',ierr)
    deallocate(row_spectral, stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','row_spectral',ierr)
    deallocate(ngwfidx_wrt_cell, stat=ierr)
    call utils_dealloc_check('bsunfold_calculate','ngwfidx_wrt_cell',ierr)
    deallocate(to_include, stat=ierr)
    call utils_dealloc_check('bsunfold_calculate', 'to_include', ierr)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving bsunfold_calculate'

  end subroutine bsunfold_calculate

end module bandstructure
