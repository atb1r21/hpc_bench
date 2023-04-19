! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=================================================================!
!                                                                 !
!       Electron energy loss spectra (EELS) Module                !
!                                                                 !
! This module calculates the electron energy loss spectra...      !
!-----------------------------------------------------------------!
! Initially written by Laura Ratcliff in October 2011.            !
! Further work by Nicholas Hine in July 2012.                     !
! Further amended by Laura Ratcliff in July 2013                  !
! and by Edward Tait between October 2014 and January 2015.       !
!=================================================================!

module eels

  use constants, only: DP

  implicit none

  private

  public :: eels_calculate

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eels_calculate(cur_spin,eigs_dens,eigen_en, &
       rep,ngwf_basis,proj_basis,core_basis,core_wvfns,mdl)

    !========================================================================!
    ! This subroutine...                                                     !
    !------------------------------------------------------------------------!
    ! Written by Laura Ratcliff in October 2011.                             !
    !========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout
    use dense, only: DEM, dense_create, dense_destroy
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use projectors, only: PROJECTOR_SET, &
         projectors_func_pos_ovlp_box, projectors_func_pos_ovlp_box_gradg
    use rundat, only: pub_paw, pub_debug_on_root, pub_eels_realspace
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_axpy, &
         sparse_transpose, sparse_copy
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_banner

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: core_basis
    type(NGWF_REP), intent(in) :: rep
    type(DEM), intent(in) :: eigs_dens
    real(kind=DP), intent(in) :: eigen_en(:) ! hamiltonian eigenvalues
    integer, intent(in) :: cur_spin
    type(PROJECTOR_SET), intent(inout) :: core_wvfns ! not needed to be inout but same as eigenstates for now
    type(MODEL), intent(in) :: mdl

    ! Local variables
    type(DEM), allocatable :: eels_mat_elements(:) ! EELS matrix elements of eigenvectors
    type(SPAM3) :: core_ngwf_pos_gradg(3), ngwf_core_pos_gradg(3)
    type(SPAM3) :: core_ngwf_pos_rspace(3), ngwf_core_pos_rspace(3)
    type(SPAM3) :: paw_correction(3), total(3)
    integer :: ierr                    ! memory allocation error flag
    integer :: xyz                     ! cartesian direction
    ! agrecocmplx
    logical :: loc_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &eels_calculate'

    ! ndmh: test for suitable settings before doing anything
    if (.not.pub_paw) then
       if (pub_on_root) then
          write(stdout, '(a)') 'WARNING in eels_calculate: eels_calculate has &
               &been set, but EELS calculations'
          write(stdout, '(a)') 'are only possible in combination with PAW core &
               &reconstruction. However, PAW is'
          write(stdout, '(a)') 'not active in the current calculation.'
       end if
       return
    end if
    if (core_basis%num==0) then
       if (pub_on_root) then
          write(stdout, '(a)') 'WARNING in eels_calculate: eels_calculate has &
               &been set, but none of the'
          write(stdout, '(a)') 'species in the calculation have been given core &
               &wavefunctions:'
          write(stdout, '(a)') 'Please add a %block species_core_wf block in &
               &your input file'
       end if
       return
    end if

    if (pub_on_root) write(stdout,'(/a)') utils_banner('=', &
         'Calculation of matrix elements between core orbitals and eigenstates')

    ! agrecocmplx
    loc_cmplx = rep%ngwfs_on_grid(1)%iscmplx

    ! lr408: Allocate spectra arrays
    allocate(eels_mat_elements(3),stat=ierr)
    call utils_alloc_check('eels_calculate','eels_mat_elements',ierr)

    do xyz=1,3
       call dense_create(eels_mat_elements(xyz),ngwf_basis%num,core_basis%num,&
            iscmplx=.true.)
    end do

    ! Calculate position matrix elements between eigenfunctions and core wvfns
    ! 1) Calculate between NGWFs and core wvfns
    ! ewt23: Grad_G method requires complex matrices.
    do xyz=1,3
       core_ngwf_pos_gradg(xyz)%structure = 'C'//rep%postfix
       call sparse_create(core_ngwf_pos_gradg(xyz),iscmplx=.true.)
       ! agrecocmplx: complex if using complex NGWFs
       core_ngwf_pos_rspace(xyz)%structure = 'C'//rep%postfix
       call sparse_create(core_ngwf_pos_rspace(xyz),iscmplx=loc_cmplx)
       paw_correction(xyz)%structure = 'C'//rep%postfix
       ! agrecocmplx
       call sparse_create(paw_correction(xyz),iscmplx=loc_cmplx)
       total(xyz)%structure = 'C'//rep%postfix
       call sparse_create(total(xyz), iscmplx=.true.)
    end do

    do xyz=1,3
       ngwf_core_pos_gradg(xyz)%structure = 'B'//rep%postfix
       call sparse_create(ngwf_core_pos_gradg(xyz),iscmplx=.true.)
       ! agrecocmplx: complex if using complex NGWFs
       ngwf_core_pos_rspace(xyz)%structure = 'B'//rep%postfix
       call sparse_create(ngwf_core_pos_rspace(xyz),iscmplx=loc_cmplx)
    end do

    ! rc2013: EMBED_FIX!
    call projectors_func_pos_ovlp_box_gradg(ngwf_core_pos_gradg, &
         rep%ngwfs_on_grid(1),ngwf_basis,core_basis,core_wvfns,mdl%fftbox,mdl%cell)

    ! rc2013: EMBED_FIX!
    call projectors_func_pos_ovlp_box(ngwf_core_pos_rspace, &
         rep%ngwfs_on_grid(1),ngwf_basis,core_basis,core_wvfns, &
         mdl%fftbox,mdl%cell)

    do xyz=1,3
       call sparse_transpose(core_ngwf_pos_gradg(xyz), &
            ngwf_core_pos_gradg(xyz))
       call sparse_transpose(core_ngwf_pos_rspace(xyz), &
            ngwf_core_pos_rspace(xyz))
    end do

    do xyz=3,1,-1
       call sparse_destroy(ngwf_core_pos_gradg(xyz))
       call sparse_destroy(ngwf_core_pos_rspace(xyz))
    end do


    ! 2) Calculate PAW correction term
    ! rc2013: EMBED_FIX!
    call eels_core_paw_corr(paw_correction, core_basis, proj_basis, &
         rep%sp_overlap%p,mdl)

    ! 3) Add the two together
    do xyz=1,3
       call sparse_copy(total(xyz),core_ngwf_pos_gradg(xyz))
       call sparse_axpy(total(xyz),paw_correction(xyz), &
            cmplx(1.0_DP,0.0_DP,kind=DP))
    end do

    ! now want to get matrix elements for the eigenvectors
    do xyz=1,3
       call eels_mat_els(eels_mat_elements(xyz),ngwf_basis%num, &
            core_basis%num,eigs_dens,total(xyz))
    end do
    call eels_print_mat_els(eels_mat_elements,ngwf_basis%num, &
         core_basis%num,eigen_en,cur_spin,rep,mdl%cell,'_grad_total')

    call eels_write_elnes_bin(mdl,rep,core_basis, &
         ngwf_basis%num,eels_mat_elements,"grad")


    ! Here we output the .elnes_bin files and
    ! .txt matrix element listings as computed
    ! using the real-space cartesian method.
    ! This is here to permit comparison with the
    ! Grad_G method - which is thought to work
    ! better.

    if (pub_eels_realspace) then
       do xyz=1,3
          call sparse_copy(total(xyz),core_ngwf_pos_rspace(xyz))
          call sparse_axpy(total(xyz),paw_correction(xyz),&
               cmplx(1.0_DP,0.0_DP,kind=DP))
       end do

       ! now we want to get matrix elements for the eigenvectors
       do xyz=1,3
          call eels_mat_els(eels_mat_elements(xyz),ngwf_basis%num, &
               core_basis%num,eigs_dens,total(xyz))
       end do
       call eels_print_mat_els(eels_mat_elements,ngwf_basis%num, &
            core_basis%num,eigen_en,cur_spin,rep,mdl%cell,'_realspace_total')

       call eels_write_elnes_bin(mdl,rep,core_basis,&
            ngwf_basis%num,eels_mat_elements,"realspace")

    end if

    do xyz=3,1,-1
       call sparse_destroy(total(xyz))
    end do

    do xyz=1,3
       call eels_mat_els(eels_mat_elements(xyz),ngwf_basis%num, &
            core_basis%num,eigs_dens,core_ngwf_pos_gradg(xyz))
    end do
    call eels_print_mat_els(eels_mat_elements,ngwf_basis%num, &
         core_basis%num,eigen_en,cur_spin,rep,mdl%cell,'_gradg_cart')

    if (pub_eels_realspace) then
       do xyz=1,3
          call eels_mat_els(eels_mat_elements(xyz),ngwf_basis%num, &
               core_basis%num,eigs_dens,core_ngwf_pos_rspace(xyz))
       end do
       call eels_print_mat_els(eels_mat_elements,ngwf_basis%num, &
            core_basis%num,eigen_en,cur_spin,rep,mdl%cell,'_rspace_cart')
    end if

    do xyz=3,1,-1
       call sparse_destroy(core_ngwf_pos_gradg(xyz))
       call sparse_destroy(core_ngwf_pos_rspace(xyz))
    end do

    ! Here we write out just the PAW corrections to the
    ! matrix elements.
    do xyz=1,3
       call eels_mat_els(eels_mat_elements(xyz),ngwf_basis%num, &
            core_basis%num,eigs_dens,paw_correction(xyz))
    end do
    call eels_print_mat_els(eels_mat_elements,ngwf_basis%num, &
         core_basis%num,eigen_en,cur_spin,rep,mdl%cell,'_paw')
    do xyz=3,1,-1
       call sparse_destroy(paw_correction(xyz))
    end do

    ! lr408: Deallocate spectra arrays
    do xyz=3,1,-1
        call dense_destroy(eels_mat_elements(xyz))
    end do
    deallocate(eels_mat_elements,stat=ierr)
    call utils_dealloc_check('eels_calculate', &
         'eels_mat_elements',ierr)

    if (pub_on_root) write(stdout,'(a)') repeat('=',80)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &eels_calculate'


  end subroutine eels_calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eels_mat_els(eels_mat_elements,ngwf_num, &
       core_num,eigs_dens,core_ngwf_pos_elements)

    !==================================================================!
    ! This subroutine ...                                              !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff in October 2011.                       !
    ! Loosely based on properties_opt_mat_els                          !
    ! Modified to use moplex matricies by Edward Tait - Jan 2015       !
    !==================================================================!

    use comms, only: comms_barrier
    use constants, only: DP, stdout
    use dense, only: DEM, dense_axpy, dense_scale, &
         dense_create, dense_destroy, dense_convert, dense_product
    use rundat, only: pub_debug_on_root
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_axpy

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: eels_mat_elements
    type(DEM), intent(in) :: eigs_dens
    type(SPAM3), intent(in) :: core_ngwf_pos_elements
    integer, intent(in) :: ngwf_num, core_num

    ! Local variables
    type(DEM) :: core_ngwf_pos_dens_cmplx
    type(DEM) :: eigs_dens_cmplx
    type(DEM) :: eigs_pos_prod
    ! agrecocmplx
    logical :: loc_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering eels_mat_els'

    ! lr408: Assume everything real for now
    ! agrecocmplx
    loc_cmplx = eigs_dens%iscmplx

    ! ndmh: Set opt_mat_elements to zero
    call dense_scale(eels_mat_elements,cmplx(0.0_DP,0.0_DP,kind=DP))

    call dense_create(core_ngwf_pos_dens_cmplx,core_num,ngwf_num,iscmplx=.true.)
    call dense_convert(core_ngwf_pos_dens_cmplx,core_ngwf_pos_elements)

    ! ndmh: Create temporary matrices
    ! agrecocmplx: only needed if original eigs_dens is real
    if (.not.loc_cmplx) then
       call dense_create(eigs_dens_cmplx,eigs_dens%nrows, &
            eigs_dens%mcols,iscmplx=.true.)
       call dense_scale(eigs_dens_cmplx,cmplx(0.0_DP,0.0_DP,kind=DP))
       call dense_axpy(eigs_dens_cmplx,eigs_dens,cmplx(1.0_DP,0.0_DP,kind=DP))
    end if
    call dense_create(eigs_pos_prod,ngwf_num,core_num,iscmplx=.true.)
    ! agrecocmplx
    ! real case: use complex version of eigs_dens
    if (.not.loc_cmplx) then
       call dense_product(eigs_pos_prod,eigs_dens_cmplx,&
            core_ngwf_pos_dens_cmplx, opA='T',opB='T')
    ! complex case: use original complex eigs_dens
    else
       call dense_product(eigs_pos_prod,eigs_dens,&
            core_ngwf_pos_dens_cmplx, opA='T',opB='T')
    end if
    call dense_axpy(eels_mat_elements,eigs_pos_prod,&
         cmplx(1.0_DP,0.0_DP,kind=DP))

    ! eae32: Synchronise procs before writing optical matrix elements to file
    call comms_barrier

    ! ndmh: Clean up temporary matrices
    call dense_destroy(eigs_pos_prod)
    call dense_destroy(core_ngwf_pos_dens_cmplx)
    ! agrecocmplx
    if (.not.loc_cmplx) then
       call dense_destroy(eigs_dens_cmplx)
    end if

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving eels_mat_els'

  end subroutine eels_mat_els


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !========================================================================!
  ! This subroutine writes out an OptaDoS compatible .elnes_bin file       !
  ! The specification for this file may be found in appendix A.4           !
  ! of the OptaDoS user guide                                              !
  !------------------------------------------------------------------------!
  ! Written by Ed Tait in October 2014.                                    !
  !========================================================================!

  subroutine eels_write_elnes_bin(mdl,rep,core_basis,ngwf_num,&
       mat_elements_in,suffix)

    use comms, only: pub_on_root, comms_barrier
    use constants, only: DP, stdout
    use dense, only: DEM, dense_get_element
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use parallel_strategy, only: PARAL_INFO
    use projectors, only: PROJECTOR_SET
    use rundat, only: pub_rootname, pub_debug_on_root, pub_num_spins
    use utils, only:  utils_unit, utils_open_unit_check, &
        utils_close_unit_check, utils_alloc_check, utils_dealloc_check, &
        utils_feature_not_supported

    implicit none

    ! Arguments
    type(MODEL), intent(in) :: mdl
    type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: core_basis
    integer, intent(in) :: ngwf_num
    type(DEM), intent(in) :: mat_elements_in(:)
    character(len=*), intent(in) :: suffix
    ! Local Variables
    real(kind=DP),parameter :: elnes_bin_file_version = 1.0_DP
    character(len=80),parameter :: elnes_bin_header = "Output from ONETEP"
    integer:: max_eigenvalues
    integer:: core_num
    integer,allocatable:: num_eigenvalues(:) !1:num_spins
    integer,allocatable:: species(:) !1:core_num
    integer,allocatable:: species_index(:) !1:core_num
    integer,allocatable:: prin_quantum_number(:) !1:core_num
    ! s=1 px=2 py=3 pz=4 etc
    integer,allocatable:: ang_mom_ident_number(:) !1:core_num
    complex(kind=DP),allocatable :: elnes_mat_els(:,:,:,:,:)
    ! Output and allocation util varibales
    integer :: output_unit, io_status,ierr
    character(len=255) :: output_file
    character(len=6) :: file_type
    integer :: current_core
    integer :: iatom
    integer :: pspc
    integer :: xyz
    integer :: final
    integer :: l,m
    integer :: elm_index,sp_ind
    integer,allocatable  :: species_ion_count(:)
    ! Number of wavefunctions we've seen on an atom so far
    integer,allocatable  :: ion_pcount(:)
    integer,allocatable  :: ion_species_index(:)
    complex(kind=DP) :: mat_el
    ! CASTEP element order array
    integer,allocatable  :: castep_order(:)

    ! How many core wavefunctions do we have?
    core_num = mdl%regions(1)%par%num_corewfs
    max_eigenvalues = ngwf_num

    if (rep%postfix=='') then
       file_type='_val'
    else if (rep%postfix=='j') then
       file_type='_joint'
    end if

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering eels_write_elnes_bin'

    if (pub_on_root) then
       write(output_file,*) &
            trim(pub_rootname)//trim(file_type)//'_'//suffix//'.elnes_bin'
       output_file = adjustl(output_file)
       output_unit = utils_unit()

       write(stdout,'(3a)',advance='no') 'Writing "',trim(output_file), '" ...'

       ! OptaDos Expects a big endian binary file
       ! try to enforce this - note the convert= keyword
       ! is non-standard...
#ifndef F2008
       ! jd: 'convert' is non-standard Fortran.
       open(unit=output_unit, form="unformatted" ,file=trim(output_file), &
            action="write",iostat=io_status &
#ifndef BIGENDIAN
            ,convert="BIG_ENDIAN")
#else
            )
#endif

#else
       ! jd: Here we want to be F2008-compliant
#ifndef BIGENDIAN
       call utils_feature_not_supported('eels_write_elnes_bin:'//&
            'OptaDos big-endian output has not been compiled in because it &
            &needs "convert", which is not part of the F2008 standard, &
            &and you asked for an F2008-compliant binary with -D2008.')
#endif

#endif
       call utils_open_unit_check('eels_write_elnes_bin','output_unit', &
            io_status)

       allocate (num_eigenvalues(pub_num_spins),stat=ierr)
       call utils_alloc_check('eels_write_elnes_bin','num_eigenvalues',ierr)

       allocate (species(core_num),stat=ierr)
       call utils_alloc_check('eels_write_elnes_bin','species',ierr)

       allocate (species_index(core_num),stat=ierr)
       call utils_alloc_check('eels_write_elnes_bin','species_index',ierr)

       allocate (prin_quantum_number(core_num),stat=ierr)
       call utils_alloc_check('eels_write_elnes_bin','prin_quantum_number', &
            ierr)

       allocate (ang_mom_ident_number(core_num),stat=ierr)
       call utils_alloc_check('eels_write_elnes_bin','ang_mom_ident_number',&
            ierr)

       allocate(elnes_mat_els(core_num,max_eigenvalues,3,1,1),stat=ierr)
       call utils_alloc_check('eels_write_elnes_bin','elnes_mat_els',&
            ierr)

       allocate(species_ion_count(mdl%regions(1)%par%num_pspecies),stat=ierr)
       call utils_alloc_check('eels_write_elnes_bin','species_ion_count',&
            ierr)

       allocate(ion_pcount(mdl%regions(1)%par%nat),stat=ierr)
       call utils_alloc_check('eels_write_elnes_bin','ion_pcount',&
            ierr)

       allocate(ion_species_index(mdl%regions(1)%par%nat),stat=ierr)
       call utils_alloc_check('eels_write_elnes_bin','ion_species_index',&
            ierr)

       allocate(castep_order(mdl%regions(1)%par%num_pspecies),stat=ierr)
       call utils_alloc_check('eels_write_elnes_bin','castep_order',&
            ierr)

       species_ion_count = 0
       ion_species_index = 0
       ion_pcount = 0

       ! Initialise the lookup table for ion species number,
       ! This way ion species number should match the order
       ! ions appear in the .dat file coords block.
       do iatom = 1,mdl%regions(1)%par%nat
          pspc = mdl%elements(iatom)%pspecies_number
          species_ion_count(pspc) = species_ion_count(pspc) + 1
          ion_species_index(iatom) = species_ion_count(pspc)
       end do


       ! A little loop here to create a lookup table for
       ! going between onetep pspecies_number and castep internal
       ! ordering. This just puts species with a label to the end
       ! otherwise preserving the order within the position block.
       castep_order = -1 ! setup
       sp_ind = 1
       do elm_index = 1, mdl%regions(1)%par%nat
          ! This can be done better, handle single char species V
          if (ichar(mdl%elements(elm_index)%species_id(3:3))&
               .eq. ichar(' ')) then
             if (castep_order(mdl%elements(elm_index)&
                  %pspecies_number) .lt. 0) then
                castep_order(mdl%elements(elm_index)&
                     %pspecies_number) = sp_ind
                sp_ind = sp_ind + 1
             end if
          end if
       end do
       do elm_index = 1, mdl%regions(1)%par%nat
          ! This can be done better, handle single char species V
          if (.not. ichar(mdl%elements(elm_index)%species_id(3:3))&
               .eq. ichar(' ')) then
             if (castep_order(mdl%elements(elm_index)%pspecies_number).lt.0) then
                castep_order(mdl%elements(elm_index)%pspecies_number) = sp_ind
                sp_ind = sp_ind + 1
             end if
          end if
       end do

    end if


    do current_core = 1,core_basis%num
       ! Get index, angular momentum, species

       ! each atom has its own projectors
       ! can look up what atom is attached to
       ! a given wave function with the FUNC_BASIS
       ! type we're passed.

       ! Get atom
       ! Get species of atom
       ! get index in species
       ! get matrix emement
       ! put each in the appropriate bit
       !  of the appropriate array
       if (pub_on_root) then

          ! Find out about this core wavefunction
          iatom = mdl%regions(1)%par%orig_atom(&
               core_basis%atom_of_func(current_core))
          pspc = mdl%elements(iatom)%pspecies_number
          ion_pcount(iatom) = ion_pcount(iatom) + 1
          species(current_core) = castep_order(pspc)
          ! Output debugging information
          if (pub_debug_on_root) then
             write(stdout,'(a,i8)') "Core --------> ", &
                  current_core
             write(stdout,'(a,a)') "Symbol         ", &
                  mdl%elements(iatom)%symbol
             write(stdout,'(a,i8)') "iatom          ", &
                  iatom
             write(stdout,'(a,i8)') "pspecies       ", &
                  pspc
             write(stdout,'(a,i8)') "Total corewfs  ", &
                  mdl%paw_sp(pspc)%n_core_wfs_tot
             write(stdout,'(a,i8)') "ion_speies_ind ", &
                  ion_species_index(iatom)
             write(stdout,'(a,i8)') "spc_ion_count  ", &
                  species_ion_count(pspc)
             write(stdout,'(a,i8)') "ion_pcount     ", &
                  ion_pcount(iatom)
             write(stdout,'(a,a)') "PAW CWF Name   ", &
                  mdl%paw_sp(pspc)%core_wf_name
             write(stdout,'(a,a)') "Elms CWF Name  ", &
                  mdl%species(mdl%elements(iatom)%&
                     species_number)%core_wf_name
          end if

          !what is our index iwthin a species?
          species_index(current_core) = ion_species_index(iatom)

          prin_quantum_number(current_core) = &
               mdl%paw_sp(pspc)%n_core_wf(&
               mdl%paw_sp(pspc)%icore_wf_tot(ion_pcount(iatom)))

          l = mdl%paw_sp(pspc)%l_core_wf_tot(ion_pcount(iatom))
          m = mdl%paw_sp(pspc)%m_core_wf_tot(ion_pcount(iatom))
          ang_mom_ident_number(current_core) = l*l+l+1+m

       end if

       ! Get matrix element for this core wvfn and eigenstate and put it in
       ! elnes_mat_els
       do final = 1,max_eigenvalues
          do xyz=1,3
             call dense_get_element(mat_el,mat_elements_in(xyz), &
                  final,current_core)
             if (pub_on_root) &
                  elnes_mat_els(current_core,final,xyz,1,1) = mat_el
          end do
       end do

    end do

    ! Now print elnes_mat_els to file
    if (pub_on_root) then

       ! Basic info checked / printed by OptaDOS
       write(output_unit) elnes_bin_file_version
       write(output_unit) elnes_bin_header
       write(output_unit) core_num
       write(output_unit) max_eigenvalues
       write(output_unit) 1 !number of kpoints
       write(output_unit) 1 !number of spins
       write(output_unit) species(1:core_num)
       write(output_unit) species_index(1:core_num)
       write(output_unit) prin_quantum_number(1:core_num)
       write(output_unit) ang_mom_ident_number(1:core_num)

       ! write out matrix elements
       write(output_unit) (((elnes_mat_els(current_core,final,xyz,1,1), &
            current_core=1,core_num),final=1,max_eigenvalues),xyz=1,3)

       close(unit=output_unit,iostat=io_status)
       call utils_close_unit_check('eels_write_elnes_bin','output_unit', &
            io_status)

       deallocate (num_eigenvalues,stat=ierr)
       call utils_dealloc_check('eels_write_elnes_bin', &
            'num_eigenvalues',ierr)

       deallocate (species,stat=ierr)
       call utils_dealloc_check('eels_write_elnes_bin','species',ierr)

       deallocate (species_index,stat=ierr)
       call utils_dealloc_check('eels_write_elnes_bin','species_index',ierr)

       deallocate (prin_quantum_number,stat=ierr)
       call utils_dealloc_check('eels_write_elnes_bin','prin_quantum_number', &
            ierr)

       deallocate (ang_mom_ident_number,stat=ierr)
       call utils_dealloc_check('eels_write_elnes_bin','ang_mom_ident_number', &
            ierr)

       deallocate (elnes_mat_els,stat=ierr)
       call utils_dealloc_check('eels_write_elnes_bin','elnes_mat_els',&
            ierr)

       deallocate(species_ion_count,stat=ierr)
       call utils_dealloc_check('eels_write_elnes_bin','species_ion_count', &
            ierr)

       deallocate(ion_pcount,stat=ierr)
       call utils_dealloc_check('eels_write_elnes_bin','ion_pcount',ierr)

       deallocate(ion_species_index ,stat=ierr)
       call utils_dealloc_check('eels_write_elnes_bin','ion_species_index', &
            ierr)

       deallocate(castep_order ,stat=ierr)
       call utils_dealloc_check('eels_write_elnes_bin','castep_order',ierr)


    end if

    call comms_barrier

    if (pub_on_root) write(stdout,'(a)') 'done'

  end subroutine eels_write_elnes_bin


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eels_print_mat_els(eels_mat_elements,ngwf_num, &
         core_num,eigen_en,cur_spin,rep,cell,tag)

    use comms, only: pub_on_root, comms_barrier
    use constants, only: DP, stdout, HARTREE_IN_EVS, UP, DN
    use dense, only: DEM, dense_get_element
    use function_basis, only: FUNC_BASIS
    use geometry, only: operator(.CROSS.), operator(.DOT.)
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_rootname, pub_spectra_print_mat_els, &
         pub_debug_on_root, pub_num_spins, pub_num_kpoints, PUB_1K
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_create, sparse_destroy, &
         sparse_axpy, sparse_copy
    use utils, only:  utils_unit, utils_open_unit_check, utils_close_unit_check,&
         utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: eels_mat_elements(:)
    real(kind=DP), intent(in) :: eigen_en(:) ! hamiltonian eigenvalues
    type(NGWF_REP), intent(in) :: rep
    type(CELL_INFO), intent(in) :: cell
    integer, intent(in) :: ngwf_num, core_num, cur_spin
    character(len=*), intent(in) :: tag

    ! Local variables
    real(kind=DP) :: cell_volume
    complex(kind=DP) :: eels_mat_el
    integer :: xyz, i, f, count, max, homo
    integer :: output_unit, io_status
    character(len=256) :: output_file  ! output file name buffer
    character(len=6) :: file_type

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering eels_print_mat_els'

    ! jme: KPOINTS_DANGER
    ! One instruction of this subroutine enforces a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine eels_calculate not ready yet for more&
         & than one k-point.')

    ! lr408: At some point, actually make use of this to scale matrix elements
    ! lr408: within ONETEP
    cell_volume = abs((cell%a1 .CROSS. cell%a2) .DOT. cell%a3)
    call utils_assert(cell_volume /= -1.0_DP,'Cell volume is zero')

    ! Determine prefix from NGWF_REP postfix - change this to get file type directly
    if (rep%postfix=='') then
       file_type='_val'
    else if (rep%postfix=='j') then
       file_type='_joint'
    end if

    homo = rep%n_occ(cur_spin,PUB_1K)
    if (homo==0) homo=1

    ! eae32: calculate number of states on all procs
    if (2*homo > ngwf_num) then
       max = ngwf_num
    else
       max = 2*homo
    end if

    max=ngwf_num

    ! only print out matrix elements if required
    if (pub_on_root .and. pub_spectra_print_mat_els) then

       write(stdout,'(a)') ''

       write(stdout,'(a)')&
            '================ Writing EELS matrix elements &
            &================'
       if (pub_num_spins == 1) then
          write(output_file,*)trim(pub_rootname)//trim(file_type)//'_EELS_MAT_ELS'//trim(tag)//'.txt'
       else if (cur_spin == UP) then
          write(output_file,*)trim(pub_rootname)//trim(file_type)//'_EELS_MAT_ELS_UP'//trim(tag)//'.txt'
       else if (cur_spin == DN) then
          write(output_file,*)trim(pub_rootname)//trim(file_type)//'_EELS_MAT_ELS_DN'//trim(tag)//'.txt'
       end if

       output_file = adjustl(output_file)

       ! cks: print output warning
       write(stdout,'(3a)',advance ='no') &
            'Writing "', trim(output_file),'" ...'

       ! cks: get a unit number that is free
       output_unit = utils_unit()

       open(unit=output_unit, form="formatted" ,file=trim(output_file), &
            action="write",iostat=io_status)
       call utils_open_unit_check('eels_print_mat_els','output_file', &
            io_status)

       ! cks: write first line
       write(output_unit,'(a)')'# initial | Energy (eV) | final  | Energy (eV) |&
           & Matrix el.  | Trans. E (eV) '
       !write(output_unit,'(a,2x,I2)') '# Spin component',cur_spin

       write(output_unit,'(I5,2x,I5)') core_num,max,1
       write(output_unit,'(F24.12)') cell_volume

    end if

    ! to do core eigenvalues properly, need to loop over species isp=1,nsp
    ! then use paw_sp(isp)%core_wf_eig(i)

    ! eae32: dense_get_element calls across all procs
    do xyz=1,3

       if (pub_on_root .and. pub_spectra_print_mat_els) then
          write(output_unit,'(a,I1)') '# Cartesian component ',xyz
       end if

       count = 1

       do i=1,core_num ! initial state
          do f=1,max ! final state

             call dense_get_element(eels_mat_el,eels_mat_elements(xyz),f,i)
             if (pub_on_root .and. pub_spectra_print_mat_els) then
                write(output_unit,'(2x,2(I5,4x,F12.6,2x),2(F24.12,2x))') &
                     i,0.0_dp,&!core_wf_eig(i)*HARTREE_IN_EVS,&
                     f,eigen_en(f)*HARTREE_IN_EVS,&
                     abs(eels_mat_el),&
                     (eigen_en(f))*HARTREE_IN_EVS
             end if
          end do
       end do

    end do ! loop over xyz

    if (pub_on_root .and. pub_spectra_print_mat_els) then

       close(unit=output_unit,iostat=io_status)
       call utils_close_unit_check('properties_spectra','output_unit', &
            io_status)

       ! cks: notify of end of output
       write(stdout,*)' done'

       write(stdout,'(a)') repeat('=',80)

    end if

    call comms_barrier

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving eels_print_mat_els'

  end subroutine eels_print_mat_els

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eels_core_paw_corr(core_ngwf_pos,core_basis,proj_basis, &
         sp_overlap,mdl)

    use comms, only: pub_my_proc_id
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use geometry, only: operator(*), operator(+)
    use model_type, only: MODEL
    use parallel_strategy, only: PARAL_INFO
    use paw, only: paw_core_pw_position_operator
    use rundat, only: pub_debug_on_root
    use sparse, only: SPAM3, sparse_create, sparse_get_element, &
         sparse_destroy, sparse_put_element, sparse_copy, sparse_get_par, &
         sparse_transpose_structure, sparse_transpose, sparse_product
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: sp_overlap
    type(SPAM3), intent(inout):: core_ngwf_pos(3)
    type(FUNC_BASIS), intent(in) :: core_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(MODEL), intent(in) :: mdl

    ! Local Variables
    real(kind=DP) :: R_atom(3),r_el,o_el
    integer       :: xyz
    integer       :: loc_iproj,iproj
    integer       :: jcwf
    integer       :: loc_iat,iat
    type(SPAM3)   :: ps_overlap
    type(SPAM3)   :: core_pw_pos(3),core_pw_ovlp
    ! agrecocmplx
    logical :: loc_cmplx
    type(SPAM3) :: temp_cmplx_mat
    type(PARAL_INFO), pointer :: par

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering eels_core_paw_corr'

    ! Start timer
    call timer_clock('eels_core_paw_corr',1)

    ! agrecocmplx
    loc_cmplx = sp_overlap%iscmplx

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, sp_overlap)

    do xyz=1,3
       ! agrecocmplx: real in any case
       core_pw_pos(xyz)%structure = 'P'
       call sparse_create(core_pw_pos(xyz))
    end do
    ! agrecocmplx: need temp complex copy in complex
    ! case in order to perform sparse_product later
    if (loc_cmplx) then
       call sparse_create(temp_cmplx_mat,core_pw_pos(1),iscmplx=loc_cmplx)
    end if

    ! agreco: this is not actually used?
    core_pw_ovlp%structure = 'P'
    call sparse_create(core_pw_ovlp)

    call sparse_transpose_structure(ps_overlap%structure,sp_overlap)
    ! agrecocmplx
    call sparse_create(ps_overlap,iscmplx=loc_cmplx)
    call sparse_transpose(ps_overlap,sp_overlap)

    ! ndmh: calculate sphere part of position operator between core wvfns
    ! ndmh: and partial waves
    call paw_core_pw_position_operator(core_pw_pos,mdl%paw_sp)

    ! ndmh: cycle over projectors on this proc, applying correction to
    ! ndmh: r_sphere to move it to the atom centre
    do loc_iproj=1,proj_basis%num_on_proc(pub_my_proc_id)
       iproj = loc_iproj + proj_basis%first_on_proc(pub_my_proc_id) - 1
       iat = proj_basis%atom_of_func(iproj)
       loc_iat = iat - par%first_atom_on_proc(pub_my_proc_id) + 1
       do jcwf=core_basis%first_on_atom(iat), &
            core_basis%first_on_atom(iat)+core_basis%num_on_atom(iat)-1

          ! Extract overlap element
          call sparse_get_element(o_el,core_pw_ovlp,jcwf,iproj)

          R_atom(1) = par%elements_on_proc(loc_iat)%centre%x
          R_atom(2) = par%elements_on_proc(loc_iat)%centre%y
          R_atom(3) = par%elements_on_proc(loc_iat)%centre%z

          ! Extract element from r_sphere and shift by R_atom*o_el
          ! ddor: get elements for one direction only if axis is specified
          do xyz=1,3
             call sparse_get_element(r_el,core_pw_pos(xyz),jcwf,iproj)
             !r_el = R_atom(xyz) * o_el + r_el
             !r_el = o_el
             call sparse_put_element(r_el,core_pw_pos(xyz),jcwf,iproj)
          end do

       end do
    end do

    ! agrecocmplx: distinguish between real and complex case
    ! in complex case use complex copy to perform product
    if (loc_cmplx) then
       do xyz=1,3
          call sparse_copy(temp_cmplx_mat,core_pw_pos(xyz))
          call sparse_product(core_ngwf_pos(xyz),temp_cmplx_mat,ps_overlap)
       end do
    ! real case
    else
       do xyz=1,3
          call sparse_product(core_ngwf_pos(xyz),core_pw_pos(xyz),ps_overlap)
       end do
    end if

    call sparse_destroy(ps_overlap)

    call sparse_destroy(core_pw_ovlp)
    if (loc_cmplx) call sparse_destroy(temp_cmplx_mat) ! agrecocmplx
    do xyz=3,1,-1
       call sparse_destroy(core_pw_pos(xyz))
    end do

    call timer_clock('eels_core_paw_corr',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving eels_core_paw_corr'

  end subroutine eels_core_paw_corr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module eels
