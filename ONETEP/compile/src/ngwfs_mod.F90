
!================================================================!
!                                                                !
!                          NGWFs  module                         !
!                                                                !
! This module initialises and manipulates NGWFs.                 !
!----------------------------------------------------------------!
! Originally written by Chris-Kriton Skylaris in May 2000.       !
! Subsequent modifications and parallelisation by                !
! Chris-Kriton Skylaris.                                         !
! Automatic NGWF initialisation added on 31/03/2006 by           !
! Chris-Kriton Skylaris.                                         !
! Reorganisation and cleanup by Nicholas Hine in September 2010  !
!================================================================!
! Reference for Cartesian NGWFs:                                 !
! [1] "The Computational Modelling of Heavy Atom Chemistry",     !
!     Chris-Kriton Skylaris (PhD thesis), Cambridge, 1999.       !
! url: www.southampton.ac.uk/assets/centresresearch/documents/   !
!      compchem/skylaris_phd.pdf                                 !
!================================================================!


module ngwfs

  use constants, only: DP
  use ion, only: RADIAL_NGWF_TYPE
  use ngwf_data, only: GTO_SET

  implicit none
  private

  ! cks: public functions and subroutines
  public :: ngwfs_initialise
  public :: ngwfs_initialise_ngwfs
  public :: ngwfs_merge_sets
  public :: ngwfs_normalise
  public :: ngwfs_initialise_from_fragments

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwfs_initialise(rep, ngwf_basis, mdl, hfxstate, &
       restart_rootname, do_not_reexpand)

    !=======================================================================!
    ! This subroutine initialises the NGWFs of an NGWF_REP.                 !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! rep (in/out)              : The NGWF_REP whose NGWFs to initialise.   !
    ! ngwf_basis (in/out)       : Function basis describing NGWFs           !
    !                             in/out because radius may be overridden   !
    !                             in ngwfs_initialise_ngwfs().              !
    ! mdl (input)               : Needed for %elements.                     !
    ! restart_rootname (in, opt): @docme                                    !
    ! do_not_reexpand (in, opt) : If present and .true., SW expansion will  !
    !                             not be performed on freshly initialised   !
    !                             NGWFs.                                    !
    !-----------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 2016.10.14.                              !
    ! Modified for embedding by Robert Charlton, 10/08/2018.                !
    !=======================================================================!

    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE, hf_exchange_dkn_indep_stage
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, ngwf_rep_register_change
    use polarisable_embedding, only: polarisable_embedding_expand_ngwf_pairs
    use rundat, only: pub_pol_emb_qmstar, pub_use_hfx, pub_use_activehfx, &
         pub_active_region

    implicit none

    ! Arguments
    type(NGWF_REP), intent(inout)   :: rep
    type(FUNC_BASIS), intent(inout) :: ngwf_basis(rep%nsub)
    type(MODEL), intent(inout), target :: mdl
    type(HFX_STATE), intent(inout), target :: hfxstate
    character(len=*), intent(in), optional :: restart_rootname
    logical, intent(in), optional   :: do_not_reexpand

    ! Locals
    integer :: isub
    logical :: local_do_not_reexpand

    ! -------------------------------------------------------------------------

    if(present(do_not_reexpand)) then
       local_do_not_reexpand = do_not_reexpand
    else
       local_do_not_reexpand = .false.
    end if

    do isub=1,mdl%nsub
       call ngwfs_initialise_ngwfs(rep%ngwfs_on_grid(isub),ngwf_basis(isub), &
            mdl,isub,restart_rootname=restart_rootname)
    end do

    ! jd: *** NGWFS changed ***
    call ngwf_rep_register_change(rep,'ngwfs_initialise')

    do isub=1,mdl%nsub
       ! jd: Re-expand new NGWFs in SWs, if needed
       if(pub_pol_emb_qmstar .and. .not. local_do_not_reexpand) then
          call polarisable_embedding_expand_ngwf_pairs(rep, ngwf_basis(isub), &
               mdl)
       end if
       ! rc2013: if this is the active region re-expand NGWFs
       if((pub_use_hfx .or. (pub_use_activehfx .and. isub==pub_active_region)) &
            .and. .not. local_do_not_reexpand) then
          call hf_exchange_dkn_indep_stage(hfxstate, mdl, isub, &
               mdl%regions(isub)%par, rep, ngwf_basis(isub))
       end if
    end do

  end subroutine ngwfs_initialise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwfs_initialise_ngwfs(ngwfs_on_grid,ngwf_basis,mdl,ireg, &
       restart_rootname)

    !=======================================================================!
    ! This subroutine initialises the NGWFs of my_proc_id ppd-storage form. !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! ngwfs_on_grid (in/out)  : NGWF data in ppds on psinc grid             !
    ! ngwf_basis (in/out)     : Function basis describing NGWFs             !
    !                           in/out because radius may be overridden.    !
    ! mdl (input)             : Needed for %elements.                       !
    ! ireg (input)            : Identifier for NGWF region.                 !
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
    ! Modified by Andrea Greco on 10/05/2015 to use complex NGWFs.          !
    ! Modified by Jacek Dziedzic on 2016.10.14 to use NGWF_REP.             !
    ! pub_par global variable removed by Robert Charlton, 09/06/2017.       !
    ! Substantially modified for embedding by Robert Charlton, 13/07/2017.  !
    ! Further modified by Joseph Prentice, May 2018                         !
    !=======================================================================!

    ! agrecocmplx
    use datatypes, only: FUNCTIONS
    use comms, only: pub_my_proc_id, pub_on_root
    use constants, only: DP, stdout, NORMAL, EDA_ISOLATED, &
        EDA_FROZENNONIDEM, EDA_FROZENIDEM, EDA_POLFRAGLOC_DEVEL, &
        EDA_CTFRAGLOC_DEVEL, EDA_FULL, EDA_INIT
    use electronic_history, only: elec_history_compose_ngwfs
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use parallel_strategy, only: PARAL_INFO
    use restart, only: restart_ngwfs_tightbox_input, &
         restart_sph_waves_input
    use rundat, only: pub_read_tightbox_ngwfs, pub_read_sw_ngwfs, &
         cond_read_tightbox_ngwfs, pub_output_detail, &
         pub_devel_code, pub_dftb, mix_ngwfs_type, pub_full_rand_ngwf, &
         md_global_restart, pub_eda_mode, &
         pub_eda_reset_ngwfs_pol, pub_eda_reset_ngwfs_ct, &
         pub_eda_read_frags, pub_frag_file_prefix, &
         pub_eda_read_super, pub_super_file_prefix, &
         pub_rand_sigma, pub_rand_conv, pub_conv_func, &
         pub_frag_counter, pub_dmft_plot_real_space, pub_dmft_optics, &
         pub_dmft_points, pub_eda_continuation, &
         pub_eda_continuation_loaded_ngwfs, pub_rootname, pub_dmft_read

    implicit none

    ! Arguments
    type(FUNCTIONS), intent(inout) :: ngwfs_on_grid
    type(FUNC_BASIS), intent(inout) :: ngwf_basis
    type(MODEL), target, intent(inout) :: mdl
    character(len=*), intent(in), optional :: restart_rootname
    integer, intent(in)             :: ireg ! region counter

    ! Local Variables
    integer :: ingwf, loc_ingwf
    integer :: shell,em      ! shell counter for atomic set
    integer :: iat, orig_iat, z_at
    integer :: isp
    integer :: ingwf_on_atom, ingwf_in_species
    integer :: iset
    integer :: ierr
    logical :: loc_cond
    logical :: initialise_new
    logical :: check, check_nabla
    ! rc2013: hack to avoid meaningless changes
    type(PARAL_INFO), pointer :: par_point
    character(len=1) :: ireg_str
    character(len=30) :: file_suffix

    ! jcap: point at correct parallel strategy
    par_point => mdl%regions(ireg)%par

    ! ebl: check for existence of files storing ham, nabla
    if (pub_dmft_points>0) then
       inquire(file=trim(pub_rootname)//'.ham1'  , exist=check)
       inquire(file=trim(pub_rootname)//'.nabla', exist=check_nabla)
       if (check .and. .not. pub_dmft_plot_real_space .and. &
            .not. (pub_dmft_optics.and..not.check_nabla) .and. pub_dmft_read) &
            return
    end if

    ! ndmh: check for cond NGWFs
    if (ngwf_basis%name=='ngwfs') then
       iset = 1
       loc_cond = .false.
    else if (ngwf_basis%name=='ngwfs_cond') then
       iset = 2
       loc_cond = .true.
    else if (ngwf_basis%name=='ngwfs_aux') then
       iset = 3
       loc_cond = .false.
    else if (ngwf_basis%name=='hub_ngwfs') then
       iset = 4
       loc_cond = .false.
    else if (ngwf_basis%name=='pdos_sws') then
       iset = 1
       loc_cond = .false.
    end if

    ! ndmh: check if loading NGWFs in from file will be required
    ! ndmh: so we know whether to run the atomsolver (must match logic below)
    initialise_new = .true.
    if (loc_cond) then
       if (cond_read_tightbox_ngwfs) then
          initialise_new = .false.
       end if
    else ! not loc_cond
       if (pub_read_tightbox_ngwfs .and. mix_ngwfs_type.eq.'NONE') then
          initialise_new = .false.
       else if (mix_ngwfs_type.ne.'NONE') then
          initialise_new = .false.
       else if ((pub_eda_read_frags .and. pub_eda_mode == EDA_ISOLATED) .or. &
            pub_eda_mode == EDA_FROZENNONIDEM) then
          initialise_new = .false.
       else if (pub_read_sw_ngwfs) then
          initialise_new = .false.
       else if (pub_read_tightbox_ngwfs .and. md_global_restart) then
          initialise_new = .false.
       end if
    end if

    ! ndmh: Check for sphere radius overrides
    do loc_ingwf=1,ngwf_basis%num_on_proc(pub_my_proc_id)

       ! ndmh: find out global NGWF index, atom index, atom index in input-file
       ! ndmh: order, species number and NGWF index on this atom
       ingwf = loc_ingwf + ngwf_basis%first_on_proc(pub_my_proc_id) - 1
       iat = ngwf_basis%atom_of_func(ingwf)
       orig_iat = mdl%regions(ireg)%par%orig_atom(iat)
       isp = mdl%regions(ireg)%elements(orig_iat)%species_number
       ingwf_on_atom = ingwf - ngwf_basis%first_on_atom(iat) + 1

       ! cks: determine the shell number and the z-angular momentum component
       ! cks: of the current function
       ingwf_in_species = 0
       if (.not. pub_dftb) then
          do shell=1,mdl%regions(ireg)%radial_ngwfs(isp,iset)%nshells
             do em=-mdl%regions(ireg)%radial_ngwfs(isp,iset)%angmom(shell), &
                  mdl%regions(ireg)%radial_ngwfs(isp,iset)%angmom(shell)
                ingwf_in_species = ingwf_in_species + 1
                if (ingwf_in_species.eq.ingwf_on_atom) exit
             enddo
             if (ingwf_in_species.eq.ingwf_on_atom) exit
          enddo
       else ! ab: in DFTB, shell number is a part of dftb_par datastructure
          z_at = mdl%elements(orig_iat)%atomic_number
          shell= mdl%dftb_par%z(z_at)%ngwf_type(ingwf_on_atom)
       end if

       ! ndmh: override sphere radius if function is smaller than sphere radius
       if (ngwf_basis%spheres(loc_ingwf)%radius > &
            mdl%regions(ireg)%radial_ngwfs(isp,iset)%rc(shell)) then
          ngwf_basis%spheres(loc_ingwf)%radius = mdl%regions(ireg)%radial_ngwfs(isp,iset)%rc(shell)
       end if

    end do

    ! rc2013: attach regional string to file name if needed
    write(file_suffix,'(a30)') trim(ngwf_basis%name)
    if(mdl%nsub .gt. 1) then
       write(ireg_str,'(i1)') ireg
       write(file_suffix,'(a30)') adjustl(trim(file_suffix)//"_"//ireg_str)
    end if
    file_suffix = adjustl(trim(file_suffix))

    ! ndmh: load NGWFs in from file if required
    initialise_new = .false.
    if ((ngwf_basis%name=='hub_ngwfs').or.(ngwf_basis%name=='ngwfs_aux')&
         .or.(ngwf_basis%name=='pdos_sws')) then
       initialise_new = .true.
    else if (loc_cond .and. cond_read_tightbox_ngwfs) then
       ! cks: initialise NGWFs by reading values stored on disk in:
       ! cks: universal tightbox representation
       call restart_ngwfs_tightbox_input(ngwfs_on_grid, &    ! in/out
            ngwf_basis, mdl%cell, mdl%fftbox, mdl%regions(ireg)%elements, &
            'tightbox_'//trim(file_suffix), mdl%regions(ireg))
    else if (.not. loc_cond .and. (pub_read_tightbox_ngwfs .or. &
                                    present(restart_rootname)) ) then
       ! cks: initialise NGWFs by reading values stored on disk in
       ! cks: universal tightbox representation
       ! kkbd: Optionally restart from a specified rootname
       call restart_ngwfs_tightbox_input(ngwfs_on_grid, &    ! in/out
            ngwf_basis, mdl%cell, mdl%fftbox, mdl%regions(ireg)%elements, &
            'tightbox_'//trim(file_suffix), mdl%regions(ireg), &
            restart_rootname=restart_rootname)
    else if (.not. loc_cond .and. (pub_eda_read_frags .and. &
          pub_eda_mode == EDA_ISOLATED)) then
       ! mjsp: initialise fragment NGWFs by reading values stored on disk in
       ! mjsp: universal tightbox representation
       call restart_ngwfs_tightbox_input(ngwfs_on_grid, &    ! in/out
            ngwf_basis, mdl%cell, mdl%fftbox, mdl%regions(ireg)%elements, &
            'tightbox_'//trim(file_suffix), mdl%regions(ireg), &
            ref_dir=pub_frag_file_prefix(pub_frag_counter), ref_is_filename=.true. )
    else if (.not. loc_cond .and. pub_eda_read_super .and. &
               .not.((pub_eda_reset_ngwfs_pol .and. &
                    (pub_eda_mode == EDA_FROZENIDEM .or. &
                     pub_eda_mode == EDA_POLFRAGLOC_DEVEL)) .or. &
                   (pub_eda_reset_ngwfs_ct .and. &
                    (pub_eda_mode == EDA_CTFRAGLOC_DEVEL) .or. &
                     pub_eda_mode == EDA_FULL)) .and. &
               .not.(pub_eda_mode == EDA_INIT)) then
       ! mjsp: initialise supermolecule NGWFs by reading values stored on
       ! mjsp: disk in universal tightbox representation
       if (pub_eda_continuation .and. .not.(pub_eda_continuation_loaded_ngwfs)) then
          ! if continuation but not yet loaded the data,
          ! then use continuation data (in working directory)
          call restart_ngwfs_tightbox_input(ngwfs_on_grid, &    ! in/out
               ngwf_basis, mdl%cell, mdl%fftbox, mdl%regions(ireg)%elements, &
               'tightbox_'//trim(file_suffix), mdl%regions(ireg))
          pub_eda_continuation_loaded_ngwfs = .true.
       else
          ! use frozen density (from preptool calculation,
          ! optionally in separate directory)
          call restart_ngwfs_tightbox_input(ngwfs_on_grid, &    ! in/out
               ngwf_basis, mdl%cell, mdl%fftbox, mdl%regions(ireg)%elements, &
               'tightbox_'//trim(file_suffix), mdl%regions(ireg), &
               pub_super_file_prefix, ref_is_filename=.true. )
       end if
    else if (.not. loc_cond .and. pub_read_sw_ngwfs) then
       ! ars: initialise NGWFs by reading values stored on disk in
       ! ars: spherical waves representation
       call restart_sph_waves_input(ngwfs_on_grid, &     ! output
            ngwf_basis, mdl%cell, mdl%fftbox, mdl%uni_tightbox, &
            'sw_'//trim(file_suffix), mdl%regions(ireg)) ! input
    else
       initialise_new = .true.
    end if

    ! ndmh: NGWFs were not loaded in, so create new ones
    if (initialise_new) then
       if (pub_on_root.and.(pub_output_detail>=NORMAL)) then
          if (pub_full_rand_ngwf) then
             write(stdout,'(a)',advance='no') 'Fully random '
          else if (pub_rand_sigma>0.0_DP) then
             write(stdout,'(a)',advance='no') 'Randomly perturbed '
          end if
          if (ngwf_basis%name=='ngwfs') then
             if (mdl%nsub.eq.1) then
                write(stdout,'(a)',advance='no') 'NGWF initialisation ...'
             else
                write(stdout,'(a,i0,a)',advance='no') 'Region '//&
                     &trim(ireg_str)//' NGWF initialisation ...'
             end if
          else if (ngwf_basis%name=='ngwfs_cond') then
             if (mdl%nsub.eq.1) then
                write(stdout,'(a)',advance='no') 'Conduction NGWF &
                     &initialisation ...'
             else
                write(stdout,'(a,i0,a)',advance='no') 'Region '//&
                     &trim(ireg_str)//' conduction NGWF initialisation ...'
             end if
          else if (ngwf_basis%name=='ngwfs_aux') then
             if (mdl%nsub.eq.1) then
                write(stdout,'(a)',advance='no') 'Auxiliary NGWF &
                     &initialisation ...'
             else
                write(stdout,'(a,i0,a)',advance='no') 'Region '//&
                     &trim(ireg_str)//' auxiliary NGWF initialisation ...'
             end if
          end if

          ! agreco: write information about convolution scheme used
          ! currently only the use of erfc is supported
          if (pub_full_rand_ngwf.and.pub_rand_conv) then
             if (pub_conv_func/='erfc') then
                write(stdout,'(a)') 'WARNING: unsupported convolution function &
                   &specified in input file: convolution of randomly &
                   &initialised NGWFs will be ignored.'
             else
                write(stdout,'(a,a)') 'Initial random NGWFs values will be convoluted &
                   &using function', pub_conv_func
             end if
          end if
       end if

       ! ndmh: Expand radial NGWFs to psinc grid
       if (index(pub_devel_code,'INIT_NGWFS_RECIP')>0) then
          ! agreco: currently randomly initialised NGWFs are not compatible
          ! with routine ngwfs_initialise_from_recip
          if (pub_full_rand_ngwf.or.(pub_rand_sigma>0.0_DP)) then
             if (pub_on_root) then
                write(stdout,'(a)') 'WARNING: randomly initialised NGWFs &
                &not yet compatible with INIT_NGWFS_RECIP: using default case...'
             end if
          end if
          call ngwfs_initialise_from_recip(ngwfs_on_grid,ngwf_basis, &
               mdl%regions(ireg)%radial_ngwfs(:,iset),mdl%regions(ireg)%elements, &
               mdl%fftbox,mdl%cell,mdl%regions(ireg)%par)
       else
          call ngwfs_initialise_from_radial(ngwfs_on_grid,ngwf_basis, &
               mdl%regions(ireg)%radial_ngwfs(:,iset),mdl%regions(ireg)%elements, &
               mdl%cell,mdl%fftbox,mdl%regions(ireg)%par,mdl%dftb_par)
       end if

       if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
            write(stdout,'(a)') '... done'

    end if


    ! smmd: Improve new NGWFs using electronic history stored in
    ! smmd: universal tightbox representation
    if (.not. loc_cond ) then
       call elec_history_compose_ngwfs(ngwfs_on_grid, ngwf_basis, &
            mdl%cell, mdl%fftbox, mdl%regions(ireg)%elements, &
            ireg, ierr, par_point)
    end if

  end subroutine ngwfs_initialise_ngwfs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_initialise_from_radial(ngwfs_on_grid, ngwf_basis, &
      radial_ngwfs,elements,cell,fftbox,par, dftb_par)

    !=======================================================================!
    ! This subroutine initialises the NGWFs of my_proc_id ppd-storage form, !
    ! according to the radial_ngwfs defined in the parent routine.          !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! ngwfs_on_grid (output)  : NGWF data in ppds on psinc grid             !
    ! ngwf_basis (input)      : function basis describing NGWFs             !
    ! elements (input)        : Array of ELEMENT type with data about ions  !
    ! radial_ngwfs (input)    : Array of RADIAL_NGWF type with radial part  !
    !                           f(r) for each species                       !
    ! dftb_par (input)        : Data structure needed in DFTB calculation.  !
    !-----------------------------------------------------------------------!
    ! Spun off from parent routine by Nicholas Hine, October 2011.          !
    ! Modified by Andrea Greco on 10/05/2015 to allow use of complex NGWFs. !
    ! Modified by Jacek Dziedzic on 2016.10.14 to use NGWF_REP.             !
    ! Extended to Cartesian NGWFs used in DFTB by Arihant Bhandari, Nov 2021!
    !=======================================================================!

    use basis, only: basis_ppd_location
    use comms, only: pub_my_proc_id, pub_on_root
    use constants, only: DP, PI, stdout
    use datatypes, only: FUNCTIONS, data_set_to_zero
    use dftb_parameters, only: DFTB_PARAMS
    use geometry, only: POINT, operator(*), operator(+), local_displacement, &
        operator(.DOT.), operator(-)
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use ion, only: element
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_full_rand_ngwf, pub_rand_seed_ngwf, pub_rand_sigma, &
        pub_rand_conv, pub_conv_func, pub_conv_region_width, &
        pub_ngwfs_phase, pub_ngwfs_rand_phase, pub_dftb, pub_dftb_cartesian_ngwfs
!$  use rundat, only: pub_threads_max
    use simulation_cell, only : CELL_INFO
    use spherical_wave, only: sw_init, sw_exit
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNCTIONS), intent(inout) :: ngwfs_on_grid
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(PARAL_INFO), intent(inout) :: par
    type(ELEMENT), intent(in) :: elements(par%nat)
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(RADIAL_NGWF_TYPE), intent(in) :: radial_ngwfs(par%num_species)
    type (DFTB_PARAMS), intent(in) :: dftb_par

    ! Local Variables
    integer :: point_counter
    integer :: ppd,ppd_loc,loc_1,loc_2,loc_3
    integer :: a1_neighbour,a2_neighbour,a3_neighbour
    integer :: a1_upper,a2_upper,a3_upper
    integer :: a1_lower,a2_lower,a3_lower
    integer :: a1_cell,a2_cell,a3_cell
    integer :: ingwf, loc_ingwf
    integer :: ingwf_on_atom, ingwf_in_species
    integer :: iat, orig_iat, z_at
    integer :: shell         ! shell counter for atomic set
    integer :: ll            ! ab: angular momentum for shell
    integer :: mx, my, mz    ! ab: components of angular momentum (ll=mx+my+mz)
    integer :: ppd_count     ! ndmh: sphere ppd list counter
    integer :: isp
    integer :: lmax
    logical :: init_sw
    real(kind=DP) :: local_1, local_2, local_3
    type(POINT) :: current_point,current_displacement,periodic_centre
    ! agreco: local variables for random initialisation
    logical :: any_rand
    integer :: rand_size, start, finish
    integer, allocatable, dimension(:) :: rand_seed
    integer :: ierr
    real(kind=DP) :: rtime
    character(len=10) :: sys_time
    real(kind=DP) :: norm_constant
    real(kind=DP) :: dist
    real(kind=DP) :: factor
    real(kind=DP) :: recip_twopi
    real(kind=DP) :: fcoord_cp(3)
    real(kind=DP) :: fcoord_centre(3)
    ! agreco
    real(kind=DP) :: rand_phase


    ! agrecocmplx: extra variables to initialise complex NGWFs
    real(kind=DP) :: ngwfs_grid_real, ngwfs_grid_imag

    ! Start timer
    call timer_clock('ngwfs_initialise_from_radial',1)

    ! agreco: initialise random seed if we have some random
    ! initialisation in the NGWFs
    any_rand = (pub_full_rand_ngwf.or.(pub_rand_sigma>0.0_DP)&
                .or.pub_ngwfs_rand_phase)

    recip_twopi = 0.5_DP / pi

    !seed the random number generator using the system time
    if (any_rand) then
       call random_seed(size=rand_size)
       allocate(rand_seed(1:rand_size),stat=ierr)
       call utils_alloc_check('ngwfs_initialise_from_radial','rand_seed',ierr)
       ! agreco: use seed specified by user in input file if present
       if (pub_rand_seed_ngwf>=0) then
          rand_seed=pub_rand_seed_ngwf
       ! agreco: if not, initialise seed from current date and time
       else
          call date_and_time(time=sys_time)
          read(sys_time,*) rtime
          rand_seed=int(rtime*1000.0_DP)
       end if
       call random_seed(put=rand_seed)
       deallocate(rand_seed,stat=ierr)
       call utils_dealloc_check('ngwfs_initialise_from_radial','rand_seed',ierr)
    end if

    ! agrecocmplx: either global or random phase
    if (pub_ngwfs_rand_phase.and.(pub_ngwfs_phase /= 0.0_DP)) then
       call utils_abort('Error: cannot use both random and non-zero NGWFs phases')
    end if

    ! agrecocmplx: ignore global phase in case of real NGWFs
    if (pub_ngwfs_phase /= 0.0_DP) then
       if (ngwfs_on_grid%iscmplx) then
          if (pub_on_root) write(stdout,'(a,a,f8.5,a)') &
             new_line('a'), 'Multiplying complex NGWFs by phase &
             &exp(i*PI*', pub_ngwfs_phase, ')...'
       else
          if (pub_on_root) write(stdout,'(a,a)') &
             new_line('a'), 'WARNING: cannot apply complex phase to real &
             &NGWFs: flag MULT_NGWF_BY_PHASE will be ignored.'
       end if
    end if

    ! agrecocmplx: ignore random phase in case of real NGWFs
    if (pub_ngwfs_rand_phase) then
       if (ngwfs_on_grid%iscmplx) then
          if (pub_on_root) write(stdout,'(a)') &
             new_line('a'), 'Multiplying complex NGWFs by random phase...'
       else
          if (pub_on_root) write(stdout,'(a,a)') &
             new_line('a'), 'WARNING: cannot apply complex phase to real &
             &NGWFs: flag MULT_NGWF_BY_RANDOM_PHASE will be ignored.'
       end if
    end if

    ! ndmh: to check for ppds included multiple times, when cell=fftbox
    ! agreco: not necessary when we initialise to fully random values
    ! but necessary when using convolution which depends on distance from
    ! the centre of the function (CHECK THIS!!!!)
    if ((fftbox%coin1).and.((.not.pub_full_rand_ngwf).or.pub_rand_conv)) then
       a1_lower = -1; a1_upper = 1
    else
       a1_lower =  0; a1_upper = 0
    end if
    if ((fftbox%coin2).and.((.not.pub_full_rand_ngwf).or.pub_rand_conv)) then
       a2_lower = -1; a2_upper = 1
    else
       a2_lower =  0; a2_upper = 0
    end if
    if ((fftbox%coin3).and.((.not.pub_full_rand_ngwf).or.pub_rand_conv)) then
       a3_lower = -1; a3_upper = 1
    else
       a3_lower =  0; a3_upper = 0
    end if

    ! ndmh: for high-l
    init_sw = .false.
    lmax = 0
    do isp=1,par%num_species
       lmax = max(lmax,maxval(radial_ngwfs(isp)%angmom(:)))
    end do
    if (lmax>=5) init_sw = .true.
    if (init_sw) call sw_init(lmax,lmax*2+1)

    !ngwfs_on_grid(:) = 0.0_DP
    ! agrecocmplx: initialise according to real or complex
    call data_set_to_zero(ngwfs_on_grid)

    ! ndmh: Loop over NGWFs on this proc

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(loc_ingwf,ingwf,iat,orig_iat,z_at,isp,ingwf_on_atom,norm_constant, &
!$OMP      ingwf_in_species,shell,ll,mx,my,mz,ppd_count,ppd,ppd_loc,a1_cell,a2_cell, &
!$OMP      a3_cell,a1_neighbour,a2_neighbour,a3_neighbour,periodic_centre, &
!$OMP      point_counter,local_1,local_2,local_3,current_displacement, &
!$OMP      current_point,start,finish,ngwfs_grid_real,ngwfs_grid_imag,factor, &
!$OMP      fcoord_cp,fcoord_centre,dist,rand_phase) &
!$OMP SHARED(a1_upper,a1_lower,a2_upper,a2_lower,a3_upper,a3_lower, &
!$OMP      radial_ngwfs,ngwf_basis,pub_my_proc_id,par,elements,dftb_par,fftbox,cell, &
!$OMP      ngwfs_on_grid,pub_full_rand_ngwf,any_rand,pub_rand_sigma,pub_threads_max, &
!$OMP      pub_rand_conv,pub_conv_func,pub_conv_region_width,recip_twopi, &
!$OMP      pub_ngwfs_phase,pub_ngwfs_rand_phase,pub_dftb,pub_dftb_cartesian_ngwfs)

!$OMP DO
    do loc_ingwf=1,ngwf_basis%num_on_proc(pub_my_proc_id)

       ! ndmh: find out global NGWF index, atom index, atom index in input-file
       ! ndmh: order, species number and NGWF index on this atom
       ingwf = loc_ingwf + ngwf_basis%first_on_proc(pub_my_proc_id) - 1
       iat = ngwf_basis%atom_of_func(ingwf)
       orig_iat = par%orig_atom(iat)
       isp = elements(orig_iat)%species_number
       ingwf_on_atom = ingwf - ngwf_basis%first_on_atom(iat) + 1

       norm_constant = 0.0_DP
       ! cks: determine the shell number and the angular momentum component
       ! cks: of the current function
       ingwf_in_species = 0
       ! ab: for Sphericals NGWFs, only z-angular momentum is relevant
       if (.not. (pub_dftb .and. pub_dftb_cartesian_ngwfs)) then
          do shell=1,radial_ngwfs(isp)%nshells
             do mz=-radial_ngwfs(isp)%angmom(shell),radial_ngwfs(isp)%angmom(shell)
                ingwf_in_species = ingwf_in_species + 1
                if (ingwf_in_species.eq.ingwf_on_atom) exit
             enddo
             if (ingwf_in_species.eq.ingwf_on_atom) exit
          enddo
       ! ab: For Cartesian NGWFs in DFTB, there are three angular momentum
       ! ab: components, mx, my and mz which sum up to the angular momentum ll
       ! ab: For each ingwf_on_atom, we determine the unique combination of
       ! ab: mx, my, mz and the corresponding shell
       else
          z_at = elements(orig_iat)%atomic_number
          do shell= 1, dftb_par%z(z_at)%nshells
             ll   = dftb_par%z(z_at)%ang(shell)
             do mx = 0, ll
                do my = 0, ll
                   do mz = 0, ll
                      if (mx+my+mz .ne. ll) then
                         cycle
                      else
                         ingwf_in_species = ingwf_in_species + 1
                      end if
                      if (ingwf_in_species.eq.ingwf_on_atom) exit
                   end do
                   if (ingwf_in_species.eq.ingwf_on_atom) exit
                end do
                if (ingwf_in_species.eq.ingwf_on_atom) exit
             end do
             if (ingwf_in_species.eq.ingwf_on_atom) exit
          end do
       end if

       ! ndmh: loop over all ppds in the NGWF sphere
       do ppd_count=1,ngwf_basis%spheres(loc_ingwf)%n_ppds_sphere

          ppd = ngwf_basis%spheres(loc_ingwf)%ppd_list(1,ppd_count)
          ppd_loc = ngwf_basis%spheres(loc_ingwf)%ppd_list(2,ppd_count)

          do a1_cell=a1_lower,a1_upper
             do a2_cell=a2_lower,a2_upper
                do a3_cell=a3_lower,a3_upper

                   ! ndmh: find which copy of the supercell this NGWF comes from
                   a1_neighbour=nint(real(ppd_loc,kind=DP)/9.0_DP)
                   a2_neighbour=nint(real(ppd_loc-9*a1_neighbour,DP)/3.0_DP)
                   a3_neighbour=ppd_loc-9*a1_neighbour-3*a2_neighbour
                   if (fftbox%coin1) a1_neighbour = a1_cell
                   if (fftbox%coin2) a2_neighbour = a2_cell
                   if (fftbox%coin3) a3_neighbour = a3_cell

                   periodic_centre=ngwf_basis%spheres(loc_ingwf)%centre &
                        + real(a1_neighbour,DP)*cell%a1 &
                        + real(a2_neighbour,DP)*cell%a2 &
                        + real(a3_neighbour,DP)*cell%a3

                   ! agreco: if partially extended NGWFs are used,
                   ! we need the distance from the NGWF centre along the
                   ! localised directions
                   if (pub_full_rand_ngwf.and.pub_rand_conv) then
                      fcoord_centre(1) = periodic_centre .dot. cell%b1
                      fcoord_centre(2) = periodic_centre .dot. cell%b2
                      fcoord_centre(3) = periodic_centre .dot. cell%b3
                      fcoord_centre = fcoord_centre * recip_twopi
                      fcoord_centre = modulo(fcoord_centre,1.0_DP)

                      if (ngwf_basis%spheres(loc_ingwf)%extended(1)) then
                         periodic_centre = periodic_centre - fcoord_centre(1)*cell%a1
                      end if
                      if (ngwf_basis%spheres(loc_ingwf)%extended(2)) then
                         periodic_centre = periodic_centre - fcoord_centre(2)*cell%a2
                      end if
                      if (ngwf_basis%spheres(loc_ingwf)%extended(3)) then
                         periodic_centre = periodic_centre - fcoord_centre(3)*cell%a3
                      end if

                   end if

                   point_counter = ngwf_basis%spheres(loc_ingwf)%offset + &
                        cell%n_pts*(ppd_count-1)

                   ! cks: loop over actual grid points in the current ppd
                   do loc_3=0,cell%n_pt3-1
                      local_3=real(loc_3,DP)*cell%d3

                      do loc_2=0,cell%n_pt2-1
                         local_2=real(loc_2,DP)*cell%d2

                         do loc_1= 0, cell%n_pt1-1
                            local_1=real(loc_1,DP)*cell%d1

                            current_displacement = local_displacement( &
                                 cell%a1_unit,cell%a2_unit, &
                                 cell%a3_unit,local_1,local_2,local_3)

                            current_point = basis_ppd_location(ppd,cell) + &
                                 current_displacement

                            ! agreco: if partially extended NGWFs are used,
                            ! we need the distance from the NGWF centre along the
                            ! localised directions
                            if (pub_full_rand_ngwf.and.pub_rand_conv) then
                               fcoord_cp(1) = current_point .dot. cell%b1
                               fcoord_cp(2) = current_point .dot. cell%b2
                               fcoord_cp(3) = current_point .dot. cell%b3
                               fcoord_cp = fcoord_cp * recip_twopi
                               fcoord_cp = modulo(fcoord_cp,1.0_DP)

                               if (ngwf_basis%spheres(loc_ingwf)%extended(1)) then
                                  current_point = current_point - fcoord_cp(1)*cell%a1
                                  end if
                               if (ngwf_basis%spheres(loc_ingwf)%extended(2)) then
                                  current_point = current_point - fcoord_cp(2)*cell%a2
                               end if
                               if (ngwf_basis%spheres(loc_ingwf)%extended(3)) then
                                  current_point = current_point - fcoord_cp(3)*cell%a3
                               end if

                            end if

                            !if (pub_full_rand_ngwf) then
                            !   ! initialise grid point to
                            !   ! random value between -0.5 and 0.5
                            !   call random_number(ngwfs_on_grid(point_counter))
                            !   ngwfs_on_grid(point_counter) = &
                            !     ngwfs_on_grid(point_counter) - 0.5_DP
                            !   ! compute norm for normalisation
                            !   norm_constant = norm_constant + &
                            !     ngwfs_on_grid(point_counter)*ngwfs_on_grid(point_counter)
                            !else
                            !   ngwfs_on_grid(point_counter) = &
                            !        ngwfs_on_grid(point_counter) + &
                            !        ngwfs_at_point(shell,em,radial_ngwfs(isp), &
                            !        periodic_centre,current_point, &
                            !        ngwf_basis%spheres(loc_ingwf)%radius)
                            !end if

                             !agrecocmplx: distinguish between real and complex
                            if (pub_full_rand_ngwf) then
                               ! initialise grid point to
                               ! random value between -0.5 and 0.5
                               call random_number(ngwfs_grid_real)
                               ngwfs_grid_real = ngwfs_grid_real - 0.5_DP
                               if (ngwfs_on_grid%iscmplx) then
                                   call random_number(ngwfs_grid_imag)
                                   ngwfs_grid_imag = ngwfs_grid_imag - 0.5_DP
                                   ngwfs_on_grid%z(point_counter) = &
                                       cmplx(ngwfs_grid_real,ngwfs_grid_imag,kind=DP)
                                   ! agreco: apply convolution if requested
                                   if (pub_rand_conv) then
                                      ngwfs_on_grid%z(point_counter) = &
                                         ngwfs_on_grid%z(point_counter)*convolution_factor( &
                                         periodic_centre, current_point, &
                                         ngwf_basis%spheres(loc_ingwf)%radius, pub_conv_func, &
                                         region_width=pub_conv_region_width)
                                   end if
                                   ! compute norm for normalisation
                                   norm_constant = norm_constant + &
                                       ngwfs_grid_real*ngwfs_grid_real + &
                                       ngwfs_grid_imag*ngwfs_grid_imag
                               else
                                   ngwfs_on_grid%d(point_counter) = ngwfs_grid_real
                                   ! agreco: apply convolution if requested
                                   if (pub_rand_conv) then
                                      ngwfs_on_grid%d(point_counter) = &
                                         ngwfs_on_grid%d(point_counter)*convolution_factor( &
                                         periodic_centre, current_point, &
                                         ngwf_basis%spheres(loc_ingwf)%radius, pub_conv_func, &
                                         region_width=pub_conv_region_width)
                                   end if
                                   norm_constant = norm_constant + &
                                       ngwfs_grid_real*ngwfs_grid_real
                               end if
                            else
                               ! agreco: perturb the initial values introducing
                               ! a small randomisation with normal distribution
                               if (pub_rand_sigma>0.0_DP) then
                                  factor = random_normal(0.0_DP,pub_rand_sigma)
                               else
                                  factor = 0.0_DP
                               end if
                               ngwfs_grid_real = &
                                    ngwfs_at_point(shell,mx,my,mz,radial_ngwfs(isp), &
                                    periodic_centre,current_point, &
                                    ngwf_basis%spheres(loc_ingwf)%radius)
                               if (ngwfs_on_grid%iscmplx) then
                                   ngwfs_on_grid%z(point_counter) = &
                                       ngwfs_on_grid%z(point_counter) + &
                                       cmplx(ngwfs_grid_real*(1.0_DP+factor),factor,kind=DP)
                                   if (pub_rand_sigma>0.0_DP) then
                                      norm_constant = norm_constant + &
                                         ngwfs_grid_real*ngwfs_grid_real + &
                                         ngwfs_grid_imag*ngwfs_grid_imag
                                   end if
                               else
                                   ngwfs_on_grid%d(point_counter) = &
                                       ngwfs_on_grid%d(point_counter) + &
                                       ngwfs_grid_real*(1.0_DP+factor)
                                   if (pub_rand_sigma>0.0_DP) then
                                      norm_constant = norm_constant + &
                                         ngwfs_grid_real*ngwfs_grid_real
                                   end if
                               end if
                            end if

                            point_counter = point_counter + 1

                         enddo  ! loc_1
                      enddo  ! loc_2
                   enddo  ! loc_3

                enddo  ! a3_cell
             enddo  ! a2_cell
          enddo  ! a1_cell

       enddo  ! ppd_count

       ! agrecocmplx: apply global phase if requested
       ! multiply by exp(i theta) where theta is real
       if (pub_ngwfs_phase /= 0.0_DP) then
          if (ngwfs_on_grid%iscmplx) then
             start = ngwf_basis%spheres(loc_ingwf)%offset
             finish = point_counter - 1
             ngwfs_on_grid%z(start:finish) = &
             ngwfs_on_grid%z(start:finish) * &
             exp(cmplx(0.0_DP,1.0_DP,kind=DP)*pub_ngwfs_phase)
          end if
       end if

       ! agrecocmplx: apply random phase to each NGWF if requested
       ! multiply by exp(i theta) where theta is real and random,
       ! between 0 and 2PI
       if (pub_ngwfs_rand_phase) then
          if (ngwfs_on_grid%iscmplx) then
             start = ngwf_basis%spheres(loc_ingwf)%offset
             finish = point_counter - 1
             ! agreco: between 0 and 1
             call random_number(rand_phase)
             ! agreco: between 0 and 2PI
             rand_phase = rand_phase * 2.0_DP * pi
             ngwfs_on_grid%z(start:finish) = &
             ngwfs_on_grid%z(start:finish) * &
             exp(cmplx(0.0_DP,1.0_DP,kind=DP)*rand_phase)
          end if
       end if

       ! agreco: (re)normalise basis set, but not when
       ! applying random phase since its modulus squared is 1
       if (any_rand.and.(.not.pub_ngwfs_rand_phase)) then
          norm_constant = sqrt(norm_constant*cell%weight)
          start = ngwf_basis%spheres(loc_ingwf)%offset
          finish = point_counter - 1
          if (norm_constant > 0.0_DP) then
             !ngwfs_on_grid(start:finish) = &
             !   ngwfs_on_grid(start:finish)/norm_constant
             ! agrecocmplx
             if (ngwfs_on_grid%iscmplx) then
                 ngwfs_on_grid%z(start:finish) = &
                 ngwfs_on_grid%z(start:finish)/norm_constant
             else
                 ngwfs_on_grid%d(start:finish) = &
                 ngwfs_on_grid%d(start:finish)/norm_constant
             end if
          else
             call utils_abort('Error in ngwfs_initialise_from_radial:&
                & Invalid norm_constant')
          end if
       end if
    end do  ! loc_ingwf
!$OMP END DO

!$OMP END PARALLEL

    if (init_sw) call sw_exit

    ! Stop timer
    call timer_clock('ngwfs_initialise_from_radial',2)

  end subroutine ngwfs_initialise_from_radial


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_initialise_from_recip(ngwfs_on_grid,ngwf_basis, &
      radial_ngwfs,elements,fftbox,cell,par)

    !=======================================================================!
    ! This subroutine initialises the NGWFs of my_proc_id ppd-storage form !
    ! according to the radial_ngwfs defined in the parent routine, by first !
    ! transforming them to reciprocal space as projectors, building them in !
    ! the reciprocal space FFTbox, and then extracting them to the ppds.    !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! ngwfs_on_grid (output)  : NGWF data in ppds on psinc grid             !
    ! ngwf_basis (input)      : function basis describing NGWFs             !
    ! elements (input)        : Array of ELEMENT type with data about ions  !
    ! radial_ngwfs (input)    : Array of RADIAL_NGWF type with radial part  !
    !                           f(r) for each species                       !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine, November 2011.                              !
    ! Modified by Andrea Greco on 10/05/2015 to allow use of complex NGWFs. !
    !=======================================================================!

    use basis, only: basis_clean_function
    use constants, only: DP
    use datatypes, only: FUNCTIONS
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use ion, only: element
    use parallel_strategy, only: PARAL_INFO
    use projectors, only: PROJECTOR_SET, projectors_allocate_set, &
         projectors_deallocate_set, projectors_create_real, &
         projectors_init_fftbox_recip, projectors_exit_fftbox_recip
    use services, only: services_regular_transform
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: ngwf_basis
    !real(kind=DP), intent(out) :: ngwfs_on_grid(ngwf_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(inout) :: ngwfs_on_grid
    type(PARAL_INFO), intent(inout) :: par
    type(ELEMENT), intent(in) :: elements(par%nat)
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(RADIAL_NGWF_TYPE), intent(in) :: radial_ngwfs(par%num_species)

    ! Local Variables
    type(PROJECTOR_SET) :: ngwf_projectors
    integer :: ingwf
    integer :: iat, orig_iat
    integer :: shell      ! shell counter for atomic set
    integer :: isp
    integer :: proj_count
    integer, parameter :: n_recip_pts = 2001
    real(kind=DP), parameter :: g_max = 50.0_DP
    ! agrecocmplx: variable for complex NGWFs
    logical :: loc_cmplx

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid%iscmplx

    ! ndmh: allocate storage for ngwf projectors
    ngwf_projectors%n_proj_species = par%num_species
    ! agrecocmplx
    call projectors_allocate_set(ngwf_projectors, &
         maxval(radial_ngwfs(:)%nshells),n_recip_pts,par=par,is_cmplx=loc_cmplx)

    ! ndmh: set species_num_proj and species_first_proj values
    proj_count = 1
    do isp=1,ngwf_projectors%n_proj_species
       ! Get an example atom of this species
       do iat=1,par%nat
          orig_iat = par%orig_atom(iat)
          if (elements(orig_iat)%species_number==isp) exit
       end do
       ! Set number of projectors of each species
       ngwf_projectors%species_num_proj(isp) = ngwf_basis%num_on_atom(iat)
       ngwf_projectors%species_first_proj(isp) = proj_count
       ngwf_projectors%gmax(isp) = g_max
       ngwf_projectors%n_rad_pts(isp) = n_recip_pts
       ngwf_projectors%num_shells(isp) = radial_ngwfs(isp)%nshells
       ngwf_projectors%ang_mom(:,isp) = 0
       ngwf_projectors%rad_proj_recip(:,:,isp) = 0.0_DP
       do shell=1,radial_ngwfs(isp)%nshells
          ngwf_projectors%ang_mom(shell,isp) = radial_ngwfs(isp)%angmom(shell)
          call services_regular_transform(radial_ngwfs(isp)%angmom(shell),2, &
               radial_ngwfs(isp)%npts,maxval(radial_ngwfs(isp)%rad), &
               n_recip_pts,g_max,radial_ngwfs(isp)%func_real(:,shell), &
               ngwf_projectors%rad_proj_recip(1:n_recip_pts,shell,isp))
       end do
       proj_count = proj_count + ngwf_projectors%species_num_proj(isp)

    end do

    ! ndmh: copy in projector centres and radii
    do iat=1,par%nat
       orig_iat = par%orig_atom(iat)
       isp = elements(orig_iat)%species_number
       ngwf_projectors%proj_centre(iat) = elements(orig_iat)%centre
       ngwf_projectors%proj_max_radius(iat) = &
            maxval(ngwf_basis%spheres(:)%radius)
       ngwf_projectors%proj_species(iat) = isp
    end do

    ! ndmh: initialise projectors in reciprocal-space fftbox representation
    call projectors_init_fftbox_recip(ngwf_projectors,cell,fftbox)

    ! ndmh: create projectors in recip space & extract from FFTbox to ppds
    !call projectors_create_real(ngwf_basis,ngwf_projectors, &
    !     fftbox,cell,ngwfs_on_grid)
    ! agrecocmplx: call modified routine in projectors_mod
    call projectors_create_real(ngwf_basis,ngwf_projectors, &
         fftbox,cell,ngwfs_on_grid,is_cmplx=loc_cmplx)

    ! ndmh: zero regions outside sphere boundary but inside nonzero ppds
    do ingwf=1,ngwf_basis%proc_num
       !call basis_clean_function(ngwfs_on_grid,ngwf_basis%spheres(ingwf), &
       !     ngwf_basis%n_ppds,cell,fftbox)
       call basis_clean_function(ngwfs_on_grid,ngwf_basis%spheres(ingwf), &
            cell,fftbox)
    end do

    call projectors_exit_fftbox_recip(ngwf_projectors)
    call projectors_deallocate_set(ngwf_projectors)

  end subroutine ngwfs_initialise_from_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_merge_sets(ngwfs_on_grid,ngwf_basis,cell, &
       ngwfs_on_grid_src1,ngwf_basis_src1,ngwfs_on_grid_src2,ngwf_basis_src2, &
       par)

    !====================================================================!
    ! This subroutine merges two sets of NGWFs together to form a joint  !
    ! set.                                                               !
    !====================================================================!
    ! Arguments:                                                         !
    !   ngwfs_on_grid (output) : Merged NGWFs in ppd storage             !
    !   ngwf_basis (input)     : Spheres etc describing merged set       !
    !   ngwfs_on_grid_src1 (in): NGWF set 1 to merge (in ppds)           !
    !   ngwf_basis_src1 (in)   : Spheres etc describing NGWF set 1       !
    !   ngwfs_on_grid_src2 (in): NGWF set 2 to merge (in ppds)           !
    !   ngwf_basis_src2 (in)   : Spheres etc describing NGWF set 2       !
    !   par (in)               : Parallel strategy for the NGWF sets.    !
    !====================================================================!
    ! Written by Nicholas Hine in April 2011.                            !
    ! Modified to remove pub_par by Joseph Prentice, August 2018         !
    !====================================================================!

    use datatypes, only: FUNCTIONS
    use comms, only: pub_my_proc_id
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: PARAL_INFO
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort, utils_assert

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNCTIONS), intent(inout) :: ngwfs_on_grid
    type(CELL_INFO), intent(in) :: cell
    type(FUNC_BASIS), intent(in) :: ngwf_basis_src1
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid_src1
    type(FUNC_BASIS), intent(in) :: ngwf_basis_src2
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid_src2
    type(PARAL_INFO), intent(in) :: par

    ! Local Variables
    integer :: ifunc_src1, ifunc_src2
    integer :: iat, loc_iat
    integer :: count, count_src1, count_src2
    integer :: start, start_src, finish, finish_src
    character(len=10) :: iat_str
    logical :: loc_cmplx

    call utils_assert((ngwfs_on_grid%iscmplx .eqv. ngwfs_on_grid_src1%iscmplx) &
         .and. (ngwfs_on_grid%iscmplx .eqv. ngwfs_on_grid_src2%iscmplx), &
         'Error in ngwfs_merge_sets: you cannot mix real and complex NGWFs.')
    loc_cmplx = ngwfs_on_grid%iscmplx

    ! ndmh: loop over all atoms on proc and copy the spheres of each one
    count = 0
    count_src1 = 0
    count_src2 = 0
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1

       if ((ngwf_basis_src1%num_on_atom(iat) + &
            ngwf_basis_src2%num_on_atom(iat)) &
            /=ngwf_basis%num_on_atom(iat)) then
          write(iat_str,'(i10)') iat
          iat_str = adjustl(trim(iat_str))
          call utils_abort('Error in ngwfs_merge_sets: &
               &Mismatching NGWF counts on atom '//iat_str)
       end if

       do ifunc_src1=1,ngwf_basis_src1%num_on_atom(iat)
          count_src1 = count_src1 + 1
          count = count + 1

          ! Find start and finish in dest array
          start = ngwf_basis%spheres(count)%offset
          finish = start - 1 + &
               ngwf_basis%spheres(count)%n_ppds_sphere*cell%n_pts

          ! Find start and finish in src array
          start_src = ngwf_basis_src1%spheres(count_src1)%offset
          finish_src = start_src - 1 + &
               ngwf_basis_src1%spheres(count_src1)%n_ppds_sphere*cell%n_pts

          ! Copy NGWF from src1 to dest array
          if (loc_cmplx) then
             ngwfs_on_grid%z(start:finish) = ngwfs_on_grid_src1%z(start_src:finish_src)
          else
             ngwfs_on_grid%d(start:finish) = ngwfs_on_grid_src1%d(start_src:finish_src)
          end if
       end do

       do ifunc_src2=1,ngwf_basis_src2%num_on_atom(iat)
          count_src2 = count_src2 + 1
          count = count + 1

          ! Find start and finish in dest array
          start = ngwf_basis%spheres(count)%offset
          finish = start - 1 + &
               ngwf_basis%spheres(count)%n_ppds_sphere*cell%n_pts

          ! Find start and finish in src array
          start_src = ngwf_basis_src2%spheres(count_src2)%offset
          finish_src = start_src - 1 + &
               ngwf_basis_src2%spheres(count_src2)%n_ppds_sphere*cell%n_pts

          ! Copy NGWF from src2 to dest array
          if (loc_cmplx) then
             ngwfs_on_grid%z(start:finish) = ngwfs_on_grid_src2%z(start_src:finish_src)
          else
             ngwfs_on_grid%d(start:finish) = ngwfs_on_grid_src2%d(start_src:finish_src)
          end if

       end do

    end do

    if (count/=ngwf_basis%proc_num) then
       call utils_abort('Error in ngwfs_merge_sets: Wrong total number of &
            &NGWFs copied into merged set')
    end if

  end subroutine ngwfs_merge_sets


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_initialise_from_fragments(ngwf_basis, rep, mdl, hfxstate)

    !=================================================================!
    ! This subroutine initialises the EDA frozen state NGWFs from the !
    ! NGWFs of the fragments (spheres and ngwfs_on_grid),             !
    ! as well as sets up the supermolecule NGWF pivot table.          !
    !-----------------------------------------------------------------!
    ! Written by Max Phipps in 2014.                                  !
    ! Modified to remove pub_par by Joseph Prentice, August 2018      !
    ! Modified further for embedding by Joseph Prentice, September 2018!
    !=================================================================!


    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: HFX_STATE, hf_exchange_dkn_indep_stage
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, ngwf_rep_register_change
    use comms, only: pub_my_proc_id, comms_barrier
    use rundat, only: pub_frag_iatm, pub_use_hfx, pub_pol_emb_qmstar
    use fragment_data, only: NGWF_INDEX_MAP_TYPE, pub_frag_data, &
         pub_super2fragid_on_proc, fragment_data_get_ngwf_index_map, &
         fragment_data_set_ngwf_index_map
    use polarisable_embedding, only: polarisable_embedding_expand_ngwf_pairs
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout)        :: ngwf_basis(1)
    type(MODEL), intent(inout), target     :: mdl
    type(NGWF_REP), intent(inout)          :: rep
    type(HFX_STATE), intent(inout), target :: hfxstate

    ! Local Variables
    ! mjsp: EDA
    integer :: it, ingwf_frag, ingwf_super
    integer :: ingwf_frag_onproc ! the index of the fragment NGWF on its proc
    integer :: ingwf_super_onproc ! the index of the supermolecule NGWF on its proc
    integer :: ingwf, ngwfs_on_grid_offset, ngwf_count
    integer :: frag_proc, super_proc
    integer :: frag_index
    integer :: num_ngwfs

    type(NGWF_INDEX_MAP_TYPE) :: ngwf_index_map
    integer, allocatable :: atom_ids_by_cum_frag_ngwf(:)
    integer, allocatable :: atom_ids_by_super_ngwf(:)

    integer :: ierr
    character(255) :: errmsg = ""

    ! mjsp: Load/overwrite the set of NGWFs with the isolated fragment NGWFs:
    !
    ! The supermolecule NGWFs arrays are not neccesarily equivalent to
    ! simply concatenating the sets of fragment NGWFs.
    ! To handle this:
    ! - The indexes of the fragment NGWFs are referenced back to the
    !   atoms they belong to in the input file,
    ! - these are then referenced back to the NGWFs in the
    !   supermolecule ngwf_basis order (pub_frag2super_ngwf_idx).
    ! This module's functionality relies on the pivot table,
    ! pub_frag2super_ngwf_idx, which has been constructed using the
    ! parallelisation strategy data. The fragment NGWFs are filled
    ! into the supermolecule NGWFs storage in the correct order using
    ! the pub_frag2super_ngwf_idx array.

    num_ngwfs = 0

    ! Determine the number of fragment NGWFs in the supermolecule.
    do frag_index = 1, pub_frag_iatm(0)
       num_ngwfs = num_ngwfs + pub_frag_data(frag_index)%ngwf_basis(1)%num
    end do

    allocate( &
         atom_ids_by_cum_frag_ngwf(num_ngwfs), &
         stat = ierr, &
         errmsg = errmsg)

    call utils_alloc_check( &
         "ngwfs_initialise_from_fragments", &
         trim(errmsg), &
         ierr)

    ! Index of cumulative fragment NGWF.
    ingwf = 1

    ! Construct the map of cumulative fragment NGWF index to atom identifier.
    do frag_index = 1, pub_frag_iatm(0)
       do ingwf_frag = 1, pub_frag_data(frag_index)%ngwf_basis(1)%num
          atom_ids_by_cum_frag_ngwf(ingwf) = &
               pub_frag_data(frag_index)%basis_ref(ingwf_frag)

          ingwf = ingwf + 1
       end do
    end do

    call fragment_data_set_ngwf_index_map( &
         atom_ids_by_cum_frag_ngwf, &
         mdl, &
         ngwf_basis(1))

    deallocate( &
         atom_ids_by_cum_frag_ngwf, &
         stat = ierr, &
         errmsg = errmsg)

    call utils_dealloc_check( &
         "ngwfs_initialise_from_fragments", &
         trim(errmsg), &
         ierr)

    ngwf_index_map = fragment_data_get_ngwf_index_map()

    ! Used for filling ngwfs_on_grid:
    ngwfs_on_grid_offset = 1 ! index offset for ngwfs_on_grid  (offset
           ! by the previous fragment(s) number of entries in the array)

    ingwf = 1 ! NGWF counter for constructing the pivot table for
              ! pub_frag2super_ngwf_idx.

    ! mjsp: Loop the spheres of the fragments and fill the
    ! supermolecule spheres with the fragment spheres using the
    ! referencing scheme of pub_frag_data(0)%basis_ref and the fragments' basis_ref
    do it=1,pub_frag_iatm(0) ! loop fragments

       do ingwf_frag=1,pub_frag_data(it)%ngwf_basis(1)%num  ! loop NGWFs on fragment

          ! mjsp: obtain our ingwf in the current supermolecular
          ! NGWF distribution:
          ingwf_super = ngwf_index_map%supers_by_cum_frag(ingwf)


          ! mjsp: we are iterating through the NGWF indexes to reconstruct the
          ! frozen density state NGWFs.  However, our NGWF sphere may not be on
          ! our proc any longer - we may need to pull the data from another proc.
          ! --> check if the sphere is on another proc:

          ! mjsp: identify (1) the proc this fragment NGWF's data is on,
          ! and (2) which proc this data should be on in our current
          ! supermolecule NGWF distribution:
          super_proc = ngwf_basis(1)%proc_of_func(ingwf_super)
          frag_proc = pub_frag_data(it)%ngwf_basis(1)%proc_of_func(ingwf_frag)


          ! mjsp: get the array index of these NGWFs after having been
          ! distributed to their procs:
          ingwf_frag_onproc = ingwf_frag - pub_frag_data(it)%ngwf_basis(1)%first_on_proc(frag_proc) + 1
          ingwf_super_onproc = ingwf_super - ngwf_basis(1)%first_on_proc(super_proc) + 1

          ! mjsp: is the fragment sphere ingwf_frag on this proc?
          if ( pub_my_proc_id == frag_proc ) then

             ! mjsp: does this fragment data want to be moved to another proc?
             if ( pub_my_proc_id .NE. super_proc ) then

                ! mjsp: we should send this fragment sphere and ngwfs_on_grid
                ! subarray to another proc!
                call internal_send_ngwf()

             else

                ! mjsp: the fragment data is on our proc already!
                call internal_reorder_ngwf()

             end if

             ! mjsp: does this proc want to receive fragment NGWF data?
          else if ( pub_my_proc_id == super_proc ) then

             if ( pub_my_proc_id .NE. frag_proc ) then

                ! mjsp: the fragment data is on another proc!
                ! so receive this fragment sphere and ngwfs_on_grid subarray from the other proc
                call internal_recv_ngwf()

             end if
             ! else sphere was already on local to our proc: captured in above
             ! nested if statement block (pub_my_proc_id == frag_proc)

          end if


          ! mjsp: populate pub_super2fragid_on_proc referencing array (proc NGWF to
          ! fragment ID):
          if ( pub_my_proc_id == super_proc ) &
               pub_super2fragid_on_proc(ingwf_super_onproc) = it

          ! mjsp: synchronize to prevent the case where
          ! two or more procs broadcast to one proc
          ! simultaneously - in such a case data mismatches can occur.
          call comms_barrier

          ingwf = ingwf + 1  ! iterate NGWF counter

       end do ! NGWF

    end do ! fragment

    ! jd: *** NGWFS changed ***
    call ngwf_rep_register_change(rep,'ngwfs_initialise_from_fragments')
    ! jd: Re-expand new NGWFs in SWs, if needed
    if(pub_pol_emb_qmstar) then
       call polarisable_embedding_expand_ngwf_pairs(rep, ngwf_basis(1), mdl)
    end if
    if(pub_use_hfx) then
       call hf_exchange_dkn_indep_stage(hfxstate, mdl, 1, &
            mdl%regions(1)%par, rep, ngwf_basis(1))
    end if

    ! mjsp: synchronize
    call comms_barrier

  contains

    subroutine internal_send_ngwf()

        ! mjsp: This subroutine sends the fragment sphere and ngwf_on_grid of
        ! ingwf_frag_onproc on this proc to the proc super_proc.
        use comms, only: comms_send, comms_wait

        implicit none

        ! Local Variables
        integer :: arr_len  ! Variables used for replacing NGWFs
                            ! from isolated monomer calculations
                            ! in frozen density representation
        integer :: frag_el_offset
        integer :: send_handles(1:9)
        integer :: send_handle_id


        ! ================================ sphere

        ! send our temporary buffer sphere (to be received in this
        ! buffer sphere):
        call comms_send(super_proc, &
          pub_frag_data(it)%ngwf_basis(1)%spheres(ingwf_frag_onproc)%centre%x, &
          return_handle=send_handles(1), add_to_stack=.true.)

        call comms_send(super_proc, &
          pub_frag_data(it)%ngwf_basis(1)%spheres(ingwf_frag_onproc)%centre%y, &
          return_handle=send_handles(2), add_to_stack=.true.)

        call comms_send(super_proc, &
          pub_frag_data(it)%ngwf_basis(1)%spheres(ingwf_frag_onproc)%centre%z, &
          return_handle=send_handles(3), add_to_stack=.true.)

        ! mjsp: sphere%radius
        call comms_send(super_proc, &
          pub_frag_data(it)%ngwf_basis(1)%spheres(ingwf_frag_onproc)%radius, &
          return_handle=send_handles(4), add_to_stack=.true.)

        ! mjsp: sphere%n_ppds_sphere
        call comms_send(super_proc, &
          pub_frag_data(it)%ngwf_basis(1)%spheres(ingwf_frag_onproc)%n_ppds_sphere, &
          return_handle=send_handles(5), add_to_stack=.true.)

        ! mjsp: sphere%extended
        call comms_send(super_proc, &
          pub_frag_data(it)%ngwf_basis(1)%spheres(ingwf_frag_onproc)%extended, &
          return_handle=send_handles(6), add_to_stack=.true.)

        ! mjsp: wait for the sends to finish
        do send_handle_id=1,6
           call comms_wait(send_handles(send_handle_id), .true.)
        end do

        ! mjsp: sphere%ppd_list
        ! mjsp: send our sphere's ppd_list array
        call comms_send(super_proc, &
          pub_frag_data(it)%ngwf_basis(1)%spheres(ingwf_frag_onproc)%ppd_list(:,:), &
          return_handle=send_handles(7), add_to_stack=.false.)
        call comms_wait(send_handles(7))

        ! mjsp: sphere%offset
!        call comms_send(super_proc, &
!          pub_frag_data(it)%ngwf_basis%spheres(ingwf_frag_onproc)%offset)



        ! ================================ ngwfs_on_grid

        ! mjsp: also send this NGWF on the grid:

        ! mjsp: the array index that the fragment NGWF begins at:
        frag_el_offset = pub_frag_data(it)%ngwf_basis(1)%spheres(ingwf_frag_onproc)%offset

        ! mjsp: prepare the neccesary metadata:
        ! mjsp: number of elements in ngwfs_on_grid that we are filling:
        arr_len = pub_frag_data(it)%ngwf_basis(1)%spheres(ingwf_frag_onproc)%n_ppds_sphere * mdl%cell%n_pts
        call comms_send(super_proc, arr_len, &
             return_handle=send_handles(8), add_to_stack=.false.)
        call comms_wait(send_handles(8))

        ! mjsp: send the actual ngwfs_on_grid:
        ! NOTE: assume real only!
        call comms_send(super_proc, &
             pub_frag_data(it)%rep%ngwfs_on_grid(1)%d((frag_el_offset):(frag_el_offset+arr_len-1)), &
             return_handle=send_handles(9), add_to_stack=.false.)
        call comms_wait(send_handles(9))


    end subroutine internal_send_ngwf

    subroutine internal_recv_ngwf()

        ! mjsp: This subroutine receives the fragment sphere
        ! from frag_proc into tempsphere, checks tempsphere
        ! is the expected sphere to be received, and copies tempsphere into
        ! ingwf_super_onproc.  This also handles receiving the
        ! ngwfs_on_grid subarray for this fragment.

        use basis, only: SPHERE, basis_copy_sphere
        use comms, only: comms_irecv, comms_wait
        use utils, only: utils_abort, utils_dealloc_check

        implicit none

        ! Local Variables
        integer :: arr_len  ! Variables used for replacing NGWFs
                            ! from isolated monomer calculations
                            ! in frozen density representation
        type(SPHERE) :: tempsphere
        integer :: recv_handles(1:9)
        integer :: recv_handle_id


        ! ================================ sphere

        ! mjsp: the other proc to receive a sphere from has sent their
        ! sphere's data to tempsphere.
        ! receive this:

        ! mjsp:sphere%centre
        call comms_irecv(frag_proc, tempsphere%centre%x, &
             handle=recv_handles(1))
        call comms_irecv(frag_proc, tempsphere%centre%y, &
             handle=recv_handles(2))
        call comms_irecv(frag_proc, tempsphere%centre%z, &
             handle=recv_handles(3))

        ! mjsp: sphere%radius
        call comms_irecv(frag_proc, tempsphere%radius, &
             handle=recv_handles(4))

        ! mjsp: sphere%n_ppds_sphere
        call comms_irecv(frag_proc, tempsphere%n_ppds_sphere, &
             handle=recv_handles(5))

        ! mjsp: sphere%extended
        call comms_irecv(frag_proc, tempsphere%extended, &
             handle=recv_handles(6))

        ! mjsp: wait for the receives to finish
        do recv_handle_id=1,6
           call comms_wait(recv_handles(recv_handle_id))
        end do

        ! mjsp: allocate using the number of ppds we have received:
        allocate(tempsphere%ppd_list(2,tempsphere%n_ppds_sphere), stat=ierr)
        call utils_alloc_check('ngwfs_initialise_from_fragments','tempsphere%ppd_list',ierr)

        ! mjsp: we are now ready to fill our sphere's array from the other
        ! proc
        call comms_irecv(frag_proc, tempsphere%ppd_list(:,:), &
             handle=recv_handles(7))

        ! mjsp: wait for the receive to finish
        call comms_wait(recv_handles(7))

!        ! mjsp: sphere%offset
!        call comms_recv(frag_proc, tempsphere%offset)


        ! mjsp: Check proc data transfer success.
        ! This code should be rigorous and should never need to abort,
        ! but included in case of small chance of modifications causing
        ! the code to fail in the future.
        if (.not.((ngwf_basis(1)%spheres(ingwf_super_onproc)%centre%x .eq. tempsphere%centre%x).and. &
                 (ngwf_basis(1)%spheres(ingwf_super_onproc)%centre%y .eq. tempsphere%centre%y).and. &
                 (ngwf_basis(1)%spheres(ingwf_super_onproc)%centre%z .eq. tempsphere%centre%z)) ) &
             call utils_abort('Error in ngwfs_initialise_from_fragments: data mismatch during &
                  &NGWF data transfer between procs.')

        ! mjsp: now we have received the sphere and have checked it is
        ! correct, copy the temporary sphere into our proc's supermolecule
        ! sphere:
        call basis_copy_sphere(ngwf_basis(1)%spheres(ingwf_super_onproc),&
             tempsphere,ngwfs_on_grid_offset)

        ! mjsp: cleanup: deallocate our tempsphere's ppd_list array
        deallocate(tempsphere%ppd_list, stat=ierr)
        call utils_dealloc_check('ngwfs_initialise_from_fragments','tempsphere%ppd_list',ierr)

        ! ================================ ngwfs_on_grid

        ! mjsp: receive the neccesary metadata:
        call comms_irecv(frag_proc, arr_len, &
             handle=recv_handles(8))
        call comms_wait(recv_handles(8))

        ! mjsp: receive the actual NGWF on grid:
        ! NOTE: assume real only!
        call comms_irecv(frag_proc, &
             rep%ngwfs_on_grid(1)%d((ngwfs_on_grid_offset):(ngwfs_on_grid_offset+arr_len-1)), &
             handle=recv_handles(9))
        call comms_wait(recv_handles(9))

        ! mjsp: The offset for the NGWFs already filled into the
        !       supermolecule ngwfs_on_grid
        ngwfs_on_grid_offset = ngwfs_on_grid_offset + arr_len

    end subroutine internal_recv_ngwf

    subroutine internal_reorder_ngwf()

        ! mjsp: This subroutine handles copying the sphere data
        ! and the ngwfs_on_grid subarray from
        ! pub_frag_data to ngwf_basis and rep on this proc.
        ! This is called for the case when the data is already on this proc.

        use basis, only: basis_copy_sphere

        implicit none

        ! Local Variables
        integer :: arr_len  ! Variables used for replacing NGWFs
                            ! from isolated monomer calculations
                            ! in frozen density representation
        integer :: frag_el_offset


        ! ================================ sphere

        ! mjsp: from function_basis_copy_spheres:
        ! copy the sphere from one basis index to its correct index
        call basis_copy_sphere(ngwf_basis(1)%spheres(ingwf_super_onproc),&
             pub_frag_data(it)%ngwf_basis(1)%spheres(ingwf_frag_onproc),ngwfs_on_grid_offset)


        ! ================================ ngwfs_on_grid

        ! the array index that the fragment NGWF begins at:
        frag_el_offset = pub_frag_data(it)%ngwf_basis(1)%spheres(ingwf_frag_onproc)%offset

        ! mjsp: number of elements in ngwfs_on_grid that we are filling:
        arr_len = ngwf_basis(1)%spheres(ingwf_super_onproc)%n_ppds_sphere * mdl%cell%n_pts

        ! mjsp: Fill ngwfs_on_grid (from ingwf_frag of fragment to
        !       pub_frag2super_ngwf_idx(ingwf_frag) of supermolecule)
        ! NOTE: assume real only!
        rep%ngwfs_on_grid(1)%d((ngwfs_on_grid_offset):(ngwfs_on_grid_offset+arr_len-1)) = &
                    pub_frag_data(it)%rep%ngwfs_on_grid(1)%d((frag_el_offset):(frag_el_offset+arr_len-1))

        ! mjsp: The offset for the NGWFs already filled into the
        !       supermolecule ngwfs_on_grid
        ngwfs_on_grid_offset = ngwfs_on_grid_offset + arr_len

    end subroutine internal_reorder_ngwf

  end subroutine ngwfs_initialise_from_fragments


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function ngwfs_at_point(shell,mx,my,mz,radial_ngwf,periodic_centre, &
       current_point,ngwf_radius)

    !=================================================================!
    ! This function returns the value of an ngwf at a given point in  !
    ! real space given a radial representation of the ngwf.           !
    !-----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.                       !
    ! f-functions added by Chris-Kriton Skylaris on 9/5/2007.         !
    ! Modified by Nicholas Hine on 03/09/2010 to use RADIAL_NGWF_TYPE.!
    ! Extended to Cartesian NGWFs by Arihant Bhandari on 30/11/2021.  !
    !=================================================================!

    use geometry,  only: point, operator(-), geometry_distance
    use services, only: services_linear_interpolation, services_locate_interp
    use spherical_wave, only: sw_real_sph_harm
    use rundat, only: pub_dftb, pub_dftb_cartesian_ngwfs

    implicit none

    ! Arguments
    real(kind=DP) ::  ngwfs_at_point
    integer, intent(in) :: shell
    integer, intent(in) :: mx, my, mz ! ab: components of ang. mom. (l=mx+my+mz)
    type(RADIAL_NGWF_TYPE),intent(in) :: radial_ngwf
    real(kind=DP), intent(in) :: ngwf_radius
    type(POINT), intent(in) :: current_point
    type(POINT), intent(in) :: periodic_centre

    ! Local variables
    real(kind=DP) :: local_distance, radial_value, angular_value
    integer :: interp_index
    type(point) :: r_vec

    ngwfs_at_point = 0.0_DP ! qoh: initialise to prevent compiler warning

    ! cks: here centre can be the centre of a ngwf function inside the
    ! cks: simulation cell but it can also be the centre of a periodic
    ! cks: image of a ngwf function which happens to have values in a
    ! cks: ppd of the simulation cell.

    local_distance = geometry_distance(periodic_centre,current_point)
    if (local_distance <= ngwf_radius) then

       ! cks: find the array index of the lowest point for the intepolation,
       ! cks: provided the point does not fall out of the region of validity
       ! cks: of the fireball ngwf function.in the opposite case set
       ! cks: radial_value=0.0_DP
       if (local_distance<radial_ngwf%rad(radial_ngwf%npts)) then

          interp_index = services_locate_interp(local_distance, &
               radial_ngwf%rad,radial_ngwf%npts)

          radial_value = services_linear_interpolation(local_distance, &
               radial_ngwf%func_real(interp_index,shell), &
               radial_ngwf%func_real(interp_index+1,shell), &
               radial_ngwf%rad(interp_index), &
               radial_ngwf%rad(interp_index+1))
       else
          radial_value = 0.0_DP
       endif

       r_vec = current_point - periodic_centre
       ! ab: i.e. for Spherical NGWFs
       if (.not. (pub_dftb .and. pub_dftb_cartesian_ngwfs)) then
          angular_value = sw_real_sph_harm(r_vec%x,r_vec%y,r_vec%z, &
               local_distance,radial_ngwf%angmom(shell),mz)
       ! ab: for Cartesian NGWFs the angular part and radial part is calculated
       ! ab: following Ref. [1]:(1.19):
       ! ab: angular_value = (x-Ax)^mx . (y-Ay)^my . (z-Az)^mz
       else
          angular_value = 1.0_DP
          if (mx > 0) then
             angular_value = angular_value * r_vec%x**mx
          end if
          if (my > 0) then
             angular_value = angular_value * r_vec%y**my
          end if
          if (mz > 0) then
             angular_value = angular_value * r_vec%z**mz
          end if
       end if

       ngwfs_at_point = radial_value*angular_value

    endif

  end function ngwfs_at_point


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_normalise(ngwf_on_grid, &   ! INPUT-OUTPUT
       ngwf_basis, cell)                       ! INPUT

    !=========================================================!
    ! Normalise a set of NGWFs to unity.                      !
    !---------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 23/10/2001.         !
    ! Modified to use function basis type by Nicholas Hine on !
    ! 14/12/2009.                                             !
    !=========================================================!

    use constants, only: DP
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort

    implicit none

    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(inout) :: ngwf_on_grid( &
         ngwf_basis%size_on_grid)
    type(CELL_INFO), intent(in) :: cell

    ! cks: internal declarations
    integer :: row, start, finish
    real(kind=DP) :: norm_constant

    do row=1,ngwf_basis%num
       start=ngwf_basis%spheres(row)%offset
       finish=start+(ngwf_basis%spheres(row)%n_ppds_sphere)*cell%n_pts -1

       norm_constant=sum(ngwf_on_grid(start:finish)*ngwf_on_grid(start:finish) )
       norm_constant=sqrt(norm_constant*cell%weight)

       if (norm_constant > 0.0_DP) then
          ngwf_on_grid(start:finish)=ngwf_on_grid(start:finish)/norm_constant
       else
          call utils_abort('Error in ngwfs_normalise(): Invalid norm_constant C&
               & for function F. The values of F and C follow. ', &
               row, opt_real_to_print1=norm_constant)
       endif

    enddo

  end subroutine ngwfs_normalise


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function convolution_factor(periodic_centre, current_point, &
               ngwf_radius, conv_function, region_width)

    !=========================================================!
    ! Function to determine the convolution factor for the    !
    ! current point, according to the convolution function.   !
    !---------------------------------------------------------!
    ! Written by Andrea Greco on 20/01/2016.                  !
    !=========================================================!

    use geometry,  only: point, operator(-), geometry_distance
    use utils, only: utils_erfc

    implicit none

    ! output
    real(kind=DP) :: convolution_factor
    ! input
    type(POINT), intent(in) :: current_point
    type(POINT), intent(in) :: periodic_centre
    real(kind=DP), intent(in)     :: ngwf_radius
    character(len=80), intent(in) :: conv_function
    ! agreco: apply convolution only in a limited region
    ! at the edge of the NGWF sphere
    real(kind=DP), optional, intent(in) :: region_width

    ! Local variables
    real(kind=DP) :: local_distance
    real(kind=DP) :: loc_region_width

    loc_region_width = ngwf_radius

    ! agreco: overwrite only if a positive value
    ! for region_width is specified; consider
    ! to add a warning message when a negative
    ! value is specified, reverting to default
    if (present(region_width)) then
       if (region_width>0.0_DP) then
          loc_region_width = region_width
       end if
    end if

    ! agreco: partially extended NGWFs already taken into account
    ! by modifying periodic_centre, current_point accordingly in
    ! parent routine
    local_distance = geometry_distance(periodic_centre,current_point)

    ! agreco: currently only the use of erfc is supported
    ! only apply convolution in the required region
    if (conv_function=='erfc'.and. &
       local_distance>(ngwf_radius-loc_region_width)) then
       convolution_factor = utils_erfc((local_distance &
          - ngwf_radius + loc_region_width)/ngwf_radius)
    else
       convolution_factor = 1.0
    end if

  end function convolution_factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function random_normal(mean,sigma)

    !=========================================================!
    ! Function to generate random numbers with normal         !
    ! distribution.
    !---------------------------------------------------------!
    ! Written by Andrea Greco on 18/01/2016.                  !
    ! Adapted from http://www.sdsc.edu/~tkaiser/f90.html      !
    !=========================================================!

    implicit none

    ! output
    real(kind=DP) :: random_normal
    ! input
    real(kind=DP), intent(in) :: mean
    real(kind=DP), intent(in) :: sigma
    ! local variables
    real(kind=DP) :: tmp
    integer, save :: flag = 0
    real(kind=DP) :: fac,rsq,r1,r2,num
    real(kind=DP), save :: gsave

    if (flag.eq.0) then
       rsq=2.0_DP
       do while(rsq.ge.1.0_DP.or.rsq.eq.0.0_DP)
          call random_number(num)
          r1=2.0_DP*num-1.0_DP
          call random_number(num)
          r2=2.0_DP*num-1.0_DP
          rsq=r1*r1+r2*r2
       end do
       fac=sqrt(-2.0_DP*log(rsq)/rsq)
       gsave=r1*fac
       tmp=r2*fac
       flag = 1
    else
       tmp=gsave
       flag=0
    end if

    random_normal=tmp*sigma+mean

    return

  end function random_normal


end module ngwfs


