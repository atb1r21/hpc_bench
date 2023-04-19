! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                  Electronic transport module                   !
!                                                                !
!----------------------------------------------------------------!
! Written by Simon M.-M. Dubois in April 2010                    !
! Revisions by Robert Bell, Nov 2013                             !
!                                                                !
! Note: parts of this module assume time-reversal symmetry       !
!================================================================!

module etrans

  use constants, only: DP
  use rundat, only: pub_debug_on_root

  implicit none

  private

  !#ifdef MKLOMP
  !#include "mkl.h"
  !#endif

  public :: etrans_calculate
  public :: etrans_init_leads

  ! Type definition

  ! Container for lead/LCR information and definitions
  type, public :: block_def
     logical            :: forward
     character(len=4)   :: type    ! lead block, or lcr block
     logical            :: ioread
     character(len=80)  :: iofile
     integer            :: num_unit_cells
     integer            :: atms(4) ! Range of atoms (original order)
     integer            :: orbs(4) ! Range of orbitals (original order)
     real(kind=dp), allocatable  :: efermi(:) ! the (lead) Fermi energy
     logical                     :: have_lead_mu
     real(kind=DP), allocatable  :: occ(:)     ! electronic occupancy of this lead
     real(kind=DP)               :: ionic_popn ! ionic charge of this lead
  end type block_def

  ! Container for data arrays for leads/LCR
  type :: hsm_type
     integer                       :: norb
     real(kind=DP)                 :: eref      ! the reference energy for the Hamiltonian
     character(len=20)             :: ham_type  ! valence/joint hamiltonian

     ! arrays for the LCR device
     real(kind=DP), allocatable    :: h_lcr(:,:)
     real(kind=DP), allocatable    :: s_lcr(:,:)

     ! arrays for the lead
     real(kind=DP), allocatable    :: s00(:,:), s01(:,:)
     real(kind=DP), allocatable    :: h00(:,:,:), h01(:,:,:)
     complex(kind=DP), allocatable :: self(:,:), self_rev(:,:)
  end type hsm_type

  ! container for transmission matrices
  type :: tmat_type
     integer                       :: norb
     complex(kind=DP), allocatable :: T(:,:) ! transmission matrix
  end type tmat_type

  ! container defining the contour
  type, public :: epath
     character(len=2)   :: type      ! TR = transmission/dos
     integer            :: nep       ! total number of points
     integer            :: seg(4)    ! (nreal, ncircle, nhline, nvline)

     real(kind=DP)      :: eref      ! system reference energy

     ! contour distribution information
     integer, allocatable          :: nepproc(:)
     integer, allocatable          :: nepptr(:)

     ! energy points and weights (not implemented)
     complex(kind=dp), allocatable :: ep(:)
     complex(kind=dp), allocatable :: ew(:)
  end type epath

  integer, public                      :: pub_nleads
  type(block_def), public, allocatable :: pub_leads(:)
  type(block_def), public              :: pub_lcr
  real(kind=DP), public, allocatable   :: pub_etrans_eigchan_en(:)
  integer, public                      :: pub_etrans_plot_num

contains

  !====================================================================!
  !====================================================================!

  subroutine etrans_calculate(ham,ham_type,rep,denskern,ngwf_basis,mdl)

    !================================================================!
    ! Main routine for calculating electronic transport              !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in April 2010                    !
    !================================================================!

    use constants, only: dp, pi, stdout, NORMAL
    use comms, only: pub_my_proc_id, pub_on_root, comms_barrier, &
         comms_reduce
    use function_basis, only: FUNC_BASIS
    use greenf, only: SPARGF, greenf_init_blocking, greenf_init_matrices, &
         greenf_prepare_inverse, trimat_deposit_atom_block, &
         greenf_alloc_blocks, greenf_destroy, greenf_dos
    use linalg, only: linalg_invert_serial
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_HAM, NGWF_REP
    use rundat, only: pub_rootname, pub_task, pub_etrans_lcr, pub_etrans_bulk, &
         pub_print_qc, pub_cond_calculate, pub_etrans_write_hs, pub_output_detail, &
         pub_etrans_num_eigchan, pub_etrans_plot_eigchan, pub_etrans_write_xyz,&
         pub_etrans_seed_lead, pub_num_spins
#ifdef MKLOMP
    use rundat, only: pub_threads_num_mkl
#endif
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_destroy
    use sparse_embed, only: SPAM3_EMBED_ARRAY, &
         sparse_embed_extract_from_array, sparse_embed_destroy_extracted_array
    use timer, only: timer_clock, timer_check_iteration_time, wrappers_etime
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(in)   :: ham
    type(NGWF_REP), intent(in)   :: rep
    type(SPAM3_EMBED_ARRAY), intent(in):: denskern
    ! rc2013: could include multiple bases
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    character(len=*), intent(in) :: ham_type
    type(MODEL), intent(in)      :: mdl

    ! Local variables
    type(epath)       :: econt
    complex(kind=DP)  :: energy
    type(hsm_type), allocatable  :: lead_hsm(:)
    type(SPARGF)      :: greenf_obj

    complex(kind=dp), allocatable :: greenf_lead(:,:)
    real(kind=DP), allocatable    :: dos_lcr(:,:,:), trc_lcr(:,:,:)
    real(kind=DP), allocatable    :: dos_lead(:,:,:), trc_lead(:,:,:)
    real(kind=DP), allocatable    :: trc_lcr_eigchan(:,:,:,:)

    integer, allocatable  :: gferr_lcr(:,:,:), gferr_lead(:,:,:)
    integer  :: norb, nleads
    integer  :: nep, ntrc
    integer  :: il, jl
    integer  :: ieloc, ienergy, is, iorb, jorb
    integer  :: iat, jat
    integer  :: ierr, info
    character(len=80)   :: filename
    logical  :: out_of_runtime
    logical  :: compute_eigchans
    real(kind=DP) :: iter_time, last_report_time
    logical  :: first_time_report

    integer :: lead_blocks(2,pub_nleads)
    ! rc2013: structure to avoid array temporaries
    type(SPAM3) :: ham_array(pub_num_spins)

    call timer_clock('etrans_calculate',1)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering etrans_calculate'

    if (pub_on_root) then
       write(stdout,'(/a)') repeat('~',80)
       write(stdout,'(30x,a)') "Electronic Transport"
       write(stdout,'(a)') repeat('~',80)
    endif


#ifdef MKLOMP
    if (pub_on_root .and. pub_output_detail >= NORMAL .and. &
      pub_threads_num_mkl .ne. 1) then
       write(stdout,'(a,i2,a)') &
           'Using MKL threading: ', pub_threads_num_mkl, ' threads'
    endif
    call mkl_set_num_threads(pub_threads_num_mkl)
#endif

    ! Initialise device setup
    if (pub_task .ne. 'NEGF') then
       call etrans_init_setup(ngwf_basis)
    endif


    ! init sparse selected inverse blocking scheme
    if (pub_etrans_lcr) then

       ! select the lead to seed the tri-diagonal partitioning
       lead_blocks(:,1) = pub_leads(pub_etrans_seed_lead)%atms(1:2)
       jl = 2
       do il = 1, pub_nleads
          if (il == pub_etrans_seed_lead) cycle
          lead_blocks(:,jl) = pub_leads(il)%atms(1:2)
          jl = jl + 1
       enddo

       ! rc2013: to get around array temporaries, copy the Hamiltonian
       ! sub-matrix into a new structure
       call sparse_embed_extract_from_array(ham_array,ham%ham)

       call greenf_init_blocking(greenf_obj,ham_array,ham_type,&
                          groups_interact=.false.,force_groups=lead_blocks,&
                          first_atom=pub_lcr%atms(1),last_atom=pub_lcr%atms(2))

       ! rc2013: now destroy this abomination
       call sparse_embed_destroy_extracted_array(ham_array)
    endif

    ! make estimate of memory required
    call internal_memory_estimate

    ! print warning if this is a conduction calculation
    if(ham_type=='joint'.and.pub_on_root) then
       write(stdout,'(5(/a)//)') &
         ' WARNING: Transport calculations using the joint (valence+conduction)', &
         '          basis may be incorrect due to numerical noise when calculating', &
         '          the lead self energy/band structure.', &
         '          Check the results below the Fermi energy agree with the', &
         '          transmission calculating using the valence NGWFs only.'
    endif

    ! write device/lead co-ordinates, if requested
    if (pub_etrans_write_xyz.and.pub_on_root) call etrans_write_xyz(mdl%elements)

    ! Initialise hsm matrices
    call etrans_init_hsm(lead_hsm,ham_type)
    if (pub_etrans_lcr) call greenf_alloc_blocks(greenf_obj)

    ! gather Lead hsm matrices
    ! LCR matrices not yet required
    call etrans_gather_lead_hsm(lead_hsm,ham,rep,ngwf_basis)

    ! Compute lead/LCR reference energies
    call etrans_compute_erefs(lead_hsm,ham,ham_type,rep,denskern,&
        & mdl%elements,mdl%cell,ngwf_basis,econt%eref)

    ! write HSM to disk, if required
    if (pub_etrans_write_hs .and. pub_on_root) then
       ! leads
       write(stdout,'(/a)',advance='no') " Writing lead H/S matrices ..."
       do il= 1, pub_nleads
          write(filename,'(2a,i2.2,a)') trim(pub_rootname), '_lead',il,'.hsm'
          call etrans_write_hsm(pub_leads(il),lead_hsm(il),trim(filename))
       enddo
       write(stdout,'(a)') " done"
    endif

    ! Initialise energy contour for transmission coefficient
    call etrans_init_energy_contour(econt,'trans')

    ! Initialise local variables
    nep = econt%nep
    nleads = pub_nleads
    ntrc = (pub_nleads*pub_nleads-pub_nleads)/2
    compute_eigchans = (pub_etrans_num_eigchan > 0)

    ! Allocate dos and transmission coefficient matrices
    if (pub_etrans_lcr) then
       allocate(trc_lcr(ntrc,pub_num_spins,nep),stat=ierr)
       call utils_alloc_check('etrans_calculate','trc_lcr',ierr)
       allocate(dos_lcr(1,pub_num_spins,nep),stat=ierr)
       call utils_alloc_check('etrans_calculate','dos_lcr',ierr)
       allocate(gferr_lcr(1,pub_num_spins,nep),stat=ierr)
       call utils_alloc_check('etrans_calculate','gferr_lcr',ierr)
       if (compute_eigchans) then
          allocate(trc_lcr_eigchan(pub_etrans_num_eigchan,pub_nleads,pub_num_spins,nep),stat=ierr)
          call utils_alloc_check('etrans_calculate','trc_lcr_eigchan',ierr)
       endif
       ! initialise all data to zero so comms can be done with a sum
       trc_lcr = 0.0_dp
       dos_lcr = 0.0_dp
       gferr_lcr = 0
       if (compute_eigchans) then
          trc_lcr_eigchan = 0.0_DP
       endif
    endif

    if (pub_etrans_bulk) then
       allocate(trc_lead(nleads,pub_num_spins,nep),stat=ierr)
       call utils_alloc_check('etrans_calculate','trc_lead',ierr)
       allocate(dos_lead(nleads,pub_num_spins,nep),stat=ierr)
       call utils_alloc_check('etrans_calculate','dos_lead',ierr)
       allocate(gferr_lead(nleads,pub_num_spins,nep),stat=ierr)
       call utils_alloc_check('etrans_calculate','gferr_lead',ierr)
       trc_lead = 0.0_dp
       dos_lead = 0.0_dp
       gferr_lead = 0
    endif


    !=================================================================
    ! Main loop over the spins and energy
    !=================================================================

    spin_loop: do is = 1, pub_num_spins

       ! Write header
       if (pub_on_root .and. pub_output_detail>=NORMAL) then

          write(stdout,'(/,a)',advance='no') &
               '=== Calculating DOS and transmission coefficients '
          if (pub_num_spins==2 .and. is == 1) then
             write(stdout,'(a,i2)') 'for spin UP'
          elseif (pub_num_spins==2 .and. is == 2) then
             write(stdout,'(a,i2)') 'for spin DN'
          else
             write(stdout,'(/)')
          endif
       endif

       ! time routine
       last_report_time = wrappers_etime()
       first_time_report = .true.


       ! gather matrices for this spin channel
       if (pub_etrans_lcr) then
          ! rc2013: to get around array temporaries, copy the Hamiltonian
          ! sub-matrix into a new structure
          call sparse_create(ham_array(is), ham%ham(is)%p)
          call sparse_copy(ham_array(is), ham%ham(is)%p)

          call greenf_init_matrices(greenf_obj,ham_array, &
               rep%overlap%p,ispin=is)
          call internal_remove_short_circuits(greenf_obj)

          ! rc2013: now destroy this abomination
          call sparse_destroy(ham_array(is))
       endif

       ! write HSM to disk, if required ! TODO
       !if (pub_etrans_write_hs .and. pub_on_root) then
       !   ! LCR
       !   write(stdout,'(a)',advance='no') " Writing LCR  H/S matrices ..."
       !   if (pub_num_spins==2 .and. is == 1) then
       !     write(filename,'(2a)') trim(pub_rootname), '_lcr_UP.hsm'
       !   elseif (pub_num_spins==2 .and. is == 2) then
       !     write(filename,'(2a)') trim(pub_rootname), '_lcr_DN.hsm'
       !   else
       !     write(filename,'(2a)') trim(pub_rootname), '_lcr.hsm'
       !   endif
       !   call etrans_write_hsm(pub_lcr,lcr_hsm,trim(filename))
       !   write(stdout,'(a/)') " done"
       !endif

       ! loop over energies
       do ieloc = 1, econt%nepproc(pub_my_proc_id)

          call timer_check_iteration_time('etrans_energy','start')

          ienergy = econt%nepptr(pub_my_proc_id) + ieloc
          energy = econt%ep(ienergy)

          call etrans_compute_self(lead_hsm,energy,is,info)

          if(info/=0) then
             gferr_lead(:,is,ienergy) = info
             gferr_lcr(1,is,ienergy) = info
             call timer_check_iteration_time('etrans_energy','cancel')
             cycle
          endif

          !-----------------------------
          ! LCR transmission -----------
          !-----------------------------
          if (pub_etrans_lcr) then

             ! deposit self energy
             do il = 1, pub_nleads
                ! add H to self energy
                lead_hsm(il)%self(:,:) = lead_hsm(il)%self(:,:) + &
                    cmplx(lead_hsm(il)%h00(:,:,is),kind=DP)
                call trimat_deposit_atom_block(lead_hsm(il)%self(:,:), &
                        greenf_obj%h(1), greenf_obj, &
                        pub_leads(il)%atms(1),pub_leads(il)%atms(2), &
                        pub_leads(il)%atms(1),pub_leads(il)%atms(2))
                ! remove H from self energy
                lead_hsm(il)%self(:,:) = lead_hsm(il)%self(:,:) - &
                    lead_hsm(il)%h00(:,:,is)
             enddo

             ! prepare sparse inverse for this energy
             call greenf_prepare_inverse(greenf_obj,energy,ierr)


             ! compute DOS
             if (ierr .eq. 0) then
                call greenf_dos(dos_lcr(1,is,ienergy),greenf_obj,ierr)
             endif

             ! compute transmission coefficients
             if (ierr .eq. 0) then
                if (compute_eigchans) then
                   call etrans_compute_transmission_channels(greenf_obj,lead_hsm, &
                       trc_lcr(:,is,ienergy), ierr, compute_eigchans, &
                       trc_lcr_eigchan(:,:,is,ienergy))
                else
                   call etrans_compute_transmission_channels(greenf_obj,lead_hsm, &
                       trc_lcr(:,is,ienergy), ierr)
                endif
             endif

             ! handle any errors in the inverse
             if (ierr .ne. 0) then
                dos_lcr(1,is,ienergy) = -huge(1.d0)
                trc_lcr(:,is,ienergy) = -huge(1.d0)
                gferr_lcr(1,is,ienergy) = ierr
             endif
          endif


          !-----------------------------
          ! Bulk transmission ----------
          !-----------------------------
          if (pub_etrans_bulk) then

             do il = 1, nleads

                ! Initialise leads_greenf^(-1)
                norb = lead_hsm(il)%norb
                allocate(greenf_lead(norb,norb),stat=ierr)
                call utils_alloc_check('etrans_calculate','greenf_lead',ierr)

                greenf_lead = energy*cmplx(lead_hsm(il)%s00(:,:),kind=dp) &
                     - cmplx(lead_hsm(il)%h00(:,:,is),kind=dp)
                greenf_lead = greenf_lead - lead_hsm(il)%self(:,:) - lead_hsm(il)%self_rev(:,:)

                ! Compute greenf
                call linalg_invert_serial(greenf_lead,ierr)

                if (ierr .eq. 0) then

                   ! Compute dos
                   do iorb = 1 , norb
                      do jorb = 1, norb
                         dos_lead(il,is,ienergy) = dos_lead(il,is,ienergy) - &
                              aimag(greenf_lead(iorb,jorb)*lead_hsm(il)%s00(jorb,iorb))
                      enddo
                   enddo
                   dos_lead(il,is,ienergy) = dos_lead(il,is,ienergy) / PI

                   ! Compute transmission coefficients
                   call etrans_compute_transmission(norb,lead_hsm(il)%self,norb,&
                        lead_hsm(il)%self_rev,greenf_lead,trc_lead(il,is,ienergy))

                   deallocate(greenf_lead,stat=ierr)
                   call utils_dealloc_check('etrans_calculate','greenf_lead',ierr)
                else
                   dos_lead(il,is,ienergy) = -huge(1.d0)
                   trc_lead(il,is,ienergy) = -huge(1.d0)
                   gferr_lead(il,is,ienergy) = ierr
                endif

             enddo
          endif


          if (pub_on_root) call internal_print_status

          ! Check if sufficient runtime remains
          call timer_check_iteration_time('etrans_energy','stop', &
               out_of_runtime,global_mpi_loop=.false.)

          ! if not, exit energy and spin loops
          if (out_of_runtime) exit spin_loop

       enddo


       if(pub_on_root) then
          write(stdout,'(a)')  ' + --------------------------------------------'
          write(stdout,'(a)')  ' + Root done.'
          write(stdout,'(a)')  ' + Waiting for other procs to complete...'
          write(stdout,'(a/)') ' + --------------------------------------------'
          write(stdout,'(a/)') '==='
       endif

       call comms_barrier

    enddo spin_loop

    ! check if any proc ran out of run time
    call comms_reduce('OR',out_of_runtime)
    if (out_of_runtime.and.pub_on_root) then
       write(stdout,'(/a)') &
          "  WARNING: At least one proc has detected insufficient run time to&
          & complete."
       write(stdout,'(a/)') &
          "           Attempting to report the results so far."
    endif
    !=================================================================

    ! destroy dense matrices in Green's function
    if (pub_etrans_lcr) then
       call greenf_destroy(greenf_obj)
    endif

    ! Communicate data
    if (pub_etrans_lcr) then
       call comms_reduce('SUM',trc_lcr)
       call comms_reduce('SUM',dos_lcr)
       call comms_reduce('SUM',gferr_lcr)
       if (compute_eigchans) then
          do il=1, nep
            call comms_reduce('SUM',trc_lcr_eigchan(:,:,:,il))
          enddo
       endif
    endif
    if (pub_etrans_bulk) then
       call comms_reduce('SUM',trc_lead)
       call comms_reduce('SUM',dos_lead)
       call comms_reduce('SUM',gferr_lead)
    endif

    ! Write transmission coefficients to file
    if (pub_on_root) then
       write(stdout,'(/a)') '=== Writing results to file'
       if (pub_etrans_lcr) then
          call etrans_write_trc('LCR',ham_type,econt%ep,trc_lcr,gferr_lcr)
          call etrans_write_dos('LCR',ham_type,econt%ep,dos_lcr,gferr_lcr)
          if (compute_eigchans) then
              call etrans_write_trc_eigchan('LCR_channels',ham_type,econt%ep, &
                   trc_lcr_eigchan,gferr_lcr)
          endif
       endif
       if (pub_etrans_bulk) then
          call etrans_write_trc('BULK',ham_type,econt%ep,trc_lead,gferr_lead)
          call etrans_write_dos('BULK',ham_type,econt%ep,dos_lead,gferr_lead)
       endif
    end if



    ! compute and plot eigenchannels
    if (pub_etrans_plot_eigchan .and. pub_etrans_plot_num>0) then
       ! select device subset of elements
       iat = pub_lcr%atms(1)
       jat = pub_lcr%atms(2)
       call etrans_compute_eigchan(ham,ham_type,rep,ngwf_basis,lead_hsm,&
            mdl)
    endif




    if (pub_print_qc .and. pub_on_root) &
         call etrans_print_qc(dos_lcr, trc_lcr, dos_lead, trc_lead, econt%eref)


    ! destroy storage of results
    if (pub_etrans_bulk) then
       deallocate(gferr_lead,stat=ierr)
       call utils_dealloc_check('etrans_calculate','gferr_lead',ierr)
       deallocate(dos_lead,stat=ierr)
       call utils_dealloc_check('etrans_calculate','dos_lead',ierr)
       deallocate(trc_lead,stat=ierr)
       call utils_dealloc_check('etrans_calculate','trc_lead',ierr)
    endif

    if (pub_etrans_lcr) then

       deallocate(gferr_lcr,stat=ierr)
       call utils_dealloc_check('etrans_calculate','gferr_lcr',ierr)
       deallocate(dos_lcr,stat=ierr)
       call utils_dealloc_check('etrans_calculate','dos_lcr',ierr)
       deallocate(trc_lcr,stat=ierr)
       call utils_dealloc_check('etrans_calculate','trc_lcr',ierr)
       if (allocated(trc_lcr_eigchan)) then
          deallocate(trc_lcr_eigchan,stat=ierr)
          call utils_dealloc_check('etrans_calculate','trc_lcr_eigchan',ierr)
       endif
    endif


    ! Reset energy contour
    call etrans_free_energy_contour(econt)
    ! do not deallocate eigenchannel energies if this is a valence run,
    ! and a joint run is to follow
    if (.not.(trim(ham_type) == 'valence' .and. pub_cond_calculate)) then
       call etrans_free_eigchan_contour
    endif

    ! Free space associated with lead hsm matrices
    call etrans_destroy_lead_hsm(lead_hsm)

    ! do not deallocate leads if this is a valence run, and a joint run is to follow
    if (.not.(trim(ham_type) == 'valence' .and. pub_cond_calculate)) then
       call etrans_destroy_leads
    endif

    if (pub_on_root) write(stdout,'(/a/)') repeat('~',80)

    call timer_clock('etrans_calculate',2)

#ifdef MKLOMP
    call mkl_set_num_threads(1)
#endif

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving etrans_calculate'

    return

    !========================================!
    !========================================!

    contains

    subroutine internal_memory_estimate

       !================================================================!
       ! Peak memory contributions                                      !
       !   LCR:        greenf_obj                                       !
       !  LEAD: (real) h00,h01,s00,s01; (complex) 2*self, GF, inv_GF    !
       ! other: (real) lcr_trc,lcr_dos, lead_trc,lead_dos               !
       !               lead_bs,
       !----------------------------------------------------------------!
       ! Written by Robert Bell in November 2013                        !
       !================================================================!

       use constants, only: real_size, cmplx_size, LONG
       use greenf, only: greenf_memory_estimate
       use utils, only: utils_report_memory_estimate

       implicit none

       ! internal variables
       integer(kind=LONG) :: lcr_mem       ! memory required for lcr data structures
       integer(kind=LONG) :: lead_mem      ! memory required for lead structures

       ! convert type sizes to longs to avoid integer overflow
       integer(kind=LONG), parameter :: cmplx_size_long=int(cmplx_size,kind=LONG)
       integer(kind=LONG), parameter :: real_size_long=int(real_size,kind=LONG)


       integer :: il, norb, max_norbs

       ! large LCR structures
       lcr_mem = 0_LONG
       if (pub_etrans_lcr) then
          lcr_mem = greenf_memory_estimate(greenf_obj)
       endif

       ! lead structures
       lead_mem = 0
       max_norbs = -1
       do il=1,pub_nleads
          norb = pub_leads(il)%orbs(2) - pub_leads(il)%orbs(1) + 1
          max_norbs = max(max_norbs,norb)
          lead_mem = lead_mem + 2*(pub_num_spins+1)*norb**2*real_size_long+& ! h/s etc
                     2*norb**2 * cmplx_size_long                     ! self
       enddo

       ! largest lead temporary matrices (max use in compute_self/transfer)
       lead_mem = lead_mem + max_norbs**2 * (14*cmplx_size_long)



       ! report
       if (pub_etrans_lcr) then
          call utils_report_memory_estimate('Electronic Transport', &
               (/'LCR arrays ','Lead arrays'/), (/lcr_mem,lead_mem/))
               ! this space ^ is required - entries must have the same length
       else
          call utils_report_memory_estimate('Electronic Transport', &
               (/'Lead arrays'/), (/lead_mem/))
       endif

    end subroutine internal_memory_estimate

    !========================================!
    !========================================!

    subroutine internal_print_status

       implicit none

       real(kind=DP), parameter :: report_time_interval = 300.0_DP ! 5 minutes


       iter_time = wrappers_etime() - last_report_time

       if (iter_time > report_time_interval.and.pub_on_root) then

          last_report_time = wrappers_etime()

          if (first_time_report) then
            first_time_report = .false.

            write(stdout,'(a)') ' + --------------------------------------------'
            write(stdout,'(a)') ' + Reporting from root proc:'
            write(stdout,'(a)') ' + --------------------------------------------'
          endif

          write(stdout,'(a,i8,a,i8,a)') ' + Completed ',ieloc, ' of ', &
               econt%nepproc(pub_my_proc_id), ' energy points'

       endif

    end subroutine internal_print_status

    !========================================!
    !========================================!

    subroutine internal_remove_short_circuits(greenf_obj)

       use constants, only: DP, cmplx_0
       use greenf, only: SPARGF, trimat_deposit_atom_block
       use parallel_strategy, only: par=>pub_par
       use utils, only: utils_alloc_check, utils_dealloc_check

       implicit none

       ! arguments
       type(SPARGF), intent(inout) :: greenf_obj

       ! internal
       complex(kind=DP), allocatable :: zero(:,:)
       integer :: ilead, col_at, row_at
       integer :: col_orb, row_orb
       integer :: col_block, row_block
       integer :: ierr, ispin
       integer :: max_ngwfs

       ! maximum number of ngwfs of these atoms
       max_ngwfs = maxval(greenf_obj%ngwfs_on_atom(:))
       allocate(zero(max_ngwfs,max_ngwfs),stat=ierr)
       call utils_alloc_check('internal_remove_short_circuits','zero',ierr)
       zero(:,:) = cmplx_0

       ispin = greenf_obj%ispin

       ! for each lead
       do ilead = 1, pub_nleads

          ! loop over all atoms in this lead
          do col_at = pub_leads(ilead)%atms(1), pub_leads(ilead)%atms(2)
             col_orb = greenf_obj%ngwfs_on_atom(par%distr_atom(col_at))
             col_block = greenf_obj%block_of_atom(par%distr_atom(col_at))

             ! loop over all other atoms
             do row_at = pub_lcr%atms(1),pub_lcr%atms(2)

                row_orb = greenf_obj%ngwfs_on_atom(par%distr_atom(row_at))
                row_block = greenf_obj%block_of_atom(par%distr_atom(row_at))


                ! skip if these atoms should interact
                if ((row_at .ge. pub_leads(ilead)%atms(1) .and. &
                     row_at .le. pub_leads(ilead)%atms(2)) .or. &
                    (row_at .ge. pub_leads(ilead)%atms(3) .and. &
                     row_at .le. pub_leads(ilead)%atms(4))) cycle

                ! skip if these atoms already have zero overlap
                if (abs(col_block-row_block) .gt. 1) cycle

                ! set the h, s elements for these atoms to zero
                call trimat_deposit_atom_block(zero(1:col_orb,1:row_orb), &
                       greenf_obj%h(ispin), greenf_obj, col_at, col_at, row_at, row_at)
                call trimat_deposit_atom_block(zero(1:col_orb,1:row_orb), &
                       greenf_obj%s, greenf_obj, col_at, col_at, row_at, row_at)

                ! if in same block, must also remove other off-diagonal
                if (col_block == row_block) then
                   call trimat_deposit_atom_block(zero(1:row_orb,1:col_orb), &
                          greenf_obj%h(ispin), greenf_obj, row_at, row_at, col_at, col_at)
                   call trimat_deposit_atom_block(zero(1:row_orb,1:col_orb), &
                          greenf_obj%s, greenf_obj, row_at, row_at, col_at, col_at)
                endif
             enddo
          enddo
       enddo

       deallocate(zero,stat=ierr)
       call utils_dealloc_check('internal_remove_short_circuits','zero',ierr)


    end subroutine internal_remove_short_circuits

  end subroutine etrans_calculate

  !====================================================================!
  !====================================================================!

  subroutine etrans_init_hsm(lead_hsm,ham_type)

    !================================================================!
    ! Allocates the hsm arrays.                                      !
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in November 2011                 !
    !================================================================!

    use comms, only: comms_barrier
    use constants, only: stdout
    use rundat, only: pub_num_spins
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(hsm_type), allocatable, intent(inout) :: lead_hsm(:)
    character(len=*), intent(in)               :: ham_type

    ! Local variables
    integer    :: norb
    integer    :: il
    integer    :: ierr

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering etrans_init_hsm'

    !--------------------------------------------------------------------------!
    ! Allocate and create dense matrices
    !--------------------------------------------------------------------------!
    allocate(lead_hsm(pub_nleads),stat=ierr)
    call utils_alloc_check('etrans_init_hsm','lead_hsm',ierr)

    do il = 1, pub_nleads
       lead_hsm(il)%norb = pub_leads(il)%orbs(2) - pub_leads(il)%orbs(1) + 1
    enddo

    do il = 1, pub_nleads
       lead_hsm(il)%ham_type = ham_type
       norb = lead_hsm(il)%norb
       allocate(lead_hsm(il)%h00(norb,norb,pub_num_spins),stat=ierr)
       call utils_alloc_check('etrans_init_hsm','lead_hsm(il)%h00',ierr)
       allocate(lead_hsm(il)%h01(norb,norb,pub_num_spins),stat=ierr)
       call utils_alloc_check('etrans_init_hsm','lead_hsm(il)%h01',ierr)
       allocate(lead_hsm(il)%s00(norb,norb),stat=ierr)
       call utils_alloc_check('etrans_init_hsm','lead_hsm(il)%s00',ierr)
       allocate(lead_hsm(il)%s01(norb,norb),stat=ierr)
       call utils_alloc_check('etrans_init_hsm','lead_hsm(il)%s01',ierr)
       allocate(lead_hsm(il)%self(norb,norb),stat=ierr)
       call utils_alloc_check('etrans_init_hsm','lead_hsm(il)%self',ierr)
       allocate(lead_hsm(il)%self_rev(norb,norb),stat=ierr)
       call utils_alloc_check('etrans_init_hsm','lead_hsm(il)%self_rev',ierr)

       if (.not. allocated(pub_leads(il)%occ)) then
          allocate(pub_leads(il)%occ(pub_num_spins),stat=ierr)
          call utils_alloc_check('etrans_init_hsm','pub_leads(il)%occ',ierr)
       endif
       if (.not. allocated(pub_leads(il)%efermi)) then
          allocate(pub_leads(il)%efermi(pub_num_spins),stat=ierr)
          call utils_alloc_check('etrans_init_hsm','pub_leads(il)%efermi',ierr)
       endif

    enddo

    call comms_barrier

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving etrans_init_hsm'

  end subroutine etrans_init_hsm

  !====================================================================!
  !====================================================================!

  subroutine etrans_gather_lead_hsm(lead_hsm,ham,rep,ngwf_basis)

    !================================================================!
    ! Routine that gathers the required Hamiltonian/overlap blocks   !
    ! for the leads.                                                 !
    !                                                                !
    !----------------------------------------------------------------!
    ! Moved here by Robert Bell in April 2014, using code originally !
    ! by Simon M.-M. Dubois in November 2011                         !
    !================================================================!

    use comms, only: pub_on_root, comms_barrier
    use constants, only: stdout, NORMAL
    use function_basis, only: FUNC_BASIS
    use ngwf_representation, only: NGWF_HAM, NGWF_REP
    use rundat, only: pub_num_spins, pub_output_detail

    implicit none

    ! Arguments
    type(hsm_type), allocatable, intent(inout) :: lead_hsm(:)
    type(NGWF_HAM), intent(in)    :: ham
    type(NGWF_REP), intent(in)    :: rep
    type(FUNC_BASIS), intent(in)  :: ngwf_basis(1)

    ! Local variables
    integer    :: norb
    integer    :: il, is
    integer    :: iat1, iat2, iat3, iat4

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering etrans_gather_lead_hsm'

    if (pub_on_root .and. pub_output_detail >= NORMAL) then
       write(stdout,'(/a)',advance='no') ' Gathering lead matrix elements '
       if (any(pub_leads(:)%num_unit_cells>1)) write(stdout,*)
    endif

    ! Dense matrices associated with LEADS
    do il = 1, pub_nleads
       iat1 = pub_leads(il)%atms(1)
       iat2 = pub_leads(il)%atms(2)
       iat3 = pub_leads(il)%atms(3)
       iat4 = pub_leads(il)%atms(4)
       norb = lead_hsm(il)%norb

       call internal_get_submatrix(iat1,iat2,iat1,iat2, &
            rep%overlap%p,lead_hsm(il)%s00(:,:),norb,norb,ngwf_basis)

       call internal_get_submatrix(iat1,iat2,iat3,iat4, &
            rep%overlap%p,lead_hsm(il)%s01(:,:),norb,norb,ngwf_basis)
       do is = 1, pub_num_spins

          call internal_get_submatrix(iat1,iat2,iat1,iat2, &
               ham%ham(is)%p,lead_hsm(il)%h00(:,:,is),norb,norb,ngwf_basis)

          call internal_get_submatrix(iat1,iat2,iat3,iat4, &
               ham%ham(is)%p,lead_hsm(il)%h01(:,:,is),norb,norb,ngwf_basis)
       enddo

       ! symmetrise the leads if defined as more than one unit cell
       if (pub_leads(il)%num_unit_cells > 1) then
          if (pub_on_root .and. pub_output_detail >= NORMAL) &
               write(stdout,'(a,i4)') ' ...symmetrising lead ', il
          call internal_symmetrise(lead_hsm(il),pub_leads(il))
       endif

    enddo

    call comms_barrier

    if (pub_on_root .and. pub_output_detail >= NORMAL) &
         write(stdout,'(a/)') ' ... done'

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving etrans_gather_lead_hsm'

    return

    contains

    subroutine internal_symmetrise(hsm,lead_def)
       !================================================================!
       ! Routine that gives the leads the correct translational symmetry!
       ! for the number of unit cells.                                  !
       ! The h00/h01, s00/s01 lead matrices are broken into blocks      !
       ! defining interactions between different unit cells. Blocks that!
       ! should be identical are averaged over, and then deposited back !
       ! into h00/h01, s00/s01.                                         !
       !                                                                !
       ! The relative parity (sign) of the NGWFs is assumed the same.   !
       !                                                                !
       !----------------------------------------------------------------!
       ! Written by Robert Bell in July 2014                            !
       !================================================================!
       use constants, only: DP
       use rundat, only: pub_num_spins
       use utils, only: utils_alloc_check, utils_dealloc_check


       implicit none

       ! arguments
       type(hsm_type), intent(inout) :: hsm
       type(block_def), intent(in)   :: lead_def

       ! internal
       real(kind=DP), allocatable :: h_tmp(:,:,:,:), s_tmp(:,:,:)
       integer :: icol, jrow, iuc, juc, norb, norb_uc
       integer :: iorb, jorb, num_ucs
       integer :: ierr

       norb = hsm%norb
       num_ucs = lead_def%num_unit_cells
       norb_uc = norb / num_ucs



       ! build matrix that contains all the symmetrized elements
       allocate(h_tmp(norb_uc,norb_uc,pub_num_spins,num_ucs*2),stat=ierr)
       call utils_alloc_check('internal_symmetrise','h_tmp',ierr)
       allocate(s_tmp(norb_uc,norb_uc,num_ucs*2),stat=ierr)
       call utils_alloc_check('internal_symmetrise','s_tmp',ierr)

       h_tmp(:,:,:,:) = 0.0_DP
       s_tmp(:,:,:)   = 0.0_DP

       ! for each element, deposit in correct block of h_tmp/s_tmp
       do icol = 1,norb
          ! number of unit cells after 1st one (plus 1) this element belongs in
          ! (once symmetrised)
          iorb = modulo(icol-1,norb_uc)+1
          iuc = (icol-1) / norb_uc + 1

          ! deal with h00, s00 first
          do jrow = icol,norb
             jorb = modulo(jrow-1,norb_uc)+1
             juc = (jrow - 1) / norb_uc + 1


             ! deposit
             h_tmp(iorb,jorb,:,juc-iuc+1) = &
                  h_tmp(iorb,jorb,:,juc-iuc+1) + hsm%h00(icol,jrow,:)
             s_tmp(iorb,jorb,juc-iuc+1)   = &
                  s_tmp(iorb,jorb,juc-iuc+1)   + hsm%s00(icol,jrow)

          enddo

          ! then with h01, s01
          do jrow = 1,norb
             jorb = modulo(jrow-1,norb_uc)+1
             juc = (norb + jrow - 1) / norb_uc + 1

             ! deposit
             if (lead_def%forward) then
               h_tmp(iorb,jorb,:,juc-iuc+1) = &
                     h_tmp(iorb,jorb,:,juc-iuc+1) + hsm%h01(icol,jrow,:)
               s_tmp(iorb,jorb,juc-iuc+1)   = &
                     s_tmp(iorb,jorb,juc-iuc+1)   + hsm%s01(icol,jrow)
             else
               h_tmp(iorb,jorb,:,juc-iuc+1) = &
                     h_tmp(iorb,jorb,:,juc-iuc+1) + hsm%h01(jrow,icol,:)
               s_tmp(iorb,jorb,juc-iuc+1)   = &
                     s_tmp(iorb,jorb,juc-iuc+1)   + hsm%s01(jrow,icol)
             endif
          enddo
       enddo

       ! must average the summed matrix elements
       do iuc = 1, num_ucs + 1
          h_tmp(:,:,:,iuc) = h_tmp(:,:,:,iuc) / num_ucs
          s_tmp(:,:,iuc)   = s_tmp(:,:,iuc)   / num_ucs
       enddo
       ! some h01/s01 elements are averaged over fewer blocks
       do iuc = num_ucs+2, 2*num_ucs-1
          h_tmp(:,:,:,iuc) = h_tmp(:,:,:,iuc) / (2*num_ucs - iuc + 1)
          s_tmp(:,:,iuc)   = s_tmp(:,:,iuc)   / (2*num_ucs - iuc + 1)
       enddo


       ! build the new lead matrix elements using these symmetrised elements
       do icol = 1,norb
          ! number of unit cells after 1st one (plus 1) this element belongs in
          ! (once symmetrised)
          iorb = modulo(icol-1,norb_uc)+1
          iuc = (icol-1) / norb_uc + 1

          ! deal with h00, s00 first
          do jrow = icol,norb
             jorb = modulo(jrow-1,norb_uc)+1
             juc = (jrow - 1) / norb_uc + 1


             ! deposit (with lower diagonal)
             hsm%h00(icol,jrow,:) = h_tmp(iorb,jorb,:,juc-iuc+1)
             hsm%h00(jrow,icol,:) = h_tmp(iorb,jorb,:,juc-iuc+1)
             hsm%s00(icol,jrow)   = s_tmp(iorb,jorb,juc-iuc+1)
             hsm%s00(jrow,icol)   = s_tmp(iorb,jorb,juc-iuc+1)

          enddo

          ! then with h01, s01
          do jrow = 1,norb
             jorb = modulo(jrow-1,norb_uc)+1
             juc = (norb + jrow - 1) / norb_uc + 1

             ! deposit
             if (lead_def%forward) then
                hsm%h01(icol,jrow,:) = h_tmp(iorb,jorb,:,juc-iuc+1)
                hsm%s01(icol,jrow)   = s_tmp(iorb,jorb,juc-iuc+1)
             else
                hsm%h01(jrow,icol,:) = h_tmp(iorb,jorb,:,juc-iuc+1)
                hsm%s01(jrow,icol)   = s_tmp(iorb,jorb,juc-iuc+1)
             endif
          enddo
       enddo


       ! tidy up
       deallocate(h_tmp,stat=ierr)
       call utils_dealloc_check('internal_symmetrise','h_tmp',ierr)
       deallocate(s_tmp,stat=ierr)
       call utils_dealloc_check('internal_symmetrise','s_tmp',ierr)


    end subroutine internal_symmetrise

  end subroutine etrans_gather_lead_hsm

  !====================================================================!
  !====================================================================!

  subroutine etrans_destroy_lead_hsm(lead_hsm)

    !================================================================!
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in November 2011                 !
    !================================================================!

    use constants, only: stdout
    use utils, only: utils_dealloc_check

    implicit none


    ! Arguments
    type(hsm_type), allocatable, intent(inout) :: lead_hsm(:)

    ! Local variables
    integer    :: il
    integer    :: ierr

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering etrans_destroy_lead_hsm'

    !--------------------------------------------------------------------------!
    ! Deallocate  dense matrices
    !--------------------------------------------------------------------------!

    do il = 1, pub_nleads
       deallocate(lead_hsm(il)%h00,stat=ierr)
       call utils_dealloc_check('etrans_destroy_lead_hsm','lead_hsm(il)%h00',ierr)
       deallocate(lead_hsm(il)%h01,stat=ierr)
       call utils_dealloc_check('etrans_destroy_lead_hsm','lead_hsm(il)%h01',ierr)
       deallocate(lead_hsm(il)%s00,stat=ierr)
       call utils_dealloc_check('etrans_destroy_lead_hsm','lead_hsm(il)%s00',ierr)
       deallocate(lead_hsm(il)%s01,stat=ierr)
       call utils_dealloc_check('etrans_destroy_lead_hsm','lead_hsm(il)%s01',ierr)
       deallocate(lead_hsm(il)%self,stat=ierr)
       call utils_dealloc_check('etrans_destroy_lead_hsm','lead_hsm(il)%self',ierr)
       deallocate(lead_hsm(il)%self_rev,stat=ierr)
       call utils_dealloc_check('etrans_destroy_lead_hsm','lead_hsm(il)%self_rev',ierr)

    enddo

    deallocate(lead_hsm,stat=ierr)
    call utils_dealloc_check('etrans_destroy_lead_hsm','lead_hsm',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving etrans_destroy_lead_hsm'

  end subroutine etrans_destroy_lead_hsm

  !====================================================================!
  !====================================================================!

  subroutine etrans_compute_erefs(lead_hsm,ham,ham_type,rep,denskern, &
       & elements,cell,ngwf_basis,eref)

    !================================================================!
    ! Routine to compute the reference energies of the calculation.  !
    ! Computes the chemical potentials of the leads, and the         !
    ! reference energy of the lcr system, which is used to define    !
    ! the energy around which the transmission is calculated.        !
    !                                                                !
    ! eref is the returned lcr reference energy.                     !
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Robert Bell in Nov 2013                             !
    ! Using original code by Simon M.-M. Dubois                      !
    !================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout, hartree_in_evs, UP, DN, NORMAL
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use ngwf_representation, only: NGWF_HAM, NGWF_REP
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_etrans_lcr, pub_etrans_eref_method, &
         pub_etrans_calc_lead_pot, pub_etrans_eref, pub_num_spins, &
         pub_output_detail, pub_etrans_plot_eigchan
    use simulation_cell, only: CELL_INFO
    use sparse_embed, only: SPAM3_EMBED_ARRAY
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(hsm_type), intent(inout) :: lead_hsm(pub_nleads)
    type(NGWF_HAM), intent(in)    :: ham
    character(len=*), intent(in)  :: ham_type
    type(NGWF_REP), intent(in)    :: rep
    type(SPAM3_EMBED_ARRAY), intent(in) :: denskern
    type(FUNC_BASIS), intent(in)  :: ngwf_basis(1)
    type(ELEMENT), intent(in)     :: elements(par%nat)
    type(CELL_INFO), intent(in)   :: cell
    real(kind=DP), intent(out)    :: eref

    ! Local variables
    integer    :: npotentials
    integer    :: il
    character(len=80) :: loc_eref_method
    ! jd: Increased^ to 80, or else length mismatch between esdf_string in rundat, through pub_etrans_eref_method
    real(kind=DP), parameter :: elec_ionic_charge_diff_tolerance = 0.1_DP

    loc_eref_method = pub_etrans_eref_method

    ! compute lead potentials
    if ( (pub_etrans_calc_lead_pot) .or. (loc_eref_method .eq. 'LEADS') ) then

       call etrans_calculate_lead_ion_occ(pub_leads, elements)

       call etrans_calculate_lead_elec_occ(pub_leads,ham_type,rep%overlap, &
               denskern,ngwf_basis)

       do il = 1, pub_nleads
          call etrans_calculate_lead_potential(pub_leads(il),lead_hsm(il),il,cell)
       enddo

       ! check if all lead potentials failed
       if ( (count(pub_leads(:)%have_lead_mu) .eq. 0) .and. &
            & (loc_eref_method .eq. 'LEADS')) then
          if (pub_on_root) write(stdout,'(a)') "Failed to determine any lead potentials, &
               &resorting to diagonalisation to determine reference energy."
          loc_eref_method = 'DIAG'
       endif
    endif


    ! read in lead matrix elements, if required, and align to chemical potential
    do il=1,pub_nleads
       if (pub_leads(il)%ioread) then
          call etrans_read_lead_hsm(pub_leads(il),lead_hsm(il),il)
          ! recalculate lead potential and bands for this lead
          call etrans_calculate_lead_potential(pub_leads(il),lead_hsm(il),il,cell)
       endif
    enddo

    ! print warning if reading .hsm and plotting eigenchannels
    ! the lead NGWFs are not changed!
    if (pub_etrans_plot_eigchan .and. any(pub_leads(:)%ioread) .and. &
         pub_on_root) then
       write(stdout,'(/a)') ' Warning: plotting eigenchannels has been &
           &requested when reading lead matrix'
       write(stdout,'(a)') ' elements. The NGWFs corresponding to the leads &
           &are not changed however.'
       write(stdout,'(a)') ' Please ensure that the current NGWFs are similar, &
           &and in the same ordering,'
       write(stdout,'(a/)') ' to those those that generated the new lead matrix &
           &elements.'
    endif

    ! determine reference energy
    eref = 0.0_DP
    select case (loc_eref_method)

       case ('LEADS')
          if (pub_output_detail > NORMAL .and. pub_on_root)  write(stdout,'(/1x,a)') &
               &"Determining reference energy by average lead potential"
          npotentials = 0
          do il = 1, pub_nleads
             if (pub_leads(il)%have_lead_mu) then
                eref = eref + sum(pub_leads(il)%efermi(:))/pub_num_spins
                npotentials = npotentials + 1
             endif
          enddo
          eref = eref / npotentials

       case ('DIAG')
          if (pub_output_detail > NORMAL .and. pub_on_root) &
               & write(stdout,'(/1x,a)') "Determining reference energy by diagonalisation"
          call etrans_eref_diagonalisation(ham,rep,ngwf_basis,eref)

       case ('REFERENCE')
          if (pub_output_detail > NORMAL .and. pub_on_root) &
               & write(stdout,'(/1x,a)') "Using input reference energy"
          eref = pub_etrans_eref

       case default
          if (pub_on_root) write(stdout,'(/1x,a)') &
               &"Energy reference method unknown, defaulting to diagonalisation"
          call etrans_eref_diagonalisation(ham,rep,ngwf_basis,eref)

    end select



    ! Write info to output file
    if (pub_on_root .and. pub_output_detail >= NORMAL) then
       write(stdout,'(/,a)') '=== Summary of system'
       do il = 1, pub_nleads

          ! lead description
          write(stdout,'(/,2x,a,i2.2,a)') 'Bulk lead (',il,') :'
          write(stdout,'(4x,a,i5.5,a,i5.5,a)') &
               '     atoms [',pub_leads(il)%atms(1),',',pub_leads(il)%atms(2),']'
          write(stdout,'(4x,a,i5.5,a,i5.5,a)') &
               '  orbitals [',pub_leads(il)%orbs(1),',',pub_leads(il)%orbs(2),']'

          ! chemical potentials
          if (pub_leads(il)%have_lead_mu) then
             if (pub_num_spins == 2) then
                write(stdout,'(4x,a,f8.3,a,f8.3,a)') '   E_fermi [', &
                     pub_leads(il)%efermi(UP)*HARTREE_IN_EVS,',',&
                     pub_leads(il)%efermi(UP)*HARTREE_IN_EVS,'] eV'
             else
                write(stdout,'(4x,a,f8.3,a)') '   E_fermi :', &
                     pub_leads(il)%efermi(1)*HARTREE_IN_EVS, ' eV'
             endif
          else
             write(stdout,'(4x,a)') '   E_fermi : N/A'
          endif

          ! electronic/ionic info
          if (pub_etrans_calc_lead_pot .or. (loc_eref_method .eq. 'LEADS')) then

             ! electron charge
             write(stdout,'(4x,a,f8.3)',advance='no') &
                  ' elec popn :', sum(pub_leads(il)%occ(:))
             if (lead_hsm(il)%ham_type/='valence') then
                write(stdout,*) ' (assumed)'
             else
                write(stdout,*)
             endif

             ! spin
             if (pub_num_spins == 2) then
                write(stdout,'(4x,a,f8.3,a)',advance='no') &
                     ' elec spin :', pub_leads(il)%occ(UP) - pub_leads(il)%occ(DN), ' hbar'
                if (lead_hsm(il)%ham_type/='valence') then
                   write(stdout,*) ' (assumed)'
                else
                   write(stdout,*)
                endif
             endif

             ! ionic charge
             write(stdout,'(4x,a,f8.3, a)') &
                  'ion charge :', pub_leads(il)%ionic_popn

             ! rab: print warning if the electronic and ionic charges don't
             ! rab: well match
             if ( abs(sum(pub_leads(il)%occ(:)) - pub_leads(il)%ionic_popn) > &
                  & elec_ionic_charge_diff_tolerance) then
                write(stdout,'(a)') ' WARNING -- lead electronic occupancy is &
                     &significantly different to the'
                write(stdout,'(12x,a)') 'lead ionic charge which may indicate &
                     &under-converged lead/scattering'
                write(stdout,'(12x,a)') 'region distance, or a poorly converged &
                     &kernel.'
             endif
          endif
       enddo

       if (pub_etrans_lcr) then
          write(stdout,'(/,2x,a)') 'Total device :'
          write(stdout,'(4x,a,i5.5,a,i5.5,a)') &
               '     atoms [',pub_lcr%atms(1),',',pub_lcr%atms(2),']'
          write(stdout,'(4x,a,i5.5,a,i5.5,a,/)') &
               '  orbitals [',pub_lcr%orbs(1),',',pub_lcr%orbs(2),']'
          write(stdout,'(1x,a,f8.3,a)') &
               'Reference energy: ', eref*HARTREE_IN_EVS, ' eV'
       endif
    endif

  end subroutine etrans_compute_erefs

  !====================================================================!
  !====================================================================!

  subroutine etrans_init_leads
    !================================================================!
    ! Allocate arrays for leads. Currently called in rundat_blocks   !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois and Robert Bell                  !
    !================================================================!

    use constants, only: stdout
    use utils, only: utils_alloc_check, utils_assert

    implicit none


    integer :: ierr, il

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering etrans_init_leads'

    call utils_assert(.not. allocated(pub_leads),'Error in etrans_init_leads: &
         &lead already allocated')

    allocate(pub_leads(pub_nleads),stat=ierr)
    call utils_alloc_check('etrans_init_leads','pub_leads',ierr)

    ! allocate associated variables
    do il=1,pub_nleads
       ! initialise flag for if chemical potential calculation succeeded
       pub_leads(il)%have_lead_mu = .false.
    enddo

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving etrans_init_leads'

  end subroutine etrans_init_leads

  !====================================================================!
  !====================================================================!

  subroutine etrans_destroy_leads
    !================================================================!
    ! Deallocate the lead arrays                                     !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois and Robert Bell                  !
    !================================================================!

    use constants, only: stdout
    use utils, only: utils_dealloc_check

    implicit none

    integer :: ierr, il

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering etrans_destroy_leads'

    do il = 1, pub_nleads
       if (allocated(pub_leads(il)%occ)) then
          deallocate(pub_leads(il)%occ,stat=ierr)
          call utils_dealloc_check('etrans_destroy_leads','leads%occ',ierr)
       endif
       if (allocated(pub_leads(il)%efermi)) then
          deallocate(pub_leads(il)%efermi,stat=ierr)
          call utils_dealloc_check('etrans_destroy_leads','leads%efermi',ierr)
       endif
    enddo

    if (allocated(pub_leads)) then
       deallocate(pub_leads,stat=ierr)
       call utils_dealloc_check('etrans_destroy_leads','pub_leads',ierr)
    endif

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving etrans_destroy_leads'

  end subroutine etrans_destroy_leads

  !====================================================================!
  !====================================================================!

  subroutine etrans_init_setup(ngwf_basis)

    !================================================================!
    ! Counts the orbitals associated with the leads and LCR region.  !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in November 2011                 !
    !================================================================!

    use constants, only: stdout, CRLF
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_etrans_lcr
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)

    ! Local variables
    integer    :: orb_count
    integer    :: iat, il, orb_lead, orb_PL

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering etrans_init_setup'

    ! Identify the range of atoms/orbitals corresponding to the leads and
    ! the etrans setup
    orb_count = 1

    do iat = 1, par%nat
       ! total device
       if (pub_etrans_lcr) then
          if (iat == pub_lcr%atms(1)) pub_lcr%orbs(1) = orb_count
          if (iat == pub_lcr%atms(2)) pub_lcr%orbs(2) = orb_count + &
               ngwf_basis(1)%num_on_atom(par%distr_atom(iat)) - 1
       endif

       ! indiviual leads
       do il = 1, pub_nleads

          if (iat == pub_leads(il)%atms(1)) pub_leads(il)%orbs(1) = orb_count
          if (iat == pub_leads(il)%atms(2)) pub_leads(il)%orbs(2) = orb_count + &
               ngwf_basis(1)%num_on_atom(par%distr_atom(iat)) - 1
          if (iat == pub_leads(il)%atms(3)) pub_leads(il)%orbs(3) = orb_count
          if (iat == pub_leads(il)%atms(4)) pub_leads(il)%orbs(4) = orb_count + &
               ngwf_basis(1)%num_on_atom(par%distr_atom(iat)) - 1
       enddo
       orb_count = orb_count + ngwf_basis(1)%num_on_atom(par%distr_atom(iat))
    enddo

    ! check that the number of orbitals in the lead and 1st PL are the same
    do il = 1, pub_nleads
       orb_lead = pub_leads(il)%orbs(2) - pub_leads(il)%orbs(1) + 1
       orb_PL   = pub_leads(il)%orbs(4) - pub_leads(il)%orbs(3) + 1
       call utils_assert(orb_lead == orb_PL, 'Error in etrans_init_setup: &
           &for lead il: lead and 1st PL have different'//CRLF//'number of &
           &orbitals, il and orbitals to follow: ',il,orb_lead,orb_PL)
    enddo

    if (pub_etrans_lcr) then
       pub_lcr%orbs(3) = pub_lcr%orbs(1)
       pub_lcr%orbs(4) = pub_lcr%orbs(2)
    endif

    do il = 1, pub_nleads
       if (pub_leads(il)%orbs(1) .le. pub_leads(il)%orbs(3)) then
          pub_leads(il)%forward = .true.
       else
          pub_leads(il)%forward = .false.
       endif
    enddo

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving etrans_init_setup'

  end subroutine etrans_init_setup

  !====================================================================!
  !====================================================================!

  subroutine etrans_write_xyz(elements)

    use constants, only: stdout
    use ion, only: ELEMENT
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_etrans_lcr, pub_rootname
    use services, only: services_write_xyz


    implicit none

    ! arguments
    type(ELEMENT), intent(in) :: elements(par%nat)

    ! local
    integer            :: il, istart, iend
    character(len=128) :: xyz_root

    ! write out geometry to .xyz file
    if (pub_etrans_lcr) then
       istart = pub_lcr%atms(1)
       iend = pub_lcr%atms(2)
       write(xyz_root,'(2a)') trim(pub_rootname),'_device'
       call services_write_xyz(elements(istart:iend),trim(xyz_root), &
            'Device geometry')
    endif

    do il = 1, pub_nleads
       istart = pub_leads(il)%atms(1)
       iend = pub_leads(il)%atms(2)
       write(xyz_root,'(2a,i2.2)') trim(pub_rootname),'_lead',il
       call services_write_xyz(elements(istart:iend),trim(xyz_root), &
            'Lead geometry')
       istart = pub_leads(il)%atms(3)
       iend = pub_leads(il)%atms(4)
       write(xyz_root,'(2a,i2.2)') trim(pub_rootname),'_lead_pl',il
       call services_write_xyz(elements(istart:iend),trim(xyz_root), &
            'Principle layer geometry')
    enddo

    write(stdout,*)


  end subroutine etrans_write_xyz

  !====================================================================!
  !====================================================================!

  subroutine internal_get_submatrix(row_start,row_stop,col_start, &
       col_stop,msparse,mdense,nrow,ncol,ngwf_basis)

    !================================================================!
    ! Collects the submatrices from the sparse H/S matrices and      !
    ! creates the corresponding dense matrices for the LCR/leads.    !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in November 2011                 !
    !================================================================!

    use comms, only: pub_my_proc_id, comms_bcast, pub_total_num_procs, &
         comms_barrier
    use constants, only: DP
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: par=>pub_par
    use sparse, only: SPAM3, sparse_get_block
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    integer, intent(in)          :: row_start, row_stop, nrow
    integer, intent(in)          :: col_start, col_stop, ncol
    type(SPAM3), intent(in)      :: msparse
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    real(kind=DP), intent(inout)     :: mdense(nrow,ncol)

    ! Internal
    real(kind=DP), allocatable   :: send_buffer(:,:)
    real(kind=DP), allocatable   :: comms_buffer(:,:)
    real(kind=DP), allocatable   :: block_tmp(:,:)
    integer, allocatable         :: norb_col_loc(:)
    integer    :: iproc
    integer    :: ierr
    integer    :: iat, jat, iat_dist, jat_dist
    integer    :: norb_col, norb_row, iat_norb, jat_norb
    integer    :: iorb_loc, jorb_loc, iorb, jorb

    call comms_barrier()


    allocate(norb_col_loc(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('internal_get_submatrix','norb_col_loc',ierr)

    norb_col = 0
    norb_col_loc = 0
    do iat = col_start, col_stop
       iat_dist = par%distr_atom(iat)
       iat_norb = ngwf_basis(1)%num_on_atom(iat_dist)
       norb_col = norb_col + iat_norb
       norb_col_loc(par%proc_of_atom(iat)) = &
            norb_col_loc(par%proc_of_atom(iat)) + iat_norb
    enddo

    norb_row = 0
    do jat = row_start, row_stop
       jat_dist = par%distr_atom(jat)
       jat_norb = ngwf_basis(1)%num_on_atom(jat_dist)
       norb_row = norb_row + jat_norb
    enddo

    ! Check consistency
    if (nrow.ne.norb_row .OR. ncol.ne.norb_col) then
       call utils_abort('Wrong matrix size in internal_get_submatrix !',&
                        nrow, norb_row, ncol,norb_col)
    endif

    if (norb_col_loc(pub_my_proc_id) .gt. 0) then
       allocate(send_buffer(norb_row,norb_col_loc(pub_my_proc_id)),stat=ierr)
       call utils_alloc_check('internal_get_submatrix','send_buffer',ierr)
    endif

    call comms_barrier()
    iorb_loc = 0
    do iat = col_start, col_stop
       iat_dist = par%distr_atom(iat)
       if (pub_my_proc_id == par%proc_of_atom(iat)) then
          iat_norb = ngwf_basis(1)%num_on_atom(iat_dist)
          jorb_loc = 0
          do jat = row_start, row_stop
             jat_dist = par%distr_atom(jat)
             jat_norb = ngwf_basis(1)%num_on_atom(jat_dist)
             allocate(block_tmp(jat_norb,iat_norb),stat=ierr)
             call utils_alloc_check('internal_get_submatrix','block_tmp',ierr)
             call sparse_get_block(block_tmp,msparse,jat_dist,iat_dist)
             send_buffer(jorb_loc+1:jorb_loc+jat_norb, &
                  iorb_loc+1:iorb_loc+iat_norb)=block_tmp
             jorb_loc = jorb_loc + jat_norb
             deallocate(block_tmp,stat=ierr)
             call utils_dealloc_check('internal_get_submatrix','block_tmp',ierr)
          enddo
       else
          iat_norb = 0
       endif
       iorb_loc = iorb_loc + iat_norb
    enddo

    call comms_barrier()

    do iproc = 0, pub_total_num_procs-1

       if (norb_col_loc(iproc) .gt. 0) then

          allocate(comms_buffer(norb_row,norb_col_loc(iproc)),stat=ierr)
          call utils_alloc_check('internal_get_submatrix','comms_buffer',ierr)

          if (iproc == pub_my_proc_id) comms_buffer = send_buffer
          call comms_bcast(iproc,comms_buffer)
          call comms_barrier()

          iorb = 0
          jorb = 0
          do iat = col_start, col_stop
             iat_dist = par%distr_atom(iat)
             iat_norb = ngwf_basis(1)%num_on_atom(iat_dist)
             if (par%proc_of_atom(iat) == iproc) then
                mdense(:,iorb+1:iorb+iat_norb) = &
                     comms_buffer(:,jorb+1:jorb+iat_norb)
                jorb = jorb + iat_norb
             endif
             iorb = iorb + iat_norb
          enddo

          deallocate(comms_buffer,stat=ierr)
          call utils_dealloc_check('internal_get_submatrix','comms_buffer',ierr)

       endif

    enddo

    if (norb_col_loc(pub_my_proc_id) .gt. 0) then
       deallocate(send_buffer,stat=ierr)
       call utils_dealloc_check('internal_get_submatrix','send_buffer',ierr)
    endif

    deallocate(norb_col_loc,stat=ierr)
    call utils_dealloc_check('internal_get_submatrix','norb_col_loc',ierr)

  end subroutine internal_get_submatrix

  !====================================================================!
  !====================================================================!

  subroutine internal_get_submatrix_dem(row_start,row_stop,col_start, &
       col_stop,msparse,mdense,nrow,ncol,ngwf_basis)

    !================================================================!
    ! Collects the submatrices from the sparse H/S matrices and      !
    ! creates the corresponding dense matrices for the LCR/leads.    !
    ! Version for DEM dense matrix.                                  !
    !----------------------------------------------------------------!
    ! Written by Robert Bell, April 2014, based on code by           !
    ! Simon M.-M. Dubois                                             !
    !================================================================!

    use comms, only: pub_my_proc_id, comms_bcast, pub_total_num_procs, &
         comms_barrier
    use constants, only: DP
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: par=>pub_par
    use sparse, only: SPAM3, sparse_get_block
    use dense, only: DEM, dense_put_col
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    integer, intent(in)          :: row_start, row_stop, nrow
    integer, intent(in)          :: col_start, col_stop, ncol
    type(SPAM3), intent(in)      :: msparse
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    type(DEM), intent(inout)     :: mdense

    ! Internal
    real(kind=DP), allocatable   :: send_buffer(:,:)
    real(kind=DP), allocatable   :: comms_buffer(:,:)
    real(kind=DP), allocatable   :: block_tmp(:,:)
    integer, allocatable         :: norb_col_loc(:)
    integer    :: iproc
    integer    :: ierr
    integer    :: iat, jat, iat_dist, jat_dist
    integer    :: norb_col, norb_row, iat_norb, jat_norb
    integer    :: iorb_loc, jorb_loc, iorb, jorb, ingwf

    call comms_barrier()


    allocate(norb_col_loc(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('internal_get_submatrix_dem','norb_col_loc',ierr)

    norb_col = 0
    norb_col_loc = 0
    do iat = col_start, col_stop
       iat_dist = par%distr_atom(iat)
       iat_norb = ngwf_basis(1)%num_on_atom(iat_dist)
       norb_col = norb_col + iat_norb
       norb_col_loc(par%proc_of_atom(iat)) = &
            norb_col_loc(par%proc_of_atom(iat)) + iat_norb
    enddo

    norb_row = 0
    do jat = row_start, row_stop
       jat_dist = par%distr_atom(jat)
       jat_norb = ngwf_basis(1)%num_on_atom(jat_dist)
       norb_row = norb_row + jat_norb
    enddo

    ! Check consistency
    if (nrow.ne.norb_row .OR. ncol.ne.norb_col) then
       call utils_abort('Wrong matrix size in internal_get_submatrix_dem !', &
                        nrow, norb_row, ncol,norb_col)
    endif

    if (norb_col_loc(pub_my_proc_id) .gt. 0) then
       allocate(send_buffer(norb_row,norb_col_loc(pub_my_proc_id)),stat=ierr)
       call utils_alloc_check('internal_get_submatrix_dem','send_buffer',ierr)
    endif

    call comms_barrier()
    iorb_loc = 0
    do iat = col_start, col_stop
       iat_dist = par%distr_atom(iat)
       if (pub_my_proc_id == par%proc_of_atom(iat)) then
          iat_norb = ngwf_basis(1)%num_on_atom(iat_dist)
          jorb_loc = 0
          do jat = row_start, row_stop
             jat_dist = par%distr_atom(jat)
             jat_norb = ngwf_basis(1)%num_on_atom(jat_dist)
             allocate(block_tmp(jat_norb,iat_norb),stat=ierr)
             call utils_alloc_check('internal_get_submatrix_dem','block_tmp',ierr)
             call sparse_get_block(block_tmp,msparse,jat_dist,iat_dist)
             send_buffer(jorb_loc+1:jorb_loc+jat_norb, &
                  iorb_loc+1:iorb_loc+iat_norb)=block_tmp
             jorb_loc = jorb_loc + jat_norb
             deallocate(block_tmp,stat=ierr)
             call utils_dealloc_check('internal_get_submatrix_dem','block_tmp',ierr)
          enddo
       else
          iat_norb = 0
       endif
       iorb_loc = iorb_loc + iat_norb
    enddo

    call comms_barrier()

    do iproc = 0, pub_total_num_procs-1

       if (norb_col_loc(iproc) .gt. 0) then

          allocate(comms_buffer(norb_row,norb_col_loc(iproc)),stat=ierr)
          call utils_alloc_check('internal_get_submatrix_dem','comms_buffer',ierr)

          if (iproc == pub_my_proc_id) comms_buffer = send_buffer
          call comms_bcast(iproc,comms_buffer)
          call comms_barrier()

          iorb = 0
          jorb = 0
          do iat = col_start, col_stop
             iat_dist = par%distr_atom(iat)
             iat_norb = ngwf_basis(1)%num_on_atom(iat_dist)
             if (par%proc_of_atom(iat) == iproc) then
                do ingwf=1,iat_norb
                   call dense_put_col(comms_buffer(:,jorb+ingwf),mdense,iorb+ingwf)
                enddo
                jorb = jorb + iat_norb
             endif
             iorb = iorb + iat_norb
          enddo

          deallocate(comms_buffer,stat=ierr)
          call utils_dealloc_check('internal_get_submatrix_dem','comms_buffer',ierr)

       endif

    enddo

    if (norb_col_loc(pub_my_proc_id) .gt. 0) then
       deallocate(send_buffer,stat=ierr)
       call utils_dealloc_check('internal_get_submatrix_dem','send_buffer',ierr)
    endif

    deallocate(norb_col_loc,stat=ierr)
    call utils_dealloc_check('internal_get_submatrix_dem','norb_col_loc',ierr)

  end subroutine internal_get_submatrix_dem

  !====================================================================!
  !====================================================================!

  subroutine etrans_init_energy_contour(econt,contype)

    !================================================================!
    ! Determines and distributes across all procs the energy points  !
    ! the transmission is to be calculated at.                       !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in November 2011                 !
    !================================================================!

    use comms, only: pub_on_root, pub_total_num_procs
    use constants, only: DP, stdout, hartree_in_evs, cmplx_1
    use rundat, only: pub_etrans_ecmplx, pub_etrans_enum , &
         pub_etrans_emin, pub_etrans_emax
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(epath), intent(inout)   :: econt
    character(len=*), intent(in)    :: contype

    ! Local variables
    real(kind=DP) :: estart, estop, estep
    integer       :: ipt, ierr, iproc
    integer       :: eblock, rem, ptr_offset

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &etrans_init_energy_contour'

    if (contype .eq. 'trans') then

       econt%type = 'TR'
       econt%seg(1) = pub_etrans_enum
       econt%seg(2:4) = 0

       econt%nep = sum(econt%seg(:))

       allocate(econt%ep(econt%nep),stat=ierr)
       call utils_alloc_check('etrans_init_energy_contour','econt%ep',ierr)
       allocate(econt%ew(econt%nep),stat=ierr)
       call utils_alloc_check('etrans_init_energy_contour','econt%ew',ierr)


       estart = pub_etrans_emin + econt%eref
       estop  = pub_etrans_emax + econt%eref
       estep  = 0.0_DP
       if (pub_etrans_enum > 1) then
          estep  = (estop-estart)/(pub_etrans_enum-1)
       endif

       do ipt = 1, pub_etrans_enum
          econt%ew(ipt) = cmplx_1
          econt%ep(ipt) = cmplx(estart+real(ipt-1,kind=DP)*estep, &
               pub_etrans_ecmplx,kind=DP)
       enddo

       ! Distribute energy contour
       allocate(econt%nepproc(0:pub_total_num_procs-1),stat=ierr)
       call utils_alloc_check('etrans_init_energy_contour','econt%nepproc',ierr)
       econt%nepproc = 0
       allocate(econt%nepptr(0:pub_total_num_procs-1),stat=ierr)
       call utils_alloc_check('etrans_init_energy_contour','econt%nepptr',ierr)
       econt%nepptr = 0

       eblock = floor(real(econt%nep)/pub_total_num_procs)
       rem = econt%nep - eblock*pub_total_num_procs

       ptr_offset = 0 ! tally the current number of energy points distributed
       do iproc =  0, rem-1
          econt%nepptr(iproc) = ptr_offset
          ptr_offset = ptr_offset + (eblock+1)
          econt%nepproc(iproc) = eblock+1
       enddo

       do iproc = rem, pub_total_num_procs-1
          econt%nepptr(iproc) = ptr_offset
          ptr_offset = ptr_offset + eblock
          econt%nepproc(iproc) = eblock
       enddo


       ! Write this information into the main output file
       if (pub_on_root) then
          write(stdout,'(1x,a,2(f8.4,a))') 'Energy range = [', &
               estart*hartree_in_evs,',',estop*hartree_in_evs,'] eV'
          write(stdout,'(1x,a,i6)') 'Number of energy points = ', econt%nep
          write(stdout,'(1x,a,2(i5,1x),a)') 'Number of energy point per &
               &proc (max, min)  = (', &
               maxval(econt%nepproc(:)),minval(econt%nepproc(:)),')'
       endif

    endif

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &etrans_init_energy_contour'

  end subroutine etrans_init_energy_contour

  !====================================================================!
  !====================================================================!

  subroutine etrans_free_energy_contour(econt)

    !================================================================!
    ! Deallocates the energy contour.                                !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in November 2011                 !
    !================================================================!

    use constants, only: stdout
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(epath), intent(inout)   :: econt

    ! Local variables
    integer       :: ierr

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering etrans_free_energy_contour'

    econt%type = ''
    econt%nep = 0
    econt%seg(:) = 0

    deallocate(econt%ep,stat=ierr)
    call utils_dealloc_check('etrans_free_energy_contour','econt%ep',ierr)
    deallocate(econt%ew,stat=ierr)
    call utils_dealloc_check('etrans_free_energy_contour','econt%ew',ierr)
    deallocate(econt%nepproc,stat=ierr)
    call utils_dealloc_check('etrans_free_energy_contour','econt%nepproc',ierr)
    deallocate(econt%nepptr,stat=ierr)
    call utils_dealloc_check('etrans_free_energy_contour','econt%nepptr',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving etrans_free_energy_contour'

  end subroutine etrans_free_energy_contour

  !====================================================================!
  !====================================================================!

  subroutine etrans_free_eigchan_contour
    !================================================================!
    ! Deallocates the eigenchannel energy contour.                   !
    !----------------------------------------------------------------!
    ! Written by Robert Bell, June 2014                              !
    !================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! internal
    integer :: ierr

    if (allocated(pub_etrans_eigchan_en)) then
       deallocate(pub_etrans_eigchan_en,stat=ierr)
       call utils_dealloc_check('etrans_free_eigchan_contour',&
            'pub_etrans_eigchan_en',ierr)
    endif

  end subroutine etrans_free_eigchan_contour

  !====================================================================!
  !====================================================================!

  subroutine etrans_compute_transmission_channels(greenf_obj,lead_hsm,trc, &
      info, compute_eigchans,trc_chan)

    !================================================================!
    ! Calculate the transmission channels for states injected by     !
    ! each lead.                                                     !
    ! If compute_eigchans = .false., only the transmission           !
    ! are calculated.                                                !
    !                                                                !
    ! Works by calculating the transmission matrix for each lead and !
    ! diagonalising to get the transmission eigenvalues:             !
    !                                                                !
    !   T_i := diag(G_a Gamma_R G_r Gamma_L)                         !
    !                                                                !
    ! The transmission matrices T = (t^dagger t)                     !
    !                            ~= (G_a Gamma_R G_r Gamma_L)        !
    ! are calculated for each lead pair, and automatically summed    !
    ! over to give the total transmission matrix for that lead. This !
    ! is diagonalised to give the eigenchannel transmissions.        !
    !----------------------------------------------------------------!
    ! Written by Robert Bell, April 2014                             !
    ! Based on code by Simon M.-M. Dubois                            !
    !================================================================!

    use constants, only: DP, stdout, cmplx_0
    use greenf, only: SPARGF, greenf_get_atom_block
    use rundat, only: pub_etrans_num_eigchan
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! arguments
    type(SPARGF), intent(in) :: greenf_obj
    type(hsm_type), intent(in)    :: lead_hsm(pub_nleads)
    real(kind=DP), intent(inout)  :: trc(pub_nleads*(pub_nleads-1))
    integer, intent(out)          :: info
    real(kind=DP), optional, intent(inout)  :: trc_chan(pub_etrans_num_eigchan,pub_nleads)
    logical, optional, intent(in)           :: compute_eigchans

    ! Local variables
    complex(kind=DP), allocatable :: gf_block(:,:)
    integer :: itrc, il, jl, ierr
    integer :: norbl, norbr
    integer :: lshift, rshift
    logical :: compute_eigchans_loc
    type(tmat_type), allocatable :: T(:) ! transmission matrices


    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering etrans_compute_transmission_channels'

    compute_eigchans_loc = .false.
    if (present(compute_eigchans)) compute_eigchans_loc = compute_eigchans

    ! allocate transmission matrix
    if (compute_eigchans_loc) then
       allocate(T(pub_nleads),stat=ierr)
       call utils_alloc_check('etrans_compute_transmission_channels','T',ierr)
       do il=1,pub_nleads
          norbl = lead_hsm(il)%norb
          T(il)%norb = norbl
          allocate(T(il)%T(norbl,norbl),stat=ierr)
          call utils_alloc_check('etrans_compute_transmission_channels','T(il)%T',ierr)
          T(il)%T = cmplx_0
       enddo
    endif

    ! Compute transmission coefficients
    itrc = 0
    do il = 1,pub_nleads-1
       norbl = lead_hsm(il)%norb
       lshift = pub_leads(il)%orbs(1) - pub_lcr%orbs(1)

       do jl = il+1,pub_nleads
          itrc = itrc + 1

          norbr = lead_hsm(jl)%norb
          rshift = pub_leads(jl)%orbs(1) - pub_lcr%orbs(1)

          allocate(gf_block(norbr,norbl),stat=ierr)
          call utils_alloc_check('etrans_compute_transmission_channels','gf_block',ierr)

          call greenf_get_atom_block(gf_block,greenf_obj, &
                 pub_leads(jl)%atms(1), pub_leads(jl)%atms(2), &
                 pub_leads(il)%atms(1), pub_leads(il)%atms(2), info)

          if (info .ne. 0) goto 100

          ! compute transmission matrix contributions
          if (compute_eigchans_loc) then
             call etrans_compute_transmission(norbl,lead_hsm(il)%self, norbr, &
                  lead_hsm(jl)%self,gf_block,trc(itrc),compute_eigchans_loc,T(il),T(jl))
          else
             call etrans_compute_transmission(norbl,lead_hsm(il)%self, norbr, &
                  lead_hsm(jl)%self,gf_block,trc(itrc),compute_eigchans_loc)
          endif

          deallocate(gf_block,stat=ierr)
          call utils_dealloc_check('etrans_compute_transmission_channels','gf_block',ierr)

       enddo
    enddo

    ! calculate eigenchannels
    if (compute_eigchans_loc) then
       do il=1,pub_nleads
         call internal_compute_eigtrans(T(il),trc_chan(:,il))
       enddo
    endif

 100 continue ! if info /= 0

    ! tidy up
    if (compute_eigchans_loc) then
       do il=1,pub_nleads
         deallocate(T(il)%T,stat=ierr)
         call utils_dealloc_check('etrans_compute_transmission_channels','T(il)%T',ierr)
       enddo
       deallocate(T,stat=ierr)
       call utils_dealloc_check('etrans_compute_transmission_channels','T',ierr)
    endif

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving etrans_compute_transmission_channels'

     contains

     subroutine internal_compute_eigtrans(T,trc_chan)
        !================================================================!
        ! Determine the eigenchannel transmission by diagonalising the   !
        ! transmission matrix for this lead.                             !
        ! On exit, T%T is overwritten and contains the (right)           !
        ! eigenvectors of T.                                             !
        !----------------------------------------------------------------!
        ! Written by Robert Bell, April 2014                             !
        !================================================================!

        use linalg, only: linalg_diag_serial
        use rundat, only: pub_etrans_num_eigchan
        use utils, only: utils_heapsort, utils_alloc_check, &
             utils_dealloc_check

        implicit none

        ! arguments
        type(tmat_type), intent(inout) :: T
        real(kind=DP), intent(out)     :: trc_chan(pub_etrans_num_eigchan)

        ! internal
        complex(kind=DP), allocatable :: eigvals_cmplx(:)
        real(kind=DP), allocatable    :: eigvals(:)
        integer :: ierr, norb, n_eigchans, ichan
        integer, allocatable :: sort_ndx(:)

        norb = T%norb
        allocate(sort_ndx(norb),stat=ierr)
        call utils_alloc_check('internal_compute_eigtrans','sort_ndx',ierr)
        allocate(eigvals_cmplx(norb),stat=ierr)
        call utils_alloc_check('internal_compute_eigtrans','eigvals_cmplx',ierr)
        allocate(eigvals(norb),stat=ierr)
        call utils_alloc_check('internal_compute_eigtrans','eigvals',ierr)

        ! compute eigenvalues of T matrix
        call linalg_diag_serial(T%T,eigvals_cmplx,norb)
        eigvals(:) = real(eigvals_cmplx(:),kind=DP)

        ! sort eigenvalues
        call utils_heapsort(norb,eigvals,sort_ndx)

        ! pick largest eigenvalues
        n_eigchans = min(norb,pub_etrans_num_eigchan)
        do ichan=1,n_eigchans
           trc_chan(ichan) = eigvals(norb-ichan+1)
        enddo


        deallocate(sort_ndx,stat=ierr)
        call utils_dealloc_check('internal_compute_eigtrans','sort_ndx',ierr)
        deallocate(eigvals_cmplx,stat=ierr)
        call utils_dealloc_check('internal_compute_eigtrans','eigvals_cmplx',ierr)
        deallocate(eigvals,stat=ierr)
        call utils_dealloc_check('internal_compute_eigtrans','eigvals',ierr)

     end subroutine internal_compute_eigtrans

  end subroutine etrans_compute_transmission_channels

  !====================================================================!
  !====================================================================!

  subroutine etrans_compute_transmission(norbl,selfl,norbr,selfr,greenf01,trc, &
      compute_Tmat,T_L,T_R)

    !================================================================!
    ! Calculate the transmission coefficient                         !
    !   T := tr(G_a Gamma_L G_r Gamma_R)                             !
    ! If compute_Tmat = .true., the transmission matrices for each   !
    ! lead are returned in T_L, T_R                                  !
    ! Assumes time-reversal symmetry                                 !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in April 2010                    !
    ! Modified to return transmission matrix, Robert Bell April 2014 !
    !================================================================!

    use constants, only: dp, stdout, cmplx_1, cmplx_i
    use linalg, only: linalg_mat_mul_serial
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in)            :: norbr, norbl
    complex(kind=DP), intent(inout):: greenf01(norbr,norbl)
    complex(kind=DP), intent(in)   :: selfl(norbl,norbl)
    complex(kind=DP), intent(in)   :: selfr(norbr,norbr)
    real(kind=dp), intent(out)     :: trc
    type(tmat_type), intent(inout), optional :: T_L
    type(tmat_type), intent(inout), optional :: T_R
    logical, intent(in), optional  :: compute_Tmat

    ! Local variables
    complex(kind=dp), allocatable  :: gammal(:,:), gammar(:,:)
    complex(kind=dp), allocatable  :: trcl(:,:), trcr(:,:)
    integer :: iorb, jorb, ierr
    logical :: compute_Tmat_loc


    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering etrans_compute_transmission'

    allocate(gammal(norbl,norbl),stat=ierr)
    call utils_alloc_check('etrans_compute_transmission','gammal',ierr)
    allocate(gammar(norbr,norbr),stat=ierr)
    call utils_alloc_check('etrans_compute_transmission','gammar',ierr)
    allocate(trcl(norbr,norbl),stat=ierr)
    call utils_alloc_check('etrans_compute_transmission','trcl',ierr)
    allocate(trcr(norbl,norbr),stat=ierr)
    call utils_alloc_check('etrans_compute_transmission','trcr',ierr)

    ! compute Gamma matrices
    gammal = cmplx_i*(selfl - conjg(transpose(selfl)))
    gammar = cmplx_i*(selfr - conjg(transpose(selfr)))

    call linalg_mat_mul_serial(trcl,greenf01,gammal)
    call linalg_mat_mul_serial(trcr,greenf01,gammar,opA='C')

    ! compute transmission
    trc = 0.0_dp
    do iorb=1,norbr
       do jorb=1,norbl
          trc = trc + real(trcl(iorb,jorb)*trcr(jorb,iorb),kind=DP)
       enddo
    enddo


    compute_Tmat_loc = .false.
    if (present(compute_Tmat)) compute_Tmat_loc = compute_Tmat

    ! if requested, calculate the transmission matrices
    if (compute_Tmat_loc .and. present(T_R)) then

       ! beta = cmplx_1 so that each contribution is summed over all leads
       call linalg_mat_mul_serial(T_R%T,trcl,trcr,beta=cmplx_1)

    endif

    ! to compute T_L = G_a Gamma_R G_r Gamma_L, a different block of the
    ! Green's function is required. Previously, G_r = G_LR, now need block
    ! G_r = G_RL.
    ! time-reversal symmetry G means that G_LR = transpose(G_RL),
    ! and Ga = G_LR^{dagger} = conjg(G_RL)
    ! Therefore trcl := G_LR^{dagger}.Gamma_L = conjg(G_RL).Gamma_L
    ! and       trcr := G_LR.Gamma_R          = tranpose(G_RL).Gamma_R
    if (compute_Tmat_loc .and. present(T_L)) then

       ! take conjg not transpose: no need to copy memory
       greenf01 = conjg(greenf01)

       ! compute trc matrices G*gamma
       call linalg_mat_mul_serial(trcr,greenf01,gammar,opA='C')
       call linalg_mat_mul_serial(trcl,greenf01,gammal)

       ! compute transmission matrix
       ! beta = cmplx_1 so that each contribution is summed over all leads
       call linalg_mat_mul_serial(T_L%T,trcr,trcl,beta=cmplx_1)

       ! undo conjugation
       greenf01 = conjg(greenf01)
    endif

    deallocate(gammal,stat=ierr)
    call utils_dealloc_check('etrans_compute_transmission','gammal',ierr)
    deallocate(gammar,stat=ierr)
    call utils_dealloc_check('etrans_compute_transmission','gammar',ierr)
    deallocate(trcl,stat=ierr)
    call utils_dealloc_check('etrans_compute_transmission','trcl',ierr)
    deallocate(trcr,stat=ierr)
    call utils_dealloc_check('etrans_compute_transmission','trcr',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving etrans_compute_transmission'

  end subroutine etrans_compute_transmission

  !====================================================================!
  !====================================================================!

  subroutine etrans_compute_self(lead_hsm,ecmp,spin,info)

    !================================================================!
    !                                                                !
    ! This subroutine compute the lead self-energies by means of     !
    ! iterative improvements of the bulk transfer matrix             !
    ! see Lopez-Sancho et al., J. Phys. F: Met.Phys. 15, 851 (1985)  !
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in April 2010                    !
    !================================================================!

    use constants, only: dp, stdout
    use linalg, only: linalg_mat_mul_serial
    use rundat, only: pub_etrans_same_leads
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(hsm_type), intent(inout):: lead_hsm(pub_nleads)
    complex(kind=dp), intent(in) :: ecmp
    integer, intent(in)          :: spin
    integer, intent(out)         :: info

    ! Internal variables
    complex(kind=dp), allocatable   :: tmp1(:,:)
    complex(kind=dp), allocatable   :: transf(:,:), transf_bar(:,:)
    complex(kind=dp), allocatable   :: self1(:,:), self2(:,:)

    integer  :: norb
    integer  :: il
    integer  :: ierr

    real(kind=dp), parameter    :: conv_tol=1.0e-7_dp

    ! smmdebug
    !integer  :: io_unit, io, jo

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering compute_self'

    !====================================================!
    ! Compute : Self_L = (e*S01-H01)^{dag} * TRANSF_bar  !
    !           Self_R = (e*S01-H01) * TRANSF            !
    !====================================================!


    if (pub_etrans_same_leads) then

       norb = lead_hsm(1)%norb
       allocate(transf(norb,norb),stat=ierr)
       call utils_alloc_check('compute_self','transf',ierr)
       allocate(transf_bar(norb,norb),stat=ierr)
       call utils_alloc_check('compute_self','transf_bar',ierr)
       allocate(self1(norb,norb),stat=ierr)
       call utils_alloc_check('compute_self','self1',ierr)
       allocate(self2(norb,norb),stat=ierr)
       call utils_alloc_check('compute_self','self2',ierr)

       ! Compute the transfer matrices
       call etrans_compute_transfer(transf,transf_bar, &
            lead_hsm(1)%h00(:,:,spin),lead_hsm(1)%h01(:,:,spin), &
            lead_hsm(1)%s00,lead_hsm(1)%s01, ecmp,norb,conv_tol,50,info)
       if(info/=0) return

       ! Compute the left and right self energies
       allocate(tmp1(norb,norb),stat=ierr)
       call utils_alloc_check('compute_self','tmp1',ierr)

       tmp1 = cmplx(lead_hsm(1)%h01(:,:,spin),kind=dp) &
           -ecmp*cmplx(lead_hsm(1)%s01,kind=dp)

       call linalg_mat_mul_serial(self1,tmp1,transf_bar,opA='C')
       call linalg_mat_mul_serial(self2,tmp1,transf)

       deallocate(tmp1,stat=ierr)
       call utils_dealloc_check('compute_self','tmp1',ierr)

       do il = 1, pub_nleads
          if (pub_leads(il)%forward .eqv. pub_leads(1)%forward) then
             lead_hsm(il)%self = self1
             lead_hsm(il)%self_rev = self2
          else
             lead_hsm(il)%self = self2
             lead_hsm(il)%self_rev = self1
          endif
       enddo

       deallocate(transf,stat=ierr)
       call utils_dealloc_check('compute_self','transf',ierr)
       deallocate(transf_bar,stat=ierr)
       call utils_dealloc_check('compute_self','transf_bar',ierr)
       deallocate(self1,stat=ierr)
       call utils_dealloc_check('compute_self','self1',ierr)
       deallocate(self2,stat=ierr)
       call utils_dealloc_check('compute_self','self2',ierr)

    else

       do il = 1, pub_nleads

          norb = lead_hsm(il)%norb
          allocate(transf(norb,norb),stat=ierr)
          call utils_alloc_check('compute_self','transf',ierr)
          allocate(transf_bar(norb,norb),stat=ierr)
          call utils_alloc_check('compute_self','transf_bar',ierr)
          allocate(self1(norb,norb),stat=ierr)
          call utils_alloc_check('compute_self','self1',ierr)
          allocate(self2(norb,norb),stat=ierr)
          call utils_alloc_check('compute_self','self2',ierr)

          ! Compute the transfer matrices
          call etrans_compute_transfer(transf,transf_bar, &
               lead_hsm(il)%h00(:,:,spin),lead_hsm(il)%h01(:,:,spin), &
               lead_hsm(il)%s00,lead_hsm(il)%s01, ecmp,norb,conv_tol,50,info)
          if(info/=0) return

          ! Compute the left and right self energies
          allocate(tmp1(norb,norb),stat=ierr)
          call utils_alloc_check('compute_self','tmp1',ierr)

          tmp1 = cmplx(lead_hsm(il)%h01(:,:,spin),kind=dp) &
              -ecmp*cmplx(lead_hsm(il)%s01,kind=dp)


          call linalg_mat_mul_serial(self1,tmp1,transf_bar,opA='C')
          call linalg_mat_mul_serial(self2,tmp1,transf)

          deallocate(tmp1,stat=ierr)
          call utils_dealloc_check('compute_self','tmp1',ierr)

          lead_hsm(il)%self = self1
          lead_hsm(il)%self_rev = self2

          deallocate(transf,stat=ierr)
          call utils_dealloc_check('compute_self','transf',ierr)
          deallocate(transf_bar,stat=ierr)
          call utils_dealloc_check('compute_self','transf_bar',ierr)
          deallocate(self1,stat=ierr)
          call utils_dealloc_check('compute_self','self1',ierr)
          deallocate(self2,stat=ierr)
          call utils_dealloc_check('compute_self','self2',ierr)

       enddo

    endif

    return

  end subroutine etrans_compute_self

  !====================================================================!
  !====================================================================!

  subroutine etrans_compute_transfer(transf,transf_bar,h_00,h_01,s_00,s_01, &
       energy,norb,convtol,maxiter,info)

    !================================================================!
    !                                                                !
    ! This subroutine compute the lead self-energies by means of     !
    ! iterative improvements of the bulk transfer matrix             !
    ! see Lopez-Sancho et al., J. Phys. F: Met.Phys. 15, 851 (1985)  !
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in April 2010                    !
    !================================================================!

    use constants, only: dp, cmplx_0, cmplx_1, stdout
    use linalg, only: linalg_mat_mul_serial
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_isnan

    implicit none

    ! Arguments
    integer, intent(in)           :: norb
    integer, intent(in)           :: maxiter
    complex(kind=dp), intent(in)  :: energy
    complex(kind=dp), intent(inout) ::  transf(norb,norb)
    complex(kind=dp), intent(inout) ::  transf_bar(norb,norb)
    real(kind=dp), intent(in)     :: h_00(norb,norb)
    real(kind=dp), intent(in)     :: h_01(norb,norb)
    real(kind=dp), intent(in)     :: s_00(norb,norb)
    real(kind=dp), intent(in)     :: s_01(norb,norb)
    real(kind=dp), intent(in)     :: convtol
    integer, intent(out)          :: info

    ! BLAS and LAPACK subroutines
    external :: zaxpy, zgesv

    ! Internal variables
    real(kind=dp)    :: conver,conver2
    complex(kind=dp), allocatable :: tprod(:,:), tprod_bar(:,:)
    complex(kind=dp), allocatable :: t0(:,:), t0_bar(:,:)
    complex(kind=dp), allocatable :: t1(:,:), t1_bar(:,:)
    complex(kind=dp), allocatable :: tmp1(:,:), tmp2(:,:), tmp3(:,:)

    real(kind=DP), parameter :: max_conv = 1.0D10
    integer  :: norb2
    integer  :: iorb, jorb, iter, ierr
    integer  :: ipiv(norb)

    ! allocate workspace
    allocate(tprod(norb,norb),stat=ierr)
    call utils_alloc_check('etrans_compute_transfer','tprod',ierr)
    allocate(tprod_bar(norb,norb),stat=ierr)
    call utils_alloc_check('etrans_compute_transfer','tprod_bar',ierr)
    allocate(t0(norb,norb),stat=ierr)
    call utils_alloc_check('etrans_compute_transfer','t0',ierr)
    allocate(t0_bar(norb,norb),stat=ierr)
    call utils_alloc_check('etrans_compute_transfer','t0_bar',ierr)
    allocate(t1(norb,norb),stat=ierr)
    call utils_alloc_check('etrans_compute_transfer','t1',ierr)
    allocate(t1_bar(norb,norb),stat=ierr)
    call utils_alloc_check('etrans_compute_transfer','t1_bar',ierr)
    allocate(tmp1(norb,norb),stat=ierr)
    call utils_alloc_check('etrans_compute_transfer','tmp1',ierr)
    allocate(tmp2(norb,norb),stat=ierr)
    call utils_alloc_check('etrans_compute_transfer','tmp2',ierr)
    allocate(tmp3(norb,norb),stat=ierr)
    call utils_alloc_check('etrans_compute_transfer','tmp3',ierr)



    !================================================================!
    !=== Compute t0 and t0_bar :
    !===   t0 =  k00^{-1} * k01^{dag}
    !===   t0_bar =  k00^{-1} * k01
    !================================================================!

    norb2 = norb*norb

    ! compute k00=(e*s00 - h00) --> tmp2
    !         k01=-(e*s01 - h01) --> tmp3
    tmp2(:,:) = energy*cmplx(s_00(:,:),kind=dp) - cmplx(h_00(:,:),kind=dp)
    tmp3(:,:) = cmplx(h_01(:,:),kind=dp) - energy*cmplx(s_01(:,:),kind=dp)


    ! invert k00 --> tmp1
    tmp1 = cmplx_0
    do iorb=1,norb
       tmp1(iorb,iorb) = cmplx_1
    end do


    call ZGESV(norb,norb,tmp2,norb,ipiv,tmp1,norb,info)

    ! exit early if inversion fails
    if (info /= 0) then
      write(stdout,'(a,i6,a)') '     Warning in compute_transfer(1): &
            &zgesv returned info=',info, '. Ignoring this energy.'
      ! make sure deallocation is done
      goto 200
    endif


    ! compute the t0 and t0_bar matrices
    t0 = cmplx_0
    t0_bar = cmplx_0


    call linalg_mat_mul_serial(t0,tmp1,tmp3,opB='C')
    call linalg_mat_mul_serial(t0_bar,tmp1,tmp3)


    ! Initialization of the bulk transfer matrices
    transf(:,:) = t0(:,:)
    tprod(:,:) = t0_bar(:,:)
    transf_bar(:,:) = t0_bar(:,:)
    tprod_bar(:,:) = t0(:,:)

    !================================================================!
    !=== Main iteration loop
    !===   transf(i) = transf(i-1) + tprod(i-1)*t(i)
    !===   transf_bar(i) = transf_bar(i-1) + tprod_bar(i-1)*t(i)
    !================================================================!

    do iter=1,maxiter

       !======================================!
       !=== Compute t1=t(i+1) from t0=t(i) ===!
       !======================================!
       call linalg_mat_mul_serial(tmp1,t0,t0_bar)
       call linalg_mat_mul_serial(tmp2,t0_bar,t0)

       tmp3(:,:) = -tmp1(:,:)-tmp2(:,:)
       do iorb=1,norb
          tmp3(iorb,iorb) = cmplx_1 + tmp3(iorb,iorb)
       end do

       tmp1 = cmplx_0
       do iorb=1,norb
          tmp1(iorb,iorb)=cmplx_1
       end do

       call ZGESV(norb,norb,tmp3,norb,ipiv,tmp1,norb,info)
       ! exit early if inversion fails
       if (info /= 0) then
          write(stdout,'(a,i6,a)') '     Warning in compute_transfer(2): &
               &zgesv returned info=',info, '. Ignoring this energy.'
          ! make sure deallocation is done
          goto 200
       endif

       call linalg_mat_mul_serial(tmp2,t0,t0)
       call linalg_mat_mul_serial(tmp3,t0_bar,t0_bar)

       call linalg_mat_mul_serial(t1,tmp1,tmp2)
       call linalg_mat_mul_serial(t1_bar,tmp1,tmp3)

       !===============================================!
       !=== Compute transf(i+1) and transf_bar(i+1) ===!
       !===============================================!
       call linalg_mat_mul_serial(tmp1,tprod,t1)
       call linalg_mat_mul_serial(tmp2,tprod,t1_bar)

       call ZAXPY(norb2,cmplx_1,tmp1,1,transf,1)
       tprod = tmp2

       call linalg_mat_mul_serial(tmp1,tprod_bar,t1_bar)
       call linalg_mat_mul_serial(tmp2,tprod_bar,t1)

       call ZAXPY(norb2,cmplx_1,tmp1,1,transf_bar,1)
       tprod_bar = tmp2

       !=================!
       !=== Update t0 ===!
       !=================!
       t0 = t1
       t0_bar = t1_bar

       !=============================!
       !=== Check the convergence ===!
       !=============================!
       conver = 0.0_dp
       conver2 = 0.0_dp

       do jorb=1,norb
          do iorb=1,norb
             conver=conver+sqrt(real(t1(iorb,jorb),dp)**2+aimag(t1(iorb,jorb))**2)
             conver2=conver2+sqrt(real(t1_bar(iorb,jorb),dp)**2+aimag(t1_bar(iorb,jorb))**2)
          end do
       end do

       ! check for nans / divergences of conver
       if(conver > max_conv .or. conver2 > max_conv .or. &
          utils_isnan(conver) .or. utils_isnan(conver2)) then
          write(stdout,'(a)') &
              '     Warning: Failed to compute transfer matrix. &
              &Ignoring this energy.'
          info = -9999

          ! make sure deallocation is done
          goto 200
       endif


       if (conver.lt.convtol .and. conver2.lt.convtol) exit

       ! print warning if maximum iterations reached
       if (iter == maxiter) then
          write(stdout,'(a,2e12.3)') '     Warning in etrans_compute_transfer: &
               &maximum iterations reached. Will continue. Residual = ', &
                conver, conver2
       endif

    end do

200 continue ! goto flag if error detected

    ! tidy up workspace
    deallocate(tprod,stat=ierr)
    call utils_dealloc_check('etrans_compute_transfer','tprod',ierr)
    deallocate(tprod_bar,stat=ierr)
    call utils_dealloc_check('etrans_compute_transfer','tprod_bar',ierr)
    deallocate(t0,stat=ierr)
    call utils_dealloc_check('etrans_compute_transfer','t0',ierr)
    deallocate(t0_bar,stat=ierr)
    call utils_dealloc_check('etrans_compute_transfer','t0_bar',ierr)
    deallocate(t1,stat=ierr)
    call utils_dealloc_check('etrans_compute_transfer','t1',ierr)
    deallocate(t1_bar,stat=ierr)
    call utils_dealloc_check('etrans_compute_transfer','t1_bar',ierr)
    deallocate(tmp1,stat=ierr)
    call utils_dealloc_check('etrans_compute_transfer','tmp1',ierr)
    deallocate(tmp2,stat=ierr)
    call utils_dealloc_check('etrans_compute_transfer','tmp2',ierr)
    deallocate(tmp3,stat=ierr)
    call utils_dealloc_check('etrans_compute_transfer','tmp3',ierr)

    return

  end subroutine etrans_compute_transfer

  !====================================================================!
  !====================================================================!

  subroutine etrans_calculate_lead_potential(lead,lead_hsm,il,cell)
    !================================================================!
    !                                                                !
    ! Approximate the lead potential by calculating the chemical     !
    ! potential, in the same manner as EDFT smearing calculation.    !
    !                                                                !
    ! Rounds non-integer electronic occupancies to an integer to     !
    ! to avoid large jumps in the chemical potential in systems with !
    ! a band gap.                                                    !
    !                                                                !
    ! The calculation is trivially parallised over kpoints.          !
    !----------------------------------------------------------------!
    ! Written by Robert Bell, June 2013                              !
    !================================================================!

    use comms, only: pub_on_root, comms_reduce, pub_total_num_procs, &
         pub_my_proc_id
    use constants, only: DP, two_pi, stdout
    use ensemble_dft, only: edft_fermi_level
    use rundat, only: pub_rootname, pub_etrans_lead_nkpoints, &
         pub_edft_smearing_width, pub_num_spins
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_heapsort

    implicit none

    ! Arguments
    type(block_def), intent(inout) :: lead
    type(hsm_type), intent(inout)  :: lead_hsm
    integer, intent(in)            :: il ! lead index
    type(CELL_INFO), intent(in)   :: cell

    ! LAPACK subroutine
    external :: zhegv

    ! internal
    integer :: nkpoints
    integer :: nbands


    real(kind=DP), allocatable :: bands(:,:,:)
    complex(kind=DP), allocatable :: hk(:,:), sk(:,:) ! h/s with bloch phases

    ! local variables
    integer :: itemp, ik, is, ierr, lwork, n_occ(pub_num_spins)
    real(kind=DP) :: integrated_Ne(1:pub_num_spins)
    complex(kind=DP), allocatable :: zwork(:)
    real(kind=DP), allocatable :: rwork(:,:)
    real(kind=DP), allocatable :: kpoints(:)
    real(kind=DP), allocatable :: weights(:)
    real(kind=DP), allocatable :: occ(:,:,:)
    real(kind=DP) :: k
    real(kind=DP) :: mu(1:pub_num_spins)  ! chemical potential
    character(len=128) :: filename
    logical :: sk_is_not_pos_def
    real(kind=DP) :: nelec

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering etrans_calculate_lead_potential'

    call timer_clock('etrans_calculate_lead_potential',1)

    nkpoints = pub_etrans_lead_nkpoints

    allocate(kpoints(0:nkpoints-1),stat=ierr)
    call utils_alloc_check('etrans_calculate_lead_potential','kpoints',ierr)
    allocate(weights(0:nkpoints-1),stat=ierr)
    call utils_alloc_check('etrans_calculate_lead_potential','weights',ierr)

    ! special case
    if (nkpoints.eq.1) then
       kpoints(:) = 0.0_DP
       weights(:) = 1.0_DP
    else
       do ik=0,nkpoints-1
          kpoints(ik) = 0.5_DP*real(ik,kind=DP)/(nkpoints-1)
          ! if k/=0 or pi/a, then use twice the weight from symmetry
          if (ik == 0 .or. ik == nkpoints-1) then
             weights(ik) = 1.0_DP
          else
             weights(ik) = 2.0_DP
          endif
       enddo
    endif
    weights(:) = weights(:) / sum(weights(:))

    ! initialise success tag
    lead%have_lead_mu = .true.

    nbands = lead_hsm%norb / lead%num_unit_cells

    ! allocate workspaces
    allocate(hk(nbands,nbands),stat=ierr)
    call utils_alloc_check('etrans_calculate_lead_potential','hk',ierr)
    allocate(sk(nbands,nbands),stat=ierr)
    call utils_alloc_check('etrans_calculate_lead_potential','sk',ierr)
    allocate(zwork(2*nbands),stat=ierr)
    call utils_alloc_check('etrans_calculate_lead_potential','zwork',ierr)
    allocate(rwork(nbands,3),stat=ierr)
    call utils_alloc_check('etrans_calculate_lead_potential','rwork',ierr)
    allocate(bands(nbands,0:nkpoints-1,pub_num_spins),stat=ierr)
    call utils_alloc_check('etrans_calculate_lead_potential','bands',ierr)
    allocate(occ(nbands,0:nkpoints,1:pub_num_spins),stat=ierr)
    call utils_alloc_check('etrans_calculate_lead_potential','occ',ierr)

    bands = 0.0_DP

    ! find workspace
    lwork = -1
    call zhegv(1,'N','L',nbands,hk,nbands,sk,nbands,bands(1,0,1),zwork,lwork,rwork(1,1),ierr)
    if (ierr /= 0) then
       lwork = 2*nbands
    else
       lwork = nint(real(zwork(1),kind=DP))
       deallocate(zwork,stat=ierr)
       call utils_dealloc_check('etrans_calculate_lead_potential','zwork',ierr)
       allocate(zwork(lwork),stat=ierr)
       call utils_alloc_check('etrans_calculate_lead_potential','zwork',ierr)
    end if

    ! account for spin degeneracy
    do is=1,pub_num_spins

       ! calculate the bands using a tight binding approach
       ! e.g. h = h00 + exp(i*k) h_01 + exp(-i*k) h_01^{dag}  (if lead is only 1 unit cell)
       kpoint_loop: &
       do ik=0,nkpoints-1

          ! trivial parellisation
          if (modulo(ik,pub_total_num_procs) .ne. pub_my_proc_id .or. &
               & pub_my_proc_id > nkpoints) cycle kpoint_loop

          k = two_pi*kpoints(ik)


          ! build hk, sk
          call internal_kpoint_matrix(hk,lead_hsm%h00(:,:,is),lead_hsm%h01(:,:,is), &
               k,lead%num_unit_cells,nbands,lead%forward)
          call internal_kpoint_matrix(sk,lead_hsm%s00(:,:),lead_hsm%s01(:,:), &
               k,lead%num_unit_cells,nbands,lead%forward)

          ! solve for the eigenvalues
          call zhegv(1,'N','L',nbands,hk,nbands,sk,nbands,bands(:,ik,is),zwork,lwork,rwork(1,1),ierr)

          ! make a note if failed due to sk non-positive definite
          if (ierr .gt. nbands) sk_is_not_pos_def = .true.

          ! if this failed, move on to next lead
          if (ierr .eq. 0) then
             lead%have_lead_mu = &
                  lead%have_lead_mu .and. .true.
          else
             lead%have_lead_mu = .false.
          end if
       end do kpoint_loop

       ! synchronise bands data, and diagonalisation success, and pos def
       call comms_reduce('SUM',bands(:,:,is))
       ! jd: I removed 'kind=DP', because it made this an INTEGER(8),
       !     which then got converted to the INTEGER(4) that n_occ is,
       !     giving a compiler warning.
       n_occ(is) = nint(0.5_DP*lead%occ(is)/lead%num_unit_cells)

    end do ! spin

    call comms_reduce('AND',lead%have_lead_mu)
    call comms_reduce('OR',sk_is_not_pos_def)

    if (lead%have_lead_mu) then
       occ = 0.0_DP
       ! kkbd: Modified call because edft_fermi_level can take fixed or free
       !       spin occupancy, and now takes in a non-integer spin channel
       !       occupancy
       nelec = sum(n_occ)
       ! ab: edft_fermi_level calculates the Fermi level for both spin channels
       !     Hence, it is called only once outside the spin loop.
       call edft_fermi_level(occ, mu, integrated_Ne, bands, &
            real(n_occ,kind=dp), nbands, nbands, nkpoints, nelec, &
            weights, pub_edft_smearing_width, -1)

       lead%efermi(:) = mu(:)
    else
       if (pub_on_root) then
          write(stdout,'(1x,a,i4,/,10x,a)') &
               &"WARNING: Unable to determine lead potential for lead ",il,&
               &"Check carefully that the calculated transmission spectrum &
               &is what is expected!"
          if(sk_is_not_pos_def) then
             write(stdout,'(10x,a)') &
                 "Non positive definite k-point overlap matrix detected"
          endif
          write(stdout,'(10x,a,i5,a,i5,/,10x,a)') &
               "ZHEGV returned error code ", ierr, '/', nbands, &
               &"Ignoring this potential."
       endif

       lead%efermi(:) = 0.0_DP

       ! make sure deallocation is done
       goto 200
    endif


    ! reference energy for this lead
    lead_hsm%eref = sum(lead%efermi(:))/pub_num_spins

    ! write out lead .bands file
    if (pub_on_root) then
       write(filename,'(2a,i2.2,a)') trim(adjustl(pub_rootname)),'_lead',il,'.bands'
       call etrans_write_lead_bands(bands,kpoints,weights,lead%efermi,&
            & lead%occ/lead%num_unit_cells,cell,filename)
    endif

200 continue ! (if error in ZHEGV, branch here. If no error, will reach here anyway)

    ! tidy up workspace
    deallocate(hk,stat=ierr)
    call utils_dealloc_check('etrans_calculate_lead_potential','hk',ierr)
    deallocate(sk,stat=ierr)
    call utils_dealloc_check('etrans_calculate_lead_potential','sk',ierr)
    deallocate(zwork,stat=ierr)
    call utils_dealloc_check('etrans_calculate_lead_potential','zwork',ierr)
    deallocate(rwork,stat=ierr)
    call utils_dealloc_check('etrans_calculate_lead_potential','rwork',ierr)
    deallocate(bands,stat=ierr)
    call utils_dealloc_check('etrans_calculate_lead_potential','bands',ierr)
    deallocate(occ,stat=ierr)
    call utils_dealloc_check('etrans_calculate_lead_potential','occ',ierr)

    deallocate(kpoints,stat=ierr)
    call utils_dealloc_check('etrans_calculate_lead_potential','kpoints',ierr)
    deallocate(weights,stat=ierr)
    call utils_dealloc_check('etrans_calculate_lead_potential','weights',ierr)

    call timer_clock('etrans_calculate_lead_potential',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving etrans_calculate_lead_potential'

    contains

    !==================================!
    !==================================!

    subroutine internal_kpoint_matrix(mk,m00,m01,k,n_unit_cells,nbands,is_fwd)
       !========================================================================!
       ! Compute the matrix m at kpoint k by adding relevant blocks of m00, m01 !
       ! with the appropriate Bloch phase.                                      !
       ! E.g. for 1 unit cell:                                                  !
       !   hk = h00 + exp(i*k) h_01 + exp(-i*k) h_01^{dag}                      !
       !------------------------------------------------------------------------!
       ! Written by Robert Bell, Aug 2013                                       !
       !========================================================================!

       use constants, only: DP, cmplx_i

       implicit none

       ! arguments
       integer,          intent(in)  :: nbands
       complex(kind=DP), intent(out) :: mk(nbands,nbands)
       real(kind=DP),    intent(in)  :: m00(:,:), m01(:,:)
       real(kind=DP),    intent(in)  :: k
       integer,          intent(in)  :: n_unit_cells
       logical,          intent(in)  :: is_fwd

       ! internal
       integer :: icell, i1, i2


       ! construct hamiltonian and overlap by adding matrix elements for subsequent
       ! unit cells, with appropriate Bloch phases
       mk = cmplx(m00(1:nbands,1:nbands),kind=DP)

       ! add m00 contributions
       do icell = 1, n_unit_cells-1
          i1 = icell*nbands+1
          i2 = (icell+1)*nbands
          mk = mk + exp(+cmplx_i*k*icell)*cmplx(m00(1:nbands,i1:i2),kind=DP) &
                  + exp(-cmplx_i*k*icell)*cmplx(transpose(m00(1:nbands,i1:i2)),kind=DP)
       enddo

       ! add m01 contributions
       do icell = n_unit_cells, 2*n_unit_cells-1
          i1 = (icell-n_unit_cells)*nbands+1
          i2 = (icell-n_unit_cells+1)*nbands
          if (is_fwd) then
             mk = mk + exp(+cmplx_i*k*icell)*cmplx(m01(1:nbands,i1:i2),kind=DP) &
                     + exp(-cmplx_i*k*icell)*cmplx(transpose(m01(1:nbands,i1:i2)),kind=DP)
          else ! FIXME
             !mk = mk + exp(+cmplx_i*k*icell)*cmplx(transpose(m01(i1:i2,1:nbands)),kind=DP) &
             !        + exp(-cmplx_i*k*icell)*cmplx(m01(i1:i2,1:nbands),kind=DP)
             mk = mk + exp(-cmplx_i*k*icell)*cmplx(m01(i1:i2,1:nbands),kind=DP) &
                     + exp(+cmplx_i*k*icell)*cmplx(transpose(m01(i1:i2,1:nbands)),kind=DP)
          endif ! FIXME
       enddo


     end subroutine internal_kpoint_matrix

  end subroutine etrans_calculate_lead_potential

  !====================================================================!
  !====================================================================!

  subroutine etrans_calculate_lead_ion_occ(leads, elements)
    !========================================================================!
    ! Routine to calculate the ionic charge of the leads.                    !
    !                                                                        !
    !------------------------------------------------------------------------!
    ! Written by Robert Bell, Aug 2013                                       !
    !========================================================================!


    use constants, only: DP
    use ion, only: ELEMENT
    use parallel_strategy, only: par=>pub_par

    implicit none

    ! arguments
    type(block_def), intent(inout) :: leads(:)
    type(ELEMENT), intent(in) :: elements(par%nat)        ! elements

    ! local
    integer :: il
    integer :: iat_start, iat_stop, iat
    real(kind=DP) :: lead_ionic_popn

    ! calculate lead ionic populations
    do il=1,pub_nleads
       iat_start = leads(il)%atms(1)
       iat_stop = leads(il)%atms(2)

       lead_ionic_popn = 0.0_DP
       do iat=iat_start,iat_stop
          lead_ionic_popn = lead_ionic_popn + elements(iat)%ion_charge
       enddo

       ! set lead information
       leads(il)%ionic_popn = lead_ionic_popn
    enddo

  end subroutine etrans_calculate_lead_ion_occ

  !====================================================================!
  !====================================================================!

  subroutine etrans_calculate_lead_elec_occ(leads, ham_type, &
       & overlap, denskern, ngwf_basis)
    !========================================================================!
    ! Routine to calculate the number of electrons associated with the leads.!
    ! Given exactly as a sum over Mulliken charges in the leads.             !
    ! Reuses code from properties_popn_analysis                              !
    !------------------------------------------------------------------------!
    ! Written by Robert Bell, Aug 2013                                       !
    !========================================================================!


    use comms, only: comms_reduce, pub_my_proc_id
    use constants, only: DP
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_num_spins, pub_spin_fac, pub_num_kpoints, &
         PUB_1K
    use sparse, only: sparse_get_block, sparse_get_par
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, &
         sparse_embed_product, sparse_embed_create, sparse_embed_destroy
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! arguments
    type(block_def), intent(inout) :: leads(:)
    character(len=*), intent(in)   :: ham_type
    type(SPAM3_EMBED), intent(in) :: overlap
    type(SPAM3_EMBED_ARRAY), intent(in) :: denskern
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)

    ! local
    type(SPAM3_EMBED) :: ks
    integer :: ierr, il, is
    integer :: iat_start, iat_stop, iat, ingwf, loc_iat
    integer :: max_ngwfs_on_atom
    real(kind=DP), allocatable :: over_block(:,:), q_atom(:,:)
    real(kind=DP) :: block_tr
    type(PARAL_INFO), pointer :: par

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine etrans_calculate_lead_elec_occ not ready yet for more&
         & than one k-point.')

    ! denskern is only in valence basis
    ! if not in this basis, assume elec occ = ionic charge
    if (ham_type /= 'valence') then
       do il=1,pub_nleads
          leads(il)%occ(:) = leads(il)%ionic_popn / pub_num_spins
       enddo

       return
    endif

    ! rc2013: get parallel strategy from SPAM3
    call sparse_get_par(par, overlap%p)

    ! Space for blocks
    max_ngwfs_on_atom = maxval(ngwf_basis(1)%num_on_atom)
    allocate(over_block(max_ngwfs_on_atom,max_ngwfs_on_atom),stat=ierr)
    call utils_alloc_check('etrans_calculate_lead_ionic_occupancies','over_block',ierr)

    ! Generate block diagonal matrix for K.S product
    ks%structure = 'D'
    call sparse_embed_create(ks)

    ! Allocate space for atomic populations
    allocate(q_atom(par%nat,pub_num_spins),stat=ierr)
    call utils_alloc_check('etrans_calculate_lead_ionic_occupancies','q_atom',ierr)
    q_atom = 0.0_DP

    ! Loop over spins
    do is=1,pub_num_spins

       ! Calculate product of density kernel and overlap
       call sparse_embed_product(ks,denskern%m(is,PUB_1K),overlap)

       iat = par%first_atom_on_proc(pub_my_proc_id)
       do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)

          ! Calculate atomic populations on this proc
          call sparse_get_block(over_block,ks%p,iat,iat)

          block_tr = 0.0_DP
          do ingwf=1,ngwf_basis(1)%num_on_atom(iat)
             block_tr = block_tr + over_block(ingwf,ingwf)
          end do
          q_atom(par%orig_atom(iat),is) = block_tr
          iat = iat + 1
       end do

    end do

    ! Sum up over all procs
    call comms_reduce('SUM',q_atom)

    ! sum occupancies for each lead
    do il=1,pub_nleads

       iat_start = leads(il)%atms(1)
       iat_stop = leads(il)%atms(2)

       leads(il)%occ(:) = 0.0_DP
       ! and then atoms in each lead
       do iat=iat_start,iat_stop
          leads(il)%occ(:) = leads(il)%occ(:) + q_atom(iat,:)
       enddo
       leads(il)%occ(:) = leads(il)%occ(:)*pub_spin_fac
    enddo

    ! Deallocate workspace
    deallocate(over_block,stat=ierr)
    call utils_dealloc_check('etrans_calculate_lead_ionic_occupancies', &
         'over_block',ierr)
    deallocate(q_atom,stat=ierr)
    call utils_dealloc_check('etrans_calculate_lead_ionic_occupancies', &
         'q_atom',ierr)
    call sparse_embed_destroy(ks)

  end subroutine etrans_calculate_lead_elec_occ

  !====================================================================!
  !====================================================================!

  subroutine etrans_eref_diagonalisation(ham,rep,ngwf_basis,efermi)

    !========================================================================!
    ! Determines the energy reference as the mid HOMO-LUMO gap of the        !
    ! periodic system. This is strictly *not*
    ! Given exactly as a sum over Mulliken charges in the leads.             !
    ! Reuses code from properties_popn_analysis                              !
    !------------------------------------------------------------------------!
    ! Written by Robert Bell, Aug 2013                                       !
    !========================================================================!

    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
         dense_eigensolve
    use function_basis, only: FUNC_BASIS
    use ngwf_representation, only: NGWF_HAM, NGWF_REP
    use rundat, only: pub_num_spins, pub_num_kpoints, PUB_1K
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! arguments
    type(NGWF_HAM), intent(in)   :: ham
    type(NGWF_REP), intent(in)   :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    real(kind=DP), intent(out)   :: efermi

    ! internal
    type(DEM)     :: sys_eigs          ! dense matrix for eigenvectors
    type(DEM)     :: sys_ham           ! dense matrix for hamiltonian
    type(DEM)     :: sys_overlap       ! dense matrix for overlap
    real(kind=DP), allocatable :: sys_eigen(:,:)     ! hamiltonian eigenvalues

    integer :: num
    integer :: is, ierr

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine etrans_calculate not ready yet for more&
         & than one k-point.')

    num = sum(ngwf_basis(:)%num)
    allocate(sys_eigen(num,pub_num_spins),stat=ierr)
    call utils_alloc_check('etrans_eref_diagonalisation','sys_eigen',ierr)

    call dense_create(sys_overlap,num,num)
    call dense_create(sys_ham,num,num)
    call dense_create(sys_eigs,num,num)
    do is = 1, pub_num_spins
       call dense_convert(sys_overlap,rep%overlap)
       call dense_convert(sys_ham,ham%ham(is))
       call dense_eigensolve(num,sys_eigen(:,is),sys_ham, &
            sys_overlap,1,sys_eigs)
    enddo
    call dense_destroy(sys_eigs)
    call dense_destroy(sys_ham)
    call dense_destroy(sys_overlap)

    efermi = 0.0_DP
    do is = 1, pub_num_spins
       if (rep%n_occ(is,PUB_1K)<ngwf_basis(1)%num.and.rep%n_occ(is,PUB_1K)>0) then
          efermi = efermi + 0.5_DP*(sys_eigen(rep%n_occ(is,PUB_1K)+1,is)+ &
               sys_eigen(rep%n_occ(is,PUB_1K),is))
       else if (rep%n_occ(is,PUB_1K) > 0) then
          efermi = efermi + sys_eigen(rep%n_occ(is,PUB_1K),is) + tiny(1.0_DP)
       else
          efermi = efermi + sys_eigen(1,is) + tiny(1.0_DP)
       end if
    end do
    efermi = efermi / pub_num_spins

    deallocate(sys_eigen,stat=ierr)
    call utils_dealloc_check('etrans_eref_diagonalisation','sys_eigen',ierr)

  end subroutine etrans_eref_diagonalisation

  !====================================================================!
  !====================================================================!

  subroutine etrans_write_hsm(ldef,hsm,filename)

    !========================================================================!
    ! Writes the matrices to unformatted file.                               !
    ! Format used is:                                                        !
    !   block_type                                                           !
    !   ham_type                                                             !
    !   number_of_spins                                                      !
    !   block_start, end, 1st_principle_layer_start, end                     !
    !   matrix_sizes (assumed to be square)                                  !
    !   reference_energy_for_block                                           !
    ! (if lead hsm:)                                                         !
    !   h00(norb,norb,spin)                                                  !
    !   h01(norb,norb,spin)                                                  !
    !   s00(norb,norb)                                                       !
    !   s01(norb,norb)                                                       !
    ! (if lcr hsm:)                                                          !
    !   h_lcr(norb,norb)                                                     !
    !   s_lcr(norb,norb)                                                     !
    !------------------------------------------------------------------------!
    ! Written by Simon Dubois                                                !
    !========================================================================!

    use rundat, only: pub_num_spins
    use utils, only: utils_unit, utils_abort

    implicit none

    ! Arguments
    type(hsm_type), intent(in)   :: hsm
    type(block_def), intent(in)  :: ldef
    character(len=*), intent(in) :: filename

    ! Local variables
    integer    :: io_unit, io, jo, is

    io_unit = utils_unit()
    open(unit=io_unit,file=filename,form='unformatted')

    write(io_unit) ldef%type
    write(io_unit) hsm%ham_type
    write(io_unit) pub_num_spins
    write(io_unit) (ldef%orbs(io),io=1,4)
    write(io_unit) hsm%norb
    write(io_unit) hsm%eref

    if(trim(ldef%type)=='lead') then
       do is = 1, pub_num_spins
          write(io_unit) ((hsm%h00(io,jo,is),io=1,hsm%norb),jo=1,hsm%norb)
       enddo
       do is = 1, pub_num_spins
          write(io_unit) ((hsm%h01(io,jo,is),io=1,hsm%norb),jo=1,hsm%norb)
       enddo

       write(io_unit) ((hsm%s00(io,jo),io=1,hsm%norb),jo=1,hsm%norb)
       write(io_unit) ((hsm%s01(io,jo),io=1,hsm%norb),jo=1,hsm%norb)
    elseif (trim(ldef%type)=='lcr') then
       write(io_unit) ((hsm%h_lcr(io,jo),io=1,hsm%norb),jo=1,hsm%norb)
       write(io_unit) ((hsm%s_lcr(io,jo),io=1,hsm%norb),jo=1,hsm%norb)
    else
       call utils_abort('Error in etrans_write_hsm: unknown hsm type '//&
           trim(ldef%type))
    endif


    close(io_unit)
  end subroutine etrans_write_hsm

  !====================================================================!
  !====================================================================!

  subroutine etrans_read_lead_hsm(lead_def,lead_hsm,ilead)

    !========================================================================!
    ! Reads the lead Hamiltonian/overlap from disk, and ensures that the     !
    ! chemical potential of the new Hamiltonian is the same as that of the   !
    ! old Hamiltonian.                                                       !
    !                                                                        !
    !------------------------------------------------------------------------!
    ! Written by Robert Bell, using code by Simon Dubois                     !
    !========================================================================!

    use constants, only: DP, stdout, HARTREE_IN_EVS
    use comms, only: pub_on_root
    use rundat, only: pub_num_spins
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, &
         utils_assert

    implicit none

    ! Argument
    type(block_def), intent(in)   :: lead_def
    type(hsm_type), intent(inout) :: lead_hsm
    integer, intent(in)           :: ilead

    ! internal
    type(hsm_type) :: buffer_hsm

    integer :: norb, is, ierr
    real(kind=DP) :: defermi

    if (pub_on_root) &
         & write(stdout,'(1x,3a)') "Reading lead H/S matrices from file '", &
         & trim(lead_def%iofile), "'"

    ! do nothing if no file to read
    if (.not. lead_def%ioread) &
         &call utils_abort('Error in etrans_read_lead_hsm: ioread is false')

    ! abort if could not determine lead potential: cannot align
    if (.not. lead_def%have_lead_mu) then
       call utils_abort('Error in etrans_read_lead_hsm: can not read lead &
            &from disk and align chemical potentials if the lead chemical &
            &potential cannot be determined')
    endif

    call utils_assert(trim(lead_def%type) == 'lead','Error in etrans_read_&
         &lead_hsm: block type must be "lead"')

    norb = lead_hsm%norb

    ! allocate buffer_hsm
    allocate(buffer_hsm%h00(norb,norb,pub_num_spins),stat=ierr)
    call utils_alloc_check('etrans_read_leads_hsm','buffer_hsm%h00',ierr)
    allocate(buffer_hsm%h01(norb,norb,pub_num_spins),stat=ierr)
    call utils_alloc_check('etrans_read_leads_hsm','buffer_hsm%h01',ierr)
    allocate(buffer_hsm%s00(norb,norb),stat=ierr)
    call utils_alloc_check('etrans_read_leads_hsm','buffer_hsm%s00',ierr)
    allocate(buffer_hsm%s01(norb,norb),stat=ierr)
    call utils_alloc_check('etrans_read_leads_hsm','buffer_hsm%s01',ierr)

    ! initialise buffer_hsm data
    buffer_hsm%ham_type = lead_hsm%ham_type
    buffer_hsm%norb = lead_hsm%norb


    ! read in file to buffer
    call etrans_read_hsm(lead_def,buffer_hsm,lead_def%iofile)

    ! copy buffer_hsm into lead_hsm, correcting reference energies if needed
    defermi = lead_hsm%eref - buffer_hsm%eref
    do is = 1, pub_num_spins
       lead_hsm%h00(:,:,is) = buffer_hsm%h00(:,:,is) + &
            defermi * buffer_hsm%s00(:,:)
       lead_hsm%h01(:,:,is) = buffer_hsm%h01(:,:,is) + &
            defermi * buffer_hsm%s01(:,:)
    enddo
    lead_hsm%s00(:,:) = buffer_hsm%s00(:,:)
    lead_hsm%s01(:,:) = buffer_hsm%s01(:,:)

    if (pub_on_root .and. abs(defermi) > tiny(1.0_DP)) then
       write(stdout,'(4x,a,i2.2,a,f10.5,a)') 'Aligned new lead (',ilead, &
            &') to old lead energy: dE_fermi = ', defermi*HARTREE_IN_EVS,' eV'
    endif


    ! deallocate
    deallocate(buffer_hsm%h00,stat=ierr)
    call utils_dealloc_check('etrans_read_leads_hsm','buffer_hsm%h00',ierr)
    deallocate(buffer_hsm%h01,stat=ierr)
    call utils_dealloc_check('etrans_read_leads_hsm','buffer_hsm%h01',ierr)
    deallocate(buffer_hsm%s00,stat=ierr)
    call utils_dealloc_check('etrans_read_leads_hsm','buffer_hsm%s00',ierr)
    deallocate(buffer_hsm%s01,stat=ierr)
    call utils_dealloc_check('etrans_read_leads_hsm','buffer_hsm%s01',ierr)

  end subroutine etrans_read_lead_hsm

  !====================================================================!
  !====================================================================!

  subroutine etrans_read_hsm(ldef,hsm,filename)

    use comms, only: pub_on_root, comms_bcast, comms_barrier, &
         pub_root_proc_id
    use constants, only: stdout
    use rundat, only: pub_num_spins
    use utils, only: utils_abort, utils_unit, utils_close_unit_check

    implicit none

    ! Arguments
    type(hsm_type), intent(inout)   :: hsm
    type(block_def), intent(in)     :: ldef
    character(len=*), intent(in)    :: filename

    ! Local variables
    character(len=4)  :: block_type
    character(len=20) :: ham_type
    integer    :: nspin, norb, orbs(4)
    integer    :: io, jo, is
    integer    :: io_unit, io_stat
    logical    :: reverse

    if (pub_on_root) then

       ! Find available unit specifier
       io_unit = utils_unit()

       open(unit=io_unit,iostat=io_stat,file=trim(filename),&
            form='unformatted',action='read')

       if (io_stat /= 0) then
          call utils_abort('ERROR : etrans_read_hsm failed &
               &to open '//trim(filename)//'. ', io_stat, io_unit)
       endif

       ! Nspin, orbs, norb
       read(io_unit) block_type
       read(io_unit) ham_type
       read(io_unit) nspin
       read(io_unit) (orbs(io),io=1,4)
       read(io_unit) norb

       if (trim(block_type) .ne. trim(ldef%type)) then
          call utils_abort('Error in etrans_read_hsm reading '//&
               &trim(filename)//': incompatible hsm block type. Found "'//&
               &trim(block_type)//'" expecting "'//trim(ldef%type)//'"')
       endif

       if (trim(ham_type) .ne. trim(hsm%ham_type)) then
          call utils_abort('Error in etrans_read_hsm reading '//&
               &trim(filename)//': incompatible hamiltonian type. &
               &Found "'//ham_type//'" expecting "'//trim(hsm%ham_type)//'"')
       endif

       if (nspin .ne. pub_num_spins) then
          ! Do not error if reading LCR .hsm: it only contains one spin
          ! even if the calculation is spin polarised, so the Hamiltonain
          ! can be reused (but print a warning)
          if (block_type == 'lcr') then
              write(stdout,'(2(a,i2),a)') 'Warning in etrans_read_hsm: matrix &
                  &elements in file'//trim(filename)//' were calculated &
                  &with ', nspin, ' spin channels, but this calculation &
                  &uses ', pub_num_spins, '. Will continue.'
          else
              call utils_abort('Error in etrans_read_hsm reading '//&
                  &trim(filename)//': incompatible number of spin &
                  &components:',nspin,pub_num_spins)
          endif
       endif
       if (norb .ne. hsm%norb) then
          call utils_abort('Error in etrans_read_hsm reading '//&
               &trim(filename)//': incompatible number of NGWFs:', &
                norb, hsm%norb)
       endif


       if (block_type=='lead') then
          if ((orbs(3).ge.orbs(1) .and. ldef%forward) .or. &
               (orbs(3).lt.orbs(1) .and. .not.ldef%forward)) then
             reverse = .false.
          else
             reverse = .true.
          endif
       endif

       read(io_unit) hsm%eref

       if (trim(block_type)=='lead') then
          do is = 1, pub_num_spins
             read(io_unit) ((hsm%h00(io,jo,is),io=1,norb),jo=1,norb)
          enddo
          do is = 1, pub_num_spins
             if (reverse) then
                read(io_unit) ((hsm%h01(jo,io,is),io=1,norb),jo=1,norb)
             else
                read(io_unit) ((hsm%h01(io,jo,is),io=1,norb),jo=1,norb)
             endif
          enddo

          read(io_unit) ((hsm%s00(io,jo),io=1,norb),jo=1,norb)
          if (reverse) then
             read(io_unit) ((hsm%s01(jo,io),io=1,norb),jo=1,norb)
          else
             read(io_unit) ((hsm%s01(io,jo),io=1,norb),jo=1,norb)
          endif
       elseif (trim(block_type)=='lcr') then
          read(io_unit) ((hsm%h_lcr(io,jo),io=1,norb),jo=1,norb)
          read(io_unit) ((hsm%s_lcr(io,jo),io=1,norb),jo=1,norb)
       else
          call utils_abort('Error in etrans_read_hsm: unknown block_type&
               & reading "'//trim(filename)//'": '//trim(block_type))
       endif


       close(unit=io_unit,iostat=io_stat)
       call utils_close_unit_check('etrans_read_hsm',filename,io_stat)

    endif

    ! broadcast data
    call comms_bcast(pub_root_proc_id,hsm%eref)

    if (trim(ldef%type)=='lead') then
       call comms_bcast(pub_root_proc_id,hsm%h00(:,:,:))
       call comms_bcast(pub_root_proc_id,hsm%h01(:,:,:))
       call comms_bcast(pub_root_proc_id,hsm%s00(:,:))
       call comms_bcast(pub_root_proc_id,hsm%s01(:,:))
    elseif (trim(ldef%type)=='lcr') then
       call comms_bcast(pub_root_proc_id,hsm%h_lcr(:,:))
       call comms_bcast(pub_root_proc_id,hsm%s_lcr(:,:))
    else
       call utils_abort('Error in etrans_read_hams: unknown block_type&
            & broadcasting "'//trim(filename)//'": '//trim(block_type))
    endif

    call comms_barrier()

  end subroutine etrans_read_hsm

  !==========================================================================!
  !==========================================================================!

  subroutine etrans_write_lead_bands(bands,kpoints,weights,efermi,n_occ,cell,filename)

    use rundat, only: pub_num_spins
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_unit, utils_abort

    implicit none

    ! arguments
    real(kind=DP), intent(in) :: bands(:,:,:)
    real(kind=DP), intent(in) :: kpoints(:)
    real(kind=DP), intent(in) :: weights(:)
    real(kind=DP), intent(in) :: efermi(pub_num_spins)
    real(kind=DP), intent(in) :: n_occ(pub_num_spins)
    type(CELL_INFO), intent(in) :: cell
    character(len=128), intent(in) :: filename

    ! local variables
    integer :: nbands, nkpts
    integer :: bands_unit, ierr, ikpt, is, iband


    nbands = size(bands,dim=1)
    nkpts = size(bands,dim=2)
    if (nkpts .ne. size(kpoints)) then
       call utils_abort('Error in etrans_write_lead_bands: inconsistent number of kpoints')
    endif
    if (pub_num_spins .ne. size(bands,dim=3)) then
       call utils_abort('Error in etrans_write_lead_bands: inconsistent &
            &number of spin components')
    endif

    ! Open output file
    bands_unit = utils_unit()
    open(bands_unit,file=trim(filename),iostat=ierr)
    if (ierr /= 0) then
       call utils_abort('Error in etrans_write_lead_bands: opening "'//&
            trim(filename)//'" failed with code ',ierr)
    end if

    write(bands_unit,'(a,i6)') 'Number of k-points',nkpts
    write(bands_unit,'(a,i2)') 'Number of spin components',pub_num_spins
    write(bands_unit,'(a,2f8.3)') 'Number of electrons ',n_occ(1:pub_num_spins)
    if (pub_num_spins == 1) then
       write(bands_unit,'(a,i8)') 'Number of eigenvalues ',nbands
       write(bands_unit,'(a,f12.6)') 'Fermi energy (in atomic units) ', &
            efermi(1)
    else
       write(bands_unit,'(a,2i8)') 'Number of eigenvalues ',nbands, nbands
       write(bands_unit,'(a,2f12.6)') 'Fermi energies (in atomic units) ', &
            efermi(1),efermi(2)
    end if
    write(bands_unit,'(a)') 'Unit cell vectors'
    write(bands_unit,'(3f12.6)') cell%a1%x,cell%a1%y,cell%a1%z
    write(bands_unit,'(3f12.6)') cell%a2%x,cell%a2%y,cell%a2%z
    write(bands_unit,'(3f12.6)') cell%a3%x,cell%a3%y,cell%a3%z

    do ikpt=1,nkpts
       write(bands_unit,'(a,i6,4f12.8)') 'K-point',ikpt,0.0,0.0, &
            kpoints(ikpt),weights(ikpt)
       do is=1,pub_num_spins
          write(bands_unit,'(a,i2)') 'Spin component',is
          do iband=1,nbands
             write(bands_unit,'(f14.8)') bands(iband,ikpt,is)
          end do
       end do
    end do

    ! Close output file
    close(bands_unit,iostat=ierr)
    if (ierr /= 0) then
       call utils_abort('Error in etrans_write_lead_bands: closing "'//&
            trim(filename)//'" failed with code ',ierr)
    end if

  end subroutine etrans_write_lead_bands


  !==========================================================================!
  !==========================================================================!

  subroutine etrans_write_trc(trc_type,ham_type,energy,trc,gf_err)

     use constants, only: DP, HARTREE_IN_EVS, stdout
     use rundat, only: pub_rootname, pub_num_spins
     use utils, only: utils_unit, utils_abort

     implicit none

       ! arguments
       character(len=*), intent(in) :: trc_type
       character(len=*), intent(in) :: ham_type
       complex(kind=DP), intent(in) :: energy(:)
       real(kind=DP), intent(in)    :: trc(:,:,:)
       integer, intent(in)          :: gf_err(:,:,:)

       ! internal
       integer :: ienergy, il, jl, nep, max_l, is, output_unit
       character(len=3)  :: spin_name
       character(len=7)  :: ham_name
       character(len=120):: filename

       nep = size(energy)
       max_l = size(trc,dim=1)

       if (trim(ham_type) == 'joint') then
          ham_name = '_joint'
       else
          ham_name = ''
       endif

       ! check for errors in inversion
       if (any(gf_err .ne. 0)) then
          write(stdout,'(3a)') &
               ' WARNING -- ', trim(trc_type), ' failed at energies (eV)'
          do ienergy=1,nep
             if (any(gf_err(:,:,ienergy).ne.0)) &
                   write(stdout,'(f16.5,a)') real(energy(ienergy))*HARTREE_IN_EVS
          enddo
          write(stdout,'(a/)') &
               ' Try increasing etrans_ecmplx if this produces bad transmission'
       endif

       do is = 1, pub_num_spins

          if (pub_num_spins == 1) then
             spin_name = ""
          else if (is == 1) then
             spin_name = "_UP"
          else
             spin_name = "_DN"
          endif

          filename = trim(pub_rootname)//'_'//trim(trc_type)//trim(spin_name)// &
               trim(ham_name)//'.TRC'

          output_unit = utils_unit()
          open(unit=output_unit, form="formatted", file=trim(filename), &
               action="write")

          ! write header
          write(output_unit,'(a,11x,a,3x)', advance='no') '#','E (eV)'
          if (trc_type == 'BULK') then
             do il = 1, pub_nleads
                write(output_unit,'(a2,7x,a4,i2,3x)',advance='no')  '| ','lead', il
             enddo
             write(output_unit,*)
          elseif (trc_type == 'LCR') then
             do il = 1, pub_nleads
                do jl = il+1, pub_nleads
                   write(output_unit,'(a,1x,i2,a,i2,a,3x)',advance='no') &
                        '| T(', il, " --> ", jl,')'
                enddo
             enddo
             write(output_unit,*)
          else
             call utils_abort('Error in etrans_write_trc: unknown trc type')
          endif

          ! write data
          do ienergy = 1, nep
             if (all(gf_err(:,is,ienergy).eq.0)) then
                write(output_unit,'(f18.8)',advance='no') &
                     & real(energy(ienergy))*HARTREE_IN_EVS
                do il = 1, max_l
                   write(output_unit,'(e18.8)',advance='no') trc(il,is,ienergy)
                enddo
                write(output_unit,*)
             endif
          enddo
          close(output_unit)

       enddo ! spins

  end subroutine etrans_write_trc

  !==========================================================================!
  !==========================================================================!

  subroutine etrans_write_trc_eigchan(trc_type,ham_type,energy,trc_eigchan,gf_err)

     use constants, only: DP, HARTREE_IN_EVS
     use rundat, only: pub_rootname, pub_etrans_num_eigchan, pub_num_spins
     use utils, only: utils_unit, utils_abort

     implicit none

     ! arguments
     character(len=*), intent(in) :: trc_type
     character(len=*), intent(in) :: ham_type
     complex(kind=DP), intent(in) :: energy(:)
     real(kind=DP), intent(in)    :: trc_eigchan(:,:,:,:)
     integer, intent(in)          :: gf_err(:,:,:)

     ! internal
     integer :: ienergy, il, nep, is, output_unit, ichan
     character(len=3)  :: spin_name
     character(len=7)  :: ham_name
     character(len=12) :: lead_name
     character(len=120):: filename

     nep = size(energy)

     if (trim(ham_type) == 'joint') then
        ham_name = '_joint'
     else
        ham_name = ''
     endif

     do il = 1, pub_nleads
       write(lead_name,"('_lead',i2.2)") il

       do is = 1, pub_num_spins

          if (pub_num_spins == 1) then
             spin_name = ""
          else if (is == 1) then
             spin_name = "_UP"
          else
             spin_name = "_DN"
          endif

          filename = trim(pub_rootname)//'_'//trim(trc_type)//trim(lead_name)//&
               trim(spin_name)//trim(ham_name)//'.TRC'

          output_unit = utils_unit()
          open(unit=output_unit, form="formatted", file=trim(filename), &
               action="write")

          ! write header
          write(output_unit,'(a,11x,a,3x)', advance='no') &
               '#','E (eV)'
          if (trc_type == 'LCR_channels') then
             do ichan = 1, pub_etrans_num_eigchan
                write(output_unit,'(a2,4x,a7,i2,3x)',advance='no')  &
                     '| ','channel', ichan
             enddo
             write(output_unit,*)
          else
             call utils_abort('Error in etrans_write_trc_eigchan: unknown trc type')
          endif

          ! write data
          do ienergy = 1, nep
             if (all(gf_err(:,is,ienergy).eq.0)) then
                write(output_unit,'(f18.8)',advance='no') &
                     & real(energy(ienergy))*HARTREE_IN_EVS
                do ichan = 1, pub_etrans_num_eigchan
                   write(output_unit,'(e18.8)',advance='no') trc_eigchan(ichan,il,is,ienergy)
                enddo
                write(output_unit,*)
             endif
          enddo
          close(output_unit)

       enddo ! spins
    enddo ! il

  end subroutine etrans_write_trc_eigchan

  !==========================================================================!
  !==========================================================================!

  subroutine etrans_write_dos(dos_type,ham_type,energy,DOS,gf_err)

     use constants, only: DP, HARTREE_IN_EVS
     use rundat, only: pub_rootname, pub_num_spins
     use utils, only: utils_unit, utils_abort

     implicit none

       ! arguments
       character(len=*), intent(in) :: dos_type
       character(len=*), intent(in) :: ham_type
       complex(kind=DP), intent(in) :: energy(:)
       real(kind=DP), intent(in)    :: DOS(:,:,:)
       integer, intent(in)          :: gf_err(:,:,:)

       ! internal
       integer :: ienergy, il, nep, max_l, is, output_unit
       character(len=3)  :: spin_name
       character(len=7)  :: ham_name
       character(len=120):: filename

       nep = size(energy)
       max_l = size(DOS,dim=1)

       if (trim(ham_type) == 'joint') then
          ham_name = '_joint'
       else
          ham_name = ''
       endif


       do is = 1, pub_num_spins

          if (pub_num_spins == 1) then
             spin_name = ""
          else if (is == 1) then
             spin_name = "_UP"
          else
             spin_name = "_DN"
          endif

          filename = trim(pub_rootname)//'_'//trim(dos_type)//trim(spin_name)// &
               trim(ham_name)//'.DOS'

          output_unit = utils_unit()
          open(unit=output_unit, form="formatted", file=trim(filename), &
               action="write")

          ! write header
          if (dos_type == 'BULK') then
             write(output_unit,'(a,11x,a,3x)', advance='no') '#','E (eV)'
             do il = 1, pub_nleads
                write(output_unit,'(a2,7x,a4,i2,3x)',advance='no')  '| ','lead', il
             enddo
             write(output_unit,*)
          elseif (dos_type == 'LCR') then
             write(output_unit,'(a,11x,a,7x,a)') '#','E (eV) |','DOS'
          else
             call utils_abort('Error in etrans_write_dos: unknown dos type')
          endif

          ! write data
          do ienergy = 1, nep
             if (all(gf_err(:,is,ienergy).eq.0)) then
                write(output_unit,'(f18.8)',advance='no') &
                     & real(energy(ienergy))*HARTREE_IN_EVS
                do il = 1, max_l
                   write(output_unit,'(e18.8)',advance='no') dos(il,is,ienergy)
                enddo
                write(output_unit,*)
             endif
          enddo
          close(output_unit)

       enddo ! spins

  end subroutine etrans_write_dos

  !==========================================================================!
  !==========================================================================!

  subroutine etrans_compute_eigchan(ham,ham_type,rep,ngwf_basis,lead_hsm,mdl)

    !================================================================!
    ! Subroutine to compute and plot the transmission eigenchannels, !
    ! calculated using the method of Paulsson and Brandbyge Phys.    !
    ! Rev. B 76, 115117 (2007).                                      !
    !                                                                !
    ! This method is cubically scaling and memory intensive, but     !
    ! uses distributed dense matrix algebra to mitigate this.        !
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Robert Bell, April 2014                             !
    !================================================================!

    ! system matrices refer to the entire system containing all the atoms,
    ! lcr refers only to atoms in the device

    use augmentation, only: augmentation_overlap
    use comms, only: pub_on_root
    use constants, only: DP, stdout, HARTREE_IN_EVS, UP, DN
    use dense, only: DEM, dense_create, dense_destroy
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_HAM, NGWF_REP
    use rundat, only: pub_etrans_num_eigchan, pub_etrans_ecmplx, pub_aug, &
        pub_print_qc, pub_num_spins
    use sparse, only: SPAM3, sparse_create, sparse_destroy
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert
    ! agrecocmplx
    use datatypes, only: FUNCTIONS, data_functions_alloc, data_functions_dealloc


    implicit none

    ! Arguments
    type(NGWF_HAM), intent(in)   :: ham
    type(NGWF_REP), intent(in)   :: rep
    type(hsm_type), intent(inout):: lead_hsm(pub_nleads)
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)
    type(MODEL), intent(in)    :: mdl
    character(len=*), intent(in) :: ham_type

    ! Local Variables
    type(DEM) :: h_dem, s_dem, gf_dem
    real(kind=DP), allocatable :: t_chan(:,:)
    complex(kind=DP), allocatable :: eigchan_coeff(:,:,:)
    real(kind=DP) :: total_trans
    ! agrecocmplx: use FUNCTIONS type to switch between
    ! real and complex case
    !real(kind=DP), allocatable :: eigchan_coeff_sys(:)
    type(FUNCTIONS) :: eigchan_coeff_sys
    real(kind=DP) :: total_dos
    real(kind=DP), allocatable :: partial_dos(:,:) ! Eigenchannel DOS contribution
    complex(kind=DP)  :: energy
    type(DEM)   :: eigs_sys, buffer_sys
    type(SPAM3) :: mo_kern(1)    ! density kernel of one molecular orbital
    type(SPAM3)  :: aug_overlap  ! Augmentation part of overlap matrix
    integer :: ichan, ilead, is, ien
    integer :: iat, iat_distr, num_ngwfs
    integer :: istart, iend, jstart, jend
    integer :: norb_sys, norb_lcr, num_energies
    integer :: ierr, info
    character(len=256) :: scalar_name, file_header

    ! agrecocmplx
    call utils_assert(.not.rep%ngwfs_on_grid(1)%iscmplx, 'Error in&
         & etrans_compute_eigchan: not ready for complex NGWFs yet.')

    norb_sys = sum(ngwf_basis(:)%num)
    norb_lcr = pub_lcr%orbs(2)-pub_lcr%orbs(1)+1
    num_energies = size(pub_etrans_eigchan_en)

    ! allocate workspace
    ! DEM matrices
    call dense_create(s_dem,norb_lcr,norb_lcr)
    call dense_create(h_dem,norb_lcr,norb_lcr)
    call dense_create(gf_dem,norb_lcr,norb_lcr,iscmplx=.true.)

    ! coefficient storage
    allocate(eigchan_coeff(norb_lcr,pub_etrans_num_eigchan,pub_nleads),stat=ierr)
    call utils_alloc_check('etrans_compute_eigchan','eigchan_coeff',ierr)
    allocate(t_chan(pub_etrans_num_eigchan,pub_nleads),stat=ierr)
    call utils_alloc_check('etrans_compute_eigchan','t_chan',ierr)
    allocate(partial_dos(pub_etrans_num_eigchan,pub_nleads),stat=ierr)
    call utils_alloc_check('etrans_compute_eigchan','partial_dos',ierr)
    ! agrecocmplx: allocate using appropriate routine
    !allocate(eigchan_coeff_sys(norb_sys),stat=ierr)
    !call utils_alloc_check('etrans_compute_eigchan','eigchan_coeff_sys',ierr)
    call data_functions_alloc(eigchan_coeff_sys,norb_sys,iscmplx=.false.)

    ! extra workspace required if augmentation charge present
    if (pub_aug) then
       ! dense system matrices
       call dense_create(eigs_sys,norb_sys,norb_sys,iscmplx=.false.)
       call dense_create(buffer_sys,norb_sys,norb_sys,iscmplx=.false.)

       mo_kern(1)%structure = 'K'//rep%postfix
       call sparse_create(mo_kern(1))

       ! ndmh: Create matrix for aug part of overlap
       call sparse_create(aug_overlap,rep%overlap%p)
       call augmentation_overlap(aug_overlap,mdl%pseudo_sp,mdl%paw_sp, &
       rep%sp_overlap%p)
    end if

    do is=1,pub_num_spins! loop over spins

       if (pub_on_root) then
          write(stdout,'(/a)') '====================== &
                         &Computing eigenchannel wavefunctions &
                         &======================'
          if (pub_num_spins == 2 .and. is==UP) then
             write(stdout,'(a)') '====================== &
                         &               UP spin               &
                         &======================'
          else if (pub_num_spins == 2 .and. is==DN) then
             write(stdout,'(a)') '====================== &
                         &               DN spin               &
                         &======================'
          endif
       endif

       ! initialise hsm DEM matrices
       call internal_gather_lcr_hsm_dem

       do ien=1,num_energies ! loop over energies
          energy = cmplx(pub_etrans_eigchan_en(ien),pub_etrans_ecmplx,kind=DP)

          ! compute the Green's function
          call etrans_compute_greenf_dem(gf_dem,h_dem,s_dem,lead_hsm,energy,is,info)
          if (info/=0) cycle

          ! compute the eigenchannels with each lead as the source
          call etrans_eigchan_coeffs_dem(gf_dem,s_dem,lead_hsm,total_trans, &
               t_chan,eigchan_coeff,total_dos,partial_dos)

          ! print out results to a table
          if (pub_on_root) then
             do ilead=1,pub_nleads
               write(stdout,*)
               write(stdout,'(a,i5,a,f11.5)') &
                  '| Source lead ',ilead,',   Energy (eV) = ', &
                  real(energy,kind=DP)*HARTREE_IN_EVS

               write(stdout,'(a)')'| Channel   transmission    Norm'
               write(stdout,'(a)')'|-------------------------------------'
               do ichan=1,pub_etrans_num_eigchan
                  write(stdout,'(a,i2,5x,f14.9,f14.5)') '| ', ichan, &
                        t_chan(ichan,ilead), partial_dos(ichan,ilead)
               enddo
               ! total_dos is less than sum(partial_dos(:,ilead))? TODO
               !write(stdout,'(a,f14.9,f14.5)') '| Total  ',total_trans,total_dos
               !write(stdout,'(a,f14.5)') '| Residual DOS ', &
               !      total_dos - sum(partial_dos(:,ilead))

             enddo

             write(stdout,*)
          endif


          ! plot the eigenchannels
          do ilead=1,pub_nleads ! loop over leads
             do ichan = 1, pub_etrans_num_eigchan ! loop over channels

                if (pub_print_qc) then
                   call internal_print_qc
                else
                   scalar_name = internal_filename(ham_type,ien,ilead,ichan,is)
                   call internal_plot_channel
                endif
             enddo ! channels
          enddo ! leads
       enddo ! energies
       if (pub_on_root) write(stdout,'(a/)')&
            '========================================================&
            &========================='
    enddo ! spin

    ! destroy workspace
    ! agrecocmplx
    !deallocate(eigchan_coeff_sys,stat=ierr)
    !call utils_dealloc_check('etrans_compute_eigchan','eigchan_coeff_sys',ierr)
    call data_functions_dealloc(eigchan_coeff_sys)
    call dense_destroy(h_dem)
    call dense_destroy(s_dem)
    deallocate(eigchan_coeff,stat=ierr)
    call utils_dealloc_check('etrans_compute_eigchan','eigchan_coeff',ierr)
    deallocate(t_chan,stat=ierr)
    call utils_dealloc_check('etrans_compute_eigchan','t_chan',ierr)
    deallocate(partial_dos,stat=ierr)
    call utils_dealloc_check('etrans_compute_eigchan','partial_dos',ierr)


    if (pub_aug) then
       call dense_destroy(eigs_sys)
       call dense_destroy(buffer_sys)
       call sparse_destroy(mo_kern(1))
       call sparse_destroy(aug_overlap)
    endif

    return

    contains

    !============================================!

    subroutine internal_gather_lcr_hsm_dem

      !================================================================!
      ! Gathers the LCR matrix elements into DEM format.               !
      !                                                                !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, April 2014, using code by Simon Dubois !
      !================================================================!

      use dense, only: dense_put_element

      implicit none

      ! internal
      integer :: iat1, iat2, il, jl
      integer :: iorb, jorb, norb
      integer :: iostart,iostop,jostart,jostop

      iat1 = pub_lcr%atms(1)
      iat2 = pub_lcr%atms(2)
      call internal_get_submatrix_dem(iat1,iat2,iat1,iat2, &
           rep%overlap%p,s_dem,norb_lcr,norb_lcr,ngwf_basis)
      call internal_get_submatrix_dem(iat1,iat2,iat1,iat2, &
           ham%ham(is)%p,h_dem,norb_lcr,norb_lcr,ngwf_basis)

      ! Get rid of the periodic boundary conditions
      do il = 1, pub_nleads
         iostart = pub_leads(il)%orbs(1) - pub_lcr%orbs(1) + 1
         iostop  = pub_leads(il)%orbs(2) - pub_lcr%orbs(1) + 1
         do jl = 1, pub_nleads

            if (il == jl) cycle

            jostart = pub_leads(jl)%orbs(1) - pub_lcr%orbs(1) + 1
            jostop  = pub_leads(jl)%orbs(2) - pub_lcr%orbs(1) + 1
             do iorb=iostart,iostop
               do jorb=jostart,jostop
                  call dense_put_element(0.0_DP,s_dem,iorb,jorb)
                  call dense_put_element(0.0_DP,h_dem,iorb,jorb)
               enddo
            enddo
         enddo
      enddo

      ! deposit lead matrices (will fix lead symmetrisation)
      do il=1,pub_nleads
         norb = lead_hsm(il)%norb
         iostart = pub_leads(il)%orbs(1) - pub_lcr%orbs(1)
         jostart = pub_leads(il)%orbs(3) - pub_lcr%orbs(1)
         do iorb=1,norb
           do jorb=1,norb
             call dense_put_element(lead_hsm(il)%s00(iorb,jorb), &
               s_dem,iorb+iostart,jorb+iostart)
             call dense_put_element(lead_hsm(il)%s01(iorb,jorb), &
               s_dem,iorb+iostart,jorb+jostart)
             call dense_put_element(lead_hsm(il)%s01(iorb,jorb), &
               s_dem,jorb+jostart,iorb+iostart)
             call dense_put_element(lead_hsm(il)%h00(iorb,jorb,is), &
               h_dem,iorb+iostart,jorb+iostart)
             call dense_put_element(lead_hsm(il)%h01(iorb,jorb,is), &
               h_dem,iorb+iostart,jorb+jostart)
             call dense_put_element(lead_hsm(il)%h01(iorb,jorb,is), &
               h_dem,jorb+jostart,iorb+iostart)
           enddo
         enddo
      enddo

    end subroutine internal_gather_lcr_hsm_dem

    !============================================!

    subroutine internal_plot_channel

      !================================================================!
      ! Plots the eigenchannel in the NGWF basis. The real and imag    !
      ! parts are printed separately.                                  !
      !                                                                !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, April 2014                             !
      !================================================================!

      use dense, only: dense_product, dense_convert, dense_put_col
      use eigenstates, only: eigenstates_plot_mo
      use parallel_strategy, only: par=>pub_par
      ! agrecocmplx
      use datatypes, only: data_set_to_zero

      implicit none

      file_header = 'Transmission eigenchannel wavefunction'

      ! ---REAL PART---------------------------
      ! deposit device coefficients into system
      ! agrecocmplx
      call data_set_to_zero(eigchan_coeff_sys)
      jstart=1
      do iat = pub_lcr%atms(1),pub_lcr%atms(2)
         iat_distr = par%distr_atom(iat)
         num_ngwfs = ngwf_basis(1)%num_on_atom(iat_distr)
         istart = ngwf_basis(1)%first_on_atom(iat_distr)
         iend = istart + num_ngwfs - 1
         jend = jstart + num_ngwfs - 1

         ! copy data
         ! agrecocmplx
         eigchan_coeff_sys%d(istart:iend) = &
              real(eigchan_coeff(jstart:jend,ichan,ilead),kind=DP)
         jstart = jend + 1
      enddo

      ! ndmh: if augmentation is active, we need to calculate the
      ! ndmh: aug part of the norm, so we need the MO 'kernel'
      if (pub_aug) then
         ! put eigenchannel in first column of a system dense matrix
         call dense_put_col(eigchan_coeff_sys%d,eigs_sys,1)

         ! ndmh: create dense 'kernel' for MO from selected eigenvector
         call dense_product(buffer_sys,eigs_sys,eigs_sys, &
              opA='N',opB='T',first_k=1,last_k=1)
         ! ndmh: convert to SPAM3
         call dense_convert(mo_kern(1),buffer_sys)
      end if

      call eigenstates_plot_mo(eigchan_coeff_sys, rep%ngwfs_on_grid(1), &
           ngwf_basis(1), mdl, file_header, trim(scalar_name)//'_real', &
           aug_overlap, mo_kern(1),print_norm=.false.)


      ! ---IMAG PART---------------------------
      ! deposit device coefficients into system
      ! agrecocmplx
      call data_set_to_zero(eigchan_coeff_sys)
      jstart=1
      do iat = pub_lcr%atms(1),pub_lcr%atms(2)
         iat_distr = par%distr_atom(iat)
         num_ngwfs = ngwf_basis(1)%num_on_atom(par%distr_atom(iat))
         istart = ngwf_basis(1)%first_on_atom(par%distr_atom(iat))
         iend = istart + num_ngwfs - 1
         jend = jstart + num_ngwfs - 1

         ! copy data
         ! agrecocmplx
         eigchan_coeff_sys%d(istart:iend) = &
              aimag(eigchan_coeff(jstart:jend,ichan,ilead))
          jstart = jend + 1
      enddo

      ! ndmh: if augmentation is active, we need to calculate the
      ! ndmh: aug part of the norm, so we need the MO 'kernel'
      if (pub_aug) then
         ! put eigenchannel in first column of a system dense matrix
         ! agrecocmplx
         call dense_put_col(eigchan_coeff_sys%d,eigs_sys,1)

         ! ndmh: create dense 'kernel' for MO from selected eigenvector
         call dense_product(buffer_sys,eigs_sys,eigs_sys, &
              opA='N',opB='T',first_k=1,last_k=1)
         ! ndmh: convert to SPAM3
         call dense_convert(mo_kern(1),buffer_sys)
      end if

      call eigenstates_plot_mo(eigchan_coeff_sys, rep%ngwfs_on_grid(1), &
           ngwf_basis(1), mdl, file_header, trim(scalar_name)//'_imag', &
           aug_overlap, mo_kern(1),print_norm=.false.)

    end subroutine internal_plot_channel

    !============================================!

    subroutine internal_print_qc

      !================================================================!
      ! Prints the eigenchannel coefficients for QC.                   !
      ! Absolute coefficients printed as may come out with different   !
      ! global phase.                                                  !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use utils, only: utils_qc_print

      implicit none

      character(len=22) :: tag
      integer :: iorb
      real(kind=DP) :: coeff
      real(kind=DP) :: t_tol = 0.0000001_DP ! minimum transmission to print

      if (pub_on_root) then
         write(tag,'(a,4(i1.1,a))') 'echan_T(',is,',',ien,',',ilead,',',ichan,')'
         call utils_qc_print(trim(tag),t_chan(ichan,ilead))
      endif

      ! do not print out if transmission for this channel is near zero
      if (t_chan(ichan,ilead) < t_tol) return

      if (pub_on_root) then
         do iorb=1,norb_sys
            write(tag,'(a,4(i1.1,a),i5.5,a)') 'echan_c(',is,',',ien,',',ilead,',',ichan,',',iorb,')'
            ! abs
            coeff = abs(eigchan_coeff(iorb,ichan,ilead))
            call utils_qc_print(trim(tag),coeff)
         enddo
      endif

    end subroutine internal_print_qc

    !============================================!

    character(len=128) function internal_filename(ham_type,ienergy,ilead,ichan,ispin)

      !================================================================!
      ! Determines the filename of the eigenchannel file               !
      ! filename = <pub_rootname>//internal_filename()                 !
      ! The filename must have real/imag appended                      !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, April 2014                             !
      !================================================================!


      use rundat, only: pub_num_spins

      implicit none

      character(len=*) :: ham_type
      integer :: ienergy, ilead, ichan, ispin

      character(len=3)  :: spin_name
      character(len=7)  :: ham_name

      if (trim(ham_type) == 'joint') then
         ham_name = '_joint'
      else
         ham_name = ''
      endif

      if (pub_num_spins == 1) then
         spin_name = ""
      else if (ispin == 1) then
         spin_name = "_UP"
      else
         spin_name = "_DN"
      endif


      write(internal_filename,'(3(a,i2.2),a)') '_LCR'//trim(ham_name)// &
                     '_E',ienergy,'_lead',ilead,'_chan',ichan,trim(spin_name)

    end function internal_filename

  end subroutine etrans_compute_eigchan

  !==========================================================================!
  !==========================================================================!

  subroutine etrans_compute_greenf_dem(gf_dem,h_dem,s_dem,lead_hsm,energy,spin,info)

    !================================================================!
    ! Computes the LCR Green's function using DEM distributed matrix !
    ! algebra.                                                       !
    !----------------------------------------------------------------!
    ! Written by Robert Bell, April 2014                             !
    !================================================================!

    use comms, only: comms_bcast, pub_root_proc_id
    use constants, only: DP, cmplx_1
    use dense, only: DEM, dense_copy, dense_scale, dense_axpy, &
         dense_get_element, dense_put_element, dense_invert

    implicit none

    ! arguments
    type(DEM), intent(in)        :: h_dem,s_dem
    type(DEM), intent(inout)     :: gf_dem
    complex(kind=DP), intent(in) :: energy
    type(hsm_type), intent(inout):: lead_hsm(pub_nleads)
    integer, intent(in)          :: spin
    integer, intent(out)         :: info

    ! internal
    complex(kind=DP) :: element
    integer :: iorb,jorb, shift, norbl, il

    ! initialise greenf as S*E - H
    call dense_copy(gf_dem,h_dem)
    call dense_scale(gf_dem,-cmplx_1)
    call dense_axpy(gf_dem,s_dem,energy)

    ! compute self energies
    call etrans_compute_self(lead_hsm,energy,spin,info)
    if (info/=0) return

    ! synchronise self energies with root proc
    do il=1,pub_nleads
      call comms_bcast(pub_root_proc_id,lead_hsm(il)%self)
    enddo

    ! deposit self energies
    do il = 1, pub_nleads
      shift = pub_leads(il)%orbs(1) - pub_lcr%orbs(1)
      norbl = lead_hsm(il)%norb
      do iorb=1,norbl
        do jorb=1,norbl
          call dense_get_element(element, gf_dem, &
               shift+iorb, shift+jorb)
          element = element - lead_hsm(il)%self(iorb,jorb)
          call dense_put_element(element, gf_dem, &
               shift+iorb,shift+jorb)
        enddo
      enddo
    enddo

    ! invert to get Green's function
    call dense_invert(gf_dem)

    return

  end subroutine etrans_compute_greenf_dem

  !==========================================================================!
  !==========================================================================!

  subroutine etrans_eigchan_coeffs_dem(G_ret,s_dem,lead_hsm,total_trans,t_chan, &
                                       eigchan_coeff,total_dos,partial_dos)
    !================================================================!
    ! Subroutine to compute the transmission eigenchannel            !
    ! coefficients in the NGWF basis. Calculated using the method of !
    ! Paulsson and Brandbyge Phys. Rev. B 76, 115117 (2007).         !
    !                                                                !
    !----------------------------------------------------------------!
    ! Written by Robert Bell, April 2014                             !
    !================================================================!

    use constants, only: stdout, DP, cmplx_i, cmplx_0, cmplx_1, PI
    use dense, only: DEM, dense_create, dense_destroy, dense_product, &
        dense_invert, dense_normal_eigensolve, dense_put_element, &
        dense_get_col, dense_copy, dense_scale, dense_get_element, &
        dense_mat_sqrt
    use rundat, only: pub_etrans_num_eigchan
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! arguments
    type(DEM), intent(in)         :: G_ret
    type(DEM), intent(in)         :: s_dem
    type(hsm_type), intent(in)    :: lead_hsm(pub_nleads)
    real(kind=DP), intent(out)    :: total_trans
    complex(kind=DP), intent(out) :: eigchan_coeff(G_ret%nrows,pub_etrans_num_eigchan,pub_nleads)
    real(kind=DP), intent(out)    :: t_chan(pub_etrans_num_eigchan,pub_nleads)
    real(kind=DP), intent(out)    :: total_dos
    real(kind=DP), intent(out)    :: partial_dos(pub_etrans_num_eigchan,pub_nleads)


    ! internal
    type(DEM), target  :: buffer1, buffer2, buffer3
    type(DEM), pointer :: gamma_L, gamma_R
    type(DEM)          :: s_inv_half, s_half
    real(kind=DP), allocatable :: eigvals(:)
    integer :: norbl,norbr,shiftl,shiftr, il, jl
    integer :: ichan, n_eigchans, ierr, norb
    integer :: iorb,jorb
    complex(kind=DP)   :: element

    if(pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering etrans_eigchan_coeffs_dem'

    call timer_clock('etrans_eigchan_coeffs_dem',1)

    norb = G_ret%nrows
    n_eigchans = min(norb,pub_etrans_num_eigchan)

    ! allocate workspace
    call dense_create(buffer3,norb,norb,iscmplx=.true.)
    call dense_create(buffer2,norb,norb,iscmplx=.true.)
    call dense_create(buffer1,norb,norb,iscmplx=.true.)
    call dense_create(s_half,norb,norb,iscmplx=.true.)
    call dense_create(s_inv_half,norb,norb,iscmplx=.true.)
    allocate(eigvals(norb),stat=ierr)
    call utils_alloc_check('Entering etrans_eigchan_coeffs_dem','eigvals',ierr)

    ! compute DOS = -1/pi * trace(gf_dem * S)
    ! cmplx(s_dem) stored in s_half for future use
    call dense_copy(s_half,s_dem)
    call dense_product(buffer3,G_ret,s_half) ! s_half = s_dem
    total_dos=0.0_DP
    do iorb=1,norb
      call dense_get_element(element,buffer3,iorb,iorb)
      total_dos = total_dos - aimag(element)
    enddo
    total_dos = total_dos / pi

    ! compute S^1/2
    call dense_mat_sqrt(s_half,buffer2,buffer1)

    ! compute S^{-1/2}
    call dense_copy(s_inv_half,s_half)
    call dense_invert(s_inv_half)

    ! loop over source leads
    do il=1,pub_nleads

      !--- Do everything to Gamma_L first to save on memory
      ! compute gamma_L, stored in buffer3
      gamma_L => buffer3

      ! set to zero
      call dense_scale(gamma_L,cmplx_0)

      ! deposit Gamma_L
      shiftl = pub_leads(il)%orbs(1) - pub_lcr%orbs(1)
      norbl = lead_hsm(il)%norb
      do iorb=1,norbl
         do jorb=1,norbl
            element = cmplx_i*( lead_hsm(il)%self(iorb,jorb) - &
                  conjg(lead_hsm(il)%self(jorb,iorb)) )
            call dense_put_element(element,gamma_L,shiftl+iorb,shiftl+jorb)
         enddo
      enddo


      ! compute buffer1 = G_adv * GammaL * G_ret
      call dense_product(buffer2,gamma_L,G_ret)
      call dense_product(buffer1,G_ret,buffer2,opA='C')

      ! finished with gamma_L: buffer3 free again
      nullify(gamma_L)

      ! orthogonalise buffer1
      call internal_orthogonalise(buffer1,s_half,buffer2)

      ! compute   ( G^adv Gamma_L G^ret )^1/2
      call dense_mat_sqrt(buffer1,buffer3,buffer2)

      !--- Now buffer3 with gamma_R
      ! compute gamma_R
      gamma_R => buffer3
      ! set to zero
      call dense_scale(gamma_R,cmplx_0)

      ! deposit all other Gamma matrices
      do jl=1,pub_nleads
        ! do not include source lead
        if (il == jl) cycle

        shiftr = pub_leads(jl)%orbs(1) - pub_lcr%orbs(1)
        norbr = lead_hsm(jl)%norb

        do iorb=1,norbr
          do jorb=1,norbr
            element = cmplx_i*( lead_hsm(jl)%self(iorb,jorb) - &
                  conjg(lead_hsm(jl)%self(jorb,iorb)) )
            call dense_put_element(element,gamma_R,shiftr+iorb,shiftr+jorb)
          enddo
        enddo
      enddo

      ! orthogonalise
      call internal_orthogonalise(gamma_R,s_inv_half,buffer2)

      ! compute buffer1^dagger . Gamma_R . buffer1
      call dense_product(buffer2,gamma_R,buffer1)
      ! (finished with gamma_R: buffer3 free again)
      nullify(gamma_R)
      call dense_product(buffer3,buffer1,buffer2,opA='C')

      ! invert sign of buffer3, so that eigenvalues come out in descending order
      call dense_scale(buffer3,-cmplx_1)

      ! diagonalize, eigenvectors are in buffer2
      ! TODO check if calculating only the highest eigenvalues is faster...
      call dense_normal_eigensolve(norb,eigvals, buffer3, buffer2)

      ! undo sign invert of eigenvalues
      eigvals(:) = -eigvals(:)

      ! convert eigenchannels to non-orthogonal basis
      call dense_product(buffer3,s_inv_half,buffer1)
      call dense_product(buffer1,buffer3,buffer2)

      ! store transmissions
      do ichan = 1,n_eigchans
         call dense_get_col(eigchan_coeff(:,ichan,il),buffer1,ichan)
         t_chan(ichan,il) = eigvals(ichan)
      enddo

      total_trans = sum(eigvals)

      ! compute partial DOS = eigchan_coeff^H * S * eigchan_coeff
      ! buffer3 = S * eigchan_coeff
      call dense_copy(buffer2,s_dem)
      call dense_product(buffer3,buffer2,buffer1)
      ! buffer2 = eigchan_coeff^H * S
      call dense_product(buffer2,buffer1,buffer3,opA='C')

      ! diagonal entries contain the partial DOS
      do ichan = 1,n_eigchans
         call dense_get_element(element,buffer2,ichan,ichan)
         partial_DOS(ichan,il) = real(element,kind=DP)
      enddo

    enddo

    ! destroy workspace
    call dense_destroy(buffer3)
    call dense_destroy(buffer2)
    call dense_destroy(buffer1)
    call dense_destroy(s_half)
    call dense_destroy(s_inv_half)
    deallocate(eigvals,stat=ierr)
    call utils_dealloc_check('Entering etrans_eigchan_coeffs_dem','eigvals',ierr)

    call timer_clock('etrans_eigchan_coeffs_dem',2)

    if(pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving etrans_eigchan_coeffs_dem'

    return

    contains

      !=====================================!

      subroutine internal_orthogonalise(mat,s_half,zbuf)

         ! zbuf = workspace. No need to allocate extra workspace

         use dense, only: DEM, dense_product

         implicit none

         ! arguments
         type(DEM), intent(inout) :: mat
         type(DEM), intent(in)    :: s_half
         type(DEM), intent(inout) :: zbuf

         ! zbuf = mat * s_half
         call dense_product(zbuf,mat,s_half)
         ! mat = s_half * zbuf
         call dense_product(mat,s_half,zbuf)

      end subroutine internal_orthogonalise

  end subroutine etrans_eigchan_coeffs_dem

  !====================================================================!
  !====================================================================!

  subroutine etrans_print_qc(dos_lcr, trc_lcr, dos_lead, trc_lead, eref)

    use constants, only: DP
    use rundat, only:  pub_etrans_bulk, pub_etrans_lcr, pub_num_spins
    use utils, only: utils_qc_print


    implicit none

    real(kind=DP), intent(in) :: dos_lcr(:,:,:), trc_lcr(:,:,:)
    real(kind=DP), intent(in) :: dos_lead(:,:,:), trc_lead(:,:,:)
    real(kind=DP), intent(in) :: eref

    integer :: ilead, jlead, is, itrc, ienergy
    integer :: enum
    character(len=9) :: lead_string
    character(len=3) :: spin_string
    integer, parameter :: estride = 1

    call utils_qc_print('eref',eref)

    enum = size(trc_lead,dim=3)

    if (pub_etrans_bulk) then
       do ilead=1,pub_nleads
          do is=1,pub_num_spins
             write(lead_string,'(a1,2(i1.1,a1),i3.3,a1)') '(', ilead,',',is,',',0,')'
             call utils_qc_print('lead'//lead_string//'_efermi',pub_leads(ilead)%efermi(is))
             do ienergy=1,enum,estride
                write(lead_string,'(a1,2(i1.1,a1),i3.3,a1)') '(', ilead,',',is,',',ienergy,')'
                call utils_qc_print('lead'//lead_string//'_dos ',dos_lead(ilead,is,ienergy))
                call utils_qc_print('lead'//lead_string//'_trc ',trc_lead(ilead,is,ienergy))
             enddo
          enddo
       enddo
    endif

    enum = size(trc_lcr,dim=3)

    if (pub_etrans_lcr) then
       do is=1,pub_num_spins
          write(spin_string,'(a1,i1,a1)') '(',is,')'
          do ienergy=1,enum,estride
             itrc = 0
             do ilead = 1,pub_nleads-1
                do jlead = ilead+1,pub_nleads
                   itrc = itrc + 1
                   write(lead_string,'(a1,2(i1,a1),i3.3,a1)') '(', itrc,',',is,',',ienergy,')'
                   call utils_qc_print('lcr'//lead_string//'_dos ',dos_lcr(1,is,ienergy))
                   call utils_qc_print('lcr'//lead_string//'_trc ',trc_lcr(itrc,is,ienergy))
                enddo
             enddo
          enddo
       enddo
    endif

  end subroutine etrans_print_qc

end module etrans
