! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                    Population analysis module                  !
!                                                                !
! This module performs population analysis using various methods,!
! including Mulliken analysis, Lowdin analysis, DDEC, NPA.       !
!----------------------------------------------------------------!
! Written by Nicholas Hine on 07/08/2015                         !
! Based on routines written by Peter Haynes, Robert Bell, Louis  !
! Lee and others, 2006-2015.                                     !
!================================================================!

module population

  use constants, only: DP
  use rundat, only: pub_debug_on_root

  implicit none

  private

  public :: population_analysis_mulliken
  public :: population_analysis_lowdin
  public :: population_analysis_ddec

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine population_analysis_mulliken(overlap, purkern, &
       mdl, ngwf_basis, partial_charge_orbital, ret_q_atom)

    !======================================================================!
    ! This subroutine performs Mulliken population analysis using the      !
    ! NGWFs as the set of local orbitals.                                  !
    !----------------------------------------------------------------------!
    ! Written by Peter Haynes on 31/07/2006.                               !
    ! Extended by Jacek Dziedzic in 2017/05 to allow simply returning      !
    ! q_atom(:,:) without any user feedback.                               !
    ! Modified for embedding by Joseph Prentice, August 2018               !
    !======================================================================!

    use comms, only: comms_bcast, pub_on_root, pub_my_proc_id, &
         pub_root_proc_id, comms_reduce
    use constants, only: DP, stdout, UP, DN, TWO_PI
    use function_basis, only: FUNC_BASIS
    use geometry, only: point, operator(*), operator(.dot.), operator(+),&
         magnitude
    use model_type, only: MODEL
    use parallel_strategy, only: OVERLAP_LIST,  &
         parallel_strategy_list_cross_overlaps, parallel_strategy_exit_overlaps
    use rundat, only : pub_popn_bond_cutoff, pub_print_qc, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_copy
    use sparse, only: sparse_get_block
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_heapsort, &
         utils_qc_print, utils_unit

    implicit none

    ! Arguments:
    type(SPAM3_EMBED), intent(in) :: overlap       ! Overlap matrix
    type(SPAM3_EMBED), intent(in) :: purkern(:)    ! Density kernel
    type(MODEL), intent(in) :: mdl                 ! MODEL container
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)  ! NGWF Basis type
    integer, intent(in) :: partial_charge_orbital  ! -1 except for in PC analysis
    real(kind=DP), intent(out), optional, allocatable :: ret_q_atom(:,:)

    ! Local variables:
    type(SPAM3_EMBED) :: loc_ks                    ! K.S product
    type(POINT) :: bond                            ! Bond
    type(OVERLAP_LIST) :: overlaps(mdl%nsub,mdl%nsub) ! Overlap list
    integer :: ierr                                ! Error flag
    integer :: is                                  ! Spin counter
    integer :: loc_iat                             ! Local atom counter
    integer :: iat,jat                             ! Atom counters
    integer :: iat_orig                            ! Atom iat in original order
    integer :: ingwf,jngwf                         ! NGWF counters
    integer :: ingwfs                              ! Number of NGWFs on iat
    integer :: iel                                 ! Element counter
    integer :: max_ngwfs_on_atom                   ! Max NGWFs on atom
    integer :: iovlap                              ! Overlap counter
    integer :: num_bonds                           ! The number of bonds
    integer :: nbond                               ! Number of bonds for output
    integer :: nb                                  ! Running index of bonds
    integer :: num_spins                           ! Number of spin channels
    integer :: isub,jsub                           ! Subsystem counter
    integer, allocatable :: bond_index(:)          ! Stores the index for
                                                   ! ordering output arrays
    real(kind=DP) :: block_tr                      ! Block trace
    real(kind=DP) :: popn                          ! Population
    real(kind=DP) :: charge,spin                   ! Atomic charge and spin
    real(kind=DP) :: bondvec(3)                    ! Bond vector
    real(kind=DP) :: frac(3,2)                     ! Fractional coordinates
    real(kind=DP), allocatable :: q_atom(:,:,:)    ! Atomic charges
    real(kind=DP), allocatable :: q_bond(:,:,:,:,:)! Bond charges
    real(kind=DP), allocatable :: over_block(:,:)  ! Block of overlap matrix
    real(kind=DP), allocatable :: kern_block(:,:)  ! Block of kernel
    ! agrecocmplx: need complex versions for complex matrices
    complex(kind=DP), allocatable :: over_block_cmplx(:,:)
    complex(kind=DP), allocatable :: kern_block_cmplx(:,:)
    real(kind=DP), allocatable :: bond_len(:)      ! Stores the length of each bond (for output)
    real(kind=DP), allocatable :: bond_pop(:)      ! Stores the population of each bond (for output)
    real(kind=DP), allocatable :: bond_spin(:)     ! Stores the spin of each bond (for output)
    character(len=26), allocatable :: bond_name(:) ! Stores the label of each bond (for output)
    character(len=12) :: iat_string
    logical :: partial_charge_only
    logical,allocatable :: orb_in_plane(:,:)
    integer             :: in1, in2, ufii, nplanes, jp
    logical             :: isorb
    ! agrecocmplx
    logical :: loc_cmplx
    integer :: max_overlaps

    ! -------------------------------------------------------------------------

    partial_charge_only = (partial_charge_orbital>0)

    ! agrecocmplx
    loc_cmplx = overlap%iscmplx

    !--------------------
    ! Atomic populations
    !--------------------

    ! Space for blocks
    max_ngwfs_on_atom=0
    do isub=1,mdl%nsub
       max_ngwfs_on_atom = max(max_ngwfs_on_atom,maxval(ngwf_basis(isub)%num_on_atom))
    end do
    ! agrecocmplx: allocate only quantities required for this case
    if (loc_cmplx) then
       allocate(over_block_cmplx(max_ngwfs_on_atom,max_ngwfs_on_atom),stat=ierr)
       call utils_alloc_check('population_analysis_mulliken','over_block_cmplx',ierr)
    else
       allocate(over_block(max_ngwfs_on_atom,max_ngwfs_on_atom),stat=ierr)
       call utils_alloc_check('population_analysis_mulliken','over_block',ierr)
    end if

    ! Generate block diagonal matrix for K.S product
    ! agrecocmplx
    loc_ks%structure='D'
    call sparse_embed_create(loc_ks,iscmplx=loc_cmplx)


    ! Allocate space for atomic populations
    num_spins = size(purkern)
    allocate(q_atom(mdl%nat,num_spins,mdl%nsub),stat=ierr)
    call utils_alloc_check('population_analysis_mulliken','q_atom',ierr)
    q_atom = 0.0_DP

    ! CW: code for heterostructures (tested but dormant)
    nplanes=0
    ufii=utils_unit()
    inquire(file='group_of_atoms_label_sorted',exist=isorb)
    if (isorb) then
       open(unit=ufii,file='group_of_atoms_label_sorted',form='unformatted')
          read(ufii) nplanes
          read(ufii) in1,in2
          if (allocated(orb_in_plane)) deallocate(orb_in_plane)
          allocate(orb_in_plane(in1,in2))
          read(ufii) orb_in_plane
       close(ufii)
    endif

    ! Loop over spins
    do is=1, num_spins

       ! Calculate product of density kernel and overlap
       call sparse_embed_product(loc_ks, purkern(is), overlap)

       do isub=1,mdl%nsub
          iat = mdl%regions(isub)%par%first_atom_on_proc(pub_my_proc_id)
          do loc_iat=1,mdl%regions(isub)%par%num_atoms_on_proc(pub_my_proc_id)

             ! Calculate atomic populations on this proc
             ! agrecocmplx: real and complex case
             if (loc_cmplx) then
                call sparse_get_block(over_block_cmplx,loc_ks%m(isub,isub),&
                     iat,iat)
             else
                call sparse_get_block(over_block,loc_ks%m(isub,isub),iat,iat)
             end if

             block_tr = 0.0_DP
             ! agrecocmplx: block_tr is real anyway
             if (loc_cmplx) then
                do ingwf=1,ngwf_basis(isub)%num_on_atom(iat)
                   block_tr = block_tr + real(over_block_cmplx(ingwf,ingwf),kind=DP)
                end do
             else
                do ingwf=1,ngwf_basis(isub)%num_on_atom(iat)
                   block_tr = block_tr + over_block(ingwf,ingwf)
                end do
             end if
             q_atom(mdl%regions(isub)%par%orig_atom(iat),is,isub) = block_tr
             iat = iat + 1
          end do
       end do

    end do

    ! Sum up over all procs
    call comms_reduce('SUM',q_atom)

    if(present(ret_q_atom)) then
       ! jd: Short-circuit to just return a copy of q_atom and exit
       allocate(ret_q_atom(mdl%nat,num_spins),stat=ierr)
       call utils_alloc_check('population_analysis_mulliken','ret_q_atom',ierr)
       ! rc2013: EMBED_FIX -- not yet compatible with multiple subsystems
       ret_q_atom(:,:) = q_atom(:,:,1)
    else
       ! Write out results
       if (pub_on_root) then

          ! ndmh: print results for first three atoms as a QC test
          ! rc2013: EMBED_FIX!
          if (pub_print_qc.and.(.not.partial_charge_only)) then
             if (num_spins == 1) then
                do iat=1,min(3,mdl%nat)
                   popn = 2.0_DP * q_atom(iat,1,1)
                   charge = mdl%elements(iat)%ion_charge - popn
                   write(iat_string,'(i10)') iat
                   iat_string = trim(adjustl(iat_string))
                   call utils_qc_print('atom_'//trim(iat_string)//'_population',popn)
                   call utils_qc_print('atom_'//trim(iat_string)//'_charge',charge)
                end do
             else
                do iat=1,min(3,mdl%nat)
                   popn = q_atom(iat,UP,1) + q_atom(iat,DN,1)
                   charge = mdl%elements(iat)%ion_charge - popn
                   spin = 0.5_DP * (q_atom(iat,UP,1) - q_atom(iat,DN,1))
                   write(iat_string,'(i10)') iat
                   iat_string = trim(adjustl(iat_string))
                   call utils_qc_print('atom_'//trim(iat_string)//'_population',popn)
                   call utils_qc_print('atom_'//trim(iat_string)//'_charge',charge)
                   call utils_qc_print('atom_'//trim(iat_string)//'_spin',spin)
                end do
             end if
          end if

          if (.not.partial_charge_only) then
             write(stdout,'(/a)') '    Mulliken Atomic Populations'
          else
             write(stdout,'(/a,i5)') '    Mulliken Partial Charges for Orbital ', &
                  partial_charge_orbital
          end if
          write(stdout,'(a)')  '    ---------------------------'
          if (.not.partial_charge_only) then
             if (num_spins == 1) then
                write(stdout,'(a)') 'Species  Ion    Total   Charge (e)'
                write(stdout,'(a)') '=================================='
             else
                write(stdout,'(a)') 'Species  Ion    Total   Charge (e)    &
                     &Spin (hbar)'
                write(stdout,'(a)') '======================================&
                     &==========='
             end if
          else
             if (num_spins == 1) then
                write(stdout,'(a)') 'Species      Ion    Partial Charge (e)'
                write(stdout,'(a)') '=================================='
             else
                write(stdout,'(a)') 'Species      Ion    Partial Charge (e)'
                write(stdout,'(a)') '======================================&
                     &==========='
             end if

          end if

          do isub=1,mdl%nsub
             if(mdl%nsub .gt. 1) then
                write(stdout,'(a,i1)') '        Region ', isub
                write(stdout,'(a)') '=================================='
             end if
             if (num_spins == 1) then
                do iat=1,mdl%regions(isub)%par%nat
                   if (.not.partial_charge_only) then
                      popn = 2.0_DP * q_atom(iat,1,isub)
                      charge = mdl%regions(isub)%elements(iat)%ion_charge - popn
                      write(stdout,'(2x,a2,1x,i5,2x,f10.3,2x,f9.3)') &
                           mdl%regions(isub)%elements(iat)%symbol, iat, popn, charge
                   else
                      popn = q_atom(iat,1,isub)
                      write(stdout,'(2x,a6,1x,i5,2x,f14.8)') &
                           mdl%regions(isub)%elements(iat)%species_id, iat, popn
                   end if
                end do
                write(stdout,'(a)') '=================================='
             else
                do iat=1,mdl%regions(isub)%par%nat
                   popn = q_atom(iat,UP,isub) + q_atom(iat,DN,isub)
                   charge = mdl%regions(isub)%elements(iat)%ion_charge - popn
                   spin = 0.5_DP * (q_atom(iat,UP,isub) - q_atom(iat,DN,isub))
                   if (.not.partial_charge_only) then
                      write(stdout,'(2x,a2,1x,i5,2x,f10.3,2x,f9.3,4x,f6.2)') &
                           mdl%regions(isub)%elements(iat)%symbol, iat, &
                           popn, charge, spin
                   else
                      write(stdout,'(2x,a6,1x,i5,2x,f14.8)') &
                           mdl%regions(isub)%elements(iat)%species_id, iat, popn
                   end if
                end do
                write(stdout,'(a)') '========================================&
                     &========='
             end if ! spins
          end do
       end if ! on root?
    end if ! early exit or printout?

    ! Deallocate workspace
    deallocate(q_atom,stat=ierr)
    call utils_dealloc_check('population_analysis_mulliken','q_atom',ierr)
    call sparse_embed_destroy(loc_ks)
    ! agrecocmplx
    if (loc_cmplx) then
       deallocate(over_block_cmplx,stat=ierr)
       call utils_dealloc_check('population_analysis_mulliken','over_block_cmplx',ierr)
    else
       deallocate(over_block,stat=ierr)
       call utils_dealloc_check('population_analysis_mulliken','over_block',ierr)
    end if

    if (partial_charge_only .or. present(ret_q_atom)) then
       return
    end if

    !------------------
    ! Bond populations
    !------------------

    ! Space for blocks
    max_ngwfs_on_atom=0
    do isub=1,mdl%nsub
       max_ngwfs_on_atom = max(max_ngwfs_on_atom,maxval(ngwf_basis(isub)%num_on_atom))
    end do
    ! agrecocmplx
    if (loc_cmplx) then
       allocate(over_block_cmplx(max_ngwfs_on_atom,max_ngwfs_on_atom),stat=ierr)
       call utils_alloc_check('population_analysis_mulliken','over_block_cmplx',ierr)
       allocate(kern_block_cmplx(max_ngwfs_on_atom,max_ngwfs_on_atom),stat=ierr)
       call utils_alloc_check('population_analysis_mulliken','kern_block_cmplx',ierr)
    else
       allocate(over_block(max_ngwfs_on_atom,max_ngwfs_on_atom),stat=ierr)
       call utils_alloc_check('population_analysis_mulliken','over_block',ierr)
       allocate(kern_block(max_ngwfs_on_atom,max_ngwfs_on_atom),stat=ierr)
       call utils_alloc_check('population_analysis_mulliken','kern_block',ierr)
    end if

    ! Obtain list of overlaps
!    call parallel_strategy_list_overlaps(overlaps,mdl%nat,mdl%elements, &
!         mdl%cell,'Fixed','Fixed',0.5_DP*pub_popn_bond_cutoff)
!
    num_bonds = 0
    max_overlaps = 0
    do jsub=1,mdl%nsub
       do isub=1,mdl%nsub
          call parallel_strategy_list_cross_overlaps(overlaps(jsub,isub), &
               mdl%regions(isub)%elements,mdl%regions(jsub)%elements, &
               mdl%cell,'Fixed','Fixed',0.5_DP*pub_popn_bond_cutoff)
          ! Allocate space for output arrays (find estimate for number of bonds first)
          if (pub_on_root) then
             do iat = 1, mdl%regions(isub)%par%nat
                num_bonds = num_bonds + overlaps(jsub,isub)%num_overlaps(iat)
             end do
          end if
          max_overlaps = max(max_overlaps, overlaps(isub,jsub)%max_overlaps)
       end do
    end do

    ! Allocate space for output arrays (find estimate for number of bonds first)
    !if (pub_on_root) then
    !   do iat = 1, mdl%nat
    !      num_bonds = num_bonds + overlaps%num_overlaps(iat)
    !   end do
    !end if
    call comms_bcast(pub_root_proc_id,num_bonds)

    ! Allocate space for atomic populations
    allocate(q_bond(max_overlaps,mdl%nat,num_spins,mdl%nsub,mdl%nsub), &
         stat=ierr)
    call utils_alloc_check('population_analysis_mulliken','q_bond',ierr)

    allocate(bond_pop(num_bonds),stat=ierr)
    call utils_alloc_check('population_analysis_mulliken','bond_pop',ierr)
    if (num_spins == 2) then
        allocate(bond_spin(num_bonds),stat=ierr)
        call utils_alloc_check('population_analysis_mulliken','bond_spin',ierr)
    end if
    allocate(bond_len(num_bonds),stat=ierr)
    call utils_alloc_check('population_analysis_mulliken','bond_len',ierr)
    allocate(bond_name(num_bonds),stat=ierr)
    call utils_alloc_check('population_analysis_mulliken','bond_name',ierr)
    allocate(bond_index(num_bonds),stat=ierr)
    call utils_alloc_check('population_analysis_mulliken','bond_index',ierr)
    q_bond = 0.0_DP

    ! Loop over atoms on this processor
    ! jcap: Loop over regions first
    do isub=1,mdl%nsub
       do jsub=1,mdl%nsub
          iat = mdl%regions(isub)%par%first_atom_on_proc(pub_my_proc_id)
          do loc_iat=1,mdl%regions(isub)%par%num_atoms_on_proc(pub_my_proc_id)
             ingwfs = ngwf_basis(isub)%num_on_atom(iat)
             iat_orig = mdl%regions(isub)%par%orig_atom(iat)
             do iovlap=1,overlaps(jsub,isub)%num_overlaps(iat_orig)
                jat = mdl%regions(jsub)%par%distr_atom(&
                      overlaps(jsub,isub)%overlap_list(iovlap,iat_orig))
                ! agrecocmplx
                if (loc_cmplx) then
                   call sparse_get_block(over_block_cmplx,overlap%m(jsub,isub),&
                        jat,iat)
                else
                   call sparse_get_block(over_block,overlap%m(jsub,isub),jat,iat)
                end if
                do is=1,num_spins
                   ! agrecocmplx: this will be real anyway
                   block_tr = 0.0_DP
                   ! agrecocmplx
                   if (loc_cmplx) then
                      call sparse_get_block(kern_block_cmplx, &
                           purkern(is)%m(jsub,isub), jat, iat)
                      do ingwf=1,ingwfs
                         do jngwf=1,ngwf_basis(jsub)%num_on_atom(jat)
                            block_tr = block_tr + real(over_block_cmplx(jngwf,ingwf) * &
                                 kern_block_cmplx(jngwf,ingwf),kind=DP)
                         end do
                      end do
                   else
                      call sparse_get_block(kern_block, purkern(is)%m(jsub,isub), &
                           jat, iat)
                      do ingwf=1,ingwfs
                         do jngwf=1,ngwf_basis(jsub)%num_on_atom(jat)
                            block_tr = block_tr + over_block(jngwf,ingwf) * &
                                 kern_block(jngwf,ingwf)
                         end do
                      end do
                   end if
                   q_bond(iovlap,iat_orig,is,jsub,isub) = block_tr
                end do
             end do
             iat = iat + 1
          end do

          ! Sum up over all procs
          call comms_reduce('SUM',q_bond(:,:,:,jsub,isub))
       end do
    end do

    ! Factor of two in definition
    q_bond = q_bond * 2.0_DP

    nbond = 0

    ! Write out results
    if (pub_on_root) then
       do jsub=1,mdl%nsub
          do isub=1,mdl%nsub
             ! rc2013: avoid double-counting regions
             if(jsub < isub) cycle
             do iat=1,mdl%regions(isub)%par%nat
                frac(1,1) = (mdl%regions(isub)%elements(iat)%centre &
                     .dot. mdl%cell%b1) / TWO_PI
                frac(2,1) = (mdl%regions(isub)%elements(iat)%centre &
                     .dot. mdl%cell%b2) / TWO_PI
                frac(3,1) = (mdl%regions(isub)%elements(iat)%centre &
                     .dot. mdl%cell%b3) / TWO_PI
                do iovlap=1,overlaps(jsub,isub)%num_overlaps(iat)
                   jat = overlaps(jsub,isub)%overlap_list(iovlap,iat)
                   if (isub == jsub .and. jat <= iat) cycle
                   frac(1,2) = (mdl%regions(jsub)%elements(jat)%centre &
                        .dot. mdl%cell%b1) / TWO_PI
                   frac(2,2) = (mdl%regions(jsub)%elements(jat)%centre &
                        .dot. mdl%cell%b2) / TWO_PI
                   frac(3,2) = (mdl%regions(jsub)%elements(jat)%centre &
                        .dot. mdl%cell%b3) / TWO_PI
                   ! Bond length
                   bondvec = abs(frac(:,1)-frac(:,2))
                   bondvec = modulo(bondvec,1.0_DP)
                   do iel=1,3
                      if (bondvec(iel) > 0.5_DP) bondvec(iel) = bondvec(iel) - 1.0_DP
                   end do
                   bond = bondvec(1) * mdl%cell%a1 + bondvec(2) * mdl%cell%a2 + &
                        bondvec(3) * mdl%cell%a3
                   nbond = nbond + 1
                   bond_len(nbond) = magnitude(bond)

                   write (bond_name(nbond),'(1x,a2,i5,a4,a2,i5)') &
                         mdl%regions(isub)%elements(iat)%symbol, iat, ' -- ', &
                         mdl%regions(jsub)%elements(jat)%symbol, jat

                   if (num_spins == 1) then
                      bond_pop(nbond)  = 2.0_DP * q_bond(iovlap,iat,1,isub,jsub)
                   else if(num_spins == 2) then
                      bond_pop(nbond)  = q_bond(iovlap,iat,UP,isub,jsub) &
                           + q_bond(iovlap,iat,DN,isub,jsub)
                      bond_spin(nbond) = 0.5_DP * (q_bond(iovlap,iat,UP,isub,jsub) &
                           - q_bond(iovlap,iat,DN,isub,jsub))
                   end if
                end do
             end do
          end do
       end do

       call utils_heapsort(nbond,bond_len,bond_index)

       if (num_spins == 1) then
          write(stdout,'(/a)') '        Bond         Population      &
               &Length (bohr)'
          write(stdout,'(a)') '======================================&
               &============'
          do nb = 1,nbond
              write(stdout,730)  bond_name(bond_index(nb)),bond_pop(bond_index(nb)),&
                                 bond_len(nb)
          end do
          write(stdout,'(a/)') '========================================&
               &=========='

       else
          write(stdout,'(/a)') '        Bond         Population      &
               &Spin       Length (bohr)'
          write(stdout,'(a)') '======================================&
               &======================='
          do nb = 1,nbond
              write(stdout,740)  bond_name(bond_index(nb)),bond_pop(bond_index(nb)),&
                                 bond_spin(bond_index(nb)),bond_len(nb)
          end do
          write(stdout,'(a/)') '========================================&
               &====================='
       end if

730    format(a23,t23,f9.2,t35,f10.5)
740    format(a23,t23,f9.2,t34,f8.2,t45,f10.5)

    end if

    !CW: code for heterostructures (tested but dormant)
    ! rc2013: EMBED_FIX -- not set up with embedding
    if (nplanes>0.and.isorb) then
       do jp=1,nplanes
          write(stdout,'(a)') '=================================='
          write(stdout,'(a)') '=================================='
          write(stdout,'(a,i4)') 'plane : ', jp
          if (pub_num_spins == 1) then
             do iat=1,mdl%nat
                if(.not.orb_in_plane(jp,iat)) cycle
                  popn = 2.0_DP * q_atom(iat,1,1)
                  charge = mdl%elements(iat)%ion_charge - popn
                  write(stdout,'(2x,a2,1x,i5,2x,f10.6,2x,f9.6)') &
                       mdl%elements(iat)%symbol, iat, popn, charge
             end do
             write(stdout,'(a)') '=================================='
          else
             do iat=1,mdl%nat
                if(.not.orb_in_plane(jp,iat)) cycle
                popn = q_atom(iat,UP,1) + q_atom(iat,DN,1)
                charge = mdl%elements(iat)%ion_charge - popn
                spin = 0.5_DP * (q_atom(iat,UP,1) - q_atom(iat,DN,1))
                write(stdout,'(2x,a2,1x,i5,2x,f10.6,2x,f9.6,4x,f6.2)') &
                     mdl%elements(iat)%symbol, iat, popn, charge, spin
             end do
             write(stdout,'(a)') '========================================&
                  &========='
          end if
       enddo
    endif

    ! Deallocate output arrays
    deallocate(bond_pop,stat=ierr)
    call utils_dealloc_check('population_analysis_mulliken','bond_pop',ierr)
    if (num_spins == 2) then
       deallocate(bond_spin,stat=ierr)
       call utils_dealloc_check('population_analysis_mulliken','bond_spin',ierr)
    end if
    deallocate(bond_len,stat=ierr)
    call utils_dealloc_check('population_analysis_mulliken','bond_len',ierr)
    deallocate(bond_name,stat=ierr)
    call utils_dealloc_check('population_analysis_mulliken','bond_name',ierr)
    deallocate(bond_index,stat=ierr)
    call utils_dealloc_check('population_analysis_mulliken','bond_index',ierr)
    if (allocated(orb_in_plane)) then
       deallocate(orb_in_plane)
       call utils_dealloc_check('population_analysis_mulliken','orb_in_plane',ierr)
    end if

    ! Destroy workspace
    do isub=1,mdl%nsub
       do jsub=1,mdl%nsub
          call parallel_strategy_exit_overlaps(overlaps(isub,jsub))
       end do
    end do
    deallocate(q_bond,stat=ierr)
    call utils_dealloc_check('population_analysis_mulliken','q_bond',ierr)
    ! agrecocmplx
    if (loc_cmplx) then
       deallocate(kern_block_cmplx,stat=ierr)
       call utils_dealloc_check('population_analysis_mulliken','kern_block_cmplx',ierr)
       deallocate(over_block_cmplx,stat=ierr)
       call utils_dealloc_check('population_analysis_mulliken','over_block_cmplx',ierr)
    else
       deallocate(kern_block,stat=ierr)
       call utils_dealloc_check('population_analysis_mulliken','kern_block',ierr)
       deallocate(over_block,stat=ierr)
       call utils_dealloc_check('population_analysis_mulliken','over_block',ierr)
    end if

  end subroutine population_analysis_mulliken


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine population_analysis_lowdin(overlap, purkern, &
       mdl, ngwf_basis)

    !======================================================================!
    ! This subroutine performs Lowdin population analysis using the        !
    ! NGWFs as the set of local orbitals.                                  !
    !----------------------------------------------------------------------!
    ! Written by Robert Bell, 21/07/2014, using corresponding Mulliken code!
    ! by Peter Haynes.                                                     !
    ! Modified for embedding by Joseph Prentice, August 2018               !
    !======================================================================!

    use comms, only: comms_bcast, pub_on_root, comms_reduce
    use constants, only: DP, stdout, UP, DN, TWO_PI
    use dense, only: DEM, dense_create, dense_convert, dense_mat_sqrt, &
         dense_destroy, dense_get_element, dense_product
    use function_basis, only: FUNC_BASIS
    use geometry, only: point, operator(*), operator(.dot.), operator(+),&
         magnitude
    use model_type, only: MODEL
    use rundat, only : pub_print_qc, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_copy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_heapsort, &
         utils_qc_print

    implicit none

    ! Arguments:
    type(SPAM3_EMBED), intent(in) :: overlap                   ! Overlap matrix
    type(SPAM3_EMBED), intent(in) :: purkern(:)                ! Density kernel
    type(MODEL), intent(in) :: mdl
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)              ! NGWF Basis type

    ! Local variables:
    type(DEM) :: s_sqrt                            ! square root of overlap
    type(DEM) :: sks                               ! S^(1/2) * K * S^(1/2)
    type(DEM) :: work                              ! workspace for dem products

    integer :: num_ngwfs                           ! number of NGWFs
    integer :: ierr                                ! Error flag
    integer :: is                                  ! Spin counter
    integer :: ingwf                               ! NGWF couter
    integer :: cumul_num_ngwfs
    integer :: isub                                ! Subsystem counter
    integer :: start                               ! first ngwf index
    integer :: iat                                 ! atom counter
    real(kind=DP) :: el_real                       ! element of sks
    real(kind=DP) :: block_tr                      ! block trace
    real(kind=DP) :: popn                          ! Population
    real(kind=DP) :: charge,spin                   ! Atomic charge and spin
    real(kind=DP), allocatable :: q_atom(:,:)      ! Atomic charges
    character(len=12) :: iat_string
    ! agrecocmplx
    logical :: loc_cmplx
    complex(kind=DP) :: el_cmplx

    integer :: start_reg


    ! start timer
    call timer_clock('population_analysis_lowdin',1)

    ! agrecocmplx
    loc_cmplx = overlap%iscmplx

    ! allocate dense purkern, overlap and work
    ! jcap: find total number of ngwfs
    num_ngwfs=0
    do isub=1,mdl%nsub
       num_ngwfs = num_ngwfs + ngwf_basis(isub)%num
    end do
    ! agrecocmplx
    call dense_create(s_sqrt,num_ngwfs,num_ngwfs,iscmplx=loc_cmplx)
    ! agrecocmplx
    call dense_create(sks,num_ngwfs,num_ngwfs,iscmplx=loc_cmplx)
    ! agrecocmplx
    call dense_create(work,num_ngwfs,num_ngwfs,iscmplx=loc_cmplx)


    !--------------------
    ! Atomic populations
    !--------------------

    ! Allocate space for atomic populations
    allocate(q_atom(mdl%nat,pub_num_spins),stat=ierr)
    call utils_alloc_check('population_analysis_lowdin','q_atom',ierr)
    q_atom = 0.0_DP


    ! convert to overlap to dense and get S^{1/2}
    call dense_convert(s_sqrt, overlap)
    call dense_mat_sqrt(s_sqrt)

    ! Loop over spins
    do is=1,pub_num_spins

       ! convert purkern to dense
       call dense_convert(sks, purkern(is))

       ! compute S^{1/2} * K * S^{1/2}, stored in sks
       call dense_product(work,sks,s_sqrt)
       call dense_product(sks,s_sqrt,work)


       ! jcap: loop over regions
       start_reg=0
       do isub=1,mdl%nsub
          cumul_num_ngwfs=0
          ! loop over atoms
          do iat=1,mdl%regions(isub)%par%nat

             block_tr = 0.0_DP

             num_ngwfs = ngwf_basis(isub)%num_on_atom(iat)
             cumul_num_ngwfs = cumul_num_ngwfs + num_ngwfs
             start = ngwf_basis(isub)%first_on_atom(iat) + start_reg
             do ingwf = start, start + num_ngwfs - 1
                ! agrecocmplx: atomic charges real in any case
                ! complex case
                if (loc_cmplx) then
                   call dense_get_element(el_cmplx, sks, ingwf, ingwf)
                   block_tr = block_tr + real(el_cmplx,kind=DP)
                   ! real case
                else
                   call dense_get_element(el_real, sks, ingwf, ingwf)
                   block_tr = block_tr + el_real
                end if
             enddo
             start_reg = start_reg + cumul_num_ngwfs

             q_atom(mdl%regions(isub)%par%orig_atom(iat),is) = block_tr

          enddo
       end do

    end do

    ! Write out results
    if (pub_on_root) then

       ! ndmh: print results for first three atoms as a QC test
       if (pub_print_qc) then
          if (pub_num_spins == 1) then
             do iat=1,min(3,mdl%nat)
                popn = 2.0_DP * q_atom(iat,1)
                charge = mdl%elements(iat)%ion_charge - popn
                write(iat_string,'(i10)') iat
                iat_string = trim(adjustl(iat_string))
                call utils_qc_print('atom_'//trim(iat_string)//'_lowdin_popn',popn)
                call utils_qc_print('atom_'//trim(iat_string)//'_lowdin_q',charge)
             end do
          else
             do iat=1,min(3,mdl%nat)
                popn = q_atom(iat,UP) + q_atom(iat,DN)
                charge = mdl%elements(iat)%ion_charge - popn
                spin = 0.5_DP * (q_atom(iat,UP) - q_atom(iat,DN))
                write(iat_string,'(i3)') iat
                iat_string = trim(adjustl(iat_string))
                call utils_qc_print('atom_'//trim(iat_string)//'_lowdin_popn',popn)
                call utils_qc_print('atom_'//trim(iat_string)//'_lowdin_q',charge)
                call utils_qc_print('atom_'//trim(iat_string)//'_lowdin_spin',spin)
             end do
          end if
       end if

       write(stdout,'(/a)') '    Lowdin Atomic Populations'
       write(stdout,'(a)')  '    -------------------------'
       if (pub_num_spins == 1) then
          write(stdout,'(a)') 'Species  Ion    Total   Charge (e)'
          write(stdout,'(a)') '=================================='
       else
          write(stdout,'(a)') 'Species  Ion    Total   Charge (e)    &
               &Spin (hbar)'
          write(stdout,'(a)') '======================================&
               &==========='
       end if

       if (pub_num_spins == 1) then
          do iat=1,mdl%nat
             popn = 2.0_DP * q_atom(iat,1)
             charge = mdl%elements(iat)%ion_charge - popn
             write(stdout,'(2x,a2,1x,i5,2x,f10.3,2x,f9.3)') &
                  mdl%elements(iat)%symbol, iat, popn, charge
          end do
          write(stdout,'(a)') '=================================='
       else
          do iat=1,mdl%nat
             popn = q_atom(iat,UP) + q_atom(iat,DN)
             charge = mdl%elements(iat)%ion_charge - popn
             spin = 0.5_DP * (q_atom(iat,UP) - q_atom(iat,DN))
             write(stdout,'(2x,a2,1x,i5,2x,f10.3,2x,f9.3,4x,f6.2)') &
                  mdl%elements(iat)%symbol, iat, popn, charge, spin
          end do
          write(stdout,'(a)') '========================================&
               &========='
       end if
    end if

    ! Deallocate workspace
    deallocate(q_atom,stat=ierr)
    call utils_dealloc_check('population_analysis_lowdin','q_atom',ierr)

    !------------------
    ! Bond populations
    !------------------

    ! Bond populations are zero by definition



    ! Destroy workspace
    call dense_destroy(s_sqrt)
    call dense_destroy(sks)
    call dense_destroy(work)
    call timer_clock('population_analysis_lowdin',2)

  end subroutine population_analysis_lowdin


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine population_analysis_ddec(denskern,rep,ngwf_basis,mdl)

    !----------------------------------------------------------------------!
    ! lpl: 17-01-2014 Subroutine to initialize and call the DDEC module    !
    !----------------------------------------------------------------------!
    ! lpl: Subroutine to initialize stuff for the DDEC module (ddec_mod)   !
    ! lpl: MOD#07 - If core density is provided by ONETEP initialize them  !
    !      at this mark                                                    !
    ! lpl: MOD#10 - Ignore shaved valence density instead of adding them   !
    !      back to the core. Should not have any effect on pspot densities !
    !      which have no sharp density spikes anyway                       !
    ! Modified for embedding by Joseph Prentice, August 2018               !
    !----------------------------------------------------------------------!

    use comms, only: comms_barrier, comms_reduce, pub_on_root
    use constants, only: stdout
    use ddec, only: ddec_main ! lpl: c,f. COMMENT#06 in ddec_mod , ddec_rmse
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use geometry, only: operator(.DOT.), operator(.CROSS.)
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_ddec_use_coredens, &
         pub_ddec_format_dens, pub_ddec_IH_frac, pub_ddec_zero_thresh, &
         pub_spin_fac, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_axpy, sparse_embed_scale, &
         sparse_embed_create, sparse_embed_destroy
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_axpy
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(NGWF_REP), intent(in) :: rep
    type(SPAM3_EMBED), intent(inout) :: denskern(:)
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(MODEL), intent(in) :: mdl

    ! lpl: Local variables
    real(kind=DP), allocatable :: density_fine(:,:,:,:), buffer_fine(:,:,:,:)
    real(kind=DP), allocatable :: pseudodensity_fine(:,:,:,:)

    real(kind=DP), target, allocatable :: density1_fine(:,:,:,:)
    real(kind=DP), target, allocatable :: density2_fine(:,:,:,:)
    real(kind=DP), target, allocatable :: density3_fine(:,:,:,:)
    real(kind=DP), target, allocatable :: density4_fine(:,:,:,:)

    ! lpl: To replace passing of null pointers to 'ddec_main' which can
    !      crash debug code
    real(kind=DP), target :: dummy_density(1,1,1,1)

    ! lpl: I heards 'allocate' is more efficient because it
    !      is almost always contiguous......anyone?
    real(kind=DP), pointer :: ref_pseudodensity_fine(:,:,:,:)
    real(kind=DP), pointer :: Yavg_pseudodensity_fine(:,:,:,:)
    real(kind=DP), pointer :: dcore_density_fine(:,:,:,:)
    real(kind=DP), pointer :: core_pseudodensity_fine(:,:,:,:)

    real(kind=DP) :: total_integrated_dens
    type(SPAM3), allocatable :: spinless_kern(:)
    integer :: isp, it1, it2, it3, old_num_spins, ierr
    integer :: isub,jsub

    logical :: include_refdens

    real(kind=DP), allocatable :: atom_charge(:)

    ! lpl: Density correction stuff as per CHARGEMOL
    ! lpl: 'ddec_max_density_val' is called 'pixel_integration_tolerance'
    !      in CHARGEMOL. Unit in electrons
    real(kind=DP), parameter :: ddec_max_density_val = 0.03_DP ! in e
    real(kind=DP), parameter :: ddec_neg_dens_thresh = 0.01_DP ! in e
    real(kind=DP), parameter :: ddec_neg_dens_thresh_frac = 0.0003_DP
    real(kind=DP) :: sum_neg_density(2), ddec_dens_zero_thresh
    real(kind=DP) :: sum_shaved_density, max_shaved_density, shaved_density
    real(kind=DP) :: vox_vol, max_density_val
    integer :: num_shaved_density

    ! lpl: This is meant as a placeholder to indicate whether the core
    !      density is available from ONETEP
    logical :: coredens_provided

    if ( pub_debug_on_root ) write(stdout,'(a)') ' DEBUG: Entering population_analysis_ddec'

    ierr = 0

!------------------------------------------------------------------------------!
! lpl: Initialization                                                          !
!------------------------------------------------------------------------------!

    if ( pub_on_root ) then
       write(stdout,'(a)') '================================================'
       write(stdout,'(a)') '                      DDEC'
       write(stdout,'(a)') '------------------------------------------------'
    end if

    ! lpl: Voxel volume
    vox_vol = mdl%fine_grid%da1.DOT.&
         (mdl%fine_grid%da2.CROSS.mdl%fine_grid%da3)

    nullify(ref_pseudodensity_fine)
    nullify(Yavg_pseudodensity_fine)
    nullify(dcore_density_fine)
    nullify(core_pseudodensity_fine)

    ! lpl: Remove this and other associated declaration if
    !      atom_charge is to be passed to other routines
    allocate(atom_charge(mdl%nat),stat=ierr)
    call utils_alloc_check('population_analysis_ddec','atom_charge',ierr)

    ! lpl: Required in all cases - real (total/valence) density and
    !      pseudodensity
    allocate(density_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,1),stat=ierr)
    call utils_alloc_check('population_analysis_ddec','density_fine',ierr)
    allocate(buffer_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,1),stat=ierr)
    call utils_alloc_check('population_analysis_ddec','buffer_fine',ierr)
    density_fine = 0.0_DP
    allocate(pseudodensity_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12,1),stat=ierr)
    call utils_alloc_check('population_analysis_ddec','pseudodensity_fine',ierr)
    pseudodensity_fine = 0.0_DP

    ! lpl: Required if chi > 0 i.e. uses reference density - reference
    !      pseudodensity and Yavg pseudodensity
    if ( pub_ddec_IH_frac > 0.0_DP ) then
       include_refdens = .true.
       allocate(density1_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12,1),stat=ierr)
       call utils_alloc_check('population_analysis_ddec','density1_fine',ierr)
       ref_pseudodensity_fine => density1_fine
       ref_pseudodensity_fine = 0.0_DP

       allocate(density2_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12,1),stat=ierr)
       call utils_alloc_check('population_analysis_ddec','density2_fine',ierr)
       Yavg_pseudodensity_fine => density2_fine
       Yavg_pseudodensity_fine = 0.0_DP
    else
       include_refdens = .false.
       ref_pseudodensity_fine => dummy_density
       Yavg_pseudodensity_fine => dummy_density
! Replaced with above since passing null pointers to 'ddec_main' could cause
! debugging to crash
!           nullify(ref_pseudodensity_fine)
!           nullify(Yavg_pseudodensity_fine)
    end if

    ! lpl: Required if there are core charges to be considered - core
    !      density and core pseudodensity
    ! lpl : Note: Be aware that there is a 'core_density_fine' defined
    !       in the main properties_mod routine
    if ( pub_ddec_use_coredens ) then
       allocate(density3_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12,1),stat=ierr)
       call utils_alloc_check('population_analysis_ddec','density3_fine',ierr)
       dcore_density_fine => density3_fine
       dcore_density_fine = 0.0_DP

       allocate(density4_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
            mdl%fine_grid%max_slabs12,1),stat=ierr)
       call utils_alloc_check('population_analysis_ddec','density4_fine',ierr)
       core_pseudodensity_fine => density4_fine
       core_pseudodensity_fine = 0.0_DP
    else
       dcore_density_fine => dummy_density
       core_pseudodensity_fine => dummy_density
! Replaced with above since passing null pointers to 'ddec_main' could cause
! debugging to crash
!           nullify(dcore_density_fine)
!           nullify(core_pseudodensity_fine)
    end if

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! lpl: Obtain valence density from density kernel (coredens goes here too)     !
!------------------------------------------------------------------------------!

    allocate(spinless_kern(1),stat=ierr)
    call utils_alloc_check('population_analysis_ddec','spinless_kern',ierr)

    ! lpl: Get total spinless density into density_fine
    !      Also fiddles with pub_num_spins for a while
    old_num_spins = pub_num_spins
    pub_num_spins = 1

    ! jcap: Loop over regions
    do jsub=1,mdl%nsub
       do isub=1,mdl%nsub
          buffer_fine=0.d0

          ! jcap: Create spinless denskern on a region-by-region basis
          call sparse_create(spinless_kern(1),denskern(1)%m(isub,jsub))

          ! lpl: Get spinless denskern
          do isp=1,pub_num_spins
             call sparse_axpy(spinless_kern(1),denskern(isp)%m(isub,jsub),pub_spin_fac)
          end do

          call density_on_grid(buffer_fine,mdl%fine_grid,mdl%dbl_grid, &
               mdl%cell, mdl%fftbox,spinless_kern,rep%ngwf_overlap%m(isub,jsub), &
               rep%ngwfs_on_grid(isub),ngwf_basis(isub),&
               rep%ngwfs_on_grid(jsub),ngwf_basis(jsub))

          density_fine=density_fine+buffer_fine

          call sparse_destroy(spinless_kern(1))
       end do
    end do

    pub_num_spins = old_num_spins

    ! lpl: Total integrated density
    total_integrated_dens = sum(density_fine)

    ! lpl: Reset densities
    call comms_barrier
    call comms_reduce('SUM',total_integrated_dens)

    deallocate(spinless_kern,stat=ierr)
    call utils_dealloc_check('population_analysis_ddec','spinless_kern',ierr)

    ! lpl: MOD#07 - If there are core densities to be provided do it at
    !      latest here, and change 'coredens_provided' to '.true.'
!        if ( coredens_provided ) dcore_density_fine = <add_core_density>
    coredens_provided = .false.

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! lpl: Format valence (and total) densities                                    !
!------------------------------------------------------------------------------!

    if ( pub_ddec_format_dens ) then
       ! lpl: Make sure any negative density 'noise' is not severe
       ddec_dens_zero_thresh = pub_ddec_zero_thresh
       sum_neg_density = 0.0_DP

       do it1=1,mdl%fine_grid%ld1
          do it2=1,mdl%fine_grid%ld2
             do it3=1,mdl%fine_grid%max_slabs12

                if ( density_fine(it1,it2,it3,1) < 0.0_DP ) then
! This is rubbish
!                       if ( abs(density_fine(it1,it2,it3,1)) &
!                           <= ddec_dens_zero_thresh ) then
                   sum_neg_density(1) = &
                        sum_neg_density(1) + density_fine(it1,it2,it3,1)
                   density_fine(it1,it2,it3,1) = 0.0_DP
!                       else
!                          write(error_string,'(e20.10)') &
!                               density_fine(it1,it2,it3,1)
!                          call utils_abort(' ERROR in population_analysis_ddec:&
!                               & density_fine < 0.0_DP ('//&
!                               trim(error_string)//')')
!                       end if
                end if

                if ( pub_ddec_use_coredens ) then
                   if ( dcore_density_fine(it1,it2,it3,1) < 0.0_DP ) then
!                          if ( abs(dcore_density_fine(it1,it2,it3,1)) &
!                               <= ddec_dens_zero_thresh ) then
                      sum_neg_density(2) = sum_neg_density(2) + &
                           dcore_density_fine(it1,it2,it3,1)
                      dcore_density_fine(it1,it2,it3,1) = 0.0_DP
!                          else
!                             write(error_string,'(e20.10)') &
!                                  dcore_density_fine(it1,it2,it3,1)
!                             call utils_abort(' ERROR in population_analysis_ddec:&
!                                  & dcore_density_fine < 0.0_DP ('//&
!                                  trim(error_string)//')')
!                          end if
                   end if
                end if

             end do
          end do
       end do

       sum_neg_density = -1.0_DP*vox_vol*sum_neg_density
       call comms_reduce('SUM',sum_neg_density)

       ddec_dens_zero_thresh = max( ddec_neg_dens_thresh, &
            ddec_neg_dens_thresh_frac*vox_vol*total_integrated_dens )

       if ( pub_on_root ) then
          write(stdout,'(a,e16.8)') ' Sum negative valence (e) : ', &
               sum_neg_density(1)
          if ( pub_ddec_use_coredens .and. coredens_provided ) &
               write(stdout,'(a,e16.8)') ' Sum negative core (e)    : ', &
                    sum_neg_density(2)

          if ( sum_neg_density(1) > ddec_dens_zero_thresh ) &
               call utils_abort(' ERROR in population_analysis_ddec:&
                    & sum_neg_density(1) > negative thresholds')
          if ( sum_neg_density(2) > ddec_dens_zero_thresh ) &
               call utils_abort(' ERROR in population_analysis_ddec:&
                    & sum_neg_density(2) > negative thresholds')
       end if

       ! lpl: Format density as per CHARGEMOL - check for density spikes etc.
       !      according to stuff done in 'format_<something>_cube_density'
       if ( pub_on_root ) write(stdout,'(a,/a,/a)') &
            " WARNING: The function 'ddec_format_dens' has", &
            ' not been thoroughly tested and is likely', &
            ' very buggy. Use at your own risk.'

       ! Get maximum voxel contribution allowed (i.e. max point value)
       max_density_val = ddec_max_density_val/vox_vol
       sum_shaved_density = 0.0_DP
       max_shaved_density = 0.0_DP
       num_shaved_density = 0

       do it1=1,mdl%fine_grid%ld1
          do it2=1,mdl%fine_grid%ld2
             do it3=1,mdl%fine_grid%max_slabs12

                if ( density_fine(it1,it2,it3,1) > max_density_val ) then

                   ! lpl: How much excess density is being shaved off
                   shaved_density = &
                        density_fine(it1,it2,it3,1) - max_density_val

                   ! lpl: Dump the excess density into the core if possible
                   ! lpl: MOD#10 - This is appararently not done, so ignore
                   !      the shaved density
! lpl: MOD#10
!                       if ( pub_ddec_use_coredens ) then
!                          dcore_density_fine(it1,it2,it3,1) = &
!                               dcore_density_fine(it1,it2,it3,1) + &
!                               shaved_density
!                       end if

                   ! lpl: Shave max voxel value to 'max_density'
                   density_fine(it1,it2,it3,1) = max_density_val
                   sum_shaved_density = sum_shaved_density + shaved_density
                   max_shaved_density = max(max_shaved_density, &
                        shaved_density)
                   num_shaved_density = num_shaved_density + 1

                end if

             end do
          end do
       end do

       sum_shaved_density = vox_vol*sum_shaved_density
       max_shaved_density = vox_vol*max_shaved_density

       call comms_reduce('SUM', sum_shaved_density)
       call comms_reduce('MAX', max_shaved_density)
       call comms_reduce('SUM', num_shaved_density)

       if ( pub_on_root ) then
          write(stdout,'(a)') ' Density formatting complete'
          if ( num_shaved_density > 0 ) then
             write(stdout,'(a)') ' Density spike(s) have been shaved:'
             write(stdout,'(a,e20.10)') '    Total (e) : ', &
                  sum_shaved_density
             write(stdout,'(a,e20.10)') '    Max   (e) : ', &
                  max_shaved_density
             write(stdout,'(a,I10)') '    Num       : ', &
                  num_shaved_density
! MOD#10
!                 if ( pub_ddec_use_coredens ) write(stdout,'(a,/a)') &
!                      ' Shaved densities have been dumped', &
!                      ' into the core density'
          end if
       end if
    end if ! END  if ( pub_ddec_format_densities )

    if ( pub_on_root ) write(stdout,'(a)') &
         '------------------------------------------------'

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! lpl: Run DDEC and finalize                                                   !
!------------------------------------------------------------------------------!

    ! lpl: Main DDEC subroutine. Main output is 'atom_charge' that
    !      could be used elsewhere
    if(mdl%nsub.gt.1) call utils_abort(&
         "DDEC is not ready to be used with more than one subsystem yet.")
    call ddec_main(atom_charge, density_fine(:,:,:,1), &
         pseudodensity_fine(:,:,:,1), ref_pseudodensity_fine(:,:,:,1), &
         Yavg_pseudodensity_fine(:,:,:,1), dcore_density_fine(:,:,:,1), &
         core_pseudodensity_fine(:,:,:,1), &
         mdl%fine_grid, mdl, total_integrated_dens, &
         include_refdens, pub_ddec_use_coredens, coredens_provided)

    if ( allocated(density4_fine) ) then
       deallocate(density4_fine, stat=ierr)
       call utils_dealloc_check('population_analysis_ddec','density4_fine',ierr)
    end if
    nullify(core_pseudodensity_fine)

    if ( allocated(density3_fine) ) then
       deallocate(density3_fine, stat=ierr)
       call utils_dealloc_check('population_analysis_ddec','density3_fine',ierr)
    end if
    nullify(dcore_density_fine)

    if ( allocated(density2_fine) ) then
       deallocate(density2_fine, stat=ierr)
       call utils_dealloc_check('population_analysis_ddec','density2_fine',ierr)
    end if
    nullify(Yavg_pseudodensity_fine)

     if ( allocated(density1_fine) ) then
! lpl: c.f. COMMENT#06 in ddec_mod
!        if ( allocated(density1_fine) .and. .not. pub_ddec_rmse ) then
       deallocate(density1_fine, stat=ierr)
       call utils_dealloc_check('population_analysis_ddec','density1_fine',ierr)
    end if
    nullify(ref_pseudodensity_fine)

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! lpl: Calculate RMSE if required and finalize                                 !
!------------------------------------------------------------------------------!
! lpl: c.f. COMMENT#06 in ddec_mod
!
!        ! lpl: Calculate RMS V(r) between DDEC charges and DFT electrostatic
!        !      potential. The 'ddec_rmse' subroutine is currently sloppy.
!        if(pub_ddec_rmse) then
!           ! lpl: Resue these grids to store Hartree potentials
!           pseudodensity_fine = 0.0_DP
!           jd: Allocates that are both commented out *and* continued are too
!               much for the alloc_checker conversion script. I put the dashes
!               in to hide them.
!           if ( .not. allocated(density1_fine) ) then
!              allocate(density1_fine(mdl%fine_grid%ld1,mdl%fine_grid%ld2, &
!                   mdl%fine_grid%max_slabs12,1),stat=ie-rr)
!              call utils_allo-c_check('population_analysis_ddec','density1_fine',ie-rr)
!           end if
!
!           ! Scale denskern for spin degeneracy
!           if (pub_num_spins == 1) then
!              call sparse_scale(denskern(1),pub_spin_fac)
!           end if
!           ! lpl: Obtain electrostatic potential into 'pseudodensity_fine'
!           !      i.e. Hartree + ion w/o xc potential
!           call hamiltonian_lhxc_calculate(pseudodensity_fine,lhxc_energy, &
!                ham%dijhat,mdl%fine_grid,localpseudo_fine,core_density_fine, &
!                rep%ngwfs_on_grid,ngwf_basis,denskern,rep%ngwf_overlap, &
!                rep%sp_overlap,add_xc_pot=.false.)
!           call visual_scalarfield(pseudodensity_fine(:,:,:,1),mdl%fine_grid, &
!                'Hartree potential fed to DDEC RMSe (in eV) for:',&
!                '_hartree_potential',elements,-HARTREE_IN_EVS)
!           ! Unscale denskern
!           if (pub_num_spins == 1) then
!              call sparse_scale(denskern(1),1.0_DP/pub_spin_fac)
!           end if
!
!           ! lpl: Calculate the RMSe V(r)
!           ! lpl: 'density1_fine' stores the point charge V(r) and is
!           !      an optional argument
!           call ddec_rmse(pseudodensity_fine,mdl%fine_grid,elements, &
!                     atom_charge,density1_fine)
!           call visual_scalarfield( &
!                density1_fine(:,:,:,1), mdl%fine_grid, &
!                'DDEC point charge potential (in eV) for:', &
!                '_ddec_potential', elements, -HARTREE_IN_EVS)
!        end if
!
!        if ( allocated(density1_fine) ) then
!           deallocate(density1_fine,stat=ierr)
!           call utils_dealloc_check('population_analysis_ddec','density1_fine',ierr)
!        end if
!        deallocate(pseudodensity_fine,stat=ierr)
!        call utils_dealloc_check('population_analysis_ddec','pseudodensity_fine',ierr)
!        ! lpl: This 'density_fine' has not been used anywhere after the
!        !      completion of ddec_main, so it could be deleted way before
!        !      the RMSE routine if deisred. It should still contain the
!        !      correct valence density generated from the density kernel
!        deallocate(density_fine,stat=ierr)
!        call utils_dealloc_check('population_analysis_ddec','density_fine',ierr)
!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! lpl: Atomic charges are here                                                 !
!------------------------------------------------------------------------------!

    ! lpl: Remove this and other associated declaration if
    !      atom_charge is to be passed to other routines
    deallocate(atom_charge,stat=ierr)
    call utils_dealloc_check('population_analysis_ddec','atom_charge',ierr)

!------------------------------------------------------------------------------!

    if ( pub_debug_on_root ) write(stdout,'(a)') ' DEBUG: Leaving population_analysis_ddec'

  end subroutine population_analysis_ddec


end module population
