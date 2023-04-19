! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                           F O R C E    T E S T                              !
!=============================================================================!
!                                                                             !
! This module tests the agreement between analytically-obtained forces and    !
! those obtained by finite difference of the energy with displacements of the !
! individual ions. Activated with TASK : FORCETEST.                           !
! Module created 31/07/2013 by Nicholas Hine, from code originally written    !
! by Arash Mostofi as part of onetep.F90.                                     !
!-----------------------------------------------------------------------------!

module forcetest

  implicit none

  private

  public :: forcetest_test_forces

contains

  subroutine forcetest_test_forces(total_energy,total_forces,mdl, &
       output_file)

    !=========================================================================!
    ! Test the calculated forces against numerical finite differences         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! total_energy -                                                          !
    ! total_forces -                                                          !
    ! elements     -                                                          !
    ! nat          - number of atoms                                          !
    ! ierr         - error flag                                               !
    ! @updateme                                                               !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! numerical_forces                                                        !
    ! analytical_forces                                                       !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   elements has been properly read and initialised                       !
    !   nat properly initialised                                              !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi, v1.0, 08/03/2005                            !
    !=========================================================================!

    use comms, only: comms_bcast, pub_on_root, pub_root_proc_id
    use constants, only: DP, stdout, periodic_table_mass
    use energy_and_force, only: energy_and_force_calculate
    use model_type, only: MODEL
    use rundat, only: pub_maxit_pen, pub_devel_code, pub_maxit_ngwf_cg, pub_edft, &
         pub_read_denskern, pub_read_tightbox_ngwfs, pub_write_denskern, &
         pub_write_tightbox_ngwfs, pub_write_hamiltonian, pub_read_hamiltonian, &
         pub_write_forces, pub_zero_total_force, pub_is_auto_solvation, &
         pub_mw_total_force
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_banner, utils_assert

    implicit none

    ! Arguments
    type(MODEL), intent(inout)   :: mdl
    real(kind=DP), intent(out)   :: total_energy
    real(kind=DP), intent(out)   :: total_forces(1:3,1:mdl%nat)
    character(len=80), intent(in) :: output_file

    ! Local Variables
    integer :: atom,i
    integer :: iat, iregion, sub_at ! jcap
    integer :: ierr
    real(kind=DP) :: average_force(1:3),fd_delta_vec(1:3),fd_delta,Eplus,Eminus,Ecentral
    real(kind=DP), allocatable, dimension(:,:) :: numerical_forces
    real(kind=DP), allocatable, dimension(:,:) :: analytical_forces
    character(len=200) :: forcetest_devel_code ! ars
    integer :: fd_type, start_pos, stop_pos, test_pos ! ars
    logical :: fd_read_dkn, fd_read_tb_ngwfs, fd_read_ham ! ars
    logical :: convgd
    integer :: atom_Z
    real(kind=DP)    :: ion_mass(mdl%nat-mdl%nat_classical)
    real(kind=DP)    :: mtot
    ! Module parameters
    real(kind=dp), parameter :: electron_mass_si = 9.1093826e-31_dp
    real(kind=dp), parameter :: avogadro_si = 6.0221415e23_dp


    if (pub_on_root) then
       write(stdout,1)
       write(stdout,2) 'Starting ONETEP Force Test'
       write(stdout,1)
    endif


    ! ars: set flags
    fd_type = 1
    fd_delta = 1e-4_dp
    fd_read_dkn = .true.
    fd_read_ham = pub_edft
    fd_read_tb_ngwfs = .true.
    if (pub_on_root) then
       forcetest_devel_code=pub_devel_code
       if (len_trim(forcetest_devel_code)>0) then
          start_pos=index(forcetest_devel_code,'FORCETEST:')
          stop_pos=index(forcetest_devel_code,':FORCETEST')
          if (stop_pos<=0) stop_pos=len_trim(forcetest_devel_code) !missing end so go to end of string
          if (start_pos>0) then

             ! ars: set finite differences scheme
             test_pos=index(forcetest_devel_code,'TYPE=')
             if (test_pos>start_pos.and.test_pos<stop_pos) then
                test_pos=test_pos+len('TYPE=')
                read(forcetest_devel_code(test_pos:test_pos+ &
                     & index(forcetest_devel_code(test_pos:stop_pos),':')-2),*) fd_type
             end if
             ! ars: set FD step
             test_pos=index(forcetest_devel_code,'DELTA=')
             if (test_pos>start_pos.and.test_pos<stop_pos) then
                test_pos=test_pos+len('DELTA=')
                read(forcetest_devel_code(test_pos:test_pos+ &
                     & index(forcetest_devel_code(test_pos:stop_pos),':')-2),*) fd_delta
             end if
             ! ars: set read_desnkern flag
             if (index(forcetest_devel_code(start_pos:stop_pos),'READ_DKN=T')>0) then
                fd_read_dkn=.true.
             else if (index(forcetest_devel_code(start_pos:stop_pos),'READ_DKN=F')>0) then
                fd_read_dkn=.false.
             end if
             ! ars: set read_desnkern flag
             if (index(forcetest_devel_code(start_pos:stop_pos),'READ_HAM=T')>0) then
                fd_read_ham=.true.
             else if (index(forcetest_devel_code(start_pos:stop_pos),'READ_HAM=F')>0) then
                fd_read_ham=.false.
             end if
             ! ars: set read_tightbox_ngwfs flag
             if (index(forcetest_devel_code(start_pos:stop_pos),'READ_TB_NGWFS=T')>0) then
                fd_read_tb_ngwfs=.true.
             else if (index(forcetest_devel_code(start_pos:stop_pos),'READ_TB_NGWFS=F')>0) then
                fd_read_tb_ngwfs=.false.
             end if

          end if
       end if
       if(fd_read_dkn .and. pub_is_auto_solvation) then
          call utils_abort('Combining forcetest and auto solvation only makes &
               &sense when the displaced configurations in the forcetest do &
               &not use restarts (i.e. READ_TB_NGWFS=F and READ_DKN=F must be &
               &specified in the DEVEL_CODE).')
       end if
    end if
    call comms_bcast(pub_root_proc_id,fd_type)
    call comms_bcast(pub_root_proc_id,fd_delta)
    call comms_bcast(pub_root_proc_id,fd_read_dkn)
    call comms_bcast(pub_root_proc_id,fd_read_tb_ngwfs)

    ! Allocate forces arrays
    allocate(numerical_forces(1:3,1:mdl%nat),stat=ierr)
    call utils_alloc_check('forcetest_test_forces',&
         'numerical_forces', ierr)
    allocate(analytical_forces(1:3,1:mdl%nat),stat=ierr)
    call utils_alloc_check('forcetest_test_forces',&
         'analytical_forces', ierr)


    ! Set flags
    if (fd_read_dkn) then
       pub_write_denskern = .TRUE.
       if (pub_on_root) write(stdout,'(a)') 'Force Test: Setting write_denskern = TRUE'
    endif
    if (fd_read_tb_ngwfs) then
       pub_write_tightbox_ngwfs = .TRUE.
       if (pub_on_root) write(stdout,'(a)') 'Force Test: Setting write_tightbox_ngwfs = TRUE'
    endif
    if (fd_read_ham) then
       pub_write_hamiltonian = .TRUE.
       if (pub_on_root) write(stdout,'(a)') 'Force Test: Setting write_hamiltonian = TRUE'
    endif

    if (.not.pub_write_forces) then
       pub_write_forces = .TRUE.
       if (pub_on_root) write(stdout,'(a)') 'Force Test: Setting write_forces = TRUE'
    endif

    ! Calculate total energy and forces for initial configuration
    call energy_and_force_calculate(Ecentral,total_forces,mdl, &
         return_converged=convgd)
    total_energy = Ecentral
    if (.not.convgd.and.(pub_maxit_ngwf_cg.gt.0)) then
       call utils_abort('Error in internal_test_forces: Initial NGWF &
            &optimisation did not converge.')
    end if

    ! Copy
    do atom=1,mdl%nat
       analytical_forces(:,atom) = total_forces(:,atom)
    enddo

    ! ars: set read dkn and ngwfs after the central calculation
    pub_write_denskern = .FALSE.
    if (pub_on_root) write(stdout,'(a)') 'Force Test: Setting write_denskern = FALSE'
    pub_write_tightbox_ngwfs = .FALSE.
    if (pub_on_root) write(stdout,'(a)') 'Force Test: Setting write_tightbox_ngwfs = FALSE'
    if(fd_read_dkn) then
       pub_maxit_pen = 0
       if (pub_on_root) write(stdout,'(a)') 'Force Test: Setting maxit_pen = 0'
       pub_read_denskern = .TRUE.
       if (pub_on_root) write(stdout,'(a)') 'Force Test: Setting read_denskern = TRUE'
    endif
    if(fd_read_ham) then
       pub_maxit_pen = 0
       if (pub_on_root) write(stdout,'(a)') 'Force Test: Setting maxit_pen = 0'
       pub_read_hamiltonian = .TRUE.
       if (pub_on_root) write(stdout,'(a)') 'Force Test: Setting read_hamiltonian = TRUE'
    endif
    if (fd_read_tb_ngwfs) then
       pub_read_tightbox_ngwfs = .TRUE.
       if (pub_on_root) write(stdout,'(a)') 'Force Test: Setting read_tightbox_ngwfs = TRUE'
    endif

    ! Calculate numerical forces
    do atom=1,mdl%nat             ! loop over ions

       ! jcap: get the region counter for copying over to region
       ! element arrays
       iregion = mdl%elements(atom)%region
       ! Loop over atoms in this region until we find the
       ! appropriate atom
       sub_at=0
       do iat = 1,mdl%regions(iregion)%par%nat
          if(mdl%regions(iregion)%elements(iat)%global_atom_number &
               &== mdl%elements(atom)%global_atom_number) then
             sub_at=iat
             exit
          end if
       end do
       call utils_assert(sub_at>0, &
            'Displaced atom cannot be found in subregion')

       do i=1,3                   ! loop over Cartesian direction

          fd_delta_vec(:)=0.0_dp
          fd_delta_vec(i)=fd_delta

          ! ars: print banner
          if(pub_on_root) then
             write(stdout,'(a)') utils_banner('~', 'FORCETEST')
             write(stdout,'(a,i6,a,i2)') "Moving atom: ", atom, &
                  " along coordinate: ", i
             write(stdout,'(a,3f10.6)') "fd_delta_vec = ", fd_delta_vec(:)
             if(fd_type.eq.1) then
                write(stdout,'(a)') "Using central finite differences"
             elseif(fd_type.eq.2) then
                write(stdout,'(a)') "Using forward finite differences"
             elseif(fd_type.eq.3) then
                write(stdout,'(a)') "Using backward finite differences"
             else ! ars: ==> abort
                call utils_abort("Unknown finite differences method for task &
                     &FORCETEST")
             end if
             write(stdout,'(a)') repeat('~',80)
          end if


          if(fd_type.eq.1) then ! ars: ==> central differences

             ! Move atom in positive direction
             mdl%elements(atom)%centre%x = mdl%elements(atom)%centre%x + fd_delta_vec(1)
             mdl%elements(atom)%centre%y = mdl%elements(atom)%centre%y + fd_delta_vec(2)
             mdl%elements(atom)%centre%z = mdl%elements(atom)%centre%z + fd_delta_vec(3)
             ! jcap: copy this into the appropriate subregion
             mdl%regions(iregion)%elements(sub_at)%centre = &
                  &mdl%elements(atom)%centre

             ! Calculate new energy
             call energy_and_force_calculate(Eplus,total_forces,mdl)

             ! Move atom in negative direction
             mdl%elements(atom)%centre%x = mdl%elements(atom)%centre%x - 2.0_dp*fd_delta_vec(1)
             mdl%elements(atom)%centre%y = mdl%elements(atom)%centre%y - 2.0_dp*fd_delta_vec(2)
             mdl%elements(atom)%centre%z = mdl%elements(atom)%centre%z - 2.0_dp*fd_delta_vec(3)

             ! jcap: copy this into the appropriate subregion
             mdl%regions(iregion)%elements(sub_at)%centre = &
                  &mdl%elements(atom)%centre

             ! Calculate new energy
             call energy_and_force_calculate(Eminus,total_forces,mdl)

             ! Calculate force by finite-difference
             numerical_forces(i,atom) = (Eminus-Eplus)/(2.0_dp*fd_delta)

             ! Move atom back to original position
             mdl%elements(atom)%centre%x = mdl%elements(atom)%centre%x + fd_delta_vec(1)
             mdl%elements(atom)%centre%y = mdl%elements(atom)%centre%y + fd_delta_vec(2)
             mdl%elements(atom)%centre%z = mdl%elements(atom)%centre%z + fd_delta_vec(3)

             ! jcap: copy this into the appropriate subregion
             mdl%regions(iregion)%elements(sub_at)%centre = &
                  &mdl%elements(atom)%centre

          elseif(fd_type.eq.2) then ! ars: ==> forward differences

             ! Move atom in positive direction
             mdl%elements(atom)%centre%x = mdl%elements(atom)%centre%x + fd_delta_vec(1)
             mdl%elements(atom)%centre%y = mdl%elements(atom)%centre%y + fd_delta_vec(2)
             mdl%elements(atom)%centre%z = mdl%elements(atom)%centre%z + fd_delta_vec(3)

             ! jcap: copy this into the appropriate subregion
             mdl%regions(iregion)%elements(sub_at)%centre = &
                  &mdl%elements(atom)%centre

             ! Calculate new energy
             call energy_and_force_calculate(Eplus,total_forces,mdl)
             if (pub_on_root) write(stdout,'(a)') "Eplus = ", Eplus
             if (pub_on_root) write(stdout,'(a)') "fd_delta = ", fd_delta

             ! Calculate force by finite-difference
             numerical_forces(i,atom) = (Ecentral-Eplus)/fd_delta

             ! Move atom back to original position
             mdl%elements(atom)%centre%x = mdl%elements(atom)%centre%x - fd_delta_vec(1)
             mdl%elements(atom)%centre%y = mdl%elements(atom)%centre%y - fd_delta_vec(2)
             mdl%elements(atom)%centre%z = mdl%elements(atom)%centre%z - fd_delta_vec(3)

             ! jcap: copy this into the appropriate subregion
             mdl%regions(iregion)%elements(sub_at)%centre = &
                  &mdl%elements(atom)%centre

          elseif(fd_type.eq.3) then ! ars: ==> backward differences

             ! Move atom in negative direction
             mdl%elements(atom)%centre%x = mdl%elements(atom)%centre%x - fd_delta_vec(1)
             mdl%elements(atom)%centre%y = mdl%elements(atom)%centre%y - fd_delta_vec(2)
             mdl%elements(atom)%centre%z = mdl%elements(atom)%centre%z - fd_delta_vec(3)

             ! jcap: copy this into the appropriate subregion
             mdl%regions(iregion)%elements(sub_at)%centre = &
                  &mdl%elements(atom)%centre

             ! Calculate new energy
             call energy_and_force_calculate(Eminus,total_forces,mdl)

             ! Calculate force by finite-difference
             numerical_forces(i,atom) = (Eminus-Ecentral)/fd_delta

             ! Move atom back to original position
             mdl%elements(atom)%centre%x = mdl%elements(atom)%centre%x + fd_delta_vec(1)
             mdl%elements(atom)%centre%y = mdl%elements(atom)%centre%y + fd_delta_vec(2)
             mdl%elements(atom)%centre%z = mdl%elements(atom)%centre%z + fd_delta_vec(3)

             ! jcap: copy this into the appropriate subregion
             mdl%regions(iregion)%elements(sub_at)%centre = &
                  &mdl%elements(atom)%centre

          end if
       enddo
    enddo

    ! Calculate average force, unless not required
    if(pub_zero_total_force) then
       average_force=0.0_dp
       do atom=1,mdl%nat
          average_force(:) = average_force(:) + numerical_forces(:,atom)/mdl%nat
       enddo

       ! Subtract average force
       do atom=1,mdl%nat
          numerical_forces(:,atom) = numerical_forces(:,atom) - average_force(:)
       enddo
    end if

    ! Calculate weighted-correction
    mtot = 0.0_DP
    if (pub_mw_total_force .and. (.not. pub_zero_total_force)) then
       do atom=1,mdl%nat - mdl%nat_classical
          average_force(:) = average_force(:) + numerical_forces(:,atom)
          atom_Z = mdl%elements(atom)%atomic_number
          ion_mass(atom) = periodic_table_mass(atom_Z)*1e-3_dp/avogadro_si/electron_mass_si
          mtot = mtot + ion_mass(atom)
       end do

       ! apply correction
       do atom=1,mdl%nat - mdl%nat_classical
          numerical_forces(:,atom) = numerical_forces(:,atom) - (ion_mass(atom)/mtot)*average_force(:)
       end do
    end if

    ! Write numerical forces to file
    call forcetest_output_forces(analytical_forces,numerical_forces,mdl%nat, &
         mdl%elements,output_file)

    ! Deallocate force arrays
    deallocate(analytical_forces,stat=ierr)
    call utils_dealloc_check('forcetest_test_forces',&
         'numerical_forces', ierr)
    deallocate(numerical_forces,stat=ierr)
    call utils_dealloc_check('forcetest_test_forces',&
         'analytical_forces', ierr)

1   format(80('='))
2   format(26('<'),1x,a,1x,26('>'))

    return

  end subroutine forcetest_test_forces


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  subroutine forcetest_output_forces(analytical_forces,numerical_forces, &
       nat,elements,output_file)

    !=========================================================================!
    ! Output calculated forces and numerical finite difference forces to file !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! elements                                                                !
    ! nat          - number of atoms                                          !
    ! output_file  - name of file to which forces will be written             !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi, v1.0, 08/03/2005                            !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: DP
    use ion, only: ELEMENT
    use utils, only: utils_unit, utils_abort

    implicit none

    ! Arguments
    integer, intent(in)       :: nat
    real(kind=DP), intent(in) :: analytical_forces(1:3,1:nat)
    real(kind=DP), intent(in) :: numerical_forces(1:3,1:nat)
    type(ELEMENT), intent(in) :: elements(nat)
    character(len=80), intent(in) :: output_file

    ! Local Variables
    integer :: atom,write_unit,ios

    if (pub_on_root) then

       ! Find available unit specifier
       write_unit = utils_unit()

       ! Open file
       open(unit=write_unit,file=output_file,iostat=ios,&
            form='FORMATTED',action='WRITE')
       if (ios/=0) call utils_abort(&
            'Failed to open file "'//trim(output_file)//'"')

       ! Write analytical forces to file
       write(write_unit,'(a)') ' '
       write(write_unit,'(a)') '******************* Analytical Forces ********************'
       write(write_unit,'(a)') '*                                                        *'
       write(write_unit,'(a)') '* Element  Atom         Cartesian components (Eh/a)      *'
       write(write_unit,'(a)') '* ------------------------------------------------------ *'
       write(write_unit,'(a)') '*                       x            y            z      *'
       write(write_unit,'(a)') '*                                                        *'
       do atom=1,nat
          write(write_unit,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
               elements(atom)%symbol,' ',atom,'   ', &
               analytical_forces(1,atom),analytical_forces(2,atom),analytical_forces(3,atom),' *'
       enddo
       write(write_unit,'(a)') '*                                                        *'
       write(write_unit,'(a)') '**********************************************************'

       ! Write numerical forces to file
       write(write_unit,'(a)') ' '
       write(write_unit,'(a)') '******************** Numerical Forces ********************'
       write(write_unit,'(a)') '*                                                        *'
       write(write_unit,'(a)') '* Element  Atom         Cartesian components (Eh/a)      *'
       write(write_unit,'(a)') '* ------------------------------------------------------ *'
       write(write_unit,'(a)') '*                       x            y            z      *'
       write(write_unit,'(a)') '*                                                        *'
       do atom=1,nat
          write(write_unit,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
               elements(atom)%symbol,' ',atom,'   ', &
               numerical_forces(1,atom),numerical_forces(2,atom),numerical_forces(3,atom),' *'
       enddo
       write(write_unit,'(a)') '*                                                        *'
       write(write_unit,'(a)') '**********************************************************'

       ! Close file
       close(write_unit,iostat=ios)
       if (ios/=0) call utils_abort(&
            'Failed to close file "'//trim(output_file)//'"')

    endif

    return

  end subroutine forcetest_output_forces

end module forcetest
