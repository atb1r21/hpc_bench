! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-


!============================================================================!
!                                                                            !
!                        Molecular Dynamics Module                           !
!                                                                            !
!----------------------------------------------------------------------------!
! Written by by Simon M.-M. Dubois (May 2010) based on the original          !
! module created by Arash A. Mostofi (May 2006)                              !
!============================================================================!

module md

  use constants, only: dp, stdout
  use rundat, only: pub_debug_on_root, md_iter_global

  implicit none

  private

  ! MD status
  real(kind=dp), save :: md_time
  real(kind=dp), save :: md_energy
  real(kind=dp), save :: md_kinetic_energy
  real(kind=dp), save :: md_potential_energy
  real(kind=dp), save :: md_actual_temp

  ! Module parameters
  real(kind=DP), parameter :: aut2fs = 0.0241888468_dp ! atomic unit of time --> fs
  real(kind=DP), parameter :: ha2k   = 3.1577462e5_dp  ! hartree --> kelvin
  character (len=30), parameter :: time_label = 'fs'
  character (len=30), parameter :: energy_label = 'Ha'
  character (len=30), parameter :: temperature_label = 'K'

  ! MTS integration parameters
  integer, save :: main_maxit_ngwf_cg
  real(kind=dp), save :: main_ngwfs_threshold
  real(kind=dp), save :: main_lnv_threshold
  real(kind=dp), save :: main_elec_energy_tol
  real(kind=dp), save :: main_elec_force_tol
  real(kind=dp), save :: main_ngwf_max_grad
  integer, save :: main_maxit_lnv
  integer, save :: main_minit_lnv
  integer, save :: main_maxit_pen

  ! List of subroutines and functions :
  public  :: md_main

  contains

!==============================================================================!
!==============================================================================!

    subroutine md_main(total_energy,forces,mdl)

      ! Modified to remove pub_par by Joseph Prentice, May 2018

      use constants,       only: dp,stdout
      use comms,           only: pub_on_root
      use md_ions,         only: md_allocate_struct, md_deallocate_struct
      use model_type,      only: MODEL
      use rundat,          only: pub_mix_bcs
      use utils,           only: utils_abort

      implicit none

      ! Arguments
      type(MODEL),      intent(inout) :: mdl
      real(kind=DP),    intent(inout) :: total_energy
      real(kind=DP),    intent(inout) :: forces(1:3,mdl%nat)

      ! Internal variables

      if (pub_mix_bcs .and. pub_on_root) then
         call utils_abort('You are attempting a calculation with &
              &inconsistent boundary conditions (neither fully periodic, nor &
              &fully open), which is not allowed in MD.')
      end if

      ! Banner
      if (pub_on_root) write(stdout,1) ' Starting ONETEP Molecular Dynamics '

      ! Allocate module arrays
      call md_allocate_struct(mdl%elements,mdl%nat)

      ! Velocity verlet integration algorithm
      call md_velocity_verlet(total_energy,forces,mdl)

      ! Deallocate module arrays
      call md_deallocate_struct()

      if (pub_on_root) write(stdout,1) '  End of ONETEP Molecular Dynamics  '

1     format(/,80('x'),/,11('x'),a47,11x,11('x'),/,80('x'),/)

      return

    end subroutine md_main

!==============================================================================!
!==============================================================================!

    subroutine md_velocity_verlet(total_energy,forces,mdl)

    !====================================================================!
    ! Velocity Verlet algorithm :                                        !
    ! Integrate the equation of motion for a system of classical         !
    ! particles, if required coupled with one of the following           !
    ! thermostat (ANDERSON/LANGEVIN/NOSE-HOOVER/NOSE-HOOVER-CHAIN)       !
    !                                                                    !
    ! Calls to thermostat subroutine are interleaved with the usual      !
    ! velocity verlet algorythm:                                         !
    !                                                                    !
    ! thermo (step1)                                                     !
    !    v(t+dt/2) = v(t) + a(t)*dt/2                                    !
    ! thermo (step2)                                                     !
    !    r(t+dt) = r(t) + v(t+dt/2)*dt                                   !
    !    call energy_and_force_calculate()                               !
    ! thermo (step3)                                                     !
    !    v(t+dt) = v(t+dt/2) + a(t+dt)*dt/2                              !
    ! thermo (step4)                                                     !
    !                                                                    !
    !--------------------------------------------------------------------!
    ! Implemented by Simon M.-M. Dubois (May 2011)                       !
    ! Modified to remove pub_par by Joseph Prentice, May 2018            !
    !====================================================================!

      use comms,            only: pub_on_root
      use constants,        only: dp,stdout
      use electronic_history, only: elec_history_initialise_method, &
           elec_history_create_storage, elec_history_destroy_storage, &
           elec_history_reset, elec_history_restore_history_info
      use md_ions,          only: a_t, v_t, r_t, ion_mass, v0_t, &
           group_num, md_initialise_groups, md_initialise_velocities, &
           md_initialise_positions, md_initialise_ion_masses, &
           kinetic_energy, md_restore_trajectory, md_atomic_to_internal_vel, &
           md_scale_velocities
      use md_thermostat,    only: md_select_thermostat, md_thermo, &
           md_restore_thermostat, md_backup_thermostat, md_write_thermostat, &
           md_velocity_verlet_thermostat, md_initialise_thermostats, &
           md_initialise_ndof
      use model_type,       only: MODEL
      use rundat,           only: md_num_iter,md_delta_t, &
           pub_write_positions, mix_dkn_type, mix_dkn_num, &
           mix_dkn_reset, mix_dkn_init_type, mix_dkn_init_num, &
           mix_ngwfs_type, mix_ngwfs_num, mix_ngwfs_reset, mix_ngwfs_init_type, &
           mix_ngwfs_init_num, mix_local_length, mix_local_smear, &
           mts_nstep, pub_write_velocities, mts_xi, md_restart, &
           md_global_restart, md_write_history, md_restart_thermo, &
           md_reset_history, pub_read_denskern, pub_read_tightbox_ngwfs, &
           mix_ngwfs_coeff, pub_read_sw_ngwfs
      use utils,            only: utils_alloc_check, utils_dealloc_check, &
           utils_assert

      implicit none

      ! Argument
      type(MODEL),      intent(inout) :: mdl
      real(kind=DP),    intent(inout) :: total_energy
      real(kind=DP),    intent(inout) :: forces(1:3,mdl%nat)

      ! Internal variables
      real(kind=DP), allocatable :: tgroup(:), new_tgroup(:)
      integer :: md_iter, md_step, mts_step
      integer :: ndof, ndof_r, new_ndof_r, md_maxstep
      integer :: iat, ig, ith
      integer :: ierr, ierr1, ierr2
      logical :: write_global

    !====================================================================!

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_velocity_verlet'


      !###############################################################
      !#### Trajectory and thermostat initialisation #################

      ! Initialise MD
      if (mts_nstep < 1) mts_nstep = 1
      md_time = 0.0_dp
      md_iter = 0
      md_iter_global = 0 ! This is overwritten if md_global_restart
      mts_step = 0
      ierr1 = 0; ierr2 = 0
      allocate(tgroup(group_num),stat=ierr)
      call utils_alloc_check('md_velocity_verlet','tgroup',ierr)
      allocate(new_tgroup(group_num),stat=ierr)
      call utils_alloc_check('md_velocity_verlet','new_tgroup',ierr)


      ! Initialise multi-time-step integration scheme
      if (pub_on_root) write(stdout,'(a)',advance='no') &
          '- MD -: Initialising multi time-step (MTS) integration scheme...'
      call init_mts_scheme()
      if (pub_on_root) write(stdout,*) 'done'

      ! Initialise ion masses
      if (pub_on_root) write(stdout,'(a)',advance='no') &
          '- MD -: Initialising ion masses...'
      call md_initialise_ion_masses(mdl%elements,mdl%nat)
      if (pub_on_root) write(stdout,*) 'done'

      ! Initialise groups
      if (pub_on_root) write(stdout,'(a)',advance='no') &
          '- MD -: Initialising ion groups...'
      call md_initialise_groups(mdl%elements,mdl%nat)
      if (pub_on_root) write(stdout,*) 'done'

      ! Read data from previous calculation
      if (md_restart) then

         ! Restore time, positions and velocities from
         ! previous calculation
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             '- MD -: Restore structure from previous calculation...'
         call md_restore_trajectory(md_time,ierr1,mdl%nat)
         if (pub_on_root) then
            if (ierr1 ==0) then
               write(stdout,*) 'done'
            else
               write(stdout,*) 'failed'
            endif
         endif

         if (md_restart_thermo) then
            if (ierr1 == 0) then
               ! Restore thermostat parameters from previous
               ! calculation
               if (pub_on_root) write(stdout,'(a)',advance='no') &
                   '- MD -: Restore thermostats from previous calculation...'
               call md_restore_thermostat(ierr2)
               if (pub_on_root) then
                  if (ierr2 ==0) then
                     write(stdout,*) 'done'
                  else
                     write(stdout,*) 'failed'
                  endif
               endif

               if (ierr2 .ne. 0) then
                  ! Initialise thermostat parameters
                  if (pub_on_root) write(stdout,'(a)',advance='no') &
                      '- MD -: Initialising thermostats...'
                  call md_initialise_thermostats()
                  if (pub_on_root) write(stdout,*) 'done'
               endif
            endif
            ! vv: Initialise tgroup
            do ig=1, group_num
               call md_select_thermostat(ith,ig,md_time)
               tgroup(ig) = md_thermo(ith)%T
            end do
         else
            ! Initialise thermostat parameters
            if (pub_on_root) write(stdout,'(a)',advance='no') &
                '- MD -: Initialising thermostats...'
            call md_initialise_thermostats()
            ! vv: Initialise tgroup
            if (pub_on_root) write(stdout,*) 'done'
            do ig=1, group_num
               call md_select_thermostat(ith,ig,md_time)
               tgroup(ig) = md_thermo(ith)%T
            end do
         end if


         ! vv: Initialise number of degrees of freedom
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             '- MD -: Initialising degrees of freedom...'
         call md_initialise_ndof(ndof_r, md_time, mdl%nat)
         ndof = 3*mdl%nat - ndof_r
         if (pub_on_root) write(stdout,*) 'done'

         ! vv: Convert user-defined atomic velocities to internal
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             '- MD -: Transforming atomic to internal velocities...'
         do ig = 1, group_num
            call md_atomic_to_internal_vel(ig,ndof_r,'forward',mdl%nat)
         end do
         if (pub_on_root) write(stdout,*) 'done'

      else if (md_global_restart) then

         ! Restore time, positions and velocities from
         ! previous calculation
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             '- MD -: Restore structure from previous restart...'
         call md_restore_trajectory(md_time,ierr1,mdl%nat)
         if (pub_on_root) then
            if (ierr1 ==0) then
               write(stdout,*) 'done'
            else
               write(stdout,*) 'failed'
            endif
         endif

         if (md_restart_thermo) then
            if (ierr1 == 0) then
               ! Restore thermostat parameters from previous
               ! calculation
               if (pub_on_root) write(stdout,'(a)',advance='no') &
                   '- MD -: Restore thermostats from previous calculation...'
               call md_restore_thermostat(ierr2)
               if (pub_on_root) then
                  if (ierr2 ==0) then
                     write(stdout,*) 'done'
                  else
                     write(stdout,*) 'failed'
                  endif
               endif

               if (ierr2 .ne. 0) then
                  ! Initialise thermostat parameters
                  if (pub_on_root) write(stdout,'(a)',advance='no') &
                      '- MD -: Initialising thermostats...'
                  call md_initialise_thermostats()
                  if (pub_on_root) write(stdout,*) 'done'
               endif
               ! vv: Initialise tgroup
               do ig=1, group_num
                  call md_select_thermostat(ith,ig,md_time)
                  tgroup(ig) = md_thermo(ith)%T
               end do
            end if
         else
            ! Initialise thermostat parameters
            if (pub_on_root) write(stdout,'(a)',advance='no') &
                '- MD -: Initialising thermostats...'
            call md_initialise_thermostats()
            ! vv: Initialise tgroup
            if (pub_on_root) write(stdout,*) 'done'
            do ig=1, group_num
               call md_select_thermostat(ith,ig,md_time)
               tgroup(ig) = md_thermo(ith)%T
            end do
         endif

         ! vv: Initialise number of degrees of freedom
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             '- MD -: Initialising degrees of freedom...'
         call md_initialise_ndof(ndof_r, md_time, mdl%nat)
         ndof = 3*mdl%nat - ndof_r
         if (pub_on_root) write(stdout,*) 'ndof =',ndof,'done'

         ! vv: Convert user-defined atomic velocities to internal
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             '- MD -: Transforming atomic to internal velocities...'
         do ig = 1, group_num
            call md_atomic_to_internal_vel(ig,ndof_r,'forward',mdl%nat)
         end do
         if (pub_on_root) write(stdout,*) 'done'

      end if


      if ((.not. md_restart .and. .not. md_global_restart) .or. (ierr1.ne.0)) then

         ! Initialise positions
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             '- MD -: Initialising positions...'
         call md_initialise_positions(mdl%elements,mdl%nat)
         if (pub_on_root) write(stdout,*) 'done'

         ! Initialise thermostat
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             '- MD -: Initialising thermostats...'
         call md_initialise_thermostats()
         if (pub_on_root) write(stdout,*) 'done'

         ! Initialise number of degrees of freedom
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             '- MD -: Initialising degrees of freedom...'
         call md_initialise_ndof(ndof_r, md_time, mdl%nat)
         ndof = 3*mdl%nat - ndof_r
         if (pub_on_root) write(stdout,*) 'done'

         ! Initialise velocities
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             '- MD -: Initialising velocities...'
         do ig = 1, group_num
            call md_select_thermostat(ith, ig, md_time)
            if (ith .gt. 0) then
               tgroup(ig) = md_thermo(ith)%T
               call md_initialise_velocities(ig, tgroup(ig), mdl%elements,ndof_r,mdl%nat)
            else
               call md_initialise_velocities(ig, 0.d0, mdl%elements,ndof_r,mdl%nat)
            endif
         enddo

      endif

      ! Adjust the local counter according to md_global_iter and user input
      call utils_assert(md_num_iter > md_iter_global,'Error: the total number of &
           &MD steps must be greater than the number of MD steps already done.')
      md_num_iter = md_num_iter - md_iter_global
      md_maxstep = int(md_num_iter/mts_nstep)

      !###############################################################
      !#### Set-up electronic history ################################

      ! Initialise the composition method for electronic degrees of
      ! freedom at subsequent time step
      if(.not. md_global_restart) then
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             '- MD -: Initialising main electronic history...'
         call elec_history_initialise_method('md',mix_dkn_type, &
                 mix_dkn_num, mix_dkn_reset, mix_dkn_init_type, &
                 mix_dkn_init_num, mix_ngwfs_type, mix_ngwfs_num, &
                 mix_ngwfs_reset, mix_ngwfs_init_type, &
                 mix_ngwfs_init_num, mix_ngwfs_coeff, &
                 mix_local_length, mix_local_smear)
         if (pub_on_root) write(stdout,*) 'done'

         if (mts_nstep .gt. 1) then

            ! Initialise the auxiliary composition method used to speed-up
            ! mts corrections to the ionic forces
            if (pub_on_root) write(stdout,'(a)',advance='no') &
                '- MD -: Initialising mts electronic history...'
            call elec_history_initialise_method('continue','REUSE',1,-2, &
               'NONE',0,'REUSE',1,-2,'NONE',0,0.d0,0.d0,0.d0)
            if (pub_on_root) write(stdout,*) 'done'

         endif

         ! Create storage for the electronic degrees of freedom to be
         ! kept in memory
         call elec_history_create_storage(mdl%nsub)
      else
         call elec_history_restore_history_info(mdl%nat,mdl%nsub)
      end if

      !###############################################################
      !#### Compute atomic forces (initial configuration) ############

        ! Update coordinates
        call md_update_coordinates(mdl)

      if(.not. md_global_restart ) then
        ! Write initial coordinates
        if (pub_write_positions) call md_write_coordinates(mdl%elements)

        ! Calculate energy and forces
        call md_energy_and_force(total_energy,forces,mdl,mts_step, &
                mts_nstep,md_iter)

        ! Calculate acceleration
        if (pub_on_root) write(stdout,*) 'Calculate accelerations... '
        do iat=1,mdl%nat
           a_t(:,iat) = forces(:,iat)/ion_mass(iat)
        enddo

        ! Update elements
        call md_update_coordinates(mdl)

        ! Write velocities
        if (pub_write_velocities) call md_write_velocities(mdl%elements)

        ! Calculate energies and temperature
        if (pub_on_root) write(stdout,*) 'Calculate energies and temperature... '
        md_potential_energy = total_energy
        md_kinetic_energy = kinetic_energy(mdl%nat)
        md_energy = md_kinetic_energy+md_potential_energy
        md_actual_temp = 2.0_dp*md_kinetic_energy/ndof

        ! Output MD info
        call md_write_md_info()

        ! Write initial state
        if (.not. md_restart ) then
           if (pub_on_root) write(stdout,'(a)',advance='no') &
                '- MD -: Writing initial state...'
           call md_write_trajectory(mdl%elements,mdl%cell)
           call md_write_thermostat(md_time)
           if (pub_on_root) write(stdout,*) 'done'
        endif

        ! Write backup files
        if (pub_on_root) write(stdout,'(a)',advance='no') &
             '- MD -: Writing MD restart files...'
        call md_backup_trajectory(mdl%nat)
        ! vv1c12: Skip writing a backup file if we are restarting a MD calculation
        ! and adding a thermostat
        if (md_restart_thermo) call md_backup_thermostat(backup_thermo=.true.)
        if (pub_on_root) write(stdout,*) 'done'

        ! Turn off reading flags
        pub_read_denskern = .false.
        pub_read_tightbox_ngwfs = .false.
        pub_read_sw_ngwfs = .false.
      end if

      !###############################################################
      !#### Compute trajectory (main MD loop) ########################

      md_loop : do md_step = 1, md_maxstep

         mts_loop : do mts_step = 1, mts_nstep


            ! Elapsed time
            md_time = md_time + md_delta_t
            md_iter = md_iter + 1
            md_iter_global = md_iter_global + 1

            if (pub_on_root) then
               write(stdout,1) 'Starting MD iteration :',&
                    &md_iter,'(',md_iter_global,')'
               if (mts_nstep .gt. 1) then
                  write(stdout,2) ' Multi-time step integration :', mts_step
               endif
            endif

            ! vv: if we switch on/off a thermostat we must re-evaluate the ndof
            !     and rescale the internal velocities
            call md_initialise_ndof(new_ndof_r,md_time,mdl%nat)
            if ((new_ndof_r /= ndof_r)) then
               ndof_r = new_ndof_r
               ndof = 3*mdl%nat - ndof_r
               do ig = 1, group_num
                 !tgroup(ig) = new_tgroup(ig)
                  call md_scale_velocities(ig,ndof,tgroup(ig))
               end do
            end if

            !if required, reset the electronic history
            if (modulo(md_iter,md_reset_history) == 0) &
               call elec_history_reset()

            ! Scale particle velocities,
            ! update thermostat velocities and positions
            if (mts_xi .or. mts_step == 1) &
               call md_velocity_verlet_thermostat(1,md_time)

            ! v(t+dt/2) = v(t) + 0.5 * dt * a(t)
            do iat=1,mdl%nat
               v_t(:,iat) = v_t(:,iat) + 0.5_dp * md_delta_t * a_t(:,iat)
            enddo

            ! Scale particle velocities,
            ! update thermostat velocities and positions
            if (mts_xi .or. mts_step == 1) &
               call md_velocity_verlet_thermostat(2,md_time)

            ! r(t+dt) = r(t) + dt * v(t+dt/2)
            call md_atomic_to_internal_vel(1,ndof_r,'backward',mdl%nat)
            do iat=1,mdl%nat
               r_t(:,iat) = r_t(:,iat) + md_delta_t * v0_t(:,iat)
            enddo

            ! Apply constraints
            ! call shake()

            ! Update elements
            call md_update_coordinates(mdl)

            ! Write coordinates
            if (pub_write_positions) call md_write_coordinates(mdl%elements)

            ! Calculate energy and forces
            call md_energy_and_force(total_energy,forces,mdl,mts_step,&
                 mts_nstep,md_iter)

            ! Calculate acceleration
            do iat=1,mdl%nat
               a_t(:,iat) = forces(:,iat)/ion_mass(iat)
            enddo

            ! Scale particle velocities,
            ! update thermostat velocities and positions
            if (mts_xi .or. mts_step == mts_nstep) &
               call md_velocity_verlet_thermostat(3,md_time)

            ! v(t+dt) = v(t+dt/2) + 0.5 * dt * a(t+dt)
            do iat=1,mdl%nat
               v_t(:,iat) = v_t(:,iat) + 0.5_dp * md_delta_t * a_t(:,iat)
            enddo
            do ig = 1, group_num
               call md_atomic_to_internal_vel(ig,ndof_r,'backward',mdl%nat)
            end do

            ! Apply constraints
            ! call shake()

            ! Scale particle velocities, update thermostat velocities and positions
            if (mts_xi .or. mts_step == mts_nstep) &
               call md_velocity_verlet_thermostat(4,md_time)

            ! Update elements
            call md_update_coordinates(mdl)

            ! Write velocities
            if (pub_write_velocities) call md_write_velocities(mdl%elements)

            ! vv: Scale velocities between two md steps if the temperature in
            ! the thermostat block has been changed
            do ig = 1, group_num
               call md_select_thermostat(ith, ig, md_time)
               if (ith .gt. 0) then
                  new_tgroup(ig) = md_thermo(ith)%T
                  if ((new_tgroup(ig) /= tgroup(ig))) then
                     tgroup(ig) = new_tgroup(ig)
                     call md_scale_velocities(ig,ndof,tgroup(ig))
                  end if
               end if
            end do

            ! Energies and temperature
            md_potential_energy = total_energy
            md_kinetic_energy = kinetic_energy(mdl%nat)
            md_energy = md_kinetic_energy + md_potential_energy
            md_actual_temp = 2.0_dp * md_kinetic_energy/ndof

            ! Output MD info
            call md_write_md_info()

            ! Write trajectory
            if (pub_on_root) write(stdout,'(a)',advance='no') &
                 ' Writing MD trajectory...'
            call md_write_trajectory(mdl%elements,mdl%cell)
            call md_write_thermostat(md_time)
            if (pub_on_root) write(stdout,*) 'done'

            ! Write backup files
            if (pub_on_root) write(stdout,'(a)',advance='no') &
                 ' Writing MD restart files...'
            call md_backup_trajectory(mdl%nat)
            ! vv: If this is the first iteration after restart and the
            ! thermostat block has changed over restarts than do not backup
            ! the old thermostat block
            if ((md_restart .or. md_global_restart) .and. (.not. md_restart_thermo) &
               .and. md_iter == 1) then
               call md_backup_thermostat(backup_thermo=.false.)
            else
               call md_backup_thermostat(backup_thermo=.true.)
            end if
            if (pub_on_root) write(stdout,*) 'done'

            ! vv: Write the global backup if necessary
            !if ((.not. md_restart) .and. (md_write_history > 0) &
            if ((md_write_history > 0) &
                .and. (md_iter_global >= md_write_history) &
                .and. mod(md_iter_global,md_write_history)==0) then
               write_global = .true.
               if (pub_on_root) write(stdout,'(a)',advance='no') &
                    ' Writing MD trajectory to global store...'
               call md_write_trajectory(mdl%elements,mdl%cell,write_global)
               call md_write_thermostat(md_time,write_global)
               call md_backup_trajectory(mdl%nat, write_global)
               call md_backup_thermostat(write_global,backup_thermo=.true.)
               if (pub_on_root) write(stdout,*) 'done'
               write_global = .false.
            end if

         enddo mts_loop

      enddo md_loop

      ! Deallocate extrapolation arrays
      call elec_history_destroy_storage()
      deallocate(tgroup,stat=ierr)
      call utils_dealloc_check('md_velocity_verlet','tgroup',ierr)
      deallocate(new_tgroup,stat=ierr)
      call utils_dealloc_check('md_velocity_verlet','new_tgroup',ierr)

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_velocity_verlet'
      return

1     format(/,80('x'),/,11('x'),a37,1x,i9,a1,i9,a1,11('x'),/,80('x'),/)
2     format(/,10x,60('x'),/,10x,'x',a37,1x,i9,11x,'x',/,10x,60('x'),/)

    end subroutine md_velocity_verlet

!==============================================================================!
!==============================================================================!

    subroutine md_energy_and_force(total_energy,forces,mdl, &
                        istep,nstep,iter)

    !====================================================================!
    ! Buffer layer between MD integration subroutines and                !
    ! energy_and_force_calculate subroutine.                             !
    !                                                                    !
    !--------------------------------------------------------------------!
    ! Implemented by Simon M.-M. Dubois (May 2011)                       !
    !====================================================================!

      use comms,            only: comms_barrier
      use constants,        only: DP, stdout
      use electronic_history, only: elec_history_create_entry, &
        elec_history_use_method, elec_history_update_method
      use energy_and_force, only: energy_and_force_calculate
      use model_type,       only: MODEL
      use md_ions,          only: r_t
      use rundat,           only: mts_mix_inc, mts_nstep, md_global_restart
      use utils,            only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Argumentas
      type(MODEL), intent(inout)     :: mdl
      real(kind=DP), intent(out)     :: total_energy
      real(kind=DP), intent(out)     :: forces(1:3,1:mdl%nat)
      integer, intent(in)            :: istep
      integer, intent(in)            :: nstep
      integer, intent(in)            :: iter

      ! Variables
      real(kind=DP), allocatable :: xforces(:,:)
      real(kind=DP) :: deltaf(3)
      integer       :: iat, ierr

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_energy_and_force'

      if (md_global_restart .and. iter == 1) then

      else
      ! Store particle coordinates (required for extrapolation
      ! of NGWFS/denskern at subsequent MD steps)
         call elec_history_create_entry(r_t)
      end if

      ! Chose the composition method for initial electronic
      ! degrees of freedom
      if (nstep .le. 1 .or. istep == 0) then
         call elec_history_use_method('md')
         call activate_md_param()
      else
         call elec_history_use_method('md')
         call activate_mts_param()
      endif

      ! Compute ionic forces and total energy
      call comms_barrier()
      call energy_and_force_calculate(total_energy,forces,mdl)

      ! Update the history of electronic
      ! degrees of freedom
      if (nstep .le. 1) then
         call elec_history_update_method('md')
      else
         if (istep == 0) then
            call elec_history_update_method('md')
         else if (istep < nstep) then
            call elec_history_update_method('md')
         else if (istep == nstep .and. .not. mts_mix_inc) then
            call elec_history_update_method('md')
            call elec_history_update_method('continue')
            call elec_history_create_entry(r_t)
            call elec_history_use_method('continue')
         else if (istep == nstep .and. mts_mix_inc) then
            call elec_history_update_method('continue')
            call elec_history_use_method('continue')
         endif
      endif

      ! Apply mts corrections to the forces
      if (nstep .gt. 1 .and. istep == nstep) then

         call activate_md_param()

         allocate(xforces(3,mdl%nat),stat=ierr)
         call utils_alloc_check('md_energy_and_force','xforces',ierr)

         call comms_barrier()
         call energy_and_force_calculate(total_energy,xforces,mdl)

         do iat = 1,mdl%nat
            deltaf(:) = mts_nstep*(xforces(:,iat) - forces(:,iat))
            forces(:,iat) = forces(:,iat) + deltaf(:)
         enddo

         deallocate(xforces,stat=ierr)
         call utils_dealloc_check('md_energy_and_force','xforces',ierr)

         if (mts_mix_inc) then
            call elec_history_update_method('md')
         endif

      endif

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_energy_and_force'

      return

    end subroutine md_energy_and_force

!==============================================================================!
!==============================================================================!


    subroutine md_update_coordinates(mdl)

      ! Modified to also update embedding structures by Joseph
      ! Prentice, May 2018

      use model_type,      only: MODEL
      use md_ions,         only: r_t, v0_t
      use services,        only: services_rationalise_coords

      implicit none

      ! Arguments
      type(MODEL), intent(inout) :: mdl

      ! Local variables
      integer :: iat,jat(mdl%nsub),iregion

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_update_coordinates'

      call services_rationalise_coords(size(mdl%elements),r_t,mdl%cell)
      do iat=1,size(mdl%elements)
         mdl%elements(iat)%centre%x =  r_t(1,iat)
         mdl%elements(iat)%centre%y =  r_t(2,iat)
         mdl%elements(iat)%centre%z =  r_t(3,iat)
         mdl%elements(iat)%ion_velocity(:) =  v0_t(:,iat)
      enddo

      ! jcap: copy the data into the appropriate subregions
      jat = 0
      do iat=1,size(mdl%elements)
         ! Get the region counter
         iregion = mdl%elements(iat)%region
         ! Count the no. of atoms allocated to each region
         jat(iregion) = jat(iregion) + 1
         ! Copy the data over
         mdl%regions(iregion)%elements(jat(iregion))%centre = &
              &mdl%elements(iat)%centre
         mdl%regions(iregion)%elements(jat(iregion))%ion_velocity = &
              &mdl%elements(iat)%ion_velocity
      end do

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_update_coordinates'

      return

    end subroutine md_update_coordinates


!==============================================================================!
!==============================================================================!


    subroutine md_write_coordinates(elements)


      use comms,           only: pub_on_root
      use constants,       only: stdout
      use ion,             only: element

      implicit none

      type(element), intent(in) :: elements(:)

      ! Local variables
      integer :: iat

      if (pub_on_root) then
         write(stdout,*)
         write(stdout,*)'                        -------------------------------'
         write(stdout,*)'                                  Cell Contents        '
         write(stdout,*)'                        -------------------------------'
         write(stdout,*)
         write(stdout,21)
         write(stdout,*)'        x  Element    Atom    Absolute co-ordinates of atoms (bohr)  x'
         write(stdout,*)'        x            Number         x           y           z        x'
         write(stdout,22)
         do iat=1,size(elements)
            write(stdout,20) elements(iat)%symbol, iat, elements(iat)%centre%x, &
                             elements(iat)%centre%y, elements(iat)%centre%z
         end do
         write(stdout,21)
         write(stdout,*)
      endif

20    format(8x,' x',a6,1x,i8,5x,3f12.6,'    x')
21    format(1x,'        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
22    format(1x,'        x------------------------------------------------------------x')

      return

    end subroutine md_write_coordinates


!==============================================================================!
!==============================================================================!


    subroutine md_write_velocities(elements)


      use comms,           only: pub_on_root
      use constants,       only: stdout
      use esdf,            only: esdf_convfac
      use ion,             only: element

      implicit none

      type(element), intent(in) :: elements(:)

      ! Local variables
      integer :: iat
      real(kind=DP)  :: convfac


      if (pub_on_root) then
         convfac = esdf_convfac('auv','bohr/fs')
         write(stdout,*)
         write(stdout,21)
         write(stdout,*)'        x  Element    Atom    Absolute velocities of atoms (bohr/fs) x'
         write(stdout,*)'        x            Number         x           y           z        x'
         write(stdout,22)
         do iat=1,size(elements)
            write(stdout,20) elements(iat)%symbol, iat, elements(iat)%ion_velocity(1), &
                             elements(iat)%ion_velocity(2), elements(iat)%ion_velocity(3)
         end do
         write(stdout,21)
         write(stdout,*)
      endif

20    format(8x,' x',a6,1x,i8,5x,3f12.6,'    x')
21    format(1x,'        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
22    format(1x,'        x------------------------------------------------------------x')

      return

    end subroutine md_write_velocities


!==============================================================================!
!==============================================================================!


    subroutine md_write_md_info()


      use comms,     only: pub_on_root
      use constants, only: stdout

      implicit none

      if (pub_on_root) then
         write(stdout,*)
         write(stdout,21)

         write(stdout,1)'MD INFO                        '
         write(stdout,1)' '
         write(stdout,2)'            time:',md_time*aut2fs,trim(time_label)
         write(stdout,1)' '
         write(stdout,2)'Potential Energy:',md_potential_energy,trim(energy_label)
         write(stdout,2)'Kinetic   Energy:',md_kinetic_energy,trim(energy_label)
         write(stdout,2)'Total     Energy:',md_energy,trim(energy_label)
         write(stdout,1)' '
         write(stdout,2)'     Temperature:',md_actual_temp*ha2k,trim(temperature_label)
         write(stdout,1)' '
         write(stdout,21)
         write(stdout,*)
      endif

1     format(11('x'),a58                ,11('x'))
2     format(11('x'),a20,f14.6,1x,a20,3x,11('x'))
21    format(80('x'))

      return

    end subroutine md_write_md_info

!==============================================================================!
!==============================================================================!

    subroutine md_backup_trajectory(nat,global)

      use constants,       only: stdout
      use comms,           only: pub_on_root
      use ion,             only: element
      use md_ions,         only: r_t, a_t, v0_t, vtot, Lcm
      use rundat,          only: pub_rootname
      use utils,           only: utils_unit, utils_open_unit_check, &
             utils_close_unit_check, utils_alloc_check, utils_dealloc_check

      implicit none

      ! Arguments
      integer           :: nat
      logical, optional :: global

      ! Local variables
      integer           :: io_unit, io_stat
      integer           :: iat,idir, iob
      character(len=80) :: filename
      real(kind=dp),allocatable :: buffer(:)
      integer           :: buf_length, ierr
      integer           :: ibuffer
      logical           :: fileexists

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_backup_trajectory'

      ! Backup the old .md.restart file to prevent any loss of data
      if(.not. present(global)) then
         ! Write .md.restart
         write(filename,'(a,a)') trim(pub_rootname),'.md.restart'
         buf_length = 6*nat+7
      elseif (present(global)) then
         if(global) then
            ! Write .md.global.restart
            write(filename,'(a,a)') trim(pub_rootname),'.md.global.restart'
            buf_length = 9*nat+7
         end if
      end if

      allocate(buffer(buf_length),stat=ierr)
      call utils_alloc_check('md_backup_trajectory','buffer',ierr)

      if (pub_on_root) then

         inquire(file=filename,exist=fileexists)

         if (fileexists) then

            io_unit = utils_unit()

            open(unit=io_unit,iostat=io_stat,file=filename, &
                 access='SEQUENTIAL',form='UNFORMATTED',position='REWIND', &
                 action='READ')
            call utils_open_unit_check('md_backup',filename,io_stat)
            ! Read md_iter_global as an integer
            read(io_unit) ibuffer
            ! Read all the rest as a real double precision
            if (present(global)) then
               if (global) then
                  do iob = 1, 9*nat+7
                     read(io_unit) buffer(iob)
                  enddo
               end if
            else
               do iob = 1, 6*nat+7
                  read(io_unit) buffer(iob)
               enddo
            end if
            close(unit=io_unit,iostat=io_stat)
            call utils_close_unit_check('md_backup',filename,io_stat)

            ! Write content in .md.restart.old
            io_unit = utils_unit()
            write(filename,'(a,a)') trim(pub_rootname),'.md.restart.old'
            open(unit=io_unit,iostat=io_stat,file=filename,status='REPLACE',&
                 form='UNFORMATTED',position='REWIND',action='WRITE')
            call utils_open_unit_check('md_backup',filename,io_stat)

            write(io_unit) ibuffer

            do iob = 1, 6*nat+7
               write(io_unit) buffer(iob)
            enddo
            close(io_unit,iostat=io_stat)
            call utils_close_unit_check('md_backup',filename,io_stat)
         endif
      endif

      ! Write a new .md.restart file
      if (pub_on_root) then

         io_unit = utils_unit()
         if(.not. present(global)) then
            ! Read .md.restart
            write(filename,'(a,a)') trim(pub_rootname),'.md.restart'
         elseif (present(global)) then
            if(global) then
               ! Read .md.global.restart
               write(filename,'(a,a)') trim(pub_rootname),'.md.global.restart'
            end if
         end if
         open(unit=io_unit,iostat=io_stat,file=filename,status='REPLACE',&
              form='UNFORMATTED',position='REWIND',action='WRITE')
         call utils_open_unit_check('md_backup',filename,io_stat)

         write(io_unit) md_iter_global

         write(io_unit) md_time

         ! vv: write the velocity of com
         do idir =1,3
            write(io_unit) vtot(idir)
         end do

         ! vv: write the angular momentum of com
         do idir =1,3
            write(io_unit) Lcm(idir)
         end do

         do iat=1,nat
            do idir =1,3
               write(io_unit) r_t(idir,iat)
            enddo
         enddo

         do iat=1,nat
            do idir =1,3
               write(io_unit) v0_t(idir,iat)
            enddo
         enddo

         if (present(global)) then
            if (global) then
               do iat=1,nat
                  do idir=1,3
                     write(io_unit) a_t(idir,iat)
                  end do
               end do
            end if
         end if

         close(io_unit,iostat=io_stat)
         call utils_close_unit_check('md_backup',filename,io_stat)

      endif

      deallocate(buffer,stat=ierr)
      call utils_dealloc_check('md_backup_trajectory','buffer',ierr)

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_backup_trajectory'

      return

    end subroutine md_backup_trajectory

!==============================================================================!
!==============================================================================!

    subroutine md_write_trajectory(elements,cell,global)


      use comms,           only: pub_on_root
      use ion,             only: element
      use md_ions,         only: r_t, v0_t, a_t, ion_mass, &
             md_ionic_polarisation, md_elec_polarisation, &
             md_total_polarisation
      use rundat,          only: pub_rootname
      use simulation_cell, only: CELL_INFO
      use utils,           only: utils_unit, utils_open_unit_check, &
             utils_close_unit_check

      implicit none

      ! Arguments
      type(ELEMENT),    intent(in) :: elements(:)
      type(CELL_INFO),  intent(in) :: cell
      logical, optional,intent(in) :: global

      ! Local variables
      integer           :: io_unit,io_stat
      integer           :: ispec,iat
      character(len=12) :: action,form,stat,position,access
      character(len=80) :: filename

      action   = 'WRITE'
      form     = 'FORMATTED'
      stat     = 'UNKNOWN'
      position = 'REWIND'
      access   = 'SEQUENTIAL'

      if (md_time==0.0_dp) then
         stat     ='REPLACE'        !create a new output file
      else
         position ='APPEND'         !append to existing output file
      end if


      if (pub_on_root) then

         ! Find available unit specifier
         io_unit = utils_unit()
         if(.not. present(global)) then
            ! Write to <filename>.md
            write(filename,'(a,a)') trim(pub_rootname),'.md'
         elseif (present(global)) then
            if(global) then
               ! Write to <filename>.md.global
               write(filename,'(a,a)') trim(pub_rootname),'.md.global'
            end if
         end if

         open(unit=io_unit,iostat=io_stat,file=filename,status=stat,&
              access=access,form=form,position=position,action=action)
         call utils_open_unit_check('md_write_trajectory',filename,io_stat)

         ! output is in atomic units except md_time (in fs) and temperature (in kelvin)
         write(io_unit,1) md_time*aut2fs
         write(io_unit,2) md_potential_energy,md_energy,md_kinetic_energy
         write(io_unit,3) md_actual_temp*ha2k

         ! cell vectors
         write(io_unit,4) cell%a1%x, cell%a1%y, cell%a1%z
         write(io_unit,4) cell%a2%x, cell%a2%y, cell%a2%z
         write(io_unit,4) cell%a3%x, cell%a3%y, cell%a3%z

         ispec=1

         ! ionic positions
         do iat=1,size(elements)
            write(io_unit,7) elements(iat)%symbol,ispec,r_t(:,iat)
         enddo

         ! ionic velocities
         do iat=1,size(elements)
            write(io_unit,8) elements(iat)%symbol,ispec,v0_t(:,iat)
         end do

         ! ionic forces
         do iat=1,size(elements)
            write(io_unit,9) elements(iat)%symbol,ispec,ion_mass(iat)*a_t(:,iat)
         end do

         ! Polarisation
         write(io_unit,10) md_elec_polarisation(:)
         write(io_unit,10) md_ionic_polarisation(:)
         write(io_unit,10) md_total_polarisation(:)

         ! blank line to signal end of iteration
         write(io_unit,*) ' '

         close(io_unit,iostat=io_stat)
         call utils_close_unit_check('md_write_trajectory',filename,io_stat)

      endif

1     format(12x,es18.8e3)
2     format(9x,3(3x,es18.8e3),'  <-- E')
3     format(12x,es18.8,T73,   '  <-- T')
4     format(9x,3(3x,es18.8e3),'  <-- h')
7     format(1x,a3,1x,i4,3(3x,es18.8e3),'  <-- R')
8     format(1x,a3,1x,i4,3(3x,es18.8e3),'  <-- V')
9     format(1x,a3,1x,i4,3(3x,es18.8e3),'  <-- F')
10    format(9x,3(3x,es18.8e3),'  <-- P')

      return

    end subroutine md_write_trajectory

!==============================================================================!
!==============================================================================!

    subroutine init_mts_scheme()

      use rundat,    only: mts_nstep, &
           ngwf_threshold_orig, pub_elec_energy_tol, &
           pub_elec_force_tol, pub_ngwf_max_grad, lnv_threshold_orig, &
           minit_lnv, maxit_lnv, pub_maxit_pen, pub_maxit_ngwf_cg, &
           mts_ngwf_threshold, mts_elec_energy_tol, mts_maxit_ngwf_cg, &
           mts_elec_force_tol, mts_ngwf_max_grad, mts_lnv_threshold, &
           mts_minit_lnv, mts_maxit_lnv, mts_maxit_pen, &
           mermin_threshold_orig, pub_mermin

      implicit none

      main_maxit_ngwf_cg    = pub_maxit_ngwf_cg
      main_ngwfs_threshold  = ngwf_threshold_orig
      !ep: mermin switch
      if (pub_mermin) then
         main_lnv_threshold    = mermin_threshold_orig
      else
         main_lnv_threshold    = lnv_threshold_orig
      end if
      !ep: mermin switch
      main_elec_energy_tol  = pub_elec_energy_tol
      main_elec_force_tol   = pub_elec_force_tol
      main_ngwf_max_grad    = pub_ngwf_max_grad
      main_minit_lnv        = minit_lnv
      main_maxit_lnv        = maxit_lnv
      main_maxit_pen        = pub_maxit_pen

      if (mts_nstep .le. 1) then
         mts_maxit_ngwf_cg    = pub_maxit_ngwf_cg
         mts_ngwf_threshold   = ngwf_threshold_orig
         !ep: mermin switch
         if (pub_mermin) then
            mts_lnv_threshold    = mermin_threshold_orig
         else
            mts_lnv_threshold    = lnv_threshold_orig
         end if
         !ep: mermin switch
         mts_elec_energy_tol  = pub_elec_energy_tol
         mts_elec_force_tol   = pub_elec_force_tol
         mts_ngwf_max_grad    = pub_ngwf_max_grad
         mts_minit_lnv        = minit_lnv
         mts_maxit_lnv        = maxit_lnv
         mts_maxit_pen        = pub_maxit_pen
      endif

      return

    end subroutine init_mts_scheme

!==============================================================================!
!==============================================================================!

    subroutine activate_mts_param()

      use rundat,  only: ngwf_threshold_orig, pub_elec_energy_tol, &
            pub_elec_force_tol, pub_ngwf_max_grad, lnv_threshold_orig, &
            minit_lnv, maxit_lnv, pub_maxit_pen, mts_maxit_ngwf_cg, &
            mts_ngwf_threshold, mts_elec_energy_tol, pub_maxit_ngwf_cg, &
            mts_elec_force_tol, mts_ngwf_max_grad, mts_lnv_threshold, &
            mts_minit_lnv, mts_maxit_lnv, mts_maxit_pen, &
            mermin_threshold_orig, pub_mermin

      implicit none

      pub_maxit_ngwf_cg = mts_maxit_ngwf_cg
      ngwf_threshold_orig = mts_ngwf_threshold
      !ep: mermin switch
      if (pub_mermin) then
          mermin_threshold_orig  = mts_lnv_threshold
      else
          lnv_threshold_orig  = mts_lnv_threshold
      end if
      !ep: mermin switch
      pub_elec_energy_tol = mts_elec_energy_tol
      pub_elec_force_tol  = mts_elec_force_tol
      pub_ngwf_max_grad   = mts_ngwf_max_grad
      minit_lnv   =  mts_minit_lnv
      maxit_lnv   =  mts_maxit_lnv
      pub_maxit_pen   =  mts_maxit_pen

    end subroutine activate_mts_param

!==============================================================================!
!==============================================================================!

    subroutine activate_md_param()

      use rundat,  only: ngwf_threshold_orig, lnv_threshold_orig, &
            pub_elec_energy_tol, pub_elec_force_tol, pub_ngwf_max_grad, &
            minit_lnv, maxit_lnv, pub_maxit_pen, pub_maxit_ngwf_cg, &
            mermin_threshold_orig, pub_mermin

      implicit none

      pub_maxit_ngwf_cg   = main_maxit_ngwf_cg
      ngwf_threshold_orig = main_ngwfs_threshold
      !ep: mermin switch
      if (pub_mermin) then
          mermin_threshold_orig  = main_lnv_threshold
      else
          lnv_threshold_orig  = main_lnv_threshold
      end if
      !ep: mermin switch
      pub_elec_energy_tol = main_elec_energy_tol
      pub_elec_force_tol  = main_elec_force_tol
      pub_ngwf_max_grad   = main_ngwf_max_grad
      minit_lnv   = main_minit_lnv
      maxit_lnv   = main_maxit_lnv
      pub_maxit_pen   = main_maxit_pen

    end subroutine activate_md_param

!==============================================================================!
!==============================================================================!

end module md
