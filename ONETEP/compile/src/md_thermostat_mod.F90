! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*- !

!============================================================================!
!                                                                            !
!                 Molcular Dynamics Thermostat Module                        !
!                                                                            !
! Thermostat subroutines for a system of classical particles coupled to an   !
! thermal bath. Includes ANDERSON/LANGEVIN/NOSE-HOOVER/NOSE-HOOVER-CHAIN/    !
! BERENDSEN/BUSSI schemes for the description of the coupling.               !
!                                                                            !
! Calls to thermostat subroutines should be interleaved with the velocity    !
! Verlet algorithm as follow :                                               !
!                                                                            !
!  thermo1                                                                   !
!     v(t+dt/2) = v(t) + a(t)*dt/2                                           !
!  thermo2                                                                   !
!     r(t+dt) = r(t) + v(t+dt/2)*dt                                          !
!     call energy_and_force_calculate()                                      !
!  thermo3                                                                   !
!     v(t+dt) = v(t+dt/2) + a(t+dt)*dt/2                                     !
!  thermo4                                                                   !
!                                                                            !
!----------------------------------------------------------------------------!
! Written by Simon M.-M. Dubois (May 2011)                                   !
!   Modified by S. M.-M. Dubois (Jul 2013)                                   !
!============================================================================!



module md_thermostat

  use constants, only: dp, stdout
  use rundat, only: pub_debug_on_root

  implicit none

  private

  ! Type definition for the thermostat.
  type, public :: thermostat

     ! General parameters
     integer         :: type        ! 0 = NONE, 1 = ANDERSEN,
                                    ! 2 = LANGEVIN, 3 = NOSE_HOOVER
                                    ! 4 = BERENDSEN, 5 = BUSSI
     real(kind=DP)   :: T           ! Thermostat temperature
     real(kind=DP)   :: Tinit       ! Thermostat temperature
     real(kind=DP)   :: Tgrad       ! Temperature gradient
     integer         :: group       ! Subset of atoms to which the thermostat
                                    ! is linked
     integer         :: start  ! Time window (expressed in md steps) in
     integer         :: stop   ! which the thermostat is defined

     ! Andersen :
     real(kind=DP)   :: mix         ! Mixing coeff for a softer rescaling

     ! Andersen and Bussi :
     real(kind=DP)   :: tau         ! Characteristic time scale of the stochastic

     ! Langevin :
     real(kind=DP)   :: damp        ! Langevin damping parameter

     ! Nose-Hoover chain :
     integer         :: nhc_length         ! #-elements in the NH chain
     integer         :: nhc_integ_nstep    ! #-steps to integrate NHC dynamics
     logical         :: nhc_upd
     real(kind=DP), pointer :: nhc_x(:)    ! Thermostat positions
     real(kind=DP), pointer :: nhc_v(:)    ! Thermostat velocities
     real(kind=DP), pointer :: nhc_g(:)    ! Thermostat accelerations
     real(kind=DP), pointer :: nhc_mass(:) ! Thermostat masses

  end type thermostat

  ! Module variables
  integer, public :: md_thermo_num
  type(thermostat), public, allocatable :: md_thermo(:)

  ! Module parameters
  real(kind=DP), parameter :: aut2fs = 0.0241888468_dp ! atomic unit of time --> fs

  ! List of subroutines
  public  :: md_velocity_verlet_thermostat
  public  :: md_select_thermostat
  public  :: md_initialise_thermostats
  public  :: md_restore_thermostat
  public  :: md_backup_thermostat
  public  :: md_write_thermostat
  public  :: md_initialise_ndof
!  private  :: md_update_thermostat
!  private  :: md_integrate_nhc

!  ^^^^^^^
! jd (retreat2014): The subroutines that are private need not be listed here,
!                   since they can't be accessed from outside this module
!                   anyway, and everything is private by default owing to the
!                   'private' clause at line 37. Normally it would not matter,
!                   but gfortran bug 54224 flags such functions as unused and
!                   generates spurious warnings. So for now I commented them out.
!                   cf. https://gcc.gnu.org/bugzilla/show_bug.cgi?id=54224

  contains


!=============================================================================!
!=============================================================================!

    subroutine md_velocity_verlet_thermostat(pos,time)

      !==============================================================!
      ! Modify the ionic forces and velocities in order to sample    !
      ! the canonical ensemble.                                      !
      ! Calls to this subroutine are interlocked with the velocity   !
      ! verlet integration of the classical equations of motion :    !
      !                                                              !
      ! [1] apply LANGEVIN, NOSEHOOVER                               !
      !         v(t+dt/2) = v(t) + a(t)*dt/2                         !
      ! [2] apply ANDERSEN, BERENDSEN, BUSSI                         !
      !         r(t+dt) = r(t) + v(t+dt/2)*dt                        !
      !         call energy_and_force_calculate()                    !
      ! [3] apply LANGEVIN                                           !
      !         v(t+dt) = v(t+dt/2) + a(t+dt)*dt/2                   !
      ! [4] apply ANDERSEN, NOSEHOOVER, BERENDSEN, BUSSI             !
      !                                                              !
      !--------------------------------------------------------------!
      ! Written by Simon M.-M. Dubois (May 2011)                     !
      !==============================================================!

      use constants,       only: DP, stdout
      use md_ions,         only: group_num
      use rundat,          only: md_delta_t,mts_nstep,mts_xi

      implicit none

      ! Argument
      integer, intent(in) :: pos
      real(kind=DP), intent(in) :: time

      ! Variables
      real(kind=DP) :: dt
      integer       :: thermo, group


      !==============================================================!

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_velocity_verlet_thermostat'

      if (mts_xi) then
         dt = md_delta_t
      else
         dt = md_delta_t*mts_nstep
      endif

      ! For each group of atoms, select the thermostat
      ! and apply appropriate rescaling
      do group = 1, group_num

         call md_select_thermostat(thermo,group,time)
         if (thermo .gt. 0) then
            call md_update_thermostat(thermo,group,time)

            ! Thermostat = NONE
            if (md_thermo(thermo)%type == 0) return

            ! Thermostat = ANDERSEN
            if (md_thermo(thermo)%type == 1) then

               if (pos == 2 .or. pos == 4) &
                  call md_andersen_scaling(thermo,group,dt)

            ! Thermostat = LANGEVIN
            elseif (md_thermo(thermo)%type == 2) then

               if (pos == 1 .or. pos == 3) &
                  call md_langevin_dynamics(thermo,group,dt)

            ! Thermostat = NOSE-HOOVER
            elseif (md_thermo(thermo)%type == 3) then

               if (pos == 1 .or. pos == 4) &
                  call md_integrate_nhc(thermo,group,dt)

            ! Thermostat = BERENDSEN
            elseif (md_thermo(thermo)%type == 4) then

               if (pos == 2 .or. pos == 4) &
                  call md_berendsen_scaling(thermo,group,dt)

            ! Thermostat = BUSSI
            elseif (md_thermo(thermo)%type == 5) then

               if (pos == 2 .or. pos == 4) &
                  call md_canonical_scaling(thermo,group,dt)

            endif
         endif
      enddo

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_velocity_verlet_thermostat'

      return


    end subroutine md_velocity_verlet_thermostat

!==============================================================================!
!==============================================================================!

    subroutine md_canonical_scaling(thermo,group,dt)

      !================================================================!
      ! Velocity rescaling for canonical sampling. For details, see:   !
      !    G. Bussi, D. Donadio, M. Parrinello,                        !
      !    J. Chem. Phys. 126, 014101 (2007)                           !
      !----------------------------------------------------------------!
      ! Implemented by Simon M.-M. Dubois (May 2013) on the basis of   !
      ! a routine written by Gilberto Teobaldi.                        !
      !================================================================!

      use comms,           only: pub_on_root, pub_root_proc_id, comms_bcast
      use constants,       only: dp, stdout
      use md_ions,         only: v_t, groups, group_kinetic_energy
      use services,        only: services_maxboltzdist, services_gammadist

      implicit none

      ! Argument
      integer, intent(in) :: thermo
      integer, intent(in) :: group
      real(kind=DP), intent(in) :: dt

      ! Variables
      real(kind=DP) :: dt2, edt
      real(kind=DP) :: R1, Ri2 !, Ri2_tmp
      real(kind=DP) :: alpha, alpha2
      real(kind=DP) :: ekin, nrkin
      integer  :: iat, jat, idir, ndof !, idof

      !====================================================================!

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering md_canonical_scaling'

      ! Propagate over Delta_t/2
      dt2 = 0.5_dp*dt

      ! Random number from Gaussian distribution
      R1 = services_maxboltzdist()

      ! Degrees of freedom
      ndof = groups(group)%num_atom*3

      ! Evaluate the sum of ndof (random_number)^2 from the corresponding
      ! Gamma distribution
      Ri2 = services_gammadist(ndof-1)
      Ri2 = Ri2 + R1*R1

      ! An alternative is to calculate Ri2 explicitly (no gamma-distribution)
      ! Ri2_tmp = 0.0_dp
      ! do idof = 2, ndof
      !   Ri = services_maxboltzdist()
      !   Ri2_tmp = Ri2_tmp + Ri*Ri
      ! enddo
      ! Ri2_tmp = Ri2_tmp + R1*R1


      ! Compute the scaling factor
      edt = exp(-dt2/md_thermo(thermo)%tau)
      ekin = group_kinetic_energy(group)
      nrkin = md_thermo(thermo)%T/(2.0_dp*ekin)

      alpha2 = edt + nrkin*(1.0_dp-edt)*Ri2 +  &
                2.0_dp*R1*SQRT(nrkin*edt*(1.0_dp-edt))
      alpha = SQRT(alpha2)
      if (pub_on_root) write(stdout,'(a,e12.5,1x,e12.5)') &
         " Canonical scaling : alpha2, alpha = ", alpha2, alpha

      ! Scale velocities
      if (pub_on_root) then
         do iat = 1, groups(group)%num_atom
            jat = groups(group)%list_atom(iat)
            do idir = 1,3
               v_t(idir,jat) = v_t(idir,jat)*alpha
            enddo
         enddo
      endif

      ! Make sure all procs have the same velocity and
      ! acceleration
      call comms_bcast(pub_root_proc_id,v_t(:,:))

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving md_canonical_scaling'

      return

    end subroutine md_canonical_scaling


!==============================================================================!
!==============================================================================!

    subroutine md_berendsen_scaling(thermo,group,dt)

      !================================================================!
      ! Berendsen thermostat. Similar implementation as for the        !
      ! canonical scaling (J.Chem.Phys. 126, 014101, 2007) but without !
      ! the stochastic terms.                                          !
      !----------------------------------------------------------------!
      ! Implemented by Simon M.-M. Dubois (May 2013)                   !
      !================================================================!

      use comms,           only: pub_on_root, pub_root_proc_id, comms_bcast
      use constants,       only: dp, stdout
      use md_ions,         only: v_t, groups, group_kinetic_energy
      use services,        only: services_maxboltzdist, services_gammadist

      implicit none

      ! Argument
      integer, intent(in) :: thermo
      integer, intent(in) :: group
      real(kind=DP), intent(in) :: dt

      ! Variables
      real(kind=DP) :: dt2, edt
      real(kind=DP) :: alpha, alpha2
      real(kind=DP) :: ekin, rkin
      integer  :: iat, jat, idir, ndof

      !====================================================================!

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering md_berendsen_scaling'

      ! Propagate over Delta_t/2
      dt2 = 0.5_dp*dt

      ! Degrees of freedom
      ndof = groups(group)%num_atom*3

      ! Compute the scaling factor
      edt = exp(-dt2/md_thermo(thermo)%tau)
      ekin = group_kinetic_energy(group)
      rkin = ndof*md_thermo(thermo)%T/(2.0_dp*ekin)

      alpha2 = edt + rkin*(1.0_dp-edt)
      alpha = SQRT(alpha2)
      if (pub_on_root) write(stdout,'(a,e12.5,1x,e12.5)') &
         " Berendsen scaling : alpha2, alpha = ", alpha2, alpha

      ! Scale velocities
      if (pub_on_root) then
         do iat = 1, groups(group)%num_atom
            jat = groups(group)%list_atom(iat)
            do idir = 1,3
               v_t(idir,jat) = v_t(idir,jat)*alpha
            enddo
         enddo
      endif

      ! Make sure all procs have the same velocity and
      ! acceleration
      call comms_bcast(pub_root_proc_id,v_t(:,:))

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving md_berendsen_scaling'

      return

    end subroutine md_berendsen_scaling


!==============================================================================!
!==============================================================================!

    subroutine md_langevin_dynamics(thermo,group,dt)

      !================================================================!
      ! Scale velocities and accelerations, and add randomly           !
      ! fluctuating accelerations according to Langevin dynamics       !
      ! equations.                                                     !
      !                                                                !
      !----------------------------------------------------------------!
      ! Implemented by Simon M.-M. Dubois (May 2011)                   !
      !================================================================!

      use comms,           only: pub_root_proc_id, comms_bcast
      use constants,       only: dp
      use md_ions,         only: v_t,a_t,groups,ion_mass
      use services,        only: services_maxboltzdist

      implicit none

      ! Argument
      integer, intent(in) :: thermo
      integer, intent(in) :: group
      real(kind=DP), intent(in) :: dt

      ! Variables
      real(kind=DP) :: dt2
      real(kind=DP) :: scale_v, scale_a
      real(kind=DP) :: scale_random, random_a
      integer  :: iat, jat, idir

      !====================================================================!

      ! Scale velocities and accelerations according to damping factor
      dt2 = 0.5_dp*dt
      scale_v = exp(-md_thermo(thermo)%damp*dt2)
      scale_a = (1-scale_v)/(md_thermo(thermo)%damp*dt2)

      do iat = 1, groups(group)%num_atom
         jat = groups(group)%list_atom(iat)

         v_t(:,jat) = scale_v*v_t(:,jat)
         a_t(:,jat) = scale_a*a_t(:,jat)
      enddo


      ! Add random fluctuating accelerations.
      scale_random = md_thermo(thermo)%T*(1-scale_v**2)/dt2**2

      do iat = 1, groups(group)%num_atom
         jat = groups(group)%list_atom(iat)

         do idir = 1,3
            random_a = sqrt(scale_random/ion_mass(jat))*services_maxboltzdist()
            a_t(idir,jat) = a_t(idir,jat) + random_a
         enddo
      enddo

      ! Make sure all procs have the same velocity and
      ! acceleration
      call comms_bcast(pub_root_proc_id,v_t(:,:))
      call comms_bcast(pub_root_proc_id,a_t(:,:))

      return

    end subroutine md_langevin_dynamics


!==============================================================================!
!==============================================================================!

    subroutine md_andersen_scaling(thermo,group,dt)

      !================================================================!
      ! Replace the velocities of some atoms with velocities derived   !
      ! from the Boltzmann distribution at the thermostat temperature  !
      !                                                                !
      !----------------------------------------------------------------!
      ! Implemented by Simon M.-M. Dubois (May 2011)                   !
      !================================================================!

      use comms,           only: pub_on_root, pub_root_proc_id, comms_bcast
      use constants,       only: dp
      use md_ions,         only: v_t, groups, ion_mass
      use services,        only: services_maxboltzdist

      implicit none

      ! Argument
      integer, intent(in) :: thermo
      integer, intent(in) :: group
      real(kind=DP), intent(in) :: dt

      ! Variables
      real(kind=DP) :: dt2, ran, sigma, pcol, mix1, mix2, vtmp
      integer  :: iat, jat, idir

      !====================================================================!

      dt2 = 0.5_dp*dt
      pcol = 1.0_dp - exp(-dt2/md_thermo(thermo)%tau)
      mix1 = sqrt(1-md_thermo(thermo)%mix**2)
      mix2 = md_thermo(thermo)%mix

      if (pub_on_root) then
         do iat = 1, groups(group)%num_atom
            jat = groups(group)%list_atom(iat)

            call random_number(ran)
            if (ran.lt.pcol) then
               sigma = sqrt(md_thermo(thermo)%T/ion_mass(jat))
               do idir = 1,3
                  vtmp = sigma*services_maxboltzdist()
                  v_t(idir,jat) = mix1*v_t(idir,jat)+mix2*vtmp
               enddo
            endif
         enddo
      endif

      ! Make sure all procs have the same velocity
      call comms_bcast(pub_root_proc_id,v_t(:,:))

      return

    end subroutine md_andersen_scaling


!==============================================================================!
!==============================================================================!

    subroutine md_integrate_nhc(thermo,group,dt)

      !================================================================!
      ! This subroutine propagates a system of particles coupled       !
      ! to a single Nose-Hoover chain of thermostat according to the   !
      ! NHC part of the extended Liouville evolution operator.         !
      !                                                                !
      ! This is part of the Nose-Hoover approach to NVT extended by    !
      ! Martyna, Tuckermann and Klein to a chain of thermostat.        !
      !                                                                !
      !----------------------------------------------------------------!
      ! Implemented by Simon M.-M. Dubois (May 2011)                   !
      !================================================================!

      use comms,           only: pub_on_root, pub_root_proc_id, comms_bcast
      use constants,       only: dp,stdout
      use md_ions,         only: v_t, groups, group_kinetic_energy

      implicit none

      ! Argument
      integer, intent(in) :: thermo
      integer, intent(in) :: group
      real(kind=DP), intent(in) :: dt

      ! Variables
      integer    :: nth
      integer    :: nstep
      integer    :: ndof
      real(kind=DP) :: wdt2, wdt4, wdt8
      real(kind=DP) :: akin
      real(kind=DP) :: scale, aa
      integer    :: is, inhc, iat, jat

      !====================================================================!


      ! Initialise NHC local parameters
      nth   = md_thermo(thermo)%nhc_length
      nstep = md_thermo(thermo)%nhc_integ_nstep
      ndof  = 3*groups(group)%num_atom

      scale = 1.0_dp
      wdt2  = dt/(2.0_dp*nstep)
      wdt4  = dt/(4.0_dp*nstep)
      wdt8  = dt/(8.0_dp*nstep)

      ! Get the group kinetic energy
      akin = 2.0_dp*group_kinetic_energy(group)

      ! Compute force on first thermostat in NHC
      md_thermo(thermo)%nhc_g(1) = (akin - ndof*md_thermo(thermo)%T)/md_thermo(thermo)%nhc_mass(1)

      ! Multiple time step procedure for the propagation
      do is = 1, md_thermo(thermo)%nhc_integ_nstep

         ! Update the thermostat velocities
         md_thermo(thermo)%nhc_v(nth) =  md_thermo(thermo)%nhc_v(nth) + wdt4*md_thermo(thermo)%nhc_g(nth)
         do inhc = 1, nth-1
            aa = exp(-wdt8*md_thermo(thermo)%nhc_v(nth+1-inhc))
            md_thermo(thermo)%nhc_v(nth-inhc) = md_thermo(thermo)%nhc_v(nth-inhc)*aa*aa &
                  + wdt4*md_thermo(thermo)%nhc_g(nth-inhc)*aa
         enddo

         ! Update the scaling factor for particle velocities
         aa = exp(-wdt2*md_thermo(thermo)%nhc_v(1))
         scale = scale*aa

         ! Update the forces on first thermostat in chain
         md_thermo(thermo)%nhc_g(1) = (scale*scale*akin - ndof * md_thermo(thermo)%T) &
                  / md_thermo(thermo)%nhc_mass(1)

         ! Update the thermostat positions
         do inhc = 1, nth
            md_thermo(thermo)%nhc_x(inhc) = md_thermo(thermo)%nhc_x(inhc) + wdt2 * md_thermo(thermo)%nhc_v(inhc)
         enddo

         ! Update the thermostat velocities and forces
         do inhc = 1, nth-1
            aa = exp(-wdt8*md_thermo(thermo)%nhc_v(inhc+1))
            md_thermo(thermo)%nhc_v(inhc) = md_thermo(thermo)%nhc_v(inhc) * aa * aa &
                  + wdt4 * md_thermo(thermo)%nhc_g(inhc) * aa
            md_thermo(thermo)%nhc_g(inhc+1) = (md_thermo(thermo)%nhc_mass(inhc) &
                  * md_thermo(thermo)%nhc_v(inhc) * md_thermo(thermo)%nhc_v(inhc) &
                  - md_thermo(thermo)%T) / md_thermo(thermo)%nhc_mass(inhc+1)
         enddo
         md_thermo(thermo)%nhc_v(nth) =  md_thermo(thermo)%nhc_v(nth) &
                  + wdt4*md_thermo(thermo)%nhc_g(nth)
      enddo

      ! Update the particle velocities according to the scaling factor
      do iat = 1, groups(group)%num_atom
         jat = groups(group)%list_atom(iat)

         v_t(:,jat) = scale  * v_t(:,jat)

      enddo

      ! Make sure all procs have the same velocity
      call comms_bcast(pub_root_proc_id,v_t(:,:))

      ! Stability check
      if (maxval(abs(md_thermo(thermo)%nhc_v)).gt. 1.0E8_dp) then
         if (pub_on_root) then
            write(stdout,*) 'Warning : dynamic of Nose-Hoover chain going unstable '
            write(stdout,*) 'Action  : increase the thermostat masses  '
            write(stdout,*) '   or/and decrease the step length for NHC integration '
            write(stdout,*) '   or/and decrease the setp lenght for MD integration. '
         endif
      endif

      return

    end subroutine md_integrate_nhc


!==============================================================================!
!==============================================================================!

    subroutine md_initialise_thermostats()

      !================================================================!
      !                                                                !
      !----------------------------------------------------------------!
      ! Implemented by Simon M.-M. Dubois (August 2011)                   !
      !================================================================!

      use constants,       only: dp
      use utils,           only: utils_alloc_check

      implicit none

      ! Local variables
      integer :: thermo, nhcl
      integer :: ierr

      !====================================================================!

      ! Loop over all thermostats
      thermo_loop : do thermo = 1, md_thermo_num

         md_thermo(thermo)%T = md_thermo(thermo)%Tinit

         if (md_thermo(thermo)%type == 3) then

            nhcl = md_thermo(thermo)%nhc_length

            ! Allocate dynamical variables
            allocate(md_thermo(thermo)%nhc_x(nhcl),stat=ierr)
            call utils_alloc_check('md_initialise_thermostats','md_thermo(thermo)%nhc_x',ierr)
            allocate(md_thermo(thermo)%nhc_v(nhcl),stat=ierr)
            call utils_alloc_check('md_initialise_thermostats','md_thermo(thermo)%nhc_v',ierr)
            allocate(md_thermo(thermo)%nhc_g(nhcl),stat=ierr)
            call utils_alloc_check('md_initialise_thermostats','md_thermo(thermo)%nhc_g',ierr)
            allocate(md_thermo(thermo)%nhc_mass(nhcl),stat=ierr)
            call utils_alloc_check('md_initialise_thermostats','md_thermo(thermo)%nhc_mass',ierr)

            ! Initialise dynamical variables
            md_thermo(thermo)%nhc_x = 0.0_dp
            md_thermo(thermo)%nhc_v = 0.0_dp
            md_thermo(thermo)%nhc_g = 0.0_dp
            md_thermo(thermo)%nhc_mass = 0.0_dp
         endif

      enddo thermo_loop

      return

    end subroutine md_initialise_thermostats

!==============================================================================!
!==============================================================================!

    subroutine md_update_thermostat(thermo,group,time)

      !================================================================!
      !                                                                !
      !----------------------------------------------------------------!
      ! Implemented by Simon M.-M. Dubois (May 2011)                   !
      !================================================================!

      use constants,       only: two_pi
      use md_ions,         only: groups
      use rundat,          only: md_delta_t
      use utils,           only: utils_alloc_check

      implicit none

      ! Argument
      integer, intent(in)        :: thermo
      integer, intent(in)        :: group
      real(kind=DP), intent(in)  :: time

      ! Local variables
      real(kind=DP)  :: time_start
      integer :: ndof

      !====================================================================!

      if (thermo .gt. 0) then

         time_start = (md_thermo(thermo)%start-1) * md_delta_t

         ! Update the thermostat temperature
         md_thermo(thermo)%T = md_thermo(thermo)%Tinit + &
            (time - time_start) * md_thermo(thermo)%Tgrad / md_delta_t


         ! If required, update the NHC masses parameters
         if (md_thermo(thermo)%type ==  3) then

            ndof = groups(group)%num_atom*3
            if (md_thermo(thermo)%nhc_upd) then
               md_thermo(thermo)%nhc_mass(:) = md_thermo(thermo)%T/(two_pi/md_thermo(thermo)%tau)**2
            else
               md_thermo(thermo)%nhc_mass(:) = md_thermo(thermo)%Tinit/(two_pi/md_thermo(thermo)%tau)**2
            endif
            md_thermo(thermo)%nhc_mass(1) = ndof*md_thermo(thermo)%nhc_mass(1)
         endif

      endif

      return

    end subroutine md_update_thermostat

!==============================================================================!
!==============================================================================!


    subroutine md_select_thermostat(thermo,group,time)

      !================================================================!
      !                                                                !
      !----------------------------------------------------------------!
      ! Implemented by Simon M.-M. Dubois (August 2011)                !
      !================================================================!

      use comms,           only: pub_on_root
      use md_ions,         only: groups
      use rundat,          only: md_delta_t
      use utils,           only: utils_abort

      implicit none

      ! Argument
      integer, intent(in)             :: group
      integer, intent(out)            :: thermo
      real(kind=DP), intent(in)       :: time

      ! Local variables
      real(kind=DP) :: time_start, time_stop
      integer       :: ith
      logical       :: found

      found = .false.

      ! Look for the thermostat corresponding to group
      do ith = 1, md_thermo_num
         time_start = (md_thermo(ith)%start-1) * md_delta_t + tiny(1.0_dp)
         time_stop = (md_thermo(ith)%stop + 0.5_dp) * md_delta_t
         if ((md_thermo(ith)%group .eq. groups(group)%label .or. &
             md_thermo(ith)%group .eq. 0 ) .and. time_start .lt. time .and. &
             time_stop .ge. time) then
            thermo = ith
            found = .true.
            exit
         endif
      enddo

      ! If not found, abort with an error message
      if (.not. found) then
         thermo = -1
         if (pub_on_root) then
            write(stdout,'(/a)') 'Error in md_select_thermostat: '
            write(stdout,'(a,i4.4,a)') 'Either no suitable thermostat found &
                 &for group ', groups(group)%label, ' (for this timestep)'
            write(stdout,'(/a)') 'or the thermostat block is missing altogether!'
            write(stdout,'(/a)') 'Please specify a thermostat block. &
                 &For an NVE calculation, please specify '
            write(stdout,'(/a)') 'in the input file a thermostat block &
                 &with "NONE" as name of the thermostat.'
            write(stdout,'(/a)') 'This is the only way in which you can &
                 &specify initial temperature for NVE.'
            call utils_abort('Error in md_select_thermostat: missing &
                 &thermostat information')
         endif
      endif

    end subroutine md_select_thermostat

!==============================================================================!
!==============================================================================!

    subroutine md_restore_thermostat(rstat)

      use constants,       only: stdout,dp
      use comms,           only: pub_on_root,pub_root_proc_id, comms_bcast
      use rundat,          only: pub_rootname, md_restart, md_global_restart
      use utils,           only: utils_unit, utils_close_unit_check, &
             utils_alloc_check, utils_dealloc_check, utils_open_unit_check

      implicit none

      ! Arguments
      integer, intent(out) :: rstat

      ! Local variables
      real(kind=dp), allocatable :: dbuf(:)
      integer           :: dbuf_length, nhcl
      integer           :: io_unit, io_stat
      integer           :: ith, il, idx
      integer           :: ierr
      character(len=80) :: filename


      if (pub_debug_on_root) write(stdout,'(/a)')  &
         'DEBUG: Entering md_restore_thermostat'

      ! Compute buffer length
      dbuf_length = 0
      do ith = 1, md_thermo_num
         dbuf_length = dbuf_length + 1
         if (md_thermo(ith)%type == 3 .and. md_thermo(ith)%nhc_length .ge. 1) then
            dbuf_length = dbuf_length + 4*md_thermo(ith)%nhc_length
         endif
      enddo

      if (dbuf_length == 0) then
         rstat = 1
         return
      endif

      ! Allocate buffer
      allocate(dbuf(dbuf_length),stat=ierr)
      call utils_alloc_check('md_restore_thermostat','dbuf',ierr)

      if (pub_on_root) then

         ! Find available unit specifier
         io_unit = utils_unit()
         if (md_restart .and. .not. md_global_restart) then
            write(filename,'(a,a)') trim(pub_rootname),'.thermo.restart'
         elseif (md_global_restart .and. .not. md_restart) then
            write(filename,'(a,a)') trim(pub_rootname),'.thermo.global.restart'
         end if

         ! Open .thermo file
         open(unit=io_unit,iostat=io_stat,file=filename,&
              access='SEQUENTIAL',form='UNFORMATTED',position='REWIND', &
              action='READ')
         call utils_open_unit_check('md_restore_thermostat',filename,io_stat)

         idx = 0
         thermo_loop1 : do il = 1, dbuf_length

            read(io_unit) dbuf(il)
            ! DEBUG
            !write(stdout,*) 'restore_thermostat ', dbuf(il)

         enddo thermo_loop1

         close(unit=io_unit,iostat=io_stat)
         call utils_close_unit_check('md_restore_thermostat',filename,io_stat)

      endif

      ! Read buffer into thermostat arrays
      call comms_bcast(pub_root_proc_id,dbuf,dbuf_length)

      idx = 0
      thermo_loop2 : do ith = 1, md_thermo_num

         idx = idx + 1
         md_thermo(ith)%T = dbuf(idx)

         nhcl = md_thermo(ith)%nhc_length
         if (md_thermo(ith)%type == 3 .and. nhcl .ge. 1) then
            allocate(md_thermo(ith)%nhc_x(nhcl),stat=ierr)
            call utils_alloc_check('md_restore_thermostats','md_thermo(ith)%nhc_x',ierr)
            allocate(md_thermo(ith)%nhc_v(nhcl),stat=ierr)
            call utils_alloc_check('md_restore_thermostats','md_thermo(ith)%nhc_v',ierr)
            allocate(md_thermo(ith)%nhc_g(nhcl),stat=ierr)
            call utils_alloc_check('md_restore_thermostats','md_thermo(ith)%nhc_g',ierr)
            allocate(md_thermo(ith)%nhc_mass(nhcl),stat=ierr)
            call utils_alloc_check('md_restore_thermostats','md_thermo(ith)%nhc_mass',ierr)

            do il = 1, nhcl
               md_thermo(ith)%nhc_x(il) = dbuf(idx+1)
               md_thermo(ith)%nhc_v(il) = dbuf(idx+2)
               md_thermo(ith)%nhc_g(il) = dbuf(idx+3)
               md_thermo(ith)%nhc_mass(il) = dbuf(idx+4)
               idx = idx + 4
            enddo
         endif

      enddo thermo_loop2

      ! Allocate buffer
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('md_restore_thermostat','dbuf',ierr)

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_restore_thermostat'

      return

    end subroutine md_restore_thermostat

!==============================================================================!
!==============================================================================!

    subroutine md_backup_thermostat(global,backup_thermo)

      use constants,       only: dp
      use comms,           only: pub_on_root
      use rundat,          only: pub_rootname
      use utils,           only: utils_unit, utils_open_unit_check, &
            utils_close_unit_check, utils_alloc_check, utils_dealloc_check

      implicit none
      logical, optional, intent(in) :: global
      logical, optional, intent(in) :: backup_thermo

      ! Local variables
      real(kind=dp), allocatable :: dbuf(:)
      integer           :: dbuf_length
      integer           :: io_unit, io_stat
      integer           :: ith, il
      integer           :: ierr
      logical           :: fileexists
      character(len=80) :: filename
      logical           :: loc_backup

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_backup_thermostat'

      loc_backup = .true.
      if (present(backup_thermo)) then
         if(.not. backup_thermo) then
            loc_backup = .false.
         end if
      end if

      ! Compute buffer length
      dbuf_length = 0
      do ith = 1, md_thermo_num
         dbuf_length = dbuf_length + 1
         if (md_thermo(ith)%type == 3 .and. md_thermo(ith)%nhc_length .ge. 1) then
            dbuf_length = dbuf_length + 4*md_thermo(ith)%nhc_length
         endif
      enddo

      if (dbuf_length == 0) return

      ! Allocate buffer
      allocate(dbuf(dbuf_length),stat=ierr)
      call utils_alloc_check('md_backup_thermostat','dbuf',ierr)

      ! Backup the old .thermo.restart file to prevent any loss
      if (pub_on_root .and. loc_backup) then

         io_unit = utils_unit()

         if(present(global)) then
           if(global) write(filename,'(a,a)') trim(pub_rootname),&
               &'.thermo.global.restart'
         else
            write(filename,'(a,a)') trim(pub_rootname),'.thermo.restart'
         end if

         inquire(file=filename,exist=fileexists)

         if (fileexists) then

            open(unit=io_unit,iostat=io_stat,file=filename,&
                 access='SEQUENTIAL',form='UNFORMATTED',position='REWIND', &
                 action='READ')

            if (io_stat == 0) then
               do il = 1, dbuf_length
                  read(io_unit) dbuf(il)
               enddo
               close(unit=io_unit,iostat=io_stat)
               call utils_close_unit_check('md_backup_thermostat',filename,io_stat)


               if(present(global)) then
                  if(global) write(filename,'(a,a)') trim(pub_rootname),&
                       &'.thermo.global.restart.old'
               else
                  write(filename,'(a,a)')trim(pub_rootname),&
                       &'.thermo.restart.old'
               end if
               open(unit=io_unit,file=filename,iostat=io_stat,status='REPLACE', &
                      form='UNFORMATTED',position='REWIND',action='WRITE')
               call utils_open_unit_check('md_backup_thermostat',filename,io_stat)

               do il = 1, dbuf_length
                  write(io_unit) dbuf(il)
                  ! DEBUG
                  !write(stdout,*) 'backup_thermostat ', dbuf(il)
               enddo
               close(unit=io_unit,iostat=io_stat)
               call utils_close_unit_check('md_backup_thermostat',filename,io_stat)

            endif
         endif
      endif

      ! Save the current variables in a new .thermo.restart file
      if (pub_on_root) then

         io_unit = utils_unit()

         if(present(global)) then
           if(global) write(filename,'(a,a)') trim(pub_rootname),&
               &'.thermo.global.restart'
         else
            write(filename,'(a,a)') trim(pub_rootname),'.thermo.restart'
         end if
         open(unit=io_unit,file=filename,iostat=io_stat,status='REPLACE', &
                form='UNFORMATTED',position='REWIND',action='WRITE')
         call utils_open_unit_check('md_backup_thermostat',filename,io_stat)

         do ith = 1, md_thermo_num
            write(io_unit) md_thermo(ith)%T
            ! DEBUG
            !write(stdout,*) 'save_thermostat ', md_thermo(ith)%T
            if (md_thermo(ith)%type == 3 .and. md_thermo(ith)%nhc_length .ge. 1) then
               do il = 1, md_thermo(ith)%nhc_length
                  write(io_unit) md_thermo(ith)%nhc_x(il)
                  write(io_unit) md_thermo(ith)%nhc_v(il)
                  write(io_unit) md_thermo(ith)%nhc_g(il)
                  write(io_unit) md_thermo(ith)%nhc_mass(il)
                  ! DEBUG
                  !write(stdout,*) 'save_thermostat ', md_thermo(ith)%nhc_x(il)
                  !write(stdout,*) 'save_thermostat ', md_thermo(ith)%nhc_v(il)
                  !write(stdout,*) 'save_thermostat ', md_thermo(ith)%nhc_g(il)
                  !write(stdout,*) 'save_thermostat ', md_thermo(ith)%nhc_mass(il)
               enddo
            endif
         enddo

         close(unit=io_unit,iostat=io_stat)
         call utils_close_unit_check('md_backup_thermostat',filename,io_stat)

      endif

      ! Dellocate buffer
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('md_backup_thermostat','dbuf',ierr)

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_backup_thermostat'

      return

    end subroutine md_backup_thermostat

!==============================================================================!
!==============================================================================!

    subroutine md_write_thermostat(time,global)

      use constants,       only: dp
      use comms,           only: pub_on_root
      use rundat,          only: pub_rootname
      use utils,           only: utils_unit, utils_open_unit_check, &
            utils_close_unit_check

      implicit none

      ! Argument
      real(kind=dp)     :: time
      logical, optional, intent(in) :: global

      ! Local variables
      character(len=80) :: filename
      character(len=12) :: tkind(6)
      character(len=80) :: position, stat
      character(len=20) :: tmp_length
      character(len=20) :: nhc_format
      integer           :: fileunit, iostatus
      integer           :: ith

      data tkind /'None','Andersen','Langevin','Nose-Hoover','Berendsen','Bussi'/

      if (time == 0.0_dp) then
         position = 'REWIND'
         stat = 'REPLACE'
      else
         position = 'APPEND'
         stat = 'UNKNOWN'
      endif

      if (pub_on_root) then
         fileunit = utils_unit()
         if(.not. present(global)) then
            ! Write to <filename>.thermo
            write(filename,'(a,a)') trim(pub_rootname),'.thermo'
         elseif (present(global)) then
            if(global) then
               ! Write to <filename>.md.global
               write(filename,'(a,a)') trim(pub_rootname),'.thermo.global'
            end if
         end if
         open(unit=fileunit,file=filename,iostat=iostatus,status=stat, &
                form='FORMATTED',position=position,action='WRITE')
         call utils_open_unit_check('md_write_thermostat',filename,iostatus)

         write(fileunit,1) time*aut2fs, 'fs'

         thermo_loop : do ith = 1, md_thermo_num
            write(fileunit,2) '       ', ith
            write(fileunit,3) '       ', tkind(md_thermo(ith)%type + 1)
            write(fileunit,4) 'T     :', md_thermo(ith)%T
            write(fileunit,6) 'Start :', md_thermo(ith)%start, 'Stop  :',md_thermo(ith)%stop
            write(fileunit,5) 'Tinit :', md_thermo(ith)%Tinit, 'Tgrad :',md_thermo(ith)%Tgrad
            write(fileunit,5) 'Tau   :', md_thermo(ith)%tau,  'Damp  :',md_thermo(ith)%damp
            write(fileunit,6) 'NHC   :', md_thermo(ith)%nhc_length, 'Int   :',md_thermo(ith)%nhc_integ_nstep

            if (md_thermo(ith)%type == 3 .and. md_thermo(ith)%nhc_length .ge. 1) then
               write(tmp_length,*) md_thermo(ith)%nhc_length
               write(nhc_format,*) '(9x,',trim(adjustl(tmp_length)),'(3x,es18.8e3))'
               write(fileunit,nhc_format) md_thermo(ith)%nhc_x(:)
               write(fileunit,nhc_format) md_thermo(ith)%nhc_v(:)
               write(fileunit,nhc_format) md_thermo(ith)%nhc_g(:)
               write(fileunit,nhc_format) md_thermo(ith)%nhc_mass(:)
            endif
         enddo thermo_loop
         write(fileunit,7)

         close(unit=fileunit,iostat=iostatus)
         call utils_close_unit_check('md_write_thermostat',filename,iostatus)

      endif

1     format(/,12x,es18.8e3,2x,a)
2     format(/,4x,a,1x,i8)
3     format(4x,a,4x,a12)
4     format(4x,a,1x,es18.8e3)
5     format(4x,2(a,1x,es18.8e3,3x))
6     format(4x,2(a,1x,i8,17x))
7     format(/,58('-'))

      return

    end subroutine md_write_thermostat

!==============================================================================!
!==============================================================================!

    subroutine md_initialise_ndof(ndof_r, time, nat)

      use md_ions,    only: group_num
      use rundat,     only: pub_open_localpseudo,pub_multigrid_hartree,&
                            pub_ii_energy_direct

      implicit none

      ! Arguments
      integer, intent(inout) :: ndof_r
      real(kind=DP), intent(in) :: time
      integer, intent(in)    :: nat

      ! Local variables
      integer :: ith,group
      logical :: found_thermo

      found_thermo = .false.

      do group = 1,group_num
         call md_select_thermostat(ith,group,time)
         ! vv: If any thermostat, i.e. NVT calculation
         if( (md_thermo(ith)%type /= 0) .and. (ith .gt. 0)) then
            found_thermo = .true.
         else
            found_thermo = .false.
         end if
      end do

      if(.not.found_thermo) then
         ! vv: Open boundary conditions ?
         if( pub_open_localpseudo .or. &
             pub_multigrid_hartree .or. pub_ii_energy_direct) then
            ! vv: Corner case: only two atoms
            if(nat == 2) then
               ndof_r = 5
            else
               ndof_r = 6
            end if
         else
            ! vv: Periodic boundary conditions ?
            ndof_r = 3
         end if
      ! vv: If thermostat
      else
         ! All the internal degrees of freedom are coupled with the external
         ! ones
         ndof_r = 0
      end if

      return

    end subroutine md_initialise_ndof

!==============================================================================!
!==============================================================================!

end module md_thermostat





