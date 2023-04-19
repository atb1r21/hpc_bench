! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-


!============================================================================!
!                                                                            !
!                        Molecular Dynamics Module                           !
!                                                                            !
!----------------------------------------------------------------------------!
! Written by by Simon M.-M. Dubois (May 2010) based on the original          !
! module created by Arash A. Mostofi (May 2006)                              !
!============================================================================!

module md_ions

  use constants, only: dp, stdout

  implicit none

  private

  ! Set of classical particles.
  type, public :: subset
     integer              :: label
     integer              :: num_atom
     integer, allocatable :: list_atom(:)
  end type subset

  ! MD global variables
  real(kind=dp), allocatable, public, save :: r_t(:,:)
  real(kind=dp), allocatable, public, save :: v_t(:,:)
  real(kind=dp), allocatable, public, save :: a_t(:,:)
  real(kind=dp), allocatable, public, save :: ion_mass(:)
  real(kind=dp), allocatable, public, save :: v0_t(:,:)
  real(kind=dp), public, save :: vtot(3)
  real(kind=dp), public, save :: Lcm(3)

  ! MD properties
  real(kind=dp), public, save :: md_volume
  real(kind=dp), public, save :: md_total_polarisation(3)
  real(kind=dp), public, save :: md_ionic_polarisation(3)
  real(kind=dp), public, save :: md_elec_polarisation(3)

  ! Subsets of atoms
  integer, public, save :: group_num
  type(subset), allocatable, public, save :: groups(:)

  ! Module parameters
  real(kind=dp), parameter :: electron_mass_si = 9.1093826e-31_dp
  real(kind=dp), parameter :: avogadro_si = 6.0221415e23_dp

  ! List of subroutines and functions :
  public :: kinetic_energy
  public :: group_kinetic_energy
  public :: md_allocate_struct
  public :: md_deallocate_struct
  public :: md_initialise_groups
  public :: md_initialise_ion_masses
  public :: md_initialise_positions
  public :: md_initialise_velocities
  public :: md_restore_trajectory
  public :: md_atomic_to_internal_vel
  public :: md_scale_velocities

  contains

!==============================================================================!
!==============================================================================!

    function kinetic_energy(nat)

      ! Modified to take nat as an argument to eliminate pub_par by
      ! Joseph Prentice, May 2018

      implicit none

      integer, intent(in) :: nat
      real(kind=dp) :: kinetic_energy
      integer :: iat, idir, ierr


      kinetic_energy = 0.0_dp

      do iat = 1, nat
         do idir = 1, 3
            kinetic_energy = kinetic_energy + ion_mass(iat)*v_t(idir,iat)**2
         enddo
      enddo
      kinetic_energy = 0.5_dp*kinetic_energy

      return

    end function kinetic_energy

!==============================================================================!
!==============================================================================!

    function group_kinetic_energy(group)

      implicit none

      integer, intent(in) :: group

      real(kind=dp) :: group_kinetic_energy
      integer       :: nat, iat, jat, idir

      nat = groups(group)%num_atom
      group_kinetic_energy = 0.0_dp
      do jat = 1, nat
         iat = groups(group)%list_atom(jat)
         do idir = 1, 3
            group_kinetic_energy = group_kinetic_energy + ion_mass(iat)*v_t(idir,iat)**2
         enddo
      enddo
      group_kinetic_energy = 0.5_dp*group_kinetic_energy

      return

    end function group_kinetic_energy

!==============================================================================!
!==============================================================================!

    subroutine md_allocate_struct(elements,nat)

      ! Modified to remove pub_par by Joseph Prentice, May 2018

      use ion,             only: element
      use utils,           only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Arguments
      integer, intent(in) :: nat
      type(ELEMENT), intent(in) :: elements(nat)

      ! Variables
      integer, allocatable :: group_nat(:)
      integer, allocatable :: group_label(:)
      integer, allocatable :: group_idx(:)
      integer :: ngroup
      integer :: iat, ig
      integer :: ierr

      ! Allocate global arrays
      allocate(r_t(3,nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct (md_main)','r_t',ierr)
      allocate(v_t(3,nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct (md_main)','v_t',ierr)
      allocate(a_t(3,nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct (md_main)','a_t',ierr)
      allocate(ion_mass(nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct (md_main)','ion_mass',ierr)
      allocate(v0_t(3,nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct (md_main)','v0_t',ierr)

      r_t=0.0_dp; v_t=0.0_dp; a_t=0.0_dp; ion_mass=0.0_dp; v0_t=0.0_dp

      ! Determine number of subset
      allocate(group_nat(nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct','group_nat',ierr)
      allocate(group_idx(nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct','group_idx',ierr)
      allocate(group_label(nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct','group_label',ierr)

      ! Determine number of groups
      ngroup = 1
      group_label(1) = elements(1)%group_id
      group_idx(1) = 1
      group_nat(1) = 1

      group_label(2:nat) = -1
      group_idx(2:nat) = -1
      group_nat(2:nat) = -1

      do iat = 2, nat
         do ig = 1, ngroup
            if (elements(iat)%group_id == group_label(ig)) then
               group_nat(ig) = group_nat(ig) + 1
               group_idx(iat) = ig
               exit
            endif
         enddo
         if (group_idx(iat) .lt. 0) then
            ngroup = ngroup + 1
            group_label(ngroup) = elements(iat)%group_id
            group_nat(ngroup) = 1
            group_idx(iat) = ngroup
         endif

      enddo
      group_num = ngroup

      ! Allocate atoms
      allocate(groups(group_num),stat=ierr)
      call utils_alloc_check('md_allocate_struct','groups',ierr)

      ! Allocate arrays for each group of atoms
      do ig = 1, group_num
         groups(ig)%label = group_label(ig)
         groups(ig)%num_atom = group_nat(ig)

         allocate(groups(ig)%list_atom(group_nat(ig)),stat=ierr)
         call utils_alloc_check('md_allocate_struct','groups(ig)%list_atom',ierr)
      enddo

      deallocate(group_nat,stat=ierr)
      call utils_dealloc_check('md_allocate_struct','group_nat',ierr)
      deallocate(group_idx,stat=ierr)
      call utils_dealloc_check('md_allocate_struct','group_idx',ierr)
      deallocate(group_label,stat=ierr)
      call utils_dealloc_check('md_allocate_struct','group_label',ierr)

      return

    end subroutine md_allocate_struct

!==============================================================================!
!==============================================================================!

    subroutine md_deallocate_struct()

      use utils,           only: utils_dealloc_check

      implicit none

      integer :: ierr, ig

      ! Deallocate global arrays
      deallocate(ion_mass,stat=ierr)
      call utils_dealloc_check('md_deallocate_struct','ion_mass',ierr)
      deallocate(a_t,stat=ierr)
      call utils_dealloc_check('md_deallocate_struct','a_t',ierr)
      deallocate(v_t,stat=ierr)
      call utils_dealloc_check('md_deallocate_struct','v_t',ierr)
      deallocate(r_t,stat=ierr)
      call utils_dealloc_check('md_deallocate_struct','r_t',ierr)
      deallocate(v0_t,stat=ierr)
      call utils_dealloc_check('md_deallocate_struct','v0_t',ierr)


      ! dellocate arrays for each group of atoms
      do ig = 1, group_num
         deallocate(groups(ig)%list_atom,stat=ierr)
         call utils_dealloc_check('md_deallocate_struct','groups(ig)%list_atom',ierr)
      enddo
      deallocate(groups,stat=ierr)
      call utils_dealloc_check('md_deallocate_struct','groups',ierr)

      return

    end subroutine md_deallocate_struct

!==============================================================================!
!==============================================================================!

    subroutine md_initialise_groups(elements,nat)

      ! Modified to remove pub_par by Joseph Prentice, May 2018

      use ion,              only: element

      implicit none

      ! Arguments
      integer, intent(in) :: nat
      type(ELEMENT), intent(in)  :: elements(nat)

      ! Local variables
      integer :: iat, jat, ig

      ! Create subset of atoms
      do ig = 1, group_num
         iat = 0
         do jat=1,nat
            if (elements(jat)%group_id == groups(ig)%label) then
               iat = iat + 1
               groups(ig)%list_atom(iat) = jat
            endif
         enddo
      enddo

    end subroutine md_initialise_groups

!==============================================================================!
!==============================================================================!

    subroutine md_initialise_ion_masses(elements,nat)

      ! Modified to remove pub_par by Joseph Prentice, May 2018

      use constants,        only: dp,periodic_table_mass
      use ion,              only: element

      implicit none

      ! Arguments
      integer, intent(in) :: nat
      type(ELEMENT), intent(in)  :: elements(nat)

      ! Local variables
      integer :: iat
      integer :: atom_Z

      ! Initialise ionic masses (convert masses to atomic units)
      do iat=1,nat
         atom_Z = elements(iat)%atomic_number
         ion_mass(iat) = periodic_table_mass(atom_Z)*1e-3_dp/avogadro_si/electron_mass_si
      enddo

    end subroutine md_initialise_ion_masses

!==============================================================================!
!==============================================================================!

    subroutine md_initialise_positions(elements,nat)

      ! Modified to remove pub_par by Joseph Prentice, May 2018

      use ion,              only: element

      implicit none

      ! Arguments
      integer, intent(in) :: nat
      type(ELEMENT), intent(in)  :: elements(nat)

      ! Local variables
      integer :: iat

      ! Initialise coordinates
      do iat=1,nat
         r_t(1,iat) = elements(iat)%centre%x
         r_t(2,iat) = elements(iat)%centre%y
         r_t(3,iat) = elements(iat)%centre%z
      enddo

    end subroutine md_initialise_positions

!==============================================================================!
!==============================================================================!

    !===============================================================!
    !                                                               !
    ! This subroutine transforms the atomic velocities into internal!
    ! velocities and viceversa.                                     !
    ! The actual transformation is dependent on the number of       !
    ! degrees of freedom.                                           !
    ! N.B. The internal velocities are used to calculate the        !
    ! correct temperature and kinetic energy.                       !
    ! However, the physical atomic velocities are printed out.      !
    !===============================================================!
    ! Written by Valerio Vitale April 2015                          !
    ! Modified to remove pub_par by Joseph Prentice, May 2018       !
    !===============================================================!

    subroutine md_atomic_to_internal_vel(group,ndor,dir,nat)

      use comms,          only: comms_bcast,pub_root_proc_id
      use constants,      only: dp
      use linalg,         only: linalg_invert_sym_matrix

      implicit none

      ! Arguments
      integer, intent(in)           :: nat
      integer, intent(in)           :: group
      integer, intent(in)           :: ndor
      character(len=*), intent(in)  :: dir

      ! Local variables
      real(kind=DP)     :: I_cm(3,3),invIL(3)
      real(kind=DP)     :: mtot,ptot(3),xtot(3)
      real(kind=DP)     :: tscale
      integer           :: ndof
      integer           :: iat,jat,icm,jcm,idir

      ! vv: Number of dofs
      ndof = 3*nat - ndor
      ! vv: Initialise variables
      mtot = 0.0_DP
      xtot(:) = 0.0_DP
      invIL(:) = 0.0_dp
      I_cm(:,:) = 0.0_dp


      ! vv: According to # of dofs perform the corresponding transformation
      if (ndor < 3) then
         if (index(dir,'backward') > 0 ) then
            v0_t(:,:) = v_t(:,:)
         elseif (index(dir,'forward') > 0 ) then
            v_t(:,:) = v0_t(:,:)
         end if

      elseif (ndor >= 3) then
         if (index(dir,'backward') > 0 ) then
            do iat = 1,groups(group)%num_atom
               jat = groups(group)%list_atom(iat)
               v0_t(:,jat) = v_t(:,jat) + vtot(:)
            end do
         else if (index(dir,'forward') > 0 ) then
            do iat = 1,groups(group)%num_atom
               jat = groups(group)%list_atom(iat)
               v_t(:,jat) = v0_t(:,jat) - vtot(:)
            end do
         end if
         if (ndor >= 4) then
            do iat = 1,groups(group)%num_atom
               jat = groups(group)%list_atom(iat)
               do idir = 1,3
                  xtot(idir) = xtot(idir) + ion_mass(jat)*r_t(idir,jat)
               end do

               mtot = mtot + ion_mass(jat)
            end do

            if (mtot .gt. tiny(1.0_dp)) then
               xtot(:) = xtot(:)/mtot
            else
               xtot = 0.0_dp
            end if

            ! vv: Inertia tensor with respect to the center of mass
            do icm = 1,3
               do jcm = 1,3
                  do iat = 1,groups(group)%num_atom
                     jat = groups(group)%list_atom(iat)
                     I_cm(icm,jcm) = I_cm(icm,jcm) + ion_mass(jat) * &
                          (r_t(icm,jat) - xtot(icm))*(r_t(jcm,jat) - xtot(jcm))
                  end do
               enddo
            end do

            ! vv: Invert the moment of inertia tensor
            call linalg_invert_sym_matrix(I_cm,3)

            do icm=1,3
               do jcm =1,3
                  invIL(icm) = invIL(icm) + I_cm(icm,jcm)*Lcm(jcm)
               end do
            end do
            if(index(dir,'backward') > 0 ) then
                ! vv: atomic velocities
                do iat = 1,groups(group)%num_atom
                   jat = groups(group)%list_atom(iat)
                   v0_t(1,jat) = v0_t(1,jat) + &
                        (invIL(2)*(r_t(3,jat)-xtot(3)) - &
                        invIL(3)*(r_t(2,jat)-xtot(2)))
                   v0_t(2,jat) = v0_t(2,jat) + &
                        (invIL(3)*(r_t(1,jat)-xtot(1)) - &
                        invIL(1)*(r_t(3,jat)-xtot(3)))
                   v0_t(3,jat) = v0_t(3,jat) + &
                        (invIL(1)*(r_t(2,jat)-xtot(2)) - &
                        invIL(2)*(r_t(1,jat)-xtot(1)))
                end do
             else if (index(dir,'forward') > 0 ) then
                do iat = 1,groups(group)%num_atom
                   jat = groups(group)%list_atom(iat)
                   v_t(1,jat) = v_t(1,jat) - &
                        (invIL(2)*(r_t(3,jat)-xtot(3)) - &
                        invIL(3)*(r_t(2,jat)-xtot(2)))
                   v_t(2,jat) = v_t(2,jat) - &
                        (invIL(3)*(r_t(1,jat)-xtot(1)) - &
                        invIL(1)*(r_t(3,jat)-xtot(3)))
                   v_t(3,jat) = v_t(3,jat) - &
                        (invIL(1)*(r_t(2,jat)-xtot(2)) - &
                        invIL(2)*(r_t(1,jat)-xtot(1)))
                end do
             end if
        end if
     end if

     ! Make sure all procs have the same velocity
     if (index(dir,'backward') > 0 ) then
        call comms_bcast(pub_root_proc_id,v0_t(:,:))
     else if (index(dir,'forward') > 0 ) then
        call comms_bcast(pub_root_proc_id,v_t(:,:))
     end if

     return

    end subroutine md_atomic_to_internal_vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine md_scale_velocities(group,ndof,temp)

     use comms,     only: pub_root_proc_id, comms_bcast
     use constants, only: DP, SAFE_DIV_EPS
     use utils,     only: utils_assert

     implicit none

     ! Arguments
     integer, intent(in) :: group
     integer, intent(in) :: ndof
     real(kind=DP), intent(in) :: temp

     ! Local variables
     real(kind=DP) :: mvtot, tscale
     integer       :: iat,jat,idir

     mvtot  = 0.0_DP
     tscale = 0.0_DP
     do iat = 1,groups(group)%num_atom
        jat = groups(group)%list_atom(iat)
        mvtot = mvtot + ion_mass(jat)*sum(v_t(:,jat)**2)
     end do

     if (ndof .gt. 0) then
        call utils_assert(abs(mvtot) > SAFE_DIV_EPS, 'md_scale_velocities: &
             &Total kinetic energy zero or practically zero.')
        tscale = sqrt(temp*ndof/mvtot)
     else
        tscale = 0.0_DP
     end if

     ! Rescale velocities
     if (temp == 0.0_DP) then
        ! no-op in this case
     else
        do iat = 1, groups(group)%num_atom
            jat = groups(group)%list_atom(iat)
            do idir = 1, 3
               v_t(idir,jat) = v_t(idir,jat)*tscale
            enddo
        end do
     end if

     ! Make sure all procs have the same velocity
     call comms_bcast(pub_root_proc_id,v0_t(:,:))

     return

    end subroutine md_scale_velocities

!==============================================================================!
!==============================================================================!

    subroutine md_initialise_velocities(group, temp, elements,ndof_r,nat)

      ! Modified to remove pub_par by Joseph Prentice, May 2018

      use comms,            only: comms_bcast, pub_root_proc_id, pub_on_root
      use constants,        only: dp, SAFE_DIV_EPS
      use ion,              only: element
      use linalg,           only: linalg_invert_sym_matrix
      use rundat,           only: md_init_velocities, pub_rand_seed
      use services,         only: services_maxboltzdist, &
           services_maxboltzdist_parallel
      use utils,            only: utils_assert, utils_alloc_check, &
           utils_dealloc_check

      implicit none

      ! Arguments
      integer, intent(in)         :: nat
      integer, intent(in)         :: group
      real(kind=dp), intent(in)   :: temp
      type(ELEMENT), intent(in)   :: elements(nat)
      integer, intent(in)         :: ndof_r

      ! Local variables
      real(kind=dp) :: I_cm(3,3),invIL(3)
      real(kind=dp) :: sigma, tscale
      real(kind=dp) :: mtot, mvtot, ptot(3), xtot(3)
      real(kind=dp), allocatable :: dist(:)
      integer :: iat, jat, idir, icm, jcm, ig, ib, istat
      integer :: ndof

      ! vv : Total degrees of freedom
      ndof = 3*nat - ndof_r
      mtot = 0.0_DP
      ptot(:) = 0.0_DP
      vtot(:) = 0.0_DP

      ! User defined atomic velocities
      if (.not. md_init_velocities) then
         do iat = 1, nat
            jat = groups(group)%list_atom(iat)
            v0_t(:,jat) = elements(iat)%ion_velocity(:)
         enddo

      ! Atomic velocities from a Maxwell-Boltzmann distribution at
      ! desired temperature
      else
      ! If we do not care about reproducibility then use the intrinsic FORTRAN
      ! function to generate random numbers in services_maxboltzdist
         if(pub_rand_seed == -1) then
            allocate(dist(1),stat=istat)
            call utils_alloc_check('md_initialise_velocity','dist',istat)
            do iat = 1, groups(group)%num_atom
               jat = groups(group)%list_atom(iat)
               sigma = sqrt(temp/ion_mass(jat))
               do idir = 1, 3
                  v0_t(idir,jat) = sigma*services_maxboltzdist()
                  ptot(idir) = ptot(idir) + ion_mass(jat)*v0_t(idir,jat)
               enddo
               mtot = mtot + ion_mass(jat)
            enddo
         else
       ! Otherwise generate the random numbers on the root proc with a linear
       ! congruental random number generator and distribute to all procs
            allocate(dist(3*nat),stat=istat)
            call utils_alloc_check('md_initialise_velocity','dist',istat)
            call services_maxboltzdist_parallel(dist,nat)
            ib = 0
            do iat = 1, groups(group)%num_atom
               jat = groups(group)%list_atom(iat)
               sigma = sqrt(temp/ion_mass(jat))
               do idir = 1, 3
                  ib = ib + 1
                  v0_t(idir,jat) = sigma*dist(ib)
                  ptot(idir) = ptot(idir) + ion_mass(jat)*v0_t(idir,jat)
               enddo
               mtot = mtot + ion_mass(jat)
            enddo
         end if

         call utils_assert(abs(mtot) > SAFE_DIV_EPS, 'md_initialise_velocities: &
              &Total mass is zero or practically zero.')
         vtot(:) = ptot(:)/mtot

         ptot(:) = 0.0_DP
         do iat = 1, groups(group)%num_atom
            jat = groups(group)%list_atom(iat)

            do idir = 1,3
               v0_t(idir,jat) = v0_t(idir,jat) - vtot(idir)
               ptot(idir) = ptot(idir) + ion_mass(jat)*v0_t(idir,jat)
            end do
         end do
         vtot(:) = ptot(:)/mtot

         deallocate(dist,stat=istat)
         call utils_dealloc_check('md_initialise_velocity','dist',istat)

      end if

      ! Make sure all procs have the same velocity
      call comms_bcast(pub_root_proc_id,v0_t(:,:))
      if (pub_on_root) write(stdout,*) 'done'

      ! vv : Transform atomic velocities to internal velocities, according to the
      !      degrees of freedom
      if (pub_on_root) write(stdout,'(a)',advance='no') &
          '- MD -: Transform atomic to internal velocities...'
      do ig = 1, group_num
         call md_atomic_to_internal_vel(ig,ndof_r,'forward',nat)
      end do
      if (pub_on_root) write(stdout,*) 'done'

      ! vv: Scale atomic velocities according to the external temperature
      if (md_init_velocities) then
         mvtot  = 0.d0
         tscale = 0.d0
         do iat = 1,groups(group)%num_atom
            jat = groups(group)%list_atom(iat)
            mvtot = mvtot + ion_mass(jat)*sum(v_t(:,jat)**2)
         end do

         call utils_assert(mvtot > SAFE_DIV_EPS, 'md_initialise_velocities: &
              &Total kinetic energy is zero or practially zero.')

         if (ndof .gt. 0) then
            tscale = sqrt(temp*ndof/mvtot)
         else
            tscale = 0.d0
         end if

         ! Rescale velocities
         do iat = 1, groups(group)%num_atom
             jat = groups(group)%list_atom(iat)
             do idir = 1, 3
                 v_t(idir,jat) = v_t(idir,jat)*tscale
             enddo
         end do
         ! Make sure all procs have the same velocity
         call comms_bcast(pub_root_proc_id,v_t(:,:))
      end if

      return

    end subroutine md_initialise_velocities


!==============================================================================!
!==============================================================================!

    subroutine md_restore_trajectory(time,rstat,nat)

      ! Modified to remove pub_par by Joseph Prentice, May 2018

      use constants,       only: stdout,dp
      use comms,           only: pub_on_root, pub_root_proc_id, comms_bcast
      use rundat,          only: md_iter_global, pub_rootname, &
           pub_debug_on_root, md_restart, md_global_restart
      use utils,           only: utils_alloc_check,  utils_unit, &
             utils_open_unit_check, utils_close_unit_check, &
             utils_dealloc_check, utils_abort

      implicit none

      ! Argument
      integer, intent(in)          :: nat
      real(kind=DP), intent(inout) :: time
      integer, intent(inout)       :: rstat

      ! Local variables
      real(kind=dp), allocatable :: dbuf(:)
      integer           :: idbuf ! vv: md_iter_global is an integer
      integer           :: dbuf_length
      integer           :: io_unit, io_stat
      integer           :: iat, idx, nr
      integer           :: ierr,read_status
      integer           :: tmp_int
      real(kind=DP)     :: tmp_dp
      character(len=80) :: filename

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_restore_trajectory'

      ! Compute buffer length
      if (md_restart .and. .not. md_global_restart) then
         dbuf_length = 6*nat + 7 ! vv: 6*N positions + velocities
                                     !     +3 for velocity of center of mass
                                     !     +3 for the angular momentum of com
                                     !     +1 md_time
                                     !     +1 for the global counter
      else if (md_global_restart .and. .not. md_restart) then
         dbuf_length = 9*nat + 7 ! accelerations too
      end if

      ! Allocate buffer
      allocate(dbuf(dbuf_length),stat=ierr)
      call utils_alloc_check('md_restore_trajectory','dbuf',ierr)

      if (pub_on_root) then

         ! Find available unit specifier
         io_unit = utils_unit()

         ! vv: Normal restart
         if (md_restart .and. .not. md_global_restart) then
            write(filename,'(a,a)') trim(pub_rootname),'.md.restart'
         ! vv: Global restart
         elseif (md_global_restart .and. .not. md_restart) then
            write(filename,'(a,a)') trim(pub_rootname),'.md.global.restart'
         end if

         open(unit=io_unit,iostat=io_stat,file=filename,&
              access='SEQUENTIAL',form='UNFORMATTED',position='REWIND', &
              action='READ')
         call utils_open_unit_check('md_restore_trajectory',filename,io_stat)

         ! vv : Make sure we have back compatibility with previous ONETEP
         !      versions
         read(io_unit) tmp_int
         nr=1
         do
            read(io_unit,iostat=read_status) tmp_dp
            if(read_status/=0) then
               exit
            end if
            nr=nr+1
         end do

         if (nr /= dbuf_length + 1) then
            write(stdout,'(a)') 'ERROR in md restore_trajectory: Since ONETEP &
                 &version 4.3.1.0'
            write(stdout,'(a)') 'an md.restart file generated with an earlier &
                 &version of ONETEP '
            write(stdout,'(a)') 'cannot be used. See documentation.'
            call utils_abort('ERROR in md_restore_trajectory: no &
                 &backward-compatibility with ONETEP versions &
                 &earlier than 4.3.1.0.')
         else
            rewind(io_unit)
         end if

        ! vv : Read md_iter_global in an integer buffer
        read(io_unit) idbuf
        ! vv : Read all the rest in a real buffer
        do idx = 1, dbuf_length
           read(io_unit) dbuf(idx)
           ! DEBUG
           !write(stdout,*) 'restore_trajectory ', dbuf(idx)
        enddo

        close(unit=io_unit,iostat=io_stat)
        call utils_close_unit_check('md_restore_trajectory',filename,io_stat)
        !else
        !   write(stdout,'(/,3a,/,a)') &
        !        '    WARNING : md_restore_trajectory failed to open "', &
        !     trim(filename),'"', &
        !     '    WARNING : Initialising positions and velocities from the &
        !     &input file!'

        !endif

      endif


      ! Let all the procs have the info
      call comms_bcast(pub_root_proc_id,idbuf)
      call comms_bcast(pub_root_proc_id,dbuf,dbuf_length)

      md_iter_global = idbuf
      time = dbuf(1)
      vtot(1:3) = dbuf(2:4)
      Lcm(1:3)  = dbuf(5:7)

      idx = 7
      do iat=1,nat
         r_t(1:3,iat) = dbuf(idx+1:idx+3)
         idx = idx +3
      enddo

      do iat=1,nat
         v0_t(1:3,iat) = dbuf(idx+1:idx+3)
         idx = idx +3
      enddo

      if (md_global_restart .and. .not. md_restart) then
         do iat=1,nat
            a_t(1:3,iat) = dbuf(idx+1:idx+3)
            idx = idx +3
         end do
      end if

      ! Deallocate buffer
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('md_restore_trajectory','dbuf',ierr)

      if (pub_debug_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_restore_trajectory'

      return

    end subroutine md_restore_trajectory

!==============================================================================!
!==============================================================================!


end module md_ions
