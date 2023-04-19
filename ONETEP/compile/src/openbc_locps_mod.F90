! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Jacek Dziedzic and Nicholas D.M. Hine
!
!   This module was created in July 2015 from bits of pseudopotentials_mod.F90
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module openbc_locps

  use constants, only : DP, PI, stdout
  use pseudopotentials, only: PSEUDO_SPECIES
  use paw, only: PAW_SPECIES
  use rundat, only: pub_paw

  public :: openbc_locps_init
  public :: openbc_locps_exit
  public :: openbc_locps_on_grid
  public :: openbc_locps_calc_forces

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! SUBROUTINES FOR PSEUDOPOTENTIAL INITIALISATION !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine openbc_locps_init(cell,pseudo_sp,paw_sp,num_pspecies)

    !========================================================================!
    ! Initialises open BC local pseudo.                                      !
    ! Can be safely called many times, does nothing if already initialised.  !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    !   None.                                                                !
    !------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 09/2013.                                  !
    !========================================================================!

    ! -----------------------------------------------------------------------

    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    integer, intent(in) :: num_pspecies
    type(CELL_INFO), intent(in) :: cell
    type(PSEUDO_SPECIES), intent(inout) :: pseudo_sp(num_pspecies)
    type(PAW_SPECIES), intent(inout) :: paw_sp(num_pspecies)

    if (pub_paw) then
       if(all(paw_sp(:)%openbc_locps_initialised)) return
    else
       if(all(pseudo_sp(:)%openbc_locps_initialised)) return
    end if

    call openbc_locps_radial(cell,pseudo_sp,paw_sp,num_pspecies)

  end subroutine openbc_locps_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine openbc_locps_exit(pseudo_sp,paw_sp,num_pspecies)
    !========================================================================!
    ! Deinitialises open BC local pseudo.                                    !
    ! Can be safely called many times, does nothing if already deinitialised !
    ! or never initialised.                                                  !
    !------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 09/2013.                                  !
    !========================================================================!

    use rundat, only: pub_forces_needed
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: num_pspecies
    type(PSEUDO_SPECIES), intent(inout) :: pseudo_sp(num_pspecies)
    type(PAW_SPECIES), intent(inout) :: paw_sp(num_pspecies)

    ! Local Variables
    integer :: ierr
    integer :: species

    ! -----------------------------------------------------------------------

    if (pub_paw) then
       if(.not.any(paw_sp(:)%openbc_locps_initialised)) return
    else
       if(.not.any(pseudo_sp(:)%openbc_locps_initialised)) return
    end if


    ! jd: Clean up
    do species = 1, num_pspecies
       if (pub_paw) then
          if (paw_sp(species)%openbc_locps_initialised) then
             deallocate(paw_sp(species)%vloc_lookup,stat=ierr)
             call utils_dealloc_check('openbc_locps_exit', &
                  'paw_sp(species)%vloc_lookup',ierr)

             if(pub_forces_needed) then
                deallocate(paw_sp(species)%vlocder_lookup,stat=ierr)
                call utils_dealloc_check('openbc_locps_exit', &
                     'paw_sp(species)%vlocder_lookup',ierr)
             end if
             paw_sp(species)%openbc_locps_initialised = .false.
          end if
       else
          if (pseudo_sp(species)%openbc_locps_initialised) then
             deallocate(pseudo_sp(species)%vloc_lookup,stat=ierr)
             call utils_dealloc_check('openbc_locps_exit', &
                  'pseudo_sp(species)%vloc_lookup',ierr)

             if(pub_forces_needed) then
                deallocate(pseudo_sp(species)%vlocder_lookup,stat=ierr)
                call utils_dealloc_check('openbc_locps_exit', &
                     'pseudo_sp(species)%vlocder_lookup',ierr)
             end if
             pseudo_sp(species)%openbc_locps_initialised = .false.
          end if
       end if

    end do

  end subroutine openbc_locps_exit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine openbc_locps_radial(cell,pseudo_sp,paw_sp,num_pspecies)
    !========================================================================!
    ! Calculates the local part of the pseudopotential on a radial grid,     !
    ! without resorting to discrete Fourier transforms, and thus realizing   !
    ! OPEN boundary conditions.                                              !
    ! This is achieved by calculating a continuous Fourier transform, to     !
    ! obtain Vloc(x) from Vloc(g), both on radial grids.                     !
    ! The integral is tricky, because:                                       !
    ! 1) for large values of x, the oscillations in sin(gx) become very      !
    !    frequent, requiring finer resolution than that of the recpot file,  !
    ! 2) the low-g part of the integral is singular,                         !
    ! 3) it must be computed to very, very high accuracy to get energies to  !
    !    an accuracy of 0.001%.                                              !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    !   None.                                                                !
    !------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2010, based on notes by Chris-Kriton  !
    ! Skylaris. Updated by Jacek Dziedzic in Sep 2010 to use the trick to    !
    ! split the potential into lr and sr parts in order to improve accuracy. !
    ! Separated into openbc_locps_radial (which only needs to be run once    !
    ! and openbc_locps_on_grid, which needs to be executed every time ions   !
    ! move, and added the calculation of the derivative, Jacek               !
    ! Dziedzic 09/2013. OMP-parallelized by Nick Hine somewhere in 2013.     !
    !========================================================================!

    use comms, only: pub_my_proc_id, pub_on_root, pub_total_num_procs, &
         comms_allgather
    use geometry, only: magnitude, operator(+), operator(-), operator(*)
    use rundat, only: pub_forces_needed, &
         pub_pspot_bc_is_periodic, &
         pub_openbc_pspot_finetune_f, pub_openbc_pspot_finetune_nptsx, &
         pub_openbc_pspot_finetune_alpha, pub_debug_on_root
!$  use rundat, only: pub_threads_max
    use services, only: services_radial_derivative
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_abort, utils_alloc_check, &
         utils_dealloc_check, utils_flush, utils_unit, utils_erf

    implicit none

    ! Arguments
    integer, intent(in) :: num_pspecies
    type(CELL_INFO), intent(in) :: cell
    type(PSEUDO_SPECIES), intent(inout) :: pseudo_sp(num_pspecies)
    type(PAW_SPECIES), intent(inout) :: paw_sp(num_pspecies)

    ! Local Variables
    real(kind=DP), allocatable :: vloc_lookup_local(:)
    ! jd: Portion of the above local to this core, for current species

    integer :: npts_x                    ! jd: Size of the radial realspace grid
    integer :: npts_x_per_core           ! jd: One core's share of the above
    integer :: my_npts_x                 ! jd: One core's share of the above
    integer :: npts_x_rem                ! ndmh: remainder of points/cores
    integer :: species, x_rad_pt         ! jd: Indices
    real(kind=DP) :: hx                  ! jd: Gridpoint separation
    real(kind=DP) :: xmax                ! jd: Last x stored in lookups
    real(kind=DP) :: x                   ! jd: Current x
    real(kind=DP) :: delta_g             ! jd: g spacing in the recpot file
    real(kind=DP) :: Z                   ! jd: Charge of current ion
    integer       :: x1, x2              ! jd: min and max x for this core
    real(kind=DP) :: coulombic           ! jd: Needed for self-test
    real(kind=DP) :: reldiff, absreldiff ! jd: Needed for self-test
    real(kind=DP) :: avgabsreldiff       ! jd: Needed for self-test
    real(kind=DP) :: avgreldiff          ! jd: Needed for self-test
    integer       :: avgpts              ! jd: Needed for self-test
    real(kind=DP) :: maxabsreldiff       ! jd: Needed for self-test
    integer       :: fineness            ! jd: fineness of the radial g-grid
    integer       :: ierr                ! jd: Error flag
    real(kind=DP) :: alpha               ! jd: Crossover between sr and lr term
    real(kind=DP) :: erf_alpha_x_over_x  ! jd: erf(alpha*x)/x
    character(len=128) :: string
    character(len=64) :: fn
    real(kind=DP), parameter :: factor = 2.0_DP / PI
    ! jd: 4pi (from ang. integration) * 4pi (onetep real/recip space convention)
    !     * 1/(2pi)**3 (Fourier transform factor) = 2/pi
    character(len=*), parameter :: myself = 'openbc_locps_radial'

    integer :: output_unit
    integer :: output_unit_der

    ! -----------------------------------------------------------------------

    call timer_clock(myself,1)

    if(pub_on_root) write(stdout,'(/a)') &
         'Calculating local pseudopotential in open BC on a radial grid ...'
    call utils_flush

    ! jd: Set up parameters
    npts_x = pub_openbc_pspot_finetune_nptsx
    alpha = pub_openbc_pspot_finetune_alpha / min( &
         magnitude(cell%a1),magnitude(cell%a2),magnitude(cell%a3))

    ! jd: Sanity check
    call utils_assert(.not.any(pub_pspot_bc_is_periodic),'Cannot use '//trim(myself)//&
         ' without having open boundary conditions in the local pseudo.')

    ! jd: Find out the maximum value of x which we might need, determine the
    !     fineness of the radial grid depending on this
    xmax = magnitude(cell%a1 + cell%a2 + cell%a3)
    hx = xmax/(npts_x-1)

    ! jd: Split the work across cores: every core will work with points x1..x2
    npts_x_per_core = npts_x / pub_total_num_procs
    npts_x_rem = mod(npts_x,pub_total_num_procs)
    my_npts_x = npts_x_per_core
    if (pub_my_proc_id<npts_x_rem) my_npts_x = my_npts_x + 1
    x1 = pub_my_proc_id*npts_x_per_core + min(pub_my_proc_id,npts_x_rem) + 1
    x2 = (pub_my_proc_id+1)*npts_x_per_core + min(pub_my_proc_id + 1,npts_x_rem)

    ! jd: Sanity check against a situation where the last core(s) is/are left
    !     with no work
    if(x1 > npts_x .or. x2 > npts_x) then
       call utils_abort('Too many cores to sensibly divide &
            &openbc_pspot_finetune_nptsx radial points in the open BC &
            &pseudopotential integration in '//trim(myself)//' Increase the &
            &number of points or reduce the number of cores.')
    end if

    ! jd: Allocate and zero lookup array.
    do species = 1, num_pspecies
       if (pub_paw) then
          allocate(paw_sp(species)%vloc_lookup(npts_x),stat=ierr)
          call utils_alloc_check(myself,'paw_sp(species)%vloc_lookup',ierr)
          paw_sp(species)%vloc_lookup = 0.0_DP

          ! jd: Same for the derivative, if forces needed
          if(pub_forces_needed) then
             allocate(paw_sp(species)%vlocder_lookup(npts_x), &
                  stat=ierr)
             call utils_alloc_check(myself,'paw_sp(species)%vlocder_lookup', ierr)
             paw_sp(species)%vlocder_lookup = 0.0_DP
          end if
       else
          allocate(pseudo_sp(species)%vloc_lookup(npts_x),stat=ierr)
          call utils_alloc_check(myself,'pseudo_sp(species)%vloc_lookup',ierr)
          pseudo_sp(species)%vloc_lookup = 0.0_DP

          ! jd: Same for the derivative, if forces needed
          if(pub_forces_needed) then
             allocate(pseudo_sp(species)%vlocder_lookup(npts_x), &
                  stat=ierr)
             call utils_alloc_check(myself,'pseudo_sp(species)%vlocder_lookup', ierr)
             pseudo_sp(species)%vlocder_lookup = 0.0_DP
          end if
       end if
    end do

    ! jd: Allocate the part of the lookup local to this core
    allocate(vloc_lookup_local(my_npts_x),stat=ierr)
    call utils_alloc_check(myself,'vloc_lookup_local', ierr)
    vloc_lookup_local = 0.0_DP

    if (pub_on_root) then
       write(stdout,'(a,f0.3,a)') 'Discrepancy between calculated realspace &
            &psloc and -Z/x on [5a0,', (npts_x-1)*hx, 'a0]:'
       write(stdout,'(a)') 'max |rel. diff| # avg |rel. diff| # &
            & avg rel. diff #  f # pseudopotential file'
    end if

    ! jd: Fill up the lookup arrays, each core deals with its own portion of x's
    do species = 1, num_pspecies

       if (pub_paw) then
          Z = paw_sp(species)%ion_charge
          delta_g = 1.0_DP / paw_sp(species)%inv_g_spacing
       else
          Z = pseudo_sp(species)%ion_charge
          delta_g = 1.0_DP / pseudo_sp(species)%inv_g_spacing
       end if

       ! jd: If the user did not choose an explicit value, select a sensible
       !     default for the fineness parameter. It is chosen such, that we
       !     still have at least 50 points to sample a full period of sin(gx)
       !     for the largest possible x, on the g-grid.
       if(pub_openbc_pspot_finetune_f == -1) then
          fineness = int(50.0_DP * delta_g / (2.0_DP*PI/xmax) + 1)
       else
          fineness = pub_openbc_pspot_finetune_f
       end if

       ! jd: Go over all of the x's that belong to this core
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(x_rad_pt,x,erf_alpha_x_over_x) &
!$OMP SHARED(vloc_lookup_local,x1,x2,species,Z,alpha,hx,fineness, &
!$OMP      pub_threads_max)
       do x_rad_pt = x1, x2
          x = (x_rad_pt-1) * hx

          if(x>1D-10) then
            erf_alpha_x_over_x = utils_erf(alpha*x)/x
          else
            ! jd: lim x->0 (erf(alpha*x)/x)
            erf_alpha_x_over_x = 2.0_DP*alpha/sqrt(PI)
          end if

          ! jd: Vloc(x) is separated into two parts: short-range and long-range,
          !     the short-range part is evaluated in internal_Is_of_x, the
          !     long-range part is -Z*erf(alpha*x)/x
          vloc_lookup_local(x_rad_pt-x1+1) = factor * &
               internal_Is_of_x(species,x,fineness,alpha) - &
               (Z*erf_alpha_x_over_x)

       end do
!$OMP END PARALLEL DO

       ! jd: Allgather the partial results vloc_lookup_local -> vloc_lookup
       if (pub_paw) then
          call comms_allgather(paw_sp(species)%vloc_lookup(:),vloc_lookup_local, &
               my_npts_x,(/-1/), (/-1/))
       else
          call comms_allgather(pseudo_sp(species)%vloc_lookup(:),vloc_lookup_local, &
               my_npts_x,(/-1/), (/-1/))
       end if

       ! jd: Calculate the derivative wrt x of Vloc(x), if forces are needed
       !     Calculating this in parallel is not straigtforward, as there are
       !     domain-boundary effects in services_radial_derivative (diff. order
       !     changes close to endpoints). We don't want the derivative to depend
       !     on ncores. Could probably be done on root and then broadcast, but
       !     this happens only once per calculation and is fast...
       if(pub_forces_needed) then
          if (pub_paw) then
             call services_radial_derivative(paw_sp(species)%vlocder_lookup(:), &
                  paw_sp(species)%vloc_lookup(:), npts_x, xmax)
          else
             call services_radial_derivative(pseudo_sp(species)%vlocder_lookup(:), &
                  pseudo_sp(species)%vloc_lookup(:), npts_x, xmax)
          end if
       end if

       ! jd: Sanity-check on the result and gather some statistics on it
       maxabsreldiff = 0.0_DP
       avgabsreldiff = 0.0_DP
       avgreldiff = 0.0_DP
       avgpts = 0
       if (pub_paw) then
          fn = paw_sp(species)%ds_name
       else
          fn = pseudo_sp(species)%pseudo_name
       end if

       do x_rad_pt = 1, npts_x

          x = (x_rad_pt-1) * hx

          ! jd: Output the calculated pseudopotential, and derivative,
          !     if calculated, to a file
          if(pub_debug_on_root) then
             if(x_rad_pt == 1) then
                output_unit = utils_unit()
                open(unit = output_unit, file = trim(fn)//'.calculated_realpot', &
                     action = "write", err=10)
                if(pub_forces_needed) then
                   output_unit_der = utils_unit()
                   open(unit = output_unit_der, file = trim(fn) &
                        //'.calculated_realpot_derivative', &
                        action = "write", err=30)
                end if
             end if
             if (pub_paw) then
                write(output_unit,'(f20.15,f20.15)') &
                     x, paw_sp(species)%vloc_lookup(x_rad_pt)
                if(pub_forces_needed) then
                   write(output_unit_der,'(f20.15,f20.15)') &
                        x, paw_sp(species)%vlocder_lookup(x_rad_pt)
                end if
             else
                write(output_unit,'(f20.15,f20.15)') &
                     x, pseudo_sp(species)%vloc_lookup(x_rad_pt)
                if(pub_forces_needed) then
                   write(output_unit_der,'(f20.15,f20.15)') &
                        x, pseudo_sp(species)%vlocder_lookup(x_rad_pt)
                end if
             end if
             if(x_rad_pt == npts_x) then
                close(output_unit,err=20)
                if(pub_forces_needed) close (output_unit_der,err=40)
             end if
          end if

          ! jd: See if the tail of the potential closely matches -Z/x, it should
          if(x > 5.0_DP) then
             coulombic = -Z/x
             if (pub_paw) then
                reldiff=(paw_sp(species)%vloc_lookup(x_rad_pt) - coulombic)/coulombic
             else
                reldiff=(pseudo_sp(species)%vloc_lookup(x_rad_pt) - coulombic)/coulombic
             end if
             absreldiff=abs(reldiff)

             if(absreldiff > maxabsreldiff) maxabsreldiff = absreldiff
             avgreldiff = avgreldiff + reldiff
             avgabsreldiff = avgabsreldiff + absreldiff
             avgpts = avgpts + 1

             if(absreldiff>0.03) then
                if(pub_on_root) then
                   if (pub_paw) then
                      write(string,'(a,f12.5,a,e12.5,a,e12.5)') 'x= ',x, &
                           ' V_Coul=', coulombic,' V_obtained=', &
                           paw_sp(species)%vloc_lookup(x_rad_pt)
                   else
                      write(string,'(a,f12.5,a,e12.5,a,e12.5)') 'x= ',x, &
                           ' V_Coul=', coulombic,' V_obtained=', &
                           pseudo_sp(species)%vloc_lookup(x_rad_pt)
                   end if
                end if
                call utils_abort('Discrepancy of more than 3% between the &
                     &calculated pseudopotential and Coulombic potential, in &
                     &the tail of the potential, detected in '//trim(myself)//&
                     '. Try increasing openbc_pspot_finetune_f or the KE &
                     &cutoff. '//trim(string))
             end if
          end if

       end do ! radial points

       ! jd: Make sure we do not divide by zero
       call utils_assert(avgpts > 0, 'Internal error in '//trim(myself)//&
            '. Check that your simulation cell is sane.')

       ! jd: Calculate and print out maximum and average relative difference
       !     from Coulombic in potential tail
       avgreldiff = avgreldiff / avgpts
       avgabsreldiff = avgabsreldiff / avgpts

       if (pub_on_root) then
          write(stdout,'(a,f14.12,a,f15.12,a,f15.12,a,i3,a,a,a)') &
               ' ',maxabsreldiff, '   ',avgabsreldiff,'  ',avgreldiff, &
               '  ',fineness,'  ''',trim(fn), ''''
       end if

       if (pub_paw) then
          paw_sp(species)%openbc_locps_initialised = .true.
       else
          pseudo_sp(species)%openbc_locps_initialised = .true.
       end if

    end do ! species

    deallocate(vloc_lookup_local,stat=ierr)
    call utils_dealloc_check(myself,'vloc_lookup_local', ierr)

    if (pub_on_root) then
       write(stdout,'(a)') '... done'
    end if

    call utils_flush

    call timer_clock(myself,2)

    return

10  call utils_abort('Could not open a debug file to write the&
         & real space pseudopotential to.')
20  call utils_abort('Could not close a debug file with the&
         & real space pseudopotential.')
30  call utils_abort('Could not open a debug file to write the&
         & real space pseudopotential derivative to.')
40  call utils_abort('Could not close a debug file with the&
         & real space pseudopotential derivative.')

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real(kind=DP) function internal_Is_of_x(s,x,fineness,alpha)
      !========================================================================!
      ! Calculates the integral of                                             !
      ! ( Vloc(g) * sin(gx)/x * g * (1-exp(-g^2/(4*alpha^2))dg from 0 to g_cut,!
      ! where g_cut is the highest g supported by the current reciprocal space !
      ! grid, which will be somewhat smaller than g_max, up to which the recpot!
      ! is defined in the recpot file.                                         !
      ! The integral is tricky in that for large x it starts to oscillate very !
      ! heavily and the g resolution in the recpot file is too low to just     !
      ! integrate over the points in rad_locpot_recip, even with the Newton-   !
      ! Cotes integrator. Thus, integration is performed on a finer grid, cf.  !
      ! the parameter 'fineness'. Unless the user chose to override it, this   !
      ! will have been set to a sensible default by the parent routine.        !
      ! A check in the parent routine checks the high-x value for adequate     !
      ! closeness to the Coulombic potential.                                  !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      ! s (input)       : number of the species whose Vloc we want.            !
      ! x (input)       : the parameter x of the integration.                  !
      ! fineness (input): the fine radial g-grid will be this much finer.      !
      ! alpha (input)   : a parameter controlling the crossover between the    !
      !                   short-range and long-range parts.                    !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in May 2010, based on notes by Chris-Kriton  !
      ! Skylaris. Updated by Jacek Dziedzic in Sep 2010 to use the trick to    !
      ! split the potential into lr and sr parts in order to improve accuracy. !
      !========================================================================!

      use services, only: services_1d_interpolation
      use rundat, only: pub_fine_grid_scale
      use utils, only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! jd: Arguments
      integer, intent(in) :: s           ! jd: Species in question
      real(kind=DP), intent(in) :: x     ! jd: Parameter in the integration
      integer, intent(in) :: fineness    ! jd: Fineness of the radial g-grid
      real(kind=DP), intent(in) :: alpha ! jd: Parameter

      ! jd: Local variables
      real(kind=DP), allocatable :: integrand(:) ! jd: Integrated quantity
      real(kind=DP) :: fac               ! jd: Newton-Cotes integration factor
      real(kind=DP) :: dg                ! jd: g spacing in the recpot file
      real(kind=DP) :: dg_fine           ! jd: g spacing in the fine radial grid
      real(kind=DP) :: g                 ! jd: Current g
      real(kind=DP) :: g_cut             ! jd: Maximum g in the recip grid
      real(kind=DP) :: sin_gx_over_x     ! jd: sin(gx)/x or limit thereof
      real(kind=DP) :: integral          ! jd: Accumulator for the integral
      real(kind=DP) :: v_loc_value       ! jd: Vloc(g) for current g
      integer :: ngpts, ngpts_fine       ! jd: Number of g points for both grids
      integer :: g_rad_pt_fine           ! jd: Index
      integer :: k                       ! jd: Counter
      integer :: ierr                    ! jd: Error flag

      ! -----------------------------------------------------------------------

      ! jd: Basic variables
      if (pub_paw) then
         Z = paw_sp(s)%ion_charge
         ngpts = paw_sp(s)%n_recip_pts
         dg = 1.0_DP / paw_sp(s)%inv_g_spacing
      else
         Z = pseudo_sp(s)%ion_charge
         ngpts = pseudo_sp(s)%n_rad_pts
         dg = 1.0_DP / pseudo_sp(s)%inv_g_spacing
      end if

      ! jd: Set up the integration domain
      dg_fine = dg/fineness

      ! jd: Ignore all g's beyond the maximum that is representable in our
      !     reciprocal space grid
      g_cut = pub_fine_grid_scale*PI/max(cell%d1,cell%d2,cell%d3)

      ! jd: Integration weight for Newton-Cotes (25.4.17, Abramowitz, Stegun)
      fac = (7.0_DP/17280.0_DP) * dg_fine

      ! jd: Allocate array holding the integrand
      allocate(integrand(ngpts*fineness),stat=ierr)
      call utils_alloc_check('internal_Is_of_x','integrand',ierr)
      integrand(:) = 0.0_DP

      ! jd: Prepare integrand by populating the array
      k=1
      do g_rad_pt_fine = 1,ngpts*fineness-1

         g = (g_rad_pt_fine-1) * dg_fine

         if(g>g_cut) exit

         if (pub_paw) then
            v_loc_value = services_1d_interpolation( &
                 paw_sp(s)%vhntzc_recip,ngpts, &
                 real(g_rad_pt_fine-1,kind=DP)/real(fineness,kind=DP),0)! jd:-1,sic
         else
            v_loc_value = services_1d_interpolation( &
                 pseudo_sp(s)%rad_locpot_recip,ngpts, &
                 real(g_rad_pt_fine-1,kind=DP)/real(fineness,kind=DP),0)! jd:-1,sic
         end if

         if(g < 1D-10) then
            v_loc_value = 0.0_DP  ! integrand(g=0) is zero
         else

            ! jd: Add back Coulombic part, subtracted earlier
            if(alpha /= 0.0_DP) then
               v_loc_value = v_loc_value -Z/(g*g)
            end if

         end if

         if(x < 1D-10) then ! jd: Avoid singularity at 0
            sin_gx_over_x = g
         else
            sin_gx_over_x = sin(g*x)/x
         end if

         if(alpha /= 0.0_DP) then
            integrand(k) = v_loc_value * (1-exp(-g*g/(4.0_DP*alpha*alpha))) * &
                 sin_gx_over_x * g
         else
            integrand(k) = v_loc_value * sin_gx_over_x * g
         end if

         k = k+1
      end do
      ngpts_fine = k-1

      ! jd: Integrate with 8-point Newton-Cotes
      integral = 0.0_DP
      do g_rad_pt_fine=1, ngpts_fine-7,7
         integral = integral + &
              751.0_DP * &
              (integrand(g_rad_pt_fine)+integrand(g_rad_pt_fine+7)) + &
              3577.0_DP* &
              (integrand(g_rad_pt_fine+1)+integrand(g_rad_pt_fine+6)) + &
              1323.0_DP* &
              (integrand(g_rad_pt_fine+2)+integrand(g_rad_pt_fine+5)) + &
              2989.0_DP* &
              (integrand(g_rad_pt_fine+3)+integrand(g_rad_pt_fine+4))
      end do

      internal_Is_of_x = integral * fac

      deallocate(integrand,stat=ierr)
      call utils_dealloc_check('internal_Is_of_x','integrand',ierr)

    end function internal_Is_of_x

  end subroutine openbc_locps_radial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine openbc_locps_on_grid(v_local_fine, grid, cell, elements, &
       pseudo_sp, paw_sp, par)
    !========================================================================!
    ! Calculates the local part of the pseudopotential on the fine grid,     !
    ! without resorting to discrete Fourier transforms, and thus realizing   !
    ! OPEN boundary conditions. Uses a radial lookup prepared in advance by  !
    ! openbc_locps_radial.                                            !
    !------------------------------------------------------------------------!
    ! Arguments                                                              !
    ! v_local_fine (output): Contains v_loc on the fine grid on output.      !
    ! grid (input): The grid on which the local pseudo is to be generated.   !
    ! elements (input): The elements array req'd to perform the calculation. !
    !------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2010, based on notes by Chris-Kriton  !
    ! Skylaris. Updated by Jacek Dziedzic in Sep 2010 to use the trick to    !
    ! split the potential into lr and sr parts in order to improve accuracy. !
    ! Separated into openbc_locps_radial (which only needs to be run         !
    ! once) and openbc_locps_on_grid, which needs to be executed every       !
    ! time ions move, and added the calculation of the derivative, Jacek     !
    ! Dziedzic 09/2013. OMP-parallelized by Nick Hine somewhere in 2013.     !
    ! Modified to remove pub_par by Joseph Prentice, May 2019                !
    !========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_my_proc_id, pub_on_root
    use geometry, only: magnitude, POINT, operator(+), operator(-), operator(*)
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_openbc_pspot_finetune_nptsx
!$  use rundat, only: pub_threads_max
    use services, only: services_1d_interpolation
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_flush

    implicit none

    ! jd: Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    real(kind=DP), intent(out)  :: v_local_fine(grid%ld1, &
         grid%ld2, grid%max_slabs12)
    type(PARAL_INFO), intent(in) :: par
    type(ELEMENT), intent(in)   :: elements(par%nat)
    type(PSEUDO_SPECIES), intent(in) :: pseudo_sp(par%num_pspecies)
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)

    integer :: npts_x                    ! jd: Size of the radial realspace grid
    integer :: species,           I      ! jd: Indices
    integer :: islab12, ipt, i1, i2, i3  ! jd: Indices
    real(kind=DP) :: accum               ! jd: Accumulator
    real(kind=DP) :: x                   ! jd: Current x
    real(kind=DP) :: xmax                ! jd: Last x stored in lookups
    type(POINT)   :: R_I                 ! jd: Position of ion I
    type(POINT)   :: r                   ! jd: Position of curr. fine gridpoint
    real(kind=DP) :: Vloc_here           ! jd: Interpolated vloc at curr. gridpt
    character(len=*), parameter :: myself = 'openbc_locps_on_grid'

    ! -----------------------------------------------------------------------

    call timer_clock(myself,1)

    if(pub_on_root) write(stdout,'(a)') 'Calculating local pseudopotential in &
         &open BC on the fine real space grid ...'
    call utils_flush

    if (pub_paw) then
       call utils_assert(all(paw_sp(:)%openbc_locps_initialised), &
            "Internal error in "//trim(myself)//". The local openbc has not been initialised.")
    else
       call utils_assert(all(pseudo_sp(:)%openbc_locps_initialised), &
            "Internal error in "//trim(myself)//". The local openbc has not been initialised.")
    end if

    ! jd: Set up parameters
    npts_x = pub_openbc_pspot_finetune_nptsx
    xmax = magnitude(cell%a1 + cell%a2 + cell%a3)

    ! jd: Fill, in a distributed fashion, the v_local_fine using the lookups
    ! jd: Loop over all my points on the fine grid
    v_local_fine = 0.0_DP ! jd: Takes care of padding between pt and ld

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,i3,r,accum,I,R_I,species,x,vloc_here) &
!$OMP SHARED(par,elements,grid,pseudo_sp,paw_sp,pub_paw, &
!$OMP      v_local_fine,xmax,npts_x,pub_my_proc_id,pub_threads_max)
    do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
       i1 = modulo(ipt-1,grid%n1) + 1
       i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
       islab12 = (ipt-(i2-1)*grid%n1-i1) / &
            (grid%n1*grid%n2) + 1
       i3 = grid%first_slab12(pub_my_proc_id) + islab12 - 1

       r = &
            real((i1-1),kind=DP) * grid%da1 + &
            real((i2-1),kind=DP) * grid%da2 + &
            real((i3-1),kind=DP) * grid%da3

       ! jd: Loop over all the ions
       accum = 0.0_DP
       do I = 1, par%nat

          R_I = elements(I)%centre
          x = magnitude(r-R_I)
          species = elements(I)%pspecies_number

          ! jd: Interpolate Vloc(x) from radial lookups
          if (pub_paw) then
             vloc_here = services_1d_interpolation( &
                  paw_sp(species)%vloc_lookup(:), &
                  npts_x, x/xmax*real(npts_x-1,kind=DP),0) ! -1, sic!
          else
             vloc_here = services_1d_interpolation( &
                  pseudo_sp(species)%vloc_lookup(:), &
                  npts_x, x/xmax*real(npts_x-1,kind=DP),0) ! -1, sic!
          end if

          ! jd: It's faster to accumulate in a variable first
          accum = accum + vloc_here

       end do ! over I

       ! jd: Add the contribution from all ions to the current point
       v_local_fine(i1,i2,islab12) = accum

    end do ! ipt
!$OMP END PARALLEL DO

    if (pub_on_root) then
       write(stdout,'(a/)') '... done'
    end if

    call utils_flush

    call timer_clock(myself,2)

  end subroutine openbc_locps_on_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine openbc_locps_calc_forces(den_slabs12, grid, cell, elements, &
       pseudo_sp, paw_sp, locps_forces, par)

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the local part of the ionic pseudopotential under open BCs.        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   Identical to pseudo_local_calculate_forces.                           !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 09/2013 using code snippets from           !
    ! pseudo_local_calculate_forces.                                          !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_my_proc_id, comms_reduce
    use constants, only: SQRT_PI
    use geometry, only: magnitude, POINT, operator(-), operator(*)
! jme: workaround for cray compiler bug 843178
#if _CRAYFTN
    use geometry, only: add_points
#else
    use geometry, only: operator(+)
#endif
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_openbc_pspot_finetune_nptsx, &
         pub_forces_needed, pub_is_smeared_ion_rep, &
         pub_is_smeared_ion_width, pub_num_spins, pub_debug_on_root
!$  use rundat, only: pub_threads_max
    use services, only: services_1d_interpolation
    use simulation_cell, only : CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_erf

    implicit none

    ! Arguments
    type(PARAL_INFO), intent(in) :: par
    type(CELL_INFO), intent(in) :: cell
    type(GRID_INFO), intent(in)  :: grid
    real(kind=DP), intent(inout) :: den_slabs12(grid%ld1, &
         grid%ld2, grid%max_slabs12, pub_num_spins)
    type(ELEMENT), intent(in)    :: elements(par%nat)
    type(PSEUDO_SPECIES), intent(in) :: pseudo_sp(par%num_pspecies)
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)
    real(kind=DP), intent(out)   :: locps_forces(1:3,par%nat)

    ! Local Variables
    integer       :: ipt, i1, i2, i3, islab12
    integer       :: npts_x
    integer       :: I
    integer       :: species
    real(kind=DP) :: x, x2
    real(kind=DP) :: xmax
    real(kind=DP) :: vlocder_here
    real(kind=DP) :: scalar_term
    type(POINT)   :: r
    type(POINT)   :: R_I
    type(POINT)   :: r_minus_R_I
    real(kind=DP) :: sigma_I
    real(kind=DP) :: Z_I
    real(kind=DP) :: term
    real(kind=DP), parameter :: two_over_sqrt_pi = 2.0_DP/SQRT_PI

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering openbc_locps_calc_forces'

    ! Start timer
    call timer_clock('openbc_locps_calc_forces',1)

    call utils_assert(pub_forces_needed, &
         'Internal error in openbc_locps_calc_forces.')

    ! jd: Set up parameters
    npts_x = pub_openbc_pspot_finetune_nptsx
! jme: workaround for cray compiler bug 843178
#if _CRAYFTN
    xmax = magnitude(add_points(add_points(cell%a1, cell%a2), cell%a3))
#else
    xmax = magnitude(cell%a1 + cell%a2 + cell%a3)
#endif

    ! If spin polarised, put total density in up spin
    if (pub_num_spins == 2) den_slabs12(:,:,:,1) = den_slabs12(:,:,:,1) + &
         den_slabs12(:,:,:,2)

    locps_forces(:,:) = 0.0_DP

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,i3,r,I,R_I,r_minus_R_I,species,x,x2,vlocder_here, &
!$OMP      scalar_term,term,Z_I,sigma_I) &
!$OMP SHARED(par,elements,grid,pseudo_sp,paw_sp,pub_paw,pub_is_smeared_ion_rep, &
!$OMP      xmax,npts_x,pub_my_proc_id,pub_threads_max,den_slabs12, &
!$OMP      pub_is_smeared_ion_width) &
!$OMP REDUCTION(+:locps_forces)
!$OMP DO

    ! jd: Loop over local grid points
    do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
       i1 = modulo(ipt-1,grid%n1) + 1
       i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
       islab12 = (ipt-(i2-1)*grid%n1-i1) / &
            (grid%n1*grid%n2) + 1
       i3 = grid%first_slab12(pub_my_proc_id) + islab12 - 1

! jme: workaround for cray compiler bug 843178
#if _CRAYFTN
       r = add_points(add_points((i1-1)*grid%da1, (i2-1)*grid%da2), &
                      (i3-1)*grid%da3)
#else
       r = (i1-1) * grid%da1 + (i2-1) * grid%da2 + (i3-1) * grid%da3
#endif

       ! jd: Loop over all the ions
       do I = 1, par%nat

          R_I = elements(I)%centre
          species = elements(I)%pspecies_number

          r_minus_R_I = r - R_I

          x = magnitude(r_minus_R_I)

          if(x >= 1E-10) then ! Usual case, atom off a gridpoint

             ! jd: Interpolate Vloc(x) from radial lookups
             if (pub_paw) then
                vlocder_here = services_1d_interpolation( &
                     paw_sp(species)%vlocder_lookup(:), &
                     npts_x, x/xmax*real(npts_x-1,kind=DP),0) ! -1, sic!
             else
                vlocder_here = services_1d_interpolation( &
                     pseudo_sp(species)%vlocder_lookup(:), &
                     npts_x, x/xmax*real(npts_x-1,kind=DP),0) ! -1, sic!
             end if

             term = vlocder_here

             ! jd: Apply smeared-ion correction, if needed. This formally
             !     belongs to is_smeared_ions_mod, but a function call here
             !     would be costly.
             if(pub_is_smeared_ion_rep) then
                x2 = x*x
                Z_I = real(elements(I)%ion_charge,kind=DP)
                sigma_I = pub_is_smeared_ion_width

                term = term + Z_I * ( &
                     - utils_erf(x/sigma_I) / x2 &
                     + two_over_sqrt_pi/sigma_I*exp(-x2/(sigma_I*sigma_I)) / x)
             end if

             scalar_term = den_slabs12(i1,i2,islab12,1) * term
             locps_forces(1,I) = locps_forces(1,I) + scalar_term*r_minus_R_I%X/x
             locps_forces(2,I) = locps_forces(2,I) + scalar_term*r_minus_R_I%Y/x
             locps_forces(3,I) = locps_forces(3,I) + scalar_term*r_minus_R_I%Z/x

          else
             ! jd: Special case -- atom exactly on a gridpoint
             !     No-op, just avoid division by x=0.
          end if

       end do ! atoms, I

    end do ! ipt
!$OMP END DO
!$OMP END PARALLEL

    ! jd: Apply grid weight and sum over grid slabs
    locps_forces = locps_forces * grid%weight
    call comms_reduce('SUM', locps_forces)

    ! If spin polarised, restore up spin density
    if (pub_num_spins == 2) den_slabs12(:,:,:,1) = den_slabs12(:,:,:,1) - &
         den_slabs12(:,:,:,2)

    ! Stop timer
    call timer_clock('openbc_locps_calc_forces',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving openbc_locps_calc_forces'

  end subroutine openbc_locps_calc_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module openbc_locps

