! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Chris-Kriton Skylaris, Arash A. Mostofi and Nicholas D.M. Hine
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!   This module was spun off in July 2015 from previously-existing
!   routines within density_mod.F90, by Nicholas Hine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module density_init

  use constants, only: DP
  use ion, only: RADIAL_DENSITY_TYPE
  use rundat, only: pub_debug_on_root

  implicit none

  private

  public :: density_init_guess_recip
  public :: density_init_guess_real
  public :: density_init_radial_init
  public :: density_init_radial_store
  public :: density_init_radial_bcast
  public :: density_init_radial_exit
  public :: density_init_check_possemidef
  public :: density_init_render_possemidef
  public :: density_on_grid_renorm

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_init_guess_recip(density, & ! output
       regions, struct_fac, grid)

    !===============================================================!
    ! This subroutine initialises an approximate (guess) charge     !
    ! density as a superposition of Gaussian functions on atoms.    !
    !---------------------------------------------------------------!
    ! Written by Arash A. Mostofi on 19/2/2004.                     !
    ! Modified by Chris-Kriton Skylaris on 22/2/2004 to work        !
    ! with data-parallel structure factor.                          !
    ! Rewritten by Peter D. Haynes on 30/6/2004 to use new Fourier  !
    ! parallelisation.                                              !
    ! Modified by Nicholas Hine on 19/09/2008 for rearranged        !
    ! structure factor array to improve cache performance           !
    ! Moved to density_mod, and dependence on p_species array       !
    ! removed so that it only uses elements array, by Nicholas Hine !
    ! on 04/11/09.                                                  !
    ! Adjusted to allow execution where number of 23-slabs is less  !
    ! than number of procs by Nicholas Hine, December 2009.         !
    !===============================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_proc_id
    use constants, only: UP, DN, DP, stdout
    use fourier, only: fourier_apply_cell_backward
    use integrals, only: integrals_trace_on_grid
    use ion, only: ELEMENT
    use model_type, only: REGION
    use rundat, only: pub_num_spins
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)   :: grid
    type(REGION), intent(in)      :: regions
    complex(kind=DP), intent(in)  :: struct_fac(regions%par%num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)
    real(kind=DP), intent(out) :: density(grid%ld1,grid%ld2, &
         grid%max_slabs12, pub_num_spins)

    ! Local variables
    integer :: ierr                            ! Error flag
    real(kind=DP), parameter :: rad_default = 1.8_DP ! Default core radius
    integer :: i3,i2,islab23,is
    integer :: atom,species
    real(kind=DP) :: ion_charge
    real(kind=DP) :: g(3),gsq,dens_value,maxrad,intdens
    real(kind=DP),allocatable :: expon(:),pref(:)
    complex(kind=DP), allocatable :: dens_recip(:,:,:)

    ! Start timer
    call timer_clock('density_init_guess_recip',1)

    ! jd: Takes care of padding between n1 and ld1 etc.
    density = 0.0_DP

    ! Allocate complex workspace for reciprocal space
    allocate(dens_recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('density_init_guess_recip','dens_recip',ierr)
    allocate(expon(regions%par%num_pspecies),stat=ierr)
    call utils_alloc_check('density_init_guess_recip','expon',ierr)
    allocate(pref(regions%par%num_pspecies),stat=ierr)
    call utils_alloc_check('density_init_guess_recip','pref',ierr)

    ! ndmh: Set up gaussian exponents and prefactors
    ! ndmh: Loop over atomic species
    do species=1,regions%par%num_pspecies

       ! Find first example of this atom in elements array
       do atom=1,regions%par%nat+1
          call utils_assert (atom<=regions%par%nat, &
               'Error in density_init_guess_recip: no atoms of at least one &
               &species were found in the elements array. Species index: ', &
               species)
          if (regions%elements(atom)%pspecies_number==species) exit
       end do

       ! Find ion charge from elements array
       ion_charge = regions%elements(atom)%ion_charge

       ! If projector radii are defined for this species, find maximum radius
       ! amongst all shells, else use default
       maxrad = regions%elements(atom)%max_core_radius
       if (maxrad < 0.01_DP) maxrad = rad_default

       ! Exponent of gaussian
       expon(species) = maxrad * maxrad * 0.25_DP
       expon(species) = expon(species) * 1.75_DP ! fudge factor

       ! Prefactor, normalised to ion charge of species
       pref(species) = ion_charge / grid%weight

    end do ! loop over species

    ! Loop over reciprocal space grid on this proc
    do islab23=1,grid%num_slabs23          ! along b1

       dens_recip(:,:,islab23) = (0.0_DP,0.0_DP)

       do i2=1,grid%n2                     ! along b2
          do i3=1,grid%n3                  ! along b3

             ! ndmh: loop over species
             do species=1,regions%par%num_pspecies

                call cell_grid_recip_pt(g,islab23 + &
                     grid%first_slab23(pub_my_proc_id) - 1,i2,i3,grid)

                gsq = sum(g(:)**2)

                dens_value = pref(species) * exp(-expon(species)*gsq)

                dens_recip(i3,i2,islab23) = dens_recip(i3,i2,islab23) + &
                     struct_fac(species,i3,i2,islab23) * dens_value

             end do ! loop over species

          end do  ! loop along b3
       end do     ! loop along b2
    end do        ! loop along b1


    ! G=0 element must be real (find first 23-slab)
    if (pub_my_proc_id==grid%proc_slab23(1)) then
       if (aimag(dens_recip(1,1,1)) /= 0.0_DP) then
          write(stdout,'(a/a)') 'WARNING in density_init_guess_recip:', &
               '  density not real - setting imaginary part of G=0 term to zero'
          dens_recip(1,1,1) = cmplx(real(dens_recip(1,1,1),kind=DP),0.0_DP, &
               kind=DP)
       end if
    end if

    ! Nyquist filter (fine grid is always going to be even)
    if (grid%num_slabs23>0) then
       dens_recip(grid%n3/2+1,:,:) = (0.0_DP,0.0_DP)
       dens_recip(:,grid%n2/2+1,:) = (0.0_DP,0.0_DP)
    end if

    ! Find last 23-slab and apply Nyquist filter
    if (pub_my_proc_id==grid%proc_slab23(grid%n1/2+1)) &
         dens_recip(:,:,grid%num_slabs23) = (0.0_DP,0.0_DP)

    ! FFT the local ionic potential from reciprocal to real space
    call fourier_apply_cell_backward(density(:,:,:,1),dens_recip,grid)

    ! Deallocate workspace
    deallocate(pref,stat=ierr)
    call utils_dealloc_check('density_init_guess_recip','pref',ierr)
    deallocate(expon,stat=ierr)
    call utils_dealloc_check('density_init_guess_recip','expon',ierr)
    deallocate(dens_recip,stat=ierr)
    call utils_dealloc_check('density_init_guess_recip','dens_recip',ierr)

    ! Density may have negative values
    call density_init_render_possemidef(density(:,:,:,1),grid)

    ! ndmh: copy up density to down density
    if (pub_num_spins == 2) then
       density(:,:,:,UP) = density(:,:,:,UP)*0.5_DP
       density(:,:,:,DN) = density(:,:,:,UP)
    end if

    ! Renormalise density (incorporating net spin)
    intdens = integrals_trace_on_grid(density(:,:,:,1),grid)
    do is=1,pub_num_spins
       call density_on_grid_renorm(density(:,:,:,is), &
            grid,intdens,regions%elements,is)
    end do

    ! Stop timer
    call timer_clock('density_init_guess_recip',2)

  end subroutine density_init_guess_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_init_guess_real(density_on_grid, & ! output
       mdl, grid, add_aug_den,&                         ! input
       active_density_on_grid)                          ! optional output

    !===============================================================!
    ! This subroutine initialises an approximate (guess) charge     !
    ! density as a superposition of Gaussian functions on atoms.    !
    !---------------------------------------------------------------!
    ! Written by Nicholas Hine in May 2011.                         !
    ! OpenMP parallelisation by Nicholas Hine in June 2014          !
    ! Edited for embedding by Robert Charlton, 30/01/2017.          !
    ! Edited for EMFT embedding by Joseph Prentice, May 2018        !
    !===============================================================!

    use cell_grid, only: GRID_INFO, cell_grid_deposit_box, &
         cell_grid_box_start_wrt_atom
    use comms, only: pub_my_proc_id, pub_total_num_procs
    use constants, only: DP, UP, DN, ANGSTROM, max_spins
    use integrals, only: integrals_trace_on_grid
    use ion, only: ELEMENT
    use model_type, only: MODEL
    use rundat, only: pub_aug, pub_devel_code, pub_num_spins, pub_emft, &
         pub_active_region
!$  use rundat, only: pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert
    use visual, only: visual_scalarfield

    implicit none

    ! Arguments
    type(MODEL), intent(inout)    :: mdl
    type(GRID_INFO), intent(in)   :: grid
    real(kind=DP), intent(out) :: density_on_grid(grid%ld1,grid%ld2, &
         grid%max_slabs12, pub_num_spins)
    logical, intent(in) :: add_aug_den
    ! jcap: optional argument for EMFT
    real(kind=DP), intent(inout), optional :: active_density_on_grid(grid%ld1,&
         grid%ld2,grid%max_slabs12, pub_num_spins)

    ! Local variables
    integer :: ierr
    integer :: iat, loc_iat
    integer :: isp
    integer :: is
    integer :: box_n1, box_n2, box_n3
    integer :: box_start1, box_start2, box_start3
    integer :: proc
    integer :: npts
    logical :: i_have_box
    real(kind=DP) :: intdens(max_spins),scale
    real(kind=DP) :: rmax
    real(kind=DP),allocatable :: density_box(:,:,:), buffer(:,:,:)
    real(kind=DP),allocatable :: rad_dens(:)
    integer :: isub
    logical :: loc_emft

    ! Start timer
    call timer_clock('density_init_guess_real',1)

    ! jcap: check to see if we are doing an emft calculation and if
    ! this matches pub_emft
    loc_emft=present(active_density_on_grid)
    call utils_assert(loc_emft.eqv.pub_emft,'Arguments of &
         &density_init_guess_real inconsistent with pub_emft:',pub_emft)

    ! jd: Takes care of padding between n1 and ld1 etc.
    density_on_grid = 0.0_DP

    ! ndmh: Loop over atomic species to check densities have been initialised
    do isub=1,mdl%nsub
        do isp=1,mdl%regions(isub)%par%num_species

            proc = modulo(isp,pub_total_num_procs)
            if (proc==pub_my_proc_id) then

                ! ndmh: check if we already have a density for this element
                call utils_assert(mdl%regions(isub)%radial_densities(isp)%present, &
                    'Error in density_init_guess_real: atom density for at least &
                    &one species has not been initialised. All species must have &
                    &atom density calculated to use real-space density calculation.')
            end if

        end do ! loop over species
    end do ! loop over regions

    ! Pick a box size: use tightbox size scaled to this grid
    scale = real(grid%n1,kind=DP)/real(mdl%cell%total_pt1,kind=DP)
    box_n1 = min(int(mdl%uni_tightbox%total_pt1*scale + 1),grid%n1)
    scale = real(grid%n2,kind=DP)/real(mdl%cell%total_pt2,kind=DP)
    box_n2 = min(int(mdl%uni_tightbox%total_pt2*scale + 1),grid%n2)
    scale = real(grid%n3,kind=DP)/real(mdl%cell%total_pt3,kind=DP)
    box_n3 = min(int(mdl%uni_tightbox%total_pt3*scale + 1),grid%n3)

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(loc_iat,iat,isp,rmax,box_start1,box_start2,box_start3, &
!$OMP      density_box,buffer,rad_dens,is,i_have_box,ierr,npts) &
!$OMP SHARED(box_n1,box_n2,box_n3,grid,pub_my_proc_id, &
!$OMP      density_on_grid,add_aug_den, &
!$OMP      pub_aug,mdl,pub_threads_num_fftboxes, &
!$OMP      pub_num_spins,loc_emft,active_density_on_grid,pub_active_region)

    ! Allocate workspace
    allocate(density_box(box_n1,box_n2,box_n3),stat=ierr)
    call utils_alloc_check('density_init_guess_real','density_box',ierr)
    allocate(buffer(box_n1,box_n2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('density_init_guess_real','buffer',ierr)

    ! Loop over atoms on this proc
    do isub=1,mdl%nsub
       ! rc2013: allocate rad_dens separately for each region
       allocate(rad_dens(maxval(mdl%regions(isub)%radial_densities(:)%npts)),stat=ierr)
       call utils_alloc_check('density_init_guess_real','rad_dens',ierr)
       do is=1,pub_num_spins
!$OMP DO
          do loc_iat=1,mdl%regions(isub)%par%max_atoms_on_proc

             if (loc_iat<=mdl%regions(isub)%par%num_atoms_on_proc(pub_my_proc_id)) then
                iat = mdl%regions(isub)%par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
                isp = mdl%regions(isub)%par%elements_on_proc(loc_iat)%species_number
                npts = mdl%regions(isub)%radial_densities(isp)%npts
                rmax = min(mdl%regions(isub)%radial_densities(isp)%rad(npts), &
                     mdl%regions(isub)%par%elements_on_proc(loc_iat)%radius)

                ! Find where box for this atom is located in simulation cell
                call cell_grid_box_start_wrt_atom( &
                     box_start1, box_start2, box_start3, &
                     mdl%regions(isub)%par%elements_on_proc(loc_iat)%centre, box_n1, box_n2, box_n3, &
                     grid,mdl%cell)

                ! Reset box to zero
                density_box(:,:,:) = 0.0_DP

                ! Set radial density profile
                rad_dens(:) = 0.0_DP
                rad_dens(1:npts) = mdl%regions(isub)%radial_densities(isp)%den(1:npts,is)

                ! Add in augmentation density if required
                if (pub_aug.and.add_aug_den) rad_dens(1:npts) = &
                     rad_dens(1:npts) + mdl%regions(isub)%radial_densities(isp)%aug_den(1:npts,is)

                ! Transfer the radial density to the 3D box
                call internal_density_to_box(density_box(:,:,:), &
                     rad_dens(1:npts),mdl%regions(isub)%radial_densities(isp)%rad, &
                     rmax,mdl%regions(isub)%radial_densities(isp)%npts, &
                     box_n1,box_n2,box_n3,box_start1,box_start2,box_start3,grid, &
                     mdl%regions(isub)%par%elements_on_proc(loc_iat)%centre)

                ! Scale by 1/2 if spin-polarised (atom solver density is
                ! always non-polarised)
                density_box(:,:,:) = density_box(:,:,:)

                i_have_box = .true.
             else
                ! Nothing to deposit on this proc
                i_have_box = .false.
             end if

             ! Deposit this box to the simulation cell if present, or just wait for
             ! data from other procs if no box
!$OMP CRITICAL
             call cell_grid_deposit_box(density_on_grid(:,:,:,is), &
                  density_box(:,:,:), buffer, grid, &
                  box_n1, box_n2, box_n3, box_n1, box_n2, &
                  box_start1, box_start2, box_start3, i_have_box, .false.)
             if (loc_emft.and.(isub==pub_active_region)) then
                ! rc2013: no need to recalculate this; just copy across
                active_density_on_grid(:,:,:,is) = density_on_grid(:,:,:,is)
             end if

!$OMP END CRITICAL
          end do
!$OMP END DO
!$OMP BARRIER
       end do
       ! rc2013: deallocate rad_dens before changing regions
       deallocate(rad_dens,stat=ierr)
       call utils_dealloc_check('density_init_guess_real','rad_dens',ierr)
    end do ! rc2013: end of regions loop

    ! Deallocate
    deallocate(density_box,stat=ierr)
    call utils_dealloc_check('density_init_guess_real','density_box',ierr)
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('density_init_guess_real','buffer',ierr)

!$OMP END PARALLEL

    ! Density may have negative values
    do is=1,pub_num_spins
       call density_init_render_possemidef(density_on_grid(:,:,:,is),grid)
       if (loc_emft) call density_init_render_possemidef(active_density_on_grid(:,:,:,is),&
            grid)
    end do

    ! Renormalise density
    do is=1,pub_num_spins
       intdens(is) = integrals_trace_on_grid(density_on_grid(:,:,:,is),grid)
       if (add_aug_den) call density_on_grid_renorm(density_on_grid(:,:,:,is), &
            grid,intdens(is),mdl%elements,is)
       if (loc_emft) then
          intdens(is) = integrals_trace_on_grid(active_density_on_grid(:,:,:,is),grid)
          if (add_aug_den) call density_on_grid_renorm(active_density_on_grid(:,:,:,is), &
               grid,intdens(is),mdl%regions(pub_active_region)%elements,is)
       end if
    end do

    ! ndmh: show density if required
    if (index(pub_devel_code,'WRITE_DENS_GUESS')>0) then
       if (pub_num_spins==2) then
          ! set spin up density to be total density
          density_on_grid(:,:,:,UP) = density_on_grid(:,:,:,UP) + &
               density_on_grid(:,:,:,DN)
       end if
       if (pub_aug.and.add_aug_den) then
          call visual_scalarfield(density_on_grid,grid,mdl%cell, &
               'Guess density (in e/ang^3) for:', '_guess_density_aug', &
               mdl%elements, ANGSTROM**3)
       else
          call visual_scalarfield(density_on_grid,grid,mdl%cell, &
               'Guess density (in e/ang^3) for:', '_guess_density', &
               mdl%elements, ANGSTROM**3)
       end if
       if (pub_num_spins==2) then
          ! restore spin up density
          density_on_grid(:,:,:,UP) = density_on_grid(:,:,:,UP) - &
               density_on_grid(:,:,:,DN)
          ! set spin down density to be spin density
          density_on_grid(:,:,:,DN) = density_on_grid(:,:,:,UP) - &
               density_on_grid(:,:,:,DN)
          if (pub_aug.and.add_aug_den) then
             call visual_scalarfield(density_on_grid(:,:,:,DN),grid,mdl%cell, &
                  'Guess density (in e/ang^3) for:', '_guess_spindensity_aug', &
                  mdl%elements, ANGSTROM**3)
          else
             call visual_scalarfield(density_on_grid(:,:,:,DN),grid,mdl%cell, &
                  'Guess density (in e/ang^3) for:', '_guess_spindensity', &
                  mdl%elements, ANGSTROM**3)
          end if
          density_on_grid(:,:,:,DN) = -density_on_grid(:,:,:,DN) + &
               density_on_grid(:,:,:,UP)
       end if
    end if

    ! Stop timer
    call timer_clock('density_init_guess_real',2)

contains

    subroutine internal_density_to_box(box,den,r,rcut,npts, &
         box_n1,box_n2,box_n3,cell_start1,cell_start2,cell_start3,grid, &
         atom_origin)

      use basis, only: basis_box_origin_to_atom
      use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
           OPERATOR(*), geometry_magnitude
      use services, only: services_locate_interp,services_linear_interpolation

      ! Arguments
      integer,intent(in) :: box_n1, box_n2, box_n3
      real(kind=DP),intent(inout) :: box(box_n1,box_n2,box_n3)
      integer,intent(in) :: npts
      real(kind=DP), intent(in) :: den(npts)
      real(kind=DP), intent(in) :: r(npts), rcut
      integer,intent(in) :: cell_start1
      integer,intent(in) :: cell_start2
      integer,intent(in) :: cell_start3
      type(GRID_INFO),intent(in) :: grid
      type(POINT),intent(in) :: atom_origin

      ! Local Variables
      integer :: i1,i2,i3
      integer :: ipt
      type(POINT) :: box_origin
      type(POINT) :: box_to_atom
      type(POINT) :: r_cell, r_sphere
      real(kind=DP) :: rmag

      call basis_box_origin_to_atom(box_to_atom,box_origin,atom_origin, &
           cell_start1,cell_start2,cell_start3,grid%da1,grid%da2,grid%da3, &
           mdl%cell)

      do i3=1,box_n3
         do i2=1,box_n2
            do i1=1,box_n1

               r_cell = box_origin + (i1-1)*grid%da1 &
                                   + (i2-1)*grid%da2 &
                                   + (i3-1)*grid%da3

               r_sphere = r_cell - atom_origin
               rmag = geometry_magnitude(r_sphere)

               if (rmag>rcut) then
                  box(i1,i2,i3) = 0.0_DP
               else
                  ipt = services_locate_interp(rmag,r,npts)
                  box(i1,i2,i3) = services_linear_interpolation(rmag,den(ipt), &
                       den(ipt+1),r(ipt),r(ipt+1))
               end if

            end do
         end do
      end do

    end subroutine internal_density_to_box

  end subroutine density_init_guess_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_init_radial_init(regions)

    !=================================================================!
    ! This subroutine initialises storage for radial densities        !
    !-----------------------------------------------------------------!
    ! Written by Nicholas D.M. Hine in May 2011.                      !
    ! Edited for embedding by Robert Charlton, 08/12/2016.            !
    !=================================================================!

    use model_type, only: REGION
    use rundat, only: pub_paw
    use utils, only: utils_alloc_check, utils_assert

    ! Arguments
    type(REGION), intent(inout) :: regions

    ! Local Variables
    integer :: ierr
    integer :: isp, ipsp, iat

    ! Allocate storage for all species
    allocate(regions%radial_densities(regions%par%num_species),stat=ierr)
    call utils_alloc_check('density_init_radial_init','regions%radial_densities',ierr)
    regions%radial_densities(:)%present = .false.

    do isp=1,regions%par%num_species

       ! ndmh: Find pspecies number for this species
       ipsp = -1
       do iat=1,regions%par%nat
          if (regions%elements(iat)%species_number==isp) then
             ipsp = regions%elements(iat)%pspecies_number
             exit
          end if
       end do
       ! rc2013: need to be careful with allocation of species labels...
       if(iat > regions%par%nat) cycle

       call utils_assert (ipsp>0, &
            'Error in density_init_radial_init: no atoms of at least one &
            &species were found in the elements array. Species index: ', &
            isp)
       if (pub_paw) then
          regions%radial_densities(isp)%npwtot = regions%paw_sp(ipsp)%npw_tot
          regions%radial_densities(isp)%ipsp   = ipsp
       end if
    end do

  end subroutine density_init_radial_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_init_radial_store(regions,isp,npts,rad,den,aug_den,rhoij0)

    !=================================================================!
    ! This subroutine stores a radial density for a given species, to !
    ! be used later in construction of the initial guess density      !
    !-----------------------------------------------------------------!
    ! Written by Nicholas D.M. Hine in May 2011.                      !
    ! Edited for embedding by Robert Charlton, 08/12/2016.            !
    !=================================================================!

    use constants, only: PI
    use model_type, only: REGION
    use paw, only: paw_store_initial_proj_kern
    use rundat, only: pub_aug, pub_paw, pub_num_spins
    use utils, only: utils_abort

    ! Arguments
    type(REGION), intent(inout) :: regions
    integer, intent(in) :: isp
    integer, intent(in) :: npts
    real(kind=DP), intent(in) :: rad(:)
    real(kind=DP), intent(in) :: den(:,:)
    real(kind=DP), intent(in) :: aug_den(:,:)
    real(kind=DP), intent(in) :: rhoij0(:,:,:)

    ! Local Variables
    integer :: ierr
    integer :: is
    character(len=10) :: isp_str

    ! Check sanity of arguments
    if ((isp<1).or.(isp>regions%par%num_species)) then
       write(isp_str,'(i9)') isp
       call utils_abort('Error in density_init_radial_store: invalid species &
            &number supplied: '//adjustl(trim(isp_str)))
    end if
    if ((npts<0).or.(npts>50000)) then
       write(isp_str,'(i9)') isp
       call utils_abort('Error in density_init_radial_store: invalid radial &
            &grid size for species '//adjustl(trim(isp_str)))
    end if
    if ((size(rad)<npts).or.(size(den,1)<npts).or.(size(aug_den,1)<npts)) then
       write(isp_str,'(i9)') isp
       call utils_abort('Error in density_init_radial_store: invalid array  &
            &sizes passed in for species '//adjustl(trim(isp_str)))
    end if

    call internal_allocate(regions%radial_densities(isp))

    if (pub_aug) then
       if (pub_paw) then
          ! ndmh: call paw_mod routine to expand rhoij0 to full m-including
          ! ndmh: array
          call paw_store_initial_proj_kern(regions%radial_densities,rhoij0, &
               isp,regions%paw_sp)
       end if
    end if

    ! Copy in arguments
    regions%radial_densities(isp)%npts = npts
    regions%radial_densities(isp)%rad(1:npts) = rad(1:npts)
    do is=1,pub_num_spins
       regions%radial_densities(isp)%den(1:npts,is) = den(1:npts,is) / 4.0_DP / PI
       if (pub_aug) regions%radial_densities(isp)%aug_den(1:npts,is) = &
            aug_den(1:npts,is) / 4.0_DP / PI
    end do

    ! Set flag to indicate this density exists
    regions%radial_densities(isp)%present = .true.

contains

    subroutine internal_allocate(rad_den)

      use utils, only: utils_alloc_check

      ! Arguments
      type(RADIAL_DENSITY_TYPE), intent(inout) :: rad_den

      ! Local Variables
      integer :: npwtot

      ! Allocate storage for this species
      allocate(rad_den%rad(npts),stat=ierr)
      call utils_alloc_check('density_init_radial_store', &
           'rad_den%rad',ierr)
      allocate(rad_den%den(npts,pub_num_spins),stat=ierr)
      call utils_alloc_check('density_init_radial_store', &
           'rad_den%den',ierr)
      if (pub_aug) then
         allocate(rad_den%aug_den(npts,pub_num_spins),stat=ierr)
         call utils_alloc_check('density_init_radial_store', &
              'rad_den%aug_den',ierr)
         if (pub_paw) then
            npwtot = rad_den%npwtot
            allocate(rad_den%rhoij0(npwtot,npwtot,pub_num_spins),stat=ierr)
            call utils_alloc_check('density_init_radial_store', &
                 'rad_den%rhoij0',ierr)
         end if
      end if

    end subroutine internal_allocate

  end subroutine density_init_radial_store


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_init_radial_bcast(rad_den,proc)

    !=================================================================!
    ! This subroutine stores a radial density for a given species, to !
    ! be used later in construction of the initial guess density      !
    !-----------------------------------------------------------------!
    ! Written by Nicholas D.M. Hine in May 2011.                      !
    !=================================================================!

    use comms, only: comms_bcast
    use rundat, only: pub_aug, pub_paw, pub_num_spins
    use utils, only: utils_alloc_check

    ! Arguments
    type(RADIAL_DENSITY_TYPE), intent(inout) :: rad_den
    integer, intent(in) :: proc

    ! Local Variables
    integer :: ierr

    ! Broadcast size of arrays
    call comms_bcast(proc,rad_den%npts)

    ! Allocate storage for this species
    if (.not.rad_den%present) then
       allocate(rad_den%rad(rad_den%npts),stat=ierr)
       call utils_alloc_check('density_init_radial_init', &
            'rad_den%rad',ierr)
       allocate(rad_den%den(rad_den%npts,pub_num_spins),stat=ierr)
       call utils_alloc_check('density_init_radial_init', &
            'rad_den%den',ierr)
       if (pub_aug) then
          allocate(rad_den%aug_den(rad_den%npts,pub_num_spins),stat=ierr)
          call utils_alloc_check('density_init_radial_init', &
               'rad_den%aug_den',ierr)
          if (pub_paw) then
             allocate(rad_den%rhoij0(rad_den%npwtot, &
                  rad_den%npwtot, pub_num_spins), stat=ierr)
             call utils_alloc_check('density_init_radial_init', &
                 'rad_den%rhoij0',ierr)
          end if
       end if
    end if

    ! Broadcast density and radial grid
    call comms_bcast(proc,rad_den%rad)
    call comms_bcast(proc,rad_den%den)
    if (pub_aug) call comms_bcast(proc,rad_den%aug_den)
    if (pub_paw) call comms_bcast(proc,rad_den%rhoij0)

    ! Set flag to indicate this density exists
    rad_den%present = .true.

  end subroutine density_init_radial_bcast


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_init_radial_exit(regions)

    !=================================================================!
    ! This subroutine deallocates storage for radial densities        !
    !-----------------------------------------------------------------!
    ! Written by Nicholas D.M. Hine in May 2011.                      !
    !=================================================================!

    use model_type, only: REGION
    use rundat, only: pub_aug, pub_paw
    use utils, only: utils_dealloc_check

    ! Arguments
    type(REGION), intent(inout) :: regions

    ! Local Variables
    integer :: isp
    integer :: ierr

    do isp=1,regions%par%num_species
       if (regions%radial_densities(isp)%present) then
          call internal_deallocate(regions%radial_densities(isp))
       end if
    end do

    deallocate(regions%radial_densities,stat=ierr)
    call utils_dealloc_check('density_init_radial_exit', &
         'regions%radial_densities',ierr)

contains

    subroutine internal_deallocate(rad_den)

      ! Arguments
      type(RADIAL_DENSITY_TYPE), intent(inout) :: rad_den

      if (pub_aug) then
         if (pub_paw) then
            deallocate(rad_den%rhoij0,stat=ierr)
            call utils_dealloc_check('density_init_radial_exit', &
                 'rad_den%rhoij0',ierr)
         end if
         deallocate(rad_den%aug_den,stat=ierr)
         call utils_dealloc_check('density_init_radial_exit', &
              'rad_den%aug_den',ierr)
      end if
      deallocate(rad_den%den,stat=ierr)
      call utils_dealloc_check('density_init_radial_exit', &
           'rad_den%den',ierr)
      deallocate(rad_den%rad,stat=ierr)
      call utils_dealloc_check('density_init_radial_exit', &
           'rad_den%rad',ierr)

    end subroutine internal_deallocate

  end subroutine density_init_radial_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_init_check_possemidef(density_fine,grid)

    !================================================================!
    ! This subroutine is checks if a given density on the fine grid  !
    ! has negative values and stops if does.                         !
    !----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.                      !
    !================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: DP
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(in) :: density_fine(grid%ld1,&
         grid%ld2,grid%n3)

    ! cks: internal definitions
    integer :: r1, r2, r3


    do r3=1,grid%n3
       do r2=1,grid%n2
          do r1=1,grid%n1

             if ( density_fine(r1,r2,r3).lt.(0.0_DP)) then
                call utils_abort('Negative density value detected in &
                     &density_check_possemidef(). Indices of grid points &
                     &and density value follow.',r1,r2,r3, &
                     opt_real_to_print1=density_fine(r1,r2,r3))
             endif

          enddo
       enddo
    enddo


  end subroutine density_init_check_possemidef


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_init_render_possemidef(density,grid)

    !=======================================================================!
    ! This subroutine renders the density of the fine grid positive         !
    ! semidefinite by zeroing out all negative values                       !
    !-----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in spring 2001.                      !
    ! Modified by Chris-Kriton Skylaris on 18/7/2004 so that it works with  !
    ! the parallel version of ONETEP.                                       !
    !=======================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: DP

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: density(grid%ld1,grid%ld2,grid%max_slabs12)

    ! cks: internal declarations
    integer :: row1, row2, islab12

    do islab12=1,grid%num_my_slabs12
       do row2=1,grid%n2
          do row1=1,grid%n1
             density(row1,row2,islab12) = &
                  max(density(row1,row2,islab12),0.0_DP)
          end do
       end do
    end do

  end subroutine density_init_render_possemidef


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_on_grid_renorm(density,grid,integrated_density,elements, &
       spin_channel)

    !======================================================================!
    ! This subroutine renormalises the density on the fine grid to the     !
    ! correct number of electrons.                                         !
    !----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2001.                            !
    ! Modified on 1/9/2004 by Chris-Kriton Skylaris to add external charge !
    !======================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use ion, only: element
    use rundat, only: pub_charge, pub_spin_polarised, pub_edft, pub_spin, pub_real_spin

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: density(grid%ld1,grid%ld2, &
         grid%max_slabs12)
    ! rc2013: remove explicit par references, extract the array size if needed
    type(ELEMENT), intent(in) :: elements(:) !size par%nat)
    real(kind=DP), intent(in) :: integrated_density
    integer, intent(in) :: spin_channel

    ! cks: internal declarations
    integer :: row, ne
    real(kind=DP) :: ne_real, spin_real
    real(kind=DP) :: factor

    ! ndmh: find spin as real
    spin_real = real(pub_spin,kind=DP)
    if (pub_edft) spin_real = pub_real_spin

    ! cks: find the total number of electrons in the simulation cell
    ! ndmh: add as reals
    ne_real = 0
    do row=1,size(elements)
       ne_real = ne_real + elements(row)%ion_charge
    end do

    if (pub_edft) then
      spin_real = pub_real_spin
    else
      spin_real = real(pub_spin, kind=dp)
    end if

    ! ndmh: calculate final target ne
    ! ndmh: spin polarisation
    if (pub_spin_polarised) then
       if (spin_channel==1) then
          ne_real = (ne_real - pub_charge + spin_real)*0.5_DP
       else
          ne_real = (ne_real - pub_charge - spin_real)*0.5_DP
       end if
    else
       ne_real = ne_real - pub_charge
    end if

    if (.not. pub_edft) ne_real = real(int(ne_real),kind=DP)

    ! cks: renormalise density
    factor = ne_real / integrated_density

    ! cks: renormalise density
    if (pub_on_root.and.((factor>1.01_DP).or.(factor<0.99_DP))) then
       if (pub_spin_polarised) then
          write(stdout,'(a,i2,a)') 'Density renormalisation, spin ', &
               spin_channel, ' ...'
       else
          write(stdout,'(a)') 'Density renormalisation ...'
       end if
       write(stdout,'(2(a,f18.12))') 'Target Charge: ',ne_real, &
            ' ; Input Charge: ', integrated_density
    end if

    density = factor * density

  end subroutine density_on_grid_renorm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module density_init

