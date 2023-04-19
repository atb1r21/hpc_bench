! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*- !
!================================================================!
!                                                                !
!                         Phonons Module                         !
!                                                                !
! This module calculates the Gamma-point vibrational frequencies !
! of the cell using the finite displacement method.              !
!----------------------------------------------------------------!
! Written by Fabiano Corsetti (May 2011)                         !
!================================================================!

module phonon

  use constants, only: DP

  implicit none

  private

  ! Conversion constants
  real(kind=DP), parameter :: electron_mass_u =5.4857990943E-4_DP
  real(kind=DP), parameter :: au2THz = 1.0_DP/2.41888468d-5
  real(kind=DP), parameter :: au2inv_cm = au2THz/2.99792458d-2
  ! Boltzmann constant in Eh/K
  real(kind=DP), parameter :: k_B = 3.1668115744561575d-06

  public :: phonon_main

contains

  subroutine phonon_main(total_energy,forces,mdl)

    !==================================================================!
    ! This subroutine calculates and outputs the phonon frequencies of !
    ! a system using a finite-displacement method.                     !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  total_energy (inout)  : Total energy of the system.             !
    !  forces (inout)        : Total force on each ion.                !
    !  mdl (inout)           : Model, containing elements list.        !
    !------------------------------------------------------------------!
    ! Written by Fabiano Corsetti in June 2011.                        !
    ! Reformatted by Nicholas Hine in October 2011.                    !
    ! Extra functionality added by Fabiano Corsetti in July 2013.      !
    ! Modified for embedding structures by Joseph Prentice, May 2018   !
    !==================================================================!

    use bibliography, only: bibliography_cite
    use comms, only: comms_bcast, pub_on_root, pub_root_proc_id
    use dense, only: DEM, dense_create, dense_destroy, dense_put_element, &
         dense_get_element, dense_eigensolve
    use constants, only: DP, TWO_PI, stdout, periodic_table_mass, ANGSTROM
    use energy_and_force, only: energy_and_force_calculate
    use model_type, only: MODEL
    use services, only: services_rationalise_coords
    use utils, only: utils_unit, utils_alloc_check, utils_dealloc_check, &
        utils_qc_print, utils_assert
    use rundat, only: pub_rootname, pub_write_denskern, pub_write_tightbox_ngwfs, &
         pub_read_denskern, pub_read_tightbox_ngwfs, pub_phonon_disp, &
         pub_phonon_fmax, pub_phonon_farming_task, pub_phonon_sampling, &
         pub_phonon_vib_free, pub_have_phonon_disp_list, &
         pub_num_disp, pub_phonon_disp_list, pub_have_phonon_except_list, &
         pub_num_except, pub_phonon_iexcept_list, pub_phonon_dexcept_list, &
         pub_phonon_tmin, pub_phonon_tmax, pub_phonon_energy_check, &
         pub_phonon_deltat, pub_phonon_min_freq, pub_print_qc, &
         pub_phonon_write_eigenvecs, pub_have_phonon_animate_list, &
         pub_phonon_animate_list, pub_num_ani, pub_phonon_animate_scale, &
         pub_geom_reuse_dk_ngwfs, pub_maxit_ngwf_cg, pub_have_supercell, &
         pub_supercell, pub_nat_unit, pub_supercell_unit_list, &
         pub_phonon_DOS, pub_phonon_grid, pub_have_phonon_qpoints, &
         pub_num_qpoints, pub_phonon_qpoints, pub_phonon_SK, &
         pub_phonon_DOS_min, pub_phonon_DOS_max, pub_phonon_DOS_delta, &
         pub_eigensolver_orfac, pub_eigensolver_abstol

    implicit none

    ! Arguments
    type(MODEL), intent(inout) :: mdl
    real(kind=DP), intent(inout) :: total_energy
    real(kind=DP), intent(inout) :: forces(1:3,mdl%nat)

    ! Local Variables
    type(DEM) :: dynamical_dem
    type(DEM) :: dynamical_symm_dem
    type(DEM) :: eigenvecs_dem
    type(DEM) :: unit_mat_dem
    character(len=1) :: dir_char
    character(len=1) :: cart_char
    character(len=256) :: output_file
    character(len=256) :: disp_number
    logical :: do_disp
    logical :: converged
    logical :: print_warning, print_warning2
    logical :: phonon_vib_free(1:3)
    logical :: found_ion_unit, found_ion
    integer :: g1, g2, g3, i, j, k, l, j_min, k_min, l_min, m, n, o, p
    integer :: T_num_iter, tot_num_disp, ncells, nqdos
    integer :: grid_loop, tot_grid_loop, output_unit
    integer :: ierr
    integer, allocatable, dimension(:) :: tot_sampling_list, temp_list
    integer, allocatable, dimension(:,:) :: tot_disp_list, tot_disp_list_back
    integer, allocatable, dimension(:,:) :: supercell_list
    integer, allocatable, dimension(:,:,:) :: q_calc_list
    integer, allocatable, dimension(:,:,:,:) :: cutoff_translate
    real(kind=DP) :: unit_cell(1:3,1:3), rec_unit_cell(1:3,1:3), shift(1:3)
    real(kind=DP) :: dist, dist_min
    real(kind=DP) :: E0
    real(kind=DP) :: dist_fact
    real(kind=DP) :: Fmax
    real(kind=DP) :: eql_pos(1:3), disp_pos(1:3)
    real(kind=DP) :: freq_conv
    real(kind=DP) :: el2
    real(kind=DP) :: el2_T
    real(kind=DP) :: el2_symm
    real(kind=DP) :: mass_n, mass_m
    real(kind=DP) :: qdos_weight
    real(kind=DP) :: q(1:3), zero_point_E, T, U, beta, sqrt_freq, exp_freq
    real(kind=DP), allocatable, dimension(:) :: tot_dist_list, freqs, F, S, C_v
    real(kind=DP), allocatable, dimension(:) :: qdos
    real(kind=DP), allocatable, dimension(:,:) :: force_consts, ion_write2
    real(kind=DP), allocatable, dimension(:,:,:) :: force_consts_m
    real(kind=DP), allocatable, dimension(:,:,:) :: abs_q_calc_list

    complex(kind=DP) :: el
    complex(kind=DP) :: el_T
    complex(kind=DP) :: el_symm
    complex(kind=DP), allocatable, dimension(:,:) :: ion_write

    complex(kind=DP), parameter :: cmplx_1=(1.0_dp,0.0_dp)
    complex(kind=DP), parameter :: cmplx_i=(0.0_dp,1.0_dp)
    complex(kind=DP), parameter :: cmplx_0=(0.0_dp,0.0_dp)

    ! jcap: further local parameters for embedding
    integer :: iat,iregion,sub_i

    call bibliography_cite('SUPERCELL_PHONON')

    ! Display Banner
    if (pub_on_root) then
       write(stdout,'()')
       write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
            &-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
       write(stdout,'(a)') 'Starting phonon calculation...'
       write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
            &-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
       write(stdout,'()')
    end if

    ! Check inputs
    call utils_assert(&
         (pub_phonon_farming_task>=0) .and. (pub_phonon_farming_task<=3), &
         'Invalid value of phonon_farming_task: ', pub_phonon_farming_task)

    call utils_assert(&
         (pub_phonon_sampling==1) .or. (pub_phonon_sampling==2), &
         'Invalid value of phonon_sampling: ', pub_phonon_sampling)

    call utils_assert(&
         (pub_phonon_vib_free>=0) .and. (pub_phonon_vib_free<=7), &
         'Invalid value of phonon_vib_free: ', pub_phonon_vib_free)

    if (pub_have_phonon_except_list) then
       call utils_assert(&
            all(pub_phonon_iexcept_list(1:2,1:pub_num_except)>=1) .and. &
            all(pub_phonon_iexcept_list(3,1:pub_num_except)>=0) .and. &
            all(pub_phonon_iexcept_list(4,1:pub_num_except)>=1) .and. &
            all(pub_phonon_iexcept_list(1,1:pub_num_except)<=mdl%nat) &
            .and. all(pub_phonon_iexcept_list(2,1:pub_num_except)<=3) .and. &
            all(pub_phonon_iexcept_list(3,1:pub_num_except)<=1) .and. &
            all(pub_phonon_iexcept_list(4,1:pub_num_except)<=2), &
            'Invalid entry in phonon_exception_list')
    end if

    if (pub_have_phonon_except_list .and. pub_have_supercell) then
       do i=1,pub_num_except
          call utils_assert(&
               any(pub_phonon_iexcept_list(1,i)==&
               pub_supercell_unit_list(1:pub_nat_unit)), &
               'Invalid entry in phonon_exception_list')
       end do
    end if

    if (pub_phonon_energy_check .and. (pub_phonon_farming_task==2)) then
       if (pub_on_root) then
          write(stdout,'(a)') 'WARNING: phonon_energy_check will not be &
               &enabled unless phonon_farming_task is set to 0'
       end if
    end if

    if (pub_have_phonon_qpoints .and. (.not. pub_have_supercell)) then
       call utils_assert(&
            all(pub_phonon_qpoints(1:3,1:pub_num_qpoints)==0.0_DP), &
            'Only Gamma point allowed in phonon_qpoints when not using a &
            &supercell')
    end if


    ! Calculate unit cell vectors
    unit_cell(1:3,1)=(/mdl%cell%a1%x,mdl%cell%a1%y,mdl%cell%a1%z/)
    unit_cell(1:3,2)=(/mdl%cell%a2%x,mdl%cell%a2%y,mdl%cell%a2%z/)
    unit_cell(1:3,3)=(/mdl%cell%a3%x,mdl%cell%a3%y,mdl%cell%a3%z/)

    ! Calculate reciprocal unit cell vectors
    rec_unit_cell(1:3,1)=(/mdl%cell%b1%x,mdl%cell%b1%y,mdl%cell%b1%z/)
    rec_unit_cell(1:3,2)=(/mdl%cell%b2%x,mdl%cell%b2%y,mdl%cell%b2%z/)
    rec_unit_cell(1:3,3)=(/mdl%cell%b3%x,mdl%cell%b3%y,mdl%cell%b3%z/)

    ! Process supercell information
    ncells=pub_supercell(1)*pub_supercell(2)*pub_supercell(3)
    if (pub_have_supercell) then

       ! Check that input file information is consistent
       call utils_assert(&
            (mdl%nat==pub_nat_unit*ncells) .and. &
            all(pub_supercell(1:3)>=1), &
            'Incorrect supercell specifications')
       call utils_assert(&
            all(pub_supercell_unit_list(1:pub_nat_unit)>=1) .and. &
            all(pub_supercell_unit_list(1:pub_nat_unit)<=mdl%nat), &
            'Incorrect supercell specifications')
       do i=1,pub_nat_unit
          do j=i+1,pub_nat_unit
             call utils_assert(&
                  pub_supercell_unit_list(i)/=pub_supercell_unit_list(j), &
                  'Incorrect supercell specifications')
          end do
       end do

       ! Calculate unit cell vectors
       unit_cell(1:3,1)=unit_cell(1:3,1)/&
            pub_supercell(1)
       unit_cell(1:3,2)=unit_cell(1:3,2)/&
            pub_supercell(2)
       unit_cell(1:3,3)=unit_cell(1:3,3)/&
            pub_supercell(3)

       ! Calculate reciprocal unit cell vectors
       rec_unit_cell(1:3,1)=rec_unit_cell(1:3,1)*&
            pub_supercell(1)
       rec_unit_cell(1:3,2)=rec_unit_cell(1:3,2)*&
            pub_supercell(2)
       rec_unit_cell(1:3,3)=rec_unit_cell(1:3,3)*&
            pub_supercell(3)

       ! Assign each ion in the supercell to a unit cell
       allocate(supercell_list(1:5,1:mdl%nat),stat=ierr)
       call utils_alloc_check('phonon','supercell_list',ierr)
       allocate(q_calc_list(1:5,1:ncells,1:pub_nat_unit),stat=ierr)
       call utils_alloc_check('phonon','q_calc_list',ierr)
       allocate(abs_q_calc_list(1:3,1:ncells,1:pub_nat_unit),stat=ierr)
       call utils_alloc_check('phonon','abs_q_calc_list',ierr)
       allocate(temp_list(1:pub_nat_unit),stat=ierr)
       call utils_alloc_check('phonon','temp_list',ierr)
       allocate(cutoff_translate(1:3,1:pub_supercell(1),1:pub_supercell(2),&
            1:pub_supercell(3)),stat=ierr)
       call utils_alloc_check('phonon','cutoff_translate',ierr)
       do j=1,pub_supercell(1)
       do k=1,pub_supercell(2)
       do l=1,pub_supercell(3)
          dist_min=99999.9_DP
          do m=-1,0
          do n=-1,0
          do o=-1,0
             shift(1:3)=unit_cell(1:3,1)*(j-1+m*(pub_supercell(1)))+&
                  unit_cell(1:3,2)*(k-1+n*(pub_supercell(2)))+&
                  unit_cell(1:3,3)*(l-1+o*(pub_supercell(3)))
             dist=shift(1)**2+shift(2)**2+shift(3)**2
             if (dist<dist_min) then
                dist_min=dist
                j_min=j-1+m*(pub_supercell(1))
                k_min=k-1+n*(pub_supercell(2))
                l_min=l-1+o*(pub_supercell(3))
             end if
          end do
          end do
          end do
          cutoff_translate(1:3,j,k,l)=(/j_min,k_min,l_min/)
       end do
       end do
       end do
       temp_list=0
       do i=1,mdl%nat
          found_ion_unit=.false.
          do j=1,pub_nat_unit
             if (pub_supercell_unit_list(j)==i) then
                supercell_list(1:5,i)=(/0,0,0,j,i/)
                temp_list(j)=temp_list(j)+1
                q_calc_list(1:5,temp_list(j),j)=(/0,0,0,i,i/)
                abs_q_calc_list(1:3,temp_list(j),j)=(/0.0_DP,0.0_DP,0.0_DP/)
                found_ion_unit=.true.
                exit
             end if
          end do
          if (.not. found_ion_unit) then
             found_ion=.false.
             do j=0,pub_supercell(1)-1
             do k=0,pub_supercell(2)-1
             do l=0,pub_supercell(3)-1
                shift(1:3)=unit_cell(1:3,1)*j+&
                     unit_cell(1:3,2)*k+&
                     unit_cell(1:3,3)*l
                if (abs(j)+abs(k)+abs(l)/=0) then
                   do m=1,pub_nat_unit
                      dist=(mdl%elements(pub_supercell_unit_list(m))&
                           %centre%x-mdl%elements(i)%centre%x+shift(1))**2+&
                           (mdl%elements(pub_supercell_unit_list(m))&
                           %centre%y-mdl%elements(i)%centre%y+shift(2))**2+&
                           (mdl%elements(pub_supercell_unit_list(m))&
                           %centre%z-mdl%elements(i)%centre%z+shift(3))**2
                      if (dist<1.0d-7) then
                         if (.not. pub_phonon_SK) then
                            j_min=cutoff_translate(1,j+1,k+1,l+1)
                            k_min=cutoff_translate(2,j+1,k+1,l+1)
                            l_min=cutoff_translate(3,j+1,k+1,l+1)
                         else
                            j_min=j
                            k_min=k
                            l_min=l
                         end if
                         supercell_list(1:5,i)=(/j_min,k_min,l_min,m,&
                              pub_supercell_unit_list(m)/)
                         temp_list(m)=temp_list(m)+1
                         q_calc_list(1:5,temp_list(m),m)=(/j_min,k_min,l_min,&
                              i,pub_supercell_unit_list(m)/)
                         abs_q_calc_list(1:3,temp_list(m),m)=&
                              unit_cell(1:3,1)*j_min+&
                              unit_cell(1:3,2)*k_min+&
                              unit_cell(1:3,3)*l_min
                         found_ion=.true.
                         exit
                      end if
                   end do
                end if
                if (found_ion) exit
             end do
                if (found_ion) exit
             end do
                if (found_ion) exit
             end do
             call utils_assert(found_ion, &
                  'Incorrect supercell specifications')
          end if
       end do
       deallocate(cutoff_translate,stat=ierr)
       call utils_dealloc_check('phonon','cutoff_translate',ierr)
       deallocate(temp_list,stat=ierr)
       call utils_dealloc_check('phonon','temp_list',ierr)

    end if

    ! Calculate and store information on the number and type of displacements
    ! to be performed
    select case(pub_phonon_vib_free)
    case(0)
       phonon_vib_free(1:3)=(/.false.,.false.,.false./)
    case(1)
       phonon_vib_free(1:3)=(/.true., .false.,.false./)
    case(2)
       phonon_vib_free(1:3)=(/.false.,.true., .false./)
    case(3)
       phonon_vib_free(1:3)=(/.true., .true., .false./)
    case(4)
       phonon_vib_free(1:3)=(/.false.,.false.,.true./)
    case(5)
       phonon_vib_free(1:3)=(/.true., .false.,.true./)
    case(6)
       phonon_vib_free(1:3)=(/.false.,.true., .true./)
    case(7)
       phonon_vib_free(1:3)=(/.true., .true., .true./)
    end select
    tot_num_disp=0
    do i=1,3
       if (pub_have_supercell) then
          if (phonon_vib_free(i)) tot_num_disp=tot_num_disp+pub_nat_unit
       else
          if (phonon_vib_free(i)) tot_num_disp=tot_num_disp+mdl%nat
       end if
    end do
    if (pub_have_phonon_except_list) then
       do i=1,pub_num_except
          if (pub_phonon_iexcept_list(3,i)==1) then
             if (.not. phonon_vib_free(pub_phonon_iexcept_list(2,i))) &
                  tot_num_disp=tot_num_disp+1
          else
             if (phonon_vib_free(pub_phonon_iexcept_list(2,i))) &
                  tot_num_disp=tot_num_disp-1
          end if
       end do
    end if
    call utils_assert(tot_num_disp/=0, &
         'Number of requested force constants cannot be 0')
    if (pub_on_root) write(stdout,'(a,i4)') 'Number of force constants: ', &
         tot_num_disp
    if ((pub_phonon_farming_task==0) .or. (pub_phonon_farming_task>1)) then
       allocate(tot_disp_list(1:2,1:tot_num_disp),stat=ierr)
       call utils_alloc_check('phonon','tot_disp_list',ierr)
       allocate(tot_disp_list_back(1:3,1:mdl%nat),stat=ierr)
       call utils_alloc_check('phonon','tot_disp_list_back',ierr)
       allocate(tot_sampling_list(1:tot_num_disp),stat=ierr)
       call utils_alloc_check('phonon','tot_sampling_list',ierr)
       allocate(tot_dist_list(1:tot_num_disp),stat=ierr)
       call utils_alloc_check('phonon','tot_dist_list',ierr)
       k=0
       tot_disp_list_back=0
       if (pub_have_supercell) then
          o=pub_nat_unit
       else
          o=mdl%nat
       end if
       do n=1,o
          if (pub_have_supercell) then
             i=pub_supercell_unit_list(n)
          else
             i=n
          end if
          do j=1,3
             m=0
             do l=1,pub_num_except
                if ((pub_phonon_iexcept_list(1,l)==i) .and. &
                     (pub_phonon_iexcept_list(2,l)==j)) then
                   m=l
                   if (pub_phonon_iexcept_list(3,l)==1) then
                      do_disp=.true.
                   else
                      do_disp=.false.
                   end if
                   exit
                end if
             end do
             if (m==0) then
                if (phonon_vib_free(j)) then
                   do_disp=.true.
                else
                   do_disp=.false.
                end if
             end if
             if (do_disp) then
                k=k+1
                tot_disp_list(1:2,k)=(/i,j/)
                tot_disp_list_back(j,i)=k
                if (m==0) then
                   tot_sampling_list(k)=pub_phonon_sampling
                   tot_dist_list(k)=pub_phonon_disp
                else
                   tot_sampling_list(k)=pub_phonon_iexcept_list(4,m)
                   tot_dist_list(k)=pub_phonon_disp*pub_phonon_dexcept_list(m)
                end if
             end if
          end do
       end do
    end if

    if (pub_have_phonon_disp_list) then
       call utils_assert(&
            all(pub_phonon_disp_list(1:pub_num_disp)>=1) .and. &
            all(pub_phonon_disp_list(1:pub_num_disp)<=tot_num_disp), &
            'Invalid displacement number in phonon_disp_list')
    end if

    !==========================================================================!
    ! STAGE 1: calculate forces for unperturbed system and check that they are !
    !          close to zero                                                   !
    !==========================================================================!

    if ((pub_phonon_farming_task==0) .or. (pub_phonon_farming_task==1)) then

       if (.not. pub_write_denskern) then
          pub_write_denskern=.true.
          if (pub_on_root) write(stdout,'(a)') 'Phonon calculation: &
               &write_denskern parameter over-ridden to TRUE'
       end if
       if (.not. pub_write_tightbox_ngwfs .and. ( pub_maxit_ngwf_cg >0)) then
          pub_write_tightbox_ngwfs=.true.
          if (pub_on_root) write(stdout,'(a)') 'Phonon calculation: &
               &write_tightbox_ngwfs parameter over-ridden to TRUE'
       end if

       call energy_and_force_calculate(total_energy,forces,mdl, &
            return_converged=converged)
       E0=total_energy

       ! Check initial calculation converged
       !call utils_assert(converged, &
       !     'Calculation of initial configuration did not converge')

       ! Calculate |F|max for system
       Fmax=0.0_DP
       do i=1,mdl%nat
          Fmax=max(Fmax,dot_product(forces(1:3,i),forces(1:3,i)))
       end do
       Fmax=sqrt(Fmax)

       ! If |F|max is too large, abort
       call utils_assert(Fmax<pub_phonon_fmax, &
            'Forces for starting configuration are too large for a meaningful &
            &phonon calculation. Please perform a geometry optimization and &
            &try again.')

    end if

    !==================================================================!
    ! STAGE 2: displace atoms in turn and save force constants to file !
    !==================================================================!

    if ((pub_phonon_farming_task==0) .or. (pub_phonon_farming_task==2)) then

       if ((.not. pub_read_denskern).and.pub_geom_reuse_dk_ngwfs) then
          pub_read_denskern=.true.
          if (pub_on_root) write(stdout,'(a)') 'Phonon calculation: &
               &read_denskern parameter over-ridden to TRUE'
       end if
       if ((.not. pub_read_tightbox_ngwfs) .and. pub_geom_reuse_dk_ngwfs .and. &
           (pub_maxit_ngwf_cg>0)) then
          pub_read_tightbox_ngwfs=.true.
          if (pub_on_root) write(stdout,'(a)') 'Phonon calculation: &
               &read_tightbox_ngwfs parameter over-ridden to TRUE'
       end if
       if (pub_write_denskern) then
          pub_write_denskern=.false.
          if (pub_on_root) write(stdout,'(a)') 'Phonon calculation: &
               &write_denskern parameter over-ridden to FALSE'
       end if
       if (pub_write_tightbox_ngwfs) then
          pub_write_tightbox_ngwfs=.false.
          if (pub_on_root) write(stdout,'(a)') 'Phonon calculation: &
               &write_tightbox_ngwfs parameter over-ridden to FALSE'
       end if

       ! Allocate array for holding a row of the force constants matrix
       if (pub_on_root) then
          allocate(force_consts(1:3,1:mdl%nat),stat=ierr)
          call utils_alloc_check('phonon','force_consts',ierr)
       end if

       ! Loop over all displacements in tot_disp_list
       do m=1,tot_num_disp

          if ((pub_phonon_farming_task==0) .or. &
               (.not. pub_have_phonon_disp_list)) then
             do_disp=.true.
          else if (any(pub_phonon_disp_list(1:pub_num_disp)==m)) then
             do_disp=.true.
          else
             do_disp=.false.
          end if

          if (do_disp) then

             i=tot_disp_list(1,m)
             j=tot_disp_list(2,m)

             if (j==1) then
                cart_char='x'
             else if (j==2) then
                cart_char='y'
             else if (j==3) then
                cart_char='z'
             end if

             if (pub_on_root) then
                write(disp_number,*) (i-1)*3+j
                write(output_file,'(a)') trim(adjustl(pub_rootname))//&
                     '.force_consts_'//trim(adjustl(disp_number))
                output_unit=utils_unit()
                open(unit=output_unit,form='unformatted',&
                     file=trim(adjustl(output_file)),action='write')
             end if

                ! Loop over +ve/-ve directions
                do k=0,2*tot_sampling_list(m)-1

                   if (k<2) then
                      dist_fact=1.0_DP
                   else
                      dist_fact=2.0_DP
                   end if

                   if ((k==0) .or. (k==2)) then
                      dir_char='+'
                   else if ((k==1) .or. (k==3)) then
                      dir_char='-'
                   end if

                   if (pub_on_root) then
                      write(stdout,'()')
                      write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-&
                           &+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
                      write(stdout,'(a,i4,a,i4,a)') &
                           'Force constant calculation ', m, ' of ', &
                           tot_num_disp, ' :'
                      write(stdout,'(a,i4,a,a1,a,a1,a,f9.6,a)') &
                           '  displacing atom ', i, ' in the ', &
                           dir_char, 've ', cart_char, '-direction by ', &
                           tot_dist_list(m)*dist_fact, ' bohr'
                      write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-&
                           &+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
                   end if

                   ! displace atomic position in required direction
                   eql_pos(1:3)=(/mdl%elements(i)%centre%x,&
                                  mdl%elements(i)%centre%y,&
                                  mdl%elements(i)%centre%z/)
                   disp_pos=eql_pos
                   disp_pos(j)=eql_pos(j)+((-1.0_DP)**k)*&
                        tot_dist_list(m)*dist_fact
                   call services_rationalise_coords(1,disp_pos,mdl%cell)
                   mdl%elements(i)%centre%x=disp_pos(1)
                   mdl%elements(i)%centre%y=disp_pos(2)
                   mdl%elements(i)%centre%z=disp_pos(3)

                   ! jcap: copy this into the appropriate subregion
                   ! Get the region counter
                   iregion = mdl%elements(i)%region
                   ! Loop over atoms in this region until we find the
                   ! appropriate atom
                   sub_i=0
                   do iat = 1,mdl%regions(iregion)%par%nat
                      if(mdl%regions(iregion)%elements(iat)%global_atom_number &
                           &== mdl%elements(i)%global_atom_number) then
                         sub_i=iat
                         exit
                      end if
                   end do
                   call utils_assert(sub_i>0, &
                        'Displaced atom cannot be found in subregion')
                   ! Copy the co-ordinates over
                   mdl%regions(iregion)%elements(sub_i)%centre = &
                           &mdl%elements(i)%centre

                   ! Calculate forces for this displacement
                   call energy_and_force_calculate(total_energy,forces,mdl)
                   if ((pub_phonon_farming_task==0) .and. &
                        pub_phonon_energy_check) then
                      call utils_assert(total_energy<=E0, &
                           'Total energy of system increases upon ionic &
                           &displacement... consider increasing &
                           &phonon_finite_disp.')
                   end if

                   ! Check last calculation converged
                   !call utils_assert(converged, &
                   !     'Calculation of last configuration did not converge')

                   ! Restore equilibrium position of displaced atom
                   mdl%elements(i)%centre%x=eql_pos(1)
                   mdl%elements(i)%centre%y=eql_pos(2)
                   mdl%elements(i)%centre%z=eql_pos(3)
                   ! jcap: copy this into the appropriate subregion
                   mdl%regions(iregion)%elements(sub_i)%centre = &
                           &mdl%elements(i)%centre

                   ! Calculate force constants and write to file
                   if (pub_on_root) then
                      if (tot_sampling_list(m)==1) then
                         if (k==0) then
                            force_consts=forces
                         else
                            force_consts=-0.5_DP*(force_consts-forces)/&
                                 tot_dist_list(m)
                            write(output_unit) force_consts
                         end if
                      else if (tot_sampling_list(m)==2) then
                         if (k==0) then
                            force_consts=forces
                         else if (k==1) then
                            force_consts=force_consts-forces
                         else if (k==2) then
                            force_consts=8.0_DP*force_consts-forces
                         else if (k==3) then
                            force_consts=-(force_consts+forces)/&
                                 (12.0_DP*tot_dist_list(m))
                            write(output_unit) force_consts
                         end if
                      end if
                   end if

                end do

                if (pub_on_root) close(output_unit)

          end if

       end do

       if (pub_on_root) then
          deallocate(force_consts,stat=ierr)
          call utils_dealloc_check('phonon','force_consts',ierr)
       end if

    end if

    !=================================================================!
    ! STAGE 3: construct dynamical matrix and find phonon frequencies !
    !=================================================================!

    if ((pub_phonon_farming_task==0) .or. (pub_phonon_farming_task==3)) then

       allocate(freqs(1:tot_num_disp),stat=ierr)
       call utils_alloc_check('phonon','freqs',ierr)
       allocate(force_consts_m(1:3,1:mdl%nat,1:tot_num_disp),stat=ierr)
       call utils_alloc_check('phonon','force_consts_m',ierr)

       ! allocate full dynamical matrix
       call dense_create(dynamical_symm_dem,tot_num_disp,tot_num_disp,.true.)
       call dense_create(dynamical_dem,tot_num_disp,tot_num_disp,.true.)
       call dense_create(eigenvecs_dem,tot_num_disp,tot_num_disp,.true.)
       call dense_create(unit_mat_dem,tot_num_disp,tot_num_disp,.true.)

       ! initialize unit matrix
       do i=1,tot_num_disp
          do j=1,tot_num_disp
             if (i==j) then
                call dense_put_element(cmplx_1,unit_mat_dem,j,i)
             else
                call dense_put_element(cmplx_0,unit_mat_dem,j,i)
             end if
          end do
       end do

       ! read force constants back in from file
       if (pub_on_root) then
          do m=1,tot_num_disp
             i=tot_disp_list(1,m)
             j=tot_disp_list(2,m)
             write(disp_number,*) (i-1)*3+j
             write(output_file,'(a)') trim(adjustl(pub_rootname))//&
                  '.force_consts_'//trim(adjustl(disp_number))
             output_unit=utils_unit()
             open(unit=output_unit,form='unformatted',&
                 file=trim(adjustl(output_file)),action='read',&
                  status='old')
             read(output_unit) force_consts_m(1:3,1:mdl%nat,m)
             close(output_unit)
          end do
       end if
       call comms_bcast(pub_root_proc_id,force_consts_m)

       ! loop over all q-points in the regular grid
       zero_point_E=0.0_DP
       T_num_iter=nint((pub_phonon_tmax-pub_phonon_tmin)/pub_phonon_deltat)
       allocate(F(1:T_num_iter),stat=ierr)
       call utils_alloc_check('phonon','F',ierr)
       allocate(S(1:T_num_iter),stat=ierr)
       call utils_alloc_check('phonon','S',ierr)
       allocate(C_v(1:T_num_iter),stat=ierr)
       call utils_alloc_check('phonon','C_v',ierr)
       F=0.0_DP
       S=0.0_DP
       C_v=0.0_DP
       print_warning=.true.
       print_warning2=.true.
       grid_loop=0
       tot_grid_loop=pub_phonon_grid(1)*pub_phonon_grid(2)*pub_phonon_grid(3)
       if (pub_phonon_DOS) then
         nqdos=int((pub_phonon_DOS_max-pub_phonon_DOS_min) / &
              pub_phonon_DOS_delta)
         allocate(qdos(0:nqdos),stat=ierr)
         call utils_alloc_check('phonon','qdos',ierr)
         qdos=0.0_DP
         qdos_weight=1.0_DP/(tot_grid_loop*pub_phonon_DOS_delta)
       end if
       do g1=0,pub_phonon_grid(1)-1
       do g2=0,pub_phonon_grid(2)-1
       do g3=0,pub_phonon_grid(3)-1

          grid_loop=grid_loop+1

          ! calculate q-point in reciprocal space
          q(1:3)=rec_unit_cell(1:3,1)*real(g1,dp)/real(pub_phonon_grid(1),dp)+&
               rec_unit_cell(1:3,2)*real(g2,dp)/real(pub_phonon_grid(2),dp)+&
               rec_unit_cell(1:3,3)*real(g3,dp)/real(pub_phonon_grid(3),dp)

          ! calculate dynamical matrix for q
          do m=1,tot_num_disp
          do n=1,tot_num_disp
             k=tot_disp_list(1,n)
             l=tot_disp_list(2,n)
             if (pub_have_supercell) then
                el=cmplx_0
                do o=1,ncells
                   p=supercell_list(4,k)
                   el=el+force_consts_m(l,q_calc_list(4,o,p),m)*&
                        exp(-cmplx_i*(q(1)*abs_q_calc_list(1,o,p)+&
                        q(2)*abs_q_calc_list(2,o,p)+&
                        q(3)*abs_q_calc_list(3,o,p)))
                end do
             else
                el=cmplx(force_consts_m(l,k,m), 0.0_DP, DP)
             end if
             mass_m = periodic_table_mass(mdl%elements( &
                  tot_disp_list(1,m))%atomic_number) / electron_mass_u
             mass_n = periodic_table_mass(mdl%elements( &
                  tot_disp_list(1,n))%atomic_number) / electron_mass_u
             el=el/sqrt(mass_m*mass_n)
             call dense_put_element(el,dynamical_dem,n,m)
          end do
          end do

          ! make dynamical matrix Hermitian
          do i=1,tot_num_disp
             call dense_get_element(el,dynamical_dem,i,i)
             call dense_put_element(el,dynamical_symm_dem,i,i)
          end do

          do i=1,tot_num_disp
             do j=i+1,tot_num_disp
                call dense_get_element(el,dynamical_dem,i,j)
                call dense_get_element(el_T,dynamical_dem,j,i)
                el_symm=(el+conjg(el_T))*0.5_DP
                call dense_put_element(el_symm,dynamical_symm_dem,i,j)
                call dense_put_element(conjg(el_symm),dynamical_symm_dem,j,i)
             end do
          end do

          ! diagonalize dynamical matrix to find phonon frequencies
          call dense_eigensolve(tot_num_disp,freqs,dynamical_symm_dem,&
               unit_mat_dem,1,eigenvecs_dem,pub_eigensolver_orfac,&
               pub_eigensolver_abstol)

          ! print out phonon frequencies at Gamma
          if ((pub_on_root) .and. (grid_loop==1)) then
             write(stdout,'()')
             write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
                  &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
             write(stdout,'(a)') 'Gamma phonon frequencies         (cm^-1)     &
                  &            (THz)'
             write(stdout,'()')
             do i=1,tot_num_disp
                freq_conv=sign(sqrt(abs(freqs(i)))/TWO_PI,freqs(i))
                write(stdout,'(12x,i5,a,2(6x,f16.11))') i, ":", &
                     freq_conv*au2inv_cm, freq_conv*au2THz
             end do
             write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
                  &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

             ! ndmh: print QC test data
             if (pub_print_qc) then
                i=1
                freq_conv=sign(sqrt(abs(freqs(i)))/TWO_PI,freqs(i))
                call utils_qc_print('phonon_freq_1',freq_conv)
                i=tot_num_disp
                freq_conv=sign(sqrt(abs(freqs(i)))/TWO_PI,freqs(i))
                call utils_qc_print('phonon_freq_3N',freq_conv)
             end if
          end if

          ! calculate and print out zero-point energy
          if (pub_on_root) then
             do i=1,tot_num_disp
                if (freqs(i)<0.0_DP) then
                   if (print_warning) then
                      print_warning=.false.
                      write(stdout,'(a)') 'WARNING: Imaginary phonon &
                           &frequencies excluded from the computation of the &
                           &zero-point energy'
                   end if
                else
                   zero_point_E=zero_point_E+sqrt(freqs(i))
                end if
             end do
             if (grid_loop==tot_grid_loop) then
                zero_point_E=0.5_DP*zero_point_E/tot_grid_loop
                write(stdout,'()')
                write(stdout,'(a,1x,f16.11,1x,a)') "Zero-point energy =", &
                     zero_point_E, "Eh"

                ! ndmh: print QC test data
                if (pub_print_qc) call utils_qc_print('zero_point_energy',&
                     zero_point_E)
             end if
          end if

          ! calculate and print out F, S, U, and C_v for specified range of T
          if (pub_on_root) then
             T=pub_phonon_tmin
             if (grid_loop==tot_grid_loop) then
                write(stdout,'()')
                write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-&
                     &+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
                     &-+-+-'
                write(stdout,'(a)') '                             Helmholtz    &
                     &                                Internal              &
                     &Specific'
                write(stdout,'(a)') ' Temperature (K)      free energy (Eh)    &
                     &    Entropy (Eh/K)           energy (Eh)           heat &
                     &(Eh/K)'
                write(stdout,'()')
             end if
             do i=0,T_num_iter
                if ((T==0.0_DP) .and. (grid_loop==tot_grid_loop)) then
                   write(stdout,'(f16.11,4(6x,f16.11))') 0.0_DP, zero_point_E, &
                        0.0_DP, zero_point_E, 0.0_DP
                else
                   beta=1.0_DP/T
                   do j=1,tot_num_disp
                      sqrt_freq=sign(sqrt(abs(freqs(j))),freqs(j))
                      if (sqrt_freq/TWO_PI < pub_phonon_min_freq) then
                         if (print_warning2) then
                            print_warning2=.false.
                            write(stdout,'(a)') 'WARNING: Low phonon &
                                 &frequencie excluded from the computation of &
                                 &thermodynamic quantities'
                         end if
                      else
                         exp_freq=exp(-beta*sqrt_freq)
                         F(i)=F(i)+log(1.0_DP-exp_freq)
                         S(i)=S(i)+sqrt_freq*exp_freq/(1.0_DP-exp_freq)
                         C_v(i)=C_v(i)+sqrt_freq**2*exp_freq/&
                              ((1.0_DP-exp_freq)**2)
                      end if
                   end do
                   if (grid_loop==tot_grid_loop) then
                      F(i)=T*F(i)/tot_grid_loop
                      S(i)=S(i)/tot_grid_loop-F(i)
                      F(i)=F(i)+zero_point_E
                      U=F(i)+S(i)
                      S(i)=k_B*S(i)/T
                      C_v(i)=k_B*C_v(i)/(tot_grid_loop*T**2)
                      write(stdout,'(f16.11,4(6x,f16.11))') T/k_B, F(i), S(i), &
                           U, C_v(i)
                   end if
                end if
                T=T+pub_phonon_deltat
             end do
             if (grid_loop==tot_grid_loop) then
                if (pub_on_root) write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+&
                  &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
                  &-+-+-+-+-+-+-+-+-+-'

                ! ndmh: print QC test data
                if (pub_print_qc) then
                   T=T-pub_phonon_deltat
                   call utils_qc_print('temperature',T/k_B)
                   call utils_qc_print('helmholtz_f_e',F(T_num_iter))
                   call utils_qc_print('entropy',S(T_num_iter))
                   call utils_qc_print('internal_e',U)
                   call utils_qc_print('specific_heat',C_v(T_num_iter))
                end if

                deallocate(C_v,stat=ierr)
                call utils_dealloc_check('phonon','C_v',ierr)
                deallocate(S,stat=ierr)
                call utils_dealloc_check('phonon','S',ierr)
                deallocate(F,stat=ierr)
                call utils_dealloc_check('phonon','F',ierr)
             end if
          end if

          ! calculate phonon DOS and write to file
          if (pub_phonon_DOS .and. pub_on_root) then
             do i=1,tot_num_disp
                j=nint((sign(sqrt(abs(freqs(i)))/TWO_PI,freqs(i))*au2inv_cm-&
                     pub_phonon_DOS_min)/pub_phonon_DOS_delta)
                if ((j>=0) .and. (j<=nqdos)) qdos(j)=qdos(j)+qdos_weight
             end do
             if (grid_loop==tot_grid_loop) then
                write(output_file,'(a)') trim(adjustl(pub_rootname))//&
                     '.qdos'
                output_unit=utils_unit()
                open(unit=output_unit,form='formatted',&
                     file=trim(adjustl(output_file)),action='write')
                write(output_unit,'(a)') '#           cm^-1              THz   &
                     &    Phonon DOS'
                do i=0,nqdos
                   dist=pub_phonon_DOS_min+i*pub_phonon_DOS_delta
                   write(output_unit,'(3(1x,f16.11))') dist, &
                        dist*au2THz/au2inv_cm, qdos(i)
                end do
                close(output_unit)
                deallocate(qdos,stat=ierr)
                call utils_dealloc_check('phonon','qdos',ierr)
             end if
          end if

       end do
       end do
       end do

       if (pub_have_phonon_qpoints) then

          ! loop over q-points in user list
          if (pub_on_root) then
             write(output_file,'(a)') trim(adjustl(pub_rootname))//&
                  '.phonon_freqs'
             output_unit=utils_unit()
             open(unit=output_unit,form='formatted',&
                  file=trim(adjustl(output_file)),action='write')
          end if
          if (pub_phonon_write_eigenvecs) then
             allocate(ion_write(1:3,1:mdl%nat),stat=ierr)
             call utils_alloc_check('phonon','ion_write',ierr)
          end if
          do g1=1,pub_num_qpoints

             q(1:3)=rec_unit_cell(1:3,1)*pub_phonon_qpoints(1,g1)+&
                  rec_unit_cell(1:3,2)*pub_phonon_qpoints(2,g1)+&
                  rec_unit_cell(1:3,3)*pub_phonon_qpoints(3,g1)

             if (pub_on_root) then
                write(output_unit,'(a,3(1x,f16.11),a)') '#', q(1:3), &
                     '    <-- QPT'
             end if

             ! work out dynamical matrix for q
             do m=1,tot_num_disp
             do n=1,tot_num_disp
                k=tot_disp_list(1,n)
                l=tot_disp_list(2,n)
                if (pub_have_supercell) then
                   el=cmplx_0
                   do o=1,ncells
                      p=supercell_list(4,k)
                      el=el+force_consts_m(l,q_calc_list(4,o,p),m)*&
                           exp(-cmplx_i*(q(1)*abs_q_calc_list(1,o,p)+&
                           q(2)*abs_q_calc_list(2,o,p)+&
                           q(3)*abs_q_calc_list(3,o,p)))
                   end do
                else
                   el=cmplx(force_consts_m(l,k,m), 0.0_DP, DP)
                end if
                mass_m = periodic_table_mass(mdl%elements( &
                     tot_disp_list(1,m))%atomic_number) / electron_mass_u
                mass_n = periodic_table_mass(mdl%elements( &
                     tot_disp_list(1,n))%atomic_number) / electron_mass_u
                el=el/sqrt(mass_m*mass_n)
                call dense_put_element(el,dynamical_dem,n,m)
             end do
             end do

             ! make dynamical matrix Hermitian
             do i=1,tot_num_disp
                call dense_get_element(el,dynamical_dem,i,i)
                call dense_put_element(el,dynamical_symm_dem,i,i)
             end do

             do i=1,tot_num_disp
                do j=i+1,tot_num_disp
                   call dense_get_element(el,dynamical_dem,i,j)
                   call dense_get_element(el_T,dynamical_dem,j,i)
                   el_symm=(el+conjg(el_T))*0.5_DP
                   call dense_put_element(el_symm,dynamical_symm_dem,i,j)
                   call dense_put_element(conjg(el_symm),dynamical_symm_dem,j,i)
                end do
             end do

             ! diagonalize dynamical matrix to find phonon frequencies
             call dense_eigensolve(tot_num_disp,freqs,dynamical_symm_dem,&
                  unit_mat_dem,1,eigenvecs_dem,pub_eigensolver_orfac,&
                  pub_eigensolver_abstol)

             ! write phonon frequencies and eigenvectors to file
             do i=1,tot_num_disp
                if (pub_on_root) &
                     write(output_unit,'(i5,1x,f16.11,a)') i, &
                     sign(sqrt(abs(freqs(i)))/TWO_PI,freqs(i))*au2inv_cm, &
                     '                                  <-- FRQ'
                if (pub_phonon_write_eigenvecs) then
                   ion_write=cmplx_0
                   do j=1,mdl%nat
                      do k=1,3
                         if (tot_disp_list_back(k,j)/=0) &
                              call dense_get_element(ion_write(k,j),&
                              eigenvecs_dem,tot_disp_list_back(k,j),i)
                      end do
                      if (pub_on_root .and. pub_phonon_write_eigenvecs) &
                           write(output_unit,'(a4,3(1x,a,2(f16.11,a)),a)') &
                           mdl%elements(j)%species_id, &
                                   '(', real(ion_write(1,j),DP), ',', &
                           aimag(ion_write(1,j)), ')', &
                           '(', real(ion_write(2,j),DP), ',', &
                           aimag(ion_write(2,j)), ')', &
                           '(', real(ion_write(3,j),DP), ',', &
                           aimag(ion_write(3,j)), ')', &
                           ' <-- VEC'
                   end do
                end if
             end do
             if (pub_phonon_write_eigenvecs .and. (g1==pub_num_qpoints)) then
                deallocate(ion_write,stat=ierr)
                call utils_dealloc_check('phonon','ion_write',ierr)
             end if

          end do

          if (pub_on_root) close(output_unit)

       end if

       call dense_destroy(unit_mat_dem)
       call dense_destroy(eigenvecs_dem)
       call dense_destroy(dynamical_dem)
       call dense_destroy(dynamical_symm_dem)

       if (pub_have_phonon_animate_list) then

          ! allocate full dynamical matrix
          call dense_create(dynamical_symm_dem,tot_num_disp,tot_num_disp,&
               .false.)
          call dense_create(dynamical_dem,tot_num_disp,tot_num_disp,.false.)
          call dense_create(eigenvecs_dem,tot_num_disp,tot_num_disp,.false.)
          call dense_create(unit_mat_dem,tot_num_disp,tot_num_disp,.false.)

          ! initialize unit matrix
          do i=1,tot_num_disp
             do j=1,tot_num_disp
                if (i==j) then
                   call dense_put_element(1.0_DP,unit_mat_dem,j,i)
                else
                   call dense_put_element(0.0_DP,unit_mat_dem,j,i)
                end if
             end do
          end do

          allocate(ion_write2(1:3,1:mdl%nat),stat=ierr)
          call utils_alloc_check('phonon','ion_write2',ierr)

          ! work out dynamical matrix at Gamma
          do m=1,tot_num_disp
          do n=1,tot_num_disp
             k=tot_disp_list(1,n)
             l=tot_disp_list(2,n)
             if (pub_have_supercell) then
                el2=0.0_DP
                do o=1,ncells
                   p=supercell_list(4,k)
                   el2=el2+force_consts_m(l,q_calc_list(4,o,p),m)
                end do
             else
                el2=force_consts_m(l,k,m)
             end if
             mass_m = periodic_table_mass(mdl%elements( &
                  tot_disp_list(1,m))%atomic_number) / electron_mass_u
             mass_n = periodic_table_mass(mdl%elements( &
                  tot_disp_list(1,n))%atomic_number) / electron_mass_u
             el2 = el2/sqrt(mass_m*mass_n)
             call dense_put_element(el2,dynamical_dem,n,m)
          end do
          end do

          ! symmetrize dynamical matrix
          do i=1,tot_num_disp
             call dense_get_element(el2,dynamical_dem,i,i)
             call dense_put_element(el2,dynamical_symm_dem,i,i)
          end do

          do i=1,tot_num_disp
             do j=i+1,tot_num_disp
                call dense_get_element(el2,dynamical_dem,i,j)
                call dense_get_element(el2_T,dynamical_dem,j,i)
                el2_symm=(el2+el2_T)*0.5_DP
                call dense_put_element(el2_symm,dynamical_symm_dem,i,j)
                call dense_put_element(el2_symm,dynamical_symm_dem,j,i)
             end do
          end do

          ! diagonalize dynamical matrix to find phonon frequencies
          call dense_eigensolve(tot_num_disp,freqs,dynamical_symm_dem,&
               unit_mat_dem,1,eigenvecs_dem,pub_eigensolver_orfac,&
               pub_eigensolver_abstol)

          ! write animation files for selected phonon frequencies
          do i=1,tot_num_disp
             if (any(pub_phonon_animate_list(1:pub_num_ani)==i)) then
                ion_write2=0.0_DP
                do j=1,mdl%nat
                   do k=1,3
                      if (tot_disp_list_back(k,j)/=0) &
                           call dense_get_element(ion_write2(k,j),&
                           eigenvecs_dem,tot_disp_list_back(k,j),i)
                   end do
                end do
                if (pub_on_root) then
                   write(disp_number,*) i
                   write(output_file,'(a)') trim(adjustl(pub_rootname))//&
                        '.phonon_'//trim(adjustl(disp_number))//'.xyz'
                   output_unit=utils_unit()
                   open(unit=output_unit,form='formatted',&
                        file=trim(adjustl(output_file)),action='write')
                   do k=0,15
                      el2=sin(k*TWO_PI*0.0625_DP)*pub_phonon_animate_scale
                      write(output_unit,'(i7)') mdl%nat
                      write(output_unit,'(f16.11,1x,a)') sign(sqrt(abs(freqs(i)))/&
                           TWO_PI,freqs(i))*au2inv_cm, 'cm^-1'
                      if (pub_have_supercell) then
                         do j=1,mdl%nat
                            write(output_unit,'(a4,3(1x,f16.11))') &
                                 mdl%elements(j)%species_id, &
                                 mdl%elements(j)%centre%x/ANGSTROM+&
                                 el2*ion_write2(1,supercell_list(5,j)), &
                                 mdl%elements(j)%centre%y/ANGSTROM+&
                                 el2*ion_write2(2,supercell_list(5,j)), &
                                 mdl%elements(j)%centre%z/ANGSTROM+&
                                 el2*ion_write2(3,supercell_list(5,j))
                         end do
                      else
                         do j=1,mdl%nat
                            write(output_unit,'(a4,3(1x,f16.11))') &
                                 mdl%elements(j)%species_id, &
                                 mdl%elements(j)%centre%x/ANGSTROM+&
                                 el2*ion_write2(1,j), &
                                 mdl%elements(j)%centre%y/ANGSTROM+&
                                 el2*ion_write2(2,j), &
                                 mdl%elements(j)%centre%z/ANGSTROM+&
                                 el2*ion_write2(3,j)
                         end do
                      end if
                   end do
                   close(output_unit)
                end if
             end if
          end do
          deallocate(ion_write2,stat=ierr)
          call utils_dealloc_check('phonon','ion_write2',ierr)

       end if

       deallocate(force_consts_m,stat=ierr)
       call utils_dealloc_check('phonon','force_consts_m',ierr)
       deallocate(freqs,stat=ierr)
       call utils_dealloc_check('phonon','freqs',ierr)

    end if

    if (pub_have_supercell) then
       deallocate(abs_q_calc_list,stat=ierr)
       call utils_dealloc_check('phonon','abs_q_calc_list',ierr)
       deallocate(q_calc_list,stat=ierr)
       call utils_dealloc_check('phonon','q_calc_list',ierr)
       deallocate(supercell_list,stat=ierr)
       call utils_dealloc_check('phonon','supercell_list',ierr)
    end if

    if ((pub_phonon_farming_task==0) .or. (pub_phonon_farming_task>1)) then
       deallocate(tot_dist_list,stat=ierr)
       call utils_dealloc_check('phonon','tot_dist_list',ierr)
       deallocate(tot_sampling_list,stat=ierr)
       call utils_dealloc_check('phonon','tot_sampling_list',ierr)
       deallocate(tot_disp_list_back,stat=ierr)
       call utils_dealloc_check('phonon','tot_disp_list_back',ierr)
       deallocate(tot_disp_list,stat=ierr)
       call utils_dealloc_check('phonon','tot_disp_list',ierr)
    end if

    if (allocated(pub_phonon_qpoints)) then
       deallocate(pub_phonon_qpoints,stat=ierr)
       call utils_dealloc_check('phonon_main','pub_phonon_qpoints',ierr)
    end if
    if (allocated(pub_supercell_unit_list)) then
       deallocate(pub_supercell_unit_list,stat=ierr)
       call utils_dealloc_check('phonon_main','pub_supercell_unit_list',ierr)
    end if
    if (allocated(pub_phonon_animate_list)) then
       deallocate(pub_phonon_animate_list,stat=ierr)
       call utils_dealloc_check('phonon_main','pub_phonon_animate_list',ierr)
    end if
    if (allocated(pub_phonon_dexcept_list)) then
       deallocate(pub_phonon_dexcept_list,stat=ierr)
       call utils_dealloc_check('phonon_main','pub_phonon_dexcept_list',ierr)
    end if
    if (allocated(pub_phonon_iexcept_list)) then
       deallocate(pub_phonon_iexcept_list,stat=ierr)
       call utils_dealloc_check('phonon_main','pub_phonon_iexcept_list',ierr)
    end if
    if (allocated(pub_phonon_disp_list)) then
       deallocate(pub_phonon_disp_list,stat=ierr)
       call utils_dealloc_check('phonon_main','pub_phonon_disp_list',ierr)
    end if

  end subroutine phonon_main

end module phonon
