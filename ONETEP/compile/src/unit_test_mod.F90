! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                      U N I T    T E S T I N G                               !
!=============================================================================!
!                                                                             !
! This module executes unit testing of the core ONETEP modules.               !
!                                                                             !
!-----------------------------------------------------------------------------!
! Originally written by Max Phipps in 2015.                                   !
!-----------------------------------------------------------------------------!

module unit_test

  implicit none

  private

  integer :: num_t_failed

  public :: unit_test_exec
  public :: unit_test_init_dummy

contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine unit_test_exec(mdl)

    !==========================================================!
    ! This subroutine is the main driver for the unit testing. !
    ! This subroutine calls the individual unit testing        !
    ! routines and handles the unit testing summary printouts. !
    !----------------------------------------------------------!
    ! Written by Max Phipps on 23/07/2015.                     !
    !==========================================================!

    use constants, only: stdout
    use comms, only: pub_on_root
    use model_type, only: MODEL
    use services, only: services_flush
    use utils, only: utils_banner

    implicit none

    ! Arguments
    type(MODEL), target, intent(inout)   :: mdl

    ! mjsp: Initialisations
    num_t_failed = 0

    ! mjsp: Print title
    if (pub_on_root) then
       write(stdout,'(/a)') repeat('#', 80)
       write(stdout,'(a)') utils_banner('#', ' U N I T    T E S T I N G ')
       write(stdout,'(a/)') repeat('#', 80)
    end if


    ! mjsp: Sparse matrix module testing
    if (pub_on_root) write(stdout,'(/4x,a)') "==>==>==>==>==>==>==>==>==>==>&
        &==>==>==>==>==>==>==>==>==>==>==>==>==>==>==>"
    if (pub_on_root) &
       write(stdout,'(4x,a)') "==> UNIT TEST: Sparse Matrix (SPAM3) Module"
    if (pub_on_root) write(stdout,'(4x,a/)') "==>==>==>==>==>==>==>==>==>==>&
        &==>==>==>==>==>==>==>==>==>==>==>==>==>==>==>"
    call services_flush()

    ! TODO: sparse module unit testing.
    if (pub_on_root) write(stdout,'(4x,a/)') "(NOT YET IMPLEMENTED)"



    ! mjsp: Fourier transform module testing
    if (pub_on_root) write(stdout,'(/4x,a)') repeat('==>', 25)
    if (pub_on_root) &
       write(stdout,'(4x,a)') "==> UNIT TEST: Fast Fourier Transform Module"
    if (pub_on_root) write(stdout,'(/4x,a)') repeat('==>', 25)

    call unit_test_fourier(mdl)


    ! mjsp: Summary printout
    if (pub_on_root) write(stdout,'(//a50)') '===== SUMMARY ===== '
    if (num_t_failed .gt. 0) then
      if (pub_on_root) &
         write(stdout,'(/a28,i3,a/)') '!!! ',num_t_failed,&
                                      ' UNIT TEST(S) FAILED !!!'
    else
      if (pub_on_root) write(stdout,'(/28x,a/)') ' ALL UNIT TESTS PASSED '
    end if


    ! mjsp: Print exit title
    if (pub_on_root) then
       write(stdout,'(/a)') repeat('#', 80)
       write(stdout,'(a)') utils_banner('#', &
            ' E N D    U N I T    T E S T I N G ')
       write(stdout,'(a/)') repeat('#', 80)
    end if

  end subroutine unit_test_exec

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine unit_test_fourier(mdl)

    !==========================================================!
    ! This subroutine handles unit testing of the FFT module.  !
    !----------------------------------------------------------!
    ! Written by Max Phipps on 23/07/2015.                     !
    !==========================================================!

    use constants, only: stdout, DP, PI
    use model_type, only: MODEL
    use rundat, only: pub_fine_grid_scale, pub_fine_is_dbl
    use fourier, only: fourier_init_cell, fourier_init_threads, &
        fourier_apply_cell_backward, fourier_apply_cell_forward, &
        fourier_exit_cell, fourier_exit_threads
    use comms, only: pub_on_root
    use cell_grid, only: cell_grid_distribute
    use services, only: services_flush
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(MODEL), target, intent(inout)   :: mdl

    ! Local Variables
    complex(kind=DP), allocatable :: recip(:,:,:)
    real(kind=DP), allocatable    :: dens_real(:,:,:)
    real(kind=DP), allocatable    :: dens_real_test(:,:,:)
    integer       :: ierr
    real(kind=DP) :: test_double, freq_mod
    integer       :: i1, i2, i3
    real          :: fft_thresh = 1E-09_DP

    if (pub_on_root) write(stdout,'(a)',advance='no') "==> Initialising FFT &
       &threading and cell..."

    ! mjsp: initialise FFT threading
    call fourier_init_threads

    ! mjsp: initialise standard grid and standard grid FFTs
    call cell_grid_distribute(mdl%fine_grid,mdl%cell, &
      pub_fine_grid_scale,'  Fine grid',.true.,.true.)
    call fourier_init_cell(mdl%fine_grid)
    pub_fine_is_dbl = .false.

    if (pub_on_root) write(stdout,'(a)') '... done'


    ! Allocate
    allocate(dens_real(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
      mdl%fine_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('unit_test_exec','dens_real',ierr)
    allocate(dens_real_test(mdl%fine_grid%ld1,mdl%fine_grid%ld2,&
      mdl%fine_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('unit_test_exec','dens_real',ierr)
    allocate(recip(mdl%fine_grid%ld3,mdl%fine_grid%ld2,&
      mdl%fine_grid%max_slabs23),stat=ierr)
    call utils_alloc_check('unit_test_exec','recip',ierr)



    ! mjsp: Zero Fourier Test:
    if (pub_on_root) &
      write(stdout,'(/a30)',advance='no') "==> Zeroes on grid ..."
    call services_flush()

    ! mjsp: Dummy density on grid
    dens_real(:,:,:) = 0.0_dp
    recip(:,:,:) = 0.0_dp
    dens_real_test(:,:,:) = 0.0_dp

    ! Fourier transform the charge density to reciprocal space
    call fourier_apply_cell_forward(dens_real(:,:,:),recip,mdl%fine_grid)

    ! FFT the local ionic potential from reciprocal to real space
    call fourier_apply_cell_backward(dens_real_test,recip,mdl%fine_grid)

    ! mjsp: Test that the difference of the back-transformed test input to the
    ! test input is within threshold.
    test_double = sum(abs(dens_real-dens_real_test))
    if (test_double .gt. fft_thresh) then
      if (pub_on_root) write(stdout,'(a,ES15.9/)') '... FAILED! Test value = ',&
                        test_double
      num_t_failed = num_t_failed + 1
    else
      if (pub_on_root) write(stdout,'(a)') '... PASSED'
    end if


    ! mjsp: Sine Fourier Test
    if (pub_on_root) write(stdout,'(a30)',advance='no') "==> Sine on grid ..."
    call services_flush()

    ! mjsp: Dummy density on grid
    dens_real(:,:,:) = 0.0_dp
    recip(:,:,:) = 0.0_dp
    dens_real_test(:,:,:) = 0.0_dp

    ! mjsp: Sine on grid
    freq_mod = 1.0_DP * PI / real(mdl%fine_grid%ld1 + mdl%fine_grid%ld2,kind=DP)
    do i1=1,mdl%fine_grid%ld1
      do i2=1,mdl%fine_grid%ld2
        do i3=1,mdl%fine_grid%max_slabs12
          dens_real(i1,i2,i3) = sin(real(i1-0.5_dp*i2,kind=DP)*freq_mod) &
                                - 0.5_DP*PI
        end do
      end do
    end do

    ! mjsp: Fourier transform the charge density to reciprocal space
    call fourier_apply_cell_forward(dens_real(:,:,:),recip,mdl%fine_grid)

    ! mjsp: FFT from reciprocal to real space
    call fourier_apply_cell_backward(dens_real_test,recip,mdl%fine_grid)

    ! mjsp: Test that the difference of the back-transformed test input to the
    ! test input is within threshold.
    test_double = sum(abs(dens_real(1:mdl%fine_grid%ld3,:,:) - &
                     dens_real_test(1:mdl%fine_grid%ld3,:,:)))
    if (test_double .gt. fft_thresh) then
      if (pub_on_root) &
         write(stdout,'(a,ES15.9/)') '... FAILED! Test value = ', test_double
      num_t_failed = num_t_failed + 1
    else
      if (pub_on_root) write(stdout,'(a)') '... PASSED'
    end if

    ! Cleanup
    call fourier_exit_cell()
    call fourier_exit_threads()

    ! Deallocate
    deallocate(recip,stat=ierr)
    call utils_dealloc_check('unit_test_exec','recip',ierr)
    deallocate(dens_real,stat=ierr)
    call utils_dealloc_check('unit_test_exec','dens_real',ierr)
    deallocate(dens_real_test,stat=ierr)
    call utils_dealloc_check('unit_test_exec','dens_real_test',ierr)

  end subroutine unit_test_fourier

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine unit_test_init_dummy(mdl)

    !==========================================================!
    ! (PLACEHOLDER) This subroutine handles initialisation of  !
    ! the unit testing dummy cell.                             !
    !----------------------------------------------------------!
    ! Written by Max Phipps on 23/07/2015.                     !
    !==========================================================!

    use constants, only: stdout, DP
    use model_type, only: MODEL
    !use simulation_cell, only: simulation_cell_initialise

    implicit none

    ! Arguments
    type(MODEL), target, intent(inout)   :: mdl

    ! mjsp: Initialise mdl%cell
    !call simulation_cell_initialise(mdl%cell,a1,a2,a3, &
    !     n_pt1,n_pt2,n_pt3,d1,d2,d3)

    !write(stdout,*) mdl%cell%a1,mdl%cell%a2,mdl%cell%a3
    !stop

  end subroutine unit_test_init_dummy

end module unit_test




