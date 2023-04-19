! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!==============================================================================!
!                                                                              !
!  ONETEP non-local exchange correlation functional module: nlxc_mod.F90       !
!                                                                              !
!  The subroutines in this file were written by Lampros Andrinopoulos          !
!  with help from Nicholas Hine, in October 2011 to November 2012.             !
!                                                                              !
!  VV10 (and VV10sol) non-local correlation-energy functionality was written by!
!  Gabriel Constantinescu in July 2014                                         !
!                                                                              !
!  The routines in this module implement the approach known as the van der     !
!  Waals Density Functional (vdW-DF) first described in                        !
!    M. Dion, H. Rydberg, E. Schroder, D.C. Langreth, and B.I. Lundqvist       !
!    Phys. Rev. Lett. 92, 246401 (2004)                                        !
!  and                                                                         !
!    T. Thonhauser, V.R. Cooper, S. Li, A. Puzder, P. Hyldgaard, and D.C.      !
!    Langreth, Phys. Rev. B 76, 125112 (2007).                                 !
!  The implementation follows the approach presented by Roman-Perez and Soler: !
!    G. Roman-Perez and J. M. Soler, Phys. Rev. Lett. 103, 096102 (2009)       !
!                                                                              !
!  For the VV10 and VV10sol parts please consult the works of Sabatini et al.  !
!  (PHYSICAL REVIEW B 87, 041108(R) 2013) and T. Bjorkman (Phys. Rev. B 86,    !
!  165109 2012), respectively                                                  !
!                                                                              !
!  Spin Polarization support in VDWDF1+2 by Jolyon Aarons, 2021.               !
!  Based on T. Thonhauser, S. Zuluaga, C. A. Arter, K. Berland, E. Schroeder   !
!  and P. Hyldgaard, Phys. Rev. Lett. 115, 136402 (2015)                       !
!  VV10 support based various statements in papers citing the VV10 paper +     !
!  the VV09 paper + own derivation.                                            !
!==============================================================================!

module nlxc

  use cell_grid, only: GRID_INFO
  use comms, only: pub_on_root
  use constants, only: stdout, DP, PI, VERBOSE
  use rundat, only: pub_output_detail, pub_num_spins
  use timer, only: timer_clock
  use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
       utils_open_unit_check, utils_close_unit_check
  use xc_funcs

  implicit none

  private

  public :: nlxc_vdw_energy
  real(kind=DP), parameter :: THIRD = 1.0_DP / 3.0_DP
  real(kind=DP), parameter :: FTHRD = 4.0_DP / 3.0_DP
  integer, parameter :: m_c = 12

contains

  subroutine nlxc_vdw_energy(nlc_energy_real,nl_pot_fine,density_fine, &
       density_grad,grid,cell,vdwdf_type,vv10_bval,vv10_Cval)

    !======================================================================!
    ! This subroutine calculates the non-local correlation energy E_c^nl   !
    ! according to Soler Eq. 8.                                            !
    !                                                                      !
    ! Note on Libxc:                                                       !
    ! For use of VV10 in combination with Libxc, use the optional          !
    ! arguments vv10_bval and vv10_Cval to pass in the values for b and C  !
    ! parameters output by Libxc. vdwdf_type should be 'VV10-LIBXC'.       !
    !----------------------------------------------------------------------!
    ! Written in Jan 2012 by Lampros Andrinopoulos.                        !
    ! OpenMP parallelised by Nicholas Hine in August 2013.                 !
    ! Support for arbitrary b and C values for use with Libxc added by     !
    ! James C. Womack, May 2016                                            !
    ! Added spin polarized support.                                        !
    ! Jolyon Aarons, 2021.                                                 !
    !======================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only : comms_reduce, pub_my_proc_id, pub_root_proc_id,&
          comms_bcast
    use constants, only: DP, PI
    use geometry, only: OPERATOR(.dot.), operator(.cross.)
    use fourier, only: fourier_apply_cell_forward, fourier_apply_cell_backward
    use rundat, only: pub_debug_on_root, pub_num_spins
!$  use rundat, only: pub_threads_max
    use simulation_cell, only: CELL_INFO
    use vdw_df_kernel, only: vdw_df_kernel_write
    use utils, only: utils_abort, utils_assert

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    real(kind=DP), intent(in) :: density_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=DP),intent(in) :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_num_spins)
    real(kind=DP), intent(out) :: nl_pot_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins) !Non-local potential
    real(kind=DP), intent(out) :: nlc_energy_real ! Non-local correlation energy
    character(*), intent(in) :: vdwdf_type
    ! JCW: Optional arguments which allow the VV10 b and C parameters to be
    ! JCW: set outside of this routine (currently only used for Libxc, when
    ! JCW: vdwdf_type = 'VV10-LIBXC')
    real(kind=DP), intent(in), optional :: vv10_bval
    real(kind=DP), intent(in), optional :: vv10_Cval

    ! Local Variables
    integer, parameter :: kernel_file = 777
    integer :: ipt,i3,i2,islab23,i1,islab12,is
    integer :: cart
    integer :: alpha, beta
    integer :: N_alpha
    integer :: N_radial
    integer :: ierr
    logical :: fileexists
    character(LEN=30)  :: double_format = "(1p4e23.14)"
    real(kind=DP) :: g(3)
    real(kind=dp) :: modk, cell_volume
    real(kind=DP) :: fac,n_2
    real(kind=dp) :: max_radius, dg, q_c, q_min
    real(kind=dp), allocatable  :: u_real(:,:,:,:),nl_pot_aux(:,:,:,:),&
         dq_drho(:,:,:,:),dq_dgradrho(:,:,:,:),pot_aux(:,:,:)
    real(kind=dp), allocatable  :: kernel_recip(:,:)
    real(kind=dp), allocatable :: q_points(:)
    real(kind=DP), allocatable :: q(:,:,:)
    real(kind=dp), allocatable :: kernel(:,:,:), kernel2(:,:,:)
    complex(kind=DP), allocatable  :: theta_recip(:,:,:,:),&
        pot_aux_recip(:,:,:), theta_recip_pt(:)
    real(kind=DP) :: nlc_energy
    real(kind=DP) :: dtheta_drho, dtheta_dgradrho,&
        A,B,C,D,E,F,h,den,p_pt,dpdq_pt
    real(kind=DP), allocatable :: poly(:), poly2(:,:)
    integer :: hi,low,bisec
    real(kind=DP), parameter :: THPI = 3.0_DP * PI
    real(kind=DP), parameter :: NPISIXR = (9.0_DP * PI)**(1.0_DP/6.0_DP)
    real(kind=DP) :: Cval,bval
    real(kind=DP) :: beta_val ! JCW: Derived from bval
    real(kind=DP), parameter :: dentol=1.0E-15_DP

    real(kind=dp) :: dummy

    ! JCW: Check that optional parameter-setting arguments for VV10 are used
    ! JCW: correctly (i.e. only when vdwdf_type == 'VV10-LIBXC').
    if (vdwdf_type == 'VV10-LIBXC') then
       call utils_assert(present(vv10_bval).and.present(vv10_Cval),&
            "Error in nlxc_vdw_energy: Optional arguments vv10_bval and &
            &vv10_Cval must be present for VV10 to be combined with a Libxc &
            &functional.")
    else
       call utils_assert(.not.(present(vv10_bval).or.present(vv10_Cval)),&
            "Error in nlxc_vdw_energy: Optional arguments vv10_bval and &
            &vv10_Cval must not be present for if vdwdf_type != 'VV10-LIBXC.")
    end if


    if(vdwdf_type == 'VV10') then
       Cval = 0.0093_DP
       bval = 6.3_DP ! Sabatini (rVV10)
    elseif(vdwdf_type == 'VV10-S') then
       Cval = 0.000001_DP
       bval = 10.25_DP ! Bjorkman parameters for AM05-VV10sol
    elseif(vdwdf_type == 'VV10BC') then
       ! JCW: Original b and C parameters used by Vydrov and Van Voorhis.
       ! JCW: See:
       ! JCW: O.A. Vydrov & T. V. Voorhis, J. Chem. Physs 133, 244103 (2010).
       Cval = 0.0093_DP
       bval = 5.9_DP
    elseif(vdwdf_type == 'VV10-B97M-V') then
       ! JCW: Mardirossian and Head-Gordon parameters for B97M-V functional
       ! JCW: (identical to omegaB97X-V parameterization).
       ! JCW: Mardirossian and Head-Gordon used b and C values for the
       ! JCW: original VV10 functional (Vydrov and Voorhis) when optimizing
       ! JCW: B97M-V, rather than Sabatini's (slightly modified) rVV10, so we
       ! JCW: would expect this to cause some discrepancy in results reported
       ! JCW: by ONETEP and other programs.
       ! JCW: See:
       ! JCW: O.A. Vydrov & T. V. Voorhis, J. Chem. Physs 133, 244103 (2010).
       ! JCW: N. Mardirossian & M. Head-Gordon,
       ! JCW:                              J. Chem. Phys. 142, 074111 (2015).
       Cval = 0.01_DP
       bval = 6.0_DP
    elseif(vdwdf_type == 'VV10-LIBXC') then
       ! JCW: Libxc provides the values of bval and Cval for each functional
       ! JCW: These must be passed in using the optional arguments vv10_Cval
       ! JCW: and vv10_bval.
       Cval = vv10_Cval
       bval = vv10_bval
    endif


    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering nlxc_vdw_energy'


    !JA: Implemented spin polarization

    ! la: open the kernel table and read it in

    ! la: check if there is already a kernel file
    fileexists = .false.
    if((vdwdf_type == 'VDWDF1').or.(vdwdf_type == 'VDWDF2')) then
       if (pub_on_root) inquire(file='vdW_df_kernel',exist=fileexists)
    else
       if (pub_on_root) inquire(file='vdW_vv10_kernel',exist=fileexists)
    endif

    call comms_bcast(pub_root_proc_id,fileexists)

    ! la: if not, generate a new one (in parallel)
    if (.not.fileexists) call vdw_df_kernel_write(vdwdf_type)

    ! la: read kernel file on root proc
    call read_kernel_file()

    ! Start timer
    call timer_clock('nlxc_vdw_energy',1)

    q_c = q_points(N_alpha)
    q_min = q_points(1)
    dg = 2.0_DP*pi/max_radius

    ! la: allocate workspace arrays for nlxc calculation
    allocate(u_real(grid%ld1,grid%ld2,grid%max_slabs12,N_alpha),stat=ierr)
    call utils_alloc_check('nlxc_vdw_energy','u_real',ierr)
    allocate(poly2(N_alpha,N_alpha),stat=ierr)
    call utils_alloc_check('nlxc_vdw_energy','poly2',ierr)
    allocate(q(grid%ld1,grid%ld2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('nlxc_vdw_energy','q',ierr)
    allocate(dq_drho(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('nlxc_vdw_energy','dq_drho',ierr)
    allocate(dq_dgradrho(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('nlxc_vdw_energy','dq_dgradrho',ierr)
    allocate(theta_recip(grid%ld3,grid%ld2,grid%max_slabs23,N_alpha),stat=ierr)
    call utils_alloc_check('nlxc_vdw_energy','theta_recip',ierr)

    ! la: initialise arrays
    nlc_energy = 0.0_DP
    nl_pot_fine(:,:,:,:) = 0.0_DP

    ! la: determine q0, theta, dq_drho and dq_dgradrho
    call nlxc_q0_theta(density_fine,density_grad,grid,q,theta_recip,&
         dq_drho,dq_dgradrho,N_alpha,q_points,u_real,vdwdf_type,&
         vv10_bval = vv10_bval, vv10_Cval = vv10_Cval)

    ! la: calculate the nonlocal energy as a recip-space sum (Eq. 8 of Soler)
!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i2,i3,islab23,g,modk,fac,theta_recip_pt,kernel_recip, &
!$OMP      alpha,beta,ierr) &
!$OMP SHARED(grid,N_alpha,theta_recip,dg,kernel,N_radial, &
!$OMP      pub_my_proc_id,pub_threads_max) &
!$OMP REDUCTION(+:nlc_energy)

    ! ndmh: NB: theta_recip is re-used as u_recip in this section
    allocate(theta_recip_pt(N_alpha),stat=ierr)
    call utils_alloc_check('nlxc_vdw_energy','theta_recip_pt',ierr)
    allocate(kernel_recip(N_alpha,N_alpha),stat=ierr)
    call utils_alloc_check('nlxc_vdw_energy','kernel_recip',ierr)
    fac = 2.0_DP
    kernel_recip(:,:) = 0.0_DP

    !$OMP DO
    do ipt=1,grid%num_slabs23*grid%n2*grid%n3

       i3 = modulo(ipt-1,grid%n3) + 1
       i2 = modulo((ipt-i3)/grid%n3,grid%n2) + 1
       islab23 = (ipt-(i2-1)*grid%n3-i3) / (grid%n3*grid%n2) + 1
       if ((islab23==1).and.(grid%proc_slab23(1) == pub_my_proc_id)) &
            fac = 1.0_DP

       call cell_grid_recip_pt(g(:),islab23 + &
            grid%first_slab23(pub_my_proc_id) - 1,i2,i3,grid)
       modk = sqrt(sum(g(:)**2))
       call interpolate(modk,kernel_recip,N_radial,N_alpha,dg,kernel)
       theta_recip_pt(1:N_alpha) = theta_recip(i3,i2,islab23,1:N_alpha)
       ! la: loop over interpolating polynomials
       do alpha=1,N_alpha
          theta_recip(i3,i2,islab23,alpha) = cmplx(0.0_DP,0.0_DP,kind=DP)
          do beta=1,N_alpha

             theta_recip(i3,i2,islab23,alpha) = &
                  theta_recip(i3,i2,islab23,alpha) + &
                  kernel_recip(alpha,beta) * &
                  theta_recip_pt(beta)

             nlc_energy = nlc_energy + &
                  real(conjg(theta_recip_pt(alpha)) * &
                  theta_recip_pt(beta) * &
                  kernel_recip(alpha,beta) * fac,kind=DP)

          enddo
       end do
       fac = 2.0_DP
    end do
!$OMP END DO

    deallocate(kernel_recip,stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','kernel_recip',ierr)
    deallocate(theta_recip_pt,stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','theta_recip_pt',ierr)

!$OMP END PARALLEL

    do alpha=1,N_alpha
       call fourier_apply_cell_backward(u_real(:,:,:,alpha), &
            theta_recip(:,:,:,alpha),grid)
    enddo

    deallocate(theta_recip,stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','theta_recip',ierr)

    allocate(nl_pot_aux(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins),stat=ierr)
    call utils_alloc_check('nlxc_vdw_energy','nl_pot_aux',ierr)

    nl_pot_aux(:,:,:,:) = 0.0_DP

    ! la: Obtain the second derivatives of the interpolating polynomials on the
    ! mesh points
    call spline(q_points,poly2)

    ! la: Real-space loop to calculate the non-local potential (Eq. 10 of Soler)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i2,i1,islab12,den,low,hi,bisec,h,A,B,C,D,E,F,alpha, &
!$OMP   poly,p_pt,dpdq_pt,dtheta_drho,dtheta_dgradrho,ierr,is) &
!$OMP SHARED(grid,nl_pot_aux,nl_pot_fine,density_fine,q,q_points,N_alpha, &
!$OMP   dq_drho,dq_dgradrho,u_real,poly2,pub_threads_max,bval,Cval,vdwdf_type, &
!$OMP   pub_num_spins)

    allocate(poly(N_alpha),stat=ierr)
    call utils_alloc_check('nlxc_vdw_energy','poly',ierr)


    if( (vdwdf_type == 'VDWDF1') .or. (vdwdf_type == 'VDWDF2') ) then

!$OMP DO
       do is=1,pub_num_spins ! JA: Spin polarized version
          do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
             i1 = modulo(ipt-1,grid%n1) + 1
             i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
             islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1

             den = density_fine(i1,i2,islab12,is)
             low = 1
             hi = N_alpha
             ! la: Obtain the values of the p(r) and dp_dq(r) at every point
             ! using a cubic spline interpolation scheme.
             ! These are then used to calculate dtheta_drho and dtheta_dgradrho
             do
                if ((hi - low) <= 1 ) exit
                bisec = (hi+low)/2
                if ( q(i1,i2,islab12) > q_points(bisec) ) then
                   low = bisec
                else
                   hi = bisec
                end if
             end do

             h = q_points(hi)-q_points(low)

             A = (q_points(hi) - q(i1,i2,islab12))/h
             B = (q(i1,i2,islab12) - q_points(low))/h
             C = (A*A*A-A)*h*h/6.0_DP
             D = (B*B*B-B)*h*h/6.0_DP
             E = (3.0_DP*A*A-1.0_DP)*h/6.0_DP
             F = (3.0_DP*B*B-1.0_DP)*h/6.0_DP

             ! la: loop over interpolating polynomials
             do alpha=1,N_alpha

                ! la: construct interpolating polynomial and its derivatives
                poly = 0
                poly(alpha) = 1
                p_pt = A * poly(low) + B * poly(hi) &
                     + ((A*A*A-A) * poly2(alpha,low) + &
                     (B*B*B-B) * poly2(alpha,hi)) * h * h / 6.0_DP
                dpdq_pt = (poly(hi)-poly(low))/h - E*poly2(alpha,low)+&
                     F*poly2(alpha,hi)

                dtheta_drho = p_pt + den * dpdq_pt * dq_drho(i1,i2,islab12,is)
                dtheta_dgradrho = den * dpdq_pt * dq_dgradrho(i1,i2,islab12,is)

                ! la: calculate nonlocal potential at this point
                nl_pot_fine(i1,i2,islab12,is) = nl_pot_fine(i1,i2,islab12,is) + &
                     u_real(i1,i2,islab12,alpha)*dtheta_drho

                ! la: calculate partial term of the potential (White & Bird style)
                if(q(i1,i2,islab12) .ne. q_points(N_alpha)) then
                   nl_pot_aux(i1,i2,islab12,is) =  nl_pot_aux(i1,i2,islab12,is) + &
                        u_real(i1,i2,islab12,alpha) * dtheta_dgradrho
                endif

             enddo
          enddo
       end do
!$OMP END DO

    else

!$OMP DO
       do is=1,pub_num_spins ! JA: Spin polarized version
          do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
             i1 = modulo(ipt-1,grid%n1) + 1
             i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
             islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1

             den = density_fine(i1,i2,islab12,is)
             low = 1
             hi = N_alpha
             ! la: Obtain the values of the p(r) and dp_dq(r) at every point
             ! using a cubic spline interpolation scheme.
             ! These are then used to calculate dtheta_drho and dtheta_dgradrho
             do
                if ((hi - low) <= 1 ) exit
                bisec = (hi+low)/2
                if ( q(i1,i2,islab12) > q_points(bisec) ) then
                   low = bisec
                else
                   hi = bisec
                end if
             end do

             h = q_points(hi)-q_points(low)

             A = (q_points(hi) - q(i1,i2,islab12))/h
             B = (q(i1,i2,islab12) - q_points(low))/h
             C = (A*A*A-A)*h*h/6.0_DP
             D = (B*B*B-B)*h*h/6.0_DP
             E = (3.0_DP*A*A-1.0_DP)*h/6.0_DP
             F = (3.0_DP*B*B-1.0_DP)*h/6.0_DP

             ! la: loop over interpolating polynomials
             do alpha=1,N_alpha

                ! la: construct interpolating polynomial and its derivatives
                poly = 0
                poly(alpha) = 1
                p_pt = A * poly(low) + B * poly(hi) &
                     + ((A*A*A-A) * poly2(alpha,low) + &
                     (B*B*B-B) * poly2(alpha,hi)) * h * h / 6.0_DP
                dpdq_pt = (poly(hi)-poly(low))/h - E*poly2(alpha,low)+&
                     F*poly2(alpha,hi)

                ! JCW: Avoid divide by zero in dtheta_drho when den = 0.0_DP.
                ! JCW: Previously, if den<dentol, nl_pot_fine and
                ! JCW: nl_pot_aux for point (i1,i2,islab12) were skipped (they
                ! JCW: were previously set to zero), so moving the conditional
                ! JCW: to enclose assignment dtheta_drho and dthata_dgradrho
                ! JCW: should not change the result.
                if(den>dentol) then
                   dtheta_drho =  (0.75_DP*p_pt*(den**(-1.0_DP/4.0_DP)) + &
                        (den**(3.0_DP/4.0_DP)) * dpdq_pt * &
                        dq_drho(i1,i2,islab12,is)) * &
                        ((THPI*bval/(2.0_DP*NPISIXR))**(-3.0_DP/2.0_DP))

                   dtheta_dgradrho = (den**(3.0_DP/4.0_DP)) * dpdq_pt * &
                        dq_dgradrho(i1,i2,islab12,is) * &
                        ((THPI*bval/(2.0_DP*NPISIXR))**(-3.0_DP/2.0_DP))

                   nl_pot_fine(i1,i2,islab12,is) = nl_pot_fine(i1,i2,islab12,is) + &
                        u_real(i1,i2,islab12,alpha)*dtheta_drho

                   if(q(i1,i2,islab12) .ne. q_points(N_alpha)) then
                      nl_pot_aux(i1,i2,islab12,is) =  nl_pot_aux(i1,i2,islab12,is) + &
                           u_real(i1,i2,islab12,alpha) * dtheta_dgradrho
                   endif
                endif

             enddo
          enddo
       end do
!$OMP END DO

    endif

    deallocate(poly,stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','poly',ierr)

!$OMP END PARALLEL

    deallocate(dq_dgradrho,stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','dq_dgradrho',ierr)
    deallocate(dq_drho,stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','dq_drho',ierr)
    deallocate(q,stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','q',ierr)
    deallocate(poly2,stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','poly2',ierr)

    allocate(pot_aux_recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('nlxc_vdw_energy','pot_aux_recip',ierr)
    allocate(pot_aux(grid%ld1,grid%ld2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('nlxc_vdw_energy','pot_aux',ierr)

    ! la: Calculation of the White and Bird-like term of the non-local potential
    ! (Eq. 10 of Soler)
    do is=1,pub_num_spins ! JA: Spin polarized version
       do cart=1,3
          pot_aux(:,:,:) = nl_pot_aux(:,:,:,is)*density_grad(:,:,:,cart,is)
          call fourier_apply_cell_forward(pot_aux,pot_aux_recip,grid)
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i2,i3,islab23,g,is) &
!$OMP SHARED(grid,cart,pot_aux_recip,pub_my_proc_id,pub_threads_max)
          do ipt=1,grid%num_slabs23*grid%n2*grid%n3
             i3 = modulo(ipt-1,grid%n3) + 1
             i2 = modulo((ipt-i3)/grid%n3,grid%n2) + 1
             islab23 = (ipt-(i2-1)*grid%n3-i3) / (grid%n3*grid%n2) + 1
             call cell_grid_recip_pt(g(:),islab23 + &
                  grid%first_slab23(pub_my_proc_id) - 1,i2,i3,grid)
             pot_aux_recip(i3,i2,islab23) = &
                  pot_aux_recip(i3,i2,islab23) * &
                  g(cart) * cmplx(0.0_DP,1.0_DP,kind=DP)
          enddo
!$OMP END PARALLEL DO
          call fourier_apply_cell_backward(pot_aux,pot_aux_recip,grid)
          nl_pot_fine(:,:,:,is) = nl_pot_fine(:,:,:,is) - pot_aux(:,:,:)
       enddo
    end do

    ! la: apply common factors and sum over all procs
    cell_volume = abs((cell%a1.cross.cell%a2).dot.cell%a3)
    n_2 = real(grid%n1*grid%n2*grid%n3,DP)*real(grid%n1*grid%n2*grid%n3,DP)
    nlc_energy=0.5_DP*nlc_energy*cell_volume/n_2
    nlc_energy_real=nlc_energy


    ! gcc32: add uniform-density compensation energy and potential
    ! for VV10 or AM05-VV10sol functinals
    ! JCW: (also the parameterization of VV10 for B97M-V)
    if( (vdwdf_type == 'VV10') .or. (vdwdf_type == 'VV10-S') .or. &
         (vdwdf_type == 'VV10-B97M-V').or.(vdwdf_type == 'VV10-LIBXC').or.&
         (vdwdf_type == 'VV10BC') ) then
       ! JCW: beta_val is determined by the bval parameter
       ! JCW:   beta_val = (1/32) * (3/bval**2)**0.75
       ! JCW: See Eq. 11 in Van Voohis and Vydrov, JCP, 2010.
       beta_val = 0.03125_DP * ((3.0_DP / (bval*bval) )**(0.75_DP))
       do is = 1,pub_num_spins
          do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
             i1 = modulo(ipt-1,grid%n1) + 1
             i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
             islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1

             nlc_energy_real = nlc_energy_real + beta_val * cell_volume * &
                  density_fine(i1,i2,islab12,is) / sqrt(n_2)
             nl_pot_fine(i1,i2,islab12,is) = nl_pot_fine(i1,i2,islab12,is) + beta_val
          end do
       end do
    endif

    call comms_reduce('SUM',nlc_energy_real)


    ! la: deallocate workspace arrays
    deallocate(pot_aux,stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','pot_aux',ierr)
    deallocate(pot_aux_recip,stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','pot_aux_recip',ierr)
    deallocate(nl_pot_aux,stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','nl_pot_aux',ierr)
    deallocate(u_real,stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','u_real',ierr)
    deallocate(kernel2,stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','kernel2',ierr)
    deallocate(kernel,stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','kernel', ierr)
    deallocate(q_points, stat=ierr)
    call utils_dealloc_check('nlxc_vdw_energy','q_points', ierr)

    ! Stop timer
    call timer_clock('nlxc_vdw_energy',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving nlxc_vdw_energy'

contains

      !=====================================================!
      ! This subroutine reads the kernel file vdW_df_kernel !
      ! which has previously been written to disk by the    !
      ! routine vdw_df_kernel_write                         !
      !-----------------------------------------------------!
      ! Written in Nov 2012 by Lampros Andrinopoulos.       !
      !=====================================================!

      subroutine read_kernel_file()

        implicit none

        integer :: ierr
        integer :: alpha,beta

        if (pub_on_root) then
           if( (vdwdf_type == 'VDWDF1') .or. (vdwdf_type == 'VDWDF2') ) then
              open(unit=kernel_file, file='vdW_df_kernel', status='old', &
                  form='formatted', action='read',iostat=ierr)
              call utils_open_unit_check('nlxc_vdw_energy','vdW_df_kernel',ierr)
           else
              open(unit=kernel_file, file='vdW_vv10_kernel', status='old', &
                  form='formatted', action='read',iostat=ierr)
              call utils_open_unit_check('nlxc_vdw_energy','vdW_vv10_kernel',ierr)
           endif

           ! la: read sizes from kernel table
           read(kernel_file, '(2i5)') N_alpha, N_radial
           N_radial = N_radial - 1
           read(kernel_file, double_format) max_radius
        end if

        ! la: broadcast kernel to all procs
        call comms_bcast(pub_root_proc_id,N_alpha)
        call comms_bcast(pub_root_proc_id,N_radial)
        call comms_bcast(pub_root_proc_id,max_radius)

        ! la: allocation for kernel table
        allocate(q_points(N_alpha), stat=ierr)
        call utils_alloc_check('nlxc_vdw_energy','q_points', ierr)
        allocate(kernel(0:N_radial,N_alpha,N_alpha),stat=ierr)
        call utils_alloc_check('nlxc_vdw_energy','kernel', ierr)
        allocate(kernel2(0:N_radial,N_alpha,N_alpha),stat=ierr)
        call utils_alloc_check('nlxc_vdw_energy','kernel2',ierr)

        if (pub_on_root) then

           ! la: read number of q_points then read kernel data
           read(kernel_file, double_format) q_points
           do alpha = 1, N_alpha
              do beta = 1, alpha
                 read(kernel_file, double_format) kernel(0:N_radial, alpha, &
                      beta)
                 kernel(0:N_radial, beta, alpha) = kernel(0:N_radial, alpha, &
                      beta)
              end do
           end do

           do alpha = 1, N_alpha
              do beta = 1, alpha
                 read(kernel_file, double_format) kernel2(0:N_radial, alpha, &
                      beta)
                 kernel2(0:N_radial, beta, alpha) = kernel2(0:N_radial, alpha, &
                      beta)
              end do
           end do

           ! la: close kernel table file
           close(kernel_file,iostat=ierr)
           if( (vdwdf_type == 'VDWDF1') .or. (vdwdf_type == 'VDWDF2') ) then
              call utils_close_unit_check('vdW_df_kernel','vdW_df_kernel',ierr)
           else
              call utils_close_unit_check('vdW_df_kernel','vdW_vv10_kernel',ierr)
           endif

        end if

        call comms_bcast(pub_root_proc_id,q_points)
        call comms_bcast(pub_root_proc_id,kernel)
        call comms_bcast(pub_root_proc_id,kernel2)

      end subroutine read_kernel_file

      !=====================================================!
      ! This subroutine does a two-dimensional cubic spline !
      ! interpolation to return the kernel value at a       !
      ! particular g-vector.                                !
      !-----------------------------------------------------!
      ! Written in Nov 2012 by Lampros Andrinopoulos.       !
      !=====================================================!

      subroutine interpolate(g,kernel_recip,N_radial,N_alpha,dg,kernel)

      ! Arguments
      real(dp), intent(in) :: g
      real(dp), intent(inout) :: kernel_recip(:,:)
      real(dp), intent(in) :: dg
      integer, intent(in) :: N_alpha,N_radial
      real(dp), intent(in) :: kernel(0:N_radial,N_alpha,N_alpha)

      ! Local Variables
      integer :: alpha, beta, i
      real(dp) :: A, B

      if ( g >= N_radial*dg ) then
         call utils_abort('Error in vdW Kernel interpolation. K-point out &
              &of range')
      end if

      kernel_recip = 0.0_DP

      i = int(g/dg) ! Use this to find the right points to be used in the
                    ! bisection

      if (mod(g,dg)==0) then
         do alpha=1,N_alpha
            do beta=1,alpha

               kernel_recip(alpha, beta) = kernel(i,alpha, beta)
               kernel_recip(beta, alpha) = kernel(i,beta, alpha)

            end do
         end do
         return
      end if

      A = (i * dg + dg - g)/dg
      B = (g - dg*i)/dg

      do alpha = 1, N_alpha
         do beta = 1, alpha

            kernel_recip(alpha, beta) = A*kernel(i, alpha, beta) + &
                 B*kernel(i+1, alpha, beta) &
                 +((A*A*A-A)*kernel2(i, alpha, beta)  &
                 + (B*B*B-B)*kernel2(i+1, alpha, beta)) * &
                 dg*dg/6.0_DP

            kernel_recip(beta, alpha) = kernel_recip(alpha, beta)

         end do
      end do

    end subroutine interpolate

  end subroutine nlxc_vdw_energy

  subroutine nlxc_q0_theta(density_fine,density_grad,grid,q,theta_recip,&
       dq_drho,dq_dgradrho,N_alpha,q_points,theta,vdwdf_type,&
       vv10_bval,vv10_Cval)

    !====================================================================!
    ! This subroutine calculates and then "saturates" the quantity q0 as !
    ! described in Soler Eq 4 and Dion Eq 11.                            !
    ! Uses u_real from parent routine as workspace (here called theta)   !
    !--------------------------------------------------------------------!
    ! Written in Jan 2012 by Lampros Andrinopoulos.                      !
    ! OpenMP parallelised by Nicholas Hine in August 2013.               !
    ! Spin polarization support by Jolyon Aarons in May 2021.            !
    !====================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP, PI, stdout
    use fourier, only: fourier_apply_cell_forward
!$  use rundat, only: pub_threads_max
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=dp), intent(in)  :: density_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=dp), intent(in)  :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_num_spins)
    real(kind=dp), intent(out) :: q(grid%ld1,grid%ld2, &
         grid%max_slabs12)
    real(kind=dp), intent(out) :: dq_drho(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    real(kind=dp), intent(out) :: dq_dgradrho(grid%ld1,grid%ld2, &
         grid%max_slabs12, pub_num_spins)
    integer, intent(in) :: N_alpha
    real(kind=DP), intent(in) :: q_points(N_alpha)
    complex(kind=DP), intent(out) :: theta_recip(grid%ld3,grid%ld2, &
         grid%max_slabs23,N_alpha)
    real(kind=DP), intent(inout) :: theta(grid%ld1,grid%ld2,grid%max_slabs12, &
         N_alpha)
    character(*), intent(in) :: vdwdf_type
    ! JCW: Optional arguments which allow the VV10 b and C parameters to be
    ! JCW: set outside of this routine (currently only used for Libxc, when
    ! JCW: vdwdf_type = 'VV10-LIBXC')
    real(kind=DP), intent(in), optional :: vv10_bval
    real(kind=DP), intent(in), optional :: vv10_Cval

    ! Local Variables
    real(kind=DP) :: den
    real(kind=DP) :: grad(3)
    real(kind=DP) :: mgd
    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient

    !JA: spin polarization
    real(kind=DP) :: den1,den2
    real(kind=DP) :: grad1(3),grad2(3)
    real(kind=dp) :: mgd1,mgd2
    real(kind=dp) :: kf1,kf2,s1,s2
    real(kind=dp) :: dq0_dq1,dq0_dq2
    real(kind=dp) :: eq7brac1,eq7brac2
    real(kind=dp) :: qx1,qx2,q0x1,q0x2
    real(kind=dp) :: dqx_dden1,dqx_dden2
    real(kind=dp) :: c_pot1,c_pot2

    real(kind=DP) :: x_energy
    real(kind=DP) :: c_energy
    real(kind=DP) :: c_pot
    real(kind=DP) :: x_pot
    real(kind=DP) :: dex_drho
    real(kind=DP) :: dec_drho
    real(kind=DP) :: dq0_dq,rs,drs_drho
    real(kind=DP), parameter :: dentol=1.0E-15_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: Z_vdwdf1 = 0.8491_DP/9.0_DP
    real(kind=DP), parameter :: Z_vdwdf2 = 1.8870_DP/9.0_DP
    real(kind=DP), parameter :: xf = -0.458165293283142893475554_DP
    real(kind=DP) :: sqrtrs,ex,ec,dummy,Z
    real(kind=DP) :: A,alpha1,beta1,beta2,beta3,beta4
    integer :: islab12,i2,i1,ipt
    integer :: alpha
    real(kind=DP) :: p_pt, q_c
    real(kind=DP) :: omega_0,omega_g_sq,omega_p_sq,kappa
    real(kind=DP) :: Cval,bval
    real(kind=DP) :: q_min
    real(kind=DP), parameter :: THPI = 3.0_DP * PI
    real(kind=DP), parameter :: FPIOT = 4.0_DP * PI / 3.0_DP
    real(kind=DP), parameter :: NPISIXR = (9.0_DP * PI)**(1.0/6.0)
    real(kind=dp), parameter :: TWPWTHD = (2.0_dp**THIRD)

    dummy = 0.0_DP
    A=0.0310907_DP
    alpha1 = 0.21370_DP
    beta1=7.5957_DP
    beta2=3.5876_DP
    beta3=1.6382_DP
    beta4=0.49294_DP
    q_c = q_points(N_alpha)
    q_min = q_points(1)

    q(:,:,:) = 0.0_DP ! jd: Takes care of padding between pt and ld

    if(vdwdf_type == 'VDWDF1') then
        Z = Z_vdwdf1
    elseif(vdwdf_type == 'VDWDF2') then
        Z = Z_vdwdf2
    elseif(vdwdf_type == 'VV10') then
       Cval = 0.0093_DP
       bval = 6.3_DP ! Sabatini (rVV10)
    elseif(vdwdf_type == 'VV10BC') then
       ! JCW: Original b and C parameters used in Vydrov and Van Voorhis'
       ! JCW: paper.
       Cval = 0.0093_DP
       bval = 5.9_DP
    elseif(vdwdf_type == 'VV10-S') then
       Cval = 0.000001_DP
       bval = 10.25_DP ! Bjorkmann parameters for AM05-VV10sol
    elseif(vdwdf_type == 'VV10-B97M-V') then
       ! JCW: Mardirossian and Head-Gordon parameters for B97M-V functional
       ! JCW: (see note at start of nlxc_vdw_energy routine)
       Cval = 0.01_DP
       bval = 6.0_DP
    elseif(vdwdf_type == 'VV10-LIBXC') then
       ! JCW: Libxc provides the values of bval and Cval for each functional
       Cval = vv10_Cval
       bval = vv10_bval
    else
        write(stdout,'(a)') 'Unrecognised vdW-DF type'
    endif



    if( (vdwdf_type == 'VDWDF1') .or. (vdwdf_type == 'VDWDF2') ) then
       if(pub_num_spins==1) then
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,den,grad,mgd,kf,s,rs,sqrtrs,drs_drho,dummy, &
!$OMP      dex_drho,dec_drho,ex,ec,x_energy,c_energy,x_pot,c_pot,dq0_dq) &
!$OMP SHARED(grid,A,alpha1,beta1,beta2,beta3,beta4,q_c,Z, &
!$OMP      q,dq_drho,dq_dgradrho,density_fine,density_grad,pub_threads_max, &
!$OMP      pub_num_spins)
          do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
             i1 = modulo(ipt-1,grid%n1) + 1
             i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
             islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1

             den = density_fine(i1,i2,islab12,1)
             grad(1:3) = density_grad(i1,i2,islab12,1:3,1)
             mgd = sqrt(grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3))

             ! ndmh: improved trap for floating point issues
             if (den>dentol/real(pub_num_spins,dp)) then
                kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
                s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient
                rs = (TONFPI / den)**THIRD
                sqrtrs = sqrt(rs)

                call xc_pw91_eq10(A,alpha1,beta1,beta2,beta3,beta4,rs,sqrtrs, &
                     dummy,dec_drho)

                drs_drho = (-THIRD)*((TONFPI)**THIRD)*(den**(-4.0_DP*THIRD))
                dex_drho=THIRD*xf*(TONFPI**(-THIRD))*(den**(-2.0_DP*THIRD))
                dec_drho = dec_drho * drs_drho

                call xc_lda_x_point(den,x_energy,x_pot)

                call xc_pw92_c_point(den,c_energy,c_pot)

                ex = x_energy / den
                ec = c_energy / den

                q(i1,i2,islab12)=(1.0_DP + ec/ex + &
                     Z * s * s) * kf

                dummy=q(i1,i2,islab12)
                call saturate(q_c,dummy,dq0_dq)
                dq_drho(i1,i2,islab12,1) = dq0_dq *( q(i1,i2,islab12)/(3.0_DP*den) &
                     + kf * (dec_drho/ex-dex_drho*ec/(ex*ex)&
                     - 2.0_DP*THIRD*Z*((THPISQ)**(-2.0_DP*THIRD))*&
                     mgd*mgd*den**(-11.0_DP*THIRD)))

                dq_dgradrho(i1,i2,islab12,1) = dq0_dq*Z/(2.0_DP*den*den*kf)

                call saturate(q_c, q(i1,i2,islab12),dummy)

             else
                kf = 0.0_DP
                s = 0.0_DP
                rs = 0.0_DP
                sqrtrs = 0.0_DP
                dummy = 0.0_DP
                dec_drho = 0.0_DP
                dex_drho = 0.0_DP
                ex = 0.0_DP
                ec = 0.0_DP
                q(i1,i2,islab12) = q_c
                dq_drho(i1,i2,islab12,1) = 0.0_DP
                dq_dgradrho(i1,i2,islab12,1) = 0.0_DP
             end if

          end do
!$OMP END PARALLEL DO
       else ! JA: Spin polarized version.
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,den,rs,sqrtrs,drs_drho,dummy, &
!$OMP      den1,den2,grad1,grad2,mgd1,mgd2,kf1,kf2,s1,s2,dq0_dq1,dq0_dq2, &
!$OMP      eq7brac1,eq7brac2,qx1,qx2,q0x1,q0x2,dqx_dden1,dqx_dden2,c_pot1,c_pot2, &
!$OMP      dex_drho,dec_drho,ex,ec,x_energy,c_energy,x_pot,c_pot,dq0_dq) &
!$OMP SHARED(grid,A,alpha1,beta1,beta2,beta3,beta4,q_c,Z, &
!$OMP      q,dq_drho,dq_dgradrho,density_fine,density_grad,pub_threads_max, &
!$OMP      pub_num_spins)
          do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
             i1 = modulo(ipt-1,grid%n1) + 1
             i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
             islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1
             den1 = density_fine(i1,i2,islab12,1)
             den2 = density_fine(i1,i2,islab12,2)
             den = den1 + den2

             if (den>dentol) then

                if (den1>dentol/real(pub_num_spins,dp)) then
                   grad1(1:3) = density_grad(i1,i2,islab12,1:3,1)
                   mgd1 = sqrt(grad1(1)*grad1(1)+grad1(2)*grad1(2)+grad1(3)*grad1(3))
                   kf1 = (THPISQ * den1)**THIRD   ! Fermi wave-vector
                   s1 = mgd1 / (2.0_DP*den1*kf1)    ! Dimensionless gradient
                   eq7brac1 = 1.0_dp + Z*(s1/TWPWTHD)**2  ! JA : bracket from equation 7, Thonhauser 2015

                   qx1 = kf1*TWPWTHD * eq7brac1
                   q0x1=qx1
                   call saturate(q_c,q0x1,dq0_dq1)

                   dqx_dden1 = ((-s1/TWPWTHD * FTHRD)*(kf1*TWPWTHD) * &
                        (2.0_dp * s1/TWPWTHD * Z) + (THIRD*kf1*TWPWTHD*eq7brac1))*dq0_dq1 + q0x1*den1/den

                   dqx_dden2 = - q0x1*den1/den

                else
                   den1 = 0.0_dp
                   grad1(1:3) = 0.0_dp
                   mgd1 = 0.0_dp
                   kf1 = 0.0_dp
                   s1 = 0.0_dp
                   eq7brac1 = 1.0_dp
                   qx1 = 0.0_dp
                   q0x1 = 0.0_dp
                   dqx_dden1 = 0.0_dp
                   dqx_dden2 = 0.0_dp
                endif

                if (den2>dentol/real(pub_num_spins,dp)) then
                   grad2(1:3) = density_grad(i1,i2,islab12,1:3,2)
                   mgd2 = sqrt(grad2(1)*grad2(1)+grad2(2)*grad2(2)+grad2(3)*grad2(3))
                   kf2 = (THPISQ * den2)**THIRD   ! Fermi wave-vector
                   s2 = mgd2 / (2.0_DP*den2*kf2)    ! Dimensionless gradient
                   eq7brac2 = 1.0_dp + Z*(s2/TWPWTHD)**2

                   qx2 = kf2*TWPWTHD * eq7brac2
                   q0x2=qx2
                   call saturate(q_c,q0x2,dq0_dq2)

                   dqx_dden2 = dqx_dden2 + ((-s2/TWPWTHD * FTHRD)*(kf2*TWPWTHD) * &
                        (2.0_dp * s2/TWPWTHD * Z) + (THIRD*kf2*TWPWTHD*eq7brac2))*dq0_dq2 + q0x2*den2/den

                   dqx_dden1 = dqx_dden1 - q0x2*den2/den
                else
                   den2 = 0.0_dp
                   grad2(1:3) = 0.0_dp
                   mgd2 = 0.0_dp
                   kf2 = 0.0_dp
                   s2 = 0.0_dp
                   eq7brac2 = 1.0_dp
                   qx1 = 0.0_dp
                   q0x1 = 0.0_dp
                end if

                ! JA: I did this slightly differently to the non-spin polarized version
                !     The exchange energy is included directly (like the paper) rather than
                !     calling x_point. This could be changed trivially if necessary.
                !     Also, rather than getting the gradients from the eq10 routine,
                !     I opted to use the pots from c_point_sp... since the spin polarized
                !     interpolation formula using eq10 would get messy.

                ! JA: correlation energy for q_c contribution
                call xc_pw92_c_point_sp(den,den1,den2,c_energy,c_pot1,c_pot2)
                ec = c_energy / den

                ! JA: sum q_x and q_c.
                q(i1,i2,islab12) = (q0x1*den1 + q0x2*den2)/den -FTHRD*pi*ec
                call saturate(q_c,q(i1,i2,islab12),dq0_dq)

                ! JA: dqx_dden holds exchange and correlation comes from xc_pw92_c_point_sp
                dq_drho(i1,i2,islab12,1) = (dqx_dden1 + FTHRD*(ec - c_pot1))*dq0_dq
                dq_drho(i1,i2,islab12,2) = (dqx_dden2 + FTHRD*(ec - c_pot2))*dq0_dq

                if (den1>dentol/real(pub_num_spins,dp)) then
                   dq_dgradrho(i1,i2,islab12,1) = 2.0_dp * dq0_dq * dq0_dq1 * den1 * kf1*TWPWTHD * &
                        2.0_dp*(s1/TWPWTHD)*Z / (2.0_dp * kf1*TWPWTHD * den1*2.0_dp)
                else
                   dq_dgradrho(i1,i2,islab12,1) = 0.0_dp
                end if

                if (den2>dentol/real(pub_num_spins,dp)) then
                   dq_dgradrho(i1,i2,islab12,2) = 2.0_dp * dq0_dq * dq0_dq2 * den2 * kf2*TWPWTHD * &
                        2.0_dp*(s2/TWPWTHD)*Z / (2.0_dp * kf2*TWPWTHD * den2*2.0_dp)
                else
                   dq_dgradrho(i1,i2,islab12,2) = 0.0_dp
                end if

             else
                kf1 = 0.0_DP
                kf2 = 0.0_DP
                s1 = 0.0_DP
                s2 = 0.0_DP
                eq7brac1 = 1.0_dp
                eq7brac2 = 1.0_dp
                qx1 = 0.0_dp
                qx2 = 0.0_dp
                q0x1 = 0.0_dp
                q0x2 = 0.0_dp
                dec_drho = 0.0_DP
                dex_drho = 0.0_DP
                ex = 0.0_DP
                ec = 0.0_DP
                q(i1,i2,islab12) = q_c
                dq_drho(i1,i2,islab12,1) = 0.0_DP
                dq_drho(i1,i2,islab12,2) = 0.0_DP
                dq_dgradrho(i1,i2,islab12,1) = 0.0_DP
                dq_dgradrho(i1,i2,islab12,2) = 0.0_DP
             end if
          end do
!$OMP END PARALLEL DO
       end if

    else ! gcc32: for VV10 and AM05-VV10sol
       !JA: Now spin polarized

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,den,grad,mgd,s,kappa,omega_g_sq,omega_0, &
!$OMP      omega_p_sq,dummy, dq0_dq, mgd1, mgd2, den1,den2,grad1,grad2) &
!$OMP SHARED(grid,q_c,bval,Cval,q,dq_drho,dq_dgradrho,density_fine, &
!$OMP      q_min,density_grad,pub_threads_max,pub_num_spins)
       do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
          i1 = modulo(ipt-1,grid%n1) + 1
          i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
          islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1

          ! JA: Spin polarization support
          den1 = density_fine(i1,i2,islab12,1)
          grad1(1:3) = density_grad(i1,i2,islab12,1:3,1)
          mgd1 = sqrt(grad1(1)*grad1(1)+grad1(2)*grad1(2)+grad1(3)*grad1(3))
          if(pub_num_spins>1) then
             den2 = density_fine(i1,i2,islab12,2)
             den = den1 + den2
             grad2(1:3) = density_grad(i1,i2,islab12,1:3,2)
             mgd2 = sqrt(grad2(1)*grad2(1)+grad2(2)*grad2(2)+grad2(3)*grad2(3))
             grad = grad1 + grad2
          else
             den = den1
             grad = grad1
          end if
          mgd = sqrt(grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3))


          ! ndmh: improved trap for floating point issues
          if (den>dentol) then
             s = mgd/den
             omega_g_sq = Cval * (s**4.0_DP)
             kappa = bval * (THPI/2.0_DP) * (den**(1.0_DP/6.0_DP)) / NPISIXR
             omega_p_sq = 4.0_DP * PI * den
             omega_0 = sqrt (omega_g_sq + omega_p_sq/3.0_DP)

             q(i1,i2,islab12) = omega_0 / kappa

             dummy=q(i1,i2,islab12)

             if(q(i1,i2,islab12) < q_min) then
                q(i1,i2,islab12) = q_min
             endif

             call saturate(q_c,dummy,dq0_dq)

             dq_drho(i1,i2,islab12,1) = (0.5_DP*(-4.0_DP*Cval*(s**4.0_DP)/den &
                  + FPIOT) / ( sqrt(Cval* (s**4.0_DP) + FPIOT*den ) * &
                  ( THPI*bval* (den**(1.0_DP/6.0_DP))/(2.0_DP*NPISIXR) ) ) - &
                  sqrt( Cval* (s**4.0_DP) + FPIOT*den) * NPISIXR * &
                  (den**(-7.0_DP/6.0_DP)) / (9.0_DP*PI*bval) ) *dq0_dq

             dq_dgradrho(i1,i2,islab12,1) = dq0_dq*2.0_DP*Cval*(s**3.0_DP) / &
                  ((THPI/(2.0_DP*NPISIXR))*bval*(den**(7.0_DP/6.0_DP)) * &
                  sqrt(Cval*(s**4.0_DP)+FPIOT*den) )

             ! JA: I think the only change to VV10 for spin polarization is to
             !     the dq_dgradrho s.
             if(pub_num_spins>1) then
                dq_drho(i1,i2,islab12,2) = dq_drho(i1,i2,islab12,1)
                dq_dgradrho(i1,i2,islab12,2) = dq_dgradrho(i1,i2,islab12,1) / mgd2
             end if
             dq_dgradrho(i1,i2,islab12,1) = dq_dgradrho(i1,i2,islab12,1) / mgd1

             ! gcc32: note that dq_dgradrho is actually dq_dgradrho /  mgd
             ! The same note applies to the vdwdf1 and vdwdf2

             call saturate(q_c, q(i1,i2,islab12),dummy)

          else
             q(i1,i2,islab12) = q_c
             dq_drho(i1,i2,islab12,1) = 0.0_DP
             dq_dgradrho(i1,i2,islab12,1) = 0.0_DP
             if(pub_num_spins>1) then
                dq_drho(i1,i2,islab12,2) = dq_drho(i1,i2,islab12,1)
                dq_dgradrho(i1,i2,islab12,2) = dq_dgradrho(i1,i2,islab12,1)
             end if
          end if
       end do
!$OMP END PARALLEL DO

    endif ! end if clause for NLC type

    ! la: get the interpolating polynomials
    call spline_interp(q_points,q,theta)

    if(pub_num_spins>1) then
       if( (vdwdf_type == 'VDWDF1') .or. (vdwdf_type == 'VDWDF2') ) then
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,den,den1,den2,p_pt,alpha) &
!$OMP SHARED(grid,theta,density_fine,N_alpha,pub_threads_max)
          do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
             i1 = modulo(ipt-1,grid%n1) + 1
             i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
             islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1
             ! la: loop over interpolating polynomials
             den1 = density_fine(i1,i2,islab12,1)
             den2 = density_fine(i1,i2,islab12,2)
             den = den1 + den2
             do alpha=1,N_alpha
                p_pt = theta(i1,i2,islab12,alpha)
                theta(i1,i2,islab12,alpha) = den * p_pt
             enddo
          enddo
!$OMP END PARALLEL DO
       else
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,den,den1,den2,p_pt,alpha) &
!$OMP SHARED(grid,theta,density_fine,N_alpha,pub_threads_max,bval)
          do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
             i1 = modulo(ipt-1,grid%n1) + 1
             i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
             islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1
             den1 = density_fine(i1,i2,islab12,1)
             den2 = density_fine(i1,i2,islab12,2)
             den = den1 + den2
             do alpha=1,N_alpha
                p_pt = theta(i1,i2,islab12,alpha)
                if(den>dentol) then
                   theta(i1,i2,islab12,alpha) = (den**(3.0_DP/4.0_DP)) * &
                        p_pt *((bval*THPI/(2.0*NPISIXR))**(-3.0_DP/2.0_DP))
                else
                   theta(i1,i2,islab12,alpha) = 0.0_DP
                endif
             enddo
          enddo
!$OMP END PARALLEL DO
       endif
    else
       if( (vdwdf_type == 'VDWDF1') .or. (vdwdf_type == 'VDWDF2') ) then
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,den,p_pt,alpha) &
!$OMP SHARED(grid,theta,density_fine,N_alpha,pub_threads_max)
          do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
             i1 = modulo(ipt-1,grid%n1) + 1
             i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
             islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1
             den = density_fine(i1,i2,islab12,1)
             ! la: loop over interpolating polynomials
             do alpha=1,N_alpha
                p_pt = theta(i1,i2,islab12,alpha)
                theta(i1,i2,islab12,alpha) = den * p_pt
             enddo
          enddo
!$OMP END PARALLEL DO
       else
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,den,p_pt,alpha) &
!$OMP SHARED(grid,theta,density_fine,N_alpha,pub_threads_max,bval)
          do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
             i1 = modulo(ipt-1,grid%n1) + 1
             i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
             islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1
             do alpha=1,N_alpha
                den = density_fine(i1,i2,islab12,1)
                p_pt = theta(i1,i2,islab12,alpha)
                if(den>dentol) then
                   theta(i1,i2,islab12,alpha) = (den**(3.0_DP/4.0_DP)) * &
                        p_pt *((bval*THPI/(2.0*NPISIXR))**(-3.0_DP/2.0_DP))
                else
                   theta(i1,i2,islab12,alpha) = 0.0_DP
                endif
             enddo
          enddo
!$OMP END PARALLEL DO
       endif
    end if

    do alpha=1,N_alpha
       call fourier_apply_cell_forward(theta(:,:,:,alpha), &
            theta_recip(:,:,:,alpha),grid)
    enddo

contains

    subroutine saturate(q_c,q,dq0_dq)

      !====================================================!
      ! This subroutine is used to saturate the q quantity !
      ! according to Soler's scheme (Eq 4).                !
      !----------------------------------------------------!
      ! Written in Jan 2012 by Lampros Andrinopoulos.      !
      !====================================================!

      implicit none

      ! Arguments
      real(kind=DP), intent(inout) :: q
      real(kind=DP), intent(in) :: q_c
      real(kind=DP), intent(out) :: dq0_dq

      ! Local Variables
      integer :: m
      real(kind=DP) :: sum_exp, prod, sum_d

      ! Initialisation
      sum_exp = 0.0_DP
      sum_d = 0.0_DP
      prod = 1.0_DP

      do m=1,m_c-1
          prod = prod * (q/q_c)
          sum_d = sum_d + prod
      enddo

      sum_d = 1.0_DP + sum_d
      prod = 1.0_DP

      do m=1,m_c
         prod = prod * (q/q_c)
         sum_exp = sum_exp + prod/real(m,kind=DP)
      end do

      dq0_dq = sum_d*exp(-sum_exp)
      q = q_c * (1.0_DP - exp(-sum_exp))

    end subroutine saturate

  end subroutine nlxc_q0_theta

  subroutine spline_interp(X,Y,Y_interp)

  !====================================================================!
  ! This subroutine performs a cubic spline interpolation used to      !
  ! calculate the p(r) from the q(r) in Soler's method.                !
  ! The method and algorithm used is from Numerical Recipes            !
  ! for Fortran.                                                       !
  ! The interpolating polynomials used are Kronecker deltas.           !
  !--------------------------------------------------------------------!
  ! Written in Jan 2012 by Lampros Andrinopoulos.                      !
  ! OpenMP parallelised by Nicholas Hine in August 2013.               !
  !====================================================================!

!$  use rundat, only: pub_threads_max

    ! Arguments
    real(dp), intent(in) :: X(:), Y(:,:,:)
    real(dp), intent(out) :: Y_interp(:,:,:,:)

    ! Local Variables
    integer :: Nx
    integer :: Ngrid_points_x
    integer :: Ngrid_points_y
    integer :: Ngrid_points_z
    integer :: low,hi,bisec,i,ipt,islab12,i2,i1,ierr
    real(dp) :: A,B,C,D,E,F,h
    real(dp), allocatable :: poly(:),poly2(:,:)

    Nx = size(X)
    Ngrid_points_x = size(Y,1)
    Ngrid_points_y = size(Y,2)
    Ngrid_points_z = size(Y,3)

    allocate(poly2(Nx,Nx),stat=ierr)
    call utils_alloc_check('spline_interp','poly2',ierr)

    call spline(X,poly2)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,low,hi,bisec,h,A,B,C,D,E,F,poly,ierr) &
!$OMP SHARED(Ngrid_points_z,Ngrid_points_y,Ngrid_points_x,Nx,X,Y,Y_interp, &
!$OMP      poly2,pub_threads_max)

    ! la: allocate workspace
    allocate(poly(Nx),stat=ierr)
    call utils_alloc_check('spline_interp','poly',ierr)

!$OMP DO
    do ipt=1,Ngrid_points_z*Ngrid_points_y*Ngrid_points_x
       i1 = modulo(ipt-1,Ngrid_points_x) + 1
       i2 = modulo((ipt-i1)/Ngrid_points_x,Ngrid_points_y) + 1
       islab12 = (ipt-(i2-1)*Ngrid_points_x-i1) / &
            (Ngrid_points_x*Ngrid_points_y) + 1

       low = 1
       hi = Nx

       do
          if ((hi - low) <= 1 ) exit

          bisec = (hi+low)/2

          if ( Y(i1,i2,islab12) > X(bisec) ) then
             low = bisec
          else
             hi = bisec
          end if

       end do

       h = X(hi)-X(low)

       A = (X(hi) - Y(i1,i2,islab12))/h
       B = (Y(i1,i2,islab12) - X(low))/h
       C = (A*A*A-A)*h*h/6.0_DP
       D = (B*B*B-B)*h*h/6.0_DP
       E = (3.0_DP*A*A-1.0_DP)*h/6.0_DP
       F = (3.0_DP*B*B-1.0_DP)*h/6.0_DP

       do i = 1, Nx

          poly = 0
          poly(i) = 1

          Y_interp(i1,i2,islab12, i) = A * poly(low) + B * poly(hi) &
               + ((A*A*A-A) * poly2(i,low) + (B*B*B-B) * poly2(i, hi)) * &
               h * h / 6.0_DP
       end do
    end do
!$OMP END DO

    ! la: deallocate workspace
    deallocate(poly,stat=ierr)
    call utils_dealloc_check('spline_interp','poly',ierr)

!$OMP END PARALLEL

    deallocate(poly2,stat=ierr)
    call utils_dealloc_check('spline_interp','poly2',ierr)

  end subroutine spline_interp

  subroutine spline(X,Y2)

    !====================================================================!
    ! This subroutine is used to calculate the second derivatives of the !
    ! interpolating polynomials, using the method and algorithm from     !
    ! Numerical Recipes for Fortran, applied to Kronecker deltas.        !
    ! These second derivatives are used in the subroutine spline_interp. !
    !--------------------------------------------------------------------!
    ! Written in Jan 2012 by Lampros Andrinopoulos.                      !
    !====================================================================!

    ! Arguments
    real(DP), intent(in)  :: X(:)
    real(DP), intent(inout) :: Y2(:,:)

    ! Local Variables
    integer :: i, j, N, ierr
    real(DP) :: S, P
    real(DP), allocatable :: U(:), Y(:)

    N = size(X)

    ! la: allocate workspace
    allocate(U(N),stat=ierr)
    call utils_alloc_check('spline','U',ierr)
    allocate(Y(N),stat=ierr)
    call utils_alloc_check('spline','Y',ierr)

    do i=1,N

       Y(:) = 0.0_DP
       Y(i) = 1.0_DP

       Y2(i,1) = 0.0_DP
       U(1) = 0.0_DP

       do j = 2, N-1

          S = (X(j)-X(j-1))/(X(j+1)-X(j-1))
          P = S * Y2(i,j-1) + 2.0_DP
          Y2(i,j) = (S - 1.0_DP)/P
          U(j) =  (6.0_DP * ((Y(j+1)-Y(j))/(X(j+1)-X(j))-(Y(j)-Y(j-1))/ &
               (X(j)-X(j-1)))/(X(j+1)-X(j-1))-S*U(j-1))/P

       enddo

       Y2(i,N) =  0.0_DP

       do j= N-1, 1, -1

          Y2(i,j) = Y2(i,j+1) * Y2(i,j) + U(j)

       enddo

    enddo

    ! la: deallocate workspace
    deallocate(U,stat=ierr)
    call utils_dealloc_check('spline','U',ierr)
    deallocate(Y,stat=ierr)
    call utils_dealloc_check('spline','Y',ierr)

  end subroutine spline

end module nlxc
