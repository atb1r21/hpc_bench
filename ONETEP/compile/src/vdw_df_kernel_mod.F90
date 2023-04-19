!-*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!==============================================================================!
!                                                                              !
!  ONETEP vdW-DF kernel module: vdW_df_kernel.F90                              !
!                                                                              !
!  The subroutines in this file were written by Lampros Andrinopoulos          !
!  with help from Nicholas Hine, in May to June 2012.                          !
!                                                                              !
!  VV10 (and VV10sol) non-local correlation-energy functionality was written by!
!  Gabriel Constantinescu in July 2014                                         !
!                                                                              !
!  The routines in this module implement the calculation of the kernel phi     !
!  described in                                                                !
!   M. Dion, H. Rydberg, E. Schroder, D.C. Langreth, and B.I. Lundqvist        !
!   Phys. Rev. Lett. 92, 246401 (2004)                                         !
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
!==============================================================================!

module vdw_df_kernel

  use constants, only : DP, pi, stdout
  use services, only : services_regular_transform
  use comms, only : pub_my_proc_id, pub_total_num_procs,&
       pub_on_root, comms_reduce
  use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort,&
       utils_open_unit_check, utils_close_unit_check

  implicit none

  integer, parameter :: radial_points = 1025
  integer, parameter :: N_alpha = 20
  integer, parameter :: file_unit = 555
  real(DP), parameter :: max_radius = 100.0_DP
  real(DP), parameter :: small = 1.0D-14
  integer, parameter :: n = 256
  real(DP) :: lo=0.0_DP
  real(DP) :: hi = 64.0_DP
  real(DP), allocatable :: x_vector(:)
  real(DP), allocatable :: wgt(:)
  real(DP) :: a, b
  integer :: i, j, l, alpha, beta, N_ab
  integer :: index, index_start, index_finish, proc, procs, ierr

  real(DP) :: d1, d2, dr, r, k, dk, test, kmax, phi_1
  real(DP) , allocatable :: phi(:,:,:), phi2(:,:,:), w_ab(:,:), phi_k(:)

  real(dp), dimension(N_alpha):: q_points = (/  1.00D-5, 0.0449420825586261D0,&
  0.0975593700991365D0, &
  0.159162633466142D0, 0.231286496836006D0, 0.315727667369529D0, &
  0.414589693721418D0, &
  0.530335368404141D0, 0.665848079422965D0, 0.824503639537924D0, &
  1.010254382520950D0, &
  1.227727621364570D0, 1.482340921174910D0, 1.780437058359530D0, &
  2.129442028133640D0, &
  2.538050036534580D0, 3.016440085356680D0, 3.576529545442460D0, &
  4.232271035198720D0, &
  5.0D0 /)

  real(dp), dimension(N_alpha):: q_points_vv = (/ 1.0D-4, 3.0D-4, &
  5.893850845618885D-4, &
  1.008103720396345D-3, 1.613958359589310D-3, 2.490584839564653D-3, &
  3.758997979748929D-3, 5.594297198907115D-3, 8.249838297569416D-3, &
  1.209220822453922D-2, 1.765183095571029D-2, 2.569619042667097D-2, &
  3.733577865542191D-2, 5.417739477463518D-2, 7.854595729872216D-2, &
  0.113805449932145D0, 0.164823306218807D0, 0.238642339497217D0, &
  0.345452975434964D0, 0.5D0 /)



contains

  subroutine vdw_df_kernel_write(vdwdf_type)

    !======================================================================!
    ! This subroutine calculates the kernel phi according to               !
    ! Dion et al. Eq. 14, Fourier transforms it and writes it to a file.   !
    !----------------------------------------------------------------------!
    ! Written in Jun 2012 by Lampros Andrinopoulos.                        !
    !======================================================================!

    use rundat, only: pub_debug_on_root

    implicit none

    character(*), intent(in) :: vdwdf_type

    logical :: df1_or_2

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering vdw_df_kernel_write'

    dr = max_radius / real(radial_points-1,DP)
    dk = 2.0_DP * pi / max_radius
    kmax = real(radial_points-1,DP) * dk
    N_ab = N_alpha * (N_alpha + 1) / 2

    allocate(phi(radial_points,N_alpha,N_alpha),stat=ierr)
    call utils_alloc_check('vdw_df_kernel_write','phi',ierr)
    allocate(phi2(radial_points,N_alpha,N_alpha),stat=ierr)
    call utils_alloc_check('vdw_df_kernel_write','phi2',ierr)
    allocate(x_vector(n),stat=ierr)
    call utils_alloc_check('vdw_df_kernel_write','x_vector',ierr)
    allocate(wgt(n),stat=ierr)
    call utils_alloc_check('vdw_df_kernel_write','wgt',ierr)
    allocate(w_ab(n,n),stat=ierr)
    call utils_alloc_check('vdw_df_kernel_write','w_ab',ierr)
    allocate(phi_k(radial_points),stat=ierr)
    call utils_alloc_check('vdw_df_kernel_write','phi_k',ierr)

    ! la: initialise arrays
    phi(:,:,:) = 0.0_DP
    phi2(:,:,:) = 0.0_DP

    if( (vdwdf_type == 'VDWDF1') .or. (vdwdf_type == 'VDWDF2') ) then

       ! la: find the integration points and weights for the Gaussian quadrature
       call gauss_quad_points_weights(atan(lo), atan(hi), x_vector, wgt, n)

       do i=1,n
          x_vector(i) = tan(x_vector(i))
          wgt(i) = wgt(i) * (1.0_DP + x_vector(i) * x_vector(i))
       enddo

       do i=1,n
          do j=1,n

             a = x_vector(i)
             b = x_vector(j)

             w_ab(i,j) = 2.0_DP *((3.0_DP-a*a)*b*cos(b)*sin(a) &
                       +(3.0_DP-b*b)*a*cos(a)*sin(b) &
                       +(a*a+b*b-3.0_DP)*sin(a)*sin(b) &
                       -3.0_DP*a*b*cos(a)*cos(b))/(a*b)
          enddo
       enddo
    endif

    ! la: calculate the starting and finishing indices for each proc
    index = 1
    procs = pub_total_num_procs
    do proc=0,procs-1
       if (proc==pub_my_proc_id) index_start = index
       index = index + N_ab/procs
       if (proc < modulo(N_ab,procs)) then
          index = index + 1
       end if
       if (proc==pub_my_proc_id) index_finish = index - 1
    end do

    ! jd: Make the string comparison ahead of the loop for efficiency
    df1_or_2 = (vdwdf_type == 'VDWDF1') .or. (vdwdf_type == 'VDWDF2')

    do index = index_start,index_finish
       alpha = int((sqrt(1.0_DP+8.0_DP*real(index-1,DP))-1)*0.5_DP)+1
       beta = index - alpha*(alpha-1)/2
       do i = 2, radial_points
          r = real(i-1,DP) * dr
          if(df1_or_2) then
             d1 = r * q_points(alpha)
             d2 = r * q_points(beta)
             call get_phi(phi_1,d1,d2,n)
          else
             d1 = 1.0_DP + r*r*q_points_vv(alpha)
             d2 = 1.0_DP + r*r*q_points_vv(beta)
             phi_1 = -1.5_DP / (d1*d2*(d1+d2))
          endif
          phi(i,alpha,beta) = phi_1
       enddo

       call services_regular_transform(0,2,radial_points,max_radius,&
            radial_points,kmax,phi(:,alpha,beta),phi_k(:))
       phi(:,alpha,beta) = phi_k(:)
       call d2phi_dx2(phi(:,alpha,beta),phi2(:,alpha,beta),dk)
    enddo

    ! la: apply factors and sum over all procs
    phi = 4.0_DP * pi * phi
    phi2 = 4.0_DP * pi * phi2
    call comms_reduce('SUM',phi(:,:,:))
    call comms_reduce('SUM',phi2(:,:,:))

    ! la: deallocate workspace
    deallocate(x_vector,stat=ierr)
    call utils_dealloc_check('vdW_df_kernel','x_vector',ierr)
    deallocate(wgt,stat=ierr)
    call utils_dealloc_check('vdW_df_kernel','wgt',ierr)
    deallocate(w_ab,stat=ierr)
    call utils_dealloc_check('vdW_df_kernel','w_ab',ierr)

    ! la: tabulate kernel values and its second derivative in the vdW_df_kernel
    ! file
    if (pub_on_root) then

       if( (vdwdf_type == 'VDWDF1') .or. (vdwdf_type == 'VDWDF2') ) then
          open(unit=file_unit,file='vdW_df_kernel',action='write',iostat=ierr)
          call utils_open_unit_check('vdW_df_kernel','vdW_df_kernel',ierr)
       else
          open(unit=file_unit,file='vdW_vv10_kernel',action='write',iostat=ierr)
          call utils_open_unit_check('vdW_df_kernel','vdW_vv10_kernel',ierr)
       endif

       write(file_unit,'(2i5)') N_alpha, radial_points
       write(file_unit, '(1p,4e23.14)') max_radius
       if( (vdwdf_type == 'VDWDF1') .or. (vdwdf_type == 'VDWDF2') ) then
          write(file_unit, '(1p,4e23.14)') q_points
       else
          write(file_unit, '(1p,4e23.14)') q_points_vv
       endif
       do alpha = 1, N_alpha
          do beta = 1, alpha
             write(file_unit,'(1p,4e23.14)') phi(:, alpha,beta)
          enddo
       enddo
       do alpha = 1, N_alpha
          do beta = 1, alpha
             write(file_unit,'(1p,4e23.14)') phi2(:, alpha,beta)
          enddo
       enddo
       close(file_unit,iostat=ierr)
       if( (vdwdf_type == 'VDWDF1') .or. (vdwdf_type == 'VDWDF2') ) then
          call utils_close_unit_check('vdW_df_kernel','vdW_df_kernel',ierr)
       else
          call utils_close_unit_check('vdW_df_kernel','vdW_vv10_kernel',ierr)
       endif

    end if

    ! la: deallocate the kernel array and the array holding its second derivatives
    deallocate(phi,stat=ierr)
    call utils_dealloc_check('vdw_df_kernel_write','phi',ierr)
    deallocate(phi2,stat=ierr)
    call utils_dealloc_check('vdw_df_kernel_write','phi2',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving vdw_df_kernel_write'

contains

!------------------------------------------------------------------------------

    pure real(kind=DP) function t(w,x,y,z)

      implicit none

      real(DP), intent(in) :: w,x,y,z

      t = (1.0_DP/(w+x) + 1.0_DP/(y+z)) &
           * (1.0_DP/((w+y)*(x+z)) + 1.0_DP/((w+z)*(y+x)))

    end function t

!------------------------------------------------------------------------------

    pure real(kind=DP) function nu(y,d)

      use constants, only: PI
      implicit none

      real(kind=DP), intent(in) :: y, d

      real(kind=DP), parameter :: gamma = 4.0_DP * PI / 9.0_DP
      real(kind=DP) :: h, y2

      y2 = y**2
      h = 1.0_DP - exp(-gamma * y2 /(d**2))
      nu = y2 /(2.0_DP * h)

    end function nu

!------------------------------------------------------------------------------

    pure subroutine get_phi(phi,d1,d2,n)

      use constants, only: PI
      implicit none

      real(kind=DP), intent(out) :: phi
      real(kind=DP), intent(in) :: d1, d2
      integer, intent(in) :: n

      real(kind=DP) :: nu_temp1, nu_temp2, nu_temp3, nu_temp4, x, y
      integer :: i, j

      phi = 0.0_DP

      if (d1 /= 0.0_DP .or. d2 /= 0.0_DP) then
         do i=1,n
            x = x_vector(i)
            nu_temp1 = nu(x,d1)
            nu_temp3 = nu(x,d2)

            do j=1,n
               y = x_vector(j)
               nu_temp2 = nu(y,d1)
               nu_temp4 = nu(y,d2)
               phi = phi + wgt(i) * wgt(j) * w_ab(i,j) * &
                    t(nu_temp1, nu_temp2, nu_temp3, nu_temp4)
            end do
         end do
         phi = phi / (PI**2)
      endif

    end subroutine get_phi

!------------------------------------------------------------------------------

    subroutine d2phi_dx2(Y, Y2, dx)

      !======================================================================!
      ! This subroutine calculates the second derivative of the kernel       !
      ! by interpolation.                                                    !
      !----------------------------------------------------------------------!
      ! Written in Jan 2012 by Lampros Andrinopoulos.                        !
      !======================================================================!

      real(DP), intent(in)  :: Y(radial_points)
      real(DP), intent(inout) :: Y2(radial_points)
      real(DP), intent(in) :: dx

      ! Local Variables
      integer :: j
      real(DP) :: S, P
      real(DP) :: U(radial_points)

      Y2 = 0.0_DP
      U = 0.0_DP

      do j = 2, radial_points-1

         S = 0.5_DP
         P = S * Y2(j-1) + 2.0_DP
         Y2(j) = (S - 1.0_DP)/P
         U(j) =  (6.0_DP * ((Y(j+1)-Y(j))/(dx)-(Y(j)-Y(j-1))/ &
              (dx))/(2.0_DP*dx)-S*U(j-1))/P
      enddo

      Y2(radial_points) =  0.0_DP

      do j= radial_points-1, 1, -1

         Y2(j) = Y2(j+1) * Y2(j) + U(j)

      enddo

    end subroutine d2phi_dx2

!------------------------------------------------------------------------------

    subroutine gauss_quad_points_weights(lo, hi, points, wgt, n)

      !======================================================================!
      ! This subroutine calculates the integration points and weights        !
      ! to be used in the Gaussian quadrature integration scheme             !
      !----------------------------------------------------------------------!
      ! Written in Jun 2012 by Lampros Andrinopoulos.                        !
      !======================================================================!

      real(dp), intent(in) :: lo, hi

      real(dp), intent(out) :: points(:), wgt(:)

      integer, intent(in) :: n

      integer :: i, j, m

      real(DP) :: temp1, temp2, temp3, x0, Dx, xf, xmid, l

      m = (n+1)/2

      xmid = 0.5_DP * (lo + hi)
      l = 0.5_DP * (hi - lo)

      do i = 1, m

         x0 = cos(pi*(i-0.25_DP)/(n+0.5_DP))

         do
            temp1 = 1.0_DP
            temp2 = 0.0_DP

            do j = 1, n

               temp3 = temp2
               temp2 = temp1
               temp1 = ((2.0_DP * j - 1.0_DP)*x0*temp2 - (j-1.0_DP)*temp3)/j

            end do

            Dx = n * (x0 * temp1 -temp2)/(x0*x0 -1.0_DP)
            xf = x0
            x0 = xf - temp1/Dx

            if (abs(x0 - xf) <= small) exit

         end do

         points(i) = xmid - l * x0
         points(n+1-i) = xmid + l * x0

         wgt(i) = 2.0_DP * l/((1.0_DP - x0*x0) * Dx *Dx)
         wgt(n+1-i) = wgt(i)

      end do

    end subroutine gauss_quad_points_weights

  end subroutine vdw_df_kernel_write

end module vdw_df_kernel
