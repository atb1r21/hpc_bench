module mg_tests
  !contains driver subroutines for various tests, called in main
  !
  !Test descriptions
  !
  ! the following 3 tests were mainly used in DL_MG development
  !
  ! 1) parabolic_rho
  !    \laplace phi = rho
  !
  ! 2) sin_eps
  !    div( (1 + a * sin(pi*px*x) * sin(pi*py*y) * sin(pi*pz*z)) grad phi ) = rho
  !
  !  where rho is set such that
  ! phi = x*(1-x)*y*(1-y)*z*(1-z) is a solution with zero BC for both cases.
  !
  ! 3) linear_phi
  !  \laplace phi = 0
  ! with BC fixed such that (a1*x+b1)*(a2*y+b2)*(a3*z+b3) is a solution.
  !
  ! 4) pbe_nosteric
  ! [2)] + lambda \sum_i c(i) q(i) exp(-\beta phi)
  !
  ! 5) pbe_point_charge
  ! \laplace phi + (labmda \beta \sum_i c(i) q(i)^2) phi = 0
  !  the field being created by a point particle outside the box
  !
  ! Some other tests may be found below, but they are deprecated.
  !
  ! Lucian Anton, July 2013
  !
  use dl_mg_params, only : wp, pi, elcharge, kB => kboltz
  implicit none

  real(wp), parameter :: twopi = 2.0_wp * pi, s2pi = sqrt(twopi), &
       fourpi = 4.0_wp * pi

  type pbe_params_t
     integer n ! number of ion species
     real(wp), allocatable :: c(:), q(:), x(:) ! concentrations, charges, neutralisation ratios
     real(wp) lam, temp! ...
     logical use_steric, linearised, use_fas
     integer neutr_method
  end type pbe_params_t

  type gauss_erf_t
     real(wp) sig, mu, eps0, d, delta, scale
  end type gauss_erf_t

  type eps_sin_params_t
     ! parameters for eps_rel model : 1.0_wp + a * sin(2*pi*px*x) * sin(2*pi*py*y) * sin(2*pi*pz*z)
     integer :: px, py, pz ! wavelength number
     real(wp) :: a ! amplitude
  end type eps_sin_params_t


  type stw_params_t
     real(wp) :: d, delta
  end type stw_params_t

  type(pbe_params_t), target :: pbe_params
  type(pbe_params_t), pointer :: pbe
  type(gauss_erf_t) ge
  type(eps_sin_params_t) eps_sin_params
  type(stw_params_t) stw_params

  logical :: periodic(3) ! true if periodic
                         ! it is a module variable for the serial code

  ! statv is the conversion factor between SI V and  and Gauss statV
  ! 1 V = 1 statv
  real(wp), parameter :: rbohr = 5.2917720859e-9_wp, statv =  0.0033356409519815_wp, epsh2o=80.0_wp


  contains

    subroutine poisson_test(mg_comm,ngx,ngy,ngz,local_shift,nlx,nly,nlz, use_cg,&
         fd_maxiters, niter, model_name,&
         tol_res_rel, tol_res_abs, tol_pot_rel, &
         tol_pot_abs, tol_newton, tol_cg, tol_mg, &
         fd_order, use_damping, errors_return, &
         neutr_method, x_ratios)

!      use dl_mg_mpi_header
!$    use omp_lib
      use dl_mg_params,  only : wp,pi, hartree, &
           dl_mg_neutralise_none , &
           dl_mg_neutralise_with_jellium_uniform, &
           dl_mg_neutralise_with_jellium_vacc, &
           dl_mg_neutralise_with_ions_fixed, &
           dl_mg_neutralise_with_ions_auto, &
           dl_mg_neutralise_with_ions_auto_linear

      use mg_utils
      use dl_mg

      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: mg_comm,ngx,ngy,ngz,local_shift(3), &
                             nlx,nly,nlz, fd_maxiters, niter
      character(len=*),intent(in) :: model_name
      real(wp), intent(in) :: tol_res_rel, tol_res_abs, &
           tol_pot_rel, tol_pot_abs, tol_newton, tol_cg, tol_mg
      logical, intent(in) :: use_damping, errors_return, use_cg
      integer, intent(in) :: fd_order
      character(len=*), intent(in) :: neutr_method
      real(wp), intent(in) :: x_ratios(:)

      integer lbx,lby,lbz,ubx,uby,ubz
      integer imin,jmin,kmin,imax,jmax,kmax
      integer myid, pdims(3),my_coords(3), ngxyz(3), ierr
      character(len=128) fn_dump
      character(len=1) version_string_short
      character(len=DL_MG_VERSION_STRING_LENGTH) version_string

      !local arrays
      real(wp),allocatable, dimension(:,:,:) :: pot,rho,res
      real(wp),allocatable, dimension(:,:,:), target :: vsteric
      real(wp),allocatable, dimension(:,:,:) :: pot_out ! for _defco_ interface
      real(wp),allocatable, dimension(:,:,:,:) :: eps_half
      real(wp),allocatable, dimension(:,:,:)   :: eps_full
      real(wp) dx,dy,dz,xl,yl,zl,pert_amplitude

      ! global grid indices
      integer sx, ex, sy, ey, sz, ez, idx_s(3), idx_e(3), boundary_condition(3)
      ! model parameters
      integer px,py,pz
      ! power for z_to_power test
      integer potpow
      real(wp) a, alpha(3), beta(3), a11, a21, a31, a32, pointsource(3), astc, bstc
      real(wp) surface_pot, inv_debye, inv_debye2, alf
      logical :: use_eps_half=.true.
      character(len=DL_MG_MAX_ERROR_STRING) err_msg

      integer status


#ifdef MPI
      call mpi_comm_rank(mg_comm,myid,ierr)
      call mpi_cart_get(mg_comm, 3, pdims, periodic, my_coords, ierr)
#else
      myid  = 0
       pdims = (/ 1, 1, 1/)
       my_coords = (/ 0, 0, 0 /)
#endif

      if(myid == 0) then
         write(*,*)
         write(*,*) 'in poisson test ', trim(model_name)
         write(*,*) 'MPI topology ', pdims, periodic
!$       write(*,*) 'OpenMP threads     ', omp_get_max_threads()
         write(*,*) 'convergence tols ', tol_res_rel, tol_res_abs, tol_pot_rel, &
              tol_pot_abs, tol_newton, tol_cg, tol_mg
!!$         write(*,*) 'pbe params', pbe_params%n, pbe_params%temp, pbe_params%lam, &
!!$              pbe_params%c(1:pbe_params%n), pbe_params%q(1:pbe_params%n), &
!!$              pbe_params%fixed_nions, &
!!$              pbe_params%use_steric, &
!!$              pbe_params%linearised, &
!!$              pbe_params%use_fas
         write(*,*)
      endif

      pbe => pbe_params

      allocate(pbe%x(pbe%n))
      select case(neutr_method)
      case('none')
         pbe%neutr_method = dl_mg_neutralise_none
      case('jellium_unif')
         pbe%neutr_method = dl_mg_neutralise_with_jellium_uniform
      case('jellium_vacc')
         pbe%neutr_method = dl_mg_neutralise_with_jellium_vacc
      case('ions_fixed')
         pbe%neutr_method = dl_mg_neutralise_with_ions_fixed
         pbe%x = x_ratios
      case('ions_auto')
         pbe%neutr_method = dl_mg_neutralise_with_ions_auto
      case('ions_auto_lin')
         pbe%neutr_method = dl_mg_neutralise_with_ions_auto_linear
      case default
         call abort('unknown neutralisation method')
      end select




   ! set start and end indices, 1d decomposition

      call set_grid_index

      idx_s = (/ sx, sy, sz /)
      idx_e = (/ ex, ey, ez /)

      if ( sum(my_coords) <= 2) then
        write(0,*) 'idx_se ', idx_s, idx_e
      endif

      ngxyz = [ngx, ngy, ngz]

      if ( myid == 0) then
         write(*,*) "grid ", ngxyz
      endif

      allocate(pot(sx:ex, sy:ey, sz:ez), &
             res(sx:ex, sy:ey, sz:ez), &
             rho(sx:ex, sy:ey, sz:ez),stat=ierr)
      if (ierr /= 0)then
         call abort('error at array allocation')
      endif

      ! initialise mg_utils module data
      call mg_utils_init(mg_comm,ngx,ngy,ngz,idx_s)

      alf = - fourpi ! rho prefactor in pbe_test

      ! some defaults are adjusted below

      select case(model_name)
      case("laplace")
         call set_bc([.false., .false., .false.], [1.0_wp, 1.0_wp, 1.0_wp], &
              ngxyz, boundary_condition, dx, dy,dz)

      case("parabolic_rho")
         xl = 1.0_wp ; yl = xl; zl = xl !104.25_wp
         call set_bc([.false., .false., .false.], [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy,dz)
         pert_amplitude = 0.01_wp

      case("sin_eps")
         xl = 1.0_wp ; yl = xl; zl = xl !104.25_wp
         call set_bc([.false., .false., .false.], [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         ! note: because of integer division the period of sin functions is not
         ! an integer number of lattice steps
         eps_sin_params = eps_sin_params_t((ngx-1)/8, (ngy-1)/8, (ngz-1)/8, 0.1_wp)
      case("linear_phi")
         xl = 1.0_wp ; yl = xl; zl = xl !104.25_wp
         call set_bc([.false., .false., .false.], [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         alpha(:) = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
         beta(:)  = (/ 0.0_wp, 1.0_wp, 1.0_wp /)
         if (any(periodic)) then
            call abort('linear_phi test works with Dirichlet BC')
         endif

      case("z_to_power")
         ! Using Dirichlet BC in z direction
         xl = 1.0_wp ; yl = xl; zl = xl !104.25_wp
         call set_bc([.true., .true., .false.], [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         potpow = 3

      case("x_to_power_pbc", "y_to_power_pbc","z_to_power_pbc")
         xl = 1.0_wp ; yl = xl; zl = xl !104.25_wp
         call set_bc([.true., .true., .true.], [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         potpow = 3

      case("pbc")
         call set_bc([.true., .true., .true.], [1.0_wp, 1.0_wp, 1.0_wp], &
              ngxyz, boundary_condition, dx, dy, dz)
         a11 = 2.0_wp * pi; a21 = 2.0_wp * pi; a31 = 2.0_wp * pi; a32 = 0.0_wp

      case("pbe_lin_pbc")
         call set_bc([.true., .true., .true.], [1.0_wp, 1.0_wp, 1.0_wp], &
              ngxyz, boundary_condition, dx, dy, dz)
         inv_debye2 = -pbe_params%lam*elcharge**2 &
              /(rbohr * pbe_params%temp*kB) &
              *sum( pbe_params%c * pbe_params%q**2)
         a11 = 2.0_wp * pi; a21 = 2.0_wp * pi
         a31 = 2.0_wp * pi; a32 = 0.0_wp

      case("pbc_eps")
         xl =1.0_wp; yl=1.0_wp; zl=1.0_wp
         call set_bc([.true., .true., .true.], [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         a11 = 2.0_wp * pi; a21 = 2.0_wp * pi; a31 = 6.0_wp * pi
         eps_sin_params = eps_sin_params_t(4, 4, 4, 0.1_wp)

      case("pbcxy")
         call set_bc([.true., .true., .false.], [1.0_wp, 1.0_wp, 1.0_wp], &
              ngxyz, boundary_condition, dx, dy, dz)
      case("sin_sum")
         xl = 1.0_wp ; yl = xl; zl = xl !104.25_wp
         call set_bc([.false., .false., .false.], [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)

      case("erf_eps")
         ! model from Fisicario et al paper
         ! Gaussian potential, sigma = 0.5
         ! eps = 1 + (eps0 -1) 1/2 (1 + erf ((r -d0)/Delta))
         ! eps0 = 78.36 , d0 = 1.7, Delta= 0.3
         xl = 10.0_wp; yl = 10.0_wp; zl = 10.0_wp
         call set_bc(periodic, [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         ge%sig = 0.5_wp; ge%mu = 5.0_wp; ge%eps0 = 78.36_wp
         ge%d = 1.7_wp; ge%delta = 0.3_wp; ge%scale=1.0_wp

      case("pbe_erf")
         xl = 30.0_wp; yl = 30.0_wp; zl = 30.0_wp
         call set_bc(periodic, [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         ge%sig = 0.5_wp; ge%mu = xl/2; ge%eps0 = 78.36_wp
         ge%d = 1.7_wp; ge%delta = 0.3_wp; ge%scale = 0.01_wp
         stw_params%d=4.0_wp; stw_params%delta = 0.3_wp
         allocate(vsteric(sx:ex,sy:ey,sz:ez)) ! set to 1 if not in use

      case("pbe_nosteric")
         ! use sin_eps for linear part
         xl = 1.0_wp ; yl = xl; zl = xl !104.25_wp
         call set_bc(periodic, [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         eps_sin_params = eps_sin_params_t(4, 4, 4, 0.1_wp)
         pert_amplitude = 0.1_wp

      case("pbez")
         ! PBC in xy
         xl = 10.0_wp ; yl = xl; zl = xl
         call set_bc([.true., .true., .false.], [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         ! steric potential need not be used with this test
         surface_pot = 200 * 1.e-3_wp *statv  ! mV -> statV
         inv_debye = sqrt(-pbe_params%lam*(elcharge**2/(epsh2o*rbohr))&
                /(pbe_params%temp*kB)*sum( pbe_params%c * pbe_params%q**2))
         if (pbe_params%use_steric) then
            call abort('steric cannot be used with pbez test')
         endif

      case("pbe_point_charge")
         inv_debye = sqrt(-pbe_params%lam*(elcharge**2/(epsh2o*rbohr))&
              /(pbe_params%temp*kB)*sum( pbe_params%c * pbe_params%q**2))
         xl = 10.0_wp ; yl = xl; zl = xl !104.25_wp
         call set_bc([.false., .false., .false.], [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         pointsource(:) = (/ (ngx/2+0.5_wp)*dx, (ngy/2+0.5_wp)*dy, (ngz/2+0.5_wp)*dz /)
         if ( pbe_params%use_steric) then
            allocate(vsteric(sx:ex,sy:ey,sz:ez))
         endif

      case("pbe_sphere")
         ! charge distribution placed at centre of the box
         ! charge density 0.75 * rmax - r
         ! i.e. neutral outside rmax
         ! steric potential astc * exp( -bstc * r)
         xl = 1.0; yl = 1.0; zl = 1.0
         call set_bc(periodic, [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         if (pbe_params%use_steric) allocate(vsteric(sx:ex,sy:ey,sz:ez))
         astc = 1.0_wp; bstc = 1.0_wp/ 0.25_wp

      case("pbe_sphere_hardcore_steric")
         ! uses pbe_sphere but with hard core steric
         xl = 6.0_wp ; yl = xl; zl = xl !104.25_wp
         call set_bc(periodic, [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         allocate(vsteric(sx:ex,sy:ey,sz:ez))

      case("pbe_molecule")
         xl = 20.0_wp ; yl = xl; zl = xl !104.25_wp
         call set_bc([.false., .false., .false.], [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         if (pbe_params%use_steric) allocate(vsteric(sx:ex,sy:ey,sz:ez))

      case("poisson_point_charge")
         xl = 1.0_wp ; yl = xl; zl = xl !104.25_wp
         call set_bc([.false., .false., .false.], [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         pointsource(:) = (/ -0.1_wp, 0.5_wp, 0.5_wp /)

      case ("wirez")
         xl = 1.0_wp ; yl = xl; zl = xl !104.25_wp
         call set_bc([.false., .false., .true.], [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         pointsource(:) = (/ (ngx/2+0.5_wp)*dx, (ngy/2+0.5_wp)*dy, 0.0_wp /)

      case("wirex")
         xl = 1.0_wp ; yl = xl; zl = xl !104.25_wp
         call set_bc([.true., .false., .false.], [xl, yl, zl], &
              ngxyz, boundary_condition, dx, dy, dz)
         pointsource(:) = (/ 0.0_wp, (ngy/2+0.5_wp)*dy, (ngz/2+0.5_wp)*dz /)

      case("onetep_t26")
! this test will work only if the ONETEP data are provided, disregarded in devel branch
         write(*,*) "ontep data"
         boundary_condition = DL_MG_BC_DIRICHLET
         call read_onetep_data

      case('jacek_pbe_pbc') !misnommer, just a tag for Jacek cases
         dx = 8.766490E-02 / 0.529177210903_wp;
         dy = dx;
         dz = dx ! A -> a0
         boundary_condition = DL_MG_BC_DIRICHLET
         !where (periodic)
         !   boundary_condition = DL_MG_BC_PERIODIC
         !end where
         if ( pbe_params%use_steric) then
            allocate(vsteric(sx:ex,sy:ey,sz:ez))
         endif

     case default
        call abort('unknown model in poisson_test :'//trim(model_name))
     end select

     call build_initial_pot
     !call write_vect(pot,lbound(pot),(/imin,jmin,kmin/),(/imax,jmax,kmax/),'initial_pot')
     call build_rho
     !call write_vect(rho,lbound(rho),(/imin,jmin,kmin/),(/imax,jmax,kmax/),'rho')

     ! test for the version string
     if (myid == 0) then
!!$        call dl_mg_version(version_string_short)
!!$        write (*,*) "DL_MG version: ", version_string_short
       call dl_mg_version(version_string)
       write (*,*) "DL_MG version: ", version_string
     endif

     !write(0,*) 'mg_test ', boundary_condition
     call dl_mg_init(ngx, ngy, ngz, &
                     dx, dy, dz, &
                     boundary_condition, &
                     idx_s, idx_e, &
                     mg_comm, &
                     98, "poisson_"//adjustl(model_name),&
                     min_thread_chunk=0, block_size=0, &
                     errors_return = errors_return, ierror=ierr)

     if (ierr /= DL_MG_SUCCESS) then
        call dl_mg_error_string(ierr, err_msg)
        write(0,*) "error in dl_mg_init: "//err_msg
        call abort("quiting...")
     endif

      select case (model_name)
      case ("pbe_nosteric","pbe_sphere", "pbe_sphere_hardcore_steric", &
           "pbe_molecule", "pbez", "jacek_pbe_pbc","pbe_lin_pbc","pbc_eps", &
           "pbe_erf")
         call pbe_test
      case ("pbe_point_charge")
         call pbe_point_charge_test
      case default
         call build_eps_half

         call dl_mg_solver(eps_full, eps_half, 1.0_wp,  rho,  pot, fd_order, &
              tol_res_rel = tol_res_rel, tol_res_abs=tol_res_abs, &
              tol_pot_rel= tol_pot_rel, tol_pot_abs=tol_pot_abs, &
              tol_cg_res_rel=tol_cg, &
              tol_vcycle_res_rel=tol_mg, use_damping=use_damping,&
              use_cg=use_cg,&
              ierror=ierr)
         if (ierr /= DL_MG_SUCCESS) then
            write(0,*) "test: error code ", ierr
            call dl_mg_error_string(ierr, err_msg)
            write(0,*) err_msg
            call abort("quitting ...")
         endif
      end select

      call dl_mg_free

      !call write_slice("xy",ngz/4, sx, sy, sz, pot, model_name)
      !call write_slice("xy",ngz/2, sx, sy, sz, pot, model_name)
      !write(22,*) pot

      call test_solution

      !call dump_pbe1d("xz", ngy/2)

      !call write_vect(pot,lbound(pot),(/imin,jmin,kmin/),(/imax,jmax,kmax/),'poisson_last_pot')
      !call write_vect(residual,lbound(residual),(/imin,jmin,kmin/),(/imax,jmax,kmax/),'poisson_last_res')

    contains

      subroutine set_grid_index
        implicit none
        integer mlz, ir, ierr
        integer zpart(0:pdims(3)-1)



#ifdef MPI

        sx = local_shift(1) + 1
        ex = local_shift(1) + nlx
        sy = local_shift(2) + 1
        ey = local_shift(2) + nly
        sz = local_shift(3) + 1
        ez = local_shift(3) + nlz

!!$        mlz = (ngz-1)/pdims(3)
!!$        ir = mod(ngz-1,pdims(3))
!!$
!!$        if (my_coords(3) < ir) mlz = mlz+1
!!$
!!$        call mpi_allreduce(mlz,ir,1,mpi_integer,MPI_SUM,mg_comm,ierr)
!!$        if ( ir /= ngz - 2) then
!!$           write(0,*) 'mg test: set grid:: problem with partition', my_coords(3),mlz
!!$        end if
!!$        call mpi_allgather(mlz,1,mpi_integer,zpart,1,mpi_integer,mg_comm,ierr)
!!$
!!$      sz = sum(zpart(0:my_coords(3)-1)) + 1 + 1
!!$      ez = sz + mlz - 1
!!$
!!$      if ( my_coords(3) == 0 ) sz = sz - 1
!!$      if ( my_coords(3) == pdims(3) - 1 ) ez = ez + 1
#else
      sx = 1; ex = ngx
      sy = 1; ey = ngy
      sz = 1; ez = ngz
#endif

    end subroutine set_grid_index


    subroutine build_eps_half
      use dl_mg_mpi_header
      implicit none

      integer i,j,k,d,ierr
      real(wp) x,y,z, dr(3,3), rmax, x0,y0,z0, r
      real(wp), allocatable :: aux(:,:,:)

      ! onetep data test allocate their arrays separately
      if ( .not. allocated(eps_half))then
         allocate(eps_half(sx:ex, sy:ey, sz:ez,3))
      endif
      if ( .not. allocated(eps_full))then
         allocate(eps_full(sx:ex, sy:ey, sz:ez))
      endif

      select case (model_name)
      case("laplace","parabolic_rho", "linear_phi", "z_to_power", "z_to_power_pbc", &
           "y_to_power_pbc", "x_to_power_pbc", "pbc",&
           "poisson_point_charge",  &
           "pbe_molecule","pbcxy","wirez", "wirex", &
           "pbe_lin_pbc")
         eps_half = 1.0_wp
         eps_full = 1.0_wp
      case ("pbe_point_charge", "pbez")
         eps_half = epsh2o
         eps_full = epsh2o

      case("sin_eps","sin_sum","pbe_nosteric", "pbc_eps")
         do k = sz, ez
            do j = sy, ey
               do i = sx, ex
                  x=(i-1.0_wp)*dx
                  y=(j-1.0_wp)*dy
                  z=(k-1.0_wp)*dz
                  eps_half(i,j,k,1) = eps_sin(eps_sin_params,(x+0.5_wp)/xl,y/yl,z/zl)
                  eps_half(i,j,k,2) = eps_sin(eps_sin_params,x/xl,(y+0.5_wp)/yl,z/zl)
                  eps_half(i,j,k,3) = eps_sin(eps_sin_params,x/xl,y/yl,(z+0.5_wp)/zl)
                  eps_full(i,j,k)   = eps_sin(eps_sin_params,x/xl,y/yl,z/zl)
               enddo
            enddo
         enddo
      case("erf_eps","pbe_erf")
         do k = sz, ez
            do j = sy, ey
               do i = sx, ex
                  x=(i-1.0_wp)*dx
                  y=(j-1.0_wp)*dy
                  z=(k-1.0_wp)*dz
                  eps_half(i,j,k,1) = erf_eps(ge%eps0,ge%d,ge%delta,ge%mu,&
                       x+0.5_wp*dx,y,z)
                  eps_half(i,j,k,2) = erf_eps(ge%eps0,ge%d,ge%delta,ge%mu,&
                       x,y+0.5_wp*dy,z)
                  eps_half(i,j,k,3) = erf_eps(ge%eps0,ge%d,ge%delta,ge%mu,&
                       x,y,z+0.5_wp*dz)
                  eps_full(i,j,k)   = erf_eps(ge%eps0,ge%d,ge%delta,ge%mu,&
                       x,y,z)
               enddo
            enddo
         enddo
      case ("pbe_sphere", "pbe_sphere_hardcore_steric")
         x0=0.5_wp *xl; y0=0.5_wp *yl; z0=0.5_wp *zl
         if (model_name == "pbe_sphere") then
           rmax=0.25_wp;
         else
           rmax=1.3_wp
         end if
         dr = reshape([ 0.5_wp*dx, 0.0_wp, 0.0_wp, &
              0.0_wp, 0.5_wp*dy, 0.0_wp, &
              0.0_wp, 0.0_wp, 0.5_wp*dz], [3,3])

         do d = 1,3
            do k = sz, ez
               do j = sy, ey
                  do i = sx, ex
                     x=(i-1.0_wp)*dx
                     y=(j-1.0_wp)*dy
                     z=(k-1.0_wp)*dz
                     r = sqrt((x+ dr(1,d)-x0)**2+(y+dr(2,d)-y0)**2+(z+dr(3,d)-z0)**2)
                     if ( r < rmax) then
                        eps_half(i,j,k,d) = 1.0_wp + (epsh2o-1.0_wp)*(r/rmax)**2 * ( 2.0_wp-(r/rmax))**2
                     else
                        eps_half(i,j,k,d) = epsh2o
                     end if
                  end do
               enddo
            enddo
         enddo

         do k = sz, ez
             do j = sy, ey
                do i = sx, ex
                   x=(i-1.0_wp)*dx
                   y=(j-1.0_wp)*dy
                   z=(k-1.0_wp)*dz
                   r = sqrt((x -x0)**2+(y-y0)**2+(z-z0)**2)
                   if ( r < rmax) then
                      eps_full(i,j,k) = 1.0_wp + (epsh2o-1.0_wp)*(r/rmax)**2 !* ( 2.0_wp-(r/rmax))**2
                   else
                      eps_full(i,j,k) = epsh2o
                   end if
                end do
             enddo
          enddo


      case ("jacek_pbe_pbc")
         allocate(aux(166,166,166))
         open(33, file='epsx.raw', status='old', action='read')
         read(33, *) aux
         eps_half(sx:ex,sy:ey,sz:ez,1) =aux(sx:ex,sy:ey,sz:ez)
         close(33)
         open(33, file='epsy.raw', status='old', action='read')
         read(33, *) aux
         eps_half(sx:ex,sy:ey,sz:ez,2)=aux(sx:ex,sy:ey,sz:ez)
         close(33)
         open(33, file='epsz.raw', status='old', action='read')
         read(33, *) aux
         eps_half(sx:ex,sy:ey,sz:ez,3)=aux(sx:ex,sy:ey,sz:ez)
         close(33)
         open(33, file='epsf.raw', status='old', action='read')
         read(33, *) aux
         eps_full(sx:ex,sy:ey,sz:ez)=aux(sx:ex,sy:ey,sz:ez)
         close(33)

!!$         do k = 2, 184
!!$            do j = 2, 184
!!$               do i = 2, 184
!!$                  eps_full(i,j,k) = 1.0_wp/6_wp * &
!!$                  ( eps_half(i-1,j,k,1) + eps_half(i,j,k,1) &
!!$                  + eps_half(i,j-1,k,2) + eps_half(i,j,k,2) &
!!$                  + eps_half(i,j,k-1,3) + eps_half(i,j,k,3))
!!$               end do
!!$            end do
!!$         end do
!!$         eps_full(1,:,:) = eps_half(1,:,:,1)

      case default
         call abort("build_eps_half: unknown model, quitting ....")
      end select

    end subroutine build_eps_half


    subroutine build_initial_pot
      use dl_mg_mpi_header
      implicit none

      integer i,j,k
      real(wp) norm_eigenfun,x,y,z,w
      real(wp), allocatable :: aux(:,:,:)

      select case (model_name)
      case("laplace")
         ! set boundary condition
         pot(sx,  :,  :) = 0.0_wp
         pot(ex,  :,  :) = 0.0_wp
         pot( :, sy,  :) = 0.0_wp
         pot( :, ey,  :) = 0.0_wp
         pot( :,  :, sz) = 0.0_wp
         pot( :,  :, ez) = 0.0_wp

         ! start with an eigenvector
         ! eigenvalue index ( imported form mg_utils)
         kx = 1; ky = 1; kz = 1


         do k = sz+1, ez-1
            do j = sy+1, ey-1
               do i = sx+1, ex-1
                  pot(i,j,k) = eigenfun_rb3d(i,j,k)
                  !pot(i,j,k) = mod(i+j+k+1,2)
               enddo
            enddo
         enddo

         !write(0,*) ' u ', u(1:nlx,1:nly,1:1)
#ifdef MPI
         call mpi_allreduce(sum(pot(sx:ex,sy:ey,sz:ez)**2),norm_eigenfun,&
              1,MPI_DOUBLE_PRECISION,MPI_SUM,mg_comm,ierr)
#else
         norm_eigenfun = sum(pot(sx:ex, sy:ey, sz:ez)**2)
#endif

         norm_eigenfun = sqrt(dx*dy*dz*norm_eigenfun)
         if (myid == 0) then
            print*, ' laplace initial norm ',  norm_eigenfun
         endif
         pot(1:nlx,1:nly,1:nlz) = pot(1:nlx,1:nly,1:nlz)/norm_eigenfun

      case("parabolic_rho","sin_eps","pbe_nosteric", "pbe_sphere",&
           "erf_eps", "pbe_erf", "pbe_sphere_hardcore_steric",&
           "pbe_molecule", "pbe_point_charge")
         ! set boundary condition
         pot(sx, :, :) = 0.d0
         pot(ex, :, :) = 0.d0
         pot(:, sy, :) = 0.d0
         pot(:, ey, :) = 0.d0
         pot(:, :, sz) = 0.d0
         pot(:, :, ez) = 0.d0

         w = 1.0_wp ! 16.0_wp
         do k=sz,ez
            do j=sy,ey
               do i=sx,ex
                  x=w*(i)*dx
                  y=w*(j)*dy
                  z=w*(k)*dz
                  pot(i,j,k) = 0.0_wp !0.5_wp * x * (1.0_wp-x) * &
                  !y * (1.0_wp-y)* z * (1.0_wp-z) + &
                  !pot_perturbation(pert_amplitude,ngx/2,ngy/2,ngz/2,x,y,z)
                  !pot(i,j,k)= sin(pi*x) * sin(pi*y) * &
                  !            sin(pi*z)/(3.0_wp*pi**2)
                  !pot(i,j,k) = mod(i+is+j+js+k+ks+1,2)
               enddo
            enddo
         enddo
      case("linear_phi")
         ! set boundary condition

         ! faces perpedicular on x direction
         if ( my_coords(1) == 0)  then
            x = 0.0_wp
            do k=sz,ez
               z=(k-1)*dz
               do j=sy,ey
                  y = (j-1)*dy
                  pot(sx, j, k) = (alpha(1)*x+beta(1))*(alpha(2)*y+beta(2))*(alpha(3)*z+beta(3))
               enddo
            enddo
         endif

         if ( my_coords(1) == pdims(1) - 1)  then
            x = xl
            do k=sz,ez
               z=(k-1)*dz
               do j=sy,ey
                  y = (j-1)*dy
                  pot(ex, j, k) = (alpha(1)*x+beta(1))*(alpha(2)*y+beta(2))*(alpha(3)*z+beta(3))
               enddo
            enddo
         endif

         ! faces perpendicular on y direction
         if ( my_coords(2) == 0 )  then
            y=0.0_wp
            do k=sz,ez
               z=(k-1)*dz
               do i=sx,ex
                  x = (i-1)*dx
                  pot(i, sy, k) = (alpha(1)*x+beta(1))*(alpha(2)*y+beta(2))*(alpha(3)*z+beta(3))
               enddo
            enddo
         endif

         if ( my_coords(2) == pdims(2) - 1)  then
            y = yl
            do k=sz,ez
               z=(k-1)*dz
               do i=sx,ex
                  x = (i-1)*dx
                  pot(i, ey, k) = (alpha(1)*x+beta(1))*(alpha(2)*y+beta(2))*(alpha(3)*z+beta(3))
               enddo
            enddo
         endif


         ! faces perpendicular on z direction
         if ( my_coords(3) == 0)  then
            z = 0.0_wp
            do j=sy,ey
               y=(j-1)*dy
               do i=sx,ex
                  x = (i-1)*dx
                  pot(i, j, sz) = (alpha(1)*x+beta(1))*(alpha(2)*y+beta(2))*(alpha(3)*z+beta(3))
               enddo
            enddo
         endif

         if ( my_coords(3) == pdims(3) - 1) then
            z = zl
            do j=sy,ey
               y=(j-1)*dy
               do i=sx,ex
                  x = (i-1)*dx
                  pot(i, j, ez) = (alpha(1)*x+beta(1))*(alpha(2)*y+beta(2))*(alpha(3)*z+beta(3))
               enddo
            enddo
         endif

         do k=sz+1,ez-1
            do j=sy+1,ey-1
               do i=sx+1,ex-1
                  pot(i,j,k) = 0.0_wp
               enddo
            enddo
         enddo

      case("z_to_power")
         ! set boundary condition
         ! bc -> Per Per Dir

         ! faces perpendicular on z direction
         if ( my_coords(3) == 0)  then
            z = 0.0_wp
            do j=sy,ey
               y=(j-1)*dy
               do i=sx,ex
                  x = (i-1)*dx
                  pot(i, j, sz) = 0.0_wp
               enddo
            enddo
         endif

         if ( my_coords(3) == pdims(3) - 1) then
            z = zl
            do j=sy,ey
               y=(j-1)*dy
               do i=sx,ex
                  x = (i-1)*dx
                  pot(i, j, ez) = 0.0_wp !z**potpow
               enddo
            enddo
         endif

         do k=sz+1,ez-1
            do j=sy+1,ey-1
               do i=sx+1,ex-1
                  pot(i,j,k) = 0.0_wp
               enddo
            enddo
         enddo

      case("poisson_point_charge")
         call set_bc_linear_pbe(pbe_params, 0.0_wp, 1.0_wp)
      case("pbc", "pbc_eps", "pbe_lin_pbc", &
           "x_to_power_pbc","y_to_power_pbc","z_to_power_pbc")
         pot(:,:,:) = 0.0_wp
      case("pbcxy")
         if (my_coords(3) == 0) then
            pot (:, :, sz) = 1.0_wp
         endif
         if (my_coords(3) == pdims(3) -1) then
            pot (:, :, ez) = 0.0_wp
         endif
      case ("pbez")
         if (my_coords(3) == 0) then
            ! x mV / (e/rB)
            pot (:, :, sz) = surface_pot &
                 / (elcharge /rbohr)
         endif
         if (my_coords(3) == pdims(3) -1) then
            pot (:, :, ez) = pbe1d(zl, surface_pot, pbe_params%temp, &
                 inv_debye)
         endif
      case("wirez")
         call set_bc_wire("z")
      case("wirex")
         call set_bc_wire("x")

      case("jacek_pbe_pbc")
         !pot = 0.0_wp
         allocate(aux(166,166,166))
         open(33,file='pot.raw', status='old', action='read')
         read(33,*) aux
         rho(sx:ex,sy:ey,sz:ez) = aux(sx:ex,sy:ey,sz:ez)
         close(33)

      case default
         write(0,*) "build_initial_pot: unknown model, quitting ...."
#ifdef MPI
         call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
#else
         stop
#endif
      end select

    end subroutine build_initial_pot


    subroutine build_rho
      use dl_mg_mpi_header
      implicit none

      integer i,j,k,ierr
      real(wp) x,y,z,w, r, r2, rmax, x0, y0, z0, q
      real(wp), allocatable :: aux(:,:,:)

      !write(0,*) ' rho  ', myid, is,js,ks,ngx,ngy,ngz

      select case(model_name)
      case("laplace")
         rho = 0.0_wp
      case("parabolic_rho")
         w = 1.0_wp !16.0_wp
         do k=sz,ez
            do j=sy,ey
               do i=sx,ex
                  x=w*(i-1.0_wp)*dx
                  y=w*(j-1.0_wp)*dy
                  z=w*(k-1.0_wp)*dz
                  rho(i,j,k) = y*(1.0_wp-y)*z*(1.0_wp-z)+ &
                       x*(1.0_wp-x)*z*(1.0_wp-z)+ &
                       x*(1.0_wp-x)*y*(1.0_wp-y)
                  !rho(i,j,k) = -n_neighbor(i+is,j+js,k+ks) * mod(i+is+j+js+k+ks,2) !
                  !rho(i,j,k) = sin(pi*x) * sin(pi*y) * sin(pi*z)
               enddo
            enddo
         enddo
      case("sin_eps","pbe_nosteric")
         do k=sz,ez
            do j=sy,ey
               do i=sx,ex
                  x=(i-1.0_wp)*dx
                  y=(j-1.0_wp)*dy
                  z=(k-1.0_wp)*dz
                  rho(i,j,k) = eps_sin(eps_sin_params,x,y,z)*y*(1.0_wp-y)*z*(1.0_wp-z)/xl**2+ &
                       dx_eps_sin(eps_sin_params,x,y,z)*(x-0.5_wp)*y*(1.0_wp-y)*z*(1.0_wp-z)/xl**2 + &
                       eps_sin(eps_sin_params,x,y,z)*x*(1.0_wp-x)*z*(1.0_wp-z)/yl**2+ &
                       dy_eps_sin(eps_sin_params,x,y,z)*(y-0.5_wp)*x*(1.0_wp-x)*z*(1.0_wp-z)/yl**2 + &
                       eps_sin(eps_sin_params,x,y,z)*x*(1.0_wp-x)*y*(1.0_wp-y)/zl**2 + &
                       dz_eps_sin(eps_sin_params,x,y,z)*(z-0.5_wp)*x*(1.0_wp-x)*y*(1.0_wp-y)/zl**2
               enddo
            enddo
         enddo
      case("erf_eps")
         do k=sz,ez
            do j=sy,ey
               do i=sx,ex
                  x=(i-1.0_wp)*dx
                  y=(j-1.0_wp)*dy
                  z=(k-1.0_wp)*dz
                  r2 = (x-ge%mu)**2 + (y-ge%mu)**2+(z-ge%mu)**2
                  r  = sqrt(r2)
                  rho(i,j,k) = gauss(ge%sig,ge%mu, ge%scale, x,y,z) &
                       / ge%sig**2 &
                       * (erf_eps(ge%eps0, ge%d, ge%delta, ge%mu, x, y, z) &
                       * (r2/ge%sig**2-3.0) &
                       - (ge%eps0 -1)*r /(ge%delta * sqrt(pi)) &
                       * exp(-((r-ge%d)/ge%delta)**2))
               enddo
            enddo
         enddo

      case("pbe_erf")

         call rho_pbe_erf(rho)

      case("linear_phi","poisson_point_charge", "pbcxy", "pbez")
         rho(:,:,:) = 0.0_wp
      case("z_to_power","z_to_power_pbc")
         do k=sz,ez
            do j=sy,ey
               do i=sx,ex
                  x=(i-1.0_wp)*dx
                  y=(j-1.0_wp)*dy
                  z=(k-1.0_wp)*dz
                  rho(i,j,k) = potpow * (potpow -1.0_wp) * z**(potpow -2) * (1.0_wp -z)**potpow &
                       -2 * potpow**2 * z**(potpow-1) * (1.0_wp -z)**(potpow-1) &
                       + potpow * (potpow -1.0_wp) * z**(potpow) * (1.0_wp -z)**(potpow-2)
               enddo
            enddo
         enddo
      case("y_to_power")
         do k=sz,ez
            do j=sy,ey
               do i=sx,ex
                  x=(i-1.0_wp)*dx
                  y=(j-1.0_wp)*dy
                  z=(k-1.0_wp)*dz
                  rho(i,j,k) = potpow * (potpow -1.0_wp) * y**(potpow -2) * (1.0_wp -y)**potpow &
                       -2 * potpow**2 * y**(potpow-1) * (1.0_wp -y)**(potpow-1) &
                       + potpow * (potpow -1.0_wp) * y**(potpow) * (1.0_wp -y)**(potpow-2)
               enddo
            enddo
         enddo

      case("x_to_power")
         do k=sz,ez
            do j=sy,ey
               do i=sx,ex
                  x=(i-1.0_wp)*dx
                  y=(j-1.0_wp)*dy
                  z=(k-1.0_wp)*dz
                  rho(i,j,k) = potpow * (potpow -1.0_wp) * x**(potpow -2) * (1.0_wp -x)**potpow &
                       -2 * potpow**2 * x**(potpow-1) * (1.0_wp -x)**(potpow-1) &
                       + potpow * (potpow -1.0_wp) * x**(potpow) * (1.0_wp -x)**(potpow-2)
               enddo
            enddo
         enddo
      case("pbe_point_charge")
         q = 1.0_wp/(8.0_wp*dx*dy*dz)
         rho(:,:,:) = 0.0_wp
         do k=ngz/2, ngz/2+1
            do j=ngy/2, ngy/2+1
               do i=ngx/2, ngx/2+1
                  if( sx <= i .and. i <= ex .and. &
                       sy <= j .and. j <= ey .and. &
                       sz <= k .and. k <= ez ) then
                     rho(i,j,k) = q
                  endif
               enddo
            enddo
         enddo
         rmax=1.3_wp;
         x0=pointsource(1); y0=pointsource(2); z0=pointsource(3)
         if ( pbe_params%use_steric) then
            do k=sz,ez
               do j=sy,ey
                  do i=sx,ex
                     x=(i-1.0_wp)*dx
                     y=(j-1.0_wp)*dy
                     z=(k-1.0_wp)*dz
                     r = sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
                     if ( r < 2.0_wp * rmax ) then
                        vsteric(i,j,k) = 0.0_wp
                     else
                        vsteric(i,j,k) = 1.0_wp
                     end if
                  enddo
               enddo
            enddo
         endif
      case("pbe_sphere")
         x0=0.5_wp; y0=0.5_wp; z0=0.5_wp
         rmax=0.25_wp;
         do k=sz,ez
            do j=sy,ey
               do i=sx,ex
                  x=(i-1.0_wp)*dx
                  y=(j-1.0_wp)*dy
                  z=(k-1.0_wp)*dz
                  r = sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
                  if ( r < rmax) then
                     rho(i,j,k) = (0.75_wp * rmax - r )
                  else
                     rho(i,j,k) = 0.0_wp
                  endif
                  if (pbe_params%use_steric) then
                     vsteric(i,j,k) = exp (- astc*exp( -bstc * r))! &
                     !* hartree / pbe_params%temp )
                  end if
               enddo
            enddo
         enddo

      case("pbe_sphere_hardcore_steric")
         x0=0.5_wp *xl; y0=0.5_wp * yl; z0=0.5_wp *zl
         rmax=1.3_wp;
         do k=sz,ez
            do j=sy,ey
               do i=sx,ex
                  x=(i-1.0_wp)*dx
                  y=(j-1.0_wp)*dy
                  z=(k-1.0_wp)*dz
                  r = sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
                  if ( r < rmax) then
                     !rho(i,j,k) = 1.0_wp / rmax**3 * (0.75_wp  - (r/rmax) )
                     ! net charge 1
                     rho(i,j,k) = 1.0_wp / rmax**3 * (15.0_wp * 0.25_wp  - (r/rmax) )/(fourpi)
                  else
                     rho(i,j,k) = 0.0_wp
                  endif
                  if ( r < 1.0_wp * rmax ) then
                     vsteric(i,j,k) = 0.0_wp
                  else
                     vsteric(i,j,k) = 1.0_wp
                  end if
               enddo
            enddo
         enddo
      case("wirez")
         ! 0 diametre wire in along z
         ! i.e. point charge in 2D
         q = 2.0_wp * pi * 1.0_wp/(4.0_wp * dx * dy)
         rho(:,:,:) = 0.0_wp
         do k=sz, ez
            do j=ngy/2, ngy/2+1
               do i=ngx/2, ngx/2+1
                  if( sx <= i .and. i <= ex .and. &
                       sy <= j .and. j <= ey ) then
                     rho(i,j,k) = q
                  endif
               enddo
            enddo
         enddo
      case("wirex")
         ! 0 diametre wire in along z
         ! i.e. point charge in 2D
         q = 2.0_wp * pi * 1.0_wp/(4.0_wp * dy * dz)
         rho(:,:,:) = 0.0_wp
         do k=ngz/2, ngz/2+1
            do j=ngy/2, ngy/2+1
               do i=sx, ex
                  if( sz <= k .and. k <= ez .and. &
                       sy <= j .and. j <= ey ) then
                     rho(i,j,k) = q
                  endif
               enddo
            enddo
         enddo
      case("pbe_molecule")
         call init_pbe_molecule

      case("pbc")
         do k=sz,ez
            do j=sy,ey
               do i=sx,ex
                  x=(i-1.0_wp)*dx
                  y=(j-1.0_wp)*dy
                  z=(k-1.0_wp)*dz
                  rho(i,j,k) = -(a11**2+a21**2+a31**2)*sin(a11*x)*sin(a21*y)&
                       *sin(a31*z)+2.3_wp
               enddo
            enddo
         enddo

         case("pbe_lin_pbc")
         do k=sz,ez
            do j=sy,ey
               do i=sx,ex
                  x=(i-1.0_wp)*dx
                  y=(j-1.0_wp)*dy
                  z=(k-1.0_wp)*dz
                  rho(i,j,k) = (a11**2+a21**2+a31**2+inv_debye2) &
                       *sin(a11*x)*sin(a21*y)*sin(a31*z)/ (4.0_wp *pi)+ 1.0
               enddo
            enddo
         enddo

      case("pbc_eps")
         do k=sz,ez
            do j=sy,ey
               do i=sx,ex
                  x=(i-1.0_wp)*dx
                  y=(j-1.0_wp)*dy
                  z=(k-1.0_wp)*dz
                  rho(i,j,k) = -1.0_wp/(4.0_wp *pi) &
                       * (-eps_sin(eps_sin_params,x,y,z)  &
                       * (a11**2+a21**2+a31**2)*sin(a11*x)*sin(a21*y)*sin(a31*z) &
                       + dx_eps_sin(eps_sin_params,x,y,z) * a11 * cos(a11*x) * sin(a21*y) * sin(a31*z) &
                       + dy_eps_sin(eps_sin_params,x,y,z) * a21 * sin(a11*x) * cos(a21*y) * sin(a31*z) &
                       + dz_eps_sin(eps_sin_params,x,y,z) * a31 * sin(a11*x) * sin(a21*y) * cos(a31*z))
               enddo
            enddo
         enddo

      case("jacek_pbe_pbc")
         allocate(aux(166,166,166))
         open(33,file='rho.raw', status='old', action='read')
         read(33,*) aux
         rho(sx:ex,sy:ey,sz:ez) = aux(sx:ex,sy:ey,sz:ez)
         close(33)
         if (pbe_params%use_steric) then
           open(33,file='gamma.raw', status='old', action='read')
           read(33,*) aux
           vsteric(sx:ex,sy:ey,sz:ez) = aux(sx:ex,sy:ey,sz:ez)
         close(33)
         end if

      case default
         call abort("build_rho: unknown model, quitting ....")
      end select

    end subroutine build_rho


      subroutine pbe_test
        implicit none

        integer np
        integer i, iter, ierror, ierr
        real(wp) c(pbe_params%n), q(pbe_params%n), &
             mu(pbe_params%n), ion_concentration(pbe_params%n), &
             xnwtr(pbe_params%n), &
             lam, temp, db, maxpot, maxpotl, gamma
        real(wp), pointer :: vsteric_pt(:,:,:)
        real(wp), allocatable, target :: dummy(:,:,:)
        logical use_pot
        character ctemp*10, caux(2)*10
        character(len=DL_MG_MAX_ERROR_STRING) :: errmsg


        caux(1)(1:) = "_lin_"
        caux(2)(1:) = "_ful_"
        lam = pbe_params%lam !-4.0_wp * pi
        np = pbe_params%n
        q(:) = pbe_params%q
        c(:) = pbe_params%c

        temp = pbe_params%temp
        use_pot =.false.

        call build_eps_half

        if (myid == 0 .and. pbe_params%lam /= 0.0 .and. &
             sum(pbe_params%c) /= 0.0) then
           write(*,*) 'pbe_test: Debye length ', &
                1.0_wp/ sqrt(-pbe_params%lam*(elcharge**2/(epsh2o*rbohr))&
                /(temp*kB)*sum( pbe_params%c * pbe_params%q**2))
        endif

        if ( pbe_params%use_steric) then
           vsteric_pt => vsteric
        else
           allocate(dummy(1:0,1:0,1:0)) ! i.e. no steric weight
           vsteric_pt => dummy
           write(0,*) 'dummy', size(dummy)
        end if

        ! apps could call the solver more than once
        do iter=1,1
           call dl_mg_solver(eps_full, eps_half, alf, &
                rho, &
                lam, temp, np, c, q, pot, fd_order, &
                tol_res_rel = tol_res_rel, tol_res_abs = tol_res_abs, &
                tol_pot_rel = tol_pot_rel, tol_pot_abs = tol_pot_abs, &
                tol_newton_res_rel=tol_newton, tol_cg_res_rel=tol_cg, &
                tol_vcycle_res_rel=tol_mg, &
                tol_vcycle_res_abs=1.e-8_wp, &
                mg_max_conv_rate =10.0_wp, &
                use_cg=use_cg,&
                linearised=pbe_params%linearised, &
                steric_weight=vsteric_pt, use_pot_in = use_pot, &
                use_damping=use_damping, use_fas=pbe_params%use_fas, &
                neutralisation_method = pbe_params%neutr_method, &
                neutralisation_ion_ratios = pbe_params%x, &
                betamu_electrostatic = mu, &
                steric_weight_average = gamma, &
                used_ion_concentrations = ion_concentration, &
                used_neutralisation_ratios = xnwtr, &
                ierror=ierror)

           if (ierror /= DL_MG_SUCCESS ) then
              if (myid == 0) then
                 write(*,*) 'PBE test failed to converge:'
                 call dl_mg_error_string(ierror,errmsg)
                 write(0,*) errmsg
              endif
           else
              if (myid == 0) then
                 print*, 'iteration', iter
                 print*, 'mu, gamma', mu, gamma
                 print*, 'returned ion concentration', ion_concentration
                 print*, 'used neutralisation fraction', xnwtr
              end if
           endif
        end do

      end subroutine pbe_test


      subroutine pbe_point_charge_test
        implicit none

        integer np
        integer i, ntemp, ierror, count
        real(wp) c(pbe_params%n), q(pbe_params%n), lam, temp, z, told
        logical usepot
        character(len=32) caux(2), ctemp

        np = pbe_params%n
        lam = pbe_params%lam
        q(:) = pbe_params%q
        c(:) = pbe_params%c
        temp = pbe_params%temp

        call build_eps_half
        ! Boundary conditions must be set here becase we need the PBEparameters

        usepot = .false.; z=0.95; count = 0; told = temp
        do while (temp >= told) ! 1 pass loop
           if (myid == 0) then
              write(*,'(/,a, E10.4,/,a,E10.4,/)') "point charge pbe, Debye lenght: ", &
                   1.0_wp / inv_debye, &
                   "                          temp: ", temp
              write(*,*) 'use linear PBE', pbe_params%linearised, &
                   'use steric', pbe_params%use_steric
           end if
           call set_bc_linear_pbe(pbe_params, potscale=1.0/epsh2o) !using BC

           if (pbe_params%use_steric) then
              call dl_mg_solver(eps_full, eps_half, -4.0_wp*pi, rho, &
                   lam, temp, np, c, q, pot, fd_order, &
                   tol_res_rel = tol_res_rel, tol_res_abs = tol_res_abs, &
                   tol_pot_rel = tol_pot_rel, tol_pot_abs = tol_pot_abs, &
                   tol_newton_res_rel=tol_newton, tol_cg_res_rel=tol_cg, &
                   tol_vcycle_res_rel=tol_mg, &
                   use_cg=use_cg, &
                   linearised=pbe_params%linearised, steric_weight=vsteric,&
                   use_pot_in = usepot, use_fas=pbe_params%use_fas, &
                   use_damping=use_damping, ierror=ierror)
           else
              call dl_mg_solver(eps_full, eps_half, -4.0_wp*pi, rho, &
                   lam, temp, np, c, q, pot, fd_order, &
                   tol_res_rel = tol_res_rel, tol_res_abs = tol_res_abs, &
                   tol_pot_rel = tol_pot_rel, tol_pot_abs = tol_pot_abs, &
                   tol_newton_res_rel=tol_newton, tol_cg_res_rel=tol_cg, &
                   tol_vcycle_res_rel=tol_mg, &
                   use_cg=use_cg, &
                   linearised=pbe_params%linearised, &
                   use_pot_in = usepot, use_fas=pbe_params%use_fas, &
                   use_damping=use_damping, ierror=ierror)
           endif

           ! testing linerisation a selected a pot
           !call dl_mg_solver_pbe(eps_half, -4.0_wp*pi, rho, pot, tol = tol, &
           !     linearised=.true., der_pot=pot, use_pot_in = usepot, ierror=ierror)


        !call dl_mg_solver_poisson(eps_half, -4.0_wp*pi, rho, pot, tol = tol, &
        !          ierror=ierror)

           if (ierror /= DL_MG_SUCCESS ) then
              if (myid == 0) then
                 write(*,*) 'linear PBE case failed to converge at beta ', beta, z
              endif
              exit
              !z = 1.0_wp -(1.0_wp-z)*0.5_wp
              !if ( bold -z*bold < 1.0_wp)  exit
              !beta = z * bold
              !if ( beta < 300.0_wp) beta = 300.0_wp
              !pbe_params%beta =  beta !needed to set BC
           else
              write(ctemp,'(E13.5)') temp
              ctemp=adjustl(ctemp)
              if (pbe_params%linearised) then
                 write(caux(1),'(a)') 'linpbe'
              else
                  write(caux(1),'(a)') 'fullpbe'
               end if
              caux(1)=adjustl(caux(1))
              if (pbe_params%use_steric) then
                 write(caux(2),'(a)') 'steric'
              else
                 write(caux(2),'(a)') 'nosteric'
              endif
              call write_slice("xy",ngz/2, sx, sy, sz, pot, trim(model_name)//"_t"//trim(ctemp)//"_"&
                   //trim(caux(1))//"_"//trim(caux(2))//"_")
              if ( temp == 300.0_wp) exit
              told = temp
              temp = z * temp
              if ( temp < 300.0_wp) temp = 300.0_wp
              pbe_params%temp =  temp !needed to set BC
           endif
           count = count+1
           usepot=.true.
        enddo

      end subroutine pbe_point_charge_test


      subroutine set_bc_linear_pbe(pbe_params, q_, potscale)
        implicit none
        type(pbe_params_t), intent(in) :: pbe_params
        real(wp), optional, intent(in) :: q_, potscale
        integer i, j, k
        real(wp) sr(3), r, x, y, z, q, scale

        if ( .not. present(q_)) then
           q = sqrt(- pbe_params%lam * (elcharge**2/(epsh2o*rbohr))/(kB *pbe_params%temp) * &
                sum( pbe_params%c(:) * pbe_params%q(:)**2))
        else
           q = q_
        end if

        if (present(potscale)) then
           scale = potscale
        else
           scale = 1.0_wp ! elcharge / rbohr
        endif

        sr = pointsource!(/ -0.1_wp, 0.5_wp, 0.5_wp /) ! point source coordinates


        ! faces perpedicular on x direction
           if ( my_coords(1) == 0)  then
              x = 0.0_wp
              do k=sz,ez
                 z=(k-1)*dz
                 do j=sy,ey
                    y = (j-1)*dy
                    r = sqrt((x-sr(1))**2+(y-sr(2))**2+(z-sr(3))**2)
                    pot(sx, j, k) = scale * exp(-q*r)/r
                 enddo
              enddo
           endif

           if ( my_coords(1) == pdims(1) - 1)  then
              x = xl
              do k=sz,ez
                 z=(k-1)*dz
                 do j=sy,ey
                    y = (j-1)*dy
                    r = sqrt((x-sr(1))**2+(y-sr(2))**2+(z-sr(3))**2)
                    pot(ex, j, k) = scale * exp(-q*r)/r
                 enddo
              enddo
           endif

           ! faces perpendicular on y direction
           if ( my_coords(2) == 0 )  then
              y=0.0_wp
              do k=sz,ez
                 z=(k-1)*dz
                 do i=sx,ex
                    x = (i-1)*dx
                    r = sqrt((x-sr(1))**2+(y-sr(2))**2+(z-sr(3))**2)
                    pot(i, sy, k) = scale * exp(-q*r)/r
                 enddo
              enddo
           endif

           if ( my_coords(2) == pdims(2) - 1)  then
              y = yl
              do k=sz,ez
                 z=(k-1)*dz
                 do i=sx,ex
                    x = (i-1)*dx
                    r = sqrt((x-sr(1))**2+(y-sr(2))**2+(z-sr(3))**2)
                    pot(i, ey, k) = scale * exp(-q*r)/r
                 enddo
              enddo
           endif

           ! faces perpendicular on z direction
           if ( my_coords(3) == 0)  then
              z = 0.0_wp
              do j=sy,ey
                 y=(j-1)*dy
                 do i=sx,ex
                    x = (i-1)*dx
                    r = sqrt((x-sr(1))**2+(y-sr(2))**2+(z-sr(3))**2)
                    pot(i, j, sz) = scale * exp(-q*r)/r
                 enddo
              enddo
           endif

           if ( my_coords(3) == pdims(3) - 1) then
              z = zl
              do j=sy,ey
                 y=(j-1)*dy
                 do i=sx,ex
                    x = (i-1)*dx
                    r = sqrt((x-sr(1))**2+(y-sr(2))**2+(z-sr(3))**2)
                    pot(i, j, ez) = scale * exp(-q*r)/r
                 enddo
              enddo
           endif

!!$           do k=sz+1,ez-1
!!$              z = (k-1)*dz
!!$              do j=sy+1,ey-1
!!$                 y = (j-1) * dy
!!$                 do i=sx+1,ex-1
!!$                    x = (i-1)*dx
!!$                    !r = sqrt((x-sr(1))**2+(y-sr(2))**2+(z-sr(3))**2)
!!$                    pot(i,j,k) = 0.0_wp !exp(-q*r)/r !0.0_wp
!!$                 enddo
!!$              enddo
!!$           enddo

           !write(*,*) 'sum pot 1:', sum(pot(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1))

         end subroutine set_bc_linear_pbe


         subroutine set_bc_wire(dir)
           implicit none
           character(len=1), intent(in) :: dir
           !type(pbe_params_t), intent(in) :: pbe_params
           !real(wp), optional, intent(in) :: q_, potscale
           integer i, j, k
           real(wp) sr(3), r, x, y, z, q, scale

!!$        if ( .not. present(q_)) then
!!$           q = sqrt(- pbe_params%lam * (elcharge**2/(epsh2o*rbohr))/(kB *pbe_params%temp) * &
!!$                sum( pbe_params%c(:) * pbe_params%q(:)**2))
!!$        else
!!$           q = q_
!!$        end if

!!$        if (present(potscale)) then
!!$           scale = potscale
!!$        else
           scale = 1.0_wp ! elcharge / rbohr
!!$        endif

        sr = pointsource!(/ -0.1_wp, 0.5_wp, 0.5_wp /) ! point source coordinates

        select case(dir)
        case("z")
           ! faces perpedicular on x direction
           if ( my_coords(1) == 0)  then
              x = 0.0_wp
              do k=sz,ez
                 z=(k-1)*dz
                 do j=sy,ey
                    y = (j-1)*dy
                    r = sqrt((x-sr(1))**2+(y-sr(2))**2)
                    pot(sx, j, k) = scale * log(r)
                 enddo
              enddo
           endif

           if ( my_coords(1) == pdims(1) - 1)  then
              x = xl
              do k=sz,ez
                 z=(k-1)*dz
                 do j=sy,ey
                    y = (j-1)*dy
                    r = sqrt((x-sr(1))**2+(y-sr(2))**2)
                    pot(ex, j, k) = scale * log(r)
                 enddo
              enddo
           endif

           ! faces perpendicular on y direction
           if ( my_coords(2) == 0 )  then
              y=0.0_wp
              do k=sz,ez
                 z=(k-1)*dz
                 do i=sx,ex
                    x = (i-1)*dx
                    r = sqrt((x-sr(1))**2+(y-sr(2))**2)
                    pot(i, sy, k) = scale * log(r)
                 enddo
              enddo
           endif

           if ( my_coords(2) == pdims(2) - 1)  then
              y = yl
              do k=sz,ez
                 z=(k-1)*dz
                 do i=sx,ex
                    x = (i-1)*dx
                    r = sqrt((x-sr(1))**2+(y-sr(2))**2)
                    pot(i, ey, k) = scale * log(r)
                 enddo
              enddo
           endif
        case("x")
           if ( my_coords(3) == 0)  then
              z = 0.0_wp
              do j=sy,ey
                 y=(j-1)*dy
                 do i=sx,ex
                    r = sqrt((y-sr(2))**2+(z-sr(3))**2)
                    pot(i, j, sz) = scale * log(r)
                 enddo
              enddo
           endif

           if ( my_coords(3) == pdims(3) - 1)  then
             z = zl
              do j=sy, ey
                 y=(j-1)*dy
                 do i=sx, ex
                    r = sqrt((y-sr(2))**2+(z-sr(3))**2)
                    pot(i, j, ez) = scale * log(r)
                 enddo
              enddo
           endif

           ! faces perpendicular on y direction
           if ( my_coords(2) == 0 )  then
              y=0.0_wp
              do k=sz,ez
                 z=(k-1)*dz
                 do i=sx,ex
                    x = (i-1)*dx
                    r = sqrt((z-sr(3))**2+(y-sr(2))**2)
                    pot(i, sy, k) = scale * log(r)
                 enddo
              enddo
           endif

           if ( my_coords(2) == pdims(2) - 1)  then
              y = yl
              do k=sz,ez
                 z=(k-1)*dz
                 do i=sx,ex
                    x = (i-1)*dx
                    r = sqrt((z-sr(3))**2+(y-sr(2))**2)
                    pot(i, ey, k) = scale * log(r)
                 enddo
              enddo
           endif
        end select

!!$           do k=sz+1,ez-1
!!$              z = (k-1)*dz
!!$              do j=sy+1,ey-1
!!$                 y = (j-1) * dy
!!$                 do i=sx+1,ex-1
!!$                    x = (i-1)*dx
!!$                    !r = sqrt((x-sr(1))**2+(y-sr(2))**2+(z-sr(3))**2)
!!$                    pot(i,j,k) = 0.0_wp !exp(-q*r)/r !0.0_wp
!!$                 enddo
!!$              enddo
!!$           enddo

           !write(*,*) 'sum pot 1:', sum(pot(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1))

         end subroutine set_bc_wire




      subroutine read_onetep_data
        implicit none

        integer i, ns, coords(3), nxi, nyi, nzi, myid, nproc, ierr ! size of the onetep input arrays

        integer, allocatable :: zl(:)
        real(wp), allocatable :: potin(:,:,:), rhoin(:,:,:), eps_halfin(:,:,:,:)

        !ngx=162;ngy=162;ngz=162
        dx=0.166698125
        dy=0.166698125
        dz=0.166698125

        nxi=161; nyi=161; nzi=161

        deallocate(pot,rho,res)

#ifdef MPI
        call mpi_comm_rank(mg_comm, myid, ierr)
#else
        myid = 0
#endif
        if ( myid == 0 ) then

           allocate (potin(nxi,nyi,nzi),rhoin(nxi,nyi,nzi),eps_halfin(nxi,nyi,nzi,3))

           write(*,*) 'before read data'
           open(55,file="onetep_t26/pot.dat",form="UNFORMATTED",status="OLD", access="STREAM")
           read(55) potin
           close(55)
           open(55,file="onetep_t26/rho.dat",form="UNFORMATTED",status="OLD", access="STREAM")
           read(55) rhoin
           close(55)
           open(55,file="onetep_t26/eps_half.dat",form="UNFORMATTED",status="OLD", access="STREAM")
           read(55) eps_halfin
           close(55)
        !write(78) eps_half(1:161,1:161,1:161,1:3)
           rhoin = -4.0_wp*pi*rhoin

        endif

#ifdef MPI
        call mpi_comm_size(mg_comm,nproc,ierr)
#else
        nproc = 1
#endif
        write(0,*) 'did read', nproc

        allocate(zl(0:nproc-1)) !z limits
#ifdef MPI
        call mpi_cart_coords(mg_comm,myid,3,coords,ierr)
#else
        coords = (/ 0, 0, 0/)
#endif
        select case(nproc)
        case (1)
           zl(0) = nzi
           allocate(pot(nxi,nyi,nzi),rho(nxi,nyi,nzi),eps_half(nxi,nyi,nzi,3),res(nxi,nyi,nzi))
        case(2)
           zl = (/ 84, 161 /)
           allocate(pot(nxi,nyi,84),rho(nxi,nyi,84),eps_half(nxi,nyi,84,3),res(nxi,nyi,84))
        case(4)
           zl(0) = 41; zl(1) = 83; zl(2) =125; zl(3) = 161
           allocate(pot(nxi,nyi,42),rho(nxi,nyi,42),eps_half(nxi,nyi,42,3),res(nxi,nyi,42))
        case default
           call abort( "onetep data test: MPI rank patition not provided")
        end select

        if ( myid == 0) then
           if ( coords(3) /= 0) then
              write(0,*) 'wrong coordinate in onetep data, please inspect!'
           endif
           pot(:,:,1:zl(0))=potin(:,:,1:zl(0))
           rho(:,:,1:zl(0))=rhoin(:,:,1:zl(0))
           eps_half(:,:,1:zl(0),:)= eps_halfin(:,:,1:zl(0),:)

#ifdef MPI
           do i=1,nproc-1
              ns = size(potin(:,:,zl(i-1)+1:zl(i)))
                 call mpi_send(potin(:,:,zl(i-1)+1:zl(i)),ns,mpi_double_precision,i,303,mg_comm,ierr)
                 call mpi_send(rhoin(:,:,zl(i-1)+1:zl(i)),ns,mpi_double_precision,i,304,mg_comm,ierr)
                 call mpi_send(eps_halfin(:,:,zl(i-1)+1:zl(i),1:3),3*ns,mpi_double_precision,i,305,mg_comm,ierr)
              enddo
#endif
              idx_s(:) = (/ 1, 1, 1 /)
              idx_e(:) = (/ ngx, ngy, zl(0) /)

              deallocate(potin,rhoin,eps_halfin)

           else
               ns = size(pot(:,:,1:zl(myid)-zl(myid-1)))
#ifdef MPI
               call mpi_recv(pot,ns,mpi_double_precision,0,303,mg_comm,mpi_status_ignore,ierr)
               call mpi_recv(rho,ns,mpi_double_precision,0,304,mg_comm,mpi_status_ignore,ierr)
               call mpi_recv(eps_half(:,:,1:zl(myid)-zl(myid-1),1:3),3*ns,mpi_double_precision,0,305,mg_comm,mpi_status_ignore,ierr)
#endif
               idx_s(:) = (/ 1, 1, zl(myid-1)+1 /)
               idx_e(:) = (/ ngx, ngy, zl(myid) /)

            endif

      end subroutine read_onetep_data


      subroutine test_solution
        use dl_mg_mpi_header
        implicit none

        integer i,j,k,myid,ierr
        real(wp) x,y,z,diff, diffl, maxl, sr(3), q, r, scale, m0, aux
        real(wp) norm_a_loc(3), norm_a(3), &
             norm2l_err, norm2_err, &
             norm2l_ref, norm2_ref, pot_av, lsum(2), gsum(2)

        character(len=100) :: msg=''
        logical unknown

        unknown = .false.

        norm2l_err =0.0_wp; norm2l_ref=0.0_wp

#ifdef MPI
        call mpi_comm_rank(mg_comm,myid,ierr)
#else
        myid = 0
#endif
        select case(model_name)
        case("parabolic_rho","sin_eps")

          diffl =0.0_wp
          maxl = 0.0_wp

          do k=sz,ez
              z=(k-1.0_wp)*dz
              do j=sy,ey
                 y = (j-1.0_wp)*dy
                 do i = sx,ex
                    x =(i-1.0_wp)*dx
                    diffl =  &
                         (-0.5_wp*x*(xl-x)*y*(yl-y)*z*(zl-z) - pot(i,j,k))
                    if (abs(diffl) > maxl) maxl = abs(diffl)
                 enddo
              enddo
           enddo

        case("linear_phi")

           diffl =0.0_wp
           maxl =0.0_wp

           do k=sz,ez
              z=(k-1)*dz
              do j=sy,ey
                 y = (j-1)*dy
                 do i = sx,ex
                    x =(i-1)*dx
                    diffl =  &
                         ((alpha(1)*x+beta(1))*(alpha(2)*y+beta(2))*(alpha(3)*z+beta(3)) - pot(i,j,k))
                    if (abs(diffl) > maxl) maxl = abs(diffl)
                    norm2l_err = norm2l_err + diffl**2
                 enddo
              enddo
           enddo
        case("z_to_power", "z_to_power_pbc")

           !           m0 = gamma(potpow+1.0_wp)**2/gamma(2*potpow+2.0_wp) ! zero mode
           if ( model_name ==  "z_to_power_pbc") then
              diffl = 0.0_wp
              do k=sz,ez
                 z=(k-1)*dz
                 do j=sy,ey
                    y = (j-1)*dy
                    do i = sx,ex
                       x =(i-1)*dx
                       diffl =  diffl + z**potpow * (1.0_wp - z)**potpow
                    enddo
                 enddo
              enddo
              diffl = dx*dy*dz*diffl
#ifdef MPI
              call mpi_allreduce(diffl,m0,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
                   mg_comm,ierr)
#else
              m0 =  diffl
#endif
           else
              m0 = 0.0_wp
           end if
           diffl = 0.0_wp
           maxl  = 0.0_wp

           do k=sz,ez
              z=(k-1)*dz
              do j=sy,ey
                 y = (j-1)*dy
                 do i = sx,ex
                    x =(i-1)*dx
                    diffl =  (z**potpow * (1.0_wp - z)**potpow - m0 - pot(i,j,k))
                    maxl = max(maxl,abs(diffl))
                 enddo
              enddo
           enddo
        case("y_to_power")

           diffl =0.0_wp
           maxl = 0.0_wp
           m0 = gamma(potpow+1.0_wp)**2/gamma(2*potpow+2.0_wp) ! zero mode

           do k=sz,ez
              z=(k-1)*dz
              do j=sy,ey
                 y = (j-1)*dy
                 do i = sx,ex
                    x =(i-1)*dx
                    diffl =  (y**potpow * (1.0_wp - y)**potpow - m0 - pot(i,j,k))
                    maxl = max(maxl,abs(diffl))
                 enddo
              enddo
           enddo

        case("x_to_power")
           diffl =0.0_wp
           maxl = 0.0_wp
           m0 = gamma(potpow+1.0_wp)**2/gamma(2*potpow+2.0_wp) ! zero mode

           do k=sz,ez
              z=(k-1)*dz
              do j=sy,ey
                 y = (j-1)*dy
                 do i = sx,ex
                    x =(i-1)*dx
                    diffl =  (x**potpow * (1.0_wp - x)**potpow - m0 - pot(i,j,k))
                    maxl = max(maxl,abs(diffl))
                 enddo
              enddo
           enddo

        case("erf_eps", "pbe_erf")
           diffl = 0.0_wp
           maxl = 0.0_wp
           pot_av = 0.0_wp
           if (all(periodic))then
              if (pbe%neutr_method == dl_mg_neutralise_with_jellium_vacc) then
                 msg = "WARNING: smart jelium is not consistent with this manufactured solution"
              end if
              !select case(model_name)
              !case(erf_eps)
              lsum = 0.0_wp
              if (allocated(vsteric) .and. pbe_params%use_steric)then
                 do k=sz,ez
                    z=(k-1)*dz
                    do j=sy,ey
                       y = (j-1)*dy
                       do i = sx,ex
                          x =(i-1)*dx
                          lsum(1) = lsum(1) + vsteric(i,j,k)*gauss(ge%sig, ge%mu, ge%scale, x,y,z)
                          lsum(2) = lsum(2) + vsteric(i,j,k)
                       enddo
                    enddo
                 enddo
              else
                do k=sz,ez
                    z=(k-1)*dz
                    do j=sy,ey
                       y = (j-1)*dy
                       do i = sx,ex
                          x =(i-1)*dx
                          lsum(1) = lsum(1) + gauss(ge%sig, ge%mu, ge%scale, x,y,z)
                          lsum(2) = lsum(2) + 1.0_wp
                       enddo
                    enddo
                 enddo
              endif

#ifdef MPI
              call mpi_allreduce(lsum, gsum, 2, MPI_DOUBLE_PRECISION, MPI_SUM,mg_comm,ierr)
#else
              gsum = lsum
#endif
              pot_av = lsum(1) / lsum(2) !(ngxyz(1) * ngxyz(2) * ngxyz(3))
              !write(*,*) 'pot av', gsum
           endif

           do k=sz,ez
              z=(k-1)*dz
              do j=sy,ey
                 y = (j-1)*dy
                 do i = sx,ex
                    x =(i-1)*dx
                    aux = gauss(ge%sig, ge%mu,ge%scale,x,y,z) - pot_av
                    diffl = aux - pot(i,j,k)
                    if (abs(diffl) > maxl) maxl = abs(diffl)
                    norm2l_err = norm2l_err + diffl**2
                    norm2l_ref = norm2l_ref + aux**2
                 enddo
              enddo
           enddo

        case("pbe_point_charge")
           sr = pointsource
           q = sqrt( - pbe_params%lam * (elcharge**2/rbohr)/(kB*pbe_params%temp)/epsh2o&
                * sum( pbe_params%c(:) * pbe_params%q(:)**2))
           scale = 1/epsh2o ! elcharge / rbohr
           diffl = 0.0_wp
           maxl = 0.0_wp
           do k=sz,ez
              z=(k-1)*dz
              do j=sy,ey
                 y = (j-1)*dy
                 do i = sx,ex
                    x =(i-1)*dx
                    r = sqrt((x-sr(1))**2+(y-sr(2))**2+(z-sr(3))**2)
                    diffl = (scale * exp(-q * r)/(r*epsh2o) - pot(i,j,k))
                    if (abs(diffl) > maxl) maxl = abs(diffl)
                 enddo
              enddo
           enddo

        case("poisson_point_charge")
           sr = pointsource
           diffl = 0.0_wp
           maxl = 0.0_wp
           do k=sz,ez
              z=(k-1)*dz
              do j=sy,ey
                 y = (j-1)*dy
                 do i = sx,ex
                    x =(i-1)*dx
                    r = sqrt((x-sr(1))**2+(y-sr(2))**2+(z-sr(3))**2)
                    diffl =  (1.0_wp/r - pot(i,j,k))
                    if (abs(diffl) > maxl) maxl = abs(diffl)
                 enddo
              enddo
           enddo
        case("wirez")
           sr = pointsource
           diffl = 0.0_wp
           maxl = 0.0_wp
           do k=sz,ez
              do j=sy,ey
                 y = (j-1)*dy
                 do i = sx,ex
                    x =(i-1)*dx
                    r = sqrt((x-sr(1))**2+(y-sr(2))**2)
                    diffl = (log(r) - pot(i,j,k))
                    if (abs(diffl) > maxl) maxl = abs(diffl)
                 enddo
              enddo
           enddo
        case("wirex")
           sr = pointsource
           diffl = 0.0_wp
           maxl = 0.0_wp
           do k=sz,ez
              z =(k-1)*dz
              do j=sy,ey
                 y = (j-1)*dy
                 do i = sx,ex
                    r = sqrt((z-sr(2))**2+(y-sr(2))**2)
                    diffl =  (log(r) - pot(i,j,k))
                    if (abs(diffl) > maxl) maxl = abs(diffl)
                 enddo
              enddo
           enddo
        case("pbc", "pbc_eps", "pbe_lin_pbc")
           diffl = 0.0_wp
           maxl = 0.0_wp
           norm2l_err = 0.0_wp
           norm2l_ref = 0.0_wp
           do k=sz,ez
              z=(k-1)*dz
              do j=sy,ey
                 y = (j-1)*dy
                 do i = sx,ex
                    x =(i-1)*dx
                    aux = sin(a11*x)*sin(a21*y)*&
                         sin(a31*z)
                    diffl = abs(aux - pot(i,j,k))
                    !write(*,*) aux, pot(i,j,k), aux/pot(i,j,k)
                    if (diffl > maxl) maxl = diffl
                    norm2l_err = norm2l_err + diffl**2
                    norm2l_ref = norm2l_ref + aux**2
                 enddo
              enddo
           enddo

        case("pbcxy")
           diffl = 0.0_wp
           maxl = 0.0_wp
           do k=sz,ez
              z=(k-1)*dz
              do j=sy,ey
                 y = (j-1)*dy
                 do i = sx,ex
                    x =(i-1)*dx
                    diffl = &
                         ((1.0_wp - z) - pot(i,j,k))
                    if (abs(diffl) > maxl) maxl = abs(diffl)
                 enddo
              enddo
           enddo

        case("pbez")
           diffl = 0.0_wp
           maxl = 0.0_wp
           norm2l_err = 0.0_wp
           norm2l_ref = 0.0_wp
           do k=sz,ez
              z=(k-1)*dz
              do j=sy,ey
                 y = (j-1)*dy
                 do i = sx,ex
                    x =(i-1)*dx
                    aux = pbe1d(z, surface_pot, pbe_params%temp, inv_debye)
                    diffl = aux - pot(i,j,k)
                    !write(300,*) i,j,k, diffl, pbe1d(z, surface_pot, pbe_params%temp, inv_debye), pot(i,j,k)
                    if (abs(diffl) > maxl) maxl = abs(diffl)
                    norm2l_err = norm2l_err + diffl**2
                    norm2l_ref = norm2l_ref + aux**2
                 enddo
              enddo
           enddo

        case default

           unknown = .true.
           if (myid == 0) then
              write(*,*) "test_solution: unknown analitic"
           endif

        end select

        if (.not. unknown) then
#ifdef MPI
           call mpi_reduce(maxl,diff,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
                0,mg_comm,ierr)
           norm_a_loc(1) = norm2l_err
           norm_a_loc(2) = norm2l_ref
           call mpi_reduce(norm_a_loc,norm_a,2,MPI_DOUBLE_PRECISION,MPI_SUM,&
                0,mg_comm,ierr)
           norm2_err = norm_a(1)
           norm2_ref = norm_a(2)
#else
           diff      =  maxl
           norm2_err = norm2l_err
           norm2_ref = norm2l_ref
#endif
           if (myid == 0) then
              write(*,'(/a/a)')    " numerical solution test vs analitical", &
                   "--------------------------------------"
              write(*,*)    "               max |error| =", diff
              write(*,*)    "             || error ||_2 =", sqrt(norm2_err*dx*dy*dz)
              if ( norm2_ref /= 0.0_wp) then
                 write(*,*) "||error||_2/||solution||_2 =",  sqrt(norm2_err/norm2_ref)
              end if
              write(*,*) trim(msg) ! possible warnings from the solutions comparison blocks
           endif
        end if
      end subroutine test_solution


      function pbe1d(x, pot, temp, inv_debye)
        implicit none
        real(wp), intent(in) :: x, pot, temp, inv_debye
        real(wp) pbe1d

        real(wp) y0, z1, z2

        y0 = 0.5_wp * elcharge * pot /( kB *temp)

        z1 = exp (y0) + 1.0_wp &
             + (exp (y0) -1.0_wp) * exp( -inv_debye * x)
        z2 = exp (y0) + 1.0_wp &
             - (exp (y0) -1.0_wp) * exp( -inv_debye * x)

        pbe1d = 2.0_wp * (kB*temp/ (elcharge**2/rbohr)) * log(z1/z2)

      end function pbe1d


      subroutine init_pbe_molecule
        implicit none

        integer i,j,k
        real(wp) x, y, z, x0, y0, z0, x1, y1, z1, ra, r0, r1
        real(wp) :: vs0, vs1, ne, dq

        ra = 1.3_wp; ! atom radius
        ne = 6.0_wp
        dq = 1.0e-5_wp

        ! fix the center of the two atoms
        x0=0.5_wp * xl - 0.9_wp * ra ; y0=0.5_wp * yl; z0=0.5_wp *zl
        x1=0.5_wp * xl + 0.9_wp * ra ; y1=0.5_wp * yl; z1=0.5_wp *zl


        do k=sz,ez
           do j=sy,ey
              do i=sx,ex
                 x=(i-1.0_wp)*dx
                 y=(j-1.0_wp)*dy
                 z=(k-1.0_wp)*dz
                 r0 = sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
                 r1 = sqrt((x-x1)**2 + (y-y1)**2 + (z-z1)**2)
                 if ( r0 < 2 * ra) then
                    vs0 = 0.0_wp
                 else
                    vs0 = 1.0_wp
                 endif
                 if ( r1 < 2 * ra) then
                    vs1 = 0.0_wp
                 else
                    vs1 = 1.0_wp
                 endif
                 vsteric(i,j,k) = vs0 * vs1
              enddo
           enddo
        enddo

        do k=sz,ez
           do j=sy,ey
              do i=sx,ex
                 x=(i-1.0_wp)*dx
                 y=(j-1.0_wp)*dy
                 z=(k-1.0_wp)*dz
                 r0 = sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
                 r1 = sqrt((x-x1)**2 + (y-y1)**2 + (z-z1)**2)
                 ! simpler model
                 if ( r0 < ra .and. r1 < ra) then
                     rho(i,j,k) = 0.0_wp
                  else if ( r0 < ra ) then
                     rho(i,j,k) = dq
                  else if ( r1 < ra ) then
                     rho (i,j,k) = -dq
                  else
                     rho(i,j,k) = 0.0_wp
                  endif


!!$                 if ( r0 < ra .and. r1 < ra) then
!!$                    ! neutral
!!$                    rho(i,j,k) =  ne / (4.0_wp/3.0_wp*pi*ra**3) &
!!$                         * (-2.0_wp +4.0_wp*(2.0_wp - 1.0_wp*(r0+r1)/ra))
!!$                 else if ( r0 < ra ) then
!!$                    rho(i,j,k) = ne / (4.0_wp/3.0_wp*pi*ra**3) &
!!$                         * (-1.0_wp -dq/ne +4.0_wp*(1.0_wp - r0/ra))
!!$                 else if ( r1 < ra) then
!!$                    rho(i,j,k) = ne / (4.0_wp/3.0_wp*pi*ra**3) &
!!$                         * (-1.0_wp + dq/ne +4.0_wp*(1.0_wp  - r1/ra))
!!$                 else
!!$                     rho(i,j,k) = 0.0_wp
!!$                 endif
               enddo
           enddo
        enddo

      end subroutine init_pbe_molecule


      subroutine dump_pbe1d(slice, level)
        implicit none

        character(len=2), intent(in) :: slice
        integer, intent(in) :: level

        integer i,j,k
        real(wp) z

        write(*,*) ' dump pb1d ',  surface_pot, pbe_params%temp, inv_debye
        do k = sz, ez
           z=(k-1.0_wp)*dz
           !write(*,*) 'z', z
           do j = sy, ey
              do i = sx, ex
                 pot(i,j,k) = pbe1d(z, surface_pot, &
                      pbe_params%temp, inv_debye)
              enddo
           enddo
        enddo

         call write_slice(slice,level, sx, sy, sz, pot, "pbe1d")

      end subroutine dump_pbe1d


      function eps_sin(params,x,y,z) result(t)
        implicit none

        type(eps_sin_params_t), intent(in) :: params
        real(wp),intent(in):: x,y,z
        real(wp) t

        integer px, py, pz
        real(wp) a

        px = params%px
        py = params%py
        pz = params%pz
        a  = params%a

        t = 1.0_wp + a * sin(twopi*px*x) * sin(twopi*py*y) * sin(twopi*pz*z)

      end function eps_sin


      function dx_eps_sin(params, x,y,z) result(t)
        implicit none

        type(eps_sin_params_t), intent(in) :: params
        real(wp),intent(in)::x,y,z
        real(wp) t

        integer px,py,pz
        real(wp) a

        px = params%px
        py = params%py
        pz = params%pz
        a  = params%a

        t = a * px * twopi * cos(twopi*px*x) * sin(twopi*py*y) * sin(twopi*pz*z)

      end function dx_eps_sin


      function dy_eps_sin(params, x, y, z) result(t)
        implicit none

        type(eps_sin_params_t),intent(in) :: params
        real(wp),intent(in):: x,y,z
        real(wp) t

        integer px,py,pz
        real(wp) a

        px = params%px
        py = params%py
        pz = params%pz
        a  = params%a

        t = a * twopi * py * sin(twopi*px*x) * cos(twopi*py*y) * sin(twopi*pz*z)

      end function dy_eps_sin


      function dz_eps_sin(params,x,y,z) result(t)
        implicit none

        type(eps_sin_params_t),intent(in) :: params
        real(wp),intent(in):: x,y,z
        real(wp) t

        integer px,py,pz
        real(wp) a

        px = params%px
        py = params%py
        pz = params%pz
        a  = params%a


        t = a * twopi * pz * sin(twopi*px*x) * sin(twopi*py*y) * cos(twopi*pz*z)

      end function dz_eps_sin


      function gauss(sig, mu, s, x, y, z)
        implicit none
        real(wp), intent(in) :: sig, mu,s, x,y,z
        real(wp) gauss

        ! small scale because of the exp term
        gauss = s * exp(-((x-mu)**2+(y-mu)**2+(z-mu)**2)/(2.0_wp*sig**2))/(sig*s2pi)**3

      end function gauss


      function erf_eps(eps0, d, delta, mu, x, y, z)
        implicit none
        real(wp), intent(in) :: eps0,d,delta,mu,x,y,z
        real(wp) erf_eps, r

        r = sqrt((x-mu)**2 + (y-mu)**2 + (z-mu)**2)
        erf_eps = 1 + (eps0-1) * 0.5_wp * (1 + erf((r-d)/delta))

      end function erf_eps


      subroutine rho_pbe_erf(rho)
        use dl_mg_nonlinear_model, only : capexp
        use dl_mg_params, only : hartree
        implicit none

        real(wp) rho(sx:, sy:, sz:)

        integer i,j,k,p
        real(wp) rbol, x,y,z,r,r2, beta, stw_sum, rho_sum
        type(pbe_params_t), pointer :: pb
        real(wp) :: mu(pbe_params%n), lsum(pbe_params%n),gsum(pbe_params%n), &
             c(pbe_params%n)
        real(wp), allocatable :: pot(:,:,:)

        logical fullpbc

        pb => pbe_params

        beta = hartree / pb%temp

        fullpbc = all(periodic)

        allocate(pot(sx:ex,sy:ey,sz:ez))

        !compute vsteric and pot
        do k=sz,ez
           do j=sy,ey
              do i=sx,ex
                 x=(i-1.0_wp)*dx
                 y=(j-1.0_wp)*dy
                 z=(k-1.0_wp)*dz
                 if (pb%use_steric) then
                    vsteric(i,j,k) = stw_erf([ge%mu, ge%mu, ge%mu], &
                         [x,y,z])
!!$                    if ( r > 2.0_wp * ge%sig) then
!!$                       vsteric(i,j,k) = 1.0_wp
!!$                    else
!!$                       vsteric(i,j,k) = 0.0_wp
!!$                    end if
                 else
                    vsteric(i,j,k) = 1.0_wp
                 endif
                 pot(i,j,k) = gauss(ge%sig,ge%mu, ge%scale, x, y, z)
              end do
           enddo
        end do

        ! compute fixed rho
         do k=sz,ez
            do j=sy,ey
               do i=sx,ex
                  x=(i-1.0_wp)*dx
                  y=(j-1.0_wp)*dy
                  z=(k-1.0_wp)*dz
                  r2 = (x-ge%mu)**2 + (y-ge%mu)**2+(z-ge%mu)**2
                  r  = sqrt(r2)
                  rho(i,j,k) = 1.0_wp/alf * (gauss(ge%sig,ge%mu, ge%scale, x,y,z) / ge%sig**2 &
                       * (erf_eps(ge%eps0, ge%d, ge%delta, ge%mu, x, y, z) &
                       * (r2/ge%sig**2-3.0) &
                       - (ge%eps0 -1)*r /(ge%delta * sqrt(pi)) &
                       * exp(-((r-ge%d)/ge%delta)**2)))

               end do
            enddo
         enddo


        ! compute mu
        mu = 0.0_wp
        if (fullpbc) then
           ! remove pot average
           lsum = 0.0
            do k=sz,ez
               do j=sy,ey
                  do i=sx,ex
                     lsum(1) = lsum(1) + vsteric(i,j,k) * pot(i,j,k)
                     lsum(2) = lsum(2) + vsteric(i,j,k)
                  end do
               enddo
            enddo

#ifdef MPI
            call mpi_allreduce(lsum, gsum, 2, MPI_DOUBLE_PRECISION, MPI_SUM,mg_comm,ierr)
#else
            gsum = lsum
#endif

            pot = pot - gsum(1)/gsum(2)
            stw_sum = gsum(2)

            !write(0,*) 'before mu', beta, pb%q, maxval(abs(pot)), maxval(abs(vsteric))
            ! mu
            lsum = 0
            do p = 1, pb%n
               do k=sz,ez
                  do j=sy,ey
                     do i=sx,ex
                        lsum(p) = lsum(p) + vsteric(i,j,k) * exp(-beta * pb%q(p) * pot(i,j,k))
                     end do
                  enddo
               enddo
            enddo

#ifdef MPI
            call mpi_allreduce(lsum, gsum, pb%n, MPI_DOUBLE_PRECISION, MPI_SUM,mg_comm,ierr)
#else
            gsum = lsum
#endif
            mu(:) = - log(gsum(:)/stw_sum)

            if (myid==0)then
               write(*,*) 'input mu', mu, stw_sum
            endif

            ! compute total rho
            lsum = 0.0_wp
            do k=sz,ez
               do j=sy,ey
                  do i=sx,ex
                     lsum(1) = lsum(1) + rho(i,j,k)
              end do
           enddo
        end do
#ifdef MPI
            call mpi_allreduce(lsum(1), rho_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM,mg_comm,ierr)
#else
            gsum = lsum
#endif
            if (myid == 0)then
               write(*,*) "rho sum", rho_sum*dx*dy*dz/(xl*yl*zl)
            endif
            c = pb%c
            ! shift the concentrations to generate a net charge 1 in \rho
            if (pbe%neutr_method /= dl_mg_neutralise_none) then
               if ( pb%q(1) > 0) then
                  c(1) = pb%c(1) + 1.0_wp/(pb%q(1)*stw_sum*dx*dy*dz)
               else
                  c(2) = pb%c(2) + 1.0_wp/(pb%q(2)*stw_sum*dx*dy*dz)
               end if

            else
               c = pb%c
               ! bump the total charge up from 0
               rho = rho - 1.0_wp/(xl*yl*zl)
            end if
            if (myid == 0)then
               write(*,*) "rho con", stw_sum*dx*dy*dz/(xl*yl*zl), &
                    sum(c*pb%q) * stw_sum*dx*dy*dz
            endif
         else
            c = pb%c
         end if

         ! subtract bolztmann rho
         if (.not. pb%linearised) then
           do k=sz,ez
              do j=sy,ey
                  do i=sx,ex
                     rbol = 0.0_wp
                     do p =1, pb%n
                        rbol = rbol + pb%lam/alf * c(p) * pb%q(p) &
                             * vsteric(i,j,k) &
                             *  exp(- beta * pb%q(p) * pot(i,j,k))
                     end do
                     rho(i,j,k) = rho(i,j,k) - rbol
                  enddo
               enddo
            enddo
         else
            do k=sz,ez
              do j=sy,ey
                  do i=sx,ex
                     rbol = 0.0_wp
                     do p =1, pb%n
                        rbol = rbol + pb%lam/alf * c(p) * pb%q(p) &
                             * vsteric(i,j,k) &
                             *  (- beta * pb%q(p) * pot(i,j,k))
                     end do
                     rho(i,j,k) = rho(i,j,k) - rbol
                  enddo
               enddo
            enddo
         endif

      end subroutine rho_pbe_erf


      function stw_erf(r0,r)
        implicit none
        real(wp), intent(in) :: r0(3), r(3)
        real(wp) stw_erf

        real(wp) y

        y = sqrt(sum((r-r0)**2))

        stw_erf = 0.5_wp * (1.0_wp + erf((y-stw_params%d)/stw_params%delta))

      end function stw_erf


      ! computes the number of nearest neighbors
      function n_neighbor(i,j,k) result (n)
        implicit none
        integer, intent(in) :: i,j,k
        integer n

        integer imin,imax,jmin,jmax,kmin,kmax

        n=0

        imin = 1; jmin = 1; kmin = 1
        imax = ngx; jmax = ngy; kmax = ngz

        if (i > imin .and. i < imax .and. &
            j > jmin .and. j< jmax .and. &
            k > kmin .and. k < kmax) then
           ! bulk
           n = 6
        else if ( (i == imin .or. i == imax) .and. &
             j > jmin .and. j< jmax .and. &
             k > kmin .and. k < kmax) then
           !side
           n = 5
        else if ( i > imin .and. i < imax .and. &
             ( j == jmin .or. j == jmax) .and. &
             k > kmin .and. k < kmax) then
           !side
           n = 5
        else if (i > imin .and. i < imax .and. &
            j > jmin .and. j < jmax .and. &
            (k == kmin .or. k == kmax)) then
           !side
           n = 5
        else if ((i == imin .or. i == imax) .and. &
            (j == jmin .or. j == jmax) .and. &
            k > kmin .and. k < kmax) then
           !edge
           n = 4
        else if ((k == kmin .or. k == kmax) .and. &
             (j == jmin .or. j == jmax) .and. &
             i > imin .and. i < imax) then
           !edge
           n = 4
        else if ((k == kmin .or. k == kmax) .and. &
             (i == imin .or. i == imax) .and. &
             j > jmin .and. j < jmax) then
           !edge
           n = 4
        else if ((i == imin .or. i == imax)  .and. &
             (j == jmin .or. j == jmax) .and. &
             (k == kmin .or. k == kmax)) then
           !corner
           n = 3
        endif

      end function n_neighbor


      function pot_perturbation(a,px,py,pz,x,y,z) result(w)
      implicit none

      integer,intent(in)   :: px,py,pz
      real(wp), intent(in) :: a,x,y,z

      real(wp) w

      w = a * ( sin(px*pi*x) + sin(py*pi*y) + sin(pz*pi*z))

    end function pot_perturbation


  end subroutine poisson_test


  !> sets boundary condition arguments for dl_mg
  !! solver using MPI topology periodicity or user provided value
  !! computes also the lattice step accordingly.
  subroutine set_bc(periodic, xyz, ng_xyz, bc, dx, dy, dz)
    use dl_mg
    implicit none

    logical, intent(in) :: periodic(3)
    real(kind=wp), intent(in) :: xyz(3)
    integer, intent(in) :: ng_xyz(3)
    integer, intent(out) :: bc(3)
    real(wp), intent(out) :: dx, dy, dz

    integer dn(3)
    real(wp)  dxyz(3)

    write(*,*) 'test periodicity', periodic

    dn = 1
    bc = DL_MG_BC_DIRICHLET
    where (periodic)
       dn = 0
       bc = DL_MG_BC_PERIODIC
    end where

    dxyz = xyz /real(ng_xyz - dn, wp)

    dx = dxyz(1)
    dy = dxyz(2)
    dz = dxyz(3)

  end subroutine set_bc


  subroutine abort(msg)
    implicit none
#ifdef MPI
      include 'mpif.h'
#endif

    character(len=*), intent(in) :: msg

    integer ierr

    write(0,*) msg

#ifdef MPI
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
    STOP
#endif

  end subroutine abort

end module mg_tests
