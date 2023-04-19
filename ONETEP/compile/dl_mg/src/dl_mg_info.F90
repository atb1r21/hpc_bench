!> \brief Subroutines for logs; also stores the solver version

module dl_mg_info
  use dl_mg_params, only : wp
  use dl_mg_mpi, only :  get_myid ! used in several functions below
  use dl_mg_common_data, only : report_unit, report_file
  implicit none
  private

  character(len=*), parameter :: version_string="3.1.0-rc.17 (18/08/2021)"

  integer, parameter, public :: info_newton_tag=222, info_mg_tag=555, info_cg_tag=777

! info solver status
  type solver_status_t
     integer status, iter, maxiter
     real(wp)  sol_norm, res_norm, res_ratio, tol
     integer res_ratio_dist(21)
  end type solver_status_t

  type defco_status_t
     integer iter, maxiter
     integer ierror
     real(wp) sol_norm, error_norm, res_norm
     real(wp) target_res_norm
     real(wp) mu, mu_prev ! stubs, should be arrays
  end type defco_status_t

  type(solver_status_t), target :: mg_status, newton_status, cg_status
  type(defco_status_t) :: defco_status

  logical :: debug_write=.false., initialised=.false.
  logical :: collect_res_hist = .false. !> enables residual ratio histogram
  logical :: print_iter_info = .true. !> if true prints residual and solution norms at each iteration

  public mg_info_init, mg_info_set,  mg_info_write, write_comms_structure, &
       mg_info_iteration_details, accumulate_residual_ratio, &
       debug_write,  write_loop_block_info, get_version,&
        write_nonlinear_params, &
       charge_neutrality_report, info_get_res_norm, &
       charge_neutralisation_info, &
       newton_status, defco_status

contains

  subroutine mg_info_init(report_file_)
!$  use omp_lib
    use dl_mg_mpi_header
    use dl_mg_common_data
    implicit none
    character(len=*) :: report_file_

    integer myid, ierr, lval
    character(len=128) val
    logical :: report_file_exists


    myid = get_myid()

    ! There my be some trouble if the report unit is set to a valid value only
    ! on one rank in the caller application, and that rank does not map to 0 in mg_comm

    ! assume that only one callee rank provides a positive report unit or
    ! all provide the same positive value
#ifdef MPI
    lval = report_unit
    call mpi_allreduce(lval, report_unit, 1, mpi_integer, mpi_max, mg_comm, ierr)
#endif

    call get_environment_variable("DL_MG_DEBUG_PRINT",value=val,length=lval)
    if (lval > 0)then
       debug_write = .true.
    endif

    ! switch collection and output of residual ratio histogram
    call get_environment_variable("DL_MG_RES_HIST",value=val,length=lval)
    if (lval > 0)then
       collect_res_hist = .true.
    endif

    call get_environment_variable("DL_MG_NO_ITER_INFO",value=val,length=lval)
    if (lval > 0)then
       print_iter_info = .false.
    endif


    report_file = report_file_
    ! reset report file via environment variable, if needed
    call get_environment_variable("DL_MG_LOG_FILE",value=val,length=lval)
    if (lval > 0) then
       report_file = val
    end if

    ! one should check that the report file is non-empty string
    ! if so, a default name should be set ( time added in ?)

    mg_status%res_ratio_dist(:)=0

    ! this is unnecessary when Newton method is not used but dl_mg_init does not know
    ! which method will be used. Anyway it's a trifle
    newton_status%res_ratio_dist(:)=0


    ! extraordinary writes to report unit by ranks > 0 go in special files
    if(myid > 0) then
       write(val, *) myid
       report_file=report_file(1:len(report_file))//".rank"//adjustl(val)
    endif

    initialised = .true.
    if (myid > 0) return
    ! JCW: Check if report file exists, if so, append, if not create new
    inquire(file=report_file,exist=report_file_exists)
    if (report_file_exists) then
      ! Exists, so append
      open(report_unit,file=report_file,status="old",position="append")
    else
      ! Does not exist, create new
      open(report_unit,file=report_file,status="new")
    end if

    write(report_unit,'(/a)') " DL_MG "//version_string
#ifndef MPI
    write(report_unit,*) "<< Compiled without MPI >>"
#else
    write(report_unit,*) "<< Compiled with MPI >>"
#endif
    write(report_unit,*) "Fine grid   : ", nx, ny, nz
    write(report_unit,*) "Grid step   : ", dx, dy, dz
    write(report_unit,*) "Coarse grid : ", nxc, nyc, nzc
    write(report_unit,*) "MG levels   : ", mg_levels
!#ifdef MPI
!    write(report_unit,*) "Full aggregation level :", full_agg_level
!    write(report_unit,*) "Mininum grid size for full aggregation  :", full_agg_size
!    if ( full_agg_size == 0) then
!       write(report_unit,*) "  If grid size for aggregation is 0 "
!       write(report_unit,*) "  full aggregation is triggered at "
!       write(report_unit,*) "  level", full_agg_level
!    endif
!#endif

    write(report_unit,*)
    !this is not working because dl_mg_init_nonlin is called after dl_mg_init
    !if(use_nonlinear) then
    !   write(report_unit,*) 'non-linear problem'
    !endif
    !write(report_unit,*)
!$  write(report_unit, *) 'using ', omp_get_max_threads(), 'OpenMP threads'
!$  write(report_unit,*)
!    flush(report_unit)
    close(report_unit)
    open(report_unit,file=report_file,status="old",position="append")

  end subroutine mg_info_init


  subroutine mg_info_set(status, iter, resnorm, res_ratio, solnorm, tol, maxiter, eq_tag)
    use dl_mg_mpi_header
    implicit none
    integer, intent(in) :: status,iter,maxiter
    real(wp), intent(in) :: resnorm, res_ratio, solnorm,tol
    integer, intent(in) :: eq_tag

    integer myid, ierr
    type(solver_status_t), pointer :: p => null()


    myid = get_myid()

    if (myid > 0) return

    select case (eq_tag)
    case (info_newton_tag)
       p => newton_status
    case ( info_mg_tag)
       p => mg_status
    case(info_cg_tag)
       p => cg_status
    end select
    p%status   = status
    p%res_norm = resnorm
    p%res_ratio= res_ratio
    p%sol_norm = solnorm
    p%iter     = iter
    p%tol      = tol
    p%maxiter  = maxiter

  end subroutine mg_info_set


  subroutine mg_info_write(eq_type)
    use dl_mg_mpi_header
    use dl_mg_common_data, only : mg_comm, DL_MG_SUCCESS
    use dl_mg_errors, only : DL_MG_ERR_NITER, DL_MG_ERR_RES_RATIO
    use dl_mg_params, only :  EQ_PBE_NEWTON,  EQ_PBE_NEWTON_NOSTERIC
    use dl_mg_nonlinear_model, only : temp
    use dl_mg_convergence_params, only : use_cg => use_cg_
    implicit none

    integer, intent(in) :: eq_type
    integer i, ip, jp, myid, nprint1, nprint2, ierr
    logical lop
    type(solver_status_t), pointer :: p => null()
    character(len=63) header

!!$     if (status /= mg_dcse_status%status) then
!!$        print*
!!$        print*, 'mg_dcse_info: wrong value of parameter status'
!!$        print*
!!$        return
!!$     endif

    myid = get_myid()

    if (myid > 0) return

    inquire (unit=report_unit,opened=lop)
    if (.not. lop) then
       open (unit=report_unit,file=report_file,status="old", position="append")
    endif

    nprint2 = 1
    if (eq_type == EQ_PBE_NEWTON &
         .or. eq_type == EQ_PBE_NEWTON_NOSTERIC) nprint2 = 2

    nprint1 = 1
    if (use_cg) nprint1 = 2
    
    do ip = 1, nprint2
       select case (ip)
       case (2)
          p => newton_status
          header="Newton iteration report"
       case(1)
          do jp = 1, nprint1
             select case(jp)
             case(2)
                p => cg_status
                header="CG iteration report"
             case(1)
                p => mg_status
                header="MG iteration report"
             end select
          end do
       end select

       write(report_unit,'(1x/a)') header
       select case(p%status)
       case (DL_MG_SUCCESS)
          write(report_unit,*) 'dl_mg_info: convergence reached after ', &
               p%iter, 'steps.'
       case(DL_MG_ERR_NITER)
          write(report_unit,*) 'dl_mg_info: failed to converge after ', &
               p%iter, 'steps.'
       case(DL_MG_ERR_RES_RATIO)
          write(report_unit,*) 'dl_mg_info: failed to converge after ', &
               p%iter, 'steps.'
          write(report_unit,'(a/a/a/a)') &
               '            Residual ratio is larger than the set limit.',&
               '            Please check the residual distribution histogram',&
               '            which can be activated by defining DL_MG_RES_HIST',&
               '            The solution tolerance parameters might need to be decreased for this problem.'
       case default
          write(report_unit,*) 'dl_mg_info: unknown status value'
          write(report_unit,*) ''
          return
       end select
       write(report_unit,*) 'residual norm ', p%res_norm
       write(report_unit,*) 'solution norm ', p%sol_norm
       write(report_unit,*) 'solution tol  ', p%tol
       !write(report_unit,*) 'max iterations', p%maxiter
       if (collect_res_hist) then
          write(report_unit,*)
          write(report_unit, *) 'residual ratio histogram in (0,1)'
          write(report_unit, '(F5.3,1x,I9)') ( (i-1.0)*0.05 + 0.025 , p%res_ratio_dist(i), i=1,20)
          write(report_unit,*)
          write(report_unit,*)'residual ratios  >= 1 ',  p%res_ratio_dist(21)

          write(report_unit,*) ''
       end if
    end do

    close(report_unit)
    open (unit=report_unit,file=report_file,status="old", position="append")
    !    flush(report_unit)

  end subroutine mg_info_write


  subroutine write_comms_structure(mg)
    use dl_mg_types
    use dl_mg_mpi_header
    use dl_mg_common_data, only : mg_comm, npx, npy, npz, mg_levels
    use dl_mg_errors
    implicit none

    type(mg_t), intent(in) :: mg(:)

#ifdef MPI
    integer a(9,mg_levels,npx*npy*npz), a0(9,mg_levels), pdims(3), &
         d, na, nactive, i, j, k, t, myid, iproc, ierr
    logical lop
    character(len=128) caux
    character(len=256):: line


    if ( .not. check_assertion(initialised)) then
       call handle_error(DL_MG_ERR_UNSPECIFIED, &
           msg = "write_comms_structure called before mg_info_init")
    endif

    a(:,:,:) = 0
    if ( 3+6+6*mg_levels > len(line) ) then
       write(0, *)' warning: write_comms_structure : line string too short!'
    endif
    do i=1, len(line)
       line(i:i) =" "
    enddo

    call mpi_comm_rank(mg_comm, myid, ierr)

    do i = mg_levels, 1, -1
       a0(1:3,i) = mg(i)%coords(:)
       if ( mg(i)%active) then
          a0(4:6,i) = (/ mg(i)%sx, mg(i)%sy, mg(i)%sz /)
          a0(7:9,i) = (/ mg(i)%ex, mg(i)%ey, mg(i)%ez /)
       else
          a0(4:6,i) = 0
          a0(7:9,i) = -1
       end if
    enddo

    call mpi_gather(a0, 9*mg_levels, mpi_integer, a, 9*mg_levels, mpi_integer, 0,&
         mg_comm,ierr)

    if (myid == 0) then
       inquire (unit=report_unit, opened=lop)
       if (.not. lop) then
          open (unit=report_unit,file=report_file,status="old",position="append")
          !write(report_unit,*) " report unit reopend, check the  calling code"
       endif
       write(report_unit,*)
       caux(1:3)='xyz'
       pdims = (/ npx, npy, npz/)

       do d=1,3
          write(report_unit,'(a)') "  ranks  right coordinate in "//caux(d:d)//" direction"
          do i = 0, pdims(d)-1

             iproc=-1
             nactive = 0
             ! find in slice i, dimesion d  a rank that
             ! have the largest number of active levels
             do k = 1, npx*npy*npz
                if ( a(d, mg_levels, k) == i) then
                   !check the number  of active levels
                   na = 0
                   do t = 1, mg_levels
                      if ( a(3+d, t, k) /= 0) then
                         na = na + 1
                      endif
                   end do
                   if ( na > nactive ) then
                      nactive = na
                      iproc = k
                   endif
                endif
             enddo
             if ( iproc < 0) then
                write(report_unit,*) "error in  write_comms_structure, iproc < 0, jumping the subroutine"
                close(report_unit)
                return
             endif
             ! indent 0 spaces'
             k=1
             write(line(k:k+5), '(I6,1x)')  a(d,mg_levels,iproc)
             k=k+6
             do t=mg_levels,1,-1
                if ( a(6+d, t, iproc) - a(3+d, t, iproc) + 1 > 0 ) then
                   write(line(k:k+5), '(I6)')   a(6+d, t, iproc)
                else
                   write(line(k:k+5), '(a)')  '      '
                endif
                k=k+6
             end do
             write(report_unit,'(a)') line(1:len_trim(line))
          end do
          write(report_unit,*)
       end do

       close(report_unit)
       open (unit=report_unit,file=report_file,status="old", position="append")
    endif
#endif

  end subroutine write_comms_structure


  subroutine write_loop_block_info(blk, mg_levels, last_active_lvl, &
       min_thread_chunk, block_size)
    use dl_mg_mpi_header
    use dl_mg_types
    use dl_mg_errors
    implicit none

    type(block_list_t), intent(in) ::  blk(:)
    integer, intent(in) :: mg_levels, last_active_lvl, min_thread_chunk
    integer, intent(in) :: block_size


    logical lop
    integer i, myid

    if ( .not. check_assertion(initialised)) then
      call handle_error(DL_MG_ERR_UNSPECIFIED, &
           msg="write_loop_block_info called before mg_info_init")
   endif

      myid = get_myid()

    if (myid == 0) then
       inquire (unit=report_unit,opened=lop)
       if (.not. lop) then
          open (unit=report_unit,file=report_file,status="old",position="append")
       endif

       write(report_unit,*) "Loop blocks info:"
       write(report_unit,*) "Smallest thread chunk in y:", min_thread_chunk
       write(report_unit,*) "Loop block size in y:      ", block_size
       write(report_unit,*) "block dimensions at each level:"
       write(report_unit,*) 'levels', (i, i = mg_levels, last_active_lvl, -1)
       write(report_unit,*) 'x     ', (blk(i)%dims(1),  i = mg_levels, last_active_lvl, -1)
       write(report_unit,*) 'y     ', (blk(i)%dims(2),  i = mg_levels, last_active_lvl, -1)
       write(report_unit,*) 'z     ', (blk(i)%dims(3),  i = mg_levels, last_active_lvl, -1)
       write(report_unit,*) 'n blk ', (size(blk(i)%start, dim=2),  i = mg_levels, last_active_lvl, -1)
       write(report_unit,*)
       close(report_unit)
       open (unit=report_unit,file=report_file,status="old", position="append")
    endif


  end subroutine write_loop_block_info


  subroutine write_nonlinear_params(eq_type, n, c, q, temp, lam)
    use dl_mg_mpi_header
    use dl_mg_params
    use dl_mg_errors
    use dl_mg_common_data, only : nx, ny, nz, fullpbc
    use dl_mg_nonlinear_model, only : steric_weight_sum
    implicit none
    integer, intent(in)  :: eq_type, n
    real(wp), intent(in) :: c(n), q(n), temp, lam


   integer myid, i
   logical lop

   if (.not. check_assertion(initialised)) then
      call handle_error(DL_MG_ERR_UNSPECIFIED, &
          msg = "write_nonlinear_params called before mg_info_init")
   endif

      myid = get_myid()

    if (myid == 0) then
       inquire (unit=report_unit,opened=lop)
       if (.not. lop) then
          open (unit=report_unit,file=report_file,status="old",position="append")
       endif

       select case(eq_type)
       case (EQ_LINEAR_PBE)
          write(report_unit, '(a)') "Solving linearised PBE with steric weight:"
       case( EQ_LINEAR_PBE_NOSTERIC)
          write(report_unit, '(a)') "Solving linearised PBE without steric weight:"
       case (EQ_PBE_FAS)
          write(report_unit, '(a)') "Solving full PBE with steric weight using FAS algorithm:"
       case( EQ_PBE_FAS_NOSTERIC)
          write(report_unit, '(a)') "Solving full PBE without steric weight using FAS algorithm:"
       case (EQ_PBE_NEWTON)
          write(report_unit, '(a)') "Solving full PBE with steric weight using Newton algorithm:"
       case( EQ_PBE_NEWTON_NOSTERIC)
          write(report_unit, '(a)') "Solving full PBE without steric weight using Newton algorithm:"
       end select
       write(report_unit, '(/,a,1x,(I4))') &
            "Ionic species:                       ", n
       write(report_unit, '(a,1x,100(F4.1,1x))') &
            "Ion charges (e units):         ", q
       write(report_unit, '(a,1x,100(E10.3,1x))') &
            "Input concentrations (rB**-3 units): ", c
       write(report_unit, '(a,1x,E10.3)') &
            "Temperature (K):                     ", temp
       write(report_unit, '(a,1x,E10.3)') &
            "Multiplicative factor lambda:        ", lam
       if (fullpbc) then
          write(report_unit, '(a,/,a37,1x,E10.3,/)') &
               "Concentration renormalisation", "Gamma:", &
               steric_weight_sum/(real(nx,wp)*ny*nz)
       end if
       write(report_unit,*)
       close(report_unit)
       open (unit=report_unit,file=report_file,status="old", position="append")
    end if

  end subroutine write_nonlinear_params

  !> this won't write level 1 data if not active
  subroutine mg_info_iteration_details(test_tag, level,iter, &
       resnorm, res_ratio, vecnorm,tol)
    use dl_mg_mpi_header
    use dl_mg_convergence_params, only : tconv_level_1, tconv_cg, &
         tconv_vcycle, tconv_newton
    implicit none
    integer, intent(in) :: test_tag, level, iter
    real(kind(0.d0)), intent(in) :: resnorm, res_ratio, vecnorm, tol

    integer myid, ierr
    logical lop
    character(len=5) tag

    if (debug_write) then
       myid = get_myid()
       if (myid == 0) then
          inquire (unit=report_unit,opened=lop)
          if (.not. lop) then
             open (unit=report_unit,file=report_file,status="old",position="append")
          endif
          if ( level == 1 ) then
             select case(test_tag)
             case (tconv_level_1)
                tag="MG1  "
             case (tconv_cg)
                tag="CG1  "
             end select
          else
             select case(test_tag)
             case (tconv_vcycle)
                tag="MGt  "
             case (tconv_newton)
                tag="NWT  "
             case (tconv_cg)
                tag="CGt  "
             end select
          end if
          write(report_unit,'(a, I3,1x,2(E22.15,1x), I3,1x, E10.2)')&
               tag, iter, resnorm, vecnorm, level, tol
       end if
!       flush(report_unit)
    end if

  end subroutine mg_info_iteration_details


  subroutine charge_neutrality_report(fixed, ions, obc_deviation, myid)
    use dl_mg_params, only : dl_mg_neutralise_with_ions_fixed, &
         dl_mg_neutralise_with_ions_auto, &
         dl_mg_neutralise_with_ions_auto_linear
    use dl_mg_nonlinear_model, only : mu, &
         neutralisation_method
    implicit none

    real(wp), intent(in) :: fixed(:), ions(:), obc_deviation(:)
    integer, intent(in), optional  :: myid

    integer idloc, ierr

    if (present(myid)) then
       idloc = myid
    else
       idloc = get_myid()
    end if

    if (idloc == 0) then
       write(report_unit,'(a/a/)') '',' Charge neutrality test: '
       write(report_unit,*)                       '                sum(c*q)       sum(c)       sum(c)/sum(N)'
       write(report_unit,*)                       '---------------------------------------------------------'
       write(report_unit,'(a,e11.3,4x,e11.3)'   ) 'fixed charge :', fixed(1:2)
       write(report_unit,'(a,e11.3,2(4x,e11.3))') 'ions         :', ions
       write(report_unit,*)                       '---------------------------------------------------------'
       write(report_unit,'(/a, 100e11.3)') 'chemical potential shifts', mu ! this will break down for more than 100 ion species
       select case (neutralisation_method)
       case (dl_mg_neutralise_with_ions_fixed, &
            dl_mg_neutralise_with_ions_auto, &
            dl_mg_neutralise_with_ions_auto_linear)
          write(report_unit,'(a,100e11.3)') &
               'deviation from the open boundary concentration:', &
               obc_deviation(1:size(obc_deviation)-1)
          write(report_unit,'(a,e11.3)') 'sum of ion mixing ratios                      :',  obc_deviation(size(obc_deviation))
       end select
    end if

  end subroutine charge_neutrality_report


  subroutine charge_neutralisation_info
    use dl_mg_params
    use dl_mg_nonlinear_model, only : neutralisation_method
    use dl_mg_neutralisation_with_ions, only : counterion_ratios
    implicit none

    integer myid
    character(len=100) line

    myid = get_myid()

    if (myid == 0) then
       write(report_unit,'(a)') "Neutralisation method info:"
       select case(neutralisation_method)
       case(dl_mg_neutralise_none)
          write(report_unit,'(a/)') '   no charge neutralisation'
       case(dl_mg_neutralise_with_jellium_uniform)
          write(report_unit,'(a/)') '   charge neutralisation with jellium'

       case(dl_mg_neutralise_with_jellium_vacc)
          write(report_unit,'(a/)') '   charge neutralisation with jellium in accessible volume'

       case(dl_mg_neutralise_with_ions_fixed, &
            dl_mg_neutralise_with_ions_auto, &
            dl_mg_neutralise_with_ions_auto_linear)
          write(report_unit,'(a)') '   using shifted concetrations'
          if (neutralisation_method == dl_mg_neutralise_with_ions_fixed) then
             write(report_unit,'(a,100e11.3/)') &
                  '      with fixed mixing fraction:', counterion_ratios
          else
             line = '      mixing fraction computed selfconsistently'
             if (neutralisation_method == dl_mg_neutralise_with_ions_auto_linear) then
                line = trim(line)//' with linear approximation'
             end if
             write(report_unit,'(a/)') line
          end if
       end select
    end if
  end subroutine charge_neutralisation_info


  subroutine accumulate_residual_ratio(eq_tag, iter, resnorm, potnorm, ratio, dmu)

    implicit none

    integer, intent(in)  :: eq_tag, iter
    real(wp), intent(in) :: resnorm, potnorm, ratio
    real(wp), intent(in) :: dmu !< max( abs(mu - mu_prev))


    integer i, myid, ierr
    character(len=5) :: hook, header*60
    type(solver_status_t), pointer :: p => null()
    logical lop

       myid = get_myid()

       if (myid == 0) then

          select case(eq_tag)
          case(info_newton_tag)
             p => newton_status
          case(info_mg_tag)
             p => mg_status
          case(info_cg_tag)
             p => cg_status
          end select


          if (iter > 1 .and. collect_res_hist) then
             do i = 1, 20
                if ( ratio < i*0.05_wp ) then
                   p%res_ratio_dist(i) = p%res_ratio_dist(i) + 1
                   exit
                endif
             enddo

             if ( ratio >= 1.0_wp)  p%res_ratio_dist(21) =  p%res_ratio_dist(21) + 1
          end if

          if (print_iter_info) then
             inquire (unit=report_unit,opened=lop)
             if (.not. lop) then
                open (unit=report_unit,file=report_file,&
                     status="old",position="append")
             end if
             select case(eq_tag)
             case(info_newton_tag)
                hook = "NWT: "
             case(info_mg_tag)
                hook = "MG:  "
             case(info_cg_tag)
                hook = "CG:  "
             end select

             if (iter == 1) then
                header = hook//" iter    residual       solution"
                if (dmu > 0.0_wp) then
                   header = trim(header)//"     max(abs(dmu))"
                end if
                write(report_unit,'(a)') header
             end if
             if ( dmu < 0.0_wp) then
                write(report_unit,100) hook, iter, resnorm, potnorm
             else
                write(report_unit,100) hook, iter, resnorm, potnorm, dmu
             end if
          end if
       end if
100    format(a, i4, 3(2x,e13.7))

   end subroutine accumulate_residual_ratio


   subroutine get_version(vs)
     implicit none
     include 'dl_mg_common_params.inc'
     character(len=DL_MG_VERSION_STRING_LENGTH), intent(out) :: vs

     vs = version_string
   end subroutine get_version


   function info_get_res_norm(key) result(norm)
     use dl_mg_errors
     implicit none

     character(len=*), intent(in) :: key

     real(wp) :: norm

     type(solver_status_t), pointer :: p

     select case(key(1:2))
     case("mg","MG")
        p => newton_status
     case ("nw", "NW")
        p => mg_status
     case default
        call handle_error(DL_MG_ERR_UNSPECIFIED, &
           msg = "Wrong key values in input of info_get_res_norm. Try mg or nw.")
    end select
     norm = p % res_norm

   end function info_get_res_norm


end module dl_mg_info
