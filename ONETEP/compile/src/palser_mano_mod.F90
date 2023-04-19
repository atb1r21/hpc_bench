! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!         Palser-Manolopoulos kernel optimisation module         !
!                                                                !
! This module implements variants of the method by               !
! Palser and Manolopoulos                                        !
!         Phys. Rev. B. 58(19), 12704 (1998)                     !
! for optimising the density kernel using purification-based     !
! algorithms.                                                    !
!----------------------------------------------------------------!
! Written by Chris-Kriton Skylaris on 24/07/2006                 !
!================================================================!

module palser_mano

  implicit none

  private

  public :: palser_mano_kernel_optimise

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine palser_mano_kernel_optimise(k_new, &
       ham, olap, s_inv, n_occ, num_iter)

    !================================================================!
    ! This subroutine implements the "Canonical purification" method !
    ! of Palser and Manolopoulos                                     !
    !         Phys. Rev. B. 58(19), 12704 (1998)                     !
    ! for obtaining the density kernel which optimises the band      !
    ! structure energy for a fixed Hamiltonian. When no truncation   !
    ! is applied to the density kernel, the results obtained are     !
    ! equivalent to the density kernel that would be obtained        !
    ! by diagonalising the Hamiltonian and building a kernel from    !
    ! the eigenvectors obtained.                                     !
    !----------------------------------------------------------------!
    ! Note: The method is suitable to sparse matrices with truncation!
    !       applied. The iterations are stopped and maximum          !
    !       convergence is assumed when either of the following      !
    !       conditions is satisfied:                                 !
    !       1) The stable fixed point value c_n is outside [0,1].    !
    !       2) The band structure energy ceases to be monotonically  !
    !          decreasing.                                           !
    !       3) The absolute change in band structure energy per elec !
    !          becomes less than the delta_e_thresh parameter.       !
    !----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 24/07/2006                 !
    ! Extracted for embedding by Robert Charlton, 28/07/17.          !
    ! It might be worth looking into new ways of building the DK.    !
    !================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout, verbose, CRLF
    use rundat, only: pub_output_detail, pub_cond_calculate, pub_debug, &
         pub_debug_on_root
    use services, only: services_flush
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
        sparse_embed_destroy, sparse_embed_extremal_eigenvalue, &
        sparse_embed_product, sparse_embed_axpy, sparse_embed_scale, &
        sparse_embed_trace, sparse_embed_copy, sparse_embed_num_rows
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout)    :: k_new   ! optimised kernel
    type(SPAM3_EMBED), intent(in)       :: ham     ! hamiltonian
    type(SPAM3_EMBED), intent(in)       :: olap    ! overlap
    type(SPAM3_EMBED), intent(in)       :: s_inv   ! inverse overlap
    integer, intent(in)           :: n_occ   ! number of up/down electrons
    integer, optional, intent(in) :: num_iter! maximum number of PM iterations

    ! Local variables
    type(SPAM3_EMBED) :: k_old    ! previous iteration kernel
    type(SPAM3_EMBED) :: ks       ! current K.S
    type(SPAM3_EMBED) :: ksk      ! current K.S.K
    type(SPAM3_EMBED) :: ksks     ! K.S.K.S structure, for calculating efermi
    type(SPAM3_EMBED) :: sks      ! current S.K.S
    type(SPAM3_EMBED) :: sinv_ham ! inverse overlap times hamiltonian
    real(kind=DP) :: min_en         ! minimum hamiltonian eigenvalue
    real(kind=DP) :: max_en         ! maximum hamiltonian eigenvalue
    real(kind=DP) :: average_en     ! everage hamiltonian eigenvalue
    real(kind=DP) :: lambda         ! PM initialisation parameter
    real(kind=DP) :: i_fac          ! PM purification formula factor
    real(kind=DP) :: h_fac          ! PM purification formula factor
    real(kind=DP) :: c_numer        ! numerator of c_n
    real(kind=DP) :: c_denom        ! denominator of c_n
    real(kind=DP) :: c_n            ! unstable fixed point of PM formula
    real(kind=DP) :: scale_fac      ! PM formula scaling factor
    real(kind=DP) :: band_energy    ! current band structure energy
    real(kind=DP) :: old_band_energy! previous band structure energy
    real(kind=DP), parameter :: delta_e_thresh =1.0E-10_DP ! energy/elec threshold
    real(kind=DP), parameter :: delta_e_thresh_small_enough = 1.0E-8_DP
    real(kind=DP), parameter :: c_numer_small_enough = 1.0E-13_DP
    real(kind=DP), parameter :: eval_prec = 0.0001_DP ! eigenvalue bound precision
    real(kind=DP) :: lbound         ! minimum kernel eigenvalue
    real(kind=DP) :: ubound         ! maximum kernel eigenvalue
    real(kind=DP) :: mid_occ        ! closest kernel eigenvalue to 0.5
    real(kind=DP) :: delta_e_per_elec
    integer       :: iter           ! iteration counter
    integer       :: loc_num_iter   ! local maximum number of PM iterations
    logical       :: quitearly      ! loop exit flag
    integer       :: quitreason     ! termination criterion satisfied
    character(len=1024) :: error_message
    real(kind=DP) :: debug_trace    ! matrix trace for debugging purposes

    if (pub_debug_on_root) then
       write(stdout,'(a)') '  '
       write(stdout,'(a)') 'DEBUG: Entering palser_mano_kernel_optimise'
    endif

    ! pdh: nothing to be done if no electrons of this spin
    if (n_occ < 1) then
       call sparse_embed_scale(k_new,0.0_DP)
       return
    end if

    ! cks: start timer
    call timer_clock("palser_mano_kernel_optimise", 1)

    ! cks: initialisations
    call sparse_embed_create(k_old, k_new)
    call sparse_embed_create(ks, k_new, olap)
    call sparse_embed_create(ksk, ks, k_new)
    call sparse_embed_create(sks, olap, ks)
    call sparse_embed_create(sinv_ham, s_inv, ham)
    old_band_energy =huge(1.0_DP)
    quitearly = .false.
    if (present(num_iter)) then
       loc_num_iter = num_iter
    else
       loc_num_iter = 1
    end if

    ! cks: --------- PRINT PRELIMINARY INFO -----------------------
    if ( (pub_output_detail >= VERBOSE) .and. pub_on_root) then
       write(stdout,'(a)')' '
       write(stdout,'(a)')'======= &
            &Optimisation of K by Palser-Manolopoulos &
            &canonical purification ======== '
       write(stdout,'(a)')'         Iteration  |  Band structure energy&
            & |    c_n   '
    end if
    ! cks: -----END PRINT PRELIMINARY INFO -----------------------


    ! cks: ======== INITIALISE DENSKERN =====================

    ! cks: S^-1.H
    call sparse_embed_product(sinv_ham, s_inv, ham)

    ! cks: maximum orbital energy
    call sparse_embed_extremal_eigenvalue(sinv_ham, olap, max_en, delta_e_thresh)
    if (pub_debug_on_root) write(stdout,'(a,f24.14)') &
         'DEBUG: Max Hamiltonian eigenvalue:', max_en

    ! cks: -S^-1.H
    call sparse_embed_scale(sinv_ham, -1.0_DP)
    ! cks: minimum orbital energy
    call sparse_embed_extremal_eigenvalue(sinv_ham, olap, min_en, delta_e_thresh)
    min_en =-min_en
    if (pub_debug_on_root) write(stdout,'(a,f24.14)') &
         'DEBUG: Min Hamiltonian eigenvalue:', min_en

    ! cks: truncate the eigenenergy bounds precision to eval_prec digits as
    ! cks: sparse_extremal_eigenvalue results can become compiler-dependent
    ! cks: if more digits are used
    max_en =max_en -mod(max_en, eval_prec) + eval_prec
    min_en =min_en -mod(min_en, eval_prec) - eval_prec
    if (pub_debug_on_root) then
       write(stdout,'(a,f24.14)') 'DEBUG: Truncated max Ham eigenvalue:', max_en
       write(stdout,'(a,f24.14)') 'DEBUG: Truncated min Ham eigenvalue:', min_en
    end if

    ! cks: S^-1.H
    call sparse_embed_scale(sinv_ham, -1.0_DP)

    ! cks: average orbital energy
    call sparse_embed_trace(average_en, s_inv, ham)
    average_en = average_en / real(sparse_embed_num_rows(ham),kind=DP)
    if (pub_debug_on_root) write(stdout,'(a,f24.14)') &
         'DEBUG: Average orbital energy:', average_en

    lambda =min( real(n_occ, kind=DP)/(max_en -average_en), &
         (sparse_embed_num_rows(ham) -real(n_occ,kind=DP))/(average_en -min_en)  )

    ! cks: (lamda*mu + N_e)/N
    i_fac =(lambda*average_en +n_occ)/real(sparse_embed_num_rows(ham),kind=DP)

    ! cks: -lambda/N
    h_fac = -lambda/real(sparse_embed_num_rows(ham),kind=DP)

    if (pub_debug_on_root) then
       write(stdout,'(a,f24.14)') 'DEBUG: lambda:', lambda
       write(stdout,'(a,f24.14)') 'DEBUG:  i_fac:', i_fac
       write(stdout,'(a,f24.14)') 'DEBUG:  h_fac:', h_fac
    end if

    ! cks: S^-1.H.S^-1
    call sparse_embed_product(k_old, sinv_ham, s_inv )

    ! cks: -lambda/N*S^-1.H.S^-1
    call sparse_embed_scale(k_old, h_fac)

    ! cks: -lambda/N*S^-1.H.S^-1 + (lamda*mu + N_e)/N*S^-1
    call sparse_embed_axpy(k_old, s_inv, i_fac)
    ! cks: ==== END INITIALISE DENSKERN =====================

    c_n = -1.0_DP ! jd: to avoid printing garbage when quitearly=T in 1st iter.

    !########################################################
    ! cks: ######### PALSER-MANO ITERATIONS #################
    pm_loop: do iter=1,loc_num_iter

       call sparse_embed_trace(band_energy, k_old, ham)

       if(pub_debug) then
          if (iter == 1) then
             call sparse_embed_trace(debug_trace, s_inv, ham)
             if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: tr[Sinv*H]:',&
                  debug_trace
             call sparse_embed_trace(debug_trace, olap, ham)
             if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: tr[S*H]:',&
                  debug_trace
             call sparse_embed_trace(debug_trace, ham, ham)
             if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: tr[H*H]:',&
                  debug_trace
             call sparse_embed_trace(debug_trace, olap, olap)
             if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: tr[S*S]:',&
                  debug_trace
             call sparse_embed_trace(debug_trace, k_old, k_old)
             if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: tr[K*K]:',&
                  debug_trace
             call sparse_embed_trace(debug_trace, k_old, olap)
             if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: tr[K*S]:',&
                  debug_trace
             call sparse_embed_trace(debug_trace, s_inv, s_inv)
             if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: tr[Sinv*Sinv]:',&
                  debug_trace

          endif ! iter==1
       end if ! pub_debug

       ! cks: K.S
       call sparse_embed_product(ks, k_old, olap)

       ! cks: K.S.K
       call sparse_embed_product(ksk, ks, k_old)

       ! cks: S.K.S
       call sparse_embed_product(sks, olap, ks)

       ! cks: -S.K.S
       call sparse_embed_scale(sks, -1.0_DP)

       ! cks: (S-S.K.S)
       call sparse_embed_axpy(sks,olap, 1.0_DP)

       ! cks: tr[K.(S -S.K.S)]
       call sparse_embed_trace(c_denom, k_old, sks)

       ! cks: tr[K.S.K.(S -S.K.S)]
       call sparse_embed_trace(c_numer, ksk, sks)

       ! cks: Check criteria to quit based on machine precision or
       ! cks: density kernel truncation and quit if any are satisfied
       ! jd: Moved these two checks here, because they don't depend on c_n
       !     and should be tested before the c_denom test.
       if (band_energy > old_band_energy) then
          quitearly =.true.
          quitreason = 2
       end if
       delta_e_per_elec = (abs(band_energy - old_band_energy)/n_occ)
       if (delta_e_per_elec < delta_e_thresh) then
          quitearly =.true.
          quitreason = 3
       end if

       ! cks: print values for current K_old
       if ( (pub_output_detail >= VERBOSE) .and. pub_on_root) then
          write(stdout,'(t12,i5,tr5,f22.14)',advance='no') iter, band_energy
       end if

       if ( .not. quitearly .and. abs(c_denom) < tiny(1.0_DP)) then
          ! jd: ndmh's observation that zero denominator is OK if c_numer tiny
          if(abs(c_numer) < tiny(1.0_DP)) then
             quitearly = .true.
             quitreason = 1
          ! jd: corrolary to the above: zero denominator is OK if c_numer small
          !     and dE/elec small
          else if(abs(c_numer) < c_numer_small_enough .and. &
               delta_e_per_elec < delta_e_thresh_small_enough ) then
             quitearly = .true.
             quitreason = 6
          else
             write(error_message, &
                  '(a,e20.10,a,e20.10,a,e20.10,a,e20.10,a,e20.10,a,e20.10)') &
                  'Error in palser_mano_kernel_optimise: &
                  &c_denom is zero, but c_numer is not.'//CRLF//'c_numer: ', &
                  c_numer, CRLF//'c_denom: ', c_denom, &
                  CRLF//'last band_energy: ', band_energy, &
                  CRLF//'last old_band_energy: ', old_band_energy, &
                  CRLF//'last dE/elec: ', &
                  delta_e_per_elec, &
                  CRLF//'last (dE/elec)_threshold: ', delta_e_thresh, &
                  CRLF//'last (dE/elec)_threshold_small_enough: ', &
                  delta_e_thresh_small_enough

             call utils_abort(error_message)
          end if
       endif

       if(.not. quitearly) c_n = c_numer/c_denom

       if ( (pub_output_detail >= VERBOSE) .and. pub_on_root) then
          write(stdout,'(tr3,f10.7)') c_n
       end if

       if (c_n > 1.0_DP) then
          quitearly =.true.
          quitreason = 4
       end if
       if (c_n < 0.0_DP) then
          quitearly =.true.
          quitreason = 5
       end if

       if (quitearly) then
          if ((pub_output_detail >= VERBOSE) .and. pub_on_root) then
             write(stdout,'(a)', advance='no') &
                  'Termination criterion satisfied: '
             if(quitreason == 1) &
                  write(stdout,'(a)') 'c_numer and c_denom practically zero.'
             if(quitreason == 2) &
                  write(stdout,'(a)') 'band energy increased during last step.'
             if(quitreason == 3) &
                  write(stdout,'(a)') 'band energy change below threshold.'
             if(quitreason == 4) &
                  write(stdout,'(a)') 'c_n > 1'
             if(quitreason == 5) &
                  write(stdout,'(a)') 'c_n < 0'
             if(quitreason == 6) &
                  write(stdout,'(a)') 'c_numer very small and band energy &
                  &change close enough to threshold.'
          end if

          ! ndmh: copy k_old to k_new if we quit on first iteration
          if (iter==1) call sparse_embed_copy(k_new,k_old)
          exit pm_loop
       end if


       ! cks: (1+cn)*I - K.S
       call sparse_embed_scale(ks, -1.0_DP, 1.0_DP+c_n)

       ! cks: (1+cn)*K.S.K - K.S.K.S.K in k_new
       call sparse_embed_product(k_new, ks, ksk)


       ! cks: use different formula to obtain K_new depending on the value of c_n
       if ( c_n >= 0.5_DP) then
          scale_fac =1.0_DP/c_n
          call sparse_embed_scale(k_new, scale_fac )
       else
          scale_fac = 1.0_DP -2.0_DP*c_n
          call sparse_embed_axpy(k_new, k_old, scale_fac)
          scale_fac = 1.0_DP/(1.0_DP -c_n)
          call sparse_embed_scale(k_new, scale_fac )
       endif

       ! cks: store new density kernel and band energy
       call sparse_embed_copy(k_old, k_new)
       old_band_energy =band_energy

    enddo pm_loop
    ! cks: ###### END PALSER-MANO ITERATIONS #################
    !#########################################################

    ! cks: print warning if convergence criteria not satisfied
    if (pub_on_root .and. .not. quitearly) write(stdout,*) &
         'WARNING: Maximum Palser-Manolopoulos iterations exceeded.'

    ! ndmh: Check eigenvalue bounds of ks

    ! ndmh: create temporary matrix for KSKS structures
    call sparse_embed_create(ksks,ks,ks)
    ! ndmh: find upper bound of occupancies
    call sparse_embed_product(ks, k_new, olap)
    call sparse_embed_extremal_eigenvalue(ks,olap,ubound)
    ! ndmh: find lower bound of occupancies
    call sparse_embed_scale(ks,-1.0_DP,1.0_DP)
    call sparse_embed_extremal_eigenvalue(ks,olap,lbound)
    lbound = 1.0_DP - lbound
    ! ndmh: find closest eigenvalue to 0.5
    call sparse_embed_product(ks, k_new, olap)
    call sparse_embed_scale(ks,-1.0_DP,0.5_DP)
    call sparse_embed_product(ksks,ks,ks)
    call sparse_embed_scale(ksks,-1.0_DP,0.0_DP)
    call sparse_embed_extremal_eigenvalue(ksks,olap,mid_occ)
    mid_occ = sqrt(abs(mid_occ)) + 0.5_DP

    ! ndmh: warn if eigenvalue suggests degeneracy at eF
    if (pub_on_root.and.(mid_occ>0.15_DP).and.(mid_occ<0.85_DP)) then
       if (.not.pub_cond_calculate) then
          write(stdout,'(/a/a)') 'WARNING in palser_mano_kernel_optimise: &
               &Possible degeneracy at fermi level'
       else
          write(stdout,'(/a/a)') 'WARNING in palser_mano_kernel_optimise: &
               &Possible degeneracy at topmost state of conduction kernel'
       end if
       ! ndmh: if the calculation was reckoned to be converged,
       ! ndmh: suggest diagonalisation instead.
       if (quitearly) write(stdout,'(a)') 'Diagonalisation &
            &(maxit_palser_mano=-1) may be required.'
    end if

    if (pub_on_root.and.(pub_output_detail >= VERBOSE)) then
       write(stdout,'(a,f20.14,a,f20.14,a)') &
            'Final occupancy bounds: [',lbound,':',ubound,']'
       write(stdout,'(a,f20.14)') &
            'Closest occupancy to 0.5: ',mid_occ
       write(stdout,'(a)')'===================================&
            &============================================='
    end if

    call sparse_embed_destroy(ksks)

    ! cks: deallocate SPAM3 memory
    call sparse_embed_destroy(sinv_ham)
    call sparse_embed_destroy(sks)
    call sparse_embed_destroy(ksk)
    call sparse_embed_destroy(ks)
    call sparse_embed_destroy(k_old)

    ! cks: stop timer
    call timer_clock("palser_mano_kernel_optimise", 2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
         &Leaving palser_mano_kernel_optimise'

    ! Flush output
    call services_flush

    return

  end subroutine palser_mano_kernel_optimise


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module palser_mano
