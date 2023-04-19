! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by Andrea Greco
!
!   January 2016
!
!   Subsequent additions by JM Escartin.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module k_points

  use constants, only: DP
  use geometry, only: POINT

  implicit none

  private

  ! Public subroutines
  public :: kpoints_generate_unique_list
  public :: kpoints_init

  ! agrecokpt: generate unique k-point from list or mesh
  interface kpoints_generate_unique_list
     module procedure kpoints_generate_unique_from_list
     module procedure kpoints_generate_unique_from_mesh
  end interface

  ! agreco: this structure defines quantities which are all associated with
  !      a particular k-point
  type KPOINT
     ! agreco: cartesian coordinates of k-point
     type(POINT)      :: centre
     ! agreco: weight of this k-point
     real(kind=DP) :: weight
     ! agreco: label for particular points in BZ, i.e. Gamma, L, M, R, ...
     character(len=5) :: label

  end type KPOINT

  ! agreco: public type definitions
  public :: KPOINT

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine kpoints_generate_unique_from_list(unique_kpoint_list, &
             num_kpoints, orig_kpoint_list, cell)

    !=========================================================================!
    ! This subroutine generates an array of unique k-points (type KPOINT) to  !
    ! be used in different routines for BZ sampling.                          !
    !                                                                         !
    ! Arguments:                                                              !
    !     orig_kpoint_list: the list of kpoints read from the input file,     !
    !     in fractional coordinates.                                          !
    !     Currently the format is orig_kpoint_list(:,kpt_index) =             !
    !     (kpt_x, kpt_y, kpt_z, weight)                                       !
    !                                                                         !
    ! Written by Andrea Greco on 06/01/2016.                                  !
    ! TO-DO: write methods to determine unique k-points from original list    !
    ! depending on the crystal structure; also include possibility of using   !
    ! a MP k-point grid instead.                                              !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout
    use geometry, only: operator(*), operator(+)
    use rundat, only: pub_gamma_point_only, pub_num_kpoints_temp
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort, utils_alloc_check, utils_banner

    implicit none

    ! Arguments
    type(KPOINT), allocatable, intent(inout) :: unique_kpoint_list(:)
    integer, intent(out) :: num_kpoints
    real(kind=DP), intent(in) :: orig_kpoint_list(:,:)
    type(CELL_INFO), intent(in) :: cell

    ! Local arguments
    integer :: ik1, ik2
    integer :: ierr
    real(kind=DP) :: kpoint_weight_sum

    ! Check arguments
    ! N.B: this check needs to be changed if we allow different formats for
    ! orig_kpoint_list
    if (size(orig_kpoint_list,2) < 1) then
       call utils_abort(&
            'Error creating unique k-point list (k_points_mod.F90):&
            & original k-point list has no k-points')
    end if
    if (size(orig_kpoint_list,1) < 4) then
       call utils_abort(&
            'Error creating unique k-point list (k_points_mod.F90):&
            & either the k-point weight or one or more coordinates &
            & are missing from the original k-point list')
    end if

    num_kpoints = size(orig_kpoint_list,2)

    ! Check for duplicate k-points in original list
    do ik1=2,num_kpoints
       do ik2=1,ik1-1
          if (all(orig_kpoint_list(1:3,ik1) == orig_kpoint_list(1:3,ik2))) then
             call utils_abort(&
                'Error creating unique k-point list (k_points_mod.F90):&
                & duplicate k-point in original k-point list')
          end if
       end do
    end do

    ! if only one k-point is specified and it is Gamma, no need to allocate
    ! unique k-point list, set pub_gamma_point_only to true
    if ((num_kpoints==1) .and. all(orig_kpoint_list(1:3,1)==0.0_DP)) then
       pub_gamma_point_only = .true.
    else
       if (.not.allocated(unique_kpoint_list)) then
          allocate(unique_kpoint_list(num_kpoints),stat=ierr)
          call utils_alloc_check('generate_unique_kpoints_list', &
               'unique_kpoint_list', ierr)
          ! TO-DO: set k-points label look-up k-point table
          ! agrecokpt: convert from fractional to cartesian coordinates
          do ik1=1,num_kpoints
             unique_kpoint_list(ik1)%centre = orig_kpoint_list(1,ik1) * &
                cell%b1 + orig_kpoint_list(2,ik1) * cell%b2 + &
                orig_kpoint_list(3,ik1) * cell%b3
             unique_kpoint_list(ik1)%weight = orig_kpoint_list(4,ik1)
          end do
          ! agreco: normalize k-points weight
          kpoint_weight_sum = sum(unique_kpoint_list(:)%weight)
          unique_kpoint_list(:)%weight = &
             unique_kpoint_list(:)%weight/kpoint_weight_sum
          pub_gamma_point_only = .false.
       else if (size(unique_kpoint_list)/=num_kpoints) then
          call utils_abort('Error in generate_unique_kpoints_list (k_points) &
            &size of unique_kpoint_list array does not match number of k-points')
       end if
    end if

    pub_num_kpoints_temp = num_kpoints

    ! agrecokpt: print summary of k-points from list with normalized weights
    if (pub_on_root) then
       write(stdout,'(/a)') utils_banner('=','K-points summary')
       write(stdout,'(a)') &
          'Number |             Fractional coordinates          |    Weight'
       do ik1=1,num_kpoints
          write(stdout,'(i3,4f16.10)') ik1, &
               orig_kpoint_list(1:3,ik1), unique_kpoint_list(ik1)%weight
       end do
       write(stdout,'(a/)') repeat('=',80)
    end if

  end subroutine kpoints_generate_unique_from_list

!-------------------------------------------------------------------------------

  subroutine kpoints_generate_unique_from_mesh(unique_kpoint_list, &
             num_kpoints, mesh, is_shift, is_time_reversal, cell, elements, nat)

    !=========================================================================!
    ! This subroutine generates an array of unique k-points (type KPOINT) to  !
    ! be used in different routines for BZ sampling. The unique k-points are  !
    ! generated from a mesh (possibly shifted), according to the symmetries   !
    ! of the crystal structure.                                               !
    !                                                                         !
    ! Arguments:                                                              !
    !     mesh: the k-point grid read from the input file                     !
    !     is_shift: the shift of the k-point grid wrt the origin              !
    !     cell: CELL_INFO type containing info about lattice vectors          !
    !     elements: array of ELEMENT type with info on the atoms              !
    !     nat: number of atoms in the simulation cell                         !
    !                                                                         !
    ! Written by Andrea Greco on 30/11/2016.                                  !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: TWO_PI, stdout
    use geometry, only: operator(.DOT.), operator(*), operator(+)
    use ion, only: ELEMENT
    use rundat, only: pub_gamma_point_only, pub_num_kpoints_temp
    use simulation_cell, only: CELL_INFO
    use spglib_f08
    use utils, only: utils_abort, utils_alloc_check, utils_banner

    implicit none

    ! Arguments
    type(KPOINT), allocatable, intent(inout) :: unique_kpoint_list(:)
    integer, intent(out) :: num_kpoints
    integer, intent(in), dimension(3) :: mesh, is_shift
    integer, intent(in) :: is_time_reversal, nat
    type(CELL_INFO), intent(in) :: cell
    type(ELEMENT), intent(in), dimension(nat) :: elements

    ! Local arguments
    integer :: iat, ik, counter, weight, ierr
    real(kind=DP), dimension(3, 3) :: lattice
    real(kind=DP), dimension(3, nat) :: posfrac
    real(kind=DP), parameter :: symprec = 1D-6
    integer, dimension(3, mesh(1)*mesh(2)*mesh(3)) :: grid_point
    integer, dimension(mesh(1)*mesh(2)*mesh(3)) :: map
    integer, dimension(nat) :: atom_types
    real(kind=DP), parameter :: recip_twopi = 1.0_DP / TWO_PI
    real(kind=DP), dimension(3) :: kptfrac
    real(kind=DP) :: weight_frac

    ! Check arguments
    if (any(mesh==0)) then
       call utils_abort(&
            'Error creating unique k-point list (k_points_mod.F90):&
            & k-point grid has no k-points along at least one dimension')
    end if

    ! if only one k-point is specified and it is Gamma, no need to allocate
    ! unique k-point list, set pub_gamma_point_only to true
    if (all(mesh==1) .and. all(is_shift==0)) then
       pub_gamma_point_only = .true.
    else

#ifdef SPG
       ! agreco: build array with cell vectors (Cartesian)
       lattice(1,:) = (/ cell%a1%x, cell%a1%y, cell%a1%z /)
       lattice(2,:) = (/ cell%a2%x, cell%a2%y, cell%a2%z /)
       lattice(3,:) = (/ cell%a3%x, cell%a3%y, cell%a3%z /)

       ! agreco: build fractional positions and atomic species
       do iat = 1,nat
          posfrac(1,iat) = elements(iat)%centre .DOT. cell%b1
          posfrac(2,iat) = elements(iat)%centre .DOT. cell%b2
          posfrac(3,iat) = elements(iat)%centre .DOT. cell%b3
          atom_types(iat) = elements(iat)%species_number
       end do

       posfrac = posfrac * recip_twopi

       posfrac = modulo(posfrac,1.0_DP)

       num_kpoints = spg_get_ir_reciprocal_mesh( grid_point, map, &
          mesh, is_shift, is_time_reversal, lattice, posfrac, &
          atom_types, nat, symprec )

       if (.not.allocated(unique_kpoint_list)) then
          allocate(unique_kpoint_list(num_kpoints),stat=ierr)
          call utils_alloc_check('generate_unique_kpoints_list', &
               'unique_kpoint_list', ierr)

          counter = 0

          if (pub_on_root) then
             write(stdout,'(/a)') utils_banner('=','K-points summary')
             write(stdout,'(a)') &
                'Number |             Fractional coordinates          |    Weight'
          end if

          do ik = 1, product(mesh)
             if (ik-1 == map(ik)) then
             ! Ad-hoc and intuitive implementation of weight. Not optimal for very large size
                weight = count(map == ik-1)
                counter = counter + 1

                kptfrac = real(is_shift + 2*grid_point(:, ik), kind=DP)/real(2*mesh, kind=DP)
                weight_frac = real(weight, kind=DP)/real(product(mesh), kind=DP)

                if (pub_on_root) then
                   write(stdout,'(i3,4f16.10)') counter, kptfrac(:), weight_frac
                end if

                unique_kpoint_list(ik)%centre = kptfrac(1) * &
                   cell%b1 + kptfrac(2) * cell%b2 + kptfrac(3) * cell%b3

                unique_kpoint_list(ik)%weight = weight_frac

             end if
          end do

          if (pub_on_root) write(stdout,'(a/)') repeat('=',80)

          pub_gamma_point_only = .false.

       else if (size(unique_kpoint_list)/=num_kpoints) then
          call utils_abort('Error in generate_unique_kpoints_list (k_points) &
            &size of unique_kpoint_list array does not match number of k-points')
       end if
#else
       num_kpoints = 0 !jme: avoid compilation warning
       call utils_abort(&
          'Error creating unique k-point list (k_points_mod.F90):&
           & mesh specification requires SPG library to work. ')
#endif

    end if

    pub_num_kpoints_temp = num_kpoints

  end subroutine kpoints_generate_unique_from_mesh

!-------------------------------------------------------------------------------

  subroutine kpoints_init(kpt_cart, cell)

    !==========================================================================!
    ! k-points initialisations and checks                                      !
    !                                                                          !
    ! Created by JM Escartin, from existing code in onetep.F90 (August 2016)   !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout
    use simulation_cell, only: CELL_INFO
    use rundat, only: pub_single_kpt, pub_gamma_point_only, pub_cmplx_ngwfs
    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: kpt_cart(3)
    type(CELL_INFO), intent(in) :: cell

    ! agrecokpt: override pub_cmplx_ngwfs value when using
    ! non Gamma k-point only; uncomment this when everything
    ! is working properly
    ! if (.not.pub_gamma_point_only) then
    !    pub_cmplx_ngwfs = .true.
    !    if (pub_on_root) write(stdout,'(a)') 'Non-Gamma &
    !       &calculation requested: complex NGWFs will be used'
    ! end if
    ! agrecokpt: currently using single kpt for testing, to be removed later
    ! use complex NGWFs when using a single non-Gamma kpt
    if (any(pub_single_kpt.ne.0.0_DP)) then
       pub_cmplx_ngwfs = .true.
       ! convert from fractional to cartesian coordinates
       kpt_cart(1) = pub_single_kpt(1) * cell%b1%x + &
          pub_single_kpt(2) * cell%b2%x + &
          pub_single_kpt(3) * cell%b3%x
       kpt_cart(2) = pub_single_kpt(1) * cell%b1%y + &
          pub_single_kpt(2) * cell%b2%y + &
          pub_single_kpt(3) * cell%b3%y
       kpt_cart(3) = pub_single_kpt(1) * cell%b1%z + &
          pub_single_kpt(2) * cell%b2%z + &
          pub_single_kpt(3) * cell%b3%z

       !pub_gamma_point_only = .false.
       if (pub_on_root) write(stdout,'(a)') 'Non-Gamma &
          &calculation requested: complex NGWFs will be used'
    ! do not set pub_cmplx_ngwfs to .false. if Gamma point only,
    ! in order to allow use of complex NGWFs with Gamma point
    ! if the user wants to do so
    else
       kpt_cart(:) = 0.0_DP
    end if
  end subroutine kpoints_init

end module k_points
