!
!> \brief  General parameters
!! Lucian Anton February 2013
module dl_mg_params
  implicit none

  integer, parameter :: single_precision = kind(0.0)
  integer, parameter :: double_precision = kind(0.d0)
  integer, parameter :: wp = double_precision ! working precision

  real(kind=wp), parameter :: pi = 4.0_wp*atan(1.0_wp)
  real(kind=wp), parameter :: elcharge = 4.80320427e-10_wp !< electron charge in CGS !
  real(kind=wp), parameter :: hartree = 3.1577464761980001721e5_wp !< hartree energy divided by kB in Gauss units,
                                                                   !< i.e. \f$ e^2/(r_B * K_B ) \f$
  real(kind=wp), parameter :: kboltz = 1.3806504e-16_wp !< Boltzmann constant in CGS

  
!> tags for equations types, used to select numerical kernels
  integer, parameter ::  EQ_POISSON = 10, &
       EQ_LINEAR_PBE = 20, EQ_LINEAR_PBE_NOSTERIC = 21, &
        EQ_LINEAR_PBE_POT = 22, &
       EQ_PBE_FAS = 30, EQ_PBE_FAS_NOSTERIC = 31, &
       EQ_PBE_NEWTON = 40, EQ_PBE_NEWTON_NOSTERIC = 41

  integer, parameter :: dl_mg_half_weight_restriction = 101
  integer, parameter :: dl_mg_full_weight_restriction = 102
  integer, parameter :: dl_mg_injection_restriction = 103


  !> parameters to be accessed from dl_mg.mod
  include  "dl_mg_common_params.inc"

  !> parameters for ion neutralisation of the solute.
  !! relevant when using PBE with PBC
  include "dl_mg_neutralisation_with_ions_params.inc"

end module dl_mg_params
