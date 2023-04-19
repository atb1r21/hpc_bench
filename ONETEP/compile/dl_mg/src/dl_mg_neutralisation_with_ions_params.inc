!! Parameters that describe the source charge neutralisation

integer, parameter :: dl_mg_neutralise_none                  = 32100

!! adds a background of opposite charge in the computational cell
!! can be used for Poisson and Poison-Boltzmann equations
integer, parameter :: dl_mg_neutralise_with_jellium_uniform  = 33100

!! adds the neutralising background only in the ion accessible region
!! to be used for Poisson-Boltzmann equation
integer, parameter :: dl_mg_neutralise_with_jellium_vacc     = 34100

!! neutralises the source charge with a user prescribed shift in
!! ion concentrations
integer, parameter :: dl_mg_neutralise_with_ions_fixed       = 35100

!! uses an asymptotic correspondence to a open system to find
!! the shift parameters for the concentrations
integer, parameter :: dl_mg_neutralise_with_ions_auto        = 36100

!! uses a linearised version of the above method
integer, parameter :: dl_mg_neutralise_with_ions_auto_linear = 37100