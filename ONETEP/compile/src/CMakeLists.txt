include_directories("${CMAKE_BINARY_DIR}/src")
include_directories("${CMAKE_BINARY_DIR}/dl_mg")

set(onetep_SRCS
alloc_checker_mod.F90
anharmonic_mod.F90
atom_mod.F90
augmentation_mod.F90
bandstructure_mod.F90
basis_mod.F90
bibliography_mod.F90
cdft_intermediate_cg_mod.F90
cell_grid_mod.F90
chebyshev_rep_mod.F90
classical_pot_mod.F90
comms_mod.F90
conduction_mod.F90
conduction_properties_mod.F90
constants_mod.F90
couplings_mod.F90
cutoff_coulomb_mod.F90
datatypes_mod.F90
ddec_mod.F90
dense_mod.F90
density_init_mod.F90
density_mod.F90
dftb_mod.F90
dftb_parameters_mod.F90
dftb_type_mod.F90
dma_mod.F90
dmft_mod.F90
eda_driver_supermol_mod.F90
eda_main_mod.F90
eels_mod.F90
eigenstates_mod.F90
eld_mod.F90
electronic_history_mod.F90
electronic_init_mod.F90
electronic_mod.F90
energy_and_force_mod.F90
ensemble_dftb_mod.F90
ensemble_dft_mod.F90
ensemble_dft_type_mod.F90
enthalpy_mod.F90
esdf_key_mod.F90
esdf_mod.F90
etrans_mod.F90
ewald_mod.F90
fft_box_mod.F90
file_handling_mod.F90
finite_differences_mod.F90
forces_mod.F90
forcetest_mod.F90
fourier_gpu_wrapper_mod.F90
fourier_mod.F90
fragment_data_mod.F90
fragment_matrix_ops_mod.F90
fragment_scfmi_mod.F90
function_basis_mod.F90
function_ops_mod.F90
gaunt_coeff_mod.F90
gaussian_ops_mod.F90
geometry_deloc_algor_mod.F90
geometry_deloc_utils_mod.F90
geometry_mod.F90
geometry_optimiser_mod.F90
greenf_mod.F90
hamiltonian_mod.F90
hartree_mod.F90
hash_table_mod.F90
hf_exchange_mod.F90
hubbard_build_mod.F90
hubbard_init_mod.F90
image_comms_mod.F90
integrals_mod.F90
ion_mod.F90
is_poisson_mod.F90
is_smeared_ions_mod.F90
is_solvation_boltzmann_mod.F90
is_solvation_mod.F90
ke_density_mod.F90
kernel_diis_mod.F90
kernel_mod.F90
kinetic_mod.F90
k_points_mod.F90
license_mod.F90
linalg_mod.F90
linear_response_mod.F90
lnv_mod.F90
lr_phonons_mod.F90
lr_tddft_mod.F90
lr_tddft_rpa_mod.F90
lr_tddft_utils_mod.F90
md_ions_mod.F90
md_mod.F90
md_thermostat_mod.F90
mermin_mod.F90
model_type_mod.F90
multigrid_methods_mod.F90
multipole_ops_mod.F90
neb_mod.F90
neighbour_list_mod.F90
ngwf_cg_mod.F90
ngwf_data_mod.F90
ngwf_gradient_mod.F90
ngwf_representation_mod.F90
ngwfs_mod.F90
ngwfs_radial_mod.F90
nho_mod.F90
nlxc_mod.F90
nonsc_forces_mod.F90
npa_mod.F90
openbc_locps_mod.F90
optics_mod.F90
palser_mano_mod.F90
parallel_strategy_mod.F90
paw_mod.F90
paw_shape_mod.F90
paw_xc_mod.F90
pbc_corrections_mod.F90
penalty_mod.F90
phonon_mod.F90
polarisable_embedding_mod.F90
population_mod.F90
potential_mod.F90
ppd_ops_mod.F90
ppd_strategy_mod.F90
projectors_mod.F90
properties_mod.F90
pseudopotentials_mod.F90
qnto_mod.F90
remote_mod.F90
restart_mod.F90
rundat_blocks_mod.F90
rundat_mod.F90
scissor_mod.F90
services_mod.F90
simulation_cell_mod.F90
smearing_operator_mod.F90
sparse_algor_mod.F90
sparse_array_mod.F90
sparse_embed_mod.F90
sparse_initialise_mod.F90
sparse_mod.F90
spglib_f08_mod.F90
spherical_wave_mod.F90
sph_harm_rotation_mod.F90
sw_expansion_mod.F90
sw_expansion_type_mod.F90
sw_resolution_of_identity_mod.F90
tddft_mod.F90
timer_mod.F90
tssearch_algor_mod.F90
tssearch_mod.F90
tssearch_utils_mod.F90
unit_test_mod.F90
utils_mod.F90
vdwcorrection_mod.F90
vdw_df_kernel_mod.F90
visual_mod.F90
xc_b97mv_mod.F90
xc_funcs_mod.F90
xc_mod.F90)

add_library(onetep ${onetep_SRCS})
if (WITH_SCALAPACK)
   set_property(TARGET onetep APPEND PROPERTY
     COMPILE_DEFINITIONS "SCALAPACK")
endif()

set_property(TARGET onetep APPEND PROPERTY
  COMPILE_DEFINITIONS ${FFT})

set_property(TARGET onetep APPEND PROPERTY
  COMPILE_DEFINITIONS ParallelFFT)

if (WITH_MPI)
  set_property(TARGET onetep APPEND PROPERTY
    COMPILE_DEFINITIONS "MPI")
endif()

if (WITH_DLMG)
  add_dependencies(onetep dl_mg)
  set_property(TARGET onetep APPEND PROPERTY
    COMPILE_DEFINITIONS "HAVE_DL_MG")
endif()

if (NOT WITH_OPENMP AND ${FFT} STREQUAL "FFTW3")
  set_property(TARGET onetep APPEND PROPERTY
    COMPILE_DEFINITIONS "FFTW3_NO_OMP")
endif()

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set_property(TARGET onetep APPEND PROPERTY
    COMPILE_DEFINITIONS "GFORTRAN")
endif()

link_directories (${LIBRARY_OUTPUT_PATH})
if(WITH_DLMG)
  target_link_libraries(onetep PRIVATE dl_mg  )
endif()

foreach(p IN LISTS EXTRA_PATHS)
  if (EXISTS ${p})
    target_link_directories(onetep PRIVATE ${p})
    message(STATUS "add linking path ${p}")
  endif()
endforeach()

foreach(p IN LISTS EXTRA_INCS)
  if (EXISTS ${p})
    message(STATUS "add include path ${p}")
    target_include_directories(onetep PRIVATE ${p})
  endif()
endforeach()

foreach(p IN LISTS EXTRA_DEFS)
  set_property(TARGET onetep APPEND PROPERTY
    COMPILE_DEFINITIONS ${p})
endforeach()

if (WITH_MKL)
  target_include_directories(onetep PRIVATE ${MKL_INCLUDE})
endif()

target_link_libraries(onetep PRIVATE ${LIBS} )
add_executable(${target_name} onetep.F90)
target_link_libraries(${target_name} PRIVATE onetep ${LIBS})

install(TARGETS ${target_name} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
