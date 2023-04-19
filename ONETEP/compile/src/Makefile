# Second-level makefile for ONETEP
#
# This version for multiple architectures by Peter Haynes
#
# System-specific changes should be made to the appropriate
#  conf.arch files
#

EXECDIR = $(ROOTDIR)/bin
INCPATH = -I$(ROOTDIR)/include
LIBPATH = -L$(ROOTDIR)/lib/$(ARCH)

UTILDIR = $(ROOTDIR)/utils

include $(ROOTDIR)/config/conf.$(ARCH)

EXTFLAGS = $(OPTFLAGS)

.SUFFIXES:
.SUFFIXES: .o .F90

#
# List of source files
#
SOURCES = alloc_checker_mod.F90 anharmonic_mod.F90 atom_mod.F90 \
          augmentation_mod.F90 \
          bandstructure_mod.F90 basis_mod.F90 bibliography_mod.F90 \
          cdft_intermediate_cg_mod.F90 cell_grid_mod.F90 chebyshev_rep_mod.F90 \
          classical_pot_mod.F90 comms_mod.F90 conduction_mod.F90 \
          conduction_properties_mod.F90 constants_mod.F90 \
          couplings_mod.F90 cutoff_coulomb_mod.F90 \
          datatypes_mod.F90 ddec_mod.F90 \
          dense_mod.F90 density_mod.F90 density_init_mod.F90 dftb_mod.F90 \
          dftb_parameters_mod.F90 dma_mod.F90 dmft_mod.F90 \
          eels_mod.F90 eigenstates_mod.F90 electronic_history_mod.F90 \
          eld_mod.F90 electronic_mod.F90 electronic_init_mod.F90 \
          eda_driver_supermol_mod.F90 eda_main_mod.F90 \
          energy_and_force_mod.F90 enthalpy_mod.F90 ensemble_dft_mod.F90 \
          ensemble_dft_type_mod.F90 ensemble_dftb_mod.F90 \
          esdf_key_mod.F90 esdf_mod.F90 etrans_mod.F90 ewald_mod.F90 \
          fft_box_mod.F90 file_handling_mod.F90 finite_differences_mod.F90 \
          forces_mod.F90 forcetest_mod.F90 fourier_mod.F90 \
          fourier_gpu_wrapper_mod.F90 \
          fragment_data_mod.F90 fragment_matrix_ops_mod.F90 \
          fragment_scfmi_mod.F90 function_basis_mod.F90 \
          function_ops_mod.F90 gaunt_coeff_mod.F90 gaussian_ops_mod.F90 \
          geometry_deloc_algor_mod.F90 geometry_deloc_utils_mod.F90 geometry_mod.F90 \
          geometry_optimiser_mod.F90 greenf_mod.F90 hamiltonian_mod.F90 \
          hartree_mod.F90 hash_table_mod.F90 hf_exchange_mod.F90 \
          hubbard_build_mod.F90 hubbard_init_mod.F90 image_comms_mod.F90 \
          integrals_mod.F90 ion_mod.F90 is_poisson_mod.F90 is_smeared_ions_mod.F90 \
          is_solvation_mod.F90 is_solvation_boltzmann_mod.F90 \
          k_points_mod.F90 ke_density_mod.F90 \
          kernel_diis_mod.F90 kernel_mod.F90 \
          kinetic_mod.F90 linalg_mod.F90 linear_response_mod.F90 \
          lnv_mod.F90 lr_phonons_mod.F90 lr_tddft_mod.F90 \
          lr_tddft_rpa_mod.F90 lr_tddft_utils_mod.F90 md_mod.F90 \
          md_ions_mod.F90 md_thermostat_mod.F90 mermin_mod.F90 model_type_mod.F90 \
          multigrid_methods_mod.F90 multipole_ops_mod.F90 neb_mod.F90 \
          neighbour_list_mod.F90 ngwf_cg_mod.F90 ngwf_data_mod.F90 \
          ngwf_gradient_mod.F90 \
          ngwf_representation_mod.F90 ngwfs_mod.F90 ngwfs_radial_mod.F90 \
          nho_mod.F90 npa_mod.F90 nlxc_mod.F90 pseudopotentials_mod.F90 \
          nonsc_forces_mod.F90 onetep.F90 openbc_locps_mod.F90 optics_mod.F90 \
          palser_mano_mod.F90 paw_mod.F90 paw_shape_mod.F90 paw_xc_mod.F90 \
          parallel_strategy_mod.F90 pbc_corrections_mod.F90 \
          penalty_mod.F90 phonon_mod.F90 potential_mod.F90 ppd_ops_mod.F90 \
          ppd_strategy_mod.F90 polarisable_embedding_mod.F90 \
          population_mod.F90 projectors_mod.F90 properties_mod.F90 \
          qnto_mod.F90 remote_mod.F90 restart_mod.F90 rundat_blocks_mod.F90 \
          rundat_mod.F90 scissor_mod.F90 services_mod.F90 simulation_cell_mod.F90 \
          smearing_operator_mod.F90 sparse_mod.F90 sparse_algor_mod.F90 \
          sparse_array_mod.F90 sparse_embed_mod.F90 sparse_initialise_mod.F90 \
          spglib_f08_mod.F90 sph_harm_rotation_mod.F90 \
          spherical_wave_mod.F90 sw_expansion_mod.F90 sw_expansion_type_mod.F90 \
          sw_resolution_of_identity_mod.F90 \
          tddft_mod.F90 timer_mod.F90 \
          tssearch_algor_mod.F90 tssearch_mod.F90 tssearch_utils_mod.F90 \
          unit_test_mod.F90 utils_mod.F90 vdwcorrection_mod.F90 \
          vdw_df_kernel_mod.F90 visual_mod.F90 xc_mod.F90 \
          xc_b97mv_mod.F90 xc_funcs_mod.F90

MODULE_SOURCES = $(filter %_mod.F90,$(SOURCES))
MODULES = $(MODULE_SOURCES:_mod.F90=.mod)
# Special case: fourier_gpu_wrapper_mod.F90 contains an inconsistent name...
MODULES := $(subst fourier_gpu_wrapper.mod,fourier_gpu_wrapper_mod.mod,$(MODULES))
OBJECTS = $(SOURCES:.F90=.o)

# Sources to be parsed by f90wrap for creation of python wrappers (currently
# a very minimal list - will expand over time)
WRAPPED_SOURCES = simulation_cell_mod.F90 constants_mod.F90 geometry_mod.F90
PP_WRAPPED_SOURCES = $(WRAPPED_SOURCES:.F90=.f90)

# If MULTIGRID if set to 'dl_mg', include obj/$(ARCH)/dl_mg*.o in the list of dependencies for onetep.
#
ifeq ($(strip $(MULTIGRID)),dl_mg)
   FFLAGS := $(FFLAGS) $(PP_PREF)-DHAVE_DL_MG
   MULTIGRID_IMPLEMENTATION := $(wildcard $(ROOTDIR)/obj/$(ARCH)/dl_mg*.o)
   OBJECTS += $(MULTIGRID_IMPLEMENTATION)
   $(info dl_mg files $(MULTIGRID_IMPLEMENTATION))
endif

# --------------------------------------------------------------------------
# End section for the DL_MG multigrid solver
# --------------------------------------------------------------------------

ifeq ($(strip $(BUILD_ELPA)),yes)
ifneq (,$(findstring -DSCALAPACK,$(FFLAGS)))
   FFLAGS := $(FFLAGS) $(PP_PREF)-DELPA2
   ELPA_OBJ := $(wildcard $(ROOTDIR)/obj/$(ARCH)/elpa*.o)
   OBJECTS += $(ELPA_OBJ)
#   $(info dl_mg files $(MULTIGRID_IMPLEMENTATION))
endif
endif

#
# General rules
#
F90_BUILD_CMD = $(F90) $(FFLAGS) $(ACCFLAGS) $(EXTFLAGS) $(INCPATH) $(CUDAINCLUDE) $(ACCLIBS) -c $<
ifeq ($(UPPERCASE_MODULE),yes)
CP_MOD = cp $$(echo $@ | tr '[a-z]' '[A-Z]' | sed 's/\.MOD/\.mod/g') $@
else
CP_MOD =
endif
%.o: %.F90
	$(F90_BUILD_CMD)
%.mod: %_mod.F90
	$(F90_BUILD_CMD)
	$(CP_MOD)
%.mod: %.F90
	$(F90_BUILD_CMD)
	$(CP_MOD)

# Preprocess to get .f90 for f90wrap
%.f90: %.F90
	$(F90_BUILD_CMD) -E > $@

#
# Specific rules
#

default:
	@echo 'No target specified!'

onetep: $(MODULES) $(MULTIGRID_IMPLEMENTATION) onetep.o
	$(F90) $(FFLAGS) $(ACCFLAGS) $(EXTFLAGS) -o $(EXECDIR)/onetep.$(ARCH) $(OBJECTS) $(LIBPATH) $(LIBS) $(CUDAINCLUDE) $(ACCLIBS)

debug:: EXTFLAGS = $(PP_PREF)-DDEBUG $(DEBUGFLAGS)
debug:: onetep

timing:: EXTFLAGS = $(EXTFLAGS) -DTIMINGS
timing:: onetep

profile:: EXTFLAGS = $(OPTFLAGS) $(PROFILEFLAGS)
profile:: onetep

warnings:: EXTFLAGS = $(OPTFLAGS) $(WARNINGFLAGS)
warnings:: onetep

libonetep: $(OBJECTS)

python: libonetep $(PP_WRAPPED_SOURCES)
	f90wrap $(PP_WRAPPED_SOURCES) -v -k ../../src/kind_map -m onetep
	# The following line is a work-in-progress and can certainly be made more
	# general to support other compilers. Currently works with Intel+MPI+OpenMP
	F90=$(shell echo $(F90) | cut -d ' ' -f 1) LDFLAGS="$(LIBS) $(shell mpif90 --show | cut -d ' ' -f 2-) -openmp"  f2py-f90wrap -c -m _onetep f90wrap*.f90 -L. -lonetep -I.

#
# Dependencies
#

alloc_checker.mod: alloc_checker_mod.F90 comms.mod constants.mod utils.mod

anharmonic.mod: anharmonic_mod.F90 bibliography.mod constants.mod comms.mod \
                dense.mod geometry.mod ion.mod linalg.mod rundat.mod \
                simulation_cell.mod utils.mod

atom.mod: atom_mod.F90 constants.mod linalg.mod ion.mod pseudopotentials.mod \
          paw.mod rundat.mod services.mod \
          spherical_wave.mod utils.mod xc.mod

augmentation.mod: augmentation_mod.F90 cell_grid.mod comms.mod constants.mod \
                  datatypes.mod function_basis.mod fft_box.mod fourier.mod \
                  geometry.mod model_type.mod parallel_strategy.mod paw.mod \
                  projectors.mod pseudopotentials.mod rundat.mod \
                  simulation_cell.mod sparse.mod spherical_wave.mod timer.mod \
                  utils.mod xc.mod \

bandstructure.mod: bandstructure_mod.F90 augmentation.mod \
                   comms.mod constants.mod dense.mod function_basis.mod \
                   geometry.mod hamiltonian.mod integrals.mod linalg.mod \
                   model_type.mod ngwf_representation.mod \
                   paw.mod projectors.mod pseudopotentials.mod \
                   rundat.mod sparse.mod sparse_embed.mod utils.mod

basis.mod: basis_mod.F90 constants.mod datatypes.mod \
           fft_box.mod geometry.mod rundat.mod \
           simulation_cell.mod timer.mod utils.mod

bibliography.mod: bibliography_mod.F90 comms.mod utils.mod

chebyshev_rep.mod: chebyshev_rep_mod.F90 comms.mod constants.mod rundat.mod \
                   utils.mod

cdft_intermediate_cg.mod: cdft_intermediate_cg_mod.F90 \
                          comms.mod constants.mod electronic.mod \
                          ensemble_dft_type.mod \
                          function_basis.mod hf_exchange.mod \
                          hubbard_build.mod hubbard_init.mod ion.mod \
                          kernel.mod linalg.mod \
                          model_type.mod ngwf_representation.mod \
                          parallel_strategy.mod \
                          rundat.mod services.mod \
                          sparse.mod timer.mod utils.mod

cell_grid.mod: cell_grid_mod.F90 comms.mod constants.mod \
               geometry.mod rundat.mod simulation_cell.mod timer.mod utils.mod

classical_pot.mod: classical_pot_mod.F90  bibliography.mod \
                   cell_grid.mod comms.mod constants.mod \
                   geometry.mod ion.mod rundat.mod \
                   services.mod timer.mod utils.mod

comms.mod: comms_mod.F90 constants.mod

conduction.mod: conduction_mod.F90 bibliography.mod comms.mod \
                conduction_properties.mod constants.mod datatypes.mod \
                electronic_init.mod ensemble_dft_type.mod function_basis.mod \
                hamiltonian.mod hf_exchange.mod hubbard_build.mod kernel.mod \
                model_type.mod projectors.mod ngwf_cg.mod ngwfs.mod \
                ngwf_representation.mod parallel_strategy.mod \
                rundat.mod services.mod sparse_embed.mod utils.mod

conduction_properties.mod: conduction_properties_mod.F90 augmentation.mod \
                           comms.mod constants.mod datatypes.mod \
                           function_basis.mod function_ops.mod \
                           hamiltonian.mod hf_exchange.mod hubbard_build.mod \
                           kernel.mod model_type.mod ngwfs.mod \
                           ngwf_representation.mod \
                           properties.mod projectors.mod rundat.mod \
                           services.mod sparse_embed.mod utils.mod

constants.mod: constants_mod.F90

couplings.mod: couplings_mod.F90 cdft_intermediate_cg.mod comms.mod \
               constants.mod datatypes.mod dense.mod function_basis.mod \
               function_ops.mod hubbard_build.mod model_type.mod \
               ngwf_representation.mod projectors.mod \
               restart.mod rundat.mod sparse.mod \
               sparse_embed.mod timer.mod utils.mod

cutoff_coulomb.mod: cutoff_coulomb_mod.F90 bibliography.mod simulation_cell.mod \
                    cell_grid.mod classical_pot.mod comms.mod \
                    constants.mod density_init.mod fourier.mod \
                    geometry.mod hartree.mod ion.mod model_type.mod \
                    pseudopotentials.mod parallel_strategy.mod \
                    paw.mod rundat.mod timer.mod utils.mod \
                    visual.mod

datatypes.mod: datatypes_mod.F90 constants.mod linalg.mod utils.mod

ddec.mod: ddec_mod.F90 cell_grid.mod comms.mod constants.mod \
          esdf.mod function_basis.mod geometry.mod ion.mod \
          model_type.mod parallel_strategy.mod rundat.mod services.mod \
          timer.mod utils.mod

dense.mod: dense_mod.F90 comms.mod constants.mod datatypes.mod image_comms.mod \
           rundat.mod sparse.mod sparse_embed.mod timer.mod utils.mod

density.mod: density_mod.F90 basis.mod simulation_cell.mod cell_grid.mod \
             comms.mod constants.mod datatypes.mod fft_box.mod fourier.mod \
             fourier_gpu_wrapper_mod.mod function_basis.mod function_ops.mod \
             geometry.mod \
             rundat.mod sparse.mod timer.mod utils.mod

density_init.mod: density_init_mod.F90 basis.mod cell_grid.mod \
                  comms.mod constants.mod fourier.mod geometry.mod integrals.mod \
                  ion.mod model_type.mod paw.mod \
                  rundat.mod services.mod simulation_cell.mod timer.mod \
                  utils.mod visual.mod

dftb.mod: dftb_mod.F90 comms.mod constants.mod dftb_parameters.mod dense.mod \
          ewald.mod function_basis.mod function_ops.mod gaussian_ops.mod \
          geometry.mod ion.mod model_type.mod \
          neighbour_list.mod ngwf_representation.mod parallel_strategy.mod \
          rundat.mod simulation_cell.mod sparse.mod sparse_embed.mod timer.mod \
          utils.mod

dftb_parameters.mod: constants.mod

dftb_type.mod: dense.mod neighbour_list.mod

dma.mod: dma_mod.F90 bibliography.mod cell_grid.mod comms.mod constants.mod \
         datatypes.mod density.mod finite_differences.mod \
         function_basis.mod geometry.mod hash_table.mod \
         ion.mod model_type.mod multigrid_methods.mod multipole_ops.mod \
         ngwf_representation.mod parallel_strategy.mod population.mod \
         rundat.mod sparse.mod sparse_array.mod sparse_embed.mod sw_expansion.mod \
         sw_expansion_type.mod sw_resolution_of_identity.mod timer.mod utils.mod

dmft.mod: dmft_mod.F90 comms.mod constants.mod dense.mod density.mod \
          file_handling.mod function_basis.mod geometry.mod \
          hubbard_build.mod hubbard_init.mod integrals.mod \
          ion.mod kernel.mod linalg.mod model_type.mod \
          ngwf_representation.mod optics.mod parallel_strategy.mod paw.mod \
          projectors.mod pseudopotentials.mod restart.mod \
          rundat.mod simulation_cell.mod sparse.mod sparse_array.mod \
          sparse_embed.mod utils.mod visual.mod

eda_driver_supermol.mod: eda_driver_supermol_mod.F90 comms.mod \
         constants.mod dense.mod density.mod ensemble_dft_type.mod \
         electronic_init.mod fragment_data.mod fragment_matrix_ops.mod \
         fragment_scfmi.mod function_basis.mod hamiltonian.mod hf_exchange.mod \
         hubbard_build.mod kernel.mod \
         model_type.mod ngwf_cg.mod ngwf_representation.mod ngwfs.mod \
         projectors.mod properties.mod restart.mod \
         rundat.mod services.mod sparse.mod sparse_embed.mod utils.mod \
         visual.mod

eda_main.mod: eda_main_mod.F90 comms.mod constants.mod cutoff_coulomb.mod \
              energy_and_force.mod fragment_data.mod function_basis.mod \
              model_type.mod parallel_strategy.mod \
              pbc_corrections.mod rundat.mod utils.mod

eels.mod: eels_mod.F90 simulation_cell.mod comms.mod constants.mod \
          dense.mod function_basis.mod geometry.mod \
          model_type.mod ngwf_representation.mod paw.mod \
          projectors.mod parallel_strategy.mod \
          rundat.mod sparse.mod timer.mod utils.mod

eigenstates.mod: eigenstates_mod.F90 augmentation.mod bibliography.mod \
                 basis.mod cell_grid.mod comms.mod constants.mod \
                 cutoff_coulomb.mod datatypes.mod \
                 dense.mod density.mod eels.mod electronic.mod \
                 electronic_init.mod ensemble_dft.mod ensemble_dft_type.mod \
                 fft_box.mod function_basis.mod function_ops.mod \
                 geometry.mod hartree.mod \
                 hubbard_build.mod hf_exchange.mod integrals.mod ion.mod \
                 is_smeared_ions.mod is_solvation.mod kernel.mod \
                 model_type.mod ngwf_representation.mod ngwfs.mod optics.mod \
                 palser_mano.mod parallel_strategy.mod \
                 population.mod projectors.mod restart.mod rundat.mod \
                 simulation_cell.mod sparse.mod sparse_embed.mod \
                 sparse_array.mod spherical_wave.mod timer.mod utils.mod \
                 visual.mod xc.mod

eld.mod: eld_mod.F90 bibliography.mod cell_grid.mod constants.mod \
         rundat.mod utils.mod services.mod

electronic.mod: electronic_mod.F90 comms.mod constants.mod \
                ensemble_dft.mod ensemble_dft_type.mod ensemble_dftb.mod \
                fragment_data.mod fragment_scfmi.mod \
                function_basis.mod hamiltonian.mod hf_exchange.mod \
                hubbard_build.mod \
                kernel.mod kernel_diis.mod lnv.mod mermin.mod \
                model_type.mod ngwf_representation.mod rundat.mod \
                sparse_embed.mod smearing_operator.mod utils.mod

electronic_init.mod: electronic_init_mod.F90 augmentation.mod \
                classical_pot.mod comms.mod constants.mod \
                cutoff_coulomb.mod datatypes.mod dense.mod density_init.mod \
                electronic_history.mod ensemble_dft.mod ensemble_dft_type.mod \
                ensemble_dftb.mod \
                fragment_matrix_ops.mod function_basis.mod geometry.mod \
                hamiltonian.mod hartree.mod hf_exchange.mod \
                hubbard_build.mod integrals.mod \
                is_smeared_ions.mod kernel.mod model_type.mod \
                ngwf_representation.mod openbc_locps.mod pseudopotentials.mod \
                palser_mano.mod \
                parallel_strategy.mod paw.mod penalty.mod potential.mod \
                projectors.mod restart.mod rundat.mod services.mod \
                sparse.mod sparse_embed.mod utils.mod xc.mod

electronic_history.mod: electronic_history_mod.F90 basis.mod \
                        simulation_cell.mod comms.mod constants.mod \
                        datatypes.mod dense.mod fft_box.mod \
                        function_basis.mod function_ops.mod geometry.mod \
                        kernel.mod ion.mod kernel.mod md_ions.mod \
                        model_type.mod ngwf_representation.mod \
                        parallel_strategy.mod projectors.mod restart.mod \
                        rundat.mod services.mod smearing_operator.mod \
                        sparse.mod sparse_embed.mod \
                        utils.mod

energy_and_force.mod: energy_and_force_mod.F90 augmentation.mod \
                      cdft_intermediate_cg.mod cell_grid.mod \
                      classical_pot.mod comms.mod conduction.mod \
                      constants.mod couplings.mod cutoff_coulomb.mod \
                      datatypes.mod dense.mod density.mod density_init.mod \
                      dftb.mod dma.mod electronic_history.mod \
                      electronic_init.mod eigenstates.mod ensemble_dft.mod \
                      ensemble_dft_type.mod ewald.mod fft_box.mod \
                      forces.mod fourier.mod fragment_data.mod \
                      fragment_matrix_ops.mod eda_driver_supermol.mod fragment_scfmi.mod \
                      function_basis.mod function_ops.mod \
                      hamiltonian.mod hf_exchange.mod hubbard_build.mod \
                      ion.mod is_poisson.mod is_smeared_ions.mod \
                      is_solvation.mod is_solvation_boltzmann.mod kernel.mod \
                      k_points.mod lnv.mod lr_tddft.mod lr_phonons.mod model_type.mod \
                      multigrid_methods.mod ngwf_cg.mod \
                      ngwf_representation.mod ngwfs.mod ngwfs_radial.mod \
                      openbc_locps.mod pseudopotentials.mod \
                      parallel_strategy.mod paw.mod \
                      pbc_corrections.mod penalty.mod polarisable_embedding.mod \
                      potential.mod \
                      ppd_strategy.mod projectors.mod properties.mod \
                      remote.mod restart.mod rundat.mod services.mod \
                      sparse.mod \
                      sparse_embed.mod sph_harm_rotation.mod sw_expansion.mod \
                      sw_resolution_of_identity.mod tddft.mod timer.mod utils.mod \
                      vdwcorrection.mod xc.mod

enthalpy.mod: enthalpy_mod.F90 cell_grid.mod comms.mod \
              constants.mod integrals.mod rundat.mod utils.mod

ensemble_dft.mod: ensemble_dft_mod.F90 comms.mod constants.mod \
                  dense.mod ensemble_dft_type.mod function_basis.mod \
                  hamiltonian.mod hf_exchange.mod hubbard_build.mod \
                  kernel.mod model_type.mod \
                  ngwf_representation.mod restart.mod rundat.mod \
                  smearing_operator.mod sparse.mod \
                  sparse_embed.mod timer.mod utils.mod

ensemble_dft_type.mod: ensemble_dft_type_mod.F90 constants.mod dense.mod \
                       sparse_embed.mod

ensemble_dftb.mod: ensemble_dftb_mod.F90 comms.mod \
                   constants.mod dense.mod dftb.mod ensemble_dft.mod \
                   ensemble_dft_type.mod function_basis.mod \
                   kernel.mod model_type.mod \
                   ngwf_representation.mod restart.mod rundat.mod \
                   smearing_operator.mod \
                   sparse_embed.mod timer.mod utils.mod

esdf_key.mod: esdf_key_mod.F90

esdf.mod: esdf_mod.F90 constants.mod esdf_key.mod utils.mod

etrans.mod: etrans_mod.F90 augmentation.mod comms.mod constants.mod \
            datatypes.mod dense.mod ensemble_dft.mod eigenstates.mod \
            function_basis.mod greenf.mod linalg.mod model_type.mod \
            ngwf_representation.mod parallel_strategy.mod \
            rundat.mod services.mod simulation_cell.mod sparse.mod \
            sparse_embed.mod timer.mod \
            utils.mod

ewald.mod: ewald_mod.F90 simulation_cell.mod comms.mod constants.mod \
           geometry.mod ion.mod rundat.mod timer.mod utils.mod

fft_box.mod: fft_box_mod.F90 comms.mod constants.mod \
             geometry.mod rundat.mod simulation_cell.mod utils.mod

file_handling.mod: file_handling_mod.F90 utils.mod

finite_differences.mod: finite_differences_mod.F90 cell_grid.mod comms.mod \
                        constants.mod geometry.mod timer.mod rundat.mod \
                        utils.mod

forces.mod: forces_mod.F90 augmentation.mod bibliography.mod \
            classical_pot.mod comms.mod constants.mod \
            cutoff_coulomb.mod datatypes.mod density.mod ewald.mod \
            function_basis.mod geometry.mod hamiltonian.mod \
            hf_exchange.mod hubbard_build.mod \
            is_smeared_ions.mod kernel_diis.mod kernel.mod lnv.mod \
            model_type.mod multigrid_methods.mod ngwf_gradient.mod \
            nonsc_forces.mod ngwf_representation.mod openbc_locps.mod \
            parallel_strategy.mod pseudopotentials.mod paw.mod \
            projectors.mod rundat.mod sparse.mod sparse_array.mod sparse_embed.mod \
            timer.mod utils.mod vdwcorrection.mod



forcetest.mod: forcetest_mod.F90 comms.mod constants.mod energy_and_force.mod \
               ion.mod model_type.mod rundat.mod utils.mod

fourier_gpu_wrapper_mod.mod: fourier_gpu_wrapper_mod.F90 constants.mod utils.mod

fourier.mod: fourier_mod.F90 cell_grid.mod comms.mod constants.mod fft_box.mod \
             fourier_gpu_wrapper_mod.mod rundat.mod timer.mod utils.mod

fragment_data.mod: fragment_data_mod.F90 \
                   constants.mod datatypes.mod dense.mod function_basis.mod \
                   model_type.mod ngwf_representation.mod sparse_embed.mod \
                   parallel_strategy.mod rundat.mod utils.mod

fragment_matrix_ops.mod: fragment_matrix_ops_mod.F90 comms.mod constants.mod \
                         datatypes.mod dense.mod fragment_data.mod \
                         function_basis.mod ngwf_representation.mod \
                         rundat.mod sparse.mod sparse_embed.mod timer.mod \
                         utils.mod

fragment_scfmi.mod: fragment_scfmi_mod.F90 \
                    comms.mod constants.mod dense.mod fragment_data.mod \
                    fragment_matrix_ops.mod kernel.mod \
                    ngwf_representation.mod rundat.mod \
                    sparse_embed.mod timer.mod utils.mod

function_basis.mod: function_basis_mod.F90 basis.mod comms.mod constants.mod \
                    datatypes.mod fft_box.mod fourier.mod geometry.mod \
                    hubbard_init.mod parallel_strategy.mod ion.mod rundat.mod \
                    simulation_cell.mod spherical_wave.mod utils.mod

function_ops.mod: function_ops_mod.F90 basis.mod comms.mod constants.mod \
                  datatypes.mod fft_box.mod function_basis.mod parallel_strategy.mod \
                  rundat.mod simulation_cell.mod sparse.mod timer.mod utils.mod

gaunt_coeff.mod: gaunt_coeff_mod.F90 constants.mod utils.mod

gaussian_ops.mod: gaussian_ops_mod.F90 constants.mod geometry.mod utils.mod

geometry.mod: geometry_mod.F90 constants.mod

geometry_deloc_algor.mod: geometry_deloc_algor_mod.F90 constants.mod \
                          geometry_deloc_utils.mod

geometry_deloc_utils.mod: geometry_deloc_utils_mod.F90 comms.mod constants.mod \
                          ion.mod rundat.mod services.mod simulation_cell.mod \
                          utils.mod

geometry_optimiser.mod: geometry_optimiser_mod.F90 simulation_cell.mod \
                        comms.mod constants.mod datatypes.mod \
                        energy_and_force.mod forces.mod \
                        geometry.mod geometry_deloc_algor.mod \
                        geometry_deloc_utils.mod hubbard_init.mod \
                        image_comms.mod ion.mod model_type.mod \
                        linalg.mod neb.mod parallel_strategy.mod \
                        rundat.mod services.mod utils.mod bibliography.mod

greenf.mod: greenf_mod.F90 comms.mod constants.mod linalg.mod neighbour_list.mod \
            parallel_strategy.mod rundat.mod sparse.mod \
            timer.mod utils.mod

hamiltonian.mod: hamiltonian_mod.F90 augmentation.mod comms.mod constants.mod \
                 cutoff_coulomb.mod datatypes.mod dense.mod density.mod dftb.mod \
                 ensemble_dft_type.mod enthalpy.mod \
                 file_handling.mod fragment_matrix_ops.mod fragment_scfmi.mod \
                 function_basis.mod function_ops.mod geometry.mod hartree.mod \
                 hf_exchange.mod hubbard_build.mod integrals.mod \
                 is_poisson.mod is_smeared_ions.mod is_solvation_boltzmann.mod \
                 ke_density.mod model_type.mod ngwf_representation.mod paw.mod \
                 polarisable_embedding.mod \
                 pseudopotentials.mod projectors.mod \
                 rundat.mod scissor.mod smearing_operator.mod sparse.mod \
                 sparse_array.mod sparse_embed.mod timer.mod utils.mod \
                 vdwcorrection.mod visual.mod xc.mod

hartree.mod: hartree_mod.F90 cell_grid.mod constants.mod fourier.mod \
             is_poisson.mod is_smeared_ions.mod is_solvation.mod \
             multigrid_methods.mod pbc_corrections.mod \
             rundat.mod simulation_cell.mod timer.mod utils.mod

hash_table.mod: hash_table_mod.F90 constants.mod rundat.mod timer.mod utils.mod

hf_exchange.mod: hf_exchange_mod.F90 basis.mod bibliography.mod comms.mod \
                 constants.mod datatypes.mod fft_box.mod fourier.mod \
                 function_basis.mod function_ops.mod geometry.mod \
                 hash_table.mod ion.mod model_type.mod neighbour_list.mod \
                 ngwf_representation.mod parallel_strategy.mod ppd_ops.mod \
                 remote.mod rundat.mod simulation_cell.mod sparse.mod \
                 sparse_embed.mod sw_expansion.mod sw_expansion_type.mod \
                 sw_resolution_of_identity.mod timer.mod utils.mod visual.mod \
                 xc.mod

hubbard_build.mod: hubbard_build_mod.F90 augmentation.mod \
                   bibliography.mod simulation_cell.mod comms.mod \
                   constants.mod datatypes.mod dense.mod fft_box.mod \
                   function_basis.mod function_ops.mod geometry.mod \
                   hubbard_init.mod integrals.mod \
                   ion.mod model_type.mod linalg.mod \
                   ngwfs.mod parallel_strategy.mod \
                   projectors.mod restart.mod rundat.mod \
                   sparse.mod timer.mod utils.mod

hubbard_init.mod: hubbard_init_mod.F90 constants.mod utils.mod

image_comms.mod: image_comms_mod.F90 rundat.mod comms.mod utils.mod \
                 constants.mod

integrals.mod: integrals_mod.F90 basis.mod cell_grid.mod \
               comms.mod constants.mod datatypes.mod \
               fft_box.mod fourier.mod function_ops.mod \
               geometry.mod kinetic.mod parallel_strategy.mod \
               potential.mod rundat.mod simulation_cell.mod sparse.mod \
               sparse_embed.mod timer.mod utils.mod function_basis.mod

ion.mod: ion_mod.F90 constants.mod geometry.mod

is_poisson.mod: is_poisson_mod.F90 cell_grid.mod \
                constants.mod simulation_cell.mod \
                utils.mod rundat.mod visual.mod

is_smeared_ions.mod: is_smeared_ions_mod.F90 cell_grid.mod \
                     comms.mod constants.mod fourier.mod geometry.mod \
                     integrals.mod ion.mod \
                     model_type.mod multigrid_methods.mod \
                     parallel_strategy.mod rundat.mod simulation_cell.mod \
                     timer.mod utils.mod

is_solvation.mod: is_solvation_mod.F90 bibliography.mod cell_grid.mod \
                  comms.mod constants.mod finite_differences.mod \
                  fourier.mod geometry.mod integrals.mod ion.mod \
                  multigrid_methods.mod is_poisson.mod is_smeared_ions.mod \
                  rundat.mod services.mod simulation_cell.mod utils.mod \
                  parallel_strategy.mod

is_solvation_boltzmann.mod: is_solvation_boltzmann_mod.F90 bibliography.mod \
                            constants.mod cell_grid.mod comms.mod geometry.mod \
                            finite_differences.mod integrals.mod ion.mod \
                            is_poisson.mod model_type.mod simulation_cell.mod \
                            timer.mod rundat.mod utils.mod visual.mod
# Add explicit dependence on DL_MG modules so that ONETEP code is recompiled
# when DL_MG is modified.
ifeq ($(strip $(MULTIGRID)),dl_mg) # !! only ifdef HAVE_DL_MG
is_solvation_boltzmann.mod: dl_mg_nonlinear_model.mod dl_mg_params.mod dl_mg.mod
endif

k_points.mod: k_points_mod.F90 comms.mod constants.mod \
              geometry.mod ion.mod rundat.mod simulation_cell.mod \
              spglib_f08.mod utils.mod

ke_density.mod: basis.mod cell_grid.mod comms.mod constants.mod datatypes.mod \
                density.mod fourier.mod function_basis.mod \
                fft_box.mod function_ops.mod integrals.mod \
                kinetic.mod rundat.mod simulation_cell.mod sparse.mod \
                timer.mod utils.mod xc.mod

kernel.mod: kernel_mod.F90 comms.mod constants.mod datatypes.mod dense.mod \
            fragment_data.mod fragment_matrix_ops.mod \
            palser_mano.mod rundat.mod sparse.mod \
            sparse_array.mod sparse_embed.mod timer.mod utils.mod

kernel_diis.mod: kernel_diis_mod.F90 comms.mod \
                 constants.mod dense.mod function_basis.mod \
                 hamiltonian.mod hf_exchange.mod \
                 hubbard_build.mod kernel.mod model_type.mod \
                 ngwf_representation.mod rundat.mod \
                 sparse.mod sparse_embed.mod timer.mod utils.mod

kinetic.mod: kinetic_mod.F90 basis.mod constants.mod datatypes.mod \
             fft_box.mod fourier.mod fourier_gpu_wrapper_mod.mod \
             function_basis.mod rundat.mod simulation_cell.mod utils.mod

linalg.mod: linalg_mod.F90 comms.mod constants.mod rundat.mod \
            timer.mod utils.mod

linear_response.mod: linear_response_mod.F90 augmentation.mod \
             cell_grid.mod constants.mod comms.mod density.mod hf_exchange.mod \
             function_basis.mod hartree.mod integrals.mod model_type.mod \
             ngwf_representation.mod paw.mod rundat.mod \
             simulation_cell.mod sparse.mod sparse_embed.mod timer.mod \
             utils.mod xc.mod

lnv.mod: lnv_mod.F90 comms.mod constants.mod \
         fragment_data.mod fragment_scfmi.mod function_basis.mod \
         hamiltonian.mod hf_exchange.mod hubbard_build.mod \
         kernel.mod model_type.mod ngwf_representation.mod restart.mod \
         rundat.mod services.mod sparse_embed.mod timer.mod utils.mod

lr_phonons.mod: lr_phonons_mod.F90 augmentation.mod basis.mod cell_grid.mod \
                comms.mod constants.mod datatypes.mod density.mod dense.mod \
                electronic_init.mod fft_box.mod fourier.mod function_basis.mod \
                function_ops.mod geometry.mod kinetic.mod hamiltonian.mod \
                hartree.mod hf_exchange.mod hubbard_build.mod \
                integrals.mod kernel.mod model_type.mod \
                ngwf_representation.mod ngwfs.mod optics.mod \
                parallel_strategy.mod paw.mod potential.mod \
                projectors.mod pseudopotentials.mod restart.mod \
                rundat.mod services.mod simulation_cell.mod \
                sparse.mod sparse_embed.mod spherical_wave.mod \
                timer.mod utils.mod visual.mod xc.mod

lr_tddft.mod: lr_tddft_mod.F90 augmentation.mod constants.mod comms.mod \
              datatypes.mod dense.mod density.mod electronic_init.mod \
              function_basis.mod hamiltonian.mod hartree.mod \
              hf_exchange.mod hubbard_build.mod \
              is_solvation.mod kernel.mod linalg.mod \
              linear_response.mod lr_tddft_rpa.mod lr_tddft_utils.mod \
              model_type.mod ngwfs.mod ngwf_representation.mod \
              projectors.mod qnto.mod rundat.mod \
              restart.mod services.mod sparse.mod \
              sparse_embed.mod timer.mod utils.mod visual.mod xc.mod

lr_tddft_rpa.mod: lr_tddft_rpa_mod.F90 \
              constants.mod comms.mod density.mod dense.mod   \
              function_basis.mod hartree.mod hf_exchange.mod \
              is_solvation.mod linalg.mod \
              linear_response.mod lr_tddft_utils.mod model_type.mod \
              ngwf_representation.mod projectors.mod rundat.mod \
              restart.mod sparse.mod sparse_embed.mod \
              timer.mod visual.mod utils.mod

lr_tddft_utils.mod: lr_tddft_utils_mod.F90 cell_grid.mod \
                    constants.mod comms.mod datatypes.mod \
                    dense.mod density.mod function_basis.mod hf_exchange.mod \
                    integrals.mod linalg.mod linear_response.mod model_type.mod \
                    ngwf_representation.mod optics.mod projectors.mod \
                    rundat.mod sparse.mod sparse_embed.mod timer.mod utils.mod \
                    visual.mod

md.mod: md_mod.F90 simulation_cell.mod comms.mod constants.mod \
        electronic_history.mod energy_and_force.mod esdf.mod \
        ion.mod md_ions.mod md_thermostat.mod model_type.mod \
        rundat.mod utils.mod services.mod

md_ions.mod: md_ions_mod.F90 comms.mod constants.mod \
             ion.mod linalg.mod rundat.mod \
             utils.mod services.mod

md_thermostat.mod: md_thermostat_mod.F90 comms.mod \
                   constants.mod md_ions.mod rundat.mod  \
                   services.mod utils.mod

mermin.mod: mermin_mod.F90 comms.mod constants.mod \
            ensemble_dft.mod function_basis.mod \
            hamiltonian.mod hf_exchange.mod hubbard_build.mod \
            kernel.mod model_type.mod ngwf_representation.mod restart.mod \
            rundat.mod services.mod smearing_operator.mod sparse_embed.mod timer.mod utils.mod

model_type.mod: model_type_mod.F90 cell_grid.mod comms.mod constants.mod \
                dftb_parameters.mod dftb_type.mod fft_box.mod ion.mod k_points.mod \
                parallel_strategy.mod paw.mod pseudopotentials.mod \
                simulation_cell.mod utils.mod

multigrid_methods.mod: multigrid_methods_mod.F90 bibliography.mod \
                       cell_grid.mod comms.mod constants.mod \
                       finite_differences.mod geometry.mod integrals.mod \
                       ion.mod is_poisson.mod is_solvation_boltzmann.mod \
                       rundat.mod simulation_cell.mod timer.mod utils.mod \
                       visual.mod

# Add explicit dependence on DL_MG modules so that ONETEP code is recompiled
# when DL_MG is modified.
ifeq ($(strip $(MULTIGRID)),dl_mg) # !! only ifdef HAVE_DL_MG
multigrid_methods.mod: dl_mg.mod
endif

neb.mod: neb_mod.F90 constants.mod comms.mod energy_and_force.mod esdf.mod \
         geometry.mod image_comms.mod model_type.mod rundat.mod \
         rundat_blocks.mod services.mod simulation_cell.mod utils.mod \
         forces.mod

multipole_ops.mod: multipole_ops_mod.F90 comms.mod constants.mod geometry.mod \
                   rundat.mod services.mod spherical_wave.mod timer.mod \
                   utils.mod

neighbour_list.mod: neighbour_list_mod.F90 comms.mod constants.mod \
                    function_basis.mod parallel_strategy.mod rundat.mod \
                    sparse.mod timer.mod utils.mod

ngwf_cg.mod: ngwf_cg_mod.F90 augmentation.mod basis.mod \
             bibliography.mod cdft_intermediate_cg.mod comms.mod constants.mod \
             datatypes.mod dense.mod electronic.mod \
             electronic_history.mod ensemble_dft_type.mod \
             forces.mod fragment_data.mod fragment_scfmi.mod \
             function_basis.mod function_ops.mod geometry.mod \
             hamiltonian.mod hf_exchange.mod hubbard_build.mod kernel.mod \
             mermin.mod model_type.mod ngwf_gradient.mod \
             ngwf_representation.mod nho.mod nonsc_forces.mod parallel_strategy.mod \
             polarisable_embedding.mod projectors.mod \
             restart.mod rundat.mod services.mod smearing_operator.mod sparse.mod \
             sparse_embed.mod timer.mod utils.mod visual.mod

ngwf_data.mod: ngwf_data_mod.F90 constants.mod

ngwf_gradient.mod: ngwf_gradient_mod.F90 augmentation.mod basis.mod \
                   simulation_cell.mod cell_grid.mod comms.mod constants.mod \
                   datatypes.mod ensemble_dft.mod fft_box.mod\
                   fourier.mod fourier_gpu_wrapper_mod.mod \
                   function_basis.mod fragment_data.mod fragment_scfmi.mod \
                   function_ops.mod geometry.mod \
                   hf_exchange.mod \
                   hubbard_build.mod kernel.mod kernel_diis.mod \
                   kinetic.mod mermin.mod model_type.mod ngwf_representation.mod \
                   polarisable_embedding.mod pseudopotentials.mod paw.mod \
                   potential.mod projectors.mod restart.mod rundat.mod \
                   scissor.mod sparse.mod sparse_array.mod sparse_embed.mod \
                   timer.mod utils.mod xc.mod

ngwf_representation.mod: ngwf_representation_mod.F90 comms.mod \
                         constants.mod datatypes.mod hash_table.mod \
                         model_type.mod remote.mod rundat.mod \
                         sparse_initialise.mod sparse_embed.mod \
                         sw_expansion_type.mod utils.mod

ngwfs.mod: ngwfs_mod.F90 basis.mod simulation_cell.mod comms.mod \
           constants.mod datatypes.mod dftb_parameters.mod electronic_history.mod \
           fragment_data.mod function_basis.mod geometry.mod hf_exchange.mod \
           fft_box.mod ion.mod model_type.mod ngwf_data.mod \
           ngwf_representation.mod \
           parallel_strategy.mod polarisable_embedding.mod \
           projectors.mod restart.mod rundat.mod services.mod \
           spherical_wave.mod timer.mod utils.mod

ngwfs_radial.mod: ngwfs_radial_mod.F90 atom.mod comms.mod \
                  constants.mod density_init.mod dftb_parameters.mod\
                  ion.mod model_type.mod ngwf_data.mod \
                  rundat.mod services.mod utils.mod

nho.mod: nho_mod.F90 simulation_cell.mod comms.mod constants.mod \
         datatypes.mod function_basis.mod geometry.mod ion.mod linalg.mod \
         parallel_strategy.mod rundat.mod sparse.mod timer.mod \
         utils.mod

nlxc.mod: nlxc_mod.F90 cell_grid.mod simulation_cell.mod comms.mod \
          constants.mod fourier.mod geometry.mod rundat.mod \
          timer.mod utils.mod vdw_df_kernel.mod xc_funcs.mod

nonsc_forces.mod: nonsc_forces_mod.F90 bibliography.mod comms.mod \
                  constants.mod datatypes.mod function_basis.mod \
                  fft_box.mod integrals.mod parallel_strategy.mod rundat.mod \
                  simulation_cell.mod \
                  sparse.mod timer.mod utils.mod

npa.mod: npa_mod.F90 bibliography.mod comms.mod constants.mod \
         dense.mod function_basis.mod ion.mod kernel.mod linalg.mod \
         ngwf_data.mod ngwf_representation.mod parallel_strategy.mod \
         rundat.mod sparse.mod sparse_embed.mod timer.mod utils.mod

optics.mod: optics_mod.F90 augmentation.mod basis.mod bibliography.mod \
            comms.mod constants.mod datatypes.mod dense.mod \
            function_basis.mod geometry.mod integrals.mod \
            ion.mod model_type.mod ngwf_representation.mod \
            pseudopotentials.mod \
            parallel_strategy.mod projectors.mod rundat.mod \
            simulation_cell.mod sparse.mod sparse_embed.mod timer.mod utils.mod

openbc_locps.mod: openbc_locps_mod.F90 cell_grid.mod comms.mod constants.mod \
                  geometry.mod ion.mod parallel_strategy.mod \
                  paw.mod pseudopotentials.mod simulation_cell.mod rundat.mod \
                  services.mod timer.mod utils.mod

onetep.o: onetep.F90 alloc_checker.mod anharmonic.mod bibliography.mod \
          comms.mod constants.mod \
          couplings.mod  eda_main.mod energy_and_force.mod esdf.mod forcetest.mod \
          geometry_optimiser.mod hubbard_init.mod image_comms.mod k_points.mod \
          md.mod model_type.mod parallel_strategy.mod phonon.mod rundat.mod \
          rundat_blocks.mod services.mod \
          spherical_wave.mod timer.mod tssearch.mod unit_test.mod \
          utils.mod vdwcorrection.mod

palser_mano.mod: palser_mano_mod.F90 comms.mod constants.mod \
                 rundat.mod services.mod sparse_embed.mod timer.mod utils.mod

parallel_strategy.mod: parallel_strategy_mod.F90 comms.mod \
                       constants.mod geometry.mod hubbard_init.mod \
                       ion.mod k_points.mod rundat.mod simulation_cell.mod utils.mod

paw.mod: paw_mod.F90 basis.mod bibliography.mod classical_pot.mod \
         cell_grid.mod constants.mod comms.mod ion.mod fft_box.mod fourier.mod \
         gaunt_coeff.mod geometry.mod image_comms.mod \
         linalg.mod paw_shape.mod paw_xc.mod parallel_strategy.mod \
         projectors.mod rundat.mod services.mod simulation_cell.mod \
         sparse.mod spherical_wave.mod timer.mod utils.mod xc.mod

paw_shape.mod: paw_shape_mod.F90 basis.mod cell_grid.mod \
               constants.mod geometry.mod services.mod \
               simulation_cell.mod spherical_wave.mod utils.mod

paw_xc.mod: paw_xc_mod.F90 constants.mod gaunt_coeff.mod \
            services.mod xc.mod

pbc_corrections.mod: pbc_corrections_mod.F90 bibliography.mod  \
                     cell_grid.mod comms.mod constants.mod \
                     fourier.mod geometry.mod ion.mod parallel_strategy.mod \
                     rundat.mod simulation_cell.mod timer.mod utils.mod

penalty.mod: penalty_mod.F90 comms.mod \
             constants.mod fragment_scfmi.mod function_basis.mod \
             hamiltonian.mod hf_exchange.mod \
             hubbard_build.mod kernel.mod model_type.mod \
             ngwf_representation.mod restart.mod rundat.mod services.mod \
             sparse.mod sparse_embed.mod timer.mod utils.mod

phonon.mod: phonon_mod.F90 bibliography.mod dense.mod constants.mod \
            energy_and_force.mod comms.mod model_type.mod \
            utils.mod rundat.mod services.mod

polarisable_embedding.mod: polarisable_embedding_mod.F90 basis.mod \
                           bibliography.mod cell_grid.mod comms.mod \
                           constants.mod datatypes.mod dma.mod fft_box.mod \
                           fourier.mod function_basis.mod function_ops.mod \
                           geometry.mod hash_table.mod integrals.mod ion.mod \
                           linalg.mod kernel.mod \
                           model_type.mod multipole_ops.mod neighbour_list.mod \
                           ngwf_representation.mod parallel_strategy.mod \
                           ppd_ops.mod remote.mod restart.mod rundat.mod \
                           simulation_cell.mod sparse.mod \
                           sparse_array.mod sparse_embed.mod sw_expansion_type.mod \
                           sw_resolution_of_identity.mod timer.mod utils.mod \
                           visual.mod

population.mod: population_mod.F90 comms.mod \
                constants.mod ddec.mod dense.mod density.mod \
                function_basis.mod geometry.mod model_type.mod \
                ngwf_representation.mod parallel_strategy.mod \
                rundat.mod sparse.mod sparse_embed.mod timer.mod utils.mod

potential.mod: potential_mod.F90 basis.mod simulation_cell.mod cell_grid.mod \
               comms.mod constants.mod datatypes.mod fft_box.mod\
               fourier.mod fourier_gpu_wrapper_mod.mod \
               function_basis.mod geometry.mod \
               geometry.mod \
               parallel_strategy.mod rundat.mod timer.mod utils.mod

ppd_ops.mod: ppd_ops_mod.F90 basis.mod constants.mod datatypes.mod \
             remote.mod rundat.mod utils.mod

ppd_strategy.mod: ppd_strategy_mod.F90 comms.mod \
                  constants.mod fft_box.mod geometry.mod \
                  ion.mod rundat.mod simulation_cell.mod utils.mod

projectors.mod: projectors_mod.F90 basis.mod simulation_cell.mod comms.mod \
                constants.mod datatypes.mod fft_box.mod fourier.mod \
                function_basis.mod function_ops.mod geometry.mod \
                parallel_strategy.mod rundat.mod services.mod sparse.mod \
                sparse_embed.mod spherical_wave.mod timer.mod utils.mod

properties.mod: properties_mod.F90 augmentation.mod \
                bandstructure.mod basis.mod cell_grid.mod comms.mod \
                constants.mod cutoff_coulomb.mod eld.mod datatypes.mod \
                dense.mod density.mod dma.mod dmft.mod \
                eigenstates.mod etrans.mod function_basis.mod \
                function_ops.mod geometry.mod hamiltonian.mod hartree.mod \
                hf_exchange.mod \
                hubbard_build.mod integrals.mod is_smeared_ions.mod \
                ke_density.mod kernel.mod linalg.mod md_ions.mod model_type.mod \
                multipole_ops.mod ngwf_representation.mod npa.mod \
                optics.mod parallel_strategy.mod paw.mod population.mod \
                projectors.mod restart.mod rundat.mod sparse.mod \
                sparse_array.mod sparse_embed.mod spherical_wave.mod \
                timer.mod utils.mod visual.mod

pseudopotentials.mod: pseudopotentials_mod.F90 cell_grid.mod \
                      classical_pot.mod comms.mod constants.mod \
                      datatypes.mod fft_box.mod fourier.mod function_basis.mod \
                      gaunt_coeff.mod geometry.mod ion.mod \
                      image_comms.mod parallel_strategy.mod \
                      pbc_corrections.mod projectors.mod rundat.mod \
                      services.mod simulation_cell.mod sparse.mod \
                      sparse_embed.mod spherical_wave.mod timer.mod \
                      utils.mod xc.mod

qnto.mod: qnto_mod.F90 augmentation.mod comms.mod constants.mod datatypes.mod \
          dense.mod eigenstates.mod function_basis.mod function_ops.mod \
          ion.mod model_type.mod ngwf_representation.mod ngwfs.mod npa.mod \
          restart.mod rundat.mod sparse.mod sparse_embed.mod \
          timer.mod utils.mod

remote.mod: remote_mod.F90 basis.mod comms.mod constants.mod datatypes.mod \
            fft_box.mod function_basis.mod geometry.mod hash_table.mod ion.mod \
            parallel_strategy.mod rundat.mod simulation_cell.mod sparse.mod \
            timer.mod utils.mod

restart.mod: restart_mod.F90 basis.mod comms.mod constants.mod \
             datatypes.mod fft_box.mod function_basis.mod geometry.mod \
             ion.mod model_type.mod parallel_strategy.mod rundat.mod simulation_cell.mod \
             spherical_wave.mod sparse_embed.mod \
             timer.mod utils.mod

rundat.mod: rundat_mod.F90 comms.mod constants.mod esdf.mod utils.mod

rundat_blocks.mod: rundat_blocks_mod.F90 classical_pot.mod \
                   comms.mod constants.mod couplings.mod dftb.mod esdf.mod \
                   fragment_data.mod \
                   geometry.mod hubbard_init.mod ion.mod is_smeared_ions.mod \
                   is_solvation_boltzmann.mod md_thermostat.mod \
                   model_type.mod parallel_strategy.mod ppd_strategy.mod \
                   polarisable_embedding.mod \
                   rundat.mod services.mod simulation_cell.mod \
                   sw_resolution_of_identity.mod utils.mod etrans.mod \
                   vdwcorrection.mod

scissor.mod: comms.mod constants.mod function_basis.mod model_type.mod \
             parallel_strategy.mod rundat.mod sparse.mod utils.mod

services.mod: services_mod.F90 cell_grid.mod comms.mod constants.mod \
                   datatypes.mod fourier.mod ion.mod rundat.mod \
                   simulation_cell.mod utils.mod

simulation_cell.mod: simulation_cell_mod.F90 comms.mod constants.mod \
                     geometry.mod ion.mod rundat.mod utils.mod

smearing_operator.mod: smearing_operator_mod.F90 comms.mod constants.mod \
                       dense.mod sparse.mod timer.mod utils.mod linalg.mod \
                       rundat.mod sparse_embed.mod

sparse.mod: sparse_mod.F90 comms.mod constants.mod datatypes.mod \
            linalg.mod parallel_strategy.mod rundat.mod timer.mod \
            utils.mod

sparse_algor.mod: sparse_algor_mod.F90 constants.mod \
                  parallel_strategy.mod sparse.mod

sparse_array.mod: sparse_array_mod.F90 sparse.mod constants.mod utils.mod

sparse_embed.mod: sparse_embed_mod.F90 comms.mod constants.mod datatypes.mod \
                  parallel_strategy.mod \
                  rundat.mod sparse.mod sparse_array.mod timer.mod utils.mod

sparse_initialise.mod: sparse_initialise_mod.F90 constants.mod comms.mod \
                       ion.mod model_type.mod parallel_strategy.mod rundat.mod \
                       simulation_cell.mod sparse.mod sparse_embed.mod \
                       timer.mod utils.mod

spglib_f08.mod: spglib_f08_mod.F90

sph_harm_rotation.mod: sph_harm_rotation_mod.F90 \
                       comms.mod constants.mod dense.mod geometry.mod \
                       rundat.mod services.mod utils.mod

spherical_wave.mod: spherical_wave_mod.F90 fft_box.mod \
                    constants.mod fourier.mod \
                    geometry.mod rundat.mod timer.mod utils.mod

sw_expansion.mod: sw_expansion_mod.F90 comms.mod constants.mod \
                  datatypes.mod fft_box.mod fourier.mod function_basis.mod \
                  geometry.mod hash_table.mod \
                  ion.mod neighbour_list.mod ngwf_representation.mod \
                  linalg.mod parallel_strategy.mod ppd_ops.mod remote.mod \
                  rundat.mod simulation_cell.mod sparse.mod \
                  sw_expansion_type.mod sw_resolution_of_identity.mod \
                  timer.mod utils.mod

sw_expansion_type.mod: sw_expansion_type_mod.F90 comms.mod constants.mod \
                       hash_table.mod rundat.mod sw_resolution_of_identity.mod \
                       utils.mod

sw_resolution_of_identity.mod: sw_resolution_of_identity_mod.F90 basis.mod \
                 bibliography.mod chebyshev_rep.mod constants.mod comms.mod \
                 geometry.mod fft_box.mod hash_table.mod ion.mod neighbour_list.mod \
                 linalg.mod parallel_strategy.mod \
                 rundat.mod simulation_cell.mod sparse.mod \
                 sph_harm_rotation.mod spherical_wave.mod \
                 timer.mod utils.mod

tddft.mod: tddft_mod.F90 cell_grid.mod comms.mod \
           constants.mod datatypes.mod function_basis.mod \
           geometry.mod hamiltonian.mod integrals.mod \
           model_type.mod linalg.mod ngwf_representation.mod properties.mod \
           rundat.mod sparse.mod sparse_embed.mod timer.mod \
           utils.mod vdwcorrection.mod xc.mod

timer.mod: timer_mod.F90 comms.mod constants.mod rundat.mod utils.mod

tssearch_algor.mod: tssearch_algor_mod.F90 tssearch_utils.mod

tssearch.mod: tssearch_mod.F90 comms.mod \
              constants.mod esdf.mod geometry_optimiser.mod ion.mod \
              neb.mod model_type.mod rundat.mod \
              rundat_blocks.mod services.mod simulation_cell.mod \
              tssearch_algor.mod tssearch_utils.mod

tssearch_utils.mod: tssearch_utils_mod.F90 simulation_cell.mod comms.mod \
                    constants.mod energy_and_force.mod forces.mod \
                    ion.mod model_type.mod utils.mod

unit_test.mod: unit_test_mod.F90 cell_grid.mod comms.mod constants.mod \
               fourier.mod model_type.mod rundat.mod services.mod \
               utils.mod

utils.mod: utils_mod.F90 comms.mod constants.mod geometry.mod

vdw_df_kernel.mod: vdw_df_kernel_mod.F90  comms.mod constants.mod \
                   rundat.mod services.mod utils.mod

vdwcorrection.mod: vdwcorrection_mod.F90 bibliography.mod comms.mod \
                   constants.mod geometry.mod ion.mod parallel_strategy.mod \
                   rundat.mod simulation_cell.mod timer.mod utils.mod

visual.mod: visual_mod.F90 basis.mod simulation_cell.mod cell_grid.mod \
            comms.mod constants.mod datatypes.mod density.mod fft_box.mod \
            finite_differences.mod fourier.mod function_basis.mod geometry.mod \
            ion.mod model_type.mod ngwf_representation.mod \
            parallel_strategy.mod rundat.mod services.mod sparse.mod timer.mod \
            utils.mod

xc_b97mv.mod: xc_b97mv_mod.F90 constants.mod

xc_funcs.mod: xc_funcs_mod.F90 constants.mod rundat.mod utils.mod xc_b97mv.mod

xc.mod: xc_mod.F90 bibliography.mod cell_grid.mod \
        comms.mod constants.mod fourier.mod nlxc.mod rundat.mod \
        services.mod simulation_cell.mod timer.mod utils.mod xc_funcs.mod

