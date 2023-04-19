cmake_minimum_required(VERSION 3.0.2)
project(onetep
  LANGUAGES Fortran
  VERSION 6.1.11.4)

set(AUTHOR "Jacek Dziedzic;Peter D. Haynes;Nicholas D. M. Hine;Arash A. Mostofi;Mike C. Payne;Chris-Kriton Skylaris")
set(AUTHOR_DETAILS "onetep.ac.uk")
set(DESCRIPTION "Onetep linear scaling electronic structure code")

option(WITH_MPI "build a MPI version" ON)
option(WITH_DLMG "build with DL_MG version" ON)
option(WITH_OPENMP "build with OpenMP version" OFF)
option(WITH_SCALAPACK "build with ScaLAPACK version" OFF)
option(WITH_BLAS "build with blas/openblas" OFF)
option(WITH_MKL "build using mkl" OFF)
option(WITH_MPI "build a MPI version" ON)
option(ENABLE_BLAS95 "enable blas95 interface" ON)
option(ENABLE_SCALAPACK "enable scalapack interface" ON)

set(FFT "FFTW3" CACHE STRING "what fft flavour to use")

set(EXTRA_LIBS "" CACHE STRING "extra libraries to be passed to linker")
set(EXTRA_LIB_PATHS "" CACHE STRING "extra library paths to be passed to linker")
set(EXTRA_INC_PATHS "" CACHE STRING "extra include paths to be passed to compiler")
set(EXTRA_DEFS "" CACHE STRING "extra defines to be passed to compiler")

set(EXE "onetep.x" CACHE STRING "binary's name")

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message(STATUS "Build release version")
endif()

set(src_dir ${CMAKE_SOURCE_DIR}/src)
set(target_name ${EXE})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake" "${CMAKE_SOURCE_DIR}/cmake/Modules")

include(GNUInstallDirs)
include(flags)

##set the output folder for libs and bins
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

if(WITH_MPI)
  ### if we want MPI check for it and set the internal kitchen
  find_package(MPI REQUIRED)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  message(STATUS "MPI_Fortran_INCLUDE_PATH ${MPI_Fortran_INCLUDE_PATH}")
  set(LIBS ${LIBS} ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
  string(REGEX REPLACE "^ " "" LIBS "${LIBS}")
  message(STATUS "Build with MPI support!")
else()
  message(STATUS "Build without MPI support!")
endif()

find_package(PkgConfig REQUIRED)
if (WITH_BLAS)
  if (WITH_SCALAPACK)
    pkg_check_modules(SCALAPACK scalapack>=2.0.2)
    if (NOT SCALAPACK_FOUND)
      pkg_check_modules(SCALAPACK scalapack-openmpi>=2.0.2)
    endif()
    if (SCALAPACK_FOUND)
      link_directories(${SCALAPACK_LIBRARY_DIRS})
      set(LIBS ${LIBS} ${SCALAPACK_LIBRARIES})
    else()
      message(WARNING "No ScaLAPACK found, please provide your own using EXTRA_LIBS!!")
    endif()
  endif()

  set(BLA_PREFER_PKGCONFIG On)
  find_package(BLAS)
  if (BLAS_FOUND)
    set(LIBS ${LIBS} ${BLAS_LIBRARIES})
  else()
    message(WARNING "requested fft ${FFT} not found used EXTRA_LIBS to provide one")
  endif()
endif()

if (WITH_OPENMP)
  find_package(OpenMP    REQUIRED)
endif()

if(WITH_MKL)
  if (WITH_SCALAPACK)
    set(ENABLE_SCALAPACK On)
  endif()
  set(MKL_INTERFACE_FULL "intel_lp64" )
  find_package(MKL CONFIG)
  if (MKL_FOUND)
    set(LIBS ${LIBS} ${MKL_IMPORTED_TARGETS})
    message(STATUS "MKL includes ${MKL_INCLUDE}")
  endif()
endif()

set(fft_not_found TRUE )
if(${FFT} STREQUAL "FFTW3")
  pkg_search_module(FFTW3 fftw3>=3.1.0 IMPORTED_TARGET)
  if (FFTW3_FOUND)
    set(fft_not_found FALSE )
    set(LIBS ${LIBS} ${FFTW3_LIBRARIES})
  endif()
elseif(${FFT} STREQUAL "MKL_FFTW3")
  if (MKL_FOUND)
    set(fft_not_found FALSE )
  endif()
endif()

if (fft_not_found)
  message(WARNING "requested fft ${FFT} not found used EXTRA_LIBS to provide one")
endif()
set(LIBS ${LIBS} ${EXTRA_LIBS})

if (WITH_DLMG)
  add_subdirectory(dl_mg)
endif()

message(STATUS "Using these libs to link: ${LIBS}")
message(STATUS "Using these extra libs to link: ${EXTRA_LIBS}")

add_subdirectory(src)

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENCE")
set(CPACK_PACKAGE_VERSION_MAJOR ${onetep_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${onetep_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${onetepe_VERSION_PATCH})
set(CPACK_PACKAGE_VERSION_TWEAK ${onetepe_VERSION_TWEAK})
set(CPACK_INSTALL_CMAKE_CONFIGURATIONS Release)

set(CPACK_SOURCE_GENERATOR "TGZ")
include(CPack)

