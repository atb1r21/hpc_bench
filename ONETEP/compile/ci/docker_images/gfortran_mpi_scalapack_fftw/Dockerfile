FROM debian:stretch-slim

LABEL maintainer="Jose M Escartin Esteban <jme52@cam.ac.uk>"
LABEL version="0.3"
LABEL os="debian-9"
LABEL compiler="gfortran"
LABEL compiler_version="6.3"
LABEL MPI="OpenMPI-2.0.2"
LABEL FFT="FFTW-3.3.5"

# We need:
#   * gawk because mawk (POSIX awk) breaks checkdeps.
#   * make.
#   * libfile-copy-recursive-perl for compile_mod.pl.
#   * gfortran (brings gcc, including cpp).
#   * OpenMPI with wrappers.
#   * libraries: blas, lapack, scalapack, fftw3.

RUN apt-get -q update && \
    apt-get -qy install --no-install-recommends \
       gawk                                     \
       make                                     \
       libfile-copy-recursive-perl              \
       gfortran=4:6.3.0-4                       \
       gfortran-6=6.3.0-18+deb9u1               \
       libopenmpi-dev=2.0.2-2                   \
       libblas3=3.7.0-2                         \
       liblapack3=3.7.0-2                       \
       libscalapack-openmpi1=1.8.0-13           \
       libfftw3-3=3.3.5-3                       \
    && rm -rf /var/lib/apt/lists/*

# Create dynamic links
# (this could also be accomplished by installing the corresponding -dev packages,
# but at the cost of also installing the static libraries, which would increase
# the size of the docker image by hundreds of MB).

RUN ln -sf libfftw3_omp.so.3.5.5              /usr/lib/x86_64-linux-gnu/libfftw3_omp.so && \
    ln -sf libfftw3.so.3.5.5                  /usr/lib/x86_64-linux-gnu/libfftw3.so     && \
    ln -sf libscalapack-openmpi.so.1          /usr/lib/libscalapack-openmpi.so          && \
    ln -sf libblacs-openmpi.so.1              /usr/lib/libblacs-openmpi.so              && \
    ln -sf libblacsCinit-openmpi.so.1         /usr/lib/libblacsCinit-openmpi.so         && \
    ln -sf /usr/lib/lapack/liblapack.so.3.7.0 /usr/lib/liblapack.so                     && \
    ln -sf /usr/lib/libblas/libblas.so.3.7.0  /usr/lib/libblas.so
