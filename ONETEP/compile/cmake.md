# build with cmake

git clone git@bitbucket.org:onetep/onetep.git
## be sure you have the needed libs in the environment

### build foss (openmpi/fftw3/openblas/sclapack)

```bash
  cmake -S ~/onetep -Bbuild-onetep -DWITH_SCALAPACK=On -DWITH_BLAS=On
  cmake --build build-onetep -j10
  #or if on linux with make
  cd build-ontep
  make -j10
```

### intel compilers and mkl with scalapack/ffw3

```bash

  cmake -S~/onetep -Bbuild-onetep-intel -DWITH_SCALAPACK=On -DWITH_MKL=On -DFFT="MKL_FFTW3"
  cmake --build build-onetep-intel -j 10
  #or if on linux with make

  cd build-onetep-intel
  make -j10
```

changing Fortran compiler flags


```bash
  FFLAGS="-O2 -g" cmake -S ~/onetep -Bbuild-onetep -DWITH_SCALAPACK=On -DWITH_BLAS=On
```
