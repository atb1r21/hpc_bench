#!/bin/bash --login

#echo $PATH

set -u
#set -x

# defaults
nx=17
ny=17
nz=17
npx=1
npy=1
npz=1
periodicity=(F F F)
agg_level=0
niter=40
nthreads=1
slice_index=-1
slice_section=-1
test_name=poisson   #laplace is deprecated
model_name=sin_eps  #laplace #sin_eps #parabolic_rho
tol_res_rel=1.e-2
tol_res_abs=1.e-2
tol_pot_rel=1.e-2
tol_pot_abs=1.e-2
tol_nwt_rel=1.e-5
tol_vcy_rel=1.e-5
tol_list=( 0 0 0 0 0 0 )
defco_order=2
use_damping=F
use_steric=F
use_linearised=F
exename=../dl_mg_test.exe
nompi=0
mpiexec=mpirun
wdir=.

while getopts "g:p:t:x:m:ZT:P:D:SLw:h" opt; do
  case "$opt" in
  g) ga=($(echo $OPTARG | sed -e 's/-/ /g'))
     nx=${ga[0]}
     ny=${ga[1]}
     nz=${ga[2]}
     ;;
  p) ga=($(echo $OPTARG | sed -e 's/-/ /g'))
     npx=${ga[0]}
     npy=${ga[1]}
     npz=${ga[2]}
     ;;
  t) nthreads=$OPTARG
     ;;
  x) exename=$OPTARG
     ;;
  m) model_name=$OPTARG
     ;;
  Z) nompi=1
     ;;
  T) ga=($(echo $OPTARG | sed -e 's/,/ /g') )
     if [ ! ${ga[0]} = 0 ] ; then tol_res_rel=${ga[0]} ; fi
     if [ ! ${ga[1]} = 0 ] ; then tol_res_abs=${ga[1]} ; fi
     if [ ! ${ga[2]} = 0 ] ; then tol_pot_rel=${ga[2]} ; fi
     if [ ! ${ga[3]} = 0 ] ; then tol_pot_abs=${ga[3]} ; fi
     if [ ! ${ga[4]} = 0 ] ; then tol_nwt_rel=${ga[4]} ; fi
     if [ ! ${ga[5]} = 0 ] ; then tol_vcy_rel=${ga[5]} ; fi
     ;;
  P) periodicity=($(echo $OPTARG | sed -e 's/-/ /g'))
     ;;
  D) defco_order=$OPTARG
     ;;
  S) use_steric=T
     ;;
  L) use_linearised=T
     ;;
  w) wdir=$OPTARG
     ;;
  h) echo "Usage: $(basename $0) -ZLS -g <nx>-<ny>-<nz> -p <px>-<py>-<pz> -t <nthreads> -x <executable> -m <model name> -P {T|F}-{T|F}-{T|F} -D <defco order> -w <work dir> -T {<tol_res_rel>|0},{<tol_res_abs>|0},{<tol_pot_rel>|0},{<tol_pot_abs>|0},{<tol_nwt_rel>|0},{<tol_vcy_rel>|0}"
     echo "-Z: no mpi"
     echo "-S: use steric weight"
     echo "-L: use linearised PBE"
     exit 0
     ;;
  \?) echo "unknown flag"; exit 1
     ;;
   :) echo "-$opt needs an argument!"; exit 1
     ;;
  esac
done

echo "grid sizes: $nx $ny $nz"
echo "periodicity ${periodicity[@]}"
echo "mpi topo  : $npx $npy $npz"

nproc=$((npx * npy * npz))

echo "nproc $nproc"

cd $wdir


cat > input << EndOfInput
# global grid sizes
$nx $ny $nz

# MPI grid
$npx $npy $npz

# periodicity
${periodicity[@]}

# aggregation level
$agg_level

# number of iterations (defect  correction, multigrid)
 $niter $niter

# write a plane in a given section
-1 -1

# test name (laplace poisson ...), model
poisson $model_name

# tol (defect correction: residual {abs,rel}, pot (abs,rel}; newton rel; mg_rel
$tol_res_rel $tol_res_abs $tol_pot_rel $tol_pot_abs $tol_nwt_rel $tol_vcy_rel

# finite diference order for defect correction; use_damping
$defco_order $use_damping

# PBE temp n fixed_nions
300.0 2 F

# PBE concetration in Mol/l
0.1 0.1

# PBE charges (in electron charge)
1.0 -1.0

# PBE lambda,  take second value only first is T (for F use default)
F 0.0

# PBE  linearized, use steric, use full approx scheme
$use_linearised $use_steric F

EndOfInput

# run the code
echo "runnig with $exename with $nproc processors"
if [ -n "$nthreads" ]
then
  export OMP_NUM_THREADS="$nthreads"
  echo " using OpenMP with $nthreads threads "
fi

echo "using $exename with $nproc MPI tasks"
echo "using mpiexec $mpiexec"

if [ "$nompi" = 0 ]
then
  $mpiexec -n $nproc $exename
else
  $exename
fi
