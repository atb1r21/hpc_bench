#!/bin/bash --login

# set of tests that cover various grid size, MPI topologies
# and number of open threads accross a set of simple problems,
#see below
#
# Lucian Anton
# July 2014
#

set -u
declare -a parallel_deviations tol_list failed_run_list failed_list

# defaults
test_list="linear_phi sin_eps pbc wirex wirez pbe_sphere_hardcore_steric pbe_sphere_hardcore_steric+fullpbc"
#test_list="linear_phi sin_eps pbc wirex wirez "
MPIEXEC=nompi
parmode=mpi+omp
use_cg=T

# grid sizes
nx=17
ny=25
nz=33

MPI_period="F F F"
agg_level=0
niter=30
tol_res_rel=1.e-2
tol_res_abs=1.e-2
tol_pot_rel=1.e-5
tol_pot_abs=1.e-5
tol_nwt_rel=1.e-5
tol_nwt_rel=1.e-5
tol_cg_rel=1.e-5
tol_vcy_rel=1.e-5
tol_list=( 0 0 0 0 0 0 0 )

defco_order=4

use_damping=F

exe=../dl_mg_test.exe

MPIOPTS=""

verbose=0

parallel_deviations=()
failed_run_list=() # list of failed runs (command return non-zero status)
failed_list=() # list of failed computation

use_steric=T

# command line args

for arg in "$@"
do
    case $arg in
	-model=*)
	    test_list="$(echo $arg | sed -e 's/-[a-z]*=//' -e 's/,/ /' )"
	    ;;
	-v) verbose=1 ;;
	-mpiexec=*)
	    MPIEXEC=${arg#*=}
	    ;;
        -mpiopts=*)
            MPIOPTS="${arg#*=}"
            ;;
	-par-mode=*)
	    parmode=${arg#*=}
	    ;;
	-defco-order=*)
	    defco_order=${arg#*=}
	    ;;
        -tol=*)
            tol_list=($(echo $arg | sed -e 's/-[a-z]*=//' -e 's/,/ /' ))
            ;;
        -exe=*)
            exe=${arg#*=}
            ;;
        -grid=*)
            ga=($(echo ${arg#*=} | sed -e 's/,/ /g' ))
            nx=${ga[0]}; ny=${ga[1]}; nz=${ga[2]}
	    ;;
	-use_cg=*)
	    use_cg=${arg#*=}
            ;;
	*) echo "Unknown flag $arg! quitting ..."; exit 1;;
    esac
done

# sanity tests

if [ ! -f "$exe" ]
then
  echo "$exe does not exist or is nort executable!"; exit 1
fi

case $parmode in
    mpi+omp)
	topolist="1-1-4 1-2-2 2-1-2 4-1-1"
	threadlist="1 2 3"
	;;
    omp)
	topolist="1-1-1"
	threadlist="1 2 3"
	;;
    serial)
	topolist="1-1-1"
	threadlist=1
	;;
    *) echo "Unknown parallel mode! quitting ..."; exit 1;;
esac

if [ ! ${tol_list[0]} = 0 ] ; then tol_res_rel=${tol_list[0]} ; fi
if [ ! ${tol_list[1]} = 0 ] ; then tol_res_abs=${tol_list[1]} ; fi
if [ ! ${tol_list[2]} = 0 ] ; then tol_pot_rel=${tol_list[2]} ; fi
if [ ! ${tol_list[3]} = 0 ] ; then tol_pot_abs=${tol_list[3]} ; fi
if [ ! ${tol_list[4]} = 0 ] ; then tol_nwt_rel=${tol_list[4]} ; fi
if [ ! ${tol_list[5]} = 0 ] ; then tol_cg_abs=${tol_list[5]} ; fi
if [ ! ${tol_list[6]} = 0 ] ; then tol_vcy_abs=${tol_list[5]} ; fi

# global test counter
counter=0
# index for failed runs ( as collected from log file)
ifail=0
#count the numner of failed  execution commands
nrunfailed=0

# tolerace for the deviation between the parallel and serial runs
tol_par=1.e-12

#shift necessary for PBC which uses n<D> - 1 grid points
dnx=0; dny=0; dnz=0

# neutralisation method for PBC, needs a test
neutralisation_method=none

echo "testing $test_list"

for model in $test_list
do
    counter_model=0
    case $model in
      pbc)
         dnx=1; dny=1; dnz=1
         MPI_period="T T T"
         ;;
      wirex)
         dnx=1; dny=0; dnz=0
         MPI_period="T F F"
         ;;
      wirez)
         dnx=0; dny=0; dnz=1
         MPI_period="F F T"
         ;;
      pbez)
         use_steric=F
         dnx=1; dny=1; dnz=0
         MPI_period="T T F"
         ;;
      pbe_sphere_hardcore_steric+fullpbc)
	 use_steric=T
         dnx=1; dny=1; dnz=1
         MPI_period="T T T"
	 neutralisation_method=ions_auto
	 model=pbe_sphere_hardcore_steric
	 ;;
       *)
         dnx=0; dny=0; dnz=0
         MPI_period="F F F"
         ;;
    esac

    for mpi_topo in $topolist
    do
        np=(${mpi_topo//-/ })

cat > input << EndOfInput
# global grid sizes
$((nx-dnx)) $((ny-dny)) $((nz-dnz))

# MPI grid
${np[0]} ${np[1]} ${np[2]}

# periodicity
$MPI_period

# use conjugate gradient
$use_cg

# number of iterations (defect  correction, multigrid)
 $niter $niter

# write a plane in a given section
-1 -1

# test name (laplace poisson ...), model
poisson $model

# tol (defect correction: residual {abs,rel}, pot (abs,rel}; newton rel; mg_rel )
$tol_res_rel $tol_res_abs $tol_pot_rel $tol_pot_abs $tol_nwt_rel $tol_cg_rel $tol_vcy_rel

# finite diference order for defect correction; use_damping
$defco_order $use_damping

# PBE temp n ions
300.0 2

# PBE concetration in Mol/l
0.1 0.1

# PBE charges (in electron charge)
1.0 -1.0

# PBE lambda,  take second value only first is T (for F uses default)
F 0.0

# PBE  linearized, use steric, errors_return, use full approx scheme
F $use_steric F F

# neutralisation method. one of the following
# none, jellium_unif, jellium_vacc ions_fixed, ions_auto, ions_auto_lin
$neutralisation_method

# ion ratios (for ion_fixed neutralisation, no need to add up to 1)
1.0 1.0

EndOfInput

	for nth in $threadlist
	do

	    counter=$(( ++counter ))
	    counter_model=$(( ++counter_model ))

	    export DL_MG_LOG_FILE=dl_mg_test_suite_"$model"_"$mpi_topo"_t"$nth".$counter

	    flog=$DL_MG_LOG_FILE
	    fout=out_dl_mg_suite_test.$counter

	    export OMP_NUM_THREADS=$nth
	    nproc=$((np[0] * np[1] * np[2]))
	    if (( verbose != 0 ))
            then
	       echo "test no $counter $counter_model"
               echo "model       $model"
               echo "fd order    $defco_order"
               echo "MPI         $mpi_topo"
               echo "OMP threads $nth"
	       echo "running $MPIEXEC -n $nproc $MPIOPTS  $exe > $fout 2>&1"
            fi
            if [ "$MPIEXEC" = nompi ]
            then
              if [ "$parmode" = mpi+omp ]
              then
                echo "inconsistentcy between MPIEXEC and parallel mode : $MPIEXEC $parmode"; exit 1
              else
                $exe  > $fout 2>&1
		cmderr=$?
              fi
            else
	      $MPIEXEC -n $nproc $MPIOPTS  $exe > $fout 2>&1
	      cmderr=$?
            fi

	    if [  "$cmderr" -ne 0  ]
	    then
		echo "failed execution with model $model, MPI $mpi_topo, threads $nth err code $cmderr"
		((++nrunfailed))
		failed_run_list[$nrunfailed]=$counter
		continue
	    fi
#collect some reference values
# was Newton iteration used ? Then we need to jump the MG report
	    have_newton=$(grep -i "Newton.*report" $flog)
	    if [ -n "$have_newton" ]
	    then
		# more elaborate version; get both values in one go
                rs_a=($(awk '/Newton.*report/,0 {if ($0 ~ /residual norm/) r=$3
                                                 if ($0 ~ /solution *norm/)  s=$3}
                              END {print r,s}' $flog))
		#residual=$(awk '/Newton.*report/,0 {if ($0 ~ /residual norm/) r= $4}' $flog)
		#solution=$(awk '/Newton.*report/,0 {if ($0 ~ /vector   norm/  print $3}' $flog)
	    else
		rs_a=($(awk '{if ($0 ~ /residual norm/) r=$3
                              if ($0 ~ /solution *norm/)  s=$3}
                             END {print r,s}' $flog))
	    fi
	    residual=${rs_a[0]}
	    solution=${rs_a[1]}
	    #residual=$(grep 'residual norm' $flog | awk '{print $4}')
	    #solution=$(grep 'vector   norm' $flog | awk '{print $3}')

# attention, failed value relies on a certain output format in the log file
	    failed=$(awk '/dl_mg_info: failed to converge/ || /DL_MG fatal error/' $flog)
	    if [ -n "$failed" ]
            then
               failed_list[ifail]="$model ${mpi_topo} $nth"
	       ((++ifail))
	       failed=''
	    fi
	    if (( counter_model == 1 ))
	    then
		residual_ref=${residual}
		solution_ref=${solution}
		res_diff=""
		sol_diff=""
		err_flag=""
	    else
		res_diff=$( echo $residual $residual_ref | awk '{printf("%20.12e", $1-$2)}' )
		sol_diff=$( echo $solution $solution_ref | awk '{printf("%20.12e", $1-$2)}' )
		ar=$( echo $res_diff | awk '{print ($1>0 ? $1 : -$1)}')
		as=$( echo $sol_diff | awk '{print ($1>0 ? $1 : -$1)}')
		err_flag=$( echo $ar $as | awk '{if ( $1 > '$tol_par' || $2 > '$tol_par' ) print "!!!"}')
		if [ -n "$err_flag" ]
		then
		    i=${#parallel_deviations[@]}
		    parallel_deviations[i]="$model ${mpi_topo} $nth"
		fi
	    fi

	    if (( verbose > 0 ))
	    then
		sol_err=$(grep 'solution error' $fout | awk '{print $3}')
		echo "error $sol_err, residual $residual $sol_diff $res_diff $err_flag"
		echo " "
	    fi
	    #we might need to keep these files
	    #rm $fout $flog
	done
    done
done

if (( ifail > 0 ))
then
    echo "The folowing computations have failed"
    for ((i=0; i<ifail;++i))
    do
	echo "${failed_list[i]}"
    done
fi

if (( nrunfailed > 0 ))
then
    echo "There were $nrunfailed failed runs:"
    echo "${failed_run_list[@]}"
fi

if (( nrunfailed == 0 && ifail == 0))
then
    echo ""
    echo "No failed runs!"
    echo ""
fi

echo "tolerance for deviations between parallel runs: $tol_par".
echo ""

n=${#parallel_deviations[@]}
if (( n > 0))
then
    echo "The following parallel runs differ significantly from the serial version"
    for ((i=0; i<n;++i))
    do echo  ${parallel_deviations[i]}
    done
fi
