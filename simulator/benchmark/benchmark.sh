#!/bin/bash
action=${1:-help}
max_threads=${2:-20}
bench=${3:-"normal"} ## another option is opt
initial_threads=${4:-${max_threads}}
malleability=${PBBCACHE_MALLEABILITY:-0} ## Malleability disabled by default
## Use one engine less so that the master process
## can use its own core
max_engines=$(( $max_threads - 1 ))
initial_engines=$(( $initial_threads - 1 ))

if [ "${PBBCACHE}" == "" ]; then
    echo "PBBCACHE envar not defined"
    env
    exit 1
fi

## Use an environment variable to force max threads 
if [ $malleability -eq 1 ] && [ "${GOMP_MAX_THREADS}" != "" ];  then
	max_engines=$(( ${GOMP_MAX_THREADS} - 1 ))
fi


## Extracted from launch-silent for initial binding
## The shell binds itself!
if [ "${HARNESS_BIND}" != "" ]; then
	b_vector=(${HARNESS_BIND})
	index=${BENCH_J_INDEX:-0} ## TO work in alone mode
	csettings=${b_vector[${index}]}
    if  [ "${csettings}" != "" ] && [ "${csettings}" != "-" ] ; then
		prefix_mask=${csettings:0:2}

		if [ "${prefix_mask}" == "0x" ]; then
			taskset -p ${csettings} $$ > /dev/null
		else
			## Interpret as list of CPUs
			taskset -cp ${csettings} $$ > /dev/null
		fi
    fi
fi

EXEC_PREFIX=""

#User-level numa settings 
if [ "${HARNESS_NUMACTL_SETTINGS}" != "" ]; then
	n_vector=(${HARNESS_NUMACTL_SETTINGS})
	index=${BENCH_J_INDEX:-0} ## TO work in alone mode
	nsettings=${n_vector[${index}]}
    if  [ "${nsettings}" != "" ] && [ "${nsettings}" != "-" ] ; then
		EXEC_PREFIX="numactl --${nsettings} -- "
    fi
fi

## Check for elasticity agent
if [ "$malleability" == "1" ] && [ ! -f ${PBBCACHE}/../util/elasticity_agent/elasticity_agent ]; then
    echo "Elasticity agent is not present. Did you forget to compile it?"
    exit 1
fi

if [ "${max_engines}" == "0" ] || [ "$bench"  == "debug" ]; then
	parallel=0
else
	parallel=1
fi


case $action in
    start) if [ $parallel -eq 1 ]; then
			echo "Starting $max_engines engines" 
	    	${EXEC_PREFIX} ${PBBCACHE}/test/run_venv3.sh ipcluster start -n ${max_engines} --daemonize
			sleep 12
      	   fi 
	   ;;
    run) 
    
	if [ "$bench"  == "debug" ]; then
		${PBBCACHE}/../util/elasticity_agent/elasticity_agent ${initial_engines} ${max_engines} sleep 20
	elif [ "$bench" == "normal" ] ||  [ "$bench" == "-" ] ; then
        if [ $parallel -eq 1 ]; then
			pflags="-P -O active_workers=${initial_engines} -O max_workers=${max_engines} -O malleability=${malleability}"
		else
			pflags=""
		fi

		## Independent "small" simulations
		${EXEC_PREFIX} ${PBBCACHE}/test/launch_sim_venv3.sh $pflags 	-a whirlpool-c,lfoc,lfoc+,kpart  \
        	-O benchmark_categories=${PBBCACHE}/data/volta_specnew_classification.csv \
        	-b 76000  -m simple  -O streaming_part_size=2  \
        	-s ${PBBCACHE}/data/volta2_turbo_off.300W.stalls.new.csv -f df benchmark/workloads20.txt 
    else
        if [ $parallel -eq 1 ]; then
				pflags="-p -O active_workers=${initial_engines} -O max_workers=${max_engines} -O malleability=${malleability}"
		else
			pflags=""
		fi

	## Big common simulation with independent tests
        ${EXEC_PREFIX} ${PBBCACHE}/test/launch_sim_venv3.sh $pflags -a optc:unf:min   \
        -O benchmark_categories=${PBBCACHE}/data/volta_specnew_classification.csv \
        -b 76000  -m simple   \
        -s ${PBBCACHE}/data/volta2_turbo_off.300W.stalls.new.csv -f df benchmark/workloads8.txt -r 1-4  ## FOur workloads only 
    fi
        ;;
    stop)  
	  if [ $parallel -eq 1 ]; then 
	    ${PBBCACHE}/test/run_venv3.sh ipcluster stop
       	  fi 
	   ;;
    help) echo "Usage: $0 [ start | run | stop | help ]"
        ;;
    *) echo "invalid option: $action"
        ;;
esac
