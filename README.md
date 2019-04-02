PBBCache is a cache-partitioning simulator that relies on offline-collected application performance data (e.g., instructions per 
cycle, memory bandwidth consumption, etc.) to determine the degree of throughput, fairness or other relevant metrics for a workload
under a particular partitioning algorithm. A key feature of the simulator is its ability to efficiently determine the optimal
cache-partitioning solution for different optimizations objectives (e.g. throughput or fairness) using a parallel branch-and-bound algorithm.

The main purpose of the simulator is twofold. Firstly, it has been designed as a tool to efficiently assess the
real potential of existing cache-partitioning approaches and to identify their limitations, 
by comparing their effectiveness with that of the optimal solution provided by the simulator. 
Secondly, the ability to study the optimal solution for various optimization metrics makes the simulator
a powerful tool to aid in the design of new cache-partitioning techniques.

## Project dependencies

The project relies on Python 2.7.9 or greater and uses the following python modules:
* matplotlib
* pandas
* numpy
* sympy
* ipyparallel
* sklearn
* datetime
* pytz

Note: After installing ipyparllel, the binaries that start the cluster are not included in the path sometimes. To do so you must find where are they installed and add them to the PATH:

	find ~ -name 'ipcluster'
	export PATH=$PATH:$HOME/.local/bin

For more information on this library, recur to the [official documentation](https://ipyparallel.readthedocs.io/en/latest/)

## Getting started

A regular use of the simulator would include the following steps:

Include the project common functions module:

	. shrc

Start a local cluster with 8 cores:

	ipcluster start -n 8 --daemonize

Execute the simulator:

	./test/simulator_generic.py -s data/example_metrics.csv data/workloads.csv -a opt-stp,ucp -f table -b 35000 -p -m stalls -r 1
	W#   Algorithm           BenchID Name                 Cluster       Mask/NR_WAYS    SlowdownNB/ANTT   STP               Slowdown         
	W1   opt-stp             1       wupwise00            1             0x1(1)          1.02462           0.97598           1.02462          
	W1   opt-stp             2       zeusmp06             2             0x2(1)          1.03593           0.96532           1.03593          
	W1   opt-stp             3       xalancbmk06          3             0x7fc(9)        1.02027           0.98013           1.02027          
	W1   opt-stp             4       OVERALL              #36           0.1527s         1.02694           2.92142           1.01534          
	W1   ucp                 1       wupwise00            1             0x1(1)          1.02462           0.97598           1.02462          
	W1   ucp                 2       zeusmp06             2             0x2(1)          1.03593           0.96532           1.03593          
	W1   ucp                 3       xalancbmk06          3             0x7fc(9)        1.02027           0.98013           1.02027          
	W1   ucp                 4       OVERALL              #1            0.0011s         1.02694           2.92142           1.01534          
          


There are many available command line options, in this example we have used the most basic ones:

* -s specifies a metrics file, that contains a table where each row store the values of various runtime metrics
gathered with performance monitoring counters when a specific application runs alone with a fixed number of cache ways.
* The following parameter (required allways) is the workloads file, which specifies the different workloads (one per line) 
that will be considered for the simulation; each workload is encoded as a comma-separated list of program names.
* -a defines the set of partitioning schemes to simulate.
* -f is the output mode
* -p activates the parallel execution of the branch and bound algorithm to calculate the optimal partiotioning.
* -m specifies the bandwidth model to be applied, in this case the one based on the stall metric.
* -b defines the maximum attainable bandwidth (MB) in the platform to explore
* -r allows to only execute a defined range of workloads. It is defined as a comma-separated list and allows to define ranges with - (e.g., if we have 10 workloads and define the flag -r 1,3,5-7,9- we will end up running the workloads 1,3,5,6,7,9,10).

To find all of them you can run the following command:

	./test/simulator_generic.py -h
	usage: simulator_generic.py [-h] [-s SUMMARY_FILE] [-b MAX_BANDWIDTH] [-p]
	                            [-H] [-C] [-m BW_MODEL] [-a ALGORITHMS]
	                            [-f FORMAT] [-d] [-r USE_RANGE] [-O key=val]
	                            WORKLOAD_FILE
	
	Test main file for simulator (UCP, Whirlpool, etc.)
	
	positional arguments:
	  WORKLOAD_FILE         a workload file
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -s SUMMARY_FILE, --summary-file SUMMARY_FILE
	                        File with the data collected offline for every
	                        benchmark
	  -b MAX_BANDWIDTH, --max-bandwidth MAX_BANDWIDTH
	                        Max bandwidth observed for the target machine
	  -p, --parallel        Enable parallel search for the optimal
	  -H, --harness         Use harness workload file as input
	  -C, --generate-chart  Enable the STP-vs-UNF chart generator
	  -m BW_MODEL, --bw-model BW_MODEL
	                        Select BW model to be applied (['simple', 'stalls'])
	  -a ALGORITHMS, --algorithms ALGORITHMS
	                        Comma-separated algorithms (['simple', 'stalls'])
	  -f FORMAT, --format FORMAT
	                        Format style for output
	  -d, --debugging       Enable debugging mode (assume normal input files) and
	                        keep workload 0
	  -r USE_RANGE, --use-range USE_RANGE
	                        Pick selected workloads only by specifying a range
	  -O key=val, --option key=val
	                        Use generic algorithm-specific options

## Project contributors

* Juan Carlos Saez Alcaide (<jcsaezal@ucm.es>)
* Jorge Casas Hernan (<jorcasas@ucm.es>)
* Adrian Garcia Garcia (<adriagar@ucm.es>)
