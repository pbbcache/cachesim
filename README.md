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

The project relies on Python 2.7.9 or greater and uses the following packages:
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

	./test/simulator_generic.py -s data/example_metrics.csv data/workoads.csv -a opt-stp,ucp,kpart -f table -b 35000 -p -m stalls -r 1

There are many available flags, in this example we have used the most basic ones:

* -s specifies a metrics file, that contains a table where each row store the values of various runtime metrics
gathered with performance monitoring counters when a specific application runs alone with a fixed number of cache ways.
* The following parameter (required allways) is the workloads file, which specifies the different workloads (one per line) 
that will be considered for the simulation; each workload is encoded as a comma-separated list of program names.
* -a defines the set of partitioning schemes to simulate.
* -f is the output mode
* -p activates the parallel execution of the branch and bound algorithm to calculate the optimal partiotioning.
* -m specifies the bandwidth model to be applied, in this case the one based on the stall metric.
* -r allows to only execute a defined range of workloads. It is defined as a comma-separated list and allows to define ranges with - (e.g., if we have 10 workloads and define the flag -r 1,3,5-7,9- we will end up running the workloads 1,3,5,6,7,9,10).

## Project contributors

* Juan Carlos Saez Alcaide (<jcsaezal@ucm.es>)
* Jorge Casas Hernan
* Adrian Garcia Garcia (<adriagar@ucm.es>)
